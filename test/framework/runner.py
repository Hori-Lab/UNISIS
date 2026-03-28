"""Test runner: orchestrates test execution across cases, modes, and levels."""

import os
import glob
import re
import shutil

from .utils import (
    run_command, ensure_dir, setup_work_dir, copy_file,
    check_prerequisites, generate_mpi_wrapper,
    generate_ref_label, write_ref_summary, resolve_reference_set,
)
from .toml_patcher import patch_toml, localize_file_paths
from .comparator import compare_out_files, compare_sampling
from .reporter import TestResult
from .builder import MODE_TO_BUILD


class TestRunner:
    """Orchestrate test execution across cases x modes x levels."""

    def __init__(self, config, builder, reporter, repo_root,
                 work_base, omp_threads=4, mpi_cmd="mpirun",
                 timeout=300, keep_work=False, verbose=False,
                 ref_label=None, condor=False, condor_memory_gb=4,
                 condor_poll_interval=30, condor_timeout=3600):
        self.config = config  # dict of test case configs
        self.builder = builder
        self.reporter = reporter
        self.repo_root = repo_root
        self.work_base = work_base
        self.omp_threads = omp_threads
        self.mpi_cmd = mpi_cmd
        self.timeout = timeout
        self.keep_work = keep_work
        self.verbose = verbose
        self.ref_label = ref_label  # None = auto-select latest
        self.run_results = {}  # (case, mode) -> work_dir for reuse
        # HTCondor settings
        self.condor = condor
        self.condor_memory_gb = condor_memory_gb
        self.condor_poll_interval = condor_poll_interval
        self.condor_timeout = condor_timeout
        self._condor_batch = None       # CondorBatch, set in condor_run_all
        self._condor_batch_mode = False  # True during batch submit pass

    def run_all(self, cases, modes, levels, print_header=True):
        """Run all selected tests.

        Parameters
        ----------
        cases : list[str]
            Test case names.
        modes : list[str]
            Parallelization modes.
        levels : list[str]
            Test levels.
        print_header : bool
            Whether to print the results table header (default True).
            Set to False when called internally to avoid duplicate headers.
        """
        if print_header:
            self.reporter.print_header()

        for case_name in cases:
            case = self.config[case_name]
            valid_modes = [m for m in modes if m in case["modes"]]

            # Check prerequisites once per case
            missing = check_prerequisites(case, self.repo_root)
            if missing:
                for level in levels:
                    mode_str = valid_modes[0] if valid_modes else "---"
                    self.reporter.add(TestResult(
                        case_name, mode_str, level, TestResult.SKIP,
                        f"missing: {', '.join(missing)}"
                    ))
                continue

            # Level: run
            if "run" in levels:
                for mode in valid_modes:
                    self._run_smoke(case_name, case, mode)

            # Level: consistency (needs multiple modes)
            if "consistency" in levels and len(valid_modes) >= 2:
                self._run_consistency(case_name, case, valid_modes)

            # Level: regression
            if "regression" in levels:
                for mode in valid_modes:
                    self._run_regression(case_name, case, mode)

            # Level: restart
            if "restart" in levels and case.get("restart_test"):
                mode = valid_modes[0] if valid_modes else None
                if mode:
                    self._run_restart(case_name, case, mode)

            # Level: sampling
            if "sampling" in levels and case.get("sampling_properties"):
                mode = valid_modes[0] if valid_modes else None
                if mode:
                    self._run_sampling(case_name, case, mode)

            # Level: gradient
            if "gradient" in levels and case.get("gradient_test"):
                mode = valid_modes[0] if valid_modes else None
                if mode:
                    self._run_gradient(case_name, case, mode)

    def _get_work_dir(self, case_name, mode, suffix=""):
        """Get or create a work directory for a test run."""
        tag = f"{case_name}_{mode}"
        if suffix:
            tag += f"_{suffix}"
        return os.path.join(self.work_base, tag)

    def _ensure_run(self, case_name, case, mode, nstep_override=None,
                    suffix="", extra_toml_patches=None, extra_args=""):
        """Ensure a test case has been run for a given mode.

        Returns (work_dir, success) tuple. Reuses previous run if available.
        """
        key = (case_name, mode, suffix)
        if key in self.run_results:
            return self.run_results[key]

        work_dir = self._get_work_dir(case_name, mode, suffix)
        ensure_dir(work_dir)

        # Setup input files
        source_dir = os.path.join(self.repo_root, case["source_dir"])

        # Collect all files to symlink
        files_to_link = list(case.get("required_files", []))
        # Also link any restart-related TOMLs
        for f in case.get("restart_input_tomls", []):
            if f not in files_to_link:
                files_to_link.append(f)

        setup_work_dir(work_dir, case["source_dir"], files_to_link, self.repo_root)

        # Copy and patch TOML
        src_toml = os.path.join(source_dir, case["input_toml"])
        dst_toml = os.path.join(work_dir, case["input_toml"])

        patches = {}
        # Set output prefix to be in work dir
        patches["prefix"] = f"./{case['output_prefix']}"

        if nstep_override is not None:
            patches["nstep"] = nstep_override
        if "nstep_save" in case and nstep_override is not None:
            patches["nstep_save"] = case["nstep_save"]
        if "nstep_save_rst" in case and nstep_override is not None:
            patches["nstep_save_rst"] = case["nstep_save_rst"]
        if extra_toml_patches:
            patches.update(extra_toml_patches)

        patch_toml(src_toml, dst_toml, patches)
        localize_file_paths(dst_toml, dst_toml)

        # Build run command
        exe = self.builder.get_exe(mode)
        if exe is None:
            self.run_results[key] = (work_dir, False)
            return work_dir, False

        if self.condor and self._condor_batch_mode:
            # Batch submit pass: submit to HTCondor without waiting.
            # Store None as a pending sentinel; analysis deferred until after
            # wait_all() updates run_results with the actual exit code.
            self._condor_submit(case, mode, exe, work_dir,
                                case["input_toml"], extra_args, key)
            self.run_results[key] = (work_dir, None)
            return work_dir, None

        success = self._execute(case, mode, exe, work_dir,
                                case["input_toml"], extra_args)
        self.run_results[key] = (work_dir, success)
        return work_dir, success

    def _build_exec_cmd_and_env(self, case, mode, exe, work_dir,
                                input_toml, extra_args=""):
        """Build the simulation command string and environment dict.

        Returns (cmd, env) without executing anything.
        """
        env = {}
        use_mpi = mode in ("mpi", "mpi_omp")

        if mode in ("serial", "omp1", "mpi"):
            env["OMP_NUM_THREADS"] = "1"
        elif mode in ("ompN", "mpi_omp"):
            env["OMP_NUM_THREADS"] = str(self.omp_threads)

        if use_mpi:
            wrapper = generate_mpi_wrapper(work_dir, exe, input_toml, extra_args)
            ranks = case.get("mpi_ranks", 4)
            cmd = f"{self.mpi_cmd} -n {ranks} {wrapper}"
        else:
            cmd = f"{exe} {input_toml}"
            if extra_args:
                cmd += f" {extra_args}"

        return cmd, env

    def _execute(self, case, mode, exe, work_dir, input_toml, extra_args=""):
        """Execute the simulation locally.

        Returns True if exit code is 0.
        """
        cmd, env = self._build_exec_cmd_and_env(
            case, mode, exe, work_dir, input_toml, extra_args)

        rc, stdout, stderr, elapsed = run_command(
            cmd, cwd=work_dir, env=env, timeout=self.timeout
        )

        if self.verbose and (stdout or stderr):
            if stdout:
                print(f"    stdout: {stdout[:500]}")
            if stderr:
                print(f"    stderr: {stderr[:500]}")

        return rc == 0

    def _condor_submit(self, case, mode, exe, work_dir,
                       input_toml, extra_args, key):
        """Submit a simulation job to HTCondor (no wait).

        Adds the resulting CondorJob to self._condor_batch.
        Returns True to signal successful submission.
        """
        from .condor import prepare_and_submit
        cmd, env = self._build_exec_cmd_and_env(
            case, mode, exe, work_dir, input_toml, extra_args)
        mpi_ranks = case.get("mpi_ranks", 4)
        job = prepare_and_submit(
            work_dir=work_dir,
            case_name=key[0],
            mode=mode,
            cmd=cmd,
            env_vars=env,
            omp_threads=self.omp_threads,
            mpi_ranks=mpi_ranks,
            memory_gb=self.condor_memory_gb,
            key=key,
            verbose=self.verbose,
        )
        self._condor_batch.add(job)
        return True

    def _condor_submit_and_wait(self, case, mode, exe, work_dir,
                                input_toml, extra_args="", label="job"):
        """Submit a single job to HTCondor and block until it completes.

        Used for sequential restart stages. Returns True if exit code is 0.
        """
        from .condor import CondorBatch, prepare_and_submit
        key = (label, mode, work_dir)
        cmd, env = self._build_exec_cmd_and_env(
            case, mode, exe, work_dir, input_toml, extra_args)
        mpi_ranks = case.get("mpi_ranks", 4)
        job = prepare_and_submit(
            work_dir=work_dir,
            case_name=label,
            mode=mode,
            cmd=cmd,
            env_vars=env,
            omp_threads=self.omp_threads,
            mpi_ranks=mpi_ranks,
            memory_gb=self.condor_memory_gb,
            key=key,
            verbose=self.verbose,
        )
        batch = CondorBatch(
            poll_interval=self.condor_poll_interval,
            timeout=self.condor_timeout,
            verbose=self.verbose,
        )
        batch.add(job)
        results = batch.wait_all()
        return results[key] == 0

    def _find_out_files(self, work_dir, case):
        """Find .out files in a work directory."""
        prefix = case["output_prefix"]
        pattern = os.path.join(work_dir, f"{prefix}*.out")
        files = sorted(glob.glob(pattern))
        # Also try just prefix.out
        single = os.path.join(work_dir, f"{prefix}.out")
        if os.path.isfile(single) and single not in files:
            files.insert(0, single)
        return files

    # ---- Level implementations ----

    def _run_smoke(self, case_name, case, mode):
        """Level 1: smoke test."""
        # Build if needed
        try:
            self.builder.build_for_mode(mode)
        except RuntimeError as e:
            self.reporter.add(TestResult(
                case_name, mode, "run", TestResult.FAIL,
                f"build failed: {e}"
            ))
            return

        nstep = case.get("short_nstep")
        work_dir, success = self._ensure_run(case_name, case, mode,
                                             nstep_override=nstep)

        if success is None:  # submitted to condor, analysis deferred
            return

        if not success:
            # Check stderr for details
            err_files = glob.glob(os.path.join(work_dir, "err.*"))
            detail = "non-zero exit"
            if err_files:
                with open(err_files[0]) as f:
                    detail = f.read().strip()[:200] or detail
            self.reporter.add(TestResult(
                case_name, mode, "run", TestResult.FAIL, detail
            ))
            return

        out_files = self._find_out_files(work_dir, case)
        if not out_files:
            self.reporter.add(TestResult(
                case_name, mode, "run", TestResult.FAIL,
                "no .out file produced"
            ))
            return

        self.reporter.add(TestResult(
            case_name, mode, "run", TestResult.PASS
        ))

    def _run_consistency(self, case_name, case, modes):
        """Level 2: consistency across modes."""
        nstep = case.get("regression_nstep") or case.get("short_nstep")

        # Ensure all modes have been run
        work_dirs = {}
        for mode in modes:
            try:
                self.builder.build_for_mode(mode)
            except RuntimeError:
                continue
            wd, ok = self._ensure_run(case_name, case, mode,
                                      nstep_override=nstep)
            if ok is None:  # submitted to condor, analysis deferred
                continue
            if ok:
                work_dirs[mode] = wd

        if len(work_dirs) < 2:
            self.reporter.add(TestResult(
                case_name, "---", "consistency", TestResult.SKIP,
                f"need >=2 modes, got {len(work_dirs)}"
            ))
            return

        # Compare pairs
        mode_list = list(work_dirs.keys())
        all_passed = True
        details = []

        for i in range(len(mode_list)):
            for j in range(i + 1, len(mode_list)):
                m1, m2 = mode_list[i], mode_list[j]
                files1 = self._find_out_files(work_dirs[m1], case)
                files2 = self._find_out_files(work_dirs[m2], case)

                if not files1 or not files2:
                    details.append(f"{m1}~{m2}: missing .out")
                    all_passed = False
                    continue

                # Compare corresponding files
                for f1, f2 in zip(sorted(files1), sorted(files2)):
                    # Use tighter tolerance for serial vs omp1
                    if {m1, m2} == {"serial", "omp1"}:
                        atol, rtol = 1e-12, 1e-12
                    else:
                        atol, rtol = 1e-8, 1e-6

                    result = compare_out_files(f1, f2, atol=atol, rtol=rtol)
                    tag = f"{m1}~{m2}"
                    if result.passed:
                        details.append(f"{tag}: rdiff={result.max_rel_diff:.1e}")
                    else:
                        details.append(f"{tag}: MISMATCH {result.summary}")
                        all_passed = False

        status = TestResult.PASS if all_passed else TestResult.FAIL
        self.reporter.add(TestResult(
            case_name, "---", "consistency", status,
            "; ".join(details)
        ))

    def _run_regression(self, case_name, case, mode):
        """Level 3: regression against reference data."""
        ref_base = os.path.join(self.repo_root, "test", "reference")
        label, ref_set_dir = resolve_reference_set(ref_base, self.ref_label)

        if ref_set_dir is None:
            self.reporter.add(TestResult(
                case_name, mode, "regression", TestResult.SKIP,
                "no reference data set found"
            ))
            return

        ref_dir = os.path.join(ref_set_dir, case_name)
        if not os.path.isdir(ref_dir):
            self.reporter.add(TestResult(
                case_name, mode, "regression", TestResult.SKIP,
                f"no reference data for this case in {label}"
            ))
            return

        nstep = case.get("regression_nstep")
        try:
            self.builder.build_for_mode(mode)
        except RuntimeError as e:
            self.reporter.add(TestResult(
                case_name, mode, "regression", TestResult.FAIL,
                f"build failed"
            ))
            return

        work_dir, success = self._ensure_run(case_name, case, mode,
                                             nstep_override=nstep)
        if success is None:  # submitted to condor, analysis deferred
            return

        if not success:
            self.reporter.add(TestResult(
                case_name, mode, "regression", TestResult.FAIL,
                "run failed"
            ))
            return

        out_files = self._find_out_files(work_dir, case)
        ref_files = sorted(glob.glob(os.path.join(ref_dir, "*.out")))

        if not out_files or not ref_files:
            self.reporter.add(TestResult(
                case_name, mode, "regression", TestResult.FAIL,
                f"out={len(out_files)}, ref={len(ref_files)} files"
            ))
            return

        tol = case["tolerances"]
        all_passed = True
        worst_rdiff = 0.0

        for actual, ref in zip(sorted(out_files), ref_files):
            result = compare_out_files(actual, ref,
                                       atol=tol["atol"], rtol=tol["rtol"])
            worst_rdiff = max(worst_rdiff, result.max_rel_diff)
            if not result.passed:
                all_passed = False

        if all_passed:
            self.reporter.add(TestResult(
                case_name, mode, "regression", TestResult.PASS,
                f"max_rdiff={worst_rdiff:.2e} (ref: {label})"
            ))
        else:
            self.reporter.add(TestResult(
                case_name, mode, "regression", TestResult.FAIL,
                f"{result.summary} (ref: {label})"
            ))

    def _run_restart(self, case_name, case, mode):
        """Level 4: restart consistency test.

        Run full simulation, then run in two halves with restart.
        Compare output from the second half.
        """
        try:
            self.builder.build_for_mode(mode)
        except RuntimeError:
            self.reporter.add(TestResult(
                case_name, mode, "restart", TestResult.FAIL, "build failed"
            ))
            return

        nstep = case.get("regression_nstep") or case.get("short_nstep", 100)
        nstep_save = case.get("nstep_save", 10)
        nstep_save_rst = case.get("nstep_save_rst", nstep // 2)
        half = nstep // 2
        # Ensure half is a multiple of nstep_save
        half = (half // nstep_save) * nstep_save
        if half < nstep_save:
            half = nstep_save

        use_mpi = mode in ("mpi", "mpi_omp")
        n_replicas = case.get("n_replicas", 1)

        # --- Stage 1: full continuous run ---
        work_full = self._get_work_dir(case_name, mode, "restart_full")
        ensure_dir(work_full)
        source_dir = os.path.join(self.repo_root, case["source_dir"])

        files_to_link = list(case.get("required_files", []))
        setup_work_dir(work_full, case["source_dir"], files_to_link, self.repo_root)

        src_toml = os.path.join(source_dir, case["input_toml"])
        dst_toml = os.path.join(work_full, case["input_toml"])
        patch_toml(src_toml, dst_toml, {
            "prefix": f"./{case['output_prefix']}",
            "nstep": nstep,
            "nstep_save": nstep_save,
            "nstep_save_rst": nstep_save_rst,
        })
        localize_file_paths(dst_toml, dst_toml)

        exe = self.builder.get_exe(mode)
        ok_full = self._execute(case, mode, exe, work_full, case["input_toml"])

        if not ok_full:
            self.reporter.add(TestResult(
                case_name, mode, "restart", TestResult.FAIL,
                "full run failed"
            ))
            return

        # --- Stage 2: first half ---
        work_s1 = self._get_work_dir(case_name, mode, "restart_s1")
        ensure_dir(work_s1)
        setup_work_dir(work_s1, case["source_dir"], files_to_link, self.repo_root)

        prefix_s1 = f"{case['output_prefix']}_s1"
        s1_toml = os.path.join(work_s1, case["input_toml"])
        patch_toml(src_toml, s1_toml, {
            "prefix": f"./{prefix_s1}",
            "nstep": half,
            "nstep_save": nstep_save,
            "nstep_save_rst": half,  # save restart at the end
        })
        localize_file_paths(s1_toml, s1_toml)

        ok_s1 = self._execute(case, mode, exe, work_s1, case["input_toml"])
        if not ok_s1:
            self.reporter.add(TestResult(
                case_name, mode, "restart", TestResult.FAIL,
                "stage 1 (first half) failed"
            ))
            return

        # --- Find and prepare restart file ---
        if use_mpi and n_replicas > 1:
            # Concatenate per-replica restart files
            rst_pattern = os.path.join(work_s1, f"{prefix_s1}_0*.rst")
            rst_files = sorted(glob.glob(rst_pattern))
            if not rst_files:
                self.reporter.add(TestResult(
                    case_name, mode, "restart", TestResult.FAIL,
                    "no restart files from stage 1"
                ))
                return
            combined_rst = os.path.join(work_s1, f"{prefix_s1}.rst")
            with open(combined_rst, "wb") as out:
                for rf in rst_files:
                    with open(rf, "rb") as inp:
                        out.write(inp.read())
            rst_path = combined_rst
        else:
            # Single restart file
            rst_pattern = os.path.join(work_s1, f"{prefix_s1}*.rst")
            rst_files = sorted(glob.glob(rst_pattern))
            if not rst_files:
                self.reporter.add(TestResult(
                    case_name, mode, "restart", TestResult.FAIL,
                    "no restart files from stage 1"
                ))
                return
            rst_path = rst_files[-1]  # last one

        # --- Stage 3: second half from restart ---
        work_s2 = self._get_work_dir(case_name, mode, "restart_s2")
        ensure_dir(work_s2)
        setup_work_dir(work_s2, case["source_dir"], files_to_link, self.repo_root)

        # Copy restart file into stage 2 work dir
        rst_basename = os.path.basename(rst_path)
        rst_in_s2 = os.path.join(work_s2, rst_basename)
        copy_file(rst_path, rst_in_s2)

        prefix_s2 = f"{case['output_prefix']}_s2"
        s2_toml = os.path.join(work_s2, case["input_toml"])
        patch_toml(src_toml, s2_toml, {
            "prefix": f"./{prefix_s2}",
            "nstep": nstep,
            "nstep_save": nstep_save,
            "nstep_save_rst": nstep,  # don't need restart from this
        })
        localize_file_paths(s2_toml, s2_toml)

        ok_s2 = self._execute(case, mode, exe, work_s2,
                              case["input_toml"],
                              extra_args=rst_basename)
        if not ok_s2:
            self.reporter.add(TestResult(
                case_name, mode, "restart", TestResult.FAIL,
                "stage 2 (restart) failed"
            ))
            return

        # --- Compare ---
        full_outs = sorted(glob.glob(
            os.path.join(work_full, f"{case['output_prefix']}*.out")
        ))
        s2_outs = sorted(glob.glob(
            os.path.join(work_s2, f"{prefix_s2}*.out")
        ))

        if not full_outs or not s2_outs:
            self.reporter.add(TestResult(
                case_name, mode, "restart", TestResult.FAIL,
                f"missing outputs: full={len(full_outs)}, s2={len(s2_outs)}"
            ))
            return

        all_passed = True
        worst_rdiff = 0.0

        for f_full, f_s2 in zip(full_outs, s2_outs):
            result = compare_out_files(
                f_s2, f_full,
                atol=1e-10, rtol=1e-10,
                align_steps=True
            )
            worst_rdiff = max(worst_rdiff, result.max_rel_diff)
            if not result.passed:
                all_passed = False

        if all_passed:
            self.reporter.add(TestResult(
                case_name, mode, "restart", TestResult.PASS,
                f"max_rdiff={worst_rdiff:.2e}"
            ))
        else:
            self.reporter.add(TestResult(
                case_name, mode, "restart", TestResult.FAIL,
                result.summary
            ))

    def _run_sampling(self, case_name, case, mode):
        """Level 5: thermodynamic sampling consistency."""
        props = case.get("sampling_properties", [])
        if not props:
            self.reporter.add(TestResult(
                case_name, mode, "sampling", TestResult.SKIP,
                "no sampling properties defined"
            ))
            return

        try:
            self.builder.build_for_mode(mode)
        except RuntimeError:
            self.reporter.add(TestResult(
                case_name, mode, "sampling", TestResult.FAIL, "build failed"
            ))
            return

        nstep = case.get("sampling_nstep", 100000)
        work_dir, success = self._ensure_run(
            case_name, case, mode, nstep_override=nstep, suffix="sampling"
        )

        if success is None:  # submitted to condor, analysis deferred
            return

        if not success:
            self.reporter.add(TestResult(
                case_name, mode, "sampling", TestResult.FAIL, "run failed"
            ))
            return

        out_files = self._find_out_files(work_dir, case)
        if not out_files:
            self.reporter.add(TestResult(
                case_name, mode, "sampling", TestResult.FAIL, "no .out file"
            ))
            return

        all_passed = True
        details = []

        for col, expected, tol_sigma in props:
            passed, detail = compare_sampling(
                out_files[0], col, expected, tol_sigma
            )
            details.append(f"col{col}: {detail}")
            if not passed:
                all_passed = False

        status = TestResult.PASS if all_passed else TestResult.FAIL
        self.reporter.add(TestResult(
            case_name, mode, "sampling", status,
            "; ".join(details)
        ))

    def _run_gradient(self, case_name, case, mode):
        """Level 6: gradient correctness (energy-force consistency)."""
        try:
            self.builder.build_for_mode(mode)
        except RuntimeError:
            self.reporter.add(TestResult(
                case_name, mode, "gradient", TestResult.FAIL, "build failed"))
            return

        work_dir, success = self._ensure_run(
            case_name, case, mode,
            suffix="gradient",
            extra_toml_patches={"type": "CHECK_FORCE"},
        )

        if success is None:  # submitted to condor, analysis deferred
            return
        if not success:
            self.reporter.add(TestResult(
                case_name, mode, "gradient", TestResult.FAIL, "run failed"))
            return

        out_files = self._find_out_files(work_dir, case)
        if not out_files:
            self.reporter.add(TestResult(
                case_name, mode, "gradient", TestResult.FAIL, "no output file"))
            return

        passed, detail = self._parse_gradient_results(out_files[0])
        status = TestResult.PASS if passed else TestResult.FAIL
        self.reporter.add(TestResult(case_name, mode, "gradient", status, detail))

    def _parse_gradient_results(self, out_file, tolerance=0.001):
        """Parse CHECK_FORCE output for energy-force differences.

        Returns (passed: bool, detail: str).
        """
        pattern = re.compile(
            r'^f_energy - f_force\s+(\d+)\s+([\d.eE+\-]+)\s+([\d.eE+\-]+)\s+([\d.eE+\-]+)')
        n_checked = 0
        failures = []
        with open(out_file) as f:
            for line in f:
                m = pattern.match(line)
                if m:
                    n_checked += 1
                    imp = int(m.group(1))
                    diffs = [abs(float(m.group(i))) for i in (2, 3, 4)]
                    if any(d > tolerance for d in diffs):
                        failures.append((imp, max(diffs)))
        if n_checked == 0:
            return False, "no gradient check lines found in output"
        if failures:
            worst_imp, worst_diff = max(failures, key=lambda x: x[1])
            return False, (
                f"{len(failures)}/{n_checked} atoms exceed tol={tolerance} "
                f"(worst: atom {worst_imp}, diff={worst_diff:.6f})")
        return True, f"all {n_checked} atoms within tol={tolerance}"

    # ---- HTCondor batch dispatch ----

    def condor_run_all(self, cases, modes, levels):
        """Run all selected tests by submitting jobs to HTCondor in batch.

        If a manifest from a prior --condor-submit-only exists, auto-detects
        completed jobs and skips re-submission. Otherwise submits all jobs,
        waits, then runs analysis.

        Parameters
        ----------
        cases : list[str]
        modes : list[str]
        levels : list[str]
        """
        from .condor import (
            CondorBatch, MANIFEST_FILENAME, read_manifest,
        )

        manifest_path = os.path.join(self.work_base, MANIFEST_FILENAME)
        manifest = read_manifest(manifest_path)

        if manifest:
            self._condor_collect_from_manifest(
                manifest, manifest_path, cases, modes, levels)
            return

        # No manifest — original behaviour: submit + wait + analyze.
        non_restart_levels = [l for l in levels if l != "restart"]
        do_restart = "restart" in levels

        self._condor_batch = CondorBatch(
            poll_interval=self.condor_poll_interval,
            timeout=self.condor_timeout,
            verbose=self.verbose,
        )

        if non_restart_levels:
            self._condor_batch_mode = True
            self.run_all(cases, modes, non_restart_levels, print_header=False)
            self._condor_batch_mode = False

            n_submitted = len(self._condor_batch.jobs)
            if n_submitted > 0:
                # Also write manifest for crash-recovery.
                self._write_manifest_from_batch(
                    manifest_path, cases, modes, levels)

                print(f"  Submitted {n_submitted} job(s) to HTCondor. "
                      f"Waiting (poll every {self.condor_poll_interval}s)...")
                results = self._condor_batch.wait_all()

                for key, exit_code in results.items():
                    work_dir, _ = self.run_results[key]
                    self.run_results[key] = (work_dir, exit_code == 0)

                # Clean up manifest after successful wait.
                if os.path.isfile(manifest_path):
                    os.remove(manifest_path)

            self.run_all(cases, modes, non_restart_levels)

        if do_restart:
            for case_name in cases:
                case = self.config[case_name]
                if not case.get("restart_test"):
                    continue
                missing = check_prerequisites(case, self.repo_root)
                if missing:
                    self.reporter.add(TestResult(
                        case_name, "---", "restart", TestResult.SKIP,
                        f"missing: {', '.join(missing)}"
                    ))
                    continue
                valid_modes = [m for m in modes if m in case["modes"]]
                mode = valid_modes[0] if valid_modes else None
                if mode:
                    self._run_restart_condor(case_name, case, mode)

    def condor_submit_only(self, cases, modes, levels):
        """Submit all jobs to HTCondor and exit immediately.

        Writes a manifest file so that a subsequent --condor run can
        collect the results without re-submitting.

        Parameters
        ----------
        cases : list[str]
        modes : list[str]
        levels : list[str]
        """
        from .condor import CondorBatch, MANIFEST_FILENAME

        non_restart_levels = [l for l in levels if l != "restart"]
        do_restart = "restart" in levels

        self._condor_batch = CondorBatch(
            poll_interval=self.condor_poll_interval,
            timeout=self.condor_timeout,
            verbose=self.verbose,
        )

        # Submit non-restart jobs.
        if non_restart_levels:
            self._condor_batch_mode = True
            self.run_all(cases, modes, non_restart_levels, print_header=False)
            self._condor_batch_mode = False

        # Submit restart stages full + s1 (s2 deferred to collect).
        restart_jobs = {}
        if do_restart:
            restart_jobs = self._submit_restart_partial(cases, modes)

        # Write manifest.
        manifest_path = os.path.join(self.work_base, MANIFEST_FILENAME)
        self._write_manifest_from_batch(
            manifest_path, cases, modes, levels, restart_jobs)

        n_total = len(self._condor_batch.jobs) + len(restart_jobs)
        print(f"  Submitted {n_total} job(s). Manifest: {manifest_path}")
        print(f"  Run with --condor (same args) to collect results.")

    def _write_manifest_from_batch(self, manifest_path, cases, modes, levels,
                                   restart_jobs=None):
        """Build a manifest from the current _condor_batch and write it."""
        from .condor import write_manifest, _key_to_str
        jobs = {}
        for job in self._condor_batch.jobs:
            key_str = _key_to_str(job.key)
            jobs[key_str] = {
                "work_dir": job.work_dir,
                "cluster_id": job.cluster_id,
                "log_path": job.log_path,
            }
        write_manifest(manifest_path, cases, modes, levels,
                       jobs, restart_jobs)

    def _submit_restart_partial(self, cases, modes):
        """Submit restart stages full + s1 (not s2) for all eligible cases.

        Returns a restart_jobs dict for the manifest.
        """
        from .condor import prepare_and_submit, _key_to_str
        restart_jobs = {}

        for case_name in cases:
            case = self.config[case_name]
            if not case.get("restart_test"):
                continue
            missing = check_prerequisites(case, self.repo_root)
            if missing:
                continue
            valid_modes = [m for m in modes if m in case["modes"]]
            mode = valid_modes[0] if valid_modes else None
            if not mode:
                continue

            try:
                self.builder.build_for_mode(mode)
            except RuntimeError:
                continue

            exe = self.builder.get_exe(mode)
            if exe is None:
                continue

            nstep = case.get("regression_nstep") or case.get("short_nstep", 100)
            nstep_save = case.get("nstep_save", 10)
            nstep_save_rst = case.get("nstep_save_rst", nstep // 2)
            half = nstep // 2
            half = (half // nstep_save) * nstep_save
            if half < nstep_save:
                half = nstep_save

            source_dir = os.path.join(self.repo_root, case["source_dir"])
            src_toml = os.path.join(source_dir, case["input_toml"])
            files_to_link = list(case.get("required_files", []))

            # --- Stage full ---
            work_full = self._get_work_dir(case_name, mode, "restart_full")
            ensure_dir(work_full)
            setup_work_dir(work_full, case["source_dir"],
                           files_to_link, self.repo_root)
            dst_toml = os.path.join(work_full, case["input_toml"])
            patch_toml(src_toml, dst_toml, {
                "prefix": f"./{case['output_prefix']}",
                "nstep": nstep,
                "nstep_save": nstep_save,
                "nstep_save_rst": nstep_save_rst,
            })
            localize_file_paths(dst_toml, dst_toml)

            key_full = (case_name, mode, "restart_full")
            cmd_full, env_full = self._build_exec_cmd_and_env(
                case, mode, exe, work_full, case["input_toml"])
            mpi_ranks = case.get("mpi_ranks", 4)
            job_full = prepare_and_submit(
                work_dir=work_full, case_name=f"{case_name}_rst_full",
                mode=mode, cmd=cmd_full, env_vars=env_full,
                omp_threads=self.omp_threads, mpi_ranks=mpi_ranks,
                memory_gb=self.condor_memory_gb, key=key_full,
                verbose=self.verbose)
            restart_jobs[_key_to_str(key_full)] = {
                "work_dir": job_full.work_dir,
                "cluster_id": job_full.cluster_id,
                "log_path": job_full.log_path,
            }

            # --- Stage s1 ---
            work_s1 = self._get_work_dir(case_name, mode, "restart_s1")
            ensure_dir(work_s1)
            setup_work_dir(work_s1, case["source_dir"],
                           files_to_link, self.repo_root)
            prefix_s1 = f"{case['output_prefix']}_s1"
            s1_toml = os.path.join(work_s1, case["input_toml"])
            patch_toml(src_toml, s1_toml, {
                "prefix": f"./{prefix_s1}",
                "nstep": half,
                "nstep_save": nstep_save,
                "nstep_save_rst": half,
            })
            localize_file_paths(s1_toml, s1_toml)

            key_s1 = (case_name, mode, "restart_s1")
            cmd_s1, env_s1 = self._build_exec_cmd_and_env(
                case, mode, exe, work_s1, case["input_toml"])
            job_s1 = prepare_and_submit(
                work_dir=work_s1, case_name=f"{case_name}_rst_s1",
                mode=mode, cmd=cmd_s1, env_vars=env_s1,
                omp_threads=self.omp_threads, mpi_ranks=mpi_ranks,
                memory_gb=self.condor_memory_gb, key=key_s1,
                verbose=self.verbose)
            restart_jobs[_key_to_str(key_s1)] = {
                "work_dir": job_s1.work_dir,
                "cluster_id": job_s1.cluster_id,
                "log_path": job_s1.log_path,
            }

        return restart_jobs

    def _condor_collect_from_manifest(self, manifest, manifest_path,
                                      cases, modes, levels):
        """Collect results from a prior condor submission using its manifest.

        Checks which jobs are done, waits for pending ones, runs analysis,
        and handles deferred restart stage s2.
        """
        from .condor import (
            CondorBatch, CondorJob, collect_results_from_manifest,
            _str_to_key,
        )

        non_restart_levels = [l for l in levels if l != "restart"]
        do_restart = "restart" in levels

        status = collect_results_from_manifest(manifest)
        n_completed = len(status["completed"])
        n_pending = len(status["pending"])
        print(f"  Manifest found: {n_completed} completed, "
              f"{n_pending} pending job(s).")

        # Populate run_results for completed jobs.
        for key_str, exit_code in status["completed"].items():
            key = _str_to_key(key_str)
            info = (manifest.get("jobs", {}).get(key_str)
                    or manifest.get("restart_jobs", {}).get(key_str))
            if info:
                self.run_results[key] = (info["work_dir"], exit_code == 0)

        # Wait for any still-pending jobs.
        if status["pending"]:
            batch = CondorBatch(
                poll_interval=self.condor_poll_interval,
                timeout=self.condor_timeout,
                verbose=self.verbose,
            )
            for key_str, info in status["pending"].items():
                key = _str_to_key(key_str)
                job = CondorJob(
                    cluster_id=info["cluster_id"],
                    work_dir=info["work_dir"],
                    log_path=info["log_path"],
                    key=key,
                )
                batch.add(job)

            print(f"  Waiting for {n_pending} pending job(s)...")
            results = batch.wait_all()
            for key, exit_code in results.items():
                info_str = "|".join(str(k) for k in key)
                info = (manifest.get("jobs", {}).get(info_str)
                        or manifest.get("restart_jobs", {}).get(info_str))
                if info:
                    self.run_results[key] = (info["work_dir"], exit_code == 0)

        # Analysis pass for non-restart levels.
        if non_restart_levels:
            self.run_all(cases, modes, non_restart_levels)

        # Handle restart: full + s1 should now be done; submit s2, compare.
        if do_restart:
            for case_name in cases:
                case = self.config[case_name]
                if not case.get("restart_test"):
                    continue
                valid_modes = [m for m in modes if m in case["modes"]]
                mode = valid_modes[0] if valid_modes else None
                if not mode:
                    continue
                self._complete_restart_from_manifest(
                    case_name, case, mode, manifest)

        # Clean up manifest.
        if os.path.isfile(manifest_path):
            os.remove(manifest_path)

    def _complete_restart_from_manifest(self, case_name, case, mode, manifest):
        """Complete a restart test whose full + s1 stages ran from manifest.

        Checks that full and s1 succeeded, copies restart file, submits s2,
        waits, and compares.
        """
        from .condor import _str_to_key

        key_full = (case_name, mode, "restart_full")
        key_s1 = (case_name, mode, "restart_s1")

        # Check that full and s1 are in run_results and succeeded.
        if key_full not in self.run_results:
            self.reporter.add(TestResult(
                case_name, mode, "restart", TestResult.SKIP,
                "full run not found in manifest"))
            return
        if key_s1 not in self.run_results:
            self.reporter.add(TestResult(
                case_name, mode, "restart", TestResult.SKIP,
                "stage 1 not found in manifest"))
            return

        work_full, ok_full = self.run_results[key_full]
        work_s1, ok_s1 = self.run_results[key_s1]

        if not ok_full:
            self.reporter.add(TestResult(
                case_name, mode, "restart", TestResult.FAIL,
                "full run failed"))
            return
        if not ok_s1:
            self.reporter.add(TestResult(
                case_name, mode, "restart", TestResult.FAIL,
                "stage 1 (first half) failed"))
            return

        try:
            self.builder.build_for_mode(mode)
        except RuntimeError:
            self.reporter.add(TestResult(
                case_name, mode, "restart", TestResult.FAIL, "build failed"))
            return

        exe = self.builder.get_exe(mode)
        nstep = case.get("regression_nstep") or case.get("short_nstep", 100)
        nstep_save = case.get("nstep_save", 10)
        use_mpi = mode in ("mpi", "mpi_omp")
        n_replicas = case.get("n_replicas", 1)
        prefix_s1 = f"{case['output_prefix']}_s1"

        # Find restart file from s1.
        if use_mpi and n_replicas > 1:
            rst_pattern = os.path.join(work_s1, f"{prefix_s1}_0*.rst")
            rst_files = sorted(glob.glob(rst_pattern))
            if not rst_files:
                self.reporter.add(TestResult(
                    case_name, mode, "restart", TestResult.FAIL,
                    "no restart files from stage 1"))
                return
            combined_rst = os.path.join(work_s1, f"{prefix_s1}.rst")
            with open(combined_rst, "wb") as out:
                for rf in rst_files:
                    with open(rf, "rb") as inp:
                        out.write(inp.read())
            rst_path = combined_rst
        else:
            rst_pattern = os.path.join(work_s1, f"{prefix_s1}*.rst")
            rst_files = sorted(glob.glob(rst_pattern))
            if not rst_files:
                self.reporter.add(TestResult(
                    case_name, mode, "restart", TestResult.FAIL,
                    "no restart files from stage 1"))
                return
            rst_path = rst_files[-1]

        # Stage s2: second half from restart.
        source_dir = os.path.join(self.repo_root, case["source_dir"])
        src_toml = os.path.join(source_dir, case["input_toml"])
        files_to_link = list(case.get("required_files", []))

        work_s2 = self._get_work_dir(case_name, mode, "restart_s2")
        ensure_dir(work_s2)
        setup_work_dir(work_s2, case["source_dir"],
                       files_to_link, self.repo_root)

        rst_basename = os.path.basename(rst_path)
        copy_file(rst_path, os.path.join(work_s2, rst_basename))

        prefix_s2 = f"{case['output_prefix']}_s2"
        s2_toml = os.path.join(work_s2, case["input_toml"])
        patch_toml(src_toml, s2_toml, {
            "prefix": f"./{prefix_s2}",
            "nstep": nstep,
            "nstep_save": nstep_save,
            "nstep_save_rst": nstep,
        })
        localize_file_paths(s2_toml, s2_toml)

        ok_s2 = self._condor_submit_and_wait(
            case, mode, exe, work_s2, case["input_toml"],
            extra_args=rst_basename, label=f"{case_name}_rst_s2")
        if not ok_s2:
            self.reporter.add(TestResult(
                case_name, mode, "restart", TestResult.FAIL,
                "stage 2 (restart) failed"))
            return

        # Compare full vs s2.
        self._compare_restart_outputs(case_name, case, mode,
                                      work_full, work_s2, prefix_s2)

    def _run_restart_condor(self, case_name, case, mode):
        """Restart consistency test using sequential HTCondor job submissions.

        Mirrors _run_restart() but replaces _execute() calls with
        _condor_submit_and_wait() so each stage runs on a cluster node.
        """
        try:
            self.builder.build_for_mode(mode)
        except RuntimeError:
            self.reporter.add(TestResult(
                case_name, mode, "restart", TestResult.FAIL, "build failed"
            ))
            return

        nstep = case.get("regression_nstep") or case.get("short_nstep", 100)
        nstep_save = case.get("nstep_save", 10)
        nstep_save_rst = case.get("nstep_save_rst", nstep // 2)
        half = nstep // 2
        half = (half // nstep_save) * nstep_save
        if half < nstep_save:
            half = nstep_save

        exe = self.builder.get_exe(mode)
        source_dir = os.path.join(self.repo_root, case["source_dir"])
        src_toml = os.path.join(source_dir, case["input_toml"])
        files_to_link = list(case.get("required_files", []))

        # --- Stage 1: full continuous run ---
        work_full = self._get_work_dir(case_name, mode, "restart_full")
        ensure_dir(work_full)
        setup_work_dir(work_full, case["source_dir"], files_to_link,
                       self.repo_root)
        dst_toml = os.path.join(work_full, case["input_toml"])
        patch_toml(src_toml, dst_toml, {
            "prefix": f"./{case['output_prefix']}",
            "nstep": nstep,
            "nstep_save": nstep_save,
            "nstep_save_rst": nstep_save_rst,
        })
        localize_file_paths(dst_toml, dst_toml)

        ok_full = self._condor_submit_and_wait(
            case, mode, exe, work_full, case["input_toml"],
            label=f"{case_name}_rst_full")
        if not ok_full:
            self.reporter.add(TestResult(
                case_name, mode, "restart", TestResult.FAIL,
                "full run failed"))
            return

        # --- Stage 2: first half ---
        work_s1 = self._get_work_dir(case_name, mode, "restart_s1")
        ensure_dir(work_s1)
        setup_work_dir(work_s1, case["source_dir"], files_to_link,
                       self.repo_root)
        prefix_s1 = f"{case['output_prefix']}_s1"
        s1_toml = os.path.join(work_s1, case["input_toml"])
        patch_toml(src_toml, s1_toml, {
            "prefix": f"./{prefix_s1}",
            "nstep": half,
            "nstep_save": nstep_save,
            "nstep_save_rst": half,
        })
        localize_file_paths(s1_toml, s1_toml)

        ok_s1 = self._condor_submit_and_wait(
            case, mode, exe, work_s1, case["input_toml"],
            label=f"{case_name}_rst_s1")
        if not ok_s1:
            self.reporter.add(TestResult(
                case_name, mode, "restart", TestResult.FAIL,
                "stage 1 (first half) failed"))
            return

        # --- Find restart file and run s2 ---
        use_mpi = mode in ("mpi", "mpi_omp")
        n_replicas = case.get("n_replicas", 1)

        if use_mpi and n_replicas > 1:
            rst_pattern = os.path.join(work_s1, f"{prefix_s1}_0*.rst")
            rst_files = sorted(glob.glob(rst_pattern))
            if not rst_files:
                self.reporter.add(TestResult(
                    case_name, mode, "restart", TestResult.FAIL,
                    "no restart files from stage 1"))
                return
            combined_rst = os.path.join(work_s1, f"{prefix_s1}.rst")
            with open(combined_rst, "wb") as out:
                for rf in rst_files:
                    with open(rf, "rb") as inp:
                        out.write(inp.read())
            rst_path = combined_rst
        else:
            rst_pattern = os.path.join(work_s1, f"{prefix_s1}*.rst")
            rst_files = sorted(glob.glob(rst_pattern))
            if not rst_files:
                self.reporter.add(TestResult(
                    case_name, mode, "restart", TestResult.FAIL,
                    "no restart files from stage 1"))
                return
            rst_path = rst_files[-1]

        work_s2 = self._get_work_dir(case_name, mode, "restart_s2")
        ensure_dir(work_s2)
        setup_work_dir(work_s2, case["source_dir"], files_to_link,
                       self.repo_root)

        rst_basename = os.path.basename(rst_path)
        copy_file(rst_path, os.path.join(work_s2, rst_basename))

        prefix_s2 = f"{case['output_prefix']}_s2"
        s2_toml = os.path.join(work_s2, case["input_toml"])
        patch_toml(src_toml, s2_toml, {
            "prefix": f"./{prefix_s2}",
            "nstep": nstep,
            "nstep_save": nstep_save,
            "nstep_save_rst": nstep,
        })
        localize_file_paths(s2_toml, s2_toml)

        ok_s2 = self._condor_submit_and_wait(
            case, mode, exe, work_s2, case["input_toml"],
            extra_args=rst_basename, label=f"{case_name}_rst_s2")
        if not ok_s2:
            self.reporter.add(TestResult(
                case_name, mode, "restart", TestResult.FAIL,
                "stage 2 (restart) failed"))
            return

        # --- Compare ---
        self._compare_restart_outputs(case_name, case, mode,
                                      work_full, work_s2, prefix_s2)

    def _compare_restart_outputs(self, case_name, case, mode,
                                 work_full, work_s2, prefix_s2):
        """Compare full run outputs against restart stage 2 outputs."""
        full_outs = sorted(glob.glob(
            os.path.join(work_full, f"{case['output_prefix']}*.out")))
        s2_outs = sorted(glob.glob(
            os.path.join(work_s2, f"{prefix_s2}*.out")))

        if not full_outs or not s2_outs:
            self.reporter.add(TestResult(
                case_name, mode, "restart", TestResult.FAIL,
                f"missing outputs: full={len(full_outs)}, s2={len(s2_outs)}"
            ))
            return

        all_passed = True
        worst_rdiff = 0.0

        for f_full, f_s2 in zip(full_outs, s2_outs):
            result = compare_out_files(
                f_s2, f_full, atol=1e-10, rtol=1e-10, align_steps=True)
            worst_rdiff = max(worst_rdiff, result.max_rel_diff)
            if not result.passed:
                all_passed = False

        if all_passed:
            self.reporter.add(TestResult(
                case_name, mode, "restart", TestResult.PASS,
                f"max_rdiff={worst_rdiff:.2e}"))
        else:
            self.reporter.add(TestResult(
                case_name, mode, "restart", TestResult.FAIL,
                result.summary))

    def generate_reference(self, cases, mode="serial", description=""):
        """Run tests and save output as reference data set.

        Creates a labeled subdirectory under test/reference/ with summary.json
        and per-case .out files.

        Parameters
        ----------
        cases : list[str]
            Test case names to generate reference for.
        mode : str
            Mode to use for reference runs.
        description : str
            Optional description for the reference set.
        """
        ref_base = os.path.join(self.repo_root, "test", "reference")
        label = generate_ref_label(self.repo_root)
        ref_set_dir = os.path.join(ref_base, label)
        ensure_dir(ref_set_dir)

        print(f"  Reference set: {label}")

        generated_cases = []

        for case_name in cases:
            case = self.config[case_name]

            if mode not in case["modes"]:
                mode_to_use = case["modes"][0]
            else:
                mode_to_use = mode

            missing = check_prerequisites(case, self.repo_root)
            if missing:
                print(f"  SKIP {case_name}: missing {', '.join(missing)}")
                continue

            try:
                self.builder.build_for_mode(mode_to_use)
            except RuntimeError as e:
                print(f"  FAIL {case_name}: build failed: {e}")
                continue

            nstep = case.get("regression_nstep")
            work_dir, success = self._ensure_run(
                case_name, case, mode_to_use,
                nstep_override=nstep, suffix="ref"
            )

            if not success:
                print(f"  FAIL {case_name}: run failed")
                continue

            out_files = self._find_out_files(work_dir, case)
            if not out_files:
                print(f"  FAIL {case_name}: no .out files")
                continue

            # Copy to reference set directory
            ref_dir = os.path.join(ref_set_dir, case_name)
            ensure_dir(ref_dir)
            for f in out_files:
                dst = os.path.join(ref_dir, os.path.basename(f))
                copy_file(f, dst)
                print(f"  Saved {dst}")

            generated_cases.append(case_name)
            print(f"  OK {case_name}: {len(out_files)} reference file(s)")

        if generated_cases:
            write_ref_summary(ref_set_dir, self.repo_root,
                              generated_cases, mode, description)
            print(f"  Summary written to {ref_set_dir}/summary.json")
        else:
            # Clean up empty directory
            shutil.rmtree(ref_set_dir, ignore_errors=True)
            print("  No reference data generated.")
