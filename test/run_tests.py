#!/usr/bin/env python3
"""UNISIS automated test runner.

Usage examples:
    # Smoke test everything
    ./test/run_tests.sh --levels run

    # Serial regression only
    ./test/run_tests.sh --levels regression --modes serial

    # Specific case + mode
    ./test/run_tests.sh --cases remd_t --levels run --modes mpi

    # Full consistency check
    ./test/run_tests.sh --levels consistency --modes serial,omp1,ompN

    # Restart validation
    ./test/run_tests.sh --cases md_simple --levels restart --modes serial

    # Generate reference data
    ./test/run_tests.sh --generate-reference --modes serial

    # Pre-built executables
    ./test/run_tests.sh --serial-exe ./build/sis --omp-exe ./build/sis
"""

import argparse
import os
import sys

# Add parent of this script to path so we can import framework and config
script_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, script_dir)

from config import TEST_CASES, ALL_MODES, ALL_LEVELS
from framework.utils import find_repo_root, ensure_dir, list_reference_sets
from framework.builder import BuildManager
from framework.reporter import TestReporter
from framework.runner import TestRunner


def parse_args():
    parser = argparse.ArgumentParser(
        description="UNISIS automated test runner",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )

    parser.add_argument(
        "--cases",
        default="all",
        help="Comma-separated test case names, or 'all' (default: all)",
    )
    parser.add_argument(
        "--modes",
        default="all",
        help="Comma-separated modes: serial,omp1,ompN,mpi,mpi_omp or 'all' (default: all)",
    )
    parser.add_argument(
        "--levels",
        default="run",
        help="Comma-separated levels: run,consistency,regression,restart,sampling or 'all' (default: run)",
    )

    # Executable paths
    parser.add_argument("--serial-exe", help="Path to serial-compiled sis")
    parser.add_argument("--omp-exe", help="Path to OpenMP-compiled sis")
    parser.add_argument("--mpi-exe", help="Path to MPI-compiled sis")

    # Build control
    build_group = parser.add_mutually_exclusive_group()
    build_group.add_argument(
        "--build", action="store_true", default=True,
        help="Build executables if not found (default)",
    )
    build_group.add_argument(
        "--no-build", action="store_true",
        help="Do not build; fail if executables missing",
    )

    # Parallelization settings
    parser.add_argument(
        "--omp-threads", type=int, default=4,
        help="Thread count for ompN/mpi_omp modes (default: 4)",
    )
    parser.add_argument(
        "--mpi-cmd", default="mpirun",
        help="MPI launcher command (default: mpirun)",
    )

    # Execution settings
    parser.add_argument(
        "--timeout", type=int, default=300,
        help="Per-test timeout in seconds (default: 300)",
    )
    parser.add_argument(
        "--work-dir",
        help="Working directory for test runs (default: test/_work)",
    )
    parser.add_argument(
        "--keep-work", action="store_true",
        help="Preserve working directories after tests",
    )
    parser.add_argument(
        "--verbose", "-v", action="store_true",
        help="Print stdout/stderr from test runs",
    )

    # HTCondor dispatch
    condor_group = parser.add_argument_group("HTCondor dispatch")
    condor_group.add_argument(
        "--condor", action="store_true",
        help="Submit jobs via HTCondor instead of running locally",
    )
    condor_group.add_argument(
        "--condor-memory", type=int, default=4,
        help="Memory per job in GB (default: 4)",
    )
    condor_group.add_argument(
        "--condor-poll-interval", type=int, default=30,
        help="Polling interval in seconds while waiting for jobs (default: 30)",
    )
    condor_group.add_argument(
        "--condor-timeout", type=int, default=3600,
        help="Maximum seconds to wait for all condor jobs (default: 3600)",
    )
    condor_group.add_argument(
        "--condor-submit-only", action="store_true",
        help="Submit jobs via HTCondor and exit (collect later with --condor)",
    )

    # Reference generation and selection
    parser.add_argument(
        "--generate-reference", action="store_true",
        help="Run tests and save output as a new reference data set",
    )
    parser.add_argument(
        "--reference",
        help="Reference data set label for regression tests (default: latest)",
    )
    parser.add_argument(
        "--ref-description", default="",
        help="Description to include in the reference set summary",
    )

    return parser.parse_args()


def resolve_list(value, all_values):
    """Resolve a comma-separated string or 'all' to a list."""
    if value == "all":
        return list(all_values)
    items = [x.strip() for x in value.split(",")]
    for item in items:
        if item not in all_values:
            print(f"Error: unknown value '{item}'. Valid: {', '.join(all_values)}")
            sys.exit(2)
    return items


def print_info():
    """Print test framework info and usage when no arguments given."""
    print("UNISIS Automated Test Framework")
    print("=" * 50)
    print()
    print(f"Test cases ({len(TEST_CASES)}):")
    for name, case in TEST_CASES.items():
        modes_str = ", ".join(case["modes"])
        print(f"  {name:<18s} {case['description']}")
        print(f"  {'':<18s} modes: {modes_str}")
    print()
    print(f"Parallelization modes: {', '.join(ALL_MODES)}")
    print(f"Test levels:           {', '.join(ALL_LEVELS)}")

    # Show available reference sets
    ref_base = os.path.join(script_dir, "reference")
    ref_sets = list_reference_sets(ref_base)
    if ref_sets:
        print()
        print(f"Reference data sets ({len(ref_sets)}):")
        for label, summary in ref_sets:
            commit = summary.get("git_commit", "?")
            date = summary.get("date", "?")
            cases = summary.get("cases", [])
            desc = summary.get("description", "")
            info = f"commit={commit}, {len(cases)} case(s)"
            if desc:
                info += f", {desc}"
            print(f"  {label:<30s} {info}")
        print(f"  (latest: {ref_sets[-1][0]})")
    else:
        print()
        print("Reference data sets: none (run --generate-reference to create)")

    print()
    print("Usage examples:")
    print("  ./test/run_tests.py --levels run --modes serial --serial-exe ./build/sis")
    print("  ./test/run_tests.py --levels run,consistency --modes serial,omp1,ompN --omp-exe ./build/sis")
    print("  ./test/run_tests.py --levels regression --cases md_simple --serial-exe ./build/sis")
    print("  ./test/run_tests.py --levels restart --cases md_simple --serial-exe ./build/sis")
    print("  ./test/run_tests.py --generate-reference --serial-exe ./build/sis")
    print("  ./test/run_tests.py --generate-reference --ref-description 'after bugfix #42' --serial-exe ./build/sis")
    print("  ./test/run_tests.py --levels regression --reference 2026-03-28_abc1234 --serial-exe ./build/sis")
    print("  ./test/run_tests.py --levels all --modes all --build")
    print()
    print("Run with --help for full option list.")


def main():
    if len(sys.argv) == 1:
        print_info()
        return

    args = parse_args()
    repo_root = find_repo_root(script_dir)

    # Resolve selections
    cases = resolve_list(args.cases, TEST_CASES.keys())
    modes = resolve_list(args.modes, ALL_MODES)
    levels = resolve_list(args.levels, ALL_LEVELS)

    # HTCondor sanity check
    if args.condor or args.condor_submit_only:
        import shutil as _shutil
        if not _shutil.which("condor_submit"):
            print("Error: condor_submit not found in PATH. Is HTCondor installed?")
            sys.exit(2)

    # Work directory
    work_base = args.work_dir or os.path.join(script_dir, "_work")
    ensure_dir(work_base)

    # Build manager
    builder = BuildManager(repo_root, verbose=args.verbose)

    # Set explicit executables if provided
    if args.serial_exe:
        builder.set_exe("serial", args.serial_exe)
    if args.omp_exe:
        builder.set_exe("omp", args.omp_exe)
    if args.mpi_exe:
        builder.set_exe("mpi", args.mpi_exe)

    # Reporter
    reporter = TestReporter(verbose=args.verbose)

    # Runner
    runner = TestRunner(
        config=TEST_CASES,
        builder=builder,
        reporter=reporter,
        repo_root=repo_root,
        work_base=work_base,
        omp_threads=args.omp_threads,
        mpi_cmd=args.mpi_cmd,
        timeout=args.timeout,
        keep_work=args.keep_work,
        verbose=args.verbose,
        ref_label=args.reference,
        condor=args.condor or args.condor_submit_only,
        condor_memory_gb=args.condor_memory,
        condor_poll_interval=args.condor_poll_interval,
        condor_timeout=args.condor_timeout,
    )

    if args.generate_reference:
        print("Generating reference data...")
        mode = modes[0] if modes else "serial"
        runner.generate_reference(cases, mode=mode,
                                  description=args.ref_description)
        print("Done.")
        return

    # Pre-build needed variants
    if not args.no_build:
        needed_modes = set()
        for case_name in cases:
            case = TEST_CASES[case_name]
            for m in modes:
                if m in case["modes"]:
                    needed_modes.add(m)
        try:
            builder.build_all_needed(list(needed_modes))
        except RuntimeError as e:
            print(f"\nBuild error: {e}")
            if args.no_build:
                sys.exit(2)

    # Run tests
    if args.condor_submit_only:
        runner.condor_submit_only(cases, modes, levels)
        return
    elif args.condor:
        runner.condor_run_all(cases, modes, levels)
    else:
        runner.run_all(cases, modes, levels)

    # Summary
    exit_code = reporter.print_summary()
    sys.exit(exit_code)


if __name__ == "__main__":
    main()
