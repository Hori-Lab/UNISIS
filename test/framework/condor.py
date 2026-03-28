"""HTCondor job submission interface for the UNISIS test framework."""

import datetime
import json
import os
import re
import subprocess
import time


# ---------------------------------------------------------------------------
# CPU count helper
# ---------------------------------------------------------------------------

def _cpu_count(mode, omp_threads, mpi_ranks):
    """Return the number of CPUs to request for a given parallelization mode."""
    if mode in ("serial", "omp1"):
        return 1
    if mode == "ompN":
        return omp_threads
    if mode == "mpi":
        return mpi_ranks
    if mode == "mpi_omp":
        return mpi_ranks * omp_threads
    raise ValueError(f"Unknown mode: {mode}")


# ---------------------------------------------------------------------------
# File writers
# ---------------------------------------------------------------------------

def write_job_script(work_dir, cmd, env_vars):
    """Write _condor_job.sh into work_dir and return its absolute path.

    Parameters
    ----------
    work_dir : str
        Directory to write the script into.
    cmd : str
        Shell command to execute (will be passed to exec).
    env_vars : dict
        Environment variables to export before running cmd.

    Returns
    -------
    str
        Absolute path to the generated script.
    """
    script_path = os.path.join(work_dir, "_condor_job.sh")
    lines = ["#!/bin/sh"]
    for key, val in env_vars.items():
        lines.append(f"export {key}={val}")
    lines.append(f"exec {cmd}")
    lines.append("")
    with open(script_path, "w") as f:
        f.write("\n".join(lines))
    os.chmod(script_path, 0o755)
    return script_path


def write_submit_file(work_dir, case_name, mode, n_cpus, memory_gb=4):
    """Write _condor_submit into work_dir and return its absolute path.

    Parameters
    ----------
    work_dir : str
        Directory to write the submit file into (also used as initialdir).
    case_name : str
        Test case name (used in batch_name).
    mode : str
        Parallelization mode (used in batch_name).
    n_cpus : int
        Number of CPUs to request.
    memory_gb : int
        Memory to request in GB.

    Returns
    -------
    str
        Absolute path to the generated submit file.
    """
    submit_path = os.path.join(work_dir, "_condor_submit")
    job_script = os.path.join(work_dir, "_condor_job.sh")
    content = (
        f"universe        = vanilla\n"
        f"batch_name      = unisis_{case_name}_{mode}\n"
        f"executable      = {job_script}\n"
        f"getenv          = True\n"
        f"initialdir      = {work_dir}\n"
        f"output          = _condor_stdout\n"
        f"error           = _condor_stderr\n"
        f"log             = _condor.log\n"
        f"request_memory  = {memory_gb} GB\n"
        f"Request_CPUs    = {n_cpus}\n"
        f"queue\n"
    )
    with open(submit_path, "w") as f:
        f.write(content)
    return submit_path


# ---------------------------------------------------------------------------
# Submit
# ---------------------------------------------------------------------------

def submit_job(submit_file_path):
    """Run condor_submit and return the cluster ID as a string.

    Parameters
    ----------
    submit_file_path : str
        Path to the HTCondor submit file.

    Returns
    -------
    str
        The cluster ID (e.g. "12345").

    Raises
    ------
    RuntimeError
        If condor_submit returns non-zero or the cluster ID cannot be parsed.
    """
    result = subprocess.run(
        ["condor_submit", submit_file_path],
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        raise RuntimeError(
            f"condor_submit failed (rc={result.returncode}):\n{result.stderr.strip()}"
        )
    m = re.search(r"submitted to cluster (\d+)", result.stdout)
    if not m:
        raise RuntimeError(
            f"Could not parse cluster ID from condor_submit output:\n{result.stdout.strip()}"
        )
    return m.group(1)


# ---------------------------------------------------------------------------
# Log parsing
# ---------------------------------------------------------------------------

def parse_log_for_exit_code(log_path):
    """Parse an HTCondor event log to find the job exit code.

    Scans for event-005 (Job terminated) and extracts the return value.

    Parameters
    ----------
    log_path : str
        Path to the _condor.log file written by HTCondor.

    Returns
    -------
    int or None
        Exit code if the job has terminated, None if still running.
        Returns -1 for abnormal termination (signal).
    """
    if not os.path.isfile(log_path):
        return None
    with open(log_path) as f:
        lines = f.readlines()
    in_term = False
    for line in lines:
        if line.startswith("005 "):
            in_term = True
            continue
        if in_term:
            m = re.search(r"Normal termination \(return value (\d+)\)", line)
            if m:
                return int(m.group(1))
            if "Abnormal termination" in line:
                return -1
    return None


# ---------------------------------------------------------------------------
# CondorJob
# ---------------------------------------------------------------------------

class CondorJob:
    """Represents a single HTCondor job submission."""

    def __init__(self, cluster_id, work_dir, log_path, key):
        """
        Parameters
        ----------
        cluster_id : str
            HTCondor cluster ID (e.g. "12345").
        work_dir : str
            Absolute path to the job's working directory.
        log_path : str
            Path to the _condor.log file for this job.
        key : tuple
            Run key (case_name, mode, suffix) identifying this run.
        """
        self.cluster_id = cluster_id
        self.work_dir = work_dir
        self.log_path = log_path
        self.key = key
        self.exit_code = None
        self.done = False


# ---------------------------------------------------------------------------
# CondorBatch
# ---------------------------------------------------------------------------

class CondorBatch:
    """Tracks a collection of submitted HTCondor jobs and waits for all to complete."""

    def __init__(self, poll_interval=30, timeout=3600, verbose=False):
        """
        Parameters
        ----------
        poll_interval : int
            Seconds between log-file polls.
        timeout : int
            Maximum total seconds to wait before giving up.
        verbose : bool
            Print progress messages.
        """
        self.poll_interval = poll_interval
        self.timeout = timeout
        self.verbose = verbose
        self.jobs = []

    def add(self, job):
        """Add a CondorJob to track."""
        self.jobs.append(job)

    def all_done(self):
        """Return True when every tracked job has completed."""
        return all(j.done for j in self.jobs)

    def _poll_once(self):
        """Check log files for all pending jobs and update their status."""
        for job in self.jobs:
            if job.done:
                continue
            rc = parse_log_for_exit_code(job.log_path)
            if rc is not None:
                job.exit_code = rc
                job.done = True

    def wait_all(self):
        """Block until all tracked jobs complete or timeout expires.

        Returns
        -------
        dict
            Mapping of job.key -> exit_code (int). Jobs that timed out
            receive exit_code -1.
        """
        t0 = time.monotonic()
        self._poll_once()
        while not self.all_done():
            elapsed = time.monotonic() - t0
            if elapsed >= self.timeout:
                pending = [j for j in self.jobs if not j.done]
                for job in pending:
                    job.exit_code = -1
                    job.done = True
                    if self.verbose:
                        print(f"  TIMEOUT: job {job.cluster_id} key={job.key}")
                break
            pending_count = sum(1 for j in self.jobs if not j.done)
            if self.verbose:
                print(f"  Waiting: {pending_count} job(s) pending "
                      f"(elapsed {int(time.monotonic() - t0)}s)...")
            time.sleep(self.poll_interval)
            self._poll_once()
        return {job.key: job.exit_code for job in self.jobs}


# ---------------------------------------------------------------------------
# Convenience entry point
# ---------------------------------------------------------------------------

def prepare_and_submit(work_dir, case_name, mode, cmd, env_vars,
                       omp_threads, mpi_ranks, memory_gb, key, verbose=False):
    """Write job files, submit to HTCondor, and return a CondorJob.

    Parameters
    ----------
    work_dir : str
        Working directory for the job (must already exist).
    case_name : str
        Test case name.
    mode : str
        Parallelization mode.
    cmd : str
        Shell command to execute inside the job.
    env_vars : dict
        Environment variables to set in the job script.
    omp_threads : int
        Number of OpenMP threads (used for CPU count calculation).
    mpi_ranks : int
        Number of MPI ranks (used for CPU count calculation).
    memory_gb : int
        Memory to request in GB.
    key : tuple
        Run key (case_name, mode, suffix) for result tracking.
    verbose : bool
        Print submission message.

    Returns
    -------
    CondorJob
    """
    n_cpus = _cpu_count(mode, omp_threads, mpi_ranks)
    write_job_script(work_dir, cmd, env_vars)
    submit_file = write_submit_file(work_dir, case_name, mode, n_cpus, memory_gb)
    cluster_id = submit_job(submit_file)
    log_path = os.path.join(work_dir, "_condor.log")
    if verbose:
        print(f"  Submitted {case_name}/{mode} -> cluster {cluster_id} "
              f"({n_cpus} CPU(s), {memory_gb} GB)")
    return CondorJob(cluster_id, work_dir, log_path, key)


# ---------------------------------------------------------------------------
# Manifest (for submit-only / collect workflow)
# ---------------------------------------------------------------------------

MANIFEST_FILENAME = "_condor_manifest.json"


def _key_to_str(key):
    """Convert a run key tuple to a manifest-safe string."""
    return "|".join(str(k) for k in key)


def _str_to_key(key_str):
    """Convert a manifest key string back to a tuple."""
    return tuple(key_str.split("|"))


def write_manifest(manifest_path, cases, modes, levels, jobs, restart_jobs=None):
    """Write a condor manifest file.

    Parameters
    ----------
    manifest_path : str
        Path to write the JSON manifest.
    cases : list[str]
        Test case names that were submitted.
    modes : list[str]
        Parallelization modes that were submitted.
    levels : list[str]
        Test levels that were submitted.
    jobs : dict
        Mapping of key_str -> {work_dir, cluster_id, log_path}.
    restart_jobs : dict, optional
        Same format as jobs, for restart stages (full, s1).
    """
    data = {
        "submitted_at": datetime.datetime.now().isoformat(),
        "cases": list(cases),
        "modes": list(modes),
        "levels": list(levels),
        "jobs": jobs,
        "restart_jobs": restart_jobs or {},
    }
    with open(manifest_path, "w") as f:
        json.dump(data, f, indent=2)


def read_manifest(manifest_path):
    """Read a condor manifest file.

    Returns
    -------
    dict or None
        The manifest data, or None if the file does not exist.
    """
    if not os.path.isfile(manifest_path):
        return None
    with open(manifest_path) as f:
        return json.load(f)


def collect_results_from_manifest(manifest):
    """Check completion status of all jobs in a manifest.

    Parameters
    ----------
    manifest : dict
        Manifest data from read_manifest().

    Returns
    -------
    dict
        With keys "completed" (key_str -> exit_code) and
        "pending" (key_str -> {work_dir, cluster_id, log_path}).
    """
    completed = {}
    pending = {}
    for section in ("jobs", "restart_jobs"):
        for key_str, info in manifest.get(section, {}).items():
            rc = parse_log_for_exit_code(info["log_path"])
            if rc is not None:
                completed[key_str] = rc
            else:
                pending[key_str] = info
    return {"completed": completed, "pending": pending}
