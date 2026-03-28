"""Utility functions for the test framework."""

import datetime
import json
import os
import platform
import subprocess
import shutil
import time


def run_command(cmd, cwd=None, env=None, timeout=300, capture=True):
    """Run a shell command and return (returncode, stdout, stderr, elapsed).

    Parameters
    ----------
    cmd : list[str] or str
        Command to run. If str, uses shell=True.
    cwd : str, optional
        Working directory.
    env : dict, optional
        Extra environment variables (merged with os.environ).
    timeout : int
        Timeout in seconds.
    capture : bool
        If True, capture stdout/stderr. If False, let them pass through.

    Returns
    -------
    tuple of (int, str, str, float)
        (returncode, stdout, stderr, elapsed_seconds)
    """
    full_env = os.environ.copy()
    if env:
        full_env.update(env)

    use_shell = isinstance(cmd, str)

    t0 = time.monotonic()
    try:
        result = subprocess.run(
            cmd,
            cwd=cwd,
            env=full_env,
            shell=use_shell,
            timeout=timeout,
            capture_output=capture,
            text=True if capture else None,
        )
        elapsed = time.monotonic() - t0
        stdout = result.stdout if capture else ""
        stderr = result.stderr if capture else ""
        return result.returncode, stdout, stderr, elapsed
    except subprocess.TimeoutExpired:
        elapsed = time.monotonic() - t0
        return -1, "", f"TIMEOUT after {timeout}s", elapsed


def ensure_dir(path):
    """Create directory if it doesn't exist."""
    os.makedirs(path, exist_ok=True)


def symlink_force(src, dst):
    """Create a symlink using a relative path, removing existing destination."""
    if os.path.islink(dst) or os.path.exists(dst):
        os.remove(dst)
    # Compute relative path from the symlink's directory to the source
    rel = os.path.relpath(os.path.abspath(src), os.path.dirname(os.path.abspath(dst)))
    os.symlink(rel, dst)


def setup_work_dir(work_dir, source_dir, input_files, repo_root):
    """Create a work directory and symlink input files into it.

    Parameters
    ----------
    work_dir : str
        Path to working directory (created if needed).
    source_dir : str
        Source directory containing test input files (relative to repo_root).
    input_files : list[str]
        Files to symlink (relative to source_dir).
    repo_root : str
        Repository root path.
    """
    ensure_dir(work_dir)

    abs_source = os.path.join(repo_root, source_dir)

    for f in input_files:
        src = os.path.normpath(os.path.join(abs_source, f))
        dst = os.path.join(work_dir, os.path.basename(f))
        if os.path.exists(src):
            symlink_force(src, dst)


def copy_file(src, dst):
    """Copy a file, overwriting destination."""
    shutil.copy2(src, dst)


def find_repo_root(start=None):
    """Walk up from start (default: this file's location) to find repo root.

    Looks for CMakeLists.txt + src/ + check/ as indicators.
    """
    if start is None:
        start = os.path.dirname(os.path.abspath(__file__))
    path = start
    for _ in range(10):
        if (os.path.isfile(os.path.join(path, "CMakeLists.txt"))
                and os.path.isdir(os.path.join(path, "src"))
                and os.path.isdir(os.path.join(path, "check"))):
            return path
        parent = os.path.dirname(path)
        if parent == path:
            break
        path = parent
    raise RuntimeError(f"Cannot find repo root from {start}")


def check_prerequisites(case_config, repo_root):
    """Check that all required files for a test case exist.

    Returns list of missing files (empty if all present).
    """
    source_dir = os.path.join(repo_root, case_config["source_dir"])
    missing = []

    # Check input TOML
    toml_path = os.path.join(source_dir, case_config["input_toml"])
    if not os.path.exists(toml_path):
        missing.append(case_config["input_toml"])

    # Check required files
    for f in case_config.get("required_files", []):
        fpath = os.path.normpath(os.path.join(source_dir, f))
        if not os.path.exists(fpath):
            missing.append(f)

    # Check restart input TOMLs if relevant
    for f in case_config.get("restart_input_tomls", []):
        fpath = os.path.join(source_dir, f)
        if not os.path.exists(fpath):
            missing.append(f)

    return missing


def generate_mpi_wrapper(work_dir, exe_path, input_toml, extra_args=""):
    """Generate an MPI wrapper script that auto-detects rank.

    Returns path to the generated script.
    """
    script_path = os.path.join(work_dir, "mpi_wrapper.sh")
    args = f"{input_toml}"
    if extra_args:
        args += f" {extra_args}"

    content = f"""#!/bin/sh
if [ -n "$OMPI_COMM_WORLD_RANK" ]; then rank=$OMPI_COMM_WORLD_RANK;
elif [ -n "$PMI_RANK" ]; then rank=$PMI_RANK;
elif [ -n "$SLURM_PROCID" ]; then rank=$SLURM_PROCID;
else rank=0; fi
suffix=$(printf "%3.3d" $rank)
{exe_path} {args} 1> ./out.$suffix 2> ./err.$suffix
"""
    with open(script_path, "w") as f:
        f.write(content)
    os.chmod(script_path, 0o755)
    return script_path


def get_git_commit(repo_root):
    """Get the current git commit hash, or 'unknown' if unavailable."""
    try:
        result = subprocess.run(
            ["git", "rev-parse", "--short", "HEAD"],
            cwd=repo_root, capture_output=True, text=True, timeout=5
        )
        if result.returncode == 0:
            return result.stdout.strip()
    except (subprocess.TimeoutExpired, FileNotFoundError):
        pass
    return "unknown"


def generate_ref_label(repo_root):
    """Generate a label for a new reference data set: YYYY-MM-DD_<commit>."""
    date_str = datetime.date.today().isoformat()
    commit = get_git_commit(repo_root)
    return f"{date_str}_{commit}"


def write_ref_summary(ref_set_dir, repo_root, cases, mode, description=""):
    """Write summary.json for a reference data set."""
    summary = {
        "date": datetime.date.today().isoformat(),
        "git_commit": get_git_commit(repo_root),
        "machine": platform.node(),
        "platform": platform.platform(),
        "mode": mode,
        "cases": cases,
        "description": description,
    }
    path = os.path.join(ref_set_dir, "summary.json")
    with open(path, "w") as f:
        json.dump(summary, f, indent=2)
    return path


def list_reference_sets(ref_base):
    """List available reference data sets, sorted by name (newest last).

    Returns list of (label, summary_dict) tuples.
    """
    sets = []
    if not os.path.isdir(ref_base):
        return sets
    for entry in os.listdir(ref_base):
        summary_path = os.path.join(ref_base, entry, "summary.json")
        if os.path.isfile(summary_path):
            with open(summary_path) as f:
                summary = json.load(f)
            sets.append((entry, summary))
    # Sort by date in summary.json (YYYY-MM-DD), with label as tiebreaker
    sets.sort(key=lambda x: (x[1].get("date", ""), x[0]))
    return sets


def resolve_reference_set(ref_base, label=None):
    """Resolve a reference set label to its directory path.

    If label is None, returns the latest (last alphabetically).
    Returns (label, path) or (None, None) if not found.
    """
    if label:
        path = os.path.join(ref_base, label)
        if os.path.isdir(path):
            return label, path
        return None, None

    sets = list_reference_sets(ref_base)
    if sets:
        latest_label = sets[-1][0]
        return latest_label, os.path.join(ref_base, latest_label)
    return None, None
