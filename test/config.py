"""Test case definitions for the UNISIS test framework.

Each test case is a dict specifying inputs, valid parallelization modes,
output expectations, and tolerance settings.

Keys
----
description : str
    Human-readable description.
source_dir : str
    Directory containing input files (relative to repo root).
input_toml : str
    Input TOML file name (relative to source_dir).
required_files : list[str]
    Files needed (relative to source_dir). Checked before running.
modes : list[str]
    Valid parallelization modes: serial, omp1, ompN, mpi, mpi_omp.
mpi_ranks : int
    Number of MPI ranks (0 = not MPI).
n_replicas : int, optional
    Number of replicas (for REMD tests).
output_prefix : str
    Prefix of output files.
has_replica_cols : bool
    Whether .out files have replica label columns.
restart_test : bool
    Whether to include in restart validation tests.
restart_input_tomls : list[str], optional
    TOML files for restart stages (stage1, stage2).
short_nstep : int or None
    Reduced nstep for smoke tests. None = use original.
regression_nstep : int or None
    nstep for regression tests. None = use original.
tolerances : dict
    {"atol": float, "rtol": float} for regression comparison.
sampling_nstep : int, optional
    nstep for sampling tests (longer runs).
sampling_properties : list[tuple], optional
    [(column_index, expected_mean, tolerance_sigma), ...] for sampling checks.
    Populated after reference values are established.
"""

TEST_CASES = {

    "md_simple": {
        "description": "Basic single-chain MD without PBC",
        "source_dir": "test/data",
        "input_toml": "input_md_simple.toml",
        "required_files": [
            "T2HP.pdb",
            "unisis.ff",
            "T2HP.fasta",
        ],
        "modes": ["serial", "omp1", "ompN"],
        "mpi_ranks": 0,
        "output_prefix": "test_md",
        "has_replica_cols": False,
        "restart_test": True,
        "short_nstep": 50,
        "regression_nstep": 100,
        "nstep_save": 10,
        "nstep_save_rst": 50,
        "tolerances": {"atol": 1e-8, "rtol": 1e-6},
        "sampling_nstep": 100000,
        "sampling_properties": [],  # to be populated
    },

}

# All valid parallelization modes
ALL_MODES = ["serial", "omp1", "ompN", "mpi", "mpi_omp"]

# All valid test levels
ALL_LEVELS = ["run", "consistency", "regression", "restart", "sampling"]
