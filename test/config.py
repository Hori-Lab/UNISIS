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
short_nstep_save : int or None, optional
    Output frequency for smoke tests. None = use original.
regression_nstep : int or None
    nstep for regression tests. None = use original.
regression_nstep_save : int or None, optional
    Output frequency for regression tests. None = use original.
tolerances : dict
    {"atol": float, "rtol": float} for regression comparison.
sampling_nstep : int, optional
    nstep for sampling tests (longer runs).
reference_mode : str, optional
    Mode to use for generating reference data (default: first in modes).
reference_extensions : list[str], optional
    Additional file extensions to copy to reference data (e.g. [".dcd", ".bp"]).
    Files matching {output_prefix}*.{ext} are copied. .out is always included.
sampling_properties : list[tuple], optional
    [(column_index, expected_mean, tolerance_sigma), ...] for sampling checks.
    Populated after reference values are established.
"""

TEST_CASES = {

    "md_simple": {
        "description": "Basic single-chain MD without PBC",
        "source_dir": "test/data",
        "input_toml": "T2HP.toml",
        "required_files": [
            "unisis.ff",
            "T2HP.fasta",
            "T2HP.db",
            "T2HP.xyz",
        ],
        "modes": ["serial", "omp1", "ompN", "mpi", "mpi_omp"],
        "reference_mode": "serial",
        "reference_extensions": [".dcd", ".bp"],
        "mpi_ranks": 4,
        "output_prefix": "test",
        "has_replica_cols": False,
        "restart_test": True,
        "short_nstep": 100,
        "short_nstep_save": 10,
        "regression_nstep": 2000,
        "regression_nstep_save": 100,
        "nstep_save_rst": 500,
        "tolerances": {"atol": 1e-8, "rtol": 1e-6},
        "sampling_nstep": 100000,
        "sampling_properties": [],  # to be populated
        "gradient_test": True,
    },

}

# All valid parallelization modes
ALL_MODES = ["serial", "omp1", "ompN", "mpi", "mpi_omp"]

# All valid test levels
ALL_LEVELS = ["run", "consistency", "regression", "restart", "sampling", "gradient"]
