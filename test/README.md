# UNISIS Test Framework

Automated test framework for verifying UNISIS across different physics options, parallelization modes, and code versions. Pure Python (standard library only) -- no external dependencies needed.

## Quick Start

From the repository root:

```bash
# Show available test cases and usage
./test/run_tests.py

# Smoke test: check that the code compiles and runs
./test/run_tests.py --levels run --modes serial --serial-exe ./build/sis

# Regression test against stored reference data
./test/run_tests.py --levels regression --modes serial --serial-exe ./build/sis

# Multiple levels and modes at once
./test/run_tests.py --levels run,consistency,regression --modes serial,omp1,ompN \
    --serial-exe ./build/sis --omp-exe ./build/sis
```

## Concepts

The framework has three orthogonal dimensions:

### Test Cases

Each test case exercises a different combination of physics options, sampling algorithm, or RNA type. Defined in `test/config.py`.

| Case | Description | Modes |
|------|-------------|-------|
| `md_simple` | Basic MD without PBC (T2HP single chain) | serial, omp1, ompN |

### Parallelization Modes

| Mode | Build | Environment | For testing |
|------|-------|-------------|-------------|
| `serial` | Without OpenMP | `OMP_NUM_THREADS=1` | Baseline, no threading |
| `omp1` | With OpenMP | `OMP_NUM_THREADS=1` | OpenMP overhead, single thread |
| `ompN` | With OpenMP | `OMP_NUM_THREADS=4` | Multi-threaded correctness |
| `mpi` | With MPI | `OMP_NUM_THREADS=1` | MPI communication |
| `mpi_omp` | With MPI | `OMP_NUM_THREADS=4` | Hybrid parallelization |

Three separate builds are required:

| Build | Directory | CMake Flags |
|-------|-----------|-------------|
| serial | `build_serial/` | `-DCMAKE_DISABLE_FIND_PACKAGE_OpenMP=TRUE` |
| omp | `build/` | (default) |
| mpi | `build_mpi/` | `-DBUILD_MPI=ON` |

You can either let the framework build them (`--build`, the default) or provide pre-built executables via `--serial-exe`, `--omp-exe`, `--mpi-exe`.

### Test Levels

Tests are organized into progressive levels. You choose which levels to run based on what you are checking.

| Level | Purpose | When to use |
|-------|---------|-------------|
| `run` | Smoke test: exits 0, output files exist | After any code change |
| `consistency` | Same case across modes gives matching results | After changing parallelization or force accumulation |
| `regression` | Output matches stored reference data | After any code change that could affect numerics |
| `restart` | Split run from checkpoint matches continuous run | After changing restart I/O or integrator state |
| `sampling` | Thermodynamic averages match expected values | After changing force field or sampling algorithm |

Typical developer workflow:

1. **Quick check while coding**: `--levels run --modes serial`
2. **Before committing**: `--levels run,regression --modes serial`
3. **Before merging**: `--levels run,consistency,regression,restart`
4. **After force field changes**: `--levels sampling`

## Command-Line Reference

```
./test/run_tests.py [OPTIONS]
```

With no arguments, prints info and usage examples.

| Option | Default | Description |
|--------|---------|-------------|
| `--cases NAMES` | `all` | Comma-separated case names or `all` |
| `--modes MODES` | `all` | Comma-separated modes or `all` |
| `--levels LEVELS` | `run` | Comma-separated levels or `all` |
| `--serial-exe PATH` | | Path to serial-compiled `sis` |
| `--omp-exe PATH` | | Path to OpenMP-compiled `sis` |
| `--mpi-exe PATH` | | Path to MPI-compiled `sis` |
| `--build` | yes | Build executables if not found |
| `--no-build` | | Do not build; fail if missing |
| `--omp-threads N` | 4 | Thread count for ompN/mpi_omp |
| `--mpi-cmd CMD` | `mpirun` | MPI launcher command |
| `--timeout SECS` | 300 | Per-test timeout |
| `--work-dir DIR` | `test/_work` | Working directory for test runs |
| `--keep-work` | | Preserve work directories after tests |
| `--verbose`, `-v` | | Print stdout/stderr from runs |
| `--generate-reference` | | Run and save output as a new reference data set |
| `--reference LABEL` | latest | Reference data set to compare against |
| `--ref-description TEXT` | | Description to include in reference set metadata |

## Managing Reference Data

Reference data for regression tests is organized under `test/reference/`, each containing a `summary.json` and per-case subdirectories with `.out` files:

```
test/reference/
  2026-01-01_abc1234/             # Auto-generated label: YYYY-MM-DD_<commit>
    summary.json
    md_simple/
      test_md.out
    nick_restraint/
      md.out
```

By default, regression tests compare against the **latest** reference set according to the date in `summary.json`. You can select a specific set with `--reference <label>`.

### Generating reference data

After confirming the code produces correct results (e.g., by independent validation or comparison with a previous trusted version):

```bash
./test/run_tests.py --generate-reference --modes serial --serial-exe ./build/sis
```

This creates a new reference set labeled `YYYY-MM-DD_<commit>` (e.g., `2026-01-01_abc1234`) with a `summary.json` recording the date, git commit, machine name, and platform.

To add a description:

```bash
./test/run_tests.py --generate-reference --ref-description "after fixing electrostatics bug" --serial-exe ./build/sis
```

To regenerate for a specific case only:

```bash
./test/run_tests.py --generate-reference --cases md_simple --modes serial --serial-exe ./build/sis
```

### Selecting a reference set for regression

```bash
# Use the latest reference set (default)
./test/run_tests.py --levels regression --serial-exe ./build/sis

# Use a specific reference set
./test/run_tests.py --levels regression --reference initial --serial-exe ./build/sis
./test/run_tests.py --levels regression --reference 2026-01-01_abc1234 --serial-exe ./build/sis
```

Run with no arguments to see all available reference sets.

### When to generate a new reference set

Generate a new reference data set when:

- **Intentional numerical changes**: A bug fix, force field correction, or algorithm change that is expected to alter output values. The old reference would cause false regression failures.
- **New test case added**: The new case has no reference data yet.
- **Different machine architecture**: You want to track expected results on a different platform (e.g., ARM vs x86).

Do **not** generate reference data to make a failing test pass without understanding why it changed. Investigate first.

Old reference sets are preserved so you can compare against any historical baseline.

### Procedure

1. Run the regression test to confirm it fails and understand the magnitude of change.
2. Verify the new results are correct (analytical check, comparison with independent code, etc.).
3. Generate a new reference set: `./test/run_tests.py --generate-reference --cases <case> --serial-exe ./build/sis`
4. Run the regression test again to confirm it passes.
5. Commit the new reference set with a message explaining why it was generated.

## Adding a New Test Case

### Step 1: Prepare input files

Place input files (TOML, FASTA, PDB, force field, etc.) in an appropriate location:

- If the test reuses existing files from `check/`, reference them via relative paths in `required_files`.
- If the test needs new files, place them in `test/data/`.

The TOML input file must use capitalized section names (`[Job]`, `[Files]`, `[MD]`, etc.) as required by the current parser.

### Step 2: Add the test case definition

Edit `test/config.py` and add an entry to `TEST_CASES`:

```python
"my_new_test": {
    "description": "Short description of what this tests",
    "source_dir": "test/data",           # or "check/some_subdir"
    "input_toml": "input_my_test.toml",
    "required_files": [
        "structure.pdb",
        "sequence.fasta",
        "../../some_force_field.ff",     # relative to source_dir
    ],
    "modes": ["serial", "omp1", "ompN"],  # which modes are valid
    "mpi_ranks": 0,                       # 0 for non-MPI tests
    "output_prefix": "my_test",           # prefix of .out files
    "has_replica_cols": False,             # True for REMD tests
    "restart_test": True,                 # include in restart tests?
    "short_nstep": 50,                    # nstep for smoke tests
    "regression_nstep": 100,              # nstep for regression tests
    "nstep_save": 10,                     # output frequency
    "nstep_save_rst": 50,                 # restart save frequency
    "tolerances": {"atol": 1e-8, "rtol": 1e-6},
},
```

Key fields explained:

- **`source_dir`**: Directory where input files live, relative to repo root.
- **`required_files`**: Paths relative to `source_dir`. Checked before running; if any are missing, the test is skipped (not failed).
- **`modes`**: Only modes listed here will be attempted. Use `["serial", "omp1", "ompN"]` for non-MPI tests, `["mpi", "mpi_omp"]` for REMD.
- **`mpi_ranks`**: Number of MPI processes. Set to 0 for non-MPI tests.
- **`short_nstep`**: Used for the `run` level smoke test. Keep it small for speed.
- **`regression_nstep`**: Used for `regression` and `restart` levels. Should be long enough to produce meaningful output but short enough to run in seconds.
- **`tolerances`**: `atol` is absolute tolerance, `rtol` is relative tolerance. A value passes if either tolerance is satisfied (OR logic).

For REMD tests, also set:

- **`n_replicas`**: Total number of replicas.
- **`has_replica_cols`**: `True` (the `.out` files have extra replica label columns).
- **`restart_input_tomls`**: List of TOML files for restart stages, if the test case provides them.

### Step 3: Verify the test runs

```bash
./test/run_tests.py --levels run --cases my_new_test --modes serial --serial-exe ./build/sis
```

### Step 4: Generate reference data

```bash
./test/run_tests.py --generate-reference --cases my_new_test --serial-exe ./build/sis
```

### Step 5: Verify regression passes

```bash
./test/run_tests.py --levels regression --cases my_new_test --serial-exe ./build/sis
```

### Step 6: Commit

Commit the new input files, the config change, and the reference data together.

## Adding Sampling Properties

The `sampling` level checks thermodynamic averages from longer runs. To enable it for a test case, add `sampling_properties` to its config:

```python
"sampling_nstep": 100000,
"sampling_properties": [
    # (column_index, expected_mean, tolerance_in_sigma)
    (3, -42.5, 4.0),   # column 3 (Ekin), expected mean -42.5, 4-sigma tolerance
    (4, -120.3, 4.0),  # column 4 (Epot)
],
```

The framework computes the mean and standard error from the `.out` file (discarding the first 20% as equilibration) and checks that `|mean - expected| < tolerance_sigma * stderr`.

## Directory Structure

```
test/
  run_tests.py              Main CLI entry point (executable)
  config.py                 Test case definitions
  README.md                 This file
  data/                     Test-specific input files
    input_md_simple.toml
    T2HP.fasta
    T2HP.pdb
  reference/                Committed reference data sets
    0000-00-00_initial/     Labeled reference set
      summary.json          Metadata (date, commit, machine)
      md_simple/
        test_md.out
  framework/                Framework modules
    __init__.py
    utils.py                Subprocess, symlinks, MPI wrapper generation
    parser.py               .out file parsing and statistics
    comparator.py           Numerical comparison (absolute + relative)
    toml_patcher.py         TOML key-value patching and path localization
    reporter.py             Color-coded pass/fail/skip reporting
    builder.py              CMake build management for 3 variants
    runner.py               Test orchestration across all levels
  _work/                    Temporary work directories (git-ignored)
```

## How It Works Internally

For each (case, mode) combination:

1. Creates an isolated work directory under `test/_work/`.
2. Symlinks input files from their source locations (using relative paths).
3. Copies and patches the TOML file: adjusts `nstep`, `prefix`, and rewrites file paths to point to local symlinks.
4. For MPI modes, generates a wrapper script that auto-detects the MPI implementation (OpenMPI, MPICH, SLURM) for rank routing.
5. Runs the executable with appropriate environment (`OMP_NUM_THREADS`, etc.).
6. Checks results according to the requested test level.

## Troubleshooting

- **All tests SKIP**: The required force field files are probably missing. Check the SKIP message for which file is needed.
- **Test fails with "non-zero exit"**: Run with `--verbose` to see stdout/stderr, or inspect files in `test/_work/<case>_<mode>/` (use `--keep-work` to preserve them).
- **Regression fails after intentional change**: See "When to regenerate" above.
- **MPI tests fail**: Check that `mpirun` (or your MPI launcher) is available. Set `--mpi-cmd` if using a different launcher (e.g., `srun`).
