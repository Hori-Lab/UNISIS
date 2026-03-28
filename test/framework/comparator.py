"""Numerical comparison of .out files."""

import math
from .parser import parse_out_file, align_by_step


class ComparisonResult:
    """Result of comparing two .out files."""

    def __init__(self):
        self.passed = True
        self.max_abs_diff = 0.0
        self.max_rel_diff = 0.0
        self.failures = []  # list of (row_idx, col_idx, actual, expected, abs_diff, rel_diff)
        self.n_compared = 0

    def add_failure(self, row, col, actual, expected, abs_diff, rel_diff):
        self.passed = False
        self.failures.append((row, col, actual, expected, abs_diff, rel_diff))

    @property
    def summary(self):
        if self.passed:
            return f"max_rdiff={self.max_rel_diff:.2e}"
        n = len(self.failures)
        first = self.failures[0]
        return (f"{n} failures; first at row {first[0]} col {first[1]}: "
                f"got {first[2]:.6e}, expected {first[3]:.6e}, "
                f"rdiff={first[5]:.2e}")


def compare_out_files(actual_path, reference_path, columns=None,
                      atol=1e-8, rtol=1e-6, align_steps=False):
    """Compare two .out files column by column.

    Parameters
    ----------
    actual_path : str
        Path to the actual output file.
    reference_path : str
        Path to the reference output file.
    columns : list[int], optional
        0-based column indices to compare. If None, compare all columns.
    atol : float
        Absolute tolerance.
    rtol : float
        Relative tolerance.
    align_steps : bool
        If True, align rows by step number (column 0) instead of by position.

    Returns
    -------
    ComparisonResult
    """
    result = ComparisonResult()

    rows_a = parse_out_file(actual_path)
    rows_b = parse_out_file(reference_path)

    if not rows_a:
        result.add_failure(0, 0, 0, 0, 0, 0)
        result.failures[-1] = (0, 0, "EMPTY", "N/A", 0, 0)
        return result

    if not rows_b:
        result.add_failure(0, 0, 0, 0, 0, 0)
        result.failures[-1] = (0, 0, "N/A", "EMPTY", 0, 0)
        return result

    if align_steps:
        rows_a, rows_b = align_by_step(rows_a, rows_b)

    n_rows = min(len(rows_a), len(rows_b))

    for i in range(n_rows):
        row_a = rows_a[i]
        row_b = rows_b[i]
        n_cols = min(len(row_a), len(row_b))

        cols_to_check = columns if columns is not None else list(range(n_cols))

        for c in cols_to_check:
            if c >= n_cols:
                continue
            a = row_a[c]
            b = row_b[c]

            abs_diff = abs(a - b)
            denom = max(abs(a), abs(b), 1e-30)
            rel_diff = abs_diff / denom

            result.n_compared += 1
            result.max_abs_diff = max(result.max_abs_diff, abs_diff)
            result.max_rel_diff = max(result.max_rel_diff, rel_diff)

            # Pass if EITHER tolerance is met
            if abs_diff <= atol or rel_diff <= rtol:
                continue

            # Both tolerances exceeded — failure
            result.add_failure(i, c, a, b, abs_diff, rel_diff)

    return result


def compare_sampling(actual_path, column, expected_mean, tolerance_sigma=4.0,
                     skip_fraction=0.2):
    """Compare a thermodynamic average against an expected value.

    Parameters
    ----------
    actual_path : str
        Path to the .out file.
    column : int
        0-based column index.
    expected_mean : float
        Expected mean value.
    tolerance_sigma : float
        Number of standard errors for the pass threshold.
    skip_fraction : float
        Fraction of initial rows to skip for equilibration.

    Returns
    -------
    tuple of (bool, str)
        (passed, detail_message)
    """
    from .parser import compute_statistics

    rows = parse_out_file(actual_path)
    if not rows:
        return False, "empty output file"

    mean, stderr, n = compute_statistics(rows, column, skip_fraction)

    if n < 10:
        return False, f"too few samples ({n})"

    if stderr < 1e-30:
        # All values identical
        diff = abs(mean - expected_mean)
        if diff < 1e-10:
            return True, f"mean={mean:.6e}, exact match"
        return False, f"mean={mean:.6e}, expected={expected_mean:.6e}, zero variance"

    z_score = abs(mean - expected_mean) / stderr
    passed = z_score < tolerance_sigma

    detail = (f"mean={mean:.6e}, expected={expected_mean:.6e}, "
              f"stderr={stderr:.2e}, z={z_score:.1f} "
              f"({'<' if passed else '>'}{tolerance_sigma:.0f}sigma)")

    return passed, detail
