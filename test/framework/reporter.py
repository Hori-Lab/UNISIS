"""Test result collection and reporting."""

import sys


class TestResult:
    """A single test result."""

    PASS = "PASS"
    FAIL = "FAIL"
    SKIP = "SKIP"

    def __init__(self, case, mode, level, status, detail="", elapsed=0.0):
        self.case = case
        self.mode = mode
        self.level = level
        self.status = status
        self.detail = detail
        self.elapsed = elapsed


class TestReporter:
    """Collect and display test results."""

    def __init__(self, verbose=False):
        self.results = []
        self.verbose = verbose

    def add(self, result):
        """Add a TestResult."""
        self.results.append(result)
        # Print immediately for real-time feedback
        self._print_result(result)

    def _print_result(self, r):
        status_str = r.status
        if r.status == TestResult.PASS:
            status_str = f"\033[32m{r.status}\033[0m"
        elif r.status == TestResult.FAIL:
            status_str = f"\033[31m{r.status}\033[0m"
        elif r.status == TestResult.SKIP:
            status_str = f"\033[33m{r.status}\033[0m"

        mode_str = r.mode if r.mode else "---"
        detail = r.detail
        if r.elapsed > 0 and not detail:
            detail = f"({r.elapsed:.1f}s)"
        elif r.elapsed > 0:
            detail = f"({r.elapsed:.1f}s) {detail}"

        print(f"  {r.case:<18s} {mode_str:<10s} {r.level:<13s} {status_str}  {detail}")

    def print_summary(self):
        """Print final summary and return exit code."""
        n_pass = sum(1 for r in self.results if r.status == TestResult.PASS)
        n_fail = sum(1 for r in self.results if r.status == TestResult.FAIL)
        n_skip = sum(1 for r in self.results if r.status == TestResult.SKIP)

        print()
        print("=" * 70)
        parts = []
        if n_pass:
            parts.append(f"\033[32m{n_pass} passed\033[0m")
        if n_fail:
            parts.append(f"\033[31m{n_fail} failed\033[0m")
        if n_skip:
            parts.append(f"\033[33m{n_skip} skipped\033[0m")
        print(f"SUMMARY: {', '.join(parts)}")
        print("=" * 70)

        if n_fail > 0:
            return 1
        return 0

    def print_header(self):
        """Print the results table header."""
        print()
        print("=" * 70)
        print("UNISIS Test Results")
        print("=" * 70)
        print(f"  {'CASE':<18s} {'MODE':<10s} {'LEVEL':<13s} STATUS  DETAILS")
        print("-" * 70)
