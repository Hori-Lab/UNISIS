"""Parse UNISIS .out energy trajectory files."""


def parse_out_file(path):
    """Parse a .out file into a list of rows.

    Each row is a list of floats. Header lines (starting with #) are skipped.
    Empty lines are skipped.

    Returns
    -------
    list[list[float]]
        Parsed rows. Each row has the same number of columns.
    """
    rows = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            tokens = line.split()
            try:
                row = [float(t) for t in tokens]
                rows.append(row)
            except ValueError:
                continue
    return rows


def extract_columns(rows, columns):
    """Extract specified columns from parsed rows.

    Parameters
    ----------
    rows : list[list[float]]
        Parsed rows from parse_out_file.
    columns : list[int]
        0-based column indices to extract.

    Returns
    -------
    list[list[float]]
        Rows with only the specified columns.
    """
    return [[row[c] for c in columns if c < len(row)] for row in rows]


def get_step_column(rows):
    """Get the step numbers (first column) from parsed rows."""
    return [int(row[0]) for row in rows if rows]


def align_by_step(rows_a, rows_b):
    """Align two sets of rows by their step column (column 0).

    Returns only rows present in both, matched by step number.

    Returns
    -------
    tuple of (list[list[float]], list[list[float]])
        Aligned rows from a and b.
    """
    step_map_b = {}
    for row in rows_b:
        step = int(row[0])
        step_map_b[step] = row

    aligned_a = []
    aligned_b = []
    for row in rows_a:
        step = int(row[0])
        if step in step_map_b:
            aligned_a.append(row)
            aligned_b.append(step_map_b[step])

    return aligned_a, aligned_b


def compute_statistics(rows, column, skip_fraction=0.2):
    """Compute mean and standard error for a column, skipping initial frames.

    Parameters
    ----------
    rows : list[list[float]]
        Parsed rows.
    column : int
        0-based column index.
    skip_fraction : float
        Fraction of initial rows to skip (equilibration).

    Returns
    -------
    tuple of (float, float, int)
        (mean, standard_error, n_samples)
    """
    n_skip = int(len(rows) * skip_fraction)
    values = [row[column] for row in rows[n_skip:] if column < len(row)]

    if len(values) < 2:
        return 0.0, float("inf"), 0

    n = len(values)
    mean = sum(values) / n
    variance = sum((v - mean) ** 2 for v in values) / (n - 1)
    stderr = (variance / n) ** 0.5

    return mean, stderr, n
