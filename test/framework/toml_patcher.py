"""Lightweight TOML key-value patcher.

Does not parse TOML fully. Matches lines of the form 'key = value'
and replaces values for specified keys. Works for the flat keys used
in UNISIS input files (nstep, nstep_save, nstep_save_rst, prefix).
"""

import os
import re


def patch_toml(input_path, output_path, patches):
    """Read a TOML file, patch specified key-value pairs, write result.

    Parameters
    ----------
    input_path : str
        Path to the original TOML file.
    output_path : str
        Path to write the patched TOML file.
    patches : dict
        Keys to patch and their new values. Values can be int, float, or str.
        String values are automatically quoted.
    """
    with open(input_path) as f:
        lines = f.readlines()

    patched_lines = []
    for line in lines:
        patched_line = line
        for key, value in patches.items():
            # Match lines like:  key = value  or  key=value
            # Allowing leading whitespace and optional comments
            pattern = rf'^(\s*{re.escape(key)}\s*=\s*)(.*)$'
            match = re.match(pattern, line)
            if match:
                prefix = match.group(1)
                if isinstance(value, str):
                    new_value = f'"{value}"'
                elif isinstance(value, bool):
                    new_value = "true" if value else "false"
                elif isinstance(value, int):
                    new_value = str(value)
                elif isinstance(value, float):
                    new_value = str(value)
                else:
                    new_value = str(value)
                patched_line = f"{prefix}{new_value}\n"
                break
        patched_lines.append(patched_line)

    with open(output_path, "w") as f:
        f.writelines(patched_lines)


def localize_file_paths(input_path, output_path):
    """Rewrite file path values in [files.in] to use local basenames.

    Converts paths like '../rna_cg1.ff' or '../../htv23_5.0_ee-17.ff'
    to './rna_cg1.ff' so they resolve to symlinks in the work directory.

    Only affects keys that look like file paths (contain '/' or '..')
    within the [files.in] or [Files.In] section.
    """
    # Keys commonly found in [files.in]
    file_keys = {"ff", "fasta", "pdb_ini", "xyz_ini", "dcd", "bpl", "bpseq",
                 "ct", "anneal", "restraint", "timed_bias_rg"}

    with open(input_path) as f:
        lines = f.readlines()

    in_files_section = False
    result_lines = []

    for line in lines:
        stripped = line.strip().lower()

        # Track if we're inside [files.in] or [Files.In]
        if stripped.startswith("["):
            in_files_section = (
                stripped in ("[files.in]", "[files.in ]")
                or "files" in stripped and "in" in stripped
                and not "out" in stripped
            )

        if in_files_section:
            for key in file_keys:
                pattern = rf'^(\s*{key}\s*=\s*)"([^"]*)"(.*)$'
                match = re.match(pattern, line, re.IGNORECASE)
                if match:
                    prefix = match.group(1)
                    path_val = match.group(2)
                    rest = match.group(3)
                    basename = os.path.basename(path_val)
                    line = f'{prefix}"./{basename}"{rest}\n'
                    break

        result_lines.append(line)

    with open(output_path, "w") as f:
        f.writelines(result_lines)
