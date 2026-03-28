"""Build manager for cmake-based compilation of UNISIS variants."""

import os
import multiprocessing
from .utils import run_command, ensure_dir


# Maps build variant to (directory_name, cmake_extra_flags)
BUILD_VARIANTS = {
    "serial": ("build_serial", ["-DCMAKE_DISABLE_FIND_PACKAGE_OpenMP=TRUE"]),
    "omp": ("build", []),
    "mpi": ("build_mpi", ["-DBUILD_MPI=ON"]),
}

# Maps parallelization mode to build variant
MODE_TO_BUILD = {
    "serial": "serial",
    "omp1": "omp",
    "ompN": "omp",
    "mpi": "mpi",
    "mpi_omp": "mpi",
}


class BuildManager:
    """Manage cmake builds for different parallelization variants."""

    def __init__(self, repo_root, build_type="Release", verbose=False):
        self.repo_root = repo_root
        self.build_type = build_type
        self.verbose = verbose
        self.executables = {}  # variant -> path

    def set_exe(self, variant, path):
        """Set an explicit executable path (from CLI --*-exe flags)."""
        if not os.path.isfile(path):
            raise FileNotFoundError(f"Executable not found: {path}")
        self.executables[variant] = os.path.abspath(path)

    def get_exe(self, mode):
        """Get path to executable for a parallelization mode.

        Returns None if the variant hasn't been built or set.
        """
        variant = MODE_TO_BUILD[mode]
        return self.executables.get(variant)

    def need_build(self, mode):
        """Check if a build is needed for the given mode."""
        variant = MODE_TO_BUILD[mode]
        if variant in self.executables:
            return False
        # Check if executable already exists from a previous build
        dir_name = BUILD_VARIANTS[variant][0]
        exe_path = os.path.join(self.repo_root, dir_name, "sis")
        if os.path.isfile(exe_path):
            self.executables[variant] = exe_path
            return False
        return True

    def build(self, variant):
        """Build a specific variant. Returns path to executable.

        Raises RuntimeError on build failure.
        """
        if variant in self.executables:
            return self.executables[variant]

        dir_name, extra_flags = BUILD_VARIANTS[variant]
        build_dir = os.path.join(self.repo_root, dir_name)
        ensure_dir(build_dir)

        # cmake
        cmake_cmd = [
            "cmake", self.repo_root,
            f"-DCMAKE_BUILD_TYPE={self.build_type}",
        ] + extra_flags

        # For MPI, try to use mpifort as the Fortran compiler
        if variant == "mpi":
            cmake_cmd.insert(1, "-DCMAKE_Fortran_COMPILER=mpifort")

        print(f"  Building {variant} in {dir_name}/ ...")
        rc, stdout, stderr, _ = run_command(cmake_cmd, cwd=build_dir, timeout=120)
        if rc != 0:
            msg = f"cmake failed for {variant}:\n{stderr}"
            if self.verbose:
                msg += f"\nstdout:\n{stdout}"
            raise RuntimeError(msg)

        # make
        njobs = min(multiprocessing.cpu_count(), 12)
        make_cmd = ["make", f"-j{njobs}"]
        rc, stdout, stderr, elapsed = run_command(make_cmd, cwd=build_dir, timeout=600)
        if rc != 0:
            msg = f"make failed for {variant}:\n{stderr}"
            if self.verbose:
                msg += f"\nstdout:\n{stdout}"
            raise RuntimeError(msg)

        exe_path = os.path.join(build_dir, "sis")
        if not os.path.isfile(exe_path):
            raise RuntimeError(f"Build succeeded but executable not found: {exe_path}")

        self.executables[variant] = exe_path
        print(f"  Built {variant}: {exe_path} ({elapsed:.1f}s)")
        return exe_path

    def build_for_mode(self, mode):
        """Build the variant needed for a parallelization mode."""
        variant = MODE_TO_BUILD[mode]
        return self.build(variant)

    def build_all_needed(self, modes):
        """Build all variants needed for the given set of modes."""
        variants_needed = set()
        for mode in modes:
            variant = MODE_TO_BUILD[mode]
            if variant not in self.executables:
                variants_needed.add(variant)

        # Check if already built on disk
        for variant in list(variants_needed):
            dir_name = BUILD_VARIANTS[variant][0]
            exe_path = os.path.join(self.repo_root, dir_name, "sis")
            if os.path.isfile(exe_path):
                self.executables[variant] = exe_path
                variants_needed.discard(variant)

        if variants_needed:
            print(f"\nBuilding required variants: {', '.join(sorted(variants_needed))}")
            for variant in sorted(variants_needed):
                self.build(variant)
