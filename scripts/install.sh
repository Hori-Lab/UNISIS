#!/usr/bin/env bash
set -euo pipefail

build_serial=1
build_mpi=0
prefix="${HOME}/.local"
build_type="${CMAKE_BUILD_TYPE:-Release}"
cmake_bin="${CMAKE:-cmake}"

usage() {
    cat <<'USAGE'
Usage: scripts/install.sh [--serial|--mpi|--both] [--prefix PREFIX] [--build-type TYPE]

Cleanly build and install UNISIS through CMake.

Options:
  --serial           Install only the serial executable (default)
  --mpi              Install only the MPI executable
  --both             Install serial and MPI executables
  --prefix PREFIX    Installation prefix, default: ~/.local
  --build-type TYPE  CMake build type, default: Release
  -h, --help         Show this help
USAGE
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --serial)
            build_serial=1
            build_mpi=0
            ;;
        --mpi)
            build_serial=0
            build_mpi=1
            ;;
        --both)
            build_serial=1
            build_mpi=1
            ;;
        --prefix)
            shift
            if [[ $# -eq 0 ]]; then
                echo "Error: --prefix requires a value" >&2
                exit 2
            fi
            prefix="$1"
            ;;
        --build-type)
            shift
            if [[ $# -eq 0 ]]; then
                echo "Error: --build-type requires a value" >&2
                exit 2
            fi
            build_type="$1"
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "Error: unknown argument: $1" >&2
            usage >&2
            exit 2
            ;;
    esac
    shift
done

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
build_dir="${repo_root}/build-install"

rm -rf "${build_dir}"

"${cmake_bin}" -S "${repo_root}" -B "${build_dir}" \
    -DCMAKE_BUILD_TYPE="${build_type}" \
    -DUNISIS_BUILD_SERIAL="$([[ ${build_serial} -eq 1 ]] && echo ON || echo OFF)" \
    -DUNISIS_BUILD_MPI="$([[ ${build_mpi} -eq 1 ]] && echo ON || echo OFF)"

"${cmake_bin}" --build "${build_dir}" -j
"${cmake_bin}" --install "${build_dir}" --prefix "${prefix}"

echo "Installed UNISIS into ${prefix}"
