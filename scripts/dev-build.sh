#!/usr/bin/env bash
set -euo pipefail

build_serial=1
build_mpi=0
create_links=0
build_type="${CMAKE_BUILD_TYPE:-Release}"
cmake_bin="${CMAKE:-cmake}"

usage() {
    cat <<'USAGE'
Usage: scripts/dev-build.sh [--serial|--mpi|--both] [--links] [--build-type TYPE]

Build UNISIS for development in ./build without installing.

Options:
  --serial           Build only the serial executable (default)
  --mpi              Build only the MPI executable
  --both             Build serial and MPI executables
  --links            Create repository-root links to built executables
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
        --links)
            create_links=1
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
build_dir="${repo_root}/build"

"${cmake_bin}" -S "${repo_root}" -B "${build_dir}" \
    -DCMAKE_BUILD_TYPE="${build_type}" \
    -DUNISIS_BUILD_SERIAL="$([[ ${build_serial} -eq 1 ]] && echo ON || echo OFF)" \
    -DUNISIS_BUILD_MPI="$([[ ${build_mpi} -eq 1 ]] && echo ON || echo OFF)"

"${cmake_bin}" --build "${build_dir}" -j

if [[ ${create_links} -eq 1 ]]; then
    if [[ ${build_serial} -eq 1 ]]; then
        ln -sfn build/bin/unisis "${repo_root}/unisis"
    fi
    if [[ ${build_mpi} -eq 1 ]]; then
        ln -sfn build/bin/unisis_mpi "${repo_root}/unisis_mpi"
    fi
fi

echo "Built executables in ${build_dir}/bin"
