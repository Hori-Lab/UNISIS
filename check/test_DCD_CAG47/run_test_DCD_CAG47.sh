#!/bin/bash

set -e

mkdir -p build

cd build
cmake -DCMAKE_BUILD_TYPE=Debug -DCMAKE_Fortran_FLAGS="-O0 -g -fbacktrace -ffpe-trap=zero,overflow,denormal,underflow -Wall -ffree-line-length-none -fbounds-check -fmax-errors=5" ../../..
make -j12
cd ..

mkdir -p build_HTN

cd build_HTN
cmake -DCMAKE_BUILD_TYPE=Debug -DCMAKE_Fortran_FLAGS="-O0 -g -fbacktrace -ffpe-trap=zero,overflow,denormal,underflow -Wall -ffree-line-length-none -fbounds-check -fmax-errors=5 -D_HTN_CONSISTENT" ../../..
make -j12
cd ..

./build/bin/unisis ./input_DCD_CAG47.toml 1> sis.out 2> sis.err
./build_HTN/bin/unisis ./input_DCD_CAG47_HTN.toml 1> sis_HTN.out 2> sis_HTN.err

gnuplot compare_DCD_CAG47.gnu
