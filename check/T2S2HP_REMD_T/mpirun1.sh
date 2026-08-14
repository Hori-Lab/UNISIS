#!/bin/sh

## OpenMPI
#rank=$OMPI_COMM_WORLD_RANK

## mpich
rank=$PMI_RANK

suffix=`printf "%3.3d" $rank`

../../build/bin/unisis_mpi input_REMD1.toml 1> ./out1.$suffix 2> ./err1.$suffix
