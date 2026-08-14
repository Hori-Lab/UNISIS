#!/bin/sh

## OpenMPI
#rank=$OMPI_COMM_WORLD_RANK

## mpich
rank=$PMI_RANK

suffix=`printf "%3.3d" $rank`

../../build/bin/unisis_mpi input_REMD.toml 1> ./out.$suffix 2> ./err.$suffix
