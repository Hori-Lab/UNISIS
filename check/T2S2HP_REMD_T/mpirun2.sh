#!/bin/sh

## OpenMPI
#rank=$OMPI_COMM_WORLD_RANK

## mpich
rank=$PMI_RANK

suffix=`printf "%3.3d" $rank`

../../build/bin/unisis_mpi input_REMD2.toml test1/test_REMD1.rst 1> ./out2.$suffix 2> ./err2.$suffix
