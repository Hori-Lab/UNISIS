#!/bin/sh

mpiver=`mpirun --version`

if [[ ${mpiver:8:8} == 'Open MPI' ]]; then
    ## OpenMPI
    rank=$OMPI_COMM_WORLD_RANK
else
    ## mpich
    rank=$PMI_RANK
fi

suffix=`printf "%3.3d" $rank`
echo $suffix

../../build/bin/unisis_mpi input_REMD2.toml test1/test_REMD1.rst 1> ./out2.$suffix 2> ./err2.$suffix
