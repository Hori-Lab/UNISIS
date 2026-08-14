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

../../build/bin/unisis_mpi input_REMD.toml 1> ./out.$suffix 2> ./err.$suffix
