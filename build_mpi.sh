set -eu

mkdir -p build_mpi
cd build_mpi

hostname=`hostname`
un=`uname`
if [[ $un == 'Linux' ]]; then
    hostnameall=`hostname -A`
else
    hostnameall=''
fi

if [[ ${hostnameall:(-14):(-1)} == 'archer2.ac.uk' ]]; then
     FC=mpifort cmake .. -DCMAKE_BUILD_TYPE=Release -DALLOW_ARGUMENT_MISMATCH=ON -DBUILD_MPI=ON
#if [[ ${hostnameall:(-14):(-1)} == 'archer2.ac.uk' ]]; then
#    FC=ftn cmake .. -DCMAKE_BUILD_TYPE=Release -DMT_USE_ALLOWINVALIDBOZ=OFF -DBP_HALT_IEEE_EXCEPTIONS=OFF -DBUILD_MPI=ON
elif [[ ${hostname:8:9} == 'sulis.hpc' ]]; then
    FC=mpiifort cmake .. -DCMAKE_BUILD_TYPE=Release -DBUILD_MPI=ON
elif [[ ${hostname:(-20)} == 'ada.nottingham.ac.uk' ]]; then
    FC=mpiifort cmake .. -DCMAKE_BUILD_TYPE=Release -DBUILD_MPI=ON
elif [[ ${hostname:13:24} == 'augusta.nottingham.ac.uk' ]]; then
    #FC=ifort cmake .. -DCMAKE_BUILD_TYPE=Release
    echo 'Edit build_mpi.sh!'
    exit 2
elif [[ ${hostname:12:5} == 'cosma' ]]; then
    FC=mpiifort cmake .. -DCMAKE_BUILD_TYPE=Release -DBUILD_MPI=ON
else
    FC=mpifort cmake .. -DCMAKE_BUILD_TYPE=Release -DBUILD_MPI=ON
fi

make -j8

cd ../
ln -sf build_mpi/sis sis.mpi
