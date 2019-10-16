# loboSHOK

try export OMPI_MCA_fs_ufs_lock_algorithm=1 
for nfs filesystems and openmpi if file writing hangs
see https://github.com/open-mpi/ompi/issues/4446

to compile with cuda need to
    export OMPI_CXX=/users/beromer/Kokkos/kokkos/bin/nvcc_wrapper
    and kokkos compiled with "--with-cuda" and "--arch=Kepler35" and "--kokkos-path=" and "--prefix="
    kokkos path needs to be in Makefile

get kokkos from github.com/kokkos/kokkos.git
get cgns-3.4 and lua-5.3 from spack
use cuda 10.0 with gcc 7.4.0
use mpich-3.3.1 on xena


example multi-gpu nvprof
mpirun -np 2 nvprof --output-profile profile.%q{OMPI_COMM_WORLD_RANK} ../loboshok.cuda input.lua --kokkos-ndevices=2
