# loboSHOK

Fast, Portable, 3D Euler Solver

try
    export OMPI_MCA_fs_ufs_lock_algorithm=1 
for nfs filesystems like in xena if file writing is slow
see https://github.com/open-mpi/ompi/issues/4446

to compile with cuda need to
    export OMPI_CXX=/users/beromer/Kokkos/kokkos/bin/nvcc_wrapper


