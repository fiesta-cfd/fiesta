# loboSHOK

try 
```
export OMPI_MCA_fs_ufs_lock_algorithm=1 
```
for nfs filesystems and openmpi if file writing hangs
see https://github.com/open-mpi/ompi/issues/4446

to compile with cuda on xena need to
```
export MPICH_CXX=/users/beromer/Kokkos/kokkos/bin/nvcc_wrapper
```

and kokkos "generate_makefile.sh" should be run with "--with-cuda" and "--arch=Kepler35" and "--kokkos-path=" and "--prefix=". kokkos path needs to be in Makefile

get kokkos from github.com/kokkos/kokkos.git
get cgns-3.4 and lua-5.3 from spack install cgns with mpich and gcc-7.4.0
```
module load gcc-7.4.0-gcc-8.1.0-j26pfmd
spack compiler find
spack install cgns%gcc@7.4.0 ^mpich@3.3.1
spack install lua@5.3%gcc@7.4.0
```

use cuda 10.0 with gcc 7.4.0 built in modules on xena


example multi-gpu nvprof
```
mpirun -np 2 nvprof --output-profile profile.%q{OMPI_COMM_WORLD_RANK} ../loboshok.cuda input.lua --kokkos-ndevices=2
```
