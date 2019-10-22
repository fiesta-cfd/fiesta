# FIESTA
**F**ast **I**nterfaces, **E**xplosions, **S**hocks and **T**ransport in the **A**tmosphere

## Install

### Dependency Installation on Xena

Compiling with Cuda support on Xena requires cgns-3.4, lua-5.3 and kokkos which are not already present on the system.
Spack can be used to install cgns and Lua on Xena with mpich and gcc-7.4.0.  Newer gcc versions cannot be used because they are not supported by Cuda.

```
module load gcc-7.4.0-gcc-8.1.0-j26pfmd
spack compiler find
spack install cgns%gcc@7.4.0 ^mpich@3.3.1
spack install lua@5.3%gcc@7.4.0
```
The system provided mpich or openmpi installationis can be used by setting their paths in .spack/packages.yaml as appropriate.  The above commands will build mpich-3.3.1 from source.

Kokkos can be obtained from github.com/kokkos/kokkos.git
```
git clone https://github.com/kokkos/kokkos.git
cd $BUILD_DIRECTORY
$KOKKOS_SOURCE_DIR/generate_makefile.sh --with-cuda --with-openmp --with-serial --arch=Kepler35 --kokkos_cuda_opt=enable_lambda --kokkos-path=$KOKKOS_SOURCE_DIR --prefix=$KOKKOS_INSTALL_DIR
```

use cuda 10.0 with gcc 7.4.0 built in modules on xena
```
export MPICH_CXX=/users/beromer/Kokkos/kokkos/bin/nvcc_wrapper
```

and kokkos "generate_makefile.sh" should be run with "--with-cuda" and "--arch=Kepler35" and "--kokkos-path=" and "--prefix=". kokkos path needs to be in Makefile


try 
```
export OMPI_MCA_fs_ufs_lock_algorithm=1 
```
for nfs filesystems and openmpi if file writing hangs. see https://github.com/open-mpi/ompi/issues/4446

example multi-gpu nvprof
```
mpirun -np 2 nvprof --output-profile profile.%q{OMPI_COMM_WORLD_RANK} ../loboshok.cuda input.lua --kokkos-ndevices=2
```
