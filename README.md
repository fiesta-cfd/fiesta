# FIESTA
**F**ast **I**nterfaces, **E**xplosions, **S**hocks and **T**ransport in the **A**tmosphere

## Installation on Xena at UNM Carc

Compiling with Cuda support on Xena requires cgns-3.4, lua-5.3 and kokkos which are not already present on the system.
Spack can be used to install cgns and Lua on Xena with mpich and gcc-7.4.0.  Newer gcc versions cannot be used because they are not supported by Cuda 10.0.

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

Now the run command scripts can be used to setup the environment.  These may need to be edited to reflect your module names.  These files load modules and export the mpi compiler environment variable

```
source mpich.rc
```

The makefile will need to be edited to reflect the kokkos paths and the installation path.  The code can then be built with 

```
make install -j
```

## Running on Xena at UNM Carc
```
qsub -I -q dualGPU -l walltime=48:00:00
fiesta.cuda input.lua --kokkos-ndevices=2
```

There is an example PBS batch script included.

## Other Commands
For NFS filesystems (like at UNM CARC) file writing may hang when using OpenMPI.  Try the following if this happens. See: https://github.com/open-mpi/ompi/issues/4446
```
export OMPI_MCA_fs_ufs_lock_algorithm=1 
```

example multi-gpu nvprof
```
mpirun -np 2 nvprof --output-profile profile.%q{OMPI_COMM_WORLD_RANK} ../fiesta.cuda input.lua --kokkos-ndevices=2
```
The profile.n files can now be viewed in nvvp.
