![](logo.jpg?raw=true)
# FIESTA
**F**ast **I**nterfac**ES** and **T**ransport in the **A**tmosphere

## Building on Kodiak

There is a new bash script and makefile system for building Fiesta.  Tested and working on Kodiak and local workstations.  Does not work on any UNM CARC machines.

Pick a working directory:
```
mkdir ~/fiesta-dev && cd ~/fiesta-dev
```

Get fiesta:
```
git clone git@github.com:beromer/fiesta.git
```

Create Build Directory:
```
mkdir build && cd build
```

Get Interactive Node:
```
salloc -N1 -t 4:00:00
```

Load Required Modules:
```
module load gcc/7.4.0 openmpi/2.1.2 cudatoolkit/10.0 cmake/3.14.6
module load hdf5-parallel/1.8.16
```

Configure and Generate Makefiles
```
../fiesta/configure.sh --device=Cuda --arch=Pascal60
```
Options for device are currently 'Serial' and 'Cuda' for cpu and gpu versions of the code.  For Kodiak use '--arch=BDW' for 'Serial' device and '--arch=Pascal60' for 'Cuda' device.  'BDW' stands for Intel Broadwell architecture.

Build the code:
```
make
```

An executable should now exist in the 'fiesta-build' subdirectory.  However, an issue with the static cgns library prevents static linking, so we have to include the dynamic library in the path using:
```
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:~/fiesta-dev/build/cgns-build/lib
```

This must also be included in any batch script. See fiesta repository root directory for example pbs and slurm batch scripts.

### Run an interactive job

Assuming that the code is built, the library path has been updated and the modules are still loaded, then running on a single interactive node is straight forward.
```
cd /your/users/scratch/directory
mkdir fiesta-test && cd fiesta-test
cp /fiesta/repository/path/test/ideal_expansion_3D/input.lua .
```

Edit the 'input.lua' file to change the number of mpi processes for 4 GPUs.  e.g.:
```
procsx = 2
procsy = 2
procsz = 1
```

Then run the code with:
```
mpirun -n 4 ~/fiesta-dev/build/fiesta-build/fiesta input.lua --kokkos-ndevices=4
```
The simulation will produce restart files and solution files.  Both are in the CGNS format and can be viewed with Paraview, Tecplot or any other mainstream visualization package.  The format is fairly well standardized.

### Run a batch job

Submit a batch script:
```
cd /your/users/scratch/directory
mkdir fiesta-test && cd fiesta-test
cp /fiesta/repository/path/fiesta.slurm .
cp /fiesta/repository/path/test/ideal_expansion_3D/input.lua .
```

Edit the 'fiesta.slurm' batch file to reflect your directories. Then submit with:
```
sbatch fiesta.slurm
```

Once the batch file executes, the simulation will produce restart files and solution files.  Both are in the CGNS format and can be viewed with Paraview, Tecplot or any other mainstream visualization package.  The format is fairly well standardized.


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
$KOKKOS_SOURCE_DIR/generate_makefile.sh --with-cuda --arch=Kepler35 --kokkos_cuda_opt=enable_lambda --kokkos-path=$KOKKOS_SOURCE_DIR --prefix=$KOKKOS_INSTALL_DIR
```

Now the run command scripts can be used to setup the environment.  These may need to be edited to reflect your module names.  These files load modules and export the mpi compiler environment variable

```
source mpich.rc
```

The makefile in the fiesta/src directory will need to be edited to reflect the kokkos paths and the installation path.  The code can then be built with 

```
make install -j
```

### Running on Xena at UNM Carc
```
qsub -I -q dualGPU -l walltime=48:00:00
fiesta.cuda input.lua --kokkos-ndevices=2
```

There is an example pbs batch script that will run on Xena.  The modules can be modified to  run on Wheeler if the '--with-serial' was used during the kokkos makefile generation step instead of --with-cuda.

There is an example PBS batch script included.

### Visualizing
Output files are written in the CFD Generalized Notation System standard (CGNS) and can be viewed with paraview, tecplot, visit, etc.  Solution files are single precision and include the grid in every solution file.  Restart files are double precision.

### Other Comments
For NFS filesystems (like at UNM CARC) file writing may hang when using OpenMPI.  Try the following if this happens. See: https://github.com/open-mpi/ompi/issues/4446
```
export OMPI_MCA_fs_ufs_lock_algorithm=1 
```

example multi-gpu nvprof
```
mpirun -np 2 nvprof --output-profile profile.%q{OMPI_COMM_WORLD_RANK} ../fiesta.cuda input.lua --kokkos-ndevices=2
```
The profile.n files can now be viewed in nvvp.
