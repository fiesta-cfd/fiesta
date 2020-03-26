![](logo.jpg?raw=true)
# FIESTA
**F**ast **I**nterfac**ES** and **T**ransport in the **A**tmosphere

## Building on Kodiak

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
```

Configure with CMAKE
```
cmake ../fiesta -DCUDA=on
```
With no options, Fiesta will be built with CPU support.  Use '-DCUDA=on' for nvidia GPUs.  Use '-DOPENMP=on' for openMP support on CPUs.

If desired, specify an installation directory with '-DCMAKE_INSTALL_PREFIX=/path/to/install'.

Build the code:
```
make -j
```
or
```
make install -j
```

The resulting executable is statically linked, so there is no need to keep the third party libraries in the path, e.g.  "export LD_LIBRARY_PATH" is not necessary. 

### Run an interactive job
Copy a sample input file to a scratch directory.

```
cd /your/users/scratch/directory
mkdir fiesta-test && cd fiesta-test
cp ~/fiesta-dev/fiesta/test/ideal_expansion_3D/fiesta.lua .
```

Edit the 'fiesta.lua' file to change the number of mpi processes for 4 GPUs.  e.g.:
```
procsx = 2
procsy = 2
procsz = 1
```

Then run the code with:
```
mpirun -n 4 ~/fiesta-dev/build/fiesta fiesta.lua --kokkos-ndevices=4
```
The simulation will produce restart files and solution files.  Both are in the CGNS format and can be viewed with Paraview, Tecplot or any other mainstream visualization package.  The format is fairly well standardized.

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
