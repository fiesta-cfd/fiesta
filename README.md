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
module load gcc/7.4.0 openmpi/2.1.2 cudatoolkit/10.1 cmake/3.14.6 hdf5-parallel/1.8.16
```

Configure with CMAKE
```
cmake ../fiesta -DCUDA=on
```
With no options, Fiesta will be built with serial CPU support.  Use '-DCUDA=on' for nvidia GPUs.  Use '-DOPENMP=on' for openMP support on CPUs.

If desired, specify an installation directory with '-DCMAKE_INSTALL_PREFIX=/path/to/install'.  If not specified, the executable will be built in the build directory.

Build the code:
```
make -j
```
or, if an installation directory was specified
```
make install -j
```

### Run an interactive job
Copy a sample input file to a scratch directory.

```
cd /lustre/scratch3/turquoise/<moniker>
mkdir fiesta-test && cd fiesta-test
cp ~/fiesta-dev/fiesta/test/ideal_expansion_2D/fiesta.lua .
```

Edit the 'fiesta.lua' file to change the number of mpi processes for 4 GPUs.  e.g.:
```
procsx = 2
procsy = 2
```

Then run the code with:
```
mpirun -n 4 ~/fiesta-dev/build/fiesta fiesta.lua --kokkos-num-devices=4
```

### Visualizing
Output files are written in the HDF5 format and can be viewed with paraview, tecplot, visit, etc with the associated xdmf file.  That is open the sol-xxxx.xmf files with your visualization application of choice.  Solution files are single precision and include the grid in every timestep.  Restart files are double precision.

### Options
* '--version' Print version and build info and exit
* '--color[=auto|on|off]' Enable color output.  Default is off.  Auto if value is omitted.  On always colors output, off never colorizes output, auto colorizes outpur only if output is detected as a tty.  Some mpi configurations report a tty, even when output is being redirected.  'less -r' can be used to view colorized logs.
* '--time-format[=n]' Format to use for time reporting.  Default is 2. Ommitted is 0. If n is zero, then use human readable formating. E.g, 4h32m37.11s.  If n is positive, us n decimal places with an exponential format. E.g. "%.ne".
* '--kokkos-num-devices=n' Number of gpu devices available per node.

### Other Comments
For NFS filesystems (like at UNM CARC) file writing may hang when using some versions of OpenMPI.  Try the following if this happens. See: https://github.com/open-mpi/ompi/issues/4446
```
export OMPI_MCA_fs_ufs_lock_algorithm=1 
```

example multi-gpu nvprof
```
mpirun -np 2 nvprof --output-profile profile.%q{OMPI_COMM_WORLD_RANK} ../fiesta.cuda input.lua --kokkos-ndevices=2
```
The profile.n files can now be viewed in nvvp.
