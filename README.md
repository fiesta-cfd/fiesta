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
With no options, Fiesta will be built with serial CPU support.  Use `-DCUDA=on` for nvidia GPUs.  Use `-DOPENMP=on` for openMP support.

If desired, specify an installation directory with `-DCMAKE_INSTALL_PREFIX=/path/to/install`.

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
cp ~/fiesta-dev/fiesta/test/idexp2d/fiesta.lua .
```

Edit the 'fiesta.lua' file to specify a domain decomposition appropriate for 4 GPUs.  e.g.:
```
procsx = 2
procsy = 2
```

Then run the code with:
```
mpirun -n 4 ~/fiesta-dev/build/fiesta fiesta.lua --color --kokkos-num-devices=4
```

### Visualizing
Output files are written in the HDF5 format and can be viewed with ParaView, tecplot, visit, etc with the associated XDMF file.  That is open the sol-xxxx.xmf files with your visualization application of choice.  Solution files are single precision and include the grid in every timestep.  Restart files are double precision.

### Options
* `--version` Print version and build info and exit
* `--color[=auto|on|off]` Enable color output.  Default is 'off'.  'auto' is used if a value is omitted.  'on' always colors output. 'off'' never colors output. 'auto' colors output only if a tty is detected.  WARNING: Some MPI configurations report a tty, even when output is being redirected.  If color is turned on in this case, then 'less -r' can be used to view colorized logs.
* `--time-format[=n]` Format to use for time reporting.  Default is 2. Omitted is 0. If n is zero, then use human readable formatting. E.g, 4h32m37.11s.  If n is positive, uses n decimal places with an exponential format to report time in seconds. E.g. "1.96e+4"
* `--kokkos-num-devices=n[,m]` Number of gpu devices available per node. Use this when running MPI jobs to assign GPU devices to processes on a node.  GPUs are assigned in a "round-robin" manner via MPI rank. An optional second argument allows for a device to be ignored which is useful for workstations with both display and compute devices.  See kokkos wiki for more information.

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
