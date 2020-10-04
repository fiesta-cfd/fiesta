## Building FIESTA on Capulin (with OpenMP)

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
salloc -N1 --qos=debug --res=debug
```

Load Modules for Build:
```
module load cmake cray-hdf5-parallel
```

### Building for MPI and HDf5
Configure with CMAKE
```
cmake ../fiesta -DCMAKE_C_COMPILER="cc" -DCMAKE_CXX_COMPILER="CC"
```
The above configures FIESTA to build with mpi and hdf5 enabled.  BEWARE, the above command is case sensitive, the c compiler is named with two lowecase c's, the cxx compiler is named with two uppercase C's.  Thanks Cray.

If desired, specify an installation directory with '-DCMAKE_INSTALL_PREFIX=/path/to/install'.  If not specified, the executable will be built in the build directory.

Build the code:
```
make -j
```
or, if an installation directory was specified
```
make install -j
```

#### Run a MPI parallel interactive job
Copy a sample input file to a scratch directory.

```
cd /lustre/scratch3/$USER
mkdir fiesta-test && cd fiesta-test
cp ~/fiesta-dev/fiesta/test/ideal_expansion_2D/fiesta.lua .
```
Modify this file so that nprocsx=2 and nprocsy=2

Then run the code with:
```
mpiexec -n 4 ~/fiesta-dev/build/fiesta fiesta.lua
```
The simulation will produce solution files in the HDF5/XDMF format and can be viewed with Paraview, Tecplot or any other mainstream visualization package.  The format is fairly well standardized.
If there are MPI errors about buffer aliasing, then try
```
export MPICH_NO_BUFFER_ALIAS_CHECK=1
```

### Building LITE with OpenMP
Configure with CMAKE
```
cmake ../fiesta -DLITE=on -DOPENMP=on
```
The above gives a FIESTA build with no device parallel or MPI, only OpenMP host parallel.

If desired, specify an installation directory with '-DCMAKE_INSTALL_PREFIX=/path/to/install'.  If not specified, the executable will be built in the build directory.

Build the code:
```
make -j
```
or, if an installation directory was specified
```
make install -j
```

#### Run an interactive job with the LITE build of fiesta
Copy a sample input file to a scratch directory.

```
cd /lustre/scratch3/$USER
mkdir fiesta-test && cd fiesta-test
cp ~/fiesta-dev/fiesta/test/ideal_expansion_2D/fiesta.lua .
```

Then run the code with:
```
~/fiesta-dev/build/fiesta fiesta.lua
```
The simulation will produce restart files and solution files.  Both are in the VTK format and can be viewed with Paraview, Tecplot or any other mainstream visualization package.  The format is fairly well standardized.
