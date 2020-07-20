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
module load gcc/8.3.0 cmake/3.17.2
```

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

### Run an interactive job
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
