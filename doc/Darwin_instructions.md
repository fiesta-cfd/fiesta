## Building FIESTA on Darwin

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

Get Interactive Node with NVidia GPU for CUDA:
```
salloc --constraint="gpu_vendor:nvidia"
```

Load Modules for Build:
```
module load cmake/3.17.3 cuda/10.2 gcc/8.3.0 openmpi/4.0.3-gcc_8.3.0
```

Configure with CMAKE
```
cmake ../fiesta -DCUDA=on
```
With no options, Fiesta will be built with CPU support.  Use '-DCUDA=on' for nvidia GPUs.  Use '-DOPENMP=on' for openMP support on CPUs.

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
cd /scratch/users/$USER
mkdir fiesta-test && cd fiesta-test
cp ~/fiesta-dev/fiesta/test/ideal_expansion_2D/fiesta.lua .
```

Then run the code with:
```
~/fiesta-dev/build/fiesta fiesta.lua
```
The simulation will produce restart files and solution files.  Both are in the CGNS format and can be viewed with Paraview, Tecplot or any other mainstream visualization package.  The format is fairly well standardized.
