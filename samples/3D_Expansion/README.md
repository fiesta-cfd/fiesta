# Sample Input Deck for a 3D Expansion Problem

## Problem Description
This input deck defines an expanding bubble of hot gas.  A 0.5 meter diameter
zone of hot gas is defined in the center of the domain with a temperature of
3000 Kelvin and a density of 1000 kg/m^3.  The surrounding air has a temperature
of 300 Kelvin and a density of 1 kg/m^3.

## Problem discretization
The computational domain is a cube 10 meters on a side.  This cube is
discretized into 1000 computational cells in each direction for a total of
1X10^9 cells.  For this problem, the data structures in FIESTA require
approximately 1.5 GB of memory per 1X10^6 cells, so the minimum memory
requirements are 1000 GB.

The MPI decomposition is controlled by the parameters `procsx`, `procsy`, and
`procsz` (Lined 44-46 in fiesta.lua).  The defaults for this problem are 4x4x4
for a total of 64 MPI ranks.  This can be arbitrarily adjusted so long as the
minimum memory requirements are met.

The computational cells are then distributed into the MPI ranks using a
Cartesian topography.  If the number of cells in one direction is not evenly
divisible by the number of ranks in that direction, then the number of cells in
each rank will not be equal.

## Problem Duration
The duration of execution is controlled by the number of time steps set
with the parameter `nt` (line 29).  `nt` should be less than 40,000.  Larger
values will result in unrealistic solutions due to boundary
interactions.

On 16nodes/64gpus (NVIDIA Tesla P100), 40,000 timesteps takes about 2.5 hours
of wall time hours wall time.

There are two halo exchanges per time step.

## Input/Output

Solution files are written in hdf5 format with a corresponding xmf file used
for visualization software.  The hdf5 writer uses parallel, contiguous
(non-chunked) hdf5 write operations.

Solution write frequency is controlled by `write_freq`.  Each solution file for
this problem are approximately 53 GB.

The repository includes a sample logfile from running this problem,
`slurm-fiesta.log`.  Comparisons can be made to this file to check correctness.

## Compiling Fiesta
A compiler and mpi implementation is required.  Cuda is required for running the
code on GPUs.  A parallel hdf5 implementation is recommended, but FIESTAs build
system will compile its own hdf5 if no system package is detected.  CMake is
required for building FIESTA.

Known working combination:

LANL Kodiak Machine
cmake/3.17.2
gcc/8.3.0
openmpi/3.1.6
cuda/10.2
hdf5-parallel/1.8.16

Create a build directory somewhere and configure the code with:

```
cmake /path/to/fiesta/ -DCUDA=on
```

Omit `-DCUDA=on` for cpu only code.

Build the code with:

```
make -j
```
This will produce an executable `fiesta`.

Information about the build configuration can be obtained with `fiesta
--version`.

## Running The Problem
mpirun -n NPROCS path/to/fiesta fiesta.lua -n NGPUS

where NPROCS = procsx*procsy*procsz
and NGPUS = number of gpus per node

The number of tasks per node should be set to the number of GPUs per node.

The input file `fiesta.lua` is a lua script.  The `os.getenv()` lua command is
useful for parametrizing jobs. For example with slurm job arrays.
