###############################################################################
Build Guide
###############################################################################

*******************************************************************************
Supported Platforms
*******************************************************************************


*******************************************************************************
Building
*******************************************************************************

===============================================================================
1. Superbuild
===============================================================================
The easiest way to install Fiesta is with the superbuild option.  This option
will install Fiesta along with several of it's dependencies including Kokkos,
HDF5 and Lua.  This method reduces the dependence on pre-installed packages.

Requirements
-------------------------------------------------------------------------------
While the superbuild option reduces the dependency requirements, several
pre-installed packages are still rquired.

* | Compiler (Required)
  | A C++ compiler which supports C++17 is required such as GNU(g++), LLVM(clang++), and Intel(icpc).
  | GNU g++ 10 or greater is recommended.  
  | Chack that g++ version 10 or later is installed with :code:`g++ --version`.

* | CMake (Required)
  | CMake version 3.12 or greater is required.
  | Check the cmake version with :code:`cmake --version`

* | Make Program (required)
  | GNU Make or Ninja-build are required.
  | Check that make version 3.82 or greater is installed with :code:`make --version` or ninja version 1.10 or later with :code:`ninja --version`

* | MPI (Required if Fiesta will be run on multiple nodes and/or GPUs)
  | MPI may be provided by OpenMPI, MPICH or Intel MPI. OpenMPI is recommended.
  | Check that OpenMPI version 3.1 or greater is installed with :code:`ompi_info --version`.
  | OpenMPI development tools must also be available, check that the following command produces output :code:`mpicc --version`.

* | CUDA (Required if Fiesta will be run on NVIDIA GPUs)
  | Check that CUDA version 11.0 or greater is available with :code:`nvcc --version`

* | HIP (Required for AMD GPUs)
  | Check that HIP version 4.0 or greater is installed with :code:`hipcc --version`

Preparation
-------------------------------------------------------------------------------
Before, configuring and installing fiesta, a directory structure must be created
and the source code must be downloaded.  The recommended method for downloading
the source code is with git clone, but pre-packaged versions of the code can be
downloaded from `<github.com/fiesta-cfd/fiesta/releases>`_

1. | Create a working directory and clone the Fiesta repository.
   | :code:`mkdir work_dir && cd work_dir` (`work_dir` is an example. Any directory name and location may be chosen)
   | :code:`git clone https://github.com/fiesta-cfd/fiesta`

2. | Create a build directory.  Out of tree builds are suggested (where the build
   | directory is outside of the source code directory) but in-tree builds are
   | supported.
   | :code:`cd ../`
   | :code:`mkdir build && cd build`

3. The resulting directory structure should look something like this:
   ::

      work_dir
      |
      |->fiesta
      |->build

  where `fiesta` contains the fiesta source code and `build` is empty.

Configuration
-------------------------------------------------------------------------------
Fiesta may now be configured.  Configuration is done with the CMAKE tool.
Configuration options include install location, single or multi-node builds and
GPU type.

To configure Fiesta for different situations, run one of the following configure
commands from the build directory.

* | For multi-gpu support for NVIDIA GPUs:
  | :code:`cmake ../fiesta -DCUDA=on`

* | For multi-gpu support for AMD GPUs:
  | :code:`cmake ../fiesta -DHIP=on`

* | For single-cpu, serial support (also known as a lite build):
  | :code:`cmake ../fiesta -DLITE=on`

* | For single-gpu support:
  | :code:`cmake ../fiesta -DLITE=on -DCUDA=on`

* | For multi-cpu support:
  | :code:`cmake ../fiesta`

* | To specify an install directory, use `-DCMAKE_INSTALL_PREFIX`. e.g.:
  | :code:`cmake ../fiesta -DCUDA=on -DCMAKE_INSTALL_PREFIX=~/.local`
  | This configures Fiesta for multi-gpu support for NVIDIA gpus and specifies that the code should be installed to the users local bin directory.

* | To use the ninja build generator, use `-G`:
  | :code:`cmake -G Ninja ../fiesta ...`

After configuration is complete, there are several ways to make changes or corrections:

1. With the `ccmake` tool, if it is installed.  Just run :code:`ccmake .` from the build directory.

2. | by running the cmake command again, but prepended with `--clean-first`.  e.g.:
   | :code:`cmake --clean-first ../fiesta -D...`

3. by deleting the contents of the build directory and running another `cmake` command.

Installation
-------------------------------------------------------------------------------
Fiesta may now be compiled and installed.

To compile with GNU Make (the default, if Ninja was not explicitly specified) run:
:code:`make -j`
The `-j` option indicates that as many CPU cores as are available will be used
to speed up compilation.  :code:`make -j 4` would limit the compilation to four
cpus for example.

If Ninja was explicitly spcified in the configure step, then, from the build directory, run:
:code:`ninja` to build in parallel or :code:`ninja -j 4` to limit the
compilation to four CPU cores (for example).

If a specific installation path was specifies with `-DCMAKE_INSTALL_PREFIX`,
then fiesta can now be installed with :code:`make install` of :code:`ninja
install`.


===============================================================================
2. Site Build
===============================================================================

Requirements
-------------------------------------------------------------------------------

Preparation
-------------------------------------------------------------------------------

Configuration
-------------------------------------------------------------------------------

Installation
-------------------------------------------------------------------------------

===============================================================================
3. Containers
===============================================================================

Requirements
-------------------------------------------------------------------------------

Preparation
-------------------------------------------------------------------------------

Configuration
-------------------------------------------------------------------------------

Installation
-------------------------------------------------------------------------------

*******************************************************************************
Testing
*******************************************************************************

*******************************************************************************
Tools
*******************************************************************************
