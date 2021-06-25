###############################################################################
Build Guide
###############################################################################

*******************************************************************************
Supported Platforms
*******************************************************************************
Fiesta has been tested on linux platforms running on ARM, AMD and Intel CPUs
with AMD and NVIDIA GPUs.

Fiesta uses Kokkos to compile for various CPU and GPU architectures.  A specific
architecture is enabled through a Kokkos "backend".  Fiesta is known to work
with the CUDA, HIP OPENMP and SERIAL backends.  Kokkos itself supports Intel
GPUs with the SYCL backend but support for Intel GPUs in Fiesta is still under
development.

The following devices have been tested and are known to work:

.. csv-table:: Fiesta Supported CPUs
   :header: "Vendor", "Architecture", "Kokkos Backend"
   :widths: 30, 30, 30

   "Intel", "Broadwell, Skylake, Cascadelake", "Serial,OpenMP"
   "AMD", "Ryzen,Epyc,Rome", "Serial,OpenMP"
   "ARM", "ThunderX2,Apple M1", "Serial,OpenMP"

.. csv-table:: Fiesta Supported GPUs
   :header: "Vendor", "Architecture", "Kokkos Backend"
   :widths: 30, 30, 30

   "NVIDIA", "Kepler, Pascal, Volta, Ampere", "CUDA"
   "AMD", "Instinct", "HIP"

*******************************************************************************
Requirements
*******************************************************************************
Several system packages are required which are not supplied by Fiesta and must be pre-installed.

* | **Compiler** (Required)
  | A C++ compiler which supports C++17 is required such as GNU(g++), LLVM(clang++), and Intel(icpc).
  | GNU g++ 9.3 or greater is recommended.  
  | Chack that g++ version 10 or later is installed with :code:`g++ --version`.

* | **CMake** (Required)
  | CMake version 3.14 or greater is required.
  | Check the cmake version with :code:`cmake --version`

* | **Make Program** (Required)
  | GNU Make or Ninja-build are required.
  | Check that make version 3.82 or greater is installed with :code:`make --version` or ninja version 1.10 or later with :code:`ninja --version`
  | GNU Make is the default on most systems, but the Ninja build tool can reduce compilation time and has an (arguably) more manageable readout.

* | **MPI** (Required if Fiesta will be run on multiple nodes and/or GPUs)
  | MPI may be provided by OpenMPI, MPICH or Intel MPI. OpenMPI is recommended.
  | Check that OpenMPI version 3.1 or greater is installed with :code:`ompi_info --version`.
  | OpenMPI development tools must also be available, check that they are installed by running :code:`mpicc --version`.

* | **CUDA** (Required if Fiesta will be run on NVIDIA GPUs)
  | Check that CUDA version 11.0 or greater is available with :code:`nvcc --version`

* | **HIP** (Required if Fiesta will be run on AMD GPUs)
  | Check that HIP version 4.0 or greater is installed with :code:`hipcc --version`

Third-Party Libraries
===============================================================================
Fiesta makes use of the following third-party libraries.  Unlike the required
packages above, Fiesta can build these packages itself.

* | **Kokkos**
  | `github.com/kokkos/kokkos <https://github.com/kokkos/kokkos>`_
  | Kokkos is a performance portability ecosystem which enables Fiesta to be run on different architectures.
  | Fiesta requires Kokkos version 3.4 or later.

* | **HDF5**
  | `www.hdfgroup.org/solutions/hdf5/ <https://www.hdfgroup.org/solutions/hdf5/>`_
  | The HDF5 library and file format are used for Fiesta's file I/O operations.
  | Fiesta requires HDF5 version 1.8 or later.

* | **Lua**
  | `lua.org <https://lua.org>`_
  | Lua is an embed-able scripting language that is used by Fiesta input decks.
  | Fiesta requires Lua version 5.3 or later.

* | **FMT**
  | `fmt.dev <https://fmt.dev>`_
  | FMT it a C++ string formatting library used for Fiesta logs and command-line output.
  | Fiesta requires FMT version 7.1.3 or later.


*******************************************************************************
Build Instructions
*******************************************************************************
Building Fiesta is done in several steps which are further detailed below.

1. | **Prepare**
   | Acquire the source code and prepare a directory structure.

2. | **Configure**
   | Inform the build system what options Fiesta should be built with and how to satisfy certain dependencies.

3. | **Compile**
   | Compile the Fiesta source code to an executable binary with the indicated configuration options.

4. | **Install** (optional)
   | Optionally install Fiesta and its supporting files to a target location.

Preparation
===============================================================================
Before configuring and installing fiesta, a directory structure must be created
and the source code must be downloaded.  The recommended method for downloading
the source code is with git clone, but pre-packaged versions of the code can be
downloaded without git from `<https://github.com/fiesta-cfd/fiesta/releases>`_.

1. | Create a working directory and clone the Fiesta repository.
   | :code:`mkdir work_dir && cd work_dir` (`work_dir` is an example. Any directory name and location may be chosen)
   | :code:`git clone https://github.com/fiesta-cfd/fiesta`

2. | Create a build directory.  Out-of-tree builds are recommended (where the build
   | directory is outside of the source code directory) but in-tree builds work too.
   | :code:`cd work_dir`
   | :code:`mkdir build && cd build`

3. The resulting directory structure for an out-of-tree build would look like the following:
   ::

      work_dir
      |
      |->fiesta
      |->build

  where `fiesta` contains the fiesta source code and `build` is empty before the
  configuration step.

Configuration
===============================================================================
Configuration is done with the CMake configuration tool.  CMake commands are
executed from the build directory and specify the location of the source code
along with configuration variables.  Configuration variables are set with the
`-D` flag.  For example :code:`cmake /project/source_code -DMY_VARIABLE=val`
will tell cmake to look for the source code at `/project/source` and set the
configuration variable `MY_VARIABLE` to the value `val`.

During configuration the backend must be specified.  The available backends
and their configuration variable are:

* | **Serial**
  | :code:`Fiesta_SERIAL`
  | This is the default option when no backend is specified.  This enables execution on one CPU cores per process.
* | **OpenMP**
  | :code:`Fiesta_OPENMP`
  | The OpenMP backend enables execution on multiple CPU cores per process.
* | **CUDA**
  | :code:`Fiesta_CUDA`
  | The CUDA backend enables execution on NVIDIA GPUs.
* | **HIP**
  | :code:`Fiesta_HIP`
  | The HIP backend enables execution on AMD GPUs.

Fiesta supports three approaches to satisfying third-party library dependencies.

Standard Build
-------------------------------------------------------------------------------
The standard build option will automatically search for the third party
libraries in default locations, or build them if they are not found. The
standard build option is the default and no special configuration variables are
required.

.. warning:: 

    When searching for Kokkos, Fiesta will look for a version which has the
    requested backend.  Configuration will fail with an informative error
    message if a version of Kokkos is found with a wrong backend.  HDF5 also
    requires parallel support.  If HDF5 is found, but without parallel support,
    an informative error message will be displayed.  Reconfigure with on of the
    custom build options below to ignore the pre-installed version of these
    packages.

Superbuild
-------------------------------------------------------------------------------
The superbuild option will compile Fiesta along with several of it's
dependencies including Kokkos, HDF5, Lua and FMT without checking for
pre-installed versions.  This method reduces the dependence on pre-installed
packages but increases compile time.

To configure a super-build, enable the :code:`Fiesta_BUILD_ALL` configuration
variable.

Custom Build
-------------------------------------------------------------------------------
The custom build option allows individual third-party libraries to be built or
found in non-standard locations.

Configuration is done with the CMAKE tool.  Configuration options include
install location, single or multi-node builds and GPU type.  The standard build
does not require any special CMake options.  

There are several configuration options to indicate whether third-party
libraries should be built or where to find them.

* | :code:`Fiesta_BUILD_KOKKOS` Do not look for Kokkos, just build it.
* | :code:`KOKKOS_ROOT` Location in which to look for kokkos.
* | :code:`Fiesta_BUILD_HDF5` Do not look for HDF5, just build it.
* | :code:`HDF5_ROOT` Location in which to look for HDF5.
* | :code:`Fiesta_BUILD_LUA` Do not look for Lua, just build it.
* | :code:`LUA_ROOT` Location in which to look for Lua.
* | :code:`Fiesta_BUILD_FMT` Do not look for FMT, just build it.
* | :code:`FMT_ROOT` Location in which to look for FMT.

Examples
-------------------------------------------------------------------------------
* | Search for pre-installed versions of all dependencies, then build Fiesta with CPU support. (Default)
  | :code:`cmake ../fiesta`

* | Build Fiesta and all it's third-party libraries with support for multiple NVIDIA GPUs.
  | :code:`cmake ../fiesta -DFiesta_BUILD_ALL=on -DFiesta_CUDA=on`

* | Build Fiesta and Kokkos with support for multiple AMD GPUs and search for pre-installed versions of all other third-party libraries.
  | :code:`cmake ../fiesta -DFiesta_BUILD_KOKKOS=on -DFiesta_HIP=on`

* | Build Fiesta with support for NVIDIA GPUs and look for a pre-installed version of Kokkos in a non-standard location.
  | :code:`cmake ../fiesta -DFiesta_CUDA=on -DKOKKOS_ROOT=/path/to/kokkos/install`

* | To use the ninja build generator, use `-G Ninja`:
  | :code:`cmake -G Ninja ../fiesta ...`

After configuration is complete, there are several ways to make changes or corrections:

1. With the `ccmake` tool, if it is installed.  Just run :code:`ccmake .` from the build directory.

2. | by running the cmake command again, but prepended with `--clean-first`.  e.g.:
   | :code:`cmake --clean-first ../fiesta -DFiesta_BUILD_ALL=on -D...`

3. By deleting the contents of the build directory and running another `cmake` command.

Compilation
===============================================================================
After configuration is complete, compile Fiesta by executing :code:`make -j`
from the build directory.  The `-j` option indicates that the maximum number of
CPU cores that are available will be used to speed up compilation. :code:`make
-j 4` would limit the compilation to four CPUs, for example.

Installation
===============================================================================
Fiesta may be run directly from the build directory once it has been compiled.
If installation is necessary (for example, to share the compiled code with
multiple users) then the command :code:`make install` (or `ninja install`) may
be executed.  If a specific installation path was specifies with
`-DCMAKE_INSTALL_PREFIX=/path/to/install`, then Fiesta will be installed there.  If not
installation directory was specified, then Fiesta will attempt to install to the
default location such as `/usr/local`.

*******************************************************************************
Testing
*******************************************************************************
The code correctness of Fiesta can be tested after compilation if the
:code:`-DFiesta_BUILD_TESTS=on` option was provided to `CMake` in the configure
step.  Enabling tests significantly increases compilation time and is most
useful during development of new Fiesta functionality.

After configuring with tests enabled and compiling, run tests from the build
directory with :code:`ctest` to see pass/fail results for each test, or
:code:`ctest -V` to also see detailed output from each test.

Testing is covered in detail in the developer guide.

*******************************************************************************
Configure Reference
*******************************************************************************
All common Fiesta CMake variables are listed below.  Consult the CMake documentation for additional options that are not specific to Fiesta.

* | :code:`Fiesta_BUILD_ALL`
  | This option will ignore pre-installed versions of all third-party libraries and build them from source.  This is equivalent to setting all of :code:`Fiesta_BUILD_KOKKOS=on` :code:`Fiesta_BUILD_HDF5=on` :code:`Fiesta_BUILD_LUA=on` :code:`Fiesta_BUILD_FMT=on`

* | :code:`Fiesta_BUILD_KOKKOS`
  | This option will ignore pre-installed versions and compile Kokkos from source with support for the indicated backend.

* | :code:`Fiesta_BUILD_HDF5`
  | This option will ignore pre-installed versions and compile HDF5 from source with parallel support.

* | :code:`Fiesta_BUILD_LUA`
  | This option will ignore pre-installed versions and compile Lua from source.

* | :code:`Fiesta_BUILD_FMT`
  | This option will ignore pre-installed versions and compile FMT from source.

* | :code:`Fiesta_BUILD_TESTS`
  | This will build ALL unit tests for Fiesta and may significantly impact compile time.

* | :code:`Fiesta_ENABLE_DEBUG`
  | This enables both Fiesta and Kokkos debugging features.

* | :code:`Fiesta_LITE`
  | This option builds Fiesta without distributed memory support (Without MPI).  This option is mostly used for development purposes, but is also useful if Fiesta will only be run on a single GPU.

* | :code:`KOKKOS_ROOT`
  | This option sets the directory in which to look for a Kokkos installation if `Fiesta_BUILD_KOKKOS` is not enabled.

* | :code:`HDF5_ROOT`
  | This option sets the directory in which to look for a HDF5 installation if `Fiesta_BUILD_HDF5` is not enabled.

* | :code:`LUA_ROOT`
  | This option sets the directory in which to look for a Lua installation if `Fiesta_BUILD_LUA` is not enabled.

* | :code:`FMT_ROOT`
  | This option sets the directory in which to look for a FMT installation if `Fiesta_BUILD_FMT` is not enabled.

* | :code:`CMAKE_INSTALL_PREFIX`
  | This option sets the installation directory.

