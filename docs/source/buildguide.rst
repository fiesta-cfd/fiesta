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
Building
*******************************************************************************
Fiesta supports three build options for compiling and installing the code.

1. **Super Build** The super build option will force Fiesta to build some of it's
   third party dependencies.  This option is best for sites which do not provide
   pre-installed version of these dependencies.  The super build option is
   simpler and requires less user interaction, however, it increases the
   compilation time and does not set special comilation options that may be
   required at a site.

2. **Standard Build** The Site-Build option first looks for properly installed
   versions of Kokkos, HDF5, Lua and FMT.  If it finds them, it will use them,
   if not, then they will be built.  The Site-build option results in reduced
   compilation times and can take advantage of site-specific configurations.

3. **Custom Build** The Custom-Build option allows control over which packages
   are built and which are searched for.

4. **Containers** Some prebuilt containers images are available for Fiesta.
   Container usage is an advanced option and their management is beyond the
   scope of this document.  However, the available containers are listed below.

Requirements
-------------------------------------------------------------------------------
Regardless of build type, several packages are required which are not supplied
by Fiesta and must be pre-installed.

* | Compiler (Required)
  | A C++ compiler which supports C++17 is required such as GNU(g++), LLVM(clang++), and Intel(icpc).
  | GNU g++ 9.3 or greater is recommended.  
  | Chack that g++ version 10 or later is installed with :code:`g++ --version`.

* | CMake (Required)
  | CMake version 3.12 or greater is required.
  | Check the cmake version with :code:`cmake --version`

* | Make Program (Required)
  | GNU Make or Ninja-build are required.
  | Check that make version 3.82 or greater is installed with :code:`make --version` or ninja version 1.10 or later with :code:`ninja --version`

* | MPI (Required if Fiesta will be run on multiple nodes and/or GPUs)
  | MPI may be provided by OpenMPI, MPICH or Intel MPI. OpenMPI is recommended.
  | Check that OpenMPI version 3.1 or greater is installed with :code:`ompi_info --version`.
  | OpenMPI development tools must also be available, check that they are installed by running :code:`mpicc --version`.

* | CUDA (Required if Fiesta will be run on NVIDIA GPUs)
  | Check that CUDA version 11.0 or greater is available with :code:`nvcc --version`

* | HIP (Required if Fiesta will be run on AMD GPUs)
  | Check that HIP version 4.0 or greater is installed with :code:`hipcc --version`

Third-Party Libraries
-------------------------------------------------------------------------------
Fiesta makes use of the following third-party libraries.  The superbuild option
will provide all these packages, but pre-installed versions may be used with the
standard and custom builds.

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
  | Lua is an embed-able scripting langage that is used by Fiesta input decks.
  | Fiesta requires Lua version 5.3 or later.

* | **FMT**
  | `fmt.dev <https://fmt.dev>`_
  | FMT it a C++ string formatting library used for Fiesta logs and command-line output.
  | Fiesta requires FMT version 7.1.3 or later.

===============================================================================
1. Superbuild
===============================================================================
The easiest way to install Fiesta is with the superbuild option.  This option
will install Fiesta along with several of it's dependencies including Kokkos,
HDF5, Lua and FMT.  This method reduces the dependence on pre-installed packages.

Preparation
-------------------------------------------------------------------------------
Before configuring and installing fiesta, a directory structure must be created
and the source code must be downloaded.  The recommended method for downloading
the source code is with git clone, but pre-packaged versions of the code can be
downloaded wirhout git from `<https://github.com/fiesta-cfd/fiesta/releases>`_.

1. | Create a working directory and clone the Fiesta repository.
   | :code:`mkdir work_dir && cd work_dir` (`work_dir` is an example. Any directory name and location may be chosen)
   | :code:`git clone https://github.com/fiesta-cfd/fiesta`

2. | Create a build directory.  Out-of-tree builds are recommended (where the build
   | directory is outside of the source code directory) but in-tree builds work too.
   | :code:`cd work_dir`
   | :code:`mkdir build && cd build`

3. The resulting directory structure should look something like this:
   ::

      work_dir
      |
      |->fiesta
      |->build

  where `fiesta` contains the fiesta source code and `build` is empty before the
  configuration step.

Configuration
-------------------------------------------------------------------------------
Fiesta may now be configured.  Configuration is done with the CMAKE tool.
Configuration options include install location, single or multi-node builds and
GPU type.  The superbuild uses the special option `Fiesta_BUILD_ALL` as
indicated below.

To configure Fiesta for different situations, run one of the following configure
commands from the build directory.

* | For multi-gpu support for NVIDIA GPUs:
  | :code:`cmake ../fiesta -DFiesta_BUILD_ALL=on -DFiesta_CUDA=on`

* | For multi-gpu support for AMD GPUs:
  | :code:`cmake ../fiesta -DFiesta_BUILD_ALL=on -DFiesta_HIP=on`

* | For single-cpu, serial support (also known as a lite build):
  | :code:`cmake ../fiesta -DFiesta_BUILD_ALL=on -DLITE=on`

* | For single-gpu support:
  | :code:`cmake ../fiesta -DFiesta_BUILD_ALL=on -DLITE=on -DFiesta_CUDA=on`

* | For multi-node multi-cpu support:
  | :code:`cmake ../fiesta -DFiesta_BUILD_ALL=on`

* | For single-node multi-cpu support:
  | :code:`cmake ../fiesta -DFiesta_BUILD_ALL=on -DFiesta_OPENMP=on`

* | To specify an install directory, use `-DCMAKE_INSTALL_PREFIX`. e.g.:
  | :code:`cmake ../fiesta -DFiesta_BUILD_ALL=on -DFiesta_CUDA=on -DCMAKE_INSTALL_PREFIX=~/.local`
  | This configures Fiesta for multi-gpu support for NVIDIA gpus and specifies that the code should be installed to the users local bin directory.

* | To use the ninja build generator, use `-G`:
  | :code:`cmake -G Ninja ../fiesta -DFiesta_BUILD_ALL=on ...`

After configuration is complete, there are several ways to make changes or corrections:

1. With the `ccmake` tool, if it is installed.  Just run :code:`ccmake .` from the build directory.

2. | by running the cmake command again, but prepended with `--clean-first`.  e.g.:
   | :code:`cmake --clean-first ../fiesta -DFiesta_BUILD_ALL=on -D...`

3. By deleting the contents of the build directory and running another `cmake` command.

Installation
-------------------------------------------------------------------------------
Fiesta may now be compiled and installed.

To compile with GNU Make (the default, if Ninja was not explicitly specified) run:
:code:`make -j`
The `-j` option indicates that as many CPU cores as are available will be used
to speed up compilation.  :code:`make -j 4` would limit the compilation to four
CPUs for example.

If Ninja was explicitly specified in the configure step, then, from the build directory, run:
:code:`ninja` to build in parallel or :code:`ninja -j N` to limit the
compilation to `N` CPU cores.

If a specific installation path was specifies with `-DCMAKE_INSTALL_PREFIX`,
then fiesta can now be installed with :code:`make install` of :code:`ninja
install`.


===============================================================================
2. Standard Build
===============================================================================
The standard build option will automatically search for the third party
libraries, then build them if they are not found.

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
The functionality and correctness of Fiesta can be tested after compilation if
the :code:`-DFiesta_BUILD_TESTS=on` option was privided to `CMake` in the
configure step.  Enabling tests significantly increases compilation time and is
most useful during development of new Fiesta functionality.

After configuring with tests enabled and compiling, run tests from the build
directory with :code:`ctest` to see pass/fail results for each test, or
:code:`ctest -V` to also see detailed output from each test.

Testing is covered in detail in the developer guide.

*******************************************************************************
Configure Reference
*******************************************************************************

* | :code:`Fiesta_BUILD_ALL`
  | Build all third-party-libraries (tpls).  This is equivalent to setting all of :code:`Fiesta_BUILD_KOKKOS=on` :code:`Fiesta_BUILD_HDF5=on` :code:`Fiesta_BUILD_LUA=on` :code:`Fiesta_BUILD_FMT=on`

