.. Fiesta documentation master file, created by
   sphinx-quickstart on Tue Feb 11 22:07:50 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

###############################################################################
Tutorials
###############################################################################

*******************************************************************************
Tutorial 1: A 3D Expansion Problem on Xena
*******************************************************************************
This tutorial walks through building Fiesta on the Xena supercomputer at the
University of New Mexico Center for Advanced Research Computing (CARC).  After
building Fiesta, a sample input file is modified in order to execute the
simulation.
    
Problem Description
===============================================================================
This input deck defines an expanding bubble of hot gas. A 0.5 meter diameter
zone of hot gas is defined in the center of the domain with a temperature of
3000 Kelvin and a density of 1000 kg/m^3. The surrounding air has a temperature
of 300 Kelvin and a density of 1 kg/m^3.

Building Fiesta
===============================================================================
After logging into Xena, create a working directory, for example
/users/<username>/fiesta.

.. code-block:: bash

    mkdir fiesta
    cd fiesta

Now clone the Fiesta repository:

.. code-block:: bash

   git clone https://github.com/fiesta-cfd/fiesta

Now create a build directory:

.. code-block:: bash

   mkdir build
   cd build

Now allocate a dual GPU interactive node for 1 hour:

.. code-block:: bash

   salloc -N 1 -p dualGPU -t1:00:00

Now load the required modules:

.. code-block:: bash

   module load gcc/9.3.0-gxji
   module load openmpi/3.1.6
   module load cuda/11.0.2
   module load cmake/3.18.2

Now configure the build:

.. code-block:: bash

    cmake ../fiesta -DFiesta_CUDA=on -DFiesta_BUILD_ALL=on

Next, compile fiesta:

.. code-block:: bash
    
    make -j


Modifying the Input File
===============================================================================

Now that Fiesta has been built, a sample problem can be run.

First, create a directory for running the test.

.. code-block:: bash

   cd ../ # now in /users/<username>/fiesta
   mkdir test
   cd test
   cp ../fiesta/samples/3D_Expansion/fiesta.lua .

Edit the input file to use the number of GPUs available on this system (2).
Change line 40 to::

    procs = {2,1,1}

Save the file and exit the text editor.  Make sure you are in the test
directory: :code:`/users/<username>/fiesta/test`.

Running Fiesta
===============================================================================

The simulation can now be run with the following command.

.. code-block:: bash

   mpirun -n 2 ../build/fiesta -n 2 -c

*******************************************************************************
Tutorial 3: Adding a Buoyancy Term
*******************************************************************************
This tutorial introduces the Fiesta development cycle by walking through the
implementation of a buoyancy term.
