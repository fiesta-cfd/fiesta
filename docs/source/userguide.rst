.. highlight:: lua

###############################################################################
User Guide
###############################################################################

Fiesta usage consists of three steps:

1. **Problem Setup** A Lua input file must be created to describe the problem to be
   simulated.

2. **Simulation** The input file formed in step 1 is then used to
   run the simulation.  This may be done interactively or in batch-mode
   depending on the system being run on.  Multiple restarts may be required to
   complete the simulation depending on job length limits.  Simulation will
   result in text output to the screen as well as a time-series of HDF5/XDMF
   solution files pairs and restart files if enabled.

3. **Visualization**  Solution files may then be visualized and post processed
   with standard tools.

*******************************************************************************
Writing Input Files
*******************************************************************************
Fortran input files are written in the Lua language and provide information on
the type of solver, grid and domain decomposition, initial conditions, file I/O,
runtime characteristics and more. Input files are full featured, executable Lua
scripts and are executed when the input file is read.  Fiesta input files may
interact with external Lua scripts and the operating system.

Lua Basics
===============================================================================
The following is a brief introduction to Lua.  Lua is an easily understand
language with a simple syntax.  The minimum language features required to
understand the Fiesta tutorials and samples is presented below with examples of
major language features.

See `<https://www.lua.org/start.html>`_ for more comprehensive resources on
Learning Lua.

Variables
-------------------------------------------------------------------------------
Lua variables or "identifiers" may only begin with uppercase of lowercase
letters or the underscore symbol '`_`' but may also contain digits in the
remainder of the name.  Many parameters that Fiesta reads are set with global
variables in the input file.  All variables are considered global unless
explicitly declares with a different scope.

Variables may take several scalar-like types including number, string, and
boolean and more complex types like table and function Tables and functions are
discussed later.  All number types in Lua are double precision floating point
numbers.

Setting scalar-like type variables (string, boolean and number) is as simple as:::

    nt = 200
    title = "Simulation Case Name"


Comments
-------------------------------------------------------------------------------
Lua scripts may be annotated with inline and block comments:::


    -- This is an inline comment

    --[[ This
         is a multiline
         block comment
    --]]


Tables
-------------------------------------------------------------------------------
Tables in Lua are analogous to arrays in other languages.  Lua tables may
contain a mixture of element types and may even contain other tables.  Tables
are accessed with either a numeric index (starting from 1) or with a key.

The following example creates a table with a string as the first element,
another table as the second element, a number as the third element and a
key/value pair.::

    --[[Create a table with a string as the first element, another table
    as the second element, a number as the third element and a key/value
    pair.--]]

    mytable = {"Apple",{2.1,2.2,2.3},3.14159,height=185}

    --Access the first element with a numerical index
    print(mytable[1]) -- prints "Apple"

    --Access the third element of the table at index 2
    print(mytable[2][3]) -- prints 2.3

    --Key/value pairs can be accessed two ways
    print(mytable["height"]) -- prints 185
    print(mytable.height)    -- prints 185


Arithmetic
-------------------------------------------------------------------------------
All arithmetic in Lua is performed with double precision.::

    -- Add subtract, multiply and divide
    c = c - (b + a)/(a*b)

    -- Modulus can be performed with `%`
    print(7%3) -- prints 1

    -- Exponentiation is done with the `^` operator
    print(2^4) -- prints 16.0

    -- Compute the sine of a number
    print(math.sin(math.pi/2)) -- prints 1.0

Flow Control
-------------------------------------------------------------------------------
Examples of Lua if blocks are below.  Relational operators are:

* :code:`==` Equal
* :code:`~=` Not equal
* :code:`>` Greater than
* :code:`<` Less than
* :code:`>=` Greater than or equal to
* :code:`<=` Less than or equal to

Logical statements can be constructed with `and`, `or` and `not` operands.

If-Else
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
::

    -- inline if statement
    if a<0 then print("a is negative") end

    -- block if/else
    if (a > -1) and (a < 1) then
        print("a is between -1 and 1")
    else
        print("a is not between -1 and 1")
    end

Loops
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Lua supports for loops, while loops and iterators.  ::

    -- while loop (prints the statement infinitely)
    while(true)
    d
        print("How long is forever?")
    end

    -- for loop (prints the digits 1 to 10)
    for i=1,10 do
        print(i)
    end

    -- Print array items and their index with the Pairs iterator
    mytable={"alpha","beta","gamma"}
    for k,v in ipairs(mytable) do
        print(k,v)
    end

Functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Functions in Lua are defined with the :code:`function` keyword below.  Arguments
may be of any type.  A value may be returned with :code:`return`.::

    -- The following function named :code`distance` returns the distance between two points

    function distance(x1,y1,x2,y2)
        return math.sqrt((x2-x1)^2+(y2-y1)^2)
    end


Environment Variables
-------------------------------------------------------------------------------
Lua can read environment variables and execute system calls.

For example to read a task id from a Slurm job array and store it in a variable
called `angle`, do the following.::

    angle = os.getenv("SLURM_ARRAY_TASK_ID")

Input File Parameters
===============================================================================
The following describe the various input file parameters and what they control.

Time
-------------------------------------------------------------------------------
The simulation duration and timestep size are controlled through the :code:`nt`
(number of time-steps) and :code:`dt` (time step size).  Currently,
time-stepping is done with a low-storage second order Runge-Kutta scheme.

By default, the physical simulation time starts from 0.0 and the time-step
indexing starts from zero.  These can be changed with the :code:`time` (start
time) and :code:`tstart` (starting time index).  These are mostly useful for
restarts as discussed below.

Restarts
-------------------------------------------------------------------------------
Fiesta restart files are written at intervals of :code:`restart_frequency` or
when the code receives a :code:`SIGUSR1` signal from the operating system.  If
:code:`restart_frequency=0`, then restart files will not be created unless
Fiesta recieves the :code:`SIGUSR1` signal.  Restart files are names
:code:`restart-####.h5` where `####` is the zero-padded timestep number when the
restart was written.  :code:`restart_path` can be used to specify the location
where restart files are written. By default restarts are written to the same
directory that Fiesta was launched from.

To run Fiesta from a restart file, :code:`restart` must be enabled, and
:code:`restart_name`, :code:`time` and :code`tstart` must be set.::

    -- Write restarts every 100 timesteps to a relative directory names restarts
    restart_frequency = 100
    restart_path = "./restarts"

    -- Restart from a restart file which was written at timestep 100.
    -- Compute the starting time based on the timestep index.
    restart=1
    restart_name=./restarts/restart-0100.h5
    tstart=100
    time=100*dt

Grid and Discretization
-------------------------------------------------------------------------------
Problem dimensions and discretization are controlled through :code:`dx`
(cell size table) :code:`ni` (cell count table) and :code:`grid_type` (grid type).

The discretization of the simulation across multiple GPUs or CPUs is controlled
by :code:`procs` (processor table).  If the number of cells in one direction is
not evenly divisible by the number of processors in that direction, then the
first N processors will have floor(ni[d]/proc[d])+1 cells if ni[d]%proc[d]=N in
direction d, while the remaining processors will have floor(ni[d]/proc[d])
cells.  For example, a problem with 101 cells in the x-direction divided amongst
3 processors will have 34 cells on the first and second processor and 33 cells
on the third.::

    --[[
        The following sets up a domain that is 5mX10mX20m divided into 50
        equally sized hexahedral computational cells in the x-direction, 100
        cells in the y-direction and 200 cells in the z-direction.  This totals
        1 million computational cells.  This computational domain is then
        divided into quadrants and assigned to four processors.
    --]]

    grid_type="cartesian"
    Lx,Ly,Lz=5,10,20    -- User variables storing the domain size
    ni={50,100,200}
    dx={Lx/ni[1],Ly/ni[2],Lz/ni[3]}
    procs={2,2,1}

The "cartesian" grid type uses equally sized rectangular hexahedral cells.  

Generalized curvilinear meshes are also supported.  For generalized meshes, the
:code:`dx` parameter need not be supplied.  Instead, the :code:`grid(i,j,k,d)` function
must be defined.  Third function mus take four parameters: `i,j,k` the index of
the cell corner, starting from 0, `d`, the direction (0 for x, 1 for y, 2 for
z).  This function must return a number indicating the position in the `d`
direction of the indicated corner.

For example, to recreate a Cartesian grid with the generalized coordinate
mechanism:

.. code-block:: lua

    grid_type="generalized"
    Lx,Ly,Lz=5,10,20    -- User variables storing the domain size
    ni={50,100,200}
    dx={Lx/ni[1],Ly/ni[2],Lz/ni[3]}
    function grid(i,j,k,d)
        if d==0 then return dx[1]*i end
        if d==1 then return dx[2]*i end
        if d==2 then return dx[3]*i end
    end
    procs={2,2,1}

Terrain following coordinates can be used with :code:`grid_type="terrain"` and
:code:`terrain_name=/path/to/height/map`.  Terrain following coordinates will
read a height map from a specially formatted HDF5 file.

Initial Conditions
-------------------------------------------------------------------------------
Initial conditions are specified with a function
:code:`initial_conditions(x,y,z,v)`.  This function takes in the coordinates of
the center of a cell and the index of a primary variable and must return the
value of the variable at that location.

In 3D, primary variables are indexed as follows (where N is the number of gas
species):

.. code-block:: lua

    0: x-momentum
    1: y-momentum
    2: z-momentum
    3: specific internal energy
    4: rho1
    5: rho2
    ...
    N+3: rhoN

Boundary Conditions
-------------------------------------------------------------------------------
Boundary conditions are controlled through the :code:`bc*min`, :code:`bc*max`
and :code:`*Per` variables, where `*` is a direction.  Periodic boundary
conditions take precedence over other boundary conditions, therefore
:code:`bcXmin` and :code:`bcXmax` will be ignored if :code:`xPer` is enabled.::
    
    --[[
    Set two reflective walls at the bottom and left of the domain, and set
    periodic conditions in the z-direction.
    --]]

    bcXmin = "reflective"
    bcXmax = "outflow"
    bcYmin = "reflective"
    bcYmax = "outflow"

    zPer = "on"

Gas Species
-------------------------------------------------------------------------------
Gas species must be identified with the :code:`species` table.  The species
table must contain at least one entry.  Each entry must contain four key value
pairs,

* :code:`name` The name of the gas species.

* :code:`gamma` The ratio of specific heats.

* :code:`M` The molecular weight of the species in :math:`\frac{kg}{mol}`.

* :code:`mu` The dynamic viscosity of the species in :math:`\frac{kg}{m\cdot s}`.

::

    --Setup two gas species:  Air and SF6
    species={
        {name="Air", gamma=1.402, M=0.028966, mu=2.928e-5},
        {name="SF6", gamma=4.095, M=0.146060, mu=1.610e-5}
    }

Physics Options
-------------------------------------------------------------------------------
By default, Fiesta solves the Euler equations.  Additional physics can be
enabled through the following variables.

* :code:`viscosity` When enabled computes the viscous stress term.  (Effectively solves the Navier-Stokes equations instead of the Euler equations.

* :code:`buoyancy` Enables gravity driven buoyancy.  A background density of 1.0 kg/m^3 is assumed.

* :code:`cequations` Enables the C-Equations.  The following parameters must also be provided.

  * :code:`kappa` The smoothness coefficient.
  * :code:`epsilon` The support coefficient.
  * :code:`alpha` The isotropic viscosity.
  * :code:`beta`  The anisotropic viscosity.
  * :code:`betae`  The anisotropic viscosity for the energy equation.

* :code:`noise` Enables the wavelet-based noise filter.  The following parameters must also be provided.

    * :code:`n_dh`  Noise indicator cutoff amplitude.
    * :code:`n_eta` Noise indicator attenuation coefficient.
    * :code:`n_coff` Noise indicator C-Equation cutoff amplitude.


Solution I/O
-------------------------------------------------------------------------------
Fiesta solution files can be described with the :code:`blocks` table.  Blocks
describe the location, size, name, resolution and frequency of output of a
region of the flow field. Solution blocks are written in single-precision.

Each block may contain several parameters:

* :code:`name` The name of the output region.

* :code:`path` The location where the solution files will be written for this region.  The default is :code:`./`, or the directory where Fiesta was launched from.

* :code:`frequency` The frequency at which to write output for this region.

* :code:`start` Index coordinates of the beginning of the region. The default is :code:`{0,0,0}`.

* :code:`limit` Index coordinates of the end of the region.  The default is :code:`{ni[1],ni[2],ni[3]}` (the full domain).

* :code:`stride` The stride at which to sample the domain. The default is :code:`{1,1,1}` (full resolution).

* :code:`average` Whether to average the solution by :code:`stride`.  Default is :code:`average=1`.  When enabled, the solution will be averaged over small blocks of size :code:`stride`, however stride blocks are limited to a maximum of :code:`{4,4,4}` when averaging to avoid extra communication.  When :code:`average=0` the solution will be sampled at :code:`stride` intervals and the stride size limit does not apply.

::

    --[[
      Write solution files for two blocks.

      1. A block named "sol" written to the "./solution" directory every 50
         timesteps.  This is a full resolution block covering the entire domain.

      2. A block named "center" writen to the "./centerline" directory every 10
         timesteps.  This block describes a plane one cell thick through the
         center of the domain.  Every 2x2 cells are averaged together to
         decrease the resolution of the output.
    --]]

    blocks={
        {name="sol", path="./solution", frequency=50},
        {name="center", path="./centerline", frequency=10,
         start={0,0,74}, limit={149,149,74}, stride={2,2,1}}
    }

Input File Examples
===============================================================================
Example input files can be found in the samples directory of he source code
tree.

*******************************************************************************
Input File Options Reference
*******************************************************************************
Parameters that take boolean values can be specified in several ways with either
a numerical value (1 or 0) or with a string specifier.  String specifiers are
not case-sensitive.  String specifiers may be surrounded by periods (like in
Fortran `.TRUE.`).  For
example the following are all equivalent:
.. code-block:: lua

  viscosity=1
  viscosity="on"
  viscosity="true"
  viscosity=".TRUE."
  viscosity="enable"
  viscosity="enabled"

Similarly, the following are also equivalent
.. code-block:: lua

  buoyancy=0
  buoyancy="off"
  buoyancy="false"
  buoyancy=".FALSE."
  buoyancy="disable"
  buoyancy="disabled"

Some parameters are defaultable.  If they are not included in the input file,
their default value is used.

* | :code:`title`
  | *string*
  | **Required**
  | An informative title for the simulation.

* | :code:`nt`
  | *integer greater than zero*
  | **Required**
  | The total number of timesteps to compute.

* | :code:`dt` [seconds]
  | *floating point number greater than zero*
  | **Required**
  | The time step.

* | :code:`R` [J/(mol**K)]
  | *floating point number greater than zero*
  | **Default: 8.314462**
  | The universal gas constant

* | :code:`tstart`
  | *integer greater than zero*
  | **Default: 0**
  | Starting time index used for restarts.  

* | :code:`time` [seconds]
  | *floating point number greater than zero*
  | **Default: 0.0**
  | The starting time of the simulation.  Used primarily for restarts.

* | :code:`cequations`
  | *boolean*
  | **Default: disabled**
  | Enable C-Equations.

* | :code:`noise`
  | *boolean*
  | **Default: disabled**
  | Enable noise indicator.

* | :code:`buoyancy`
  | *boolean*
  | **Default: disabled**
  | Enable buoyancy.

* | :code:`viscosity`
  | *boolean*
  | **Default: disabled**
  | Enable viscosity.

* | :code:`ndim`
  | *2 or 3*
  | **Default: 2**
  | Dimensionality of problem.

* | :code:`xPer`
  | *boolean*
  | **Default: false**
  | Periodicity in the x-direction.  Periodic boundaries have higher precedence than other boundary conditions.

* | :code:`yPer`
  | *boolean*
  | **Default: false**
  | Periodicity in the y-direction.  Periodic boundaries have higher precedence than other boundary conditions.

* | :code:`zPer`
  | *boolean*
  | **Default: false**
  | Periodicity in the z-direction.  Periodic boundaries have higher precedence than other boundary conditions.

* | :code:`bcXmin`
  | *0: Freeflow*
  | *1: Reflective*
  | *2: No-slip*
  | *3: Hydrostatic*
  | **Default: 0**
  | Boundary condition type for the minimum x face.

* | :code:`bcXmax`
  | *0: Freeflow*
  | *1: Reflective*
  | *2: No-slip*
  | *3: Hydrostatic*
  | **Default: 0**
  | Boundary condition type for the maximum x face.

* | :code:`restart`
  | *boolean*
  | **Default: 0**
  | Enable starting the simulation from a restart file.
  | Requires setting :code:`restart_name`, :code:`tstart`, and :code:`time`.

* | :code:`restart_name`
  | *string*
  | **Required if** :code:`restart` **is enabled**
  | The name and location of the restart file with which to start the simulation if restarts are enabled.

* | :code:`restart_path`
  | *string*
  | **Default: ./**
  | The location where restart files will be written.

* | :code:`terrain_name`
  | *string*
  | **Default: terrain.h5**
  | Location and name of terrain data file.

* | :code:`progress_frequency`
  | *integer greater than zero*
  | **Default: 0**
  | Frequency at which to write out messages showing current time step number and simulation elapsed time.
  | A value of `0` disables these messages.

* | :code:`write_frequency`
  | *integer greater than zero*
  | **Default: 0**
  | Frequency at which to write default solution files.  Default solution files are full resolution, single precision files with primary and secondary flow variables.
  | A value of `0` disables writing of default solution files.

* | :code:`restart_frequency`
  | *integer greater than zero*
  | **Default: 0**
  | Frequency at which restart files are written.  Restart files may be very large and can impact performance if written high frequencies.
  | A value of `0` disables these restart dumps.

* | :code:`status_frequency`
  | *integer greater than zero*
  | **Default: 0**
  | Frequency at which to report simulation status.  Status reports consist of estimated time to completion and minimum and maximum values of flow variables. Status reports have a small impact on performance and which may become significant at high frequencies.
  | A value of `0` disables these status reports.

* | :code:`advection_scheme`
  | *'weno5', 'quick', 'centered4'*
  | **Default: 'weno5'**
  | Which scheme to use for advective terms.

* | :code:`grid_type`
  | *'cartesian', 'generalized', 'terrain'*
  | **Default: 'cartesian'**
  | Grid type.

* | :code:`dx`
  | *2 or 3 element array*
  | **Required**
  | Array of cell sizes.

* | :code:`ni`
  | *2 or 3 element array*
  | **Required**
  | Array of cell counts describing the problem size.

* | :code:`procs`
  | *2 or 3 element array*
  | **Required**
  | Array of process counts in each direction.

* | :code:`kappa`
  | *floating point number*
  | **Required if** :code:`cequation` **is enabled**
  | C-Equation parameter

* | :code:`epsilon`
  | *floating point number*
  | **Required if** :code:`cequation` **is enabled**
  | C-Equation parameter

* | :code:`alpha`
  | *floating point number*
  | **Required if** :code:`cequation` **is enabled**
  | C-Equation parameter

* | :code:`beta`
  | *floating point number*
  | **Required if** :code:`cequation` **is enabled**
  | C-Equation parameter

* | :code:`betae`
  | *floating point number*
  | **Required if** :code:`cequation` **is enabled**
  | C-Equation parameter

* | :code:`n_dh`
  | *floating point number*
  | **Required if** :code:`noise` **is enabled**
  | Noise indicator cutoff amplitude.

* | :code:`n_eta`
  | *floating point number*
  | **Required if** :code:`noise` **is enabled**
  | Noise indicator attenuation coefficient.

* | :code:`n_coff`
  | *floating point number*
  | **Required if** :code:`noise` **is enabled**
  | Noise indicator C-Equation cutoff amplitude.

* | :code:`block`
  | *table of I/O block specifiers*
  | **Optional** :code:`noise` **is enabled**
  | Table of parameters specifying I/O blocks for solution date output.  See `Output and Visualization`_ for more information.

* | :code:`species`
  | *table of species names and data*
  | **Minimum one entry required**
  | Table of parameters specifying a species name and physical properties.  See `Gas Species`_ for more information.


*******************************************************************************
Running Fiesta
*******************************************************************************

Once Fiesta has been built and a problem has been setup through a Lua input
file, then the simulation may be executed.

Fiesta can be run in interactive mode and in batch-mode.  Batch-mode runs depend
on the scheduling software installed on the cluster and is not covered here.
See your cluster documentation for setting up batch jobs.

In interactive mode, Fiesta can be executed on a system with four GPUs with
color text output with the following command.

.. code-block:: bash

   mpirun -n /path/to/fiesta -c

By default fiesta looks for an input file names :code:`fiesta.lua` in the
current directory. An input file with a different name can be specified on the
command line.

.. code-block:: bash

   mpirun -n /path/to/fiesta /path/to/input.lua -c

The verbosity level can be changed to show more or less information.  By default
Fiesta will output all errors, warnings and messages except for debugging
messages.  See the command line reference for a description of the different
verbosity levels.

The version of Fiesta can be checked with the :code:`--version` flag.  This will
print out the current version and configuration settings of Fiesta and exit
without reading any input files or running any simulation.  During execution of
a simulation, this version information is also included in the standard output.

Command Line Options Reference
===============================================================================

* | :code:`-n <N>`
  | This flag indicates the number of GPUs per node for GPU builds or the number of OpenMP threads per rank for OpenMP builds.

* | :code:`-c, --color[=when]`
  | Colorize output.  Default: on
  | on: always colorize output
  | off: never colorize output
  | auto: colorize output only if the detected output device supports color

* | :code:`-v, --verbosity[=level]`
  | Set verbosity level. Default: 4
  | 1: Errors only
  | 2: Errors and warnings
  | 3: Errors, warnings and status messages
  | 4: Errors, warnings, status messages and informational messages
  | 5: Errors, warnings, status messages, informational messages and debug messages

* | :code:`-V, --version`
  | Display version and build information.

*******************************************************************************
Output and Visualization
*******************************************************************************

For each I/O block, fiesta outputs a time series of HDF5/XDMF file pairs.  The
XDMF file is used to describe the layout of the HDF5 file to visualization
software.  The XDMF files can be opened in ParaView, Tecplot, VisIt or any other
visualization package that supports XDMF.

The HDF5 files may also be post processed through various scripting languages.
There are examples of post processing scripts in the samples directory using
Python and Julia.

