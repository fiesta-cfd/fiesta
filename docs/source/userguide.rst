###############################################################################
User Guide
###############################################################################

*******************************************************************************
Writing Input Files
*******************************************************************************
Fortran input files are written in the Lua language and provide information on
the type of solver, grid and domain decomposition, initial conditions, file I/O,
runtime characteristics and more. Input files are full featured Lua scritps and
are executed when the input file is read.  Fiesta input files may interact with
external Lua scripts and the operating system through system commands.

Lua Basics
===============================================================================
What follows is a "crash-course" in Lua.  Lua is an easy to understand language
with a simple syntax.  The minimum language features required to understand with
the Fiesta tutorials and sample input files is presented below mostly with
examples and brief explanations.

See the resources at `<https://www.lua.org/start.html>`_ for more comprehensive
resources on Learning Lua.

Variables
-------------------------------------------------------------------------------
Lua variables of "identifiers" may only begin with uppercase of lowercase
letters or the undescore symbol '`_`' but may also contain digits in the
remainder of the name.  Many parameters that Fiesta reads are set with global
variables in the input file.  All variables are considered global unless
explicitly declares with a different scope.

Variables may take several scalar types including number, string, and boolean
and more complex types like table and function Tables and functions are
discussed later.  All number types in Lua are double precision floating point
numbers.

Setting scalar type variables (string, boolean and number) is as simple as:
::

    nt = 200
    title = "Simulation Case Name"


Comments
-------------------------------------------------------------------------------
Lua supports inline and block comments:

::

    -- This is an inline comment

::

    --[[ This
         is a multiline
         block comment
    --]]


Tables
-------------------------------------------------------------------------------
Tables in Lua are analogous to arrays in other languages.  Lua tables may
contain a mixture of element types and may even contain other tables.  Tables
are accessed with either a numeric index (starting from 1) or with a key.

Create a table with a string as the first element, another table
as the second element, a number as the third element and a key/value pair.

::

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
All arithmetic in Lua is performed with double precision.

::

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
* :code:`<=` Lesss than or equal to

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
Lua supports for loops, while loops and iterators.

::
    
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


Environment Variables
-------------------------------------------------------------------------------
Lua can read environment variables and execute system calls.

For example to read a task id from a slurm job array and store it in a variable
called `angle`, do the following.

::

    angle = os.getenv("SLURM_ARRAY_TASK_ID")


Time
===============================================================================
The simulation duration and timestep size are controlled through the :code:`nt`
(number of time-steps) and :code:`dt` (time step size).  Currently,
time-stepping is done with a low-storage second order Runge-Kutta scheme.

by default, the physical simulation time starts from 0.0 and the time-step
indexing starts from zero.  These can be changed with the :code:`time` (start
time) and :code:`tstart` (starting time index).  These are mostly useful for
restarts as discussed below.

Restarts
===============================================================================

Grid and Discretization
===============================================================================
Problem dimensions and discretization are controlled through :code:`dx`
(cell size table) :code:`ni` (cell count table) and :code:`grid_type` (grid type).

The discretization of the simulation across multiple GPUs or CPUs is controlled
by :code:`procs` (processor table).  If the number of cells in one direction is
not evenly divisible by the number of processors in that direction, then the
first N processors will have floor(ni[d]/proc[d])+1 cells if ni[d]%proc[d]=N in
direction d, while the remaining processors will have floor(ni[d]/proc[d])
cells.  For example, a problem with 101 cells in the x-direction divided amongst
3 processors will have 34 cells on the first and second processor and 33 cells
on the third.

Example:
::

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

Genarlized curvilinear meshes are also supported.  For generalized meshes, the
:code:`dx` parameter need not be supplied.  Instead, the :code:`grid(i,j,k,d)` function
must be defined.  Thid function mus take four parameters: `i,j,k` the index of
the cell corner, starting from 0, `d`, the direction (0 for x, 1 for y, 2 for
z).  This function must return a number indicating the position in the `d`
direction of the indicated corner.

For example, to recreate a cartesian grid with the generalized coordinate
mechanism:
::

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

InitialConditions
===============================================================================
Initial conditions are specified with a function
:code:`initial_conditions(x,y,z,v)`.  This function takes in the coordinates of
the center of a cell and the index of a primary variable and must return the
value of the variable at that location.

In 3D, primary variables are indexed as follows (where N is the number of gas
species):
::

    0: x-momentum
    1: y-momentum
    2: z-momentum
    3: specific internal energy
    4: rho1
    5: rho2
    ...
    N+3: rhoN

Boundary Conditions
===============================================================================

Gas Species
===============================================================================

Solution I/O
===============================================================================

Input File Options Reference
===============================================================================
Parameters that take boolean values can be specified in several ways with either
a numerical value (1 or 0) or with a string specifier.  Srting specifiers are
not case-sensitive.  String specifiers may be surrounded by periods (like in
fortran `.TRUE.`).  For
example the following are all equivalent:
::

  viscosity=1
  viscosity="on"
  viscosity="true"
  viscosity=".TRUE."
  viscosity="enable"
  viscosity="enabled"

Similariy, the following are also equivalent
::

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
  | Enable cequations.

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
  | Frequency at which to report simulation status.  Status reports concist of estimated time to completion and minimum and maximum values of flow variables. Status reports have a small impact on performance and which may become significant at high frequencies.
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
  | Cequation parameter

* | :code:`epsilon`
  | *floating point number*
  | **Required if** :code:`cequation` **is enabled**
  | Cequation parameter

* | :code:`alpha`
  | *floating point number*
  | **Required if** :code:`cequation` **is enabled**
  | Cequation parameter

* | :code:`beta`
  | *floating point number*
  | **Required if** :code:`cequation` **is enabled**
  | Cequation parameter

* | :code:`betae`
  | *floating point number*
  | **Required if** :code:`cequation` **is enabled**
  | Cequation parameter

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
  | Noise indicator c-equation cutoff amplitude.

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

Command Line Options
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

Examples
===============================================================================

*******************************************************************************
Output and Visualization
*******************************************************************************

