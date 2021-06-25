###############################################################################
User Guide
###############################################################################

*******************************************************************************
Command Line Options
*******************************************************************************
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
Input File
*******************************************************************************

Options
===============================================================================
Parameters that take boolean values can be specified in several ways with either
a numerical value (1 or 0) or with a string specifier.  Srting specifiers are
not case-sensitive.  String specifiers may be surrounded by periods (like in
fortram `.TRUE.`).  For
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
  | Table of parameters specifying I/O blocks for solution date output.  See Solution I/O section for more information.

Grid
===============================================================================

InitialConditions
===============================================================================

Boundary Conditions
===============================================================================

Lua Basics
===============================================================================

*******************************************************************************
Output and Visualization
*******************************************************************************

Supported Formats
===============================================================================

ASCII Output
===============================================================================

*******************************************************************************
Examples
*******************************************************************************

*******************************************************************************
Issues
*******************************************************************************

Documentation for Fiesta

.. math::
    (a + b)^2  &=  (a + b)(a + b) \\
               &=  a^2 + 2ab + b^2


Input options are as follows

.. .. tabularcolumns:: |l|l|l|

+-------+-------+----------------------------------------+
|Option |Default|Description                             |
+=======+=======+========================================+
|visc   |0      |Enable(1) or Disable(0) the viscous term|
+-------+-------+----------------------------------------+
|ceq    |1      |Enable(1) or Disable(0) the C-Equations |
+-------+-------+----------------------------------------+

visc
ceq
