Logger
======

Prints timestamped log messages to a file.  A program logger is created near the
beginning of execution and is stored in the configuration object ``cf.log``.
Additional logger objects can be created if needed.  See the logger constructor
documentation below. Timestamps are in seconds since the logger was created. For
the program logger ``cf.log`` this is essentially the number of seconds since
program execution began.

Colorized logs can be viewed with the command ``less -R fiesta.log``.
Sufficiently modern versions of ``cat`` will also interpret color.

Command Line Options
--------------------

* ``-v [level]`` (lower case v) Verbosity Level

  1. Errors
  2. Errors and Warnings
  3. Errors, Warnings and Messages
  4. Errors, Warnings, Messages and Info
  5. Errors, Warnings, Messages, Info and Debug Messages

* ``-l`` (lower case L) Colorize Log Output
* ``--log-file-name=[filename]`` Name of file with logger output.
  Default: ``fiesta.log``

Usage
-----

Example (in code):
::

  cf.log->error("This is an error ","with data: ",42);
  cf.log->warning("This is a warning ","with data: ",3.14," And more date: ",159);
  cf.log->message("This is a message ","with data: ",0xAF);
  cf.log->info("This is info ","with data: ",1.675);
  cf.log->debug("This is a debug message ","with data: ",42);

Output (contents of ``fiesta.log``):
::

  [       2.5395367]    Error: This is an error with data: 42
  [       2.5395501]  Warning: This is a warning with data: 3.14 And more date: 159
  [       2.5395559]  Message: This is a message with data:175
  [       2.5395617]     Info: This is info with data: 1.675
  [       2.5395666]    Debug: This is a debug message with data: 42


When colorized log files are enabled (with ``-l``), the following nominal colors are specified.
The actual colors used will depend on your terminal configuration.

* Debug: magenta
* Info: none(white)
* Message: blue
* Warning: yellow
* Error: red
  
Logger Constructor
------------------

Additional logger objects can be created if necessary.  They will not obey the
command line options above unless the appropriate values of ``cArgs`` are used
in the constructor.
::

  Logger myLogger = Logger(int verbosity, int colorMode, int colorLogs, int mpiRank, string logName);

* ``verbosity`` Verbosity level, see level descriptions above.  Use ``cArgs.verbosity`` to retrieve
  the log level requested on the command line.  
* ``colorMode`` Colormode to use if ``colorLogs`` is true.

  * ``1`` Force Color
  * ``2`` Use Color only if output it to a tty(e.g. /dev/s
  * Use ``cArgs.colorFlag`` to
    retrieve the colormode requested on the command line.tdout).

* ``colorLogs`` Whether or not to colorize logs.  Use ``cArgs.colorLogs`` to
  retrieve the logger color flag requested on the command line.
* ``mpiRank`` Rank of current process.  Only rank ``0`` will actually write to
  log files.
* ``logName`` Name of file where logs will be written.

