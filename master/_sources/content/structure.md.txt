# Structure of TSMP #

After cloning TSMP into a directory, for example `TSMP`, you obtain a fixed top level folder structure. This folder structured is discussed in this section.

## Component model directories and PDAF directory ##

TSMP essentially consists of an interface which couples dedicated versions of the Consortium for Small-scale Modeling (COSMO,
<http://www.cosmo-model.org>) atmospheric model in NWP or climate mode, the Community Land Model (CLM, <http://www.cesm.ucar.edu/models/clm/>), and the hydrologic model ParFlow (<https://www.parflow.org>) through the OASIS3-MCT coupler (<https://portal.enes.org/oasis>, <https://www.mcs.anl.gov/research/projects/mct/>).


### TSMP framework 

In the install directory `TSMP` the `bldsva` can be found, which contain all files of the `TSMP` framework (details below):
``` text
TSMP/bldsva
```

### TSMP component directory 

The directories are called after the respectivel model component and ends with the version number of the model compoennt. 

``` text
TSMP/clm3_5
TSMP/clm4_0
TSMP/cosmo4_21
TSMP/cosmo5_1
TSMP/icon2_1
TSMP/icon2_622
TSMP/parflow
TSMP/parflow3_0
```

The directories `clm3_5`, `clm4_0`, `cosmo4_21`, `cosmo5_1`, `parflow`, and `parflow3_0` contain the model source codes of the respective model components of TSMP.

### OASIS3-MCT and PDAF directory

The[ `bldsva`](#directory-bldsva-and-subdirectories) directory includes the scripts for configuring and compiling TSMP.

The directories `oasis3` and `oasis3-mct` contain the source code of the OASIS3 and OASIS3-MCT coupler, resepectivly. `oasis3`  (without MCT) is deprecated and is no longer used in current versions of TSMP. 

``` text
TSMP/oasis3
TSMP/oasis3-mct
```

After successfully downloading the three parts of TSMP-PDAF (TSMP-repository, PDAF, component models), rename PDAF-source folder
called `PDAF-D_V1.13.2` to `pdaf1_1`.

The directory `pdaf1_1` contains the PDAF source code.

``` text
TSMP/pdaf1_1
```

### Example

The following directory structure should be present in the install directory `TSMP` in case you want to run a fully coupled `TSMP` with COSMO5_1, CLM3_5 and ParFlow:
``` text
TSMP/bldsva
TSMP/clm3_5
TSMP/cosmo5_1
TSMP/parflow
TSMP/oasis3-mct
```

## Directory `bldsva` and subdirectories ##

``` text
TSMP/bldsva/data_oas3
TSMP/bldsva/intf_DA
TSMP/bldsva/intf_oas3
TSMP/bldsva/machines
TSMP/bldsva/setups
```

### data\_oas3

The directory `data_oas3` holds files with static parameters, namely all
variants of the `namcouple` file and `cf_name_table.txt`.

### intf\_oas3

Holds the interfaces of the component models and scripts to build them.
This folder again holds sub-folders with one for each component model
version that is available.

``` text
TSMP/bldsva/intf_oas3/clm3_5
TSMP/bldsva/intf_oas3/clm3_5-icon
TSMP/bldsva/intf_oas3/clm4_0
TSMP/bldsva/intf_oas3/cosmo4_21
TSMP/bldsva/intf_oas3/cosmo5_1
TSMP/bldsva/intf_oas3/icon2-1
TSMP/bldsva/intf_oas3/icon2-622
TSMP/bldsva/intf_oas3/oasis3
TSMP/bldsva/intf_oas3/oasis3-mct
TSMP/bldsva/intf_oas3/parflow
TSMP/bldsva/intf_oas3/parflow3_0
```

Each component model folder contains itself a subset of the following
subfolders, as an example take `clm3_5`:

``` text
TSMP/bldsva/intf_oas3/clm3_5/arch
TSMP/bldsva/intf_oas3/clm3_5/mct
TSMP/bldsva/intf_oas3/clm3_5/oas3
TSMP/bldsva/intf_oas3/clm3_5/pfile  (only clm3_5, possibly deprecated; currently not used)
TSMP/bldsva/intf_oas3/clm3_5/tsmp
```

-   `arch`: This folder contains a sub-folder for every machine on which
    the component model is supported, f.e. `JURECA` and `JUWELS`.

    Each of the machine-folders then contains 2 optional folders:

    1.  `config`, which includes makefiles, configure scripts,
        machine-dependent build files etc, and

    2.  `src` which contains machine specific source code modifications.

-   `tsmp`: This folder holds source code modifications that need to be
    made in the component model for every machine (not only for specific
    machines as the modifications in `/arch/*/src/`).

-   `oas3`: This folder holds source code of the interface between the
    component model and oasis.

-   `mct`: This folder holds source code that is only needed when using
    Oasis3-MCT.

### machines

Holds scripts with machine specifications and also scripts to load
modules among other configurations. This folder has one sub-folder for
all supported machines.

### setups

Holds scripts to specify experiment setups and namelist templates for
this particular setup. This folder has one sub-folder for each supported
setup

### `intf_DA`

See [Structure of DA inside TSMP](./intf_da.md)

## Shell Scripts in `bldsva` ##

The following files reside under `$tsmp_root`/bldsva/:

-   [`build_tsmp.ksh`](#shell-script-build_tsmpksh)

-   [`setup_tsmp.ksh`](#shell-script-setup_tsmpksh)

-   `supported_versions.ksh`

-   `download_data_for_test_test_cases.ksh`

The two important scripts to automatically build TSMP and setup an
experiment are

``` text
build_tsmp.ksh
setup_tsmp.ksh
```

There will be no explicit user documentation, since both scripts include
a man-page styled output by calling the individual script with the
option `–man`.

``` bash
./build_tsmp.ksh --man
./setup_tsmp.ksh --man
```

The shell script `supported_versions.ksh` that contains the supported
versions of the TSMP instance is described in the following
subsection.

The script `download_data_for_test_test_cases.ksh` simplifies working
through the `README.md` of TSMP by automating the download of the
CORDEX test case input data.

### `supported_versions.ksh`

This script is a source list of supported machines and configurations.
It is called by running the main scripts with option `-a` / `–avail`:

``` bash
./build_tsmp.ksh -a
./setup_tsmp.ksh -a
```

The output includes:

-   for `./build_tsmp.ksh`: A list of supported machines including the
    information which version is available on which machine; a list of
    versions including the information of which combinations of
    component models is available for each version and the directory
    names of the component model versions of each TSMP version.

-   for `./setup_tsmp.ksh`: The information from before and additionally
    the available setups for each TSMP-version and a list of setups with
    short documentation. Setups are automatically build test cases for
    TSMP.

The dictionaries containing all information in `supported_versions.ksh`
are implemented as `ksh-associative-array`. These dictionaries are used
by other scripts as well, for example for sanity checks, as possible
user choices in the interactive mode and for naming conventions.

Whenever you pass a machine-name, version, combination, setup, etc. as a
flag or need to know the folder name of a component-model make sure it
is identically written as in this script (case sensitive).

Some array entries are string lists and are converted to a list later
on. Make sure these lists start and end with a space (see comments in
this file).

        # IMPORTANT: add a leading and trailing " "(space)

The spaces are necessary because strings might be parsed with spaces,
i.e. as ` value`. If a default is necessary the first entry of a string
list will be taken.

With versions you have a variety of options to constrain your support.
For example you could create a special version that is only supported on
a certain machine. This version then consists of special
component-model-versions and might only allow certain combinations (for
example only standalone - if MPMD is not supported).

One special naming convention I introduced is for Oasis3-MCT coupling.
Oasis3-MCT is used if \"MCT\" is somehow part of the version-name. (One
could think about adding this information as a flag, but for now I
decided not to do so.)

## Shell Scripts in subdirectories ##

### `intf_oas3/common_build_interface.ksh`

This file has the only intention to reduce duplicate code. It is sourced
by the main-scripts and thus available in all subscripts. I outsourced
code to this file that is the same amongst different versions of a
component-model and amongst different machines. Since many things are
similar between model-versions/machines this lowers the complexity of
adding a new model a lot.

If you add a new machine/model-version and encounter inconsistencies
between this file and your commands you are responsible to solve it.
Either by removing the inconsistent commands from the
`common_build_interface` to the \"real\" interfaces of the form
`intf_oas3/MODEL/arch/ARCH/build_interface_MODEL_ARCH.ksh` (f.e.
`intf_oas3/clm3_5/arch/JUWELS/build_interface_clm3_5_JUWELS.ksh`). If a
lot is inconsistent you can simply not make use of the common interface.

### `intf_oas3/MODEL/arch/ARCH/build_interface_MODEL_ARCH.ksh`

This is the \"real\" interface to compile a specific model on a specific
machine and to handle namelist substitutions.

As you see by the location in the directory tree
(`intf_oas3/MODEL/arch/ARCH/`), there must be a script for every
supported model-platform combination (`MODEL` and `ARCH` ). The script
`build_interface_MODEL_ARCH.ksh` must implement 5 interface routines:

-   `always_MOD()`,

-   `configure_MOD()`,

-   `make_MOD()`,

-   `substitution_MOD()`,

-   `setup_MOD()`, where `MOD`
    $\in \{\mathtt{clm},\mathtt{cos},\mathtt{oas},\mathtt{pfl}\}$

`always_MOD()`: in this routine stuff is handled that always needs to be
done even if compilation is skipped. An example is declaring paths to
oasis libs in [always\_oas()]{.roman}.

`configure_MOD()`: this routine handles everything that is needed to
configure or to create Makefiles. It also edits Makefile-templates after
creation. It is important that also config file-templates are always
copied from the config folder in order to have a fresh template before
editing.

`make_MOD()`: this routine calls `make` and starts compiling and moves
the created executables to a bindir.

`substitution_MOD()`: This routine substitutes source code that is
unchanged for this configuration and is independent from
making/configuring (other than in `configure_MOD()`). This is usually
the whole `oas3` and `tsmp` folder or the `src` folder.

`setup_MOD()`: This routine handles the creation of necessary files for
the run. In especially editing the namelist template. And also execute
the namelist.

The ParFlow `tcl` file and the Cosmo `lmrun_uc` always needed to get
executed. But now also the clm namelist needs an execution (compare with
`lnd.stdin` in nrw setup). This was necessary for `clm4_0` and also for
ensembles.

Most of these routines call a function with the same name but with
leading `c_` for example `c_configure_oas()`. This is the according
routine in the
[`common_build_interface.ksh`](#intf_oas3common_build_interfaceksh)
and handles stuff that is the same throughout versions/machines.

### `machines/ARCH/build_interface_ARCH.ksh`

This script must implement 3 interface routines:

-   `getMachineDefaults()`,

-   `finalizeMachine()`

-   and `createRunscript()`.

In `getMachineDefaults()` some variables need to be defined, basically
loading modules (or executing a load-script) but also the location of
the libraries that are necessary for compiling. Plus some additional
things like default optimization or profiling on this machine (not yet
supported). This routine is used by the interactive mode to pre-load the
values. Because of that this routine might get called several times.

For this reason it is separated from `finalizeMachine()`. It was
supposed to do fixed steps, after a selection was made. Originally it
was planned to load the modules there, but it turned out that this was
mostly also necessary for the first step in order to get library-paths.
It is still used by some machines.

In `createRunscript()` the runscript for a job scheduler on this machine
is created. It calculates processor-distribution and also creates
map-files if needed. It also creates the `instanceMap.txt` which is used
for ensemble runs.

Ensemble runs are currently only supported if Oasis3-MCT is used as
coupler. Ensembles also increase the used processors (as given) by the
factor of ensemble members.

### `setups/SETUP/SETUP_ARCH_setup.ksh`

This script must implement two interface routines:

-   `initSetup()`

-   and `finalizeSetup()`

In `initSetup()` some variables need to be defined, basically values
that are needed to fill the namelist-templates. These might be optional
if they are not needed for the template. But others are mandatory like
namelist-location, forcing-dir etc. Some variables only serve as default
and might get overwritten by the main-scripts (usually start with
\"default\").

This routine is used by the interactive mode to pre-load the values.
Because of that this routine might get called several times.

For this reason it is separated from `finalizeSetup()` where
setup-specific edits to the run-directory is done. For example copying
files like rmp-files or parflow slope-files and other resource files
that are not common to all experiments.

## Shell script `build_tsmp.ksh` ##

This script does everything that is needed to build (create executables)
a given version of TSMP for on certain machine with a certain
model-combination.

### Routines

#### `getDefaults()`

This section serves as hardcoded defaults. All necessary variables are
here defined as a default version `def_...`. They will be used if they
are not overwritten by flags or interactive modifications.

This section is intended to be modified also by users. It is useful if
you are always using a specific setup and don't want to type a lengthy
command every time. Note that by running the script without argument
flags it will start the interactive session. But with providing the flag
`-b` the batch mode is forced. Thus, with `./build_tsmp.ksh -b` only the
hardcoded defaults of this section will be used.

#### `setDefaults()`

This section copies the hardcoded default of the previous session into
the \"real\" data structure/variables.

Some variables need a decision. If they are not hardcoded in the
previous section and not given as flag/interactively a decision is made
in this routine (really hard coded). Most importantly `platform=CLUMA`
and `version=1.1.0MCT`.

#### `clearMachineSelection()`

This is only necessary for handling the interactive mode. If the
platform is switched during the interactive session some variables needs
to be cleared.

#### `clearPathSelection()`

This is only necessary for handling the interactive mode. If some paths
(fore example root) are switched during the interactive session some
variables needs to be cleared.

#### `setSelection()`

If a new machine is chosen (at the beginning or during interactive
session) machine dependent default variables are getting read from
`machines/ARCH/build_interface_ARCH.ksh`.

Note: They will only loaded if this values were not set in any other
selection (like interactive/flag/hardcoded).

#### `finalizeSelection()`

Doing stuff after selections were made. Currently only creating a
bindir.

Questionable if this routine is necessary anymore.

#### `setCombination()`

This routine is setting the `withMOD` variables depending on the chosen
combination. This routine was designed to be flexible with the
combination options. But other routines are hindering this flexibility
anyways.

This routine needs a revision.

#### `compileClm()`

Compiles CLM. It will first source the `build_script` for the selected
machine and then calls the interface. Which interface routines are
called is dependent of the chosen build option (`skip`, `fresh`,
`build`, `configure`, `make`). If `fresh` (default) is chosen, also a
backup is made.

#### `compileCosmo()`

Analog to `compileClm()`

#### `compileOasis()`

Analog to `compileClm()`

#### `compileParflow()`

Analog to `compileClm()`

#### `runCompilation()`

Determines which models are part of the selected combination and above
functions are called.

#### `interactive()`

Handles the interactive session. Reads selections from the console and
pushes return values to the desired variables. Some Variables need extra
handling, in especially if dependencies follow. For example `machine`,
`version`, `rootdir`, etc.

#### `printState()`

Prints all selectable variables. This is used by the interactive mode
and as debug information for the log.

#### `check()`

Helper function to check the `$?` return value of each command.

#### `terminate()`

Clean exit handling.

#### `comment()`

Prints to log and terminal. Usual approach for every system call is to
comment the action and append a `check()` to make sure no error
happened.

#### `route()`

An call that should be made as first and last action of a routine to get
output in which routine the script currently is (debugging).

#### `warning()`

If some sanity checks trigger it is asked for a decision to continue on
own responsibility.

####  `hardSanityCheck()`

Exits if incompatible decisions were made.

####  `softSanityCheck()`

Gives warning if decisions are not supported but could work if it was
taken care of. For example setups/combinations that are not in the
`supported_machine.ksh`.

####  `listAvailabilities()`

Returns a human readable output of supported machines/combinations etc
basically everything that is given in `supported_machines.ksh`. Script
will exit without modifications afterwards.

#### `listTutorial()`

Returns a tutorial. Probably not sufficient enough. Script will exit
without modifications afterwards.

#### `getRoot()`

Determines the root-directory out of pwd and call-command. Working on
usual Linux.

#### `MAIN()`

Starts with determining root and loading defaults. Then handling the
command line flags. Then doing sanity checks and sourcing machine-script
based on flags. Then starting the interactive session. Then finalizing
the selection and starting the compilation. At the end doing some
cleanup and logging.

### Features

#### Layers of defaults and selections

We have two kinds of defaults: Hardcoded defaults from
[`getDefaults`](#getdefaults) (starting with `def_`) and some machine-
and setup specific defaults (starting with `default...`) in the setup-
and machine-interfaces. Then there is the batch-mode which overwrites
existing values with command-line flags and the interactive mode which
will be started if no flags are given.

The hierarchy of overwriting is:

    defaults... < def_ < batch-mode < interactive-mode.

Even if providing flags, the interactive mode can be forced with `-i`
flag. Additional flags will then be pre-selected in the interactive
mode.

#### Logging and debugging

Three log files are created during the build and setup script.

`log_all_$DATE` -- stores all std-output that is generated by system
calls and executing external scripts (like `make` etc.) This output is
not shown while scripts run.

`err_all_$DATE` -- stores all err-output that is generated by system
calls and executing external scripts (like `make` etc.) This output is
not shown while scripts run.

`stdout_all_$DATE` -- stores the output of the script (which you also
see on screen when scripts run)

Usually, a comment output is written to `stdout_*` before system
calls. 

After execution, the system call is usually checked by
[`check()`](#check). If errors are found, the build is terminated.

Additionally, every function entry and exit is logged. 

**Debugging Tip**: If a check fails the script stops and it is
immediately visible in which routine and call it failed. Thus, the
output in `stdout_*` is a good starting point for finding the location
in `TSMP`, where an error occurred during building.
	
</details>


#### Sanity checks

2 types of sanity checks exist, hard and soft.

Hard: The script need to know the machine, version and setup. Otherwise
it does not know what files to source. Thus, if an unsupported
version/machine/setup is given the script will exit.

Soft: Combinations and model versions might work even if not supported/
or under development. An warning will be shown with the option to quit
or continue. Furthermore, during the interactive session only supported
selections will be given as option.

#### Build options

Currently the following build options are implemented to provide
flexibility for debugging/porting/tuning/source modifications.

`skip`: the model is part of the combination, but it is already compiled
and further compilation can be skipped.

`fresh`: the model needs to be build completely and is also backed up to
a new folder. The source folder remains untouched. Also the model source
must lie under `$root/model` (as in the catalog).

`build`: the model needs to be build completely but in the original
folder/or the folder-name that is given to `-wxyz` flags.

`configure`: the model is only configured. New set of
makefile(-templates) make clean configure, etc.

`make`: the model is only compiled (make). No make clean etc., thus, can
resume a started make.

#### Directories

Paths like `bin-dir`, `run-dir`, `forcing-dirs`, `namelist-dirs` are
arbitrary. Even the `root-dir`, which means, that theoretically a
different version of `bldsva` or component-models under `root-dir` could
be used. If this makes sense is questionable (error prone).

## Shell script `setup_tsmp.ksh` ##

### Routines

Most of the routines are identical or have the same purpose as those in
`build_tsmp.ksh`.

The following will only list significantly different.

#### finalizeSelection()

Additionally to the identical stuff it also calculates some processor
counts for the models.

#### softSanityCheck()

Additionally to the identical stuff it also tests for multi-instance
support (is not supported with Oasis3).

#### MAIN()

Starts with determining root and loading defaults. Then handling the
command line flags. Then doing sanity checks and sourcing machine-script
and setup-script based on flags. Then starting the interactive session.
And finalizing the selection afterwards.

Then a loop loops over the number of instances and searches for a
individual namelist for that instance. Then a sub-run-dir is created for
each instance. To create the run the `setup_MOD()` interface routine is
called for every model and executables are copied to the `run-dir`.
After the loop the run-script is created.

At the end doing some cleanup and logging.

### Feature

#### Layers of defaults and selections

See section [Layers of defaults and
selections.](#layers-of-defaults-and-selections).

#### Logging and debugging

See section [Logging and Debugging](#logging-and-debugging).

#### Sanity checks

See section [Sanity checks](#sanity-checks).

#### Directories

See section [Directories](#directories)

For creating run-dirs some extra safety provisions were made. An
experiment-ID (`-I` flag) is appended to the run-dir. This is per
default the date. Also, if rundirs already exist, a backup is created to
avoid potential data losses.

#### Multi instance

Running TSMP in multiple instances is supported by splitting the
`MPI_COMM_WORLD` in sub-communicators inside the models (standalone) or
in Oasis3-MCT (coupled). Inside the run directory sub-directories are
created for every instance. But all start the same executable and will
then `chdir` depending on the rank. The information is given in an
instance-mapfile which is created in `createRunscript()`.

Every instance is looking for an individual namelist which needs to be
named like, for example, `coup_oas.tcl_0` for the 0th member. If no such
namelist is found it will take the base (`coup_oas.tcl`). In this way
not all component models need namelists for every instance if no
perturbation is applied.

#### Restarts

The restarts functionality is basically handled in the namelists itself
in the model-typical way. But two things need to be given to the script
in order to perform a restart.

Restart files: Every restarted model needs a restart file to start from.
This must be provided by the flags `-jkl`.

Initial date: For CLM and ParFlow this only regulates the output naming,
restart files are sufficient for this two. But for Cosmo also the
`start_hours` is determined by `startDate`-`initDate`.

Performing an automated spinup or chain-jobs needs to be handled by the
user (example scripts might be available in the future). This basically
needs a loop over some restart periods and run the setup-script which
creates a specific setup in every iteration.
