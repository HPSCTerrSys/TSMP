# Structure of DA inside TSMP #

The directory `intf_DA` holds the interfaces of the component models
to data assimilation frameworks (`dart`, `kenda` (currently not in
`master`) and `pdaf`) alongside scripts for building the data
assimilation frameworks.

The directory `intf_DA/pdaf` contains build interface scripts and
the interface code for coupling and running PDAF along with model
component `clm3_5`, `parflow`, `parflow3_0`, `parflow3_2`, `oasis`,
`oasis-mct` and `cosmo4_21`.  It has the following internal structure:

``` text
TSMP/bldsva/intf_DA/pdaf/arch
TSMP/bldsva/intf_DA/pdaf/framework
TSMP/bldsva/intf_DA/pdaf/model
TSMP/bldsva/intf_DA/pdaf/tsmp
```

-   [`arch`](#arch) contains the configuration and build scripts for
    different computers/clusters.
	- Most important file(s): `build_interface_pdaf.ksh`

-   [`framework`](#framework) contains the interface to PDAF. 
	- Most important file: `pdaf_terrsysmp.F90`.

-   [`model`](#model) contains the wrapper/interface to the component
	models.  In the subdirectories, there is component-model-specific
	wrapper functions.
	- Most important files: `wrapper_tsmp.c`, `Makefile`

-   [`tsmp`](#tsmp) contains modified source files for `clm3_5`,
    `parflow`, `parflow3_0`, `parflow3_2`, `oasis`, `oasis-mct` and
    `cosmo4_21`. These source files are copied into the sources of the
    component models, when PDAF is invoked and on top of the files
    copied from `intf_oas3`.
	- Most important file: `Makefile`

## `arch` ##

`arch` currently contains the PDAF build scripts tested for the
following two computers:

``` text
├── JURECA
└── JUWELS
```

The directory looks as follows:

``` text
├── build_interface_pdaf.ksh
└── config
    ├── linux_gfortran_openmpi.h
    └── linux_ifort.h

1 directory, 3 files
```

The build-script `build_interface_pdaf.ksh` contains the following
functions. The script is source and the functions are called in
`/TSMP/bldsva/build_tsmp.ksh` in the function `compileDA()`.

- `always_da()`: 
  - always executed, when DA is used (`$withDA=="true"`) ()
  - up to now: dummy function, no effect
  - `always_*` only used for OASIS3 for setting some directories
- `substitutions_da()`
  - make directory `$dadir/interface` (inside the newly generated
    PDAF-directory-backup with build-suffix)
  - copy directories [`model`](#model) (component model interface) and [`framework`](#framework) (PDAF interface)
    into `interface`
  - make directory `$dadir/lib` (make sure the empty `/lib` exists inside PDAF, this is where the static PDAF library `libpdaf-d.a` is written after compilation)
- `configure_da()`
  - configure three Makefiles (three-step process: (1)copying Makfiles
    into compilation directories, (2) setting variables and then
    overwriting dummy placeholders in the template makefiles), (3)
    `make clean` for all Makefiles)
	- `make.arch/linux_*.h`: PDAF machine-specific configuration
	- `framework/Makefile`
    - `model/Makefile`
- `make_da()`
  - three `make` commands
	- PDAF source `src`
    - PDAF `interface/model`
    - PDAF `interface/framework`
- `setup_da()`


## `framework` ##

## `model` ##

## `tsmp` ##


