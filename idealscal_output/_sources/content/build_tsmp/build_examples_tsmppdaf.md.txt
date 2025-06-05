# TSMP-PDAF Build Examples #

A collection of build examples for TSMP-PDAF.

All builds involving more than one coupled component model use the
coupling software Oasis3-MCT. Builds in the TSMP-PDAF branch follows the syntax
shown here. Building PDAF in the master branch has a syntax without the internal
 version number:
 			
 			./build_tsmp.ksh -m JUWELS -c clm3-cos5-pfl-pdaf -O Intel

- [`Compile fully coupled TSMP-PDAF`](#compile-fully-coupled-tsmp-pdaf-on-jureca)
- [`Compile ParFlow standalone with PDAF`](#compile-parflow-and-clm-with-pdaf)
- [`Compile ParFlow and CLM with PDAF`](#compile-parflow-standalone-with-pdaf-on-juwels-without-compiler-optimization)

## Compile fully coupled TSMP-PDAF (on JURECA) ##

The first example runs the build script of TSMP-PDAF, specifies the
machine `JURECA`, specifies the component models CLM, COSMO and
ParFlow and specifies version `1.1.0MCTPDAF`.

      ./build_tsmp.ksh -m JURECA -c clm-cos-pfl -v 1.1.0MCTPDAF

The `-c` and `-v` flag are optional in this case because the inputs for
these flags on `-m JURECA` are equivalent to the corresponding default
inputs.

The `-v` flag is optional in this case because `-c clm-cos-pfl` is the
default input for this flag on `-m JURECA` and for version `-v
1.1.0MCTPDAF`.

By default the `build_tsmp.ksh` script will make a copy of the original
model source folder (named as in the catalog, f.e. `clm3_5`,
`cosmo4_21`, and `parflow`) and back it up to
`MODEL_ARCH_VERSION_COMBINATION`.

With the flags `-wxyz` (see man-page) you can specify your own directory
where it is backed up to. If you like to use your own folder structure
for the models you can again select your own model-path with `-wxyz` and
avoid making a backup with `-WXYZ` option set to \"`build`\".

## Compile Parflow and CLM with PDAF ##

### Build commands ###

Build commands for different machines.

Remarks:
- For building TSMP-PDAF with DA, use branch `TSMP_pdaf`
- Parflow-3.9 is supported using version tag `3.1.0MCTPDAF`. For older
  versions, use `3.0.0MCTPDAF` (>=3.2, <3.7), `1.1.0MCTPDAF` (<3.2).

#### JURECA ####

      ./build_tsmp.ksh -m JURECA -c clm-pfl -v 3.1.0MCTPDAF -O Intel

#### JUWELS ####

      ./build_tsmp.ksh -m JUWELS -c clm-pfl -v 3.1.0MCTPDAF -O Intel

#### Remote build ####

``` shell
	./build_tsmp_remote -m JUWELS -v 3.1.0MCTPDAF -O Intel -c clm-pfl -branch TSMP_pdaf -host juwels &
	./build_tsmp_remote -m JURECA -v 3.1.0MCTPDAF -O Intel -c clm-pfl -branch TSMP_pdaf -host jureca &
```

### Component Models ###

The build commands were last tested for these component models

- **TSMP**
  - current revision of branch `TSMP_pdaf` from
    <https://github.com/HPSCTerrSys/TSMP>

- **clm3\_5**
  - Version 3.5, `clm3_5/Copyright`, `share3_070321` in
    `clm3_5/src/csm_share/ChangeLog`
    <https://icg4geo.icg.kfa-juelich.de/ModelSystems/tsmp_src/clm3.5_fresh>
    revision `801b5304179f0a8cbe3dc2c50b584a6bfee387b0`
- **oasis3-mct**
  - (svn revision r1506, `svn info` (caused error);
  - new: Gitlab repo
	<https://icg4geo.icg.kfa-juelich.de/ModelSystems/tsmp_src/oasis3-mct.git>,
	branch `oasis3-MCT_2.0`, revision `bc58342`)
- **parflow3\_0** (`1.0.0MCTPDAF`)
  - Version v3.2.0 in `VERSION`, no `CMakeLists.txt` yet, revision
    `faaee2c` in `parflow_vTdz` from
    <https://git.meteo.uni-bonn.de/git/parflow>
- **parflow3\_2** (`3.0.0MCTPDAF`)
  - Version v3.2.0 in `VERSION`, no `CMakeLists.txt` yet, revision
    `98a87011` in `master` from
    <https://github.com/parflow/parflow.git>
- **parflow** (`3.1.0MCTPDAF`)
  - Version Tag `v3.9.0`, revision `bc80e3ac` from
    <https://github.com/parflow/parflow.git>
- **pdaf**
  - Version 2.0 in `/pdaf/src/PDAF_print_version.F90`, new version

Not used:
- **cosmo4\_21**
  - Gitlab REPO
    <https://icg4geo.icg.kfa-juelich.de/ModelSystems/tsmp_src/cosmo4.21_fresh.git>,
    `c81de76`
	- Makefile different to other **cosmo4\_21** versions

## Compile ParFlow standalone with PDAF (on JUWELS without compiler optimization) ##

      ./build_tsmp.ksh -m JUWELS -c pfl -v 1.1.0MCTPDAF -o -O0

The `-c` and `-v` flag are optional in this case because the inputs
for these flags on `-m JUWELS` are equivalent to the corresponding
default inputs.

      ./build_tsmp.ksh -m JUWELS -c pfl -v 1.1.0MCTPDAF -o -g

The `-g` flag, is build to produce debugging information.

## Compile CLM5 with PDAF ##

CLM5-PDAF compilation, currently available only on branch
`TSMP_pdaf-clm5` (January 2023).

``` bash
	./build_tsmp.ksh -v 4.4.0MCTPDAF -c clm -m JURECA -O Intel
```

### Prerequisite1: CLM5 preparation ###

For obtaining and preparing the component model `clm5_0`, you
have to run the following commands

``` bash
	git clone git@github.com:ESCOMP/CTSM.git clm5_0
	cd clm5_0
	git checkout release-clm5.0
	./manage_externals/checkout_externals
```

### Prerequisite2: Path to cesm ###

The path `CESMDATAROOT` needs to be changed in
`TSMP/bldsva/intf_oas3/clm5_0/arch/JURECA/config/softwarepaths.ksh`. (analogous
for `JUWELS`).

Currently you need to have access to compute project `cjicg41` and use the line:
```bash
	export CESMDATAROOT=/p/scratch/cjicg41/<user>/cesm
```
For `<user>`, please contact someone on the inside.
	
