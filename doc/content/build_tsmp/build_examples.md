# Building TSMP #

These are example builds of TSMP that are frequently tested.
	
## Build script of TSMP

To get information about TSMP build options for different version
using standalone or different combination of model components on
different machines inside the folder `TSMP/bldsva` execute

``` bash
	./build_tsmp.ksh -a
```

you can also use `--man` or `--help` to get terrsysmp build options

``` bash
	./build_tsmp.ksh --man ./build_tsmp.ksh --help
```

for building `TSMP` with COSMO-CLM3.5-ParFlow the following combination is recommended:
-   `-c clm3-cos5-pfl ` (recommended).

for building `TSMP-PDAF` the following combinations can be used:
-   `-c clm3-cos4-pfl-pdaf`
-   `-c clm3-cos5-pfl-pdaf` (recommended).

### Building on GPU

For building ParFlow 3.12 on GPU:

```shell
cd $TSMP_DIR/bldsva
./build_tsmp.ksh -c clm3-cos5-pfl -m JUWELS -O Intel -A GPU
```

## Build Examples

```{toctree} 
---
maxdepth: 3
---
build_examples_tsmp.md
build_examples_tsmppdaf.md
```

## Patching the original source code 

In order to prepare the original Cosmo5_1 source code for coupling, the Cosmo source code is patched using the diff files which are located in `$TSMP_DIR/bldsva/intf_oas3/cosmo5_1/pfile`.
The diff files are generated using the script `fpatch.sh` which is located in the same directory. Note that the patching is done only for cosmo5_1 source code.
If the user aims to modify the coupling source files, should modify the files located in `$TSMP_DIR/bldsva/intf_oas3/cosmo5_1/tsmp` and run again the script `fpatch.sh` in order to generate the the new diff files accordingly.
The necessary changes for making the ParFlow3.12 ready for coupling are already included in the official ParFlow3.12 release.
The changed source files for CLM3.5 are located in  `$TSMP_DIR/bldsva/intf_oas3/clm3_5/tsmp` and will be copied by TSMP scripts to the user's source code directory of CLM3.5 '(`bld/usr.src`)

## Automatic Porting of TSMP on x86 machines

The users who want to port TSMP on GENERIC_X86 Linux, the TSMP team provided a script to install all the necessary libraries (Netcdf, GRIBAPI, OpenMPI, HDF5, TCL, Hypre and Silo) automatically in TSMP root directory. Please run the script "lib_install.sh" located in bldsva directory to install the libraries. Note that if you exported already one of the libraries Netcdf, HDF5, GRIBAPI, Silo, Hypre and TCL in the .baschrc or .profile, you need to comment them out in order to not mess up the installation via the script lib_install.sh.

The instruction on how to build and how to configure TSMP for the cordex test case ([EURO-CORDEX test case experiment](./../gettingstarted.md/#the-fully-coupled-pan-european-euro-cordex-evaluation-experiment-with-tsmp)) are given in [step 5 ](./../gettingstarted.md/#step-5-build-tsmp-interface-and-component-models) and [step 6 ](./../gettingstarted.md/#step-6-setup-and-configuration-of-the-respective-usage-and-test-case) in the getting started section, respectively.
