# Debugging TSMP #

See also [Debug Tips](build_examples_tsmppdaf.md#debug-tips) for building
TSMP!

If an error appears and the source can not be localized by the default
outputs, it might make sense to build TSMP-PDAF again, but with more
debug-flags enabled.

For `JUWELS` these flags can be set in

- `bldsva/intf_DA/pdaf1_1/arch/JUWELS/build_interface_pdaf1_1_JUWELS.ksh`
- `bldsva/machines/JUWELS/build_interface_JUWELS.ksh`


In `build_interface_JUWELS.ksh`, there is `defaultOptC`, which can be
set, f.e. as

``` ksh
defaultOptC="-g -O0 -xHost -traceback" # Intel
```

For TSMP-PDAF: In `build_interface_pdaf1_1_JUWELS.ksh`, there is the general
`importFlags` that could for example be set as

```ksh
importFlags="-g -traceback"
```

Watch out here, to get the right condition (Intel / Fortran
compilers)! The examples are for the Intel compiler, gcc-flags are
slightly different.

## Debugging with TotalView

It is possible to debug TSMP with multiple coupled components with the mean of the TotalView debugger. 

Please contact the SDLTS for support for the initial set up TotalView on JSC machines. 

1) Start interactive session. E.g. on JUWELS machine 
```sh
salloc --partition=devel --nodes=5 --account=PROJECT --time=01:30:00
```

2) Setup and modifiy the submission script of TSMP
```sh
totalview -args srun --multi-prog slm_multiprog_mapping.conf
```
