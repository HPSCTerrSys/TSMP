# Debugging flags

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

