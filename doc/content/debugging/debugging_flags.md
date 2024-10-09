# Debugging flags

If an error appears and the source can not be located by the default
outputs, it might make sense to build TSMP-PDAF with more debug-flags
enabled.

Quick debug build: Add the option `-D true` to your original
`./build_tsmp.ksh` command. E.g. 

```
./build_tsmp.ksh -c clm5-pdaf -m JURECA -O Intel
```

becomes

```
./build_tsmp.ksh -c clm5-pdaf -D true -m JURECA -O Intel
```

## More detailed information

Debugging flags are set in

- `bldsva/intf_DA/pdaf/arch/build_interface_pdaf.ksh`
- `bldsva/machines/config_JUWELS.ksh`


The following examples are for the Intel compiler, gcc-flags are
slightly different.

In `config_JUWELS.ksh`, there is `defaultOptC`, which can be
set, f.e. as

``` ksh
defaultOptC="-g -O0 -xHost -traceback" # Intel
```

For TSMP-PDAF: In `build_interface_pdaf.ksh`, there is the general
`importFlags` that could for example be set as

```ksh
importFlags="-g -traceback"
```

