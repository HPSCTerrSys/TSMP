# Input files for ParFlow #

```{danger}
**DEPRECATED**: This documentation is deprecated and no longer maintained. 
Please refer to the updated documentation for current information.

   - TSMP2: <https://hpscterrsys.github.io/TSMP2>
   - TSMP-PDAF: <https://hpscterrsys.github.io/pdaf>
```

Similar to the input for CLM, the used has to create individual
ParFlow data base files for each model realisation which share the
same file name prefix. For example, if the chosen ParFlow file name
prefix is `pfinput`, the file names for the first three model
realisation are:

``` text
pfinput_00000.pfidb
pfinput_00001.pfidb
pfinput_00002.pfidb
```

As for the CLM input, the timing information contained in the ParFlow
input files needs to be consistent among the different realisation.
Other input information like initial conditions, parameters, etc. can
differ between the realisations.

