# Building TSMP #

These are example builds of TSMP that are frequently tested.
	
**Build script of TSMP**

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

for building `TSMP` with COSMO-CLM3.5-ParFlow the following internal version is recommended:
-   `3.1.0MCT` (recommended).

for building `TSMP-PDAF` the following internal versions can be used:
-   `1.1.0MCTPDAF`
-   `3.0.0MCTPDAF`
-   `3.1.0MCTPDAF` (recommended).

```{toctree} 
---
maxdepth: 3
caption: Build Examples for TSMP-PDAF
---
build_examples_tsmp.md
build_examples_tsmppdaf.md
build_environment_variables.md
```
