# TSMP Build Examples #

**Attention**: Not all of this examples were tested and paths are
mostly made up.

## Builds ##

### Fully coupled TSMP with Oasis3-MCT ###

        ./build_tsmp.ksh -m JUWELS -c clm-cos-pfl -v 3.1.0MCT -O Intel
        ./build_tsmp.ksh -m JURECA -c clm-cos-pfl -v 3.1.0MCT -O Intel

`-O Intel`: Currently, on `JUWELS`, on stage `2020`, the compilation
with Intel-Compiler is recommended. It can explicitly be called by the
flag `-O Intel`.

The `-c` and `-v` flag are optional in this case because they are
default.

    -   **TSMP**

        -   revision `1b69ac4` on `master` Github repository `TSMP`

    -   **clm3\_5**

        -   Version 3.5, `clm3_5/Copyright`, `share3_070321` in
            `clm3_5/src/csm_share/ChangeLog`

        -   Gitlab repo
            <https://icg4geo.icg.kfa-juelich.de/ModelSystems/tsmp_src/clm3.5_fresh.git>,
            revision `801b5304`

    -   **cosmo5\_1**

        -   Gitlab repo
            <https://icg4geo.icg.kfa-juelich.de/ModelSystems/tsmp_src/cosmo5.01_fresh.git>,
            `f407b9b`

    -   **oasis3-mct**

        -   similar: svn revision r1506, `svn info` (caused error)

        -   now: Gitlab repo
            <https://icg4geo.icg.kfa-juelich.de/ModelSystems/tsmp_src/oasis3-mct.git>,
            tag `svn-r1506`, revision `bc58342`)

    -   **parflow3\_2**

        -   Version v3.2.0 in `VERSION`

        -   no `CMakeLists.txt` yet

        -   NOW: repo
            <https://icg4geo.icg.kfa-juelich.de/ModelSystems/tsmp_src/parflow3.2_legacy.git>,
            revision `98a8701` Remark: This is the same repository as the official ParFlow 
			repository on <https://github.com/parflow/parflow>

    -   **NO pdaf1\_1**

        -   Version 1.15.0 in `/pdaf1_1/src/PDAF-D_print_version.F90`,
            new version

### Compile ParFlow standalone on CLUMA2 without compiler optimization ###

      ./build_tsmp.ksh -m CLUMA2 -c pfl -v 1.1.0MCT -o "O0"

The `-m` and `-v` flag are optional in this case because they are
default.

### Compile new CLM4.0 + Cosmo5.1 on JURECA ###

      ./build_tsmp.ksh -m JURECA -c clm-cos -v 2.1.0MCT -C false

The new models are currently only supported for clm-cos and with
alternative coupling-scheme.

### Compile standalone eCLM ###

Branch `master`. Building eCLM. Currently, only standalone.

``` bash
	./build_tsmp.ksh -v 5.0.0 -c clm -m JURECA -O Intel
	./build_tsmp.ksh -v 5.0.0MCT -c clm -m JURECA -O Intel	
```

