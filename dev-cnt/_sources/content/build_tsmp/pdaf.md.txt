# Structure of PDAF inside TSMP

```{toctree}
---
maxdepth: 3
---
build_preprocessor_variables.md
```

## Call Graph for TSMP-PDAF

Here is an exemplary call-graph from TSMP-PDAF's main program
`pdaf_terrsysmp()` inside `pdaf_terrsysmp.F90` up to the function,
where the enkf update is carried out:
`PDAF_enkf_analysis_rlm`/`PDAF_enkf_analysis_rsm`.

Inside TSMP, `TSMP/bldsva/intf_DA/pdaf/framework`:

``` text
pdaf_terrsysmp.F90: pdaf_terrsysmp
-> assimilate_pdaf.F90: assimilate_pdaf
-> PDAF_assimilate_enkf
```

Inside PDAF, `pdaf/src/`:

``` text
PDAF-D_assimilate_enkf.F90: PDAF_assimilate_enkf
-> PDAF-D_put_state_enkf.F90: PDAF_put_state_enkf
-> PDAF-D_enkf_update.F90: PDAF_enkf_update

-> PDAF-D_enkf_analysis_rlm.F90: PDAF_enkf_analysis_rlm
or
-> PDAF-D_enkf_analysis_rsm.F90: PDAF_enkf_analysis_rsm
```

## Important Functions of TSMP-PDAF

### `pdaf_terrsysmp` (main program)

Subroutines:

- `init_parallel_pdaf()`
  1. `MPI_Init` if not yet done (should be done)
  2. Set `npes_world`, `mype_world`
  3. Read command line input: `n_modeltasks`
  4. Initialize communicators for ensemble evaluations
	- consistency checks `n_modeltasks > npes_world` and (not used
since `dim_ens` is set to dummy-zero) `n_modeltasks > dim_ens`
  5. `COMM_ensemble` (local copy of `MPI_COMM_WORLD`)
  6. `COMM_model`: Generate communicators for model runs
	 - Set `local_npes_model` depending on complete number
       `npes_world` and `n_modeltasks`
	 - Example: `mype_world=24` and 4 PEs per model, then mype_world
       will be model task/color `9 = floor(24/4)+1` and rank `0`
       (`25->1`, `26->2`, `27->3`).
  7. `COMM_filter`: Same as `COMM_model`, but only first task used.
  8. `COMM_couple`: For each rank in `COMM_model` (example 4 ranks),
     there is a task in `COMM_couple` with `n_modeltasks` ranks
     (example 48)
  9. Model equivalents of `COMM_model` are set: `COMM_model_oas`,
     `COMM_model_pfl`, `COMM_model_clm`
  10. Model equivalent of `COMM_couple` is set: `COMM_model_clm`

- [`initialize_tsmp()`](#initialize_tsmp)
  1. read parameter file `enkfpf.par`
  2. Set model number for a process
  3. Set input function names
  4. Call model-specific input functions
	- `clm_init`
	- `enkfparflowinit` (some more allocations in ParFlow part)
    - `cosmo_init`

### `initialize_tsmp()`

File: `./bldsva/intf_DA/pdaf/model/wrapper_tsmp.c`

Subroutines:
- `clm_init`
  1. Setting `nlfilename` (essentially the same as `finname`?)
  2. Initialize Communicators depending on mode
  3. Initialize ESMF (needed for time-manager)
  4. Initialize timing library, and set full path to namelist
  5. Initialize Orbital parameters (`shr_orb_params`)
  6. Initialize land model (`clm_init0, clm_init1, clm_init2`)
  7. Initialize "external" atmospheric forcing (`atmdrv_init`)
- `enkfparflowinit`
  1. 

## List of interesting variables of TSMP-PDAF

- `local_npes_model`: PEs per ensemble
- `n_modeltasks`: Number of parallel model tasks
- `npes_world`: Number of MPI processes in `MPI_COMM_WORLD`
- `mype_world`: Rank of calling process in `MPI_COMM_WORLD`
- `pf_statevec`: ParFlow-Statevector in C-wrapper
  - from `enkf_parflow.h` and set in `initialize_tsmp` from the
    `wrapper_tsmp.c`
- `pf_statevec_fortran`: ParFlow-Statevector in Fortran-part of
  TSMP-PDAF
  - defined in `init_pdaf.F90` by `C_F_POINTER` subroutine (Doc from
    Gnu-Compiler:
    <https://gcc.gnu.org/onlinedocs/gfortran/C_005fF_005fPOINTER.html>)
  - takes values from `pf_statevec`
- `subvec_p`: state vector pressure (from `enkf_parflow.c:371`)
  - defined: `enkf_parflow.h`, allocated: `wrapper_tsmp.c`
- `tag_model_clm`, (`0`), `tag_model_parflow` (`1`) and
  `tag_model_cosmo` (`2`) identifiers for the component models
    - defined `mod_tsmp.F90`
    - used throughout TSMP-PDAF, for example in
      `collect_state_pdaf.F90`

### Model Communicator `COMM_model`

Variable names for the model communicator:

-   ParFlow standalone:
    - `COMM_model` (`init_parallel_pdaf.F90`)
	- `COMM_model_pfl` (`init_parallel_pdaf.F90`)
	- C: `COMM_model_pfl` (`enkf_parflow.h`)
	- C: `pfcomm` (`enkfparflowinit` in `enkf_parflow.c`)
	- C: handed to `amps_EmbeddedInitComm`

-   CLM (3.5) standalone:
    - `COMM_model` (`init_parallel_pdaf.F90`)
	- `COMM_model_clm` (`init_parallel_pdaf.F90`)
	- `COMM_model_clm` (`clm_init` in `enkf_clm.F90`)
	- handed to `spmd_init` and `mct_world_init`


-   Coupled simulation using OASIS:
    - `COMM_model` (`init_parallel_pdaf.F90`)
	- `COMM_model_oas` (`init_parallel_pdaf.F90`)
	- handed over to OASIS to `oasis_init_comp` in
      `mod_oasis_method.F90`

The two call-sequences of loading `COMM_model_oas` in coupled
simulations:
- ParFlow: Sets `amps_CommWorld`
  - `enkfparflowinit` in `enkf_parflow.c`
  - `amps_init` in `oas3/amps_init.c` (entering ParFlow)
  - `oas_pfl_init` in `oas3/oas_pfl_init`
  - `prism_init_comp_proto` (`mod_prism`) / `oasis_init_comp` in
    `mod_oasis_method.F90` (entering OASIS)
- CLM (3.5): Sets `kl_comm`
  - `clm_init` in `enkf_clm.F90`
  - `oas_clm_init` in `oas_clm_init.F90` (entering CLM)
  - `prism_init_comp_proto` (`mod_prism`)/ `oasis_init_comp` in
    `mod_oasis_method.F90` (entering OASIS)


### AMPS used by ParFlow in TSMP-PDAF

Information on AMPS (Another Message Passing System):
<https://github.com/parflow/parflow/blob/master/pfsimulator/amps/oas3/intro.dxx>

The directory for AMPS used in TSMP-PDAF is set using the
CMAKE-variable `PARFLOW_AMPS_LAYER` in `build_interface_parflow.ksh`.

AMPS directory / `PARFLOW_AMPS_LAYER` for
-   ParFlow standalone: `da` (`mpi1`)
-   Coupled ParFlow: `oas3`

#### AMPS-Layer `da` / `mpi1` ####

`da` is a new AMPS-Layer introduced for ParFlow with TSMP-PDAF.

`da` is based on the existing AMPS-Layer `mpi1`.

The changes between `da` and `mpi1` are small, therefore in
TSMP2-PDAF, `da` will be removed and replaced by a pre-patched version
of `mpi1`.

Changes:
- `dacomm` communicator defined and suspected to be
  unused. ParFlow-standalone testmodel needed to confirm.
- commented out `MPI_Finalize`

#### AMPS-Layer `oas3` ####

`oas3` is an existing AMPS-Layer.

`oas3` uses the communicator setup defined in
`prism_init_comp_proto`/`oasis_init_comp`, where `COMM_model_oas` is
inserted.

### CMEM ###

If you want to use the CMEM operator, you have to download the source
code under https://confluence.ecmwf.int/display/LDAS/CMEM+Download
(Tested: Version 5.1) and copy it into the `framework` directory of
TSMP-PDAF.
