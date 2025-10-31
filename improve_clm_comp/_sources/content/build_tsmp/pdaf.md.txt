# Structure of PDAF inside TSMP

```{toctree}
---
maxdepth: 3
---
build_environment_variables.md
```

## Call Graph for TSMP-PDAF

Here is an exemplary call-graph from TSMP-PDAF's main program
`pdaf_terrsysmp()` inside `pdaf_terrsysmp.F90` up to the function,
where the enkf update is carried out:
`PDAF_enkf_analysis_rlm`/`PDAF_enkf_analysis_rsm`.

Inside TSMP, `TSMP/bldsva/intf_DA/pdaf1_1/framework`:

``` text
pdaf_terrsysmp.F90: pdaf_terrsysmp
-> assimilate_pdaf.F90: assimilate_pdaf
-> PDAF_assimilate_enkf
```

Inside PDAF, `pdaf1_1/src/`:

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
  9. Finally `da_comm` is set to `COMM_model`; and
     `cosmo_input_suffix` is set

- [`initialize_tsmp()`](#initialize_tsmp)
  1. read parameter file `enkfpf.par`
  2. Set model number for a process
  3. Set input function names
  4. Call model-specific input functions
	- `clm_init`
	- `enkfparflowinit` (some more allocations in ParFlow part)
    - `cosmo_init`

### `initialize_tsmp()`

File: `./bldsva/intf_DA/pdaf1_1/model/wrapper_tsmp.c`

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
    https://gcc.gnu.org/onlinedocs/gfortran/C_005fF_005fPOINTER.html)
  - takes values from `pf_statevec`
- `subvec_p`: state vector pressure (from `enkf_parflow.c:371`)
  - defined: `enkf_parflow.h`, allocated: `wrapper_tsmp.c`
- `tag_model_clm`, (`0`), `tag_model_parflow` (`1`) and
  `tag_model_cosmo` (`2`) identifiers for the component models
    - defined `mod_tsmp.F90`
    - used throughout TSMP-PDAF, for example in
      `collect_state_pdaf.F90`

### Model Communicators

AMPS directory:
-   ParFlow standalone: amps directory `da/` is used,
-   Coupled ParFlow-CLM component models: amps directory `oas3/` is used.

Variable names for the model communicator:
-   ParFlow standalone:
    -   `comm_model` -> C: `comm_model_pdaf` -> `pfcomm` (enkf_parflow/enkfparflowinit) -> `amps_CommWorld` (da/amps_init.c)
-   Coupled with OASIS:
    -   `comm_model` -> `da_comm` (init_parallel_pdaf) -> C: `fsubcomm` -> C-version not used
    -   `comm_model` -> `da_comm` (init_parallel_pdaf) -> `mpi_comm_global` (mod_oasis_method)
