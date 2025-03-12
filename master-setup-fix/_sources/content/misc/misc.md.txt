# Misc #

Various topics related to TSMP-PDAF.

## Gitlab: TSMP-PDAF group ##

The main repository of TSMP-PDAF is found on GitHub
([https://github.com/HPSCTerrSys/TSMP](https://github.com/HPSCTerrSys/TSMP)).

Additional development of TSMP-PDAF is organized on Gitlab in the
group TSMP-PDAF:

[https://icg4geo.icg.kfa-juelich.de/ExternalRepos/tsmp-pdaf](https://icg4geo.icg.kfa-juelich.de/ExternalRepos/tsmp-pdaf)

### Repositories ###

The Gitlab-group TSMP-PDAF currently contains the following
repositories

- `TSMP`: main repository, clone of GitHub TSMP from
  [https://github.com/HPSCTerrSys/TSMP](https://github.com/HPSCTerrSys/TSMP)
- `TSMP-PDAF Manual`: Repository for this manual.
- `TSMP-PDAF Build Scripts`: Scripts for (remote) building of TSMP
- `TSMP-PDAF Component Models`: Component Models used for building TSMP-PDAF
- `TSMP-PDAF Testcases`: Testcases for TSMP-PDAF
- `TSMP-PDAF Presentations`: Presentation slides for TSMP-PDAF

## Open loop compilation of (TSMP-)PDAF ##

For compiling TSMP-PDAF with PDAF but without data assimilation, you
can set the flag `PDAF_NO-UPDATE` in Makefiles under
`./TSMP/bldsva/intf_DA/pdaf/arch/<machine>/config/<compiler>.h.`
and in the line, where `CPP_DEFS` are set. Usually only `-DUSE_PDAF`
is set, here you can add `-DPDAF_NO_UPDATE`.

If you compile once with this flag and once without, you should obtain
two build-versions of TSMP-PDAF, one for open-loop runs and one for
data assimilation runs.

## Scripts and other resources ##

- Script for calculating discharge for ParFlow along the river network
by Carina Furusho. Discharge can be calculated for each grid cell and
each time step:
<https://icg4geo.icg.kfa-juelich.de/ModelSystems/realpep_p4/blob/master/tools_postprocessing/functions_output_UNIVERSAL.py>

