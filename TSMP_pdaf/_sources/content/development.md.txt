# Developing TSMP #

## Supporting a new machine ##

First of all make sure that the new machine is introduced to the
`supported_machines.ksh` script.

Create the necessary folder structure under `intf_oas3` and create a
script for this machine and every component-model that will be
supported and implement the interface.

If necessary add `src/config` folders if machine specific source-code or
`config/make` files are needed.

Add a folder with the machine name under `/machines`. Also add the
machine-interface and implement the interface. Also add files to load
modules if necessary.

For both cases it is easiest to copy a similar machine and change names
and specifics accordingly.

## Setting up a new experiment ##

First of all make sure that the new setup is introduced to the
`supported_machines.ksh` script.

Add a folder with the setup name under `/setups`.

Also add the setup-interface for all machines that are supported and
implement the interface.

Also add namelists for all component-models.

The namelist name must match to the according entry in the setup
interface.

namelists for `cosmo5.1` and `clm4.0` must have the version as suffix
even if declared otherwise in the setup interface. For example
`lmrun_uc5_1` and `lnd.stdin4_0`. But in the setup-interface it must
be written `lmrun_uc` and `lnd.stdin`.

Also for ensembles in TSMP-PDAF the ensemble namelists must then be named like
`lmrun_uc5_1_4711` or `lnd.stdin4_0_4711` for ensemble member 4711 even
if `lmrun_uc` and `lnd.stdin` is specified by the setup interface or
flags. Again it is easiest to copy a similar setup and change names and
specifics accordingly.

