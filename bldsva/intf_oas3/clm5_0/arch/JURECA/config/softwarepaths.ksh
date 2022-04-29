#! /bin/ksh

# Specific software paths:
# Determined based on loaded modules.
# which 'specific config tool for that netcdf path' 
# gives the full path including /bin/config 
# sed 's#\(bin\).*#\1#' selects everything before and including "bin"
# sed 's#\(.*\)...#\1#' selects everything except the last 3 characters i.e. bin
 export NETCDF=$(which nc-config | sed 's#\(bin\).*#\1#' | sed 's#\(.*\)...#\1#')
 export NETCDF_FORTRAN=$(which nf-config | sed 's#\(bin\).*#\1#' | sed 's#\(.*\)...#\1#')
 export NETCDF_C=$(which ncxx4-config | sed 's#\(bin\).*#\1#' | sed 's#\(.*\)...#\1#')
 export PNETCDF=$(which pnetcdf-config | sed 's#\(bin\).*#\1#' | sed 's#\(.*\)...#\1#')

# User defined data paths
# Path to the root directory for all CESM / CLM input files
# TODO (change to shared with other users)
 export CESMDATAROOT=$SCRATCH/cesm
# # Path to the CSM specific input files (usually subfolder of CESMDATAROOT)
 export CSMDATA=$CESMDATAROOT/inputdata
#
