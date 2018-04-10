#!/bin/csh
# This script is intended to update the DART namelist
# Uses, git_to_dart_XXXX for syncing with the DART source codes...
# Usage 
#
# USER SETUP ----DART NAMELIST PARAMETERS--------
set cutoff_radius = 0.0001 #radians
set ver_norm_hgt = 5000.  # meters
set horiz_dist_only = .false.  #logical
# USER SETUP ENDS---------------------------------
#
set ensemble_size = $1
set DART_DIR = $2

##
set coswork = ${DART_DIR}/cosmo/work
set clmwork = ${DART_DIR}/clm/work
set pflwork = ${DART_DIR}/parflow/work

# Update the dart namelist
./git_to_dart_cosmo.csh
./git_to_dart_clm.csh
./git_to_dart_parflow.csh
cd ${coswork}
echo $coswork
sed "s,__nens__,${ensemble_size}," -i input.nml
sed "s,__cutrad__,${cutoff_radius}," -i input.nml
sed "s,__vnhgt__,${ver_norm_hgt}," -i input.nml 
sed "s,__hdonly__,${horiz_dist_only}," -i input.nml 
cd ${clmwork}
sed "s,__nens__,${ensemble_size}," -i input.nml
sed "s,__cutrad__,${cutoff_radius}," -i input.nml
sed "s,__vnhgt__,${ver_norm_hgt}," -i input.nml
sed "s,__hdonly__,${horiz_dist_only}," -i input.nml
cd ${pflwork}
sed "s,__nens__,${ensemble_size}," -i input.nml
sed "s,__cutrad__,${cutoff_radius}," -i input.nml
sed "s,__vnhgt__,${ver_norm_hgt}," -i input.nml 
sed "s,__hdonly__,${horiz_dist_only}," -i input.nml 
exit 0 
