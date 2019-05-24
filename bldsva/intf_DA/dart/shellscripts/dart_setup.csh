#!/bin/csh
# This script is intended to update the DART namelist
# Uses, git_to_dart_XXXX for syncing with the DART source codes...
# Usage ./dart_setup.csh 20 /home/pshrestha/DART/lanai/models/terrsysmp 
#
# USER SETUP ----DART NAMELIST PARAMETERS--------
# soil-vegetation
set sv_cutoff_radius = 0.0001 #radians
set sv_ver_norm_hgt = 5000.  # meters
set sv_horiz_dist_only = .false.  #logical

# atmosphere
set a_cutoff_radius = 0.0005 #radians
set a_ver_norm_hgt = 200000.  # meters
set a_horiz_dist_only = .false.  #logical
# USER SETUP ENDS---------------------------------
#
set ensemble_size = $1
set DART_DIR = $2

##
set coswork = ${DART_DIR}/cosmo/work
set clmwork = ${DART_DIR}/clm/work
set pflwork = ${DART_DIR}/parflow/work

# Update the dart namelist
cd ${coswork}
echo $coswork
sed "s,__nens__,${ensemble_size}," -i input.nml
sed "s,__cutrad__,${a_cutoff_radius}," -i input.nml
sed "s,__vnhgt__,${a_ver_norm_hgt}," -i input.nml 
sed "s,__hdonly__,${a_horiz_dist_only}," -i input.nml 
cd ${clmwork}
sed "s,__nens__,${ensemble_size}," -i input.nml
sed "s,__cutrad__,${sv_cutoff_radius}," -i input.nml
sed "s,__vnhgt__,${sv_ver_norm_hgt}," -i input.nml
sed "s,__hdonly__,${sv_horiz_dist_only}," -i input.nml
cd ${pflwork}
sed "s,__nens__,${ensemble_size}," -i input.nml
sed "s,__cutrad__,${sv_cutoff_radius}," -i input.nml
sed "s,__vnhgt__,${sv_ver_norm_hgt}," -i input.nml 
sed "s,__hdonly__,${sv_horiz_dist_only}," -i input.nml 
exit 0 
