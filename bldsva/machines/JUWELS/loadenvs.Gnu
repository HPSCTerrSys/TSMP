####!/bin/tcsh
#2020-11-26_Agh
#USAGE="source $0"
#for bld and work w/ TerrSysMP on JUWELS
#for standard user

module --force purge
module use $OTHERSTAGES
module load Stages/2019a
module load GCC/8.3.0
module load ParaStationMPI/5.2.2-1

module load Hypre/2.15.1
module load Silo/4.10.2
# module load Tcl/8.6.9 # this version is not compatible with ParFlow3.2

module load HDF5/1.10.5
module load netCDF/4.6.3
module load netCDF-Fortran/4.4.5
module load parallel-netcdf
module load ecCodes/2.12.0
module load cURL/7.64.1.lua
module load Szip/.2.1.1
module load Python/2.7.16
module load CMake
module load NCO
#module load CDO/1.9.5
export EBROOTTCL="/p/project/cslts/local/juwels/tcl8.6.8_gcc803"
export PATH="/p/project/cslts/local/juwels/tcl8.6.8_gcc803/bin:$PATH"
export LD_LIBRARY_PATH="/p/project/cslts/local/juwels/tcl8.6.8_gcc803/lib:$LD_LIBRARY_PATH"
module li
