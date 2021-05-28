#2020-01-28_Agh
#USAGE="source $0"
#for standard user

module --force purge
module use $OTHERSTAGES
module load Stages/2020
module load GCC/9.3.0
module load ParaStationMPI/5.4.7-1

module load ecCodes/2.18.0
module load Hypre/2.20.0
module load Silo/4.10.2
# module load Tcl/8.6.10

module load HDF5/1.10.6
module load netCDF/4.7.4 netCDF-Fortran/4.5.3
module load parallel-netcdf/1.12.1
module load cURL/7.71.1
module load Szip/.2.1.1
ml Python/3.8.5
module load NCO/4.9.5
module load CMake/3.18.0
#
export EBROOTTCL="/p/project/cslts/local/juwels/tcl8.6.8_gcc903"
export PATH="/p/project/cslts/local/juwels/tcl8.6.8_gcc903/bin:$PATH"
export LD_LIBRARY_PATH="/p/project/cslts/local/juwels/tcl8.6.8_gcc903/lib:$LD_LIBRARY_PATH"
#
#module load NCO
#module load CDO
module li
