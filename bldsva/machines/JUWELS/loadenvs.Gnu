#!/bin/tcsh
#2019-01-28_SBr
#USAGE="source $0"
#for bld and work w/ TerrSysMP on JUWELS
#for standard user

module --force purge

# Only necessary if modules are missing after maintenence
module use /gpfs/software/juwels/otherstages/

module load Stages/2018a
#
# standard modules on Juwels
module load StdEnv
module load jscslurm/.17.02.11
module load jsctools/.0.1
module load .juwels-env

# Only necessary if modules are missing after maintenence
module load GCC/7.3.0
module load ParaStationMPI/5.2.1-1

# libs for TerrSysMP
module load Tcl/.8.6.8
module load Hypre/2.14.0
module load Silo/4.10.2

#grib1 dwd
module load netCDF/4.6.1
module load netCDF-Fortran/4.4.4
module load parallel-netcdf/1.9.0


# standard analysis and vis. libs and tools 
module load NCO/4.7.4
#module load NCL/6.5.0
#module load CDO/1.9.3
module load ncview/2.1.7
module load Python/2.7.14
#module load matplotlib/1.3.1-Python-2.7.9
#module load netcdf4-python/1.4.2-Python-3.6.6
module load SciPy-Stack/2018a-Python-2.7.14
#module load basemap/1.0.7-Python-2.7.14
module load h5py/2.7.1-Python-2.7.14
#module load petsc4py/3.9.0-Python-3.6.5
module load mpi4py/3.0.0-Python-2.7.14
#image magick / feh ... some standard viewer is missing

module li
