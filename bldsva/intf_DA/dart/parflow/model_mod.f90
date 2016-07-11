! This code is not necessarily under the DART copyright ...
!

! DART $Id: $

module model_mod

! This module contains routines to work with ParFlow Binary data 
! in the DART framework

! Based on DART template directory

! Modules that are absolutely required for use are listed
use        types_mod, only : r8, MISSING_R8

use time_manager_mod, only : time_type, set_time, set_date, get_date, get_time, &
                             print_time, print_date, set_calendar_type

use     location_mod, only : location_type,      get_close_maxdist_init,        &
                             get_close_obs_init, get_close_obs, set_location,   &
                             set_location,                                      &
                             vert_is_height,       VERTISHEIGHT

use    utilities_mod, only : register_module, error_handler, nc_check, &
                             get_unit, open_file, close_file, E_ERR, E_MSG, &
                             nmlfileunit, do_output, do_nml_file, do_nml_term,  &
                             find_namelist_in_file, check_namelist_read
use netcdf

implicit none
private

! required by DART code - will be called from filter and other
! DART executables.  interfaces to these routines are fixed and
! cannot be changed in any way.
public :: get_model_size,         &
          adv_1step,              &
          get_state_meta_data,    &
          model_interpolate,      &
          get_model_time_step,    &
          end_model,              &
          static_init_model,      &
          init_time,              &
          init_conditions,        &
          nc_write_model_atts,    &
          nc_write_model_vars,    &
          pert_model_state,       &
          get_close_maxdist_init, &
          get_close_obs_init,     &
          get_close_obs,          &
          ens_mean_for_model

! not required by DART but for larger models can be useful for
! utility programs that are tightly tied to the other parts of
! the model_mod code.
public :: get_state_vector, &
          get_parflow_filename, &
          dart_vector_to_model_file

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL: https://proxy.subversion.ucar.edu/DAReS/DART/releases/Lanai/models/template/model_mod.f90 $"
character(len=32 ), parameter :: revision = "$Revision: 6256 $"
character(len=128), parameter :: revdate  = "$Date: 2013-06-12 18:19:10 +0200 (Wed, 12 Jun 2013) $"

character(len=256) :: string1, string2, string3
logical, save :: module_initialized = .false.

! Dimension and grid resolution from PFB file
integer(kind=4)                ::  &
                   nx,             &     ! Longitude dimension
                   ny,             &     ! Latitude dimension
                   nz                    ! Vertical dimension
real(r8)                       ::  &
                   dx,             &     ! Grid Resolution in m
                   dy,             &     ! Grid Resoultion in m
                   dz                    ! Grid Resolution scale, check parflow namelist

! EXAMPLE: define model parameters here
integer                          :: model_size
type(time_type)                  :: time_step
type(location_type), allocatable :: state_loc(:)

! Everything needed to describe a variable

type progvartype
   private
   character(len=NF90_MAX_NAME) :: varname
   character(len=NF90_MAX_NAME) :: long_name
   character(len=NF90_MAX_NAME) :: units
   integer  :: numdims
   integer  :: numEW
   integer  :: numNS
   integer  :: numZ
   integer  :: varsize     ! prod(dimlens(1:numdims))
   integer  :: index1      ! location in dart state vector of first occurrence
   integer  :: indexN      ! location in dart state vector of last  occurrence
   integer  :: rangeRestricted
   real(r8) :: minvalue
   real(r8) :: maxvalue
   character(len=256) :: kind_string
   logical  :: update
end type progvartype

type(progvartype)     :: progvar

real(r8), allocatable ::   lon(:,:)      ! longitude degrees_east
real(r8), allocatable ::   lat(:,:)      ! latitude  degrees_north
real(r8), allocatable ::   vcoord(:)     ! vertical co-ordinate in m
type(time_type)       ::   parflow_time  ! parflow ouptut time

! EXAMPLE: perhaps a namelist here 
character(len=256) :: parflow_file                 = 'parflow_file'
character(len=256) :: pfidb_file                   = 'pfidb_file'
character(len=256) :: grid_file                    = 'grid_file'
logical            :: output_state_vector_as_1D    = .false.
integer            :: assimilation_period_days     = 0
integer            :: assimilation_period_seconds  = 21400
integer            :: debug = 0 

namelist /model_nml/                &
      parflow_file,                 &
      pfidb_file,                   &
      grid_file,                    &
      output_state_vector_as_1D,    &
      assimilation_period_days,     &
      assimilation_period_seconds,  &
      debug

contains


!==================================================================


subroutine static_init_model()
!------------------------------------------------------------------
!
! Called to do one time initialization of the model. As examples,
! might define information about the model size or model timestep.
! In models that require pre-computed static data, for instance
! spherical harmonic weights, these would also be computed here.
! Can be a NULL INTERFACE for the simplest models.

real(r8) :: x_loc
integer  :: i
integer  :: iunit, io

if ( module_initialized ) return ! only need to do this once.

! Since this routine calls other routines that could call this routine
! we'll say we've been initialized pretty dang early.
module_initialized = .true.

! Print module information to log file and stdout.
call register_module(source, revision, revdate)

! This is where you would read a namelist, for example.
call find_namelist_in_file("input.nml", "model_nml", iunit)
read(iunit, nml = model_nml, iostat = io)
call check_namelist_read(iunit, io, "model_nml")

! Record the namelist values used for the run ...
if (do_nml_file()) write(nmlfileunit, nml=model_nml)
if (do_nml_term()) write(     *     , nml=model_nml)

call set_calendar_type('Gregorian')

! Get the dimensions of the grid and the grid variables from the netCDF file.
call pfread_dim(parflow_file)

if (debug > 0 .and. do_output()) write(*,*) '   ... pfb dimensions are ' , nx, ny, nz
if (debug > 0 .and. do_output()) write(*,'(A,3(1X,F9.2))') '    ... pfb grid resoluion are ', dx, dy, dz

model_size = nx*ny*nz

! Get the vertical co-ordinate
allocate(vcoord(nz))
call pfidb_read(pfidb_file)

! Get the geo-locaion
allocate(lon(nx,ny))
allocate(lat(nx,ny))

call grid_read(grid_file)
if (debug > 0 .and. do_output()) then
   write(*, '(A,2(1X,F9.2))') '    ...parflow lon', minval(lon), maxval(lon)
   write(*, '(A,2(1X,F9.2))') '    ...parflow lat', minval(lat), maxval(lat)
end if

!CPS call error_handler(E_ERR,'static_init_model','routine not tested',source, revision,revdate)

! Create storage for locations
allocate(state_loc(model_size))

! Define the locations of the model state variables
! naturally, this can be done VERY differently for more complicated models.
! set_location() is different for 1D vs. 3D models, not surprisingly.
do i = 1, model_size
   x_loc = (i - 1.0_r8) / model_size
   ! must do one of these:
   !state_loc(i) =  set_location(x_loc)
   !state_loc(i) =  set_location(x_loc,y_loc,v_loc,v_type)
end do

! The time_step in terms of a time type must also be initialized.
time_step = set_time(assimilation_period_seconds, assimilation_period_days)

end subroutine static_init_model




subroutine init_conditions(x)
!------------------------------------------------------------------
! subroutine init_conditions(x)
!
! Returns a model state vector, x, that is some sort of appropriate
! initial condition for starting up a long integration of the model.
! At present, this is only used if the namelist parameter 
! start_from_restart is set to .false. in the program perfect_model_obs.
! If this option is not to be used in perfect_model_obs, or if no 
! synthetic data experiments using perfect_model_obs are planned, 
! this can be a NULL INTERFACE.

real(r8), intent(out) :: x(:)

if ( .not. module_initialized ) call static_init_model

call error_handler(E_ERR,'init_conditions','routine not tested',source, revision,revdate)
x = MISSING_R8

end subroutine init_conditions



subroutine adv_1step(x, time)
!------------------------------------------------------------------
! subroutine adv_1step(x, time)
!
! Does a single timestep advance of the model. The input value of
! the vector x is the starting condition and x is updated to reflect
! the changed state after a timestep. The time argument is intent
! in and is used for models that need to know the date/time to 
! compute a timestep, for instance for radiation computations.
! This interface is only called if the namelist parameter
! async is set to 0 in perfect_model_obs of filter or if the 
! program integrate_model is to be used to advance the model
! state as a separate executable. If one of these options
! is not going to be used (the model will only be advanced as
! a separate model-specific executable), this can be a 
! NULL INTERFACE.

real(r8),        intent(inout) :: x(:)
type(time_type), intent(in)    :: time

if ( .not. module_initialized ) call static_init_model

write(string1,*) 'Cannot advance ParFlow with a subroutine call; async cannot equal 0'
call error_handler(E_ERR,'adv_1step',string1,source,revision,revdate)

end subroutine adv_1step



function get_model_size()
!------------------------------------------------------------------
!
! Returns the size of the model as an integer. Required for all
! applications.

integer :: get_model_size

if ( .not. module_initialized ) call static_init_model

get_model_size = model_size

end function get_model_size



subroutine init_time(time)
!------------------------------------------------------------------
!
! Companion interface to init_conditions. Returns a time that is somehow 
! appropriate for starting up a long integration of the model.
! At present, this is only used if the namelist parameter 
! start_from_restart is set to .false. in the program perfect_model_obs.
! If this option is not to be used in perfect_model_obs, or if no 
! synthetic data experiments using perfect_model_obs are planned, 
! this can be a NULL INTERFACE.

type(time_type), intent(out) :: time

if ( .not. module_initialized ) call static_init_model

call error_handler(E_ERR,'init_time','routine not tested',source, revision,revdate)
! for now, just set to 0
time = set_time(0,0)

end subroutine init_time



subroutine model_interpolate(x, location, itype, obs_val, istatus)
!------------------------------------------------------------------
!
! Given a state vector, a location, and a model state variable type,
! interpolates the state variable field to that location and returns
! the value in obs_val. The istatus variable should be returned as
! 0 unless there is some problem in computing the interpolation in
! which case an alternate value should be returned. The itype variable
! is a model specific integer that specifies the type of field (for
! instance temperature, zonal wind component, etc.). In low order
! models that have no notion of types of variables, this argument can
! be ignored. For applications in which only perfect model experiments
! with identity observations (i.e. only the value of a particular
! state variable is observed), this can be a NULL INTERFACE.

real(r8),            intent(in) :: x(:)
type(location_type), intent(in) :: location
integer,             intent(in) :: itype
real(r8),           intent(out) :: obs_val
integer,            intent(out) :: istatus

if ( .not. module_initialized ) call static_init_model

call error_handler(E_ERR,'model_interpolate','routine not tested',source, revision,revdate)
! This should be the result of the interpolation of a
! given kind (itype) of variable at the given location.
obs_val = MISSING_R8

! The return code for successful return should be 0. 
! Any positive number is an error.
! Negative values are reserved for use by the DART framework.
! Using distinct positive values for different types of errors can be
! useful in diagnosing problems.
istatus = 1

end subroutine model_interpolate



function get_model_time_step()
!------------------------------------------------------------------
!
! Returns the the time step of the model; the smallest increment
! in time that the model is capable of advancing the state in a given
! implementation. This interface is required for all applications.

type(time_type) :: get_model_time_step

if ( .not. module_initialized ) call static_init_model

call error_handler(E_ERR,'get_model_time_step','routine not tested',source, revision,revdate)
get_model_time_step = time_step

end function get_model_time_step

!------------------------------------------------------------------------
!> Given an index into the DART state vector, return the location and
!> (optionally) what kind of variable. This  optional argument kind
!> can be returned if the model has more than one type of field (for
!> instance temperature and zonal wind component). This interface is
!> required for all filter applications as it is required for computing
!> the distance between observations and state variables.

subroutine get_state_meta_data(index_in, location, var_type)

integer,             intent(in)            :: index_in
type(location_type), intent(out)           :: location
integer,             intent(out), optional :: var_type

real(r8) :: mylon,mylat,vloc
integer  :: iloc, jloc, kloc
integer  :: local_ind 

if ( .not. module_initialized ) call static_init_model

local_ind = index_in

kloc =  local_ind / (nx*ny) + 1
jloc = (local_ind - (kloc-1)*nx*ny)/nx + 1
iloc =  local_ind - (kloc-1)*nx*ny - (jloc-1)*nx + 1

write(*,*)'.. index_in ',index_in, ' and dereferences to ',iloc,jloc,kloc

! Now that we know the i,j,k we have to index the right set of
! coordinate arrays
mylon = lon(iloc,jloc)
mylat = lat(iloc,jloc)
vloc  = vcoord(kloc)

call error_handler(E_ERR,'get_state_meta_data','routine not tested',source, revision,revdate)
! these should be set to the actual location and obs kind
location = set_location(mylon, mylat, vloc, VERTISHEIGHT) ! meters
if (present(var_type)) var_type = 0  

end subroutine get_state_meta_data



subroutine end_model()
!------------------------------------------------------------------
!
! Does any shutdown and clean-up needed for model. Can be a NULL
! INTERFACE if the model has no need to clean up storage, etc.

if ( .not. module_initialized ) call static_init_model

call error_handler(E_ERR,'end_model','routine not tested',source, revision,revdate)
! good style ... perhaps you could deallocate stuff (from static_init_model?).
! deallocate(state_loc)

end subroutine end_model



function nc_write_model_atts( ncFileID ) result (ierr)
!------------------------------------------------------------------
! TJH 24 Oct 2006 -- Writes the model-specific attributes to a netCDF file.
!     This includes coordinate variables and some metadata, but NOT
!     the model state vector. We do have to allocate SPACE for the model
!     state vector, but that variable gets filled as the model advances.
!
! As it stands, this routine will work for ANY model, with no modification.
!
! The simplest possible netCDF file would contain a 3D field
! containing the state of 'all' the ensemble members. This requires
! three coordinate variables -- one for each of the dimensions 
! [model_size, ensemble_member, time]. A little metadata is useful, 
! so we can also create some 'global' attributes. 
! This is what is implemented here.
!
! Once the simplest case is working, this routine (and nc_write_model_vars)
! can be extended to create a more logical partitioning of the state vector,
! fundamentally creating a netCDF file with variables that are easily 
! plotted. The bgrid model_mod is perhaps a good one to view, keeping
! in mind it is complicated by the fact it has two coordinate systems. 
! There are stubs in this template, but they are only stubs.
!
! TJH 29 Jul 2003 -- for the moment, all errors are fatal, so the
! return code is always '0 == normal', since the fatal errors stop execution.
!
! assim_model_mod:init_diag_output uses information from the location_mod
!     to define the location dimension and variable ID. All we need to do
!     is query, verify, and fill ...
!
! Typical sequence for adding new dimensions,variables,attributes:
! NF90_OPEN             ! open existing netCDF dataset
!    NF90_redef         ! put into define mode 
!    NF90_def_dim       ! define additional dimensions (if any)
!    NF90_def_var       ! define variables: from name, type, and dims
!    NF90_put_att       ! assign attribute values
! NF90_ENDDEF           ! end definitions: leave define mode
!    NF90_put_var       ! provide values for variable
! NF90_CLOSE            ! close: save updated netCDF dataset

use typeSizes
use netcdf

integer, intent(in)  :: ncFileID      ! netCDF file identifier
integer              :: ierr          ! return value of function

integer :: nDimensions, nVariables, nAttributes, unlimitedDimID

integer :: StateVarDimID   ! netCDF pointer to state variable dimension (model size)
integer :: MemberDimID     ! netCDF pointer to dimension of ensemble    (ens_size)
integer :: TimeDimID       ! netCDF pointer to time dimension           (unlimited)

integer :: StateVarVarID   ! netCDF pointer to state variable coordinate array
integer :: StateVarID      ! netCDF pointer to 3D [state,copy,time] array

character(len=129)    :: errstring

! we are going to need these to record the creation date in the netCDF file.
! This is entirely optional, but nice.

character(len=8)      :: crdate      ! needed by F90 DATE_AND_TIME intrinsic
character(len=10)     :: crtime      ! needed by F90 DATE_AND_TIME intrinsic
character(len=5)      :: crzone      ! needed by F90 DATE_AND_TIME intrinsic
integer, dimension(8) :: values      ! needed by F90 DATE_AND_TIME intrinsic
character(len=NF90_MAX_NAME) :: str1

integer :: i

if ( .not. module_initialized ) call static_init_model

call error_handler(E_ERR,'nc_write_model_atts','routine not tested',source, revision,revdate)
!-------------------------------------------------------------------------------
! make sure ncFileID refers to an open netCDF file, 
! and then put into define mode.
!-------------------------------------------------------------------------------

ierr = -1 ! assume things go poorly

call nc_check(nf90_inquire(ncFileID,nDimensions,nVariables,nAttributes,unlimitedDimID), &
                     "nc_write_model_atts", "inquire")
call nc_check(nf90_redef(ncFileID), "nc_write_model_atts", "redef")

!-------------------------------------------------------------------------------
! We need the dimension ID for the number of copies/ensemble members, and
! we might as well check to make sure that Time is the Unlimited dimension. 
! Our job is create the 'model size' dimension.
!-------------------------------------------------------------------------------

call nc_check(nf90_inq_dimid(ncid=ncFileID, name="copy", dimid=MemberDimID), &
                            "nc_write_model_atts", "inq_dimid copy")
call nc_check(nf90_inq_dimid(ncid=ncFileID, name="time", dimid= TimeDimID), &
                            "nc_write_model_atts", "inq_dimid time")

if ( TimeDimID /= unlimitedDimId ) then
   write(errstring,*)"Time Dimension ID ",TimeDimID, &
                     " should equal Unlimited Dimension ID",unlimitedDimID
   call error_handler(E_ERR,"nc_write_model_atts", errstring, source, revision, revdate)
endif

!-------------------------------------------------------------------------------
! Define the model size / state variable dimension / whatever ...
!-------------------------------------------------------------------------------
call nc_check(nf90_def_dim(ncid=ncFileID, name="StateVariable",  &
                           len=model_size, dimid=StateVarDimID), &
                           "nc_write_model_atts", "def_dim state")

!-------------------------------------------------------------------------------
! Write Global Attributes 
!-------------------------------------------------------------------------------

call DATE_AND_TIME(crdate,crtime,crzone,values)
write(str1,'(''YYYY MM DD HH MM SS = '',i4,5(1x,i2.2))') &
                  values(1), values(2), values(3), values(5), values(6), values(7)

call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "creation_date" ,str1), &
                          "nc_write_model_atts", "put_att creation_date")
call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "model_source"  ,source), &
                          "nc_write_model_atts", "put_att model_source")
call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "model_revision",revision), &
                          "nc_write_model_atts", "put_att model_revision")
call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "model_revdate" ,revdate), &
                          "nc_write_model_atts", "put_att model_revdate")
call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "model","template"), &
                          "nc_write_model_atts", "put_att model")

!-------------------------------------------------------------------------------
! Here is the extensible part. The simplest scenario is to output the state vector,
! parsing the state vector into model-specific parts is complicated, and you need
! to know the geometry, the output variables (PS,U,V,T,Q,...) etc. We're skipping
! complicated part.
!-------------------------------------------------------------------------------

if ( output_state_vector_as_1D ) then

   !----------------------------------------------------------------------------
   ! Create a variable for the state vector
   !----------------------------------------------------------------------------

  ! Define the state vector coordinate variable and some attributes.
   call nc_check(nf90_def_var(ncid=ncFileID,name="StateVariable", xtype=NF90_INT, &
                              dimids=StateVarDimID, varid=StateVarVarID), &
                             "nc_write_model_atts", "def_var StateVariable")
   call nc_check(nf90_put_att(ncFileID, StateVarVarID,"long_name","State Variable ID"), &
                             "nc_write_model_atts", "put_att StateVariable long_name")
   call nc_check(nf90_put_att(ncFileID, StateVarVarID, "units",     "indexical"), &
                             "nc_write_model_atts", "put_att StateVariable units")
   call nc_check(nf90_put_att(ncFileID, StateVarVarID, "valid_range", (/ 1, model_size /)), &
                             "nc_write_model_atts", "put_att StateVariable valid_range")

   ! Define the actual (3D) state vector, which gets filled as time goes on ... 
   call nc_check(nf90_def_var(ncid=ncFileID, name="state", xtype=NF90_REAL, &
                 dimids = (/ StateVarDimID, MemberDimID, unlimitedDimID /), &
                 varid=StateVarID), "nc_write_model_atts", "def_var state")
   call nc_check(nf90_put_att(ncFileID, StateVarID, "long_name", "model state or fcopy"), &
                             "nc_write_model_atts", "put_att state long_name")

   ! Leave define mode so we can fill the coordinate variable.
   call nc_check(nf90_enddef(ncfileID),"nc_write_model_atts", "state_vector enddef")

   ! Fill the state variable coordinate variable
   call nc_check(nf90_put_var(ncFileID, StateVarVarID, (/ (i,i=1,model_size) /)), &
                                    "nc_write_model_atts", "put_var state")

else

   !----------------------------------------------------------------------------
   ! We need to process the prognostic variables.
   !----------------------------------------------------------------------------

   ! This block is a stub for something more complicated.
   ! Usually, the control for the execution of this block is a namelist variable.
   ! Take a peek at the bgrid model_mod.f90 for a (rather complicated) example.

   call nc_check(nf90_enddef(ncfileID), "nc_write_model_atts", "prognostic enddef")

endif

!-------------------------------------------------------------------------------
! Flush the buffer and leave netCDF file open
!-------------------------------------------------------------------------------
call nc_check(nf90_sync(ncFileID),"nc_write_model_atts", "sync")

ierr = 0 ! If we got here, things went well.

end function nc_write_model_atts



function nc_write_model_vars( ncFileID, statevec, copyindex, timeindex ) result (ierr)         
!------------------------------------------------------------------
! TJH 24 Oct 2006 -- Writes the model variables to a netCDF file.
!
! TJH 29 Jul 2003 -- for the moment, all errors are fatal, so the
! return code is always '0 == normal', since the fatal errors stop execution.
!
! For the lorenz_96 model, each state variable is at a separate location.
! that's all the model-specific attributes I can think of ...
!
! assim_model_mod:init_diag_output uses information from the location_mod
!     to define the location dimension and variable ID. All we need to do
!     is query, verify, and fill ...
!
! Typical sequence for adding new dimensions,variables,attributes:
! NF90_OPEN             ! open existing netCDF dataset
!    NF90_redef         ! put into define mode
!    NF90_def_dim       ! define additional dimensions (if any)
!    NF90_def_var       ! define variables: from name, type, and dims
!    NF90_put_att       ! assign attribute values
! NF90_ENDDEF           ! end definitions: leave define mode
!    NF90_put_var       ! provide values for variable
! NF90_CLOSE            ! close: save updated netCDF dataset

use typeSizes
use netcdf

integer,                intent(in) :: ncFileID      ! netCDF file identifier
real(r8), dimension(:), intent(in) :: statevec
integer,                intent(in) :: copyindex
integer,                intent(in) :: timeindex
integer                            :: ierr          ! return value of function

integer :: nDimensions, nVariables, nAttributes, unlimitedDimID

integer :: StateVarID

!-------------------------------------------------------------------------------
! make sure ncFileID refers to an open netCDF file, 
!-------------------------------------------------------------------------------

if ( .not. module_initialized ) call static_init_model
call error_handler(E_ERR,'nc_write_model_vars','routine not tested',source, revision,revdate)

ierr = -1 ! assume things go poorly

call nc_check(nf90_inquire(ncFileID,nDimensions,nVariables,nAttributes,unlimitedDimID), &
                          "nc_write_model_vars", "inquire")

if ( output_state_vector_as_1D ) then

   call nc_check(nf90_inq_varid(ncFileID, "state", StateVarID), &
                               "nc_write_model_vars", "inq_varid state" )
   call nc_check(nf90_put_var(ncFileID, StateVarID, statevec,  &
                              start=(/ 1, copyindex, timeindex /)), &
                             "nc_write_model_vars", "put_var state")                   

else

   !----------------------------------------------------------------------------
   ! We need to process the prognostic variables.
   !----------------------------------------------------------------------------

   ! This block is a stub for something more complicated.
   ! Usually, the control for the execution of this block is a namelist variable.
   ! Take a peek at the bgrid model_mod.f90 for a (rather complicated) example.
   !
   ! Generally, it is necessary to take the statevec and decompose it into 
   ! the separate prognostic variables. In this (commented out) example,
   ! global_Var is a user-defined type that has components like:
   ! global_Var%ps, global_Var%t, ... etc. Each of those can then be passed
   ! directly to the netcdf put_var routine. This may cause a huge storage
   ! hit, so large models may want to avoid the duplication if possible.

   ! call vector_to_prog_var(statevec, get_model_size(), global_Var)

   ! the 'start' array is crucial. In the following example, 'ps' is a 2D
   ! array, and the netCDF variable "ps" is a 4D array [lat,lon,copy,time]

   ! call nc_check(nf90_inq_varid(ncFileID, "ps", psVarID), &
   !                             "nc_write_model_vars",  "inq_varid ps")
   ! call nc_check(nf90_put_var( ncFileID, psVarID, global_Var%ps, &
   !                             start=(/ 1, 1, copyindex, timeindex /)), &
   !                            "nc_write_model_vars", "put_var ps")

endif

!-------------------------------------------------------------------------------
! Flush the buffer and leave netCDF file open
!-------------------------------------------------------------------------------

call nc_check(nf90_sync(ncFileID), "nc_write_model_vars", "sync")

ierr = 0 ! If we got here, things went well.

end function nc_write_model_vars



subroutine pert_model_state(state, pert_state, interf_provided)
!------------------------------------------------------------------
!
! Perturbs a model state for generating initial ensembles.
! The perturbed state is returned in pert_state.
! A model may choose to provide a NULL INTERFACE by returning
! .false. for the interf_provided argument. This indicates to
! the filter that if it needs to generate perturbed states, it
! may do so by adding an O(0.1) magnitude perturbation to each
! model state variable independently. The interf_provided argument
! should be returned as .true. if the model wants to do its own
! perturbing of states.  The returned pert_state should in any
! case be valid, since it will be read by filter even if 
! interf_provided is .false.

real(r8), intent(in)  :: state(:)
real(r8), intent(out) :: pert_state(:)
logical,  intent(out) :: interf_provided

if ( .not. module_initialized ) call static_init_model
call error_handler(E_ERR,'pert_model_state','routine not tested',source, revision,revdate)

pert_state      = state
interf_provided = .false.

end subroutine pert_model_state




subroutine ens_mean_for_model(ens_mean)
!------------------------------------------------------------------
! Not used in low-order models

real(r8), intent(in) :: ens_mean(:)

if ( .not. module_initialized ) call static_init_model
call error_handler(E_ERR,'ens_mean_for_model','routine not tested',source, revision,revdate)

end subroutine ens_mean_for_model


!==================================================================
! PUBLIC interfaces that aren't required by the DART code but are
! generally useful for other related utility programs.
! (less necessary for small models; generally used for larger models
! with predefined file formats and control structures.)
!==================================================================

!------------------------------------------------------------------
!> Reads the current time and state variables from a model data
!> file and packs them into a dart state vector.

subroutine get_state_vector(sv, model_time)

real(r8),         intent(inout)           :: sv(1:model_size)
type(time_type),  intent(out), optional   :: model_time
real(r8),allocatable                      :: pfbdata(:,:,:)
 
!

if ( .not. module_initialized ) call static_init_model

allocate(pfbdata(nx,ny,nz))

call pfread_var(parflow_file,pfbdata) 

sv(:) = reshape(pfbdata,(/ (nx*ny*nz) /))

model_time = parflow_time

if (debug > 0 .and. do_output()) write(*,*) '   ... data written to state_vector'
   
deallocate(pfbdata)

end subroutine get_state_vector

!------------------------------------------------------------------
!> Writes the current time and state variables from a dart state
!> vector (1d array) into a ncommas netcdf restart file.

subroutine dart_vector_to_model_file(state_vector, filename, statedate)

real(r8),         intent(in) :: state_vector(:)
character(len=*), intent(in) :: filename
type(time_type),  intent(in) :: statedate

if ( .not. module_initialized ) call static_init_model
call error_handler(E_ERR,'dart_vector_to_model_file','routine not tested',source, revision,revdate)

! code goes here

end subroutine dart_vector_to_model_file


!------------------------------------------------------------------
!> Reads the current time and state variables from a model data
!> file and packs them into a dart state vector.


function get_parflow_filename()

character(len=256) :: get_parflow_filename

if ( .not. module_initialized ) call static_init_model

get_parflow_filename = trim(parflow_file)

end function get_parflow_filename

!------------------------------------------------------------------
!> Reads pfidb ascii file to extract time and vertical co-ordinate

subroutine pfidb_read(filename)

character(len=*), intent(in)   :: filename

integer(kind=4)                :: nudat, izerr, iz
character(len=256)             :: errmsg, fmt
real(r8)                       :: pfb_dz(nz)
integer(kind=4)                :: yyyy, mm, dd, hh, mn, ss,  ts

! code starts here
  nudat   = get_unit()
  write(*,*) filename
  open(nudat, file=trim(filename),status='old')
 
  read(nudat,*,iostat=izerr) yyyy, mm, dd, ts 
  if (izerr < 0) call error_handler(E_ERR,'pfidb_read','error time', source, revision, revdate)
  if (debug > 0 .and. do_output()) write(*,*) '  ...parflow time ',yyyy, mm, dd, ts 

  hh = int(ts/3600._r8)
  mn = int(ts - hh*3600)
  ss = mod(ts, 60)  

  parflow_time = set_date(yyyy,mm,dd,hh,mn,ss) 

  do iz = 1, nz
    read(nudat,'(F5.2)',iostat=izerr) pfb_dz(iz)
    if (izerr < 0) then
       write(errmsg,*) 'ERROR READING pfidb_dz file ', iz
       call error_handler(E_ERR,'pfidb_read', errmsg, source, revision, revdate)
    endif
  end do

  do iz = nz, 1, -1
    if (iz == nz) then
      vcoord(iz) = 0.5_r8 * pfb_dz(iz) 
    else
      vcoord(iz) = vcoord(iz+1) + 0.5_r8 *(pfb_dz(iz+1) + pfb_dz(iz))
    end if
    if (debug > 0 .and. do_output()) write(*,'(A,1X,I2,1X,F5.2,1X,F6.3)')   '... pfidb ' ,&
                                         iz, pfb_dz(iz), vcoord(iz)
  end do

  close(nudat)

end subroutine pfidb_read

!------------------------------------------------------------------
!> Reads the oasis grids.nc file to extract parflow geo-location 

subroutine grid_read(filename)

! gpfl.lon(y_gpfl,x_gpfl)
! gpfl.lat(y_gpfl,x_gpfl)

character(len=*), intent(in)   :: filename
integer :: ncid, io, nlon, nlat, dimid(2), varid(2) 

io = nf90_open(filename, NF90_NOWRITE, ncid)
if (io /= nf90_noerr) call error_handler(E_ERR,'grid_read','ERR opening',source, revision,revdate) 

io = nf90_inq_dimid(ncid, "x_gpfl", dimid(1))
if (io /= nf90_noerr) call error_handler(E_ERR,'grid_read','ERR inq nlon',source, revision,revdate)

io = nf90_inq_dimid(ncid, "y_gpfl", dimid(2))
if (io /= nf90_noerr) call error_handler(E_ERR,'grid_read','ERR inq nlat',source, revision,revdate)

io = nf90_inq_varid(ncid, "gpfl.lon", varid(1))
if (io /= nf90_noerr) call error_handler(E_ERR,'grid_read','ERR inq lon',source, revision,revdate)

io = nf90_inq_varid(ncid, "gpfl.lat", varid(2))
if (io /= nf90_noerr) call error_handler(E_ERR,'grid_read','ERR inq lat',source, revision,revdate)

!
io = nf90_inquire_dimension(ncid, dimid(1), string1, nlon)
if (io /= nf90_noerr) call error_handler(E_ERR,'grid_read','ERR get nlon',source, revision,revdate)

io = nf90_inquire_dimension(ncid, dimid(2), string1, nlat)
if (io /= nf90_noerr) call error_handler(E_ERR,'grid_read','ERR get nlat',source, revision,revdate)

if (nlon /= nx .or. nlat /=ny) then
  write(*, *) 'Dimensions ...', nlon, nx, nlat, ny
  call error_handler(E_ERR,'grid_read','ERR dimension mismatch',source, revision,revdate)
end if

io = nf90_get_var(ncid, varid(1), lon)
if (io /= nf90_noerr) call error_handler(E_ERR,'grid_read','ERR get lon',source, revision,revdate)

io = nf90_get_var(ncid, varid(2), lat)
if (io /= nf90_noerr) call error_handler(E_ERR,'grid_read','ERR get lat',source, revision,revdate)

where(lon <   0.0_r8) lon = lon + 360.0_r8
where(lat < -90.0_r8) lat = -90.0_r8
where(lat >  90.0_r8) lat =  90.0_r8

io = nf90_close(ncid)

end subroutine grid_read
!------------------------------------------------------------------
!> Reads pbf dimensions. based on tr32-z4-tools work of P. Shrestha

subroutine pfread_dim(filename)

character(len=*), intent(in)   :: filename

real(r8)                       :: x1, y1 ,z1
integer(kind=4)                :: nudat, izerr 
character(len=256)             :: errmsg

! code starts here
  nudat   = get_unit() 
  open(nudat,file=trim(filename),form='unformatted',access='stream' , &
                    convert='BIG_ENDIAN',status='old')         ! gfortran

  !read in header infor
  read(nudat, iostat=izerr) x1 !X
  if (izerr /= 0) then
    errmsg   = "unable to read X"
    call error_handler(E_ERR,'pfread_dim',errmsg,source,revision,revdate)
  endif
  read(nudat, iostat=izerr) y1 !Y
  if (izerr /= 0) then
    errmsg   = "unable to read Y"
    call error_handler(E_ERR,'pfread_dim',errmsg,source,revision,revdate)
  endif
  read(nudat, iostat=izerr) z1 !Z
  if (izerr /= 0) then
    errmsg   = "unable to read Z"
    call error_handler(E_ERR,'pfread_dim',errmsg,source,revision,revdate)
  endif

  read(nudat, iostat=izerr) nx !NX
  if (izerr /= 0) then
    errmsg   = "unable to read NX"
    call error_handler(E_ERR,'pfread_dim',errmsg,source,revision,revdate)
  endif
  read(nudat, iostat=izerr) ny !NY
  if (izerr /= 0) then
    errmsg   = "unable to read NY"
    call error_handler(E_ERR,'pfread_dim',errmsg,source,revision,revdate)
  endif
  read(nudat, iostat=izerr) nz !NZ
  if (izerr /= 0) then
    errmsg   = "unable to read NZ"
    call error_handler(E_ERR,'pfread_dim',errmsg,source,revision,revdate)
  endif

  read(nudat, iostat=izerr) dx !dX
  if (izerr /= 0) then
    errmsg   = "unable to read dx"
    call error_handler(E_ERR,'pfread_dim',errmsg,source,revision,revdate)
  endif
  read(nudat, iostat=izerr) dy !dY
  if (izerr /= 0) then
    errmsg   = "unable to read dy"
    call error_handler(E_ERR,'pfread_dim',errmsg,source,revision,revdate)
  endif
  read(nudat, iostat=izerr) dz !dZ
  if (izerr /= 0) then
    errmsg   = "unable to read dz"
    call error_handler(E_ERR,'pfread_dim',errmsg,source,revision,revdate)
  endif

  close(nudat)
end subroutine pfread_dim

subroutine pfread_var(filename,pfvar)
!------------------------------------------------------------------
! Reads pbf dimensions. based on tr32-z4-tools work of P. Shrestha

character(len=*), intent(in)   :: filename
real(r8),         intent(out)  ::  &
                   pfvar(nx,ny,nz)       ! ParFlow pressure files 

real(r8)                       :: x1, y1, z1 

real(r8)                       ::  &
                   dummyRes              ! Grid Resolution in m

integer(kind=4)                :: nudat, dummy, izerr
integer(kind=4)                ::  &
                   i, j, k,        &
                   ix, iy, iz,     &
                   is, ns,         &
                   rx, ry, rz,     &
                   nnx, nny, nnz
character(len=256)             :: errmsg

! code starts here
  nudat   = get_unit()
  open(nudat,file=trim(filename),form='unformatted',access='stream' , &
                    convert='BIG_ENDIAN',status='old')         ! gfortran

  !read in header infor
  read(nudat, iostat=izerr) x1 !X
  if (izerr /= 0) then
    errmsg   = "unable to read X"
    call error_handler(E_ERR,'pfread_var',errmsg,source,revision,revdate)
  endif
  read(nudat, iostat=izerr) y1 !Y
  if (izerr /= 0) then
    errmsg   = "unable to read Y"
    call error_handler(E_ERR,'pfread_var',errmsg,source,revision,revdate)
  endif
  read(nudat, iostat=izerr) z1 !Z
  if (izerr /= 0) then
    errmsg   = "unable to read Z"
    call error_handler(E_ERR,'pfread_var',errmsg,source,revision,revdate)
  endif

  read(nudat, iostat=izerr) dummy !NX
  if (izerr /= 0) then
    errmsg   = "unable to read NX"
    call error_handler(E_ERR,'pfread_var',errmsg,source,revision,revdate)
  endif
  read(nudat, iostat=izerr) dummy !NY
  if (izerr /= 0) then
    errmsg   = "unable to read NY"
    call error_handler(E_ERR,'pfread_var',errmsg,source,revision,revdate)
  endif
  read(nudat, iostat=izerr) dummy !NZ
  if (izerr /= 0) then
    errmsg   = "unable to read NZ"
    call error_handler(E_ERR,'pfread_var',errmsg,source,revision,revdate)
  endif

  read(nudat, iostat=izerr) dummyRes !dX
  if (izerr /= 0) then
    errmsg   = "unable to read dx"
    call error_handler(E_ERR,'pfread_var',errmsg,source,revision,revdate)
  endif
  read(nudat, iostat=izerr) dummyRes !dY
  if (izerr /= 0) then
    errmsg   = "unable to read dy"
    call error_handler(E_ERR,'pfread_var',errmsg,source,revision,revdate)
  endif
  read(nudat, iostat=izerr) dummyRes !dZ
  if (izerr /= 0) then
    errmsg   = "unable to read dz"
    call error_handler(E_ERR,'pfread_var',errmsg,source,revision,revdate)
  endif

  read(nudat, iostat=izerr) ns !num_subgrids
  if (izerr /= 0) then
    errmsg   = "unable to read ns"
    call error_handler(E_ERR,'pfread_var',errmsg,source,revision,revdate)
  endif
!
! Start loop over number of subgrids
  do is = 0, ns-1
! Start: reading subgrid spatial information
! ix, iy, iz
   read(nudat, iostat=izerr) ix
   if (izerr /= 0) then
     errmsg   = "unable to read ix"
     call error_handler(E_ERR,'pfread_var',errmsg,source,revision,revdate)
   endif
   read(nudat, iostat=izerr) iy
   if (izerr /= 0) then
     errmsg   = "unable to read ix"
     call error_handler(E_ERR,'pfread_var',errmsg,source,revision,revdate)
   endif
   read(nudat, iostat=izerr) iz
   if (izerr /= 0) then
     errmsg   = "unable to read ix"
     call error_handler(E_ERR,'pfread_var',errmsg,source,revision,revdate)
   endif

! nnx. nny, nnz 
   read(nudat, iostat=izerr) nnx
   if (izerr /= 0) then
     errmsg   = "unable to read nnx"
     call error_handler(E_ERR,'pfread_var',errmsg,source,revision,revdate)
   endif
   read(nudat, iostat=izerr) nny
   if (izerr /= 0) then
     errmsg   = "unable to read nny"
     call error_handler(E_ERR,'pfread_var',errmsg,source,revision,revdate)
   endif
   read(nudat, iostat=izerr) nnz
   if (izerr /= 0) then
     errmsg   = "unable to read nnz"
     call error_handler(E_ERR,'pfread_var',errmsg,source,revision,revdate)
   endif

! rx,ry,rz
   read(nudat, iostat=izerr) rx
   if (izerr /= 0) then
     errmsg   = "unable to read rx"
     call error_handler(E_ERR,'pfread_var',errmsg,source,revision,revdate)
   endif
   read(nudat, iostat=izerr) ry
   if (izerr /= 0) then
     errmsg   = "unable to read ry"
     call error_handler(E_ERR,'pfread_var',errmsg,source,revision,revdate)
   endif
   read(nudat, iostat=izerr) rz
   if (izerr /= 0) then
     errmsg   = "unable to read rz"
     call error_handler(E_ERR,'pfread_var',errmsg,source,revision,revdate)
   endif

! Start reading in the data finally for each individual subgrid
   ! Start: Read in data from each individual subgrid
   do  k=iz +1 , iz + nnz
   do  j=iy +1 , iy + nny
   do  i=ix +1 , ix + nnx
     read(nudat, iostat=izerr) pfvar(i,j,k)
   end do
   end do
   end do
   ! End: Read in data from each individual subgrid

  end do
! End: loop over number of sub grids

  close(nudat)
  return
end subroutine pfread_var

subroutine pfwrite_var(filename,nx,ny,nz,dx,dy,dz)
!------------------------------------------------------------------
! Reads pbf dimensions. based on tr32-z4-tools work of P. Shrestha

character(len=*), intent(in)   :: filename
integer(kind=4),  intent(out)  ::  &
                   nx,             &     ! Longitude dimension
                   ny,             &     ! Latitude dimension
                   nz                    ! Vertical dimension
real(r8),         intent(out)  ::  &
                   dx,             &     ! Grid Resolution in m
                   dy,             &     ! Grid Resoultion in m
                   dz                    ! Grid Resolution scale, check parflow namelist

real(r8)                       :: x1, y1, z1 
integer(kind=4)                :: nudat, izerr
character(len=256)             :: errmsg

! code starts here
  nudat   = get_unit()
  open(nudat,file=trim(filename),form='unformatted',access='stream' , &
                    convert='BIG_ENDIAN',status='old')         ! gfortran

  !read in header infor
  read(nudat, iostat=izerr) x1 !X
  if (izerr /= 0) then
    errmsg   = "unable to read X"
    call error_handler(E_ERR,'pfread_dim',errmsg,source,revision,revdate)
  endif
  read(nudat, iostat=izerr) y1 !Y
  if (izerr /= 0) then
    errmsg   = "unable to read Y"
    call error_handler(E_ERR,'pfread_dim',errmsg,source,revision,revdate)
  endif
  read(nudat, iostat=izerr) z1 !Z
  if (izerr /= 0) then
    errmsg   = "unable to read Z"
    call error_handler(E_ERR,'pfread_dim',errmsg,source,revision,revdate)
  endif

  read(nudat, iostat=izerr) nx !NX
  if (izerr /= 0) then
    errmsg   = "unable to read NX"
    call error_handler(E_ERR,'pfread_dim',errmsg,source,revision,revdate)
  endif
  read(nudat, iostat=izerr) ny !NY
  if (izerr /= 0) then
    errmsg   = "unable to read NY"
    call error_handler(E_ERR,'pfread_dim',errmsg,source,revision,revdate)
  endif
  read(nudat, iostat=izerr) nz !NZ
  if (izerr /= 0) then
    errmsg   = "unable to read NZ"
    call error_handler(E_ERR,'pfread_dim',errmsg,source,revision,revdate)
  endif

  read(nudat, iostat=izerr) dx !dX
  if (izerr /= 0) then
    errmsg   = "unable to read dx"
    call error_handler(E_ERR,'pfread_dim',errmsg,source,revision,revdate)
  endif
  read(nudat, iostat=izerr) dy !dY
  if (izerr /= 0) then
    errmsg   = "unable to read dy"
    call error_handler(E_ERR,'pfread_dim',errmsg,source,revision,revdate)
  endif
  read(nudat, iostat=izerr) dz !dZ
  if (izerr /= 0) then
    errmsg   = "unable to read dz"
    call error_handler(E_ERR,'pfread_dim',errmsg,source,revision,revdate)
  endif

  close(nudat)
  return
end subroutine pfwrite_var
!===================================================================
! End of model_mod
!===================================================================
end module model_mod

! <next few lines under version control, do not edit>
! $URL: https://proxy.subversion.ucar.edu/DAReS/DART/releases/Lanai/models/template/model_mod.f90 $
! $Id: model_mod.f90 6256 2013-06-12 16:19:10Z thoar $
! $Revision: 6256 $
! $Date: 2013-06-12 18:19:10 +0200 (Wed, 12 Jun 2013) $
