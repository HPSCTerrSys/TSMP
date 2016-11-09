! This code is not necessarily under the DART copyright ...
!

! DART $Id: $

module model_mod

! This module contains routines to work with ParFlow Binary data 
! in the DART framework

! Based on DART template directory

! Modules that are absolutely required for use are listed
use        types_mod, only : r8, obstypelength, MISSING_R8

use time_manager_mod, only : time_type, set_time, set_date, get_date, get_time, &
                             print_time, print_date, set_calendar_type

use     location_mod, only : location_type,      get_close_maxdist_init,        &
                             get_close_obs_init, get_close_obs, set_location,   &
                             get_location,                       &
                             vert_is_height,       VERTISHEIGHT

use    utilities_mod, only : register_module, error_handler, nc_check, &
                             get_unit, open_file, close_file, E_ERR, E_MSG, &
                             nmlfileunit, do_output, do_nml_file, do_nml_term,  &
                             find_namelist_in_file, file_exist, check_namelist_read

use     obs_kind_mod, only : KIND_SOIL_MOISTURE,          &
                             paramname_length,            &
                             get_raw_obs_kind_index,      &
                             get_raw_obs_kind_name

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
          write_state_times,      &
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
integer, parameter :: max_state_variables = 2   !ParFlow has two ouputs 

type progvartype
   private
   character(len=NF90_MAX_NAME) :: varname
   character(len=NF90_MAX_NAME) :: long_name
   character(len=NF90_MAX_NAME) :: units
   character(len=obstypelength) :: dimnames(NF90_MAX_VAR_DIMS)
   integer  :: dimlens(NF90_MAX_VAR_DIMS)
   integer  :: numdims
   integer  :: numEW
   integer  :: numNS
   integer  :: numZ
   integer  :: index1      ! location in dart state vector of first occurrence
   integer  :: varsize     ! prod(dimlens(1:numdims))
   integer  :: rangeRestricted
   integer  :: pfb_kind
   real(r8) :: minvalue
   real(r8) :: maxvalue
   character(len=256) :: kind_string
   logical  :: update
end type progvartype

type(progvartype), dimension(max_state_variables)     :: progvar

real(r8), allocatable ::   lon(:,:)      ! longitude degrees_east
real(r8), allocatable ::   lat(:,:)      ! latitude  degrees_north
real(r8), allocatable ::   vcoord(:)     ! vertical co-ordinate in m
type(time_type)       ::   parflow_time  ! parflow ouptut time

! EXAMPLE: perhaps a namelist here 
character(len=256) :: parflow_press_file           = 'parflow_press_file'
character(len=256) :: parflow_satur_file           = 'parflow_satur_file'
character(len=256) :: pfidb_file                   = 'pfidb_file'
character(len=256) :: grid_file                    = 'grid_file'
character(len=256) :: clm_file                     = 'clmoash0_file'
logical            :: output_1D_state_vector       = .false.
integer            :: assimilation_period_days     = 0
integer            :: assimilation_period_seconds  = 21400
integer            :: debug = 0 

namelist /model_nml/                &
      parflow_press_file,           &
      parflow_satur_file,           &
      pfidb_file,                   &
      grid_file,                    &
      clm_file,                     &
      output_1D_state_vector,       &
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
integer  :: ivar

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
call pfread_dim(parflow_press_file)

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
if (debug > 1 .and. do_output()) then
   write(*, '(A,2(1X,F9.2))') '    ...parflow lon', minval(lon), maxval(lon)
   write(*, '(A,2(1X,F9.2))') '    ...parflow lat', minval(lat), maxval(lat)
   write(*, *) '-------------------------------'
end if

!CPS call error_handler(E_ERR,'static_init_model','routine not tested',source, revision,revdate)
! TODO, var_type could switch between 1 and 2 for pressure and sat
ivar = 1
progvar(ivar)%varname         = "psi"
progvar(ivar)%long_name       = "Pressure Head"
progvar(ivar)%units           = "m"
progvar(ivar)%numdims         = 3
progvar(ivar)%dimnames(1)     = 'pfl_lon'
progvar(ivar)%dimnames(2)     = 'pfl_lat'
progvar(ivar)%dimnames(3)     = 'level'
progvar(ivar)%dimlens(1)      = nx
progvar(ivar)%dimlens(2)      = ny
progvar(ivar)%dimlens(3)      = nz
progvar(ivar)%pfb_kind        = KIND_SOIL_MOISTURE   !Need pedotransfer function 
progvar(ivar)%index1          = 1

! Create storage for locations
allocate(state_loc(model_size))

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
write(string1,*) 'input.nml:start_from_restart cannot be FALSE'
call error_handler(E_ERR,'init_conditions',string1,source, revision,revdate)
x = MISSING_R8

return
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


subroutine write_state_times(iunit, statetime)

integer,                   intent(in) :: iunit
type(time_type),           intent(in) :: statetime

integer           :: iyear, imonth, iday, ihour, imin, isec 
integer           :: ndays, nhours, nmins, nsecs,nfreq
type(time_type)   :: interval

call get_date(statetime, iyear, imonth, iday, ihour, imin, isec)
write(iunit, '(''defaultInitDate '',I4.4,2(''-'',I2.2),1x,i2.2)') iyear, imonth, iday, ihour

return 
end subroutine write_state_times

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

write(string1,*)'input.nml:start_from_restart cannot be FALSE'
call error_handler(E_ERR,'init_time',string1,source, revision,revdate)
! for now, just set to 0
time = set_time(0,0)

end subroutine init_time

!------------------------------------------------------------------------
!>Given a state vector, a location, and a model state variable type,
!>interpolates the state variable field to that location and returns
!>the value in obs_val. The istatus variable should be returned as
!>0 unless there is some problem in computing the interpolation in
!>which case an alternate value should be returned. The obs_type variable
!>is a model specific integer that specifies the type of field .

subroutine model_interpolate(x, location, obs_type, interp_val, istatus)
!------------------------------------------------------------------
!
! Error codes:
! istatus = 99 : unknown error
! istatus = 10 : observation type is not in state vector
! istatus = 15 : observation lies outside the model domain (horizontal)
! istatus = 16 : observation lies outside the model domain (vertical)
! istatus = 19 : observation vertical coordinate is not supported

real(r8),            intent(in) :: x(:)
type(location_type), intent(in) :: location
integer,             intent(in) :: obs_type
real(r8),           intent(out) :: interp_val
integer,            intent(out) :: istatus

integer, parameter :: NO_KIND_IN_STATE     = 10
integer, parameter :: VERTICAL_UNSUPPORTED = 19

! Local Variables

integer  :: pfb_kind
integer  :: i, ivar
integer  :: geo_inds(4), hgt_inds(2)
real(r8) :: geo_wgts(4), hgt_wgt
real(r8) :: point_coords(1:3)
real(r8) :: geo_lon, geo_lat, geo_hgt
real(r8) :: interp_h_var(2)

! Variable Definition Complete

if ( .not. module_initialized ) call static_init_model

interp_val = MISSING_R8     ! the DART bad value flag
istatus = 99                ! unknown error

pfb_kind = obs_type

if ( pfb_kind < 0 ) then
   interp_val = x(-1*pfb_kind)
   istatus = 0
   return
endif

! Determine if this pfb_kind exists in ParFlow state vector
! Unnecessary but who knows if ivar will increase in future

ivar = -1
FoundIt: do i = 1, max_state_variables     ! we do not yet define nfields in input_nml
   if (progvar(i)%pfb_kind == pfb_kind) then
      ivar = i
      exit FoundIt
   endif
enddo FoundIt

! Get the geo-location
point_coords(1:3) = get_location(location)
geo_lon = point_coords(1) ! degrees East
geo_lat = point_coords(2) ! degrees North
geo_hgt  = point_coords(3) ! depth in meters 

call get_corners(geo_lon, geo_lat, nx, ny, lon(1:nx,1), lat(1,1:ny), geo_inds, geo_wgts, istatus) 

if (istatus /= 0) return

call get_level_indices(geo_hgt, hgt_inds, hgt_wgt, istatus)

if (istatus /= 0) return

! Horizontal interpolation at two levels given by hgt_inds
call horizontal_interpolate(x, ivar, hgt_inds, geo_inds, geo_wgts, interp_h_var, istatus)

interp_val = interp_h_var(1)*hgt_wgt + interp_h_var(2)*(1.0_r8 - hgt_wgt)

istatus = 0

return 

end subroutine model_interpolate



function get_model_time_step()
!------------------------------------------------------------------
!
! Returns the the time step of the model; the smallest increment
! in time that the model is capable of advancing the state in a given
! implementation. This interface is required for all applications.

type(time_type) :: get_model_time_step

if ( .not. module_initialized ) call static_init_model

get_model_time_step = time_step

return
!call error_handler(E_ERR,'get_model_time_step','routine not tested',source, revision,revdate)
!get_model_time_step = time_step

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
integer  :: ivar

if ( .not. module_initialized ) call static_init_model

local_ind = index_in - 1    !CPS offset for algorithm below

kloc =  local_ind / (nx*ny) + 1
jloc = (local_ind - (kloc-1)*nx*ny)/nx + 1
iloc =  local_ind - (kloc-1)*nx*ny - (jloc-1)*nx + 1

if ((debug > 5) .and. do_output()) then
  write(*,*)'.. index_in ',index_in, ' and dereferences to i,j,k ',iloc,jloc,kloc
end if
! Now that we know the i,j,k we have to index the right set of
! coordinate arrays
mylon = lon(iloc,jloc)
mylat = lat(iloc,jloc)
vloc  = vcoord(kloc)

!CPS call error_handler(E_ERR,'get_state_meta_data','routine not tested',source, revision,revdate)
! these should be set to the actual location and obs kind
location = set_location(mylon, mylat, vloc, VERTISHEIGHT) ! meters

if (present(var_type)) var_type = 0  

end subroutine get_state_meta_data



subroutine end_model()
!------------------------------------------------------------------
!
! Does any shutdown and clean-up needed for model. Can be a NULL
! INTERFACE if the model has no need to clean up storage, etc.

deallocate(lon,lat)
deallocate(vcoord)
deallocate(state_loc)
!if ( allocated(ens_mean) ) deallocate(ens_mean)

end subroutine end_model


!------------------------------------------------------------------------
!> 

function nc_write_model_atts( ncFileID ) result (ierr)
!------------------------------------------------------------------
! TJH 24 Oct 2006 -- Writes the model-specific attributes to a netCDF file.
!     This includes coordinate variables and some metadata, but NOT
!     the model state vector. We do have to allocate SPACE for the model
!     state vector, but that variable gets filled as the model advances.
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
integer :: vcoordVarID
integer :: VarID

integer :: StateVarVarID   ! netCDF pointer to state variable coordinate array
integer :: StateVarID      ! netCDF pointer to 3D [state,copy,time] array

integer :: lonDimID
integer :: latDimID
integer :: levelDimID

character(len=129)    :: errstring

character(len=256)   :: filename

! we are going to need these to record the creation date in the netCDF file.
! This is entirely optional, but nice.

character(len=8)      :: crdate      ! needed by F90 DATE_AND_TIME intrinsic
character(len=10)     :: crtime      ! needed by F90 DATE_AND_TIME intrinsic
character(len=5)      :: crzone      ! needed by F90 DATE_AND_TIME intrinsic
integer, dimension(8) :: values      ! needed by F90 DATE_AND_TIME intrinsic
character(len=NF90_MAX_NAME) :: str1

integer, dimension(NF90_MAX_VAR_DIMS) :: mydimids
character(len=NF90_MAX_NAME) :: varname

integer :: io, ndims, ivar, i

if ( .not. module_initialized ) call static_init_model

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

if ( output_1D_state_vector ) then

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

   io = nf90_def_dim(ncid=ncFileID, name='pfl_lon', len=nx, dimid = lonDimID)
   call nc_check(io, 'nc_write_model_atts', 'pfl_lon def_dim '//trim(filename))

   io = nf90_def_dim(ncid=ncFileID, name='pfl_lat', len=ny, dimid = latDimID)
   call nc_check(io, 'nc_write_model_atts', 'pfl_lat def_dim '//trim(filename))

   io = nf90_def_dim(ncid=ncFileID, name='level', len=nz, dimid = levelDimID)
   call nc_check(io, 'nc_write_model_atts', 'level def_dim '//trim(filename))

! Standard Grid Longitudes
   call nc_check(nf90_def_var(ncFileID,name='lon', xtype=nf90_real, &
                 dimids=(/ lonDimID, latDimID /), varid=VarID),&
                 'nc_write_model_atts', 'lon def_var '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'long_name', 'longitude'), &
                 'nc_write_model_atts', 'lon long_name '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'cartesian_axis', 'X'),  &
                 'nc_write_model_atts', 'lon cartesian_axis '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'units', 'degrees_east'), &
                 'nc_write_model_atts', 'lon units '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'valid_range', (/ 0.0_r8, 360.0_r8 /)), &
                 'nc_write_model_atts', 'lon valid_range '//trim(filename))

! Standard Grid Latitudes
   call nc_check(nf90_def_var(ncFileID,name='lat', xtype=nf90_real, &
                 dimids=(/ lonDimID, latDimID /), varid=VarID),&
                 'nc_write_model_atts', 'lat def_var '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'long_name', 'latitudes of grid'), &
                 'nc_write_model_atts', 'lat long_name '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'cartesian_axis', 'Y'),  &
                 'nc_write_model_atts', 'lat cartesian_axis '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'units', 'degrees_east'), &
                 'nc_write_model_atts', 'lat units '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'valid_range', (/ -90.0_r8, 90.0_r8 /)), &
                 'nc_write_model_atts', 'lat valid_range '//trim(filename))
! vcoord

   io = nf90_def_var(ncFileID,name='vcoord', xtype=nf90_real, &
                 dimids=(/ levelDimID /), varid=vcoordVarID)
   call nc_check(io, 'nc_write_model_atts', 'vcoord def_var '//trim(filename))
   io = nf90_put_att(ncFileID, vcoordVarID,'long_name', 'vertical height')
   call nc_check(io, 'nc_write_model_atts', 'vcoord long_name '//trim(filename))

   io = nf90_put_att(ncFileID, vcoordVarID, 'units', 'm ')
   call nc_check(io, 'nc_write_model_atts', 'vcoord units '//trim(filename))

   !----------------------------------------------------------------------------
   ! Create the (empty) Prognostic Variables and the Attributes
   !----------------------------------------------------------------------------

   do ivar=1, 1 

      varname = trim(progvar(ivar)%varname)
      string1 = trim(filename)//' '//trim(varname)

      ! match shape of the variable to the dimension IDs
      call define_var_dims(ivar, ncFileID, MemberDimID, unlimitedDimID, ndims, mydimids)

      io = nf90_def_var(ncid=ncFileID, name=trim(varname), xtype=nf90_double, &
                    dimids = mydimids(1:ndims), varid=VarID)
      call nc_check(io, 'nc_write_model_atts', trim(string1)//' def_var' )

      io = nf90_put_att(ncFileID, VarID, 'long_name', trim(progvar(ivar)%long_name))
      call nc_check(io, 'nc_write_model_atts', trim(string1)//' put_att long_name' )

      io = nf90_put_att(ncFileID, VarID, 'units', trim(progvar(ivar)%units))
      call nc_check(io, 'nc_write_model_atts', trim(string1)//' put_att units' )

   enddo

   call nc_check(nf90_enddef(ncfileID), "nc_write_model_atts", "prognostic enddef")

   !----------------------------------------------------------------------------
   ! Fill the coordinate variables - reshape the 1D arrays to the 2D shape
   !----------------------------------------------------------------------------

   io = nf90_inq_varid(ncFileID, 'lon', VarID)
   call nc_check(io, 'nc_write_model_atts', 'lon inq_varid '//trim(filename))
   io = nf90_put_var(ncFileID, VarID, lon)
   call nc_check(io, 'nc_write_model_atts', 'lon put_var '//trim(filename))

   io = nf90_inq_varid(ncFileID, 'lat', VarID)
   call nc_check(io, 'nc_write_model_atts', 'lat inq_varid '//trim(filename))
   io = nf90_put_var(ncFileID, VarID, lat )
   call nc_check(io, 'nc_write_model_atts', 'lat put_var '//trim(filename))

   io = nf90_inq_varid(ncFileID, 'vcoord', VarID)
   call nc_check(io, 'nc_write_model_atts', 'vcoord inq_varid '//trim(filename))
   io = nf90_put_var(ncFileID, VarID, vcoord )
   call nc_check(io, 'nc_write_model_atts', 'vcoord put_var '//trim(filename))
endif

!-------------------------------------------------------------------------------
! Flush the buffer and leave netCDF file open
!-------------------------------------------------------------------------------

call nc_check(nf90_sync(ncFileID),"nc_write_model_atts", "sync")

ierr = 0 ! If we got here, things went well.
!CPS call error_handler(E_ERR,'nc_write_model_atts','routine not tested',source, revision,revdate)

end function nc_write_model_atts


!------------------------------------------------------------------------
!>
!> Writes the model variables to a netCDF file.
!> All errors are fatal, so the return code is always '0 == normal'.

function nc_write_model_vars( ncFileID, statevec, copyindex, timeindex ) result (ierr)         
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

character(len=NF90_MAX_NAME) :: varname
integer ::  dimIDs(NF90_MAX_VAR_DIMS)
integer :: ncstart(NF90_MAX_VAR_DIMS)
integer :: nccount(NF90_MAX_VAR_DIMS)

integer :: io, i, ivar, VarID, ndims, dimlen
integer :: TimeDimID, CopyDimID

character(len=256) :: filename

integer :: StateVarID

!-------------------------------------------------------------------------------
! make sure ncFileID refers to an open netCDF file, 
!-------------------------------------------------------------------------------

if ( .not. module_initialized ) call static_init_model

ierr = -1 ! assume things go poorly

write(filename,*) 'ncFileID', ncFileID

! make sure ncFileID refers to an open netCDF file,

io = nf90_inq_dimid(ncFileID, 'copy', dimid=CopyDimID)
call nc_check(io, 'nc_write_model_vars', 'inq_dimid copy '//trim(filename))

io = nf90_inq_dimid(ncFileID, 'time', dimid=TimeDimID)
call nc_check(io, 'nc_write_model_vars', 'inq_dimid time '//trim(filename))

if ( output_1D_state_vector ) then

   io = nf90_inq_varid(ncFileID, 'state', VarID)
   call nc_check(io, 'nc_write_model_vars', 'state inq_varid '//trim(filename))

   io = nf90_put_var(ncFileID,VarID,statevec,start=(/1,copyindex,timeindex/))
   call nc_check(io, 'nc_write_model_vars', 'state put_var '//trim(filename))

else

   !----------------------------------------------------------------------------
   ! We need to process the prognostic variables.
   !----------------------------------------------------------------------------
   do ivar = 1,1

      varname = trim(progvar(ivar)%varname)
      string2 = trim(filename)//' '//trim(varname)

      call nc_check(nf90_inq_varid(ncFileID, varname, VarID), &
            'nc_write_model_vars', 'inq_varid '//trim(string2))

      call nc_check(nf90_inquire_variable(ncFileID,VarID,dimids=dimIDs,ndims=ndims), &
            'nc_write_model_vars', 'inquire '//trim(string2))

      ncstart = 1   ! These are arrays, actually
      nccount = 1
      DimCheck : do i = 1,progvar(ivar)%numdims

         write(string1,'(a,i2,A)') 'inquire dimension ',i,trim(string2)
         call nc_check(nf90_inquire_dimension(ncFileID, dimIDs(i), len=dimlen), &
               'nc_write_model_vars', trim(string1))

         if (progvar(ivar)%dimnames(i) == 'time') cycle DimCheck

         if ( dimlen /= progvar(ivar)%dimlens(i) ) then
            write(string1,*)trim(string2),' dim/dimlen ',i,dimlen, &
                            ' not ',progvar(ivar)%dimlens(i)
            write(string2,*)' but it should be.'
            call error_handler(E_ERR, 'nc_write_model_vars', trim(string1), &
                            source, revision, revdate, text2=trim(string2))
         endif

         nccount(i) = dimlen

      enddo DimCheck

      where(dimIDs == CopyDimID) ncstart = copyindex
      where(dimIDs == CopyDimID) nccount = 1
      where(dimIDs == TimeDimID) ncstart = timeindex
      where(dimIDs == TimeDimID) nccount = 1

      if ((debug > 10) .and. do_output()) then
         write(*,*)'nc_write_model_vars '//trim(varname)//' start is ',ncstart(1:ndims)
         write(*,*)'nc_write_model_vars '//trim(varname)//' count is ',nccount(1:ndims)
      endif

      ! revelation - no need to reshape the state vector before nf90_put_var

      call nc_check(nf90_put_var(ncFileID, VarID, statevec,          &
                   start = ncstart(1:ndims), count=nccount(1:ndims)), &
                   'nc_write_model_vars', 'put_var '//trim(string2))
   enddo


endif

!-------------------------------------------------------------------------------
! Flush the buffer and leave netCDF file open
!-------------------------------------------------------------------------------

call nc_check(nf90_sync(ncFileID), "nc_write_model_vars", "sync")

ierr = 0 ! If we got here, things went well.

!CPS call error_handler(E_ERR,'nc_write_model_vars','routine not tested',source, revision,revdate)

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

!not needed
! allocate(ens_mean(1:model_size))
! ens_mean(:) = filter_ens_mean(:)
!call error_handler(E_ERR,'ens_mean_for_model','routine not tested',source, revision,revdate)

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
integer                                   :: ncid 
!

if ( .not. module_initialized ) call static_init_model

allocate(pfbdata(nx,ny,nz))

call pfread_var(parflow_press_file,pfbdata) 

sv(:) = reshape(pfbdata,(/ (nx*ny*nz) /))

if ( .not. file_exist(clm_file) ) then 
   write(string1,*) 'cannot open file ', trim(clm_file),' for reading.'
   call error_handler(E_ERR,'restart_file_to_sv',string1,source,revision,revdate)
endif

call nc_check(nf90_open(trim(clm_file), NF90_NOWRITE, ncid), &
              'restart_file_to_sv','open '//trim(clm_file))

model_time = get_state_time_ncid(ncid)

!CPS model_time = parflow_time, THIS IS DUMMY FROM PFIDB_DZ FILE

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
!not needed at the moment
!call error_handler(E_ERR,'dart_vector_to_model_file','routine not tested',source, revision,revdate)

! code goes here

end subroutine dart_vector_to_model_file


!------------------------------------------------------------------
!> Reads the current time and state variables from a model data
!> file and packs them into a dart state vector.


function get_parflow_filename()

character(len=256) :: get_parflow_filename

if ( .not. module_initialized ) call static_init_model

get_parflow_filename = trim(parflow_press_file)

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
  if (debug > 1 .and. do_output()) write(*,*) filename
  open(nudat, file=trim(filename),status='old')
 
  read(nudat,*,iostat=izerr) yyyy, mm, dd, ts 
  if (izerr < 0) call error_handler(E_ERR,'pfidb_read','error time', source, revision, revdate)
  if (debug > 10 .and. do_output()) write(*,*) ' DUMMY ...parflow time ',yyyy, mm, dd, ts 

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
    if (debug > 10 .and. do_output()) write(*,'(A,1X,I2,1X,F5.2,1X,F6.3)')   '... pfidb ' ,&
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
  if (nx > 9999 ) then
     errmsg   = "problem readng NX (NX>9999), check pfb file"
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

!> Reads pfb binary file, based on tr32-z4-tools work of P. Shrestha
!>
subroutine pfread_var(filename,pfvar)

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

!> Writes pfb binary file, based on tr32-z4-tools work of P. Shrestha
!>
subroutine pfwrite_var(filename,nx,ny,nz,dx,dy,dz,xd,yd,zd,nxs,nys,pfvar)

character(len=*), intent(in)   :: filename
integer(kind=4),  intent(in)   ::  &
                   nx,             &     ! Longitude dimension
                   ny,             &     ! Latitude dimension
                   nz,             &                    ! Vertical dimension
                   nxs, nys
real(r8),         intent(in)   ::  &
                   dx,             &     ! Grid Resolution in m
                   dy,             &     ! Grid Resoultion in m
                   dz,             &    ! Grid Resolution scale, check parflow namelist
                   xd,yd,zd
real(r8),         intent(in)   ::  &
                   pfvar(nx,ny,nz)       ! ParFlow pressure files 
integer(kind=4)                :: nudat, izerr
character(len=256)             :: errmsg
integer(kind=4)                :: i,j,k, ix, iy, iz, is, ns,       &
                                  rx, ry, rz, nnx, nny, nnz,       &
                                  ixs, iys

! code starts here
  nudat   = get_unit()
  open(nudat,file=trim(filename),form='unformatted',access='stream' , &
                    convert='BIG_ENDIAN',status='new')         ! gfortran

  !read in header infor
  write(nudat, iostat=izerr) xd !X
  if (izerr /= 0) then
    errmsg   = "unable to write  X"
    call error_handler(E_ERR,'pfwrite_var',errmsg,source,revision,revdate)
  endif
  write(nudat, iostat=izerr) yd !Y
  if (izerr /= 0) then
    errmsg   = "unable to write Y"
    call error_handler(E_ERR,'pfwrite_var',errmsg,source,revision,revdate)
  endif
  write(nudat, iostat=izerr) zd !Z
  if (izerr /= 0) then
    errmsg   = "unable to write  Z"
    call error_handler(E_ERR,'pfwrite_var',errmsg,source,revision,revdate)
  endif

  write(nudat, iostat=izerr) nx !NX
  if (izerr /= 0) then
    errmsg   = "unable to write NX"
    call error_handler(E_ERR,'pfwrite_var',errmsg,source,revision,revdate)
  endif
  write(nudat, iostat=izerr) ny !NY
  if (izerr /= 0) then
    errmsg   = "unable to write NY"
    call error_handler(E_ERR,'pfwrite_var',errmsg,source,revision,revdate)
  endif
  write(nudat, iostat=izerr) nz !NZ
  if (izerr /= 0) then
    errmsg   = "unable to write NZ"
    call error_handler(E_ERR,'pfwrite_var',errmsg,source,revision,revdate)
  endif

  write(nudat, iostat=izerr) dx !dX
  if (izerr /= 0) then
    errmsg   = "unable to write  dx"
    call error_handler(E_ERR,'pfwrite_var',errmsg,source,revision,revdate)
  endif
  write(nudat, iostat=izerr) dy !dY
  if (izerr /= 0) then
    errmsg   = "unable to write dy"
    call error_handler(E_ERR,'pfwrite_var',errmsg,source,revision,revdate)
  endif
  write(nudat, iostat=izerr) dz !dZ
  if (izerr /= 0) then
    errmsg   = "unable to write  dz"
    call error_handler(E_ERR,'pfwrite_var',errmsg,source,revision,revdate)
  endif

  ns = int(nxs*nys)
  write(nudat,iostat=izerr)  ns !num_subgrids
  if (izerr /= 0) then
    errmsg   = "unable to write  ns"
    call error_handler(E_ERR,'pfwrite_var',errmsg,source,revision,revdate)
  endif
! End: Writing of domain spatial information

! Start: loop over number of sub grids
  nnx = int(nx/nxs)
  nny = int(ny/nys)
  nnz = nz
  iz = 0
!
  do iys = 0, nys-1
  do ixs = 0, nxs-1

! Start: Writing of sub-grid spatial information

    ix = int(nnx*ixs)
    iy = int(nny*iys)

    write(nudat,iostat=izerr) ix
    if (izerr /= 0) then
      errmsg   = "unable to write  ix"
      call error_handler(E_ERR,'pfwrite_var',errmsg,source,revision,revdate)
    endif

    write(nudat,iostat=izerr) iy
    if (izerr /= 0) then
      errmsg   = "unable to write  iy"
      call error_handler(E_ERR,'pfwrite_var',errmsg,source,revision,revdate)
    endif

    write(nudat,iostat=izerr) iz
    if (izerr /= 0) then
      errmsg   = "unable to write  iz"
      call error_handler(E_ERR,'pfwrite_var',errmsg,source,revision,revdate)
    endif

    write(nudat,iostat=izerr) nnx
    if (izerr /= 0) then
      errmsg   = "unable to write  nnx"
      call error_handler(E_ERR,'pfwrite_var',errmsg,source,revision,revdate)
    endif

    write(nudat,iostat=izerr) nny
    if (izerr /= 0) then
      errmsg   = "unable to write  nny"
      call error_handler(E_ERR,'pfwrite_var',errmsg,source,revision,revdate)
    endif

    write(nudat,iostat=izerr) nnz
    if (izerr /= 0) then
      errmsg   = "unable to write  nnz"
      call error_handler(E_ERR,'pfwrite_var',errmsg,source,revision,revdate)
    endif

    rx = 0; ry = 0; rz=0
    write(nudat,iostat=izerr) rx
    if (izerr /= 0) then
      errmsg   = "unable to write  rx"
      call error_handler(E_ERR,'pfwrite_var',errmsg,source,revision,revdate)
    endif

    write(nudat,iostat=izerr) ry
    if (izerr /= 0) then
      errmsg   = "unable to write  ry"
      call error_handler(E_ERR,'pfwrite_var',errmsg,source,revision,revdate)
    endif

    write(nudat,iostat=izerr) rz
    if (izerr /= 0) then
      errmsg   = "unable to write  rz"
      call error_handler(E_ERR,'pfwrite_var',errmsg,source,revision,revdate)
    endif

! End: Writing of sub-grid spatial information

    ! Start: Write in data from each individual subgrid
    do  k=iz +1 , iz + nnz
    do  j=iy +1 , iy + nny
    do  i=ix +1 , ix + nnx
      write(nudat,iostat=izerr) pfvar(i,j,k)
      if (izerr /= 0) then
        errmsg   = "unable to write  pfvar"
        call error_handler(E_ERR,'pfwrite_var',errmsg,source,revision,revdate)
      endif
    end do
    end do
    end do 

! End: Write in data from each individual subgrid
  end do 
  end do 
! End: loop over number of sub grids

  close(nudat)
  return
end subroutine pfwrite_var

!------------------------------------------------------------------------
!>

!>  define_var_dims() takes the N-dimensional variable and appends the DART
!>  dimensions of 'copy' and 'time'. If the variable initially had a 'time'
!>  dimension, it is ignored because (by construction) it is a singleton
!>  dimension.

subroutine define_var_dims(ivar, ncid, memberdimid, unlimiteddimid, ndims, dimids)

integer,               intent(in)  :: ivar, ncid, memberdimid, unlimiteddimid
integer,               intent(out) :: ndims
integer, dimension(:), intent(out) :: dimids

character(len=NF90_MAX_NAME),dimension(NF90_MAX_VAR_DIMS) :: dimnames

integer :: i, io, mydimid

ndims = 0

DIMLOOP : do i = 1,progvar(ivar)%numdims

   if (progvar(ivar)%dimnames(i) == 'time') cycle DIMLOOP

   io = nf90_inq_dimid(ncid=ncid, name=progvar(ivar)%dimnames(i), dimid=mydimid)
   call nc_check(io, 'define_var_dims','inq_dimid '//trim(progvar(ivar)%dimnames(i)))

   ndims         = ndims + 1
   dimids(ndims) = mydimid
   dimnames(ndims) = progvar(ivar)%dimnames(i)

enddo DIMLOOP

! The last two dimensions are always 'copy' and 'time'
ndims           = ndims + 1
dimids(ndims)   = memberdimid
dimnames(ndims) = 'copy'
ndims           = ndims + 1
dimids(ndims)   = unlimitedDimid
dimnames(ndims) = 'time'

return
end subroutine define_var_dims

!------------------------------------------------------------------------
!> Horizontal interpolation at two height indices

subroutine horizontal_interpolate(x, ivar, kk_inds, ij_inds, ij_wgts, interp_hv, istatus)

real(r8), intent(in)  :: x(:)
integer,  intent(in)  :: ivar
integer,  intent(in)  :: kk_inds(2)
integer,  intent(in)  :: ij_inds(4)     ! ileft, iright, jbot, jtop
real(r8), intent(in)  :: ij_wgts(2)     ! ifrac, jfrac
real(r8), intent(out) :: interp_hv(2)  ! at two height levels 
integer,  intent(out) :: istatus

integer  :: k, ll, lr, ur, ul
real(r8) :: ij_rwgts(2)

! Location of dart_indices
ll = ijk_to_dart(ivar, ij_inds(1), ij_inds(3), kk_inds(1))
lr = ijk_to_dart(ivar, ij_inds(2), ij_inds(3), kk_inds(1))
ur = ijk_to_dart(ivar, ij_inds(2), ij_inds(4), kk_inds(1))
ul = ijk_to_dart(ivar, ij_inds(1), ij_inds(4), kk_inds(1))

! Interpolate
ij_rwgts(1) = 1.0_r8 - ij_wgts(1)  !ifrac
ij_rwgts(2) = 1.0_r8 - ij_wgts(2)  !jfrac

do k = 1, 2 

 interp_hv(k) = ij_rwgts(2) * ( ij_rwgts(1) * x(ll) + ij_wgts(1) * x(lr)) + &
                 ij_wgts(2)  * ( ij_rwgts(1) * x(ul) + ij_wgts(1) * x(ur))

end do

if (debug > 99 .and. do_output()) then
   write(*,*)
   write(*,*)' ileft, jbot, level_index decompose to ', ll, x( ll), ij_rwgts(1)  
   write(*,*)'iright, jbot, level_index decompose to ',lr, x(lr), ij_wgts(1)
   write(*,*)'iright, jtop, level_index decompose to ',ur, x(ur), ij_rwgts(2)  
   write(*,*)' ileft, jtop, level_index decompose to ', ul, x( ul),ij_rwgts(2) 
   write(*,*)' horizontal value is ', kk_inds(1), interp_hv(1)
   write(*,*)' horizontal value is ', kk_inds(2), interp_hv(2)
endif

return

end subroutine horizontal_interpolate

!------------------------------------------------------------------------
!>  Adapted from COSMO

function ijk_to_dart(ivar, i, j, k) result(dartindex)

integer, intent(in) :: ivar
integer, intent(in) :: i
integer, intent(in) :: j
integer, intent(in) :: k
integer             :: dartindex

if (progvar(ivar)%numdims == 3) then
   dartindex = progvar(ivar)%index1 + (i-1) + &
               (j-1) * (progvar(ivar)%dimlens(1)) + &
               (k-1) * (progvar(ivar)%dimlens(1) * progvar(ivar)%dimlens(2))

elseif (progvar(ivar)%numdims == 2) then
   dartindex = progvar(ivar)%index1 + (i-1) + &
               (j-1) * (progvar(ivar)%dimlens(1))

elseif (progvar(ivar)%numdims == 1) then
   dartindex = progvar(ivar)%index1 + (i-1)

endif

end function

!------------------------------------------------------------------------
!> Get model corners for the observed location
!> Adapted from cosmo interp routines

subroutine get_corners(lon, lat, nx, ny, gridlons, gridlats, ij_inds, &
                        ij_wgts, istatus)

real(r8), intent(in)  :: lon
real(r8), intent(in)  :: lat
integer,  intent(in)  :: nx
integer,  intent(in)  :: ny
real(r8), intent(in)  :: gridlons(nx)
real(r8), intent(in)  :: gridlats(ny)

integer,  intent(out) :: ij_inds(4) ! ileft, iright, jbot, jtop
real(r8), intent(out) :: ij_wgts(2)    ! ifrac, jfrac

!Local Variables
integer               :: istatus, indarr(1)
integer, parameter    :: OUTSIDE_HORIZONTALLY = 15
real(r8)              :: dlon, dlonT, dlat, dlatT

! Initialize
ij_inds    = -1
ij_wgts    = 0.0_r8
istatus    = 99

! check to make sure the location is within the domain

if ( lon < gridlons(1) .or.  lon > gridlons(nx) .or. &
     lat < gridlats(1) .or.  lat > gridlats(ny) )then
   istatus = OUTSIDE_HORIZONTALLY
   return
endif

!    | ............ dlonT ............. |
!    | .. dlon .. |
!    |------------|---------------------|
! lonleft       lon                  lonright
! ij_inds(1)                         ij_inds(2)

if (lon == gridlons(nx)) then
   ij_inds(1)  = nx-1
   ij_inds(2)  = nx
   ij_wgts(1)  = 0.0_r8  ! fractional distance to min
else
   indarr = maxloc(gridlons , gridlons <= lon)
   ij_inds(1)  = indarr(1)
   ij_inds(2)  = ij_inds(1) + 1
   dlon        =     lon          - gridlons(ij_inds(1))
   dlonT       = gridlons(ij_inds(2)) - gridlons(ij_inds(1))
   ij_wgts(1)  = dlon / dlonT
endif

! latitudes are arranged 'south' to 'north', i.e. -90 to 90  (albeit rotated)
!
!    | ............ dlatT ............. |
!    | .. dlat .. |
!    |------------|---------------------|
! latbot         lat                  lattop
! ij_inds(3)                         ij_inds(4)

if (lat == gridlats(ny)) then
   ij_inds(3) = ny-1
   ij_inds(4) = ny
   ij_wgts(2) = 0.0_r8
else
   indarr = maxloc( gridlats,  gridlats <= lat)
   ij_inds(3)   = indarr(1)
   ij_inds(4)   = ij_inds(3) + 1
   dlat         =     lat        - gridlats(ij_inds(3))
   dlatT        = gridlats(ij_inds(4)) - gridlats(ij_inds(3))
   ij_wgts(2)   = dlat / dlatT
endif

istatus = 0

if (debug > 5 .and. do_output()) then
   write(*,*)
   write(*,*)'longitude index to the  west and  east are ',ij_inds(1:2)
   write(*,*)'latitude  index to the south and north are ',ij_inds(3:4) 
   write(*,*)'lon  west, lon, lon  east',gridlons(ij_inds(1)),lon,gridlons(ij_inds(2))
   write(*,*)'lat south, lat, lat north',gridlats(ij_inds(3)), lat,gridlats(ij_inds(4))
endif

end subroutine get_corners


!------------------------------------------------------------------------
!> Get bottom and top model levels  for the observed location
!> Adapted from cosmo interp routines

subroutine get_level_indices(height, kk_inds, kk_wgt, istatus)

! The height for vcoord array is from bottom to top
! vcoord (1) = bottom of the subsurface
! vcoord (nz) = top of the subsurface

!         bottom <--------------------------------------------| surface
! layer     1 ... 10                                 11  ... nz 
!                 | ............ dztot ............. |
!                 | ... dz ... |
!                 |------------|---------------------|
!               ztop         height                 zbot
!              iabove                              ibelow
!              kk_inds(1)                          kk_inds(2) 


real(r8), intent(in)  :: height
integer,  intent(out) :: kk_inds(2) 
real(r8), intent(out) :: kk_wgt 
integer,  intent(out) :: istatus

integer, parameter :: OUTSIDE_VERTICALLY   = 16

integer :: indarr(1), i
real(r8) :: dz, dztot

! Variable definition complete
if (height > vcoord(1) .or. height < vcoord(nz)) then
  istatus = OUTSIDE_VERTICALLY
  return
endif

if (height == vcoord(1)) then
  kk_inds(1) = 1
  kk_inds(2) = 2
  kk_wgt     = 0._r8
else
  indarr  = minloc( vcoord,  vcoord >= height )
  !write(*,*) "CPS",vcoord
  !write(*,*) "CPS", height, indarr(1), vcoord(indarr(1))
  kk_inds(1) = indarr(1)
  kk_inds(2) = kk_inds(1) + 1
  dz     = vcoord(kk_inds(1)) - height
  dztot  = vcoord(kk_inds(1)) - vcoord(kk_inds(2))
  kk_wgt = dz / dztot
endif

if (debug > 99 .and. do_output()) then
  write(*,*)
  write(*,*)'vertical levels for ParFlow '
  write(*,*)'iabove, levelfrac, ibelow ', kk_inds(1), kk_wgt, kk_inds(2)
  write(*,*)'iabove, height, ibelow ', vcoord(kk_inds(1)), height, vcoord(kk_inds(2))
  write(*,*)
endif

istatus = 0

end subroutine get_level_indices

!------------------------------------------------------------------
function get_state_time_ncid( ncid )
!------------------------------------------------------------------
! The restart netcdf files have the time of the state.

type(time_type) :: get_state_time_ncid
integer, intent(in) :: ncid 

integer :: VarID
integer :: rst_curr_ymd, rst_curr_tod, leftover
integer :: year, month, day, hour, minute, second

if ( .not. module_initialized ) call static_init_model

call nc_check(nf90_inq_varid(ncid, 'timemgr_rst_curr_ymd', VarID),'get_state_time_ncid', &
&  'inq_varid timemgr_rst_curr_ymd'//trim(clm_file))
call nc_check(nf90_get_var(  ncid, VarID,   rst_curr_ymd),'get_state_time_ncid', &
&            'get_var rst_curr_ymd'//trim(clm_file))

call nc_check(nf90_inq_varid(ncid, 'timemgr_rst_curr_tod', VarID),'get_state_time_ncid', &
&  'inq_varid timemgr_rst_curr_tod'//trim(clm_file))
call nc_check(nf90_get_var(  ncid, VarID,   rst_curr_tod),'get_state_time_ncid', &
&            'get_var rst_curr_tod'//trim(clm_file))

year     = rst_curr_ymd/10000
leftover = rst_curr_ymd - year*10000
month    = leftover/100
day      = leftover - month*100

hour     = rst_curr_tod/3600
leftover = rst_curr_tod - hour*3600
minute   = leftover/60
second   = leftover - minute*60

get_state_time_ncid = set_date(year, month, day, hour, minute, second)

end function get_state_time_ncid

!------------------------------------------------------------------


!===================================================================
! End of model_mod
!===================================================================
end module model_mod

! <next few lines under version control, do not edit>
! $URL: https://proxy.subversion.ucar.edu/DAReS/DART/releases/Lanai/models/template/model_mod.f90 $
! $Id: model_mod.f90 6256 2013-06-12 16:19:10Z thoar $
! $Revision: 6256 $
! $Date: 2013-06-12 18:19:10 +0200 (Wed, 12 Jun 2013) $
