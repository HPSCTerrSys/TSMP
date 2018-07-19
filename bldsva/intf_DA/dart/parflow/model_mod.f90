! This code is not necessarily under the DART copyright ...
!
! DART $Id: model_mod.f90 Wed Apr 11 20:26:43 CEST 2018 $

module model_mod

! This module contains routines to work with ParFlow Binary data 
! in the DART framework

! Based on DART template directory

! Modules that are absolutely required for use are listed
use        types_mod, only : r8, obstypelength, MISSING_R8

use time_manager_mod, only : time_type, set_time, get_time, print_time, &
                                        set_date, get_date, print_date, &
                             set_calendar_type, set_time_missing,       &
                             operator(*),  operator(+), operator(-),    &    
                             operator(>),  operator(<), operator(/),    &    
                             operator(/=), operator(<=), operator(==)

use     location_mod, only : location_type,      get_close_maxdist_init,        &
                             get_close_obs_init, get_close_obs, set_location,   &
                             get_location, vert_is_height, VERTISHEIGHT

use    utilities_mod, only : register_module, error_handler, nc_check, &
                             get_unit, open_file, close_file, E_ERR, E_MSG, &
                             nmlfileunit, do_output, do_nml_file, do_nml_term,  &
                             find_namelist_in_file, file_exist, to_upper, &
                             check_namelist_read, logfileunit

use     obs_kind_mod, only : paramname_length, &
                             get_raw_obs_kind_index, &
                             get_raw_obs_kind_name, &
                             KIND_3D_PARAMETER, &
                             KIND_SOIL_SATURATION, &
                             KIND_SOIL_WATER_CONTENT, &
                             KIND_SOIL_MOISTURE

use netcdf

implicit none
private

! version controlled file description for error handling, do not edit
character(len=*), parameter :: source   = "$URL: model_mod.f90 $"
character(len=*), parameter :: revision = "$Revision: Bonn $"
character(len=*), parameter :: revdate  = "$Date: Wed Apr 11 2018 $"

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
public :: get_state_vector,     &
          get_parflow_filename, &
          get_parflow_id,       &
          write_parflow_file,   &
          PRESSURE_HEAD, SATURATION


character(len=512) :: string1, string2, string3
logical, save :: module_initialized = .false.

! Dimension and grid resolution from (parflow binary) PFB file
integer(kind=4)                ::  &
                   nx,             &     ! Longitude dimension
                   ny,             &     ! Latitude dimension
                   nz                    ! Vertical dimension
real(r8)                       ::  &
                   dx,             &     ! Grid Resolution in m
                   dy,             &     ! Grid Resoultion in m
                   dz                    ! Grid Resolution scale, check parflow namelist
real(r8)                       :: xdpfb, ydpfb, zdpfb            ! ParFlow PFB file

integer          :: model_size
type(time_type)  :: time_step
type(time_type)  :: start_date      !Needed for Parflow nos

type(location_type), allocatable :: state_loc(:)

! Codes for restricting the range of a variable
integer, parameter :: BOUNDED_NONE  = 0 ! ... unlimited range
integer, parameter :: BOUNDED_BELOW = 1 ! ... minimum, but no maximum
integer, parameter :: BOUNDED_ABOVE = 2 ! ... maximum, but no minimum
integer, parameter :: BOUNDED_BOTH  = 3 ! ... minimum and maximum

integer, parameter :: PRESSURE_HEAD = 1
integer, parameter :: SATURATION = 2

! Everything needed to describe a variable
integer, parameter :: max_state_variables = 2   !ParFlow has press and satur 

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
integer,  allocatable ::   soil_type(:,:,:)     ! pfl soil indicator 
real(r8), allocatable ::   soil_parameters(:,:) ! soil attributes (Ss,Sr,a,N,Phi)
type(time_type)       ::   parflow_time  ! parflow ouptut time
integer               ::   ipfb ! determine from  parflow_assim_variable

! run-time options - the namelist and the defaults
character(len=256) :: parflow_press_file           = 'parflow_press_file'
character(len=256) :: parflow_satur_file           = 'parflow_satur_file'
character(len=256) :: grid_file                    = 'grid_file'
character(len=256) :: clm_file                     = 'clm_restart.nc'
character(len=256) :: clm_file_s                   = 'clm_restart_s.nc'
character(len=256) :: soilInd_file                 = 'pfl_soil_file'
character(len=256) :: parflow_assim_variable       = 'saturation'
integer            :: assimilation_period_days     = 0
integer            :: assimilation_period_seconds  = 21400
integer            :: debug = 0 


namelist /model_nml/                &
      parflow_press_file,           &
      parflow_satur_file,           &
      grid_file,                    &
      clm_file,                     &
      clm_file_s,                   &
      soilInd_file,                 &
      parflow_assim_variable,       &
      assimilation_period_days,     &
      assimilation_period_seconds,  &
      debug

contains


!==================================================================


!------------------------------------------------------------------
!> Called to do one time initialization of the model. As examples,
!> might define information about the model size or model timestep.
!> In models that require pre-computed static data, for instance
!> spherical harmonic weights, these would also be computed here.
!> Can be a NULL INTERFACE for the simplest models.

subroutine static_init_model()

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
start_date = set_time_missing()

! Determine which variable will be part of the DART state
! PRESSURE_HEAD is the prognostic variable in parflow, so it will always be updated
! no matter what is in the DART state.
write(string3,*)'read "'//trim(parflow_assim_variable)//'"'
call to_upper(parflow_assim_variable)
select case (parflow_assim_variable)
   case ('SATURATION')
      ipfb = SATURATION
   case ('PRESSURE_HEAD')
      ipfb = PRESSURE_HEAD
   case default
      write(string1,*)'unsupported value for "model_nml:parflow_assim_variable"' 
      write(string2,*)'can be either (case-insensitive) "pressure_head" or "saturation"'
      call error_handler(E_ERR, 'static_init_model', string1, &
                 source, revision, revdate, text2=string2,text3=string3)
end select   

! Get the dimensions of the grid and the grid variables from the netCDF file.
call pfread_dim(parflow_press_file)

if (debug > 0 .and. do_output()) write(*,*) '   ... pfb      dimensions are ' , nx, ny, nz
if (debug > 0 .and. do_output()) write(*,'(A,3(1X,F9.2))') '    ... pfb grid resolution are ', dx, dy, dz

model_size = nx*ny*nz

! Get the geo-locaion
allocate(lon(nx,ny))
allocate(lat(nx,ny))

call grid_read(grid_file)

if (debug > 1 .and. do_output()) then
   write(*, '(A,2(1X,F20.14))') '    ... parflow lon range is ', minval(lon), maxval(lon)
   write(*, '(A,2(1X,F20.14))') '    ... parflow lat range is ', minval(lat), maxval(lat)
   write(*, *) '-------------------------------'
end if

! Get the vertical co-ordinate and soil indicators
call read_soil_table(soilInd_file)

! TODO, var_type could switch between 1 and 2 for pressure and sat

ivar                          = PRESSURE_HEAD
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
progvar(ivar)%pfb_kind        = KIND_3D_PARAMETER   !Need pedotransfer function 
progvar(ivar)%index1          = 1
progvar(ivar)%rangeRestricted = BOUNDED_NONE        
progvar(ivar)%minvalue        = MISSING_R8 
progvar(ivar)%maxvalue        = MISSING_R8

ivar                          = SATURATION
progvar(ivar)%varname         = "sat"
progvar(ivar)%long_name       = "Relative Saturation"
progvar(ivar)%units           = "-"
progvar(ivar)%numdims         = 3  
progvar(ivar)%dimnames(1)     = 'pfl_lon'
progvar(ivar)%dimnames(2)     = 'pfl_lat'
progvar(ivar)%dimnames(3)     = 'level'
progvar(ivar)%dimlens(1)      = nx 
progvar(ivar)%dimlens(2)      = ny 
progvar(ivar)%dimlens(3)      = nz 
progvar(ivar)%pfb_kind        = KIND_SOIL_SATURATION
progvar(ivar)%index1          = 1  
progvar(ivar)%rangeRestricted = BOUNDED_BOTH      
progvar(ivar)%minvalue        = 0._r8                != f(soil index) 
progvar(ivar)%maxvalue        = 1._r8 

! Create storage for locations
allocate(state_loc(model_size))

! The time_step in terms of a time type must also be initialized.
time_step = set_time(assimilation_period_seconds, assimilation_period_days)

end subroutine static_init_model




!------------------------------------------------------------------
!> Returns a model state vector, x, that is some sort of appropriate
!> initial condition for starting up a long integration of the model.
!> At present, this is only used if the namelist parameter 
!> start_from_restart is set to .false. in the program perfect_model_obs.
!> If this option is not to be used in perfect_model_obs, or if no 
!> synthetic data experiments using perfect_model_obs are planned, 
!> this can be a NULL INTERFACE.

subroutine init_conditions(x)

real(r8), intent(out) :: x(:)

if ( .not. module_initialized ) call static_init_model

write(string1,*) 'input.nml:start_from_restart cannot be FALSE'
call error_handler(E_ERR,'init_conditions',string1,source, revision,revdate)

x = MISSING_R8

return
end subroutine init_conditions



!------------------------------------------------------------------
!> Does a single timestep advance of the model. The input value of
!> the vector x is the starting condition and x is updated to reflect
!> the changed state after a timestep. The time argument is intent
!> in and is used for models that need to know the date/time to 
!> compute a timestep, for instance for radiation computations.
!> This interface is only called if the namelist parameter
!> async is set to 0 in perfect_model_obs of filter or if the 
!> program integrate_model is to be used to advance the model
!> state as a separate executable. If one of these options
!> is not going to be used (the model will only be advanced as
!> a separate model-specific executable), this can be a 
!> NULL INTERFACE.

subroutine adv_1step(x, time)

real(r8),        intent(inout) :: x(:)
type(time_type), intent(in)    :: time

if ( .not. module_initialized ) call static_init_model

write(string1,*) 'Cannot advance ParFlow with a subroutine call; async cannot equal 0'
call error_handler(E_ERR,'adv_1step',string1,source,revision,revdate)

end subroutine adv_1step


!------------------------------------------------------------------
!> writes the time details to a text file for parflow

subroutine write_state_times(iunit, statetime)

integer,                   intent(in) :: iunit
type(time_type),           intent(in) :: statetime

integer           :: iyear, imonth, iday, ihour, imin, isec 
integer           :: ndays, nhours, nmins, nsecs,nfreq
type(time_type)   :: interval
type(time_type)   :: missing

missing = set_time_missing()

if (start_date == missing) then
   write(string1,*)'TJH start_date was not set'
   call error_handler(E_ERR,'write_state_time',string1,source,revision,revdate)
endif

call get_date(statetime, iyear, imonth, iday, ihour, imin, isec)
nsecs = (ihour*60 + imin)*60 + isec 
write(iunit, '(''clmext '',I4.4,2(''-'',I2.2),''-'',i5.5)') iyear, imonth, iday, nsecs
write(iunit, '(''defaultInitDate '',I4.4,2(''-'',I2.2),1x,i2.2)') iyear, imonth, iday, ihour 

interval = statetime - start_date

call get_time(interval, nsecs, ndays)

nhours = nsecs / (60*60)
nsecs  = nsecs - (nhours * 60*60)
nmins  = nsecs / 60 
nsecs  = nsecs - (nmins * 60)
nfreq  = 3   !PFL SPECIFIC

write(iunit, '(''cosrbin '',''lrff'',4(I2.2),''o'')') ndays, nhours, nmins, nsecs
write(iunit, '(''coshist '',''lfff'',4(I2.2),''.nc'')') ndays, nhours, nmins, nsecs
write(iunit, '(''pflhist '',I5.5)') (ndays*24 + nhours)/nfreq

call get_date(start_date, iyear, imonth, iday, ihour, imin, isec)
write(iunit, '(''defaultStartDate '',I4.4,2(''-'',I2.2),1x,i2.2)') iyear,imonth, iday, ihour 

return

end subroutine write_state_times


!------------------------------------------------------------------
!> Returns the size of the model as an integer. Required for all
!> applications.

function get_model_size()

integer :: get_model_size

if ( .not. module_initialized ) call static_init_model

get_model_size = model_size

end function get_model_size


!------------------------------------------------------------------
!> Companion interface to init_conditions. Returns a time that is somehow 
!> appropriate for starting up a long integration of the model.
!> At present, this is only used if the namelist parameter 
!> start_from_restart is set to .false. in the program perfect_model_obs.
!> If this option is not to be used in perfect_model_obs, or if no 
!> synthetic data experiments using perfect_model_obs are planned, 
!> this can be a NULL INTERFACE.

subroutine init_time(time)

type(time_type), intent(out) :: time

if ( .not. module_initialized ) call static_init_model

write(string1,*)'input.nml:start_from_restart cannot be FALSE'
call error_handler(E_ERR,'init_time',string1,source, revision,revdate)
! for now, just set to 0
time = set_time(0,0)

end subroutine init_time


!------------------------------------------------------------------------
!> Given a state vector, a location, and a model state variable type,
!> interpolates the state variable field to that location and returns
!> the value in obs_val. The istatus variable should be returned as
!> 0 unless there is some problem in computing the interpolation in
!> which case an alternate value should be returned. The obs_type variable
!> is a model specific integer that specifies the type of field .
 
subroutine model_interpolate(x, location, obs_type, interp_val, istatus)

! Error codes:
! istatus = 99 : unknown error
! istatus = 10 : observation type is not in state vector
! istatus = 11 : observation type cannot be computed
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

! Identity observation support

if ( obs_type < 0 ) then
   interp_val = x(-1*obs_type)
   istatus = 0
   return
endif

! Get the geo-location
point_coords(1:3) = get_location(location)
geo_lon = point_coords(1) ! degrees East
geo_lat = point_coords(2) ! degrees North
geo_hgt = point_coords(3) ! depth in meters 

! Determine if the state has soil saturation

ivar = -1
FoundIt: do i = 1, max_state_variables
   if (progvar(i)%pfb_kind == KIND_SOIL_SATURATION) then
      ivar = i
      exit FoundIt
   endif
enddo FoundIt
if (ivar < 0) then
   istatus = NO_KIND_IN_STATE
   return
endif

! Get the corners, weights, etc. needed by every method
call get_corners(geo_lon, geo_lat, nx, ny, lon(1:nx,1), lat(1,1:ny), &
                 geo_inds, geo_wgts, istatus) 
if (istatus /= 0) return

call get_level_indices(geo_hgt, hgt_inds, hgt_wgt, istatus)
if (istatus /= 0) return

if (obs_type == KIND_SOIL_SATURATION) then

   ! Horizontal interpolation at two levels given by hgt_inds
   call horizontal_interpolate(x, ivar, hgt_inds, geo_inds, geo_wgts, &
                            interp_h_var, istatus)
   if (istatus /= 0) return

   ! Vertical interpolation to scalar estimate.
   interp_val = interp_h_var(2)*hgt_wgt + interp_h_var(1)*(1.0_r8 - hgt_wgt)

else if (obs_type == KIND_SOIL_MOISTURE) then  ! cm^3/cm^3

   call calculate_soil_moisture(x, ivar, hgt_inds, geo_inds, geo_wgts, interp_h_var, istatus)
   if (istatus /= 0) return
   interp_val   = interp_h_var(2)*hgt_wgt + interp_h_var(1)*(1.0_r8 - hgt_wgt)

elseif (obs_type == KIND_SOIL_WATER_CONTENT) then

   call calculate_soil_moisture(x, ivar, hgt_inds, geo_inds, geo_wgts, interp_h_var, istatus)
   if (istatus /= 0) return
   interp_val = interp_h_var(2)*hgt_wgt + interp_h_var(1)*(1.0_r8 - hgt_wgt)

   interp_val = interp_val * 100.0_r8

endif

end subroutine model_interpolate


!------------------------------------------------------------------
!> Returns the the time step of the model; the smallest increment
!> in time that the model is capable of advancing the state in a given
!> implementation. This interface is required for all applications.

function get_model_time_step()

type(time_type) :: get_model_time_step

if ( .not. module_initialized ) call static_init_model

get_model_time_step = time_step

return

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

local_ind = index_in - 1    !CPS offset for algorithm below

kloc =  local_ind / (nx*ny) + 1
jloc = (local_ind - (kloc-1)*nx*ny)/nx + 1
iloc =  local_ind - (kloc-1)*nx*ny - (jloc-1)*nx + 1

if (debug > 99 .and. do_output()) then
  write(*,*)'.. index_in ',index_in, ' and dereferences to i,j,k ',iloc,jloc,kloc
endif

! Now that we know the i,j,k we have to index the right set of
! coordinate arrays
mylon = lon(iloc,jloc)
mylat = lat(iloc,jloc)
vloc  = vcoord(kloc)

location = set_location(mylon, mylat, vloc, VERTISHEIGHT) ! meters

if (present(var_type)) var_type = progvar(ipfb)%pfb_kind 

end subroutine get_state_meta_data


!------------------------------------------------------------------
!> Does any shutdown and clean-up needed for model. Can be a NULL
!> INTERFACE if the model has no need to clean up storage, etc.

subroutine end_model()

deallocate(lon,lat)
deallocate(vcoord)
deallocate(state_loc)
deallocate(soil_type)
deallocate(soil_parameters)
!if ( allocated(ens_mean) ) deallocate(ens_mean)

end subroutine end_model


!------------------------------------------------------------------------
!> 
!> Writes the model-specific attributes to a netCDF file.
!> This includes coordinate variables and some metadata, but NOT
!> the actual model state vector. We do have to allocate SPACE for the model
!> state vector, but that variable gets filled as the model advances.

function nc_write_model_atts( ncFileID ) result (ierr)

use typeSizes
use netcdf

integer, intent(in)  :: ncFileID      ! netCDF file identifier
integer              :: ierr          ! return value of function

! All errors are fatal, so the return code is always '0 == normal'.

integer :: nDimensions, nVariables, nAttributes, unlimitedDimID

integer :: StateVarDimID   ! netCDF pointer to state variable dimension (model size)
integer :: MemberDimID     ! netCDF pointer to dimension of ensemble    (ens_size)
integer :: TimeDimID       ! netCDF pointer to time dimension           (unlimited)
integer :: vcoordVarID
integer :: VarID

integer :: lonDimID
integer :: latDimID
integer :: levelDimID

character(len=512) :: errstring
character(len=256) :: filename

! we are going to need these to record the creation date in the netCDF file.

character(len=8)      :: crdate      ! needed by F90 DATE_AND_TIME intrinsic
character(len=10)     :: crtime      ! needed by F90 DATE_AND_TIME intrinsic
character(len=5)      :: crzone      ! needed by F90 DATE_AND_TIME intrinsic
integer, dimension(8) :: values      ! needed by F90 DATE_AND_TIME intrinsic

integer, dimension(NF90_MAX_VAR_DIMS) :: mydimids
character(len=NF90_MAX_NAME) :: varname
character(len=NF90_MAX_NAME) :: str1

integer :: io, ndims, ivar

if ( .not. module_initialized ) call static_init_model

! make sure ncFileID refers to an open netCDF file, 
! and then put into define mode.

ierr = 0 ! all errors are fatal, so this is always 0 ... kinda stupid on my part (TJH)

call nc_check(nf90_inquire(ncFileID,nDimensions,nVariables,nAttributes,unlimitedDimID), &
                     "nc_write_model_atts", "inquire")
call nc_check(nf90_redef(ncFileID), "nc_write_model_atts", "redef")

! We need the dimension ID for the number of copies/ensemble members, and
! we might as well check to make sure that Time is the Unlimited dimension. 
! Our job is create the 'model size' dimension.

call nc_check(nf90_inq_dimid(ncid=ncFileID, name="copy", dimid=MemberDimID), &
                            "nc_write_model_atts", "inq_dimid copy")
call nc_check(nf90_inq_dimid(ncid=ncFileID, name="time", dimid= TimeDimID), &
                            "nc_write_model_atts", "inq_dimid time")

if ( TimeDimID /= unlimitedDimId ) then
   write(errstring,*)"Time Dimension ID ",TimeDimID, &
                     " should equal Unlimited Dimension ID",unlimitedDimID
   call error_handler(E_ERR,"nc_write_model_atts", errstring, source, revision, revdate)
endif

! Define the model size / state variable dimension / whatever ...

call nc_check(nf90_def_dim(ncid=ncFileID, name="StateVariable",  &
                           len=model_size, dimid=StateVarDimID), &
                           "nc_write_model_atts", "def_dim state")

! Write Global Attributes 

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

! Declare space for the variables we want in the diagnostic files.

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

! Create the (empty) Variables and the Attributes

do ivar=ipfb, ipfb

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

call nc_check(nf90_enddef(ncfileID), "nc_write_model_atts", "enddef")

! Fill the coordinate variables - reshape the 1D arrays to the 2D shape

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

! Flush the buffer and leave netCDF file open

call nc_check(nf90_sync(ncFileID),"nc_write_model_atts", "sync")

end function nc_write_model_atts


!------------------------------------------------------------------------
!> Writes the model variables to a netCDF file.
!> All errors are fatal, so the return code is always '0 == normal'.

function nc_write_model_vars( ncFileID, statevec, copyindex, timeindex ) result (ierr)         
use typeSizes
use netcdf

integer,  intent(in) :: ncFileID      ! netCDF file identifier
real(r8), intent(in) :: statevec(:)
integer,  intent(in) :: copyindex
integer,  intent(in) :: timeindex
integer              :: ierr          ! return value of function

character(len=NF90_MAX_NAME) :: varname
integer ::  dimIDs(NF90_MAX_VAR_DIMS)
integer :: ncstart(NF90_MAX_VAR_DIMS)
integer :: nccount(NF90_MAX_VAR_DIMS)

integer :: io, i, ivar, VarID, ndims, dimlen
integer :: TimeDimID, CopyDimID

character(len=256) :: filename

if ( .not. module_initialized ) call static_init_model

ierr = 0 ! all errors are fatal, so this is always 0 ... kinda stupid on my part (TJH)

write(filename,*) 'ncFileID', ncFileID

io = nf90_inq_dimid(ncFileID, 'copy', dimid=CopyDimID)
call nc_check(io, 'nc_write_model_vars', 'inq_dimid copy '//trim(filename))

io = nf90_inq_dimid(ncFileID, 'time', dimid=TimeDimID)
call nc_check(io, 'nc_write_model_vars', 'inq_dimid time '//trim(filename))

do ivar = ipfb, ipfb

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

   if (debug > 99 .and. do_output()) then
      write(*,*)'nc_write_model_vars '//trim(varname)//' start is ',ncstart(1:ndims)
      write(*,*)'nc_write_model_vars '//trim(varname)//' count is ',nccount(1:ndims)
   endif

   ! revelation - no need to reshape the state vector before nf90_put_var

   call nc_check(nf90_put_var(ncFileID, VarID, statevec,          &
                start = ncstart(1:ndims), count=nccount(1:ndims)), &
                'nc_write_model_vars', 'put_var '//trim(string2))
enddo

! Flush the buffer and leave netCDF file open

call nc_check(nf90_sync(ncFileID), "nc_write_model_vars", "sync")

ierr = 0 ! If we got here, things went well.

end function nc_write_model_vars


!------------------------------------------------------------------
!> Perturbs a model state for generating initial ensembles.
!> The perturbed state is returned in pert_state.
!> A model may choose to provide a NULL INTERFACE by returning
!> .false. for the interf_provided argument. This indicates to
!> the filter that if it needs to generate perturbed states, it
!> may do so by adding an O(0.1) magnitude perturbation to each
!> model state variable independently. The interf_provided argument
!> should be returned as .true. if the model wants to do its own
!> perturbing of states.  The returned pert_state should in any
!> case be valid, since it will be read by filter even if 
!> interf_provided is .false.

subroutine pert_model_state(state, pert_state, interf_provided)

real(r8), intent(in)  :: state(:)
real(r8), intent(out) :: pert_state(:)
logical,  intent(out) :: interf_provided

if ( .not. module_initialized ) call static_init_model
call error_handler(E_ERR,'pert_model_state','routine not tested',source, revision,revdate)

pert_state      = state
interf_provided = .false.

end subroutine pert_model_state


!------------------------------------------------------------------
!> Needed for vertical conversion. Since we are not
!> converting from some hybrid vertical coordinates ... not needed.

subroutine ens_mean_for_model(ens_mean)

real(r8), intent(in) :: ens_mean(:)

if ( .not. module_initialized ) call static_init_model

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

subroutine get_state_vector(sv, sv_id, model_time)

real(r8),         intent(inout)           :: sv(1:model_size)
integer,          intent(in)              :: sv_id     !to choose what to get
type(time_type),  intent(out), optional   :: model_time

character(len=*), parameter :: routine = 'get_state_vector'
real(r8),      allocatable  :: pfbdata(:,:,:)

if ( .not. module_initialized ) call static_init_model

allocate(pfbdata(nx,ny,nz))

if (sv_id == PRESSURE_HEAD) then

   call pfread_var(parflow_press_file, pfbdata) 
   write(string1,'(A)')'reading from "'//trim(parflow_press_file)//'"'

elseif (sv_id == SATURATION) then

   call pfread_var(parflow_satur_file, pfbdata)
   write(string1,'(A)')'reading from "'//trim(parflow_satur_file)//'"'

else

   write(string1,*)'sv_id should be 1 or 2'
   write(string2,*)'was ',sv_id
   call error_handler(E_ERR, routine, string1, &
              source, revision, revdate, text2=string2)  
endif

if (debug > 0 .and. do_output()) call error_handler(E_MSG,routine,string1)

sv(:) = reshape(pfbdata,(/ (nx*ny*nz) /))

!>@todo What is the logic here ...
!> if the sv_id to read is the one we are using to assimilate ... 
!> we define the start_date, if not -   !>@todo check

! if (sv_id == ipfb) then  !>@todo check

   !>@todo should these two calls be in static_init_model()?
   if (present(model_time)) model_time = get_model_time(clm_file)

   start_date = get_start_time(clm_file_s)

   ! HAVE TO USE THE START FILE TO GET THE CORRECT START DATE :(((
   !CPS model_time = parflow_time, THIS IS DUMMY FROM PFIDB_DZ FILE

   if (debug > 0 .and. do_output()) then
      if (present(model_time)) then
      call print_date( model_time,' get_state_vector model date')
      call print_time( model_time,' get_state_vector model time')
      call print_date( model_time,' get_state_vector model date',logfileunit)
      call print_time( model_time,' get_state_vector model time',logfileunit)
      endif

      call print_date( start_date,' get_state_vector start date')
      call print_time( start_date,' get_state_vector start time')
      call print_date( start_date,' get_state_vector start date',logfileunit)
      call print_time( start_date,' get_state_vector start time',logfileunit)
   endif

! endif  !>@todo check

deallocate(pfbdata)

end subroutine get_state_vector


!------------------------------------------------------------------
!> Writes the current time and state variables from a dart state
!> vector (1d array) into a parflow binary file.

subroutine write_parflow_file(sv, pfb_state, dart_file, newfile)

real(r8),         intent(in) :: sv(:)           ! the DART posterior
real(r8),         intent(in) :: pfb_state(:)    ! diagnostic vector
character(len=*), intent(in) :: dart_file       ! the filename
character(len=*), intent(in) :: newfile         ! the name of the new parflow restart

character(len=*), parameter :: routine = 'write_parflow_file'
real(r8)                     :: rbuf(nx*ny*nz) ! data to be read
logical                      :: desiredG  = .false.
integer                      :: ivar, izerr, soilID

!pfb
real(r8), parameter          :: max_press_head = 0.005_r8
real(r8)                     :: pfvar(nx,ny,nz)       ! array->DART state 
real(r8)                     :: pfarr(nx,ny,nz)       ! array->ParFlow state 
real(r8)                     :: a_vG, N_vG, Ssat, Sres, Sw  ! vG param
integer(kind=4)              :: iunit

integer(kind=4)              ::  i, j, k,               &
                                 ix, iy, iz, ixs, iys,  &
                                 ns, nnx, nny, nnz
integer(kind=4), parameter   :: nxs = 1
integer(kind=4), parameter   :: nys = 1
integer(kind=4), parameter   :: rx  = 0
integer(kind=4), parameter   :: ry  = 0
integer(kind=4), parameter   :: rz  = 0

if ( .not. module_initialized ) call static_init_model

if (debug > 99 .and. do_output()) then
   write(string1,*) 'The DART posterior file is "'//trim(dart_file)//'"'
   write(string2,*) 'The new (posterior) parflow restart file is "'//trim(newfile)//'"'
   call error_handler(E_MSG, routine, string1, &
              source, revision, revdate, text2=string2)
endif

if ( .not. file_exist(dart_file) ) then 
   write(string1,*) 'cannot open file "', trim(dart_file),'" for reading.'
   call error_handler(E_ERR,routine,string1,source,revision,revdate)
endif

if (desiredG) rbuf = apply_clamping(ivar, sv)    ! For global clamping CPS

!vector to array conversion
pfarr = RESHAPE(pfb_state,(/nx,ny,nz/))

!vector to array conversion
pfvar = RESHAPE(sv,(/nx,ny,nz/))

iunit  = get_unit() 
open(iunit, file=newfile, status='unknown', access='stream', &
     convert='BIG_ENDIAN', form='unformatted', iostat=izerr)
if (izerr /= 0) then
   write(string1,*) 'cannot create file "', trim(newfile),'".'
   call error_handler(E_ERR,routine,string1,source,revision,revdate)
endif
rewind(iunit)

write(iunit) xdpfb !X
write(iunit) ydpfb !Y
write(iunit) zdpfb !Z

write(iunit) nx !NX
write(iunit) ny !NY
write(iunit) nz !NZ

write(iunit) dx !DX
write(iunit) dy !DY
write(iunit) dz !DZ

ns = INT(nxs*nys)
write(iunit) ns !num_subgrids
! End: Writing of domain spatial information

! Start: loop over number of sub grids
nnx = INT(nx/nxs)
nny = INT(ny/nys)
nnz = nz
iz = 0

do iys = 0, nys-1
do ixs = 0, nxs-1

! Start: Writing of sub-grid spatial information

  ix = INT(nnx*ixs)
  iy = INT(nny*iys)

  write(iunit) ix
  write(iunit) iy
  write(iunit) iz

  write(iunit) nnx
  write(iunit) nny
  write(iunit) nnz

  write(iunit) rx
  write(iunit) ry
  write(iunit) rz

! End: Writing of sub-grid spatial information

! Start: Write in data from each individual subgrid
! Follow the kji looping always
  do  k=iz +1 , iz + nnz
  do  j=iy +1 , iy + nny
  do  i=ix +1 , ix + nnx

    soilID = soil_type(i,j,k)
    a_vG   = soil_parameters(1,soilID)
    N_vG   = soil_parameters(2,soilID)
    Sres   = soil_parameters(3,soilID)
    Ssat   = soil_parameters(4,soilID)

    !CPS psi has to approach infinity to Sw to reach Sres
    Sw   = min(1._r8,max(pfvar(i,j,k),Sres + 0.001_r8))

    if (k .eq. nz .and. Sw .eq. 1._r8) then ! saturated
      pfarr(i,j,k) = max(pfarr(i,j,k),0._r8) 

    elseif (Sw.lt.1._r8 ) then  !UnsatZ
      pfarr(i,j,k) = -(1._r8/a_vG) * &
                ((((Ssat-Sres)/(Sw-Sres))**(N_vG/(N_vG-1._r8))) - 1._r8)**(1._r8/N_vG) 
    else
       write(string1,*)'Should not be possible.'
       call error_handler(E_ERR,routine,string1,source,revision,revdate)
    endif

    write(iunit) pfarr(i,j,k)
  end do
  end do
  end do
! End: Write in data from each individual subgrid

end do
end do
! End: loop over number of sub grids

close(iunit)

return

end subroutine write_parflow_file


!------------------------------------------------------------------
!> replace the posterior values that are outside physical limits with
!> the limit 

function apply_clamping(ivar, posterior) result (slab)

integer,  intent(in) :: ivar 
real(r8), intent(in) :: posterior(:)
real(r8)             :: slab(size(posterior))

slab = posterior

if ((progvar(ivar)%rangeRestricted == BOUNDED_ABOVE ) .or. &
    (progvar(ivar)%rangeRestricted == BOUNDED_BOTH )) then 
   where ((slab /= MISSING_R8) .and. &
          (slab > progvar(ivar)%maxvalue)) &
           slab = progvar(ivar)%maxvalue
endif

if ((progvar(ivar)%rangeRestricted == BOUNDED_BELOW ) .or. &
    (progvar(ivar)%rangeRestricted == BOUNDED_BOTH )) then 
   where ((slab /= MISSING_R8) .and. &
          (slab < progvar(ivar)%minvalue)) &
           slab = progvar(ivar)%minvalue
endif

end function apply_clamping


!------------------------------------------------------------------
!> returns the integer code defining what is being used for the DART state.
!> the id is set by interpreting the model_nml:parflow_assim_variable

function get_parflow_id()

integer :: get_parflow_id

if ( .not. module_initialized ) call static_init_model

get_parflow_id = ipfb

end function get_parflow_id


!------------------------------------------------------------------
!> Reads the parflow file name

function get_parflow_filename(pfbid)

integer, intent(in) :: pfbid 
character(len=256)  :: get_parflow_filename

if ( .not. module_initialized ) call static_init_model

if (pfbid .eq. 1) then
  get_parflow_filename = trim(parflow_press_file)
elseif (pfbid .eq. 2) then
  get_parflow_filename = trim(parflow_satur_file)
endif

end function get_parflow_filename


!-----------------------------------------------------------------
!> Reads the soilInd file to extract the parflow soil parameters
!> The soil_type matrix is an integer table of soil classifications.
!> The value in the matrix is an index into a 2D array of soil
!> parameters. 
!> int sID(lev, lat, lon) ;
!>     sID:units = "-" ;
!>     sID:longname = "ParFlow soil ID" ;
!>     sID:_FillValue = -9999 ;
!> double sparmval(sID, parameters) ;
!>     sparmval:units = "-" ;
!>     sparmval:longname = "Soil Parameters" ;
!>     sparmval:parameters = "a_vG N_vG Sres Sat Phi" ;
!>     sparmval:_FillValue = -9999.
!>
!> routine fills in module variables:
!>    soil_type(nx,ny,nz) .... integer soil types at each gridcell
!>    soil_parameters(nparams,:) ... 5 parameters for each integer type
!>    vcoord(nz) ... depth (m)

subroutine read_soil_table(filename)

character(len=*), intent(in) :: filename

character(len=*), parameter :: routine = 'read_soil_table'

character(len=obstypelength) :: dimnames(NF90_MAX_VAR_DIMS)

integer :: dimids(NF90_MAX_VAR_DIMS)
integer :: dimlens(NF90_MAX_VAR_DIMS)
integer :: ncid, io, ndims, varid 
integer :: i, iz
real(r8):: pfb_dz(nz) 

io = nf90_open(filename, NF90_NOWRITE, ncid)
call nc_check(io, routine, 'cannot open "'//trim(filename)//'"')

io = nf90_inq_varid(ncid, 'sID', varid)
call nc_check(io, routine, 'inq_varid sID')

io = nf90_inquire_variable(ncid, varid, dimids=dimids, ndims=ndims)
do i = 1,ndims
   write(string1,*)'inquire dimension ',i,' of "sID"'
   io = nf90_inquire_dimension(ncid, dimids(i), name=dimnames(i), len=dimlens(i))
   call nc_check(io, routine, string1)
enddo

if (dimlens(1) /= nx .or. dimlens(2) /=ny .or. dimlens(3) /=nz) then
  write(string1, *) 'dimension mismatch ... nlon ', dimlens(1), ' /= ', nx, ' or '
  write(string2, *) 'nlat ', dimlens(2), ' /= ', ny, ' or '
  write(string3, *) 'nlev ', dimlens(3), ' /= ', nz
  call error_handler(E_ERR,routine,string1, &
             source, revision, revdate, text2=string2, text3=string3)
endif

! Now we know the variable is shaped the same as the parflow state.

allocate(soil_type(nx,ny,nz))

io = nf90_get_var(ncid, varid, soil_type)
call nc_check(io, routine, 'get_var soil_type "'//trim(filename)//'"')

! move on to the lookup table of soil properties

io = nf90_inq_varid(ncid, 'sparmval', varid)
call nc_check(io, routine, 'inq_varid sparmval "'//trim(filename)//'"')

io = nf90_inquire_variable(ncid, varid, dimids=dimids, ndims=ndims)
do i = 1,ndims
   write(string1,*)'inquire dimension ',i,' of sparmval "'//trim(filename)//'"'
   io = nf90_inquire_dimension(ncid, dimids(i), name=dimnames(i), len=dimlens(i))
   call nc_check(io, routine, string1)
enddo

allocate(soil_parameters(dimlens(1),dimlens(2)))

io = nf90_get_var(ncid, varid, soil_parameters)
call nc_check(io, routine, 'get_var soil_parameters "'//trim(filename)//'"')

if (debug > 0 .and. do_output()) then
   do i = 1,size(soil_parameters,1)
     write(*,'(A,1x,i8,A,5(1X,F5.3))') 'parameters for soil type ',i, &
                                       ' are', soil_parameters(:,i)
   enddo
endif

! move on to the lookup table of soil properties

io = nf90_inq_varid(ncid, 'dz', varid)
call nc_check(io, routine, 'inquire_varid dz "'//trim(filename)//'"')

io = nf90_get_var(ncid, varid, pfb_dz)
call nc_check(io, routine, 'get_var dz "'//trim(filename)//'"')

io = nf90_close(ncid)
call nc_check(io, routine, 'close "'//trim(filename)//'"')

allocate(vcoord(nz))

do iz = nz, 1, -1
  if (iz == nz) then 
    vcoord(iz) = 0.5_r8 * pfb_dz(iz) 
  else 
    vcoord(iz) = vcoord(iz+1) + 0.5_r8 *(pfb_dz(iz+1) + pfb_dz(iz))
  end if
  if (debug > 10 .and. do_output()) &
    write(*,'(A,1X,I2,1X,F5.2,1X,F6.3)') '...dz_sID ', iz, pfb_dz(iz), vcoord(iz)
end do

end subroutine read_soil_table


!------------------------------------------------------------------
!> Reads the oasis grids.nc file to extract parflow geo-location 

subroutine grid_read(filename)

! gpfl.lon(y_gpfl,x_gpfl)
! gpfl.lat(y_gpfl,x_gpfl)

character(len=*), intent(in)   :: filename

character(len=*), parameter :: routine = 'grid_read'
integer :: ncid, io, nlon, nlat, dimid(2), varid(2) 

io = nf90_open(filename, NF90_NOWRITE, ncid)
call nc_check(io, routine, 'opening "'//trim(filename)//'"')

io = nf90_inq_dimid(ncid, "x_gpfl", dimid(1))
call nc_check(io, routine, 'inq_dimid x_gpfl "'//trim(filename)//'"')

io = nf90_inq_dimid(ncid, "y_gpfl", dimid(2))
call nc_check(io, routine, 'inq_dimid y_gpfl "'//trim(filename)//'"')

io = nf90_inq_varid(ncid, "gpfl.lon", varid(1))
call nc_check(io, routine, 'inq_varid gpfl.lon "'//trim(filename)//'"')

io = nf90_inq_varid(ncid, "gpfl.lat", varid(2))
call nc_check(io, routine, 'inq_varid gpfl.lat'//trim(filename)//'"')

io = nf90_inquire_dimension(ncid, dimid(1), string1, len=nlon)
call nc_check(io, routine, 'inquire_dimension nlon "'//trim(filename)//'"')

io = nf90_inquire_dimension(ncid, dimid(2), string1, len=nlat)
call nc_check(io, routine, 'inquire_dimension nlat "'//trim(filename)//'"')

if (nlon /= nx .or. nlat /=ny) then
  write(*, *) 'Dimensions ...', nlon, nx, nlat, ny
  call error_handler(E_ERR,routine,'ERR dimension mismatch',source, revision,revdate)
end if

io = nf90_get_var(ncid, varid(1), lon)
call nc_check(io, routine, 'get_var lon "'//trim(filename)//'"')

io = nf90_get_var(ncid, varid(2), lat)
call nc_check(io, routine, 'get_var lat "'//trim(filename)//'"')

where(lon <   0.0_r8) lon = lon + 360.0_r8
where(lat < -90.0_r8) lat = -90.0_r8
where(lat >  90.0_r8) lat =  90.0_r8

if (debug > 99 .and. do_output()) then
  write(*, *) 'grid_read:Longitude range ...', minval(lon), maxval(lon)
  write(*, *) 'grid_read:Latitude  range ...', minval(lat), maxval(lat)
endif

io = nf90_close(ncid)
call nc_check(io,routine,'closing "'//trim(filename)//'"')

end subroutine grid_read


!-----------------------------------------------------------------------
!> Reads pbf dimensions. based on tr32-z4-tools work of P. Shrestha

subroutine pfread_dim(filename)

character(len=*), intent(in)   :: filename

character(len=*), parameter :: routine = 'pfread_dim'
integer(kind=4)             :: nudat, izerr
character(len=256)          :: errmsg

! code starts here
  nudat   = get_unit() 

  open(nudat,file=trim(filename),form='unformatted',access='stream' , &
                    convert='BIG_ENDIAN',status='old', iostat=izerr)         ! gfortran
  if (izerr /= 0) then
     write(string1,*) 'cannot open "', trim(filename),'" for reading.'
     call error_handler(E_ERR,routine,string1,source,revision,revdate)
  endif

  !read in header infor
  read(nudat, iostat=izerr) xdpfb !X
  if (izerr /= 0) then
    errmsg   = "unable to read X"
    call error_handler(E_ERR,routine,errmsg,source,revision,revdate)
  endif
  read(nudat, iostat=izerr) ydpfb !Y
  if (izerr /= 0) then
    errmsg   = "unable to read Y"
    call error_handler(E_ERR,routine,errmsg,source,revision,revdate)
  endif
  read(nudat, iostat=izerr) zdpfb !Z
  if (izerr /= 0) then
    errmsg   = "unable to read Z"
    call error_handler(E_ERR,routine,errmsg,source,revision,revdate)
  endif

  read(nudat, iostat=izerr) nx !NX
  if (izerr /= 0) then
    errmsg   = "unable to read NX"
    call error_handler(E_ERR,routine,errmsg,source,revision,revdate)
  endif
  if (nx > 9999 ) then
     errmsg   = "problem readng NX (NX>9999), check pfb file"
     call error_handler(E_ERR,routine,errmsg,source,revision,revdate)
  endif
  read(nudat, iostat=izerr) ny !NY
  if (izerr /= 0) then
    errmsg   = "unable to read NY"
    call error_handler(E_ERR,routine,errmsg,source,revision,revdate)
  endif
  read(nudat, iostat=izerr) nz !NZ
  if (izerr /= 0) then
    errmsg   = "unable to read NZ"
    call error_handler(E_ERR,routine,errmsg,source,revision,revdate)
  endif

  read(nudat, iostat=izerr) dx !dX
  if (izerr /= 0) then
    errmsg   = "unable to read dx"
    call error_handler(E_ERR,routine,errmsg,source,revision,revdate)
  endif
  read(nudat, iostat=izerr) dy !dY
  if (izerr /= 0) then
    errmsg   = "unable to read dy"
    call error_handler(E_ERR,routine,errmsg,source,revision,revdate)
  endif
  read(nudat, iostat=izerr) dz !dZ
  if (izerr /= 0) then
    errmsg   = "unable to read dz"
    call error_handler(E_ERR,routine,errmsg,source,revision,revdate)
  endif

  close(nudat)
end subroutine pfread_dim


!-----------------------------------------------------------------------
!> Reads pfb binary file, based on tr32-z4-tools work of P. Shrestha
!>

subroutine pfread_var(filename,pfvar)

character(len=*), intent(in)  :: filename
real(r8),         intent(out) ::  pfvar(nx,ny,nz)       ! ParFlow pressure files 

character(len=*), parameter :: routine = 'pfread_var'
real(r8)                       :: dummyRes              ! Grid Resolution in m 
integer(kind=4)                :: nudat, dummy, izerr
integer(kind=4)                :: i, j, k,            &
                                  ix, iy, iz, is,     &
                                  ns, nnx, nny, nnz,  &
                                  rx, ry, rz
character(len=256)             :: errmsg

! code starts here
  nudat   = get_unit()
  open(nudat,file=trim(filename),form='unformatted',access='stream' , &
                    convert='BIG_ENDIAN',status='old', iostat=izerr) ! gfortran
   
  if (izerr /= 0) then
     write(string1,*) 'cannot open "', trim(filename),'" for reading.'
     call error_handler(E_ERR,routine,string1,source,revision,revdate)
  endif

  !read in header infor
  read(nudat, iostat=izerr) dummyRes  !X
  if (izerr /= 0) then
    errmsg   = "unable to read X"
    call error_handler(E_ERR,routine,errmsg,source,revision,revdate)
  endif
  read(nudat, iostat=izerr) dummyRes  !Y
  if (izerr /= 0) then
    errmsg   = "unable to read Y"
    call error_handler(E_ERR,routine,errmsg,source,revision,revdate)
  endif
  read(nudat, iostat=izerr) dummyRes  !Z
  if (izerr /= 0) then
    errmsg   = "unable to read Z"
    call error_handler(E_ERR,routine,errmsg,source,revision,revdate)
  endif

  read(nudat, iostat=izerr) dummy !NX
  if (izerr /= 0) then
    errmsg   = "unable to read NX"
    call error_handler(E_ERR,routine,errmsg,source,revision,revdate)
  endif
  read(nudat, iostat=izerr) dummy !NY
  if (izerr /= 0) then
    errmsg   = "unable to read NY"
    call error_handler(E_ERR,routine,errmsg,source,revision,revdate)
  endif
  read(nudat, iostat=izerr) dummy !NZ
  if (izerr /= 0) then
    errmsg   = "unable to read NZ"
    call error_handler(E_ERR,routine,errmsg,source,revision,revdate)
  endif

  read(nudat, iostat=izerr) dummyRes !dX
  if (izerr /= 0) then
    errmsg   = "unable to read dx"
    call error_handler(E_ERR,routine,errmsg,source,revision,revdate)
  endif
  read(nudat, iostat=izerr) dummyRes !dY
  if (izerr /= 0) then
    errmsg   = "unable to read dy"
    call error_handler(E_ERR,routine,errmsg,source,revision,revdate)
  endif
  read(nudat, iostat=izerr) dummyRes !dZ
  if (izerr /= 0) then
    errmsg   = "unable to read dz"
    call error_handler(E_ERR,routine,errmsg,source,revision,revdate)
  endif

  read(nudat, iostat=izerr) ns !num_subgrids
  if (izerr /= 0) then
    errmsg   = "unable to read ns"
    call error_handler(E_ERR,routine,errmsg,source,revision,revdate)
  endif
!
! Start loop over number of subgrids
  do is = 0, ns-1
! Start: reading subgrid spatial information
! ix, iy, iz
   read(nudat, iostat=izerr) ix
   if (izerr /= 0) then
     errmsg   = "unable to read ix"
     call error_handler(E_ERR,routine,errmsg,source,revision,revdate)
   endif
   read(nudat, iostat=izerr) iy
   if (izerr /= 0) then
     errmsg   = "unable to read iy"
     call error_handler(E_ERR,routine,errmsg,source,revision,revdate)
   endif
   read(nudat, iostat=izerr) iz
   if (izerr /= 0) then
     errmsg   = "unable to read iz"
     call error_handler(E_ERR,routine,errmsg,source,revision,revdate)
   endif

! nnx. nny, nnz 
   read(nudat, iostat=izerr) nnx
   if (izerr /= 0) then
     errmsg   = "unable to read nnx"
     call error_handler(E_ERR,routine,errmsg,source,revision,revdate)
   endif
   read(nudat, iostat=izerr) nny
   if (izerr /= 0) then
     errmsg   = "unable to read nny"
     call error_handler(E_ERR,routine,errmsg,source,revision,revdate)
   endif
   read(nudat, iostat=izerr) nnz
   if (izerr /= 0) then
     errmsg   = "unable to read nnz"
     call error_handler(E_ERR,routine,errmsg,source,revision,revdate)
   endif

! rx,ry,rz
   read(nudat, iostat=izerr) rx
   if (izerr /= 0) then
     errmsg   = "unable to read rx"
     call error_handler(E_ERR,routine,errmsg,source,revision,revdate)
   endif
   read(nudat, iostat=izerr) ry
   if (izerr /= 0) then
     errmsg   = "unable to read ry"
     call error_handler(E_ERR,routine,errmsg,source,revision,revdate)
   endif
   read(nudat, iostat=izerr) rz
   if (izerr /= 0) then
     errmsg   = "unable to read rz"
     call error_handler(E_ERR,routine,errmsg,source,revision,revdate)
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


!-----------------------------------------------------------------------
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
!> Horizontal interpolation for each of two heights.
!> Result is two values - one at each height.

subroutine horizontal_interpolate(x, ivar, kk_inds, ij_inds, ij_wgts, interp_hv, istatus)

real(r8), intent(in)  :: x(:)
integer,  intent(in)  :: ivar
integer,  intent(in)  :: kk_inds(2)
integer,  intent(in)  :: ij_inds(4)    ! ileft, iright, jbot, jtop
real(r8), intent(in)  :: ij_wgts(2)    ! ifrac, jfrac
real(r8), intent(out) :: interp_hv(2)  ! at two height levels 
integer,  intent(out) :: istatus

integer  :: k, ll, lr, ur, ul
real(r8) :: ij_rwgts(2)

! Interpolate
ij_rwgts(1) = 1.0_r8 - ij_wgts(1)  !ifrac
ij_rwgts(2) = 1.0_r8 - ij_wgts(2)  !jfrac

do k = 1, 2
  ! Location of dart_indices
  ll = ijk_to_dart(ivar, ij_inds(1), ij_inds(3), kk_inds(k))
  lr = ijk_to_dart(ivar, ij_inds(2), ij_inds(3), kk_inds(k))
  ur = ijk_to_dart(ivar, ij_inds(2), ij_inds(4), kk_inds(k))
  ul = ijk_to_dart(ivar, ij_inds(1), ij_inds(4), kk_inds(k))

  interp_hv(k) = ij_rwgts(2) * ( ij_rwgts(1) * x(ll) + ij_wgts(1) * x(lr)) + &
                 ij_wgts(2)  * ( ij_rwgts(1) * x(ul) + ij_wgts(1) * x(ur))

end do

if (debug > 99 .and. do_output()) then
   write(*,*)
   write(*,*)' ileft, jbot, level_index decompose to ', ll, x(ll), ij_rwgts(1)  
   write(*,*)'iright, jbot, level_index decompose to ', lr, x(lr), ij_wgts(1)
   write(*,*)'iright, jtop, level_index decompose to ', ur, x(ur), ij_wgts(1)  
   write(*,*)' ileft, jtop, level_index decompose to ', ul, x(ul), ij_rwgts(1) 
   write(*,*)' bottom, top weights are ', ij_rwgts(2), ij_rwgts(1) 
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

end function ijk_to_dart


!------------------------------------------------------------------------
!> Get model corners for the observed location
!> Adapted from cosmo interp routines

subroutine get_corners(lon, lat, nx, ny, gridlons, gridlats, &
                       ij_inds, ij_wgts, istatus)

real(r8), intent(in)  :: lon
real(r8), intent(in)  :: lat
integer,  intent(in)  :: nx
integer,  intent(in)  :: ny
real(r8), intent(in)  :: gridlons(nx)
real(r8), intent(in)  :: gridlats(ny)
integer,  intent(out) :: ij_inds(4) ! ileft, iright, jbot, jtop
real(r8), intent(out) :: ij_wgts(2) ! ifrac, jfrac

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
! 'dlon' is ij_wghts(1)

if (lon == gridlons(nx)) then
   ij_inds(1)  = nx-1
   ij_inds(2)  = nx
   ij_wgts(1)  = 0.0_r8  ! weight for the lonleft
else
   indarr      = minloc(abs(gridlons-lon))   
   if (gridlons(indarr(1)) .gt. lon) then
     ij_inds(1)  = indarr(1) - 1
   else
     ij_inds(1)  = indarr(1)
   endif
   ij_inds(2)  = ij_inds(1) + 1
   dlon        =     lon              - gridlons(ij_inds(1)) ! see figure
   dlonT       = gridlons(ij_inds(2)) - gridlons(ij_inds(1))
   ij_wgts(1)  = dlon / dlonT
endif

! latitudes are arranged 'south' to 'north', i.e. -90 to 90 (geographical) 
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
   indarr      = minloc(abs(gridlats-lat))
    if (gridlats(indarr(1)) .gt. lat) then 
     ij_inds(3)  = indarr(1) - 1
   else 
     ij_inds(3)  = indarr(1)
   endif
   ij_inds(4)   = ij_inds(3) + 1
   dlat         =           lat        - gridlats(ij_inds(3))
   dlatT        = gridlats(ij_inds(4)) - gridlats(ij_inds(3))
   ij_wgts(2)   = dlat / dlatT
endif

istatus = 0

if (debug > 5 .and. do_output()) then
   write(*,*)
   write(*,*)'longitude index to the  west and  east are ',ij_inds(1:2)
   write(*,*)'latitude  index to the south and north are ',ij_inds(3:4) 
   write(*,*)'lonW, lon, lonE',gridlons(ij_inds(1)),lon,gridlons(ij_inds(2))
   write(*,*)'latS, lat, latN',gridlats(ij_inds(3)), lat,gridlats(ij_inds(4))
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
!               ideeper      height                 ishallower
!              kk_inds(1)                          kk_inds(2) 

real(r8), intent(in)  :: height
integer,  intent(out) :: kk_inds(2) 
real(r8), intent(out) :: kk_wgt 
integer,  intent(out) :: istatus

integer, parameter :: OUTSIDE_VERTICALLY   = 16

integer :: indarr(1)
real(r8) :: dz, dztot

! Variable definition complete
if (height > vcoord(1) .or. height < vcoord(nz)) then
  istatus = OUTSIDE_VERTICALLY
  return
endif

if (height == vcoord(1)) then
  kk_inds(1) = 1
  kk_inds(2) = 1
  kk_wgt     = 0._r8
elseif (height == vcoord(nz)) then
  kk_inds(1) = nz
  kk_inds(2) = nz   
  kk_wgt     = 0._r8
else
  indarr  = minloc( vcoord,  vcoord > height )
  kk_inds(1) = indarr(1)       ! Deeper layer
  kk_inds(2) = kk_inds(1) + 1  ! Shallower layer
  dz     = vcoord(kk_inds(1)) - height
  dztot  = vcoord(kk_inds(1)) - vcoord(kk_inds(2))
  kk_wgt = dz / dztot
endif

if (debug > 99 .and. do_output()) then
  write(*,*)
  write(*,'(A30,f20.15,1x,i2)')'height at shallower layer',vcoord(kk_inds(2)), kk_inds(2)
  write(*,*)'frac',kk_wgt
  write(*,'(A30,f20.15)')'desired        height ',height
  write(*,*)'frac',1.0_r8 - kk_wgt
  write(*,'(A30,f20.15,1x,i2)')'height at deeper    layer',vcoord(kk_inds(1)), kk_inds(1)
  write(*,*)
endif

istatus = 0

end subroutine get_level_indices


!------------------------------------------------------------------
!> The CLM restart netcdf files have the time of the parflow state.

function get_model_time( filename )

type(time_type) :: get_model_time
character(len=*), intent(in) :: filename

character(len=*), parameter :: routine = 'get_model_time'
integer :: io, ncid, VarID
integer :: rst_curr_ymd, rst_curr_tod, leftover
integer :: year, month, day, hour, minute, second

if ( .not. file_exist(filename) ) then 
  write(string1,*) 'cannot open "', trim(filename),'" for reading.'
  call error_handler(E_ERR,routine,string1,source,revision,revdate)
endif

io = nf90_open(filename, NF90_NOWRITE, ncid)
call nc_check(io, routine, 'open "'//trim(filename)//'"')

io = nf90_inq_varid(ncid, 'timemgr_rst_curr_ymd', VarID)
call nc_check(io, routine, 'inq_varid timemgr_rst_curr_ymd "'//trim(filename)//'"')

io = nf90_get_var(ncid, VarID, rst_curr_ymd)
call nc_check(io, routine, 'get_var rst_curr_ymd "'//trim(filename)//'"')

io = nf90_inq_varid(ncid, 'timemgr_rst_curr_tod', VarID)
call nc_check(io, routine, 'inq_varid timemgr_rst_curr_tod "'//trim(filename)//'"')

io = nf90_get_var(ncid, VarID, rst_curr_tod)
call nc_check(io, routine, 'get_var rst_curr_tod "'//trim(filename)//'"')

io = nf90_close(ncid)
call nc_check(io, routine, 'close "'//trim(filename)//'"')

year     = rst_curr_ymd/10000
leftover = rst_curr_ymd - year*10000
month    = leftover/100
day      = leftover - month*100

hour     = rst_curr_tod/3600
leftover = rst_curr_tod - hour*3600
minute   = leftover/60
second   = leftover - minute*60

get_model_time = set_date(year, month, day, hour, minute, second)

end function get_model_time


!------------------------------------------------------------------
!>Get start date from CLM restart netcdf file

function get_start_time(filename)
character(len=*), intent(in) :: filename
type(time_type) :: get_start_time

character(len=*), parameter :: routine = 'get_start_time'
integer :: io, ncid, VarID
integer :: rst_start_ymd, rst_start_tod, leftover
integer :: year, month, day, hour, minute, second

if ( .not. file_exist(filename) ) then 
   write(string1,*) 'cannot open "', trim(filename),'" to read start date.'
   call error_handler(E_ERR,routine,string1,source,revision,revdate)
endif

io = nf90_open(trim(filename), NF90_NOWRITE, ncid)
call nc_check(io, routine,'open "'//trim(filename)//'"')

io = nf90_inq_varid(ncid, 'timemgr_rst_ref_ymd', VarID)
call nc_check(io, routine, 'inq_varid timemgr_rst_ref_ymd "'//trim(filename)//'"')

io = nf90_get_var(ncid, VarID, rst_start_ymd)
call nc_check(io, routine, 'get_var rst_ref_ymd "'//trim(filename)//'"')

io = nf90_inq_varid(ncid, 'timemgr_rst_ref_tod', VarID)
call nc_check(io, routine, 'inq_varid timemgr_rst_ref_tod "'//trim(filename)//'"')

io = nf90_get_var(  ncid, VarID,   rst_start_tod)
call nc_check(io, routine, 'get_var rst_ref_tod "'//trim(filename)//'"')

io = nf90_close(ncid)
  call nc_check(io, routine,'close "'//trim(filename)//'"')

year     = rst_start_ymd/10000
leftover = rst_start_ymd - year*10000
month    = leftover/100
day      = leftover - month*100

hour     = rst_start_tod/3600
leftover = rst_start_tod - hour*3600
minute   = leftover/60
second   = leftover - minute*60

get_start_time = set_date(year, month, day, hour, minute, second)

end function get_start_time


!-----------------------------------------------------------------------
!> routine to compute the soil moisture at all the surrounding corners
!> and interpolate them to two locations - above and below 
!> Returns soil moisture in cm^3 / cm^3
!>
!> The porosity is required at the grid locations to convert soil saturations
!> to soil moisture. The soil moisture at the locations is then linearly
!> interpolated to the desired locations. The porosity at the gridcells can
!> be obtained from the file referenced in the &model_nml:soilInd_file variable.
!> Necessary module variables containing this metadata are:
!>
!> soil_type
!> soil_parameters

subroutine calculate_soil_moisture(x, ivar, kk_inds, ij_inds, ij_wgts, interp_hv, istatus)

real(r8), intent(in)  :: x(:)            ! soil saturation (the state)
integer,  intent(in)  :: ivar
integer,  intent(in)  :: kk_inds(2)      ! layer indices
integer,  intent(in)  :: ij_inds(4)      ! ileft, iright, jbot, jtop
real(r8), intent(in)  :: ij_wgts(2)      ! ifrac, jfrac
real(r8), intent(out) :: interp_hv(2)    ! at two levels
integer,  intent(out) :: istatus

integer, parameter :: POROSITY = 5
integer  :: soilID
integer  :: k, ll, lr, ur, ul
real(r8) :: ij_rwgts(2)
real(r8) :: saturation(4)   ! at each corner
real(r8) :: porosities(4)   ! at each corner

ij_rwgts(1) = 1.0_r8 - ij_wgts(1)  !ifrac
ij_rwgts(2) = 1.0_r8 - ij_wgts(2)  !jfrac

! calculate the soil moisture at each corner ... for each layer

do k = 1, 2

  ! Location of dart_indices
  ll = ijk_to_dart(ivar, ij_inds(1), ij_inds(3), kk_inds(k))
  lr = ijk_to_dart(ivar, ij_inds(2), ij_inds(3), kk_inds(k))
  ur = ijk_to_dart(ivar, ij_inds(2), ij_inds(4), kk_inds(k))
  ul = ijk_to_dart(ivar, ij_inds(1), ij_inds(4), kk_inds(k))

  saturation = (/ x(ll), x(lr), x(ur), x(ul) /)

  soilID = soil_type(ij_inds(1), ij_inds(3), kk_inds(k)) 
  porosities(1) = soil_parameters(POROSITY,soilID) ! lower left

  soilID = soil_type(ij_inds(2), ij_inds(3), kk_inds(k)) 
  porosities(2) = soil_parameters(POROSITY,soilID) ! lower right
  
  soilID = soil_type(ij_inds(2), ij_inds(4), kk_inds(k)) 
  porosities(3) = soil_parameters(POROSITY,soilID) ! upper right
  
  soilID = soil_type(ij_inds(1), ij_inds(4), kk_inds(k)) 
  porosities(4) = soil_parameters(POROSITY,soilID) ! upper left

  saturation = saturation * porosities ! now soil moisture, actually
 
  interp_hv(k) = ij_rwgts(2) * (ij_rwgts(1)*saturation(1) + ij_wgts(1)*saturation(2)) + &
                 ij_wgts(2)  * (ij_rwgts(1)*saturation(4) + ij_wgts(1)*saturation(3))
enddo

if (debug > 99 .and. do_output()) then
   write(*,*)
   write(*,*)' ileft, jbot, level_index decompose to ', ll, x(ll), ij_rwgts(1)  
   write(*,*)'iright, jbot, level_index decompose to ', lr, x(lr), ij_wgts(1)
   write(*,*)'iright, jtop, level_index decompose to ', ur, x(ur), ij_rwgts(2)  
   write(*,*)' ileft, jtop, level_index decompose to ', ul, x(ul), ij_rwgts(2) 
   write(*,*)' bottom, top weights are ', ij_rwgts(2), ij_rwgts(1) 
   write(*,*)' horizontal value is ', kk_inds(1), interp_hv(1)
   write(*,*)' horizontal value is ', kk_inds(2), interp_hv(2)
endif

istatus = 0

end subroutine calculate_soil_moisture


!===================================================================
end module model_mod
!===================================================================
