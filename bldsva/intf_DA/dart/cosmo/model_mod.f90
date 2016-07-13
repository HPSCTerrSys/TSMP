! DART software - Copyright 2004 - 2016 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! DART $Id: $

module model_mod

! This module provides routines to work with COSMO data
! files in the DART framework

use        types_mod, only : i4, r4, r8, digits12, SECPERDAY, MISSING_R8,          &
                             rad2deg, deg2rad, PI, obstypelength

use time_manager_mod, only : time_type, set_time, set_date, get_date, get_time,&
                             print_time, print_date, set_calendar_type,        &
                             operator(*),  operator(+), operator(-),           &
                             operator(>),  operator(<), operator(/),           &
                             operator(/=), operator(<=)

use     location_mod, only : location_type, get_dist, query_location,          &
                             get_close_maxdist_init, get_close_type,           &
                             set_location, get_location, horiz_dist_only,      &
                             vert_is_undef,        VERTISUNDEF,                &
                             vert_is_surface,      VERTISSURFACE,              &
                             vert_is_level,        VERTISLEVEL,                &
                             vert_is_pressure,     VERTISPRESSURE,             &
                             vert_is_height,       VERTISHEIGHT,               &
                             vert_is_scale_height, VERTISSCALEHEIGHT,          &
                             get_close_obs_init, get_close_obs

use    utilities_mod, only : register_module, error_handler,                   &
                             E_ERR, E_WARN, E_MSG, logfileunit, get_unit,      &
                             nc_check, do_output, to_upper,                    &
                             find_namelist_in_file, check_namelist_read,       &
                             open_file, close_file, file_exist,                &
                             find_textfile_dims, file_to_text,                 &
                             do_nml_file, do_nml_term

use     obs_kind_mod, only : KIND_U_WIND_COMPONENT,       &
                             KIND_V_WIND_COMPONENT,       &
                             KIND_VERTICAL_VELOCITY,      &
                             KIND_TEMPERATURE,            &
                             KIND_PRESSURE,               &
                             KIND_PRESSURE_PERTURBATION,  &
                             KIND_SPECIFIC_HUMIDITY,      &
                             KIND_CLOUD_LIQUID_WATER,     &
                             KIND_CLOUD_ICE,              &
                             KIND_SURFACE_ELEVATION,      &
                             KIND_SURFACE_GEOPOTENTIAL,   &
                             paramname_length,            &
                             get_raw_obs_kind_index,      &
                             get_raw_obs_kind_name

use    random_seq_mod, only: random_seq_type, init_random_seq, random_gaussian

use          sort_mod, only: index_sort

use netcdf

implicit none
private

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = "$URL: model_mod.f90 $"
character(len=32 ), parameter :: revision = "$Revision: none $"
character(len=128), parameter :: revdate  = "$Date: none $"

character(len=256) :: string1, string2, string3
logical, save :: module_initialized = .false.

! these routines must be public and you cannot change
! the arguments - they will be called *from* the DART code.

public :: get_model_size,         &
          adv_1step,              &
          get_state_meta_data,    &
          model_interpolate,      &
          get_model_time_step,    &
          static_init_model,      &
          end_model,              &
          init_time,              &
          init_conditions,        &
          nc_write_model_atts,    &
          nc_write_model_vars,    &
          pert_model_state,       &
          get_close_maxdist_init, &
          get_close_obs_init,     &
          get_close_obs,          &
          ens_mean_for_model

! generally useful routines for various support purposes.
! the interfaces here can be changed as appropriate.

public :: get_state_time,     &
          get_state_vector,   &
          write_cosmo_file,    &
          get_cosmo_filename, &
          write_state_times

! TODO FIXME  ultimately the dart_variable_info type will be removed
type dart_variable_info
   character(len=16)    :: varname_short
   character(len=256)   :: varname_long
   character(len=32)    :: units
   logical              :: is_present
   integer              :: nx
   integer              :: ny
   integer              :: nz
   real(r8),allocatable :: vertical_level(:)
   integer              :: vertical_coordinate
   integer              :: horizontal_coordinate
   integer,allocatable  :: state_vector_sindex(:) ! starting index in state vector for every vertical level
   integer,allocatable  :: cosmo_state_index(:)   ! index in cosmo state of every vertical level
end type dart_variable_info

! Codes for restricting the range of a variable
integer, parameter :: BOUNDED_NONE  = 0 ! ... unlimited range
integer, parameter :: BOUNDED_BELOW = 1 ! ... minimum, but no maximum
integer, parameter :: BOUNDED_ABOVE = 2 ! ... maximum, but no minimum
integer, parameter :: BOUNDED_BOTH  = 3 ! ... minimum and maximum

integer :: nfields
integer, parameter :: max_state_variables = 80
integer, parameter :: num_state_table_columns = 7
character(len=obstypelength) :: variable_table(max_state_variables, num_state_table_columns)

! Codes for interpreting the columns of the variable_table
integer, parameter :: VT_GRIBVERSIONINDX = 1 ! ... the version of the grib table being used
integer, parameter :: VT_GRIBVARINDX     = 2 ! ... variable name
integer, parameter :: VT_VARNAMEINDX     = 3 ! ... netcdf variable name
integer, parameter :: VT_KINDINDX        = 4 ! ... DART kind
integer, parameter :: VT_MINVALINDX      = 5 ! ... minimum value if any
integer, parameter :: VT_MAXVALINDX      = 6 ! ... maximum value if any
integer, parameter :: VT_STATEINDX       = 7 ! ... update (state) or not

integer, parameter :: NPDS  = 321 ! dimension of product definition section
integer, parameter :: NGDS  = 626 ! dimension of grid description section

! Everything needed to describe a variable
!>@ TODO FIXME remove the unused netcdf bits ... we're working with binary only

type progvartype
   private
   character(len=NF90_MAX_NAME) :: varname
   character(len=NF90_MAX_NAME) :: long_name
   character(len=NF90_MAX_NAME) :: units
   character(len=obstypelength) :: dimnames(NF90_MAX_VAR_DIMS)
   integer  :: dimlens(NF90_MAX_VAR_DIMS)
   integer  :: tableID     ! grib table version
   integer  :: variableID  ! variable ID in grib table
   integer  :: levtypeID   ! the kind of vertical coordinate system
   integer  :: numdims
   integer  :: numEW
   integer  :: numNS
   integer  :: numZ
   integer  :: varsize     ! prod(dimlens(1:numdims))
   integer  :: index1      ! location in dart state vector of first occurrence
   integer  :: indexN      ! location in dart state vector of last  occurrence
   integer  :: slab1
   integer  :: slabN
   integer  :: dart_kind
   integer  :: rangeRestricted
   real(r8) :: minvalue
   real(r8) :: maxvalue
   character(len=paramname_length) :: kind_string
   logical  :: update
end type progvartype

type(progvartype), dimension(max_state_variables) :: progvar

! Dimensions (from the netCDF file) specifying the grid shape, etc.
integer :: ntime
integer :: nbnds
integer :: nrlon
integer :: nrlat
integer :: nsrlon
integer :: nsrlat
integer :: nlevel
integer :: nlevel1
integer :: nsoil
integer :: nsoil1

!> @TODO FIXME ... check to make sure the storage order in netCDF is the
!> same as the storage order of the binary variables.
!> vcoord has units, etc. and will require a transformation

real(r8), allocatable ::   lon(:,:) !                  longitude degrees_east
real(r8), allocatable ::   lat(:,:) !                  latitude  degrees_north
real(r8), allocatable :: slonu(:,:) ! staggered U-wind longitude degrees_east
real(r8), allocatable :: slatu(:,:) ! staggered U-wind latitude  degrees_north
real(r8), allocatable :: slonv(:,:) ! staggered V-wind longitude degrees_east
real(r8), allocatable :: slatv(:,:) ! staggered V-wind latitude  degrees_north

! variables pertaining to the computational grid
real(r8), allocatable ::  rlon(:) !           rotated longitude degrees
real(r8), allocatable ::  rlat(:) !           rotated latitude  degrees
real(r8), allocatable :: srlon(:) ! staggered rotated longitude degrees
real(r8), allocatable :: srlat(:) ! staggered rotated latitude  degrees
real(r8) :: north_pole_latitude
real(r8) :: north_pole_longitude

! just mimic what is in the netCDF variable
!       float vcoord(level1) ;
!               vcoord:long_name = "Height-based hybrid Gal-Chen coordinate" ;
!               vcoord:units = "Pa" ;
!               vcoord:ivctype = 2 ;
!               vcoord:irefatm = 2 ;
!               vcoord:p0sl = 100000. ;
!               vcoord:t0sl = 300. ;
!               vcoord:dt0lp = 42. ;
!               vcoord:vcflat = 11000. ;
!               vcoord:delta_t = 75. ;
!               vcoord:h_scal = 10000. ;

type verticalobject
   private
   real(r8), allocatable :: level1(:)
   integer               :: nlevel1
   real(r8), allocatable :: level(:)
   integer               :: nlevel
   character(len=128)    :: long_name = "Height-based hybrid Gal-Chen coordinate"
   character(len=32)     :: units = "Pa"
   integer               :: ivctype = 2
   integer               :: irefatm = 2
   real(r8)              :: p0sl = 100000.
   real(r8)              :: t0sl = 300.
   real(r8)              :: dt0lp = 42.
   real(r8)              :: vcflat = 11000.
   real(r8)              :: delta_t = 75.
   real(r8)              :: h_scal = 10000.
end type verticalobject

type(verticalobject), save :: vcoord  ! some compilers require save if initializing a structure

! things which can/should be in the model_nml

character(len=256) :: cosmo_restart_file           = "cosmo_restart_file"
character(len=256) :: cosmo_netcdf_file            = "cosmo_netcdf_file"
integer            :: assimilation_period_days     = 0
integer            :: assimilation_period_seconds  = 60
integer            :: model_dt                     = 40
logical            :: output_1D_state_vector       = .FALSE.
real(r8)           :: model_perturbation_amplitude = 0.1
integer            :: debug                        = 0
character(len=obstypelength) :: variables(max_state_variables*num_state_table_columns) = ' '

namelist /model_nml/             &
   cosmo_restart_file,           &
   cosmo_netcdf_file,            &
   assimilation_period_days,     &
   assimilation_period_seconds,  &
   model_dt,                     &
   model_perturbation_amplitude, &
   output_1D_state_vector,       &
   debug,                        &
   variables

integer         :: model_size = 0
type(time_type) :: model_timestep ! smallest time to adv model

integer, parameter       :: n_max_kinds=400
type(dart_variable_info) :: state_vector_vars(1:n_max_kinds)

real(r8), allocatable :: ens_mean(:)

type(random_seq_type) :: random_seq

!> TODO ... do we need these
type(time_type) :: cosmo_fc_time
type(time_type) :: cosmo_an_time

INTERFACE get_grid_var
      MODULE PROCEDURE get_1d_grid_var
      MODULE PROCEDURE get_2d_grid_var
END INTERFACE


contains

!------------------------------------------------------------------------
!>


function get_model_size()

integer :: get_model_size

if ( .not. module_initialized ) call static_init_model

get_model_size = model_size

end function get_model_size


!------------------------------------------------------------------------
!> Called to do one-time initialization of the model.
!>
!> All the grid information comes from the COSMOS netCDF file
!>@ TODO FIXME All the variable information comes from all over the place
!> Not actually reading in the state, that is done in get_state_vector()


subroutine static_init_model()

! Local variables - all the important ones have module scope

integer :: io, iunit, ivar 

if ( module_initialized ) return ! only need to do this once.

module_initialized = .TRUE.

! read the DART namelist for this model
call find_namelist_in_file('input.nml', 'model_nml', iunit)
read(iunit, nml = model_nml, iostat = io)
call check_namelist_read(iunit, io, 'model_nml')

! Record the namelist values used for the run
if (do_nml_file()) write(logfileunit, nml=model_nml)
if (do_nml_term()) write(     *     , nml=model_nml)

call set_calendar_type('Gregorian')

! Get the dimensions of the grid and the grid variables from the netCDF file.
call get_cosmo_grid(cosmo_netcdf_file)

! rectify the user input for the variables to include in the DART state vector
call parse_variable_table(variables, nfields, variable_table)

! Read the input file to decode and confirm shapes of the variables.
! Without a second argument, read_binary_file just collects sizes.

!>@ do we need the model time here ... 
do ivar = 1,nfields
    call read_binary_file(ivar)
enddo

! calculate how we store what we have in the DART state vector

call set_variable_layout()

if (debug > 5 .and. do_output()) call progvar_summary()

! ens_mean_for_model, that sort of thing

return

end subroutine static_init_model


!------------------------------------------------------------------------
!> Given an index into the DART state vector, return the location and
!> (optionally) what kind of variable ... KIND_TEMPERATURE, KIND_U_WIND_COMPONENT, etc.


subroutine get_state_meta_data(index_in, location, var_type)

integer,             intent(in)  :: index_in
type(location_type), intent(out) :: location
integer, optional,   intent(out) :: var_type

real(r8) :: mylon,mylat,vloc
integer  :: nx, ny, nz
integer  :: local_ind
integer  :: iloc, jloc, kloc
integer  :: ivar

if (.not. module_initialized) call static_init_model()

ivar = Find_Variable_by_index(index_in,'get_state_meta_data')

! Determine the i,j,k location of interest

local_ind = index_in - progvar(ivar)%index1

nx = progvar(ivar)%dimlens(1)
ny = progvar(ivar)%dimlens(2)
nz = progvar(ivar)%dimlens(3)  ! may or may not be needed

if     (progvar(ivar)%numdims == 2) then
   kloc =  1
   jloc =  local_ind / nx + 1
   iloc =  local_ind - (jloc-1)*nx + 1
elseif (progvar(ivar)%numdims == 3) then
   kloc =  local_ind / (nx*ny) + 1
   jloc = (local_ind - (kloc-1)*nx*ny)/nx + 1
   iloc =  local_ind - (kloc-1)*nx*ny - (jloc-1)*nx + 1
else

   write(string1,*) 'Can not calculate indices for variable ', &
      trim(progvar(ivar)%varname), ' ndims = ', progvar(ivar)%numdims
   call error_handler(E_ERR, 'get_state_meta_data',string1)

endif

if (debug > 99 .and. do_output()) &
write(*,*)'gsmd: index_in ',index_in,' local_ind is ',local_ind + 1, ' and dereferences to ',iloc,jloc,kloc

! Now that we know the i,j,k we have to index the right set of
! coordinate arrays

if     (progvar(ivar)%dart_kind == KIND_U_WIND_COMPONENT) then
   mylon = slonu(iloc,jloc)
   mylat = slatu(iloc,jloc)
elseif (progvar(ivar)%dart_kind == KIND_V_WIND_COMPONENT) then
   mylon = slonv(iloc,jloc)
   mylat = slatv(iloc,jloc)
else
   mylon = lon(iloc,jloc)
   mylat = lat(iloc,jloc)
endif

if (progvar(ivar)%dart_kind == KIND_VERTICAL_VELOCITY) then
   vloc = vcoord%level1(kloc)
else
   vloc = vcoord%level(kloc)
endif

location = set_location(mylon, mylat, vloc, VERTISHEIGHT) ! meters

if (present(var_type)) then
   var_type = progvar(ivar)%dart_kind
endif

return

end subroutine get_state_meta_data


!------------------------------------------------------------------------
!> Returns the smallest increment of time that we want to advance the model.
!> This defines the minimum assimilation interval.
!> It is NOT the dynamical timestep of the model.


function get_model_time_step()

type(time_type) :: get_model_time_step

call error_handler(E_ERR,'get_model_time_step','routine not written',source,revision,revdate)

if ( .not. module_initialized ) call static_init_model

model_timestep      = set_time(model_dt)
get_model_time_step = model_timestep

return

end function get_model_time_step


!------------------------------------------------------------------------
!>


subroutine model_interpolate(x, location, obs_type, interp_val, istatus)

! Error codes:
! istatus = 99 : unknown error
! istatus = 10 : observation type is not in state vector
! istatus = 15 : observation lies outside the model domain (horizontal)
! istatus = 16 : observation lies outside the model domain (vertical)
! istatus = 19 : observation vertical coordinate is not supported
!
! The third argument would be more accurately called 'DART_KIND', but
! cannot/should not be changed as this is a mandatory interface.
! consequently, the first thing we will do is to assign the value
! to a variable that is more appropriately named ... 'dart_kind'

! Passed variables

real(r8),            intent(in)  :: x(:)
type(location_type), intent(in)  :: location
integer,             intent(in)  :: obs_type
real(r8),            intent(out) :: interp_val
integer,             intent(out) :: istatus

integer, parameter :: NO_KIND_IN_STATE     = 10
integer, parameter :: VERTICAL_UNSUPPORTED = 19

! Local storage

integer :: dart_kind
integer :: ivar  ! index into the progvar structure for this dart_kind

real(r8) :: point_coords(1:3)

real(r8), parameter :: polgam = 0.0_r8 ! angle between the north poles of the systems

real(r8) :: geo_lat, geo_lon, height
real(r8) :: rotated_lat, rotated_lon

integer  :: i, iabove, ibelow
integer  :: ileft, iright, jbot, jtop
real(r8) :: ifrac, jfrac, levelfrac

real(r8) :: value_above, value_below

IF ( .not. module_initialized ) call static_init_model

interp_val = MISSING_R8     ! the DART bad value flag
istatus = 99                ! unknown error

dart_kind = obs_type  ! compensate for a poor original choice of varname

! If identity observation (dart_kind < 0), then no need to interpolate
! identity observation -> -(dart_kind) = DART state vector index
! obtain "interpolated" value directly from state using (abs(dart_kind))

if ( dart_kind < 0 ) then
   interp_val = x(-1*dart_kind)
   istatus = 0
   return
endif

! Determine if this dart_kind exists in the DART state vector
! Sometimes special action needs to happen to interpolate to model levels.
! To determine the number of model levels, one can try to interpolate
! KIND_GEOPOTENTIAL_HEIGHT, KIND_GEOMETRIC_HEIGHT ... which are not 
! normally part of the state, so this next test will have to be modified.

ivar = -1
FoundIt: do i = 1,nfields
   if (progvar(i)%dart_kind == dart_kind) then
      ivar = i
      exit FoundIt
   endif
enddo FoundIt

if ( ivar < 0 ) then
   istatus = NO_KIND_IN_STATE
   return
endif

! Carry on

point_coords(1:3) = get_location(location)
geo_lon = point_coords(1) ! degrees East
geo_lat = point_coords(2) ! degrees North
height  = point_coords(3) ! whatever

! Determine the rotated lat & lon of the desired location

rotated_lat = phi2phirot( geo_lat, geo_lon, north_pole_latitude, north_pole_longitude)
rotated_lon = rla2rlarot( geo_lat, geo_lon, north_pole_latitude, north_pole_longitude, polgam)

if (debug > 99 .and. do_output()) then
   write(*,*)
   write(*,*)'The geographic longitude then rotated longitude is ',geo_lon, rotated_lon
   write(*,*)'The geographic  latitude then rotated  latitude is ',geo_lat, rotated_lat
endif

if (    dart_kind == KIND_U_WIND_COMPONENT ) then
   call get_corners(rotated_lon, rotated_lat, nsrlon,  nrlat, srlon,  rlat, ileft, iright, ifrac, jbot, jtop, jfrac, istatus)
elseif ( dart_kind == KIND_V_WIND_COMPONENT ) then
   call get_corners(rotated_lon, rotated_lat,  nrlon, nsrlat,  rlon, srlat, ileft, iright, ifrac, jbot, jtop, jfrac, istatus)
else
   call get_corners(rotated_lon, rotated_lat,  nrlon,  nrlat,  rlon,  rlat, ileft, iright, ifrac, jbot, jtop, jfrac, istatus)
endif

if (istatus /= 0) return   ! and pass on the failed istatus of get_corners()

call get_level_indices(height, dart_kind, ibelow, iabove, levelfrac, istatus)

if (istatus /= 0) return   ! and pass on the failed istatus of get_level_indices()

call horizontal_interpolate(x, ivar, iabove, ileft, iright, ifrac, jbot, jtop, jfrac, value_above, istatus) 
call horizontal_interpolate(x, ivar, ibelow, ileft, iright, ifrac, jbot, jtop, jfrac, value_below, istatus) 

! vertically interpolate the layers to the desired height

if (debug > 99 .and. do_output()) then
   ! to test, use model_mod_check_nml and pick a level that is a known distance
   value_above = 10.0_r8
   value_below = 20.0_r8
   
   interp_val = value_below*levelfrac + value_above*(1.0_r8 - levelfrac)
   
   write(*,*)
   write(*,*)' level above, level, level below ',vcoord%level(iabove), height, vcoord%level(ibelow)
   write(*,*)' distances   ',vcoord%level(iabove) - height,  height - vcoord%level(ibelow)
   write(*,*)' fractions   ',levelfrac, 1.0_r8 - levelfrac
   write(*,*)' layervalues ',value_above, value_below
   write(*,*)' final interpolated value is ',interp_val
   write(*,*)
endif

interp_val = value_below*levelfrac + value_above*(1.0_r8 - levelfrac)
istatus    = 0

return

end subroutine model_interpolate


!------------------------------------------------------------------------
!>
!> Returns a model state vector, x, that is some sort of appropriate
!> initial condition for starting up a long integration of the model.
!> At present, this is only used if the namelist parameter
!> start_from_restart is set to .false. in the program perfect_model_obs.

subroutine init_conditions(x)

real(r8), intent(out) :: x(:)

if ( .not. module_initialized ) call static_init_model

write(string1,*)'Cannot initialize COSMO time via subroutine call.'
write(string2,*)'input.nml:start_from_restart cannot be FALSE'
call error_handler(E_ERR, 'init_conditions', string1, &
           source, revision, revdate, text2=string2)

x = 0.0_r8  ! suppress compiler warnings about unused variables

return

end subroutine init_conditions


!------------------------------------------------------------------------
!> Companion interface to init_conditions. Returns a time that is somehow
!> appropriate for starting up a long integration of the model.
!> At present, this is only used if the namelist parameter
!> start_from_restart is set to .false. in the program perfect_model_obs.


subroutine init_time(time)

type(time_type), intent(out) :: time

if ( .not. module_initialized ) call static_init_model

write(string1,*)'Cannot initialize COSMO time via subroutine call.'
write(string2,*)'input.nml:start_from_restart cannot be FALSE'
call error_handler(E_ERR, 'init_time', string1, &
           source, revision, revdate, text2=string2)

time = set_time(0,0) ! suppress compiler warnings about unused variables

return

end subroutine init_time


!------------------------------------------------------------------------
!> As COSMO can only be advanced as a separate executable,
!> this is a NULL INTERFACE and is a fatal error if invoked.

subroutine adv_1step(x, time)

real(r8),        intent(inout) :: x(:)
type(time_type), intent(in)    :: time

if ( .not. module_initialized ) call static_init_model

if (do_output()) then
    call print_time(time,'NULL interface adv_1step (no advance) DART time is')
    call print_time(time,'NULL interface adv_1step (no advance) DART time is',logfileunit)
endif

write(string1,*) 'Cannot advance COSMO with a subroutine call; async cannot equal 0'
call error_handler(E_ERR,'adv_1step',string1,source,revision,revdate)

return

end subroutine adv_1step


!------------------------------------------------------------------------
!>


subroutine end_model()

deallocate(lon,lat,slonu,slatu,slonv,slatv)
deallocate(vcoord%level1, vcoord%level) 
deallocate(rlon,rlat,srlon,srlat)

return

end subroutine end_model


!------------------------------------------------------------------------
!>


function nc_write_model_atts( ncFileID ) result (ierr)

integer, intent(in)  :: ncFileID      ! netCDF file identifier
integer              :: ierr          ! return value of function

integer  :: nDimensions, nVariables, nAttributes, unlimitedDimID, TimeDimID
integer  :: StateVarDimID
integer  :: MemberDimID
integer  :: LineLenDimID
integer  :: StateVarVarID, StateVarID, VarID, vcoordVarID

integer  ::   lonDimID,  srlonDimID
integer  ::   latDimID,  srlatDimID
integer  :: levelDimID, level1DimID

integer  :: io, ivar, ndims, i


character(len=256)   :: filename

character(len=8)      :: crdate      ! needed by F90 DATE_AND_TIME intrinsic
character(len=10)     :: crtime      ! needed by F90 DATE_AND_TIME intrinsic
character(len=5)      :: crzone      ! needed by F90 DATE_AND_TIME intrinsic
integer, dimension(8) :: values      ! needed by F90 DATE_AND_TIME intrinsic

integer, dimension(NF90_MAX_VAR_DIMS) :: mydimids
character(len=NF90_MAX_NAME) :: varname

if ( .not. module_initialized ) call static_init_model

ierr = -1 ! assume things go poorly

write(filename,*) 'ncFileID', ncFileID

!-------------------------------------------------------------------------------
! make sure ncFileID refers to an open netCDF file,
! and then put into define mode.
!-------------------------------------------------------------------------------

call nc_check(nf90_Inquire(ncFileID,nDimensions,nVariables,nAttributes,unlimitedDimID),&
                                   'nc_write_model_atts', 'inquire '//trim(filename))
call nc_check(nf90_Redef(ncFileID),'nc_write_model_atts',   'redef '//trim(filename))

!-------------------------------------------------------------------------------
! We need the dimension ID for the number of copies/ensemble members, and
! we might as well check to make sure that Time is the Unlimited dimension.
! Our job is create the 'model size' dimension.
!-------------------------------------------------------------------------------

call nc_check(nf90_inq_dimid(ncid=ncFileID, name='NMLlinelen', dimid=LineLenDimID), &
 'nc_write_model_atts','inq_dimid NMLlinelen')
call nc_check(nf90_inq_dimid(ncid=ncFileID, name='copy', dimid=MemberDimID), &
 'nc_write_model_atts', 'copy dimid '//trim(filename))
call nc_check(nf90_inq_dimid(ncid=ncFileID, name='time', dimid=  TimeDimID), &
 'nc_write_model_atts', 'time dimid '//trim(filename))

if ( TimeDimID /= unlimitedDimId ) then
   write(string1,*)'Time Dimension ID ',TimeDimID, &
         ' must equal Unlimited Dimension ID ',unlimitedDimID
   call error_handler(E_ERR,'nc_write_model_atts', string1, source, revision, revdate)
endif

!-------------------------------------------------------------------------------
! Define the model size / state variable dimension / whatever ...
!-------------------------------------------------------------------------------
call nc_check(nf90_def_dim(ncid=ncFileID, name='StateVariable', len=model_size, &
 dimid = StateVarDimID),'nc_write_model_atts', 'state def_dim '//trim(filename))

!-------------------------------------------------------------------------------
! Write Global Attributes
!-------------------------------------------------------------------------------

call DATE_AND_TIME(crdate,crtime,crzone,values)
write(string1,'(''YYYY MM DD HH MM SS = '',i4,5(1x,i2.2))') &
 values(1), values(2), values(3), values(5), values(6), values(7)

call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, 'creation_date' ,string1 ), &
              'nc_write_model_atts', 'creation put '//trim(filename))
call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, 'model_source'  ,source  ), &
              'nc_write_model_atts', 'source put '//trim(filename))
call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, 'model_revision',revision), &
              'nc_write_model_atts', 'revision put '//trim(filename))
call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, 'model_revdate' ,revdate ), &
              'nc_write_model_atts', 'revdate put '//trim(filename))
call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, 'model',  'cosmo' ), &
              'nc_write_model_atts', 'model put '//trim(filename))

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
   call nc_check(nf90_def_var(ncid=ncFileID,name='StateVariable', xtype=nf90_int, &
                 dimids=StateVarDimID, varid=StateVarVarID), 'nc_write_model_atts', &
                 'statevariable def_var '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,StateVarVarID,'long_name','State Variable ID'),&
                 'nc_write_model_atts','statevariable long_name '//trim(filename))
   call nc_check(nf90_put_att(ncFileID, StateVarVarID, 'units','indexical'), &
                 'nc_write_model_atts', 'statevariable units '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,StateVarVarID,'valid_range',(/ 1,model_size /)),&
                 'nc_write_model_atts', 'statevariable valid_range '//trim(filename))

   ! Define the actual (3D) state vector, which gets filled as time goes on ...
   call nc_check(nf90_def_var(ncid=ncFileID, name='state', xtype=nf90_real, &
                 dimids=(/StateVarDimID,MemberDimID,unlimitedDimID/),varid=StateVarID),&
                 'nc_write_model_atts','state def_var '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,StateVarID,'long_name','model state or fcopy'),&
                 'nc_write_model_atts', 'state long_name '//trim(filename))

   ! Leave define mode so we can fill the coordinate variable.
   call nc_check(nf90_enddef(ncfileID),'nc_write_model_atts','state enddef '//trim(filename))

   ! Fill the state variable coordinate variable
   call nc_check(nf90_put_var(ncFileID, StateVarVarID, (/ (i,i=1,model_size) /) ), &
                 'nc_write_model_atts', 'state put_var '//trim(filename))

else

   !----------------------------------------------------------------------------
   ! We need to output the prognostic variables.
   !----------------------------------------------------------------------------
   ! Define the new dimensions IDs
   ! TJH note: we are following the COSMO use of dimensions names. Even though
   ! we are not using the rotated coordinates, they are still using the
   ! rotated coordinate dimensions.
   !----------------------------------------------------------------------------

   io = nf90_def_dim(ncid=ncFileID, name='rlon', len=nrlon, dimid = lonDimID)
   call nc_check(io, 'nc_write_model_atts', 'rlon def_dim '//trim(filename))

   io = nf90_def_dim(ncid=ncFileID, name='rlat', len=nrlat, dimid = latDimID)
   call nc_check(io, 'nc_write_model_atts', 'rlat def_dim '//trim(filename))

   io = nf90_def_dim(ncid=ncFileID, name='srlon', len=nsrlon, dimid = srlonDimID)
   call nc_check(io, 'nc_write_model_atts', 'srlon def_dim '//trim(filename))

   io = nf90_def_dim(ncid=ncFileID, name='srlat', len=nsrlat, dimid = srlatDimID)
   call nc_check(io, 'nc_write_model_atts', 'srlat def_dim '//trim(filename))

   io = nf90_def_dim(ncid=ncFileID, name='level', len=nlevel, dimid = levelDimID)
   call nc_check(io, 'nc_write_model_atts', 'level def_dim '//trim(filename))

   io = nf90_def_dim(ncid=ncFileID, name='level1', len=nlevel1, dimid = level1DimID)
   call nc_check(io, 'nc_write_model_atts', 'level1 def_dim '//trim(filename))

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


   ! staggered U Grid Longitudes
   call nc_check(nf90_def_var(ncFileID,name='slonu', xtype=nf90_real, &
                 dimids=(/ srlonDimID, latDimID /), varid=VarID),&
                 'nc_write_model_atts', 'slonu def_var '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'long_name', 'staggered U-wind longitude'), &
                 'nc_write_model_atts', 'slonu long_name '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'cartesian_axis', 'X'),  &
                 'nc_write_model_atts', 'slonu cartesian_axis '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'units', 'degrees_east'), &
                 'nc_write_model_atts', 'slonu units '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'valid_range', (/ 0.0_r8, 360.0_r8 /)), &
                 'nc_write_model_atts', 'slonu valid_range '//trim(filename))

   ! staggered U Grid Latitudes
   call nc_check(nf90_def_var(ncFileID,name='slatu', xtype=nf90_real, &
                 dimids=(/ lonDimID, latDimID /), varid=VarID),&
                 'nc_write_model_atts', 'slatu def_var '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'long_name', 'staggered U-wind latitude'), &
                 'nc_write_model_atts', 'slatu long_name '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'cartesian_axis', 'Y'),  &
                 'nc_write_model_atts', 'slatu cartesian_axis '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'units', 'degrees_east'), &
                 'nc_write_model_atts', 'slatu units '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'valid_range', (/ -90.0_r8, 90.0_r8 /)), &
                 'nc_write_model_atts', 'slatu valid_range '//trim(filename))

   ! staggered V Grid Longitudes
   call nc_check(nf90_def_var(ncFileID,name='slonv', xtype=nf90_real, &
                 dimids=(/ lonDimID, latDimID /), varid=VarID),&
                 'nc_write_model_atts', 'slonv def_var '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'long_name', 'staggered V-wind longitude'), &
                 'nc_write_model_atts', 'slonv long_name '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'cartesian_axis', 'X'),  &
                 'nc_write_model_atts', 'slonv cartesian_axis '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'units', 'degrees_east'), &
                 'nc_write_model_atts', 'slonv units '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'valid_range', (/ 0.0_r8, 360.0_r8 /)), &
                 'nc_write_model_atts', 'slonv valid_range '//trim(filename))

   ! staggered V Grid Latitudes
   call nc_check(nf90_def_var(ncFileID,name='slatv', xtype=nf90_real, &
                 dimids=(/ lonDimID, latDimID /), varid=VarID),&
                 'nc_write_model_atts', 'slatv def_var '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'long_name', 'staggered V-wind latitude'), &
                 'nc_write_model_atts', 'slatv long_name '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'cartesian_axis', 'Y'),  &
                 'nc_write_model_atts', 'slatv cartesian_axis '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'units', 'degrees_east'), &
                 'nc_write_model_atts', 'slatv units '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'valid_range', (/ -90.0_r8, 90.0_r8 /)), &
                 'nc_write_model_atts', 'slatv valid_range '//trim(filename))

   ! vcoord 

   io = nf90_def_var(ncFileID,name='vcoord', xtype=nf90_real, &
                 dimids=(/ level1DimID /), varid=vcoordVarID)
   call nc_check(io, 'nc_write_model_atts', 'vcoord def_var '//trim(filename))
   io = nf90_put_att(ncFileID, vcoordVarID,'long_name', trim(vcoord%long_name))
   call nc_check(io, 'nc_write_model_atts', 'vcoord long_name '//trim(filename))

   io = nf90_put_att(ncFileID, vcoordVarID, 'units', trim(vcoord%units))
   call nc_check(io, 'nc_write_model_atts', 'vcoord units '//trim(filename))

   io = nf90_put_att(ncFileID, vcoordVarID, 'ivctype', vcoord%ivctype)
   call nc_check(io, 'nc_write_model_atts', 'vcoord ivctype '//trim(filename))

   io = nf90_put_att(ncFileID, vcoordVarID, 'irefatm', vcoord%irefatm)
   call nc_check(io, 'nc_write_model_atts', 'vcoord irefatm '//trim(filename))

   io = nf90_put_att(ncFileID, vcoordVarID, 'p0sl', vcoord%p0sl)
   call nc_check(io, 'nc_write_model_atts', 'vcoord p0sl '//trim(filename))

   io = nf90_put_att(ncFileID, vcoordVarID, 't0sl', vcoord%t0sl)
   call nc_check(io, 'nc_write_model_atts', 'vcoord t0sl '//trim(filename))

   io = nf90_put_att(ncFileID, vcoordVarID, 'dt0lp', vcoord%dt0lp)
   call nc_check(io, 'nc_write_model_atts', 'vcoord dt0lp '//trim(filename))

   io = nf90_put_att(ncFileID, vcoordVarID, 'vcflat', vcoord%vcflat )
   call nc_check(io, 'nc_write_model_atts', 'vcoord vcflat '//trim(filename))

   io = nf90_put_att(ncFileID, vcoordVarID, 'delta_t', vcoord%delta_t)
   call nc_check(io, 'nc_write_model_atts', 'vcoord delta_t '//trim(filename))

   io = nf90_put_att(ncFileID, vcoordVarID, 'h_scal', vcoord%h_scal)
   call nc_check(io, 'nc_write_model_atts', 'vcoord h_scal '//trim(filename))

! TJH    ! Standard Z Levels
! TJH    call nc_check(nf90_def_var(ncFileID,name='LEV', xtype=nf90_real, &
! TJH                  dimids=(/ levDimID /), varid=VarID),&
! TJH                  'nc_write_model_atts', 'LEV def_var '//trim(filename))
! TJH    call nc_check(nf90_put_att(ncFileID,  VarID, 'long_name', 'standard hybrid model levels'), &
! TJH                  'nc_write_model_atts', 'LEV long_name '//trim(filename))
! TJH    call nc_check(nf90_put_att(ncFileID,  VarID, 'cartesian_axis', 'Z'),  &
! TJH                  'nc_write_model_atts', 'LEV cartesian_axis '//trim(filename))
! TJH    call nc_check(nf90_put_att(ncFileID,  VarID, 'units', 'model level'), &
! TJH                  'nc_write_model_atts', 'LEV units '//trim(filename))
! TJH    call nc_check(nf90_put_att(ncFileID,  VarID, 'valid_range', (/ 1._r8,float(nz)+1._r8 /)), &
! TJH                  'nc_write_model_atts', 'LEV valid_range '//trim(filename))
! TJH 
! TJH    ! W-wind Z Levels
! TJH    call nc_check(nf90_def_var(ncFileID,name='WLEV', xtype=nf90_real, &
! TJH                  dimids=(/ wlevDimID /), varid=VarID),&
! TJH                  'nc_write_model_atts', 'WLEV def_var '//trim(filename))
! TJH    call nc_check(nf90_put_att(ncFileID,  VarID, 'long_name', 'standard model levels for W-wind'), &
! TJH                  'nc_write_model_atts', 'WLEV long_name '//trim(filename))
! TJH    call nc_check(nf90_put_att(ncFileID,  VarID, 'cartesian_axis', 'Z'),  &
! TJH                  'nc_write_model_atts', 'WLEV cartesian_axis '//trim(filename))
! TJH    call nc_check(nf90_put_att(ncFileID,  VarID, 'units', 'model level'), &
! TJH                  'nc_write_model_atts', 'WLEV units '//trim(filename))
! TJH    call nc_check(nf90_put_att(ncFileID,  VarID, 'valid_range', (/ 1._r8,float(nz)+1._r8 /)), &
! TJH                  'nc_write_model_atts', 'WLEV valid_range '//trim(filename))

   !----------------------------------------------------------------------------
   ! Create the (empty) Prognostic Variables and the Attributes
   !----------------------------------------------------------------------------

   do ivar=1, nfields

      varname = trim(progvar(ivar)%varname)
      string1 = trim(filename)//' '//trim(varname)

      ! match shape of the variable to the dimension IDs

      call define_var_dims(ivar, ncFileID, MemberDimID, unlimitedDimID, ndims, mydimids)

      ! define the variable and set the attributes

      io = nf90_def_var(ncid=ncFileID, name=trim(varname), xtype=nf90_double, &
                    dimids = mydimids(1:ndims), varid=VarID)
      call nc_check(io, 'nc_write_model_atts', trim(string1)//' def_var' )

      io = nf90_put_att(ncFileID, VarID, 'long_name', trim(progvar(ivar)%long_name))
      call nc_check(io, 'nc_write_model_atts', trim(string1)//' put_att long_name' )

      io = nf90_put_att(ncFileID, VarID, 'DART_kind', trim(progvar(ivar)%kind_string))
      call nc_check(io, 'nc_write_model_atts', trim(string1)//' put_att dart_kind' )

      io = nf90_put_att(ncFileID, VarID, 'units', trim(progvar(ivar)%units))
      call nc_check(io, 'nc_write_model_atts', trim(string1)//' put_att units' )

   enddo

   !----------------------------------------------------------------------------
   ! Leave define mode so we can fill the coordinate variable.
   !----------------------------------------------------------------------------

   io = nf90_enddef(ncfileID)
   call nc_check(io ,'nc_write_model_atts','prognostic enddef '//trim(filename))

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

   io = nf90_inq_varid(ncFileID, 'slonu', VarID)
   call nc_check(io, 'nc_write_model_atts', 'slonu inq_varid '//trim(filename))
   io = nf90_put_var(ncFileID, VarID, slonu )
   call nc_check(io, 'nc_write_model_atts', 'slonu put_var '//trim(filename))

   io = nf90_inq_varid(ncFileID, 'slatu', VarID)
   call nc_check(io, 'nc_write_model_atts', 'slatu inq_varid '//trim(filename))
   io = nf90_put_var(ncFileID, VarID, slatu )
   call nc_check(io, 'nc_write_model_atts', 'slatu put_var '//trim(filename))

   io = nf90_inq_varid(ncFileID, 'slonv', VarID)
   call nc_check(io, 'nc_write_model_atts', 'slonv inq_varid '//trim(filename))
   io = nf90_put_var(ncFileID, VarID, slonv )
   call nc_check(io, 'nc_write_model_atts', 'slonv put_var '//trim(filename))

   io = nf90_inq_varid(ncFileID, 'slatv', VarID)
   call nc_check(io, 'nc_write_model_atts', 'slatv inq_varid '//trim(filename))
   io = nf90_put_var(ncFileID, VarID, slatv )
   call nc_check(io, 'nc_write_model_atts', 'slatv put_var '//trim(filename))

   io = nf90_put_var(ncFileID, vcoordVarID, vcoord%level1)
   call nc_check(io, 'nc_write_model_atts', 'vcoord put_var '//trim(filename))

endif

!-------------------------------------------------------------------------------
! Flush the buffer and leave netCDF file open
!-------------------------------------------------------------------------------
call nc_check(nf90_sync(ncFileID), 'nc_write_model_atts', 'atts sync')

ierr = 0 ! If we got here, things went well.

end function nc_write_model_atts


!------------------------------------------------------------------------
!>
!> Writes the model variables to a netCDF file.
!> All errors are fatal, so the return code is always '0 == normal'.


function nc_write_model_vars( ncFileID, state_vec, copyindex, timeindex ) result (ierr)

integer,  intent(in) :: ncFileID
real(r8), intent(in) :: state_vec(:)
integer,  intent(in) :: copyindex
integer,  intent(in) :: timeindex
integer              :: ierr          ! return value of function


character(len=NF90_MAX_NAME) :: varname
integer ::  dimIDs(NF90_MAX_VAR_DIMS)
integer :: ncstart(NF90_MAX_VAR_DIMS)
integer :: nccount(NF90_MAX_VAR_DIMS)

integer :: io, i, ivar, VarID, ndims, dimlen
integer :: index1, indexN
integer :: TimeDimID, CopyDimID


character(len=256) :: filename

if ( .not. module_initialized ) call static_init_model

ierr = -1 ! assume things go poorly

! we only have a netcdf handle here so we do not know the filename
! or the fortran unit number.  but construct a string with at least
! the netcdf handle, so in case of error we can trace back to see
! which netcdf file is involved.

write(filename,*) 'ncFileID', ncFileID

! make sure ncFileID refers to an open netCDF file,

io = nf90_inq_dimid(ncFileID, 'copy', dimid=CopyDimID)
call nc_check(io, 'nc_write_model_vars', 'inq_dimid copy '//trim(filename))

io = nf90_inq_dimid(ncFileID, 'time', dimid=TimeDimID)
call nc_check(io, 'nc_write_model_vars', 'inq_dimid time '//trim(filename))


if ( output_1D_state_vector ) then

   ! blast out the DART vector as a 1D array.

   io = nf90_inq_varid(ncFileID, 'state', VarID)
   call nc_check(io, 'nc_write_model_vars', 'state inq_varid '//trim(filename))

   io = nf90_put_var(ncFileID,VarID,state_vec,start=(/1,copyindex,timeindex/))
   call nc_check(io, 'nc_write_model_vars', 'state put_var '//trim(filename))

else

   ! write out the DART vector as the native variables and shapes

   do ivar = 1,nfields  ! Very similar to loop in sv_to_restart_file

      varname = trim(progvar(ivar)%varname)
      string2 = trim(filename)//' '//trim(varname)

      ! Ensure netCDF variable is conformable with progvar quantity.
      ! The TIME and Copy dimensions are intentionally not queried
      ! by looping over the dimensions stored in the progvar type.

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

      if ((debug > 99) .and. do_output()) then
         write(*,*)'nc_write_model_vars '//trim(varname)//' start is ',ncstart(1:ndims)
         write(*,*)'nc_write_model_vars '//trim(varname)//' count is ',nccount(1:ndims)
      endif

      ! revelation - no need to reshape the state vector before nf90_put_var

      index1 = progvar(ivar)%index1
      indexN = progvar(ivar)%indexN
      call nc_check(nf90_put_var(ncFileID, VarID, state_vec(index1:indexN), &
                   start = ncstart(1:ndims), count=nccount(1:ndims)), &
                   'nc_write_model_vars', 'put_var '//trim(string2))
   enddo

endif

return

end function nc_write_model_vars


!------------------------------------------------------------------------
!>


  subroutine pert_model_state(state, pert_state, interf_provided)
    !------------------------------------------------------------------
    ! Perturbs a model state for generating initial ensembles.
    ! The perturbed state is returned in pert_state.
    ! A model may choose to provide a NULL INTERFACE by returning
    ! .false. for the interf_provided argument. This indicates to
    ! the filter that if it needs to generate perturbed states, it
    ! may do so by adding a perturbation to each model state
    ! variable independently. The interf_provided argument
    ! should be returned as .true. if the model wants to do its own
    ! perturbing of states.
    !------------------------------------------------------------------
    ! Currently only implemented as rondom perturbations
    !------------------------------------------------------------------

    real(r8), intent(in)  :: state(:)
    real(r8), intent(out) :: pert_state(:)
    logical,  intent(out) :: interf_provided

    real(r8)              :: stddev,mean

    integer               :: ikind,ilevel,i,istart,iend
    logical, save         :: random_seq_init = .false.

  call error_handler(E_ERR,'pert_model_state','routine not written',source,revision,revdate)
    if ( .not. module_initialized ) call static_init_model

    interf_provided = .true.

    ! Initialize my random number sequence (no seed is submitted here!)
    if(.not. random_seq_init) then
      call init_random_seq(random_seq)
      random_seq_init = .true.
    endif

    ! add some uncertainty to every state vector element
    do ikind=1,size(state_vector_vars)
      if (state_vector_vars(ikind)%is_present) then
        do ilevel=1,state_vector_vars(ikind)%nz
          istart=state_vector_vars(ikind)%state_vector_sindex(ilevel)
          iend=istart+(state_vector_vars(ikind)%nx*state_vector_vars(ikind)%ny)-1

          mean=sum(abs(state(istart:iend)))/float(iend-istart+1)
          stddev=sqrt(sum((state(istart:iend)-mean)**2))/float(iend-istart+1)

          do i=istart,iend
            pert_state(i) = random_gaussian(random_seq, state(i),model_perturbation_amplitude*stddev)
          enddo
          if ((ikind==KIND_SPECIFIC_HUMIDITY) .or. &
              (ikind==KIND_CLOUD_LIQUID_WATER) .or. &
              (ikind==KIND_CLOUD_ICE)) then
            where (pert_state(istart:iend)<0.)
              pert_state(istart:iend)=0.
            end where
          endif
        enddo
      endif
    enddo

    return

  end subroutine pert_model_state


!------------------------------------------------------------------------
!>


subroutine ens_mean_for_model(filter_ens_mean)

real(r8), dimension(:), intent(in) :: filter_ens_mean

call error_handler(E_ERR,'ens_mean_for_model','routine not written',source,revision,revdate)

if ( .not. module_initialized ) call static_init_model

allocate(ens_mean(1:model_size))
ens_mean(:) = filter_ens_mean(:)

!  write(string1,*) 'COSMO has no ensemble mean in storage.'
!  call error_handler(E_ERR,'ens_mean_for_model',string1,source,revision,revdate)

return

end subroutine ens_mean_for_model


!------------------------------------------------------------------------
!>


function get_state_time() result (time)
type(time_type) :: time

call error_handler(E_ERR,'get_state_time','routine not written',source,revision,revdate)

if ( .not. module_initialized ) call static_init_model
time=cosmo_fc_time

return

end function get_state_time


!------------------------------------------------------------------------
!>


subroutine get_state_vector(statevector, model_time)

real(r8),        intent(out) :: statevector(1:model_size)
type(time_type), intent(out) :: model_time

integer  :: ivar

if ( .not. module_initialized ) call static_init_model

do ivar=1,nfields
   call read_binary_file(ivar, statevector, model_time)
enddo

return

end subroutine get_state_vector


!------------------------------------------------------------------------
!>


subroutine write_cosmo_file(sv, dart_file, newfile)

real(r8),         intent(in) :: sv(:)      ! the DART posterior
character(len=*), intent(in) :: dart_file  ! where the posterior came from
character(len=*), intent(in) :: newfile    ! the name of the new cosmo restart

integer(i4) :: ipdsbuf(NPDS)     ! pds: product definition section
integer(i4) :: igdsbuf(NGDS)     ! gds: grid definition section
real(r8)    :: rbuf(nrlon*nrlat) ! data to be read

logical :: desired = .false.
integer :: index1, indexN
integer :: iunit_in, iunit_out

integer :: izerr      ! error status
integer :: iz_countl  ! read counter for binary data

if ( .not. module_initialized ) call static_init_model

write(string1,*)'..  The (untouched) input cosmos restart file is "'//trim(cosmo_restart_file)//'"'
write(string2,*)    'The DART posterior file is "'//trim(dart_file)//'"'
write(string3,*)    'The new (posterior) cosmos restart file is "'//trim(newfile)//'"'
call error_handler(E_MSG, 'write_cosmo_file', string1, &
           source, revision, revdate, text2=string2, text3=string3)

! The idea is to read the original input file and check to see if the record data needs
! to be replaced with the DART posterior or not. We can use the progvar(ivar)%slab1, etc.
! information to do this.

! prepare all data from the input binary file
iunit_in  = open_file(cosmo_restart_file, form='unformatted', action='read')
iunit_out = open_file(newfile,            form='unformatted', action='write')

! Can call rw_binary_header with ivar = 1 just to get past the header.
! With the third argument, whatever is read in the first file will
! be written to the output file.

call rw_binary_header(iunit_in, 1, iunit_out)

iz_countl = 0        ! keep track of the 'slab index'
desired   = .false.  ! presume we do not want this slab

COPY_LOOP: do

   read(iunit_in, iostat=izerr) ipdsbuf, igdsbuf, rbuf   ! read a 'slab'
   if (izerr < 0) then
      exit COPY_LOOP
   elseif (izerr > 0) then
      write(string1,*) 'reading input cosmos file around data record ',iz_countl
      call error_handler(E_ERR,'write_cosmo_file', string1, source, revision, revdate)
   endif

   iz_countl = iz_countl + 1

   call get_dart_indices(iz_countl, index1, indexN, desired)

   if (desired) rbuf = sv(index1:indexN)

   write(iunit_out, iostat=izerr) ipdsbuf, igdsbuf, rbuf   ! write a 'slab'
   if (izerr /= 0) then
      write(string1,*) 'writing posterior cosmos file around data record ',iz_countl
      call error_handler(E_ERR,'write_cosmo_file', string1, source, revision, revdate)
   endif

enddo COPY_LOOP

call close_file(iunit_in)
call close_file(iunit_out)

return

end subroutine write_cosmo_file


!------------------------------------------------------------------------
!>


function get_cosmo_filename(filetype)
character(len=*), optional, intent(in) :: filetype
character(len=256) :: get_cosmo_filename

character(len=256) :: lj_filename

lj_filename = adjustl(cosmo_restart_file)

if (present(filetype)) then
   if (trim(filetype) == 'netcdf') then
      lj_filename = adjustl(cosmo_netcdf_file)
   endif
endif

get_cosmo_filename = trim(lj_filename)

end function get_cosmo_filename


!------------------------------------------------------------------------
!>


  subroutine write_state_times(iunit, statetime, advancetime)
    integer,         intent(in) :: iunit
    type(time_type), intent(in) :: statetime, advancetime

    character(len=32) :: timestring
    integer           :: iyear, imonth, iday, ihour, imin, isec
    integer           :: ndays, nhours, nmins, nsecs
    type(time_type)   :: interval

    call error_handler(E_ERR,'write_state_times','routine not written',source,revision,revdate)

    call get_date(statetime, iyear, imonth, iday, ihour, imin, isec)
    write(timestring, "(I4,5(1X,I2))") iyear, imonth, iday, ihour, imin, isec
    write(iunit, "(A)") trim(timestring)

    call get_date(advancetime, iyear, imonth, iday, ihour, imin, isec)
    write(timestring, "(I4,5(1X,I2))") iyear, imonth, iday, ihour, imin, isec
    write(iunit, "(A)") trim(timestring)

    interval = advancetime - statetime
    call get_time(interval, nsecs, ndays)
    nhours = nsecs / (60*60)
    nsecs  = nsecs - (nhours * 60*60)
    nmins  = nsecs / 60
    nsecs  = nsecs - (nmins * 60)

    write(timestring, "(I4,3(1X,I2))") ndays, nhours, nmins, nsecs
    write(iunit, "(A)") trim(timestring)

    return

  end subroutine write_state_times


!------------------------------------------------------------------------
!>


subroutine get_cosmo_grid(filename)

character(len=*), intent(in) :: filename

! Our little test case had these dimensions - only the names are important
!       time = UNLIMITED ; // (1 currently)
!       bnds = 2 ;
!       rlon = 30 ;
!       rlat = 20 ;
!       srlon = 30 ;
!       srlat = 20 ;
!       level = 50 ;
!       level1 = 51 ;
!       soil = 7 ;
!       soil1 = 8 ;

integer :: ncid, io, dimid, VarID
integer :: i, j, k, iii, iunit, iunitm

real(r8) :: rotated_lon, rotated_lat

io = nf90_open(filename, NF90_NOWRITE, ncid)
call nc_check(io, 'get_cosmo_grid','open "'//trim(filename)//'"')

write(string1,*) 'time: '//trim(filename)
io = nf90_inq_dimid(ncid, 'time', dimid)
call nc_check(io, 'get_cosmo_grid','inq_dimid '//trim(string1))
io = nf90_inquire_dimension(ncid, dimid, len=ntime)
call nc_check(io, 'get_cosmo_grid','inquire_dimension '//trim(string1))

write(string1,*) 'bnds: '//trim(filename)
io = nf90_inq_dimid(ncid, 'bnds', dimid)
call nc_check(io, 'get_cosmo_grid','inq_dimid '//trim(string1))
io = nf90_inquire_dimension(ncid, dimid, len=nbnds)
call nc_check(io, 'get_cosmo_grid','inquire_dimension '//trim(string1))

write(string1,*) 'rlon: '//trim(filename)
io = nf90_inq_dimid(ncid, 'rlon', dimid)
call nc_check(io, 'get_cosmo_grid','inq_dimid '//trim(string1))
io = nf90_inquire_dimension(ncid, dimid, len=nrlon)
call nc_check(io, 'get_cosmo_grid','inquire_dimension '//trim(string1))

write(string1,*) 'rlat: '//trim(filename)
io = nf90_inq_dimid(ncid, 'rlat', dimid)
call nc_check(io, 'get_cosmo_grid','inq_dimid '//trim(string1))
io = nf90_inquire_dimension(ncid, dimid, len=nrlat)
call nc_check(io, 'get_cosmo_grid','inquire_dimension '//trim(string1))

write(string1,*) 'srlon: '//trim(filename)
io = nf90_inq_dimid(ncid, 'srlon', dimid)
call nc_check(io, 'get_cosmo_grid','inq_dimid '//trim(string1))
io = nf90_inquire_dimension(ncid, dimid, len=nsrlon)
call nc_check(io, 'get_cosmo_grid','inquire_dimension '//trim(string1))

write(string1,*) 'srlat: '//trim(filename)
io = nf90_inq_dimid(ncid, 'srlat', dimid)
call nc_check(io, 'get_cosmo_grid','inq_dimid '//trim(string1))
io = nf90_inquire_dimension(ncid, dimid, len=nsrlat)
call nc_check(io, 'get_cosmo_grid','inquire_dimension '//trim(string1))

write(string1,*) 'level: '//trim(filename)
io = nf90_inq_dimid(ncid, 'level', dimid)
call nc_check(io, 'get_cosmo_grid','inq_dimid '//trim(string1))
io = nf90_inquire_dimension(ncid, dimid, len=nlevel)
call nc_check(io, 'get_cosmo_grid','inquire_dimension '//trim(string1))

write(string1,*) 'level1: '//trim(filename)
io = nf90_inq_dimid(ncid, 'level1', dimid)
call nc_check(io, 'get_cosmo_grid','inq_dimid '//trim(string1))
io = nf90_inquire_dimension(ncid, dimid, len=nlevel1)
call nc_check(io, 'get_cosmo_grid','inquire_dimension '//trim(string1))

write(string1,*) 'soil: '//trim(filename)
io = nf90_inq_dimid(ncid, 'soil', dimid)
call nc_check(io, 'get_cosmo_grid','inq_dimid '//trim(string1))
io = nf90_inquire_dimension(ncid, dimid, len=nsoil)
call nc_check(io, 'get_cosmo_grid','inquire_dimension '//trim(string1))

write(string1,*) 'soil1: '//trim(filename)
io = nf90_inq_dimid(ncid, 'soil1', dimid)
call nc_check(io, 'get_cosmo_grid','inq_dimid '//trim(string1))
io = nf90_inquire_dimension(ncid, dimid, len=nsoil1)
call nc_check(io, 'get_cosmo_grid','inquire_dimension '//trim(string1))

if ( debug > 99 .and. do_output() ) then

   write(string1,*)'time   has dimension ',ntime
   call error_handler(E_MSG,'get_cosmo_grid',string1)

   write(string1,*)'bnds   has dimension ',nbnds
   call error_handler(E_MSG,'get_cosmo_grid',string1)

   write(string1,*)'rlon   has dimension ',nrlon
   call error_handler(E_MSG,'get_cosmo_grid',string1)

   write(string1,*)'rlat   has dimension ',nrlat
   call error_handler(E_MSG,'get_cosmo_grid',string1)

   write(string1,*)'srlon  has dimension ',nsrlon
   call error_handler(E_MSG,'get_cosmo_grid',string1)

   write(string1,*)'srlat  has dimension ',nsrlat
   call error_handler(E_MSG,'get_cosmo_grid',string1)

   write(string1,*)'level  has dimension ',nlevel
   call error_handler(E_MSG,'get_cosmo_grid',string1)

   write(string1,*)'level1 has dimension ',nlevel1
   call error_handler(E_MSG,'get_cosmo_grid',string1)

   write(string1,*)'soil   has dimension ',nsoil
   call error_handler(E_MSG,'get_cosmo_grid',string1)

   write(string1,*)'soil1  has dimension ',nsoil1
   call error_handler(E_MSG,'get_cosmo_grid',string1)

endif

! Get the geographic grid information

call get_grid_var(ncid,  'lon' ,  nrlon,  nrlat, filename)
call get_grid_var(ncid,  'lat' ,  nrlon,  nrlat, filename)
call get_grid_var(ncid, 'slonu', nsrlon,  nrlat, filename)
call get_grid_var(ncid, 'slatu', nsrlon,  nrlat, filename)
call get_grid_var(ncid, 'slonv',  nrlon, nsrlat, filename)
call get_grid_var(ncid, 'slatv',  nrlon, nsrlat, filename)

where(lon <   0.0_r8) lon = lon + 360.0_r8
where(lat < -90.0_r8) lat = -90.0_r8
where(lat >  90.0_r8) lat =  90.0_r8

call get_vcoord(ncid, filename)

! Get the rotated grid information

call get_grid_var(ncid, 'rlon',   nrlon, filename)
call get_grid_var(ncid, 'rlat',   nrlat, filename)
call get_grid_var(ncid, 'srlon', nsrlon, filename)
call get_grid_var(ncid, 'srlat', nsrlat, filename)

io = nf90_inq_varid(ncid, 'rotated_pole', VarID)
call nc_check(io, 'get_cosmo_grid', 'rotated_pole inq_varid '//trim(filename))

io = nf90_get_att(ncid, VarID, 'grid_north_pole_latitude', north_pole_latitude)
call nc_check(io, 'get_cosmo_grid', 'grid_north_pole_latitude get_att '//trim(filename))

io = nf90_get_att(ncid, VarID, 'grid_north_pole_longitude', north_pole_longitude)
call nc_check(io, 'get_cosmo_grid', 'grid_north_pole_longitude get_att '//trim(filename))

! close up

call nc_check(nf90_close(ncid),'get_cosmo_grid','close "'//trim(filename)//'"' )

! this should only be done once, if at all
if (debug > 99 .and. do_output()) then

   iunit = open_file('exhaustive_grid_table.txt',form='formatted')
   write(iunit,'(''This is for the U grid only.'')')

   iunitm = open_file('exhaustive_grid.m',form='formatted')

   iii = 0

   do k = 1,nlevel
   do j = 1,nrlat
   do i = 1,nrlon
      iii = iii + 1
      
      rotated_lat = phi2phirot( lat(i,j), slonu(i,j), north_pole_latitude, north_pole_longitude)
      rotated_lon = rla2rlarot( lat(i,j), slonu(i,j), north_pole_latitude, north_pole_longitude, 0.0_r8)

      write(iunit,100) i, j, k, iii, slonu(i,j), lat(i,j), vcoord%level(k), rotated_lon, rotated_lat
      write(iunitm,200) i, j, k, iii, slonu(i,j), lat(i,j), vcoord%level(k), rotated_lon, rotated_lat
      
   enddo
   enddo
   enddo

   close(iunit)
   close(iunitm)

 100 format(3(1x,i3),1x,i6,1x,'geo lon,lat,vert',2(1x,f14.9),1x,f12.4,' rlon,rlat ',2(1x,f14.9))
 200 format(3(1x,i3),1x,i6,1x,2(1x,f14.9),1x,f12.4,2(1x,f16.11))

endif

return

end subroutine get_cosmo_grid


!------------------------------------------------------------------------
!> get the 1D grid variables in the netCDF file
!> this does not handle scale, offset, missing_value, _FillValue etc.


subroutine get_1d_grid_var(ncid,varstring,expected_dim1,filename)
integer,          intent(in) :: ncid
character(len=*), intent(in) :: varstring
integer,          intent(in) :: expected_dim1
character(len=*), intent(in) :: filename

integer :: dimIDs(NF90_MAX_VAR_DIMS)
integer :: dimlens
integer :: io, VarID, ndims
integer :: i

write(string3,*)trim(varstring)//' '//trim(filename)

call nc_check(nf90_inq_varid(ncid, trim(varstring), VarID), &
         'get_1d_grid_var', 'inq_varid '//trim(string3))

call nc_check(nf90_inquire_variable(ncid, VarID, dimids=dimIDs, ndims=ndims), &
         'get_1d_grid_var', 'inquire_variable '//trim(string3))

!> @TODO more informative error message
if (ndims /= 1) then
   call error_handler(E_ERR,'get_1d_grid_var','wrong shape for '//string3, &
              source, revision, revdate)
endif

write(string1,'(''inquire dimension'',i2,A)') i,trim(string3)

call nc_check(nf90_inquire_dimension(ncid, dimIDs(1), len=dimlens), &
                                    'get_1d_grid_var', string1)

! Check that the variable actual size matches the expected size

if (dimlens .ne. expected_dim1 ) then
   write(string1,*)'expected dimension 1 to be ',expected_dim1, 'was', dimlens
   call error_handler(E_ERR, 'get_1d_grid_var', string1, &
              source, revision, revdate, text2=string2)
endif

! finally allocate and fill the desired variable

select case (trim(varstring))
   case ("rlon")
      allocate( rlon(dimlens) )
      io = nf90_get_var(ncid, VarID, rlon)
      call nc_check(io, 'get_1d_grid_var', 'get_var '//trim(string3))

   case ("rlat")
      allocate( rlat(dimlens) )
      io = nf90_get_var(ncid, VarID, rlat)
      call nc_check(io, 'get_1d_grid_var', 'get_var '//trim(string3))

   case ("srlon")
      allocate( srlon(dimlens) )
      io = nf90_get_var(ncid, VarID, srlon)
      call nc_check(io, 'get_1d_grid_var', 'get_var '//trim(string3))

   case ("srlat")
      allocate( srlat(dimlens) )
      io = nf90_get_var(ncid, VarID, srlat)
      call nc_check(io, 'get_1d_grid_var', 'get_var '//trim(string3))

   case default
      write(string1,*)'unsupported grid variable '
      call error_handler(E_ERR,'get_1d_grid_var', string1, &
                 source, revision, revdate, text2=string3)

end select

return

end subroutine get_1d_grid_var

!------------------------------------------------------------------------
!> get the 2D grid variables in the netCDF file
!> this does not handle scale, offset, missing_value, _FillValue etc.


subroutine get_2d_grid_var(ncid,varstring,expected_dim1,expected_dim2,filename)
integer,          intent(in) :: ncid
character(len=*), intent(in) :: varstring
integer,          intent(in) :: expected_dim1
integer,          intent(in) :: expected_dim2
character(len=*), intent(in) :: filename

integer, dimension(NF90_MAX_VAR_DIMS) :: dimIDs
integer, dimension(NF90_MAX_VAR_DIMS) :: dimlens
integer :: io, VarID, ndims
integer :: i

write(string3,*)trim(varstring)//' '//trim(filename)

call nc_check(nf90_inq_varid(ncid, trim(varstring), VarID), &
         'get_2d_grid_var', 'inq_varid '//trim(string3))

call nc_check(nf90_inquire_variable(ncid, VarID, dimids=dimIDs, ndims=ndims), &
         'get_2d_grid_var', 'inquire_variable '//trim(string3))

!> @TODO more informative error message
if (ndims /= 2) then
   call error_handler(E_ERR,'get_2d_grid_var','wrong shape for '//string3, &
              source, revision, revdate)
endif

DimensionLoop : do i = 1,ndims

   write(string1,'(''inquire dimension'',i2,A)') i,trim(string3)

   call nc_check(nf90_inquire_dimension(ncid, dimIDs(i), len=dimlens(i)), &
                                       'get_2d_grid_var', string1)

enddo DimensionLoop

! Check that the variable actual sizes match the expected sizes,

if (dimlens(1) .ne. expected_dim1 .or. &
    dimlens(2) .ne. expected_dim2 ) then
   write(string1,*)'expected dimension 1 to be ',expected_dim1, 'was', dimlens(1)
   write(string2,*)'expected dimension 2 to be ',expected_dim2, 'was', dimlens(2)
   call error_handler(E_ERR, 'get_2d_grid_var', string1, &
              source, revision, revdate, text2=string2)
endif

select case (trim(varstring))
   case ("lon")
      allocate( lon(dimlens(1), dimlens(2)) )
      io = nf90_get_var(ncid, VarID, lon)
      call nc_check(io, 'get_2d_grid_var', 'get_var '//trim(string3))

   case ("lat")
      allocate( lat(dimlens(1), dimlens(2)) )
      io = nf90_get_var(ncid, VarID, lat)
      call nc_check(io, 'get_2d_grid_var', 'get_var '//trim(string3))

   case ("slonu")
      allocate( slonu(dimlens(1), dimlens(2)) )
      io = nf90_get_var(ncid, VarID, slonu)
      call nc_check(io, 'get_2d_grid_var', 'get_var '//trim(string3))

   case ("slatu")
      allocate( slatu(dimlens(1), dimlens(2)) )
      io = nf90_get_var(ncid, VarID, slatu)
      call nc_check(io, 'get_2d_grid_var', 'get_var '//trim(string3))

   case ("slonv")
      allocate( slonv(dimlens(1), dimlens(2)) )
      io = nf90_get_var(ncid, VarID, slonv)
      call nc_check(io, 'get_2d_grid_var', 'get_var '//trim(string3))

   case ("slatv")
      allocate( slatv(dimlens(1), dimlens(2)) )
      io = nf90_get_var(ncid, VarID, slatv)
      call nc_check(io, 'get_2d_grid_var', 'get_var '//trim(string3))

   case default
      write(string1,*)'unsupported grid variable '
      call error_handler(E_ERR,'get_2d_grid_var', string1, &
                 source, revision, revdate, text2=string3)

end select

return

end subroutine get_2d_grid_var


!------------------------------------------------------------------------
!>


subroutine get_vcoord(ncid, filename)
integer,          intent(in) :: ncid
character(len=*), intent(in) :: filename

integer, dimension(NF90_MAX_VAR_DIMS) :: dimIDs
integer, dimension(NF90_MAX_VAR_DIMS) :: dimlens
integer :: io, VarID, ndims, i

write(string3,*)' vcoord from '//trim(filename)

call nc_check(nf90_inq_varid(ncid, 'vcoord', VarID), &
         'get_vcoord', 'inq_varid '//trim(string3))

call nc_check(nf90_inquire_variable(ncid, VarID, dimids=dimIDs, ndims=ndims), &
         'get_vcoord', 'inquire_variable '//trim(string3))

!> @TODO more informative error message
if (ndims /= 1) then
   call error_handler(E_ERR,'get_vcoord','wrong shape for '//string3, &
              source, revision, revdate)
endif

io = nf90_inquire_dimension(ncid, dimIDs(1), len=dimlens(1))
call nc_check(io, 'get_vcoord', 'inquire_dimension '//trim(string3))

! Check that the variable actual sizes match the expected sizes,

if (dimlens(1) .ne. nlevel1 ) then
   write(string1,*)'expected dimension to be ',nlevel1, 'was', dimlens(1)
   call error_handler(E_ERR, 'get_vcoord', string1, source, revision, revdate)
endif

vcoord%nlevel1 = nlevel1
vcoord%nlevel  = nlevel 

allocate( vcoord%level1(nlevel1) )
allocate( vcoord%level( nlevel ) )

io = nf90_get_var(ncid, VarID, vcoord%level1 )
call nc_check(io, 'get_vcoord', 'get_var '//trim(string3))

do i = 1,nlevel
   vcoord%level(i) = (vcoord%level1(i) + vcoord%level1(i+1))/2.0_r8
enddo

! write(*,*)'debug vcoord%level1 ',vcoord%level1
! write(*,*)'debug vcoord%level ',vcoord%level

io = nf90_get_att(ncid, VarID, 'long_name', vcoord%long_name)
call nc_check(io, 'get_vcoord', 'get_att long_name '//trim(string3))

io = nf90_get_att(ncid, VarID, 'units', vcoord%units)
call nc_check(io, 'get_vcoord', 'get_att units '//trim(string3))

io = nf90_get_att(ncid, VarID, 'ivctype', vcoord%ivctype)
call nc_check(io, 'get_vcoord', 'get_att ivctype '//trim(string3))

io = nf90_get_att(ncid, VarID, 'irefatm', vcoord%irefatm)
call nc_check(io, 'get_vcoord', 'get_att irefatm '//trim(string3))

io = nf90_get_att(ncid, VarID, 'p0sl', vcoord%p0sl)
call nc_check(io, 'get_vcoord', 'get_att p0sl '//trim(string3))

io = nf90_get_att(ncid, VarID, 't0sl', vcoord%t0sl)
call nc_check(io, 'get_vcoord', 'get_att t0sl '//trim(string3))

io = nf90_get_att(ncid, VarID, 'dt0lp', vcoord%dt0lp)
call nc_check(io, 'get_vcoord', 'get_att dt0lp '//trim(string3))

io = nf90_get_att(ncid, VarID, 'vcflat', vcoord%vcflat)
call nc_check(io, 'get_vcoord', 'get_att vcflat '//trim(string3))

io = nf90_get_att(ncid, VarID, 'delta_t', vcoord%delta_t)
call nc_check(io, 'get_vcoord', 'get_att delta_t '//trim(string3))

io = nf90_get_att(ncid, VarID, 'h_scal', vcoord%h_scal)
call nc_check(io, 'get_vcoord', 'get_att h_scal '//trim(string3))

! Print a summary
if (debug > 99 .and. do_output()) then
   write(logfileunit,*)
   write(logfileunit,*)'vcoord long_name: ',trim(vcoord%long_name)
   write(logfileunit,*)'vcoord     units: ',vcoord%units
   write(logfileunit,*)'vcoord   ivctype: ',vcoord%ivctype
   write(logfileunit,*)'vcoord   irefatm: ',vcoord%irefatm
   write(logfileunit,*)'vcoord      p0sl: ',vcoord%p0sl
   write(logfileunit,*)'vcoord      t0sl: ',vcoord%t0sl
   write(logfileunit,*)'vcoord     dt0lp: ',vcoord%dt0lp
   write(logfileunit,*)'vcoord    vcflat: ',vcoord%vcflat
   write(logfileunit,*)'vcoord   delta_t: ',vcoord%delta_t
   write(logfileunit,*)'vcoord    h_scal: ',vcoord%h_scal

   write(*,*)
   write(*,*)'vcoord long_name: ',trim(vcoord%long_name)
   write(*,*)'vcoord     units: ',vcoord%units
   write(*,*)'vcoord   ivctype: ',vcoord%ivctype
   write(*,*)'vcoord   irefatm: ',vcoord%irefatm
   write(*,*)'vcoord      p0sl: ',vcoord%p0sl
   write(*,*)'vcoord      t0sl: ',vcoord%t0sl
   write(*,*)'vcoord     dt0lp: ',vcoord%dt0lp
   write(*,*)'vcoord    vcflat: ',vcoord%vcflat
   write(*,*)'vcoord   delta_t: ',vcoord%delta_t
   write(*,*)'vcoord    h_scal: ',vcoord%h_scal
endif

end subroutine get_vcoord


!------------------------------------------------------------------------
!>  This routine checks the user input against the variables available in the
!>  input netcdf file to see if it is possible to construct the DART state vector
!>  specified by the input.nml:model_nml:clm_variables  variable.
!>  Each variable must have 6 entries.
!>  1: GRIB table version number
!>  2: variable name
!>  3: DART KIND
!>  4: minimum value - as a character string - if none, use 'NA'
!>  5: maximum value - as a character string - if none, use 'NA'
!>  6: does the variable get updated in the restart file or not ...
!>     all variables will be updated INTERNALLY IN DART
!>     'UPDATE'       => update the variable in the restart file
!>     'NO_COPY_BACK' => do not copy the variable back to the restart file


subroutine parse_variable_table( state_variables, ngood, table )

character(len=*), dimension(:),   intent(in)  :: state_variables
integer,                          intent(out) :: ngood
character(len=*), dimension(:,:), intent(out) :: table

integer  :: nrows, ncols, ivar, dart_kind, ios, idummy
real(r8) :: minvalue, maxvalue

character(len=NF90_MAX_NAME) :: gribtableversion
character(len=NF90_MAX_NAME) :: gribvar
character(len=NF90_MAX_NAME) :: varname
character(len=NF90_MAX_NAME) :: dartstr
character(len=NF90_MAX_NAME) :: minvalstring
character(len=NF90_MAX_NAME) :: maxvalstring
character(len=NF90_MAX_NAME) :: state_or_aux

nrows = size(table,1)
ncols = size(table,2)

! This loop just repackages the 1D array of values into a 2D array.
! We can do some miniminal checking along the way.
! Determining which file to check is going to be more complicated.

ngood = 0
MyLoop : do ivar = 1, nrows

   gribtableversion = trim(state_variables(ncols*ivar - 6))
   gribvar          = trim(state_variables(ncols*ivar - 5))
   varname          = trim(state_variables(ncols*ivar - 4))
   dartstr          = trim(state_variables(ncols*ivar - 3))
   minvalstring     = trim(state_variables(ncols*ivar - 2))
   maxvalstring     = trim(state_variables(ncols*ivar - 1))
   state_or_aux     = trim(state_variables(ncols*ivar    ))

   call to_upper(state_or_aux)

   table(ivar,VT_GRIBVERSIONINDX) = trim(gribtableversion)
   table(ivar,VT_GRIBVARINDX)     = trim(gribvar)
   table(ivar,VT_VARNAMEINDX)     = trim(varname)
   table(ivar,VT_KINDINDX)        = trim(dartstr)
   table(ivar,VT_MINVALINDX)      = trim(minvalstring)
   table(ivar,VT_MAXVALINDX)      = trim(maxvalstring)
   table(ivar,VT_STATEINDX)       = trim(state_or_aux)

   ! If the first element is empty, we have found the end of the list.
   if ( table(ivar,1) == ' ' ) exit MyLoop

   ! Any other condition is an error.
   if ( any(table(ivar,:) == ' ') ) then
      string1 = 'input.nml &model_nml:variables not fully specified'
      string2 = 'must be 7 entries per variable. Last known variable name is'
      string3 = '['//trim(table(ivar,1))//'] ... (without the [], naturally)'
      call error_handler(E_ERR, 'parse_variable_table', string1, &
         source, revision, revdate, text2=string2, text3=string3)
   endif

   ! Make sure DART kind is valid
   dart_kind = get_raw_obs_kind_index(dartstr)

   if( dart_kind < 0 ) then
      write(string1,'(''there is no obs_kind <'',a,''> in obs_kind_mod.f90'')') trim(dartstr)
      call error_handler(E_ERR,'parse_variable_table',string1,source,revision,revdate)
   endif

   ! Now that we have the variable table, fill what we can in the progvar structure

   progvar(ivar)%long_name       = 'unknown'
   progvar(ivar)%units           = 'unknown'
   progvar(ivar)%tableID         = 0
   progvar(ivar)%variableID      = 0
   progvar(ivar)%varname         = trim(table(ivar,VT_VARNAMEINDX))
   progvar(ivar)%kind_string     = trim(table(ivar,VT_KINDINDX))
   progvar(ivar)%dart_kind       = dart_kind
   progvar(ivar)%numEW           = 0
   progvar(ivar)%numNS           = 0
   progvar(ivar)%numZ            = 0
   progvar(ivar)%dimnames        = 'unknown'
   progvar(ivar)%dimlens         = 0
   progvar(ivar)%numdims         = 0
   progvar(ivar)%rangeRestricted = BOUNDED_NONE
   progvar(ivar)%minvalue        = MISSING_R8
   progvar(ivar)%maxvalue        = MISSING_R8
   progvar(ivar)%update          = .false.

   read(table(ivar,VT_GRIBVERSIONINDX),*,iostat=ios) idummy
   if (ios == 0) progvar(ivar)%tableID = idummy

   read(table(ivar,VT_GRIBVARINDX),*,iostat=ios) idummy
   if (ios == 0) progvar(ivar)%variableID = idummy

   read(table(ivar,VT_MINVALINDX),*,iostat=ios) minvalue
   if (ios == 0) progvar(ivar)%minvalue = minvalue

   read(table(ivar,VT_MAXVALINDX),*,iostat=ios) maxvalue
   if (ios == 0) progvar(ivar)%maxvalue = maxvalue

   if (table(ivar,VT_STATEINDX) == 'UPDATE') progvar(ivar)%update = .true.

   ! rangeRestricted == BOUNDED_NONE  == 0 ... unlimited range
   ! rangeRestricted == BOUNDED_BELOW == 1 ... minimum, but no maximum
   ! rangeRestricted == BOUNDED_ABOVE == 2 ... maximum, but no minimum
   ! rangeRestricted == BOUNDED_BOTH  == 3 ... minimum and maximum

   if (   (progvar(ivar)%minvalue /= MISSING_R8) .and. &
          (progvar(ivar)%maxvalue /= MISSING_R8) ) then
      progvar(ivar)%rangeRestricted = BOUNDED_BOTH

   elseif (progvar(ivar)%maxvalue /= MISSING_R8) then
      progvar(ivar)%rangeRestricted = BOUNDED_ABOVE

   elseif (progvar(ivar)%minvalue /= MISSING_R8) then
      progvar(ivar)%rangeRestricted = BOUNDED_BELOW

   else
      progvar(ivar)%rangeRestricted = BOUNDED_NONE

   endif

   ! Check to make sure min is less than max if both are specified.

   if ( progvar(ivar)%rangeRestricted == BOUNDED_BOTH ) then
      if (maxvalue < minvalue) then
         write(string1,*)'&model_nml state_variable input error for ',trim(progvar(ivar)%varname)
         write(string2,*)'minimum value (',minvalue,') must be less than '
         write(string3,*)'maximum value (',maxvalue,')'
         call error_handler(E_ERR,'parse_variable_table',string1, &
            source,revision,revdate,text2=trim(string2),text3=trim(string3))
      endif
   endif

   if (debug > 99 .and. do_output()) then
      write(logfileunit,*) &
         'variable ',ivar,' is ',trim(table(ivar,1)), ' ', trim(table(ivar,2)),' ', &
                                 trim(table(ivar,3)), ' ', trim(table(ivar,4)),' ', &
                                 trim(table(ivar,5)), ' ', trim(table(ivar,6)),' ', trim(table(ivar,7))
      write(     *     ,*) &
         'variable ',ivar,' is ',trim(table(ivar,1)), ' ', trim(table(ivar,2)),' ', &
                                 trim(table(ivar,3)), ' ', trim(table(ivar,4)),' ', &
                                 trim(table(ivar,5)), ' ', trim(table(ivar,6)),' ', trim(table(ivar,7))
   endif

   ngood = ngood + 1

enddo MyLoop

if (ngood == nrows) then
   string1 = 'WARNING: There is a possibility you need to increase ''max_state_variables'''
   write(string2,'(''WARNING: you have specified at least '',i4,'' perhaps more.'')')ngood
   call error_handler(E_MSG,'parse_variable_table',string1,text2=string2)
endif

return

end subroutine parse_variable_table



!------------------------------------------------------------------------
!> Read the binary restart file and record the locations and shapes of the
!> variables of interest. Since we don't know the order of the desired 
!> variables in the binary file, each file must be a rewind. We can assume
!> that the variables are stored contiguously so that if we encounter a 
!> new variable type we can stop looking.

subroutine read_binary_file(dartid, statevector, model_time)

integer,                   intent(in)  :: dartid
real(r8),        optional, intent(out) :: statevector(:)
type(time_type), optional, intent(out) :: model_time

! Table to decode the record contents of ipdsbuf
!> @TODO no seconds?
integer, parameter :: indx_gribver   =   2, &
                      indx_var       =   7, &
                      indx_zlevtyp   =   8, &
                      indx_zlevtop   =   9, &
                      indx_zlevbot   =  10

type(time_type) :: slab_time

integer(i4) :: ipdsbuf(NPDS)     ! pds: product definition section
integer(i4) :: igdsbuf(NGDS)     ! gds: grid definition section
real(r8)    :: rbuf(nrlon*nrlat) ! data to be read

integer :: izerr      ! error status
integer :: iz_countl  ! read counter for binary data
integer :: ivar       ! variable reference number based on iver
integer :: iver       ! version number of GRIB1 indicator table
integer :: ilevtyp    ! type of vertical coordinate system

integer :: fid, nx, ny, ilev, ilevp1

real(r8) :: lat1, latN, lon1, lonN

logical :: desired
logical :: filling

!-----------------------------------------------------------------------

if (present(statevector) .and. model_size == 0) then
   string1 = 'model_size unknown, cannot fill state vector yet.'
   string2 = 'must call read_binary_file once to define model_size.'
   call error_handler(E_ERR,'read_binary_file',string1, &
              source, revision, revdate, text2=string2)
endif

if (present(statevector)) then
   filling = .true.
else
   filling = .false.
endif

!-----------------------------------------------------------------------

fid = open_file(cosmo_restart_file, form='unformatted', action='read')

if (debug > 99 .and. do_output()) then
   write(string1,*)'searching for '//trim(progvar(dartid)%varname), &
                   progvar(dartid)%tableID, progvar(dartid)%variableID
   call error_handler(E_MSG,'read_binary_file', string1, source, revision, revdate)
endif

call rw_binary_header(fid, dartid)

!------------------------------------------------------------------------------
!Section 3: READ ALL RECORDS
!------------------------------------------------------------------------------

iz_countl = 0        ! keep track of the 'slab index'
desired   = .false.  ! presume we do not want this slab

read_loop: DO

   read(fid, iostat=izerr) ipdsbuf, igdsbuf, rbuf   ! read a 'slab'

   if (izerr < 0) then
     exit read_loop
   elseif (izerr > 0) then
      write(string1,*) 'ERROR READING RESTART FILE around data record ',iz_countl
      call error_handler(E_ERR,'read_binary_file', string1, source, revision, revdate)
   endif

   iz_countl = iz_countl + 1

   ! Determine table ID, variable ID, the vertical coordinate system, etc.

   iver    = ipdsbuf(indx_gribver)
   ivar    = ipdsbuf(indx_var)
   ilevtyp = ipdsbuf(indx_zlevtyp)
   ilev    = ipdsbuf(indx_zlevtop)
   ilevp1  = ipdsbuf(indx_zlevbot)

   call decode_location(igdsbuf, nx, ny, lat1, lon1, latN, lonN)

   ! Since Fortran binary reads have no way to know if they fail or not
   ! we are going to compare our input with the sizes encoded in the file.
   ! nrlon and nrlat come from the netCDF file. They are our expectation.

   if (nx .ne. nrlon .or. ny .ne. nrlat) then
      write(string1,*) 'reading variable ',trim(progvar(dartid)%varname), ' from '//trim(cosmo_restart_file)
      write(string2,*) '(file) nx /= nrlon ',nx,nrlon, ' or '
      write(string3,*) '(file) ny /= nrlat ',ny,nrlat
      call error_handler(E_ERR,'read_binary_file', string1, &
                 source, revision, revdate, text2=string2, text3=string3)
   endif

   ! Check to see if the slab/variable/etc. is one we want
   ! if so - count up the slabs, add to the variable size and keep moving
   !>@ TODO FIXME the special cases will be difficult. 
   !> What if we want KIND_SURFACE_TEMPERATURE as well as KIND_TEMPERATURE ...
   !> can we key on DART KIND instead of ilevtyp 

   if (progvar(dartid)%tableID == iver .and. progvar(dartid)%variableID == ivar ) then

      ! special case - specific humidity occurs twice - one with
      ! ilevtyp = 1 (surface) and again with ilevtyp = 110 - a 3D variable
      if (progvar(dartid)%tableID == 2 .and. progvar(dartid)%variableID == 51 .and.  ilevtyp < 109 ) then
         cycle read_loop
      endif

      ! special case - temperature occurs twice - one with
      ! ilevtyp = 1 (surface) and again with ilevtyp = 110 - a 3D variable
      if (progvar(dartid)%tableID == 2 .and. progvar(dartid)%variableID == 11 .and.  ilevtyp < 109 ) then
         cycle read_loop
      endif

      if (.not. desired) then
         ! This is the first slab of the variable we want
         progvar(dartid)%slab1 = iz_countl
         call decode_time(ipdsbuf,slab_time)
         if (present(model_time)) model_time = slab_time
      endif

      desired = .true.

      if (debug > 99 .and. do_output()) &
      write (*,'(A,2x,8(1x,i4),2(1xf15.8))')'wanted slab ',iz_countl, &
        iver, ivar, nx, ny, ilev, ilevp1, ilevtyp, MINVAL(rbuf), MAXVAL(rbuf)

      if ( filling ) then
         call insert_slab_in_state(dartid, iz_countl, nx, ny, rbuf, statevector)
      else
         call record_sizes(dartid, nx, ny, iz_countl, ilevtyp, ilev, ilevp1)
      endif

   elseif (desired) then
      ! If we wanted the last slab but the new slab is not of the same variable - we are done.
      exit read_loop
   endif

enddo read_loop

close(fid, iostat=izerr)
if (izerr /= 0) then
   write(string1,*) 'closing '//trim(cosmo_restart_file)//' while looking for ',trim(progvar(dartid)%varname)
   call error_handler(E_ERR,'read_binary_file', string1, source, revision, revdate)
endif

if (.not. desired) then
   write(string1,*) 'never found '//trim(progvar(dartid)%varname)//' in '//trim(cosmo_restart_file)
   call error_handler(E_ERR,'read_binary_file', string1, source, revision, revdate)
endif

return

end subroutine read_binary_file


!------------------------------------------------------------------------
!> determine where to unwrap each of the variables into a DART 1D vector


subroutine set_variable_layout()

integer :: ivar

model_size = 0

do ivar = 1,nfields
   progvar(ivar)%index1 = model_size + 1
   progvar(ivar)%indexN = model_size + progvar(ivar)%varsize
   model_size = progvar(ivar)%indexN
enddo

end subroutine set_variable_layout


!------------------------------------------------------------------------
!> utility to report on the metadata for the composition of the DART state 


subroutine progvar_summary

integer :: ivar,i

! only have 1 task write the report
if ( .not. do_output() ) return

do ivar = 1,nfields
   write(*,*)
   write(*,*)'variable ',ivar,' is ',trim(progvar(ivar)%varname)
   write(*,*)'   long_name       ',  trim(progvar(ivar)%long_name)
   write(*,*)'   units           ',  trim(progvar(ivar)%units)
   write(*,*)'   tableID         ',       progvar(ivar)%tableID
   write(*,*)'   variableID      ',       progvar(ivar)%variableID
   write(*,*)'   levtypeID       ',       progvar(ivar)%levtypeID
   write(*,*)'   numdims         ',       progvar(ivar)%numdims
   write(*,*)'   numEW           ',       progvar(ivar)%numEW
   write(*,*)'   numNS           ',       progvar(ivar)%numNS
   write(*,*)'   numZ            ',       progvar(ivar)%numZ
   write(*,*)'   varsize         ',       progvar(ivar)%varsize
   write(*,*)'   slab1           ',       progvar(ivar)%slab1
   write(*,*)'   slabN           ',       progvar(ivar)%slabN
   write(*,*)'   index1          ',       progvar(ivar)%index1
   write(*,*)'   indexN          ',       progvar(ivar)%indexN
   write(*,*)'   dart_kind       ',       progvar(ivar)%dart_kind
   write(*,*)'   rangeRestricted ',       progvar(ivar)%rangeRestricted
   write(*,*)'   minvalue        ',       progvar(ivar)%minvalue
   write(*,*)'   maxvalue        ',       progvar(ivar)%maxvalue
   write(*,*)'   kind_string     ',  trim(progvar(ivar)%kind_string)
   write(*,*)'   update          ',       progvar(ivar)%update
   do i = 1,progvar(ivar)%numdims
      write(*,*)'   dim ',i, 'has length = ', progvar(ivar)%dimlens(i), &
                ' and is ', trim(progvar(ivar)%dimnames(i)) 
   enddo
   write(*,*)

   write(logfileunit,*)
   write(logfileunit,*)'variable ',ivar,' is ',trim(progvar(ivar)%varname)
   write(logfileunit,*)'   long_name       ',  trim(progvar(ivar)%long_name)
   write(logfileunit,*)'   units           ',  trim(progvar(ivar)%units)
   write(logfileunit,*)'   tableID         ',       progvar(ivar)%tableID
   write(logfileunit,*)'   variableID      ',       progvar(ivar)%variableID
   write(logfileunit,*)'   levtypeID       ',       progvar(ivar)%levtypeID
   write(logfileunit,*)'   numdims         ',       progvar(ivar)%numdims
   write(logfileunit,*)'   numEW           ',       progvar(ivar)%numEW
   write(logfileunit,*)'   numNS           ',       progvar(ivar)%numNS
   write(logfileunit,*)'   numZ            ',       progvar(ivar)%numZ
   write(logfileunit,*)'   varsize         ',       progvar(ivar)%varsize
   write(logfileunit,*)'   slab1           ',       progvar(ivar)%slab1
   write(logfileunit,*)'   slabN           ',       progvar(ivar)%slabN
   write(logfileunit,*)'   index1          ',       progvar(ivar)%index1
   write(logfileunit,*)'   indexN          ',       progvar(ivar)%indexN
   write(logfileunit,*)'   dart_kind       ',       progvar(ivar)%dart_kind
   write(logfileunit,*)'   rangeRestricted ',       progvar(ivar)%rangeRestricted
   write(logfileunit,*)'   minvalue        ',       progvar(ivar)%minvalue
   write(logfileunit,*)'   maxvalue        ',       progvar(ivar)%maxvalue
   write(logfileunit,*)'   kind_string     ',  trim(progvar(ivar)%kind_string)
   write(logfileunit,*)'   update          ',       progvar(ivar)%update
   write(logfileunit,*)
enddo

end subroutine progvar_summary


!------------------------------------------------------------------------
!>


subroutine decode_time(ipdsbuf, model_time)

integer,         intent(in)  :: ipdsbuf(:)
type(time_type), intent(out) :: model_time

integer  :: icc, iyy, imm, idd, ihh, imin, iccyy
integer  :: istartstep, iendstep, nztri
type(time_type) :: base_time
type(time_type) :: run_length

! Table to decode the record contents of ipdsbuf
!> @TODO no seconds?
integer, parameter :: indx_year      =  11, &
                      indx_mm        =  12, &
                      indx_dd        =  13, &
                      indx_hh        =  14, &
                      indx_min       =  15, &
                      indx_startstep =  17, &
                      indx_endstep   =  18, &
                      indx_nztri     =  19, &
                      indx_cc        =  22

icc   = ipdsbuf(indx_cc)-1
iyy   = ipdsbuf(indx_year)
iccyy = iyy + icc*100
imm   = ipdsbuf(indx_mm)
idd   = ipdsbuf(indx_dd)
ihh   = ipdsbuf(indx_hh)
imin  = ipdsbuf(indx_min)

base_time = set_date(iccyy,imm,idd,ihh,imin,0)

istartstep = ipdsbuf(indx_startstep) ! number of hours into the run
iendstep   = ipdsbuf(indx_endstep)   ! not needed by DART, I think
nztri      = ipdsbuf(indx_nztri)

!> @TODO confirm that if nztri /= 0, we do nothing
if (nztri == 0) then
   run_length = set_time(istartstep*60*60,0)
else
   ! In the absence of more information, we are going to error out.
   write(string1,*)'cannot decode model time. expected nztri = 0, it was ',nztri
   call error_handler(E_ERR, 'decode_time', string1, source, revision, revdate)
endif

model_time = base_time + run_length

if (debug > 99 .and. do_output()) then
   call print_time(model_time,'model time is ')
   call print_date(model_time,'model date is ')

!  write(*,'(A,6(1x,i2),3(1x,i4))') 'time as read = ', &
!         icc,iyy,imm,idd,ihh,imin, nztri,istartstep,iendstep
endif

end subroutine decode_time


!------------------------------------------------------------------------
!>


subroutine rw_binary_header(rfid, ivar, wfid)

! some of these might be useful, but not at the moment.
! just need to skip the header to get to the records.

integer,           intent(in) :: rfid    ! The file handle to read from
integer,           intent(in) :: ivar
integer, optional, intent(in) :: wfid    ! the file handle to write to

integer, parameter :: KHMAX = 250

real(r8)    :: psm0                ! initial value for mean surface pressure ps
real(r8)    :: dsem0               ! initial value for mean dry static energy
real(r8)    :: msem0               ! initial value for mean moist static energy
real(r8)    :: kem0                ! initial value for mean kinetic energy
real(r8)    :: qcm0                ! initial value for mean cloudwater content
integer(i4) :: ntke                ! time level for TKE

integer(i4) :: izvctype_read       ! check vertical coordinate type in restarts
real(r8)    :: refatm_p0sl         ! constant reference pressure on sea-level
real(r8)    :: refatm_t0sl         ! constant reference temperature on sea-level
real(r8)    :: refatm_dt0lp        ! d (t0) / d (ln p0)
real(r8)    :: vcoord_vcflat       ! coordinate where levels become flat
real(r8)    :: zvc_params(KHMAX)   ! height levels

real(r8) :: refatm_delta_t ! temperature difference between sea level and stratosphere (for irefatm=2)
real(r8) :: refatm_h_scal  ! scale height (for irefatm=2)
real(r8) :: refatm_bvref   ! constant Brund-Vaisala-frequency for irefatm=3

integer :: izerr

! first record
read(rfid, iostat=izerr) psm0, dsem0, msem0, kem0, qcm0, ntke

if (izerr /= 0) then
   write(string1,*)'unable to read first record while searching for '//trim(progvar(ivar)%varname)
   call error_handler(E_ERR,'rw_binary_header', string1, source, revision, revdate)
endif

if (present(wfid)) then
   write(wfid, iostat=izerr) psm0, dsem0, msem0, kem0, qcm0, ntke

   if (izerr /= 0) then
      write(string1,*)'unable to write first record while searching for '//trim(progvar(ivar)%varname)
      call error_handler(E_ERR,'rw_binary_header', string1, source, revision, revdate)
   endif
endif

! second record
read(rfid, iostat=izerr) izvctype_read, refatm_p0sl, refatm_t0sl,   &
                        refatm_dt0lp, vcoord_vcflat, zvc_params
if (izerr /= 0) then
   write(string1,*)'unable to read second record while searching for '//trim(progvar(ivar)%varname)
   call error_handler(E_ERR,'rw_binary_header', string1, source, revision, revdate)
endif

if (present(wfid)) then
   write(wfid, iostat=izerr) izvctype_read, refatm_p0sl, refatm_t0sl,   &
                           refatm_dt0lp, vcoord_vcflat, zvc_params
   if (izerr /= 0) then
      write(string1,*)'unable to write second record while searching for '//trim(progvar(ivar)%varname)
      call error_handler(E_ERR,'rw_binary_header', string1, source, revision, revdate)
   endif
endif

! auxiliary record
if     ( izvctype_read > 0 .and. izvctype_read <= 100 ) then
   !write(*,*) "izvctype_read = ", izvctype_read
   continue

elseif ( (izvctype_read > 100) .and. (izvctype_read <= 200) ) then

   read(rfid, iostat=izerr) refatm_delta_t, refatm_h_scal
   if (izerr /= 0) then
      write(string1,*)'izvctype_read is ',izvctype_read, 'requiring us to read "refatm_delta_t, refatm_h_scal"'
      write(string2,*)'unable to read record while searching for '//trim(progvar(ivar)%varname)
      call error_handler(E_ERR,'rw_binary_header', string1, source, revision, revdate, text2=string2)
   endif

   if (present(wfid)) then
      write(wfid, iostat=izerr) refatm_delta_t, refatm_h_scal
      if (izerr /= 0) then
         write(string1,*)'izvctype_read is ',izvctype_read, 'requiring us to write "refatm_delta_t, refatm_h_scal"'
         write(string2,*)'unable to write record while searching for '//trim(progvar(ivar)%varname)
         call error_handler(E_ERR,'rw_binary_header', string1, source, revision, revdate, text2=string2)
      endif
   endif

elseif ( (izvctype_read > 200) .and. (izvctype_read <= 300) ) then

   read(rfid, iostat=izerr) refatm_bvref
   if (izerr /= 0) then
      write(string1,*)'izvctype_read is ',izvctype_read, 'requiring us to read "refatm_bvref"'
      write(string2,*)'unable to read record while searching for '//trim(progvar(ivar)%varname)
      call error_handler(E_ERR,'rw_binary_header', string1, source, revision, revdate, text2=string2)
   endif

   if (present(wfid)) then
      write(wfid, iostat=izerr) refatm_bvref
      if (izerr /= 0) then
         write(string1,*)'izvctype_read is ',izvctype_read, 'requiring us to write "refatm_bvref"'
         write(string2,*)'unable to write record while searching for '//trim(progvar(ivar)%varname)
         call error_handler(E_ERR,'rw_binary_header', string1, source, revision, revdate, text2=string2)
      endif
   endif

else
   write(string1,*) 'izvctype_read is ',izvctype_read,' is unsupported.'
   call error_handler(E_ERR,'rw_binary_header', string1, source, revision, revdate)
endif

end subroutine rw_binary_header


!------------------------------------------------------------------------
!>


subroutine decode_location(igdsbuf, nx, ny, lat1, lon1, latN, lonN)

integer(i4), intent(in)  :: igdsbuf(:)
integer,     intent(out) :: nx
integer,     intent(out) :: ny
real(r8),    intent(out) :: lat1
real(r8),    intent(out) :: lon1
real(r8),    intent(out) :: latN
real(r8),    intent(out) :: lonN

! Table to decode the record contents of igdsbuf
integer, parameter :: indx_numEW    =  5, &
                      indx_numNS    =  6, &
                      indx_startlat =  7, &
                      indx_startlon =  8, &
                      indx_endlat   = 10, &
                      indx_endlon   = 11

nx   = igdsbuf(indx_numEW)
ny   = igdsbuf(indx_numNS)

lat1 = real(igdsbuf(indx_startlat),r8) * 0.001_r8
lon1 = real(igdsbuf(indx_startlon),r8) * 0.001_r8
latN = real(igdsbuf(indx_endlat  ),r8) * 0.001_r8
lonN = real(igdsbuf(indx_endlon  ),r8) * 0.001_r8

end subroutine decode_location


!------------------------------------------------------------------------
!>


subroutine record_sizes(ivar, nx, ny, iz_countl, ilevtyp, ilev, ilevp1)

integer, intent(in) :: ivar
integer, intent(in) :: nx
integer, intent(in) :: ny
integer, intent(in) :: iz_countl
integer, intent(in) :: ilevtyp
integer, intent(in) :: ilev
integer, intent(in) :: ilevp1

if (    ilevtyp == 109) then
   progvar(ivar)%numZ        = ilevp1
   progvar(ivar)%numdims     = 3
elseif (ilevtyp == 110) then
   progvar(ivar)%numZ        = ilev
   progvar(ivar)%numdims     = 3
elseif (ilevtyp == 1) then
   progvar(ivar)%numZ        = 1
   progvar(ivar)%numdims     = 2
else
   write(string1,*)'unsupported vertical coordinate system of ',ilevtyp
   write(string2,*)'trying to find ',progvar(ivar)%varname
   write(string3,*)'grib table ',progvar(ivar)%tableID, &
                   ' grib variable ', progvar(ivar)%variableID
   call error_handler(E_ERR, 'record_sizes', string1, &
              source, revision, revdate, text2=string2, text3=string3)
endif

progvar(ivar)%levtypeID  = ilevtyp

if     (progvar(ivar)%dart_kind == KIND_U_WIND_COMPONENT) then
   progvar(ivar)%dimnames(1) = 'srlon'
   progvar(ivar)%dimnames(2) = 'rlat'
   progvar(ivar)%dimnames(3) = 'level'
   progvar(ivar)%dimlens(1)  = nrlon
   progvar(ivar)%dimlens(2)  = nrlat
   progvar(ivar)%dimlens(3)  = nlevel

elseif (progvar(ivar)%dart_kind == KIND_V_WIND_COMPONENT) then
   progvar(ivar)%dimnames(1) = 'rlon'
   progvar(ivar)%dimnames(2) = 'srlat'
   progvar(ivar)%dimnames(3) = 'level'
   progvar(ivar)%dimlens(1)  = nrlon
   progvar(ivar)%dimlens(2)  = nrlat
   progvar(ivar)%dimlens(3)  = nlevel

elseif (progvar(ivar)%dart_kind == KIND_VERTICAL_VELOCITY) then
   progvar(ivar)%dimnames(1) = 'rlon'
   progvar(ivar)%dimnames(2) = 'rlat'
   progvar(ivar)%dimnames(3) = 'level1'
   progvar(ivar)%dimlens(1)  = nrlon
   progvar(ivar)%dimlens(2)  = nrlat
   progvar(ivar)%dimlens(3)  = nlevel1

else
   progvar(ivar)%dimnames(1) = 'rlon'
   progvar(ivar)%dimnames(2) = 'rlat'
   progvar(ivar)%dimnames(3) = 'level'
   progvar(ivar)%dimlens(1)  = nrlon
   progvar(ivar)%dimlens(2)  = nrlat
   progvar(ivar)%dimlens(3)  = nlevel
endif

progvar(ivar)%numEW      = nx
progvar(ivar)%numNS      = ny
progvar(ivar)%slabN      = iz_countl
progvar(ivar)%varsize    = progvar(ivar)%numEW * &
                           progvar(ivar)%numNS * &
                           progvar(ivar)%numZ

return

end subroutine record_sizes


!------------------------------------------------------------------------
!>


subroutine insert_slab_in_state(ivar, slabID, nx, ny, rbuf, statevector)

integer,  intent(in)    :: ivar
integer,  intent(in)    :: slabID      ! index of slice we have in rbuf
integer,  intent(in)    :: nx
integer,  intent(in)    :: ny
real(r8), intent(in)    :: rbuf(nx,ny)
real(r8), intent(inout) :: statevector(:)

integer :: islab, istart, iend

if ( slabID < progvar(ivar)%slab1 .or. slabID > progvar(ivar)%slabN ) then
   write(string1,*)'slab out-of-range for ',trim(progvar(ivar)%varname)
   write(string2,*)'slab1 = ', progvar(ivar)%slab1, ', slab = ',slabID, &
                 ', slabN = ', progvar(ivar)%slabN
   call error_handler(E_ERR,'insert_slab_in_state', string1, &
              source, revision, revdate, text2=string2)
endif

islab  = slabID - progvar(ivar)%slab1 + 1
istart = progvar(ivar)%index1 + (islab -1)*nx*ny
iend   = istart + nx*ny - 1

if (debug > 99 .and. do_output()) &
   write(*,*)'stuffing slab ',slabID,' into ',istart, iend

statevector(istart:iend) = reshape(rbuf, (/ nx*ny /))

return

end subroutine insert_slab_in_state


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

if (debug > 99 .and. do_output()) then

   write(logfileunit,*)
   write(logfileunit,*)'define_var_dims knowledge'

   write(logfileunit,*)trim(progvar(ivar)%varname),' has original     dimnames: ', &
                   (/( trim(progvar(ivar)%dimnames(i))//' ',i=1,progvar(ivar)%numdims) /)
   write(logfileunit,*)trim(progvar(ivar)%varname),' repackaging into dimnames: ', &
                       (/ (trim(dimnames(i))//' ',i=1,ndims) /)

   write(logfileunit,*)'thus dimids ',dimids(1:ndims)
   write(     *     ,*)
   write(     *     ,*)'define_var_dims knowledge'
   write(     *     ,*)trim(progvar(ivar)%varname),' has original     dimnames: ', &
                   (/( trim(progvar(ivar)%dimnames(i))//' ',i=1,progvar(ivar)%numdims) /)
   write(     *     ,*)trim(progvar(ivar)%varname),' repackaging into dimnames: ', &
                       (/ (trim(dimnames(i))//' ',i=1,ndims) /)
   write(     *     ,*)'thus dimids ',dimids(1:ndims)

endif

return

end subroutine define_var_dims


!------------------------------------------------------------------------
!>


function Find_Variable_by_index(myindx, msgstring)
! Given an index into the DART state vector, return the index of metadata
! variable 'progvar' responsible for this portion of the state vector
integer,          intent(in) :: myindx
character(len=*), intent(in) :: msgstring
integer                      :: Find_Variable_by_index

integer :: ivar

Find_Variable_by_index = -1

FindIndex : do ivar = 1,nfields
   if ((myindx >= progvar(ivar)%index1)  .and. &
       (myindx <= progvar(ivar)%indexN)) then
      Find_Variable_by_index = ivar
      exit FindIndex
   endif
enddo FindIndex

if (Find_Variable_by_index < 0) then
   write(string1,*)'index ',myindx,' is out of range of all variables.'
   write(string2,*)'model size is ',model_size
   call error_handler(E_ERR, 'Find_Variable_by_index'//trim(msgstring), string1, &
                      source, revision, revdate, text2=string2 )
endif

end function Find_Variable_by_index


!>@ TODO FIXME put these in a separate module to indicate they came from COSMO
!##############################################################################
!> Both phi2phirot and rla2rlarot are part of the 
!> COSMO utilities (5.21)
!##############################################################################

FUNCTION  phi2phirot ( phi, rla, polphi, pollam )

!------------------------------------------------------------------------------
! Description:
!   This routine converts phi from the real geographical system to phi
!   in the rotated system.
!
! Method:
!   Transformation formulas for converting between these two systems.
!
! Originally part of geo2rotated_cosmo.f90 ...  
!------------------------------------------------------------------------------

integer, parameter :: wp=r8   ! to save on changing DART kind to COSMO kind

! Parameter list:
REAL (KIND=wp),     INTENT (IN)      ::        &
  polphi,  & ! latitude of the rotated north pole
  pollam,  & ! longitude of the rotated north pole
  phi,     & ! latitude in the geographical system
  rla        ! longitude in the geographical system

REAL (KIND=wp)                       ::        &
  phi2phirot ! latitude in the rotated system

! Local variables
REAL (KIND=wp)                           ::    &
  zsinpol, zcospol, zlampol, zphi, zrla, zarg1, zarg2, zrla1

REAL (KIND=wp),     PARAMETER            ::    &
  zrpi18 = 57.2957795_wp,                      & !
  zpir18 = 0.0174532925_wp

!------------------------------------------------------------------------------

! Begin function phi2phirot

  zsinpol  = SIN (zpir18 * polphi)
  zcospol  = COS (zpir18 * polphi)
  zlampol  =      zpir18 * pollam
  zphi     =      zpir18 * phi
  IF (rla > 180.0_wp) THEN
    zrla1  = rla - 360.0_wp
  ELSE
    zrla1  = rla
  ENDIF
  zrla     = zpir18 * zrla1

  zarg1    = SIN (zphi) * zsinpol
  zarg2    = COS (zphi) * zcospol * COS (zrla - zlampol)

  phi2phirot = zrpi18 * ASIN (zarg1 + zarg2)

END FUNCTION phi2phirot

!==============================================================================
!==============================================================================

! polgam = 0.0

FUNCTION  rla2rlarot ( phi, rla, polphi, pollam, polgam )

!------------------------------------------------------------------------------
!
! Description:
!   This routine converts lambda from the real geographical system to lambda 
!   in the rotated system.
!
! Method:
!   Transformation formulas for converting between these two systems.
!
!------------------------------------------------------------------------------

integer, parameter :: wp=r8   ! to save on changing DART kind to COSMO kind

! Parameter list:
REAL (KIND=wp),     INTENT (IN)      ::        &
  polphi,  & ! latitude of the rotated north pole
  pollam,  & ! longitude of the rotated north pole
  phi,     & ! latitude in geographical system
  rla        ! longitude in geographical system

REAL (KIND=wp),     INTENT (IN)      ::        &
  polgam      ! angle between the north poles of the systems

REAL (KIND=wp)                       ::        &
  rla2rlarot ! longitude in the the rotated system

! Local variables
REAL (KIND=wp)                           ::    &
  zsinpol, zcospol, zlampol, zphi, zrla, zarg1, zarg2, zrla1

REAL (KIND=wp),     PARAMETER            ::    &
  zrpi18 = 57.2957795_wp,                      & !
  zpir18 = 0.0174532925_wp

!------------------------------------------------------------------------------

! Begin function rla2rlarot

  zsinpol  = SIN (zpir18 * polphi)
  zcospol  = COS (zpir18 * polphi)
  zlampol  =      zpir18 * pollam
  zphi     =      zpir18 * phi
  IF (rla > 180.0_wp) THEN
    zrla1  = rla - 360.0_wp
  ELSE
    zrla1  = rla
  ENDIF
  zrla     = zpir18 * zrla1

  zarg1    = - SIN (zrla-zlampol) * COS(zphi)
  zarg2    = - zsinpol * COS(zphi) * COS(zrla-zlampol) + zcospol * SIN(zphi)

  IF (zarg2 == 0.0_wp) zarg2 = 1.0E-20_wp

  rla2rlarot = zrpi18 * ATAN2 (zarg1,zarg2)

  IF (polgam /= 0.0_wp) THEN
    rla2rlarot = polgam + rla2rlarot
    IF (rla2rlarot > 180._wp) rla2rlarot = rla2rlarot -360._wp
  ENDIF

END FUNCTION rla2rlarot

!##############################################################################
! Both phi2phirot and rla2rlarot are part of the 
! COSMO utilities (5.21)
!##############################################################################

!------------------------------------------------------------------------
!>


subroutine get_corners(lon, lat, nx, ny, gridlons, gridlats,  &
                ileft, iright, ifrac, jbot, jtop, jfrac, istatus)

real(r8), intent(in)  :: lon
real(r8), intent(in)  :: lat
integer,  intent(in)  :: nx
integer,  intent(in)  :: ny
real(r8), intent(in)  :: gridlons(nx)
real(r8), intent(in)  :: gridlats(ny)

integer,  intent(out) :: ileft
integer,  intent(out) :: iright
real(r8), intent(out) :: ifrac
integer,  intent(out) :: jbot
integer,  intent(out) :: jtop
real(r8), intent(out) :: jfrac
integer,  intent(out) :: istatus

integer, parameter :: OUTSIDE_HORIZONTALLY = 15

real(r8) :: dlon1, dlonT
real(r8) :: dlat1, dlatT

integer :: indarr(1)

! set output to facilitate early failed returns. 
! If all goes well, these get replaced.
ileft    = -1
iright   = -1
ifrac    = 0.0_r8
jbot     = -1
jtop     = -1
jfrac    = 0.0_r8
istatus  = 99

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

if (lon == gridlons(nx)) then
   ileft  = nx-1
   iright = nx
   ifrac  = 0.0_r8  ! fractional distance to min
else
   indarr = maxloc(gridlons , gridlons <= lon)
   ileft  = indarr(1)
   iright = ileft + 1
   dlon1  =     lon          - gridlons(ileft)
   dlonT  = gridlons(iright) - gridlons(ileft)
   ifrac  = dlon1 / dlonT 
endif

! latitudes are arranged 'south' to 'north', i.e. -90 to 90  (albeit rotated)
!
!    | ............ dlatT ............. |
!    | .. dlat .. |
!    |------------|---------------------|
! latbot         lat                  lattop

if (lat == gridlats(ny)) then
   jbot  = ny-1
   jtop  = ny
   jfrac = 0.0_r8
else
   indarr = maxloc( gridlats,  gridlats <= lat)
   jbot   = indarr(1)
   jtop   = jbot + 1
   dlat1  =     lat        - gridlats(jbot)
   dlatT  = gridlats(jtop) - gridlats(jbot)
   jfrac  = dlat1 / dlatT
endif

istatus = 0

if (debug > 99 .and. do_output()) then
   write(*,*)
   write(*,*)'longitude index to the  west and  east are ',ileft, iright
   write(*,*)'latitude  index to the south and north are ',jbot, jtop
   write(*,*)'lon  west, lon, lon  east',gridlons(ileft),lon,gridlons(iright)
   write(*,*)'lat south, lat, lat north',gridlats(jbot), lat,gridlats(jtop)
endif

return

end subroutine get_corners


!------------------------------------------------------------------------
!>


subroutine horizontal_interpolate(x, ivar, level_index, ileft, iright, ifrac, jbot, jtop, jfrac, &
                      layervalue, istatus) 

real(r8), intent(in) :: x(:)
integer,  intent(in) :: ivar
integer,  intent(in) :: level_index
integer,  intent(in) :: ileft
integer,  intent(in) :: iright
real(r8), intent(in) :: ifrac
integer,  intent(in) :: jbot
integer,  intent(in) :: jtop
real(r8), intent(in) :: jfrac

real(r8), intent(out) :: layervalue
integer,  intent(out) :: istatus

integer :: lowerleft
integer :: lowerright
integer :: upperright
integer :: upperleft
real(r8) :: ifracrem
real(r8) :: jfracrem

! figure out the dart state vector indices of the locations of interest

lowerleft  = ijk_to_dart(ivar, ileft,  jbot, level_index) 
lowerright = ijk_to_dart(ivar, iright, jbot, level_index)
upperright = ijk_to_dart(ivar, iright, jtop, level_index) 
upperleft  = ijk_to_dart(ivar, ileft,  jtop, level_index)

! apply the weights

jfracrem = 1.0_r8 - jfrac
ifracrem = 1.0_r8 - ifrac

layervalue =  jfracrem * ( ifracrem*x(lowerleft) + ifrac*x(lowerright) ) + &
              jfrac    * ( ifracrem*x(upperleft) + ifrac*x(upperright) )

if (debug > 99 .and. do_output()) then
   write(*,*)
   write(*,*)' ileft, jbot, level_index decompose to ', lowerleft, x( lowerleft), ifracrem
   write(*,*)'iright, jbot, level_index decompose to ',lowerright, x(lowerright), ifrac
   write(*,*)'iright, jtop, level_index decompose to ',upperright, x(upperright), jfracrem
   write(*,*)' ileft, jtop, level_index decompose to ', upperleft, x( upperleft), jfrac
   write(*,*)' horizontal value is ',layervalue
endif

return

end subroutine horizontal_interpolate


!------------------------------------------------------------------------
!>

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
!>


subroutine get_level_indices(height, dart_kind, ibelow, iabove, levelfrac, istatus)

! The heights in the vcoord array are from the top of the atmosphere down.
! vcoord%level1(   1   ) = top of atmosphere
! vcoord%level1(nlevel1) = bottom level of model, earth surface
! Only W is staggered to use all nlevel1 dimensions.
! All other variables use the level midpoints, stored in vcoord%level
!
!
!         space <--------------------------------------------| surface
! layer     0 ... 10                                 11 .... 51
!                 | ............ dztot ............. |
!                 | ... dz ... |
!                 |------------|---------------------|
!               ztop         height                 zbot
!              iabove                              ibelow
!              big numbers                      small numbers

!         space <--------------------------------------------| surface
! layer     0 ... 1                                  2 ..... 51
!                 | ............ dztot ............. |
!                 | ... dz ... |
!                 |------------|---------------------| ...... 0
!               22000        height                21000
!               iabove                             ibelow
!
! A note about the use of the minloc() function. I struggled with the logic
! on this one and could not get maxloc() to work despite the fact I thought
! it was the right choice. Using minloc() _worked_ despite the fact I thought
! it was the wrong choice. Perhaps the fact that the array is descending
! rather than ascending made it invert the logic. TJH

real(r8), intent(in)  :: height
integer,  intent(in)  :: dart_kind
integer,  intent(out) :: ibelow
integer,  intent(out) :: iabove
real(r8), intent(out) :: levelfrac
integer,  intent(out) :: istatus

integer, parameter :: OUTSIDE_VERTICALLY   = 16

integer :: indarr(1)
real(r8) :: dz, dztot

if ( dart_kind == KIND_VERTICAL_VELOCITY ) then

   if (height > vcoord%level1(1) .or. height < vcoord%level1(nlevel1)) then
      istatus = OUTSIDE_VERTICALLY
      return
   endif

   if (height == vcoord%level1(nlevel1)) then
      iabove = nlevel1
      ibelow = nlevel1
      levelfrac = 1.0_r8
   else
      indarr = minloc( vcoord%level1,  vcoord%level1 >= height )
      iabove = indarr(1)
      ibelow = iabove + 1
      dz     = vcoord%level1(iabove) - height
      dztot  = vcoord%level1(iabove) - vcoord%level1(ibelow)
      levelfrac = dz / dztot
   endif

else

   if (height > vcoord%level(1) .or. height < vcoord%level(nlevel)) then
      istatus = OUTSIDE_VERTICALLY
      return
   endif

   if (height == vcoord%level(nlevel)) then
      iabove = nlevel
      ibelow = nlevel
      levelfrac = 1.0_r8
   else
      indarr = minloc( vcoord%level,  vcoord%level >= height )
      iabove = indarr(1)
      ibelow = iabove + 1
      dz     = vcoord%level(iabove) - height
      dztot  = vcoord%level(iabove) - vcoord%level(ibelow)
      levelfrac = dz / dztot
   endif

   if (debug > 99 .and. do_output()) then
      write(*,*)
      write(*,*)'computine vertical levels for DART kind ',dart_kind
      write(*,*)'iabove, levelfrac, ibelow ', iabove, levelfrac, ibelow
      write(*,*)'above, height, below ', vcoord%level(iabove), height, vcoord%level(ibelow)
      write(*,*)
   endif
endif

istatus = 0

return

end subroutine get_level_indices


!------------------------------------------------------------------------
!>


subroutine get_dart_indices(slab_index, index1, indexN, desired)
integer, intent(in)  :: slab_index
integer, intent(out) :: index1
integer, intent(out) :: indexN
logical, intent(out) :: desired

integer :: ivar, kindex

desired = .false.
index1  = -1
indexN  = -1

IndexLoop : do ivar = 1,nfields

   if ( slab_index < progvar(ivar)%slab1 .or. slab_index > progvar(ivar)%slabN ) cycle IndexLoop

   desired = .true.
   kindex  = slab_index - progvar(ivar)%slab1 + 1

   index1  = progvar(ivar)%index1 + (kindex -1) * progvar(ivar)%dimlens(1)*progvar(ivar)%dimlens(2)
   indexN  =               index1 +               progvar(ivar)%dimlens(1)*progvar(ivar)%dimlens(2) - 1

!  write(*,*)'slab ',slab_index,' is level ',kindex,' of variable ',ivar,' index1,N are ',index1, indexN

enddo IndexLoop

return

end subroutine get_dart_indices


!------------------------------------------------------------------------
!>


end module model_mod

! <next few lines under version control, do not edit>
! $URL: $
! $Id: $
! $Revision: $
! $Date: $

