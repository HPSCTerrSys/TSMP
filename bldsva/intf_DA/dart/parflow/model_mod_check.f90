! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! DART $Id: model_mod_check.f90 Wed Apr 11 20:26:43 CEST 2018 $

program model_mod_check

!----------------------------------------------------------------------
! purpose: test routines.  this version for models with oned locations.
!----------------------------------------------------------------------

use        types_mod, only : r8, digits12, metadatalength, MISSING_R8

use    utilities_mod, only : initialize_utilities, finalize_utilities, nc_check, &
                             open_file, close_file, find_namelist_in_file, &
                             check_namelist_read, E_MSG, E_ERR, error_handler, &
                             do_output

use     location_mod, only : location_type, set_location, write_location, get_dist, &
                             query_location, LocationDims, get_location,            &
                             VERTISUNDEF, VERTISSURFACE, VERTISLEVEL,               &
                             VERTISPRESSURE, VERTISHEIGHT, VERTISSCALEHEIGHT

use     obs_kind_mod, only : get_raw_obs_kind_name, get_raw_obs_kind_index

use  assim_model_mod, only : open_restart_read, open_restart_write, close_restart, &
                             aread_state_restart, awrite_state_restart, &
                             netcdf_file_type, aoutput_diagnostics, &
                             init_diag_output, finalize_diag_output, static_init_assim_model

use time_manager_mod, only : time_type, set_calendar_type, GREGORIAN, &
                             read_time, get_time, set_time,  &
                             print_time, print_date, write_time, operator(-)

use        model_mod, only : static_init_model, get_model_size, get_state_meta_data, &
                             get_state_vector, get_parflow_filename, model_interpolate, &
                             get_parflow_id
use netcdf
use typesizes

implicit none

! version controlled file description for error handling, do not edit
character(len=*), parameter :: source   = "$URL: model_mod_check.f90 $"
character(len=*), parameter :: revision = "$Revision: Bonn $"
character(len=*), parameter :: revdate  = "$Date: Wed Apr 11 2018 $"

character(len=512) :: string1, string2, string3

!------------------------------------------------------------------
! The namelist variables
!------------------------------------------------------------------

character(len=256)            :: dart_input_file       = 'dart_ics'
character(len=256)            :: output_file           = 'check_me'
logical                       :: advance_time_present  = .FALSE.
logical                       :: verbose               = .FALSE.
integer                       :: test1thru             = -1
integer                       :: x_ind                 = -1
real(r8), dimension(3)        :: loc_of_interest       = -1.0_r8
character(len=metadatalength) :: kind_of_interest      = 'ANY'
character(len=metadatalength) :: interp_test_vertcoord = 'VERTISHEIGHT'

real(r8) :: interp_test_lonrange(2)  = (/ 320.0, 90.0 /)
real(r8) :: interp_test_dlon         = 1.0
real(r8) :: interp_test_latrange(2)  = (/ 0.0, 75.0 /)
real(r8) :: interp_test_dlat         = 1.0
real(r8) :: interp_test_vertrange(2) = (/ 1.0, 2.0 /)
real(r8) :: interp_test_dvert        = 1.0

namelist /model_mod_check_nml/ dart_input_file, output_file,        &
                        advance_time_present, test1thru, x_ind,     &
                        loc_of_interest, kind_of_interest, verbose, &
                        interp_test_vertcoord,                      &
                        interp_test_lonrange, interp_test_dlon,     &
                        interp_test_latrange, interp_test_dlat,     &
                        interp_test_vertrange, interp_test_dvert

!----------------------------------------------------------------------
! integer :: numlons, numlats, numlevs

integer :: ios_out, iunit, io
integer :: x_size
integer :: myobsIndex           !Observation Type to use in interpolation

type(time_type)       :: model_time, adv_to_time
real(r8), allocatable :: statevector(:)

integer :: pfbid, num_failed

character(len=metadatalength) :: state_meta(1)
character(len=256)            :: parflow_input_file
type(netcdf_file_type) :: ncFileID
type(location_type) :: loc

real(r8) :: interp_val

!----------------------------------------------------------------------
! This portion checks the geometry information. 
!----------------------------------------------------------------------

call initialize_utilities(progname='model_mod_check')
call set_calendar_type(GREGORIAN)

write(*,*)
write(*,*)'Reading the namelist to get the input filename.'

call find_namelist_in_file("input.nml", "model_mod_check_nml", iunit)
read(iunit, nml = model_mod_check_nml, iostat = io)
call check_namelist_read(iunit, io, "model_mod_check_nml")

loc = set_location(loc_of_interest(1), loc_of_interest(2), loc_of_interest(3), VERTISHEIGHT)
myobsIndex = get_raw_obs_kind_index(kind_of_interest)

if (test1thru < 1) goto 999

write(*,*)
write(*,*)'model_mod_check ....1'
write(*,*)'static_init_model test STARTING ...'
call static_init_model()
write(*,*)'static_init_model test COMPLETE ...'

if (test1thru < 2) goto 999

write(*,*)
write(*,*)'model_mod_check ....2'
write(*,*)'get_model_size test STARTING ...'
x_size = get_model_size()
write(*,*)'get_model_size test : state vector has length',x_size
write(*,*)'get_model_size test COMPLETE ...'

!----------------------------------------------------------------------
! Write a supremely simple restart file. Most of the time, I just use
! this as a starting point for a Matlab function that replaces the 
! values with something more complicated.
!----------------------------------------------------------------------

if (test1thru < 3) goto 999

allocate(statevector(x_size))
write(*,*)
write(*,*)'model_mod_check ....3'
write(*,*)'initialize statevector and set_time'
statevector = 1.0_r8;
model_time  = set_time(21600, 149446)   ! 06Z 4 March 2010
write(*,*)'set_time complete ...'

!----------------------------------------------------------------------
! Open a test DART initial conditions file.
! Reads the valid time, the state, and (possibly) a target time.
!----------------------------------------------------------------------
if (test1thru < 4) goto 999

write(*,*)
write(*,*)'model_mod_check ....4'
write(*,*)'Reading '//trim(dart_input_file)

iunit = open_restart_read(dart_input_file)
if ( advance_time_present ) then
   call aread_state_restart(model_time, statevector, iunit, adv_to_time)
else
   call aread_state_restart(model_time, statevector, iunit)
endif

call close_restart(iunit)
call print_time( model_time,'model_mod_check:model time')
call print_time( model_time,'model_mod_check:model time')

!----------------------------------------------------------------------
! Output the state vector to a netCDF file ...
! This is the same procedure used by 'perfect_model_obs' & 'filter'
! init_diag_output()
! aoutput_diagnostics()
! finalize_diag_output()
!----------------------------------------------------------------------

if (test1thru < 5) goto 999

write(*,*)
write(*,*)'model_mod_check ....5'
write(*,*)'Exercising the netCDF routines.'
write(*,*)'Creating '//trim(output_file)//'.nc'

state_meta(1) = 'restart test'
ncFileID = init_diag_output(trim(output_file),'just testing a restart', 1, state_meta)

call aoutput_diagnostics(ncFileID, model_time, statevector, 1)

call nc_check( finalize_diag_output(ncFileID), 'model_mod_check:main', 'finalize')


!----------------------------------------------------------------------
! Checking get_state_meta_data (and get_state_indices, get_state_kind)
!----------------------------------------------------------------------

if (test1thru < 6) goto 999

write(*,*)
write(*,*)'model_mod_check ....6'
write(*,*)'Checking metadata routines for dart index ....'

call check_meta_data( x_ind )

!----------------------------------------------------------------------
! Trying to find the state vector index closest to a particular ...
! Checking for valid input is tricky ... we don't know much. 
!----------------------------------------------------------------------

if (test1thru < 7) goto 999

write(*,*)
write(*,*)'model_mod_check ....7'
write(*,*)'Find location of state vector...'

if ( loc_of_interest(1) >= 0.0_r8 ) call find_closest_gridpoint( loc_of_interest )

!----------------------------------------------------------------------
! Check the interpolation - print initially to STDOUT
!----------------------------------------------------------------------

if (test1thru < 8) goto 999

write(*,*)
write(*,*)'model_mod_check ....8'
write(*,*)'Testing model_interpolate ...'

if ( myobsIndex < 0 ) then
   write(*,*)'WARNING'
   write(*,*)'WARNING - input.nml:model_mod_check does not have a known "kind_of_interest"'
   write(*,*)'WARNING - skipping the model_interpolate tests.'
   write(*,*)'WARNING'
   goto 200 
endif

call model_interpolate(statevector, loc, myobsIndex , interp_val, ios_out)

if ( ios_out == 0 ) then 
   write(*,*)'model_interpolate SUCCESS: The interpolated value is ',interp_val
else
   write(*,*)'model_interpolate ERROR: model_interpolate failed with error code ',ios_out
endif

!----------------------------------------------------------------------
! Exhaustive test of model_interpolate.
!----------------------------------------------------------------------

if (test1thru < 9 ) goto 999

   write(*,*)
   write(*,*)'Exhaustive test of model interpolate ... please be patient.'

   num_failed = test_interpolate_range( statevector, myobsIndex, &
                                        interp_test_dlon,        &
                                        interp_test_dlat,        &
                                        interp_test_dvert,       &
                                        interp_test_vertcoord,   &
                                        interp_test_lonrange,    &
                                        interp_test_latrange,    &
                                        interp_test_vertrange,   &
                                        kind_of_interest,        &
                                        verbose )
 200 continue

!----------------------------------------------------------------------
! Writing dart output file.
!----------------------------------------------------------------------

if (test1thru < 10) goto 999

pfbid = get_parflow_id()
parflow_input_file = get_parflow_filename(pfbid)

write(*,*)
write(*,*)'Reading restart files from  '//trim(parflow_input_file)

call get_state_vector(statevector, 1, model_time)

write(*,*)
write(*,*)'Writing data into '//trim(output_file)
iunit = open_restart_write(output_file)
call awrite_state_restart(model_time, statevector, iunit)
call close_restart(iunit)

!----------------------------------------------------------------------
! Open a test DART initial conditions file.
! Reads the valid time, the state, and (possibly) a target time.
!----------------------------------------------------------------------

if (test1thru < 11) goto 999

write(*,*)
write(*,*)'Reading '//trim(output_file)

iunit = open_restart_read(output_file)
if ( advance_time_present ) then
   call aread_state_restart(model_time, statevector, iunit, adv_to_time)
else
   call aread_state_restart(model_time, statevector, iunit)
endif
call close_restart(iunit)

call print_date( model_time,'model_mod_check:model date')
call print_time( model_time,'model_mod_check:model time')

 999 continue

call finalize_utilities(progname='model_mod_check')

! end of main program

contains

!---------------------------------------------------------
!> Checking meta data for the given location in state vector

subroutine check_meta_data( iloc )

integer, intent(in) :: iloc
type(location_type) :: loc
integer             :: var_type
character(len=129)  :: string1

write(*,*)
write(*,*)'Checking metadata routines.'

call get_state_meta_data( iloc, loc, var_type)

call write_location(0, loc, charstring=string1)
write(*,*)' indx ',iloc,' is type ',var_type,trim(string1)

end subroutine check_meta_data


!---------------------------------------------------------
!>Simple exhaustive search to find the indices into the
!>state vector of a particular location.

subroutine find_closest_gridpoint( loc_of_interest )

real(r8), dimension(:), intent(in) :: loc_of_interest

type(location_type) :: loc0, loc1
integer  :: i, which_vert, var_type
real(r8) :: plon, plat, plev, vals(3)
real(r8) :: closest
real(r8), allocatable, dimension(:) :: thisdist
logical :: matched

write(*,*)
write(*,'(''Checking for the index in the state vector that is closest to '')')
write(*,'(''lon/lat/lev '',3(1x,f20.14))')loc_of_interest(1:LocationDims)

allocate( thisdist(get_model_size()) )
thisdist  = 9999999999.9_r8         ! really far away 
matched   = .false.

plon = loc_of_interest(1)
plat = loc_of_interest(2)
plev = loc_of_interest(3)

! Since there can be multiple variables with
! identical distances, we will just cruise once through 
! the array and come back to find all the 'identical' values.
do i = 1,get_model_size()

   ! REALLY INEFFICIENT, but grab the 'which_vert' from the
   ! grid and set our target location to have the same.
   ! Then, compute the distance and compare.

   call get_state_meta_data(i, loc1, var_type)
   which_vert  = nint( query_location(loc1) )
   loc0        = set_location(plon, plat, plev, which_vert)
   thisdist(i) = get_dist( loc1, loc0, no_vert=.false.)
   matched    = .true.
enddo

if (.not. matched) then
   write(*,*)'No state vector elements of type '//trim(kind_of_interest)
   return
endif

! Now that we know  ... report 

closest = minval(thisdist)
if (closest == 9999999999.9_r8) then
   write(*,*)'No closest gridpoint found'
   return
endif

matched = .false.
do i = 1,get_model_size()

   if ( thisdist(i) == closest ) then
      call get_state_meta_data(i, loc1)
      vals = get_location(loc1)
      write(*,'(A,1x,3(1x,f20.14),i10)') 'matched at ',vals, i
      matched = .true.
   endif

enddo

if ( .not. matched ) then
   write(*,*)'Nothing matched the closest gridpoint'
endif


deallocate( thisdist )

end subroutine find_closest_gridpoint


!-------------------------------------------------------------------------------
!> Interpolate over a range of lat, lon, and vert values.
!> Returns the number of failures.
!> Exercises model_mod:model_interpolate().
!> This will result in a netCDF file with all salient metadata.

function test_interpolate_range( statevector,           &
                                 obsIndex,              &
                                 interp_test_dlon,      &
                                 interp_test_dlat,      &
                                 interp_test_dvert,     &
                                 interp_test_vertcoord, &
                                 interp_test_lonrange,  &
                                 interp_test_latrange,  &
                                 interp_test_vertrange, &
                                 quantity_string,       &
                                 verbose )

real(r8)              , intent(in) :: statevector(:)
integer               , intent(in) :: obsIndex
real(r8)              , intent(in) :: interp_test_dlon
real(r8)              , intent(in) :: interp_test_dlat
real(r8)              , intent(in) :: interp_test_dvert
character(len=*)      , intent(in) :: interp_test_vertcoord
real(r8), dimension(2), intent(in) :: interp_test_latrange
real(r8), dimension(2), intent(in) :: interp_test_lonrange
real(r8), dimension(2), intent(in) :: interp_test_vertrange
character(len=*),       intent(in) :: quantity_string
logical               , intent(in) :: verbose

integer :: test_interpolate_range

! Local variables

character(len=*), parameter :: routine = 'test_interpolate_range'

real(r8), allocatable :: lon(:), lat(:), vert(:)
real(r8), allocatable :: field(:,:,:,:)
integer,  allocatable :: all_ios_out(:,:)
real(r8) :: lonrange_top
integer :: nlon, nlat, nvert
integer :: ilon, jlat, kvert, nfailed
character(len=128) :: ncfilename, txtfilename

character(len=8)      :: crdate      ! needed by F90 DATE_AND_TIME intrinsic
character(len=10)     :: crtime      ! needed by F90 DATE_AND_TIME intrinsic
character(len=5)      :: crzone      ! needed by F90 DATE_AND_TIME intrinsic
integer, dimension(8) :: values      ! needed by F90 DATE_AND_TIME intrinsic

integer, parameter :: ens_size = 1
integer :: ncid, nlonDimID, nlatDimID, nvertDimID
integer :: VarID(ens_size), lonVarID, latVarID, vertVarID

character(len=256)  :: output_file = 'check_me'
character(len=32)   :: field_name
type(location_type) :: loc
integer :: iunit, ios_out(ens_size), imem
integer :: vertcoord

test_interpolate_range = 0

if ((interp_test_dlon < 0.0_r8) .or. (interp_test_dlat < 0.0_r8)) then
   if ( do_output() ) then
      write(*,'(A)')    'Skipping the rigorous interpolation test because one of'
      write(*,'(A)')    'interp_test_dlon,interp_test_dlat are < 0.0'
      write(*,'(A,I2)') 'interp_test_dlon  = ',interp_test_dlon
      write(*,'(A,I2)') 'interp_test_dlat  = ',interp_test_dlat
      write(*,'(A,I2)') 'interp_test_dvert = ',interp_test_dvert
   endif
   return
endif

vertcoord = get_location_index(interp_test_vertcoord)

write( ncfilename,'(a,a)')trim(output_file),'_interptest.nc'
write(txtfilename,'(a,a)')trim(output_file),'_interptest.m'

! for longitude, allow wrap.
lonrange_top = interp_test_lonrange(2)
if (interp_test_lonrange(2) < interp_test_lonrange(1)) &
   lonrange_top = interp_test_lonrange(2) + 360.0_r8

! round down to avoid exceeding the specified range
nlat  = aint(( interp_test_latrange(2) - interp_test_latrange(1))  / interp_test_dlat) + 1
nlon  = aint((            lonrange_top - interp_test_lonrange(1))  / interp_test_dlon) + 1
nvert = aint((interp_test_vertrange(2) - interp_test_vertrange(1)) / interp_test_dvert) + 1

iunit = open_file(trim(txtfilename), action='write')
write(iunit,'(''missingvals = '',f12.4,'';'')')MISSING_R8
write(iunit,'(''nlon = '',i8,'';'')')nlon
write(iunit,'(''nlat = '',i8,'';'')')nlat
write(iunit,'(''nvert = '',i8,'';'')')nvert
write(iunit,'(''interptest = [ ... '')')

allocate(lon(nlon), lat(nlat), vert(nvert), field(nlon,nlat,nvert,ens_size))
allocate(all_ios_out(nlon*nlat*nvert,ens_size))

all_ios_out = 0 ! assume successful interpolation for every grid location, all members.
nfailed = 0

do ilon = 1, nlon
   lon(ilon) = interp_test_lonrange(1) + real(ilon-1,r8) * interp_test_dlon
   if (lon(ilon) >= 360.0_r8) lon(ilon) = lon(ilon) - 360.0_r8
   if (lon(ilon) <    0.0_r8) lon(ilon) = lon(ilon) + 360.0_r8
   do jlat = 1, nlat
      lat(jlat) = interp_test_latrange(1) + real(jlat-1,r8) * interp_test_dlat
      do kvert = 1, nvert
         vert(kvert) = interp_test_vertrange(1) + real(kvert-1,r8) * interp_test_dvert
         loc = set_location(lon(ilon), lat(jlat), vert(kvert), vertcoord)

         call model_interpolate(statevector, loc, myobsIndex, &
                                field(ilon,jlat,kvert,1), ios_out(1))

         write(iunit,*) field(ilon,jlat,kvert,:)
         if (any(ios_out /= 0)) then

            nfailed    = nfailed + 1
            ! don't really care which location was causing the failure
            all_ios_out(nfailed,:) = ios_out

            if (verbose) then
               write(string1,*) 'interpolation return code was', ios_out
               write(string2,'(''ilon,jlat,kvert,lon,lat,vert'',3(1x,i6),3(1x,f14.6))') &
                                 ilon,jlat,kvert,lon(ilon),lat(jlat),vert(kvert)
               call error_handler(E_MSG, routine, string1, &
                                  source, revision, revdate, text2=string2)
            endif

         endif
      enddo
   enddo
enddo

write(iunit,'(''];'')')
write(iunit,'(''datmat = reshape(interptest,nvert,nlat,nlon,nens);'')')
write(iunit,'(''datmat = permute(datmat,[4,1,2,3]);'')')
write(iunit,'(''datmat(datmat == missingvals) = NaN;'')')
call close_file(iunit)

if ( do_output() ) then
   write(*,'(A)')     '-------------------------------------------------------------'
   write(*,'(A,I10)') 'total  interpolations : ', nlon*nlat*nvert
   write(*,'(A,I10)') 'failed interpolations : ', nfailed
   write(*,'(A)')     '-------------------------------------------------------------'
endif

! Write out the netCDF file for easy exploration.

call DATE_AND_TIME(crdate,crtime,crzone,values)
write(string1,'(''YYYY MM DD HH MM SS = '',i4,5(1x,i2.2))') &
                  values(1), values(2), values(3), values(5), values(6), values(7)

call nc_check( nf90_create(path=trim(ncfilename), cmode=NF90_clobber, ncid=ncid), &
                  routine, 'open '//trim(ncfilename))
call nc_check( nf90_put_att(ncid, NF90_GLOBAL, 'creation_date' ,trim(string1) ), &
                  routine, 'creation put '//trim(ncfilename))

! Define dimensions

call nc_check(nf90_def_dim(ncid=ncid, name='lon', len=nlon, &
        dimid = nlonDimID),routine, 'nlon def_dim '//trim(ncfilename))

call nc_check(nf90_def_dim(ncid=ncid, name='lat', len=nlat, &
        dimid = nlatDimID),routine, 'nlat def_dim '//trim(ncfilename))

call nc_check(nf90_def_dim(ncid=ncid, name='vert', len=nvert, &
        dimid = nvertDimID),routine, 'nvert def_dim '//trim(ncfilename))

! Define variables

call nc_check(nf90_def_var(ncid=ncid, name='lon', xtype=nf90_double, &
        dimids=nlonDimID, varid=lonVarID), routine, &
                 'lon def_var '//trim(ncfilename))
call nc_check(nf90_put_att(ncid, lonVarID, 'range', interp_test_lonrange), &
           routine, 'put_att lonrange '//trim(ncfilename))
call nc_check(nf90_put_att(ncid, lonVarID, 'cartesian_axis', 'X'),   &
           routine, 'lon cartesian_axis '//trim(ncfilename))

call nc_check(nf90_def_var(ncid=ncid, name='lat', xtype=nf90_double, &
        dimids=nlatDimID, varid=latVarID), routine, &
                 'lat def_var '//trim(ncfilename))
call nc_check(nf90_put_att(ncid, latVarID, 'range', interp_test_latrange), &
           routine, 'put_att latrange '//trim(ncfilename))
call nc_check(nf90_put_att(ncid, latVarID, 'cartesian_axis', 'Y'),   &
           routine, 'lat cartesian_axis '//trim(ncfilename))

call nc_check(nf90_def_var(ncid=ncid, name='vert', xtype=nf90_double, &
        dimids=nvertDimID, varid=vertVarID), routine, &
                 'vert def_var '//trim(ncfilename))
call nc_check(nf90_put_att(ncid, vertVarID, 'range', interp_test_vertcoord), &
           routine, 'put_att vertrange '//trim(ncfilename))
call nc_check(nf90_put_att(ncid, vertVarID, 'cartesian_axis', 'Z'),   &
           routine, 'vert cartesian_axis '//trim(ncfilename))

! loop over ensemble members
do imem = 1, ens_size
   if ( ens_size > 1) then
      write(field_name,'(A,I2)') "field_",imem
   else
      field_name = "field"
   endif
   call nc_check(nf90_def_var(ncid=ncid, name=field_name, xtype=nf90_double, &
           dimids=(/ nlonDimID, nlatDimID, nvertDimID /), varid=VarID(imem)), routine, &
                    'field def_var '//trim(ncfilename))
   call nc_check(nf90_put_att(ncid, VarID(imem), 'long_name', quantity_string), &
              routine, 'put_att field long_name '//trim(ncfilename))
   call nc_check(nf90_put_att(ncid, VarID(imem), '_FillValue', MISSING_R8), &
              routine, 'put_att field FillValue '//trim(ncfilename))
   call nc_check(nf90_put_att(ncid, VarID(imem), 'missing_value', MISSING_R8), &
              routine, 'put_att field missing_value '//trim(ncfilename))
   call nc_check(nf90_put_att(ncid, VarID(imem), 'interp_test_vertcoord', interp_test_vertcoord ), &
              routine, 'put_att field interp_test_vertcoord '//trim(ncfilename))
enddo

! Leave define mode so we can fill the variables.
call nc_check(nf90_enddef(ncid), &
              routine,'field enddef '//trim(ncfilename))

! Fill the variables
call nc_check(nf90_put_var(ncid, lonVarID, lon), &
              routine,'lon put_var '//trim(ncfilename))
call nc_check(nf90_put_var(ncid, latVarID, lat), &
              routine,'lat put_var '//trim(ncfilename))
call nc_check(nf90_put_var(ncid, vertVarID, vert), &
              routine,'vert put_var '//trim(ncfilename))

do imem = 1, ens_size
   call nc_check(nf90_put_var(ncid, VarID(imem), field(:,:,:,imem)), &
                 routine,'field put_var '//trim(ncfilename))
enddo

! tidy up
call nc_check(nf90_close(ncid), routine,'close '//trim(ncfilename))

deallocate(lon, lat, vert, field)
deallocate(all_ios_out)

test_interpolate_range = nfailed

end function test_interpolate_range


!-------------------------------------------------------------------------------
!> need to convert the character string for the test vertical coordinate into
!> the corresponding dart index.

function  get_location_index(test_vertcoord)

character(len=*) , intent(in) :: test_vertcoord

integer :: get_location_index

select case (test_vertcoord)
   case ('VERTISUNDEF')
      get_location_index = VERTISUNDEF
   case ('VERTISSURFACE')
      get_location_index = VERTISSURFACE
   case ('VERTISLEVEL')
      get_location_index = VERTISLEVEL
   case ('VERTISPRESSURE')
      get_location_index = VERTISPRESSURE
   case ('VERTISHEIGHT')
      get_location_index = VERTISHEIGHT
   case ('VERTISSCALEHEIGHT')
      get_location_index = VERTISSCALEHEIGHT
   case default
      get_location_index = VERTISUNDEF
end select

end function  get_location_index



end program model_mod_check
