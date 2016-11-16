! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id: model_mod_check.f90 6311 2013-07-17 22:04:54Z thoar $

program model_mod_check

!----------------------------------------------------------------------
! purpose: test routines.  this version for models with oned locations.
!----------------------------------------------------------------------

use        types_mod, only : r8, digits12, metadatalength

use    utilities_mod, only : initialize_utilities, finalize_utilities, nc_check, &
                             open_file, close_file, find_namelist_in_file, &
                             check_namelist_read
use     location_mod, only : location_type, set_location, write_location, get_dist, &
                             query_location, LocationDims, get_location,            &
                             VERTISHEIGHT

use     obs_kind_mod, only : get_raw_obs_kind_name, get_raw_obs_kind_index

use  assim_model_mod, only : open_restart_read, open_restart_write, close_restart, &
                             aread_state_restart, awrite_state_restart, &
                             netcdf_file_type, aoutput_diagnostics, &
                             init_diag_output, finalize_diag_output, static_init_assim_model

use time_manager_mod, only : time_type, set_calendar_type, GREGORIAN, &
                             read_time, get_time, set_time,  &
                             print_time, print_date, write_time, operator(-)

use        model_mod, only : static_init_model, get_model_size, get_state_meta_data, &
                             get_state_vector, get_parflow_filename, model_interpolate
use netcdf
use typesizes

implicit none

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL: https://proxy.subversion.ucar.edu/DAReS/DART/releases/Lanai/models/template/model_mod_check.f90 $"
character(len=32 ), parameter :: revision = "$Revision: 6311 $"
character(len=128), parameter :: revdate  = "$Date: 2013-07-18 00:04:54 +0200 (Thu, 18 Jul 2013) $"

!------------------------------------------------------------------
! The namelist variables
!------------------------------------------------------------------

character (len = 256) :: dart_input_file               = 'dart_ics'
character (len = 129) :: output_file                   = 'check_me'
logical               :: advance_time_present          = .FALSE.
logical               :: verbose                       = .FALSE.
integer                :: test1thru                    = -1
integer               :: x_ind                         = -1
real(r8), dimension(3) :: loc_of_interest              = -1.0_r8
character(len=metadatalength) :: kind_of_interest      = 'ANY'
character(len=metadatalength) :: interp_test_vertcoord = 'VERTISHEIGHT'

namelist /model_mod_check_nml/ dart_input_file, output_file,        &
                        advance_time_present, test1thru, x_ind,     &
                        loc_of_interest, kind_of_interest, verbose, &
                        interp_test_vertcoord

!----------------------------------------------------------------------
! integer :: numlons, numlats, numlevs

integer :: in_unit, out_unit, ios_out, iunit, io, offset
integer :: x_size
integer :: year, month, day, hour, minute, second
integer :: secs, days
integer :: myobsIndex           !Observation Type to use in interpolation

type(time_type)       :: model_time, adv_to_time
real(r8), allocatable :: statevector(:)

integer :: i, skip

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
write(*,*)'Exhaustive test of model interpolate not written yet ...'

 200 continue

!----------------------------------------------------------------------
! Writing dart output file.
!----------------------------------------------------------------------

if (test1thru < 10) goto 999

parflow_input_file = get_parflow_filename()

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

call finalize_utilities()

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
character(len=129)  :: string1
real(r8), allocatable, dimension(:) :: thisdist
logical :: matched

write(*,*)
write(*,'(''Checking for the index in the state vector that is closest to '')')
write(*,'(''lon/lat/lev'',3(1x,f14.5))')loc_of_interest(1:LocationDims)

allocate( thisdist(get_model_size()) )
thisdist  = 9999999999.9_r8         ! really far away 
matched   = .false.
!

plon = loc_of_interest(1)
plat = loc_of_interest(2)
plev = loc_of_interest(3)


! Since there can be multiple variables with
! identical distances, we will just cruise once through 
! the array and come back to find all the 'identical' values.
do i = 1,get_model_size()

   ! Really inefficient, but grab the 'which_vert' from the
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
      write(*,'(3(1x,f14.5),i3)') vals, i
      matched = .true.
   endif

enddo

if ( .not. matched ) then
   write(*,*)'Nothing matched the closest gridpoint'
endif


deallocate( thisdist )

end subroutine find_closest_gridpoint


end program model_mod_check

! <next few lines under version control, do not edit>
! $URL: https://proxy.subversion.ucar.edu/DAReS/DART/releases/Lanai/models/template/model_mod_check.f90 $
! $Id: model_mod_check.f90 6311 2013-07-17 22:04:54Z thoar $
! $Revision: 6311 $
! $Date: 2013-07-18 00:04:54 +0200 (Thu, 18 Jul 2013) $
