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
                             print_time, write_time, operator(-)

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

character (len = 129) :: dart_input_file               = 'dart_ics'
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

type(time_type)       :: model_time, adv_to_time
real(r8), allocatable :: statevector(:)

integer :: i, skip

character(len=metadatalength) :: state_meta(1)
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

!CPSloc = set_location(loc_of_interest(1), loc_of_interest(2), loc_of_interest(3), VERTISHEIGHT)
!CPSmykindindex = get_raw_obs_kind_index(kind_of_interest)

if (test1thru < 1) goto 999

write(*,*)
write(*,*)'static_init_model test STARTING ...'
call static_init_model()
write(*,*)'static_init_model test COMPLETE ...'

if (test1thru < 2) goto 999

write(*,*)
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

statevector = 1.0_r8;
model_time  = set_time(21600, 149446)   ! 06Z 4 March 2010

!----------------------------------------------------------------------
! Open a test DART initial conditions file.
! Reads the valid time, the state, and (possibly) a target time.
!----------------------------------------------------------------------
if (test1thru < 4) goto 999

write(*,*)
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
write(*,*)'Exercising the netCDF routines.'
write(*,*)'Creating '//trim(output_file)//'.nc'

state_meta(1) = 'restart test'
ncFileID = init_diag_output(trim(output_file),'just testing a restart', 1, state_meta)

call aoutput_diagnostics(ncFileID, model_time, statevector, 1)

call nc_check( finalize_diag_output(ncFileID), 'model_mod_check:main', 'finalize')


!----------------------------------------------------------------------
! Checking get_state_meta_data (and get_state_indices, get_state_kind)
!----------------------------------------------------------------------

write(*,*)
write(*,*)'Checking metadata routines.'

if (test1thru < 6) goto 999

skip = 1000000

do i = 1, x_size, skip
   if ( i > 0 .and. i <= x_size ) call check_meta_data( i )
enddo

!----------------------------------------------------------------------
! Trying to find the state vector index closest to a particular ...
! Checking for valid input is tricky ... we don't know much. 
!----------------------------------------------------------------------

if (test1thru < 7) goto 999

if ( loc_of_interest(1) >= 0.0_r8 ) call find_closest_gridpoint( loc_of_interest )

!----------------------------------------------------------------------
! Check the interpolation - print initially to STDOUT
!----------------------------------------------------------------------

write(*,*)
write(*,*)'Testing model_interpolate ...'

call model_interpolate(statevector, loc, 1 , interp_val, ios_out)

if ( ios_out == 0 ) then 
   write(*,*)'model_interpolate SUCCESS: The interpolated value is ',interp_val
else
   write(*,*)'model_interpolate ERROR: model_interpolate failed with error code ',ios_out
endif

 999 continue

call finalize_utilities()

! end of main program

contains


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



subroutine find_closest_gridpoint( loc_of_interest )
! Simple exhaustive search to find the indices into the 
! state vector of a particular location.
real(r8), intent(in) :: loc_of_interest(:)

type(location_type) :: loc0, loc1
integer  :: i, indx(1)
real(r8) :: closest
character(len=129)  :: string1
real(r8), allocatable, dimension(:) :: thisdist

loc0 = set_location(loc_of_interest)

write(*,*)
write(*,'(''Checking for the index in the state vector that is closest to '')')
call write_location(0, loc0, charstring=string1)
write(*,*) trim(string1)

allocate( thisdist(get_model_size()) )
thisdist  = 9999999999.9_r8         ! really far away 



! Since there can be multiple variables with
! identical distances, we will just cruise once through 
! the array and come back to find all the 'identical' values.
do i = 1,get_model_size()

   ! Really inefficient, but grab the 'which_vert' from the
   ! grid and set our target location to have the same.
   ! Then, compute the distance and compare.

   call get_state_meta_data(i, loc1)
   thisdist(i) = get_dist( loc1, loc0)

enddo

indx = minloc(thisdist)

! Now that we know  ... report 

write(*, *) 'closest to the given location is index ', indx(1)

deallocate( thisdist )

end subroutine find_closest_gridpoint


end program model_mod_check

! <next few lines under version control, do not edit>
! $URL: https://proxy.subversion.ucar.edu/DAReS/DART/releases/Lanai/models/template/model_mod_check.f90 $
! $Id: model_mod_check.f90 6311 2013-07-17 22:04:54Z thoar $
! $Revision: 6311 $
! $Date: 2013-07-18 00:04:54 +0200 (Thu, 18 Jul 2013) $
