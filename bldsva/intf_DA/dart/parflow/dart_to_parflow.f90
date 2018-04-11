! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! DART $Id: dart_to_parflow.f90 Wed Apr 11 20:26:43 CEST 2018 $

program dart_to_parflow

!----------------------------------------------------------------------
! method: Read DART state vector and overwrite values in the parflow restart file.
!         If the DART state vector has an 'advance_to_time' present, a
!         file called "dart_posterior_times.txt" is created to help 
!         parflow to advance to the requested time.
!
!         The dart_to_pfb_nml namelist setting for advance_time_present 
!         determines whether or not the input file has an 'advance_to_time'.
!         Typically, only temporary files like 'assim_model_state_ic' have
!         an 'advance_to_time'.
!----------------------------------------------------------------------

use        types_mod, only : r8

use    utilities_mod, only : initialize_utilities, finalize_utilities, &
                             find_namelist_in_file, check_namelist_read, &
                             logfileunit, open_file, close_file, E_ERR, E_MSG, &
                             error_handler

use  assim_model_mod, only : open_restart_read, aread_state_restart, close_restart

use time_manager_mod, only : time_type, print_time, print_date, operator(-), &
                             get_time, get_date

use        model_mod, only : static_init_model, write_parflow_file, &
                             get_parflow_id, PRESSURE_HEAD, SATURATION,  &
                             get_model_size, get_state_vector, write_state_times

implicit none

! version controlled file description for error handling, do not edit
character(len=*), parameter :: source   = "$URL: dart_to_parflow.f90 $"
character(len=*), parameter :: revision = "$Revision: Bonn $"
character(len=*), parameter :: revdate  = "$Date: Wed Apr 11 2018 $"

character(len=512) :: string1, string2, string3

!------------------------------------------------------------------
! The namelist variables
!------------------------------------------------------------------

character(len=256) :: pfb_restart_filename   = 'model_restartfile'
character(len=256) :: dart_to_pfb_input_file = 'dart_restart'
logical            :: advance_time_present   = .false.
logical            :: verbose                = .false.

namelist /dart_to_pfb_nml/  dart_to_pfb_input_file, &
                            advance_time_present,   &
                            pfb_restart_filename,   &
                            verbose

!----------------------------------------------------------------------

integer               :: iunit, io, x_size
type(time_type)       :: model_time, adv_to_time, base_time
real(r8), allocatable :: statevector(:)
real(r8), allocatable :: pfb_state(:)
integer               :: pfid, jpfb

!----------------------------------------------------------------------

call initialize_utilities(progname='dart_to_parflow')

! model_mod:static_init_model() reads the model namelist and set grid sizes, etc.

call static_init_model()

x_size = get_model_size()
allocate(statevector(1:x_size))

! Read the namelist to get the input filename. 

call find_namelist_in_file("input.nml", "dart_to_pfb_nml", iunit)
read(iunit, nml = dart_to_pfb_nml, iostat = io)
call check_namelist_read(iunit, io, "dart_to_pfb_nml")

!----------------------------------------------------------------------
! Reads the valid time, the state, and the target time.
!----------------------------------------------------------------------

iunit = open_restart_read(dart_to_pfb_input_file)

if ( advance_time_present ) then
   call aread_state_restart(model_time, statevector, iunit, adv_to_time)
else
   call aread_state_restart(model_time, statevector, iunit)
endif
call close_restart(iunit)

call print_date( model_time,' dart_to_parflow model date')
call print_time( model_time,' dart_to_parflow model time')
call print_date( model_time,' dart_to_parflow model date',logfileunit)
call print_time( model_time,' dart_to_parflow model time',logfileunit)

!----------------------------------------------------------------------
! update the current model state vector
! Convey the amount of time to integrate the model ...
! time_manager_nml: stop_option, stop_count increments
!----------------------------------------------------------------------

allocate( pfb_state(1:x_size) )

pfid = get_parflow_id()

!>@todo I do not understand why we switch ... don't we always update the pressure
if (pfid .eq. PRESSURE_HEAD) then
   jpfb               = SATURATION 
   write(string3,'(A)')'SATURATION'
elseif (pfid .eq. SATURATION) then
   jpfb               = PRESSURE_HEAD
   write(string3,'(A)')'PRESSURE_HEAD'
endif

write(string1,'(A)')'converting DART file "'//trim(dart_to_pfb_input_file)//'"'
write(string2,'(A)')'using parflow state "'//trim(string3)//'"'
write(string3,'(A)')'updating "'//trim(pfb_restart_filename)//'"'
call error_handler(E_MSG,'dart_to_parflow',string1,text2=string2,text3=string3)

call get_state_vector(pfb_state, jpfb)

call write_parflow_file(statevector, pfb_state, &
             dart_to_pfb_input_file, pfb_restart_filename)

!>@todo ... dart_posterior_times.txt and dart_prior_time.txt will always
!>          have the same information ...

iunit = open_file('dart_posterior_times.txt', action='write')
call write_state_times(iunit, model_time)
call close_file(iunit)

!----------------------------------------------------------------------
! Log what we think we're doing, and exit.
!----------------------------------------------------------------------

!>@ the adv_to_time is never used ... we can skip this if need be ...

if ( advance_time_present ) then
call print_time(adv_to_time,' dart_to_parflow advance_to time')
call print_date(adv_to_time,' dart_to_parflow advance_to date')
call print_time(adv_to_time,' dart_to_parflow advance_to time',logfileunit)
call print_date(adv_to_time,' dart_to_parflow advance_to date',logfileunit)
endif

call finalize_utilities(progname='dart_to_parflow')

end program dart_to_parflow
