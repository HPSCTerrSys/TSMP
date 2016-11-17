! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id: dart_to_model.f90 6311 2013-07-17 22:04:54Z thoar $

program dart_to_parflow

!----------------------------------------------------------------------
! purpose: interface between DART and the model model
!
! method: Read DART state vector and overwrite values in a model restart file.
!         If the DART state vector has an 'advance_to_time' present, a
!         file called model_in.DART is created with a time_manager_nml namelist 
!         appropriate to advance model to the requested time.
!
!         The dart_to_pfb_nml namelist setting for advance_time_present 
!         determines whether or not the input file has an 'advance_to_time'.
!         Typically, only temporary files like 'assim_model_state_ic' have
!         an 'advance_to_time'.
!
! author: Tim Hoar 25 Jun 09, revised 12 July 2010
!----------------------------------------------------------------------

use        types_mod, only : r8
use    utilities_mod, only : initialize_utilities, finalize_utilities, &
                             find_namelist_in_file, check_namelist_read, &
                             logfileunit, open_file, close_file
use  assim_model_mod, only : open_restart_read, aread_state_restart, close_restart
use time_manager_mod, only : time_type, print_time, print_date, operator(-), &
                             get_time, get_date
use        model_mod, only : static_init_model, write_parflow_file, &
                             get_model_size, get_state_vector, write_state_times

implicit none

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL: https://proxy.subversion.ucar.edu/DAReS/DART/releases/Lanai/models/template/dart_to_model.f90 $"
character(len=32 ), parameter :: revision = "$Revision: 6311 $"
character(len=128), parameter :: revdate  = "$Date: 2013-07-18 00:04:54 +0200 (Thu, 18 Jul 2013) $"

!------------------------------------------------------------------
! The namelist variables
!------------------------------------------------------------------

character (len = 128) :: dart_to_pfb_input_file = 'dart_restart'
logical               :: advance_time_present     = .false.
character(len=256)    :: pfb_restart_filename   = 'model_restartfile'

namelist /dart_to_pfb_nml/  dart_to_pfb_input_file, &
                            advance_time_present,    &
                            pfb_restart_filename

!----------------------------------------------------------------------

integer               :: iunit, io, x_size, diff1, diff2
type(time_type)       :: model_time, adv_to_time, base_time
real(r8), allocatable :: statevector(:)
real(r8), allocatable :: sv_sat(:)         !diagnostic state vector
logical               :: verbose              = .FALSE.

!----------------------------------------------------------------------

call initialize_utilities(progname='dart_to_parflow', output_flag=verbose)

!----------------------------------------------------------------------
! Call model_mod:static_init_model() which reads the model namelists
! to set grid sizes, etc.
!----------------------------------------------------------------------

call static_init_model()

x_size = get_model_size()
allocate(statevector(1:x_size))

! Read the namelist to get the input filename. 

call find_namelist_in_file("input.nml", "dart_to_pfb_nml", iunit)
read(iunit, nml = dart_to_pfb_nml, iostat = io)
call check_namelist_read(iunit, io, "dart_to_pfb_nml")

write(*,*)
write(*,*) 'dart_to_parflow: converting DART file ', "'"//trim(dart_to_pfb_input_file)//"'"
write(*,*) 'to model restart files named        ', "'"//trim(pfb_restart_filename)//"'" 

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

!----------------------------------------------------------------------
! update the current model state vector
! Convey the amount of time to integrate the model ...
! time_manager_nml: stop_option, stop_count increments
!----------------------------------------------------------------------

write(*,*) 'Obtaining diagnostic state vector for clamping'
allocate( sv_sat(1:x_size) )
call get_state_vector(sv_sat, 2 )         ! For clamping

write(*,*) 'Calling write_parflow_file to restart file'
call write_parflow_file(statevector,sv_sat, dart_to_pfb_input_file, pfb_restart_filename)

iunit = open_file('dart_posterior_times.txt', action='write')
call write_state_times(iunit, model_time)
call close_file(iunit)

!----------------------------------------------------------------------
! Log what we think we're doing, and exit.
!----------------------------------------------------------------------

call print_date( model_time,'dart_to_parflow:model model date')
call print_time( model_time,'dart_to_parflow:DART model time')
call print_date( model_time,'dart_to_parflow:model model date',logfileunit)
call print_time( model_time,'dart_to_parflow:DART model time',logfileunit)

if ( advance_time_present ) then
call print_time(adv_to_time,'dart_to_parflow:advance_to time')
call print_date(adv_to_time,'dart_to_parflow:advance_to date')
call print_time(adv_to_time,'dart_to_parflow:advance_to time',logfileunit)
call print_date(adv_to_time,'dart_to_parflow:advance_to date',logfileunit)
endif

call finalize_utilities()

end program dart_to_parflow

! <next few lines under version control, do not edit>
! $URL: https://proxy.subversion.ucar.edu/DAReS/DART/releases/Lanai/models/template/dart_to_model.f90 $
! $Id: dart_to_model.f90 6311 2013-07-17 22:04:54Z thoar $
! $Revision: 6311 $
! $Date: 2013-07-18 00:04:54 +0200 (Thu, 18 Jul 2013) $
