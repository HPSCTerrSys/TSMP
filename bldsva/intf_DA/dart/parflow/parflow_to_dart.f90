! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! DART $Id: parflow_to_dart.f90 Wed Apr 11 20:26:43 CEST 2018 $

program parflow_to_dart

!----------------------------------------------------------------------
! purpose: interface between parflow and DART
!
! method: Read pfb files of ParFlow model state.
!         Reform fields into a DART state vector (control vector).
!         Write out state vector in "proprietary" format for DART.
!         The output is a "DART restart file" format.
!
! 
! USAGE:  The parflow pfb filename is read from the input.nml
!         <edit pfb_to_dart_output_file in input.nml:pfb_to_dart_nml>
!         ./parflow_to_dart
!
! author: Tim Hoar 6/24/09
!         P. Shrestha 7/7/16 Update the template for ParFlow
!----------------------------------------------------------------------

use        types_mod, only : r8

use    utilities_mod, only : initialize_utilities, finalize_utilities, &
                             find_namelist_in_file, check_namelist_read, &
                             open_file, close_file,E_ERR, E_MSG, error_handler

use        model_mod, only : get_model_size, get_state_vector, &
                             write_state_times, get_parflow_filename, &
                             get_parflow_id, static_init_model

use  assim_model_mod, only : awrite_state_restart, open_restart_write, close_restart

use time_manager_mod, only : time_type, print_time, print_date

implicit none

! version controlled file description for error handling, do not edit
character(len=*), parameter :: source   = "$URL: parflow_to_dart.f90 $"
character(len=*), parameter :: revision = "$Revision: Bonn $"
character(len=*), parameter :: revdate  = "$Date: Wed Apr 11 2018 $"

!-----------------------------------------------------------------------
! namelist parameters with default values.
!-----------------------------------------------------------------------

logical            :: verbose = .TRUE.
character(len=128) :: pfb_to_dart_output_file  = 'dart_ics'

namelist /pfb_to_dart_nml/ pfb_to_dart_output_file, verbose

!----------------------------------------------------------------------
! global storage
!----------------------------------------------------------------------

character(len=512)    :: string1, string2
character(len=256)    :: parflow_filename
integer               :: io, iunit, x_size
type(time_type)       :: model_time
real(r8), allocatable :: x(:)       ! statevector
integer               :: pfid       ! ParFlow ID (PRESSURE_HEAD -or- SATURATION)

!======================================================================

call initialize_utilities(progname='parflow_to_dart')

call static_init_model() ! to initialize the model_mod

!----------------------------------------------------------------------
! Read the namelist to get the output filename.
!----------------------------------------------------------------------

call find_namelist_in_file("input.nml", "pfb_to_dart_nml", iunit)
read(iunit, nml = pfb_to_dart_nml, iostat = io)
call check_namelist_read(iunit, io, "pfb_to_dart_nml") ! closes, too.

pfid = get_parflow_id()

parflow_filename = get_parflow_filename(pfid)

write(string1,*) 'converting parflow file "'//trim(parflow_filename)//'"'
write(string2,*) ' to DART file "'//trim(pfb_to_dart_output_file)//'"'
call error_handler(E_MSG, 'parflow_to_dart', string1, text2=string2)

!----------------------------------------------------------------------
! get to work
!----------------------------------------------------------------------

x_size = get_model_size()
allocate( x(x_size) )

call get_state_vector(x, pfid, model_time)

iunit = open_restart_write(pfb_to_dart_output_file)

call awrite_state_restart(model_time, x, iunit)
call close_restart(iunit)

deallocate(x)

iunit = open_file('parflow_prior_time.txt',form='formatted')
call write_state_times(iunit, model_time)
call close_file(iunit)

!----------------------------------------------------------------------
! finish up
!----------------------------------------------------------------------

call print_date(model_time, str=' parflow_to_dart model date')
call print_time(model_time, str=' parflow_to_dart model time')
call finalize_utilities(progname='parflow_to_dart')

end program parflow_to_dart
