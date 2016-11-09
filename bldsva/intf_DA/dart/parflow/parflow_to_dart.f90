! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id: model_to_dart.f90 6311 2013-07-17 22:04:54Z thoar $

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
!         <edit parflow_to_dart_output_file in input.nml:model_to_dart_nml>
!         model_to_dart
!
! author: Tim Hoar 6/24/09
!         P. Shrestha 7/7/16 Update the template for ParFlow
!----------------------------------------------------------------------

use        types_mod, only : r8

use    utilities_mod, only : initialize_utilities, finalize_utilities, &
                             find_namelist_in_file, check_namelist_read, &
                             open_file, close_file,E_ERR, E_MSG, error_handler

use        model_mod, only : get_model_size, get_state_vector, &
                             write_state_times,get_parflow_filename, static_init_model

use  assim_model_mod, only : awrite_state_restart, open_restart_write, close_restart

use time_manager_mod, only : time_type, print_time, print_date

implicit none

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL: https://proxy.subversion.ucar.edu/DAReS/DART/releases/Lanai/models/template/model_to_dart.f90 $"
character(len=32 ), parameter :: revision = "$Revision: 6311 $"
character(len=128), parameter :: revdate  = "$Date: 2013-07-18 00:04:54 +0200 (Thu, 18 Jul 2013) $"

!-----------------------------------------------------------------------
! namelist parameters with default values.
!-----------------------------------------------------------------------

character(len=128) :: pfb_to_dart_output_file  = 'dart_ics'

namelist /pfb_to_dart_nml/    &
     pfb_to_dart_output_file

!----------------------------------------------------------------------
! global storage
!----------------------------------------------------------------------

character(len=256)    :: parflow_filename
logical               :: verbose = .TRUE.
integer               :: io, iunit, x_size
type(time_type)       :: model_time
real(r8), allocatable :: x(:)               !statevector
character(len=512)   :: string1, string2

!======================================================================

call initialize_utilities(progname='parflow_to_dart')

!----------------------------------------------------------------------
! Read the namelist to get the output filename.
!----------------------------------------------------------------------

call find_namelist_in_file("input.nml", "pfb_to_dart_nml", iunit)
read(iunit, nml = pfb_to_dart_nml, iostat = io)
call check_namelist_read(iunit, io, "pfb_to_dart_nml") ! closes, too.

parflow_filename = get_parflow_filename()

write(string1,*) 'converting parflow file "'//trim(parflow_filename)//'"'
write(string2,*) ' to DART file "'//trim(pfb_to_dart_output_file)//'"'
call error_handler(E_MSG,'parflow_to_dart',string1,text2=string2)

!----------------------------------------------------------------------
! get to work
!----------------------------------------------------------------------

call static_init_model()

x_size = get_model_size()
allocate( x(x_size) )

call get_state_vector(x, model_time) 

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

call print_date(model_time, str='parflow_to_dart:model model date')
call print_time(model_time, str='parflow_to_dart:DART model time')
call finalize_utilities()

end program parflow_to_dart

! <next few lines under version control, do not edit>
! $URL: https://proxy.subversion.ucar.edu/DAReS/DART/releases/Lanai/models/template/model_to_dart.f90 $
! $Id: model_to_dart.f90 6311 2013-07-17 22:04:54Z thoar $
! $Revision: 6311 $
! $Date: 2013-07-18 00:04:54 +0200 (Thu, 18 Jul 2013) $
