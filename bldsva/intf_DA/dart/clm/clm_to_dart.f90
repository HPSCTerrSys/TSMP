! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id: clm_to_dart.f90 7195 2014-10-03 17:01:29Z thoar $

program clm_to_dart

!----------------------------------------------------------------------
! purpose: interface between clm and DART
!
! method: Read clm "restart" files of model state
!         Reform fields into a DART state vector (control vector).
!         Write out state vector in "proprietary" format for DART.
!         The output is a "DART restart file" format.
! 
! USAGE:  The clm filename is read from the clm_in namelist
!         <edit clm_to_dart_output_file in input.nml:clm_to_dart_nml>
!         clm_to_dart
!
! author: Tim Hoar 12 July 2011
!----------------------------------------------------------------------

use        types_mod, only : r8
use    utilities_mod, only : initialize_utilities, finalize_utilities, &
                             open_file, close_file, find_namelist_in_file, check_namelist_read
use        model_mod, only : get_model_size, clm_to_dart_state_vector, &
                             write_state_times
use  assim_model_mod, only : awrite_state_restart, open_restart_write, close_restart
use time_manager_mod, only : time_type, print_time, print_date

implicit none

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL: https://svn-dares-dart.cgd.ucar.edu/DART/trunk/models/clm/clm_to_dart.f90 $"
character(len=32 ), parameter :: revision = "$Revision: 7195 $"
character(len=128), parameter :: revdate  = "$Date: 2014-10-03 11:01:29 -0600 (Fri, 03 Oct 2014) $"

!-----------------------------------------------------------------------
! namelist parameters with default values.
!-----------------------------------------------------------------------

character(len=512) :: clm_to_dart_output_file  = 'dart_ics'

namelist /clm_to_dart_nml/ clm_to_dart_output_file

!----------------------------------------------------------------------
! global storage
!----------------------------------------------------------------------

integer               :: io, iunit, x_size
type(time_type)       :: model_time
real(r8), allocatable :: statevector(:)

!======================================================================

call initialize_utilities(progname='clm_to_dart')

!----------------------------------------------------------------------
! Read the namelist to get the output filename.
!----------------------------------------------------------------------

call find_namelist_in_file("input.nml", "clm_to_dart_nml", iunit)
read(iunit, nml = clm_to_dart_nml, iostat = io)
call check_namelist_read(iunit, io, "clm_to_dart_nml") ! closes, too.

!----------------------------------------------------------------------
! get to work
!----------------------------------------------------------------------

x_size = get_model_size()
allocate(statevector(x_size))

! Each variable specifies its own file of origin.
call clm_to_dart_state_vector(statevector, model_time) 

iunit = open_restart_write(clm_to_dart_output_file)

call awrite_state_restart(model_time, statevector, iunit)
call close_restart(iunit)

call print_date(model_time, str='clm_to_dart:clm  model date')
call print_time(model_time, str='clm_to_dart:DART model time')

!CPS need this with terrsysmp
iunit = open_file('clm_prior_time.txt',form='formatted')
call write_state_times(iunit, model_time)
call close_file(iunit)

call finalize_utilities('clm_to_dart')

end program clm_to_dart

! <next few lines under version control, do not edit>
! $URL: https://svn-dares-dart.cgd.ucar.edu/DART/trunk/models/clm/clm_to_dart.f90 $
! $Id: clm_to_dart.f90 7195 2014-10-03 17:01:29Z thoar $
! $Revision: 7195 $
! $Date: 2014-10-03 11:01:29 -0600 (Fri, 03 Oct 2014) $
