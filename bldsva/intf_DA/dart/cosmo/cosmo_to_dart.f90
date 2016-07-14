! DART software - Copyright 2004 - 2011 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

program cosmo_to_dart

! DART $Id: cosmo_to_dart.f90 $

!----------------------------------------------------------------------
! purpose: interface between cosmo and DART
!
! method: Read grib files of cosmo model state.
!         Reform fields into a DART state vector.
!         Write out state vector in "proprietary" format for DART.
!         The output is a "DART restart file" format.
!         The cosmo filename is read from input.nml:model_nml:cosmo_filename
!         DART filename is specified in 
! 
! USAGE:  
!         cosmo_to_dart
!
!----------------------------------------------------------------------

use        types_mod, only : r8

use    utilities_mod, only : initialize_utilities, finalize_utilities, &
                             find_namelist_in_file, check_namelist_read, &
                             open_file, E_MSG, E_ERR, error_handler

use        model_mod, only : get_model_size, get_state_vector, &
                             get_cosmo_filename, static_init_model, &
                             write_state_times

use  assim_model_mod, only : awrite_state_restart, open_restart_write, close_restart

use time_manager_mod, only : time_type, print_time, print_date

implicit none

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = "$URL: cosmo_to_dart.f90 $"
character(len=32 ), parameter :: revision = "$Revision: none $"
character(len=128), parameter :: revdate  = "$Date: none $"

!-----------------------------------------------------------------------
! namelist parameters with default values.
!-----------------------------------------------------------------------

character(len=256) :: cosmo_to_dart_output_file  = 'dart.ud'

namelist /cosmo_to_dart_nml/ cosmo_to_dart_output_file

!----------------------------------------------------------------------
! global storage
!----------------------------------------------------------------------

integer              :: io, iunit, x_size
real(r8),allocatable :: x(:)
type(time_type)      :: model_time
character(len=256)   :: cosmo_filename
character(len=512)   :: string1, string2, string3

!======================================================================

call initialize_utilities(progname='cosmo_to_dart')

!----------------------------------------------------------------------
! Read the namelist to get the output filename.
!----------------------------------------------------------------------

call find_namelist_in_file("input.nml", "cosmo_to_dart_nml", iunit)
read(iunit, nml = cosmo_to_dart_nml, iostat = io)
call check_namelist_read(iunit, io, "cosmo_to_dart_nml") ! closes, too.

!----------------------------------------------------------------------
! Call model_mod:static_init_model() which reads the model namelists
! to set grid sizes, etc.
!----------------------------------------------------------------------

call static_init_model()

cosmo_filename = get_cosmo_filename(filetype='restart')

write(string1,*)'converting cosmo file "'//trim(cosmo_filename)//'"'
write(string2,*)' to DART file "'//trim(cosmo_to_dart_output_file)//'"'
call error_handler(E_MSG,'cosmo_to_dart',string1,text2=string2)

x_size = get_model_size()

allocate( x(x_size) )

call get_state_vector(x, model_time)

iunit = open_restart_write(cosmo_to_dart_output_file)

call awrite_state_restart(model_time, x, iunit)
call close_restart(iunit)

deallocate(x)

iunit = open_file('cosmo_prior_time.txt',form='formatted')
call write_state_times(iunit, model_time)

call print_date(model_time, str='cosmo_to_dart:cosmo model date')
call print_time(model_time, str='cosmo_to_dart:DART  model time')
call finalize_utilities()

end program cosmo_to_dart

! <next few lines under version control, do not edit>
! $URL: $
! $Id: cosmo_to_dart.f90 $
! $Revision: $
! $Date: $
