! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id: column_rand.f90 6256 2013-06-12 16:19:10Z thoar $
 
program column_rand

! Allows creation of input file for generating a set of randomly located
! observation stations with full column of obs for b-grid model. Should be
! nearly identical to similar thing for CAM, etc.

use      types_mod, only : r8, PI
use random_seq_mod, only : random_seq_type, init_random_seq, random_uniform
use  utilities_mod, only : get_unit
use   location_mod, only : VERTISSURFACE, VERTISLEVEL

implicit none

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL: https://proxy.subversion.ucar.edu/DAReS/DART/releases/Lanai/models/bgrid_solo/column_rand.f90 $"
character(len=32 ), parameter :: revision = "$Revision: 6256 $"
character(len=128), parameter :: revdate  = "$Date: 2013-06-12 18:19:10 +0200 (Wed, 12 Jun 2013) $"

integer  :: level, num_cols, num_levs,num_dates, i,k, iunit
real(r8) :: lat, lon, t_err_var, uv_err_var, ps_err_var
real(r8) :: lat_data(10), lon_data(10), level_data(5)
integer  :: yyyy, mm, dd_data(6), hh, mn, ss

data  lat_data/49.93, 49.90, 49.87, 49.90, 49.93, 49.90, 49.87, 40.90, 49.93, 49.90/
data  lon_data/ 5.43,  5.45,  5.47,  5.49,  5.51,  5.53,  5.55,  5.57,  5.59,  5.60/
data  level_data/10.0, 100.0, 200.0, 500., 1000./
data  dd_data/8, 9, 10, 11, 12, 13/

!Pre Specified Date
yyyy = 2008
mm   = 5
hh   = 0
mn   = 0
ss   = 0
! Open an output file and write header info
iunit = get_unit()
open(unit = iunit, file = 'column_rand.out')

write(*, *) 'input the number of columns:10'
read(*, *) num_cols

write(*, *) 'input the number of model levels:5'
read(*, *) num_levs

write(*, *) 'input the number of dates:6'
read(*, *) num_dates

! Output the total number of obs
write(*, *) 'Max num is + 10 ', num_cols*(num_levs*num_dates)+10
write(iunit, *) num_cols*(num_levs*num_dates)+10     !CPS_write

! Get error variance for t
write(*, *) 'Input error VARIANCE for T obs'
read(*, *) t_err_var

! No values or qc
write(iunit, *) 0                          !CPS_write
write(iunit, *) 0                          !CPS_write

! Loop through each column
do k = 1, num_dates
do i = 1, num_cols
   ! Longitude 
   lon = lon_data(i) 

   ! Latitude 
   lat = lat_data(i)

   ! Loop through each observation in the column
   do level = 1, num_levs

      write(iunit, *) 0                   !CPS_write
      ! Write out the t observation
      write(iunit, *) 5                   !CPS_write 'RADIOSONDE_TEMPERATURE'
      write(iunit, *) 3                   !CPS_write height
      write(iunit, *) level_data(level)   !CPS_write
      write(iunit, *) lon                 !CPS_write
      write(iunit, *) lat                 !CPS_write
      write(iunit, *) yyyy,mm,dd_data(k),hh,mn,ss        !CPS_write
      write(iunit, *) t_err_var           !CPS_write

   end do
end do
end do

write(iunit, *) -1                        !CPS_write
write(iunit, *) 'set_def.out'             !CPS_write

end program column_rand

! <next few lines under version control, do not edit>
! $URL: https://proxy.subversion.ucar.edu/DAReS/DART/releases/Lanai/models/bgrid_solo/column_rand.f90 $
! $Id: column_rand.f90 6256 2013-06-12 16:19:10Z thoar $
! $Revision: 6256 $
! $Date: 2013-06-12 18:19:10 +0200 (Wed, 12 Jun 2013) $
