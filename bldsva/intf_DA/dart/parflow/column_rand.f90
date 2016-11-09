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
real(r8) :: lat, lon, sw_err_var, uv_err_var, ps_err_var
real(r8) :: lat_data(10), lon_data(10), level_data(5)
integer  :: yyyy, mm, dd_data(13), hh, mn, ss
character(len=1024) :: filename

data  lat_data/49.87, 49.95, 49.87, 49.95, 49.87, 49.95, 49.87, 49.95, 49.87, 49.95/
data  lon_data/ 5.48,  5.48,  5.55,  5.55,  5.62,  5.62,  5.69,  5.69,  5.75,  5.75/
data  level_data/0.06, 0.13, 0.54 /    !depth from surface [meters]
data  dd_data/9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21/

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

write(*, *) 'input the number of model levels:3'
read(*, *) num_levs

write(*, *) 'input the number of dates:13'
read(*, *) num_dates

! Output the total number of obs
write(*, *) 'Max num is + 10 ', num_cols*(num_levs*num_dates)+10
write(iunit, *) num_cols*(num_levs*num_dates)+10 !CPS_write upperbound on obs

! Get error variance for t
write(*, *) 'Input error VARIANCE for soil moisture obs'
read(*, *) sw_err_var

! No values or qc
write(iunit, *) 0                          !CPS_write number of copies
write(iunit, *) 0                          !CPS_write number of qc flags

! Loop through each column
do k = 1, num_dates
do i = 1, num_cols
   ! Longitude 
   lon = lon_data(i) 

   ! Latitude 
   lat = lat_data(i)

   ! Loop through each observation in the column
   do level = 1, num_levs

      write(iunit, *) 0                   !CPS_write -1 if no obs
      ! Write out the sw observation
      write(iunit, *) 3                   !CPS_write 'SOIL_MOISTURE'
      write(iunit, *) 3                   !CPS_write vertical co-ordinate 
      write(iunit, *) level_data(level)   !CPS_write height
      write(iunit, *) lon                 !CPS_write lon
      write(iunit, *) lat                 !CPS_write lat
      write(iunit, *) yyyy,mm,dd_data(k),hh,mn,ss        !CPS_write date
      write(iunit, *) sw_err_var           !CPS_write error variance

   end do
end do
end do

write(iunit, *) -1                        !CPS_write -1 if no obs
write(iunit, *) 'set_def.out'             !CPS_write filename for sequence
close(iunit)

!CPS added for export
do k = 1, num_dates
  iunit = get_unit()
  if (k .le. 9) then
    write (filename, "(A19,I1)") "column_export_out_0",k 
  else
    write (filename, "(A18,I2)") "column_export_out_",k
  endif
  open(unit = iunit, file = trim(filename))
  write(iunit, *) 'set_def.out'             !CPS_write
  write(iunit, *) 1                         !CPS_write
  write(iunit, *) 1                         !CPS_write
  write(iunit, *) yyyy,mm,dd_data(k),hh,mn,ss        !CPS_write
  write(iunit, *) 0,0                         !CPS_write
  write(iunit, *) 'obs_seq.in'             !CPS_write
  close(iunit)
end do

end program column_rand

! <next few lines under version control, do not edit>
! $URL: https://proxy.subversion.ucar.edu/DAReS/DART/releases/Lanai/models/bgrid_solo/column_rand.f90 $
! $Id: column_rand.f90 6256 2013-06-12 16:19:10Z thoar $
! $Revision: 6256 $
! $Date: 2013-06-12 18:19:10 +0200 (Wed, 12 Jun 2013) $
