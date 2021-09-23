! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id: column_rand.f90 Wed Apr 11 20:26:43 CEST 2018 $
 
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
character(len=*), parameter :: source   = "$URL: column_rand.f90 $"
character(len=*), parameter :: revision = "$Revision: Bonn $"
character(len=*), parameter :: revdate  = "$Date: Wed Apr 11 2018 $"

integer  :: level, num_cols, num_levs,num_dates, i,k, iunit
real(r8) :: lat, lon, t_err_var, uv_err_var, ps_err_var
real(r8) :: lat_data(10), lon_data(10), level_data(7)
integer  :: yyyy, daysinmonth(12), mm(60), dd_data(60), hh, mn, ss
character(len=1024) :: filename

data  lat_data/49.87, 49.95, 49.87, 49.95, 49.87, 49.95, 49.87, 49.95, 49.87, 49.95/
data  lon_data/ 5.48,  5.48,  5.55,  5.55,  5.62,  5.62,  5.69,  5.69,  5.75,  5.75/
data  level_data/10.0, 100.0, 200.0, 500., 1000.,3000.,5000./
data daysinmonth/31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/
!data  dd_data/ 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 24, &
!              25, 26, 27, 28, 29, 30, 31,  1,  2,  3,  4,  5,  6,  7,  8/

!Pre Specified Date
yyyy = 2008
hh   = 0
mn   = 0
ss   = 0
mm(1)      = 6
dd_data(1) = 7
! Open an output file and write header info
iunit = get_unit()
open(unit = iunit, file = 'column_rand.out')

write(*, *) 'input the number of columns:10'
read(*, *) num_cols

write(*, *) 'input the number of model levels:7'
read(*, *) num_levs

write(*, *) 'input the number of dates:30'
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
! Allows day and month update for a year only
if (k .gt. 1) then
  dd_data(k) = dd_data(k-1) + 1
  if (dd_data(k) .le. daysinmonth(mm(k-1))) then
    mm(k) = mm(k-1)
  else
    dd_data(k) = 1    !start of new month
    mm(k) = mm(k-1) + 1
  end if
end if
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
      write(iunit, *) yyyy,mm(k),dd_data(k),hh,mn,ss        !CPS_write
      write(iunit, *) t_err_var           !CPS_write

   end do
end do
end do

write(iunit, *) -1                        !CPS_write
write(iunit, *) 'set_def.out'             !CPS_write
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
  write(iunit, *) yyyy,mm(k),dd_data(k),hh,mn,ss        !CPS_write
  write(iunit, *) 0,0                         !CPS_write
  write(iunit, *) 'obs_seq.in'             !CPS_write
  close(iunit)
end do

end program column_rand
