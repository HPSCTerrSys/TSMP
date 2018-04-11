! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! DART $Id: column_rand.f90 Wed Apr 11 20:26:43 CEST 2018 $
 
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
character(len=*), parameter :: source   = & "$URL: column_rand.f90 $"
character(len=*), parameter :: revision = "$Revision: Bonn $"
character(len=*), parameter :: revdate  = "$Date: Wed Apr 11 2018 $"

integer  :: level, num_cols, num_levs,num_dates, i,k, iunit
real(r8) :: lat, lon, t_err_var, uv_err_var, ps_err_var
real(r8) :: lat_data(10), lon_data(10), level_data(7)
integer  :: yyyy, mm, dd_data(13), hh, mn, ss
character(len=512) :: filename

!This is based in CLM grids
data  lat_data/49.87, 49.92, 49.87, 49.92, 49.87, 49.92, 49.87, 49.92, 49.87, 49.92/
data  lon_data/ 5.45,  5.45,  5.48,  5.48,  5.51,  5.51,  5.55,  5.55, 5.59,  5.59/
data  level_data/ 0.02, 0.06, 0.10, 0.20, 0.30, 0.50, 0.80 /
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

write(*, *) 'input the number of model levels:7'
read(*, *) num_levs

write(*, *) 'input the number of dates:13'
read(*, *) num_dates

! Output the total number of obs
write(*, *) 'Max num is + 10 ', num_cols*(num_levs*num_dates)+10
write(iunit, *) num_cols*(num_levs*num_dates)+10     !CPS_write

! Get error variance for t
write(*, *) 'Input error VARIANCE for T soil obs'
read(*, *) t_err_var

write(iunit, *) 0  ! number of copies of data (0 for just a definition)  
write(iunit, *) 0  ! number of quality control values per field (0 or greater)

! Loop through each column
do k = 1, num_dates
do i = 1, num_cols
   ! Longitude 
   lon = lon_data(i) 

   ! Latitude 
   lat = lat_data(i)

   ! Loop through each observation in the column
   do level = 1, num_levs

      write(iunit, *) 0                   ! input a -1 if there are no more obs
      write(iunit, *) 2                   ! 'SOIL_TEMPERATURE' - maybe
      write(iunit, *) 3                   ! 3  --> height
      write(iunit, *) level_data(level)
      write(iunit, *) lon
      write(iunit, *) lat
      write(iunit, *) yyyy,mm,dd_data(k),hh,mn,ss ! TJH: this time will be overwritten in next phase
      write(iunit, *) t_err_var           ! observation error variance

   end do
end do
end do

write(iunit, *) -1                        ! there are no more obs
write(iunit, *) 'set_def.out'             ! name of template for next phase
close(iunit)

!CPS added for export
do k = 1, num_dates
  iunit = get_unit()
  write (filename, "(A18,I2.2)") "column_export_out_",k
  open(unit = iunit, file = trim(filename))
  write(iunit, *) 'set_def.out'
  write(iunit, *) 1
  write(iunit, *) 1
  write(iunit, *) yyyy,mm,dd_data(k),hh,mn,ss
  write(iunit, *) 0,0
  write(iunit, *) 'obs_seq.in'
  close(iunit)
end do

end program column_rand
