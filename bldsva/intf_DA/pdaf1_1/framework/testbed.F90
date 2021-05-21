subroutine testbed

use clm4cmem, only: SATELLITE
implicit none

character(len=300) :: clmin, inpar
type(SATELLITE) :: SAT

clmin = "/homeb/hbn29/hbn29q/command/clmoas.clm2.h0.2015-07-01-00900.nc"
inpar = "input"
call cmem_main(trim(clmin),trim(inpar),SAT)

print*,SAT%name
print*,SAT%TBSAT_HV(1:10,1,1,1)
print*,size(SAT%TBSAT_HV)

END subroutine testbed
