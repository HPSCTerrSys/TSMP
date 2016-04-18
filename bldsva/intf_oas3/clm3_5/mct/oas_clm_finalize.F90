SUBROUTINE oas_clm_finalize

!---------------------------------------------------------------------
! Description:
!  This routine finalizes the coupling. If MPI_init has not been
!  called explicitly before oas_clm_init it will also close
!  MPI communication. 
!
! References:
!  CEREFACS/ETH: E. Maisonnave, Edoward Davin
!
! Current Code Owner: TR32, Z4: Prabhakar Shrestha
!    phone: 0228733453
!    email: pshrestha@uni-bonn.de
!
! History:
! Version    Date       Name
! ---------- ---------- ----
! 1.1        2011/11/28 Prabhakar Shrestha 
!   Modfied and Implemented in CLM3.5, Initial release
! @VERSION@    @DATE@     <Your name>
!  <Modification comments>         
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!==============================================================================

! Declarations:
!
! Modules used:

USE oas_clm_vardef
USE spmdMod      ,       ONLY : masterproc

!==============================================================================

IMPLICIT NONE

!==============================================================================

!------------------------------------------------------------------------------
!- End of header
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!- Begin Subroutine oas_clm_finalize 
!------------------------------------------------------------------------------



 WRITE(6,*) "oasclm: oas_clm_finalize:  prsim terminate"
 CALL prism_terminate_proto ( nerror )

!------------------------------------------------------------------------------
!- End of the Subroutine
!------------------------------------------------------------------------------

END SUBROUTINE oas_clm_finalize
