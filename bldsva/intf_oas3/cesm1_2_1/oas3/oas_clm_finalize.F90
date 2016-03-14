SUBROUTINE oas_clm_finalize

!---------------------------------------------------------------------
! Description:
!  This routine finalizes the coupling. If MPI_init has not been
!  called explicitly before oas_clm_init it will also close
!  MPI communication. 
!
! Current Code Owner: TR32, Z4: Prabhakar Shrestha
!    phone: 0228733453
!    email: pshrestha@uni-bonn.de
!
! History:
! Version    Date       Name
! ---------- ---------- ----
! 2.1.0        2016/02/29 Prabhakar Shrestha
! Implementation for CESM 1.2.1
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

!==============================================================================

IMPLICIT NONE

!==============================================================================

!------------------------------------------------------------------------------
!- End of header
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!- Begin Subroutine oas_clm_finalize 
!------------------------------------------------------------------------------

 DEALLOCATE(exfld)

 WRITE(6,*) "oasclm: oas_clm_finalize:  prsim terminate"
 CALL prism_terminate_proto ( nerror )

!------------------------------------------------------------------------------
!- End of the Subroutine
!------------------------------------------------------------------------------

END SUBROUTINE oas_clm_finalize
