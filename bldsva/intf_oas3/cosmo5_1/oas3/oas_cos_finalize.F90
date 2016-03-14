SUBROUTINE oas_cos_finalize

!---------------------------------------------------------------------
! Description:
!  This routine finalizes the coupling. If MPI_init has not been
!  called explicitly before oas_cos_init it will also close
!  MPI communication.
!
! Current Code Owner: TR32, Z4: Prabhakar Shrestha
!    phone: 0228733453
!    email: pshrestha@uni-bonn.de
!
! History:
! Version    Date       Name
! ---------- ---------- ----
! 1.1.1        2011/11/28 Prabhakar Shrestha 
!   Modfied and Implemented in COSMO4.11, Initial release
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

USE oas_cos_vardef

!==============================================================================

IMPLICIT NONE

!==============================================================================

!------------------------------------------------------------------------------
!- End of header
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!- Begin Subroutine oas_cos_finalize 
!------------------------------------------------------------------------------


 WRITE(6,*) 'oascos: oas_cos_finalize: prism_terminate '
 CALL flush(6)
 IF ( lpe_cpl ) DEALLOCATE(exfld)
 CALL prism_terminate_proto ( nerror )         

!------------------------------------------------------------------------------
!- End of the Subroutine
!------------------------------------------------------------------------------

END SUBROUTINE oas_cos_finalize
