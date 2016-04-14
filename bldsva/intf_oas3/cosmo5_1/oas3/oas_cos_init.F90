SUBROUTINE oas_cos_init 

!---------------------------------------------------------------------
! Description:
!  This routine initializes coupler to get the MPI communicator
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

! Local Variable

CHARACTER(LEN=6)   :: MODNAME = 'oascos'    ! Name of the model

!------------------------------------------------------------------------------
!- End of header
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!- Begin Subroutine oas_cos_init 
!------------------------------------------------------------------------------

 ! 1st Initialize the PRISM system for the component
 CALL prism_init_comp_proto( ncomp_id, MODNAME, nerror )
 IF( nerror /= PRISM_Success )   CALL prism_abort_proto( ncomp_id, 'oas_cos_init', 'Failure in prism_init_comp' )

 ! 2nd Get an MPI communicator fr CLM local communication
 CALL prism_get_localcomm_proto( kl_comm, nerror )
 IF( nerror /= PRISM_Success )   CALL prism_abort_proto( ncomp_id, 'oas_cos_init', 'Failure in prism_get_localcomm' )
      
 WRITE(6,*) "oascos: oas_cos_init: prism_init" 

!------------------------------------------------------------------------------
!- End of the Subroutine
!------------------------------------------------------------------------------

END SUBROUTINE oas_cos_init
