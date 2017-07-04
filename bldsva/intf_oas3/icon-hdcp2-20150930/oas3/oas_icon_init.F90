SUBROUTINE oas_icon_init 

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
! 1.1.1        2017/07/ Slavko Brdar 
!   Modfied and Implemented in icon-hdcp2-20150930
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

USE oas_icon_vardef

!==============================================================================

IMPLICIT NONE

!==============================================================================

! Local Variable

CHARACTER(LEN=6)   :: MODNAME = 'oas_icon'    ! Name of the model

!------------------------------------------------------------------------------
!- End of header
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!- Begin Subroutine oas_icon_init
!------------------------------------------------------------------------------

 ! 1st Initialize the PRISM system for the component
 CALL prism_init_comp_proto( ncomp_id, MODNAME, nerror )
 IF( nerror /= PRISM_Success )   CALL prism_abort_proto( ncomp_id, 'oas_icon_init', 'Failure in prism_init_comp' )

 ! 2nd Get an MPI communicator fr CLM local communication
 CALL prism_get_localcomm_proto( kl_comm, nerror )
 IF( nerror /= PRISM_Success )   CALL prism_abort_proto( ncomp_id, 'oas_icon_init', 'Failure in prism_get_localcomm' )
      
 WRITE(6,*) "oas_icon: oas_icon_init: prism_init" 

!------------------------------------------------------------------------------
!- End of the Subroutine
!------------------------------------------------------------------------------

END SUBROUTINE oas_icon_init
