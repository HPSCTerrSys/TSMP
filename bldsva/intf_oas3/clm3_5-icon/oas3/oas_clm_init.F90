SUBROUTINE oas_clm_init 

!---------------------------------------------------------------------
! Description:
!  This routine initializes coupler to get the MPI communicator
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

!==============================================================================

IMPLICIT NONE

!==============================================================================

! Local Variable
   
CHARACTER(LEN=6)   :: CLMODNAME = 'oasclm'  ! Name of the mode

!------------------------------------------------------------------------------
!- End of header
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!- Begin Subroutine oas_clm_init 
!------------------------------------------------------------------------------


 ! 1st Initialize the PRISM system for the application
  CALL prism_init_comp_proto ( ncomp_id, CLMODNAME, nerror )
  IF ( nerror /= PRISM_Ok ) &
    CALL prism_abort_proto (ncomp_id, 'oas_clm_init', 'Failure in prism_init_comp_proto')

 ! 2nd Get an MPI communicator fr CLM local communication
  CALL prism_get_localcomm_proto ( kl_comm, nerror )
  IF ( nerror /= PRISM_Ok ) &
    CALL prism_abort_proto (ncomp_id, 'oas_clm_init','Failure in prism_get_localcomm_proto' )

  WRITE(6,*) "oasclm: oas_clm_init:  prism_init"

!------------------------------------------------------------------------------
!- End of the Subroutine
!------------------------------------------------------------------------------

END SUBROUTINE oas_clm_init
