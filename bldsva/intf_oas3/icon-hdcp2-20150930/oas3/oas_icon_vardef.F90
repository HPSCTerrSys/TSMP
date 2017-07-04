MODULE oas_icon_vardef

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
! 1.2.1        2012/10/15 Markus Uebel, Prabhakar Shrestha
!   Modfied and Implemented in COSMO4.21
!              2013/09/17 P. Shrestha
!   Added cpl_scheme to TAB to choose COSMO-CLM coupling scheme
!   Two Options: True, transfer coeffiecent sent from CLM direclty
!                False, transfer coefficient inverted from CLM fluxes
! 1.2.2        2017/06/03 Slavko Brdar
!   Modfied and Implemented in icon-hdcp2-20150930
!
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

USE mo_kind,            ONLY:  wp
!USE mpi
USE mod_prism_proto              ! OASIS3 prism module

!==============================================================================

IMPLICIT NONE

!==============================================================================

SAVE

! Debug level of OASIS
!     0 : Minimum debugging
!     1 : Debugging
!     2 : Perfs measurement
!     3 : OASIS restart production

! Variables

INTEGER :: IOASISDEBUGLVL = 0 

! Variable ids

INTEGER            :: ncomp_id           ! id returned by prism_init_comp
INTEGER            :: kl_comm            ! Local communicator 
INTEGER            :: nerror             ! return error code
INTEGER, PUBLIC    :: OASIS_Rcv  = 1     ! return code if received field
INTEGER, PUBLIC    :: OASIS_idle = 0     ! return code if nothing done by oasis
INTEGER, PUBLIC                    :: PRISM_Success = 0  ! return code if no error in oasis

INTEGER, PARAMETER :: nmaxfld=40    ! Maximum number of coupling fields
INTEGER            :: ksnd=17,     &! Number of send coupling fields                !MU: changed to 17 due to CO2 coupling
                                      krcv=15       ! Number of received coupling fields            !MU: changed to 15 due to CO2 coupling

TYPE, PUBLIC     ::   FLD_CPL                 ! Type for coupling field information
  LOGICAL                 ::        laction   ! To be coupled or not
  CHARACTER(len = 8)      ::        clname    ! Name of the coupling field
  CHARACTER(len = 1)      ::        clgrid    ! Grid type
  INTEGER(KIND=iintegers) ::        nid       ! Id of the field
END TYPE FLD_CPL

TYPE(FLD_CPL), DIMENSION(nmaxfld), PUBLIC      :: srcv, ssnd   ! Coupling fields

REAL(KIND=wp), DIMENSION(:,:), ALLOCATABLE :: exfld        ! Temporary buffer for receiving
REAL(KIND=wp), DIMENSION(:,:), ALLOCATABLE :: zmask        ! CPS moved here

INTEGER                        :: nldi,nlei, nldj,nlej ! halo limits on a local subdomain
INTEGER                        :: jih, jjh             ! subdomain limit

LOGICAL :: lpe_cpl = .FALSE.
LOGICAL :: cpl_scheme                          !Coupling Scheme with CLM, now set in oas_icon_define  
LOGICAL :: fhalo =.True.                       !Removes Infinity with inversion of tcm/tch 1/0. 

REAL(KIND=wp), DIMENSION(:,:,:),ALLOCATABLE  ::   frcv        ! all fields recieved from soil model

INTEGER                        :: cplfreq         ! Coupling frequency of COSMO with clm
END MODULE oas_icon_vardef
