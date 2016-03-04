MODULE oas_clm_vardef

!---------------------------------------------------------------------
! Description:
!  Definition and variables for OASIS3 communications
!
! Current Code Owner: TR32, Z4: Prabhakar Shrestha
!    phone: 0228733453
!    email: pshrestha@uni-bonn.de
!
! History:
! Version    Date       Name
! ---------- ---------- ----
! 1.1.1        2011/11/28 Prabhakar Shrestha 
!   Modfied and Implemented in CLM3.5, Initial release
! 1.2.1        2012/09/13 Markus Uebel
!   Modified due to CO2 coupling
!   Added cpl_scheme to TAB to choose COSMO-CLM coupling scheme
!   Two Options: True, transfer coeffiecent sent from CLM direclty
!                False, transfer coefficient inverted from CLM fluxes
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

USE shr_kind_mod , only : r8 => shr_kind_r8
USE mpi
USE mod_prism_proto 

!==============================================================================

IMPLICIT NONE

!==============================================================================

SAVE

! Debug level of OASIS
!     0 : Minimum debugging
!     1 : Debugging
!     2 : Perfs measurement


! Variables

INTEGER            :: IOASISDEBUGLVL = 0 
! Variable ids
INTEGER            :: ncomp_id           ! id returned by prism_init_comp
INTEGER            :: kl_comm            ! Local communicator 
INTEGER            :: nerror             ! return error code
INTEGER, PUBLIC    :: OASIS_Rcv  = 1     ! return code if received field
INTEGER, PUBLIC    :: OASIS_idle = 0     ! return code if nothing done by oasis
INTEGER, PUBLIC    :: PRISM_Success = 0  ! return code if no error in oasis

INTEGER, PARAMETER :: nmaxfld=200          ! Maximum number of coupling fields
INTEGER            :: ksnd=15,  krcv=17    ! Number of send/received coupling fields between COSMO and CLM      !MU: changed to ksnd=15, krcv=17 for CO2 coupling
INTEGER            :: vsnd=10, vrcv=20    ! Number of send/received coupling fields between CLM and ParFlow


TYPE, PUBLIC       ::   FLD_CPL          ! Type for coupling field information
    LOGICAL            ::   laction   ! To be coupled or not
    CHARACTER(len =12) ::   clname    ! Name of the coupling field
    CHARACTER(len =3)  ::   ref       ! Type of the coupling field
    CHARACTER(len = 1) ::   clgrid    ! Grid type
    INTEGER            ::   nid       ! Id of the field
    INTEGER            ::   level     ! # of soil layer
END TYPE FLD_CPL

TYPE(FLD_CPL), DIMENSION(nmaxfld), PUBLIC    :: srcv, ssnd   ! Coupling fields

INTEGER                :: ndlon, ndlat

LOGICAL, ALLOCATABLE   :: llmask(:,:,:)
LOGICAL                :: cpl_scheme = .True.                 !Coupling Scheme with COSMO, False for inversion technique
INTEGER                :: start1d, length1d, pe_loc1d     ! Parallel Decomp

END MODULE oas_clm_vardef
