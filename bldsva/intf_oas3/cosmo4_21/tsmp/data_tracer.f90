! $RCSfile: data_tracer.f90,v $
! $Revision: 4.11 $ $Date: 2011/06/27
! Data module for all global meteorological fields used for tracer transport
!------------------------------------------------------------------------------

MODULE data_tracer

!------------------------------------------------------------------------------
!
! Description:
! ===========
!  This module declares all meteorological fields used for tracer transport
!  that have to reside in the long term storage, i.e. that are used in
!  more than one module.
!
!  All fields are declared as allocatable arrays. They are allocated in the
!  setup of the model and deallocated in the cleanup at the end of the
!  program.
!
! Current Code Owner: Markus Uebel, MIUB
!  phone:  +49  228 73 5186
!  fax:    +49  228 73 5188
!  email:  muebel@uni-bonn.de
!
! History:
! Version    Date       Name
! ---------- ---------- ----
! V_4.11     2011/12/12 Markus Uebel
!
! @VERSION@    @DATE@     <Your name>
!  <Modification comments>         
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!==============================================================================
!
! Declarations:
! -------------
!
USE data_parameters, ONLY :   &
    ireals,    & ! KIND-type parameter for real variables
    iintegers    ! KIND-type parameter for standard integer variables

!==============================================================================

IMPLICIT NONE

!==============================================================================

! Global (i.e. public) Declarations:

  REAL  (KIND=ireals), TARGET, ALLOCATABLE ::           &
! Note: 4 dim. variables: var1(i,j,k,ntracer) 
!       5 dim. variables: var2(i,j,k,nstep,ntracer)
    tracer(:,:,:,:,:),&       ! tracer concentration                  ("var")
    tracertens(:,:,:,:),&     ! tracer tendency                       ("var"/s)
    tracer_bd(:,:,:,:,:),&    ! boundary field for tracer             ("var")
    tracert_conv(:,:,:,:),&   ! tracer tendency due to convection     ("var"/s)
    tracer_source(:,:,:,:,:)  ! tracer source                         ("var"/s)

  INTEGER (KIND=iintegers), TARGET, ALLOCATABLE::           &
! Note: 2 dim. switch: ltracer(processes(i.e. advection, turbulence, ...),ntracer)
    ltracer(:,:),&            ! switch for tracer
    ltracer_rad(:)            ! switch for radiation and aerosol type    

  INTEGER (KIND=iintegers) ::          &
    ntracer,&                 ! number of tracers
    ntracer_adv,&             ! number of tracers with advection
    ntracer_conv,&            ! number of tracers with convection
    ntracer_hdiff,&           ! number of tracers with horizontal diffusion
    ntracer_rad,&             ! number of tracers with radiation
    ntracer_dim,&             ! number of tracers for dimensioning of fields (min = 1)
    nsrc1,nsrc2,&             ! indices for permutation of tracer source data  
    hincsource,&              ! time increment for tracer source data
    nincsource,&              ! time step increment for tracer source data
    nlastsource               ! time step for last tracer source data

! needed for CO2 coupling
  REAL (KIND=ireals) ::      &
    molmass_da  = 28.954_ireals, &    ! molar mass of dry air                        (g/mol)
    molmass_w   = 18.015_ireals, &    ! molar mass of water                          (g/mol)
    molmass_co2 = 44.010_ireals       ! molar mass of CO2                            (g/mol)

!CPS#ifdef COUP_OAS_COS
  REAL  (KIND=ireals), TARGET, ALLOCATABLE ::           &
    co2fl_s(:,:), &                   ! CO2 tendency due to land-atmosphere exchange ( 1/s )
    psn_tens(:,:), &                  ! CO2 tendency due to photosynthesis           ( 1/s )
    plres_tens(:,:)                   ! CO2 tendency due to plant respiration        ( 1/s )
!CPS#endif

! variables needed for radiation 
  REAL (KIND=ireals) ::      &
    seanew(8,8,4),&           ! extinction coefficients for aerosol types
    saanew(8,8,4),&           ! absorption coefficients for aerosol types
    ganew(8,8,4),&            ! asymmetry factors for aerosol types
    feux(8)                   ! reference relative humidities
  REAL  (KIND=ireals), TARGET, ALLOCATABLE ::           &
    dz_aer(:,:), &            ! dz (for calculation of radiation)
    rh_aer(:,:), &            ! relative humidity (for calculation of radiation) 
    n_aer(:,:,:)              ! two-dimensional tracer number density (1/m^3)
  LOGICAL  ::    laer_cosmo   ! switch for cosmo aerosol on/off                  
 
END MODULE data_tracer
