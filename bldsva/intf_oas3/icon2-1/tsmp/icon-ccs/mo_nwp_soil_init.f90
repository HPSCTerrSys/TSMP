!>
!! Initialization routine for multi-layer soil model
!!
!! <Describe the concepts of the procedures and algorithms used in the module.>
!! <Details of procedures are documented below with their definitions.>
!! <Include any applicable external references inline as module::procedure,>
!! <external_procedure(), or by using @see.>
!! <Don't forget references to literature.>
!!
!! @author Juergen Helmert, DWD
!! @author Ekaterina Machulskaya, DWD
!! @author Guenther Zaengl, DWD
!! @author Daniel Reinert, DWD
!!
!!
!! @par Revision History
!! Moved to separate module from mo_soil_ml.f90 by Daniel Reinert, DWED (2014-04-01)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_nwp_soil_init


  USE mo_kind,               ONLY: ireals=>wp, &
    &                              iintegers=>i4
  USE mo_run_config,         ONLY: msg_level
  USE mo_physical_constants, ONLY: t0_melt => tmelt,& ! absolute zero for temperature
    &                              rho_w => rhoh2o, & ! density of liquid water (kg/m^3)
    &                              lh_f  => alf   , & ! latent heat of fusion
    &                              g     => grav      ! acceleration due to gravity

  USE mo_lnd_nwp_config,     ONLY: lmulti_snow, lana_rho_snow, &
    &                              lmelt, lmelt_var,           &
    &                              max_toplaydepth, l2lay_rho_snow
  USE mo_phyparam_soil,      ONLY: cdw0, cdw1, ckw0, ckw1, cik2, &
    &                              cporv, cpwp, cadp, cfcap,     & 
    &                              crock, crhoc, cala0, cala1,   &
    &                              csandf, cclayf, csnow_tmin,   &
   &                               crhosmax_tmin, crhosmax_ml


  IMPLICIT NONE

  PRIVATE



  PUBLIC  :: terra_multlay_init
  PUBLIC  :: get_wsnow

CONTAINS

  !>
  !! Performs initialization of soil model TERRA
  !!
  !! Performs initialization of soil model TERRA. Depending on the initialization mode 
  !! chosen, a coldstart- , warmstart-  or warmstart with snow increments (IAU)-initialization
  !! is performed
  !!
  !!
  !! @par Revision History
  !! Initial Revision by Juergen Helmert, DWD (2011)
  !!
  SUBROUTINE terra_multlay_init (                & 
                  init_mode,    & ! 1 = coldstart, 2 = warmstart, 3 = warmstart with snow increments (IAU)
                  ie,           & ! array dimensions
                  istartpar,    & ! start index for computations in the parallel program
                  iendpar,      & ! end index for computations in the parallel program
! 4. variables for the time discretization and related variables
!                  ke,           & ! nsubs0=1 for single tile, nsubs0=2 for multi-tile
!                  nsubs0ubs1,& ! nsubs1=1 for single tile, nsubs1=#tiles+1 for multi-tile
                  ke_soil, ke_snow , &
                  czmls            , & ! processing soil level structure 
                  soiltyp_subs     , & ! type of the soil (keys 0-9)                     --
                  rootdp           , & ! depth of the roots                            ( m  )
                  plcov            , & ! fraction of surface covered by plants         ( -  )
                  t_snow_now       , & ! temperature of the snow-surface               (  K  )
                  t_snow_mult_now  , & ! temperature of the snow-surface               (  K  )
                  t_s_now          , & ! temperature of the ground surface             (  K  )
                  t_s_new          , & ! temperature of the ground surface             (  K  )
                  w_snow_now       , & ! water content of snow                         (m H2O)
                  h_snow           , & ! snow depth                                   (m H2O)
                  rho_snow_now     , & ! snow density                                  (kg/m**3)
                  rho_snow_mult_now, & ! snow density                                  (kg/m**3)
                  t_so_now         , & ! soil temperature (main level)                 (  K  )
                  t_so_new         , & ! soil temperature (main level)                 (  K  )
                  w_so_now         , & ! total water conent (ice + liquid water)       (m H20)
                  w_so_new         , & ! total water conent (ice + liquid water)       (m H20)
                  w_so_ice_now     , & ! ice content                                   (m H20)
                  w_so_ice_new     , & ! ice content                                   (m H20)
                  wliq_snow_now    , & ! liquid water content in the snow              (m H2O)
                  wtot_snow_now    , & ! total (liquid + solid) water content of snow  (m H2O)
                  dzh_snow_now       & ! layer thickness between half levels in snow   (  m  )
                                  )

                                      

!-------------------------------------------------------------------------------
! Declarations for ICON (USE statements for COSMO)
!-------------------------------------------------------------------------------

  IMPLICIT NONE


  INTEGER                 , INTENT(IN)  ::  &
                  init_mode    ! 1 = coldstart, 2 = warmstart, 3 = warmstart with snow increments (IAU)
  INTEGER (KIND=iintegers), INTENT(IN)  ::  &
                  ie,           & ! array dimensions
                  istartpar,    & ! start index for computations in the parallel program
                  iendpar,      & ! end index for computations in the parallel program
!                  ke,           &
!                 nsubs0ubs1 , &
                  ke_soil, ke_snow      
  REAL    (KIND = ireals), DIMENSION(ke_soil+1), INTENT(IN) :: &
                  czmls           ! processing soil level structure 

  INTEGER (KIND=iintegers),DIMENSION(ie), INTENT(IN) :: & 
                  soiltyp_subs      ! type of the soil (keys 0-9)                     --
  REAL    (KIND = ireals), DIMENSION(ie), INTENT(IN) :: & 
                  rootdp           ! depth of the roots                            ( m  )
  REAL    (KIND = ireals), DIMENSION(ie), INTENT(IN) :: & 
                  plcov           ! plant coverage                            ( -  )
  REAL    (KIND = ireals), DIMENSION(ie), INTENT(INOUT) :: &
                  t_snow_now              ! temperature of the snow-surface (K)
  REAL    (KIND = ireals), DIMENSION(ie,0:ke_snow), INTENT(INOUT) :: &
                  t_snow_mult_now      ! temperature of the snow-surface               (  K  )
  REAL    (KIND = ireals), DIMENSION(ie), INTENT(INOUT) :: &
                  t_s_now              ! temperature of the ground surface             (  K  )
  REAL    (KIND = ireals), DIMENSION(ie), INTENT(OUT) :: &
                  t_s_new              ! temperature of the ground surface             (  K  )
  REAL    (KIND = ireals), DIMENSION(ie), INTENT(INOUT) :: &
                  w_snow_now       , & ! water content of snow                         (m H2O)
                  rho_snow_now         ! snow density                                  (kg/m**3)
  REAL    (KIND = ireals), DIMENSION(ie), INTENT(INOUT) :: &
                  h_snow               ! snow depth                                   (m H2O)
  REAL    (KIND = ireals), DIMENSION(ie,ke_snow), INTENT(INOUT) :: &
                  rho_snow_mult_now    ! snow density                                  (kg/m**3)
  REAL    (KIND = ireals), DIMENSION(ie,0:ke_soil+1), INTENT(INOUT) :: &
                  t_so_now             ! soil temperature (main level)                 (  K  )
  REAL    (KIND = ireals), DIMENSION(ie,0:ke_soil+1), INTENT(OUT) :: &
                  t_so_new             ! soil temperature (main level)                 (  K  )
  REAL    (KIND = ireals), DIMENSION(ie,ke_soil+1), INTENT(INOUT) :: &
                  w_so_now         , & ! total water conent (ice + liquid water)       (m H20)
                  w_so_ice_now         ! ice content                                   (m H20)
  REAL    (KIND = ireals), DIMENSION(ie,ke_soil+1), INTENT(OUT) :: &
                  w_so_new         , & ! total water conent (ice + liquid water)       (m H20)
                  w_so_ice_new         ! ice content                                   (m H20)
  REAL    (KIND = ireals), DIMENSION(ie,ke_snow), INTENT(INOUT) :: &
                  wliq_snow_now    , & ! liquid water content in the snow              (m H2O)
                  wtot_snow_now        ! total (liquid + solid) water content of snow  (m H2O)
  REAL    (KIND = ireals), DIMENSION(ie,ke_snow), INTENT(INOUT) :: &
                  dzh_snow_now         ! layer thickness between half levels in snow   (  m  )

!--------------------------------------------------------------------------------
! TERRA Declarations

! New declaration for ICON

!------------------------------------------------------------------------------
! Subroutine arguments: None
! --------------------

!
! Local parameters:
! ----------------

  REAL    (KIND=ireals   ), PARAMETER ::  &
    zepsi  = 1.0E-6_ireals ! security constant

  REAL    (KIND=ireals   ), PARAMETER ::  &
    T_ref_ice = 0.1_ireals , & !degC Soil ice parameterization
    T_star_ice= 0.01_ireals, & !degC according to K. Schaefer and Jafarov, E.,2016, doi:10.5194/bg-13-1991-2016
    b_clay= -0.3_ireals , & 
    b_silt= -0.5_ireals , &
    b_sand= -0.9_ireals , &
    b_org = -1.0_ireals 

! Local scalars:
! -------------

  INTEGER (KIND=iintegers) ::  &
    kso            , & ! loop index for soil moisture layers           
    ksn,k          , & ! loop index for snow layers
    i,ic           , & ! loop index in x-direction              
    icount_soil    , & ! "true" soil
    icount_rockice , & ! rock and ice points
    mstyp          , & ! soil type index
    istarts        , & ! start index for x-direction      
    iends              ! end   index for x-direction

!  REAL    (KIND=ireals   ) ::  &
!    zdel_t_so          ! auxiliary variable

  REAL    (KIND=ireals   ) ::  &
    zpsi0=.01_ireals    ! air entry potential at water saturation (m)

  REAL    (KIND=ireals   ) ::  &
    zroota    (ie) ! root density profile parameter (1/m)


  REAL    (KIND=ireals   ) ::  &
    zaa, zh_snow, zzz, zd, zd1, zd2, zd3, zd4        ! utility variables

  REAL    (KIND = ireals) :: &
!
    zmls     (ke_soil+1)  , & ! depth of soil main level
    zzhls    (ke_soil+1)  , & ! depth of the half level soil layers in m
    zdzhs    (ke_soil+1)  , & ! layer thickness between half levels
    zdzms    (ke_soil+1)  , & ! distance between main levels

!
! External parameters
!
    zbwt     (ie)      , & ! root depth (with artificial minimum value)
    zsandf   (ie)      , & ! mean fraction of sand (weight percent)
    zclayf   (ie)      , & ! mean fraction of clay (weight percent)
    zsiltf   (ie)      , & ! mean fraction of silt (weight percent)
    zb_por   (ie)      , & ! pore size distribution index
    zpsis    (ie)      , & ! air entry potential (m)
    zw_m_org           , &  ! maximum of  liquid water content   organic
    zw_m_soil          , &  ! maximum of  liquid water content   mineral soil
    zw_m_up            , &  ! maximum of  liquid water content   at temp -3 degC
    zw_m_low           , &  ! maximum of  liquid water content   at temp -40 degC
    zw_m_org_low       , &  ! maximum of  liquid water content   organic at temp -40 degC
    zw_m_soil_low      , &  ! maximum of  liquid water content   mineral soil at temp -40 degC
    t_zw_up            , &  ! temp -3 degC
    t_zw_low           , &  ! temp -40 degC
    zw_m     (ie)          ! maximum of liquid water content  (m)

  INTEGER  (KIND=iintegers ) ::  &
    m_styp   (ie)      , & ! soil type
    zicount1 (ie)      , &
    zicount2 (ie)      , &    
    soil_list(ie)      , & ! list of "true" soil points
                           ! i.e. no ice no rock       
    rockice_list(ie)       ! list of rock and ice points
  REAL    (KIND=ireals   ) ::  & ! Soil and plant parameters
!
    zadp     (ie)      , & ! air dryness point
    zporv    (ie)      , & ! pore volume (fraction of volume)
    zedb     (ie)          ! utility variable

  REAL    (KIND=ireals   ) ::  &
    zw_snow_old(ie)      ,& !
    zrho_snow_old(ie)    ,&
    h_snow_fg(ie)        ,&
    h_snow_incr(ie)      ,&
    zhh_snow (ie,ke_snow),& ! depth of the half level snow layers
    zhm_snow (ie,ke_snow),& ! depth of snow main levels
    sum_weight(ie)      , &
    t_new  (ie,ke_snow) , &
    rho_new(ie,ke_snow) , &
    wl_new (ie,ke_snow) , &
    z_old  (ie,ke_snow) , & 
    dz_old (ie,ke_snow) , &
    weight, wso_ice_tolerance, wso_ice_equil, rwso_ice_tolerance

  LOGICAL :: l_redist(ie)

!- End of header
!==============================================================================


!------------------------------------------------------------------------------
! Begin Subroutine terra_multlay              
!------------------------------------------------------------------------------
!==============================================================================
!  Computation of the diagnostic part I of the soil parameterization scheme
!  In this part, evaporation from the surface and transpiration is calculated.
!  A multi-layer soil water distribution is calculated by simple bulk
!  parameterisation or a Penman-Monteith-type version of evaporation
!  and transpiration which can be used alternatively.
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Section I.1: Initializations
!------------------------------------------------------------------------------

! Horizontal domain for computation
  istarts = istartpar
  iends   = iendpar

! grids for temperature and water content

  zzhls(1) = 2._ireals*czmls(1)   !depth of first half level
  zdzhs(1) = zzhls(1)      !layer thickness betw. half levels of uppermost layer
  zmls(1)  = czmls(1)      !depth of 1st main level
  zdzms(1) = czmls(1)      !layer thickness between soil surface and main level
                           ! of uppermost layer
 
  DO kso = 2,ke_soil+1
    zzhls(kso)  = zzhls(kso-1) + 2._ireals*(czmls(kso) -zzhls(kso-1))
    zdzhs(kso) = zzhls(kso) - zzhls(kso-1) ! layer thickness betw. half levels
    zmls(kso)  = czmls(kso)                ! depth of main levels
    zdzms(kso) = zmls(kso) - zmls(kso-1)   ! layer thickness betw. main levels
  ENDDO




! Prepare basic surface properties (for land-points only)
  icount_soil=0
  icount_rockice=0
  DO i = istarts, iends
    mstyp       = soiltyp_subs(i)        ! soil type
    m_styp(i) = mstyp                     ! array for soil type
    IF (mstyp >= 3) THEN
      icount_soil=icount_soil+1
      soil_list(icount_soil)=i
    ELSE
      icount_rockice=icount_rockice+1
      rockice_list(icount_rockice)=i
    END IF
    zporv(i)  = cporv(mstyp)              ! pore volume
    zadp (i)  = cadp (mstyp)              ! air dryness point
    zbwt(i)   = MAX(0.001_ireals,rootdp(i))! Artificial minimum value
                                           ! for root depth
    zroota(i) = 3._ireals/zbwt(i)       ! root density profile parameter (1/m)
                                        ! zroota=0. creates the original TERRA_LM
                                        ! version with constant root density
                                        ! for root depth

    ! Arrays for soil water freezing/melting
    zsandf(i)   = csandf(mstyp)
    zclayf(i)   = cclayf(mstyp)
    zsiltf(i)   = 100._ireals -csandf(mstyp)-cclayf(mstyp) ! Residuum of sand and clay
    zpsis(i)    = -zpsi0 * 10._ireals**(1.88_ireals-0.013_ireals*zsandf(i))
    zb_por(i)   = 2.91_ireals + .159_ireals*zclayf(i)
    zedb(i)     = 1._ireals/zb_por(i)

! print*,i,zporv(i),zclayf(i),zpsis(i),zb_por(i)

  ENDDO
 
  ! Further parameters for soil water freezing/melting
  t_zw_up  = 270.15_ireals ! temp -3 degC
  t_zw_low = 233.15_ireals ! temp -40 degC
  zd = LOG((T_ref_ice-(t_zw_low-t0_melt))/T_star_ice)
  zd1 = EXP(b_sand*zd)
  zd2 = EXP(b_clay*zd)
  zd3 = EXP(b_silt*zd)
  zd4 = EXP(b_org*zd)

  ! Tolerances for the allowed deviation of w_so_ice from its equilibrium value
  wso_ice_tolerance  = 1.05_ireals ! 5% deviation allowed
  rwso_ice_tolerance = 1._ireals/wso_ice_tolerance


! For ntstep=0 : Some preparations
! ================================

  w_so_new(:,:) = w_so_now(:,:)

! SBr, CHa: DO NOT set minimal values to w_so if the threshold is passed
!! Provide for a soil moisture 1 % above air dryness point, reset soil
!! moisture to zero in case of ice and rock
!  DO kso   = 1,ke_soil+1
!    DO ic = 1, icount_soil
!      i=soil_list(ic)
!      w_so_now (i,kso) = MAX(w_so_now(i,kso),                     &
!                         1.01_ireals*zadp(i)*zdzhs(kso) )
!      w_so_new (i,kso) = MAX(w_so_new(i,kso),                     &
!                         1.01_ireals*zadp(i)*zdzhs(kso) )
!    END DO
!    DO ic = 1, icount_rockice
!      i=rockice_list(ic)
!      w_so_now(i,kso) = 0.0_ireals
!      w_so_new(i,kso) = 0.0_ireals
!    END DO
!  END DO


! adjust temperature profile in lower soil layers, if temperature of first soil
! layer was reduced due to the detection of snow (e.g. in the analysis)
! loop over grid points
  DO i = istarts, iends
    IF (w_snow_now(i) <= zepsi) THEN
      ! spurious snow is removed
      t_snow_now(i)=t_so_now(i,0)
!jh   t_snow_new(i)=t_so_now(i,0)
      w_snow_now(i)= 0.0_ireals
      IF(lmulti_snow) THEN
        t_snow_mult_now(i,0) = t_so_now(i,0)
        DO ksn = 1, ke_snow
          t_snow_mult_now(i,ksn) = t_so_now(i,0)
          wliq_snow_now(i,ksn) = 0.0_ireals
          wtot_snow_now(i,ksn) = 0.0_ireals
          rho_snow_mult_now(i,ksn) = 0.0_ireals
          dzh_snow_now(i,ksn) = 0.0_ireals
        END DO
      END IF
! GZ: this soil temperature adjustment tends to cause more damage than benefit in practice
! because if there are traces of snow over warm ground, soil heat fluxes can become so large
! that the scheme becomes numerically unstable if the soil top temperature is artificially 
! adjusted to freezing temperature
    ELSE
! !           adjust soil temperatures in presence of snow
!             zdel_t_so = MAX (0._ireals,t_so_now(i,0)-t0_melt)
!             t_so_now(i,0)=MIN( t0_melt,t_so_now(i,0) )
!             t_so_new(i,0)=MIN( t0_melt,t_so_now(i,0) )
!             DO kso=1,ke_soil
!               IF ( t_so_now(i,kso) > t0_melt) THEN
!                  t_so_now(i,kso) = MAX( t0_melt,                  &
!                                    t_so_now(i,kso) -zdel_t_so    &
!                                    *(zmls(ke_soil+1)-zmls(kso))   &
!                                    /(zmls(ke_soil+1)-zmls( 1 )) )
!                  t_so_new(i,kso) = MAX( t0_melt,                  &
!                                    t_so_now(i,kso) -zdel_t_so    &
!                                    *(zmls(ke_soil+1)-zmls(kso))   &
!                                    /(zmls(ke_soil+1)-zmls( 1 )) )
!               ENDIF
!               IF ( t_so_now(i,kso) <= t0_melt) EXIT
!             ENDDO
      t_snow_now(i) = MIN( t0_melt, t_snow_now(i))
      IF(lmulti_snow) THEN
        DO ksn = 0, ke_snow
          t_snow_mult_now(i,ksn) = MIN( t0_melt, t_snow_mult_now(i,ksn))
        END DO
      END IF 
    ENDIF

    t_s_now(i) = t_so_now(i,0)
    t_s_new(i) = t_so_now(i,0)

!   Set level 1 to level 0 for t_so for every landpoint
    t_so_now(i,1) = t_so_now(i,0)
    t_so_new(i,1) = t_so_now(i,0)
    t_so_new(i,0) = t_so_now(i,0)

  END DO


  ! Note that a few parts need to be shifted to the .NOT. is_coldstart branch.
  ! I.e. those which are dealing with the contributions from the analysis.
  ! Currently, the multi layer soil model starts from the FG.
  IF ( init_mode == 1 ) THEN

    ! Initialization of 
    ! - snow density (if necessary)
    ! - multi layer snow variables (if necessary)
    !   * wtot_snow, dzh_snow, wliq_snow
    ! - soil ice water content w_so_ice
    ! - snow height h_snow

!   Initialization of snow density, if necessary
!   --------------------------------------------
    IF (.NOT. lana_rho_snow) THEN
       IF(lmulti_snow) THEN
         ksn = 1
         DO i = istarts, iends
           rho_snow_mult_now(i,ksn) = 250.0_ireals    ! average initial density
           t_snow_mult_now  (i,ksn) = t_snow_now(i)
           wliq_snow_now    (i,ksn) = 0.0_ireals
           ! Limit top layer to max_toplaydepth
           dzh_snow_now     (i,ksn) = MIN(max_toplaydepth, w_snow_now(i)/REAL(ke_snow,ireals) &
                  &                            *rho_w/rho_snow_mult_now(i,ksn))
           wtot_snow_now    (i,ksn) = dzh_snow_now(i,ksn)*rho_snow_mult_now(i,ksn)/rho_w
           t_snow_mult_now  (i,0  ) = t_snow_now(i)
         ENDDO
         k = MIN(2, ke_snow-1)
         DO ksn = 2, ke_snow
           DO i = istarts, iends
             rho_snow_mult_now(i,ksn) = 250.0_ireals    ! average initial density
             t_snow_mult_now  (i,ksn) = t_snow_now(i)
             wliq_snow_now    (i,ksn) = 0.0_ireals
             ! distribute the remaining snow equally among the layers
             IF (k == ksn) THEN
               dzh_snow_now   (i,ksn) = MIN(8._ireals*max_toplaydepth, (w_snow_now(i)-wtot_snow_now(i,1)) &
                    &                            /REAL(ke_snow-1,ireals)*rho_w/rho_snow_mult_now(i,ksn) )
             ELSE
               dzh_snow_now   (i,ksn) = (w_snow_now(i)-wtot_snow_now(i,k))/REAL(ke_snow-k,ireals) &
                    &                            *rho_w/rho_snow_mult_now(i,ksn)
             ENDIF
             wtot_snow_now    (i,ksn) = dzh_snow_now(i,ksn)*rho_snow_mult_now(i,ksn)/rho_w
           ENDDO
         ENDDO
       ELSE
         rho_snow_now(istarts:iends) = 250.0_ireals    ! average initial density
         IF (l2lay_rho_snow) rho_snow_mult_now(istarts:iends,:) = 250.0_ireals
       ENDIF
    ELSE
      IF(lmulti_snow) THEN
!The sum of wtot_snow, i.e. the "old" total water equivalent depth
        zw_snow_old  (:) = 0.0_ireals
        zrho_snow_old(:) = 0.0_ireals
        zicount1(:)      = 0
        zicount2(:)      = 0
        sum_weight(:)    = 0.0_ireals 
        DO ksn = 1, ke_snow
          DO i = istarts, iends
            zw_snow_old(i) = zw_snow_old(i) + wtot_snow_now(i,ksn)
            zrho_snow_old(i) = zrho_snow_old(i) + dzh_snow_now(i,ksn)
            IF (wtot_snow_now(i,ksn) > zepsi) zicount1(i) = zicount1(i)+1
            IF (dzh_snow_now(i,ksn) > zepsi)  zicount2(i) = zicount2(i)+1
          END DO
        END DO
!The "old" snow density
        DO i = istarts, iends
          zrho_snow_old(i) = zw_snow_old(i) / MAX(zrho_snow_old(i),1.0E-09_ireals) *rho_w
        END DO 
        k = MIN(2, ke_snow-1)
        DO ksn = 1, ke_snow
          DO i = istarts, iends
            IF(zicount1(i).EQ.ke_snow .AND. zicount2(i).EQ.ke_snow) THEN
              rho_snow_mult_now(i,ksn) = wtot_snow_now(i,ksn)/dzh_snow_now(i,ksn)*rho_w
              wtot_snow_now(i,ksn) = wtot_snow_now(i,ksn)*w_snow_now(i)/zw_snow_old(i)
              dzh_snow_now (i,ksn) = dzh_snow_now (i,ksn)*w_snow_now(i)/zw_snow_old(i)
              wliq_snow_now(i,ksn) = wliq_snow_now(i,ksn)*w_snow_now(i)/zw_snow_old(i)
            ELSE
              IF(rho_snow_now(i) .EQ. 0._ireals) rho_snow_now(i) = 250._ireals
              IF (ksn == 1) THEN ! Limit top layer to max_toplaydepth
                IF (t_snow_now(i) > 268.15_ireals) THEN
                  rho_snow_mult_now(i,ksn) = MIN(400._ireals,rho_snow_now(i))
                ELSE ! at very low temperatures, the near-surface snow is usually not subject to strong densification
                  rho_snow_mult_now(i,ksn) = MIN(rho_snow_now(i),400._ireals+6._ireals* &
                    (MAX(228.15_ireals,t_snow_now(i))-268.15_ireals)                    )
                ENDIF
                dzh_snow_now(i,ksn) = MIN(max_toplaydepth, w_snow_now(i)/REAL(ke_snow,ireals) &
                  &                  *rho_w/rho_snow_mult_now(i,ksn))
              ELSE ! distribute the remaining snow equally among the layers
                rho_snow_mult_now(i,ksn) = rho_snow_now(i)
                IF (ksn == k) THEN
                  dzh_snow_now(i,ksn) = MIN(8._ireals*max_toplaydepth, (w_snow_now(i)-wtot_snow_now(i,1)) &
                      &                 /REAL(ke_snow-1,ireals)*rho_w/rho_snow_mult_now(i,ksn) )
                ELSE
                  dzh_snow_now(i,ksn) = (w_snow_now(i)-sum_weight(i))/REAL(ke_snow-k,ireals) &
                      &                 *rho_w/rho_snow_mult_now(i,ksn)
                ENDIF
              ENDIF
              wtot_snow_now    (i,ksn) = dzh_snow_now(i,ksn)*rho_snow_mult_now(i,ksn)/rho_w
              sum_weight(i) = sum_weight(i) + wtot_snow_now(i,ksn)
              wliq_snow_now    (i,ksn) = 0.0_ireals
            END IF
          END DO
        END DO
        DO i = istarts, iends
          t_snow_mult_now  (i,0  ) = t_snow_now(i)
        END DO
      ELSE
        DO i = istarts, iends
          IF(rho_snow_now(i) .EQ. 0._ireals) rho_snow_now(i) = 250._ireals
          IF (l2lay_rho_snow) THEN
            IF(rho_snow_mult_now(i,1) == 0._ireals) rho_snow_mult_now(i,1) = rho_snow_now(i)
            IF(rho_snow_mult_now(i,2) == 0._ireals) rho_snow_mult_now(i,2) = rho_snow_now(i)
          ENDIF
        END DO
      END IF
    ENDIF


!   Initialization of soil ice content
!   ----------------------------------
!   Generally the combination lmelt=.true., lmelt_var=.true. has to be used.
!   lmelt_var=.false. (soil water freezing at 273.15K) should be used only
!   for special investigations because at model start a partial frozen layer
!   is only considered if the soil ice water content is read in from
!   a pre run.
!   ----------------------------------
    IF (lmelt .AND. .NOT. lmelt_var) THEN
      DO kso   = 1,ke_soil+1
        DO i = istarts, iends
          w_so_ice_now(i,kso) = 0.0_ireals
          w_so_ice_new(i,kso) = 0.0_ireals

          IF (t_so_now(i,kso) < t0_melt) THEN
            w_so_ice_now(i,kso) = w_so_now(i,kso)
            w_so_ice_new(i,kso) = w_so_now(i,kso)
          END IF ! t_so(kso) < t0_melt
        END DO
      END DO
    END IF           ! lmelt .AND. .NOT. lmelt_var
    IF(lmelt .AND. lmelt_var) THEN
      DO kso   = 1,ke_soil+1
        DO i = istarts, iends
          IF (t_so_now(i,kso) < (t0_melt-zepsi)) THEN 

            zaa    = g*zpsis(i)/lh_f

            ! J. Helmert: Soil ice parameterization according to K. Schaefer and Jafarov, E.,2016,
            !                                                    doi:10.5194/bg-13-1991-2016

            ! Liq. water content at -3 degC
            zw_m_up =  zporv(i)*zdzhs(kso)*EXP(-zedb(i)*LOG((t_zw_up - t0_melt)/(t_zw_up*zaa)) )

            ! Determine liq. water content at -40 degC
            zw_m_soil = 0.01_ireals*(zsandf(i)*zd1 + zclayf(i)*zd2 + zsiltf(i)*zd3)

            ! Scale soil ice content with organic soil horizon.
            ! should decrease the root zone liquid water content of frozen soil for low temperatures significantly!
            IF(zmls(kso) < rootdp(i)) THEN
              zzz = plcov(i)*(rootdp(i)-zmls(kso))/rootdp(i)
              zw_m_low = zporv(i)*zdzhs(kso)*(zzz*zd4 + (1._ireals-zzz)*zw_m_soil)
            ELSE
              zw_m_low = zporv(i)*zdzhs(kso)*zw_m_soil
            END IF
   
            IF (t_so_now(i,kso) < t_zw_low) THEN
              zw_m(i) = zw_m_low
            ELSE IF (t_so_now(i,kso) < t_zw_up) THEN ! Logarithmic Interpolation between -3 degC and -40 degC 
              zw_m(i) = zw_m_low*EXP((t_so_now(i,kso) - t_zw_low)*(LOG(zw_m_up) - LOG(zw_m_low))/(t_zw_up-t_zw_low))
            ELSE
              zw_m(i) = zporv(i)*zdzhs(kso)*EXP(-zedb(i)*LOG((t_so_now(i,kso)-t0_melt)/(t_so_now(i,kso)*zaa)) )
            END IF

            w_so_ice_now(i,kso) = MAX (0.0_ireals,w_so_now(i,kso) - zw_m(i))
            w_so_ice_new(i,kso) = MAX (0.0_ireals,w_so_now(i,kso) - zw_m(i))
          ELSE ! ensure that w_so_ice is zero
            w_so_ice_now(i,kso) = 0.0_ireals
            w_so_ice_new(i,kso) = 0.0_ireals
          END IF
        END DO
      END DO
    END IF           ! lmelt .AND. lmelt_var


!   Initialization of the local array of the grid in snow
!   -----------------------------------------------------
    IF(lmulti_snow) THEN
      sum_weight(:)    = 0.0_ireals 
      DO i = istarts, iends
        h_snow(i) = 0.0_ireals
        DO ksn = 1,ke_snow
          h_snow(i) = h_snow(i) + dzh_snow_now(i,ksn)
        END DO
      END DO
      DO ksn = 1,ke_snow
        DO i = istarts, iends
          IF(h_snow(i) > 0.0_ireals) THEN
            t_snow_mult_now(i,ksn) = t_snow_now(i) + &
              (t_s_now(i)-t_snow_now(i))*(sum_weight(i)+dzh_snow_now(i,ksn)*0.5)/h_snow(i)
            sum_weight(i) = sum_weight(i) + dzh_snow_now(i,ksn)
          END IF
        END DO
      END DO
    ELSE
      DO i = istarts, iends
        h_snow(i) = w_snow_now(i)/rho_snow_now(i)*rho_w
      END DO
    END IF

  ELSE  ! init_mode > 1, i.e. warmstart (i.e. assimilation cycle and forecast)
        ! in this case, w_snow and rho_snow are read from the first guess, and h_snow is updated by the snow analysis;
        ! in the IAU mode, the snow depth increment (h_snow_incr) is directly provided by the snow analysis


    IF ( init_mode == 2 ) THEN  ! snow depth increment needs to be computed
      DO i = istarts, iends
        h_snow_fg(i)   = w_snow_now(i)/rho_snow_now(i)*rho_w
        ! Self-consistency fix needed because of possible GRIB truncation errors
        IF (lmulti_snow .AND. h_snow_fg(i) > zepsi) THEN
          dzh_snow_now(i,2:ke_snow) = MAX(dzh_snow_now(i,1),dzh_snow_now(i,2:ke_snow))
          zh_snow       = SUM(dzh_snow_now(i,1:ke_snow))
          w_snow_now(i) = w_snow_now(i) * zh_snow/h_snow_fg(i)
          h_snow_fg(i)  = zh_snow
        ELSE IF (lmulti_snow) THEN
          dzh_snow_now(i,:) = 0._ireals
        ENDIF
        h_snow_incr(i) = h_snow(i) - h_snow_fg(i)
      ENDDO
    ELSE IF ( init_mode == 3 ) THEN  ! snow depth increment is provided from snow analysis
      IF (lmulti_snow) THEN
        DO i = istarts, iends
          ! Self-consistency fix needed because of possible GRIB truncation errors
          IF (dzh_snow_now(i,1) > zepsi) THEN
            dzh_snow_now(i,2:ke_snow) = MAX(dzh_snow_now(i,1),dzh_snow_now(i,2:ke_snow))
            h_snow_fg(i) = SUM(dzh_snow_now(i,1:ke_snow))
          ELSE
            dzh_snow_now(i,:) = 0._ireals
            h_snow_fg(i)      = 0._ireals
          ENDIF
          h_snow_incr(i) = h_snow(i) - h_snow_fg(i)
        ENDDO
      ELSE
        DO i = istarts, iends
          h_snow_fg(i) = h_snow(i)
          h_snow_incr(i) = 0._ireals
        ENDDO
      ENDIF
    ENDIF

    ! Provisional fix to get rid of isolated extremely cold points in parallel routine (2014-10-16)
    ! This is applied only to the lowest prognostic soil layer, i.e. ke_soil
    DO i = istarts, iends
      IF (MIN(t_so_now(i,ke_soil+1),t_so_now(i,ke_soil-1)) - t_so_now(i,ke_soil) > 5._ireals) THEN
        t_so_now(i,ke_soil) = MIN(t_so_now(i,ke_soil+1),t_so_now(i,ke_soil-1)) - 5._ireals
      ENDIF
      IF (t_so_now(i,ke_soil) - MAX(t_so_now(i,ke_soil+1),t_so_now(i,ke_soil-1)) > 5._ireals) THEN
        t_so_now(i,ke_soil) = MAX(t_so_now(i,ke_soil+1),t_so_now(i,ke_soil-1)) + 5._ireals
      ENDIF
    ENDDO

    ! Ensure that w_so_ice stays within 5% of its equilibrium value
    ! In addition, it must not exceed w_so
    IF(lmelt .AND. lmelt_var) THEN 
      DO kso   = 1,ke_soil+1
        DO i = istarts, iends
          IF (t_so_now(i,kso) < (t0_melt-zepsi)) THEN 

            zaa    = g*zpsis(i)/lh_f

            ! J. Helmert: Soil ice parameterization according to K. Schaefer and Jafarov, E.,2016,
            !                                                    doi:10.5194/bg-13-1991-2016

            ! Liq. water content at -3 degC
            zw_m_up =  zporv(i)*zdzhs(kso)*EXP(-zedb(i)*LOG((t_zw_up - t0_melt)/(t_zw_up*zaa)) )

            ! Determine liq. water content at -40 degC
            zw_m_soil = 0.01_ireals*(zsandf(i)*zd1 + zclayf(i)*zd2 + zsiltf(i)*zd3)

            ! Scale soil ice content with organic soil horizon.
            ! should decrease the root zone liquid water content of frozen soil for low temperatures significantly!
            IF(zmls(kso) < rootdp(i)) THEN
              zzz = plcov(i)*(rootdp(i)-zmls(kso))/rootdp(i)
              zw_m_low = zporv(i)*zdzhs(kso)*(zzz*zd4 + (1._ireals-zzz)*zw_m_soil)
            ELSE
              zw_m_low = zporv(i)*zdzhs(kso)*zw_m_soil
            END IF
   
            IF (t_so_now(i,kso) < t_zw_low) THEN
              zw_m(i) = zw_m_low
            ELSE IF (t_so_now(i,kso) < t_zw_up) THEN ! Logarithmic Interpolation between -3 degC and -40 degC 
              zw_m(i) = zw_m_low*EXP((t_so_now(i,kso) - t_zw_low)*(LOG(zw_m_up) - LOG(zw_m_low))/(t_zw_up-t_zw_low))
            ELSE
              zw_m(i) = zporv(i)*zdzhs(kso)*EXP(-zedb(i)*LOG((t_so_now(i,kso)-t0_melt)/(t_so_now(i,kso)*zaa)) )
            END IF

            wso_ice_equil = MAX (0.0_ireals,w_so_now(i,kso) - zw_m(i))
            w_so_ice_now(i,kso) = MIN(w_so_ice_now(i,kso),  wso_ice_tolerance*wso_ice_equil, w_so_now(i,kso))
            w_so_ice_now(i,kso) = MAX(w_so_ice_now(i,kso), rwso_ice_tolerance*wso_ice_equil)
          ELSE ! ensure that w_so_ice is zero
            w_so_ice_now(i,kso) = 0.0_ireals
          END IF
        END DO
      END DO
    END IF

    ! For the single-layer snow case, the prognostic variables now need to be updated. This requires appropriate
    ! assumptions about the new snow density

    IF(lmulti_snow) THEN ! update of prognostic variables of multi-layer snow model
      l_redist(:) = .FALSE.
      DO i = istarts, iends
        IF (h_snow(i) <= zepsi) THEN  ! snow has disappeared or was already absent
          rho_snow_now(i) = 250._ireals
          w_snow_now(i)   = 0._ireals

          dzh_snow_now (i,1:ke_snow) = 0._ireals
          wtot_snow_now(i,1:ke_snow) = 0._ireals
          rho_snow_mult_now(i,1:ke_snow) = 0._ireals
          wliq_snow_now(i,1:ke_snow) = 0._ireals
          t_snow_mult_now(i,1:ke_snow) = t_snow_now(i)
        ELSE IF (h_snow_fg(i) <= zepsi .AND. h_snow(i) > zepsi) THEN 

          ! snow analysis has either reestablished a snow cover that had been erroneously melted away by the model
          ! or generated a new snow cover that was missed by the model due to a lack of precipitation (or wrong phase of precip)
          ! ** snow density should be distinguished between new snow and reestablished old snow **
          rho_snow_now(i) = 250._ireals ! needs to be improved
          w_snow_now(i)   = h_snow(i)*rho_snow_now(i)/rho_w
          t_snow_now(i)   = MIN(t_snow_now(i),t0_melt)

          rho_snow_mult_now(i,1:ke_snow) = rho_snow_now(i)
          dzh_snow_now (i,1:ke_snow) = h_snow(i)/REAL(ke_snow,ireals)
          wtot_snow_now(i,1:ke_snow) = w_snow_now(i)/REAL(ke_snow,ireals)
          wliq_snow_now(i,1:ke_snow) = 0._ireals
          t_snow_mult_now(i,1:ke_snow) = t_snow_now(i)
          IF (h_snow(i) > REAL(ke_snow,ireals)*max_toplaydepth) l_redist(i) = .TRUE.
        ELSE IF (h_snow_fg(i) > zepsi .AND. h_snow_incr(i) < 0._ireals) THEN

          ! a snow cover is present but reduced by the snow analysis increment
          ! in this case, we retain the snow densities but adjust the water equivalents
          l_redist(i) = .TRUE.
          dz_old(i,1) = MAX(dzh_snow_now (i,1),zepsi/REAL(ke_snow,ireals))
          dzh_snow_now (i,1) = dzh_snow_now (i,1) + h_snow_incr(i)
          wtot_snow_now(i,1) = wtot_snow_now(i,1)*dzh_snow_now (i,1)/dz_old(i,1)
        ELSE
          ! a snow cover is present and increased by the snow analysis
          ! ** snow density should be distinguished between new snow and reestablished old snow **
          ! IF (w_snow_incr(i) > 0._ireals) THEN ! assume temperature-dependent fresh snow density - maybe air temp should be used here?!
          !   rho_new(i,1) = MAX(50._ireals,125._ireals-10._ireals*(t0_melt-t_snow_mult_now(i,1)))
          ! ELSE
          !   rho_new(i,1) = rho_snow_mult_now(i,1)
          ! ENDIF
          l_redist(i) = .TRUE.
          rho_new(i,1) = rho_snow_mult_now(i,1)
          wtot_snow_now(i,1) = wtot_snow_now(i,1) + h_snow_incr(i)*rho_new(i,1)/rho_w
          dzh_snow_now (i,1) = dzh_snow_now (i,1) + h_snow_incr(i)
          rho_snow_mult_now(i,1) = rho_w*wtot_snow_now(i,1)/dzh_snow_now (i,1)
        ENDIF
      ENDDO

      ! Redistribute snow cover among layers where necessary
      sum_weight(:) = 0._ireals
      DO ksn = ke_snow,1,-1
        DO i = istarts, iends
          IF (l_redist(i)) THEN
            dz_old(i,ksn) = dzh_snow_now(i,ksn)
            z_old(i,ksn) = -sum_weight(i) - dzh_snow_now(i,ksn)/2._ireals
            sum_weight(i) = sum_weight(i) + dzh_snow_now(i,ksn)
          ENDIF
        END DO
      END DO

      k = MIN(2, ke_snow-1)
      DO ksn = 1,ke_snow
        DO i = istarts, iends
          IF (l_redist(i)) THEN

            IF (ksn == 1) THEN ! Limit top layer to max_toplaydepth
              zhh_snow(i,ksn) = -MAX( h_snow(i)-max_toplaydepth, h_snow(i)/ke_snow*(ke_snow-ksn) )
              zhm_snow(i,ksn) = (-h_snow(i) + zhh_snow(i,ksn))/2._ireals   
              dzh_snow_now(i,ksn) = zhh_snow(i,ksn) + h_snow(i)    !layer thickness betw. half levels of uppermost snow layer
            ELSE IF (ksn == 2 .AND. ke_snow > 2) THEN ! Limit second layer to 8*max_toplaydepth
              zhh_snow(i,ksn) = MIN( 8._ireals*max_toplaydepth+zhh_snow(i,1), zhh_snow(i,1)/(ke_snow-1)*(ke_snow-ksn) )
            ELSE ! distribute the remaining snow equally among the layers
              zhh_snow(i,ksn) = zhh_snow(i,k)/(ke_snow-k)*(ke_snow-ksn)
            ENDIF

            t_new  (i,ksn) = 0.0_ireals
            rho_new(i,ksn) = 0.0_ireals
            wl_new (i,ksn) = 0.0_ireals

            IF(dz_old(i,ksn) > 0._ireals .AND. rho_snow_mult_now(i,ksn) > 0._ireals) THEN
              wliq_snow_now(i,ksn) = wliq_snow_now(i,ksn)/dz_old(i,ksn)
            END IF

          ENDIF
        END DO
      END DO

      DO ksn = 2, ke_snow
        DO i = istarts, iends
          IF(l_redist(i)) THEN
            zhm_snow(i,ksn) = (zhh_snow(i,ksn) + zhh_snow(i,ksn-1))/2._ireals
            dzh_snow_now(i,ksn) = zhh_snow(i,ksn) - zhh_snow(i,ksn-1) ! layer thickness betw. half levels
          END IF
        END DO
      END DO

      ! mass-weighted recomputation of temperatures and densities
      DO ksn = ke_snow,1,-1 
        DO k = ke_snow,1,-1
          DO i = istarts, iends
            IF(l_redist(i)) THEN
         
              weight = MAX(MIN(z_old(i,k)+dz_old(i,k)/2._ireals,zhm_snow(i,ksn) + dzh_snow_now(i,ksn)/2._ireals)-   &
                       MAX(z_old(i,k)-dz_old(i,k)/2._ireals, zhm_snow(i,ksn)-dzh_snow_now(i,ksn)/2._ireals),0._ireals) &
                       &/dzh_snow_now(i,ksn)
      
              t_new  (i,ksn) = t_new  (i,ksn) + t_snow_mult_now  (i,k)*weight
              rho_new(i,ksn) = rho_new(i,ksn) + rho_snow_mult_now(i,k)*weight
              wl_new (i,ksn) = wl_new (i,ksn) + wliq_snow_now    (i,k)*weight
            END IF
          END DO
        END DO
      END DO

      DO ksn = ke_snow,1,-1
        DO i = istarts, iends
          IF(l_redist(i)) THEN
            t_snow_mult_now  (i,ksn) = t_new  (i,ksn)
            rho_snow_mult_now(i,ksn) = rho_new(i,ksn)
            wtot_snow_now    (i,ksn) = rho_new(i,ksn)*dzh_snow_now(i,ksn)/rho_w
            wliq_snow_now    (i,ksn) = wl_new (i,ksn)*dzh_snow_now(i,ksn)
          END IF
        END DO   
      END DO


      ! Re-diagnose integrated/averaged snow-fields, in order to 
      ! have a consistent state.
      !
      h_snow      (istarts:iends) = 0.0_ireals
      w_snow_now  (istarts:iends) = 0.0_ireals
      rho_snow_now(istarts:iends) = 0.0_ireals
      t_snow_now  (istarts:iends) = 0.0_ireals
      ! Re-diagnose h_snow, w_snow, rho_snow, t_snow
      DO ksn = ke_snow,1,-1
        DO i = istarts, iends
          h_snow(i)    = h_snow(i)     + dzh_snow_now(i,ksn)
          w_snow_now(i)= w_snow_now(i) + wtot_snow_now(i,ksn)
        END DO
      ENDDO

      DO ksn = ke_snow,1,-1
        DO i = istarts, iends
          IF (h_snow(i) > zepsi) THEN
            rho_snow_now(i) = rho_snow_now(i) + rho_snow_mult_now(i,ksn) * dzh_snow_now(i,ksn)/h_snow(i)
          ELSE
            rho_snow_now(i) = 250._ireals
          ENDIF
        END DO
      ENDDO

      DO i = istarts, iends
        t_snow_now(i) = t_snow_mult_now(i,1)
      ENDDO


    ELSE  ! single snow layer

      ! diagnose w_snow and rho_snow from analyzed h_snow
      ! h_snow is taken from analysis, while up to this point 
      ! rho_snow_now and w_snow_now contain the first guess.
      CALL get_wsnow(h_snow,         &  ! in
        &            rho_snow_now,   &  ! inout
        &            t_snow_now,     &  ! in
        &            istarts, iends, &  ! in
        &            soiltyp_subs,   &  ! in
        &            w_snow_now      )  ! out

      DO i = istarts, iends
        IF (w_snow_now(i) <= zepsi) THEN
          ! spurious snow is removed
          t_snow_now(i)= t_so_now(i,0)
          w_snow_now(i)= 0.0_ireals
        ELSE
          t_snow_now(i) = MIN( t0_melt, t_snow_now(i))
        ENDIF
      ENDDO

    ENDIF ! multi snow

  ENDIF  ! init_mode



! End of timestep 0 preparations
! ==============================

  END SUBROUTINE terra_multlay_init



  !>
  !! Transform snow depth to snow water equivalent
  !!
  !! Description:
  !!  Snow depth observations are reported in m and this is also the analysis
  !!  variable. The prognostic model variables describing the evolution of
  !!  snow are snow water eqivalent and snow density.
  !!  This subroutine calculates snow water equivalent from snow depth.
  !!  For non-permanent ice points, hsnow is limited to 4m. w_snow is limited 
  !!  accordingly, using rho_snow from the first guess.
  !!
  !! Method:
  !!  Simple linear relation between snow water and snow depth using snow density
  !!
  !! @par Revision History
  !! Adapted to ICON by Daniel Reinert, DWD (2014-04-01)
  !! Original implementation by Michael Buchhold (see below)
  !!
  SUBROUTINE get_wsnow(h_snow, rho_snow, t_snow, istart, iend, soiltyp, w_snow)

    REAL(ireals), INTENT(INOUT) :: h_snow(:)    ! snow depth            [m]
    REAL(ireals), INTENT(INOUT) :: rho_snow(:)  ! snow density          [kg/m**3]
    REAL(ireals), INTENT(IN)    :: t_snow(:)    ! snow temperature      [K]
    REAL(ireals), INTENT(INOUT) :: w_snow(:)    ! snow water equivalent [m H2O]
    INTEGER,  INTENT(IN)    :: istart       ! start index
    INTEGER,  INTENT(IN)    :: iend         ! end index
    INTEGER,  INTENT(IN)    :: soiltyp(:)   ! type of soil (1-9)

    ! local
    INTEGER :: jc         ! loop index

    REAL(ireals), PARAMETER :: rho_snw_default=250._ireals
  !--------------------------------------------------------------------

    ! calculate water equivalent from snow depth
    !
    DO jc = istart, iend
      IF (rho_snow(jc) < 1._ireals) rho_snow(jc)=rho_snw_default      ! average initial density
      IF (h_snow(jc) > 0._ireals) THEN
        ! multiply analysed snow depth [m] by (first-guess) density 
        ! to get water equivalent in [m H2O]
        IF (soiltyp(jc) /= 1) THEN       ! 1=ice
          ! limit snow depth to 4m for non-glacier points
          h_snow(jc) = MIN(h_snow(jc), 4._ireals)
          w_snow(jc) = h_snow(jc) * rho_snow(jc)/rho_w
        ELSE
          ! limit snow depth to at least 1 m and set snow density to the equilibrium value for the current snow temperature
          h_snow(jc)   = MAX(h_snow(jc), 1._ireals)
          rho_snow(jc) = crhosmax_tmin + MAX(0._ireals,(t_snow(jc)-csnow_tmin)/(t0_melt-csnow_tmin)) &
                         * (crhosmax_ml-crhosmax_tmin)
          w_snow(jc)   = h_snow(jc) * rho_snow(jc)/rho_w
        ENDIF  ! soiltyp
      ENDIF
    ENDDO


  END SUBROUTINE get_wsnow




!!$!
!!$!+ Transform snow depth to snow water equivalent and vice versa
!!$!
!!$      SUBROUTINE snw_tran (snw,rho,ie,je,imode)
!!$!
!!$! Description:
!!$!  Snow depth observations are reported in m and this is also the analysis
!!$!  variable. The prognostic model variables describing the evolution of
!!$!  snow are  snow water eqivalent and snow density.
!!$!  This subroutine calculates snow depth from water eqivalent and vice versa.
!!$!
!!$! Method:
!!$!  Simple linear relation between snow water and snow depth using snow density
!!$!
!!$! Input  files: none
!!$! Output files: none
!!$!
!!$! Current Code Owner: DWD, Michael Buchhold                          
!!$!  phone:  +49  69  8062 2726
!!$!  fax:    +49  69  8062 3721
!!$!  email:  michael.buchhold@dwd.de
!!$!
!!$! History:
!!$! Version    Date       Name
!!$! ---------- ---------- ----
!!$! 1.4        2003/07/08 Michael Buchhold
!!$!  Initial release
!!$! 1.13       2005/04/18 Michael Buchhold
!!$!  Replacement of common blocks by data modules
!!$! 1.17       2005/12/01 Michael Buchhold
!!$!  Snow density is used for transformation between snow water
!!$!  and snow depth
!!$! V1_26        2010/04/27 Michael Gertz
!!$!  Adaptions to SVN
!!$!
!!$! Code Description:
!!$! Language: Fortran 90.
!!$! Software Standards: "European Standards for Writing and
!!$! Documenting Exchangeable Fortran 90 Code".
!!$!=======================================================================
!!$!
!!$! Declarations:
!!$!
!!$! Modules used:
!!$      USE data_all, ONLY:
!!$     & rhosnw
!!$
!!$      IMPLICIT NONE
!!$!=======================================================================
!!$!
!!$! Subroutine / Function arguments
!!$! Scalar arguments with intent(in):
!!$      INTEGER, INTENT(IN) ::  ie, je   ! Dimensions of snws
!!$      INTEGER, INTENT(IN) ::  imode  
!!$!     imode= 1: calculate water equivalent from snow depth
!!$!     imode=-1: calculate snow depth from water equivalent
!!$
!!$! Array arguments with intent(inout):
!!$      REAL, INTENT(INOUT) ::   snw (ie, je), rho(ie,je)
!!$!
!!$! Local Type Definitions:
!!$
!!$!
!!$! Local Parameters:
!!$
!!$!
!!$! Local Scalars:
!!$      INTEGER i,j
!!$      REAL zcrit, zwcrit
!!$!
!!$! Local Arrays:
!!$
!!$!=======================================================================
!!$! Include statements: none
!!$!
!!$!- End of header
!!$!------------------------------------------------------------------------
!!$      if (imode == 1) then
!!$!     calculate water equivalent from snow depth
!!$        do j = 1,je
!!$        do i = 1,ie
!!$          if (snw(i,j) > 0.) then
!!$             if (rho(i,j) < 1.) rho(i,j)=rhosnw
!!$!           scale analysed snow depth [mm] to get in m and multiply it by density 
!!$            snw(i,j) = snw(i,j)/1000.*rho(i,j)
!!$          endif
!!$        enddo
!!$        enddo
!!$
!!$      else if (imode ==-1) then
!!$!     calculate  snow depth from water eqivalent
!!$        do j = 1,je
!!$        do i = 1,ie
!!$          if (snw(i,j) >= 0.) then
!!$            if (rho(i,j) < 1.) rho(i,j)=rhosnw
!!$            snw(i,j) = snw(i,j)/rho(i,j)
!!$          else
!!$            rho(i,j) = 0.
!!$          endif
!!$!         scale snow depth to have it in mm 
!!$          snw(i,j) = snw(i,j)*1000.
!!$        enddo
!!$        enddo
!!$      endif
!!$
!!$      end


END MODULE mo_nwp_soil_init
