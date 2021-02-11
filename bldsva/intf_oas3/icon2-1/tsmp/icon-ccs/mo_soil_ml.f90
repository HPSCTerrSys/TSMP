!>
!! Soil Vegetation Atmosphere Transfer (SVAT) scheme TERRA
!! Source Module  "mo_soil_ml.f90"
!! "Nihil in TERRA sine causa fit." (Cicero)
!!------------------------------------------------------------------------------
!!-----------------------------------------------------------------------------
!!
!! @par Description:
!!   The module "mo_soil_ml_v413.f90" performs calculations related to the
!!   parameterization of soil processes. It contains the soubroutine
!!   terra_multlay which is the combination of the former parts
!!   terra1_multlay.incf and terra2_multlay.incf of the LM.
!!   All parametric scalar and array data for this soil model routine are
!!   defined in the data module data_soil.f90.
!!
!!   All global variables of the model that are used by the soil model routine
!!   terra_multlay.incf is imported by USE statements below.
!!
!!   The parameterization package has been provided by B. Ritter in a
!!   Plug-compatible Fortran77-Version, which is based on the EM/DM soil model
!!   by E. Heise. Some technical modifications have been done for the
!!   F90 and the parallel Version:
!!   Internal communication by common-blocks is replaced by module parameters,
!!   scalars and arrays defined in module data_soil. The plug compatible
!!   I/O lists of the subroutines have been replaced by the Module interface
!!   defined by Use lists below.
!!
!! @author E. Heise, R. Schrodin, B. Ritter
!! @author E. Machulskaya, F. Ament, J. Helmert
!!
!! @par reference   This is an adaption of subroutine hydci_pp in file src_gscp.f90
!!  of the lm_f90 library (COSMO code). Equation numbers refer to
!!  Doms, Foerstner, Heise, Herzog, Raschendorfer, Schrodin, Reinhardt, Vogel
!!    (September 2005): "A Description of the Nonhydrostatic Regional Model LM",
!!
!!
!! @par Revision History
!! implemented into ICON by K. Froehlich, E. Machulskaya, and J. Helmert (2010-11-XX)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
!!==============================================================================
!
!
!
! Current Code Owner: DWD, Reinhold Schrodin
!  phone:  +49  69  8062 2709
!  fax:    +49  69  8062 3721
!  email:  reinhold.schrodin@dwd.de
!
! History:
! Version    Date       Name
! ---------- ---------- ----
! 2.17       2002/05/08 Reinhold Schrodin
!  Initial release
! 2.18       2002/07/16 Reinhold Schrodin
!  Corrections and further developments
! 3.6        2003/12/11 Reinhold Schrodin
!  Reformulation of bare soil evaporation, transpiration, snow density,
!  snow height, snow melting, subsoil runoff calculation
! 3.7        2004/02/18 Reinhold Schrodin / Ulrich Schaettler
!  Adapted use of nold to 2 time level scheme and some corrections
! 3.13       2004/12/03 Reinhold Schrodin, Bodo Ritter, Erdmann Heise
!  Snow albedo:   Snow age indicator calculation (freshsnow)
!  Transpiration: 2m temperature and 10m wind considered instead temperature
!                 and wind at lowest atmospheric main level
!  Runoff:        Reformulation of soil water runoff (gravitational part)
!  Soil water:    Soil water calculations restricted to ke_soil_hy,
!                 determined by soil level nearest 2.5 m, combined with
!                 a reformulation of the tridiagonal equation system
!  Snow melting:  Correction in case of pore volume overshooting of
!                 uppermost soil layer by infiltration of snow water
!  Precautions to avoid
!   - soil water content below air dryness point
!   - soil surface temperature problems at lateral LM boundaries.
!   - excessive turbulent fluxes at the soil surface due to mismatch
!     of atmospheric and soil variables in the initial state
!   - evaporation if soil water content is below air dryness point
! 3.16       2005/07/22 Reinhold Schrodin
!   - Modification of evaporation treatment at plant wilting point
!   - Soil water transport: Soil water diffusion term driven by the gradient
!                        of the total soil water concentration (water + ice)
!   - Combining terra1_multlay and terra2_multlay in terra_multlay
! 3.17       2005/12/12 Reinhold Schrodin
!   LME adaptations to GME soil model modifications (B. Ritter):
!   - Constant snow density (250 kg/m**3) replaced by prognostic snow density
!     formulation (min = 50 kg/m**3, max = 400 kg/m**3)
!   - Adjustment of soil temperatures in presence of snow
!   - Extended formulation of transfer coefficient limitation
! 3.18       2006/03/03 Ulrich Schaettler / Erdmann Heise
!  Adaptations to do some initializations also for restart runs
!  Field h_snow is now time dependent, because of use in the FLake model
!  Increase infiltration and reduce surface runoff (bugfix) (Erdmann Heise)
! 3.19       2006/04/25 Erdmann Heise
!  Remove spurious snow and set consistent t_snow and w_snow at the beginning
! 3.21       2006/12/04 Ulrich Schaettler
!  crsmin, rat_lam put to data_soil; eliminated variables that are not used
!  Use dt for ntstep=0 instead of dt/2
!  Additional use of graupel, if present (Thorsten Reinhardt)
!  Modifications to avoid soil water content below air dryness point
!  Gravitational runoff calculation modified in case of frozen soil layer
!                                                      (Reinhold Schrodin)
! V3_23        2007/03/30 Matthias Raschendorfer
!  Introduction of tfv to consider different laminar resistance for heat
!  and water vapour; thus 'rat_lam' is not used here any longer.
! V4_4         2008/07/16 Ulrich Schaettler
!  Splitting of a loop in Section I.4.3b (m_styp is not defined for sea points
!  and must not occur together with llandmask in the same IF-clause)
! V4_7         2008/12/12 Ulrich Schaettler
!  There were still some loops left with llandmask and m_styp in one line
! V4_9         2009/07/16 Ulrich Schaettler, Christian Bollmann
!  Inserted Call to collapse loops
! V4_10        2009/09/11 Christian Bollmann
!  Added compiler directive to use option _on_adb for NEC
! V4_11        2009/11/30 Ekaterina Machulskaya, Juergen Helmert, Lucio Torrisi
!  Implementation of multi-layer snow model (EM)
!  Use of an external parameter field for stomata resistance (JH)
!  Implementation of ground water as lower boundary of soil column and
!    soil moisture dependent heat conductivity of the soil (JH)
!  Save additional fluxes and stomata resistance to global memory for output (LT)
! V4_12        2010/05/11 Ulrich Schaettler, Ekaterina Machulskaya
!  Renamed t0 to t0_melt because of conflicting names
!  Renamed prs_min to rsmin2d because of conflicting names
!  Update of the new snow model (EM)
! V4_13        2010/05/11 Michael Gertz
!  Adaptions to SVN
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!==============================================================================

MODULE mo_soil_ml
!
! Declarations:
!
! Modules used:

#ifdef __COSMO__

USE data_parameters, ONLY :   &
    ireals,    & ! KIND-type parameter for real variables
    iintegers    ! KIND-type parameter for standard integer variables

!------------------------------------------------------------------------------

USE data_modelconfig, ONLY :   &
! 2. horizontal and vertical sizes of the fields and related variables
! --------------------------------------------------------------------
    ie,           & ! number of grid points in zonal direction
    ke,           & ! number of grid points in vertical direction
    ke_soil,      & ! number of layers in multi-layer soil model
    ke_snow,      & ! number of layers in multi-layer soil model
    czmls,        & ! depth of the main level soil layers in m
                    ! (see organize_physics.f90)
! 3. start- and end-indices for the computations in the horizontal layers
! -----------------------------------------------------------------------
!    These variables give the start- and the end-indices of the
!    forecast for the prognostic variables in a horizontal layer.
!    Note, that the indices for the wind-speeds u and v differ from
!    the other ones because of the use of the staggered Arakawa-C-grid.
!
    istartpar,    & ! start index for computations in the parallel program
    iendpar,      & ! end index for computations in the parallel program
! 4. variables for the time discretization and related variables
! --------------------------------------------------------------
    dt              ! long time-step
!   dt2             ! dt*2.

! end of data_modelconfig

!------------------------------------------------------------------------------

USE data_constants  , ONLY :   &

! 1. physical constants and related variables
! -------------------------------------------

    t0_melt,      & ! absolute zero for temperature
    r_d,          & ! gas constant for dry air
    rdv,          & ! r_d / r_v
    o_m_rdv,      & ! 1 - r_d/r_v
    rvd_m_o,      & ! r_v/r_d -1
    cp_d,         & ! specific heat of dry air at constant pressure
    rdocp,        & ! r_d / cp_d
    lh_v,         & ! latent heat of vapourization
    lh_f,         & ! latent heat of fusion
    lh_s,         & ! latent heat of sublimation
    g,            & ! acceleration due to gravity
    sigma,        & ! Boltzmann-constant

! 2. constants for parametrizations
! ---------------------------------
    b1,           & ! variables for computing the saturation vapour pressure
    b2w,          & ! over water (w) and ice (i)
    b2i,          & !               -- " --
    b3,           & !               -- " --
    b4w,          & !               -- " --
    b4i,          & !               -- " --
    rho_w           ! density of liquid water

! end of data_constants

!------------------------------------------------------------------------------

USE data_fields     , ONLY :   &

! 1. constant fields for the reference atmosphere                     (unit)
! -----------------------------------------------
    p0         ,    & ! base state pressure                           (Pa)

! 2. external parameter fields                                        (unit)
! ----------------------------
    soiltyp_subs    ,    & ! type of the soil (keys 0-9)                     --
    plcov      ,    & ! fraction of plant cover                         --
    rootdp     ,    & ! depth of the roots                            ( m  )
    sai        ,    & ! surface area index                              --
    tai        ,    & ! transpiration area index                        --
    eai        ,    & ! earth area (evaporative surface area) index     --
    llandmask  ,    & ! landpoint mask                                  --
    rsmin2d    ,    & ! minimum stomata resistance                    ( s/m )


! 3. prognostic variables                                             (unit)
! -----------------------
    u          ,    & ! zonal wind speed                              ( m/s )
    v          ,    & ! meridional wind speed                         ( m/s )
    t          ,    & ! temperature                                   (  k  )
    qv         ,    & ! specific water vapor content                  (kg/kg)
    pp         ,    & ! deviation from the reference pressure         ( pa  )

! 5. fields for surface values and soil model variables               (unit )
! -----------------------------------------------------
    ps        ,     & ! surface pressure                              ( pa  )
    t_snow    ,     & ! temperature of the snow-surface               (  K  )
    t_snow_mult,    & ! temperature of the snow-surface               (  K  )
    t_s       ,     & ! temperature of the ground surface             (  K  )
    t_g       ,     & ! weighted surface temperature                  (  K  )
    qv_s      ,     & ! specific humidity at the surface              (kg/kg)
    w_snow    ,     & ! water content of snow                         (m H2O)
    rho_snow  ,     & ! snow density                                  (kg/m**3)
    rho_snow_mult,  & ! snow density                                  (kg/m**3)
    h_snow    ,     & ! snow height                                   (  m
    w_i       ,     & ! water content of interception water           (m H2O)
    t_so      ,     & ! soil temperature (main level)                 (  K  )
    w_so      ,     & ! total water conent (ice + liquid water)       (m H20)
    swi       ,     & ! soil wetness index                            (  1  )
    w_so_ice  ,     & ! ice content                                   (m H20)
!   t_2m      ,     & ! temperature in 2m                             (  K  )
!   u_10m     ,     & ! zonal wind in 10m                             ( m/s )
!   v_10m     ,     & ! meridional wind in 10m                        ( m/s )
    uv_low     ,     & ! wind speed at lowest model level              ( m/s )
    freshsnow ,     & ! indicator for age of snow in top of snow layer(  -  )
    wliq_snow ,     & ! liquid water content in the snow              (m H2O)
    wtot_snow ,     & ! total (liquid + solid) water content of snow  (m H2O)
    dzh_snow          ! layer thickness between half levels in snow   (  m  )

USE data_fields     , ONLY :   &

! 6. fields that are computed in the parametrization and dynamics     (unit )
! ---------------------------------------------------------------
!   fields of convective and grid-scale precipitation
    prr_con     ,   & ! precipitation rate of rain, convective        (kg/m2*s)
    prs_con     ,   & ! precipitation rate of snow, convective        (kg/m2*s)
    prr_gsp     ,   & ! precipitation rate of rain, grid-scale        (kg/m2*s)
    prs_gsp     ,   & ! precipitation rate of snow, grid-scale        (kg/m2*s)
    prg_gsp     ,   & ! precipitation rate of graupel, grid-scale     (kg/m2*s)

!   fields of the turbulence scheme
    tch         ,   & ! turbulent transfer coefficient for heat       ( -- )
    tcm         ,   & ! turbulent transfer coefficient for momentum   ( -- )
    tfv         ,   & ! laminar reduction factor for evaporation      ( -- )

!   fields of the radiation
    sobs        ,   & ! solar radiation at the ground                 ( W/m2)
    thbs        ,   & ! thermal radiation at the ground               ( W/m2)
    pabs        ,   & ! photosynthetic active radiation               ( W/m2)

! 7. fields for model output and diagnostics                          (unit )
! ---------------------------------------------------------------
    runoff_s    ,   & ! surface water runoff; sum over forecast      (kg/m2)
    runoff_g    ,   & ! soil water runoff; sum over forecast         (kg/m2)
    rstom       ,   & ! stomata resistance                           ( s/m )
    lhfl_bs     ,   & ! average latent heat flux from bare soil evap.( W/m2)
    lhfl_pl           ! average latent heat flux from plants         ( W/m2)

! end of data_fields

!------------------------------------------------------------------------------

USE data_runcontrol , ONLY :   &

! 1. start and end of the forecast
! --------------------------------
    nstart,       & ! first time step of the forecast
    ntstep,       & ! actual time step
                    ! indices for permutation of three time levels
    nold  ,       & ! corresponds to ntstep - 1
    nnow  ,       & ! corresponds to ntstep
    nnew  ,       & ! corresponds to ntstep + 1
    lmulti_snow,  & ! run the multi-layer snow model

! 3. controlling the physics
! --------------------------
    itype_trvg,   & ! type of vegetation transpiration parameterization
    itype_evsl,   & ! type of parameterization of bare soil evaporation
    itype_root,   & ! type of root density distribution
    itype_heatcond,&! type of soil heat conductivity
    itype_hydbound,&! type of hydraulic lower boundary condition
    lstomata,     & ! map of minimum stomata resistance

! 5. additional control variables
! --------------------------
     l2tls           ! forecast with 2-TL integration scheme

! end of data_runcontrol

!------------------------------------------------------------------------------

USE data_parallel,      ONLY:  &
    my_cart_id        ! rank of this subdomain in the global communicator

!------------------------------------------------------------------------------

USE data_soil       ! All variables from data module "data_soil are used by
                    ! this module. These variables start with letter "c"

!------------------------------------------------------------------------------

! External routines used by this module

USE environment,     ONLY : collapse
USE meteo_utilities, ONLY : tgcom

#endif

#ifdef __ICON__
!
USE mo_kind,               ONLY: ireals=>wp,    &
                                 iintegers=>i4
USE mo_math_constants    , ONLY: pi
!
USE mo_physical_constants, ONLY: t0_melt => tmelt,& ! absolute zero for temperature
                                 r_d   => rd    , & ! gas constant for dry air
                                 rvd_m_o=>vtmpc1, & ! r_v/r_d - 1
                                 o_m_rdv        , & ! 1 - r_d/r_v
                                 rdv            , & ! r_d / r_v
                                 lh_v  => alv   , & ! latent heat of vapourization
                                 lh_s  => als   , & ! latent heat of sublimation
                                 lh_f  => alf   , & ! latent heat of fusion
                                 cp_d  => cpd   , & ! specific heat of dry air at constant press
                                 g     => grav  , & ! acceleration due to gravity
                                 sigma => stbo  , & ! Boltzmann-constant
                                 rho_w => rhoh2o, & ! density of liquid water (kg/m^3)
                                 rdocp => rd_o_cpd  ! r_d / cp_d
!
USE mo_convect_tables,     ONLY: b1    => c1es  , & !! constants for computing the sat. vapour
                                 b2w   => c3les , & !! pressure over water (l) and ice (i)
                                 b2i   => c3ies , & !!               -- " --
                                 b4w   => c4les , & !!               -- " --
                                 b4i   => c4ies , & !!               -- " --
                                 b234w => c5les     !!               -- " --

!
USE mo_phyparam_soil
!
USE mo_lnd_nwp_config,     ONLY: lmulti_snow, l2lay_rho_snow,     &
  &                              itype_trvg, itype_evsl,          &
  &                              itype_root, itype_heatcond,      &
  &                              itype_hydbound, lstomata,        &
  &                              max_toplaydepth, itype_interception, &
  &                              cwimax_ml
!
!
USE mo_exception,          ONLY: message, message_text
USE mo_run_config,         ONLY: msg_level
USE mo_impl_constants,     ONLY: iedmf
#endif



!------------------------------------------------------------------------------
! Declarations
!------------------------------------------------------------------------------

IMPLICIT NONE

PRIVATE


!------------------------------------------------------------------------------
! Public subroutines
!------------------------------------------------------------------------------

PUBLIC :: terra_multlay

!------------------------------------------------------------------------------
! Public variables
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
! Parameters and variables which are global in this module
!------------------------------------------------------------------------------

#ifdef COSMO
CHARACTER(132) :: message_text = ''

!AD:Setting this identifier as global in this module
!   was the only way for now. In case iedmf is redefined
!   in mo_impl_constants, it should be chagned in here too
INTEGER, PARAMETER :: iedmf   =  3
#endif


CONTAINS

#ifdef COSMO
SUBROUTINE message (name, text, all_print)
IMPLICIT NONE

CHARACTER (*) :: name, text
LOGICAL, INTENT(in), OPTIONAL :: all_print

LOGICAL :: lprint

IF (PRESENT(all_print)) THEN
  lprint = all_print
ELSE
  lprint = .FALSE.
ENDIF

IF (lprint .OR. my_cart_id==0) &
  WRITE (*,*) TRIM(name), TRIM(text)

END SUBROUTINE message

#endif



!==============================================================================
!+ Computation of the first part of the soil parameterization scheme
!------------------------------------------------------------------------------

!option! -pvctl _on_adb


  SUBROUTINE terra_multlay (         &
                  icant            , & ! canopy type
                  ie               , & ! array dimensions
                  istartpar        , & ! start index for computations in the parallel program
                  iendpar          , & ! end index for computations in the parallel program
                  ldiag_tg         , & ! if true: diagnose t_g and snow-cover fraction with tgcom
                  ke_soil, ke_snow , &
                  ke_soil_hy       , & ! number of active soil moisture layers
                  czmls            , & ! processing soil level structure
                  inwp_turb        , & ! turbulence scheme number
                  nclass_gscp      , & ! number of hydrometeor classes of grid scale microphysics
                  dt               , & ! time step
!
                  soiltyp_subs     , & ! type of the soil (keys 0-9)                     --
                  plcov            , & ! fraction of plant cover                         --
                  rootdp           , & ! depth of the roots                            ( m  )
                  sai              , & ! surface area index                              --
                  tai              , & ! transpiration area index                        --
                  eai              , & ! earth area (evaporative surface area) index     --
!                 llandmask        , & ! landpoint mask                                  --
                  rsmin2d          , & ! minimum stomata resistance                    ( s/m )
!
                  u                , & ! zonal wind speed                              ( m/s )
                  v                , & ! meridional wind speed                         ( m/s )
                  t                , & ! temperature                                   (  k  )
                  qv               , & ! specific water vapor content                  (kg/kg)
                  p0               , & !!!! base state pressure                        ( Pa  )
!                 pp               , & ! deviation from the reference pressure         ( Pa  )
                  ps               , & ! surface pressure                              ( Pa  )
!
                  t_snow_now       , & ! temperature of the snow-surface               (  K  )
                  t_snow_new       , & ! temperature of the snow-surface               (  K  )
!
                  t_snow_mult_now  , & ! temperature of the snow-surface               (  K  )
                  t_snow_mult_new  , & ! temperature of the snow-surface               (  K  )
!
                  t_s_now          , & ! temperature of the ground surface             (  K  )
                  t_s_new          , & ! temperature of the ground surface             (  K  )
!
                  t_g              , & ! weighted surface temperature                  (  K  )
                  qv_s             , & ! specific humidity at the surface              (kg/kg)
                  w_snow_now       , & ! water content of snow                         (m H2O)
                  w_snow_new       , & ! water content of snow                         (m H2O)
!
                  rho_snow_now     , & ! snow density                                  (kg/m**3)
                  rho_snow_new     , & ! snow density                                  (kg/m**3)
!
                  rho_snow_mult_now, & ! snow density                                  (kg/m**3)
                  rho_snow_mult_new, & ! snow density                                  (kg/m**3)
!
                  h_snow           , & ! snow depth                                   (  m  )
                  h_snow_gp        , & ! grid-point averaged snow depth               (  m  )
                  meltrate         , & ! snow melting rate                             (kg/(m**2*s))
!
                  w_i_now          , & ! water content of interception water           (m H2O)
                  w_i_new          , & ! water content of interception water           (m H2O)
!
                  w_p_now          , & ! water content of pond interception water     (m H2O)
                  w_p_new          , & ! water content of pond interception water     (m H2O)
!
                  w_s_now          , & ! water content of interception snow           (m H2O)
                  w_s_new          , & ! water content of interception snow           (m H2O)
!
                  t_so_now         , & ! soil temperature (main level)                 (  K  )
                  t_so_new         , & ! soil temperature (main level)                 (  K  )
!
                  w_so_now         , & ! total water conent (ice + liquid water)       (m H20)
                  w_so_new         , & ! total water conent (ice + liquid water)       (m H20)
!
                  w_so_ice_now     , & ! ice content                                   (m H20)
                  w_so_ice_new     , & ! ice content                                   (m H20)
!
!                 t_2m             , & ! temperature in 2m                             (  K  )
                  u_10m            , & ! zonal wind in 10m                             ( m/s )
                  v_10m            , & ! meridional wind in 10m                        ( m/s )
                  freshsnow        , & ! indicator for age of snow in top of snow layer(  -  )
                  zf_snow          , & ! snow-cover fraction                           (  -  )
!
                  wliq_snow_now    , & ! liquid water content in the snow              (m H2O)
                  wliq_snow_new    , & ! liquid water content in the snow              (m H2O)
!
                  wtot_snow_now    , & ! total (liquid + solid) water content of snow  (m H2O)
                  wtot_snow_new    , & ! total (liquid + solid) water content of snow  (m H2O)
!
                  dzh_snow_now     , & ! layer thickness between half levels in snow   (  m  )
                  dzh_snow_new     , & ! layer thickness between half levels in snow   (  m  )
!
                  prr_con          , & ! precipitation rate of rain, convective        (kg/m2*s)
                  prs_con          , & ! precipitation rate of snow, convective        (kg/m2*s)
                  conv_frac        , & ! convective area fraction as assumed in convection scheme
                  prr_gsp          , & ! precipitation rate of rain, grid-scale        (kg/m2*s)
                  prs_gsp          , & ! precipitation rate of snow, grid-scale        (kg/m2*s)
                  prg_gsp          , & ! precipitation rate of graupel, grid-scale     (kg/m2*s)
!
                  tch              , & ! turbulent transfer coefficient for heat       ( -- )
                  tcm              , & ! turbulent transfer coefficient for momentum   ( -- )
                  tfv              , & ! laminar reduction factor for evaporation      ( -- )
!
                  sobs             , & ! solar radiation at the ground                 ( W/m2)
                  thbs             , & ! thermal radiation at the ground               ( W/m2)
                  pabs             , & !!!! photosynthetic active radiation            ( W/m2)
!
                  runoff_s         , & ! surface water runoff; sum over forecast       (kg/m2)
                  runoff_g         , & ! soil water runoff; sum over forecast          (kg/m2)
!
                  zshfl_s          , & ! sensible heat flux soil/air interface         (W/m2)
                  zlhfl_s          , & ! latent   heat flux soil/air interface         (W/m2)
                  zshfl_snow       , & ! sensible heat flux snow/air interface         (W/m2)
                  zlhfl_snow       , & ! latent   heat flux snow/air interface         (W/m2)
                  lhfl_bs          , & ! latent heat flux from bare soil evap.         (W/m2)
                  lhfl_pl          , & ! latent heat flux from plants                  (W/m2)
                  rstom            , & ! stomatal resistance                           ( s/m )
                  zshfl_sfc        , & ! sensible heat flux surface interface          (W/m2)
                  zlhfl_sfc        , & ! latent   heat flux surface interface          (W/m2)
                  zqhfl_sfc          & ! moisture      flux surface interface          (kg/m2/s)
                                     )

!-------------------------------------------------------------------------------
! Declarations for ICON (USE statements for COSMO)
!-------------------------------------------------------------------------------


  INTEGER (KIND=iintegers), INTENT(IN)  ::  &
                  icant,             & ! canopy type
                  ie,                & ! array dimensions
                  istartpar,         & ! start index for computations in the parallel program
                  iendpar,           & ! end index for computations in the parallel program
                  ke_soil, ke_snow,  &
                  ke_soil_hy           ! number of active soil moisture layers
  REAL    (KIND = ireals), DIMENSION(ke_soil+1), INTENT(IN) :: &
                  czmls                ! processing soil level structure
  LOGICAL, INTENT(IN) :: ldiag_tg      ! if .TRUE., use tgcom to diagnose t_g and snow-cover fraction
  INTEGER (KIND=iintegers), INTENT(IN)  :: &
                  inwp_turb            ! turbulence scheme number
  INTEGER (KIND=iintegers), INTENT(IN)  :: &
                  nclass_gscp          ! number of hydrometeor classes of grid scale microphysics
  REAL    (KIND = ireals), INTENT(IN)  ::  &
                  dt                   ! time step

  INTEGER (KIND=iintegers),DIMENSION(ie), INTENT(IN) :: &
                  soiltyp_subs         ! type of the soil (keys 0-9)                     --
  REAL    (KIND = ireals), DIMENSION(ie), INTENT(IN) :: &
                  plcov            , & ! fraction of plant cover                         --
                  rootdp           , & ! depth of the roots                            ( m  )
                  sai              , & ! surface area index                              --
                  tai              , & ! transpiration area index                        --
                  eai                  ! earth area (evaporative surface area) index     --
!  LOGICAL                , DIMENSION(ie), INTENT(IN) :: &
!                  llandmask           ! landpoint mask                                  --
  REAL    (KIND = ireals), DIMENSION(ie), INTENT(IN) :: &
                 rsmin2d               ! minimum stomata resistance                    ( s/m )
  REAL    (KIND = ireals), DIMENSION(ie), INTENT(IN) :: &
                  u                , & ! zonal wind speed                              ( m/s )
                  v                , & ! meridional wind speed                         ( m/s )
                  t                , & ! temperature                                   (  k  )
                  qv                   ! specific water vapor content                  (kg/kg)
  REAL    (KIND = ireals), DIMENSION(ie), INTENT(IN) :: &
                  p0                   !!!! base state pressure                        ( Pa )
  REAL    (KIND = ireals), DIMENSION(ie), INTENT(IN) ::    &
                  ps                   ! surface pressure                              ( pa  )
  REAL    (KIND = ireals), DIMENSION(ie), INTENT(INOUT) :: &
                  t_snow_now           ! temperature of the snow-surface (K)
  REAL    (KIND = ireals), DIMENSION(ie), INTENT(OUT) :: &
                  t_snow_new
  REAL    (KIND = ireals), DIMENSION(ie,0:ke_snow), INTENT(INOUT) :: &
                  t_snow_mult_now      ! temperature of the snow-surface               (  K  )
  REAL    (KIND = ireals), DIMENSION(ie,0:ke_snow), INTENT(OUT) :: &
                  t_snow_mult_new      ! temperature of the snow-surface               (  K  )
  REAL    (KIND = ireals), DIMENSION(ie), INTENT(INOUT) :: &
                  t_s_now              ! temperature of the ground surface             (  K  )
  REAL    (KIND = ireals), DIMENSION(ie), INTENT(OUT) :: &
                  t_s_new              ! temperature of the ground surface             (  K  )
  REAL    (KIND = ireals), DIMENSION(ie), INTENT(INOUT) ::&
                  t_g              , & ! weighted surface temperature                  (  K  )
                  qv_s                 ! specific humidity at the surface              (kg/kg)
  REAL    (KIND = ireals), DIMENSION(ie), INTENT(INOUT) :: &
                  w_snow_now       , & ! water content of snow                         (m H2O)
                  rho_snow_now         ! snow density                                  (kg/m**3)
  REAL    (KIND = ireals), DIMENSION(ie), INTENT(OUT) :: &
                  w_snow_new       , & ! water content of snow                         (m H2O)
                  rho_snow_new         ! snow density                                  (kg/m**3)
  REAL    (KIND = ireals), DIMENSION(ie,ke_snow), INTENT(INOUT) :: &
                  rho_snow_mult_now    ! snow density                                  (kg/m**3)
  REAL    (KIND = ireals), DIMENSION(ie,ke_snow), INTENT(OUT) :: &
                  rho_snow_mult_new    ! snow density                                  (kg/m**3)
  REAL    (KIND = ireals), DIMENSION(ie), INTENT(INOUT) :: &
                  h_snow               ! snow depth
  REAL    (KIND = ireals), DIMENSION(ie), INTENT(IN) :: &
                  h_snow_gp            ! grid-point averaged snow depth
  REAL    (KIND = ireals), DIMENSION(ie), INTENT(OUT) :: &
                  meltrate             ! snow melting rate
  REAL    (KIND = ireals), DIMENSION(ie), INTENT(INOUT) :: &
                  w_i_now              ! water content of interception water           (m H2O)
  REAL    (KIND = ireals), DIMENSION(ie), INTENT(OUT) :: &
                  w_i_new              ! water content of interception water           (m H2O)

  REAL    (KIND = ireals), DIMENSION(ie), INTENT(INOUT) :: &
                  w_p_now              ! water content of interception water pond      (m H2O)
  REAL    (KIND = ireals), DIMENSION(ie), INTENT(OUT) :: &
                  w_p_new              ! water content of interception water pond      (m H2O)
  REAL    (KIND = ireals), DIMENSION(ie), INTENT(INOUT) :: &
                  w_s_now              ! water content of interception snow water      (m H2O)
  REAL    (KIND = ireals), DIMENSION(ie), INTENT(OUT) :: &
                  w_s_new              ! water content of interception snow water

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
  REAL    (KIND = ireals), DIMENSION(ie), INTENT(IN) :: &
!                 t_2m             , & ! temperature in 2m                             (  K  )
                  u_10m            , & ! zonal wind in 10m                             ( m/s )
                  v_10m                ! meridional wind in 10m                        ( m/s )
  REAL    (KIND = ireals), DIMENSION(ie), INTENT(INOUT) :: &
                  freshsnow        , & ! indicator for age of snow in top of snow layer(  -  )
                  zf_snow              ! snow-cover fraction
  REAL    (KIND = ireals), DIMENSION(ie,ke_snow), INTENT(INOUT) :: &
                  wliq_snow_now    , & ! liquid water content in the snow              (m H2O)
                  wtot_snow_now        ! total (liquid + solid) water content of snow  (m H2O)
  REAL    (KIND = ireals), DIMENSION(ie,ke_snow), INTENT(OUT) :: &
                  wliq_snow_new    , & ! liquid water content in the snow              (m H2O)
                  wtot_snow_new        ! total (liquid + solid) water content of snow  (m H2O)
  REAL    (KIND = ireals), DIMENSION(ie,ke_snow), INTENT(INOUT) :: &
                  dzh_snow_now         ! layer thickness between half levels in snow   (  m  )
  REAL    (KIND = ireals), DIMENSION(ie,ke_snow), INTENT(OUT) :: &
                  dzh_snow_new         ! layer thickness between half levels in snow   (  m  )
  REAL    (KIND = ireals), DIMENSION(ie), INTENT(IN) ::    &
                  prr_con          , & ! precipitation rate of rain, convective        (kg/m2*s)
                  prs_con          , & ! precipitation rate of snow, convective        (kg/m2*s)
                  conv_frac        , & ! convective area fraction
                  prr_gsp          , & ! precipitation rate of rain, grid-scale        (kg/m2*s)
                  prs_gsp          , & ! precipitation rate of snow, grid-scale        (kg/m2*s)
                  prg_gsp              ! precipitation rate of graupel, grid-scale     (kg/m2*s)

  REAL    (KIND = ireals), DIMENSION(ie), INTENT(INOUT) :: &
                  tch              , & ! turbulent transfer coefficient for heat       ( -- )
                  tcm              , & ! turbulent transfer coefficient for momentum   ( -- )
                  tfv                  ! laminar reduction factor for evaporation      ( -- )

  REAL    (KIND = ireals), DIMENSION(ie), INTENT(IN) :: &
                  sobs             , & ! solar radiation at the ground                 ( W/m2)
                  thbs             , & ! thermal radiation at the ground               ( W/m2)
                  pabs                 !!!! photosynthetic active radiation            ( W/m2)

  REAL    (KIND = ireals), DIMENSION(ie), INTENT(INOUT) :: &
                  runoff_s         , & ! surface water runoff; sum over forecast       (kg/m2)
                  runoff_g             ! soil water runoff; sum over forecast          (kg/m2)
  REAL    (KIND = ireals), DIMENSION(ie), INTENT(OUT) :: &
                  zshfl_s          , & ! sensible heat flux soil/air interface         (W/m2)
                  zlhfl_s          , & ! latent   heat flux soil/air interface         (W/m2)
                  zshfl_snow       , & ! sensible heat flux snow/air interface         (W/m2)
                  zlhfl_snow           ! latent   heat flux snow/air interface         (W/m2)

  REAL    (KIND = ireals), DIMENSION(ie), INTENT(OUT) :: &
                  rstom            ! stomata resistance                                ( s/m )
  REAL    (KIND = ireals), DIMENSION(ie), INTENT(OUT) :: &
                  lhfl_bs          ! latent heat flux from bare soil evap.             ( W/m2)
  REAL    (KIND = ireals), DIMENSION(ie,ke_soil+1), INTENT(OUT) :: &
                  lhfl_pl          ! average latent heat flux from plants              ( W/m2)
  REAL    (KIND = ireals), DIMENSION(ie), OPTIONAL, INTENT(OUT) :: &
                  zshfl_sfc        , & ! sensible heat flux surface interface          (W/m2)
                  zlhfl_sfc        , & ! latent   heat flux surface interface          (W/m2)
!DR start
                  zqhfl_sfc            ! latent   heat flux surface interface          (W/m2)
!DR end

!--------------------------------------------------------------------------------
! TERRA Declarations

! New declaration for ICON

!------------------------------------------------------------------------------
! Subroutine arguments: None
! --------------------

  INTEGER (KIND=iintegers)              ::  &
    ierror                        ! error status variable

  CHARACTER (LEN=80)                    ::  &
    yerror
!
! Local parameters:
! ----------------

  REAL    (KIND=ireals   ), PARAMETER ::  &
    zepsi  = 1.0E-6_ireals , & ! security constant
    zalfa  = 1.0_ireals    , & ! degree of impliciteness (1: full implicit,
                               !    (0.5: Cranck-Nicholson)
    T_ref_ice = 0.1_ireals , & !degC Soil ice parameterization
    T_star_ice= 0.01_ireals, & !degC according to K. Schaefer and Jafarov, E.,2016, doi:10.5194/bg-13-1991-2016
    b_clay= -0.3_ireals , & 
    b_silt= -0.5_ireals , &
    b_sand= -0.9_ireals , &
    b_org = -1.0_ireals , &
! 
    rho_i  = 910._ireals       ! density of solid ice (soil model)  (kg/m**3)


! Local scalars:
! -------------

  INTEGER (KIND=iintegers) ::  &
!
!   Indices
!
!   nx             , & ! time-level for integration
    kso            , & ! loop index for soil moisture layers
    ksn            , & ! loop index for snow layers
    k              , & ! loop index for snow layers
    i,ic           , & ! loop index in x-direction
    icount_snow    , & ! Counter for snow
    icount_soil    , & ! "true" soil
    icount_rockice , & ! rock and ice points
#ifndef __ICON__
!    im1, jm1       , & ! i-1, j-1   ! must be removed completely
#endif
    jb             , & ! loop index for soil-type
    mstyp          , & ! soil type index
    msr_off        , & ! number of layers contributing to surface run off
    istarts        , & ! start index for x-direction
    iends          , & ! end   index for x-direction
    k10cm          , & ! index of half level closest to 0.1 m
    k100cm             ! index of half level closest to 1.0 m

  REAL    (KIND=ireals   ) ::  &
!
!   Timestep parameters
!
    zdt            , & ! integration time-step [s]
!   zdel_t_so      , & ! auxiliary variable
    zdtdrhw        , & ! timestep / density of water
    zrhwddt        , & ! density of water / timestep
    zroffdt        , & ! timestep for runoff calculation
    z1d2dt         , & ! 1./2*timestep
!
!   Connection to the atmosphere
!
    zuv            , & ! wind velocity in lowest atmospheric layer
    ztvs           , & ! virtual temperature at the surface
    zplow          , & ! pressure of lowest atmospheric layer
    zqvlow         , & ! specific humidity of lowest atmospheric layer
!em    zrss           , & ! ice covered part of grid element
    zrww           , & ! water covered part of grid element
!
!   Snow parameters
!
    zsnow_rate     , & ! rate of snow fall            [kg/m**2 s]
    zrain_rate     , & ! rate of rain fall            [kg/m**2 s]
    zrime          , & ! ground riming rate
    zdsn_new       , & ! snow age refresh increment   [-]
    zdsn_old       , & ! snow age decay increment     [-]
    zdz_snow(ie)    ! snow depth (not for snow heat content but snow
                       ! density calculation)

  REAL    (KIND=ireals   ) ::  &
!
!   Multi-snow layer parameters
    zsn_porosity   , & ! snow porosity                    (   -   )
    zp1            , &
    zfukt          , &
    zq0            , &
    zqbase(ie)  , &
    zdzh_old       , &
    zrefr(ie)   , & ! rate of liquid water refreezing
    zmelt(ie)   , & ! rate of snow melting
    ze_in          , &
    ze_out(ie)  , &
    zadd_dz        , &
    zrho_dry_old(ie)   , &
    zeta           , &
    zdens_old      , &
    zp(ie,ke_snow), &
    zcounter(ie), &
    ze_rad(ie)  , &
    zswitch(ie) , &
    tmp_num(ie) , &
    sum_weight(ie)      , &
    t_new  (ie,ke_snow) , &
    rho_new(ie,ke_snow) , &
    wl_new (ie,ke_snow) , &
    z_old  (ie,ke_snow) , &
    dz_old (ie,ke_snow) , &
    t_so_free_new(ie,0:ke_soil+1), &
    t_so_snow_new(ie,0:ke_soil+1), &
    sn_frac(ie), &
    zf_snow_lim(ie), &
    weight         , &
!
!   Plant parameters
    zbeta          , & ! reduction factor for evaporation
    zevap          , & ! auxiliary variable
    zalpha              ! NP89 bare soil evaporation

  REAL    (KIND=ireals   ) ::  &
!
!   Local scalars for BATS-scheme
!
    zpsi0=.01_ireals,& ! air entry potential at water saturation (m)
    zbf1           , & ! auxiliary variable
    zbf2           , & ! auxiliary variable
    zdmax          , & ! auxiliary variable
    zd             , & ! auxiliary variable
    zck            , & ! auxiliary variable
    z1             , & ! depth of half level closest to 0.1 m
    znull          , & ! depth of half level closest to 1.0 m
    zfqmax         , & ! maximum sustainable flux in uppermost soil layer
    zevapor        , & ! evaporation rate of bare soil in BATS-scheme
    zrla           , & ! atmospheric resistance
    zcatm          , & ! transfer function CA
    zpar           , & ! PAR (interim value)
    zf_wat         , & ! soil water function for stomatal resistance
    zf_tem         , & ! temperature function for stomatal resistance
    zf_sat         , & ! saturation deficit function for stomatal resistance
    zepsat         , & ! saturation vapour pressure at near surface temperature
    zepke          , & ! near surface water vapour pressure
    zrstom         , & ! stomata resistance
    zrveg          , & ! sum of atmospheric and stomata resistance
    zedrstom       , & ! 1./zrstom
    zustar         , & ! friction velocity
    ztrabpf        , & ! area average transpiration rate
    ! variables for non-uniform root-distribution
    ztrfr          , & ! fraction of total transpiration associated with each layer
    zrootfr        , & ! normalized root fraction in each layer (=1 for z_soil=0.0)
    zrootdz        , & ! normalized root fraction * min(layer thickness,rootlength-top_of_layer)
!
!   Soil parameters
!
    zice           , & ! indicator of soil type ice
    zwqg           , & ! mean of fcap and pwp
    z4wdpv         , & ! 4*zwqg/porv
!   zroot          , & ! fraction of layer filled by roots (Bucket scheme)
    zropart        , & ! fraction of layer filled by roots (Dickinson scheme)
    zrootfc        , & ! distributes transpiration to soil layers
!   zr_root        , & ! zroot/zbwt (Bucket scheme)
    ze_sum             ! sum of all contributions to evapotranspiration

  REAL    (KIND=ireals   ) ::  &
!
!   Thermal parameters
!
!    ztgt0          , & ! Indicator T_g > T_0
    zgstr          , & ! downward longwave radiation
    zalas          , & ! heat conductivity of snow
    zrnet_snow     , & ! net radiation at snow surface
    zfak           , & ! utility variable for implicit snow temperature forecast
    ztsnow_im      , & ! utility variable for implicit snow temperature forecast
    ztsnew         , & ! preliminary value of snow surface temperature
    ztsnownew      , & ! preliminary value of snow surface temperature
    zfor_snow      , & ! total forcing at snow surface
    zfr_melt       , & ! melting snow fraction
!    zdwsnm         , & ! utility variable for snow melt determination
    zdwgme         , & ! utility variable for snow melt determination
    zdelt_s        , & ! utility variable for snow melt determination
    ze_avail       , & ! utility variable for snow melt determination
    ze_total           ! utility variable for snow melt determination

  ! Some of these values have to be fields for the multi-layer snow model
  REAL    (KIND=ireals   ) ::  &
    zalas_mult    (ie,ke_snow),    & ! heat conductivity of snow
    ztsnownew_mult(ie,0:ke_snow),  & ! preliminary value of snow surface temperature
    zextinct      (ie,ke_snow),    & ! solar radiation extinction coefficient in snow (1/m)
    zfor_snow_mult(ie)               ! total forcing at snow surface

! REAL    (KIND=ireals   ) ::  &
!   zfor_total    (ie)               ! total forcing at surface

  REAL    (KIND=ireals   ) ::  &
!
!   Interception variables
!
    zwinstr  (ie  )       , & ! preliminary value of interception store
    zinfmx   (ie  )       , & ! maximum infiltration rate
    zwimax   (ie  )       , & ! maximum interception store
    zvers    (ie  )       , & ! water supply for infiltration
    zwisnstr (ie  )      , & ! water content of snow interception store (t+1) (mH2O)
    zwpnstr  (ie  )      , & ! water content of pond store (t+1) (mH2O)
    zfd_wi   (ie  )      , & ! surface fraction of intercepted water
    zf_pd    (ie  )      , & ! surface fraction covered by pond
    zwisn    (ie  )      , & ! water content of snow interception store (mH2O)
    zwpn     (ie  )      , & ! water content of pond
    zewi (ie  )   ,       & ! evaporation from interception store
    zepd (ie  )   ,       & ! evaporation from pond
    zept (ie  )   ,       & ! potential evaporation
    zesn (ie  )   ,       & ! evaporation from snow
    zdrr (ie  )   ,       & ! formation of dew
    zrrs (ie  )             ! formation of rime
!
  REAL    (KIND=ireals   ) ::  &
!   Hydraulic parameters
    zalf           , & ! utility variable
    zro_inf        , & ! surface runoff
    zdwsndtt       , & ! time rate of change of snow store
    zdwisndtt      , & ! time rate of change of snow store
    zdwpdtt        , & ! time rate of change of pond store
    zwsnstr        , & ! preliminary value of snow store
    zdwseps        , & ! artificial change of small amount of snow store
    zdwidtt        , & ! time rate of change of interception store
    zdwieps        , & ! artificial change of small amount of interception store
    zro_wi         , & ! surface runoff due to limited infiltration capacity
    zro_sfak       , & ! utility variable
    zro_gfak       , & ! utility variable
    zfmb_fak       , & ! utility variable
    zdwg           , & ! preliminary change of soil water content
    zwgn           , & ! preliminary change of soil water content
    zredfu         , & ! utility variable for runoff determination
    zro            , & ! utility variable for runoff determination
    zro2           , & ! utility variable for runoff determination
    zkorr          , & ! utility variable for runoff determination
    zfr_ice        , & ! reduction factor for water transport
    zfr_ice_free   , & ! reduction factor for water transport
    zwso_new       , & ! preliminary value of soil water content
!    zwsnew         , & ! preliminary value of snow water equivalent
    zw_ovpv        , & ! utility variable
    zfcorr_wi      , & ! utility variable
    zpercmax       , & ! utility variable
!
!   Implicit solution of thermal and hydraulic equations
!
    zakb           , & ! utility variable
    zakb1          , & ! utility variable
    zakb2          , & ! utility variable
    zzz            , & ! utility variable
    z1dgam1        , & ! utility variable
    zredm          , & ! utility variable
    zredm05        , & ! utility variable
    zredp05        , & ! utility variable
    zgam2m05       , & ! utility variable
    zgam2p05           ! utility variable

  REAL    (KIND=ireals   ) ::  &
!
!   Statement functions
!
!!$    zsf_heav       , & ! Statement function: Heaviside function
!!$    zsf_psat_iw    , & ! Saturation water vapour pressure over ice or water
!!$                       ! depending on temperature "zstx"
!!$    zsf_qsat       , & ! Specific humidity at saturation pressure
!!$                       ! (depending on the saturation water vapour pressure
!!$                       !  "zspsatx" and the air pressure "zspx")
!!$    zsf_dqvdt_iw   , & ! Statement function: First derivative of specific
!!$                       ! saturation humidity
!!$                       ! with respect to temperature (depending on temperature
!!$                       ! "zstx" and saturation specific humidity pressure
!!$                       ! "zsqsatx")
!!$    watrcon_RT, watrcon_BC , & !
!!$    watrdiff_RT,watrdiff_BC, & !
!!$    ks,ds, lw,tw,m,a,kw1,dw1,pv,adp ,&
!!$    zstx           , & ! dummy argument for Stmt. function
!!$    zspx           , & ! dummy argument for Stmt. function
!!$    zspsatx        , & ! dummy argument for Stmt. function
!!$    zsqsatx        , & ! dummy argument for Stmt. function
    z2iw           , & ! dummy argument for Stmt. function
    z4iw           , & ! dummy argument for Stmt. function
    z234iw         , & ! dummy argument for Stmt. function
    zqs            , & ! saturation humidty at T_s and ps
    zdqs           , & ! derivative of zqs with respect to T_s
    zqsnow         , & ! saturation humidty at T_snow and ps
    zdqsnow        , & ! derivative of zqs with respect to T_snow
!   Local scalars for transfer coefficient limitation calculations and for
!   treatment of possible problems at lateral boundaries
    zch_soil = 1.400E06_ireals ,& ! approximate average heat capacity of soil (J/m**3 K)
    zlim_dtdt = 2.5_ireals   ,&! maximum allowed temperature increment (K) in
                               ! one time step in  uppermost soil layer
    zdT_snow                            ,& ! maximum possible increment in snow layer before
                                           ! melting sets in
    zg1       (ie                  ) ,& ! heat flux between two uppermost soil layers (W/m**2)
    zlhfl     (ie                  ) ,& ! estimated       latent    heat flux surface (W/m**2)
    zshfl     (ie                  ) ,& ! estimated       sensible  heat flux surface (W/m**2)
    zthfl     (ie                  ) ,& ! estimated total turbulent heat flux surface (W/m**2)
    zradfl    (ie                  ) ,& ! total radiative flux at surface             (W/m**2)
    ze_melt   (ie                  ) ,& ! energy required for melting of snow         (J/m**2)
    zch_snow  (ie                  ) ,& ! heat capacity of snow * w_snow * Rho_w      (J/m**2 K)
    zeb1      (ie                  ) ,& ! estimated energy budget of first soil layer (W/m**2)
    ztchv     (ie                  ) ,& ! transfer coefficient multiplied by wind velocity
    ztchv_max (ie                  ) ,& ! dito, but upper limit
    zrho_atm  (ie                  ) ,& ! air density of lowest atmospheric layer     (kg/m**3)
    zdt_atm   (ie                  ) ,& ! temperature difference between lowest
    zdq_atm   (ie                  ) ,& ! dito for moisture (kg/kg)
    !additional variables for root distribution
    zroota    (ie                  ) ,& ! root density profile parameter (1/m)
    zwrootdz  (ie         , ke_soil) ,& ! mean water content over root depth weighted by root density
    zrootdz_int (ie                ) ,& ! parameter needed to initialize the root density profile integral
    zwrootdz_int(ie                ) ,& ! parameter needed to initialize the root water content integral
!DR start
    zqhfl_s    (ie                 ) ,& ! moisture flux at soil/air interface
    zqhfl_snow (ie                 )    ! moisture flux at snow/air interface
!DR end


  INTEGER m_limit                          ! counter for application of limitation
  CHARACTER *7 yhc                         ! heating or cooling indicator


  REAL    (KIND=ireals   ) ::  &
!
!   Local scalars for water content dependent freezing/melting
!
    zliquid        , & ! utility variable
    zxx            , & ! utility variable
    znen           , & ! utility variable
    zargu          , & ! utility variable
!
!   Water transport
!
    zice_fr_ksom1  , & ! fractional ice content of actual layer - 1
    zice_fr_kso    , & ! fractional ice content of actual layer
    zice_fr_ksop1  , & ! fractional ice content of actual layer + 1
    zlw_fr_ksom1   , & ! fractional liquid water content of actual layer - 1
    zlw_fr_ksom1_new,& ! fractional liquid water content of actual layer -1
    zlw_fr_kso     , & ! fractional liquid water content of actual layer
    zlw_fr_kso_new , & ! fractional liquid water content of actual layer
    zlw_fr_ksop1   , & ! fractional liquid water content of actual layer + 1
    zlw_fr_ksom05  , & ! fractional liquid water content of actual layer - 1/2
    zdlw_fr_ksom05 , & ! hydraulic diffusivity coefficient at half level above
    zklw_fr_ksom05 , & ! hydraulic conductivity coefficient at half level above
    zklw_fr_kso_new, & ! hydraulic conductivity coefficient at main level
                       !    for actual runoff_g
    zdlw_fr_kso    , & ! hydraulic diff coefficient at main level
    zklw_fr_kso    , & ! hydraulic conductivity coefficient at main level
    zklw_fr_ksom1  , & ! hydraulic conductivity coefficient at main level above
    zlw_fr_ksop05  , & ! fractional liquid water content of actual layer + 1/2
    zdlw_fr_ksop05 , & ! hydraulic diffusivity coefficient at half level below
    zklw_fr_ksop05 , & ! hydraulic conductivity coefficient at half level below
    zinf           , & ! infiltration
!
!   Snow density
    ztau_snow      , & ! 'ageing constant' for snow density          (-)
    zrhosmax       , & ! temperature-dependent target density for snow ageing
    zrho_snowe     , & ! updated density of existing snow ('ageing') kg/m**3
    zrho_snowf     , & ! density of fresh snow                       kg/m**3
    znorm          , & ! normalisation factor for weighted snow density mH2O
!
!
!   Freezing/melting of soil water/ice:
!
    zdelwice       , & ! amount of melted soil ice/frozen soil water
    zdwi_scal      , & ! time scale parameter for freezing/melting soil water
    ztx            , & ! water content dependent freezing/melting temperature
    zd1, zd2, zd3, zd4 ! auxiliary variables

! Local (automatic) arrays:
! -------------------------

  REAL    (KIND = ireals) :: &
!
! Two-time level variables exchange with interface
!
    h_snow_now (ie)    , & ! snow height  (m)
    h_snow_new (ie)    , & ! snow height  (m)
!
!
! Model geometry
!
    zmls     (ke_soil+1)  , & ! depth of soil main level
    zzhls    (ke_soil+1)  , & ! depth of the half level soil layers in m
    zdzhs    (ke_soil+1)  , & ! layer thickness between half levels
    zdzms    (ke_soil+1)  , & ! distance between main levels
    zdz_snow_fl(ie)    , & ! snow depth for snow temperature gradient
!
! Multi-layer snow model
    zhh_snow (ie,ke_snow)    , & ! depth of the half level snow layers
    zhm_snow (ie,ke_snow)    , & ! depth of snow main levels
    zdzh_snow(ie,ke_snow)    , & ! layer thickness between half levels
    zdzm_snow(ie,ke_snow)    , & ! distance between main levels
!
! External parameters
!
    zbwt     (ie)      , & ! root depth (with artificial minimum value)
    zrock    (ie)      , & ! ice/rock-indicator: 0 for ice and rock
    zsandf   (ie)      , & ! mean fraction of sand (weight percent)
    zclayf   (ie)      , & ! mean fraction of clay (weight percent)
    zsiltf   (ie)      , & ! mean fraction of clay (weight percent)
    zb_por   (ie)      , & ! pore size distribution index
    zpsis    (ie)      , & ! air entry potential (m)
    zw_m_org  (ie)     , &  ! maximum of  liquid water content organic
    zw_m_soil (ie)     , &  ! maximum of  liquid water content   mineral soil
    zw_m_up   (ie)       , &  ! maximum of  liquid water content   at temp -3 degC
    zw_m_low  (ie,ke_soil+1), &  ! maximum of  liquid water content   at temp -40 degC
    t_zw_up          , &  ! temp -3 degC
    t_zw_low         , &  ! temp -40 degC
    zw_m     (ie)         ! maximum of liquid water content  (m)

  INTEGER  (KIND=iintegers ) ::  &
    m_styp   (ie)      , & ! soil type
    melt_list(ie)      , & ! list of melting snow points
    soil_list(ie)      , & ! list of "true" soil points
                           ! i.e. no ice no rock
    rockice_list(ie)       ! list of rock and ice points

  REAL    (KIND=ireals   ) ::  &
!
! Connection to the atmosphere
!
    zrr      (ie)      , & ! total rain rate including formation of dew
    zrs      (ie)      , & ! total snow rate including formation of rime
    zesoil   (ie)      , & ! evaporation from bare soil
    zrhoch   (ie)      , & ! transfer coefficient*rho*g
    zth_low  (ie)      , & ! potential temperature of lowest layer
    zf_wi    (ie)      , & ! surface fraction covered by interception water
    ztmch    (ie)      , & ! heat transfer coefficient*density*velocity
    zep_s    (ie)      , & ! potential evaporation for t_s
    zep_snow (ie)      , & ! potential evaporation for t_snow
    zverbo   (ie)      , & ! total evapotranspiration
    zversn   (ie)      , & ! total evaporation from snow surface
    zthsoi   (ie)      , & ! thermal flux at soil surface
    zthsnw   (ie)      , & ! thermal flux at snow surface
    zfor_s   (ie)      , & ! total forcing at soil surface
    zgsb     (ie)      , & ! heat-flux through snow
    zrnet_s  (ie)      , & ! net radiation
    zsprs    (ie)      , & ! utility variable
!
! Tendencies
!
    zdwidt   (ie)      , & ! interception store tendency
    zdwsndt  (ie)      , & ! snow water content tendency
    zdtsdt   (ie)      , & ! tendency of zts
    zdtsnowdt(ie)      , & ! tendency of ztsnow
    zdtsnowdt_mult(ie,0:ke_snow)      , & ! tendency of ztsnow
    zdwgdt   (ie,ke_soil)  ! tendency of water content [kg/(m**3 s)]

  REAL    (KIND=ireals   ) ::  &
!
! Soil and plant parameters
!
    zroc(ie,ke_soil+1) , & ! heat capacity
    zfcap    (ie,ke_soil+1)      , & ! field capacity of soil
    zadp     (ie,ke_soil+1)      , & ! air dryness point
    zporv    (ie,ke_soil+1)      , & ! pore volume (fraction of volume)
    zdlam    (ie)      , & ! heat conductivity parameter
    zdw      (ie,ke_soil+1)      , & ! hydrological diff.parameter
    zdw1     (ie,ke_soil+1)      , & ! hydrological diff.parameter
    zkw      (ie,ke_soil+1)      , & ! hydrological cond.parameter
    zkw1     (ie,ke_soil+1)      , & ! hydrological cond.parameter
    zik2     (ie)      , & ! minimum infiltration rate
    zpwp     (ie,ke_soil+1)      , & ! plant wilting point  (fraction of volume)
    ztlpmwp  (ie)      , & ! turgor-loss-point minus plant wilting point
    zedb     (ie)      , & ! utility variable
    zaa      (ie)      , & ! utility variable
!
! Hydraulic variables
!
    ztrang(ie,ke_soil) , & ! transpiration contribution by the different layers
    ztrangs(ie)        , & ! total transpiration (transpiration from all
                              !    soil layers)
    zwin     (ie)      , & ! water cont. of interception store   (m H20)
    zwsnow   (ie)      , & ! snow water equivalent               (m H20)
    zwsnew   (ie)      , & ! snow water equivalent               (m H20)
    zdwsnm   (ie)      , & ! utility variable for snow melt determination
    zw_fr    (ie,ke_soil+1),&!fractional total water content of soil layers
    zinfil   (ie)      , & ! infiltration rate
    zlw_fr   (ie,ke_soil+1), & ! fractional liqu. water content of soil layer
    ziw_fr   (ie,ke_soil+1), & ! fractional ice content of soil layer
    zwsnn    (ie)          , & ! new value of zwsnow
    zflmg    (ie,ke_soil+1), & ! flux of water at soil layer interfaces
    zrunoff_grav(ie,ke_soil+1) ! main level water gravitation
                                  !    flux (for runoff calc.)

  REAL    (KIND=ireals   ) ::  &
!
! Local arrays for BATS-scheme
!
    zk0di    (ie)      , & ! surface type dependent parameter
    zbedi    (ie)      , & ! surface type dependent parameter
    zsnull   (ie)      , & ! mean relative(zporv) water fraction of the
                              !    so-called total active soil layer (above 1m)
    zs1      (ie)      , & ! mean relative (zporv) water fraction  of
                              !    layers above 0.1m
    zf_rad   (ie)      , & ! radiation function for stomata resistance
    ztraleav (ie)      , & ! transpiration rate of dry leaves
!
! Root functions
!
    zwroot   (ie)      , & ! mean water content over root depth
    zropartw(ie,ke_soil)   ! fraction of layer filled by roots * w_fr

  REAL    (KIND=ireals   ) ::  &
!
! Thermal variables
!
    zts      (ie)      , & ! soil surface temperature
    ztsnow   (ie)      , & ! snow surface temperaure
    ztsnow_mult   (ie,0:ke_snow)      , & ! snow surface temperaure
! HEATCOND (soil moisture dependent heat conductivity of the soil)
    zalamtmp (ie,ke_soil),& ! heat conductivity
    zalam    (ie,ke_soil),&! heat conductivity
    zrocg    (ie,ke_soil+1),& ! total volumetric heat capacity of soil
    zrocg_soil(ie,ke_soil+1),& ! volumetric heat capacity of bare soil
    zrocs    (ie)      , & ! heat capacity of snow
    ztsn     (ie)      , & ! new value of zts
    ztsnown  (ie)      , & ! new value of ztsnow
    ztsnown_mult  (ie,0:ke_snow)      , & ! new value of ztsnow
    znlgw1f  (ke_soil)    , & ! utility variable
    zaga(ie,0:ke_soil+ke_snow+1),& ! utility variable
    zagb(ie,0:ke_soil+ke_snow+1),& ! utility variable
    zagc(ie,0:ke_soil+ke_snow+1),& ! utility variable
    zagd(ie,0:ke_soil+ke_snow+1),& ! utility variable
    zage(ie,0:ke_soil+ke_snow+1),& ! utility variable
!
! Auxiliary variables
!
!    zw_snow_old(ie)    , & !
    zdqvtsnow(ie)      , & ! first derivative of saturation specific humidity
                              !    with respect to t_snow
    zrho_snow(ie)      , & ! snow density used for computing heat capacity and conductivity
    zts_pm   (ie)      , & ! indicator zts > < T_melt
    ztfunc   (ie)      , & ! smoothed transition function between T_melt and T_melt+2K
    ztsnow_pm(ie)          ! indicator ztsnow > < T_melt


  LOGICAL     :: &
!
    ldebug                , & !
    limit_tch (ie)         ! indicator for flux limitation problem

  ! ground water as lower boundary of soil column
  REAL    (KIND=ireals) ::  &
    zdelta_sm, zdhydcond_dlwfr ! zklw_fr_kso, zklw_fr_ksom1, zdlw_fr_kso

  ! HEATCOND
  REAL    (KIND=ireals) ::          &
    hzalam   (ie,ke_soil+1),     & ! heat conductivity
    zthetas, zlamli, zlamsat, zlams, rsandf, zlamq, zlam0, zrhod, zlamdry,  &
    zlamdry_soil, zsri, zKe, zthliq, zlamic, zlamdry_c1, zlamdry_c2

  ! For performance improvement
  REAL    (KIND=ireals) :: ln_2, ln_3, ln_10, ln_006

!  INTEGER (KIND=iintegers) :: i_loc, isub

!- End of header
!==============================================================================

  zaga = 1.0_ireals
  zagb = 1.0_ireals
  zagc = 1.0_ireals

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


! Declaration of STATEMENT-FUNCTIONS

!!$  zsf_heav     (zstx                    ) = 0.5_ireals+SIGN( 0.5_ireals, zstx )
!!$  zsf_psat_iw  (zstx,z2iw   ,z4iw       )                                     &
!!$                   = b1*EXP(z2iw*(zstx - b3)/(zstx - z4iw))

!!$  zsf_qsat     (zspsatx, zspx           )                                     &
!!$                   = rdv*zspsatx/(zspx-o_m_rdv*zspsatx)

!!$  zsf_dqvdt_iw (zstx,zsqsatx,z4iw,z234iw)                                     &
!!$                   = z234iw*(1._ireals+rvd_m_o*zsqsatx)*zsqsatx/(zstx-z4iw)**2

!!$  watrcon_RT   (ks,lw,kw1,pv,adp   )  =  ks*EXP(kw1*(pv-lw)/(pv-adp))
!!$
!!$  watrdiff_RT  (ds,lw,dw1,pv,adp   )  =  ds*EXP(dw1*(pv-lw)/(pv-adp))
!!$
!!$  watrcon_BC(ks,tw,m) = ks*tw**(5._ireals/2._ireals+2._ireals/ &
!!$                           ((1._ireals/(1._ireals-m))-1._ireals))
!!$  watrdiff_BC(ks,tw,m,a,pv,adp)=ks/(a*((1._ireals/(1._ireals-m))-1._ireals)*(pv-adp))*&
!!$            tw**(3._ireals/2._ireals + 1._ireals/((1._ireals/(1._ireals-m))-1._ireals))


!------------------------------------------------------------------------------
! Section I.1: Initializations
!------------------------------------------------------------------------------

! Horizontal domain for computation
  istarts = istartpar
  iends   = iendpar

!>JH
!  prg_gsp=0._ireals ! graupel not implemented yet
!<JH

  ierror = 0
  yerror = '        '

! select timelevel and timestep for calculations
! In the new formulation a 2 timelevel scheme is used for the soil model.
! Therefore it is not necessary any more to distinguish between Leapfrog
! and Runge-Kutta scheme.
! But the first timestep must not be done with dt/2, as in the Leapfrog
! scheme
!  nx  = nnow

!!$  IF ( (ntstep == 0) .AND. (.NOT. l2tls) ) THEN
!!$    ! use the original dt and not dt/2
!!$    zdt = 2._ireals * dt
!!$  ELSE
  zdt = dt
!!$  ENDIF

  ! time step for run-off computation
  zroffdt = zdt


! Computation of derived constants

  zrhwddt = rho_w/zdt     ! density of liquid water/timestep
  zdtdrhw = zdt/rho_w     ! timestep/density of liquid water

  zdwi_scal = zdt/1800._ireals ! time scale parameter for freezing/melting soil water

! time constant for infiltration of water from interception store
  ctau_i        = MAX(ctau_i,zdt)

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

  ln_10 = LOG(10._ireals)


!---loop over tiles---
!DO ns=nsubs0,nsubs1
!---------------------


! Prepare basic surface properties (for land-points only)

  icount_soil     =0
  icount_rockice  =0
  soil_list(:)    =0
  rockice_list(:) =0
  melt_list(:)    =0
  meltrate(:)     = 0._ireals


! Temperaturedifference for liquid water content in frozen soil at -40 degC
!  J. Helmert: Soil ice parameterization according to K. Schaefer and Jafarov, E.,2016,
!                                                    doi:10.5194/bg-13-1991-2016
  t_zw_up  = 270.15_ireals ! temp -3 degC
  t_zw_low = 233.15_ireals ! temp -40 degC


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
    ! ensure that glaciers are covered with at least 1 m of snow
    IF (mstyp == 1) h_snow(i) = MAX(1._ireals, h_snow(i))
    zdw   (i,:)  = cdw0  (mstyp)
    zdw1  (i,:)  = cdw1  (mstyp)
    ! SBr, CHa: In order no to allow moisture sink to lower levels we have to
    ! set hydraulic conductivity to zero
    !zkw   (i,:)  = ckw0  (mstyp)
    zkw   (i,:)  = ckw0  (mstyp) * 0.0e-4_ireals
    zkw1  (i,:)  = ckw1  (mstyp)
    zik2  (i)  = cik2  (mstyp)
    zporv(i,:)  = cporv(mstyp)              ! pore volume
    zpwp (i,:)  = cpwp (mstyp)              ! plant wilting point
    zadp (i,:)  = cadp (mstyp)              ! air dryness point
    zfcap(i,:)  = cfcap(mstyp)              ! field capacity
    zrock(i)  = crock(mstyp)              ! EQ 0 for Ice and Rock EQ 1 else
    zrocg(i,:)  = crhoc(mstyp)              ! heat capacity
    zrocg_soil(i,:)  = crhoc(mstyp)              ! heat capacity
    zalam(i,:)  = cala0(mstyp)              ! heat conductivity parameter
    zdlam(i)  = cala1(mstyp)-cala0(mstyp) ! heat conductivity parameter
    zbwt(i)   = MAX(0.001_ireals,rootdp(i))! Artificial minimum value
                                                ! for root depth
    zroota(i) = 3._ireals/zbwt(i)       ! root density profile parameter (1/m)
                                                ! zroota=0. creates the original TERRA_LM
                                                ! version with constant root density

    ! New arrays for BATS-scheme
    zk0di(i)  = ck0di(mstyp)              !
    zbedi(i)  = cbedi(mstyp)              !
  ENDDO

  ! Arrays for soil water freezing/melting
  zd = LOG((T_ref_ice-(t_zw_low-t0_melt))/T_star_ice)
  zd1 = EXP(b_sand*zd)
  zd2 = EXP(b_clay*zd)
  zd3 = EXP(b_silt*zd)
  zd4 = EXP(b_org*zd)
  DO i = istarts, iends
    mstyp       = soiltyp_subs(i)        ! soil type
    zsandf(i)   = csandf(mstyp)
    zclayf(i)   = cclayf(mstyp)
    zsiltf(i)   = 100._ireals -csandf(mstyp)-cclayf(mstyp) ! Residuum of sand and clay
    zpsis(i)  = -zpsi0 * EXP(ln_10*(1.88_ireals-0.013_ireals*zsandf(i)))
    zb_por(i) = 2.91_ireals + .159_ireals*zclayf(i)
    zedb(i)   = 1._ireals/zb_por(i)
    zaa(i)    = g*zpsis(i)/lh_f
    ! Liq. water content at -3 degC
    zw_m_up(i) = EXP(-zedb(i)*LOG((t_zw_up - t0_melt)/(t_zw_up*zaa(i))) ) ! Without zporv(i,kso)*zdzhs(kso)!!
    ! Determine liq. water content at -40 degC
    !  J. Helmert: Soil ice parameterization according to K. Schaefer and Jafarov, E.,2016,
    !                                                    doi:10.5194/bg-13-1991-2016
    zw_m_soil(i) = 0.01_ireals*(zsandf(i)*zd1 + zclayf(i)*zd2 + zsiltf(i)*zd3)
    zw_m_org(i) = zd4
  ENDDO


  ! Set three-dimensional variables
  DO kso = 1, ke_soil
      DO i = istarts, iends
        !fc=2 1/m Exponential Ksat-profile decay parameter,see Decharme et al. (2006)
        zkw   (i,kso) = zkw   (i,kso)*EXP(-2._ireals*(zmls(kso)-rootdp(i)))

        ! Scale soil heat capacity with organic fraction -> Chadburn et al., 2015
        IF(zmls(kso) < rootdp(i)) THEN
          zzz = plcov(i)*(rootdp(i)-zmls(kso))/rootdp(i)
          zrocg(i,kso)=(1._ireals-zzz)*zrocg_soil(i,kso)+zzz*0.58E+06_ireals
          !  J. Helmert: Soil ice parameterization according to K. Schaefer and Jafarov, E.,2016,
          !  Organic fraction                                        doi:10.5194/bg-13-1991-2016
          zw_m_low(i,kso) = zporv(i,kso)*zdzhs(kso)*(zzz*zw_m_org(i) + (1._ireals-zzz)*zw_m_soil(i))
        ELSE
          zw_m_low(i,kso) = zporv(i,kso)*zdzhs(kso)*zw_m_soil(i)
        END IF
      ENDDO
  END DO


! For ntstep=nstart : Some preparations
! =====================================

!>JH  IF (ntstep == nstart) THEN
!   Determine constants clgk0 for BATS-scheme
  DO jb       = 1, 10
    clgk0(jb) = LOG(MAX(zepsi,ck0di(jb)/ckrdi))/ln_10
  END DO
!<JH  ENDIF


  DO i = istarts, iends
      ztrangs(i)      = 0.0_ireals
      zsnull(i)       = 0.0_ireals
      zs1(i)          = 0.0_ireals
      zwroot(i)       = 0.0_ireals
      zdwidt (i)      = 0.0_ireals         ! Initialisation of all
      zdwsndt(i)      = 0.0_ireals         ! evaporation quantities
      zesoil (i)      = 0.0_ireals         ! to be equal to zero
      zrr    (i)      = 0.0_ireals         ! in first part formation of dew
      zrs    (i)      = 0.0_ireals         ! in first part formation of rime
      zw_fr(i,ke_soil+1)  = w_so_now(i,ke_soil+1)/zdzhs(ke_soil+1)
      lhfl_bs(i)      = 0.0_ireals
      lhfl_pl(i,:)    = 0.0_ireals
      rstom  (i)      = 0.0_ireals
  END DO

  DO kso   = 1, ke_soil
      DO i = istarts, iends
        zw_fr(i,kso)    = w_so_now(i,kso)/zdzhs(kso)
        ztrang (i,kso)  = 0._ireals
        zropartw(i,kso) = 0._ireals
      END DO
  END DO

  !>JH WRITE(*,*) 'zw_fr: ',zw_fr(1,1,:)

! Determine the layer indices and total depth for water content averaging
! (in soil evaporation after Dickinson)
  k10cm  = 1
  k100cm = 1
  z1     = zzhls(1)
  znull  = zzhls(1)
  DO kso = 1, ke_soil
    IF (zmls(kso).le.0.1_ireals) THEN
      z1      = zzhls(kso)
      k10cm   = kso
    END IF
    IF (zmls(kso).le.1.0_ireals) THEN
      znull   = zzhls(kso)
      k100cm  = kso
    END IF
  END DO

! Determine average soil water content over 10 and 100 cm, respectively.
  DO kso   = 1, k100cm
      DO i = istarts, iends
!        IF (llandmask(i)) THEN   ! for land-points only
          IF (kso.le.k10cm) zs1(i) = zs1(i) + w_so_now(i,kso)
          zsnull(i)   = zsnull(i) + w_so_now(i,kso)
!        END IF
      END DO
  END DO

! JH No soil moisture for Ice and Rock
  DO kso   = 1, ke_soil+1
!CDIR NODEP,VOVERTAKE,VOB
    DO ic = 1, icount_rockice
      i=rockice_list(ic)
      w_so_now(i,kso)         = 0._ireals
      w_so_ice_now(i,kso)     = 0._ireals
      w_so_new(i,kso)         = w_so_now(i,kso)
      w_so_ice_new(i,kso)     = w_so_ice_now(i,kso)
    END DO

!CDIR NODEP,VOVERTAKE,VOB
    DO ic = 1, icount_soil
      i=soil_list(ic)
      w_so_new(i,kso)         = w_so_now(i,kso)
      w_so_ice_new(i,kso)     = w_so_ice_now(i,kso)
    END DO
  END DO

  ! Decide which snow density is used for computing the heat capacity
  IF (l2lay_rho_snow) THEN
    zrho_snow(istarts:iends) = rho_snow_mult_now(istarts:iends,1)
  ELSE
    zrho_snow(istarts:iends) = rho_snow_now(istarts:iends)
  ENDIF

  IF (itype_heatcond >= 2) THEN

! heat conductivity dependent on actual soil water content
! based on Peters-Lidard et al. (1998) and Johansen (1975),
! see also Block, Alexander (2007), Dissertation BTU Cottbus

    zlamli  = LOG(0.57_ireals)       ! LOG(thermal conductivity of water)
    zlamic  = LOG(2.2_ireals)        ! LOG(thermal conductivity of ice)
    zlamq   = LOG(7.7_ireals)        ! LOG(thermal conductivity of quartz)
    ln_2    = LOG(2._ireals)
    ln_3    = LOG(3._ireals)
    ln_006  = LOG(0.06_ireals)

    ! tuning constants for dry thermal conductivity formula
    IF (itype_evsl == 4) THEN
      zlamdry_c1 = 437._ireals
      zlamdry_c2 = 0.901_ireals
    ELSE
      zlamdry_c1 = 64.7_ireals
      zlamdry_c2 = 0.947_ireals
    ENDIF

    DO kso = 1, ke_soil+1
      DO i = istarts, iends
        zthetas = zporv(i,kso)                                 ! porosity
        zthliq  = zthetas - w_so_ice_now(i,kso)/zdzhs(kso) ! unfrozen volume fraction

        rsandf = zsandf(i)/100._ireals                     ! quartz content

        if (rsandf >= 0.2_ireals)  zlam0 = ln_2      ! LOG(thermal conductivity non-quartz)
        if (rsandf <  0.2_ireals)  zlam0 = ln_3

  ! saturated thermal conductivity

        zlams = zlamq*rsandf + zlam0*(1._ireals-rsandf)  ! LOG(solids thermal conductivity)

        zlamsat = EXP(zlams*(1.0_ireals-zthetas) + zlamic*(zthetas-zthliq) + zthliq*zlamli)

  ! dry thermal conductivity

        zrhod   = 2700.0_ireals*(1.0_ireals-zthetas)       ! dry density
        zlamdry_soil = ( 0.135_ireals*zrhod + zlamdry_c1 )                     &
                / ( 2700.0_ireals - zlamdry_c2*zrhod )
        ! missing: crushed rock formulation for dry thermal conductivity (see PL98)
  ! Scale zlamdry with organic fraction
        IF(zmls(kso) < rootdp(i)) THEN
          zzz = plcov(i)*(rootdp(i)-zmls(kso))/rootdp(i)
          zlamdry = EXP(LOG(zlamdry_soil)*(1._ireals-zzz)+ln_006*zzz) ! Chadburn et al.,2015, Dankers et al., 2011
        ELSE
          zzz = 0._ireals
          zlamdry=zlamdry_soil
        END IF

  ! Kersten number

        zsri = MIN(1.0_ireals, w_so_now(i,kso)/zdzhs(kso) / zthetas) ! degree of saturation

        IF ( t_so_now(i,kso) < t0_melt) THEN                         ! frozen
          zKe = zsri
        ELSE                                                         ! unfrozen
          zKe = 0.0_ireals
          IF ( soiltyp_subs(i) == 3 .or. soiltyp_subs(i) == 4 ) THEN ! coarse soil
            IF ( zsri >= 0.05_ireals ) THEN
              zKe = 0.7_ireals*LOG(zsri)/ln_10 + 1.0_ireals
            ENDIF
          ELSE                                                       ! fine soil (other)
            IF ( zsri >= 0.1_ireals ) THEN
              zKe = LOG(zsri)/ln_10 + 1.0_ireals
            ENDIF
          ENDIF
        ENDIF
        zKe = MAX(0.0_ireals, (MIN(1.0_ireals, zKe)) )

  ! thermal conductivity

        ! tuning factor to indirectly account for the impact of vegetation, which does not depend on soil moisture
        IF(itype_heatcond == 3 .AND. zmls(kso) < 0.075_ireals) THEN
          zxx = 12.5_ireals*(0.075_ireals-zmls(kso))*zzz
        ELSE
          zxx = 0._ireals
        ENDIF
        hzalam(i,kso) = (zKe*(zlamsat - zlamdry) + zlamdry)*(1._ireals-zxx) + zxx*0.06_ireals

      ENDDO
    ENDDO

    DO kso = 1, ke_soil
      DO i = istarts, iends
        ! mean heat conductivity
        zalam(i,kso) = hzalam(i,kso)*hzalam(i,kso+1)*(zmls(kso+1)-zmls(kso))    &
                       / ( hzalam(i,kso)*(zmls(kso+1)-zzhls(kso))                   &
                       +   hzalam(i,kso+1)*(zzhls(kso)-zmls(kso)) )
      ENDDO
    ENDDO

  ELSE

! heat conductivity based on assumption of a soil water content which is equal to the
! average between wilting point and field capacity
    DO kso = 1, ke_soil
    DO i = istarts, iends
        zwqg         = 0.5_ireals*(zfcap(i,kso) + zpwp(i,kso))
        z4wdpv       = 4._ireals*zwqg/zporv(i,kso)
        ! heat conductivity
        zalamtmp(i,kso) =              zdlam(i)                         &
                      * (0.25_ireals + 0.30_ireals*zdlam(i)           &
                      / (1._ireals+0.75_ireals*zdlam(i)))             &
                      * MIN (z4wdpv, 1.0_ireals + (z4wdpv-1.0_ireals)   &
                      *(1.0_ireals+0.35_ireals*zdlam(i))              &
                      /(1.0_ireals+1.95_ireals*zdlam(i)))
    ENDDO
    ENDDO

    DO kso = 1, ke_soil
        DO i = istarts, iends
            zalam(i,kso) = zalam(i,kso) + zalamtmp(i,kso)
            hzalam(i,kso) = zalam(i,kso)
        ENDDO
    ENDDO

  ENDIF


! Initialisations and conversion of tch to tmch

  limit_tch(:) = .false.  !  preset problem indicator
  ldebug       = .false.

  DO i = istarts, iends
!      IF (llandmask(i)) THEN     ! for land-points only
#ifdef __ICON__
        zuv        = SQRT ( u(i)**2 + v(i)**2 )
#else
        im1        = MAX( 1, i-1)
        zuv        = 0.5_ireals*SQRT ( (u(i) + u(im1))**2 &
                               +(v(i) + v(im1))**2 )
#endif
        ztvs       = t_g (i)*(1.0_ireals + rvd_m_o*qv_s(i))

        !  'potential' temperature of lowest atmospheric layer
        zplow          = p0(i) ! + pp(i)
        zth_low (i)  =  t(i) * EXP(rdocp*LOG(ps(i)/zplow))

        zdt_atm (i)  =  zth_low(i)-t_g(i)
        zdq_atm (i)  =  qv(i)-qv_s(i)

!   introduce an artificical upper boundary on transfer coefficients in cases
!   where an extreme cooling/heating of topmost soil layer may be caused due
!   to excessive sensible/latent heat fluxes (e.g. after creation of unbalanced
!   structures in the data assimilation or following massive changes in radiative
!   forcing due to infrequent radiation computations)

!   estimate current energy budget of topmost soil layer

!   heat flux between layers 1&2 based on current temperature profile
        zg1(i)= zalam(i,1)*(t_so_now(i,1)-t_so_now(i,2))/zdzms(2)

!   estimates of sensible and latent heat flux
        zrho_atm(i)=ps(i)/(r_d*ztvs)
        zshfl(i) = tch(i)*zuv*zrho_atm(i)*cp_d*zdt_atm(i)
        zlhfl(i) = tch(i)*zuv*zrho_atm(i)*lh_v*zdq_atm(i)

        IF (zshfl(i)*zlhfl(i) >= 0._ireals) THEN
          zthfl(i) = zshfl(i) + zlhfl(i)
        ELSE IF (ABS(zshfl(i)) > ABS(zlhfl(i))) THEN
          zthfl(i) = zshfl(i) + SIGN(MIN(500._ireals,ABS(zlhfl(i))),zlhfl(i))
        ELSE
          zthfl(i) = zlhfl(i) + SIGN(MIN(500._ireals,ABS(zshfl(i))),zshfl(i))
        ENDIF

        IF (ABS(zthfl(i)) <= zepsi) zthfl(i)=SIGN(zepsi,zthfl(i))

!   net radiative fluxes at surface
        zradfl(i) = sobs(i)+thbs(i)

!   unconstrained estimated energy budget of topmost soil layer
        zeb1(i) = zthfl(i)+zradfl(i)-zg1(i)

!   energy required to melt existing snow
        ze_melt(i)=w_snow_now(i)*rho_w*lh_f     ! (J/m**2)
!   heat capacity of snow layer, limited to a snow depth of 1.5 m for consistency with subsequent calculations
        zch_snow(i)=MIN(w_snow_now(i),1.5_ireals*rho_snow_now(i)/rho_w)*rho_w*chc_i   ! (J/(m**2 K))

!   constrain transfer coefficient, if energy budget  of topmost soil layer is:
!   a) negative & surface layer is unstable (i.e   upward directed turbulent heat flux)
!   b) positive & surface layer is stable   (i.e downward directed turbulent heat flux)

        IF (zeb1(i)<0.0_ireals .AND. zthfl(i)<0.0_ireals) THEN
    ! cooling of 1st soil layer&upward SHF+LHF

          ztchv_max(i) = ( zlim_dtdt*(zch_soil*zdzhs(1)+zch_snow(i))/zdt    &
                       + ABS(zg1(i)-zradfl(i)) ) / ABS(zthfl(i)) * tch(i)*zuv
        ELSEIF (zeb1(i)>0.0_ireals .AND. zthfl(i)>0.0_ireals) THEN
    ! heating of 1st soil layer & downward SHF+LHF
    !   Note: The heat capacity of snow is only relevant for the difference
    !         between the actual temperature and the melting point. The mean
    !         snow temperature is set to the average of t_snow & t_so(1)
          IF(lmulti_snow) THEN
            zdT_snow=MIN(0._ireals, &
              &          t0_melt-0.5_ireals*(t_snow_mult_now(i,1)+t_so_now(i,1)))
          ELSE
            zdT_snow=MIN(0._ireals,t0_melt-0.5_ireals*(t_snow_now(i)+t_so_now(i,1)))
          ENDIF
          ztchv_max(i) = ( (zlim_dtdt*zch_soil*zdzhs(1)+zdT_snow*zch_snow(i)+ze_melt(i))/zdt  &
                       + ABS(zg1(i)-zradfl(i)) ) / zthfl(i) * tch(i)*zuv
        ELSE
          ! unlimited transfer coefficient
          ztchv_max(i) = HUGE(1._ireals)
        ENDIF
                                                        ! required constraint as non-turbulent
                                                        ! energy budget components alone may
        ztchv_max(i) =MAX( ztchv_max(i), zepsi)         ! exceed the upper limit in the energy
                                                        ! budget of the 1st soil layer

        ! Additional limitation for better numerical stability at long time steps
        ztchv_max(i) = MIN(ztchv_max(i),(4._ireals*zlim_dtdt*(zch_soil*zdzhs(1)+zch_snow(i))/zdt &
                      +ABS(zradfl(i)))/ABS(zthfl(i))*tch(i)*zuv)

        ztchv(i)    = tch(i)*zuv  ! transfer coefficient * velocity

        IF ( inwp_turb /= iedmf ) THEN
          LIM: IF (ztchv(i) > ztchv_max(i)) THEN
            tch(i)=ztchv_max(i)/MAX(zuv,1.E-06_ireals)
!           IF (ntstep > 10) THEN          ! control print only after initial adaptation
!                                          ! to analysis state
            limit_tch(i) = .true.          ! set switch for later use
          END IF LIM
        ENDIF

        ztmch(i) = tch(i)*zuv*g*zrho_atm(i)
!     ENDIF
  ENDDO

! counter for limitation of transfer coefficients
  m_limit = COUNT( limit_tch(:) )


! In debugging mode and if transfer coefficient occured for at least one grid point
  IF (m_limit > 0 .AND. ldebug .AND. msg_level >= 19) THEN
    WRITE(*,'(1X,A,/,2(1X,A,F10.2,A,/),1X,A,F10.2,/,1X,A,F10.3,/)')                  &
           'terra1: transfer coefficient had to be constrained',                     &
           'model time step                                 :', zdt     ,' seconds', &
           'max. temperature increment allowed per time step:',zlim_dtdt,' K',       &
           'upper soil model layer thickness                :', zdzhs(1)

      DO i = istarts, iends
#ifdef __ICON__
      IF (limit_tch(i)) THEN
        zuv        = SQRT (u(i)**2 + v(i)**2 )
#else
      im1 = MAX(1,i-1)
      IF (limit_tch(i)) THEN
        zuv        = 0.5_ireals*SQRT ( (u(i) + u(im1))**2 &
                               +(v(i) + v(im1))**2 )
#endif
        yhc       ='COOLING'
        IF (zeb1(i) > 0._ireals) Yhc='HEATING'
        END IF
      END DO
  ENDIF


! Update indicator for age of snow in top of snow layer

  DO i = istarts, iends
!    IF(llandmask(i)) THEN        ! for land-points only
      IF (w_snow_now(i) <=0.0_ireals) THEN
        ! if no snow exists, reinitialize age indicator
        freshsnow(i) = 1.0_ireals
      ELSE
        IF ( nclass_gscp >= 6 ) THEN
          zsnow_rate = prs_gsp(i)+prs_con(i)+prg_gsp(i) ! [kg/m**2 s]
        ELSE
          zsnow_rate = prs_gsp(i)+prs_con(i)              ! [kg/m**2 s]
        ENDIF

        zrain_rate = prr_gsp(i)+prr_con(i)  ! [kg/m**2 s]

        ! temperature-dependent aging timescale: 2 days at freezing point, 28 days below -15 deg C
        ztau_snow = 86400._ireals*MIN(28.0_ireals,2._ireals+1.733_ireals*(t0_melt-MIN(t0_melt,t_snow_now(i))))

        ! wind-dependent snow aging: a thin snow cover tends to get broken under strong winds, which reduces the albedo
        ! an offset is added in order to ensure moderate aging for low snow depths
        zuv = MIN(300._ireals, u_10m(i)**2 + v_10m(i)**2 + 12._ireals )
        ztau_snow = MIN(ztau_snow,MAX(86400._ireals,2.e8_ireals*MAX(0.05_ireals,h_snow_gp(i))/zuv))

        ! decay rate for fresh snow including contribution by rain (full aging after 10 mm of rain)
        zdsn_old   = zdt/ztau_snow + zdt*zrain_rate*0.1_ireals

        ! linear growth rate equals 1.0 in 1 day for a temperature-dependent snow rate between
        ! 10 mmH2O (kg/m**2) per day (0.1) and 5 mmH2O (kg/m**2) per day (0.2)
        zdsn_new   = zdt*zsnow_rate*(0.1_ireals + MIN(0.1_ireals,0.02_ireals*(t0_melt-t(i))))

        ! reduce decay rate, if new snow is falling and as function of snow
        ! age itself
        zdsn_old   = (zdsn_old - zdsn_new)*freshsnow(i)
        zdsn_old   = MAX(zdsn_old,0._ireals)

        freshsnow(i) = freshsnow(i) + zdsn_new-zdsn_old

        freshsnow(i) = MIN(1._ireals,MAX(0._ireals,freshsnow(i)))

     END IF
!   END IF
  ENDDO

!------------------------------------------------------------------------------
! Section I.2: temperatures, water contents (in mH2O), surface pressure,
!------------------------------------------------------------------------------

  IF(lmulti_snow) THEN

    DO ksn = 0,ke_snow
        DO i = istarts, iends
!          IF (llandmask(i)) THEN     ! for land-points only
            IF (w_snow_now(i) > 0.0_ireals) THEN
              ! existence of snow
              t_snow_mult_now(i,ksn) = MIN (t0_melt - zepsi, t_snow_mult_now(i,ksn) )
            ELSE IF (t_snow_mult_now(i,ke_snow) >= t0_melt) THEN
              ! no snow and t_snow >= t0_melt --> t_s > t0_melt and t_snow = t_s
              t_snow_mult_now(i,ksn) = MAX (t0_melt + zepsi, t_s_now(i) )
            ELSE
              ! no snow and  t_snow < t0_melt
              ! --> t_snow = t_s
              t_snow_mult_now(i,ksn) = MIN (t0_melt - zepsi, t_s_now(i) )
            END IF
!          END IF
        ENDDO
    ENDDO

    DO i = istarts, iends
!      IF (llandmask(i)) THEN     ! for land-points only
        IF (w_snow_now(i) > 0.0_ireals) THEN
          ! existence of snow
          ! --> no water in interception store and t_snow < t0_melt
          w_snow_now(i) = w_snow_now(i) + w_i_now(i)
          wtot_snow_now(i,1) = wtot_snow_now(i,1) + w_i_now(i)
          dzh_snow_now(i,1)  = dzh_snow_now(i,1)  + w_i_now(i) &
            &                      /rho_snow_mult_now(i,1)*rho_w
          w_i_now(i) = 0.0_ireals
          h_snow_now(i) = h_snow(i)
        ELSE IF (t_snow_mult_now(i,ke_snow) >= t0_melt) THEN
          ! no snow and t_snow >= t0_melt --> t_s > t0_melt and t_snow = t_s
          t_s_now   (i) = t_snow_mult_now(i,ke_snow)
          h_snow_now(i) = 0.0_ireals
        ELSE
          ! no snow and  t_snow < t0_melt
          ! --> t_snow = t_s and no water w_i in interception store
          t_s_now   (i) = t_snow_mult_now(i,ke_snow)
          w_i_now(i) = 0.0_ireals
          h_snow_now(i) = 0.0_ireals
        END IF
!      END IF
    ENDDO
  ELSE ! no multi-layer snow

    DO i = istarts, iends
!      IF (llandmask(i)) THEN     ! for land-points only
        IF (w_snow_now(i) > 0.0_ireals) THEN
          ! existence of snow
          ! --> no water in interception store and t_snow < t0_melt
          ! GZ: this effectively suppresses rime formation because deposition rates per time
          ! step are usually less than 1e-6 m (zepsi)
      !!!    w_snow_now(i) = w_snow_now(i) + w_i_now(i)
      !!!    w_i_now(i) = 0.0_ireals
          t_snow_now(i) = MIN (t0_melt - zepsi, t_snow_now(i) )
        ELSE IF (t_snow_now(i) >= t0_melt) THEN
          ! no snow and t_snow >= t0_melt --> t_s > t0_melt and t_snow = t_s
          t_s_now   (i) = MAX (t0_melt + zepsi, t_s_now(i) )
          t_snow_now(i) = t_s_now(i)
        ELSE
          ! no snow and  t_snow < t0_melt
          ! --> t_snow = t_s and no water w_i in interception store
          t_s_now   (i) = MIN (t0_melt - zepsi, t_s_now(i) )
          t_snow_now(i) = t_s_now(i)
     !!!     w_i_now(i) = 0.0_ireals
        END IF
!      END IF
    ENDDO

  ENDIF

  !Inizializations for the next sections
  IF (lmulti_snow) THEN
    ! some preparations for ksn==0 and ksn==1
    DO i = istarts, iends
!      IF (llandmask(i)) THEN   ! for land-points only
        ztsnow_mult   (i,0) = t_snow_mult_now(i,0)
        ztsnow_mult   (i,1) = t_snow_mult_now(i,1)

        zhh_snow(i,1) =  -h_snow_now(i) + dzh_snow_now(i,1)
        zhm_snow(i,1) = (-h_snow_now(i) + zhh_snow(i,1))/2._ireals

        zdzh_snow(i,1) = dzh_snow_now(i,1)
        zextinct (i,1) = 0.13_ireals*rho_snow_mult_now(i,1)+3.4_ireals

        ! set ztsnow to ztsnow_mult(ksn=1)
        ztsnow   (i) = t_snow_mult_now(i,1)

        ! initialize zdz_snow (from Section I.3) to zdzh_snow(i,1)
        zdz_snow (i) = dzh_snow_now(i,1) ! zdzh_snow
        zwsnow   (i) = w_snow_now(i)
!      END IF
    ENDDO

    DO  ksn = 2,ke_snow
        DO i = istarts, iends
!          IF (llandmask(i)) THEN   ! for land-points only
            ztsnow_mult(i,ksn) = t_snow_mult_now(i,ksn)

            zhh_snow   (i,ksn) = zhh_snow(i,ksn-1) + dzh_snow_now(i,ksn)
            zhm_snow   (i,ksn) = (zhh_snow(i,ksn) + zhh_snow(i,ksn-1))/2._ireals

            zdzh_snow  (i,ksn) = dzh_snow_now(i,ksn)
            zextinct   (i,ksn) = 0.13_ireals*rho_snow_mult_now(i,ksn)+3.4_ireals

            ! build sum over all layers for zdzh_snow in zdz_snow
            zdz_snow   (i)     = zdz_snow(i) + zdzh_snow(i,ksn)
 !         END IF
        ENDDO
    ENDDO
  ELSE

    ! set ztsnow to t_snow
    DO i = istarts, iends
!      IF (llandmask(i)) THEN   ! for land-points only
        ztsnow   (i) = t_snow_now(i)
        zwsnow   (i) = w_snow_now(i)
        zdz_snow (i) = zwsnow(i)*rho_w/rho_snow_now(i)
!      END IF
    ENDDO
  ENDIF

  DO i = istarts, iends
!    IF (llandmask(i)) THEN   ! for land-points only

    ! ztsnow is now set according to lmulti_snow (see above);
    ! therefore we need not take care of it in this loop
    ! ztsnow   (i) = t_snow(i,nx)
    zts      (i) = t_s_now   (i)
    zts_pm   (i) = zsf_heav(zts   (i) - t0_melt)
    ztfunc   (i) = MAX(0._ireals,1._ireals - MAX(0._ireals,0.5_ireals*(zts(i)-t0_melt)))
    ztsnow_pm(i) = zsf_heav(ztsnow(i) - t0_melt)
   IF (itype_interception == 1) THEN
    zwin     (i) = w_i_now(i)
    w_p_now (i) =  0._ireals
    w_p_new (i) =  0._ireals
    w_s_now (i) =  0._ireals
    w_s_new (i) =  0._ireals
   ELSE IF (itype_interception == 2) THEN
    zwin     (i) = w_i_now   (i)
    zwpn     (i) = w_p_now   (i)
    zwisn    (i) = w_s_now   (i)
   END IF


    ! moisture and potential temperature of lowest atmospheric layer
    zplow          = p0(i) ! + pp(i)
    zqvlow         = qv(i)
    zth_low (i)  =  t(i) * EXP(rdocp*LOG(ps(i)/zplow))

    ! density*transfer coefficient*wind velocity
    zrhoch(i)    = ztmch(i)*(1._ireals/g) + zepsi

    ! saturation specific humidity for t_s and t_snow and first derivative
    z2iw        = zts_pm(i)*b2w + (1._ireals - zts_pm(i))*b2i
    z4iw        = zts_pm(i)*b4w + (1._ireals - zts_pm(i))*b4i
    z234iw      = z2iw*(b3 - z4iw)
    zqs         = zsf_qsat( zsf_psat_iw(zts(i), z2iw,z4iw), ps(i) )
    zdqs        = zqvlow - zqs
    IF (ABS(zdqs).LT.0.01_ireals*zepsi) zdqs = 0._ireals
    z2iw        = ztsnow_pm(i)*b2w + (1._ireals - ztsnow_pm(i))*b2i
    z4iw        = ztsnow_pm(i)*b4w + (1._ireals - ztsnow_pm(i))*b4i
    z234iw      = z2iw*(b3 - z4iw)
    zqsnow      = zsf_qsat(zsf_psat_iw(ztsnow(i),z2iw,z4iw), ps(i))
    zdqvtsnow(i)= zsf_dqvdt_iw(ztsnow(i), zqsnow, z4iw,z234iw)
    zdqsnow     = zqvlow - zqsnow
    IF (ABS(zdqsnow).LT.0.01_ireals*zepsi) zdqsnow = 0._ireals

    ! potential evaporation at T_snow and Ts
    zep_snow(i) = (1._ireals-ztsnow_pm(i))* tfv(i)*zrhoch(i)*zdqsnow
    zep_s   (i) =                           tfv(i)*zrhoch(i)*zdqs
!    END IF
  ENDDO



!------------------------------------------------------------------------------
! Section I.3: heat conductivity, frozen fraction, snow and
!            water covered fraction, snow density and  height,
!            volumetric heat content
!------------------------------------------------------------------------------

  DO i = istarts, iends
!      IF (llandmask(i)) THEN          ! land-points only
        ! snow and water covered fraction
!em        zrss = MAX( 0.01_ireals, MIN(1.0_ireals,zwsnow(i)/cf_snow) )
        zzz  = MAX(0.25*cf_w,0.4*cwimax_ml*tai(i))
        zrww = MAX( 0.01_ireals, 1.0_ireals -                           &
                                 EXP(MAX( -5.0_ireals, - zwin(i)/zzz) ) )
!em        zf_snow(i) = zrss*zsf_heav(zwsnow(i) - zepsi)
  IF (itype_interception == 1) THEN
        zf_wi  (i) = zrww*zsf_heav(zwin  (i) - 1.e-4_ireals*zepsi)
  ELSE IF (itype_interception == 2) THEN
        zf_wi  (i) = plcov (i) !Fraction of interception store on grid box area scales with plcov
  END IF

! BR 7/2005 prognostic  snow density
! US DMironov: for the FLake Model, h_snow has to be a prognostic variable but
!              here it is used only diagnostically. So the values are put
!              to every time level
        ! zdz_snow has been computed in the initializations above depending
        ! on lmulti_snow
        ! zdz_snow(i)=zwsnow(i)*rho_w/rho_snow(i,nx)
        h_snow_new(i) = zdz_snow(i)
        h_snow_now(i) = zdz_snow(i)
!        IF (.NOT. l2tls) h_snow(i,nold) = zdz_snow(i)

!       constrain snow depth and consider this constraint for the computation
!       of average snow density of snow layer
        zf_snow_lim(i) = MAX(0.01_ireals,0.1_ireals*freshsnow(i),zf_snow(i))
        zdz_snow (i) =  zdz_snow(i)/zf_snow_lim(i)
        zdz_snow (i) =  MAX(cdsmin,zdz_snow(i))

!       limitation of snow depth to 1.5m for snow cover heat transfer
        zdz_snow_fl(i) = MIN(1.5_iREALs, zdz_snow(i))
        IF(.NOT. lmulti_snow) &
          zrocs(i) = chc_i*zdz_snow_fl(i)*zrho_snow(i)
!BR 7/2005 End
!      END IF
  ENDDO



!------------------------------------------------------------------------------
! Section I.4: Hydrology, 1.Section
!------------------------------------------------------------------------------

  !----------------------------------------------------------------------------
  ! Section I.4.1: Evaporation from interception store and from snow cover,

  !----------------------------------------------------------------------------
  ! Evaporation and transpiration are negative, dew and rime
  ! positive quantities, since positive sign indicates a flux
  ! directed towards the earth's surface!

  DO i = istarts, iends
!      IF (llandmask(i)) THEN             ! land points only
        ! Evaporation from interception store if it contains water (wi>0) and
        ! if zep_s<0 indicates potential evaporation for temperature Ts
        ! amount of water evaporated is limited to total content of store
        zzz = (1._ireals + 0.5_ireals*ztfunc(i))/3._ireals
        zdwidt(i) = zsf_heav(-zep_s(i)) * MAX(-zrhwddt*zwin(i),              &
          zzz*zf_wi(i)*zep_s(i), -MAX(300._ireals,0.75_ireals*zradfl(i))/lh_v)
        ! Evaporation of snow, if snow exists (wsnow>0) and if zep_snow<0
        ! indicates potential evaporation for temperature t_snow
        zdwsndt(i) = zsf_heav(-zep_snow(i))  &
                       * MAX(-zrhwddt*zwsnow(i), zf_snow(i)*zep_snow(i))
        ! Formation of dew or rime, if zep_s > 0 . distinction between
        ! dew or rime is only controlled by sign of surface temperature
        ! and not effected by presence of snow !
        zrr(i)=zsf_heav(zep_s   (i))*zep_s   (i)*            zts_pm(i)
        zrs(i)=zsf_heav(zep_snow(i))*zep_snow(i)*(1.0_ireals-zts_pm(i))
!      END IF
  ENDDO

    IF (itype_interception == 2) THEN

  DO i = istarts, iends
        zwimax(i) = MAX(1.E-6_ireals,4.E-4_ireals * sai(i)) ! Security min. value 1E-6 m
        zpercmax = 2.E-3_ireals
        zfd_wi(i)=MIN(1._ireals,MAX(0.0_ireals, EXP((2._ireals/3._ireals)*LOG(zwin(i)/zwimax(i))) ))
        zf_pd(i) = MIN(1._ireals,MAX(0.0_ireals, EXP((2._ireals/3._ireals)*LOG(zwpn(i)/zpercmax)) ))

        zewi(i)=zsf_heav(-zep_s(i)) * zf_wi(i)*zfd_wi(i)*zep_s(i) ! canopy covered part
        zepd(i)=zsf_heav(-zep_s(i)) * (1._ireals - zf_wi(i))*zf_pd(i)*zep_s(i) ! bare soil part
        zept(i)=zsf_heav(-zep_s(i)) * zep_s(i) ! potential evaporation
        zesn(i)=zsf_heav(-zep_snow(i)) * zep_snow(i) ! Snow evaporation
        zdrr(i) = zsf_heav(zep_s   (i))*zep_s   (i)*     zts_pm(i)
        zrrs(i) = zsf_heav(zep_snow(i))*zep_snow(i)*(1.0_ireals-zts_pm(i))

  END DO

  END IF

  !----------------------------------------------------------------------------
  ! Section I.4.2b: Bare soil evaporation, BATS version
  !----------------------------------------------------------------------------

  IF (itype_evsl.EQ.2) THEN
    ! Calculation of bare soil evaporation after Dickinson (1984)
    ! Determination of mean water content relative to volume of voids
      DO i = istarts, iends
!        IF (llandmask(i)) THEN       ! land points only
          IF (zep_s(i) < 0.0_ireals) THEN   ! upwards directed potential
                                              ! evaporation
            zsnull(i) = zsnull(i)/(znull*zporv(i,1))
            ! Treatment of ice (m_styp=1) and rocks (m_styp=2)
            zice   = zsf_heav(1.5_ireals - REAL(m_styp(i),ireals)) ! 1 only for ice
            zevap  = zrock(i) + zice                               ! 1 for all soil types
                                                                     ! but rock and ice (=0)
            zbeta  = 0.0_ireals
            IF (m_styp(i).ge.3) THEN ! Computations not for ice and rocks
              ! auxiliary quantities
              zbf1   = 5.5_ireals - 0.8_ireals* zbedi(i)*                &
                      (1.0_ireals + 0.1_ireals*(zbedi(i) - 4.0_ireals)*  &
                       clgk0(m_styp(i)) )
              zbf2   = (zbedi(i) - 3.7_ireals + 5.0_ireals/zbedi(i))/  &
                      (5.0_ireals + zbedi(i))
              zdmax  = zbedi(i)*cfinull*zk0di(i)/crhowm
              zs1(i)  = zs1(i)/(z1*zporv(i,1))
              zd     = 1.02_ireals*zdmax*EXP( (zbedi(i)+2._ireals)*LOG(zs1(i)) ) * &
                                         EXP( zbf1*LOG(zsnull(i)/zs1(i)) )
              zck    = (1.0_ireals + 1550.0_ireals*cdmin/zdmax)*zbf2
              ! maximum sustainable moisture flux in the uppermost surface
              ! layer in kg/(s*m**2)
              zfqmax = - rho_w*zck*zd*zsnull(i)/SQRT(znull*z1)
              zevapor= MAX(zep_s(i),zfqmax)
              IF(zw_fr(i,1)+zevapor*(1.0_ireals - zf_wi(i))   &
                                     *(1.0_ireals - zf_snow(i)) &
                                     *eai(i)/sai(i)* zdtdrhw/zdzhs(1) &
                                     .LE.zadp(i,1)) zevapor = 0._ireals
              zbeta  = zevapor/MIN(zep_s(i),-zepsi)
            END IF ! Computations not for ice and rocks
            zbeta  = zbeta + (1.0_ireals - zbeta)*zice
            ! zbeta=1 (ice), zbeta=0 (rocks), zbeta unchanged for all other
            ! soil types
            ! consideration of plant or snow/water cover
         IF (itype_interception == 1) THEN
            zesoil(i) = zevap*zbeta*zep_s(i)       & ! evaporation
            !!!              *(1.0_ireals - zf_wi  (i)) & ! not water covered
                          *(1.0_ireals - zf_snow(i)) & ! not snow covered
                          * eai(i)/sai(i) ! relative source surface
                                              ! of the bare soil

         ELSE IF (itype_interception == 2) THEN
            zesoil(i) = zevap*zbeta*zep_s(i)       & ! evaporation
                          *(1.0_ireals - plcov(i))   & ! plant cover weighting
                          *(1.0_ireals - zf_snow(i)) & ! not snow covered
                          *(1.0_ireals - zf_pd(i))     ! not pond covered
         END IF ! interception
            lhfl_bs(i) = lh_v * zesoil(i)
        END IF  ! upwards directed potential evaporation
!        END IF    ! land points
      END DO
  END IF ! BATS version


  !----------------------------------------------------------------------------
  ! Section I.4.2c: Bare soil evaporation, Noilhan and Planton, 1989
  !----------------------------------------------------------------------------

  IF (itype_evsl.EQ.3) THEN
      DO i = istarts, iends
!        IF (llandmask(i)) THEN       ! land points only
          IF (zep_s(i) < 0.0_ireals) THEN   ! upwards directed potential
                                              ! evaporation
            zsnull(i) = zsnull(i)/(znull*zporv(i,1))
            ! Treatment of ice (m_styp=1) and rocks (m_styp=2)
            zice   = zsf_heav(1.5_ireals - REAL(m_styp(i),ireals)) ! 1 only for ice
            zevap  = zrock(i) + zice                  ! 1 for all soil types
                                                        ! but rock and ice (=0)
            zbeta  = 0.0_ireals
            IF (m_styp(i).ge.3) THEN ! Computations not for ice and rocks

               IF (zw_fr(i,1)> zfcap(i,1)) THEN
                  zalpha = 1.0_ireals
               ELSE
                  zalpha = 0.5_ireals * (1.0_ireals - COS ( pi *                  &
                           (zw_fr(i,1) - zadp(i,1)) / ( zfcap(i,1) - zadp(i,1)) ) )
               ENDIF
               z2iw   = ztsnow_pm(i)*b2w + (1._ireals - ztsnow_pm(i))*b2i
               z4iw   = ztsnow_pm(i)*b4w + (1._ireals - ztsnow_pm(i))*b4i
               zqs    = zsf_qsat( zsf_psat_iw(zts(i), z2iw,z4iw), ps(i) )
               zevapor= MIN(0.0_ireals,zrhoch(i)*(qv(i)-zalpha*zqs))

               zbeta  = zevapor/MIN(zep_s(i),-zepsi)
            END IF ! Computations not for ice and rocks
            zbeta  = zbeta + (1.0_ireals - zbeta)*zice
            ! zbeta=1 (ice), zbeta=0 (rocks), zbeta unchanged for all other
            ! soil types
            ! consideration of plant or snow/water cover

         IF (itype_interception == 1) THEN
            zesoil(i) = zevap*zbeta*zep_s(i)       & ! evaporation
            !!!              *(1.0_ireals - zf_wi  (i)) & ! not water covered
                          *(1.0_ireals - zf_snow(i)) & ! not snow covered
                          * eai(i)/sai(i) ! relative source surface
                                              ! of the bare soil

         ELSE IF (itype_interception == 2) THEN
            zesoil(i) = zevap*zbeta*zep_s(i)       & ! evaporation
                          *(1.0_ireals - plcov(i))   & ! plant cover weighting
                          *(1.0_ireals - zf_snow(i)) & ! not snow covered
                          *(1.0_ireals - zf_pd(i))     ! not pond covered
         END IF ! interception

            lhfl_bs(i) = lh_v * zesoil(i)
          END IF  ! upwards directed potential evaporation
!        END IF    ! land points
      END DO
  END IF ! NP89

  !----------------------------------------------------------------------------
  ! Section I.4.2d: Bare soil evaporation, resistance version
  !----------------------------------------------------------------------------

  IF (itype_evsl.EQ.4) THEN   ! Resistance version
    ! Calculation of bare soil evaporation using a resistance formulation.
    ! For a review see Schulz et al. (1998) 
    DO i = istarts, iends
      IF (zep_s(i) < 0.0_ireals) THEN   ! upwards directed potential evaporation
        ! Treatment of ice (m_styp=1) and rocks (m_styp=2)
        zice  = zsf_heav(1.5_ireals - REAL(m_styp(i),ireals)) ! 1 only for ice
        zevap = zrock(i) + zice                               ! 1 for all, but rock
        zbeta = 0.0_ireals

        IF (m_styp(i).GE.3) THEN ! Computations not for ice and rocks
           zalpha = MAX( 0.0_ireals, MIN( 1.0_ireals,                     &
                    (zw_fr(i,1) - zadp(i,1)) / (zfcap(i,1) - zadp(i,1)) ) )
           zalpha = 50.0_ireals / (zalpha + zepsi)
           zbeta  = 1.0_ireals                                &
                  / (1.0_ireals + zrhoch(i)*zalpha/zrho_atm(i))
        END IF ! Computations not for ice and rocks

        ! zbeta=1 (ice), zbeta=0 (rocks), zbeta unchanged for all other soil types
        zbeta  = zbeta + (1.0_ireals - zbeta)*zice

        ! Consideration of plant or snow/water cover
        IF (itype_interception == 1) THEN         ! Interception
           zesoil(i) = zevap*zbeta*zep_s(i)     & ! evaporation
           !!!            *(1.0_ireals - zf_wi  (i)) & ! not water covered
                     *(1.0_ireals - zf_snow(i)) & ! not snow covered
                     * eai(i)/sai(i)              ! relative source surface
                                                  ! of the bare soil
        ELSE IF (itype_interception == 2) THEN
           zesoil(i) = zevap*zbeta*zep_s(i)     & ! evaporation
                     *(1.0_ireals - plcov(i))   & ! plant cover weighting
                     *(1.0_ireals - zf_snow(i)) & ! not snow covered
                     *(1.0_ireals - zf_pd(i))     ! not pond covered
        END IF ! Interception

        lhfl_bs(i) = lh_v * zesoil(i)
      END IF  ! upwards directed potential evaporation
    END DO
  END IF      ! Resistance version


  !----------------------------------------------------------------------------
  ! Section I.4.3b: transpiration by plants, BATS version
  !----------------------------------------------------------------------------

  IF (itype_trvg.EQ.2) THEN   ! BATS version
    ! This version is based on Dickinson's (1984) BATS scheme, simplified by
    ! neglecting the water and energy transports between the soil and the plant
    ! canopy. This leads to a Monteith combination formula for the computation
    ! of plant transpiration.


  ! Root distribution

    IF (itype_root == 2) THEN
        DO i = istarts, iends
          zrootdz_int (i)= 0.0_ireals   !initialize the root density profile integral
          zwrootdz_int(i)= 0.0_ireals   !initialize the root water content integral
        END DO

      DO kso   = 1,ke_soil
!CDIR NODEP,VOVERTAKE,VOB
         DO ic=1,icount_soil
            i=soil_list(ic)
!!$          DO i = istarts, iends
!!$!            IF (llandmask(i)) THEN ! land points only,
!!$              IF (m_styp(i).ge.3) THEN ! neither ice or rocks
                IF (zep_s(i) < 0.0_ireals) THEN  ! upwards directed potential evaporation
                  ! consider the effect of root depth & root density
                  zrootfr = EXP (-zroota(i)*zmls(kso)) ! root density
                  zrootdz = zrootfr*MIN(zdzhs(kso),MAX(0.0_ireals, zbwt(i)-(zmls(kso) &
                    &       -0.5_ireals*zdzhs(kso)) ) )
                  zrootdz_int(i)=zrootdz_int(i) + zrootdz
                  zwrootdz(i,kso)=zrootdz*(zw_fr(i,kso)-w_so_ice_now(i,kso)/zdzhs(kso))
                  zwrootdz_int(i)=zwrootdz_int(i) + zwrootdz(i,kso)
                END IF  ! negative potential evaporation only
!!$              END IF  ! neither ice or rocks
!            END IF    ! land-points only
          END DO
      END DO


       ! Compute root zone integrated average of liquid water content

          DO i = istarts, iends
             zwrootdz_int(i)=zwrootdz_int(i)/MAX(zrootdz_int(i),zepsi)
          END DO

    ELSE

      DO kso   = 1,ke_soil
!CDIR NODEP,VOVERTAKE,VOB
         DO ic=1,icount_soil
            i=soil_list(ic)
!!$          DO i = istarts, iends
!!$!            IF (llandmask(i)) THEN ! land points only,
!!$              IF (m_styp(i).ge.3) THEN ! neither ice or rocks
                IF (zep_s(i) < 0.0_ireals) THEN  ! upwards directed potential
                                                   ! evaporation
                  zropart  = MIN ( zdzhs(kso), MAX(0.0_ireals,                     &
                             zbwt(i) - (zmls(kso) - 0.5_ireals*zdzhs(kso))))
                  zropartw(i,kso) = zropart*(zw_fr(i,kso)-w_so_ice_now(i,kso)/zdzhs(kso))
                  zwroot(i) = zwroot(i) + zropartw(i,kso)/MAX(zepsi,zbwt(i))
                END IF  ! negative potential evaporation only
!!$              END IF  ! neither ice or rocks
!            END IF    ! land-points only
          END DO
      END DO

    ENDIF


    ! Determination of the transfer functions CA, CF, and CV
!CDIR NODEP,VOVERTAKE,VOB
         DO ic=1,icount_soil
            i=soil_list(ic)
!!$      DO i = istarts, iends
!!$!        IF (llandmask(i)) THEN ! land points only,
!!$          IF (m_styp(i).ge.3) THEN ! neither ice or rocks
            IF (zep_s(i) < 0.0_ireals) THEN  ! upwards directed potential evaporation
              zuv        = SQRT (u_10m(i) **2 + v_10m(i)**2 )
              zcatm      = tch(i)*zuv           ! Function CA

              IF(icant == 1) THEN !additional laminar canopy resistance in case of Louis-transfer-scheme
                zustar     = zuv*SQRT(tcm(i))
                zrla       = 1.0_ireals/MAX(cdash*SQRT(zustar),zepsi)
              ELSE !in case of Raschendorfer-transfer-scheme a laminar canopy resistance is already considered
                zrla       = 0._ireals
              ENDIF


              ! to compute CV, first the stomatal resistance has to be determined
              ! this requires the determination of the F-functions:
              ! Radiation function
              zpar       = pabs(i)  !  PAR
              zf_rad(i)= MAX(0.0_ireals,MIN(1.0_ireals,zpar/cparcrit))
              ztlpmwp(i) = (zfcap(i,1) - zpwp(i,1))*(0.81_ireals +       &
                 0.121_ireals*ATAN(-86400._ireals*zep_s(i) - 4.75_ireals))

              ! Soil water function
              IF (itype_root == 2) THEN
                zf_wat     = MAX(0.0_ireals,MIN(1.0_ireals,(zwrootdz_int(i) -  &
                                                zpwp(i,1))/ztlpmwp(i)))
              ELSE
                zf_wat     = MAX(0.0_ireals,MIN(1.0_ireals,(zwroot(i) -  &
                                                zpwp(i,1))/ztlpmwp(i)))
              ENDIF

              ! Temperature function
!              IF (ntstep .EQ. 0 .AND. icant .NE. 2) THEN
!                t_2m(i)=t(i)
!              ENDIF
!             zf_tem     = MAX(0.0_ireals,MIN(1.0_ireals,4.0_ireals*     &
!                          (t_2m(i)-t0_melt)*(ctend-t_2m(i))/(ctend-t0_melt)**2))
              ! T at lowest model level used (approximation of leaf height)
              zf_tem     = MAX(0.0_ireals,MIN(1.0_ireals,4.0_ireals*     &
                           (t(i)-t0_melt)*(ctend-t(i))/(ctend-t0_melt)**2))

              ! Saturation deficit function (function not used, but computations
              ! necessary for determination of  slope of the saturation curve)
              z2iw       = zts_pm(i)*b2w + (1._ireals - zts_pm(i))*b2i
              z4iw       = zts_pm(i)*b4w + (1._ireals - zts_pm(i))*b4i
              zepsat     = zsf_psat_iw(t(i),z2iw,z4iw)
              zepke      = qv(i)*ps(i)/                   &
                                        (rdv + o_m_rdv*qv(i))
              zf_sat     = MAX(0.0_ireals,MIN(1.0_ireals,1.0_ireals -  &
                                              (zepsat - zepke)/csatdef))
              ! zf_sat paralysed:
              zf_sat     = 1.0_ireals

              IF (lstomata) THEN
                zedrstom   = 1.0_ireals/crsmax + (1.0_ireals/MAX(40._ireals,rsmin2d(i)) -     &
                             1.0_ireals/crsmax)*zf_rad(i)*zf_wat*zf_tem*zf_sat
              ELSE
                zedrstom   = 1.0_ireals/crsmax + (1.0_ireals/crsmin -     &
                             1.0_ireals/crsmax)*zf_rad(i)*zf_wat*zf_tem*zf_sat
              END IF

              zrstom     = 1.0_ireals/zedrstom              ! stomatal resistance
              rstom(i) = zrstom
              zrveg      = zrla + zrstom
              ! Transpiration rate of dry leaves:
              ztraleav(i)=zep_s(i)*tai(i)/(sai(i)+zrveg*zcatm)
            END IF  ! upwards directed potential evaporation only
!!$          END IF    ! m_styp > 2
!        END IF      ! land points
      END DO

    ! Consideration of water and snow coverage, distribution to the different
    ! soil layers

IF (itype_interception == 1) THEN
    DO     kso       = 1,ke_soil
!CDIR NODEP,VOVERTAKE,VOB
         DO ic=1,icount_soil
            i=soil_list(ic)
!!$        DO i         = istarts, iends
!!$!          IF (llandmask(i)) THEN ! land points only,
!!$            IF (m_styp(i).ge.3) THEN ! neither ice or rocks
              IF (zep_s(i) < 0.0_ireals) THEN    ! upwards potential evaporation
                ztrabpf  = ztraleav(i)*                   & ! plant covered part
                !!!           (1.0_ireals - zf_wi(i))*       & ! not water covered
                           (1.0_ireals - zf_snow(i))        ! not snow covered

                ! for root distribution
                IF (itype_root == 2) THEN
                  ztrfr    = zwrootdz(i,kso)/(zrootdz_int(i)*zwrootdz_int(i))
                  ztrang(i,kso) = ztrabpf*ztrfr
                ELSE
                  zrootfc = zropartw(i,kso)/(zwroot(i) + zepsi)
                  ztrang(i,kso) = ztrabpf*zrootfc/MAX(zepsi,zbwt(i))
                ENDIF

                ! Limit evaporation such that the soil water content does not fall beyond the wilting point
                IF(zw_fr(i,kso)+ztrang(i,kso)*zdtdrhw/zdzhs(kso) < zpwp(i,kso)) &
                  ztrang(i,kso) = MIN(0._ireals,(zpwp(i,kso)-zw_fr(i,kso))*zdzhs(kso)/zdtdrhw)

                lhfl_pl(i,kso)= lh_v * ztrang(i,kso)
                ztrangs(i)    = ztrangs(i) + ztrang(i,kso)
              END IF  ! upwards directed potential evaporation only
!!$            END IF    ! m_styp > 2
!          END IF      ! land points
        END DO
    END DO          ! loop over soil layers
ELSE          IF (itype_interception == 2) THEN
    DO     kso       = 1,ke_soil
!CDIR NODEP,VOVERTAKE,VOB
         DO ic=1,icount_soil
            i=soil_list(ic)
!!$        DO i         = istarts, iends
!!$!          IF (llandmask(i)) THEN ! land points only,
!!$            IF (m_styp(i).ge.3) THEN ! neither ice or rocks
              IF (zep_s(i) < 0.0_ireals) THEN    ! upwards potential evaporation
                ztrabpf  = ztraleav(i)*                   & ! plant covered part
                           (1.0_ireals - zfd_wi(i))*       & ! not water covered
                           (1.0_ireals - zf_snow(i))        ! not snow covered

                ! for root distribution
                IF (itype_root == 2) THEN
                  ztrfr    = zwrootdz(i,kso)/(zrootdz_int(i)*zwrootdz_int(i))
                  ztrang(i,kso) = ztrabpf*ztrfr
                ELSE
                  zrootfc = zropartw(i,kso)/(zwroot(i) + zepsi)
                  ztrang(i,kso) = ztrabpf*zrootfc/MAX(zepsi,zbwt(i))
                  IF(zw_fr(i,kso)+ztrang(i,kso)*zdtdrhw/zdzhs(kso) &
                                    .LT.zpwp(i,kso)) ztrang(i,kso) = 0._ireals
                ENDIF
                lhfl_pl(i,kso)= lh_v * ztrang(i,kso)
                ztrangs(i)    = ztrangs(i) + ztrang(i,kso)
              END IF  ! upwards directed potential evaporation only
!!$            END IF    ! m_styp > 2
!          END IF      ! land points
        END DO
    END DO          ! loop over soil layers
 END IF
  END IF ! BATS version


  !----------------------------------------------------------------------------
  ! Section I.4.4: total evapotranspiration and
  !              associated ficticious soil humidity qv_s
  !----------------------------------------------------------------------------

  ! Ensure that the sum of the evaporation terms does not exceed the potential evaporation
  DO i = istarts, iends
    ze_sum = zdwsndt(i) + zdwidt(i) + zesoil(i) + ztrangs(i)
    IF (zep_s(i) < 0._ireals .AND. ze_sum < zep_s(i)) THEN
      zzz = zep_s(i)/ze_sum
      zdwsndt(i) = zdwsndt(i)*zzz
      zdwidt(i)  = zdwidt(i) *zzz
      zesoil(i)  = zesoil(i) *zzz
      ztrangs(i) = ztrangs(i)*zzz
    ENDIF
  ENDDO

  IF (itype_interception == 1) THEN
    DO i = istarts, iends

        ze_sum = zdwsndt(i  )  & ! evaporation of snow
               + zdwidt (i  )  & ! evaporation from interception store
               + zesoil (i  )  & ! evaporation from bare soil
               + ztrangs(i  )  & ! transpiration from all soil layers
               + zrr    (i  )  & ! formation of dew
               + zrs    (i  )    ! formation of rime
        qv_s(i) = qv (i) - ze_sum /(zrhoch(i) + zepsi)
!JH     qv_s(i,nnew) = qv_s(i,nx)


    END DO
  ELSE          IF (itype_interception == 2) THEN
    DO i = istarts, iends
           ze_sum = zesn   (i) &
                  + zepd   (i) &
                  + zewi   (i) &
                  + zesoil (i) &
                  + ztrangs(i) &
                  + zdrr   (i) &
                  + zrrs   (i)
        qv_s(i) = qv (i) - ze_sum /(zrhoch(i) + zepsi)
!JH     qv_s(i,nnew) = qv_s(i,nx)

    END DO
  END IF


!------------------------------------------------------------------------------
! End of former module procedure terra1_multlay
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
! Former SUBROUTINE terra2_multlay (yerror, ierror)
!------------------------------------------------------------------------------

!   In the prognostic part II the equation of heat conduction and water
!   transport is solved for a multi-layer soil using the same vertical grid
!   Freezing/melting of soil water/ice is accounted for (optionally). A
!   simple one-layer snow model provides the snow surface temperature and
!   the snow water equivalent.


!------------------------------------------------------------------------------
! Section II.1: Initializations
!------------------------------------------------------------------------------

  ! Computation of derived constants
  z1d2dt = 1._ireals/zdt      ! 1./2*timestep

  ! Number of soil layers contributing to surface run-off
  msr_off  = 0

  zdtdrhw = zdt/rho_w   ! timestep/density of liquid water

  ! time constant for infiltration of water from interception store
  ! must not be less than 2*time step
  ctau_i  = MAX( ctau_i, zdt )

  ! Utility variable to avoid IF-constructs
  DO kso     = 1,ke_soil
    znlgw1f(kso)  = 0.0_ireals
  END DO
  znlgw1f(1)         = 1.0_ireals

! Initialisations
  IF (lmulti_snow) THEN
    DO ksn = 0, ke_snow
        DO i = istarts, iends
!          IF (llandmask(i)) THEN                 ! land-points
            zdtsnowdt_mult(i,ksn)  = 0.0_ireals
            zdtsdt(i)     = 0.0_ireals
!          END IF
        END DO
    END DO
  ELSE
      DO i = istarts, iends
!        IF (llandmask(i)) THEN                 ! land-points
          zdtsnowdt(i)  = 0.0_ireals
          zdtsdt(i)     = 0.0_ireals
!        END IF
      END DO
  END IF



!------------------------------------------------------------------------------
! Section II.2: Prepare basic surface properties and create some local
!               arrays of surface related quantities (for land-points only)
!               Initialise some fields
!------------------------------------------------------------------------------

  DO i = istarts, iends
!      IF (llandmask(i)) THEN                 ! land-points
        mstyp        = soiltyp_subs(i)
        m_styp(i)  = mstyp
        zsandf(i)  = csandf(mstyp)
        zclayf(i)  = cclayf(mstyp)
        zgsb   (i) = 0.0_ireals
        ztrangs(i) = 0.0_ireals
!      END IF
  END DO


  DO   kso = 1,ke_soil+1
      DO i = istarts, iends
!        IF (llandmask(i)) THEN     ! land-points only
          ziw_fr(i,kso) = w_so_ice_now(i,kso)/zdzhs(kso)   ! ice frac.
          zlw_fr(i,kso) = zw_fr(i,kso) - ziw_fr(i,kso)  ! liquid water frac.
          zroc(i,kso)   = zrocg(i,kso) + rho_w*zlw_fr(i,kso)*chc_w +          &
                                         rho_w*ziw_fr(i,kso)*chc_i
          IF (kso<=ke_soil) THEN
            ztrangs(i) = ztrangs(i) + ztrang(i,kso)
            zdwgdt(i,kso) = 0.0_ireals
          END IF
          zflmg(i,kso)  = 0.0_ireals
          zrunoff_grav(i,kso)  = 0.0_ireals
!        END IF
      END DO
  END DO      !soil layers



!------------------------------------------------------------------------------
! Section II.3: Estimate thermal surface fluxes
!------------------------------------------------------------------------------
  IF (itype_interception == 1) THEN
  DO i = istarts, iends
!      IF(llandmask(i))THEN     ! land-points only

        ! store forcing terms due to evapotranspiration, formation of dew
        ! and rime for later use
        zverbo(i) = zdwidt(i) + zesoil(i) + ztrangs(i) +              &
                        (1._ireals-zf_snow(i))*(zrr(i) + zrs(i))

        zversn(i) = zdwsndt(i) + zrs(i)                                 &
                                  + zsf_heav (zwsnow(i) - zepsi) * zrr(i)

        ! add grid scale and convective precipitation (and graupel, if present)
        ! to dew and rime
        zrr(i) = zrr(i) + prr_con(i) + prr_gsp(i)
        zrime  = zrs(i)
        zrs(i) = zrs(i) + prs_con(i) + prs_gsp(i)
        IF ( nclass_gscp >= 6 ) zrs(i) = zrs(i) + prg_gsp(i)

        ! Decide whether riming is added to interception store or to snow cover
        IF (zrs(i) >= 1.05_ireals*zrime .OR. zf_snow(i) >= 0.9_ireals) zrime = 0._ireals
        zrs(i) = zrs(i) - zrime

        ! infiltration and surface run-off

        ! ice free fraction of first soil layer scaled by pore volume
        ! is used as reduction factor for infiltration rate
        zfr_ice_free     = 1._ireals-ziw_fr(i,1)/zporv(i,1)

        ! subtract evaporation from interception store to avoid negative
        ! values due to sum of evaporation+infiltration
        zwinstr(i) = zwin(i) + zdwidt(i)*zdtdrhw
        zwinstr(i) = MAX(0.0_ireals,zwinstr(i))

        ! maximum infiltration rate of the soil (rock/ice/water-exclusion
   !     zinfmx(i) = zrock(i)*zfr_ice_free*csvoro &
   !              *( cik1*MAX(0.5_ireals,plcov(i))*MAX(0.0_ireals,           &
   !              zporv(i,1)-zw_fr(i,1))/zporv(i,1) + zik2(i) )

        zinfmx(i) = zrock(i)*zfr_ice_free*csvoro*zkw(i,1)*rho_w


        ! to avoid pore volume water excess of the uppermost layer by
        ! infiltration
        zinfmx(i) = MIN(zinfmx(i), (zporv(i,1) - zw_fr(i,1))*zdzhs(1)*zrhwddt)

        ! to avoid infiltration at snow covered parts of soil surface
        zinfmx(i) = zinfmx(i)*(1._ireals - zf_snow(i))

        zwimax(i) = cwimax_ml*(1._ireals+ztfunc(i))*MAX(ztfunc(i), zepsi, tai(i))
        zalf   = SQRT(MAX(0.0_ireals,1.0_ireals - zwinstr(i)/zwimax(i)))

        ! water supply from interception store (if Ts above freezing)
        zuv    = SQRT ( u_10m(i)**2 + v_10m(i)**2 )
        zzz    = MAX(0.1_ireals, 0.4_ireals - 0.05_ireals*zuv)
        zinf   = MAX(0._ireals,zwinstr(i)-zzz*zwimax(i))*rho_w/ctau_i*                       &
                 (1._ireals+0.75_ireals*MAX(0._ireals,zuv-1._ireals))*(1._ireals-ztfunc(i))**2

        ! possible contribution of rain to infiltration
    !    IF (zrr(i)-zepsi > 0.0_ireals) THEN
    !      zalf = MAX( zalf,                                                   &
    !             (zrhwddt*MAX(0.0_ireals, zwimax(i)-zwinstr(i)) + zinf)/zrr(i) )
    !      zalf = MAX( 0.01_ireals, MIN(1.0_ireals, zalf) )

          ! if rain falls onto snow, all rain is considered for infiltration
          ! as no liquid water store is considered in the snowpack
          IF (zwsnow(i) > 0.0_ireals) zalf = 0.0_ireals
    !    END IF
        ! Increase infiltration and reduce surface runoff (bugfix)
    !    zalf = 0.0_ireals ! this deactivates filling the interception store!!
        ! rain freezes on the snow surface
        IF (lmulti_snow .AND. zwsnow(i) > 0.0_ireals) zalf = 1.0_ireals

        ! interception store; convective precip is taken into account with a fractional area passed from the convection scheme
        zdwidtt  = zalf*(zrr(i)+(conv_frac(i)-1._ireals)*prr_con(i)) + zrime + zdwidt(i)-zinf
        zwinstr(i)  = zwin(i) + zdwidtt*zdtdrhw
        zwinstr(i)  = MAX(0.0_ireals, zwinstr(i)) !avoid negative values (security)
        zdwieps  = 0.0_ireals
        IF (zwinstr(i) > 0.0_ireals .AND. zwinstr(i) < 1.e-4_ireals*zepsi) THEN
          zdwieps    = zwinstr(i)*zrhwddt
          runoff_s(i)= runoff_s(i) + zdwieps*zroffdt
          zdwidtt    = zdwidtt - zdwieps
          zwinstr(i)    = 0.0_ireals
        END IF
        ! add excess over zwimax(i) to infiltration
        zro_wi       = zrhwddt*MAX( 0.0_ireals, zwinstr(i)-zwimax(i) )
        zdwidtt      = zdwidtt - zro_wi
        zdwidt(i)    = zdwidtt
        zinf         = zinf + zro_wi
        IF (zts(i) <= t0_melt) THEN
          ! add excess rime to snow
          zrs(i) = zrs(i) + zinf
          zinf = 0._ireals
        ENDIF

        ! add rain contribution to water supply for infiltration
        zvers(i) = zinf + (1._ireals - zalf)*zrr(i) + (1._ireals-conv_frac(i))*zalf*prr_con(i)

        ! compute surface runoff; accounting for the fractional area of convective precip proved to be disadvantageous
        ! here because of excessive drying of the soil in longer-term runs
        zro_inf  = MAX(0._ireals,zvers(i)-zinfmx(i))

        ! final infiltration rate
        zinfil(i) = zvers(i) - zro_inf

        ! surface run-off (residual of potential minus actual infiltration)
        runoff_s(i) = runoff_s(i) + zro_inf*zroffdt

        ! change of snow water and interception water store
        ! (negligible residuals are added to the run-off)

        ! snow store
        zdwsndtt = zrs(i) + zdwsndt(i)
        zwsnstr  = zwsnow(i) + zdwsndtt*zdtdrhw
        zwsnstr  = MAX(0.0_ireals, zwsnstr) ! avoid negative values (security)
        zdwseps  = 0.0_ireals
        IF (zwsnstr > 0.0_ireals .AND. zwsnstr < zepsi) THEN ! shift marginal snow amounts to interception storage
!         IF (ztsnow_pm(i) > 0.0_ireals) THEN
            zwinstr(i) = zwinstr(i) + zwsnstr
            zdwseps    = zwsnstr*zrhwddt
!            runoff_s(i) = runoff_s(i) + zdwseps*zroffdt ! previous implementation
            zdwsndtt   = zdwsndtt - zdwseps
            zdwidt(i)  = zdwidt(i) + zdwseps
!         END IF
        END IF
        zdwsndt(i) = zdwsndtt

!      END IF            ! land-points only
  END DO
ELSE   IF (itype_interception == 2) THEN
  DO i = istarts, iends
!       Compute forcing terms after (???) interception, infiltration calculation --> below
!       store forcing terms due to evapotranspiration, formation of dew
!       and rime for later use
        zverbo(i) = zewi(i) + zesoil(i) + zepd(i) + ztrangs(i) +              &
                        (1._ireals-zf_snow(i))*(zdrr(i) + zrrs(i))

        zversn(i) = zdwsndt(i) + zrrs(i)                                 &
                                  + zsf_heav (zwsnow(i) - zepsi) * zdrr(i)

!       ice free fraction of first soil layer scaled by pore volume
!       is used as reduction factor for infiltration rate (previously zts_pm switch bobo)
        zfr_ice_free     = 1._ireals-ziw_fr(i,1)/zporv(i,1)

!       add grid scale and convective precipitation (and graupel, if present) to dew and rime
        zrr(i) = zdrr(i) + prr_con(i) + prr_gsp(i)

        IF ( nclass_gscp >= 6 ) THEN
          zrs(i) = zrrs(i) + prs_con(i) + prs_gsp(i) + prg_gsp(i)
        ELSE
          zrs(i) = zrrs(i) + prs_con(i) + prs_gsp(i)
        ENDIF

!       Preliminary interception store budget
!       According to Wang et al. (Evaluation of canopy interception schemes in land surface models, J. of Hydrology,2007)
!       the water balance equation for the interception store is   zdwidtt  = Ic + Ew - Dr, where Ic is the canopy
!       interception rate, Dr the canopy drip rate, and Ew is the evaporation rate from wet foliage.
!       This equation is independent of the approach used to model the canopy hydrological processes.

!       Comparable to the community land model (CLM), the precipitation intercepted by vegetation canopy Ic is
!       considered as an exponential function of canopy density using a canopy cover fraction
!       zf_wi  (j1,j2) = 1._ireals - exp(-0.5_ireals * plai(j1,j2)).  Therefore, similar to plai,
!       the  canopy cover fraction shows an annual cycle.
!       The canopy capacity zwimax(i) is considered linearly related to leaf area index
!       (van Dijk and Bruijnzeel, 2001; Wang et al., 2007)
!       The canopy dripping occurs only when canopy water storage winstr exceeds the water holding capacity zwimax(i).
!       This interception runoff is part of the soil infiltration.
!       The melting of snow on the leafs, if T > T_melt is also considered, but from a frozen interception store only
!       evaporation is allowed.




!       Interception of rain water
        IF (ztsnow(i).gt.t0_melt) THEN
!         snow fall on leafs/soil with T>T_melt, snow water content increases
!         interception store water content
          zfcorr_wi=1._ireals ! Start value

          zdwidtt  = zf_wi(i)*zrr(i) + zf_wi(i)*zrs(i) + zewi(i)
          zwinstr(i)  = zwin (i) + zdtdrhw * zdwidtt
          zewi(i) = zewi(i) - MIN(0._ireals,zwinstr(i) * zrhwddt) ! Partial evaporation


!         Final calculation of the interception store budget
          zdwidtt  = zf_wi(i)*zrr(i) + zf_wi(i)*zrs(i) + zewi(i)
          zwinstr(i)  = zwin (i) + zdtdrhw * zdwidtt

        ELSE   IF (ztsnow(i).le.t0_melt) THEN!  T =< t0_melt

          zdwidtt = zewi(i) ! only evaporation from frozen interception store allowed
          zwinstr(i)  =  zwin (i) + zdtdrhw * zdwidtt
          zewi(i) = zewi(i) - MIN(0._ireals,zwinstr(i) * zrhwddt) ! Partial evaporation


!         Final calculation of the interception store budget
          zdwidtt  =  zewi(i) ! only evaporation from frozen interception store allowed
          zwinstr(i)  = zwin (i) + zdtdrhw * zdwidtt

        END IF ! Temperature

        zro_wi        = MAX( 0.0_ireals, zwinstr(i)-zwimax(i) )   ! excess over zwimax -> runoff and infiltration
        zwinstr(i)  = MAX(0.0_ireals, zwinstr(i))                 ! avoid negative values (security)
        zwinstr(i)  = MIN(zwinstr(i),zwimax(i))                 ! correction of interc. store due to runoff
        zdwidt(i)   = zdwidtt - zro_wi

        zuv        = SQRT ( u(i)**2 + v(i)**2 )

!       SNOW Interception Model (Roesch et al., 2001)
! Need forest_fraction
!        zdwisndtt=(for_e(i)+for_d(i))*zrs(i) + (for_e(i)+for_d(i))*zesn(i) - zwisn(i)* &
!        ((t(i,j,ke,nx)-270.15)/1.87E5_ireals + SQRT(u(i,j,ke,nx)*u(i,j,ke,nx)+v(i,j,ke,nx)*v(i,j,ke,nx))/1.56E5_ireals)
! Just for testing, using plcov
        zdwisndtt=plcov(i)*zrs(i) + plcov(i)*zesn(i) - zwisn(i)* &
        ((t(i)-270.15)/1.87E5_ireals + zuv/1.56E5_ireals)
        zwisnstr(i)  = zwisn (i) + zdtdrhw * zdwisndtt
        zwisnstr(i)  = MAX(0.0_ireals, zwisnstr(i)) !avoid negative values (security)



!       Calculation of infiltration
!       add rain and interception excess to water supply for infiltration
        zvers(i) =  zrhwddt*zro_wi + (1._ireals - zf_wi(i))*zrr(i) ! Test  Excess over zwimax will infiltrated

!       maximum infiltration rate of the soil (rock/ice/water-exclusion
  !      zinfmx(i) = zrock(i)*zfr_ice_free *csvoro &
  !       *( cik1*MAX(0.5_ireals,plcov(i))*MAX(0.0_ireals,           &
  !          zporv(i,1)-zw_fr(i,1))/zporv(i,1) + zik2(i) )

        zinfmx(i) = zrock(i)*zfr_ice_free*csvoro*zkw(i,1)*rho_w

!       to avoid pore volume water excess of the uppermost layer by infiltration
        zinfmx(i) = MIN(zinfmx(i), (zporv(i,1) - zw_fr(i,1))*zdzhs(1)*zrhwddt)

!       final infiltration rate limited by maximum value
        zinfil(i) = MIN(zinfmx(i),zvers(i)+ zwpn(i)*zrhwddt ) !GME

!       surface run-off (residual of potential minus actual infiltration)
        zro_inf     = MAX(0._ireals, zvers(i) - zinfil(i))

!       Pond interception
        zdwpdtt  = zvers(i) - zinfil(i)  + zepd(i)
        zwpnstr(i)  = zwpn   (i) + zdtdrhw * zdwpdtt
        IF (zvers(i) - zinfil(i).ge.0._ireals) then
          zepd(i) = zepd(i) - MIN(0._ireals,zwpnstr(i) * zrhwddt) ! Partial evaporation
        ELSE
          zepd(i) = zepd(i) - MIN(0._ireals,zwpn(i) * zrhwddt) ! Partial evaporation
        END IF

!       Final calculation
        zdwpdtt  = zvers(i) - zinfil(i)  + zepd(i)
        zwpnstr(i)  = MAX(0._ireals,zwpn (i) + zdtdrhw * zdwpdtt )
        zro_inf = MAX(0._ireals, (zwpnstr(i)  - zpercmax)*zrhwddt )   ! Pond overflow?
        zwpnstr(i)  = MIN(zwpnstr(i),zpercmax) ! Correction of interc. store due to runoff

!       Calculation of surface runoff
!       change of snow water and interception water store
!       (negligible residuals are added to the run-off)

!       snow store
        zdwsndtt = zrs(i) + zdwsndt(i)
        zwsnstr  = w_s_now(i) + zdwsndtt*zdtdrhw
        zwsnstr  = MAX(0.0_ireals, zwsnstr) ! avoid negative values (security)
        zdwseps  = 0.0_ireals
        IF (zwsnstr > 0.0_ireals .AND. zwsnstr < zepsi) THEN
          zdwseps    = zwsnstr*zrhwddt
          runoff_s(i) = runoff_s(i) + zdwseps*zroffdt
          zdwsndtt   = zdwsndtt - zdwseps
        END IF
        zdwsndt(i) = zdwsndtt
        runoff_s(i)= runoff_s(i) + zro_wi*zroffdt

     END DO
  END IF


!------------------------------------------------------------------------------
! Section II.4: Soil water transport and runoff from soil layers
!------------------------------------------------------------------------------

! uppermost layer, kso = 1
!CDIR NODEP,VOVERTAKE,VOB

      DO ic = 1, icount_soil
        i=soil_list(ic)
      ! sedimentation and capillary transport in soil
      ! Note: The fractional liquid water content (concentration)  of each layer
      !       is normalized by the ice free fraction of each layer in order to
      !       obtain a representative concentration of liquid water in the
      !       'active' part of each soil layer
      !       Hydraulic diffusivity and conductivity coefficients are multiplied
      !       by a reduction factor depending on the maximum ice fraction of the
      !       adjacent layers in order to avoid the transport of liquid water
      !       in to the frozen part of the adjacent layer
          zice_fr_kso   = ziw_fr(i,1)
          zice_fr_ksop1 = ziw_fr(i,2)
          zlw_fr_kso    = zlw_fr(i,1)
          zlw_fr_ksop1  = zlw_fr(i,2)

          ! compute reduction factor for transport coefficients
          zfr_ice  = max (zice_fr_kso,zice_fr_ksop1)
          zredp05  = 1._ireals-zfr_ice/MAX(zlw_fr_kso+zice_fr_kso,zlw_fr_ksop1+zice_fr_ksop1)

          ! interpolated scaled liquid water fraction at layer interface
          zlw_fr_ksop05  = 0.5_ireals*(zdzhs(2)*zlw_fr_kso+zdzhs(1)*zlw_fr_ksop1) &
                                               /zdzms(2)
!!$          zdlw_fr_ksop05 = zredp05*zdw(i,1)*EXP(zdw1(i,1)*                        &
!!$                           (zporv(i,1)-zlw_fr_ksop05)/(zporv(i,1)-zadp(i,1)) )
!!$          zklw_fr_ksop05 = zredp05*zkw(i,1)*EXP(zkw1(i,1)*                        &
!!$                           (zporv(i,1)-zlw_fr_ksop05)/(zporv(i,1)-zadp(i,1)) )

         zdlw_fr_ksop05= zredp05*watrdiff_RT(zdw(i,1),zlw_fr_ksop05,&
                                  zdw1(i,1),zporv(i,1),zadp(i,1))
         zklw_fr_ksop05= zredp05*watrcon_RT(zkw(i,1),zlw_fr_ksop05,zkw1(i,1),zporv(i,1),zadp(i,1))

          ! coefficients for implicit flux computation
          z1dgam1     = zdt/zdzhs(1)
          zgam2p05    = zdlw_fr_ksop05/zdzms(2)
          zaga(i,1) = 0._ireals
          zagb(i,1) = 1._ireals+zalfa*zgam2p05*z1dgam1
          zagc(i,1) = -zalfa * zgam2p05*z1dgam1
          zagd(i,1) = zlw_fr(i,1) + zinfil(i)*z1dgam1/rho_w  &
                       -zklw_fr_ksop05*z1dgam1                     &
                       +(1._ireals - zalfa)* zgam2p05*z1dgam1*(zlw_fr_ksop1 - zlw_fr_kso)  &
                       +                     zgam2p05*z1dgam1*(zice_fr_ksop1-zice_fr_kso)

          ! explicit part of soil surface water flux:
          zflmg (i,1) = - zinfil(i)! boundary value for soil water transport
  END DO


! inner layers 2 <=kso<=ke_soil_hy-1
  DO kso =2,ke_soil_hy-1
!CDIR NODEP,VOVERTAKE,VOB
!DIR$ IVDEP, PREFERVECTOR
      DO ic = 1, icount_soil
        i=soil_list(ic)
! sedimentation and capillary transport in soil
            zice_fr_ksom1 = ziw_fr(i,kso-1)
            zice_fr_kso   = ziw_fr(i,kso  )
            zice_fr_ksop1 = ziw_fr(i,kso+1)
            zlw_fr_ksom1  = zlw_fr(i,kso-1)
            zlw_fr_kso    = zlw_fr(i,kso  )
            zlw_fr_ksop1  = zlw_fr(i,kso+1)
            ! interpolated scaled liquid water content at interface to layer
            ! above and below
            zlw_fr_ksom05 = 0.5_ireals*(zdzhs(kso-1)*zlw_fr_kso+   &
                                   zdzhs(kso)*zlw_fr_ksom1)/zdzms(kso)
            zlw_fr_ksop05 = 0.5_ireals*(zdzhs(kso+1)*zlw_fr_kso+   &
                                   zdzhs(kso)*zlw_fr_ksop1)/zdzms(kso+1)

            ! compute reduction factor for coefficients
            zfr_ice          = max (zice_fr_kso,zice_fr_ksom1)
            zredm05 = 1._ireals-zfr_ice/max (zlw_fr_kso+zice_fr_kso,zlw_fr_ksom1+zice_fr_ksom1)
            zfr_ice          = max (zice_fr_kso,zice_fr_ksop1)
            zredp05 = 1._ireals-zfr_ice/max (zlw_fr_kso+zice_fr_kso,zlw_fr_ksop1+zice_fr_ksop1)
!!$            zdlw_fr_ksom05= zredm05*zdw(i,kso)*EXP( zdw1(i,kso)*   &
!!$                               (zporv(i,kso)-zlw_fr_ksom05)/(zporv(i,kso)-zadp(i,kso)) )
!!$            zdlw_fr_ksop05= zredp05*zdw(i,kso)*EXP( zdw1(i,kso)*   &
!!$                               (zporv(i,kso)-zlw_fr_ksop05)/(zporv(i,kso)-zadp(i,kso)) )
!!$            zklw_fr_ksom05= zredm05*zkw(i,kso)*EXP( zkw1(i,kso)*   &
!!$                               (zporv(i,kso)-zlw_fr_ksom05)/(zporv(i,kso)-zadp(i,kso)) )
!!$            zklw_fr_ksop05= zredp05*zkw(i,kso)*EXP( zkw1(i,kso)*   &
!!$                               (zporv(i,kso)-zlw_fr_ksop05)/(zporv(i,kso)-zadp(i,kso)) )

            zdlw_fr_ksom05= zredm05*watrdiff_RT(zdw(i,kso),zlw_fr_ksom05,&
                                  zdw1(i,kso),zporv(i,kso),zadp(i,kso))

            zdlw_fr_ksop05= zredp05*watrdiff_RT(zdw(i,kso),zlw_fr_ksop05,&
                                  zdw1(i,kso),zporv(i,kso),zadp(i,kso))

            zklw_fr_ksom05= zredm05*watrcon_RT(zkw(i,kso),zlw_fr_ksom05,zkw1(i,kso),zporv(i,kso),zadp(i,kso))

            zklw_fr_ksop05= zredp05*watrcon_RT(zkw(i,kso),zlw_fr_ksop05,zkw1(i,kso),zporv(i,kso),zadp(i,kso))


            ! coefficients for implicit flux computation
            z1dgam1 = zdt/zdzhs(kso)
            zgam2m05  = zdlw_fr_ksom05/zdzms(kso)
            zgam2p05  = zdlw_fr_ksop05/zdzms(kso+1)
            zaga (i,kso) = -zalfa*zgam2m05*z1dgam1
            zagc (i,kso) = -zalfa*zgam2p05*z1dgam1
            zagb (i,kso) = 1._ireals +zalfa*(zgam2m05+zgam2p05)*z1dgam1
            zagd (i,kso) = zlw_fr(i,kso)+                               &
                                  z1dgam1*(-zklw_fr_ksop05+zklw_fr_ksom05)+ &
                                  (1._ireals-zalfa)*z1dgam1*                &
                                  (zgam2p05*(zlw_fr_ksop1-zlw_fr_kso  )     &
                                  -zgam2m05*(zlw_fr_kso  -zlw_fr_ksom1)   ) &
                                 +z1dgam1*                                  &
                                  (zgam2p05*(zice_fr_ksop1-zice_fr_kso  )   &
                                  -zgam2m05*(zice_fr_kso-zice_fr_ksom1)   )

            !soil water flux, explicit part, for soil water flux investigations
            ! only)
            zflmg(i,kso) = rho_w &
              &              *(zdlw_fr_ksom05*(zlw_fr_kso+zice_fr_kso-zlw_fr_ksom1-zice_fr_ksom1) &
              &              / zdzms(kso) - zklw_fr_ksom05)

            IF(kso==ke_soil_hy-1) THEN
              zflmg(i,kso+1)=rho_w &
                &       *(zdlw_fr_ksop05*(zlw_fr_ksop1+zice_fr_ksop1-zlw_fr_kso-zice_fr_kso) &
                &       / zdzms(kso+1) - zklw_fr_ksop05)
            ENDIF
      END DO
   END DO

!CDIR NODEP,VOVERTAKE,VOB
!DIR$ IVDEP
     DO ic = 1, icount_soil
        i=soil_list(ic)
          ! lowest active hydrological layer ke_soil_hy-1
          zice_fr_ksom1 = ziw_fr(i,ke_soil_hy-1)
          zice_fr_kso   = ziw_fr(i,ke_soil_hy  )
          zlw_fr_ksom1  = zlw_fr(i,ke_soil_hy-1)
          zlw_fr_kso    = zlw_fr(i,ke_soil_hy  )
          zlw_fr_ksom05 = 0.5_ireals*(zdzhs(ke_soil_hy-1)*zlw_fr_kso+ &
                              zdzhs(ke_soil_hy)*zlw_fr_ksom1)/zdzms(ke_soil_hy)

          zfr_ice          = max (zice_fr_kso,zice_fr_ksom1)
          zredm05 = 1._ireals-zfr_ice/max (zlw_fr_kso+zice_fr_kso,zlw_fr_ksom1+zice_fr_ksom1)
!!$          zdlw_fr_ksom05= zredm05*zdw(i,ke_soil_hy)*EXP( zdw1(i,ke_soil_hy)* &
!!$                            (zporv(i,ke_soil_hy)-zlw_fr_ksom05)/(zporv(i,ke_soil_hy)-zadp(i,ke_soil_hy)) )

          zdlw_fr_ksom05= zredm05*watrdiff_RT(zdw(i,ke_soil_hy),zlw_fr_ksom05,&
                                  zdw1(i,ke_soil_hy),zporv(i,ke_soil_hy),zadp(i,ke_soil_hy))


          z1dgam1 = zdt/zdzhs(ke_soil_hy)
          zgam2m05  = zdlw_fr_ksom05/zdzms(ke_soil_hy)
!!$   zklw_fr_ksom05= zredm05*zkw(i,ke_soil_hy)*EXP( zkw1(i,ke_soil_hy)* &
!!$               (zporv(i,ke_soil_hy)-zlw_fr_ksom05)/(zporv(i,ke_soil_hy)-zadp(i,ke_soil_hy)) )
!!$

          zklw_fr_ksom05= zredm05*watrcon_RT(zkw(i,ke_soil_hy),zlw_fr_ksom05,&
               zkw1(i,ke_soil_hy),zporv(i,ke_soil_hy),zadp(i,ke_soil_hy))

          zaga(i,ke_soil_hy) = -zalfa* zgam2m05*z1dgam1
          zagb(i,ke_soil_hy) = 1._ireals+ zalfa*zgam2m05*z1dgam1
          zagc(i,ke_soil_hy) = 0.0_ireals
          zagd(i,ke_soil_hy) = zlw_fr(i,ke_soil_hy)+z1dgam1*zklw_fr_ksom05 &
                            +(1._ireals-zalfa)*z1dgam1*                              &
                             zgam2m05*(zlw_fr_ksom1  - zlw_fr_kso)            &
                            +z1dgam1*                                         &
                             zgam2m05*(zice_fr_ksom1-zice_fr_kso )
  END DO

!CDIR NODEP,VOVERTAKE,VOB
!DIR$ IVDEP
  DO ic = 1, icount_soil
        i=soil_list(ic)
        ! generalized upper boundary condition
          zagc(i,1) = zagc(i,1)/zagb(i,1)
          zagd(i,1) = zagd(i,1)/zagb(i,1)
  END DO

  DO kso=2,ke_soil_hy-1
!CDIR NODEP,VOVERTAKE,VOB
!DIR$ IVDEP
  DO ic = 1, icount_soil
        i=soil_list(ic)
            zzz = 1._ireals/(zagb(i,kso) - zaga(i,kso)*zagc(i,kso-1))
            zagc(i,kso) = zagc(i,kso) * zzz
            zagd(i,kso) = (zagd(i,kso) - zaga(i,kso)*zagd(i,kso-1)) * zzz
      END DO
  END DO                ! soil layers

!CDIR NODEP,VOVERTAKE,VOB
!DIR$ IVDEP
  DO ic = 1, icount_soil
        i=soil_list(ic)
           zage(i,ke_soil_hy) = (zagd(i,ke_soil_hy)-zaga(i,ke_soil_hy)*  &
                             zagd(i,ke_soil_hy-1))/                          &
                            (zagb(i,ke_soil_hy) - zaga(i,ke_soil_hy)*      &
                             zagc(i,ke_soil_hy-1))
  END DO

  DO kso = ke_soil_hy-1,1,-1
!CDIR NODEP,VOVERTAKE,VOB
!DIR$ IVDEP
  DO ic = 1, icount_soil
        i=soil_list(ic)
            zage(i,kso)     = zagd(i,kso) - zagc(i,kso)*zage(i,kso+1)
            ! compute implicit part of new liquid water content and add existing
            ! ice content
            w_so_new(i,kso) = zage(i,kso)*zdzhs(kso) + w_so_ice_now(i,kso)
      END DO
  END DO                ! soil layers

!lowest active hydrological level
!CDIR NODEP,VOVERTAKE,VOB
!DIR$ IVDEP
  DO ic = 1, icount_soil
        i=soil_list(ic)
          ! boundary values ensure that the calculation below leaves the climate
          ! layer water contents unchanged compute implicit part of new liquid
          ! water content and add existing ice content
          w_so_new(i,ke_soil_hy) = zage(i,ke_soil_hy)*zdzhs(ke_soil_hy) + &
                          w_so_ice_now(i,ke_soil_hy)
  END DO

! to ensure vertical constant water concentration profile beginning at
! layer ke_soil_hy for energetic treatment only
! soil water climate layer(s)

  IF (itype_hydbound == 3) THEN
    ! ground water as lower boundary of soil column
    DO kso = ke_soil_hy+1,ke_soil+1
!CDIR NODEP,VOVERTAKE,VOB
  DO ic = 1, icount_soil
        i=soil_list(ic)
              w_so_new(i,kso) = zporv(i,kso)*zdzhs(kso)
        END DO
    END DO
  ELSE
    DO kso = ke_soil_hy+1,ke_soil+1
!CDIR NODEP,VOVERTAKE,VOB
  DO ic = 1, icount_soil
        i=soil_list(ic)
             w_so_new(i,kso) = w_so_new(i,kso-1)*zdzhs(kso)/zdzhs(kso-1)
       END DO
    END DO
  ENDIF

! combine implicit part of sedimentation and capillary flux with explicit part
! (for soil water flux investigations only)
  DO kso = 2,ke_soil+1
!CDIR NODEP,VOVERTAKE,VOB
!DIR$ IVDEP
     DO ic = 1, icount_soil
        i=soil_list(ic)
            zice_fr_ksom1 = ziw_fr(i,kso-1)
            zice_fr_kso   = ziw_fr(i,kso)
            zlw_fr_ksom1_new= w_so_new(i,kso-1)/zdzhs(kso-1) - zice_fr_ksom1
            zlw_fr_kso_new  = w_so_new(i,kso  )/zdzhs(kso  ) - zice_fr_kso
            zlw_fr_ksom1  = w_so_now(i,kso-1)/zdzhs(kso-1) - zice_fr_ksom1
            zlw_fr_kso    = w_so_now(i,kso  )/zdzhs(kso  ) - zice_fr_kso
            !... additionally for runoff_g at lower level of lowest active water
            ! layer calculated with (upstream) main level soil water content
            ! compute reduction factor for transport coefficients
            zfr_ice          = max (zice_fr_kso,zice_fr_ksom1)
            zredm05 = 1._ireals-zfr_ice/max (zlw_fr_kso+zice_fr_kso,zlw_fr_ksom1+zice_fr_ksom1)

          ! interpolated liquid water content at interface to layer above
            zlw_fr_ksom05 =0.5_ireals*(zdzhs(kso)*zlw_fr_ksom1+zdzhs(kso-1)*zlw_fr_kso) &
                              /zdzms(kso)
!!$            zdlw_fr_ksom05= zredm05*zdw(i,kso)*EXP(zdw1(i,kso) *  &
!!$                            (zporv(i,kso)-zlw_fr_ksom05)/(zporv(i,kso)-zadp(i,kso)) )
!!$            zklw_fr_ksom05= zredm05*zkw(i,kso) * EXP(zkw1(i,kso)* &
!!$                            (zporv(i,kso)-zlw_fr_ksom05)/(zporv(i,kso)-zadp(i,kso)) )
!!$

           zdlw_fr_ksom05= zredm05*watrdiff_RT(zdw(i,kso),zlw_fr_ksom05,&
                                  zdw1(i,kso),zporv(i,kso),zadp(i,kso))

           zklw_fr_ksom05= zredm05*watrcon_RT(zkw(i,kso),zlw_fr_ksom05,zkw1(i,kso),&
                             zporv(i,kso),zadp(i,kso))


            IF (kso> ke_soil_hy) zdlw_fr_ksom05=0.0_ireals   ! no flux gradient
                                                      ! contribution below 2.5m
            IF (kso> ke_soil_hy) zklw_fr_ksom05=0.0_ireals   ! no gravitation flux below 2.5m
            zflmg(i,kso) =                              &
                     (1._ireals-zalfa) * zflmg(i,kso) + & ! explicit flux component
                                zalfa  * rho_w *          & ! implicit flux component
                   (zdlw_fr_ksom05 * (zlw_fr_kso_new+zice_fr_kso-zlw_fr_ksom1_new-zice_fr_ksom1) &
                     /zdzms(kso) - zklw_fr_ksom05)
            zredm = 1._ireals-zice_fr_kso/(zlw_fr_kso+zice_fr_kso)
!!$            zklw_fr_kso_new = zredm*zkw(i,kso) * EXP(zkw1(i,kso)* &
!!$                              (zporv(i,kso) - zlw_fr_kso_new)/(zporv(i,kso) - zadp(i,kso)) )

            zklw_fr_kso_new = zredm*watrcon_RT(zkw(i,kso),&
                 zlw_fr_kso_new,zkw1(i,kso),zporv(i,kso),zadp(i,kso))

            ! actual gravitation water flux
            IF (w_so_new(i,kso).LT.1.01_ireals*zadp(i,kso)*zdzhs(kso)) THEN
              zklw_fr_kso_new=0._ireals
            ENDIF
            zrunoff_grav(i,kso) =  - rho_w * zklw_fr_kso_new

            ! ground water as lower boundary of soil column
            IF ((kso == ke_soil_hy+1).and.(itype_hydbound == 3)) THEN
               zdelta_sm=( zlw_fr_kso_new - zlw_fr_ksom1_new )

!!$               zdlw_fr_kso = zredm05*zdw(i,kso)*EXP(zdw1(i,kso) *  &
!!$                    (zporv(i,kso)-zlw_fr_kso_new)/(zporv(i,kso)-zadp(i,kso)) )
!!$               zklw_fr_kso = zredm05*zkw(i,kso) * EXP(zkw1(i,kso)* &
!!$                    (zporv(i,kso)-zlw_fr_kso_new)/(zporv(i,kso)-zadp(i,kso)) )
!!$               zklw_fr_ksom1 = zredm05*zkw(i,kso) * EXP(zkw1(i,kso)* &
!!$                    (zporv(i,kso)-zlw_fr_ksom1_new)/(zporv(i,kso)-zadp(i,kso)) )

               zdlw_fr_kso = zredm05*watrdiff_RT(zdw(i,kso),zlw_fr_kso_new,&
                                  zdw1(i,kso),zporv(i,kso),zadp(i,kso))
               zklw_fr_kso = zredm05*watrcon_RT(zkw(i,kso),&
                 zlw_fr_kso_new,zkw1(i,kso),zporv(i,kso),zadp(i,kso))
               zklw_fr_ksom1 = zredm05*watrcon_RT(zkw(i,kso),&
                 zlw_fr_ksom1_new,zkw1(i,kso),zporv(i,kso),zadp(i,kso))

               zdhydcond_dlwfr=( zklw_fr_kso - zklw_fr_ksom1 ) / zdelta_sm
               zrunoff_grav(i,ke_soil_hy)=zrunoff_grav(i,ke_soil_hy)+ &
                  zdhydcond_dlwfr / &
                  (1.0_ireals-exp(-zdhydcond_dlwfr/zdlw_fr_kso*0.5_ireals*zdzms(ke_soil_hy+1)))* &
                  zdelta_sm
            ENDIF
      END DO
  END DO


  DO  kso = 1,ke_soil
    ! utility variables used to avoid if-constructs in following loops
    zro_sfak = zsf_heav(0.5_ireals + REAL(msr_off - kso,ireals))  ! 1.0 for 'surface runoff'
    zro_gfak = 1._ireals - zro_sfak                               ! 1.0 for 'ground runoff'

    ! - Compute subsoil runoff (runoff_g) as drainage flux through bottom
    !   of layer ke_soil_hy (as suggested by the Rhone-Aggregation
    !   Experiment)
    ! - soil moisture gradient related flux is switched off below
    !   (i.e. only sedimentation flux allowed between ke_soil_hy and ke_soil_hy+1)
    zfmb_fak = MERGE(1.0_ireals, 0.0_ireals, kso==ke_soil_hy)

     ! sedimentation and capillary transport in soil
!CDIR NODEP,VOVERTAKE,VOB
!DIR$ IVDEP
      DO ic = 1, icount_soil
         i=soil_list(ic)
            ! first runoff calculation without consideration of
            ! evapotranspiration
            !zdwg   =  zflmg(i,kso+1) - zflmg(i,kso)
            !zdwg calculated above by flux divergence has to be aequivalent with
            zdwg =  (w_so_new(i,kso)/zdzhs(kso)-zw_fr(i,kso))*zdzhs(kso) &
                                                                   /zdtdrhw
            zdwg =  zdwg + zrunoff_grav(i,kso)*zfmb_fak
            zredfu =  MAX( 0.0_ireals, MIN( 1.0_ireals,(zw_fr(i,kso) -     &
                       zfcap(i,kso))/MAX(zporv(i,kso) - zfcap(i,kso),zepsi)) )
            zredfu = zsf_heav(zdwg)*zredfu
            zro    = zdwg*zredfu
            zdwg   = zdwg*(1._ireals - zredfu)

            ! add evaporation (znlgw1f: first layer only)
            ! and transpiration (for each layer)
            zdwg   = zdwg + znlgw1f(kso) * zesoil(i) + ztrang (i,kso)
            zwgn   = zw_fr(i,kso) + zdtdrhw*zdwg/zdzhs(kso)
            zro2   = zrhwddt*zdzhs(kso)*MAX(0.0_ireals, zwgn - zporv(i,kso))
            zkorr  = zrhwddt*zdzhs(kso)*MAX(0.0_ireals, zadp(i,kso) - zwgn )

            ! SBr, CHa: prevent the subroutine terra_multlay from avoiding soil
            ! water content below air dryness point by setting zkorr=0
            zkorr  = 0.0_ireals
            zdwgdt(i,kso)= zdwg + zkorr - zro2
            ! SBr, CHa: prevent the subroutine terra_multlay from avoiding soil
            ! water content changing with time by setting tendency zdwgdt=0.
            zdwgdt(i,kso) = 0.0_ireals
            zro    = zro      + zro2
            runoff_s(i) = runoff_s(i) + zro*zro_sfak*zroffdt
            runoff_g(i) = runoff_g(i) + zro*zro_gfak*zroffdt
            ! runoff_g reformulation:
            runoff_g(i) = runoff_g(i) - (zrunoff_grav(i,kso) * zfmb_fak &
                                          + zkorr) * zroffdt
      END DO
   END DO         ! end loop over soil layers



!------------------------------------------------------------------------------
! Section II.5: Soil surface heat flux (thermal forcing)
!------------------------------------------------------------------------------

  IF (lmulti_snow) THEN

    DO i = istarts, iends
        ! Estimate thermal surface fluxes:
        ! Estimate thermal surface fluxes over snow covered and snow free
        ! part of surface based on area mean values calculated in radiation
        ! code (positive = downward)
        zgstr =   sigma*(1._ireals - Ctalb) * ( (1._ireals - zf_snow(i))* &
                  zts(i) + zf_snow(i)*ztsnow_mult(i,1) )**4 + thbs(i)
        zthsnw(i) = - sigma*(1._ireals - Ctalb)*ztsnow_mult(i,1)**4 + zgstr
        zthsoi(i) = - sigma*(1._ireals - Ctalb)*zts(i)**4 + zgstr
        ! the estimation of the solar component would require the availability
        ! of the diffuse and direct components of the solar flux
        !
        ! Forcing for snow-free soil:
        ! (evaporation, transpiration, formation of dew and rime are already
        !  weighted by correspondind surface fraction)
        ! net radiation, sensible and latent heat flux

        zrnet_s(i) = sobs(i) + zthsoi(i)
        zshfl_s(i) = cp_d*zrhoch(i) * (zth_low(i) - zts(i))
        zlhfl_s(i) = (zts_pm(i)*lh_v + (1._ireals-zts_pm(i))*lh_s)*zverbo(i) &
                     / MAX(zepsi,(1._ireals - zf_snow(i)))  ! take out (1-f) scaling
!DR start
        zqhfl_s(i) = zverbo(i)/ MAX(zepsi,(1._ireals - zf_snow(i)))  ! take out (1-f) scaling
!DR end
    END DO

    DO i = istarts, iends
        ! thawing of snow falling on soil with Ts > T0
        IF (ztsnow_pm(i)*zrs(i) > 0.0_ireals) THEN
          ! snow fall on soil with T>T0, snow water content increases
          ! interception store water content
          zsprs  (i) = - lh_f*zrs(i)
          zdwidt (i) = zdwidt (i) + zrs(i)
          zdwsndt(i) = zdwsndt(i) - zrs(i)

          ! avoid overflow of interception store, add possible excess to
          ! surface run-off
          zwinstr(i)      = zwin(i) + zdwidt(i)*zdtdrhw
          IF (zwinstr(i) > zwimax(i)) THEN  ! overflow of interception store
            zro        = (zwinstr(i) - zwimax(i))*zrhwddt
            zdwidt(i)= zdwidt(i) - zro
            runoff_s(i)  = runoff_s(i) + zro*zroffdt
          ENDIF                       ! overflow of interception store

        ! freezing of rain falling on soil with Ts < T0  (black-ice !!!)
        ELSEIF ((1._ireals-zts_pm(i))*zrr(i) > 0.0_ireals) THEN
          zsprs  (i) = lh_f*zrr(i)
          ! keep freezing rain in interception storage rather than shifting it to snow
         ! zdwidt (i) = zdwidt (i) - zrr(i)
         ! zdwsndt(i) = zdwsndt(i) + zrr(i)
        ELSE
          zsprs  (i) = 0.0_ireals
        END IF

!       Influence of heatflux through snow on total forcing:
        zdwsndt(i) = zdwsndt(i)*zsf_heav(zdwsndt(i) - zepsi/zdtdrhw)
        zwsnew(i)  = zwsnow(i) + zdwsndt(i)*zdtdrhw
        IF (zwsnew(i).GT.zepsi) THEN

          zrho_snowf = crhosminf+(crhosmaxf-crhosminf)* (zth_low(i)-csnow_tmin) &
                                                  /(t0_melt          -csnow_tmin)
          zrho_snowf = MAX(crhosminf,MIN(crhosmaxf,zrho_snowf))

          IF(zdwsndt(i)-zrs(i)-zrr(i).GT.0.0_ireals) THEN

            wtot_snow_now(i,1) = max(wtot_snow_now(i,1) + zdwsndt(i)*zdtdrhw, &
              &                           0.0_ireals)

            zhm_snow(i,1) = zhm_snow(i,1) - (zdwsndt(i)-zrs(i)- &
                              zrr(i))*zdt/rho_i/2._ireals-  &
              zrs(i)*zdt/zrho_snowf/2._ireals- zrr(i)*zdt/rho_i/2._ireals
            zdzh_snow(i,1) = zdzh_snow(i,1) + (zdwsndt(i)-zrs(i)-zrr(i))*zdt/rho_i +  &
              zrs(i)*zdt/zrho_snowf + zrr(i)*zdt/rho_i

            rho_snow_mult_now(i,1) = max(wtot_snow_now(i,1)*rho_w/zdzh_snow(i,1), &
              &                              0.0_ireals)
          ELSE

            wtot_snow_now(i,1) = max(wtot_snow_now(i,1) + (zrs(i)+zrr(i))*zdtdrhw, &
              &                          0.0_ireals)

            zhm_snow(i,1)  = zhm_snow(i,1) - zrs(i)*zdt/zrho_snowf/2._ireals- &
                               zrr(i)*zdt/rho_i/2._ireals
            zdzh_snow(i,1) = zdzh_snow(i,1) + zrs(i)*zdt/zrho_snowf + zrr(i)*zdt/rho_i

            IF(wtot_snow_now(i,1) .GT. 0._ireals) THEN
              rho_snow_mult_now(i,1) = max(wtot_snow_now(i,1)*rho_w/zdzh_snow(i,1), &
                &                              0.0_ireals)

              wtot_snow_now(i,1) = max(wtot_snow_now(i,1) &
                &                      + (zdwsndt(i)-zrs(i)-zrr(i))*zdtdrhw,0.0_ireals)

              zhm_snow(i,1)  = zhm_snow(i,1) - (zdwsndt(i)-zrs(i)-zrr(i)) &
                &                *zdt/rho_snow_mult_now(i,1)/2._ireals
              zdzh_snow(i,1) = zdzh_snow(i,1) + (zdwsndt(i)-zrs(i)-zrr(i)) &
                &                *zdt/rho_snow_mult_now(i,1)
            ELSE
              rho_snow_mult_now(i,1) = 0.0_ireals
              zdzh_snow(i,1) = 0.0_ireals
            END IF

          END IF
        END IF
        h_snow_now(i) = 0.0_ireals
        sum_weight(i) = 0.0_ireals
!      END IF          ! land-points only
    END DO

!!$IF (msg_level >= 14) THEN
!!$  DO i = istarts, iends
!!$  IF (soiltyp_subs(i) == 1) THEN  !1=glacier and Greenland
!!$    IF ( ABS( zshfl_snow(i) )  >  500.0  .OR. &
!!$         ABS( zlhfl_snow(i) )  > 2000.0 ) THEN
!!$      write(*,*) 'hello mo_soil_ml 1: ', zshfl_s(i),cp_d, zrhoch(i),zth_low(i),zts(i), &
!!$        '  ......  ', zlhfl_s(i),zts_pm(i),lh_v,          lh_s,zverbo(i),zf_snow(i), &
!!$        '  ......  ', tch(i), tcm(i)
!!$    ENDIF
!!$  ENDIF
!!$  END DO
!!$ENDIF

    DO ksn = 1,ke_snow
      DO i = istarts, iends
!        IF (llandmask(i)) THEN          ! land-points only
          IF (zwsnew(i).GT.zepsi) THEN
            h_snow_now(i) = h_snow_now(i) + zdzh_snow(i,ksn)
          END IF
!        END IF          ! land-points only
      END DO
    END DO

    k = MIN(2,ke_snow-1)
    DO ksn = 1,ke_snow
      DO i = istarts, iends
        IF(zwsnew(i) .GT. zepsi) THEN
          IF(zwsnow(i) .GT. zepsi) THEN
            IF (ksn == 1) THEN ! Limit top layer to max_toplaydepth
              zhh_snow(i,ksn) = -MAX( h_snow_now(i)-max_toplaydepth, h_snow_now(i)/ke_snow*(ke_snow-ksn) )
            ELSE IF (ksn == 2 .AND. ke_snow > 2) THEN ! Limit second layer to 8*max_toplaydepth
              zhh_snow(i,ksn) = MIN( 8._ireals*max_toplaydepth+zhh_snow(i,1), zhh_snow(i,1)/(ke_snow-1)*(ke_snow-ksn) )
            ELSE ! distribute the remaining snow equally among the layers
              zhh_snow(i,ksn) = zhh_snow(i,k)/(ke_snow-k)*(ke_snow-ksn)
            ENDIF
          ELSE ! a newly generated snow cover will not exceed max_toplaydepth
            zhh_snow(i,ksn) = -h_snow_now(i)/ke_snow*(ke_snow-ksn)
          END IF
        END IF
      END DO
    END DO

    DO ksn = ke_snow,1,-1
      DO i = iends, istarts, -1
        IF(zwsnew(i) .GT. zepsi .AND. zwsnow(i) .GT. zepsi) THEN
          dz_old(i,ksn) = zdzh_snow(i,ksn)
          z_old(i,ksn) = -sum_weight(i) - zdzh_snow(i,ksn)/2._ireals
          sum_weight(i) = sum_weight(i) + zdzh_snow(i,ksn)
        END IF
      END DO
    END DO

    DO i = istarts, iends
!      IF (llandmask(i)) THEN  ! for landpoints only
        IF(zwsnew(i) .GT. zepsi) THEN
          zhm_snow (i,1) = (-h_snow_now(i) + zhh_snow(i,1))/2._ireals
          zdzh_snow(i,1) = zhh_snow(i,1) + h_snow_now(i)            !layer thickness betw. half levels of uppermost snow layer
          zdzm_snow(i,1) = zhm_snow(i,1) + h_snow_now(i)            !layer thickness between snow surface and main level of uppermost layer
          IF(zwsnow(i) .GT. zepsi) THEN
            IF(dz_old(i,1).ne.0..and.rho_snow_mult_now(i,1).ne.0.) THEN
              wliq_snow_now(i,1) = wliq_snow_now(i,1)/dz_old(i,1)
            END IF
          END IF
        END IF
!      END IF          ! land-points only
    END DO
    DO ksn = 2,ke_snow
      DO i = istarts, iends
!        IF (llandmask(i)) THEN  ! for landpoints only
          IF(zwsnew(i) .GT. zepsi) THEN
            zhm_snow(i,ksn) = (zhh_snow(i,ksn) + zhh_snow(i,ksn-1))/2._ireals
            zdzh_snow(i,ksn) = zhh_snow(i,ksn) - zhh_snow(i,ksn-1) ! layer thickness betw. half levels
            zdzm_snow(i,ksn) = zhm_snow(i,ksn) - zhm_snow(i,ksn-1) ! layer thickness betw. main levels
            IF(zwsnow(i) .GT. zepsi) THEN
              IF(dz_old(i,ksn).ne.0..and.rho_snow_mult_now(i,ksn).ne.0.) THEN
                wliq_snow_now(i,ksn) = wliq_snow_now(i,ksn)/dz_old(i,ksn)
              END IF
            END IF
          END IF
!        END IF          ! land-points only
      END DO
    END DO

    DO ksn = ke_snow,1,-1
      DO i = iends, istarts, -1
        t_new  (i,ksn) = 0.0_ireals
        rho_new(i,ksn) = 0.0_ireals
        wl_new (i,ksn) = 0.0_ireals
      END DO

      DO k = ke_snow,1,-1
        DO i = iends, istarts, -1
!          IF (llandmask(i)) THEN  ! for landpoints only
            IF(zwsnew(i) .GT. zepsi .AND. zwsnow(i) .GT. zepsi) THEN

              weight = MIN(dz_old(i, k), &
                   &       z_old(i, k) + dz_old(i, k)*0.5_ireals &
                   &         - zhm_snow(i, ksn) + zdzh_snow(i, ksn)*0.5_ireals,&
                   &       zhm_snow(i, ksn) + zdzh_snow(i,ksn)*0.5_ireals &
                   &         - z_old(i, k) + dz_old(i, k)*0.5_ireals, &
                   &       zdzh_snow(i,ksn))

              weight = (weight + ABS(weight)) * 0.5_ireals
              weight = weight / zdzh_snow(i,ksn)
              t_new  (i,ksn) = t_new  (i,ksn) + ztsnow_mult  (i,k      )*weight
              rho_new(i,ksn) = rho_new(i,ksn) + rho_snow_mult_now(i,k)*weight
              wl_new (i,ksn) = wl_new (i,ksn) + wliq_snow_now(i,k)*weight
            END IF
!          END IF          ! land-points only
        END DO
      END DO
    END DO

                   ! MIN(z_old(i, k) + dz_old(i, k)/2._ireals, &
                   ! &   zhm_snow(i, ksn) + zdzh_snow(i,ksn)/2._ireals) &
                   ! - MAX(z_old(i, k) - dz_old(i, k)/2._ireals, &
                   ! &     zhm_snow(i, ksn) - zdzh_snow(i, ksn)/2._ireals), &

    DO ksn = ke_snow,1,-1
      DO i = iends, istarts, -1
!        IF (llandmask(i)) THEN  ! for landpoints only
          IF(zwsnew(i) .GT. zepsi) THEN
            IF(zwsnow(i) .GT. zepsi) THEN
              ztsnow_mult  (i,ksn      ) = t_new  (i,ksn)
              rho_snow_mult_now(i,ksn) = rho_new(i,ksn)
              wtot_snow_now    (i,ksn) = rho_new(i,ksn)*zdzh_snow(i,ksn)/rho_w
              wliq_snow_now    (i,ksn) = wl_new (i,ksn)*zdzh_snow(i,ksn)
            ELSE ! Remark: if there was now snow in the previous time step, snow depth will not exceed the limit for equipartitioning
              ztsnow_mult  (i,ksn      ) = t_s_now(i)
              rho_snow_mult_now(i,ksn) = rho_snow_mult_now(i,1)
              wtot_snow_now    (i,ksn) = zwsnew(i)/ke_snow
            END IF
          END IF
!        END IF          ! land-points only
      END DO
    END DO

    ! heat conductivity of snow as funtion of water content
    DO ksn = 1, ke_snow
        DO i = istarts, iends
!          IF (llandmask(i)) THEN  ! for landpoints only
            IF (zwsnew(i).GT.zepsi) THEN
              zalas_mult(i,ksn) = 2.22_ireals*EXP(1.88_ireals*LOG(rho_snow_mult_now(i,ksn)/rho_i))
            END IF
!          END IF          ! land-points only
        END DO
    END DO

    DO i = istarts, iends
!      IF (llandmask(i)) THEN  ! for landpoints only
        IF (zwsnew(i).GT.zepsi) THEN
           zgsb(i) = ((zalas_mult(i,ke_snow)*(-zhm_snow(i,ke_snow))+hzalam(i,1)*zdzms(1))/ &
                       (-zhm_snow(i,ke_snow)+zdzms(1)) * &
                       (ztsnow_mult(i,ke_snow) - t_so_now(i,1))/(-zhm_snow(i,ke_snow) &
                       +zdzms(1)))*zf_snow(i)

!          ! GZ: use formulation of single-layer snow model, which is numerically more stable
!          zgsb(i) = zalas_mult(i,ke_snow)*(ztsnow_mult(i,ke_snow) - t_so_now(i,1))/ &
!            MAX(-zhm_snow(i,ke_snow),cdsmin)

        END IF

        ! total forcing for uppermost soil layer
!<em new solution
!        zfor_s(i) = ( zrnet_s(i) + zshfl_s(i) + zlhfl_s(i) + zsprs(i) ) * (1._ireals - zf_snow(i)) &
!                    + (1._ireals-ztsnow_pm(i)) * zgsb(i)
        zfor_s(i) = ( zrnet_s(i) + zshfl_s(i) + zlhfl_s(i) + zsprs(i) ) * (1._ireals - zf_snow(i))
!em>

        IF(zwsnew(i) .GT. zepsi) THEN
          IF(zextinct(i,1).gt.0.0_ireals) THEN
            zrnet_snow = zthsnw(i)
          ELSE
            zrnet_snow = sobs(i) + zthsnw(i)
          END IF
        ELSE
          zrnet_snow = sobs(i) + zthsnw(i)
        END IF
        zshfl_snow(i) = zrhoch(i)*cp_d*(zth_low(i) - ztsnow_mult(i,1))
        zlhfl_snow(i) = lh_s*zversn(i)
        zqhfl_snow(i) = zversn(i)
        zfor_snow_mult(i)  = (zrnet_snow + zshfl_snow(i) + zlhfl_snow(i) + lh_f*zrr(i))*zf_snow(i)

!      END IF          ! land-points only
    END DO

  ELSE  ! single-layer snow model

    DO i = istarts, iends
!      IF (llandmask(i)) THEN          ! land-points only
        ! Estimate thermal surface fluxes:
        ! Estimate thermal surface fluxes over snow covered and snow free
        ! part of surface based on area mean values calculated in radiation
        ! code (positive = downward)
        zgstr =   sigma*(1._ireals - Ctalb) * ( (1._ireals - zf_snow(i))* &
                  zts(i) + zf_snow(i)*ztsnow(i) )**4 + thbs(i)
        zthsnw(i) = - sigma*(1._ireals - Ctalb)*ztsnow(i)**4 + zgstr
        zthsoi(i) = - sigma*(1._ireals - Ctalb)*zts(i)**4 + zgstr
        ! the estimation of the solar component would require the availability
        ! of the diffuse and direct components of the solar flux
        !
        ! Forcing for snow-free soil:
        ! (evaporation, transpiration, formation of dew and rime are already
        !  weighted by correspondind surface fraction)
        ! net radiation, sensible and latent heat flux

        zrnet_s(i) = sobs(i) + zthsoi(i)
        zshfl_s(i) = cp_d*zrhoch(i) * (zth_low(i) - zts(i))
        zlhfl_s(i) = (zts_pm(i)*lh_v + (1._ireals-zts_pm(i))*lh_s)*zverbo(i) &
                     / MAX(zepsi,(1._ireals - zf_snow(i)))  ! take out (1-f) scaling
!DR start
        zqhfl_s(i) = zverbo(i)/ MAX(zepsi,(1._ireals - zf_snow(i)))  ! take out (1-f) scaling
!DR end
        zsprs  (i) = 0.0_ireals
        ! thawing of snow falling on soil with Ts > T0
        IF (ztsnow_pm(i)*zrs(i) > 0.0_ireals) THEN
          ! snow fall on soil with T>T0, snow water content increases
          ! interception store water content
          zsprs  (i) = - lh_f*zrs(i)
          zdwidt (i) = zdwidt (i) + zrs(i)
          zdwsndt(i) = zdwsndt(i) - zrs(i)

          ! avoid overflow of interception store, add possible excess to
          ! surface run-off
          zwinstr(i) = zwin(i) + zdwidt(i)*zdtdrhw
          IF (zwinstr(i) > zwimax(i)) THEN  ! overflow of interception store
            zro         = (zwinstr(i) - zwimax(i))*zrhwddt
            zdwidt(i)   = zdwidt(i) - zro
            runoff_s(i) = runoff_s(i) + zro*zroffdt
          ENDIF                       ! overflow of interception store

        ! freezing of rain falling on soil with Ts < T0  (black-ice !!!)
        ELSEIF (zwsnow(i) == 0.0_ireals .AND.                            &
               (1._ireals-ztsnow_pm(i))*zrr(i) > 0.0_ireals) THEN
          zsprs  (i) = MIN(lh_f*zrr(i),(t0_melt-t_s_now(i))*zroc(i,1)*zdzhs(1)/zdt)
          ! keep freezing rain in interception storage rather than shifting it to snow
         ! zdwidt (i) = zdwidt (i) - zrr(i)
         ! zdwsndt(i) = zdwsndt(i) + zrr(i)
        END IF

!       Influence of heatflux through snow on total forcing:
        zwsnew(i)   = zwsnow(i) + zdwsndt(i)*zdtdrhw
        IF (zwsnew(i).GT.zepsi) THEN
!         heat conductivity of snow as funtion of water content
! BR      zalas  = MAX(calasmin,MIN(calasmax, calasmin + calas_dw*zwsnow(i)))
!
! BR 7/2005 Introduce new dependency of snow heat conductivity on snow density
!
          zalas  = 2.22_ireals*EXP(1.88_ireals*LOG(zrho_snow(i)/rho_i))

! BR 11/2005 Use alternative formulation for heat conductivity by Sun et al., 1999
!            The water vapour transport associated conductivity is not included.

!        zalas   = 0.023_ireals+(2.290_ireals-0.023_ireals)* &
!                               (7.750E-05_ireals*rho_snow(i,nx) + &
!                                1.105E-06_ireals*prho_snow(i,nx)**2)

!          zgsb(i) = zalas*(ztsnow(i) - zts(i))/zdz_snow_fl(i)
          zgsb(i) = zalas*(ztsnow(i) - zts(i))/zdz_snow_fl(i)
        END IF

        ! total forcing for uppermost soil layer
        zfor_s(i) = ( zrnet_s(i) + zshfl_s(i) + zlhfl_s(i) ) &
                         * (1._ireals - zf_snow(i)) + zsprs(i) &
                    + zf_snow(i) * (1._ireals-ztsnow_pm(i)) * zgsb(i)

!      END IF          ! land-points only
    END DO

IF (msg_level >= 19) THEN
  DO i = istarts, iends
  IF (soiltyp_subs(i) == 1) THEN  !1=glacier and Greenland
    IF ( ABS( zshfl_s(i) )  >  500.0_ireals  .OR. &
         ABS( zlhfl_s(i) )  > 2000.0_ireals ) THEN
      write(*,*) 'hello mo_soil_ml 2: ', zshfl_s(i), zrhoch(i),zth_low(i),t(i),zts(i), &
        '  ...LHF...  ',                 zlhfl_s(i), zts_pm(i),zverbo(i),zf_snow(i),qv(i),qv_s(i), &
        '  ...CH,CM...  ', tch(i), tcm(i)
    ENDIF
  ENDIF
  END DO
ENDIF

  ENDIF ! lmulti_snow



!------------------------------------------------------------------------------
! Section II.6: Solution of the heat conduction equation, freezing/melting
!               of soil water/ice (optionally)
!------------------------------------------------------------------------------

  ! EM: If the single-layer snow model is used, nothing changes;
  ! if the multi-layer snow model is used and zf_snow(i) == 1,
  ! then the heat conduction equation for the whole column "soil + snow" is solved
  ! (see below, after the statement IF (lmulti_snow) THEN);
  ! if the multi-layer snow model is used, but zf_snow(i) < 1._ireals,
  ! then the two partial temperature updates are computed: first, the heat conduction equation
  ! for the snow-free surface is solved, second, the heat conduction equation
  ! for the whole column "soil + snow" for snow-covered surface is solved.
  ! Then, the two updates are merged together.

  DO i = istarts, iends
    sn_frac(i) = zf_snow(i)
    IF (zwsnow(i) < zepsi .AND. zwsnew(i) >= zepsi) THEN
      sn_frac(i) = 0.01_ireals
    ELSEIF (zwsnow(i) >= zepsi .AND. zwsnew(i) < zepsi) THEN
      sn_frac(i) = 0._ireals
    ENDIF
  END DO

  DO kso = 1, ke_soil+1
    DO i = istarts, iends
      t_so_free_new(i,kso) = t_so_now(i,kso)
      t_so_snow_new(i,kso) = t_so_now(i,kso)
    END DO
  END DO        ! soil layers

  DO kso = 2, ke_soil
    DO i = istarts, iends
      IF (.NOT. lmulti_snow .OR. sn_frac(i) < 1._ireals) THEN
        ! for heat conductivity: zalam is now 3D
        zakb1 = zalam(i,kso-1)/zroc(i,kso)
        zakb2 = zalam(i,kso  )/zroc(i,kso)
        zaga(i,kso) = -zalfa*zdt*zakb1/(zdzhs(kso)*zdzms(kso))
        zagc(i,kso) = -zalfa*zdt*zakb2/(zdzhs(kso)*zdzms(kso+1))
        zagb(i,kso) = 1._ireals - zaga(i,kso) - zagc(i,kso)
        zagd(i,kso) = t_so_now(i,kso) +                                     &
               (1._ireals - zalfa)*( - zaga(i,kso)/zalfa*t_so_now(i,kso-1)+ &
               (zaga(i,kso)/zalfa + zagc(i,kso)/zalfa)*t_so_now(i,kso) -  &
                zagc(i,kso)/zalfa*t_so_now(i,kso+1)  )
      END IF
    END DO
  END DO        ! soil layers

  DO i = istarts, iends
!    IF (.NOT. lmulti_snow .OR. (zf_snow(i) < 1._ireals .AND. zwsnew(i).GT.zepsi)) THEN
    IF (.NOT. lmulti_snow .OR. sn_frac(i) < 1._ireals) THEN
      ! for heat conductivity: zalam is now 3D: here we need layer 1
      zakb1 = hzalam(i,1)/zroc(i,1)
      zakb2 =  zalam(i,1)/zroc(i,1)
      zaga(i,1) = -zalfa*zdt*zakb1/(zdzhs(1)*zdzms(1))
      zagc(i,1) = -zalfa*zdt*zakb2/(zdzhs(1)*zdzms(2))
      zagb(i,1) = 1._ireals - zaga(i,1) - zagc(i,1)
      zagd(i,1) = t_so_now(i,1) + (1._ireals - zalfa)* (                 &
                      - zaga(i,1)/zalfa * t_s_now(i) +                     &
                      (zaga(i,1) + zagc(i,1))/zalfa * t_so_now(i,1) -    &
                       zagc(i,1)/zalfa * t_so_now(i,2)   )
      zaga(i,0) = 0.0_ireals
      zagb(i,0) = zalfa
      zagc(i,0) = -zalfa
      ! EM: In the case of multi-layer snow model, zfor_s(i) does not include the heat conductivity flux
      ! between soil and snow (zgsb). It will be accounted for at the next semi-step, see below.
      zagd(i,0)    = zdzms(1) * zfor_s(i)/hzalam(i,1)+(1._ireals-zalfa)* &
                      (t_so_now(i,1) - t_s_now(i))
      zaga(i,ke_soil+1) = 0.0_ireals
      zagb(i,ke_soil+1) = 1.0_ireals
      zagc(i,ke_soil+1) = 0.0_ireals
      zagd(i,ke_soil+1) = t_so_now(i,ke_soil+1)
    END IF
  END DO

  DO i = istarts, iends
    IF (.NOT. lmulti_snow .OR. sn_frac(i) < 1._ireals) THEN
      zagc(i,0) = zagc(i,0)/zagb(i,0)
      zagd(i,0) = zagd(i,0)/zagb(i,0)
    END IF
  END DO

  DO kso = 1, ke_soil
    DO i = istarts, iends
      IF (.NOT. lmulti_snow .OR. sn_frac(i) < 1._ireals) THEN
        zzz = 1._ireals/(zagb(i,kso) - zaga(i,kso)*zagc(i,kso-1))
        zagc(i,kso) = zagc(i,kso) * zzz
        zagd(i,kso) = (zagd(i,kso) - zaga(i,kso)*zagd(i,kso-1)) * zzz
      END IF
    END DO
  END DO                ! soil layers

  DO i = istarts, iends
    IF (.NOT. lmulti_snow .OR. sn_frac(i) < 1._ireals) THEN
      zage(i,ke_soil+1) = (zagd(i,ke_soil+1) - zaga(i,ke_soil+1)*       &
                             zagd(i,ke_soil))/                              &
                            (zagb(i,ke_soil+1) - zaga(i,ke_soil+1)*       &
                             zagc(i,ke_soil))
      t_so_free_new(i,ke_soil+1) = zage(i,ke_soil+1) ! climate value, unchanged
    END IF
  END DO

  DO kso = ke_soil,0,-1
    DO i = istarts, iends
      IF (.NOT. lmulti_snow .OR. sn_frac(i) < 1._ireals) THEN
        zage(i,kso)     = zagd(i,kso) - zagc(i,kso)*zage(i,kso+1)
        ! The surface temperature computed by t_so(i,0,nnew)=zage(i,0) is
        ! presently unused
        t_so_free_new(i,kso) = zage(i,kso)
      END IF
    END DO
  END DO                ! soil layers

  IF (lmulti_snow) THEN

    ! If there is snow, the solution of the heat conduction equation
    ! goes through the whole column "soil+snow"
    DO i = istarts, iends
      IF (sn_frac(i) > 0._ireals) THEN
        ! Uppermost snow layer, Neumann boundary condition
        zrocs(i) = (wliq_snow_now(i,1)/wtot_snow_now(i,1)*chc_w + &
              (wtot_snow_now(i,1)-wliq_snow_now(i,1))/wtot_snow_now(i,1)*chc_i)*rho_snow_mult_now(i,1)
        zakb      = zalas_mult(i,1)/zrocs(i)
        zaga(i,1) = -zalfa*zdt*zakb/(zdzh_snow(i,1)*zdzm_snow(i,1))
        zagc(i,1) = -zalfa*zdt*zakb/(zdzh_snow(i,1)*zdzm_snow(i,2))
        zagb(i,1) = 1._ireals - zaga(i,1) - zagc(i,1)
        zagd(i,1) = ztsnow_mult(i,1) + (1._ireals - zalfa)* (-zaga(i,1)/zalfa * t_s_now(i) + &
                    (zaga(i,1) + zagc(i,1))/zalfa * t_so_now(i,1) - zagc(i,1)/zalfa * t_so_now(i,2) )
        zaga(i,0) = 0.0_ireals
        zagb(i,0) = zalfa
        zagc(i,0) = -zalfa
        zagd(i,0) = zdzm_snow(i,1) * zfor_snow_mult(i)/zalas_mult(i,1)+(1._ireals-zalfa)* &
                   (t_so_now(i,1) - t_s_now(i))
        zagc(i,0) = zagc(i,0)/zagb(i,0)
        zagd(i,0) = zagd(i,0)/zagb(i,0)

        ! Lowermost soil layer, Dirichlet boundary condition
        zaga(i,ke_snow + ke_soil+1) = 0.0_ireals
        zagb(i,ke_snow + ke_soil+1) = 1.0_ireals
        zagc(i,ke_snow + ke_soil+1) = 0.0_ireals
        zagd(i,ke_snow + ke_soil+1) = t_so_now(i,ke_soil+1)

        ! Lowermost snow layer, special treatment
        zrocs(i) = (wliq_snow_now(i,ke_snow)/wtot_snow_now(i,ke_snow)*chc_w +  &
          & (wtot_snow_now(i,ke_snow)-wliq_snow_now(i,ke_snow))/                    &
          & wtot_snow_now(i,ke_snow)*chc_i)*rho_snow_mult_now(i,ke_snow)
        zakb = zalas_mult(i,ke_snow)/zrocs(i)
        zaga(i,ke_snow) = -zalfa*zdt*zakb/(zdzh_snow(i,ke_snow)*zdzm_snow(i,ke_snow))
        zakb = (zalas_mult(i,ke_snow)/zrocs(i)*(-zhm_snow(i,ke_snow))+hzalam(i,1)/ &
                zroc(i,1)*zmls(1))/(zmls(1)-zhm_snow(i,ke_snow))
        zagc(i,ke_snow) = -zalfa*zdt*zakb/(zdzh_snow(i,ke_snow)*(zmls(1)-zhm_snow(i,ke_snow)))
        zagb(i,ke_snow) = 1.0_ireals - zaga(i,ke_snow) - zagc(i,ke_snow)
        zagd(i,ke_snow) = ztsnow_mult(i,ke_snow) + &
          &     (1._ireals - zalfa)*( - zaga(i,ke_snow)/zalfa*ztsnow_mult(i,ke_snow-1) + &
          &     (zaga(i,ke_snow)/zalfa + zagc(i,ke_snow)/zalfa)*ztsnow_mult(i,ke_snow) - &
          &     zagc(i,ke_snow)/zalfa*t_so_now(i,1)  )

        ! Uppermost soil layer, special treatment
        zakb = (zalas_mult(i,ke_snow)/zrocs(i)*(-zhm_snow(i,ke_snow))+     &
                hzalam(i,1)/zroc(i,1)*zmls(1))/(zmls(1)-zhm_snow(i,ke_snow))
        zaga(i,ke_snow+1) = -zalfa*zdt*zakb/(zdzhs(1)*(zmls(1)-zhm_snow(i,ke_snow)))
        zakb = hzalam(i,1)/zroc(i,1)
        zagc(i,ke_snow+1) = -zalfa*zdt*zakb/(zdzhs(1)*zdzms(2))
        zagb(i,ke_snow+1) = 1._ireals - zaga(i,ke_snow+1) - zagc(i,ke_snow+1)
        zagd(i,ke_snow+1) = t_so_now(i,1) + &
          &    (1._ireals - zalfa)*( - zaga(i,ke_snow+1)/zalfa*ztsnow_mult(i,ke_snow) + &
          &    (zaga(i,ke_snow+1)/zalfa + zagc(i,ke_snow+1)/zalfa)*t_so_now(i,1)      - &
          &    zagc(i,ke_snow+1)/zalfa*t_so_now(i,2)  )
      END IF
    END DO

    ! Snow layers
    DO ksn = 2, ke_snow-1
      DO i = istarts, iends
        IF (sn_frac(i) > 0._ireals) THEN
          zrocs(i) = (wliq_snow_now(i,ksn)/wtot_snow_now(i,ksn)*chc_w + &
            & (wtot_snow_now(i,ksn)-wliq_snow_now(i,ksn))/ &
            & wtot_snow_now(i,ksn)*chc_i)*rho_snow_mult_now(i,ksn)
          zakb = zalas_mult(i,ksn)/zrocs(i)
          zaga(i,ksn) = -zalfa*zdt*zakb/(zdzh_snow(i,ksn)*zdzm_snow(i,ksn))
          zakb = (zalas_mult(i,ksn)/zrocs(i)*zdzm_snow(i,ksn)+zalas_mult(i,ksn+1)/zrocs(i)*zdzm_snow(i,ksn+1))/ &
                 (zdzm_snow(i,ksn)+zdzm_snow(i,ksn+1))
          zagc(i,ksn) = -zalfa*zdt*zakb/(zdzh_snow(i,ksn)*zdzm_snow(i,ksn+1))
          zagb(i,ksn) = 1._ireals - zaga(i,ksn) - zagc(i,ksn)
          zagd(i,ksn) = ztsnow_mult(i,ksn) + &
            &     (1._ireals - zalfa)*( - zaga(i,ksn)/zalfa*ztsnow_mult(i,ksn-1) + &
            &     (zaga(i,ksn)/zalfa + zagc(i,ksn)/zalfa)*ztsnow_mult(i,ksn)     - &
            &     zagc(i,ksn)/zalfa*ztsnow_mult(i,ksn+1)  )
        END IF
      END DO
    END DO                ! snow layers

    ! Soil layers
    DO ksn = ke_snow+2, ke_snow + ke_soil
      DO i = istarts, iends
        IF (sn_frac(i) > 0._ireals) THEN
          kso = ksn - ke_snow
          zakb1 = zalam(i,kso-1)/zroc(i,kso)
          zakb2 = zalam(i,kso  )/zroc(i,kso)
          zaga(i,ksn) = -zalfa*zdt*zakb1/(zdzhs(kso)*zdzms(kso))
          zagc(i,ksn) = -zalfa*zdt*zakb2/(zdzhs(kso)*zdzms(kso+1))
          zagb(i,ksn) = 1._ireals - zaga(i,ksn) - zagc(i,ksn)
          zagd(i,ksn) = t_so_now(i,kso) + &
            &    (1._ireals - zalfa)*( - zaga(i,ksn)/zalfa*t_so_now(i,kso-1) +  &
            &    (zaga(i,ksn)/zalfa + zagc(i,ksn)/zalfa)*t_so_now(i,kso)     -  &
            &    zagc(i,ksn)/zalfa*t_so_now(i,kso+1)  )
        END IF
      END DO
    END DO                ! soil layers

    DO kso = 1, ke_snow + ke_soil
      DO i = istarts, iends
        IF (sn_frac(i) > 0._ireals) THEN
          zzz = 1._ireals/(zagb(i,kso) - zaga(i,kso)*zagc(i,kso-1))
          zagc(i,kso) = zagc(i,kso) * zzz
          zagd(i,kso) = (zagd(i,kso) - zaga(i,kso)*zagd(i,kso-1)) * zzz
        END IF
      END DO
    END DO                ! snow + soil layers

    ! Back substitution, lowermost soil layer
    DO i = istarts, iends
      IF (sn_frac(i) > 0._ireals) THEN
        zage(i,ke_snow+ke_soil+1) =  &
          &    (zagd(i,ke_snow+ke_soil+1) - zaga(i,ke_snow+ke_soil+1) * zagd(i,ke_snow+ke_soil))/ &
          &    (zagb(i,ke_snow+ke_soil+1) - zaga(i,ke_snow+ke_soil+1) * zagc(i,ke_snow+ke_soil))
        t_so_snow_new(i,ke_soil+1) = zage(i,ke_snow+ke_soil+1) ! climate value, unchanged
      END IF
    END DO

    ! Back substitution, soil layers
    DO kso = ke_snow+ke_soil, ke_snow+1, -1
      DO i = istarts, iends
        IF (sn_frac(i) > 0._ireals) THEN
          zage(i,kso) = zagd(i,kso) - zagc(i,kso)*zage(i,kso+1)
          t_so_snow_new(i,kso-ke_snow) = zage(i,kso)
        END IF
      END DO
    END DO                ! soil layers

    ! Back substitution, snow layers
    DO ksn = ke_snow,1,-1
      DO i = istarts, iends
        IF (sn_frac(i) > 0._ireals) THEN
          zage(i,ksn) = zagd(i,ksn) - zagc(i,ksn)*zage(i,ksn+1)
          ztsnown_mult(i,ksn) = zage(i,ksn)
        END IF
      END DO
    END DO                ! snow layers

!in case of thin snowpack (less than zswitch), apply single-layer snow model
    DO i = istarts, iends
      IF (sn_frac(i) > 0._ireals) THEN
        zgsb(i) = ((zalas_mult(i,ke_snow)*(-zhm_snow(i,ke_snow))+hzalam(i,1)*zdzms(1))/ &
                  (-zhm_snow(i,ke_snow)+zdzms(1)) * &
                  (ztsnown_mult(i,ke_snow) - t_so_now(i,1))/(-zhm_snow(i,ke_snow) &
                  +zdzms(1)))*zf_snow(i)
        zrocs(i) = (wliq_snow_now(i,1)/wtot_snow_now(i,1)*chc_w + &
          (wtot_snow_now(i,1)-wliq_snow_now(i,1))/wtot_snow_now(i,1)*chc_i)*rho_snow_mult_now(i,1)
        zswitch(i) = (-zfor_snow_mult(i)+zgsb(i))/50./zrocs(i)*zdt*ke_snow
        zswitch(i) = MAX(zswitch(i),1.E-03_ireals)
      ELSE
        zswitch(i) = 0.0_ireals
      END IF
    END DO
    DO i = istarts, iends
      IF(zwsnew(i) .LT. zswitch(i) .AND. sn_frac(i) > 0._ireals) THEN
        ztsn  (i) = t_so_now(i,1)
        tmp_num(i) = ztsnow_mult(i,1) + zdt*2._ireals*(zfor_snow_mult(i) - zgsb(i))  &
                           /zrocs(i)/(zswitch(i)/rho_snow_mult_now(i,1)*rho_w) &
                           &- ( ztsn(i) - zts(i) )
        zalas  = 2.22_ireals*EXP(1.88_ireals*LOG(rho_snow_mult_now(i,1)/rho_i))

        ztsnow_im    = - zrhoch(i) * (cp_d + zdqvtsnow(i) * lh_s) - zalas/h_snow_now(i)
        zfak  = MAX(zepsi,1.0_ireals - zdt*zalfa*ztsnow_im/zrocs(i)/h_snow_now(i))
        tmp_num(i) = ztsnow_mult(i,1) + (tmp_num(i)-ztsnow_mult(i,1))/zfak
      END IF
    END DO
    DO ksn = 1, ke_snow
      DO i = istarts, iends
        IF(sn_frac(i) > 0._ireals) THEN
          IF(zwsnew(i) .LT. zswitch(i)) THEN
!            ztsnown_mult(i,ksn) = tmp_num(i)
          END IF
          ztsnown_mult(i,0) = ztsnown_mult(i,1)
        END IF
      END DO
    END DO

    DO ksn = 1, ke_snow
        DO i = istarts, iends
            IF(zwsnew(i) .GT. zepsi .and. zwsnew(i) .LT. zswitch(i)) THEN

              IF((zfor_snow_mult(i)-zgsb(i))*zdt > zwsnew(i)*rho_w*lh_f) THEN
                ztsnown_mult(i,ksn) = t_so_now(i,0)
              ELSE IF(zfor_snow_mult(i)-zgsb(i) .GT. 0._ireals) THEN
                ztsnown_mult(i,ksn) = ztsnow_mult(i,ksn) + &
                  (zfor_snow_mult(i)-zgsb(i))*zdt/(chc_i*wtot_snow_now(i,ksn))/rho_w/ke_snow
              END IF
            END IF
        END DO
    END DO

    DO i = istarts, iends
        IF(zwsnew(i) .GT. zepsi .and. zwsnew(i) .LT. zswitch(i)) THEN

          IF((zfor_snow_mult(i)-zgsb(i))*zdt > zwsnew(i)*rho_w*lh_f) THEN
            zdwsndt(i) = zdwsndt(i) - zwsnew(i)*rho_w/zdt
            zwsnew(i)  = 0._ireals
            ztsnown_mult(i,0) = t_so_now(i,0)
          ELSE IF((zfor_snow_mult(i)-zgsb(i)) .GT. 0._ireals) THEN
            ztsnown_mult(i,0) = ztsnown_mult(i,1)
          END IF

        END IF
    END DO

  END IF

  ! Combining the two partial updates of soil temperature
  IF(lmulti_snow) THEN
    DO kso = 1, ke_soil+1
      DO i = istarts, iends
        t_so_new(i,kso) = t_so_snow_new(i,kso)*sn_frac(i) + t_so_free_new(i,kso)*(1._ireals - sn_frac(i))
      END DO
    END DO
  ELSE
    DO kso = 1, ke_soil+1
      DO i = istarts, iends
        t_so_new(i,kso) = t_so_free_new(i,kso)
      END DO
    END DO
  END IF



!  IF(lmelt) THEN ! + lmelt_var
      DO kso = 1,ke_soil
!CDIR NODEP,VOVERTAKE,VOB
!DIR$ IVDEP
         DO ic=1,icount_soil
            i=soil_list(ic)
                ztx      = t0_melt
                zw_m(i)     = zporv(i,kso)*zdzhs(kso)
                IF(t_so_new(i,kso).LT.(t0_melt-zepsi)) THEN
!                  zw_m(i) = zw_m(i)*EXP(-zedb(i,kso)*LOG((t_so_new(i,kso) - t0_melt)/(t_so_new(i,kso)*zaa)) )

                  IF (t_so_new(i,kso) < t_zw_low) THEN
                    zw_m(i) = zw_m_low(i,kso)
                  ELSE IF (t_so_new(i,kso) < t_zw_up) THEN ! Logarithmic Interpolation between -3 degC and -40 degC 
                    zw_m(i) = zw_m_low(i,kso)*EXP((t_so_new(i,kso) - t_zw_low)*                          &
                      (LOG(zporv(i,kso)*zdzhs(kso)*zw_m_up(i)) - LOG(zw_m_low(i,kso)))/(t_zw_up-t_zw_low))
                  ELSE
                    zw_m(i) = zw_m(i)*EXP(-zedb(i)*LOG((t_so_new(i,kso) - t0_melt)/(t_so_new(i,kso)*zaa(i))) )
                  END IF
       
                  zliquid= MAX(zepsi,w_so_now(i,kso) -  w_so_ice_now(i,kso))
                  znen   = 1._ireals-zaa(i)*EXP(zb_por(i)*LOG(zporv(i,kso)*zdzhs(kso)/zliquid))
                  ztx    = t0_melt/znen
                ENDIF
                ztx      = MIN(t0_melt,ztx)
                zfak     = zroc(i,kso)*zdzhs(kso)/(lh_f*rho_w)
                zdelwice = - zfak*(t_so_new(i,kso)-ztx)
                zwso_new  = w_so_now(i,kso) + zdt*zdwgdt(i,kso)/rho_w
                zargu = zwso_new - zw_m(i) - w_so_ice_now(i,kso)
                IF (t_so_new(i,kso) > t0_melt .AND. w_so_ice_now(i,kso) > 0._ireals) THEN
                  ! melting point adjustment (time scale 30 min)
                  zdelwice = - MIN(w_so_ice_now(i,kso), zdwi_scal*(t_so_new(i,kso)-t0_melt)*zfak)
                ELSE IF (zdelwice < 0.0_ireals) THEN
                  zdelwice = - MIN( - zdelwice,MIN(-zargu,w_so_ice_now(i,kso)))
                  ! limit latent heat consumption due to melting to half the temperature increase since last time step
                  ! or 2.5 K within 30 min
                  zdelwice = - MIN( - zdelwice,MAX(2.5_ireals*zdwi_scal,0.5_ireals*(t_so_new(i,kso)-t_so_now(i,kso)))*zfak)
                ELSE
                  zdelwice = MIN(zdelwice,MAX(zargu,0.0_ireals))
                  ! limit latent heat release due to freezing to half the differene from the melting point
                  zdelwice = MIN(zdelwice,0.5_ireals*(t0_melt-t_so_new(i,kso))*zfak)
                ENDIF
                w_so_ice_new(i,kso) = w_so_ice_now(i,kso) + zdelwice
                t_so_new(i,kso) = t_so_new(i,kso) + zdelwice/zfak
          END DO
      ENDDO
!  END IF ! lmelt


!------------------------------------------------------------------------------
! Section II.7: Energy budget and temperature prediction at snow-surface
!------------------------------------------------------------------------------

  DO i = istarts, iends
!      IF (llandmask(i)) THEN          ! land-points only
        ! next line has to be changed if a soil surface temperature is
        ! predicted by the heat conduction equation
        zdtsdt (i) = (t_so_new(i,1) - zts(i))*z1d2dt
        ztsn   (i) =  t_so_new(i,1)
        IF(.NOT. lmulti_snow)  &
          ztsnown(i) = ztsn(i)       ! default setting
        zwsnew(i)       = zwsnow(i) + zdwsndt(i)*zdtdrhw

        ! forcing contributions for snow formation of dew and rime are
        ! contained in ze_ges, heat fluxes must not be multiplied by
        ! snow covered fraction

        IF(.NOT. lmulti_snow) THEN

          zrnet_snow    = sobs(i) + zthsnw(i)
          zshfl_snow(i) = zrhoch(i)*cp_d*(zth_low(i) - ztsnow(i))
          zlhfl_snow(i) = lh_s*zversn(i)
!DR start
          zqhfl_snow(i) = zversn(i)
!DR end
          zfor_snow     = zrnet_snow + zshfl_snow(i) + zlhfl_snow(i)

          ! forecast of snow temperature Tsnow
          IF (ztsnow(i) < t0_melt .AND. zwsnew(i) > zepsi) THEN
            ztsnown(i) = ztsnow(i) + zdt*2._ireals*(zfor_snow - zgsb(i))  &
                           /zrocs(i) - ( ztsn(i) - zts(i) )

            ! implicit formulation
! BR        zalas  = MAX(calasmin,MIN(calasmax, calasmin + calas_dw*zwsnew(i)))
! BR 7/2005 Introduce new dependency of snow heat conductivity on snow density
!
            zalas  = 2.22_ireals*EXP(1.88_ireals*LOG(zrho_snow(i)/rho_i))

            ztsnow_im    = - zrhoch(i) * (cp_d + zdqvtsnow(i) * lh_s)       &
                                         - zalas/zdz_snow_fl(i)
            zfak  = MAX(zepsi,1.0_ireals - zdt*zalfa*ztsnow_im/zrocs(i))
            ztsnown(i) = ztsnow(i) + (ztsnown(i)-ztsnow(i))/zfak
          END IF

          zdtsnowdt(i) = (ztsnown(i) - ztsnow(i))*z1d2dt
        ENDIF
!      END IF          ! land-points only
  END DO

IF (msg_level >= 19) THEN
  DO i = istarts, iends
  IF (soiltyp_subs(i) == 1) THEN  !1=glacier and Greenland
    IF ( ABS( zshfl_snow(i) )  >  500.0_ireals  .OR. &
         ABS( zlhfl_snow(i) )  > 2000.0_ireals ) THEN
      write(*,*) 'soil: ', zshfl_snow(i), zlhfl_snow(i), '....', &
        zth_low(i), ztsnow(i), '....', &
        zwsnow(i), zrr(i), zrs(i), zdwsndt(i)
    ENDIF
  ENDIF
  END DO
ENDIF

  IF (lmulti_snow) THEN
    DO ksn = 0, ke_snow
        DO i = istarts, iends
!          IF (llandmask(i)) THEN          ! land-points only
            IF (zwsnew(i) > zepsi) THEN
              zdtsnowdt_mult(i,ksn) = (ztsnown_mult(i,ksn) - ztsnow_mult(i,ksn))*z1d2dt
            END IF
!          ENDIF
        END DO
    ENDDO
  ENDIF


!------------------------------------------------------------------------------
! Section II.8: Melting of snow ,infiltration and surface runoff of snow water
!------------------------------------------------------------------------------

! If the soil surface temperature predicted by the equation of heat conduction
! is used instead of using T_s = T_so(1), the following section has to be
! adjusted accordingly.
!
! Basically this snow model uses heat fluxes to either heat the uppermost soil
! layer, if the snow surface temperature exceeds t0_melt, or to heat the snow, if
! the temperature of the uppermost soil layer exceeds t0_melt. Melting is considered
! after this process, if the snow temperature equals t0_melt AND the temperature of
! the uppermost soil layer exceeds t0_melt. The excess heat (t_so(1)-t0_melt) is used for
! melting snow. In cases this melting may be postponed to the next time step.

  IF (.NOT. lmulti_snow) THEN

      DO i = istarts, iends
!        IF (llandmask(i)) THEN          ! land-points only

          zwsnn  (i)  = zwsnow(i) + zdtdrhw*zdwsndt(i)
          zwsnew (i)  = zwsnn(i)
          ztsnownew     = ztsnown(i)
          ze_avail      = 0.0_ireals
          ze_total      = 0.0_ireals
          zfr_melt      = 0.0_ireals

          IF (zwsnew(i) > zepsi) THEN        ! points with snow cover only
            ! first case: T_snow > t0_melt: melting from above
            ! ----------
            IF (ztsnown(i) > t0_melt .AND. t_so_new(i,1) < t0_melt ) THEN
              ! Limit max w_snow in melting conditions for consistency with heat capacity calculation
              zdwsnm(i)    = MIN(1.5_ireals*rho_snow_now(i)/rho_w,zwsnew(i))* &
                    .5_ireals*(ztsnown(i) - (t0_melt - zepsi))/ &
                   (.5_ireals* (zts(i) - (t0_melt - zepsi)) - lh_f/chc_i)
              zdwsnm(i)    = zdwsnm(i)*z1d2dt*rho_w
              zdwsndt(i)   = zdwsndt (i) + zdwsnm(i)
              meltrate(i)  = - zdwsnm(i)
              ztsnownew    = t0_melt - zepsi
              zdtsnowdt(i) = zdtsnowdt(i) + (ztsnownew - ztsnown(i))*z1d2dt
              runoff_s (i) = runoff_s(i) - zdwsnm(i)*zroffdt
            ENDIF ! melting from above

            IF (t_so_new(i,1) >= t0_melt) THEN
              !second case:  temperature of uppermost soil layer > t0_melt. First a
              !-----------   heat redistribution is performed. As a second step,
              !              melting of snow is considered.
              ! a) Heat redistribution
              ztsnew = t0_melt + zepsi
              ztsnownew      = ztsnown(i) + zf_snow_lim(i)*(ztsn(i) - ztsnew) +  &
                   2._ireals*(ztsn(i) - ztsnew)*zroc(i,1)*zdzhs(1)/zrocs(i)
              zdtsdt(i)    = zdtsdt(i) + zf_snow_lim(i)*(ztsnew - t_so_new(i,1))*z1d2dt
              zdtsnowdt(i) = zdtsnowdt(i) + (ztsnownew - ztsnown(i))*z1d2dt
              ! b) Melting of snow (if possible)
              IF (ztsnownew > t0_melt) THEN
                ze_avail     = 0.5_ireals*(ztsnownew - t0_melt)*zrocs(i)*zf_snow_lim(i)
                ze_total     = lh_f*zwsnew(i)*rho_w
                zfr_melt     = MIN(1.0_ireals,ze_avail/ze_total)
                zdtsnowdt(i)= zdtsnowdt(i) + (t0_melt - ztsnownew)*z1d2dt
                zdelt_s      = MAX(0.0_ireals,(ze_avail - ze_total)/(zroc(i,1)* &
                                                                        zdzhs(1)))
                zdtsdt(i)  = zdtsdt(i) + zf_snow_lim(i)*zdelt_s*z1d2dt

                ! melted snow is allowed to penetrate the soil (up to field
                ! capacity), if the soil type is neither ice nor rock (zrock = 0);
                ! else it contributes to surface run-off;
                ! fractional water content of the first soil layer determines
                ! a reduction factor which controls additional run-off
                zdwsnm(i)   = zfr_melt*zwsnew(i)*z1d2dt*rho_w  ! available water
                zdwsndt(i)  = zdwsndt (i) - zdwsnm(i)
                meltrate(i) = meltrate(i) + zdwsnm(i)
                zdwgme        = zdwsnm(i)*zrock(i)             ! contribution to w_so
                zro           = (1._ireals - zrock(i))*zdwsnm(i)      ! surface runoff
                zredfu        = MAX( 0.0_ireals,  MIN( 1.0_ireals, (zw_fr(i,1) -  &
                                zfcap(i,1))/MAX(zporv(i,1)-zfcap(i,1), zepsi)))
                zdwgdt(i,1) = zdwgdt(i,1) + zdwgme*(1._ireals - zredfu)
                zro           = zro + zdwgme*zredfu    ! Infiltration not possible
                                                       ! for this fraction

                ! zro-, zdw_so_dt-correction in case of pore volume overshooting
                zw_ovpv = MAX(0._ireals, zw_fr(i,1)* zdzhs(1) * zrhwddt +  &
                           zdwgdt(i,1) - zporv(i,1) * zdzhs(1) * zrhwddt)
                zro = zro + zw_ovpv
                zdwgdt(i,1)= zdwgdt(i,1) - zw_ovpv


                IF (zfr_melt > 0.9999_ireals) zdwsndt(i)= -zwsnow(i)*zrhwddt
                runoff_s (i)= runoff_s(i) + zro*zroffdt

              END IF   ! snow melting
            END IF     ! snow and/or soil temperatures
          END IF       ! points with snow cover only
!        END IF         ! land-points only
       END DO

  ELSE       ! new snow scheme

    DO i = istarts, iends
!      IF (llandmask(i)) THEN          ! land-points only

        zwsnew(i) = zwsnow(i) + zdtdrhw*zdwsndt(i)
        zdwsnm(i) = 0.0_ireals

        ze_out  (i) = 0.0_ireals
        zqbase  (i) = 0.0_ireals
        zcounter(i) = 0.0_ireals

        ze_rad(i) = 0.0_ireals
        IF(zextinct(i,1).gt.0.0_ireals) ze_rad(i) = zf_snow(i) * sobs(i)

        ztsnownew_mult(i,0) = ztsnown_mult(i,0)
!      END IF         ! land-points only
    END DO

    DO ksn = 1,ke_snow

        DO i = istarts, iends
!          IF (llandmask(i)) THEN          ! land-points only
            IF (zwsnew(i) > zepsi) THEN        ! points with snow cover only

              zrefr(i) = 0.0_ireals
              zmelt(i) = 0.0_ireals
              ztsnownew_mult(i,ksn) = ztsnown_mult(i,ksn)

              IF(zdzh_snow(i,ksn) - wliq_snow_now(i,ksn).GT.zepsi .OR. &
                wtot_snow_now(i,ksn) - wliq_snow_now(i,ksn).GT.zepsi) THEN
                zrho_dry_old(i) = MAX(wtot_snow_now(i,ksn)-wliq_snow_now(i,ksn), &
                  &                     zepsi)                                               &
                  &                 *rho_w/(zdzh_snow(i,ksn) - wliq_snow_now(i,ksn))
              ELSE
                zrho_dry_old(i) = rho_w
              END IF

              ztsnownew_mult(i,ksn) = (ztsnown_mult(i,ksn)*wtot_snow_now(i,ksn) &
                &                       + t0_melt*zqbase(i)*zdt)/(zqbase(i)*zdt     &
                &                       + wtot_snow_now(i,ksn))

              IF(zextinct(i,ksn).eq.0.0_ireals) THEN
                ze_in = ze_out(i)
              ELSE
                IF(ksn.eq.ke_snow) THEN     ! all the rest of radiation is absorbed by the lowermost snow layer
                  ze_in = ze_out(i) + ze_rad(i)
                ELSE
                  zcounter(i) = EXP (-zextinct(i,ksn)*zdzh_snow(i,ksn))
                  ze_in = ze_out(i) + ze_rad(i) * (1._ireals - zcounter(i))
                  ze_rad(i) = ze_rad(i) * zcounter(i)
                END IF
              END IF

              ztsnownew_mult(i,ksn) = ztsnownew_mult(i,ksn) &
                &                       + ze_in*zdt/(chc_i*wtot_snow_now(i,ksn))/rho_w
              wtot_snow_now(i,ksn) = wtot_snow_now(i,ksn) + zqbase(i)*zdt
              wliq_snow_now(i,ksn) = wliq_snow_now(i,ksn) + zqbase(i)*zdt

              zdzh_old = zdzh_snow(i,ksn)
              zdzh_snow(i,ksn) = zdzh_snow(i,ksn) + zqbase(i)*zdt

              rho_snow_mult_now(i,ksn) = MAX(wtot_snow_now(i,ksn)*&
                &                                rho_w/zdzh_snow(i,ksn), &
                &                                0.0_ireals)

              IF(ztsnownew_mult(i,ksn) .GT. t0_melt) THEN

                IF(wtot_snow_now(i,ksn) .LE. wliq_snow_now(i,ksn)) THEN
                  ze_out(i) = chc_i*wtot_snow_now(i,ksn)*(ztsnownew_mult(i,ksn)-t0_melt) &
                    &      *z1d2dt*rho_w
                  zmelt(i) = 0.0_ireals
                ELSEIF(chc_i*wtot_snow_now(i,ksn)*(ztsnownew_mult(i,ksn)-t0_melt)/lh_f <= &
                  wtot_snow_now(i,ksn)-wliq_snow_now(i,ksn)) THEN
                  zmelt(i) = chc_i*wtot_snow_now(i,ksn)*(ztsnownew_mult(i,ksn)-t0_melt) &
                    &          *z1d2dt/lh_f
                  ze_out(i) = 0.0_ireals
                  wliq_snow_now(i,ksn) = wliq_snow_now(i,ksn) + zmelt(i)*zdt
                ELSE
                  zmelt(i) = (wtot_snow_now(i,ksn)-wliq_snow_now(i,ksn))*z1d2dt
                  ze_out(i) = chc_i*wtot_snow_now(i,ksn)*(ztsnownew_mult(i,ksn)-t0_melt) &
                    &      *z1d2dt*rho_w - zmelt(i)*lh_f*rho_w
                  wliq_snow_now(i,ksn) = wliq_snow_now(i,ksn) + zmelt(i)*zdt
                END IF
                ztsnownew_mult(i,ksn) = t0_melt

              ELSE
!          T<0
                IF(wliq_snow_now(i,ksn) .GT. -chc_i*wtot_snow_now(i,ksn) &
                  & *(ztsnownew_mult(i,ksn) - t0_melt)/lh_f) THEN
                  zrefr(i) = -chc_i*wtot_snow_now(i,ksn)*(ztsnownew_mult(i,ksn) &
                    &          - t0_melt)*z1d2dt/lh_f
                  ztsnownew_mult(i,ksn)   = t0_melt
                  wliq_snow_now(i,ksn) = wliq_snow_now(i,ksn) - zrefr(i)*zdt
                ELSE
                  zrefr(i) = wliq_snow_now(i,ksn)*z1d2dt
                  wliq_snow_now(i,ksn) = 0.0_ireals
                  ztsnownew_mult(i,ksn)   = ztsnownew_mult(i,ksn) + zrefr(i)*zdt*lh_f &
                    &                         /(chc_i*wtot_snow_now(i,ksn))
                END IF
                ze_out(i) = 0.0_ireals

              END IF

              zdtsnowdt_mult(i,ksn) = zdtsnowdt_mult(i,ksn) + &
                                        (ztsnownew_mult(i,ksn) - ztsnown_mult(i,ksn))*z1d2dt
              IF(wtot_snow_now(i,ksn) .LE. wliq_snow_now(i,ksn)) THEN
                zqbase(i)           = wliq_snow_now(i,ksn)*z1d2dt
                wliq_snow_now(i,ksn) = 0.0_ireals
                wtot_snow_now(i,ksn) = 0.0_ireals
                zdzh_snow(i,ksn)    = 0.0_ireals
                rho_snow_mult_now(i,ksn)  = 0.0_ireals
              ELSE
                IF(zrefr(i) .GT. 0.0_ireals .OR. zmelt(i) .GT. 0.0_ireals) THEN
                  zadd_dz = 0.0_ireals
                  zadd_dz = MAX(zrefr(i),0._ireals)*(-1.0_ireals + 1.0_ireals/rho_i*rho_w)*zdt
                  zadd_dz = MAX(zmelt(i),0._ireals)*(-1.0_ireals/zrho_dry_old(i)*rho_w &
                    &       + 1.0_ireals)*zdt
                  zdzh_snow(i,ksn)   = zdzh_snow(i,ksn) + zadd_dz
                  rho_snow_mult_now(i,ksn) = MAX(wtot_snow_now(i,ksn)*rho_w &
                    &                            /zdzh_snow(i,ksn),0.0_ireals)
                  IF(wtot_snow_now(i,ksn) .LE. 0.0_ireals) zdzh_snow(i,ksn) = 0.0_ireals
                  IF(rho_snow_mult_now(i,ksn) .GT. rho_w) THEN
                    zdzh_snow(i,ksn)   = zdzh_snow(i,ksn)*rho_snow_mult_now(i,ksn)/rho_w
                    rho_snow_mult_now(i,ksn) = rho_w
                  END IF
                END IF

                zsn_porosity = 1._ireals - (rho_snow_mult_now(i,ksn)/rho_w -  &
                               wliq_snow_now(i,ksn)/zdzh_snow(i,ksn))/rho_i*rho_w - &
                               wliq_snow_now(i,ksn)/zdzh_snow(i,ksn)
                zsn_porosity = MAX(zsn_porosity,cwhc + 0.1_ireals)
                zp1 = zsn_porosity - cwhc

                IF (wliq_snow_now(i,ksn)/zdzh_snow(i,ksn) .GT. cwhc) THEN
                  zfukt             = (wliq_snow_now(i,ksn)/zdzh_snow(i,ksn) - cwhc)/zp1
                  zq0               = chcond * zfukt**3
                  zqbase(i)       = MIN(zq0*zdt,wliq_snow_now(i,ksn))
                  wliq_snow_now(i,ksn) = wliq_snow_now(i,ksn) - zqbase(i)
                  wtot_snow_now(i,ksn) = wtot_snow_now(i,ksn) - zqbase(i)

                  zdzh_old = zdzh_snow(i,ksn)
                  zdzh_snow(i,ksn) = zdzh_snow(i,ksn) - zqbase(i)
                  zqbase(i)        = zqbase(i)*z1d2dt

                  IF(zdzh_snow(i,ksn) .LT. zepsi*0.01_ireals) THEN
                    wliq_snow_now(i,ksn) = 0.0_ireals
                    wtot_snow_now(i,ksn) = 0.0_ireals
                    zdzh_snow(i,ksn)     = 0.0_ireals
                    rho_snow_mult_now(i,ksn)  = 0.0_ireals
                  ELSE
                    rho_snow_mult_now(i,ksn) = MAX(wtot_snow_now(i,ksn)*rho_w &
                      &                            /zdzh_snow(i,ksn),0.0_ireals)
                    IF(wtot_snow_now(i,ksn) .LE. 0.0_ireals) zdzh_snow(i,ksn) = 0.0_ireals
                    IF(rho_snow_mult_now(i,ksn) .GT. rho_w) THEN
                      zdzh_snow(i,ksn)   = zdzh_snow(i,ksn)*rho_snow_mult_now(i,ksn)/rho_w
                      rho_snow_mult_now(i,ksn) = rho_w
                    END IF
                  END IF
                ELSE
                  zqbase(i) = 0.0_ireals
                END IF
              END IF

            END IF       ! points with snow cover only
!          END IF         ! land-points only
        END DO
    END DO        ! snow layers

    DO i = istarts, iends
!      IF (llandmask(i)) THEN          ! land-points only
        IF (zwsnew(i) > zepsi) THEN        ! points with snow cover only
          zdwsnm(i) = zqbase(i)*rho_w       ! ksn == ke_snow
        END IF       ! points with snow cover only
!      END IF         ! land-points only
    END DO

    DO i = istarts, iends
!       IF (llandmask(i)) THEN          ! land-points only
         IF (zwsnew(i) > zepsi) THEN        ! points with snow cover only
           zdwsndt(i)  = zdwsndt(i) - zdwsnm(i)

             ! melted snow is allowed to penetrate the soil (up to field
             ! capacity), if the soil type is neither ice nor rock (zrock = 0);
             ! else it contributes to surface run-off;
             ! fractional water content of the first soil layer determines
             ! a reduction factor which controls additional run-off

           zdwgme        = zdwsnm(i)*zrock(i)             ! contribution to w_so
           zro           = (1._ireals - zrock(i))*zdwsnm(i)      ! surface runoff
           zredfu        = MAX( 0.0_ireals,  MIN( 1.0_ireals, (zw_fr(i,1) -  &
                           zfcap(i,1))/MAX(zporv(i,1)-zfcap(i,1), zepsi)))
           zdwgdt(i,1) = zdwgdt(i,1) + zdwgme*(1._ireals - zredfu)
           zro           = zro + zdwgme*zredfu    ! Infiltration not possible
                                                      ! for this fraction

           ! zro-, zdw_so_dt-correction in case of pore volume overshooting
           zw_ovpv = MAX(0._ireals, zw_fr(i,1)* zdzhs(1) * zrhwddt +  &
                     zdwgdt(i,1) - zporv(i,1) * zdzhs(1) * zrhwddt)
           zro = zro + zw_ovpv
           zdwgdt(i,1)= zdwgdt(i,1) - zw_ovpv

           runoff_s(i) = runoff_s(i) + zro*zroffdt
         END IF       ! points with snow cover only
!       END IF         ! land-points only
    END DO

! snow densification due to gravity and metamorphism

    DO ksn = 2, ke_snow
      zp(:,ksn) = 0.0_ireals                         ! gravity, Pa
      DO k = ksn,1,-1
          DO i = istarts, iends
!            IF (llandmask(i)) THEN          ! land-points only
              zp(i,ksn) = zp(i,ksn) + rho_snow_mult_now(i,k)*g*zdzh_snow(i,ksn)
!            END IF         ! land-points only
          END DO
      END DO
    END DO

    DO ksn = 2, ke_snow
        DO i = istarts, iends
!          IF (llandmask(i)) THEN          ! land-points only
            IF (zwsnew(i) > zepsi) THEN        ! points with snow cover only
              IF(rho_snow_mult_now(i,ksn) .LT. 600._ireals .AND. &
                rho_snow_mult_now(i,ksn) .NE. 0.0_ireals) THEN
                zdens_old = rho_snow_mult_now(i,ksn)
                zeta =         &! compactive viscosity of snow
                  ca2*EXP(19.3_ireals*rho_snow_mult_now(i,ksn)/rho_i)* &
                  EXP(67300._ireals/8.31_ireals/ztsnownew_mult(i,ksn))
                rho_snow_mult_now(i,ksn) = rho_snow_mult_now(i,ksn) + &
                  zdt*rho_snow_mult_now(i,ksn)*(csigma+zp(i,ksn))/zeta
                rho_snow_mult_now(i,ksn) = MIN(rho_snow_mult_now(i,ksn),rho_i)
                zdzh_snow(i,ksn)   = zdzh_snow(i,ksn) * zdens_old/rho_snow_mult_now(i,ksn)
              END IF
            END IF       ! points with snow cover only
!          END IF         ! land-points only
        END DO
    END DO

    DO i = istarts, iends
!      IF (llandmask(i)) THEN          ! land-points only
        IF (zwsnew(i) > zepsi) THEN        ! points with snow cover only

          IF(ztsnownew_mult(i,0) .GT. t0_melt) THEN
            ztsnownew_mult(i,0) = t0_melt
            zdtsnowdt_mult(i,0) = zdtsnowdt_mult(i,0) +     &
                                    (ztsnownew_mult(i,0) - ztsnown_mult(i,0))*z1d2dt
          END IF
        END IF       ! points with snow cover only
!      END IF         ! land-points only
    END DO

  END IF

!------------------------------------------------------------------------------
! Section II.9: Final updating of prognostic values
!------------------------------------------------------------------------------

  IF (lmulti_snow) THEN
    ! First for ksn == 0
    DO i = istarts, iends
!      IF (llandmask(i)) THEN  ! for landpoints only
        t_snow_mult_new  (i,0) = t_snow_mult_now(i,0) + zdt*zdtsnowdt_mult(i,0)
        t_snow_new(i) = t_snow_mult_new (i,0)
!      ENDIF
    ENDDO

    DO ksn = 1, ke_snow
        DO i = istarts, iends
!          IF (llandmask(i)) THEN  ! for landpoints only
            t_snow_mult_new  (i,ksn) = t_snow_mult_now(i,ksn) + &
              &                              zdt*zdtsnowdt_mult(i,ksn)
            dzh_snow_new     (i,ksn) = zdzh_snow(i,ksn)
            wtot_snow_new    (i,ksn) = wtot_snow_now(i,ksn)
            rho_snow_mult_new(i,ksn) = rho_snow_mult_now(i,ksn)
            wliq_snow_new    (i,ksn) = wliq_snow_now(i,ksn)
!          ENDIF
        ENDDO
    ENDDO
  ELSE
    DO i = istarts, iends
!      IF (llandmask(i)) THEN  ! for landpoints only
        t_snow_new(i)  = t_snow_now(i) + zdt*zdtsnowdt(i)
!      ENDIF
    ENDDO
  ENDIF

  DO i = istarts, iends
!      IF (llandmask(i)) THEN  ! for landpoints only
        ! t_snow is computed above
        ! t_snow(i,nnew)  = t_snow(i,nx) + zdt*zdtsnowdt(i)
        t_so_new(i,1)  = t_so_now(i,1) + zdt*zdtsdt   (i)         ! (*)

        ! Next line has to be changed, if the soil surface temperature
        ! t_so(i,0,nnew) predicted by the heat conduction equation is used
        t_s_new   (i)    = t_so_new(i,1)
        t_so_new  (i,0)  = t_so_new(i,1)
        w_snow_new(i)  = w_snow_now(i) + zdt*zdwsndt  (i)/rho_w
 IF (itype_interception == 1) THEN
        w_i_new   (i)  = w_i_now(i) + zdt*zdwidt   (i)/rho_w
 ELSE IF (itype_interception == 2) THEN
      w_i_new (i)  = w_i_now(i) + zdt*zdwidt   (i)/rho_w
      w_p_new (i)  = zwpnstr(i)
      w_s_new (i)  = zwisnstr(i)
 END IF


!      END IF          ! land-points only
     END DO

  ! Reset t_snow_new to t_so(0) if no snow was present at the beginning of the time step
  ! The heat balance calculation is incomplete in this case and sometimes yields unreasonable results
  DO i = istarts, iends
    IF (w_snow_now(i) < zepsi .AND. w_snow_new(i) >= zepsi) THEN
      t_snow_new(i) = MIN(t0_melt,t_so_new(i,0))
      IF (lmulti_snow) THEN
        t_snow_mult_new(i,:) = t_snow_new(i)
      ENDIF
    ENDIF
  ENDDO

!>JH New solution of heat conduction for snow points which melted completly
!    during time step
!------------------------------------------------------------------------------
! Section II.6n: Solution of the heat conduction equation, freezing/melting
!               of soil water/ice (optionally)
!------------------------------------------------------------------------------

! Index list of completely melting snow points
     icount_snow=0
     DO i = istarts, iends
        IF (w_snow_now(i) > zepsi .AND. w_snow_new(i) < zepsi) THEN ! Snow vanished during time step
         icount_snow=icount_snow+1
         melt_list(icount_snow)=i
         zfor_s(i)=0._ireals ! no soil forcing is needed at this step
                             ! only distribution of heat
       END IF
     END DO

!     New update of soil heat capacity
      DO   kso = 1,ke_soil+1
!CDIR NODEP,VOVERTAKE,VOB
        DO ic=1,icount_snow
          i=melt_list(ic)
          ziw_fr(i,kso) = w_so_ice_new(i,kso)/zdzhs(kso)   ! ice frac.
          zlw_fr(i,kso) = w_so_new(i,kso)/zdzhs(kso) - ziw_fr(i,kso)  ! liquid water frac.
          zroc(i,kso)   = zrocg(i,kso) + rho_w*zlw_fr(i,kso)*chc_w +          &
                                         rho_w*ziw_fr(i,kso)*chc_i
       END DO      !soil layers
      END DO


      DO kso = 2,ke_soil
!CDIR NODEP,VOVERTAKE,VOB
        DO ic=1,icount_snow
          i=melt_list(ic)
!        IF (llandmask(i)) THEN          ! land-points only
          ! for heat conductivity: zalam is now 3D
          zakb1 = zalam(i,kso-1)/zroc(i,kso)
          zakb2 = zalam(i,kso  )/zroc(i,kso)
          zaga(i,kso) = -zalfa*zdt*zakb1/(zdzhs(kso)*zdzms(kso))
          zagc(i,kso) = -zalfa*zdt*zakb2/(zdzhs(kso)*zdzms(kso+1))
          zagb(i,kso) = 1._ireals - zaga(i,kso) - zagc(i,kso)
          zagd(i,kso) = t_so_new(i,kso) +                                     &    ! distribute heat in (*)
                 (1._ireals - zalfa)*( - zaga(i,kso)/zalfa*t_so_new(i,kso-1)+ &
                 (zaga(i,kso)/zalfa + zagc(i,kso)/zalfa)*t_so_new(i,kso) -  &
                  zagc(i,kso)/zalfa*t_so_new(i,kso+1)  )
!        END IF  ! land-points only
         END DO        ! soil layers
       END DO


!CDIR NODEP,VOVERTAKE,VOB
     DO ic=1,icount_snow
        i=melt_list(ic)
!    IF (llandmask(i)) THEN          ! land-points only
      ! for heat conductivity: zalam is now 3D: here we need layer 1
      zakb1 = hzalam(i,1)/zroc(i,1)
      zakb2 =  zalam(i,1)/zroc(i,1)
      zaga(i,  1) = -zalfa*zdt*zakb1/(zdzhs(1)*zdzms(1))
      zagc(i,  1) = -zalfa*zdt*zakb2/(zdzhs(1)*zdzms(2))
      zagb(i,  1) = 1._ireals - zaga(i,1) - zagc(i,1)
      zagd(i,  1) = t_so_new(i,1) + (1._ireals - zalfa)* (                 &
                      - zaga(i,1)/zalfa * t_s_new(i) +                     &
                      (zaga(i,1) + zagc(i,1))/zalfa * t_so_new(i,1) -    &
                       zagc(i,1)/zalfa * t_so_new(i,2)   )
      zaga(i,0)    = 0.0_ireals
      zagb(i,0)    = zalfa
      zagc(i,0)    = -zalfa
      zagd(i,0)    = zdzms(1) * zfor_s(i)/hzalam(i,1)+(1._ireals-zalfa)* &
                      (t_so_new(i,1) - t_s_new(i))
      zaga(i,ke_soil+1) = 0.0_ireals
      zagb(i,ke_soil+1) = 1.0_ireals
      zagc(i,ke_soil+1) = 0.0_ireals
      zagd(i,ke_soil+1) = t_so_new(i,ke_soil+1)

!    END IF          ! land-points only
  END DO

!CDIR NODEP,VOVERTAKE,VOB
     DO ic=1,icount_snow
        i=melt_list(ic)
!    IF (llandmask(i)) THEN          ! land-points only
      zagc(i,0) = zagc(i,0)/zagb(i,0)
      zagd(i,0) = zagd(i,0)/zagb(i,0)
!    END IF          ! land-points only
   END DO


     DO kso=1,ke_soil
!CDIR NODEP,VOVERTAKE,VOB
       DO ic=1,icount_snow
         i=melt_list(ic)

!        IF (llandmask(i)) THEN          ! land-points only
         zzz = 1._ireals/(zagb(i,kso) - zaga(i,kso)*zagc(i,kso-1))
         zagc(i,kso) = zagc(i,kso) * zzz
         zagd(i,kso) = (zagd(i,kso) - zaga(i,kso)*zagd(i,kso-1)) * zzz
!        END IF          ! land-points only
       END DO                ! soil layers
     END DO


     DO ic=1,icount_snow
        i=melt_list(ic)
!    IF (llandmask(i)) THEN          ! land-points only
      zage(i,ke_soil+1) = (zagd(i,ke_soil+1) - zaga(i,ke_soil+1)*       &
                             zagd(i,ke_soil))/                              &
                            (zagb(i,ke_soil+1) - zaga(i,ke_soil+1)*       &
                             zagc(i,ke_soil))
!    END IF          ! land-points only
    END DO


    DO kso = ke_soil,0,-1
!CDIR NODEP,VOVERTAKE,VOB
      DO ic=1,icount_snow
        i=melt_list(ic)
!        IF (llandmask(i)) THEN          ! land-points only
          zage(i,kso)     = zagd(i,kso) - zagc(i,kso)*zage(i,kso+1)
        ! The surface temperature computed by t_so(i,0,nnew)=zage(i,0) is
        ! presently unused
          t_so_new(i,kso) = zage(i,kso)
!        END IF          ! land-points only
      END DO                ! soil layers
    END DO


     DO ic=1,icount_snow
       i=melt_list(ic)
!    IF (llandmask(i)) THEN          ! land-points only
      t_so_new(i,ke_soil+1) = zage(i,ke_soil+1) ! climate value, unchanged
!    END IF          ! land-points only
   END DO



!  IF(lmelt) THEN ! + melt_var!
     DO kso = 1,ke_soil
!CDIR NODEP,VOVERTAKE,VOB
       DO ic=1,icount_snow
         i=melt_list(ic)
!            IF (llandmask(i)) THEN ! land points only,
              IF (m_styp(i).ge.3) THEN ! neither ice or rocks
                ztx      = t0_melt
                zw_m(i)     = zporv(i,kso)*zdzhs(kso)
                IF(t_so_new(i,kso).LT.(t0_melt-zepsi)) THEN
!                  zw_m(i) = zw_m(i)*EXP(-zedb(i,kso)*LOG((t_so_new(i,kso) - t0_melt)/(t_so_new(i,kso)*zaa)) )

                  IF (t_so_new(i,kso) < t_zw_low) THEN
                    zw_m(i) = zw_m_low(i,kso)
                  ELSE IF (t_so_new(i,kso) < t_zw_up) THEN ! Logarithmic Interpolation between -3 degC and -40 degC 
                    zw_m(i) = zw_m_low(i,kso)*EXP((t_so_new(i,kso) - t_zw_low)*                          &
                      (LOG(zporv(i,kso)*zdzhs(kso)*zw_m_up(i)) - LOG(zw_m_low(i,kso)))/(t_zw_up-t_zw_low))
                  ELSE
                    zw_m(i) = zw_m(i)*EXP(-zedb(i)*LOG((t_so_new(i,kso) - t0_melt)/(t_so_new(i,kso)*zaa(i))) )
                  END IF

                  zliquid= MAX(zepsi,w_so_now(i,kso) -  w_so_ice_now(i,kso))
                  znen   = 1._ireals-zaa(i)*EXP(zb_por(i)*LOG(zporv(i,kso)*zdzhs(kso)/zliquid))
                  ztx    = t0_melt/znen
                ENDIF
                ztx      = MIN(t0_melt,ztx)
                zfak     = zroc(i,kso)*zdzhs(kso)/(lh_f*rho_w)
                zdelwice = - zfak*(t_so_new(i,kso)-ztx)
                zwso_new  = w_so_now(i,kso) + zdt*zdwgdt(i,kso)/rho_w
                zargu = zwso_new - zw_m(i) - w_so_ice_now(i,kso)
                IF (t_so_new(i,kso) > t0_melt .AND. w_so_ice_now(i,kso) > 0._ireals) THEN
                  ! melting point adjustment (time scale 30 min)
                  zdelwice = - MIN(w_so_ice_now(i,kso), zdwi_scal*(t_so_new(i,kso)-t0_melt)*zfak)
                ELSE IF (zdelwice < 0.0_ireals) THEN
                  zdelwice = - MIN( - zdelwice,MIN(-zargu,w_so_ice_now(i,kso)))
                  ! limit latent heat consumption due to melting to half the temperature increase since last time step
                  ! or 2.5 K within 30 min
                  zdelwice = - MIN( - zdelwice,MAX(2.5_ireals*zdwi_scal,0.5_ireals*(t_so_new(i,kso)-t_so_now(i,kso)))*zfak)
                ELSE
                  zdelwice = MIN(zdelwice,MAX(zargu,0.0_ireals))
                  ! limit latent heat release due to freezing to half the differene from the melting point
                  zdelwice = MIN(zdelwice,0.5_ireals*(t0_melt-t_so_new(i,kso))*zfak)
                ENDIF
                w_so_ice_new(i,kso) = w_so_ice_now(i,kso) + zdelwice
                t_so_new(i,kso) = t_so_new(i,kso) + zdelwice/zfak

             END IF                   ! m_stpy > 2
!            END IF                    ! land-points only
          END DO
       ENDDO
! End of heat transfer

  DO i = istarts, iends
!      IF (llandmask(i)) THEN  ! for landpoints only
        ! Next line has to be changed, if the soil surface temperature
        ! t_so(i,0,nnew) predicted by the heat conduction equation is used
        t_s_new   (i)    = t_so_new(i,1)
        t_so_new  (i,0)  = t_so_new(i,1)
        w_snow_new(i)  = w_snow_now(i) + zdt*zdwsndt  (i)/rho_w
 IF (itype_interception == 1) THEN
      w_i_new   (i)  = w_i_now(i) + zdt*zdwidt   (i)/rho_w
 ELSE IF (itype_interception == 2) THEN
      w_i_new   (i)  = w_i_now(i) + zdt*zdwidt   (i)/rho_w
      w_p_new (i)  = zwpnstr(i)
      w_s_new (i)  = zwisnstr(i)
 END IF


!      END IF          ! land-points only
     END DO



  DO i = istarts, iends
!!$!GZ:    ! *** Provisional fix for numerical instabilities in the presence of very
!!$        !     small amounts of snow *** !
!!$        ! melting-point adjustment of snow: if snow temp is above freezing and a non-negligible
!!$        ! amount of snow is available, then melt as much snow as needed to get snow temp
!!$        ! back to t0_melt while conserving energy
!!$        IF (t_snow_new(i) > t0_melt .AND. w_snow_new(i) > zepsi) THEN
!!$          w_snow_new(i) = MAX(w_snow_new(i)*(1._ireals-(t_snow_new(i)-t0_melt) &
!!$            *chc_i/lh_f), 0.0_ireals)
!!$          t_snow_new(i) = t0_melt
!!$        ELSE IF (w_snow_new(i) <= zepsi .OR. w_snow_new(i) > zepsi .AND. t_s_new(i) > t0_melt+15.0_ireals .AND. &
!!$          t(i) > t0_melt+15.0_ireals .AND. t_so_new(i,2) > t0_melt+15.0_ireals) THEN
!!$          ! if the amount of snow is negligible, or the environment is very warm, then just remove it
!!$          w_snow_new(i) = 0.0_ireals
!!$          ! If the snow has just melted, limit soil temperature increment to 2.5 deg C
!!$          ! in order to avoid nonsensically large temperature jumps
!!$
!!$          IF (w_snow_now(i) > zepsi) THEN
!!$!CDIR BEGIN EXPAND=7
!!$            t_so_new(i,1:7) = MIN(t_so_now(i,1:7)+2.5_ireals,t_so_new(i,1:7))
!!$            t_so_new(i,1:7) = MAX(t_so_now(i,1:7)-2.5_ireals,t_so_new(i,1:7))
!!$!CDIR END
!!$            t_so_new(i,0) = t_so_new(i,1)
!!$            t_s_new(i)    = t_so_new(i,1)
!!$          ENDIF
!!$          t_snow_new(i) = t_so_new(i,0)
!!$       ENDIF
!!$        ! *** End of provisional stability fix *** !
        !

        ! *** original code *** !
         IF (w_snow_new(i) <= zepsi) THEN
           w_i_new(i)    = w_i_new(i) + w_snow_new(i)
           w_snow_new(i) = 0.0_ireals
           t_snow_new(i) = t_so_new(i,0)
         ENDIF
        IF (w_i_new(i) <= 1.e-4_ireals*zepsi) w_i_new(i) = 0.0_ireals
     end DO
!<JH


  ! Eliminate snow for multi-layer snow model, if w_snow = 0
  IF (lmulti_snow) THEN
    DO ksn = 1, ke_snow
        DO i = istarts, iends
!          IF (llandmask(i)) THEN  ! for landpoints only
            IF (w_snow_new(i) <= zepsi) THEN
              t_snow_mult_new(i,ksn) = t_so_new(i,0)
              wliq_snow_new(i,ksn) = 0.0_ireals
              wtot_snow_new(i,ksn) = 0.0_ireals
              rho_snow_mult_new(i,ksn) = 0.0_ireals
              dzh_snow_new(i,ksn) = 0.0_ireals
            ENDIF
!          END IF          ! land-points only
        END DO
    END DO
  ENDIF


  IF(.NOT. lmulti_snow) THEN

      DO i = istarts, iends

!BR 7/2005 Update snow density
!
!     a) aging of existing snow
!
!     temperature dependence of relaxation/ageing constant

         zzz       = (t_snow_new(i)-csnow_tmin)/(t0_melt-csnow_tmin)
         ztau_snow = crhosmint+(crhosmaxt-crhosmint)*zzz
         ztau_snow = MAX(0.05_ireals,MIN(crhosmaxt,ztau_snow)) ! use 20 days in combination with temperature-dependent equilibrium density
         zrhosmax  = crhosmax_tmin+MAX(0._ireals,zzz)*(crhosmax_ml-crhosmax_tmin)
         zrho_snowe= MAX(rho_snow_now(i),zrhosmax+(rho_snow_now(i)-zrhosmax)* &
                     EXP(-ztau_snow*zdt/86400._ireals) )
!
!     b) density of fresh snow
!
         zrho_snowf= crhosminf+(crhosmaxf-crhosminf)* (zth_low(i)-csnow_tmin) &
                                                     /(t0_melt   -csnow_tmin)
         zrho_snowf= MAX(crhosminf,MIN(crhosmaxf,zrho_snowf))
!
!     c) new snow density is computed by adding depths of exisiting and new snow
!
         IF ( nclass_gscp >= 6 ) THEN
           zzz = (prs_gsp(i)+prs_con(i)+prg_gsp(i))*zdtdrhw
         ELSE
           zzz = (prs_gsp(i)+prs_con(i))*zdtdrhw
         ENDIF
         ! prevent accumulation of new snow if the air temperature is above 1 deg C with
         ! linear transition between 0.5 and 1 deg C
         IF (zzz > 0.5_ireals*zepsi .AND. zth_low(i) > t0_melt + 0.5_ireals) THEN
           ! part of the new snow that accumulates on the ground
           zxx = MAX(0._ireals, zzz*(t0_melt + 1._ireals - zth_low(i))*2._ireals)
           !
           ! the rest is transferred to the interception storage, soil moisture or runoff:
           w_i_new(i) = w_i_new(i) + zzz-zxx
           IF (w_i_new(i) > zwimax(i)) THEN  ! overflow of interception store
             zinf = w_i_new(i) - zwimax(i)
             w_i_new(i) = zwimax(i)
           ELSE
             zinf = 0._ireals
           ENDIF 
           zdwgme        = zinf*zrock(i)                    ! contribution to w_so
           zro           = (1._ireals - zrock(i))*zinf      ! surface runoff
           zredfu        = MAX( 0.0_ireals,  MIN( 1.0_ireals, (zw_fr(i,1) -  &
                           zfcap(i,1))/MAX(zporv(i,1)-zfcap(i,1), zepsi)))
           zdwgdt(i,1) = zdwgdt(i,1) + zdwgme*(1._ireals - zredfu)
           zro           = zro + zdwgme*zredfu    ! Infiltration not possible
                                                  ! for this fraction

           ! zro-, zdw_so_dt-correction in case of pore volume overshooting
           zw_ovpv = MAX(0._ireals, zw_fr(i,1)* zdzhs(1) * zrhwddt +  &
                     zdwgdt(i,1) - zporv(i,1) * zdzhs(1) * zrhwddt)
           zro = zro + zw_ovpv
           zdwgdt(i,1)= zdwgdt(i,1) - zw_ovpv

           runoff_s(i) = runoff_s(i) + zro*zroffdt

           ! correct SWE for immediately melted new snow
           zzz = zxx
           w_snow_new(i) = MIN(w_snow_new(i),w_snow_now(i)+zzz)
         ENDIF
         rho_snow_new(i)  = (w_snow_now(i)+zzz) / &
          ( MAX(w_snow_now(i),zepsi)/zrho_snowe + zzz/zrho_snowf )

    ! previous code based on weighted averaging of rho_snow
    !       znorm=MAX(w_snow_now(i)+(prs_gsp(i)+prs_con(i)+prg_gsp(i))      &
    !                 *zdtdrhw,zepsi)
    !       rho_snow_new(i)  = ( zrho_snowe*w_snow_now(i) + &
    !                           zrho_snowf*(prs_gsp(i)+prs_con(i)+prg_gsp(i)) &
    !                              *zdtdrhw )    /znorm
    !     ELSE
    !       znorm=MAX(w_snow_now(i)+(prs_gsp(i)+prs_con(i) )      &
    !                 *zdtdrhw,zepsi)
    !       rho_snow_new(i)  = ( zrho_snowe*w_snow_now(i) + &
    !                          zrho_snowf*(prs_gsp(i)+prs_con(i) ) &
    !                             *zdtdrhw )    /znorm
    !     ENDIF

         rho_snow_new(i) = MIN(crhosmax_ml,MAX(crhosmin_ml, rho_snow_new(i)))
!        END IF          ! land-points only

! New calculation of snow height for single layer snow model
         zdz_snow(i)=w_snow_new(i)*rho_w/rho_snow_new(i)
         h_snow_new(i) = zdz_snow(i)

         ! Calculation of top-layer snow density for two-layer snow density scheme
         IF (l2lay_rho_snow) THEN
           zrho_snowe = MAX(rho_snow_mult_now(i,1),zrhosmax+(rho_snow_mult_now(i,1)-zrhosmax)* &
                        EXP(-ztau_snow*zdt/86400._ireals) )
           zwsnow(i)  = MIN(max_toplaydepth,h_snow_gp(i))*rho_snow_mult_now(i,1)/rho_w
           rho_snow_mult_new(i,1) = (zwsnow(i)+zzz) / ( MAX(zwsnow(i),zepsi)/zrho_snowe + zzz/zrho_snowf )
           rho_snow_mult_new(i,1) = MIN(crhosmax_ml,MAX(crhosmin_ml, rho_snow_mult_new(i,1)))
           rho_snow_mult_new(i,2) = rho_snow_new(i)
         ENDIF

      END DO

  ELSE   ! new snow scheme

    DO i = istarts, iends
!      IF (llandmask(i)) THEN  ! for landpoints only
        h_snow_new(i) = 0.0_ireals
        sum_weight(i) = 0.0_ireals
!      END IF          ! land-points only
    END DO
    DO ksn = 1,ke_snow
        DO i = istarts, iends
!          IF (llandmask(i)) THEN  ! for landpoints only
            IF(w_snow_new(i) .GT. zepsi) THEN
              h_snow_new(i) = h_snow_new(i) + zdzh_snow(i,ksn)
            END IF
!          END IF          ! land-points only
        END DO
    END DO

    k = MIN(2,ke_snow-1)
    DO ksn = 1,ke_snow
      DO i = istarts, iends
        IF (w_snow_new(i) .GT. zepsi) THEN
          IF (ksn == 1) THEN ! Limit top layer to max_toplaydepth
            zhh_snow(i,ksn) = -MAX( h_snow_new(i)-max_toplaydepth, h_snow_new(i)/ke_snow*(ke_snow-ksn) )
          ELSE IF (ksn == 2 .AND. ke_snow > 2) THEN ! Limit second layer to 8*max_toplaydepth
            zhh_snow(i,ksn) = MIN( 8._ireals*max_toplaydepth+zhh_snow(i,1), zhh_snow(i,1)/(ke_snow-1)*(ke_snow-ksn) )
          ELSE ! distribute the remaining snow equally among the layers
            zhh_snow(i,ksn) = zhh_snow(i,k)/(ke_snow-k)*(ke_snow-ksn)
          ENDIF
        ENDIF
      END DO
    END DO

    DO ksn = ke_snow,1,-1
      DO i = istarts, iends
        IF (w_snow_new(i) .GT. zepsi) THEN
          dz_old(i,ksn) = dzh_snow_new(i,ksn)
          z_old(i,ksn) = -sum_weight(i) - dzh_snow_new(i,ksn)/2._ireals
          sum_weight(i) = sum_weight(i) + dzh_snow_new(i,ksn)
        END IF
      END DO
    END DO

    DO i = istarts, iends
!      IF (llandmask(i)) THEN  ! for landpoints only
        IF(w_snow_new(i) .GT. zepsi) THEN
          zhm_snow(i,1) = (-h_snow_new(i) + zhh_snow(i,1))/2._ireals
          dzh_snow_new(i,1) = zhh_snow(i,1) + h_snow_new(i)            !layer thickness betw. half levels of uppermost snow layer
          zdzm_snow(i,1     ) = zhm_snow(i,1) + h_snow_new(i)            !layer thickness between snow surface and main level of uppermost layer
          IF(dz_old(i,1).ne.0..and.rho_snow_mult_new(i,1).ne.0.) THEN
            wliq_snow_new(i,1) = wliq_snow_new(i,1)/dz_old(i,1)
          END IF
        END IF
!      END IF          ! land-points only
    END DO
    DO ksn = 2,ke_snow
        DO i = istarts, iends
!          IF (llandmask(i)) THEN  ! for landpoints only
            IF(w_snow_new(i) .GT. zepsi) THEN
              zhm_snow(i,ksn) = (zhh_snow(i,ksn) + zhh_snow(i,ksn-1))/2._ireals
              dzh_snow_new(i,ksn) = zhh_snow(i,ksn) - zhh_snow(i,ksn-1) ! layer thickness betw. half levels
              zdzm_snow(i,ksn     ) = zhm_snow(i,ksn) - zhm_snow(i,ksn-1) ! layer thickness betw. main levels
              IF(dz_old(i,ksn).ne.0..and.rho_snow_mult_new(i,ksn).ne.0.) THEN
                wliq_snow_new(i,ksn) = wliq_snow_new(i,ksn)/dz_old(i,ksn)
              END IF
            END IF
!          END IF          ! land-points only
        END DO
    END DO

    DO ksn = ke_snow,1,-1
      DO i = iends, istarts, -1
        t_new  (i,ksn) = 0.0_ireals
        rho_new(i,ksn) = 0.0_ireals
        wl_new (i,ksn) = 0.0_ireals
      END DO

      DO k = ke_snow,1,-1
          DO i = iends, istarts, -1
!            IF (llandmask(i)) THEN  ! for landpoints only
              IF(w_snow_new(i) .GT. zepsi) THEN

                weight = MIN(&
                     dz_old(i,k), &
                     z_old(i,k) + dz_old(i,k)/2._ireals &
                     - zhm_snow(i,ksn) + dzh_snow_new(i,ksn)/2._ireals , &
                     zhm_snow(i,ksn) + dzh_snow_new(i,ksn)/2._ireals &
                     - z_old(i,k) + dz_old(i,k)/2._ireals, &
                     dzh_snow_new(i,ksn))

                weight = (weight + ABS(weight)) * 0.5_ireals / dzh_snow_new(i,ksn)

                t_new  (i,ksn) = t_new  (i,ksn) + t_snow_mult_new  (i,k)*weight
                rho_new(i,ksn) = rho_new(i,ksn) + rho_snow_mult_new(i,k)*weight
                wl_new (i,ksn) = wl_new (i,ksn) + wliq_snow_new(i,k)*weight
              END IF
!            END IF          ! land-points only
          END DO
      END DO
    END DO

    DO ksn = ke_snow,1,-1
        DO i = istarts, iends
!          IF (llandmask(i)) THEN  ! for landpoints only
            IF(w_snow_new(i) > zepsi) THEN
              t_snow_mult_new  (i,ksn) = t_new  (i,ksn)
              rho_snow_mult_new(i,ksn) = rho_new(i,ksn)
              wtot_snow_new    (i,ksn) = rho_new(i,ksn)*dzh_snow_new(i,ksn)/rho_w
              wliq_snow_new    (i,ksn) = wl_new (i,ksn)*dzh_snow_new(i,ksn)
            END IF
!          END IF          ! land-points only
        END DO
    END DO

    DO i = istarts, iends
!      IF (llandmask(i)) THEN  ! for landpoints only
        IF(w_snow_new(i) > zepsi) THEN
          rho_snow_new(i) = w_snow_new(i)/h_snow_new(i)*rho_w
        ELSE !JH
          rho_snow_new(i) = 250._ireals ! workaround need to be inspected!!
        END IF
        IF(w_snow_new(i) > zepsi) THEN
          ! linear extrapolation from t_snow_mult_new(i,2) and t_snow_mult_new(i,1) to t_snow_mult_new(i,0)
          t_snow_mult_new(i,0) = (t_snow_mult_new(i,1)*(2._ireals*dzh_snow_new(i,1)+dzh_snow_new(i,2))- &
                                t_snow_mult_new(i,2)*dzh_snow_new(i,1))/(dzh_snow_new(i,1)+dzh_snow_new(i,2))
          ! limiter to prevent unphysical values and/or numerical instabilities
          t_snow_mult_new(i,0) = MIN(273.15_ireals,t_snow_mult_new(i,0),t_snow_mult_new(i,1)+5.0_ireals)
          t_snow_mult_new(i,0) = MAX(t_snow_mult_new(i,0),t_snow_mult_new(i,1)-5.0_ireals)
        END IF
        t_snow_new(i) = t_snow_mult_new (i,0)
!      END IF          ! land-points only
    END DO

 ENDIF ! lmulti_snow


  DO kso = 1,ke_soil
      DO i = istarts, iends
!        IF (llandmask(i)) THEN  ! for landpoints only
          w_so_new(i,kso) = w_so_now(i,kso) + zdt*zdwgdt(i,kso)/rho_w
!        END IF  ! land-points only
      END DO
  END DO        ! soil layers

  ! Update of two-time level interface variables
  DO i = istarts, iends
!    IF (llandmask(i)) THEN  ! for landpoints only
     h_snow(i) = h_snow_new(i)
!    END IF
  END DO

!---loop over tiles---
!END DO
!---------------------

!  DO ns = nsubs0, nsubs1
  !em>

  ! computation of the temperature at the boundary soil/snow-atmosphere
  IF (ldiag_tg) THEN
    CALL tgcom ( t_g(:), t_snow_new(:), t_s_new(:),       &
                 w_snow_new(:), zf_snow(:), ie, cf_snow, &
                 istarts, iends )
  ENDIF
!  END DO

  ! computation of the weighted turbulent fluxes at the boundary surface-atmosphere
  IF(PRESENT(zshfl_sfc)) THEN
    DO i = istarts, iends
      zshfl_sfc(i) = zshfl_s(i)*(1._ireals - zf_snow(i)) + zshfl_snow(i)*zf_snow(i)
      zlhfl_sfc(i) = zlhfl_s(i)*(1._ireals - zf_snow(i)) + zlhfl_snow(i)*zf_snow(i)

!DR start
      zqhfl_sfc(i) = zqhfl_s(i)*(1._ireals - zf_snow(i)) + zqhfl_snow(i)*zf_snow(i)
!DR end

!        zlhfl_s(i) = (zts_pm(i)*lh_v + (1._ireals-zts_pm(i))*lh_s)*zverbo(i) &
!                     / MAX(zepsi,(1._ireals - zf_snow(i)))  ! take out (1-f) scaling
!        zlhfl_snow(i) = lh_s*zversn(i)
!      zlhfl_sfc(i) = zverbo(i) + zversn(i)*zf_snow(i)
    END DO
  END IF

! Debug messages


#ifdef __ICON__
  IF (msg_level >= 19) THEN
!     DO ic=1,icount_snow
!       i=melt_list(ic)
      DO i = istarts, iends

        IF (t_snow_new(i) > 355._ireals .OR. t_s_now(i) > 355._ireals .OR. t_s_new(i) > 355._ireals .OR. &
            t_snow_new(i) < 190._ireals .OR. t_s_now(i) < 200._ireals .OR. t_s_new(i) < 200._ireals) THEN

!        IF ((t_snow_new(i) > t0_melt .AND. w_snow_new(i) > zepsi).OR.&
!   (w_snow_new(i) <= zepsi .OR. w_snow_new(i) > zepsi .AND. t_s_new(i) > t0_melt+15.0_ireals .AND. &
!         t(i) > t0_melt+15.0_ireals .AND. t_so_new(i,2) > t0_melt+15.0_ireals)) THEN

              write(0,*) "SFC-DIAGNOSIS TERRA ",i !,icount_snow!,ntstep
              write(0,*) "t",i, t(i)
              write(0,*) "p0",i,p0(i)
              write(0,*) "qv",i,qv(i)
              write(0,*) "ps",i,ps(i)
              write(0,*) "t_g",i,t_g(i)
              write(0,*) "t_s",i,t_s_now(i),t_s_new(i)
              write(0,*) "t_snow",i,t_snow_now(i),t_snow_new(i)
              write(0,*) "w_snow",i,w_snow_now(i),w_snow_new(i)
              write(0,*) "h_snow",i,h_snow_now(i),h_snow_new(i)
              write(0,*) "rho_snow",i,rho_snow_now(i),rho_snow_new(i)
              write(0,*) "qv_s",i,qv_s(i)
              write(0,*) "t_so_now",i,t_so_now(i,:)
              write(0,*) "t_so_new",i,t_so_new(i,:)
              write(0,*) "w_so_now",i,w_so_now(i,:)
              write(0,*) "w_so_new",i,w_so_new(i,:)
              IF (lmulti_snow) THEN
                write(0,*) "t_snow_mult",i,t_snow_mult_now(i,:),t_snow_mult_new(i,:)
                write(0,*) "t_snow_mult_now",i,t_snow_mult_now(i,0),t_snow_mult_now(i,1),t_snow_mult_now(i,2)
                write(0,*) "t_snow_mult_new",i,t_snow_mult_new(i,0),t_snow_mult_new(i,1),t_snow_mult_new(i,2)
                write(0,*) "rho_snow_mult",i,rho_snow_mult_now(i,:),rho_snow_mult_new(i,:)
                write(0,*) "wliq_snow",i,wliq_snow_now(i,:),wliq_snow_new(i,:)
                write(0,*) "wtot_snow",i,wtot_snow_now(i,:),wtot_snow_new(i,:)
                write(0,*) "dzh_snow",i,dzh_snow_now(i,:),dzh_snow_new(i,:)
              ENDIF
              write(0,*) "tch_t",i,tch(i)
              write(0,*) "tcm_t",i,tcm(i)
!              write(0,*) " tfv_t",i,tfv(i)
              write(0,*) "zshfl,zlhfl,zradfl,zg1",i,zshfl(i),zlhfl(i),zradfl(i),zg1(i)
              write(0,*) "zshfl_s,zlhfl_s",i, zshfl_s(i), zlhfl_s(i)
              write(0,*) "zshfl_snow,zlhfl_snow,zf_snow,zgsb",i,zshfl_snow(i),zlhfl_snow(i), zf_snow(i),zgsb(i)
              write(0,*) "soiltyp_t",i,soiltyp_subs(i)
              write(0,*) "plcov_t",i,  plcov(i)
!              write(0,*) "rootdp_t",i, rootdp(i)
!              write(0,*) "sai_t",i,   sai(i)
!              write(0,*) "tai_t",i,   tai(i)
!              write(0,*) "eai_t",i,   eai(i)
!             write(0,*) "t_2m_t",i,  t_2m(i)
!              write(0,*) "u_10m_t",i, u_10m(i)
!              write(0,*) "v_10m_t",i, v_10m(i)
              write(0,*) "sobs_t",i,  sobs(i)
              write(0,*) "thbs_t",i,  thbs(i)
              write(0,*) "pabs_t",i,  pabs(i)
              write(0,*) "prr_gsp, prr_con, prs_gsp, prs_con",i,prr_gsp(i), prr_con(i), prs_gsp(i), prs_con(i)
              write(0,*) "zrr, zrs",i, zrr(i), zrs(i)
              write(0,*) "zsprs, ztchv, ztchv_max",i, zsprs(i), ztchv(i), ztchv_max(i), ztchv(i)/MAX(SQRT(u(i)**2+v(i)**2),zepsi)
!              write(0,*) "llandmask_t",i,llandmask(i)
           END IF
         END DO
  ENDIF
#endif



!------------------------------------------------------------------------------
! End of module procedure terra_multlay
!------------------------------------------------------------------------------


END SUBROUTINE terra_multlay


!==============================================================================


SUBROUTINE tgcom (tg, ts, tb, ws, snowfrac, ie, cf_snow, istart, iend)

!-------------------------------------------------------------------------------
!
! Description:
!   Computation of the temperature tg at the boundary layer between the ground
!   and the atmosphere. Only 2-dimensional arrays can be passed to tgcom. It
!   must be called using the desired time level.
!
! Method:
!   For grid points above water and for grid points on land that are not
!   covered with snow:   tg = ground surface temperature tb
!   For snow covered land points, tg is a function of the temperature of the
!   the snow surface ts and the ground surface temperature tb:
!       tg = ts + exp( -rhde*ws ) * (tb-ts)
!   from Version 2.18 on replaced by
!       tg = ts + ( 1. - MIN(1.,ws/cf_snow)) * (tb -ts)
!
!-------------------------------------------------------------------------------

! Parameter list:

INTEGER (KIND=iintegers), INTENT (IN)    ::    &
  ie,              & ! dimensions of the fields
  istart, iend       ! start and end-indices of the computation

REAL (KIND=ireals), INTENT (INOUT)       ::    &
  tg (ie)    ! temperature at the boundary between ground and atmosphere

REAL (KIND=ireals), INTENT (IN)          ::    &
!DR  ts (ie), & ! temperature of the snow surface
  tb (ie), & ! temperature of the ground surface
  ws (ie)    ! water content of snow

REAL (KIND=ireals), INTENT (INOUT)          ::    &
  ts (ie),   & ! temperature of the snow surface
  snowfrac(ie) ! snow-cover fraction

REAL (KIND=ireals), INTENT (IN)          ::    &
  cf_snow       ! factor for the computation

INTEGER :: i

!-------------------------------------------------------------------------------

! Begin subroutine tgcom

  DO i = istart, iend
    snowfrac(i) = MIN(1.0_ireals,ws(i)/cf_snow)
    tg(i) = ts(i) + (1.0_ireals - snowfrac(i))*(tb(i) - ts(i))
  ENDDO

END SUBROUTINE tgcom


!==============================================================================

REAL (KIND = ireals) FUNCTION zsf_heav (zstx) !Heaviside function
   implicit none
   real (KIND=ireals   ), intent(in)  :: zstx
   zsf_heav      = 0.5_ireals+SIGN( 0.5_ireals, zstx )
END FUNCTION zsf_heav

REAL (KIND = ireals) FUNCTION zsf_psat_iw  (zstx,z2iw   ,z4iw)! Saturation water vapour pressure over ice or water
                                                              ! depending on temperature "zstx"
   implicit none
   real (KIND=ireals   ), intent(in)  :: zstx,z2iw   ,z4iw
   zsf_psat_iw   = b1*EXP(z2iw*(zstx - b3)/(zstx - z4iw))
END FUNCTION zsf_psat_iw


REAL (KIND = ireals) FUNCTION zsf_qsat  (zspsatx, zspx     ) ! Specific humidity at saturation pressure
                                                             ! (depending on the saturation water vapour pressure
   implicit none
   real (KIND=ireals   ), intent(in)  :: zspsatx, zspx
  zsf_qsat      = rdv*zspsatx/(zspx-o_m_rdv*zspsatx)
END FUNCTION zsf_qsat

REAL (KIND = ireals) FUNCTION zsf_dqvdt_iw (zstx,zsqsatx,z4iw,z234iw) ! First derivative of specific saturation humidity
                                                                      ! with respect to temperature (depending on temperature
                                                                      ! "zstx" and saturation specific humidity pressure
                                                                      ! "zsqsatx")
   implicit none
   real (KIND=ireals   ), intent(in)  :: zstx,zsqsatx,z4iw,z234iw
  zsf_dqvdt_iw  = z234iw*(1._ireals+rvd_m_o*zsqsatx)*zsqsatx/(zstx-z4iw)**2
END FUNCTION zsf_dqvdt_iw


REAL (KIND = ireals) FUNCTION watrcon_RT(ks,lw,kw1,pv,adp)
   implicit none
   real (KIND=ireals   ), intent(in)  :: ks,lw,kw1,pv,adp
   watrcon_RT=ks*EXP(kw1*(pv-lw)/(pv-adp))
END FUNCTION watrcon_RT

REAL (KIND = ireals) FUNCTION watrdiff_RT(ds,lw,dw1,pv,adp)
   implicit none
   real (KIND=ireals   ), intent(in)  ::  ds,lw,dw1,pv,adp
   watrdiff_RT = ds*EXP(dw1*(pv-lw)/(pv-adp))
END FUNCTION watrdiff_RT


!------------------------------------------------------------------------------
! End of module src_soil_multlay
!------------------------------------------------------------------------------




END MODULE mo_soil_ml

