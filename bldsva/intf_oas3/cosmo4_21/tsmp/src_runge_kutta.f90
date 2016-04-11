! $RCSfile: src_runge_kutta.f90,v $
! $Revision: 4.21 $ $Date: 2012/03/30 $
!+ Source module for 2-timelevel Runge-Kutta version of the dynamics
!------------------------------------------------------------------------------

MODULE src_runge_kutta

!------------------------------------------------------------------------------
!
! Description:
!   The module "src_runge_kutta" performs one time step of the integration
!   of the spatially discretized thermodynamical equations using one of  
!   several selectable 2-timelevel Runge-Kutta schemes.
!   Normal RK-variants of 2nd- or 3rd-order (in principle it is also possible
!   to use a simple 1st-order Euler forward scheme) as well as a slightly
!   more sophisticated 3rd-order TVD-Runge-Kutta (total variation diminishing)
!   scheme are implemented.
!   The type of scheme is set using the namelist-Parameter irunge_kutta
!   (0: Almuts scheme, 1: normal RK-scheme, 2: TVD-RK-scheme) and
!   the order of the scheme (makes only sense for case 1) is set with
!   the namelist-Parameter irk_order (value: 1, 2 or 3).
!   In addition, the order of the operator for the horizontal advection of
!   the variables u, v, w, pp an T is choosen via the namelist-Parameter
!   iadv_order (value: 3, 4, 5 or 6).
!   Best choice imho: irunge_kutta=2, irk_order=3, iadv_order=5.
!   In contrary to Almuts scheme the RK-loop contains the effects
!   of the complete slow tendencies - especially the whole 3D advection is
!   computed in each RK-step - and the integration of the small time steps
!   done in the routine "fast_waves_runge_kutta".
!   The option to use vertical explicit advection is prepared
!   (lva_impl_dyn=.FALSE.) but not yet working correctly.
!   Driving routine is the model procedure "org_runge_kutta", which
!   calls the routines required. 
!
! Current Code Owner: DWD, Jochen Foerstner    
!  phone:  +49  69  678667 35
!  fax:    +49  69  8062 3721
!  email:  jochen.foerstner@dwd.de
!
! History:
! Version    Date       Name
! ---------- ---------- ----
! 3.7        2004/02/18 Jochen Foerstner
!  Initial release
! 3.8        2004/03/23 Jochen Foerstner
!  Technical changes and some optimizations
! 3.13       2004/12/03 Jochen Foerstner
!  Adaptations for calls to boundary exchange; more modifications
!  Modifications to run with latent heat nudging (Klaus Stephan et.al.)
! 3.14       2005/01/25 Jochen Foerstner
!  Enhanced formulation of moist turbulence for 3D turbulence scheme
! 3.16       2005/07/22 Jochen Foerstner
!  New formulation of advection for the moisture quantities
!  (Courant number independent Eulerian and Semi-Lagrangian variants)
!  Optimizations
! 3.18       2006/03/03 Jochen Foerstner / Klaus Stephan
!  Conservation form of Bott Scheme is now standard;
!  Changed interfaces and some additional optimizations
!  LHN namelist parameter moved to data_lheat_nudge to avoid to many dependencies
! 3.19       2006/04/25 Michael Baldauf
!  Corrections for computing small time step
! 3.21       2006/12/04 Jochen Foerstner / Ulrich Schaettler
!  More updates
!  Added additional terms for deep atmosphere (R. Petrik)
!  Put subroutines for slow tendencies (to src_slow_tendencies_rk) and 
!  advection (to src_advection_rk) to extra modules
! V3_23        2007/03/30 Jochen Foerstner, Michael Baldauf, Lucio Torrisi
!  Bug correction for computing CFL criterion in SR check_cfl_horiz_advection
!  Apply Rayleigh damping only, where necessary
!  Eliminated SR exchange_boundaries_nnow
!  Adaptations to work with the DFI (by Lucio Torrisi)
! V4_1         2007/12/04 Ulrich Schaettler
!  Compute MAXVAL for zwmax only for istartpar:iendpar,jstartpar:jendpar
!  to get reproducible results
! V4_4         2008/07/16 Ulrich Schaettler, Michael Baldauf
!  Adapted interface of get_timings
!  Adaptations to changes in src_advection_rk
!  Replaced lkainfri by itype_conv
!  Introduced tendencies of sub-grid scale orography scheme (J-P Schulz)
! V4_5         2008/09/10 Guenther Zaengl, Dmitrii Mironov, Michael Baldauf
!  Adaptations for new reference atmosphere (Guenther Zaengl)
!  Deliver tendencies of qi and qc of the new convection scheme for the
!    cloud microphysics scheme. (D. Mironov, M. Baldauf)
! V4_8         2009/02/16 Michael Baldauf
!  Stable integration of coriolis terms (is still hardcoded to false; 
!  for using it, set l_Coriolis_every_RK_substep=.TRUE.)
!  Bug fix in the computation of dt0dz for irefatm=2
!    (hhl has to be used on main levels)
! V4_9         2009/07/16 Guenther Zaengl
!  Option for using potential temperature as prognostic variable in RK part
!  Check calls to get_timings with ltime (U. Schaettler)
!  Inserted compiler directives
!  Inserted call to SR collapse to change loop indices
! V4_11        2009/11/30 Guenther Zaengl
!  Implemented optional use of irefatm=3 with constant BV frequency
! V4_11        2010/01/18 Volker Kuell
!  Additional fields for HYMACS convection scheme, inclusion of convective pressure tendencies
! V4_11        2011/06/08 Volker Kuell
!  all global HYMACS variables exported to data_hymacs.f90 and all HYMACS related code put into
!  ifdef blocks
! V4_11      2011/06/20 Volker Kuell
!  All HYMACS related global Variables renamed to hymacs_...
! V4_12        2010/05/11 Oli Fuhrer
!  Replaced DOUBLE PRECISION by KIND-definition in SR coriolis
!  Adapted call to SR comp_hori_diffusion and eliminated option itype_hdiff=3
! V4_13        2010/05/11 Michael Gertz
!  Adaptions to SVN
! V4_15        2010/11/19 Michael Baldauf
!  Possibility to call a 3rd order implicit vertical advection.
!  New subroutine 'calc_wcon_sqrtg'
! V4_17        2011/02/24 Ulrich Blahak
!  Adapted interface of exchg_boundaries; corrected kzdims(1:20) -> kzdims(1:24);
!  eliminated my_peri_neigh; added possibility to output max. V_h and min./max. W
!  to stdout for diagnostic purposes, if (ldebug_dyn .AND. idbg_level > 3).
! V4_18        2011/05/26 Ulrich Schaettler, Michael Baldauf
!  Introduced conditional compilation for Nudging
!  Bug fix: Subr. 'coriolis' is now called with the correct timelevels
!     in the case l_Coriolis_every_RK_substep=.TRUE. (Michael Baldauf)
!  Implementation of COSMO-ART interfaces (Christoph Knote)
! V4_20        2011/08/31 Michael Baldauf
!  Modified number of small time steps for the Runge-Kutta substeps
!  Call SR coriolis only, if lcori is true
!  Initialize some variables in any case, not only for itheta_avd >= 1 (Oli Fuhrer)
!  Print debug output only in case of ldebug_dyn
! V4_21        2012/03/15 Volker Kuell
!  Additional fields for HYMACS convection scheme, inclusion of convective pressure tendencies
! V4_21        2012/03/30 Alexander Kelbch
!  Inclusion of tracer (adopted from Markus Uebel)
! V4_21        2012/09/19 Markus Uebel
!  Inclusion of CO2 coupling with Community Land Model (CLM) (ifdef COUP_OAS_COS) 
! V4_21        2013/02/03 Prabhakar Shrestha
!  Inclusion of transfer coefficeint for moisture
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!==============================================================================
!
! Declarations:
!
! Modules used:

USE data_parameters , ONLY :   &
  ireals,    & ! KIND-type parameters for real variables
  iintegers    ! KIND-type parameter for standard integer variables

!------------------------------------------------------------------------------

USE data_modelconfig, ONLY :   &
  
  irefatm,      & ! 1: old reference atm. based on dT/dlnp = const
                  ! 2: new reference atm. with asymptotically isothermal stratosphere
  dt0lp,        & ! d (t0) / d (ln p0)
  h_scal,       & ! scale height for irefatm = 2
  delta_t,      & ! temp. difference between sea level and stratosphere (for irefatm = 2)
  bvref,        & ! reference Brunt-Vaisala-frequency for irefatm=3
                  ! (available for idealized simulations only)

  ! 2. horizontal and vertical sizes of the fields and related variables
  ! --------------------------------------------------------------------
  ie,           & ! number of grid points in zonal direction
  je,           & ! number of grid points in meridional direction
  ke,           & ! number of grid points in vertical direction
  ke1,          & ! ke+1
  
  ! 3. start- and end-indices for the computations in the horizontal layers
  ! -----------------------------------------------------------------------
  !    These variables give the start- and the end-indices of the 
  !    forecast for the prognostic variables in a horizontal layer.
  !    Note, that the indices for the wind-speeds u and v differ from 
  !    the other ones because of the use of the staggered Arakawa-B-grid.
  !    
  !   zonal direction
  istart,       & ! start index for the forecast of w, t, qv, qc, tracer (AK) and pp
  iend,         & ! end index for the forecast of w, t, qv, qc, tracer (AK) and pp
  istartu,      & ! start index for the forecast of u
  iendu,        & ! end index for the forecast of u
  istartv,      & ! start index for the forecast of v
  iendv,        & ! end index for the forecast of v
  istartpar,    & ! start index for computations in the parallel program
  iendpar,      & ! end index for computations in the parallel program
  
  !   meridional direction
  jstart,       & ! start index for the forecast of w, t, qv, qc, tracer (AK) and pp
  jend,         & ! end index for the forecast of w, t, qv, qc, tracer (AK) and pp
  jstartu,      & ! start index for the forecast of u
  jendu,        & ! start index for the forecast of u
  jstartv,      & ! start index for the forecast of v
  jendv,        & ! end index for the forecast of v
  jstartpar,    & ! start index for computations in the parallel program
  jendpar,      & ! end index for computations in the parallel program
  
  ! 4. constants for the horizontal rotated grid and related variables
  ! ------------------------------------------------------------------

  dlon,         & ! grid point distance in zonal direction (in degrees)
  dlat,         & ! grid point distance in meridional direction (in degrees)

  ! 5. variables for the time discretization and related variables
  ! --------------------------------------------------------------

  dt,           & ! long time-step
  betasw,       & ! beta-variable for treatment of soundwaves    in w
  betagw,       & ! beta-variable for treatment of gravity-waves in w
  beta2sw,      & ! beta-variable for treatment of soundwaves    in p*, T*
  beta2gw,      & ! beta-variable for treatment of gravity-waves in p*, T*
  kflat           ! level-index where half-levels bcome flat


! end of data_modelconfig

!------------------------------------------------------------------------------

USE data_constants  , ONLY :   &
  
  pi,           & ! circle constant

  ! 2. physical constants and related variables
  ! -------------------------------------------

  r_d,          & ! gas constant for dry air
  rdv,          & ! r_d / r_v
  o_m_rdv,      & ! 1 - r_d/r_v
  rvd_m_o,      & ! r_v/r_d - 1
  cp_d,         & ! specific heat of dry air at constant pressure
  cpdr,         & ! 1 / cp_d
  gamma,        & ! = cp_d / cv_d
  lh_v,         & ! latent heat of vapourization
  g,            & ! acceleration due to gravity
  r_earth,      & ! mean radius of the earth
  

  ! 3. constants for parametrizations
  ! ---------------------------------

  b1,           & ! variables for computing the saturation vapour pressure
  b2w,          & ! over water (w) and ice (i)
  b3,           & !               -- " --
  b4w,          & !               -- " --
  b234w,        & ! b2w * (b3 - b4w)
  aks4            ! variable for horizontal diffusion of fourth order

! end of data_constants

!------------------------------------------------------------------------------

USE data_fields     , ONLY :   &
  
  ! 1. constant fields for the reference atmosphere                 (unit)
  ! -----------------------------------------------
  rho0       ,    & ! reference density at the full model levels    (kg/m3)
  p0         ,    & ! reference pressure at main levels             ( Pa  )
  dt0dz      ,    & ! temperature gradient of reference atmosphere  ( K/m )
  t0         ,    & ! reference temperature                         ( K   )
  hhl        ,    & ! geometical height of half model levels        ( m )
  
  ! 2. external parameter fields                                    (unit)
  ! ----------------------------
  fc         ,    & ! coriolis-parameter                            ( 1/s )
  fccos      ,    & ! coriolis-parameter mit cosinus                ( 1/s )
  crlat      ,    & ! cosine of transformed latitude
  
  ! 3. prognostic variables                                         (unit)
  ! -----------------------
  u          ,    & ! zonal wind speed                              ( m/s )
  v          ,    & ! meridional wind speed                         ( m/s )
  w          ,    & ! vertical wind speed (defined on half levels)  ( m/s )
  t          ,    & ! temperature                                   (  k  )
  qv         ,    & ! specific water vapor content                  (kg/kg)
  qc         ,    & ! specific cloud water content                  (kg/kg)
  qi         ,    & ! specific cloud ice   content                  (kg/kg)
  qr         ,    & ! specific rain water  content                  (kg/kg)
  qs         ,    & ! specific snow        content                  (kg/kg)
  qg         ,    & ! specific graupel     content                  (kg/kg)
  tke        ,    & ! turbulent kinetic energy (on half levels)     (m2/s2)
  pp         ,    & ! deviation from the reference pressure         ( pa  )
  
  ! 4. tendency fields for the prognostic variables                 (unit )
  ! -----------------------------------------------
  !    timely deviation  by diabatic and adiabatic processes 
  !    without sound-wave terms
  utens,          & ! u-tendency without sound-wave terms           ( m/s2)
  vtens,          & ! v-tendency without sound-wave terms           ( m/s2)
  wtens,          & ! w-tendency without sound-wave terms           ( m/s2)
  ttens,          & ! t-tendency without sound-wave terms           ( K/s )
  qvtens,         & ! qv-tendency                                   ( 1/s )
  qctens,         & ! qc-tendency                                   ( 1/s )
  qitens,         & ! qi-tendency                                   ( 1/s )
  pptens            ! pp-tendency without sound-wave terms          ( Pa/s )

USE data_fields     , ONLY :   &
  
  ! 6. fields that are computed in the parametrization and dynamics (unit )
  ! ---------------------------------------------------------------
  qvt_diff  ,   & ! humidity    tendency  due to diffusion          ( 1/s )
  tinc_lh,      & ! temperature increment due to latent heat        (  K  )

  !   turbulent coefficients in the atmosphere
  !   (defined on half levels)
                    ! vertical   turbulent diffusion coefficients
  !   fields of the radiation
  sohr,         & ! rate of solar heating                          ( K/s )
  thhr,         & ! rate of thermal heating                        ( K/s )
  
  !   fields of the precipitation
  tt_conv,      & ! temperature tendency due to convection        (K /s )
  qvt_conv,     & ! humidity    tendency due to convection        ( s^-1)
  qct_conv,     & ! qc-tendency tendency due to convection        ( s^-1)
  qit_conv,     & ! qi-tendency tendency due to convection        ( s^-1)
  ut_conv,      & ! u-tendency due to convection                  (ms^-2)
  vt_conv,      & ! v-tendency due to convection                  (ms^-2)
  
  !   fields from the sub-grid scale orography scheme
  ut_sso   ,    & ! u-tendency due to SSO                         ( m/s2)
  vt_sso   ,    & ! v-tendency due to SSO                         ( m/s2)
  tt_sso   ,    & ! temperature tendency due to SSO               ( K/s )

  !   fields that are computed in the dynamics
  wcon     ,    & ! contravariant vertical velocity
  uadvt    ,    & ! advective tendency of u
  vadvt    ,    & ! advective tendency of v
  wadvt    ,    & ! advective tendency of w
  ppadvt   ,    & ! advective tendency of pp
  tadvt           ! advective tendency of t

! end of data_fields

!AK (30.03.2012)
USE data_tracer     , ONLY :   &
  tracer,       & ! tracer concentration                          ( "var" )
  tracertens  , & ! tracer tendency                               ("var"/s)
  tracert_conv, & ! tracer tendency due to convection             ("var"/s)
!MU (12.04.13)
#ifdef COUP_OAS_COS
  co2fl_s,      & ! CO2 tendency due to land-atmosphere exchange  (  1/s  )
  psn_tens,     & ! CO2 tendency due to photosynthesis            (  1/s  )
  plres_tens,   & ! CO2 tendency due to plant respiration         (  1/s  )
#endif
!MU (12.04.13)
  ltracer,      & ! switch for tracer
  ntracer_dim,  & ! number of tracers for dimensioning of fields
  ntracer_adv,  & ! number of tracers with advection
  ntracer         ! number of tracer
!AK (30.03.2012)

#ifdef HYMACS
! VK 2012/03/15
USE data_hymacs     , ONLY :   &
  hymacs_ppt_conv ! pressure tendency due to convection           (Pa/s )
#endif

!------------------------------------------------------------------------------

USE data_runcontrol , ONLY :   &
  
  ! 1. start and end of the forecast
  ! --------------------------------
  nstart,        & ! first time step of the forecast
  ntstep,        & ! actual time step
  ! indices for permutation of three time levels
  nnow,          & ! corresponds to ntstep
  nnew,          & ! corresponds to ntstep + 1
  
  ! 3. controlling the physics
  ! --------------------------
  lphys,         & ! forecast with physical parametrizations
  lprog_qi,      & ! if .TRUE., running with cloud ice
  lprogprec,     & ! forecast with prognostic rain, snow and graupel
  lsso,          & ! forecast with sub-grid scale orography scheme
  itype_gscp,    & ! type of microphys. parametrization
  itype_turb,    & ! type of turbulent diffusion parametrization
  itype_conv,    & ! type of convection parametrization
  l3dturb,       & ! 3D-turbulence (additional horizontal diffusion)
  lprog_tke,     & ! prognostic treatment of TKE (for itype_turb=5/7)
  l_cosmo_art,   & ! if .TRUE., run the COSMO_ART
  l_pollen,      & ! if .TRUE., run the Pollen

  ! 4. controlling the dynamics
  ! ---------------------------
  lcori,         & ! lartif_data=.TRUE.:  with Coriolis force
  lcori_deep,    & ! if =.TRUE.: take cos(phi) coriolis terms into account

  ! 6. controlling the upper boundary condition
  ! -------------------------------------------
  lrubc            ! with radiative upper boundary condition
  
USE data_runcontrol , ONLY :   &

  ! 7. additional control variables
  ! -------------------------------
  lcond,         & ! forecast with condensation/evaporation
  ldiabf_lh,     & ! include diabatic forcing due to latent heat in RK-scheme
  ldiabf_satad,  & ! include diabatic forcing due to saturation adjustment
  lperi_x,       & ! if lartif_data=.TRUE.: periodic boundary conditions in x-dir.
                   ! or with Davies conditions (.FALSE.)
  lperi_y,       & ! if lartif_data=.TRUE.: periodic boundary conditions in y-dir.
                   ! or with Davies conditions (.FALSE.)
  l2dim,         & ! 2 dimensional runs
  ltime,         & ! detailed timings of the program are given
  y_scalar_advect,&! type of scalar advection scheme
                   ! "SL3_MF", "SL3_SFD", "Bott2", "Bott4"
                   ! "Bott2_Strang", "Bott4_Strang", "vanLeer", "PPM"
  y_vert_adv_dyn,& ! switch to choose between explicit and several implicit schemes
                   ! for the vertical advection of dynamic variables
  lreproduce,    & ! the results are reproducible in parallel mode
  lhordiff,      & ! running with horizontal diffusion
  nbl_exchg,     & ! number of boundlines to exchange
  irunge_kutta,  & ! =0: use Almuts scheme, =1: use new RK scheme,
                   ! =2: use new TVD-RK scheme
  irk_order,     & ! order of the Runge-Kutta time-integration scheme
  iadv_order,    & ! order of the horizontal advection scheme in dyn. core
  ieva_order,    & ! order of the explicit vertical adv. scheme in dyn. core
  itheta_adv,    & ! =0: use T' (perturbation temperature) for advection
                   ! =1: use theta' (perturbation potential temperature)
                   ! =2: use theta (full potential temperature)
  xkd,           & ! coefficient for divergence damping
  itype_hdiff,   & ! type of horizontal diffusion (=1: 4th order linear),
                   ! =2: 4th order linear monotonic with orographic limit)
  hd_corr_q_bd,  & ! correction factor for horizontal diffusion fluc of qv,qc in boundary zone
  hd_corr_q_in,  & ! correction factor for horizontal diffusion fluc of qv,qc in domain

! 12. controlling verbosity of debug output
! -----------------------------------------
  idbg_level,    & ! to control the verbosity of debug output
  ldebug_dyn,    & ! if .TRUE., debug output for dynamics
  lprintdeb_all    ! .TRUE.:  all tasks print debug output
                   ! .FALSE.: only task 0 prints debug output
! end of data_runcontrol 

!------------------------------------------------------------------------------

USE data_parallel,      ONLY :  &
  num_compute,   & ! number of compute PEs
  nboundlines,   & ! number of boundary lines of the domain for which
                   ! no forecast is computed = overlapping boundary
                   ! lines of the subdomains
  ldatatypes,    & ! if .TRUE.: use MPI-Datatypes for some communications
  ltime_barrier, & ! if .TRUE.: use additional barriers for determining the
                   ! load-imbalance
  ncomm_type,    & ! type of communication
  my_cart_id,    & ! rank of this subdomain in the cartesian communicator
  my_cart_neigh, & ! neighbors of this subdomain in the cartesian grid
  icomm_cart,    & ! communicator for the virtual cartesian topology
                   ! that can be used by MPI_WAIT to identify the send
  imp_reals,     & ! determines the correct REAL type used in the model
                   ! for MPI
  imp_integers,  &
  sendbuf,       & ! sending buffer for boundary exchange:
                   ! 1-4 are used for sending, 5-8 are used for receiving
  isendbuflen      ! length of one column of sendbuf

!------------------------------------------------------------------------------

! coefficients used in complete_tendencies
! ----------------------------------------

USE src_slow_tendencies_rk, ONLY:  &
  zsqrtgrho_r_s,        & ! reciprocal square root of G * rho at skalar points
  zsqrtgrho_r_u,        & ! reciprocal square root of G * rho at u points
  zsqrtgrho_r_v,        & ! reciprocal square root of G * rho at v points
  zsqrtgrho_r_w,        & ! reciprocal square root of G * rho at w points
  za1t,                 & !
  za2t,                 & !
  zpia,                 & !
  zpianf,               & !
  ztheta,               & !
  ztheta_l,             & !
  ztmkvm,               & !
  ztmkvh,               & !
  ztch,                 & !
  ztcm,                 & !
#ifdef COUP_OAS_COS
  ztcw,                 & !
#endif
  zkh,                  & !
  ztmkvw,               & !
  ztmkvtke,             & !
  lvertdiff_w,          & ! .TRUE. for vertical diffusion of w
  lmoist_turb,          & ! .TRUE. if moist turb. param. (Postdam) is used
  lmassf_diffusion        !

USE src_slow_tendencies_rk, ONLY:  &
  complete_tendencies_init,        &
  complete_tendencies_uvwtpp,      &
  complete_tendencies_uvwtpp_eva,  &
  complete_tendencies_qvqcqi_tke,  &
#ifdef COSMOART
  complete_tendencies_cosmo_art,   &
#endif
#ifdef POLLEN
  complete_tendencies_pollen,      &
#endif
!AK (02.04.2012)
  complete_tendencies_tracer,      &
!AK (02.04.2012)
  complete_tend_uvwtpp_CN3Crow,    &
  explicit_horizontal_diffusion,   &
  implicit_vert_diffusion_uvwt

! coefficients used in advection
! ------------------------------

USE src_advection_rk, ONLY:  &
  lef_adv_qx_notpd,     & ! .TRUE. if Euler forward adv. scheme for qx
                          ! is NOT positive definite
  lqrqsqg_trilin,       & ! .TRUE. if trilin. interpol. is used for qr, qs, qg
  lqx_conserv_form,     & ! .TRUE. if qx transport in conservation form
  lcalrho_advprog,      & ! .TRUE. if rho is advected prognostically
  implicit_sedim,       & !
  advection,            & !
  advection_pd            !

!------------------------------------------------------------------------------

USE environment,              ONLY :  &
  exchg_boundaries,       & ! performs the boundary exchange between 
                            ! neighboring processors
  comm_barrier,           & ! 
  model_abort, collapse

!------------------------------------------------------------------------------

USE time_utilities,           ONLY :  get_timings, i_dyn_computations,      &
     i_slow_tendencies, i_horizontal_diffusion, i_horizontal_advection,     &
     i_fast_waves, i_fast_waves_barrier, i_fast_waves_comm,                 &
     i_add_tend_moisture, i_communications_dyn, i_barrier_waiting_dyn

USE meteo_utilities,          ONLY :  &
  satad                     !

USE parallel_utilities,       ONLY :  &
  i_global, j_global,     & !
  global_values             !

USE hori_diffusion,           ONLY :  &
  comp_hori_diff

USE numeric_utilities_rk,     ONLY :  &
  clipping,               & !
  init_bott_coeffs          !

USE fast_waves_rk,            ONLY :  &
  fast_waves_runge_kutta    ! fast waves solver for Runge-Kutta core

!------------------------------------------------------------------------------

#ifdef NUDGING
USE data_lheat_nudge,    ONLY :  &
    llhn,         & ! main switch for latent heat nudging (lhn)
    llhnverif,    & ! main switch for latent heat nudging (lhn)
    tt_lheat                ! profile of t-increments due to latent heating   ( K/s )
                            ! (stored for current and previous timestep)

!------------------------------------------------------------------------------

USE src_lheating,             ONLY :  &
    get_gs_lheating         ! storage of grid scale latent heating for lhn
#endif

!------------------------------------------------------------------------------

! CK 20101204 ART RK interfaces
#ifdef POLLEN
USE data_pollen,         ONLY :     &
    isp_pollen   , & ! number of pollen species
    cpollen      , & ! pollen concentration                          (# m-3 )
    cpollentens  , & ! cpolle-tendency without sound-wave terms     (# m-3/s)
    cpollent_conv, & !
    vsed_p           ! sedimentation velocity aerosols                  (m/s)
#endif

#ifdef COSMOART
USE data_cosmo_art,      ONLY :     &
    lgas         , & ! of gases
    laero        , & ! of aerosols
    isp_gastrans , & ! number of transported gas phase species
    isp_aerotrans, & ! number of transported aerosol variables
    cgas         , & ! gas phase concentration                        ( ppm )
    caero        , & ! aerosol concentration                        (ug m-3 )
    cgastens     , & ! cgas-tendency without sound-wave terms         (ppm/s)
    cgast_conv   , & !
    caerotens    , & ! caero-tendency without sound-wave terms     (ug m-3/s)
    caerot_conv  , &
    vd           , & ! deposition velocity gas phase                    (m/s)
    vdepa        , & ! deposition velocity aerosols                     (m/s)
    vseda            ! sedimentation velocity aerosols                  (m/s)
#endif

!==============================================================================

IMPLICIT NONE

!==============================================================================

!------------------------------------------------------------------------------
! Declarations
!------------------------------------------------------------------------------

REAL (KIND = ireals)              :: &
  alpha_rk,       & ! 1. parameter for RK-scheme
  beta_rk,        & ! 2. parameter for RK-scheme
  gamma_rk,       & ! 3. parameter for RK-scheme
  dts,            & ! small time step for time splitting
  dtsmall,        & ! small time step in last step
  dtsmall_1,      & ! small time step first intermediate step
  dtsmall_2,      & ! small time step first intermediate step
  cfl_eva,        & ! CFL-value (for explicit vertical advection)
  epsray            ! Rayleigh damping coefficient for vhmx_vol > vhmx_cfl

INTEGER (KIND=iintegers)          :: &
  ismtstep,       & ! number of small time steps in last step
  ismtstep_1,     & ! number of small time steps in first intermediate step
  ismtstep_2,     &   ! number of small time steps in second intermediate step
!AK (30.03.2012)
  iig               ! loop index for ntracer
!AK (30.03.2012)

INTEGER (KIND=iintegers)          :: &
  j2dim,          & ! middle j index for 2D-runs
  ny_2dim,        & ! number of gridpoints in y-direction for 2D-runs
  im, ip,         & ! index-boundaries for advection stencil
  imipmx            ! maximum of the stencil for explicit vertical advection

!==============================================================================

CONTAINS

!==============================================================================
!+ Module procedure in "src_runge_kutta" for organizing the time stepping  
!------------------------------------------------------------------------------

SUBROUTINE org_runge_kutta

!------------------------------------------------------------------------------
!
! Description:
!   This module procedure of module "src_runge_kutta" is the driving routine 
!   of the module, i.e. it acts as interface to the organizing routine
!   organize_dynamics
!
! Method:
!
!------------------------------------------------------------------------------

! Local scalars:
! -------------
INTEGER (KIND=iintegers) ::  &
  i,j,k,irk, isp,      & !  Loop indices
  kzdims(24),          & !  vertical dimensions for exchg_boundaries
  kitpro,              & !  number of iterations in the saturation adjustment
  izstata,             & !  error status at allocation
  izstatd,             & !  error status at deallocation
  izerror,             & !
  izs,                 & ! sum of error status
  ismts,               & ! actual number of small time steps
  nadv,                & ! actual time level to use in advection
  izdebug

REAL    (KIND=ireals   ) ::  &
  zdtr

CHARACTER (LEN=80)       ::  &
  yzerrmsg

LOGICAL :: lapply_Rayleigh_damping

! Local allocatable arrays:
! -------------------------
REAL (KIND = ireals), ALLOCATABLE :: &  
  u_half(:,:,:)  ,& ! 1. velocity component after 
                    ! half of the small time steps
  v_half(:,:,:)  ,& ! 2. velocity component ...
  w_half(:,:,:)  ,& ! 3. velocity component ...
  epsray(:,:,:)     ! Rayleigh-damping-coefficient

#ifdef COSMOART
REAL (KIND = ireals), ALLOCATABLE :: &   !  COSMO_ART
  zcaero(:,:,:,:)
#endif

#ifdef POLLEN
REAL (KIND = ireals), ALLOCATABLE :: &   !  COSMO_ART
  zcpollen(:,:,:,:)
#endif

! Local (automatic) arrays:
! -------------------------
REAL    (KIND=ireals   ) ::  &
  zwmax(ke)            ! maximum vertical velocity in a horizontal layer

! for fast_waves
REAL    (KIND=ireals   ) ::  &
  zubdt_west (je,ke),& ! u-boundary tendency at west boundary (total domain)
  zubdt_east (je,ke),& ! u-boundary tendency at east boundary (total domain)
  zvbdt_south(ie,ke),& ! v-boundary tendency at south boundary (total domain)
  zvbdt_north(ie,ke)   ! v-boundary tendency at north boundary (total domain)

LOGICAL :: l_Coriolis_every_RK_substep

! For potential temperature advection
REAL    (KIND=ireals   ) ::  &
  rovcp, rp00, govcp, exner(ie,je,ke)

! UB>>
! For debug output of MAX(U) and MIN/MAX(W) every timestep to stdout:
REAL    (KIND=ireals   ) ::  &
     vhmx_loc, wmx_loc, wmn_loc
! UB<<

! End of header
!==============================================================================

!------------------------------------------------------------------------------
! Begin Subroutine org_runge_kutta
!------------------------------------------------------------------------------

! Initialize, whether debug output shall be done
  IF (ldebug_dyn) THEN
    IF (lprintdeb_all) THEN
      izdebug = idbg_level
    ELSE
      IF (my_cart_id == 0) THEN
        izdebug = idbg_level
      ELSE
        izdebug = 0
      ENDIF
    ENDIF
  ELSE
    izdebug = 0
  ENDIF

! l_Coriolis_every_RK_substep = .TRUE.      ! stable variant
  l_Coriolis_every_RK_substep = .FALSE.

  SELECT CASE( TRIM(y_scalar_advect) )
  CASE( "PPM", "PPM_STRANG" )
    lef_adv_qx_notpd = .TRUE.
  CASE DEFAULT
    lef_adv_qx_notpd = .FALSE.
  END SELECT

  ! Tri-linear or tri-cubic interpolation for qr, qs and qg
  lqrqsqg_trilin = .FALSE.

  ! Transport of qx in conservation form on or off
  lqx_conserv_form = .TRUE.

  ! Compute advection of rho prognostically simultaneously to qx transport
  lcalrho_advprog = .TRUE.

  ! Vertical diffusion of vertical velocity on or off
  lvertdiff_w = .TRUE.

  ! Vertical mass redistribution scheme for cloud water on or off
  lmassf_diffusion = .TRUE. ! after calculation of the diffusion

  ! Include diabatic forcing due to latent heat of saturation adjustment
  ldiabf_satad = .FALSE.


  ! Prog. advection of rho is only an option if conservation form for qx
  ! is choosen.
  IF ( .NOT.lqx_conserv_form ) lcalrho_advprog = .FALSE.

  ! Set ldiabf_satad according to NAMELIST parameter ldiabf_lh
  IF ( .NOT.ldiabf_lh ) ldiabf_satad = .FALSE.

  ! Set constants for potential temperature advection
  rovcp = r_d/cp_d
  rp00  = 1.E-5_ireals
  govcp = g/cp_d

!------------------------------------------------------------------------------
! Section 1: Some preparations at the initial time step(s):
!------------------------------------------------------------------------------

  IF ( ntstep == nstart ) THEN

    SELECT CASE( TRIM(y_scalar_advect) )
    CASE( "BOTT2", "BOTT4", "BOTT2_STRANG", "BOTT4_STRANG" )
      CALL init_bott_coeffs
    END SELECT

    CALL calc_small_timestep( dt, ismtstep, ismtstep_1, ismtstep_2,         &
                              dtsmall, dtsmall_1, dtsmall_2 )

    ! some preparations to simplify the choice of
    ! different advection shemes 
    ! setting of the indices of the advected stencil
    SELECT CASE(iadv_order)
    CASE(1)
      im = -1
      ip = 1
    CASE(2)
      im = -1
      ip = 1
    CASE(3)
      im = -2
      ip = 2
    CASE(4)
      im = -2
      ip = 2
    CASE(5)
      im = -3
      ip = 3
    CASE(6)
      im = -3
      ip = 3
    END SELECT
    IF ( y_vert_adv_dyn == "expl" ) THEN
      cfl_eva = 1.0_ireals
      SELECT CASE(ieva_order)
      CASE(1)
        imipmx = 1
      CASE(2)
        imipmx = 1
      CASE(3)
        imipmx = 2
        IF ( irk_order == 3 ) cfl_eva = 1.4_ireals
      CASE(4)
        imipmx = 2
        cfl_eva = 0.9_ireals
      CASE(5)
        imipmx = 3
        cfl_eva = 0.95_ireals
      CASE(6)
        imipmx = 3
        cfl_eva = 0.75_ireals
      END SELECT
    END IF

    j2dim = nboundlines + 1
    ny_2dim = nboundlines + j2dim

    IF (my_cart_id == 0) THEN
      PRINT *,' LEVELINDEX WHERE LEVELS BECOME FLAT, KFLAT = ', kflat
      PRINT *,' BETA_SW  = ', betasw,  ' BETA_GW  = ', betagw
      PRINT *,' BETA2_SW = ', beta2sw, ' BETA2_GW = ', beta2gw
      PRINT *,' XKD = ', xkd
      IF (lrubc) THEN
        PRINT *,' CAUTION!!: RADIATIVE UPPER BOUNDARY CONDITION NOT &
          &IMPLEMENTED'
        PRINT *,' CAUTION!!: MODEL RUN USES RIGID UPPER LID'
      ENDIF
    ENDIF

    ! calculate reference temperature t0 and gradient of t0
    IF ((itheta_adv == 0) .AND. (irefatm == 1)) THEN
      DO  k = 1, ke
        t0(:,:,k)    = p0(:,:,k) / ( r_d * rho0(:,:,k) )
        dt0dz(:,:,k) = - dt0lp * g * rho0(:,:,k) / p0(:,:,k)
      ENDDO
    ELSE IF ((itheta_adv == 0) .AND. (irefatm == 2)) THEN
      DO  k = 1, ke
        t0(:,:,k)    = p0(:,:,k) / ( r_d * rho0(:,:,k) )
        dt0dz(:,:,k) =  -delta_t/h_scal*EXP(-0.5_ireals*(hhl(:,:,k)+ &
                         hhl(:,:,k+1))/h_scal)
      ENDDO
    ELSE IF ((itheta_adv == 0) .AND. (irefatm == 3)) THEN
      DO  k = 1, ke
        t0(:,:,k)    = p0(:,:,k) / ( r_d * rho0(:,:,k) )
        dt0dz(:,:,k) = t0(:,:,k)*bvref**2/g - govcp
      ENDDO
    ! when potential temperature is used, t0 and dt0dz carry
    ! theta_0 and dtheta_0/dz
    ELSE IF ((itheta_adv >= 1) .AND. (irefatm == 1)) THEN
      DO  k = 1, ke
        t0(:,:,k)    = p0(:,:,k) / ( r_d * rho0(:,:,k) )* &
          EXP(-rovcp*LOG(p0(:,:,k)*rp00))
        dt0dz(:,:,k) = (govcp - dt0lp * g * rho0(:,:,k) / p0(:,:,k))* &
          EXP(-rovcp*LOG(p0(:,:,k)*rp00))
      ENDDO
    ELSE IF ((itheta_adv >= 1) .AND. (irefatm == 2)) THEN
      DO  k = 1, ke
        t0(:,:,k)    = p0(:,:,k) / ( r_d * rho0(:,:,k) )* &
          EXP(-rovcp*LOG(p0(:,:,k)*rp00))
        dt0dz(:,:,k) =  (govcp - delta_t/h_scal*EXP(-0.5_ireals*(hhl(:,:,k)+ &
          hhl(:,:,k+1))/h_scal)) * EXP(-rovcp*LOG(p0(:,:,k)*rp00))
      ENDDO
    ELSE IF ((itheta_adv >= 1) .AND. (irefatm == 3)) THEN
      DO  k = 1, ke
        t0(:,:,k)    = p0(:,:,k) / ( r_d * rho0(:,:,k) )* &
          EXP(-rovcp*LOG(p0(:,:,k)*rp00))
        dt0dz(:,:,k) = t0(:,:,k)*bvref**2/g
      ENDDO
    ENDIF

  ENDIF

! UB>>
! Determine the min/.max. W in the model domain and the max |vh| and write to stdout
! for diagnostic purposes:
  IF (ldebug_dyn .AND. idbg_level > 3) THEN
    vhmx_loc = SQRT(MAXVAL( u(:,:,:,nnow)**2 + v(:,:,:,nnow)**2 ))
    wmx_loc  = MAXVAL(w(:,:,:,nnow))
    wmn_loc  = MINVAL(w(:,:,:,nnow))
    IF (num_compute > 1) THEN
      CALL global_values (vhmx_loc, 1, 'MAX', imp_reals, icomm_cart, -1,    &
           yzerrmsg, izerror)
      CALL global_values (wmx_loc, 1, 'MAX', imp_reals, icomm_cart, -1,    &
           yzerrmsg, izerror)
      CALL global_values (wmn_loc, 1, 'MIN', imp_reals, icomm_cart, -1,    &
           yzerrmsg, izerror)
    ENDIF
    IF ( my_cart_id == 0 ) THEN
      WRITE (*,'(a,es12.5)') 'Maximum V_h : ', vhmx_loc
      WRITE (*,'(a,es12.5,a,es12.5)') 'Maximum W : ', wmx_loc, '     Minimum W : ', wmn_loc
    END IF
  END IF
! UB<<

  ALLOCATE( epsray(1:ie, 1:je, 1:ke), STAT=izstata )
  IF ( izstata /= 0 ) THEN
    yzerrmsg="allocation of epsray"
    CALL model_abort (my_cart_id, 100+izerror, yzerrmsg, 'org_runge_kutta')
  END IF

  CALL check_cfl_horiz_advection ( lapply_Rayleigh_damping, epsray )

!------------------------------------------------------------------------------
! Section 2: Allocate and store boundary tendencies for u and v in case of
!            non-periodic boundary conditions.
!------------------------------------------------------------------------------

  ! These tendencies are required to compute the wind divergence along the
  ! boundary lines in the fast-wave solver.
  ! The boundary tendencies are premultiplied by the small time step
  ! The nnew-values are the boundary values provided in the beginning of the
  ! time step.

  zubdt_west(:,:) = 0.0_ireals
  zubdt_east(:,:) = 0.0_ireals
  zvbdt_south(:,:) = 0.0_ireals
  zvbdt_north(:,:) = 0.0_ireals

  IF ( .NOT.lperi_x ) THEN
    IF (my_cart_neigh(1) == -1) THEN
      DO k = 1, ke
        DO j = jstartu, jendu
          zubdt_west(j,k) = (u(istartu-1,j,k,nnew)-u(istartu-1,j,k,nnow))/dt
        ENDDO
      ENDDO
    ENDIF
    IF (my_cart_neigh(3) == -1) THEN
      DO k = 1, ke
        DO j = jstartu, jendu
          zubdt_east(j,k) = (u(iendu+1,j,k,nnew) - u(iendu+1,j,k,nnow))/dt
        ENDDO
      ENDDO
    ENDIF
  END IF
  IF ( .NOT.lperi_y ) THEN
    IF (my_cart_neigh(4) == -1) THEN
      DO k = 1, ke
        DO i = istartv, iendv
          zvbdt_south(i,k) = (v(i,jstartv-1,k,nnew)-v(i,jstartv-1,k,nnow))/dt
        ENDDO
      ENDDO
    ENDIF
    IF (my_cart_neigh(2) == -1) THEN
      DO k = 1, ke
        DO i = istartv, iendv
          zvbdt_north(i,k) = (v(i,jendv+1,k,nnew) - v(i,jendv+1,k,nnow))/dt
        ENDDO
      ENDDO
    ENDIF
  ENDIF


!------------------------------------------------------------------------------
! Section 3: Do forecast using the two time level scheme
!------------------------------------------------------------------------------

  ! Allocation of the separate tendency arrays for advection
  izs = 0
  ALLOCATE ( uadvt (ie,je,ke),  STAT=izstata ); izs = izs + izstata
  ALLOCATE ( vadvt (ie,je,ke),  STAT=izstata ); izs = izs + izstata
  ALLOCATE ( wadvt (ie,je,ke1), STAT=izstata ); izs = izs + izstata
  ALLOCATE ( ppadvt(ie,je,ke),  STAT=izstata ); izs = izs + izstata
  ALLOCATE ( tadvt (ie,je,ke),  STAT=izstata ); izs = izs + izstata

  IF ( izs /= 0 ) THEN
    yzerrmsg="allocation of xadvt"
    CALL model_abort (my_cart_id, 100+izerror, yzerrmsg, 'org_runge_kutta')
  END IF

  uadvt (:,:,:) = 0.0_ireals
  vadvt (:,:,:) = 0.0_ireals
  wadvt (:,:,:) = 0.0_ireals
  ppadvt(:,:,:) = 0.0_ireals
  tadvt (:,:,:) = 0.0_ireals

  ALLOCATE ( wcon (ie,je,ke1), STAT=izstata )
  wcon(:,:,:) = 0.0_ireals

  IF (ltime) CALL get_timings (i_dyn_computations, ntstep, dt, izerror)

!------------------------------------------------------------------------------
! Section 3a: Complete the slow tendencies from physical process by some
!             adiabatic forcings (Coriolis force, water loading etc.)
!             and add them to the advective tendencies.
!------------------------------------------------------------------------------

  ! Completion of temperature tendency due to radiative forcing.
  ! Completion of tendencies due to convective forcings.
  ! Add Rayleigh friction to the velocity field.
  ! ---------------------------------------------------------------------------

  ! Set the time level for Rayleigh friction and add damping terms

  IF ( lapply_Rayleigh_damping ) THEN
    DO  k = 1 , ke
      DO   j = jstartu, jendu
        DO i = istartu, iendu
          utens(i,j,k) = utens(i,j,k) - epsray(i,j,k) * u(i,j,k,nnow)
        ENDDO
      ENDDO
      DO   j = jstartv, jendv
        DO i = istartv, iendv
          vtens(i,j,k) = vtens(i,j,k) - epsray(i,j,k) * v(i,j,k,nnow)
        ENDDO
      ENDDO
    ENDDO
  END IF

  DEALLOCATE( epsray, STAT=izstatd )
  IF ( izstatd /= 0 ) THEN
    yzerrmsg="deallocation of epsray"
    CALL model_abort (my_cart_id, 100+izerror, yzerrmsg, 'org_runge_kutta')
  END IF

  ! Add physical forcings if required
  IF (lphys) THEN
    DO  k = 1 , ke
      DO    j = jstart, jend 
        DO  i = istart, iend
          ttens (i,j,k) = ttens (i,j,k) + sohr(i,j,k) + thhr(i,j,k)
          ttens (i,j,k) = ttens (i,j,k) + tt_conv (i,j,k)
          ! ... substract explicit diffusion part from qv-tendency
          !     for correct consistent calculation of dqvdt
          qvtens(i,j,k) = qvtens(i,j,k) + qvt_conv(i,j,k) - qvt_diff(i,j,k)
!AK (30.03.2012)
          DO iig=1, ntracer
            tracertens (i,j,k,iig) = tracertens (i,j,k,iig) + tracert_conv (i,j,k,iig)
          ENDDO
!AK (30.03.2012)
        ENDDO
      ENDDO

      IF (lsso) THEN
        DO    j = jstart, jend
          DO  i = istart, iend
            ttens(i,j,k) = ttens(i,j,k) + tt_sso(i,j,k)
          ENDDO
        ENDDO
      ENDIF

      IF (itype_conv == 1) THEN
        DO    j = jstart, jend
          DO  i = istart, iend
            qctens(i,j,k) = qctens(i,j,k) + qct_conv(i,j,k)
          ENDDO
        ENDDO
      ENDIF

#ifndef HYMACS
      IF ( (itype_conv == 0) .OR. (itype_conv == 2) ) THEN
#endif
#ifdef HYMACS
! VK 2012/03/15
      IF ( (itype_conv == 0) .OR. (itype_conv == 2) .OR. (itype_conv == 4)) THEN
#endif
        ! Add convective tendencies of qc and qi computed by Tiedtke(0) or
        ! Bechtold Scheme(2)
        DO    j = jstart, jend
          DO  i = istart, iend
            qctens(i,j,k) = qctens(i,j,k) + qct_conv(i,j,k)
            qitens(i,j,k) = qitens(i,j,k) + qit_conv(i,j,k)
          ENDDO
        ENDDO
      ENDIF

#ifdef HYMACS
! VK 2012/03/15
      IF (itype_conv == 4) THEN
        DO    j = jstart, jend
          DO  i = istart, iend
            pptens(i,j,k) = pptens(i,j,k) + hymacs_ppt_conv(i,j,k)
          ENDDO
        ENDDO
      ENDIF
#endif

#ifdef COSMOART
      ! CK 20101204 ART RK interfaces
      IF (l_cosmo_art) THEN
        IF ( (itype_conv == 0) ) THEN

          IF (lgas) THEN
            DO isp = 1 , isp_gastrans
              DO j = jstart, jend
                DO i = istart, iend
                  cgastens(i,j,k,isp) = cgastens(i,j,k,isp)   &
                                      + cgast_conv(i,j,k,isp)
                ENDDO
              ENDDO
            ENDDO
          ENDIF

          IF (laero) THEN
            DO isp = 1 , isp_aerotrans
              DO j = jstart, jend
                DO i = istart, iend
                  caerotens(i,j,k,isp) = caerotens(i,j,k,isp)   &
                                       + caerot_conv(i,j,k,isp)
                ENDDO
              ENDDO
            ENDDO
          ENDIF
        ENDIF
      ENDIF
#endif

#ifdef POLLEN
      IF (l_pollen) THEN
        IF ( (itype_conv == 0) ) THEN

          DO isp = 1 , isp_pollen
            DO j = jstart, jend
              DO i = istart, iend
                cpollentens(i,j,k,isp) = cpollentens(i,j,k,isp)   &
                                       + cpollent_conv(i,j,k,isp)
              ENDDO
            ENDDO
          ENDDO

        ENDIF
      ENDIF
#endif

      DO   j = jstartu, jendu
        DO i = istartu, iendu
          utens(i,j,k) = utens(i,j,k) + ut_conv(i,j,k)
        ENDDO
      ENDDO
      IF (lsso) THEN
        DO   j = jstartu, jendu
          DO i = istartu, iendu
            utens(i,j,k) = utens(i,j,k) + ut_sso(i,j,k)
          ENDDO
        ENDDO
      ENDIF
      DO   j = jstartv, jendv
        DO i = istartv, iendv
          vtens(i,j,k) = vtens(i,j,k) + vt_conv(i,j,k)
        ENDDO
      ENDDO
      IF (lsso) THEN
        DO   j = jstartv, jendv
          DO i = istartv, iendv
            vtens(i,j,k) = vtens(i,j,k) + vt_sso(i,j,k)
          ENDDO
        ENDDO
      ENDIF
    ENDDO

!MU (19.09.2012)
#ifdef COUP_OAS_COS
    DO iig=1, ntracer
      IF ( ltracer(7,iig) == 1) THEN
        DO    j = jstart, jend 
          DO  i = istart, iend
            tracertens (i,j,ke,iig) = tracertens (i,j,ke,iig) + co2fl_s(i,j)
!MU (12.04.13)
!            tracertens (i,j,ke,iig) = tracertens (i,j,ke,iig) + psn_tens(i,j) + plres_tens(i,j)
!MU (12.04.13)
          ENDDO
        ENDDO
      ENDIF
    ENDDO
#endif
!MU (19.09.2012)

    ! add estimated temperature tendency due to latent heat
    ! release in previous time step
    IF ( ldiabf_lh ) THEN
      zdtr = 1.0_ireals / dt
      DO k = 1, ke
        DO j = jstart, jend 
          DO i = istart, iend
            ttens(i,j,k) = ttens(i,j,k) + zdtr * tinc_lh(i,j,k)
          ENDDO
        ENDDO
      ENDDO
    ENDIF
  ENDIF

  IF ( lcori .AND. (.NOT. l_Coriolis_every_RK_substep) ) THEN
    CALL coriolis( utens, vtens, wtens, nnow ) 
  END IF

  ! precalculate some coefficients for complete_tendencies_...
  ! ----------------------------------------------------------------

  ALLOCATE ( zsqrtgrho_r_s(ie,je,ke), zsqrtgrho_r_u(ie,je,ke),           &
    zsqrtgrho_r_v(ie,je,ke), STAT=izstata )
  IF ( izstata /= 0 ) THEN
    yzerrmsg="allocation of zsqrtgrho_r_s, ..."
    CALL model_abort (my_cart_id, 100+izerror, yzerrmsg, 'org_runge_kutta')
  END IF

#ifdef COUP_OAS_COS
  ALLOCATE ( za1t(ke1), za2t(ke1), zpia(ie,je,ke), zpianf(ie,je,ke1),    &
    ztmkvm(ie,je,ke), ztmkvh(ie,je), ztch(ie,je), ztcm(ie,je),           &
    ztcw(ie,je), ztheta(ie,je,ke), zkh(ie,je,ke), STAT=izstata )
#else
  ALLOCATE ( za1t(ke1), za2t(ke1), zpia(ie,je,ke), zpianf(ie,je,ke1),    &
    ztmkvm(ie,je,ke), ztmkvh(ie,je), ztch(ie,je), ztcm(ie,je),           &
    ztheta(ie,je,ke), zkh(ie,je,ke), STAT=izstata )
#endif
  IF ( izstata /= 0 ) THEN
    yzerrmsg="allocation of za1t"
    CALL model_abort (my_cart_id, 100+izerror, yzerrmsg, 'org_runge_kutta')
  END IF

  IF ( lvertdiff_w ) THEN
    ALLOCATE ( ztmkvw(ie,je,ke),                        &
      zsqrtgrho_r_w(ie,je,ke), STAT=izstata )
    IF ( izstata /= 0 ) THEN
      yzerrmsg="allocation of ztmkvw, ..."
      CALL model_abort (my_cart_id, 100+izerror, yzerrmsg, 'org_runge_kutta')
    END IF
  END IF

  IF ( lprog_tke ) THEN
    ALLOCATE ( ztmkvtke(ie,je,ke), STAT=izstata )
    IF ( izstata /= 0 ) THEN
      yzerrmsg="allocation of ztmkvtke"
      CALL model_abort (my_cart_id, 100+izerror, yzerrmsg, 'org_runge_kutta')
    END IF
  END IF

  IF ( itype_turb == 6 .OR. itype_turb == 8 ) THEN
    lmoist_turb = .TRUE.
    IF ( l3dturb ) THEN
      ALLOCATE ( ztheta_l(ie,je,ke), STAT=izstata )
      IF ( izstata /= 0 ) THEN
        yzerrmsg="allocation of ztheta_l"
        CALL model_abort (my_cart_id, 100+izerror, yzerrmsg, 'org_runge_kutta')
      END IF
    END IF
  ELSE
    lmoist_turb = .FALSE.
  END IF

  CALL complete_tendencies_init

  IF (ltime) CALL get_timings (i_slow_tendencies, ntstep, dt, izerror)

  ! Horizontal diffusion
  ! ----------------------------------------------------------------
  IF ( l3dturb ) CALL explicit_horizontal_diffusion

  ! Vertical diffusion
  ! ----------------------------------------------------------------
  CALL implicit_vert_diffusion_uvwt

  IF (ltime) CALL get_timings (i_horizontal_diffusion, ntstep, dt, izerror)

  IF (itheta_adv == 0) THEN
    ! calculate temperature perturbation tp=t-t0
    ! (temporarily override array t)
!CDIR COLLAPSE
    DO k = 1, ke
      DO j = 1, je
        DO i = 1, ie
          t(i,j,k,nnow) = t(i,j,k,nnow) - t0(i,j,k)
!          t(i,j,k,nnew) = t(i,j,k,nnew) - t0(i,j,k) ! not needed
        END DO
      END DO
    END DO
  ELSE IF (itheta_adv == 1) THEN
    ! Convert temperature into potential temperature and compute perturbation
!CDIR COLLAPSE
    DO k = 1, ke
      DO j = 1, je
        DO i = 1, ie
          exner(i,j,k) = EXP(-rovcp*LOG((p0(i,j,k)+pp(i,j,k,nnow))*rp00))
          t(i,j,k,nnow) = t(i,j,k,nnow)*exner(i,j,k) - t0(i,j,k)
          ttens(i,j,k)  = ttens(i,j,k)*exner(i,j,k)
!           t(i,j,k,nnew) = t(i,j,k,nnew)*EXP(-rovcp*LOG((p0(i,j,k)+ &
!                           pp(i,j,k,nnew))*rp00)) - t0(i,j,k) ! not needed
        END DO
      END DO
    END DO
  ELSE IF (itheta_adv == 2) THEN
    ! Convert temperature into potential temperature
!CDIR COLLAPSE
    DO k = 1, ke
      DO j = 1, je
        DO i = 1, ie
          exner(i,j,k) = EXP(-rovcp*LOG((p0(i,j,k)+pp(i,j,k,nnow))*rp00))
          t(i,j,k,nnow) = t(i,j,k,nnow)*exner(i,j,k)
          ttens(i,j,k)  = ttens(i,j,k)*exner(i,j,k)
!           t(i,j,k,nnew) = t(i,j,k,nnew)*EXP(-rovcp*LOG((p0(i,j,k)+ &
!                           pp(i,j,k,nnew))*rp00)) ! not needed
        END DO
      END DO
    END DO
  ENDIF

  !----------------------------------------------------------------------------
  ! Runge-Kutta time integration ...
  !----------------------------------------------------------------------------

  ! set time level for advection in first RK step
  nadv = nnow

  alpha_rk = 1.0_ireals
  beta_rk  = 0.0_ireals

  runge_kutta_loop: DO irk = 1, irk_order

    !--------------------------------------------------------------------------
    ! Section 3b: Compute slow tendencies by first evaluating the advective
    !             and then providing slow tendencies including physical ones.
    !--------------------------------------------------------------------------

    IF (irk_order == 3) THEN
      SELECT CASE(irunge_kutta)
      CASE(1)
        SELECT CASE(irk)
        CASE(1)
          gamma_rk = 1.0_ireals/3.0_ireals
          ismts = ismtstep_1   ! set number of small time steps
          dts   = dtsmall_1
        CASE(2)
          gamma_rk = 0.5_ireals
          ismts = ismtstep_2   ! set number of small time steps
          dts   = dtsmall_2
        CASE(3)
          gamma_rk = 1.0_ireals
          ismts = ismtstep     ! set number of small time steps
          dts   = dtsmall
        END SELECT
      CASE(2)
        SELECT CASE(irk)
        CASE(1)
          alpha_rk = 1.0_ireals
          beta_rk  = 0.0_ireals
          gamma_rk = 1.0_ireals
          ismts = ismtstep_1   ! set number of small time steps
          dts   = dtsmall_1
        CASE(2)
          alpha_rk = 0.75_ireals
          beta_rk  = 0.25_ireals
          gamma_rk = 0.25_ireals
          ismts = ismtstep_2   ! set number of small time steps
          dts   = dtsmall_2
        CASE(3)
          alpha_rk = 1.0_ireals/3.0_ireals
          beta_rk  = 2.0_ireals/3.0_ireals
          gamma_rk = 2.0_ireals/3.0_ireals
          ismts = ismtstep     ! set number of small time steps
          dts   = dtsmall
        END SELECT
      END SELECT
    ELSE
      IF (irk == irk_order) THEN
        gamma_rk = 1.0_ireals
        ismts = ismtstep     ! set number of small time steps
        dts   = dtsmall
      ELSE
        gamma_rk = 0.5_ireals
        ismts = ismtstep_1   ! set number of small time steps
        dts   = dtsmall_1
      END IF
    END IF

    !--------------------------------------------------------------------------
    ! Section 3b1: Compute the time tendencies due to horizontal advection (and Coriolis terms)
    !              and store them on separate tendency fields (uadvt,...)
    !--------------------------------------------------------------------------

    ! Compute advection for dynamic variables
    CALL advection( nadv, irk, im, ip, imipmx, cfl_eva, j2dim, ny_2dim)

    IF ( lcori .AND. l_Coriolis_every_RK_substep ) THEN
      ! Compute Coriolisterms for u, v (and perhaps w)
      CALL coriolis( uadvt, vadvt, wadvt, nadv )
    END IF

    IF (ltime) CALL get_timings (i_horizontal_advection, ntstep, dt, izerror)

    !--------------------------------------------------------------------------
    ! Section 3b2: Compute the slow tendencies
    !              and store them on separate tendency fields (uadvt,...)
    !--------------------------------------------------------------------------

    IF ( y_vert_adv_dyn == "impl2" ) THEN

      IF ( irk_order == 3 .AND. irunge_kutta == 2 ) THEN
        CALL complete_tendencies_uvwtpp( nadv, gamma_rk*dt, alpha_rk, beta_rk )
      ELSE
        CALL complete_tendencies_uvwtpp( nadv, gamma_rk*dt )
      END IF

    ELSE IF ( y_vert_adv_dyn == "expl" ) THEN

      CALL complete_tendencies_uvwtpp_eva

    ELSE IF ( y_vert_adv_dyn == "impl3" ) THEN  
       ! tendencies are not added in this implicit scheme, therefore do it here explicitely:

       uadvt(:,:,:)  = uadvt(:,:,:)  + utens(:,:,:)
       vadvt(:,:,:)  = vadvt(:,:,:)  + vtens(:,:,:)
       wadvt(:,:,:)  = wadvt(:,:,:)  + wtens(:,:,:)
       ppadvt(:,:,:) = ppadvt(:,:,:) + pptens(:,:,:)
       tadvt(:,:,:)  = tadvt(:,:,:)  + ttens(:,:,:)
    ELSE
      yzerrmsg="false value for y_vert_adv_dyn"
      CALL model_abort (my_cart_id, 100+izerror, yzerrmsg, 'org_runge_kutta')
    END IF


    IF (ltime) CALL get_timings (i_slow_tendencies, ntstep, dt, izerror)

    ! set time level for advection in remaining RK step(s)
    nadv = nnew

    !premultiply tendencies with small time step
    uadvt  = uadvt * dts
    vadvt  = vadvt * dts
    wadvt  = wadvt * dts
    tadvt  = tadvt * dts
    ppadvt = ppadvt * dts

    !--------------------------------------------------------------------------
    ! Runge-Kutta step(s) ...
    ! Section 3c: Compute the intermediate values of the dynamic fields
    !             (using the time level nnew for storage) by integrating
    !             the sound and gravity wave terms, including the slow
    !             tendencies with small time steps for the time interval
    !             (dt/4,) (dt/2,) dt with the split-explicit scheme.
    !             Exchange all boundary fields afterwards.
    !--------------------------------------------------------------------------

    CALL fast_waves_runge_kutta( ismts, dts, xkd,                          &
      uadvt, vadvt, wadvt, tadvt, ppadvt,                                  &
      zubdt_west*dts, zubdt_east*dts, zvbdt_north*dts, zvbdt_south*dts,    &
      irk, alpha_rk, beta_rk )

    IF ( irk == irk_order .AND. ldiabf_lh ) THEN
      ! subtract estimated temperature increment due to latent heat 
      ! release in previous time step
      IF (itheta_adv == 0) THEN
        t(:,:,:,nnew) = t(:,:,:,nnew) - tinc_lh(:,:,:)
      ELSE
        t(:,:,:,nnew) = t(:,:,:,nnew) - tinc_lh(:,:,:) * &
          exner(:,:,:)
      ENDIF
    ENDIF
    

    IF (ltime) THEN
      CALL get_timings (i_fast_waves, ntstep, dt, izerror)
      IF (ltime_barrier) THEN
        CALL comm_barrier (icomm_cart, izerror, yzerrmsg)
        CALL get_timings (i_fast_waves_barrier, ntstep, dt, izerror)
      ENDIF
    ENDIF

!!$    IF (num_compute > 1) THEN
      kzdims(1:20)=(/ke,ke,ke1,ke,ke,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
      CALL exchg_boundaries                                                &
        (52+nnew, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute, ie, je,     &
        kzdims, jstartpar, jendpar, nbl_exchg, nboundlines, my_cart_neigh, &
        lperi_x, lperi_y, l2dim, &
        20000+ntstep, ldatatypes, ncomm_type, izerror, yzerrmsg,           &
        u (:,:,:,nnew), v (:,:,:,nnew), w (:,:,:,nnew), t (:,:,:,nnew),    &
        pp(:,:,:,nnew))
!!$    ENDIF
!!$    IF ( l2dim ) THEN
!!$      DO j = 1, ny_2dim
!!$        u (:,j,:,nnew) = u (:,j2dim,:,nnew)
!!$        v (:,j,:,nnew) = v (:,j2dim,:,nnew)
!!$        w (:,j,:,nnew) = w (:,j2dim,:,nnew)
!!$        t (:,j,:,nnew) = t (:,j2dim,:,nnew)
!!$        pp(:,j,:,nnew) = pp(:,j2dim,:,nnew)
!!$      ENDDO
!!$    ENDIF

    IF (ltime) CALL get_timings (i_fast_waves_comm, ntstep, dt, izerror)

  END DO runge_kutta_loop


  IF ( y_vert_adv_dyn == "impl3" ) THEN
    ! complete operatorsplitting between hor-Adv./RK  and vertical advection

    CALL calc_wcon_sqrtg( u(:,:,:,nnew), v(:,:,:,nnew), w(:,:,:,nnew), wcon(:,:,:) )
    ! attention: now  wcon  contains  'sqrtg * contravar. vertical velocity' !
 
    CALL complete_tend_uvwtpp_CN3Crow( dt, nnew )

    u (:,:,:,nnew) = u (:,:,:,nnew) + dt * uadvt (:,:,:)
    v (:,:,:,nnew) = v (:,:,:,nnew) + dt * vadvt (:,:,:)
    w (:,:,:,nnew) = w (:,:,:,nnew) + dt * wadvt (:,:,:)
    pp(:,:,:,nnew) = pp(:,:,:,nnew) + dt * ppadvt(:,:,:)
    t (:,:,:,nnew) = t (:,:,:,nnew) + dt * tadvt (:,:,:)

  END IF


  IF (itheta_adv == 0) THEN
    ! calculate absolute temperature t=tp+t0
!CDIR COLLAPSE
    DO k = 1, ke
      DO j = 1, je
        DO i = 1, ie
          t(i,j,k,nnow) = t(i,j,k,nnow) + t0(i,j,k)
          t(i,j,k,nnew) = t(i,j,k,nnew) + t0(i,j,k)
        END DO
      END DO
    END DO
  ELSE IF (itheta_adv == 1) THEN
    ! convert theta perturbation back to temperature
!CDIR COLLAPSE
    DO k = 1, ke
      DO j = 1, je
        DO i = 1, ie
          t(i,j,k,nnow) = (t(i,j,k,nnow) + t0(i,j,k)) / exner(i,j,k)
          t(i,j,k,nnew) = (t(i,j,k,nnew) + t0(i,j,k)) * &
                          EXP(rovcp*LOG((p0(i,j,k)+pp(i,j,k,nnew))*rp00))
        END DO
      END DO
    END DO
  ELSE IF (itheta_adv == 2) THEN
    ! convert theta back to temperature
!CDIR COLLAPSE
    DO k = 1, ke
      DO j = 1, je
        DO i = 1, ie
          t(i,j,k,nnow) = t(i,j,k,nnow) / exner(i,j,k)
          t(i,j,k,nnew) = t(i,j,k,nnew) * &
                          EXP(rovcp*LOG((p0(i,j,k)+pp(i,j,k,nnew))*rp00))
        END DO
      END DO
    END DO
  ENDIF

  ! Deallocate the private arrays for the advective tendencies
  izs = 0
  DEALLOCATE ( uadvt , STAT=izstatd ); izs = izs + izstatd
  DEALLOCATE ( vadvt , STAT=izstatd ); izs = izs + izstatd
  DEALLOCATE ( wadvt , STAT=izstatd ); izs = izs + izstatd
  DEALLOCATE ( ppadvt, STAT=izstatd ); izs = izs + izstatd
  DEALLOCATE ( tadvt , STAT=izstatd ); izs = izs + izstatd
  IF ( izstatd /= 0 ) THEN
    yzerrmsg="deallocation of uadvt, ..."
    CALL model_abort (my_cart_id, 100+izerror, yzerrmsg, 'org_runge_kutta')
  END IF

  !----------------------------------------------------------------------------
  ! Section 4a: Advection of scalars: qv, qc, qi, qr, qs, qg and tke
  !----------------------------------------------------------------------------

  ! Allocate working arrays
  ALLOCATE ( u_half(ie,je,ke), v_half(ie,je,ke), w_half(ie,je,ke1), &
             STAT=izstata )
  IF ( izstata /= 0 ) THEN
    yzerrmsg="allocation of u_half, ..."
    CALL model_abort (my_cart_id, 100+izerror, yzerrmsg, 'org_runge_kutta')
  END IF

  ! Calculate estimates of the velocities for time n+1/2 
  ! -------------------------------------------------------------------------
  ! Integrate the velocities
!CDIR COLLAPSE
  u_half(:,:,:)  = 0.5_ireals*( u(:,:,:,nnow)+u(:,:,:,nnew) )
!CDIR COLLAPSE
  v_half(:,:,:)  = 0.5_ireals*( v(:,:,:,nnow)+v(:,:,:,nnew) )
!CDIR COLLAPSE
  w_half(:,:,:)  = 0.5_ireals*( w(:,:,:,nnow)+w(:,:,:,nnew) )

  !
  ! Perform positive-definite advection of moisture variables + TKE
  !
!AK (30.03.2012)
    ! to avoid allocatable variables for tracer in the advection scheme and to
    ! guarantee that these variables have the default value 1 if ntracer_adv = 0
  ntracer_dim = MAX(ntracer_adv,1)
!AK (30.03.2012)
  CALL advection_pd( u_half(:,:,:), v_half(:,:,:), w_half(:,:,:), nnow, dt, &
                     im, ip, j2dim, ny_2dim)

  ! Deallocate working arrays
  DEALLOCATE ( u_half, v_half, w_half, STAT=izstatd )
  IF ( izstatd /= 0 ) THEN
    yzerrmsg="deallocation of u_half, ..."
    CALL model_abort (my_cart_id, 100+izerror, yzerrmsg, 'org_runge_kutta')
  END IF

  IF (ltime) CALL get_timings (i_add_tend_moisture, ntstep, dt, izerror)

#ifdef COSMOART
  ! CK 20101204 ART RK interfaces
  ! Sedimentation of aerosols
  IF (l_cosmo_art) THEN
    IF (laero) THEN
      ALLOCATE ( zcaero(ie,je,ke,isp_aerotrans),STAT=izstata )
      DO isp=1,isp_aerotrans

        DO k = 1, ke
          DO j = 1, je
            DO i = 1, ie
              zcaero(i,j,k,isp)=caero(i,j,k,isp,nnew)
              vseda(i,j,k,isp) = (-1.)*vseda(i,j,k,isp)
            ENDDO
          ENDDO
        ENDDO
        CALL implicit_sedim( caero(:,:,:,isp,nnew), zcaero(:,:,:,isp), &
                             vseda(:,:,:,isp), vseda(:,:,:,isp))
      ENDDO
      DEALLOCATE ( zcaero )
    ENDIF
  ENDIF
#endif

#ifdef POLLEN
  ! Sedimentation of pollen
  IF (l_pollen) THEN
    ALLOCATE ( zcpollen(ie,je,ke,isp_pollen),STAT=izstata )
    zcpollen(:,:,:,:)=cpollen(:,:,:,:,nnew)
    DO isp=1,isp_pollen
      DO k = 1,ke
        DO j = 1,je
          DO i = 1,ie
            vsed_p(i,j,k,isp) = -1.0_ireals * vsed_p(i,j,k,isp)
          ENDDO
        ENDDO
      ENDDO
      CALL implicit_sedim( cpollen(:,:,:,isp,nnew), zcpollen(:,:,:,isp), &
                           vsed_p(:,:,:,isp), vsed_p(:,:,:,isp) )
    ENDDO
    DEALLOCATE ( zcpollen )
  ENDIF
#endif

  !----------------------------------------------------------------------------
  ! Section 4b: Horizontal diffusion and tendencies of qv, qc, qi, tke and tracer (AK)
  !----------------------------------------------------------------------------

#ifdef COSMOART
  ! CK 20101117 there have been calls to comp_hori_diff here,
  ! but in the current version the diffusion and tendencies are
  ! calculated directly in the comp_hori_diff call (20 lines below)
  ! so I assume this is unnecessary...
#endif

  IF (ltime) CALL get_timings (i_horizontal_diffusion, ntstep, dt, izerror)

  ! Complete timestep, if required by vertical diffusion
  CALL complete_tendencies_qvqcqi_tke

#ifdef COSMOART
  IF (l_cosmo_art) CALL complete_tendencies_cosmo_art
#endif

#ifdef POLLEN
  IF (l_pollen) CALL complete_tendencies_pollen
#endif
!AK (02.04.2012)
  CALL complete_tendencies_tracer
!AK (02.04.2012)
  DEALLOCATE ( wcon, STAT=izstatd )
  DEALLOCATE ( zsqrtgrho_r_s, zsqrtgrho_r_u, zsqrtgrho_r_v, STAT=izstatd )
#ifdef COUP_OAS_COS
  DEALLOCATE ( za1t, za2t, zpia, zpianf, ztmkvm, ztmkvh, ztch, ztcm, zkh,  &
               ztcw, ztheta, STAT=izstatd )
#else
  DEALLOCATE ( za1t, za2t, zpia, zpianf, ztmkvm, ztmkvh, ztch, ztcm, zkh,  &
               ztheta, STAT=izstatd )
#endif
  IF ( lvertdiff_w )  DEALLOCATE ( ztmkvw, zsqrtgrho_r_w, STAT=izstatd )
  IF ( lprog_tke )    DEALLOCATE ( ztmkvtke, STAT=izstatd )
  IF ( lmoist_turb .AND. l3dturb )  DEALLOCATE ( ztheta_l, STAT=izstatd )

  IF (ltime) CALL get_timings (i_add_tend_moisture, ntstep, dt, izerror)

  !----------------------------------------------------------------------------
  ! Section 5: Horizontal diffusion at end of time step
  !----------------------------------------------------------------------------

  IF ( lhordiff ) THEN

!MU (16.04.2012)
    DO iig = 1, ntracer
      IF ((my_cart_id==1) .AND. (itype_hdiff==1) .AND. (ltracer(2,iig)==1) .AND. (ntstep<2)) THEN
        PRINT *, 'WARNING: No horizontal diffusion of tracer if itype_hdiff == 1!'
      ENDIF
    ENDDO
!MU (16.04.2012)
    
    CALL comp_hori_diff( itype_hdiff )
    IF ( hd_corr_q_bd /= 0.0_ireals .OR. hd_corr_q_in /= 0.0_ireals ) THEN
!!$      IF (num_compute > 1) THEN
        kzdims(1:24)=(/ke,ke,ke1,ke,ke,ke,ke,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
        CALL exchg_boundaries                                                 &
          (54+nnew, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute, &
           ie, je, kzdims, jstartpar, jendpar, 2, nboundlines, my_cart_neigh, &
          lperi_x, lperi_y, l2dim,                                            &
          20000+ntstep, ldatatypes, ncomm_type, izerror, yzerrmsg,            &
          u (:,:,:,nnew), v (:,:,:,nnew), w (:,:,:,nnew), t (:,:,:,nnew),     &
          qv(:,:,:,nnew), qc(:,:,:,nnew), pp(:,:,:,nnew))
!!$      ENDIF
!!$      IF ( l2dim ) THEN
!!$        DO j = 1, ny_2dim
!!$          u (:,j,:,nnew) = u (:,j2dim,:,nnew)
!!$          v (:,j,:,nnew) = v (:,j2dim,:,nnew)
!!$          w (:,j,:,nnew) = w (:,j2dim,:,nnew)
!!$          t (:,j,:,nnew) = t (:,j2dim,:,nnew)
!!$          pp(:,j,:,nnew) = pp(:,j2dim,:,nnew)
!!$          qv(:,j,:,nnew) = qv(:,j2dim,:,nnew)
!!$          qc(:,j,:,nnew) = qc(:,j2dim,:,nnew)
!!$        ENDDO
!!$      ENDIF
      CALL clipping( qv(:,:,:,nnew), ie, je, ke )
      CALL clipping( qc(:,:,:,nnew), ie, je, ke )

#ifdef COSMOART
      ! CK 20101222 adapted COSMO-ART to changes in hordiff
      IF (l_cosmo_art) THEN
        IF (lgas) THEN
          DO isp=1,isp_gastrans
            kzdims(1:24)=(/ke,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
            CALL exchg_boundaries                                               &
              (0, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute,     &
              ie, je, kzdims, jstartpar, jendpar, 2, nboundlines, my_cart_neigh,&
              lperi_x, lperi_y, l2dim, 20000+ntstep, .FALSE., ncomm_type,       &
              izerror, yzerrmsg, cgas(:,:,:,isp,nnew))
          ENDDO
          DO isp=1,isp_gastrans
            CALL clipping( cgas(:,:,:,isp,nnew), ie, je, ke )
          ENDDO
        ENDIF

        IF (laero) THEN
          DO isp=1,isp_aerotrans
            kzdims(1:24)=(/ke,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
            CALL exchg_boundaries                                               &
              (0, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute,     &
              ie, je, kzdims, jstartpar, jendpar, 2, nboundlines, my_cart_neigh,&
              lperi_x, lperi_y, l2dim, 20000+ntstep, .FALSE., ncomm_type,       &
              izerror, yzerrmsg, caero(:,:,:,isp,nnew))
          ENDDO
          DO isp=1,isp_aerotrans
            CALL clipping( caero(:,:,:,isp,nnew), ie, je, ke )
          ENDDO
        ENDIF

      ENDIF
#endif

#ifdef POLLEN
      IF (l_pollen) THEN
        DO isp=1,isp_pollen
          kzdims(1:24)=(/ke,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
          CALL exchg_boundaries                                               &
            (0, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute,     &
            ie, je, kzdims, jstartpar, jendpar, 2, nboundlines, my_cart_neigh,&
            lperi_x, lperi_y, l2dim, 20000+ntstep, .FALSE., ncomm_type,       &
            izerror, yzerrmsg, cpollen(:,:,:,isp,nnew))
        ENDDO
        DO isp=1,isp_pollen
          CALL clipping( cpollen(:,:,:,isp,nnew), ie, je, ke )
        ENDDO
      ENDIF
#endif

!MU (13.04.2012)
      DO iig = 1, ntracer
        IF (ltracer(2,iig) == 1) THEN
            IF (num_compute > 1) THEN
              kzdims(1:24)=(/ke,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
              CALL exchg_boundaries                                              &
                (0, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute,    &
                                                                      ie, je,    &
                kzdims, jstartpar, jendpar, 2, nboundlines, my_cart_neigh,       &
                          lperi_x, lperi_y, l2dim,                               &
                20000+ntstep, .FALSE., ncomm_type, izerror, yzerrmsg,         &
                tracer (:,:,:,nnew,iig))
            ENDIF
            IF (ltracer(2,iig) == 1) THEN
              CALL clipping( tracer(:,:,:,nnew,iig), ie, je, ke )
            ENDIF
        ENDIF
      ENDDO
!MU (13.04.2012)

    ELSE
!!$      IF (num_compute > 1) THEN
        kzdims(1:24)=(/ke,ke,ke1,ke,ke,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
        CALL exchg_boundaries                                                 &
          ( 0, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute,      &
           ie, je, kzdims, jstartpar, jendpar, 2, nboundlines, my_cart_neigh, &
          lperi_x, lperi_y, l2dim,                                            &
          20000+ntstep, .FALSE., ncomm_type, izerror, yzerrmsg,               &
          u (:,:,:,nnew), v (:,:,:,nnew), w (:,:,:,nnew), t (:,:,:,nnew),     &
          pp(:,:,:,nnew))
!!$      ENDIF
!!$      IF ( l2dim ) THEN
!!$        DO j = 1, ny_2dim
!!$          u (:,j,:,nnew) = u (:,j2dim,:,nnew)
!!$          v (:,j,:,nnew) = v (:,j2dim,:,nnew)
!!$          w (:,j,:,nnew) = w (:,j2dim,:,nnew)
!!$          t (:,j,:,nnew) = t (:,j2dim,:,nnew)
!!$          pp(:,j,:,nnew) = pp(:,j2dim,:,nnew)
!!$        ENDDO
!!$      ENDIF
    ENDIF
    
  ENDIF

  IF (ltime) CALL get_timings (i_horizontal_diffusion, ntstep, dt, izerror)

  !----------------------------------------------------------------------------
  ! Section 6: First saturation adjustment for t, qv, qc at time level nnew
  !----------------------------------------------------------------------------

  IF ( ldiabf_lh ) THEN
    ! initialize temperature increment due to latent heat
    IF ( ldiabf_satad ) THEN
      tinc_lh(:,:,:) = - t(:,:,:,nnew)
    ELSE
      tinc_lh(:,:,:) = 0.0_ireals
    END IF
  END IF
  
  IF (lcond) THEN

#ifdef NECSX
  CALL collapse(.TRUE., ie, je, istartpar, iendpar, jstartpar, jendpar)
#endif

    ! The array wcon is allocated again here but just as intermediate storage
    ! for routine satad.
    ALLOCATE (wcon(ie,je,9), STAT=izstata )
    IF ( izstata /= 0 ) THEN
      yzerrmsg="allocation of wcon"
      CALL model_abort (my_cart_id, 100+izerror, yzerrmsg, 'org_runge_kutta')
    END IF

    ! Determine the maxima of w(:,:,k,:) for reproducible results in the
    ! parallel program
    DO k = 1,ke
      zwmax(k) = MAXVAL ( ABS(w(istartpar:iendpar,jstartpar:jendpar,k,nnew)) )
    ENDDO

    IF ( (lreproduce) .AND. (num_compute > 1) ) THEN
      IF (ltime) THEN
        CALL get_timings (i_dyn_computations, ntstep, dt, izerror)
        IF (ltime_barrier) THEN
          CALL comm_barrier (icomm_cart, izerror, yzerrmsg)
          CALL get_timings (i_barrier_waiting_dyn, ntstep, dt, izerror)
        ENDIF
      ENDIF

      CALL global_values (zwmax, ke, 'MAX', imp_reals, icomm_cart, -1,    &
        yzerrmsg, izerror)

      IF (ltime) CALL get_timings (i_communications_dyn, ntstep, dt, izerror)
    ENDIF

#ifdef NUDGING
    IF (llhn .OR. llhnverif) THEN
      ! For calculation of latent heating rate 
      ! T has to be stored before saturation adjustment
      IF (itype_gscp>=3) THEN
        ! in case of itype_gscp==3 (hydci or hydci_pp) get latent heat rate 
        ! of the last time step and store it on the next time step. Both 
        ! Routines are runnig at the end of the time step and therefore after 
        ! the latent heat nudging.
        tt_lheat(:,:,:,nnew)=tt_lheat(:,:,:,nnow)
      ENDIF

      CALL get_gs_lheating ('add', 1, ke)
    ENDIF
#endif

    DO  k = 1, ke
      DO j = jstartpar, jendpar
        DO i = istartpar, iendpar
          wcon(i,j,1) =  p0(i,j,k) + pp(i,j,k,nnew)
        ENDDO
      ENDDO
      kitpro       = 1
      IF (zwmax(k) >=  2.0) kitpro = 2
      IF (zwmax(k) >= 10.0) kitpro = 3

      CALL satad ( kitpro, t(:,:,k,nnew), qv(:,:,k,nnew),       &
        qc(:,:,k,nnew), t(:,:,k,nnow), wcon(:,:,1),             &
        wcon(:,:,2), wcon(:,:,3), wcon(:,:,4), wcon(:,:,5),     &
        wcon(:,:,6), wcon(:,:,7), wcon(:,:,8), wcon(:,:,9),     &
        b1, b2w, b3, b4w, b234w, rdv, o_m_rdv,                  &
        rvd_m_o, lh_v, cpdr, cp_d,                              &
        ie, je, istartpar, iendpar, jstartpar, jendpar  )
    ENDDO

    DEALLOCATE (wcon, STAT=izstatd)

#ifdef NUDGING
    ! calculate gridscale latent heating rate from the saturation adjustment
    IF (llhn .OR. llhnverif) CALL get_gs_lheating ('inc', 1, ke)
#endif

#ifdef NECSX
  CALL collapse(.FALSE., ie, je, istartpar, iendpar, jstartpar, jendpar)
#endif

  ENDIF

  IF ( ldiabf_satad ) THEN
    ! compute temperature increment due to latent heat
    tinc_lh(:,:,:) = tinc_lh(:,:,:) + t(:,:,:,nnew)
  END IF
  
  IF (ltime) CALL get_timings (i_dyn_computations, ntstep, dt, izerror)

  !----------------------------------------------------------------------------
  ! End of module procedure org_runge_kutta
  !----------------------------------------------------------------------------

END SUBROUTINE org_runge_kutta

!==============================================================================
!+ Module procedure in "src_runge_kutta" for computing small time step
!------------------------------------------------------------------------------

SUBROUTINE calc_small_timestep( dt, ismtstep, ismtstep_1, ismtstep_2, &
                                dtsmall, dtsmall_1, dtsmall_2 )

!----------------------------------------------------------------------------
!
! Description:
!   Calculation of the small time step sizes dtsmall, _1, _2 and the number of
!   small steps (ismtstep, _1, _2) for each Runge-Kutta-substep
!
! Method:
!
!----------------------------------------------------------------------------

  ! Subroutine arguments:
  ! ---------------------

  REAL    (KIND=ireals),    INTENT(IN)  :: dt
  INTEGER (KIND=iintegers), INTENT(OUT) :: ismtstep, ismtstep_1, ismtstep_2
  REAL    (KIND=ireals),    INTENT(OUT) :: dtsmall,  dtsmall_1,  dtsmall_2

  ! Local scalars:
  ! -------------

  REAL    (KIND=ireals) ::  &
    z_cs, z_dtsmax, z_dx, z_dy, z_crlat_min, z_dt

  INTEGER (KIND=iintegers) :: &
    z_ismtstep_max, z_ismtstep_sum, izerror

  CHARACTER (LEN=80) ::  &
    yzerrmsg

  INTEGER (KIND=iintegers) :: icalc_version  ! 1 = as before, 0 = Standard RK

  !----------------------------------------------------------------------------
  ! Begin Subroutine calc_small_timestep
  !----------------------------------------------------------------------------

  icalc_version = 0
  !WRITE(*,*) "Subr [calc_small_timestep]: icalc_version =", icalc_version

  z_crlat_min = MINVAL( crlat(:,1) )

  IF (num_compute > 1) THEN
    IF (ltime) THEN
      CALL get_timings (i_dyn_computations, ntstep, dt, izerror)
      IF (ltime_barrier) THEN
        CALL comm_barrier (icomm_cart, izerror, yzerrmsg)
        CALL get_timings (i_barrier_waiting_dyn, ntstep, dt, izerror)
      ENDIF
    ENDIF

    CALL global_values (z_crlat_min, 1, 'MIN', imp_reals, icomm_cart, -1, &
        yzerrmsg, izerror)

    IF (ltime) CALL get_timings (i_communications_dyn, ntstep, dt, izerror)
  ENDIF

  z_dt = ABS(dt)   ! in DFI backward stepping, the time step is negative
  z_dx = r_earth * 3.1415927_ireals / 180.0_ireals * dlon * z_crlat_min
  z_dy = r_earth * 3.1415927_ireals / 180.0_ireals * dlat

  z_cs        = SQRT( gamma * r_d * 303.0_ireals )
  z_dtsmax    = 0.9_ireals / z_cs / SQRT( 1.0/z_dx**2 + 1.0/z_dy**2 )
  z_dtsmax    = MIN( z_dtsmax, 32.0_ireals )      !MB: reason?

  SELECT CASE(irk_order)
  CASE(1)
    ! = Euler-Forward-scheme
    ismtstep = INT( z_dt/z_dtsmax, iintegers) + 1
    IF ( icalc_version == 1 ) THEN
      IF ( ismtstep == 1 ) ismtstep = 2
      dtsmall  = z_dt/ismtstep
    ELSE
      dtsmall  = z_dt/ismtstep
    ENDIF
    z_ismtstep_max = ismtstep
    z_ismtstep_sum = ismtstep
  CASE(2)
    ismtstep_1 = INT( ( z_dt/2.0_ireals) / z_dtsmax, iintegers) + 1
    IF ( ismtstep_1 == 1 ) ismtstep_1 = 2
    ismtstep   = INT( z_dt/z_dtsmax, iintegers) + 1
    IF ( ismtstep == 1 ) ismtstep = 2
    dtsmall_1  = (z_dt/2.0_ireals) / ismtstep_1
    dtsmall    = z_dt / ismtstep
    z_ismtstep_max = ismtstep
    z_ismtstep_sum = ismtstep_1 + ismtstep
  CASE(3)
    SELECT CASE(irunge_kutta)
    CASE(1)
      IF ( icalc_version == 1 ) THEN
        ismtstep_1 = INT( (z_dt/3.0_ireals) / z_dtsmax, iintegers) + 1
        IF ( ismtstep_1 == 1 ) ismtstep_1 = 2
        ismtstep_2 = INT( (z_dt/2.0_ireals) / z_dtsmax, iintegers) + 1
        IF ( ismtstep_2 == 1 ) ismtstep_2 = 2
        ismtstep   = INT( z_dt/z_dtsmax, iintegers) + 1
        IF ( ismtstep == 1 ) ismtstep = 2
        dtsmall_1  = (z_dt/3.0_ireals) / ismtstep_1
        dtsmall_2  = (z_dt/2.0_ireals) / ismtstep_2
        dtsmall    = z_dt / ismtstep
        z_ismtstep_max = ismtstep
      ELSE
        ! generate a ratio dt/dtsmall so, that numbers 2/3/6, or 4/6/12, ...
        ! of small steps in the 3 RK substeps are used
        ! (remark: 2/3/6 is more efficient than higher numbers,
        ! therefore it could even be advisable to reduce the BIG timestep dt
        ! a bit, to achieve these lower numbers)

        ismtstep = INT( z_dt / z_dtsmax, iintegers)
        ! round up to multiples of 6:
        ismtstep = 6 * INT( ismtstep / 6 + 1, iintegers )

        dtsmall = z_dt / ismtstep

        ismtstep_1 = ismtstep / 3
        dtsmall_1 = dtsmall

        ismtstep_2 = ismtstep / 2
        dtsmall_2 = dtsmall

        z_ismtstep_max = ismtstep
      ENDIF
    CASE(2)
      ismtstep_1 = INT( z_dt / z_dtsmax, iintegers) + 1
      IF ( ismtstep_1 == 1 ) ismtstep_1 = 2
      ismtstep_2 = INT( ( z_dt/4.0_ireals) / z_dtsmax, iintegers) + 1
      IF ( ismtstep_2 == 1 ) ismtstep_2 = 2
      ismtstep   = INT( (2.0_ireals * z_dt/3.0_ireals) / z_dtsmax, iintegers) + 1
      IF ( ismtstep == 1 ) ismtstep = 2
      dtsmall_1  = z_dt / ismtstep_1
      dtsmall_2  = (z_dt/4.0_ireals) / ismtstep_2
      dtsmall    = (2.0_ireals*z_dt/3.0_ireals) / ismtstep
      z_ismtstep_max = ismtstep_1
    END SELECT
    z_ismtstep_sum = ismtstep_1 + ismtstep_2 + ismtstep
  END SELECT

  ! Adapt the sign for the small steps (for DFI)
  IF (dt < 0.0_ireals) THEN
    dtsmall   = - dtsmall
    dtsmall_1 = - dtsmall_1
    dtsmall_2 = - dtsmall_2
  ENDIF

  IF (my_cart_id == 0) THEN
    PRINT *,"CS=",z_cs," DT=",dt," DTSMAX=" ,z_dtsmax
    IF      ( irk_order == 1 ) THEN
      PRINT *,"SHORT TIME STEP ", dtsmall, ", number of steps: ", ismtstep
    ELSE IF ( irk_order == 2 ) THEN
      PRINT *,"SHORT TIME STEP of 1. RK step ", dtsmall_1, ", number of steps: ", ismtstep_1
      PRINT *,"SHORT TIME STEP of 2. RK step ", dtsmall,   ", number of steps: ", ismtstep
    ELSE IF ( irk_order == 3 ) THEN
      PRINT *,"SHORT TIME STEP of 1. RK step ", dtsmall_1, ", number of steps: ", ismtstep_1
      PRINT *,"SHORT TIME STEP of 2. RK step ", dtsmall_2, ", number of steps: ", ismtstep_2
      PRINT *,"SHORT TIME STEP of 3. RK step ", dtsmall,   ", number of steps: ", ismtstep
    END IF
    PRINT *,"ismtstep_sum = ",z_ismtstep_sum," ismtstep_max = ", z_ismtstep_max
  ENDIF

END SUBROUTINE calc_small_timestep

!==============================================================================
!==============================================================================

SUBROUTINE check_cfl_horiz_advection ( lapply_Rayleigh_damping, epsray )

!------------------------------------------------------------------------------
!
! Description:
! 1.) Check, if the CFL-criterion for 2-dim. horizontal advection at
!     timelevel nnow is violated.
! 2.) If this is the case then set the Rayleigh damping coefficients
!
!------------------------------------------------------------------------------

!============================================================================

! Declarations:

! Subroutine arguments:
! ---------------------

REAL(KIND=ireals), INTENT(OUT) :: epsray(1:ie, 1:je, 1:ke)
      ! Rayleigh damping coefficient, if CFL-crit. is violated
LOGICAL, INTENT(OUT) ::  lapply_Rayleigh_damping

! Local scalars:
! -------------

INTEGER           :: i, j, k
REAL(KIND=ireals) :: cfl_limit_1D     ! maximum CFL in the 1-dim. case
REAL(KIND=ireals) :: cfl, cfl_max
REAL(KIND=ireals) :: dt_dx, dt_dy
REAL(KIND=ireals) :: epsray0
INTEGER           :: izerror, izdebug
CHARACTER(LEN=80) :: yzerrmsg
LOGICAL           :: lz_crit_point_already_printed

!----------------------------------------------------------------------------
! Begin Subroutine check_cfl_advection
!----------------------------------------------------------------------------

! Initialize, whether debug output shall be done
  IF (ldebug_dyn) THEN
    IF (lprintdeb_all) THEN
      izdebug = idbg_level
    ELSE
      IF (my_cart_id == 0) THEN
        izdebug = idbg_level
      ELSE
        izdebug = 0
      ENDIF
    ENDIF
  ELSE
    izdebug = 0
  ENDIF

  dt_dx = ABS(dt) / ( Pi * r_earth * dlon / 180.0_ireals ) 
          ! correction with 'crlat' see below
  dt_dy = ABS(dt) / ( Pi * r_earth * dlat / 180.0_ireals )

  ! 1D-CFL-criteria for different advection schemes
  ! (e.g. Wicker, Skamarock (2002)
  SELECT CASE ( irk_order )
  CASE ( 3 )

    SELECT CASE ( iadv_order )
    CASE(1)
      cfl_limit_1D = 1.24_ireals
    CASE(2)
      cfl_limit_1D = 1.72_ireals
    CASE(3)
      cfl_limit_1D = 1.61_ireals
    CASE(4)
      cfl_limit_1D = 1.25_ireals
    CASE(5)
      cfl_limit_1D = 1.42_ireals
    CASE(6)
      cfl_limit_1D = 1.08_ireals
    CASE DEFAULT
      cfl_limit_1D = 1.0_ireals
    END SELECT

  CASE default
    cfl_limit_1D = 1.0_ireals
  END SELECT

  ! calculate CFL-criterion
  cfl_max = 0.0_ireals
  DO k = 1,ke
    DO j = jstartv,jendv
      DO i = istartu,iendu
        cfl =   ABS( u(i,j,k,nnow) ) * dt_dx / crlat(j,1)  &
              + ABS( v(i,j,k,nnow) ) * dt_dy
        cfl_max = MAX( cfl_max, cfl )
      ENDDO
    ENDDO
  ENDDO

  IF (num_compute > 1) THEN
    IF (ltime) THEN
      CALL get_timings (i_dyn_computations, ntstep, dt, izerror)
      IF (ltime_barrier) THEN
        CALL comm_barrier (icomm_cart, izerror, yzerrmsg)
        CALL get_timings (i_barrier_waiting_dyn, ntstep, dt, izerror)
      ENDIF
    ENDIF

    CALL global_values ( cfl_max, 1, 'MAX', imp_reals, icomm_cart, -1,  &
      yzerrmsg, izerror)

    IF (ltime) CALL get_timings (i_communications_dyn, ntstep, dt, izerror)
  ENDIF

  IF (izdebug > 2) THEN
    PRINT *, 'CFL=', cfl_max, ', Max (|u|+|v|) ~ ', cfl_max / dt_dx
    ! this is only a rough estimation for the velocity,
    ! assumption: dt_dx = dt_dy = const. is used
  END IF

  IF ( cfl_max > 0.95_ireals * cfl_limit_1D ) THEN

    IF ( my_cart_id == 0 ) THEN
      PRINT*,' !!!!*** WARNING ***!!! CFL-criterion for horizontal &
               &advection is violated'
      PRINT*,' !!!! Max ( |cfl_x| + |cfl_y| ) = ',cfl_max,  &
                            " ( > 0.95 * ", cfl_limit_1D, " )"
    ENDIF

    ! calculate the Rayleigh-damping-coefficient
    ! and search again the point of the CFL-violation for debugging

    lapply_Rayleigh_damping = .TRUE.

    epsray0 = 1.0_ireals / (1000.0_ireals * ABS(dt))
    epsray0 = epsray0 * (cfl_max - 0.95_ireals * cfl_limit_1D) /  &
                        (0.05_ireals * cfl_limit_1D)
    ! stability constraint
    epsray0 = SIGN(1._ireals,dt) * MIN( epsray0, 1.0_ireals / (2.0 * ABS(dt) ) )

    lz_crit_point_already_printed = .FALSE.

    epsray(:,:,:)=0.0_ireals

    DO k = 1,ke
      DO j = jstartv,jendv
        DO i = istartu,iendu

          cfl =   ABS( u(i,j,k,nnow) ) * dt_dx / crlat(j,1)  &
                + ABS( v(i,j,k,nnow) ) * dt_dy

          ! reduce only high velocities:
          IF ( cfl > 0.90_ireals * cfl_limit_1D ) THEN
            epsray(i,j,k) = epsray0
          ELSE
            epsray(i,j,k) = 0.0_ireals
          END IF

          IF ( .NOT.(lz_crit_point_already_printed) .AND.   &
               (cfl >= 0.9999_ireals * cfl_max )   ) THEN
            PRINT *, ' global indices: i=', i_global(i),                  &
                                    ', j=', j_global(j), ', k=', k
            PRINT *, "|cfl_x| + |cfl_y| = ", cfl, ", u=", u(i,j,k,nnow),  &
                                                  ", v=", v(i,j,k,nnow)
            PRINT *, "epsray0 = ", epsray0
            lz_crit_point_already_printed = .TRUE.
               ! printing of only one grid point with
               ! extremely high velocities is sufficient
          END IF

        ENDDO
      ENDDO
    ENDDO

  ELSE

    lapply_Rayleigh_damping = .FALSE.
    !epsray(:,:,:) = 0.0_ireals

  ENDIF

END SUBROUTINE check_cfl_horiz_advection



!==============================================================================


SUBROUTINE coriolis( utens, vtens, wtens, nn ) 

!-----------------------------------------------------
!
! Description:
! Add Coriolis force to u- and v-tendencies
!
! Input:
!   u(:,:,:,nn), v(:,:,:,nn), w(:,:,:,nn)  at timelevel nn
!   fc, fccos
!
! Input- and Output:
!   utens, vtens, wtens
!
!------------------------------------------------------
  REAL(KIND=ireals), INTENT(INOUT) :: utens(ie,je,ke), vtens(ie,je,ke)
  REAL(KIND=ireals), INTENT(INOUT) :: wtens(ie,je,ke1)
  INTEGER (KIND=iintegers), INTENT(IN) :: nn

  INTEGER          :: i,j,k
  REAL(KIND=ireals):: z_fv_north, z_fv_south, zfq
  REAL(KIND=ireals):: z_fcw_west, z_fcw_east, zfcw 
  REAL(KIND=ireals):: z_fu_west,  z_fu_east
  REAL(KIND=ireals):: zfcu

  DO k = 1 , ke
    DO j = jstartu, jendu
!CDIR ON_ADB(fc)
!CDIR ON_ADB(v)
      DO i = istartu, iendu
        z_fv_north = fc(i,j)   * ( v(i,j  ,k,nn) + v(i+1,j  ,k,nn) )
        z_fv_south = fc(i,j-1) * ( v(i,j-1,k,nn) + v(i+1,j-1,k,nn) )
        zfq = 0.25_ireals * ( z_fv_north + z_fv_south )
        utens(i,j,k) = utens(i,j,k) + zfq
      ENDDO
    ENDDO

    DO j = jstartv, jendv
!CDIR ON_ADB(fc)
!CDIR ON_ADB(u)
      DO i = istartv, iendv
        z_fu_east = fc(i,j)   * ( u(i  ,j,k,nn) + u(i  ,j+1,k,nn) )
        z_fu_west = fc(i-1,j) * ( u(i-1,j,k,nn) + u(i-1,j+1,k,nn) )
        zfq = 0.25_ireals * ( z_fu_east + z_fu_west )
        vtens(i,j,k) = vtens(i,j,k) - zfq 
      ENDDO
    ENDDO
  END DO

  IF (lcori_deep) THEN
    ! additional 2*Omega*cos(phi) * w -term in the momentum equation for u
    DO k = 1 , ke
      DO j = jstartu, jendu
        DO i = istartu, iendu
          z_fcw_west = fccos(i  ,j) * ( w(i  ,j,k,nn) + w(i  ,j,k+1,nn) )
          z_fcw_east = fccos(i+1,j) * ( w(i+1,j,k,nn) + w(i+1,j,k+1,nn) )
          zfcw = 0.25_ireals * ( z_fcw_west + z_fcw_east )
          utens(i,j,k) = utens(i,j,k) - zfcw
        ENDDO
      ENDDO
    END DO

    ! additional 2*Omega*cos(phi) * u -term in the momentum equation for w
    DO k = 2 , ke
      ! no calculations for wtens(:,:,1/ke) due to boundary conditions
      DO j = jstart, jend
        DO i = istart, iend
          zfcu = 0.25_ireals * fccos(i,j)                   &
               * ( u(i-1,j,k-1,nn) + u(i,j,k-1,nn) +        &
                   u(i-1,j,k  ,nn) + u(i,j,k  ,nn) )
          wtens(i,j,k) = wtens(i,j,k) + zfcu
        ENDDO
      ENDDO
    ENDDO
  ENDIF

END SUBROUTINE coriolis


SUBROUTINE calc_wcon_sqrtg( u, v, w, wcon_sqrtg )

!-------------------------------------------------------------
! Description:
!   calculate the product
!   contravariant vertical velocity wcon * sqrt(G)
!   from the spherical wind components u, v, w 
! 
! Input:
!   u(:,:,:), v(:,:,:), w(:,:,:)
!
! Output: 
!   wcon_sqrtg(:,:,:)
!-------------------------------------------------------------

  USE src_advection_rk, ONLY: horiz_adv_driver

  USE data_constants, ONLY:   &
    r_earth           ! mean radius of the earth
    
  USE data_fields, ONLY:   &
    hhl        ,    & ! geometical height of half model levels        (  m  )
    acrlat            ! 1 / ( crlat * radius of the earth )           ( 1/m )
    ! sqrtg_r_w         ! reciprocal square root of G at w points       ( 1/m )

  USE fast_waves_rk, ONLY:   wgtfac_u, wgtfac_v   ! weighting factors for linear interpolation from full to half levels

  IMPLICIT NONE 

  REAL  (KIND=ireals), INTENT(IN):: &
    &    u(1:ie, 1:je, 1:ke), &
    &    v(1:ie, 1:je, 1:ke), &
    &    w(1:ie, 1:je, 1:ke+1)

  REAL  (KIND=ireals), INTENT(OUT):: &
    &    wcon_sqrtg(1:ie, 1:je, 1:ke+1)


  REAL  (KIND=ireals), ALLOCATABLE ::  &
    zu   (:,:,:),    & !
    zv   (:,:,:)

  REAL  (KIND=ireals)       :: zsign

  INTEGER  (KIND=iintegers) :: i, j, k

  INTEGER  (KIND=iintegers) :: izstata, izstatd, izerror
  CHARACTER(LEN=80)         :: yzerrmsg

  REAL    (KIND=ireals   )  :: r_earth_recip


  ALLOCATE( zu(1:ie, 1:je, 1:ke1),  &
    &       zv(1:ie, 1:je, 1:ke1),  &
    &       STAT=izstata)
  IF ( izstata /= 0 ) THEN
    yzerrmsg="allocation of zu, zv"
    CALL model_abort (my_cart_id, 100+izerror, yzerrmsg, 'calc_wcon_sqrtg')
  END IF

  zsign = SIGN(1._ireals,dt)

  r_earth_recip = 1.0_ireals / r_earth

  ! calc.  u= d lambda/dt  and  v= d phi/dt and average to w-position:
  !DO  k = 2, ke
  !  DO j = jstart, jend
  !    DO i = istart, iend
  !      zu(i,j,k) = 0.25_ireals * ( u(i-1,j  ,k-1) + u(i,j,k-1)    &
  !        &                       + u(i-1,j  ,k  ) + u(i,j,k  ) )  &
  !        &                 * acrlat(j,1)
  !      zv(i,j,k) = 0.25_ireals * ( v(i  ,j-1,k-1) + v(i,j,k-1)    &
  !        &                       + v(i  ,j-1,k  ) + v(i,j,k  ) )  &
  !        &                 * r_earth_recip
  !    END DO
  !  END DO
  !END DO

  ! calc.  u= d lambda/dt  and  v= d phi/dt and average to w-position:
  DO  k = 2, ke
    DO j = jstart, jend
      DO i = istart, iend
        zu(i,j,k) = (                wgtfac_u(i-1,j,k)   * u(i-1,j,k  )    &
          &         + ( 1.0_ireals - wgtfac_u(i-1,j,k) ) * u(i-1,j,k-1)    &
          &         +                wgtfac_u(i  ,j,k)   * u(i,  j,k  )    &
          &         + ( 1.0_ireals - wgtfac_u(i  ,j,k) ) * u(i,  j,k-1) )  &
          &                 * 0.5_ireals * acrlat(j,1)
        zv(i,j,k) = (                wgtfac_v(i,j-1,k)   * v(i,j-1,k  )    &
          &         + ( 1.0_ireals - wgtfac_v(i,j-1,k) ) * v(i,j-1,k-1)    &
          &         +                wgtfac_v(i,j  ,k)   * v(i,j  ,k  )    &
          &         + ( 1.0_ireals - wgtfac_v(i,j  ,k) ) * v(i,j  ,k-1) )  &
          &                 * 0.5_ireals * r_earth_recip
      END DO
    END DO
  END DO

  ! Calculation of the contravariant vertical velocity (wcon):

  wcon_sqrtg(:,:,:) = 0.0_ireals

  ! (negative) tendency of horizontal advection of z:

  CALL horiz_adv_driver( hhl, zu, zv, wcon_sqrtg, zsign,  &
    &             istart, iend, jstart, jend, 2, ke, ke1, iadv_order )

  ! die folgenden beiden Zeilen sollten NICHT noetig sein:
  wcon_sqrtg(:,:,1)    = 0.0_ireals
  wcon_sqrtg(:,:,ke+1) = 0.0_ireals

  DO k = 2, ke
    DO j = jstart, jend
      DO  i = istart, iend
        !wcon(i,j,k) = sqrtg_r_w(i,j,k) * ( -wcon(i,j,k) - w(i,j,k) )
        wcon_sqrtg(i,j,k) = -wcon_sqrtg(i,j,k) - w(i,j,k)
      ENDDO
    ENDDO
  ENDDO

  DEALLOCATE( zu, zv, STAT=izstatd )
  IF ( izstatd /= 0 ) THEN
    yzerrmsg = 'deallocation of zu, zv'
    CALL model_abort (my_cart_id, 100+izerror, yzerrmsg, 'calc_wcon_sqrtg')
  END IF

END SUBROUTINE  calc_wcon_sqrtg



END MODULE src_runge_kutta
