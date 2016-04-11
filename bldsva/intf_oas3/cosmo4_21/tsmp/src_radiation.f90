!+ Source module  "src_radiation"
!------------------------------------------------------------------------------

MODULE src_radiation

!------------------------------------------------------------------------------
!
! Description:
!   The module "radiation" performs calculations related to the 
!   parameterization of radiative transfer.
!   Driving routine is the model procedure "organize_radiation", which
!   calls the main routine "fesft" of the radiation package.
!   Additionally, some diagnostics and gridpoint output is done.
!
!   All global variables of the model that are used by the interface routine
!   "organize_radiation" are imported by USE statements below.
!
!   The program package for the parameterization of radiative transfer consists
!   of following routines:
!         fesft, coe_so, coe_th, inv_so, ivs_th, opt_so and opt_th.
!
!   All parametric data that are required by these routines are defined in the
!   data module data_radiation.
!
!   Additionally, the routine init_radiation has to be called once before the 
!   first call of the driving routine fesft for the radiation package in the 
!   driving routine organize_radiation. Aerdis is a small help routine to 
!   receive some parameters for the vertical distribution of background 
!   aerosol in the driving routine.
!
!   The parameterization package has been provided by B. Ritter in a
!   Plug-compatible Fortran77-Version. Some technical modifications have been
!   done for the F90-Version:
!   Internal communication by common-blocks is replaced by module parameters,
!   scalars and arrays defined in module data_radiation. 
!
! Current Code Owner: DWD, Bodo Ritter
!  phone:  +49  69  8062 2703
!  fax:    +49  69  8062 3721
!  email:  Bodo.Ritter@dwd.de
!
! History:
! Version    Date       Name
! ---------- ---------- ----
! 1.1        1998/03/11 Guenter Doms     
!  Initial release
! 1.4        1998/05/22 Guenter Doms
!  Inclusion of control parameter l2tls to select time levels
!  according to the time integration scheme used.
! 1.8        1998/08/03 Ulrich Schaettler
!  Elmination of dependency from module data_io
! 1.20       1999/01/07 Guenther Doms
!  Renaming of some global variables
! 1.24       1999/03/01 Guenther Doms
!  Inclusion of the new prognostic 3-D array 'qi'.
! 1.29       1999/05/11 Ulrich Schaettler
!  Adapted interfaces to utility-modules
! 1.30       1999/06/24 Matthias Raschendofer
!  Use from module data_constants: rvd_m_o, lhop, b234w.
!  Use form module data_fields: rcld.
!  Use form module data-runcontrol: itype_wcld, icldm_rad.
!  Use form module meteo_utilities: cloud_diag
! 1.32       1999/08/24 Guenther Doms
!  top_con and bas_con removed from the use-list
! 1.34       1999/12/10 Ulrich Schaettler
!  Use new timing routines
! 1.39       2000/05/03 Ulrich Schaettler
!  Use variables klv800 and klv400 from data_modelconfig.
! 2.18       2002/07/16 Reinhold Schrodin
!  Eliminated variable rhde, use cf_snow instead
!  Use new variable for soil moisture in case of multi-layer soil model
! 3.5        2003/09/02 Ulrich Schaettler
!  Avoid global communication by providing the fields phi_tot, rla_tot
! 3.6        2003/12/11 Reinhold Schrodin
!  Adaptations for multi-layer soil model
! 3.7        2004/02/18 Ulrich Schaettler
!  Replace local logical variable lcl_ice by global variable lprog_qi
!  Renamed alb (alb_rad), phi (rlat), rla (rlon)
! 3.16       2005/07/22 Reinhold Schrodin
!  Use variables for_e and for_d from data_fields, csalb_snow_fd and
!  csalb_snow_fe from data_soil
! 3.18       2006/03/03 Ulrich Schaettler
!  Really included all the *.incf files in src_radiation.f90
!  Add variables and fields used by the lake model FLake
!  Some additionals for the Climate LM Version
! 3.21       2006/12/04 Burkhardt Rockel, Ulrich Schaettler
!  Added A2 scenarios for CO2; renamed variable sodwdir to sodwddm
!  Use alternatively t_g for determining soil type of grid cell
!  Use new NL variables from data_constants: clc_diag, q_crit
! 3.22       2007/01/24 Axel Seifert
!  Adapted a constant for computation of cloud cover of ice clouds
! V3_23        2007/03/30 Matthias Raschendorfer, Matteo Buzzi
!  Moving 'clc_diag' and 'q_crit' to MODULE data_turbulence.
!  Calculation of topographical correction if lradtopo=.TRUE.
! V4_1         2007/12/04 Ulrich Schaettler
!  Bug correction (found out by Catherine Meissner, IMK Karlsruhe):
!  The downward and upward component have to be changed in call to SR fesft
! V4_3         2008/02/25 Ulrich Schaettler
!  There were also other downward and upward component in the wrong order
! V4_4         2008/07/16 Ulrich Schaettler
!  Eliminated timing variables which are unused
!  Changed NL parameter lyear_360 to itype_calendar, to have several options
! V4_8         2009/02/16 Ulrich Schaettler; Guenther Zaengl
!  Included additional (output) fields (for CLM: Uli)
!  Use p0hl (reference pressure at half levels) for full consistency with
!  new reference atmosphere implementation
! V4_9         2009/07/16 Ulrich Schaettler
!  Removed boundary exchange for radiation averaging (now in organize_physics)
!  Introduced full 3D fields of local variables in organize_radiation for a
!   better vectorization of the interpolation for radiation averaging
!  Implemented COSMO-ART features (with ifdef)
! V4_10        2009/09/11 Matthias Raschendorfer, Jan-Peter Schulz
!  Bug correction for the case "icldm_rad==2" when coud ice is present.
!  Modifications for the new seaice scheme (Jan-Peter Schulz)
!  Additional vector optimizations by NEC
! V4_11        2009/11/30 Ekaterina Machulskaya, Juergen Helmert
!  Adaptations for running with multi-layer snow model
!  Implemented options for aerosol distribution and use of an external 
!      emissivity map (JH)
! V4_12        2010/05/11 Ulrich Schaettler
!  Renamed t0 to t0_melt because of conflicting names
!  Bug fix in call to exchg_boundaries (wrong vertical dimension) (Lucio Torrisi)
!  Init and use unified variables for all aerosol-types (JH)
!  Compute additional fields for sunshine duration (Oli Fuhrer)
!  Adaptations in the COSMO-ART part
! V4_13        2010/05/11 Michael Gertz
!  Adaptions to SVN
! V4_17        2011/02/24 Ulrich Blahak
!  Adapted interface of exchg_boundaries; corrected kzdims(1:20) -> kzdims(1:24);
!  eliminated my_peri_neigh
! V4_18        2011/05/26 Ulrich Schaettler
!  Use local 2D array zskyview and use 1D slice of it (because global array is
!   allocated only for lradtopo, but was used in SR interfaces)
!  Smoothing field swdifu_s (for output) after coarse grid computation
!  Use of lgas and call to SR calcjval of COSMO-ART only for lgas (Christoph Knote)
!  Bug fix for computing qc_rad, qi_rad in case of a coarser grid (Victor Venema)
! V4_20        2011/08/31 Juergen Helmert / Victor Venema
!  Bug correction in the sub grid scale cloudiness
!  Replaced arguments p0 and pp in SR cloud_diag by p_tot (Uli Schaettler)
!  Bug fix in setting full variables after coarse radiation step (TR)
!
! CPS: Updated by Markus Übel(C4) based on the TerrSysMPV1.0 formulation in cosmo4.11
! Sep 05 2013: Fixed sun angle and aerosol for periodic boundary condition (as discussion
!              with Uli Blahak) : P. Shrestha
! Jan 21 2014: Update for masked simulation on calculation of parallel albedo
! for sea points : P.Shrestha (see zpalp estimation for ocean)
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

USE data_parameters, ONLY :   &
    ireals,    & ! KIND-type parameter for real variables
    iintegers    ! KIND-type parameter for standard integer variables

!------------------------------------------------------------------------------

USE data_modelconfig, ONLY :   &

! 1. vertical coordinate parameters and related variables
! -------------------------------------------------------
    sigmr,        & ! sigma-coordinate refering to PMSL             half level

! 2. horizontal and vertical sizes of the fields and related variables
! --------------------------------------------------------------------
    ie,           & ! number of grid points in zonal direction
    ie_tot,       & ! the same for the total field
    je,           & ! number of grid points in meridional direction
    ke,           & ! number of grid points in vertical direction
    ke1,          & ! ke+1
    czmls,        & ! depth of the soil layers in meters

! 3. start- and end-indices for the computations in the horizontal layers
! -----------------------------------------------------------------------
!    These variables give the start- and the end-indices of the
!    forecast for the prognostic variables in a horizontal layer.
!    Note, that the indices for the wind-speeds u and v differ from
!    the other ones because of the use of the staggered Arakawa-C-grid.

!   zonal direction
    istart,       & ! start index for the forecast of w, t, qd, qw and pp
    iend,         & ! end index for the forecast of w, t, qd, qw and pp
    istartpar,    & ! start index for computations in the parallel program
    iendpar,      & ! end index for computations in the parallel program

!   meridional direction
    jstart,       & ! start index for the forecast of w, t, qd, qw and pp
    jend,         & ! end index for the forecast of w, t, qd, qw and pp
    jstartpar,    & ! start index for computations in the parallel program
    jendpar,      & ! end index for computations in the parallel program

! 4. variables for the time discretization and related variables
! --------------------------------------------------------------
    dt,           & ! long time-step
!CPS00
    degrad,       & ! factor for transforming degree to rad

! 5. constants for the horizontal rotated grid and related variables
! ------------------------------------------------------------------

    pollon,       & ! longitude of the rotated north pole (in degrees, E>0)
    pollat,       & ! latitude of the rotated north pole (in degrees, N>0)

!CPS00

! 7. Layer index corresponding to a specified pressure
! ----------------------------------------------------

    klv800,       & ! k index of the LM-mainlevel, on 800 HPa
    klv500,       & ! k index of the LM-mainlevel, on 500 HPa
    klv400          ! k index of the LM-mainlevel, on 400 HPa

! end of data_modelconfig

!------------------------------------------------------------------------------

USE data_constants  , ONLY :   &

! 1. mathematical constants
! -------------------------
    pi,           & ! circle constant

! 2. physical constants and related variables
! -------------------------------------------

    t0_melt,      & ! melting temperature of ice
    r_v,          & ! gas constant for water vapor
    rdv,          & ! r_d / r_v
    o_m_rdv,      & ! 1 - r_d/r_v
    rvd_m_o,      & ! r_v/r_d - 1
    cp_d,         & ! specific heat of dry air at constant pressure
    cpdr,         & ! 1 / cp_d
    rdocp,        & ! r_d / cp_d
    lh_v,         & ! latent heat of vapourization
    lhocp,        & ! lh_v/cp_d
    g,            & ! acceleration due to gravity
    sigma,        & ! Boltzmann-constant
    solc,         & ! solar constant

! 3. constants for parametrizations
! ---------------------------------
    b1,           & ! variables for computing the saturation vapour pressure
    b2w,          & ! over water (w) and ice (i)
    b2i,          & !               -- " --
    b3,           & !               -- " --
    b4w,          & !               -- " --
    b4i,          & !               -- " --
    b234w,        & ! b2w * (b3 - b4w)
    uc1,          & ! variable for computing the rate of cloud cover in 
    uc2,          & ! the unsaturated case
    ucl             !               -- " --

! end of data_constants
!------------------------------------------------------------------------------

USE data_turbulence     , ONLY :   &

! 1. tuning constants for statistical cloud scheme:
! ------------------------------------------------------------

    clc_diag,     & ! cloud cover at saturation in statistical cloud diagnostic
    q_crit          ! critical value for normalized over-saturation

!------------------------------------------------------------------------------

USE data_fields     , ONLY :   &

! 1. constant fields for the reference atmosphere                     (unit)
! -----------------------------------------------
    p0         ,    & ! reference pressure at full levels             ( pa  )
    p0hl       ,    & ! reference pressure at half levels             ( Pa  )
    dp0        ,    & ! pressure thickness
    depth_lk   ,    & ! lake depth                                    (  m  )

! 2. external parameter fields                                        (unit)
! ----------------------------
    soiltyp    ,    & ! type of the soil (keys 0-9)                     --
    vio3       ,    & ! vertical integrated ozone contents            (pa O3)
    hmo3       ,    & ! ozone maximum                                 ( pa  )
    rlat       ,    & ! geographical latitude                         ( rad )
    rlon       ,    & ! geographical longitude                        ( rad )
    rlattot    ,    & ! geographical latitude                         ( rad )
    rlontot    ,    & ! geographical longitude                        ( rad )
    aer_su    ,    & ! monthly aerosol climatology sulfate drops               (0 - 1)
    aer_du    ,    & ! monthly aerosol climatology total dust                  (0 - 1)
    aer_or    ,    & ! monthly aerosol climatology organic (water sol.)        (0 - 1)
    aer_bc    ,    & ! monthly aerosol climatology black carbon                (0 - 1)
    aer_ss    ,    & ! monthly aerosol climatology sea salt                    (0 - 1)
    emis_rad  ,    & ! thermal surface emissivity                              (0 - 1)
    aerlan     ,    & ! aerosol-distribution for rural areas            --
    aerurb     ,    & ! aerosol-distribution for urban areas            --
    aerdes     ,    & ! aerosol-distribution for desert areas           --
    aersea     ,    & ! aerosol-distribution for sea                    --
    plcov      ,    & ! fraction of plant cover                         --
    llandmask  ,    & ! landpoint mask
    for_e      ,    & ! ground fraction covered by evergreen forest     --
    for_d      ,    & ! ground fraction covered by deciduous forest     --

! 3. prognostic variables                                             (unit)
! -----------------------
    t          ,    & ! temperature                                   (  k  )
    qv         ,    & ! specific water vapor content                  (kg/kg)
    qc         ,    & ! specific cloud water content                  (kg/kg)
    qi         ,    & ! specific cloud ice content                    (kg/kg)
    pp         ,    & ! deviation from the reference pressure         ( pa  )

! 5. fields for surface values and soil model variables               (unit )
! -----------------------------------------------------
    ps        ,     & ! surface pressure                              ( pa  )
    t_s       ,     & ! surface temperature                           (  K  )
    t_snow    ,     & ! temperature of the snow-surface               (  k  )
    t_snow_mult,    & ! temperature of the snow-surface               (  k  )
    t_g       ,     & ! weighted surface temperature                  (  k  )
    w_snow    ,     & ! water content of snow                         (m H2O)
    w_g1      ,     & ! water content of the upper soil layer         (m H2O)
    w_so              ! multi-layer soil moisture                     (m H2O)

USE data_fields     , ONLY :   &

!   fields for prognostic variables of the lake model FLake or ocean
!   variables
    t_ice      ,    & ! temperature of ice/water surface              (  K  )
    h_ice      ,    & ! lake/sea ice thickness                        (  m  )

! 6. fields that are computed in the parametrization and dynamics     (unit )
! ---------------------------------------------------------------
!   turbulence statistics in the atmosphere
!   (defined on full levels)
    rcld        ,   & ! standard deviation of the saturation deficit    --

!   fields of the radiation
    sohr        ,   & ! rate of solar heating                         ( k/s )
    thhr        ,   & ! rate of thermal heating                       ( k/s )
    clc_sgs     ,   & ! subgrid-scale stratiform cloud cover            --
    alb_rad     ,   & ! albedo of the ground                            --
    sobs        ,   & ! solar radiation at the ground                 ( w/m2)
    thbs        ,   & ! thermal radiation at the ground               ( w/m2)
    pabs        ,   & ! photosynthetic active radiation at the ground ( w/m2)
    sobt        ,   & ! solar radiation at the upper boundary         ( w/m2)
                      ! of the atmosphere
    thbt        ,   & ! thermal radiation at the upper boundary       ( w/m2)
                      ! of the atmosphere
    clch        ,   & ! cloud cover with high clouds                    --   
    clcm        ,   & ! cloud cover with medium clouds                  --   
    clcl        ,   & ! cloud cover with low clouds                     --   
    clct        ,   & ! total cloud cover                               --   
    freshsnow   ,   & ! weighting function indicating 'freshness' of snow in
                      ! upper few centimeters of snow cover            ( -- )
    sun_el      ,   & ! sun elevation angle                           (deg  )
    sun_azi     ,   & ! sun azimuth  angle                            (deg  )

    ! and for the Climate-LM Version
    sodwddm     ,   & ! downward direct solar radiative flux / smu0   ( W/m2)
    qc_rad      ,   & ! subgrid-scale specific cloud water content    (kg/kg)
    qi_rad      ,   & ! subgrid-scale specific ice water content      (kg/kg)

!   fields of the convection
    clc_con     ,   & ! cloud cover due to convection                   --     

!   fields for the radiation correction scheme
    ! these are actual values
    swdir_s     ,   & ! direct comp. of solar radiative flux at surface ( W/m2)
    swdifd_s    ,   & ! diffuse downward comp. of short wave rad. flux  ( W/m2)
    swdifu_s    ,   & ! diffuse upward   comp. of short wave rad. flux  ( W/m2)
    lwd_s       ,   & !         downward comp. of long  wave rad. flux  ( W/m2)
    lwu_s       ,   & !         upward   comp. of long  wave rad. flux  ( W/m2)

    ! this is the essential correction factor
    swdir_cor   ,   & ! direct short wave radiation correction factor actual value

    ! these are topographic parameters
    skyview     ,   & ! sky view
    slo_asp     ,   & ! slope aspect
    slo_ang     ,   & ! slope angle
    horizon     ,   & ! horizon

! 7. fields for model output and diagnostics                          (unit )
! ------------------------------------------
     sod_t      ,   & ! solar downward radiation at top of atmosphere (     )
    asod_t            ! averaged solar downward radiation at top      (     )

! end of data_fields

!------------------------------------------------------------------------------

USE data_runcontrol , ONLY :   &

! 1. start and end of the forecast
! --------------------------------
    ntstep,       & ! actual time step
                    ! indices for permutation of three time levels
    nold  ,       & ! corresponds to ntstep - 1
    nnow  ,       & ! corresponds to ntstep 
    nnew  ,       & ! corresponds to ntstep + 1

! 3. controlling the physics
! --------------------------
    lconv,        & ! forecast with convection
    itype_aerosol,& ! type of aerosol map internal/external
    lemiss,       & ! external emissivity map
    lgsp,         & ! switch for grid-scale cloud and precipitation scheme
    lforest,      & ! if .true., run with forest (evergreen and deciduous)
    lsoil,        & ! forecast with soil model
    lseaice,      & ! forecast with sea ice model
    llake,        & ! forecst with lake model FLake
    l_cosmo_art,  & ! if .TRUE., run the COSMO_ART
    nincrad,      & ! time step increment for running the radiation
    nradcoarse,   & ! radiation coarse-grid number of gpts per hor. direction !T.R.
    lradf_avg,    & ! switch for filtering of radiative increments !T.R.
    nlgw,         & ! number of prognostic soil water levels
    itype_wcld,   & ! type of water cloud diagnosis
    icldm_rad,    & ! mode of cloud representation in radiation  parametr.
    lmulti_layer, & ! run multi-layer soil model
    lmulti_snow , & ! run multi-layer snow model
    lradtopo,     & ! if .TRUE., calculate topographic correction of radiation
    nhori,        & ! number of sectors for the horizont array by the topographic
                    ! correction of the radiation

    ! and for the Climate-LM Version
    ico2_rad,     & ! type of CO2 concentration in radiation parameterization

! 5. additional control variables
! -------------------------------
    ltime,        & ! 
    lreproduce,   & ! the results are reproducible in parallel mode
    l2tls     ,   & ! forecast with 1-TL integration scheme        
    lprog_qi,     & ! if .TRUE., running with cloud ice
    lperi_x,      & ! if lgen=.TRUE.: periodic boundary conditions (.TRUE.) in x-dir.
                    ! or with Davies conditions (.FALSE.)
    lperi_y,      & ! if lgen=.TRUE.: periodic boundary conditions (.TRUE.) in y-dir.
                    ! or with Davies conditions (.FALSE.)
    l2dim,        & ! 2 dimensional runs

! 9. Variables for Ascii file handling, time measuring, ...
! ---------------------------------------------------------
    itype_calendar,&! for specifying the calendar used

! 12. controlling verbosity of debug output
! -----------------------------------------
    idbg_level,   & ! to control the verbosity of debug output
    ldebug_rad,   & ! if .TRUE., debug output for radiation
    lprintdeb_all   ! .TRUE.:  all tasks print debug output
                    ! .FALSE.: only task 0 prints debug output

! end of data_runcontrol 

!------------------------------------------------------------------------------

USE data_soil , ONLY :   &
    csalb     , & !
    csalbw    , & !
    csalb_p   , & !
    csalb_snow, & !
    csalb_snow_min, & ! min. solar albedo of snow for forest free surfaces
    csalb_snow_max, & ! max. solar albedo of snow for forest free surfaces
    csalb_snow_fd , & ! solar albedo of snow for surfaces with deciduous forest
    csalb_snow_fe , & ! solar albedo of snow for surfaces with evergreen forest
    ctalb     , & !
    cdzw12    , & !
    cdzw13    , & !
    cf_snow

! end of data_soil       

!------------------------------------------------------------------------------

USE data_flake, ONLY : &
  ! flake_parameters
  h_Ice_min_flk             , & ! Minimum ice thickness [m]
  tpl_T_f                   , & ! Fresh water freezing point [K]

  ! flake_albedo_ref
  albedo_whiteice_ref       , & ! White ice
  albedo_blueice_ref        , & ! Blue ice
  c_albice_MR                   ! Constant in the interpolation formula for
                                ! the ice albedo (Mironov and Ritter 2004)

!------------------------------------------------------------------------------

USE data_parallel,      ONLY :  &
  nprocx,isubpos, &
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
  imp_integers,  & ! determines the correct INTEGER type used in the model
                   ! for MPI
  sendbuf,       & ! sending buffer for boundary exchange:
                   ! 1-4 are used for sending, 5-8 are used for receiving
  isendbuflen      ! length of one column of sendbuf

!------------------------------------------------------------------------------

USE environment,              ONLY :  &
  exchg_boundaries!,       & ! performs the boundary exchange between
                             ! neighboring processors
!------------------------------------------------------------------------------

USE data_radiation, ONLY : &
    idim_rad,        & ! ie-dimension of the coarser grid
    istartrad,       & ! start- and end-indices for computing the radiation
    iendrad,         & !   (when running on a coarser grid, the input values for
    jstartrad,       & !    fesft are computed on all grid points, to compute an
    jendrad,         & !    average input over several grid points)
    iendparrad,      & ! end-index just for preparations
    jendparrad,      & ! end-index just for preparations

! Imported array data variables with intent (in) for init_radiation:
  jpgas,                   & ! Number of gases
  coali,                   & !
  cobti,                   & !
  pgas,                    & !
  tgas,                    & !
  jpsol   , & ! Number of solar spectral intervals
  jpther  , & ! Number of thermal spectral intervals
  jpspec  , & ! =jpsol+jpther (Total number of spectral intervals)

  solant  , & ! Fraction of solar energy at TOA contained in individual
              ! solar spectral intervals

  planck  , & ! coefficients for the description of the fraction of the total
              ! black body radiation contained in an thermal spectral interval
              ! as a function (2.order polynomial) of temperature
  zketypr , & ! e-type continuum-coefficient for all spectral intervals 
              ! (PA (H2O)**-2) at 296 K
  zketypa , & ! (r) following ROBERTS ET AL. 1976,(a) implicitly derived 
              ! from the AFGL spectral data
  ztetypr , & ! constant for the temperature dependancy of e-type absorption 
              ! for all intervals
  ztetypa , & ! (r) following ROBERTS ET AL. 1976,(a) implicitly derived 
              ! from the AFGL spectral data
  zteref  , & ! reference temperature

  grenze  , & ! Limits of spectral intervals (for information only)

  ! absorption properties of atmospheric gases
  ncgas   , & ! number of coefficients for each interval and gas (maximum=7)
  nfast   , & ! control variable for choice between ESFT/FESFT method in 
              ! each interval
  coai    , & ! weigthing coefficients
  cobi    , & ! absorption coefficients

! array data variables for opt_th and opt_so with intent (in):
  zaea, zaes, zaeg, zaef,         & !
  zlwe, zlwemn, zlwemx,           & !
  zlww, zlwg,                     & !
  ziwe, ziwemn, ziwemx,           & !
  ziww, ziwg,                     & !
  zrsc,                           & !

! for albedo
  rad_csalbw


!------------------------------------------------------------------------------

USE meteo_utilities,    ONLY :  cloud_diag
USE utilities,          ONLY :  get_utc_date

!------------------------------------------------------------------------------

#ifdef COSMOART
USE data_cosmo_art,     ONLY :     &
    Eup          , & ! upward flux (3. band of GRAALS)            (W m-2)
    Edown        , & ! downward flux (3. band of GRAALS)          (W m-2)
    Edir         , & ! direkt intensity (3. band of GRAALS)       (W m-2)
    mmy          , & ! cosinus of sun-zenith-angle in radiation   (1)
    ! CK 20101204 lgas necessary to check if jvals need to be calculated
    lgas         , & ! with gas phase chemistry
    lrad_dust    , & ! mineral dust aerosols
    lrad_seas    , & ! sea salt aerosols
    lrad_aero    , & ! anthropogenic aerosols
    asym_ges     , & !
    asym_seas    , & !
    asym_aero     , & !
    tau_abs_dust  , & !
    tau_streu_dust, & !
    tau_abs_seas  , & !
    tau_streu_seas, & !
    tau_abs_aero  , & !
    tau_streu_aero    !

USE art_papa,           ONLY :  calcjval
#endif

!==============================================================================

IMPLICIT NONE

!==============================================================================

CONTAINS

!==============================================================================
!+ Module procedure in "Radiation" for initializing necessary variables
!------------------------------------------------------------------------------

SUBROUTINE init_radiation

!------------------------------------------------------------------------------
!
! Description:
! The module procedure init_radiation initializes the data necessary for the
! the radiation scheme. It provides the constant aerosol arrays 
! (aersea, aerlan, aerurb and aerdes) and processes the data variables 
! provided by module data_radiation for the calculation of the absorption 
! properties of atmospheric gases.
! The routine is called once at the beginning of a model run in order to
! convert the raw coffeicients from the data variables to those entities used
! in the radiation code.
!
! Method:
! - Computation of the inverse transformations (Legendre and Fourier) of a T10
!   representation of the four aerosol fields
! - scaling of absorption coefficients from the individual reference
!   temperature and pressure to unified conditions, i.e.
!   reference temperature = 273.15 K
!   reference pressure    = 1013.25 hPa
!
!------------------------------------------------------------------------------

! Subroutine arguments: None
! --------------------

! Local arrays and scalars:
! -------------------------
  INTEGER (KIND=iintegers)  ::  &
    i  , j   , jg  , js , jc ,  & ! loop indices
    jzj, jzm1, jzm2, jzm, jzn,  & ! indices for Legendre coefficients
    imn, imnc, imns, jmm, jnn,  & !
    ist

  REAL (KIND=ireals)        ::  &
                               ! arrays for the T10 distrubution of
    zaesc(66) , zaess (55) , & ! sea    type aerosols                     
    zaelc(66) , zaels (55) , & ! land   type aerosols
    zaeuc(66) , zaeus (55) , & ! urban  type aerosols
    zaedc(66) , zaeds (55) , & ! desert type aerosols
    zfaes(21) , zfael (21) , & ! coefficients for spectral
    zfaeu(21) , zfaed (21) , & ! expansion
    zalp (66) ,              & !
    zsinphi   , zcosphi    , & !
    zm, z2m, zre1, ze1, ze2, & !
    zf1m, zf2m, zn, zn2,     & !
    zsin1, zsin2, zsin3, zsin4, zsin5,  & ! 
    zsin6, zsin7, zsin8, zsin9, zsin10, & !
    zcos1, zcos2, zcos3, zcos4, zcos5,  & ! 
    zcos6, zcos7, zcos8, zcos9, zcos10    !

  REAL (KIND=ireals)        ::  &
    zdzwb
 
!- End of header
!==============================================================================

!------------------------------------------------------------------------------
! Section 0: Data for the Fourier coefficients of the four aerosol types          
!------------------------------------------------------------------------------

  DATA zaesc/ &
     +.6688E+00,-.1172E+00,-.1013E+00,+.1636E-01,-.3699E-01,+.1775E-01,      &
                -.9635E-02,+.1290E-02,+.4681E-04,-.9106E-04,+.9355E-04,      &
     -.7076E-01,-.1782E-01,+.1856E-01,+.1372E-01,+.8210E-04,+.2149E-02,      &
                           +.4856E-03,+.2231E-03,+.1824E-03,+.1960E-05,      &
     +.2057E-01,+.2703E-01,+.2424E-01,+.9716E-02,+.1312E-02,-.8846E-03,      &
                                      -.3347E-03,+.6231E-04,+.6397E-04,      &
     -.3341E-02,-.1295E-01,-.4598E-02,+.3242E-03,+.8122E-03,-.2975E-03,      &
                                                 -.7757E-04,+.7793E-04,      &
     +.4455E-02,-.1584E-01,-.2551E-02,+.1174E-02,+.1335E-04,+.5112E-04,      &
                                                            +.5605E-04,      &
     +.7412E-04,+.1857E-02,-.1917E-03,+.4460E-03,+.1767E-04,-.5281E-04,      &
     -.5043E-03,+.2467E-03,-.2497E-03,-.2377E-04,-.3954E-04, &
     +.2666E-03,-.8186E-03,-.1441E-03,-.1904E-04, &
     +.3337E-03,-.1696E-03,-.2503E-04, &
     +.1239E-03,-.9983E-04, &
     -.5283E-04_ireals  /

  DATA zaess/                                                                &
     -.3374E-01,-.3247E-01,-.1012E-01,+.6002E-02,+.5190E-02,+.7784E-03,      &
                           -.1090E-02,+.3294E-03,+.1719E-03,-.5866E-05,      &
     -.4124E-03,-.3742E-01,-.5054E-02,+.3430E-02,+.5513E-03,-.6235E-03,      &
                                      +.2892E-03,-.9730E-04,+.7078E-04,      &
     -.3300E-01,+.5104E-03,-.2156E-02,-.3194E-02,-.5079E-03,-.5517E-03,      &
                                                 +.4632E-04,+.5369E-04,      &
     -.2731E-01,+.5126E-02,+.2241E-02,-.5789E-03,-.3048E-03,-.1774E-03,      &
                                                            +.1946E-05,      &
     -.8247E-02,+.2338E-02,+.1021E-02,+.1575E-04,+.2612E-05,+.1995E-04,      &
     -.1319E-02,+.1384E-02,-.4159E-03,-.2337E-03,+.5764E-04,                 &
     +.1495E-02,-.3727E-03,+.6075E-04,-.4642E-04,                            &
     +.5368E-03,-.7619E-04,+.3774E-04,                                       &
     +.1206E-03,-.4104E-06,                                                  &
     +.2158E-04  /

  DATA zaelc/                                                                &
     +.1542E+00,+.8245E-01,-.1879E-03,+.4864E-02,-.5527E-02,-.7966E-02,      &
                -.2683E-02,-.2011E-02,-.8889E-03,-.1058E-03,-.1614E-04,      &
     +.4206E-01,+.1912E-01,-.9476E-02,-.6780E-02,+.1767E-03,-.5422E-03,      &
                           -.7753E-03,-.2106E-03,-.9870E-04,-.1721E-04,      &
     -.9536E-02,-.9580E-02,-.1050E-01,-.5747E-02,-.1282E-02,+.2248E-03,      &
                                      +.1694E-03,-.4782E-04,-.2441E-04,      &
     +.5781E-03,+.6212E-02,+.1921E-02,-.1102E-02,-.8145E-03,+.2497E-03,      &
                                                 +.1539E-03,-.2538E-04,      &
     -.3993E-02,+.9777E-02,+.4837E-03,-.1304E-02,+.2417E-04,-.1370E-04,      &
                                                            -.3731E-05,      &
     +.1922E-02,-.5167E-03,+.4295E-03,-.1888E-03,+.2427E-04,+.4012E-04,      &
     +.1529E-02,-.2120E-03,+.8166E-04,+.2579E-04,+.3488E-04,                 &
     +.2140E-03,+.2274E-03,-.3447E-05,-.1075E-04,                            &
     -.1018E-03,+.2864E-04,+.3442E-04,                                       &
     -.1002E-03,+.7117E-04,                                                  &
     +.2045E-04  /

  DATA zaels/                                                                &
     +.1637E-01,+.1935E-01,+.1080E-01,+.2784E-02,+.1606E-03,+.1860E-02,      &
                           +.1263E-02,-.2707E-03,-.2290E-03,-.9761E-05,      &
     -.7317E-02,+.2465E-01,+.6799E-02,-.1913E-02,+.1382E-02,+.6691E-03,      &
                                      +.1414E-03,+.3527E-04,-.5210E-04,      &
     +.1873E-01,+.2977E-02,+.4650E-02,+.2509E-02,+.3680E-03,+.1481E-03,      &
                                                 -.6594E-04,-.5634E-04,      &
     +.1592E-01,-.1875E-02,-.1093E-02,+.3022E-03,+.2625E-03,+.3252E-04,      &
                                                            -.3803E-04,      &
     +.4218E-02,-.1843E-02,-.1351E-02,-.2952E-03,-.8171E-05,-.1473E-04,      &
     +.9076E-03,-.1057E-02,+.2676E-03,+.1307E-03,-.3628E-04,                 &
     -.9158E-03,+.4335E-03,+.2927E-04,+.6602E-04,                            &
     -.3570E-03,+.5760E-04,-.3465E-04,                                       &
     -.8535E-04,-.2011E-04,                                                  &
     +.6612E-06  /  

  DATA zaeuc/ &
     +.8005E-01,+.7095E-01,+.2014E-01,-.1412E-01,-.2425E-01,-.1332E-01,      &
                -.2904E-02,+.5068E-03,+.9369E-03,+.4114E-03,+.7549E-04,      &
     +.1922E-01,+.2534E-01,+.2088E-01,+.1064E-01,+.1063E-02,-.2526E-02,      &
                           -.2091E-02,-.9660E-03,-.2030E-03,+.3865E-04,      &
     -.9900E-02,-.5964E-02,+.2223E-02,+.4941E-02,+.3277E-02,+.1038E-02,      &
                                      -.1480E-03,-.2844E-03,-.1208E-03,      &
     +.3999E-02,+.6282E-02,+.2813E-02,+.1475E-02,+.4571E-03,-.1349E-03,      &
                                                 -.9011E-04,-.1936E-04,      &
     +.1994E-02,+.3540E-02,+.8837E-03,+.1992E-03,+.3092E-04,-.7979E-04,      &
                                                            -.2664E-04,      &
     -.5006E-04,+.6447E-03,+.5550E-03,+.1197E-03,+.6657E-04,+.1488E-04,      &
     -.9141E-04,-.2896E-03,-.1561E-03,-.6524E-04,-.1559E-04,                 &
     -.1082E-03,-.4126E-03,-.1732E-03,-.8286E-04,                            &
     -.1993E-04,+.3850E-04,+.2870E-04,                                       &
     +.4493E-04,+.4721E-04,                                                  &
     +.1338E-04  /

  DATA zaeus/                                                                &
     +.6646E-02,+.8373E-02,+.5463E-02,+.4554E-02,+.3301E-02,+.5725E-03,      &
                           -.7482E-03,-.6222E-03,-.2603E-03,-.5127E-04,      &
     -.3849E-04,+.9741E-02,+.8190E-02,+.5712E-02,+.3039E-02,+.5290E-03,      &
                                      -.2044E-03,-.2309E-03,-.1160E-03,      &
     +.9160E-02,+.1286E-01,+.1170E-01,+.5491E-02,+.1393E-02,-.6288E-04,      &
                                                 -.2715E-03,-.1047E-03,      &
     +.4873E-02,+.3545E-02,+.3069E-02,+.1819E-02,+.6947E-03,+.1416E-03,      &
                                                            -.1538E-04,      &
     -.4351E-03,-.1907E-02,-.5774E-03,-.2247E-03,+.5345E-04,+.9052E-04,      &
     -.3972E-04,-.9665E-04,+.7912E-04,-.1094E-04,-.6776E-05,                 &
     +.2724E-03,+.1973E-03,+.6837E-04,+.4313E-04,                            &
     -.7174E-05,+.8527E-05,-.2160E-05,                                       &
     -.7852E-04,+.3453E-06,                                                  &
     -.2402E-05  /

  DATA zaedc/                                                                &
     +.2840E-01,+.1775E-01,-.1069E-01,-.1553E-01,-.3299E-02,+.3583E-02,      &
                +.2274E-02,+.5767E-04,-.3678E-03,-.1050E-03,+.2133E-04,      &
     +.2326E-01,+.1566E-01,-.3130E-02,-.8253E-02,-.2615E-02,+.1247E-02,      &
                           +.1059E-02,+.1196E-03,-.1303E-03,-.5094E-04,      &
     +.1185E-01,+.7238E-02,-.1562E-02,-.3665E-02,-.1182E-02,+.4678E-03,      &
                                      +.4448E-03,+.8307E-04,-.3468E-04,      &
     +.5273E-02,+.3037E-02,-.4014E-03,-.1202E-02,-.4647E-03,+.5148E-04,      &
                                                 +.1014E-03,+.2996E-04,      &
     +.2505E-02,+.1495E-02,+.2438E-03,-.1223E-03,-.7669E-04,-.1638E-04,      &
                                                            +.1869E-05,      &
     +.1094E-02,+.6131E-03,+.1508E-03,+.1765E-04,+.1360E-05,-.7998E-06,      &
     +.4475E-03,+.2737E-03,+.6430E-04,-.6759E-05,-.6761E-05,                 &
     +.1992E-03,+.1531E-03,+.4828E-04,+.5103E-06,                            &
     +.7454E-04,+.5917E-04,+.2152E-04,                                       &
     +.9300E-05,+.9790E-05,                                                  &
     -.8853E-05 /

  DATA zaeds/                                                                &
     +.9815E-02,+.8436E-02,+.1087E-02,-.2717E-02,-.1755E-02,-.1559E-03,      &
                           +.2367E-03,+.8808E-04,+.2001E-05,-.1244E-05,      &
     +.1041E-01,+.8039E-02,+.1005E-02,-.1981E-02,-.1090E-02,+.1595E-05,      &
                                      +.1787E-03,+.4644E-04,-.1052E-04,      &
     +.6593E-02,+.3983E-02,-.1527E-03,-.1235E-02,-.5078E-03,+.3649E-04,      &
                                                 +.1005E-03,+.3182E-04,      &
     +.3225E-02,+.1672E-02,-.7752E-04,-.4312E-03,-.1872E-03,-.1666E-04,      &
                                                            +.1872E-04,      &
     +.1133E-02,+.5643E-03,+.7747E-04,-.2980E-04,-.2092E-04,-.8590E-05,      &
     +.2988E-03,+.6714E-04,-.6249E-05,+.1052E-04,+.8790E-05,                 &
     +.1569E-03,-.1175E-04,-.3033E-04,-.9777E-06,                            &
     +.1101E-03,+.6827E-05,-.1023E-04,                                       &
     +.4231E-04,+.4905E-05,                                                  &
     +.6229E-05  /

 
!------------------------------------------------------------------------------
! Begin Subroutine init_radiation             
!------------------------------------------------------------------------------
 
!------------------------------------------------------------------------------
! Section 1: Compute the start- and end-indices
!------------------------------------------------------------------------------

IF (nradcoarse > 1) THEN
  ! Set indices: 
  !  iendparrad, jendparrad                end of i-/j-index for preparations
  !  istartrad/jstartrad/iendrad/jendrad   begin/end of i-/j-index for averaging

  IF (my_cart_neigh(1) == -1) THEN
    istartrad = 1
  ELSE
    istartrad = istartpar-MOD(isubpos(my_cart_id,1)-1,nradcoarse)
  ENDIF
  iendrad=iend+ABS(MIN(nradcoarse,nboundlines)-MOD(isubpos(my_cart_id,1)          &
                   +iend-MAX(nboundlines,nradcoarse)-1,nradcoarse))
  IF (iend == iendpar) iendrad=iend+MOD(iendrad-iend,nradcoarse)

  IF (my_cart_neigh(4) == -1) THEN
    jstartrad = 1
  ELSE
    jstartrad = jstartpar-MOD(isubpos(my_cart_id,2)-1,nradcoarse)
  ENDIF
  jendrad=jend+ABS(MIN(nradcoarse,nboundlines)-MOD(isubpos(my_cart_id,2)          &
                     +jend-MAX(nboundlines,nradcoarse)-1,nradcoarse))
  IF (jend == jendpar) jendrad=jend+MOD(jendrad-jend,nradcoarse)

  IF (my_cart_neigh(3) == -1) THEN
    iendparrad=iendpar
  ELSE
    iendparrad=iendrad
  ENDIF

  IF (my_cart_neigh(2) == -1) THEN
    jendparrad=jendpar
  ELSE
    jendparrad=jendrad
  ENDIF

  idim_rad = (iendrad-istartrad+nradcoarse)/nradcoarse
  IF (iendpar > iend) idim_rad = idim_rad + 1

ENDIF

!------------------------------------------------------------------------------
! Section 2: Calculation of the inverse Legendre and Fourier transformation  
!------------------------------------------------------------------------------

! loops in i and j over model domain

IF (itype_aerosol == 1) THEN

  zaea=RESHAPE((/0.0477, 0.0875,  0.1198, 0.0458, &
                 0.0387, 0.0439,  0.0599, 0.0396, &
                 0.0381, 0.0129,  0.0130, 0.1304, &
                 0.1757, 0.0949,  0.0653, 0.0795, &
                 0.0962, 0.2046,  0.4116, 0.0169, &
                 0.0204, 0.0263,  0.0348, 0.0361, &
                 0.0030, 0.0271,  0.0613, 0.0118, &
                 0.0160, 0.0231,  0.0287, 0.0127, &
                 0.0103, 0.000016,0.0000, 0.0087, &
                 0.0238, 0.0511,  0.0734, 0.0809/),(/8,5/))

  zaes=RESHAPE((/0.1407, 0.4256,  1.0066, 0.0279, &
                 0.0391, 0.0445,  0.0485, 0.0362, &
                 0.6746, 0.8761,  1.0139, 0.0443, &
                 0.0624, 0.0921,  0.1491, 0.2327, &
                 0.0605, 0.2761,  0.7449, 0.0023, &
                 0.0034, 0.0051,  0.0065, 0.0045, &
                 0.0284, 0.5524,  0.9683, 0.0001, &
                 0.0004, 0.0024,  0.0049, 0.0030, &
                 0.0467, 0.3854,  1.1008, 0.0000, &
                 0.00005,0.0004,  0.0006, 0.0006/),(/8,5/))

  zaeg=RESHAPE((/0.6989, 0.6329,  0.6418, 0.6243, &
                 0.7299, 0.7430,  0.7086, 0.8569, &
                 0.7833, 0.7575,  0.7456, 0.4997, &
                 0.6130, 0.7440,  0.7426, 0.7590, &
                 0.5753, 0.5867,  0.5957, 0.6027, &
                 0.6766, 0.6117,  0.5439, 0.6905, &
                 0.5170, 0.6674,  0.7004, 0.0340, &
                 0.0570, 0.1289,  0.1597, 0.1906, &
                 0.3751, 0.6353,  0.7259, 0.0037, &
                 0.0083, 0.0177,  0.0201, 0.0332/),(/8,5/))

  zaef(:,:)= 0.0_ireals

  DO j = 1, je
    DO i = 1, ie
 
!CPS01
     ! Calculation of the values zalp for the sine of latitude (zsinphi) of the
     ! normalized Legendre associated functions. The limit wave number is 10.
     IF (lperi_x) THEN
       ! In case of periodic BCs, set a constant reference
       ! point for the aerosol distribution to make it equal everywhere
       ! in the domain. This reference point is chosen to be the reference point
       ! of the model domain as determined by pollon, pollat:
       zsinphi  = SIN (degrad*(90.0_ireals-ABS(pollat)))
     ELSE
       zsinphi  = SIN(rlat(i,j) )
     END IF

!CPS01
!CPS  zsinphi  = SIN(rlat(i,j) )
     zcosphi  = SQRT(1.-zsinphi**2)
     jzj      = 2
     zf1m     = SQRT(3.0)
     zalp (1) = 1.0
     zalp (2) = zf1m*zsinphi
     wave_number_loop : DO jzm1 = 1, 11 
       jzm  = jzm1-1
       zm   = jzm
       z2m  = zm + zm
       zre1 = SQRT(z2m+3.0)
       ze1  = 1.0/zre1
       IF (jzm.NE.0) THEN     
          zf2m      = zf1m*zcosphi/SQRT(z2m)
          zf1m      = zf2m*zre1
          jzj       = jzj + 1
          zalp(jzj) = zf2m
          IF(jzm ==10) CYCLE wave_number_loop
          jzj       = jzj + 1
          zalp(jzj) = zf1m*zsinphi
          IF(jzm1==10) CYCLE wave_number_loop
       ENDIF  
       jzm2 = jzm+2
       DO jzn = jzm2, 10
          zn        = jzn
          zn2       = zn**2
          ze2       = SQRT( (4.0*zn2-1.0)/(zn2-zm**2) )
          jzj       = jzj+1
          zalp(jzj) = ze2*(zsinphi*zalp(jzj-1)-ze1*zalp(jzj-2))
          ze1       = 1.0/ze2
       ENDDO  
     ENDDO wave_number_loop

     ! Legendre transform of aerosols
     zfaes(:) = 0.0_ireals
     zfael(:) = 0.0_ireals
     zfaeu(:) = 0.0_ireals
     zfaed(:) = 0.0_ireals

     imn  = 0
     imnc = 0
     imns = 0

     DO jmm = 1, 11
        imn  = imn  + 1
        DO jnn = jmm, 11
           imnc       = imnc + 1
           zfaes(imn) = zfaes(imn)+zalp(imnc)*zaesc(imnc)
           zfael(imn) = zfael(imn)+zalp(imnc)*zaelc(imnc)
           zfaeu(imn) = zfaeu(imn)+zalp(imnc)*zaeuc(imnc)
           zfaed(imn) = zfaed(imn)+zalp(imnc)*zaedc(imnc)
        ENDDO    
        IF(jmm.NE.1) THEN
           imn  = imn+1
           DO jnn = jmm, 11
              imns       = imns + 1
              zfaes(imn) = zfaes(imn)+zalp(imns+11)*zaess(imns)
              zfael(imn) = zfael(imn)+zalp(imns+11)*zaels(imns)
              zfaeu(imn) = zfaeu(imn)+zalp(imns+11)*zaeus(imns)
              zfaed(imn) = zfaed(imn)+zalp(imns+11)*zaeds(imns)
           ENDDO  
        ENDIF
     ENDDO   
 
     ! Inverse Fourier transformation   
!CPS02
     IF (lperi_y .OR. l2dim) THEN
       ! In case of periodic BCs, set a constant reference
       ! point for the aerosol distribution to make it equal everywhere
       ! in the domain. This reference point is chosen to be the reference point
       ! of the model domain as determined by pollon, pollat:
       zcos1   = COS(degrad*(pollon-SIGN(1.0_ireals,pollon)*180.0_ireals))
       zsin1   = sin(degrad*(pollon-SIGN(1.0_ireals,pollon)*180.0_ireals))
     ELSE
       zcos1   = COS(rlon(i,j) )
       zsin1   = SIN(rlon(i,j) )
     END IF

!CPS02
!CPS     zcos1   = COS(rlon(i,j) )
!CPS     zsin1   = SIN(rlon(i,j) )
     zcos2   = zcos1*zcos1 - zsin1*zsin1
     zsin2   = zsin1*zcos1 + zcos1*zsin1
     zcos3   = zcos2*zcos1 - zsin2*zsin1
     zsin3   = zsin2*zcos1 + zcos2*zsin1
     zcos4   = zcos3*zcos1 - zsin3*zsin1
     zsin4   = zsin3*zcos1 + zcos3*zsin1
     zcos5   = zcos4*zcos1 - zsin4*zsin1
     zsin5   = zsin4*zcos1 + zcos4*zsin1
     zcos6   = zcos5*zcos1 - zsin5*zsin1
     zsin6   = zsin5*zcos1 + zcos5*zsin1
     zcos7   = zcos6*zcos1 - zsin6*zsin1
     zsin7   = zsin6*zcos1 + zcos6*zsin1
     zcos8   = zcos7*zcos1 - zsin7*zsin1
     zsin8   = zsin7*zcos1 + zcos7*zsin1
     zcos9   = zcos8*zcos1 - zsin8*zsin1
     zsin9   = zsin8*zcos1 + zcos8*zsin1
     zcos10  = zcos9*zcos1 - zsin9*zsin1
     zsin10  = zsin9*zcos1 + zcos9*zsin1
 
     aersea(i,j) = zfaes(1) + 2.* ( zfaes(2 ) * zcos1 + zfaes(3 ) * zsin1    &
                                  + zfaes(4 ) * zcos2 + zfaes(5 ) * zsin2    &
                                  + zfaes(6 ) * zcos3 + zfaes(7 ) * zsin3    &
                                  + zfaes(8 ) * zcos4 + zfaes(9 ) * zsin4    &
                                  + zfaes(10) * zcos5 + zfaes(11) * zsin5    &
                                  + zfaes(12) * zcos6 + zfaes(13) * zsin6    &
                                  + zfaes(14) * zcos7 + zfaes(15) * zsin7    &
                                  + zfaes(16) * zcos8 + zfaes(17) * zsin8    &
                                  + zfaes(18) * zcos9 + zfaes(19) * zsin9    &
                                  + zfaes(20) * zcos10+ zfaes(21) * zsin10 )

     aerlan(i,j) = zfael(1) + 2.* ( zfael(2 ) * zcos1 + zfael(3 ) * zsin1    &
                                  + zfael(4 ) * zcos2 + zfael(5 ) * zsin2    &
                                  + zfael(6 ) * zcos3 + zfael(7 ) * zsin3    &
                                  + zfael(8 ) * zcos4 + zfael(9 ) * zsin4    &
                                  + zfael(10) * zcos5 + zfael(11) * zsin5    &
                                  + zfael(12) * zcos6 + zfael(13) * zsin6    &
                                  + zfael(14) * zcos7 + zfael(15) * zsin7    &
                                  + zfael(16) * zcos8 + zfael(17) * zsin8    &
                                  + zfael(18) * zcos9 + zfael(19) * zsin9    &
                                  + zfael(20) * zcos10+ zfael(21) * zsin10 )
     aerurb(i,j) = zfaeu(1) + 2.* ( zfaeu(2 ) * zcos1 + zfaeu(3 ) * zsin1    &
                                  + zfaeu(4 ) * zcos2 + zfaeu(5 ) * zsin2    &
                                  + zfaeu(6 ) * zcos3 + zfaeu(7 ) * zsin3    &
                                  + zfaeu(8 ) * zcos4 + zfaeu(9 ) * zsin4    &
                                  + zfaeu(10) * zcos5 + zfaeu(11) * zsin5    &
                                  + zfaeu(12) * zcos6 + zfaeu(13) * zsin6    &
                                  + zfaeu(14) * zcos7 + zfaeu(15) * zsin7    &
                                  + zfaeu(16) * zcos8 + zfaeu(17) * zsin8    &
                                  + zfaeu(18) * zcos9 + zfaeu(19) * zsin9    &
                                  + zfaeu(20) * zcos10+ zfaeu(21) * zsin10 )
     aerdes(i,j) = zfaed(1) + 2.* ( zfaed(2 ) * zcos1 + zfaed(3 ) * zsin1    &
                                  + zfaed(4 ) * zcos2 + zfaed(5 ) * zsin2    &
                                  + zfaed(6 ) * zcos3 + zfaed(7 ) * zsin3    &
                                  + zfaed(8 ) * zcos4 + zfaed(9 ) * zsin4    &
                                  + zfaed(10) * zcos5 + zfaed(11) * zsin5    &
                                  + zfaed(12) * zcos6 + zfaed(13) * zsin6    &
                                  + zfaed(14) * zcos7 + zfaed(15) * zsin7    &
                                  + zfaed(16) * zcos8 + zfaed(17) * zsin8    &
                                  + zfaed(18) * zcos9 + zfaed(19) * zsin9    &
                                  + zfaed(20) * zcos10+ zfaed(21) * zsin10 )

     aersea(i,j) = MAX( 0.0_ireals, MIN( 1.0_ireals, aersea(i,j) ) )
     aerlan(i,j) = MAX( 0.0_ireals, MIN( 1.0_ireals, aerlan(i,j) ) )
     aerurb(i,j) = MAX( 0.0_ireals, MIN( 1.0_ireals, aerurb(i,j) ) )
     aerdes(i,j) = MAX( 0.0_ireals, MIN( 1.0_ireals, aerdes(i,j) ) )
 
    ! end of loops over model domain      
    ENDDO 
  ENDDO       
ENDIF ! itype_aerosol = 1
 
IF (itype_aerosol == 2) THEN
  zaea=RESHAPE((/0.0345_ireals, 0.0511_ireals,  0.0847_ireals, 0.0336_ireals, &
                 0.0499_ireals, 0.0364_ireals,  0.0382_ireals, 0.0260_ireals, &
                 0.0457_ireals, 0.0018_ireals,  0.0015_ireals, 0.1361_ireals, &
                 0.2346_ireals, 0.1177_ireals,  0.0684_ireals, 0.0808_ireals, &
                 0.0707_ireals, 0.0689_ireals,  0.1557_ireals, 0.1258_ireals, &
                 0.1588_ireals, 0.1973_ireals,  0.2766_ireals, 0.1134_ireals, &
                 0.0597_ireals, 0.1077_ireals,  0.2095_ireals, 0.0299_ireals, &
                 0.0456_ireals, 0.0358_ireals,  0.0377_ireals, 0.0304_ireals, &
                 0.0103_ireals, 0.000016_ireals,0.0000_ireals, 0.0087_ireals, &
                 0.0238_ireals, 0.0511_ireals,  0.0734_ireals, 0.0809_ireals/),(/8,5/))

  zaes=RESHAPE((/0.1030_ireals, 0.3977_ireals,  1.0680_ireals, 0.0084_ireals, &
                 0.0142_ireals, 0.0191_ireals,  0.0234_ireals, 0.0140_ireals, &
                 0.7894_ireals, 0.9734_ireals,  1.0110_ireals, 0.0307_ireals, &
                 0.0531_ireals, 0.0546_ireals,  0.0839_ireals, 0.2142_ireals, &
                 0.7157_ireals, 0.8698_ireals,  0.8604_ireals, 0.0645_ireals, &
                 0.0781_ireals, 0.1256_ireals,  0.2317_ireals, 0.1409_ireals, &
                 0.0859_ireals, 0.3442_ireals,  0.9496_ireals, 0.0067_ireals, &
                 0.0113_ireals, 0.0153_ireals,  0.0187_ireals, 0.0113_ireals, &
                 0.0467_ireals, 0.3854_ireals,  1.1008_ireals, 0.0000_ireals, &
                 0.00005_ireals,0.0004_ireals,  0.0006_ireals, 0.0006_ireals/),(/8,5/))

  zaeg=RESHAPE((/0.6562_ireals, 0.6614_ireals,  0.7109_ireals, 0.5043_ireals, &
                 0.6486_ireals, 0.6814_ireals,  0.6489_ireals, 0.7799_ireals, &
                 0.8105_ireals, 0.7906_ireals,  0.7947_ireals, 0.4374_ireals, &
                 0.5203_ireals, 0.7076_ireals,  0.7246_ireals, 0.7535_ireals, &
                 0.6932_ireals, 0.6962_ireals,  0.7402_ireals, 0.4029_ireals, &
                 0.5587_ireals, 0.5618_ireals,  0.4520_ireals, 0.7120_ireals, &
                 0.6462_ireals, 0.6510_ireals,  0.6955_ireals, 0.5041_ireals, &
                 0.6482_ireals, 0.6805_ireals,  0.6477_ireals, 0.7753_ireals, &
                 0.3751_ireals, 0.6353_ireals,  0.7259_ireals, 0.0037_ireals, &
                 0.0083_ireals, 0.0177_ireals,  0.0201_ireals, 0.0332_ireals/),(/8,5/))

  zaef(:,:) = 0.0_ireals
ENDIF ! itype_aerosol = 2

!------------------------------------------------------------------------------
! Section 3: Data for radiative transfer calculations
!------------------------------------------------------------------------------

IF (my_cart_id == 0) THEN
  PRINT *,'*****************************************************'
  PRINT *,'*   Radiative transfer calculations employ data     *'
  PRINT *,'*        provided in routine rad_aibi               *'
  PRINT *,'*****************************************************'
ENDIF

! Include reference pressure and temperature in *cobi*

  DO jg = 1, jpgas
    DO js = 1, jpspec
      DO jc = 1, ncgas(js,jg)
         cobi(jc,js,jg) = cobi(jc,js,jg) * (1./pgas(js,jg))**coali(jc,js,jg) &
                                         * (   tgas(js,jg))**cobti(jc,js,jg)
      ENDDO
    ENDDO
  ENDDO

! security settings, if a gas shall not be considered in an intervall
! where the esft will be used.

  DO jg = 1, jpgas
    DO js = 1, jpspec
      IF ( nfast(js) == 0 .AND. ncgas(js,jg) == 0 ) THEN
         ncgas     (js,jg) = 1
         coai    (1,js,jg) = 1.00
         cobi    (1,js,jg) = 0.00
         coali   (1,js,jg) = 1.00
         cobti   (1,js,jg) = 1.00
      END IF
    ENDDO
  ENDDO

!------------------------------------------------------------------------------
! Section 4: Precalculation of albedo of soil type as function of soil water
!            content and depth of upper soil layer
!------------------------------------------------------------------------------

  ! Albedo of soil type (without vegetation) as function of soil water content
  ! (in mH2O) and depth of the upper soil layer      

  zdzwb    = 0.0_ireals 
  IF (lmulti_layer) THEN
     zdzwb = 1.0_ireals / (2.0_ireals * czmls(1))
  ELSE
    IF ( nlgw == 2 ) THEN
      zdzwb = 1.0_ireals / cdzw12
    ELSE
      zdzwb = 1.0_ireals / cdzw13
    ENDIF
  ENDIF
  DO ist = 1, 10
    rad_csalbw(ist) = csalbw(ist) * zdzwb
  ENDDO

!------------------------------------------------------------------------------
! End of the subroutine 
!------------------------------------------------------------------------------

END SUBROUTINE init_radiation 

!==============================================================================
!==============================================================================
!+ Module procedure in "Radiation" 
!------------------------------------------------------------------------------

SUBROUTINE organize_radiation (ydate_ini)

!------------------------------------------------------------------------------
!
! Description:
!
! The module procedure organize_radiations forms the interface between
! the model and the radiation code adapted from the global model gm_e.
!
! Method:
!
! All variables that are required for the radiation code (i.e. input arrays
! and scalars) are provided or calculated from the model variabels.
! The results are stored as solar and thermal heating rates on the 
! corresponding global arrays sohr and thhr.
!     
!------------------------------------------------------------------------------
!
! T.R.:
! Die Variablen zfltf und zflsf (pfltf pflsf in fesft) wurden entfernt, da sie
! nicht weiterverwemdet werden. lcrf and ldebug wurden von Variablen zu
! Parametern.
! Methode fuer Strahlungsrechnung auf groeberem Gitter:
! Die Eingangsgroessen fuer die Routine fesft werden auf groeberes Gitter gemittelt
! (Variablennamen der gemittelten Groessen: XXX_rn).
! Dann Aufruf von fesft mit _rn-Variablen.
! Dann Zurueckspeichern auf LM-Gitter.
! Dann Abziehen der jeweils auf dem groben Gitter berechneten therm. Ausstrahlung und
! reflektierten Solarstrahlung.
! Dann Filtern der Strahlungsinkremente mit diskretem Filter (analog lconf_avg),
! falls lradf_avg==.TRUE.
! Dann Hinzufuegen der therm. Ausstrahlung (gemaess tg auf feinem Gitter) und der
! solaren Rueckstrahlung (gemaess Albedo auf feinem Gitter).
! Kurzerklaerung zu einigen Variablen:
! nradcoarse: Zahl der GP (je x-/y-Richtung, ueber die gemittelt wird
!          1 = wie bisher, 2 = 2-mal-2-Gebiet usw.
! lradave = interner Schalter, der abfragt, ob ueberhaupt gemittelt werden soll
!           (ist true fuer nradcoarse > 1) ==> bei nradcoarse==1 wird nichts gemacht
! alb_rad am Ende, wie es ausgegeben wird, ist die originale Albedo!
! qc_rad und qi_rad sind ebenfalls originale Felder (urspruengliches LM-Gitter)!
! Version auk9_neu8.
!
! The variables zfltf and zflsf (pfltf,pflsf inside fesft) were removed since they
! are not used. lcrf and ldebug were made parameters.
! Method for radiation calculation on coarser grid:
! Input fields for subroutine fesft are averaged onto coarser grid (variable names
! of coarse-grid variables: XXX_rn).
! Then subroutine fesft is called with _rn variables.
! Then coarse-grid fields are stored back onto original (LM) grid.
! Then surface outgoing thermal radiation (calculated with coarse-grid t_g) and
! reflected solar radiation (calculated with coaerse-grid albedo) are substracted from
! thbs and sobs, resp.
! Then the radiative increments are filtered with discrete filter (analogous to
! lconf_avg), if lradf_avg == .true.
! Then surface outgoing thermal radiation (calculated with coarse-grid t_g) and
! reflected solar radiation (calculated with coaerse-grid albedo) are added to thbs
! and sobs, resp.
!
! nradcoarse: number of gridpoints (per x/y direction) to be averaged
!          1: as hitherto; 2: 2-times-2 area and so on
! lradave = internal switch, whether radiation calculation on coarser grid is applied
! at all (.true. for nradcoarse > 1) ==> if nradcoarse==1, then lradave=.false. and
! nothing new is done
! alb_rad at end, as in output, is original (LM fine-scale) albedo!
! qc_rad and qi_rad are also original (LM fine-scale) fields!

!==============================================================================

! Parameterlist
! -------------

CHARACTER (LEN=10), INTENT(IN)     ::   &
     ydate_ini   ! start of the forecast yyyymmddhh (year, month, day, hour)

! Local parameters: 
! ----------------

  LOGICAL, PARAMETER :: &
    lcrf  = .FALSE. , &
    ldebug= .FALSE.

  REAL    (KIND=ireals   ), PARAMETER ::  &
     zcent  = 0.2500_ireals, & ! centre weight in a nine point stencil !T.R.
     zside  = 0.1250_ireals, & ! weight for side points !T.R.
     zedge  = 0.0625_ireals, & ! weight for edge points !T.R.
     zepclc = 1.0E-8_ireals, & ! avoids cloud cover =1.0 and = 0.0 
     zeph2o = 1.0E-9_ireals, & ! minimum value for specific humidity
     zepemu = 1.0E-9_ireals, & ! avoids cosine of zenith angle = 0.0
     zclwcm = 1.0E-9_ireals, & ! avoids cloud water content = 0.0
     rtod  = 57.2957795_ireals ! conversion from radians to degrees

! the former parameter zqco2 =  0.5014E-3_ireals is now a variable 
! (for specifying different co2 scenarios in the Climate-LM Version).
! It is set later in this subroutine dependent on the setting of the Namelist
! Parameter ico2_rad

! Local scalars:
! -------------

! Input for the radiation routine fesft
! -------------------------------------
  INTEGER (KIND=iintegers)  ::  &
     ki1sd,       & ! start index for first  array dimension
     ki1ed,       & ! end   index for first  array dimension
     ki2sd,       & ! start index for second array dimension
     ki2ed,       & ! end   index for second array dimension
     ki3sd,       & ! start index for third  array dimension
     ki3ed,       & ! end   index for third  array dimension

   ! and the same for the computations
     ki1sc,       & ! start index for first  array computation
     ki1ec,       & ! end   index for first  array computation
     ki2sc,       & ! start index for second array computation
     ki2ec,       & ! end   index for second array computation
     ki3sc,       & ! start index for third  array computation
     ki3ec          ! end   index for third  array computation

  LOGICAL                  ::  &
     lradave,     & ! internal switch whether radiation is coarse-grid ! T.R.
     lrady,       & ! for radiation on coaerse-grid !T.R.
     lsolar(je)      ! control switch for solar calculations

  REAL    (KIND=ireals   ) ::  &
     zstb,        & ! Stefan-Boltzman constant 
     zsct           ! solar constant (at time of year)

  INTEGER (KIND=iintegers) ::  &
    kzdims(24),            &
    j_rn,nradcoarse_y,     & ! for radiation on coaerse grid! T.R.
    izz,ii,jz1, n,         & ! for radiation on coaerse grid! T.R.
    i, j, js, k, i_ld,     & ! loop indices over spatial dimensions
    nzx     , i_std, i_etd,& ! time level of prognostic variables
    nzrad,                 & !
    jj   , itaja,          & ! output from routine get_utc_dat
    ist  ,                 & ! loop index for soil type  
    nztgpk, nzp, nzpa,     & ! loop index for gridpoint output
    izstata,               & ! error status at allocation !T.R.
    izstatd,               & ! error status at deallocation !T.R.
    izerror, izdebug         ! for error status
 
  REAL    (KIND=ireals   ) ::  &
    zalbfak,                   & ! albedo correction factor !T.R.
    zfactor,zfactor_b,         & ! for radiation on coaerse grid !T.R.
    zalbradtopo,               & ! buz/T.R.
    zemissfac,                 & ! buz/T.R.
    zemissivity,               & ! buz/T.R.
    ztrbga, zvobga, zstbga,    & ! output from routine aerdis
    zaeops, zaeopl, zaeopu,    & !   "
    zaeopd, ztrpt , zaeadm,    & !   "
    zstunde,                   & ! output from routine get_utc_dat
    ztwo  , ztho  ,            & !
    zdtzgl, zdek  ,            & !
    zsocof, zeit0 ,            & !
    zdeksin,zdekcos,           & !
    zmaxmu0(je)    ,zeitrad,   & !
    zsinphi,zcosphi, zcosthi,  & !
            zsnow  , zvege  ,  & !
    zph    ,zsigma , zdthdz ,  & !
    zpio   ,zpiu   , zpim   ,  & !
    zpnf   ,zphf   , zphfo  ,  & !
    zthvo  ,zthvu  , zuc    ,  & !
    zclwfac,zclwcmn,           & !
    zclwfs ,zclwcs ,           & !
    zclics ,zclws  ,           & !
    zclick ,zclwck , zclwk,    & !
    zclwfk ,                   & !
    zcs    ,zck    , zfac   ,  & !
    zt_ice1,zt_ice2,           & !
    fgew   ,fgee   , fgqv   ,  & ! name of statement functions
    ztt    ,zzpv   , zzpa   ,  & ! dummy arguments of stat. func.
    zrealdiff, zsmu0_loc    ,  & ! dummy arguments of stat. func.
    zsalb_snow, zsnow_alb   ,  & !
    zqdw, zsex, zf_ice, zdpo,  & !
    zdpn, x1, x2, phi_s

  CHARACTER (LEN=10)   yrad1         ! output from routine get_utc_dat
  CHARACTER (LEN=22)   yrad2         ! output from routine get_utc_dat
  CHARACTER (LEN=80)   yzerrmsg      ! for error message
 
  REAL    (KIND=ireals   ) ::  &
    zyear, zqco2                     ! for specifying different CO2 scenarios
                                     ! in the Climate-LM Version

! Local (automatic) arrays:
! ------------------------

  REAL    (KIND=ireals   ) ::  &

  ! Input for the radiation routine fesft
! zqdw   (ie,ke )     , & ! Total water (qv+qc)
! zsex   (ie,ke )     , & ! 
! zf_ice (ie,ke )     , & !
! zphl   (ie,ke1)     , & ! Pressure at half evels       s
! zdpr   (ie,ke )     , & ! Pressure thickness of layers

  zti    (ie,je,ke1)     , & ! Tempeature at layer boundaries
  zclc   (ie,je,ke )     , & ! Cloud cover in each layer
  zwv    (ie,je,ke )     , & ! Water vapour mixing ratio
  zsw    (ie,je,ke )     , & ! Saturation water vapour mixing ratio over water
  zse    (ie,je,ke )     , & ! Saturation water vapour mixing ratio over ice
  zclwc  (ie,je,ke )     , & ! liquid water mixing ratio
  zciwc  (ie,je,ke )     , & ! ice mixing ratio
  zduco2f(ie,je,ke )     , & ! CO2 content in layer
  zduo3f (ie,je,ke )     , & ! O3  content in layer
  zaeq1  (ie,je,ke )     , & ! Type1-Aerosole optical depth at 0.55 micrometer
  zaeq2  (ie,je,ke )     , & ! Type2  "
  zaeq3  (ie,je,ke )     , & ! Type3  "
  zaeq4  (ie,je,ke )     , & ! Type4  "
  zaeq5  (ie,je,ke )     , & ! Type5  "
  zapre  (ie,je )     , & ! Surface pressure
  zsmu0  (ie,je )     , & ! Cosine of zenith angle 
  zalth  (ie,je )     , & ! Thermal surface albedo 
  zalso  (ie,je )     , & ! Solar surface albedo
#ifdef COUP_OAS_COS
  zpalp  (ie,je )     , & ! Direct albedo
#endif

  ! other values for intermediate storage
  zclcmax(ie,je,ke )     , & !
  zclcmin(ie,je,ke )     , & !
  zclcm1 (ie,je)         , & !
! zclx   (ie,ke )     , & !

  ! Output from the radiation routine fesft
  zflt   (ie,ke1)     , & ! Thermal radiative flux at layer boundary
  zfls   (ie,ke1)     , & ! Solar radiative flux at layer boundary
  zflsdir(ie,ke1)     , & ! solar direct downward radiative flux at 
                          ! layer boundary

  ! surface flux of photosynthetic active radiation and components
  zflpar    (ie )     , & ! surface flux of photosynthetic acive radiation
  zflsp_par (ie )     , & ! direct component
  zflsd_par (ie )     , & ! diffuse downward component
  zflsu_par (ie )     , & ! diffuse upward   component

  ! 2D fields for averaging and distribution, if working on a coarse grid
  zzflsp_par(ie,je)   , & ! direct component
  zzflsd_par(ie,je)   , & ! diffuse downward component
  zzflsu_par(ie,je)       ! diffuse upward component

  REAL    (KIND=ireals   ) ::  &

  ! corrected solar and thermal fluxes at layer boundary and components
  zfls_s    (ie)      , & ! Corrected solar
  zflt_s    (ie)      , & !         thermal
  zflsp     (ie)      , & ! direct component of solar radiative flux
  zflsd     (ie)      , & ! diffuse downward component of solar flux
  zflsu     (ie)      , & ! diffuse upward   component of solar flux
  zfltd     (ie)      , & ! diffuse downward component of thermal flux
  zfltu     (ie)      , & ! diffuse upward   component of thermal flux
  zskyview  (ie,je)       ! used as argument for SR fesft
! zfcor     (ie,je)   , & !
! zslo_ang  (ie)      , & !
! zslo_asp  (ie)      , & !
! zhori     (ie,nhori), & !
! zzsmu0    (ie)      , & !
! zrlat     (ie)      , & !
! zrlon     (ie)

  REAL    (KIND=ireals   ) ::  &

  ! Other local utility arrays
  zqcfo    (ie,je), zqcfn ,                 &
  zo3h     (ie,je),                         &
  zqofo    (ie,je), zqofn,                  &
  zaeqdo   (ie,je), zaeqdn,                 &
  zaequo   (ie,je), zaequn,                 &
  zaeqlo   (ie,je), zaeqln,                 &
  zaeqso   (ie,je), zaeqsn,                 &
  zaetr_top(ie,je), zaetr_bot, zaetr,       &

  ! Constants for vertical distribution of aerosole      
  ! (output from routine aerdis)
  zsign(ke1), zvdaes(ke1), zvdael(ke1),      &
  zvdaeu(ke1), zvdaed(ke1),                  &
  zaeadk(3  ), t_test

! in case of nradcoarse > 1: fields on coarse grid: !T.R.

  REAL (KIND = ireals), ALLOCATABLE :: &
  tg_rn     (:  )      , & ! ground temperature
  tg_ra     (:,:)      , & ! ground temperature

  ! Input for the radiation routine fesft
  zti_rn    (:,:)      , & ! Tempeature at layer boundaries
  zdpr_rn   (:,:)      , & ! Pressure thickness of layers
  zclc_rn   (:,:)      , & ! Cloud cover in each layer
  zwv_rn    (:,:)      , & ! Water vapour mixing ratio
  zsw_rn    (:,:)      , & ! Saturation water vapour mixing ratio over water
  zclwc_rn  (:,:)      , & ! liquid water mixing ratio
  zciwc_rn  (:,:)      , & ! ice mixing ratio
  zduco2f_rn(:,:)      , & ! CO2 content in layer
  zduo3f_rn (:,:)      , & ! O3  content in layer
  zaeq1_rn  (:,:)      , & ! Type1-Aerosole optical depth at 0.55 micrometer
  zaeq2_rn  (:,:)      , & ! Type2  "
  zaeq3_rn  (:,:)      , & ! Type3  "
  zaeq4_rn  (:,:)      , & ! Type4  "
  zaeq5_rn  (:,:)      , & ! Type5  "
  zapre_rn  (:  )      , & ! Surface pressure
  zsmu0_rn  (:  )      , & ! Cosine of zenith angle
  zalth_rn  (:  )      , & ! Thermal surface albedo
#ifdef COUP_OAS_COS
  zpalp_rn  (:  )      , & ! Solar diffuse albedo
#endif
  zalso_rn  (:  )          ! Solar surface albedo

  REAL (KIND = ireals), ALLOCATABLE :: &

  ! Output from the radiation routine fesft
  zflt_rn      (:,:)   , & ! Thermal radiative flux at layer boundary
  zfls_rn      (:,:)   , & ! Solar radiative flux at layer boundary
  zflsdir_rn   (:,:)   , & ! solar direct downward radiative flux at layer boundary

  ! surface flux of photosynthetic active radiation and components
  zflpar_rn    (:  )   , & ! surface flux of photosynthetic acive radiation
  zflsp_par_rn (:  )   , & ! direct component
  zflsu_par_rn (:  )   , & ! diffuse upward component
  zflsd_par_rn (:  )   , & ! diffuse downward component

  ! corrected solar and thermal fluxes at layer boundary and components
  zfls_s_rn    (:  )   , & ! corrected solar radiative flux
  zflt_s_rn    (:  )   , & ! thermal flux
  zflsp_rn     (:  )   , & ! direct component of solar flux at surface
  zflsd_rn     (:  )   , & ! diffuse downward component of solar flux
  zflsu_rn     (:  )   , & ! diffuse upward   component of solar flux
  zfltd_rn     (:  )   , & ! diffuse downward component of thermal flux
  zfltu_rn     (:  )   , & ! diffuse upward   component of thermal flux
  zskyv_rn     (:  )   , & ! diffuse upward   component of thermal flux

  zalb_rad_ori(:,: )   , & ! albedo of original grid
  zsohr      (:,:,:)   , & ! rate of solar heating                     ( k/s )
  zthhr      (:,:,:)   , & ! rate of thermal heating                   ( k/s )
  zsobs        (:,:)   , & ! solar radiation at ground                 ( w/m2)
  zsobt        (:,:)   , & ! solar radiation at upper boundary of atmosphere
  zthbs        (:,:)   , & ! thermal radiation at ground               ( w/m2)
  zthbt        (:,:)   , & ! thermal radiation at upper boundary of atmosphere
  zpabs        (:,:)   , & ! photosynthetic active radiation at ground ( w/m2)
  zsodwddm     (:,:)   , &

  ! 2D fields for averaging and distribution, if working on a coarse grid
  z_zzfltd     (:,:)   , & ! diffuse downward component of thermal flux
  z_zzfltu     (:,:)   , & ! diffuse upward   component of thermal flux
  z_zzflsp     (:,:)   , & ! direct component of solar radiative flux
  z_zzflsd     (:,:)   , & ! diffuse downward component of solar flux
  z_zzflsu     (:,:)   , & ! diffuse upward   component of solar flux

  ! photosynthetic active radiation at the ground ( w/m2): components
  z_zzflsp_par (:,:)   , & ! direct component
  z_zzflsd_par (:,:)   , & ! diffuse downward component
  z_zzflsu_par (:,:)       ! diffuse upward   component

!- End of header
!==============================================================================
 
!------------------------------------------------------------------------------
! Begin Subroutine organize_radiation  
!------------------------------------------------------------------------------
      
! statement function to calculate saturation vapour pressure over water
  fgew(ztt)       = b1 * EXP( b2w*(ztt - b3)/(ztt - b4w) ) ! ztt: temperature

! statement function to calculate saturation vapour pressure over ice
  fgee(ztt)       = b1 * EXP( b2i*(ztt - b3)/(ztt - b4i) ) ! ztt: temperature

! statement function to calculate specific humitdity  
  fgqv(zzpv,zzpa) = rdv*zzpv/(zzpa - o_m_rdv*zzpv)   ! zzpv: vapour pressure
                                                     ! zzpa: total air pressure
 
!------------------------------------------------------------------------------
! Section 1: Initializations
!------------------------------------------------------------------------------

  IF (ldebug_rad) THEN
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

  !----------------------------------------------------------------------------
  ! Section 1.1:  Some preparations and calculation of local utility variables
  !----------------------------------------------------------------------------

  IF (izdebug > 10) THEN
    PRINT *, '        START organize_radiation'
  ENDIF

  ! select time level according to the integration scheme used
  IF ( l2tls ) THEN
    nzx  = nnow
  ELSE
    nzx  = nold
  ENDIF

  !----------------------------------------------------------------------------
  ! Section 1.1a: calculate zenith angle
  !----------------------------------------------------------------------------

  IF (izdebug > 10) THEN
    PRINT *, '           calculate zenith angle'
  ENDIF

  ! Calculation of zenith angle and related quantities                   
  nzrad   = ntstep + nincrad/2-1
  CALL get_utc_date ( nzrad, ydate_ini, dt, itype_calendar, yrad1, yrad2,  &
                      itaja, zstunde )
  READ (yrad1(1:4),'(I4)') jj
  ztwo    = 0.681 + 0.2422*(jj-1949)-(jj-1949)/4
  ztho    = 2.*pi*( REAL(itaja, ireals) -1.0 + ztwo )/365.2422
  zdtzgl  = 0.000075 + 0.001868*COS(   ztho) - 0.032077*SIN(   ztho) &
                     - 0.014615*COS(2.*ztho) - 0.040849*SIN(2.*ztho)
  zdek    = 0.006918 - 0.399912*COS(   ztho) + 0.070257*SIN(   ztho) &
                     - 0.006758*COS(2.*ztho) + 0.000907*SIN(2.*ztho) &
                     - 0.002697*COS(3.*ztho) + 0.001480*SIN(3.*ztho)
  zsocof  = 1.000110 + 0.034221*COS(   ztho) + 0.001280*SIN(   ztho) &
                     + 0.000719*COS(2.*ztho) + 0.000077*SIN(2.*ztho)
  zeit0   = pi*(zstunde-12.)/12. + zdtzgl
  zdeksin = SIN (zdek)
  zdekcos = COS (zdek)
 
  !----------------------------------------------------------------------------
  ! Section 1.1b: choose CO2 scenario
  !----------------------------------------------------------------------------

  IF (izdebug > 10) THEN
    PRINT *, '           choose CO2 scenario   '
  ENDIF

  ! Now set zqco2, dependent on the chosen CO2-Scenario
  zyear = REAL(jj,ireals) + REAL(itaja,ireals)/365.0_ireals

  SELECT CASE (ico2_rad)
  CASE (0)
    ! specific CO2 content of the atmosphere (=330 PPM) (default for DWD)
    zqco2  = 0.5014E-3_ireals

  CASE (1)
    ! time dependent CO2 content (fits of IPCC scenario values, taken from ECHAM5)
    ! A1B scenario (for 1950 <= zyear <= 2100)
    !   only CO2
    zqco2 = (- 2.2915249519766070E+07_ireals                &
             + 45714.032150104744_ireals      * zyear       &
             - 34.178190922262594_ireals      * zyear*zyear &
             + 0.01134997524110402_ireals     * zyear**3    &
             - 1.4124678138498344E-06_ireals  * zyear**4) * 1.519E-06_ireals

  CASE (2)
    ! A1B scenario (for 1950 <= zyear <= 2100)
    !   eff. CO2 (i.e. CO2 & CH4 & N2O)
    zqco2 = ( -2.131843263017098E+07_ireals                &
             + 42697.69425574343_ireals      * zyear       &
             - 32.04969808544885_ireals      * zyear*zyear &
             + 0.010685253016710392_ireals   * zyear**3    &
             - 1.3349801856070718E-06_ireals * zyear**4) * 1.519E-06_ireals

  CASE (3)
    ! B1 scenario (for 1950 <= zyear <= 2100)
    !   only CO2
    zqco2 = (- 1.0401357268181011E+07_ireals               &
             + 21152.707545487563_ireals     * zyear       &
             - 16.116691528852456_ireals     * zyear*zyear &
             + 0.005452554505141226_ireals   * zyear**3    &
             - 6.910849734430986E-07_ireals  * zyear**4) * 1.519E-06_ireals

  CASE (4)
    ! B1 scenario (for 1950 <= zyear <= 2100)
    !   eff. CO2 (i.e. CO2 & CH4 & N2O)
    zqco2 = (- 7.716609874947305E+06_ireals                &
             + 15881.335647116388_ireals     * zyear       &
             - 12.239258629216023_ireals     * zyear*zyear &
             + 0.0041862325463834565_ireals  * zyear**3    &
             - 5.361489502050553E-07         * zyear**4) * 1.519E-06_ireals

  CASE (5)
    ! A2 scenario (for 1950 <= zyear <= 2100)
    !   only CO2
    zqco2 = (  3.682237592956851E06_ireals                 &
             - 7547.069807360021_ireals      * zyear       &
             + 5.8133367065151145_ireals     * zyear*zyear &
             - 0.001994454601121309_ireals   * zyear**3    &
             + 2.571600007798381E-07_ireals  * zyear**4 ) * 1.519E-06_ireals

  CASE (6)
    ! A2 scenario (for 1950 <= zyear <= 2100)
    !   eff. CO2 (i.e. CO2 & CH4 & N2O)
    zqco2 = ( - 340960.0590212098_ireals                    &
              + 403.20639583857496_ireals     * zyear       &
              - 0.074859345260926_ireals      * zyear*zyear &
              - 0.00005743139714985962_ireals * zyear**3    &
              + 1.837122734626407E-08         * zyear**4) * 1.519E-06_ireals

  END SELECT

  !----------------------------------------------------------------------------
  ! Section 1.1c:  initialize background aerosol (aerdis)
  !----------------------------------------------------------------------------

  IF (izdebug > 10) THEN
    PRINT *, '           initialize background aerosol (aerdis)'
  ENDIF

  ! The routine aerdis is called to recieve some parameters for the vertical 
  ! distribution of background aerosol.
  zsign(1) = 0
  DO k = 2, ke1
    zsign(k) = sigmr(k)
  ENDDO
  CALL aerdis ( zsign, zvdaes, zvdael, zvdaeu, zvdaed, ke1,            &
                      ztrbga, zvobga, zstbga, zaeops, zaeopl, zaeopu, &
                      zaeopd, ztrpt , zaeadk, zaeadm)

  !----------------------------------------------------------------------------
  ! Section 1.1d:  setting of boundaries and lradave
  !----------------------------------------------------------------------------

  ! Setting of array boundaries and constant scalar input for routine fesft 
  ki2sd = 1
  ki2ed = 1
  ki3sd = 1
  ki3ed = ke
  ki2sc = 1
  ki2ec = 1
  ki3sc = 1
  ki3ec = ke

  zstb  = sigma  
  zsct  = zsocof*solc

  lradave = nradcoarse > 1

  IF (lradave) THEN

    IF (izdebug > 10) THEN
      PRINT *, '           calculations for radiation averaging  ', nradcoarse, lradave
    ENDIF

    IF (izdebug > 10) THEN
      PRINT *, '           memory allocation'
    ENDIF

    ! Allocate the fields for the coarser grid with idim_rad
    ! (has been computed in init_radiation)
    ALLOCATE ( zti_rn    (idim_rad,ke1) , STAT=izstata )
    ALLOCATE ( zdpr_rn   (idim_rad,ke ) , STAT=izstata )
    ALLOCATE ( zclc_rn   (idim_rad,ke ) , STAT=izstata )
    ALLOCATE ( zwv_rn    (idim_rad,ke ) , STAT=izstata )
    ALLOCATE ( zsw_rn    (idim_rad,ke ) , STAT=izstata )
    ALLOCATE ( zclwc_rn  (idim_rad,ke ) , STAT=izstata )
    ALLOCATE ( zciwc_rn  (idim_rad,ke ) , STAT=izstata )
    ALLOCATE ( zduco2f_rn(idim_rad,ke ) , STAT=izstata )
    ALLOCATE ( zduo3f_rn (idim_rad,ke ) , STAT=izstata )
    ALLOCATE ( zaeq1_rn  (idim_rad,ke ) , STAT=izstata )
    ALLOCATE ( zaeq2_rn  (idim_rad,ke ) , STAT=izstata )
    ALLOCATE ( zaeq3_rn  (idim_rad,ke ) , STAT=izstata )
    ALLOCATE ( zaeq4_rn  (idim_rad,ke ) , STAT=izstata )
    ALLOCATE ( zaeq5_rn  (idim_rad,ke ) , STAT=izstata )
    ALLOCATE ( zapre_rn  (idim_rad    ) , STAT=izstata )
    ALLOCATE ( zsmu0_rn  (idim_rad    ) , STAT=izstata )
    ALLOCATE ( zalth_rn  (idim_rad    ) , STAT=izstata )
    ALLOCATE ( zalso_rn  (idim_rad    ) , STAT=izstata )
#ifdef COUP_OAS_COS
    ALLOCATE ( zpalp_rn  (idim_rad    ) , STAT=izstata )
#endif
    ALLOCATE ( zflt_rn   (idim_rad,ke1) , STAT=izstata )
    ALLOCATE ( zfls_rn   (idim_rad,ke1) , STAT=izstata )
    ALLOCATE ( zflsdir_rn(idim_rad,ke1) , STAT=izstata )
    ALLOCATE ( zflpar_rn (idim_rad    ) , STAT=izstata )
    ALLOCATE ( zflsu_par_rn(idim_rad  ) , STAT=izstata )
    ALLOCATE ( zflsd_par_rn(idim_rad  ) , STAT=izstata )
    ALLOCATE ( zflsp_par_rn(idim_rad  ) , STAT=izstata )
    ALLOCATE ( tg_rn     (ie)           , STAT=izstata )
    ALLOCATE ( tg_ra     (ie,je)        , STAT=izstata )
    ALLOCATE ( zalb_rad_ori(ie,je)      , STAT=izstata )
    ALLOCATE ( zfls_s_rn (idim_rad    ) , STAT=izstata )
    ALLOCATE ( zflt_s_rn (idim_rad    ) , STAT=izstata )
    ALLOCATE ( zflsp_rn  (idim_rad    ) , STAT=izstata )
    ALLOCATE ( zflsd_rn  (idim_rad    ) , STAT=izstata )
    ALLOCATE ( zflsu_rn  (idim_rad    ) , STAT=izstata )
    ALLOCATE ( zfltd_rn  (idim_rad    ) , STAT=izstata )
    ALLOCATE ( zfltu_rn  (idim_rad    ) , STAT=izstata )
    ALLOCATE ( zskyv_rn  (idim_rad    ) , STAT=izstata )

    zti_rn(:,:)     = 0.0_ireals
    zdpr_rn(:,:)    = 0.0_ireals
    zclc_rn(:,:)    = 0.0_ireals
    zwv_rn(:,:)     = 0.0_ireals
    zsw_rn(:,:)     = 0.0_ireals
    zclwc_rn(:,:)   = 0.0_ireals
    zciwc_rn(:,:)   = 0.0_ireals
    zduco2f_rn(:,:) = 0.0_ireals
    zduo3f_rn(:,:)  = 0.0_ireals
    zaeq1_rn(:,:)   = 0.0_ireals
    zaeq2_rn(:,:)   = 0.0_ireals
    zaeq3_rn(:,:)   = 0.0_ireals
    zaeq4_rn(:,:)   = 0.0_ireals
    zaeq5_rn(:,:)   = 0.0_ireals
    zalso_rn(:)     = 0.0_ireals
#ifdef COUP_OAS_COS
    zpalp_rn(:)     = 0.0_ireals
#endif
    zalth_rn(:)     = 0.0_ireals
    zapre_rn(:)     = 0.0_ireals
    zsmu0_rn(:)     = 0.0_ireals

    zfls_s_rn(:)    = 0.0_ireals
    zflt_s_rn(:)    = 0.0_ireals
    zflsp_rn (:)    = 0.0_ireals
    zflsd_rn (:)    = 0.0_ireals
    zflsu_rn (:)    = 0.0_ireals
    zfltd_rn (:)    = 0.0_ireals
    zfltu_rn (:)    = 0.0_ireals
    zskyv_rn (:)    = 0.0_ireals

    ! Setting of array boundaries for routine fesft
    ki1sc=1
    ki1ed=idim_rad
    ki1sd=1

  ELSE !.NOT.lradave:

  IF (izdebug > 10) THEN
    PRINT *, '           settings for no radiation averaging'
  ENDIF

    ! Set zapre for the interface to fesft
!CDIR COLLAPSE
    zapre(:,:) = p0hl(:,:,ke+1)

    istartrad = istartpar
    iendparrad   = iendpar
    jstartrad = jstartpar
    jendparrad   = jendpar

    ! Setting of array boundaries for routine fesft
    ki1sd = 1
    ki1ed = ie
    ki1sc = istartpar
    ki1ec = iendpar

  ENDIF !lradave

! maximum (in-)cloud water content: 0.5% of specific humidity at saturation
  zclwfs=0.005 !0.5% of specific humidity at saturation in stratiform clouds
  zclwfk=0.010 !1.0% of specific humidity at saturation in convective clouds

  !----------------------------------------------------------------------------
  ! Section 1.2: Determine sunshine-condition for every latitude (j-row)
  !----------------------------------------------------------------------------

  IF (izdebug > 10) THEN
    PRINT *, '           determine sunshine conditions'
  ENDIF

  IF ((nprocx > 1) .AND. (lreproduce)) THEN
    zmaxmu0(:) = 0.0_ireals

    i_std = isubpos(my_cart_id,1) - nboundlines
    i_etd = isubpos(my_cart_id,3) + nboundlines

!    DO  js = jstartpar, jendpar
 DO  js = jstartrad, jendparrad
      i_ld  = 0
      DO  i = 1, ie_tot
!CPS03
      IF (lperi_y .OR. l2dim) THEN
          ! Similar thing for lperi_y=.true. or l2dim=.true.:
          ! Use the geogr. latitude of the model reference point:
          zsinphi    = SIN (degrad*(90.0_ireals-ABS(pollat)))
      ELSE
         zsinphi     = SIN (rlattot(i,js))     
      END IF
!
      IF (lperi_x) THEN
      ! In case of lperi_x=.true., use the solar time of the model reference point
      ! (as implied by pollon,pollat) to avoid boundary disturbances.
      ! The problem is that the "true" sun is not periodic, but has to be
      ! "artificially forced" to be periodic for periodic BCs.
        IF (pollat >= 0.0_ireals) THEN
           zeitrad   = zeit0 + degrad*(pollon-SIGN(1.0_ireals,pollon)*180.0_ireals)
        ELSE
           zeitrad   = zeit0 + degrad*pollon
        END IF
     ELSE
          zeitrad      = zeit0 + rlontot(i,js)
     END IF
     zcosphi      = SQRT(1.0_ireals - zsinphi**2)
!CPS03
!CPS         zsinphi      = SIN (rlattot(i,js))
!CPS         zcosphi      = SQRT(1.0_ireals - zsinphi**2)
!CPS         zeitrad      = zeit0 + rlontot(i,js)
         zcosthi      = zdeksin * zsinphi + zdekcos * zcosphi * COS(zeitrad)
         zsmu0_loc    = MAX (zcosthi, zepemu)
         IF ( (i >= i_std) .AND. (i <= i_etd) ) THEN
           i_ld = i_ld + 1
           zsmu0(i_ld,js)  = zsmu0_loc

           ! Sun azimuth and sun elevation (for computation of relative sunshine duration) buz
           sun_el(i_ld,js) = ASIN(zsmu0_loc)
           x1 = zdekcos * SIN(zeitrad) / COS(sun_el(i_ld,js))
           x2 = ( SIN(rlat(i_ld,js)) * zdekcos * COS(zeitrad) - &
               COS(rlat(i_ld,js)) * zdeksin ) / COS(sun_el(i_ld,js))
           IF (x2 < -1.0) x2 = -1.0_ireals
           IF (x2 >  1.0) x2 = 1.0_ireals
           phi_s = ACOS(x2)
           IF (x1 < 0) phi_s = - phi_s
           sun_azi(i_ld,js) = rtod*(phi_s + pi)
           sun_el(i_ld,js) = rtod*sun_el(i_ld,js)
         ENDIF
         zmaxmu0(js)  = MAX (zsmu0_loc, zmaxmu0(js))
      ENDDO
    ENDDO
  ELSE
    zmaxmu0(:) = 0.0_ireals

!    DO  js = jstartpar, jendpar
!      DO  i = istartpar, iendpar
    DO  js = jstartrad, jendparrad
      DO  i = istartrad, iendparrad

!CPS For single processors-------------------------------------------------
     IF (lperi_y .OR. l2dim) THEN
          ! Similar thing for lperi_y=.true. or l2dim=.true.:
          ! Use the geogr. latitude of the model reference point:
          zsinphi    = SIN (degrad*(90.0_ireals-ABS(pollat)))
      ELSE
         zsinphi     = SIN (rlat(i,js))
      END IF
!
      IF (lperi_x) THEN
      ! In case of lperi_x=.true., use the solar time of the model reference point
      ! (as implied by pollon,pollat) to avoid boundary disturbances.
      ! The problem is that the "true" sun is not periodic, but has to be
      ! "artificially forced" to be periodic for periodic BCs.
        IF (pollat >= 0.0_ireals) THEN
           zeitrad   = zeit0 + degrad*(pollon-SIGN(1.0_ireals,pollon)*180.0_ireals)
        ELSE
           zeitrad   = zeit0 + degrad*pollon
        END IF
     ELSE
          zeitrad      = zeit0 + rlon(i,js)
     END IF
     zcosphi      = SQRT(1.0_ireals - zsinphi**2)
!CPS------------------------------------------------

!CPS         zsinphi      = SIN (rlat(i,js))
!CPS         zcosphi      = SQRT(1.0_ireals - zsinphi**2)
!CPS         zeitrad      = zeit0 + rlon(i,js)
         zcosthi      = zdeksin * zsinphi + zdekcos * zcosphi * COS(zeitrad)
         zsmu0(i,js)  = MAX (zcosthi, zepemu)
         zmaxmu0(js)  = MAX (zsmu0(i,js), zmaxmu0(js))
      ENDDO
    ENDDO

    DO  js = jstartrad, jendparrad
      DO  i = istartrad, iendparrad
         ! Sun azimuth and sun elevation (for computation of relative sunshine duration) buz
         sun_el(i,js) = ASIN(zsmu0(i,js))
         x1 = zdekcos * SIN(zeitrad) / COS(sun_el(i,js))
         x2 = ( SIN(rlat(i,js)) * zdekcos * COS(zeitrad) - &
              COS(rlat(i,js)) * zdeksin ) / COS(sun_el(i,js))
         IF (x2 < -1.0) x2 = -1.0_ireals
         IF (x2 >  1.0) x2 = 1.0_ireals
         phi_s = ACOS(x2)
         IF (x1 < 0) phi_s = - phi_s
         sun_azi(i,js) = rtod*(phi_s + pi)
         sun_el(i,js) = rtod*sun_el(i,js)
      ENDDO
    ENDDO

  ENDIF

#ifdef COSMOART
  IF(l_cosmo_art) THEN
    DO  js = jstartpar, jendpar
      DO  i = istartpar, iendpar
        mmy(i,js) = zsmu0(i,js)
      ENDDO
    ENDDO
  ENDIF
#endif

  !----------------------------------------------------------------------------
  ! Section 1.3: Start of loop over the model domain from south to north
  !----------------------------------------------------------------------------

  IF (izdebug > 10) THEN
    PRINT *, '           computation loop over model domain'
  ENDIF

DO j = jstartrad, jendparrad   ! jstartpar, jendpar
  IF (zmaxmu0(j) > zepemu) THEN
    lsolar(j) = .TRUE.
  ELSE
    lsolar(j) = .FALSE.
  ENDIF
ENDDO
  
!------------------------------------------------------------------------------
! Section 2:  Calculation of surface albedo taking soil type,              
!             vegetation and snow/ice conditions into account
!------------------------------------------------------------------------------

DO j = jstartrad, jendparrad      ! jstartpar, jendpar
  DO  i = istartrad, iendparrad   ! istartpar, iendpar
     IF (lemiss) THEN
        zalth(i,j) = 1._ireals-emis_rad(i,j)  ! geographical dependent thermal albedo
     ELSE
       zalth(i,j) = ctalb
     ENDIF

     ist = 10

     ! In the following IF statement, t_snow has been used up to now.
     ! In NetCDF files, t_snow is undefined (-1E20) where no snow exists.
     ! This leads to ice-points over the whole sea. t_g could be used instead,
     ! but this changes the results and has to be tested more intensively.
     ! As an intermediate solution, we use t_snow, where it is defined,
     ! otherwise t_g (in grib-files, t_snow is defined as t_s, where no snow
     ! exists.

     IF(lmulti_snow) THEN
       IF (t_snow_mult(i,j,1,nzx) < 0.0_ireals) THEN
         t_test = t_g   (i,j,nzx)
       ELSE
         t_test = t_snow_mult(i,j,1,nzx)
       ENDIF
     ELSE
       IF (t_snow(i,j,nzx) < 0.0_ireals) THEN
         t_test = t_g   (i,j,nzx)
       ELSE
         t_test = t_snow(i,j,nzx)
       ENDIF
     ENDIF

     IF ( llandmask(i,j) .OR. t_test >= t0_melt - 1.7 ) THEN
       ist = NINT(soiltyp(i,j)) ! water (ist=9) and sea ice (ist=10) included
     ENDIF
     zalso(i,j) = csalb(ist)
     IF (lsoil .AND. llandmask(i,j)) THEN
       IF(lmulti_layer) THEN
         zalso(i,j) = csalb(ist) - rad_csalbw(ist)*w_so(i,j,1,nzx)
       ELSE
         zalso(i,j) = csalb(ist) - rad_csalbw(ist)*w_g1(i,j,nzx)
       ENDIF
     ENDIF  ! lsoil, llandmask
  ENDDO       

  IF (lseaice) THEN
    DO i = istartrad, iendparrad
      ! In case the sea ice model is used AND water point AND ice is present,
      ! compute ice albedo for water points with an empirical formula taken from GME.
      ! The ice albedo is the lower the warmer, and therefore wetter, the ice is.
      ! Use ice temperature at time level nnow (2-time level scheme in sea ice model).
  
      IF ((.NOT. llandmask(i,j)) .AND. (h_ice(i,j,nnow) > 0.0_ireals))               &
      zalso(i,j) = (1.0_ireals-0.3846_ireals*EXP(-0.35_ireals*(t0_melt-t_ice(i,j,nnow)))) &
                 * csalb(10)
    ENDDO
  ENDIF

  IF (llake) THEN
    DO  i = istartrad, iendparrad
      IF((depth_lk(i,j)      >  0.0_ireals) .AND.    &
         (h_ice   (i,j,nnow) >= h_Ice_min_flk) ) THEN
      !  In case the lake model FLake is used AND lake point AND ice is present,
      !  compute ice albedo for lake points with an empirical formulation 
      !  proposed by Mironov and Ritter (2004) for use in GME 
      !  [ice_albedo=function(ice_surface_temperature)].
      !  Use surface temperature at time level "nnow".

        zalso(i,j) = EXP(-c_albice_MR*(tpl_T_f-t_s(i,j,nnow))/tpl_T_f)
        zalso(i,j) = albedo_whiteice_ref * (1._ireals-zalso(i,j)) +      &
                     albedo_blueice_ref  * zalso(i,j)
      ENDIF
    ENDDO
  ENDIF

  ! Snow cover and vegetation
  ! -------------------------

  IF (lsoil) THEN
    DO  i = istartrad, iendparrad    ! istartpar, iendpar
      zvege= 0.0_ireals
      zsnow= 0.0_ireals
      IF (llandmask(i,j)) THEN 
        IF (lmulti_layer) THEN
          ! consider effects of aging on solar snow albedo
          zsalb_snow = csalb_snow_min + &
                       freshsnow(i,j)*(csalb_snow_max-csalb_snow_min)
          IF (lforest) THEN
            zsnow_alb = zsalb_snow*(1._ireals-for_e(i,j)-for_d(i,j))       &
                        + csalb_snow_fe * for_e(i,j)                       &
                        + csalb_snow_fd * for_d(i,j)
          ELSE
            zsnow_alb = zsalb_snow
          ENDIF
        ELSE
          zsnow_alb = csalb_snow
        ENDIF

        ! account for snow cover and plant cover and compute final solar
        ! snow albedo
        zvege = plcov(i,j)
        IF (w_snow(i,j,nzx) > 0.0_ireals)                              &
            zsnow = MIN(1.0_ireals, w_snow(i,j,nzx)/cf_snow)
        zalso(i,j) = zsnow * zsnow_alb +                               &
          (1.0_ireals - zsnow) * (zvege * csalb_p + (1.0_ireals - zvege) * zalso(i,j))
      ENDIF !   llandmask
    ENDDO
  ENDIF
ENDDO

IF (lradave) THEN
  DO j = jstartrad, jendparrad      ! jstartpar, jendpar
    DO  i = istartrad, iendparrad
      zalb_rad_ori (i,j) = zalso (i,j) !T.R. fuer Albedokorrektur
    ENDDO
  ENDDO
ENDIF

#if defined COUP_OAS_COS
!------------------------------------------------------------------------------
! provide Community Land Model surface albedo instead of albedo above
! Only for land points (llandmask)
!------------------------------------------------------------------------------
  DO j = jstartrad, jendparrad      ! jstartpar, jendpar
    DO  i = istartrad, iendparrad
      IF (llandmask(i,j)) THEN
         zpalp (i,j) = alb_rad(i,j,1)
         zalso (i,j) = alb_rad(i,j,2)
         if (lradave) zalb_rad_ori(i,j) = alb_rad(i,j,2)
      ELSE
!CPS need this for parallel radiation
        zpalp  (i,j) = (1.0_ireals +                                           &
         0.5_ireals * ((1.0_ireals/zsmu0(i,j)) * (1.0_ireals/zalso(i,j) - 1.0_ireals))) &
      / (1.0_ireals + ((1.0_ireals/zsmu0(i,j)) * (1.0_ireals/zalso(i,j) - 1.0_ireals)))**2
!CPS need this for parallel radiation
      ENDIF !   llandmask
    ENDDO
  ENDDO
#endif

!------------------------------------------------------------------------------
! Section 3:  Set cloudiness and humidity on input for fesft; 
!             Store cloud cover on corresponding global arrays
!------------------------------------------------------------------------------

  !----------------------------------------------------------------------------
  ! Section 3.1: Calculate water vapour saturation mixing ratios of 
  !              over water and over ice 
  !----------------------------------------------------------------------------

  zt_ice1= t0_melt -  5.0_ireals
  zt_ice2= t0_melt - 25.0_ireals

  DO k = 1, ke
    DO j = jstartrad, jendparrad       ! jstartpar, jendpar
      DO  i = istartrad, iendparrad    ! istartpar, iendpar

        ! specific humidity (zwv) specific total water content (zqdw), 
        ! specific humidity at saturation
        ! over water (zsw ) and ice (zse)
        zph       = p0(i,j,k) + pp(i,j,k,nzx)
        zse (i,j,k) = fgqv ( fgee(t(i,j,k,nzx)), zph)
        zsw (i,j,k) = fgqv ( fgew(t(i,j,k,nzx)), zph)
      ENDDO
    ENDDO
  ENDDO
 
  !----------------------------------------------------------------------------
  ! Section 3.2: Calculate stratiform cloud cover (non-convective)
  !----------------------------------------------------------------------------

  IF     ( icldm_rad == 0 ) THEN

    ! a) No interpretation of clouds at all for radiative calculations
    !-----------------------------------------------------------------

    DO  k = 1, ke
      DO j = jstartrad, jendparrad       ! jstartpar, jendpar
        DO  i = istartrad, iendparrad    ! istartpar, iendpar
          zclwc(i,j,k) = 0.0_ireals
          zciwc(i,j,k) = 0.0_ireals
          zclc (i,j,k) = 0.0_ireals
        ENDDO
      ENDDO
    ENDDO

  ELSEIF ( icldm_rad == 1 ) THEN

    ! b) Only grid-sale water clouds are passed to the radiation routine
    !-------------------------------------------------------------------
   
    DO  k = 1, ke
      DO j = jstartrad, jendparrad       ! jstartpar, jendpar
        DO  i = istartrad, iendparrad    ! istartpar, iendpar
          zclwc(i,j,k) = qc(i,j,k,nzx)
          IF (lprog_qi) THEN
            IF ( qc(i,j,k,nzx)+qi(i,j,k,nnow) > 0.0_ireals ) THEN
              zclc(i,j,k) = 1.0_ireals
            ELSE
              zclc(i,j,k) = 0.0_ireals
            END IF
            zciwc(i,j,k) = qi(i,j,k,nnow)
          ELSE
            IF ( qc(i,j,k,nzx) > 0.0_ireals ) THEN
              zclc(i,j,k) = 1.0_ireals
            ELSE
              zclc(i,j,k) = 0.0_ireals
            END IF
            zciwc(i,j,k) = 0.0_ireals
          ENDIF
          clc_sgs(i,j,k) = zclc(i,j,k)
        ENDDO
      ENDDO
    ENDDO

  ELSEIF (icldm_rad == 2) THEN

    ! c) Cloud cover and water content from statistical diagnosis
    !-------------------------------------------------------------------

    CALL cloud_diag(zclc,zclwc,                                      &
         !istartpar,iendpar,js,js,1,ke,                              &
         istartrad,iendparrad,jstartrad,jendparrad,1,ke,             &
         1,ie, 1,je,1,ke,                                            &
         ie,je,ke,ke1,                                               &
         rdv,o_m_rdv,rvd_m_o,lhocp,t0_melt,                          &
         b1,b2w,b3,b4w,b234w,b2i,b4i,                                &
         uc1,uc2,ucl, clc_diag, q_crit,                              &
         t(:,:,:,nzx),qv(:,:,:,nzx),qc(:,:,:,nzx),                   &
         pp(:,:,:,nzx)+p0(:,:,:),rcld,ps(:,:,nzx),                   &
         itype_wcld)

    DO  k = 1, ke
      DO j = jstartrad, jendparrad       ! jstartpar, jendpar
        DO  i = istartrad, iendparrad    ! istartpar, iendpar
          ! convective (in-)cloud water content 
          ! as a function of specific humidity at saturation
          IF ( t(i,j,k,nzx) >= t0_melt )  THEN
             zclwck = zsw(i,j,k)*zclwfk  !                           
          ELSE
             zclwck = zse(i,j,k)*zclwfk   
          ENDIF

          ! cloud cover of the non convective part of the grid box and cloud ice
          zcs = zclc(i,j,k)
          zciwc(i,j,k) = 0.0_ireals

          IF (lprog_qi) THEN
            ! if there is a grid scale cloud with cloud ice,
            ! even there might has been diagnosed subgrid scale water clouds,
            ! their water is thought to be distributed over the 
            ! whole grid volume:
            IF ( qi(i,j,k,nnow) > 0.0_ireals ) THEN
              zcs = 1.0_ireals
            ENDIF
            zciwc(i,j,k) = qi(i,j,k,nnow)
          ENDIF
          clc_sgs(i,j,k) = zcs

          ! convective cloud cover
          zck = clc_con(i,j,k) 

          ! grid scale cloud cover and water contend
          zclc (i,j,k) = zcs + zck*(1.0-zcs)
          zclwc(i,j,k) = zclwc(i,j,k)*(1.0-zck) + zclwck*zck
        ENDDO
      ENDDO        
    ENDDO        

  ELSEIF ( icldm_rad == 4 .OR. icldm_rad == 3 ) THEN 

    ! a) Standard diagnosis
    ! ---------------------

    DO  k = 1, ke
      DO j = jstartrad, jendparrad       ! jstartpar, jendpar
        DO  i = istartrad, iendparrad    ! istartpar, iendpar
          ! Critical relative humidity as function of thermal stability
          zph      = p0(i,j,k) + pp(i,j,k,nzx)
          zsigma   = zph / ps(i,j,nzx)
          zdthdz   = 0.0

          zsex      = zsw(i,j,k)
          zqdw      = qv(i,j,k,nzx) + qc(i,j,k,nzx)
          IF (lprog_qi) THEN
            zf_ice      = 1.0_ireals - MIN( 1.0_ireals, MAX( 0.0_ireals, &
                          (t(i,j,k,nzx)-zt_ice2)/(zt_ice1-zt_ice2) ) )
            zqdw        = zqdw      + qi(i,j,k,nnow)
            zsex        = zsw(i,j,k) * (1.0_ireals - zf_ice) + zse(i,j,k)*zf_ice
          ENDIF

          IF(k == ke) THEN
            zpio    = ( 1.0E-5 *( p0(i,j,k)+pp(i,j,k,nzx) ) )**rdocp
            zpiu    = ( 1.0E-5 * ps(i,j,nzx)  )**rdocp
            zpim    = 0.5*(zpio+zpiu)
            zthvo   = t  (i,j,k  ,nzx)/zpio
            zthvu   = t_g(i,j,    nzx)/zpiu
            zdthdz  = zthvo - zthvu 
          ELSE IF(zsigma.GT.0.95) THEN
            zpio    = ( 1.0E-5 *( p0(i,j,k  )+pp(i,j,k  ,nzx) ) )**rdocp
            zpiu    = ( 1.0E-5 *( p0(i,j,k+1)+pp(i,j,k+1,nzx) ) )**rdocp
            zpim    = 0.5*(zpio+zpiu)
            zthvo   = t(i,j,k  ,nzx)/zpio
            zthvu   = t(i,j,k+1,nzx)/zpiu
            zdthdz  = zthvo - zthvu + (lh_v*cpdr/zpim)*(qv(i,j,k,nzx)-qv(i,j,k+1,nzx)) 
          ENDIF
 
          ! grid scale cloud cover as function of relative humidity
          zuc     = 0.95 - uc1*zsigma*(1.-zsigma)*(1.+uc2*(zsigma-0.5)) 
          zcs     = MAX ( 0.0_ireals,                &
                    MIN ( 1.0_ireals, (zqdw/zsex-zuc)/(ucl-zuc) ) )**2

          ! Corrections and limitations
          IF ( (zsigma > 0.95_ireals) .AND. (zdthdz < 0.0_ireals) ) THEN
            zcs = 0.0_ireals  ! no cloud cover in unstable stratification
          ENDIF
          IF ( qc(i,j,k,nzx) > 0.0_ireals ) THEN  ! grid scale clouds
            IF ( llandmask(i,j) .AND. k < ke ) zcs = 1.0_ireals
          ENDIF
          IF (lprog_qi) THEN
            IF (qi(i,j,k,nnow) > 1.0E-7_ireals) THEN
              zcs = 1.0_ireals ! grid scale clouds with cloud ice
            ENDIF
          ENDIF

          ! store grid-scale cloud cover on global array
          clc_sgs(i,j,k) = zcs

          ! Maximum in-cloud water content:  1.0% of specific humidity at saturation
          !                                  except for convective clouds (fixed)
          ! Standard diagnosis

          IF (lprog_qi) THEN
            zclws  = 0.005*zsex
            zclwcs = zclws*(1.0_ireals-zf_ice)
            zclics = zclws*zf_ice

            ! Check for grid-scale water or ice-clouds
            ! Now change in zclwcs only if qc(i,j,k,nzx) > 0.0
            IF ( qc(i,j,k,nzx) > 0.0_ireals ) THEN      ! grid scale cloud water
              zclwcs = MAX( zclwcs, 0.5_ireals*qc(i,j,k,nzx) )
            ENDIF
            ! Now change in zclics only if qi(i,j,k,nzx) > 1.0E-7
            IF ( qi(i,j,k,nnow) > 1.0E-7_ireals ) THEN  ! grid scale cloud ice
              zclics = MAX( zclics, 0.5_ireals*qi(i,j,k,nnow) )
            ENDIF

            ! Convective cloud water / ice content
            zclwk  = MAX( 2.0_ireals*zclws, 0.0002_ireals )
            zclwck = zclwk*(1.0_ireals-zf_ice)
            zclick = zclwk*zf_ice

            ! Reduce the cloud cover of ice clouds in the upper troposphere
            ! for the diagnosis of clch and clct
            IF ((k <= klv500) .AND. (zclwcs <= 1.0E-10_ireals) .AND. &
                                    (zclics  > 0.0_ireals) ) THEN
              clc_sgs(i,j,k) = clc_sgs(i,j,k)*MIN( 1._ireals, MAX(0.2_ireals, &
                                 ( LOG(zclics)       - LOG(1.E-7_ireals) )/ &
                                 ( LOG(5.E-5_ireals) - LOG(1.E-7_ireals) )) )
            ENDIF

            ! set area-average cloud water/ice content
            zclwc(i,j,k) = zclwck*clc_con(i,j,k) + &
                           zclwcs*clc_sgs(i,j,k)*(1.0_ireals-clc_con(i,j,k))
            zciwc(i,j,k) = zclick*clc_con(i,j,k) + &
                           zclics*clc_sgs(i,j,k)*(1.0_ireals-clc_con(i,j,k))

          ELSE

            zclwcs = 0.005*zsw(i,j,k)
            zclwck = MAX( zclwcs, 0.0002_ireals )
            IF ( qc(i,j,k,nzx) > 0.0 ) THEN  ! grid scale clouds
               zclwcs = MAX( zclwcs, 0.5*qc(i,j,k,nzx) )
            ENDIF

            ! set area-average cloud water/ice content
            zclwc(i,j,k) = zclwck*clc_con(i,j,k) + &
                           zclwcs*clc_sgs(i,j,k)*(1.0_ireals-clc_con(i,j,k))

            ! set area average cloud ice content (in-cloud)
            zciwc(i,j,k) = 0.0_ireals

          ENDIF

          ! calculate combined cloud cover
          zclc (i,j,k) = clc_sgs(i,j,k) + &
                         clc_con(i,j,k)*( 1.0_ireals - clc_sgs(i,j,k) )

        ENDDO
      ENDDO        
    ENDDO        

  ENDIF  ! icldm_rad

  ! Restrictions for radiative calculations
  ! ------------------------------------

  DO  k = 1, ke
    DO j = jstartrad, jendparrad       ! jstartpar, jendpar
      DO  i = istartrad, iendparrad    ! istartpar, iendpar
        zwv  (i,j,k) = MIN( MAX(zeph2o,qv(i,j,k,nzx)), zsw(i,j,k) )
        zclc (i,j,k) = MAX( zepclc, MIN(1.0_ireals-zepclc,zclc(i,j,k)) )
        zclwc(i,j,k) = MAX( zclwcm, zclwc(i,j,k) )
        zciwc(i,j,k) = MAX( zclwcm, zciwc(i,j,k) )

        ! set qc_rad, qi_rad
        qc_rad(i,j,k) = zclwc(i,j,k)
        qi_rad(i,j,k) = zciwc(i,j,k)
      ENDDO
    ENDDO
  ENDDO        

  !----------------------------------------------------------------------------
  ! Section 3.3:  Calculate and store total cloud cover for 3 integral layers
  !               (high, medium, low)
  !----------------------------------------------------------------------------

  DO j = jstartpar, jendpar
    DO i = istartpar, iendpar
      zclcm1(i,j) = 1.0 - zclc(i,j,1)
    ENDDO
  ENDDO
 
  DO j = jstartpar, jendpar
!CDIR UNROLL=10
    DO  k  = 2, ke
      DO i = istartpar, iendpar
        zclcmax(i,j,k) = 1.0_ireals - MAX(zclc(i,j,k), zclc(i,j,k-1))
        zclcmin(i,j,k) = 1.0_ireals / (1.0_ireals - zclc(i,j,k-1))
      ENDDO
    ENDDO
  ENDDO

  DO j = jstartpar, jendpar
!CDIR UNROLL=10
    DO  k  = 2, klv400
      DO i = istartpar, iendpar
        zclcm1(i,j) = zclcm1(i,j) * zclcmax(i,j,k) * zclcmin(i,j,k)
      ENDDO
    ENDDO
  ENDDO
 
  DO j = jstartpar, jendpar
    DO i = istartpar, iendpar
      clch  (i,j) = 1.0 - zclcm1(i,j) - zepclc
    ENDDO
  ENDDO
 
  DO j = jstartpar, jendpar
!CDIR UNROLL=10
    DO  k  = klv400+1, ke
      DO i = istartpar, iendpar
        zclcm1(i,j) = zclcm1(i,j)* zclcmax(i,j,k) * zclcmin(i,j,k)
      ENDDO
    ENDDO   
  ENDDO   

  DO j = jstartpar, jendpar
    DO i = istartpar, iendpar
      clct  (i,j) = 1.0 - zclcm1(i,j) - zepclc
      zclcm1(i,j) = 1.0 - zclc(i,j,klv400+1)
    ENDDO
  ENDDO

  DO j = jstartpar, jendpar
!CDIR UNROLL=10
    DO  k  = klv400+2,klv800
      DO i = istartpar, iendpar
        zclcm1(i,j) = zclcm1(i,j) * zclcmax(i,j,k) * zclcmin(i,j,k)
      ENDDO
    ENDDO   
  ENDDO   

  DO j = jstartpar, jendpar
    DO i = istartpar, iendpar
      clcm  (i,j) = 1.0 - zclcm1(i,j) - zepclc
      zclcm1(i,j) = 1.0 - zclc(i,j,klv800+1)
    ENDDO
  ENDDO

  DO j = jstartpar, jendpar
!CDIR UNROLL=10
    DO  k  = klv800+2,ke
      DO i= istartpar, iendpar
        zclcm1(i,j) = zclcm1(i,j) * zclcmax(i,j,k) * zclcmin(i,j,k)
      ENDDO
    ENDDO   
  ENDDO   

  DO j = jstartpar, jendpar
    DO i= istartpar, iendpar
      clcl  (i,j) = 1.0 - zclcm1(i,j) - zepclc
    ENDDO
  ENDDO
 
!------------------------------------------------------------------------------
! Section 4:  Set pressure and temperature on input for fesft; 
!------------------------------------------------------------------------------

  ! Surface pressure, half level pressure and pressure thickness. 
  ! At present, pressure is replaced by model reference pressure 
  ! for radiation calculations

  ! The following local variables have been replaced by the global variables
  ! DO  i = istartrad, iendparrad         ! istartpar, iendpar
  !   zapre(i) = p0hl(i,js,ke+1)
  !   zphl(i,ke+1) = zapre(i)
  ! ENDDO
  ! DO k = 1, ke
  !   DO  i = istartrad, iendparrad         ! istartpar, iendpar
  !    zdpr(i,k) = dp0(i,js,k)
  !    zphl(i,k) = p0hl(i,js,k)
  !   ENDDO
  ! ENDDO

  ! Temperatures at layer boundaries
  DO  k = 2, ke
    DO j = jstartrad, jendparrad       ! jstartpar, jendpar
      DO  i = istartrad, iendparrad    ! istartpar, iendpar
        zpnf     = p0hl(i,j,k  )
        zphf     = p0  (i,j,k  )
        zphfo    = p0  (i,j,k-1)
          zti(i,j,k) = (t(i,j,k-1,nzx)*zphfo*(zphf - zpnf )      &
                      + t(i,j,k  ,nzx)*zphf *(zpnf - zphfo) )    &
                                 * (1./(zpnf *(zphf - zphfo)))
      ENDDO
    ENDDO
  ENDDO         
 
  DO j = jstartrad, jendparrad       ! jstartpar, jendpar
    DO  i = istartrad, iendparrad       ! istartpar, iendpar
      zpnf    = p0hl  (i,j,2)
      zphf    = p0    (i,j,1)
      zti(i,j,ke1) = t_g(i,j,nzx)
      zti(i,j,  1) = t  (i,j,1,nzx) - zphf*(t(i,j,1,nzx)-zti(i,j,2))/(zphf - zpnf)
    ENDDO
  ENDDO
 
!------------------------------------------------------------------------------
! Section 5:  Calculate amounts of absorbers (CO2, O3, Aerosol)
!------------------------------------------------------------------------------

#ifdef COSMOART
! Change climatology for dust
  IF (l_cosmo_art .AND. lrad_dust) THEN
    DO j = jstartrad, jendparrad       ! jstartpar, jendpar
      DO  i = istartrad, iendparrad       ! istartpar, iendpar
        IF (itype_aerosol == 1) THEN
          aerdes(i,j) = 0.0_ireals
        ELSEIF (itype_aerosol == 2) THEN
          aer_du(i,j) = 0.0_ireals
        ENDIF
      ENDDO
    ENDDO
  ENDIF
#endif

  IF (itype_aerosol == 1) THEN
    DO j = jstartrad, jendparrad       ! jstartpar, jendpar
      DO  i = istartrad, iendparrad       ! istartpar, iendpar
        zdpo           = p0hl(i,j,1)
        zo3h     (i,j) = SQRT(hmo3(i,j))**3
        zqcfo    (i,j) = zqco2*zdpo
        zaeqso   (i,j) = zaeops*aersea(i,j)*zvdaes(1)
        zaeqlo   (i,j) = zaeopl*aerlan(i,j)*zvdael(1)
        zaequo   (i,j) = zaeopu*aerurb(i,j)*zvdaeu(1)
        zaeqdo   (i,j) = zaeopd*aerdes(i,j)*zvdaed(1)
        zaetr_top(i,j) = 1.0_ireals
        zqofo    (i,j) = vio3(i,j)*SQRT(zdpo**3)/(SQRT(zdpo**3) + zo3h(i,j))
      ENDDO
    ENDDO
  ELSEIF (itype_aerosol == 2) THEN
    ! new Tegen aerosol climatology: no multiplication with tau(max) as climatology not normalised!
    DO j = jstartrad, jendparrad       ! jstartpar, jendpar
      DO  i = istartrad, iendparrad       ! istartpar, iendpar
        zdpo           = p0hl(i,j,1)
        zo3h     (i,j) = SQRT(hmo3(i,j))**3
        zqcfo    (i,j) = zqco2*zdpo
        zaeqso   (i,j) =  aer_ss(i,j)              *zvdaes(1)
        zaeqlo   (i,j) =( aer_or(i,j)+aer_su(i,j) )*zvdael(1)
        zaequo   (i,j) =  aer_bc(i,j)              *zvdaeu(1)
        zaeqdo   (i,j) =  aer_du(i,j)              *zvdaed(1)
        zaetr_top(i,j) =  1.0_ireals
        zqofo    (i,j) = vio3(i,j)*SQRT(zdpo**3)/(SQRT(zdpo**3) + zo3h(i,j))
      ENDDO
    ENDDO
  ENDIF

  IF (itype_aerosol == 1) THEN
    DO j = jstartrad, jendparrad       ! jstartpar, jendpar
!CDIR UNROLL=10
      DO  k = 1, ke
        DO  i = istartrad, iendparrad       ! istartpar, iendpar
          zdpn           = p0hl(i,j,k+1)
          zaeqsn         = zaeops*aersea(i,j)*zvdaes(k+1)
          zaeqln         = zaeopl*aerlan(i,j)*zvdael(k+1)
          zaequn         = zaeopu*aerurb(i,j)*zvdaeu(k+1)
          zaeqdn         = zaeopd*aerdes(i,j)*zvdaed(k+1)
          zaetr_bot      = zaetr_top(i,j) * ( MIN (1.0_ireals, zti(i,j,k)/zti(i,j,k+1)) )**ztrpt
          zqcfn          = zqco2 * zdpn
          zqofn          = vio3(i,j)*SQRT(zdpn**3)/(SQRT(zdpn**3) + zo3h(i,j))
          zduco2f(i,j,k) = zqcfn-zqcfo(i,j)
          zduo3f (i,j,k) = zqofn-zqofo(i,j)
          zaetr          = SQRT(zaetr_bot*zaetr_top(i,j))
          zaeq1(i,j,k)   = (1.-zaetr) * (ztrbga*dp0(i,j,k)+zaeqln-zaeqlo(i,j)+zaeqdn-zaeqdo(i,j))
          zaeq2(i,j,k)   = (1.-zaetr) * ( zaeqsn-zaeqso(i,j) )
          zaeq3(i,j,k)   = (1.-zaetr) * ( zaequn-zaequo(i,j) )
          zaeq4(i,j,k)   =     zaetr  *   zvobga*dp0(i,j,k)
          zaeq5(i,j,k)   =     zaetr  *   zstbga*dp0(i,j,k)

          zqcfo(i,j)     = zqcfn
          zqofo(i,j)     = zqofn
          zaetr_top(i,j) = zaetr_bot
          zaeqso(i,j)    = zaeqsn
          zaeqlo(i,j)    = zaeqln
          zaequo(i,j)    = zaequn
          zaeqdo(i,j)    = zaeqdn
        ENDDO
      ENDDO
    ENDDO
  ELSEIF (itype_aerosol == 2) THEN
    ! new Tegen aerosol climatology: no multiplication with tau(max) as climatology not normalised!
    DO j = jstartrad, jendparrad       ! jstartpar, jendpar
!CDIR UNROLL=10
      DO  k = 1, ke
        DO  i = istartrad, iendparrad       ! istartpar, iendpar
          zdpn           = p0hl(i,j,k+1)
          zaeqsn         =   aer_ss(i,j)              * zvdaes(k+1)
          zaeqln         =  (aer_or(i,j)+aer_su(i,j)) * zvdael(k+1)
          zaequn         =   aer_bc(i,j)              * zvdaeu(k+1)
          zaeqdn         =   aer_du(i,j)              * zvdaed(k+1)
          zaetr_bot      = zaetr_top(i,j) * ( MIN (1.0_ireals, zti(i,j,k)/zti(i,j,k+1)) )**ztrpt
          zqcfn          = zqco2 * zdpn
          zqofn          = vio3(i,j)*SQRT(zdpn**3)/(SQRT(zdpn**3) + zo3h(i,j))
          zduco2f(i,j,k) = zqcfn-zqcfo(i,j)
          zduo3f (i,j,k) = zqofn-zqofo(i,j)
          zaetr          = SQRT(zaetr_bot*zaetr_top(i,j))

          zaeq1(i,j,k)   = (1.0_ireals-zaetr)*( ztrbga*dp0(i,j,k) + zaeqln - zaeqlo(i,j) )
          zaeq2(i,j,k)   = (1.0_ireals-zaetr)*(zaeqsn-zaeqso(i,j))
          zaeq3(i,j,k)   = (1.0_ireals-zaetr)*(zaeqdn-zaeqdo(i,j))
          zaeq4(i,j,k)   = (1.0_ireals-zaetr)*(zaequn-zaequo(i,j))
          zaeq5(i,j,k) =    zaetr * zstbga*dp0(i,j,k)

          zqcfo(i,j)     = zqcfn
          zqofo(i,j)     = zqofn
          zaetr_top(i,j) = zaetr_bot
          zaeqso(i,j)    = zaeqsn
          zaeqlo(i,j)    = zaeqln
          zaequo(i,j)    = zaequn
          zaeqdo(i,j)    = zaeqdn
        ENDDO
      ENDDO
    ENDDO
  ENDIF

!------------------------------------------------------------------------------
! Section 6:  Correction factors for radiation in complex topography
!------------------------------------------------------------------------------

  IF (lradtopo) THEN
    IF (izdebug > 10) THEN
      PRINT *,'        organize_radiation with lradtopo = ', lradtopo
    ENDIF

    !US:  Careful: this does NOT work with nradcoarse > 1!!!
    CALL calc_rad_corrections (slo_ang, slo_asp, horizon, zsmu0,       &
              rlat, rlon, zdeksin, zdekcos, zeit0, swdir_cor,          &
              ie, je, nhori, istartpar, iendpar, jstartpar, jendpar,   &
              izdebug)

    ! and set zskyview (1-dimensional slice)
    DO j = jstartrad, jendparrad       ! jstartpar, jendpar
      DO  i = istartrad, iendparrad       ! istartpar, iendpar
        zskyview(i,j) = skyview(i,j)
      ENDDO
    ENDDO
  ELSE
    ! Set default value for skyview
    DO j = jstartrad, jendparrad       ! jstartpar, jendpar
      DO  i = istartrad, iendparrad       ! istartpar, iendpar
        zskyview(i,j) = 1.0_ireals
      ENDDO
    ENDDO
  ENDIF

!------------------------------------------------------------------------------
! Section 7:  Average input values for fesft on nradcoarse**2 gridpoints
!------------------------------------------------------------------------------

loop_south_north: DO  js = jstartrad, jendparrad       ! jstartpar, jendpar

  IF (lradave) THEN

    IF (jendpar > jend .AND. js > jendrad)  THEN
      lrady=.FALSE.
      IF ( js == jendpar ) THEN
        lrady=.TRUE.
        nradcoarse_y=jendpar-jendrad
      ENDIF
    ELSE
      j_rn=MOD(-jstartrad+js+1,nradcoarse)
      lrady=.FALSE.
      IF (MOD(j_rn,nradcoarse) == 0) THEN
        IF (js > nradcoarse-1) THEN
          lrady=.TRUE.
          nradcoarse_y=nradcoarse
        ENDIF
      ENDIF
    ENDIF

    IF (lrady) THEN
      zfactor=1.0_ireals/(REAL(nradcoarse)*REAL(nradcoarse_y))
      IF ( iendpar > iend .AND. iendpar > iendrad ) THEN
        zfactor_b=1.0_ireals/(REAL(nradcoarse_y)*REAL(iendpar-iendrad))
      ENDIF
    ENDIF

    IF (lrady) THEN

      ! Atmosphere
      izz = (iendrad-istartrad+nradcoarse) / nradcoarse

!US   IF (js < jendpar) THEN
      IF (nradcoarse_y == 2) THEN
        DO k=1,ke
          DO n = 1, izz
            i = istartrad + (n-1) * nradcoarse
            zti_rn    (n,k) = zti_rn    (n,k) + zti    (i,js-1,k) + zti    (i+1,js-1,k) + zti    (i,js,k) + zti    (i+1,js,k)
            zdpr_rn   (n,k) = zdpr_rn   (n,k) + dp0    (i,js-1,k) + dp0    (i+1,js-1,k) + dp0    (i,js,k) + dp0    (i+1,js,k)
            zclc_rn   (n,k) = zclc_rn   (n,k) + zclc   (i,js-1,k) + zclc   (i+1,js-1,k) + zclc   (i,js,k) + zclc   (i+1,js,k)
            zwv_rn    (n,k) = zwv_rn    (n,k) + zwv    (i,js-1,k) + zwv    (i+1,js-1,k) + zwv    (i,js,k) + zwv    (i+1,js,k)
            zsw_rn    (n,k) = zsw_rn    (n,k) + zsw    (i,js-1,k) + zsw    (i+1,js-1,k) + zsw    (i,js,k) + zsw    (i+1,js,k)
            zclwc_rn  (n,k) = zclwc_rn  (n,k) + zclwc  (i,js-1,k) + zclwc  (i+1,js-1,k) + zclwc  (i,js,k) + zclwc  (i+1,js,k)
            zciwc_rn  (n,k) = zciwc_rn  (n,k) + zciwc  (i,js-1,k) + zciwc  (i+1,js-1,k) + zciwc  (i,js,k) + zciwc  (i+1,js,k)
            zduco2f_rn(n,k) = zduco2f_rn(n,k) + zduco2f(i,js-1,k) + zduco2f(i+1,js-1,k) + zduco2f(i,js,k) + zduco2f(i+1,js,k)
            zduo3f_rn (n,k) = zduo3f_rn (n,k) + zduo3f (i,js-1,k) + zduo3f (i+1,js-1,k) + zduo3f (i,js,k) + zduo3f (i+1,js,k)
            zaeq1_rn  (n,k) = zaeq1_rn  (n,k) + zaeq1  (i,js-1,k) + zaeq1  (i+1,js-1,k) + zaeq1  (i,js,k) + zaeq1  (i+1,js,k)
            zaeq2_rn  (n,k) = zaeq2_rn  (n,k) + zaeq2  (i,js-1,k) + zaeq2  (i+1,js-1,k) + zaeq2  (i,js,k) + zaeq2  (i+1,js,k)
            zaeq3_rn  (n,k) = zaeq3_rn  (n,k) + zaeq3  (i,js-1,k) + zaeq3  (i+1,js-1,k) + zaeq3  (i,js,k) + zaeq3  (i+1,js,k)
            zaeq4_rn  (n,k) = zaeq4_rn  (n,k) + zaeq4  (i,js-1,k) + zaeq4  (i+1,js-1,k) + zaeq4  (i,js,k) + zaeq4  (i+1,js,k)
            zaeq5_rn  (n,k) = zaeq5_rn  (n,k) + zaeq5  (i,js-1,k) + zaeq5  (i+1,js-1,k) + zaeq5  (i,js,k) + zaeq5  (i+1,js,k)
          ENDDO ! n
        ENDDO   ! k
      ELSE ! nradcoarse_y == 1
        DO k=1,ke
          DO n = 1, izz
            i = istartrad + (n-1) * nradcoarse
            zti_rn    (n,k) = zti_rn    (n,k) + zti    (i,js,k) + zti    (i+1,js,k)
            zdpr_rn   (n,k) = zdpr_rn   (n,k) + dp0    (i,js,k) + dp0    (i+1,js,k)
            zclc_rn   (n,k) = zclc_rn   (n,k) + zclc   (i,js,k) + zclc   (i+1,js,k)
            zwv_rn    (n,k) = zwv_rn    (n,k) + zwv    (i,js,k) + zwv    (i+1,js,k)
            zsw_rn    (n,k) = zsw_rn    (n,k) + zsw    (i,js,k) + zsw    (i+1,js,k)
            zclwc_rn  (n,k) = zclwc_rn  (n,k) + zclwc  (i,js,k) + zclwc  (i+1,js,k)
            zciwc_rn  (n,k) = zciwc_rn  (n,k) + zciwc  (i,js,k) + zciwc  (i+1,js,k)
            zduco2f_rn(n,k) = zduco2f_rn(n,k) + zduco2f(i,js,k) + zduco2f(i+1,js,k)
            zduo3f_rn (n,k) = zduo3f_rn (n,k) + zduo3f (i,js,k) + zduo3f (i+1,js,k)
            zaeq1_rn  (n,k) = zaeq1_rn  (n,k) + zaeq1  (i,js,k) + zaeq1  (i+1,js,k)
            zaeq2_rn  (n,k) = zaeq2_rn  (n,k) + zaeq2  (i,js,k) + zaeq2  (i+1,js,k)
            zaeq3_rn  (n,k) = zaeq3_rn  (n,k) + zaeq3  (i,js,k) + zaeq3  (i+1,js,k)
            zaeq4_rn  (n,k) = zaeq4_rn  (n,k) + zaeq4  (i,js,k) + zaeq4  (i+1,js,k)
            zaeq5_rn  (n,k) = zaeq5_rn  (n,k) + zaeq5  (i,js,k) + zaeq5  (i+1,js,k)
          ENDDO ! n
        ENDDO   ! k
      ENDIF ! nradcoarse_y
      DO k = 1, ke
        DO n = 1, izz
          zti_rn    (n,k) = zti_rn    (n,k) * zfactor
          zdpr_rn   (n,k) = zdpr_rn   (n,k) * zfactor
          zclc_rn   (n,k) = zclc_rn   (n,k) * zfactor
          zwv_rn    (n,k) = zwv_rn    (n,k) * zfactor
          zsw_rn    (n,k) = zsw_rn    (n,k) * zfactor
          zclwc_rn  (n,k) = zclwc_rn  (n,k) * zfactor
          zciwc_rn  (n,k) = zciwc_rn  (n,k) * zfactor
          zduco2f_rn(n,k) = zduco2f_rn(n,k) * zfactor
          zduo3f_rn (n,k) = zduo3f_rn (n,k) * zfactor
          zaeq1_rn  (n,k) = zaeq1_rn  (n,k) * zfactor
          zaeq2_rn  (n,k) = zaeq2_rn  (n,k) * zfactor
          zaeq3_rn  (n,k) = zaeq3_rn  (n,k) * zfactor
          zaeq4_rn  (n,k) = zaeq4_rn  (n,k) * zfactor
          zaeq5_rn  (n,k) = zaeq5_rn  (n,k) * zfactor
        ENDDO ! n
      ENDDO   ! k

      ! Treatment at the eastern boundary
      IF (iendpar > iend) THEN
        i   = iendpar
        izz = izz+1
        IF (nradcoarse_y == 2) THEN
          DO k = 1, ke
            zti_rn    (izz,k) = zti_rn    (izz,k) + zti    (i,js-1,k) + zti    (i,js,k)
            zdpr_rn   (izz,k) = zdpr_rn   (izz,k) + dp0    (i,js-1,k) + dp0    (i,js,k)
            zclc_rn   (izz,k) = zclc_rn   (izz,k) + zclc   (i,js-1,k) + zclc   (i,js,k)
            zwv_rn    (izz,k) = zwv_rn    (izz,k) + zwv    (i,js-1,k) + zwv    (i,js,k)
            zsw_rn    (izz,k) = zsw_rn    (izz,k) + zsw    (i,js-1,k) + zsw    (i,js,k)
            zclwc_rn  (izz,k) = zclwc_rn  (izz,k) + zclwc  (i,js-1,k) + zclwc  (i,js,k)
            zciwc_rn  (izz,k) = zciwc_rn  (izz,k) + zciwc  (i,js-1,k) + zciwc  (i,js,k)
            zduco2f_rn(izz,k) = zduco2f_rn(izz,k) + zduco2f(i,js-1,k) + zduco2f(i,js,k)
            zduo3f_rn (izz,k) = zduo3f_rn (izz,k) + zduo3f (i,js-1,k) + zduo3f (i,js,k)
            zaeq1_rn  (izz,k) = zaeq1_rn  (izz,k) + zaeq1  (i,js-1,k) + zaeq1  (i,js,k)
            zaeq2_rn  (izz,k) = zaeq2_rn  (izz,k) + zaeq2  (i,js-1,k) + zaeq2  (i,js,k)
            zaeq3_rn  (izz,k) = zaeq3_rn  (izz,k) + zaeq3  (i,js-1,k) + zaeq3  (i,js,k)
            zaeq4_rn  (izz,k) = zaeq4_rn  (izz,k) + zaeq4  (i,js-1,k) + zaeq4  (i,js,k)
            zaeq5_rn  (izz,k) = zaeq5_rn  (izz,k) + zaeq5  (i,js-1,k) + zaeq5  (i,js,k)
          ENDDO ! k
        ELSE  ! nradcoarse_y == 1
          DO k = 1, ke
            zti_rn    (izz,k) = zti_rn    (izz,k) + zti    (i,js,k)
            zdpr_rn   (izz,k) = zdpr_rn   (izz,k) + dp0    (i,js,k)
            zclc_rn   (izz,k) = zclc_rn   (izz,k) + zclc   (i,js,k)
            zwv_rn    (izz,k) = zwv_rn    (izz,k) + zwv    (i,js,k)
            zsw_rn    (izz,k) = zsw_rn    (izz,k) + zsw    (i,js,k)
            zclwc_rn  (izz,k) = zclwc_rn  (izz,k) + zclwc  (i,js,k)
            zciwc_rn  (izz,k) = zciwc_rn  (izz,k) + zciwc  (i,js,k)
            zduco2f_rn(izz,k) = zduco2f_rn(izz,k) + zduco2f(i,js,k)
            zduo3f_rn (izz,k) = zduo3f_rn (izz,k) + zduo3f (i,js,k)
            zaeq1_rn  (izz,k) = zaeq1_rn  (izz,k) + zaeq1  (i,js,k)
            zaeq2_rn  (izz,k) = zaeq2_rn  (izz,k) + zaeq2  (i,js,k)
            zaeq3_rn  (izz,k) = zaeq3_rn  (izz,k) + zaeq3  (i,js,k)
            zaeq4_rn  (izz,k) = zaeq4_rn  (izz,k) + zaeq4  (i,js,k)
            zaeq5_rn  (izz,k) = zaeq5_rn  (izz,k) + zaeq5  (i,js,k)
          ENDDO ! k
        ENDIF
        DO k = 1, ke
          zti_rn    (izz,k) = zti_rn    (izz,k) * zfactor_b
          zdpr_rn   (izz,k) = zdpr_rn   (izz,k) * zfactor_b
          zclc_rn   (izz,k) = zclc_rn   (izz,k) * zfactor_b
          zwv_rn    (izz,k) = zwv_rn    (izz,k) * zfactor_b
          zsw_rn    (izz,k) = zsw_rn    (izz,k) * zfactor_b
          zclwc_rn  (izz,k) = zclwc_rn  (izz,k) * zfactor_b
          zciwc_rn  (izz,k) = zciwc_rn  (izz,k) * zfactor_b
          zduco2f_rn(izz,k) = zduco2f_rn(izz,k) * zfactor_b
          zduo3f_rn (izz,k) = zduo3f_rn (izz,k) * zfactor_b
          zaeq1_rn  (izz,k) = zaeq1_rn  (izz,k) * zfactor_b
          zaeq2_rn  (izz,k) = zaeq2_rn  (izz,k) * zfactor_b
          zaeq3_rn  (izz,k) = zaeq3_rn  (izz,k) * zfactor_b
          zaeq4_rn  (izz,k) = zaeq4_rn  (izz,k) * zfactor_b
          zaeq5_rn  (izz,k) = zaeq5_rn  (izz,k) * zfactor_b
        ENDDO ! k
      ENDIF ! Rand_x

      ! Surface
      izz = (iendrad-istartrad+nradcoarse) / nradcoarse

      IF (nradcoarse_y == 2) THEN
        DO n = 1, izz
          i = istartrad + (n-1) * nradcoarse
          zti_rn  (n,ke+1) = zti_rn  (n,ke+1) + zti  (i,js-1,ke+1) + zti  (i+1,js-1,ke+1) + zti  (i,js,ke+1) + zti  (i+1,js,ke+1)
          zalso_rn(n)      = zalso_rn(n)      + zalso(i,js-1)      + zalso(i+1,js-1)      + zalso(i,js)      + zalso(i+1,js)
#ifdef COUP_OAS_COS
          zpalp_rn(n)      = zpalp_rn(n)      + zpalp(i,js-1)      + zpalp(i+1,js-1)      + zpalp(i,js)      + zpalp(i+1,js)
#endif
          zalth_rn(n)      = zalth_rn(n)      + zalth(i,js)        + zalth(i+1,js-1)      + zalth(i,js)      + zalth(i+1,js)
          zapre_rn(n)      = zapre_rn(n)      + p0hl (i,js,ke+1)   + p0hl (i+1,js-1,ke+1) + p0hl (i,js,ke+1) + p0hl (i+1,js,ke+1)
          zsmu0_rn(n)      = zsmu0_rn(n)      + zsmu0(i,js)        + zsmu0(i+1,js-1)      + zsmu0(i,js)      + zsmu0(i+1,js)
        ENDDO ! n
      ELSE  ! nradcoarse_y == 1
        DO n = 1, izz
          i = istartrad + (n-1) * nradcoarse
          zti_rn  (n,ke+1) = zti_rn  (n,ke+1) + zti  (i,js,ke+1) + zti  (i+1,js,ke+1)
          zalso_rn(n)      = zalso_rn(n)      + zalso(i,js)      + zalso(i+1,js)
#ifdef COUP_OAS_COS
          zpalp_rn(n)      = zpalp_rn(n)      + zpalp(i,js)      + zpalp(i+1,js)
#endif
          zalth_rn(n)      = zalth_rn(n)      + zalth(i,js)      + zalth(i+1,js)
          zapre_rn(n)      = zapre_rn(n)      + p0hl (i,js,ke+1) + p0hl (i+1,js,ke+1)
          zsmu0_rn(n)      = zsmu0_rn(n)      + zsmu0(i,js)      + zsmu0(i+1,js)
        ENDDO ! n
      ENDIF ! js < jendpar
      DO n = 1, izz
        zti_rn  (n,ke+1) = zti_rn  (n,ke+1) * zfactor
        zalso_rn(n)      = zalso_rn(n)      * zfactor
#ifdef COUP_OAS_COS
        zpalp_rn(n)      = zpalp_rn(n)      * zfactor
#endif
        zalth_rn(n)      = zalth_rn(n)      * zfactor
        zapre_rn(n)      = zapre_rn(n)      * zfactor
        zsmu0_rn(n)      = zsmu0_rn(n)      * zfactor
      ENDDO ! n

      ! Treatment at the eastern boundary
      IF (iendpar > iend) THEN
        izz = izz+1
        i   = iendpar
        IF (nradcoarse_y == 2) THEN
          zti_rn  (izz,ke+1) = zti_rn  (izz,ke+1) + zti  (i,js-1,ke+1) + zti  (i,js,ke+1)
          zalso_rn(izz)      = zalso_rn(izz)      + zalso(i,js-1)      + zalso(i,js)
#ifdef COUP_OAS_COS
          zpalp_rn(izz)      = zpalp_rn(izz)      + zpalp(i,js-1)      + zpalp(i,js)
#endif
          zalth_rn(izz)      = zalth_rn(izz)      + zalth(i,js-1)      + zalth(i,js)
          zapre_rn(izz)      = zapre_rn(izz)      + p0hl (i,js-1,ke+1) + p0hl (i,js,ke+1)
          zsmu0_rn(izz)      = zsmu0_rn(izz)      + zsmu0(i,js-1)      + zsmu0(i,js)
        ELSE  ! nradcoarse_y == 1
          zti_rn  (izz,ke+1) = zti_rn  (izz,ke+1) + zti  (i,js,ke+1)
          zalso_rn(izz)      = zalso_rn(izz)      + zalso(i,js)
#ifdef COUP_OAS_COS
          zpalp_rn(izz)      = zpalp_rn(izz)      + zpalp(i,js)
#endif
          zalth_rn(izz)      = zalth_rn(izz)      + zalth(i,js)
          zapre_rn(izz)      = zapre_rn(izz)      + p0hl (i,js,ke+1)
          zsmu0_rn(izz)      = zsmu0_rn(izz)      + zsmu0(i,js)
        ENDIF !
        zti_rn  (izz,ke+1) = zti_rn  (izz,ke+1) * zfactor_b
        zalso_rn(izz)      = zalso_rn(izz)      * zfactor_b
#ifdef COUP_OAS_COS
        zpalp_rn(izz)      = zpalp_rn(izz)      * zfactor_b
#endif
        zalth_rn(izz)      = zalth_rn(izz)      * zfactor_b
        zapre_rn(izz)      = zapre_rn(izz)      * zfactor_b
        zsmu0_rn(izz)      = zsmu0_rn(izz)      * zfactor_b
      ENDIF !Rand_x

      ! Setting of first-dimension array boundary for routine fesft
      ki1ec=izz

      ! set default value for sykview; lradtopo must not be chosen with nradcoarse > 1
      zskyv_rn(:) = 1.0_ireals

#ifdef COUP_OAS_COS  !MU: change to 4.11 (zskyv_rn instead of skyview)
      CALL fesft                                                            &
          (zti_rn,    zdpr_rn,     zclc_rn,     zwv_rn,    zsw_rn,          &
           zclwc_rn,  zciwc_rn,    zduco2f_rn,  zduo3f_rn,                  &
           zaeq1_rn,  zaeq2_rn,    zaeq3_rn,    zaeq4_rn,  zaeq5_rn,        &
           zapre_rn,  zsmu0_rn,    zalso_rn,    zalth_rn,  zskyv_rn,        &
           swdir_cor, zstb,        zsct,        zpalp_rn,                   &
           ki1sd,     ki1ed,       ki2sd,       ki2ed,     ki3sd,    ki3ed, &
           ki1sc,     ki1ec,       ki2sc,       ki2ec,     ki3sc,    ki3ec, &
           lsolar(js),lcrf,        .FALSE.,     izdebug,   js,              &
           zflt_rn,   zfls_rn,     zflt_s_rn,   zfls_s_rn, zflsdir_rn,      &
           zfltd_rn,  zfltu_rn,    zflsd_rn,    zflsu_rn,  zflsp_rn,        &
           zflpar_rn, zflsu_par_rn, zflsd_par_rn, zflsp_par_rn)
#else
      CALL fesft                                                            &
          (zti_rn,    zdpr_rn,     zclc_rn,     zwv_rn,    zsw_rn,          &
           zclwc_rn,  zciwc_rn,    zduco2f_rn,  zduo3f_rn,                  &
           zaeq1_rn,  zaeq2_rn,    zaeq3_rn,    zaeq4_rn,  zaeq5_rn,        &
           zapre_rn,  zsmu0_rn,    zalso_rn,    zalth_rn,  zskyv_rn,        &
           swdir_cor, zstb,        zsct,                                    &
           ki1sd,     ki1ed,       ki2sd,       ki2ed,     ki3sd,    ki3ed, &
           ki1sc,     ki1ec,       ki2sc,       ki2ec,     ki3sc,    ki3ec, &
           lsolar(js),lcrf,        .FALSE.,     izdebug,   js,              &
           zflt_rn,   zfls_rn,     zflt_s_rn,   zfls_s_rn, zflsdir_rn,      &
           zfltd_rn,  zfltu_rn,    zflsd_rn,    zflsu_rn,  zflsp_rn,        &
           zflpar_rn, zflsu_par_rn, zflsd_par_rn, zflsp_par_rn)
#endif

      ! Store back results from fesft
      DO k=1,ke
        DO ii=0,nradcoarse-1
          izz=0
          DO i=istartrad,iendrad,nradcoarse
            izz=izz+1
            zflt(i+ii,k) = zflt_rn(izz,k)
            zfls(i+ii,k) = zfls_rn(izz,k)
!           zdpr(i+ii,k) = zdpr_rn(izz,k)
          ENDDO !i
        ENDDO !ii
      ENDDO !k
      IF (iendpar > iend) THEN
!CDIR NOVECTOR
        DO i=iendrad+1,iendpar,nradcoarse
          izz=izz+1
!CDIR NOVECTOR
          DO ii=0,(iendpar-iendrad)-1
            DO k=1,ke
              zflt(i+ii,k) = zflt_rn(izz,k)
              zfls(i+ii,k) = zfls_rn(izz,k)
!             zdpr(i+ii,k) = zdpr_rn(izz,k)
            ENDDO !k
          ENDDO !i
        ENDDO !ii
      ENDIF

      ! for the Climate-LM Version: solar direct radiation
      DO k=1,ke1
        DO ii=0,nradcoarse-1
          izz=0
          DO i=istartrad,iendrad,nradcoarse
            izz=izz+1
            zflsdir(i+ii,k) = zflsdir_rn(izz,k)
          ENDDO !i
        ENDDO !ii
      ENDDO !k
      IF (iendpar > iend) THEN
!CDIR NOVECTOR
        DO i=iendrad+1,iendpar,nradcoarse
          izz=izz+1
!CDIR NOVECTOR
          DO ii=0,(iendpar-iendrad)-1
            DO k=1,ke1
              zflsdir(i+ii,k) = zflsdir_rn(izz,k)
            ENDDO !k
          ENDDO !ii
        ENDDO !i
      ENDIF

      DO ii=0,nradcoarse-1
        izz=0
        DO i=istartrad,iendrad,nradcoarse
          izz=izz+1
          zflt     (i+ii,ke1) = zflt_rn     (izz,ke1)
          zfls     (i+ii,ke1) = zfls_rn     (izz,ke1)
          zflt_s   (i+ii)     = zflt_s_rn   (izz)
          zfls_s   (i+ii)     = zfls_s_rn   (izz)
          zalso    (i+ii,js)  = zalso_rn    (izz)
#ifdef COUP_OAS_COS
          zpalp    (i+ii,js)  = zpalp_rn    (izz)
#endif
          tg_rn    (i+ii)     = zti_rn      (izz,ke1)
          zfltu    (i+ii)     = zfltu_rn    (izz)
          zfltd    (i+ii)     = zfltd_rn    (izz)
          zflsu    (i+ii)     = zflsu_rn    (izz)
          zflsd    (i+ii)     = zflsd_rn    (izz)
          zflsp    (i+ii)     = zflsp_rn    (izz)
          zflpar   (i+ii)     = zflpar_rn   (izz)
          zflsu_par(i+ii)     = zflsu_par_rn(izz)
          zflsd_par(i+ii)     = zflsd_par_rn(izz)
          zflsp_par(i+ii)     = zflsp_par_rn(izz)
        ENDDO !i
      ENDDO

      IF (iendpar > iend) THEN
        DO i=iendrad+1,iendpar,nradcoarse
          izz=izz+1
          DO ii=0,(iendpar-iendrad)-1
            zflt     (i+ii,ke1) = zflt_rn     (izz,ke1)
            zfls     (i+ii,ke1) = zfls_rn     (izz,ke1)
            zflt_s   (i+ii)     = zflt_s_rn   (izz)
            zfls_s   (i+ii)     = zfls_s_rn   (izz)
            zalso    (i+ii,js)  = zalso_rn    (izz)
#ifdef COUP_OAS_COS
            zpalp    (i+ii,js)  = zpalp_rn    (izz)
#endif
            tg_rn    (i+ii)     = zti_rn      (izz,ke1)
            zfltu    (i+ii)     = zfltu_rn    (izz)
            zfltd    (i+ii)     = zfltd_rn    (izz)
            zflsu    (i+ii)     = zflsu_rn    (izz)
            zflsd    (i+ii)     = zflsd_rn    (izz)
            zflsp    (i+ii)     = zflsp_rn    (izz)
            zflpar   (i+ii)     = zflpar_rn   (izz)
            zflsu_par(i+ii)     = zflsu_par_rn(izz)
            zflsd_par(i+ii)     = zflsd_par_rn(izz)
            zflsp_par(i+ii)     = zflsp_par_rn(izz)
          ENDDO !ii
        ENDDO !i
      ENDIF

      !------------------------------------------------------------------------
      !
      ! Section 6:  Heating rates and radiation budget at surface
      !
      !------------------------------------------------------------------------

      DO  k = 1, ke
        DO jz1 = js-nradcoarse_y+1, js
          DO  i = istartrad, iendparrad
!           zfac = g/(cp_d*dp0 (i,js,k))
            zfac = g/(cp_d*dp0 (i,jz1,k))
!US         zfac = g/(cp_d*zdpr(i,k))
            sohr(i,jz1,k) = 0.0_ireals
!           IF (zsmu0(i,js) > zepemu) THEN
            IF (zsmu0(i,jz1) > zepemu) THEN
              sohr(i,jz1,k) = zfac * (zfls(i,k)-zfls(i,k+1))
            ENDIF
            thhr(i,jz1,k)   = zfac * (zflt(i,k)-zflt(i,k+1))
          ENDDO
        ENDDO
      ENDDO

      ! for the Climate-LM Version: solar direct radiation
      DO jz1 = js-nradcoarse_y+1, js
        DO  i = istartrad, iendparrad
          sodwddm(i,jz1) = 0.0_ireals
          IF (zsmu0(i,jz1) > zepemu) THEN
            sodwddm(i,jz1) = zflsdir(i,ke1) / zsmu0(i,jz1)
          ENDIF
        ENDDO

        DO  i = istartrad, iendparrad
          sobs(i,jz1) = 0.0_ireals
          pabs(i,jz1) = 0.0_ireals
          sobt(i,jz1) = 0.0_ireals

          swdir_s (i,jz1) = zflsp(i)
          swdifd_s(i,jz1) = zflsd(i)
          swdifu_s(i,jz1) = zflsu(i)
          lwd_s   (i,jz1) = zfltd(i)
          lwu_s   (i,jz1) = zfltu(i)

          zzflsp_par  (i,jz1) = zflsp_par(i)
          zzflsd_par  (i,jz1) = zflsd_par(i)
          zzflsu_par  (i,jz1) = zflsu_par(i)

!         IF (zsmu0(i,js) > zepemu) THEN
          IF (zsmu0(i,jz1) > zepemu) THEN
            sobs(i,jz1) = zfls_s(i)
            pabs(i,jz1) = zflpar(i)
            sobt(i,jz1) = zfls  (i,  1)
          ENDIF
          thbs    (i,jz1) = zflt_s(i)
          tg_ra   (i,jz1) = tg_rn (i)
          thbt    (i,jz1) = zflt  (i,  1)
#ifdef COUP_OAS_COS
          alb_rad (i,jz1,1) = zpalp (i,js)
          alb_rad (i,jz1,2) = zalso (i,js)
#else
          alb_rad (i,jz1) = zalso (i,js)
#endif
        ENDDO
      ENDDO

      zti_rn(:,:)     = 0.0_ireals
      zdpr_rn(:,:)    = 0.0_ireals
      zclc_rn(:,:)    = 0.0_ireals
      zwv_rn(:,:)     = 0.0_ireals
      zsw_rn(:,:)     = 0.0_ireals
      zclwc_rn(:,:)   = 0.0_ireals
      zciwc_rn(:,:)   = 0.0_ireals
      zduco2f_rn(:,:) = 0.0_ireals
      zduo3f_rn(:,:)  = 0.0_ireals
      zaeq1_rn(:,:)   = 0.0_ireals
      zaeq2_rn(:,:)   = 0.0_ireals
      zaeq3_rn(:,:)   = 0.0_ireals
      zaeq4_rn(:,:)   = 0.0_ireals
      zaeq5_rn(:,:)   = 0.0_ireals
      zalso_rn(:)     = 0.0_ireals
#ifdef COUP_OAS_COS
      zpalp_rn(:)     = 0.0_ireals
#endif
      zalth_rn(:)     = 0.0_ireals
      zapre_rn(:)     = 0.0_ireals
      zsmu0_rn(:)     = 0.0_ireals

      zfls_s_rn(:)    = 0.0_ireals
      zflt_s_rn(:)    = 0.0_ireals
      zflsp_rn (:)    = 0.0_ireals
      zflsd_rn (:)    = 0.0_ireals
      zflsu_rn (:)    = 0.0_ireals
      zfltd_rn (:)    = 0.0_ireals
      zfltu_rn (:)    = 0.0_ireals

    ENDIF !lrady

  ELSE ! .NOT. lradave:

!------------------------------------------------------------------------------
! Section 8:  Calculation of radiation fluxes in routine fesft
!------------------------------------------------------------------------------

#ifdef COUP_OAS_COS     !MU: change to 4.11 (zskyview(:,js) instead of skyview)
    CALL fesft                                                                             &
          (zti  (:,js,:),   dp0  (:,js,:), zclc   (:,js,:), zwv   (:,js,:), zsw  (:,js,:), &
           zclwc(:,js,:),   zciwc(:,js,:), zduco2f(:,js,:), zduo3f(:,js,:),                &
           zaeq1(:,js,:),   zaeq2(:,js,:), zaeq3  (:,js,:), zaeq4 (:,js,:), zaeq5(:,js,:), &
           zapre(:,js),     zsmu0(:,js),   zalso  (:,js),   zalth (:,js), zskyview(:,js),  &
           swdir_cor(:,js), zstb,          zsct,   zpalp  (:,js),                          &
           ki1sd,      ki1ed,       ki2sd,       ki2ed,     ki3sd,    ki3ed, &
           ki1sc,      ki1ec,       ki2sc,       ki2ec,     ki3sc,    ki3ec, &
           lsolar(js), lcrf,        lradtopo,    izdebug,   js,              &
           zflt,       zfls,        zflt_s,      zfls_s,    zflsdir,         &
           zfltd,      zfltu,       zflsd,       zflsu,     zflsp,           &
           zflpar,     zflsu_par,   zflsd_par,   zflsp_par)
#else
    CALL fesft                                                                             &
          (zti  (:,js,:),   dp0  (:,js,:), zclc   (:,js,:), zwv   (:,js,:), zsw  (:,js,:), &
           zclwc(:,js,:),   zciwc(:,js,:), zduco2f(:,js,:), zduo3f(:,js,:),                &
           zaeq1(:,js,:),   zaeq2(:,js,:), zaeq3  (:,js,:), zaeq4 (:,js,:), zaeq5(:,js,:), &
           zapre(:,js),     zsmu0(:,js),   zalso  (:,js),   zalth (:,js), zskyview(:,js),  &
           swdir_cor(:,js), zstb,          zsct,                                           &
           ki1sd,      ki1ed,       ki2sd,       ki2ed,     ki3sd,    ki3ed, &
           ki1sc,      ki1ec,       ki2sc,       ki2ec,     ki3sc,    ki3ec, &
           lsolar(js), lcrf,        lradtopo,    izdebug,   js,              &
           zflt,       zfls,        zflt_s,      zfls_s,    zflsdir,         &
           zfltd,      zfltu,       zflsd,       zflsu,     zflsp,           &
           zflpar,     zflsu_par,   zflsd_par,   zflsp_par)
#endif

!------------------------------------------------------------------------------
! Section 9:  Heating rates and radiation budget at surface
!------------------------------------------------------------------------------

  DO  k = 1, ke 
    DO  i = istartpar, iendpar
      zfac = g/(cp_d*dp0 (i,js,k))
      sohr(i,js,k) = 0.0
      IF (zsmu0(i,js) > zepemu) THEN
        sohr(i,js,k) = zfac * (zfls(i,k)-zfls(i,k+1))
      ENDIF
      thhr(i,js,k)   = zfac * (zflt(i,k)-zflt(i,k+1))
    ENDDO     
  ENDDO         
 
  DO  i = istartpar, iendpar
    sobs    (i,js) = 0.0
    pabs    (i,js) = 0.0
    sobt    (i,js) = 0.0
    swdir_s (i,js) = 0.0
    swdifd_s(i,js) = 0.0
    swdifu_s(i,js) = 0.0
    lwd_s   (i,js) = 0.0
    lwu_s   (i,js) = 0.0
    IF (zsmu0(i,js) > zepemu) THEN
      sobs  (i,js) = zfls_s(i)
      pabs  (i,js) = zflpar(i)
      sobt  (i,js) = zfls  (i,  1)
    END IF
    thbs    (i,js) = zflt_s(i)
    thbt    (i,js) = zflt  (i,  1)
#ifdef COUP_OAS_COS
    alb_rad (i,js,1) = zpalp (i,js)
    alb_rad (i,js,2) = zalso (i,js)
#else
    alb_rad (i,js) = zalso (i,js)
#endif
    swdir_s (i,js) = zflsp(i)
    swdifd_s(i,js) = zflsd(i)
    swdifu_s(i,js) = zflsu(i)
    lwd_s   (i,js) = zfltd(i)
    lwu_s   (i,js) = zfltu(i)

    zzflsp_par(i,js) = zflsp_par(i)
    zzflsd_par(i,js) = zflsd_par(i)
    zzflsu_par(i,js) = zflsu_par(i)
  ENDDO         
 
  ! for the Climate-LM Version: solar direct parallel radiation at the surface
  DO i = istartpar, iendpar
    sodwddm(i,js) = 0.0_ireals
    IF (zsmu0(i,js) > zepemu) THEN
      sodwddm(i,js) = zflsdir(i,ke1) / zsmu0(i,js)
    ENDIF
  ENDDO

  ENDIF !lradave

!-------------------------------------------------------------------------------
! end of loop over the model domain from south to north
!-------------------------------------------------------------------------------

ENDDO loop_south_north

  IF (lradave) THEN

    DEALLOCATE ( zti_rn       , STAT=izstatd )
    DEALLOCATE ( zdpr_rn      , STAT=izstatd )
    DEALLOCATE ( zclc_rn      , STAT=izstatd )
    DEALLOCATE ( zwv_rn       , STAT=izstatd )
    DEALLOCATE ( zsw_rn       , STAT=izstatd )
    DEALLOCATE ( zclwc_rn     , STAT=izstatd )
    DEALLOCATE ( zciwc_rn     , STAT=izstatd )
    DEALLOCATE ( zduco2f_rn   , STAT=izstatd )
    DEALLOCATE ( zduo3f_rn    , STAT=izstatd )
    DEALLOCATE ( zaeq1_rn     , STAT=izstatd )
    DEALLOCATE ( zaeq2_rn     , STAT=izstatd )
    DEALLOCATE ( zaeq3_rn     , STAT=izstatd )
    DEALLOCATE ( zaeq4_rn     , STAT=izstatd )
    DEALLOCATE ( zaeq5_rn     , STAT=izstatd )
    DEALLOCATE ( zapre_rn     , STAT=izstatd )
    DEALLOCATE ( zsmu0_rn     , STAT=izstatd )
    DEALLOCATE ( zalth_rn     , STAT=izstatd )
    DEALLOCATE ( zalso_rn     , STAT=izstatd )
#ifdef COUP_OAS_COS
    DEALLOCATE ( zpalp_rn     , STAT=izstatd )
#endif
    DEALLOCATE ( zflt_rn      , STAT=izstatd )
    DEALLOCATE ( zfls_rn      , STAT=izstatd )
    DEALLOCATE ( zflpar_rn    , STAT=izstatd )
    DEALLOCATE ( zflsu_par_rn , STAT=izstatd )
    DEALLOCATE ( zflsd_par_rn , STAT=izstatd )
    DEALLOCATE ( zflsp_par_rn , STAT=izstatd )
    DEALLOCATE ( zflsdir_rn   , STAT=izstatd )
    DEALLOCATE ( tg_rn        , STAT=izstatd )
    DEALLOCATE ( zfls_s_rn    , STAT=izstatd )
    DEALLOCATE ( zflt_s_rn    , STAT=izstatd )
    DEALLOCATE ( zflsp_rn     , STAT=izstatd )
    DEALLOCATE ( zflsd_rn     , STAT=izstatd )
    DEALLOCATE ( zflsu_rn     , STAT=izstatd )
    DEALLOCATE ( zfltd_rn     , STAT=izstatd )
    DEALLOCATE ( zfltu_rn     , STAT=izstatd )

    DO js=jstartpar,jendpar
      DO i = istartpar,iendpar
        thbs (i,js)  = thbs(i,js) + zstb*(1._ireals - ctalb)*(tg_ra(i,js)**4)

        ! this was eliminated in testsuite 3.4
        ! but keep it for the moment to reproduce results with Version 3.22
#ifdef COUP_OAS_COS
        zalbfak      = 1._ireals/(1._ireals-alb_rad(i,js,2))
#else
        zalbfak      = 1._ireals/(1._ireals-alb_rad(i,js))
#endif
        sobs (i,js)  = sobs(i,js) * zalbfak
        pabs (i,js)  = pabs(i,js) * zalbfak
      ENDDO
    ENDDO

    DEALLOCATE ( tg_ra, STAT=izstatd )

    IF (lradf_avg) THEN

!!$      IF (num_compute > 1) THEN
        kzdims(1:24)=(/ke,ke,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0/)
#ifdef COUP_OAS_COS        !MU: some changes compared to 4.11
        CALL exchg_boundaries                                              &
          (9, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute, ie, je,         &
           kzdims, jstartpar, jendpar, 1, nboundlines, my_cart_neigh,      &
           lperi_x, lperi_y, l2dim, &
           20000+ntstep, .FALSE.,    ncomm_type, izerror, yzerrmsg,        &
           thhr(:,:,:),sohr(:,:,:),thbs(:,:),thbt(:,:),                    &
           sobs(:,:),sobt(:,:),pabs(:,:),alb_rad(:,:,:), sodwddm(:,:),       &
           swdir_s(:,:), swdifd_s(:,:), swdifu_s(:,:), lwd_s(:,:),           &
           lwu_s(:,:), zzflsp_par(:,:),zzflsd_par(:,:),zzflsu_par(:,:) )
#else
        CALL exchg_boundaries                                              &
          (9, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute, ie, je,         &
           kzdims, jstartpar, jendpar, 1, nboundlines, my_cart_neigh,      &
           lperi_x, lperi_y, l2dim, &
           20000+ntstep, .FALSE.,    ncomm_type, izerror, yzerrmsg,        &
           thhr(:,:,:),sohr(:,:,:),thbs(:,:),thbt(:,:),                    &
           sobs(:,:),sobt(:,:),pabs(:,:),alb_rad(:,:), sodwddm(:,:),       &
           swdir_s(:,:), swdifd_s(:,:), swdifu_s(:,:), lwd_s(:,:),         &
           lwu_s(:,:), zzflsp_par(:,:),zzflsd_par(:,:),zzflsu_par(:,:) )
#endif
!!$      ENDIF

      ALLOCATE ( zsohr       (ie,je,ke), STAT=izstata )
      ALLOCATE ( zthhr       (ie,je,ke), STAT=izstata )
      ALLOCATE ( zsobs       (ie,je)   , STAT=izstata )
      ALLOCATE ( zsobt       (ie,je)   , STAT=izstata )
      ALLOCATE ( zthbs       (ie,je)   , STAT=izstata )
      ALLOCATE ( zthbt       (ie,je)   , STAT=izstata )
      ALLOCATE ( zpabs       (ie,je)   , STAT=izstata )
      ALLOCATE ( zsodwddm    (ie,je)   , STAT=izstata )
      ALLOCATE ( z_zzflsp    (ie,je)   , STAT=izstata )
      ALLOCATE ( z_zzflsd    (ie,je)   , STAT=izstata )
      ALLOCATE ( z_zzflsu    (ie,je)   , STAT=izstata )
      ALLOCATE ( z_zzflsp_par(ie,je)   , STAT=izstata )
      ALLOCATE ( z_zzflsd_par(ie,je)   , STAT=izstata )
      ALLOCATE ( z_zzflsu_par(ie,je)   , STAT=izstata )
      ALLOCATE ( z_zzfltd    (ie,je)   , STAT=izstata )
      ALLOCATE ( z_zzfltu    (ie,je)   , STAT=izstata ) 

      zthhr        (:,:,:) = thhr       (:,:,:)
      zsohr        (:,:,:) = sohr       (:,:,:)
      zthbs        (:,:)   = thbs       (:,:)
      zthbt        (:,:)   = thbt       (:,:)
      zsobs        (:,:)   = sobs       (:,:)
      zsobt        (:,:)   = sobt       (:,:)
      zpabs        (:,:)   = pabs       (:,:)
      zsodwddm     (:,:)   = sodwddm    (:,:)
      z_zzflsp     (:,:)   = swdir_s    (:,:)
      z_zzflsd     (:,:)   = swdifd_s   (:,:)
      z_zzflsu     (:,:)   = swdifu_s   (:,:)
      z_zzfltd     (:,:)   = lwd_s      (:,:)
      z_zzfltu     (:,:)   = lwu_s      (:,:)
      z_zzflsp_par (:,:)   = zzflsp_par (:,:)
      z_zzflsd_par (:,:)   = zzflsd_par (:,:)
      z_zzflsu_par (:,:)   = zzflsu_par (:,:)

      DO k=1,ke
        DO js=jstart,jend
          DO i = istart,iend
            thhr(i,js,k)  =  ( zcent*zthhr(i,js,k)                  &
                                + zside*( zthhr(i-1,js,k  ) + zthhr(i+1,js,k  )    &
                                        + zthhr(i,js-1,k  ) + zthhr(i,js+1,k  )  ) &
                                + zedge*( zthhr(i-1,js-1,k) + zthhr(i+1,js-1,k)    &
                                +         zthhr(i-1,js+1,k) + zthhr(i+1,js+1,k)  ) )
            sohr(i,js,k)  =  ( zcent*zsohr(i,js,k)                  &
                                + zside*( zsohr(i-1,js,k  ) + zsohr(i+1,js,k  )    &
                                        + zsohr(i,js-1,k  ) + zsohr(i,js+1,k  )  ) &
                                + zedge*( zsohr(i-1,js-1,k) + zsohr(i+1,js-1,k)    &
                                        + zsohr(i-1,js+1,k) + zsohr(i+1,js+1,k)  ) )
          ENDDO
        ENDDO
      ENDDO !k
      DO k=1,ke1
        DO js=jstart,jend
          DO i = istart,iend
            sodwddm(i,js)  = ( zcent*zsodwddm(i,js)                                &
                               + zside*( zsodwddm(i-1,js  ) + zsodwddm(i+1,js  )   &
                                       + zsodwddm(i,js-1  ) + zsodwddm(i,js+1  ) ) &
                               + zedge*( zsodwddm(i-1,js-1) + zsodwddm(i+1,js-1)   &
                               +         zsodwddm(i-1,js+1) + zsodwddm(i+1,js+1) ) )
          ENDDO
        ENDDO
      ENDDO !k
      DO js=jstart,jend
        DO i = istart,iend
          thbs(i,js)  =  ( zcent*zthbs(i,js)                  &
                            + zside*( zthbs(i-1,js  ) + zthbs(i+1,js  )    &
                                    + zthbs(i,js-1  ) + zthbs(i,js+1  )  ) &
                            + zedge*( zthbs(i-1,js-1) + zthbs(i+1,js-1)    &
                                    + zthbs(i-1,js+1) + zthbs(i+1,js+1)  ) )
          thbt(i,js)  =  ( zcent*zthbt(i,js)                  &
                            + zside*( zthbt(i-1,js  ) + zthbt(i+1,js  )    &
                                    + zthbt(i,js-1  ) + zthbt(i,js+1  )  ) &
                            + zedge*( zthbt(i-1,js-1) + zthbt(i+1,js-1)    &
                                    + zthbt(i-1,js+1) + zthbt(i+1,js+1)  ) )
          sobs(i,js)  =  ( zcent*zsobs(i,js)                  &
                            + zside*( zsobs(i-1,js  ) + zsobs(i+1,js  )    &
                                    + zsobs(i,js-1  ) + zsobs(i,js+1  )  ) &
                            + zedge*( zsobs(i-1,js-1) + zsobs(i+1,js-1)    &
                                    + zsobs(i-1,js+1) + zsobs(i+1,js+1)  ) )
          sobt(i,js)  =  ( zcent*zsobt(i,js)                  &
                            + zside*( zsobt(i-1,js  ) + zsobt(i+1,js  )    &
                                    + zsobt(i,js-1  ) + zsobt(i,js+1  )  ) &
                            + zedge*( zsobt(i-1,js-1) + zsobt(i+1,js-1)    &
                                    + zsobt(i-1,js+1) + zsobt(i+1,js+1)  ) )
          pabs(i,js)  =  ( zcent*zpabs(i,js)                  &
                            + zside*( zpabs(i-1,js  ) + zpabs(i+1,js  )    &
                                    + zpabs(i,js-1  ) + zpabs(i,js+1  )  ) &
                            + zedge*( zpabs(i-1,js-1) + zpabs(i+1,js-1)    &
                                    + zpabs(i-1,js+1) + zpabs(i+1,js+1)  ) )
          swdir_s(i,js) =  ( zcent*z_zzflsp(i,js)                  &
                            + zside*( z_zzflsp(i-1,js  ) + z_zzflsp(i+1,js  )    &
                                    + z_zzflsp(i,js-1  ) + z_zzflsp(i,js+1  )  ) &
                            + zedge*( z_zzflsp(i-1,js-1) + z_zzflsp(i+1,js-1)    &
                                    + z_zzflsp(i-1,js+1) + z_zzflsp(i+1,js+1)  ) )
          swdifd_s(i,js)  =  ( zcent*z_zzflsd(i,js)                  &
                            + zside*( z_zzflsd(i-1,js  ) + z_zzflsd(i+1,js  )    &
                                    + z_zzflsd(i,js-1  ) + z_zzflsd(i,js+1  )  ) &
                            + zedge*( z_zzflsd(i-1,js-1) + z_zzflsd(i+1,js-1)    &
                                    + z_zzflsd(i-1,js+1) + z_zzflsd(i+1,js+1)  ) )
          swdifu_s(i,js)  =  ( zcent*z_zzflsu(i,js)                  &
                            + zside*( z_zzflsu(i-1,js  ) + z_zzflsu(i+1,js  )    &
                                    + z_zzflsu(i,js-1  ) + z_zzflsu(i,js+1  )  ) &
                            + zedge*( z_zzflsu(i-1,js-1) + z_zzflsu(i+1,js-1)    &
                                    + z_zzflsu(i-1,js+1) + z_zzflsu(i+1,js+1)  ) )
          lwd_s(i,js)  =  ( zcent*z_zzfltd(i,js)                  &
                            + zside*( z_zzfltd(i-1,js  ) + z_zzfltd(i+1,js  )    &
                                    + z_zzfltd(i,js-1  ) + z_zzfltd(i,js+1  )  ) &
                            + zedge*( z_zzfltd(i-1,js-1) + z_zzfltd(i+1,js-1)    &
                                    + z_zzfltd(i-1,js+1) + z_zzfltd(i+1,js+1)  ) )
          lwu_s(i,js)  =  ( zcent*z_zzfltu(i,js)                  &
                            + zside*( z_zzfltu(i-1,js  ) + z_zzfltu(i+1,js  )    &
                                    + z_zzfltu(i,js-1  ) + z_zzfltu(i,js+1  )  ) &
                            + zedge*( z_zzfltu(i-1,js-1) + z_zzfltu(i+1,js-1)    &
                                    + z_zzfltu(i-1,js+1) + z_zzfltu(i+1,js+1)  ) )
          zzflsp_par(i,js)  =  ( zcent*z_zzflsp_par(i,js)                  &
                            + zside*( z_zzflsp_par(i-1,js  ) + z_zzflsp_par(i+1,js  )    &
                                    + z_zzflsp_par(i,js-1  ) + z_zzflsp_par(i,js+1  )  ) &
                            + zedge*( z_zzflsp_par(i-1,js-1) + z_zzflsp_par(i+1,js-1)    &
                                    + z_zzflsp_par(i-1,js+1) + z_zzflsp_par(i+1,js+1)  ) )
          zzflsd_par(i,js)  =  ( zcent*z_zzflsd_par(i,js)                  &
                            + zside*( z_zzflsd_par(i-1,js  ) + z_zzflsd_par(i+1,js  )    &
                                    + z_zzflsd_par(i,js-1  ) + z_zzflsd_par(i,js+1  )  ) &
                            + zedge*( z_zzflsd_par(i-1,js-1) + z_zzflsd_par(i+1,js-1)    &
                                    + z_zzflsd_par(i-1,js+1) + z_zzflsd_par(i+1,js+1)  ) )
          zzflsu_par(i,js)  =  ( zcent*z_zzflsu_par(i,js)                  &
                            + zside*( z_zzflsu_par(i-1,js  ) + z_zzflsu_par(i+1,js  )    &
                                    + z_zzflsu_par(i,js-1  ) + z_zzflsu_par(i,js+1  )  ) &
                            + zedge*( z_zzflsu_par(i-1,js-1) + z_zzflsu_par(i+1,js-1)    &
                                    + z_zzflsu_par(i-1,js+1) + z_zzflsu_par(i+1,js+1)  ) )
        ENDDO !i
      ENDDO !js

      DEALLOCATE ( zsohr       , STAT=izstatd )
      DEALLOCATE ( zthhr       , STAT=izstatd )
      DEALLOCATE ( zsobs       , STAT=izstatd )
      DEALLOCATE ( zsobt       , STAT=izstatd )
      DEALLOCATE ( zthbs       , STAT=izstatd )
      DEALLOCATE ( zthbt       , STAT=izstatd )
      DEALLOCATE ( zpabs       , STAT=izstatd )
      DEALLOCATE ( zsodwddm    , STAT=izstatd )
      DEALLOCATE ( z_zzflsp    , STAT=izstatd )
      DEALLOCATE ( z_zzflsd    , STAT=izstatd )
      DEALLOCATE ( z_zzflsu    , STAT=izstatd )
      DEALLOCATE ( z_zzflsp_par, STAT=izstatd )
      DEALLOCATE ( z_zzflsd_par, STAT=izstatd )
      DEALLOCATE ( z_zzflsu_par, STAT=izstatd )
      DEALLOCATE ( z_zzfltd    , STAT=izstatd )
      DEALLOCATE ( z_zzfltu    , STAT=izstatd )

    ENDIF !lradf_avg

    DO js=jstartpar,jendpar
      DO i = istartpar,iendpar
        thbs (i,js) = thbs(i,js) - zstb*(1._ireals - ctalb)*(t_g(i,js,nzx)**4)

        ! such it has been tested in testsuite 3.4
        ! but keep it for the moment to reproduce results with Version 3.22
        zalbfak = (1.0_ireals-zalb_rad_ori(i,js))
        sobs (i,js) = sobs(i,js) * zalbfak
        pabs (i,js) = pabs(i,js) * zalbfak

        ! And this seems to be the better version
        ! swdifu_s(i,js) = ( swdir_s(i,js) + swdifd_s(i,js) ) * zalb_rad_ori(i,js)
        ! sobs(i,js)     =   swdir_s(i,js) + swdifd_s(i,js)   - swdifu_s(i,js)
        ! pabs(i,js)     = zzflsp_par(i,js)+zzflsd_par(i,js)  -  zzflsu_par(i,js)

      ENDDO !i
    ENDDO !js
#ifdef COUP_OAS_COS
    alb_rad(:,:,2) = zalb_rad_ori(:,:)
#else
    alb_rad(:,:) = zalb_rad_ori(:,:)
#endif

    DEALLOCATE ( zalb_rad_ori, STAT=izstatd )

    IF (lradtopo) THEN
      ! Storage of individual flux components for thermal radiative surface flux

      IF (.NOT. lemiss) THEN
        zemissivity =         1.0_ireals - ctalb ! surface emissivity
        zemissfac   = zstb * (1.0_ireals - ctalb)
      ENDIF

      DO js = jstartpar,jendpar
        DO i = istartpar,iendpar

          ! Recompute surface thermal flux components based on 
          ! lower boundary condition
          IF (lemiss) THEN
             lwd_s(i,js) = (thbs(i,js) + zstb*emis_rad(i,js)*t_g(i,js,nzx)**4) / emis_rad(i,js)
          ELSE
             lwd_s(i,js) = (thbs(i,js) + zemissfac*t_g(i,js,nzx)**4) / zemissivity
          ENDIF

          ! correction for thermal fluxes
          lwu_s(i,js) = lwd_s(i,js) - thbs(i,js)
          lwd_s(i,js) = lwd_s(i,js) *               skyview(i,js) +        &
                        lwu_s(i,js) * (1.0_ireals - skyview(i,js))

          thbs(i,js)  = lwd_s(i,js) - lwu_s(i,js)

          IF (zsmu0(i,js) > zepemu) THEN
            ! correction for solar fluxes
            ! direct down corrected
            swdir_s(i,js) = swdir_cor(i,js) * swdir_s(i,js)

            ! diffuse down corrected
            swdifd_s(i,js) = swdifd_s(i,js) *             skyview(i,js)    &
                           + swdifu_s(i,js) * (1.0_ireals-skyview(i,js))

            ! diffuse up adapted to new other components
#ifdef COUP_OAS_COS
            swdifu_s(i,js) = (swdir_s(i,js) + swdifd_s(i,js)) * alb_rad(i,js,2)
#else
            swdifu_s(i,js) = (swdir_s(i,js) + swdifd_s(i,js)) * alb_rad(i,js)
#endif
            sobs    (i,js) =  swdir_s(i,js) + swdifd_s(i,js)  - swdifu_s(i,js)

            ! correction for solar fluxes of photosynthetic active radiation
            IF ( (zzflsp_par(i,js) > 0.0_ireals) .OR.                      &
                 (zzflsd_par(i,js) > 0.0_ireals)             ) THEN
              zalbradtopo =  zzflsu_par(i,js) /                            &
                            (zzflsp_par(i,js)+zzflsd_par(i,js))

              ! direct down corrected
              zzflsp_par(i,js) = swdir_cor(i,js) * zzflsp_par(i,js)

              ! diffuse down corrected
              zzflsd_par(i,js) = zzflsd_par(i,js) *             skyview(i,js) &
                               + zzflsu_par(i,js) * (1.0_ireals-skyview(i,js))

              ! diffuse up adapted to new other components
              zzflsu_par(i,js) = (zzflsp_par(i,js) + zzflsd_par(i,js))      &
                                 * zalbradtopo
              pabs(i,js)       = zzflsp_par(i,js) + zzflsd_par(i,js)        &
                               - zzflsu_par(i,js)
            ELSE
              pabs(i,js)       = 0.0_ireals
            ENDIF
          ELSE
            sobs(i,js) = 0.0_ireals
            pabs(i,js) = 0.0_ireals
          ENDIF
        ENDDO
      ENDDO

    ENDIF !lradtopo

  ENDIF !lradave

  ! Set some additional output parameters
  DO js = jstartpar,jendpar
    DO i = istartpar,iendpar
      sod_t(i,js) = zsct * zsmu0(i,js)
    ENDDO
  ENDDO

#ifdef COSMOART
  IF(l_cosmo_art) THEN
    ! CK 20101204 Photolysis rates only need to be calculated if gases are ON
    If (lgas) CALL calcjval
  ENDIF
#endif

!-------------------------------------------------------------------------------
! End of the subroutine
!-------------------------------------------------------------------------------

END SUBROUTINE organize_radiation

!==============================================================================
!==============================================================================
!+ Module procedure in "Radiation" 
!------------------------------------------------------------------------------

#ifdef COUP_OAS_COS
!CPS we receive albedo for direct radiation also (palp)
SUBROUTINE fesft(                                                        &
       pti   , pdp  , pclc  , pwv  , psw  , pqlwc, pqiwc, pduco2, pduo3, &
       paeq1 , paeq2, paeq3 , paeq4, paeq5,                              &
       papre , psmu0, palso , palth, pskyview, pfcor,                    &
       psig  , psct , palp  ,                                            &
       ki1sd , ki1ed, ki2sd , ki2ed, ki3sd, ki3ed,                       &
       ki1sc , ki1ec, ki2sc , ki2ec, ki3sc, ki3ec,                       &
       lsolar, lcrf , lradtopo, idebug, jindex,                          &
       pflt  , pfls , pflt_s, pfls_s, pflsdir,                           &
       pfltd , pfltu, pflsd , pflsu, pflsp,                              &
       pflpar, pflsu_par, pflsd_par, pflsp_par)
#else
SUBROUTINE fesft(                                                        &
       pti   , pdp  , pclc  , pwv  , psw  , pqlwc, pqiwc, pduco2, pduo3, &
       paeq1 , paeq2, paeq3 , paeq4, paeq5,                              &
       papre , psmu0, palso , palth, pskyview, pfcor,                    &
       psig  , psct ,                                                    &
       ki1sd , ki1ed, ki2sd , ki2ed, ki3sd, ki3ed,                       &
       ki1sc , ki1ec, ki2sc , ki2ec, ki3sc, ki3ec,                       &
       lsolar, lcrf , lradtopo, idebug, jindex,                          &
       pflt  , pfls , pflt_s, pfls_s, pflsdir,                           &
       pfltd , pfltu, pflsd , pflsu, pflsp,                              &
       pflpar, pflsu_par, pflsd_par, pflsp_par)
#endif

!------------------------------------------------------------------------------
!
! Description:
!
!   The module procedure fesft organizes the radiative transfer calculations.
!
! Method:
!
! This routine organizes the radiative transfer calculations by
! calling a set of dedicated routines for the calculation of
! basic optical properties (opt_th/opt_so), the derivation of
! layer coefficients (coe_th/coe_so) for an implicit delta-two-
! stream scheme (cf.Ritter and Geleyn, 1992) and the inversion
! (inv_th/inv_so) of the corresponding system matrix. These 
! operations are performed seperately for thermal and solar parts
! of the spectrum and are embedded in loops over the various
! spectral intervals. Within each interval, a data-controlled
! decision is taken, whether the so-called ESFT or FESFT approach
! is used for the handling of gaseous absorption (cf. Ritter and
! Geleyn, 1992).
! Controlled by the logical input variable LCRF, the calculation
! of radiative fluxes in cloud-free conditions can be done in
! addition to the results for the given atmospheric cloud structure.
! (not implemented yet)
! Before the actual flux calculation starts, some preliminary steps
! provide utility arrays which are applicable to all spectral inter-
! vals (e.g. cloud geometry factors, integrated layer water content, etc.)
!
!------------------------------------------------------------------------------

! Subroutine arguments:
! --------------------

! Input data
! ----------
  INTEGER (KIND=iintegers), INTENT (IN) ::  &
    ! indices for array dimensioning
    ki1sd,       & ! start index for first  array dimension
    ki1ed,       & ! end   index for first  array dimension
    ki2sd,       & ! start index for second array dimension
    ki2ed,       & ! end   index for second array dimension
    ki3sd,       & ! start index for third  array dimension
    ki3ed,       & ! end   index for third  array dimension

    ! and the same for the computations
    ki1sc,       & ! start index for first  array computation
    ki1ec,       & ! end   index for first  array computation
    ki2sc,       & ! start index for second array computation
    ki2ec,       & ! end   index for second array computation
    ki3sc,       & ! start index for third  array computation
    ki3ec,       & ! end   index for third  array computation

    ! for activating debug output
    jindex,      & ! actual j-index computed
    idebug         ! debug control switch       

  LOGICAL                 , INTENT (IN) ::  &
    lcrf,        & ! control switch for cloud-free calcul.
    lsolar,      & ! control switch for solar calculations
    lradtopo       ! control switch for topographic corrections

  REAL    (KIND=ireals   ), INTENT (IN) ::  &
    psig,        & ! Stefan-Boltzman constant 
    psct           ! solar constant (at time of year)

  REAL    (KIND=ireals   ), INTENT (IN) ::  &
    ! Temperature at layer boundaries, pressure thickness of layers,
    ! cloud cover in each layer, water vapour and saturation water
    ! vapour mixing ratio, liquid water and ice mixing ratio
    ! layer CO2 content and O3 content
    pti    (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed+1)    , & ! (K)
    pdp    (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed)      , & ! (Pa)
    pclc   (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed)      , & ! (1)
    pwv    (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed)      , & ! (kg/kg)
    psw    (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed)      , & ! (kg/kg)
    pqlwc  (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed)      , & ! (kg/kg)
    pqiwc  (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed)      , & ! (kg/kg)
    pduco2 (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed)      , & ! (Pa CO2)
    pduo3  (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed)      , & ! (Pa O3)

    ! Aerorsole optical depth at 0.55 micrometer for five types
    paeq1  (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed)      , & ! (1)
    paeq2  (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed)      , & ! (1)
    paeq3  (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed)      , & ! (1)
    paeq4  (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed)      , & ! (1)
    paeq5  (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed)      , & ! (1)

    ! Cosine of zenith angle, thermal and solar surface albedo
    psmu0  (ki1sd:ki1ed,ki2sd:ki2ed)                  , & ! (1)
    palth  (ki1sd:ki1ed,ki2sd:ki2ed)                  , & ! (1)
    palso  (ki1sd:ki1ed,ki2sd:ki2ed)                  , & ! (1)
#ifdef COUP_OAS_COS
    palp   (ki1sd:ki1ed,ki2sd:ki2ed)                  , & ! (1)
         ! Solar surface albedo for parallel radiation
#endif
  
    ! External data and radiation corrections
    pskyview(ki1sd:ki1ed,ki2sd:ki2ed)                 , & !
    pfcor   (ki1sd:ki1ed,ki2sd:ki2ed)

  REAL    (KIND=ireals   ), INTENT (INOUT) ::  &
    ! Surface pressure
    papre  (ki1sd:ki1ed,ki2sd:ki2ed)                      ! (Pa)

! Output data
! -----------
  REAL    (KIND=ireals   ), INTENT (OUT) ::  &
    ! Thermal and solar radiative fluxes at each layer boundary
    ! dito for cloud-free conditions (TOA and surface only)
    pflt      (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed+1) , & ! (W/m**2)
    pfls      (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed+1) , & ! (W/m**2)

    ! corrected thermal and solar surface flux
    ! (if not lradtopo, just pflt(ke1) and pfls(ke1)
    pflt_s    (ki1sd:ki1ed,ki2sd:ki2ed)               , & ! (W/m**2)
    pfls_s    (ki1sd:ki1ed,ki2sd:ki2ed)               , & ! (W/m**2)

    ! and for the Climate LM Version: solar direct downward radiative flux 
    pflsdir   (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed+1) , & ! (W/m**2)

    ! components of solar and thermal fluxes at surface
    ! (influenced by lradtopo for topographic corrections)
    pfltd     (ki1sd:ki1ed,ki2sd:ki2ed)               , & ! (W/m**2)
    pfltu     (ki1sd:ki1ed,ki2sd:ki2ed)               , & ! (W/m**2)
    pflsd     (ki1sd:ki1ed,ki2sd:ki2ed)               , & ! (W/m**2)
    pflsu     (ki1sd:ki1ed,ki2sd:ki2ed)               , & ! (W/m**2)
    pflsp     (ki1sd:ki1ed,ki2sd:ki2ed)               , & ! (W/m**2)

    ! surface flux of photosynthetic active radiation and components
    pflpar    (ki1sd:ki1ed,ki2sd:ki2ed)               , & ! (W/m**2)
    pflsu_par (ki1sd:ki1ed,ki2sd:ki2ed)               , & ! (W/m**2)
    pflsd_par (ki1sd:ki1ed,ki2sd:ki2ed)               , & ! (W/m**2)
    pflsp_par (ki1sd:ki1ed,ki2sd:ki2ed)                   ! (W/m**2)

! Local parameters: 
! ----------------
  INTEGER (KIND=iintegers), PARAMETER ::  &
    j1b    = 1,           & ! debug point index first dimension
    j2b    = 1              ! debug point index second dimension

  REAL    (KIND=ireals   ), PARAMETER ::  &
    zepflx = 1.0E-8_ireals,     & ! Minimum 'grey' flux to avoid 1./0.
    zrd    = 287.05_ireals,     & ! Ra (gas constant of dry air)

    zrvdm1 = 461.51_ireals/287.05_ireals-1.0_ireals,  & ! Rv/Ra - 1
    zrvd   = 461.51_ireals/287.05_ireals,             & ! Rv/Ra
    zepai  = 0.0_ireals           ! Could be used to save computing time for 
                                  ! 'unimportant' gaseous absorption coefficients

! Local scalars:
! -------------
  INTEGER (KIND=iintegers) ::  &
    icrf, igase,           & !
    igasm1, igasz,         & !
    j1,j2,j3,              & ! loop indices over spatial dimensions
    jc,jh2o,jco2,jo3,      & ! loop indices over gaseous coefficients
    jspec,jspect,          & ! loop indices over spectrum
    jg,jjg ,               & ! loop indices over gases      
    icc,ih2o,ico2,io3        ! loop limit for gaseous absorption

  LOGICAL                  ::  &
    ldebug_th      ,       & ! debug control switch for thermal     
    ldebug_so      ,       & ! debug control switch for solar       
    ldebug_opt_th  ,       & ! debug control switch for opt_th      
    ldebug_opt_so  ,       & ! debug control switch for opt_so      
    ldebug_inv_th  ,       & ! debug control switch for inv_th      
    ldebug_inv_so            ! debug control switch for inv_so      

  REAL    (KIND=ireals   ) ::  &
    zet ,zaiprod,          & !
    zcoai,zcobi,           & !
    zemissivity, zalbedo     !

! Local arrays:
! -------------
  INTEGER (KIND=iintegers) ::  &
    icgas (3)                !

  REAL    (KIND=ireals   ) ::  &
    zketyp (jpther) ,      & !
    ztetyp (jpther)          !

! Local (automatic) arrays:
! ------------------------
! Arrays local to *fesft* or required for communication with
! subroutines called from *fesft*

  REAL    (KIND=ireals   ) ::  &
  ! 'Grey' and gaseous fluxes for individual spectral intervals
  ! "_c" means: corrected if lradtopo
  zflux    (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed+1), & ! (W/m**2)
  zflux_c  (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed+1), & ! (W/m**2)
  zfluxi   (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed+1), & ! 1./(W/m**2)
  zfluxu   (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed+1), & ! (W/m**2)
  zfluxu_c (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed+1), & ! (W/m**2)
  zfluxui  (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed+1), & ! 1./(W/m**2)
  zfluxd   (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed+1), & ! (W/m**2)
  zfluxd_c (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed+1), & ! (W/m**2)
  zfluxdi  (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed+1), & ! 1./(W/m**2)
  zfgas    (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed+1), & ! (W/m**2)
  zfgasu   (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed+1), & ! (W/m**2)
  zfgasd   (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed+1)    ! (W/m**2)

  REAL    (KIND=ireals   ) ::  &
  pbbr     (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed+1), &
         ! (W/m**2) Black body radiation at layer boundaries
  pflpt    (ki1sd:ki1ed,ki2sd:ki2ed)              , &
         ! Solar flux at TOA
#if ! defined COUP_OAS_COS
  palp     (ki1sd:ki1ed,ki2sd:ki2ed)              , &
         ! Solar surface albedo for parallel radiation
#endif
  pqsmu0   (ki1sd:ki1ed,ki2sd:ki2ed)              , &
         ! Inverse of cosine of zenith angle
  palogt   (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed)  , &
         ! ln T
  palogp   (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed)  , &
         ! ln p
  papra    (ki1sd:ki1ed,ki2sd:ki2ed)              , &
         ! (Pa) pressure at one level
  pduh2oc  (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed)  , &
         ! layer water vapour content (Pa), cloudy
  pduh2of  (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed)  , &
         ! layer water vapour content (Pa), cloud-free)
  pdulwc   (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed)  , &
         ! (Pa H2O-liquid) layer
         !incloud liquid water content
  pduiwc   (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed)  , &
         ! (Pa H2O-ice) layer incloud ice content
  prholwc  (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed)  , &
         ! (kg/m**3)
  prhoiwc  (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed)  , &
         ! (kg/m**3)
  zduetpc  (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed)  , &
         ! water vapour e-type contribution (cloudy)
  zduetpf  (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed)  , &
         ! water vapour e-type contribution (cloud-free)
 
  ! layer mean temperature, water vapour mixing ratio, utility arrays
  ztm      (ki1sd:ki1ed,ki2sd:ki2ed)              , & !
  zzwv     (ki1sd:ki1ed,ki2sd:ki2ed)              , & !
  zcpo     (ki1sd:ki1ed,ki2sd:ki2ed)              , & !
  zcpn     (ki1sd:ki1ed,ki2sd:ki2ed)              , & !
  zcmo     (ki1sd:ki1ed,ki2sd:ki2ed)              , & !
  zcmn     (ki1sd:ki1ed,ki2sd:ki2ed)                  !

  REAL    (KIND=ireals   ) ::  &
  ! Output data from opt_th/opt_so
  podac   (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed)   , & ! absorption optical depth
  podaf   (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed)   , & ! in cloudy and free part
  podsc   (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed)   , & ! scattering optical depth
  podsf   (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed)   , & ! in cloudy and free part
  pbsfc   (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed)   , & ! backscattering fraction 
  pbsff   (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed)   , & ! in cloudy and free part
  pusfc   (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed)   , & ! upscattering   fraction 
  pusff   (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed)   , & ! in cloudy and free part

  ! cloud geometry factors
  pca1    (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed)   , & !
  pcb1    (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed)   , & !
  pcc1    (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed)   , & !
  pcd1    (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed)   , & !
  pca2    (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed)   , & !
  pcb2    (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed)   , & !
  pcc2    (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed)   , & !
  pcd2    (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed)   , & !

  !fluxes calculated in inv_th/inv_so
  pflfd   (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed+1) , & ! (W/m**2)
  pflfu   (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed+1) , & ! (W/m**2)
  pflfp   (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed+1) , & ! (W/m**2)
  pflcd   (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed+1) , & ! (W/m**2)
  pflcu   (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed+1) , & ! (W/m**2)
  pflcp   (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed+1)     ! (W/m**2)

!- End of header
!==============================================================================
 
!------------------------------------------------------------------------------
! Begin Subroutine fesft               
!------------------------------------------------------------------------------
      
!------------------------------------------------------------------------------
! Section 1: Initializations
!------------------------------------------------------------------------------

  icrf  = 0
 
  ! Debug switches for lower level subroutines
  ldebug_th     = .FALSE.
  ldebug_so     = .FALSE.
  ldebug_opt_th = .FALSE.
  ldebug_opt_so = .FALSE.
  ldebug_inv_th = .FALSE.
  ldebug_inv_so = .FALSE.

  IF (idebug > 15) THEN
    PRINT *,' **** FESFT  *********************** debug point : ',j1b,j2b
  ENDIF 

  ! Preset output arrays
  DO j3 = ki3sc, ki3ec+1
    DO j2 = ki2sc, ki2ec
      DO j1 = ki1sc, ki1ec
        pflt   (j1,j2,j3) = 0.0_ireals
        pfls   (j1,j2,j3) = 0.0_ireals
        ! and for Climate-LM Version
        pflsdir(j1,j2,j3) = 0.0_ireals
      ENDDO
    ENDDO
  ENDDO

#ifdef COSMOART
  IF(l_cosmo_art) THEN
    DO j3 = ki3sc, ki3ec+1
      DO j1 = ki1sc, ki1ec  
        Edir (j1,jindex,j3) =  0.0_ireals
        Edown(j1,jindex,j3) =  0.0_ireals
        Eup  (j1,jindex,j3) =  0.0_ireals
      ENDDO
    ENDDO
  ENDIF
#endif

  DO j2 = ki2sc, ki2ec
    DO j1 = ki1sc, ki1ec
      pflpar   (j1,j2) = 0.0_ireals
      pflsu_par(j1,j2) = 0.0_ireals
      pflsd_par(j1,j2) = 0.0_ireals
      pflsp_par(j1,j2) = 0.0_ireals
      pfls_s   (j1,j2) = 0.0_ireals
      pflt_s   (j1,j2) = 0.0_ireals
      pfltu    (j1,j2) = 0.0_ireals
      pfltd    (j1,j2) = 0.0_ireals
      pflsu    (j1,j2) = 0.0_ireals
      pflsd    (j1,j2) = 0.0_ireals
      pflsp    (j1,j2) = 0.0_ireals
    ENDDO
  ENDDO
 
  ! Choice of e-type-absorption and temperature correction coefficients
  DO jspec = 1, jpther     
     zketyp(jspec) = zketypa(jspec)
     ztetyp(jspec) = ztetypa(jspec)
  ENDDO

  ! cloud geometry factors
 
  ! first part for top layer
  DO j2 = ki2sc, ki2ec
    DO j1 = ki1sc, ki1ec
      pca2(j1,j2,ki3sc)   = 1.0_ireals - pclc(j1,j2,ki3sc)
      pcd2(j1,j2,ki3sc)   = 1.0_ireals
      zcpn(j1,j2)         = MAX(pclc(j1,j2,ki3sc),pclc(j1,j2,ki3sc+1))
      zcmn(j1,j2)         = MIN(pclc(j1,j2,ki3sc),pclc(j1,j2,ki3sc+1))
      pca2(j1,j2,ki3sc+1) = (1.0_ireals-zcpn(j1,j2))/pca2(j1,j2,ki3sc)
      pcd2(j1,j2,ki3sc+1) =      zcmn(j1,j2) /pclc(j1,j2,ki3sc)
    ENDDO
  ENDDO
 
  ! first part for inner layers
 
  DO j3 = ki3sc+1, ki3ec-1
    DO j2 = ki2sc, ki2ec
      DO j1 = ki1sc, ki1ec
        zcpo(j1,j2)      = zcpn(j1,j2)
        zcmo(j1,j2)      = zcmn(j1,j2)
        zcpn(j1,j2)      = MAX (pclc(j1,j2,j3),pclc(j1,j2,j3+1))
        zcmn(j1,j2)      = MIN (pclc(j1,j2,j3),pclc(j1,j2,j3+1))
        pca2(j1,j2,j3+1) = (1.-zcpn(j1,j2))/(1.-pclc(j1,j2,j3))
        pca1(j1,j2,j3-1) = (1.-zcpo(j1,j2))/(1.-pclc(j1,j2,j3))
        pcd2(j1,j2,j3+1) = zcmn(j1,j2)/pclc(j1,j2,j3)
        pcd1(j1,j2,j3-1) = zcmo(j1,j2)/pclc(j1,j2,j3)
      ENDDO
    ENDDO
  ENDDO
 
  ! first part for lowest layer
  DO j2 = ki2sc, ki2ec
    DO j1 = ki1sc, ki1ec
      pca1(j1,j2,ki3ec-1) = (1.-zcpn(j1,j2))/(1.-pclc(j1,j2,ki3ec))
      pcd1(j1,j2,ki3ec-1) = zcmn(j1,j2)/pclc(j1,j2,ki3ec)
      pca1(j1,j2,ki3ec  ) = 1.0_ireals
      pcd1(j1,j2,ki3ec  ) = 1.0_ireals
    ENDDO
  ENDDO
 
  ! second part of geometry factors
  DO j3 = ki3sc, ki3ec
    DO j2 = ki2sc, ki2ec
      DO j1 = ki1sc, ki1ec
        pcb1(j1,j2,j3) = 1.-pca1(j1,j2,j3)
        pcc1(j1,j2,j3) = 1.-pcd1(j1,j2,j3)
        pcb2(j1,j2,j3) = 1.-pca2(j1,j2,j3)
        pcc2(j1,j2,j3) = 1.-pcd2(j1,j2,j3)
      ENDDO
    ENDDO
  ENDDO
 
  ! Optically relevant layer constituents
  ! (Note: CO2 and O3 amounts are provided by calling routine)
  DO j2 = ki2sc, ki2ec
    DO j1 = ki1sc, ki1ec
      papra(j1,j2) = papre(j1,j2)   ! surface pressure
    ENDDO
  ENDDO
 
  ! water vapour, liquid water and ice content, logarithm of layer
  ! mean temperature and pressure, absorber amount for e-type absorption
 
  DO j3 = ki3ec, ki3sc,-1    ! Bottom to top
    DO j2 = ki2sc, ki2ec
      DO j1 = ki1sc, ki1ec
        ztm   (j1,j2)    = 0.5*(pti(j1,j2,j3)+pti(j1,j2,j3+1))
        papra (j1,j2)    = papra(j1,j2) - 0.5*pdp(j1,j2,j3)
        palogt(j1,j2,j3) = LOG  (ztm  (j1,j2))
        palogp(j1,j2,j3) = LOG  (papra(j1,j2))
    
        ! cloud-free:  water vapour and e-type absorber amount
        zzwv   (j1,j2   ) = MAX( (pwv(j1,j2,j3)-pclc(j1,j2,j3)*psw(j1,j2,j3)) &
                              /(1.0_ireals-pclc(j1,j2,j3)) , 0.0_ireals)
        pduh2of(j1,j2,j3) = pdp(j1,j2,j3)*zzwv(j1,j2)
        zduetpf(j1,j2,j3) = pduh2of(j1,j2,j3)*pduh2of(j1,j2,j3) &
                           *papra(j1,j2)*zrvd/pdp(j1,j2,j3)
    
        ! cloudy:  water vapour, e-type absorber amount, liquid water and ice
        pdulwc (j1,j2,j3) = pdp(j1,j2,j3) * (pqlwc(j1,j2,j3)/pclc(j1,j2,j3))
        pdulwc (j1,j2,j3) = MAX( pdulwc(j1,j2,j3), 0.0_ireals )
        pduiwc (j1,j2,j3) = pdp(j1,j2,j3) * (pqiwc(j1,j2,j3)/pclc(j1,j2,j3))
        pduiwc (j1,j2,j3) = MAX( pduiwc(j1,j2,j3), 0.0_ireals )
        pduh2oc(j1,j2,j3) = pdp(j1,j2,j3) * psw(j1,j2,j3)
        zduetpc(j1,j2,j3) = pduh2oc(j1,j2,j3) * pduh2oc(j1,j2,j3) *           &
                                  papra(j1,j2) * zrvd / pdp(j1,j2,j3)
        prholwc(j1,j2,j3) = (pqlwc(j1,j2,j3) / pclc(j1,j2,j3)) * papra(j1,j2) &
                              / (zrd*ztm(j1,j2) * (1.+zrvdm1*psw(j1,j2,j3)))
        prhoiwc(j1,j2,j3) = (pqiwc(j1,j2,j3) / pclc(j1,j2,j3)) * papra(j1,j2) &
                              / (zrd*ztm(j1,j2) * (1.+zrvdm1*psw(j1,j2,j3)))
        ! Secure minium for ice density for use in empirical function with ALOG
        prhoiwc(j1,j2,j3) = MAX (prhoiwc(j1,j2,j3),1.0E-06_ireals) 
         
        papra  (j1,j2   ) = papra(j1,j2) - 0.5_ireals * pdp(j1,j2,j3)
      ENDDO
    ENDDO
  ENDDO     ! End of vertical loop

  ! Identify *papre* with top of model pressure (for Rayleigh scattering)
  DO j2 = ki2sc, ki2ec
    DO j1 = ki1sc, ki1ec
      papre(j1,j2) = papra(j1,j2)
    ENDDO
  ENDDO
 
#ifdef COSMOART
! aerosol optical properties for online radiation feedback
  IF(l_cosmo_art) THEN
    IF (lrad_dust) CALL rad_dust(jindex)
    IF (lrad_seas) CALL rad_seas(jindex)
    IF (lrad_aero) CALL rad_aero(jindex)
  ENDIF
#endif

1 CONTINUE  ! Address for backward jump to perform cloud-free calculations

!------------------------------------------------------------------------------
! Section 2: Thermal radiative flux calculations 
!------------------------------------------------------------------------------

  ! Loop over thermal spectral intervals
  !================================================================
  thermal_spectral_loop:  DO jspec=jpsol+1,jpspec
  !================================================================
      
  !----------------------------------------------------------------------------
  ! Section 2.1: Initializations
  !----------------------------------------------------------------------------

    jspect = jspec - jpsol
 
    ! Black body radiation at layer boundaries in spectral interval
    DO j3 = ki3sc, ki3ec+1
      DO j2 = ki2sc, ki2ec
        DO j1 = ki1sc, ki1ec
          pbbr(j1,j2,j3)= ( planck(1,jspect) + pti(j1,j2,j3)                 &
                   * ( planck(2,jspect) + pti(j1,j2,j3)*planck(3,jspect) ) ) &
                   * psig * (pti(j1,j2,j3)**2)**2
        ENDDO
      ENDDO
    ENDDO
 
    ! Optical properties of non-gaseous constituents        
    IF (idebug > 10 ) THEN
       print *,' FESFT    Call to opt_th for jspec: ',jspec
    ENDIF 

    CALL opt_th( prholwc,pdulwc,prhoiwc,pduiwc,              &
                 paeq1  ,paeq2 ,paeq3  ,paeq4 , paeq5 ,      &
                 ki1sd ,ki1ed,ki2sd,ki2ed,ki3sd,ki3ed,       &
                 jspec  ,ki1sc ,ki1ec  ,ki2sc , ki2ec ,      &
                         ki3sc ,ki3ec  ,ldebug_opt_th,       &
                 podac  ,podaf ,podsc  ,podsf , pbsfc ,pbsff )

    ! Addition of e-type contribution
    IF (zketyp(jspect) /= 0.) THEN
      zet   = 1./EXP(ztetyp(jspect)/zteref)
      DO j3 = ki3sc,ki3ec
        DO j2 = ki2sc, ki2ec
          DO j1 = ki1sc, ki1ec
            ztm(j1,j2)      = 0.5*(pti(j1,j2,j3)+pti(j1,j2,j3+1))
            podaf(j1,j2,j3) = podaf(j1,j2,j3) + zduetpf(j1,j2,j3) &
                      * zet * EXP(ztetyp(jspect)/ztm(j1,j2)) * zketyp(jspect)
            podac(j1,j2,j3) = podac(j1,j2,j3) + zduetpc(j1,j2,j3) &
                      * zet * EXP(ztetyp(jspect)/ztm(j1,j2)) * zketyp(jspect)
          ENDDO  
        ENDDO  
      ENDDO  
    ENDIF

  !----------------------------------------------------------------------------
  ! Section 2.2: Selection of ESFT or FESFT method for interval considered
  !----------------------------------------------------------------------------

    !--------------------------------------------------------------
    IF (nfast(jspec) == 0) THEN           ! ESFT method
    !--------------------------------------------------------------

    DO j3 = ki3sc, ki3ec+1
      DO j2 = ki2sc, ki2ec
        DO j1 = ki1sc, ki1ec
          zflux(j1,j2,j3) = 0.0_ireals  ! Preset flux in spectral interval
        ENDDO
      ENDDO
    ENDDO
 
      ! Loop over various absorption coefficients of each of the three gases
 
      ih2o = ncgas(jspec,1)
      ico2 = ncgas(jspec,2)
      io3  = ncgas(jspec,3)
 
      DO jh2o=    1,ih2o    ! Loop over H2O coefficients
        DO jco2=  1,ico2    ! Loop over CO2 coefficients
          DO jo3= 1,io3     ! Loop over O3  coefficients
 
          zaiprod = coai(jh2o,jspec,1)*coai(jco2,jspec,2)*coai(jo3,jspec,3)
 
          IF (icrf.eq.0) THEN      ! partially cloudy atmosphere
            call       inv_th (                                              &
              pclc   ,pca1   ,pca2  ,pcb1  ,pcb2  ,pcc1  ,pcc2 ,pcd1 ,pcd2 , &
              pduh2oc,pduh2of,pduco2,pduo3 ,palogp,palogt,                   &
              podsc  ,podsf  ,podac ,podaf ,pbsfc ,pbsff ,                   &
              pbbr   ,palth,                                                 &
              jspec  ,jh2o   ,jco2  ,jo3   ,                                 &
              ki1sd  ,ki1ed  ,ki2sd ,ki2ed ,ki3sd ,ki3ed ,                   &
              ki1sc  ,ki1ec  ,ki2sc ,ki2ec ,ki3sc ,ki3ec , ldebug_inv_th ,   &
              pflcu  ,pflfu  ,pflcd ,pflfd)
          ELSE                     ! 'cloud-free' atmosphere    
            print *,' CRF not yet implemented'
          END IF
 
          ! Incrementation of flux in spectral interval

          DO j3 = ki3sc, ki3ec+1
            DO j2 = ki2sc, ki2ec
              DO j1 = ki1sc, ki1ec
                zflux(j1,j2,j3) = zflux(j1,j2,j3)                            &
                  + zaiprod * ( pflfu(j1,j2,j3) + pflcu(j1,j2,j3)            &
                              - pflfd(j1,j2,j3) - pflcd(j1,j2,j3) )
              ENDDO
            ENDDO
          ENDDO

          ENDDO       ! Loop over O3 absorption coefficients
        ENDDO         ! Loop over CO2 absorption coefficients
      ENDDO           ! Loop over H2O absorption coefficients

    !--------------------------------------------------------------
    ELSE                                  ! FESFT method
    !--------------------------------------------------------------

      igase = 0
      DO jg=1,3
        icgas (jg) = 0
        IF (ncgas(jspec,jg).GT.1) THEN
          igase     = igase + 1
        ENDIF
      ENDDO
      igasm1    = igase -1
 
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF (igase.le.1) THEN  !(no 'grey' fluxes required)
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        DO j3 = ki3sc, ki3ec+1
          DO j2 = ki2sc, ki2ec
            DO j1 = ki1sc, ki1ec
              zfluxi(j1,j2,j3) = 1.0_ireals
            ENDDO
          ENDDO
        ENDDO

      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ELSE   ! more than 1 gas --> 'grey' fluxes required
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        IF (icrf.EQ.0) THEN    ! partially cloudy atmosphere
          call       inv_th (                                               &
             pclc   ,pca1   ,pca2  ,pcb1  ,pcb2  ,pcc1  ,pcc2 ,pcd1 ,pcd2 , &
             pduh2oc,pduh2of,pduco2,pduo3 ,palogp,palogt,                   &
             podsc  ,podsf  ,podac ,podaf ,pbsfc ,pbsff ,                   &
             pbbr   ,palth,                                                 &
             jspec  ,0      ,0     ,0     ,                                 &
             ki1sd  ,ki1ed  ,ki2sd ,ki2ed ,ki3sd ,ki3ed ,                   &
             ki1sc  ,ki1ec  ,ki2sc ,ki2ec ,ki3sc ,ki3ec , ldebug_inv_th ,   &
             pflcu  ,pflfu  ,pflcd ,pflfd)
        ELSE                     ! 'cloud-free' atmosphere    
          print *,' CRF not yet implemented'
        ENDIF

        ! Storage of 'grey' fluxes and their inverse (**igasm1)
        DO j3 = ki3sc, ki3ec+1
          DO j2 = ki2sc, ki2ec
            DO j1 = ki1sc, ki1ec
              zflux (j1,j2,j3) = pflfu(j1,j2,j3) + pflcu(j1,j2,j3)          &
                               - pflfd(j1,j2,j3) - pflcd(j1,j2,j3)
              zfluxi(j1,j2,j3) = 1.0_ireals                                 &
                                 / MIN( -zepflx, zflux(j1,j2,j3) )**igasm1
            ENDDO
          ENDDO
        ENDDO

      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ENDIF    ! No.of relevant gases
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      igasz = 0
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      DO jg = 3, 1, -1   !     Loop over gases
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 
        IF (ncgas(jspec,jg).GT.1) THEN   ! include gas only, if necessary
          igasz = igasz + 1
 
          DO jjg = 1,3
            icgas(jjg) = 0      !   Set absorption coefficient index for all
          ENDDO                 !   gases to zero
 
          DO j3 = ki3sc, ki3ec+1
            DO j2 = ki2sc, ki2ec
              DO j1 = ki1sc, ki1ec
                zfgas(j1,j2,j3) = 0.0_ireals    ! Preset 'gaseous' flux
              ENDDO
            ENDDO
          ENDDO
 
          icc = ncgas(jspec,jg) ! No.of relevant coefficients    
 
          ! -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
          DO jc = icc,1,-1        ! Loop over absorption coefficients
          ! -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
            zcoai = coai(jc,jspec,jg)
            zcobi = cobi(jc,jspec,jg)
 
            ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            IF ( ((zcoai.GE.zepai).AND.(zcobi.GT.0.0)) .OR. (igase.EQ.1) ) THEN   
            ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
              ! Solve linear system, if necessary
              icgas(jg) = jc
              IF (icrf.EQ.0) THEN  ! partially cloudy atmosphere
                call       inv_th (                                              & 
                  pclc   ,pca1   ,pca2  ,pcb1  ,pcb2  ,pcc1  ,pcc2 ,pcd1 ,pcd2 , &
                  pduh2oc,pduh2of,pduco2,pduo3 ,palogp,palogt,                   &
                  podsc  ,podsf  ,podac ,podaf ,pbsfc ,pbsff ,                   &
                  pbbr   ,palth,                                                 &
                  jspec  ,icgas(1),icgas(2),icgas(3),                            &
                  ki1sd  ,ki1ed  ,ki2sd ,ki2ed ,ki3sd ,ki3ed ,                   &
                  ki1sc  ,ki1ec  ,ki2sc ,ki2ec ,ki3sc ,ki3ec , ldebug_inv_th ,   &
                  pflcu  ,pflfu  ,pflcd ,pflfd)
              ELSE                     ! 'cloud-free' atmosphere    
                print *,' CRF not yet implemented'
              ENDIF
            ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            ELSE               ! use 'grey' fluxes directly
            ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
              DO j3 = ki3sc, ki3ec+1
                DO j2 = ki2sc, ki2ec
                  DO j1 = ki1sc, ki1ec
                    pflfu(j1,j2,j3) = zflux(j1,j2,j3)
                    pflfd(j1,j2,j3) = 0.0_ireals
                    pflcu(j1,j2,j3) = 0.0_ireals
                    pflcd(j1,j2,j3) = 0.0_ireals
                  ENDDO
                ENDDO
              ENDDO
            ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            ENDIF            ! Necessity to calculate fluxes
            ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 
            DO j3 = ki3sc, ki3ec+1
              DO j2 = ki2sc, ki2ec
                DO j1 = ki1sc, ki1ec
                  zfgas(j1,j2,j3) = zfgas(j1,j2,j3) + zcoai * ( pflfu(j1,j2,j3) &
                      + pflcu(j1,j2,j3) - pflfd(j1,j2,j3) - pflcd(j1,j2,j3) )
                ENDDO
              ENDDO
            ENDDO

            !?????????????????????????????????????????????????????????????????
            IF (ldebug_th) THEN
               print *,' FESFT in debug mode for thermal fluxes'
               print *,' only one interval/coefficient considered '
               print *,'zfgas(j1b,j2b,ki3sc): ',zfgas(j1b,j2b,ki3sc)
               EXIT thermal_spectral_loop
            ENDIF 
            !?????????????????????????????????????????????????????????????????
          ! -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
          ENDDO      ! Loop over absorption coefficients
          ! -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
 
          ! Combination of inverse of 'grey' fluxes and actual gaseous flux
          DO j3 = ki3sc, ki3ec+1
            DO j2 = ki2sc, ki2ec
              DO j1 = ki1sc, ki1ec
                zfluxi(j1,j2,j3) = zfluxi(j1,j2,j3)*zfgas(j1,j2,j3)
              ENDDO
            ENDDO
          ENDDO
 
          IF (igasz.eq.igasm1) THEN    ! Avoid unphysical pseudo-transmission
            DO j3 = ki3sc, ki3ec+1
              DO j2 = ki2sc, ki2ec
                DO j1 = ki1sc, ki1ec
                  zfluxi(j1,j2,j3) = MIN( 1.0_ireals,                         &
                                          MAX( 0.0_ireals, zfluxi(j1,j2,j3) ) )
                ENDDO
              ENDDO
            ENDDO
          ENDIF
    
        ENDIF                   ! Test, whether gas needs to be included
    
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      END DO          ! End of loop over gases
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 
      ! Store FESFT result in zflux
      DO j3 = ki3sc, ki3ec+1
        DO j2 = ki2sc, ki2ec
          DO j1 = ki1sc, ki1ec
            zflux(j1,j2,j3) = zfluxi(j1,j2,j3)
          ENDDO
        ENDDO
      ENDDO

    !--------------------------------------------------------------
    END IF                             ! ESFT/FESFT-Selection 
    !--------------------------------------------------------------

  !----------------------------------------------------------------------------
  ! Section 2.3: Addition of flux for spectral interval to total thermal flux
  !----------------------------------------------------------------------------
 
    IF (icrf.EQ.0) THEN     ! Add flux at all levels
      DO j3 = ki3sc, ki3ec+1
        DO j2 = ki2sc, ki2ec
          DO j1 = ki1sc, ki1ec
            pflt(j1,j2,j3) = pflt(j1,j2,j3) + zflux(j1,j2,j3)
          ENDDO
        ENDDO
      ENDDO
    END IF
 
  ! End of spectral loop
  ! ===============================================================
  ENDDO  thermal_spectral_loop
  ! ===============================================================

  !----------------------------------------------------------------------------
  ! Section 2.4: Storage of components for thermal radiative surface flux
  !----------------------------------------------------------------------------

  DO j2 = ki2sc, ki2ec
    DO j1 = ki1sc, ki1ec

      ! Recompute surface thermal flux components based on lower boundary 
      ! condition (cf. Ritter and Geleyn (1992))
      zemissivity   =  1.0_ireals-palth(j1,j2) ! surface emissivity
      pfltd (j1,j2) = (pflt(j1,j2,ki3ec+1) +                               &
                             zemissivity * psig * pti(j1,j2,ki3ec+1)**4)   &
                        / zemissivity
      pfltu (j1,j2) = pfltd(j1,j2) - pflt(j1,j2,ki3ec+1)
      pflt_s(j1,j2) = pflt(j1,j2,ki3ec+1)
    ENDDO
  ENDDO

!NEC_CB Moved Debug-Prints of out the loop
  IF (idebug > 15) THEN
    IF ((j2b>=ki2sc).AND.(j2b<=ki2ec).AND.(j1b>=ki1sc).AND.(j1b<=ki1ec)) THEN
      j1=j1b
      j2=j2b
      WRITE (*,'(A32,2F16.6)') 'FESFT: zemissivity, palth',              &
           zemissivity, palth(j1,j2)
      WRITE (*,'(A60, F16.6)')                                           &
          'FESFT: thermal fluxes before and after correction, skyview',  &
           pskyview(j1,j2)
      WRITE (*,'(A32,2I4,3F16.6)') 'FESFT th before: net, down, up:',    &
           j1, j2, pflt(j1,j2,ki3ec+1), pfltd(j1,j2), pfltu(j1,j2)
    ENDIF
  ENDIF

  IF (lradtopo) THEN
    DO j2 = ki2sc, ki2ec
      DO j1 = ki1sc, ki1ec
        ! corrected thermal balance
        ! correction as in Mueller and Scherrer (2005)
        pfltd (j1,j2) = pfltd(j1,j2) *               pskyview(j1,j2)  +   &
                        pfltu(j1,j2) * (1.0_ireals - pskyview(j1,j2))
        pflt_s(j1,j2) = pfltd(j1,j2)-pfltu(j1,j2)

      ENDDO
    ENDDO
!NEC_CB Moved Debug-Prints of out the loop
    IF (idebug > 15) THEN
      IF ((j2b>=ki2sc).AND.(j2b<=ki2ec).AND.(j1b>=ki1sc).AND.(j1b<=ki1ec)) THEN
        j1=j1b
        j2=j2b
        WRITE (*,'(A32,2I4,3F16.6)') 'FESFT th after: net,down,up:',     &
                j1, j2, pflt_s(j1,j2), pfltd(j1,j2), pfltu(j1,j2)
      END IF
    END IF
  END IF

!------------------------------------------------------------------------------
! Section 3: Solar flux calculations, if required
!------------------------------------------------------------------------------

  IF (lsolar) THEN

  !----------------------------------------------------------------------------
  ! Section 3.1: Initializations
  !----------------------------------------------------------------------------

  ! Inverse of cosine of zenith angle and surface albedo for
  ! parallel radiation

  DO j2 = ki2sc, ki2ec
    DO j1 = ki1sc, ki1ec
      pqsmu0(j1,j2) = 1.0_ireals / psmu0(j1,j2)
#if ! defined COUP_OAS_COS
      palp  (j1,j2) = (1.0_ireals +                                           &
        0.5_ireals * (psmu0(j1,j2) * (1.0_ireals/palso(j1,j2) - 1.0_ireals))) &
     / (1.0_ireals + (psmu0(j1,j2) * (1.0_ireals/palso(j1,j2) - 1.0_ireals)))**2
#endif
    ENDDO
  ENDDO

  ! Loop over solar spectral intervals
  ! ===============================================================
  solar_spectral_loop:  DO jspec = 1, jpsol
  ! ===============================================================

    ! Preset flux in spectral interval
    DO j3 = ki3sc, ki3ec+1
      DO j2 = ki2sc, ki2ec
        DO j1 = ki1sc, ki1ec
          zflux (j1,j2,j3) = 0.0_ireals
          zfluxd(j1,j2,j3) = 0.0_ireals
          zfluxu(j1,j2,j3) = 0.0_ireals
        ENDDO
      ENDDO
    ENDDO
 
    ! Upper boundary condition and reference pressure for Rayleigh sc.
    DO j2 = ki2sc, ki2ec
      DO j1 = ki1sc, ki1ec
        pflpt(j1,j2) = psct * solant(jspec) * psmu0(j1,j2)
        papra(j1,j2) = papre(j1,j2) * pqsmu0(j1,j2)
      ENDDO
    ENDDO
 
    ! Optical properties of non-gaseous constituents        
 
    IF (idebug > 10) THEN
       print *,' FESFT    Call to opt_so for jspec: ',jspec
    ENDIF   

    CALL opt_so ( prholwc,pdulwc,prhoiwc,pduiwc,               &
                  paeq1  ,paeq2 ,paeq3  ,paeq4 , paeq5 ,       &
                  pdp    ,papra ,psmu0  ,pqsmu0,               &
                  ki1sd  ,ki1ed,ki2sd,ki2ed,ki3sd,ki3ed,       &
                  jspec  ,ki1sc ,ki1ec  ,ki2sc , ki2ec ,       &
                          ki3sc ,ki3ec  ,ldebug_opt_so,        &
                  podac  ,podaf ,podsc  ,podsf , pbsfc ,pbsff ,&
                  pusfc  ,pusff  )

  !----------------------------------------------------------------------------
  ! Section 3.2: Selection of ESFT or FESFT method for interval considered
  !----------------------------------------------------------------------------

    !--------------------------------------------------------------
    IF (nfast(jspec).eq.0) THEN           ! ESFT method
    !--------------------------------------------------------------

      ih2o = ncgas(jspec,1)
      ico2 = ncgas(jspec,2)
      io3  = ncgas(jspec,3)
 
      DO jh2o    = 1, ih2o    ! Loop over H2O coefficients
        DO jco2  = 1, ico2    ! Loop over CO2 coefficients
          DO jo3 = 1, io3     ! Loop over O3  coefficients
 
            zaiprod=coai(jh2o,jspec,1)*coai(jco2,jspec,2)*coai(jo3,jspec,3)
 
            IF (icrf.eq.0) THEN      ! partially cloudy atmosphere
              CALL inv_so (                                                     &
                 pclc   ,pca1   ,pca2  ,pcb1  ,pcb2  ,pcc1  ,pcc2 ,pcd1 ,pcd2 , &
                 pflpt  ,psmu0  ,pqsmu0,palp  ,palso ,                          &
                 pduh2oc,pduh2of,pduco2,pduo3 ,palogp,palogt,                   &
                 podsc  ,podsf  ,podac ,podaf ,pbsfc ,pbsff ,pusfc,pusff,       &
                 jspec  ,jh2o   ,jco2  ,jo3   ,                                 &
                 ki1sd  ,ki1ed  ,ki2sd ,ki2ed ,ki3sd ,ki3ed,                    &
                 ki1sc  ,ki1ec  ,ki2sc ,ki2ec ,ki3sc ,ki3ec , ldebug_inv_so ,   &
                 pflcu  ,pflfu  ,pflcd ,pflfd ,pflcp ,pflfp)
 
            ELSE                     ! cloud-free calculation
              print *,' CRF-Code not implemented yet'
            ENDIF
 
            ! Incrementation of flux in spectral interval
            DO j3 = ki3sc, ki3ec+1
              DO j2 = ki2sc, ki2ec
                DO j1 = ki1sc, ki1ec
                  zflux (j1,j2,j3) = zflux (j1,j2,j3)                           &
                            + zaiprod * (pflfp(j1,j2,j3) + pflcp(j1,j2,j3))
                  zfluxd(j1,j2,j3) = zfluxd(j1,j2,j3)                           &
                            + zaiprod * (pflfd(j1,j2,j3) + pflcd(j1,j2,j3))
                  zfluxu(j1,j2,j3) = zfluxu(j1,j2,j3)                           &
                            + zaiprod * (pflfu(j1,j2,j3) + pflcu(j1,j2,j3))
                ENDDO
              ENDDO
            ENDDO
 
          ENDDO       ! Loop over O3 absorption coefficients
        ENDDO         ! Loop over CO2 absorption coefficients
      ENDDO           ! Loop over H2O absorption coefficients

    !--------------------------------------------------------------
    ELSE                                  ! FESFT method
    !--------------------------------------------------------------

      igase = 0
      DO jg = 1, 3
         icgas(jg) = 0
         IF (ncgas(jspec,jg).GT.1) THEN
            igase = igase + 1
         ENDIF
      ENDDO

      igasm1    = igase -1
 
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF (igase.le.1) THEN  !(no 'grey' fluxes required)
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        DO j3 = ki3sc, ki3ec+1
          DO j2 = ki2sc, ki2ec
            DO j1 = ki1sc, ki1ec
              zflux  (j1,j2,j3) = 1.0_ireals
              zfluxd (j1,j2,j3) = 1.0_ireals
              zfluxu (j1,j2,j3) = 1.0_ireals
              zfluxi (j1,j2,j3) = 1.0_ireals
              zfluxdi(j1,j2,j3) = 1.0_ireals
              zfluxui(j1,j2,j3) = 1.0_ireals
            ENDDO
          ENDDO
        ENDDO

      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ELSE   ! more than 1 gas --> 'grey' fluxes required
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        IF (icrf.eq.0) THEN      ! partially cloudy atmosphere
          CALL inv_so (                                                    &
             pclc   ,pca1   ,pca2  ,pcb1  ,pcb2  ,pcc1  ,pcc2 ,pcd1 ,pcd2 ,&
             pflpt  ,psmu0  ,pqsmu0,palp  ,palso ,                         &
             pduh2oc,pduh2of,pduco2,pduo3 ,palogp,palogt,                  & 
             podsc  ,podsf  ,podac ,podaf ,pbsfc ,pbsff ,pusfc,pusff,      &
             jspec  ,0      ,0     ,0     ,                                &
             ki1sd  ,ki1ed  ,ki2sd ,ki2ed ,ki3sd ,ki3ed,                   &
             ki1sc  ,ki1ec  ,ki2sc ,ki2ec ,ki3sc ,ki3ec , ldebug_inv_so ,  &
             pflcu  ,pflfu  ,pflcd ,pflfd ,pflcp ,pflfp)
        ELSE                     ! cloud-free calculation
          print *,' CRF-Code not implemented yet !!!'
          stop 'no-crf'
        ENDIF

        ! Storage of 'grey' fluxes and their inverse (**igasm1)
        DO j3 = ki3sc, ki3ec+1
          DO j2 = ki2sc, ki2ec
            DO j1 = ki1sc, ki1ec
              zfluxi (j1,j2,j3) = 1.0_ireals                                  &
                   / MAX (pflfp(j1,j2,j3)+pflcp(j1,j2,j3), zepflx) **  igasm1
              zfluxdi(j1,j2,j3) = 1.0_ireals                                  &
                   / MAX (pflfd(j1,j2,j3)+pflcd(j1,j2,j3), zepflx) ** igasm1
              zfluxui(j1,j2,j3) = 1.0_ireals                                  &
                   / MAX (pflfu(j1,j2,j3)+pflcu(j1,j2,j3), zepflx) ** igasm1
              zflux  (j1,j2,j3) = pflfp(j1,j2,j3) + pflcp(j1,j2,j3)
              zfluxd (j1,j2,j3) = pflfd(j1,j2,j3) + pflcd(j1,j2,j3)
              zfluxu (j1,j2,j3) = pflfu(j1,j2,j3) + pflcu(j1,j2,j3)
            ENDDO
          ENDDO
        ENDDO
      
        IF (ldebug_so) THEN
           print *,' FESFT in debug mode for solar   fluxes'
           print *,' Grey fluxes   '
           DO j3=ki3sc,ki3ec+1
              print *,'par/down/up    : ',zflux (j1b,j2b,j3) &
                                         ,zfluxd(j1b,j2b,j3) &
                                         ,zfluxu(j1b,j2b,j3),j3
           ENDDO
        ENDIF 
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ENDIF    ! No.of relevant gases
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      igasz = 0
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      DO jg = 3, 1, -1   !     Loop over gases
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 
      IF (ncgas(jspec,jg).GT.1) THEN   ! include gas only, if necessary
   
        igasz = igasz + 1
 
        DO jjg = 1,3
           icgas(jjg) = 0
        ENDDO 
 
        ! Initialize 'gaseous' fluxes
        DO j3 = ki3sc, ki3ec+1
          DO j2 = ki2sc, ki2ec
            DO j1 = ki1sc, ki1ec
              zfgas (j1,j2,j3) = 0.0_ireals
              zfgasd(j1,j2,j3) = 0.0_ireals
              zfgasu(j1,j2,j3) = 0.0_ireals
            ENDDO
          ENDDO
        ENDDO
 
        ! -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
        icc = ncgas(jspec,jg)
        DO jc = icc,1,-1        ! Loop over absorption coefficients
        ! -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
          zcoai = coai(jc,jspec,jg)
          zcobi = cobi(jc,jspec,jg)
          ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
          IF ( ((zcoai.GE.zepai).AND.(zcobi.GT.0.0)) .OR. (igase.EQ.1) ) THEN 
          ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            ! Solve linear system, if necessary
            icgas(jg) = jc
            IF (icrf.eq.0) THEN      ! partially cloudy atmosphere
              CALL inv_so (                                                     &
                 pclc   ,pca1   ,pca2  ,pcb1  ,pcb2  ,pcc1  ,pcc2 ,pcd1 ,pcd2 , &
                 pflpt  ,psmu0  ,pqsmu0,palp  ,palso ,                          &
                 pduh2oc,pduh2of,pduco2,pduo3 ,palogp,palogt,                   &
                 podsc  ,podsf  ,podac ,podaf ,pbsfc ,pbsff ,pusfc,pusff,       &
                 jspec  ,icgas(1),icgas(2),icgas(3),                            &
                 ki1sd  ,ki1ed  ,ki2sd ,ki2ed ,ki3sd ,ki3ed,                    &
                 ki1sc  ,ki1ec  ,ki2sc ,ki2ec ,ki3sc ,ki3ec , ldebug_inv_so ,   &
                 pflcu  ,pflfu  ,pflcd ,pflfd ,pflcp ,pflfp)
            ELSE                     ! cloud-free clculations
              print *,'crf-code not yet implemented'
              stop 'no-crf'
            ENDIF

          ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
          ELSE               ! use 'grey' fluxes directly
          ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            DO j3 = ki3sc, ki3ec+1
              DO j2 = ki2sc, ki2ec
                DO j1 = ki1sc, ki1ec
                  pflfp(j1,j2,j3) = zflux (j1,j2,j3)
                  pflcp(j1,j2,j3) = 0.0_ireals
                  pflfd(j1,j2,j3) = zfluxd(j1,j2,j3)
                  pflcd(j1,j2,j3) = 0.0_ireals
                  pflfu(j1,j2,j3) = zfluxu(j1,j2,j3)
                  pflcu(j1,j2,j3) = 0.0_ireals
                ENDDO
              ENDDO
            ENDDO

          ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
          ENDIF            ! Necessity to calculate fluxes
          ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 
          DO j3 = ki3sc, ki3ec+1
            DO j2 = ki2sc, ki2ec
              DO j1 = ki1sc, ki1ec
                zfgas (j1,j2,j3) = zfgas (j1,j2,j3)                          &
                                 + zcoai * (pflfp(j1,j2,j3) + pflcp(j1,j2,j3))
                zfgasd(j1,j2,j3) = zfgasd(j1,j2,j3)                          &
                                 + zcoai * (pflfd(j1,j2,j3) + pflcd(j1,j2,j3))
                zfgasu(j1,j2,j3) = zfgasu(j1,j2,j3)                          &
                                 + zcoai * (pflfu(j1,j2,j3) + pflcu(j1,j2,j3))
              ENDDO
            ENDDO
          ENDDO
 
          !?????????????????????????????????????????????????????????????????
          IF (ldebug_so) THEN
             print *,' FESFT in debug mode for solar   fluxes'
             print *,' only one interval/coefficient considered '
             print *,' zcoai = ',zcoai
             DO j3=ki3sc,ki3ec+1
                print *,'zfgas(j1b,j2b,j3): ',zfgas(j1b,j2b,j3)
             ENDDO 
             EXIT solar_spectral_loop
          ENDIF 
          !?????????????????????????????????????????????????????????????????
        ! -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
        ENDDO      ! Loop over absorption coefficients
        ! -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
 
        ! Combination of inverse of 'grey' fluxes and actual gaseous flux
        DO j3 = ki3sc, ki3ec+1
          DO j2 = ki2sc, ki2ec
            DO j1 = ki1sc, ki1ec
              zfluxi (j1,j2,j3) = zfluxi (j1,j2,j3) * zfgas (j1,j2,j3)
              zfluxdi(j1,j2,j3) = zfluxdi(j1,j2,j3) * zfgasd(j1,j2,j3)
              zfluxui(j1,j2,j3) = zfluxui(j1,j2,j3) * zfgasu(j1,j2,j3)
            ENDDO
          ENDDO
        ENDDO
 
        IF (igasz.eq.igasm1) THEN    !     Avoid unphysical pseudo-transmission
          DO j3 = ki3sc, ki3ec+1
            DO j2 = ki2sc, ki2ec
              DO j1 = ki1sc, ki1ec
                zfluxi (j1,j2,j3) = MIN( 1.0_ireals,                           &
                                         MAX( 0.0_ireals, zfluxi (j1,j2,j3)) )
                zfluxdi(j1,j2,j3) = MIN( 1.0_ireals,                           &
                                         MAX( 0.0_ireals, zfluxdi(j1,j2,j3)) )
                zfluxui(j1,j2,j3) = MIN( 1.0_ireals,                           &
                                         MAX( 0.0_ireals, zfluxui(j1,j2,j3)) )
              ENDDO
            ENDDO
          ENDDO
        END IF

      END IF                   ! Test, whether gas needs to be included
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      END DO          ! End of loop over gases
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 
      ! Store FESFT result in zflux
      DO j3 = ki3sc, ki3ec+1
        DO j2 = ki2sc, ki2ec
          DO j1 = ki1sc, ki1ec
            zflux (j1,j2,j3) = zfluxi (j1,j2,j3) 
            zfluxd(j1,j2,j3) = zfluxdi(j1,j2,j3)
            zfluxu(j1,j2,j3) = zfluxui(j1,j2,j3)
          ENDDO
        ENDDO
      ENDDO
 
    !--------------------------------------------------------------
    END IF                             ! ESFT/FESFT-Selection 
    !--------------------------------------------------------------

  !----------------------------------------------------------------------------
  ! Section 3.3: Compute corrected fluxes, if lradtopo
  !----------------------------------------------------------------------------

    IF (lradtopo) THEN
      IF (idebug > 15) THEN
        WRITE (*,'(A60,2F16.6)')                                             &
                  'FESFT: solar fluxes before and after, skyview, fcor',     &
                          pskyview(j1b,j2b), pfcor(j1b,j2b)
        WRITE (*,'(A32,2I4,3F16.6)')  'FESFT: zfluxd, zflux, zfluxu',        &
                j1b, j2b, zfluxd(j1b,j2b,ki3ec+1), zflux(j1b,j2b,ki3ec+1),   &
                          zfluxu(j1b,j2b,ki3ec+1)
        zalbedo = zfluxu(j1b,j2b,ki3ec+1) /                                  &
                          (zflux(j1b,j2b,ki3ec+1)+zfluxd(j1b,j2b,ki3ec+1))
        WRITE (*,'(A32,2F16.6)') 'albedo, palso', zalbedo, palso(j1b,j2b)
      ENDIF

      DO j2 = ki2sc, ki2ec
        DO j1 = ki1sc, ki1ec

          zalbedo = zfluxu(j1,j2,ki3ec+1) /                                  &
                        (zflux(j1,j2,ki3ec+1)+zfluxd(j1,j2,ki3ec+1))

          ! direct down corrected
          zflux_c (j1,j2,ki3ec+1) = pfcor(j1,j2) * zfluxi(j1,j2,ki3ec+1)

          ! diffuse down corrected
          zfluxd_c(j1,j2,ki3ec+1) =                                          &
                    zfluxdi(j1,j2,ki3ec+1) *             pskyview(j1,j2)     &
                  + zfluxui(j1,j2,ki3ec+1) * (1.0_ireals-pskyview(j1,j2))

          ! diffuse up adapted to new other components
          zfluxu_c(j1,j2,ki3ec+1) =                                          &
                   (zflux_c(j1,j2,ki3ec+1) + zfluxd_c(j1,j2,ki3ec+1)) * zalbedo

        ENDDO
      ENDDO

      IF (idebug > 15) THEN
        WRITE (*,'(A42,2I4,3F16.6)')                                         &
                      'FESFT corrected: zfluxd, zflux, zfluxu', j1b, j2b,    &
                      zfluxd_c(j1b,j2b,ki3ec+1), zflux_c(j1b,j2b,ki3ec+1),   &
                      zfluxu_c(j1b,j2b,ki3ec+1)
      ENDIF
    ENDIF

  !----------------------------------------------------------------------------
  ! Section 3.4: Addition of flux for spectral interval to total solar flux
  !----------------------------------------------------------------------------
 
    IF (icrf == 0) THEN     ! Add flux at all levels

      DO j3 = ki3sc, ki3ec+1
        DO j2 = ki2sc, ki2ec
          DO j1 = ki1sc, ki1ec
            pfls(j1,j2,j3) = pfls (j1,j2,j3)                                 &
                  + zflux(j1,j2,j3) + zfluxd(j1,j2,j3) - zfluxu(j1,j2,j3)
            ! for the Climate-LM Version
            IF ( (.NOT. lradtopo) .OR. (pfcor(j1,j2) /= 0.0_ireals) ) THEN
              pflsdir(j1,j2,j3) = pflsdir (j1,j2,j3)     + zflux(j1,j2,j3)
            ! the else part just lets pflsdir untouched
            ENDIF
          ENDDO
        ENDDO
      ENDDO

      ! Store individual components of solar flux at surface
      IF (lradtopo) THEN
        DO j2 = ki2sc, ki2ec
          DO j1 = ki1sc, ki1ec
            pflsu (j1,j2) = pflsu (j1,j2) + zfluxu_c(j1,j2,ki3ec+1)
            pflsd (j1,j2) = pflsd (j1,j2) + zfluxd_c(j1,j2,ki3ec+1)
            pflsp (j1,j2) = pflsp (j1,j2) + zflux_c (j1,j2,ki3ec+1)
            pfls_s(j1,j2) = pfls_s(j1,j2) + zflux_c (j1,j2,ki3ec+1)          &
                      + zfluxd_c(j1,j2,ki3ec+1) - zfluxu_c(j1,j2,ki3ec+1)
          ENDDO
        ENDDO
      ELSE
        DO j2 = ki2sc, ki2ec
          DO j1 = ki1sc, ki1ec
            pflsu (j1,j2) = pflsu (j1,j2) + zfluxu  (j1,j2,ki3ec+1)
            pflsd (j1,j2) = pflsd (j1,j2) + zfluxd  (j1,j2,ki3ec+1)
            pflsp (j1,j2) = pflsp (j1,j2) + zflux   (j1,j2,ki3ec+1)
            pfls_s(j1,j2) = pfls_s(j1,j2) + zflux   (j1,j2,ki3ec+1)          &
                      + zfluxd  (j1,j2,ki3ec+1) - zfluxu  (j1,j2,ki3ec+1)
          ENDDO
        ENDDO
      ENDIF
      IF (idebug > 15) THEN
        WRITE (*,'(A40,3F16.6)') 'FESFT: diff_up, diff_down, dir_down',      &
            pflsu(j1b,j2b), pflsd(j1b,j2b), pflsp(j1b,j2b)
      ENDIF

      IF (jspec == 3) THEN   ! Photosynthetic active radiation
        IF (lradtopo) THEN   ! T.R.
          DO j2 = ki2sc, ki2ec
            DO j1 = ki1sc, ki1ec
              pflsu_par(j1,j2) = pflsu_par (j1,j2) + zfluxu_c(j1,j2,ki3ec+1)
              pflsd_par(j1,j2) = pflsd_par (j1,j2) + zfluxd_c(j1,j2,ki3ec+1)
              pflsp_par(j1,j2) = pflsp_par (j1,j2) + zflux_c (j1,j2,ki3ec+1)
              pflpar   (j1,j2) = pflpar    (j1,j2) + zflux_c (j1,j2,ki3ec+1) &
                      + zfluxd_c(j1,j2,ki3ec+1) - zfluxu_c(j1,j2,ki3ec+1)
            ENDDO
          ENDDO
        ELSE
          DO j2 = ki2sc, ki2ec
            DO j1 = ki1sc, ki1ec
              pflsu_par(j1,j2) = pflsu_par (j1,j2) + zfluxu(j1,j2,ki3ec+1)
              pflsd_par(j1,j2) = pflsd_par (j1,j2) + zfluxd(j1,j2,ki3ec+1)
              pflsp_par(j1,j2) = pflsp_par (j1,j2) + zflux (j1,j2,ki3ec+1)
              pflpar   (j1,j2) = pflpar    (j1,j2) + zflux (j1,j2,ki3ec+1)   &
                      + zfluxd  (j1,j2,ki3ec+1) - zfluxu  (j1,j2,ki3ec+1)
            ENDDO
          ENDDO
        ENDIF
      ENDIF !(jspec == 3)

    ENDIF
 
#ifdef COSMOART
  IF(l_cosmo_art) THEN
    IF (jspec == 3) THEN
      IF (lradtopo) THEN
        DO j3 = ki3sc, ki3ec+1
          DO j2 = ki2sc, ki2ec
            DO j1 = ki1sc, ki1ec
              Edir (j1,jindex,j3) =  zflux_c(j1,j2,j3)
              Edown(j1,jindex,j3) =  zfluxd_c(j1,j2,j3)
              Eup  (j1,jindex,j3) =  zfluxu_c(j1,j2,j3)
            ENDDO
          ENDDO
        ENDDO
      ELSE
        DO j3 = ki3sc, ki3ec+1
          DO j2 = ki2sc, ki2ec
            DO j1 = ki1sc, ki1ec
              Edir (j1,jindex,j3) =  zflux(j1,j2,j3)
              Edown(j1,jindex,j3) =  zfluxd(j1,j2,j3)
              Eup  (j1,jindex,j3) =  zfluxu(j1,j2,j3)
            ENDDO
          ENDDO
        ENDDO
      ENDIF
    ENDIF
  ENDIF
#endif

  ! End of solar spectral loop
  ! ===============================================================
  ENDDO solar_spectral_loop 
  ! ===============================================================
 
  ENDIF     ! Test, whether solar calculation or not
 
!------------------------------------------------------------------------------
! Section 4: Repeat calculations for cloud-free fluxes
!------------------------------------------------------------------------------

  ! Repeat calculations for cloud-free fluxes if switch for CRF
  ! is set to .true. and cloud-free fluxes have not yet been
  ! computed

  IF (.NOT. lcrf) THEN
!T.R.: pfltf & pflsf removed
!    DO j2 = ki2sc, ki2ec
!      DO j1 = ki1sc, ki1ec
!        pflsf(j1,j2,1)  = pfls(j1,j2,ki3sc  )
!        pflsf(j1,j2,2)  = pfls(j1,j2,ki3ec+1)
!        pfltf(j1,j2,1)  = pflt(j1,j2,ki3sc  )
!        pfltf(j1,j2,2)  = pflt(j1,j2,ki3ec+1)
!      ENDDO
!    ENDDO
  ELSE IF (icrf.eq.0) THEN  ! Branch to cloud-free calculations only once
    icrf = 1     
    GO TO 1
  ENDIF
 
!------------------------------------------------------------------------------
! End of the subroutine
!------------------------------------------------------------------------------

END SUBROUTINE fesft 

!==============================================================================
!==============================================================================
!+ Module procedure in "Radiation"
!------------------------------------------------------------------------------

SUBROUTINE coe_th (                                                    &
       pduh2oc,pduh2of,pduco2 ,pduo3  ,palogp ,palogt ,                &
       podsc  ,podsf  ,podac  ,podaf  ,pbsfc  ,pbsff  ,                &
       ki3    ,kspec  ,kh2o   ,kco2   ,ko3    ,                        &
       ki1sd  ,ki1ed  ,ki2sd  ,ki2ed  ,ki3sd  ,ki3ed,                  &
       ki1sc  ,ki1ec  ,ki2sc  ,ki2ec  ,ki3sc  ,ki3ec,                  &
       ldebug ,                                                        &
       pa1c   ,pa1f   ,pa2c   ,pa2f   ,pa3c   ,pa3f)

!------------------------------------------------------------------------------
!
! Description:
!
!   The module procedure coe_th calculates the optical effects of atmospheric 
!   layers on thermal radiation based on basic optical properties of non-gaseous 
!   constituents and gaseous absorption coefficients selected through the 
!   corresponding control variables in the argument list.
!   This routine computes layer effects (transmissivity, reflectivity
!   and emmisivity) in the thermal part of the radiative spectrum
!   both for the cloud-free and the cloudy part of a model layer.
!   The calculation is based on the implicit delt-two-stream equations
!   (cf. Ritter and Geleyn, 1992) and uses basic optical properties
!   (i.e. absorption and scattering optical depth and backscattered
!   fraction for non-gaseous atmospheric constituents as well as 
!   gaseous absorption properties) as input. 
!
! Method:
!
! - addition of individual gaseous absorption effects to the optical
!   properties of the non-gaseous constituents
! - determination of layer effects (cf. Zdunkowski et al., 1982, 1986
!   and Ritter and Geleyn, 1992)
!     
!------------------------------------------------------------------------------

! Subroutine arguments:
! --------------------

! Input data
! ----------
  INTEGER (KIND=iintegers), INTENT (IN) ::  &
     ki1sd,       & ! start index for first  array dimension
     ki1ed,       & ! end   index for first  array dimension
     ki2sd,       & ! start index for second array dimension
     ki2ed,       & ! end   index for second array dimension
     ki3sd,       & ! start index for third  array dimension
     ki3ed,       & ! end   index for third  array dimension

   ! and the same for the computations
     ki1sc,       & ! start index for first  array computation
     ki1ec,       & ! end   index for first  array computation
     ki2sc,       & ! start index for second array computation
     ki2ec,       & ! end   index for second array computation
     ki3sc,       & ! start index for third  array computation
     ki3ec,       & ! end   index for third  array computation
     ki3  ,       & ! vertical layer considered    
     kspec,       & ! spectral interval considered
     kh2o ,       & ! table index for h2o absorption properties
     kco2 ,       & ! table index for co2 absorption properties
     ko3            ! table index for o3  absorption properties

  LOGICAL                 , INTENT (IN) ::  &
     ldebug         ! debug control switch       

  REAL    (KIND=ireals   ), INTENT (IN) ::  &

     ! opticall relevant gas quantities (Pa)
     pduh2oc(ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed), & ! h2o inside cloud
     pduh2of(ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed), & ! h2o out of cloud
     pduco2 (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed), & ! co2 content 
     pduo3  (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed), & ! o3  content 

     ! Logarithm of layer mean temperature and pressure
     palogt (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed), & ! ln T
     palogp (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed), & ! ln p

     ! Optical properties of non-gaseous constituents (..c=cloudy; ..f=free)  
     podsc  (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed), & ! 
     podsf  (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed), & ! 
     podac  (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed), & ! 
     podaf  (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed), & ! 
     pbsfc  (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed), & ! 
     pbsff  (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed)    ! 

! Output data
! -----------
  REAL    (KIND=ireals   ), INTENT (OUT) ::  &
     pa1c  (ki1sd:ki1ed,ki2sd:ki2ed), & ! transmissivity in cloud   
     pa1f  (ki1sd:ki1ed,ki2sd:ki2ed), & ! transmissivity cloud-free  
     pa2c  (ki1sd:ki1ed,ki2sd:ki2ed), & ! reflectivity in cloud    
     pa2f  (ki1sd:ki1ed,ki2sd:ki2ed), & ! reflectivity cloud-free      
     pa3c  (ki1sd:ki1ed,ki2sd:ki2ed), & ! emissivity in cloud    
     pa3f  (ki1sd:ki1ed,ki2sd:ki2ed)    ! emissivity cloud-free       

! Local parameters: 
! ----------------
  REAL    (KIND=ireals   ), PARAMETER ::  &
     zargli  = 80.0     , &  ! argument limit for EXP 
     ztsec   = 1.0E-35  , &  ! (=exp(-zargli) avoids ALOG(0.0)
     zodmax  = 1.0E+6   , &  ! maximum allowed optical depth
     zudiff  = 2.0      , &  ! Diffusivity factors for gases and other constituents
     zangfa  = 1.648721271   ! exp(0.5)

  INTEGER (KIND=iintegers), PARAMETER ::  &
     j1b    = 1,           & ! debug point index first dimension
     j2b    = 1              ! debug point index second dimension

! Local scalars:
! -------------
  INTEGER (KIND=iintegers) ::  &
    j1,j2,j3                 ! loop indices over spatial dimensions

  REAL    (KIND=ireals   ) ::  &
    zeps, ztau, zrho, zodgf, zodgc, zod1, zod2
 
! End of header
!==============================================================================
 
!------------------------------------------------------------------------------
! Begin Subroutine coe_th              
!------------------------------------------------------------------------------

  j3     = ki3

  IF (ldebug) THEN
     print *,'**** coe_th ******************************'
     print *,'**** debug point : ',j1b,j2b
     print *,'**** coe_th kspec=',kspec
     print *,'**** coe_th j3   =',j3   
     print *,'**** coe_th kh2o =',kh2o 
     print *,'**** coe_th kco2 =',kco2 
     print *,'**** coe_th ko3  =',ko3  
     print *,'**** pduh2of(j1b,j2b,j3)=',pduh2of(ki1sc,j2b,j3)
     print *,'**** pduh2oc(j1b,j2b,j3)=',pduh2oc(ki1sc,j2b,j3)
     print *,'**** pduco2 (j1b,j2b,j3)=',pduco2 (ki1sc,j2b,j3)
     print *,'**** pduo3  (j1b,j2b,j3)=',pduo3  (ki1sc,j2b,j3)
     print *,'**** palogp (j1b,j2b,j3)=',palogp (ki1sc,j2b,j3)
     print *,'**** palogt (j1b,j2b,j3)=',palogt (ki1sc,j2b,j3)
     print *,'**** cobi (kh2o,kspec,1)    =',cobi (kh2o,kspec,1)    
     print *,'**** cobi (kco2,kspec,2)    =',cobi (kco2,kspec,2)    
     print *,'**** cobi (ko3 ,kspec,3)    =',cobi (ko3 ,kspec,3)    
     print *,'**** coali(kh2o,kspec,1)    =',coali (kh2o,kspec,1)    
     print *,'**** coali(kco2,kspec,2)    =',coali (kco2,kspec,2)    
     print *,'**** coali(ko3 ,kspec,3)    =',coali (ko3 ,kspec,3)    
     print *,'**** cobti(kh2o,kspec,1)    =',cobti(kh2o,kspec,1)    
     print *,'**** cobti(kco2,kspec,2)    =',cobti(kco2,kspec,2)    
     print *,'**** cobti(ko3 ,kspec,3)    =',cobti(ko3 ,kspec,3)    
  ENDIF 

  ! Optical depth of gases

  DO j2 = ki2sc, ki2ec
    DO j1 = ki1sc, ki1ec
      zodgf = 0.0     ! Initialisation
 
      IF (kco2.ne.0) then     ! Include CO2 contribution
        zodgf = zodgf + pduco2(j1,j2,j3) * (cobi(kco2,kspec,2) &
                    * EXP ( coali(kco2,kspec,2) * palogp(j1,j2,j3)    &
                           -cobti(kco2,kspec,2) * palogt(j1,j2,j3)))
      ENDIF                  ! CO2
      !US IF (ldebug) print *,'**** zodgf(CO2)        =',zodgf
 
      IF (ko3 /= 0) THEN  ! Include O3 contribution
        zodgf = zodgf + pduo3 (j1,j2,j3) * (cobi(ko3 ,kspec,3)* &
                    EXP ( coali(ko3 ,kspec,3) * palogp(j1,j2,j3)      &
                         -cobti(ko3 ,kspec,3) * palogt(j1,j2,j3)))
      ENDIF
      !US IF (ldebug) print *,'**** zodgf(CO2+O3)     =',zodgf
 
      ! Cloudy = cloud free for CO2 and O3 :
      zodgc = zodgf
 
      IF (kh2o /= 0) THEN  ! Include H2O contribution
        zodgf = zodgf + pduh2of(j1,j2,j3)* (cobi(kh2o,kspec,1)* &
                    EXP ( coali(kh2o,kspec,1) * palogp(j1,j2,j3)      &
                         -cobti(kh2o,kspec,1) * palogt(j1,j2,j3)))
        zodgc = zodgc + pduh2oc(j1,j2,j3)* (cobi(kh2o,kspec,1)* &
                    EXP ( coali(kh2o,kspec,1) * palogp(j1,j2,j3)      &
                         -cobti(kh2o,kspec,1) * palogt(j1,j2,j3)))
      ENDIF
!------------------------------------------------------------------------------
      !US IF (ldebug) print *,'**** zodgf(CO2+O3+H2O) =',zodgf
      !US IF (ldebug) print *,'**** zodgc(CO2+O3+H2O) =',zodgc
 
      zodgf = MIN (zodgf, zodmax)
      zodgc = MIN (zodgc, zodmax)

      !US IF (ldebug) print *,'**** nach securit auf optical depth '
      !US IF (ldebug) print *,'**** zodgf(CO2+O3+H2O) =',zodgf
      !US IF (ldebug) print *,'**** zodgc(CO2+O3+H2O) =',zodgc
 
      ! Pseudo-optical depth in cloud-free part of layer

      zod2 = zudiff * pbsff(j1,j2,j3) * podsf(j1,j2,j3)
      zod1 = zod2 + zudiff * podaf(j1,j2,j3)
      zod1 = zod1 + zangfa * zodgf

      !US IF (ldebug) THEN
      !US   print *,'**** cloud-free zod1 (j1b,j2b)=',zod1
      !US   print *,'**** cloud-free zod2 (j1b,j2b)=',zod2
      !US ENDIF 
 
      ! Layer coefficients in cloud-free part of layer
 
      zeps=SQRT(zod1*zod1-zod2*zod2)
      IF (zeps.LT.zargli) THEN
        ztau = EXP  (-zeps)
      ELSE
        ztau = ztsec
      END IF
      zrho = zod2/(zod1+zeps)
      pa1f(j1,j2)=ztau*(1.-(zrho**2))*(1./(1.-(zrho**2)*(ztau**2)))
      pa2f(j1,j2)=zrho*(1.-(ztau**2))*(1./(1.-(zrho**2)*(ztau**2)))
      pa3f(j1,j2)=(1.-pa1f(j1,j2)+pa2f(j1,j2))/(zod1+zod2)

      !US IF (ldebug) THEN
      !US   print *,'**** cloud-free pa1f (j1b,j2b)=',pa1f (j1b,j2b)
      !US   print *,'**** cloud-free pa2f (j1b,j2b)=',pa2f (j1b,j2b)
      !US   print *,'**** cloud-free pa3f (j1b,j2b)=',pa3f (j1b,j2b)
      !US ENDIF
 
      ! Pseudo-optical depth in cloudy part of layer
      zod2 = zudiff * pbsfc(j1,j2,j3) * podsc(j1,j2,j3)
      zod1 = zod2 + zudiff * podac(j1,j2,j3)
      zod1 = zod1 + zangfa * zodgc
 
      ! Layer coefficients in cloudy part of layer
 
      zeps=SQRT(zod1*zod1-zod2*zod2)
      IF (zeps.LT.zargli) THEN
        ztau = EXP  (-zeps)
      ELSE
        ztau = ztsec
      END IF
      zrho = zod2/(zod1+zeps)
      pa1c(j1,j2)=ztau*(1.-(zrho**2))*(1./(1.-(zrho**2)*(ztau**2)))
      pa2c(j1,j2)=zrho*(1.-(ztau**2))*(1./(1.-(zrho**2)*(ztau**2)))
      pa3c(j1,j2)=(1.-pa1c(j1,j2)+pa2c(j1,j2))/(zod1+zod2)
    ENDDO
  ENDDO
 
!------------------------------------------------------------------------------
! End of the subroutine
!------------------------------------------------------------------------------

END SUBROUTINE coe_th

!==============================================================================
!==============================================================================
!+ Module procedure in "Radiation"
!------------------------------------------------------------------------------

SUBROUTINE coe_so (                                                     &
       pduh2oc,pduh2of,pduco2 ,pduo3  ,palogp ,palogt ,                 &
       podsc  ,podsf  ,podac  ,podaf  ,pbsfc  ,pbsff  ,pusfc ,pusff ,   &
       psmu0  ,pqsmu0 ,                                                 &
       ki3    ,kspec  ,kh2o   ,kco2   ,ko3    ,                         &
       ki1sd  ,ki1ed  ,ki2sd  ,ki2ed  ,ki3sd  ,ki3ed,                   &
       ki1sc  ,ki1ec  ,ki2sc  ,ki2ec  ,ki3sc  ,ki3ec,                   &
       ldebug ,                                                         &
       pa1c   ,pa1f   ,pa2c   ,pa2f   ,pa3c   ,pa3f ,                   &
       pa4c   ,pa4f   ,pa5c   ,pa5f )

!------------------------------------------------------------------------------
!
! Description:
!
!   The module procedure coe_so calculates the optical effects of atmospheric
!   layers on solar radiation based on basic optical properties of non-gaseous
!   constituents and gaseous absorption coefficients selected through the
!   corresponding control variables.
!   This routine computes layer effects (transmissivity, reflectivity)
!   for diffuse and direct solar radiation both for the cloudy and the
!   cloud-free part of a model layer.
!   The calculation is based on the implicit delt-two-stream equations
!   (cf. Ritter and Geleyn, 1992) and uses basic optical properties
!   (i.e. absorption and scattering optical depth, backscattered and
!   upscattered fraction for non-gaseous atmospheric constituents and
!   gaseous absorption properties) as input.
!
! Method:
!
! - addition of individual gaseous absorption effects to the optical
!   properties of the non-gaseous constituents
!   (optical depth multiplied by alpha1 to alpha4)
!
! - determination of layer effects (cf. Zdunkowski et al., 1982, 1986
!   and Ritter and Geleyn, 1992)
! - the resonance case for those effects related to the direct solar
!   radiation is avoided by a small displacement of the local inverse
!   of the cosine of the zenith angle (if necessary)
!
!------------------------------------------------------------------------------

! Subroutine arguments:
! --------------------

! Input data
! ----------
  INTEGER (KIND=iintegers), INTENT (IN) ::  &
     ki1sd,       & ! start index for first  array dimension
     ki1ed,       & ! end   index for first  array dimension
     ki2sd,       & ! start index for second array dimension
     ki2ed,       & ! end   index for second array dimension
     ki3sd,       & ! start index for third  array dimension
     ki3ed,       & ! end   index for third  array dimension

   ! and the same for the computations
     ki1sc,       & ! start index for first  array computation
     ki1ec,       & ! end   index for first  array computation
     ki2sc,       & ! start index for second array computation
     ki2ec,       & ! end   index for second array computation
     ki3sc,       & ! start index for third  array computation
     ki3ec,       & ! end   index for third  array computation
     ki3  ,       & ! vertical layer considered
     kspec,       & ! spectral interval considered
     kh2o ,       & ! table index for h2o absorption properties
     kco2 ,       & ! table index for co2 absorption properties
     ko3            ! table index for o3  absorption properties

  LOGICAL                 , INTENT (IN) ::  &
     ldebug         ! debug control switch

  REAL    (KIND=ireals   ), INTENT (IN) ::  &

     ! opticall relevant gas quantities (Pa)
     pduh2oc(ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed), & ! h2o inside cloud
     pduh2of(ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed), & ! h2o out of cloud
     pduco2 (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed), & ! co2 content
     pduo3  (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed), & ! o3  content

     ! Logarithm of layer mean temperature and pressure
     palogt (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed), & ! ln T
     palogp (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed), & ! ln p

     ! Optical properties of non-gaseous constituents (..c=cloudy; ..f=free)
     podsc  (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed), & !
     podsf  (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed), & !
     podac  (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed), & !
     podaf  (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed), & !
     pbsfc  (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed), & !
     pbsff  (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed), & !
     pusfc  (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed), & !
     pusff  (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed), & !

     psmu0  (ki1sd:ki1ed,ki2sd:ki2ed)            , & ! cosine of zenith angle
     pqsmu0 (ki1sd:ki1ed,ki2sd:ki2ed)                ! inverse of cosine ...

! Output data
! -----------
  REAL    (KIND=ireals   ), INTENT (OUT) ::  &
     pa1c  (ki1sd:ki1ed,ki2sd:ki2ed), & ! direct radiation transmis-
     pa1f  (ki1sd:ki1ed,ki2sd:ki2ed), & ! sivity cloudy/cloud-free
     pa2c  (ki1sd:ki1ed,ki2sd:ki2ed), & ! direct radition downward
     pa2f  (ki1sd:ki1ed,ki2sd:ki2ed), & ! scattering cloudy/cloud-free
     pa3c  (ki1sd:ki1ed,ki2sd:ki2ed), & ! direct radiation upward
     pa3f  (ki1sd:ki1ed,ki2sd:ki2ed), & ! scattering cloudy/cloud-free
     pa4c  (ki1sd:ki1ed,ki2sd:ki2ed), & ! diffuse flux transmissivity
     pa4f  (ki1sd:ki1ed,ki2sd:ki2ed), & ! cloudy/cloud-free
     pa5c  (ki1sd:ki1ed,ki2sd:ki2ed), & ! diffuse flux reflectivity
     pa5f  (ki1sd:ki1ed,ki2sd:ki2ed)    ! cloudy/cloud-free

! Local parameters:
! ----------------
  REAL    (KIND=ireals   ), PARAMETER ::  &
     zargli  = 80.0     , &  ! argument limit for EXP
     ztsec   = 1.0E-35  , &  ! (=exp(-zargli) avoids ALOG(0.0)
     zodmax  = 1.0E+6   , &  ! maximum allowed optical depth
     zepres  = 1.0E-7   , &  ! for resonance case avoidance
                             ! 32bit-accuracy (1.E-14 for 64bit machine)
     zudiff  = 2.0      , &  ! Diffusivity factors for gases and other constituents
     zangfa  = 1.648721271   ! exp(0.5)

  INTEGER (KIND=iintegers), PARAMETER ::  &
     j1b    = 1,           & ! debug point index first dimension
     j2b    = 1              ! debug point index second dimension

! Local scalars:
! -------------
  INTEGER (KIND=iintegers) ::  &
    j1,j2,j3                 ! loop indices over spatial dimensions

  REAL    (KIND=ireals   ) ::  &
    zeps,                      & !
    ze,zm,zg1,zg2,ze1mwf,zmu0if  !

  REAL    (KIND=ireals   ) ::  &
     zodgf, zodgc, zod1, zod2, zod3, zod4, zod5
 
!- End of header
!==============================================================================

!------------------------------------------------------------------------------
! Begin Subroutine coe_so
!------------------------------------------------------------------------------

  j3     = ki3

  IF (ldebug) THEN
     print *,'**** coe_so ******************************'
     print *,'**** debug point index : ',j1b,j2b
     print *,'**** coe_so kspec=',kspec
     print *,'**** coe_so j3   =',j3
     print *,'**** coe_so kh2o =',kh2o
     print *,'**** coe_so kco2 =',kco2
     print *,'**** coe_so ko3  =',ko3
     print *,'**** pduh2of(j1b,j2b,j3)=',pduh2of(j1b,j2b,j3)
     print *,'**** pduh2oc(j1b,j2b,j3)=',pduh2oc(j1b,j2b,j3)
     print *,'**** pduco2 (j1b,j2b,j3)=',pduco2 (j1b,j2b,j3)
     print *,'**** pduo3  (j1b,j2b,j3)=',pduo3  (j1b,j2b,j3)
     print *,'**** psmu0  (j1b,j2b)   =',psmu0  (j1b,j2b)
  ENDIF

  ! Optical depth of gases

  DO j2 = ki2sc , ki2ec
    DO j1 = ki1sc, ki1ec
      zodgf = 0.0_ireals        ! Initialisation

      IF (kco2 /= 0) THEN     ! Include CO2 contribution
        zodgf = zodgf + pduco2(j1,j2,j3)* (cobi(kco2,kspec,2)*       &
                    EXP ( coali(kco2,kspec,2) * palogp(j1,j2,j3)     &
                         -cobti(kco2,kspec,2) * palogt(j1,j2,j3)))
      ENDIF                  ! CO2
      !US IF (ldebug) print *,'**** zodgf(CO2)        =',zodgf

      IF (ko3 /= 0) THEN     ! Include O3 contribution
        zodgf = zodgf + pduo3 (j1,j2,j3)* (cobi(ko3 ,kspec,3)*       &
                    EXP ( coali(ko3 ,kspec,3) * palogp(j1,j2,j3)     &
                         -cobti(ko3 ,kspec,3) * palogt(j1,j2,j3)))
      ENDIF
      !US IF (ldebug) print *,'**** zodgf(CO2+O3)     =',zodgf

      ! Cloudy = cloud free for CO2 and O3 :
      zodgc = zodgf

      IF (kh2o /= 0) THEN    ! Include H2O contribution
        zodgf = zodgf + pduh2of(j1,j2,j3)* (cobi(kh2o,kspec,1)*       &
                    EXP ( coali(kh2o,kspec,1) * palogp(j1,j2,j3)      &
                         -cobti(kh2o,kspec,1) * palogt(j1,j2,j3)))
        zodgc = zodgc + pduh2oc(j1,j2,j3)* (cobi(kh2o,kspec,1)*       &
                    EXP ( coali(kh2o,kspec,1) * palogp(j1,j2,j3)      &
                         -cobti(kh2o,kspec,1) * palogt(j1,j2,j3)))
      ENDIF
      !US IF (ldebug) print *,'**** zodgf(CO2+O3+H2O) =',zodgf
      !US IF (ldebug) print *,'**** zodgc(CO2+O3+H2O) =',zodgc

      zodgf = MIN (zodgf, zodmax)
      zodgc = MIN (zodgc, zodmax)
      !US IF (ldebug) print *,'**** nach securit auf optical depth '
      !US IF (ldebug) print *,'**** zodgf(CO2+O3+H2O) =',zodgf
      !US IF (ldebug) print *,'**** zodgc(CO2+O3+H2O) =',zodgc

      ! Pseudo-optical depth in cloud-free part of layer

      zod2 = zudiff * pbsff(j1,j2,j3) * podsf(j1,j2,j3)
      zod1 = zod2 + zudiff * podaf(j1,j2,j3)
      zod3 = pusff(j1,j2,j3) * podsf(j1,j2,j3)
      zod4 = podsf(j1,j2,j3) - zod3
      zod5 = podsf(j1,j2,j3) + podaf(j1,j2,j3)
      zod1 = zod1 + zangfa * zodgf
      zod5 = zod5 + zodgf
      !US IF (ldebug) THEN
      !US   print *,'**** cloud-free zod1 (j1b,j2b)    =',zod1
      !US   print *,'**** cloud-free zod2 (j1b,j2b)    =',zod2
      !US   print *,'**** cloud-free zod3 (j1b,j2b)    =',zod3
      !US   print *,'**** cloud-free zod4 (j1b,j2b)    =',zod4
      !US   print *,'**** cloud-free zod5 (j1b,j2b)    =',zod5
      !US ENDIF

      ! Layer coefficients in cloud-free part of layer

      zeps=SQRT(zod1*zod1-zod2*zod2)
      IF (zeps.LT.zargli) THEN
        ze = EXP  (-zeps)
      ELSE
        ze = ztsec
      END IF
      zm = zod2/(zod1+zeps)
      pa4f(j1,j2)=ze*(1.-(zm**2))*(1./(1.-(zm**2)*(ze**2)))
      pa5f(j1,j2)=zm*(1.-(ze**2))*(1./(1.-(zm**2)*(ze**2)))

      ze1mwf = zeps / zod5
      zmu0if = ze1mwf + SIGN ( MAX(ABS(pqsmu0(j1,j2)-ze1mwf),zepres) &
                              ,(pqsmu0(j1,j2)-ze1mwf) )
      zod3 = zod3 * zmu0if
      zod4 = zod4 * zmu0if
      zod5 = zod5 * zmu0if
      IF (zod5.LT.zargli) THEN
        pa1f(j1,j2) = EXP  (-zod5)
      ELSE
        pa1f(j1,j2) = ztsec
      END IF
      zg1 = ( zod3*(zod5-zod1) -zod2*zod4) /(zod5*zod5 - zeps*zeps)
      zg2 =-( zod4*(zod5+zod1) +zod2*zod3) /(zod5*zod5 - zeps*zeps)
      pa2f(j1,j2) = zg2*(pa1f(j1,j2)-pa4f(j1,j2)) -zg1*pa5f(j1,j2)*pa1f(j1,j2)
      pa3f(j1,j2) = zg1*(1.-pa4f(j1,j2)*pa1f(j1,j2)) -zg2*pa5f(j1,j2)

      !US IF (ldebug) THEN
      !US   print *,'**** cloud-free pa1f (j1b,j2b)    =',pa1f (j1b,j2b)
      !US   print *,'**** cloud-free pa2f (j1b,j2b)    =',pa2f (j1b,j2b)
      !US   print *,'**** cloud-free pa3f (j1b,j2b)    =',pa3f (j1b,j2b)
      !US   print *,'**** cloud-free pa4f (j1b,j2b)    =',pa4f (j1b,j2b)
      !US   print *,'**** cloud-free pa5f (j1b,j2b)    =',pa5f (j1b,j2b)
      !US ENDIF

      ! Pseudo-optical depth in cloudy part of layer

      zod2 = zudiff * pbsfc(j1,j2,j3) * podsc(j1,j2,j3)
      zod1 = zod2 + zudiff * podac(j1,j2,j3)
      zod3 = pusfc(j1,j2,j3) * podsc(j1,j2,j3)
      zod4 = podsc(j1,j2,j3) - zod3
      zod5 = podsc(j1,j2,j3) + podac(j1,j2,j3)
      zod1 = zod1 + zangfa * zodgc
      zod5 = zod5 + zodgc

      !US IF (ldebug) THEN
      !US   print *,'**** cloudy     zod1 (j1b,j2b)    =',zod1
      !US   print *,'**** cloudy     zod2 (j1b,j2b)    =',zod2
      !US   print *,'**** cloudy     zod3 (j1b,j2b)    =',zod3
      !US   print *,'**** cloudy     zod4 (j1b,j2b)    =',zod4
      !US   print *,'**** cloudy     zod5 (j1b,j2b)    =',zod5
      !US ENDIF

      ! Layer coefficients in cloudy part of layer

      zeps=SQRT(zod1*zod1-zod2*zod2)
      IF (zeps.LT.zargli) THEN
        ze = EXP  (-zeps)
      ELSE
        ze = ztsec
      END IF
      zm = zod2/(zod1+zeps)
      pa4c(j1,j2)=ze*(1.-(zm**2))*(1./(1.-(zm**2)*(ze**2)))
      pa5c(j1,j2)=zm*(1.-(ze**2))*(1./(1.-(zm**2)*(ze**2)))

      ze1mwf = zeps / zod5
      zmu0if = ze1mwf + SIGN ( MAX(ABS(pqsmu0(j1,j2)-ze1mwf),zepres) &
                              ,(pqsmu0(j1,j2)-ze1mwf) )
      zod3 = zod3 * zmu0if
      zod4 = zod4 * zmu0if
      zod5 = zod5 * zmu0if
      IF (zod5.LT.ZARGLI) THEN
         pa1c(j1,j2) = EXP  (-zod5)
      ELSE
         pa1c(j1,j2) = ztsec
      END IF
      zg1 = ( zod3*(zod5-zod1) -zod2*zod4) /(zod5*zod5 - zeps*zeps)
      zg2 =-( zod4*(zod5+zod1) +zod2*zod3) /(zod5*zod5 - zeps*zeps)
      pa2c(j1,j2) = zg2*(pa1c(j1,j2)-pa4c(j1,j2)) -zg1*pa5c(j1,j2)*pa1c(j1,j2)
      pa3c(j1,j2) = zg1*(1.-pa4c(j1,j2)*pa1c(j1,j2)) -zg2*pa5c(j1,j2)
    ENDDO
  ENDDO

  IF (ldebug) THEN
      print *,'**** cloudy     pa1c (j1b,j2b)    =',pa1c (j1b,j2b)
      print *,'**** cloudy     pa2c (j1b,j2b)    =',pa2c (j1b,j2b)
      print *,'**** cloudy     pa3c (j1b,j2b)    =',pa3c (j1b,j2b)
      print *,'**** cloudy     pa4c (j1b,j2b)    =',pa4c (j1b,j2b)
      print *,'**** cloudy     pa5c (j1b,j2b)    =',pa5c (j1b,j2b)
  ENDIF

!-------------------------------------------------------------------------------
! End of the subroutine
!-------------------------------------------------------------------------------

END SUBROUTINE coe_so

!==============================================================================
!==============================================================================
!+ Module procedure in "Radiation" 
!------------------------------------------------------------------------------

SUBROUTINE inv_th (                                                   &
       pclc   ,pca1   ,pca2  ,pcb1  ,pcb2  ,pcc1  ,pcc2 ,pcd1 ,pcd2 , &
       pduh2oc,pduh2of,pduco2,pduo3 ,palogp,palogt,                   &
       podsc  ,podsf  ,podac ,podaf ,pbsfc ,pbsff ,                   &
       pbbr   ,palth,                                                 &
       kspec  ,kh2o   ,kco2  ,ko3   ,                                 &
       ki1sd  ,ki1ed  ,ki2sd ,ki2ed ,ki3sd ,ki3ed,                    &
       ki1sc  ,ki1ec  ,ki2sc ,ki2ec ,ki3sc ,ki3ec , ldebug ,          &
       pflcu  ,pflfu  ,pflcd ,pflfd)

!------------------------------------------------------------------------------
!
! Description:
!
!   The module procedure inv_th solves a linear equation system for thermal 
!   fluxes using a Gaussian elimination-backsubstitution algorithm dedicated 
!   to the specific structure of the system matrix.
!
! Method:
!
! - setting of the RHS of the system using the layer boundary black
!   body radiation and allowing for partial cloud cover in each layer
! - solution of the equation system including the lower boundary
!   condition
! - matrix coefficients are calculated in the course of the elimination
!   step for one layer at a time through a call to routine *coe_th*
! - the final result, i.e. the so-called black body flux differences
!   (cf.Ritter and Geleyn, 1992) are stored seperately for cloudy and
!   cloud-free part of each layer boundary
!
!------------------------------------------------------------------------------

! Subroutine arguments:
! --------------------

! Input data
! ----------
  INTEGER (KIND=iintegers), INTENT (IN) ::  &
     ki1sd,       & ! start index for first  array dimension
     ki1ed,       & ! end   index for first  array dimension
     ki2sd,       & ! start index for second array dimension
     ki2ed,       & ! end   index for second array dimension
     ki3sd,       & ! start index for third  array dimension
     ki3ed,       & ! end   index for third  array dimension

   ! and the same for the computations
     ki1sc,       & ! start index for first  array computation
     ki1ec,       & ! end   index for first  array computation
     ki2sc,       & ! start index for second array computation
     ki2ec,       & ! end   index for second array computation
     ki3sc,       & ! start index for third  array computation
     ki3ec,       & ! end   index for third  array computation
     kspec,       & ! spectral interval considered
     kh2o ,       & ! table index for h2o absorption properties
     kco2 ,       & ! table index for co2 absorption properties
     ko3            ! table index for o3  absorption properties

  LOGICAL                 , INTENT (IN) ::  &
     ldebug         ! debug control switch       

  REAL    (KIND=ireals   ), INTENT (IN) ::  &
     pclc  (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed),   & ! cloud cover
     pca1  (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed),   & ! cloud geometry factor  
     pca2  (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed),   & ! cloud geometry factor  
     pcb1  (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed),   & ! cloud geometry factor  
     pcb2  (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed),   & ! cloud geometry factor  
     pcc1  (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed),   & ! cloud geometry factor  
     pcc2  (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed),   & ! cloud geometry factor  
     pcd1  (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed),   & ! cloud geometry factor  
     pcd2  (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed),   & ! cloud geometry factor  
     pbbr  (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed+1), & ! black body radiation   
     palth (ki1sd:ki1ed,ki2sd:ki2ed)                  ! surface albedo

  ! Input data to be passed to *coe_th*
  REAL    (KIND=ireals   ), INTENT (IN) ::  &
 
     ! layer gas contents (cloudy and cloud-free, if distinction necessary)
     pduh2oc(ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed), & ! h2o-vapour cloudy      
     pduh2of(ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed), & ! h2o-vapour cloud-free  
     pduco2 (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed), & ! co2
     pduo3  (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed), & ! o3 
     ! optical properties of 'grey' constituents (cloudy and cloud-free)
     podsc  (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed), & ! scattering optical depth
     podsf  (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed), & ! scattering optical depth
     podac  (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed), & ! absorption optical depth
     podaf  (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed), & ! absorption optical depth
     pbsfc  (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed), & ! backscatter fraction
     pbsff  (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed), & ! backscatter fraction
      
     palogp (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed), & ! ln(p)
     palogt (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed)    ! ln(T)
 
! Output data
! -----------
  REAL    (KIND=ireals   ), INTENT (OUT) ::  &
     pflcu (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed+1), & ! flux up   cloudy
     pflfu (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed+1), & ! flux up   cloud-free 
     pflcd (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed+1), & ! flux down cloudy     
     pflfd (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed+1)    ! flux down cloud-free 

! Local parameters: 
! ----------------
  INTEGER (KIND=iintegers), PARAMETER ::  &
     j1b    = 1,           & ! debug point index first dimension
     j2b    = 1              ! debug point index second dimension

! Local scalars:
! -------------
  INTEGER (KIND=iintegers) ::  &
    j1,j2,j3                 ! loop indices over spatial dimensions

  LOGICAL                  ::  &
    ldebug_coe_th            ! debug switch for *coe_th*

  REAL    (KIND=ireals   ) ::  &
    ztd1 ,ztd2 ,ztd3 ,ztd4 ,ztd5 ,ztd6 , ztd7,  & !
    ztds1,ztds2,ztds3,ztus1                       !
 
! Local (automatic) arrays:
! ------------------------
  REAL    (KIND=ireals   ) ::  &

    ! layer properties calculated in *coe_th*
    pa1c   (ki1sd:ki1ed,ki2sd:ki2ed), & !
    pa1f   (ki1sd:ki1ed,ki2sd:ki2ed), & ! 
    pa2c   (ki1sd:ki1ed,ki2sd:ki2ed), & ! 
    pa2f   (ki1sd:ki1ed,ki2sd:ki2ed), & ! 
    pa3c   (ki1sd:ki1ed,ki2sd:ki2ed), & ! 
    pa3f   (ki1sd:ki1ed,ki2sd:ki2ed), & ! 

    ! Utility arrays
    ztu1 (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed), & !
    ztu2 (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed), & !
    ztu3 (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed), & !
    ztu4 (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed), & !
    ztu5 (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed), & !
    ztu6 (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed), & !
    ztu7 (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed), & !
    ztu8 (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed), & !
    ztu9 (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed)    !
 
!- End of header
!==============================================================================
 
!------------------------------------------------------------------------------
! Begin Subroutine inv_th              
!------------------------------------------------------------------------------

  ldebug_coe_th = .FALSE. 

! Upper boundary condition
 
  DO j3 = ki3sc, ki3ec+1
    DO j2 = ki2sc, ki2ec
      DO j1 = ki1sc, ki1ec
        pflfd(j1,j2,j3) = pbbr(j1,j2,j3)
        pflcd(j1,j2,j3) = 0.0_ireals
      ENDDO
    ENDDO
  ENDDO

  IF (ldebug) THEN
     print *,' *** INV_TH **************************'
     print *,' *** debug point : ',j1b,j2b
     print *,'pflfd(j1b,j2b,ki3sc) : ',pflfd(j1b,j2b,ki3sc)
     print *,'pflcd(j1b,j2b,ki3sc) : ',pflcd(j1b,j2b,ki3sc)
  ENDIF
 
! Determine effects of first layer in *coe_th*
  CALL coe_th ( pduh2oc,pduh2of,pduco2 ,pduo3  ,palogp ,palogt , &
                podsc  ,podsf  ,podac  ,podaf  ,pbsfc  ,pbsff  , &
                ki3sc  ,kspec  ,kh2o   ,kco2   ,ko3    ,         &
                ki1sd  ,ki1ed  ,ki2sd  ,ki2ed  ,ki3sd  ,ki3ed,   &
                ki1sc  ,ki1ec  ,ki2sc  ,ki2ec  ,ki3sc  ,ki3ec,   &
                ldebug_coe_th ,                                  &
                pa1c   ,pa1f   ,pa2c   ,pa2f   ,pa3c   ,pa3f)
 
! Set RHS
  DO j2 = ki2sc, ki2ec
    DO j1 = ki1sc, ki1ec
      pflfu(j1,j2,ki3sc)   = (1.-pclc(j1,j2,ki3sc))*pa3f(j1,j2)* &
                              (pbbr(j1,j2,ki3sc)-pbbr(j1,j2,ki3sc+1))
      pflcu(j1,j2,ki3sc)   =     pclc(j1,j2,ki3sc) *pa3c(j1,j2)* &
                              (pbbr(j1,j2,ki3sc)-pbbr(j1,j2,ki3sc+1))
      pflfd(j1,j2,ki3sc+1) = -pflfu(j1,j2,ki3sc)
      pflcd(j1,j2,ki3sc+1) = -pflcu(j1,j2,ki3sc)
    ENDDO
  ENDDO
 
! Elimination for first layer
  DO j2 = ki2sc, ki2ec
    DO j1 = ki1sc, ki1ec
      pflfu(j1,j2,ki3sc)   = pflfu(j1,j2,ki3sc  ) + pa2f(j1,j2) *            &
                            (pca2 (j1,j2,ki3sc)   * pflfd(j1,j2,ki3sc))
      pflfd(j1,j2,ki3sc+1) = pflfd(j1,j2,ki3sc+1) + pa1f(j1,j2) *            &
                            (pca2 (j1,j2,ki3sc)   * pflfd(j1,j2,ki3sc))
      pflcu(j1,j2,ki3sc)   = pflcu(j1,j2,ki3sc  ) + pa2c(j1,j2) *            &
                            (pcb2 (j1,j2,ki3sc)   * pflfd(j1,j2,ki3sc))
      pflcd(j1,j2,ki3sc+1) = pflcd(j1,j2,ki3sc+1) + pa1c(j1,j2) *            &
                            (pcb2 (j1,j2,ki3sc)   * pflfd(j1,j2,ki3sc))
    ENDDO
  ENDDO

  IF (ldebug) THEN
     print *,' after elimination'
     print *,'pflfd(j1b,j2b,ki3sc+1) : ',pflfd(j1b,j2b,ki3sc+1)
  ENDIF     
 
! Store some utitlity variables for first layer

  DO j2 = ki2sc, ki2ec
    DO j1 = ki1sc, ki1ec
      ztu1(j1,j2,ki3sc) = 0.0_ireals
      ztu2(j1,j2,ki3sc) = pca1(j1,j2,ki3sc)*pa1f(j1,j2)
      ztu3(j1,j2,ki3sc) = pcc1(j1,j2,ki3sc)*pa1f(j1,j2)
      ztu4(j1,j2,ki3sc) = pcb1(j1,j2,ki3sc)*pa1c(j1,j2)
      ztu5(j1,j2,ki3sc) = pcd1(j1,j2,ki3sc)*pa1c(j1,j2)
      ztu6(j1,j2,ki3sc) = pca1(j1,j2,ki3sc)*pa2f(j1,j2)
      ztu7(j1,j2,ki3sc) = pcc1(j1,j2,ki3sc)*pa2f(j1,j2)
      ztu8(j1,j2,ki3sc) = pcb1(j1,j2,ki3sc)*pa2c(j1,j2)
      ztu9(j1,j2,ki3sc) = pcd1(j1,j2,ki3sc)*pa2c(j1,j2)
    ENDDO
  ENDDO
 
! Vertical loop
 
  DO j3 = ki3sc+1, ki3ec
 
    ! Determine effect of the layer in *coe_th*

    CALL coe_th ( pduh2oc,pduh2of,pduco2 ,pduo3  ,palogp ,palogt , &
                  podsc  ,podsf  ,podac  ,podaf  ,pbsfc  ,pbsff  , &
                  j3     ,kspec  ,kh2o   ,kco2   ,ko3    ,         &
                  ki1sd  ,ki1ed  ,ki2sd  ,ki2ed  ,ki3sd  ,ki3ed,   &
                  ki1sc  ,ki1ec  ,ki2sc  ,ki2ec  ,ki3sc  ,ki3ec,   &
                  ldebug_coe_th ,                                  &
                  pa1c   ,pa1f   ,pa2c   ,pa2f   ,pa3c   ,pa3f)
 
    ! Set RHS
    DO j2 = ki2sc, ki2ec
      DO j1 = ki1sc, ki1ec
        pflfu(j1,j2,j3  ) = (1.0_ireals - pclc(j1,j2,j3)) * pa3f(j1,j2)      &
                                       * (pbbr(j1,j2,j3)  - pbbr(j1,j2,j3+1))
        pflcu(j1,j2,j3  ) =   pclc(j1,j2,j3) * pa3c(j1,j2)                   &
                                       * (pbbr(j1,j2,j3)  - pbbr(j1,j2,j3+1))
        pflfd(j1,j2,j3+1) = - pflfu(j1,j2,j3)
        pflcd(j1,j2,j3+1) = - pflcu(j1,j2,j3)
      ENDDO
    ENDDO
 
    IF (ldebug) THEN
       print *,' in vertical loop j3=',j3
       print *,'pflfd(j1b,j2b,j3+1) : ',pflfd(j1b,j2b,j3+1)
    ENDIF 

    ! Elimination and storage of utility variables
 
    DO j2 = ki2sc, ki2ec
      DO j1 = ki1sc, ki1ec
         ztd1 = 1./(1.-pa2f(j1,j2)*( pca2(j1,j2,j3)*ztu6(j1,j2,j3-1) &
                                    +pcc2(j1,j2,j3)*ztu8(j1,j2,j3-1)))
         pflfu(j1,j2,j3) = ztd1*( pflfu(j1,j2,j3) +                  & 
                       pa2f(j1,j2)*( pca2(j1,j2,j3)*pflfd(j1,j2,j3)  &
                                    +pcc2(j1,j2,j3)*pflcd(j1,j2,j3)))
         ztu1 (j1,j2,j3) = ztd1*                                     &
                       pa2f(j1,j2)*( pca2(j1,j2,j3)*ztu7(j1,j2,j3-1) &
                                    +pcc2(j1,j2,j3)*ztu9(j1,j2,j3-1))
         ztu2 (j1,j2,j3) = ztd1*pa1f(j1,j2)*pca1(j1,j2,j3)
         ztu3 (j1,j2,j3) = ztd1*pa1f(j1,j2)*pcc1(j1,j2,j3)
         ztd2 = pa2c(j1,j2)*(  pcb2(j1,j2,j3)*ztu6(j1,j2,j3-1)  &
                             + pcd2(j1,j2,j3)*ztu8(j1,j2,j3-1))
         ztd3 = 1./(1.-pa2c(j1,j2)*( pcb2(j1,j2,j3)*ztu7(j1,j2,j3-1)  &
                                    +pcd2(j1,j2,j3)*ztu9(j1,j2,j3-1)) &
                               -ztd2*ztu1(j1,j2,j3))
         pflcu(j1,j2,j3) = ztd3*( pflcu(j1,j2,j3) +                   & 
                       pa2c(j1,j2)*( pcb2(j1,j2,j3)*pflfd(j1,j2,j3)   &
                                    +pcd2(j1,j2,j3)*pflcd(j1,j2,j3))  &
                             + ztd2*pflfu(j1,j2,j3))
         ztu4 (j1,j2,j3) = ztd3*( pa1c(j1,j2)*pcb1(j1,j2,j3)+ztd2*ztu2(j1,j2,j3))
         ztu5 (j1,j2,j3) = ztd3*( pa1c(j1,j2)*pcd1(j1,j2,j3)+ztd2*ztu3(j1,j2,j3))
         ztd4 = pa1f(j1,j2)*( pca2(j1,j2,j3)*ztu6(j1,j2,j3-1) &
                             +pcc2(j1,j2,j3)*ztu8(j1,j2,j3-1))
         ztd5 = pa1f(j1,j2)*( pca2(j1,j2,j3)*ztu7(j1,j2,j3-1) &
                             +pcc2(j1,j2,j3)*ztu9(j1,j2,j3-1))
         pflfd(j1,j2,j3+1) = pflfd(j1,j2,j3+1)                 &
                +pa1f(j1,j2)*( pca2(j1,j2,j3)*pflfd(j1,j2,j3)  &
                              +pcc2(j1,j2,j3)*pflcd(j1,j2,j3)) &
                + ztd4*pflfu(j1,j2,j3) + ztd5*pflcu(j1,j2,j3)
         ztu6 (j1,j2,j3) = pa2f(j1,j2)*pca1(j1,j2,j3) &
                          +ztd4*ztu2(j1,j2,j3)+ztd5*ztu4(j1,j2,j3)
         ztu7 (j1,j2,j3) = pa2f(j1,j2)*pcc1(j1,j2,j3) &
                          +ztd4*ztu3(j1,j2,j3)+ztd5*ztu5(j1,j2,j3)
         ztd6 = pa1c(j1,j2)*( pcb2(j1,j2,j3)*ztu6(j1,j2,j3-1) &
                             +pcd2(j1,j2,j3)*ztu8(j1,j2,j3-1))
         ztd7 = pa1c(j1,j2)*( pcb2(j1,j2,j3)*ztu7(j1,j2,j3-1) &
                             +pcd2(j1,j2,j3)*ztu9(j1,j2,j3-1))
         pflcd(j1,j2,j3+1) = pflcd(j1,j2,j3+1)                 &
                +pa1c(j1,j2)*( pcb2(j1,j2,j3)*pflfd(j1,j2,j3)  &
                              +pcd2(j1,j2,j3)*pflcd(j1,j2,j3)) &
                + ztd6*pflfu(j1,j2,j3) + ztd7*pflcu(j1,j2,j3)
         ztu8(j1,j2,j3) = pa2c(j1,j2)*pcb1(j1,j2,j3) &
                         +ztd6*ztu2(j1,j2,j3)+ztd7*ztu4(j1,j2,j3)
         ztu9(j1,j2,j3) = pa2c(j1,j2)*pcd1(j1,j2,j3) &
                         +ztd6*ztu3(j1,j2,j3)+ztd7*ztu5(j1,j2,j3)
      ENDDO
    ENDDO
    IF (ldebug) THEN
       print *,' after elimination in vertical loop j3=',j3
       print *,'pflfd(j1b,j2b,j3+1) : ',pflfd(j1b,j2b,j3+1)
    ENDIF  
 
  ENDDO     ! End of vertical loop over layers
 
  ! Elimination and backsubstitution at the surface

  DO j2 = ki2sc, ki2ec
    DO j1 = ki1sc, ki1ec
       ztds1    =1./(1.-palth(j1,j2)*ztu6(j1,j2,ki3ec))
       pflfu(j1,j2,ki3ec+1)= ztds1 *palth(j1,j2)*pflfd(j1,j2,ki3ec+1)
       ztus1    =ztds1 *palth(j1,j2)*ztu7(j1,j2,ki3ec)
       ztds2    =palth(j1,j2)*ztu8(j1,j2,ki3ec)
       ztds3    =1./(1.-palth(j1,j2)*ztu9(j1,j2,ki3ec)-ztds2*ztus1)
       pflcu(j1,j2,ki3ec+1)=ztds3*( palth(j1,j2)*pflcd(j1,j2,ki3ec+1) &
                                   +ztds2       *pflfu(j1,j2,ki3ec+1))
       pflfu(j1,j2,ki3ec+1)=pflfu(j1,j2,ki3ec+1)+ztus1*pflcu(j1,j2,ki3ec+1)
    ENDDO
  ENDDO

!     Layer-by-layer backsubstitution
 
  DO j3 =ki3ec,ki3sc,-1
    DO j2 = ki2sc, ki2ec
      DO j1 = ki1sc, ki1ec
        pflcd(j1,j2,j3+1) = pflcd(j1,j2,j3+1) + ztu8 (j1,j2,j3)               &
                   * pflfu(j1,j2,j3+1) + ztu9 (j1,j2,j3) * pflcu(j1,j2,j3+1)
        pflfd(j1,j2,j3+1) = pflfd(j1,j2,j3+1) + ztu6 (j1,j2,j3)               &
                   * pflfu(j1,j2,j3+1) + ztu7 (j1,j2,j3) * pflcu(j1,j2,j3+1)
        pflcu(j1,j2,j3  ) = pflcu(j1,j2,j3  ) + ztu4 (j1,j2,j3)               &
                   * pflfu(j1,j2,j3+1) + ztu5 (j1,j2,j3) * pflcu(j1,j2,j3+1)
        pflfu(j1,j2,j3  ) = pflfu(j1,j2,j3  ) + ztu2 (j1,j2,j3)               &
                   * pflfu(j1,j2,j3+1) + ztu3 (j1,j2,j3) * pflcu(j1,j2,j3+1)  &
                                       + ztu1 (j1,j2,j3) * pflcu(j1,j2,j3)
      ENDDO
    ENDDO

     IF (ldebug) THEN
        print *,' after backsubst.  in vertical loop j3=',j3
        print *,'pflfd(j1b,j2b,j3+1) : ',pflfd(j1b,j2b,j3+1)
     ENDIF
  ENDDO
 
!------------------------------------------------------------------------------
! End of the subroutine
!------------------------------------------------------------------------------

END SUBROUTINE inv_th

!==============================================================================
!==============================================================================
!+ Module procedure in "Radiation" 
!------------------------------------------------------------------------------

SUBROUTINE inv_so (                                                    &
       pclc   ,pca1   ,pca2  ,pcb1  ,pcb2  ,pcc1  ,pcc2 ,pcd1 ,pcd2 ,  &
       pflpt  ,psmu0  ,pqsmu0,palp  ,palso ,                           &
       pduh2oc,pduh2of,pduco2,pduo3 ,palogp,palogt,                    &
       podsc  ,podsf  ,podac ,podaf ,pbsfc ,pbsff ,pusfc,pusff,        &
       kspec  ,kh2o   ,kco2  ,ko3   ,                                  &
       ki1sd  ,ki1ed  ,ki2sd ,ki2ed ,ki3sd ,ki3ed ,                    &
       ki1sc  ,ki1ec  ,ki2sc ,ki2ec ,ki3sc ,ki3ec , ldebug ,           &
       pflcu  ,pflfu  ,pflcd ,pflfd ,pflcp ,pflfp)

!------------------------------------------------------------------------------
!
! Description:
!
!   The module procedure inv_so solves the linear system of equations for 
!   solar fluxes.
!   The routine solves a linear equation system for solar fluxes using
!   a Gaussian elimination-backsubstitution algorithm dedicated to the
!   specific structure of the system matrix.
!
! Method:
!
! - setting of the RHS of the system using the parallel solar radiation
!   at the top of the atmosphere and allowing for partial cloud cover
! - solution of the equation system including the lower boundary
!   condition
! - matrix coefficients are calculated in the course of the elimination
!   step for one layer at a time through a call to routine *coe_so*
! - the final result, i.e. upward and downward diffuse and parallel 
!   solar fluxes are stored seperately for cloudy and cloud-free parts
!   of each layer boundary
!     
!------------------------------------------------------------------------------

! Subroutine arguments:
! --------------------

! Input data
! ----------
  INTEGER (KIND=iintegers), INTENT (IN) ::  &
     ki1sd,       & ! start index for first  array dimension
     ki1ed,       & ! end   index for first  array dimension
     ki2sd,       & ! start index for second array dimension
     ki2ed,       & ! end   index for second array dimension
     ki3sd,       & ! start index for third  array dimension
     ki3ed,       & ! end   index for third  array dimension

   ! and the same for the computations
     ki1sc,       & ! start index for first  array computation
     ki1ec,       & ! end   index for first  array computation
     ki2sc,       & ! start index for second array computation
     ki2ec,       & ! end   index for second array computation
     ki3sc,       & ! start index for third  array computation
     ki3ec,       & ! end   index for third  array computation
     kspec,       & ! spectral interval considered
     kh2o ,       & ! table index for h2o absorption properties
     kco2 ,       & ! table index for co2 absorption properties
     ko3            ! table index for o3  absorption properties

  LOGICAL                 , INTENT (IN) ::  &
     ldebug         ! debug control switch       

  REAL    (KIND=ireals   ), INTENT (IN) ::  &
     pclc  (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed), & ! cloud cover
     pca1  (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed), & ! cloud geometry factor  
     pca2  (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed), & ! cloud geometry factor  
     pcb1  (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed), & ! cloud geometry factor  
     pcb2  (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed), & ! cloud geometry factor  
     pcc1  (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed), & ! cloud geometry factor  
     pcc2  (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed), & ! cloud geometry factor  
     pcd1  (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed), & ! cloud geometry factor  
     pcd2  (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed), & ! cloud geometry factor  

     pflpt (ki1sd:ki1ed,ki2sd:ki2ed), &  ! parallel solar flux at TOA
     palp  (ki1sd:ki1ed,ki2sd:ki2ed), &  ! surface albedo for parallel
     palso (ki1sd:ki1ed,ki2sd:ki2ed)     ! and for diffuse radiation  

  ! Input data to be passed to *coe_so*
  REAL    (KIND=ireals   ), INTENT (IN) ::  &
 
     ! layer gas contents (cloudy and cloud-free, if distinction necessary)
     pduh2oc(ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed), & ! h2o-vapour cloudy      
     pduh2of(ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed), & ! h2o-vapour cloud-free  
     pduco2 (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed), & ! co2
     pduo3  (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed), & ! o3 
     ! optical properties of 'grey' constituents (cloudy and cloud-free)
     podsc  (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed), & ! scattering optical depth
     podsf  (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed), & ! scattering optical depth
     podac  (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed), & ! absorption optical depth
     podaf  (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed), & ! absorption optical depth
     pbsfc  (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed), & ! backscatter fraction
     pbsff  (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed), & ! backscatter fraction
     pusfc  (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed), & ! upscatter   fraction
     pusff  (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed), & ! upscatter   fraction
      
     palogp (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed), & ! ln(p)
     palogt (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed), & ! ln(T)
 
     psmu0 (ki1sd:ki1ed,ki2sd:ki2ed), & ! cosine of zenith angle
     pqsmu0(ki1sd:ki1ed,ki2sd:ki2ed)    ! 1./cosine of zenith angle

! Output data
! -----------
  REAL    (KIND=ireals   ), INTENT (OUT) ::  &
     pflcu (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed+1), & ! flux up   cloudy
     pflfu (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed+1), & ! flux up   cloud-free 
     pflcd (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed+1), & ! flux down cloudy     
     pflfd (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed+1), & ! flux down cloud-free 
     pflcp (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed+1), & ! flux par. cloudy     
     pflfp (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed+1)    ! flux par. cloud-free 

! Local parameters: 
! ----------------
  INTEGER (KIND=iintegers), PARAMETER ::  &
     j1b    = 1,           & ! debug point index first dimension
     j2b    = 1              ! debug point index second dimension

! Local scalars:
! -------------
  INTEGER (KIND=iintegers) ::  &
    j1,j2,j3                 ! loop indices over spatial dimensions

  LOGICAL                  ::  &
    ldebug_coe_so            ! debug switch for *coe_so*

  REAL    (KIND=ireals   ) ::  &
    ztd1 ,ztd2 ,ztd3 ,ztd4 ,ztd5 ,ztd6 , ztd7,  & !
    ztds1,ztds2,ztds3,ztus1                       !
 
! Local (automatic) arrays:
! ------------------------
  REAL    (KIND=ireals   ) ::  &

    ! layer properties calculated in *coe_so*
    pa1c   (ki1sd:ki1ed,ki2sd:ki2ed), & !
    pa1f   (ki1sd:ki1ed,ki2sd:ki2ed), & ! 
    pa2c   (ki1sd:ki1ed,ki2sd:ki2ed), & ! 
    pa2f   (ki1sd:ki1ed,ki2sd:ki2ed), & ! 
    pa3c   (ki1sd:ki1ed,ki2sd:ki2ed), & ! 
    pa4f   (ki1sd:ki1ed,ki2sd:ki2ed), & ! 
    pa4c   (ki1sd:ki1ed,ki2sd:ki2ed), & ! 
    pa5f   (ki1sd:ki1ed,ki2sd:ki2ed), & ! 
    pa5c   (ki1sd:ki1ed,ki2sd:ki2ed), & ! 
    pa3f   (ki1sd:ki1ed,ki2sd:ki2ed), & ! 

    ! Utility arrays
    ztu1 (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed), & !
    ztu2 (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed), & !
    ztu3 (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed), & !
    ztu4 (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed), & !
    ztu5 (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed), & !
    ztu6 (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed), & !
    ztu7 (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed), & !
    ztu8 (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed), & !
    ztu9 (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed)    !
 
!- End of header
!==============================================================================
 
!------------------------------------------------------------------------------
! Begin Subroutine inv_so              
!------------------------------------------------------------------------------

  ldebug_coe_so = .FALSE. 

  ! Upper boundary condition
 
  DO j2 = ki2sc, ki2ec
    DO j1 = ki1sc, ki1ec
      pflfp(j1,j2,ki3sc) = pflpt(j1,j2)
      pflcp(j1,j2,ki3sc) = 0.0_ireals
      pflfd(j1,j2,ki3sc) = 0.0_ireals
      pflcd(j1,j2,ki3sc) = 0.0_ireals
    ENDDO
  ENDDO

  IF (ldebug) THEN
     print *,' *** INV_SO **************************'
     print *,' *** Debug point: ',j1b,j2b
     print *,'pflfp(j1b,j2b,ki3sc) : ',pflfp(j1b,j2b,ki3sc)
     print *,'pflcp(j1b,j2b,ki3sc) : ',pflcp(j1b,j2b,ki3sc)
     print *,'pflfp(j1b,j2b,ki3sc) : ',pflfp(j1b,j2b,ki3sc)
     print *,'pflcd(j1b,j2b,ki3sc) : ',pflcd(j1b,j2b,ki3sc)
  ENDIF 

  ! Determine effects of first layer in *coe_so*

  CALL  coe_so (                                                     &
      pduh2oc,pduh2of,pduco2 ,pduo3  ,palogp ,palogt ,               &
      podsc  ,podsf  ,podac  ,podaf  ,pbsfc  ,pbsff  ,pusfc ,pusff , &
      psmu0  ,pqsmu0 ,                                               &
      ki3sc  ,kspec  ,kh2o   ,kco2   ,ko3    ,                       &
      ki1sd  ,ki1ed  ,ki2sd  ,ki2ed  ,ki3sd  ,ki3ed,                 &
      ki1sc  ,ki1ec  ,ki2sc  ,ki2ec  ,ki3sc  ,ki3ec,                 &
      ldebug_coe_so ,                                                &
      pa1c   ,pa1f   ,pa2c   ,pa2f   ,pa3c   ,pa3f ,                 &
      pa4c   ,pa4f   ,pa5c   ,pa5f )
 
  ! Top layer elimination
 
  DO j2 = ki2sc, ki2ec
    DO j1 = ki1sc, ki1ec
      pflfu(j1,j2,ki3sc  ) = pa3f(j1,j2) * pca2(j1,j2,ki3sc) *pflfp(j1,j2,ki3sc)
      pflfp(j1,j2,ki3sc+1) = pa1f(j1,j2) * pca2(j1,j2,ki3sc) *pflfp(j1,j2,ki3sc)
      pflfd(j1,j2,ki3sc+1) = pa2f(j1,j2) * pca2(j1,j2,ki3sc) *pflfp(j1,j2,ki3sc)
      pflcu(j1,j2,ki3sc  ) = pa3c(j1,j2) * pcb2(j1,j2,ki3sc) *pflfp(j1,j2,ki3sc)
      pflcp(j1,j2,ki3sc+1) = pa1c(j1,j2) * pcb2(j1,j2,ki3sc) *pflfp(j1,j2,ki3sc)
      pflcd(j1,j2,ki3sc+1) = pa2c(j1,j2) * pcb2(j1,j2,ki3sc) *pflfp(j1,j2,ki3sc)
    ENDDO
  ENDDO

  IF (ldebug) THEN
     print *,' *** INV_SO **************************'
     print *,'pflfu(j1b,j2b,ki3sc)  : ',pflfu(j1b,j2b,ki3sc)
     print *,'pflcu(j1b,j2b,ki3sc)  : ',pflcu(j1b,j2b,ki3sc)
     print *,'pflfd(j1b,j2b,ki3sc+1): ',pflfd(j1b,j2b,ki3sc+1)
     print *,'pflcd(j1b,j2b,ki3sc+1): ',pflcd(j1b,j2b,ki3sc+1)
     print *,'pa1f (j1b,j2b)        : ',pa1f (j1b,j2b)        
     print *,'pa1c (j1b,j2b)        : ',pa1c (j1b,j2b)        
     print *,'pa2f (j1b,j2b)        : ',pa2f (j1b,j2b)        
     print *,'pa2c (j1b,j2b)        : ',pa2c (j1b,j2b)        
     print *,'pa3f (j1b,j2b)        : ',pa3f (j1b,j2b)        
     print *,'pa3c (j1b,j2b)        : ',pa3c (j1b,j2b)        
  ENDIF 
 
  ! Storage of utility arrays for the top layer

  DO j2 = ki2sc, ki2ec
    DO j1 = ki1sc, ki1ec
      ztu1(j1,j2,1) = 0.0_ireals
      ztu2(j1,j2,1) = pca1(j1,j2,1) * pa4f(j1,j2)
      ztu3(j1,j2,1) = pcc1(j1,j2,1) * pa4f(j1,j2)
      ztu4(j1,j2,1) = pcb1(j1,j2,1) * pa4c(j1,j2)
      ztu5(j1,j2,1) = pcd1(j1,j2,1) * pa4c(j1,j2)
      ztu6(j1,j2,1) = pca1(j1,j2,1) * pa5f(j1,j2)
      ztu7(j1,j2,1) = pcc1(j1,j2,1) * pa5f(j1,j2)
      ztu8(j1,j2,1) = pcb1(j1,j2,1) * pa5c(j1,j2)
      ztu9(j1,j2,1) = pcd1(j1,j2,1) * pa5c(j1,j2)
    ENDDO
  ENDDO
 
  ! Suczessive layer-by-layer elimination
 
  DO j3 = ki3sc+1, ki3ec         ! Loop over vertical

     ! Determine effects of layer in *coe_so*
      CALL  coe_so (                                                     &
          pduh2oc,pduh2of,pduco2 ,pduo3  ,palogp ,palogt ,               &
          podsc  ,podsf  ,podac  ,podaf  ,pbsfc  ,pbsff  ,pusfc ,pusff , &
          psmu0  ,pqsmu0 ,                                               &
          j3     ,kspec  ,kh2o   ,kco2   ,ko3    ,                       &
          ki1sd  ,ki1ed  ,ki2sd  ,ki2ed  ,ki3sd  ,ki3ed,                 &
          ki1sc  ,ki1ec  ,ki2sc  ,ki2ec  ,ki3sc  ,ki3ec,                 &
          ldebug_coe_so ,                                                &
          pa1c   ,pa1f   ,pa2c   ,pa2f   ,pa3c   ,pa3f ,                 &
          pa4c   ,pa4f   ,pa5c   ,pa5f )
 
     ! Elimination
 
     DO j2 = ki2sc, ki2ec
       DO j1 = ki1sc, ki1ec
         ztd1 = 1./(1.-pa5f(j1,j2)*(pca2(j1,j2,j3)*ztu6(j1,j2,j3-1)    &
                                   +pcc2(j1,j2,j3)*ztu8(j1,j2,j3-1)))
         pflfu(j1,j2,j3) = ztd1*(                                       &
                         pa3f(j1,j2)*( pca2(j1,j2,j3)*pflfp(j1,j2,j3)   &
                                      +pcc2(j1,j2,j3)*pflcp(j1,j2,j3) ) &
                        +pa5f(j1,j2)*( pca2(j1,j2,j3)*pflfd(j1,j2,j3)   &
                                      +pcc2(j1,j2,j3)*pflcd(j1,j2,j3) )   )
         ztu1 (j1,j2,j3) = ztd1*pa5f(j1,j2)*                            &
                                     ( pca2(j1,j2,j3)*ztu7(j1,j2,j3-1)  &
                                      +pcc2(j1,j2,j3)*ztu9(j1,j2,j3-1))
         ztu2(j1,j2,j3)  =  ztd1*pa4f(j1,j2)*pca1(j1,j2,j3)
         ztu3(j1,j2,j3)  =  ztd1*pa4f(j1,j2)*pcc1(j1,j2,j3)
         ztd2            = pa5c(j1,j2)*( pcb2(j1,j2,j3)*ztu6(j1,j2,j3-1) &
                                        +pcd2(j1,j2,j3)*ztu8(j1,j2,j3-1) )
         ztd3            = 1./( 1.-pa5c(j1,j2)*(pcb2(j1,j2,j3)*ztu7(j1,j2,j3-1) &
                                               +pcd2(j1,j2,j3)*ztu9(j1,j2,j3-1))&
                                - ztd2*ztu1(j1,j2,j3) )
         pflcu(j1,j2,j3) = ztd3 *(                                               &
                                  pa3c(j1,j2)*( pcb2(j1,j2,j3)*pflfp(j1,j2,j3)   &
                                               +pcd2(j1,j2,j3)*pflcp(j1,j2,j3) ) &
                                 +pa5c(j1,j2)*( pcb2(j1,j2,j3)*pflfd(j1,j2,j3)   &
                                               +pcd2(j1,j2,j3)*pflcd(j1,j2,j3) ) &
                                 +ztd2*pflfu(j1,j2,j3)   )
         ztu4(j1,j2,j3)  = ztd3 *( pa4c(j1,j2)*pcb1(j1,j2,j3)+ztd2*ztu2(j1,j2,j3) )
         ztu5(j1,j2,j3)  = ztd3 *( pa4c(j1,j2)*pcd1(j1,j2,j3)+ztd2*ztu3(j1,j2,j3) )
         pflfp(j1,j2,j3+1) = pa1f(j1,j2)*(pca2(j1,j2,j3)*pflfp(j1,j2,j3) &
                                         +pcc2(j1,j2,j3)*pflcp(j1,j2,j3))
         pflcp(j1,j2,j3+1) = pa1c(j1,j2)*(pcb2(j1,j2,j3)*pflfp(j1,j2,j3) &
                                         +pcd2(j1,j2,j3)*pflcp(j1,j2,j3))
         ztd4 = pa4f(j1,j2)*( pca2(j1,j2,j3)*ztu6(j1,j2,j3-1) &
                             +pcc2(j1,j2,j3)*ztu8(j1,j2,j3-1) )
         ztd5 = pa4f(j1,j2)*( pca2(j1,j2,j3)*ztu7(j1,j2,j3-1) &
                             +pcc2(j1,j2,j3)*ztu9(j1,j2,j3-1) )
         pflfd(j1,j2,j3+1) = pa2f(j1,j2)*(pca2(j1,j2,j3)*pflfp(j1,j2,j3)  &
                                         +pcc2(j1,j2,j3)*pflcp(j1,j2,j3)) &
                            +pa4f(j1,j2)*(pca2(j1,j2,j3)*pflfd(j1,j2,j3)  &
                                         +pcc2(j1,j2,j3)*pflcd(j1,j2,j3)) &
                            +ztd4*pflfu(j1,j2,j3) + ztd5*pflcu(j1,j2,j3)
         ztu6(j1,j2,j3) = pa5f(j1,j2)*pca1(j1,j2,j3) &
                          + ztd4*ztu2(j1,j2,j3) + ztd5*ztu4(j1,j2,j3)
         ztu7(j1,j2,j3) = pa5f(j1,j2)*pcc1(j1,j2,j3) &
                          + ztd4*ztu3(j1,j2,j3) + ztd5*ztu5(j1,j2,j3)
         ztd6 = pa4c(j1,j2)*( pcb2(j1,j2,j3)*ztu6(j1,j2,j3-1) &
                             +pcd2(j1,j2,j3)*ztu8(j1,j2,j3-1) )
         ztd7 = pa4c(j1,j2)*( pcb2(j1,j2,j3)*ztu7(j1,j2,j3-1) &
                             +pcd2(j1,j2,j3)*ztu9(j1,j2,j3-1) )
         pflcd(j1,j2,j3+1) =  pa2c(j1,j2)*(pcb2(j1,j2,j3)*pflfp(j1,j2,j3)  &
                                          +pcd2(j1,j2,j3)*pflcp(j1,j2,j3)) &
                            + pa4c(j1,j2)*(pcb2(j1,j2,j3)*pflfd(j1,j2,j3)  &
                                          +pcd2(j1,j2,j3)*pflcd(j1,j2,j3)) &
                            +ztd6*pflfu(j1,j2,j3) + ztd7*pflcu(j1,j2,j3)
         ztu8(j1,j2,j3) = pa5c(j1,j2)*pcb1(j1,j2,j3) &
                          + ztd6*ztu2(j1,j2,j3) + ztd7*ztu4(j1,j2,j3)
         ztu9(j1,j2,j3) = pa5c(j1,j2)*pcd1(j1,j2,j3) &
                          + ztd6*ztu3(j1,j2,j3) + ztd7*ztu5(j1,j2,j3)
       ENDDO
     ENDDO
     IF (ldebug) THEN
       print *,' inv_so j3=',j3
       print *,'pflfu(j1b,j2b,j3)  : ',pflfu(j1b,j2b,j3)
       print *,'pflcu(j1b,j2b,j3)  : ',pflcu(j1b,j2b,j3)
       print *,'pflfd(j1b,j2b,j3+1): ',pflfd(j1b,j2b,j3+1)
       print *,'pflcd(j1b,j2b,j3+1): ',pflcd(j1b,j2b,j3+1)
       print *,'pflfp(j1b,j2b,j3+1): ',pflfp(j1b,j2b,j3+1)
       print *,'pflcp(j1b,j2b,j3+1): ',pflcp(j1b,j2b,j3+1)
       print *,'ztu1 (j1b,j2b,j3)  : ',ztu1 (j1b,j2b,j3)
       print *,'ztu2 (j1b,j2b,j3)  : ',ztu2 (j1b,j2b,j3)
       print *,'ztu3 (j1b,j2b,j3)  : ',ztu3 (j1b,j2b,j3)
       print *,'ztu4 (j1b,j2b,j3)  : ',ztu4 (j1b,j2b,j3)
       print *,'ztu5 (j1b,j2b,j3)  : ',ztu5 (j1b,j2b,j3)
       print *,'ztu6 (j1b,j2b,j3)  : ',ztu6 (j1b,j2b,j3)
       print *,'ztu7 (j1b,j2b,j3)  : ',ztu7 (j1b,j2b,j3)
       print *,'ztu8 (j1b,j2b,j3)  : ',ztu8 (j1b,j2b,j3)
       print *,'ztu9 (j1b,j2b,j3)  : ',ztu9 (j1b,j2b,j3)
       print *,' .....'
     ENDIF 

  END DO       ! Vertical loop
 
  ! Elimination and back-substitution at surface
 
  DO j2 = ki2sc, ki2ec
    DO j1 = ki1sc, ki1ec
       ztds1  = 1./(1.-palso(j1,j2)*ztu6(j1,j2,ki3ec))
       pflfu(j1,j2,ki3ec+1)= ztds1 *(palp (j1,j2)*pflfp(j1,j2,ki3ec+1) &
                                    +palso(j1,j2)*pflfd(j1,j2,ki3ec+1))
       ztus1  =  ztds1*palso(j1,j2)*ztu7(j1,j2,ki3ec)
       ztds2  =        palso(j1,j2)*ztu8(j1,j2,ki3ec)
       ztds3  = 1./(1.-palso(j1,j2)*ztu9(j1,j2,ki3ec)-ztds2*ztus1)
       pflcu(j1,j2,ki3ec+1) = ztds3 *(palp (j1,j2)*pflcp(j1,j2,ki3ec+1) &
                                     +palso(j1,j2)*pflcd(j1,j2,ki3ec+1) &
                                     +ztds2       *pflfu(j1,j2,ki3ec+1))
       pflfu(j1,j2,ki3ec+1) = pflfu(j1,j2,ki3ec+1) +ztus1*pflcu(j1,j2,ki3ec+1)
    ENDDO
  ENDDO  

  IF (ldebug) THEN
     print *,' inv_so surface ------------------------------'
     print *,'pflfu(j1b,j2b,ki3ec+1) : ',pflfu(j1b,j2b,ki3ec+1)
     print *,'pflcu(j1b,j2b,ki3ec+1) : ',pflcu(j1b,j2b,ki3ec+1)
     print *,'palp (j1b,j2b): ',palp (j1b,j2b)     
     print *,'palso(j1b,j2b): ',palso(j1b,j2b)     
     print *,'ztds1               ',ztds1                  
     print *,'ztds2               ',ztds2                  
     print *,'ztds3               ',ztds3                  
     print *,'ztus1               ',ztus1                  
     print *,' .....'
  ENDIF    

  ! Layer-by-layer backsubstitution

  IF (ldebug) print *,' inv_so BACKSUBSTITUTION'

!------------------------------------------------------------------------------
  DO j3 = ki3ec, ki3sc, -1
    DO j2 = ki2sc, ki2ec
      DO j1 = ki1sc, ki1ec
        pflcd(j1,j2,j3+1) =  pflcd(j1,j2,j3+1)                                &
                                   + ztu8(j1,j2,j3)*pflfu(j1,j2,j3+1)         &
                                   + ztu9(j1,j2,j3)*pflcu(j1,j2,j3+1)
        pflfd(j1,j2,j3+1) =  pflfd(j1,j2,j3+1)                                &
                                   + ztu6(j1,j2,j3)*pflfu(j1,j2,j3+1)         &
                                   + ztu7(j1,j2,j3)*pflcu(j1,j2,j3+1)
        pflcu(j1,j2,j3  ) =  pflcu(j1,j2,j3  )                                &
                                   + ztu4(j1,j2,j3)*pflfu(j1,j2,j3+1)         &
                                   + ztu5(j1,j2,j3)*pflcu(j1,j2,j3+1)
        pflfu(j1,j2,j3  ) =  pflfu(j1,j2,j3  )                                &
                                   + ztu2(j1,j2,j3)*pflfu(j1,j2,j3+1)         &
                                   + ztu3(j1,j2,j3)*pflcu(j1,j2,j3+1)         &
                                   + ztu1(j1,j2,j3)*pflcu(j1,j2,j3)
      ENDDO
    ENDDO

    IF (ldebug) THEN
      print *,' inv_so j3=',j3
      print *,'pflfu(j1b,j2b,j3)  : ',pflfu(j1b,j2b,j3)
      print *,'pflcu(j1b,j2b,j3)  : ',pflcu(j1b,j2b,j3)
      print *,'pflfd(j1b,j2b,j3+1): ',pflfd(j1b,j2b,j3+1)
      print *,'pflcd(j1b,j2b,j3+1): ',pflcd(j1b,j2b,j3+1)
    ENDIF  
  ENDDO  
 
!------------------------------------------------------------------------------
! End of the subroutine
!------------------------------------------------------------------------------

END SUBROUTINE inv_so

!==============================================================================
!==============================================================================
!+ Module procedure in "Radiation" 
!------------------------------------------------------------------------------

SUBROUTINE opt_th(prholwc,pdulwc,prhoiwc,pduiwc,                &
                  paeq1  ,paeq2 ,paeq3  ,paeq4 , paeq5 ,        &
                  ki1sd  ,ki1ed ,ki2sd  ,ki2ed , ki3sd ,ki3ed,  &
                  kspec  ,ki1sc ,ki1ec  ,ki2sc , ki2ec ,        &
                          ki3sc ,ki3ec  ,ldebug,                &
                  podac  ,podaf ,podsc  ,podsf , pbsfc ,pbsff  )

!------------------------------------------------------------------------------
!
! Description:
!
!   The module procedure opt_th calculates the optical properties of the 
!   non-gaseous constituents for one spectral interval in the thermal part 
!   of the spectrum.
!   Purpose is the calculation of (absorption and scattering) optical
!   depth and backward scattered fraction of diffuse radiation (excluding 
!   the contribution by gaseous constituents).
!
! Method:
!
! - determination of optical properies (i.e. extinction coefficient,
!   single scattering albedo and asymetry factor of the phase function)
!   using approximate relations for cloud water and cloud ice and
!   combination of five type of aerosols
!
! - calculation of optical depth (scattering and absorption) and back-
!   scattered fraction suitable for use in an implicit delta-two-stream
!   scheme
!
!------------------------------------------------------------------------------

! Subroutine arguments:
! --------------------

! Input data
! ----------
  INTEGER (KIND=iintegers), INTENT (IN) ::  &
     ki1sd,       & ! start index for first  array dimension
     ki1ed,       & ! end   index for first  array dimension
     ki2sd,       & ! start index for second array dimension
     ki2ed,       & ! end   index for second array dimension
     ki3sd,       & ! start index for third  array dimension
     ki3ed,       & ! end   index for third  array dimension
     kspec,       & ! selected spectral interval

   ! and the same for the computations
     ki1sc,       & ! start index for first  array computation
     ki1ec,       & ! end   index for first  array computation
     ki2sc,       & ! start index for second array computation
     ki2ec,       & ! end   index for second array computation
     ki3sc,       & ! start index for third  array computation
     ki3ec          ! end   index for third  array computation

  LOGICAL                 , INTENT (IN) ::  &
     ldebug         ! debug control switch       

  REAL    (KIND=ireals   ), INTENT (IN) ::  &
           !  Liquid and ice water density and content within for the cloudy
           !  part of each layer
     prholwc(ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed), &
     prhoiwc(ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed), &
     pdulwc (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed), &
     pduiwc (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed), &

           !  Aerosole contents (optical depths at 0.55 micrometer) for 5 types
     paeq1  (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed), &
     paeq2  (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed), &
     paeq3  (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed), &
     paeq4  (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed), &
     paeq5  (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed)   

! Output data
! -----------
  REAL    (KIND=ireals   ), INTENT (OUT) ::  &
     podac (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed), & ! absorption optical depth
     podaf (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed), & ! in cloudy and free part
     podsc (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed), & ! scattering optical depth
     podsf (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed), & ! in cloudy and free part
     pbsfc (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed), & ! backscattering fraction 
     pbsff (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed)    ! in cloudy and free part

! Local parameters: 
! ----------------
  REAL    (KIND=ireals   ), PARAMETER ::  &
     z1dg   = 1.0/9.80665, & ! 1./g
     z1d8   = 0.125      , & ! 1./8 
     zepopd = 1.0E-6     , & ! Security constant for optical depth
     zepssa = 1.0E-6         ! Security constant for single scattering albedo

  INTEGER (KIND=iintegers), PARAMETER ::  &
     j1b    = 22,           & ! debug point index first dimension
     j2b    = 1              ! debug point index second dimension

! Local scalars:
! -------------
  INTEGER (KIND=iintegers) ::  &
    j1,j2,j3,              & ! loop indices over spatial dimensions
    ja                       ! local loop index

  REAL    (KIND=ireals   ) ::  &
    ! individual optical properties of liquid water and ice 
    z_lwe, z_iwe,          & ! extinction coefficient
    z_lww, z_iww,          & ! single scattering coefficient
    z_lwg, z_iwg,          & ! asymetry factor
    z_lwf, z_iwf,          & ! forward scattered fraction 
    zzg
 
! Local (automatic) arrays:
! ------------------------
!    optical properties (absorption, scattering, backscatter fraction)
!    for liquid water, ice and total aerosole

  REAL    (KIND=ireals   ) ::  &
     zlwoda (ki1sd:ki1ed,ki2sd:ki2ed), & !
     zlwods (ki1sd:ki1ed,ki2sd:ki2ed), & ! 
     zlwb0  (ki1sd:ki1ed,ki2sd:ki2ed), & ! 
     ziwoda (ki1sd:ki1ed,ki2sd:ki2ed), & ! 
     ziwods (ki1sd:ki1ed,ki2sd:ki2ed), & ! 
     ziwb0  (ki1sd:ki1ed,ki2sd:ki2ed), & ! 
     zaeoda (ki1sd:ki1ed,ki2sd:ki2ed), & ! 
     zaeods (ki1sd:ki1ed,ki2sd:ki2ed), & ! 
     zaeb0  (ki1sd:ki1ed,ki2sd:ki2ed)    ! 
 
!- End of header
!==============================================================================
 
!------------------------------------------------------------------------------
! Begin Subroutine opt_th              
!------------------------------------------------------------------------------


  IF (ldebug) THEN
     print *,' **** opt-th   start ********************'
     print *,' **** debug point : ',j1b,j2b
  ENDIF     
 
  DO ja=1,5
    zaef(kspec,ja)  = zaeg(kspec,ja)**2 ! forward sc.fraction f.aerosols
  ENDDO

  IF (ldebug) THEN
     DO ja=1,5
       print *,'ja, zaef(kspec,ja): ',ja,zaef(kspec,ja)
     ENDDO
  ENDIF 

  IF (ldebug) print *,' In opt-th   vor vertical loop'

! Vertical loop
! ------------- 

  DO j3 = ki3sc, ki3ec
 
    IF (ldebug) print *,' In opt-th   j3 = ',j3  

    ! Optical properties of liquid water and ice as function of the specific
    ! liquid water and ice content
 
    DO j2 = ki2sc, ki2ec
      DO j1 = ki1sc, ki1ec

      ! Liquid water

      z_lwg    = zlwg(1,kspec) + zlwg(2,kspec)*prholwc(j1,j2,j3)
      z_lwg    = MAX (0.0_ireals,MIN(1.0_ireals,z_lwg))
      z_lwf    = z_lwg*z_lwg
      z_lww    = zlww(1,kspec) + zlww(2,kspec)*prholwc(j1,j2,j3)
      z_lww    = MAX(zepssa,MIN(1.0_ireals,z_lww)) 
      z_lwe    = z1dg * (zlwe(1,kspec) + zlwe(2,kspec)/ &
                 (zlwe(3,kspec)*prholwc(j1,j2,j3)+zlwe(4,kspec)))
      z_lwe    = MAX(zlwemn(kspec),MIN(zlwemx(kspec),z_lwe))

      zlwoda(j1,j2)= z_lwe *pdulwc(j1,j2,j3) * (1.-z_lww) 
      zlwods(j1,j2)= z_lwe *pdulwc(j1,j2,j3) *     z_lww  * (1.-z_lwf)
      zlwb0 (j1,j2)= z1d8*(4.+z_lwg)/(1.+z_lwg)
 
      ! Ice
 
      z_iwg    = ziwg(1,kspec) + ziwg(2,kspec)* LOG(prhoiwc(j1,j2,j3))
      z_iwg    = MAX(0.0_ireals,MIN(1.0_ireals,z_iwg))
      z_iwf    = z_iwg*z_iwg
      z_iww    = ziww(1,kspec) + ziww(2,kspec)* LOG(prhoiwc(j1,j2,j3))
      z_iww    = MAX(1.E-12_ireals , MIN (1.0_ireals , z_iww) )
      z_iwe    = z1dg * (ziwe(1,kspec) + ziwe(2,kspec)/  &
                 (ziwe(3,kspec)*prhoiwc(j1,j2,j3)+ziwe(4,kspec)))
      z_iwe    = MAX(ziwemn(kspec),MIN(ziwemx(kspec),z_iwe  ))
 
      ziwoda(j1,j2) = z_iwe * pduiwc(j1,j2,j3)*(1.0-z_iww)
      ziwods(j1,j2) = z_iwe * pduiwc(j1,j2,j3)*     z_iww  *(1.-z_iwf)
      ziwb0 (j1,j2) = z1d8*(4.+z_iwg)/(1.+z_iwg)
      END DO
    END DO

    IF (ldebug) THEN
      print *,' prholwc (j1b,j2b) :',prholwc (j1b,j2b,j3)
      print *,' pdulwc  (j1b,j2b) :',pdulwc  (j1b,j2b,j3)
      print *,' zlwoda  (j1b,j2b) :',zlwoda  (j1b,j2b)
      print *,' zlwods  (j1b,j2b) :',zlwods  (j1b,j2b)
      print *,' zlwb0   (j1b,j2b) :',zlwb0   (j1b,j2b)
      print *,' z_lwg                 :',z_lwg                
      print *,' z_lwf                 :',z_lwf                
      print *,' z_lww                 :',z_lww                
      print *,' z_lwe                 :',z_lwe                
      print *,' prhoiwc (j1b,j2b) :',prhoiwc (j1b,j2b,j3)
      print *,' pduiwc  (j1b,j2b) :',pduiwc  (j1b,j2b,j3)
      print *,' ziwoda  (j1b,j2b) :',ziwoda  (j1b,j2b)
      print *,' ziwods  (j1b,j2b) :',ziwods  (j1b,j2b)
      print *,' ziwb0   (j1b,j2b) :',ziwb0   (j1b,j2b)
    ENDIF   
 
    ! Aerosoles
 
    DO j2 = ki2sc, ki2ec
      DO j1 = ki1sc, ki1ec
        zaeoda(j1,j2) =  paeq1(j1,j2,j3)*zaea(kspec,1) &
                       + paeq2(j1,j2,j3)*zaea(kspec,2) &
                       + paeq3(j1,j2,j3)*zaea(kspec,3) &
                       + paeq4(j1,j2,j3)*zaea(kspec,4) &
                       + paeq5(j1,j2,j3)*zaea(kspec,5)
        zaeods(j1,j2) =                                                &
                  ( paeq1(j1,j2,j3)*zaes(kspec,1)*(1.-zaef(kspec,1)) ) &
                 +( paeq2(j1,j2,j3)*zaes(kspec,2)*(1.-zaef(kspec,2)) ) &
                 +( paeq3(j1,j2,j3)*zaes(kspec,3)*(1.-zaef(kspec,3)) ) &
                 +( paeq4(j1,j2,j3)*zaes(kspec,4)*(1.-zaef(kspec,4)) ) &
                 +( paeq5(j1,j2,j3)*zaes(kspec,5)*(1.-zaef(kspec,5)) )
        zzg=                                                                   &
            ((paeq1(j1,j2,j3)*zaes(kspec,1)*(1.-zaef(kspec,1)))*zaeg(kspec,1)  &
            +(paeq2(j1,j2,j3)*zaes(kspec,2)*(1.-zaef(kspec,2)))*zaeg(kspec,2)  &
            +(paeq3(j1,j2,j3)*zaes(kspec,3)*(1.-zaef(kspec,3)))*zaeg(kspec,3)  &
            +(paeq4(j1,j2,j3)*zaes(kspec,4)*(1.-zaef(kspec,4)))*zaeg(kspec,4)  &
            +(paeq5(j1,j2,j3)*zaes(kspec,5)*(1.-zaef(kspec,5)))*zaeg(kspec,5)) &
            / MAX( zaeods(j1,j2),zepopd)
        zaeb0 (j1,j2) = z1d8*(4.+zzg)/(1.+zzg)
      ENDDO
    ENDDO

#ifdef COSMOART
  IF(l_cosmo_art) THEN
    IF ((lrad_dust) .AND. (.NOT. lrad_seas) .AND. (.NOT. lrad_aero)) THEN
      DO j2 = ki2sc, ki2ec
        DO j1 = ki1sc, ki1ec
          zaeoda(j1,j2) =   paeq1(j1,j2,j3)*zaea(kspec,1) + tau_abs_dust(j1,j3,kspec)   &
                          + paeq2(j1,j2,j3)*zaea(kspec,2)                         &
                          + paeq3(j1,j2,j3)*zaea(kspec,3)                         &
                          + paeq4(j1,j2,j3)*zaea(kspec,4)                         &
                          + paeq5(j1,j2,j3)*zaea(kspec,5)

          zaeods(j1,j2) = paeq1(j1,j2,j3)*zaes(kspec,1)*(1.-zaef(kspec,1))+       &     
                  tau_streu_dust(j1,j3,kspec)*(1.-(asym_ges(j1,j3,kspec)**2))     &     
                        + paeq2(j1,j2,j3)*zaes(kspec,2)*(1.-zaef(kspec,2))        &
                        + paeq3(j1,j2,j3)*zaes(kspec,3)*(1.-zaef(kspec,3))        &   
                        + paeq4(j1,j2,j3)*zaes(kspec,4)*(1.-zaef(kspec,4))        &   
                        + paeq5(j1,j2,j3)*zaes(kspec,5)*(1.-zaef(kspec,5))

          zzg=( paeq1(j1,j2,j3)*zaes(kspec,1)*(1.-zaef(kspec,1))*zaeg(kspec,1) +  &
                          tau_streu_dust(j1,j3,kspec)*                            &
                          (1.-(asym_ges(j1,j3,kspec)**2))*asym_ges(j1,j3,kspec)   &     
                        + paeq2(j1,j2,j3)*zaes(kspec,2)*(1.-zaef(kspec,2))*zaeg(kspec,2)  &         
                        + paeq3(j1,j2,j3)*zaes(kspec,3)*(1.-zaef(kspec,3))*zaeg(kspec,3)  &         
                        + paeq4(j1,j2,j3)*zaes(kspec,4)*(1.-zaef(kspec,4))*zaeg(kspec,4)  &         
                        + paeq5(j1,j2,j3)*zaes(kspec,5)*(1.-zaef(kspec,5))*zaeg(kspec,5)) &
                          / MAX( zaeods(j1,j2),zepopd)
          zaeb0 (j1,j2) = z1d8*(4.+zzg)/(1.+zzg)
        ENDDO
      ENDDO

    ELSEIF ((lrad_dust) .AND. (.NOT. lrad_seas) .AND. (lrad_aero)) THEN

      DO j2 = ki2sc, ki2ec
        DO j1 = ki1sc, ki1ec
          zaeoda(j1,j2) = paeq1(j1,j2,j3)*zaea(kspec,1) + tau_abs_dust(j1,j3,kspec)    &
                + paeq2(j1,j2,j3)*zaea(kspec,2)                                &
                + paeq3(j1,j2,j3)*zaea(kspec,3) + tau_abs_aero(j1,j3,kspec)    &
                + paeq4(j1,j2,j3)*zaea(kspec,4)                                &
                + paeq5(j1,j2,j3)*zaea(kspec,5)

          zaeods(j1,j2) = paeq1(j1,j2,j3)*zaes(kspec,1)*(1.-zaef(kspec,1)) +     &       
                  tau_streu_dust(j1,j3,kspec)*(1.-(asym_ges(j1,j3,kspec)**2))    &
                 +paeq2(j1,j2,j3)*zaes(kspec,2)*(1.-zaef(kspec,2))               &
                 +paeq3(j1,j2,j3)*zaes(kspec,3)*(1.-zaef(kspec,3)) +             &
                  tau_streu_aero(j1,j3,kspec)*(1.-(asym_aero(j1,j3,kspec)**2))   &
                 +paeq4(j1,j2,j3)*zaes(kspec,4)*(1.-zaef(kspec,4))               &
                 +paeq5(j1,j2,j3)*zaes(kspec,5)*(1.-zaef(kspec,5))

          zzg= ( paeq1(j1,j2,j3)*zaes(kspec,1)*(1.-zaef(kspec,1))*zaeg(kspec,1) +                   & 
               tau_streu_dust(j1,j3,kspec)*(1.-(asym_ges(j1,j3,kspec)**2))*asym_ges(j1,j3,kspec)    &
             + paeq2(j1,j2,j3)*zaes(kspec,2)*(1.-zaef(kspec,2))*zaeg(kspec,2)                       &
             + paeq3(j1,j2,j3)*zaes(kspec,3)*(1.-zaef(kspec,3))*zaeg(kspec,3)  +                    &
               tau_streu_aero(j1,j3,kspec)*(1.-(asym_aero(j1,j3,kspec)**2.))*asym_aero(j1,j3,kspec) &
             + paeq4(j1,j2,j3)*zaes(kspec,4)*(1.-zaef(kspec,4))*zaeg(kspec,4)                       &
             + paeq5(j1,j2,j3)*zaes(kspec,5)*(1.-zaef(kspec,5))*zaeg(kspec,5) )                     &
               / MAX( zaeods(j1,j2),zepopd)
          zaeb0 (j1,j2) = z1d8*(4.+zzg)/(1.+zzg)
        ENDDO
      ENDDO

    ELSEIF ((.NOT. lrad_dust) .AND. (.NOT. lrad_seas) .AND. (lrad_aero)) THEN

      DO j2 = ki2sc, ki2ec
        DO j1 = ki1sc, ki1ec
          zaeoda(j1,j2) = paeq1(j1,j2,j3)*zaea(kspec,1) &
                + paeq2(j1,j2,j3)*zaea(kspec,2) &
                + paeq3(j1,j2,j3)*zaea(kspec,3) + tau_abs_aero(j1,j3,kspec) &
                + paeq4(j1,j2,j3)*zaea(kspec,4) &
                + paeq5(j1,j2,j3)*zaea(kspec,5)

          zaeods(j1,j2) = paeq1(j1,j2,j3)*zaes(kspec,1)*(1.-zaef(kspec,1))              &
                 +paeq2(j1,j2,j3)*zaes(kspec,2)*(1.-zaef(kspec,2))              &
                 +paeq3(j1,j2,j3)*zaes(kspec,3)*(1.-zaef(kspec,3))+             &
                  tau_streu_aero(j1,j3,kspec)*(1.-(asym_aero(j1,j3,kspec)**2.)) &
                 +paeq4(j1,j2,j3)*zaes(kspec,4)*(1.-zaef(kspec,4))              &
                 +paeq5(j1,j2,j3)*zaes(kspec,5)*(1.-zaef(kspec,5))

          zzg=( paeq1(j1,j2,j3)*zaes(kspec,1)*(1.-zaef(kspec,1))*zaeg(kspec,1)                       &
               +paeq2(j1,j2,j3)*zaes(kspec,2)*(1.-zaef(kspec,2))*zaeg(kspec,2)                       &
               +paeq3(j1,j2,j3)*zaes(kspec,3)*(1.-zaef(kspec,3))*zaeg(kspec,3)  +                    &
                tau_streu_aero(j1,j3,kspec)*(1.-(asym_aero(j1,j3,kspec)**2.))*asym_aero(j1,j3,kspec) &
               +paeq4(j1,j2,j3)*zaes(kspec,4)*(1.-zaef(kspec,4))*zaeg(kspec,4)                       &
               +paeq5(j1,j2,j3)*zaes(kspec,5)*(1.-zaef(kspec,5))*zaeg(kspec,5) )                     &
                   / MAX( zaeods(j1,j2),zepopd)
          zaeb0 (j1,j2) = z1d8*(4.+zzg)/(1.+zzg)
        ENDDO
      ENDDO

    ELSEIF ((lrad_dust) .AND. (lrad_seas) .AND. (.NOT.lrad_aero)) THEN

      DO j2 = ki2sc, ki2ec
        DO j1 = ki1sc, ki1ec
          zaeoda(j1,j2) =   paeq1(j1,j2,j3)*zaea(kspec,1) + tau_abs_dust(j1,j3,kspec)   &
                  + tau_abs_seas(j1,j3,kspec)&
                  + paeq3(j1,j2,j3)*zaea(kspec,3) &
                  + paeq4(j1,j2,j3)*zaea(kspec,4) &
                  + paeq5(j1,j2,j3)*zaea(kspec,5)

          zaeods(j1,j2) = paeq1(j1,j2,j3)*zaes(kspec,1)*(1.-zaef(kspec,1))+                  &
                  tau_streu_dust(j1,j3,kspec)*(1.-(asym_ges(j1,j3,kspec)**2))        &
                + tau_streu_seas(j1,j3,kspec)*(1.-(asym_seas(j1,j3,kspec)**2))       &
                + paeq3(j1,j2,j3)*zaes(kspec,3)*(1.-zaef(kspec,3))                   &
                + paeq4(j1,j2,j3)*zaes(kspec,4)*(1.-zaef(kspec,4))                   &
                + paeq5(j1,j2,j3)*zaes(kspec,5)*(1.-zaef(kspec,5))

          zzg=( paeq1(j1,j2,j3)*zaes(kspec,1)*(1.-zaef(kspec,1))*zaeg(kspec,1) +                          &
                tau_streu_dust(j1,j3,kspec)*(1.-(asym_ges(j1,j3,kspec)**2))*asym_ges(j1,j3,kspec)         &
               +tau_streu_seas(j1,j3,kspec)*(1.-(asym_seas(j1,j3,kspec)**2))*asym_seas(j1,j3,kspec)       &
               +paeq3(j1,j2,j3)*zaes(kspec,3)*(1.-zaef(kspec,3))*zaeg(kspec,3)                            &
               +paeq4(j1,j2,j3)*zaes(kspec,4)*(1.-zaef(kspec,4))*zaeg(kspec,4)                            &
               +paeq5(j1,j2,j3)*zaes(kspec,5)*(1.-zaef(kspec,5))*zaeg(kspec,5) )                          &
             / MAX( zaeods(j1,j2),zepopd)
          zaeb0 (j1,j2) = z1d8*(4.+zzg)/(1.+zzg)
        ENDDO
      ENDDO

    ELSEIF ((.NOT. lrad_dust) .and. (lrad_seas) .and. (.NOT. lrad_aero)) THEN

      DO j2 = ki2sc, ki2ec
        DO j1 = ki1sc, ki1ec
          zaeoda(j1,j2) = paeq1(j1,j2,j3)*zaea(kspec,1)                                &
                 + tau_abs_seas(j1,j3,kspec)    &
                 + paeq3(j1,j2,j3)*zaea(kspec,3)                                &
                 + paeq4(j1,j2,j3)*zaea(kspec,4)                                &
                 + paeq5(j1,j2,j3)*zaea(kspec,5)

          zaeods(j1,j2) =  paeq1(j1,j2,j3)*zaes(kspec,1)*(1.-zaef(kspec,1))            &
                  + tau_streu_seas(j1,j3,kspec)*(1.-(asym_seas(j1,j3,kspec)**2))&
                  + paeq3(j1,j2,j3)*zaes(kspec,3)*(1.-zaef(kspec,3))            &
                  + paeq4(j1,j2,j3)*zaes(kspec,4)*(1.-zaef(kspec,4))            &
                  + paeq5(j1,j2,j3)*zaes(kspec,5)*(1.-zaef(kspec,5))

          zzg=   ( paeq1(j1,j2,j3)*zaes(kspec,1)*(1.-zaef(kspec,1))*zaeg(kspec,1)                       &
           +tau_streu_seas(j1,j3,kspec)*(1.-(asym_seas(j1,j3,kspec)**2))*asym_seas(j1,j3,kspec)  &
           +paeq3(j1,j2,j3)*zaes(kspec,3)*(1.-zaef(kspec,3))*zaeg(kspec,3)                       &
           +paeq4(j1,j2,j3)*zaes(kspec,4)*(1.-zaef(kspec,4))*zaeg(kspec,4)                       &
           +paeq5(j1,j2,j3)*zaes(kspec,5)*(1.-zaef(kspec,5))*zaeg(kspec,5) )                     &
           / MAX( zaeods(j1,j2),zepopd)
          zaeb0 (j1,j2) = z1d8*(4.+zzg)/(1.+zzg)
        ENDDO
      ENDDO

    ELSEIF ((lrad_dust) .AND. (lrad_seas) .AND. (lrad_aero)) THEN

      DO j2 = ki2sc, ki2ec
        DO j1 = ki1sc, ki1ec
          zaeoda(j1,j2) = paeq1(j1,j2,j3)*zaea(kspec,1) + tau_abs_dust(j1,j3,kspec)    &
                + tau_abs_seas(j1,j3,kspec)    &
                + paeq3(j1,j2,j3)*zaea(kspec,3) + tau_abs_aero(j1,j3,kspec)    &
                + paeq4(j1,j2,j3)*zaea(kspec,4)                                &
                + paeq5(j1,j2,j3)*zaea(kspec,5)

          zaeods(j1,j2) = paeq1(j1,j2,j3)*zaes(kspec,1)*(1.-zaef(kspec,1)) +      &      
                  tau_streu_dust(j1,j3,kspec)*(1.-(asym_ges(j1,j3,kspec)**2))    &
                 +tau_streu_seas(j1,j3,kspec)*(1.-(asym_seas(j1,j3,kspec)**2))   &
                 +paeq3(j1,j2,j3)*zaes(kspec,3)*(1.-zaef(kspec,3)) +             &
                  tau_streu_aero(j1,j3,kspec)*(1.-(asym_aero(j1,j3,kspec)**2))   &
                 +paeq4(j1,j2,j3)*zaes(kspec,4)*(1.-zaef(kspec,4))               &
                 +paeq5(j1,j2,j3)*zaes(kspec,5)*(1.-zaef(kspec,5))

          zzg= ( paeq1(j1,j2,j3)*zaes(kspec,1)*(1.-zaef(kspec,1))*zaeg(kspec,1) +                     &
               tau_streu_dust(j1,j3,kspec)*(1.-(asym_ges(j1,j3,kspec)**2))*asym_ges(j1,j3,kspec)    &
              +tau_streu_seas(j1,j3,kspec)*(1.-(asym_seas(j1,j3,kspec)**2))*asym_seas(j1,j3,kspec)  &
              +paeq3(j1,j2,j3)*zaes(kspec,3)*(1.-zaef(kspec,3))*zaeg(kspec,3)  +                    &
               tau_streu_aero(j1,j3,kspec)*(1.-(asym_aero(j1,j3,kspec)**2.))*asym_aero(j1,j3,kspec) &
              +paeq4(j1,j2,j3)*zaes(kspec,4)*(1.-zaef(kspec,4))*zaeg(kspec,4)                       &
              +paeq5(j1,j2,j3)*zaes(kspec,5)*(1.-zaef(kspec,5))*zaeg(kspec,5) )                     &
             / MAX( zaeods(j1,j2),zepopd)
          zaeb0 (j1,j2) = z1d8*(4.+zzg)/(1.+zzg)
        ENDDO
      ENDDO

    ELSEIF ((.NOT. lrad_dust) .AND. (lrad_seas) .AND. (lrad_aero)) THEN

      DO j2 = ki2sc, ki2ec
        DO j1 = ki1sc, ki1ec
          zaeoda(j1,j2) = paeq1(j1,j2,j3)*zaea(kspec,1) &
                + tau_abs_seas(j1,j3,kspec) &
                + paeq3(j1,j2,j3)*zaea(kspec,3) + tau_abs_aero(j1,j3,kspec) &
                + paeq4(j1,j2,j3)*zaea(kspec,4) &
                + paeq5(j1,j2,j3)*zaea(kspec,5)

          zaeods(j1,j2) = paeq1(j1,j2,j3)*zaes(kspec,1)*(1.-zaef(kspec,1))               &
                 +tau_streu_seas(j1,j3,kspec)*(1.-(asym_seas(j1,j3,kspec)**2.))  &
                 +paeq3(j1,j2,j3)*zaes(kspec,3)*(1.-zaef(kspec,3)) +             &
                  tau_streu_aero(j1,j3,kspec)*(1.-(asym_aero(j1,j3,kspec)**2.))  &
                 +paeq4(j1,j2,j3)*zaes(kspec,4)*(1.-zaef(kspec,4))               &
                 +paeq5(j1,j2,j3)*zaes(kspec,5)*(1.-zaef(kspec,5))

          zzg=( paeq1(j1,j2,j3)*zaes(kspec,1)*(1.-zaef(kspec,1))*zaeg(kspec,1)                       &
             +tau_streu_seas(j1,j3,kspec)*(1.-(asym_seas(j1,j3,kspec)**2.))*asym_seas(j1,j3,kspec) &
             +paeq3(j1,j2,j3)*zaes(kspec,3)*(1.-zaef(kspec,3))*zaeg(kspec,3) +                     &
              tau_streu_aero(j1,j3,kspec)*(1.-(asym_aero(j1,j3,kspec)**2.))*asym_aero(j1,j3,kspec) &
             +paeq4(j1,j2,j3)*zaes(kspec,4)*(1.-zaef(kspec,4))*zaeg(kspec,4)                       &
             +paeq5(j1,j2,j3)*zaes(kspec,5)*(1.-zaef(kspec,5))*zaeg(kspec,5) )                     &
                 / MAX( zaeods(j1,j2),zepopd)
          zaeb0 (j1,j2) = z1d8*(4.+zzg)/(1.+zzg)
        ENDDO
      ENDDO

    ENDIF
  ENDIF
#endif

    IF (ldebug) THEN
      print *,' zaeoda  (j1b,j2b) :',zaeoda  (j1b,j2b)
      print *,' zaeods  (j1b,j2b) :',zaeods  (j1b,j2b)
      print *,' zaeb0   (j1b,j2b) :',zaeb0   (j1b,j2b)
    ENDIF    
 
    ! Combined effects
 
    ! a) cloud free part of each layer

    DO j2 = ki2sc, ki2ec
      DO j1 = ki1sc, ki1ec
        podaf(j1,j2,j3) = MAX(zaeoda(j1,j2), zepopd)   
        podsf(j1,j2,j3) = MAX(zaeods(j1,j2), zepopd)
        pbsff(j1,j2,j3) =     zaeb0 (j1,j2)
      ENDDO
    ENDDO

    IF (ldebug) THEN
      print *,' podaf   (j1b,j2b,j3) :',podaf   (j1b,j2b,j3)
      print *,' podsf   (j1b,j2b,j3) :',podsf   (j1b,j2b,j3)
      print *,' pbsff   (j1b,j2b,j3) :',pbsff   (j1b,j2b,j3)
    ENDIF
 
    ! b) cloudy part of each layer

    DO j2 = ki2sc, ki2ec
      DO j1 = ki1sc, ki1ec
        podac(j1,j2,j3) = MAX( zlwoda (j1,j2) + ziwoda (j1,j2) + zaeoda(j1,j2),&
                               zepopd)
        podsc(j1,j2,j3) = MAX( zlwods (j1,j2) + ziwods (j1,j2) + zaeods(j1,j2),&
                               zepopd)
        podsc(j1,j2,j3) = MIN( podsc(j1,j2,j3),                                &
                   (1.0_ireals-zepssa) * (podac(j1,j2,j3)+podsc(j1,j2,j3)))
        pbsfc(j1,j2,j3) = (  zlwb0 (j1,j2)*zlwods (j1,j2)                      &
                           + ziwb0 (j1,j2)*ziwods (j1,j2)                      &
                           + zaeb0 (j1,j2)*zaeods (j1,j2) ) / podsc(j1,j2,j3)
      ENDDO
    ENDDO

    IF (ldebug) THEN
      print *,' podac   (j1b,j2b,j3) :',podac   (j1b,j2b,j3)
      print *,' podsc   (j1b,j2b,j3) :',podsc   (j1b,j2b,j3)
      print *,' pbsfc   (j1b,j2b,j3) :',pbsfc   (j1b,j2b,j3)
    ENDIF 
 
  ! End of vertical loop
  ! --------------------
  ENDDO       
 
!------------------------------------------------------------------------------
! End of the subroutine
!------------------------------------------------------------------------------

END SUBROUTINE opt_th

!==============================================================================
!==============================================================================
!+ Module procedure in "Radiation" 
!------------------------------------------------------------------------------

SUBROUTINE opt_so(prholwc,pdulwc,prhoiwc,pduiwc,                &
                  paeq1  ,paeq2 ,paeq3  ,paeq4 , paeq5 ,        &
                  pdp    ,papra ,psmu0  ,pqsmu0,                &
                  ki1sd  ,ki1ed ,ki2sd  ,ki2ed , ki3sd , ki3ed, &
                  kspec  ,ki1sc ,ki1ec  ,ki2sc , ki2ec ,        &
                          ki3sc ,ki3ec  ,ldebug,                &
                  podac  ,podaf ,podsc  ,podsf , pbsfc ,pbsff , &
                  pusfc  ,pusff  )

!------------------------------------------------------------------------------
!
! Description:
!
!   The module procedure opt_so calculates the optical properties of the 
!   non-gaseous constituents for one spectral interval in the solar part 
!   of the spectrum.
!   Purpose is the calculation of (absorption and scattering) optical
!   depth and backward scattered fraction for diffuse and upward scattered 
!   fraction for direct solar radiation (excluding the contribution by 
!   gaseous constituents).
!
! Method:
!
! - determination of optical properies (i.e. extinction coefficient,
!   single scattering albedo and asymetry factor of the phase function)
!   using approximate relations for Rayleigh effects, cloud water,
!   cloud ice and a combination of five type of aerosols
!
! - calculation of optical depth (scattering and absorption), back-
!   scattered and upscattered fraction of radiation suitable for use
!   in an implicit delta-two-stream scheme
!     
!------------------------------------------------------------------------------

! Subroutine arguments:
! --------------------

! Input data
! ----------
  INTEGER (KIND=iintegers), INTENT (IN) ::  &
     ki1sd,       & ! start index for first  array dimension
     ki1ed,       & ! end   index for first  array dimension
     ki2sd,       & ! start index for second array dimension
     ki2ed,       & ! end   index for second array dimension
     ki3sd,       & ! start index for third  array dimension
     ki3ed,       & ! end   index for third  array dimension

   ! and the same for the computations
     ki1sc,       & ! start index for first  array computation
     ki1ec,       & ! end   index for first  array computation
     ki2sc,       & ! start index for second array computation
     ki2ec,       & ! end   index for second array computation
     ki3sc,       & ! start index for third  array computation
     ki3ec,       & ! end   index for third  array computation
     kspec          ! selected spectral interval

  LOGICAL                 , INTENT (IN) ::  &
     ldebug         ! debug control switch       

  REAL    (KIND=ireals   ), INTENT (IN) ::  &
           !  Liquid and ice water density and content within for the cloudy
           !  part of each layer
     prholwc(ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed), &
     prhoiwc(ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed), &
     pdulwc (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed), &
     pduiwc (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed), &

           !  Aerosole contents (optical depths at 0.55 micrometer) for 5 types
     paeq1  (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed), &
     paeq2  (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed), &
     paeq3  (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed), &
     paeq4  (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed), &
     paeq5  (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed), &

           !  pressure thickness of layers 
     pdp    (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed), &
 
           !  zenith angle and it's inverse
     psmu0  (ki1sd:ki1ed,ki2sd:ki2ed)            , &
     pqsmu0 (ki1sd:ki1ed,ki2sd:ki2ed) 

  REAL    (KIND=ireals   ), INTENT (INOUT) ::  &
           ! mean layer pressure (TOA on input)
     papra  (ki1sd:ki1ed,ki2sd:ki2ed)                
  
! Output data
! -----------
  REAL    (KIND=ireals   ), INTENT (OUT) ::  &
     podac (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed), & ! absorption optical depth
     podaf (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed), & ! in cloudy and free part
     podsc (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed), & ! scattering optical depth
     podsf (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed), & ! in cloudy and free part
     pbsfc (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed), & ! backscattering fraction 
     pbsff (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed), & ! in cloudy and free part
     pusfc (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed), & ! upward scattered fraction
     pusff (ki1sd:ki1ed,ki2sd:ki2ed,ki3sd:ki3ed)    ! in cloudy and free part

! Local parameters: 
! ----------------
  REAL    (KIND=ireals   ), PARAMETER ::  &
     z1dg   = 1.0/9.80665, & ! 1./g
     z1d8   = 0.125      , & ! 1./8 
     zepopd = 1.0E-6     , & ! Security constant for optical depth
     zepssa = 1.0E-6         ! Security constant for single scattering albedo

  INTEGER (KIND=iintegers), PARAMETER ::  &
     j1b    = 1,           & ! debug point index first dimension
     j2b    = 1              ! debug point index second dimension

! Local scalars:
! -------------
  INTEGER (KIND=iintegers) ::  &
    j1,j2,j3,              & ! loop indices over spatial dimensions
    ja      ,              & ! local loop index
    isp                      ! (=kspec, but shorter for notation purposes)

  REAL    (KIND=ireals   ) ::  &
    ! individual optical properties of liquid water and ice 
    z_lwe, z_iwe,          & ! extinction coefficient
    z_lww, z_iww,          & ! single scattering coefficient
    z_lwg, z_iwg,          & ! asymetry factor
    z_lwf, z_iwf,          & ! forward scattered fraction 
    zzg
 
! Local (automatic) arrays:
! ------------------------
!    optical properties (absorption, scattering, back- and upscatter fraction)
!    for liquid water, ice and total aerosole
!    for Rayleigh scattering, only the optical depth is stored as array, since
!    back- and upscatter fractions are constant (=0.5)

  REAL    (KIND=ireals   ) ::  &
     zlwoda (ki1sd:ki1ed,ki2sd:ki2ed), & !
     zlwods (ki1sd:ki1ed,ki2sd:ki2ed), & ! 
     zlwb0  (ki1sd:ki1ed,ki2sd:ki2ed), & ! 
     zlwb   (ki1sd:ki1ed,ki2sd:ki2ed), & ! 
     ziwoda (ki1sd:ki1ed,ki2sd:ki2ed), & ! 
     ziwods (ki1sd:ki1ed,ki2sd:ki2ed), & ! 
     ziwb0  (ki1sd:ki1ed,ki2sd:ki2ed), & ! 
     ziwb   (ki1sd:ki1ed,ki2sd:ki2ed), & ! 
     zaeoda (ki1sd:ki1ed,ki2sd:ki2ed), & ! 
     zaeods (ki1sd:ki1ed,ki2sd:ki2ed), & ! 
     zaeb0  (ki1sd:ki1ed,ki2sd:ki2ed), & ! 
     zaeb   (ki1sd:ki1ed,ki2sd:ki2ed), & ! 
     zraods (ki1sd:ki1ed,ki2sd:ki2ed)    ! 
 
!- End of header
!==============================================================================
 
!------------------------------------------------------------------------------
! Begin Subroutine opt_so              
!------------------------------------------------------------------------------
 

  isp      = kspec
  IF (ldebug) THEN
      print *,' **** opt-so   start **********************'
      print *,' **** debug point : ',j1b,j2b
      print *,' **** interval    : ',isp    
  ENDIF 
 
  DO ja=1,5
    zaef(isp,ja)  = zaeg(isp,ja)**2 ! forward sc.fraction f.aerosols
  ENDDO

  IF (ldebug) THEN
     DO ja=1,5
       print *,'ja, zaef(isp,ja): ',ja,zaef(isp,ja)
     ENDDO 
  ENDIF 

  IF (ldebug) print *,' In opt-so   vor vertical loop'
 
! Vertical loop
! ------------ 

  DO j3=ki3sc,ki3ec
 
    if (ldebug) print *,' In opt-so   j3 = ',j3  

    !     Optical properties of liquid water and ice as function of the specific
    !     liquid water and ice content
 
    DO j2=ki2sc,ki2ec
      DO j1=ki1sc,ki1ec
 
        ! liquid water effects
 
        z_lwg = zlwg(1,isp) + zlwg(2,isp)*prholwc(j1,j2,j3)
        z_lwg = MAX(0.0_ireals,MIN(1.0_ireals,z_lwg))
        z_lwf = z_lwg*z_lwg 
        z_lww = zlww(1,isp) + zlww(2,isp)*prholwc(j1,j2,j3)
        z_lww = MAX(zepssa,MIN(1.0-zepssa,z_lww))
        z_lwe = z1dg * (zlwe(1,isp) + zlwe(2,isp)/ &
                (zlwe(3,isp)*prholwc(j1,j2,j3)+zlwe(4,isp)))
        z_lwe = MAX(zlwemn(isp),MIN(zlwemx(isp),z_lwe))
 
        zlwoda(j1,j2) = z_lwe*pdulwc(j1,j2,j3)*(1.0-z_lww)
        zlwods(j1,j2) = z_lwe*pdulwc(j1,j2,j3)*     z_lww *(1.-z_lwf)
        zlwb0 (j1,j2) = z1d8*(4.+z_lwg)/(1.+z_lwg) 
        zlwb  (j1,j2) = 0.5-0.75*psmu0(j1,j2)*z_lwg/(1.+z_lwg)
 
        ! ice water effects
 
        z_iwg = ziwg(1,isp) + ziwg(2,isp)* LOG(prhoiwc(j1,j2,j3))
        z_iwg = MAX(0.0_ireals,MIN(1.0_ireals,z_iwg))
        z_iwf = z_iwg*z_iwg
        z_iww = ziww(1,isp) + ziww(2,isp)* LOG(prhoiwc(j1,j2,j3))
        z_iww = MAX(zepssa,MIN(1.0_ireals-zepssa,z_iww))
        z_iwe = z1dg * (ziwe(1,isp) + ziwe(2,isp)/ &
                   (ziwe(3,isp)*prhoiwc(j1,j2,j3)+ziwe(4,isp)))
        z_iwe = MAX(ziwemn(isp),MIN(ziwemx(isp),z_iwe  ))
 
        ziwoda(j1,j2) = z_iwe*pduiwc(j1,j2,j3)*(1.0-z_iww)
        ziwods(j1,j2) = z_iwe*pduiwc(j1,j2,j3)*     z_iww *(1.-z_iwf)
        ziwb0 (j1,j2) = z1d8*(4.+z_iwg)/(1.+z_iwg)
        ziwb  (j1,j2) = 0.5-0.75*psmu0(j1,j2)*z_iwg/(1.+z_iwg)
      END DO
    END DO

    IF (ldebug) THEN
      print *,' prholwc (j1b,j2b) :',prholwc (j1b,j2b,j3)   
      print *,' pdulwc  (j1b,j2b) :',pdulwc  (j1b,j2b,j3)
      print *,' zlwoda  (j1b,j2b) :',zlwoda  (j1b,j2b)
      print *,' zlwods  (j1b,j2b) :',zlwods  (j1b,j2b)
      print *,' zlwb0   (j1b,j2b) :',zlwb0   (j1b,j2b)
      print *,' zlwb    (j1b,j2b) :',zlwb    (j1b,j2b)
      print *,' z_lwg             :',z_lwg                
      print *,' z_lwf             :',z_lwf                
      print *,' z_lww             :',z_lww                
      print *,' z_lwe             :',z_lwe                
      print *,' prhoiwc (j1b,j2b) :',prhoiwc (j1b,j2b,j3)
      print *,' pduiwc  (j1b,j2b) :',pduiwc  (j1b,j2b,j3)
      print *,' ziwoda  (j1b,j2b) :',ziwoda  (j1b,j2b)
      print *,' ziwods  (j1b,j2b) :',ziwods  (j1b,j2b)
      print *,' ziwb0   (j1b,j2b) :',ziwb0   (j1b,j2b)
      print *,' ziwb    (j1b,j2b) :',ziwb    (j1b,j2b)
    ENDIF 

    ! Optical properties of five aerosol types combined
 
    DO j2=ki2sc,ki2ec
      DO j1=ki1sc,ki1ec

        zaeoda(j1,j2)=                                             &
          paeq1(j1,j2,j3)*zaea(isp,1)+paeq2(j1,j2,j3)*zaea(isp,2)+ &
          paeq3(j1,j2,j3)*zaea(isp,3)+paeq4(j1,j2,j3)*zaea(isp,4)+ &
          paeq5(j1,j2,j3)*zaea(isp,5)

        zaeods(j1,j2)=                                        &
            ( paeq1(j1,j2,j3)*zaes(isp,1)*(1.-zaef(isp,1)) )  &
           +( paeq2(j1,j2,j3)*zaes(isp,2)*(1.-zaef(isp,2)) )  &
           +( paeq3(j1,j2,j3)*zaes(isp,3)*(1.-zaef(isp,3)) )  &
           +( paeq4(j1,j2,j3)*zaes(isp,4)*(1.-zaef(isp,4)) )  &
           +( paeq5(j1,j2,j3)*zaes(isp,5)*(1.-zaef(isp,5)) )
 
        zzg=(                                                          &
           (paeq1(j1,j2,j3)*zaes(isp,1)*(1.-zaef(isp,1)))*zaeg(isp,1)  &
          +(paeq2(j1,j2,j3)*zaes(isp,2)*(1.-zaef(isp,2)))*zaeg(isp,2)  &
          +(paeq3(j1,j2,j3)*zaes(isp,3)*(1.-zaef(isp,3)))*zaeg(isp,3)  &
          +(paeq4(j1,j2,j3)*zaes(isp,4)*(1.-zaef(isp,4)))*zaeg(isp,4)  &
          +(paeq5(j1,j2,j3)*zaes(isp,5)*(1.-zaef(isp,5)))*zaeg(isp,5)) &
          / MAX( zaeods(j1,j2),zepopd)

        zaeb0(j1,j2) = z1d8*(4.+zzg)/(1.+zzg)
        zaeb (j1,j2) = 0.5-0.75*psmu0(j1,j2)*zzg/(1.+zzg)

      ENDDO
    ENDDO

#ifdef COSMOART
  IF(l_cosmo_art) THEN
    IF ((lrad_dust) .AND. (.NOT.lrad_seas) .AND. (.NOT. lrad_aero)) THEN    

      DO j2 = ki2sc, ki2ec
        DO j1 = ki1sc, ki1ec
          zaeoda(j1,j2) =   paeq1(j1,j2,j3)*zaea(kspec,1) + tau_abs_dust(j1,j3,kspec)   &
                  + paeq2(j1,j2,j3)*zaea(kspec,2) &
                  + paeq3(j1,j2,j3)*zaea(kspec,3) &
                  + paeq4(j1,j2,j3)*zaea(kspec,4) &
                  + paeq5(j1,j2,j3)*zaea(kspec,5)

          zaeods(j1,j2) = paeq1(j1,j2,j3)*zaes(kspec,1)*(1.-zaef(kspec,1))+             &
                  tau_streu_dust(j1,j3,kspec)*(1.-(asym_ges(j1,j3,kspec)**2))        &
                + paeq2(j1,j2,j3)*zaes(kspec,2)*(1.-zaef(kspec,2))              &
                + paeq3(j1,j2,j3)*zaes(kspec,3)*(1.-zaef(kspec,3))                   &
                + paeq4(j1,j2,j3)*zaes(kspec,4)*(1.-zaef(kspec,4))                   &
                + paeq5(j1,j2,j3)*zaes(kspec,5)*(1.-zaef(kspec,5))

          zzg=( paeq1(j1,j2,j3)*zaes(kspec,1)*(1.-zaef(kspec,1))*zaeg(kspec,1) +                          &
                tau_streu_dust(j1,j3,kspec)*(1.-(asym_ges(j1,j3,kspec)**2))*asym_ges(j1,j3,kspec)         &
               +paeq2(j1,j2,j3)*zaes(kspec,2)*(1.-zaef(kspec,2))*zaeg(kspec,2)                            &
               +paeq3(j1,j2,j3)*zaes(kspec,3)*(1.-zaef(kspec,3))*zaeg(kspec,3)                            &
               +paeq4(j1,j2,j3)*zaes(kspec,4)*(1.-zaef(kspec,4))*zaeg(kspec,4)                            &
               +paeq5(j1,j2,j3)*zaes(kspec,5)*(1.-zaef(kspec,5))*zaeg(kspec,5) )                          &
             / MAX( zaeods(j1,j2),zepopd)
          zaeb0 (j1,j2) = z1d8*(4.+zzg)/(1.+zzg)
        ENDDO
      ENDDO

    ELSEIF ((lrad_dust) .AND. (.NOT. lrad_seas) .AND. (lrad_aero)) THEN

      DO j2 = ki2sc, ki2ec
        DO j1 = ki1sc, ki1ec
          zaeoda(j1,j2) = paeq1(j1,j2,j3)*zaea(kspec,1) + tau_abs_dust(j1,j3,kspec)    &
                + paeq2(j1,j2,j3)*zaea(kspec,2)                                &
                + paeq3(j1,j2,j3)*zaea(kspec,3) + tau_abs_aero(j1,j3,kspec)    &
                + paeq4(j1,j2,j3)*zaea(kspec,4)                                &
                + paeq5(j1,j2,j3)*zaea(kspec,5)

          zaeods(j1,j2) = paeq1(j1,j2,j3)*zaes(kspec,1)*(1.-zaef(kspec,1)) +             &
                  tau_streu_dust(j1,j3,kspec)*(1.-(asym_ges(j1,j3,kspec)**2))    &
                + paeq2(j1,j2,j3)*zaes(kspec,2)*(1.-zaef(kspec,2))               &
                + paeq3(j1,j2,j3)*zaes(kspec,3)*(1.-zaef(kspec,3)) +             &
                  tau_streu_aero(j1,j3,kspec)*(1.-(asym_aero(j1,j3,kspec)**2))   &
                + paeq4(j1,j2,j3)*zaes(kspec,4)*(1.-zaef(kspec,4))               &
                + paeq5(j1,j2,j3)*zaes(kspec,5)*(1.-zaef(kspec,5))

          zzg= (paeq1(j1,j2,j3)*zaes(kspec,1)*(1.-zaef(kspec,1))*zaeg(kspec,1) +                       &
                  tau_streu_dust(j1,j3,kspec)*(1.-(asym_ges(j1,j3,kspec)**2))*asym_ges(j1,j3,kspec)    &
                + paeq2(j1,j2,j3)*zaes(kspec,2)*(1.-zaef(kspec,2))*zaeg(kspec,2)                       &
                + paeq3(j1,j2,j3)*zaes(kspec,3)*(1.-zaef(kspec,3))*zaeg(kspec,3)  +                    &
                  tau_streu_aero(j1,j3,kspec)*(1.-(asym_aero(j1,j3,kspec)**2.))*asym_aero(j1,j3,kspec) &
                + paeq4(j1,j2,j3)*zaes(kspec,4)*(1.-zaef(kspec,4))*zaeg(kspec,4)                       &
                + paeq5(j1,j2,j3)*zaes(kspec,5)*(1.-zaef(kspec,5))*zaeg(kspec,5) )                     &
                / MAX( zaeods(j1,j2),zepopd)
          zaeb0 (j1,j2) = z1d8*(4.+zzg)/(1.+zzg)
        ENDDO
      ENDDO

    ELSEIF ((.NOT. lrad_dust) .AND. (.NOT. lrad_seas) .AND. (lrad_aero)) THEN

      DO j2 = ki2sc, ki2ec
        DO j1 = ki1sc, ki1ec
          zaeoda(j1,j2) = paeq1(j1,j2,j3)*zaea(kspec,1) &
                + paeq2(j1,j2,j3)*zaea(kspec,2) &
                + paeq3(j1,j2,j3)*zaea(kspec,3) + tau_abs_aero(j1,j3,kspec) &
                + paeq4(j1,j2,j3)*zaea(kspec,4) &
                + paeq5(j1,j2,j3)*zaea(kspec,5)

          zaeods(j1,j2) = paeq1(j1,j2,j3)*zaes(kspec,1)*(1.-zaef(kspec,1))              &
                 +paeq2(j1,j2,j3)*zaes(kspec,2)*(1.-zaef(kspec,2))              &
                 +paeq3(j1,j2,j3)*zaes(kspec,3)*(1.-zaef(kspec,3))+             &
                  tau_streu_aero(j1,j3,kspec)*(1.-(asym_aero(j1,j3,kspec)**2.)) &
                 +paeq4(j1,j2,j3)*zaes(kspec,4)*(1.-zaef(kspec,4))              &
                 +paeq5(j1,j2,j3)*zaes(kspec,5)*(1.-zaef(kspec,5))

          zzg=( paeq1(j1,j2,j3)*zaes(kspec,1)*(1.-zaef(kspec,1))*zaeg(kspec,1)                       &
               +paeq2(j1,j2,j3)*zaes(kspec,2)*(1.-zaef(kspec,2))*zaeg(kspec,2)                       &
               +paeq3(j1,j2,j3)*zaes(kspec,3)*(1.-zaef(kspec,3))*zaeg(kspec,3)  +                    &
                tau_streu_aero(j1,j3,kspec)*(1.-(asym_aero(j1,j3,kspec)**2.))*asym_aero(j1,j3,kspec) &
               +paeq4(j1,j2,j3)*zaes(kspec,4)*(1.-zaef(kspec,4))*zaeg(kspec,4)                       &
               +paeq5(j1,j2,j3)*zaes(kspec,5)*(1.-zaef(kspec,5))*zaeg(kspec,5) )                     &
                   / MAX( zaeods(j1,j2),zepopd)
          zaeb0 (j1,j2) = z1d8*(4.+zzg)/(1.+zzg)
        ENDDO
      ENDDO
   
    ELSEIF ((lrad_dust) .AND. (lrad_seas) .AND. (.NOT. lrad_aero)) THEN   

      DO j2 = ki2sc, ki2ec
        DO j1 = ki1sc, ki1ec
          zaeoda(j1,j2) =   paeq1(j1,j2,j3)*zaea(kspec,1) + tau_abs_dust(j1,j3,kspec)   &
                  + tau_abs_seas(j1,j3,kspec)&
                  + paeq3(j1,j2,j3)*zaea(kspec,3) &
                  + paeq4(j1,j2,j3)*zaea(kspec,4) &
                  + paeq5(j1,j2,j3)*zaea(kspec,5)

          zaeods(j1,j2) = paeq1(j1,j2,j3)*zaes(kspec,1)*(1.-zaef(kspec,1))+                  &
                  tau_streu_dust(j1,j3,kspec)*(1.-(asym_ges(j1,j3,kspec)**2))        &
                + tau_streu_seas(j1,j3,kspec)*(1.-(asym_seas(j1,j3,kspec)**2))       &
                + paeq3(j1,j2,j3)*zaes(kspec,3)*(1.-zaef(kspec,3))                   &
                + paeq4(j1,j2,j3)*zaes(kspec,4)*(1.-zaef(kspec,4))                   &
                + paeq5(j1,j2,j3)*zaes(kspec,5)*(1.-zaef(kspec,5))

          zzg=( paeq1(j1,j2,j3)*zaes(kspec,1)*(1.-zaef(kspec,1))*zaeg(kspec,1) +                          &
                tau_streu_dust(j1,j3,kspec)*(1.-(asym_ges(j1,j3,kspec)**2))*asym_ges(j1,j3,kspec)         &
               +tau_streu_seas(j1,j3,kspec)*(1.-(asym_seas(j1,j3,kspec)**2))*asym_seas(j1,j3,kspec)       &
               +paeq3(j1,j2,j3)*zaes(kspec,3)*(1.-zaef(kspec,3))*zaeg(kspec,3)                            &
               +paeq4(j1,j2,j3)*zaes(kspec,4)*(1.-zaef(kspec,4))*zaeg(kspec,4)                            &
               +paeq5(j1,j2,j3)*zaes(kspec,5)*(1.-zaef(kspec,5))*zaeg(kspec,5) )                          &
             / MAX( zaeods(j1,j2),zepopd)
          zaeb0 (j1,j2) = z1d8*(4.+zzg)/(1.+zzg)
        ENDDO
      ENDDO

    ELSEIF ((.NOT.lrad_dust) .AND. (lrad_seas) .AND. (.NOT. lrad_aero)) THEN

      DO j2 = ki2sc, ki2ec
        DO j1 = ki1sc, ki1ec
          zaeoda(j1,j2) = paeq1(j1,j2,j3)*zaea(kspec,1)                                &
                 + tau_abs_seas(j1,j3,kspec)    &
                 + paeq3(j1,j2,j3)*zaea(kspec,3)                                &
                 + paeq4(j1,j2,j3)*zaea(kspec,4)                                &
                 + paeq5(j1,j2,j3)*zaea(kspec,5)

          zaeods(j1,j2) =  paeq1(j1,j2,j3)*zaes(kspec,1)*(1.-zaef(kspec,1))            &
                  + tau_streu_seas(j1,j3,kspec)*(1.-(asym_seas(j1,j3,kspec)**2))&
                  + paeq3(j1,j2,j3)*zaes(kspec,3)*(1.-zaef(kspec,3))            &
                  + paeq4(j1,j2,j3)*zaes(kspec,4)*(1.-zaef(kspec,4))            &
                  + paeq5(j1,j2,j3)*zaes(kspec,5)*(1.-zaef(kspec,5))

          zzg=   ( paeq1(j1,j2,j3)*zaes(kspec,1)*(1.-zaef(kspec,1))*zaeg(kspec,1)                       &
                  +tau_streu_seas(j1,j3,kspec)*(1.-(asym_seas(j1,j3,kspec)**2))*asym_seas(j1,j3,kspec)  &
                  +paeq3(j1,j2,j3)*zaes(kspec,3)*(1.-zaef(kspec,3))*zaeg(kspec,3)                       &
                  +paeq4(j1,j2,j3)*zaes(kspec,4)*(1.-zaef(kspec,4))*zaeg(kspec,4)                       &
                  +paeq5(j1,j2,j3)*zaes(kspec,5)*(1.-zaef(kspec,5))*zaeg(kspec,5) )                     &
                  / MAX( zaeods(j1,j2),zepopd)
          zaeb0 (j1,j2) = z1d8*(4.+zzg)/(1.+zzg)
        ENDDO
      ENDDO

    ELSEIF ((lrad_dust) .AND. (lrad_seas) .AND. (lrad_aero)) THEN

      DO j2 = ki2sc, ki2ec
        DO j1 = ki1sc, ki1ec
          zaeoda(j1,j2) = paeq1(j1,j2,j3)*zaea(kspec,1) + tau_abs_dust(j1,j3,kspec)    &
                + tau_abs_seas(j1,j3,kspec)    &
                + paeq3(j1,j2,j3)*zaea(kspec,3) + tau_abs_aero(j1,j3,kspec)    &
                + paeq4(j1,j2,j3)*zaea(kspec,4)                                &
                + paeq5(j1,j2,j3)*zaea(kspec,5)

          zaeods(j1,j2) = paeq1(j1,j2,j3)*zaes(kspec,1)*(1.-zaef(kspec,1)) +             &
                  tau_streu_dust(j1,j3,kspec)*(1.-(asym_ges(j1,j3,kspec)**2))    &
                 +tau_streu_seas(j1,j3,kspec)*(1.-(asym_seas(j1,j3,kspec)**2))   &
                 +paeq3(j1,j2,j3)*zaes(kspec,3)*(1.-zaef(kspec,3)) +             &
                  tau_streu_aero(j1,j3,kspec)*(1.-(asym_aero(j1,j3,kspec)**2))   &
                 +paeq4(j1,j2,j3)*zaes(kspec,4)*(1.-zaef(kspec,4))               &
                 +paeq5(j1,j2,j3)*zaes(kspec,5)*(1.-zaef(kspec,5))

          zzg= ( paeq1(j1,j2,j3)*zaes(kspec,1)*(1.-zaef(kspec,1))*zaeg(kspec,1) +                     &
                  tau_streu_dust(j1,j3,kspec)*(1.-(asym_ges(j1,j3,kspec)**2))*asym_ges(j1,j3,kspec)    &
                 +tau_streu_seas(j1,j3,kspec)*(1.-(asym_seas(j1,j3,kspec)**2))*asym_seas(j1,j3,kspec)  &
                 +paeq3(j1,j2,j3)*zaes(kspec,3)*(1.-zaef(kspec,3))*zaeg(kspec,3)  +                    &
                  tau_streu_aero(j1,j3,kspec)*(1.-(asym_aero(j1,j3,kspec)**2.))*asym_aero(j1,j3,kspec) &
                 +paeq4(j1,j2,j3)*zaes(kspec,4)*(1.-zaef(kspec,4))*zaeg(kspec,4)                       &
                 +paeq5(j1,j2,j3)*zaes(kspec,5)*(1.-zaef(kspec,5))*zaeg(kspec,5) )                     &
                / MAX( zaeods(j1,j2),zepopd)
          zaeb0 (j1,j2) = z1d8*(4.+zzg)/(1.+zzg)
        ENDDO
      ENDDO

    ELSEIF ((.NOT. lrad_dust) .AND. (lrad_seas) .AND. (lrad_aero)) THEN

      DO j2 = ki2sc, ki2ec
        DO j1 = ki1sc, ki1ec
          zaeoda(j1,j2) = paeq1(j1,j2,j3)*zaea(kspec,1) &
                + tau_abs_seas(j1,j3,kspec) &
                + paeq3(j1,j2,j3)*zaea(kspec,3) + tau_abs_aero(j1,j3,kspec) &
                + paeq4(j1,j2,j3)*zaea(kspec,4) &
                + paeq5(j1,j2,j3)*zaea(kspec,5)

          zaeods(j1,j2) = paeq1(j1,j2,j3)*zaes(kspec,1)*(1.-zaef(kspec,1))               &
                 +tau_streu_seas(j1,j3,kspec)*(1.-(asym_seas(j1,j3,kspec)**2.))  &
                 +paeq3(j1,j2,j3)*zaes(kspec,3)*(1.-zaef(kspec,3)) +             &
                  tau_streu_aero(j1,j3,kspec)*(1.-(asym_aero(j1,j3,kspec)**2.))  &
                 +paeq4(j1,j2,j3)*zaes(kspec,4)*(1.-zaef(kspec,4))               &
                 +paeq5(j1,j2,j3)*zaes(kspec,5)*(1.-zaef(kspec,5))

          zzg=( paeq1(j1,j2,j3)*zaes(kspec,1)*(1.-zaef(kspec,1))*zaeg(kspec,1)                       &
               +tau_streu_seas(j1,j3,kspec)*(1.-(asym_seas(j1,j3,kspec)**2.))*asym_seas(j1,j3,kspec) &
               +paeq3(j1,j2,j3)*zaes(kspec,3)*(1.-zaef(kspec,3))*zaeg(kspec,3) +                     &
                tau_streu_aero(j1,j3,kspec)*(1.-(asym_aero(j1,j3,kspec)**2.))*asym_aero(j1,j3,kspec) &
               +paeq4(j1,j2,j3)*zaes(kspec,4)*(1.-zaef(kspec,4))*zaeg(kspec,4)                       &
               +paeq5(j1,j2,j3)*zaes(kspec,5)*(1.-zaef(kspec,5))*zaeg(kspec,5) )                     &
                   / MAX( zaeods(j1,j2),zepopd)
          zaeb0 (j1,j2) = z1d8*(4.+zzg)/(1.+zzg)
        ENDDO
      ENDDO
   
    ENDIF
      
  ENDIF
#endif

    IF (ldebug) THEN
      print *,' zaeoda  (j1b,j2b) :',zaeoda  (j1b,j2b)
      print *,' zaeods  (j1b,j2b) :',zaeods  (j1b,j2b)
      print *,' zaeb0   (j1b,j2b) :',zaeb0   (j1b,j2b)
      print *,' zaeb    (j1b,j2b) :',zaeb    (j1b,j2b)
    ENDIF   

    ! Optical thickness for Rayleigh scattering

    DO j2 = ki2sc, ki2ec
      DO j1 = ki1sc, ki1ec
        zraods(j1,j2) = zrsc(isp) * pdp(j1,j2,j3) / (1.0_ireals + zrsc(isp) * &
               (papra(j1,j2) + 0.5_ireals * pdp(j1,j2,j3) * pqsmu0(j1,j2)))
        papra (j1,j2) = papra(j1,j2) + pdp(j1,j2,j3) * pqsmu0(j1,j2)
      ENDDO
    ENDDO

    IF (ldebug) THEN
      !  print *,' Rayleigh coefficient:',zrsc(isp)
      !  print *,' Papra  (j1b,j2b)    :',papra   (j1b,j2b)
      !  print *,' Pqsmu0 (j1b,j2b)    :',pqsmu0  (j1b,j2b)
      !  print *,' Pdp    (j1b,j2b)    :',pdp     (j1b,j2b,j3)
         print *,' zraods (j1b,j2b)    :',zraods  (j1b,j2b)
    ENDIF   

    !-----------------------------------------------------------------------
 
    ! Linear combination of individual contributions 
 

    ! a) cloud free part of layer

    DO j2 = ki2sc, ki2ec
      DO j1 = ki1sc, ki1ec
        podaf(j1,j2,j3) = MAX ( zaeoda(j1,j2)                , zepopd)
        podsf(j1,j2,j3) = MAX ( zaeods(j1,j2) + zraods(j1,j2), zepopd)
        podsf(j1,j2,j3) = MIN ( podsf (j1,j2,j3),                             &
                   (1.0_ireals-zepssa) * (podaf(j1,j2,j3) + podsf(j1,j2,j3)) )
        pbsff(j1,j2,j3) = (zaeb0(j1,j2) * zaeods(j1,j2)                       &
                           + 0.5_ireals * zraods(j1,j2)) / podsf(j1,j2,j3)
        pusff(j1,j2,j3) = (zaeb (j1,j2) * zaeods(j1,j2)                       &
                           + 0.5_ireals * zraods(j1,j2)) / podsf(j1,j2,j3)
      ENDDO
    ENDDO
 
!-------------------------------------------------------------------------------
    ! b) cloudy part of layer

    DO j2 = ki2sc, ki2ec
      DO j1 = ki1sc, ki1ec
        podac(j1,j2,j3) = zlwoda(j1,j2) + ziwoda(j1,j2) + zaeoda(j1,j2)
        podsc(j1,j2,j3) = zlwods(j1,j2) + ziwods(j1,j2)                        &
                                        + zaeods(j1,j2) + zraods(j1,j2)
        podac(j1,j2,j3) = MAX( podac(j1,j2,j3), zepopd)
        podsc(j1,j2,j3) = MAX( podsc(j1,j2,j3), zepopd)
        podsc(j1,j2,j3) = MIN( podsc(j1,j2,j3),                                &
                  (1.0_ireals-zepssa) * (podac(j1,j2,j3) + podsc(j1,j2,j3)))

        pbsfc(j1,j2,j3)= (zlwb0(j1,j2) * zlwods(j1,j2)                         &
                 + ziwb0(j1,j2) * ziwods(j1,j2) + zaeb0(j1,j2) * zaeods(j1,j2) &
                 + 0.5 * zraods(j1,j2)) / podsc(j1,j2,j3)
        pusfc(j1,j2,j3)= (zlwb (j1,j2) * zlwods(j1,j2)                         &
                 + ziwb (j1,j2) * ziwods(j1,j2) + zaeb (j1,j2) * zaeods(j1,j2) &
                 + 0.5 * zraods(j1,j2)) / podsc(j1,j2,j3)
      ENDDO
    ENDDO
 
  ! End of vertical loop
  END DO  
 
  IF (ldebug) THEN
      print *,' podaf   (j1b,j2b,1) :',podaf   (j1b,j2b,1)
      print *,' podsf   (j1b,j2b,1) :',podsf   (j1b,j2b,1)
      print *,' pbsff   (j1b,j2b,1) :',pbsff   (j1b,j2b,1)
      print *,' pusff   (j1b,j2b,1) :',pusff   (j1b,j2b,1)
      print *,' podac   (j1b,j2b,1) :',podac   (j1b,j2b,1)
      print *,' podsc   (j1b,j2b,1) :',podsc   (j1b,j2b,1)
      print *,' pbsfc   (j1b,j2b,1) :',pbsfc   (j1b,j2b,1)
      print *,' pusfc   (j1b,j2b,1) :',pusfc   (j1b,j2b,1)
      print *,' -------------------------------------------------'
  ENDIF     
 
!------------------------------------------------------------------------------
! End of the subroutine
!------------------------------------------------------------------------------

END SUBROUTINE opt_so

!==============================================================================
!==============================================================================
!+ Module procedure in "Radiation" 
!------------------------------------------------------------------------------

SUBROUTINE aerdis ( petah,  pvdaes, pvdael, pvdaeu, pvdaed, klevp1, &
                    ptrbga, pvobga, pstbga, paeops, paeopl, paeopu, &
                    paeopd, ptrpt , paeadk, paeadm )

!------------------------------------------------------------------------------
!
! Description:
!
! The module procedure aerdis provides parameters for the vertical distribution
! of aerosols (based on the original code of J.F. Geleyn (ECMWF, 4.11.82).
!
! The routine computes the values PVDAE* (* = s, l, u or d for sea, land
! urban or desert) of a surfach-normalised vertical distribution of aerosols'
! optical depth from the argument petah (vertical coordinate) at klevp1 levels.
! It also sets values for non-geograpically weighted total optical depths (at
! 55 micrometer wavelength) paeopn for the same four types and similar optical
! depths diveded by pressure for bachground well-mixed aerosols of three types
! p**bga (** = tr, vo or st for tropospheric, volcanic (stratosperic ashes) or
! stratosperic (sulfuric type)). It finally sets values for the power to be
! applied to a temperature ratio smaller than two in order to obtain an index
! one in the stratosphere and zero in the troposphere with a relatively smooth
! transistion (ptrpt), as well as for adsorption coefficients fo water to the
! three type of troposperic aerosols (paeadk) with a minimum value ( in the 
! whole atmosphere) for the sum of the products paeadk by the optical depths
! divided by pressure thickness: paeadm. 
!
! Method:
!
! Straightforward, equivalent heights are given in meters (8434 for the
! atmosphere) and tropospheric and stratospheric pressure boundary values
! are set at 101325 and 19330 Pascal. 
!
!------------------------------------------------------------------------------

! Subroutine arguments:
! --------------------

! Input data
! ----------
  INTEGER (KIND=iintegers), INTENT (IN) ::  &
     klevp1           ! number of model layer interfaces      

  REAL    (KIND=ireals   ), INTENT (IN) ::  &
     petah(klevp1)    ! normalized vertical coordinate at half levels

! Output data
! -----------
  REAL    (KIND=ireals   ), INTENT (OUT) ::  &
    pvdaes(klevp1), & ! normalized vertical distribution (sea)
    pvdael(klevp1), & ! normalized vertical distribution (land)
    pvdaeu(klevp1), & ! normalized vertical distribution (urban)
    pvdaed(klevp1), & ! normalized vertical distrubution (desert)
    ptrbga        , & ! b. optical depths div. by pressure (tropospheric)
    pvobga        , & ! b. optical depths div. by pressure (volcanic)
    pstbga        , & ! b. optical depths div. by pressure (stratospheric)
    paeops        , & ! total opt. depths for ver. varying aerosols (sea)
    paeopl        , & ! total opt. depths for ver. varying aerosols (land)
    paeopu        , & ! total opt. depths for ver. varying aerosols (urban)
    paeopd        , & ! total opt. depths for ver. varying aerosols (desert)
    ptrpt         , & ! temperature exponent for the stratosperic definition
    paeadk(3)     , & ! constants for definition of the quantity of water 
    paeadm            ! vapour that will be adsorbed to the dry aerosols to
                      ! form moist aerosols

! Local parameters:
! -------------
  REAL (KIND=ireals), PARAMETER  ::  &
    zhss = 8434.0_ireals/1000.0_ireals ,  & !
    zhsl = 8434.0_ireals/1000.0_ireals ,  & !
    zhsu = 8434.0_ireals/1000.0_ireals ,  & !
    zhsd = 8434.0_ireals/3000.0_ireals      !

 
!- End of header
!==============================================================================

!------------------------------------------------------------------------------
! Begin Subroutine aerdis              
!------------------------------------------------------------------------------

  pvdaes(1) = 0.0_ireals
  pvdael(1) = 0.0_ireals
  pvdaeu(1) = 0.0_ireals
  pvdaed(1) = 0.0_ireals

  IF(petah(1).NE.0.) THEN
     pvdaes(1) = petah(1)**zhss
     pvdael(1) = petah(1)**zhsl
     pvdaeu(1) = petah(1)**zhsu
     pvdaed(1) = petah(1)**zhsd
  END IF

  pvdaes(2:klevp1) = petah(2:klevp1)**zhss
  pvdael(2:klevp1) = petah(2:klevp1)**zhsl
  pvdaeu(2:klevp1) = petah(2:klevp1)**zhsu
  pvdaed(2:klevp1) = petah(2:klevp1)**zhsd

  ptrbga = 0.03_ireals  / (101325.0_ireals - 19330.0_ireals)
  pvobga = 0.007_ireals /  19330.0_ireals
  pstbga = 0.045_ireals /  19330.0_ireals

  paeops = 0.05_ireals
  paeopl = 0.2_ireals
  paeopu = 0.1_ireals
  paeopd = 1.9_ireals
  ptrpt  = 30.0_ireals

  paeadk(1) = 0.3876E-03_ireals
  paeadk(2) = 0.6693E-02_ireals
  paeadk(3) = 0.8563E-03_ireals
  paeadm    = 2.6000E-10_ireals

 
!------------------------------------------------------------------------------
! End of the subroutine 
!------------------------------------------------------------------------------

END SUBROUTINE aerdis  

!==============================================================================
!==============================================================================

SUBROUTINE calc_rad_corrections (thetain, phi, horizonte, smu0, rlati,     &
                                 rloni, deksini, dekcosi, zeit0i, fcor,    &
                                 idim, jdim, nhordim, isc, iec, jsc, jec, idbg)

!------------------------------------------------------------------------------
!
! Description:
!   Compute parameters swdir_cor needed for the grid scale topographic 
!   correction of direct solar radiation at the surface.
!   If correction option is chosen, this subroutine is called before fesft.
!   Following Mueller and Sherrer (2005), MWR
!
!------------------------------------------------------------------------------

! Variables of the parameter list
INTEGER (KIND=iintegers),  INTENT (IN)  ::   &
  idim, jdim, nhordim,           & ! field dimensions
  isc, iec, jsc, jec,            & ! start and end indices for computation
  idbg                             ! for debug output

REAL    (KIND=ireals   ),  INTENT (IN)     ::   &
  thetain  (idim,jdim),          & ! slope angle
  phi      (idim,jdim),          & ! slope aspect
  horizonte(idim,jdim,nhordim),  & ! horizont
  rlati    (idim,jdim),          & ! latitude (geogr) [rad]
  rloni    (idim,jdim),          & ! longitude(geogr) [rad]
  smu0     (idim,jdim)             ! sun zenith angle

REAL    (KIND=ireals   ),  INTENT (IN)     ::   &
  deksini,                       & ! sin of sun declination angle
  dekcosi,                       & ! cos of sun declination angle
  zeit0i                           !T.R.

! Output of the routine
REAL    (KIND=ireals     ),  INTENT (OUT) ::   &
  fcor     (idim,jdim)

! Local parameters and variables
REAL (KIND=ireals), PARAMETER :: &
  zepemu = 1.0E-07,              &
  rtod = 57.2957795                ! radiantas to degrees

REAL    (KIND=ireals)   ::       &
  zeitrad,                       & ! T.R.
  phi_s,                         & !
  phi_sun  (idim,jdim),          & ! sun azimuth angle [rad]
  theta_sun(idim,jdim),          & ! sun elevation angle [rad]
  theta    (idim,jdim),          & ! theta for computation
  x1,x2,ha_sun                     !


INTEGER (KIND=iintegers) ::        &
  ii,shadow,i,k,j

LOGICAL                  ::        &
  lshade, lslope_aspect            !switches

!------------------------------------------------------------------------------

lshade        = .TRUE.
lslope_aspect = .TRUE.

DO j = jsc, jec
  DO i = isc, iec

    ! sun elevation angle
    theta_sun(i,j) = ASIN(smu0(i,j))

    ! sun azimuth angle
    zeitrad = zeit0i + rloni(i,j)  !T.R.
    x1 = dekcosi * SIN(zeitrad) / COS(theta_sun(i,j))
    x2 = ( SIN(rlati(i,j)) * dekcosi * COS(zeitrad) - &
                        COS(rlati(i,j)) * deksini ) / COS(theta_sun(i,j))
    IF (x2 < -1) x2 = -1
    IF (x2 >  1) x2 = 1
    phi_s = COS(x2)
    IF (x1 < 0) phi_s = - phi_s
    phi_sun(i,j) = phi_s + pi

    ! sun elevation angle corrected by refraction (enpiric formula, A.H.Siemer(1988))
    theta_sun(i,j) = theta_sun(i,j) + (1.569000_ireals - theta_sun(i,j)) / &
                                    (185.5_ireals + 3620.0_ireals * theta_sun(i,j))

    ! night or day
    IF (theta_sun(i,j) < 0.0) THEN
       theta_sun(i,j) = 0.0_ireals
    ENDIF

    ! compute shadow
    ! the horizon has a spatial resolution of 360/nhordim degrees. a distance weighted
    ! linear interpolation is done between the two neighboring points.
    ii = INT(rtod*phi_sun(i,j)/(360.0_ireals/nhordim))
    IF (ii >= nhordim) THEN
      ii = nhordim - 1
    ENDIF

    k  = MOD(ii+1,24)

    IF (horizonte(i,j,k+1) > 90.0_ireals .OR. horizonte(i,j,k+1) < 0.0_ireals) THEN
      PRINT *,'!!ERROR!!, horizon_angle > 90deg or < 0deg',horizonte(i,j,k+1)
      !US there shall be no stop in a parallel program!!!  STOP
    ENDIF

    ha_sun = (horizonte(i,j,k+1)*(rtod*phi_sun(i,j)-15*ii)+ &
              horizonte(i,j,ii+1)*(15*(ii+1)-rtod*phi_sun(i,j)))/15.0

    ! compute shadowmask
    shadow = 1
    IF (rtod*theta_sun (i,j) < ha_sun .AND. lshade) THEN
      shadow = 0
    ENDIF

    ! compute fcor
    ! slope angle and aspect switched off
    IF (.NOT. lslope_aspect) THEN
      theta(i,j) = 0.0_ireals
    ELSE
      theta(i,j) = thetain(i,j)
    ENDIF
 
    IF (theta_sun(i,j) > 0.01_ireals) THEN
      ! Mueller and Scherrer (2005) formula (MWR)
      ! fcor(i,j) = shadow * (1 + ( tan(theta(i,j)) / tan(theta_sun(i,j)) )*&
      !              cos( phi_sun(i,j) - phi(i,j)) )

      ! New formula (lower correction, theoretically correct derived)
      fcor(i,j) = shadow * ( COS(theta(i,j)) + (SIN(theta(i,j))/TAN(theta_sun(i,j)) )*&
                    COS( phi_sun(i,j) - phi(i,j)) )
    ELSE
      fcor(i,j) = 1.0
    ENDIF

    ! Consistency check to avoid negative corrections:
    ! active in situations with low sun elevation (slope angles > sun elevation)
    ! when the slope aspect is greater than the daily maxima or smaller than the
    ! daily minima of the sun azimuth angle (during the sunshine time, a kind 
    ! of self shading effect).
    IF (fcor(i,j) < 0.0_ireals) THEN
      fcor(i,j) = 0.0
    ENDIF

    IF ( (idbg > 20) .AND. (i == 4) ) THEN
      PRINT *, '   calc_rad_corrections:  debug point:  ', i, j
      PRINT *, '   deksini, dekcosi, zeitrad = ', deksini, dekcosi, zeitrad
      PRINT *, '   rlat, rlon, theta, phi, smu0 = ', rlati(i,j), rloni(i,j), &
                                                      theta(i,j), phi(i,j), smu0(i,j)
      PRINT *, '   ii, k, horizon = ', ii, k, horizonte(i,j,ii+1), horizonte(i,j,k+1)
      PRINT *, '   ha_sun, theta_sun, phi_sun = ', ha_sun, rtod*theta_sun(i,j), &
                                                           rtod*phi_sun(i,j)

      IF ( (shadow == 0) .AND. (theta_sun(i,j) > 0.01_ireals) ) THEN
        PRINT *, '    calc_rad_corrections:  ', 'DAY-SHADOW'
      ENDIF
      IF ( (shadow == 1) .AND. (theta_sun(i,j) > 0.01_ireals) ) THEN
        PRINT *, '    calc_rad_corrections:  ', 'DAY-SUN'
      ENDIF
      IF (theta_sun(i,j) <= 0.01_ireals ) then
        PRINT *, '    calc_rad_corrections:  ', 'NIGHT'
      ENDIF
      PRINT *, '    calc_rad_corrections: fcor  ', fcor(i,j)
    ENDIF

  ENDDO
ENDDO


END SUBROUTINE calc_rad_corrections

!==============================================================================

END MODULE src_radiation
