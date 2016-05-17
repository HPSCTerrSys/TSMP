!+ Data module for all global meteorological fields
!------------------------------------------------------------------------------

MODULE data_fields

!------------------------------------------------------------------------------
!
! Description:
!  This module declares all meteorological fields that have to reside in 
!  the long term storage, i.e. that are used in more than one module.
!  Fields included are
!    - constant fields defining the reference atmosphere
!    - external parameter fields
!    - prognostic variables
!    - tendency fields for the prognostic variables
!    - fields for surface values
!    - fields that are computed in the parametrization packages 
!      or in the dynamics
!    - fields for model-output and diagnostics
!    - fields for the boundary values
!
!  All fields are declared as allocatable arrays. They are allocated in the
!  setup of the model and deallocated in the cleanup at the end of the
!  program.
!
! Current Code Owner: DWD, Ulrich Schaettler
!  phone:  +49  69  8062 2739
!  fax:    +49  69  8062 3721
!  email:  ulrich.schaettler@dwd.de
!
! History:
! Version    Date       Name
! ---------- ---------- ----
! 1.1        1998/03/11 Ulrich Schaettler
!  Initial release
! 1.3        1998/04/15 Guenther Doms
!  Definition of new tendency arrays for convection
! 1.7        1998/07/16 Guenther Doms
!  Removal of global array 'rrssk'.
! 1.20       1999/01/07 Guenhter Doms
!  Renaming of some global variables.
! 1.24       1999/03/01 Guenther Doms
!  Declaration of a new 3-D prognostic variable (qi).
! 1.30       1999/06/24 Matthias Raschendorfer
!  Declaration of a 3-D prognostic var. (tke) and its tendenci field (tketens).
!  Declaration of a 3-D array (rcld).
!  Declaration of 5 3-D arrays for canopy layers (c_big,c_sml, r_air, t_e,qv_e).
!  Declaration of a 2-D arrays (tfh,tfm) and (h_can,d_pat).
!  Declaration of 4 2-D variables for new convection closures
! 1.33       1999/10/14 Reinhold Hess
!  Declaration of new global 2-D arrays idiv_hum and aevap_s for diagnosis
!  of the model water budget
!  Declaration of a new 2-D array 'sai' for surface area index (M.Raschendorfer)
! 1.34       1999/12/10 Ulrich Schaettler
!  Named all boundary fields with "_bd" (for consistency)
! 1.39       2000/05/03 Ulrich Schaettler
!  Add declaration of a boundary field for w (used for interactive nesting)
! 2.2        2000/08/18 Matthias Raschendorfer
!  Declaration of the 2-D arrays 'eai' and 'tai'.
! 2.4        2001/01/29 Christoph Schraff
!  Declaration of the 2-D array 'prne_con', for humidity balancing at T-nudging.
! 2.8        2001/07/06 Ulrich Schaettler
!  Added new fields for multi-layer soil model and surface fluxes
! 2.11       2001/09/28 Ulrich Schaettler
!  Added new fields for lateral values of cloud ice
! 2.17       2002/05/08 Ulrich Schaettler
!  New fields for Kain-Fritsch convection scheme
! 2.18       2002/07/16 Ulrich Schaettler
!  New fields for specific rain and snow content;
!  included declaration of a1t, a2t from src_leapfrog
! 3.5        2003/09/02 Ulrich Schaettler
!  New fields phi_tot, rla_tot to avoid global communication in the radiation
! 3.6        2003/12/11 Reinhold Schrodin
!  New field freshsnow for new multi-layer soil model added
! 3.7        2004/02/18 Ulrich Schaettler
!  New fields for computing synthetic satellite images (synme5-7, synmsg)
!  New field for storing convective cloud water (clw_con)
!  Renamed alb (alb_rad), idiv_hum (tdiv_hum), phi (rlat), rla (rlon),
!      cphi (crlat), acphir (acrlat), tgphi (tgrlat) (for consistency with GME
! 3.13       2004/12/03 Ulrich Schaettler
!  New fields for graupel scheme (qg, prg_gsp, grau_gsp): Thorsten Reinhardt
!  New fields for 3D turbulence  (tkhm, tkhh): Jochen Foerstner
!  New fields for fast_waves_rk: 1/SQRT(G): sqrtg_r_(s,u,v,w): Jochen Foerstner
!  New fields for external parameters (not yet used: for_e, for_d) R. Schrodin
!  Renamed w_ice to w_so_ice (to be consistent with GME)
! 3.17       2005/12/12 Reinhold Schrodin
!  New fields rho_snow (prognostic snow density) and h_snow (snow height) added
! 3.18       2006/03/03 Ulrich Schaettler
!  New fields for the CLM version              (by CLM Community)
!  New fields for the lake model FLake         (Dmitrii Mironov)
!  New field rh_2m for relative humidity in 2m (Matthias Raschendorfer)
!  New field tinc_lh to gather temperature increments due to latent heat 
!   conversion in microphysics (mainly saturation adjustment; Jochen Foerstner)
! 3.19       2006/04/25 Ulrich Schaettler
!  New field T_S_LAKE to save lake values from Nudging cycle
! 3.21       2006/12/04 Jochen Foerstner, Burkhardt Rockel, Christoph Schraff
!  New variables (LMK): t0, dt0dz, hd_mask_dcoeff, qvt_diff
!  Put declaration of zwcon here, because it is used in several modules
!  New boundary variables: qr_bd, qs_bd, qg_bd
!  Renamed variable sunshhrs, sodwdir to dursun, sodwddm
!  Additional variable introduced for deep atmosphere (Ronny Petrik)
!  Additional fields for Bechtold convection scheme (MeteoSwiss)
!  New fields for time-integrated analysis increments from nudging (section 10)
! V3_23        2007/03/30 Matthias Raschendorfer, Matteo Buzzi, Jochen Foerstner
!  Added field 'tfv' containing the laminar reduction factor for evaporation
!  Added fields for computations of topographic radiation correction (Matteo Buzzi)
!  Added fields for new relaxation (Jochen Foerstner)
! V4_3         2008/02/25 Matthias Raschendorfer
!  Introduction of a 3D diagnostic field 'edr' for the eddy dissipotion rate
! V4_4         2008/07/16 Jan-Peter Schulz
!  Added external parameters (sso_stdh, sso_gamma, sso_theta, sso_sigma) and
!  fields (ut_sso, vt_sso, tt_sso, ustr_sso, vstr_sso, vdis_sso,
!  austr_sso, avstr_sso, avdis_sso) for SSO scheme
! V4_5         2008/09/10 Ulrich Schaettler
!  Moved field for damping coefficients from src_relaxation to data_fields 
!  for global use
! V4_8         2009/02/16 Ulrich Schaettler
!  New fields for convective and dynamical gust in 10m
!  New fields for additional output of radiation values (for CLM)
!  Add global field for reference pressure at half levels (p0hl) (G. Zaengl)
! V4_9         2009/07/16 Ulrich Schaettler
!  New fields for diffusion masks for t-, q- and u-fields
! V4_10        2009/09/11 Ulrich Schaettler
!  New field for snow melt
!  New field for sea-mask (Jan-Peter Schulz for sea-ice model)
! V4_11        2009/11/30 Ekaterina Machulskaya, Juergen Helmert, Lucio Torrisi
!  Additional fields for the snow model (EM)
!  Additional external parameters for aerosol distributions, surface emissivity
!   and stomata resistance;  (JH)
!  Additional fields in output for SMA (LT)
! V4_12        2010/05/11 Michael Baldauf, Ulrich Schaettler
!  New fields dzeta_dlam, dzeta_dphi
!  Renamed hd_mask_dcoeff to hd_mask_dcoeff_p for diffusion of pressure
! V4_13        2010/05/11 Michael Gertz
!  Adaptions to SVN
! V4_15        2010/11/19 Ulrich Schaettler
!  The field t_s_lake is removed again after adaptations in the SST Analysis
! V4_20        2011/08/31 Ulrich Schaettler
!  tgrlat needs 2 dimensions because of v-point dependence 
!    (reported by Andreas Will)
!  Introduction of additional 3D-arrays 'tket_(conv, sso, hshr)' for additional
!   TKE source by the action of other sub grid scale flow patterns.
!    (by Matthias Raschendorfer)
!  V4.21        2013/02/12 Prabhakar Shrestha
!  Inclusion of transfer coefficient for moisture to be used with CLM
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

!==============================================================================

IMPLICIT NONE

!==============================================================================

! Global (i.e. public) Declarations:

! 1. constant fields for the reference atmosphere                     (unit)
! -----------------------------------------------

  REAL  (KIND=ireals), TARGET, ALLOCATABLE ::           &
    rho0 (:,:,:),     & ! reference density at the full model levels  (kg/m3)
    dp0  (:,:,:),     & ! reference pressure thickness of layers      ( Pa  )
    p0   (:,:,:),     & ! reference pressure at full levels           ( Pa  )
    p0hl (:,:,:),     & ! reference pressure at half levels           ( Pa  )
    dt0dz(:,:,:),     & ! temperature grad. of reference atmosphere   ( K/m )
    t0   (:,:,:),     & ! reference temperature                       ( K   )
    hhl  (:,:,:),     & ! geometrical height of half levels           ( m   )
    sqrtg_r_s(:,:,:), & ! 1 / square root of G at scalar points       ( 1/m )
    sqrtg_r_u(:,:,:), & ! 1 / square root of G at u points            ( 1/m )
    sqrtg_r_v(:,:,:), & ! 1 / square root of G at v points            ( 1/m )
    sqrtg_r_w(:,:,:), & ! 1 / square root of G at w points            ( 1/m )
    dzeta_dlam(:,:,:),& ! d zeta / d lambda (for constant phi,    z)
                        ! at the scalar position                      ( 1   )
    dzeta_dphi(:,:,:)   ! d zeta / d phi    (for constant lambda, z)
                        ! at the scalar position                      ( 1   )

! 2. external parameter fields                                        (unit)
! ----------------------------

  REAL  (KIND=ireals), TARGET, ALLOCATABLE ::           &
    hsurf  (:,:),   & ! height of surface topography                  ( m   )
    sso_stdh (:,:), & ! standard deviation of sub-grid scale orography( m   )
    sso_gamma(:,:), & ! anisotropy of sub-grid scale orography          --
    sso_theta(:,:), & ! angle betw. principal axis of orography and E ( rad )
    sso_sigma(:,:), & ! mean slope of sub-grid scale orography          --
    aer_su (:,:),   & ! monthly aerosol climatology sulfate drops     (0 - 1)
    aer_du (:,:),   & ! monthly aerosol climatology total dust        (0 - 1)
    aer_or (:,:),   & ! monthly aerosol climatology organic (water sol.)(0 - 1)
    aer_bc (:,:),   & ! monthly aerosol climatology black carbon      (0 - 1)
    aer_ss (:,:),   & ! monthly aerosol climatology sea salt          (0 - 1)
    emis_rad (:,:), & ! external thermal emissivity map               (0 - 1)
    rsmin2d(:,:)  , & ! minimum stomata resistance                    ( s/m )
    swi    (:,:,:), & ! soil wetness index                            (0 - 1)
    gz0    (:,:),   & ! surface roughness * g                         (m2/s2)
    fr_land(:,:),   & ! fraction of land in a grid element              --
    soiltyp(:,:),   & ! type of the soil (keys 0-9)                     --
    vio3   (:,:),   & ! vertical integrated ozone contents            (Pa O3)
    hmo3   (:,:),   & ! ozone maximum                                 ( Pa  )
    rlat   (:,:),   & ! geographical latitude                         ( rad )
    rlon   (:,:),   & ! geographical longitude                        ( rad )
    rlattot(:,:),   & ! geographical latitude                         ( rad )
    rlontot(:,:),   & ! geographical longitude                        ( rad )
    fccos  (:,:),   & ! horizontal coriolis-parameter                 ( 1/s )
    fc     (:,:),   & ! coriolis-parameter                            ( 1/s )
    rmy    (:,:,:), & ! Davis-parameter for relaxation (mass, qv, qc)   --
    rmyq   (:,:),   & ! Davis-parameter for relaxation (qr, qs, qg)     --
    hd_mask_dcoeff_p(:,:,:), & ! 3D-domain mask for horizontal diffusion  --
    hd_mask_dcoeff_t(:,:,:), & ! 3D-domain mask for horizontal diffusion  --
    hd_mask_dcoeff_q(:,:,:), & ! 3D-domain mask for horizontal diffusion  --
    hd_mask_dcoeff_u(:,:,:), & ! 3D-domain mask for horizontal diffusion  --
    ofa_hdx(:,:,:), & !
    ofa_hdy(:,:,:), & !
    hd_mask(:,:,:)    ! 3D-domain mask for horizontal diffusion         --

  REAL  (KIND=ireals), TARGET, ALLOCATABLE ::           &
    crlat  (:,:),   & ! cosine of transformed latitude                  --
    acrlat (:,:),   & ! 1 / ( crlat * radius of the earth )           ( 1/m )
    tgrlat (:,:),   & ! tangens of transformed latitude                 --
    aerlan (:,:),   & ! aerosol-distribution on rural areas             --
    aerurb (:,:),   & ! aerosol-distribution on urban areas             --
    aerdes (:,:),   & ! aerosol-distribution on desert areas            --
    aersea (:,:),   & ! aerosol-distribution on the sea                 --
    plcov  (:,:),   & ! fraction of plant cover                         --
    lai    (:,:),   & ! leaf area index of plants                       --
    tai    (:,:),   & ! transpiration area index                        --
    sai    (:,:),   & ! surface area index                              --
    eai    (:,:),   & ! (evaporative) earth area index                  --
    rootdp (:,:),   & ! depth of the roots                            ( m  )
    for_e  (:,:),   & ! ground fraction covered by evergreen forest     --
    for_d  (:,:),   & ! ground fraction covered by deciduous forest     --

    h_can (:,:),    & ! hight of the vertically resolved canopy       ( m )
    d_pat (:,:),    & ! horizontal pattern length scale               ( m )
    c_big (:,:,:),  & ! effective drag coefficient of canopy elements
                      ! larger than or equal to the tubulent length
                      ! scale                                         (1/m)
    c_sml (:,:,:),  & ! effective drag coefficient of canopy elements
                      ! smaller than the tubulent length scale        (1/m)
    r_air (:,:,:)     ! air containing fraction of a gridbox inside
                      ! the canopy                                    ( 1 )

  LOGICAL, ALLOCATABLE ::           &
    least_lbdz(:,:), & ! mask for eastern  lateral boundary zone
    lwest_lbdz(:,:), & ! mask for western  lateral boundary zone
    lnorth_lbdz(:,:),& ! mask for northern lateral boundary zone
    lsouth_lbdz(:,:)   ! mask for southern lateral boundary zone

  ! external parameter fields for the lake model FLake
  REAL  (KIND=ireals), TARGET, ALLOCATABLE ::           &
    fr_lake(:,:),   & ! lake fraction in a grid element [0,1]         (  -  )
    depth_lk(:,:),  & ! lake depth                                    (  m  )
    fetch_lk(:,:),  & ! wind fetch over lake                          (  m  )
    dp_bs_lk(:,:),  & ! depth of the thermally active layer
                      ! of bottom sediments                           (  m  )
    t_bs_lk (:,:),  & ! climatological temperature at the bottom of
                      ! the thermally active layer of sediments       (  K  )
    gamso_lk(:,:)     ! attenuation coefficient for
                      ! solar radiation in lake water                 ( 1/m )


  LOGICAL, ALLOCATABLE ::           &
    llandmask(:,:), & ! landpoint mask
    lseamask(:,:)     ! ocean point mask, i.e. water but not lake

! 3. prognostic variables                                             (unit)
! -----------------------

  REAL  (KIND=ireals), TARGET, ALLOCATABLE ::           &
    u (:,:,:,:),    & ! zonal wind speed                              ( m/s )
    v (:,:,:,:),    & ! meridional wind speed                         ( m/s )
    w (:,:,:,:),    & ! vertical wind speed (defined on half levels)  ( m/s )
    t (:,:,:,:),    & ! temperature                                   (  k  )
    qv(:,:,:,:),    & ! specific water vapor content                  (kg/kg)
    qc(:,:,:,:),    & ! specific cloud water content                  (kg/kg)
    qi(:,:,:,:),    & ! specific cloud ice content                    (kg/kg)
    qr(:,:,:,:),    & ! specific rain content                         (kg/kg)
    qs(:,:,:,:),    & ! specific snow content                         (kg/kg)
    qg(:,:,:,:),    & ! specific graupel content                      (kg/kg)
    pp(:,:,:,:),    & ! deviation from the reference pressure         ( pa  )

! fields of the turbulent scheme defined on half-levels:

   tke(:,:,:,:),    & ! SQRT(2 * turbulent kinetik energy)            ( m/s )
   edr(:,:,:),      & ! eddy dissipation rate of TKE (EDR)            (m2/s3)
   tketens(:,:,:),  & ! tke-tendency (defined on half-levels)         ( m/s )
   tket_conv(:,:,:),& ! TKE-tendency due to convective buoyancy       ( m2/s3)
   tket_hshr(:,:,:),& ! TKE-tendency due to (sep.) horiz. shear       ( m2/s3)
   tket_sso (:,:,:)   ! TKE-tendency due to SSO wake production       ( m2/s3)

! 4. tendency fields for the prognostic variables                     (unit )
! -----------------------------------------------
!    time tendencies  by diabatic and adiabatic processes
!    without sound-wave terms

  REAL  (KIND=ireals), ALLOCATABLE ::           &
    utens (:,:,:),  & ! u-tendency without sound-wave terms           ( m/s2)
    vtens (:,:,:),  & ! v-tendency without sound-wave terms           ( m/s2)
    wtens (:,:,:),  & ! w-tendency without sound-wave terms           ( m/s2)
                      ! (defined on half levels)
    ttens (:,:,:),  & ! t-tendency without sound-wave terms           ( K/s )
    qvtens(:,:,:),  & ! qv-tendency                                   ( 1/s )
    qctens(:,:,:),  & ! qc-tendency                                   ( 1/s )
    qitens(:,:,:),  & ! qi-tendency                                   ( 1/s )
    pptens(:,:,:)     ! pp-tendency without sound-wave terms          (Pa/s )


! 5. fields for surface values and soil/canopy model variables        (unit )
! -----------------------------------------------------

  REAL  (KIND=ireals), TARGET, ALLOCATABLE ::           &
    ps       (:,:,:),  & ! surface pressure                           ( pa  )
    t_snow   (:,:,:),  & ! temperature of the snow-surface            (  K  )
    t_snow_mult(:,:,:,:),& ! temperature of the snow-surface       (  K  )
    t_s      (:,:,:),  & ! temperature of the ground surface (soil)   (  K  )
    t_g      (:,:,:),  & ! weighted surface temperature               (  K  )
    qv_s     (:,:,:),  & ! specific water vapor content at the surface(kg/kg)
    t_m      (:,:,:),  & ! temperature between upper and medium
                         ! soil layer                                 (  K  )
    t_cl     (:,:),    & ! temperature between medium and lower
                         ! soil layer (climatology)                   (  K  )
    t_so     (:,:,:,:),& ! multi-layer soil temperature               (  K  )
    w_snow   (:,:,:),  & ! water content of snow                      (m H2O)
    w_i      (:,:,:),  & ! water content of interception water        (m H2O)
    w_g1     (:,:,:),  & ! water content of the upper soil layer      (m H2O)
    w_g2     (:,:,:),  & ! water content of the medium soil layer     (m H2O)
    w_g3     (:,:,:),  & ! water content of the lower soil layer      (m H2O)
                         ! (if nlgw=3, unused otherwise)
    w_so     (:,:,:,:),& ! multi-layer soil moisture                  (m H2O)
    w_so_ice (:,:,:,:),& ! multi-layer soil ice                       (m H2O)
    w_cl     (:,:),    & ! climatological water content               (m H2O)
    freshsnow(:,:),    & ! weighting function indicating 'freshness' of snow
    snow_melt(:,:),    & ! snow melt amount                           (kg/m2)
    rho_snow (:,:,:),  & ! prognostic snow density                    (kg/m3)
    wliq_snow(:,:,:,:),& ! liquid water content in snow              (m H2O)
    wtot_snow(:,:,:,:),& ! total (liquid + solid) water content of snow  (m H2O)
    lev_snow (:,:,:,:),& ! vertical grid in snowpack                     (  m  )
    dzh_snow (:,:,:,:),& ! layer thickness between half levels in snow   (  m  )
    rho_snow_mult (:,:,:,:),& ! prognostic snow density                  (kg/m3)
    h_snow   (:,:,:),  & ! snow height                                (  m  )
    t_e      (:,:,:),  & ! surface temperature of the canopy elements (  K  )
    qv_e     (:,:,:)     ! surface value of qv of the canopy elements (Kg/Kg)

  ! fields for prognostic variables of the lake model FLake or ocean
  ! variables
  REAL  (KIND=ireals), TARGET, ALLOCATABLE ::           &
    fr_ice(:,:),    & ! ice fraction for ocean/lake surfaces          (  -  )
    t_ice (:,:,:),  & ! temperature of ice/water surface              (  K  )
    h_ice (:,:,:)     ! lake/sea ice thickness                        (  m  )

  REAL  (KIND=ireals), TARGET, ALLOCATABLE ::           &
    t_mnw_lk(:,:,:),& ! mean temperature of the water column          (  K  )
    t_wml_lk(:,:,:),& ! mixed-layer temperature                       (  K  )
    t_bot_lk(:,:,:),& ! temperature at the water-bottom sediment
                      ! interface                                     (  K  )
    t_b1_lk (:,:,:),& ! temperature at the bottom of the upper layer
                      ! of the sediments                              (  K  )
    c_t_lk  (:,:,:),& ! shape factor with respect to the
                      ! temperature profile in lake thermocline       (  -  )
    h_ml_lk (:,:,:),& ! thickness of the mixed-layer                  (  m  )
    h_b1_lk (:,:,:)   ! thickness of the upper layer
                      ! of bottom sediments                           (  m  )

! 6. fields that are computed in the parametrization and dynamics     (unit )
! ---------------------------------------------------------------

  REAL  (KIND=ireals), TARGET, ALLOCATABLE ::           &

    qvt_diff(:,:,:),& ! humidity    tendency  due to diffusion        ( 1/s )
    tinc_lh(:,:,:), & ! temperature increment due to latent heat      (  K  )

!   air density at present time level (main levels)       
    rho (:,:,:),    & ! total density of air                          (kg/m3)

!   coefficients for turbulent diffusion in the atmosphere
!   (defined on half levels)
                      ! vertical   turbulent diffusion coefficients
    tkvm(:,:,:),    & ! ... for momentum                              (m2/s)
    tkvh(:,:,:),    & ! ... for heat and moisture                     (m2/s)
                      ! horizontal turbulent diffusion coefficients
    tkhm(:,:,:),    & ! ... for momentum                              (m2/s)
    tkhh(:,:,:),    & ! ... for heat and moisture                     (m2/s)

!   vertical varying implicitness of vertical diffusion
!   (def. at half levels)
    a1t(:),         & !                                               ( -- )
    a2t(:),         & !                                               ( -- )

    ! Rayleigh damping coefficient
    rdcoef(:),      & !

!   turbulence statistics in the atmosphere
!   (defined on full levels)
    rcld(:,:,:),    & ! standard deviation of the saturation deficit    --

!   turbulent coefficients at the surface
    tcm (:,:),      & ! transfer coefficient for momentum             ( -- )
    tch (:,:),      & ! transfer coefficient for heat and moisture    ( -- )
#ifdef COUP_OAS_COS
    tcw (:,:),      & ! transfer coefficient for moisture             ( -- )   !CPS
#endif
    tfm (:,:),      & ! factor of laminar transfer of momentum           --
    tfh (:,:),      & ! factor of laminar transfer of scalars            --
    tfv (:,:),      & ! laminar reduction factor for evaporation         --
 
!   fields from the radiation scheme
    sohr   (:,:,:), & ! rate of solar heating                         ( K/s )
    thhr   (:,:,:), & ! rate of thermal heating                       ( K/s )
    sodwddm(:,:),   & ! downward direct solar radiative flux / smu0   ( W/m2)
    qc_rad (:,:,:), & ! subgrid-scale specific cloud liq. water cont. (kg/kg)
    qi_rad (:,:,:), & ! subgrid-scale specific ice water              (kg/kg)
    clc_sgs(:,:,:), & ! subgrid-scale stratiform cloud cover            --
#ifdef COUP_OAS_COS
    alb_rad(:,:,:), & ! direct and diffuse albedo of the ground
#else
    alb_rad(:,:),   & ! albedo of the ground                            --
#endif
    sobs   (:,:),   & ! solar radiation at the ground                 ( w/m2)
    thbs   (:,:),   & ! thermal radiation at the ground               ( w/m2)
    pabs   (:,:),   & ! photosynthetic active radiation at the ground ( w/m2)
    sobt   (:,:),   & ! solar radiation at the upper boundary         ( w/m2)
                      ! of the atmosphere
    thbt   (:,:),   & ! thermal radiation at the upper boundary       ( w/m2)
                      ! of the atmosphere
    clch   (:,:),   & ! cloud cover with high clouds                    --   
    clcm   (:,:),   & ! cloud cover with medium clouds                  --   
    clcl   (:,:),   & ! cloud cover with low clouds                     --   
    clct   (:,:),   & ! total cloud cover                               --   
    sun_el (:,:),   & ! sun elevation angle                           (deg  )
    sun_azi(:,:)      ! sun azimuth  angle                            (deg  )

  REAL  (KIND=ireals), TARGET, ALLOCATABLE ::           &
!   fields from the convection scheme
    clc_con(:,:,:), & ! cloud cover due to convection                   --     
    clw_con(:,:,:), & ! convective cloud liquid water
    prr_con(:,:),   & ! precipitation rate of rain, convective        (kg/m2*s)
    prs_con(:,:),   & ! precipitation rate of snow, convective        (kg/m2*s)
    prne_con(:,:),  & ! precipitation rate, no evaporat., convective  (kg/m2*s)
    bas_con(:,:),   & ! level index of convective cloud base            --
    top_con(:,:),   & ! level index of convective cloud top             --
    tt_conv (:,:,:),& ! temperature tendency due to convection        ( K/s  )
    qvt_conv(:,:,:),& ! humidity    tendency due to convection        ( 1/s  )
    qct_conv(:,:,:),& ! qc tendency due to convection                 ( 1/s  )
    qit_conv(:,:,:),& ! qi tendency due to convection                 ( 1/s  )
    qrt_conv(:,:,:),& ! qr tendency due to convection                 ( 1/s  )
    qst_conv(:,:,:),& ! qs tendency due to convection                 ( 1/s  )
    ut_conv (:,:,:),& ! u-tendency due to convection                  ( m/s^2)
    vt_conv (:,:,:),& ! v-tendency due to convection                  ( m/s^2)
    mflx_con(:,:),  & ! cloud base massflux                           (kg/m2*s)
    cape_con(:,:),  & ! convective available energy                   (   J/kg)
    tke_con (:,:),  & ! convective turbulent energy                   (   J/kg)
    qcvg_con(:,:),  & ! moisture convergence for Kuo-type closure     (    1/s)
    w0avg   (:,:,:)   ! running average of w

  REAL  (KIND=ireals), TARGET, ALLOCATABLE ::           &
!   fields for Kain-Fritsch/Bechtold convection scheme
    ma_usl_kmin(:,:)   , & ! level index of lower boundary of updraft source layer
    ma_usl_kmax(:,:)   , & ! level index of upper boundary of updraft source layer
    ma_lfs     (:,:)   , & ! level of free sink
    ma_etl_k   (:,:)   , & ! equilibrium temperature level
    ma_ml      (:,:)   , & ! melting level
    ma_ddt     (:,:)   , & ! top of detrainment layer
    ma_usl_buoy(:,:)   , & ! buoyancy at start of ascent (T only)
    ma_nb_k    (:,:)   , & ! level index of neutral buoyancy below the LCL (T only)
    ma_nb_k_min(:,:)   , & ! minimum of ma_nb_k in last hour (T only)
    ma_lcl_k   (:,:)   , & ! level index of lifting condensation level (LCL)
    ma_lcl_t   (:,:)   , & ! parcel temperature at the LCL
    ma_lcl_dt  (:,:)   , & ! temperature perturbation at the LCL
    ma_lcl_tenv(:,:)   , & ! environment temperature at the LCL
    ma_trg     (:,:)   , & ! trigger criteria (degrees)
    ma_trg_max (:,:)   , & ! maximum of ma_trg in last hour
    ma_top_k   (:,:)   , & ! level index of cloud top
    ma_type    (:,:)       ! type of convection (1=deep, 2=shallow, 3=mid-level)

  REAL  (KIND=ireals), TARGET, ALLOCATABLE ::           &
    ma_umf   (:,:,:)   , &  ! updraft mass flux
    ma_udr   (:,:,:)   , &  ! updraft detrainment rate
    ma_uer   (:,:,:)   , &  ! updraft entrainment rate
    ma_urv   (:,:,:)   , &  ! water vapour in updraft
    ma_urci  (:,:,:)   , &  ! total condensat in updraft
    ma_ls_rad(:,:)     , &  ! updraft area
    ma_uw    (:,:)     , &  ! velocity in updraft
    ma_wsub  (:,:,:)   , &  ! compensating mass flux in environment
    ma_dmf   (:,:,:)   , &  ! downdraft mass flux
    ma_der   (:,:,:)   , &  ! downdraft entrainment rate
    ma_ddr   (:,:,:)   , &  ! downdraft detrainment rate
    ma_drw   (:,:,:)   , &  ! total downdraft water
    ma_prlflx(:,:,:)   , &  ! liquid precipitation flux
    ma_prsflx(:,:,:)   , &  ! solid precipitation flux
    ma_urr   (:,:,:)   , &  ! liquid precipitation produced in model layer
    ma_urs   (:,:,:)        ! solid precipitation produced in model layer

  INTEGER (KIND=iintegers), TARGET, ALLOCATABLE ::           &
    nca     (:,:)     !

  REAL  (KIND=ireals), TARGET, ALLOCATABLE ::           &
!   fields from the grid-scale precipitation scheme
    qrs    (:,:,:), & ! precipitation water content (water loading)   (kg/kg)
    prr_gsp(:,:),   & ! precipitation rate of rain, grid-scale        (kg/m2*s)
    prs_gsp(:,:),   & ! precipitation rate of snow, grid-scale        (kg/m2*s)
    prg_gsp(:,:)      ! precipitation rate of graupel, grid-scale     (kg/m2*s)

  REAL  (KIND=ireals), TARGET, ALLOCATABLE ::           &
!   fields from the sub-grid scale orography scheme
    ut_sso (:,:,:), & ! u-tendency due to SSO                         ( m/s2)
    vt_sso (:,:,:), & ! v-tendency due to SSO                         ( m/s2)
    tt_sso (:,:,:), & ! temperature tendency due to SSO               ( K/s )
    ustr_sso (:,:), & ! u-stress (surface momentum flux) due to SSO   ( N/m2)
    vstr_sso (:,:), & ! v-stress (surface momentum flux) due to SSO   ( N/m2)
    vdis_sso (:,:), & ! vert. int. dissipation of kin. en. due to SSO ( W/m2)
    austr_sso(:,:), & ! average of ustr_sso                           ( N/m2)
    avstr_sso(:,:), & ! average of vstr_sso                           ( N/m2)
    avdis_sso(:,:)    ! average of vdis_sso                           ( W/m2)

  REAL  (KIND=ireals), TARGET, ALLOCATABLE ::           &
!   fields for the radiation correction scheme
    ! these are actual values
    swdir_s  (:,:), & ! direct comp. of solar radiative flux at surface ( W/m2)
    swdifd_s (:,:), & ! diffuse downward comp. of short wave rad. flux  ( W/m2)
    swdifu_s (:,:), & ! diffuse upward   comp. of short wave rad. flux  ( W/m2)
    lwd_s    (:,:), & !         downward comp. of long  wave rad. flux  ( W/m2)
    lwu_s    (:,:), & !         upward   comp. of long  wave rad. flux  ( W/m2)

    ! these are accumulated values
    aswdir_s (:,:), & ! direct comp. of solar radiative flux at surface ( W/m2)
    aswdifd_s(:,:), & ! diffuse downward comp. of short wave rad. flux  ( W/m2)
    aswdifu_s(:,:), & ! diffuse upward   comp. of short wave rad. flux  ( W/m2)
    alwd_s   (:,:), & !         downward comp. of long  wave rad. flux  ( W/m2)
    alwu_s   (:,:), & !         upward   comp. of long  wave rad. flux  ( W/m2)

    ! this is the essential correction factor
    swdir_cor(:,:), & ! direct short wave radiation correction factor actual value

    ! these are topographic parameters
    skyview  (:,:), & ! sky view
    slo_asp  (:,:), & ! slope aspect
    slo_ang  (:,:), & ! slope angle
    horizon(:,:,:), & ! horizon

!   fields that are computed in the dynamics
    dqvdt  (:,:,:), & ! threedimensional moisture convergence         ( 1/s )
    qvsflx (:,:),   & ! surface flux of water vapour                  (kg/m2s)
    dpsdt  (:,:),   & ! tendency of the surface pressure              ( pa/s)
    umfl_s (:,:),   & ! u-momentum flux (surface)                     ( N/m2)
    vmfl_s (:,:),   & ! v-momentum flux (surface)                     ( N/m2)
    shfl_s (:,:),   & ! sensible heat flux (surface)                  ( W/m2)
    lhfl_s (:,:),   & ! latent heat flux (surface)                    ( W/m2)
    aumfl_s(:,:),   & ! average u-momentum flux (surface)             ( N/m2)
    avmfl_s(:,:),   & ! average v-momentum flux (surface)             ( N/m2)
    ashfl_s(:,:),   & ! average sensible heat flux (surface)          ( W/m2)
    alhfl_s(:,:),   & ! average latent heat flux (surface)            ( W/m2)
    rstom  (:,:),   & ! stomata resistance                            ( s/m )
    lhfl_bs(:,:),   & ! latent heat flux from bare soil evap.         ( W/m2)
    lhfl_pl(:,:,:), & ! latent heat flux from plants                  ( W/m2)
    alhfl_bs(:,:),  & ! average latent heat flux from bare soil evap. ( W/m2)
    alhfl_pl(:,:,:)   ! average latent heat flux from plants          ( W/m2)

  ! fields used in the Runge-Kutta scheme (only per time step)
  REAL (KIND = ireals), ALLOCATABLE :: &
    wcon(:,:,:),    & ! contravariant vertical velocity
    uadvt(:,:,:),   & ! advective tendency of u
    vadvt(:,:,:),   & ! advective tendency of v
    wadvt(:,:,:),   & ! advective tendency of w
    ppadvt(:,:,:),  & ! advective tendency of pp
    tadvt(:,:,:)      ! advective tendency of t

! 7. fields for model output and diagnostics                          (unit )
! ------------------------------------------

  REAL  (KIND=ireals), TARGET, ALLOCATABLE ::           &
    t_2m    (:,:),  & ! temperature in 2m                             (  K  )
    t_2m_av (:,:),  & ! time mean temperature in 2m                   (  K  )
    qv_2m   (:,:),  & ! specific water vapor content in 2m            (kg/kg)
    td_2m   (:,:),  & ! dew-point in 2m                               (  K  )
    td_2m_av(:,:),  & ! time mean dew-point in 2m                     (  K  )
    rh_2m   (:,:),  & ! relative humidity in 2m                       (  %  )
    u_10m   (:,:),  & ! zonal wind in 10m                             ( m/s )
    u_10m_av(:,:),  & ! time_mean zonal wind in 10m                   ( m/s )
    v_10m   (:,:),  & ! meridional wind in 10m                        ( m/s )
    v_10m_av(:,:),  & ! time mean meridional wind in 10m              ( m/s )
    tmin_2m (:,:),  & ! minimum temperature in 2m                     (  K  )
    tmax_2m (:,:),  & ! maximum temperature in 2m                     (  K  )
    vmax_10m(:,:),  & ! maximal wind gust in 10m                      ( m/s )
    vgust_dyn(:,:), & ! maximal dynamical wind gust in 10m            ( m/s )
    vgust_con(:,:), & ! maximal convective wind gust in 10m           ( m/s )
    asob_s  (:,:),  & ! average solar radiation budget (surface)      ( W/m2)
    athb_s  (:,:),  & ! average thermal radiation budget (surface)    ( W/m2)
    apab_s  (:,:),  & ! average photosynthetic active radiation (sfc) ( W/m2)
    asob_t  (:,:),  & ! average solar radiation budget (model top)    ( W/m2) 
    athb_t  (:,:),  & ! average thermal radiation budget (model top)  ( W/m2)
     sod_t  (:,:),  & ! solar downward radiation at top of atmosphere (     )
    asod_t  (:,:),  & ! averaged solar downward radiation at top      (     )
    dursun  (:,:),  & ! sunshine duration                             (  s  )
    dursun_m(:,:),  & ! maximum possible sunshine duration            (  s  )
    dursun_r(:,:),  & ! relative sunshine duration                    (  s  )
    rain_gsp(:,:),  & ! amount of rain from grid-scale precip. (sum)  (kg/m2)
    snow_gsp(:,:),  & ! amount of snow from grid-scale precip. (sum)  (kg/m2)
    grau_gsp(:,:),  & ! amount of graupel from grid-scale precip. (sum) (kg/m2)
    rain_con(:,:),  & ! amount of rain from convective precip. (sum)  (kg/m2)
    snow_con(:,:),  & ! amount of snow from convective precip. (sum)  (kg/m2)
    runoff_s(:,:),  & ! surface water runoff; sum over forecast       (kg/m2)
    runoff_g(:,:),  & ! soil water runoff; sum over forecast          (kg/m2)
    tdiv_hum(:,:),  & ! vertical sum for  divergence of humidity      (kg/m2)
    aevap_s (:,:)     ! accumulated surface moisture flux             (kg/m2)

! 8. fields for the boundary values                                   (unit )
! ---------------------------------

  REAL  (KIND=ireals), TARGET, ALLOCATABLE ::           &
    u_bd   (:,:,:,:), & ! boundary field for u                        ( m/s )
    v_bd   (:,:,:,:), & ! boundary field for v                        ( m/s )
    w_bd   (:,:,:,:), & ! boundary field for w                        ( m/s )
    t_bd   (:,:,:,:), & ! boundary field for t                        (  K  )
    qv_bd  (:,:,:,:), & ! boundary field for qv                       (kg/kg)
    qc_bd  (:,:,:,:), & ! boundary field for qc                       (kg/kg)
    qi_bd  (:,:,:,:), & ! boundary field for qi                       (kg/kg)
    qr_bd  (:,:,:,:), & ! boundary field for qr                       (kg/kg)
    qs_bd  (:,:,:,:), & ! boundary field for qs                       (kg/kg)
    qg_bd  (:,:,:,:), & ! boundary field for qg                       (kg/kg)
    pp_bd  (:,:,:,:), & ! boundary field for pp                       (  pa )
    qv_s_bd  (:,:,:), & ! boundary field for qv_s                     (kg/kg)
    t_snow_bd(:,:,:), & ! boundary field for t_snow                   (  K  )
    t_s_bd   (:,:,:), & ! boundary field for t_s                      (  K  )
    t_m_bd   (:,:,:), & ! boundary field for t_m                      (  K  )
    w_snow_bd(:,:,:), & ! boundary field for w_snow                   (m H2O)
    w_g1_bd  (:,:,:), & ! boundary field for w_g1                     (m H2O)
    w_g2_bd  (:,:,:), & ! boundary field for wb2                      (m H2O)
    w_g3_bd  (:,:,:), & ! boundary field for wb3                      (m H2O)

    ! and for the CLM Version
    hmo3_bd  (:,:,:), & ! boundary field for hmo3                     (m    )
    vio3_bd  (:,:,:), & ! boundary field for vio3                     (pa O3)
    w_cl_bd  (:,:,:), & ! boundary field for w_cl                     (m H2O)
    t_cl_bd  (:,:,:), & ! boundary field for t_cl                     (  K  )
    lai_bd   (:,:,:), & ! boundary field for lai                      ( --  )
    rootdp_bd(:,:,:), & ! boundary field for rootdp                   (m    )
    plcov_bd (:,:,:)    ! boundary field for plcov                    ( --  )

! 9. fields for the synthetic satellite images
! --------------------------------------------

  REAL  (KIND=ireals), TARGET, ALLOCATABLE ::           &
    synme5 (:,:,:)  , & ! Meteosat 5
    synme6 (:,:,:)  , & ! Meteosat 6
    synme7 (:,:,:)  , & ! Meteosat 7
    synmsg (:,:,:)      ! Meteosat Second Generation

! 10. analysis increment fields
! -----------------------------

  REAL  (KIND=ireals), TARGET, ALLOCATABLE ::           &
    ff_anai (:,:,:) , & ! wind velocity                               ( m/s )
    dd_anai (:,:,:) , & ! wind direction                              ( rad )
    t_anai  (:,:,:) , & ! temperature                                 (  k  )
    p_anai  (:,:,:) , & ! deviation from the reference pressure       ( Pa  )
    qv_anai (:,:,:) , & ! specific water vapor content                (kg/kg)
    qc_anai (:,:,:)     ! specific cloud water content (via saturation adjustm)
!   fi_anai (:,:,:) , & ! geopotential
!   pmsl_anai (:,:) , & ! mean sea level pressure
!   tqv_anai  (:,:) , & ! wind velocity
!   tqc_anai  (:,:)     ! wind velocity

!==============================================================================

END MODULE data_fields
