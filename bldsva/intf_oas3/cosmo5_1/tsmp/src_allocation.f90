!+ Module for the allocation/deallocation of meteorological LM fields
!------------------------------------------------------------------------------

MODULE src_allocation

!------------------------------------------------------------------------------
!
! Description:
!   This module performs the allocation and deallocations of the meteorological
!   fields of the model which are declared in the data module data_fields.
!
! Current Code Owner: DWD, Ulrich Schaettler
!  phone:  +49  69  8062 2739
!  fax:    +49  69  8062 3721
!  email:  ulrich.schaettler@dwd.de
!
! History:
! Version    Date       Name
! ---------- ---------- ----
! 1.34       1999/12/10 Ulrich Schaettler
!  Initial release
! 1.39       2000/05/03 Ulrich Schaettler
!  Put subroutine constant_fields to src_setup and introduced subroutine
!  organize_allocation (needed for the nesting version).
! 2.2        2000/08/18 Matthias Raschendorfer
!  Allocation of the 2-D fields 'eai' and 'tai'.
! 2.4        2001/01/29 Christoph Schraff
!  Allocation of the 2-D field 'prne_con' added.
! 2.8        2001/07/06 Ulrich Schaettler
!  Added new fields for multi-layer soil model
! 2.11       2001/09/28 Ulrich Schaettler
!  Some corrections for allocating fields for the multi-layer soil model
!  Allocation of boundary fields for cloud ice
! 2.17       2002/05/08 Ulrich Schaettler
!  Allocation of new fields for the Kain-Fritsch convection scheme
! 2.18       2002/07/16 Ulrich Schaettler
!  Allocation and deallocation of new variables qr and qs
! 3.5        2003/09/02 Ulrich Schaettler
!  Allocation and deallocation of new variables phi_tot, rla_tot
! 3.6        2003/12/11 Reinhold Schrodin
!  Allocation and deallocation of variables for multi-layer soil model
! 3.7        2004/02/18 Ulrich Schaettler
!  Renamed alb (alb_rad), idiv_hum (tdiv_hum), phi (rlat), rla (rlon)
!  Renamed cphi (crlat), acphir (acrlat), tgphi (tgrlat)
! 3.13       2004/12/03 Ulrich Schaettler
!  Use and allocation/deallocation of new fields for 
!  multi-layer soil model (Reinhold Schrodin), Graupel scheme (Thorsten 
!  Reinhardt) and Runge-Kutta scheme (Jochen Foerstner)
! 3.14       2005/01/25 Ulrich Schaettler
!  Corrections for allocation of soil variables for old soil model
! 3.16       2005/07/22 Reinhold Schrodin
!  Allocation and deallocation of new variables for_e and for_d
! 3.17       2005/12/12 Reinhold Schrodin
!  Add prognostic snow density rho_snow and snow heigth h_snow
! 3.18       2006/03/03 Ulrich Schaettler / Klaus Stephan
!  Added fields for CLM, LHN, FLake
!  LHN namelist parameter moved to data_lheat_nudge to avoid to many dependencies
!  Added field for 2m relative humidity (rh_2m) and temperature increment due
!  to latent heat (tinc_lh)
! 3.19       2006/04/25 Ulrich Schaettler
!  Added field T_S_LAKE for saving the lake temperatures in the Nudging cycle
! 3.21       2006/12/04 Ulrich Schaettler, et.al.
!  Renamed sunshhrs, sodwdir to dursun, sodwddm (Burkhardt Rockel)
!  Added fields for deep atmosphere (Ronny Petrik)
!  Additional fields for KainFritsch/Bechtold convection schemes (MeteoSwiss)
!  Time integrated analysis increment fields introduced (Christoph Schraff)
! V3_23        2007/03/30 Matthias Raschendorfer, Matteo Buzzi, Jochen Foerstner
!  Added field 'tfv' containing the laminar reduction factor for evaporation
!  Added fields for topographic radiation correction (Matteo Buzzi)
!  Added fields for new relaxation (Jochen Foerstner)
! V3_24        2007/04/26 Ulrich Schaettler
!  Added prints for verbosity of output
! V3_26        2007/05/29 Ulrich Schaettler
!  More debug level output
! V4_3         2008/02/25 Matthias Raschendorfer
!  Introductinon of a 3D diagnostic field 'edr' for the eddy dissipotion rate
! V4_4         2008/07/16 Jan-Peter Schulz
!  Added external parameters (sso_stdh, sso_gamma, sso_theta, sso_sigma) and
!  fields (ut_sso, vt_sso, tt_sso, ustr_sso, vstr_sso, vdis_sso,
!  austr_sso, avstr_sso, avdis_sso) for SSO scheme
!  Replaced lkainfri by itype_conv (Ulrich Schaettler)
! V4_8         2009/02/16 Ulrich Schaettler
!  Added fields for convective and dynamical gust in 10m
!  Added fields for additional output of radiation values (for CLM)
!  Use p0hl (reference pressure at half levels) for full consistency with
!  the new reference atmosphere implementation (Guenther Zaengl)
! V4_9         2009/07/16 Ulrich Schaettler
!  Allocate additional mask fields for t-, q- and u-fields
!  Allocate sqrtg_r_x fields in any case
! V4_10        2009/09/11 Ulrich Schaettler, Jan-Peter Schulz
!  Allocation / Deallocation of snow_melt
!  Allocation / Deallocation of sea-ice fields
! V4_11        2009/11/30 Ekaterina Machulskaya, Juergen Helmert, Lucio Torrisi
!  Allocation / Deallocation of fields for multi layer snow model (EM)
!  Allocation / Deallocation of additional fields for radiation / soil (JH)
!  Allocation / Deallocation of additional fields for output (LT)
! V4_12        2010/05/11 Michael Baldauf, Ulrich Schaettler
!  Allocate dzeta_dlam, dzeta_dphi
!  Renamed hd_mask_dcoeff to hd_mask_dcoeff_p for diffusion of pressure
!  Allocation of additional fields for sunshine duration (Oli Fuhrer)
! V4_13        2010/05/11 Michael Gertz
!  Adaptions to SVN
! V4_15        2010/11/19 Ulrich Schaettler
!  The field t_s_lake is removed again after adaptations in the SST Analysis
! V4_17        2011/02/24 Ulrich Blahak
!  Fixed allocation of rcld and edr in case of itype_tran=2;
!  Fixec allocation of tke and tketens for the case lphys=.false. and
!  at the same time RK-scheme iype_turb between 5 and 8 as an intermediate solution. 
!  This is a pathological case, but
!  if tke and tketens are not allocated, the model wants to diffuse and transport
!  them anyways, because lphys=.false. is not implemented correctly into the RK-core. 
!  Fixed initialization of dzeta_dphi and dzeta_dlam (cut-and-paste error)
! V4_19        2011/08/01 Ulrich Schaettler
!  Introduced conditional compilation for Nudging
! V4_20        2011/08/31 Matthias Raschendorfer
!  Introduction of additional 3D-arrays 'tket_(conv, sso, hshr)' for additional
!   TKE source by the action of other sub grid scale flow patterns.
!  Removed action=canopy fully and move allocation of canopy fields to init_canopy
!   from src_turbdiff
!  tgrlat needs 2 dimensions because of v-point dependence
! V4_22        2012/01/31 Thorsten Reinhardt
!  Additional fields used to update the solar zenith angle in the radiation.
! V4_23        2012/05/10 Michael Baldauf, Ulrich Schaettler, Oliver Fuhrer, CLM
!                         Lucio Torrisi
!  Eliminate (de)allocation of sqrtg_r_*, dzeta_*; eliminate l_dzeta_d_needed
!  Allocate some special fields anyhow, not depending on settings, because they
!    are used in interfaces: depth_lk, rcld, edr, ut_sso, vt_sso
!  Eliminated switch lprogprec (US)
!  Eliminated field qvt_diff (OF)
!  Introduction of fields for surface albedo (alb_dry, alb_sat: CLM; alb_dif: JH)
!  Introduction of diagnostic variable vabsmx_10m (max wind speed in 10 meters) (CLM)
!  Implementation of time dependent boundary values for aerosol optical depths (CLM)
!  Take care that t_s_bd is allocated, if lbdsst=.TRUE. (LT)
! V4_24        2012/06/22 Michael Baldauf
!  Added treatment of t0hl
! V4_25        2012/09/28 Anne Roches, Oliver Fuhrer, Ulrich Blahak,
!                         Ekaterina Machulskaya
!  Removed allocation/deallocation of qx and their tendencies and boundary variables
!  Added tracer variables for the 2mom microphysics scheme. (UB)
!  Added a new field tg_radstep: ground temperature at the last call to radiation
! V4_26        2012/12/06 Ulrich Schaettler
!  Adapted variable names of multi-layer snow model to corresponding
!   short names for I/O
! V4_27        2013/03/19 Astrid Kerkweg, Ulrich Schaettler
!  MESSy interface introduced
! V4_28        2013/07/12 Ulrich Schaettler
!  Eliminated use of vcoord, sigmr, hhlr (now used from vgrid_refatm_utils)
! V4_29        2013/10/04 Ulrich Schaettler, Astrid Kerkweg, Ulrich Blahak
!  Added 3 global fields: fr_snow, fr_wi, ustar_fv (which were local to
!   the soil model before): for a better handling in MESSY
!  Allocate memory for synthetic satellite images here for long term storage
!  Allocation of tkvh, tkvm, tkhm, tkhh for 1:ke1 in the vertical
!  Added some missing lalloc_t_s_bd=.true. and  a missing l_alloc_t_cl=.true. (UB)
! V5_1         2014-11-28 Ulrich Schaettler, Ulrich Blahak, Oliver Fuhrer
!  Allocate fr_lake anyhow, because it is used during NetCDF I/O
!  Added a 1D field hhl_prof, to save a special hhl_profile (to replace vert_coord) (US)
!  Removed unnecessary TKE allocation catch in case of itpye_turb 5,6,7,8. This
!   is now fixed by checking lprog_tke in organize_physics more thoroughly.
!  Added allocation of some fields for 2-moment scheme (KIT)
!  Added extra advective tendency of TKE for itype_turb=3.
!  Added "_wp" to all real variables (OF)
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
!
USE data_parameters, ONLY :   &
    wp,        & ! KIND-type parameter for real variables
    iintegers    ! KIND-type parameter for standard integer variables

!------------------------------------------------------------------------------

USE data_modelconfig, ONLY :   &

! 1. vertical coordinate parameters and related variables
! -------------------------------------------------------

    hhl_prof,     & ! a special hhl_profile

! 2. horizontal and vertical sizes of the fields and related variables
! --------------------------------------------------------------------

    ie,           & ! number of grid points in zonal direction
    ie_tot,       & ! the same for the total domain
    je,           & ! number of grid points in meridional direction
    ke,           & ! number of grid points in vertical direction
    ke_soil,      & ! number of layers in the multi-layer soil model
    ke_snow,      & ! number of layers in the multi-layer snow model
    ke1,          & ! KE+1
    kcm,          & ! index of the lowest model layer higher than the canopy

! 7. Layer index corresponding to a specified pressure
!    and other organizational variables
! ----------------------------------------------------

    lalloc_h_ice,   & !
    lalloc_w_g3,    & !
    lalloc_t_cl,    & !
    lalloc_t_s_bd,  & !
    lalloc_w_g3_bd, & !
    lalloc_prr_gsp, & !
    lalloc_prr_con, & !
    lalloc_prs_gsp, & !
    lalloc_prs_con, & !
    lalloc_prg_gsp, & !
    lalloc_tke

!------------------------------------------------------------------------------

USE data_runcontrol , ONLY :   &

! 3. controlling the physics
! --------------------------
    itype_gscp,   & ! type of grid-scale precipitation physics
    itype_turb,   & ! type of turbulent diffusion parametrization
    itype_tran,   & ! type of turbulent diffusion parametrization
    itype_aerosol,& ! type of aerosol map
    lemiss,       & ! surface emissivity map
    nradcoarse,   & ! number of horiz. gridpoints for radiation on coarser grid
    lstomata,     & ! minimum stomata resistance
    nlgw,         & ! number of prognostic soil water levels
    nlgw_ini,     & ! number of prognostic soil water levels in initial data
    nlgw_bd,      & ! number of prognostic soil water levels in boundary data
    ltur,         & ! forecast with vertical diffusion
    l3dturb,      & ! 3D-turbulence (additional horizontal diffusion)
    lprog_tke,    & ! prognostic treatment of TKE (for itype_turb=5/7)
    lphys,        & ! forecast with physical parametrizations
    lseaice,      & ! forecast with sea ice model
    llake,        & ! forecast with lake model
    lsso,         & ! forecast with sub-grid scale orography scheme
    lforest,      & ! if .true., run with forest (evergreen and deciduous)
    lgsp,         & ! forecast with grid scale precipitation
    lmulti_layer, & ! run multi-layer soil model
    lmulti_snow,  & ! run multi-layer snow model
    itype_albedo, & ! type of surface albedo treatment

! 4. controlling the dynamics
! ---------------------------

    lcori_deep,    & ! if =.TRUE.: take cos(phi) coriolis terms into account

! 7. additional control variables
! -------------------------------

    l2tls,        & ! time integration by two timelevel RK-scheme (.TRUE.)
                    ! or by default three-time level KW-scheme (.FALSE.)
    ldiabf_lh,    & ! include diabatic forcing due to latent heat in RK-scheme
    lw_freeslip,  & ! if .TRUE.: with free slip lateral boundary condition and
                    ! if .FALSE. specified lateral boundary values for w
    lreproduce,   & ! the results are reproducible in parallel mode
    lradtopo,     & ! if .TRUE., calculate topographic correction of radiation
    luse_rttov,   & ! if rttov-library is used

! 8. diagnostic calculations
! --------------------------

    ldiagnos,     & ! perform diagnostic calculations

! 9. Other variables
! ------------------

    lout_anai,    & ! allocate fields to enable writing analysis increments
    nhori,        & ! number of sectors for the horizont array by the topographic
                    ! correction of the radiation

! 12. controlling verbosity of debug output
! -----------------------------------------

    idbg_level,   & ! to control the verbosity of debug output
                    ! Some basic output is written anyhow
                    ! but for most components the output must be
                    ! activated with the logical switches below
    lprintdeb_all   ! .TRUE.:  all tasks print debug output
                    ! .FALSE.: only task 0 prints debug output

!------------------------------------------------------------------------------

USE data_fields     , ONLY :   &

! 1. constant fields for the reference atmosphere                      (unit)
! -----------------------------------------------

    rho0       ,    & ! reference density at the full model levels     (kg/m3)
    dp0        ,    & ! pressure thickness of model layers             ( Pa  )
    p0         ,    & ! reference pressure at full levels              ( Pa  )
    p0hl       ,    & ! reference pressure at falf levels              ( Pa  )
    t0hl       ,    & ! reference temperature at half levels           ( K   )
    dt0dz      ,    & ! temperature gradient of reference atmosphere   ( K/m )
    t0         ,    & ! reference temperature                          ( K   )
    hhl        ,    & ! geometical height of half levels               (  m  )
    hsurf      ,    & ! geometical heigt of surface topography         (  m  )
    sso_stdh   ,    & ! standard deviation of sub-grid scale orography ( m   )
    sso_gamma  ,    & ! anisotropy of sub-grid scale orography           --
    sso_theta  ,    & ! angle betw. principal axis of orography and E  ( rad )
    sso_sigma  ,    & ! mean slope of sub-grid scale orography           --
    aer_su     ,    & ! monthly aerosol climatology sulfate drops               (0 - 1)
    aer_du     ,    & ! monthly aerosol climatology total dust                  (0 - 1)
    aer_or     ,    & ! monthly aerosol climatology organic (water sol.)        (0 - 1)
    aer_bc     ,    & ! monthly aerosol climatology black carbon                (0 - 1)
    aer_ss     ,    & ! monthly aerosol climatology sea salt                    (0 - 1)
    emis_rad   ,    & ! surface emissivity                                      (0 - 1)
    rsmin2d    ,    & ! minimum stomata resistance                              ( s/m )
    swi        ,    & ! soil wetness index                                      (0 - 1)
    gz0        ,    & ! surface roughness  * g                         (m2/s2)
    fr_land    ,    & ! fraction of land in a grid element             (  -- )
    soiltyp    ,    & ! type of the soil (keys 0-9)                    (  -- )
    vio3       ,    & ! vertical integrated ozone contents             (pa O3)
    hmo3       ,    & ! ozone maximum                                  ( pa  )
    rlat       ,    & ! geographical latitude                          ( rad )
    rlon       ,    & ! geographical longitude                         ( rad )
    rlattot    ,    & ! geographical latitude                          ( rad )
    rlontot    ,    & ! geographical longitude                         ( rad )
    fc         ,    & ! coriolis-parameter                             ( 1/s )
    fccos      ,    & ! coriolis-parameter mit cos(phi)                ( 1/s )
    rmy        ,    & ! Davis-parameter for relaxation (mass, qv, qc)    --
    rmyq              ! Davis-parameter for relaxation (qr, qs, qg)      --

USE data_fields     , ONLY :   &

! 2. external parameter fields                                         (unit)
! ----------------------------

    hd_mask_dcoeff_p, & ! 3D-domain mask for horizontal diffusion * dcoeff --
    hd_mask_dcoeff_t, & ! 3D-domain mask for horizontal diffusion * dcoeff --
    hd_mask_dcoeff_q, & ! 3D-domain mask for horizontal diffusion * dcoeff --
    hd_mask_dcoeff_u, & ! 3D-domain mask for horizontal diffusion * dcoeff --
    ofa_hdx    ,    & ! 
    ofa_hdy    ,    & ! 
    hd_mask    ,    & ! 3D-domain mask for horizontal diffusion * dcoeff --
    least_lbdz ,    & ! mask for eastern  lateral boundary zone
    lwest_lbdz ,    & ! mask for western  lateral boundary zone
    lnorth_lbdz,    & ! mask for northern lateral boundary zone
    lsouth_lbdz,    & ! mask for southern lateral boundary zone
    crlat      ,    & ! cosine of transformed latitude
    acrlat     ,    & ! 1 / ( crlat * radius of the earth )            ( 1/m )
    tgrlat     ,    & ! tangens of transformed latitude                  --
    aerlan     ,    & ! aerosol-distribution for rural areas             --
    aerurb     ,    & ! aerosol-distribution for urban areas             --
    aerdes     ,    & ! aerosol-distribution for desert areas            --
    aersea     ,    & ! aerosol-distribution for sea                     --
    plcov      ,    & ! fraction of plant cover                          --
    lai        ,    & ! leaf area index of plants                        --
    tai        ,    & ! transpiration area index                         --
    sai        ,    & ! surface area index                               --
    eai        ,    & ! (evaporative) earth area index                   --
    rootdp     ,    & ! depth of the roots                             (  m  )
    llandmask  ,    & ! landpoint mask
    lseamask          ! ocean point mask, i.e. water but not lake

USE data_fields     , ONLY :   &
    for_e      ,    & ! ground fraction covered by evergreen forest      --
    for_d      ,    & ! ground fraction covered by deciduous forest      --
    h_can      ,    & ! hight of the vertically resolved canopy        (  m  )
    d_pat      ,    & ! horizontal pattern length scale                (  m  )
    c_big      ,    & ! effective drag coefficient of canopy elements
                      ! larger than or equal to the tubulent length
                      ! scale                                          ( 1/m )
    c_sml      ,    & ! effective drag coefficient of canopy elements
                      ! smaller than the tubulent length scale         ( 1/m )
    r_air      ,    & ! air containing fraction of a gridbox inside
                      ! the canopy                                     (  1  )
    fr_lake    ,    & ! lake fraction in a grid element [0,1]          (  -- )
    depth_lk   ,    & ! lake depth                                     (  m  )
    fetch_lk   ,    & ! wind fetch over lake                           (  m  )
    dp_bs_lk   ,    & ! depth of the thermally active layer
                      ! of bottom sediments                            (  m  )
    t_bs_lk    ,    & ! climatological temperature at the bottom of
                      ! the thermally active layer of sediments        (  K  )
    gamso_lk   ,    & ! attenuation coefficient for
                      ! solar radiation in lake water                  ( 1/m )
    alb_dry    ,    & ! surface albedo field for dry soil
    alb_sat    ,    & ! surface albedo field for saturated soil
    alb_dif           ! diffuse albedo field 

USE data_fields     , ONLY :   &

! 3. prognostic variables                                              (unit)
! -----------------------

    u          ,    & ! zonal wind speed                               ( m/s )
    v          ,    & ! meridional wind speed                          ( m/s )
    w          ,    & ! vertical wind speed (defined on half levels)   ( m/s )
    t          ,    & ! temperature                                    (  k  )
    pp         ,    & ! deviation from the reference pressure          ( pa  )

! Fields of the turbulence scheme defined on half levels:

    tke        ,    & ! SQRT(2*TKE); TKE='turbul. kin. energy'        ( m/s )
    tketens    ,    & ! non-advective tendency of SQRT(2*TKE)         ( m/s2)
    tket_adv   ,    & ! pure advective tendency of SQRT(2*TKE)        ( m/s2)
    edr        ,    & ! eddy dissipation rate of TKE (EDR)            ( m2/s3)
    tket_conv  ,    & ! TKE-tendency due to convective buoyancy       ( m2/s3)
    tket_hshr  ,    & ! TKE-tendency due to (sep.) horiz. shear       ( m2/s3)
    tket_sso   ,    & ! TKE-tendency due to SSO wake production       ( m2/s3)

! 4. tendency fields for the prognostic variables                     (unit )
! -----------------------------------------------
!    timely deviation  by diabatic and adiabatic processes 
!    without sound-wave terms

    utens        ,  & ! u-tendency without sound-wave terms            ( m/s2)
    vtens        ,  & ! v-tendency without sound-wave terms            ( m/s2)
    wtens        ,  & ! w-tendency without sound-wave terms            ( m/s2)
                        ! (defined on half levels )
    ttens        ,  & ! t-tendency without sound-wave terms            ( m/s2)
    pptens            ! pp-tendency without sound-wave terms           ( m/s2)

USE data_fields     , ONLY :   &

! 5. fields for surface values and soil model variables                (unit )
! -----------------------------------------------------

    ps        ,     & ! surface pressure                               ( pa  )
    t_snow    ,     & ! temperature of the snow-surface                (  k  )
    t_snow_mult,    & ! temperature of the snow-surface                (  k  )
    dzh_snow_mult,  & !
    w_snow_mult,    & ! total water content of snow                    (m H2O)
    wliq_snow ,     & ! liquid water content of snow                   (m H2O)
    t_s       ,     & ! temperature of the ground surface              (  k  )
    t_g       ,     & ! weighted surface temperature                   (  k  )
    tg_radstep,     & ! ground temperature at the last call            (  K  )
                      ! of the radiation routine
    qv_s      ,     & ! specific water vapor content on the surface    (kg/kg)
    t_m       ,     & ! temperature between upper and medium 
                      ! soil layer                                     (  k  )
    t_cl      ,     & ! temperature between medium and lower
                      ! soil layer                                     (  k  )
    t_so      ,     & ! multi-layer soil temperature                   (  k  )

    w_snow    ,     & ! water content of snow                          (m H2O)
    w_i       ,     & ! water content of interception water            (m H2O)
    w_g1      ,     & ! water content of the upper soil layer          (m H2O)
    w_g2      ,     & ! water content of the medium soil layer         (m H2O)
    w_g3      ,     & ! water content of the lower soil layer          (m H2O)
                      ! (if NLWB=3, unused otherwise)
    w_so      ,     & ! multi-layer soil moisture                      (m H2O)
    w_so_ice  ,     & ! multi-layer soil ice                           (m H2O)
    w_cl      ,     & ! climatological water content                   (m H2O) 
    freshsnow ,     & ! weighting function indicating 'freshness' of snow in
                      ! upper few centimeters of snow cover            ( -- )
    snow_melt ,     & ! snow melt amount                               (kg/m2)
    rho_snow  ,     & ! prognostic density of snow                     (kg/m3)
    rho_snow_mult,  & ! prognostic density of snow                     (kg/m3)
    h_snow    ,     & ! snow height                                    (  m  )
    fr_snow   ,     & ! surface fraction covered by snow               (  -  )
    fr_wi     ,     & ! surface fraction covered by interception water (  -  )
    ustar_fv  ,     & ! friction velocity (ustar)                      ( m/s )
    t_e       ,     & ! surface temperature of the canopy elements     (  k  )
    qv_e      ,     & ! surface value of qd of the canopy elements     (Kg/Kg)
    t_ice     ,     & ! temperature at the snow-ice or
                      ! air-ice interface                              (  K  )
    t_mnw_lk  ,     & ! mean temperature of the water column           (  K  )
    t_wml_lk  ,     & ! mixed-layer temperature                        (  K  )
    t_bot_lk  ,     & ! temperature at the water-bottom sediment
                      ! interface                                      (  K  )
    t_b1_lk   ,     & ! temperature at the bottom of the upper layer
                      ! of the sediments                               (  K  )
    c_t_lk    ,     & ! shape factor with respect to the
                      ! temperature profile in lake thermocline        (  -  )
    h_ice     ,     & ! ice thickness                                  (  m  )
    h_ml_lk   ,     & ! thickness of the mixed-layer                   (  m  )
    h_b1_lk           ! thickness of the upper layer
                      ! of bottom sediments                            (  m  )

!------------------------------------------------------------------------------

USE data_fields     , ONLY :   &

! 6. fields that are computed in the parametrization and dynamics      (unit )
! ---------------------------------------------------------------

    tinc_lh    ,    & ! temperature increment due to latent heat       (  K  )

!   density of moist air 
    rho        ,    & ! total density of moist air                     (kg/m3)

!   coefficients for turbulent diffusion in the atmosphere
!   (defined on half levels)
                      ! vertical   turbulent diffusion coefficients
    tkvm     ,      & ! ... for momentum                               (m2/s)
    tkvh     ,      & ! ... for heat and moisture                      (m2/s)
                      ! horizontal turbulent diffusion coefficients
    tkhm     ,      & ! ... for momentum                               (m2/s)
    tkhh     ,      & ! ... for heat and moisture                      (m2/s)

!   turbulence statistics in the atmosphere
!   (defined on full levels)
    rcld       ,    & ! standard deviation of the saturation deficit      --

!   turbulent coefficients at the surface 
    tcm      ,      & ! transfer coefficient for momentum               ( -- )
    tch      ,      & ! transfer coefficient for heat and moisture      ( -- )
#ifdef COUP_OAS_COS
    tcw      ,      & ! transfer coefficient for          moisture      ( -- )
    trad_clm ,      & ! radiative temperature in CLM as bottom BC       ( K  )
#endif
    tfm      ,      & ! factor of laminar transfer of momentum            --
    tfh      ,      & ! factor of laminar transfer of scalars             --
    tfv      ,      & ! laminar reduction factor for evaporation        ( -- )

!   fields from the radiation scheme
    sohr      ,     & ! rate of solar heating                          ( K/s )
    sotr,sotr_par,  & ! solar transmissivity
    thhr      ,     & ! rate of thermal heating                        ( K/s )
    clc_sgs   ,     & ! subgrid-scale stratiform cloud cover              --
    alb_rad   ,     & ! albedo of the ground                              --
    alb_rad_coarse, & ! albedo of the ground, coarse grid if nradcoarse>1  --
    sobs      ,     & ! solar radiation at the ground                  ( w/m2)
    thbs      ,     & ! thermal radiation at the ground                ( w/m2)
    pabs      ,     & ! photosynthetic active radiation at the ground  ( w/m2)
    sobt      ,     & ! solar radiation at the upper boundary          ( w/m2)
                      ! of the atmosphere
    thbt      ,     & ! thermal radiation at the upper boundary        ( w/m2)
                      ! of the atmosphere
    clch      ,     & ! cloud cover with high clouds                      --   
    clcm      ,     & ! cloud cover with medium clouds                    --   
    clcl      ,     & ! cloud cover with low clouds                       --   
    clct      ,     & ! total cloud cover                                 --   
    sun_el    ,     & ! sun elevation angle                            ( deg )
    sun_azi   ,     & ! sun azimuth  angle                             ( deg )

    ! and used in the Climate-LM Version
    sodwddm   ,     & ! downward direct solar radiative flux / smu0    ( W/m2)
    qc_rad    ,     & ! subgrid-scale specific cloud water             (kg/kg)
    qi_rad            ! subgrid-scale specific ice water               (kg/kg)

USE data_fields     , ONLY :   &
    !   fields for the radiation correction scheme
    swdir_s  ,      & ! direct comp. of solar radiative flux at surface ( W/m2)
    swdifd_s ,      & ! diffuse downward comp. of short wave rad. flux  ( W/m2)
    swdifu_s ,      & ! diffuse upward   comp. of short wave rad. flux  ( W/m2)
    swtrdir_s   ,   & ! direct comp. of solar radiative transmiss. at surface
    swtrdifd_s  ,   & ! diffuse downward comp. of short wave rad. transmiss.
    swtrdifu_s  ,   & ! diffuse upward   comp. of short wave rad. transmiss.
    lwd_s    ,      & !         downward comp. of long  wave rad. flux  ( W/m2)
    lwu_s    ,      & !         upward   comp. of long  wave rad. flux  ( W/m2)

    ! these are accumulated values
    aswdir_s ,      & ! direct comp. of solar radiative flux at surface ( W/m2)
    aswdifd_s,      & ! diffuse downward comp. of short wave rad. flux  ( W/m2)
    aswdifu_s,      & ! diffuse upward   comp. of short wave rad. flux  ( W/m2)
    alwd_s   ,      & !         downward comp. of long  wave rad. flux  ( W/m2)
    alwu_s   ,      & !         upward   comp. of long  wave rad. flux  ( W/m2)

    ! this is the essential correction factor
    swdir_cor,      & ! direct short wave radiation correction factor actual value

    ! these are topographic parameters
    skyview  ,      & ! sky view
    slo_asp  ,      & ! slope aspect
    slo_ang  ,      & ! slope angle
    horizon           ! horizon

!------------------------------------------------------------------------------

USE data_fields     , ONLY :   &

!   fields from the convection scheme
    clc_con    ,    & ! cloud cover due to convection                     --
    clw_con    ,    & ! cloud liquid water due to convection              --
    prr_con    ,    & ! precipitation rate of rain, convective        (kg/m2s)
    prs_con    ,    & ! precipitation rate of snow, convective        (kg/m2s)
    prne_con   ,    & ! precipitation rate, no evaporat., convective  (kg/m2s)
    bas_con    ,    & ! level index of convective cloud base            --
    top_con    ,    & ! level index of convective cloud top             --
    tt_conv    ,    & ! temperature tendency due to convection        ( K/s  )
    qvt_conv   ,    & ! humidity    tendency due to convection        ( 1/s  )
    qct_conv   ,    & ! qc-tendency tendency due to convection        ( 1/s  )
    qit_conv   ,    & ! qi-tendency tendency due to convection        ( 1/s  )
    qrt_conv   ,    & ! qr-tendency tendency due to convection        ( 1/s  )
    qst_conv   ,    & ! qs-tendency tendency due to convection        ( 1/s  )
    ut_conv    ,    & ! u-tendency due to convection                  ( m/s^2)
    vt_conv    ,    & ! v-tendency due to convection                  ( m/s^2)
    mflx_con   ,    & ! convective massflux                           (kg/m2s)
    cape_con   ,    & ! convective available energy                   ( J/kg )
    tke_con    ,    & ! convective turbulent kinetic energy           ( J/kg )
    qcvg_con   ,    & ! moisture convergence for Kuo-type closure     ( 1/s  )
    w0avg      ,    &

!   fields of the precipitation
#ifdef TWOMOM_SB
    reffc_out ,     & ! effective radius of cloud droplets            ( m    )
    reffi_out ,     & ! effective radius of ice particles             ( m    )
    odepthw_so ,    & ! optical depth of cloud droplets
    odepthi_so ,    & ! optical depth of ice particles
    odepthw_th ,    & ! optical depth of cloud droplets
    odepthi_th ,    & ! optical depth of ice particles
#endif
    qrs        ,    & ! precipitation water (water loading)           (kg/kg )
    prr_gsp    ,    & ! precipitation rate of rain, grid-scale        (kg/m2s)
    prs_gsp    ,    & ! precipitation rate of snow, grid-scale        (kg/m2s)
    prg_gsp    ,    & ! precipitation rate of graupel, grid-scale     (kg/m2s)
    prh_gsp           ! precipitation rate of hail, grid-scale        (kg/m2s)

USE data_fields     , ONLY :   &

!   fields from the sub-grid scale orography scheme
    ut_sso     ,    & ! u-tendency due to SSO                         ( m/s2)
    vt_sso     ,    & ! v-tendency due to SSO                         ( m/s2)
    tt_sso     ,    & ! temperature tendency due to SSO               ( K/s )
    ustr_sso   ,    & ! u-stress (surface momentum flux) due to SSO   ( N/m2)
    vstr_sso   ,    & ! v-stress (surface momentum flux) due to SSO   ( N/m2)
    vdis_sso   ,    & ! vert. int. dissipation of kin. en. due to SSO ( W/m2)
    austr_sso  ,    & ! average of ustr_sso                           ( N/m2)
    avstr_sso  ,    & ! average of vstr_sso                           ( N/m2)
    avdis_sso  ,    & ! average of vdis_sso                           ( W/m2)

!   fields that are computed in the dynamics
    dqvdt      ,    & ! threedimensional moisture convergence         ( 1/s  )
    qvsflx     ,    & ! surface flux of water vapour                  (kg/m2s)
    dpsdt      ,    & ! tendency of the surface pressure              ( pa/s )
    umfl_s     ,    & ! u-momentum flux (surface)                     ( N/m2 )
    vmfl_s     ,    & ! v-momentum flux (surface)                     ( N/m2 )
    shfl_s     ,    & ! sensible heat flux (surface)                  ( W/m2 )
    lhfl_s     ,    & ! latent heat flux (surface)                    ( W/m2 )
    aumfl_s    ,    & ! average u-momentum flux (surface)             ( N/m2 )
    avmfl_s    ,    & ! average v-momentum flux (surface)             ( N/m2 )
    ashfl_s    ,    & ! average sensible heat flux (surface)          ( W/m2 )
    alhfl_s           ! average latent heat flux (surface)            ( W/m2 )

!------------------------------------------------------------------------------

USE data_fields     , ONLY :   &

! 7. fields for model output and diagnostics                          (unit)
! ------------------------------------------

    t_2m       ,    & ! temperature in 2m                             (  K   )
    t_2m_av    ,    & ! time mean temperature in 2m                   (  K   )
    qv_2m      ,    & ! specific water vapor content in 2m            (kg/kg )
    td_2m      ,    & ! dew-point in 2m                               (  K   )
    td_2m_av   ,    & ! time mean dew-point in 2m                     (  K   )
    rh_2m      ,    & ! relative humidity in 2m                       (  %   )
    u_10m      ,    & ! zonal wind in 10m                             ( m/s  )
    u_10m_av   ,    & ! time mean zonal wind in 10m                   ( m/s  )
    v_10m      ,    & ! meridional wind in 10m                        ( m/s  )
    v_10m_av   ,    & ! time mean meridional wind in 10m              ( m/s  )
    tmin_2m    ,    & ! minimum temperature in 2m                     (  K   )
    tmax_2m    ,    & ! maximum temperature in 2m                     (  K   )
    vmax_10m   ,    & ! maximal windspeed in 10m                      ( m/s  )
    vabsmx_10m ,    & ! maximal windspeed in 10m without gust         ( m/s  )
    vgust_dyn  ,    & ! maximal dynamical wind gust in 10m            ( m/s )
    vgust_con  ,    & ! maximal convective wind gust in 10m           ( m/s )
    asob_s     ,    & ! average solar radiation budget (surface)      ( W/m2 )
    athb_s     ,    & ! average thermal radiation budget (surface)    ( W/m2 )
    apab_s     ,    & ! average photosynthetic active radiation (sfc) ( W/m2 )
    asob_t     ,    & ! average solar radiation budget (model top)    ( W/m2 )
    athb_t     ,    & ! average thermal radiation budget (model top)  ( W/m2 )
     sod_t     ,    & ! solar downward radiation at top of atmosphere (     )
    asod_t     ,    & ! averaged solar downward radiation at top      (     )
    dursun     ,    & ! sunshine duration                             (  s   )
    dursun_m   ,    & ! maximum possible sunshine duration            (  s  )
    dursun_r          ! relative sunshine duration                    (  s  )

USE data_fields     , ONLY :   &
    rain_gsp   ,    & ! amount of rain from grid-scale precip. (sum)  (kg/m2 )
    snow_gsp   ,    & ! amount of snow from grid-scale precip. (sum)  (kg/m2 )
    grau_gsp   ,    & ! amount of graupel from grid-scale prec. (sum) (kg/m2 )
    hail_gsp   ,    & ! amount of hail from grid-scale prec.    (sum) (kg/m2 )
    rain_con   ,    & ! amount of rain from convective precip. (sum)  (kg/m2 )
    snow_con   ,    & ! amount of snow from convective precip. (sum)  (kg/m2 )
    runoff_s   ,    & ! surface water runoff; sum over forecast       (kg/m2 )
    runoff_g   ,    & ! soil water runoff; sum over forecast          (kg/m2 )
    rstom      ,    & ! stomata resistance                            ( s/m )
    lhfl_bs    ,    & ! latent heat flux from bare soil evap.         ( W/m2)
    lhfl_pl    ,    & ! latent heat flux from plants                  ( W/m2)
    alhfl_bs   ,    & ! average latent heat flux from bare soil evap. ( W/m2)
    alhfl_pl   ,    & ! average latent heat flux from plants          ( W/m2)
    tdiv_hum   ,    & ! vertical integral divergence of humidity      (kg/m2 )
    aevap_s           ! accumulated surface moisture flux             (kg/m2 )

!------------------------------------------------------------------------------

USE data_fields     , ONLY :   &

! 8. fields for the boundary values                                   (unit)
! ---------------------------------

    u_bd          , & ! boundary field for u                          ( m/s  )
    v_bd          , & ! boundary field for v                          ( m/s  )
    w_bd          , & ! boundary field for w                          ( m/s  )
    t_bd          , & ! boundary field for t                          (  k   )
    pp_bd         , & ! boundary field for pp                         (  pa  )
    qv_s_bd       , & ! boundary field for qv_s                       (kg/kg )
    t_snow_bd     , & ! boundary field for t_snow                     (  k   )
    t_s_bd        , & ! boundary field for t_s                        (  k   )
    t_m_bd        , & ! boundary field for t_m                        (  k   )
    w_snow_bd     , & ! boundary field for w_snow                     (m H2O )
    w_g1_bd       , & ! boundary field for w_g1                       (m H2O )
    w_g2_bd       , & ! boundary field for w_g2                       (m H2O )
    w_g3_bd       , & ! boundary field for w_g3                       (m H2O )
    hmo3_bd       , & ! boundary field for hmo3                       (m    )
    vio3_bd       , & ! boundary field for vio3                       (pa O3)
    w_cl_bd       , & ! boundary field for w_cl                       (m H2O)
    t_cl_bd       , & ! boundary field for t_cl                       (  K  )
    lai_bd        , & ! boundary field for lai                        ( --  )
    rootdp_bd     , & ! boundary field for rootdp                     (m    )
    plcov_bd      , & ! boundary field for plcov                      ( --  )
    aer_su_bd     , & ! boundary field for aer_su                     (0 - 1)
    aer_bc_bd     , & ! boundary field for aer_bc                     (0 - 1)
    aer_du_bd     , & ! boundary field for aer_du                     (0 - 1)
    aer_or_bd     , & ! boundary field for aer_or                     (0 - 1)
    aer_ss_bd     , & ! boundary field for aer_ss                     (0 - 1)

! 9. fields for the synthetic satellite images
! --------------------------------------------

    synme7        , & ! Meteosat 7
    synmsg        , & ! Meteosat Second Generation

! 10. analysis increment fields
! -----------------------------

    ff_anai       , & ! wind velocity                                 ( m/s )
    dd_anai       , & ! wind direction                                ( rad )
    t_anai        , & ! temperature                                   (  k  )
    p_anai        , & ! deviation from the reference pressure         ( Pa  )
    qv_anai       , & ! specific water vapor content                  (kg/kg)
    qc_anai           ! specific cloud water content (via saturation adjustm)

! end of data_fields

!------------------------------------------------------------------------------

USE data_io,        ONLY :   &
    lbdclim,      & ! boundary data in climate model
    lbdsst          ! T_S boundary data are used only over the sea
                    ! (SST is not maintained constant during the integration)

!------------------------------------------------------------------------------

USE data_parallel,  ONLY :   &
    my_cart_id,   & !
    nprocx          ! number of processors in x-direction

!------------------------------------------------------------------------------

#ifdef NUDGING
USE data_lheat_nudge,  ONLY :   &
    llhn         ,& ! on/off switch for latent heat nudging (lhn)
    llhnverif    ,& ! on/off switch for latent heat nudging (lhn)
    lhn_qrs      ,& ! use integrated precipitation flux as reference
    tt_lheat     ,& ! profile of t-increments due to latent heating   ( K/s  )
    qrsflux         ! total precipitation flux

!------------------------------------------------------------------------------

USE data_lhn_diag, ONLY  :  &
!   fields for lh-diagnostics
    tinc_lhn_o     ,& ! temperature increments due to lhn               ( K/s  )
    tt_lheat_o     ,& ! array for cumulated latent heating (grid scale + conv)( K )
    ttm_cv_o          ! array for test output of diverse 2D fields
#endif

!------------------------------------------------------------------------------

#ifdef MESSY
 USE messy_main_channel_bi,  ONLY: messy_COSMO_create_channel &
                                 , messy_COSMO_dealloc_meteofields
#endif

!------------------------------------------------------------------------------

IMPLICIT NONE

!------------------------------------------------------------------------------

!==============================================================================
! Module Procedures in src_allocation
!==============================================================================

CONTAINS

!==============================================================================
!+ organizing the allocation and deallocation of the meteorological fields
!------------------------------------------------------------------------------

SUBROUTINE organize_allocation ( yaction, ierrstat )
!------------------------------------------------------------------------------
!
! Description:
!   This external subroutine organizes the allocations and deallocations
!   of the global meteorological fields for all grids of the model run.
!
!   The actions are:
!     - 'default' :  allocation of the default fields before first data input
!     - 'canopy'  :  allocation of the extra fields for turbulence and canopy
!     - 'dealloc' :  deallocation of all allocated fields
!
! Method:
!   IF-construction for checking yaction
!
!==============================================================================

! Parameter list:
  CHARACTER (LEN= *),       INTENT(IN)            ::                        &
    yaction       ! action to be performed

  INTEGER (KIND=iintegers), INTENT(OUT)           ::                        &
    ierrstat      ! for global error code

!==============================================================================

  ierrstat = 0

  IF ( yaction == 'default' .OR. yaction == 'canopy' ) THEN

#ifndef MESSY
    CALL alloc_meteofields( yaction, ierrstat )
#else
    CALL messy_COSMO_create_channel( yaction, ierrstat )
#endif

  ELSEIF ( yaction == 'dealloc' ) THEN

#ifndef MESSY
    CALL dealloc_meteofields( ierrstat )
#else
    CALL messy_COSMO_dealloc_meteofields
#endif

  ENDIF

END SUBROUTINE organize_allocation

!==============================================================================
!+ Allocates the meteorological fields and initializes them with 0
!------------------------------------------------------------------------------
#ifndef MESSY
SUBROUTINE alloc_meteofields (yaction, ist)

!------------------------------------------------------------------------------
!
! Description:
!   This routine allocates space for the meteorological fields that are used 
!   in the long term storage and initializes them with 0.
!
! Method:
!   All ALLOCATABLE fields are allocated with the ALLOCATE statement.
!   The status of the allocation is checked. In case of an error, the
!   error variable is set: ist = 1
!
!==============================================================================
!
! Subroutine / Function arguments
INTEGER (KIND=iintegers), INTENT(OUT)   ::       &
  ist                ! for global error-code

! Parameter list:
CHARACTER (LEN= *),       INTENT(IN)    ::       &
  yaction            ! action to be performed

! Local variables
INTEGER (KIND=iintegers)                ::       &
  izl      ,       & ! for local error-code
  izdebug  ,       & ! for local debug output
  nztlev             ! number of time-levels for prognostic variables

!------------------------------------------------------------------------------
!- End of header -
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!- Begin SUBROUTINE alloc_fields
!------------------------------------------------------------------------------

  ist = 0
  izl = 0

  IF (lprintdeb_all) THEN
    izdebug = idbg_level
  ELSE
    IF (my_cart_id == 0) THEN
      izdebug = idbg_level
    ELSE
      izdebug = 0
    ENDIF
  ENDIF

! The default number of time-levels is 3 for the Klemp-Wilhelmson/Leapfrog
! scheme. If the two time level RK-schem is used (l2tls = .TRUE), nztlev
! is set to 2 and less space is allocated for the prognostic variables.
  IF ( l2tls ) THEN
    nztlev = 2
  ELSE
    nztlev = 3
  ENDIF

!------------------------------------------------------------------------------
! Section 1:  default allocation
!------------------------------------------------------------------------------

IF ( yaction == 'default' ) THEN

! constant fields for the reference atmosphere and the vertical coordinate
! ------------------------------------------------------------------------

  ALLOCATE ( hhl_prof(0:ke+1)   , STAT=izl ) ; hhl_prof=0.0_wp  ; ist = ist + izl
  ALLOCATE ( rho0   (ie,je,ke)  , STAT=izl ) ; rho0   = 0.0_wp  ; ist = ist + izl
  ALLOCATE ( p0     (ie,je,ke)  , STAT=izl ) ; p0     = 0.0_wp  ; ist = ist + izl
  ALLOCATE ( p0hl   (ie,je,ke+1), STAT=izl ) ; p0hl   = 0.0_wp  ; ist = ist + izl
  ALLOCATE ( dp0    (ie,je,ke)  , STAT=izl ) ; dp0    = 0.0_wp  ; ist = ist + izl
  ALLOCATE ( hhl    (ie,je,ke+1), STAT=izl ) ; hhl    = 0.0_wp  ; ist = ist + izl
  ALLOCATE ( dt0dz  (ie,je,ke)  , STAT=izl ) ; dt0dz  = 0.0_wp  ; ist = ist + izl
  ALLOCATE ( t0     (ie,je,ke)  , STAT=izl ) ; t0     = 0.0_wp  ; ist = ist + izl
  ALLOCATE ( t0hl   (ie,je,ke+1), STAT=izl ) ; t0hl   = 0.0_wp  ; ist = ist + izl

! external parameter fields
! -------------------------

  ALLOCATE ( hsurf  (ie,je)  , STAT=izl ) ; hsurf  = 0.0_wp   ; ist = ist + izl
  IF (lsso) THEN
    ALLOCATE (sso_stdh(ie,je) ,STAT=izl ) ; sso_stdh = 0.0_wp ; ist = ist + izl
    ALLOCATE (sso_gamma(ie,je),STAT=izl ) ; sso_gamma= 0.0_wp ; ist = ist + izl
    ALLOCATE (sso_theta(ie,je),STAT=izl ) ; sso_theta= 0.0_wp ; ist = ist + izl
    ALLOCATE (sso_sigma(ie,je),STAT=izl ) ; sso_sigma= 0.0_wp ; ist = ist + izl
  ENDIF

  IF (itype_aerosol == 2) THEN
    ALLOCATE (aer_su (ie,je),  STAT=izl ) ; aer_su   = 0.0_wp ; ist = ist + izl
    ALLOCATE (aer_du (ie,je),  STAT=izl ) ; aer_du   = 0.0_wp ; ist = ist + izl
    ALLOCATE (aer_or (ie,je),  STAT=izl ) ; aer_or   = 0.0_wp ; ist = ist + izl
    ALLOCATE (aer_bc (ie,je),  STAT=izl ) ; aer_bc   = 0.0_wp ; ist = ist + izl
    ALLOCATE (aer_ss (ie,je),  STAT=izl ) ; aer_ss   = 0.0_wp ; ist = ist + izl
  ENDIF

  IF (lemiss) THEN
    ALLOCATE (emis_rad(ie,je), STAT=izl ) ; emis_rad = 0.0_wp ; ist = ist + izl
  ENDIF

  IF (lstomata) THEN
    ALLOCATE (rsmin2d (ie,je), STAT=izl ) ; rsmin2d  = 0.0_wp ; ist = ist + izl
  ENDIF

  ALLOCATE ( gz0    (ie,je)  , STAT=izl ) ; gz0    = 0.0_wp   ; ist = ist + izl
  ALLOCATE ( fr_land(ie,je)  , STAT=izl ) ; fr_land= 0.0_wp   ; ist = ist + izl
  ALLOCATE ( soiltyp(ie,je)  , STAT=izl ) ; soiltyp= 0.0_wp   ; ist = ist + izl
  ALLOCATE ( vio3   (ie,je)  , STAT=izl ) ; vio3   = 0.0_wp   ; ist = ist + izl
  ALLOCATE ( hmo3   (ie,je)  , STAT=izl ) ; hmo3   = 0.0_wp   ; ist = ist + izl
  ALLOCATE ( rlat   (ie,je)  , STAT=izl ) ; rlat   = 0.0_wp   ; ist = ist + izl
  ALLOCATE ( rlon   (ie,je)  , STAT=izl ) ; rlon   = 0.0_wp   ; ist = ist + izl
  IF ((nprocx > 1) .AND. (lreproduce)) THEN
    ALLOCATE (rlattot(ie_tot,je),STAT=izl); rlattot= 0.0_wp   ; ist = ist + izl
    ALLOCATE (rlontot(ie_tot,je),STAT=izl); rlontot= 0.0_wp   ; ist = ist + izl
  ENDIF
  ALLOCATE ( fc     (ie,je)  , STAT=izl ) ; fc     = 0.0_wp   ; ist = ist + izl
  IF (lcori_deep) THEN
    ALLOCATE (fccos (ie,je)  , STAT=izl ) ; fccos  = 0.0_wp   ; ist = ist + izl
  ENDIF
  ALLOCATE ( rmy    (ie,je,3), STAT=izl ) ; rmy    = 0.0_wp   ; ist = ist + izl
  ALLOCATE ( rmyq   (ie,je)  , STAT=izl ) ; rmyq   = 0.0_wp   ; ist = ist + izl
  ALLOCATE ( hd_mask_dcoeff_p(ie,je,ke), STAT=izl ) ; hd_mask_dcoeff_p = 0.0_wp ; ist = ist + izl
  ALLOCATE ( hd_mask_dcoeff_t(ie,je,ke), STAT=izl ) ; hd_mask_dcoeff_t = 0.0_wp ; ist = ist + izl
  ALLOCATE ( hd_mask_dcoeff_q(ie,je,ke), STAT=izl ) ; hd_mask_dcoeff_q = 0.0_wp ; ist = ist + izl
  ALLOCATE ( hd_mask_dcoeff_u(ie,je,ke), STAT=izl ) ; hd_mask_dcoeff_u = 0.0_wp ; ist = ist + izl
  ALLOCATE ( ofa_hdx(ie,je,ke),STAT=izl ) ; ofa_hdx= 1.0_wp   ; ist = ist + izl
  ALLOCATE ( ofa_hdy(ie,je,ke),STAT=izl ) ; ofa_hdy= 1.0_wp   ; ist = ist + izl
  ALLOCATE ( hd_mask(ie,je,ke),STAT=izl ) ; hd_mask= 1.0_wp   ; ist = ist + izl
  ALLOCATE ( least_lbdz(ie,je),STAT=izl ) ; least_lbdz =.FALSE.; ist=ist+izl
  ALLOCATE ( lwest_lbdz(ie,je),STAT=izl ) ; lwest_lbdz =.FALSE.; ist=ist+izl
  ALLOCATE ( lnorth_lbdz(ie,je),STAT=izl) ; lnorth_lbdz=.FALSE.; ist=ist+izl
  ALLOCATE ( lsouth_lbdz(ie,je),STAT=izl) ; lsouth_lbdz=.FALSE.; ist=ist+izl
  ALLOCATE ( crlat  (je,2)   , STAT=izl ) ; crlat  = 0.0_wp   ; ist = ist + izl
  ALLOCATE ( acrlat (je,2)   , STAT=izl ) ; acrlat = 0.0_wp   ; ist = ist + izl
  ALLOCATE ( tgrlat (je,2)   , STAT=izl ) ; tgrlat = 0.0_wp   ; ist = ist + izl
  ALLOCATE ( aerlan (ie,je)  , STAT=izl ) ; aerlan = 0.0_wp   ; ist = ist + izl
  ALLOCATE ( aerurb (ie,je)  , STAT=izl ) ; aerurb = 0.0_wp   ; ist = ist + izl
  ALLOCATE ( aerdes (ie,je)  , STAT=izl ) ; aerdes = 0.0_wp   ; ist = ist + izl
  ALLOCATE ( aersea (ie,je)  , STAT=izl ) ; aersea = 0.0_wp   ; ist = ist + izl
  ALLOCATE ( plcov  (ie,je)  , STAT=izl ) ; plcov  = 0.0_wp   ; ist = ist + izl
  ALLOCATE ( lai    (ie,je)  , STAT=izl ) ; lai    = 0.0_wp   ; ist = ist + izl
  ALLOCATE ( tai    (ie,je)  , STAT=izl ) ; tai    = 0.0_wp   ; ist = ist + izl
  ALLOCATE ( sai    (ie,je)  , STAT=izl ) ; sai    = 0.0_wp   ; ist = ist + izl
  ALLOCATE ( eai    (ie,je)  , STAT=izl ) ; eai    = 0.0_wp   ; ist = ist + izl
  ALLOCATE ( rootdp (ie,je)  , STAT=izl ) ; rootdp = 0.0_wp   ; ist = ist + izl
  ALLOCATE ( llandmask(ie,je), STAT=izl ) ; llandmask = .TRUE.; ist = ist + izl

  IF (lseaice) THEN
    ALLOCATE (lseamask(ie,je), STAT=izl ) ; lseamask = .FALSE.; ist = ist + izl
  ENDIF

  ALLOCATE ( h_can  (ie,je)  , STAT=izl ) ; h_can  = 0.0_wp   ; ist = ist + izl
  ALLOCATE ( d_pat  (ie,je)  , STAT=izl ) ; d_pat  = 0.0_wp   ; ist = ist + izl

  IF (lforest) THEN
    ALLOCATE (for_e (ie,je)  , STAT=izl ) ; for_e  = 0.0_wp   ; ist = ist + izl
    ALLOCATE (for_d (ie,je)  , STAT=izl ) ; for_d  = 0.0_wp   ; ist = ist + izl
  ENDIF

  ! Because these fields are used at other parts of the code, allocate anyhow
  ALLOCATE (depth_lk(ie,je), STAT=izl ) ; depth_lk = 0.0_wp ; ist = ist + izl
  ALLOCATE (fr_lake (ie,je), STAT=izl ) ; fr_lake  = 0.0_wp ; ist = ist + izl
  IF (llake) THEN
    ALLOCATE (fetch_lk(ie,je), STAT=izl ) ; fetch_lk = 0.0_wp ; ist = ist + izl
    ALLOCATE (dp_bs_lk(ie,je), STAT=izl ) ; dp_bs_lk = 0.0_wp ; ist = ist + izl
    ALLOCATE (t_bs_lk (ie,je), STAT=izl ) ; t_bs_lk  = 0.0_wp ; ist = ist + izl
    ALLOCATE (gamso_lk(ie,je), STAT=izl ) ; gamso_lk = 0.0_wp ; ist = ist + izl
  ENDIF

  IF     (itype_albedo == 2) THEN
    ALLOCATE (alb_dry (ie,je), STAT=izl ) ; alb_dry  = 0.0_wp ; ist = ist + izl
    ALLOCATE (alb_sat (ie,je), STAT=izl ) ; alb_sat  = 0.0_wp ; ist = ist + izl
  ELSEIF (itype_albedo == 3) THEN
    ALLOCATE (alb_dif (ie,je), STAT=izl ) ; alb_dif  = 0.0_wp ; ist = ist + izl
  ENDIF

  IF (izdebug > 10) THEN
    PRINT *, '   ALLOCATED constant and external parameter fields:  ', ist
  ENDIF

! prognostic variables
! --------------------

  ALLOCATE ( u  (ie,je,ke  ,nztlev), STAT=izl ) ; u  = 0.0_wp  ; ist = ist + izl
  ALLOCATE ( v  (ie,je,ke  ,nztlev), STAT=izl ) ; v  = 0.0_wp  ; ist = ist + izl
  ALLOCATE ( w  (ie,je,ke+1,nztlev), STAT=izl ) ; w  = 0.0_wp  ; ist = ist + izl
  ALLOCATE ( t  (ie,je,ke  ,nztlev), STAT=izl ) ; t  = 0.0_wp  ; ist = ist + izl
  ALLOCATE ( pp (ie,je,ke  ,nztlev), STAT=izl ) ; pp = 0.0_wp  ; ist = ist + izl

! tendency fields for the prognostic variables
! --------------------------------------------

  ALLOCATE ( utens  (ie,je,ke)  , STAT=izl ) ; utens  = 0.0_wp ; ist = ist + izl
  ALLOCATE ( vtens  (ie,je,ke)  , STAT=izl ) ; vtens  = 0.0_wp ; ist = ist + izl
  ALLOCATE ( wtens  (ie,je,ke+1), STAT=izl ) ; wtens  = 0.0_wp ; ist = ist + izl
  ALLOCATE ( ttens  (ie,je,ke)  , STAT=izl ) ; ttens  = 0.0_wp ; ist = ist + izl
  ALLOCATE ( pptens (ie,je,ke)  , STAT=izl ) ; pptens = 0.0_wp ; ist = ist + izl

  IF (izdebug > 10) THEN
    PRINT *, '   ALLOCATED prognostic variables and tendencies:     ', ist
  ENDIF

! fields for surface values and soil model variables
! --------------------------------------------------

  ALLOCATE ( ps  (ie,je,nztlev) , STAT=izl ) ; ps     = 0.0_wp  ; ist = ist + izl
  ALLOCATE ( t_snow(ie,je,nztlev),STAT=izl ) ; t_snow = 0.0_wp  ; ist = ist + izl
  ALLOCATE ( t_s (ie,je,nztlev) , STAT=izl ) ; t_s    = 0.0_wp  ; ist = ist + izl
  ALLOCATE ( t_g (ie,je,nztlev) , STAT=izl ) ; t_g    = 0.0_wp  ; ist = ist + izl
  ALLOCATE ( tg_radstep(ie,je)  , STAT=izl ) ; tg_radstep=0.0_wp; ist = ist + izl
  ALLOCATE ( qv_s(ie,je,nztlev) , STAT=izl ) ; qv_s   = 0.0_wp  ; ist = ist + izl
  ALLOCATE ( w_snow(ie,je,nztlev),STAT=izl ) ; w_snow = 0.0_wp  ; ist = ist + izl
  ALLOCATE ( w_i (ie,je,nztlev) , STAT=izl ) ; w_i    = 0.0_wp  ; ist = ist + izl
  ALLOCATE ( h_snow(ie,je,nztlev),STAT=izl ) ; h_snow = 0.0_wp  ; ist = ist + izl
  ALLOCATE ( snow_melt(ie,je),    STAT=izl ) ; snow_melt=0.0_wp ; ist = ist + izl
  ALLOCATE ( fr_snow  (ie,je),    STAT=izl ) ; fr_snow= 0.0_wp  ; ist = ist + izl
  ALLOCATE ( fr_wi    (ie,je),    STAT=izl ) ; fr_wi  = 0.0_wp  ; ist = ist + izl
  ALLOCATE ( ustar_fv (ie,je),    STAT=izl ) ; ustar_fv=0.0_wp  ; ist = ist + izl

  IF (lmulti_layer) THEN
    ALLOCATE (t_so      (ie,je,0:ke_soil+1,nztlev), STAT=izl);
                                                 t_so = 0.0_wp  ; ist = ist + izl
    ALLOCATE (w_so      (ie,je,  ke_soil+1,nztlev), STAT=izl);
                                                 w_so = 0.0_wp  ; ist = ist + izl
    ALLOCATE (w_so_ice  (ie,je,  ke_soil+1,nztlev), STAT=izl);
                                              w_so_ice= 0.0_wp  ; ist = ist + izl
    ALLOCATE (rho_snow  (ie,je,            nztlev), STAT=izl);
                                              rho_snow = 0.0_wp ; ist = ist + izl
    IF(lmulti_snow) THEN
      ALLOCATE (t_snow_mult(ie,je,0:ke_snow,nztlev), STAT=izl );
                                           t_snow_mult = 0.0_wp  ; ist = ist + izl
      ALLOCATE (dzh_snow_mult(ie,je,  ke_snow,nztlev), STAT=izl );
                                        dzh_snow_mult  = 0.0_wp  ; ist = ist + izl
      ALLOCATE (w_snow_mult  (ie,je,  ke_snow,nztlev), STAT=izl );
                                           w_snow_mult = 0.0_wp  ; ist = ist + izl
      ALLOCATE (wliq_snow    (ie,je,  ke_snow,nztlev), STAT=izl );
                                             wliq_snow = 0.0_wp  ; ist = ist + izl
      ALLOCATE (rho_snow_mult(ie,je,ke_snow,nztlev), STAT=izl);
                                         rho_snow_mult = 0.0_wp ; ist = ist + izl
    ENDIF

    ALLOCATE (freshsnow (ie,je),                    STAT=izl);
                                             freshsnow = 1.0_wp ; ist = ist + izl
    ! in case of climate runs, also allocate climatological deep soil variables
    IF (lbdclim) THEN
      ALLOCATE (t_cl(ie,je     ), STAT=izl ) ; t_cl   = 0.0_wp  ; ist = ist + izl
      lalloc_t_cl = .TRUE.
      ALLOCATE (w_cl(ie,je     ), STAT=izl ) ; w_cl   = 0.0_wp  ; ist = ist + izl
    ENDIF
  ELSE
    ALLOCATE (t_m (ie,je,nztlev), STAT=izl ) ; t_m    = 0.0_wp  ; ist = ist + izl
    ALLOCATE (t_cl(ie,je       ), STAT=izl ) ; t_cl   = 0.0_wp  ; ist = ist + izl
    lalloc_t_cl = .TRUE.
    ALLOCATE (w_g1(ie,je,nztlev), STAT=izl ) ; w_g1   = 0.0_wp  ; ist = ist + izl
    ALLOCATE (w_g2(ie,je,nztlev), STAT=izl ) ; w_g2   = 0.0_wp  ; ist = ist + izl
    IF ( (nlgw_ini == 3) .OR. (nlgw == 3) ) THEN
      ALLOCATE (w_g3(ie,je,nztlev), STAT=izl); w_g3   = 0.0_wp  ; ist = ist + izl
      lalloc_w_g3 = .TRUE.
    ENDIF
    ALLOCATE (w_cl(ie,je       ), STAT=izl ) ; w_cl   = 0.0_wp  ; ist = ist + izl
  ENDIF

  IF (lseaice .OR. llake) THEN
    ALLOCATE (t_ice   (ie,je,nztlev), STAT=izl ) ; t_ice    = 0.0_wp ; ist = ist + izl
    ALLOCATE (h_ice   (ie,je,nztlev), STAT=izl ) ; h_ice    = 0.0_wp ; ist = ist + izl
    lalloc_h_ice = .TRUE.
  ENDIF

  IF (llake) THEN
    ALLOCATE (t_mnw_lk(ie,je,nztlev), STAT=izl ) ; t_mnw_lk = 0.0_wp ; ist = ist + izl
    ALLOCATE (t_wml_lk(ie,je,nztlev), STAT=izl ) ; t_wml_lk = 0.0_wp ; ist = ist + izl
    ALLOCATE (t_bot_lk(ie,je,nztlev), STAT=izl ) ; t_bot_lk = 0.0_wp ; ist = ist + izl
    ALLOCATE (t_b1_lk (ie,je,nztlev), STAT=izl ) ; t_b1_lk  = 0.0_wp ; ist = ist + izl
    ALLOCATE (c_t_lk  (ie,je,nztlev), STAT=izl ) ; c_t_lk   = 0.0_wp ; ist = ist + izl
    ALLOCATE (h_ml_lk (ie,je,nztlev), STAT=izl ) ; h_ml_lk  = 0.0_wp ; ist = ist + izl
    ALLOCATE (h_b1_lk (ie,je,nztlev), STAT=izl ) ; h_b1_lk  = 0.0_wp ; ist = ist + izl
  ENDIF

  IF (izdebug > 10) THEN
    PRINT *, '   ALLOCATED soil fields:                             ', ist
  ENDIF

! fields that are computed in the parametrization and dynamics
! ------------------------------------------------------------

  IF ( ldiabf_lh ) THEN
    ALLOCATE ( tinc_lh(ie,je,ke), STAT=izl ) ; tinc_lh = 0.0_wp  ; ist = ist + izl
  ENDIF
  ALLOCATE ( rho    (ie,je,ke)  , STAT=izl ) ; rho     = 0.0_wp  ; ist = ist + izl
  ALLOCATE ( tkvm   (ie,je,ke1) , STAT=izl ) ; tkvm    = 0.0_wp  ; ist = ist + izl
  ALLOCATE ( tkvh   (ie,je,ke1) , STAT=izl ) ; tkvh    = 0.0_wp  ; ist = ist + izl
!US  IF ( l3dturb ) THEN
    ALLOCATE ( tkhm (ie,je,ke1) , STAT=izl ) ; tkhm    = 0.0_wp  ; ist = ist + izl
    ALLOCATE ( tkhh (ie,je,ke1) , STAT=izl ) ; tkhh    = 0.0_wp  ; ist = ist + izl
!US  ENDIF
  ALLOCATE ( tcm    (ie,je)     , STAT=izl ) ; tcm     = 0.0_wp  ; ist = ist + izl
  ALLOCATE ( tch    (ie,je)     , STAT=izl ) ; tch     = 0.0_wp  ; ist = ist + izl
  ALLOCATE ( tfm    (ie,je)     , STAT=izl ) ; tfm     = 1.0_wp  ; ist = ist + izl
  ALLOCATE ( tfh    (ie,je)     , STAT=izl ) ; tfh     = 1.0_wp  ; ist = ist + izl
  ALLOCATE ( tfv    (ie,je)     , STAT=izl ) ; tfv     = 1.0_wp  ; ist = ist + izl
  ALLOCATE ( sohr   (ie,je,ke)  , STAT=izl ) ; sohr    = 0.0_wp  ; ist = ist + izl
  ALLOCATE ( sotr   (ie,je,ke1) , STAT=izl ) ; sotr    = 0.0_wp  ; ist = ist + izl
  ALLOCATE ( sotr_par(ie,je)    , STAT=izl ) ; sotr_par= 0.0_wp  ; ist = ist + izl
  ALLOCATE ( thhr   (ie,je,ke)  , STAT=izl ) ; thhr    = 0.0_wp  ; ist = ist + izl
  ALLOCATE ( clc_sgs(ie,je,ke)  , STAT=izl ) ; clc_sgs = 0.0_wp  ; ist = ist + izl
#ifdef COUP_OAS_COS
  ALLOCATE ( tcw    (ie,je)     , STAT=izl ) ; tcw     = 0.0_wp  ; ist = ist + izl
  ALLOCATE ( trad_clm(ie,je)    , STAT=izl ) ; trad_clm= 0.0_wp  ; ist = ist + izl
  ALLOCATE ( alb_rad(ie,je,2)   , STAT=izl ) ; alb_rad = 0.0_wp  ; ist = ist + izl
#else
  ALLOCATE ( alb_rad(ie,je)     , STAT=izl ) ; alb_rad = 0.0_wp  ; ist = ist + izl
#endif
  IF ( nradcoarse > 1 ) THEN
    ALLOCATE ( alb_rad_coarse(ie,je) , STAT=izl ) ; alb_rad_coarse = 0.0_wp  ; ist = ist + izl
  ENDIF
  ALLOCATE ( sobs   (ie,je)     , STAT=izl ) ; sobs    = 0.0_wp  ; ist = ist + izl
  ALLOCATE ( thbs   (ie,je)     , STAT=izl ) ; thbs    = 0.0_wp  ; ist = ist + izl
  ALLOCATE ( pabs   (ie,je)     , STAT=izl ) ; pabs    = 0.0_wp  ; ist = ist + izl
  ALLOCATE ( sobt   (ie,je)     , STAT=izl ) ; sobt    = 0.0_wp  ; ist = ist + izl
  ALLOCATE ( thbt   (ie,je)     , STAT=izl ) ; thbt    = 0.0_wp  ; ist = ist + izl
  ALLOCATE ( clch   (ie,je)     , STAT=izl ) ; clch    = 0.0_wp  ; ist = ist + izl
  ALLOCATE ( clcm   (ie,je)     , STAT=izl ) ; clcm    = 0.0_wp  ; ist = ist + izl
  ALLOCATE ( clcl   (ie,je)     , STAT=izl ) ; clcl    = 0.0_wp  ; ist = ist + izl
  ALLOCATE ( clct   (ie,je)     , STAT=izl ) ; clct    = 0.0_wp  ; ist = ist + izl
  ALLOCATE ( sodwddm(ie,je)     , STAT=izl ) ; sodwddm = 0.0_wp  ; ist = ist + izl
  ALLOCATE ( qc_rad (ie,je,ke)  , STAT=izl ) ; qc_rad  = 0.0_wp  ; ist = ist + izl
  ALLOCATE ( qi_rad (ie,je,ke)  , STAT=izl ) ; qi_rad  = 0.0_wp  ; ist = ist + izl
  ALLOCATE ( sun_el (ie,je)     , STAT=izl ) ; sun_el  = 0.0_wp  ; ist = ist + izl
  ALLOCATE ( sun_azi(ie,je)     , STAT=izl ) ; sun_azi = 0.0_wp  ; ist = ist + izl

  IF (lradtopo) THEN
    ALLOCATE (skyview  (ie,je)   , STAT=izl ); skyview = 0.0_wp  ; ist = ist + izl
    ALLOCATE (slo_asp  (ie,je)   , STAT=izl ); slo_asp = 0.0_wp  ; ist = ist + izl
    ALLOCATE (slo_ang  (ie,je)   , STAT=izl ); slo_ang = 0.0_wp  ; ist = ist + izl
    ALLOCATE (horizon  (ie,je,nhori), STAT=izl ); horizon = 0.0_wp  ; ist = ist + izl
  ENDIF
  ! the other fields are calculated in any case in the radiation scheme
  ALLOCATE (swdir_s  (ie,je)    , STAT=izl ) ; swdir_s = 0.0_wp  ; ist = ist + izl
  ALLOCATE (swdifd_s (ie,je)    , STAT=izl ) ; swdifd_s= 0.0_wp  ; ist = ist + izl
  ALLOCATE (swdifu_s (ie,je)    , STAT=izl ) ; swdifu_s= 0.0_wp  ; ist = ist + izl
  ALLOCATE (swtrdir_s(ie,je)    , STAT=izl ) ; swtrdir_s = 0.0_wp; ist = ist + izl
  ALLOCATE (swtrdifd_s(ie,je)   , STAT=izl ) ; swtrdifd_s= 0.0_wp; ist = ist + izl
  ALLOCATE (swtrdifu_s(ie,je)   , STAT=izl ) ; swtrdifu_s= 0.0_wp; ist = ist + izl
  ALLOCATE (lwd_s    (ie,je)    , STAT=izl ) ; lwd_s   = 0.0_wp  ; ist = ist + izl
  ALLOCATE (lwu_s    (ie,je)    , STAT=izl ) ; lwu_s   = 0.0_wp  ; ist = ist + izl
  ALLOCATE (aswdir_s (ie,je)    , STAT=izl ) ; aswdir_s= 0.0_wp  ; ist = ist + izl
  ALLOCATE (aswdifd_s(ie,je)    , STAT=izl ) ; aswdifd_s=0.0_wp  ; ist = ist + izl
  ALLOCATE (aswdifu_s(ie,je)    , STAT=izl ) ; aswdifu_s=0.0_wp  ; ist = ist + izl
  ALLOCATE (alwd_s   (ie,je)    , STAT=izl ) ; alwd_s  = 0.0_wp  ; ist = ist + izl
  ALLOCATE (alwu_s   (ie,je)    , STAT=izl ) ; alwu_s  = 0.0_wp  ; ist = ist + izl
  ALLOCATE (swdir_cor(ie,je)    , STAT=izl ) ; swdir_cor=0.0_wp  ; ist = ist + izl

#ifdef TWOMOM_SB
  ALLOCATE (reffc_out (ie,je,ke), STAT=izl ) ; reffc_out = 10.0E-6_wp; ist = ist + izl
  ALLOCATE (reffi_out (ie,je,ke), STAT=izl ) ; reffi_out = 20.0E-6_wp; ist = ist + izl
  ALLOCATE (odepthw_so(ie,je,ke), STAT=izl ) ; odepthw_so = 0.0_wp   ; ist = ist + izl
  ALLOCATE (odepthi_so(ie,je,ke), STAT=izl ) ; odepthi_so = 0.0_wp   ; ist = ist + izl
  ALLOCATE (odepthw_th(ie,je,ke), STAT=izl ) ; odepthw_th = 0.0_wp   ; ist = ist + izl
  ALLOCATE (odepthi_th(ie,je,ke), STAT=izl ) ; odepthi_th = 0.0_wp   ; ist = ist + izl
#endif

  ALLOCATE (  clc_con(ie,je,ke) , STAT=izl ) ; clc_con = 0.0_wp  ; ist = ist + izl
  ALLOCATE (  clw_con(ie,je,ke) , STAT=izl ) ; clw_con = 0.0_wp  ; ist = ist + izl
  ALLOCATE (  prr_con(ie,je)    , STAT=izl ) ; prr_con = 0.0_wp  ; ist = ist + izl
  lalloc_prr_con = .TRUE.
  ALLOCATE (  prs_con(ie,je)    , STAT=izl ) ; prs_con = 0.0_wp  ; ist = ist + izl
  lalloc_prs_con = .TRUE.
  ALLOCATE ( prne_con(ie,je)    , STAT=izl ) ; prne_con= 0.0_wp  ; ist = ist + izl
  ALLOCATE (  bas_con(ie,je)    , STAT=izl ) ; bas_con = 0.0_wp  ; ist = ist + izl
  ALLOCATE (  top_con(ie,je)    , STAT=izl ) ; top_con = 0.0_wp  ; ist = ist + izl
  ALLOCATE (  tt_conv(ie,je,ke) , STAT=izl ) ; tt_conv = 0.0_wp  ; ist = ist + izl
  ALLOCATE ( qvt_conv(ie,je,ke) , STAT=izl ) ; qvt_conv= 0.0_wp  ; ist = ist + izl
  ALLOCATE ( qct_conv(ie,je,ke) , STAT=izl ) ; qct_conv= 0.0_wp  ; ist = ist + izl
  ALLOCATE ( qit_conv(ie,je,ke) , STAT=izl ) ; qit_conv= 0.0_wp  ; ist = ist + izl
  ALLOCATE ( qrt_conv(ie,je,ke) , STAT=izl ) ; qrt_conv= 0.0_wp  ; ist = ist + izl
  ALLOCATE ( qst_conv(ie,je,ke) , STAT=izl ) ; qst_conv= 0.0_wp  ; ist = ist + izl
  ALLOCATE (  ut_conv(ie,je,ke) , STAT=izl ) ; ut_conv = 0.0_wp  ; ist = ist + izl
  ALLOCATE (  vt_conv(ie,je,ke) , STAT=izl ) ; vt_conv = 0.0_wp  ; ist = ist + izl
  ALLOCATE (tket_conv(ie,je,ke1), STAT=izl ) ;tket_conv= 0.0_wp  ; ist = ist + izl
  ALLOCATE (tket_hshr(ie,je,ke1), STAT=izl ) ;tket_hshr= 0.0_wp  ; ist = ist + izl
  ALLOCATE (tket_sso (ie,je,ke1), STAT=izl ) ;tket_sso = 0.0_wp  ; ist = ist + izl
  ALLOCATE ( mflx_con(ie,je)    , STAT=izl ) ; mflx_con= 0.0_wp  ; ist = ist + izl
  ALLOCATE ( cape_con(ie,je)    , STAT=izl ) ; cape_con= 0.0_wp  ; ist = ist + izl
  ALLOCATE (  tke_con(ie,je)    , STAT=izl ) ; tke_con = 0.0_wp  ; ist = ist + izl
  ALLOCATE ( qcvg_con(ie,je)    , STAT=izl ) ; qcvg_con= 0.0_wp  ; ist = ist + izl
  ALLOCATE (    w0avg(ie,je,ke) , STAT=izl ) ; w0avg   = 0.0_wp  ; ist = ist + izl
  ALLOCATE (  qrs    (ie,je,ke) , STAT=izl ) ; qrs     = 0.0_wp  ; ist = ist + izl
  ALLOCATE (  prr_gsp(ie,je)    , STAT=izl ) ; prr_gsp = 0.0_wp  ; ist = ist + izl
  lalloc_prr_gsp = .TRUE.
  ALLOCATE (  prs_gsp(ie,je)    , STAT=izl ) ; prs_gsp = 0.0_wp  ; ist = ist + izl
  lalloc_prs_gsp = .TRUE.
  IF (itype_gscp >= 4) THEN
    ALLOCATE ( prg_gsp(ie,je)   , STAT=izl ) ; prg_gsp = 0.0_wp  ; ist = ist + izl
    lalloc_prg_gsp = .TRUE.
    ALLOCATE ( grau_gsp(ie,je)  , STAT=izl ) ; grau_gsp = 0.0_wp ; ist = ist + izl
    IF (itype_gscp >= 100) THEN 
      ALLOCATE ( prh_gsp(ie,je)   , STAT=izl ); prh_gsp = 0.0_wp ; ist = ist + izl
      ALLOCATE ( hail_gsp(ie,je)  , STAT=izl ); hail_gsp = 0.0_wp; ist = ist + izl
    ENDIF     
  ENDIF
  ! because these 2 fields are used in the turbulence interface
  ! they are allocated anyhow
  ALLOCATE ( ut_sso(ie,je,ke) , STAT=izl ) ; ut_sso  = 0.0_wp  ; ist = ist + izl
  ALLOCATE ( vt_sso(ie,je,ke) , STAT=izl ) ; vt_sso  = 0.0_wp  ; ist = ist + izl
  IF (lsso) THEN
    ALLOCATE ( tt_sso(ie,je,ke) , STAT=izl ) ; tt_sso  = 0.0_wp  ; ist = ist + izl
    ALLOCATE ( ustr_sso (ie,je) , STAT=izl ) ; ustr_sso= 0.0_wp  ; ist = ist + izl
    ALLOCATE ( vstr_sso (ie,je) , STAT=izl ) ; vstr_sso= 0.0_wp  ; ist = ist + izl
    ALLOCATE ( vdis_sso (ie,je) , STAT=izl ) ; vdis_sso= 0.0_wp  ; ist = ist + izl
    ALLOCATE ( austr_sso(ie,je) , STAT=izl ) ; austr_sso=0.0_wp  ; ist = ist + izl
    ALLOCATE ( avstr_sso(ie,je) , STAT=izl ) ; avstr_sso=0.0_wp  ; ist = ist + izl
    ALLOCATE ( avdis_sso(ie,je) , STAT=izl ) ; avdis_sso=0.0_wp  ; ist = ist + izl
  ENDIF
  ALLOCATE ( dqvdt  (ie,je,ke)  , STAT=izl ) ; dqvdt   = 0.0_wp  ; ist = ist + izl
  ALLOCATE ( qvsflx (ie,je)     , STAT=izl ) ; qvsflx  = 0.0_wp  ; ist = ist + izl
  ALLOCATE ( dpsdt  (ie,je)     , STAT=izl ) ; dpsdt   = 0.0_wp  ; ist = ist + izl
  ALLOCATE ( umfl_s (ie,je)     , STAT=izl ) ; umfl_s  = 0.0_wp  ; ist = ist + izl
  ALLOCATE ( vmfl_s (ie,je)     , STAT=izl ) ; vmfl_s  = 0.0_wp  ; ist = ist + izl
  ALLOCATE ( shfl_s (ie,je)     , STAT=izl ) ; shfl_s  = 0.0_wp  ; ist = ist + izl
  ALLOCATE ( lhfl_s (ie,je)     , STAT=izl ) ; lhfl_s  = 0.0_wp  ; ist = ist + izl
  ALLOCATE ( aumfl_s(ie,je)     , STAT=izl ) ; aumfl_s = 0.0_wp  ; ist = ist + izl
  ALLOCATE ( avmfl_s(ie,je)     , STAT=izl ) ; avmfl_s = 0.0_wp  ; ist = ist + izl
  ALLOCATE ( ashfl_s(ie,je)     , STAT=izl ) ; ashfl_s = 0.0_wp  ; ist = ist + izl
  ALLOCATE ( alhfl_s(ie,je)     , STAT=izl ) ; alhfl_s = 0.0_wp  ; ist = ist + izl

  IF (izdebug > 10) THEN
    PRINT *, '   ALLOCATED fields for dynamics/parameterizations:   ', ist
  ENDIF

! fields for model output and diagnostics
! ---------------------------------------

  ALLOCATE ( t_2m    (ie,je) , STAT=izl ) ; t_2m     = 0.0_wp   ; ist = ist + izl
  ALLOCATE ( qv_2m   (ie,je) , STAT=izl ) ; qv_2m    = 0.0_wp   ; ist = ist + izl
  ALLOCATE ( td_2m   (ie,je) , STAT=izl ) ; td_2m    = 0.0_wp   ; ist = ist + izl
  ALLOCATE ( rh_2m   (ie,je) , STAT=izl ) ; rh_2m    = 0.0_wp   ; ist = ist + izl
  ALLOCATE ( u_10m   (ie,je) , STAT=izl ) ; u_10m    = 0.0_wp   ; ist = ist + izl
  ALLOCATE ( v_10m   (ie,je) , STAT=izl ) ; v_10m    = 0.0_wp   ; ist = ist + izl
  IF (lbdclim) THEN
    ALLOCATE (t_2m_av (ie,je), STAT=izl ) ; t_2m_av  = 0.0_wp   ; ist = ist + izl
    ALLOCATE (td_2m_av(ie,je), STAT=izl ) ; td_2m_av = 0.0_wp   ; ist = ist + izl
    ALLOCATE (u_10m_av(ie,je), STAT=izl ) ; u_10m_av = 0.0_wp   ; ist = ist + izl
    ALLOCATE (v_10m_av(ie,je), STAT=izl ) ; v_10m_av = 0.0_wp   ; ist = ist + izl
  ENDIF
  ALLOCATE ( tmin_2m (ie,je) , STAT=izl ) ; tmin_2m  = 500.0_wp ; ist = ist + izl
  ALLOCATE ( tmax_2m (ie,je) , STAT=izl ) ; tmax_2m  = 0.0_wp   ; ist = ist + izl
  ALLOCATE ( vmax_10m(ie,je) , STAT=izl ) ; vmax_10m = 0.0_wp   ; ist = ist + izl
  ALLOCATE ( vabsmx_10m(ie,je),STAT=izl ) ; vabsmx_10m = 0.0_wp ; ist = ist + izl
  ALLOCATE ( vgust_dyn(ie,je), STAT=izl ) ; vgust_dyn= 0.0_wp   ; ist = ist + izl
  ALLOCATE ( vgust_con(ie,je), STAT=izl ) ; vgust_con= 0.0_wp   ; ist = ist + izl
  ALLOCATE ( asob_s  (ie,je) , STAT=izl ) ; asob_s   = 0.0_wp   ; ist = ist + izl
  ALLOCATE ( athb_s  (ie,je) , STAT=izl ) ; athb_s   = 0.0_wp   ; ist = ist + izl
  ALLOCATE ( apab_s  (ie,je) , STAT=izl ) ; apab_s   = 0.0_wp   ; ist = ist + izl
  ALLOCATE ( asob_t  (ie,je) , STAT=izl ) ; asob_t   = 0.0_wp   ; ist = ist + izl
  ALLOCATE ( athb_t  (ie,je) , STAT=izl ) ; athb_t   = 0.0_wp   ; ist = ist + izl
  ALLOCATE (  sod_t  (ie,je) , STAT=izl ) ;  sod_t   = 0.0_wp   ; ist = ist + izl
  ALLOCATE ( asod_t  (ie,je) , STAT=izl ) ; asod_t   = 0.0_wp   ; ist = ist + izl
  ALLOCATE ( dursun  (ie,je) , STAT=izl ) ; dursun   = 0.0_wp   ; ist = ist + izl
  ALLOCATE ( dursun_m(ie,je) , STAT=izl ) ; dursun_m = 0.0_wp   ; ist = ist + izl
  ALLOCATE ( dursun_r(ie,je) , STAT=izl ) ; dursun_r = 0.0_wp   ; ist = ist + izl
  ALLOCATE ( rain_gsp(ie,je) , STAT=izl ) ; rain_gsp = 0.0_wp   ; ist = ist + izl
  ALLOCATE ( snow_gsp(ie,je) , STAT=izl ) ; snow_gsp = 0.0_wp   ; ist = ist + izl
  ALLOCATE ( rain_con(ie,je) , STAT=izl ) ; rain_con = 0.0_wp   ; ist = ist + izl
  ALLOCATE ( snow_con(ie,je) , STAT=izl ) ; snow_con = 0.0_wp   ; ist = ist + izl
  ALLOCATE ( runoff_s(ie,je) , STAT=izl ) ; runoff_s = 0.0_wp   ; ist = ist + izl
  ALLOCATE ( runoff_g(ie,je) , STAT=izl ) ; runoff_g = 0.0_wp   ; ist = ist + izl
  ALLOCATE ( rstom   (ie,je) , STAT=izl ) ; rstom    = 0.0_wp   ; ist = ist + izl
  ALLOCATE ( lhfl_bs(ie,je)  , STAT=izl ) ; lhfl_bs = 0.0_wp    ; ist = ist + izl
  ALLOCATE ( lhfl_pl(ie,je,ke_soil) , STAT=izl ) ; lhfl_pl = 0.0_wp   ; ist = ist + izl
  ALLOCATE ( alhfl_bs(ie,je) , STAT=izl ) ; alhfl_bs = 0.0_wp   ; ist = ist + izl
  ALLOCATE ( alhfl_pl(ie,je,ke_soil) , STAT=izl ) ; alhfl_pl = 0.0_wp   ; ist = ist + izl
  IF (ldiagnos) THEN
    ALLOCATE (tdiv_hum(ie,je), STAT=izl ) ; tdiv_hum = 0.0_wp   ; ist = ist + izl
    ALLOCATE (aevap_s (ie,je), STAT=izl ) ; aevap_s  = 0.0_wp   ; ist = ist + izl
  ENDIF

  IF (luse_rttov) THEN
    ALLOCATE (synme7(ie,je,8), synmsg(ie,je,32))
  ENDIF

  IF (izdebug > 10) THEN
    PRINT *, '   ALLOCATED fields for model output / diagnostics:   ', ist
  ENDIF

! fields for the boundary values
! ------------------------------

  ALLOCATE ( u_bd   (ie,je,ke,2), STAT=izl ) ; u_bd     = 0.0_wp; ist = ist + izl
  ALLOCATE ( v_bd   (ie,je,ke,2), STAT=izl ) ; v_bd     = 0.0_wp; ist = ist + izl
  IF ( .NOT. lw_freeslip ) THEN
    ALLOCATE ( w_bd   (ie,je,ke1,2), STAT=izl ) ; w_bd  = 0.0_wp; ist = ist + izl
  ENDIF
  ALLOCATE ( t_bd   (ie,je,ke,2), STAT=izl ) ; t_bd     = 0.0_wp; ist = ist + izl
  ALLOCATE ( pp_bd  (ie,je,ke,2), STAT=izl ) ; pp_bd    = 0.0_wp; ist = ist + izl
  ALLOCATE ( qv_s_bd   (ie,je,2), STAT=izl ) ; qv_s_bd  = 0.0_wp; ist = ist + izl
  ALLOCATE ( t_snow_bd (ie,je,2), STAT=izl ) ; t_snow_bd= 0.0_wp; ist = ist + izl
  ALLOCATE ( w_snow_bd (ie,je,2), STAT=izl ) ; w_snow_bd= 0.0_wp; ist = ist + izl
  IF (.NOT. lmulti_layer) THEN
    ALLOCATE (t_s_bd   (ie,je,2), STAT=izl ) ; t_s_bd   = 0.0_wp; ist = ist + izl
    lalloc_t_s_bd = .TRUE.
    ALLOCATE (t_m_bd   (ie,je,2), STAT=izl ) ; t_m_bd   = 0.0_wp; ist = ist + izl
    ALLOCATE (w_g1_bd  (ie,je,2), STAT=izl ) ; w_g1_bd  = 0.0_wp; ist = ist + izl
    ALLOCATE (w_g2_bd  (ie,je,2), STAT=izl ) ; w_g2_bd  = 0.0_wp; ist = ist + izl
    IF ( (nlgw_bd == 3) .OR. (nlgw == 3) ) THEN
      ALLOCATE (w_g3_bd(ie,je,2), STAT=izl ) ; w_g3_bd  = 0.0_wp; ist = ist + izl
      lalloc_w_g3_bd = .TRUE.
    ENDIF
  ENDIF
  IF (lbdclim) THEN
    ALLOCATE (plcov_bd (ie,je,2), STAT=izl ) ; plcov_bd = 0.0_wp; ist = ist + izl
    ALLOCATE (lai_bd   (ie,je,2), STAT=izl ) ; lai_bd   = 0.0_wp; ist = ist + izl
    ALLOCATE (rootdp_bd(ie,je,2), STAT=izl ) ; rootdp_bd= 0.0_wp; ist = ist + izl
    ALLOCATE (vio3_bd  (ie,je,2), STAT=izl ) ; vio3_bd  = 0.0_wp; ist = ist + izl
    ALLOCATE (hmo3_bd  (ie,je,2), STAT=izl ) ; hmo3_bd  = 0.0_wp; ist = ist + izl
    ALLOCATE (t_cl_bd  (ie,je,2), STAT=izl ) ; t_cl_bd  = 0.0_wp; ist = ist + izl
    ALLOCATE (w_cl_bd  (ie,je,2), STAT=izl ) ; w_cl_bd  = 0.0_wp; ist = ist + izl
    IF (lmulti_layer) THEN
      ALLOCATE (t_s_bd (ie,je,2), STAT=izl ) ; t_s_bd   = 0.0_wp; ist = ist + izl
      lalloc_t_s_bd = .TRUE.
    ENDIF
    IF (itype_aerosol == 2) THEN
      ALLOCATE (aer_su_bd(ie,je,2),STAT=izl ); aer_su_bd = 0.0_wp; ist = ist + izl
      ALLOCATE (aer_du_bd(ie,je,2),STAT=izl ); aer_du_bd = 0.0_wp; ist = ist + izl
      ALLOCATE (aer_bc_bd(ie,je,2),STAT=izl ); aer_bc_bd = 0.0_wp; ist = ist + izl
      ALLOCATE (aer_or_bd(ie,je,2),STAT=izl ); aer_or_bd = 0.0_wp; ist = ist + izl
      ALLOCATE (aer_ss_bd(ie,je,2),STAT=izl ); aer_ss_bd = 0.0_wp; ist = ist + izl
    ENDIF
  ELSE
    IF (lbdsst .AND. lmulti_layer) THEN
      ! otherwise t_s_bd is already allocated
      ALLOCATE (t_s_bd (ie,je,2), STAT=izl ) ; t_s_bd   = 0.0_wp; ist = ist + izl
      lalloc_t_s_bd = .TRUE.
    ENDIF
  ENDIF

  IF (izdebug > 10) THEN
    PRINT *, '   ALLOCATED boundary fields:                         ', ist
  ENDIF

! turbulence variables
! --------------------

  IF ( lphys .AND. ltur ) THEN
    ! these 2 fields were allocated only in special cases; but now they are 
    ! used in the turbulence interface, so allocate them anyhow
    ALLOCATE(rcld(ie,je,ke1), STAT=izl);           rcld=0.0_wp; ist = ist + izl
    ALLOCATE(edr(ie,je,ke+1), STAT=izl);            edr=0.0_wp; ist = ist + izl
!US allocate it in any case
    ALLOCATE(tket_adv(ie,je,ke+1), STAT=izl);  tket_adv=0.0_wp; ist = ist + izl

    SELECT CASE( itype_turb )
    CASE( 3 )
      ALLOCATE(tke(ie,je,ke+1,nztlev), STAT=izl );    tke=0.0_wp; ist = ist + izl
      ALLOCATE(tketens(ie,je,ke+1), STAT=izl);    tketens=0.0_wp; ist = ist + izl
!US      ALLOCATE(tket_adv(ie,je,ke+1), STAT=izl);  tket_adv=0.0_wp; ist = ist + izl
    CASE( 5:8 )
      IF ( lprog_tke ) THEN
        ALLOCATE(tke(ie,je,ke+1,nztlev), STAT=izl);   tke=0.0_wp; ist = ist + izl
        ALLOCATE(tketens (ie,je,ke+1), STAT=izl); tketens=0.0_wp; ist = ist + izl
      ELSE
        ALLOCATE(tke(ie,je,ke+1,1), STAT=izl );       tke=0.0_wp; ist = ist + izl
      END IF
    CASE default
      ALLOCATE(tke(ie,je,ke+1,1), STAT=izl ) ;        tke=0.0_wp; ist = ist + izl
    END SELECT
    lalloc_tke=.TRUE.
  ENDIF

#ifdef NUDGING
! fields for lhn
! --------------

  IF (llhn .OR. llhnverif) THEN
    ALLOCATE ( tt_lheat(ie,je,ke,nztlev),STAT=izl);tt_lheat = 0.0_wp ; ist=ist+izl
    IF (lhn_qrs) THEN
       ALLOCATE ( qrsflux(ie,je,ke),STAT=izl);qrsflux = 0.0_wp ; ist = ist + izl
    ENDIF
    ALLOCATE ( tt_lheat_o (ie,je,ke) , STAT=izl ) ; tt_lheat_o  = 0.0_wp ; ist=ist+izl
    ALLOCATE ( tinc_lhn_o (ie,je,ke) , STAT=izl ) ; tinc_lhn_o  = 0.0_wp ; ist=ist+izl
    ALLOCATE ( ttm_cv_o   (ie,je,ke) , STAT=izl ) ; ttm_cv_o    = 0.0_wp ; ist=ist+izl
  ENDIF

! time integrated analysis increment fields
! -----------------------------------------

  IF (lout_anai) THEN
    ALLOCATE ( ff_anai(ie,je,ke) , STAT=izl ) ; ff_anai = 0.0_wp ; ist = ist + izl
    ALLOCATE ( dd_anai(ie,je,ke) , STAT=izl ) ; dd_anai = 0.0_wp ; ist = ist + izl
    ALLOCATE ( t_anai (ie,je,ke) , STAT=izl ) ; t_anai  = 0.0_wp ; ist = ist + izl
    ALLOCATE ( p_anai (ie,je,ke) , STAT=izl ) ; p_anai  = 0.0_wp ; ist = ist + izl
    ALLOCATE ( qv_anai(ie,je,ke) , STAT=izl ) ; qv_anai = 0.0_wp ; ist = ist + izl
    ALLOCATE ( qc_anai(ie,je,ke) , STAT=izl ) ; qc_anai = 0.0_wp ; ist = ist + izl
  ENDIF
#endif

  IF (izdebug > 10) THEN
    PRINT *, '   ALLOCATED additional fields:                       ', ist
  ENDIF

!US    !------------------------------------------------------------------------------
!US    ! Section 2:  allocation for canopy fields
!US    !------------------------------------------------------------------------------
!US    
!US    ELSEIF ( yaction == 'canopy' ) THEN
!US    
!US    ! external canopy fields
!US    ! ----------------------
!US    
!US      IF ( ltur ) THEN
!US    
!US        ! Allocation with kcm-1, because kcm was increased by 1 in init_canopy
!US        ALLOCATE ( c_big(ie,je,kcm-1:ke+1), STAT=izl ); c_big = 0.0_wp; ist = ist+izl
!US        ALLOCATE ( c_sml(ie,je,kcm-1:ke+1), STAT=izl ); c_sml = 0.0_wp; ist = ist+izl
!US        ALLOCATE ( r_air(ie,je,kcm-1:ke+1), STAT=izl ); r_air = 0.0_wp; ist = ist+izl
!US        ALLOCATE ( t_e  (ie,je,kcm-1:ke  ), STAT=izl ); t_e   = 0.0_wp; ist = ist+izl
!US        ALLOCATE ( qv_e (ie,je,kcm-1:ke  ), STAT=izl ); qv_e  = 0.0_wp; ist = ist+izl
!US    
!US      ENDIF

ENDIF

!------------------------------------------------------------------------------
!  End of the Subroutine
!------------------------------------------------------------------------------

END SUBROUTINE alloc_meteofields


!+ Module Procedure in "src_allocations"
!+ performing the deallocation of the meteorological fields
!------------------------------------------------------------------------------

!==============================================================================
!+ Deallocates the allocated meteorological fields
!------------------------------------------------------------------------------

SUBROUTINE dealloc_meteofields  (istat)

!------------------------------------------------------------------------------
!
! Description:
!   This routine deallocates the allocated meteorological fields.
!
! Method:
!   DEALLOCATE statement
!
!------------------------------------------------------------------------------
!
! Declarations must be of the form:
! <type> <VariableName> ! Description/ purpose of variable
!
! Subroutine / Function arguments
INTEGER (KIND=iintegers), INTENT(OUT)   ::       &
  istat              ! for local error-code

INTEGER (KIND=iintegers)                ::       &
  izdebug            ! for local debug output

!- End of header
!==============================================================================

!------------------------------------------------------------------------------
!- Begin SUBROUTINE dealloc_meteofields
!------------------------------------------------------------------------------

istat = 0

IF (lprintdeb_all) THEN
  izdebug = idbg_level
ELSE
  IF (my_cart_id == 0) THEN
    izdebug = idbg_level
  ELSE
    izdebug = 0
  ENDIF
ENDIF

! constant fields for the reference atmosphere
! --------------------------------------------

  DEALLOCATE ( hhl_prof, STAT=istat )
  DEALLOCATE ( rho0    , STAT=istat )
  DEALLOCATE ( dp0     , STAT=istat )
  DEALLOCATE ( p0      , STAT=istat )
  DEALLOCATE ( p0hl    , STAT=istat )
  DEALLOCATE ( t0hl    , STAT=istat )
  DEALLOCATE ( hhl     , STAT=istat )
  DEALLOCATE ( dt0dz   , STAT=istat )
  DEALLOCATE ( t0      , STAT=istat )

! external parameter fields
! -------------------------

  DEALLOCATE ( hsurf  , STAT=istat )
  IF (lsso) THEN
    DEALLOCATE ( sso_stdh , STAT=istat )
    DEALLOCATE ( sso_gamma, STAT=istat )
    DEALLOCATE ( sso_theta, STAT=istat )
    DEALLOCATE ( sso_sigma, STAT=istat )
  ENDIF
  IF (itype_aerosol == 2) THEN
    DEALLOCATE ( aer_su   , STAT=istat )
    DEALLOCATE ( aer_du   , STAT=istat )
    DEALLOCATE ( aer_or   , STAT=istat )
    DEALLOCATE ( aer_bc   , STAT=istat )
    DEALLOCATE ( aer_ss   , STAT=istat )
  ENDIF
  IF (lemiss) THEN
    DEALLOCATE ( emis_rad , STAT=istat )
  ENDIF
  IF (lstomata) THEN
    DEALLOCATE ( rsmin2d  , STAT=istat )
  ENDIF
  DEALLOCATE ( gz0    , STAT=istat )
  DEALLOCATE ( fr_land, STAT=istat )
  DEALLOCATE ( soiltyp, STAT=istat )
  DEALLOCATE ( vio3   , STAT=istat )
  DEALLOCATE ( hmo3   , STAT=istat )
  DEALLOCATE ( rlat   , STAT=istat )
  DEALLOCATE ( rlon   , STAT=istat )
  IF (ALLOCATED (rlattot)) THEN
    DEALLOCATE ( rlattot, STAT=istat )
    DEALLOCATE ( rlontot, STAT=istat )
  ENDIF
  DEALLOCATE ( fc     , STAT=istat )
  IF (ALLOCATED (fccos)) THEN
    DEALLOCATE ( fccos  , STAT=istat )
  ENDIF
  DEALLOCATE ( rmy    , STAT=istat )
  DEALLOCATE ( rmyq   , STAT=istat )
  DEALLOCATE ( hd_mask_dcoeff_p, STAT=istat )
  DEALLOCATE ( hd_mask_dcoeff_t, STAT=istat )
  DEALLOCATE ( hd_mask_dcoeff_q, STAT=istat )
  DEALLOCATE ( hd_mask_dcoeff_u, STAT=istat )
  DEALLOCATE ( ofa_hdx, STAT=istat )
  DEALLOCATE ( ofa_hdy, STAT=istat )
  DEALLOCATE ( hd_mask, STAT=istat )
  DEALLOCATE ( least_lbdz, STAT=istat )
  DEALLOCATE ( lwest_lbdz, STAT=istat )
  DEALLOCATE ( lnorth_lbdz,STAT=istat )
  DEALLOCATE ( lsouth_lbdz,STAT=istat )
  DEALLOCATE ( crlat  , STAT=istat )
  DEALLOCATE ( acrlat , STAT=istat )
  DEALLOCATE ( tgrlat , STAT=istat )
  DEALLOCATE ( aerlan , STAT=istat )
  DEALLOCATE ( aerurb , STAT=istat )
  DEALLOCATE ( aerdes , STAT=istat )
  DEALLOCATE ( aersea , STAT=istat )
  DEALLOCATE ( plcov  , STAT=istat )
  DEALLOCATE ( lai    , STAT=istat )
  DEALLOCATE ( tai    , STAT=istat )
  DEALLOCATE ( sai    , STAT=istat )
  DEALLOCATE ( eai    , STAT=istat )
  DEALLOCATE ( rootdp , STAT=istat )
  DEALLOCATE ( llandmask , STAT=istat )

  IF (lseaice) THEN
    DEALLOCATE ( lseamask, STAT=istat )
  ENDIF

  DEALLOCATE ( h_can   , STAT=istat )
  DEALLOCATE ( d_pat   , STAT=istat )

  IF (ALLOCATED (for_e)) THEN
    DEALLOCATE ( for_e   , STAT=istat )
    DEALLOCATE ( for_d   , STAT=istat )
  ENDIF

  DEALLOCATE ( depth_lk, STAT=istat )
  DEALLOCATE ( fr_lake , STAT=istat )
  IF (llake) THEN
    DEALLOCATE ( fetch_lk, STAT=istat )
    DEALLOCATE ( dp_bs_lk, STAT=istat )
    DEALLOCATE ( t_bs_lk , STAT=istat )
    DEALLOCATE ( gamso_lk, STAT=istat )
  ENDIF

  IF     (itype_albedo == 2) THEN
    DEALLOCATE ( alb_dry , STAT=istat )
    DEALLOCATE ( alb_sat , STAT=istat )
  ELSEIF (itype_albedo == 3) THEN
    DEALLOCATE ( alb_dif , STAT=istat )
  ENDIF

  IF (izdebug > 10) THEN
    PRINT *, '   DE-ALLOCATED constant and external parameter fields:  ', istat
  ENDIF

! prognostic variables
! --------------------

  DEALLOCATE ( u      , STAT=istat )
  DEALLOCATE ( v      , STAT=istat )
  DEALLOCATE ( w      , STAT=istat )
  DEALLOCATE ( t      , STAT=istat )
  DEALLOCATE ( pp     , STAT=istat )
  DEALLOCATE ( tke    , STAT=istat )

! tendency fields for the prognostic variables
! --------------------------------------------

  DEALLOCATE ( utens  , STAT=istat )
  DEALLOCATE ( vtens  , STAT=istat )
  DEALLOCATE ( wtens  , STAT=istat )
  DEALLOCATE ( ttens  , STAT=istat )
  DEALLOCATE ( pptens , STAT=istat )
  IF ( ALLOCATED(tketens)) THEN
     DEALLOCATE ( tketens , STAT=istat )
  ENDIF
  IF ( ALLOCATED(tket_adv)) THEN
     DEALLOCATE ( tket_adv , STAT=istat )
  ENDIF

  IF (izdebug > 10) THEN
    PRINT *, '   DE-ALLOCATED prognostic variables and tendencies:  ', istat
  ENDIF

! fields for surface values and soil model variables
! --------------------------------------------------

  DEALLOCATE ( ps        , STAT=istat )
  DEALLOCATE ( t_snow    , STAT=istat )
  DEALLOCATE ( t_s       , STAT=istat )
  DEALLOCATE ( t_g       , STAT=istat )
  DEALLOCATE ( tg_radstep,STAT=istat )
  DEALLOCATE ( qv_s      , STAT=istat )
  DEALLOCATE ( w_snow    , STAT=istat )
  DEALLOCATE ( w_i       , STAT=istat )
  DEALLOCATE ( h_snow    , STAT=istat )
  DEALLOCATE ( snow_melt , STAT=istat )
  DEALLOCATE ( fr_snow   , STAT=istat )
  DEALLOCATE ( fr_wi     , STAT=istat )
  DEALLOCATE ( ustar_fv  , STAT=istat )

  IF (lmulti_layer) THEN
    DEALLOCATE ( t_so     , STAT=istat )
    DEALLOCATE ( w_so     , STAT=istat )
    DEALLOCATE ( w_so_ice , STAT=istat )
    DEALLOCATE ( freshsnow, STAT=istat )
    DEALLOCATE ( rho_snow , STAT=istat )
    IF(lmulti_snow) THEN
      DEALLOCATE ( t_snow_mult  , STAT=istat )
      DEALLOCATE ( dzh_snow_mult, STAT=istat )
      DEALLOCATE ( w_snow_mult  , STAT=istat )
      DEALLOCATE ( wliq_snow    , STAT=istat )
      DEALLOCATE ( rho_snow_mult, STAT=istat )
    ENDIF
    IF (lbdclim) THEN
      DEALLOCATE ( t_cl   , STAT=istat )
      DEALLOCATE ( w_cl   , STAT=istat )
    ENDIF
  ELSE
    DEALLOCATE ( t_m    , STAT=istat )
    DEALLOCATE ( t_cl   , STAT=istat )
    DEALLOCATE ( w_g1   , STAT=istat )
    DEALLOCATE ( w_g2   , STAT=istat )
    IF ( (nlgw_ini == 3) .OR. (nlgw == 3) ) THEN
      DEALLOCATE ( w_g3   , STAT=istat )
    ENDIF
    DEALLOCATE ( w_cl   , STAT=istat )
  ENDIF

  IF (lseaice .OR. llake) THEN
    DEALLOCATE ( t_ice    , STAT=istat )
    DEALLOCATE ( h_ice    , STAT=istat )
  ENDIF

  IF (llake) THEN
    DEALLOCATE ( t_mnw_lk , STAT=istat )
    DEALLOCATE ( t_wml_lk , STAT=istat )
    DEALLOCATE ( t_bot_lk , STAT=istat )
    DEALLOCATE ( t_b1_lk  , STAT=istat )
    DEALLOCATE ( c_t_lk   , STAT=istat )
    DEALLOCATE ( h_ml_lk  , STAT=istat )
    DEALLOCATE ( h_b1_lk  , STAT=istat )
  ENDIF

  IF (izdebug > 10) THEN
    PRINT *, '   DE-ALLOCATED prognostic variables and tendencies:  ', istat
  ENDIF

! fields that are computed in the parametrization and dynamics
! ------------------------------------------------------------

  IF ( ldiabf_lh ) THEN
    DEALLOCATE ( tinc_lh, STAT=istat )
  ENDIF
  DEALLOCATE ( rho    , STAT=istat )

  ! coefficients for turbulent diffusion in the atmosphere
  ! (defined on half levels)
  DEALLOCATE ( tkvm   , STAT=istat )
  DEALLOCATE ( tkvh   , STAT=istat )
  IF ( l3dturb ) THEN
    DEALLOCATE ( tkhm   , STAT=istat )
    DEALLOCATE ( tkhh   , STAT=istat )
  ENDIF
  IF ( ALLOCATED(rcld)) THEN
     DEALLOCATE (rcld   , STAT=istat )
  ENDIF 
  IF (ALLOCATED(edr)) THEN
     DEALLOCATE ( edr    , STAT=istat )
  ENDIF
  DEALLOCATE ( tcm    , STAT=istat )
  DEALLOCATE ( tch    , STAT=istat )
#ifdef COUP_OAS_COS
  DEALLOCATE ( tcw    , STAT=istat )
  DEALLOCATE (trad_clm, STAT=istat )
#endif
  DEALLOCATE ( tfm    , STAT=istat )
  DEALLOCATE ( tfh    , STAT=istat )
  DEALLOCATE ( tfv    , STAT=istat )

  ! fields from the radiation scheme
  DEALLOCATE ( sohr   , STAT=istat )
  DEALLOCATE ( sotr   , STAT=istat )
  DEALLOCATE ( sotr_par,STAT=istat )
  DEALLOCATE ( thhr   , STAT=istat )
  DEALLOCATE ( clc_sgs, STAT=istat )
  DEALLOCATE ( alb_rad, STAT=istat )
  IF ( nradcoarse > 1 ) THEN
    DEALLOCATE ( alb_rad_coarse, STAT=istat )
  ENDIF
  DEALLOCATE ( sobs   , STAT=istat )
  DEALLOCATE ( thbs   , STAT=istat )
  DEALLOCATE ( pabs   , STAT=istat )
  DEALLOCATE ( sobt   , STAT=istat )
  DEALLOCATE ( thbt   , STAT=istat )
  DEALLOCATE ( clch   , STAT=istat )
  DEALLOCATE ( clcm   , STAT=istat )
  DEALLOCATE ( clcl   , STAT=istat )
  DEALLOCATE ( clct   , STAT=istat )
  DEALLOCATE ( sodwddm, STAT=istat )
  DEALLOCATE ( qc_rad , STAT=istat )
  DEALLOCATE ( qi_rad , STAT=istat )

  IF (lradtopo) THEN
    DEALLOCATE (skyview  , STAT=istat )
    DEALLOCATE (slo_asp  , STAT=istat )
    DEALLOCATE (slo_ang  , STAT=istat )
    DEALLOCATE (horizon  , STAT=istat )
  ENDIF
  DEALLOCATE ( swdir_s  , STAT=istat )
  DEALLOCATE ( swdifd_s , STAT=istat )
  DEALLOCATE ( swdifu_s , STAT=istat )
  DEALLOCATE ( swtrdir_s, STAT=istat )
  DEALLOCATE ( swtrdifd_s,STAT=istat )
  DEALLOCATE ( swtrdifu_s,STAT=istat )
  DEALLOCATE ( lwd_s    , STAT=istat )
  DEALLOCATE ( lwu_s    , STAT=istat )
  DEALLOCATE ( aswdir_s , STAT=istat )
  DEALLOCATE ( aswdifd_s, STAT=istat )
  DEALLOCATE ( aswdifu_s, STAT=istat )
  DEALLOCATE ( alwd_s   , STAT=istat )
  DEALLOCATE ( alwu_s   , STAT=istat )
  DEALLOCATE ( swdir_cor, STAT=istat )
  DEALLOCATE ( sun_el , STAT=istat )
  DEALLOCATE ( sun_azi, STAT=istat )

#ifdef TWOMOM_SB
  DEALLOCATE ( reffc_out , STAT=istat )
  DEALLOCATE ( reffi_out , STAT=istat )
  DEALLOCATE ( odepthw_so, STAT=istat )
  DEALLOCATE ( odepthi_so, STAT=istat )
  DEALLOCATE ( odepthw_th, STAT=istat )
  DEALLOCATE ( odepthi_th, STAT=istat )
#endif

  ! fields from the convection scheme
  DEALLOCATE ( clc_con, STAT=istat )
  DEALLOCATE ( clw_con, STAT=istat )
  DEALLOCATE ( prr_con, STAT=istat )
  DEALLOCATE ( prs_con, STAT=istat )
  DEALLOCATE (prne_con, STAT=istat )
  DEALLOCATE ( bas_con, STAT=istat )
  DEALLOCATE ( top_con, STAT=istat )
  DEALLOCATE ( tt_conv, STAT=istat )
  DEALLOCATE (qvt_conv, STAT=istat )
  DEALLOCATE (qct_conv, STAT=istat )
  DEALLOCATE (qit_conv, STAT=istat )
  DEALLOCATE (qrt_conv, STAT=istat )
  DEALLOCATE (qst_conv, STAT=istat )
  DEALLOCATE ( ut_conv, STAT=istat )
  DEALLOCATE ( vt_conv, STAT=istat )
  DEALLOCATE (mflx_con, STAT=istat )
  DEALLOCATE (cape_con, STAT=istat )
  DEALLOCATE ( tke_con, STAT=istat )
  DEALLOCATE (qcvg_con, STAT=istat )
  DEALLOCATE (   w0avg, STAT=istat )

  ! fields from the grid-scale precipitation scheme
  DEALLOCATE ( qrs    , STAT=istat )
  DEALLOCATE ( prr_gsp, STAT=istat )
  DEALLOCATE ( prs_gsp, STAT=istat )
  IF (itype_gscp >= 4) THEN
    DEALLOCATE ( prg_gsp, STAT=istat )
    DEALLOCATE ( grau_gsp, STAT=istat )
    IF (itype_gscp >= 100) THEN
      DEALLOCATE ( prh_gsp, STAT=istat )
      DEALLOCATE ( hail_gsp, STAT=istat )
    ENDIF
  ENDIF

  ! because these 2 fields are used in the turbulence interface
  DEALLOCATE ( ut_sso   , STAT=istat )
  DEALLOCATE ( vt_sso   , STAT=istat )
  IF (lsso) THEN
    DEALLOCATE ( tt_sso   , STAT=istat )
    DEALLOCATE ( ustr_sso , STAT=istat )
    DEALLOCATE ( vstr_sso , STAT=istat )
    DEALLOCATE ( vdis_sso , STAT=istat )
    DEALLOCATE ( austr_sso, STAT=istat )
    DEALLOCATE ( avstr_sso, STAT=istat )
    DEALLOCATE ( avdis_sso, STAT=istat )
  ENDIF
  ! fields that are computed in the dynamics
  DEALLOCATE ( dqvdt  , STAT=istat )
  DEALLOCATE ( qvsflx , STAT=istat )
  DEALLOCATE ( dpsdt  , STAT=istat )
  DEALLOCATE ( umfl_s , STAT=istat )
  DEALLOCATE ( vmfl_s , STAT=istat )
  DEALLOCATE ( shfl_s , STAT=istat )
  DEALLOCATE ( lhfl_s , STAT=istat )
  DEALLOCATE ( aumfl_s, STAT=istat )
  DEALLOCATE ( avmfl_s, STAT=istat )
  DEALLOCATE ( ashfl_s, STAT=istat )
  DEALLOCATE ( alhfl_s, STAT=istat )

  ! fields of the canopy
  IF (ALLOCATED(c_big))    DEALLOCATE ( c_big  , STAT=istat )
  IF (ALLOCATED(c_sml))    DEALLOCATE ( c_sml  , STAT=istat )
  IF (ALLOCATED(r_air))    DEALLOCATE ( r_air  , STAT=istat )
  IF (ALLOCATED(t_e  ))    DEALLOCATE ( t_e    , STAT=istat )
  IF (ALLOCATED(qv_e ))    DEALLOCATE ( qv_e   , STAT=istat )

  IF (izdebug > 10) THEN
    PRINT *, '   DE-ALLOCATED fields for dynamics/parameterizations: ', istat
  ENDIF

! fields for model output and diagnostics
! ---------------------------------------

  DEALLOCATE ( t_2m    , STAT=istat )
  DEALLOCATE ( qv_2m   , STAT=istat )
  DEALLOCATE ( td_2m   , STAT=istat )
  DEALLOCATE ( rh_2m   , STAT=istat )
  DEALLOCATE ( u_10m   , STAT=istat )
  DEALLOCATE ( v_10m   , STAT=istat )
  IF (lbdclim) THEN
    DEALLOCATE ( t_2m_av , STAT=istat )
    DEALLOCATE ( td_2m_av, STAT=istat )
    DEALLOCATE ( u_10m_av, STAT=istat )
    DEALLOCATE ( v_10m_av, STAT=istat )
  ENDIF
  DEALLOCATE ( tmin_2m , STAT=istat )
  DEALLOCATE ( tmax_2m , STAT=istat )
  DEALLOCATE ( vmax_10m, STAT=istat )
  DEALLOCATE ( vabsmx_10m,STAT=istat)
  DEALLOCATE ( vgust_dyn,STAT=istat )
  DEALLOCATE ( vgust_con,STAT=istat )
  DEALLOCATE ( asob_s  , STAT=istat )
  DEALLOCATE ( athb_s  , STAT=istat )
  DEALLOCATE ( apab_s  , STAT=istat )
  DEALLOCATE ( asob_t  , STAT=istat )
  DEALLOCATE ( athb_t  , STAT=istat )
  DEALLOCATE (  sod_t  , STAT=istat )
  DEALLOCATE ( asod_t  , STAT=istat )
  DEALLOCATE ( dursun  , STAT=istat )
  DEALLOCATE ( dursun_m , STAT=istat)
  DEALLOCATE ( dursun_r , STAT=istat)
  DEALLOCATE ( rain_gsp, STAT=istat )
  DEALLOCATE ( snow_gsp, STAT=istat )
  DEALLOCATE ( rain_con, STAT=istat )
  DEALLOCATE ( snow_con, STAT=istat )
  DEALLOCATE ( runoff_s, STAT=istat )
  DEALLOCATE ( runoff_g, STAT=istat )
  DEALLOCATE ( rstom   , STAT=istat )
  DEALLOCATE ( lhfl_bs , STAT=istat )
  DEALLOCATE ( lhfl_pl , STAT=istat )
  DEALLOCATE ( alhfl_bs, STAT=istat )
  DEALLOCATE ( alhfl_pl, STAT=istat )
  IF (ldiagnos) THEN
    DEALLOCATE ( tdiv_hum, STAT=istat )
    DEALLOCATE ( aevap_s , STAT=istat )
  ENDIF

  IF (luse_rttov) THEN
    DEALLOCATE (synme7, synmsg)
  ENDIF

  IF (izdebug > 10) THEN
    PRINT *, '   DE-ALLOCATED fields for model output / diagnostics: ', istat
  ENDIF

! fields for the boundary values
! ------------------------------

  DEALLOCATE ( u_bd     , STAT=istat )
  DEALLOCATE ( v_bd     , STAT=istat )
  IF ( .NOT. lw_freeslip ) THEN
    DEALLOCATE ( w_bd   , STAT=istat )
  ENDIF
  DEALLOCATE ( t_bd     , STAT=istat )
  DEALLOCATE ( pp_bd    , STAT=istat )
  DEALLOCATE ( qv_s_bd  , STAT=istat )
  DEALLOCATE ( t_snow_bd, STAT=istat )
  DEALLOCATE ( w_snow_bd, STAT=istat )

  IF (.NOT. lmulti_layer) THEN
    DEALLOCATE ( t_s_bd   , STAT=istat )
    DEALLOCATE ( t_m_bd   , STAT=istat )
    DEALLOCATE ( w_g1_bd  , STAT=istat )
    DEALLOCATE ( w_g2_bd  , STAT=istat )
    IF ( (nlgw_bd == 3) .OR. (nlgw == 3) ) THEN
      DEALLOCATE ( w_g3_bd , STAT=istat )
    ENDIF
  ENDIF
  IF (lbdclim) THEN
    DEALLOCATE ( plcov_bd , STAT=istat )
    DEALLOCATE ( lai_bd   , STAT=istat )
    DEALLOCATE ( rootdp_bd, STAT=istat )
    DEALLOCATE ( vio3_bd  , STAT=istat )
    DEALLOCATE ( hmo3_bd  , STAT=istat )
    DEALLOCATE ( t_cl_bd  , STAT=istat )
    DEALLOCATE ( w_cl_bd  , STAT=istat )
    IF (lmulti_layer) THEN
      DEALLOCATE ( t_s_bd   , STAT=istat )
    ENDIF
    IF (itype_aerosol == 2) THEN
      DEALLOCATE ( aer_su_bd, STAT=istat )
      DEALLOCATE ( aer_du_bd, STAT=istat )
      DEALLOCATE ( aer_bc_bd, STAT=istat )
      DEALLOCATE ( aer_or_bd, STAT=istat )
      DEALLOCATE ( aer_ss_bd, STAT=istat )
    ENDIF
  ELSE
    IF (lbdsst .AND. lmulti_layer) THEN
      DEALLOCATE (t_s_bd  , STAT=istat )
    ENDIF
  ENDIF

  IF (izdebug > 10) THEN
    PRINT *, '   DE-ALLOCATED boundary fields:                      ', istat
  ENDIF

#ifdef NUDGING
! fields for latent heating
! -------------------------

  IF (llhn .OR. llhnverif) THEN
    DEALLOCATE ( tt_lheat , STAT=istat )
    IF (lhn_qrs) THEN
       DEALLOCATE ( qrsflux , STAT=istat )
    ENDIF
    DEALLOCATE ( tt_lheat_o , STAT=istat )
    DEALLOCATE ( tinc_lhn_o , STAT=istat )
    DEALLOCATE ( ttm_cv_o   , STAT=istat )
  ENDIF

! time integrated analysis increment fields
! -----------------------------------------

  IF (lout_anai) THEN
    DEALLOCATE ( ff_anai , STAT=istat )
    DEALLOCATE ( dd_anai , STAT=istat )
    DEALLOCATE ( t_anai  , STAT=istat )
    DEALLOCATE ( p_anai  , STAT=istat )
    DEALLOCATE ( qv_anai , STAT=istat )
    DEALLOCATE ( qc_anai , STAT=istat )
  ENDIF
#endif

  IF (izdebug > 10) THEN
    PRINT *, '   DE-ALLOCATED additional fields:                    ', istat
  ENDIF

!------------------------------------------------------------------------------
!  End of the Subroutine
!------------------------------------------------------------------------------

END SUBROUTINE dealloc_meteofields
#endif
!==============================================================================

END MODULE src_allocation

