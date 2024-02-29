!+ Data module for variables controlling the run of the model
!------------------------------------------------------------------------------

MODULE data_runcontrol

!------------------------------------------------------------------------------
!
! Description:
!  This module contains variables for running and controlling the forecast.
!  Concerned are the organization of the forecast and the grib I/O.
!  The variables are divided into several groups:
!    
!    - start and end of the forecast
!    - boundary definition and update
!    - controlling the physics
!    - controlling the dynamics
!    - controlling the nudging
!    - controlling the upper boundary condition
!    - additional control variables
!    - controlling the grib I/O
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
! 1.2        1998/03/30 Ulrich Schaettler
!  Introduction of Namelist variable lcond for switching on/off condensation
! 1.3        1998/04/15 Guenther Doms
!  Introduction of Namelist variabel nincconv (timestep increment for convection)
! 1.4        1998/05/22 Guenther Doms
!  Introduction of Namelist variable l2tls for switching on/off 
!  the two timelevel RK-integration scheme
! 1.5        1998/06/29 Guenther Doms
!  Introduction of Namelist variable rdheight and nrdtaur to control
!  the Rayleigh damping layer at the upper boundary.
! 1.9        1998/09/16 Guenther Doms
!  The namelist variable 'nincmxn' specifying the time averaging intervall
!  for certain output fields has been replaces by 'nincmxt' and 'nincmxu'
!  to allow for different averaging periods.
! 1.10       1998/09/29 Ulrich Schaettler
!  Added new control variables for semi-imp. time stepping, nudging, llm.
! 1.11       1998/10/13 Michael Buchhold
!  Specification of initial fields to be checked for the time range indicator.
! 1.17       1998/11/17 Ulrich Schaettler
!  New control variables for reading and writing ready files
! 1.19       1998/12/11 Christoph Schraff
!  Initial field to be checked corrected from VMO3 to HMO3.
! 1.21       1999/01/25 Guenhter Doms
!  Character variables in list 'yunaman' corrected to new names
! 1.29       1999/05/11 Reinhold Hess
!  Check soil water levels for additional element number
! 1.30       1999/06/24 Matthias Raschendofer
!  Introduction of the time index ntke belonging to the TKE-field.
!  Introduction of 9 INTEGER-namelist-parameters controlling the physics:
!  These are: itype_(wcld, tran, turb, synd), imode_(tran, turb), icldm_(rad, tran, turb).
!  Introduction of 4 LOGICAL-namelist-parameters controlling the physics:
!  These are: lturhor, lexpcor, lnonloc, lcpfluc.
!  Introduction of 3 REAL-namelist-parameters controlling the physics:
!  These are: lam_h, lam_m, pat_len.
!  Introduction of a LOGICAL namelist-parameter controlling the output (lgpspec)
!  Introduction of 3 LOGICAL Namelist Parameters for the convection closure and
!  the soil model :lcape, lctke, lbats
! 1.32       1999/08/24 Guenther Doms
!  New logical control variable 'l2dim' for 2D-runs added.
! 1.33       1999/10/14 Matthias Raschendorfer
!  Introduction of 2 LOGICAL namelist-parameter controlling the physics 
!  (ltmpcor,lprfcor).
!  Removal of a LOGICAL namelist-parameter (lbats).
!  Introduction of 2 INTEGER-namelist-parameters controlling the evaporation:
!  These are: itype_(trvg, evsl).
!  Introduction of a REAL-namelist-parameter (crsmin) to control transpiration.
!  Introduction of 7 REAL-namelist-parameter controlling the turbulence:
!  (tur_len, a_heat, d_heat, a_mom, d_mom, c_diff, rat_lam, rat_can, 
!   c_lnd, c_see)
! 1.34       1999/12/10 Ulrich Schaettler
!  Added variables for unit numbers for output files
! 1.39       2000/05/03 Ulrich Schaettler
!  Removed some organizational data (--> now in organize_diagnosis)
!  Included switch ldfi             (in data_filter before)
!  Introduced switch lw_freeslip    (for treatment of w in nesting)
! 2.2        2000/08/18 Guenther Doms
!  Introduction of the logical switch 'lconf_avg' on Namelist input to
!  enable (default) or disable a horizontal averaging of the convective
!  forcing functions. Also, two new REAL namelist input parameters
!  'c_soil' and 'e_surf' to control surface fluxes have been introduced. 
! 2.8        2001/07/06 Ulrich Schaettler
!  Added new variables for multi-layer soil model and moved some variables to 
!  data_io.f90. Added variable nvers for documenting purposes in Grib-Code
! 2.9        2001/07/16 Guenther Doms
!  Introduction of new global contol parameters (for Namelist input in
!  group INPUT_DYN to control horizontal diffusion: itype_hdiff, hd_corr_t,
!  hd_corr_q and hd_dhmax.
! 2.11       2001/09/28 Ulrich Schaettler
!  Added another variable lmelt_var for multi-layer soil model
! 2.15       2002/03/20 Matthias Raschendorfer
!  Introduction of a REAL namelist-parameter controlling the physics (z0m_dia)
! 2.17       2002/05/08 Ulrich Schaettler
!  New Namelist parameters lkainfri, ltiedtke for choosing a convection scheme
! 3.7        2004/02/18 Ulrich Schaettler
!  New Namelist parameters for computing synthetic satellite images,
!    for turning on/off prognostic precipitation
!    for control variables for 2 tl scheme         and
!    for the ratio of laminar scaling factors for heat over sea and land
! 3.10       2004/06/23 Ulrich Schaettler
!  Added variable isynsat_stat for the status of synthetic satellite images
! 3.12       2004/09/15 Christoph Schraff
!  New Namelist variable ldiniprec for initialising prognostic rain / snow.
! 3.13       2004/12/03 Ulrich Schaettler
!  New Namelist variable for 3D turbulence and prognostic treatment of TKE
!  in Runge-Kutta scheme (Jochen Foerstner)
!  New Namelist variable for explicit formulation of lateral boundary
!  relaxation (Jochen Foerstner)
!  New Namelist variable for option of Rayleigh damping (Lucio Torrisi)
! 3.16       2005/07/22 Matthias Raschendorfer, Jochen Foerstner
!  New Namelist parameters
!   - for Turbulence scheme: tkesmot, wichfakt, securi, tkhmin, tkmmin
!   - for Convection:        lshallow   (for shallow convection)
!                            lconv_inst (output of top_con/bas_con)
!   - for Radiation:         lforest (for use of external fields for_e, for_d)
!   - for Runge-Kutta scheme: lsl_adv_qx, lva_impl_dyn, lva_impl_qvqc
!                             (lvertad_impl eliminated)
!                             yef_adv_qx, nbl_exchg, ieva_order, intcr_max
! 3.18       2006/03/03 Ulrich Schaettler
!  Restart: moved some variables from src_meanvalues, to save them here
!  New variables and control switches for CLM Version
!  New switch "llake" for the lake model FLake (Dmitrii Mironov)
!  New switches for very high resolution: l3dturb_metr, ldyn_bbc, ldiabf_lh
!      ldiabf_satad (only internally) (Jochen Foerstner, et.al.)
!      eliminated lva_impl_qvqc
! 3.21       2006/12/04 Ulrich Schaettler (et al)
!  Introduced control variables for Ensemble Prediction   (C. Gebhardt)
!  Introduced additional control variable for FLake model (D. Mironov)
!  Removed some tuning variables and put them to data_turbulence (U. Schaettler)
!  New Namelist variables for radiation on coarser grid: nradcoarse, lradf_avg
!                                                         (T. Reinhardt)
!  New Namelist variables for dynamics: crltau, lhdiff_mask (J. Foerstner)
!  Variables (lcori_deep, ladv_deep) introduced for deep atmosphere (R.Petrik)
!  New Namelist variable lbechtol for Bechtold convection scheme (Meteoswiss)
!  New variable for diagnostics of nudging: lout_anai (Ch. Schraff)
!  New Namelist variable nfi_spubc2 (for itype_spubc=2) (Lucio Torrisi)
! V3_23        2007/03/30 Ulrich Schaettler
!  Added idbg_level and logical ldebug_xxx for controlling verbosity of output
!  xxx stands for one of the components of the model
!  Added some variables for controlling boundary updates
!  Added NL variables lradtopo, nhori for topographic radiation correction
!     (by Matteo Buzzi)
!  Added NL variable itype_lbcqx for lateral boundary treatment of qr,qs,qg
!     (by Jochen Foerstner)
! V3_24        2007/04/26 Ulrich Schaettler
!  Eliminated nincmxu, nincmxt and introduced control as for other increments
! V4_4         2008/07/16 Ulrich Schaettler
!  Eliminated ltime_mean, ltime_proc; replaced by itype_timing
!  New Namelist variables itype_diag_t2m, itype_diag_gusts (DIACTL), lsso (PHYCTL)
!  Changed lyear_360 to itype_calendar (RUNCTL)
!  Replaced logical switches for convection (ltiedtke...) by itype_conv
! V4_5         2008/09/10 Ulrich Schaettler
!  Added parameter nincsso
!  Add namelist options to facilitate idealized simulations (G. Zaengl)
!  Add namelist parameters for modifying values of lai, plcov, rootdp in
!   ensemble mode (Christoph Gebhardt)
! V4_8         2009/02/16 Ulrich Schaettler
!  Renamed logical switch lgen to lartif_data (similar to GME)
!  Added NL switch linit_fields: whether initialization should be done or not
!  New variable nfinalstop for end of total simulation
! V4_9         2009/07/16 Ulrich Schaettler
!  New NL switches l_cosmo_art, ldebug_art, l_pollen
!  New NL switch itheta_adv (introduced option for theta advection)
!  Renamed switch itype_lbcqx to more meaningful itype_outflow_qrsg 
!  New NL switch itype_lbc_qrsg for lateral boundary treatment of qr, qs, qg
! V4_10        2009/09/11 Matthias Raschendorfer, Jan-Peter Schulz
!  Introduction of LOGICALs limpltkediff, ltkesso and removing lturhor.
!  Introduction of INTEGER itype_sher
!  Introduction of LOGICAL lseaice for sea-ice model
! V4_11        2009/11/30 Ekaterina Machulskaya, Juergen Helmert, Guenther Zaengl
!  Introduced switch lmulti_snow to run multi-layer snow model (EM)
!  Introduced switches: itype_aeorsol, itype_root, itype_heatcond, itype_hydbound
!     and lemiss, lstomata for use of additional external parameter fields
!  Introduced itype_lbc_w (for type of lower boundary condition for w)
!     and ltadv_limiter (to use a limiter for temperature advection; but only
!     in case of itheta_adv=2) (GZ)
! V4_12        2010/05/11 Ulrich Schaettler, Michael Baldauf, Oli Fuhrer
!  Renamed itype_lbc_w itype_bbc_w because of "bottom boundary condition"
!  New flag l_dzeta_d_needed
!  Eliminated lhdiff_mask; introduced additional switches for treating
!  horizontal diffusion different in the interior and the boundary zones
!  and also different for special variables
! V4_13        2010/05/11 Michael Gertz
!  Adaptions to SVN
! V4_15        2010/11/19 Michael Baldauf
!  Replace Namelist Variable lva_impl_dyn by y_vert_adv_dyn
! V4_17        2011/02/24 Ulrich Blahak
!  Introduced lperi_x and lperi_y; eliminated lperi and lfreeslip_sfc.
! V4_18        2011/05/26 Michael Baldauf
!  Introduced new NL switch y_scalar_advect (which replaces lsl_adv_qx and yef_adv_qx)
! V4_20        2011/08/31 Matthias Raschendorfer
!  Introduction of  the switch 'ltkecon' for considering the TKE source from 
!  convection
! V4_21        2011/12/06 Michael Baldauf
!  Introduced new Namelist parameter l_diff_Smag for Smagorinsky diffusion
! V4_23        2012/05/10 Ulrich Schaettler, Oli Fuhrer, CLM
!  Removed switches lprogprec, ltrans_prec: these are eliminated from the code
!  Removed src_2timelevel and modified comments accordingly
!  Editorial Changes (by Oli Fuhrer)
!  New Namelist parameters iy_co2_stab and lco2_stab                       (CLM)
!  New Namelist parameter itype_albedo (for presribed surface albedo)      (CLM)
!  New Namelist parameter nincsn for calling frequency of spectral nudging (CLM)
! V4_24        2012/06/22 Michael Baldauf, Hendrik Reich
!  Introduced new Namelist variable itype_fast_waves (MB)
!  Enlarged number of digits for date variables yakdat1, yakdat2 to include
!    minutes and seconds (HR)
! V4_25        2012/09/28 Ulrich Blahak, Carlos Osuna
!  Implemented two internal switches:
!   - l_2mom:  to indicate that the 2-moment scheme is running
!              (users control of 2-moment scheme is by itype_gscp)
!   - l2mom_satads to switch on extra saturation adjustments outside the microphysics.
!  Introduce indices for current output step, used in logic of asynchronous I/O
!    netcdf (OS)
! V4_26        2012/12/06 Anne Roches
!  Replacement of hd_corr_q_XXX by hd_corr_trcr_XXX in order to be consistent
!  also with the naming of other switches (e.g. ltrcr_trilin, lef_adv_trcr_notpd)
! V4_27        2013/03/19 Michael Baldauf, Astrid Kerkweg, Ulrich Schaettler
!  New Namelist parameter divdamp_slope for new fast-waves solver
!  MESSy interface introduced (AK)
! V4_28        2013/07/12 Ulrich Schaettler
!  Added new Namelist parameter lroutine (to indicate an operational job and to
!   set GRIB2 meta data correct)
! V5_1         2014-11-28 Michael Baldauf, Matthias Raschendorfer, Ulrich Blahak
!                         Oliver Fuhrer, Anne Roches, Xavier Lapillonne, Lucio Torrisi
!  Introduced new namelist parameters:
!   - lhor_pgrad_Mahrer, l_3D_div_damping:             (MB)
!   - lscm (to activate single column model)           (MR)
!   - luse_radarfwo to activate radar forward operator (UB)
!  Removed Namelist-Parameter lexpl_lbc; Replaced crltau by crltau_inv (MB)
!  Replaced ireals by wp (working precision) (OF)
!  Introduced new namelist parameter ltraj, to run the online trajectory module (AR)
!  Introduced new namelist parameter nproma and variables nlastproma, nblock (XL)
!  Introduced new namelist parameter lsuper_coolw for switching on/off super-cooled
!    liquid water in microphysics (US)
!  Added new namelist parameters lsppt (LT)
!  Section 12 added for controlling the random number field generation for
!  the stochastic perturbation of physics tendencies scheme (SPPT) (LT)
!  Introduced new switch for 2-moment scheme: iradpar_cloud (KIT)
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
    wp,        & ! KIND-type parameter for real variables
    iintegers    ! KIND-type parameter for standard integer variables

!==============================================================================

IMPLICIT NONE

!==============================================================================

! Global (i.e. public) Declarations:

! 1. start and end of the forecast
! --------------------------------

  INTEGER   (KIND=iintegers), TARGET       ::           &
    nstart,       & ! first time step of the forecast
    nstop,        & ! last time step of the forecast period
    nfinalstop,   & ! last time step of the total forecast
                    ! (necessary, if simulation is splitted into periods)
    ntstep,       & ! actual time step
                    ! indices for permutation of three time levels
    nold,         & ! corresponds to ntstep - 1
    nnow,         & ! corresponds to ntstep
    nnew,         & ! corresponds to ntstep + 1
                    ! indices for permutation of two time levels
    ntke            ! corresponds to ntstep

! 2. boundary definition and update
! ---------------------------------

  INTEGER   (KIND=iintegers)       ::           &
    nlastbound,   & ! time step of the last boundary update
    nnextbound,   & ! time step of the next boundary update
    nincbound,    & ! time step increment of boundary update
    ndiff_ini_bd, & ! difference between start date and date of boundary data
    nbd1,         & ! indices for permutation of the 
    nbd2,         & ! two boundary time levels
    newbc,        & ! number of times that boundary update is analysis after 1 h
    newbcdt,      & ! time step increment of boundary update being derived from
                    ! the (latest) analysis (rather than forecast) fields
    nincboufac      ! factor to 'nincbound' when new boundary update is analysis

  REAL      (KIND=wp)              ::           &
    hstart,       & ! start of the forecast in full hours
    hstop,        & ! end of the forecast in hours
    hlastbound,   & ! last hour with boundary update
    hincbound,    & ! hour increment for reading boundary update
    hnextbound      ! next hour with boundary update

! 3. controlling the physics
! --------------------------

  INTEGER   (KIND=iintegers)       ::           &
    nincrad,      & ! time step increment for running the radiation
    nextrad,      & ! next step for running the radiation
    nradcoarse,   & ! number of horiz. gridpoints for radiation on coarser grid
    ninctura,     & ! time step increment for running the vertical diffusion
    nincconv,     & ! time step increment for running the convection scheme 
    nincsso,      & ! time step increment for running the SSO scheme

    itype_trvg,   & ! type of vegetation transpiration parameterization
    itype_evsl,   & ! type of parameterization of bare soil evaporation

    itype_gscp,   & ! type of grid-scale precipitaiton physics

    itype_wcld,   & ! type of water cloud diagnosis
    itype_tran,   & ! type of surface-atmosphere transfer
    itype_turb,   & ! type of turbulent diffusion parametrization
    itype_synd,   & ! type of diagnosis of synop. station values

    imode_tran,   & ! mode of surface-atmosphere transfer
    imode_turb,   & ! mode of turbulent diffusion parametrization

    ico2_rad,     & ! type of CO2 concentration in radiation parametrization
    iy_co2_stab,  & ! year when CO2 gets constant
#ifdef TWOMOM_SB
    iradpar_cloud,& ! method of calculating cloud optical properties >2 depends on effective radii
#endif
    icldm_rad,    & ! mode of cloud representation in radiation parametr.
    icldm_tran,   & ! mode of cloud representation in transfer parametr.
    icldm_turb,   & ! mode of cloud representation in turbulence parametr.

    itype_conv,   & ! type of convection parameterization:
                    !  0:  Tiedtke scheme (Default)
                    !  1:  Kain-Fritsch scheme
                    !  2:  Bechtold scheme (not yet implemented)
                    !  3:  Shallow convection based on Tiedtke scheme
    itype_aerosol,& ! type of aerosol map
                    !  1:  Tanre (Default)
                    !  2:  Tegen (exp)
    itype_root   ,& ! type of root distribution
                    ! 1: uniform (Default)
                    ! 2: exponential (exp)
    itype_heatcond,&! type of soil heat conductivity
                    ! 1: use average soil moisture
                    ! 2: take into account soil moisture/soil ice
    itype_hydbound,&! type of hydraulic lower boundary
                    ! 1: allow for drainage but not diffusion
                    ! 2: rigid lid
                    ! 3: ground water with drainage and diffusion
    itype_sher,   & ! type of shear production for TKE
                    !  1:  only vertical shear
                    !  2:  full isotropic 3D-shear
                    !  3:  vertical shear and separted horizontal shear mode

    itype_albedo, & ! type of surface albedo treatment
                    ! 1 (default) surface albedo is a function of soiltype
                    ! 2 surface albedo is prescribed by external fields
                    ! 3 background albedo is prescribed by external fields
                    ! 4 vegetation albedo is modified by forest fraction

    nlgw_ini,     & ! number of prognostic soil water levels in initial data
    nlgw_bd,      & ! number of prognostic soil water levels in boundary data
    nlgw,         & ! number of prognostic soil water levels

    ibot_w_so,    & ! number of hydrological active soil layers
    nhori,        & ! number of sectors for the horizont array by the topographic
                    ! correction of the radiation
    nproma,       & ! block size for physical parameterizations
    nlastproma,   & ! size of last block
    nblock          ! number of blocks

  LOGICAL                          ::           &
    lphys,        & ! forecast with physical parametrizations
    lrad,         & ! forecast with radiation
    lforest,      & ! if .true., run with forest (evergreen and deciduous)
    lsso,         & ! process parameters for sso scheme
    lemiss,       & ! external surface emissivity map
    lstomata,     & ! minimum stomata resistance map
    ltur,         & ! forecast with vertical diffusion
    l3dturb,      & ! 3D-turbulence: CALL explicit_horizontal_diffusion (RK)
    l3dturb_metr, & ! switch on/off additional metric terms for 3D-turbulence
    l_dzeta_d_needed = .FALSE., &
                    ! metric coeff. dzeta_dlam, dzeta_dphi are needed
    lprog_tke,    & ! prognostic treatment of TKE (for itype_turb=5/7)
    limpltkediff, & ! use semi-implicit TKE diffusion
    lconv,        & ! forecast with convection
    lconv_inst,   & ! output of instantaneous values of top_con/bas_con
                    ! instead of min/max for an output interval
    lgsp,         & ! forecast with grid scale precipitation
    lsuper_coolw, & ! switch for supercooled liquid water (work from Felix Rieper)
    ldiniprec,    & ! diagnostic initialisation of prognostic precip (qr, qs)
    lsoil,        & ! forecast with soil model
    lmelt,        & ! soil model with melting process
    lmelt_var,    & ! freezing temperature dependent on water content
    lmulti_layer, & ! run multi-layer soil model
    lmulti_snow,  & ! run multi-layer snow model
    lseaice,      & ! forecast with sea ice model
    llake,        & ! forecast with lake model
    lconf_avg,    & ! average convective forcings in case of massflux closure
    lradf_avg,    & ! average radiative forcings if radiation is calculated on coarser grid
    lcape,        & ! convection with CAPE closure
    lctke,        & ! convection with turbulent convective energy closure
                    ! warning: lctke not yet fully implemented
    ltkesso,      & ! calculation SSO-wake turbulence production for TKE
    ltkecon,      & ! consider convective buoyancy production for TKE
    lexpcor,      & ! explicit corrections of the implicit calculated
                    ! turbulent diffusion (only if itype_turb=3)
    ltmpcor,      & ! consideration of thermal TKE-sources in the 
                    ! enthalpy budget
    lprfcor,      & ! using the profile values of the lowest main level instead
                    ! of the mean value of the lowest layer for surface flux
                    ! calulations
    lnonloc,      & ! nonlocal calculation of vertical gradients used
                    ! for turbulent diffusion (only if itype_turb=3)
#ifdef tkvfilter
    ltkvfilter,   & ! vertical filtering of tkvh/tkvm 
#endif
    lcpfluc,      & ! consideration of fluctuations of the heat capacity of air
    lco2_stab,    & ! enable CO2 stabilisation
    lradtopo        ! if .TRUE., calculate topographic correction of radiation

  REAL (KIND=wp)     ::           &
    hincrad,      & ! increment for running the radiation in hours
    hnextrad,     & ! next step for running the radiation in hours
    czbot_w_so      ! depth of bottom of last hydrological active soil layer [m]

! 4. controlling the dynamics
! ---------------------------

  LOGICAL                          ::           &
    l2tls,             & ! time integration by two timelevel RK-scheme (.TRUE.)
                         ! or by default three-time level KW-scheme (.FALSE.)
    ltadv_limiter,     & ! use limiter for temperature advection (itheta_adv=2 only)
    lsemi_imp,         & ! if .TRUE.,  running with semi-implicit scheme,
                         ! else with split-explicit scheme (only for l2tls=FALSE!)
    ldyn_bbc,          & ! dynamical bottom boundary condition
    lcori_deep,        & ! if =.TRUE.: take cos(phi) coriolis terms into account
    ladv_deep,         & ! if =.TRUE.: use all metric advective terms
    lhor_pgrad_Mahrer, & ! if =.TRUE., horizontal p-gradients are calculated 
                         ! analogous to Mahrer (1984) (only for itype_fast_waves==2)
    l_3D_div_damping     ! if =.TRUE., the fully 3D (=isotropic) divergence damping
                         ! is used (only for itype_fast_waves==2)

  CHARACTER (LEN=10)               ::           &
    y_vert_adv_dyn       ! ="expl": explicit vertical advection ...
                         ! ="impl2": current implicit vertical advection in RK-scheme
                         ! ="impl3": new implicit method outside of the RK-scheme

  CHARACTER (LEN=20)               ::           &
    y_scalar_advect      ! type of scalar advection scheme
                         ! "SL3_MF", "SL3_SFD", "Bott2", "Bott4"
                         ! "Bott2_Strang", "Bott4_Strang", "vanLeer", "PPM"

  INTEGER   (KIND=iintegers)       ::           &
    nbl_exchg,         & ! number of boundlines to exchange: because this varies
                         ! especially for the Runge-Kutta scheme, it is set by the
                         ! program, depending on iadv_order, intcr_max
    irunge_kutta,      & ! =1: use new RK scheme
                         ! =2: use new TVD-RK scheme
    irk_order,         & ! order of the Runge-Kutta scheme
    iadv_order,        & ! order of the horizontal advection scheme for dynamics
    ieva_order,        & ! order of the explicit vertical advection scheme
    itheta_adv,        & ! =0: use T' (perturbation temperature) for advection
                         ! =1: use theta' (perturbation potential temperature)
                         ! =2: use theta (full potential temperature)
    itype_fast_waves,  & ! Type of fast waves solver for Runge-Kutta dynamics
                         ! = 1: old scheme (from module fast_waves_rk.f90)
                         ! = 2: new scheme (from module fast_waves_sc.f90)
    itype_bbc_w,       & ! Bottom boundary condition for vertical wind
                         ! =0/1: RK-like method following iadv_order
                         ! =2/3: differencing following iadv_order without RK stepping
                         ! =4/5: Fourth-order centered differences
                         ! 0/2/4: include linear extrapolation of horizontal wind to sfc
                         ! 1/3/5: no extrapolation of horizontal wind to sfc
    intcr_max,         & ! max. integer courant number in cr-independent advection
    itype_outflow_qrsg,& ! type of relaxation treatment for qr, qs, qg
                         ! =1: (default) same treatment as all variables
                         ! =2: no relaxation for qr,qs,qg at outflow boundary points
    itype_lbc_qrsg,    & ! type of lateral boundary treatment for qr, qs, qg, 
                         ! =1: (default) zero-gradient condition 
                         ! =2: set qr,qs,qg to 0.0 at the boundaries
                         ! =3: no presetting at the boundaries 
                         !     (must not be chosen for Leapfrog applications)
    itype_spubc,       & ! type of Rayleigh damping in the upper levels
    nfi_spubc2,        & ! Number of applications of smoother for the determination
                         !  of the large scale field used in the Rayleigh damping 
                         !  with itype_spubc=2
    ikrylow_si,        & ! dimension of the Krylow space used in the elliptic
                         ! solver for the semi-implicit scheme
    maxit_si,          & ! maximum number of iterations for the elliptic solver
    iprint_si            ! to control whether statistics of the solver are printed

  REAL      (KIND=wp)              ::           &
    xkd,               & ! coefficient for divergence damping
    divdamp_slope,     & ! exceed the theoretical slope stability criterion of
                         ! the divergence damping (only for itype_fast_waves=2)
    eps_si               ! precision limit for the elliptic solver

! 5. controlling the observation processing
! -----------------------------------------

  LOGICAL                          ::           &
    luseobs         ! on - off switch for using observational data for:
                    ! - nudging (of conventional data)
                    ! - latent heat nudging (not implemented yet)
                    ! - 2-dim. analyses (2m-Temperature, 2m-Humidity, precipit.)
                    ! - verification of model data against observations

! 6. controlling the upper and lateral boundary conditions
! --------------------------------------------------------

  LOGICAL                          ::           &
    lspubc,       & ! with Rayleigh damping in the upper levels
    lrubc           ! radiative upper boundary condition

  REAL      (KIND=wp)              ::           &
    rdheight,     & ! bottom height of Rayleigh damping layer
    crltau_inv,   & ! factor for relaxation time 1/tau_r = crltau_inv * 1/dt
    rlwidth,      & ! width of relaxation layer
    relax_fac       ! reduction factor for strength of lateral boundary relaxation
                    ! (relevant for radiative lateral boundary conditions)

  INTEGER   (KIND=iintegers)       ::           &
    nrdtau          ! number of time steps in Rayleigh damping time scale     

! 7. additional control variables
! -------------------------------

  LOGICAL                          ::           &
    ltraj,        & ! if .TRUE., compute trajectories
    l_cosmo_art,  & ! if .TRUE., run the COSMO_ART
    l_pollen,     & ! if .TRUE., run the Pollen component
    lroutine,     & ! if .TRUE., run an operational forecast
    llm,          & ! if .TRUE., running with a lowered upper boundary
    lprog_qi,     & ! if .TRUE., running with cloud ice
                    !   (this is set internally by the program,
                    !    depending on itype_gscp)
    lcond,        & ! forecast with condensation/evaporation
    ldiabf_lh,    & ! include diabatic forcing due to latent heat in RK-scheme
    ldiabf_satad, & ! include diabatic forcing due to saturation adjustment
    l_2mom,       & ! to indicate, whether the 2-moment scheme is running
    l2mom_satads, & ! in case of 2-moment scheme, do all the satads
                    ! (like for the 1-moment schemes), not just the
                    ! satad after the microphysics at the end of the timestep.
    lclock,       & ! system clock is present
    ltime,        & ! detailed timings of the program are given
    lreproduce,   & ! the results are reproducible in parallel mode
    lhordiff,     & ! running with horizontal diffusion
    ldump_ascii,  & ! for flushing (close and re-open) the ASCII files
    lrerun,       & ! not implemented
    lrout,        & ! routine-forecast of the model
    lartif_data,  & ! forecast with self-defined artificial data
    lperi_x,      & ! lartif_data=.TRUE.:  periodic boundary conditions in x-dir.
                    !            =.FALSE.: with Davies conditions
    lperi_y,      & ! lartif_data=.TRUE.:  periodic boundary conditions in y-dir.
                    !            =.FALSE.: with Davies conditions
    l2dim,        & ! lartif_data=.TRUE.:  2dimensional model version
                    !            =.FALSE.: full 3dimensional version
    lcori,        & ! lartif_data=.TRUE.:  with Coriolis force
                    !            =.FALSE.: or without Coriolis force
    lmetr,        & ! lartif_data=.TRUE.:  with metric terms
                    !            =.FALSE.: or without metric terms
    lradlbc,      & ! lartif_data=.TRUE.:  radiative lateral boundary conditions
                    !            =.FALSE.: or with Davies conditions
    lw_freeslip,  & ! if .TRUE.: with free slip lateral boundary condition and
                    ! if .FALSE. specified lateral boundary values for w
    ldfi,         & ! switch for initialization by digital filtering
                    ! (if .TRUE. do apply dfi)
    luse_rttov,   & ! if rttov-library is used
    lscm,         & ! if Single Column Model is used (default: FALSE)
    luse_radarfwo   ! if the radar forward operator is switched on

  REAL      (KIND=wp)              ::           &
    hlastmxu,     & ! last hour when vbmax was "nullified"
    hnextmxu,     & ! next hour when vbmax will be "nullified"
    hincmxu,      & ! increment that can be specified via Namelist
    hlastmxt,     & ! last hour when tmin, tmax were "nullified"
    hnextmxt,     & ! next hour when tmin, tmax will be "nullified"
    hincmxt         ! increment that can be specified via Namelist

  INTEGER   (KIND=iintegers)       ::           &
    nlastmxu,     & ! last step when vbmax was "nullified"
    nnextmxu,     & ! next step when vbmax will be "nullified"
    nlastmxt,     & ! last step when tmin, tmax were "nullified"
    nnextmxt,     & ! next step when tmin, tmax will be "nullified"
    nvers,        & ! version number of experiment for documentation
    isynsat_stat, & ! for global status of synthetic satellite images
    itype_hdiff     ! type of horizontal diffusion (=1: 4th order linear),
                    ! =2: 4th order linear monotonic with orographic limiter)

  REAL      (KIND=wp)              ::           &
    ! Correction factor for horizontal diffusion fluxes of:
    hd_corr_u_bd,    & ! u,v,w         in boundary zone
    hd_corr_t_bd,    & ! t             in boundary zone
    hd_corr_trcr_bd, & ! tracers       in boundary zone
    hd_corr_p_bd,    & ! p             in boundary zone
    hd_corr_u_in,    & ! u,v,w         in domain
    hd_corr_t_in,    & ! t             in domain
    hd_corr_trcr_in, & ! tracers       in domain
    hd_corr_p_in,    & ! p             in domain
    hd_dhmax           ! maximum gridpoint height difference for applying
                       ! horizontal diffusion fluxes between them

  LOGICAL                          ::           & 
    l_diff_Smag     ! use Smagorinsky-diffusion for u and v

! 8. diagnostic calculations
! --------------------------

  LOGICAL                          ::           &
    ldiagnos           ! perform diagnostic calculations

  INTEGER   (KIND=iintegers)       ::           &
    itype_diag_t2m,  & ! type of T_2M diagnostics
    itype_diag_gusts   ! type of gusts diagnostics

  ! for meanvalues, that have to be saved for restart runs
#ifndef MESSY
  REAL    (KIND=wp)                      ::             &
#else
  REAL    (KIND=wp),     POINTER         ::             &
#endif
    psm0,   & ! initial value for mean surface pressure ps
    dsem0,  & ! initial value for mean dry static energy
    msem0,  & ! initial value for mean moist static energy
    kem0,   & ! initial value for mean kinetic energy
    qcm0      ! initial value for mean cloudwater content

! 9. Variables for Ascii file handling, time measuring, ...
! ---------------------------------------------------------

  INTEGER   (KIND=iintegers)       ::           &
    itype_calendar   ! for specifying the calendar used
                     !  = 0: gregorian calendar (default)
                     !    (but this needs a bug fix in get_utc_date,
                     !    because up to now we only have julian calendar)
                     !  = 1: every year has 360 days
                     !  = 2: every year has 365 days

  CHARACTER (LEN=14)               ::           &
    yakdat1   ! actual date (ydate_ini+ntstep/dt) in the form
              ! yyyymmddhhmmss (year, month, day, hour, minutes, seconds)
  CHARACTER (LEN=28)               ::           &
    yakdat2   ! actual date (ydate_ini+ntstep/dt) in the form
              ! wd dd.mm.yyyy  hh mm ss UTC  (weekday, ...)

  INTEGER   (KIND=iintegers)       ::           &
    nusolver,     & ! unit number for file YUSOLVER
    nudebug,      & ! unit number for file YUDEBUG
    nuspecif        ! unit number for file YUSPECIF

  INTEGER   (KIND=iintegers)       ::           &
    itype_timing    ! determines, how to handle the timing
                    !  1,2: timings per processor, hourly (1) or summed up (2)
                    !  3,4: mean values for timings, hourly (3) or summed up (4)

  CHARACTER (LEN= 8) :: yusolver='YUSOLVER'
  CHARACTER (LEN= 7) :: yudebug ='YUDEBUG'
  CHARACTER (LEN= 8) :: yuspecif='YUSPECIF'

  LOGICAL                          ::           &
    lout_anai       ! allocate fields to enable writing analysis increments

! 10. Variables for Spectral Nudging
! ----------------------------------
! This feature has been introduced by the Climate-LM Community.
! (Implemented by GKSS: Rockel)

! Some of the structures defined here are also defined in data_io.
! But we did not want to have dependencies to the I/O here

  INTEGER (KIND=iintegers), PARAMETER :: &
    nzmxsn =  30    ! maximum number of variables for spectral nudging

  CHARACTER (LEN= 10) :: yvarsn (nzmxsn) ! list of fields for spectral nudging
                                         ! (must be a subset of yvarbd)

  INTEGER   (KIND=iintegers)       ::           &
    nyvar_sn,     & ! number of variables for spectral nudging
                    !  (must be <= nyvar_b)
    isc_sn,       & ! spectral nudging in i-direction
    jsc_sn,       & ! spectral nudging in j-direction
    nincsn          ! time step increment for calling spectral nudging

  REAL   (KIND=wp)           ::           &
    pp_sn,        & ! lowest pressure level for spectral nudging
    alpha_sn        ! amplification factor for spectral nudging
                    ! (0. <= alpha_sn <= 1.)

  LOGICAL                          ::           &
    lspecnudge      ! spectral nudging of boundary data  ! GKSS (Rockel)

  ! Type for the list descriptions
  TYPE list_description
    CHARACTER (LEN=10)             :: name                ! name of variable
    INTEGER (KIND=iintegers)       :: iloc1, iloc2, iloc3 ! location in vartab
    INTEGER (KIND=iintegers)       :: idimvert            ! vertical dimension
  END TYPE

  ! Variables for describing the different lists
  TYPE (list_description), ALLOCATABLE :: list_sn(:)

! 11. controlling ensemble mode (EPS) and stochastic physics (SPPT)
! -----------------------------------------------------------------

  LOGICAL                          ::           &
    leps             ! switch ensemble mode on/off

  INTEGER   (KIND=iintegers)       ::           &
    iepsmem,              & ! ID of the member in the ensemble (ID >= 0)
    iepstot,              & ! total number of ensemble members (>=0)
    iepstyp                 ! ID of the ensemble generation type (ID >= 0)

  REAL      (KIND=wp)              ::           &
    fac_plcov,            & ! modification factor for PLCOV
    rmin_plcov,           & ! lower limit of PLCOV
    rmax_plcov,           & ! upper limit of PLCOV
    fac_rootdp,           & ! modification factor for ROOTDP
    rmin_rootdp,          & ! lower limit of ROOTDP
    rmax_rootdp,          & ! upper limit of ROOTDP
    fac_lai,              & ! modification factor for LAI
    rmin_lai,             & ! lower limit of LAI
    rmax_lai                ! upper limit of LAI


! 12. controlling random number generation for SPPT (stochastic physics)
! -------------------------------------------------

  LOGICAL                          ::           &
    lsppt,              & ! switch, if .true., perturb the physical tendencies
    lhorint_rn,         & ! random numbers (defined on a rn horiz. coarse grid)
                          !   horizontally interpolated on model grid (otherwise
                          !   model grid points contained in the same rn coarse
                          !   grid point have the same random number value)
    ltimeint_rn,        & ! random numbers (defined on a rn horiz. coarse grid)
                          !   are interpolated in time
    lgauss_rn             ! use a gaussian distribution of random numbers
                          !   (otherwise a uniform distribution is used)

  INTEGER   (KIND=iintegers), PARAMETER  ::     &
    npatmax = 5           ! max. number of random number patterns (with
                          !   different correlation lengths (horiz'ly, in time)

  INTEGER   (KIND=iintegers)       ::           &
    npattern_rn,        & ! number of rn patterns with different scale/correl.
!   ninc_rn  (npatmax), & ! timestep increment for drawing a new field of
                          !   random numbers (used if (imode_rn == 0))
    nseed_rn (npatmax), & ! external part of seed for random number generation
    nseed_rn2(npatmax), & ! ext. seed to generate rn at ntstep>0 (if imode_rn>0)
    imode_rn,           & ! = 0: use only 1 stream of rn for all rn time steps
                          ! = 1: use a new stream of rn for every rn time step,
                          !      this enables temporal corrations in DA cycles
    itype_vtaper_rn,    & ! type of vertical tapering near surface and/or
                          !   in stratosphere
    itype_qxpert_rn,    & ! define which hum variables tend. are perturbed 
    itype_qxlim_rn,     & ! type of reduction/removal of the perturbation 
                          !   in case of negative (qv, qc, qi) or 
                          !   supersaturated (qv) values
    n1_rn    (npatmax), & ! / indices for permutation of the
    n2_rn    (npatmax), & ! \ two random number time levels
    ie_rn    (npatmax), & ! number of horiz. coarse grid points in zonal direct.
    je_rn    (npatmax)    ! number of horiz. coarse grid points in merid. dir.
                          !   where random numbers are defined

  REAL      (KIND=wp)              ::           &
    hinc_rn  (npatmax), & ! time increment (in [hrs]) for drawing a new field
                          !   of random numbers (used if (imode_rn > 0))
    hstep_rn (npatmax), & ! time [hrs] for which the latest rn field is valid
    dlat_rn  (npatmax), & ! random number coarse grid point distance in
                          !   meridional direction (in degrees)
    dlon_rn  (npatmax), & ! random number coarse grid point distance in
                          !   zonal direction (in degrees)
    stdv_rn  (npatmax), & ! standard deviation of the gaussian distribution of
                          !   random numbers
    range_rn (npatmax)    ! max magnitude of random numbers

  REAL (KIND=wp),     ALLOCATABLE          ::           &
    rvtaper_rn (:)        ! externally specified function for vertical tapering
                          !   of the random numbers

! 13. controlling verbosity of debug output and miscellaneous
! -----------------------------------------------------------

  INTEGER   (KIND=iintegers)       ::           &
    idbg_level      ! to control the verbosity of debug output
                    ! Some basic output is written anyhow
                    ! but for most components the output must be
                    ! activated with the logical switches below

  ! the following switches can activate additional debug output for
  ! special components (in combination with idbg_level)
  LOGICAL                          ::           &
    ldebug_dyn,   & ! if .TRUE., debug output for dynamics
    ldebug_gsp,   & ! if .TRUE., debug output for grid scale precipitation
    ldebug_rad,   & ! if .TRUE., debug output for radiation
    ldebug_tur,   & ! if .TRUE., debug output for turbulence
    ldebug_con,   & ! if .TRUE., debug output for convection
    ldebug_soi,   & ! if .TRUE., debug output for soil model
    ldebug_io ,   & ! if .TRUE., debug output for I/O
    ldebug_mpe,   & ! if .TRUE., debug output for mpe_io
    ldebug_dia,   & ! if .TRUE., debug output for diagnostics
    ldebug_art,   & ! if .TRUE., debug output for COSMO_ART
    ldebug_ass,   & ! if .TRUE., debug output for assimilation
    ldebug_lhn      ! if .TRUE., debug output for latent heat nudging

  LOGICAL                          ::           &
    lprintdeb_all,& ! .TRUE.:  all tasks print debug output
                    ! .FALSE.: only task 0 prints debug output
    linit_fields    ! to initialize also local variables with a default value
                    ! (some compilers seem to have problems, if fields are not
                    !  initialized, even if values are never used;
                    !  also, debugging could be easier then)

! 14. controlling asynchronous output
! -----------------------------------

INTEGER (KIND=iintegers) ::  &
    cur_outstep,        &   ! current output time step
    cur_outstep_idx,    &   ! index of current output time step
    cur_gribout_idx         ! index of current gribout section

!==============================================================================

END MODULE data_runcontrol
