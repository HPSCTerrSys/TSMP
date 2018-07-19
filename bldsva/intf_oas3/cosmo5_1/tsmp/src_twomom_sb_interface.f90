!+ Module for interface routine for the Seifert-Beheng twomoment microphysics
!------------------------------------------------------------------------------

MODULE src_twomom_sb_interface

!------------------------------------------------------------------------------
!
! Description:
!
! This module contains the interface subroutine "seifert_pp" for the
! Seifert/Beheng two-moment bulk cloud microphysical scheme.
!
! The routines for the scheme itself can be found in the file
! "src_twomom_sb.f90", which you can obtain via contacting
! the current code owner (see below).
!
! This interface is necessary because the original scheme has
! been developed originally for the Karlsruhe Atmospheric Mesoscale Model
! Verison 2 (KAMM2), which used different conventions for the
! hydrometeor tracers (mass- and number densities instead of 
! mass specific quantities).  
!
! The scheme has been originally developed by Axel Seifert and
! Klaus Beheng at the Institute for Meteorology and Climate Reserach,
! Karlsruhe Institute of Technology (KIT), Germany.
! In its original version, it comprised the 5 hydrometeor species
! cloud droplets, rain, cloud ice, snow, and graupel.
! Ulrich Blahak and Heike Noppel (both also at KIT at that time)
! extended the scheme by an additional hail class to enabble
! more realistic simulations of hail processes in deep convective
! clouds at high spatial resolution.
!
! More comments on using the scheme can be found in the header
! information of subroutine "seifert_pp", in particular how to set
! the namelist parameter "itype_gscp" accordingly to configure
! the scheme to one's needs.
!
! NOTE: ice nucleation needs a clean-up! At the moment,
!       one can either use the new scheme based on works of
!       V. Phillips/B. Kaercher/U. Lohmann with its 9 suboptions
!       or the old scheme as originally implemented by A. Seifert
!       with also its 9 suboptions (comprising most of the "classical"
!       ice nucleation schemes)
!
!
! Preprocessor flags necessary for the Seifert/Beheng 2-moment scheme 
!   (set in Makefile / Fopts):
!
! MPI_CLOUD:             Define it, if you want to enable the MPI-code
!                        necessary to use load balancing.
!                        If you then really want to use load balancing,
!                        set the below parameter
!                        LOGICAL, PARAMETER :: MPI_LOAD_BALANCING = .true.
!                        Load balancing speeds up computations on
!                        scalar many-core architectures.
!
! SEDI_VECTORIZED:       If defined, vectorizable versions of the sedimenation routines
!                        are used. You can find the code in src_seifert.f90. This is recommended
!                        on vector machines. On scalar machines, you can use
!                        both versions, and the original routines should be
!                        slightly faster.
!
! CGP_SEARCH_VEC:        If defined, use vectorizable version of cloudy 
!                        grid point search.
!                        Is more efficient on vector machines.
!                        On scalar machines, this switch is not necessary.
!                        
!-------------------------------------
!
! Current Code Owner: DWD, Ulrich Blahak
!  phone:  +49  69  8062 2393
!  fax:    +49  69  8062 3721
!  email:  ulrich.blahak@dwd.de
!
! History:
! Version    Date       Name
! ---------- ---------- ----
! 1.1        2012/07/25  Ulrich Blahak     
!  Initial release
! @VERSION@    @DATE@    Ulrich Blahak
!  Bugfix for l2mom_satads=.true.: cloud base w (w_cb) is now
!   determined using qc instead of ssw, because ssw is 0 if
!   the satad after the dynamics and before the microphysics is
!   made.
!
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


  ! COSMO modules

  USE data_parameters , ONLY :   &
       wp, dp, sp,   & ! KIND-type parameters for real variables
       iintegers       ! kind-type parameter for "normal" integer variables

  USE data_runcontrol,   ONLY :   &
       ntstep,       & ! actual time step
       nstart,       & ! first time step of the forecast
       nnew,         & ! corresponds to ntstep + 1
       itype_gscp,   & ! type of grid-scale precipitation physics
       l2tls,        & ! forecast with 2-TL integration scheme
       ldiabf_lh,    & ! include diabatic forcing due to latent heat in RK-scheme
       ltime, &
       l2mom_satads    ! in case of 2-moment scheme, do all the satads
                       ! (like for the 1-moment schemes), not just the
                       ! satad after the microphysics at the end of the timestep.
  
  USE data_modelconfig,  ONLY :   &
       ie,           & ! number of grid points in zonal direction
       je,           & ! number of grid points in meridional direction
       ke,           & ! number of grid points in vertical direction
       ie_tot,       & ! number of grid points in zonal direction
       je_tot,       & ! number of grid points in meridional direction
       istartpar,    & ! start index for computations in the parallel program
       iendpar,      & ! end index for computations in the parallel program
       jstartpar,    & ! start index for computations in the parallel program
       jendpar,      & ! end index for computations in the parallel program
       dt,           & ! timestep
       dt2             ! 2 * dt

  USE data_constants  , ONLY :   &
       pi,           & ! 

 ! 2. physical constants and related variables
 ! -------------------------------------------
       t0_melt,      & ! melting temperature of ice    
       r_d,          & ! gas constant for dry air
       r_v,          & ! gas constant for water vapour
       rdv,          & ! r_d / r_v
       rcpv,         & ! cp_d / cp_v - 1
       rcpl,         & ! cp_d / cp_l - 1
       o_m_rdv,      & ! 1 - r_d/r_v
       rvd_m_o,      & ! r_v/r_d - 1
       cp_d,         & ! specific heat of dry air
       cpdr,         & ! 1 / cp_d
       lh_v,         & ! latent heat of vapourization
       lh_f,         & ! latent heat of fusion
       lh_s,         & ! latent heat of sublimation
       g,            & ! acceleration due to gravity

  ! 3. constants for parametrizations
  ! ---------------------------------
       b1,           & ! variables for computing the saturation vapour pressure
       b2w,          & ! over water (w) and ice (i)
       b2i,          & !               -- " --
       b3,           & !               -- " --
       b4w,          & !               -- " --
       b234w,        & !               -- " --
       b4i             !               -- " --

  USE data_fields,       ONLY :   &
       hhl        ,    & ! height of model half levels                   (m )
       rho0       ,    & ! reference density at the full model levels    (kg/m3)
       p0         ,    & ! reference pressure at full levels             (Pa)
       w          ,    & ! vertical wind speed (defined on half levels)  ( m/s )
       t          ,    & ! temperature                                   (  k  )
       pp         ,    & ! deviation from the reference pressure         ( pa  )
       tinc_lh    ,    & ! temperature increment due to latent heat      (  K  )
       rho        ,    & ! total density of moist air                    (kg/m3)
       qrs        ,    & ! precipitation water (water loading)           (kg/kg )
       prr_gsp    ,    & ! precipitation rate of rain, grid-scale        (kg/m2s)
       prs_gsp    ,    & ! precipitation rate of snow, grid-scale        (kg/m2s)
       prg_gsp    ,    & ! precipitation rate of graupel, grid-scale     (kg/m2s)
       prh_gsp           ! precipitation rate of hail, grid-scale        (kg/m2s)

  USE data_parallel,      ONLY :  &
       num_compute,     & ! number of compute PEs
       nboundlines,     & ! number of boundary lines of the domain for which
                          ! no forecast is computed = overlapping boundary
       my_world_id,     & ! rank of this subdomain in the global communicator
       my_cart_id,      & ! rank of this subdomain in the cartesian communicator
       icomm_cart,      & ! communicator for the virtual cartesian topology
       imp_reals,       & ! determines the correct REAL type used in the model
                          ! for MPI
       imp_integers       ! determines the correct INTEGER type used in the model

  USE parallel_utilities, ONLY :   &
       global_values      ! collects max/min/sum values from all nodes

  USE meteo_utilities,          ONLY :  &
       calrho, satad

  USE pp_utilities, ONLY : gamma_fct

  USE time_utilities, ONLY:  get_timings, &
       i_2mom_start, i_2mom_init, i_2mom_todens, i_2mom_clouds, i_2mom_tospecif, &
       i_2mom_sedi,  i_2mom_cleanup, i_2mom_init_0, i_2mom_cgps

  USE environment,              ONLY : model_abort, get_free_unit, release_unit

  USE src_tracer,               ONLY : trcr_get, trcr_errorstr

#ifdef NUDGING
  USE data_lheat_nudge,           ONLY :  &
       llhn,         & ! main switch for latent heat nudging
       llhnverif,    & ! main switch for latent heat nudging
       lhn_qrs,      & ! use integrated precipitaion flux as reference
       tt_lheat,     & ! latent heat release of the model
       qrsflux         ! total precipitation flux
  
  USE src_lheating,             ONLY :  &
       get_gs_lheating            ! storage of grid scale latent heating for lhn
#endif

  ! 2MOM modules from src_twomom_sb.f90:

  USE wolken_driver,     ONLY: cloud_type => wolke_typ,          &
       &                       loc_ix, loc_iy, loc_iz,           &
       &                       dim_ix, dim_iy, dim_iz,           &
       &                       dt_2mom => dt,                    &
       &                       T_2mom => T,                      &
       &                       p_2mom => p,                      &
       &                       rho_2mom => rho,                  &
       &                       w_2mom => w,                      &
       &                       q,                                &
       &                       q_cloud, n_cloud,                 &
       &                       q_ice, q_rain, q_snow, q_graupel, &
       &                       n_ice, n_rain, n_snow, n_graupel, &
       &                       n_hail, q_hail,                   &
       &                       alloc_driver, dealloc_driver, cp, &
       &                       S_w, S_i,dSwdz, dSidz, dT0dz,     &
       &                       zml_k,                            &
       &                       w_cb,                             &
       &                       speichere_dqdt,       &
       &                       speichere_precipstat,             &
       &                       dqdt, cond_neu_sb, evap_neu_sb,   &
       &                       nrates, nmicrorates,nsedifluxdiv, &
       &                       dmin_wg_g, pvec_wg_g, Tvec_wg_g,  &
       &                       qwvec_wg_g, qivec_wg_g, &
       &                       anzp_wg, anzT_wg, anzi_wg, anzw_wg, &
       &                       ltabdminwgg


#ifdef COSMOART
  USE wolken_driver,     ONLY: calc_i,calc_j,calc_k,cloudact
#endif

  USE parallele_umgebung, ONLY:                                  &
       &                       isIO,                    &
       &                       global_maxval,global_minval,      &
       &                       global_maxval_2d
  USE wolken,            ONLY: alloc_wolken, dealloc_wolken,     &
       &                       clouds,                           &
       &                       clip_and_limit_ireals,            &
       &                       clip_only_ireals

  USE wolken_sedi,       ONLY:                                   &
       &                       rain_sedimentation_lm,            &
       &                       ice_sedimentation_lm,             &
       &                       snow_sedimentation_lm,            &
       &                       graupel_sedimentation_lm,         &
       &                       hail_sedimentation_lm


  USE wolken_eis, ONLY : init_dmin_wetgrowth, &
                         init_dmin_wg_gr_ltab_equi, het_nuc_phillips

  USE wolken_konstanten, ONLY: init_seifert,                     &
!!$ we use now COSMO-constants>   L_ew,L_ed,R_d,R_l,                &
                               cloud,rain,ice,snow,graupel,hail, &
                               rho0_2mom => rho0,               &
                               e_ws,e_es,                &
                               e_ws_vec,e_es_vec,                &
                               T_nuc, T_f, satad_nach_mikrophysik, & 
                               graupel_shedding, hail_shedding, tshift_cloudfreeze,    &
                               qnc_const_seifert => qnc_const, nuc_c_typ, &
                               z_ncn_segal_seifert => z_ncn_segal, &
                               z1oe_ncn_segal_seifert => z1oe_ncn_segal




!==============================================================================

IMPLICIT NONE

!==============================================================================

! Look-up table for Phillips et al. nucleation
INTEGER(KIND=iintegers), PARAMETER :: &
  ttmax  = 30,      &  ! sets limit for temperature in look-up table
  ssmax  = 60          ! sets limit for ice supersaturation in look-up table
! UB_20100421>> von integer nach real(kind=wp) geaendert:
REAL(KIND=wp), PARAMETER :: &
  ttstep = 2.0_wp,       &  ! increment for temperature in look-up table
  ssstep = 1.0_wp           ! increment for ice supersaturation in look-up table

REAL(KIND=wp), DIMENSION(0:100,0:100), SAVE :: &
  afrac_dust, &  ! look-up table of activated fraction of dust particles acting as ice nuclei
  afrac_soot, &  ! ... of soot particles
  afrac_orga     ! ... of organic material

INCLUDE 'phillips_nucleation.incf'


!==============================================================================

!==============================================================================
! Module procedures
!==============================================================================

CONTAINS

!==============================================================================

!------------------------------------------------------------------------------
! Begin of function section
!------------------------------------------------------------------------------

!==============================================================================
!
! Some utility functions for the Phillips-nucleation:
!
!==============================================================================

REAL (KIND=wp) ELEMENTAL FUNCTION diff_vapor(temp,pres)
IMPLICIT NONE

REAL (KIND=wp), INTENT(IN)  :: temp, pres

REAL (KIND=wp), PARAMETER   :: a = 8.7602e-5_wp
REAL (KIND=wp), PARAMETER   :: b = 1.81_wp

diff_vapor = a * EXP( b * LOG(temp)) / pres

!diff_vapor = 0.211e-4 * (tg/t0)**1.94 * (1013.25e2/ppg)

END FUNCTION diff_vapor

!==============================================================================

ELEMENTAL FUNCTION sat_pres_water(temp)
IMPLICIT NONE

REAL (KIND=wp)              :: sat_pres_water
REAL (KIND=wp), INTENT(IN)  :: temp

sat_pres_water = b1*EXP( b2w*(temp-b3)/(temp-b4w) )

END FUNCTION sat_pres_water

!==============================================================================

ELEMENTAL FUNCTION sat_pres_ice(temp)
IMPLICIT NONE

REAL (KIND=wp)              :: sat_pres_ice
REAL (KIND=wp), INTENT(IN)  :: temp

sat_pres_ice = b1*EXP( b2i*(temp-b3)/(temp-b4i) )

END FUNCTION sat_pres_ice

!==============================================================================

ELEMENTAL FUNCTION spec_humi(pvap,pres)
IMPLICIT NONE

REAL (KIND=wp)              :: spec_humi
REAL (KIND=wp), INTENT(IN)  :: pres,pvap

spec_humi = rdv*pvap/( pres - o_m_rdv*pvap )

END FUNCTION spec_humi

!==============================================================================

ELEMENTAL FUNCTION rain_mue_Dm(Dm,qc)
IMPLICIT NONE

REAL (KIND=wp)              :: rain_mue_Dm
REAL (KIND=wp), INTENT(IN)  :: Dm,qc
REAL    (KIND=wp   ), PARAMETER ::  &
  zrmu0  = 2.00e+00_wp, & ! Seifert (2008) mue-Dmv-relation 
  zrmu1  = 2.00e+01_wp, & ! for the raindrop size distribution
  zrmu2  = 1.00e+03_wp, & ! 
  zrmu3  = 1.10e-03_wp, & ! 
  zrmu4  = 0.00e+00_wp

IF (Dm.LE.zrmu3) THEN    
  !      zrmue(i,j) = zrmu0*TANH((4.*zrmu2*(zDmr-zrmu3))**2) + zrmu4
  rain_mue_Dm = zrmu0*(1.-Dm/zrmu3)**2
ELSE
  !      zrmue(i,j) = zrmu1*TANH((1.*zrmu2*(zDmr-zrmu3))**2) + zrmu4
  rain_mue_Dm = zrmu0*(Dm/zrmu3-1.)**2
ENDIF
rain_mue_Dm = MAX(0.0d0,rain_mue_Dm*MIN(1.0d0,1.0d0-1.e3*qc))

END FUNCTION rain_mue_Dm

!==============================================================================

!***********************************************************************
! dimensionless growth time scale (ECHAM Code)
! B.Kaercher and S. Solomon, JGR 104(D22), 27441-27459, 1999 
!***********************************************************************

REAL(KIND=wp) ELEMENTAL FUNCTION dep_growth_timescale(b,y,x_0) 
IMPLICIT NONE
REAL(KIND=wp), INTENT(in) :: b,y,x_0

REAL(KIND=wp), PARAMETER  :: &
  SQ31 = 0.577350269_wp,       &
  SIX1 = 0.166666666_wp

REAL(KIND=wp) ::  X,F1,F10,F2,F20
    
IF (y.LE.x_0) THEN
  dep_growth_timescale = 0.
ELSE
  X     = MIN( y,DBLE(0.9999999) )
  F1    = SIX1 * LOG( (1.+X +X**2)  / (1.-X)**2 ) 
  F10   = SIX1 * LOG( (1.+x_0+x_0**2) / (1.-x_0)**2 ) 
  F2    = SQ31 * ATAN( SQ31*(1.+2.*X) )
  F20   = SQ31 * ATAN( SQ31*(1.+2.*x_0) )      
  dep_growth_timescale = (b+1.)*(F1-F10) + (b-1.)*(F2-F20)            
END IF

END FUNCTION dep_growth_timescale

!==============================================================================

!==============================================================================
!
! Functions for setting initial and boundary conditions for NCX in case
! these fields are not present in input and/ or boundary data sets. Will
! be called in subroutine organize_input() from src_input.f90:
!
!==============================================================================

! For Seifert/Beheng scheme:
! Set NCCLOUD in a way that it leads to an average mass corresponding to a diameter of 10 microns:
! Input: QC in kg/m^3; Output: NCC in 1/m^3
! Function is called from organize_input in organize_data ('init'). At this point,
! the structures cloud, rain, ice, snow, graupel and hail are not initialized yet.
! So they cannot be used here.

FUNCTION set_qnc_from_qc_sb(qcin)

  IMPLICIT NONE

  REAL(kind=wp) :: set_qnc_from_qc_sb(ie, je, ke)
  REAL(kind=wp), INTENT(in) :: qcin(ie, je, ke)

  ! Diameter of mean particle mass:
  REAL(kind=wp), PARAMETER :: Dmean = 10e-6_wp
!  REAL(kind=wp), PARAMETER :: Dmean = 30e-6_wp

  REAL(kind=wp), PARAMETER :: rhowater = 1000.0_wp
  REAL(kind=wp), PARAMETER :: pi = 3.14159265358979323844_wp

  set_qnc_from_qc_sb = &
       qcin * 6.0_wp / (pi * rhowater * Dmean**3.0_wp)

END FUNCTION set_qnc_from_qc_sb

! For Seifert/Beheng scheme:
! Set NCICE in a way that it leads to an average mass corresponding to a diameter of 100 microns:
! Input: QI in kg/m^3; Output: NCICE in 1/m^3
FUNCTION set_qni_from_qi_sb(qiin)

  IMPLICIT NONE

  REAL(kind=wp) :: set_qni_from_qi_sb(ie, je, ke)
  REAL(kind=wp), INTENT(in) :: qiin(ie, je, ke)

  ! Diameter of mean particle mass:
  REAL(kind=wp), PARAMETER :: Dmean = 100e-6_wp
  ! Mass-size-relation parameter of "iceHK"
  REAL(kind=wp), PARAMETER :: ageo  = 0.217_wp
  REAL(kind=wp), PARAMETER :: bgeo  = 0.302115_wp

  set_qni_from_qi_sb = &
       qiin / 1e-10
     !  qiin / ( ( Dmean / ageo) ** (1.0_wp / bgeo) )

END FUNCTION set_qni_from_qi_sb

! For Seifert/Beheng scheme:
! Same procedure, but for a single scalar value instead of a field:
FUNCTION set_qni_from_qi_sb_scalar(qiin)

  IMPLICIT NONE

  REAL(kind=wp) :: set_qni_from_qi_sb_scalar
  REAL(kind=wp), INTENT(in) :: qiin

  ! Diameter of mean particle mass:
  REAL(kind=wp), PARAMETER :: Dmean = 100e-6_wp
  ! Mass-size-relation parameter of "iceHK"
  REAL(kind=wp), PARAMETER :: ageo  = 0.217_wp
  REAL(kind=wp), PARAMETER :: bgeo  = 0.302115_wp

  set_qni_from_qi_sb_scalar = &
       qiin / 1e-10
     !         qiin / ( ( Dmean / ageo) ** (1.0_wp / bgeo) )

END FUNCTION set_qni_from_qi_sb_scalar

! For Seifert/Beheng scheme:
! Set NCRAIN in a way that it is consistent to the assumptions in itype_gscp == 4:
! Input: QR in kg/m^3; Output: NCRAIN in 1/m^3
FUNCTION set_qnr_from_qr_sb(qrin)

  IMPLICIT NONE

  REAL(kind=wp) :: set_qnr_from_qr_sb(ie, je, ke)
  REAL(kind=wp), INTENT(in) :: qrin(ie, je, ke)

  ! Diameter of mean particle mass:
  REAL(kind=wp), PARAMETER :: N0r = 8000.0e3_wp
  REAL(kind=wp), PARAMETER :: rhowater = 1000.0_wp
  REAL(kind=wp), PARAMETER :: pi = 3.14159265358979323844_wp

  set_qnr_from_qr_sb = &
       N0r * ( qrin * 6.0_wp / (pi * rhowater * N0r * gamma_fct(4.0_wp)))**(0.25_wp)

END FUNCTION set_qnr_from_qr_sb

! For Seifert/Beheng scheme:
! Set NCSNOW in a way that it is consistent to the assumptions in itype_gscp == 4:
! Input: QS in kg/m^3; Output: NCSNOW in 1/m^3
FUNCTION set_qns_from_qs_sb(qsin)

  IMPLICIT NONE

  REAL(kind=wp) :: set_qns_from_qs_sb(ie, je, ke)
  REAL(kind=wp), INTENT(in) :: qsin(ie, je, ke)

  ! Diameter of mean particle mass:
  REAL(kind=wp), PARAMETER :: Dmean = 100e-6_wp
  REAL(kind=wp), PARAMETER :: N0s = 800.0e3_wp
  REAL(kind=wp), PARAMETER :: ams = 0.038_wp
  REAL(kind=wp), PARAMETER :: bms = 2.0_wp

  set_qns_from_qs_sb = &
       N0s * ( qsin / ( ams * N0s * gamma_fct(bms+1.0_wp)))**( 1.0_wp/(1.0_wp+bms) )
!  set_qns_from_qs_sb = &
!       qsin / ( ( Dmean / ams) ** (1.0_wp / bms) )

END FUNCTION set_qns_from_qs_sb

! For Seifert/Beheng scheme:
! Set NCGRAUPEL in a way that it is consistent to the assumptions in itype_gscp == 4:
! Input: QG in kg/m^3; Output: NCGRAUPEL in 1/m^3
FUNCTION set_qng_from_qg_sb(qgin)

  IMPLICIT NONE

  REAL(kind=wp) :: set_qng_from_qg_sb(ie, je, ke)
  REAL(kind=wp), INTENT(in) :: qgin(ie, je, ke)

  ! Diameter of mean particle mass:
  REAL(kind=wp), PARAMETER :: Dmean = 1000e-6_wp
  REAL(kind=wp), PARAMETER :: N0g = 4000.0e3_wp
  REAL(kind=wp), PARAMETER :: amg = 169.6_wp
  REAL(kind=wp), PARAMETER :: bmg = 3.1_wp

  set_qng_from_qg_sb = &
       N0g * ( qgin / ( amg * N0g * gamma_fct(bmg+1.0_wp)))**( 1.0_wp/(1.0_wp+bmg) )

!  set_qng_from_qg_sb = &
!       qgin / ( ( Dmean / amg) ** (1.0_wp / bmg) )

END FUNCTION set_qng_from_qg_sb


!==============================================================================
!==============================================================================
!==============================================================================


SUBROUTINE seifert_pp ()


!==============================================================================
!                                              
! Two-moment mixed-phase bulk microphysics for COSMO
!       
! original version by Axel Seifert, May 2003
! with modifications by Ulrich Blahak, 2007-2012
!                                                                       
! Description:
!
! The subroutine is the interface between COSMO and the original KAMM2 modules, 
! which are provided by src_twomom_sb.f90. A major difference between COSMO and 
! KAMM2 is that KAMM2 uses mass densities instead of mixing ratios, thus
! q_cloud=rho*qc. 
! Temporary allocation of memory to the KAMM2 variables is done by the 
! subroutines ALLOC_DRIVER and ALLOC_WOLKEN.
! All microphysical source terms e.g. nucleation, condensation, coagulation,
! freezing and melting are then calculated and time integrated within the 
! subroutine CLOUDS.
! Sedimentation is done using explicit upstream advection within the four
! sedimentation subroutines, e.g. GRAUPEL_SEDIMENTATION (sedimenation of 
! cloud droplets is neglected in this version). 
!
! Various different parameterizations can be choosen using ITYPE_GSCP,
! currently the following parameters are recommended:
!
!     2    with ice phase processes turned ON, including hail 
!     4    Phillips et al. (2008) parameterization of heterogenous nucleation,
!     8    'continental' aerosol assumed for nucleation of cloud droplets
!          using the parameterization of Segal and Khain (2006)
!     3    use autoconversion/accretion scheme by Seifert&Beheng (2001) 
!          including breakup and shape effects (see Seifert 2007 for details)
!
! other possible values are
!
!  *    -> Switch for ice categories 
!  2723    with hail category
!  1723    without hail category, only ice/snow/graupel
!  0723    without ice phase processes, pure warm rain two-moment scheme
!
!   *   -> Switch for ice nucleation (IN):
!  2323    for Meyers formula
!  2423    for Phillips scheme 
!  2123    homogeneous nucleation only (do not use for real cases)
!
!    *  -> Switch for droplet nucleation (CCN):
!  2413    for 'maritime' aerosol using Eq. (17) of SB2006
!  2433    for some 'intermediate' aerosol using Eq. (17) of SB2006
!  2443    for old treatment of 'continental' CCN (Seifert 2002)
!  2463    for 'maritime' CCN using Segal and Khain parameterization
!  2473    for 'intermediate' CCN using Segal and Khain parameterization
!  2483    for 'continental' CCN using Segal and Khain parameterization
!  2493    for 'continental polluted' CCN using Segal and Khain parameterization
!
!     * -> Switch for warm rain scheme (autoconversion, accretion, breakup, ...)
!  2423    currently recommended version of the SB2001 scheme
!  0100    simple one-moment warm rain Kessler-scheme 
!
! details can be found in src_twomom_sb.F90 and in the following references
!
! Seifert und Beheng (2001): Atmos. Res. 59-60, 265-281    
! Seifert, Axel (2002): Dissertation, Uni Karlsruhe
! Seifert and Beheng (2006): Meteor. Atmos. Phys. 92, 45-66
! Seifert (2008): J. Atmos. Sci.
! Noppel et al. (2010): Atmos. Res.
!
!
! Praeprozessor-Flags: (set in Makefile or define below):
!
! MPI_CLOUD:             Define it, if you want to enable the MPI-code
!                        necessary to use load balancing.
!                        If you then really want to use load balancing,
!                        set the below parameter
!                        LOGICAL, PARAMETER :: MPI_LOAD_BALANCING = .true.
!                        Load balancing speeds up computations on
!                        scalar many-core architectures.
!
! SEDI_VECTORIZED:       If defined, vectorizable versions of the sedimenation routines
!                        are used. You can find the code in src_twomom_sb.f90. This is recommended
!                        on vector machines. On scalar machines, you can use
!                        both versions, and the original routines should be
!                        slightly faster.
!
! CGP_SEARCH_VEC:        If defined, use vectorizable version of cloudy 
!                        grid point search.
!                        Is more efficient on vector machines.
!                        On scalar machines, this switch is not necessary.
!                        
!
!==============================================================================

!!! comment out the items, which are not desired;
!!! comment all items out, if their definition is left to the Makefile!
!!! Otherwise, they take precedence over the Makefile definitions and
!!! cause compiler warnings about "redefined macros".
!#define MPI_CLOUD
!#define SEDI_VECTORIZED
!#define CGP_SEARCH_VEC

!!!!! Satad wird bei mir nicht ueberall gemacht, sondern nur an Wolkenpunkten
!!!!! Vielleicht meine Routine erst nach dem Rueckkopieren anwenden?
!!!!!
!!!!! Residuale Uebersaettigung muss nach dem Rückverteilen nur an
!!!!! Wolkenpunkten gerechnet werden, weil sie sich bis dahin anderswo nicht
!!!!! geaendert hat. Allgemein: dort Rechnen, wo Phasenuebergaenge
!!!!! stattfinden. Wenn also z.b. die satad ueberall gerechnet wird,
!!!!! dann auch die ssw, ssi ueberall nachberechnen.

  ! LM modules

  IMPLICIT NONE

  !===============================================================================
  ! .. Set configuration parameters and switches:
  !===============================================================================

  REAL(KIND=wp), PARAMETER :: qnc_const = 1e8_wp    ! Assumed constant NCCLOUD (1/m^3) in case of nuc_c_typ = 0 
  REAL(KIND=wp), PARAMETER :: lewfaktor = 1.0_wp
  LOGICAL, PARAMETER :: lgshed_ascm = .false.
  LOGICAL, PARAMETER :: lhshed_ascm = .false.
  REAL(KIND=wp), PARAMETER :: tshift_clfr_sbcm = 0.0_wp
  LOGICAL, PARAMETER :: speichere_umwandlungsraten = .false.
  LOGICAL, PARAMETER :: lwrite_gspstatistics = .false.
  LOGICAL, PARAMETER :: lwrite_dqdt_gspstatistics = .false.
  REAL(KIND=wp), PARAMETER :: z_ncn_segal = 2000.0
  REAL(KIND=wp), PARAMETER :: z1oe_ncn_segal = 2000.0
  LOGICAL, PARAMETER :: lincloud = .true.

  !===============================================================================
  !===============================================================================

  ! ... Local Variables 

  INTEGER(KIND=iintegers) :: its,ite,jts,jte,kts,kte
  INTEGER(KIND=iintegers) :: ims,ime,jms,jme,kms,kme
  INTEGER(KIND=iintegers) :: nx
  INTEGER(KIND=iintegers) :: i,j,k,ii,jj,kk,n,nt_2mom,nt_sedi,igridpoints,izstat,kkk,kkp1
  INTEGER, SAVE  :: firstcall, firstcall_init
  REAL(KIND=wp) :: dzmin, xmax(20)
  REAL(KIND=wp) :: q_vap_new,q_vap_old
  REAL(KIND=wp) :: q_liq_new,q_liq_old
  REAL(KIND=wp) :: q_ice_new,q_ice_old
  REAL(KIND=wp) :: q_v, rho_v, x_v, hlp, dt_sedi,e_v,T_a,qerr, mrho

  INTEGER, DIMENSION(0:ie*je*ke)           :: ilm,jlm,klm
  REAL(KIND=wp), DIMENSION(ie,je,ke) :: rain_r,rain_g,rain_s,rain_i,rain_h
  REAL(KIND=wp), DIMENSION(ie,je,ke) :: adz,rhocorr,ssi,ssw,zml

  REAL(KIND=wp),    PARAMETER :: eps  = 1e-15
  LOGICAL, PARAMETER :: debug = .FALSE.
  LOGICAL, PARAMETER :: CGP_SEARCH = .TRUE.
  LOGICAL, PARAMETER :: MPI_LOAD_BALANCING = .FALSE.
  LOGICAL, PARAMETER :: MPI_DEBUG = .FALSE.
  LOGICAL, PARAMETER :: debug_maxval = .true.  ! Output von Max./Min.-Werten der QX auf stdout
! UB_20090814>>
! neuer Schalter, ermoeglicht getrennte Spezifizierung, ob Output von MAX(RR) auf stdout.
! Voher war das an den Output der Max./Min.-Werte der QX gekoppelt.
  LOGICAL, PARAMETER :: debug_maxrain = .true.

! UB_20090814>>
! replaces debug_maxval at these code lines below, where all processors
! come by to do calculations. If debug_maxval = .TRUE., it is set
! to true if at least one processor has got some cloudy grid points:
  LOGICAL            :: debug_maxvalloc
  LOGICAL            :: debug_maxrainloc
! Similar thing, but it is intended for these code parts where only
! processors possessing cloudy grid points come by. If debug_maxval = .TRUE., it is set
! to true if all processors have some cloudy grid points, either because
! of load balancing or because this is just so:
  LOGICAL            :: debug_maxvalglob
! UB_20090814<<

  INTEGER (KIND=iintegers) :: izerror

  ! ... MPI stuff

  INTEGER(KIND=iintegers) :: MPI_COMM, MPI_ERR, MPI_PROCS, MPI_MYRANK, MPI_REAL, MPI_INTEGER
  INTEGER(KIND=iintegers) :: lb_total_num_cgp,lb_arrays_per_gp,lb_sendlength,&
    &                        lb_recvlength,igpmin_lbalance
  REAL(KIND=wp)       :: lb_avg_num_cgp
  INTEGER(KIND=iintegers), ALLOCATABLE   :: lb_proc_num_cgp(:),lb_proc_num_bgp(:),&
       &                                    lb_pp_num_bgp(:,:),lb_cgp_to_proc(:), &
       &                                    sendcnt(:),recvcnt(:),&
       &                                    sdispl(:),rdispl(:)
  REAL(KIND=wp), ALLOCATABLE         :: sendbuf(:),recvbuf(:)
  INTEGER(KIND=iintegers), ALLOCATABLE   :: isendbuf(:),irecvbuf(:)
  LOGICAL :: lb_imbalance

! field for inst. microphys. conv. rates     (kg/kg*s)
  REAL(KIND=wp), ALLOCATABLE         :: &
       dqdt_cl(:,:,:), heizrate(:,:,:), qalt(:,:,:), qnalt(:,:,:)      
  INTEGER(KIND=iintegers), ALLOCATABLE         :: &
       dqdt_cl_flag(:,:,:,:)  ! flag-field for inst. microphys. conv. rates     ( - )
  INTEGER(KIND=iintegers) :: iiii, unitnr
  REAL(KIND=wp) :: dum1, dum2, dum3, dum4, dum5, dum6, dum7

  REAL(KIND=wp), PARAMETER :: wcb_min=0.1, scb_min=0.0e0, scb_max=-0.0e0


! UB_20090901>> Working arrays for vectorized grid point search:
  INTEGER(kind=iintegers), DIMENSION(ie,je*ke) :: jo_cgp, ko_cgp
  INTEGER(kind=iintegers), DIMENSION(ie)       :: njko_cgp
! UB_20090901<<

  ! timestep used for micrphysics integration:
  REAL(KIND=wp) :: zdt

  DOUBLE PRECISION :: dummy3(3)

  ! local constants for t and pp equation:
  REAL    (KIND=wp   ) ::  &
    zlh_s, zlh_f, &  ! Temperature dependent lh_s and lh_f
    zttfac, zpptfac, zlhTfac, &
    cv_d, cv_v, cv_l, cv_i

  INTEGER (KIND=iintegers) ::  &
    satad_errstat,  &  ! error status flag for satad_v_3d
    satad_count        ! number of actually performed Newton iterations during satad_v_3d

  ! Fields for satad from meteo_utilities:
  REAL    (KIND=wp   ) ::  &
       zpres(ie,je), zdummy(ie,je,9)

  CHARACTER (LEN=80)          :: yzerrmsg
  CHARACTER (LEN=80)          :: yzroutine

! Tracer pointers
!----------------
  REAL (KIND=wp), POINTER :: &
    qv  (:,:,:) => NULL(),         & ! QV at nx
    qc  (:,:,:) => NULL(),         & ! QC at nx
    qr  (:,:,:) => NULL(),         & ! QR at nx
    qi  (:,:,:) => NULL(),         & ! QI at nx
    qs  (:,:,:) => NULL(),         & ! QS at nx
    qg  (:,:,:) => NULL(),         & ! QG at nx
    qh  (:,:,:) => NULL(),         & ! QH at nx
    qnc (:,:,:) => NULL(),         & ! NCCLOUD at nx
    qnr (:,:,:) => NULL(),         & ! NCRAIN at nx
    qni (:,:,:) => NULL(),         & ! NCICE at nx
    qns (:,:,:) => NULL(),         & ! NCSNOW at nx
    qng (:,:,:) => NULL(),         & ! NCGRAUPEL at nx
    qnh (:,:,:) => NULL()            ! NCHAIL at nx

!------------ End of header ---------------------------------------------------

!------------------------------------------------------------------------------
! Begin Subroutine seifert_pp
!------------------------------------------------------------------------------

  yzroutine(:) = ' '
  yzroutine = 'seifert_pp'

  IF (isIO()) WRITE (*,*) TRIM(yzroutine)//": start"

  ims = 1 
  ime = ie
  jms = 1 
  jme = je
  kms = 1
  kme = ke
  its = istartpar  ! 1
  ite = iendpar    ! ie
  jts = jstartpar  ! 1
  jte = jendpar    ! je
  kts = 1
  kte = ke

  ! Dimensions for field allocations:
  dim_ix = (jte-jts+1)*(kte-kts+1)*(ite-its+1)
  dim_iy = 1
  dim_iz = 1

  ! Actual computing dimensions (will be set below, here just initialization):
  loc_ix = -999
  loc_iy = -999
  loc_iz = -999

  IF (ltime) CALL get_timings (i_2mom_start, ntstep, dt, izerror)

  graupel_shedding = lgshed_ascm
  hail_shedding    = lhshed_ascm
  tshift_cloudfreeze = tshift_clfr_sbcm
  qnc_const_seifert  = qnc_const
  z_ncn_segal_seifert = z_ncn_segal
  z1oe_ncn_segal_seifert = z1oe_ncn_segal

  ! select timelevel and timestep for calculations
  nx    = nnew
  IF ( l2tls ) THEN
    zdt   = dt
  ELSE
    zdt   = dt2
  ENDIF

  ! retrieve the required microphysics tracers (at corresponding timelevel nnew)
  yzerrmsg(:) = ' '
  CALL trcr_get(izerror, 'QV', ptr_tlev = nx, ptr = qv)
  IF (izerror /= 0) THEN
    yzerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
  ENDIF
  CALL trcr_get(izerror, 'QC', ptr_tlev = nx, ptr = qc)
  IF (izerror /= 0) THEN
    yzerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
  ENDIF
  CALL trcr_get(izerror, 'QR', ptr_tlev = nx, ptr = qr)
  IF (izerror /= 0) THEN
    yzerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
  ENDIF
  CALL trcr_get(izerror, 'QI', ptr_tlev = nx, ptr = qi)
  IF (izerror /= 0) THEN
    yzerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
  ENDIF
  CALL trcr_get(izerror, 'QS', ptr_tlev = nx, ptr = qs)
  IF (izerror /= 0) THEN
    yzerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
  ENDIF
  CALL trcr_get(izerror, 'QG', ptr_tlev = nx, ptr = qg)
  IF (izerror /= 0) THEN
    yzerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
  ENDIF
  CALL trcr_get(izerror, 'QH', ptr_tlev = nx, ptr = qh)
  IF (izerror /= 0) THEN
    yzerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
  ENDIF
  CALL trcr_get(izerror, 'NCCLOUD', ptr_tlev = nx, ptr = qnc)
  IF (izerror /= 0) THEN
    yzerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
  ENDIF
  CALL trcr_get(izerror, 'NCRAIN', ptr_tlev = nx, ptr = qnr)
  IF (izerror /= 0) THEN
    yzerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
  ENDIF
  CALL trcr_get(izerror, 'NCICE', ptr_tlev = nx, ptr = qni)
  IF (izerror /= 0) THEN
    yzerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
  ENDIF
  CALL trcr_get(izerror, 'NCSNOW', ptr_tlev = nx, ptr = qns)
  IF (izerror /= 0) THEN
    yzerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
  ENDIF
  CALL trcr_get(izerror, 'NCGRAUPEL', ptr_tlev = nx, ptr = qng)
  IF (izerror /= 0) THEN
    yzerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
  ENDIF
  CALL trcr_get(izerror, 'NCHAIL', ptr_tlev = nx, ptr = qnh)
  IF (izerror /= 0) THEN
    yzerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
  ENDIF


  ! Delete values of the precipitation rate from last timestep:
  prr_gsp = 0.0_wp
  prs_gsp = 0.0_wp
  prg_gsp = 0.0_wp
  prh_gsp = 0.0_wp

  speichere_precipstat = .FALSE.

  ! Initialisierung der Aufzeichnung der gesamten Kondensationsmenge im Zeitschritt:
  IF (lwrite_dqdt_gspstatistics  .AND. lwrite_gspstatistics) THEN
    IF (my_cart_id == 0) THEN
      WRITE (*,*) TRIM(yzroutine)//' *** lwrite_gspstatistics = .true. and lwrite_dqdt_gspstatistics = .true., ',&
           'but preprocessor flag SAVE_CONVERSIONRATES has not been defined!'
    END IF
    CALL model_abort (my_world_id, 11124, TRIM(yzroutine)//': Preprozessor flag error, SAVE_CONVERSIONRATES', &
         TRIM(yzroutine)//': Preprozessor flag error')      
  END IF


  MPI_COMM    = icomm_cart
  MPI_REAL    = imp_reals
  MPI_INTEGER = imp_integers
  MPI_PROCS   = num_compute      !.. CALL mpi_comm_size(MPI_COMM,MPI_PROCS,MPI_ERR)
  MPI_MYRANK  = my_cart_id       !.. CALL mpi_comm_rank(MPI_COMM,MPI_MYRANK,MPI_ERR)
  IF (MPI_DEBUG) THEN
    WRITE (*,'(2(A,I5))') TRIM(yzroutine)//": MPI_PROCS   = ",MPI_PROCS,  "    MPI_MYRANK = ",MPI_MYRANK
  ENDIF

  IF (isIO()) WRITE (*,*) TRIM(yzroutine)//": itype_gscp = ",itype_gscp

  cloud_type = itype_gscp

!!! Wird in organize_physics.f90 gemacht:
!!! CALL init_seifert( cloud_type )
#ifdef NUDGING
  ! add part of latent heating calculated in subroutine seifert_pp to model latent
  ! heating field: subtract temperature from model latent heating field
  IF ((llhn .OR. llhnverif)) THEN
    IF (lhn_qrs) THEN
      qrsflux(:,:,:) = 0.0_wp
    ENDIF
    CALL get_gs_lheating ('add',1,ke)
  ENDIF
#endif
  IF ( ldiabf_lh ) THEN
    ! initialize temperature increment due to latent heat (for Jochen)
    tinc_lh(:,:,:) = tinc_lh(:,:,:) - t(:,:,:,nx)
  END IF

  ! Other initializations:
  IF (firstcall_init /= 1) THEN

    ! Initialize the lookup table for wet growth of graupel -> hail:
    ! This is done only on one processor, because the lookup table file
    ! is read from hard disk:
    IF (isIO()) WRITE (*,*) TRIM(yzroutine)//": init_dmin_wetgrowth"
    IF (MPI_MYRANK == 0) THEN
      CALL get_free_unit(unitnr)
      CALL init_dmin_wetgrowth('dmin_'//TRIM(ADJUSTL(graupel%name))//'_wetgrowth_lookup.dat', unitnr)
      CALL release_unit(unitnr)
    END IF

    !.. In MPI-run, broadcast the lookup tables for wet growth to all other nodes:
    IF (MPI_PROCS > 1) THEN
      CALL MPI_BCAST(anzp_wg, 1, mpi_integer, 0, MPI_COMM, MPI_err)
      CALL MPI_BCAST(anzT_wg, 1, mpi_integer, 0, MPI_COMM, MPI_err)
      CALL MPI_BCAST(anzw_wg, 1, mpi_integer, 0, MPI_COMM, MPI_err)
      CALL MPI_BCAST(anzi_wg, 1, mpi_integer, 0, MPI_COMM, MPI_err)
      
      IF (MPI_MYRANK > 0) THEN
        ALLOCATE(pvec_wg_g(anzp_wg))
        ALLOCATE(Tvec_wg_g(anzT_wg))
        ALLOCATE(qwvec_wg_g(anzw_wg))
        ALLOCATE(qivec_wg_g(anzi_wg))
        ALLOCATE(dmin_wg_g(anzp_wg,anzT_wg,anzw_wg,anzi_wg))
      END IF
      ALLOCATE(sendbuf(anzp_wg*anzT_wg*anzw_wg*anzi_wg))
      
      IF (MPI_MYRANK == 0) sendbuf(1:anzp_wg) = REAL(pvec_wg_g,kind=wp)
      CALL MPI_BCAST(sendbuf(1:anzp_wg), anzp_wg, mpi_real, 0, MPI_COMM, MPI_err)
      IF (MPI_MYRANK > 0) pvec_wg_g = sendbuf(1:anzp_wg)
      IF (MPI_MYRANK == 0) sendbuf(1:anzT_wg) = REAL(Tvec_wg_g,kind=wp)
      CALL MPI_BCAST(sendbuf(1:anzT_wg), anzT_wg, mpi_real, 0, MPI_COMM, MPI_err)
      IF (MPI_MYRANK > 0) Tvec_wg_g = sendbuf(1:anzT_wg)
      IF (MPI_MYRANK == 0) sendbuf(1:anzw_wg) = REAL(qwvec_wg_g,kind=wp)
      CALL MPI_BCAST(sendbuf(1:anzw_wg), anzw_wg, mpi_real, 0, MPI_COMM, MPI_err)
      IF (MPI_MYRANK > 0) qwvec_wg_g = sendbuf(1:anzw_wg)
      IF (MPI_MYRANK == 0) sendbuf(1:anzi_wg) = REAL(qivec_wg_g,kind=wp)
      CALL MPI_BCAST(sendbuf(1:anzi_wg), anzi_wg, mpi_real, 0, MPI_COMM, MPI_err)
      IF (MPI_MYRANK > 0) qivec_wg_g = sendbuf(1:anzi_wg)
      
      IF (MPI_MYRANK == 0) THEN
        jj = 0
        DO ii=1, anzi_wg
          DO k=1, anzw_wg
            DO j=1, anzT_wg
              DO i=1, anzp_wg
                jj = jj + 1
                sendbuf(jj) = dmin_wg_g(i,j,k,ii)
              END DO
            END DO
          END DO
        END DO
      END IF
      
      CALL MPI_BCAST(jj, 1, mpi_integer, 0, MPI_COMM, MPI_err)
      CALL MPI_BCAST(sendbuf, jj, mpi_real, 0, MPI_COMM, MPI_err)
      
      IF (MPI_MYRANK > 0) THEN
        jj = 0
        DO ii=1, anzi_wg
          DO k=1, anzw_wg
            DO j=1, anzT_wg
              DO i=1, anzp_wg
                jj = jj + 1
                dmin_wg_g(i,j,k,ii) = sendbuf(jj)
              END DO
            END DO
          END DO
        END DO
      END IF
      
      DEALLOCATE(sendbuf)
    
    END IF

    IF (MPI_MYRANK == 0) WRITE (*,*) TRIM(yzroutine)//': init_dmin_wetgrowth done!'

! UB_20090227>>
!     Initialize an equidistant lookup table for better vectorization of 
!     the table lookup itself:
    IF (isIO()) WRITE (*,*) TRIM(yzroutine)//": init_dmin_wg_gr_ltab_equi"
    IF (MPI_MYRANK == 0) THEN
      CALL get_free_unit(unitnr)
      CALL init_dmin_wg_gr_ltab_equi('dmin_'//TRIM(ADJUSTL(graupel%name))//'_wetgrowth_lookup.dat', &
           unitnr, 61, ltabdminwgg)
      CALL release_unit(unitnr)
    END IF
    
    !.. In MPI-run, broadcast the lookup table for wet growth to all other nodes:
    IF (MPI_PROCS > 1) THEN
      
      !.. 1) Broadcast scalars of the struct ltabdminwgg:
      
      ALLOCATE(isendbuf(4))
      IF (MPI_MYRANK == 0) THEN
        isendbuf(1) = ltabdminwgg%n1
        isendbuf(2) = ltabdminwgg%n2
        isendbuf(3) = ltabdminwgg%n3
        isendbuf(4) = ltabdminwgg%n4
      END IF
      CALL MPI_BCAST(isendbuf, 4, mpi_integer, 0, MPI_COMM, MPI_err)
      IF (MPI_MYRANK > 0) THEN
        ltabdminwgg%n1 = isendbuf(1)
        ltabdminwgg%n2 = isendbuf(2)
        ltabdminwgg%n3 = isendbuf(3)
        ltabdminwgg%n4 = isendbuf(4)
      END IF
      DEALLOCATE(isendbuf)

      ALLOCATE(sendbuf(8))
      IF (MPI_MYRANK == 0) THEN
        sendbuf(1) = ltabdminwgg%dx1
        sendbuf(2) = ltabdminwgg%odx1
        sendbuf(3) = ltabdminwgg%dx2
        sendbuf(4) = ltabdminwgg%odx2
        sendbuf(5) = ltabdminwgg%dx3
        sendbuf(6) = ltabdminwgg%odx3
        sendbuf(7) = ltabdminwgg%dx4
        sendbuf(8) = ltabdminwgg%odx4
      END IF
      CALL MPI_BCAST(sendbuf,  8, mpi_real, 0, MPI_COMM, MPI_err)
      IF (MPI_MYRANK > 0) THEN
        ltabdminwgg%dx1 = sendbuf(1)
        ltabdminwgg%odx1 = sendbuf(2)
        ltabdminwgg%dx2 = sendbuf(3)
        ltabdminwgg%odx2 = sendbuf(4)
        ltabdminwgg%dx3 = sendbuf(5)
        ltabdminwgg%odx3 = sendbuf(6)
        ltabdminwgg%dx4 = sendbuf(7)
        ltabdminwgg%odx4 = sendbuf(8)
      END IF
      DEALLOCATE(sendbuf)

      !.. 2) Allocate/Broadcast table vectors:

      IF (MPI_MYRANK > 0) THEN
        ALLOCATE( ltabdminwgg%x1(ltabdminwgg%n1) )
        ALLOCATE( ltabdminwgg%x2(ltabdminwgg%n2) )
        ALLOCATE( ltabdminwgg%x3(ltabdminwgg%n3) )
        ALLOCATE( ltabdminwgg%x4(ltabdminwgg%n4) )
      END IF

      ALLOCATE(sendbuf(ltabdminwgg%n1*ltabdminwgg%n2*ltabdminwgg%n3*ltabdminwgg%n4))

      IF (MPI_MYRANK == 0) sendbuf(1:ltabdminwgg%n1) = REAL(ltabdminwgg%x1)
      CALL MPI_BCAST(sendbuf(1:ltabdminwgg%n1), ltabdminwgg%n1, mpi_real, 0, MPI_COMM, MPI_err)
      IF (MPI_MYRANK > 0) ltabdminwgg%x1 = sendbuf(1:ltabdminwgg%n1)
      IF (MPI_MYRANK == 0) sendbuf(1:ltabdminwgg%n2) = REAL(ltabdminwgg%x2)
      CALL MPI_BCAST(sendbuf(1:ltabdminwgg%n2), ltabdminwgg%n2, mpi_real, 0, MPI_COMM, MPI_err)
      IF (MPI_MYRANK > 0) ltabdminwgg%x2 = sendbuf(1:ltabdminwgg%n2)
      IF (MPI_MYRANK == 0) sendbuf(1:ltabdminwgg%n3) = REAL(ltabdminwgg%x3)
      CALL MPI_BCAST(sendbuf(1:ltabdminwgg%n3), ltabdminwgg%n3, mpi_real, 0, MPI_COMM, MPI_err)
      IF (MPI_MYRANK > 0) ltabdminwgg%x3 = sendbuf(1:ltabdminwgg%n3)
      IF (MPI_MYRANK == 0) sendbuf(1:ltabdminwgg%n4) = REAL(ltabdminwgg%x4)
      CALL MPI_BCAST(sendbuf(1:ltabdminwgg%n4), ltabdminwgg%n4, mpi_real, 0, MPI_COMM, MPI_err)
      IF (MPI_MYRANK > 0) ltabdminwgg%x4 = sendbuf(1:ltabdminwgg%n4)


      !.. 3) Allocate/Broadcast the table:

      IF (MPI_MYRANK > 0) THEN
        ALLOCATE( ltabdminwgg%ltable(ltabdminwgg%n1,ltabdminwgg%n2,ltabdminwgg%n3,ltabdminwgg%n4) )
      END IF
      
      IF (MPI_MYRANK == 0) THEN
        jj = 0
        DO ii=1, ltabdminwgg%n4
          DO k=1, ltabdminwgg%n3
            DO j=1, ltabdminwgg%n2
              DO i=1, ltabdminwgg%n1
                jj = jj + 1
                sendbuf(jj) = ltabdminwgg%ltable(i,j,k,ii)
              END DO
            END DO
          END DO
        END DO
      END IF
      
      CALL MPI_BCAST(jj, 1, mpi_integer, 0, MPI_COMM, MPI_err)
      CALL MPI_BCAST(sendbuf, jj, mpi_real, 0, MPI_COMM, MPI_err)
      
      IF (MPI_MYRANK > 0) THEN
        jj = 0
        DO ii=1, ltabdminwgg%n4
          DO k=1, ltabdminwgg%n3
            DO j=1, ltabdminwgg%n2
              DO i=1, ltabdminwgg%n1
                 jj = jj + 1
                ltabdminwgg%ltable(i,j,k,ii) = sendbuf(jj)
              END DO
            END DO
          END DO
        END DO
      END IF
     
      DEALLOCATE(sendbuf)

    END IF
! UB_20090227<<

    IF (MPI_MYRANK == 0) WRITE (*,*) TRIM(yzroutine)//': init_dmin_wg_gr_ltab_equi done!'

    firstcall_init = 1

    IF (ltime) CALL get_timings (i_2mom_init_0, ntstep, dt, izerror)

  END IF

  ! Limit QX and NCX after the advection routines at all gridpoints:
  ! This is to prevent artifacts and model crashes when there is qnx but no qx, so that
  ! at such gridpoints, no cloud microphysics is computed but the stuff gets advected.
  ! For example, in situations with high NCCLOUD at gridpoints with very low QC, the advection routine
  ! led to ever increasing values of NCCLOUD above 10^34 with nothing to stop it, which is the maximum value which fits into
  ! the grib data routines of DWD. This inevitably leads to a model crash for no apparent physical reason.
  !
  ! Clip values of QX < 0 and then clip NCX values to their allowed range:

  qv(:,:,:)  = MAX(qv(:,:,:),  0.0_wp)
  qc(:,:,:)  = MAX(qc(:,:,:),  0.0_wp)
  qnc(:,:,:) = MAX(qnc(:,:,:), 0.0_wp)
  qr(:,:,:)  = MAX(qr(:,:,:),  0.0_wp)
  qnr(:,:,:) = MAX(qnr(:,:,:), 0.0_wp)
  IF (itype_gscp >= 1000) THEN
    qi(:,:,:)  = MAX(qi(:,:,:),  0.0_wp)
    qni(:,:,:) = MAX(qni(:,:,:), 0.0_wp)
    qs(:,:,:)  = MAX(qs(:,:,:),  0.0_wp)
    qns(:,:,:) = MAX(qns(:,:,:), 0.0_wp)
    qg(:,:,:)  = MAX(qg(:,:,:),  0.0_wp)
    qng(:,:,:) = MAX(qng(:,:,:), 0.0_wp)
  END IF
  IF (itype_gscp >= 2000) THEN
    qh(:,:,:)  = MAX(qh(:,:,:),  0.0_wp)
    qnh(:,:,:) = MAX(qnh(:,:,:), 0.0_wp)
  END IF

  qnc(:,:,:) = MAX(qnc(:,:,:), qc(:,:,:) / cloud%x_max)
  qnc(:,:,:) = MIN(qnc(:,:,:), qc(:,:,:) / cloud%x_min)
  qnr(:,:,:) = MAX(qnr(:,:,:), qr(:,:,:) / rain%x_max)
  qnr(:,:,:) = MIN(qnr(:,:,:), qr(:,:,:) / rain%x_min)
  IF (itype_gscp >= 1000) THEN
    qni(:,:,:) = MAX(qni(:,:,:), qi(:,:,:) / ice%x_max)
    qni(:,:,:) = MIN(qni(:,:,:), qi(:,:,:) / ice%x_min)
    qns(:,:,:) = MAX(qns(:,:,:), qs(:,:,:) / snow%x_max)
    qns(:,:,:) = MIN(qns(:,:,:), qs(:,:,:) / snow%x_min)
    qng(:,:,:) = MAX(qng(:,:,:), qg(:,:,:) / graupel%x_max)
    qng(:,:,:) = MIN(qng(:,:,:), qg(:,:,:) / graupel%x_min)
  END IF
  IF (itype_gscp >= 2000) THEN
    qnh(:,:,:) = MAX(qnh(:,:,:), qh(:,:,:) / hail%x_max)
    qnh(:,:,:) = MIN(qnh(:,:,:), qh(:,:,:) / hail%x_min)
  END IF



  ! calculate total density rho(:,:,:):
  DO kk = kts, kte
!CDIR COLLAPSE
    qrs(:,:,kk) = qr(:,:,kk) + &
         qi(:,:,kk) + qs(:,:,kk) + &
         qg(:,:,kk) + qh(:,:,kk)
  END DO
  CALL calrho ( t(:,:,:,nx), pp(:,:,:,nx), &
       qv(:,:,:), qc(:,:,:),&
       qrs(:,:,:), p0(:,:,:), &
       rho(:,:,:), ie, je, ke, r_d, r_v/r_d-1.0_wp)

#ifdef COSMOART
  ALLOCATE( calc_i(0:dim_ix), calc_j(0:dim_ix), calc_k(0:dim_ix))  
#endif

  IF (debug.AND.isIO()) WRITE (*,*) TRIM(yzroutine)//": allocation of hydrometeor variables, dim_ix = ", dim_ix
  ! ... Allocate memory to 2MOM variables  
  CALL alloc_driver()
  CALL alloc_wolken()

  IF (debug.AND.isIO()) WRITE (*,*) TRIM(yzroutine)//": calculating supersaturations after dynamics/advection"
  !..calculate supersaturations after dynamics/advection
  DO kk = kts, kte

     ! supersaturation w.r.t. ice
    ssi(:,:,kk) = e_es_vec(DBLE(T(:,:,kk,nx)), ie, je)
    ssi(:,:,kk) = r_v & 
          * rho(:,:,kk) * qv(:,:,kk) &
          * T(:,:,kk,nx) / ssi(:,:,kk) - 1.0_wp

     ! supersaturation w.r.t. liquid water
    ssw(:,:,kk) = e_ws_vec(DBLE(T(:,:,kk,nx)), ie, je)
    ssw(:,:,kk) = r_v & 
          * rho(:,:,kk) * qv(:,:,kk) &
          * T(:,:,kk,nx) / ssw(:,:,kk) - 1.0_wp

  ENDDO



  IF (ltime) CALL get_timings (i_2mom_init, ntstep, dt, izerror)

  !..search for cloudy grid points and store locations
#ifdef CGP_SEARCH_VEC

! UB_20090901>> vectorized version:  
  IF (CGP_SEARCH .OR. MPI_LOAD_BALANCING) THEN

    ! working arrays:
    jo_cgp  = -1
    ko_cgp  = -1
    njko_cgp =  0
    
    ! grid point search part with "intelligent" storage 
    ! of locations in the working arrays:
    DO kk = kts, kte
      DO jj = jts, jte
!CDIR NODEP
        DO ii = its, ite           
          IF ( ssi(ii,jj,kk) >= eps .or. &
               & ssw(ii,jj,kk) >= eps .or. &
               & qc(ii,jj,kk)  >= eps .or. &
               & qrs(ii,jj,kk)  >= eps ) THEN
            njko_cgp(ii) = njko_cgp(ii) + 1
            jo_cgp (ii,njko_cgp(ii)) = jj
            ko_cgp (ii,njko_cgp(ii)) = kk
          ENDIF
        ENDDO
      ENDDO
    ENDDO

    ! non-vectorized loop; the inner field operations are vectorized, however:
    i = 0
!CDIR NOVECTOR
    DO ii=its, ite
      IF (njko_cgp(ii) > 0) THEN
        ilm(i:i+njko_cgp(ii)-1) = ii
        jlm(i:i+njko_cgp(ii)-1) = jo_cgp(ii,1:njko_cgp(ii))
        klm(i:i+njko_cgp(ii)-1) = ko_cgp(ii,1:njko_cgp(ii))
        i = i + njko_cgp(ii)
      END IF
    END DO
    igridpoints = i - 1
  ELSE
    DO kk = kts, kte
      DO jj = jts, jte
        i = (kk-kts) * (jte-jts+1) * (ite-its+1) + (jj-jts) * (ite-its+1)
!CDIR NODEP
        DO ii = its, ite
          ilm(i+ii-its) = ii
        END DO
        jlm(i:i+ite-its) = jj
        klm(i:i+ite-its) = kk
      ENDDO
    ENDDO
    igridpoints = (kte-kts+1) * (jte-jts+1) * (ite-its+1) - 1
    ! igridpoints is the end index of the gridpoint loop, which below starts with 0, not 1 !
  ENDIF

#else

! UB_20090901>> old version, not safely vectorizable because of dependency on i:
  i = -1
  IF (CGP_SEARCH .OR. MPI_LOAD_BALANCING) THEN
    DO kk = kts, kte
      DO jj = jts, jte
        DO ii = its, ite           
          IF ( ssi(ii,jj,kk) >= eps .or. &
               & ssw(ii,jj,kk) >= eps .or. &
               & qc(ii,jj,kk)  >= eps .or. &
               & qrs(ii,jj,kk)  >= eps ) THEN
            i = i+1
            ilm(i) = ii     ! they start with i=0
            jlm(i) = jj
            klm(i) = kk
          ENDIF
        ENDDO
      ENDDO
    ENDDO
  ELSE
    DO kk = kts, kte
      DO jj = jts, jte
        DO ii = its, ite           
          i = i+1
          ilm(i) = ii
          jlm(i) = jj
          klm(i) = kk
        ENDDO
      ENDDO
    ENDDO
  ENDIF
  igridpoints = i

#endif

  IF (ltime) CALL get_timings (i_2mom_cgps, ntstep, dt, izerror)

  lb_imbalance = .FALSE.
  ALLOCATE(lb_proc_num_cgp(MPI_PROCS))
  IF (MPI_DEBUG) WRITE (*,*) "calling mpi_allgather, MPI_PROCS = ", MPI_PROCS
  IF (MPI_PROCS > 1) THEN
    CALL mpi_allgather(igridpoints+1,1,mpi_integer,lb_proc_num_cgp,1,mpi_integer,mpi_comm,MPI_err)
  ELSE
    lb_proc_num_cgp(1) = igridpoints+1
  ENDIF
#ifdef MPI_CLOUD
  IF (MPI_LOAD_BALANCING) THEN
    lb_total_num_cgp = sum(lb_proc_num_cgp)
    lb_avg_num_cgp = lb_total_num_cgp/MPI_PROCS
    IF (MPI_DEBUG) print *,'lb_my_num_cgp =',igridpoints+1
    IF (MPI_DEBUG) print *,'lb_total_num_cgp =',lb_total_num_cgp
    IF (MPI_DEBUG) print *,'lb_avg_num_cgp =',lb_avg_num_cgp
    ! Mindestens ein Prozessor muss mehr als MPI_PROCS Wolkenpunkte aufweisen, sonst
    ! passiert im weiteren Verlauf ein Segmentation fault! Des weiteren sollte die Wolkenpunktzahl
    ! dieses Prozessors deutlich groesser sein als die Prozessorzahl MPI_PROCS, sonst
    ! lohnt sich das load balancing nicht.
    igpmin_lbalance = 10*MPI_PROCS
    IF (MPI_PROCS > 1 .AND. MAXVAL(lb_proc_num_cgp) > igpmin_lbalance &
         & .AND. REAL(MAXVAL(lb_proc_num_cgp)) > 2.00 * REAL(MINVAL(lb_proc_num_cgp))) &
         &  lb_imbalance = .TRUE.
    IF (debug.and.isIO()) print *,'lb_proc_num_cgp =',lb_proc_num_cgp,lb_imbalance
  END IF
  IF (isIO()) WRITE (*,*) TRIM(yzroutine)//": lb_imbalance = ",lb_imbalance
#endif
  IF (isIO() .AND. (CGP_SEARCH .OR. MPI_LOAD_BALANCING)) &
       WRITE (*,'(a,i10,a,f8.4,a)') TRIM(yzroutine)//": cloudy grid points = ", SUM(lb_proc_num_cgp), &
!!$       ' (', SUM(lb_proc_num_cgp) * 100.0 / ((ie_tot-nboundlines)*(je_tot-nboundlines)*ke),' % of total Domain)'
       ' (', SUM(lb_proc_num_cgp) * 100.0 / (ie_tot*je_tot*ke),' % of total Domain)'

! UB_20090814>>
  IF (debug_maxval .and. MAXVAL(lb_proc_num_cgp) > 0 ) THEN
    ! only if debug_maxval=.true. and at least one processor has some cloudy grid points
    debug_maxvalloc = .TRUE.
  ELSE
    debug_maxvalloc = .FALSE.
  END IF
  IF (debug_maxrain .and. MAXVAL(lb_proc_num_cgp) > 0 ) THEN
    ! only if debug_maxval=.true. and at least one processor has some cloudy grid points
    debug_maxrainloc = .TRUE.
  ELSE
    debug_maxrainloc = .FALSE.
  END IF

  IF (debug_maxval .and. (MINVAL(lb_proc_num_cgp).gt.0 .or. lb_imbalance)) THEN
    ! only if debug_maxval=.true. and all processors have some cloudy grid points,
    ! either because of load balancing or just because it is so.
    debug_maxvalglob = .TRUE.
  ELSE
    debug_maxvalglob = .FALSE.
  END IF
! UB_20090814<<

  ! Regenratenberechnung und -ausgabe vorbereiten:
  rain_r = 0.0
  rain_g = 0.0
  rain_h = 0.0
  rain_s = 0.0
  rain_i = 0.0

!.. temporarily here, but later in src_setup.f90:
cv_d     =  cp_d - r_d
cv_v     =  (rcpv + 1.0_wp) * cp_d - r_v
cv_l     =  (rcpl + 1.0_wp) * cp_d
cv_i     =  2060.0_wp

!!$ UB>>
! Precalculated coefficients of t and pp tendencies (DO NOT TOUCH!):
!!$  IF (ldiab_isochoric) THEN
!!$    ! .. 1. possibility: take into account the phase change terms Qh and Qm in the pressure equation,
!!$    !       which in turn give an additional contribution in the T-equation due to the pressure term
!!$    !       (this effectively leads to changing cp_d to cv_d and is the "isochoric" formulation):
!!$    zttfac  = 1.0_wp / cv_d  ! isochoric
!!$    zpptfac = 1.0_wp / cv_d  ! isochoric
!!$    zlhTfac = 1.0_wp         ! T-dependent lh_v, lh_s and lh_f and the Qm-term
!!$  ELSE
!!$    ! .. 2. possibility: disregard the phase change terms Qh and Qm in the pressure equation
!!$    !       and their additional contribution (through the p-term) to the T-equation.
!!$    !       This is the original formulation of the COSMO equation system:
    zttfac  = cpdr
    zpptfac = 0.0_wp  ! zero out the pp-tendency
    zlhTfac = 0.0_wp  ! zero out the T-dependence of lh_v, lh_s and lh_f and the Qm-term


!!$  END IF

  IF (debug_maxvalloc) THEN
    xmax(1)  = MAXVAL(w(its:ite,jts:jte,kts:kte,nx))
    xmax(2)  = MAXVAL(qv(its:ite,jts:jte,kts:kte))
    xmax(3)  = MAXVAL(qc(its:ite,jts:jte,kts:kte))
    xmax(4)  = MAXVAL(qr(its:ite,jts:jte,kts:kte))
    xmax(5)  = MAXVAL(qi(its:ite,jts:jte,kts:kte))
    xmax(6)  = MAXVAL(qs(its:ite,jts:jte,kts:kte))
    xmax(7)  = MAXVAL(qg(its:ite,jts:jte,kts:kte))  
    xmax(8)  = MAXVAL(qh(its:ite,jts:jte,kts:kte))  
    xmax(9)  = MAXVAL(qnc(its:ite,jts:jte,kts:kte))  
    xmax(10) = MAXVAL(qni(its:ite,jts:jte,kts:kte))  
    xmax(11) = MAXVAL(t(its:ite,jts:jte,kts:kte,nx))  
    IF (num_compute > 1) THEN
      CALL global_values (xmax,11,'MAX',imp_reals,icomm_cart,-1,yzerrmsg,izerror)
    END IF
    IF (isIO()) THEN
      WRITE (*,*) TRIM(yzroutine)//": output 1"
      WRITE(*,'(10x,a)') 'Maximum Values:'
      WRITE(*,'(A10,11A11)')   '   ', 'w','qv','qc','qr','qi','qs','qg','qh','nc','ni','tmax'
      WRITE(*,'(A10,11ES11.3)') '   ', (xmax(i), i=1,11)
    ENDIF
    xmax(1)  = MINVAL(w(its:ite,jts:jte,kts:kte,nx))
    xmax(2)  = MINVAL(qv(its:ite,jts:jte,kts:kte))
    xmax(3)  = MINVAL(qc(its:ite,jts:jte,kts:kte))
    xmax(4)  = MINVAL(qr(its:ite,jts:jte,kts:kte))
    xmax(5)  = MINVAL(qi(its:ite,jts:jte,kts:kte))
    xmax(6)  = MINVAL(qs(its:ite,jts:jte,kts:kte))
    xmax(7)  = MINVAL(qg(its:ite,jts:jte,kts:kte))  
    xmax(8)  = MINVAL(qh(its:ite,jts:jte,kts:kte))  
    xmax(9)  = MINVAL(qnc(its:ite,jts:jte,kts:kte))  
    xmax(10) = MINVAL(qni(its:ite,jts:jte,kts:kte))  
    xmax(11) = MINVAL(t(its:ite,jts:jte,kts:kte,nx))  
    IF (num_compute > 1) THEN
      CALL global_values (xmax,11,'MIN',imp_reals,icomm_cart,-1,yzerrmsg,izerror)
    END IF
    IF (isIO()) THEN
      WRITE(*,'(10x,a)') 'Minimum Values:'
      WRITE(*,'(A10,11A11)')   '   ', 'w','qv','qc','qr','qi','qs','qg','qh','nc','ni','tmax'
      WRITE(*,'(A10,11ES11.3)') '   ', (xmax(i), i=1,11)
    END IF
  ENDIF

  !..return now, if no clouds are found
  IF (.NOT.lb_imbalance .AND. igridpoints == -1) THEN 

    IF (debug) WRITE (*,*) TRIM(yzroutine)//": no clouds found!"

#ifdef COSMOART
    DEALLOCATE(calc_i,calc_j,calc_k)
#endif


    DEALLOCATE(lb_proc_num_cgp)
    IF (debug) WRITE (*,*) TRIM(yzroutine)//": deallocated"

  ELSE ! cloudy points have been found

    IF (isIO() .AND. MPI_LOAD_BALANCING) THEN
      IF (lb_imbalance) THEN
        WRITE (*,'(a)') '   '//TRIM(yzroutine)//' *** WITH *** load balancing ...'
      ELSE
        WRITE (*,'(a)') '   '//TRIM(yzroutine)//' *** NO *** load balancing ...'
      END IF
    END IF

    ! ... ... reciprocal vertical grid
    DO kk = kts, kte
      DO jj = jts, jte
        DO ii = its, ite           
          adz(ii,jj,kk) = 1.0/(hhl(ii,jj,kk) - hhl(ii,jj,kk+1))   ! reciprocal vertical grid
          zml(ii,jj,kk) = 0.5*(hhl(ii,jj,kk) + hhl(ii,jj,kk+1))   ! model levels
        ENDDO
      ENDDO
    ENDDO

    IF (lb_imbalance) THEN
#ifdef MPI_CLOUD
      ! ... LOAD BALANCING
      allocate(lb_proc_num_bgp(MPI_PROCS))
      allocate(lb_pp_num_bgp(MPI_PROCS,MPI_PROCS))
      ALLOCATE(lb_cgp_to_proc(lb_total_num_cgp), stat=izstat)
      IF (izstat /= 0) THEN
        WRITE (*,*) lb_total_num_cgp
        CALL model_abort (my_world_id, 11125, 'SEIFERT: Allocation error, lb_cgp_to_proc', &
             'SEIFERT: Allocation error')      
      END IF
      ! ... Round robin: determines the communication pattern
      j = 0
      lb_pp_num_bgp = 0.0
      IF (.TRUE.) THEN
        IF (MPI_debug) print *,' robin 2, mpi_myrank =',MPI_MYRANK
        ALLOCATE(isendbuf(MPI_PROCS))
        ALLOCATE(irecvbuf(MPI_PROCS*MPI_PROCS))
        isendbuf = 0
        k = MPI_MYRANK+1
        IF (MPI_debug) print *,'lb_proc_num_cgp(k) = ',lb_proc_num_cgp(k)
        DO i=1,lb_proc_num_cgp(k)
          j = j+1
          isendbuf(j) = isendbuf(j)+1  ! MY PART OF THE COMMUNICATION MATRIX
          if (j == MPI_PROCS) j=0
        ENDDO
        IF (MPI_debug) print *,'lb_pp_num_bgp (isendbuffer)= ',isendbuf
        call mpi_allgather(isendbuf,MPI_PROCS,mpi_integer,irecvbuf,&
          &                MPI_PROCS,mpi_integer,mpi_comm,MPI_err)
        i=1
        IF (MPI_debug) print *,'lb_pp_num_bgp (irecvbuffer) = ',irecvbuf
        do k=1,MPI_PROCS
          do j=1,MPI_PROCS
            lb_pp_num_bgp(k,j) = irecvbuf(i)
            i=i+1
          enddo
        enddo
        IF (MPI_debug) print *,'lb_pp_num_bgp(k) = ',lb_pp_num_bgp(k,:)
        deallocate(isendbuf)
        deallocate(irecvbuf)
      ENDIF
      IF (MPI_debug) print *,'lb_avg_num_cgp = ',lb_avg_num_cgp
      IF (MPI_debug) print *,'lb_pp_num_bgp = '
      IF (MPI_PROCS.le.10) THEN
        do i=1,MPI_PROCS
          IF (MPI_debug) WRITE (*,'(10i7)') lb_pp_num_bgp(i,:)
        enddo
      ENDIF
      ! ... Built send-recv buffers for all2all communication
      lb_arrays_per_gp = 22
      allocate(sendcnt(0:MPI_PROCS-1))
      allocate(recvcnt(0:MPI_PROCS-1))
      allocate(sdispl(0:MPI_PROCS-1))
      allocate(rdispl(0:MPI_PROCS-1))
      do j=1,MPI_PROCS     
        lb_proc_num_bgp(j) = sum(lb_pp_num_bgp(:,j)) 
        sendcnt(j-1) = lb_arrays_per_gp * lb_pp_num_bgp(MPI_MYRANK+1,j)
        recvcnt(j-1) = lb_arrays_per_gp * lb_pp_num_bgp(j,MPI_MYRANK+1)
      enddo
      IF (MPI_debug) print *,'lb_proc_num_cgp = ',lb_proc_num_cgp
      IF (MPI_debug) print *,'lb_proc_num_bgp = ',lb_proc_num_bgp
      IF (MPI_debug) print *,'lb_total_num_cgp = ',sum(lb_proc_num_cgp)
      IF (MPI_debug) print *,'lb_total_num_bgp = ',sum(lb_proc_num_bgp)
      IF (MPI_debug) WRITE (*,'(A,20i8)') ' sendcnt = ',sendcnt
      IF (MPI_debug) WRITE (*,'(A,20i8)') ' recvcnt = ',recvcnt
      if (sum(lb_proc_num_cgp).ne.sum(lb_proc_num_bgp)) then
        WRITE (*,'(A,4i8)') ' lb error stop'
        stop
      endif
      sdispl(0) = 0
      rdispl(0) = 0
      do j=2,MPI_PROCS
        sdispl(j-1) = sdispl(j-2)+sendcnt(j-2)
        rdispl(j-1) = rdispl(j-2)+recvcnt(j-2)
      enddo
      IF (MPI_debug) WRITE (*,'(A,20i8)') ' sdispl  =  ',sdispl
      IF (MPI_debug) WRITE (*,'(A,20i8)') ' rdispl  =  ',rdispl
      lb_sendlength = lb_arrays_per_gp*lb_proc_num_cgp(MPI_MYRANK+1)
      lb_recvlength = lb_arrays_per_gp*lb_proc_num_bgp(MPI_MYRANK+1)
      ALLOCATE(sendbuf(lb_sendlength), stat=izstat)
      IF (izstat /= 0) THEN
        WRITE (*,*) lb_sendlength
        CALL model_abort (my_world_id, 11125, 'SEIFERT: Allocation error, lb_sendlength', &
             'SEIFERT: Allocation error')      
      END IF
      ALLOCATE(recvbuf(lb_recvlength), stat=izstat)
      IF (izstat /= 0) THEN
        WRITE (*,*) lb_recvlength
        CALL model_abort (my_world_id, 11125, 'SEIFERT: Allocation error, lb_recvlength', &
             'SEIFERT: Allocation error')      
      END IF
      IF (MPI_debug) WRITE (*,'(A,20i8)') ' slength =  ',lb_sendlength
      IF (MPI_debug) WRITE (*,'(A,20i8)') ' rlength =  ',lb_recvlength
      ! ... fill sendbuffer with variables used in cloud module     
      if(lb_proc_num_cgp(MPI_MYRANK+1).ne.0) then
!!! ub>> alt        j = 1
        do i=1,lb_proc_num_cgp(MPI_MYRANK+1)
!!! ub>> neu
          j = (i-1) * lb_arrays_per_gp + 1
          ! ... these are the local WRF gripoints
          ii = ilm(i-1)
          jj = jlm(i-1)
          kk = klm(i-1)
          ! ... copy cloud stuff to sendbuffer
          sendbuf(j)   = qv(ii,jj,kk)        ! qv  -> q
          sendbuf(j+1) = qc(ii,jj,kk)        ! qc  -> q_cloud
          sendbuf(j+2) = qr(ii,jj,kk)        ! qr  -> q_rain
          sendbuf(j+3) = qi(ii,jj,kk)        ! qi  -> q_ice
          sendbuf(j+4) = qs(ii,jj,kk)        ! qs  -> q_snow
          sendbuf(j+5) = qg(ii,jj,kk)        ! qg  -> q_graupel
          sendbuf(j+6) = qh(ii,jj,kk)        ! qh  -> q_hail
          sendbuf(j+7) = qnc(ii,jj,kk)       ! qnc -> n_cloud
          sendbuf(j+8) = qnr(ii,jj,kk)       ! qnr -> n_rain
          sendbuf(j+9) = qni(ii,jj,kk)       ! qni -> n_ice
          sendbuf(j+10) = qns(ii,jj,kk)       ! qns -> n_snow
          sendbuf(j+11) = qng(ii,jj,kk)      ! qng -> n_graupel           
          sendbuf(j+12) = qnh(ii,jj,kk)      ! qnh -> n_hail
          sendbuf(j+13) = T(ii,jj,kk,nx)        ! Th*pii -> T_2mom
          sendbuf(j+14) = p0(ii,jj,kk)+pp(ii,jj,kk,nx)                ! p -> p_2mom    
          sendbuf(j+15) = rho(ii,jj,kk)                ! rho -> rho_2mom    
          sendbuf(j+16) = (ssw(ii,jj,MAX(kk-1,kts))-ssw(ii,jj,MAX(kk,kts))) * adz(ii,jj,kk)     &
               &             *MAX(w(ii,jj,kk,nx),w(ii,jj,MAX(kk-1,kts),nx)) ! -> dSwdz
          sendbuf(j+17) = (ssi(ii,jj,MAX(kk-1,kts))-ssi(ii,jj,MAX(kk,kts))) * adz(ii,jj,kk)     &
               &             *MAX(w(ii,jj,kk,nx),w(ii,jj,MAX(kk-1,kts),nx)) ! -> dSidz
          sendbuf(j+18) = zml(ii,jj,kk)        ! zml -> zml_k

! detect cloud base and assign cloud base updraft w_cb, incloud nucleation is possible depending on switch "lincloud":  

!MB is loadbalancing possible in case of aerosol interaction??

          IF (l2mom_satads) THEN

            ! In this case, we cannot use grid scale supersaturation for searching the cloud base,
            ! because there is no grid scale supersaturation! So we use qc (mass specific quantity)
            ! as precursor for previously supersaturated conditions:

          
            kkp1 = MIN(kk+1,kte)

            ! UNTESTED!
            IF ( lincloud .AND. (w(ii,jj,kkp1,nx) > wcb_min .AND. qc(ii,jj,kk) >= eps &
                 & .AND.  qc(ii,jj,kk) > qc(ii,jj,kkp1)))  THEN
              sendbuf(j+19) = w(ii,jj,kkp1,nx)
            ELSE IF (.NOT.lincloud .AND. (w(ii,jj,kkp1,nx) > wcb_min .AND. qc(ii,jj,kk) >= eps .AND. &
                 qc(ii,jj,kkp1) < eps)) THEN
              sendbuf(j+19) = w(ii,jj,kkp1,nx)
            ELSE
              sendbuf(j+19) = 0.0
            ENDIF

          ELSE

            ! In this case there is grid scale supersaturation and we can use it
            ! to find the cloud base:

            kkp1 = MIN(kk+1,kte)
          
            IF ( lincloud .AND. (w(ii,jj,kkp1,nx) > wcb_min .AND. ssw(ii,jj,kk) >= scb_min &
                 & .AND.  ssw(ii,jj,kk) > ssw(ii,jj,kkp1)))  THEN
              sendbuf(j+19) = w(ii,jj,kkp1,nx)
            ELSE IF (.NOT.lincloud .AND. (w(ii,jj,kkp1,nx) > wcb_min .AND. ssw(ii,jj,kk) >= scb_min .AND. &
                 ssw(ii,jj,kkp1) < scb_max)) THEN
              sendbuf(j+19) = w(ii,jj,kkp1,nx)
            ELSE
              sendbuf(j+19) = 0.0
            ENDIF

          END IF


          sendbuf(j+20) = ( (T(ii,jj,MAX(kk-1,kts),nx)-T(ii,jj,kk,nx)) * adz(ii,jj,kk)     &
               &             * MAX(w(ii,jj,kk,nx),w(ii,jj,MAX(kk-1,kts),nx)) ) ! -> w
          sendbuf(j+21) = w(ii,jj,kk,nx)       ! w -> w_2mom



!!! ub>> alt          j = j+lb_arrays_per_gp
        enddo
      else
        sendbuf = 0.0
      endif
      ! ... do ALL2ALL communication
      !WRITE (*,*) ' MPI_ALLTOALLV 1, send = ',lb_proc_num_cgp(MPI_MYRANK+1)
      CALL MPI_ALLTOALLV(sendbuf,sendcnt,sdispl,mpi_real,&
           &             recvbuf,recvcnt,rdispl,mpi_real,mpi_comm,mpi_err)
      !WRITE (*,*) ' MPI_ALLTOALLV 1, recv = ',lb_proc_num_bgp(MPI_MYRANK+1)
      ! ... copy cloud variables from recvbuffer into cloud arrays
      IF(lb_proc_num_bgp(MPI_MYRANK+1) /= 0) THEN
        ! .. define grid for cloud module
        loc_ix = lb_proc_num_bgp(MPI_MYRANK+1)-1
        loc_iy = 1
        loc_iz = 1
        ! ... copy the recvbuffer and convert to densities
!!! ub>> alt        j = 1
        do i=0,loc_ix
!!! ub>> neu
          j = i * lb_arrays_per_gp + 1
          T_2mom(i,1,1)       = DBLE(recvbuf(j+13))  ! Th*pii -> T_2mom
          p_2mom(i,1,1)       = DBLE(recvbuf(j+14))  ! p   -> p_2mom    
          x_v              = DBLE(recvbuf(j)) 
          rho_2mom(i,1,1)     = DBLE(recvbuf(j+15))
          hlp              = rho_2mom(i,1,1)
          q(i,1,1)         = hlp*DBLE(recvbuf(j))    ! qv  -> q
          q_cloud(i,1,1)   = hlp*DBLE(recvbuf(j+1))  ! qc  -> q_cloud
          q_rain(i,1,1)    = hlp*DBLE(recvbuf(j+2))  ! qr  -> q_rain
          q_ice(i,1,1)     = hlp*DBLE(recvbuf(j+3))  ! qi  -> q_ice
          q_snow(i,1,1)    = hlp*DBLE(recvbuf(j+4))  ! qs  -> q_snow
          q_graupel(i,1,1) = hlp*DBLE(recvbuf(j+5))  ! qg  -> q_graupel
          q_hail(i,1,1)    = hlp*DBLE(recvbuf(j+6))  ! qh  -> q_hail
          n_cloud(i,1,1)   = hlp*DBLE(recvbuf(j+7))  ! qnc -> n_cloud
          n_rain(i,1,1)    = hlp*DBLE(recvbuf(j+8))  ! qnr -> n_rain
          n_ice(i,1,1)     = hlp*DBLE(recvbuf(j+9))  ! qni -> n_ice
          n_snow(i,1,1)    = hlp*DBLE(recvbuf(j+10)) ! qns -> n_snow
          n_graupel(i,1,1) = hlp*DBLE(recvbuf(j+11)) ! qng -> n_graupel
          n_hail(i,1,1)    = hlp*DBLE(recvbuf(j+12)) ! qnh -> n_snow

! Berechnung des Dampfdruckes unter Vernachl. des Covolumens der Hydrometeore:
! (Erklaerung: Bei Anwendung der Gasgleichung muesste eigentlich die abs. Feuchte
! als Verh. von Masse Dampf zu Volumen Dampf+Luft stehen, so dass der Druck
! den Druck im Gasanteil darstellt; das hier verwendete q(i,1,1) ist aber
! Masse Dampf pro Volumen Luft+Dampf+Hydrometeore)
          e_v              = q(i,1,1)*T_2mom(i,1,1)*r_v
          S_w(i,1,1)       = e_v/e_ws(T_2mom(i,1,1)) - 1.0  
          S_i(i,1,1)       = e_v/e_es(T_2mom(i,1,1)) - 1.0
          dSwdz(i,1,1)     = DBLE(recvbuf(j+16))     !     -> dSwdz
          dSidz(i,1,1)     = DBLE(recvbuf(j+17))     !     -> dSidz
          zml_k(i,1,1)     = DBLE(recvbuf(j+18))     !     -> zml
          w_cb(i,1,1)      = DBLE(recvbuf(j+19))     !     -> w_cb
          dT0dz(i,1,1)     = DBLE(recvbuf(j+20))     !     -> dT0dz
          w_2mom(i,1,1)   = DBLE(recvbuf(j+21))     !     -> w_2mom

!!! ub>> alt          j = j+lb_arrays_per_gp
        enddo
      endif

      DEALLOCATE(sendbuf, recvbuf)

#endif
    ELSE
      ! ... NO LOAD BALANCING -> COPY ARRAYS DIRECTLY
      ! ... Grid
      loc_ix = igridpoints
      loc_iy = 1
      loc_iz = 1
      IF (debug.and.isIO()) THEN
        WRITE (*,*) TRIM(yzroutine)//": grid"
        WRITE (*,*) "      its = ",its
        WRITE (*,*) "      ite = ",ite
        WRITE (*,*) "      jts = ",jts
        WRITE (*,*) "      jte = ",jte
        WRITE (*,*) "      kts = ",kts
        WRITE (*,*) "      kte = ",kte
        WRITE (*,*) "      loc_ix = ",loc_ix
        WRITE (*,*) "      loc_iy = ",loc_iy
        WRITE (*,*) "      loc_iz = ",loc_iz
      ENDIF
      IF (debug.and.isIO()) WRITE (*,*) TRIM(yzroutine)//": trafo 1, no load balancing", loc_ix
      ! ... transpose to one-dimensional array and variables used in cloud module
      j = 1
      k = 1
!CDIR    NODEP                                      
      DO i=0,loc_ix
        ! ... WRF grid points
        ii = ilm(i)
        jj = jlm(i)
        kk = klm(i)
        ! ... dynamics
        T_2mom(i,j,k)      = T(ii,jj,kk,nx)
        p_2mom(i,j,k)      = p0(ii,jj,kk)+pp(ii,jj,kk,nx)
        ! ... use own moist density
        x_v             = qv(ii,jj,kk)
        rho_2mom(i,j,k)    = rho(ii,jj,kk)
        ! .. the supersaturations and their vertical gradients * vertical velocity
        S_w(i,j,k)   = ssw(ii,jj,kk) 
        S_i(i,j,k)   = ssi(ii,jj,kk)
        dSwdz(i,j,k) = (ssw(ii,jj,max(kk-1,kts))-ssw(ii,jj,max(kk,kts))) * adz(ii,jj,kk)     &
             &             *MAX(w(ii,jj,kk,nx),w(ii,jj,max(kk-1,kts),nx))
        dSidz(i,j,k) = (ssi(ii,jj,max(kk-1,kts))-ssi(ii,jj,max(kk,kts))) * adz(ii,jj,kk)     &
             &             *MAX(w(ii,jj,kk,nx),w(ii,jj,max(kk-1,kts),nx))
        ! ... mass concentrations --> number densities
        n_cloud(i,j,k)   = rho_2mom(i,j,k) * DBLE( qnc(ii,jj,kk) )
        n_rain(i,j,k)    = rho_2mom(i,j,k) * DBLE( qnr(ii,jj,kk) )
        n_ice(i,j,k)     = rho_2mom(i,j,k) * DBLE( qni(ii,jj,kk) )
        n_snow(i,j,k)    = rho_2mom(i,j,k) * DBLE( qns(ii,jj,kk) )
        n_graupel(i,j,k) = rho_2mom(i,j,k) * DBLE( qng(ii,jj,kk) )
        n_hail(i,j,k)    = rho_2mom(i,j,k) * DBLE( qnh(ii,jj,kk) )

        ! ... mixing ratios -> mass densities
        q(i,j,k)         = rho_2mom(i,j,k) * DBLE( x_v )
        q_cloud(i,j,k)   = rho_2mom(i,j,k) * DBLE( qc(ii,jj,kk) )
        q_rain(i,j,k)    = rho_2mom(i,j,k) * DBLE( qr(ii,jj,kk) ) 
        q_ice(i,j,k)     = rho_2mom(i,j,k) * DBLE( qi(ii,jj,kk) )
        q_snow(i,j,k)    = rho_2mom(i,j,k) * DBLE( qs(ii,jj,kk) )
        q_graupel(i,j,k) = rho_2mom(i,j,k) * DBLE( qg(ii,jj,kk) )           
        q_hail(i,j,k)    = rho_2mom(i,j,k) * DBLE( qh(ii,jj,kk) )           

! detect cloud base and assign cloud base updraft w_cb, incloud nucleation is possible depending on switch "lincloud":
!MB Need to implement a cloud diagnosis e.g. cloudact=1 new cloud; cloudact=2 cloud base; cloudact=3 incloud 

        IF (l2mom_satads) THEN

          ! In this case, we cannot use grid scale supersaturation for searching the cloud base,
          ! because there is no grid scale supersaturation! So we use qc (mass specific quantity)
          ! as precursor for previously supersaturated conditions:
          
          kkp1 = MIN(kk+1,kte)

          ! UNTESTED!
          IF ( lincloud .AND. (w(ii,jj,kkp1,nx) > wcb_min .AND. qc(ii,jj,kk) >= eps &
               & .AND.  qc(ii,jj,kk) > qc(ii,jj,kkp1)))  THEN

            w_cb(i,j,k) = w(ii,jj,kkp1,nx)

#ifdef COSMOART
            IF(qc(ii,jj,kkp1) < eps)    THEN   ! cloud base
              cloudact(i,j,k)=2 
            ELSE IF(qc(ii,jj,kk) > qc(ii,jj,kkp1)) THEN  ! incloud activation
              cloudact(i,j,k)=3
            ENDIF
#endif

          ELSE IF (.NOT.lincloud .AND. (w(ii,jj,kkp1,nx) > wcb_min .AND. qc(ii,jj,kk) >= eps .AND. &
               qc(ii,jj,kkp1) < eps)) THEN
            
            w_cb(i,j,k) = w(ii,jj,kkp1,nx)
            
#ifdef COSMOART
            cloudact(i,j,k)=2 
#endif
            
          ELSE
            w_cb(i,j,k) = 0.0
          ENDIF

        ELSE

          ! In this case there is grid scale supersaturation and we can use it
          ! to find the cloud base:

          kkp1 = MIN(kk+1,kte)

          IF ( lincloud .AND. (w(ii,jj,kkp1,nx) > wcb_min .AND. ssw(ii,jj,kk) >= scb_min &
               & .AND.  ssw(ii,jj,kk) > ssw(ii,jj,kkp1)))  THEN

            w_cb(i,j,k) = w(ii,jj,kkp1,nx)

#ifdef COSMOART
            IF(ssw(ii,jj,kkp1) < scb_max)    THEN   ! cloud base
              cloudact(i,j,k)=2 
            ELSE IF(ssw(ii,jj,kk) > ssw(ii,jj,kkp1)) THEN  ! incloud activation
              cloudact(i,j,k)=3
            ENDIF
#endif

          ELSE IF (.NOT.lincloud .AND. (w(ii,jj,kkp1,nx) > wcb_min .AND. ssw(ii,jj,kk) >= scb_min .AND. &
               ssw(ii,jj,kkp1) < scb_max)) THEN
            
            w_cb(i,j,k) = w(ii,jj,kkp1,nx)
            
#ifdef COSMOART
            cloudact(i,j,k)=2                
#endif
            
          ELSE
            w_cb(i,j,k) = 0.0
          ENDIF
          
        END IF

#ifdef COSMOART
        !IF  n_cloud 'too low' do activation as if new
        IF( n_cloud(i,j,k) <= 10.0e6 )          THEN ! new cloud
          cloudact(i,j,k)=1
        ENDIF
#endif
        zml_k(i,j,k) = zml(ii,jj,kk) 

        dT0dz(i,j,k) = ( (T(ii,jj,MAX(kk-1,kts),nx)-T(ii,jj,kk,nx)) * adz(ii,jj,kk)     &
             &             * MAX(w(ii,jj,kk,nx),w(ii,jj,MAX(kk-1,kts),nx)) ) ! -> dT0dz

        w_2mom(i,j,k)   = w(ii,jj,kk,nx)

#ifdef COSMOART
        calc_i(i)=ii
        calc_j(i)=jj
        calc_k(i)=kk
#endif
      enddo
    ENDIF

    ! ... timestep
    dt_2mom = zdt

    IF (debug.and.isIO()) WRITE (*,*) TRIM(yzroutine)//": calling clouds"

    IF (ltime) CALL get_timings (i_2mom_todens, ntstep, dt, izerror)

    ! .. this subroutine calculates all the microphysical sources and sinks
    CALL clouds ()

    IF (ltime) CALL get_timings (i_2mom_clouds, ntstep, dt, izerror)

    IF (debug.and.isIO()) WRITE (*,*) TRIM(yzroutine)//": beeing back from clouds"

    ! ... output coeffs and density correction for fall velocity of precip particles
    IF (isIO().AND.firstcall.NE.1.AND.debug) THEN
      WRITE (*,*) TRIM(yzroutine)//": particles and PSD shapes"
      WRITE (*,*) "          type = ", cloud%name  
      WRITE (*,*) "          nu   = ", cloud%nu
      WRITE (*,*) "          mu   = ", cloud%mu
      WRITE (*,*) "          type = ", rain%name  
      WRITE (*,*) "          nu   = ", rain%nu
      WRITE (*,*) "          mu   = ", rain%mu
      WRITE (*,*) "          type = ", ice%name  
      WRITE (*,*) "          nu   = ", ice%nu
      WRITE (*,*) "          mu   = ", ice%mu
      WRITE (*,*) "          type = ", snow%name  
      WRITE (*,*) "          nu   = ", snow%nu
      WRITE (*,*) "          mu   = ", snow%mu
      WRITE (*,*) "          type = ", graupel%name  
      WRITE (*,*) "          nu   = ", graupel%nu
      WRITE (*,*) "          mu   = ", graupel%mu
      WRITE (*,*) "          type = ", hail%name  
      WRITE (*,*) "          nu   = ", hail%nu
      WRITE (*,*) "          mu   = ", hail%mu
      WRITE (*,*) 
      firstcall = 1     
    ENDIF

    IF (debug) THEN
      IF (MINVAL(q) < 0.0) THEN
        WRITE (*,*) ' SEIFERT: q < 0 ENCOUNTERED, FILLED WITH 0.0'
        WHERE (q < 0.0) q = 0.0
      ENDIF
      IF (MINVAL(q_cloud) < 0.0) THEN
        write (*,*) ' SEIFERT: q_cloud < 0, STOPPED AFTER CLOUDS 1'
        stop
      ENDIF
      IF (MINVAL(q_rain) < 0.0) THEN
        write (*,*) ' SEIFERT: q_rain < 0, STOPPED AFTER CLOUDS 1'
        stop
      ENDIF
      IF (MINVAL(q_ice) < 0.0) THEN
        write (*,*) ' SEIFERT: q_ice < 0, STOPPED AFTER CLOUDS 1'
        stop
      ENDIF
      IF (MINVAL(q_snow) < 0.0) THEN
        write (*,*) ' SEIFERT: q_snow < 0, STOPPED AFTER CLOUDS 1'
        stop
      ENDIF
      IF (MINVAL(q_graupel) < 0.0) THEN
        WRITE (*,*) ' SEIFERT: q_graupel < 0, STOPPED AFTER CLOUDS 1'
        stop
      ENDIF
      IF (MINVAL(q_hail) < 0.0) THEN
        WRITE (*,*) ' SEIFERT: q_hail < 0, STOPPED AFTER CLOUDS 1'
        stop
      ENDIF
      IF (MINVAL(n_cloud) < 0.0) THEN
        write (*,*) ' SEIFERT: n_cloud < 0, STOPPED AFTER CLOUDS 1'
        stop
      ENDIF
      IF (MINVAL(n_rain) < 0.0) THEN
        write (*,*) ' SEIFERT: n_rain < 0, STOPPED AFTER CLOUDS 1'
        stop
      ENDIF
      IF (MINVAL(n_ice) < 0.0) THEN
        write (*,*) ' SEIFERT: n_ice < 0, STOPPED AFTER CLOUDS 1'
        stop
      ENDIF
      IF (MINVAL(n_snow) < 0.0) THEN
        write (*,*) ' SEIFERT: n_snow < 0, STOPPED AFTER CLOUDS 1'
        stop
      ENDIF
      IF (MINVAL(n_graupel) < 0.0) THEN
        write (*,*) ' SEIFERT: n_graupel < 0, STOPPED AFTER CLOUDS 1'
        stop
      ENDIF
      IF (MINVAL(n_hail) < 0.0) THEN
        write (*,*) ' SEIFERT: n_hail < 0, STOPPED AFTER CLOUDS 1'
        stop
      ENDIF
    ENDIF

    ! ... Transformation of variables back to driving model and latent heat equation
    IF (lb_imbalance) THEN
#ifdef MPI_CLOUD
      ! ... copy cloud variables from cloud arrays to recvbuffer which is now the sendbuffer
      ! ... Built send-recv buffers for all2all communication

      lb_arrays_per_gp = 14

      IF (speichere_precipstat) lb_arrays_per_gp = lb_arrays_per_gp + 2
      lb_sendlength = lb_arrays_per_gp*lb_proc_num_cgp(MPI_MYRANK+1)
      lb_recvlength = lb_arrays_per_gp*lb_proc_num_bgp(MPI_MYRANK+1)
      ALLOCATE(sendbuf(lb_sendlength), stat=izstat)
      IF (izstat /= 0) THEN
        WRITE (*,*) lb_sendlength
        CALL model_abort (my_world_id, 11125, 'SEIFERT: Allocation error, lb_sendlength', &
             'SEIFERT: Allocation error')      
      END IF
      ALLOCATE(recvbuf(lb_recvlength), stat=izstat)
      IF (izstat /= 0) THEN
        WRITE (*,*) lb_recvlength
        CALL model_abort (my_world_id, 11125, 'SEIFERT: Allocation error, lb_recvlength', &
             'SEIFERT: Allocation error')      
      END IF
      
      DO j=1,MPI_PROCS     
        sendcnt(j-1) = lb_arrays_per_gp * lb_pp_num_bgp(MPI_MYRANK+1,j)
        recvcnt(j-1) = lb_arrays_per_gp * lb_pp_num_bgp(j,MPI_MYRANK+1)
      ENDDO
      IF (MPI_debug) WRITE (*,'(A,20i8)') ' sendcnt = ',sendcnt
      IF (MPI_debug) WRITE (*,'(A,20i8)') ' recvcnt = ',recvcnt
      sdispl(0) = 0
      rdispl(0) = 0
      DO j=2,MPI_PROCS
        sdispl(j-1) = sdispl(j-2)+sendcnt(j-2)
        rdispl(j-1) = rdispl(j-2)+recvcnt(j-2)
      ENDDO
      IF (MPI_debug) WRITE (*,'(A,20i8)') ' sdispl  =  ',sdispl
      IF (MPI_debug) WRITE (*,'(A,20i8)') ' rdispl  =  ',rdispl
      IF (MPI_debug) WRITE (*,'(A,20i8)') ' slength =  ',lb_sendlength
      IF (MPI_debug) WRITE (*,'(A,20i8)') ' rlength =  ',lb_recvlength
      ! ... copy into recvbuffer and convert to WRF variables
      IF(lb_proc_num_bgp(MPI_MYRANK+1).NE.0) THEN
!!! ub>> alt        j = 1
!CDIR    NODEP                                      
        DO i=0,loc_ix
!!! ub>> neu
          j = i * lb_arrays_per_gp + 1
          ! ... vapor mixing ratio
          recvbuf(j)    = q(i,1,1)
          recvbuf(j+1)  = q_cloud(i,1,1)   ! qc  -> q_cloud
          recvbuf(j+2)  = q_rain(i,1,1)    ! qr  -> q_rain
          recvbuf(j+3)  = q_ice(i,1,1)     ! qi  -> q_ice
          recvbuf(j+4)  = q_snow(i,1,1)    ! qs  -> q_snow
          recvbuf(j+5)  = q_graupel(i,1,1) ! qg  -> q_graupel
          recvbuf(j+6)  = q_hail(i,1,1)    ! qh  -> q_hail
          recvbuf(j+7)  = n_cloud(i,1,1)         ! qnc -> n_cloud
          recvbuf(j+8)  = n_rain(i,1,1)          ! qnr -> n_rain
          recvbuf(j+9)  = n_ice(i,1,1)           ! qni -> n_ice
          recvbuf(j+10) = n_snow(i,1,1)          ! qns -> n_snow
          recvbuf(j+11) = n_graupel(i,1,1)       ! qng -> n_graupel
          recvbuf(j+12) = n_hail(i,1,1)       ! qnh -> n_hail
          recvbuf(j+13) = rho_2mom(i,1,1)                 ! rho_2mom


!!! ub>>          j = j+lb_arrays_per_gp
        ENDDO
      ELSE
        recvbuf = 0.0
      ENDIF
      ! ... the ALL2ALL communication: now backwards
      CALL MPI_ALLTOALLV(recvbuf,recvcnt,rdispl,mpi_real,&
           &             sendbuf,sendcnt,sdispl,mpi_real,mpi_comm,mpi_err)
      ! ... copy cloud variables from sendbuffer, which is now the recvbuffer, to WRF
      !     and calculate the latent heat release


      IF(lb_proc_num_cgp(MPI_MYRANK+1).NE.0) THEN

!!! ub>>        j = 1
!CDIR    NODEP                                      
        DO i=0,lb_proc_num_cgp(MPI_MYRANK+1)-1
!!! ub>> neu
          j = i * lb_arrays_per_gp + 1
          ! ... the old WRF grid points
          ii = ilm(i)
          jj = jlm(i)
          kk = klm(i)
          ! ... latent heat for temperature equation        

!!$          rho(ii,jj,kk) = sendbuf(j+13)
          hlp = 1.0 / rho(ii,jj,kk)
          q_vap_new = hlp * sendbuf(j)
          q_liq_new = hlp * (sendbuf(j+1)+sendbuf(j+2))
          q_ice_new = hlp * (sendbuf(j+3)+sendbuf(j+4)+sendbuf(j+5)+sendbuf(j+6))
          ! ... ... old moisture variables
          q_vap_old = qv(ii,jj,kk)
          q_liq_old = qc(ii,jj,kk) + qr(ii,jj,kk)
          q_ice_old = qi(ii,jj,kk) + qs(ii,jj,kk) &
               + qg(ii,jj,kk) + qh(ii,jj,kk)   

          ! ... temperature equation
          zlh_s = lh_s - zlhTfac * ( cv_i - r_v-cv_v ) * ( T(ii,jj,kk,nx) - t0_melt )
          zlh_f = lh_s - lh_v - zlhTfac * (cv_i-cv_l)  * ( T(ii,jj,kk,nx) - t0_melt )
          T(ii,jj,kk,nx) = T(ii,jj,kk,nx) - lewfaktor * zttfac * ( &
                       ( zlh_s - zlhTfac*r_v*T(ii,jj,kk,nx) ) * (q_vap_new - q_vap_old) &
               &     +   zlh_f                                  * (q_liq_new - q_liq_old) )


          ! ... pressure equation:
 ! UB_20111222>>
!!$ .. The pressure tendency terms
!!$     due to phase changes as they appear in the original unapproximated
!!$     pressure equation:

!!$ currently deactivated because of zpptfac = 0.0 above:
          pp(ii,jj,kk,nx) = pp(ii,jj,kk,nx) - rho(ii,jj,kk) * zpptfac * ( &
               (q_liq_new - q_liq_old) * ( r_d*zlh_f )                                    + & 
               (q_vap_new - q_vap_old) * ( r_d*zlh_s - zlhTfac*cp_d*r_v*T(ii,jj,kk,nx) ) )
               
! UB_20111222<<

          qv(ii,jj,kk) = hlp * sendbuf(j)
          qc(ii,jj,kk) = hlp * sendbuf(j+1)
          qr(ii,jj,kk) = hlp * sendbuf(j+2)
          qi(ii,jj,kk) = hlp * sendbuf(j+3)
          qs(ii,jj,kk) = hlp * sendbuf(j+4)
          qg(ii,jj,kk) = hlp * sendbuf(j+5)
          qh(ii,jj,kk) = hlp * sendbuf(j+6)
          ! ... number concentrations
          qnc(ii,jj,kk) = hlp * sendbuf(j+7)
          qnr(ii,jj,kk) = hlp * sendbuf(j+8)
          qni(ii,jj,kk) = hlp * sendbuf(j+9)
          qns(ii,jj,kk) = hlp * sendbuf(j+10)
          qng(ii,jj,kk) = hlp * sendbuf(j+11)
          qnh(ii,jj,kk) = hlp * sendbuf(j+12)


        ENDDO

      END IF



      DEALLOCATE(sendbuf, recvbuf)

#endif


    ELSE

      ! NO LOAD BALANCING
      IF (debug.and.isIO()) WRITE (*,*) TRIM(yzroutine)//": trafo 2, no load balancing"

      j = 1
      k = 1
!CDIR    NODEP                                      
      DO i=0,loc_ix
        ii = ilm(i)
        jj = jlm(i)
        kk = klm(i)

!!$        rho(ii,jj,kk) = rho_2mom(i,j,k)
        hlp = 1.0 / rho(ii,jj,kk)
        q_v = hlp * q(i,j,k)

        ! ... latent heat for temperature equation        
        ! ... ... new variables
        q_vap_new = q_v
        q_liq_new = hlp * ( q_cloud(i,j,k)+q_rain(i,j,k) )
        q_ice_new = hlp * ( q_ice(i,j,k)+q_snow(i,j,k)+q_graupel(i,j,k)+q_hail(i,j,k) )
        ! ... ... old variables
        q_vap_old = qv(ii,jj,kk)
        q_liq_old = qc(ii,jj,kk) + qr(ii,jj,kk)
        q_ice_old = qi(ii,jj,kk) + qs(ii,jj,kk) + qg(ii,jj,kk) + qh(ii,jj,kk)

        ! ... temperature equation
        zlh_s = lh_s - zlhTfac * ( cv_i - r_v-cv_v ) * ( T(ii,jj,kk,nx) - t0_melt )
        zlh_f = lh_s - lh_v - zlhTfac * (cv_i-cv_l)  * ( T(ii,jj,kk,nx) - t0_melt )
        T(ii,jj,kk,nx) = T(ii,jj,kk,nx) - lewfaktor * zttfac * ( &
                       ( zlh_s - zlhTfac*r_v*T(ii,jj,kk,nx) ) * (q_vap_new - q_vap_old) &
                     +   zlh_f                                  * (q_liq_new - q_liq_old) )



          ! ... pressure equation:
 ! UB_20111222>>
!!$ .. The pressure tendency terms
!!$     due to phase changes as they appear in the original unapproximated
!!$     pressure equation:

!!$ currently deactivated because of zpptfac = 0.0 above:
        pp(ii,jj,kk,nx) = pp(ii,jj,kk,nx) - rho(ii,jj,kk) * zpptfac * ( &
             (q_liq_new - q_liq_old) * ( r_d*zlh_f )                                    + & 
             (q_vap_new - q_vap_old) * ( r_d*zlh_s - zlhTfac*cp_d*r_v*T(ii,jj,kk,nx) ) )
               
! UB_20111222<<

        qv(ii,jj,kk) = hlp * q(i,j,k)
        qc(ii,jj,kk) = hlp * q_cloud(i,j,k)
        qr(ii,jj,kk) = hlp * q_rain(i,j,k)
        qi(ii,jj,kk) = hlp * q_ice(i,j,k)
        qs(ii,jj,kk) = hlp * q_snow(i,j,k)
        qg(ii,jj,kk) = hlp * q_graupel(i,j,k)
        qh(ii,jj,kk) = hlp * q_hail(i,j,k)
        ! ... mass concentrations
        qnc(ii,jj,kk) = hlp * n_cloud(i,j,k)
        qnr(ii,jj,kk) = hlp * n_rain(i,j,k)
        qni(ii,jj,kk) = hlp * n_ice(i,j,k)
        qns(ii,jj,kk) = hlp * n_snow(i,j,k)
        qng(ii,jj,kk) = hlp * n_graupel(i,j,k)
        qnh(ii,jj,kk) = hlp * n_hail(i,j,k)


      ENDDO



    ENDIF   ! no load balancing

#ifdef COSMOART
    deallocate(calc_i,calc_j,calc_k)    
#endif

    ! Final saturation adjustment using satad from meteo_utilities at all grid points,
    ! because the saturation_adjust_h2o_qv above might have been only on 
    ! cloudy gridpoints:
    DO k=1,ke

      ! Do a final saturation adjustment for new values of t, qv and qc
      zpres(:,:) = p0(:,:,k) + pp(:,:,k,nx)
      
      zdummy(:,:,9) = t(:,:,k,nx)

      CALL satad ( 1, t(:,:,k,nx), qv(:,:,k),          &
           qc(:,:,k), zdummy(:,:,9), zpres,      &
           zdummy(:,:,1),zdummy(:,:,2),zdummy(:,:,3), &
           zdummy(:,:,4),zdummy(:,:,5),zdummy(:,:,6), &
           zdummy(:,:,7),zdummy(:,:,8),               &
           b1, b2w, b3, b4w, b234w, rdv, o_m_rdv,     &
           rvd_m_o, lh_v, cpdr, cp_d,                 &
           ie, je, istartpar, iendpar, jstartpar, jendpar)

    END DO

    ! Temperature, qv and qc has been changed, so rediagnose density:
    DO kk = kts, kte
!CDIR COLLAPSE
      qrs(:,:,kk) = qr(:,:,kk) + &
           qi(:,:,kk) + qs(:,:,kk) + &
           qg(:,:,kk) + qh(:,:,kk)
    END DO

! UB_20100525>>  For now deactivated, but this will be the future saturation adjustment interface:
!!$    CALL satad_driver (calling_routine = TRIM(yzroutine), &
!!$         lisochoric = ldiab_isochoric,                           &
!!$         te=t(:,:,:,nx), qve=qv(:,:,:),              &
!!$         qce=qc(:,:,:), ppe=pp(:,:,:,nx),            &
!!$         qle=qr(:,:,:), qie=qrs(:,:,:)-qr(:,:,:), p0e=p0,    &
!!$         its=its, ite=ite, jts=jts, jte=jte, kts=kts, kte=kte,   &
!!$         ims=1, ime=ie, jms=1, jme=je, kms=1, kme=ke,                             &
!!$         b1=b1, b2w=b2w, b3=b3, b4w=b4w, b234w=b234w,              &
!!$         lwd=lh_v, liw=lh_f, r_d=r_d, r_v=r_v,  rvd_m_o=rvd_m_o,   &
!!$         cvd=cv_d, cvv=cv_v, cl=cv_l, ci=cv_i,                     &
!!$         rdv=rdv, o_m_rdv=o_m_rdv, cpdr=cpdr, cp_d=cp_d,           & ! for old satad
!!$         satad_count=satad_count, satad_errstat=satad_errstat)
!!$    
!!$    ! Check error status of saturation adjustment:
!!$    SELECT CASE(satad_errstat)
!!$    CASE (0)
!!$      ! no error occured, so do nothing ...
!!$    CASE (1)
!!$      ! Allocation error for working arrays, abort:
!!$      CALL model_abort(my_cart_id, 10101,                                    &
!!$           'SATAD_DRIVER '//TRIM(yzroutine)//': Error allocating work arrays', 'seifert_pp')
!!$    CASE (2)
!!$      ! Convergence problem, report a warning but continue model run:
!!$      WRITE (*,'(a)') 'SATAD_DRIVER '//TRIM(yzroutine)//': WARNING: ' // &
!!$           'Convergence problem at some gridpoint. CONTINUE WITH CAUTION!'
!!$    CASE default
!!$      ! other error states are currently not defined, do nothing
!!$    END SELECT
    

    CALL calrho ( t(:,:,:,nx), pp(:,:,:,nx), &
         qv(:,:,:), qc(:,:,:),&
         qrs(:,:,:), p0(:,:,:), &
         rho(:,:,:), ie, je, ke, r_d, rvd_m_o)
! UB_20090804<<

    IF (nuc_c_typ .EQ. 0) THEN
      IF (debug.and.isIO()) WRITE (*,*) '  ... force constant cloud droplet number conc. after SATAD'
! UB_20090901>> Bugfix: qnr --> qnc !
!CDIR COLLAPSE
      WHERE (qc(:,:,:) > 0.0d0) qnc(:,:,:) = qnc_const / rho(:,:,:)
    END IF

!!! Diese Abfragen mit any() in Feldform schreiben, damit vektorisieren!
    IF (debug) THEN
      DO kk = kts, kte
        DO jj = jts, jte
          DO ii = its, ite           
            if (qc(ii,jj,kk).lt.0) then
              write (*,*) ' SEIFERT: qc < 0, STOPPED AFTER CLOUDS AND RECOPYING TO LM-VARIABLES'
              CALL model_abort (my_world_id, 999, 'negative values in qc',TRIM(yzroutine))      
            endif
            if (qr(ii,jj,kk).lt.0) then
              write (*,*) ' SEIFERT: qr < 0, STOPPED AFTER CLOUDS AND RECOPYING TO LM-VARIABLES'
              CALL model_abort (my_world_id, 999, 'negative values in qr',TRIM(yzroutine))      
            endif
            if (qi(ii,jj,kk).lt.0) then
              write (*,*) ' SEIFERT: qi < 0, STOPPED AFTER CLOUDS AND RECOPYING TO LM-VARIABLES'
              CALL model_abort (my_world_id, 999, 'negative values in qi',TRIM(yzroutine))      
            endif
            if (qs(ii,jj,kk).lt.0) then
              write (*,*) ' SEIFERT: qs < 0, STOPPED AFTER CLOUDS AND RECOPYING TO LM-VARIABLES'
              CALL model_abort (my_world_id, 999, 'negative values in qs',TRIM(yzroutine))      
            endif
            if (qg(ii,jj,kk).lt.0) then
              write (*,*) ' SEIFERT: qg < 0, STOPPED AFTER CLOUDS AND RECOPYING TO LM-VARIABLES'
              CALL model_abort (my_world_id, 999, 'negative values in qg',TRIM(yzroutine))      
            endif
            if (qh(ii,jj,kk).lt.0) then
              write (*,*) ' SEIFERT: qh < 0, STOPPED AFTER CLOUDS AND RECOPYING TO LM-VARIABLES'
              CALL model_abort (my_world_id, 999, 'negative values in qh',TRIM(yzroutine))      
            endif
            if (qnc(ii,jj,kk).lt.0) then
              write (*,*) ' SEIFERT: qnc < 0, STOPPED AFTER CLOUDS AND RECOPYING TO LM-VARIABLES'
              CALL model_abort (my_world_id, 999, 'negative values in qnc',TRIM(yzroutine))      
            endif
            if (qnr(ii,jj,kk).lt.0) then
              write (*,*) ' SEIFERT: qnr < 0, STOPPED AFTER CLOUDS AND RECOPYING TO LM-VARIABLES'
              CALL model_abort (my_world_id, 999, 'negative values in qc',TRIM(yzroutine))      
            endif
            if (qni(ii,jj,kk).lt.0) then
              write (*,*) ' SEIFERT: qni < 0, STOPPED AFTER CLOUDS AND RECOPYING TO LM-VARIABLES'
              CALL model_abort (my_world_id, 999, 'negative values in qni',TRIM(yzroutine))      
            endif
            if (qns(ii,jj,kk).lt.0) then
              write (*,*) ' SEIFERT: qns < 0, STOPPED AFTER CLOUDS AND RECOPYING TO LM-VARIABLES'
              CALL model_abort (my_world_id, 999, 'negative values in qns',TRIM(yzroutine))      
            endif
            if (qng(ii,jj,kk).lt.0) then
              write (*,*) ' SEIFERT: qng < 0, STOPPED AFTER CLOUDS AND RECOPYING TO LM-VARIABLES'
              CALL model_abort (my_world_id, 999, 'negative values in qng',TRIM(yzroutine))      
            endif
            if (qnh(ii,jj,kk).lt.0) then
              write (*,*) ' SEIFERT: qnh < 0, STOPPED AFTER CLOUDS AND RECOPYING TO LM-VARIABLES'
              CALL model_abort (my_world_id, 999, 'negative values in qnh',TRIM(yzroutine))      
            endif
!            if (rho(ii,jj,kk).lt.0) then
!              WRITE (*,*) ' SEIFERT: rho < 0, STOPPED AFTER CLOUDS AND RECOPYING TO LM-VARIABLES'
!              CALL model_abort (my_world_id, 999, 'negative values in rho',TRIM(yzroutine))      
!            endif
            if (T(ii,jj,kk,nx).lt.0) then
              WRITE (*,*) ' SEIFERT: T < 0, STOPPED AFTER CLOUDS AND RECOPYING TO LM-VARIABLES'
              CALL model_abort (my_world_id, 999, 'negative values in T',TRIM(yzroutine))      
            endif
          ENDDO
        ENDDO
      ENDDO
    END IF

    IF (debug.and.debug_maxvalglob) THEN
      ! This only works if all processors got some cloudy grid points
      xmax(1)  = MINVAL(w(its:ite,jts:jte,kts:kte,nx))
      xmax(2)  = MINVAL(qv(its:ite,jts:jte,kts:kte))
      xmax(3)  = MINVAL(qc(its:ite,jts:jte,kts:kte))
      xmax(4)  = MINVAL(qr(its:ite,jts:jte,kts:kte))
      xmax(5)  = MINVAL(qi(its:ite,jts:jte,kts:kte))
      xmax(6)  = MINVAL(qs(its:ite,jts:jte,kts:kte))
      xmax(7)  = MINVAL(qg(its:ite,jts:jte,kts:kte))  
      xmax(8)  = MINVAL(qh(its:ite,jts:jte,kts:kte))  
      xmax(9)  = MINVAL(qnc(its:ite,jts:jte,kts:kte))  
      xmax(10) = MINVAL(qni(its:ite,jts:jte,kts:kte))  
      xmax(11) = MINVAL(t(its:ite,jts:jte,kts:kte,nx))  
      IF (num_compute > 1) THEN
        CALL global_values (xmax,11,'MIN',imp_reals,icomm_cart,-1,yzerrmsg,izerror)
      END IF
      IF (isIO()) THEN
        WRITE (*,*) TRIM(yzroutine)//": output 1"
        WRITE(*,'(10x,a)') 'Minimum Values: after trafo 2'
        WRITE(*,'(A10,11A11)')   '   ', 'w','qv','qc','qr','qi','qs','qg','qh','nc','ni','tmax'
        WRITE(*,'(A10,11ES11.3)') '   ', (xmax(i), i=1,11)
      ENDIF
    END IF

    IF (ltime) CALL get_timings (i_2mom_tospecif, ntstep, dt, izerror)

    ! ... Sedimentation on the full LM grid
    IF (debug.AND.isIO()) WRITE (*,*) TRIM(yzroutine)//": rhocorr, rainrates allocated"

    ! ... ... calculate density correction and reciprocal vertical grid
    DO kk = kts, kte
      DO jj = jts, jte
        DO ii = its, ite           
          ! density correction for fall velocities
          rhocorr(ii,jj,kk) = rho0_2mom / rho(ii,jj,kk)
        ENDDO
      ENDDO
    ENDDO

    dzmin = 1.0e10
    dzmin = MIN(MINVAL(hhl(its:ite,jts:jte,kts:kte)-hhl(its:ite,jts:jte,kts+1:kte+1)),dzmin)

    ! ... ... Iteration of vertical advection (sedimentation)
    !dt_sedi = MIN(dt,4.7 * dzmin/20.0) ! MIN(dt,0.5 * dzmin/20.0)
    nt_sedi = 1 ! CEILING(dt/dt_sedi)      ! 1
    dt_sedi = zdt/nt_sedi

    IF (debug.and.isIO()) &
      & WRITE (*,'(A,i6,6(f10.2))') TRIM(yzroutine)//": calling sedi rain", nt_sedi &
      & , 20.0/dzmin*zdt/nt_sedi,dzmin

    DO n=1,nt_sedi  ! time loop with smaller timestep for rain
      ! This sedimentation routine uses the 1D scheme after Lin and Rood (1996) 
      ! and is quite robust even for large Courant Numbers (at least in 1D); 
      ! so we use a relaxed Courant number limitation. Courant number is limited at all because
      ! of numerical diffusion: in a 2-moment-scheme we need some degree of diffusion 
      ! for a good sedimentation parameterization, so the time step should be small enough 
      ! to enable lots of sedimentation steps, since the more time steps are needed for
      ! a given transport length, the more diffusion the scheme produces.
      ! (the same is used for graupel and hail)
      CALL rain_sedimentation_lm (             &
           qr     (its:ite,jts:jte,kts:kte), &
           qnr    (its:ite,jts:jte,kts:kte), &
           qc     (its:ite,jts:jte,kts:kte), &
           rain_r (its:ite,jts:jte,kts:kte),      &
           rho    (its:ite,jts:jte,kts:kte),      &
           rhocorr(its:ite,jts:jte,kts:kte),      &
           adz    (its:ite,jts:jte,kts:kte),      &
           dt_sedi,                               &
           its,ite,jts,jte,kts,kte)
    ENDDO


    ! Calculate sedimentation of ice species only when it is needed:
    IF (itype_gscp >= 1000) THEN


      !dt_sedi = MIN(dt,4.7 * dzmin/30.0)  ! MIN(dt,0.9 * dzmin/20.0)
      nt_sedi = 1 ! CEILING(dt/dt_sedi)       ! 1
      dt_sedi = zdt/nt_sedi

      IF (debug.AND.isIO()) &
           & WRITE (*,'(A,i6,4(f10.2))') TRIM(yzroutine)//": calling sedi graupel", nt_sedi &
           & , 30.0/dzmin*zdt/nt_sedi 

      DO n=1,nt_sedi  ! time loop with smaller timestep for graupel
        CALL graupel_sedimentation_lm (             &
             qg     (its:ite,jts:jte,kts:kte), &
             qng    (its:ite,jts:jte,kts:kte), &
             rain_g (its:ite,jts:jte,kts:kte),      &
             rho    (its:ite,jts:jte,kts:kte),      &
             rhocorr(its:ite,jts:jte,kts:kte),      &
             adz    (its:ite,jts:jte,kts:kte),      &
             dt_sedi,                               &
             its,ite,jts,jte,kts,kte)
      ENDDO



    ENDIF

    IF (itype_gscp >= 2000) THEN

      !dt_sedi = MIN(dt,4.7 * dzmin/30.0)  ! MIN(dt,0.9 * dzmin/20.0)
      nt_sedi = 1 ! CEILING(dt/dt_sedi)       ! 1
      dt_sedi = zdt/nt_sedi

      IF (debug.and.isIO()) &
        & WRITE (*,'(A,i6,4(f10.2))') TRIM(yzroutine)//": calling sedi hail", nt_sedi &
        & , 30.0/dzmin*zdt/nt_sedi 

      DO n=1,nt_sedi  ! time loop with smaller timestep for graupel
        CALL hail_sedimentation_lm (                &
             qh     (its:ite,jts:jte,kts:kte), &
             qnh    (its:ite,jts:jte,kts:kte), &
             rain_h (its:ite,jts:jte,kts:kte),      &
             rho    (its:ite,jts:jte,kts:kte),      &
             rhocorr(its:ite,jts:jte,kts:kte),      &
             adz    (its:ite,jts:jte,kts:kte),      &
             dt_sedi,                               &
             its,ite,jts,jte,kts,kte)
      ENDDO


    END IF


    IF (itype_gscp >= 1000) THEN


      !dt_sedi = MIN(dt,0.9 * dzmin/3.0)  ! MIN(dt,0.5 * dzmin/5.0)
      nt_sedi = 1 ! CEILING(dt/dt_sedi)      
      dt_sedi = zdt!/nt_sedi

      IF (debug.AND.isIO()) &
           & WRITE (*,'(A,i6,4(f10.2))') TRIM(yzroutine)//": calling sedi snow", nt_sedi

      DO n=1,nt_sedi  ! time loop with smaller timestep for snow
        CALL snow_sedimentation_lm (                &
             qs     (its:ite,jts:jte,kts:kte), &
             qns    (its:ite,jts:jte,kts:kte), &
             rain_s (its:ite,jts:jte,kts:kte),      &
             rho    (its:ite,jts:jte,kts:kte),      &
             rhocorr(its:ite,jts:jte,kts:kte),      &
             adz    (its:ite,jts:jte,kts:kte),      &
             dt_sedi,                               &
             its,ite,jts,jte,kts,kte)
      END DO



    ENDIF

    IF (itype_gscp >= 1000) THEN


      nt_sedi = 1 ! ... assume this works.
      dt_sedi = zdt!/nt_sedi

      IF (debug.AND.isIO()) &
           & WRITE (*,'(A,i6,4(f10.2))') TRIM(yzroutine)//": calling sedi ice", nt_sedi

      DO n=1,nt_sedi  ! time loop with small timestep
        CALL ice_sedimentation_lm (                 &
             qi     (its:ite,jts:jte,kts:kte), &
             qni    (its:ite,jts:jte,kts:kte), &
             rain_i (its:ite,jts:jte,kts:kte),      &
             rho    (its:ite,jts:jte,kts:kte),      &
             rhocorr(its:ite,jts:jte,kts:kte),      &
             adz    (its:ite,jts:jte,kts:kte),      &
             dt_sedi,                               &
             its,ite,jts,jte,kts,kte)
      END DO



    ENDIF

    IF (ltime) CALL get_timings (i_2mom_sedi, ntstep, dt, izerror)


    ! Finally, clip values of QX < 0 and then clip NCX values to their allowed range:

    qv(:,:,:)  = MAX(qv(:,:,:),  0.0_wp)
    qc(:,:,:)  = MAX(qc(:,:,:),  0.0_wp)
    qnc(:,:,:) = MAX(qnc(:,:,:), 0.0_wp)
    qr(:,:,:)  = MAX(qr(:,:,:),  0.0_wp)
    qnr(:,:,:) = MAX(qnr(:,:,:), 0.0_wp)
    IF (itype_gscp >= 1000) THEN
      qi(:,:,:)  = MAX(qi(:,:,:),  0.0_wp)
      qni(:,:,:) = MAX(qni(:,:,:), 0.0_wp)
      qs(:,:,:)  = MAX(qs(:,:,:),  0.0_wp)
      qns(:,:,:) = MAX(qns(:,:,:), 0.0_wp)
      qg(:,:,:)  = MAX(qg(:,:,:),  0.0_wp)
      qng(:,:,:) = MAX(qng(:,:,:), 0.0_wp)
    END IF
    IF (itype_gscp >= 2000) THEN
      qh(:,:,:)  = MAX(qh(:,:,:),  0.0_wp)
      qnh(:,:,:) = MAX(qnh(:,:,:), 0.0_wp)
    END IF
    
    qnc(:,:,:) = MAX(qnc(:,:,:), qc(:,:,:) / cloud%x_max)
    qnc(:,:,:) = MIN(qnc(:,:,:), qc(:,:,:) / cloud%x_min)
    qnr(:,:,:) = MAX(qnr(:,:,:), qr(:,:,:) / rain%x_max)
    qnr(:,:,:) = MIN(qnr(:,:,:), qr(:,:,:) / rain%x_min)
    IF (itype_gscp >= 1000) THEN
      qni(:,:,:) = MAX(qni(:,:,:), qi(:,:,:) / ice%x_max)
      qni(:,:,:) = MIN(qni(:,:,:), qi(:,:,:) / ice%x_min)
      qns(:,:,:) = MAX(qns(:,:,:), qs(:,:,:) / snow%x_max)
      qns(:,:,:) = MIN(qns(:,:,:), qs(:,:,:) / snow%x_min)
      qng(:,:,:) = MAX(qng(:,:,:), qg(:,:,:) / graupel%x_max)
      qng(:,:,:) = MIN(qng(:,:,:), qg(:,:,:) / graupel%x_min)
    END IF
    IF (itype_gscp >= 2000) THEN
      qnh(:,:,:) = MAX(qnh(:,:,:), qh(:,:,:) / hail%x_max)
      qnh(:,:,:) = MIN(qnh(:,:,:), qh(:,:,:) / hail%x_min)
    END IF
    
!!$ Further clipping of values < eps to 0 will be done in 
!!$ lmorg.f90, nullify_tracers(), at the end of the timestep.
!!$ However, the following cross-clipping of qnX depending
!!$ on qX is not done there and is most easily done here:

!!$ MB: for consistent output of the qnX:
!!$     if there is no mass, there should also be no
!!$     number concentration in the output and eff. radii
!!$     should be 0.0:
    
    WHERE (qc < eps) qnc = 0.0_wp
    WHERE (qr < eps) qnr = 0.0_wp
    IF (itype_gscp >= 1000) THEN
      WHERE (qi < eps) qni = 0.0_wp
      WHERE (qs < eps) qns = 0.0_wp
      WHERE (qg < eps) qng = 0.0_wp
    END IF
    IF (itype_gscp >= 2000) THEN
      WHERE (qh < eps) qnh = 0.0_wp
    END IF


    rain_r = MAX(rain_r, 0.0_wp)
    rain_i = MAX(rain_i, 0.0_wp)
    rain_s = MAX(rain_s, 0.0_wp)
    rain_g = MAX(rain_g, 0.0_wp)
    rain_h = MAX(rain_h, 0.0_wp)


    ! ... ... surface rainrate in kg/(m2*s) as calculated in sedimentation -> mm/timestep
!CDIR COLLAPSE
    DO j = 1, je
      DO i = 1, ie
        ! ... surface rainrate in kg/(m2*s)
        prr_gsp(i,j) = rain_r(i,j,ke)
        prs_gsp(i,j) = rain_s(i,j,ke) + rain_i(i,j,ke)
        prg_gsp(i,j) = rain_g(i,j,ke)
        prh_gsp(i,j) = rain_h(i,j,ke)
      ENDDO
    ENDDO
    DO k = kts, kte
!CDIR COLLAPSE
      DO j = 1, je
        DO i = 1, ie           
          qrs(i,j,k) = qr(i,j,k)+qi(i,j,k)+qs(i,j,k)+qg(i,j,k)+qh(i,j,k)
        ENDDO
      ENDDO
    ENDDO
    ! for the latent heat nudging
#ifdef NUDGING
    IF ((llhn .OR. llhnverif) .AND. lhn_qrs) THEN
      DO k = kts, kte
!CDIR COLLAPSE
        DO j = 1, je
          DO i = 1, ie           
            !qrsflux(i,j,k) = rain_i(i,j,ke)+rain_r(i,j,ke)+rain_s(i,j,ke)+rain_g(i,j,ke)+rain_h(i,j,ke)
            qrsflux(i,j,k) = rain_i(i,j,k)+rain_r(i,j,k)+rain_s(i,j,k)+rain_g(i,j,k)+rain_h(i,j,k)
          ENDDO
        ENDDO
      ENDDO
    ENDIF
#endif

    !... The qx have been changed by sedimentation, so rediagnose density again:
    CALL calrho ( t(:,:,:,nx), pp(:,:,:,nx), &
         qv(:,:,:), qc(:,:,:),&
         qrs(:,:,:), p0(:,:,:), &
         rho(:,:,:), ie, je, ke, r_d, rvd_m_o)
    
#ifdef MPI_CLOUD
    IF (MPI_DEBUG) WRITE (*,*) TRIM(yzroutine)//": dealloc_mpi"
    IF (lb_imbalance) THEN
      DEALLOCATE(lb_proc_num_bgp,lb_pp_num_bgp,lb_cgp_to_proc, &
           &     sendcnt,recvcnt,sdispl,rdispl)

    ENDIF
#endif
    DEALLOCATE(lb_proc_num_cgp)


    IF (debug.and.isIO()) WRITE (*,*) TRIM(yzroutine)//": end"

  END IF ! ... This ends the loooong if-block 'cloudy points present'

  !... deallocation
!!$    IF (debug.and.isIO()) WRITE (*,*) TRIM(yzroutine)//": calling deallocation_driver"
!!$    CALL dealloc_driver()
!!$    IF (debug.and.isIO()) WRITE (*,*) TRIM(yzroutine)//": calling deallocation_wolken"
!!$    CALL dealloc_wolken()
!!$    IF (debug.and.isIO()) WRITE (*,*) TRIM(yzroutine)//": deallocation"


  IF (debug_maxrainloc) THEN
    xmax(1) = MAXVAL(prr_gsp+prs_gsp+prg_gsp+prh_gsp)
    xmax(2) = MAXVAL(rain_r(:,:,ke))
    xmax(3) = MAXVAL(rain_i(:,:,ke))
    xmax(4) = MAXVAL(rain_s(:,:,ke))
    xmax(5) = MAXVAL(rain_g(:,:,ke))
    xmax(6) = MAXVAL(rain_h(:,:,ke))
    IF (num_compute > 1) THEN
      CALL global_values (xmax,6,'MAX',imp_reals,icomm_cart,-1,yzerrmsg,izerror)
    END IF
    IF (isIO()) THEN
      WRITE (*,*) TRIM(yzroutine)//": rainrate = ", xmax(1)*3600
      WRITE (*,*) "            rain_r   = ", xmax(2)*3600
      WRITE (*,*) "            rain_g   = ", xmax(3)*3600
      WRITE (*,*) "            rain_h   = ", xmax(4)*3600
      WRITE (*,*) "            rain_s   = ", xmax(5)*3600
      WRITE (*,*) "            rain_i   = ", xmax(6)*3600
    ENDIF
  ENDIF
    
  IF (debug_maxvalloc) THEN
    xmax(1)  = MAXVAL(w(its:ite,jts:jte,kts:kte,nx))
    xmax(2)  = MAXVAL(qv(its:ite,jts:jte,kts:kte))
    xmax(3)  = MAXVAL(qc(its:ite,jts:jte,kts:kte))
    xmax(4)  = MAXVAL(qr(its:ite,jts:jte,kts:kte))
    xmax(5)  = MAXVAL(qi(its:ite,jts:jte,kts:kte))
    xmax(6)  = MAXVAL(qs(its:ite,jts:jte,kts:kte))
    xmax(7)  = MAXVAL(qg(its:ite,jts:jte,kts:kte))  
    xmax(8)  = MAXVAL(qh(its:ite,jts:jte,kts:kte))  
    xmax(9)  = MAXVAL(qnc(its:ite,jts:jte,kts:kte))  
    xmax(10) = MAXVAL(qni(its:ite,jts:jte,kts:kte))  
    xmax(11) = MAXVAL(t(its:ite,jts:jte,kts:kte,nx))  
    IF (num_compute > 1) THEN
      CALL global_values (xmax,11,'MAX',imp_reals,icomm_cart,-1,yzerrmsg,izerror)
    END IF
    IF (isIO()) THEN
      WRITE (*,*) TRIM(yzroutine)//": output after microphysics"
      WRITE(*,'(10x,a)') 'Maximum Values:'
      WRITE(*,'(A10,11A11)')   '   ', 'w','qv','qc','qr','qi','qs','qg','qh','nc','ni','tmax'
      WRITE(*,'(A10,11ES11.3)') '   ', (xmax(i), i=1,11)
    ENDIF
    xmax(1)  = MINVAL(w(its:ite,jts:jte,kts:kte,nx))
    xmax(2)  = MINVAL(qv(its:ite,jts:jte,kts:kte))
    xmax(3)  = MINVAL(qc(its:ite,jts:jte,kts:kte))
    xmax(4)  = MINVAL(qr(its:ite,jts:jte,kts:kte))
    xmax(5)  = MINVAL(qi(its:ite,jts:jte,kts:kte))
    xmax(6)  = MINVAL(qs(its:ite,jts:jte,kts:kte))
    xmax(7)  = MINVAL(qg(its:ite,jts:jte,kts:kte))  
    xmax(8)  = MINVAL(qh(its:ite,jts:jte,kts:kte))  
    xmax(9)  = MINVAL(qnc(its:ite,jts:jte,kts:kte))  
    xmax(10) = MINVAL(qni(its:ite,jts:jte,kts:kte))  
    xmax(11) = MINVAL(t(its:ite,jts:jte,kts:kte,nx))  
    IF (num_compute > 1) THEN
      CALL global_values (xmax,11,'MIN',imp_reals,icomm_cart,-1,yzerrmsg,izerror)
    END IF
    IF (isIO()) THEN
      WRITE(*,'(10x,a)') 'Minimum Values:'
      WRITE(*,'(A10,11A11)')   '   ', 'w','qv','qc','qr','qi','qs','qg','qh','nc','ni','tmax'
      WRITE(*,'(A10,11ES11.3)') '   ', (xmax(i), i=1,11)
    END IF
  ENDIF

  ! add part of latent heating calculated in subroutine hydci_pp_gr to model latent
  ! heating field: add temperature to model latent heating field
#ifdef NUDGING
  IF (llhn .OR. llhnverif) &
    CALL get_gs_lheating ('inc',1,ke)
#endif
  ! compute temperature increment due to latent heat (for Jochen)
  IF ( ldiabf_lh ) THEN
    tinc_lh(:,:,:) = tinc_lh(:,:,:) + t(:,:,:,nx)
  END IF

  IF (ltime) CALL get_timings (i_2mom_cleanup, ntstep, dt, izerror)

  IF (isIO()) WRITE (*,*) TRIM(yzroutine)//": end"

  RETURN

CONTAINS

  SUBROUTINE saturation_adjust_h2o_qv(i, j, k, dq37, dq38, dq39, dq40, tempaenderung)
    
    
    USE wolken_konstanten, ONLY: T_f, T_3, L_wd, A_w, B_w, e_3
    USE wolken_driver, ONLY: cp
    
    IMPLICIT NONE
    
    INTEGER(kind=iintegers), INTENT(IN)          :: i, j, k
    INTEGER(kind=iintegers)                      :: ice_typ
    REAL(kind=wp) :: a1, b1, dq_d, q_sat, T_a, p_a, hlp, rlrd, rlrdm1, e_sat
    REAL(kind=wp), INTENT(OUT)   :: dq37, dq38, dq39, dq40, tempaenderung
    REAL(kind=wp) :: t_, e_ws_

! statement function for saturation waper pressure w.r.t. water:
    e_ws_(t_)  = e_3 * EXP (A_w * (t_ - T_3) / (t_ - B_w))
! end of statement function

    ice_typ   = cloud_type/1000
    
    rlrd  = r_d / r_v
    rlrdm1 = rlrd - 1.0_wp
    
    !.. Initialisieren der Ouput-Umwandlungsraten:
    dq37 = 0.0
    dq38 = 0.0
    dq39 = 0.0
    dq40 = 0.0
    tempaenderung = 0.0

    T_a   = t(i,j,k,nx)
    
    IF (T_a > T_f .OR. ice_typ == 0) THEN

      p_a   = p0(i,j,k) + pp(i,j,k,nx)
      e_sat = e_ws_ (T_a) * 1.00
      
      !...Berechnung der spezifischen Saettigungsfeuchte, Dotzek (4.33)
      q_sat = rlrd / (p_a / e_sat + rlrdm1)
      
      !...Berechnung der Adjustierungs-Terme, Dotzek, S. 35
      a1 = A_w * (T_3 - B_w) / (T_a - B_w)**2
      b1 = a1 * lewfaktor * L_wd / cp * q_sat
      
      !...Berechnung des Phasenuebergangs, Dotzek, S. 36
      dq_d = (qv(i,j,k) - q_sat) / (1.0d0 + b1)
      
      !...Verdunstung bis zur Saettigung oder bis kein Kondensat mehr da ist
      dq_d = MIN( MAX (dq_d, -qc(i,j,k)), qv(i,j,k))
      
            
      !...Berechnung der neuen H_2O-Komponenten
      qv(i,j,k)       = qv(i,j,k)       - dq_d
      qc(i,j,k)       = qc(i,j,k)       + dq_d
      

      ! Vervollstaendigung der Heizratenspeicherung:
      tempaenderung =  lewfaktor * L_wd / cp * dq_d

      !...Adjustierung der Temperatur:
      t(i,j,k,nx) = t(i,j,k,nx) + tempaenderung

      !...Begrenzung von nc so, dass der mittlere Radius nicht zu gross wird:
!      qnc(i,j,k) = MIN(qnc(i,j,k), qc(i,j,k)/cloud%x_min)
      qnc(i,j,k) = MAX(qnc(i,j,k), qc(i,j,k)/cloud%x_max)

    END IF
    
  END SUBROUTINE saturation_adjust_h2o_qv


  SUBROUTINE saturation_adjust_h2o_q(T_a, p_a, q_cloud_a, n_cloud_a, q_a, &
    & rho_a, lewfakt, dq37, dq38, dq39, dq40, tempaenderung)
    
    USE wolken_konstanten, ONLY: T_f, T_3, L_wd, A_w, B_w, e_3
    USE wolken_driver, ONLY: cp
    
    IMPLICIT NONE

    REAL(kind=wp), INTENT(in)    :: rho_a, lewfakt, p_a
    REAL(kind=wp), INTENT(inout) :: T_a, q_cloud_a, n_cloud_a, q_a
    INTEGER(kind=iintegers)          :: ice_typ
    REAL(kind=wp)                :: a1, b1, dq_d, q_sat, hlp, rlrd, rlrdm1, e_sat
    REAL(kind=wp), INTENT(OUT)   :: dq37, dq38, dq39, dq40, tempaenderung
    REAL(kind=wp) :: t_, e_ws_

! statement function for saturation waper pressure w.r.t. water:
    e_ws_(t_)  = e_3 * EXP (A_w * (t_ - T_3) / (t_ - B_w))
! end of statement function

    ice_typ   = cloud_type/1000
    
    rlrd  = r_d / r_v
    rlrdm1 = rlrd - 1.0_wp

    !.. Initialisieren der Ouput-Umwandlungsraten:
    dq37 = 0.0_wp
    dq38 = 0.0_wp
    dq39 = 0.0_wp
    dq40 = 0.0_wp
    tempaenderung = 0.0_wp
    
    IF (T_a > T_f .OR. ice_typ == 0) THEN

      e_sat = e_ws_ (T_a) * 1.00_wp
      
      !...Berechnung der spezifischen Saettigungsfeuchte, Dotzek (4.33)
      q_sat = rlrd / (p_a / e_sat + rlrdm1)
      
      !...Berechnung der Adjustierungs-Terme, Dotzek, S. 35
      a1 = A_w * (T_3 - B_w) / (T_a - B_w)**2
      b1 = a1 * lewfakt * L_wd / cp * q_sat
      
      !...Berechnung des Phasenuebergangs, Dotzek, S. 36
      dq_d = (q_a - q_sat*rho_a) / (1.0d0 + b1)
      
      !...Verdunstung bis zur Saettigung oder bis kein Kondensat mehr da ist
      dq_d = MIN( MAX (dq_d, -q_cloud_a), q_a)
      
      !...Berechnung der Uebersaettigung vor der Adjustierung:
            
      !...Berechnung der neuen H_2O-Komponenten
      q_a       = q_a       - dq_d
      q_cloud_a = q_cloud_a + dq_d
      

      ! Vervollstaendigung der Heizratenspeicherung:
      tempaenderung =  lewfaktor * L_wd / cp * dq_d / rho_a

      !...Adjustierung der Temperatur:
      T_a = T_a + tempaenderung

      !...Begrenzung von nc so, dass der mittlere Radius nicht zu klein wird:
      n_cloud_a = MAX(n_cloud_a, q_cloud_a/cloud%x_max)

    END IF
    
  END SUBROUTINE saturation_adjust_h2o_q
  
!------------------------------------------------------------------------------
! End of module procedure seifert_pp
!------------------------------------------------------------------------------

END SUBROUTINE seifert_pp

!==============================================================================

!------------------------------------------------------------------------------
! Subroutines for ice nucleation 
!------------------------------------------------------------------------------

SUBROUTINE het_nuc_phillips_table(temp,ssi,frac,imax)

IMPLICIT NONE

INTEGER(KIND=iintegers),              INTENT(in)  :: imax
REAL(KIND=wp), DIMENSION(3,imax), INTENT(out) :: frac
REAL(KIND=wp), DIMENSION(imax),   INTENT(in)  :: temp,ssi

REAL(KIND=wp) :: xt,xs,tc2, s_i, s_w, e_v, e_vw, e_vi, q_v

INTEGER :: i, tt, ss

DO i=1,imax

  IF (temp(i).GT.273.0 .OR. ssi(i).LE. 1.0) THEN

    frac(i,1) = 0.0
    frac(i,2) = 0.0
    frac(i,3) = 0.0

  ELSE

    ! calculate indices used for look-up tables
    xt = (274. - temp(i))  / ttstep
    xs = (100*(ssi(i)-1.)) / ssstep
    
    xt = MIN(xt,ttmax-DBLE(1.))
    xs = MIN(xs,ssmax-DBLE(1.))          
    tt = INT(xt)
    ss = INT(xs)
    
    ! bi-linear interpolation in look-up tables
    frac(i,1) = (tt+1-xt) * (ss+1-xs) * afrac_dust(tt,ss  ) +  (xt-tt)*(ss+1-xs) * afrac_dust(tt+1,ss  ) &
              + (tt+1-xt) * (xs-ss)   * afrac_dust(tt,ss+1) +  (xt-tt)*(xs-ss)   * afrac_dust(tt+1,ss+1)
    frac(i,2) = (tt+1-xt) * (ss+1-xs) * afrac_soot(tt,ss  ) +  (xt-tt)*(ss+1-xs) * afrac_soot(tt+1,ss  ) &
              + (tt+1-xt) * (xs-ss)   * afrac_soot(tt,ss+1) +  (xt-tt)*(xs-ss)   * afrac_soot(tt+1,ss+1)
    frac(i,3) = (tt+1-xt) * (ss+1-xs) * afrac_orga(tt,ss  ) +  (xt-tt)*(ss+1-xs) * afrac_orga(tt+1,ss  ) &
              + (tt+1-xt) * (xs-ss)   * afrac_orga(tt,ss+1) +  (xt-tt)*(xs-ss)   * afrac_orga(tt+1,ss+1)

  END IF
ENDDO

END SUBROUTINE het_nuc_phillips_table

!==============================================================================

SUBROUTINE het_nuc_phillips_table_init

IMPLICIT NONE

REAL(KIND=wp), PARAMETER :: pref=500e2_wp

REAL(KIND=wp) :: xt,xs,tc2, s_i, s_w, e_v, e_vw, e_vi, q_v

! UB_20100421>>
! REAL(KIND=8), DIMENSION(3) :: fractions
DOUBLE PRECISION, DIMENSION(3) :: fractions
! UB_20100421<<

INTEGER :: tt, ss, nunit

LOGICAL, PARAMETER :: debug = .true.

CHARACTER(LEN=30) :: FMT  = "(A,I2,A,9(E9.3,','),E9.3,A)"
CHARACTER(LEN=30) :: FMT2 = "(A,I2,30E12.3)"

WRITE (*,*) ' INITIALIZE het_nuc_phillips_table, please be patient'

afrac_dust = 0.0
afrac_soot = 0.0
afrac_orga = 0.0

! condensation/immersion freezing at water saturation

IF (debug) WRITE(*,*) " Immersion freezing"

DO tt=1,ttmax

  tc2 = (274.0 - tt*ttstep) 
  ss  = 99

! UB_20100421>>
!  CALL het_nuc_phillips(tc2,500d2,-999.d0,fractions)
  CALL het_nuc_phillips(DBLE(tc2),DBLE(pref),-999.d0,fractions)
! UB_20100421<<
  afrac_dust(tt,ss) = fractions(1)
  afrac_soot(tt,ss) = fractions(2)
  afrac_orga(tt,ss) = fractions(3)

  IF (debug) THEN
    WRITE (*,'(2x,2i5,3(2x,a,f6.2),4(2x,a,f8.4))') tt,ss,&
      ' T = ',tc2,  &
      'Si = ',S_i, &
      'Sw = ',S_w, &
      'IFRC_dust = ',100.*afrac_dust(tt,ss), &
      'IFRC_soot = ',100.*afrac_soot(tt,ss), &
      'IFRC_orga = ',100.*afrac_orga(tt,ss)
  END IF
END DO

! deposition freezing below water saturation

IF (debug) WRITE(*,*)
IF (debug) WRITE(*,*) " Deposition freezing"

DO tt=1,ttmax
  DO ss=1,ssmax

    tc2 = (274.0 - tt*ttstep) 

    ! water vapour pressure [Pa]
    e_vw = sat_pres_water(tc2)
    e_vi = sat_pres_ice(tc2)

    S_i = ss/100.0_wp * ssstep
    S_w = e_vi/e_vw * (S_i+1.0_wp)
    S_w = MIN(S_w,1.0_wp)

    e_v = e_vi * (S_i+1.)                

    IF (s_w.LT.1.0 .OR. ss.EQ.1) THEN
      q_v = spec_humi(e_v,pref)
! UB_20100421>>
!      CALL het_nuc_phillips(tc2,500d2,q_v,fractions)
      CALL het_nuc_phillips(DBLE(tc2),DBLE(pref),DBLE(q_v),fractions)
! UB_20100421<<
      afrac_dust(tt,ss) = fractions(1)
      afrac_soot(tt,ss) = fractions(2)
      afrac_orga(tt,ss) = fractions(3)
    ELSE
      afrac_dust(tt,ss) = 0.0_wp
      afrac_soot(tt,ss) = 0.0_wp
      afrac_orga(tt,ss) = 0.0_wp
    END IF
    
    IF (debug) THEN
      WRITE (*,'(2x,2i5,3(2x,a,f6.2),4(2x,a,f8.4))') tt,ss,&
        ' T = ',tc2,  &
        'Si = ',S_i, &
        'Sw = ',S_w, &
        'IFRC_dust = ',100.*afrac_dust(tt,ss), &
        'IFRC_soot = ',100.*afrac_soot(tt,ss), &
        'IFRC_orga = ',100.*afrac_orga(tt,ss)
    END IF
  END DO

END DO

! dump out DATA statements for include file
IF (my_cart_id == 0) THEN
  CALL get_free_unit(nunit)
  OPEN(nunit, file='phillips_nucleation.incf', status='new', form='formatted')
  DO ss=1,ssmax
    WRITE (nunit,FMT) 'DATA afrac_dust( 1:10,',ss,') /', afrac_dust( 1:10,ss),    '/'
    WRITE (nunit,FMT) 'DATA afrac_dust(11:20,',ss,') /', afrac_dust(11:20,ss),    '/'
    WRITE (nunit,FMT) 'DATA afrac_dust(21:30,',ss,') /', afrac_dust(21:ttmax,ss), '/'
  END DO
  DO ss=1,ssmax
    WRITE (nunit,FMT) 'DATA afrac_soot( 1:10,',ss,') /', afrac_soot( 1:10,ss),     '/'
    WRITE (nunit,FMT) 'DATA afrac_soot(11:20,',ss,') /', afrac_soot(11:20,ss),     '/'
    WRITE (nunit,FMT) 'DATA afrac_soot(21:30,',ss,') /', afrac_soot(21:ttmax,ss), '/'
  END DO
  DO ss=1,ssmax
    WRITE (nunit,FMT) 'DATA afrac_orga( 1:10,',ss,') /', afrac_orga( 1:10,ss),     '/'
    WRITE (nunit,FMT) 'DATA afrac_orga(11:20,',ss,') /', afrac_orga(11:20,ss),     '/'
    WRITE (nunit,FMT) 'DATA afrac_orga(21:30,',ss,') /', afrac_orga(21:ttmax,ss), '/'
  END DO
  DO ss=99,99
    WRITE (nunit,FMT) 'DATA afrac_dust( 1:10,',ss,') /', afrac_dust( 1:10,ss),     '/'
    WRITE (nunit,FMT) 'DATA afrac_dust(11:20,',ss,') /', afrac_dust(11:20,ss),     '/'
    WRITE (nunit,FMT) 'DATA afrac_dust(21:30,',ss,') /', afrac_dust(21:ttmax,ss), '/'
  END DO
  DO ss=99,99
    WRITE (nunit,FMT) 'DATA afrac_soot( 1:10,',ss,') /', afrac_soot( 1:10,ss),     '/'
    WRITE (nunit,FMT) 'DATA afrac_soot(11:20,',ss,') /', afrac_soot(11:20,ss),     '/'
    WRITE (nunit,FMT) 'DATA afrac_soot(21:30,',ss,') /', afrac_soot(21:ttmax,ss), '/'
  END DO
  DO ss=99,99
    WRITE (nunit,FMT) 'DATA afrac_orga( 1:10,',ss,') /', afrac_orga( 1:10,ss),     '/'
    WRITE (nunit,FMT) 'DATA afrac_orga(11:20,',ss,') /', afrac_orga(11:20,ss),     '/'
    WRITE (nunit,FMT) 'DATA afrac_orga(21:30,',ss,') /', afrac_orga(21:ttmax,ss), '/'
  END DO
  CLOSE(nunit)
  CALL release_unit(nunit)
END IF

! dump out code for NCL
IF (my_cart_id == 0) THEN
  CALL get_free_unit(nunit)
  OPEN(nunit, file='phillips_nucleation.ncl', status='new', form='formatted')
  DO ss=1,ssmax
    WRITE (nunit,FMT2) '  afrac_dust ',ss,afrac_dust(1:ttmax,ss)
  END DO
  DO ss=1,ssmax
    WRITE (nunit,FMT2) '  afrac_soot ',ss,afrac_soot(1:ttmax,ss)
  END DO
  DO ss=1,ssmax
    WRITE (nunit,FMT2) '  afrac_orga ',ss,afrac_orga(1:ttmax,ss)
  END DO
  DO ss=99,99
    WRITE (nunit,FMT2) '  afrac_dust ',ss,afrac_dust(1:ttmax,ss)
  END DO
  DO ss=99,99
    WRITE (nunit,FMT2) '  afrac_soot ',ss,afrac_soot(1:ttmax,ss)
  END DO
  DO ss=99,99
    WRITE (nunit,FMT2) '  afrac_orga ',ss,afrac_orga(1:ttmax,ss)
  END DO
  CLOSE(nunit)
  CALL release_unit(nunit)
END IF

IF(debug) WRITE (*,*) ' ... INITIALIZED TABLE.'

END SUBROUTINE het_nuc_phillips_table_init

!==============================================================================
!==============================================================================

END MODULE src_twomom_sb_interface
