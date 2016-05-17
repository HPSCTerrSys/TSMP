!+ External procedure for organizing the calls to the physics packages
!------------------------------------------------------------------------------

SUBROUTINE organize_physics (yaction, ierror, yerrmsg)

!------------------------------------------------------------------------------
!
! Description:
!   This procedure is the driving routine for calling the physical
!   parametrizations. It is called first at the beginning of the program just
!   to initialize the physical packages. Later it is called during the 
!   time-stepping in the initialization and the main program.
!
! Method:
!
! Current Code Owner: DWD, Ulrich Schaettler
!  phone:  +49  69  8062 2739
!  fax:    +49  69  8062 3721
!  email:  ulrich.schaettler@dwd.de
!
! History:
! Version    Date       Name
! ---------- ---------- ----
! 1.29       1999/05/11 Ulrich Schaettler
!  Initial release
! 1.30       1999/06/24 Matthias Raschendorfer
!  Including the call of the new turbulence routines: turbdiff, turbtran,
!  canopy_source, turbdiff_hori                 
! 1.33       1999/10/14 Matthias Raschendorfer
!  Introduction of the LOGICAL lstfnct, controlling the calculation of a
!  new stability function in sub. turbdiff.
! 1.34       1999/12/10 Ulrich Schaettler
!  Changed Calls to timing routines and included Call to hydci (lmorg before)
! 1.39       2000/05/03 Ulrich Schaettler
!  Included subroutine for Namelist input and some technical changes.
! 2.2        2000/08/18 Guenther Doms
!  Introduction of the logical switch 'lconf_avg' on Namelist input to
!  enable (default) or disable a horizontal averaging of the convective
!  forcing functions. Also, two new REAL namelist input parameters
!  'c_soil' and 'e_surf' to control surface fluxes have been introduced. 
!  An error in the check of namelist input parameters was corrected.
! 2.3        2000/11/15 Guenther Doms
!  An error in the call of the TKE-subroutines 'turbtran.incf' and
!  'turbdiff.incf' was corrected.
! 2.8        2001/07/06 Ulrich Schaettler
!  Eliminated non-necessary variables from the USE-lists
!  Introduced new NAMELIST variables for multi-layer soil-model 
!         (ke_soil, isoillevels, lmulti_layer, lmelt)
! 2.11       2001/09/28 Ulrich Schaettler
!  Renamed a variable for the multi-layer soil model and added another one
! 2.15       2002/03/20 Matthias Raschendorfer
!  Introduction of the roughness length for a typical SYNOP-station (z0m_dia)
! 2.17       2002/05/08 Ulrich Schaettler
!  Adaptations to use the new parameterization schemes
! 2.18       2002/07/16 Ulrich Schaettler
!  Added error variables in call to terra2_multlay; extended range for
!  itype_gscp for prognostic rain.
! 2.19       2002/10/24 Ulrich Schaettler
!  Eliminated call to vertical_diffusion_impl in case of 2 time level scheme
! 3.2        2003/02/07 Ulrich Schaettler
!  Moved some communications from the Convection scheme to organize_physics.
! 3.5        2003/09/02 Ulrich Schaettler
!  Adaptations for the interface of exchg_boundaries
! 3.6        2003/12/11 Reinhold Schrodin
!  Adaptations for the multi-layer soil model;
!  Modifications for checking the IOSTAT-value when reading the NAMELIST
! 3.7        2004/02/18 Ulrich Schaettler
!  Namelist Input for new variable lprogprec; ltrans_prec; rat_sea
!  Changed selection of precipitation module
! 3.8        2004/03/23 Jochen Foerstner
!  For the 2 timelevel scheme and prognostic precipitation, the precipitation
!  scheme (renamed to hydci_pp) is now also called after the dynamics.
! 3.12       2004/09/15 Christoph Schraff
!  New Namelist variable 'ldiniprec'.
! 3.13       2004/12/03 Jochen Foerstner / Thorsten Reinhardt
!  New Namelist variables for 3D turbulence
! 3.15       2005/03/03 Ulrich Schaettler
!  Editorial changes
! 3.16       2005/07/22 M. Raschendorfer, J. Foerstner, R. Schrodin
!  New Namelist variables for Physics, Turbulence scheme and shallow convection
!  Multi layer soil model completely called before convection
! 3.18       2006/03/03 Ulrich Schaettler, Dmitrii Mironov
!  Added variables for CLM version (ico2_rad, czbot_w_so, ibot_w_so)
!  Added switch llake and calls to the routines of the lake model FLake
!  Added switch l3dturb_metr for using metric terms for 3D-turbulence
!  Determination of ntke in case of restarts
! 3.19       2006/04/25 Ulrich Schaettler
!  Provide the depths of soil half levels for NetCDF IO
! 3.21       2006/12/04 Ulrich Schaettler, Burkhardt Rockel
!  Bug correction in distributing Namelist input
!  Maximum value for ico2_rad set to 6
!  In case of climate runs, call init_canopy every time step
!  Removed several NL variables and put them to TUNING (in src_setup)
!  Added NL variable lbechtold and call to Bechtold scheme (MeteoSwiss)
! 3.22       2007/01/24 Jochen Foerstner
!  Corrections for Runge-Kutta Restart: TKE has to be put to "nnew"
! V3_23        2007/03/30 Ulrich Schaettler
!  Introduction of idbg_level;
!  Eliminated hincconv, hinctura (these values can only be set in time steps)
!  Changed determination when to call the radiation scheme.
!  Changed default settings of some namelist variables
!  Added new NL variables lradtopo, nhori (Matteo Buzzi)
! V3_26        2007/05/29 Ulrich Schaettler
!  Bug correction when running with nincrad = 1
! V4_4         2008/07/16 Ulrich Schaettler, Jan-Peter Schulz
!  Adapted interface of get_timings
!  Added NL variable lsso
!  Replaced logical switches for convection (ltiedtke...) by itype_conv
!  Initialize nincrad, if only hincrad is given in NL input
! V4_5         2008/09/10 Jan-Peter Schulz
!  Activation of SSO scheme with additional Namelist parameter nincsso
!  Add possibility to use turbulence scheme without surface friction (G. Zaengl)
! V4_9         2009/07/16 Ulrich Schaettler
!  Corrected a comment for CO2
!  Check calls to get_timings
!  Adapted a check for ltrans_prec and l2tls
!  Included boundary exchange necessary for radiation averaging
! V4_10        2009/09/11 Matthias Raschendorfer, Jan-Peter Schulz
!  Introduction of LOGICALs limpltkediff, ltkesso and removing lturhor.
!  Introduction of INTEGER itype_shar
!  Moving the call of 'organize_sso' to be executed before the call of turbulence.
!  Introduction of sea-ice model
! V4_11        2009/11/30 Ekaterina Machulskaya, Juergen Helmert
!  Read Namelist switches lmulti_snow, ke_snow  (EM)
!  Read Namelist switches itype_aerosol, itype_root, itype_heatcond, 
!       itype_hydbound, lemiss, lstomata
! V4_12        2010/05/11 Michael Baldauf
!  set flag l_dzeta_d_needed, if 3D Turbulence is used
! V4_13        2010/05/11 Michael Gertz
!  Adaptions to SVN
! V4_17        2011/02/24 Ulrich Blahak
!  - Adapted interface of exchg_boundaries;
!  - corrected kzdims(1:20) -> kzdims(1:24);
!  - eliminated my_peri_neigh;
!  - eliminated lfreeslip_surf (free-slip BC and/or no-surface-heat/moisture-flux
!    conditions may now be imposed by new switches lnosurffluxes_m/h in namelist IDEAL);
!  - introduced new itype_turb=100 to specify the fixed diffusion coefficients
!    tkvhfix, tkhhfix, tkvmfix, tkhmfix in namelist IDEAL (mainly intended for idealized runs, but may also
!    be used for real cases); 
!  - added some missing communication for tkhm, tkhh, tke, tketens.
! V4_18        2011/05/26 Michael Baldauf
!  Introduced new NL switch y_scalar_advect (which replaces lsl_adv_qx and yef_adv_qx)
!  Restrict use of lradtopo to nradcoarse=1 (US)
! V4_20        2011/08/31 Ulrich Schaettler
!  Introduced call to new subroutine organize_turbulence which organizes the
!    calls to the packages for turbulence parameterization (therefore removed
!    calls to routines from src_turbulence and src_turbdiff)
!  Re-organized 'init'-phase for the turbulence (now calls init_volume_canopy
!    and init_surface_canopy
!  Introduction of NAMELIST parameter 'ltkecon' (Raschendorfer)
! Setting 'lctke=F' in case of 'ltur=F'.
! V4_21        2011/12/06 Ulrich Schaettler
!  Additional debug output
! V4_22        2012/01/31 Thorsten Reinhardt
!  Adaptations for solar zenith angle updating in every timestep.
! V4_23        2012/05/10 Ulrich Schaettler, Oli Fuhrer, CLM
!  Removed src_2timelevel and related stuff
!  Removed switches lprogprec, ltrans_prec
!  Removed Kain-Fritsch convection scheme 
!    only possible at the moment: itype_conv = 0 / 3
!  Added computation of total physical tendencies at the end of action='compute'.
!    These computations have been moved from the dynamical cores, to clean up
!    the code.
!  New Namelist parameters iy_co2_stab, lco2_stab; extended range of ico2_rad to 10
!  Introduction of prescribed soil albedo: new namelist parameter itype_albedo
!     (CLM)
! V4_25        2012/09/28 Anne Roches, Oliver Fuhrer, Ulrich Blahak
!                         Hans-Juergen Panitz (CLM)
!  Added part for the initialization of the tracer module
!  Replaced qx-variables by using them from the tracer module
!  UB:
!  Added additional tracers for the 2-moment scheme.
!  For 2-moment microphysics: added calls to init_seifert and seifert_pp;
!    added consistency checks with relevant namelist switches.
!  Implemented internal switch "l2mom_satads" to be able to switch on the 
!    extra saturation adjustments outside the microphysics parts.
!  Added plausibility checks regarding lrad and itype_aerosols in combination
!   with periodic BCs.
!  Prohibit concurrent usage of nboundlines > 3 and nradcoarse > 2, 
!    which is not working (Uli Schaettler)
!  Implement IFS convection for CLM with conditional compilation (ifdef CLM)
!  HJP for CLM:
!  Allow a new value for itype_conv (=2) to choose IFS convection
!  Introduced a warning in case itype_conv=2 and lconf_avg=TRUE (will be set
!    to FALSE then)  
!  Introduced nexch_tag for MPI boundary exchange tag to replace ntstep (HJP)
! V4_26        2012/12/06 Hans-Juergen Panitz, Ulrich Schaettler, Anne Roches
!  Computation of physical forcings (tendencies):
!   Delete lines related to itpye_conv=1 since itype_conv=1 is not valid anymore
!  Read ntke in case of restarts from restart file instead of recomputing
!  Changes and technical adaptations to the tracer handling (by Anne Roches)
! V4_27        2013/03/19 Astrid Kerkweg
!  MESSy interface introduced
!  bugfix for tke vertical dimension in exchg_boundaries call in case of l3dturb
!  Modified default values of some tuning constants to reflect settings of COSMO-EU (US)
! V4_28        2013/07/12 Ulrich Schaettler, KIT
!  Compute depth of soil half layers (czhls) in all subdomains (US)
!  Changes to adapt COSMO-ART to new tracer module
! V4_29        2013-10-02 Astrid Kerkweg, Ulrich Schaettler
!  Unification of MESSy interfaces and COSMO Tracer structure
!  Allocation and boundary exchange for tkvm, tkvh, tkhm, tkhh for 1:ke1 
!   in the vertical
! V4_30        2013/11/08 Lucio Torrisi
!  Set dt2, which is needed for a diagnostic initialization calling hydci
! V5_1         2014-11-28 Hans-Juergen Panitz, Ulrich Blahak, Matthias Raschendorfer
!                         Oliver Fuhrer, Anne Roches, Xavier Lapillonne
!  Corrected check for ltkecon (HJP)
!  Changed the format of some YUSPECIF entries for the CLM namelist tool. (UB)
!  Implemented F2003 IOMSG-mechanism for better namelist error messages.
!  Enabled lprog_tke=.true. also for the turbulence scheme with itype_turb=3
!  Catched setting of lprog_tke in case of ltur=.false. or lphys=.false.
!  Introduction of LOGICAL lscm for Single Column runs.
!  Do not interpolate and store ut_conv, vt_conv, ut_sso, vt_sso to u- or v-grid points, 
!  because they are defined on mass grid points (and at least the sso-fields are used 
!  as such in the turbulence scheme).
!  Replaced ireals by wp (working precision) (OF)
!  Removed hacks for the tracers (AR)
!  Implemented block data structure and copying from and to blocked fields (XL)
!  Included COSMO-ICON Version of microphysics  (XL)
!  For 2-moment scheme: read NL parameter iradpar_cloud (Tobias Schad)
!  For 2-moment scheme: renamed number concentrations to new official grib-shortnames
!    QNCLOUD --> NCCLOUD, etc.
!  Exchange of t_g instead of t_s after physics, because dynamics now uses
!    t_g instead of t_s for flux calculations (UB).
!  
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!==============================================================================

USE data_parameters,    ONLY:  wp, iintegers

USE data_modelconfig,   ONLY:                                                &
       ie, je, ke, ke_soil, jstartpar, jendpar, dt, dt2, czmls, msoilgrib,   &
       czhls, istartu, iendu, jstartu, jendu, istartv, iendv, jstartv, jendv,&
       ke_snow, ke1, kcm, istartpar, iendpar, istart, iend, jstart, jend,    &
       idt_qv, idt_qc,  idt_qi,  idt_qr,  idt_qs,  idt_qg,  idt_qh,          &
               idt_qnc, idt_qnr, idt_qni, idt_qns, idt_qng, idt_qnh

USE data_runcontrol,    ONLY:                                                &
       lphys, lgsp, lrad, ltur, lconv, itype_conv, lsoil, lmelt, lmelt_var,  &
       lmulti_layer, ldiniprec, lprog_qi, nradcoarse, lsuper_coolw,          &
       lradf_avg, nincrad, nextrad, hincrad, hnextrad, ninctura, nincconv,   &
       itype_tran, itype_turb, itype_gscp, itype_synd, itype_wcld,           &
       imode_tran, imode_turb, icldm_rad, icldm_tran, icldm_turb, nstart,    &
       ntstep, l2tls, ltime, nold, nnow, nnew, nuspecif,                     &
       lexpcor, lnonloc, lcpfluc, lconf_avg, lcape, lctke, ltmpcor, lprfcor, &
       itype_trvg, itype_evsl, nlgw, nbl_exchg, l3dturb,                     &
       lprog_tke, lconv_inst, lforest, luse_rttov, ntke, lseaice, llake,     &
       ico2_rad, ibot_w_so, czbot_w_so, l3dturb_metr, idbg_level,            &
       lprintdeb_all, lradtopo, nhori, hstart, lsso, nincsso,                &
       itype_sher, limpltkediff, ltkesso, ltkecon, lmulti_snow, lemiss,      &
       lstomata, itype_aerosol, itype_root, itype_heatcond, itype_hydbound,  &
       l_dzeta_d_needed, lperi_x, lperi_y, l2dim, nbl_exchg, lscm,           &
       lartif_data, y_scalar_advect, ldiabf_lh, lco2_stab, iy_co2_stab,      &
       itype_albedo, itype_lbc_qrsg, itype_outflow_qrsg, lspubc,             &
       itype_spubc, l2mom_satads, l_cosmo_art, l_pollen, l_2mom, lsppt

USE data_parallel,      ONLY:                                                &
       num_compute, my_cart_id, icomm_cart, sendbuf, isendbuflen, iexch_req, &
       imp_reals, imp_integers, imp_logical, nboundlines, my_cart_neigh,     &
       nproc, my_world_id, icomm_world, intbuf, realbuf, logbuf, ncomm_type, &
       ldatatypes, ltime_barrier, nexch_tag

USE data_fields,        ONLY:                                                &
       tkvm, tkvh, tcm, tch, qrs, t_s, t_g, t_snow, w_snow, freshsnow, qv_s, &
       w_so, w_g1, ut_conv, vt_conv, clc_con, ut_sso, vt_sso, fr_land, d_pat,&
       c_big, c_sml, r_air, plcov, sai, tai, lai, eai, hhl, h_can,           &
       tkhm, tkhh, tke, tketens, ttens, utens, vtens,                        &
       sohr, thhr, tt_conv, qvt_conv, qct_conv, qit_conv, ut_conv, vt_conv,  &
       ut_sso, vt_sso, tt_sso, tinc_lh

USE data_io,            ONLY:  ydate_ini, lbdclim, lana_qi, llb_qi, lana_qg, &
                               llb_qg, lana_qr_qs, llb_qr_qs

!CPS
#ifdef COUP_OAS_COS
USE data_fields,        ONLY:                                                &
                               tcw, trad_clm
#endif

USE data_tracer,        ONLY:  T_ADV_OFF , T_ADV_ON  , T_DIFF_OFF, T_DIFF_ON, &
                               T_TURB_OFF, T_TURB_1D , T_CONV_OFF, T_CONV_ON, &
                               T_INI_ZERO, T_INI_FILE, T_LBC_ZERO, T_LBC_FILE,&
                               T_LBC_ZEROGRAD        , T_LBC_CST             ,&
                               T_BBC_ZEROFLUX        , T_BBC_ZEROVAL         ,&
                               T_BBC_SURF_VAL        , T_RELAX_OFF           ,&
                               T_RELAX_FULL          , T_RELAX_INFLOW        ,&
                               T_DAMP_OFF            , T_DAMP_ON             ,&
                               T_CLP_OFF             , T_CLP_ON              ,&
                               T_ERR_NOTFOUND        , T_TURB_3D

USE environment,        ONLY:  exchg_boundaries, comm_barrier, model_abort
USE parallel_utilities, ONLY:  distribute_values
USE time_utilities,     ONLY:  get_timings, i_precipitation, i_radiation,    &
                               i_turbulence, i_convection, i_soil_model,     &
                               i_communications_phy, i_barrier_waiting_phy,  &
                               i_add_computations, i_sso

USE src_gscp,           ONLY:  kessler_pp, hydor_pp, hydci !, hydci_pp, hydci_pp_gr
USE gscp_interface,     ONLY:  gscp_organize, gscp_init, gscp_init_copy
USE src_radiation,      ONLY:  init_radiation, organize_radiation
USE turbulence_utilities, ONLY: init_volume_canopy, init_surface_canopy
USE turbulence_interface, ONLY: organize_turbulence 
USE src_conv_tiedtke,   ONLY:  organize_conv_tiedtke
!not yet:  USE src_conv_bechtold,  ONLY:  organize_conv_becht
USE src_conv_shallow,   ONLY:  organize_conv_shallow
USE src_soil,           ONLY:  terra1, terra2
USE src_soil_multlay,   ONLY:  terra_multlay
USE src_flake,          ONLY:  flake_init, flake_interface
USE src_sso,            ONLY:  organize_sso
USE src_seaice,         ONLY:  seaice
USE src_artifdata,      ONLY:  tkvhfix, tkhhfix, tkvmfix, tkhmfix, &
                               lnosurffluxes_m, lnosurffluxes_h

USE src_tracer,         ONLY: trcr_new, trcr_errorstr, trcr_meta_define,      &
                              trcr_meta_set, trcr_get

#ifdef TWOMOM_SB
USE src_twomom_sb_interface,   ONLY:  seifert_pp
USE wolken_konstanten,         ONLY:  init_seifert
USE data_runcontrol,           ONLY:  iradpar_cloud
#endif

#ifdef CLM
USE src_conv_ifs,       ONLY:  organize_conv_ifs
#endif

!==============================================================================

IMPLICIT NONE

!==============================================================================
! Parameter list:
CHARACTER (LEN= *),       INTENT(IN)            ::                      &
  yaction      ! action to be performed

INTEGER (KIND=iintegers), INTENT(OUT)           ::                      &
  ierror       ! error status

CHARACTER (LEN= *),       INTENT(OUT)           ::                      &
  yerrmsg      ! error message

! Local variables: 
INTEGER (KIND=iintegers)   ::   izerrstat, izerr,                       &
                                nuin, izerror, nx, i, j, k, kzdims(24), &
                                izdebug, ist, izl, izic_qi, izlbc_qi, izic_qg,&
                                izlbc_qg, izic_qr, izlbc_qr, izic_qs,         &
                                izlbc_qs, izrelax_qr, izrelax_qs, izrelax_qg, &
                                izdamp_qv, izdamp_qc, izdamp_qi, izdamp_qr,   &
                                izdamp_qs, izdamp_qg

#ifdef TWOMOM_SB
INTEGER (KIND=iintegers)   ::   izdamp_qh, izic_qh, izlbc_qh
#endif

LOGICAL                    ::   lzconv     ! if convection is computed

CHARACTER (LEN=  9)        ::   yinput             ! Namelist INPUT file
CHARACTER (LEN=255)        ::   yzerror
CHARACTER (LEN= 25)        ::   yzroutine

REAL (KIND=wp)             ::   zdt, zdt_orig

! Tracer pointers:
REAL (KIND=wp),     POINTER:: &
  qv_tens(:,:,:) => NULL()  , &
  qc_tens(:,:,:) => NULL()  , &
  qi_tens(:,:,:) => NULL()

!------------------------------------------------------------------------------
!- End of header
!------------------------------------------------------------------------------
 
!------------------------------------------------------------------------------
!- Begin Subroutine organize_physics
!------------------------------------------------------------------------------

yzroutine = 'organize_physics'
ierror    = 0_iintegers
izerr     = 0_iintegers
izerrstat = 0_iintegers
lzconv    = .FALSE.
kzdims(:) = 0_iintegers

! Initialize, whether debug output shall be done
IF (lprintdeb_all) THEN
  izdebug = idbg_level
ELSE
  IF (my_cart_id == 0) THEN
    izdebug = idbg_level
  ELSE
    izdebug = 0
  ENDIF
ENDIF

!------------------------------------------------------------------------------
! Section 1: Input of the Namelist
!------------------------------------------------------------------------------

IF (yaction == 'input') THEN

  ! Open NAMELIST-INPUT file
  IF (my_world_id == 0) THEN
    IF (idbg_level > 0) THEN
      PRINT *,'    INPUT OF THE NAMELISTS FOR PHYSICS'
    ENDIF
    yinput   = 'INPUT_PHY'
    nuin     =  1
    OPEN(nuin   , FILE=yinput  , FORM=  'FORMATTED', STATUS='UNKNOWN',  &
         IOSTAT=izerrstat)
    IF(izerrstat /= 0) THEN
      yerrmsg  = ' ERROR    *** Error while opening file INPUT_PHY *** '
      ierror   = 2001
      RETURN
    ENDIF
  ENDIF

  ! Read all NAMELIST-groups
  CALL input_phyctl (nuspecif, nuin, izerrstat)

  IF (my_world_id == 0) THEN
    ! Close file for input of the NAMELISTS
    CLOSE (nuin    , STATUS='KEEP')
  ENDIF

  IF (izerrstat < 0) THEN
    yerrmsg = 'ERROR *** while reading NAMELIST Group /PHYCTL/ ***'
    ierror  = 2003
  ELSEIF (izerrstat > 0) THEN
    yerrmsg = 'ERROR *** Wrong values occured in NAMELIST INPUT_PHY ***'
    ierror  = 2004
  ENDIF

!------------------------------------------------------------------------------
! Section 1.5: Definition of required tracers by all physics modules
!------------------------------------------------------------------------------

ELSEIF (yaction == 'tracer') THEN

  ! init error value
  izerr = 0_iintegers

  ! Define all tracers required by the physics
  ! ------------------------------------------

  ! define QV
  IF (lspubc) THEN
    izdamp_qv = T_DAMP_ON
  ELSE
    izdamp_qv = T_DAMP_OFF
  ENDIF 

  CALL trcr_new(                                                              &
    ierr           = izerr,                                                   &
    yshort_name    = 'QV',                                                    &
    igribparam     = 51,                                                      &
    igribtable     = 1,                                                       &
    yparent        = 'organize_physics',                                      &
    yunits         = 'kg kg-1',                                               &
    ystandard_name = 'specific_humidity',                                     &
    ylong_name     = 'specific humidity',                                     &
    itype_adv      = T_ADV_ON,                                                &
    itype_diff     = T_DIFF_ON,                                               &
    itype_turbmix  = T_TURB_1D,                                               &
    itype_passconv = T_CONV_OFF,                                              &
    itype_ini      = T_INI_FILE,                                              &
    itype_lbc      = T_LBC_FILE,                                              &
    itype_bbc      = T_BBC_SURF_VAL,                                          &
    itype_relax    = T_RELAX_FULL,                                            & 
    itype_damp     = izdamp_qv,                                               &
    itype_clip     = T_CLP_ON,                                                &
    idx_trcr       = idt_qv)

  IF (izerr /= 0_iintegers) THEN
    ierror  = izerr
    yerrmsg = trcr_errorstr(izerr)
    RETURN
  ENDIF

  ! define QC
  IF (lspubc) THEN
    izdamp_qc = T_DAMP_ON
  ELSE
    izdamp_qc = T_DAMP_OFF
  ENDIF 

  CALL trcr_new(                                                              &
    ierr           = izerr,                                                   &
    yshort_name    = 'QC',                                                    &
    igribparam     = 31,                                                      &
    igribtable     = 2,                                                       &
    yparent        = 'organize_physics',                                      &
    yunits         = 'kg kg-1',                                               &
    ystandard_name = 'mass_fraction_of_cloud_liquid_water_in_air',            &
    ylong_name     = 'specific cloud liquid water content',                   &
    itype_adv      = T_ADV_ON,                                                &
    itype_diff     = T_DIFF_ON,                                               &
    itype_turbmix  = T_TURB_1D,                                               &
    itype_passconv = T_CONV_OFF,                                              &
    itype_ini      = T_INI_FILE,                                              &
    itype_lbc      = T_LBC_FILE,                                              &
    itype_bbc      = T_BBC_ZEROFLUX,                                          &
    itype_relax    = T_RELAX_FULL,                                            &
    itype_damp     = izdamp_qc,                                               &
    itype_clip     = T_CLP_ON,                                                &
    idx_trcr       = idt_qc)

  IF (izerr /= 0_iintegers) THEN
    ierror  = izerr
    yerrmsg = trcr_errorstr(izerr)
    RETURN
  ENDIF

  ! define QI 
  IF (lprog_qi) THEN
    IF (lspubc) THEN
      izdamp_qi = T_DAMP_ON
    ELSE
      izdamp_qi = T_DAMP_OFF
    ENDIF 
    IF (lana_qi) THEN
      izic_qi = T_INI_FILE
    ELSE
      izic_qi = T_INI_ZERO
    ENDIF
    IF (llb_qi) THEN
      izlbc_qi = T_LBC_FILE
    ELSE
      izlbc_qi = T_LBC_ZERO     ! this might be changed by the user 
                                ! to T_LBC_CST or T_LBC_ZEROGRAD
    ENDIF

    CALL trcr_new(                                                            &
      ierr           = izerr,                                                 &
      yshort_name    = 'QI',                                                  &
      igribparam     = 33,                                                    &
      igribtable     = 2,                                                     &
      yparent        = 'organize_physics',                                    &
      yunits         = 'kg kg-1',                                             &
      ystandard_name = 'mass_fraction_of_cloud_ice_in_air',                   &
      ylong_name     = 'specific cloud ice content',                          &
      itype_adv      = T_ADV_ON,                                              &
      itype_diff     = T_DIFF_OFF,                                            &
      itype_turbmix  = T_TURB_1D,                                             &
      itype_passconv = T_CONV_OFF,                                            &
      itype_ini      = izic_qi,                                               &
      itype_lbc      = izlbc_qi,                                              &
      itype_bbc      = T_BBC_ZEROFLUX,                                        &
      itype_relax    = T_RELAX_FULL,                                          &
      itype_damp     = izdamp_qi,                                             &
      itype_clip     = T_CLP_ON,                                              &
      idx_trcr       = idt_qi)

    IF (izerr /= 0_iintegers) THEN
      ierror  = izerr
      yerrmsg = trcr_errorstr(izerr)
      RETURN
    ENDIF
  ENDIF

  ! define QG
  IF (itype_gscp >= 4) THEN
    IF (lspubc) THEN
      IF (itype_spubc == 1_iintegers) THEN
        izdamp_qg = T_DAMP_OFF
      ELSE
        izdamp_qg = T_DAMP_ON
      ENDIF
    ELSE
      izdamp_qg = T_DAMP_OFF
    ENDIF 
    IF (lana_qg) THEN
      izic_qg = T_INI_FILE
    ELSE
      izic_qg = T_INI_ZERO
    ENDIF

    IF (llb_qg) THEN
      izlbc_qg = T_LBC_FILE
      IF (itype_outflow_qrsg == 1) THEN
        izrelax_qg = T_RELAX_FULL
      ELSEIF (itype_outflow_qrsg == 2) THEN
        izrelax_qg = T_RELAX_INFLOW
      ENDIF
    ELSE
      izrelax_qg = T_RELAX_OFF
      IF (itype_lbc_qrsg == 1) THEN 
        izlbc_qg = T_LBC_ZEROGRAD 
      ELSEIF (itype_lbc_qrsg == 2) THEN
        izlbc_qg = T_LBC_ZERO
      ENDIF
    ENDIF
 
    CALL trcr_new(                                                            &
      ierr           = izerr,                                                 &
      yshort_name    = 'QG',                                                  &
      igribparam     = 39,                                                    &
      igribtable     = 2,                                                     &
      yparent        = 'organize_physics',                                    &
      yunits         = 'kg kg-1',                                             &
      ystandard_name = 'mass_fraction_of_graupel_in_air',                     &
      ylong_name     = 'specific graupel content',                            &
      itype_adv      = T_ADV_ON,                                              &
      itype_diff     = T_DIFF_OFF,                                            &
      itype_turbmix  = T_TURB_OFF,                                            &
      itype_passconv = T_CONV_OFF,                                            &
      itype_ini      = izic_qg,                                               &
      itype_lbc      = izlbc_qg,                                              &
      itype_bbc      = T_BBC_ZEROFLUX,                                        &
      itype_relax    = izrelax_qg,                                            &
      itype_damp     = izdamp_qg,                                             &
      itype_clip     = T_CLP_ON,                                              &
      idx_trcr       = idt_qg)

    ! check for errors
    IF (izerr /= 0_iintegers) THEN
      ierror  = izerr
      yerrmsg = trcr_errorstr(izerr)
      RETURN
    ENDIF
  ENDIF

  ! define QR 
  IF (lspubc) THEN
    IF (itype_spubc == 1_iintegers) THEN
      izdamp_qr = T_DAMP_OFF
    ELSE
      izdamp_qr = T_DAMP_ON
    ENDIF
  ELSE
    izdamp_qr = T_DAMP_OFF
  ENDIF 
  IF (lana_qr_qs) THEN
    izic_qr = T_INI_FILE
  ELSE
    izic_qr = T_INI_ZERO
  ENDIF
  IF (llb_qr_qs) THEN
    izlbc_qr = T_LBC_FILE
    IF (itype_outflow_qrsg == 1) THEN
      izrelax_qr = T_RELAX_FULL
    ELSEIF (itype_outflow_qrsg == 2) THEN
      izrelax_qr = T_RELAX_INFLOW
    ENDIF
  ELSE
    izrelax_qr = T_RELAX_OFF
    IF (itype_lbc_qrsg == 1) THEN 
      izlbc_qr = T_LBC_ZEROGRAD 
    ELSEIF (itype_lbc_qrsg == 2) THEN
      izlbc_qr = T_LBC_ZERO
    ENDIF
  ENDIF

  CALL trcr_new(                                                              &
    ierr           = izerr,                                                   &
    yshort_name    = 'QR',                                                    &
    igribparam     = 35,                                                      &
    igribtable     = 2,                                                       &
    yparent        = 'organize_physics',                                      &
    yunits         = 'kg kg-1',                                               &
    ystandard_name = 'mass_fraction_of_rain_in_air',                          &
    ylong_name     = 'specific rain content',                                 &
    itype_adv      = T_ADV_ON,                                                &
    itype_diff     = T_DIFF_OFF,                                              &
    itype_turbmix  = T_TURB_OFF,                                              &
    itype_passconv = T_CONV_OFF,                                              &
    itype_ini      = izic_qr,                                                 &
    itype_lbc      = izlbc_qr,                                                &
    itype_bbc      = T_BBC_ZEROFLUX,                                          &
    itype_relax    = izrelax_qr,                                              &
    itype_damp     = izdamp_qr,                                               &
    itype_clip     = T_CLP_ON,                                                &
    idx_trcr       = idt_qr)

  IF (izerr /= 0_iintegers) THEN
    ierror  = izerr
    yerrmsg = trcr_errorstr(izerr)
    RETURN
  ENDIF

  ! define QS
  IF (itype_gscp > 1_iintegers) THEN
    IF (lspubc) THEN
      IF (itype_spubc == 1_iintegers) THEN
        izdamp_qs = T_DAMP_OFF
      ELSE
        izdamp_qs = T_DAMP_ON
      ENDIF
    ELSE
      izdamp_qs = T_DAMP_OFF
    ENDIF 
    IF (lana_qr_qs) THEN
      izic_qs = T_INI_FILE
    ELSE
      izic_qs = T_INI_ZERO
    ENDIF
    IF (llb_qr_qs) THEN
      izlbc_qs = T_LBC_FILE
      IF (itype_outflow_qrsg == 1) THEN
        izrelax_qs = T_RELAX_FULL
      ELSEIF (itype_outflow_qrsg == 2) THEN
        izrelax_qs = T_RELAX_INFLOW
      ENDIF
    ELSE
      izrelax_qs = T_RELAX_OFF
      IF (itype_lbc_qrsg == 1) THEN 
        izlbc_qs = T_LBC_ZEROGRAD 
      ELSEIF (itype_lbc_qrsg == 2) THEN
        izlbc_qs = T_LBC_ZERO
      ENDIF
    ENDIF

    CALL trcr_new(                                                            &
      ierr           = izerr,                                                 &
      yshort_name    = 'QS',                                                  &
      igribparam     = 36,                                                    &
      igribtable     = 2,                                                     &
      yparent        = 'organize_physics',                                    &
      yunits         = 'kg kg-1',                                             &
      ystandard_name = 'mass_fraction_of_snow_in_air',                        &
      ylong_name     = 'specific snow content',                               &
      itype_adv      = T_ADV_ON,                                              &
      itype_diff     = T_DIFF_OFF,                                            &
      itype_turbmix  = T_TURB_OFF,                                            &
      itype_passconv = T_CONV_OFF,                                            &
      itype_ini      = izic_qs,                                               &
      itype_lbc      = izlbc_qs,                                              &
      itype_bbc      = T_BBC_ZEROFLUX,                                        &
      itype_relax    = izrelax_qs,                                            &
      itype_damp     = izdamp_qs,                                             &
      itype_clip     = T_CLP_ON,                                              &
      idx_trcr       = idt_qs)

    IF (izerr /= 0_iintegers) THEN
      ierror  = izerr
      yerrmsg = trcr_errorstr(izerr)
      RETURN
    ENDIF
  ENDIF

#ifdef TWOMOM_SB
  IF ( itype_gscp >= 100 ) THEN

    ! define NCCLOUD
   
    ! Overtake value of
    !   izdamp_qc
    ! from the definition of QC above to do consistent things!
    
    CALL trcr_new(                                                                 &
         ierr           = izerr,                                                   &
         yshort_name    = 'NCCLOUD',                                               &
         igribparam     = 221,                                                     &
         igribtable     = 2,                                                       &
         yparent        = 'organize_physics',                                      &
         yunits         = 'kg-1',                                                  &
         ystandard_name = 'specif_number_of_cloud_droplets_in_air',                &
         ylong_name     = 'specific cloud droplet number',                         &
         itype_adv      = T_ADV_ON,                                                &
         itype_diff     = T_DIFF_ON,                                               &
         itype_turbmix  = T_TURB_3D,                                               &
         itype_passconv = T_CONV_OFF,                                              &
         itype_ini      = T_INI_FILE,                                              &
         itype_lbc      = T_LBC_FILE,                                              &
         itype_bbc      = T_BBC_ZEROFLUX,                                          &
         itype_relax    = T_RELAX_FULL,                                            &
         itype_damp     = izdamp_qc,                                               &
         itype_clip     = T_CLP_ON,                                                &
         idx_trcr       = idt_qnc)

    IF (izerr /= 0_iintegers) THEN
      ierror  = izerr
      yerrmsg = trcr_errorstr(izerr)
      RETURN
    ENDIF

    ! define NCICE 
    ! ------------

    IF (lprog_qi) THEN

      ! Overtake values of
      !   izdamp_qi
      !   izic_qi
      !   izlbc_qi
      ! from the definition of QI above to do consistent things!

      CALL trcr_new(                                                               &
           ierr           = izerr,                                                 &
           yshort_name    = 'NCICE',                                               &
           igribparam     = 223,                                                   &
           igribtable     = 2,                                                     &
           yparent        = 'organize_physics',                                    &
           yunits         = 'kg-1',                                                &
           ystandard_name = 'specif_number_of_cloud_ice_in_air',                   &
           ylong_name     = 'specific cloud ice number content',                   &
           itype_adv      = T_ADV_ON,                                              &
           itype_diff     = T_DIFF_OFF,                                            &
           itype_turbmix  = T_TURB_3D,                                             &
           itype_passconv = T_CONV_OFF,                                            &
           itype_ini      = izic_qi,                                               &
           itype_lbc      = izlbc_qi,                                              &
           itype_bbc      = T_BBC_ZEROFLUX,                                        &
           itype_relax    = T_RELAX_FULL,                                          &
           itype_damp     = izdamp_qi,                                             &
           itype_clip     = T_CLP_ON,                                              &
           idx_trcr       = idt_qni)
 
      IF (izerr /= 0_iintegers) THEN
        ierror  = izerr
        yerrmsg = trcr_errorstr(izerr)
        RETURN
      ENDIF

    ENDIF

    ! define NCRAIN
    ! -------------

    ! Overtake values of
    !    izdamp_qr
    !    izic_qr
    !    izlbc_qr
    !    izrelax_qr
    ! from the definition of QR above to do consistent things!

    CALL trcr_new(                                                                 &
         ierr           = izerr,                                                   &
         yshort_name    = 'NCRAIN',                                                &
         igribparam     = 222,                                                     &
         igribtable     = 2,                                                       &
         yparent        = 'organize_physics',                                      &
         yunits         = 'kg-1',                                                  &
         ystandard_name = 'specif_number_of_rain_drops_in_air',                    &
         ylong_name     = 'specific rain drop number',                             &
         itype_adv      = T_ADV_ON,                                                &
         itype_diff     = T_DIFF_OFF,                                              &
         itype_turbmix  = T_TURB_OFF,                                              &
         itype_passconv = T_CONV_OFF,                                              &
         itype_ini      = izic_qr,                                                 &
         itype_lbc      = izlbc_qr,                                                &
         itype_bbc      = T_BBC_ZEROFLUX,                                          &
         itype_relax    = izrelax_qr,                                              &
         itype_damp     = izdamp_qr,                                               &
         itype_clip     = T_CLP_ON,                                                &
         idx_trcr       = idt_qnr)

    IF (izerr /= 0_iintegers) THEN
      ierror  = izerr
      yerrmsg = trcr_errorstr(izerr)
      RETURN
    ENDIF

    ! define NCSNOW
    ! -------------

    ! Overtake values of
    !    izdamp_qg
    !    izic_qs
    !    izlbc_qs
    !    izrelax_qs
    ! from the definition of QS above to do consistent things!

    CALL trcr_new(                                                               &
         ierr           = izerr,                                                 &
         yshort_name    = 'NCSNOW',                                              &
         igribparam     = 224,                                                   &
         igribtable     = 2,                                                     &
         yparent        = 'organize_physics',                                    &
         yunits         = 'kg-1',                                                &
         ystandard_name = 'specif_number_of_snow_flakes_in_air',                 &
         ylong_name     = 'specific snow flake number',                          &
         itype_adv      = T_ADV_ON,                                              &
         itype_diff     = T_DIFF_OFF,                                            &
         itype_turbmix  = T_TURB_OFF,                                            &
         itype_passconv = T_CONV_OFF,                                            &
         itype_ini      = izic_qs,                                               &
         itype_lbc      = izlbc_qs,                                              &
         itype_bbc      = T_BBC_ZEROFLUX,                                        &
         itype_relax    = izrelax_qs,                                            &
         itype_damp     = izdamp_qs,                                             &
         itype_clip     = T_CLP_ON,                                              &
         idx_trcr       = idt_qns)

    IF (izerr /= 0_iintegers) THEN
      ierror  = izerr
      yerrmsg = trcr_errorstr(izerr)
      RETURN
    ENDIF

    ! define NCGRAUPEL
    ! ------------

    ! Overtake values of
    !    izdamp_qg
    !    izic_qg
    !    izlbc_qg
    !    izrelax_qg
    ! from the definition of QG above to do consistent things!

    CALL trcr_new(                                                               &
         ierr           = izerr,                                                 &
         yshort_name    = 'NCGRAUPEL',                                           &
         igribparam     = 225,                                                   &
         igribtable     = 2,                                                     &
         yparent        = 'organize_physics',                                    &
         yunits         = 'kg-1',                                                &
         ystandard_name = 'specif_number_of_graupel_in_air',                     &
         ylong_name     = 'specific graupel number',                             &
         itype_adv      = T_ADV_ON,                                              &
         itype_diff     = T_DIFF_OFF,                                            &
         itype_turbmix  = T_TURB_OFF,                                            &
         itype_passconv = T_CONV_OFF,                                            &
         itype_ini      = izic_qg,                                               &
         itype_lbc      = izlbc_qg,                                              &
         itype_bbc      = T_BBC_ZEROFLUX,                                        &
         itype_relax    = izrelax_qg,                                            &
         itype_damp     = izdamp_qg,                                             &
         itype_clip     = T_CLP_ON,                                              &
         idx_trcr       = idt_qng)

    ! check for errors
    IF (izerr /= 0_iintegers) THEN
      ierror  = izerr
      yerrmsg = trcr_errorstr(izerr)
      RETURN
    ENDIF

    ! define QH (in any case, not only for itype_gscp >= 2000)
    ! --------------------------------------------------------

    CALL trcr_new(                                                               &
         ierr           = izerr,                                                 &
         yshort_name    = 'QH',                                                  &
         igribparam     = 98,                                                    &
         igribtable     = 2,                                                     &
         yparent        = 'organize_physics',                                    &
         yunits         = 'kg kg-1',                                             &
         ystandard_name = 'mass_fraction_of_hail_in_air',                        &
         ylong_name     = 'specific hail content',                               &
         itype_adv      = T_ADV_ON,                                              &
         itype_diff     = T_DIFF_OFF,                                            &
         itype_turbmix  = T_TURB_OFF,                                            &
         itype_passconv = T_CONV_OFF,                                            &
         itype_ini      = T_INI_ZERO,                                            &
         itype_lbc      = T_LBC_ZERO,                                            &
         itype_bbc      = T_BBC_ZEROFLUX,                                        &
         itype_relax    = T_RELAX_OFF,                                           &
         itype_damp     = T_DAMP_OFF,                                            &
         itype_clip     = T_CLP_ON,                                              &
         idx_trcr       = idt_qh)

    ! check for errors
    IF (izerr /= 0_iintegers) THEN
      ierror  = izerr
      yerrmsg = trcr_errorstr(izerr)
      RETURN
    ENDIF

    ! define NCHAIL
    ! -------------

    ! Overtake values of
    !    izdamp_qg
    !    izic_qg
    !    izlbc_qg
    !    izrelax_qg
    ! from the definition of QG above to do consistent things!

    CALL trcr_new(                                                               &
         ierr           = izerr,                                                 &
         yshort_name    = 'NCHAIL',                                              &
         igribparam     = 226,                                                   &
         igribtable     = 2,                                                     &
         yparent        = 'organize_physics',                                    &
         yunits         = 'kg-1',                                                &
         ystandard_name = 'specif_number_of_hail_in_air',                        &
         ylong_name     = 'specific hail number',                                &
         itype_adv      = T_ADV_ON,                                              &
         itype_diff     = T_DIFF_OFF,                                            &
         itype_turbmix  = T_TURB_OFF,                                            &
         itype_passconv = T_CONV_OFF,                                            &
         itype_ini      = T_INI_ZERO,                                            &
         itype_lbc      = T_LBC_ZERO,                                            &
         itype_bbc      = T_BBC_ZEROFLUX,                                        &
         itype_relax    = T_RELAX_OFF,                                           &
         itype_damp     = T_DAMP_OFF,                                            &
         itype_clip     = T_CLP_ON,                                              &
         idx_trcr       = idt_qnh)

    ! check for errors
    IF (izerr /= 0_iintegers) THEN
      ierror  = izerr
      yerrmsg = trcr_errorstr(izerr)
      RETURN
    ENDIF


  END IF
#endif

  ! Define and set all additional metadata required by the microphysics
  ! -------------------------------------------------------------------

  ! This list of metadata is necessary to produce good quality results with the 
  ! Leapfrog dycore. These metadata should however be used with caution for other tracers.
  ! Testing recommended!

  ! SP_ADV_LF: Define advection schemes and time stepping scheme combinations
  !            in case of LF core (src_leapfrog.f90, slow_tendencies.f90,
  !            src_relaxation.f90)
  !            1: hadv_cd2 (centered differences, LF time step)
  !               and implicit advection in the vertical
  !            2: hadv_pd (positive definite, Euler forward)
  !               and vadv_pd (explicit vertical advection, positive definite)
  !            3: interpol_sl_trilin (Semi-Lagrange, LF timestep)
  ! WARNING: SP_ADV_LF must be 1 (default) for QV !!!
  ! WARNING: SP_ADV_LF = 2 or 3 is not compatible with T_ADV_OFF !!!

  IF ( .NOT. l2tls ) THEN
    CALL trcr_meta_define(izerr, 'SP_ADV_LF', 1_iintegers)
    IF (izerr /= 0_iintegers) THEN
      ierror  = izerr
      yerrmsg = trcr_errorstr(izerr)
      RETURN
    ENDIF
    CALL trcr_meta_set(izerr, idt_qi, 'SP_ADV_LF', 2_iintegers)
    IF (izerr /= 0_iintegers .AND. izerr /= T_ERR_NOTFOUND) THEN
      ierror  = izerr
      yerrmsg = trcr_errorstr(izerr)
      RETURN
    ENDIF
    CALL trcr_meta_set(izerr, idt_qr, 'SP_ADV_LF', 3_iintegers)
    IF (izerr /= 0_iintegers .AND. izerr /= T_ERR_NOTFOUND) THEN
      ierror  = izerr
      yerrmsg = trcr_errorstr(izerr)
      RETURN
    ENDIF
    CALL trcr_meta_set(izerr, idt_qs, 'SP_ADV_LF', 3_iintegers)
    IF (izerr /= 0_iintegers .AND. izerr /= T_ERR_NOTFOUND) THEN
      ierror  = izerr
      yerrmsg = trcr_errorstr(izerr)
      RETURN
    ENDIF
    CALL trcr_meta_set(izerr, idt_qg, 'SP_ADV_LF', 3_iintegers)
    IF (izerr /= 0_iintegers .AND. izerr /= T_ERR_NOTFOUND) THEN
      ierror  = izerr
      yerrmsg = trcr_errorstr(izerr)
      RETURN
    ENDIF
#ifdef TWOMOM_SB
    CALL trcr_meta_set(izerr, idt_qh, 'SP_ADV_LF', 3_iintegers)
    IF (izerr /= 0_iintegers .AND. izerr /= T_ERR_NOTFOUND) THEN
      ierror  = izerr
      yerrmsg = trcr_errorstr(izerr)
      RETURN
    ENDIF
    CALL trcr_meta_set(izerr, idt_qnc, 'SP_ADV_LF', 2_iintegers)
    IF (izerr /= 0_iintegers .AND. izerr /= T_ERR_NOTFOUND) THEN
      ierror  = izerr
      yerrmsg = trcr_errorstr(izerr)
      RETURN
    ENDIF
    CALL trcr_meta_set(izerr, idt_qnr, 'SP_ADV_LF', 3_iintegers)
    IF (izerr /= 0_iintegers .AND. izerr /= T_ERR_NOTFOUND) THEN
      ierror  = izerr
      yerrmsg = trcr_errorstr(izerr)
      RETURN
    ENDIF
    CALL trcr_meta_set(izerr, idt_qni, 'SP_ADV_LF', 3_iintegers)
    IF (izerr /= 0_iintegers .AND. izerr /= T_ERR_NOTFOUND) THEN
      ierror  = izerr
      yerrmsg = trcr_errorstr(izerr)
      RETURN
    ENDIF
    CALL trcr_meta_set(izerr, idt_qns, 'SP_ADV_LF', 3_iintegers)
    IF (izerr /= 0_iintegers .AND. izerr /= T_ERR_NOTFOUND) THEN
      ierror  = izerr
      yerrmsg = trcr_errorstr(izerr)
      RETURN
    ENDIF
    CALL trcr_meta_set(izerr, idt_qng, 'SP_ADV_LF', 3_iintegers)
    IF (izerr /= 0_iintegers .AND. izerr /= T_ERR_NOTFOUND) THEN
      ierror  = izerr
      yerrmsg = trcr_errorstr(izerr)
      RETURN
    ENDIF
    CALL trcr_meta_set(izerr, idt_qnh, 'SP_ADV_LF', 3_iintegers)
    IF (izerr /= 0_iintegers .AND. izerr /= T_ERR_NOTFOUND) THEN
      ierror  = izerr
      yerrmsg = trcr_errorstr(izerr)
      RETURN
    ENDIF
#endif
  ENDIF

  ! BD_SET_FORCED: Force application of a boundary condition (zero gradiend or zero
  !                value) after dynamics also if boundary data from a file is available
  !                (set_trcr_special_bc() in lmorg.f90). This is to reproduce the original
  !                behaviour of QR, QS, QG.
  !                0 : do not do a forced application of a boundary condition
  !                1 : do a forced boundary condition with zero gradiend
  !                2 : do a forced boundary condition with zero value
  CALL trcr_meta_define(izerr, 'BD_SET_FORCED', 0_iintegers)
  IF (izerr /= 0_iintegers) THEN
    ierror  = izerr
    yerrmsg = trcr_errorstr(izerr)
    RETURN
  ENDIF
  IF (itype_lbc_qrsg == 1 .OR. itype_lbc_qrsg == 2) THEN
    CALL trcr_meta_set(izerr, idt_qr, 'BD_SET_FORCED', itype_lbc_qrsg)   
    IF (izerr /= 0_iintegers .AND. izerr /= T_ERR_NOTFOUND) THEN
      ierror  = izerr
      yerrmsg = trcr_errorstr(izerr)
      RETURN
    ENDIF
    CALL trcr_meta_set(izerr, idt_qs, 'BD_SET_FORCED', itype_lbc_qrsg)
    IF (izerr /= 0_iintegers .AND. izerr /= T_ERR_NOTFOUND) THEN
      ierror  = izerr
      yerrmsg = trcr_errorstr(izerr)
      RETURN
    ENDIF
    CALL trcr_meta_set(izerr, idt_qg, 'BD_SET_FORCED', itype_lbc_qrsg)
    IF (izerr /= 0_iintegers .AND. izerr /= T_ERR_NOTFOUND) THEN
      ierror  = izerr
      yerrmsg = trcr_errorstr(izerr)
      RETURN
    ENDIF
#ifdef TWOMOM_SB
    CALL trcr_meta_set(izerr, idt_qh, 'BD_SET_FORCED', itype_lbc_qrsg)
    IF (izerr /= 0_iintegers .AND. izerr /= T_ERR_NOTFOUND) THEN
      ierror  = izerr
      yerrmsg = trcr_errorstr(izerr)
      RETURN
    ENDIF
    CALL trcr_meta_set(izerr, idt_qnr, 'BD_SET_FORCED', itype_lbc_qrsg)
    IF (izerr /= 0_iintegers .AND. izerr /= T_ERR_NOTFOUND) THEN
      ierror  = izerr
      yerrmsg = trcr_errorstr(izerr)
      RETURN
    ENDIF
    CALL trcr_meta_set(izerr, idt_qns, 'BD_SET_FORCED', itype_lbc_qrsg)
    IF (izerr /= 0_iintegers .AND. izerr /= T_ERR_NOTFOUND) THEN
      ierror  = izerr
      yerrmsg = trcr_errorstr(izerr)
      RETURN
    ENDIF
    CALL trcr_meta_set(izerr, idt_qng, 'BD_SET_FORCED', itype_lbc_qrsg)
    IF (izerr /= 0_iintegers .AND. izerr /= T_ERR_NOTFOUND) THEN
      ierror  = izerr
      yerrmsg = trcr_errorstr(izerr)
      RETURN
    ENDIF
    CALL trcr_meta_set(izerr, idt_qnh, 'BD_SET_FORCED', itype_lbc_qrsg)
    IF (izerr /= 0_iintegers .AND. izerr /= T_ERR_NOTFOUND) THEN
      ierror  = izerr
      yerrmsg = trcr_errorstr(izerr)
      RETURN
    ENDIF
#endif
  ENDIF

  ! MASSFLX_CLP: Massflux correction clipping in both LF and RK core (after
  !              turbulent diffusion of tracers in slow_tendencies.f90 and
  !              src_slow_tendencies_rk.f90). This is to reproduce the original
  !              behaviour of QV and QC.
  !              NOTE: this is only active if T_CLP is set to T_CLP_ON
  !              NOTE: this cannot be switched on at the same time as CLP_10E-12

  IF ( .NOT. l2tls ) THEN
    CALL trcr_meta_define(izerr, 'MASSFLX_CLP', .FALSE.)
    IF (izerr /= 0_iintegers) THEN
      ierror  = izerr
      yerrmsg = trcr_errorstr(izerr)
      RETURN
    ENDIF
    CALL trcr_meta_set(izerr, idt_qv, 'MASSFLX_CLP', .TRUE.)
    IF (izerr /= 0_iintegers) THEN
      ierror  = izerr
      yerrmsg = trcr_errorstr(izerr)
      RETURN
    ENDIF
    CALL trcr_meta_set(izerr, idt_qc, 'MASSFLX_CLP', .TRUE.)
    IF (izerr /= 0_iintegers) THEN
      ierror  = izerr
      yerrmsg = trcr_errorstr(izerr)
      RETURN
    ENDIF
  ENDIF


#ifdef TWOMOM_SB
  ! Initialize hydrometeor parameters for the 2-Moment cloud microphysical
  ! scheme of Seifert & Beheng
  IF (itype_gscp >= 100 ) THEN
    CALL init_seifert(itype_gscp)
  END IF

  ! ONLY EFFECTIVE IF itype_gscp >= 100:
  ! Determine if all the satads should be done
  ! (like for the 1-moment schemes), not just the
  ! satad after the microphysics at the end of the timestep:
  l2mom_satads = .FALSE.

  ! However, l2mom_do_extra_satads = .TRUE. can only be used
  ! for itype_gscp = ??[6-9]? (e.g., 2767, but not 2737),
  ! i.e., for the Segal-Khain parameterization of cloud nucleation.
  ! For the other cloud nucleation types, we need supersaturation
  ! after the dynamics.
  IF (l2mom_satads .AND. itype_gscp >= 100) THEN
    IF (MODULO(itype_gscp/10,10) < 6) THEN
      yerrmsg = 'ERROR ** l2mom_satads=.TRUE. only possible for itype_gscp=??[6-9]? **'
      ierror  = 2005
    END IF
  END IF
#endif

!------------------------------------------------------------------------------
! Section 2: Initialization of the packages for the first call
!------------------------------------------------------------------------------

ELSEIF (yaction == 'init') THEN

  IF (izdebug > 0) THEN
    PRINT *, '    PHYSICAL PACKAGES'
  ENDIF

  ! Initialize values for grid scale precipitation
  IF (lgsp) CALL gscp_init

  ! Initialize values for the radiation and the canopy fields
  IF (lrad) CALL init_radiation
  IF (lco2_stab) THEN
    IF (izdebug > 1) THEN
      WRITE (*,'(A44,I4,A4)') '           stabilised GHG forcing from year ', iy_co2_stab, ' on!'
    ENDIF
  ENDIF

  IF (ltur .OR. lsoil) THEN
     ! Determine the index of the highest layer which still is a canopy layer

     ! h_can is a primary external parameter. The initial values of h_can are 0.
     ! If we don't change this, no canopy will be resolved in the vertical direction.
     ! Because there is no h_can available at the moment, the next lines 
     ! are not in affect

     kcm=ke
     DO WHILE (MAXVAL( h_can - hhl(:,:,kcm) + hhl(:,:,ke1) ) > 0.0_wp)
        kcm=kcm-1
     END DO

     ! Now kcm points to the lowest layer being not a canopy layer.
     ! Therefore add +1
     kcm=kcm+1

     ! Allocate the canopy fields c_big, c_sml, r_air:
     ist = 0
     izl = 0
     ALLOCATE ( c_big(ie,je,kcm:ke1), STAT=izl ); c_big = 0.0_wp; ist = ist+izl
     ALLOCATE ( c_sml(ie,je,kcm:ke1), STAT=izl ); c_sml = 0.0_wp; ist = ist+izl
     ALLOCATE ( r_air(ie,je,kcm:ke1), STAT=izl ); r_air = 0.0_wp; ist = ist+izl

     IF (ist /= 0) THEN
       yerrmsg = 'ERROR *** Allocation of canopy fields failed ***'
       ierror  = 2014
       RETURN
     ENDIF

     CALL init_volume_canopy  (ie, je, ke, ke1, kcm,                   &
                               istartpar, iendpar, jstartpar, jendpar, &
                               fr_land, d_pat, c_big, c_sml, r_air)
     CALL init_surface_canopy (ie, je, itype_tran,                     &
                               istartpar, iendpar, jstartpar, jendpar, &
                               fr_land, plcov, lai, sai, tai, eai)
  ENDIF

  ! initialize the time stepping for the TKE 
  ! (necessary in any case, but absolutely for restarts!!)
  IF (nstart == 0) THEN
     ntke  = 0
! do not recompute, but read it from restart-file in src_input
! ELSE
!   IF (l2tls) THEN
!     ntke = nnew  !???
!   ELSE
!     ntke = nnow  !???
!   ENDIF
  ENDIF

  !Initialize varaiables of the lake model FLake
  IF (llake) CALL flake_init

  ! Diagnostic initialisation of prognostic rain / snow
  IF (ldiniprec) THEN
    ! set the timestep for the initialization phase
    zdt_orig = dt
    IF (.NOT. l2tls) THEN
      ! for Leapfrog only dt/2 is used in the first step
      dt = 0.5_wp * dt
    ENDIF
    dt2 = 2.0_wp * dt

    CALL hydci ( ldiniprec, yzerror, izerr )
    IF (izerr /= 0_iintegers) THEN
      ierror  = 2021
      yerrmsg = yzerror
    ENDIF

    ! reset the timestep again
    IF (.NOT. l2tls) THEN
      dt = zdt_orig
    ENDIF
    dt2 = 2.0_wp * dt
  ENDIF

!------------------------------------------------------------------------------
! Section 2.1: Initialization of the copy to block
!------------------------------------------------------------------------------

ELSEIF (yaction == 'init_copy') THEN


   IF (lgsp) CALL gscp_init_copy

!------------------------------------------------------------------------------
! Section 3: Physical packages at the beginning of the time stepping
!------------------------------------------------------------------------------

ELSEIF (yaction == 'compute') THEN

  ! get the correct timelevel
  IF (l2tls) THEN
    nx = nnow
    zdt = dt
  ELSE
    nx = nold
    zdt = dt2
  ENDIF

  IF (lbdclim) THEN
    ! the canopy layer has to be initialized every step, because the
    ! leaf area index, the plant cover and the root depth are changing
    IF (ltur .OR. lsoil)                                                &
      CALL init_surface_canopy (ie, je, itype_tran,                     &
                                istartpar, iendpar, jstartpar, jendpar, &
                                fr_land, plcov, lai, sai, tai, eai)
  ENDIF

  !----------------------------------------------------------------------------
  ! Section 3.1: Precipitation
  !----------------------------------------------------------------------------

  IF (ltime) CALL get_timings (i_add_computations, ntstep, dt, izerror)

! IF (lgsp) THEN
!   IF (izdebug > 5) THEN
!     PRINT *, '      GRID SCALE PRECIPITATION'
!   ENDIF

!   CALL organize_gscp
!   IF (ltime) CALL get_timings (i_precipitation, ntstep, dt, izerror)
! ENDIF

  !----------------------------------------------------------------------------
  ! Section 3.2: Radiation
  !----------------------------------------------------------------------------

  IF (lrad) THEN
!    IF ( (ntstep < 2) .OR. (ntstep == nextrad) ) THEN
      IF (izdebug > 5) THEN
        PRINT *, '      RADIATION'
      ENDIF

      CALL organize_radiation (ydate_ini)

!!$      IF ((ntstep >= 1) .OR. (nincrad == 1)) THEN
!!$        IF (hincrad /= 0.0) THEN
!!$          hnextrad = hnextrad + hincrad
!!$          nextrad  = NINT (hnextrad * 3600.0_wp / dt)
!!$        ELSE
!!$          nextrad  = nextrad + nincrad
!!$        ENDIF
!!$      ENDIF
      IF (ltime) CALL get_timings (i_radiation, ntstep, dt, izerror)
!    ENDIF
    
  ENDIF !lrad

  !----------------------------------------------------------------------------
  ! Section 3.3: Sub-grid scale orography
  !----------------------------------------------------------------------------
       
  IF (lsso) THEN
    IF ( (ntstep <= 10) .OR. (MOD(ntstep+1,nincsso) == 0) ) THEN
      IF (izdebug > 5) THEN
        PRINT *, '      SUB-GRID SCALE OROGRAPHY'
      ENDIF

      CALL organize_sso

      IF (ltime) CALL get_timings (i_sso, ntstep, dt, izerror)
    ENDIF
  ENDIF

  !----------------------------------------------------------------------------
  ! Section 3.4: Turbulence
  !----------------------------------------------------------------------------

  IF (ltur) THEN

    CALL organize_turbulence (izerror, yerrmsg)

    IF (ltime) CALL get_timings (i_turbulence, ntstep, dt, izerror)

  ENDIF

  !----------------------------------------------------------------------------
  ! Section 3.5: Soil Model  ( call of 1st part of 2-layer soil model or the
  !                            complete multi layer soil model respectively)
  !              Sea ice model
  !              Lake model
  !----------------------------------------------------------------------------

  IF (lsoil) THEN
#ifdef COUP_OAS_COS
  CALL send_fld_2clm
  CALL receive_fld_2clm
#else
    IF (lmulti_layer) THEN
      IF (izdebug > 5) THEN
        PRINT *, '      MULTI LAYER SOIL MODEL'
      ENDIF

      CALL terra_multlay (yzerror, izerr)
    ELSE
      IF (izdebug > 5) THEN
        PRINT *, '      TWO LAYER SOIL MODEL; 1st part'
      ENDIF

      CALL terra1
    ENDIF
    IF (izerr /= 0_iintegers) THEN
      ierror  = 2009
      yerrmsg = yzerror
      RETURN
    ENDIF
#endif

    IF (lseaice) THEN
      IF (izdebug > 5) THEN
        PRINT *, '      SEA ICE MODEL'
      ENDIF

      CALL seaice
    ENDIF

    IF (llake) THEN
      IF (izdebug > 5) THEN
        PRINT *, '      FLAKE MODEL'
      ENDIF

      CALL flake_interface
    ENDIF

    IF (ltime) CALL get_timings (i_soil_model, ntstep, dt, izerror)
  ENDIF
  
  !----------------------------------------------------------------------------
  ! Section 3.6: Convection
  !----------------------------------------------------------------------------

  IF (lconv) THEN
    IF ( (ntstep < 2) .OR. (MOD(ntstep+1,nincconv) == 0) ) THEN
      lzconv = .TRUE.

      SELECT CASE (itype_conv)

      CASE (0) ! Tiedtke scheme
        IF (izdebug > 5) THEN
          PRINT *, '      TIEDTKE CONVECTION SCHEME'
        ENDIF
        CALL organize_conv_tiedtke

!no more:      CASE (1) ! Kain-Fritsch scheme
!no more:        IF (izdebug > 5) THEN
!no more:          PRINT *, '      KAIN-FRITSCH CONVECTION SCHEME'
!no more:        ENDIF
!no more:        CALL organize_conv_kainfri

#ifdef CLM
      CASE (2) ! ECMWF IFS scheme
        IF (izdebug > 5) THEN
          PRINT *, '      ECMWF IFS CONVECTION SCHEME'
        ENDIF
        CALL organize_conv_ifs
#endif

      CASE (3) ! Shallow convection scheme
        IF (izdebug > 5) THEN
          PRINT *, '      SHALLOW CONVECTION SCHEME'
        ENDIF
        CALL organize_conv_shallow

      CASE DEFAULT
        ierror = 2008
        yerrmsg = 'No valid convection scheme'
        RETURN

      END SELECT

    ENDIF

#ifdef MESSY
    CALL messy_convec
#endif

    IF (ltime) CALL get_timings (i_convection, ntstep, dt, izerror)
  ENDIF
  
  !----------------------------------------------------------------------------
  ! Section 3.7: Soil Model  ( 2nd part of the 2-layer soil model)
  !----------------------------------------------------------------------------

  IF (lsoil) THEN
#ifdef COUP_OAS_COS

#else
    IF (.NOT. lmulti_layer) THEN
      IF (izdebug > 5) THEN
        PRINT *, '      TWO LAYER SOIL MODEL; 2nd part'
      ENDIF

      CALL terra2
    ENDIF
#endif

    IF (ltime) CALL get_timings (i_soil_model, ntstep, dt, izerror)
  ENDIF

  !----------------------------------------------------------------------------
  ! Section 3.8: Boundary exchange for the parallel program
  !----------------------------------------------------------------------------

  IF (izdebug > 5) THEN
    PRINT *, '      BOUNDARY EXCHANGE AFTER PHYSICS'
  ENDIF

  IF (ltime) THEN
    IF (ltime_barrier) THEN
      CALL comm_barrier (icomm_cart, izerror, yerrmsg)
      CALL get_timings (i_barrier_waiting_phy, ntstep, dt, izerror)
    ENDIF
  ENDIF

  IF (lsso) THEN
    IF (.NOT. lzconv) THEN
      kzdims(1:24)=(/ke1,ke1,ke,ke,ke,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
      CALL exchg_boundaries                                                  &
         (nx+73,  sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute,  &
          ie, je, kzdims, jstartpar, jendpar, nbl_exchg, nboundlines,        &
          my_cart_neigh, lperi_x, lperi_y, l2dim,                            &
          10000+nexch_tag, ldatatypes, ncomm_type, izerror, yerrmsg,         &
          tkvm(:,:,:), tkvh(:,:,:), qrs(:,:,:), ut_sso(:,:,:), vt_sso(:,:,:),&
          t_g(:,:,nx), qv_s(:,:,nx), tcm(:,:), tch(:,:))
    ELSE
      kzdims(1:24)=(/ke1,ke1,ke,ke,ke,ke,ke,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0/)
      CALL exchg_boundaries                                                  &
         (nx+70,  sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute,  &
          ie, je, kzdims, jstartpar, jendpar, nbl_exchg, nboundlines,        &
          my_cart_neigh, lperi_x, lperi_y, l2dim,                            &
          10000+nexch_tag, ldatatypes, ncomm_type, izerror, yerrmsg,         &
          tkvm(:,:,:), tkvh(:,:,:), qrs(:,:,:), ut_conv(:,:,:),              &
          vt_conv(:,:,:), ut_sso(:,:,:), vt_sso(:,:,:), t_g(:,:,nx),         &
          qv_s(:,:,nx), tcm(:,:), tch(:,:))
    ENDIF
  ELSE
    IF (.NOT. lzconv) THEN
      kzdims(1:24)=(/ke1,ke1,ke,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
      CALL exchg_boundaries                                                  &
         (nx+13,  sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute,  &
          ie, je, kzdims, jstartpar, jendpar, nbl_exchg, nboundlines,        &
          my_cart_neigh, lperi_x, lperi_y, l2dim,                            &
          10000+nexch_tag, ldatatypes, ncomm_type, izerror, yerrmsg,         &
          tkvm(:,:,:), tkvh(:,:,:), qrs(:,:,:), t_g(:,:,nx), qv_s(:,:,nx),   &
          tcm(:,:), tch(:,:))
    ELSE
!CPS Added tcw for exchange
#ifdef COUP_OAS_COS
      kzdims(1:24)=(/ke1,ke1,ke,ke,ke,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0/)
      CALL exchg_boundaries                                                  &
         (nx+10,  sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute,  &
          ie, je, kzdims, jstartpar, jendpar, nbl_exchg, nboundlines,        &
          my_cart_neigh, lperi_x, lperi_y, l2dim,                            &
          10000+nexch_tag, ldatatypes, ncomm_type, izerror, yerrmsg,         &
          tkvm(:,:,:), tkvh(:,:,:), qrs(:,:,:), ut_conv(:,:,:),              &
          vt_conv(:,:,:),t_g(:,:,nx),qv_s(:,:,nx),tcm(:,:),tch(:,:),         &
          tcw(:,:),trad_clm(:,:))
#else
      kzdims(1:24)=(/ke1,ke1,ke,ke,ke,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
      CALL exchg_boundaries                                                  &
         (nx+10,  sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute,  &
          ie, je, kzdims, jstartpar, jendpar, nbl_exchg, nboundlines,        &
          my_cart_neigh, lperi_x, lperi_y, l2dim,                            &
          10000+nexch_tag, ldatatypes, ncomm_type, izerror, yerrmsg,         &
          tkvm(:,:,:), tkvh(:,:,:), qrs(:,:,:), ut_conv(:,:,:),              &
          vt_conv(:,:,:), t_g(:,:,nx), qv_s(:,:,nx), tcm(:,:), tch(:,:))
#endif
    ENDIF
  ENDIF

  IF ( ltur ) THEN

    IF ( l3dturb ) THEN
      kzdims(1:24)=(/ke1,ke1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
      CALL exchg_boundaries                                                  &
           (0,  sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute,    &
           ie, je, kzdims, jstartpar, jendpar, nbl_exchg, nboundlines,       &
           my_cart_neigh, lperi_x, lperi_y, l2dim,                           &
           15000+nexch_tag, .FALSE., ncomm_type, izerror, yerrmsg,           &
           tkhm(:,:,:), tkhh(:,:,:) )
    END IF
    
    SELECT CASE ( itype_turb )
    CASE (3)
      IF ( l3dturb .OR. lprog_tke ) THEN   ! Not sure if necessary in all cases or just in 3D-Turb.
        kzdims(1:24)=(/ke1,ke1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
        CALL exchg_boundaries                                                &
             (0, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute,   &
              ie, je, kzdims, jstartpar, jendpar, nbl_exchg, nboundlines,    &
              my_cart_neigh, lperi_x, lperi_y, l2dim,                        &
              15000+nexch_tag, .FALSE., ncomm_type, izerror, yerrmsg,        &
              tke(:,:,:,ntke), tketens(:,:,:) )
      END IF
    CASE (5:8)
      IF ( lprog_tke ) THEN
        kzdims(1:24)=(/ke1,ke1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
        CALL exchg_boundaries                                                &
             (0, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute,   &
              ie, je, kzdims, jstartpar, jendpar, nbl_exchg, nboundlines,    &
              my_cart_neigh, lperi_x, lperi_y, l2dim,                        &
              15000+nexch_tag, .FALSE., ncomm_type, izerror, yerrmsg,        &
              tke(:,:,:,nnow), tketens(:,:,:) )
      ELSE
        IF ( l3dturb ) THEN   ! Not sure if necessary in all cases or just in 3D-Turb.
          kzdims(1:24)=(/ke1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
          CALL exchg_boundaries                                              &
               (0, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute, &
                ie, je, kzdims, jstartpar, jendpar, nbl_exchg, nboundlines,  &
                my_cart_neigh, lperi_x, lperi_y, l2dim,                      &
                15000+nexch_tag, .FALSE., ncomm_type, izerror, yerrmsg,      &
                tke(:,:,:,1) )
        END IF
      END IF
    END SELECT

  END IF
    
  IF (nradcoarse > 1) THEN
    ! radiation computations on a coarser grid need additional
    ! boundary exchange

    IF ( lsoil ) THEN
      IF(lmulti_layer) THEN
        kzdims(1:24)=(/1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
       CALL exchg_boundaries                                                 &
         (70+nnew, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute, &
          ie, je, kzdims, jstartpar, jendpar, nradcoarse-1, nboundlines,     &
          my_cart_neigh, lperi_x, lperi_y, l2dim,                            &
          17000+nexch_tag, ldatatypes, ncomm_type, izerror, yerrmsg,         &
          t_g(:,:,nnew),w_so(:,:,1,nnew),t_snow(:,:,nnew),w_snow(:,:,nnew),  &
          freshsnow(:,:) )
      ELSE
        kzdims(1:24)=(/1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
        CALL exchg_boundaries                                                &
          (73+nnew, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute,&
           ie, je, kzdims, jstartpar, jendpar, nradcoarse-1, nboundlines,    &
           my_cart_neigh, lperi_x, lperi_y, l2dim,                           &
           17000+nexch_tag, ldatatypes, ncomm_type, izerror, yerrmsg,        &
           t_g(:,:,nnew),w_g1(:,:,nnew),t_snow(:,:,nnew),w_snow(:,:,nnew))
      ENDIF
    ELSE ! .NOT.lsoil:
        kzdims(1:24)=(/1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
        ! exchange of t_snow probably not necessary
        CALL exchg_boundaries                                                &
          (76+nnew, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute,&
           ie, je, kzdims, jstartpar, jendpar, nradcoarse-1, nboundlines,    &
           my_cart_neigh, lperi_x, lperi_y, l2dim,                           &
           17000+nexch_tag, ldatatypes, ncomm_type, izerror, yerrmsg,        &
           t_g(:,:,nnew),t_snow(:,:,nnew) )
    ENDIF

    IF (lzconv) THEN
      kzdims(1:24)=(/ke,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
      CALL exchg_boundaries                                                  &
        (2, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute,        &
         ie, je, kzdims, jstartpar, jendpar, nradcoarse-1, nboundlines,      &
         my_cart_neigh, lperi_x, lperi_y, l2dim,                             &
         19000+nexch_tag, ldatatypes, ncomm_type, izerror, yerrmsg,          &
         clc_con(:,:,:))
    ENDIF

  ENDIF   ! nradcoarse > 1

  IF (ltime) CALL get_timings (i_communications_phy, ntstep, dt, izerror)

  !----------------------------------------------------------------------------
  ! Section 3.9: Computation of total physical forcings
  !              (this has been in the dynamical cores before Version 4.23)
  !----------------------------------------------------------------------------

  ! retrieve the required microphysics tracers
  CALL trcr_get(ierror, idt_qv, ptr_tens=qv_tens)
  IF (izerror /= 0) THEN
    yerrmsg = trcr_errorstr(ierror)
    CALL model_abort(my_cart_id, ierror, yerrmsg, yzroutine)
  ENDIF
  CALL trcr_get(ierror, idt_qc, ptr_tens=qc_tens)
  IF (ierror /= 0) THEN
    yerrmsg = trcr_errorstr(ierror)
    CALL model_abort(my_cart_id, ierror, yerrmsg, yzroutine)
  ENDIF
  CALL trcr_get(ierror, idt_qi, ptr_tens=qi_tens)
  IF (ierror /= 0 .AND. ierror /= T_ERR_NOTFOUND) THEN
    yerrmsg = trcr_errorstr(ierror)
    CALL model_abort(my_cart_id, ierror, yerrmsg, yzroutine)
  ENDIF
  
  ! Reset error
  ierror = 0_iintegers

  ! Compute total physical forcings (tendencies)
  DO  k = 1 , ke
#ifdef SCLM
   IF (lscm .AND. lradtnd) THEN
#endif       
    DO    j = jstart, jend 
      DO  i = istart, iend
        ttens  (i,j,k) = ttens  (i,j,k) + sohr(i,j,k) + thhr(i,j,k)
        ttens  (i,j,k) = ttens  (i,j,k) + tt_conv (i,j,k)
        qv_tens(i,j,k) = qv_tens(i,j,k) + qvt_conv(i,j,k)
      ENDDO
    ENDDO
#ifdef SCLM
   ENDIF
#endif

    IF (lsso) THEN
      DO    j = jstart, jend
        DO  i = istart, iend
          ttens(i,j,k) = ttens(i,j,k) + tt_sso(i,j,k)
        ENDDO
      ENDDO
    ENDIF

!US this is not active at the moment, but might be reconsidered in the future
!US therefore we leave these lines in the code
!   IF (itype_conv == 1) THEN
!     DO    j = jstart, jend
!       DO  i = istart, iend
!         qc_tens(i,j,k) = qc_tens(i,j,k) + qct_conv(i,j,k)
!       ENDDO
!     ENDDO
!   ENDIF

    IF ( (itype_conv == 0) .OR. (itype_conv == 2) ) THEN
      ! Add convective tendencies of qc and qi computed by Tiedtke(0) or
      ! IFS scheme (2)
      DO    j = jstart, jend
        DO  i = istart, iend
          qc_tens(i,j,k) = qc_tens(i,j,k) + qct_conv(i,j,k)
        ENDDO
      ENDDO
      IF ( ASSOCIATED(qi_tens) ) THEN
        DO    j = jstart, jend
          DO  i = istart, iend
            qi_tens(i,j,k) = qi_tens(i,j,k) + qit_conv(i,j,k)
          ENDDO
        ENDDO
      ENDIF
    ENDIF

    ! SSO and convective tendencies for u and v are defined on the mass grid points
    ! and have to remain there for the output. And ut_sso, vt_sso are used later on
    ! also in the turbulence scheme on the mass grid points!
    ! only for adding the tendencies to utens, vtens, they are interpolated to the
    ! u- and v-grid points

    IF (lconv) THEN
      DO   j = jstartu, jendu
        DO i = istartu, iendu
          utens(i,j,k) = utens(i,j,k) + 0.5_wp*( ut_conv(i+1,j,k) + ut_conv(i,j,k) )
        ENDDO
      ENDDO
      DO   j = jstartv, jendv
        DO i = istartv, iendv
          vtens(i,j,k) = vtens(i,j,k) + 0.5_wp*( vt_conv(i,j+1,k) + vt_conv(i,j,k) )
        ENDDO
      ENDDO
    ENDIF
    IF (lsso) THEN
      DO   j = jstartu, jendu
        DO i = istartu, iendu
          utens(i,j,k) = utens(i,j,k) + 0.5_wp*( ut_sso(i+1,j,k) + ut_sso(i,j,k) )
        ENDDO
      ENDDO
      DO   j = jstartv, jendv
        DO i = istartv, iendv
          vtens(i,j,k) = vtens(i,j,k) + 0.5_wp*( vt_sso(i,j+1,k) + vt_sso(i,j,k) )
        ENDDO
      ENDDO
    ENDIF

  ENDDO

  IF (izdebug > 5) THEN
    PRINT *, '      PHYSICAL COMPUTATIONS before dynamics finished '
  ENDIF

!------------------------------------------------------------------------------
! Section 4: Physical packages at the end of the time stepping
!------------------------------------------------------------------------------

ELSEIF (yaction == 'finish_compute') THEN

  IF (lgsp) THEN
    IF (izdebug > 5) THEN
      PRINT *, '      GRID SCALE PRECIPITATION after dynamics'
    ENDIF

    IF     (itype_gscp <= 4) THEN
      ! blocked version of microphysics: routines are called within the interface
      CALL gscp_organize (ierror, yerrmsg)
      IF (ierror /= 0) RETURN
#ifdef TWOMOM_SB
    ELSEIF (itype_gscp >= 100 ) THEN
      ! this is still in the old data structure
      CALL seifert_pp
#endif
    ENDIF

!US old versions of the microphysics, still present in src_gscp
! can be activated again just by un-commenting and commening the above calls
!   IF     (itype_gscp == 4) THEN
!     CALL hydci_pp_gr
!   ELSEIF (itype_gscp == 3) THEN
!     CALL hydci_pp
!   ELSEIF (itype_gscp == 2) THEN
!     CALL hydor_pp
!   ELSEIF (itype_gscp == 1) THEN
!     CALL kessler_pp
!ifdef TWOMOM_SB
!   ELSEIF (itype_gscp >= 100 ) THEN
!     CALL seifert_pp
!endif
!   ENDIF

    IF (izdebug > 5) THEN
      PRINT *, '      PHYSICAL COMPUTATIONS after dynamics finished '
    ENDIF

    IF (ltime) CALL get_timings (i_precipitation, ntstep, dt, izerror)
  ENDIF

!------------------------------------------------------------------------------
! Section 5: All other actions are wrong
!------------------------------------------------------------------------------

ELSE

  ierror  = 1
  yerrmsg = 'ERROR *** No valid action for the physics ***'

ENDIF

!------------------------------------------------------------------------------
! Internal Procedures
!------------------------------------------------------------------------------

CONTAINS

!==============================================================================
!+ Module procedure in "setup" for the input of NAMELIST phyctl
!------------------------------------------------------------------------------

SUBROUTINE input_phyctl (nuspecif, nuin, ierrstat)

!------------------------------------------------------------------------------
! Description:
!   This subroutine organizes the input of the NAMELIST-group phyctl. 
!   The group phyctl contains variables for the organization of the physics.
!   These are logical variables whether a certain package has to be performed
!   and organizational variables that determine how often it is performed.
!
! Method:
!   All variables are initialized with default values and then read in from
!   the file INPUT. The input values are checked for errors and for
!   consistency. If wrong input values are detected the program prints
!   an error message. The program is not stopped in this routine but an
!   error code is returned to the calling routine that aborts the program after
!   reading in all other namelists.
!   In parallel mode, the variables are distributed to all nodes with the
!   environment-routine distribute_values.
!   Both, default and input values are written to the file YUSPECIF
!   (specification of the run).
!
!------------------------------------------------------------------------------

! Subroutine / Function arguments

  INTEGER   (KIND=iintegers),   INTENT (IN)      ::        &
    nuspecif,     & ! Unit number for protocolling the task
    nuin            ! Unit number for Namelist INPUT file

  INTEGER   (KIND=iintegers),   INTENT (INOUT)   ::        &
    ierrstat        ! error status variable

! Variables for default values
  REAL (KIND=wp)             ::       &
    hincrad_d,    & ! hour increment for running the radiation
    czbot_w_so_d, & ! depth of bottom of last hydrological active soil layer
    czml_soil_d(20), & ! depth of the main level of the soil layers (defaults)
    czml_soil  (20)    ! depth of the main level of the soil layers (read in)


  INTEGER (KIND=iintegers)   ::       &
    nradcoarse_d, & ! number of horiz. gridpoints for radiation on coarser grid
    nincrad_d,    & ! time step increment for running the radiation
    nincsso_d,    & ! time step increment for running the radiation
    ninctura_d,   & ! time step increment for running the vertical diffusion
    nincconv_d,   & ! time step increment for running the convection scheme 

    itype_trvg_d, & ! type of vegetation transpiration parameterization
    itype_evsl_d, & ! type of parameterization of bare soil evaporation

    itype_gscp_d, & ! type of grid-scale precipitation physics               

    itype_sher_d, & ! type of shear production for TKE
    itype_wcld_d, & ! type of water cloud diagnosis
    itype_tran_d, & ! type of surface-atmosphere transfer
    itype_turb_d, & ! type of turbulent diffusion parametrization
    itype_synd_d, & ! type of diagnosis of synop. station values

    imode_tran_d, & ! mode of surface-atmosphere transfer
    imode_turb_d, & ! mode of turbulent diffusion parametrization

    ico2_rad_d,   & ! type of CO2 concentration in radiation parameterization
                    ! (used in the CLM)
    iy_co2_stab_d,& ! default stabilisation year
#ifdef TWOMOM_SB
    iradpar_cloud_d, & ! type of parameterization of radiative transfer parameters (ext., sing. alb., asym.)
#endif
    icldm_rad_d,  & ! mode of cloud representation in radiation  parametr.
    icldm_tran_d, & ! mode of cloud representation in transfer  parametr.
    icldm_turb_d, & ! mode of cloud representation in turbulence parametr.
    itype_conv_d, & ! type of convection parameterization

    itype_aerosol_d, &! type of aerosol map
    itype_root_d,    &! type of root density distribution
    itype_heatcond_d,&! type of heat conductivity
    itype_hydbound_d,&! type of hydraulic lower boundary

    nhori_d,      & ! number of sectors for the horizont by the topographic radcorr
    ke_soil_d,    & ! number of layers in multi-layer soil model
    ke_snow_d,    & ! number of layers in multi-layer snow model
    nlgw_d,       & ! number of prognostic soil water levels
    itype_albedo_d  ! type of soil albedo treatment

  LOGICAL                    ::       &
    lrad_d,       & ! forecast with radiation
    lforest_d,    & ! to run with forest data (evergreen and deciduous)
    ltur_d,       & ! forecast with vertical diffusion
    lconv_d,      & ! forecast with convection
    lconv_inst_d, & ! output of instantaneous values of top_con/bas_con
                    ! instead of min/max for an output interval
    lgsp_d,       & ! forecast with grid scale precipitation
    lsuper_coolw_d,&! switch for supercooled liquid water (work from Felix Rieper)
    ldiniprec_d,  & ! diagnostic initialisation of prognostic precip (qr, qs)
    l3dturb_d,    & ! 3D-turbulence: CALL explicit_horizontal_diffusion (RK)
    l3dturb_metr_d,&! switch on/off additional metric terms for 3D-turbulence
    lprog_tke_d,  & ! prognostic treatment of TKE (for itype_turb=5/7)
    limpltkediff_d,&! use semi_implicit TKE-diffusion
    lsoil_d,      & ! forecast with soil model
    lmelt_d,      & ! soil model with melting process
    lmelt_var_d,  & ! freezing temperature dependent on water content
    lmulti_layer_d,&! run multi-layer soil model
    lmulti_snow_d, &! run multi-layer snow model
    lseaice_d,    & ! forecast with sea ice model
    llake_d,      & ! forecast with lake model
    lsso_d,       & ! forecast with sub-grid scale orography scheme

    lemiss_d,     & ! external surface emissivity map
    lstomata_d,   & ! external minimum stomata resistance

    ltkesso_d,    & ! calculation SSO-wake turbulence production for TKE
    ltkecon_d,    & ! consider convective buoyancy production for TKE
    lexpcor_d,    & ! explicit corrections of the implicit calculated
                    ! turbulent diffusion (only if itype_turb=3)
    ltmpcor_d,    & ! consideration of thermal TKE-sourcel in the enthalpy budget
    lprfcor_d,    & ! using the profile values of the lowest main level instead of
                    ! the mean value of the lowest layer for surface flux calulations
    lnonloc_d,    & ! nonlocal calculation of vertical gradients used
                    ! for turbulent diffusion (only if itype_turb=3)
    lcpfluc_d,    & ! consideration of fluctuations of the heat capacity of air
    lconf_avg_d,  & ! average convective forcings in case of massflux closure
    lradf_avg_d,  & ! average radiative forcings when running on coarser grid
    lcape_d,      & ! convection with CAPE closure
    lctke_d,      & ! convection with turbulent convective energy closure
                    ! warning: lctke not yet fully implemented
    lradtopo_d,   & ! uses topographic correction of radiation
    lco2_stab_d     ! a  year for CO2 stabilisation will be defined

  INTEGER (KIND=iintegers)   :: i, k, invar, ierr, iz_err
  LOGICAL                    :: lzequiv

  CHARACTER(LEN=250)         :: iomsg_str

! Define the namelist group

  NAMELIST /phyctl/ lrad, ltur, lconv, itype_conv, lgsp, lsuper_coolw,        &
                    lsoil, lmelt, lmelt_var, lmulti_layer, lexpcor,           &
                    ltmpcor, lprfcor, lnonloc, lcpfluc,lcape,lctke, lconf_avg,&
                    nincrad, hincrad, ninctura, nincconv,                     &
                    ldiniprec, itype_trvg, itype_evsl,                        &
                    itype_gscp, itype_wcld, itype_tran, itype_turb,itype_synd,&
                    icldm_rad, icldm_tran, icldm_turb, imode_tran, imode_turb,&
                    ke_soil, czml_soil, nlgw, l3dturb, lprog_tke,             &
                    lforest, lconv_inst, lseaice, llake, ico2_rad, czbot_w_so,&
                    l3dturb_metr, nradcoarse, lradf_avg, lradtopo, nhori,     &
                    lsso, nincsso, limpltkediff, ltkesso, ltkecon,            &
                    itype_sher, lmulti_snow, ke_snow, lemiss, lstomata,       &
                    itype_aerosol, itype_root, itype_heatcond, itype_hydbound,&
#ifdef TWOMOM_SB
                    iradpar_cloud,                                            &
#endif
                    lco2_stab, iy_co2_stab, itype_albedo

!------------------------------------------------------------------------------
!- End of header -
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!- Begin SUBROUTINE input_phyctl
!------------------------------------------------------------------------------

ierrstat = 0_iintegers
iz_err   = 0_iintegers

IF (my_world_id == 0) THEN

!------------------------------------------------------------------------------
!- Section 1: Initialize the default variables
!------------------------------------------------------------------------------

  lrad_d         = .TRUE.
  lforest_d      = .TRUE.
  ltur_d         = .TRUE.
  lconv_d        = .TRUE.
  lconv_inst_d   = .FALSE.
  lgsp_d         = .TRUE.
  lsuper_coolw_d = .FALSE.
  ldiniprec_d    = .FALSE.
  l3dturb_d      = .FALSE.
  l3dturb_metr_d = .TRUE.
  lprog_tke_d    = .FALSE.
  limpltkediff_d = .TRUE.
  lsoil_d        = .TRUE.
  lmelt_d        = .TRUE.
  lmelt_var_d    = .TRUE.
  lmulti_layer_d = .TRUE.
  lmulti_snow_d  = .FALSE.
  lseaice_d      = .TRUE.
  llake_d        = .FALSE.
  lemiss_d       = .FALSE.
  lstomata_d     = .FALSE.

  lsso_d         = .TRUE.
  ltkesso_d      = .FALSE.
  ltkecon_d      = .FALSE.
  lexpcor_d      = .FALSE.
  ltmpcor_d      = .FALSE.
  lprfcor_d      = .FALSE.
  lnonloc_d      = .FALSE.
  lcpfluc_d      = .FALSE.
  lconf_avg_d    = .TRUE. 
  lradf_avg_d    = .FALSE.
  lcape_d        = .FALSE.
  lctke_d        = .FALSE.
  lradtopo_d     = .FALSE.
  lco2_stab_d    = .FALSE.

  nhori_d      = 24
  nradcoarse_d = 1
  nincrad_d    = 0
  hincrad_d    = 1.0_wp
  ninctura_d   = 1
  nincconv_d   = 4
  nincsso_d    = 5

  itype_trvg_d = 2 ! vegetation transpiration using the BATS-approach
  itype_evsl_d = 2 ! bare soil evaporation using the BATS-approach

  itype_gscp_d = 3 ! grid-scale-cloud-precipitation with 'cloudice'

  itype_sher_d = 1 ! only vertical shear
  itype_wcld_d = 2 ! water cloud diagnosis (in SUB cloud_diag) using statistical condensation
  itype_tran_d = 2 ! new surface-atmosphere transfer using SUB turb_tran
  itype_turb_d = 3 ! new moist scheme with prognostic tke-equation using SUB turb_diff
  itype_synd_d = 2 ! new method with resistance formulation using SUB turb_tran

  imode_tran_d = 1 ! diagnostik TKE in 'turbtran'
  imode_turb_d = 1 ! implicite calculation of turbulent diffusion using the explicit
                   ! surface flux densities as lower boundary condition

  ico2_rad_d   = 0 ! constant CO2 concentration (330 ppm)
  iy_co2_stab_d = 2001 ! year of CO2 stabilisation

#ifdef TWOMOM_SB
  iradpar_cloud_d = 1 ! old method, optical properties of clouds only dep. on LWC / IWC
#endif

  icldm_rad_d  = 4 ! special (old) cloud diagnosis (using relative humidity) for radiation
  icldm_tran_d = 0 ! surface-atmosphere-transfer considers no water clouds
  icldm_turb_d = 2 ! turbulence considers water clouds calculated in SUB cloud_diag
                   ! dependent on itype_wcld.

  itype_conv_d = 0 ! Tiedtke convection scheme

  itype_aerosol_d   = 1 ! fixed aerosol map
  itype_root_d      = 1 ! uniform root density distribution
  itype_heatcond_d  = 1 ! average soil moisture for heat conductivity
  itype_hydbound_d  = 1 ! drainage no diffusion

  ke_soil_d      = 7
  czml_soil_d(1)  = 0.005_wp
  czml_soil_d(2)  = 0.02_wp
  czml_soil_d(3)  = 0.06_wp
  czml_soil_d(4)  = 0.18_wp
  czml_soil_d(5)  = 0.54_wp
  czml_soil_d(6)  = 1.62_wp
  czml_soil_d(7)  = 4.86_wp
  czml_soil_d(8)  =14.58_wp
  czml_soil_d(9:) = -99.9999_wp
  czbot_w_so_d    = 2.5_wp

  ke_snow_d      = 2
  nlgw_d         = 2

  itype_albedo_d = 1

!------------------------------------------------------------------------------
!- Section 2: Initialize variables with defaults
!------------------------------------------------------------------------------

  lrad           = lrad_d
  lforest        = lforest_d
  ltur           = ltur_d
  lconv          = lconv_d
  lconv_inst     = lconv_inst_d
  lgsp           = lgsp_d
  lsuper_coolw   = lsuper_coolw_d
  ldiniprec      = ldiniprec_d
  l3dturb        = l3dturb_d
  l3dturb_metr   = l3dturb_metr_d
  lprog_tke      = lprog_tke_d
  limpltkediff   = limpltkediff_d
  lsoil          = lsoil_d
  lmelt          = lmelt_d
  lmelt_var      = lmelt_var_d
  lmulti_layer   = lmulti_layer_d
  lmulti_snow    = lmulti_snow_d
  lseaice        = lseaice_d
  llake          = llake_d
  lemiss         = lemiss_d
  lstomata       = lstomata_d

  lsso           = lsso_d
  ltkesso        = ltkesso_d
  ltkecon        = ltkecon_d
  lexpcor        = lexpcor_d
  ltmpcor        = ltmpcor_d
  lprfcor        = lprfcor_d
  lnonloc        = lnonloc_d
  lcpfluc        = lcpfluc_d
  lconf_avg      = lconf_avg_d
  lradf_avg      = lradf_avg_d
  lcape          = lcape_d  
  lctke          = lctke_d  
  lradtopo       = lradtopo_d
  lco2_stab      = lco2_stab_d !HJP 2011-12-19

  nhori          = nhori_d
  nradcoarse     = nradcoarse_d
  nincrad        = nincrad_d
  ninctura       = ninctura_d
  nincconv       = nincconv_d
  nincsso        = nincsso_d
  hincrad        = hincrad_d

  itype_trvg     = itype_trvg_d
  itype_evsl     = itype_evsl_d

  itype_gscp     = itype_gscp_d

  itype_sher     = itype_sher_d
  itype_wcld     = itype_wcld_d
  itype_tran     = itype_tran_d
  itype_turb     = itype_turb_d
  itype_synd     = itype_synd_d

  imode_tran     = imode_tran_d
  imode_turb     = imode_turb_d

  ico2_rad       = ico2_rad_d
  iy_co2_stab    = iy_co2_stab_d

#ifdef TWOMOM_SB
  iradpar_cloud  = iradpar_cloud_d
#endif

  icldm_rad      = icldm_rad_d
  icldm_tran     = icldm_tran_d
  icldm_turb     = icldm_turb_d

  itype_conv     = itype_conv_d

  itype_aerosol  =  itype_aerosol_d
  itype_root     =  itype_root_d
  itype_heatcond =  itype_heatcond_d
  itype_hydbound =  itype_hydbound_d

  ke_soil        = ke_soil_d
  ke_snow        = ke_snow_d
  czml_soil(:)   = -1.0_wp
  czbot_w_so     = czbot_w_so_d
  nlgw           = nlgw_d

  itype_albedo   = itype_albedo_d

!------------------------------------------------------------------------------
!- Section 3: Input of the namelist values
!------------------------------------------------------------------------------

  iomsg_str(:) = ' '
  READ (nuin, phyctl, IOSTAT=iz_err, IOMSG=iomsg_str)

  IF (iz_err /= 0) WRITE (*,'(A,A)') 'Namelist-ERROR PHYCTL: ', TRIM(iomsg_str)
ENDIF

IF (nproc > 1) THEN
  ! distribute error status to all processors
  CALL distribute_values  (iz_err, 1, 0, imp_integers,  icomm_world, ierr)
ENDIF

IF (iz_err /= 0) THEN
  ierrstat = -1
  RETURN
ENDIF

IF (my_world_id == 0) THEN

!------------------------------------------------------------------------------
!- Section 4: Check values for errors and consistency
!------------------------------------------------------------------------------

  ! Check whether the values for the increments are given in hours 
  ! and calculate the values in time steps
  IF ( nincrad /= nincrad_d ) THEN
    ! hour values take priority over time step values
    IF ( hincrad /= hincrad_d ) THEN
      IF (hincrad <= 0.0_wp) THEN
        PRINT *, ' ERROR  *** Wrong value for hincrad: ', hincrad, ' *** '
        PRINT *, '        ***   must be > 0.0 *** '
        ierrstat = 1002
      ENDIF
      nincrad = NINT ( hincrad * 3600.0_wp/dt)
    ELSE
      hincrad = 0.0_wp
    ENDIF
  ELSE  
    IF ( hincrad /= hincrad_d ) THEN
      IF (hincrad <= 0.0_wp) THEN
        PRINT *, ' ERROR  *** Wrong value for hincrad: ', hincrad, ' *** '
        PRINT *, '        ***   must be > 0.0 *** '
        ierrstat = 1002
      ENDIF
    ENDIF
    nincrad = NINT ( hincrad * 3600.0_wp/dt)
  ENDIF
#ifndef MESSY
  hnextrad = hstart
#else
  IF (hstart == 0.0_wp) THEN
     hnextrad = hstart
  ELSE
     ! NOTE: hnextrad is restart attribute to ensure radiation trigger
     !       at exactly the same point in time
  ENDIF
#endif
  nextrad  = NINT ( hnextrad * 3600.0_wp / dt)


  IF ( (nradcoarse > nboundlines+1) .OR. (nradcoarse < 1) ) THEN
    PRINT *, ' ERROR  *** Wrong value for nradcoarse: ', nradcoarse, ' *** '
    PRINT *, '        ***   must be >= 1 and <= nboundlines+1 !*** '
    ierrstat = 1002
  ENDIF
  IF ( (nradcoarse > 1) .AND. (nboundlines > 3) ) THEN
    PRINT *, ' ERROR  *** nradcoarse > 1 and nboundlines > 3 is not working  *** '
    ierrstat = 1002
  ENDIF
  IF ( (nradcoarse < 2) .AND. lradf_avg ) THEN
    PRINT *, ' ERROR: If nradcoarse<2, lradf_avg not implemented yet. '
    ierrstat = 1002
  ELSEIF ( (nradcoarse > 1) .AND. .NOT.lradf_avg ) THEN
    PRINT *, ' WARNING: IF (nradcoarse > 1), lradf_avg==.TRUE. recommended.'
  ENDIF
  IF ( (nradcoarse > 1) .AND. lradtopo ) THEN
    PRINT *, ' ERROR: If nradcoarse > 1, lradtopo is not possible. '
    ierrstat = 1002
  ENDIF
  IF ( lperi_x .AND. lrad ) THEN
    PRINT *, ' WARNING: lrad=.true. and lperi_x=.true. => constant aerosol and zenith angle in X-dir.,'
    PRINT *, '          taken from domain reference point (as defined via pollon, pollat)!'
  ENDIF
  IF ( lperi_y .AND. lrad ) THEN
    PRINT *, ' WARNING: lrad=.true. and lperi_y=.true. => constant aerosol and zenith angle in Y-dir.,'
    PRINT *, '          taken from domain reference point (as defined via pollon, pollat)!'
  ENDIF

#ifndef TWOMOM_SB
  l_2mom = .FALSE.
  IF ( (itype_gscp < 1) .OR. (itype_gscp > 4) ) THEN
    PRINT *, ' ERROR  *** Wrong value for itype_gscp: ', itype_gscp, ' *** '
    PRINT *, '        ***   must be >= 1 and <= 4 !!!                  *** '
    PRINT *, ' NOTE:  *** No 2-moment scheme compiled in this binary   *** '
    ierrstat = 1002
  ENDIF
#else
  IF ( (itype_gscp < 1) .OR. ( itype_gscp > 4 .and. itype_gscp < 100 ) ) THEN
    PRINT *, ' ERROR  *** Wrong value for itype_gscp: ', itype_gscp, ' *** '
    PRINT *, '        ***   must be >= 1 and <= 4 or >= 100 !!!        *** '
    PRINT *, ' NOTE:  *** 2-moment scheme is compiled in this binary   *** '
    ierrstat = 1002
  ENDIF
  IF (itype_gscp >= 100) THEN
    l_2mom = .TRUE.
    IF (dt > 30.0_wp) THEN
      itype_gscp = itype_gscp_d
      PRINT *, ' ERROR  *** itype_gscp >= 100 (2mom microphysics) should only       ***'
      PRINT *, '        ***   run with dt <= 30 seconds!                            *** '
      ierrstat = 1002
    END IF
    IF (luse_rttov) THEN
      PRINT  *, ' WARNING  *** Satellite images together with 2mom microphysics',          &
                                ' not correctly implemented yet! ***'
    END IF
  ELSE
    l_2mom = .FALSE.
  END IF

  IF ( .NOT. lrad .AND. iradpar_cloud /= iradpar_cloud_d ) THEN
    PRINT *, ' WARNING *** lrad=.false., so iradpar_cloud has been       ***'
    PRINT *, '         *** reset to its default (=1) to save comp. time  ***'
    iradpar_cloud = iradpar_cloud_d
  END IF
  IF (iradpar_cloud < 1 .OR. iradpar_cloud > 3) THEN
    PRINT *, ' ERROR  *** iradpar_cloud can only be 1, 2, or 3      ***'
    ierrstat = 1002
  END IF
  IF ( .NOT. l_2mom .AND. iradpar_cloud /= iradpar_cloud_d  ) THEN
    PRINT *, ' WARNING *** l_2mom=.false., so iradpar_cloud has been       ***'
    PRINT *, '         *** reset to its default (=1) to save comp. time  ***'
    PRINT *, '         *** and has no effect on cloud properties  ***'
    iradpar_cloud = iradpar_cloud_d
  ENDIF
#endif

  IF ( (.NOT. lgsp) .AND. (itype_gscp /= itype_gscp_d) ) THEN
    itype_gscp = itype_gscp_d
    PRINT *, ' WARNING:  *** itype_gscp is set to the default *** ',    &
             '           *** because precipitation is turned off ***'
  ENDIF
  IF ( (itype_gscp < 3) .AND. ldiniprec) THEN
    PRINT  *, ' ERROR  *** Initialization of prognostic precipitation ***'
    PRINT  *, '        *** only possible for itype_gscp >= 3          *** '
    ierrstat = 1002
  ENDIF
  IF ( (itype_gscp < 3) .AND. luse_rttov) THEN
    PRINT  *, ' ERROR  *** Without cloud ice, satellite images ',          &
                                           'cannot be computed ***'
    ierrstat = 1002
  ENDIF
  IF ((itype_gscp < 3) .AND. lsppt) THEN
    PRINT *, ' ERROR  *** itype_gscp = 1/2 is not possible with lsppt  *** '
    PRINT *, '        *** lsppt not implemented for kessler and hydor  *** '
    ierrstat = 1002
  ENDIF
  IF (ldiniprec)                                                              &
    PRINT *,' NOTE: *** ldiniprec will be set FALSE if lana_qr_qs is true *** '

  IF (itype_sher.LT.1 .OR.  itype_sher.GT.3) THEN
     itype_sher = itype_sher_d
     PRINT *,' WARNING: *** itype_sher is set to default again *** '
  END IF
  IF (itype_wcld.LT.1 .OR.  itype_wcld.GT.2) THEN
     itype_wcld = itype_wcld_d
     PRINT *,' WARNING: *** itype_wcld is set to default again *** '
  END IF
  IF (itype_tran.LT.1 .OR.  itype_tran.GT.2) THEN
     itype_tran = itype_tran_d
     PRINT *,' WARNING: *** itype_tran is set to default again *** '
  END IF

  IF (.NOT. (ltur .AND. lphys) .AND. lprog_tke) THEN
    lprog_tke  = .FALSE.
    PRINT *,' WARNING: *** lprog_tke=.true. not possible for ltur=.false. or lphys=.false.! *** '
    PRINT *,'          *** lprog_tke is set to .false.!         *** '
  END IF

  SELECT CASE( itype_turb )
  CASE( 1 )
    IF (lprog_tke) THEN
      lprog_tke  = .FALSE.
      PRINT *,' WARNING: *** lprog_tke=.true. not possible for itype_turb = 1! *** '
      PRINT *,'          *** lprog_tke is set to .false.!         *** '
    END IF
  CASE( 3 )
    IF (.NOT. l2tls .AND. lprog_tke) THEN
      lprog_tke  = .FALSE.
      PRINT *,' WARNING: *** lprog_tke=.true. only for RK-Scheme! *** '
      PRINT *,'          *** lprog_tke is set to .false.!         *** '
    END IF
  CASE( 5:8 )
    IF (.NOT. l2tls .AND. lprog_tke) THEN
      lprog_tke  = .FALSE.
      PRINT *,' WARNING: *** lprog_tke=.true. only for RK-Scheme! *** '
      PRINT *,'          *** lprog_tke is set to .false.!         *** '
    END IF
    IF (.NOT. l2tls) THEN
      itype_turb = itype_turb_d
      PRINT *,' WARNING: *** itype_turb = ',itype_turb,' only available for RK-Scheme! *** '
      PRINT *,'          *** itype_turb is set to default again *** '
    END IF
  CASE( 100 )
    lprog_tke  = .FALSE.
  CASE default
    lprog_tke  = .FALSE.
    itype_turb = itype_turb_d
    PRINT *,' WARNING: *** no valid value of itype_turb! *** '
    PRINT *,'          *** itype_turb is set to default and lprog_tke to .false. *** '
  END SELECT

  IF (itype_synd.LT.1 .OR.  itype_synd.GT.2) THEN
     itype_synd = itype_synd_d
     PRINT *,' WARNING: *** itype_synd is set to default again *** '
  END IF
  
  IF (imode_tran.LT.1 .OR.  imode_tran.GT.2) THEN
     imode_tran = imode_tran_d
     PRINT *,' WARNING: *** imode_tran is set to default again *** '
  END IF
  IF (imode_turb.LT.0 .OR.  imode_turb.GT.4) THEN
     imode_turb = imode_turb_d
     PRINT *,' WARNING: *** imode_turb is set to default again *** '
  END IF

  IF (icldm_rad.LT.0  .OR. icldm_rad.GT.4 )  THEN
     icldm_rad  = icldm_rad_d
     PRINT *,' WARNING: *** icldm_rad is set to default again *** '
  END IF
     
  IF (icldm_tran.LT.0 .OR. icldm_tran.GT.2)  THEN
     icldm_tran = icldm_tran_d
     PRINT *,' WARNING: *** icldm_tran is set to default again *** '
  END IF
  ! if icldm_turb = -1: real dry scheme: qc is not taken into account
  IF (icldm_turb.LT.-1 .OR. icldm_turb.GT.2) THEN
     icldm_turb = icldm_turb_d
     PRINT *,' WARNING: *** icldm_turb is set to default again *** '
  END IF

  IF (ico2_rad > 10)  ico2_rad  = ico2_rad_d

  ! the prognostic TKE-scheme of Matthias Raschendorfer is not used
  IF (itype_turb.NE.3) THEN
    imode_turb = 0
    itype_tran = MIN(itype_tran,2)
    itype_wcld = 1
    itype_sher = 1
    icldm_turb = 0
  END IF
  ! the Louis-scheme for transfer is used
  IF (itype_tran.EQ.1) THEN
    imode_tran = 1
    icldm_tran = 0
    itype_synd = 1
  ENDIF

  IF (.NOT. lmulti_layer) THEN
    PRINT *,' ERROR  *** lmulti_layer=.false. no longer supported    *** '
    PRINT *,' ERROR  *** If lmulti_layer=.false. is desired,         *** '
    PRINT *,' ERROR  *** you will have to switch off this error      *** '
    PRINT *,' ERROR  *** check manually in organize_physics.f90!     *** '
    ierrstat = 1099
  END IF
  
  IF (ltkesso .AND. .NOT.lsso) THEN
    PRINT *,' WARNING: *** ltkesso cannot be active, since SSO scheme is not running *** '
  END IF

  IF (ltkecon) THEN
    IF (.NOT. lconv) THEN
      PRINT *,' WARNING: *** ltkecon cannot be active, since convection is not running *** '
    ELSEIF (itype_conv == 2) THEN
      PRINT *,' WARNING: *** ltkecon is not yet supported by convection according to IFS *** '
    END IF
  END IF

  IF (lmulti_layer) THEN
    ALLOCATE (czmls(1:ke_soil+1), STAT=ierr)
    ALLOCATE (czhls(0:ke_soil+1), STAT=ierr)
    ALLOCATE (msoilgrib(0:ke_soil+1), STAT=ierr)
                    ! (careful: the level k=1 will be coded with 1,
                    !           but is in the depth of 0.5 cm!)

    ! Check, how many soil levels have been specified
    invar = COUNT(czml_soil(:) /= -1.0_wp)
    IF (invar == 0) THEN
      ! no level specifications have been read
      IF (ke_soil == ke_soil_d) THEN
        ! use the default
        PRINT *,'  *** Default specifications of soil main levels are used    *** '
        czmls(1:ke_soil+1) = czml_soil_d(1:ke_soil+1)
      ELSE
        PRINT *,' ERROR  *** no specifications of soil levels,    *** '
        PRINT *,' ERROR  *** but using wrong default              *** ', &
                  ke_soil_d, ke_soil
        ierrstat = 1002
      ENDIF
    ELSE
      IF (ke_soil+1 == invar) THEN
        lzequiv = .TRUE.
        DO k = 1, ke_soil+1
          IF (czml_soil(k) /= czml_soil_d(k)) THEN
            lzequiv = .FALSE.
          ENDIF
        ENDDO
        IF (lzequiv) THEN
          PRINT *,'  *** Default specifications of soil main levels are used *** '
          czmls(1:ke_soil+1) = czml_soil_d(1:ke_soil+1)
        ELSE
          PRINT *,'  *** WARNING: Own specifications of soil main levels are used *** '
          PRINT *,'  ***          These have to correspond to the levels  of the  *** '
          PRINT *,'  ***          coarse grid model!                              *** '
          czmls(1:ke_soil+1) = czml_soil(1:ke_soil+1)
        ENDIF
      ELSE
        PRINT *,' ERROR  *** wrong number of specifications ',           &
                'for soil levels  *** ', ke_soil, invar
        ierrstat = 1002
      ENDIF
    ENDIF

    IF (ierrstat == 0) THEN
      ! compute grib coded values of the depth of main soil levels
      msoilgrib(0) = 0_iintegers
      DO i = 1, ke_soil+1
        msoilgrib(i) = NINT (100 * czmls(i)+1.0E-7_wp)
      ENDDO
    ENDIF

    ! determine depth of half levels out of czmls
    czhls(0) = 0.0_wp
    DO k = 1, ke_soil+1
      czhls(k) = 2.0_wp * czmls(k) - czhls(k-1)
    ENDDO

    czbot_w_so = MIN(czbot_w_so, czhls(ke_soil))
    ibot_w_so = ke_soil
    DO i=1,ke_soil+1
      IF (czhls(i) <= czbot_w_so) ibot_w_so=i
    ENDDO
  ENDIF

  IF (lconv) THEN
#ifdef CLM
    IF ( (itype_conv /= 0) .AND. (itype_conv /= 2) .AND. (itype_conv /= 3) ) THEN
      PRINT  *, ' ERROR  *** Wrong type of convection scheme: ', itype_conv
      PRINT  *, '        *** Must be 0, 2 or 3!'
      ierrstat = 1002
    ENDIF
    IF ( itype_conv == 2 .and. lconf_avg ) THEN
      lconf_avg = .FALSE.
      PRINT  *
      PRINT *, '*** WARNING: Namelist parameter lconf_avg has a different effect ***'
      PRINT *, '***          on the IFS scheme than with the Tiedtke scheme.     ***'
      PRINT *, '***          Only lconf_avg=.FALSE. has yet been tested.         ***'
      PRINT *, '***          Thus, lconf_avg has been set to .FALSE.             ***'
    ENDIF
#else
    IF ( (itype_conv /= 0) .AND. (itype_conv /= 3) ) THEN
      PRINT  *, ' ERROR  *** Wrong type of convection scheme: ', itype_conv
      PRINT  *, '        *** Must be 0 or 3!'
      ierrstat = 1002
    ENDIF
#endif
    IF (.NOT.ltur) THEN
       lctke=.FALSE.
    END IF
  ENDIF

  IF ( (itype_albedo  < 1) .OR. (itype_albedo  > 4) ) THEN
    PRINT  *, ' ERROR  *** Wrong type of albedo scheme: ', itype_albedo
    PRINT  *, '        *** Must be between 1 and 4!'
    ierrstat = 1002
  ENDIF

  IF ( (itype_albedo == 4) .AND. (.NOT. lforest) ) THEN
    PRINT  *, ' ERROR  *** albedo scheme 4 only possible with lforest=.TRUE.!'
    ierrstat = 1002
  ENDIF

  IF ( (itype_aerosol < 1) .OR. (itype_aerosol > 2) ) THEN
    PRINT  *, ' ERROR  *** Wrong type of aerosol scheme: ', itype_aerosol
    PRINT  *, '        *** Must be between 1 and 2!'
    ierrstat = 1002
  ENDIF

  IF ( (lperi_x .OR. lperi_y) .AND. itype_aerosol > 1 ) THEN
    PRINT  *, ' ERROR  *** Wrong type of aerosol scheme: ', itype_aerosol
    PRINT  *, '        *** For periodic BCs only itype_aerosol = 1 possible!'
    ierrstat = 1002
  ENDIF

  IF ( (itype_root    < 1) .OR. (itype_root    > 2) ) THEN
    PRINT  *, ' ERROR  *** Wrong type of root distribution: ', itype_root
    PRINT  *, '        *** Must be between 1 and 2!'
    ierrstat = 1002
  ENDIF

  IF ( (itype_heatcond < 1) .OR. (itype_heatcond > 2) ) THEN
    PRINT  *, ' ERROR  *** Wrong type of soil heat conductivity: ', itype_heatcond
    PRINT  *, '        *** Must be between 1 and 2!'
    ierrstat = 1002
  ENDIF

  IF ( (itype_hydbound /= 1) .AND. (itype_hydbound /= 3) ) THEN
    PRINT  *, ' ERROR  *** Wrong type of hydraulic lower boundary: ', itype_hydbound
    PRINT  *, '        *** Must be 1 or 3! (2 is not implemented yet)'
    ierrstat = 1002
  ENDIF

  IF (ltur .AND. l3dturb) THEN
    ! Add a warning, because field qvt_diff has been eliminated in Version 4.23!!
    PRINT *,'  *** WARNING: Explicit turbulent mixing tendencies from       *** '
    PRINT *,'  ***          horizontal mixing (l3dturb=TRUE) not accounted  *** '
    PRINT *,'  ***          for moisture divergence (dqvdt)!                *** '
  ENDIF

ENDIF

!------------------------------------------------------------------------------
!- Section 5: Distribute variables to all nodes
!------------------------------------------------------------------------------

IF (nproc > 1) THEN

  IF (my_world_id == 0) THEN
    intbuf ( 1) = nincrad
    intbuf ( 2) = ninctura
    intbuf ( 3) = nincconv
    intbuf ( 4) = nlgw

    intbuf ( 5) = itype_trvg
    intbuf ( 6) = itype_evsl

    intbuf ( 7) = itype_gscp

    intbuf ( 8) = itype_sher
    intbuf ( 9) = itype_wcld
    intbuf (10) = itype_tran
    intbuf (11) = itype_turb
    intbuf (12) = itype_synd

    intbuf (13) = imode_tran
    intbuf (14) = imode_turb

    intbuf (15) = icldm_rad
    intbuf (16) = icldm_tran
    intbuf (17) = icldm_turb
    intbuf (18) = ke_soil
    intbuf (19) = ico2_rad
    intbuf (20) = ibot_w_so
    intbuf (21) = nradcoarse
    intbuf (22) = nhori
    intbuf (23) = nextrad
    intbuf (24) = itype_conv
    intbuf (25) = nincsso
    intbuf (26) = ke_snow
    intbuf (27) = itype_aerosol
    intbuf (28) = itype_root
    intbuf (29) = itype_heatcond
    intbuf (30) = itype_hydbound
    intbuf (31) = iy_co2_stab
    intbuf (32) = itype_albedo
#ifdef TWOMOM_SB
    intbuf (33) = iradpar_cloud
#endif

    logbuf ( 1) = lrad
    logbuf ( 2) = ltur
    logbuf ( 3) = lconv
    logbuf ( 4) = lgsp
    logbuf ( 5) = lsoil
    logbuf ( 6) = ltkecon
    logbuf ( 7) = lexpcor
    logbuf ( 8) = ltmpcor
    logbuf ( 9) = lprfcor
    logbuf (10) = lnonloc
    logbuf (11) = lcpfluc
    logbuf (12) = lcape
    logbuf (13) = lctke
    logbuf (14) = lconf_avg
    logbuf (15) = lmelt
    logbuf (16) = lmelt_var
    logbuf (17) = lmulti_layer
    logbuf (18) = ldiniprec
    logbuf (19) = l3dturb
    logbuf (20) = lprog_tke
    logbuf (21) = limpltkediff
    logbuf (22) = lforest
    logbuf (23) = lconv_inst
    logbuf (24) = llake
    logbuf (25) = l3dturb_metr
    logbuf (26) = lradf_avg
    logbuf (27) = lradtopo
    logbuf (28) = lsso
    logbuf (29) = ltkesso
    logbuf (30) = lseaice
    logbuf (31) = lmulti_snow
    logbuf (32) = lemiss
    logbuf (33) = lstomata
    logbuf (34) = lco2_stab
    logbuf (35) = l_2mom
    logbuf (36) = lsuper_coolw

    realbuf( 1) = czbot_w_so
    realbuf( 2) = hincrad
    realbuf( 3) = hnextrad

    IF (lmulti_layer) THEN
      DO i = 1, ke_soil+1
        realbuf(3+i) = czmls(i)
      ENDDO
    ELSE
      realbuf( 4:40) = 0.0_wp
    ENDIF
    
  ENDIF

  CALL distribute_values (intbuf, 33, 0, imp_integers, icomm_world, ierr)
  CALL distribute_values (logbuf, 36, 0, imp_logical,  icomm_world, ierr)
  CALL distribute_values (realbuf,40, 0, imp_reals,    icomm_world, ierr)

  IF (my_world_id /= 0) THEN
    nincrad      = intbuf ( 1)
    ninctura     = intbuf ( 2)
    nincconv     = intbuf ( 3)
    nlgw         = intbuf ( 4)

    itype_trvg   = intbuf ( 5)
    itype_evsl   = intbuf ( 6)

    itype_gscp   = intbuf ( 7)

    itype_sher   = intbuf ( 8)
    itype_wcld   = intbuf ( 9)
    itype_tran   = intbuf (10)
    itype_turb   = intbuf (11)
    itype_synd   = intbuf (12)

    imode_tran   = intbuf (13)
    imode_turb   = intbuf (14)

    icldm_rad    = intbuf (15)
    icldm_tran   = intbuf (16)
    icldm_turb   = intbuf (17)
    ke_soil      = intbuf (18)
    ico2_rad     = intbuf (19)
    ibot_w_so    = intbuf (20)
    nradcoarse   = intbuf (21)
    nhori        = intbuf (22)
    nextrad      = intbuf (23)
    itype_conv   = intbuf (24)
    nincsso      = intbuf (25)
    ke_snow      = intbuf (26)
    itype_aerosol= intbuf (27)
    itype_root   = intbuf (28)
    itype_heatcond=intbuf (29)
    itype_hydbound=intbuf (30)
    iy_co2_stab  = intbuf (31)
    itype_albedo = intbuf (32)
#ifdef TWOMOM_SB
    iradpar_cloud= intbuf (33)
#endif

    lrad         = logbuf ( 1)
    ltur         = logbuf ( 2)
    lconv        = logbuf ( 3)
    lgsp         = logbuf ( 4)
    lsoil        = logbuf ( 5)

    ltkecon      = logbuf ( 6)
    lexpcor      = logbuf ( 7)
    ltmpcor      = logbuf ( 8)
    lprfcor      = logbuf ( 9)
    lnonloc      = logbuf (10)
    lcpfluc      = logbuf (11)

    lcape        = logbuf (12)
    lctke        = logbuf (13)
    lconf_avg    = logbuf (14)
    lmelt        = logbuf (15)
    lmelt_var    = logbuf (16)
    lmulti_layer = logbuf (17)
    ldiniprec    = logbuf (18)
    l3dturb      = logbuf (19)
    lprog_tke    = logbuf (20)
    limpltkediff = logbuf (21)
    lforest      = logbuf (22)
    lconv_inst   = logbuf (23)
    llake        = logbuf (24)
    l3dturb_metr = logbuf (25)
    lradf_avg    = logbuf (26)
    lradtopo     = logbuf (27)
    lsso         = logbuf (28)
    ltkesso      = logbuf (29)
    lseaice      = logbuf (30)
    lmulti_snow  = logbuf (31)
    lemiss       = logbuf (32)
    lstomata     = logbuf (33)
    lco2_stab    = logbuf (34)
    l_2mom       = logbuf (35)
    lsuper_coolw = logbuf (36)

    czbot_w_so   = realbuf( 1)
    hincrad      = realbuf( 2)
    hnextrad     = realbuf( 3)

    IF (lmulti_layer) THEN
      ALLOCATE (czmls(1:ke_soil+1), STAT=ierr)
      ALLOCATE (czhls(0:ke_soil+1), STAT=ierr)
      ALLOCATE (msoilgrib(0:ke_soil+1), STAT=ierr)

      ! determine depth of half levels out of czmls
      msoilgrib(0) = 0_iintegers
      czhls    (0) = 0.0_wp
      DO i = 1, ke_soil+1
        czmls(i) = realbuf(3+i)
        czhls(i) = 2.0_wp * czmls(i) - czhls(i-1)
        msoilgrib(i) = NINT (100 * czmls(i)+1.0E-7_wp)
      ENDDO
    ENDIF
  ENDIF

ENDIF

! Set flag l_dzeta_d_needed
IF (l3dturb_metr) THEN
  l_dzeta_d_needed = .TRUE.
END IF

! Set lprog_qi now internally, depending on itype_gscp
IF (itype_gscp < 3) THEN
  lprog_qi = .FALSE.
ELSE
  lprog_qi = .TRUE.
ENDIF

!------------------------------------------------------------------------------
!- Section 6: Output of the namelist variables and their default values
!------------------------------------------------------------------------------

IF (my_world_id == 0) THEN

  WRITE (nuspecif, '(A2)')  '  '
  WRITE (nuspecif, '(A23)') '0     NAMELIST:  phyctl'
  WRITE (nuspecif, '(A23)') '      -----------------'
  WRITE (nuspecif, '(A2)')  '  '
  WRITE (nuspecif, '(T7,A,T33,A,T52,A,T70,A)') 'Variable', 'Actual Value',   &
                                               'Default Value', 'Format'

  WRITE (nuspecif, '(T8,A,T33,L12  ,T52,L12  ,T71,A3)')                      &
                                           'lgsp'   ,lgsp   ,lgsp_d,   ' L '
  WRITE (nuspecif, '(T8,A,T33,L12  ,T52,L12  ,T71,A3)')                      &
                   'lsuper_coolw'   , lsuper_coolw  ,lsuper_coolw_d,   ' L '
  WRITE (nuspecif, '(T8,A,T33,L12  ,T52,L12  ,T71,A3)')                      &
                            'ldiniprec'   ,ldiniprec   ,ldiniprec_d,   ' L '
  WRITE (nuspecif, '(T8,A,T33,L12  ,T52,L12  ,T71,A3)')                      &
                                           'lrad',   lrad,   lrad_d,   ' L '
  WRITE (nuspecif, '(T8,A,T33,I12  ,T52,I12  ,T71,A3)')                      &
                       'itype_aerosol', itype_aerosol, itype_aerosol_d,' I '
  WRITE (nuspecif, '(T8,A,T33,L12  ,T52,L12  ,T71,A3)')                      &
                                     'lemiss   ',lemiss   ,lemiss_d   ,' L '
  WRITE (nuspecif, '(T8,A,T33,L12  ,T52,L12  ,T71,A3)')                      &
                                           'lforest',lforest,lforest_d,' L '
  WRITE (nuspecif, '(T8,A,T33,L12  ,T52,L12  ,T71,A3)')                      &
                                           'ltur   ',ltur   ,ltur_d   ,' L '
  WRITE (nuspecif, '(T8,A,T33,L12  ,T52,L12  ,T71,A3)')                      &
                                      'l3dturb', l3dturb, l3dturb_d,   ' L '
  WRITE (nuspecif, '(T8,A,T33,L12  ,T52,L12  ,T71,A3)')                      &
                       'l3dturb_metr', l3dturb_metr, l3dturb_metr_d,   ' L '
  WRITE (nuspecif, '(T8,A,T33,L12  ,T52,L12  ,T71,A3)')                      &
                                'lprog_tke', lprog_tke, lprog_tke_d,   ' L '
  WRITE (nuspecif, '(T8,A,T33,L12  ,T52,L12  ,T71,A3)')                      &
                       'limpltkediff', limpltkediff, limpltkediff_d,   ' L '
  WRITE (nuspecif, '(T8,A,T33,L12  ,T52,L12  ,T71,A3)')                      &
                                           'lconv  ',lconv  ,lconv_d  ,' L '
  WRITE (nuspecif, '(T8,A,T33,I12  ,T52,I12  ,T71,A3)')                      &
                                'itype_conv', itype_conv, itype_conv_d,' I '
  WRITE (nuspecif, '(T8,A,T33,L12  ,T52,L12  ,T71,A3)')                      &
                                  'lconv_inst',lconv_inst,lconv_inst_d,' L '
  WRITE (nuspecif, '(T8,A,T33,L12  ,T52,L12  ,T71,A3)')                      &
                                           'lsoil  ',lsoil  ,lsoil_d  ,' L '
  WRITE (nuspecif, '(T8,A,T33,L12  ,T52,L12  ,T71,A3)')                      &
                                           'lseaice',lseaice,lseaice_d,' L '
  WRITE (nuspecif, '(T8,A,T33,L12  ,T52,L12  ,T71,A3)')                      &
                                           'llake  ',llake  ,llake_d  ,' L '
  WRITE (nuspecif, '(T8,A,T33,L12  ,T52,L12  ,T71,A3)')                      &
                                           'lsso   ',lsso   ,lsso_d   ,' L '
  WRITE (nuspecif, '(T8,A,T33,L12  ,T52,L12  ,T71,A3)')                      &
                                           'lmelt  ',lmelt  ,lmelt_d  ,' L '
  WRITE (nuspecif, '(T8,A,T33,L12  ,T52,L12  ,T71,A3)')                      &
                                 'lmelt_var',lmelt_var  ,lmelt_var_d  ,' L '
  WRITE (nuspecif, '(T8,A,T33,L12  ,T52,L12  ,T71,A3)')                      &
                         'lmulti_layer ',lmulti_layer,lmulti_layer_d  ,' L '

  IF (lmulti_layer) THEN
    ! Write specification for soil levels
    WRITE (nuspecif, '(T8,A,T33,I12  ,T52,I12  ,T71,A3)')                    &
                                        'ke_soil ',ke_soil ,ke_soil_d, ' I '
    WRITE (nuspecif, '(T8,A)') 'Main levels of the soil layers:'
!!$    WRITE (nuspecif, '(T10,A)') '                  (m)           (cm)'
!!$    WRITE (nuspecif, '(A,I12)') '              0:               ', msoilgrib(0)
!!$    DO i = 1, ke_soil+1
!!$      WRITE (nuspecif, '(I15,A,F12.4,I12)') i, ':   ',czmls(i), msoilgrib(i)
!!$    ENDDO

    DO i = 1, MAX(ke_soil+1,ke_soil_d+1)
      WRITE (nuspecif, '(T8,A,I3.3,A,T33,F12.4,T52,F12.4,T71,A)')            &
           'czml_soil(',i,')',czmls(i),  czml_soil_d(i), ' R '
    ENDDO

    WRITE (nuspecif, '(T8,A,T33,F12.4,T52,F12.4,T71,A3)')                    &
                                  'czbot_w_so',czbot_w_so,czbot_w_so_d,' R '
  ENDIF

  WRITE (nuspecif, '(T8,A,T33,L12  ,T52,L12  ,T71,A3)')                      &
                         'lmulti_snow  ',lmulti_snow,  lmulti_snow_d  ,' L '
  WRITE (nuspecif, '(T8,A,T33,I12  ,T52,I12  ,T71,A3)')                      &
                                        'ke_snow ',ke_snow ,ke_snow_d, ' I '
  WRITE (nuspecif, '(T8,A,T33,I12  ,T52,I12  ,T71,A3)')                      &
              'itype_heatcond',itype_heatcond  ,itype_heatcond_d   ,   ' I '
  WRITE (nuspecif, '(T8,A,T33,I12  ,T52,I12  ,T71,A3)')                      &
              'itype_hydbound',itype_hydbound  ,itype_hydbound_d   ,   ' I '
  WRITE (nuspecif, '(T8,A,T33,I12  ,T52,I12  ,T71,A3)')                      &
                          'itype_root',itype_root  ,itype_root_d   ,   ' I '
  WRITE (nuspecif, '(T8,A,T33,L12  ,T52,L12  ,T71,A3)')                      &
                                'lstomata   ',lstomata   ,lstomata_d  ,' L '
  WRITE (nuspecif, '(T8,A,T33,L12  ,T52,L12  ,T71,A3)')                      &
                                           'ltkesso',ltkesso,ltkesso_d,' L '
  WRITE (nuspecif, '(T8,A,T33,L12  ,T52,L12  ,T71,A3)')                      &
                                           'ltkecon',ltkecon,ltkecon_d,' L '
  WRITE (nuspecif, '(T8,A,T33,L12  ,T52,L12  ,T71,A3)')                      &
                                           'lexpcor',lexpcor,lexpcor_d,' L '
  WRITE (nuspecif, '(T8,A,T33,L12  ,T52,L12  ,T71,A3)')                      &
                                           'ltmpcor',ltmpcor,ltmpcor_d,' L '
  WRITE (nuspecif, '(T8,A,T33,L12  ,T52,L12  ,T71,A3)')                      &
                                           'lprfcor',lprfcor,lprfcor_d,' L '
  WRITE (nuspecif, '(T8,A,T33,L12  ,T52,L12  ,T71,A3)')                      &
                                           'lnonloc',lnonloc,lnonloc_d,' L '
  WRITE (nuspecif, '(T8,A,T33,L12  ,T52,L12  ,T71,A3)')                      &
                                           'lcpfluc',lcpfluc,lcpfluc_d,' L '
  WRITE (nuspecif, '(T8,A,T33,L12  ,T52,L12  ,T71,A3)')                      &
                                   'lconf_avg',lconf_avg,lconf_avg_d  ,' L '
  WRITE (nuspecif, '(T8,A,T33,L12  ,T52,L12  ,T71,A3)')                      &
                                   'lradf_avg',lradf_avg,lradf_avg_d  ,' L '
  WRITE (nuspecif, '(T8,A,T33,L12  ,T52,L12  ,T71,A3)')                      &
                                           'lcape  ',lcape  ,lcape_d  ,' L '
  WRITE (nuspecif, '(T8,A,T33,L12  ,T52,L12  ,T71,A3)')                      &
                                           'lctke  ',lctke  ,lctke_d  ,' L '

  WRITE (nuspecif, '(T8,A,T33,F12.4,T52,F12.4,T71,A3)')                      &
                                           'hincrad',hincrad,hincrad_d,' R '
  WRITE (nuspecif, '(T8,A,T33,I12  ,T52,I12  ,T71,A3)')                      &
                                           'nincrad',nincrad,nincrad_d,' I '
  WRITE (nuspecif, '(T8,A,T33,I12  ,T52,I12  ,T71,A3)')                      &
                                           'nradcoarse',nradcoarse,nradcoarse_d,' I '
  WRITE (nuspecif, '(T8,A,T33,I12  ,T52,I12  ,T71,A3)')                      &
                                        'ninctura',ninctura,ninctura_d,' I '
  WRITE (nuspecif, '(T8,A,T33,I12  ,T52,I12  ,T71,A3)')                      &
                                        'nincconv',nincconv,nincconv_d,' I '
  WRITE (nuspecif, '(T8,A,T33,I12  ,T52,I12  ,T71,A3)')                      &
                                           'nincsso',nincsso,nincsso_d,' I '

  WRITE (nuspecif, '(T8,A,T33,I12  ,T52,I12  ,T71,A3)')                      &
                                    'itype_trvg',itype_trvg,itype_trvg_d,' I '
  WRITE (nuspecif, '(T8,A,T33,I12  ,T52,I12  ,T71,A3)')                      &
                                    'itype_evsl',itype_evsl,itype_evsl_d,' I '

  WRITE (nuspecif, '(T8,A,T33,I12  ,T52,I12  ,T71,A3)')                      &
                                    'itype_gscp',itype_gscp,itype_gscp_d,' I '

  WRITE (nuspecif, '(T8,A,T33,I12  ,T52,I12  ,T71,A3)')                      &
                                    'itype_sher',itype_sher,itype_sher_d,' I '
  WRITE (nuspecif, '(T8,A,T33,I12  ,T52,I12  ,T71,A3)')                      &
                                    'itype_wcld',itype_wcld,itype_wcld_d,' I '
  WRITE (nuspecif, '(T8,A,T33,I12  ,T52,I12  ,T71,A3)')                      &
                                    'itype_tran',itype_tran,itype_tran_d,' I '
  WRITE (nuspecif, '(T8,A,T33,I12  ,T52,I12  ,T71,A3)')                      &
                                    'itype_turb',itype_turb,itype_turb_d,' I '
  WRITE (nuspecif, '(T8,A,T33,I12  ,T52,I12  ,T71,A3)')                      &
                                    'itype_synd',itype_synd,itype_synd_d,' I '

  WRITE (nuspecif, '(T8,A,T33,I12  ,T52,I12  ,T71,A3)')                      &
                                    'imode_tran',imode_tran,imode_tran_d,' I '
  WRITE (nuspecif, '(T8,A,T33,I12  ,T52,I12  ,T71,A3)')                      &
                                    'imode_turb',imode_turb,imode_turb_d,' I '
  WRITE (nuspecif, '(T8,A,T33,I12  ,T52,I12  ,T71,A3)')                      &
                                    'ico2_rad',ico2_rad,ico2_rad_d,      ' I '
  WRITE (nuspecif, '(T8,A,T33,I12  ,T52,I12  ,T71,A3)')                      &
                                    'iy_co2_stab',iy_co2_stab,iy_co2_stab_d,' I '
  WRITE (nuspecif, '(T8,A,T33,L12  ,T52,L12  ,T71,A3)')                      &
                                    'lco2_stab',lco2_stab,lco2_stab_d,   ' L '
#ifdef TWOMOM_SB
  WRITE (nuspecif, '(T8,A,T33,I12  ,T52,I12  ,T71,A3)')                      &
                           'iradpar_cloud',iradpar_cloud,iradpar_cloud_d,' I '
#endif
  WRITE (nuspecif, '(T8,A,T33,I12  ,T52,I12  ,T71,A3)')                      &
                                    'icldm_rad',icldm_rad,icldm_rad_d,   ' I '
  WRITE (nuspecif, '(T8,A,T33,I12  ,T52,I12  ,T71,A3)')                      &
                                    'icldm_tran',icldm_tran,icldm_tran_d,' I '
  WRITE (nuspecif, '(T8,A,T33,I12  ,T52,I12  ,T71,A3)')                      &
                                    'icldm_turb',icldm_turb,icldm_turb_d,' I '

  WRITE (nuspecif, '(T8,A,T33,I12  ,T52,I12  ,T71,A3)')                      &
                                        'nlgw    ',nlgw    ,nlgw_d,    ' I '
  WRITE (nuspecif, '(T8,A,T33,L12  ,T52,L12  ,T71,A3)')                      &
                                      'lradtopo',lradtopo,lradtopo_d  ,' L '
  WRITE (nuspecif, '(T8,A,T33,I12  ,T52,I12  ,T71,A3)')                      &
                                        'nhori', nhori, nhori_d,       ' I '
  WRITE (nuspecif, '(T8,A,T33,I12  ,T52,I12  ,T71,A3)')                      &
                       'itype_albedo', itype_albedo, itype_albedo_d,' I '
  WRITE (nuspecif, '(A2)')  '  '

ENDIF

!------------------------------------------------------------------------------
!- End of the Subroutine
!------------------------------------------------------------------------------

END SUBROUTINE input_phyctl

!==============================================================================

!------------------------------------------------------------------------------
! End of module procedure organize_physics
!------------------------------------------------------------------------------

END SUBROUTINE organize_physics
