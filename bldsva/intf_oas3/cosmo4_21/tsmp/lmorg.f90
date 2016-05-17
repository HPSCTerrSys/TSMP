! $RCSfile: lmorg.f90,v $
! $Revision: 4.21 $ $Date: 2012/03/20 $
!+ Main program for the LM
!------------------------------------------------------------------------------

PROGRAM lmorg

!------------------------------------------------------------------------------
!
! Description:
!  Lmorg is the main program that organizes the forecast done by the LM.
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
!  Introduction of digital filtering for data initialization
! 1.4        1998/05/22 Guenther Doms
!  Adaptions for the two time-level time integration scheme
! 1.5        1998/06/29 Guenther Doms
!  Call of new routine 'print_vertcoord' for printing the vertical coordinates
! 1.8        1998/08/03 Ulrich Schaettler
!  Added parameterlist to init_meanvalues
! 1.9        1998/09/16 Guenther Doms
!  Call of routine 'vertical_diffusion' in case of two_time_level integration.
! 1.10       1998/09/29 Ulrich Schaettler
!  Change parameterlist of routine collect_timings.
! 1.12       1998/10/19 Ulrich Schaettler
!  Introduced possibility of boundary values from different runs.
! 1.13       1998/10/22 Christoph Schraff
!  Introduction of data assimilation by the nudging method.
! 1.17       1998/11/17 Ulrich Schaettler
!  Writing ready files after output steps, if required.
! 1.19       1998/12/11 Christoph Schraff
!  Bug fix for writing ready files for analyses.
! 1.20       1999/01/07 Guenhter Doms
!  Renaming of some global variables.
! 1.22       1999/02/08 Michael Buchhold
!  Introduction of 2D analyses of near surface parameters (e.g. t2m, rh2m)
! 1.24       1999/03/01 Guenther Doms
!  Call of the cloud ice parameterization schema (optionally).
! 1.29       1999/05/11 Ulrich Schaettler
!  Adapted interfaces to utility-modules and excluded several internal 
!  procedures
! 1.33       1999/10/14 Reinhold Hess
!  Call to budget diagnosis 'organize_diabudget', additional use of ldiag.
! 1.34       1999/12/10 Ulrich Schaettler
!  Introduced organization routines to most components which changes most
!  interfaces.
! 1.39       2000/05/03 Ulrich Schaettler
!  Adaptations to include asynchronous IO and the interactive nesting.
! 1.40       2000/05/23 Ulrich Schaettler
!  Correction for using MPE_IO on sequential platforms.
! 2.6        2001/06/12 Guenther Doms
!  Correction in an data exchange related to the cloud-ice scheme
! 2.8        2001/07/06 Ulrich Schaettler
!  Eliminated non-necessary variables from the USE-lists
!  Inquire the ASCII files now by UNIT-number instead of filename
! 2.11       2001/09/28 Ulrich Schaettler
!  Changed treatment of cloud ice, if initial and boundary fields are provided
! 2.17       2002/05/08 Ulrich Schaettler
!  Initializations of new time level for multi-layer soil model variables
! 2.18       2002/07/16 Ulrich Schaettler
!  Eliminated variable rhde, use cf_snow instead (by Reinhold Schrodin)
!  Exchange of variables qr, qs for 2 time level scheme (by Almut Gassmann)
!  Changed initialization of new time level, if working with boundary data
!  defined on frames (by Lucio Torrisi).
!  Use lmulti_layer to define boundary values for old soil model correctly
! 2.19       2002/10/24 Ulrich Schaettler
!  Moved re-initialization of vmax_10m from near_surface to initialize_loop.
! 3.2        2003/02/07 Ulrich Schaettler
!  Moved some communications from the Convection scheme to lmorg.
! 3.5        2003/09/02 Ulrich Schaettler
!  Adaptations for the interface of exchg_boundaries
!  Corrections for the digital filtering
! 3.6        2003/12/11 Ulrich Schaettler
!  Changed interface for calling the digital filtering dfi_initialization
! 3.7        2004/02/18 Ulrich Schaettler
!  Introduced CALLs to new component for computing synthetic satellite images.
!  Adaptations in the boundary exchange for prognostic precipitation.
!  Treatment of humidity-variables (are set to 0.0, if they are too small)
! 3.8        2004/03/23 Ulrich Schaettler
!  Correction in initialize_loop for the 2 timelevel scheme.
! 3.11       2004/07/21 Michael Gertz
!  FPE trapping for IBM machines introduced
! 3.13       2004/12/03 Ulrich Schaettler
!  Introduced CALL for Latent Heat Nudging (LHN) within organize_assimilation
!                                                      (Klaus Stephan)
!  Adaptations in boundary exchange for new graupel scheme; 
!  Modification for 2D runs (l2dim) (Thorsten Reinhardt)
! 3.14       2005/01/25 Ulrich Schaettler
!  Adapted treatment of t_s for new multi-layer soil model 
! 3.15       2005/03/03 Ulrich Schaettler
!  Replaced FLOAT by REAL
! 3.16       2005/07/22 Ulrich Schaettler
!  Adaptations for the 2 time level Runge-Kutta scheme
! 3.18       2006/03/03 Ulrich Schaettler
!  Changed treatment of ASCII files to introduce restart possibility
!  Additional boundary variables in case of lbdclim=.TRUE.
!  Added get_timing routine for LHN   (Klaus Stephan)
!  Changed computation of surface specific humidity for lakes (D. Mironov)
!  Added treatment of vertical velocity w for lw_freeslip=.FALSE.
! 3.21       2006/12/04 Ulrich Schaettler, Jochen Foerstner
!  Read Namelist Input for physics in any case (also for "dry" runs)
!  Read Namelist Input for EPS mode, if leps=.TRUE. 
!  Call near_surface with time level nnow
!  Set new time level for qr, qs, qg 
! V3_23        2007/03/30 Ulrich Schaettler
!  Adapting print outs to idbg_level and flexible dt
! V3_24        2007/04/26 Ulrich Schaettler
!  Bug correction in call to dfi_initialization
!  Eliminated nincmxu and introduced control as for other increments
! V4_3         2008/02/25 Ulrich Schaettler 
!  Correction for computing zforecasttime in case of restart runs
!  Adapted call to init_timings
! V4_4         2008/07/16 Ulrich Schaettler, Dmitrii Mironov
!  Use qitens so that cloud ice can be treated in the same way as cloud water (DM)
!  Eliminated timing variables which are unused
!  Adapted interface of get_timings
!  Changed NL parameter lyear_360 to itype_calendar, to have several options
! V4_8         2009/02/16 Ulrich Schaettler
!  Renamed switch for artificial data from lgen to lartif_data
!  Reset also variables for dynamical and convective gusts
! V4_9         2009/07/16 Ulrich Schaettler, Heike Vogel
!  Implement a first version of COSMO_ART
!  Introduced internal subroutine set_qrqsqg_boundaries. Set these boundaries
!   for the total field to 0, if no boundary fields are given.
!  Replace hardcoded 0 in call to get_timings for cleanup by nstop
!  Read diagnosis namelist even if ldiagnos=.FALSE.; (some settings are needed)
!  Inserted NEC Compiler directives
!  Eliminate reduction of qv in case lana_qi, llb_qi = .TRUE.
! V4_11        2009/11/30 Ekaterina Machulskaya
!  Adaptations to run with the multi-layer snow model
! V4_12        2010/05/11 Jan-Peter Schulz
!  Adapted interface parameter t_s for SR tgcom to account for seaice (JPS)
! V4_13        2010/05/11 Michael Gertz
!  Adaptions to SVN
!  New information module info_lm_f90 added
! V4_14        2010/06/14 Ulrich Schaettler
!  Putting the call for the info-module to the right place
! V4_17        2011/02/24 Ulrich Blahak
!  Adapted interface of exchg_boundaries; corrected kzdims(1:20) -> kzdims(1:24);
!  eliminated my_peri_neigh; added calls for artificial disturbances on
!  initial temperature (atmosphere and/or soil) or artificial heating rates
!  (atmosphere and/or soil) or artificial t_s disturbances in case of lsoil=.false.;
! V4_18        2011/05/26 Ulrich Schaettler
!  Introduced conditional compilation for Nudging and synthetic satellite images
!  from Christoph Knote (for COSMO-ART)
!  lphot removed, not necessary anymore
!  lho, lho2 used instead of hardcoded cgas species index
!  ART boundaries use and calculate their own relaxation weights
! V4_20        2011/08/31 Ulrich Schaettler
!  Implemented interface to OASIS coupler using conditional compilation with -DCOUP_OAS
!   (by CLM Community)
! V4_21        2011/12/06 Oliver Fuhrer
!  Bug Fix in call to exchg_boundaries for POLLEN runs
! V4_21        2012/03/20 Alexander Kelbch
!  Inclusion of tracer (adopted from Markus Uebel)
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!==============================================================================
! Declarations:
! Modules used:

! Load the library information data:
USE info_lm_f90,         ONLY: info_define, info_readnl, info_print

USE data_parameters,     ONLY:   ireals, iintegers
USE data_constants,      ONLY:   b1, b2w, b3, b4w, b2i, b4i, rdv, o_m_rdv,   &
                                 rvd_m_o, r_d
USE data_soil,           ONLY:   cf_snow
USE data_fields,         ONLY:                                               &
       dp0, p0, rho, rho0, qrs, llandmask, t_g, vmax_10m, dqvdt, qvsflx,     &
       u, v, w, t, qv, qc, qi, qr, qs, qg, tke, pp, ps, u_bd, v_bd, w_bd,    &
       t_bd, qv_bd, qc_bd, qi_bd, qr_bd, qs_bd, qg_bd, pp_bd,                &
       t_snow,    t_s,    qv_s,    t_m,    w_snow,    w_g1,    w_g2,         &
       t_snow_bd, t_s_bd, qv_s_bd, t_m_bd, w_snow_bd, w_g1_bd, w_g2_bd,      &
       plcov_bd, lai_bd, rootdp_bd, vio3_bd, hmo3_bd, t_cl_bd, w_cl_bd,      &
       plcov,    lai,    rootdp,    vio3,    hmo3,    t_cl,    w_cl,         &
       utens, vtens, wtens, ttens, qvtens, qctens, pptens, w_g3, w_g3_bd,    &
       fr_lake, depth_lk, h_ice, qvt_diff, qitens, vgust_con, vgust_dyn,     &
       prr_gsp, prr_con, prne_con, bas_con, t_snow_mult, t_ice

!AK (20.03.12)
USE data_tracer,         ONLY:   tracer, ntracer, ltracer, tracer_bd,        &
                                 tracertens
!AK (20.03.12)

USE data_modelconfig,    ONLY:   ie, je, ke, ke1, jstartpar, jendpar,        &
                                 istartpar, iendpar, ivctype, vcoord,        &
                                 dt, dt2, ed2dt, dtdeh, nehddt, hhlr, sigmr, &
                                 vcflat, kflat, p0sl, t0sl, dt0lp,           &
                                 jstartu, jendu, jstartv, jendv,             &
                                 istart, iend, jstart, jend, ie_tot, je_tot

USE data_runcontrol,     ONLY:                                               &
       nstart, nstop, ntstep, nlastbound, nincbound, lmulti_layer,           &
       ltime, itype_timing, newbcdt, newbc, nincboufac, nlgw, itype_gscp,    &
       nbd1, nbd2, nold, nnow, nnew, lartif_data, ldiagnos, luse_rttov,      &
       lphys, lconv, luseobs, l2tls, lsemi_imp, lgsp, lprog_qi, lprogprec,   &
       yakdat1, yakdat2, nuspecif, ldfi, nincconv, lconf_avg, l2dim,         &
       irunge_kutta, nbl_exchg, lprog_tke, lspecnudge, lseaice, llake,       &
       lw_freeslip, leps, idbg_level, lprintdeb_all, itype_calendar,         &
       hlastmxu, hnextmxu, hincmxu, nlastmxu, nnextmxu, l_cosmo_art,         &
       l_pollen, hstart, itype_lbc_qrsg, lmulti_snow, lperi_x, lperi_y, l2dim

USE data_parallel,       ONLY:                                               &
       num_compute, icomm_cart, my_cart_id, sendbuf, isendbuflen, iexch_req, &
       imp_reals, nboundlines, my_cart_neigh, nprocx, nprocy, my_world_id,   &
       lcompute_pe, lasync_io, ncomm_type, ldatatypes, ltime_barrier

USE data_io,             ONLY:   ydate_ini, pp_nl, root, lbdclim,            &
                                 lana_qi, llb_qi, llb_qr_qs, llb_qg,         &
                                 lbd_frame, undef

USE data_flake,          ONLY:   h_Ice_min_flk

USE mpe_io,              ONLY:   mpe_io_node, mpe_io_shutdown

USE environment,         ONLY:   exchg_boundaries, comm_barrier,             &
                                 final_environment, model_abort, get_free_unit
USE meteo_utilities,     ONLY:   calrho, calps, tgcom
USE time_utilities,      ONLY:   nutiming, init_timings, get_timings,        &
                  i_initializations, i_add_computations, i_phy_computations, &
                  i_dyn_computations,i_lhn_computations, i_spectr_nudging,   &
                  i_relaxation, i_input, i_output, i_barrier_waiting_dyn,    &
                  i_communications_dyn, i_cleanup, collect_timings
USE utilities,           ONLY:   get_utc_date

USE src_setup,           ONLY:   organize_setup, constant_fields

USE src_allocation,      ONLY:   organize_allocation

USE src_artifdata,       ONLY:   artif_heatrate_dist, artif_heatrate_dist_tso, set_tempdist, &
                                 set_tempdist_tso, set_tempdist_bbc_ts, add_noise_tw

#ifdef COSMOART
USE data_cosmo_art,      ONLY:   lgas, laero, ldust, lseas,                  &
                                 lrad_dust, lrad_aero,                       &
                                 lgasbd, laerobd, lgasini, laeroini,         &
                                 cgas, caero, cgas_bd, caero_bd, cgastens,   &
                                 caerotens, isp_gas, isp_aero,               &
                                 nlastbound_art, nincbound_art,              &
                                 nbd1_art, nbd2_art
USE art_utilities,       ONLY:   gm3ppm_a, ppmgm3_a
USE art_gas_const,       ONLY:   lho, lho2
USE art_aerosol_const,   ONLY:   aerostart
#endif

#ifdef POLLEN
USE data_pollen,         ONLY:   cpollen, cpollentens, cpollen_bd,           &
                                 isp_pollen, dtpollen
#endif

!AK (20.03.12)
USE src_tracer_supply,   ONLY:   organize_tracer_init, organize_tracer,       &
                                 organize_tracer_bound, organize_tracer_source
!AK (20.03.12) 

!==============================================================================

IMPLICIT NONE

!==============================================================================
!
! Declarations must be of the form:
! <type> <VariableName> ! Description/ purpose of variable
!
! Local Scalars:

INTEGER  (KIND=iintegers) ::                          &
  ierrstat,      & ! error status variable
  izerror,       & !
  nzhours,       & ! for recording a forecast hour
  nzdays,        & ! for recording a forecast day
  i,j,k, kzdims(24), izdebug, nzdiv, isp, kstart,     &
!AK (20.03.2012)
  iig,           & ! loop index for ltracer
  iprog,         & ! loop index for prognostic processes of tracers
  nprog            ! number of prognostic processes of tracer
!AK (20.03.2012)

REAL (KIND=ireals)        ::                          &
  zgrids_dt(1),  & ! structure for interface to organize_data
                   ! (would be necessary only for the nesting version)
  zforecasttime, & ! for recording a forecast hour
  zrealdiff        ! for time measurement

LOGICAL                   ::                          &
  lopen,         & ! check whether YUCHKDAT is opened
  lzconv           ! to determine whether extra communication for convection
                   ! is needed

CHARACTER (LEN=46)        :: ynote
CHARACTER (LEN=80)        :: yzerrmsg

!==============================================================================

!------------------------------------------------------------------------------
!- End of Header
!------------------------------------------------------------------------------
#ifdef FPEABORT
! Floating point exception trapping
      CALL initialize_fpe_trap(.TRUE., ierrstat)
      IF ( ierrstat /= 0 ) THEN
        STOP 'Error initializing ftp_trap'
      END IF
#endif

  ierrstat = 0_iintegers
  izerror  = 0_iintegers

!------------------------------------------------------------------------------
!- Section 1: Setup of the model and Namelist Input for all components
!------------------------------------------------------------------------------

  ! Section 1 has to be computed by all PEs (Compute and IO-PEs)
  ! After Section 1 the PEs for computing and IO are splitted

  CALL organize_setup

  IF (my_cart_id == 0) THEN
    ! Print the default information to stdout:
    CALL info_define ('lmparbin')          ! Pre-define the program name
    CALL info_readnl ('INPUT_COSMO')       ! Read additional information from namelist file
    CALL info_print ()                     ! Print the information to stdout
  ENDIF

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

  ! Input of the dynamics namelist
  CALL organize_dynamics ('input', izerror, yzerrmsg, dt, .FALSE.)
  IF (izerror /= 0_iintegers) THEN
    CALL model_abort (my_world_id, 100+izerror, yzerrmsg,               &
                                   'organize_dynamics: input')
  ENDIF

  ! Input of the physics namelist
  CALL organize_physics ('input', izerror, yzerrmsg)
  IF (izerror /= 0_iintegers) THEN
    CALL model_abort (my_world_id, 100+izerror, yzerrmsg,             &
                                   'organize_physics: input')
  ENDIF

#ifdef COSMOART
  ! Input of the COSMO_ART namelist
  IF (l_cosmo_art) THEN
    CALL organize_cosmo_art ('input', ydate_ini, izerror, yzerrmsg)
    IF (izerror /= 0_iintegers) THEN
      CALL model_abort (my_world_id, 100+izerror, yzerrmsg,             &
                                   'organize_cosmo_art: input')
    ENDIF
  ENDIF
#endif

#ifdef POLLEN
  ! Input of the Pollen namelist
  IF (l_pollen) THEN
    CALL organize_pollen ('input', ydate_ini, izerror, yzerrmsg)
    IF (izerror /= 0_iintegers) THEN
      CALL model_abort (my_world_id, 100+izerror, yzerrmsg,             &
                                   'organize_pollen: input')
    ENDIF
  ENDIF
#endif

  ! Input of the diagnosis namelist
  CALL organize_diagnosis ('input', izerror, yzerrmsg)
  IF (izerror /= 0_iintegers) THEN
    CALL model_abort (my_world_id, 100+izerror, yzerrmsg,             &
                                   'organize_diagnosis: input')
  ENDIF

#ifdef NUDGING
  ! Input of the assimilation namelist
  IF (luseobs) THEN
    CALL organize_assimilation ('input', izerror, yzerrmsg)
    IF (izerror /= 0_iintegers) THEN
      CALL model_abort (my_world_id, 100+izerror, yzerrmsg,             &
                                     'organize_assimilation: input')
    ENDIF
  ENDIF
#endif

  ! Input of the EPS namelist
  IF (leps) THEN
    CALL organize_eps ('input', izerror, yzerrmsg)
    IF (izerror /= 0_iintegers) THEN
      CALL model_abort (my_world_id, 100+izerror, yzerrmsg,             &
                                      'organize_eps: input')
    ENDIF
  ENDIF

  ! Input of the namelists for the I/O-package
  zgrids_dt(1) = dt
  CALL organize_data ('input-init', 0, 1, 1, zgrids_dt,    &
                       izerror, yzerrmsg)
  IF (izerror /= 0_iintegers) THEN
    CALL model_abort (my_world_id, 100+izerror, yzerrmsg,               &
                                   'organize_data: input-init')
  ENDIF

#if defined RTTOV7 || defined RTTOV9 || defined RTTOV10
  IF (luse_rttov) THEN
    ! Input of the namelists for the I/O-package
    CALL organize_satellites ('input', izerror, yzerrmsg)
    IF (izerror /= 0_iintegers) THEN
      CALL model_abort (my_world_id, 100+izerror, yzerrmsg,             &
                                     'organize_satellites: input')
    ENDIF
  ENDIF
#endif

!AK (20.03.12)
  ! Input of the namelists for tracer transport
  CALL organize_tracer(izerror, yzerrmsg)
    IF (izerror /= 0_iintegers) THEN
      CALL model_abort (my_cart_id, 100+izerror, yzerrmsg,              &
                                     'organize_tracer')
    ENDIF
!AK (20.03.12)

  ! Initialize the timings
  CALL get_free_unit (nutiming)
  CALL init_timings (nstart, nstop, dt, itype_timing, ldfi,             &
                     lphys, luseobs, l2tls, lsemi_imp, l_cosmo_art, izerror)
  IF (izerror /= 0) THEN
    ! no system clock present
    ltime = .FALSE.
  ENDIF

!------------------------------------------------------------------------------
! Section 2: Allocation of space and computation of constant fields
!------------------------------------------------------------------------------

! Now comes the part for the compute PEs. The part for the IO-PEs is in
! the ELSE-part at the end of the program.

comp_pe: IF (lcompute_pe) THEN

  ! allocate space
  IF (izdebug > 0) THEN
    PRINT *,'    ALLOCATE SPACE'
  ENDIF

  ! fields for the meteorological variables
  CALL organize_allocation ('default', izerror)

  IF (izerror /= 0) THEN
    ierrstat = 1004
    yzerrmsg  = ' ERROR    *** Allocation of space for meteofields failed ***'
    CALL model_abort (my_cart_id, ierrstat, yzerrmsg, 'allocation: default')
  ENDIF
  
#ifdef COSMOART
  IF (l_cosmo_art) THEN
    CALL organize_cosmo_art ('allocate', ydate_ini, izerror, yzerrmsg)

    IF (izerror /= 0) THEN
      ierrstat = 1005
      yzerrmsg  = ' ERROR    *** Allocation of space for COSMO_ART failed ***'
      CALL model_abort (my_cart_id, ierrstat, yzerrmsg,             &
                                    'organize_cosmo_art: allocate')
    ENDIF
  ENDIF
#endif

#ifdef POLLEN
  IF (l_pollen) THEN
    CALL organize_pollen ('allocate', ydate_ini, izerror, yzerrmsg)

    IF (izerror /= 0) THEN
      ierrstat = 1005
      yzerrmsg  = ' ERROR    *** Allocation of space for Pollen failed ***'
      CALL model_abort (my_cart_id, ierrstat, yzerrmsg,             &
                                    'organize_pollen: allocate')
    ENDIF
  ENDIF
#endif

  CALL constant_fields

!------------------------------------------------------------------------------
!- Section 3: Input of first data sets
!------------------------------------------------------------------------------

  ! Read or generate initial data and the first boundary data sets
  zgrids_dt(1) = dt
  CALL organize_data ('start', 0, 1, 1, zgrids_dt, izerror, yzerrmsg)
  IF (izerror /= 0_iintegers) THEN
    CALL model_abort (my_cart_id, 100+izerror, yzerrmsg,                &
                                   'start: input-init')
  ENDIF

!------------------------------------------------------------------------------
!- Section 4: Initializations and allocation of extra space
!------------------------------------------------------------------------------

  IF (izdebug > 0) THEN
    PRINT *, '  INITIALIZATIONS'
  ENDIF

  ! 3.1:  Initialization of different packages
  CALL organize_dynamics ('init', izerror, yzerrmsg, dt, .FALSE.)
  IF (izerror /= 0_iintegers) THEN
    CALL model_abort (my_cart_id, 100+izerror, yzerrmsg,                &
                                   'organize_dynamics: init')
  ENDIF

!MU (17.10.12)
  ! Initialization of tracer
  CALL organize_tracer_init(izerror,yzerrmsg)
  IF (izerror /= 0_iintegers) THEN
      CALL model_abort (my_cart_id, 100+izerror, yzerrmsg,              &
                                     'organize_tracer_init: init')
  ENDIF
!MU (17.10.12)

  ! Initialization of the physics
  IF (lphys) THEN
    CALL organize_physics ('init', izerror, yzerrmsg)
    IF (izerror /= 0_iintegers) THEN
      CALL model_abort (my_cart_id, 100+izerror, yzerrmsg,              &
                                     'organize_physics: init')
    ENDIF

    ! allocate fields for meteorological variables concerned to the canopy
    CALL organize_allocation ('canopy', izerror)
    IF (izerror /= 0) THEN
      ierrstat = 1004
      yzerrmsg  = ' ERROR    *** Allocation of extra space failed ***'
      CALL model_abort (my_cart_id, ierrstat, yzerrmsg, 'allocation: canopy')
    ENDIF
  ENDIF
  
#ifdef COSMOART
  ! Initialization of COSMO_ART
  IF (l_cosmo_art) THEN
    CALL organize_cosmo_art ('init', ydate_ini, izerror, yzerrmsg)
    IF (izerror /= 0_iintegers) THEN
      CALL model_abort (my_cart_id, 100+izerror, yzerrmsg,              &
                                    'organize_cosmo_art: init')
    ENDIF

    ! Initial profiles of gas phase species of LM_ART
    IF (lgas) THEN
      CALL organize_cosmo_art ('start_gas', ydate_ini, izerror, yzerrmsg)
      IF (izerror /= 0_iintegers) THEN
        CALL model_abort (my_cart_id, 100+izerror, yzerrmsg,              &
                                    'organize_cosmo_art: start_gas')
      ENDIF
    ENDIF
  ENDIF
#endif

#ifdef POLLEN
  ! Initialization of the pollen
  IF (l_pollen) THEN
    CALL organize_pollen ('init', ydate_ini, izerror, yzerrmsg)
    IF (izerror /= 0_iintegers) THEN
      CALL model_abort (my_cart_id, 100+izerror, yzerrmsg,              &
                                     'organize_pollen: init')
    ENDIF
  ENDIF
#endif

  ! Initialization of the diagnosis
  IF (ldiagnos) THEN
    CALL organize_diagnosis ('init', izerror, yzerrmsg)
    IF (izerror /= 0_iintegers) THEN
      CALL model_abort (my_cart_id, 100+izerror, yzerrmsg,              &
                                     'organize_diagnosis: init')
    ENDIF
  ENDIF

#ifdef NUDGING
  ! Initialization of the assimilation
  IF (luseobs) THEN
    CALL organize_assimilation ('init', izerror, yzerrmsg)
    IF (izerror /= 0_iintegers) THEN
      CALL model_abort (my_cart_id, 100+izerror, yzerrmsg,              &
                                     'organize_assimilation: init')
    ENDIF
  ENDIF
#endif

#if defined RTTOV7 || defined RTTOV9 || defined RTTOV10
  IF (luse_rttov) THEN
    ! initialization of variables for the synthetic satellite computations
    CALL organize_satellites ('init', izerror, yzerrmsg)
    IF (izerror /= 0_iintegers) THEN
      CALL model_abort (my_cart_id, 100+izerror, yzerrmsg,              &
                                     'organize_satellites: init')
    ENDIF
  ENDIF
#endif

!------------------------------------------------------------------------------
! Section 5: Model initialization by digital filtering
!------------------------------------------------------------------------------

  IF (ldfi) THEN
    IF (izdebug > 0) THEN
      PRINT *, '    DIGITAL FILTERING'
    ENDIF
    IF(l2tls) THEN
      IF (izdebug > 0) THEN
        PRINT *, ' **** CAUTION **** CAUTION **** CAUTION **** CAUTION ****'
        PRINT *, ' ****   DIGITAL FILTERING not tested for 2TL-scheme  ****'
        PRINT *, ' **** **** ****  Proceed on your own risk  **** **** ****'
      ENDIF
    ENDIF
    CALL dfi_initialization (lana_qi, llb_qi, llb_qr_qs, llb_qg,      &
                             lbd_frame, undef, izerror, yzerrmsg)
    IF (izerror /= 0_iintegers) THEN
      CALL model_abort (my_cart_id, 100+izerror, yzerrmsg,            &
                                     'dfi_initialization')
    ENDIF
  ENDIF
  
  ! Close file for control output
  IF (my_cart_id == 0) THEN
    CLOSE (nuspecif, STATUS='KEEP')
  ENDIF

  ! halve the time step, if ntstep = 0
  ! (in src_setup: nstep = nstart)
  IF (.NOT. l2tls) THEN
    IF (ntstep == 0) THEN
      zforecasttime = 0.0
      dt = 0.5 * dt
    ELSE
      zforecasttime = 0.5*dt + (ntstep-1)*dt
      nzdiv         = INT (zforecasttime / 3600.0_ireals, iintegers)
      zforecasttime = zforecasttime - nzdiv * 3600_ireals
    ENDIF
  ELSE
    IF (ntstep == 0) THEN
      zforecasttime = 0.0
    ELSE
      zforecasttime = ntstep*dt
      nzdiv         = INT (zforecasttime / 3600.0_ireals, iintegers)
      zforecasttime = zforecasttime - nzdiv * 3600_ireals
    ENDIF
  ENDIF

  IF (lbdclim) THEN
    ynote         = '...... FORECAST TIME IS NOW xxxxxx DAYS ......'
  ELSE
    ynote         = '...... FORECAST TIME IS NOW xxx HOURS   ......'
  ENDIF

  IF (ltime) THEN
    CALL get_timings (i_initializations, ntstep, dt, izerror)
  ENDIF

  !------------------------------------------------------------------------------
  ! 5a. Temperature disturbance(s) (either in the air or at the soil surface
  !    or within the soil) in the initial conditions
  !------------------------------------------------------------------------------

  ! There are different types of possible disturbances, see the documentation
  ! of the corresponding namelist parameters in INPUT_IDEAL, and there is the possibility
  ! to specify more than one disturbance (up to 50 right now).
  
  IF (lartif_data) THEN
    CALL set_tempdist(nnew)
    ! Initial condition on t_so (takes only effect if lsoil=.true.)
    CALL set_tempdist_tso(nnew)
  END IF

#if defined COUP_OAS_COS
  CALL oas_cos_define
!OASIS4 only
!  CALL oas_cos_update_time(0)
#endif

!------------------------------------------------------------------------------
!- Section 6: Time stepping
!------------------------------------------------------------------------------

  IF (izdebug > 0) THEN
    PRINT *, '  TIME STEPPING'
  ENDIF

  timeloop: DO ntstep = nstart , nstop

    IF ( (izdebug > 1) .AND. (.NOT. lbdclim)) THEN
      PRINT *, '    STEP ',ntstep
    ENDIF

    !--------------------------------------------------------------------------
    !- Section 6.1: Initialization of this time step
    !--------------------------------------------------------------------------

!AK (20.03.12)
    tracertens(:,:,:,:) = 0._ireals
    CALL organize_tracer_bound(izerror,yzerrmsg)
    IF (izerror /= 0_iintegers) THEN
      CALL model_abort (my_cart_id, 100+izerror, yzerrmsg,               &
                                     'organize_tracer_bound: compute')
    ENDIF

    CALL organize_tracer_source(izerror,yzerrmsg)
    IF (izerror /= 0_iintegers) THEN
      CALL model_abort (my_cart_id, 100+izerror, yzerrmsg,               &
                                       'organize_tracer_source: compute')
    ENDIF
!AK (20.03.12)

    CALL initialize_loop (ntstep, nbd1, nbd2, nold, nnow, nnew)
 
    IF (ltime) CALL get_timings (i_add_computations, ntstep, dt, izerror)

    !--------------------------------------------------------------------------
    !- Section 6.2.1: physics 
    !--------------------------------------------------------------------------

#ifdef COSMOART
    ! CK 20101204 unit conversion necessary before physics
    ! more universal approach: get a general injection point for ART
    IF (l_cosmo_art) THEN
      CALL organize_cosmo_art ('prepare_physics', ydate_ini, izerror, yzerrmsg)
      IF (izerror /= 0_iintegers) THEN
        CALL model_abort (my_cart_id, 100+izerror, yzerrmsg,               &
                                       'organize_cosmo_art: prepare_physics')
      ENDIF
    ENDIF
#endif

#ifdef POLLEN
    ! Preparations for Pollen
    IF (l_pollen) THEN
      CALL organize_pollen ('prepare_physics', ydate_ini, izerror, yzerrmsg)
      IF (izerror /= 0_iintegers) THEN
        CALL model_abort (my_cart_id, 100+izerror, yzerrmsg,               &
                                       'organize_pollen: prepare_physics')
      ENDIF
    ENDIF
#endif

    IF (lphys) CALL organize_physics ('compute', izerror, yzerrmsg)
    IF (izerror /= 0_iintegers) THEN
      CALL model_abort (my_cart_id, 100+izerror, yzerrmsg,               &
                                     'organize_physics: compute')
    ENDIF

    IF (ltime) CALL get_timings (i_phy_computations, ntstep, dt, izerror)

    !--------------------------------------------------------------------------
    !- Section 6.2.1a: After the call to the soil model, add extra
    !  artificial heating rate disturbances on the soil temperature if desired
    !  (formulated as Euler forward updates of the temperatures)
    !--------------------------------------------------------------------------

    IF (lartif_data) THEN
      ! Set possible artificial heating rate disturbance(s) in the soil 
      ! (affects t_so or t_s/t_m/t_cl depending on soil model
      ! and takes effect only IF lsoil=.TRUE.). 
      ! Because the soil model has already done the time integration,
      ! the artificial disturbances have to be imposed on
      ! timelevel nnew:
      CALL artif_heatrate_dist_tso(nnew)
    END IF

#ifdef COSMOART
    ! CK 20101204 unit conversion necessary after physics
    ! more universal approach: get a general injection point for ART
    IF (l_cosmo_art) THEN
      CALL organize_cosmo_art ('finalize_physics', ydate_ini, izerror, yzerrmsg)
      IF (izerror /= 0_iintegers) THEN
        CALL model_abort (my_cart_id, 100+izerror, yzerrmsg,               &
                                       'organize_cosmo_art: finalize_physics')
      ENDIF
    ENDIF
#endif

#ifdef POLLEN
    ! Clean up after Pollen
    IF (l_pollen) THEN
      CALL organize_pollen ('finalize_physics', ydate_ini, izerror, yzerrmsg)
      IF (izerror /= 0_iintegers) THEN
        CALL model_abort (my_cart_id, 100+izerror, yzerrmsg,               &
                                       'organize_pollen: finalize_physics')
      ENDIF
    ENDIF
#endif

#ifdef COSMOART
    !--------------------------------------------------------------------------
    !- Section 6.2.2: emissions for COSMO_ART
    !--------------------------------------------------------------------------

    IF (l_cosmo_art) THEN
      CALL organize_cosmo_art ('emiss', ydate_ini, izerror, yzerrmsg)
      IF (izerror /= 0_iintegers) THEN
        CALL model_abort (my_cart_id, 100+izerror, yzerrmsg,               &
                                       'organize_cosmo_art: emiss')
      ENDIF
    ENDIF
#endif

#ifdef POLLEN
    IF (l_pollen) THEN
      CALL organize_pollen ('emiss', ydate_ini, izerror, yzerrmsg)
      IF (izerror /= 0_iintegers) THEN
        CALL model_abort (my_cart_id, 100+izerror, yzerrmsg,               &
                                       'organize_pollen: emiss')
      ENDIF
    ENDIF
#endif

    !--------------------------------------------------------------------------
    !- Section 6.3.1: dynamics
    !--------------------------------------------------------------------------

#ifdef COSMOART
    ! CK 20101204 unit conversion necessary before dynamics
    ! more universal approach: get a general injection point for ART
    IF (l_cosmo_art) THEN
      CALL organize_cosmo_art ('prepare_dynamics', ydate_ini, izerror, yzerrmsg)
      IF (izerror /= 0_iintegers) THEN
        CALL model_abort (my_cart_id, 100+izerror, yzerrmsg,               &
                                       'organize_cosmo_art: prepare_dynamics')
      ENDIF
    ENDIF
#endif

#ifdef POLLEN
    ! Preparations for Pollen
    IF (l_pollen) THEN
      CALL organize_pollen ('prepare_dynamics', ydate_ini, izerror, yzerrmsg)
      IF (izerror /= 0_iintegers) THEN
        CALL model_abort (my_cart_id, 100+izerror, yzerrmsg,               &
                                       'organize_pollen: prepare_dynamics')
      ENDIF
    ENDIF
#endif

    CALL organize_dynamics ('compute', izerror, yzerrmsg, dt, .FALSE.)
    IF (izerror /= 0_iintegers) THEN
      CALL model_abort (my_cart_id, 100+izerror, yzerrmsg,               &
                                    'organize_dynamics: compute')
    ENDIF

#ifdef COSMOART
    ! CK 20101204 unit conversion necessary after dynamics
    ! more universal approach: get a general injection point for ART
    IF (l_cosmo_art) THEN
      CALL organize_cosmo_art ('finalize_dynamics', ydate_ini, izerror, yzerrmsg)
      IF (izerror /= 0_iintegers) THEN
        CALL model_abort (my_cart_id, 100+izerror, yzerrmsg,               &
                                       'organize_cosmo_art: finalize_dynamics')
      ENDIF
    ENDIF
#endif

#ifdef POLLEN
    ! Clean up after Pollen and washout
    IF (l_pollen) THEN
      CALL organize_pollen ('finalize_dynamics', ydate_ini, izerror, yzerrmsg)
      IF (izerror /= 0_iintegers) THEN
        CALL model_abort (my_cart_id, 100+izerror, yzerrmsg,               &
                                       'organize_pollen: finalize_dynamics')
      ENDIF
    ENDIF
#endif

    CALL set_qrqsqg_boundaries

    IF (ltime) CALL get_timings (i_dyn_computations, ntstep, dt, izerror)

#ifdef COSMOART
    !--------------------------------------------------------------------------
    !- Section 6.3a: chemistry and aerosol dynamics, COSMO_ART; Pollen
    !--------------------------------------------------------------------------

    IF (l_cosmo_art) THEN

      IF (laero) THEN
        CALL organize_cosmo_art ('init_aero', ydate_ini,izerror, yzerrmsg)
      ENDIF

      CALL organize_cosmo_art ('compute', ydate_ini, izerror, yzerrmsg)

      IF (izerror /= 0_iintegers) THEN
        CALL model_abort (my_cart_id, 100+izerror, yzerrmsg,               &
                                       'organize_cosmo_art: chemie')
      ENDIF

    ENDIF
#endif

#ifdef POLLEN
    IF (l_pollen) THEN
      WHERE (cpollen(:,:,:,:,nnew) < 1.0E-15_ireals)
        cpollen(:,:,:,:,nnew) = 1.0E-14_ireals
      ENDWHERE
    ENDIF
#endif

    !--------------------------------------------------------------------------
    !- Section 6.4: nudging
    !--------------------------------------------------------------------------

#ifdef NUDGING
    IF (luseobs) CALL organize_assimilation ('nudge', izerror, yzerrmsg)
    IF (izerror /= 0_iintegers) THEN
      CALL model_abort (my_cart_id, 100+izerror, yzerrmsg,               &
                                    'organize_assimilation: nudge')
    ENDIF
#endif

    !--------------------------------------------------------------------------
    !- Section 6.4a: latent heat nudging (LHN)
    !--------------------------------------------------------------------------

#ifdef NUDGING
    IF (luseobs) CALL organize_assimilation ('lhn', izerror, yzerrmsg)
    IF (izerror /= 0_iintegers) THEN
      CALL model_abort (my_cart_id, 100+izerror, yzerrmsg,               &
                                    'organize_assimilation: lhn')
    ENDIF

    IF (ltime) CALL get_timings (i_lhn_computations, ntstep, dt, izerror)
#endif

    !--------------------------------------------------------------------------
    !- Section 6.5: water budget
    !--------------------------------------------------------------------------
 
    IF (ldiagnos .AND. (l2tls .OR. (ntstep > 0))) THEN
      ! for the leapfrog scheme the summations in diagbudget must not be done
      ! in the first intermediate step ntstep==0. These calculations are done
      ! again for ntstep==1.
      CALL organize_diagnosis ('diagbudget', izerror, yzerrmsg)
      IF (izerror /= 0_iintegers) THEN
        CALL model_abort (my_cart_id, 100+izerror, yzerrmsg,             &
                                     'organize_diagnosis: diagbudget')
      ENDIF
    ENDIF

    IF (ltime) CALL get_timings (i_add_computations, ntstep, dt, izerror)

    !--------------------------------------------------------------------------
    !- Section 6.6: spectral nudging and relaxation
    !--------------------------------------------------------------------------

    IF (lspecnudge) THEN
      CALL organize_dynamics ('specnudge', izerror, yzerrmsg, dt, .FALSE.)
      IF (izerror /= 0_iintegers) THEN
        CALL model_abort (my_cart_id, 100+izerror, yzerrmsg,               &
                                      'organize_dynamics: specnudge')
      ENDIF
      IF (ltime) CALL get_timings (i_spectr_nudging, ntstep, dt, izerror)
    ENDIF

#ifdef COSMOART
    ! CK 20101204 unit conversion necessary before relaxation
    ! more universal approach: get a general injection point for ART
    IF (l_cosmo_art) THEN
      CALL organize_cosmo_art ('prepare_relaxation', ydate_ini, izerror, yzerrmsg)
      IF (izerror /= 0_iintegers) THEN
        CALL model_abort (my_cart_id, 100+izerror, yzerrmsg,               &
                                       'organize_cosmo_art: prepare_relaxation')
      ENDIF
    ENDIF
#endif

#ifdef POLLEN
    ! Preparations for Pollen
    IF (l_pollen) THEN
      CALL organize_pollen ('prepare_relaxation', ydate_ini, izerror, yzerrmsg)
      IF (izerror /= 0_iintegers) THEN
        CALL model_abort (my_cart_id, 100+izerror, yzerrmsg,               &
                                       'organize_pollen: prepare_relaxation')
      ENDIF
    ENDIF
#endif

    CALL organize_dynamics ('relaxation', izerror, yzerrmsg, dt, .FALSE.)
    IF (izerror /= 0_iintegers) THEN
      CALL model_abort (my_cart_id, 100+izerror, yzerrmsg,               &
                                    'organize_dynamics: relaxation')
    ENDIF
    IF (ltime) CALL get_timings (i_relaxation, ntstep, dt, izerror)

#ifdef COSMOART
    ! CK 20101204 unit conversion necessary after dynamics
    ! more universal approach: get a general injection point for ART
    IF (l_cosmo_art) THEN
      CALL organize_cosmo_art ('finalize_relaxation', ydate_ini, izerror, yzerrmsg)
      IF (izerror /= 0_iintegers) THEN
        CALL model_abort (my_cart_id, 100+izerror, yzerrmsg,               &
                                       'organize_cosmo_art: finalize_relaxation')
      ENDIF
    ENDIF
#endif

#ifdef POLLEN
    ! Clean up after Pollen
    IF (l_pollen) THEN
      CALL organize_pollen ('finalize_relaxation', ydate_ini, izerror, yzerrmsg)
      IF (izerror /= 0_iintegers) THEN
        CALL model_abort (my_cart_id, 100+izerror, yzerrmsg,               &
                                       'organize_pollen: finalize_relaxation')
      ENDIF
    ENDIF
#endif

    ! Final update of temperature and humidity variables due to
    ! cloud microphysics in case of the cloud ice scheme
    IF (lphys) CALL organize_physics ('finish_compute', izerror, yzerrmsg)

    IF (izerror /= 0_iintegers) THEN
      CALL model_abort (my_cart_id, 100+izerror, yzerrmsg,               &
                                     'organize_physics: finish_compute')
    ENDIF

    CALL nullify_tracers

#ifdef COSMOART
    ! CK 20101204 setting minima for tracers
    ! more universal approach: get a general injection point for ART
    IF (l_cosmo_art) THEN
      CALL organize_cosmo_art ('finalize', ydate_ini, izerror, yzerrmsg)
      IF (izerror /= 0_iintegers) THEN
        CALL model_abort (my_cart_id, 100+izerror, yzerrmsg,               &
                                       'organize_cosmo_art: finalize')
      ENDIF
    ENDIF
#endif

    IF (ltime) CALL get_timings (i_phy_computations, ntstep, dt, izerror)

    !--------------------------------------------------------------------------
    !- Section 6.7: Exchange of boundary data
    !--------------------------------------------------------------------------
 
      ! The calls to exchg_boundaries have to be here in any case, even for a
      ! sequential version, because of possible periodic boundary conditions 
      ! the check, which kind of communication is necessary, is done within 
      ! the subroutine now.

      ! Check, whether additional communication for the convection is
      ! necessary
      lzconv = lconv .AND. lconf_avg .AND.                         &
                ((ntstep < 1) .OR. (MOD(ntstep+2,nincconv)==0))

      IF (ltime_barrier) THEN
        CALL comm_barrier (icomm_cart, ierrstat, yzerrmsg)
        IF (ltime) CALL get_timings (i_barrier_waiting_dyn, ntstep, dt, izerror)
      ENDIF

      IF     ( l2tls .AND. irunge_kutta /= 0 ) THEN
        CALL exchange_runge_kutta
      ELSEIF ( l2tls .AND. irunge_kutta == 0 ) THEN
        CALL exchange_2timelevel
      ELSE ! Leapfrog:
        CALL exchange_leapfrog
      ENDIF

      IF (ltime) CALL get_timings (i_communications_dyn, ntstep, dt, izerror)

    !--------------------------------------------------------------------------
    !- Section 6.8: diagnostics
    !--------------------------------------------------------------------------

    CALL near_surface (nnow)

    !   Analysis of near surface parameters
    !   -----------------------------------

#ifdef NUDGING
    IF (luseobs) CALL organize_assimilation ('surface', izerror, yzerrmsg)
    IF (izerror /= 0_iintegers) THEN
      CALL model_abort (my_cart_id, 100+izerror, yzerrmsg,               &
                                    'organize_assimilation: surface')
    ENDIF
#endif

    IF (ldiagnos) THEN
      CALL organize_diagnosis ('compute', izerror, yzerrmsg)
      IF (izerror /= 0_iintegers) THEN
        CALL model_abort (my_cart_id, 100+izerror, yzerrmsg,             &
                                       'organize_diagnosis: compute')
      ENDIF
    ENDIF

#if defined RTTOV7 || defined RTTOV9 || defined RTTOV10
    IF (luse_rttov) THEN
      CALL organize_satellites ('compute', izerror, yzerrmsg)
      IF (izerror /= 0_iintegers) THEN
        CALL model_abort (my_cart_id, 100+izerror, yzerrmsg,             &
                                       'organize_satellites: compute')
      ENDIF
    ENDIF
#endif

    IF (ltime) CALL get_timings (i_add_computations, ntstep, dt, izerror)

    !--------------------------------------------------------------------------
    !- Section 6.9: output of results
    !--------------------------------------------------------------------------

    zgrids_dt(1) = dt
    CALL organize_data ('result', ntstep, 1, 1, zgrids_dt,  &
                         izerror, yzerrmsg)
    IF (izerror /= 0_iintegers) THEN
      CALL model_abort (my_cart_id, 100+izerror, yzerrmsg,               &
                                     'result: input-init')
    ENDIF

    IF (ltime) CALL get_timings (i_output, ntstep, dt, izerror)

    !--------------------------------------------------------------------------
    !- Section 6.10: Finalization of this time step
    !--------------------------------------------------------------------------
 
#if defined RTTOV7 || defined RTTOV9 || defined RTTOV10
    ! deallocate the satellite variables
    IF (luse_rttov) THEN
      CALL organize_satellites ('dealloc', izerror, yzerrmsg)
      IF (izerror /= 0_iintegers) THEN
        CALL model_abort (my_cart_id, 100+izerror, yzerrmsg,             &
                                       'organize_satellites: dealloc')
      ENDIF
    ENDIF
#endif

    IF (lbdclim) THEN
      ! record a forecast day
      IF (my_cart_id ==  0) THEN
        zforecasttime = zforecasttime + dt
        IF (zforecasttime >= 86400.0) THEN
          nzdays  = NINT ((ntstep+1)*dt) / 86400
          WRITE (ynote(29:34),'(I6.6)') nzdays
          PRINT *, ynote
          zforecasttime = zforecasttime - 86400.0_ireals
        ENDIF
      ENDIF
    ELSE
      ! record a forecast hour
      IF (my_cart_id ==  0) THEN
        zforecasttime = zforecasttime + dt
        IF (zforecasttime >= 3600.0) THEN
          nzhours = NINT ((ntstep+1)*dt) / 3600
          WRITE (ynote(29:31),'(I3.3)') nzhours
          PRINT *, ynote
          zforecasttime = zforecasttime - 3600.0_ireals
        ENDIF
      ENDIF
    ENDIF

    ! Reset the time step for leapfrog integration
    IF ( ntstep == 0 .AND. (.NOT.l2tls) ) THEN
      dt = 2.0 * dt
    ENDIF

#if defined COUP_OAS_COS
! OASIS4 only
!    CALL oas_cos_update_time(ntstep+1)
#endif

  ENDDO timeloop

  IF (izdebug > 0) THEN
    PRINT *, 'END OF TIME STEPPING'
  ENDIF

!------------------------------------------------------------------------------
!- Section 7: Final clean up
!------------------------------------------------------------------------------

  IF (izdebug > 0) THEN
    PRINT *, 'CLEAN UP'
  ENDIF

  CALL organize_allocation ('dealloc', ierrstat)
  IF (ldiagnos) THEN
    CALL organize_diagnosis ('dealloc', izerror, yzerrmsg)
    IF (izerror /= 0_iintegers) THEN
      CALL model_abort (my_cart_id, 100+izerror, yzerrmsg,             &
                                     'organize_diagnosis: dealloc')
    ENDIF
  ENDIF

#if defined RTTOV7 || defined RTTOV9 || defined RTTOV10
  IF (luse_rttov) THEN
    CALL organize_satellites ('cleanup', izerror, yzerrmsg)
    IF (izerror /= 0_iintegers) THEN
       CALL model_abort (my_cart_id, 100+izerror, yzerrmsg,             &
          'organize_satellites: cleanup')
    ENDIF
  ENDIF
#endif

#ifdef COSMOART
  IF (l_cosmo_art) THEN
    CALL organize_cosmo_art ('deallocate', ydate_ini, izerror, yzerrmsg)

    IF (izerror /= 0) THEN
      ierrstat = 1005
      yzerrmsg  = ' ERROR    *** Deallocation of space for COSMO_ART failed ***'
      CALL model_abort (my_cart_id, ierrstat, yzerrmsg,             &
                                    'organize_cosmo_art: deallocate')
    ENDIF
  ENDIF
#endif

#ifdef POLLEN
  IF (l_pollen) THEN
    CALL organize_pollen ('deallocate', ydate_ini, izerror, yzerrmsg)

    IF (izerror /= 0) THEN
      ierrstat = 1005
      yzerrmsg  = ' ERROR    *** Deallocation of space for Pollen failed ***'
      CALL model_abort (my_cart_id, ierrstat, yzerrmsg,             &
                                    'organize_pollen: deallocate')
    ENDIF
  ENDIF
#endif

  IF (ltime) THEN
    CALL get_timings (i_cleanup, nstop, dt, izerror)
    CALL collect_timings
  ENDIF

  IF (lasync_io .OR. (num_compute > 1) ) THEN
    CALL mpe_io_shutdown()
  ENDIF

!------------------------------------------------------------------------------
!- Section 8: Part of the IO-PEs
!------------------------------------------------------------------------------

ELSE comp_pe

   CALL mpe_io_node()

ENDIF comp_pe

!------------------------------------------------------------------------------
!- Section 9: Final MPI-cleanup
!------------------------------------------------------------------------------

  CALL final_environment (ierrstat, yzerrmsg)

!------------------------------------------------------------------------------
!- End of the main program
!------------------------------------------------------------------------------

!==============================================================================
! Internal procedures in lmorg
!==============================================================================

CONTAINS

!==============================================================================
!+ Internal procedure in "lmorg" for initializing each time step
!------------------------------------------------------------------------------

!option! -pvctl ifopt
SUBROUTINE initialize_loop (ntstep, nbd1, nbd2, nold, nnow, nnew)

!------------------------------------------------------------------------------
!
! Description:
!   This routine initializes each time step. It checks whether certain 
!   actions have to be performed and sets the logical variables from
!   the parameterlist. Organizational variables are updated.
!
! Method:
!   - input of new boundary values, if necessary
!   - check whether diagnostic computations have to be performed
!   - check whether outputs have to be done
!   - update organiziational variables
!   - initialize the new time level with boundary values
!   - calculate the density of moist air (rho) for this time level
!
!==============================================================================
!
! Declarations must be of the form:
! <type> <VariableName> ! Description/ purpose of variable
!
! Subroutine / Function arguments
! Scalar arguments with intent(in):
INTEGER (KIND=iintegers), INTENT (IN)     ::       &
  ntstep             ! actual time step

! Scalar arguments with intent(inout):
INTEGER (KIND=iintegers), INTENT (INOUT)  ::       &
  nbd1, nbd2,      & ! indices for the boundary levels
  nold, nnow, nnew   ! indices for the time levels

!------------------------------------------------------------------------------
! Local variables

INTEGER (KIND=iintegers)   ::  &
  nsp, nx0, i, j, k, i1,       & ! support variable
  nzjulianday                    ! day of the year

REAL (KIND=ireals)         ::        &
  tact, tlbc, tivl,            & ! temporal quantities for control output
  z1, z2,                      & ! factors for time interpolation
  fpvsw, fpvsi, fqvs, zst, zsge, zsp,  & ! for the statement functions
  zacthour                               ! actual hour of the forecast

#ifdef COSMOART
! CK 20101204
! ART bd update frequency does not need to be the same as meteo.
! therefore, the weights calculated for meteo might be wrong.

REAL (KIND=ireals)         ::        &
  z1_art, z2_art                       ! factors for time interpolation
#endif

REAL (KIND=ireals)         ::        &
  zt_s(ie,je)                    ! = t_s   on land and sea
                                 ! = t_ice on sea ice (if present)

!------------------------------------------------------------------------------
!- End of header -
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!- Begin SUBROUTINE initialize_loop
!------------------------------------------------------------------------------

! Statement functions
!--------------------
  fpvsw(zst)       = b1 * EXP( b2w * (zst-b3)/(zst-b4w) )
  fpvsi(zst)       = b1 * EXP( b2i * (zst-b3)/(zst-b4i) )
  fqvs (zsge, zsp) = rdv * zsge / ( zsp - o_m_rdv * zsge )

! get new actual utc date
  CALL get_utc_date(ntstep, ydate_ini, dt, itype_calendar, yakdat1, yakdat2, &
                    nzjulianday, zacthour)

!------------------------------------------------------------------------------
!- Section 1: input of new boundary values, if necessary
!------------------------------------------------------------------------------

  IF (ltime) CALL get_timings (i_add_computations, ntstep, dt, izerror)

  ! get new boundary data file or generate artificial data
  zgrids_dt(1) = dt
  CALL organize_data ('boundary', ntstep, 1, 1, zgrids_dt, &
                       izerror, yzerrmsg)
  IF (izerror /= 0_iintegers) THEN
    CALL model_abort (my_cart_id, 100+izerror, yzerrmsg,                &
                                   'boundary: input-init')
  ENDIF

  IF (ltime) CALL get_timings (i_input, ntstep, dt, izerror)

!------------------------------------------------------------------------------
!- Section 2: update organizational variables
!------------------------------------------------------------------------------

  ! cyclic changing of the indices for the time levels
  IF (l2tls) THEN
    nnow = 3 - nnow
    nnew = 3 - nnew
    nx0  = nnow
  ELSE
    nsp    = nold
    nold   = nnow
    nnow   = nnew
    nnew   = nsp
    nx0    = nold
  ENDIF

  ! variables concerned with the time step
  dt2    = 2.0 * dt
  ed2dt  = 1.0 / dt2
  dtdeh  = dt / 3600.0
  nehddt = NINT ( 3600.0_ireals / dt )

!------------------------------------------------------------------------------
!- Section 3: Initialize new time level with boundary values
!------------------------------------------------------------------------------

  ! factors for linear time interpolation
  z2 = REAL (ntstep+1-nlastbound, ireals) / REAL (nincbound, ireals)
  z2 = MIN ( 1.0_ireals , z2 )
  z1 = 1.0 - z2

  ! fields of the atmosphere
  IF (lbd_frame) THEN
!CDIR COLLAPSE
    WHERE (t_bd (:,:,:,nbd2) /= undef)
!CDIR COLLAPSE
      u (:,:,:,nnew) = z1 * u_bd (:,:,:,nbd1) + z2 * u_bd (:,:,:,nbd2)
!CDIR COLLAPSE
      v (:,:,:,nnew) = z1 * v_bd (:,:,:,nbd1) + z2 * v_bd (:,:,:,nbd2)
!CDIR COLLAPSE
      t (:,:,:,nnew) = z1 * t_bd (:,:,:,nbd1) + z2 * t_bd (:,:,:,nbd2)
!CDIR COLLAPSE
      pp(:,:,:,nnew) = z1 * pp_bd(:,:,:,nbd1) + z2 * pp_bd(:,:,:,nbd2)
!CDIR COLLAPSE
      qv(:,:,:,nnew) = z1 * qv_bd(:,:,:,nbd1) + z2 * qv_bd(:,:,:,nbd2)
!CDIR COLLAPSE
      qc(:,:,:,nnew) = z1 * qc_bd(:,:,:,nbd1) + z2 * qc_bd(:,:,:,nbd2)
    ELSEWHERE
!CDIR COLLAPSE
      u (:,:,:,nnew) = u (:,:,:,nnow)
!CDIR COLLAPSE
      v (:,:,:,nnew) = v (:,:,:,nnow)
!CDIR COLLAPSE
      t (:,:,:,nnew) = t (:,:,:,nnow)
!CDIR COLLAPSE
      pp(:,:,:,nnew) = pp(:,:,:,nnow)
!CDIR COLLAPSE
      qv(:,:,:,nnew) = qv(:,:,:,nnow)
!CDIR COLLAPSE
      qc(:,:,:,nnew) = qc(:,:,:,nnow)
!CDIR COLLAPSE
    ENDWHERE

!AK (20.03.12)
    DO iig=1, ntracer
      WHERE (t_bd (:,:,:,nbd2) /= undef)
        tracer(:,:,:,nnew,iig) = z1 * tracer_bd(:,:,:,nbd1,iig) + z2 * tracer_bd(:,:,:,nbd2,iig)
      ELSEWHERE
        tracer(:,:,:,nnew,iig) = tracer(:,:,:,nnow,iig)
      ENDWHERE
    ENDDO
!AK (20.03.12)

    IF (.NOT. lw_freeslip) THEN
!CDIR COLLAPSE
      WHERE (t_bd (:,:,:,nbd2) /= undef)
!CDIR COLLAPSE
        w(:,:,:,nnew) = z1 *  w_bd(:,:,:,nbd1) + z2 *  w_bd(:,:,:,nbd2)
      ELSEWHERE
!CDIR COLLAPSE
        w(:,:,:,nnew) =  w(:,:,:,nnow)
      ENDWHERE
    ENDIF
  ELSE
!CDIR COLLAPSE
    u (:,:,:,nnew) = z1 * u_bd (:,:,:,nbd1) + z2 * u_bd (:,:,:,nbd2)
!CDIR COLLAPSE
    v (:,:,:,nnew) = z1 * v_bd (:,:,:,nbd1) + z2 * v_bd (:,:,:,nbd2)
!CDIR COLLAPSE
    t (:,:,:,nnew) = z1 * t_bd (:,:,:,nbd1) + z2 * t_bd (:,:,:,nbd2)
!CDIR COLLAPSE
    pp(:,:,:,nnew) = z1 * pp_bd(:,:,:,nbd1) + z2 * pp_bd(:,:,:,nbd2)
!CDIR COLLAPSE
    qv(:,:,:,nnew) = z1 * qv_bd(:,:,:,nbd1) + z2 * qv_bd(:,:,:,nbd2)
!CDIR COLLAPSE
    qc(:,:,:,nnew) = z1 * qc_bd(:,:,:,nbd1) + z2 * qc_bd(:,:,:,nbd2)

!AK (20.03.12)
    DO iig=1, ntracer
      tracer(:,:,:,nnew,iig) = z1 * tracer_bd(:,:,:,nbd1,iig) + z2 * tracer_bd(:,:,:,nbd2,iig)
    ENDDO
!AK (20.03.12)

 IF (.NOT. lw_freeslip) THEN
!CDIR COLLAPSE
      w (:,:,:,nnew) = z1 * w_bd (:,:,:,nbd1) + z2 * w_bd (:,:,:,nbd2)
    END IF
  ENDIF

#ifdef COSMOART
  IF (l_cosmo_art)  THEN
    ! CK 20101204 ART calculates its own interpolation weights
    IF(lgasbd .OR. laerobd) THEN
      ! factors for linear time interpolation
      z2_art = REAL (ntstep+1-nlastbound_art, ireals) / &
               REAL (nincbound_art, ireals)
      z2_art = MIN ( 1.0_ireals , z2_art )
      z1_art = 1.0 - z2_art
    ENDIF
    IF (lgas) THEN
      IF (lgasbd)  THEN
!CDIR COLLAPSE
        ! CK 20101204 ART uses its own interpolation weights
        cgas(:,:,:,:,nnew) = z1_art*cgas_bd (:,:,:,:,nbd1_art) + &
                             z2_art*cgas_bd (:,:,:,:,nbd2_art)
      ELSE
!CDIR COLLAPSE
        cgas(:,:,:,:,nnew)   = cgas(:,:,:,:,nnow)
      ENDIF
!CDIR COLLAPSE
      ! CK 20101204 no more hardcoded references to certain species
      cgas(:,:,:,lho,nnew)   = cgas(:,:,:,lho,nnow)
!CDIR COLLAPSE
      cgas(:,:,:,lho2,nnew)   = cgas(:,:,:,lho2,nnow)
    ENDIF

    IF (laero)  THEN
      IF (laerobd)  THEN
!CDIR COLLAPSE
        ! CK 20101204 ART uses its own interpolation weights
        caero(:,:,:,:,nnew) = z1_art*caero_bd (:,:,:,:,nbd1_art) + &
                              z2_art*caero_bd (:,:,:,:,nbd2_art)
      ELSE
!CDIR COLLAPSE
        caero(:,:,:,:,nnew)  = caero(:,:,:,:,nnow)
      ENDIF
    ENDIF
  ENDIF
#endif

#ifdef POLLEN
!CDIR COLLAPSE
  IF (l_pollen) cpollen(:,:,:,:,nnew) = cpollen(:,:,:,:,nnow)
#endif

  ! field for the cloud ice scheme
  IF (lprog_qi) THEN
    IF (llb_qi) THEN
      IF (lbd_frame) THEN
!CDIR COLLAPSE
        WHERE (qi_bd (:,:,:,nbd2) /= undef)
!CDIR COLLAPSE
          qi(:,:,:,nnew) = z1 * qi_bd(:,:,:,nbd1) + z2 * qi_bd(:,:,:,nbd2)
        ELSEWHERE
!CDIR COLLAPSE
          qi(:,:,:,nnew) = qi(:,:,:,nnow)
        END WHERE
      ELSE
!CDIR COLLAPSE
        qi(:,:,:,nnew) = z1 * qi_bd(:,:,:,nbd1) + z2 * qi_bd(:,:,:,nbd2)
      ENDIF
    ELSE
      ! Boundary values of cloud ice are interpreted from qc and
      ! qv is recalculated from relative humidity over ice below
      ! a threshold temperature. 
      DO k = 1, ke
!CDIR COLLAPSE
        qi(:,:,k,nnew) = 0.0_ireals
!CDIR COLLAPSE
        IF ( MINVAL(t(:,:,k,nnew)) < 248.15_ireals ) THEN
          DO j = 1, je
!CDIR COLLAPSE
            DO i = 1, ie
              IF ( t(i,j,k,nnew) < 248.15_ireals ) THEN
                qi(i,j,k,nnew) = qc(i,j,k,nnew)
                qc(i,j,k,nnew) = 0.0_ireals
! this reduction is only useful for GME data without cloud ice
! put all these computations into the INT2LM later
!               qv(i,j,k,nnew) = qv(i,j,k,nnew) &
!                      *fqvs( fpvsi(t(i,j,k,nnew)), p0(i,j,k)+pp(i,j,k,nnew) ) &
!                      /fqvs( fpvsw(t(i,j,k,nnew)), p0(i,j,k)+pp(i,j,k,nnew) )
              ENDIF
            ENDDO
          ENDDO
        ENDIF
      ENDDO
    ENDIF
  ENDIF

  ! the initial values are reformulated similarily
  IF (lprog_qi) THEN
    IF ( (ntstep == 0) .AND. (.NOT. lana_qi) ) THEN
      DO k = 1, ke
!CDIR COLLAPSE
        IF ( MINVAL(t(:,:,k,nnow)) < 248.15_ireals ) THEN
          DO j = 1, je
!CDIR COLLAPSE
            DO i = 1, ie
              IF ( t(i,j,k,nnow) < 248.15_ireals ) THEN
                qi(i,j,k,nnow) = qc(i,j,k,nnow)
                qi(i,j,k,nx0)  = qi(i,j,k,nnow)
                qc(i,j,k,nnow) = 0.0_ireals
                qc(i,j,k,nx0)  = 0.0_ireals
! this reduction is only useful for GME data without cloud ice
! put all these computations into the INT2LM later
!               qv(i,j,k,nnow) = qv(i,j,k,nnow) &
!                      *fqvs( fpvsi(t(i,j,k,nnow)), p0(i,j,k)+pp(i,j,k,nnow) ) &
!                      /fqvs( fpvsw(t(i,j,k,nnow)), p0(i,j,k)+pp(i,j,k,nnow) )
!               qv(i,j,k,nx0)  = qv(i,j,k,nnow)
              ENDIF
            ENDDO
          ENDDO
        ENDIF
      ENDDO
    ENDIF
  ENDIF

  ! prognostic precipitation fields
  IF (llb_qr_qs) THEN
    IF (lbd_frame) THEN
!CDIR COLLAPSE
      WHERE (qr_bd (:,:,:,nbd2) /= undef)
!CDIR COLLAPSE
        qr(:,:,:,nnew) = z1 * qr_bd(:,:,:,nbd1) + z2 * qr_bd(:,:,:,nbd2)
!CDIR COLLAPSE
        qs(:,:,:,nnew) = z1 * qs_bd(:,:,:,nbd1) + z2 * qs_bd(:,:,:,nbd2)
      ELSEWHERE
!CDIR COLLAPSE
        qr(:,:,:,nnew) = qr(:,:,:,nnow)
!CDIR COLLAPSE
        qs(:,:,:,nnew) = qs(:,:,:,nnow)
      END WHERE
    ELSE
!CDIR COLLAPSE
      qr(:,:,:,nnew) = z1 * qr_bd(:,:,:,nbd1) + z2 * qr_bd(:,:,:,nbd2)
!CDIR COLLAPSE
      qs(:,:,:,nnew) = z1 * qs_bd(:,:,:,nbd1) + z2 * qs_bd(:,:,:,nbd2)
    ENDIF
  ENDIF
  IF (llb_qg) THEN
    IF (lbd_frame) THEN
!CDIR COLLAPSE
      WHERE (qg_bd (:,:,:,nbd2) /= undef)
!CDIR COLLAPSE
        qg(:,:,:,nnew) = z1 * qg_bd(:,:,:,nbd1) + z2 * qg_bd(:,:,:,nbd2)
      ELSEWHERE
!CDIR COLLAPSE
        qg(:,:,:,nnew) = qg(:,:,:,nnow)
      END WHERE
    ELSE
!CDIR COLLAPSE
      qg(:,:,:,nnew) = z1 * qg_bd(:,:,:,nbd1) + z2 * qg_bd(:,:,:,nbd2)
    ENDIF
  ENDIF

  ! initialize tendency fields with 0.0
  utens  (:,:,:) = 0.0_ireals
  vtens  (:,:,:) = 0.0_ireals
  wtens  (:,:,:) = 0.0_ireals
  ttens  (:,:,:) = 0.0_ireals
  qvtens (:,:,:) = 0.0_ireals
  qctens (:,:,:) = 0.0_ireals
!WS (01.06.2012) tracertens is already set to zero before call of organize_tracer_source
!  DO iig = 1, ntracer
!    tracertens(:,:,:,iig) = 0.0_ireals
!  ENDDO
!WS (01.06.2012)

  ! Initialise qi tendency so that cloud ice can be treated
  ! in the same way as cloud water.
  qitens (:,:,:) = 0.0_ireals
  pptens (:,:,:) = 0.0_ireals
  ! ... also for humidity tendency due to diffusion
  qvt_diff (:,:,:) = 0.0_ireals

#ifdef COSMOART
  IF (l_cosmo_art) THEN
    IF (lgas)   cgastens (:,:,:,:) = 0.0_ireals
    IF (laero)  caerotens(:,:,:,:) = 0.0_ireals
  ENDIF
#endif

#ifdef POLLEN
  IF (l_pollen) cpollentens(:,:,:,:) = 0.0_ireals
#endif

  ! Calculate the surface pressure ps for the new time level nnew
  CALL calps ( ps(:,:   ,nnew), pp(:,:,ke,nnew), t(:,:,ke,nnew),     &
               qv(:,:,ke,nnew), qc(:,:,ke,nnew), qrs(:,:,ke)   ,     &
               rho0(:,:,ke), p0(:,:,ke), dp0(:,:,ke),                &
               ie, je, rvd_m_o, r_d,                                 &
               istartpar, iendpar, jstartpar, jendpar )

  ! surface fields
  IF (lbdclim) THEN
    DO j = 1,je
!CDIR COLLAPSE
      DO i = 1,ie
        vio3(i,j) = z1 * vio3_bd(i,j,nbd1) + z2 * vio3_bd(i,j,nbd2)
        hmo3(i,j) = z1 * hmo3_bd(i,j,nbd1) + z2 * hmo3_bd(i,j,nbd2)

        IF ( llandmask(i,j) .EQV. .TRUE. ) THEN
          ! in this case no distinction between undef/defined (for frames)
          ! is done, because we assume that for climate simulations always
          ! the full field is given
          lai   (i,j) = z1 * lai_bd   (i,j,nbd1) + z2 * lai_bd   (i,j,nbd2)
          rootdp(i,j) = z1 * rootdp_bd(i,j,nbd1) + z2 * rootdp_bd(i,j,nbd2)
          plcov (i,j) = z1 * plcov_bd (i,j,nbd1) + z2 * plcov_bd (i,j,nbd2)
          IF (.NOT. lmulti_layer) THEN
            t_cl  (i,j) = z1 * t_cl_bd  (i,j,nbd1) + z2 * t_cl_bd  (i,j,nbd2)
            w_cl  (i,j) = z1 * w_cl_bd  (i,j,nbd1) + z2 * w_cl_bd  (i,j,nbd2)
          ENDIF
        ENDIF
      ENDDO
    ENDDO
  ENDIF

  DO j = 1,je
!CDIR COLLAPSE
    DO i = 1,ie
      IF ( llandmask(i,j) .EQV. .TRUE. ) THEN

        IF (t_snow_bd   (i,j,nbd2) /= undef) THEN

          IF(lmulti_layer .AND. lmulti_snow) THEN
            t_snow_mult(i,j,1,nnew) = z1*t_snow_bd(i,j,nbd1) + z2*t_snow_bd(i,j,nbd2)
          ELSE
            t_snow(i,j,nnew) = z1*t_snow_bd(i,j,nbd1) + z2*t_snow_bd(i,j,nbd2)
          ENDIF

          qv_s  (i,j,nnew) = z1*qv_s_bd  (i,j,nbd1) + z2*qv_s_bd  (i,j,nbd2)
          w_snow(i,j,nnew) = z1*w_snow_bd(i,j,nbd1) + z2*w_snow_bd(i,j,nbd2)
          IF(.NOT.lmulti_layer) THEN
            t_s   (i,j,nnew) = z1*t_s_bd   (i,j,nbd1) + z2*t_s_bd (i,j,nbd2)
            t_m   (i,j,nnew) = z1*t_m_bd   (i,j,nbd1) + z2*t_m_bd (i,j,nbd2)
            w_g1  (i,j,nnew) = z1*w_g1_bd  (i,j,nbd1) + z2*w_g1_bd(i,j,nbd2)
            w_g2  (i,j,nnew) = z1*w_g2_bd  (i,j,nbd1) + z2*w_g2_bd(i,j,nbd2)
            IF ( nlgw == 3 ) THEN
              w_g3(i,j,nnew) = z1*w_g3_bd  (i,j,nbd1) + z2*w_g3_bd(i,j,nbd2)
            ENDIF
          ELSE
            t_s   (i,j,nnew) = t_s   (i,j,nx0) ! should only be used, if lphys=.F.
          ENDIF
        ELSE

          IF(lmulti_layer .AND. lmulti_snow) THEN
            t_snow_mult(i,j,1,nnew) = t_snow_mult(i,j,1,nnow)
          ELSE
            t_snow(i,j,nnew) = t_snow(i,j,nnow)
          ENDIF

          qv_s  (i,j,nnew) = qv_s  (i,j,nnow)
          w_snow(i,j,nnew) = w_snow(i,j,nnow)
          IF(.NOT.lmulti_layer) THEN
            t_s   (i,j,nnew) = t_s (i,j,nnow)
            t_m   (i,j,nnew) = t_m (i,j,nnow)
            w_g1  (i,j,nnew) = w_g1(i,j,nnow)
            w_g2  (i,j,nnew) = w_g2(i,j,nnow)
            IF ( nlgw == 3 ) THEN
              w_g3(i,j,nnew) = w_g3(i,j,nnow)
            ENDIF
          ELSE
            t_s   (i,j,nnew) = t_s   (i,j,nx0) ! should only be used, if lphys=.F.
          ENDIF
        ENDIF
      ENDIF
    ENDDO
  ENDDO

  DO j = 1,je
!CDIR COLLAPSE
    DO i = 1,ie
      IF ( llandmask(i,j) .EQV. .FALSE. ) THEN

        ! For lakes: compute saturation specific humidity over the lake
        ! surface at "nnow".
        ! This is required to compute surface fluxes in "flake_interface".

        IF (llake) THEN
          IF (depth_lk(i,j) > 0.0_ireals) THEN
            ! Lake model is used and this is a lake point
            t_s   (i,j,nnew) = t_s   (i,j,nnow)
            IF(lmulti_layer .AND. lmulti_snow) THEN
              t_snow_mult(i,j,1,nnew) = t_snow_mult(i,j,1,nnow)
            ELSE
              t_snow(i,j,nnew) = t_snow(i,j,nnow)
            ENDIF
            IF (h_ice(i,j,nnow) < h_Ice_min_flk) THEN
              ! Water surface
              qv_s(i,j,nnow) = fqvs ( fpvsw ( t_s(i,j,nnow) ) , ps (i,j,nnow) )
            ELSE
              ! Ice surface
              qv_s(i,j,nnow) = fqvs ( fpvsi ( t_s(i,j,nnow) ) , ps (i,j,nnow) )
            END IF
            qv_s  (i,j,nnew) = qv_s  (i,j,nnow)
          ELSE
            ! This is a sea point
            IF (lbdclim) THEN
              ! in this case the sea surface temperature must not be constant!!
              IF(lmulti_layer .AND. lmulti_snow) THEN
                t_snow_mult(i,j,1,nnew) = z1*t_snow_bd(i,j,nbd1) + z2*t_snow_bd(i,j,nbd2)
              ELSE
                t_snow(i,j,nnew) = z1*t_snow_bd(i,j,nbd1) + z2*t_snow_bd(i,j,nbd2)
              ENDIF
              t_s   (i,j,nnew) = z1*t_s_bd   (i,j,nbd1) + z2*t_s_bd   (i,j,nbd2)
            ELSE
              IF(lmulti_layer .AND. lmulti_snow) THEN
                t_snow_mult(i,j,1,nnew) = t_snow_mult(i,j,1,nx0)
              ELSE
                t_snow(i,j,nnew) = t_snow(i,j,nx0)
              ENDIF
              t_s   (i,j,nnew) = t_s   (i,j,nx0)
            ENDIF
            IF (lseaice) THEN
              IF (h_ice(i,j,nnow) > 0.0_ireals) THEN
                ! Ice surface
                qv_s(i,j,nnow) = fqvs ( fpvsi ( t_ice(i,j,nnow) ) , ps (i,j,nnow) )
              ELSE
                ! Water surface
                qv_s(i,j,nnow) = fqvs ( fpvsw ( t_s(i,j,nnow) ) , ps (i,j,nnow) )
              ENDIF
            ELSE
              qv_s  (i,j,nx0 ) = fqvs ( fpvsw ( t_s(i,j,nx0) ) , ps (i,j,nx0) )
              qv_s  (i,j,nnow) = qv_s  (i,j,nx0)
            ENDIF
            qv_s  (i,j,nnew) = qv_s  (i,j,nnow)
          ENDIF
        ELSE
          ! Lake model is not used
          IF (lbdclim) THEN
            ! in this case the sea surface temperature must not be constant!!
            IF(lmulti_layer .AND. lmulti_snow) THEN
              t_snow_mult(i,j,1,nnew) = z1*t_snow_bd(i,j,nbd1) + z2*t_snow_bd(i,j,nbd2)
            ELSE
              t_snow(i,j,nnew) = z1*t_snow_bd(i,j,nbd1) + z2*t_snow_bd(i,j,nbd2)
            ENDIF
            t_s   (i,j,nnew) = z1*t_s_bd   (i,j,nbd1) + z2*t_s_bd   (i,j,nbd2)
          ELSE
            IF(lmulti_layer .AND. lmulti_snow) THEN
              t_snow_mult(i,j,1,nnew) = t_snow_mult(i,j,1,nx0)
            ELSE
              t_snow(i,j,nnew) = t_snow(i,j,nx0)
            ENDIF
            t_s   (i,j,nnew) = t_s   (i,j,nx0)
          ENDIF
          IF (lseaice) THEN
            IF (h_ice(i,j,nnow) > 0.0_ireals) THEN
              ! Ice surface
              qv_s(i,j,nnow) = fqvs ( fpvsi ( t_ice(i,j,nnow) ) , ps (i,j,nnow) )
            ELSE
              ! Water surface
              qv_s(i,j,nnow) = fqvs ( fpvsw ( t_s(i,j,nnow) ) , ps (i,j,nnow) )
            ENDIF
          ELSE
            qv_s  (i,j,nx0 ) = fqvs ( fpvsw ( t_s(i,j,nx0) ) , ps (i,j,nx0) )
            qv_s  (i,j,nnow) = qv_s  (i,j,nx0)
          ENDIF
          qv_s  (i,j,nnew) = qv_s  (i,j,nnow)
        ENDIF
      ENDIF
    ENDDO
  ENDDO

  ! compute the temperature at the boundary soil-atmosphere
  DO j = 1,je
!CDIR COLLAPSE
    DO i = 1,ie
      zt_s(i,j) = t_s(i,j,nnew)
      IF (lseaice) THEN
        IF (.NOT. llandmask(i,j) .AND. h_ice(i,j,nx0) > 0.0_ireals) THEN
          zt_s(i,j) = t_ice(i,j,nx0)
        ENDIF
      ENDIF
    ENDDO
  ENDDO

#if !defined COUP_OAS_COS
  IF(lmulti_layer .AND. lmulti_snow) THEN
    CALL tgcom ( t_g (:,:,nnew), t_snow_mult(:,:,1,nnew), &
                 zt_s(:,:)     , w_snow(:,:,nnew), &
                 llandmask(:,:) , ie, je, cf_snow, &
                 istartpar, iendpar, jstartpar, jendpar )
  ELSE
    CALL tgcom ( t_g (:,:,nnew), t_snow(:,:,nnew), &
                 zt_s(:,:)     , w_snow(:,:,nnew), &
                 llandmask(:,:) , ie, je, cf_snow, &
                 istartpar, iendpar, jstartpar, jendpar )
  ENDIF
#endif
  ! compute density of moist air for time-level nnow        
  CALL calrho ( t(:,:,:,nnow), pp(:,:,:,nnow), qv(:,:,:,nnow), qc(:,:,:,nnow),&
                qrs, p0, rho, ie, je, ke, r_d, rvd_m_o)

! Set one or more heating rate disturbances (up to 50). Time(s), location(s) and intensity(ies) are
! determined by namelist parameters of namelist IDEAL
  IF (lartif_data) THEN
    ! Set possible heating rate disturbance(s) in the atmosphere (affects ttens and qtens):
    CALL artif_heatrate_dist(nnow)
    ! Set possible disturbance(s) in the bottom boundary condition for t_s (takes effect if lsoil=.false.):
    ! in case the soil model is turned off:
    CALL set_tempdist_bbc_ts( )
    ! Add white noise to the T and W - fields:
    CALL  add_noise_tw(nnow)
  END IF

!-------------------------------------------------------------------------------
!  Section 4: Reinitialize vmax_10m
!-------------------------------------------------------------------------------

  IF (ntstep-1 == nnextmxu) THEN
    vmax_10m (:,:) =   0.0_ireals
    vgust_dyn(:,:) =   0.0_ireals
    vgust_con(:,:) =   0.0_ireals

    ! Determine next step for re-initializing
    hlastmxu = hnextmxu
    hnextmxu = hlastmxu + hincmxu
    nlastmxu = NINT (hlastmxu * 3600.0_ireals / dt)
    nnextmxu = NINT (hnextmxu * 3600.0_ireals / dt)
  ENDIF

!------------------------------------------------------------------------------
!- End of the Subroutine
!------------------------------------------------------------------------------

END SUBROUTINE initialize_loop

!==============================================================================
!==============================================================================

SUBROUTINE exchange_leapfrog

  IF (lprog_qi .AND. lzconv .AND. .NOT. lprogprec) THEN
    kzdims(1:24) =                                                          &
       (/ke,ke,ke,ke,ke1,ke1,ke,ke,ke,ke,ke,ke,                             &
         ke,ke,ke,ke,ke,ke,1,0,0,0,0,0/)
    CALL exchg_boundaries                                                   &
       (nnew+39, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute,  &
        ie, je, kzdims, jstartpar, jendpar,                                 &
        nbl_exchg, nboundlines, my_cart_neigh, lperi_x, lperi_y, l2dim,     &
        20000+ntstep, ldatatypes, ncomm_type, izerror, yzerrmsg,            &
        u (:,:,:,nnow), u (:,:,:,nnew), v (:,:,:,nnow), v (:,:,:,nnew),     &
        w (:,:,:,nnow), w (:,:,:,nnew), t (:,:,:,nnow), t (:,:,:,nnew),     &
        qv(:,:,:,nnow), qv(:,:,:,nnew), qc(:,:,:,nnow), qc(:,:,:,nnew),     &
        qi(:,:,:,nnow), qi(:,:,:,nnew), pp(:,:,:,nnow), pp(:,:,:,nnew),     &
        qrs(:,:,:)    , dqvdt(:,:,:)  , qvsflx(:,:) )
  ELSEIF (lprog_qi .AND. lzconv .AND. lprogprec) THEN
    kzdims(1:24) =                                                          &
       (/ke,ke,ke,ke,ke1,ke1,ke,ke,ke,ke,ke,ke,                             &
         ke,ke,ke,ke,ke,ke,ke,ke,ke,ke,1,0/)
    CALL exchg_boundaries                                                   &
       (nnew+39, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute,  &
        ie, je, kzdims, jstartpar, jendpar,                                 &
        nbl_exchg, nboundlines, my_cart_neigh,                              &
        lperi_x, lperi_y, l2dim,                                            &
        20000+ntstep, ldatatypes, ncomm_type, izerror, yzerrmsg,            &
        u (:,:,:,nnow), u (:,:,:,nnew), v (:,:,:,nnow), v (:,:,:,nnew),     &
        w (:,:,:,nnow), w (:,:,:,nnew), t (:,:,:,nnow), t (:,:,:,nnew),     &
        qv(:,:,:,nnow), qv(:,:,:,nnew), qc(:,:,:,nnow), qc(:,:,:,nnew),     &
        qi(:,:,:,nnow), qi(:,:,:,nnew), qr(:,:,:,nnow), qr(:,:,:,nnew),     &
        qs(:,:,:,nnow), qs(:,:,:,nnew), pp(:,:,:,nnow), pp(:,:,:,nnew),     &
        qrs(:,:,:)    , dqvdt(:,:,:)  , qvsflx(:,:) )
    IF (itype_gscp==4) THEN
      kzdims(1:24) =(/ke,ke,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
      CALL exchg_boundaries                                                 &
       ( 0, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute,       &
        ie, je, kzdims, jstartpar, jendpar,                                 &
        nbl_exchg, nboundlines, my_cart_neigh,                              &
        lperi_x, lperi_y, l2dim,                                            &
        20000+ntstep, .FALSE., ncomm_type, izerror, yzerrmsg,               &
        qg (:,:,:,nnow), qg (:,:,:,nnew) )
    ENDIF
  ELSEIF (lprog_qi .AND. .NOT. lzconv .AND. .NOT. lprogprec) THEN
    kzdims(1:24) =                                                          &
       (/ke,ke,ke,ke,ke1,ke1,ke,ke,ke,ke,ke,ke,                             &
         ke,ke,ke,ke,ke,0,0,0,0,0,0,0/)
    CALL exchg_boundaries                                                   &
       (nnew+36, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute,  &
        ie, je, kzdims, jstartpar, jendpar,                                 &
        nbl_exchg, nboundlines, my_cart_neigh,                              &
        lperi_x, lperi_y, l2dim,                                            &
        20000+ntstep, ldatatypes, ncomm_type, izerror, yzerrmsg,            &
        u (:,:,:,nnow), u (:,:,:,nnew), v (:,:,:,nnow), v (:,:,:,nnew),     &
        w (:,:,:,nnow), w (:,:,:,nnew), t (:,:,:,nnow), t (:,:,:,nnew),     &
        qv(:,:,:,nnow), qv(:,:,:,nnew), qc(:,:,:,nnow), qc(:,:,:,nnew),     &
        qi(:,:,:,nnow), qi(:,:,:,nnew), pp(:,:,:,nnow), pp(:,:,:,nnew),     &
        qrs(:,:,:)    )
  ELSEIF (lprog_qi .AND. .NOT. lzconv) THEN
    kzdims(1:24) =                                                          &
       (/ke,ke,ke,ke,ke1,ke1,ke,ke,ke,ke,ke,ke,                             &
         ke,ke,ke,ke,ke,ke,ke,ke,ke,0,0,0/)
    CALL exchg_boundaries                                                   &
       (nnew+36, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute,  &
        ie, je, kzdims, jstartpar, jendpar,                                 &
        nbl_exchg, nboundlines, my_cart_neigh,                              &
        lperi_x, lperi_y, l2dim,                                            &
        20000+ntstep, ldatatypes, ncomm_type, izerror, yzerrmsg,            &
        u (:,:,:,nnow), u (:,:,:,nnew), v (:,:,:,nnow), v (:,:,:,nnew),     &
        w (:,:,:,nnow), w (:,:,:,nnew), t (:,:,:,nnow), t (:,:,:,nnew),     &
        qv(:,:,:,nnow), qv(:,:,:,nnew), qc(:,:,:,nnow), qc(:,:,:,nnew),     &
        qi(:,:,:,nnow), qi(:,:,:,nnew), qr(:,:,:,nnow), qr(:,:,:,nnew),     &
        qs(:,:,:,nnow), qs(:,:,:,nnew), pp(:,:,:,nnow), pp(:,:,:,nnew),     &
        qrs(:,:,:)    )
    IF (itype_gscp==4) THEN
      kzdims(1:24) =(/ke,ke,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
      CALL exchg_boundaries                                                 &
       ( 0, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute,       &
        ie, je, kzdims, jstartpar, jendpar,                                 &
        nbl_exchg, nboundlines, my_cart_neigh,                              &
        lperi_x, lperi_y, l2dim,                                            &
        20000+ntstep, .FALSE., ncomm_type, izerror, yzerrmsg,               &
        qg (:,:,:,nnow), qg (:,:,:,nnew) )
    ENDIF
  ELSEIF (.NOT. lprog_qi .AND. lzconv) THEN
    kzdims(1:24) =                                                          &
       (/ke,ke,ke,ke,ke1,ke1,ke,ke,ke,ke,ke,ke,                             &
         ke,ke,ke,1,0,0,0,0,0,0,0,0/)
    CALL exchg_boundaries                                                   &
       (nnew+33, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute,  &
        ie, je, kzdims, jstartpar, jendpar,                                 &
        nbl_exchg, nboundlines, my_cart_neigh,                              &
        lperi_x, lperi_y, l2dim,                                            &
        20000+ntstep, ldatatypes, ncomm_type, izerror, yzerrmsg,            &
        u (:,:,:,nnow), u (:,:,:,nnew), v (:,:,:,nnow), v (:,:,:,nnew),     &
        w (:,:,:,nnow), w (:,:,:,nnew), t (:,:,:,nnow), t (:,:,:,nnew),     &
        qv(:,:,:,nnow), qv(:,:,:,nnew), qc(:,:,:,nnow), qc(:,:,:,nnew),     &
        pp(:,:,:,nnow), pp(:,:,:,nnew), dqvdt(:,:,:)  , qvsflx(:,:) )
    IF (lprogprec) THEN
      IF (itype_gscp > 1) THEN
        kzdims(1:24) =                                                      &
           (/ke,ke,ke,ke,ke,0,0,0,0,0,0,0,                                  &
             0,0,0,0,0,0,0,0,0,0,0,0/)
        CALL exchg_boundaries                                               &
         ( 0, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute,     &
          ie, je, kzdims, jstartpar, jendpar,                               &
          nbl_exchg, nboundlines, my_cart_neigh,                            &
          lperi_x, lperi_y, l2dim,                                          &
          20000+ntstep, .FALSE., ncomm_type, izerror, yzerrmsg,             &
          qr(:,:,:,nnow), qr(:,:,:,nnew),                                   &
          qs(:,:,:,nnow), qs(:,:,:,nnew), qrs(:,:,:) )
      ELSE
        kzdims(1:24) =                                                      &
           (/ke,ke,ke,0,0,0,0,0,0,0,0,0,                                    &
             0,0,0,0,0,0,0,0,0,0,0,0/)
        CALL exchg_boundaries                                               &
         ( 0, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute,     &
          ie, je, kzdims, jstartpar, jendpar,                               &
          nbl_exchg, nboundlines, my_cart_neigh,                            &
          lperi_x, lperi_y, l2dim,                                          &
          20000+ntstep, .FALSE., ncomm_type, izerror, yzerrmsg,             &
          qr(:,:,:,nnow), qr(:,:,:,nnew), qrs(:,:,:) )
      ENDIF
    ENDIF
  ELSE
    kzdims(1:24) =                                                          &
       (/ke,ke,ke,ke,ke1,ke1,ke,ke,ke,ke,ke,ke,                             &
         ke,ke,0,0,0,0,0,0,0,0,0,0/)
    CALL exchg_boundaries                                                   &
       (nnew+30, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute,  &
        ie, je, kzdims, jstartpar, jendpar,                                 &
        nbl_exchg, nboundlines, my_cart_neigh,                              &
        lperi_x, lperi_y, l2dim,                                            &
        20000+ntstep, ldatatypes, ncomm_type, izerror, yzerrmsg,            &
        u (:,:,:,nnow), u (:,:,:,nnew), v (:,:,:,nnow), v (:,:,:,nnew),     &
        w (:,:,:,nnow), w (:,:,:,nnew), t (:,:,:,nnow), t (:,:,:,nnew),     &
        qv(:,:,:,nnow), qv(:,:,:,nnew), qc(:,:,:,nnow), qc(:,:,:,nnew),     &
        pp(:,:,:,nnow), pp(:,:,:,nnew))
    IF (lprogprec) THEN
      IF (itype_gscp > 1) THEN
        kzdims(1:24) =                                                      &
           (/ke,ke,ke,ke,ke,0,0,0,0,0,0,0,                                  &
             0,0,0,0,0,0,0,0,0,0,0,0/)
        CALL exchg_boundaries                                               &
         ( 0, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute,     &
          ie, je, kzdims, jstartpar, jendpar,                               &
          nbl_exchg, nboundlines, my_cart_neigh,                            &
          lperi_x, lperi_y, l2dim,                                          &
          20000+ntstep, .FALSE., ncomm_type, izerror, yzerrmsg,             &
          qr(:,:,:,nnow), qr(:,:,:,nnew),                                   &
          qs(:,:,:,nnow), qs(:,:,:,nnew), qrs(:,:,:) )
      ELSE
        kzdims(1:24) =                                                      &
           (/ke,ke,ke,0,0,0,0,0,0,0,0,0,                                    &
             0,0,0,0,0,0,0,0,0,0,0,0/)
        CALL exchg_boundaries                                               &
         ( 0, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute,     &
          ie, je, kzdims, jstartpar, jendpar,                               &
          nbl_exchg, nboundlines, my_cart_neigh,                            &
          lperi_x, lperi_y, l2dim,                                          &
          20000+ntstep, .FALSE., ncomm_type, izerror, yzerrmsg,             &
          qr (:,:,:,nnow), qr (:,:,:,nnew), qrs(:,:,:) )
      ENDIF
    ENDIF
  ENDIF

#ifdef COSMOART
  IF (l_cosmo_art) THEN
    IF (lgas) THEN
      DO isp = 1,isp_gas
        kzdims(1:24) =                                                      &
           (/ke,ke,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
        CALL exchg_boundaries                                               &
           ( 0, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute,   &
            ie, je, kzdims, jstartpar, jendpar,                             &
            nbl_exchg, nboundlines, my_cart_neigh,                          &
            lperi_x, lperi_y, l2dim,                                        &
            20000+ntstep, .FALSE., ncomm_type, izerror, yzerrmsg,           &
            cgas(:,:,:,isp,nnow), cgas(:,:,:,isp,nnew))
      ENDDO
    ENDIF
    IF (laero) THEN
      DO isp = 1,isp_aero
        kzdims(1:24) =                                                      &
           (/ke,ke,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
        CALL exchg_boundaries                                               &
           ( 0, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute,   &
            ie, je, kzdims, jstartpar, jendpar,                             &
            nbl_exchg, nboundlines, my_cart_neigh,                          &
            lperi_x, lperi_y, l2dim,                                        &
            20000+ntstep, .FALSE., ncomm_type, izerror, yzerrmsg,           &
            caero(:,:,:,isp,nnow), caero(:,:,:,isp,nnew))
      ENDDO
    ENDIF
  ENDIF
#endif

#ifdef POLLEN
  IF (l_pollen) THEN
    DO isp = 1,isp_pollen
      kzdims(1:24) =                                                        &
         (/ke,ke,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
      CALL exchg_boundaries                                                 &
         ( 0, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute,     &
          ie, je, kzdims, jstartpar, jendpar,                               &
          nbl_exchg, nboundlines, my_cart_neigh,                            &
          lperi_x, lperi_y, l2dim,                                          &
          20000+ntstep, .FALSE., ncomm_type, izerror, yzerrmsg,             &
          cpollen(:,:,:,isp,nnow), cpollen(:,:,:,isp,nnew))
    ENDDO
  ENDIF
#endif

!AK (20.03.2012)
  DO iig=1, ntracer
    nprog = 0_iintegers
    DO iprog=1, 7
      nprog = nprog + ltracer(iprog,iig)
    ENDDO
    IF (nprog .GE. 1) THEN
      kzdims(1:24)=(/ke,ke,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
      CALL exchg_boundaries                                                &
        (0, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute,      &
                                                               ie, je,     &
        kzdims, jstartpar, jendpar, nbl_exchg, nboundlines, my_cart_neigh, &
        lperi_x, lperi_y, l2dim,                                           &
        20000+ntstep, .FALSE., ncomm_type, izerror, yzerrmsg,              &
        tracer (:,:,:,nnow,iig), tracer (:,:,:,nnew,iig) )
    ENDIF
  ENDDO
!AK (20.03.2012)

END SUBROUTINE exchange_leapfrog

!==============================================================================
!==============================================================================

SUBROUTINE exchange_runge_kutta
  
  IF (lprog_qi) THEN
    IF (lprogprec) THEN
      ! this is former itype_gscp = 5
      IF (itype_gscp == 3) THEN
        kzdims(1:24)=(/ke,ke,ke1,ke,ke,ke,ke,ke,ke,ke,ke,0,0,0,0,0,0,0,0,0,0,0,0,0/)
        CALL exchg_boundaries                                                  &
         (50+nnew, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute,   &
          ie, je, kzdims, jstartpar, jendpar,                                  &
          nbl_exchg, nboundlines, my_cart_neigh,                               &
          lperi_x, lperi_y, l2dim,                                             &
          20000+ntstep, ldatatypes, ncomm_type, izerror, yzerrmsg,             &
          u (:,:,:,nnew), v (:,:,:,nnew), w (:,:,:,nnew), t (:,:,:,nnew),      &
          qv(:,:,:,nnew), qc(:,:,:,nnew), qi(:,:,:,nnew), qr(:,:,:,nnew),      &
          qs(:,:,:,nnew), pp(:,:,:,nnew), qrs(:,:,:) )
      END IF
      IF (itype_gscp == 4) THEN
        kzdims(1:24)=(/ke,ke,ke1,ke,ke,ke,ke,ke,ke,ke,ke,ke,0,0,0,0,0,0,0,0,0,0,0,0/)
        CALL exchg_boundaries                                                  &
         (50+nnew, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute,   &
          ie, je, kzdims, jstartpar, jendpar,                                  &
          nbl_exchg, nboundlines, my_cart_neigh,                               &
          lperi_x, lperi_y, l2dim,                                             &
          20000+ntstep, ldatatypes, ncomm_type, izerror, yzerrmsg,             &
          u (:,:,:,nnew), v (:,:,:,nnew), w (:,:,:,nnew), t (:,:,:,nnew),      &
          qv(:,:,:,nnew), qc(:,:,:,nnew), qi(:,:,:,nnew), qr(:,:,:,nnew),      &
          qs(:,:,:,nnew), qg(:,:,:,nnew), pp(:,:,:,nnew), qrs(:,:,:) )
      ENDIF
    ELSE ! .NOT. lprogprec:
      ! this is former itype_gscp = 3
      kzdims(1:24)=(/ke,ke,ke1,ke,ke,ke,ke,ke,ke,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
      CALL exchg_boundaries                                                    &
       (50+nnew, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute,     &
        ie, je, kzdims, jstartpar, jendpar,                                    &
        nbl_exchg, nboundlines, my_cart_neigh,                                 &
        lperi_x, lperi_y, l2dim,                                               &
        20000+ntstep, ldatatypes, ncomm_type, izerror, yzerrmsg,               &
        u (:,:,:,nnew), v (:,:,:,nnew), w (:,:,:,nnew), t (:,:,:,nnew),        &
        qv(:,:,:,nnew), qc(:,:,:,nnew), qi(:,:,:,nnew), pp(:,:,:,nnew),        &
        qrs(:,:,:) )
    ENDIF
  ELSE ! .NOT. lprog_qi:
    IF (lprogprec) THEN
      IF (itype_gscp > 1) THEN        
        ! this is former itype_gscp = 4
        kzdims(1:24)=(/ke,ke,ke1,ke,ke,ke,ke,ke,ke,ke,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
        CALL exchg_boundaries                                                  &
         (50+nnew, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute,   &
          ie, je, kzdims, jstartpar, jendpar,                                  &
          nbl_exchg, nboundlines, my_cart_neigh,                               &
          lperi_x, lperi_y, l2dim,                                             &
          20000+ntstep, ldatatypes, ncomm_type, izerror, yzerrmsg,             &
          u (:,:,:,nnew), v (:,:,:,nnew), w (:,:,:,nnew), t (:,:,:,nnew),      &
          qv(:,:,:,nnew), qc(:,:,:,nnew), qr(:,:,:,nnew), qs(:,:,:,nnew),      &
          pp(:,:,:,nnew), qrs(:,:,:) )
      ELSE ! kessler_pp:
        kzdims(1:24)=(/ke,ke,ke1,ke,ke,ke,ke,ke,ke,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
        CALL exchg_boundaries                                                  &
         (50+nnew, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute,   &
          ie, je, kzdims, jstartpar, jendpar,                                  &
          nbl_exchg, nboundlines, my_cart_neigh,                               &
          lperi_x, lperi_y, l2dim,                                             &
          20000+ntstep, ldatatypes, ncomm_type, izerror, yzerrmsg,             &
          u (:,:,:,nnew), v (:,:,:,nnew), w (:,:,:,nnew), t (:,:,:,nnew),      &
          qv(:,:,:,nnew), qc(:,:,:,nnew), qr(:,:,:,nnew), pp(:,:,:,nnew),      &
          qrs(:,:,:) )          
      ENDIF          
    ELSE ! .NOT. lprogprec:
      kzdims(1:24)=(/ke,ke,ke1,ke,ke,ke,ke,ke,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
      CALL exchg_boundaries                                                    &
       (50+nnew, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute,     &
        ie, je, kzdims, jstartpar, jendpar,                                    &
        nbl_exchg, nboundlines, my_cart_neigh,                                 &
        lperi_x, lperi_y, l2dim,                                               &
        20000+ntstep, ldatatypes, ncomm_type, izerror, yzerrmsg,               &
        u (:,:,:,nnew), v (:,:,:,nnew), w (:,:,:,nnew), t (:,:,:,nnew),        &
        qv(:,:,:,nnew), qc(:,:,:,nnew), pp(:,:,:,nnew), qrs(:,:,:) )
    ENDIF
  END IF
  
  IF ( lzconv ) THEN
    IF ( lprog_tke ) THEN
      kzdims(1:24)=(/ke,1,ke1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
      CALL exchg_boundaries                                                    &
       (0, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute,           &
        ie, je, kzdims, jstartpar, jendpar,                                    &
        nbl_exchg, nboundlines, my_cart_neigh,                                 &
        lperi_x, lperi_y, l2dim,                                               &
        20000+ntstep, .FALSE., ncomm_type, izerror, yzerrmsg,                  &
        dqvdt(:,:,:), qvsflx(:,:), tke(:,:,:,nnew) )
    ELSE
      kzdims(1:24)=(/ke,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
      CALL exchg_boundaries                                                    &
       (0, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute,           &
        ie, je, kzdims, jstartpar, jendpar,                                    &
        nbl_exchg, nboundlines, my_cart_neigh,                                 &
        lperi_x, lperi_y, l2dim,                                               &
        20000+ntstep, .FALSE., ncomm_type, izerror, yzerrmsg,                  &
        dqvdt(:,:,:), qvsflx(:,:) )
    END IF
  ELSE
    IF ( lprog_tke ) THEN
      kzdims(1:24)=(/ke1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
      CALL exchg_boundaries                                                    &
       (0, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute,           &
        ie, je, kzdims, jstartpar, jendpar,                                    &
        nbl_exchg, nboundlines, my_cart_neigh,                                 &
        lperi_x, lperi_y, l2dim,                                               &
        20000+ntstep, .FALSE., ncomm_type, izerror, yzerrmsg,                  &
        tke(:,:,:,nnew) )
    END IF
  END IF

#ifdef COSMOART
! CK 20101204 boundary exchange also for ART
  IF (l_cosmo_art) THEN
    IF (lgas) THEN
      DO isp = 1,isp_gas
        kzdims(1:24) =                                                         &
           (/ke,ke,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
        CALL exchg_boundaries                                                  &
           (2, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute,       &
            ie, je, kzdims, jstartpar, jendpar,                                &
            nbl_exchg, nboundlines, my_cart_neigh,                             &
            lperi_x, lperi_y, l2dim,                                           &
            20000+ntstep, .FALSE., ncomm_type, izerror, yzerrmsg,              &
            cgas(:,:,:,isp,nnew))
      ENDDO
    ENDIF
    IF (laero) THEN
      DO isp = 1,isp_aero
        kzdims(1:24) =                                                         &
           (/ke,ke,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
        CALL exchg_boundaries                                                  &
           (2, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute,       &
            ie, je, kzdims, jstartpar, jendpar,                                &
            nbl_exchg, nboundlines, my_cart_neigh,                             &
            lperi_x, lperi_y, l2dim,                                           &
            20000+ntstep, .FALSE., ncomm_type, izerror, yzerrmsg,              &
            caero(:,:,:,isp,nnew))
      ENDDO
    ENDIF
  ENDIF
#endif

#ifdef POLLEN
  IF (l_pollen) THEN
    DO isp = 1,isp_pollen
      kzdims(1:24) =                                                           &
         (/ke,ke,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
      CALL exchg_boundaries                                                    &
         (2, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute,         &
          ie, je, kzdims, jstartpar, jendpar,                                  &
          nbl_exchg, nboundlines, my_cart_neigh,                               &
          lperi_x, lperi_y, l2dim,                                             &
          20000+ntstep, .FALSE., ncomm_type, izerror, yzerrmsg,                &
          cpollen(:,:,:,isp,nnew))
    ENDDO
  ENDIF
#endif

!AK (20.03.2012)
  DO iig=1, ntracer
    nprog = 0_iintegers
    DO iprog=1, 7
      nprog = nprog + ltracer(iprog,iig)
    ENDDO
    IF (nprog .GE. 1) THEN
      kzdims(1:24)=(/ke,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
      CALL exchg_boundaries                                                   &
        (2, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute,ie, je,  &
        kzdims, jstartpar, jendpar, nbl_exchg, nboundlines, my_cart_neigh,    &
        lperi_x, lperi_y, l2dim,                                              &
        20000+ntstep, .FALSE., ncomm_type, izerror, yzerrmsg,                 &
        tracer (:,:,:,nnew,iig) )
    ENDIF
  ENDDO
!AK (20.03.2012)


END SUBROUTINE exchange_runge_kutta

!==============================================================================
!==============================================================================

SUBROUTINE exchange_2timelevel

  IF (lprog_qi .AND. lzconv) THEN
    kzdims(1:24) =(/ke,ke,ke1,ke,ke,ke,ke,ke,ke,ke,1,0,0,0,0,0,0,0,0,0,0,0,0,0/)
    CALL exchg_boundaries                                                      &
       (nnew+39, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute,     &
        ie, je, kzdims, jstartpar, jendpar,                                    &
        nbl_exchg, nboundlines, my_cart_neigh,                                 &
        lperi_x, lperi_y, l2dim,                                               &
        20000+ntstep, ldatatypes, ncomm_type, izerror, yzerrmsg,               &
        u (:,:,:,nnew), v (:,:,:,nnew), w (:,:,:,nnew), t (:,:,:,nnew),        &
        qv(:,:,:,nnew), qc(:,:,:,nnew), qi(:,:,:,nnew), pp(:,:,:,nnew),        &
        qrs(:,:,:)    , dqvdt(:,:,:)  , qvsflx(:,:) )
  ELSEIF (lprog_qi .AND. .NOT. lzconv) THEN
    kzdims(1:24) =(/ke,ke,ke1,ke,ke,ke,ke,ke,ke,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
    CALL exchg_boundaries                                                      &
       (nnew+36, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute,     &
        ie, je, kzdims, jstartpar, jendpar,                                    &
        nbl_exchg, nboundlines, my_cart_neigh,                                 &
        lperi_x, lperi_y, l2dim,                                               &
        20000+ntstep, ldatatypes, ncomm_type, izerror, yzerrmsg,               &
        u (:,:,:,nnew), v (:,:,:,nnew), w (:,:,:,nnew), t (:,:,:,nnew),        &
        qv(:,:,:,nnew), qc(:,:,:,nnew), qi(:,:,:,nnew), pp(:,:,:,nnew),        &
        qrs(:,:,:)    )
  ELSEIF (.NOT. lprog_qi .AND. lzconv) THEN
    kzdims(1:24) =(/ke,ke,ke1,ke,ke,ke,ke,ke,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
    CALL exchg_boundaries                                                      &
       (nnew+33, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute,     &
        ie, je, kzdims, jstartpar, jendpar,                                    &
        nbl_exchg, nboundlines, my_cart_neigh,                                 &
        lperi_x, lperi_y, l2dim,                                               &
        20000+ntstep, ldatatypes, ncomm_type, izerror, yzerrmsg,               &
        u (:,:,:,nnew), v (:,:,:,nnew), w (:,:,:,nnew), t (:,:,:,nnew),        &
        qv(:,:,:,nnew), qc(:,:,:,nnew), pp(:,:,:,nnew),                        &
        dqvdt(:,:,:)  , qvsflx(:,:) )
  ELSE
    kzdims(1:24) =(/ke,ke,ke1,ke,ke,ke,ke,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
    CALL exchg_boundaries                                                      &
       (nnew+30, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute,     &
        ie, je, kzdims, jstartpar, jendpar,                                    &
        nbl_exchg, nboundlines, my_cart_neigh,                                 &
        lperi_x, lperi_y, l2dim,                                               &
        20000+ntstep, ldatatypes, ncomm_type, izerror, yzerrmsg,               &
        u (:,:,:,nnew), v (:,:,:,nnew), w (:,:,:,nnew), t (:,:,:,nnew),        &
        qv(:,:,:,nnew), qc(:,:,:,nnew), pp(:,:,:,nnew))
  ENDIF

!AK (20.03.2012)
  DO iig=1, ntracer
    nprog = 0_iintegers
    DO iprog=1, 7
      nprog = nprog + ltracer(iprog,iig)
    ENDDO
    IF (nprog .GE. 1) THEN
      kzdims(1:24)=(/ke,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
      CALL exchg_boundaries                                                &
        (nnew+30, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute, ie, je,&
        kzdims, jstartpar, jendpar, nbl_exchg, nboundlines, my_cart_neigh,   &
        lperi_x, lperi_y, l2dim,                                               &
        20000+ntstep, .FALSE., ncomm_type, izerror, yzerrmsg,           &
        tracer (:,:,:,nnew,iig) )
    ENDIF
  ENDDO
!AK (20.03.2012)

  IF (lprogprec) THEN
    IF (itype_gscp == 4) THEN
      kzdims(1:24) =                                                           &
         (/ke,ke,ke,ke,0,0,0,0,0,0,0,0,                                        &
           0,0,0,0,0,0,0,0,0,0,0,0/)
      CALL exchg_boundaries                                                    &
       ( 0, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute,          &
        ie, je, kzdims, jstartpar, jendpar,                                    &
        nbl_exchg, nboundlines, my_cart_neigh,                                 &
        lperi_x, lperi_y, l2dim,                                               &
       20000+ntstep, .FALSE., ncomm_type, izerror, yzerrmsg,                   &
        qr(:,:,:,nnew), qs(:,:,:,nnew), qg(:,:,:,nnew), qrs(:,:,:) )
    ELSEIF (itype_gscp > 1) THEN
      kzdims(1:24) =                                                           &
         (/ke,ke,ke,0,0,0,0,0,0,0,0,0,                                         &
           0,0,0,0,0,0,0,0,0,0,0,0/)
      CALL exchg_boundaries                                                    &
       ( 0, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute,          &
        ie, je, kzdims, jstartpar, jendpar,                                    &
        nbl_exchg, nboundlines, my_cart_neigh,                                 &
        lperi_x, lperi_y, l2dim,                                               &
       20000+ntstep, .FALSE., ncomm_type, izerror, yzerrmsg,                   &
        qr(:,:,:,nnew), qs(:,:,:,nnew), qrs(:,:,:) )
    ELSE
      kzdims(1:24) =                                                           &
         (/ke,ke,0,0,0,0,0,0,0,0,0,0,                                          &
           0,0,0,0,0,0,0,0,0,0,0,0/)
      CALL exchg_boundaries                                                    &
       ( 0, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute,          &
        ie, je, kzdims, jstartpar, jendpar,                                    &
        nbl_exchg, nboundlines, my_cart_neigh,                                 &
        lperi_x, lperi_y, l2dim,                                               &
        20000+ntstep, .FALSE., ncomm_type, izerror, yzerrmsg,                  &
        qr(:,:,:,nnew), qrs(:,:,:) )
    ENDIF
  ENDIF

END SUBROUTINE exchange_2timelevel

!==============================================================================
!==============================================================================

SUBROUTINE exchange_l2dim

  DO k = 1, ke
    DO j = 1,nboundlines
          
      t (:,jstart-j,k,nnew) = t (:,jend  +1-j,k,nnew)
      pp(:,jstart-j,k,nnew) = pp(:,jend  +1-j,k,nnew)
      qv(:,jstart-j,k,nnew) = qv(:,jend  +1-j,k,nnew)
      qc(:,jstart-j,k,nnew) = qc(:,jend  +1-j,k,nnew)
      qrs(:,jstart-j,k)     = qrs(:,jend  +1-j,k)
      t (:,jend  +j,k,nnew) = t (:,jstart-1+j,k,nnew)
      pp(:,jend  +j,k,nnew) = pp(:,jstart-1+j,k,nnew)
      qv(:,jend  +j,k,nnew) = qv(:,jstart-1+j,k,nnew)
      qc(:,jend  +j,k,nnew) = qc(:,jstart-1+j,k,nnew)
      qrs (:,jend  +j,k)    = qrs (:,jstart-1+j,k)

!AK (20.03.2012)
      DO iig=1, ntracer
        tracer(:,jstart-j,k,nnew,iig) = tracer(:,jend  +1-j,k,nnew,iig)
        tracer(:,jend  +j,k,nnew,iig) = tracer(:,jstart-1+j,k,nnew,iig)
      ENDDO
!AK (20.03.2012)

      IF (lprogprec) THEN
        qr(:,jstart-j,k,nnew) = qr(:,jend  +1-j,k,nnew)
        qr(:,jend  +j,k,nnew) = qr(:,jstart-1+j,k,nnew)
        IF (itype_gscp > 1) THEN
          qs(:,jstart-j,k,nnew) = qs(:,jend  +1-j,k,nnew)
          qs(:,jend  +j,k,nnew) = qs(:,jstart-1+j,k,nnew)
          IF (itype_gscp == 4) THEN
            qg(:,jstart-j,k,nnew) = qg(:,jend  +1-j,k,nnew)
            qg(:,jend  +j,k,nnew) = qg(:,jstart-1+j,k,nnew)
          ENDIF
        ENDIF
      ENDIF
      IF (lprog_qi) THEN
        qi(:,jstart-j,k,nnew) = qi(:,jend  +1-j,k,nnew)
        qi(:,jend  +j,k,nnew) = qi(:,jstart-1+j,k,nnew)
      ENDIF

      IF (lzconv) THEN
        dqvdt(:,jstart-j,k)   = dqvdt(:,jend  +1-j,k)
        dqvdt(:,jend  +j,k)   = dqvdt(:,jstart-1+j,k)
      ENDIF

      u(:,jstartu-j,k,nnew) = u(:,jendu  +1-j,k,nnew)
      v(:,jstartv-j,k,nnew) = v(:,jendv  +1-j,k,nnew)
      u(:,jendu  +j,k,nnew) = u(:,jstartu-1+j,k,nnew)
      v(:,jendv  +j,k,nnew) = v(:,jstartv-1+j,k,nnew)
      
    ENDDO
  ENDDO

  DO k = 1, ke1
    DO j = 1,nboundlines
      w (:,jstart-j,k,nnew) = w (:,jend  +1-j,k,nnew)
      w (:,jend  +j,k,nnew) = w (:,jstart-1+j,k,nnew)
    ENDDO
  ENDDO

  IF (lzconv) THEN
    DO j = 1,nboundlines
      qvsflx (:,jstart-j) = qvsflx (:,jend  +1-j)
      qvsflx (:,jend  +j) = qvsflx (:,jstart-1+j)
    ENDDO
  ENDIF
     
  IF ( .NOT.l2tls ) THEN

    DO k = 1, ke
      DO j = 1,nboundlines
        
        t (:,jstart-j,k,nnow) = t (:,jend  +1-j,k,nnow)
        pp(:,jstart-j,k,nnow) = pp(:,jend  +1-j,k,nnow)
        qv(:,jstart-j,k,nnow) = qv(:,jend  +1-j,k,nnow)
        qc(:,jstart-j,k,nnow) = qc(:,jend  +1-j,k,nnow)
        t (:,jend  +j,k,nnow) = t (:,jstart-1+j,k,nnow)
        pp(:,jend  +j,k,nnow) = pp(:,jstart-1+j,k,nnow)
        qv(:,jend  +j,k,nnow) = qv(:,jstart-1+j,k,nnow)
        qc(:,jend  +j,k,nnow) = qc(:,jstart-1+j,k,nnow)

        IF (lprogprec) THEN
          qr(:,jstart-j,k,nnow) = qr(:,jend  +1-j,k,nnow)
          qr(:,jend  +j,k,nnow) = qr(:,jstart-1+j,k,nnow)
          IF (itype_gscp > 1) THEN
            qs(:,jstart-j,k,nnow) = qs(:,jend  +1-j,k,nnow)
            qs(:,jend  +j,k,nnow) = qs(:,jstart-1+j,k,nnow)
            IF (itype_gscp == 4) THEN
              qg(:,jstart-j,k,nnow) = qg(:,jend  +1-j,k,nnow)
              qg(:,jend  +j,k,nnow) = qg(:,jstart-1+j,k,nnow)
            ENDIF
          ENDIF
        ENDIF
        IF (lprog_qi) THEN
          qi(:,jstart-j,k,nnow) = qi(:,jend  +1-j,k,nnow)
          qi(:,jend  +j,k,nnow) = qi(:,jstart-1+j,k,nnow)
        ENDIF

        u(:,jstartu-j,k,nnow) = u(:,jendu  +1-j,k,nnow)
        v(:,jstartv-j,k,nnow) = v(:,jendv  +1-j,k,nnow)
        u(:,jendu  +j,k,nnow) = u(:,jstartu-1+j,k,nnow)
        v(:,jendv  +j,k,nnow) = v(:,jstartv-1+j,k,nnow)
        
      ENDDO
    ENDDO

    DO k = 1, ke1
      DO j = 1,nboundlines
        w (:,jstart-j,k,nnow) = w (:,jend  +1-j,k,nnow)
        w (:,jend  +j,k,nnow) = w (:,jstart-1+j,k,nnow)
      ENDDO
    ENDDO

  ENDIF
      
END SUBROUTINE exchange_l2dim

!==============================================================================

SUBROUTINE set_qrqsqg_boundaries

  ! Now we have to set the nnew values for qr and qs in a consistent way:
  ! this is an intermediate solution, as long as no better treatment of 
  ! the boundary values is found

  ! Treatment of rain and snow
  ! --------------------------

!US-JF:  IF (.NOT.llb_qr_qs) THEN
    IF     (itype_lbc_qrsg == 1) THEN
      ! set boundary value to first interior row

      ! qr
      IF (ALLOCATED(qr)) THEN
        ! western boundary
        IF (my_cart_neigh(1) == -1) THEN
          DO i = 1, nboundlines
            DO  k = 1, ke
!CDIR NOLOOPCHG
              DO  j = jstart, jend
                qr(i,j,k,nnew) = qr(istart,j,k,nnew)
              ENDDO
            ENDDO
          ENDDO
        ENDIF
        ! eastern boundary
        IF (my_cart_neigh(3) == -1) THEN
          DO i = ie-nboundlines+1, ie
            DO  k = 1, ke
!CDIR NOLOOPCHG
              DO  j = jstart, jend
                qr(i,j,k,nnew) = qr(iend  ,j,k,nnew)
              ENDDO
            ENDDO
          ENDDO
        ENDIF
        ! southern boundary
        IF (my_cart_neigh(4) == -1) THEN
          DO  k = 1, ke
            DO  j = 1, nboundlines
              qr(:,j,k,nnew) = qr(:,jstart,k,nnew)
            ENDDO
          ENDDO
        ENDIF
        ! northern boundary
        IF (my_cart_neigh(2) == -1) THEN
          DO  k = 1, ke
            DO  j = je-nboundlines+1, je
              qr(:,j,k,nnew) = qr(:,jend  ,k,nnew)
            ENDDO
          ENDDO
        ENDIF
      ENDIF

      ! qs
      IF (ALLOCATED(qs)) THEN
        ! western boundary
        IF (my_cart_neigh(1) == -1) THEN
          DO  k = 1, ke
            DO i = 1, nboundlines
!CDIR NOLOOPCHG
              DO  j = jstart, jend
                qs(i,j,k,nnew) = qs(istart,j,k,nnew)
              ENDDO
            ENDDO
          ENDDO
        ENDIF
        ! eastern boundary
        IF (my_cart_neigh(3) == -1) THEN
          DO  k = 1, ke
            DO i = ie-nboundlines+1, ie
!CDIR NOLOOPCHG
              DO  j = jstart, jend
                qs(i,j,k,nnew) = qs(iend  ,j,k,nnew)
              ENDDO
            ENDDO
          ENDDO
        ENDIF
        ! southern boundary
        IF (my_cart_neigh(4) == -1) THEN
          DO  k = 1, ke
            DO  j = 1, nboundlines
              qs(:,j,k,nnew) = qs(:,jstart,k,nnew)
            ENDDO
          ENDDO
        ENDIF
        ! northern boundary
        IF (my_cart_neigh(2) == -1) THEN
          DO  k = 1, ke
            DO  j = je-nboundlines+1, je
              qs(:,j,k,nnew) = qs(:,jend  ,k,nnew)
            ENDDO
          ENDDO
        ENDIF
      ENDIF

    ELSEIF (itype_lbc_qrsg == 2) THEN

      ! set all values to 0

      ! qr
      IF (ALLOCATED(qr)) THEN
        ! western boundary
        IF (my_cart_neigh(1) == -1) THEN
          DO i = 1, nboundlines
            DO  k = 1, ke
!CDIR NOLOOPCHG
              DO  j = jstart, jend
                qr(i,j,k,nnew) = 0.0_ireals
              ENDDO
            ENDDO
          ENDDO
        ENDIF
        ! eastern boundary
        IF (my_cart_neigh(3) == -1) THEN
          DO i = ie-nboundlines+1, ie
            DO  k = 1, ke
!CDIR NOLOOPCHG
              DO  j = jstart, jend
                qr(i,j,k,nnew) = 0.0_ireals
              ENDDO
            ENDDO
          ENDDO
        ENDIF
        ! southern boundary
        IF (my_cart_neigh(4) == -1) THEN
          DO  k = 1, ke
            DO  j = 1, nboundlines
              qr(:,j,k,nnew) = 0.0_ireals
            ENDDO
          ENDDO
        ENDIF
        ! northern boundary
        IF (my_cart_neigh(2) == -1) THEN
          DO  k = 1, ke
            DO  j = je-nboundlines+1, je
              qr(:,j,k,nnew) = 0.0_ireals
            ENDDO
          ENDDO
        ENDIF
      ENDIF

      ! qs
      IF (ALLOCATED(qs)) THEN
        ! western boundary
        IF (my_cart_neigh(1) == -1) THEN
          DO  k = 1, ke
            DO i = 1, nboundlines
!CDIR NOLOOPCHG
              DO  j = jstart, jend
                qs(i,j,k,nnew) = 0.0_ireals
              ENDDO
            ENDDO
          ENDDO
        ENDIF
        ! eastern boundary
        IF (my_cart_neigh(3) == -1) THEN
          DO  k = 1, ke
            DO i = ie-nboundlines+1, ie
!CDIR NOLOOPCHG
              DO  j = jstart, jend
                qs(i,j,k,nnew) = 0.0_ireals
              ENDDO
            ENDDO
          ENDDO
        ENDIF
        ! southern boundary
        IF (my_cart_neigh(4) == -1) THEN
          DO  k = 1, ke
            DO  j = 1, nboundlines
              qs(:,j,k,nnew) = 0.0_ireals
            ENDDO
          ENDDO
        ENDIF
        ! northern boundary
        IF (my_cart_neigh(2) == -1) THEN
          DO  k = 1, ke
            DO  j = je-nboundlines+1, je
              qs(:,j,k,nnew) = 0.0_ireals
            ENDDO
          ENDDO
        ENDIF
      ENDIF

!   ELSEIF (itype_lbc_qrsg == 3) THEN
!     ! nothing is done for this option
    ENDIF
!US-JF:  ENDIF      ! .NOT.llb_qr_qs

  ! Treatment of graupel
  ! --------------------

!US-JF:  IF (.NOT.llb_qg) THEN
    IF (ALLOCATED(qg)) THEN
      IF     (itype_lbc_qrsg == 1) THEN

        ! western boundary
        IF (my_cart_neigh(1) == -1) THEN
          DO  k = 1, ke
            DO i = 1, nboundlines
!CDIR NOLOOPCHG 
              DO  j = jstart, jend
                qg(i,j,k,nnew) = qg(istart,j,k,nnew)
              ENDDO
            ENDDO
          ENDDO
        ENDIF
        ! eastern boundary 
        IF (my_cart_neigh(3) == -1) THEN
          DO  k = 1, ke
            DO i = ie-nboundlines+1, ie
!CDIR NOLOOPCHG
              DO  j = jstart, jend
                qg(i,j,k,nnew) = qg(iend  ,j,k,nnew)
              ENDDO
            ENDDO
          ENDDO
        ENDIF
        ! southern boundary
        IF (my_cart_neigh(4) == -1) THEN
          DO  k = 1, ke
            DO  j = 1, nboundlines
              qg(:,j,k,nnew) = qg(:,jstart,k,nnew)
            ENDDO
          ENDDO
        ENDIF
          ! northern boundary
        IF (my_cart_neigh(2) == -1) THEN
          DO  k = 1, ke
            DO  j = je-nboundlines+1, je
              qg(:,j,k,nnew) = qg(:,jend  ,k,nnew)
            ENDDO
          ENDDO
        ENDIF

      ELSEIF (itype_lbc_qrsg == 2) THEN

        ! set all values to 0
        ! western boundary
        IF (my_cart_neigh(1) == -1) THEN
          DO  k = 1, ke
            DO i = 1, nboundlines
!CDIR NOLOOPCHG 
              DO  j = jstart, jend
                qg(i,j,k,nnew) = 0.0_ireals
              ENDDO
            ENDDO
          ENDDO
        ENDIF
        ! eastern boundary 
        IF (my_cart_neigh(3) == -1) THEN
          DO  k = 1, ke
            DO i = ie-nboundlines+1, ie
!CDIR NOLOOPCHG
              DO  j = jstart, jend
                qg(i,j,k,nnew) = 0.0_ireals
              ENDDO
            ENDDO
          ENDDO
        ENDIF
        ! southern boundary
        IF (my_cart_neigh(4) == -1) THEN
          DO  k = 1, ke
            DO  j = 1, nboundlines
              qg(:,j,k,nnew) = 0.0_ireals
            ENDDO
          ENDDO
        ENDIF
        ! northern boundary
        IF (my_cart_neigh(2) == -1) THEN
          DO  k = 1, ke
            DO  j = je-nboundlines+1, je
              qg(:,j,k,nnew) = 0.0_ireals
            ENDDO
          ENDDO
        ENDIF

!     ELSEIF (itype_lbc_qrsg == 3) THEN
!       ! nothing has to be done for this option
      ENDIF
    ENDIF
!US-JF:  ENDIF

END SUBROUTINE set_qrqsqg_boundaries

!==============================================================================

SUBROUTINE nullify_tracers

    ! If the humidity variables are rather small, set them to 0
    IF (ALLOCATED(qr)) THEN
      WHERE (qr (:,:,:,nnew) < 1.0E-15_ireals)
        qr (:,:,:,nnew) = 0.0_ireals
      ENDWHERE
    ENDIF
    IF (ALLOCATED(qs)) THEN
      WHERE (qs (:,:,:,nnew) < 1.0E-15_ireals)
        qs (:,:,:,nnew) = 0.0_ireals
      ENDWHERE
    ENDIF
    IF (ALLOCATED(qg)) THEN
      WHERE (qg (:,:,:,nnew) < 1.0E-15_ireals)
        qg (:,:,:,nnew) = 0.0_ireals
      ENDWHERE
    ENDIF
    IF (ALLOCATED(qi)) THEN
      WHERE (qi (:,:,:,nnew) < 1.0E-15_ireals)
        qi (:,:,:,nnew) = 0.0_ireals
      ENDWHERE
    ENDIF
    WHERE (qc (:,:,:,nnew) < 1.0E-15_ireals)
      qc (:,:,:,nnew) = 0.0_ireals
    ENDWHERE
    WHERE (qv (:,:,:,nnew) < 1.0E-15_ireals)
      qv (:,:,:,nnew) = 0.0_ireals
    ENDWHERE
    WHERE (qrs(:,:,:) < 1.0E-15_ireals)
      qrs(:,:,:) = 0.0_ireals
    ENDWHERE

#ifdef POLLEN
    IF (l_pollen) THEN
      WHERE (cpollen(:,:,:,:,nnew) < 1.0E-15_ireals)
        cpollen(:,:,:,:,nnew) = 1.0E-14_ireals
      ENDWHERE
    ENDIF
#endif

END SUBROUTINE nullify_tracers

!==============================================================================
! End of the program
!==============================================================================

END PROGRAM lmorg

!------------------------------------------------------------------------------
