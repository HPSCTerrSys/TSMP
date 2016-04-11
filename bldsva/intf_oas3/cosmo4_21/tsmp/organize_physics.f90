! $RCSfile: organize_physics.f90,v $
! $Revision: 4.21 $ $Date: 2012/04/03 $
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
! V4_11        2010/01/18 Volker Kuell
!  Subroutine call and additional fields for HYMACS convection scheme
! V4_11        2011/06/08 Volker Kuell
!  all global HYMACS variables exported to data_hymacs.f90 and all HYMACS related code put into
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
! V4_21        2012/03/15 Volker Kuell
!  Subroutine call and additional fields for HYMACS convection scheme
! V4_21        2012/04/03 Alexander Kelbch
!  Inclusion of tracer (adopted from Markus Uebel)
! V4_21        2013/02/03 Prabhakar Shrestha
!  Inclusion of transfere coefficeint for moisture 

!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!==============================================================================

USE data_parameters,    ONLY:  ireals, iintegers

USE data_modelconfig,   ONLY:                                                &
       ie, je, ke, ke_soil, jstartpar, jendpar, dt, dt2, czmls, msoilgrib,   &
       czhls, istartu, iendu, jstartu, jendu, istartv, iendv, jstartv, jendv,&
       ke_snow, ke1, kcm, istartpar, iendpar

USE data_runcontrol,    ONLY:                                                &
       lgsp, lrad, ltur, lconv, itype_conv, lsoil, lmelt, lmelt_var,         &
       lmulti_layer, lprogprec, ltrans_prec, ldiniprec, lprog_qi, nradcoarse,&
       lradf_avg, nincrad, nextrad, hincrad, hnextrad, ninctura, nincconv,   &
       itype_tran, itype_turb, itype_gscp, itype_synd, itype_wcld,           &
       imode_tran, imode_turb, icldm_rad, icldm_tran, icldm_turb, nstart,    &
       ntstep, l2tls, irunge_kutta, ltime, nold, nnow, nnew, nuspecif,       &
       lexpcor, lnonloc, lcpfluc, lconf_avg, lcape, lctke, ltmpcor, lprfcor, &
       itype_trvg, itype_evsl, nlgw, nbl_exchg, l3dturb,                     &
       lprog_tke, lconv_inst, lforest, luse_rttov, ntke, lseaice, llake,     &
       ico2_rad, ibot_w_so, czbot_w_so, l3dturb_metr, idbg_level,            &
       lprintdeb_all, lradtopo, nhori, hstart, lsso, nincsso,                &
       itype_sher, limpltkediff, ltkesso, ltkecon, lmulti_snow, lemiss,      &
       lstomata, itype_aerosol, itype_root, itype_heatcond, itype_hydbound,  &
       l_dzeta_d_needed, lperi_x, lperi_y, l2dim, nbl_exchg,                 &
       lartif_data, y_scalar_advect

USE data_parallel,      ONLY:                                                &
       num_compute, my_cart_id, icomm_cart, sendbuf, isendbuflen, iexch_req, &
       imp_reals, imp_integers, imp_logical, nboundlines, my_cart_neigh,     &
       nproc, my_world_id, icomm_world, intbuf, realbuf, logbuf, ncomm_type, &
       ldatatypes, ltime_barrier

USE data_fields,        ONLY:                                                &
       tkvm, tkvh, tcm, tch, qrs, t_s, t_g, t_snow, w_snow, freshsnow, qv_s, &
       w_so, w_g1, ut_conv, vt_conv, clc_con, ut_sso, vt_sso, fr_land, d_pat,&
       c_big, c_sml, r_air, plcov, sai, tai, lai, eai, hhl, h_can,           &
       tkhm, tkhh, tke, tketens

!CPS
#ifdef COUP_OAS_COS
USE data_fields,        ONLY:                                                &
       tcw
#endif

!VK 2012/03/15
USE data_tracer,        ONLY:  ntracer_conv, ntracer_dim
!VK 2012/03/15

USE data_io,            ONLY:  ydate_ini, lbdclim

USE environment,        ONLY:  exchg_boundaries, comm_barrier
USE parallel_utilities, ONLY:  distribute_values
USE time_utilities,     ONLY:  get_timings, i_precipitation, i_radiation,    &
                               i_turbulence, i_convection, i_soil_model,     &
                               i_communications_phy, i_barrier_waiting_phy,  &
                               i_add_computations, i_sso

USE src_gscp,           ONLY:  organize_gscp, hydci, hydci_pp, hydci_pp_gr,  &
                               kessler_pp, hydor_pp
USE src_radiation,      ONLY:  init_radiation, organize_radiation
!US  USE src_turbulence,     ONLY:  partura, parturs, prankolmo_rk,               &
!US                                 horizontal_diffcoeffs
!US  USE src_turbdiff,       ONLY:  init_canopy_old,                              &
!US                                 turbtran_old, turbdiff_old, turbdiff_hori_old
USE turbulence_utilities, ONLY: init_volume_canopy, init_surface_canopy
USE turbulence_interface, ONLY: organize_turbulence 
USE src_conv_tiedtke,   ONLY:  organize_conv_tiedtke
USE src_conv_kainfri,   ONLY:  organize_conv_kainfri
!not yet:  USE src_conv_bechtold,  ONLY:  organize_conv_becht
USE src_conv_shallow,   ONLY:  organize_conv_shallow

#ifdef HYMACS
!VK 2012/03/15
USE src_conv_hymacs,    ONLY:  organize_conv_hymacs
#endif

USE src_soil,           ONLY:  terra1, terra2
USE src_soil_multlay,   ONLY:  terra_multlay
USE src_flake,          ONLY:  flake_init, flake_interface
USE src_sso,            ONLY:  organize_sso
USE src_seaice,         ONLY:  seaice
USE src_artifdata,      ONLY:  tkvhfix, tkhhfix, tkvmfix, tkhmfix, &
                               lnosurffluxes_m, lnosurffluxes_h

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
                                izdebug, ist, izl
LOGICAL                    ::   lstfnct, & ! calculation of a new stability
                                           ! function in sub. turbdiff
                                lzconv     ! if convection is computed

CHARACTER (LEN= 9)         ::   yinput             ! Namelist INPUT file
CHARACTER (LEN=80)         ::   yzerror

REAL (KIND=ireals)         ::   zdt

!------------------------------------------------------------------------------
!- End of header
!------------------------------------------------------------------------------
 
!------------------------------------------------------------------------------
!- Begin Subroutine organize_physics
!------------------------------------------------------------------------------

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
! Section 2: Initialization of the packages for the first call
!------------------------------------------------------------------------------

ELSEIF (yaction == 'init') THEN

  IF (izdebug > 0) THEN
    PRINT *, '    PHYSICAL PACKAGES'
  ENDIF

  ! Initialize values for the radiation and the canopy fields
  IF (lrad) CALL init_radiation
  IF (ltur .OR. lsoil) THEN
     ! Determine the index of the highest layer which still is a canopy layer

     ! h_can is a primary external parameter. The initial values of h_can are 0.
     ! If we don't change this, no canopy will be resolved in the vertical direction.
     ! Because there is no h_can available at the moment, the next lines 
     ! are not in affect

     kcm=ke
     DO WHILE (MAXVAL( h_can - hhl(:,:,kcm) + hhl(:,:,ke1) ) > 0.0_ireals)
        kcm=kcm-1
     END DO

     ! Now kcm points to the lowest layer being not a canopy layer.
     ! Therefore add +1
     kcm=kcm+1

     ! Allocate the canopy fields c_big, c_sml, r_air:
     ist = 0
     izl = 0
     ALLOCATE ( c_big(ie,je,kcm:ke1), STAT=izl ); c_big = 0.0; ist = ist+izl
     ALLOCATE ( c_sml(ie,je,kcm:ke1), STAT=izl ); c_sml = 0.0; ist = ist+izl
     ALLOCATE ( r_air(ie,je,kcm:ke1), STAT=izl ); r_air = 0.0; ist = ist+izl

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
  ELSE
    IF (l2tls) THEN
      ntke = nnew  !???
    ELSE
      ntke = nnow  !???
    ENDIF
  ENDIF

  !Initialize varaiables of the lake model FLake
  IF (llake) CALL flake_init

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

  IF (lgsp) THEN
    IF (izdebug > 5) THEN
      PRINT *, '      GRID SCALE PRECIPITATION'
    ENDIF

    CALL organize_gscp
    IF (ltime) CALL get_timings (i_precipitation, ntstep, dt, izerror)
  ENDIF

  !----------------------------------------------------------------------------
  ! Section 3.2: Radiation
  !----------------------------------------------------------------------------

  IF (lrad) THEN
    IF ( (ntstep < 2) .OR. (ntstep == nextrad) ) THEN
      IF (izdebug > 5) THEN
        PRINT *, '      RADIATION'
      ENDIF

      CALL organize_radiation (ydate_ini)

      IF ((ntstep >= 1) .OR. (nincrad == 1)) THEN
        IF (hincrad /= 0.0) THEN
          hnextrad = hnextrad + hincrad
          nextrad  = NINT (hnextrad * 3600.0_ireals / dt)
        ELSE
          nextrad  = nextrad + nincrad
        ENDIF
      ENDIF
      IF (ltime) CALL get_timings (i_radiation, ntstep, dt, izerror)
    ENDIF
  ENDIF

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
!AK (03.04.2012)
        ! to avoid allocatable variables for tracer in convection and to guarantee
        ! that these variables have the default value 1 if ntracer_conv = 0
        ntracer_dim = MAX(1,ntracer_conv)
!AK (03.04.2012)
        CALL organize_conv_tiedtke

      CASE (1) ! Kain-Fritsch scheme
        IF (izdebug > 5) THEN
          PRINT *, '      KAIN-FRITSCH CONVECTION SCHEME'
        ENDIF
!AK (03.04.2012)
        IF (ntracer_conv /= 0) THEN
          IF ( (my_cart_id == 1) .AND. (ntstep < 2) ) THEN
            PRINT *, 'WARNING: Tracer not implemented in Kain-Fritsch!'
          ENDIF
        ENDIF
!AK (03.04.2012)
        CALL organize_conv_kainfri

!not yet:      CASE (2) ! Bechtold scheme
!not yet:        IF (izdebug > 5) THEN
!not yet:          PRINT *, '      BECHTOLD CONVECTION SCHEME'
!not yet:        ENDIF
!not yet:        CALL organize_conv_becht

      CASE (3) ! Shallow convection scheme
        IF (izdebug > 5) THEN
          PRINT *, '      SHALLOW CONVECTION SCHEME'
        ENDIF
!AK (03.04.2012)
        ! to avoid allocatable variables for tracer in convection and to guarantee
        ! that these variables have the default value 1 if ntracer_conv = 0
        ntracer_dim = MAX(1,ntracer_conv)
!AK (03.04.2012)
        CALL organize_conv_shallow

#ifdef HYMACS
! VK 2012/03/15
      CASE (4) ! HYMACS
        IF (izdebug > 5) THEN
          PRINT *, '      HYMACS CONVECTION SCHEME'
        ENDIF
!AK (03.04.2012)
        ! to avoid allocatable variables for tracer in convection and to guarantee
        ! that these variables have the default value 1 if ntracer_conv = 0
        ntracer_dim = MAX(1,ntracer_conv)
!AK (03.04.2012)
        CALL organize_conv_hymacs
#endif

      CASE DEFAULT
        ierror = 2008
        yerrmsg = 'No valid convection scheme'
        RETURN

      END SELECT

      IF (ltime) CALL get_timings (i_convection, ntstep, dt, izerror)
    ENDIF
  ENDIF
!AK (03.04.2012)
  IF (.NOT. lconv .AND. (ntracer_conv /= 0) .AND. (ntstep < 2) .AND. (my_cart_id == 1) ) THEN
    PRINT *, 'WARNING: lconv = .FALSE. (no convective transport!), ltracer(4,*) must be "0"!'
  ENDIF
!AK (03.04.2012)
  
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

!!$  IF (num_compute > 1) THEN
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
        kzdims(1:24)=(/ke,ke,ke,ke,ke,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
        CALL exchg_boundaries                                                  &
           (nx+73,  sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute, ie, je,       &
            kzdims, jstartpar, jendpar, nbl_exchg, nboundlines, my_cart_neigh, &
            lperi_x, lperi_y, l2dim, &
            10000+ntstep, ldatatypes, ncomm_type, izerror, yerrmsg,            &
            tkvm, tkvh, qrs, ut_sso, vt_sso, t_s(:,:,nx),                      &
            qv_s(:,:,nx), tcm, tch )
      ELSE
        kzdims(1:24)=(/ke,ke,ke,ke,ke,ke,ke,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0/)
        CALL exchg_boundaries                                                  &
           (nx+70,  sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute, ie, je,       &
            kzdims, jstartpar, jendpar, nbl_exchg, nboundlines, my_cart_neigh, &
            lperi_x, lperi_y, l2dim, &
            10000+ntstep, ldatatypes, ncomm_type, izerror, yerrmsg,            &
            tkvm, tkvh, qrs, ut_conv, vt_conv, ut_sso, vt_sso, t_s(:,:,nx),    &
            qv_s(:,:,nx), tcm, tch )
      ENDIF
    ELSE
      IF (.NOT. lzconv) THEN
        kzdims(1:24)=(/ke,ke,ke,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
        CALL exchg_boundaries                                                  &
           (nx+13,  sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute, ie, je,       &
            kzdims, jstartpar, jendpar, nbl_exchg, nboundlines, my_cart_neigh, &
            lperi_x, lperi_y, l2dim, &
            10000+ntstep, ldatatypes, ncomm_type, izerror, yerrmsg,            &
            tkvm, tkvh, qrs, t_s(:,:,nx), qv_s(:,:,nx), tcm, tch)
      ELSE
!CPS Added tcw for exchange
#ifdef COUP_OAS_COS
        kzdims(1:24)=(/ke,ke,ke,ke,ke,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
        CALL exchg_boundaries                                                  &
           (nx+10,  sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute, ie, je,       &
            kzdims, jstartpar, jendpar, nbl_exchg, nboundlines, my_cart_neigh, &
            lperi_x, lperi_y, l2dim, &
            10000+ntstep, ldatatypes, ncomm_type, izerror, yerrmsg,            &
            tkvm, tkvh, qrs, ut_conv, vt_conv, t_s(:,:,nx), qv_s(:,:,nx),      &
            tcm, tch, tcw)
#else
        kzdims(1:24)=(/ke,ke,ke,ke,ke,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
        CALL exchg_boundaries                                                  &
           (nx+10,  sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute, ie, je,       &
            kzdims, jstartpar, jendpar, nbl_exchg, nboundlines, my_cart_neigh, &
            lperi_x, lperi_y, l2dim, &
            10000+ntstep, ldatatypes, ncomm_type, izerror, yerrmsg,            &
            tkvm, tkvh, qrs, ut_conv, vt_conv, t_s(:,:,nx), qv_s(:,:,nx),      &
            tcm, tch)
#endif
      ENDIF
    ENDIF

! ub>> ! Added exchg. of tkhm, tkhh, tke, tketens when necessary:
    IF ( ltur ) THEN

      IF ( l3dturb ) THEN

        kzdims(1:24)=(/ke,ke,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
        CALL exchg_boundaries                                                  &
             (0,  sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute, ie, je,       &
             kzdims, jstartpar, jendpar, nbl_exchg, nboundlines, my_cart_neigh, &
             lperi_x, lperi_y, l2dim, &
             10000+ntstep, .FALSE., ncomm_type, izerror, yerrmsg,            &
             tkhm, tkhh )
      END IF
      
      SELECT CASE (itype_turb )
      CASE (3)
        IF ( l3dturb ) THEN   ! Not sure if necessary in all cases or just in 3D-Turb.
          kzdims(1:24)=(/ke,ke,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
          CALL exchg_boundaries                                                  &
               (0, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute, ie, je,       &
               kzdims, jstartpar, jendpar, nbl_exchg, nboundlines, my_cart_neigh, &
               lperi_x, lperi_y, l2dim, &
               10000+ntstep, .FALSE., ncomm_type, izerror, yerrmsg,            &
               tke, tketens )
        END IF
      CASE (5:8)
        IF ( lprog_tke ) THEN
          kzdims(1:24)=(/ke,ke,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
          CALL exchg_boundaries                                                  &
               (0, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute, ie, je,       &
               kzdims, jstartpar, jendpar, nbl_exchg, nboundlines, my_cart_neigh, &
               lperi_x, lperi_y, l2dim, &
               10000+ntstep, .FALSE., ncomm_type, izerror, yerrmsg,            &
               tke, tketens )
        ELSE
          IF ( l3dturb ) THEN   ! Not sure if necessary in all cases or just in 3D-Turb.
            kzdims(1:24)=(/ke,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
            CALL exchg_boundaries                                                  &
                 (0, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute, ie, je,       &
                 kzdims, jstartpar, jendpar, nbl_exchg, nboundlines, my_cart_neigh, &
                 lperi_x, lperi_y, l2dim, &
                 10000+ntstep, .FALSE., ncomm_type, izerror, yerrmsg,            &
                 tke )
          END IF
        END IF
      END SELECT

    END IF
! ub<<
      
    IF (nradcoarse > 1) THEN
      ! radiation computations on a coarser grid need additional
      ! boundary exchange

      IF ( lsoil ) THEN
        IF(lmulti_layer) THEN
          kzdims(1:24)=(/1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
         CALL exchg_boundaries                                               &
           (70+nnew, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute, ie, je,    &
            kzdims, jstartpar, jendpar, nradcoarse-1, nboundlines, my_cart_neigh,  &
            lperi_x, lperi_y, l2dim, &
            20000+ntstep, ldatatypes, ncomm_type, izerror, yerrmsg,          &
            t_g(:,:,nnew),w_so(:,:,1,nnew),t_snow(:,:,nnew),w_snow(:,:,nnew),&
            freshsnow(:,:) )
        ELSE
          kzdims(1:24)=(/1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
          CALL exchg_boundaries                                              &
            (73+nnew, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute, ie, je,   &
            kzdims, jstartpar, jendpar, nradcoarse-1, nboundlines, my_cart_neigh,  &
            lperi_x, lperi_y, l2dim, &
            20000+ntstep, ldatatypes, ncomm_type, izerror, yerrmsg,          &
            t_g(:,:,nnew),w_g1(:,:,nnew),t_snow(:,:,nnew),w_snow(:,:,nnew))
        ENDIF
      ELSE ! .NOT.lsoil:
          kzdims(1:24)=(/1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
          CALL exchg_boundaries  & !exchange of t_snow probably not necessary
            (76+nnew, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute, ie, je,   &
            kzdims, jstartpar, jendpar, nradcoarse-1, nboundlines, my_cart_neigh,  &
            lperi_x, lperi_y, l2dim, &
            20000+ntstep, ldatatypes, ncomm_type, izerror, yerrmsg,          &
            t_g(:,:,nnew),t_snow(:,:,nnew) )
      ENDIF

      IF (lzconv) THEN
        kzdims(1:24)=(/ke,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
        CALL exchg_boundaries                                               &
          (2, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute, ie, je,          &
          kzdims, jstartpar, jendpar, nradcoarse-1, nboundlines, my_cart_neigh,  &
          lperi_x, lperi_y, l2dim, &
          20000+ntstep, ldatatypes, ncomm_type, izerror, yerrmsg,           &
          clc_con(:,:,:))
      ENDIF

    ENDIF   ! nradcoarse > 1

    IF (ltime) CALL get_timings (i_communications_phy, ntstep, dt, izerror)
!!$  ENDIF

  IF (lzconv) THEN
    ! convective tendencies for u and v have been exchanged above and
    ! have to be interpolated to the u- and v-grid points
    DO  k = 1, ke
      DO    j = jstartu, jendu
        DO  i = istartu, iendu
          ut_conv(i,j,k) = 0.5*( ut_conv(i+1,j,k) + ut_conv(i,j,k) )
        ENDDO
      ENDDO
      DO    j = jstartv, jendv
        DO  i = istartv, iendv
          vt_conv(i,j,k) = 0.5*( vt_conv(i,j+1,k) + vt_conv(i,j,k) )
        ENDDO
      ENDDO
    ENDDO
  ENDIF

  ! SSO tendencies for u and v have been exchanged above and
  ! have to be interpolated to the u- and v-grid points
  IF (lsso) THEN
    DO  k = 1, ke
      DO    j = jstartu, jendu
        DO  i = istartu, iendu
          ut_sso(i,j,k) = 0.5*( ut_sso(i+1,j,k) + ut_sso(i,j,k) )
        ENDDO
      ENDDO
      DO    j = jstartv, jendv
        DO  i = istartv, iendv
          vt_sso(i,j,k) = 0.5*( vt_sso(i,j+1,k) + vt_sso(i,j,k) )
        ENDDO
      ENDDO
    ENDDO
  ENDIF

  IF (izdebug > 5) THEN
    PRINT *, '      PHYSICAL COMPUTATIONS before dynamics finished '
  ENDIF

!------------------------------------------------------------------------------
! Section 4: Physical packages at the end of the time stepping
!------------------------------------------------------------------------------

ELSEIF (yaction == 'finish_compute') THEN

  IF (lgsp .AND. .NOT.(l2tls .AND. irunge_kutta == 0))  THEN
    IF (izdebug > 5) THEN
      PRINT *, '      GRID SCALE PRECIPITATION after dynamics'
    ENDIF

    IF (lprogprec) THEN
      IF     (itype_gscp == 4) THEN
        CALL hydci_pp_gr
      ELSEIF (itype_gscp == 3) THEN
        CALL hydci_pp
      ELSEIF (itype_gscp == 2) THEN
        CALL hydor_pp
      ELSEIF (itype_gscp == 1) THEN
        CALL kessler_pp
      ENDIF
    ELSE
      IF (lprog_qi) CALL hydci
    ENDIF

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
  REAL (KIND=ireals)         ::       &
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

    icldm_rad_d,  & ! mode of cloud representation in radiation  parametr.
    icldm_tran_d, & ! mode of cloud representation in transfer  parametr.
    icldm_turb_d, & ! mode of cloud representation in turbulence parametr.
    itype_conv_d, & ! type of convection parameterization

    itype_aerosol_d, &! type of aerosol map
    itype_root_d,   & ! type of root density distribution
    itype_heatcond_d,&! type of heat conductivity
    itype_hydbound_d,&! type of hydraulic lower boundary

    nhori_d,      & ! number of sectors for the horizont by the topographic radcorr
    ke_soil_d,    & ! number of layers in multi-layer soil model
    ke_snow_d,    & ! number of layers in multi-layer snow model
    nlgw_d          ! number of prognostic soil water levels

  LOGICAL                    ::       &
    lrad_d,       & ! forecast with radiation
    lforest_d,    & ! to run with forest data (evergreen and deciduous)
    ltur_d,       & ! forecast with vertical diffusion
    lconv_d,      & ! forecast with convection
    lconv_inst_d, & ! output of instantaneous values of top_con/bas_con
                    ! instead of min/max for an output interval
    lgsp_d,       & ! forecast with grid scale precipitation
    lprogprec_d,  & ! forecast with prognostic rain and snow (qr, qs)
    ltrans_prec_d,& ! forecast with prognostic rain and snow (qr, qs)
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
    lradtopo_d      ! uses topographic correction of radiation

  INTEGER (KIND=iintegers)   :: i, k, invar, ierr, iz_err
  LOGICAL                    :: lzequiv

! Define the namelist group

  NAMELIST /phyctl/ lrad, ltur, lconv, itype_conv, lgsp,                      &
                    lsoil, lmelt, lmelt_var, lmulti_layer, lexpcor,           &
                    ltmpcor, lprfcor, lnonloc, lcpfluc,lcape,lctke, lconf_avg,&
                    nincrad, hincrad, ninctura, nincconv,                     &
                    lprogprec, ltrans_prec, ldiniprec, itype_trvg, itype_evsl,&
                    itype_gscp, itype_wcld, itype_tran, itype_turb,itype_synd,&
                    icldm_rad, icldm_tran, icldm_turb, imode_tran, imode_turb,&
                    ke_soil, czml_soil, nlgw, l3dturb, lprog_tke,             &
                    lforest, lconv_inst, lseaice, llake, ico2_rad, czbot_w_so,&
                    l3dturb_metr, nradcoarse, lradf_avg, lradtopo, nhori,     &
                    lsso, nincsso, limpltkediff, ltkesso, ltkecon,            &
                    itype_sher, lmulti_snow, ke_snow, lemiss, lstomata,       &
                    itype_aerosol, itype_root, itype_heatcond, itype_hydbound

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
  lprogprec_d    = .TRUE.
  ltrans_prec_d  = .TRUE.
  ldiniprec_d    = .FALSE.
  l3dturb_d      = .FALSE.
  l3dturb_metr_d = .TRUE.
  lprog_tke_d    = .FALSE.
  limpltkediff_d = .FALSE.
  lsoil_d        = .TRUE.
  lmelt_d        = .TRUE.
  lmelt_var_d    = .TRUE.
  lmulti_layer_d = .TRUE.
  lmulti_snow_d  = .FALSE.
  lseaice_d      = .FALSE.
  llake_d        = .FALSE.
  lemiss_d       = .FALSE.
  lstomata_d     = .FALSE.

  lsso_d         = .FALSE.
  ltkesso_d      = .FALSE.
  ltkecon_d      = .FALSE.
  lexpcor_d      = .TRUE.
  ltmpcor_d      = .FALSE.
  lprfcor_d      = .FALSE.
  lnonloc_d      = .FALSE.
  lcpfluc_d      = .FALSE.
  lconf_avg_d    = .TRUE. 
  lradf_avg_d    = .FALSE.
  lcape_d        = .FALSE.
  lctke_d        = .FALSE.
  lradtopo_d     = .FALSE.

  nhori_d      = 24
  nradcoarse_d = 1
  nincrad_d    = 0
  hincrad_d    = 1.0_ireals
  ninctura_d   = 1
  nincconv_d   = 10
  nincsso_d    = 5

  itype_trvg_d = 2 ! vegetation transpiration using the BATS-approach
  itype_evsl_d = 2 ! bare soil evaporation using the BATS-approach

  itype_gscp_d = 3 ! grid-scale-cloud-precipitation with 'hydor'

  itype_sher_d = 1 ! only vertical shear
  itype_wcld_d = 2 ! water cloud diagnosis (in SUB cloud_diag) using statistical condensation
  itype_tran_d = 2 ! new surface-atmosphere transfer using SUB turb_tran
  itype_turb_d = 3 ! new moist scheme with prognostic tke-equation using SUB turb_diff
  itype_synd_d = 2 ! new method with resistance formulation using SUB turb_tran

  imode_tran_d = 1 ! diagnostik TKE in 'turbtran'
  imode_turb_d = 1 ! implicite calculation of turbulent diffusion using the explicit
                   ! surface flux densities as lower boundary condition

  ico2_rad_d   = 0 ! constant CO2 concentration (330 ppm)

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
  czml_soil_d(1)  = 0.005_ireals
  czml_soil_d(2)  = 0.02_ireals
  czml_soil_d(3)  = 0.06_ireals
  czml_soil_d(4)  = 0.18_ireals
  czml_soil_d(5)  = 0.54_ireals
  czml_soil_d(6)  = 1.62_ireals
  czml_soil_d(7)  = 4.86_ireals
  czml_soil_d(8)  =14.58_ireals
  czbot_w_so_d    = 2.5_ireals

  ke_snow_d      = 2
  nlgw_d         = 2

!------------------------------------------------------------------------------
!- Section 2: Initialize variables with defaults
!------------------------------------------------------------------------------

  lrad           = lrad_d
  lforest        = lforest_d
  ltur           = ltur_d
  lconv          = lconv_d
  lconv_inst     = lconv_inst_d
  lgsp           = lgsp_d
  lprogprec      = lprogprec_d
  ltrans_prec    = ltrans_prec_d
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
  czml_soil(:)   = -1.0_ireals
  czbot_w_so     = czbot_w_so_d
  nlgw           = nlgw_d

!------------------------------------------------------------------------------
!- Section 3: Input of the namelist values
!------------------------------------------------------------------------------

  READ (nuin, phyctl, IOSTAT=iz_err)
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
      IF (hincrad <= 0.0_ireals) THEN
        PRINT *, ' ERROR  *** Wrong value for hincrad: ', hincrad, ' *** '
        PRINT *, '        ***   must be > 0.0 *** '
        ierrstat = 1002
      ENDIF
      nincrad = NINT ( hincrad * 3600.0_ireals/dt)
    ELSE
      hincrad = 0.0_ireals
    ENDIF
  ELSE  
    IF ( hincrad /= hincrad_d ) THEN
      IF (hincrad <= 0.0_ireals) THEN
        PRINT *, ' ERROR  *** Wrong value for hincrad: ', hincrad, ' *** '
        PRINT *, '        ***   must be > 0.0 *** '
        ierrstat = 1002
      ENDIF
    ENDIF
    nincrad = NINT ( hincrad * 3600.0_ireals/dt)
  ENDIF
  hnextrad = hstart
  nextrad  = NINT ( hnextrad * 3600.0_ireals / dt)


  IF ( (nradcoarse > nboundlines+1) .OR. (nradcoarse < 1) ) THEN
    PRINT *, ' ERROR  *** Wrong value for nradcoarse: ', nradcoarse, ' *** '
    PRINT *, '        ***   must be >= 1 and <= nboundlines+1 !*** '
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

  IF ( (itype_gscp < 1) .OR. (itype_gscp > 4) ) THEN
    itype_gscp = itype_gscp_d
    PRINT *, ' ERROR  *** Wrong value for itype_gscp: ', itype_gscp, ' *** '
    PRINT *, '        ***   must be >= 1 and <= 4 !!!                  *** '
    ierrstat = 1002
  ENDIF
  IF ( (.NOT. lgsp) .AND. (itype_gscp /= itype_gscp_d) ) THEN
    itype_gscp = itype_gscp_d
    PRINT *, ' WARNING:  *** itype_gscp is set to the default *** ',    &
             '           *** because precipitation is turned off ***'
  ENDIF
  IF ( (.NOT. lgsp) .AND. (lprogprec) ) THEN
    lprogprec = .FALSE.
    PRINT *, ' WARNING:  *** prognostic precipitation is turned off *** '
    PRINT *, '           *** because no precipitation is computed *** '
  ENDIF
  IF ( (itype_gscp == 1) .AND. l2tls .AND. lprogprec .AND. (irunge_kutta == 0)) THEN
    PRINT  *, ' ERROR  *** Kessler scheme with prognostic precipitation'
    PRINT  *, 'only implemented for leapfrog or 2 timelevel RK *** '
    ierrstat = 1002
  ENDIF
  IF ( (itype_gscp == 1) .AND. ldiniprec) THEN
    PRINT  *, ' ERROR  *** Initialization of prognostic precipitation ***'
    PRINT  *, '        *** only possible for itype_gscp > 1           *** '
    ierrstat = 1002
  ENDIF
  IF (.NOT. ltrans_prec .AND. l2tls .AND.             &
        .NOT. ( (y_scalar_advect=="SL3_MF") .OR.      &
                (y_scalar_advect=="SL3_SFD")     ) ) THEN
    PRINT  *, ' ERROR  *** Transport of Precipitation cannot be switched ***'
    PRINT  *, '        ***   off for the 2 timelevel Runge-Kutta scheme  ***'
    ierrstat = 1002
  ENDIF
  IF ( (itype_gscp == 4) .AND. .NOT.lprogprec)  THEN
    PRINT  *, ' ERROR  *** Graupel scheme only with ',          &
                         'prognostic precipitation *** '
    ierrstat = 1002
  ENDIF
  IF ( (itype_gscp == 4) .AND. l2tls .AND. (irunge_kutta == 0)) THEN
    PRINT  *, ' ERROR  *** Graupel scheme only implemented for ',          &
                         'leapfrog or 2 timelevel RK *** '
    ierrstat = 1002
  ENDIF
  IF ( (itype_gscp < 3) .AND. luse_rttov) THEN
    PRINT  *, ' ERROR  *** Without cloud ice, satellite images ',          &
                                           'cannot be computed ***'
    ierrstat = 1002
  ENDIF
  IF ( (.NOT. lprogprec) .AND. (ltrans_prec) ) THEN
    ltrans_prec=.FALSE.
    PRINT *,' WARNING:  *** transport of precipitation has been switched  *** '
    PRINT *,'           *** off because precipitation is not prognostic   *** '
  END IF
  IF ( (.NOT. lprogprec) .AND. (ldiniprec) ) THEN
    ldiniprec = .FALSE.
    PRINT *,' WARNING:  *** diagnostic initialisation of precip is turned *** '
    PRINT *,'           *** off because precipitation is not prognostic   *** '
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

  SELECT CASE( itype_turb )
  CASE( 1, 3 )
    lprog_tke  = .FALSE.
  CASE( 5:8 )
    IF ( .NOT.( l2tls .AND. irunge_kutta >= 1) ) THEN
      lprog_tke  = .FALSE.
      itype_turb = itype_turb_d
      PRINT *,' WARNING: *** itype_turb is set to default again *** '
    END IF
  CASE( 100 )
    lprog_tke  = .FALSE.
  CASE default
    lprog_tke  = .FALSE.
    itype_turb = itype_turb_d
    PRINT *,' WARNING: *** itype_turb is set to default again *** '
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

  IF (ico2_rad.GT.6 )  ico2_rad  = ico2_rad_d

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
      PRINT *,' WARNING: *** ltkeconv cannot be active, since convection is not running *** '
    ELSEIF (itype_conv == 2 .OR. itype_conv == 3) THEN
      PRINT *,' WARNING: *** ltkeconv is not yet supported by convection according Kain-Frisch or Bechthold *** '
    END IF
  END IF

  IF (lmulti_layer) THEN
    ALLOCATE (czmls(1:ke_soil+1), STAT=ierr)
    ALLOCATE (czhls(0:ke_soil+1), STAT=ierr)
    ALLOCATE (msoilgrib(0:ke_soil+1), STAT=ierr)
                    ! (careful: the level k=1 will be coded with 1,
                    !           but is in the depth of 0.5 cm!)

    ! Check, how many soil levels have been specified
    invar = COUNT(czml_soil(:) /= -1.0_ireals)
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
        msoilgrib(i) = NINT (100 * czmls(i)+1.0E-7)
      ENDDO
    ENDIF

    ! determine depth of half levels out of czmls
    czhls(0) = 0.0_ireals
    DO k = 1, ke_soil+1
      czhls(k) = 2.0_ireals * czmls(k) - czhls(k-1)
    ENDDO

    czbot_w_so = MIN(czbot_w_so, czhls(ke_soil))
    ibot_w_so = ke_soil
    DO i=1,ke_soil+1
      IF (czhls(i) <= czbot_w_so) ibot_w_so=i
    ENDDO
  ENDIF

  IF (lconv) THEN

#ifndef HYMACS
    IF ( (itype_conv < 0) .OR. (itype_conv > 3) ) THEN
      PRINT  *, ' ERROR  *** Wrong type of convection scheme: ', itype_conv
      PRINT  *, '        *** Must be between 0 and 3!'
      ierrstat = 1002
    ENDIF
#endif
#ifdef HYMACS
! VK 2012/03/15
    IF ( (itype_conv < 0) .OR. (itype_conv > 4) ) THEN
      PRINT  *, ' ERROR  *** Wrong type of convection scheme: ', itype_conv
      PRINT  *, '        *** Must be between 0 and 4!'
      ierrstat = 1002
    ENDIF
#endif

     IF (.NOT.ltur) THEN
       lctke=.FALSE.
    END IF
  ENDIF

  IF ( (itype_aerosol < 1) .OR. (itype_aerosol > 2) ) THEN
    PRINT  *, ' ERROR  *** Wrong type of aerosol scheme: ', itype_aerosol
    PRINT  *, '        *** Must be between 1 and 2!'
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
    logbuf (18) = lprogprec
    logbuf (19) = ltrans_prec
    logbuf (20) = ldiniprec
    logbuf (21) = l3dturb
    logbuf (22) = lprog_tke
    logbuf (23) = limpltkediff
    logbuf (24) = lforest
    logbuf (25) = lconv_inst
    logbuf (26) = llake
    logbuf (27) = l3dturb_metr
    logbuf (28) = lradf_avg
    logbuf (29) = lradtopo
    logbuf (30) = lsso
    logbuf (31) = ltkesso
    logbuf (32) = lseaice
    logbuf (33) = lmulti_snow
    logbuf (34) = lemiss
    logbuf (35) = lstomata

    realbuf( 1) = czbot_w_so
    realbuf( 2) = hincrad
    realbuf( 3) = hnextrad

    IF (lmulti_layer) THEN
      DO i = 1, ke_soil+1
        realbuf(3+i) = czmls(i)
      ENDDO
    ELSE
      realbuf( 4:40) = 0.0_ireals
    ENDIF
    
  ENDIF

  CALL distribute_values (intbuf, 30, 0, imp_integers, icomm_world, ierr)
  CALL distribute_values (logbuf, 35, 0, imp_logical,  icomm_world, ierr)
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
    lprogprec    = logbuf (18)
    ltrans_prec  = logbuf (19)
    ldiniprec    = logbuf (20)
    l3dturb      = logbuf (21)
    lprog_tke    = logbuf (22)
    limpltkediff = logbuf (23)
    lforest      = logbuf (24)
    lconv_inst   = logbuf (25)
    llake        = logbuf (26)
    l3dturb_metr = logbuf (27)
    lradf_avg    = logbuf (28)
    lradtopo     = logbuf (29)
    lsso         = logbuf (30)
    ltkesso      = logbuf (31)
    lseaice      = logbuf (32)
    lmulti_snow  = logbuf (33)
    lemiss       = logbuf (34)
    lstomata     = logbuf (35)

    czbot_w_so   = realbuf( 1)
    hincrad      = realbuf( 2)
    hnextrad     = realbuf( 3)

    IF (lmulti_layer) THEN
      ALLOCATE (czmls(1:ke_soil+1), STAT=ierr)
      ALLOCATE (msoilgrib(0:ke_soil+1), STAT=ierr)
      msoilgrib(0) = 0_iintegers
      DO i = 1, ke_soil+1
        czmls(i) = realbuf(3+i)
        msoilgrib(i) = NINT (100 * czmls(i)+1.0E-7)
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
  WRITE (nuspecif, '(T7,A,T21,A,T39,A,T58,A)') 'Variable', 'Actual Value',   &
                                               'Default Value', 'Format'

  WRITE (nuspecif, '(T8,A,T21,L12  ,T40,L12  ,T59,A3)')                      &
                                           'lgsp'   ,lgsp   ,lgsp_d,   ' L '
  WRITE (nuspecif, '(T8,A,T21,L12  ,T40,L12  ,T59,A3)')                      &
                            'lprogprec'   ,lprogprec   ,lprogprec_d,   ' L '
  WRITE (nuspecif, '(T8,A,T21,L12  ,T40,L12  ,T59,A3)')                      &
                      'ltrans_prec'   ,ltrans_prec   ,ltrans_prec_d,   ' L '
  WRITE (nuspecif, '(T8,A,T21,L12  ,T40,L12  ,T59,A3)')                      &
                            'ldiniprec'   ,ldiniprec   ,ldiniprec_d,   ' L '
  WRITE (nuspecif, '(T8,A,T21,L12  ,T40,L12  ,T59,A3)')                      &
                                           'lrad',   lrad,   lrad_d,   ' L '
  WRITE (nuspecif, '(T8,A,T21,I12  ,T40,I12  ,T59,A3)')                      &
                       'itype_aerosol', itype_aerosol, itype_aerosol_d,' I '
  WRITE (nuspecif, '(T8,A,T21,L12  ,T40,L12  ,T59,A3)')                      &
                                     'lemiss   ',lemiss   ,lemiss_d   ,' L '
  WRITE (nuspecif, '(T8,A,T21,L12  ,T40,L12  ,T59,A3)')                      &
                                           'lforest',lforest,lforest_d,' L '
  WRITE (nuspecif, '(T8,A,T21,L12  ,T40,L12  ,T59,A3)')                      &
                                           'ltur   ',ltur   ,ltur_d   ,' L '
  WRITE (nuspecif, '(T8,A,T21,L12  ,T40,L12  ,T59,A3)')                      &
                                      'l3dturb', l3dturb, l3dturb_d,   ' L '
  WRITE (nuspecif, '(T8,A,T21,L12  ,T40,L12  ,T59,A3)')                      &
                       'l3dturb_metr', l3dturb_metr, l3dturb_metr_d,   ' L '
  WRITE (nuspecif, '(T8,A,T21,L12  ,T40,L12  ,T59,A3)')                      &
                                'lprog_tke', lprog_tke, lprog_tke_d,   ' L '
  WRITE (nuspecif, '(T8,A,T21,L12  ,T40,L12  ,T59,A3)')                      &
                       'limpltkediff', limpltkediff, limpltkediff_d,   ' L '
  WRITE (nuspecif, '(T8,A,T21,L12  ,T40,L12  ,T59,A3)')                      &
                                           'lconv  ',lconv  ,lconv_d  ,' L '
  WRITE (nuspecif, '(T8,A,T21,I12  ,T40,I12  ,T59,A3)')                      &
                                'itype_conv', itype_conv, itype_conv_d,' I '
  WRITE (nuspecif, '(T8,A,T21,L12  ,T40,L12  ,T59,A3)')                      &
                                  'lconv_inst',lconv_inst,lconv_inst_d,' L '
  WRITE (nuspecif, '(T8,A,T21,L12  ,T40,L12  ,T59,A3)')                      &
                                           'lsoil  ',lsoil  ,lsoil_d  ,' L '
  WRITE (nuspecif, '(T8,A,T21,L12  ,T40,L12  ,T59,A3)')                      &
                                           'lseaice',lseaice,lseaice_d,' L '
  WRITE (nuspecif, '(T8,A,T21,L12  ,T40,L12  ,T59,A3)')                      &
                                           'llake  ',llake  ,llake_d  ,' L '
  WRITE (nuspecif, '(T8,A,T21,L12  ,T40,L12  ,T59,A3)')                      &
                                           'lsso   ',lsso   ,lsso_d   ,' L '
  WRITE (nuspecif, '(T8,A,T21,L12  ,T40,L12  ,T59,A3)')                      &
                                           'lmelt  ',lmelt  ,lmelt_d  ,' L '
  WRITE (nuspecif, '(T8,A,T21,L12  ,T40,L12  ,T59,A3)')                      &
                                 'lmelt_var',lmelt_var  ,lmelt_var_d  ,' L '
  WRITE (nuspecif, '(T8,A,T21,L12  ,T40,L12  ,T59,A3)')                      &
                         'lmulti_layer ',lmulti_layer,lmulti_layer_d  ,' L '

  IF (lmulti_layer) THEN
    ! Write specification for soil levels
    WRITE (nuspecif, '(T8,A,T21,I12  ,T40,I12  ,T59,A3)')                    &
                                        'ke_soil ',ke_soil ,ke_soil_d, ' I '
    WRITE (nuspecif, '(T10,A)') 'Main levels of the soil layers'
    WRITE (nuspecif, '(T10,A)') '                  (m)           (cm)'
    WRITE (nuspecif, '(A,I12)') '              0:               ', msoilgrib(0)
    DO i = 1, ke_soil+1
      WRITE (nuspecif, '(I15,A,F12.4,I12)') i, ':   ',czmls(i), msoilgrib(i)
    ENDDO
    WRITE (nuspecif, '(T8,A,T21,F12.4,T40,F12.4,T59,A3)')                    &
                                  'czbot_w_so',czbot_w_so,czbot_w_so_d,' R '
  ENDIF

  WRITE (nuspecif, '(T8,A,T21,L12  ,T40,L12  ,T59,A3)')                      &
                         'lmulti_snow  ',lmulti_snow,  lmulti_snow_d  ,' L '
  WRITE (nuspecif, '(T8,A,T21,I12  ,T40,I12  ,T59,A3)')                      &
                                        'ke_snow ',ke_snow ,ke_snow_d, ' I '
  WRITE (nuspecif, '(T8,A,T27,I6   ,T40,I12  ,T59,A3)')                      &
              'itype_heatcond',itype_heatcond  ,itype_heatcond_d   ,   ' I '
  WRITE (nuspecif, '(T8,A,T27,I6   ,T40,I12  ,T59,A3)')                      &
              'itype_hydbound',itype_hydbound  ,itype_hydbound_d   ,   ' I '
  WRITE (nuspecif, '(T8,A,T21,I12  ,T40,I12  ,T59,A3)')                      &
                          'itype_root',itype_root  ,itype_root_d   ,   ' I '
  WRITE (nuspecif, '(T8,A,T21,L12  ,T40,L12  ,T59,A3)')                      &
                                'lstomata   ',lstomata   ,lstomata_d  ,' L '
  WRITE (nuspecif, '(T8,A,T21,L12  ,T40,L12  ,T59,A3)')                      &
                                           'ltkesso',ltkesso,ltkesso_d,' L '
  WRITE (nuspecif, '(T8,A,T21,L12  ,T40,L12  ,T59,A3)')                      &
                                           'ltkecon',ltkecon,ltkecon_d,' L '
  WRITE (nuspecif, '(T8,A,T21,L12  ,T40,L12  ,T59,A3)')                      &
                                           'lexpcor',lexpcor,lexpcor_d,' L '
  WRITE (nuspecif, '(T8,A,T21,L12  ,T40,L12  ,T59,A3)')                      &
                                           'ltmpcor',ltmpcor,ltmpcor_d,' L '
  WRITE (nuspecif, '(T8,A,T21,L12  ,T40,L12  ,T59,A3)')                      &
                                           'lprfcor',lprfcor,lprfcor_d,' L '
  WRITE (nuspecif, '(T8,A,T21,L12  ,T40,L12  ,T59,A3)')                      &
                                           'lnonloc',lnonloc,lnonloc_d,' L '
  WRITE (nuspecif, '(T8,A,T21,L12  ,T40,L12  ,T59,A3)')                      &
                                           'lcpfluc',lcpfluc,lcpfluc_d,' L '
  WRITE (nuspecif, '(T8,A,T21,L12  ,T40,L12  ,T59,A3)')                      &
                                   'lconf_avg',lconf_avg,lconf_avg_d  ,' L '
  WRITE (nuspecif, '(T8,A,T21,L12  ,T40,L12  ,T59,A3)')                      &
                                   'lradf_avg',lradf_avg,lradf_avg_d  ,' L '
  WRITE (nuspecif, '(T8,A,T21,L12  ,T40,L12  ,T59,A3)')                      &
                                           'lcape  ',lcape  ,lcape_d  ,' L '
  WRITE (nuspecif, '(T8,A,T21,L12  ,T40,L12  ,T59,A3)')                      &
                                           'lctke  ',lctke  ,lctke_d  ,' L '

  WRITE (nuspecif, '(T8,A,T21,F12.4,T40,F12.4,T59,A3)')                      &
                                           'hincrad',hincrad,hincrad_d,' R '
  WRITE (nuspecif, '(T8,A,T21,I12  ,T40,I12  ,T59,A3)')                      &
                                           'nincrad',nincrad,nincrad_d,' I '
  WRITE (nuspecif, '(T8,A,T21,I12  ,T40,I12  ,T59,A3)')                      &
                                           'nradcoarse',nradcoarse,nradcoarse_d,' I '
  WRITE (nuspecif, '(T8,A,T21,I12  ,T40,I12  ,T59,A3)')                      &
                                        'ninctura',ninctura,ninctura_d,' I '
  WRITE (nuspecif, '(T8,A,T21,I12  ,T40,I12  ,T59,A3)')                      &
                                        'nincconv',nincconv,nincconv_d,' I '
  WRITE (nuspecif, '(T8,A,T21,I12  ,T40,I12  ,T59,A3)')                      &
                                           'nincsso',nincsso,nincsso_d,' I '

  WRITE (nuspecif, '(T8,A,T21,I12  ,T40,I12  ,T59,A3)')                      &
                                    'itype_trvg',itype_trvg,itype_trvg_d,' I '
  WRITE (nuspecif, '(T8,A,T21,I12  ,T40,I12  ,T59,A3)')                      &
                                    'itype_evsl',itype_evsl,itype_evsl_d,' I '

  WRITE (nuspecif, '(T8,A,T21,I12  ,T40,I12  ,T59,A3)')                      &
                                    'itype_gscp',itype_gscp,itype_gscp_d,' I '

  WRITE (nuspecif, '(T8,A,T21,I12  ,T40,I12  ,T59,A3)')                      &
                                    'itype_sher',itype_sher,itype_sher_d,' I '
  WRITE (nuspecif, '(T8,A,T21,I12  ,T40,I12  ,T59,A3)')                      &
                                    'itype_wcld',itype_wcld,itype_wcld_d,' I '
  WRITE (nuspecif, '(T8,A,T21,I12  ,T40,I12  ,T59,A3)')                      &
                                    'itype_tran',itype_tran,itype_tran_d,' I '
  WRITE (nuspecif, '(T8,A,T21,I12  ,T40,I12  ,T59,A3)')                      &
                                    'itype_turb',itype_turb,itype_turb_d,' I '
  WRITE (nuspecif, '(T8,A,T21,I12  ,T40,I12  ,T59,A3)')                      &
                                    'itype_synd',itype_synd,itype_synd_d,' I '

  WRITE (nuspecif, '(T8,A,T21,I12  ,T40,I12  ,T59,A3)')                      &
                                    'imode_tran',imode_tran,imode_tran_d,' I '
  WRITE (nuspecif, '(T8,A,T21,I12  ,T40,I12  ,T59,A3)')                      &
                                    'imode_turb',imode_turb,imode_turb_d,' I '
  WRITE (nuspecif, '(T8,A,T21,I12  ,T40,I12  ,T59,A3)')                      &
                                    'ico2_rad',ico2_rad,ico2_rad_d,      ' I '

  WRITE (nuspecif, '(T8,A,T21,I12  ,T40,I12  ,T59,A3)')                      &
                                    'icldm_rad',icldm_rad,icldm_rad_d,   ' I '
  WRITE (nuspecif, '(T8,A,T21,I12  ,T40,I12  ,T59,A3)')                      &
                                    'icldm_tran',icldm_tran,icldm_tran_d,' I '
  WRITE (nuspecif, '(T8,A,T21,I12  ,T40,I12  ,T59,A3)')                      &
                                    'icldm_turb',icldm_turb,icldm_turb_d,' I '

  WRITE (nuspecif, '(T8,A,T21,I12  ,T40,I12  ,T59,A3)')                      &
                                        'nlgw    ',nlgw    ,nlgw_d,    ' I '
  WRITE (nuspecif, '(T8,A,T21,L12  ,T40,L12  ,T59,A3)')                      &
                                      'lradtopo',lradtopo,lradtopo_d  ,' L '
  WRITE (nuspecif, '(T8,A,T21,I12  ,T40,I12  ,T59,A3)')                      &
                                        'nhori', nhori, nhori_d,       ' I '
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
