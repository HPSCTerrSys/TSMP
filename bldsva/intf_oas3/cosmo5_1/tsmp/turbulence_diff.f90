!+ Source module for computing diffusion coefficients
!------------------------------------------------------------------------------

 MODULE turbulence_diff

!------------------------------------------------------------------------------
!
! Description:
!
! The module turbulence_diff calculates the tendencies for turbulent
! vertial transport of momentum and heat and the coefficients
! for turbulent diffusion as well.
!
! The clousure is made on lever 2.5 (Mellor/Yamada) using a prognostik
! TKE-equation and includes the formulation of a flow through a porous 
! medium (roughnes layer)
!
! The turbulence model (with some Prandtl-layer approximations is used 
! for the calculation of turbulent transfer between atmosphere and the
! lower boundary too.
!
!------------------------------------------------------------------------------
!
! Current Code Owner: DWD, Matthias Raschendorfer
!  phone:  +49  69  8062 2708
!  fax:    +49  69  8062 3721
!  email:  matthias.raschendorfer@dwd.de
!
! History:
! Version    Date       Name
! ---------- ---------- ----
! 1.30       1999/06/24 Matthias Raschendorfer
!  Initial release
! 1.33       1999/10/14 Matthias Raschendorfer
!  Introduction of the LOGICAL parameters ltmpcor, lstfnct.
!  Rearranging the code in order to make it run faster.
!  Reformulation of the Charnock-Formula using the additional parameter alpha1.
!  Using the redefined parameter c_g (including the factor 3).
! 1.34       1999/12/10 Matthias Raschendofer
!  Introduction of minimal Diffusion Coefficients.
!  Modification in the Formulation of subgrid scale condensation.
!  Consideration of partial cloudiness on the TKE production due to thermal circulations.
! 1.39       2000/05/03 Ulrich Schaettler
!  Changed some variable names.
! 2.2        2000/08/18 Matthias Raschendorfer
!  tkv(h,m) limited by the molecular diff.coef. con_(h,m) for small values.
! 2.3        2000/11/15 Guenther Doms
!  Some local variable have been redefined to allow for reproducible
!  results on MPP platforms.
! 2.15       2002/03/20 Matthias Raschendorfer
!  Some formal modifications to make the code more efficient on vector machines
! 2.17       2002/05/08 Ulrich Schaettler
!  Some more optimizations for vector machines
! 2.19       2002/10/24 Ulrich Schaettler
!  Adaptations to 2-timelevel scheme (use of ntlev)
! 3.5        2003/09/02 Ulrich Schaettler
!  Adaptation of the interface for exchg_boundaries
! 3.6        2003/12/11 Ulrich Schaettler
!  Only editorial changes
! 3.7        2004/02/18 Ulrich Schaettler
!  Increased dimension of kzdims
! 3.13       2004/12/03 Thorsten Reinhardt
!  Replaced SIGN-Function by IF-statements
! 3.16       2005/07/22 Matthias Raschendorfer
!  Introduction of the smoothing parameter 'tsmot'; removal of the parameter 'q_max'.
!  Changing the restriction in order to avoid a singularity in the stability function.
!  Modification of vertical diffusion in the TKE-equation
! 3.18       2006/03/03 Ulrich Schaettler
!  Eliminated call to exchg_boundaries for wind-tendencies, which is
!  not necessary here
! 3.19       2006/04/25 Matthias Raschendorfer / Jochen Foerstner
!  Correction of a bug related to the Exner-factor in the explicit correction
!  and some formal modifications
! 3.21       2006/12/04 Matthias Raschendorfer, Dmitrii Mironov, Ulrich Schaettler
!  Numerical security check for the second term of explicit vertical TKE diffusion
!  Introduction of the total dry option "icldm_turb.EQ.-1".
!  Changes to use the FLake model.
!  Changed interface to cloud_diag from meteo_utilities.
! V3_23        2007/03/30 Matthias Raschendorfer
!  Introduction of some output variables for SCLM.
!  Change of some parameter names.
!  Eliminated stab_funct.incf and turb_param.incf
!  (Substitution of file 'stab_funct.incf' by a SUBROUTINE)
!  New file 'statement_functs.incf' containing statement functions
! V3_24        2007/04/26 Ulrich Schaettler
!  Introduced call to exchg_boundaries for wind-tendencies again;
!  it is necessary for imode_turb=2/3!
! V4_3         2008/02/25 Matthias Raschendorfer
!  Calculation of the 3D diagnostic field 'edr' for the eddy dissipotion rate.
!  Changing the treatment of parameter field dd(:,:,:) to achiev better vectorisation
!  in SUBROUTINE 'stab_funct'.
! V4_9         2009/07/16 Ulrich Schaettler, Christian Bollmann
!   Eliminated small loops to get the compiler vectorize over the correct loop
! V4_10        2009/09/11 Matthias Raschendorfer, Jan-Peter Schulz
!  Introduction of implicit vertical diffusion also for TKE and 
!   correction of a bug related with the explicite TKE diffusion (by Oliver Fuhrer).
!  Removing the horizontal loops for SCLM treatment.
!  Adoptions with respect to the implicit vertical TKE diffusion.  
!  Introduction of 3D and horizontal corrections for windshear production incl. metric terms.
!  Introduction of a separate horizontal shere mode
!   and wake turbulence terms due to the SSO-scheme.
!  Introduction of a stablity correction for turbulent length scale.
!  Modifications for seaice model: eliminated l_ls_ice, introduced lseaice
!   (Jan-Peter Schulz)
! V4_12        2010/05/11 Ulrich Schaettler
!  Renamed t0 to t0_melt because of conflicting names
! V4_13        2010/05/11 Michael Gertz
!  Adaptions to SVN
! V4_17        2011/02/24 Ulrich Blahak
!  Adapted interface of exchg_boundaries, introduced lperi_x, lperi_y
! V4_20        2011/08/31 Ulrich Schaettler
!  Restructured SR turbdiff now as an extra module with argument list
!  for future unified COSMO-ICON physics
! V4_23        2012/05/10 Ulrich Schaettler, Oliver Fuhrer
!  Bug fix in call to SR cloud_diag: the same field is given twice in different
!    parameters. Now a local copy is done before.
!  Eliminated qvt_diff
!  Use necessary SR exchg_boundaries for ifdef __COSMO__
! V4_25        2012/09/28 Hans-Juergen Panitz
!  Introduced nexch_tag for MPI boundary exchange tag to replace ntstep (HJP)
! V4_26        2012/12/06 Matthias Raschendorfer
!  Saving scale interaction TKE source terms for output.
! V4_29        2013/10/04 Astrid Kerkweg, Ulrich Schaettler
!  Allocation of tkvm, tkvh for 1:ke1 in the vertical
! V5_1         2014-11-28 Ulrich Blahak, Matthias Raschendorfer, Oliver Fuhrer
!  Added extra advective tendency of TKE, to avoid exponential time filtering
!   of the advection, which would slow down the advection speed
!  Some preparations for COSMO-SC.
!  Introduction of horizontal diffusion coefficients. These are 0 for
!   itype_sher<2, isotropic for itype_sher=2 and have an additional
!   horizontal mode for itype_sher=3.
!  Replaced ireals by wp (working precision) (OF)
!  Enforce minimal allowed value vel_min for tke in lowest level in case of itype_tran=1 (UB)
!
! Code Description:
! Language: Fortran 90 
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!==============================================================================

! Modules used:
!-------------------------------------------------------------------------------

USE data_constants,  ONLY:   &
    t0_melt, grav=>g, con_m, con_h, cp_d, rvd_m_o, r_d,  &
    lhocp, lh_v, rcpv, rcpl, uc1, uc2, ucl,              &
    rdv, o_m_rdv, b1, b2w, b3, b4w, b234w, rdocp, b2i, b4i ! for statement_functs

USE data_flake,      ONLY:   &
    h_ice_min_flk

USE data_modelconfig,ONLY:   &
    eddlon, edadlat

USE data_fields,     ONLY:   &
    acrlat

USE data_parameters, ONLY:   &
    wp,           & ! KIND-type parameter for real variables
    iintegers       ! KIND-type parameter for standard integer variables

USE data_runcontrol, ONLY:   &
    lseaice, llake, lsso, lconv, ltkesso, ltkecon,                           &
    icldm_tran, icldm_turb, itype_turb, imode_turb, itype_tran, itype_wcld,  &
    itype_sher,                                                              &
    ltmpcor, lexpcor, lnonloc, lcpfluc, limpltkediff, l3dturb,               &
#ifdef tkvfilter
    ltkvfilter,                                                              &
#endif
    ! but only for __COSMO__
    ntstep, lperi_x, lperi_y, l2dim, lscm

USE data_turbulence, ONLY:   &
    zt_ice, z0_ice, len_min, tet_g, rim, akt, l_scal, alpha0,                &
    alpha1, rlam_mom, rlam_heat, rat_lam, it_end, vel_min, tkesmot,          &
    c_scld, a_stab, tkhmin, tkmmin, c_tke, l_hori, securi, wichfakt, c_diff, &
    clc_diag, q_crit, a_hshr, pat_len,                                       &
    a_3, a_5, a_6,           b_1, b_2,                                       &
    d_1, d_2, d_3, d_4, d_5, d_6,                                            &
    c_g,                                                                     &
    d_m=>d_mom, a_h=>a_heat, a_m=>a_mom, d_h=>d_heat

USE meteo_utilities,       ONLY:   cloud_diag
USE turbulence_utilities,  ONLY:   stab_funct, turb_param

#ifdef __COSMO__
USE data_parallel,      ONLY :  &
    num_compute,     & ! number of compute PEs
    my_cart_neigh,   & ! neighbors of this subdomain in the cartesian grid
                       ! load-imbalance
    ncomm_type,      & ! type of communication
    icomm_cart,      & ! communicator for the virtual cartesian topology
    imp_reals,       & ! determines the correct REAL type used in the model
                       ! for MPI
    nboundlines,     & ! number of boundary lines of the domain for which
                       ! no forecast is computed = overlapping boundary
                       ! lines of the subdomains
    nexch_tag,       & ! tag to be used for MPI boundary exchange
                       !  (in calls to exchg_boundaries)
    sendbuf,         & ! sending buffer for boundary exchange:
                       ! 1-4 are used for sending, 5-8 are used for receiving
    isendbuflen        ! length of one column of sendbuf

USE environment,           ONLY:   exchg_boundaries
#endif

!SCLM---------------------------------------------------------------------------
USE data_1d_global, ONLY : &

    i_cal, im, jm, &

    UUA, VVA, UWA, VWA, WWA, UST, TTA, TWA, SHF, LHF, &
    TKE_SCLM=>TKE, BOYPR, SHRPR, DISSI, TRANP
!SCLM---------------------------------------------------------------------------

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

!==============================================================================

IMPLICIT NONE

!==============================================================================

REAL (KIND=wp),     PARAMETER :: &
    z0=0.0_wp,&
    z1=1.0_wp,&
    z2=2.0_wp,&
    z3=3.0_wp,&
    z4=4.0_wp,&
    z5=5.0_wp,&
    z6=6.0_wp,&
    z7=7.0_wp,&
    z8=8.0_wp,&
    z9=9.0_wp,& 
    z10=10.0_wp

REAL (KIND=wp)     :: &
    z1d2=z1/z2,&
    z1d3=z1/z3,&
    z2d3=z2/z3,&
    z3d2=z3/z2

INTEGER (KIND=iintegers) :: &
    istat=0, ilocstat=0

LOGICAL :: lerror=.FALSE.

!==============================================================================

CONTAINS

!==============================================================================
!+ Module procedure for computing the tendencies for vertical diffusion

SUBROUTINE turbdiff (dt_var, dt_tke, lstfnct, iini, ntur, ntim, nvor,         &
                     ie, je, ke, ke1, kcm, vst,                               &
                     istartpar, iendpar, jstartpar, jendpar,                  &
                     istart   , iend   , jstart   , jend   ,                  &
                     istartu  , iendu  , jstartu  , jendu  ,                  &
                     istartv  , iendv  , jstartv  , jendv  ,                  &
                     hhl, fr_land, depth_lk, d_pat, c_big, c_sml, r_air,      &
                     dp0,                                                     &
                     u, v, w, t, qv, qc, prs, ps, qv_s, t_g, h_ice,           &
                     gz0, tcm, tch, tfm, tfh, tfv,                            &
                     tke, rcld,                                               &
                     tkvm, tkvh, tkhm, tkhh,                                  &
                     epr, rho,                                                &
                     u_tens, v_tens, t_tens, qv_tens, qc_tens, tketens,       &
                     tketens_adv,                                             &
                     ut_sso, vt_sso, tket_conv, tket_sso, tket_hshr,          &
                     edr,                                                     &
                     yerrormsg, ierrstat)

!------------------------------------------------------------------------------
!
! Description:
!
!     Es werden die Diffusionskoeffizienten berechnet und ggf. Anteile 
!     der zeitlichen Tendenzen der turbulenten Diffusion bestimmt 
!     und zu den Tendenzfeldern hinzuaddiert. 
!     Optional wird eine explizite oder (teil-)implizite Berechnung der
!     Diffusionstendenzen oder aber nur eine Berechnung der Diffusions-
!     koeffizienten durchgefuehrt. Im letzten Fall wird dann ein 
!     implizit zu berechnender Anteil der Diffusionstendenzen an 
!     anderer Stelle (slow_tendencies) bestimmt. 
!     Allerdings koennen dann zusaetzliche explizite Korrekturtendenzen
!     hier in tubdiff bestimmt werden.
!
! Method:
!
!     Die Berechnung basiert auf einer Schliessung 2-ter Ordnung auf 
!     dem level 2.5 (nach Mellor/Yamada). Demnach wird also eine 
!     prognostische Gleichung fuer die TKE geloest. 
!     Ausser der TKE-Advektion, die zusammen mit den Advektionstendenzen
!     der anderen prognostischen Variablen an anderer Stelle berechnet 
!     wird, geschieht die gesamte TKE-Prognose in diesem Unterprogramm.

!     Die Formulierung des Schemas erfolgt mit thermodynamischen 
!     Variablen, die bei feuchtadiabatischen vertikalen Verrueckungen
!     erhalten bleiben (pot. Fluessigw.temp. und Gesamtwassergehalt), 
!     so dass der Kondesationseffekt auf subskalige Vertikalbewegungen
!     beruecksichtigt wird.
!     Die turbulenten Flussdichten der Erhaltungsgroessen werden in 
!     solche der Modellvariablen konvertiert, so dass die thermodyn. 
!     Kopplung der Flussdichten richtig erhalten bleibt.
!
!     Angeschlossen ist auch ein optionales statistisches Wolkenschema
!     (nach Sommeria und Deardorff), sub turb_cloud, welches auch
!     subskalige Bewoelkung mit Hilfe der ueber das Feld rcld ausge-
!     gebenen Standardabweichung des Saettigungsdefizites berechnet.
!
!     Das Turbulenzschema wurde so verallgemeinert, dass es auch bei
!     einer vertikal aufgeloesten Bestandesschicht gueltig ist, indem
!     idealisierend von der Durchstroemung eines poroesen Mediums
!     ausgegangen wird. Die Bilanzgleichungen 1-ter und 2-ter Ordnung
!     enthalten dann zusaetzliche Terme, welche die Wechselwirkungen
!     mit der Bestandes-Matrix beschreiben. Dies wirkt sich zum einen
!     auf die Stabilitaetsfunktionen und zum anderen vor allem auf die
!     TKE-Gleichung aus, welche einen auf den Formwiderstand der 
!     Bestandeselemente zurueckzufuehrenden zusaetzlichen Quellterm
!     (Nachlaufturbulenz) enthaelt. Ausserdem werden die turbulenten
!     Flussdichtedivergenzen noch um einen Zusatzterm, welcher der 
!     Reduktion des lufterfuellten Volumens im Gitterelement Rechnung
!     traegt, erweitert. Der Effekt des Formwiderstandes in der 
!     Impulsgleichung ist ebenfalls beruecksichtigt. Die zusaetzlichen
!     Tendenzterme, die auf die Flussdichten zwischen Bestandes-Matrix
!     und umgebender Luft zurueckzufuehren sind (Bestandesquellen),
!     muessen noch in einem separaten Bestandesmodell parametrisiert
!     werden und sind nicht Gegenstand des Turbulenzmodells.
!     Schliesslich wird auch der Effekt der Transformation von Turbulenz
!     auf der dominierenden Skala in kleinskalige dissipative Turbulenz
!     durch Wirbelbrechen an Koerpern mit 
!          Laengenskalen der Abmessungen << turbulente Laengenskala 
!     beruecksichtigt; was sich durch eine (von der Laengenskala und 
!     Volumendichte jener sehr kleinen Bestandeselemente aubhaengige)
!     Modifikation der Modellkonstanten ausdruecken laesst.
!
!     Es wird auch versucht den Effekt thermisch induzierter 
!     Zirkulationen auf die TKE-Produktion zu beruecksichtigen, was 
!     durch eine Parametrisierung des Drucktransport-Termes erfolgt.
!     Hierdurch wird (vor allem) der Austauch in der naechtlichen 
!     Grenzschicht erhoeht, was der Tendenz des alten Schemas, 
!     in Bodennaehe zu kalte und nicht schnell genug anwachsende 
!     Inversionen zu produzieren, entgegenwirkt.   
!
!     Optional kann die Berechnung der vertikalen Gradienten durch eine
!     nicht-lokale Variante erfolgen. Hierbei werden die Gradienten
!     mit Profilen gebildet, die mit einem ueber die stabilitaets-
!     abhaengige Laengenskala gebildeten gleitenden Mittel behandelt
!     wurden.
!
!     Die Bildung der Anteile der Diffusionstendenzen, die numerisch 
!     durch die Multiplikation mit einer Tridiagonalmatrix ausdrueckbar 
!     sind, kann (neben der expliziten Variante) auch implizit erfolgen
!     (im Falle der Berechnung von nich-lokalen Gradienten und fuer die
!     TKE allerdings nur explizit).
!     Bei expliziter Rechnung ist, um auch bei Zeitschritten von 
!     mehreren Minuten numerisch stabil zu bleiben, eine Limitierung
!     der Groesse der Diffusionskoeffezienten und eine im vertikalen Integral
!     quellenfreie numerische Glaettung der Vertikalprofile der Diffusions-
!     tendenzen erforderlich, sowie eine teilimplizite Behandlung der
!     Diffusionstendenzen in der untersten Modellschicht, erforderlich.
!
!     Die unteren Randwerte der turbulenten Flussdichten werden ueber 
!     die Transferkoeffizienten zwischen Erdboden und unterster 
!     Modellschicht (tcm und tch) bestimmt.
!     Optional koennen die Transferkoeffizienten auch mit diesem
!     Unterprogramm bestimmt werden, indem das Turbulenzmodell auch auf
!     das Niveau d+z0 angewandt wird, wobei vertikale Gradienten in 
!     diesem Niveau mit Hilfe der Prandtl-Schicht-Hypothese berechnet
!     werden.
!     In diesem Zusammenhang wird auch die Wirkung der laminaren 
!     Grenzschicht behandelt.
!
!     Turbulente Horizontaldiffusion (um ein 3-d-Schema zu erhalten)
!     ist noch nicht enthalten, kann aber integriert werden. 
!     Uebergabevariablen:
!
!------------------------------------------------------------------------------

! Declarations:

!------------------------------------------------------------------------------

! Arguments with intent(in):

  REAL (KIND=wp),           INTENT (IN) :: &
    dt_var,  & ! time step for turbulent diffusion (all prognostic var but TKE)
    dt_tke     ! time step for TKE stepping

  INTEGER (KIND=iintegers), INTENT(IN) :: &
    iini,    & ! type of initialization (0: no, 1: separate before the time loop
               !                              , 2: within the first time step)
    ntim       ! number of tke time levels 

  INTEGER (KIND=iintegers), INTENT(INOUT) :: &
    nvor,    & ! actual time step index of TKE (also used in turbdiff)
    ntur       ! current new time step index of tke 

  INTEGER (KIND=iintegers), INTENT (IN) :: &
    ie, je, ke, ke1,                       & ! dimensions of the array
    istartpar, iendpar, jstartpar, jendpar,& ! start- and end-indices for computation
    istart   , iend   , jstart   , jend   ,&
    istartu  , iendu  , jstartu  , jendu  ,&
    istartv  , iendv  , jstartv  , jendv

  INTEGER (KIND=iintegers), INTENT (INOUT) :: &
    kcm        ! level index of the upper canopy bound
               !  (might be modified in init_canopy!!)

  INTEGER (KIND=iintegers), INTENT (IN) :: &
    vst             ! velocity component staggering

  LOGICAL,                  INTENT (IN) :: &
    lstfnct

!------------------------------------------------------------------------------

  REAL (KIND=wp),           INTENT (IN) :: &

    ! External parameter fields:
    ! --------------------------

    hhl     (ie,je,ke1), & ! height of model half levels                   ( m )
    dp0     (ie,je,ke),  & ! pressure thickness of layer                   (pa )
    fr_land (ie,je),     & ! land portion of a grid point area             ( 1 )
    depth_lk(ie,je)        ! lake depth                                    ( m )
 
  REAL (KIND=wp),           INTENT (IN), TARGET           :: &
!
    d_pat(ie,je),        & ! horizontal pattern length scale
    c_big(ie,je,kcm:ke1),& ! effective drag coefficient of canopy elements
                           ! larger than or equal to the turbulent length scale (1/m)
    c_sml(ie,je,kcm:ke1),& ! effective drag coefficient of canopy elements
                           ! smaller than the turbulent length scale            (1/m)
    r_air(ie,je,kcm:ke1)   ! log of air containing fraction of a gridbox inside
                           ! the canopy                                          (1)

!------------------------------------------------------------------------------
    
  REAL (KIND=wp),           INTENT(IN)    :: &

    ! Fields for surface values and soil/canopy model variables:
    ! ----------------------------------------------------------

    h_ice(ie,je),        & ! ice thickness                                 (  m  )
    ps   (ie,je),        & ! surface pressure                              ( pa  )
    qv_s (ie,je),        & ! specific water vapor content on the surface   (kg/kg)
    t_g  (ie,je)           ! specific water vapor content on the surface   (kg/kg)

!------------------------------------------------------------------------------
    
  REAL (KIND=wp),           INTENT(INOUT), TARGET :: &

    ! Atmospheric model variables:
    ! ----------------------------

!US the latest news from ICON: not only the microphysics, but also the turbulence
!   scheme (the fast physics) should update the prognostic variables, not the
!   tendencies

    u  (ie,je,ke),       & ! zonal wind speed                              ( m/s )
    v  (ie,je,ke),       & ! meridional wind speed                         ( m/s )
    w  (ie,je,ke1),      & ! vertical wind speed (defined on half levels)  ( m/s )
    t  (ie,je,ke),       & ! temperature                                   (  k  )
    qv (ie,je,ke),       & ! specific water vapor content                  (kg/kg)
    qc (ie,je,ke),       & ! specific cloud water content                  (kg/kg)
    prs(ie,je,ke)          ! full pressure                                 ( pa  )

!------------------------------------------------------------------------------

  REAL (KIND=wp),           INTENT(INOUT) :: &

    ! Atmospheric variables of the turbulence model:
    ! ----------------------------------------------

    tke(ie,je,ke1,ntim), & ! SQRT(2*TKE); TKE='turbul. kin. energy'        ( m/s )
                           ! (defined on half levels)
    tkvm(ie,je,ke1),     & ! turbulent diffusion coefficients for momentum (m/s2 )
    tkvh(ie,je,ke1),     & ! turbulent diffusion coefficients for heat     (m/s2 )
                           ! and moisture
    rcld(ie,je,ke1),     & ! standard deviation of the saturation deficit
                           ! (as input and output)
                           ! fractional cloud cover (in turbdiff)            --
    tkhm(ie,je,ke),      & ! turbulent diffusion coefficients for momentum (m/s2 )
    tkhh(ie,je,ke)         ! turbulent diffusion coefficients for heat     (m/s2 )

!------------------------------------------------------------------------------

  REAL (KIND=wp),           INTENT(IN), TARGET, OPTIONAL :: &

     rho(ie,je,ke),      & ! total density of air                          (kg/m3)
     epr(ie,je,ke)         ! exner pressure                                 (1)

!------------------------------------------------------------------------------

  REAL (KIND=wp),           INTENT(INOUT), TARGET, OPTIONAL :: &

    ! Tendency fields for the prognostic variables:
    ! ---------------------------------------------

    u_tens (ie,je,ke),  & ! u-tendency                                    ( m/s2)
    v_tens (ie,je,ke),  & ! v-tendency                                    ( m/s2)
    t_tens (ie,je,ke),  & ! t-tendency                                    ( K/s )
    qv_tens(ie,je,ke),  & ! qd-tendency                                   ( 1/s )
    qc_tens(ie,je,ke)     ! qw-tendency                                   ( 1/s )

  REAL (KIND=wp),           INTENT(IN),  OPTIONAL :: &
    ut_sso(ie,je,ke),   & ! u-tendency due to the SSO-Scheme              ( 1/s )
    vt_sso(ie,je,ke),   & ! v-tendency due to the SSO-Scheme              ( 1/s )
    tket_conv(ie,je,ke1)  ! TKE-tendency due to convective buoyancy       (m2/s3)

  REAL (KIND=wp),           INTENT(OUT), OPTIONAL :: &
    edr     (ie,je,ke1),& ! eddy dissipation rate of TKE (EDR)            (m2/s3)
    tket_hshr(ie,je,ke1),&! TKE-tendency due to (sep.) horiz. shear       (m2/s3)
    tket_sso(ie,je,ke1)   ! TKE-tendency due to SSO wake production       (m2/s3)

  REAL (KIND=wp),           INTENT(INOUT) :: &

    tketens(ie,je,ke1)    ! tendency of SQRT(2*TKE)                       ( m/s2)

  REAL (KIND=wp),           INTENT(IN) :: &
    tketens_adv(ie,je,ke1)  ! pure advective tendency of SQRT(2*TKE)      ( m/s2)

!------------------------------------------------------------------------------

  REAL (KIND=wp),           INTENT(INOUT) :: &

    ! Diagnostic surface variable of the turbulence model:
    ! ----------------------------------------------------

    gz0(ie,je),          & ! roughness length * g of the vertically not
                           ! resolved canopy                               (m2/s2)
    tcm(ie,je),          & ! turbulent transfer coefficients for momentum    --
    tch(ie,je),          & ! turbulent transfer coefficients for heat        --
    tfm(ie,je),          & ! factor of laminar transfer of momentum          --
    tfh(ie,je),          & ! factor of laminar transfer of scalars           --
    tfv(ie,je)             ! laminar reduction factor for evaporation        --
!------------------------------------------------------------------------------

!US  REAL (KIND=wp),           INTENT(OUT), OPTIONAL   :: &
!US
!US     shfl_s(ie,je),      & ! sensible heat flux at the surface             (W/m2) (positive upward)
!US     lhfl_s(ie,je)         ! latent   heat flux at the surface             (W/m2) (positive upward)

!------------------------------------------------------------------------------
  
  INTEGER (KIND=iintegers), INTENT(OUT)   :: &
    ierrstat     ! error status
  
  CHARACTER (LEN= *),       INTENT(OUT)   :: &
    yerrormsg    ! error message

!------------------------------------------------------------------------------
!-----------------------------------------------------------------------
!
! Local parameters:

! Local scalars:

!     Lokale logical Variablen:

      LOGICAL limplizit,& !(teil-)implizite Berechnung der Diffus.tend.
              lexplizit,& !(teil-)explizite Berechnung der Diffus.tend.
              linitial !Initialisierungsdurchgang bei Modellstart

!     Lokale Integer-Parameter:
     
      INTEGER (KIND=iintegers), PARAMETER :: &

!             Indexgrenzen:

              nscal=3,     & !Skalare Groessen
              ninv=2,      & !daraus abgeleitete gegenueber vertikalen 
                             !(feuchtadiabatischen) Verrueckungen 
                             !invarianten Groessen 
              nvel=2,      & !Geschwindigkeitskomponenten
              ndiff=nscal+nvel, &
              nred=ninv+nvel,   &
              ntmax=3,     & !max. Anzahl der Zeitebenen fuer die TKE

!             Zeiger fuer die Variablen :

              u_m=1,       & !zonale Geschw.komp. im Massenzentrum
              v_m=2,       & !meridionale  ,,      ,,     ,, 
              tet_l=3,     & !feucht-potentielle Temperatur
              h2o_g=4,     & !Gesamtwasseergehalt

              u_s=1,       & !zonale Geschw.komp. im gestag. Gitter 
              v_s=2,       & !meridionale  ,,     ,,   ,,       ,, 
              tem=3,       & !Temperatur                  
              tet=3,       & !pot.Temperatur                  
              vap=4,       & !Wasserdampfmischungsverh.
              liq=5          !Fluessigwasser  ,,              

      INTEGER (KIND=iintegers)            :: &
              izerror,     & !
              ntlev          !tatsaechliche Anzahl der Zeitebenen fuer die TKE

!     Beachte:u_m,v_m, sowie u_s,v_s muessen in [1,nvel] liegen;
!             aber: tet_l,h2o_g in [nvel+1,nred]
!             und   tem (tet),vap,liq in [nvel+1,ndiff]

!     Lokale Integer-Hilfsvariablen:

      INTEGER (KIND=iintegers) :: &

              izloc,jzloc,  & !Startindices bei Windinterpolation
              i,j,k,        & !Diskretis.index fuer lamda, phi, sigma
              kzdims(24),   & !dimensions of fields to be exchanged
              n,m,ii,jj,kk, & !Indices fuer diverse Schleifen
              ku,ko,k1,k2,  & !unterer und oberer Schicht-Index
              nkorr,nk1,nk2,& !Startindices der explizit zu behandelnden
                              !Variablen
              kem,          & !ke oder ke1
              it_start,     & !Startindex der Iterationen
              it_durch        !Durchgangsindex der Iterationen

      INTEGER (KIND=iintegers) :: &

              tkestep(ntmax),&!Zeitstufen des tke-Feldes

!             Eingrenzende Hoehenvieaus bei der Berechnung der 
!             gemittelten Profile bei nicht-lokaler Gradient-Berechnung:

              lev(2) 


!     Lokale real Variablen:

      REAL (KIND=wp)     :: &

!          Hilfsvariablen:

           wert,& ! Platzhalter fuer beliebige Zwischenergebnisse
           fakt,& !  ,,         ,,     ,,     Faktoren
           sum,&  !  ,,         ,,     ,,     Summen

!     Platzh. fuer pot. Temp., relat. Wolkenanteil, Exner- und rezp. 
!     virt. Faktor, sowie fuer die Temperaturtendenz des Wasser-
!     saettigungsmisch.verh.:

           teta,rc,ex,virt,alf,&

!     Platzh. fuer therm. und mech. Antrieb der Turbulenz in (1/s)**2 
!     (fh2,fm2), die Diffusionskoeff. (bzw. Transferkoeff.) fuer 
!     Skalare und Impuls in m^2/s (kh,km), Stabilitaetsfunkt. (sh,sm),
!     imodifizierte Transferkoeff. (tmch,tmcm) in Kg/m^2/s:

           fh2,fm2,kh,km,khc,sh,sm,tmch,tmcm,&

!     Platzh. fuer horiz. Geschw.-Komponenten und Geschw.betrag:

           vl1,vl2,vel,&

!     Platzh. fuer SQRT(2*TKE)-Werte:

           q0,q1,q2,q3,&  
           q0vec(ie,je), &

!     Platzh. fuer die Hoehe ueber Grund, Hoehendifferenzen, obere und
!     untere Hoehe, turbulente Laengaenskalen, Kohaerenzlaenge,
!     Dicke der laminaren Grenzschicht,sowie eine Laenge allgemein:

           h,hu,ho,l_turb,ld,lh,lm,kohae_len,l_lam,len,&

           edh,& ! Kehrwert von Schichtdicken

!     Platzhalter fuer Druckdifferenzen

           dpo,dpu,dp2,&

!     Zwischenspeicher fuer rcpv, rcpl und lhocp

           zrcpv,zrcpl,zlhocp,&

!     Zwischenspeicher fuer

           thermik,& !(negative) Auftriebsproduktion an TKE
           phasdif,& !Temperaturtendenz durch Phasendiffusion

           liqfak ,& !Faktor zur Beruecksichtigung von Fluessigwasser

!     sonstiges:

           tmp1,tmp2,tmp3,aa,bb,&
           alfa,beta,gama,x0,x1,x2,x3,xx

      REAL (KIND=wp)     :: flukon33, flukon43, flukon53, &
                            flukon34, flukon44, flukon54, flukon55, &
                            flux_3,   flux_4,   flux_5, &
                            zs11, zs22, zs33, zs12, zs21, zs13, zs23

! Local arrays:

      LOGICAL            :: lo_ice(ie,je) ! logical sea ice indicator

!     Lokale Hilfsfelder:

      REAL (KIND=wp)     :: &
           rc_2d(ie,je)    !

      REAL (KIND=wp)     :: &

           dicke(ie,je,ke1),& !(effektive) Dicke einer Modellschicht 
                              !bzw. sonstiges Hilfsfeld, bei k=ke1 steht
                              !eine effekt. Dicke der Modell-Prandtlsch.

           wind(ie,je,ke1,nvel),& !horizontale Windkomponenten im
                                  !Massenzentrum, bzw. deren Tendenzen

           rhon(ie,je,ke1),& !Luftdichte auf Nebenflaechen
                             !einschl. des skin-layers

!     Hilfsfeld fuer die korrigierten effektiven Gradienten der
!     ndiff diffundierenden Modellvariablen:

           a(ie,je,ke1,ndiff), &

!          a() enthaelt zu beginn die Werte der 5 thermodynamischen
!          Koeffizienten dQs/dT,ex_fakt,cp_fakt,A_tet,A_q 
!          auf den Positionen 1 bis 5 des letzten Index.

!          Im Falle der impliziten Berechnung von Flussdivergenzen
!          enthaelt a() hierfuer benoetigte Hilfsgroessen

!     3-d Hilfsfelder

           hlp (ie,je,ke1), &
           hlp2(ie,je,ke1), &
#ifdef tkvfilter
           tkvm_filter(ie,je,ke1), &
           tkvh_filter(ie,je,ke1), &
#endif

!     Flusskonversionsmatrix, um von den turb. Flussdichten der 2
!     scalaren Erhaltungsgroessen (bzgl. feuchtadiab. vert. Verrueck.)
!     auf diejenigen der nscal scalaren Modellvariablen zu gelangen:

           flukon(nvel+1:ndiff,nvel+1:nred), &

!     Hilfsfelder fuer eine und zwei Variablenschichten:

           lay(ie,je), lays(ie,je,2), &

           pr(ie,je),   & !Druck             (am atm. Modellunterrand)
           tp(ie,je),   & !Temperatur              ,,       
           qd(ie,je),   & !Wasserdampfgehalt       ,,
           ql(ie,je),   & !Fluessigwassergehalt    ,,
           z0m(ie,je),  & !Rauhigkeitslaenge (fuer Impuls)
           l_pat(ie,je),& !Laengenskala der therm. Inhomogenitaeten
                          !der Erdbodenoberflaeche

           frc(ie,je),  & !forcing
           src(ie,je),  & !source term
           src2(ie,je), & !source term

!     Vertikale Gradienten und turbulente Flussdichten verschiedener 
!     thermodynamischer Variablen:

           grad(ndiff), flux(ndiff), varc(ndiff), &

           hig(2) ! obere und untere Referenzhoehe bei der Bildung
                  ! nicht-lokaler Gradienten

      LOGICAL                  :: &
        lini

CHARACTER (LEN=80)             :: &
  yzerrmsg      ! error message


REAL (KIND=wp),     TARGET :: &
   len_scale(ie,je,ke1),       & !
        vari(ie,je,ke1,ndiff), & !
          dd(ie,je,0:7)          ! local derived turbulence parameter

!     ------------------------------------------------------------------
!     Festlegung der Formelfunktionen fuer die Turbulenzparametrisierung:

!     ------------------------------------------------------------------
      INCLUDE 'statement_functs.incf' 
!     ------------------------------------------------------------------

! - end of header
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!     nur in LM-Umgebung:

      istat=0
      izerror  = 0
      yzerrmsg = '   '

!     abhaengig, ob das 2- oder das 3-Zeitebenenverfahren angewendet
!     wird, muss ntlev=2 bzw. =3 gesetzt werden:
!US      IF (l2tls) THEN
!US        ntlev = 2
!US      ELSE
!US        ntlev = 3
!US      ENDIF
      ntlev = ntim

!     Fuer die Turb.par. benutzter Variablensatz auf Hauptflaechen:
!     Bei k=ke1 stehen die unteren Randwerte der Prandtlschicht
!     (skin-layer-Werte) 

!     Der letzte Index bezeichnet die physik. Bedeutung der Variablen
!     und hat einen der Werte u_m,v_m,tet_l,h2o_g,liq;
!                        bzw. u_s,v_s,tem,vap,liq.
!     Der letzte Index bezeichnet die physik. Bedeutung der Variablen
!     Bei k=ke1 stehen die unteren Randwerte der Prandtlschicht
!     (skin-layer-Werte)
!     Am Ende des Programms wird vari mit den (Co-Varianzen) der
!     Geschwindigkeitskomponenten ueberschrieben, die fuer die
!     Berechung der Horizontaldiffusion benoetigt werden.
!     vari() enthaelt spaeter auch die ndiff (nichtlokalen) vertikalen
!     Gradienten und auch die durch Wirkung der subskaligen Kondensation
!     veraenderten (effektiven) vertikalen Gradienten.
!     Zum Schluss enthaelt vari() fuer die turbulente Horizontaldiff.
!     benoetigte Komponenten des turbulenten Spannungstensors.

!-----------------------------------------------------------------------


! 1)  Vorbereitungen:

 IF (iini.GT.0) THEN !an initialization run
!   lini=.TRUE.
    linitial=.TRUE.
    IF (iini.EQ.1) THEN !separate initialization before the time loop
       it_start=1 !only 'it_end' iterations for initialization
                  !and an additional one at the first time loop
    ELSE !initialization within the first time step
       it_start=0 !"it_end+1" iterations for initializatzion
    END IF
 ELSE !not an initialization run
!   lini=.FALSE.
    linitial=.FALSE.
    it_start=it_end !only one iteration
 END IF

!     Setzen von Steuerparametern:

      IF (imode_turb.EQ.3) THEN
         limplizit=.TRUE. 
      ELSE
         limplizit=.FALSE.
      END IF

      IF (imode_turb.EQ.2) THEN
         lexplizit=.TRUE.
         nkorr=1
      ELSE
         IF (lexpcor) THEN
            lexplizit=.TRUE.
            IF (lnonloc) THEN
               nkorr=1
            ELSE
               nkorr=nvel+1
            END IF
         ELSE
            lexplizit=.FALSE.
            nkorr=ndiff+1
         END IF
      END IF

      nk1=nkorr
      nk2=MAX(nkorr,nvel+1)

      IF (itype_tran.EQ.3) THEN
         kem=ke1
      ELSE
         kem=ke
      END IF   

!US   IF (ntke.EQ.0) THEN
!US      linitial=.true.
!US      it_start=0
!US   ELSE
!US      linitial=.false.
!US      it_start=it_end
!US   END IF

      IF (lcpfluc) THEN
         zrcpv=rcpv
         zrcpl=rcpl
      ELSE
         zrcpv=z0  
         zrcpl=z0  
      END IF

!     Berechnung abgeleiteter Parameter:

      CALL turb_param (istartpar, iendpar, jstartpar, jendpar, grav, cp_d, dd)

!US!     Bestimmung der Zeitebenen des TKE-Feldes:
!US
!US      DO n=1,ntlev
!US         tkestep(n)=MOD(ntke+n,ntlev)+1
!US      END DO
!US
!US      nvor=tkestep(ntlev-1)
!US      ntke=tkestep(ntlev)

!     Bestimmung der initialen Werte fuer die laminaren Reduktions-
!     koeffizienten :

      IF (linitial.AND.itype_tran.NE.2) THEN 
         DO j=jstartpar,jendpar
         DO i=istartpar,iendpar    
            tfh(i,j)=z1
            tfm(i,j)=z1
         END DO
         END DO
      END IF

!     Berechnung der effektiven Bedeckung mit nichtkonvektiven
!     Wasserwolken:

      IF (icldm_turb.EQ.-1) THEN       
!        Keine Wolkenberuecksichtigung;
!        Wolkenwasser wird ignoriert (vollkommen trockenes Schema):

         DO k=1,ke
            DO j=jstartpar,jendpar
            DO i=istartpar,iendpar    
               rcld(i,j,k)=z0
               hlp(i,j,k)=z0
            END DO
            END DO
         END DO

         zlhocp=z0 !no condensation
         liqfak=z0 !total water without liquid water
      
      ELSE

      zlhocp=lhocp; liqfak=z1

      IF (icldm_turb.EQ.0) THEN       

!        Keine Wolkenberuecksichtigung;
!        alles Wolkenwasser wird verdunstet:

         DO k=1,ke
            DO j=jstartpar,jendpar
            DO i=istartpar,iendpar    
               rcld(i,j,k)=z0
               hlp(i,j,k)=-qc(i,j,k)
            END DO
            END DO
         END DO

      ELSEIF (icldm_turb.EQ.1) THEN
!       Wolkenwasser wird als skalige Wolke interpretiert:

        DO k=1,ke
          DO j=jstartpar,jendpar
            DO i=istartpar,iendpar
              IF ( qc(i,j,k) > z0 ) THEN
                rcld(i,j,k) = z1
              ELSE
                rcld(i,j,k) = z0
              END IF
              hlp(i,j,k)=z0
            ENDDO
          ENDDO
        ENDDO

      ELSEIF (icldm_turb.EQ.2) THEN
!       spezielle Diagnose von Wasserwolken:

        hlp2(:,:,:) = rcld(:,:,:)
        CALL cloud_diag(rcld, hlp,                                     &
           istartpar,iendpar,jstartpar,jendpar,1,ke,                   &
           1,ie,1,je,1,ke1,                                            &
           ie,je,ke,ke1,                                               &
           rdv,o_m_rdv,rvd_m_o,lhocp,t0_melt,                          &
           b1,b2w,b3,b4w,b234w,b2i,b4i,                                &
           uc1,uc2,ucl,clc_diag,q_crit,                                &
           t(:,:,:), qv(:,:,:), qc(:,:,:), prs(:,:,:), hlp2, ps(:,:),  &
           itype_wcld)

         DO k=1,ke
            DO j=jstartpar,jendpar
            DO i=istartpar,iendpar    
               hlp(i,j,k)=hlp(i,j,k)-qc(i,j,k)
            END DO
            END DO
         END DO

      END IF
      END IF
      
      DO j=jstartpar,jendpar
      DO i=istartpar,iendpar    
         hlp(i,j,ke1)=hlp(i,j,ke)
      END DO
      END DO

      IF (icldm_tran.EQ.0) THEN
         DO j=jstartpar,jendpar
         DO i=istartpar,iendpar    
            rcld(i,j,ke1)=z0
         END DO
         END DO
      ELSEIF (icldm_tran.EQ.1) THEN
         DO j=jstartpar,jendpar
           DO i=istartpar,iendpar    
             IF ( qc(i,j,ke) > z0 ) THEN
               rcld(i,j,ke1) = z1
             ELSE
               rcld(i,j,ke1) = z0
             END IF
           END DO
         END DO
      ELSE
         DO j=jstartpar,jendpar
         DO i=istartpar,iendpar    
            rcld(i,j,ke1)=rcld(i,j,ke)
         END DO
         END DO
      END IF


!     Berechnung der thermodynamischen Hilfsfelder:

!     Bem.:  Am unteren Modellrand gibt es eine skalige Erdboden-
!            oberflaeche mit einer vertikal nicht aufgeloesten 
!            Rauhigkeit der Laenge d+z0m. Diese von Rauhigkeits-
!            elementen durchsetzte 'Mikro'-Rauhigkeitsschicht zaehlt
!            nicht mehr zum atm. Modellgebiet! Sie stellt einen 
!            Mikrobestand der Hoehe d+z0m dar, in dem die konstante
!            Schliessungslaenge akt*z0m gelten soll.

!            Es wird davon ausgegangen, dass das Niveau der Hoehe
!            hhl(i,j,ke1) (Erdbodenniveau) gerade der Untergrenze der
!            Prandtlschicht entspricht, also um d+z0m hoeher als die
!            Hoehe des eigentlichen (festen) Erdbodens liegt, (wobei
!            d die Verdraengungshoehe der Rauhigkeitselemente und 
!            z0m die effektive Rauhigkeitslaenge ist.

!            Hierbei wird angenommen, dass es am unteren Modellrand
!            keine turb. Fluessigwasserfluesse gibt, dass also der
!            Randwert fuer liq dem Wert der untersten Modellschicht
!            entspricht. Die Flussdichte fuer tet_l entspricht dann
!            der fuer tet und die fuer h2o_g der fuer vap.

      DO k=1,ke1
         IF (k.LT.ke1) THEN
            DO j=jstartpar,jendpar
            DO i=istartpar,iendpar    
               tp(i,j)=t(i,j,k)+zlhocp*hlp(i,j,k)
               qd(i,j)=qv(i,j,k)-hlp(i,j,k)
               ql(i,j)=qc(i,j,k)+hlp(i,j,k)
               pr(i,j)=prs(i,j,k)
            END DO   
            END DO   
         ELSE
!           Beachte:
!           tp(), qd(), ql() und pr() sind vom Durchgang k=ke
!           bereits vorhanden.

!           Es werden diese Variablen fuer den atm. Unterrand bestimmt.
!           Dabei wird ueber die Reduktionskoeff. tfh und tfm 
!           die Wirkung der laminaren Grenzschicht zwischen der 
!           Erdbodenumsatzflaeche (im Abstand einer eventuellen
!           Verdraengungshoehe von der festen Erdbodenoberflaeche)
!           und diesem Niveau beruecksichtigt:
 
            DO j=jstartpar,jendpar
            DO i=istartpar,iendpar    
               tp(i,j)=tp(i,j)*(z1-tfh(i,j)) &
!                     +(t_g(i,j)+zlhocp*hlp(i,j,ke1))*tfh(i,j)
                      +t_g(i,j)*tfh(i,j)
               qd(i,j)=qd(i,j)*(z1-tfh(i,j)) &
!                     +(qv_s(i,j)-hlp(i,j,ke1))*tfh(i,j)
                      +qv_s(i,j)*tfh(i,j)
               virt=(z1+rvd_m_o*qd(i,j)-ql(i,j))
               pr(i,j)=ps(i,j) 
            END DO
            END DO
         END IF

         DO j=jstartpar,jendpar
         DO i=istartpar,iendpar    
            fakt=c_scld/(z1+rcld(i,j,k)*(c_scld-z1))
            rc=fakt*rcld(i,j,k)
            ex=zexner(pr(i,j)) !Exner-Faktor
            virt=z1/(z1+rvd_m_o*qd(i,j)-ql(i,j)) !virtueller Faktor
            teta=tp(i,j)/ex !pot. Temp

!           Temp.tendenz der Saettigungsfeuchte:
            alf=zdqsdt(tp(i,j),zqvap(zpsat_w(tp(i,j)),pr(i,j)))

            fakt=(zlhocp/tp(i,j)-(z1+rvd_m_o)*virt)/(z1+alf*zlhocp)

            a(i,j,k,1)=z1+zrcpv*qd(i,j)+zrcpl*ql(i,j) !Cp/Cpd
            a(i,j,k,2)=ex                 !Exner-Faktor
            a(i,j,k,3)=alf                !dQs/dT

!           Berechn. der thermodyn. Faktoren A_tet und A_q:
            a(i,j,k,4)=grav*(z1/teta-rc*fakt*ex*alf)
            a(i,j,k,5)=grav*(rvd_m_o*virt+rc*fakt)

!           Berechnung der feucht-potentiellen Temp. und des
!           Gesamtwassergehaltes:
            vari(i,j,k,tet_l)=teta-ql(i,j)*zlhocp/ex
            vari(i,j,k,h2o_g)=qd(i,j)+liqfak*ql(i,j)
            vari(i,j,k,liq)=ql(i,j)

            rhon(i,j,k)=virt*pr(i,j)/(r_d*tp(i,j)) !Luftdichte

!           Die thermodynamischen Hilfsgroessen wurden hier
!           unter Beruecksichtigung der diagnostizierten
!           Kondensationskorrektur gebildet, indem die entspr.
!           korrigierten Werte fuer tem, qv und ql benutzt wurden.

         END DO
         END DO
      END DO

!     Beachte:
!     tp(), qd(), ql() und pr() gehoeren jetzt zu k=ke1 

!     Berechnung der horizontalen Windgeschwindigkeiten
!     im Massenzentrum der Gitterbox:

      DO j=jstartpar,jendpar
         jj=MAX(j-1,1)
         DO k=1,ke
            DO i=istartpar,iendpar
               ii=MAX(i-1,1)
               vari(i,j,k,u_m)=(u(i,j,k)+u(ii,j,k))*z1d2 
               vari(i,j,k,v_m)=(v(i,j,k)+v(i,jj,k))*z1d2
            END DO
         END DO
         DO i=istartpar,iendpar
            vari(i,j,ke1,u_m)=vari(i,j,ke,u_m)*(z1-tfm(i,j))
            vari(i,j,ke1,v_m)=vari(i,j,ke,v_m)*(z1-tfm(i,j))
         END DO
      END DO         

!     Sichern der auf Massenpunkte interpolierten Windkomponenten:

      DO n=1,nvel
         DO k=1,ke1
            DO j=jstartpar,jendpar
            DO i=istartpar,iendpar    
               wind(i,j,k,n)=vari(i,j,k,n)
            END DO
            END DO
         END DO
      END DO    

!     Berechnung der Modellschichtdicken:
        
      DO k=1,ke
         DO j=jstartpar,jendpar
         DO i=istartpar,iendpar    
            dicke(i,j,k)=hhl(i,j,k)-hhl(i,j,k+1)
         END DO
         END DO
      END DO   

!     Interpolation der thermodyn. Hilfsgroessen im Feld a(),
!     der Wolkendichte und der Luftdichte auf Nebenflaechen:

      DO k=ke,2,-1
         DO j=jstartpar,jendpar
         DO i=istartpar,iendpar
            hlp(i,j,k)=dp0(i,j,k)+dp0(i,j,k-1)
         END DO
         END DO
      END DO
      DO ii=1,ndiff
         DO k=ke,2,-1
            DO j=jstartpar,jendpar
            DO i=istartpar,iendpar
               a(i,j,k,ii)=(a(i,j,k,ii)*dp0(i,j,k-1)+a(i,j,k-1,ii)*dp0(i,j,k))/hlp(i,j,k)
            END DO
            END DO
         END DO
      END DO
      DO k=ke,2,-1
         DO j=jstartpar,jendpar
         DO i=istartpar,iendpar
            rcld(i,j,k)=(rcld(i,j,k)*dp0(i,j,k-1)+rcld(i,j,k-1)*dp0(i,j,k))/hlp(i,j,k)
            rhon(i,j,k)=(rhon(i,j,k)*dicke(i,j,k-1)  &
                        +rhon(i,j,k-1)*dicke(i,j,k)) &
                       /(hhl(i,j,k-1)-hhl(i,j,k+1))

!           Beachte die andere Behandlung der Dichte!
         END DO
         END DO
      END DO

!     Bestimmun der initialen Werte fuer die Rauhigkeitslaenge
!     ueber Wasserpunkten:

      ! Set the logical mask lo_ice to distinguish between ice covered
      ! and open water sea or lake grid points.

      DO j=jstartpar,jendpar
        DO i=istartpar,iendpar

          IF (fr_land(i,j).LT.z1d2) THEN
            ! Water point.
            IF (.NOT. lseaice) THEN
              ! Sea ice model is not used.
              ! Ice surface if SST is less than the salt water freezing temperature.
              lo_ice(i,j) = t_g(i,j) < t0_melt + zt_ice
            ELSE
              ! Sea ice model is used.
              ! Ice surface if ice is present.
              lo_ice(i,j) = h_ice(i,j) > 0.0_wp
            END IF
            IF (llake) THEN
              ! Lake model is used.
              ! Ice surface if this is a lake point AND ice is present.
              IF ((depth_lk(i,j) > 0.0_wp) .AND. (h_ice(i,j) >= h_Ice_min_flk)) &
              lo_ice(i,j) = .TRUE.
            END IF
          END IF

        END DO
      END DO

      IF (linitial.AND.itype_tran.EQ.3) THEN 

         DO j=jstartpar,jendpar
         DO i=istartpar,iendpar    
            IF (fr_land(i,j).LT.z1d2) THEN

               ! Use ice surface roughness or open-water surface roughness
               ! according to lo_ice
               IF ( lo_ice(i,j) ) THEN
                  gz0(i,j)=grav*z0_ice
               ELSE
!                 Bei von Schubspannung abhaengiger Wellenhoehe:

!                 Einfachste Schaetzung der Schubspannung als Impusls-
!                 flussdichte durch die Nebenflaeche ke mit Hilfe
!                 einer diagnostischen TKE und Bestimmung der Rauhig-
!                 keitslaenge unter Anwendung der Charnockformel:

                  l_turb=dicke(i,j,ke)
                  l_turb=akt*MAX(len_min,l_turb/(z1+l_turb/l_scal))
                  edh=z2/(hhl(i,j,ke-1)-hhl(i,j,ke1))

                  DO n=1,nred
                     grad(n)=(vari(i,j,ke-1,n)-vari(i,j,ke,n))*edh
                  END DO

                  fh2=a(i,j,k,4)*grad(tet_l)+a(i,j,k,5)*grad(h2o_g)
                  fm2=grad(u_m)**2+grad(v_m)**2

                  ! Vereinfachte Loesung mit Rf=Ri:
                  IF (fh2.GE.(z1-rim)*fm2) THEN
                     ! Die krit. Ri-Zahl wird ueberschritten und lm, sowie lh
                     ! werden durch lm bei der krit. Ri-Zahl angenaehert:
                     fakt=z1/rim-z1
                     lm=l_turb*(b_2-(a_6+a_3)*fakt)
                     lh=lm
                  ELSE
                     fakt=fh2/(fm2-fh2)
                     lm=l_turb*(b_2-(a_6+a_3)*fakt)
                     lh=l_turb*(b_1-a_5*fakt)
                  END IF

                  wert=lm*fm2
                  wert=MAX(wert-lh*fh2,rim*wert)

                  q0=MAX(vel_min,SQRT(d_m*l_turb*wert))

                  wert=lm*q0*SQRT(fm2) 
                  gz0(i,j)=MAX(grav*len_min,alpha0*wert+alpha1*grav*con_m/SQRT(wert))
               END IF
            END IF   
         END DO   
         END DO   
      END IF     

!     Bestimmung der Rauhigkeitslaenge 
!     und der effektiven Dicke der Modell-Prandtlschicht:

      DO j=jstartpar,jendpar
      DO i=istartpar,iendpar    
         z0m(i,j)=gz0(i,j)/grav
         dicke(i,j,ke1)=z0m(i,j) &
                       *log(z1d2*dicke(i,j,ke)/z0m(i,j)+z1)
      END DO
      END DO


!     Berechnung der tubulenten Laengenscalen:

      DO j=jstartpar,jendpar
      DO i=istartpar,iendpar
         len_scale(i,j,ke1)=z0m(i,j)
      END DO
      END DO
      DO k=ke,kcm,-1 !Innerhalb des Bestandesmodells
         DO j=jstartpar,jendpar
         DO i=istartpar,iendpar
            IF (c_big(i,j,k).gt.z0) THEN

!              Die turbulente Laengenskala wird durch die Laengen-
!              skala der lufterfuellten Zwischenraeume limitiert:

               l_turb=z1/(c_big(i,j,k)*sqrt(z1/exp(r_air(i,j,k))-z1))
               len_scale(i,j,k)=MIN(dicke(i,j,k)+len_scale(i,j,k+1) &
                                   ,l_turb)
            ELSE
               len_scale(i,j,k)=dicke(i,j,k)+len_scale(i,j,k+1)
            END IF
         END DO
         END DO
      END DO   
      DO k=kcm-1,1,-1
         DO j=jstartpar,jendpar
         DO i=istartpar,iendpar
            len_scale(i,j,k)=dicke(i,j,k)+len_scale(i,j,k+1)
         END DO
         END DO
      END DO   

!     Uebergang von der maximalen turbulenten Laengenskala zur
!     effektiven turbulenten Laengenskala:

      DO k=1,ke1
         DO j=jstartpar,jendpar
         DO i=istartpar,iendpar
               len_scale(i,j,k)=akt*MAX(len_min,&
               len_scale(i,j,k)/(z1+len_scale(i,j,k)/l_scal))
         END DO
         END DO
      END DO       

!     Initialisierung der Felder fuer tke,tkvh,tkvm:

      IF (linitial) THEN  !nur beim allerersten Durchgang

         DO k=2,kem
            DO j=jstartpar,jendpar
            DO i=istartpar,iendpar

!              der Einfachheit halber nur lokale Berechnung der
!              vertikalen Gradienten:

               len=len_scale(i,j,k)

               IF (k.EQ.ke1) THEN
                  edh=z1/dicke(i,j,ke1)
               ELSE
                  edh=z2/(hhl(i,j,k+1)-hhl(i,j,k-1))
               END IF

!NEC_CB               DO n=1,nred
!NEC_CB                  grad(n)=(vari(i,j,k,n)-vari(i,j,k-1,n))*edh
!NEC_CB               END DO
!              ! nred=4 in all cases:
               ! thus these expressions are vectorized correctly
               grad(1)=(vari(i,j,k,1)-vari(i,j,k-1,1))*edh
               grad(2)=(vari(i,j,k,2)-vari(i,j,k-1,2))*edh
               grad(3)=(vari(i,j,k,3)-vari(i,j,k-1,3))*edh
               grad(4)=(vari(i,j,k,4)-vari(i,j,k-1,4))*edh

               fh2=a(i,j,k,4)*grad(tet_l)+a(i,j,k,5)*grad(h2o_g)
               fm2=grad(u_m)**2+grad(v_m)**2

               ! Vereinfachte Loesung mit Rf=Ri:
               IF (fh2.GE.(z1-rim)*fm2) THEN
                  ! Die krit. Ri-Zahl wird ueberschritten und lm, sowie lh
                  ! werden durch lm bei der krit. Ri-Zahl angenaehert:
                  fakt=z1/rim-z1
                  lm=len*(b_2-(a_6+a_3)*fakt)
                  lh=lm
               ELSE
                  fakt=fh2/(fm2-fh2)
                  lm=len*(b_2-(a_6+a_3)*fakt)
                  lh=len*(b_1-a_5*fakt)
               END IF

               wert=lm*fm2
               wert=MAX(wert-lh*fh2,rim*wert)

               q0vec(i,j)=MAX(vel_min,SQRT(d_m*len*wert))

               tkvm(i,j,k)=lm
               tkvh(i,j,k)=lh

!              Am Anfang konnte noch keine Advektion oder 
!              Diffusion von SQRT(2*TKE) berechnet werden:

               tketens(i,j,k)=z0
            END DO 
            END DO        
!NEC_CB Separated for vectorization
            DO n=1,ntlev
              DO j=jstartpar,jendpar
                DO i=istartpar,iendpar
                  tke(i,j,k,n)=q0vec(i,j)
                END DO
              END DO
            END DO
            DO n=1,ntlev
               DO j=jstartpar,jendpar
               DO i=istartpar,iendpar
                  tke(i,j,1,n)=tke(i,j,2,n)
               END DO
               END DO
            END DO
         END DO    
      ELSE
         DO k=2,kem
            DO j=jstartpar,jendpar
            DO i=istartpar,iendpar
               tkvh(i,j,k)=tkvh(i,j,k)/tke(i,j,k,nvor)
               tkvm(i,j,k)=tkvm(i,j,k)/tke(i,j,k,nvor)
            END DO
            END DO
         END DO    
      END IF

!     tkvh und tkvm enthalten jetzt die stabilitaetsabhaengigen 
!     Laengenmasse, nicht die Diffusionskoeffizienten!


! 2)  Berechnung der benoetigten vertikalen Gradienten und
!     Abspeichern auf vari():

!     Am unteren Modellrand:

      DO j=jstartpar,jendpar
      DO i=istartpar,iendpar
         hlp(i,j,ke1)=z1/dicke(i,j,ke1)
         vari(i,j,ke1,liq)=z0
      END DO
      END DO
 
      DO n=1,nred 
         DO j=jstartpar,jendpar
         DO i=istartpar,iendpar
            vari(i,j,ke1,n)=(vari(i,j,ke,n)-vari(i,j,ke1,n))*hlp(i,j,ke1)
         END DO
         END DO
      END DO   

!     An den darueberliegenden Nebenflaechen:

      IF (lnonloc) THEN

!        Berechnung nicht-lokaler Gradienten:

         DO n=1,ndiff

!           Berechnung vertikalen Integralfunktionen in hlp():
      
            DO j=jstartpar,jendpar
            DO i=istartpar,iendpar
               hlp(i,j,ke1)=z0
            END DO 
            END DO 

            DO k=ke,2,-1
               DO j=jstartpar,jendpar
               DO i=istartpar,iendpar
                  hlp(i,j,k)=hlp(i,j,k+1)+vari(i,j,k,n) &
                                         *(hhl(i,j,k)-hhl(i,j,k+1))
               END DO
               END DO
            END DO

            k1=1
            k2=2
            DO j=jstartpar,jendpar
            DO i=istartpar,iendpar
               lays(i,j,k1)=hhl(i,j,1)-hhl(i,j,2)
            END DO
            END DO
            DO k=2,ke   

!              Berechnung der nicht-lokalen Gradienten als mittlere
!              Differenzenquotienten ueber die stabilitaetsabhaengige
!              Laengenskala in tkvh() bzw tkvm():

!              Speichern der stab.abh. Laengenskala unter lay():

               IF (n.LE.nvel) THEN 
                  DO j=jstartpar,jendpar
                  DO i=istartpar,iendpar
                     lay(i,j)=tkvm(i,j,k)
                  END DO
                  END DO
               ELSE   
                  DO j=jstartpar,jendpar
                  DO i=istartpar,iendpar
                     lay(i,j)=tkvh(i,j,k)
                  END DO
                  END DO
               END IF
                  
!              Bestimmung der nicht-lokalen Gradienten und 
!              Zwischenspeichern derselben auf dem Feld dicke():

               DO j=jstartpar,jendpar
               DO i=istartpar,iendpar


                  lays(i,j,k2)=hhl(i,j,k)-hhl(i,j,k+1)

                  IF (lay(i,j).LE. &
                      z1d2*MIN(lays(i,j,k1),lays(i,j,k2))) THEN

!                    Die vertikalen Diffusionswege schneiden weder
!                    eine untere noch eine obere Hauptflaeche. Es 
!                    koennen somit die lokalen Gradienten genommen 
!                    werden. Bei sehr kleinen Diffusionslaengen, muss 
!                    aus num. Gruenden sogar die lokale Berechnung 
!                    gewaehlt werden:

                     dicke(i,j,k)=z2*(vari(i,j,k-1,n)-vari(i,j,k,n)) &
                                    /(hhl(i,j,k-1)-hhl(i,j,k+1))
                  ELSE

!                    Berechn. der benoetigten Referenzhoehen und -level:

                     h=hhl(i,j,k)
                     hu=MAX(h-lay(i,j),hhl(i,j,ke1))
                     hig(1)=hu+lay(i,j)
                     hig(2)=h+lay(i,j)
 
                     kk=k
111                  IF (hhl(i,j,kk).GT.hu) THEN
                        kk=kk+1
                        GOTO 111
                     END IF
                     ku=kk
                     DO ii=1,2
112                     IF (kk.GT.1) THEN
                           IF (hhl(i,j,kk).LE.hig(ii)) THEN
                              kk=kk-1
                              GOTO 112
                           END IF
                        END IF
                        lev(ii)=kk+1
                     END DO

!                    Berechnung der gemittelten Differenzenquotienten
!                    als Ausdruck fuer die nicht-lokalen Gradienten: 

                     wert=hlp(i,j,ku)-hlp(i,j,k) &
                         +hlp(i,j,lev(2))-hlp(i,j,lev(1)) &
                         +vari(i,j,ku-1,n)*(hu-hhl(i,j,ku)) &
                         -vari(i,j,lev(1)-1,n)*(hig(1)-hhl(i,j,lev(1)))&
                         +vari(i,j,lev(2)-1,n)*(hig(2)-hhl(i,j,lev(2)))

                     dicke(i,j,k)=wert/(lay(i,j)*(h-hu))
                  END IF   
               END DO    
               END DO    
               kk=k1
               k1=k2
               k2=kk

            END DO   

!           Sichern der nicht-lokalen Gradienten im Feld vari():

            DO k=2,ke
               DO j=jstartpar,jendpar
               DO i=istartpar,iendpar
                  vari(i,j,k,n)=dicke(i,j,k)
               END DO
               END DO
            END DO

         END DO

!        Belegung von dicke() mit den Schichtdicken
!        bzgl. Nebenflaechen:

         DO k=2,ke   
            DO j=jstartpar,jendpar
            DO i=istartpar,iendpar
               dicke(i,j,k)=(hhl(i,j,k-1)-hhl(i,j,k+1))*z1d2
            END DO
            END DO
         END DO      

      ELSE       

!        Berechnung lokaler Gradienten:
         
         DO k=ke,2,-1
            DO j=jstartpar,jendpar
            DO i=istartpar,iendpar
               dicke(i,j,k)=(hhl(i,j,k-1)-hhl(i,j,k+1))*z1d2
               hlp(i,j,k)=z1/dicke(i,j,k)
            END DO
            END DO
         END DO   

         DO n=1,ndiff
            DO k=ke,2,-1
               DO j=jstartpar,jendpar
               DO i=istartpar,iendpar
                  vari(i,j,k,n)=(vari(i,j,k-1,n)-vari(i,j,k,n))*hlp(i,j,k)
               END DO
               END DO
            END DO   
         END DO      

      END IF

!------------------------------------------------------------------------------------

     IF (l3dturb .OR. PRESENT(tket_hshr)) THEN
         !3D-turbulent shear or separate horizontal shear mode to be considered as
         !TKE source or at least to be calculated for output:

!        Berechnung der 3D-Korrektur des Scherungsterms der TKE incl. metr. Terme:

         DO k=2,ke

            !Neigung der Modellflaeche:
            DO j = jstart, jend
            DO i = istart, iend
               lay(i,j) = ( hhl(i+1,j,k) - hhl(i-1,j,k) ) * z1d2
               frc(i,j) = ( hhl(i,j+1,k) - hhl(i,j-1,k) ) * z1d2
            END DO
            END DO

            !Scherung durch horizontale Konfluenz:
            DO j = jstart, jend
            DO i = istart, iend
               zs11 = ( u(i,j,k  ) - u(i-1,j,k  )   +       &
                        u(i,j,k-1) - u(i-1,j,k-1) ) * z1d2
               zs22 = ( v(i,j,k  ) - v(i,j-1,k  )   +       &
                        v(i,j,k-1) - v(i,j-1,k-1) ) * z1d2

               zs11 = ( zs11 - lay(i,j)*vari(i,j,k,u_m) ) * eddlon * acrlat(j,1)
               zs22 = ( zs22 - frc(i,j)*vari(i,j,k,v_m) ) * edadlat

               hlp(i,j,k) = z2 * ( zs11**2 + zs22**2 )
            END DO
            END DO

            !Horizontale Scherung:
            DO j = jstart, jend
            DO i = istart-1, iend
               lays(i,j,1) = ( v(i  ,j,k  ) + v(i  ,j-1,k  ) +     &
                               v(i+1,j,k  ) + v(i+1,j-1,k  ) +     &
                               v(i  ,j,k-1) + v(i  ,j-1,k-1) +     &
                               v(i+1,j,k-1) + v(i+1,j-1,k-1) )     &
                               * 0.125_wp
            END DO
            END DO
            DO j = jstart-1, jend
            DO i = istart, iend
               lays(i,j,2) = ( u(i,j  ,k  ) + u(i-1,j  ,k  ) +     &
                               u(i,j+1,k  ) + u(i-1,j+1,k  ) +     &
                               u(i,j  ,k-1) + u(i-1,j  ,k-1) +     &
                               u(i,j+1,k-1) + u(i-1,j+1,k-1) )     &
                               * 0.125_wp
            END DO
            END DO
            DO j = jstart, jend
            DO i = istart, iend
               zs12 = lays(i,j,1) - lays(i-1,j,1)
               zs21 = lays(i,j,2) - lays(i,j-1,2)

               zs12 = ( zs12 - lay(i,j)*vari(i,j,k,v_m) ) * acrlat(j,1) * eddlon
               zs21 = ( zs21 - frc(i,j)*vari(i,j,k,u_m) ) * edadlat

               hlp(i,j,k) = hlp(i,j,k) + ( zs12 + zs21 )**2

               IF (l3dturb) THEN
                 ! tkvm was divided by tke above, so we need multiplications 
                 ! with tke to get the diffusion coefficients (isotropic part)
                 tkhm(i,j,k)=tkvm(i,j,k)*tke(i,j,k,nvor) ! for momentum
                 tkhh(i,j,k)=tkvh(i,j,k)*tke(i,j,k,nvor) ! for heat (scalar)
               END IF

            END DO
            END DO

            IF (itype_sher.GE.2 .OR. PRESENT(tket_hshr)) THEN
               !Separate horizontale Scherungsmode soll berechnet werden:

               DO j = jstart, jend
               DO i = istart, iend
                  ! related horizontal diffusion coefficient for momentum
                  src(i,j)=(a_hshr*l_hori)**2 * hlp(i,j,k)**z1d2 
               END DO
               END DO

               IF (PRESENT(tket_hshr)) THEN
                  DO j = jstart, jend
                  DO i = istart, iend
                     tket_hshr(i,j,k)=src(i,j)*hlp(i,j,k)
                  END DO
                  END DO
               END IF

               IF (itype_sher.GE.2 .AND. l3dturb) THEN
                 !Korrektur durch separate horizontal Scherungsmode:
                 DO j = jstart, jend
                 DO i = istart, iend
                   ! tkvm was divided by tke above, so we need multiplications 
                   ! with tke to get the diffusion coefficients
                   hlp(i,j,k) = hlp(i,j,k) + src(i,j)*hlp(i,j,k)/(tkvm(i,j,k)*tke(i,j,k,nvor))
                 END DO
                 END DO
                 IF (itype_sher.EQ.3) THEN
                    !Inlcuding related horizontal diffusion:
                    DO j = jstart, jend
                    DO i = istart, iend
                       ! total horizontal diffusion coefficient for
                       tkhm(i,j,k)= tkhm(i,j,k)+          src(i,j) ! momentum
                       tkhh(i,j,k)= tkhh(i,j,k)+(b_1/b_2)*src(i,j) ! heat (scalar)
                    END DO
                    END DO
                 END IF
               END IF
            END IF

            !Vertikale Scherungskorrektur:
            DO j = jstart, jend
            DO i = istart-1, iend
               lays(i,j,1) = z1d2 * ( w(i,j,k) + w(i+1,j,k) )
            END DO
            END DO
            DO j = jstart-1, jend
            DO i = istart, iend
               lays(i,j,2) = z1d2 * ( w(i,j,k) + w(i,j+1,k) )
            END DO
            END DO
            DO j = jstart, jend
            DO i = istart, iend
               zs13 = lays(i,j,1) - lays(i-1,j,1)
               zs23 = lays(i,j,2) - lays(i,j-1,2)

               zs33 = ( w(i,j,k+1) - w(i,j,k-1) ) &
                        / ( hhl(i,j,k-1) - hhl(i,j,k+1) )

               zs13  = ( zs13 - lay(i,j)*zs33 ) * acrlat(j,1) * eddlon
               zs23  = ( zs23 - frc(i,j)*zs33 ) * edadlat

               hlp(i,j,k) = hlp(i,j,k) + zs13 * ( z2*vari(i,j,k,u_m) + zs13) &
                                       + zs23 * ( z2*vari(i,j,k,v_m) + zs23) &
                                       + z2 * zs33**2
            END DO
            END DO

         END DO

         !Neigungskorrektur der vertikalen Scherung am unteren Modellrand
         !bei 3D-Scherung :

         DO j = jstart, jend
         DO i = istart, iend
            !Neigung der Erdoberflaeche:
            lay(i,j) = ( hhl(i+1,j,ke1) - hhl(i-1,j,ke1) )*z1d2
            frc(i,j) = ( hhl(i,j+1,ke1) - hhl(i,j-1,ke1) )*z1d2

            zs33 = w(i,j,ke) / hhl(i,j,ke)
            zs13 = lay(i,j)*zs33
            zs23 = frc(i,j)*zs33

            hlp(i,j,ke1) = ( frc(i,j)*vari(i,j,ke1,u_m) - lay(i,j)*vari(i,j,ke1,v_m) )**2 &
                         + zs13 * ( z2*vari(i,j,ke1,u_m) + zs13 ) &
                         + zs23 * ( z2*vari(i,j,ke1,v_m) + zs23 )
         END DO
         END DO

         !Beachte: Von der gesamten 3D-Scherung wurde die reine vertikale Scherung
         !         des Horizontalwindes ausgespart. Diese wird i.f. bestimmt.

      END IF

!__COSMO__---------------------------------------------------------------------------

! 3)  Hauptschleife: Bestimmung der TKE und der Stabilitaetsfunkt.:

      DO it_durch=it_start,it_end

         IF (kcm.LE.kem) THEN

!           Es gibt einen Rauhigkeitsbestand:

            DO j=jstartpar,jendpar
            DO i=istartpar,iendpar

!              Belegung mit Werten gueltig ausserhalb des Bestandes:

               dd(i,j,0)=d_m

               dd(i,j,1)=d_1
               dd(i,j,2)=d_2
               dd(i,j,3)=d_3
               dd(i,j,4)=d_4
               dd(i,j,5)=d_5
               dd(i,j,6)=d_6

               dd(i,j,7)=rim
            END DO
            END DO

         END IF

         DO k=2,kem

!           Berechnung der atmosphaerischen Vertikal-Antriebe:

            DO j=jstartpar,jendpar
            DO i=istartpar,iendpar

               lays(i,j,1)=vari(i,j,k,u_m)**2+vari(i,j,k,v_m)**2
               lays(i,j,2)=a(i,j,k,4)*vari(i,j,k,tet_l) &
                          +a(i,j,k,5)*vari(i,j,k,h2o_g)
 
               wert=tkvm(i,j,k)*lays(i,j,1)

               frc(i,j)=MAX(dd(i,j,7)*wert,wert-tkvh(i,j,k)*lays(i,j,2))

!              Durch die MAX-Funktion wird frc so nach unten 
!              beschraenkt, dass im level-2 Gleichgewicht die
!              krit. Rizahl nie ueberschritten wird. In diesem
!              Gleichgewicht entsrpicht dies dann einem von der
!              Windscherung abhaengigen Minimalwert von TKE.

!              Die Stabilitaetsfunktionen werden aber unabhaengig von
!              dieser Beschraenkung berechnet!

!              Stabilitaetskorrektur der turbulenten Laengenskala
!              bei stabilier Schichtung:

               vel=tke(i,j,k,nvor)
               wert=a_stab*SQRT(MAX(z0,lays(i,j,2)))

               len_scale(i,j,k)=vel*len_scale(i,j,k)/(vel+wert*len_scale(i,j,k))

            END DO
            END DO

!           Berechnung des Antriebs durch Nachlaufproduktion:

            IF (k.GE.kcm) THEN !Innerhalb des Bestandes:

               IF (k.EQ.ke1) THEN
                  DO j=jstartpar,jendpar
                  DO i=istartpar,iendpar
                     lay(i,j)=SQRT(wind(i,j,ke1,u_m)**2 &
                                  +wind(i,j,ke1,v_m)**2)
                  END DO
                  END DO
               ELSE
                  DO j=jstartpar,jendpar
                  DO i=istartpar,iendpar
                     dpu=dp0(i,j,k)
                     dpo=dp0(i,j,k-1)
                     dp2=dpu+dpo

                     vl1=(wind(i,j,k,u_m)*dpo+wind(i,j,k-1,u_m)*dpu)/dp2
                     vl2=(wind(i,j,k,v_m)*dpo+wind(i,j,k-1,v_m)*dpu)/dp2
                     lay(i,j)=SQRT(vl1**2+vl2**2+w(i,j,k)**2)
                  END DO
                  END DO
               END If

               DO j=jstartpar,jendpar
               DO i=istartpar,iendpar

!                 Hilfsgroesse zur Berechn. der reduzierten Konstanten:

                  wert=z3*c_sml(i,j,k)*len_scale(i,j,k) &
                                      *lay(i,j)/tke(i,j,k,nvor)

!                 Berechnung der modifizierten Modellparameter:

                  dd(i,j,0)=d_m/(z1+d_m*wert)

                  dd(i,j,1)=a_h/(z1+a_h*wert)
                  dd(i,j,2)=a_m/(z1+z2*a_m*wert)

                  dd(i,j,3)=z9*dd(i,j,1)
                  dd(i,j,4)=z6*dd(i,j,2)
                  dd(i,j,5)=z3*(d_h+dd(i,j,4))
                  dd(i,j,6)=dd(i,j,3)+z3*dd(i,j,4)
                  dd(i,j,1)=z1/dd(i,j,1)
                  dd(i,j,2)=z1/dd(i,j,2)

                  dd(i,j,7)=z1/(z1+(dd(i,j,0)-dd(i,j,4))/dd(i,j,5))

!                 TKE-Forcing durch Addition der Nachlaufproduktion zum vertikalen Forcing:

                  lay(i,j)=frc(i,j)+c_big(i,j,k)*lay(i,j)**3/tke(i,j,k,nvor)

               END DO
               END DO

            ELSE !Ausserhalb des Rauhigkeitsbestandes

!              TKE-Forcing gleich vertikales Forcing:

               DO j=jstartpar,jendpar
               DO i=istartpar,iendpar
                  lay(i,j)=frc(i,j)
               END DO 
               END DO 

            END IF

            IF (.NOT.(linitial.OR.k.eq.ke1)) THEN !nicht bei der Initialisierung oder am Unterrand

               IF (lsso .AND. ltkesso .AND. PRESENT(ut_sso) .AND. PRESENT(vt_sso)) THEN
                  !SSO-Schema ist aktiv und SSO-Tendenzen des Windes sind vorhanden:

!                 Berechnung der TKE-Tendenz durch Nachlaufwirbelproduktion aus SSO_Tendenzen:               

                  DO j=jstart,jend
                  DO i=istart,iend
                     dpu=dp0(i,j,k)
                     dpo=dp0(i,j,k-1)
                     dp2=dpu+dpo

                     vl1=-(ut_sso(i,j,k)  *wind(i,j,k  ,u_m) + vt_sso(i,j,k)  *wind(i,j,k  ,v_m))*dpo
                     vl2=-(ut_sso(i,j,k-1)*wind(i,j,k-1,u_m) + vt_sso(i,j,k-1)*wind(i,j,k-1,v_m))*dpu

                     !US: the two fields are only to ensure for the moment bit-identical 
                     !    results to the former version 4.25
                     src(i,j)=  MAX( z0, (vl1+vl2)/ dp2                  )
                     src2(i,j)= MAX( z0, (vl1+vl2)/(dp2*tke(i,j,k,nvor)) )

!                    Beachte:
!                    Die SSO-Tendenzen beziehen sich tatsaechlich auf Massenpunkte, sie werden
!                    erst spaeter in SUB 'organize_physics' auf die Zwischenpositionen interpoliert!
!                    Obwohl vl1 und vl2 immer positiv sein muessten, wird zur Sicherheit MAX benutzt!
                  END DO
                  END DO

                  IF (PRESENT(tket_sso)) THEN
                     DO j=jstart,jend
                     DO i=istart,iend
                        tket_sso(i,j,k)=src(i,j)
                     END DO
                     END DO
                  END IF

                  IF (ltkesso) THEN !Nachlaufwirbeltendenzen sollen beruecksichtigt werden
                     DO j=jstart,jend
                     DO i=istart,iend
                        lay(i,j)=lay(i,j) + src2(i,j)
                     END DO
                     END DO
                  END IF

               END IF

               IF (lconv .AND. ltkecon .AND. PRESENT(tket_conv)) THEN
                  !Konvektionsschema ist aktiv, soll mit Turbulenz interagieren und conv. TKE-Tend. ist vorhanden:

!                 Addition die TKE-Quelle durch konvektive Aktivitaet:

                  DO j=jstart,jend
                  DO i=istart,iend
                     lay(i,j)=lay(i,j) + MAX( z0, tket_conv(i,j,k)/tke(i,j,k,nvor) )
                  END DO
                  END DO

!                 Beachte:  Obwohl tket_conv immer positiv sein muesste, wird zur Sicherheit MAX benutzt!
               END IF

            END IF

!           Beruecksichtigung der 3D_Scherungskorrektur:

            IF (itype_sher.GT.1) THEN !3D-turbulence or separate horizontal shear

!              Addition der 3D-Scherungskorrektur zum bisherigen TKE-Forcing:

               DO j = jstart, jend
               DO i = istart, iend
                  lay(i,j)=lay(i,j)+tkvm(i,j,k)*hlp(i,j,k)
               END DO
               END DO

            END IF

!           Berechnung der neuen TKE-Werte:

            DO j=jstartpar,jendpar
            DO i=istartpar,iendpar

!              Berechnung einiger Hilfsgroessen:

               q0=tke(i,j,k,nvor) + tketens_adv(i,j,k)*dt_tke
               l_turb=len_scale(i,j,k)
               ld=dd(i,j,0)*l_turb

!              q0 wird explizit advektionskorrigiert, damit die untenstehende
!              zeitliche Glaettung (tkesmot) nicht die Transportgeschindigkeit
!              der Advektion beeinflusst (verlangsamt), aber trotzdem
!              verwendet werden kann, um Instabilitaeten durch andere Prozesse
!              in der Zeitintegration zu unterdrücken.

!              SQRT(2*TKE)-Prognose:

               q1=ld/dt_tke
               q2=MAX(z0,q0+tketens(i,j,k)*dt_tke)+lay(i,j)*dt_tke

!              In tketens steht die zuvor berechnete 
!              Diffusionstendenz (einschl. Drucktransport).
!              Die MAX-Funkt. verhindert, dass Transportterme die TKE
!              negativ machen koennen.

               q3=q1*(sqrt(z1+z4*q2/q1)-z1)*z1d2

               q1=SQRT(l_turb*frc(i,j)*dd(i,j,4)/(z1-c_g))
               tke(i,j,k,ntur)=MAX(vel_min,q1,q3*(z1-tkesmot)+q0*tkesmot)

!              Die MAX-Funktion und die zeitliche Glaettung
!              sind (vor allem wegen der Wirkung des Zirkulationstermes)
!              nur in Ausnahmefaellen noetig und koennen dann
!              vermeiden, dass die TKE zu klein wird.
!              q1 ist ein Minimalwert, um eine moegliche Singuaritaet bei der
!              Berechnung der Stabilitaetsfunktion fuer labile Schichtung
!              zu vermeiden.

!              Sichern des thermischen TKE-Antriebs:

               hlp(i,j,k)=lays(i,j,2)

!              Sichern der TKE-Dissipation u.a. als thermische Quelle:
              !edr(i,j,k)=q3**3/ld
               edr(i,j,k)=tke(i,j,k,ntur)**3/ld !mit gefiltertem tke-Feld

            END DO    
            END DO    

!           Berechnung der neuen stabilitaetsabhangigen Laengenskalen:

            IF (lstfnct) THEN

!US            CALL stab_funct_old(lays(:,:,1), lays(:,:,2), frc, &
!US                            istartpar,iendpar, jstartpar,jendpar, k, ntur)

               CALL stab_funct(sm=tkvm(:,:,k), sh=tkvh(:,:,k), fm2=lays(:,:,1), fh2=lays(:,:,2), &
                               frc=frc, tvs=tke(:,:,k,ntur), tls=len_scale(:,:,k), dd=dd,        &
                               i_st=istartpar,i_en=iendpar, j_st=jstartpar,j_en=jendpar)


               DO j=jstartpar,jendpar
               DO i=istartpar,iendpar
                  tkvm(i,j,k)=len_scale(i,j,k)*tkvm(i,j,k)
                  tkvh(i,j,k)=len_scale(i,j,k)*tkvh(i,j,k)
               END DO
               END DO

            END IF

!SCLM------------------------------------------------------------------
#ifdef SCLM
            IF (lscm .AND. it_durch.EQ.it_end) THEN

               q3=tke(im,jm,k,ntur); q2=q3**2

               fakt=len_scale(im,jm,k)/q3; x0=dd(im,jm,0)

               km=tkvm(im,jm,k)*q3; kh=tkvh(im,jm,k)*q3

               flux(u_m)  =-km*vari(im,jm,k,u_m)
               flux(v_m)  =-km*vari(im,jm,k,v_m)
               flux(tet_l)=-kh*vari(im,jm,k,tet_l)

               x1=flux(u_m)*vari(im,jm,k,u_m)
               x2=flux(v_m)*vari(im,jm,k,v_m)
               x3=-kh*lays(im,jm,2)

               !Achtung: TKE%mod(k)%val zeigt z.Z. noch auf tke(im,jm,k,nvor),
               !somit wird also der alte tke-Wert mit dem neuen ueberschrieben,
               !was aber ohne Bedeutung ist, weil ab jetzt der alte tke-Wert nicht
               !mehr benoetigt wird. Beachte, dass die Modelleinheit hier [m/s] ist.

               TKE_SCLM%mod(k)%val=q3; TKE_SCLM%mod(k)%vst=i_cal

               BOYPR%mod(k)%val=x3      ; BOYPR%mod(k)%vst=i_cal
               SHRPR%mod(k)%val=-(x1+x2); SHRPR%mod(k)%vst=i_cal
               DISSI%mod(k)%val=q3**3/(x0*len_scale(im,jm,k))
                                          DISSI%mod(k)%vst=i_cal
               TRANP%mod(k)%val=q3*tketens(im,jm,k)
                                          TRANP%mod(k)%vst=i_cal

               varc(u_m)  =z1d3*q2+x0*fakt*(-z4*x1+z2*x2-z2*x3)
               varc(v_m)  =z1d3*q2+x0*fakt*(+z2*x1-z4*x2-z2*x3)
               varc(tet_l)=-d_h*fakt*flux(tet_l)*vari(im,jm,k,tet_l)

               UWA%mod(k)%val=flux(u_m)  ; UWA%mod(k)%vst=i_cal
               VWA%mod(k)%val=flux(v_m)  ; VWA%mod(k)%vst=i_cal
               TWA%mod(k)%val=flux(tet_l); TWA%mod(k)%vst=i_cal
               UST%mod(k)%val=SQRT(km*SQRT(lays(im,jm,1)))
                                           UST%mod(k)%vst=i_cal

               IF (varc(u_m).GE.z0 .AND. varc(v_m).GE.z0 .AND. &
                   varc(u_m)+varc(v_m).LE.q2) THEN

                  UUA%mod(k)%val=varc(u_m); UUA%mod(k)%vst=i_cal
                  VVA%mod(k)%val=varc(v_m); VVA%mod(k)%vst=i_cal
                  WWA%mod(k)%val=q2-UUA%mod(k)%val-VVA%mod(k)%val
                                            WWA%mod(k)%vst=i_cal
               END IF

               TTA%mod(k)%val=varc(tet_l); TTA%mod(k)%vst=i_cal

            END IF   
#endif
!SCLM------------------------------------------------------------------

            IF (ltmpcor) THEN

!              Sichern der TKE-Dissipation als thermische Quelle:

               DO j=jstartpar,jendpar
               DO i=istartpar,iendpar
                  tketens(i,j,k)=tke(i,j,k,ntur)**3/(dd(i,j,0)*len_scale(i,j,k))
               END DO
               END DO

            END IF

         END DO 

         IF (it_durch.LT.it_end) THEN

!           Permutation der TKE-Zeitebenen als Vorbereitung des 
!           naechsten Iterationsschrittes:

            DO n=1,ntlev
               tkestep(n)=mod(ntur+n,ntlev)+1
            END DO

            nvor=tkestep(ntlev-1)
            ntur=tkestep(ntlev)

            IF (itype_tran.EQ.2) THEN

!              tke(i,j,ke1,ntur) wurde schon in turbtran bestimmt
!              und auf der Zeitstufe ntur vor der Permutation
!              abgelegt. Dies ist jetzt die Zeitstufe nvor:

               DO j=jstartpar,jendpar
               DO i=istartpar,iendpar

                  tke(i,j,ke1,ntur)=tke(i,j,ke1,nvor)

               END DO
               END DO
            END IF
         END IF   

      END DO !Iterationen ueber it_durch

      DO j=jstartpar,jendpar
      DO i=istartpar,iendpar
         tke(i,j,1,ntur)=tke(i,j,2,ntur)
      END DO
      END DO

!  4) Berechnung der effektiven turbulenten Vertikalgradienten
!     und weiterer Turbulenzgroessen:

         DO k=2,kem

            IF (ltmpcor) THEN

!              Berechnung des vert. Temp.grad. fuer den Phasendiffusionsterm:
!              
               DO j=jstartpar,jendpar
               DO i=istartpar,iendpar
                  lay(i,j)=a(i,j,k,2)*vari(i,j,k,tet_l)-tet_g &
                          +zlhocp*vari(i,j,k,liq) ! vert. Temperaturgradient
               END DO   
               END DO   

!              Dies geschieht schon hier, weil im naechsten Schritt das Feld vari()
!              durch die effiktiven Gradienten ueberschrieben wird.
            END IF    

            DO j=jstartpar,jendpar
            DO i=istartpar,iendpar

               x2=a(i,j,k,3)*a(i,j,k,2)

               rc_2d(i,j) = (c_scld * rcld(i,j,k)) / (z1+rcld(i,j,k)*(c_scld-z1))

               rcld(i,j,k)=sqrt(len_scale(i,j,k)*tkvh(i,j,k)*d_h)* &
                           abs(x2*vari(i,j,k,tet_l)-vari(i,j,k,h2o_g))

!              Unter rcld steht jetzt die (geschaetzte) Standardabw.
!              des Saettigungsdefizites!

!              Berechn. der Diffusionskoeffizienten:

               tkvh(i,j,k)=MAX(con_h,tkhmin,tkvh(i,j,k)*tke(i,j,k,ntur))
               tkvm(i,j,k)=MAX(con_m,tkmmin,tkvm(i,j,k)*tke(i,j,k,ntur))

!              tkvh und tkvm enthalten jetzt nicht mehr Diffusions-
!              laengen, sondern Diffusionskoeffizienten in m^2/s!
            END DO
            END DO

            DO j=jstartpar,jendpar
            DO i=istartpar,iendpar

!              Berech. der thermodynamischen Flusskonversionsmatrix
!              und der Standardabweichnung der Verteilung des
!              Saettigungsdefizites:

               x1=z1/(z1+zlhocp*a(i,j,k,3))
               x2=a(i,j,k,3)*a(i,j,k,2)

               a(i,j,k,3)=rc_2d(i,j)

               flukon33         =z1-rc_2d(i,j)*(z1-x1)
               flukon34         =x1*zlhocp/a(i,j,k,2)*rc_2d(i,j)

               flukon53         =-x1*x2*rc_2d(i,j)
               flukon54         =x1*rc_2d(i,j)
               flukon55         =z1-liqfak

               flukon43         =-flukon53
               flukon44         =z1-flukon54

               flux_3   =a(i,j,k,2)*a(i,j,k,1) & !ex*Cp/Cpd
                        *(flukon33         *vari(i,j,k,tet_l) &
                         +flukon34         *vari(i,j,k,h2o_g))
               flux_4   =(flukon43         *vari(i,j,k,tet_l) &
                         +flukon44         *vari(i,j,k,h2o_g))
               flux_5   =(flukon53         *vari(i,j,k,tet_l) &
                         +flukon54         *vari(i,j,k,h2o_g) &
                         +flukon55         *vari(i,j,k,liq)) 

               vari(i,j,k,tem)=flux_3 !ex*(Cp/Cpd)*grad(tet)
               vari(i,j,k,vap)=flux_4
               vari(i,j,k,liq)=flux_5

               !Achtung: 
               !'vari(i,j,k,tem)' ist der in der Temperaturgleichung benoetigte
               !effective Gradient und entspricht ex*(Cp/Cpd)*tet_flux_dens.
               !Die Indices 'tet', 'tem' und 'tet_l' haben den gleichen Wert
               !und sollen nur den jeweiligen Inhalt verdeutlichen.

               !Im Falle "icldm_turb.NE.-1" ist "liqfak.EQ.z1" und 'flukon55' verschwindet.
               !Im Falle "icldm_turb.EQ.-1" dagegen verschwinden 'flukon53' und 'flukon54'.
               !Dann ist aber "flukon55.EQ.z1". Ferner verschwinden in diesem Fall auch
               !'flukon34' und 'flukon43', wobei die physikalischen Groesen 'tet_l' und 'tet'
               !sowie 'h2o_g' und 'vap' identisch sind.
               
!-----------------------------------------------------------------------
!     nur im alten 1-d_Modell:

!              Berechnung der der Flussdichten fuer Impuls,
!              fuehlbare und latente Energie (zum Erdboden hin posit.): 

!              km=-rhon(i,j,k)*tkvm(i,j,k)
!              kh=-rhon(i,j,k)*tkvh(i,j,k)

!              flimpuls(i,j,k)=km*sqrt(vari(i,j,k,u_s)**2  &
!                                      +vari(i,j,k,v_s)**2)
!              flfbwae(i,j,k)=cp_d*kh*vari(i,j,k,tem)
!              fllawae(i,j,k)=lh_v*kh*vari(i,j,k,vap)

!        Beachte: flimpuls,flfbwae,fllawael
!                 werden nur fuer Diagnosezwecke berechnet.
!-----------------------------------------------------------------------

            END DO
            END DO
           
            IF (ltmpcor) THEN

!              Berechnung der Temperaturtendenzen durch TKE-Quellen
!              auf Nebenflaechen (ausser der Divergenz des 
!              Drucktransportes):

               DO j=jstartpar,jendpar
               DO i=istartpar,iendpar
                  thermik=tkvh(i,j,k)*hlp(i,j,k)
                  phasdif=tkvh(i,j,k)*lay(i,j) &
                         *(zrcpv*vari(i,j,k,vap)+zrcpl*vari(i,j,k,liq))  

                  tketens(i,j,k)=len_scale(i,j,k)/a(i,j,k,1) &
                          *((tketens(i,j,k)+thermik)/cp_d+phasdif)
               END DO   
               END DO   

!              Beachte:
!              In tketens() steht bisher die TKE-Dissipation.
!              Wegen der spaeteren Interpolation auf Hauptflaechen,
!              wird mit der turbulenten Laengenskala multipliziert.
            END IF
   
         END DO


!  5) Untere Randwerte der Turbulenzgroessen:

         IF (kem.EQ.ke1) THEN

            DO j=jstartpar,jendpar
            DO i=istartpar,iendpar
 
!              Dicke der laminaren Grenzschicht:
               l_lam=z0m(i,j)*SQRT(con_m/tkvm(i,j,ke1))

!              Reduktionsfaktor:
               fakt=(l_lam/dicke(i,j,ke1))*(tkvh(i,j,ke1)/con_h)
               tfh(i,j)=z1/(z1+rlam_heat*fakt)
               tfm(i,j)=z1/(z1+rlam_mom *fakt)

               vel=MAX(vel_min, &
                       SQRT(wind(i,j,ke,u_m)**2+wind(i,j,ke,v_m)**2))

!              Neue Transferkoeffizienten:

               fakt=z1/(vel*dicke(i,j,ke1))
               tcm(i,j)=fakt*tkvm(i,j,ke1)*tfm(i,j)
               tch(i,j)=fakt*tkvh(i,j,ke1)*tfh(i,j)

!              Berechnung der neuen Rauhigkeitslaenge ueber Meer:

               IF (fr_land(i,j).LT.z1d2) THEN

                  ! Use ice surface roughness or open-water surface roughness
                  ! according to lo_ice
                  IF ( lo_ice(i,j) ) THEN
                     gz0(i,j)=grav*z0_ice
                  ELSE
!                    Berechnung der Schubspannung:

                     wert=tcm(i,j)*vel*SQRT(vel**2+tke(i,j,ke1,ntur)**2)
                     
!                    grav*z0 mittels Charnock-Formel: 

                     gz0(i,j)=MAX(grav*len_min,alpha0*wert+alpha1*grav*con_m/SQRT(wert))
                  END IF
               END IF

            END DO   
            END DO   

         ELSE

!           Untere Randwerte der Turbulenzgroessen mit hilfe der
!           bereits vorher (in turbtran oder partur(b,s) berechneten
!           Transferkoeffizienten:

            IF (itype_tran.EQ.1) THEN

!              Es wurde noch keine unteren Randwerte fuer 
!              tke(), sowie tkvm(), tkvh() berechnet:

               DO j=jstartpar,jendpar
               DO i=istartpar,iendpar
                  vel=MAX(vel_min, &
                          SQRT(wind(i,j,ke,u_m)**2+wind(i,j,ke,v_m)**2))
                  wert=vel*dicke(i,j,ke1)
                  tkvh(i,j,ke1)=wert*tch(i,j)
                  tkvm(i,j,ke1)=wert*tcm(i,j)

!                 Weil in partur(b,s) keine laminare Grenzschicht
!                 beruecksichtigt wird, sind die Reduktionskoeff.
!                 immer gleich 1:

!                 tfh(i,j)=z1
!                 tfm(i,j)=z1

!                 Der SQRT(2*TKE)-Wert am unteren Modellrand wird mit 
!                 Hilfe von c_tke*Ustar bestimmt: 

                  tke(i,j,ke1,ntur)=MAX(vel_min, c_tke*vel*SQRT(tcm(i,j)))
               END DO   
               END DO   

            END IF
   
            DO j=jstartpar,jendpar
            DO i=istartpar,iendpar

!              Thermischer Antrieb (fh2):

               hlp(i,j,ke1)=a(i,j,ke1,4)*vari(i,j,ke1,tet_l) &
                           +a(i,j,ke1,5)*vari(i,j,ke1,h2o_g)

!              Thermodynamische Korrektur des effektiven vertikalen 
!              Gradienten der pot. Temperatur am unteren Modellrand:

               vari(i,j,ke1,tem)=vari(i,j,ke1,tem) & !grad(tet)
                                *a(i,j,ke1,2)*a(i,j,ke1,1) !ex*Cp/Cpd

               !Beachte:
               !Am Unterrand soll grad(liq) = 0 sein. Daher ist dort
               !grad(tet_l) = grad(tet).
               !'vari(i,j,ke1,tem)' ist der effective in der Temperaturgleichung
               !benoetigte Gradient und entspricht ex*(Cp/Cpd)*tet_flux_dens.

               a(i,j,ke1,3)=rcld(i,j,ke1)
            END DO   
            END DO   

            IF (ltmpcor) THEN  

!              Bestimmung der unteren Randwerte der Enthalpieproduktion
!              durch Wechselwirkung mit der TKE (diss und thermik),
!              sowie durch Phasendiffusion:
          
               DO j=jstartpar,jendpar
               DO i=istartpar,iendpar

                  wert=(t(i,j,ke)-tp(i,j))/dicke(i,j,ke1) !T-Grad.
                  thermik=tkvh(i,j,ke1)*hlp(i,j,ke1)
                  phasdif=tkvh(i,j,ke1)*wert &
                         *(zrcpv*vari(i,j,ke1,vap)+zrcpl*vari(i,j,ke1,liq))

                  wert=tke(i,j,ke1,ntur)**3/(d_m*len_scale(i,j,ke1))
                  tketens(i,j,ke1)=len_scale(i,j,ke1)/a(i,j,ke1,1) &
                         *((wert+thermik)/cp_d+phasdif)

!test             tketens(i,j,ke1)=tketens(i,j,ke)

               END DO   
               END DO   
                
            END IF    
         
         END IF 

#ifdef tkvfilter
! SPo: vertical filtering for diffusion coefficients
         IF (ltkvfilter) THEN
            DO k=2,ke
               IF (k==2 .OR. k==ke) THEN

                  tkvh_filter(:,:,k) = 0.5*tkvh(:,:,k) + 0.25*(tkvh(:,:,k-1)+tkvh(:,:,k+1))
                  tkvm_filter(:,:,k) = 0.5*tkvm(:,:,k) + 0.25*(tkvm(:,:,k-1)+tkvm(:,:,k+1))

               ELSE

                  tkvh_filter(:,:,k) = 0.5*tkvh(:,:,k) + 0.2*(tkvh(:,:,k-1)+tkvh(:,:,k+1)) + 0.05*(tkvh(:,:,k-2)+tkvh(:,:,k+2))
                  tkvm_filter(:,:,k) = 0.5*tkvm(:,:,k) + 0.2*(tkvm(:,:,k-1)+tkvm(:,:,k+1)) + 0.05*(tkvm(:,:,k-2)+tkvm(:,:,k+2))

               END IF
            END DO

! no filter of lowermost layer to persits energy/mass conservation
!             tkvh_filter(:,:,ke1) = 0.9*tkvh(:,:,ke1) + 0.1*tkvh(:,:,ke)
!             tkvm_filter(:,:,ke1) = 0.9*tkvm(:,:,ke1) + 0.1*tkvm(:,:,ke)

            ! overwrite unfiltered with filtered values
            tkvh(:,:,2:ke) = tkvh_filter(:,:,2:ke)
            tkvm(:,:,2:ke) = tkvm_filter(:,:,2:ke)
         END IF
#endif

! 6)  Berechnung der zu TKE-Quellen gehoerigen Temperatur-
!     tendenzen ausser der Divergenz des Drucktransportes:

      IF (ltmpcor) THEN   
         DO j=jstart,jend
         DO i=istart,iend
            t_tens(i,j,1)=t_tens(i,j,1) &
                        +tketens(i,j,2) &
                         /(len_scale(i,j,1)+len_scale(i,j,2))
         END DO
         END DO
         DO k=2,ke
            DO j=jstart,jend
            DO i=istart,iend
               t_tens(i,j,k)=t_tens(i,j,k) &
                           +(tketens(i,j,k)+tketens(i,j,k+1)) &
                            /(len_scale(i,j,k)+len_scale(i,j,k+1))
            END DO
            END DO
         END DO
      END IF   


! 7)  Bestimmung des Drucktransporttermes:

      IF (pat_len.GT.z0) THEN

         DO j=jstartpar,jendpar
         DO i=istartpar,iendpar
            l_pat(i,j)=l_hori*d_pat(i,j)/(l_hori+d_pat(i,j))
         END DO
         END DO

         DO k=2,ke1

            IF (k.LT.ke1) THEN
!              Interpolation des Druckes auf die naechst hoehere 
!              Nebenflaeche:

               DO j=jstartpar,jendpar
               DO i=istartpar,iendpar
                  lay(i,j)=((prs(i,j,k  ))*dp0(i,j,k-1)&
                           +(prs(i,j,k-1))*dp0(i,j,k  ))&
                          /(dp0(i,j,k)+dp0(i,j,k-1))
               END DO   
               END DO   
            ELSE
               DO j=jstartpar,jendpar
               DO i=istartpar,iendpar
                  lay(i,j)=pr(i,j)
               END DO
               END DO
            END IF

            DO j=jstartpar,jendpar
            DO i=istartpar,iendpar

               fakt=z1-z2*ABS(a(i,j,k,3)-z1d2)
               len=MAX(l_pat(i,j),SQRT(fakt*len_scale(i,j,k)*l_hori))

!              Berechnung der lokalen Kohaerenzlaenge fuer die 
!              Drucktransport-Parameterisierung:
!              kohae_len=l_pat*ex_fakt*grad(tet_v)*r/g:

               fakt=hlp(i,j,k)*lay(i,j)/(rhon(i,j,k)*grav**2)
               kohae_len=len*SIGN(z1,fakt)*MIN(ABS(fakt),z1)

!              Belegung von tketens mit dem zur vert. turb. 
!              TKE-Flussdichte mittels Drucktransport gehoerigen
!              TKE-Gradienten:

               tketens(i,j,k)=kohae_len*hlp(i,j,k)

!              Die Divergenz des zugehoerigen Flusses ist gleichzeitig
!              eine weitere Quelle fuer thermische Energie.

            END DO
            END DO

         END DO   

      ELSE

         DO k=2,ke1
            DO j=jstartpar,jendpar
            DO i=istartpar,iendpar
               tketens(i,j,k)=z0
            END DO
            END DO
         END DO   
         
      END IF   


! 8)  Berechnung der korrigierten Flussdichten, die fuer die 
!     expliziten Berechnungen (Korrekturen) benoetigt werden 
!     und Ablegen derselben im Feld a(), sowie weiter Vorbereitungen:

      IF (lexplizit.or.limplizit) THEN

         wert=securi/dt_var

         DO j=jstartpar,jendpar
         DO i=istartpar,iendpar
            lay(i,j)=(rhon(i,j,ke)+rhon(i,j,ke1)) &
                    *(hhl(i,j,ke)-hhl(i,j,ke1))*z1d2
         END DO
         END DO

         IF (nk2.LE.ndiff) THEN
!           Berechnung der korrigierten Gradienten der thermodyn.
!           (skalaren) Variablen:

!           Beachte:
!           Cp/Cpd bleibt noch unter a(i,j,k,1) gespeichert.

            IF (imode_turb.EQ.2) THEN
!              Die gesamte Rechnung ist explizit, daher werden die
!              unkorrigierten effektiven Gradienten benutzt:

               DO n=nk2,ndiff
                  DO k=2,ke1
                     DO j=jstartpar,jendpar
                     DO i=istartpar,iendpar
                        a(i,j,k,n)=vari(i,j,k,n)
                     END DO
                     END DO
                  END DO
               END DO
            ELSE   
!              Es wird ein Teil der Flussdichtedivergenzen implizit
!              berechnet. Daher werden nur korrigierte Gradienten
!              benutzt, so dass das explizite Verfahren nur den bisher 
!              nicht berechneten Flussdivergenzanteil uebernimmt:

!              Am unteren Rand ist keine explizite Korrektur noetig,
!              weil grad(liq)=0 ist, und hier kein Unterschied 
!              zwischen lokaler und nicht-lokaler Gradientberechnung
!              besteht:

               DO n=nk2,ndiff
                  DO j=jstartpar,jendpar
                  DO i=istartpar,iendpar
                     a(i,j,ke1,n)=z0
                  END DO
                  END DO
               END DO

!              Fuer die darueberliegenden Schichten werden die lokalen
!              Gradienten der eigentlichen thermodyn. Modellvariablen
!              benoetigt:

               DO k=2,ke
                  DO j=jstartpar,jendpar
                  DO i=istartpar,iendpar
                     edh=z1/dicke(i,j,k)
                     a(i,j,k,tem)=(t(i,j,k-1)-t(i,j,k))*edh+tet_g
                     a(i,j,k,vap)=(qv(i,j,k-1)-qv(i,j,k))*edh
                     a(i,j,k,liq)=(qc(i,j,k-1)-qc(i,j,k))*edh
                  END DO
                  END DO
               END DO

               IF (imode_turb.EQ.3) THEN
!                 Die impliziten Anteile werden mit den thermodyn.
!                 Erhaltungsgroessen bestimmt.
!                 Fuer die explizite Korrektur werden die Differenzen
!                 zwischen den effektiven Gradienten und den lokalen
!                 Gradienten der zugehoerigen Erhaltungsgroessen
!                 benoetigt:

                  DO k=2,ke
                     DO j=jstartpar,jendpar
                     DO i=istartpar,iendpar
                        a(i,j,k,tem)=vari(i,j,k,tem) & !effectiver Temp.grad.
                         -(a(i,j,k,tem)-zlhocp*a(i,j,k,liq))*a(i,j,k,1) !eff. Grad. der Fluess.wass.temp.
                        a(i,j,k,vap)=vari(i,j,k,vap) &
                         -(a(i,j,k,vap)+liqfak*a(i,j,k,liq))
                        a(i,j,k,liq)=vari(i,j,k,liq)
                     END DO   
                     END DO   
                  END DO   
               ELSE   
!                 Die imliziten Anteile werden mit den eigentlichen
!                 thermodyn. Modellvariablen bestimmt.
!                 Fuer die explizite Korrektur werden die Differenzen
!                 zwischen den effektiven und den lokalen Gradienten
!                 dieser Variablen benoetigt:

                  DO n=nk2,ndiff
                     DO k=2,ke
                        DO j=jstartpar,jendpar
                        DO i=istartpar,iendpar
                           a(i,j,k,n)=vari(i,j,k,n)-a(i,j,k,n)
                        END DO
                        END DO
                     END DO
                  END DO
               END IF   
            END IF

!           Addition des Gradienten, welcher zur Temperatur-
!           flussdichte durch Drucktransport gehoert:

            IF (ltmpcor.AND.pat_len.GT.z0) THEN
               DO k=2,ke1 
                  DO j=jstartpar,jendpar
                  DO i=istartpar,iendpar
                     a(i,j,k,tem)=a(i,j,k,tem) &
                                  -tketens(i,j,k)/(cp_d*a(i,j,k,1))
                  END DO
                  END DO
               END DO
            END IF

!           Berechnung der zugehoerigen Flussdichten (mit aus numer.
!           Gruenden reduzierten Diffusionskoeff.):

            DO k=2,ke
               DO j=jstartpar,jendpar
               DO i=istartpar,iendpar
                  lays(i,j,1)=rhon(i,j,k) &
                              *MIN(wert*dicke(i,j,k)**2,tkvh(i,j,k))
               END DO
               END DO
               DO n=nk2,ndiff
                  DO j=jstartpar,jendpar
                  DO i=istartpar,iendpar
                     a(i,j,k,n)=-lays(i,j,1)*a(i,j,k,n)
                  END DO 
                  END DO 
               END DO 
            END DO 

!           Quasi impliziete Behandlung der untersten Schicht (Prandtl-
!           schicht):

            DO j=jstartpar,jendpar
            DO i=istartpar,iendpar
               lays(i,j,1)=rhon(i,j,ke1)*tkvh(i,j,ke1)
               lays(i,j,2)=lay(i,j) &
                   /(lays(i,j,1)*dt_var/(dicke(i,j,ke1)*securi)+lay(i,j))
            END DO 
            END DO 
            DO n=nk2,ndiff
               DO j=jstartpar,jendpar
               DO i=istartpar,iendpar

!                 Durch Anpassung der Flussdichte durch die ke-te
!                 Nebenflaeche, kann als untere Randbedingung der
!                 Fluss aus den Boden erhalten bleiben:

                  a(i,j,ke1,n)=-lays(i,j,1)*a(i,j,ke1,n)
                  a(i,j,ke,n)=a(i,j,ke1,n)*(z1-lays(i,j,2)) &
                             +a(i,j,ke,n)*lays(i,j,2)
               END DO
               END DO 
            END DO 

         END IF   

         IF (nk1.LE.nvel) THEN

!           Berechnung der korrigierten Gradienten der 
!           Impulskomponenten:

            IF (imode_turb.EQ.2) THEN
!              Die gesamte Rechnung ist explizit, daher werden die
!              unkorrigierten Gradienten benutzt:

               DO n=nk1,nvel
                  DO k=2,ke1
                     DO j=jstartpar,jendpar
                     DO i=istartpar,iendpar
                        a(i,j,k,n)=vari(i,j,k,n)
                     END DO
                     END DO
                  END DO
               END DO
            ELSE   
!              Es wird ein Teil der Flussdichtedivergenzen implizit
!              berechnet. Daher werden nur korrigierte Gradienten
!              benutzt, so dass das explizite Verfahren nur den bisher 
!              nicht berechneten Flussdivergenzanteil uebernimmt.
!              Dies sind bei den Impulskomponenten nur die Effekte
!              der nicht-lokalen Gradienten

!              Am unteren Rand ist keine explizite Korrektur noetig,
!              weil hier kein Unterschied zwischen lokaler und 
!              nicht-lokaler Gradientberechnung besteht:

               DO n=1,nvel
                  DO j=jstartpar,jendpar
                  DO i=istartpar,iendpar
                     a(i,j,ke1,n)=z0
                  END DO
                  END DO
               END DO

               DO n=nk1,nvel
                  DO k=2,ke
                     DO j=jstartpar,jendpar
                     DO i=istartpar,iendpar
                        a(i,j,k,n)=vari(i,j,k,n) &
                        -(wind(i,j,k-1,n)-wind(i,j,k,n))/dicke(i,j,k)
                     END DO
                     END DO
                  END DO
               END DO
            END IF   

!           Berechnung der zugehoerigen Flussdichten (mit aus numer.
!           Gruenden reduzierten Diffusionskoeff.):

            DO k=2,ke 
               DO j=jstartpar,jendpar
               DO i=istartpar,iendpar
                  lays(i,j,1)=rhon(i,j,k) &
                              *MIN(wert*dicke(i,j,k)**2,tkvm(i,j,k))
               END DO   
               END DO   
               DO n=nk1,nvel
                  DO j=jstartpar,jendpar
                  DO i=istartpar,iendpar
                     a(i,j,k,n)=-lays(i,j,1)*a(i,j,k,n)
                  END DO
                  END DO 
               END DO 
            END DO 

!           Quasi impliziete Behandlung der untersten Schicht (Prandtl-
!           schicht):

            DO j=jstartpar,jendpar
            DO i=istartpar,iendpar
               lays(i,j,1)=rhon(i,j,ke1)*tkvm(i,j,ke1) 
               lays(i,j,2)=lay(i,j) &
                   /(lays(i,j,1)*dt_var/(dicke(i,j,ke1)*securi)+lay(i,j))
            END DO   
            END DO   
            DO n=nk1,nvel 
               DO j=jstartpar,jendpar
               DO i=istartpar,iendpar

!                 Durch Anpassung der Flussdichte durch die ke-te
!                 Nebenflaeche, kann als untere Randbedingung der
!                 Fluss aus den Boden erhalten bleiben:

                  a(i,j,ke1,n)=-lays(i,j,1)*a(i,j,ke1,n)
                  a(i,j,ke,n)=a(i,j,ke1,n)*(z1-lays(i,j,2)) &
                             +a(i,j,ke,n)*lays(i,j,2)
               END DO
               END DO 
            END DO 

         END IF   

      END IF


!     Effektive Schichtdicke multipliziert mit der Dichte:

      IF (lexplizit.OR.limplizit) THEN
         kk=1
      ELSE
         kk=kcm
      END IF

      DO k=kk,ke
         DO j=jstartpar,jendpar
         DO i=istartpar,iendpar
            dicke(i,j,k)=(rhon(i,j,k)+rhon(i,j,k+1)) &
                        *(hhl(i,j,k)-hhl(i,j,k+1))*z1d2
         END DO
         END DO
      END DO


! 9)  Explizite Berechnung von Bestandtendenzen:

      IF (kcm.LE.ke) THEN

!        Berechnung des Formwiderstandes im Bestand:

         ku=1
         ko=2
         DO j=jstartpar,jendpar
         DO i=istartpar,iendpar
            lays(i,j,ku)=z0
         END DO
         END DO
         DO k=ke,kcm,-1   
            DO j=jstartpar,jendpar
            DO i=istartpar,iendpar
               lays(i,j,ko)=w(i,j,k)
               wert=z1d2*(c_big(i,j,k)  +c_sml(i,j,k) &
                         +c_big(i,j,k+1)+c_sml(i,j,k+1)) &
                   *SQRT(wind(i,j,k,u_m)**2+wind(i,j,k,v_m)**2 &
                        +(z1d2*(lays(i,j,ku)+lays(i,j,ko)))**2)

!              Durch die min-Funktion, wird verhindert, dass durch
!              die Wirkung der Formreibung eine Richtungsumkehr
!              der Geschwindigkeit erfolgen kann:

               wert=MIN(z1/dt_var,wert)

               wind(i,j,k,u_m)=-wert*wind(i,j,k,u_m)
               wind(i,j,k,v_m)=-wert*wind(i,j,k,v_m)
            END DO   
            END DO   
            kk=ku
            ku=ko
            ko=kk
         END DO   

!        Berechnung des Volumeneffektes im Bestand:

         DO n=1,ndiff
            ku=1
            ko=2
            IF (n.LE.nvel) THEN
               DO j=jstartpar,jendpar
               DO i=istartpar,iendpar
                  lays(i,j,ku)=-rhon(i,j,ke1)*tkvm(i,j,ke1) &
                                             *vari(i,j,ke1,n)
               END DO
               END DO
            ELSE
               DO j=jstartpar,jendpar
               DO i=istartpar,iendpar
                  lays(i,j,ku)=-rhon(i,j,ke1)*tkvh(i,j,ke1) &
                                             *vari(i,j,ke1,n)
               END DO
               END DO
            END IF

            DO k=ke,kcm,-1
               IF (kcm.EQ.1) THEN
                  DO j=jstartpar,jendpar
                  DO i=istartpar,iendpar
                     lays(i,j,ku)=lays(i,j,ku)*z1d2 &
                         *(r_air(i,j,k)-r_air(i,j,k+1))/dicke(i,j,k)
                  END DO
                  END DO
               ELSE   
                  DO j=jstartpar,jendpar
                  DO i=istartpar,iendpar
                     lays(i,j,ko)=-rhon(i,j,k)*tkvh(i,j,k)*vari(i,j,k,n)
                     lays(i,j,ku)=(lays(i,j,ko)+lays(i,j,ku))*z1d2 &
                         *(r_air(i,j,k)-r_air(i,j,k+1))/dicke(i,j,k)
                  END DO
                  END DO
               END IF
               IF (n.LE.nvel) THEN
                  DO j=jstartpar,jendpar
                  DO i=istartpar,iendpar
                     wind(i,j,k,n)=wind(i,j,k,n)-lays(i,j,ku)
                  END DO
                  END DO
               ELSEIF (n.EQ.tem) THEN
                  DO j=jstart,jend
                  DO i=istart,iend
                     t_tens(i,j,k)=t_tens(i,j,k)-lays(i,j,ku)
                  END DO
                  END DO
               ELSEIF (n.EQ.vap) THEN
                  DO j=jstart,jend
                  DO i=istartpar,iendpar
                     qv_tens(i,j,k)=qv_tens(i,j,k)-lays(i,j,ku)
                  END DO
                  END DO
               ELSEIF (n.EQ.liq) THEN
                  DO j=jstart,jend
                  DO i=istart,iend
                     qc_tens(i,j,k)=qc_tens(i,j,k)-lays(i,j,ku)
                  END DO
                  END DO
               END IF
               kk=ku
               ku=ko
               ko=kk 
            END DO    
         END DO    

      END IF   

! 10) Berechnung der Flussdichtedivergenzen:

! 10a) Vorbereitungen: 

      IF (lexplizit.OR.limplizit) THEN

!        Im Feld wind() stehen jetzt innerhalb der Bestandes die dort
!        anfallenden Tendenzen der horizontalen Windkomponenten.
!        I.f. werden weitere Wind-Tendenzen auf dieses Feld addiert.
!        Daher muss es ausserhalb des Bestandes auf 0 gesetzt werden:
 
         DO k=1,kcm-1
            DO j=jstartpar,jendpar
            DO i=istartpar,iendpar
               wind(i,j,k,u_m)=z0     
               wind(i,j,k,v_m)=z0     
            END DO
            END DO
         END DO

      END IF    

#ifdef SCLM
      IF (lscm) THEN

!        Nun werden die unter vari() abgespeicherten effektiven
!        Gradienten nur noch fuer die Berechnung der Komponenten
!        des turbulenten Spannungstensors benoetigt, welche u.U.
!        auch als Eingangsfeld fuer die Bestimmung der turbulenten
!        Horizontaldiffusion gebraucht werden:
         
         DO k=2,ke1
            DO j=jstartpar,jendpar
            DO i=istartpar,iendpar
               flux(u_s)=-tkvm(i,j,k)*vari(i,j,k,u_s)
               flux(v_s)=-tkvm(i,j,k)*vari(i,j,k,v_s)

               x0=z1d3*tke(i,j,k,ntur)**2
               fakt=d_m*len_scale(i,j,k)/tke(i,j,k,ntur)
               x1=flux(u_s)*vari(i,j,k,u_m)
               x2=flux(v_s)*vari(i,j,k,v_m)
               x3=-tkvh(i,j,k)*hlp(i,j,k)

               vari(i,j,k,1)=flux(u_s)
               vari(i,j,k,2)=flux(v_s)
               vari(i,j,k,3)=x0+fakt*(-z4*x1+z2*x2-z2*x3)
               vari(i,j,k,4)=x0+fakt*(+z2*x1-z4*x2-z2*x3)
            END DO
            END DO
         END DO

!        Jetzt werden die in hlp() gespeicherten fh2-Werte nicht mehr
!        benoetigt und hlp() kann anderweitig benutzt werden.
         
      END IF     
#endif


! 10b) Berechnung des expliziten Anteiles der Flussdichtedivergenzen:

       IF (lexplizit) THEN

!         Explizite Berechnung der Flussdichtedivergenzen:

            DO n=nkorr,ndiff
              
               DO j=jstartpar,jendpar
               DO i=istartpar,iendpar
                  hlp(i,j,1)=-a(i,j,2,n)
               END DO
               END DO
               DO k=2,ke
                  DO j=jstartpar,jendpar
                  DO i=istartpar,iendpar
                     hlp(i,j,k)=a(i,j,k,n)-a(i,j,k+1,n)
                  END DO
                  END DO
               END DO

!              Die folgende (quellenfreie) vertikale Glaettung ist
!              bei groesseren Zeitschritten (> ca. 30sec) in der
!              expliziten Variante noetig:

               DO j=jstartpar,jendpar
               DO i=istartpar,iendpar
                  a(i,j,1,n)=-(hlp(i,j,1)+wichfakt*(-hlp(i,j,1) &
                              +hlp(i,j,2)))/dicke(i,j,1)
               END DO
               END DO
               DO k=2,ke-1
                  DO j=jstartpar,jendpar
                  DO i=istartpar,iendpar
                     a(i,j,k,n)= &
                      -(hlp(i,j,k)+wichfakt*(hlp(i,j,k-1)-2*hlp(i,j,k) &
                                         +hlp(i,j,k+1)))/dicke(i,j,k)
                  END DO
                  END DO
               END DO
               DO j=jstartpar,jendpar
               DO i=istartpar,iendpar
                  a(i,j,ke,n)=-(hlp(i,j,ke)+wichfakt*(hlp(i,j,ke-1) &
                               -hlp(i,j,ke)))/dicke(i,j,ke)
               END DO
               END DO

!              Beruecksichtigung der vertikalen Cp-Variationen bei
!              der Berechnung der Temperaturtendenzen

               IF (n.eq.tem) THEN
                  DO k=1,ke
                     DO j=jstartpar,jendpar
                     DO i=istartpar,iendpar
                        a(i,j,k,tem)=a(i,j,k,tem) &
                         /(z1+zrcpv*qv(i,j,k)+zrcpl*qc(i,j,k))
                     END DO
                     END DO
                  END DO
!                 Beachte:
!                 Der Einfachheit halber wurde hier der Cp-Faktor
!                 ohne Berucksichtigung der subskaligen Wolken-
!                 wasserkorrektur berechnet!
               END IF   

            END DO   

!           Abspeichern der Tendenzen in den Tendenzfeldern, damit
!           das Feld a() fuer die ev. folgende implizite Rechnung
!           verfuegbar ist:

            DO n=nk1,nvel
               DO k=1,ke
                  DO j=jstartpar,jendpar
                  DO i=istartpar,iendpar
                     wind(i,j,k,n)=wind(i,j,k,n)+a(i,j,k,n)
                  END DO
                  END DO
               END DO
            END DO  

            DO k=1,ke
               DO j=jstart,jend
               DO i=istart,iend
                  t_tens(i,j,k)=t_tens(i,j,k)+a(i,j,k,tem)
                  qv_tens(i,j,k)=qv_tens(i,j,k)+a(i,j,k,vap)
                  qc_tens(i,j,k)=qc_tens(i,j,k)+a(i,j,k,liq)
               END DO
               END DO
            END DO

      END IF   


! 10c) Berechnung des impliziten Anteiles der Flussdichtedivergenzen:

      IF (limplizit) THEN

            DO k=1,ke
               DO j=jstartpar,jendpar
               DO i=istartpar,iendpar
                  dicke(i,j,k)=dicke(i,j,k)/dt_var
               END DO
               END DO
            END DO

            DO j=jstartpar,jendpar
            DO i=istartpar,iendpar
               a(i,j,1,1)=z0
            END DO
            END DO

            DO n=1,nvel

!              Bestimmung der Tendenzen der horiz. Windgeschw.:

               DO k=2,ke 
                  DO j=jstartpar,jendpar
                  DO i=istartpar,iendpar
                     a(i,j,k,1)=tkvm(i,j,k) &
                            *z2*rhon(i,j,k)/(hhl(i,j,k-1)-hhl(i,j,k+1))
                  END DO
                  END DO
               END DO
               DO j=jstartpar,jendpar
               DO i=istartpar,iendpar
                  a(i,j,ke1,1)=tkvm(i,j,ke1) &
                              *rhon(i,j,ke1)/dicke(i,j,ke1)
               END DO
               END DO

               IF (n.EQ.u_m) THEN
                  DO k=1,ke
                     DO j=jstartpar,jendpar
                     DO i=istartpar,iendpar
                        ii=MAX(i-1,1)
                        hlp(i,j,k)=(u(i,j,k)+u(ii,j,k))*z1d2
                     END DO
                     END DO
                  END DO   
               ELSEIF (n.EQ.v_m) THEN
                  DO k=1,ke
                     DO j=jstartpar,jendpar
                     DO i=istartpar,iendpar
                        jj=MAX(j-1,1)
                        hlp(i,j,k)=(v(i,j,k)+v(i,jj,k))*z1d2
                     END DO
                     END DO
                  END DO   
               END IF   
               DO j=jstartpar,jendpar
               DO i=istartpar,iendpar
                  hlp(i,j,ke1)=hlp(i,j,ke)*(z1-tfm(i,j))
               END DO
               END DO

               DO k=1,ke-1
                  DO j=jstartpar,jendpar
                  DO i=istartpar,iendpar
                     a(i,j,k,3)=a(i,j,k+1,1)
                     a(i,j,k,2)=-(a(i,j,k,1)+a(i,j,k,3)+dicke(i,j,k))
                     a(i,j,k,4)=-dicke(i,j,k)*hlp(i,j,k)
                  END DO
                  END DO
               END DO
               DO j=jstartpar,jendpar
               DO i=istartpar,iendpar
                  a(i,j,ke,2)=-(a(i,j,ke,1)+a(i,j,ke1,1) &
                                           +dicke(i,j,ke))
                  a(i,j,ke,4)=-dicke(i,j,ke)*hlp(i,j,ke)
               END DO
               END DO

               DO k=2,ke
                  DO j=jstartpar,jendpar
                  DO i=istartpar,iendpar
                     a(i,j,k,1)=a(i,j,k,1)/a(i,j,k-1,2)
                     a(i,j,k,2)=a(i,j,k,2)-a(i,j,k,1)*a(i,j,k-1,3)
                     a(i,j,k,4)=a(i,j,k,4)-a(i,j,k,1)*a(i,j,k-1,4)
                  END DO
                  END DO
               END DO

               DO j=jstartpar,jendpar
               DO i=istartpar,iendpar
                  a(i,j,ke,4)=a(i,j,ke,4)/a(i,j,ke,2)
               END DO
               END DO
               DO k=ke-1,1,-1
                  DO j=jstartpar,jendpar
                  DO i=istartpar,iendpar
                     a(i,j,k,4)=(a(i,j,k,4)-a(i,j,k,3)*a(i,j,k+1,4)) &
                                                       /a(i,j,k,2) 
                  END DO
                  END DO
               END DO

               DO k=1,ke
                  DO j=jstartpar,jendpar
                  DO i=istartpar,iendpar
                     wind(i,j,k,n)=wind(i,j,k,n) &
                                  +(a(i,j,k,4)-hlp(i,j,k))/dt_var
                  END DO
                  END DO
               END DO

            END DO   

            DO n=nred,nvel+1,-1

!              Bestimmung der Tendenzen der reduzierten scalaren Var..
!              Dabei soll tet_l als letztes behandelt werden:

               DO k=2,ke1
                  DO j=jstartpar,jendpar
                  DO i=istartpar,iendpar
                     a(i,j,k,1)=tkvh(i,j,k) &
                            *z2*rhon(i,j,k)/(hhl(i,j,k-1)-hhl(i,j,k+1))
                  END DO
                  END DO
               END DO
               DO j=jstartpar,jendpar
               DO i=istartpar,iendpar
                  a(i,j,ke1,1)=tkvh(i,j,ke1) &
                              *rhon(i,j,ke1)/dicke(i,j,ke1)
               END DO
               END DO

               IF (n.EQ.tet_l) THEN
                  ko=1
                  ku=2
                  DO j=jstartpar,jendpar
                  DO i=istartpar,iendpar
                     fakt=zexner(prs(i,j,1))
                     wert=z1+zrcpv*qv(i,j,1)+zrcpl*qc(i,j,1)

                     hlp(i,j,1)=(t(i,j,1) &
                                -zlhocp*qc(i,j,1))/fakt
                     dicke(i,j,1)=dicke(i,j,1)*wert
                     lays(i,j,ko)=fakt*wert
                  END DO
                  END DO

                  DO k=2,ke   
                     DO j=jstartpar,jendpar
                     DO i=istartpar,iendpar
                        fakt=zexner(prs(i,j,k))
                        wert=z1 &
                            +zrcpv*qv(i,j,k)+zrcpl*qc(i,j,k)

                        hlp(i,j,k)=(t(i,j,k) &
                                   -zlhocp*qc(i,j,k))/fakt
                        dicke(i,j,k)=dicke(i,j,k)*wert
                        lays(i,j,ku)=fakt*wert

                        wert=(lays(i,j,ku)*dp0(i,j,k-1) &
                             +lays(i,j,ko)*dp0(i,j,k)) &
                            /(dp0(i,j,k-1)+dp0(i,j,k))
                        a(i,j,k,1)=a(i,j,k,1)*wert
                     END DO   
                     END DO   
                     kk=ko
                     ko=ku
                     ku=kk
                  END DO   
                  DO j=jstartpar,jendpar
                  DO i=istartpar,iendpar
!                    Beachte:
!                    tp(), qd(), ql() und pr() sind fuer die 
!                    Schicht k=ke1 bereits vorhanden.
     
                     fakt=zexner(pr(i,j))
                     wert=z1+zrcpv*qd(i,j)+zrcpl*ql(i,j)

                     hlp(i,j,ke1)=(tp(i,j)-zlhocp*ql(i,j))/fakt

                     a(i,j,ke1,1)=a(i,j,ke1,1)*fakt*wert
                  END DO
                  END DO
!                 Weil tet_l als letztes behandelt wird, konnte
!                 dicke() mit der speziellen Korrektur bzgl. 
!                 cp_fakt ueberschrieben werden!
               ELSEIF (n.EQ.h2o_g) THEN
                  DO k=1,ke
                     DO j=jstartpar,jendpar
                     DO i=istartpar,iendpar
                        hlp(i,j,k)=qv(i,j,k)+liqfak*qc(i,j,k)
                     END DO
                     END DO
                  END DO   
                  DO j=jstartpar,jendpar
                  DO i=istartpar,iendpar
!                    Beachte:
!                    qd() und ql() sind fuer die Schicht k=ke1
!                    bereits vorhanden.

                     hlp(i,j,ke1)=qd(i,j)+liqfak*ql(i,j)
                  END DO
                  END DO
               END IF    
               DO k=1,ke-1
                  DO j=jstartpar,jendpar
                  DO i=istartpar,iendpar
                     a(i,j,k,3)=a(i,j,k+1,1)
                     a(i,j,k,2)=-(a(i,j,k,1)+a(i,j,k,3)+dicke(i,j,k))
                     a(i,j,k,4)=-dicke(i,j,k)*hlp(i,j,k)
                  END DO
                  END DO
               END DO

               DO j=jstartpar,jendpar
               DO i=istartpar,iendpar
                  a(i,j,ke,2)=-(a(i,j,ke,1)+dicke(i,j,ke))
                  a(i,j,ke,4)=-dicke(i,j,ke)*hlp(i,j,ke) &
                             +a(i,j,ke1,1)*(hlp(i,j,ke)-hlp(i,j,ke1))
               END DO
               END DO

               DO k=2,ke
                  DO j=jstartpar,jendpar
                  DO i=istartpar,iendpar
                     a(i,j,k,1)=a(i,j,k,1)/a(i,j,k-1,2)
                     a(i,j,k,2)=a(i,j,k,2)-a(i,j,k,1)*a(i,j,k-1,3)
                     a(i,j,k,4)=a(i,j,k,4)-a(i,j,k,1)*a(i,j,k-1,4)
                  END DO
                  END DO
               END DO

               DO j=jstartpar,jendpar
               DO i=istartpar,iendpar
                  a(i,j,ke,4)=a(i,j,ke,4)/a(i,j,ke,2)
               END DO
               END DO
               DO k=ke-1,1,-1
                  DO j=jstartpar,jendpar
                  DO i=istartpar,iendpar
                     a(i,j,k,4)=(a(i,j,k,4)-a(i,j,k,3)*a(i,j,k+1,4)) &
                                                       /a(i,j,k,2) 
                  END DO
                  END DO
               END DO

               IF (n.EQ.tet_l) THEN
                  DO k=1,ke
                     DO j=jstart,jend
                     DO i=istart,iend
                        t_tens(i,j,k)=t_tens(i,j,k) &
                                     +(a(i,j,k,4)-hlp(i,j,k))/dt_var
                     END DO
                     END DO
                  END DO
               ELSEIF (n.EQ.h2o_g) THEN
                  DO k=1,ke
                     DO j=jstart,jend
                     DO i=istart,iend
                        qv_tens(i,j,k)=qv_tens(i,j,k) &
                                     +(a(i,j,k,4)-hlp(i,j,k))/dt_var
                     END DO
                     END DO
                  END DO
               END IF   
                
            END DO   

      END IF   


! 11) Berechnung der Tendenzen infloge horizontaler turb. Diffusion
!     (und aufsummieren auf die Tendenzfelder):

!        ** wird erst spaeter eingefuehrt **


! 12) Aktualisierung der Tendenzfelder fuer den horizontalen Wind:

      IF (limplizit .OR. nk1.EQ.1 .OR. kcm.LE.ke) THEN

#ifdef __COSMO__
!        Interpolation der Diffusionstendenzen von horizontalen
!        Impulskomponenten auf die 'gestaggerten' Positionen:

         kzdims=(/ke,ke,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
         CALL exchg_boundaries                                             &
              (  0  ,  sendbuf, isendbuflen, imp_reals, icomm_cart,        &
              num_compute, ie, je, kzdims, jstartpar, jendpar, 2,          &
              nboundlines, my_cart_neigh, lperi_x, lperi_y, l2dim,         &
              2200+nexch_tag, .FALSE.   , ncomm_type, izerror, yzerrmsg,   &
              wind(:,:,:,u_m),wind(:,:,:,v_m))
#endif

         IF (limplizit.OR.(nk1.EQ.1)) THEN
            kk=1
         ELSE
            kk=kcm
         END IF

         DO k=kk,ke   
            DO j=jstartu,jendu     
            DO i=istartu,iendu   
               u_tens(i,j,k)=u_tens(i,j,k) &
                           +(wind(i,j,k,u_m)+wind(i+1,j,k,u_m))*z1d2
            END DO
            END DO
            DO j=jstartv,jendv     
            DO i=istartv,iendv   
               v_tens(i,j,k)=v_tens(i,j,k) &
                           +(wind(i,j,k,v_m)+wind(i,j+1,k,v_m))*z1d2
            END DO
            END DO
         END DO

      END IF
                  
! 13) Berechnung der Diffusionstendenz von SQRT(2*TKE):

!        Interpolationen auf Hauptflaechen fuer die Standardabweichnung
!        des Saettigungsdefizites und den Drucktransport:

         IF (pat_len.GT.z0) THEN
            DO k=2,ke1
               DO j=jstartpar,jendpar
               DO i=istartpar,iendpar  
                  tketens(i,j,k)=rhon(i,j,k)*tkvh(i,j,k)/tke(i,j,k,ntur) &
                                            *tketens(i,j,k)*len_scale(i,j,k)
               END DO
               END DO
            END DO
            DO k=2,ke
               DO j=jstartpar,jendpar
               DO i=istartpar,iendpar  
                  tketens(i,j,k)=(tketens(i,j,k)+tketens(i,j,k+1)) &
                                /(len_scale(i,j,k)+len_scale(i,j,k+1))
               END DO
               END DO
            END DO
         END IF

         DO j=jstartpar,jendpar
         DO i=istartpar,iendpar  
            rcld(i,j,1)=rcld(i,j,2)
         END DO
         END DO
         DO k=2,kem-1
            DO j=jstartpar,jendpar
            DO i=istartpar,iendpar  
               rcld(i,j,k)=(rcld(i,j,k)+rcld(i,j,k+1))*z1d2
            END DO
            END DO
         END DO
!        Fuer die unterste Hauptflaeche (k=ke) wird bei kem=ke
!        der Wert auf der entspr. Nebenflaeche beibehalten.

!        tketens wird hiernach als Teilbetrag der SQRT(2*TKE)-Flussdichte
!        behandelt. 

!        Berechnung der vertikalen SQRT(2*TKE)-Flussdichten
!        und Ablegen der selben im Feld tketens:
          
!        Weil die TKE-Gleichung eine Gleichung 2-ter Ordnung ist und
!        die Parametrisierung der Diffusion in dieser Gleichung ohne
!        Zuhilfenahme von Gleichungen hoeherer Ordnung quasi ad hoc 
!        erfolgt, sind die Genauhigkeitsansprueche der numerischen
!        Berechnung ensprechend geringer als bei der Berechnung der
!        Diffusion in den Gleichungen 1-ter Ordnung. Es wird daher
!        hier auf nicht-lokale Gradientbildung und implizite Behandlung
!        generell verzichtet.
!
!        If limpltkediff=.false. we compute the original explicit diffusion
!        of TKE. Note that in order to achieve stability with this
!        explicit scheme, the following condition 0 <= securi < 0.5
!        should be met. Also, for very shallow model levels and TKE
!        values on the order of 1 m2/s2, the diffusion of TKE is
!        strongly limited by stability constraints, and an implicit
!        solution (by choosing limpltkediff=.true.) should be favored.
!
!        Nevertheless, we need for the present at least an exlicit correction
!        in order to express the circulation term and roughness layer terms.

!        Es soll keinen turb. TKE-FLuss aus der Atm. hinaus geben:

         DO j=jstartpar,jendpar
         DO i=istartpar,iendpar
            tketens(i,j,1)=z0
         END DO
         END DO

         IF (.not.limpltkediff) THEN

!           Bei der reinen explizite Diffusion von TKE muss die explizite
!           Flussdichte berechnet werden:

            DO k=ke,2,-1
               DO j=jstartpar,jendpar
               DO i=istartpar,iendpar
                  dicke(i,j,k)=hhl(i,j,k)-hhl(i,j,k+1)
   
                  q3=(tke(i,j,k,ntur)+tke(i,j,k+1,ntur))*z1d2
                  kh=q3*(len_scale(i,j,k)+len_scale(i,j,k+1))*z1d2
                  kh=(rhon(i,j,k)+rhon(i,j,k+1))*z1d2 &
                    *MIN(securi*dicke(i,j,k)**2/dt_tke,c_diff*kh)
                     
!-----------------------------------------------------------------------
!                    nur im 1-d-Modell:
   
!                    if (ldruck) then
!                       write (7,*) 'k,tke,p_flux,tke_flux' &
!                      ,k,tke(i,j,k,ntur),tketens(i,j,k)    &
!                      ,kh*q3*(tke(i,j,k+1,ntur)-tke(i,j,k,ntur)) &
!                            /dicke(i,j,k)
!                       end if
!-----------------------------------------------------------------------

                  tketens(i,j,k)=tketens(i,j,k) &
                                +kh*(tke(i,j,k+1,ntur)-tke(i,j,k,ntur))&
                                   /dicke(i,j,k)
               END DO
               END DO
            END DO

!           Bem.: Der viskosen Zusatzterm kann auch entfallen!

!           Die Diffusionskoeffizienten wurden mit Hilfe der MIN-Funktion
!           auf das numerisch ertraegliche Mass reduziert.

!           Zuvor wurde tketens bereits mit dem Drucktransport belegt,
!           so dass tketens also die Summe aus der eigentlichen TKE-
!           Flussdichte und dem Drucktrnasport enthaelt!

!           Positiv definite Diffusionskorrektur in der Gleichung fuer SQRT(2*TKE)
!           wobei die Laengenskala c_diff*len_scale beschraenkt werden muss:

            DO k=2,ke
               DO j=jstartpar,jendpar
               DO i=istartpar,iendpar
                  len=hhl(i,j,k-1)-hhl(i,j,k+1)
                  hlp(i,j,k)=(tketens(i,j,k-1)-tketens(i,j,k)) &
                     +MIN(securi*len**2/(4*tke(i,j,k,ntur)*dt_tke), c_diff*len_scale(i,j,k)) &
                      *( (tke(i,j,k-1,ntur)-tke(i,j,k+1,ntur))/len )**2
 
               END DO
               END DO
            END DO

         ELSE

!           Bei der impliziten Diffusion von TKE muss nur die Flussdichte bzueglich des
!           Zirkulationstermes explizit behandelt werden:

            DO k=2,ke
               DO j=jstartpar,jendpar
               DO i=istartpar,iendpar
                  hlp(i,j,k)=tketens(i,j,k-1)-tketens(i,j,k) 
               END DO
               END DO
            END DO

         END IF    

!        Volumenkorrektur in der Rauhigkeitsschicht:

         DO k=MAX(2,kcm),ke
            DO j=jstartpar,jendpar
            DO i=istartpar,iendpar
               hlp(i,j,k)=hlp(i,j,k) &
                  +z1d2*(tketens(i,j,k)*dp0(i,j,k-1)+tketens(i,j,k-1)*dp0(i,j,k)) &
                       /(dp0(i,j,k)+dp0(i,j,k-1))*(r_air(i,j,k-1)-r_air(i,j,k+1))
            END DO
            END DO
         END DO

         DO j=jstartpar,jendpar
         DO i=istartpar,iendpar
            hlp(i,j,1)=hlp(i,j,2)
            hlp(i,j,ke1)=hlp(i,j,ke)
         END DO
         END DO

!        Die folgende (quellenfreie) vertikale Glaettung (mittels 'wichfakt')
!        ist bei groesseren Zeitschritten (> ca. 30sec) in der expliziten
!        Variante noetig:

         DO k=2,ke
            DO j=jstartpar,jendpar
            DO i=istartpar,iendpar
               dicke(i,j,k)=z1d2*(hhl(i,j,k-1)-hhl(i,j,k+1))*rhon(i,j,k)
               tketens(i,j,k)= &
                     -(hlp(i,j,k)+wichfakt*(hlp(i,j,k-1)-2*hlp(i,j,k) &
                                        +hlp(i,j,k+1)))/dicke(i,j,k)
            END DO
            END DO
         END DO

         IF (limpltkediff) THEN
!        *************** semi-implicit TKE diffusion tendency *************

         ! define implicit/explicit weight
         !          explicit : beta = 0.0
         !   Crank-Nicholson : beta = 0.5
         !          implicit : beta = 1.0
         beta=0.85_wp

         ! compute level thickness at full level height
         DO k=1,ke
            DO j=jstartpar,jendpar
            DO i=istartpar,iendpar
               hlp(i,j,k)=hhl(i,j,k)-hhl(i,j,k+1)
            END DO
            END DO
         END DO

         ! compute level thickness at half level height
         DO j=jstartpar,jendpar
         DO i=istartpar,iendpar
            dicke(i,j,1)=hlp(i,j,1)*rhon(i,j,1)
         END DO
         END DO
         DO k=2,ke
            DO j=jstartpar,jendpar
            DO i=istartpar,iendpar
               dicke(i,j,k)=z1d2*(hhl(i,j,k-1)-hhl(i,j,k+1))*rhon(i,j,k)
            END DO
            END DO
         END DO
         DO j=jstartpar,jendpar
         DO i=istartpar,iendpar
            dicke(i,j,ke1)=hlp(i,j,ke)*rhon(i,j,ke1)
         END DO
         END DO

         ! compute diffusion constant at full level height
         ! note: kh = c_diff * rhon * len_scale * q**2
         DO k=1,ke
            DO j=jstartpar,jendpar
            DO i=istartpar,iendpar
               a(i,j,k,1)=c_diff * z1d2 * (rhon(i,j,k)+rhon(i,j,k+1)) &
                                 * z1d2 * (len_scale(i,j,k)+len_scale(i,j,k+1)) &
                                 *(z1d2 * (tke(i,j,k,ntur)+tke(i,j,k+1,ntur)) )**2
               ! limiter to mimic old scheme
!               a(i,j,k,1)=min( a(i,j,k,1), &
!                    securi*hlp(i,j,k)**2/dt_tke/(z1d2*(rhon(i,j,k)+rhon(i,j,k+1))) )
            END DO
            END DO
         END DO

         ! Top layer (k=1)
         DO j=jstartpar,jendpar
         DO i=istartpar,iendpar
            ! setup matrix (a,b,c) and vector (d)
            tmp1       = dt_tke / tke(i,j,1,ntur) / dicke(i,j,1)
            tmp3       = tmp1 * a(i,j,1,1) / hlp(i,j,1)
            aa         = z0
            bb         = z0 !- tmp3
            a(i,j,1,3) = z0 !+ tmp3
            a(i,j,1,4) = + (z1+(z1-beta)*bb        )*tke(i,j,1,ntur) &
                         + (   (z1-beta)*a(i,j,1,3))*tke(i,j,2,ntur)
            bb         = z1 - beta*bb
            a(i,j,1,3) =    - beta*a(i,j,1,3)
            ! forward elimination
            a(i,j,1,3) = a(i,j,1,3)/bb
            a(i,j,1,4) = a(i,j,1,4)/bb
            a(i,j,1,2) = bb
         END DO
         END DO

         ! The layers from k=2 to k=ke1-1
         DO k=2,ke1-1
            DO j=jstartpar,jendpar
            DO i=istartpar,iendpar
               ! setup matrix (a,b,c) and vector (d)
               tmp1       = dt_tke / tke(i,j,k,ntur) / dicke(i,j,k)
               tmp2       = tmp1 * a(i,j,k-1,1) / hlp(i,j,k-1)
               tmp3       = tmp1 * a(i,j,k  ,1) / hlp(i,j,k  )
               aa         = + tmp2
               bb         = - tmp2 - tmp3
               a(i,j,k,3) = + tmp3
               a(i,j,k,4) = + (   (z1-beta)*aa        )*tke(i,j,k-1,ntur) &
                            + (z1+(z1-beta)*bb        )*tke(i,j,k  ,ntur) &
                            + (   (z1-beta)*a(i,j,k,3))*tke(i,j,k+1,ntur)
               aa         =    - beta*aa
               bb         = z1 - beta*bb
               a(i,j,k,3) =    - beta*a(i,j,k,3)
               ! forward elimination
               tmp1       = aa/a(i,j,k-1,2)
               a(i,j,k,4) = a(i,j,k,4) - tmp1*a(i,j,k-1,4)*a(i,j,k-1,2)
               a(i,j,k,2) = bb         - tmp1*a(i,j,k-1,3)*a(i,j,k-1,2)
               a(i,j,k,3) = a(i,j,k,3)/a(i,j,k,2)
               a(i,j,k,4) = a(i,j,k,4)/a(i,j,k,2)
            END DO
            END DO
         END DO

         ! The bottom layer (k=ke1)
         DO j=jstartpar,jendpar
         DO i=istartpar,iendpar
            ! setup matrix (a,b,c) and vector (d)
            tmp1         = dt_tke / tke(i,j,ke1,ntur) / dicke(i,j,ke1)
            tmp2         = tmp1 * a(i,j,ke1-1,1) / hlp(i,j,ke1-1)
            aa           = z0 !+ tmp2
            bb           = z0 !- tmp2
            a(i,j,ke1,3) = z0
            a(i,j,ke1,4) = + (   (z1-beta)*aa       )*tke(i,j,ke1-1,ntur) &
                           + (z1+(z1-beta)*bb       )*tke(i,j,ke1  ,ntur)
            aa           =    - beta*aa
            bb           = z1 - beta*bb
            ! forward elimination
            tmp1         = aa/a(i,j,ke1-1,2)
            a(i,j,ke1,4) = a(i,j,ke1,4) - tmp1*a(i,j,ke1-1,4)*a(i,j,ke1-1,2)
            a(i,j,ke1,2) = bb           - tmp1*a(i,j,ke1-1,3)*a(i,j,ke1-1,2)
            a(i,j,ke1,4) = a(i,j,ke1,4)/a(i,j,ke1,2)
         END DO
         END DO

         ! Backsubstitution and storage of new SQRT(2*TKE) field
         DO j=jstartpar,jendpar
         DO i=istartpar,iendpar
            hlp(i,j,ke1) = a(i,j,ke1,4)
         END DO
         END DO
         DO k=ke1-1,1,-1
            DO j=jstartpar,jendpar
            DO i=istartpar,iendpar
               hlp(i,j,k) = a(i,j,k,4) - a(i,j,k,3) * hlp(i,j,k+1)
            END DO
            END DO
         END DO

         ! Computation of SQRT(2*TKE) tendencies
         DO k=2,ke
            DO j=jstartpar,jendpar
            DO i=istartpar,iendpar
               tketens(i,j,k)=tketens(i,j,k)+ &
                              (hlp(i,j,k)-tke(i,j,k,ntur))/dt_tke
            END DO
            END DO
         END DO

!        *************** end of TKE diffusion tendency *************
         END IF

!        Im skin-layer soll es keine TKE-Tendenzen infolge
!        Diffusion oder Drucktransport geben:

         DO j=jstartpar,jendpar
         DO i=istartpar,iendpar
            tketens(i,j,ke1)=z0
         END DO
         END DO

!        Die folgende MAX-Funktion soll sicherstellen, dass allein durch die Diffusionstendenz
!        von SQRT(2*TKE) (einschl. des Zirkulationstermes) die TKE nicht negativ werden kann.

         DO k=2,ke
            DO j=jstartpar,jendpar
            DO i=istartpar,iendpar
               tketens(i,j,k)=MAX(-tke(i,j,k,ntur)/dt_tke, tketens(i,j,k))
            END DO
            END DO
         END DO

!        Jetzt enthaelt tketens die SQRT(2*TKE)-Tend. infolge Diffusion
!        und Drucktransport. Zu tketens wird nun an anderer Stelle noch
!        die Advektionstendenz hinzuaddiert und das ganze dann im
!        naechsten Zeitschritt fuer die SQRT(2*TKE)-Prognose benutzt.

! 14) Berechnung der turb. horizontalen TKE-Flussdichtedivergenzen
!     und Aufsummieren auf das SQRT(2*TKE)-Tendenzfeld:

!        ** wird erst spaeter eingefuehrt **

! 15) Deallocierung lokaler Felder

END SUBROUTINE turbdiff

!==============================================================================

END MODULE turbulence_diff
