! $RCSfile: src_slow_tendencies_rk.f90,v $
! $Revision: 4.21 $ $Date: 2012/03/30 $
!+ Source module for computing slow tendencies of Runge-Kutta dynamics
!------------------------------------------------------------------------------

MODULE src_slow_tendencies_rk

!------------------------------------------------------------------------------
!
! Description:
!   This module contains subroutines which compute the slow tendencies for
!   the Runge-Kutta version (irunge_kutta=1/2). 
!
!   These routines have been in module src_runge_kutta.f90 before, which 
!   has been splitted.
!
!   Routines currenctly included are:
!    - complete_tendencies_init
!    - complete_tendencies_uvwtpp
!    - complete_tendencies_uvwtpp_eva
!    - complete_tendencies_qvqcqi_tke
!    - explicit_horizontal_diffusion
!    - implicit_vert_diffusion_uvwt
!
! Current Code Owner: DWD, Jochen Foerstner    
!  phone:  +49  69  678667 35
!  fax:    +49  69  8062 3721
!  email:  jochen.foerstner@dwd.de
!
! History:
! Version    Date       Name
! ---------- ---------- ----
! 3.21       2006/12/04 Ulrich Schaettler
!  Initial release
! V3_23        2007/03/30 Jochen Foerstner
!  Corrections for l3dturb_metr
! V4_5         2008/09/10 Guenther Zaengl
!  Possibility to switch off the surface momentum fluxes in idealized cases
!  by setting lfreeslip_sfc=.TRUE.
! V4_9         2009/07/16 Ulrich Schaettler, Christian Bollmann (NEC)
!  Implemented NEC ON_ADB Directives; changed some j- and k-loops
! V4_10        2009/09/11 Ulrich Schaettler, Christian Bollmann (NEC)
!  Changed syntax of on_adb option for NEC because of new compiler
! V4_11        2009/11/30 Guenther Zaengl
!  Implemented use of cloud ice tendencies from convection scheme
! V4_12        2010/05/11 Michael Baldauf
!  Replace local variables zMlambda_s, zMphi_s by new metric coefficients
!  dzeta_dlam, dzeta_dphi, resp.
! V4_13        2010/05/11 Michael Gertz
!  Adaptions to SVN
! V4_15        2010/11/19 Michael Baldauf
!  new SUBROUTINE complete_tend_uvwtpp_CN3Crow
! V4_17        2011/02/24 Ulrich Blahak
!  Eliminated lfreeslip_surf (free-slip BC and/or no-surface-heat/moisture-flux
!    conditions can now be imposed by switches lnosurffluxes_m/h in namelist ARTIFCTL);
!  fixed use of lkge2_prog_tke in horizontal turbulent diffusion of TKE;
!  fixed various bugs in loop index boundaries in horizontal turbulent diffusion;
!  correction for ztch: introduced lower velocity limit vel_min (COMMENTED OUT FOR NOW)
! V4_18        2011/05/26 Ulrich Schaettler
!  Implementation of new SR complete_tendencies_cosmo_art, complete_tendencies_pollen
!  for COSMO-ART and Pollen (Christoph Knote)
!  Added comments for ztvb, ztch regarding the use of t_g instead of t_s and
!   the clipping of zvbke by vel_min. These things have to be tested further
!   and made consistent in the rest of the model in the future. (by Uli Blahak)
! V4_20        2011/08/31 Ulrich Schaettler
!  Activated lower velocity limit vel_min
! V4_21        2012/04/03 Alexander Kelbch / Markus Uebel
!  Inclusion of tracer tendencies (adopted from Markus Uebel)
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

USE data_parameters , ONLY :   &
  ireals,    & ! KIND-type parameters for real variables
  iintegers    ! KIND-type parameter for standard integer variables

!------------------------------------------------------------------------------

USE data_modelconfig, ONLY :   &
  
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
  
  !   meridional direction
  jstart,       & ! start index for the forecast of w, t, qv, qc, tracer (AK) and pp
  jend,         & ! end index for the forecast of w, t, qv, qc, tracer (AK) and pp
  jstartu,      & ! start index for the forecast of u
  jendu,        & ! start index for the forecast of u
  jstartv,      & ! start index for the forecast of v
  jendv,        & ! end index for the forecast of v
  
  ! 4. constants for the horizontal rotated grid and related variables
  ! ------------------------------------------------------------------

  eddlon,       & ! 1 / dlon
  eddlat,       & ! 1 / dlat

  ! 5. variables for the time discretization and related variables
  ! --------------------------------------------------------------

  dt              ! long time-step

! end of data_modelconfig

!------------------------------------------------------------------------------

USE data_constants  , ONLY :   &
  
  ! 2. physical constants and related variables
  ! -------------------------------------------

  r_d,          & ! gas constant for dry air
  rvd_m_o,      & ! r_v/r_d - 1
  cp_d,         & ! specific heat of dry air at constant pressure
  rdocp,        & ! r_d / cp_d
  lh_v,         & ! latent heat of vapourization
  g,            & ! acceleration due to gravity (for COSMO_ART)
  r_earth         ! mean radius of the earth

! end of data_constants

!------------------------------------------------------------------------------

USE data_fields     , ONLY :   &
  
  ! 1. constant fields for the reference atmosphere                 (unit)
  ! -----------------------------------------------
  p0         ,    & ! reference pressure at main levels             ( Pa  )
  hhl        ,    & ! geometical height of half model levels        (  m  )
  sqrtg_r_s  ,    & ! reciprocal square root of G at skalar points  ( 1/m )
  sqrtg_r_u  ,    & ! reciprocal square root of G at u points       ( 1/m )
  sqrtg_r_v  ,    & ! reciprocal square root of G at v points       ( 1/m )
  sqrtg_r_w  ,    & ! reciprocal square root of G at w points       ( 1/m )
  dzeta_dlam ,    & ! metric coefficient                            ( 1   )
  dzeta_dphi ,    & ! metric coefficient                            ( 1   )

  ! 2. external parameter fields                                    (unit)
  ! ----------------------------
  crlat      ,    & ! cosine of transformed latitude
  acrlat     ,    & ! 1 / ( crlat * radius of the earth )           ( 1/m )

  ! 3. prognostic variables                                         (unit)
  ! -----------------------
  u          ,    & ! zonal wind speed                              ( m/s )
  v          ,    & ! meridional wind speed                         ( m/s )
  w          ,    & ! vertical wind speed (defined on half levels)  ( m/s )
  t          ,    & ! temperature                                   (  k  )
  qv         ,    & ! specific water vapor content                  (kg/kg)
  qc         ,    & ! specific cloud water content                  (kg/kg)
  qi         ,    & ! specific cloud ice   content                  (kg/kg)
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
  tketens,        & ! tke-tendency                                  (m2/s3)
  pptens            ! pp-tendency without sound-wave terms          (Pa/s )

USE data_fields     , ONLY :   &

  ! 5. fields for surface values and soil model variables           (unit )
  ! -----------------------------------------------------
  ps        ,     & ! surface pressure                              ( pa  )
  t_s       ,     & ! temperature of the ground surface             (  k  )
  t_g       ,     & ! weighted surface temperature                  (  k  )
  qv_s      ,     & ! specific water vapor content on the surface   (kg/kg)
  
  ! 6. fields that are computed in the parametrization and dynamics (unit )
  ! ---------------------------------------------------------------
  rho,          & ! density of moist air
  qvt_diff  ,   & ! humidity    tendency  due to diffusion          ( 1/s )
  a1t, a2t  ,   & ! implicit weight of vertical diffusion
  
  !   turbulent coefficients in the atmosphere
  !   (defined on half levels)
  !   vertical turbulent diffusion coefficients
  tkvm     ,      & ! ... for momentum                               (m2/s)
  tkvh     ,      & ! ... for heat and moisture                      (m2/s)
  ! horizontal turbulent diffusion coefficients
  tkhm     ,      & ! ... for momentum                               (m2/s)
  tkhh     ,      & ! ... for heat and moisture                      (m2/s)

  !   turbulent coefficients at the surface
  tcm      ,    & ! turbulent diffusion coefficients for momentum   --
  tch      ,    & ! turbulent diffusion coefficients for heat       --
                  ! and moisture
#ifdef COUP_OAS_COS
  tcw      ,    & ! turbulent diffusion coefficients for moisture   --  !CPS
#endif
  !   fields that are computed in the dynamics
  dqvdt    ,    & ! threedimensional moisture convergence         (  1/s)
  umfl_s   ,    & ! average u-momentum flux (surface)             ( N/m2)
  vmfl_s   ,    & ! average v-momentum flux (surface)             ( N/m2)
  shfl_s   ,    & ! average sensible heat flux (surface)          ( W/m2)
  qvsflx   ,    & ! surface flux of water vapour                  (kg/m2s)
  lhfl_s   ,    & ! average latent heat flux (surface)            ( W/m2)
  wcon     ,    & ! contravariant vertical velocity
  uadvt    ,    & ! advective tendency of u
  vadvt    ,    & ! advective tendency of v
  wadvt    ,    & ! advective tendency of w
  ppadvt   ,    & ! advective tendency of pp
  tadvt           ! advective tendency of t

! end of data_fields
!------------------------------------------------------------------------------

!AK (30.03.2012)
USE data_tracer     , ONLY:   &
  tracer,       & ! tracer concentration                          ("var")
  tracertens,   & ! tracer tendency                               ("var"/s)
  ltracer,      & ! switch for tracer
  ntracer         ! number of tracers

USE data_runcontrol , ONLY: ntstep

USE data_parallel   , ONLY: my_cart_id
!AK (30.03.2012)

!------------------------------------------------------------------------------

USE data_runcontrol , ONLY :   &
  nnow,          & ! corresponds to ntstep
  nnew,          & ! corresponds to ntstep + 1
  lprog_qi,      & ! if .TRUE., running with cloud ice
  ltur,          & ! forecast with turbulent diffusion
  itype_turb,    & ! type of turbulent diffusion parametrization
  imode_turb,    & ! mode of turbulent diffusion parametrization
  l3dturb,       & ! 3D-turbulence (additional horizontal diffusion)
  l3dturb_metr,  & ! switch on/off additional metric terms for 3D-turbulence
  lprog_tke,     & ! prognostic treatment of TKE (for itype_turb=5/7)
  l_cosmo_art,   & ! if .TRUE., run the COSMO_ART
  l_pollen         ! of pollen

!------------------------------------------------------------------------------

! UB>>
USE data_turbulence , ONLY :   &
  vel_min          ! minimal velocity scale [m/s]
! UB<<

!------------------------------------------------------------------------------

USE numeric_utilities_rk,     ONLY :  &
  clipping                    !

!------------------------------------------------------------------------------

#ifdef POLLEN
USE data_pollen,         ONLY :     &
    isp_pollen   , & ! number of pollen species
    cpollen      , & ! pollen concentration                          (# m-3 )
    cpollen_s    , & ! surface pollen concentration                  (# m-3 )
    cpollentens  , & ! cpolle-tendency without sound-wave terms     (# m-3/s)
    cpollent_diff, & ! tendency by horizontal diffusion
    vsed_p           ! sedimentation velocity aerosols                  (m/s)
#endif

#ifdef COSMOART
USE data_cosmo_art,      ONLY :     &
    lgas         , & ! of gases
    laero        , & ! of aerosols
    isp_gastrans , & ! number of transported gas phase species
    isp_aerotrans, & ! number of transported aerosol variables
    cgas         , & ! gas phase concentration                        ( ppm )
    cgas_s       , & ! surface gas phase concentration                ( ppm )
    caero        , & ! aerosol concentration                          (ug m-3 )
    caero_s      , & ! surface aerosol concentration                  (ug m-3 )
    cgastens     , & ! cgas-tendency without sound-wave terms         (ppm/s)
    caerotens    , & ! caero-tendency without sound-wave terms     (ug m-3/s)
    cgast_diff   , & ! tendency by horizontal diffusion
    caerot_diff  , & ! tendency by horizontal diffusion
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

! coefficients used in complete_tendencies
! ----------------------------------------
REAL    (KIND=ireals   ) ::  &
  zfac_qc,              & !
  za1t_surf, za2t_surf    !

REAL (KIND=ireals), ALLOCATABLE ::  &
  zsqrtgrho_r_s(:,:,:), & ! reciprocal square root of G * rho at skalar points
  zsqrtgrho_r_u(:,:,:), & ! reciprocal square root of G * rho at u points
  zsqrtgrho_r_v(:,:,:), & ! reciprocal square root of G * rho at v points
  zsqrtgrho_r_w(:,:,:), & ! reciprocal square root of G * rho at w points
!  zqit_hd(:,:,:),       & !
  za1t(:),              & !
  za2t(:),              & !
  zpia(:,:,:),          & !
  zpianf(:,:,:),        & !
  ztheta(:,:,:),        & !
  ztheta_l(:,:,:),      & !
  ztmkvm(:,:,:),        & !
  ztmkvh(:,:),          & !
  ztch(:,:),            & !
#ifdef COUP_OAS_COS
  ztcw(:,:),            & ! CPS
#endif
  ztcm(:,:),            & !
  zkh(:,:,:),           & !
  ztmkvw(:,:,:),        & !
  ztmkvtke(:,:,:)         !

LOGICAL                  ::  &
  lvertdiff,            & ! .TRUE. if vertical diffusion is calculated here
  lvertdiff_w,          & ! .TRUE. for vertical diffusion of w
  lmoist_turb,          & ! .TRUE. if moist turb. param. (Postdam) is used
  lmassf_diffusion        !

!==============================================================================

CONTAINS

!==============================================================================
!+ init procedure for completing the time step precalculate some coeffs.
!------------------------------------------------------------------------------

SUBROUTINE complete_tendencies_init

!------------------------------------------------------------------------------
!
! Description:
!   This procedure provides some coefficients needed in the computation of 
!   the final slow tendencies for the dynamic variables.
!
!------------------------------------------------------------------------------

! Declarations:

! Local scalars:
! -------------
INTEGER (KIND=iintegers) ::  &
  i,  j,  k              !  Loop indices in lon., lat. and vert. direction


REAL    (KIND=ireals   ) ::  &
  ztmkvtke_tmp(ie,je), & !
  ztvb, zvbke,         & !
  zdr, zpp,            & !
  zfpi, zppa,          & !
  zlhvocp,             & !
  ztkvz, ztkvl

! Statement function zfpi for calculation of the exner function
! where the dummy argument zppa is pressure

zfpi(zppa) = (1.E-5_ireals*zppa)**rdocp

! End of header
!==============================================================================

!------------------------------------------------------------------------------
! Begin Subroutine init_complete_tendencies
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Section 1: Some preparations
!------------------------------------------------------------------------------

! Precalculate some help variables to avoid divisions in subsequent code
! zedsqrtgrho = 1/(sqrt(G)*rho)
!------------------------------------------------------------------------------

  DO  k = 1, ke
    DO j = jstart, jend
      DO i = istart, iend
        zsqrtgrho_r_s(i,j,k) = sqrtg_r_s(i,j,k) / rho(i,j,k)
      ENDDO
    ENDDO
  ENDDO

  DO  k = 1, ke
    DO j = jstartu, jendu
      DO i = istartu, iendu
        zsqrtgrho_r_u(i,j,k) = 2.0_ireals*sqrtg_r_u(i,j,k)              &
          / ( rho(i,j,k)+rho(i+1,j,k) )
      ENDDO
    ENDDO
  ENDDO

  DO  k = 1, ke
    DO j = jstartv, jendv
      DO i = istartv, iendv
        zsqrtgrho_r_v(i,j,k) = 2.0_ireals*sqrtg_r_v(i,j,k)              &
          / ( rho(i,j,k)+rho(i,j+1,k) )
      ENDDO
    ENDDO
  ENDDO

  IF ( lvertdiff_w ) THEN
    DO  k = 2, ke
      DO j = jstart, jend
        DO i = istart, iend
          zsqrtgrho_r_w(i,j,k) = 2.0_ireals*sqrtg_r_w(i,j,k)            &
            / ( rho(i,j,k)+rho(i,j,k-1) )
        ENDDO
      ENDDO
    ENDDO
  END IF

  ! Setting of parameters for horizontal averaging of vertical diffusion
  ! coefficients: ztkvz, ztkvl
  ztkvz = 0.9
  ztkvl = (1.-ztkvz)*0.25

  ! In order to switch off the calculation of vertical diffusion,
  ! we use a local copy for the implicit weights a1t and a2t with values
  ! set to zero. This has to be done if the tendencies have alreay been
  ! calculated in turbdiff. An exception is the cloud ice, the tendencies
  ! have not beec calculated in turbdiff.
  ! This will be removed when a final form of time-integration has been found.

  IF ( ltur .AND. (itype_turb /= 3 .OR. imode_turb < 2) ) THEN
    DO k = 1, ke1
      za1t(k) = a1t(k)
      za2t(k) = a2t(k)
    ENDDO
    lvertdiff = .TRUE.
  ELSE
    DO k = 1, ke1
      za1t(k) = 0.0_ireals
      za2t(k) = 0.0_ireals
    ENDDO
    lvertdiff = .FALSE.
  END IF

  ! Selection of lower boundary conditions

  IF (imode_turb == 0) THEN
    ! condition of the lower boundary is the surface concentration
    za1t_surf = za1t(ke1)  ! implicit weight for ke-value
    za2t_surf = za2t(ke1)  ! explicit weight for ke-value
    zfac_qc   = 1.0_ireals ! zero-qc boundary condition
  ELSE
    ! condition of the lower boundary is the explicit surface mass flux
    za1t_surf = 0.0_ireals ! no implicit weight
    za2t_surf = 1.0_ireals ! full explicit weight
    zfac_qc   = 0.0_ireals ! zero qc-flux boundary condition
  ENDIF

  ! Calculation of Exner-pressure and potential temperature
  ! -------------------------------------------------------
  DO  k = 1 , ke
    DO  j = jstart-1, jend+1
      DO  i = istart-1, iend+1
        zpp           = p0(i,j,k) + pp(i,j,k,nnow)
        zpia(i,j,k)   = zfpi( zpp )
        ztheta(i,j,k) = t(i,j,k,nnow) / zpia(i,j,k)
      ENDDO
    ENDDO
  ENDDO

  ! Some preparations for moist turbulence formulation of sensible heat flux
  ! ------------------------------------------------------------------------
  IF ( lmoist_turb .AND. l3dturb ) THEN
    
    zlhvocp = lh_v / cp_d
    
    ! liquid water potential temperature
    DO  k = 1 , ke
      DO  j = jstart-1, jend+1
        DO  i = istart-1, iend+1
          ztheta_l(i,j,k) = ztheta(i,j,k)                               &
                          - zlhvocp * qc(i,j,k,nnow) / zpia(i,j,k)
        ENDDO
      ENDDO
    ENDDO
    
  END IF
  
  ! Calculation of modified transfer coefficients 
  ! ---------------------------------------------
  DO   j  =  jstart, jend+1
    DO i  =  istart, iend+1
      zvbke      = 0.5*SQRT ( (u(i,j,ke,nnow) + u(i-1,j,ke,nnow))**2    &
                            + (v(i,j,ke,nnow) + v(i,j-1,ke,nnow))**2 )
! UB>> muss hier nicht t_g anstelle von t_s hin???
      ztvb       = t_s (i,j,nnow)*(1.0 + rvd_m_o*qv_s(i,j,nnow))
!      ztvb       = t_g (i,j,nnow)*(1.0 + rvd_m_o*qv_s(i,j,nnow))
! UB<<
      ztcm(i,j)  = tcm(i,j)*zvbke*ps(i,j,nnow)/(r_d*ztvb)
! UB>> limit the horizontal velocity zvbke to a max. of vel_min, to prevent 0 heat/moisture fluxes at 0 windspeed:
!  (NOT ACTIVATED FOR NOW) --> should also be done consistently in the rest
!  of the code, i.e., src_soil_multlay.f90
!!!      ztch(i,j)  = tch(i,j)*zvbke*ps(i,j,nnow)/(r_d*ztvb)
      ztch(i,j)  = tch(i,j)*MAX(zvbke,vel_min)*ps(i,j,nnow)/(r_d*ztvb)
#ifdef COUP_OAS_COS
      ztcw(i,j)  = tcw(i,j)*MAX(zvbke,vel_min)*ps(i,j,nnow)/(r_d*ztvb)  !CPS
#endif
! UB<<
    ENDDO
  ENDDO


  ! Calculation of the modified turbulent diffusion coefficients
  ! ------------------------------------------------------------

  DO  k = 2 , ke
    DO  j = jstart-1, jend+1
      DO  i = istart-1, iend+1
        zdr    = 0.5_ireals*( rho(i,j,k) + rho(i,j,k-1) ) * sqrtg_r_w(i,j,k)
        ztmkvh(i,j  ) = tkvh(i,j,k) * zdr
        ztmkvm(i,j,k) = tkvm(i,j,k) * zdr
        zpianf(i,j,k) = 0.5_ireals * ( zpia(i,j,k) + zpia(i,j,k-1) )
      ENDDO
    ENDDO
    DO    j = jstart, jend
      DO  i = istart, iend
        zkh(i,j,k) = ztkvz*ztmkvh(i,j)                       &
                   + ztkvl*( ztmkvh(i-1,j) + ztmkvh(i+1,j)   &
                           + ztmkvh(i,j-1) + ztmkvh(i,j+1) )
        ztmkvw  (i,j,k) = ztkvz*ztmkvm(i,j,k)                                 &
                        + ztkvl*( ztmkvm(i+1,j,k) + ztmkvm(i,j+1,k)           &
                                + ztmkvm(i-1,j,k) + ztmkvm(i,j-1,k) )
      ENDDO
    ENDDO
  ENDDO

  IF ( lprog_tke ) THEN
    DO  k = 2 , ke
      DO  j = jstart-1, jend+1
        DO  i = istart-1, iend+1
          ztmkvtke_tmp(i,j) = tkvm(i,j,k) * sqrtg_r_w(i,j,k)
        ENDDO
      ENDDO
      DO    j = jstart, jend
        DO  i = istart, iend
          ztmkvtke(i,j,k) = ztkvz*ztmkvtke_tmp(i,j)                             &
                          + ztkvl*( ztmkvtke_tmp(i+1,j) + ztmkvtke_tmp(i,j+1)   &
                                  + ztmkvtke_tmp(i-1,j) + ztmkvtke_tmp(i,j-1) )
        ENDDO
      ENDDO
    ENDDO
  END IF
  
  DO    j = jstart-1, jend+1
    DO  i  = istart-1, iend+1
      zpianf(i,j,ke1) = zfpi( ps(i,j,nnow) )
    ENDDO
  ENDDO

END SUBROUTINE complete_tendencies_init

!==============================================================================
!==============================================================================
!+ procedure for completing the time step for u, v, w, t and pp 
!+ in combination with calculation of the diffusion at the beginning
!------------------------------------------------------------------------------

SUBROUTINE complete_tendencies_uvwtpp( nadv, dtadv, opt_ark, opt_brk )

!------------------------------------------------------------------------------
!
! Description:
!   This procedure calculates the vertical tendencies for the 
!   prognostic variables u,v,w,pp and T.
!
! Method:
!   Using the previously calculated tendencies from horizontal advection,
!   which are stored on the tendency arrays, and tendencies due to 
!   physics and adiabatic processes of the dynamic variables, 
!   the vertical advection is solved by a vertically implicit 
!   scheme (modified Crank-Nicolson).
!
!------------------------------------------------------------------------------

! Declarations:

! Subroutine arguments:
! ---------------------

INTEGER (KIND=iintegers), INTENT(IN) ::  &
  nadv

REAL (KIND=ireals), INTENT(IN) ::  &
  dtadv

REAL (KIND=ireals), INTENT(IN), OPTIONAL ::  &
  opt_ark, opt_brk


! Local scalars:
! -------------
INTEGER (KIND=iintegers) ::  &
  i,  j,  k              !  Loop indices in lon., lat. and vert. direction


REAL    (KIND=ireals   ) ::  &
  zbetav, zbetp, zbetm,& !
  zgav  , zgcv ,       & !
  zag, zas,            & !
  zcg, zcs,            & !
  zbg, zdg,            & !
  zd1g, zd2g,          & !
  znew, zz,            & !
  zdr, zzdtr,          & !
  zark, zbrk             !

! Local (automatic) arrays:
! ------------------------
REAL    (KIND=ireals   ) ::  &
  zgavx   (ie,je,ke),    & !
  zgcvx   (ie,je,ke),    & !
  zc      (ie,je,ke),    & ! Upper main diagonal of ...
  zd1     (ie,je,ke),    & ! Right hand side of ...
  zd2     (ie,je,ke),    & ! Right hand side of ...
  ze1     (ie,je,ke),    & ! Soluton vector of ...
  ze2     (ie,je,ke)       ! Soluton vector of ...

! End of header
!==============================================================================

!------------------------------------------------------------------------------
! Begin Subroutine complete_tendencies_uvwtpp
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Section 1: Some preparations
!------------------------------------------------------------------------------

  ! test if optional parameters are present -
  ! used in combination with TVD-variant of Runge-Kutta scheme
  IF ( PRESENT(opt_ark) .AND. PRESENT(opt_brk) ) THEN
    zark = opt_ark
    zbrk = opt_brk
  ELSE
    zark = 1.0_ireals
    zbrk = 0.0_ireals
  END IF

  ! Setting of parameters for implicit calculation of vertical advection       
  ! (zbetav: weight for the n+1 time level           
  !          zbetav=0: centered, +1: full implicit,-1: explicit)
  zbetav = 0.0
  zbetp  = 0.5*( 1.0 + zbetav )
  zbetm  = 0.5*( 1.0 - zbetav )

  ! Setting of reciprocal time step
  zzdtr = 1.0_ireals/dtadv

! Precalculate some help variables to avoid divisions in subsequent code
!-----------------------------------------------------------------------

!NEC_CB Included below
!  DO  k = 1, ke
!    DO j = jstart, jend
!      DO i = istart, iend
!        zgavx (i,j,k) = - 0.5_ireals*wcon(i,j,k  )
!        zgcvx (i,j,k) =   0.5_ireals*wcon(i,j,k+1)
!      ENDDO
!    ENDDO
!  ENDDO

!------------------------------------------------------------------------------
! Section 2: Setup of tridiagonal matrix systems resulting from the implicit
!            numerical formulation of advection.
!            -  pressure perturbation and temperature
!------------------------------------------------------------------------------

  DO j = jstart, jend

    ! Top layer       
!!CDIR ON_ADB(pp)
!!CDIR ON_ADB(t)
!!CDIR ON_ADB(zc)
!CDIR ON_ADB(zd1)
!CDIR ON_ADB(zd2)
    DO i = istart, iend
      zgcv =   0.5_ireals*wcon(i,j,2)
      zcg  = zgcv*zbetp
      zcs  = zgcv*zbetm
      zbg  = zzdtr - zcg
      zd1g = zzdtr * ( zark*pp(i,j,1,nnow) + zbrk*pp(i,j,1,nadv) )   &
           + pptens(i,j,1) + ppadvt(i,j,1)                           &
           - zcs * ( pp(i,j,2,nadv) - pp(i,j,1,nadv) )
      zd2g = zzdtr * ( zark*t (i,j,1,nnow) + zbrk*t (i,j,1,nadv) )   &
           + ttens (i,j,1) + tadvt (i,j,1)                           &
           - zcs * ( t (i,j,2,nadv) - t (i,j,1,nadv) )
      zd1(i,j,1) = zd1g / zbg
      zd2(i,j,1) = zd2g / zbg
      zc(i,j,1)  = zcg  / zbg
    ENDDO

    ! The layers from k=2 to k=ke-1
    DO k = 2, ke-1
!CDIR ON_ADB(pp)
!CDIR ON_ADB(t)
!CDIR ON_ADB(zc)
!CDIR ON_ADB(zd1)
!CDIR ON_ADB(zd2)
      DO i = istart, iend
        zgav = - 0.5_ireals*wcon(i,j,k  )
        zgcv =   0.5_ireals*wcon(i,j,k+1)
        zag  = zgav*zbetp
        zas  = zgav*zbetm
        zcg  = zgcv*zbetp
        zcs  = zgcv*zbetm
        zbg  = zzdtr - zag - zcg
        zd1g = zzdtr * ( zark*pp(i,j,k,nnow) + zbrk*pp(i,j,k,nadv) )   &
             + pptens(i,j,k) + ppadvt(i,j,k)                           &
             - zas * ( pp(i,j,k-1,nadv) - pp(i,j,k,nadv) )             &
             - zcs * ( pp(i,j,k+1,nadv) - pp(i,j,k,nadv) )
        zd2g = zzdtr * ( zark*t (i,j,k,nnow) + zbrk*t (i,j,k,nadv) )   &
             + ttens (i,j,k) + tadvt (i,j,k)                           &
             - zas * ( t (i,j,k-1,nadv) - t (i,j,k,nadv) )             &
             - zcs * ( t (i,j,k+1,nadv) - t (i,j,k,nadv) )
        zz   = 1.0_ireals / ( zbg - zag*zc(i,j,k-1) )
        zc (i,j,k) = zcg * zz
        zd1(i,j,k) = ( zd1g -zag*zd1(i,j,k-1) ) * zz
        zd2(i,j,k) = ( zd2g -zag*zd2(i,j,k-1) ) * zz
      ENDDO
    ENDDO

    ! The bottom layer
!CDIR ON_ADB(pp)
!CDIR ON_ADB(t)
!CDIR ON_ADB(zc)
!CDIR ON_ADB(zd1)
!CDIR ON_ADB(zd2)
!CDIR ON_ADB(ze1)
!CDIR ON_ADB(ze2)
    DO i = istart, iend
      zgav = - 0.5_ireals*wcon(i,j,ke)
      zag  = zgav*zbetp 
      zas  = zgav*zbetm 
      zbg  = zzdtr - zag
      zd1g = zzdtr * ( zark*pp(i,j,ke,nnow) + zbrk*pp(i,j,ke,nadv) )   &
           + pptens(i,j,ke) + ppadvt(i,j,ke)                           &
           - zas * ( pp(i,j,ke-1,nadv) - pp(i,j,ke,nadv) )
      zd2g = zzdtr * ( zark*t (i,j,ke,nnow) + zbrk*t (i,j,ke,nadv) )   &
           + ttens (i,j,ke) + tadvt (i,j,ke)                           &
           - zas * ( t (i,j,ke-1,nadv) - t (i,j,ke,nadv) )
      zz   = 1.0_ireals / ( zbg - zag*zc(i,j,ke-1) )
      znew = ( zd1g - zag*zd1(i,j,ke-1) ) * zz
      ppadvt(i,j,ke) =  &
        ( znew - zark*pp(i,j,ke,nnow) - zbrk*pp(i,j,ke,nadv) ) * zzdtr
      ze1(i,j,ke) = znew
      znew = ( zd2g - zag*zd2(i,j,ke-1) ) * zz
      tadvt (i,j,ke) =  &
        ( znew - zark*t (i,j,ke,nnow) - zbrk*t (i,j,ke,nadv) ) * zzdtr
      ze2(i,j,ke) = znew
    ENDDO

    ! Backsubstitution and storage of the complete slow tendencies

    DO k = ke-1, 1, -1
!CDIR ON_ADB(ze1)
!CDIR ON_ADB(ze2)
      DO i = istart, iend
        ze1(i,j,k)     = zd1(i,j,k) - zc(i,j,k)*ze1(i,j,k+1)
        ppadvt(i,j,k) =  &
          ( ze1(i,j,k) - zark*pp(i,j,k,nnow) - zbrk*pp(i,j,k,nadv) ) * zzdtr
        ze2(i,j,k)     = zd2(i,j,k) - zc(i,j,k)*ze2(i,j,k+1)
        tadvt (i,j,k) =  &
          ( ze2(i,j,k) - zark*t (i,j,k,nnow) - zbrk*t (i,j,k,nadv) ) * zzdtr
      ENDDO
    ENDDO

  ENDDO ! j-loop

!------------------------------------------------------------------------------
! Section 3: Setup of tridiagonal matrix systems resulting from the implicit
!            numerical formulation of advection.
!            -  vertical velocity
!------------------------------------------------------------------------------

!NEC_CB Included below
!  DO  k = 2, ke
!    ! Precalculate some help variables to avoid divisions in subsequent code
!    DO j = jstart, jend
!      DO i = istart, iend
!        zgavx(i,j,k) =  &
!          - 0.25_ireals*( wcon(i,j,k) + wcon(i,j,k-1) )
!        zgcvx(i,j,k) =  &
!          + 0.25_ireals*( wcon(i,j,k) + wcon(i,j,k+1) )
!      ENDDO
!    ENDDO
!  ENDDO

  ! Top layer       
  DO j = jstart, jend
!CDIR ON_ADB(w)
!CDIR ON_ADB(zc)
!CDIR ON_ADB(zd1)
    DO i = istart, iend
      zgav       = - 0.25_ireals*( wcon(i,j,2) + wcon(i,j,1) )
      zgcv       = + 0.25_ireals*( wcon(i,j,2) + wcon(i,j,3) )
      zcg        = zgcv*zbetp
      zbg        = zzdtr - zcg - zgav*zbetp
      zdg        = zzdtr * ( zark*w(i,j,2,nnow) + zbrk*w(i,j,2,nadv) )   &
        + wtens(i,j,2) + wadvt(i,j,2)                                    &
        - zgcv*zbetm * ( w(i,j,3,nadv) - w(i,j,2,nadv) )         &
        + zgav*zbetm*w(i,j,2,nadv) - zgav*w(i,j,1,nadv)
      zd1(i,j,2) = zdg / zbg
      zc (i,j,2) = zcg / zbg
    ENDDO

    ! The layers from k=3 to k=ke-1
    DO k = 3, ke-1
!CDIR ON_ADB(w)
!CDIR ON_ADB(zc)
!CDIR ON_ADB(zd1)
      DO i = istart, iend
        zgav       = - 0.25_ireals*( wcon(i,j,k) + wcon(i,j,k-1) )
        zgcv       = + 0.25_ireals*( wcon(i,j,k) + wcon(i,j,k+1) )
        zag        = zgav*zbetp
        zas        = zgav*zbetm
        zcg        = zgcv*zbetp
        zcs        = zgcv*zbetm
        zbg        = zzdtr - zag - zcg
        zdg        = zzdtr * ( zark*w(i,j,k,nnow) + zbrk*w(i,j,k,nadv) )  &
          + wtens(i,j,k) + wadvt(i,j,k)                                   &
          - zas * ( w(i,j,k-1,nadv) - w(i,j,k,nadv) )                     &
          - zcs * ( w(i,j,k+1,nadv) - w(i,j,k,nadv) )
        zz         = 1.0_ireals/( zbg - zag*zc(i,j,k-1) )
        zc (i,j,k) = zcg * zz
        zd1(i,j,k) = ( zdg - zag*zd1(i,j,k-1) ) * zz
      ENDDO
    ENDDO

    ! The bottom layer
!CDIR ON_ADB(ze1)
    DO i = istart, iend
      zgav       = - 0.25_ireals*( wcon(i,j,ke) + wcon(i,j,ke-1) )
      zgcv       = + 0.25_ireals*( wcon(i,j,ke) + wcon(i,j,ke+1) )
      zag        = zgav*zbetp 
      zbg        = zzdtr - zag - zgcv*zbetp
      zdg        = zzdtr * ( zark*w(i,j,ke,nnow) + zbrk*w(i,j,ke,nadv) )  &
        + wtens(i,j,ke) + wadvt(i,j,ke)                                   &
        - zgav*zbetm*( w(i,j,ke-1,nadv) - w(i,j,ke,nadv) )       &
        + zgcv*zbetm*w(i,j,ke,nadv)                              &
        - zgcv*w(i,j,ke1,nadv)
      znew       = ( zdg -zag*zd1(i,j,ke-1) ) / ( zbg - zag*zc(i,j,ke-1) )
      wadvt(i,j,ke) =  &
        ( znew - zark*w(i,j,ke,nnow) - zbrk*w(i,j,ke,nadv) ) * zzdtr
      ze1(i,j,ke) = znew
    ENDDO

    ! Backsubstitution and storage of the complete slow tendencies

    DO k = ke-1, 2, -1
!CDIR ON_ADB(ze1)
      DO i = istart, iend
        ze1(i,j,k)     = zd1(i,j,k) - zc(i,j,k)*ze1(i,j,k+1)
        wadvt (i,j,k) =  &
          ( ze1(i,j,k) - zark*w(i,j,k,nnow) - zbrk*w(i,j,k,nadv) ) * zzdtr
      ENDDO
    ENDDO
  ENDDO ! j-loop

!------------------------------------------------------------------------------
! Section 4: Setup of tridiagonal matrix systems resulting from the implicit
!            numerical formulation of advection.
!            -  horizontal wind velocity u
!------------------------------------------------------------------------------

  ! Top layer:   k = 1
  DO j = jstartu, jendu
!CDIR ON_ADB(wcon)
!CDIR ON_ADB(u)
!CDIR ON_ADB(zc)
!CDIR ON_ADB(zd1)
    DO i = istartu, iendu
      zgcv       = 0.25_ireals*( wcon(i,j,2)+wcon(i+1,j,2) )
      zcg        = zgcv*zbetp
      zbg        = zzdtr - zcg
      zdg        = zzdtr * ( zark*u(i,j,1,nnow) + zbrk*u(i,j,1,nadv) )  &
        + utens(i,j,1) + uadvt(i,j,1)                                   &
        - zgcv*zbetm * ( u(i,j,2,nadv) - u(i,j,1,nadv) )
      zd1(i,j,1) = zdg / zbg
      zc (i,j,1) = zcg / zbg
    ENDDO

    ! The layers from k=2 to k=ke-1
    DO k = 2, ke-1
!CDIR ON_ADB(wcon)
!CDIR ON_ADB(u)
!CDIR ON_ADB(zc)
!CDIR ON_ADB(zd1)
      DO i = istartu, iendu
        zgav       =  &
          -0.25_ireals*( wcon(i,j,k  )+wcon(i+1,j,k  ) )
        zgcv       =  &
           0.25_ireals*( wcon(i,j,k+1)+wcon(i+1,j,k+1) )
        zag        = zgav*zbetp
        zas        = zgav*zbetm
        zcg        = zgcv*zbetp
        zcs        = zgcv*zbetm
        zbg        = zzdtr - zag - zcg
        zdg        = zzdtr * ( zark*u(i,j,k,nnow) + zbrk*u(i,j,k,nadv) )  &
          + utens(i,j,k) + uadvt(i,j,k)                                   &
          - zas * ( u(i,j,k-1,nadv) - u(i,j,k,nadv) )                     &
          - zcs * ( u(i,j,k+1,nadv) - u(i,j,k,nadv) )
        zz         = 1.0_ireals / ( zbg - zag*zc(i,j,k-1) )
        zc (i,j,k) = zcg * zz
        zd1(i,j,k) = ( zdg -zag*zd1(i,j,k-1) ) * zz
      ENDDO
    ENDDO

    ! The bottom layer:  k = ke
!CDIR ON_ADB(ze1)
    DO i = istartu, iendu
      zgav       = - 0.25_ireals*( wcon(i,j,ke)+wcon(i+1,j,ke) )
      zag        = zgav*zbetp
      zas        = zgav*zbetm
      zbg        = zzdtr - zag
      zdg        = zzdtr * ( zark*u(i,j,ke,nnow) + zbrk*u(i,j,ke,nadv) )  &
        + utens(i,j,ke) + uadvt(i,j,ke)                                   &
        - zas * ( u(i,j,ke-1,nadv) - u(i,j,ke,nadv) )
      znew       = ( zdg -zag*zd1(i,j,ke-1) ) / ( zbg - zag*zc(i,j,ke-1) )
      uadvt (i,j,ke) =  &
        ( znew - zark*u(i,j,ke,nnow) - zbrk*u(i,j,ke,nadv) ) * zzdtr
      ze1(i,j,ke) = znew
    ENDDO

    ! Backsubstitution and storage of the complete slow tendencies

    DO k = ke-1, 1, -1
!CDIR ON_ADB(ze1)
      DO i = istartu, iendu
        ze1(i,j,k) =  zd1(i,j,k) - zc(i,j,k)*ze1(i,j,k+1)
        uadvt(i,j,k) =  &
          ( ze1(i,j,k) - zark*u(i,j,k,nnow) - zbrk*u(i,j,k,nadv) ) * zzdtr
      ENDDO
    ENDDO
  ENDDO

!------------------------------------------------------------------------------
! Section 5: Setup of tridiagonal matrix systems resulting from the implicit
!            numerical formulation of advection.
!            -  horizontal wind velocity v
!------------------------------------------------------------------------------

  ! Top layer  k=1
  DO j = jstartv, jendv
!CDIR ON_ADB(wcon)
!CDIR ON_ADB(v)
!CDIR ON_ADB(zc)
!CDIR ON_ADB(zd1)
    DO i = istartv, iendv
      zgcv       = 0.25_ireals*( wcon(i,j,2)+wcon(i,j+1,2) )
      zcg        = zgcv*zbetp
      zbg        = zzdtr - zcg
      zdg        = zzdtr * ( zark*v(i,j,1,nnow) + zbrk*v(i,j,1,nadv) )  &
        + vtens(i,j,1) + vadvt(i,j,1)                                   &
        - zgcv*zbetm * ( v(i,j,2,nadv) - v(i,j,1,nadv) )
      zd1(i,j,1) = zdg / zbg
      zc (i,j,1) = zcg / zbg
    ENDDO

    ! The layers from k=2 to k=ke-1
    DO  k = 2, ke-1
!CDIR ON_ADB(wcon)
!CDIR ON_ADB(v)
!CDIR ON_ADB(zc)
!CDIR ON_ADB(zd1)
      DO i = istartv, iendv
        zgav       = -0.25_ireals*( wcon(i,j,k  )+wcon(i,j+1,k  ) )
        zgcv       =  0.25_ireals*( wcon(i,j,k+1)+wcon(i,j+1,k+1) )
        zag        = zgav*zbetp
        zas        = zgav*zbetm
        zcg        = zgcv*zbetp
        zcs        = zgcv*zbetm 
        zbg        = zzdtr - zag - zcg
        zdg        = zzdtr * ( zark*v(i,j,k,nnow) + zbrk*v(i,j,k,nadv) )  &
          + vtens(i,j,k) + vadvt(i,j,k)                                   &
          - zas * ( v(i,j,k-1,nadv) - v(i,j,k,nadv) )                     &
          - zcs * ( v(i,j,k+1,nadv) - v(i,j,k,nadv) )
        zz         = 1.0_ireals / ( zbg - zag*zc(i,j,k-1) )
        zc (i,j,k) = zcg * zz
        zd1(i,j,k) = ( zdg -zag*zd1(i,j,k-1) ) * zz
      ENDDO
    ENDDO

    ! The bottom layer k=ke
!CDIR ON_ADB(ze1)
    DO i = istartv, iendv
      zgav       = - 0.25_ireals*( wcon(i,j,ke)+wcon(i,j+1,ke) )
      zag        = zgav*zbetp
      zas        = zgav*zbetm
      zbg        = zzdtr - zag
      zdg        = zzdtr * ( zark*v(i,j,ke,nnow) + zbrk*v(i,j,ke,nadv) )  &
        + vtens(i,j,ke) + vadvt(i,j,ke)                                   &
        - zas * ( v(i,j,ke-1,nadv) - v(i,j,ke,nadv) )
      znew       = ( zdg -zag*zd1(i,j,ke-1) ) / ( zbg - zag*zc(i,j,ke-1) )
      vadvt (i,j,ke) =  &
        ( znew - zark*v(i,j,ke,nnow) - zbrk*v(i,j,ke,nadv) ) * zzdtr
      ze1(i,j,ke) = znew
    ENDDO

    ! Backsubstitution and storage of the complete slow tendencies

    DO k = ke-1, 1, -1
!CDIR ON_ADB(ze1)
      DO i = istartv, iendv
        ze1(i,j,k) = zd1(i,j,k) - zc(i,j,k)*ze1(i,j,k+1)
        vadvt(i,j,k) =  &
          ( ze1(i,j,k) - zark*v(i,j,k,nnow) - zbrk*v(i,j,k,nadv) ) * zzdtr
      ENDDO
    ENDDO
  ENDDO

!------------------------------------------------------------------------------
! End of subroutine complete_tendencies_uvwtpp
!------------------------------------------------------------------------------

END SUBROUTINE complete_tendencies_uvwtpp

!==============================================================================
!==============================================================================
!+ procedure for completing the time step for u, v, w, t and pp
!+ in combination with vertical explicit advection and
!+ calculation of the vertical diffusion at the beginning
!------------------------------------------------------------------------------

!option! -pvctl _on_adb
SUBROUTINE complete_tendencies_uvwtpp_eva

!------------------------------------------------------------------------------
!
! Description:
!   This procedure provides the final slow tendencies for the dynamic 
!   variables.
!
!------------------------------------------------------------------------------

! Declarations:

! Local scalars:
! -------------
INTEGER (KIND=iintegers) ::  &
  i,  j,  k              !  Loop indices in lon., lat. and vert. direction

! End of header
!==============================================================================

!------------------------------------------------------------------------------
! Begin Subroutine complete_tendencies_uvwtpp_eva
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Section : Pressure perturbation
!------------------------------------------------------------------------------

  DO k = 1, ke
    DO j = jstart, jend
      DO i = istart, iend
        ppadvt(i,j,k) = ppadvt(i,j,k) + pptens(i,j,k)
      ENDDO
    ENDDO
  ENDDO

!------------------------------------------------------------------------------
! Section : Vertical velocity
!------------------------------------------------------------------------------

  DO k = 2, ke
    DO j = jstart, jend
      DO i = istart, iend
        wadvt(i,j,k) = wadvt(i,j,k) + wtens(i,j,k)
      ENDDO
    ENDDO
  ENDDO

!------------------------------------------------------------------------------
! Section : Temperature 
!------------------------------------------------------------------------------

  DO k = 1, ke
    DO j = jstart, jend
      DO i = istart, iend
        tadvt(i,j,k) = tadvt(i,j,k) + ttens(i,j,k)
      ENDDO
    ENDDO
  ENDDO

!------------------------------------------------------------------------------
! Section : Horizontal wind velocity u 
!------------------------------------------------------------------------------

  DO k = 1, ke
    DO j = jstartu, jendu
      DO i = istartu, iendu
        uadvt(i,j,k) = uadvt(i,j,k) + utens(i,j,k)
      ENDDO
    ENDDO
  ENDDO

!------------------------------------------------------------------------------
! Section : Horizontal wind velocity v
!------------------------------------------------------------------------------

  DO k = 1, ke
    DO j = jstartv, jendv
      DO i = istartv, iendv
        vadvt(i,j,k) = vadvt(i,j,k) + vtens(i,j,k)
      ENDDO
    ENDDO
  ENDDO

!------------------------------------------------------------------------------
! End of subroutine complete_tendencies_uvwtpp_eva
!------------------------------------------------------------------------------

END SUBROUTINE complete_tendencies_uvwtpp_eva

!==============================================================================
!==============================================================================
!+ procedure for completing the time step for qv, qc, qi and tke
!------------------------------------------------------------------------------

SUBROUTINE complete_tendencies_qvqcqi_tke

!------------------------------------------------------------------------------
!
! Description:
!   This procedure is only used in src_runge_kutta and calculates the
!   final updates for the prognostic moisture variables qv, qc, tracer (AK) and qi
!   at time level n+1 (nnew). 
!   The vertical diffusion as the last slow tendency 
!   is computed here for all variables the vertical diffusion acts on. 
!
! Method:
!   Using the previously calculated updates (on the timelevel nnew) for the
!   moisture variables from advection, physics and adiabatic processes, the 
!   vertical diffusion is solved by a vertically implicit scheme (modified 
!   Crank-Nicolson). The water vapour, cloud water and cloud ice variables
!   are updated directly. Also, the surface fluxes of moisture are stored.
!
!------------------------------------------------------------------------------

! Declarations:

! Local scalars:
! -------------
INTEGER (KIND=iintegers) ::  &
  i,  j,  k                   !  Loop indices in lon., lat. and vert. direction

INTEGER (KIND=iintegers) :: &
  km1, kp1
  
REAL    (KIND=ireals   ) ::  &
  zbetav, zbetp, zbetm,  & !
  zgat, zgct,            & !
  zag, zas,              & !
  zcg, zcs,              & !
  zbg, zdg,              & !
  zd1g, zd2g,            & !
  zd3g, zd4g,            & !
  znew, zz,              & !
  zd, zzdtr                !

REAL    (KIND=ireals   ) :: help1,help2

! Local (automatic) arrays:
! ------------------------
REAL    (KIND=ireals   ) ::  &
  zgavx   (ie,je,ke),    & !
  zgcvx   (ie,je,ke),    & !
  zwcdz   (ie,je,ke),    & !
  zgatz   (ie,je,ke),    & !
  zgctz   (ie,je,ke),    & !
  zqvcor  (ie,je   ),    & !
  zqccor  (ie,je   ),    & !
  zc      (ie,je,ke),    & ! Upper main diagonal of ...
  zd1     (ie,je,ke),    & ! Right hand side of ...
  zd2     (ie,je,ke),    & ! Right hand side of ...
  zd3     (ie,je,ke),    & ! Right hand side of ...
  zd4     (ie,je,ke),    & ! Right hand side of ...
  ze      (ie,je,ke)       ! Soluton vector of ...

! End of header
!==============================================================================

!------------------------------------------------------------------------------
! Begin Subroutine complete_tendencies_qvqcqi_tke  
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Section 1: Some preparations
!------------------------------------------------------------------------------

  ! Setting of reciprocal time step
  zzdtr = 1.0 / dt

  IF ( lvertdiff ) THEN
  
    !--------------------------------------------------------------------------
    ! Section 2a: Setup of tridiagonal matrix systems resulting from the
    !             implicit numerical formulation of diffusion of qc and qv
    !             and qi (if required).
    !--------------------------------------------------------------------------

    ! First, the matrix elements a(k), b(k) and c(k) of the coefficient
    ! matrix are set (these are the same for qv, qc and qi).
    ! The right hand side is stored on d(k).

    IF ( lprog_qi ) THEN
      
      ! Top layer
      DO j = jstart, jend
!CDIR ON_ADB(qv)
!CDIR ON_ADB(qc)
!CDIR ON_ADB(qi)
!CDIR ON_ADB(ze)
!CDIR ON_ADB(zc)
!CDIR ON_ADB(zd1)
!CDIR ON_ADB(zd2)
!CDIR ON_ADB(zd3)
!CDIR ON_ADB(zd4)
        DO i = istart, iend
          zgct       = - zkh(i,j,2)*zsqrtgrho_r_s(i,j,1)
          zcg        = zgct*za1t(2)
          zcs        = zgct*za2t(2)
          zbg        = zzdtr - zcg
          zd3g       = zzdtr * qv(i,j,1,nnew) + qvt_diff(i,j,1)    &
                     - zcs * ( qv(i,j,2,nnow) - qv(i,j,1,nnow) )
          zd2g       = zzdtr * qc(i,j,1,nnew) + qctens(i,j,1)      &
                     - zcs * ( qc(i,j,2,nnow) - qc(i,j,1,nnow) )
          zd1g       = zd3g + qvtens(i,j,1)
          zd4g       = zzdtr * qi(i,j,1,nnew) + qitens(i,j,1)      &
                     - zcs * ( qi(i,j,2,nnow) - qi(i,j,1,nnow) )
          zc (i,j,1) = zcg / zbg
          zd1(i,j,1) = zd1g / zbg
          zd2(i,j,1) = zd2g / zbg
          zd3(i,j,1) = zd3g / zbg
          zd4(i,j,1) = zd4g / zbg
        ENDDO

        ! The layers from k=2 to k=ke-1
        DO k = 2, ke-1
!CDIR ON_ADB(qv)
!CDIR ON_ADB(qc)
!CDIR ON_ADB(qi)
!CDIR ON_ADB(ze)
!CDIR ON_ADB(zc)
!CDIR ON_ADB(zd1)
!CDIR ON_ADB(zd2)
!CDIR ON_ADB(zd3)
!CDIR ON_ADB(zd4)
          DO i = istart, iend
            zgat       = - zkh(i,j,k  )*zsqrtgrho_r_s(i,j,k)
            zgct       = - zkh(i,j,k+1)*zsqrtgrho_r_s(i,j,k)
            zag        = zgat*za1t(k)
            zas        = zgat*za2t(k)
            zcg        = zgct*za1t(k+1)
            zcs        = zgct*za2t(k+1)
            zbg        = zzdtr - zag - zcg
            zd3g       = zzdtr * qv(i,j,k,nnew) + qvt_diff(i,j,k)      &
                       - zas * ( qv(i,j,k-1,nnow) - qv(i,j,k,nnow) )   &
                       - zcs * ( qv(i,j,k+1,nnow) - qv(i,j,k,nnow) )
            zd2g       = zzdtr * qc(i,j,k,nnew) + qctens(i,j,k)        &
                       - zas * ( qc(i,j,k-1,nnow) - qc(i,j,k,nnow) )   &
                       - zcs * ( qc(i,j,k+1,nnow) - qc(i,j,k,nnow) )
            zd1g       = zd3g + qvtens(i,j,k)
            zd4g       = zzdtr * qi(i,j,k,nnew) + qitens(i,j,k)        &
                       - zas * ( qi(i,j,k-1,nnow) - qi(i,j,k,nnow) )   &
                       - zcs * ( qi(i,j,k+1,nnow) - qi(i,j,k,nnow) )
            zz         = 1.0_ireals / ( zbg - zag*zc(i,j,k-1) )
            zc (i,j,k) = zcg * zz
            zd1(i,j,k) = ( zd1g -zag*zd1(i,j,k-1) ) * zz
            zd2(i,j,k) = ( zd2g -zag*zd2(i,j,k-1) ) * zz
            zd3(i,j,k) = ( zd3g -zag*zd3(i,j,k-1) ) * zz
            zd4(i,j,k) = ( zd4g -zag*zd4(i,j,k-1) ) * zz
          ENDDO
        ENDDO

        ! the bottom layer
!CDIR ON_ADB(qv)
!CDIR ON_ADB(qc)
!CDIR ON_ADB(qi)
!CDIR ON_ADB(ze)
        DO i = istart, iend
          zgat       = - zkh (i,j,ke)*zsqrtgrho_r_s(i,j,ke)
#ifdef COUP_OAS_COS
          zgct       = - ztcw(i,j   )*zsqrtgrho_r_s(i,j,ke)
#else
          zgct       = - ztch(i,j   )*zsqrtgrho_r_s(i,j,ke)
#endif
          zag        = zgat*za1t(ke)
          zas        = zgat*za2t(ke)
          zcg        = zgct*za1t_surf
          zcs        = zgct*za2t_surf
          zbg        = zzdtr - zag - zcg
          zd3g       = zzdtr * qv(i,j,ke,nnew) + qvt_diff(i,j,ke)       &
                     - zas * ( qv(i,j,ke-1,nnow) - qv(i,j,ke,nnow) )    &
                     + zcs *  qv(i,j,ke,nnow) - zgct * qv_s(i,j,nnow)
          zd2g       = zzdtr * qc(i,j,ke,nnew) + qctens(i,j,ke)         &
                     - zas * ( qc(i,j,ke-1,nnow) - qc(i,j,ke,nnow) )    &
                     + zcs * zfac_qc * qc(i,j,ke,nnow)
          zd1g       = zd3g + qvtens(i,j,ke)
          zd4g       = zzdtr * qi(i,j,ke,nnew) + qitens(i,j,ke)         &
                     - zas * ( qi(i,j,ke-1,nnow) - qi(i,j,ke,nnow) )    &
                     + zcs * zfac_qc * qi(i,j,ke,nnow)
          zz         = 1.0_ireals / ( zbg  - zag*zc(i,j,ke-1) )
          qv(i,j,ke,nnew) = ( zd1g - zag*zd1(i,j,ke-1) ) * zz
          qc(i,j,ke,nnew) = ( zd2g - zag*zd2(i,j,ke-1) ) * zz
          znew            = ( zd3g - zag*zd3(i,j,ke-1) ) * zz
          dqvdt(i,j,ke)   =  &
            ( MAX(0.0_ireals, znew) - qv(i,j,ke,nnow) )*zzdtr
          ze(i,j,ke)      = znew
          qi(i,j,ke,nnew) = ( zd4g - zag*zd4(i,j,ke-1) ) * zz
        ENDDO

        ! Backsubstitution and storage of the complete slow tendencies
!CDIR ON_ADB(qv)
!CDIR ON_ADB(qc)
!CDIR ON_ADB(qi)
!CDIR ON_ADB(ze)
        DO k = ke-1, 1, -1
          DO i = istart, iend
            qv(i,j,k,nnew) = zd1(i,j,k) - zc(i,j,k)*qv(i,j,k+1,nnew)
            qc(i,j,k,nnew) = zd2(i,j,k) - zc(i,j,k)*qc(i,j,k+1,nnew)
            ze(i,j,k)      = zd3(i,j,k) - zc(i,j,k)*ze(i,j,k+1)
            dqvdt(i,j,k)   =  &
              ( MAX(0.0_ireals,ze(i,j,k)) - qv(i,j,k,nnow) ) * zzdtr
            qi(i,j,k,nnew) = zd4(i,j,k) - zc(i,j,k)*qi(i,j,k+1,nnew)
          ENDDO
        ENDDO
      ENDDO

    ELSE

      ! Top layer
      DO j = jstart, jend
        DO i = istart, iend
          zgct       = - zkh(i,j,2)*zsqrtgrho_r_s(i,j,1)
          zcg        = zgct*za1t(2)
          zcs        = zgct*za2t(2)
          zbg        = zzdtr - zcg
          zd3g       = zzdtr * qv(i,j,1,nnew) + qvt_diff(i,j,1)    &
                     - zcs * ( qv(i,j,2,nnow) - qv(i,j,1,nnow) )
          zd2g       = zzdtr * qc(i,j,1,nnew) + qctens(i,j,1)      &
                     - zcs * ( qc(i,j,2,nnow) - qc(i,j,1,nnow) )
          zd1g       = zd3g + qvtens(i,j,1)
          zc (i,j,1) = zcg / zbg
          zd1(i,j,1) = zd1g / zbg
          zd2(i,j,1) = zd2g / zbg
          zd3(i,j,1) = zd3g / zbg
        ENDDO
      ENDDO

      ! The layers from k=2 to k=ke-1
      DO k = 2, ke-1
        DO j = jstart, jend
          DO i = istart, iend
            zgat       = - zkh(i,j,k  )*zsqrtgrho_r_s(i,j,k)
            zgct       = - zkh(i,j,k+1)*zsqrtgrho_r_s(i,j,k)
            zag        = zgat*za1t(k)
            zas        = zgat*za2t(k)
            zcg        = zgct*za1t(k+1)
            zcs        = zgct*za2t(k+1)
            zbg        = zzdtr - zag - zcg
            zd3g       = zzdtr * qv(i,j,k,nnew) + qvt_diff(i,j,k)      &
                       - zas * ( qv(i,j,k-1,nnow) - qv(i,j,k,nnow) )   &
                       - zcs * ( qv(i,j,k+1,nnow) - qv(i,j,k,nnow) )
            zd2g       = zzdtr * qc(i,j,k,nnew) + qctens(i,j,k)        &
                       - zas * ( qc(i,j,k-1,nnow) - qc(i,j,k,nnow) )   &
                       - zcs * ( qc(i,j,k+1,nnow) - qc(i,j,k,nnow) )
            zd1g       = zd3g + qvtens(i,j,k)
            zz         = 1.0_ireals / ( zbg - zag*zc(i,j,k-1) )
            zc (i,j,k) = zcg * zz
            zd1(i,j,k) = ( zd1g -zag*zd1(i,j,k-1) ) * zz
            zd2(i,j,k) = ( zd2g -zag*zd2(i,j,k-1) ) * zz
            zd3(i,j,k) = ( zd3g -zag*zd3(i,j,k-1) ) * zz
          ENDDO
        ENDDO
      ENDDO

      ! the bottom layer
      DO j = jstart, jend
        DO i = istart, iend
          zgat       = - zkh (i,j,ke)*zsqrtgrho_r_s(i,j,ke)
#ifdef COUP_OAS_COS
          zgct       = - ztcw(i,j   )*zsqrtgrho_r_s(i,j,ke)
#else
          zgct       = - ztch(i,j   )*zsqrtgrho_r_s(i,j,ke)
#endif
          zag        = zgat*za1t(ke)
          zas        = zgat*za2t(ke)
          zcg        = zgct*za1t_surf
          zcs        = zgct*za2t_surf
          zbg        = zzdtr - zag - zcg
          zd3g       = zzdtr * qv(i,j,ke,nnew) + qvt_diff(i,j,ke)       &
                     - zas * ( qv(i,j,ke-1,nnow) - qv(i,j,ke,nnow) )    &
                     + zcs *  qv(i,j,ke,nnow) - zgct * qv_s(i,j,nnow)
          zd2g       = zzdtr * qc(i,j,ke,nnew) + qctens(i,j,ke)         &
                     - zas * ( qc(i,j,ke-1,nnow) - qc(i,j,ke,nnow) )    &
                     + zcs * zfac_qc * qc(i,j,ke,nnow)
          zd1g       = zd3g + qvtens(i,j,ke)
          zz         = 1.0_ireals / ( zbg  - zag*zc(i,j,ke-1) )
          qv(i,j,ke,nnew) = ( zd1g - zag*zd1(i,j,ke-1) ) * zz
          qc(i,j,ke,nnew) = ( zd2g - zag*zd2(i,j,ke-1) ) * zz
          znew            = ( zd3g - zag*zd3(i,j,ke-1) ) * zz
          dqvdt(i,j,ke)   =  &
            ( MAX(0.0_ireals, znew) - qv(i,j,ke,nnow) )*zzdtr
          ze(i,j,ke)      = znew
        ENDDO
      ENDDO

      ! Backsubstitution and storage of the complete slow tendencies
      DO k = ke-1, 1, -1
        DO j = jstart, jend
          DO i = istart, iend
            qv(i,j,k,nnew) = zd1(i,j,k) - zc(i,j,k)*qv(i,j,k+1,nnew)
            qc(i,j,k,nnew) = zd2(i,j,k) - zc(i,j,k)*qc(i,j,k+1,nnew)
            ze(i,j,k)      = zd3(i,j,k) - zc(i,j,k)*ze(i,j,k+1)
            dqvdt(i,j,k)   =  &
              ( MAX(0.0_ireals,ze(i,j,k)) - qv(i,j,k,nnow) ) * zzdtr
          ENDDO
        ENDDO
      ENDDO

    END IF

  ELSE
    
    DO k = 1, ke
      DO j = jstart, jend
        DO i = istart, iend
          qv(i,j,k,nnew) = qv(i,j,k,nnew) + dt * qvt_diff(i,j,k)
          dqvdt(i,j,k)   =  &
            ( MAX(0.0_ireals, qv(i,j,k,nnew)) - qv(i,j,k,nnow) )*zzdtr
          qv(i,j,k,nnew) = qv(i,j,k,nnew) + dt * qvtens(i,j,k)
          qc(i,j,k,nnew) = qc(i,j,k,nnew) + dt * qctens(i,j,k)
        ENDDO
      ENDDO
    ENDDO
    
  END IF

  !--------------------------------------------------------------------------
  ! Section 2c: Do vertical diffusion for cloud ice (qi), if required.
  !             setup of tridiagonal matrix systems resulting from the
  !             implicit numerical formulation of diffusion.
  !--------------------------------------------------------------------------

  IF ( lprog_qi .AND. .NOT.lvertdiff ) THEN
    
    ! First, the matrix elements a(k), b(k) and c(k) of the coefficient matrix
    ! are set (for vert. diff. only becaus vertical advection has already
    ! been computed.  The right hand side is stored on zd.

    ! Top layer
    DO j = jstart, jend
      DO i = istart, iend
        zgct       = - zkh(i,j,2)*zsqrtgrho_r_s(i,j,1)
        zcg        = zgct*a1t(2)
        zcs        = zgct*a2t(2)
        zbg        = zzdtr - zcg
        zdg        = zzdtr * qi(i,j,1,nnew) + qitens(i,j,1)   &
          - zcs * (qi(i,j,2,nnow) - qi(i,j,1,nnow))
        zd1(i,j,1) = zdg/zbg
        zc (i,j,1) = zcg/zbg
      ENDDO
    ENDDO

    ! The layers from k=2 to k=ke-1
    DO k = 2, ke-1
      DO j = jstart, jend
        DO i = istart, iend
          zgat       = - zkh(i,j,k  )*zsqrtgrho_r_s(i,j,k)
          zgct       = - zkh(i,j,k+1)*zsqrtgrho_r_s(i,j,k)
          zag        = zgat*a1t(k)
          zas        = zgat*a2t(k)
          zcg        = zgct*a1t(k+1)
          zcs        = zgct*a2t(k+1)
          zbg        = zzdtr - zag - zcg
          zdg        = zzdtr * qi(i,j,k,nnew) + qitens(i,j,k)   &
            - zas* (qi(i,j,k-1,nnow) - qi(i,j,k,nnow))          &
            - zcs* (qi(i,j,k+1,nnow) - qi(i,j,k,nnow))
          zz         = 1.0_ireals/( zbg - zag*zc(i,j,k-1) )
          zc (i,j,k) = zcg*zz
          zd1(i,j,k) = ( zdg -zag*zd1(i,j,k-1) )*zz
        ENDDO
      ENDDO
    ENDDO

    ! the bottom layer
    DO j = jstart, jend
      DO i = istart, iend
        zgat       = - zkh(i,j,ke)*zsqrtgrho_r_s(i,j,ke)
#ifdef COUP_OAS_COS
        zgct       = - ztcw(i,j )*zsqrtgrho_r_s(i,j,ke)
#else
        zgct       = - ztch(i,j )*zsqrtgrho_r_s(i,j,ke)
#endif
        zag        = zgat*za1t(ke)
        zas        = zgat*za2t(ke)
        zcg        = zgct*za1t_surf
        zcs        = zgct*za2t_surf*zfac_qc   
        zbg        = zzdtr - zag - zcg
        zdg        = zzdtr * qi(i,j,ke,nnew) + qitens(i,j,ke)   &
          - zas * (qi(i,j,ke-1,nnow) - qi(i,j,ke,nnow))         &
          + zcs * qi(i,j,ke,nnow)
        zz         = 1.0_ireals / ( zbg  - zag*zc(i,j,ke-1) )
        qi(i,j,ke,nnew) = ( zdg - zag*zd1(i,j,ke-1) ) * zz
      ENDDO
    ENDDO

    ! Backsubstitution 
    DO k = ke-1, 1, -1
      DO j = jstart, jend
        DO i = istart, iend
          qi(i,j,k,nnew) = zd1(i,j,k) - zc(i,j,k)*qi(i,j,k+1,nnew)
        ENDDO
      ENDDO
    ENDDO

  ENDIF

  !----------------------------------------------------------------------------
  ! Section 3: Calculation of the surface moisture flux 'qvsflx'
  !            The latent heat flux is integrated in time
  !----------------------------------------------------------------------------

  DO j =jstart, jend
    DO i = istart , iend
#ifdef COUP_OAS_COS
      qvsflx(i,j)  =  - ztcw(i,j) * &
        ( za2t_surf*(qv_s(i,j,nnow) - qv(i,j,ke,nnow)) +  &
          za1t_surf*(qv_s(i,j,nnew) - qv(i,j,ke,nnew)) )
      lhfl_s(i,j) = lh_v*qvsflx(i,j)
#else
      qvsflx(i,j)  =  - ztch(i,j) * &
        ( za2t_surf*(qv_s(i,j,nnow) - qv(i,j,ke,nnow)) +  &
          za1t_surf*(qv_s(i,j,nnew) - qv(i,j,ke,nnew)) )
      lhfl_s(i,j) = lh_v*qvsflx(i,j)
#endif
    ENDDO
  ENDDO

  !----------------------------------------------------------------------------
  ! Section 4: Massflux correction scheme for removing negative qv and qc
  !            values (if required)
  !----------------------------------------------------------------------------

  IF ( lmassf_diffusion ) THEN

    DO j = jstart, jend
      DO i = istart, iend
        zqvcor(i,j) = 0.0_ireals
        zqccor(i,j) = 0.0_ireals
      ENDDO
    ENDDO
    DO k = 1, ke
      DO j = jstart, jend
        DO i = istart, iend
          help1 = qv(i,j,k,nnew) + zqvcor(i,j)*zsqrtgrho_r_s(i,j,k)
          IF ( help1 < 0.0_ireals ) THEN
            help2 = 0.0_ireals
            help1 = help1/zsqrtgrho_r_s(i,j,k)
          ELSE
            help2 = help1
            help1 = 0.0_ireals
          ENDIF
          zqvcor(i,j)    = help1
          qv(i,j,k,nnew) = help2

          help1 = qc(i,j,k,nnew) + zqccor(i,j)*zsqrtgrho_r_s(i,j,k)
          IF ( help1 < 0.0_ireals ) THEN
            help2 = 0.0_ireals
            help1 = help1/zsqrtgrho_r_s(i,j,k)
          ELSE
            help2 = help1 
            help1 = 0.0_ireals
          ENDIF
          zqccor(i,j)    = help1
          qc(i,j,k,nnew) = help2
        ENDDO
      ENDDO
    ENDDO

  ELSE
    
    CALL clipping( qv(:,:,:,nnew), ie, je, ke )
    CALL clipping( qc(:,:,:,nnew), ie, je, ke )
    
  ENDIF ! end mass-flux correction

  IF ( lprog_qi ) CALL clipping( qi(:,:,:,nnew), ie, je, ke )
  
  !----------------------------------------------------------------------------
  ! Section 5: Do vertical diffusion for turbulent kinetic energy (tke),
  !            if required.
  !----------------------------------------------------------------------------

  IF ( lprog_tke ) THEN
    
    DO  k = 2, ke
      
      km1 = MAX( 2, k-1 )
      kp1 = MIN( ke, k+1 )
      
      ! Precalculate some help variables to avoid divisions in subsequent code
      DO j = jstart, jend
        DO i = istart, iend
          zgatz(i,j,k) = - sqrtg_r_w(i,j,k)   &
                         * ( ztmkvtke(i,j,k)+ztmkvtke(i,j,km1) )
          zgctz(i,j,k) = - sqrtg_r_w(i,j,k)   &
                         * ( ztmkvtke(i,j,kp1)+ztmkvtke(i,j,k) )
        ENDDO
      ENDDO
      
    ENDDO

    ! Top layer       
    DO j = jstart, jend
      DO i = istart, iend
        zag        = 0.5*(a1t(2)+a1t(1)) * zgatz(i,j,2)
        zas        = 0.5*(a2t(2)+a2t(1)) * zgatz(i,j,2)
        zcg        = 0.5*(a1t(3)+a1t(2)) * zgctz(i,j,2)
        zcs        = 0.5*(a2t(3)+a2t(2)) * zgctz(i,j,2)
        zbg        = zzdtr - zag - zcg
        zdg        = zzdtr * tke(i,j,2,nnew) + tketens(i,j,2)                &
                     + zas * tke(i,j,2,nnow)                                 &
                     - zcs * ( tke(i,j,3,nnow) - tke(i,j,2,nnow) )
        zd1(i,j,2) = zdg/zbg
        zc (i,j,2) = zcg/zbg
      ENDDO
    ENDDO

    ! The layers from k=3 to k=ke-1
    DO k = 3, ke-1
      DO j = jstart, jend
        DO i = istart, iend
          zag        = 0.5*(a1t(k)+a1t(k-1)) * zgatz(i,j,k)
          zas        = 0.5*(a2t(k)+a2t(k-1)) * zgatz(i,j,k)
          zcg        = 0.5*(a1t(k+1)+a1t(k)) * zgctz(i,j,k)
          zcs        = 0.5*(a2t(k+1)+a2t(k)) * zgctz(i,j,k)
          zbg        = zzdtr - zag - zcg
          zdg        = zzdtr * tke(i,j,k,nnew) + tketens(i,j,k)              &
                     - zas * ( tke(i,j,k-1,nnow) - tke(i,j,k,nnow) )         &
                     - zcs * ( tke(i,j,k+1,nnow) - tke(i,j,k,nnow) )
          zz         = 1.0_ireals/( zbg - zag*zc(i,j,k-1) )
          zc (i,j,k) = zcg*zz
          zd1(i,j,k) = ( zdg - zag*zd1(i,j,k-1) )*zz
        ENDDO
      ENDDO
    ENDDO

    ! The bottom layer
    DO j = jstart, jend
      DO i = istart, iend
        zag        = 0.5*(a1t(ke)+a1t(ke-1)) * zgatz(i,j,ke)
        zas        = 0.5*(a2t(ke)+a2t(ke-1)) * zgatz(i,j,ke)
        zcg        = 0.5*(a1t(ke1)+a1t(ke))  * zgctz(i,j,ke)
        zcs        = 0.5*(a2t(ke1)+a2t(ke))  * zgctz(i,j,ke)
        zbg        = zzdtr - zag - zcg
        zdg        = zzdtr * tke(i,j,ke,nnew) + tketens(i,j,ke)               &
                     - zas * ( tke(i,j,ke-1,nnow) - tke(i,j,ke,nnow) )        &
                     + zcs * tke(i,j,ke,nnow)
        zz         = 1.0_ireals / ( zbg  - zag*zc(i,j,ke-1) )
        tke(i,j,ke,nnew) = ( zdg -zag*zd1(i,j,ke-1) ) * zz
      ENDDO
    ENDDO

    ! Backsubstitution
    DO k = ke-1, 2, -1
      DO j = jstart, jend
        DO i = istart, iend
          tke(i,j,k,nnew) = zd1(i,j,k) - zc(i,j,k)*tke(i,j,k+1,nnew)
        ENDDO
      ENDDO
    ENDDO

    CALL clipping( tke(:,:,1:ke,nnew), ie, je, ke )

  END IF

  !----------------------------------------------------------------------------
  ! End of subroutine complete_tendencies_qvqcqi_tke 
  !----------------------------------------------------------------------------

END SUBROUTINE complete_tendencies_qvqcqi_tke

!==============================================================================

! CK 20101204 ART RK interfaces
#ifdef COSMOART
!==============================================================================
!+ COSMO_ART procedure for completing the time step for cgas,caero
!------------------------------------------------------------------------------

SUBROUTINE complete_tendencies_cosmo_art

!------------------------------------------------------------------------------
!
! Description:
!   This procedure is only used in src_runge_kutta and calculates the
!   final updates for the prognostic variables cgas and caero
!   at time level n+1 (nnew). 
!   The vertical diffusion as the last slow tendency 
!   is computed here for all variables the vertical diffusion acts on. 
!
! Method:
!   Using the previously calculated updates (on the timelevel nnew) cgas and caero
!   from advection, physics and adiabatic processes, the 
!   vertical diffusion is solved by a vertically implicit scheme (modified 
!   Crank-Nicolson). 
!
!------------------------------------------------------------------------------

! Declarations:

! Local scalars:
! -------------
INTEGER (KIND=iintegers) ::  &
  i,  j,  k, isp           !  Loop indices in lon., lat. and vert. direction
INTEGER (KIND=iintegers) ::  &
  izstata                  !  Loop indices in lon., lat. and vert. direction

INTEGER (KIND=iintegers) :: &
  km1, kp1
  
REAL    (KIND=ireals   ) ::  &
  zbetav, zbetp, zbetm,  & !
  zgat, zgct,            & !
  zag, zas,              & !
  zcg, zcs,              & !
  zbg, zdg,              & !
  zd1g, zd2g,            & !
  zd3g, zd4g,            & !
  znew, zz,              & !
  zd, zzdtr,             & !
  zd1ggas,zd1gaero,      & ! COSMO_ART
  ztvb                !

! Local (automatic) arrays:
! ------------------------
REAL    (KIND=ireals   ) ::  &
  zgavx   (ie,je,ke),    & !
  zgcvx   (ie,je,ke),    & !
  zwcdz   (ie,je,ke),    & !
  zgatz   (ie,je,ke),    & !
  zgctz   (ie,je,ke),    & !
  zqvcor  (ie,je   ),    & !
  zqccor  (ie,je   ),    & !
  zc      (ie,je,ke),    & ! Upper main diagonal of ...
  zd1     (ie,je,ke),    & ! Right hand side of ...
  zd2     (ie,je,ke),    & ! Right hand side of ...
  zd3     (ie,je,ke),    & ! Right hand side of ...
  zd4     (ie,je,ke),    & ! Right hand side of ...
  ze      (ie,je,ke)       ! Soluton vector of ...

REAL (KIND = ireals), ALLOCATABLE :: &
  zd1gas    (:,:,:,:),         & ! Right hand side of ...
  zd1aero   (:,:,:,:),         & !
  zd3gas    (:,:,:,:),         & ! Right hand side of ...
  zd3aero   (:,:,:,:)

! End of header
!==============================================================================

!------------------------------------------------------------------------------
! Begin Subroutine complete_tendencies_cosmo_art
!------------------------------------------------------------------------------

  IF (lgas) THEN
    ALLOCATE (zd1gas(ie,je,ke,isp_gastrans),    &
              zd3gas(ie,je,ke,isp_gastrans),STAT=izstata)
  ENDIF

  IF (laero) THEN
    ALLOCATE (zd1aero(ie,je,ke,isp_aerotrans),    &
              zd3aero(ie,je,ke,isp_aerotrans),STAT=izstata)
  ENDIF

!------------------------------------------------------------------------------
! Section 1: Some preparations
!------------------------------------------------------------------------------

  ! Setting of reciprocal time step
  zzdtr = 1.0 / dt

  IF ( lvertdiff ) THEN
  
    !--------------------------------------------------------------------------
    ! Section 2a: Setup of tridiagonal matrix systems resulting from the
    !             implicit numerical formulation of diffusion of cgas and caero.
    !--------------------------------------------------------------------------

    ! First, the matrix elements a(k), b(k) and c(k) of the coefficient
    ! matrix are set.
    ! The right hand side is stored on d(k).

    ! Top layer
    IF(lgas) THEN
      DO isp = 1,isp_gastrans
        DO j = jstart, jend
          DO i = istart, iend
            zgct       = - zkh(i,j,2)*zsqrtgrho_r_s(i,j,1)
            zcg        = zgct*za1t(2)
            zcs        = zgct*za2t(2)
            zbg        = zzdtr - zcg
            zd3g       = zzdtr * cgas(i,j,1,isp,nnew) + cgast_diff(i,j,1,isp)    &
                       - zcs * ( cgas(i,j,2,isp,nnow) - cgas(i,j,1,isp,nnow) )
            zd1g       = zd3g + cgastens(i,j,1,isp)
            zc (i,j,1) = zcg / zbg
            zd1gas(i,j,1,isp) = zd1g / zbg
            zd3gas(i,j,1,isp) = zd3g / zbg
          ENDDO
        ENDDO
      ENDDO
    ENDIF

    IF(laero) THEN
      DO isp = 1,isp_aerotrans
        DO j = jstart, jend
          DO i = istart, iend
            zgct       = - zkh(i,j,2)*zsqrtgrho_r_s(i,j,1)
            zcg        = zgct*za1t(2)
            zcs        = zgct*za2t(2)
            zbg        = zzdtr - zcg
            zd3g       = zzdtr * caero(i,j,1,isp,nnew) + caerot_diff(i,j,1,isp)    &
                       - zcs * ( caero(i,j,2,isp,nnow) - caero(i,j,1,isp,nnow) )
            zd1g       = zd3g + caerotens(i,j,1,isp)
            zc (i,j,1) = zcg / zbg
            zd1aero(i,j,1,isp) = zd1g / zbg
            zd3aero(i,j,1,isp) = zd3g / zbg
          ENDDO
        ENDDO
      ENDDO
    ENDIF

    ! The layers from k=2 to k=ke-1

    IF(lgas) THEN
      DO isp = 1,isp_gastrans
        DO k = 2, ke-1
          DO j = jstart, jend
            DO i = istart, iend
              zgat       = - zkh(i,j,k  )*zsqrtgrho_r_s(i,j,k)
              zgct       = - zkh(i,j,k+1)*zsqrtgrho_r_s(i,j,k)
              zag        = zgat*za1t(k)
              zas        = zgat*za2t(k)
              zcg        = zgct*za1t(k+1)
              zcs        = zgct*za2t(k+1)
              zbg        = zzdtr - zag - zcg
              zd3g       = zzdtr * cgas(i,j,k,isp,nnew) + cgast_diff(i,j,k,isp)      &
                         - zas * ( cgas(i,j,k-1,isp,nnow) - cgas(i,j,k,isp,nnow) )   &
                         - zcs * ( cgas(i,j,k+1,isp,nnow) - cgas(i,j,k,isp,nnow) )
              zd1g       = zd3g + cgastens(i,j,k,isp)
              zz         = 1.0_ireals / ( zbg - zag*zc(i,j,k-1) )
              zc (i,j,k) = zcg * zz
              zd1gas(i,j,k,isp) = ( zd1g -zag*zd1gas(i,j,k-1,isp) ) * zz
              zd3gas(i,j,k,isp) = ( zd3g -zag*zd3gas(i,j,k-1,isp) ) * zz
            ENDDO
          ENDDO
        ENDDO
      ENDDO
    ENDIF

    IF(laero) THEN
      DO isp = 1,isp_aerotrans
        DO k = 2, ke-1
          DO j = jstart, jend
            DO i = istart, iend
              zgat       = - zkh(i,j,k  )*zsqrtgrho_r_s(i,j,k)
              zgct       = - zkh(i,j,k+1)*zsqrtgrho_r_s(i,j,k)
              zag        = zgat*za1t(k)
              zas        = zgat*za2t(k)
              zcg        = zgct*za1t(k+1)
              zcs        = zgct*za2t(k+1)
              zbg        = zzdtr - zag - zcg
              zd3g       = zzdtr * caero(i,j,k,isp,nnew) + caerot_diff(i,j,k,isp)      &
                         - zas * ( caero(i,j,k-1,isp,nnow) - caero(i,j,k,isp,nnow) )   &
                         - zcs * ( caero(i,j,k+1,isp,nnow) - caero(i,j,k,isp,nnow) )
              zd1g       = zd3g + caerotens(i,j,k,isp)
              zz         = 1.0_ireals / ( zbg - zag*zc(i,j,k-1) )
              zc (i,j,k) = zcg * zz
              zd1aero(i,j,k,isp) = ( zd1g -zag*zd1aero(i,j,k-1,isp) ) * zz
              zd3aero(i,j,k,isp) = ( zd3g -zag*zd3aero(i,j,k-1,isp) ) * zz
            ENDDO
          ENDDO
        ENDDO
      ENDDO
    ENDIF

    ! the bottom layer

    IF (lgas) THEN
      DO isp = 1,isp_gastrans
        DO j = jstart, jend
          DO i = istart, iend
            zgat       = - zkh (i,j,ke)*zsqrtgrho_r_s(i,j,ke)
            zgct       = - ztch(i,j   )*zsqrtgrho_r_s(i,j,ke)
            zag        = zgat*za1t(ke)
            zas        = zgat*za2t(ke)
            zcg        = zgct*za1t_surf
            zcs        = zgct*za2t_surf
            zbg        = zzdtr - zag - zcg
            zd3g       = zzdtr * cgas(i,j,ke,isp,nnew) + cgast_diff(i,j,ke,isp)       &
                       - zas * ( cgas(i,j,ke-1,isp,nnow) - cgas(i,j,ke,isp,nnow) )    &
                       + zcs *  cgas(i,j,ke,isp,nnow) - zgct * cgas_s(i,j,isp,nnow)
!                      + zcs *  cgas(i,j,ke,isp,nnow)
            zd1g       = zd3g + cgastens(i,j,ke,isp)
            zz         = 1.0_ireals / ( zbg  - zag*zc(i,j,ke-1) )
            cgas(i,j,ke,isp,nnew) = ( zd1g - zag*zd1gas(i,j,ke-1,isp) ) * zz
          ENDDO
        ENDDO

        ! Backsubstitution and storage of the complete slow tendencies
        DO k = ke-1, 1, -1
          DO j = jstart, jend
            DO i = istart, iend
              cgas(i,j,k,isp,nnew) = zd1gas(i,j,k,isp) - zc(i,j,k)*cgas(i,j,k+1,isp,nnew)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
    ENDIF

    IF (laero) THEN
      DO isp = 1,isp_aerotrans
        DO j = jstart, jend
          DO i = istart, iend
            zgat       = - zkh (i,j,ke)*zsqrtgrho_r_s(i,j,ke)
            zgct       = - ztch(i,j   )*zsqrtgrho_r_s(i,j,ke)
            zag        = zgat*za1t(ke)
            zas        = zgat*za2t(ke)
            zcg        = zgct*za1t_surf
            zcs        = zgct*za2t_surf
            zbg        = zzdtr - zag - zcg
            zd3g       = zzdtr * caero(i,j,ke,isp,nnew) + caerot_diff(i,j,ke,isp)       &
                       - zas * ( caero(i,j,ke-1,isp,nnow) - caero(i,j,ke,isp,nnow) )    &
                       + zcs *  caero(i,j,ke,isp,nnow) - zgct * caero_s(i,j,isp,nnow)
!                      + zcs *  caero(i,j,ke,isp,nnow)
            zd1g       = zd3g + caerotens(i,j,ke,isp)
            zz         = 1.0_ireals / ( zbg  - zag*zc(i,j,ke-1) )
            caero(i,j,ke,isp,nnew) = ( zd1g - zag*zd1aero(i,j,ke-1,isp) ) * zz
          ENDDO
        ENDDO

        ! Backsubstitution and storage of the complete slow tendencies
        DO k = ke-1, 1, -1
          DO j = jstart, jend
            DO i = istart, iend
              caero(i,j,k,isp,nnew) = zd1aero(i,j,k,isp) - zc(i,j,k)*caero(i,j,k+1,isp,nnew)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
    ENDIF

  ELSE      ! not lvertdiff

    IF (lgas) THEN
      DO isp = 1,isp_gastrans
        DO k = 1, ke
          DO j = jstart, jend
            DO i = istart, iend
              cgas(i,j,k,isp,nnew) = cgas(i,j,k,isp,nnew) + dt * cgast_diff(i,j,k,isp)
              cgas(i,j,k,isp,nnew) = cgas(i,j,k,isp,nnew) + dt * cgastens(i,j,k,isp)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
    ENDIF

    IF (laero) THEN
      DO isp = 1,isp_aerotrans
        DO k = 1, ke
          DO j = jstart, jend
            DO i = istart, iend
              caero(i,j,k,isp,nnew) = caero(i,j,k,isp,nnew) + dt * caerot_diff(i,j,k,isp)
              caero(i,j,k,isp,nnew) = caero(i,j,k,isp,nnew) + dt * caerotens(i,j,k,isp)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
    ENDIF

  ENDIF       ! lvertdiff

  IF (lgas) THEN
    DEALLOCATE (zd1gas, zd3gas)
  ENDIF

  IF (laero) THEN
    DEALLOCATE (zd1aero, zd3aero)
  ENDIF

  !----------------------------------------------------------------------------
  ! End of subroutine complete_tendencies_cosmo_art
  !----------------------------------------------------------------------------

END SUBROUTINE complete_tendencies_cosmo_art
#endif

#ifdef POLLEN
!==============================================================================
!+ COSMO_ART procedure for completing the time step for cpollen
!------------------------------------------------------------------------------

SUBROUTINE complete_tendencies_pollen

!------------------------------------------------------------------------------
!
! Description:
!   This procedure is only used in src_runge_kutta and calculates the
!   final updates for the prognostic variable cpollen
!   at time level n+1 (nnew). 
!   The vertical diffusion as the last slow tendency 
!   is computed here for all variables the vertical diffusion acts on. 
!
! Method:
!   Using the previously calculated updates (on the timelevel nnew) for 
!   cpollen from advection, physics and adiabatic processes, the 
!   vertical diffusion is solved by a vertically implicit scheme (modified 
!   Crank-Nicolson). 
!
!------------------------------------------------------------------------------

! Declarations:

! Local scalars:
! -------------
INTEGER (KIND=iintegers) ::  &
  i,  j,  k, isp           !  Loop indices in lon., lat. and vert. direction
INTEGER (KIND=iintegers) ::  &
  izstata                  !  Loop indices in lon., lat. and vert. direction

INTEGER (KIND=iintegers) :: &
  km1, kp1
  
REAL    (KIND=ireals   ) ::  &
  zbetav, zbetp, zbetm,  & !
  zgat, zgct,            & !
  zag, zas,              & !
  zcg, zcs,              & !
  zbg, zdg,              & !
  zd1g, zd2g,            & !
  zd3g, zd4g,            & !
  znew, zz,              & !
  zd, zzdtr,             & !
  zd1gpolle,ztvb                !

! Local (automatic) arrays:
! ------------------------
REAL    (KIND=ireals   ) ::  &
  zgavx   (ie,je,ke),    & !
  zgcvx   (ie,je,ke),    & !
  zwcdz   (ie,je,ke),    & !
  zgatz   (ie,je,ke),    & !
  zgctz   (ie,je,ke),    & !
  zqvcor  (ie,je   ),    & !
  zqccor  (ie,je   ),    & !
  zc      (ie,je,ke),    & ! Upper main diagonal of ...
  zd1     (ie,je,ke),    & ! Right hand side of ...
  zd2     (ie,je,ke),    & ! Right hand side of ...
  zd3     (ie,je,ke),    & ! Right hand side of ...
  zd4     (ie,je,ke),    & ! Right hand side of ...
  ze      (ie,je,ke)       ! Soluton vector of ...

REAL (KIND = ireals), ALLOCATABLE :: &
  zd1polle  (:,:,:,:),         & ! Right hand side of ...
  zd3polle  (:,:,:,:)            ! Right hand side of ...


! End of header
!==============================================================================

!------------------------------------------------------------------------------
! Begin Subroutine complete_tendencies_pollen
!------------------------------------------------------------------------------

  ALLOCATE (zd1polle(ie,je,ke,isp_pollen),     &
            zd3polle(ie,je,ke,isp_pollen),STAT=izstata)

!------------------------------------------------------------------------------
! Section 1: Some preparations
!------------------------------------------------------------------------------

  ! Setting of reciprocal time step
  zzdtr = 1.0 / dt

  IF ( lvertdiff ) THEN
  
    !--------------------------------------------------------------------------
    ! Section 2a: Setup of tridiagonal matrix systems resulting from the
    !             implicit numerical formulation of diffusion of cpollen
    !--------------------------------------------------------------------------

    ! First, the matrix elements a(k), b(k) and c(k) of the coefficient
    ! matrix are set.
    ! The right hand side is stored on d(k).

    ! Top layer
    IF(l_pollen) THEN
      DO isp = 1,isp_pollen
        DO j = jstart, jend
          DO i = istart, iend
            zgct       = - zkh(i,j,2)*zsqrtgrho_r_s(i,j,1)
            zcg        = zgct*za1t(2)
            zcs        = zgct*za2t(2)
            zbg        = zzdtr - zcg
            zd3g       = zzdtr * cpollen(i,j,1,isp,nnew) + cpollent_diff(i,j,1,isp)    &
                       - zcs * ( cpollen(i,j,2,isp,nnow) - cpollen(i,j,1,isp,nnow) )
            zd1g       = zd3g + cpollentens(i,j,1,isp)
            zc (i,j,1) = zcg / zbg
            zd1polle(i,j,1,isp) = zd1g / zbg
            zd3polle(i,j,1,isp) = zd3g / zbg
          ENDDO
        ENDDO
      ENDDO

      ! The layers from k=2 to k=ke-1
      DO isp = 1,isp_pollen
        DO k = 2, ke-1
          DO j = jstart, jend
            DO i = istart, iend
              zgat       = - zkh(i,j,k  )*zsqrtgrho_r_s(i,j,k)
              zgct       = - zkh(i,j,k+1)*zsqrtgrho_r_s(i,j,k)
              zag        = zgat*za1t(k)
              zas        = zgat*za2t(k)
              zcg        = zgct*za1t(k+1)
              zcs        = zgct*za2t(k+1)
              zbg        = zzdtr - zag - zcg
              zd3g       = zzdtr * cpollen(i,j,k,isp,nnew) + cpollent_diff(i,j,k,isp)      &
                         - zas * ( cpollen(i,j,k-1,isp,nnow) - cpollen(i,j,k,isp,nnow) )   &
                         - zcs * ( cpollen(i,j,k+1,isp,nnow) - cpollen(i,j,k,isp,nnow) )
              zd1g       = zd3g + cpollentens(i,j,k,isp)
              zz         = 1.0_ireals / ( zbg - zag*zc(i,j,k-1) )
              zc (i,j,k) = zcg * zz
              zd1polle(i,j,k,isp) = ( zd1g -zag*zd1polle(i,j,k-1,isp) ) * zz
              zd3polle(i,j,k,isp) = ( zd3g -zag*zd3polle(i,j,k-1,isp) ) * zz
            ENDDO
          ENDDO
        ENDDO
      ENDDO

      ! the bottom layer
      DO isp = 1,isp_pollen
        DO j = jstart, jend
          DO i = istart, iend
            zgat       = - zkh (i,j,ke)*zsqrtgrho_r_s(i,j,ke)
            zgct       = - ztch(i,j   )*zsqrtgrho_r_s(i,j,ke)
            zag        = zgat*za1t(ke)
            zas        = zgat*za2t(ke)
            zcg        = zgct*za1t_surf
            zcs        = zgct*za2t_surf
            zbg        = zzdtr - zag - zcg
            zd3g       = zzdtr * cpollen(i,j,ke,isp,nnew) + cpollent_diff(i,j,ke,isp)       &
                       - zas * ( cpollen(i,j,ke-1,isp,nnow) - cpollen(i,j,ke,isp,nnow) )    &
                       + zcs *  cpollen(i,j,ke,isp,nnow) - zgct * cpollen_s(i,j,isp,nnow)
!                      + zcs *  cpollen(i,j,ke,isp,nnow)
            zd1g       = zd3g + cpollentens(i,j,ke,isp)
            zz         = 1.0_ireals / ( zbg  - zag*zc(i,j,ke-1) )
            cpollen(i,j,ke,isp,nnew) = ( zd1g - zag*zd1polle(i,j,ke-1,isp) ) * zz
          ENDDO
        ENDDO

        ! Backsubstitution and storage of the complete slow tendencies
        DO k = ke-1, 1, -1
          DO j = jstart, jend
            DO i = istart, iend
              cpollen(i,j,k,isp,nnew) = zd1polle(i,j,k,isp) - zc(i,j,k)*cpollen(i,j,k+1,isp,nnew)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
    ENDIF

  ELSE        ! not lvertdiff

    IF(l_pollen) THEN
      DO isp = 1,isp_pollen
        DO k = 1, ke
          DO j = jstart, jend
            DO i = istart, iend
              cpollen(i,j,k,isp,nnew) = cpollen(i,j,k,isp,nnew) + dt * cpollent_diff(i,j,k,isp)
              cpollen(i,j,k,isp,nnew) = cpollen(i,j,k,isp,nnew) + dt * cpollentens(i,j,k,isp)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
    ENDIF
    
  ENDIF

  !----------------------------------------------------------------------------
  ! End of subroutine complete_tendencies_pollen
  !----------------------------------------------------------------------------

  DEALLOCATE (zd1polle, zd3polle,STAT=izstata)

END SUBROUTINE complete_tendencies_pollen
#endif

!MU (17.04.2012) --- BEGIN ---

!==============================================================================
!+ COSMO procedure for completing the time step for tracer
!------------------------------------------------------------------------------

SUBROUTINE complete_tendencies_tracer

!------------------------------------------------------------------------------
!
! Description:
!   This procedure is only used in src_runge_kutta and calculates the
!   final updates for the prognostic variable tracer
!   at time level n+1 (nnew). 
!   The vertical diffusion as the last slow tendency 
!   is computed here for all tracers the vertical diffusion acts on. 
!
! Method:
!   Using the previously calculated updates (on the timelevel nnew) for 
!   tracer from advection, physics and adiabatic processes, the 
!   vertical diffusion is solved by a vertically implicit scheme (modified 
!   Crank-Nicolson). 
!
!------------------------------------------------------------------------------

! Declarations:

! Local scalars:
! -------------
INTEGER (KIND=iintegers) ::  &
  i,  j,  k,                 & ! Loop indices in lon., lat. and vert. direction
  iig                          ! Loop index for tracer

REAL    (KIND=ireals   ) ::  &
  zgat, zgct,            & !
  zag, zas,              & !
  zcg, zcs,              & !
  zbg, zdg,              & !
  zd5g,                  & !
  znew, zz,              & !
  zd, zzdtr                !

REAL    (KIND=ireals   ) :: help1,help2

! Local (automatic) arrays:
! ------------------------
REAL    (KIND=ireals   ) ::  &
  ztracercor(ie,je,ntracer),    & !
  zc        (ie,je,ke),         & ! Upper main diagonal of ...
  zd5       (ie,je,ke)            ! Right hand side of ...

! End of header
!==============================================================================

!------------------------------------------------------------------------------
! Begin Subroutine complete_tendencies_tracer  
!------------------------------------------------------------------------------

  !----------------------------------------------------------------------------
  ! Section 1: Some preparations
  !----------------------------------------------------------------------------

  ! Setting of reciprocal time step
  zzdtr = 1.0 / dt

  !----------------------------------------------------------------------------
  ! Section 2a: Setup of tridiagonal matrix systems resulting from the
  !             implicit numerical formulation of tracer.
  !----------------------------------------------------------------------------

!    For tracer the switch lvertdiff is replaced by ltracer(3,ntracer) to switch
!    turbulence on/off independently from the other prognostic variables.
!    Turbulent transport of tracer (ltracer(3,*) == 1) .AND. lvertdiff == .FALSE.
!    is not allowed.
!    Note: lvertdiff = (ltur .AND. (itype_turb /= 3 .OR. imode_turb < 2))
!    The procedure is the same as for the other prognostic variables:
!
!    First, the matrix elements a(k), b(k) and c(k) of the coefficient
!    matrix are set.
!    The right hand side is stored on d(k).

  DO iig = 1, ntracer

    IF ( ltracer(3,iig) == 1 .AND. ltur .AND. (itype_turb /= 3 .OR. imode_turb < 2)  ) THEN

      ! Top layer
      DO j = jstart, jend
        DO i = istart, iend
          zgct       = - zkh(i,j,2)*zsqrtgrho_r_s(i,j,1)
          zcg        = zgct*za1t(2)
          zcs        = zgct*za2t(2)
          zbg        = zzdtr - zcg
          zd5g       = zzdtr * tracer(i,j,1,nnew,iig) + tracertens(i,j,1,iig)    &
                       - zcs * ( tracer(i,j,2,nnow,iig) - tracer(i,j,1,nnow,iig) )
          zc (i,j,1) = zcg / zbg
          zd5(i,j,1) = zd5g / zbg
        ENDDO
        ! The layers from k=2 to k=ke-1
        DO k = 2, ke-1
          DO i = istart, iend
            zgat       = - zkh(i,j,k  )*zsqrtgrho_r_s(i,j,k)
            zgct       = - zkh(i,j,k+1)*zsqrtgrho_r_s(i,j,k)
            zag        = zgat*za1t(k)
            zas        = zgat*za2t(k)
            zcg        = zgct*za1t(k+1)
            zcs        = zgct*za2t(k+1)
            zbg        = zzdtr - zag - zcg
            zd5g       = zzdtr * tracer(i,j,k,nnew,iig) + tracertens(i,j,k,iig)      &
                         - zas * ( tracer(i,j,k-1,nnow,iig) - tracer(i,j,k,nnow,iig) )   &
                         - zcs * ( tracer(i,j,k+1,nnow,iig) - tracer(i,j,k,nnow,iig) )
            zz         = 1.0_ireals / ( zbg - zag*zc(i,j,k-1) )
            zc (i,j,k) = zcg * zz
            zd5(i,j,k) = ( zd5g -zag*zd5(i,j,k-1) ) * zz
          ENDDO
        ENDDO
        ! the bottom layer
        DO i = istart, iend
          zgat       = - zkh (i,j,ke)*zsqrtgrho_r_s(i,j,ke)
          zgct       = - ztch(i,j   )*zsqrtgrho_r_s(i,j,ke)
          zag        = zgat*za1t(ke)
          zas        = zgat*za2t(ke)
          zcg        = zgct*za1t_surf
          zcs        = zgct*za2t_surf
          zbg        = zzdtr - zag - zcg
          zd5g       = zzdtr * tracer(i,j,ke,nnew,iig) + tracertens(i,j,ke,iig)       &
                       - zas * ( tracer(i,j,ke-1,nnow,iig) - tracer(i,j,ke,nnow,iig) )    &
                       + zcs * zfac_qc * tracer(i,j,ke,nnow,iig)
          zz         = 1.0_ireals / ( zbg  - zag*zc(i,j,ke-1) )
          tracer(i,j,ke,nnew,iig) = ( zd5g - zag*zd5(i,j,ke-1) ) * zz
        ENDDO
        ! Backsubstitution and storage of the complete slow tendencies
        DO k = ke-1, 1, -1
          DO i = istart, iend
            tracer(i,j,k,nnew,iig) = zd5(i,j,k) - zc(i,j,k)*tracer(i,j,k+1,nnew,iig)
          ENDDO
        ENDDO
      ENDDO

    ELSE

      DO k = 1, ke
        DO j = jstart, jend
          DO i = istart, iend
            tracer(i,j,k,nnew,iig) = tracer(i,j,k,nnew,iig) + dt * tracertens(i,j,k,iig)
          ENDDO
        ENDDO
      ENDDO

    END IF

  ENDDO

  DO iig = 1, ntracer
    IF (.NOT. lvertdiff .AND. (my_cart_id==1) .AND. (ltracer(3,iig)==1) .AND. (ntstep<2) ) THEN
      PRINT *, 'WARNING: No turbulent tracer transport although ltracer(3,*) == 1!'
    ENDIF
  ENDDO

  !----------------------------------------------------------------------------
  ! Section 4: Massflux correction scheme for removing negative tracer values
  !            (if required)
  !----------------------------------------------------------------------------

  DO iig = 1, ntracer

    IF ( ltracer(6,iig) == 1 ) THEN

      IF ( lmassf_diffusion ) THEN

        DO j = jstart, jend
          DO i = istart, iend
            ztracercor(i,j,iig) = 0.0_ireals
          ENDDO
        ENDDO
        DO k = 1, ke
          DO j = jstart, jend
            DO i = istart, iend
              help1 = tracer(i,j,k,nnew,iig) + ztracercor(i,j,iig)*zsqrtgrho_r_s(i,j,k)
              IF ( help1 < 0.0_ireals ) THEN
                help2 = 0.0_ireals
                help1 = help1/zsqrtgrho_r_s(i,j,k)
              ELSE
                help2 = help1
                help1 = 0.0_ireals
              ENDIF
              ztracercor(i,j,iig)    = help1
              tracer(i,j,k,nnew,iig) = help2
            ENDDO
          ENDDO
        ENDDO

      ELSE

        CALL clipping( tracer(:,:,:,nnew,iig), ie, je, ke)
    
      ENDIF ! end mass-flux correction

    ENDIF

  ENDDO

!------------------------------------------------------------------------------
! End of subroutine complete_tendencies_tracer
!------------------------------------------------------------------------------

END SUBROUTINE complete_tendencies_tracer

!MU (17.04.2012) --- END ---

!==============================================================================

!option! -pvctl _on_adb
SUBROUTINE explicit_horizontal_diffusion

! USE data_runcontrol, ONLY: l3dturb_metr

  IMPLICIT NONE

  ! Local scalars:
  !---------------
  REAL    (KIND=ireals   ) ::  &
    zhfdx, zhfdy, zhfkh,    & !
    zcrlatr, zarhor

  INTEGER (KIND=iintegers) :: &
    i, j, k, klow, kup, &
    iig !AK

  LOGICAL                  :: &
    lkge2_prog_tke

  ! Local arrays:
  !--------------
  REAL    (KIND=ireals   ) ::  &
    zrhokvm (ie,je   ),     & !
    zrhokhm (ie,je   ),     & !
    zrhokhh (ie,je   ),     & !
    ztau11  (ie,je,ke),     & !
    ztau12  (ie,je,ke),     & !
    ztau13  (ie,je,ke),     & !
    ztau22  (ie,je,ke),     & !
    ztau23  (ie,je,ke),     & !
    ztaud13 (ie,je,ke),     & !
    ztaud23 (ie,je,ke),     & !
    zth1    (ie,je,ke),     & !
    zth2    (ie,je,ke),     & !
    zqvh1   (ie,je,ke),     & !
    zqch1   (ie,je,ke),     & !
    zqih1   (ie,je,ke),     & !
    zqvh2   (ie,je,ke),     & !
    zqch2   (ie,je,ke),     & !
    zqih2   (ie,je,ke),     & !
!AK (30.03.2012)
!    ztracerh1 (ie,je,ke,ntracer), & ! local variable to handle with tracer
!    ztracerh2 (ie,je,ke,ntracer), & ! local variable to handle with tracer
!AK (30.03.2012)
    ztkeh1  (ie,je,ke),     & !
    ztkeh2  (ie,je,ke)

  ! End of header
  !==============================================================================

  ! -------- (1) calc. of the fluxes -------------------------------------
  DO k = 1, ke

    kup  = MAX( 2, k )
    klow = MIN( ke, k+1 )

    IF ( k >= 2 .AND. lprog_tke ) THEN
      lkge2_prog_tke = .TRUE.
    ELSE
      lkge2_prog_tke = .FALSE.
    END IF

! UB>>
!!$    DO j = jstart-1, jend+1
!!$      DO  i = istart-1, iend+1
    DO j = jstart-2, jend+2
      DO  i = istart-2, iend+2
! UB<<

        zrhokhm(i,j) = rho(i,j,k) * 0.5_ireals*(tkhm(i,j,kup)+tkhm(i,j,klow))
        zrhokhh(i,j) = rho(i,j,k) * 0.5_ireals*(tkhh(i,j,kup)+tkhh(i,j,klow))
        zrhokvm(i,j) = tkvm(i,j,kup) * 0.5_ireals*(rho(i,j,kup-1)+rho(i,j,k))

      ENDDO
    ENDDO

! UB>>
    IF ( lprog_tke ) THEN
      ztkeh1(:,:,1) = 0.0_ireals
      ztkeh2(:,:,1) = 0.0_ireals
    END IF
! UB<<

    DO j = jstart, jendv+1
      DO  i = istart, iendu+1
        ztau11(i,j,k) = -2.0 * zrhokhm(i,j) * acrlat(j,1) *                   &
                      ( u(i,j,k,nnow) - u(i-1,j,k,nnow) ) * eddlon
        ztau22(i,j,k) = -2.0 * zrhokhm(i,j) / r_earth *                       &
                      ( v(i,j,k,nnow) - v(i,j-1,k,nnow) ) * eddlat
      ENDDO
    ENDDO

! UB>>
!!$    DO j = jstart-1, jend
!!$      DO i = istart-1, iend
    DO j = jstart-1, jend+1
      DO i = istart-1, iend+1
! UB<<
        ztau12(i,j,k) = -0.25_ireals *                                        &
                           ( zrhokhm(i,j  )  + zrhokhm(i+1,j  )               &
                           + zrhokhm(i,j+1)  + zrhokhm(i+1,j+1) )             &
                         * acrlat(j,2) *                                      &
                           ( ( v(i+1,j,k,nnow) - v(i,j,k,nnow) ) * eddlon     &
              + crlat(j,2) * ( u(i,j+1,k,nnow) - u(i,j,k,nnow) ) * eddlat )
      ENDDO
    ENDDO

    DO j = jstart, jend
      DO i = istart-1, iend

        ! part of T13 which is explicitly handled in this subr.:
        ztaud13(i,j,k) = -0.5_ireals * (zrhokvm(i,j)+zrhokvm(i+1,j))          &
                          * acrlat(j,1) *                                     &
                          ( w(i+1,j,k,nnow) - w(i,j,k,nnow) ) * eddlon
        ! part of T13 which is implicitly handled for u-eq:
        ztau13(i,j,k)  =  ztaud13(i,j,k) +                                    &
                          0.5_ireals * (ztmkvm(i,j,kup)+ztmkvm(i+1,j,kup)) *  &
                                  ( u(i,j,k,nnow) - u(i,j,kup-1,nnow) )

        IF ( lmoist_turb ) THEN
          zhfdx        = ( ztheta_l(i+1,j,k) - ztheta_l(i,j,k) ) * eddlon
        ELSE
          zhfdx        = ( ztheta(i+1,j,k) - ztheta(i,j,k) ) * eddlon
        END IF
        zhfkh        = 0.5_ireals * ( zpia(i  ,j,k) * zrhokhh(i  ,j)          &
                                    + zpia(i+1,j,k) * zrhokhh(i+1,j) )        &
                                    * acrlat(j,1)
        zth1(i,j,k)  = - zhfkh * zhfdx

        zhfkh        = 0.5_ireals * ( zrhokhh(i,j) + zrhokhh(i+1,j) )         &
                             * acrlat(j,1)

        zhfdx        = ( qv(i+1,j,k,nnow) - qv(i,j,k,nnow) ) * eddlon
        zqvh1(i,j,k) = - zhfkh * zhfdx

        zhfdx        = ( qc(i+1,j,k,nnow) - qc(i,j,k,nnow) ) * eddlon
        zqch1(i,j,k) = - zhfkh * zhfdx
!AK (30.03.2012)
!        DO iig = 1, ntracer
!          zhfdx = ( tracer(i+1,j,k,nnow,iig) - tracer(i,j,k,nnow,iig) ) * eddlon
!          ztracerh1(i,j,k,iig) = - zhfkh * zhfdx
!        ENDDO
!AK (30.03.2012)

        IF ( lprog_qi ) THEN
          zhfdx        = ( qi(i+1,j,k,nnow) - qi(i,j,k,nnow) ) * eddlon
          zqih1(i,j,k) = - zhfkh * zhfdx
        END IF

        IF ( lkge2_prog_tke ) THEN
          zhfkh         = ( tkhm(i,j,k) + tkhm(i+1,j,k) ) * acrlat(j,1)
          zhfdx         = ( tke(i+1,j,k,nnow) - tke(i,j,k,nnow) ) * eddlon
          ztkeh1(i,j,k) = - zhfkh * zhfdx
        END IF

      ENDDO
    ENDDO

    DO j = jstart-1, jend
! UB>>
!!$      DO i = istart, iend
      DO i = istart-1, iend
! UB<<

        ! part of T23 which is explicitly handled:
        ztaud23(i,j,k) = -0.5_ireals * (zrhokvm(i,j)+ zrhokvm(i,j+1))          &
                          / r_earth *                                          &
                              ( w(i,j+1,k,nnow) - w(i,j,k,nnow) ) * eddlat
        ! part of T23 which is implicitly handled for v-eq:
        ztau23(i,j,k)  =  ztaud23(i,j,k) +                                     &
                          0.5_ireals * (ztmkvm(i,j,kup)+ztmkvm(i,j+1,kup)) *   &
                               ( v(i,j,k,nnow) - v(i,j,kup-1,nnow) )

        IF ( lmoist_turb ) THEN
          zhfdy        = ( ztheta_l(i,j+1,k) - ztheta_l(i,j,k) ) * eddlat
        ELSE
          zhfdy        = ( ztheta(i,j+1,k) - ztheta(i,j,k) ) * eddlat
        END IF
        zhfkh        = 0.5_ireals * ( zpia(i,j  ,k) * zrhokhh(i,j  )           &
                                    + zpia(i,j+1,k) * zrhokhh(i,j+1) )         &
                                    / r_earth
        zth2(i,j,k)  = - zhfkh * zhfdy

        zhfkh        = 0.5_ireals * ( zrhokhh(i,j) + zrhokhh(i,j+1) )          &
                              / r_earth

        zhfdy        = ( qv(i,j+1,k,nnow) - qv(i,j,k,nnow) ) * eddlat
        zqvh2(i,j,k) = - zhfkh * zhfdy

        zhfdy        = ( qc(i,j+1,k,nnow) - qc(i,j,k,nnow) ) * eddlat
        zqch2(i,j,k) = - zhfkh * zhfdy

!AK (30.03.2012)
!        DO iig = 1, ntracer
!          zhfdy      = ( tracer(i,j+1,k,nnow,iig) - tracer(i,j,k,nnow,iig) ) * eddlat
!          ztracerh2(i,j,k,iig) = - zhfkh * zhfdy
!        ENDDO
!AK (30.03.2012)

        IF ( lprog_qi ) THEN
          zhfdy        = ( qi(i,j+1,k,nnow) - qi(i,j,k,nnow) ) * eddlat
          zqih2(i,j,k) = - zhfkh * zhfdy
        END IF

        IF ( lkge2_prog_tke ) THEN
          zhfkh         = ( tkhm(i,j,k) + tkhm(i,j+1,k) ) / r_earth
          zhfdy         = ( tke(i,j+1,k,nnow) - tke(i,j,k,nnow) ) * eddlat
          ztkeh2(i,j,k) = - zhfkh * zhfdy
        END IF

      ENDDO
    ENDDO


    ! metric terms for 3D turbulence
    IF ( l3dturb_metr ) THEN

      IF ( ( k>=2 ) .AND. ( k<=ke-1 ) ) THEN

        DO j = jstart, jend
          DO  i = istart-1, iend

            ztaud13(i,j,k) = ztaud13(i,j,k)                                   &
                    -0.5_ireals * (zrhokvm(i,j)+zrhokvm(i+1,j))               &
                    * acrlat(j,1)                                             &
                    * 0.25 * ( dzeta_dlam(i,j,k  ) + dzeta_dlam(i+1,j,k  )    &
                             + dzeta_dlam(i,j,k-1) + dzeta_dlam(i+1,j,k-1) )  &
                    * 0.25 * ( ( w(i  ,j,k+1,nnow) - w(i  ,j,k-1,nnow) )      &
                             + ( w(i+1,j,k+1,nnow) - w(i+1,j,k-1,nnow) ) )

            ztaud23(i,j,k) = ztaud23(i,j,k)                                   &
                    -0.5_ireals * (zrhokvm(i,j)+ zrhokvm(i,j+1))              &
                        / r_earth                                             &
                    * 0.25 * ( dzeta_dlam(i,j,k  ) + dzeta_dlam(i,j+1,k  )    &
                             + dzeta_dlam(i,j,k-1) + dzeta_dlam(i,j+1,k-1) )  &
                    * 0.25 * ( ( w(i,j  ,k+1,nnow) - w(i,j  ,k-1,nnow) )      &
                             + ( w(i,j+1,k+1,nnow) - w(i,j+1,k-1,nnow) ) )

          ENDDO
        ENDDO

      ELSE IF ( k==1 ) THEN

        DO j = jstart, jend
          DO  i = istart-1, iend

            ztaud13(i,j,k) = ztaud13(i,j,k)                                   &
                    - 0.5_ireals * (zrhokvm(i,j)+zrhokvm(i+1,j))              &
                    * acrlat(j,1)                                             &
                    * 0.5  * ( dzeta_dlam(i,j,k  ) + dzeta_dlam(i+1,j,k  ) )  &
                    * 0.5  * ( ( w(i  ,j,k+1,nnow) - w(i  ,j,k,nnow) )        &
                             + ( w(i+1,j,k+1,nnow) - w(i+1,j,k,nnow) ) )

            ztaud23(i,j,k) = ztaud23(i,j,k)                                   &
                    - 0.5_ireals * (zrhokvm(i,j)+ zrhokvm(i,j+1))             &
                    / r_earth                                                 &
                    * 0.5  * ( dzeta_dlam(i,j,k  ) + dzeta_dlam(i,j+1,k  ) )  &
                    * 0.5  * ( ( w(i,j  ,k+1,nnow) - w(i,j  ,k,nnow) )        &
                             + ( w(i,j+1,k+1,nnow) - w(i,j+1,k,nnow) ) )

          ENDDO
        ENDDO

      ELSE IF ( k==ke ) THEN

        DO j = jstart, jend
          DO  i = istart-1, iend

            ztaud13(i,j,k) = ztaud13(i,j,k)                                   &
                    - 0.5_ireals * (zrhokvm(i,j)+zrhokvm(i+1,j))              &
                    * acrlat(j,1)                                             &
                    * 0.25 * ( dzeta_dlam(i,j,k  ) + dzeta_dlam(i+1,j,k  )    &
                             + dzeta_dlam(i,j,k-1) + dzeta_dlam(i+1,j,k-1) )  &
                    * 0.5  * ( ( w(i  ,j,k,nnow) - w(i  ,j,k-1,nnow) )        &
                             + ( w(i+1,j,k,nnow) - w(i+1,j,k-1,nnow) ) )

            ztaud23(i,j,k) = ztaud23(i,j,k)                                   &
                    - 0.5_ireals * (zrhokvm(i,j)+ zrhokvm(i,j+1))             &
                    / r_earth                                                 &
                    * 0.25 * ( dzeta_dlam(i,j,k  ) + dzeta_dlam(i,j+1,k  )    &
                             + dzeta_dlam(i,j,k-1) + dzeta_dlam(i,j+1,k-1) )  &
                    * 0.5  * ( ( w(i,j  ,k,nnow) - w(i,j  ,k-1,nnow) )        &
                             + ( w(i,j+1,k,nnow) - w(i,j+1,k-1,nnow) ) )

          ENDDO
        ENDDO

      END IF

      IF ( ( k>=2 ) .AND. ( k<=ke-1 ) ) THEN

        DO j = jstart, jendv+1
          DO  i = istart, iendu+1
            ztau11(i,j,k) = ztau11(i,j,k)                                      &
                    - 2.0 * zrhokhm(i,j) * acrlat(j,1)                         &
                     *  dzeta_dlam(i,j,k)                                      &
                     * 0.25 * ( ( u(i-1,j,k+1,nnow) - u(i-1,j,k-1,nnow) )      &
                              + ( u(i  ,j,k+1,nnow) - u(i  ,j,k-1,nnow) ) )

            ztau12(i,j,k) =  ztau12(i,j,k)                                     &
                    -0.25_ireals *                                             &
                         ( zrhokhm(i,j  )  + zrhokhm(i+1,j  )                  &
                         + zrhokhm(i,j+1)  + zrhokhm(i+1,j+1) )                &
                     * acrlat(j,2)                                             &
                     * 0.25 * ( dzeta_dlam(i,j  ,k) + dzeta_dlam(i+1,j  ,k)    &
                              + dzeta_dlam(i,j+1,k) + dzeta_dlam(i+1,j+1,k) )  &
                     * 0.25 * ( ( v(i  ,j,k+1,nnow) - v(i  ,j,k-1,nnow) )      &
                              + ( v(i+1,j,k+1,nnow) - v(i+1,j,k-1,nnow) ) )

            ztau22(i,j,k) = ztau22(i,j,k)                                      &
                    -2.0 * zrhokhm(i,j) / r_earth                              &
                     *  dzeta_dphi(i,j,k)                                      &
                     *  0.25 * ( ( v(i,j  ,k+1,nnow) - v(i,j  ,k-1,nnow) )     &
                               + ( v(i,j-1,k+1,nnow) - v(i,j-1,k-1,nnow) ) )

          ENDDO
        ENDDO

        ! scalar variables, 1. flux component
        DO j = jstart, jend
          DO i = istart-1, iend

            IF ( lmoist_turb ) THEN
              zhfdx  = 0.25 * ( ( ztheta_l(i  ,j,k+1) - ztheta_l(i  ,j,k-1) )  &
                              + ( ztheta_l(i+1,j,k+1) - ztheta_l(i+1,j,k-1) ) )
            ELSE
              zhfdx  = 0.25 * ( ( ztheta  (i  ,j,k+1) - ztheta  (i  ,j,k-1) )  &
                              + ( ztheta  (i+1,j,k+1) - ztheta  (i+1,j,k-1) ) )
            END IF
            zhfkh        = 0.5_ireals * ( zpia(i  ,j,k) * zrhokhh(i  ,j)       &
                                        + zpia(i+1,j,k) * zrhokhh(i+1,j) )     &
                                      * acrlat(j,1)                            &
                            * 0.5 * ( dzeta_dlam(i,j,k) + dzeta_dlam(i+1,j,k) )
            zth1(i,j,k)  = zth1(i,j,k) - zhfkh * zhfdx

            zhfkh        = 0.5_ireals * ( zrhokhh(i,j) + zrhokhh(i+1,j) )      &
                         * acrlat(j,1)                                         &
                         * 0.5 * ( dzeta_dlam(i,j,k) + dzeta_dlam(i+1,j,k) )

            zhfdx = 0.25 * ( ( qv(i  ,j,k+1,nnow) - qv(i  ,j,k-1,nnow) )       &
                           + ( qv(i+1,j,k+1,nnow) - qv(i+1,j,k-1,nnow) ) )
            zqvh1(i,j,k) = zqvh1(i,j,k) - zhfkh * zhfdx

            zhfdx = 0.25 * ( ( qc(i  ,j,k+1,nnow) - qc(i  ,j,k-1,nnow) )       &
                           + ( qc(i+1,j,k+1,nnow) - qc(i+1,j,k-1,nnow) ) )
            zqch1(i,j,k) = zqch1(i,j,k) - zhfkh * zhfdx

!AK (30.03.2012)
!            DO iig = 1, ntracer
!              zhfdx = 0.25 * ( ( tracer(i  ,j,k+1,nnow,iig) - tracer(i  ,j,k-1,nnow,iig) )       &
!                             + ( tracer(i+1,j,k+1,nnow,iig) - tracer(i+1,j,k-1,nnow,iig) ) )
!              ztracerh1(i,j,k,iig) = ztracerh1(i,j,k,iig) - zhfkh * zhfdx
!            ENDDO
!AK (30.03.2012)

            IF ( lprog_qi ) THEN
              zhfdx = 0.25 * ( ( qi(i  ,j,k+1,nnow) - qi(i  ,j,k-1,nnow) )     &
                             + ( qi(i+1,j,k+1,nnow) - qi(i+1,j,k-1,nnow) ) )
              zqih1(i,j,k) = zqih1(i,j,k) - zhfkh * zhfdx
            END IF

! UB>>
            IF ( lkge2_prog_tke ) THEN
!!$            IF ( lprog_tke ) THEN
              zhfkh = ( tkhm(i,j,k) + tkhm(i+1,j,k) ) * acrlat(j,1)
              zhfdx = 0.25 * ( ( tke(i  ,j,k+1,nnow) - tke(i  ,j,k-1,nnow) )   &
                             + ( tke(i+1,j,k+1,nnow) - tke(i+1,j,k-1,nnow) ) )
              ztkeh1(i,j,k) = ztkeh1(i,j,k) - zhfkh * zhfdx
            END IF

          ENDDO
        ENDDO

        ! scalar variables, 2. flux component
!UB>>
!!$        DO j = jstart, jend
!!$          DO i = istart-1, iend
        DO j = jstart-1, jend
          DO i = istart, iend
!UB<<

            IF ( lmoist_turb ) THEN
              zhfdy = 0.25 * ( ( ztheta_l(i,j  ,k+1) - ztheta_l(i,j  ,k-1) )   &
                             + ( ztheta_l(i,j+1,k+1) - ztheta_l(i,j+1,k-1) ) )
            ELSE
              zhfdy = 0.25 * ( ( ztheta  (i,j  ,k+1) - ztheta  (i,j  ,k-1) )   &
                             + ( ztheta  (i,j+1,k+1) - ztheta  (i,j+1,k-1) ) )
            END IF
            zhfkh        = 0.5_ireals * ( zpia(i,j  ,k) * zrhokhh(i,j  )       &
                                        + zpia(i,j+1,k) * zrhokhh(i,j+1) )     &
                                      / r_earth                                &
                            * 0.5 * ( dzeta_dphi(i,j,k) + dzeta_dphi(i+1,j,k) )
            zth2(i,j,k)  = zth2(i,j,k) - zhfkh * zhfdy

            zhfkh        = 0.5_ireals * ( zrhokhh(i,j) + zrhokhh(i,j+1) )      &
                         / r_earth                                             &
                         * 0.5 * ( dzeta_dphi(i,j,k) + dzeta_dphi(i,j+1,k) )

            zhfdy = 0.25 * ( ( qv(i,j  ,k+1,nnow) - qv(i,j  ,k-1,nnow) )       &
                           + ( qv(i,j+1,k+1,nnow) - qv(i,j+1,k-1,nnow) ) )
            zqvh2(i,j,k) = zqvh2(i,j,k) - zhfkh * zhfdy

            zhfdy = 0.25 * ( ( qc(i,j  ,k+1,nnow) - qc(i,j  ,k-1,nnow) )       &
                           + ( qc(i,j+1,k+1,nnow) - qc(i,j+1,k-1,nnow) ) )
            zqch2(i,j,k) = zqch2(i,j,k) - zhfkh * zhfdy

!AK (30.03.2012)
!            DO iig = 1, ntracer
!              zhfdy = 0.25 * ( ( tracer(i,j  ,k+1,nnow,iig) - tracer(i,j  ,k-1,nnow,iig) )       &
!                             + ( tracer(i,j+1,k+1,nnow,iig) - tracer(i,j+1,k-1,nnow,iig) ) )
!              ztracerh2(i,j,k,iig) = ztracerh2(i,j,k,iig) - zhfkh * zhfdy
!            ENDDO
!AK (30.03.2012)

            IF ( lprog_qi ) THEN
              zhfdy = 0.25 * ( ( qi(i,j  ,k+1,nnow) - qi(i,j  ,k-1,nnow) )     &
                             + ( qi(i,j+1,k+1,nnow) - qi(i,j+1,k-1,nnow) ) )
              zqih2(i,j,k) = zqih2(i,j,k) - zhfkh * zhfdy
            END IF

! UB>>            
            IF ( lkge2_prog_tke ) THEN
!!$            IF ( lprog_tke ) THEN
              zhfkh         = ( tkhm(i,j,k) + tkhm(i+1,j,k) ) * acrlat(j,1)
              zhfdy = 0.25 * ( ( tke(i,j  ,k+1,nnow) - tke(i,j  ,k-1,nnow) )   &
                             + ( tke(i,j+1,k+1,nnow) - tke(i,j+1,k-1,nnow) ) )
              ztkeh2(i,j,k) = ztkeh2(i,j,k) - zhfkh * zhfdy
            END IF

          ENDDO
        ENDDO

      ELSE IF ( k==1 ) THEN

        DO j = jstart, jendv+1
          DO  i = istart, iendu+1
            ztau11(i,j,k) = ztau11(i,j,k)                                      &
                    - 2.0 * zrhokhm(i,j) * acrlat(j,1)                         &
                     *  dzeta_dlam(i,j,k)                                      &
                     * 0.5 * ( ( u(i-1,j,k+1,nnow) - u(i-1,j,k,nnow) )         &
                             + ( u(i  ,j,k+1,nnow) - u(i  ,j,k,nnow) ) )

            ztau12(i,j,k) =  ztau12(i,j,k)                                     &
                    -0.25_ireals *                                             &
                         ( zrhokhm(i,j  )  + zrhokhm(i+1,j  )                  &
                         + zrhokhm(i,j+1)  + zrhokhm(i+1,j+1) )                &
                     * acrlat(j,2)                                             &
                     * 0.25 * ( dzeta_dlam(i,j  ,k) + dzeta_dlam(i+1,j  ,k)    &
                              + dzeta_dlam(i,j+1,k) + dzeta_dlam(i+1,j+1,k) )  &
                     * 0.5  * ( ( v(i  ,j,k+1,nnow) - v(i  ,j,k,nnow) )        &
                              + ( v(i+1,j,k+1,nnow) - v(i+1,j,k,nnow) ) )

            ztau22(i,j,k) = ztau22(i,j,k)                                      &
                    -2.0 * zrhokhm(i,j) / r_earth                              &
                     *  dzeta_dphi(i,j,k)                                      &
                     *  0.5 * ( ( v(i,j  ,k+1,nnow) - v(i,j  ,k,nnow) )        &
                              + ( v(i,j-1,k+1,nnow) - v(i,j-1,k,nnow) ) )

          ENDDO
        ENDDO

        ! scalar variables, 1. flux component
        DO j = jstart, jend
          DO i = istart-1, iend

            IF ( lmoist_turb ) THEN
              zhfdx  = 0.5 * ( ( ztheta_l(i  ,j,k+1) - ztheta_l(i  ,j,k) )     &
                             + ( ztheta_l(i+1,j,k+1) - ztheta_l(i+1,j,k) ) )
            ELSE
              zhfdx  = 0.5 * ( ( ztheta  (i  ,j,k+1) - ztheta  (i  ,j,k) )     &
                             + ( ztheta  (i+1,j,k+1) - ztheta  (i+1,j,k) ) )
            END IF
            zhfkh = 0.5_ireals * ( zpia(i  ,j,k) * zrhokhh(i  ,j)              &
                                 + zpia(i+1,j,k) * zrhokhh(i+1,j) )            &
                                 * acrlat(j,1)                                 &
                         * 0.5 * ( dzeta_dlam(i,j,k) + dzeta_dlam(i+1,j,k) )
            zth1(i,j,k)  = zth1(i,j,k) - zhfkh * zhfdx

            zhfkh = 0.5_ireals * ( zrhokhh(i,j) + zrhokhh(i+1,j) )             &
                                 * acrlat(j,1)                                 &
                         * 0.5 * ( dzeta_dlam(i,j,k) + dzeta_dlam(i+1,j,k) )

            zhfdx        = 0.5 * ( ( qv(i  ,j,k+1,nnow) - qv(i  ,j,k,nnow) )   &
                                 + ( qv(i+1,j,k+1,nnow) - qv(i+1,j,k,nnow) ) )
            zqvh1(i,j,k) = zqvh1(i,j,k) - zhfkh * zhfdx

            zhfdx        = 0.5 * ( ( qc(i  ,j,k+1,nnow) - qc(i  ,j,k,nnow) )   &
                                 + ( qc(i+1,j,k+1,nnow) - qc(i+1,j,k,nnow) ) )
            zqch1(i,j,k) = zqch1(i,j,k) - zhfkh * zhfdx

!AK (30.03.2012)
!            DO iig = 1, ntracer
!              zhfdx        = 0.5 * ( ( tracer(i  ,j,k+1,nnow,iig) - tracer(i  ,j,k,nnow,iig) )   &
!                                   + ( tracer(i+1,j,k+1,nnow,iig) - tracer(i+1,j,k,nnow,iig) ) )
!              ztracerh1(i,j,k,iig) = ztracerh1(i,j,k,iig) - zhfkh * zhfdx
!            ENDDO
!AK (30.03.2012)

            IF ( lprog_qi ) THEN
              zhfdx = 0.5 * ( ( qi(i  ,j,k+1,nnow) - qi(i  ,j,k,nnow) )        &
                            + ( qi(i+1,j,k+1,nnow) - qi(i+1,j,k,nnow) ) )
              zqih1(i,j,k) = zqih1(i,j,k) - zhfkh * zhfdx
            END IF

! UB>> unnecessary ?
!!$            IF ( lkge2_prog_tke ) THEN
!!$              zhfkh         = ( tkhm(i,j,k) + tkhm(i+1,j,k) ) * acrlat(j,1)
!!$              zhfdx = 0.5 * ( ( tke(i  ,j,k+1,nnow) - tke(i  ,j,k,nnow) )      &
!!$                            + ( tke(i+1,j,k+1,nnow) - tke(i+1,j,k,nnow) ) )
!!$              ztkeh1(i,j,k) = ztkeh1(i,j,k) - zhfkh * zhfdx
!!$            END IF
! UB<<

          ENDDO
        ENDDO

        ! scalar variables, 2. flux component
! UB>>
!!$        DO j = jstart, jend
!!$          DO i = istart-1, iend
        DO j = jstart-1, jend
          DO i = istart, iend
! UB<<

            IF ( lmoist_turb ) THEN
              zhfdy = 0.5 * ( ( ztheta_l(i,j  ,k+1) - ztheta_l(i,j  ,k) )      &
                            + ( ztheta_l(i,j+1,k+1) - ztheta_l(i,j+1,k) ) )
            ELSE
              zhfdy = 0.5 * ( ( ztheta  (i,j  ,k+1) - ztheta  (i,j  ,k) )      &
                            + ( ztheta  (i,j+1,k+1) - ztheta  (i,j+1,k) ) )
            END IF
            zhfkh = 0.5_ireals * ( zpia(i,j  ,k) * zrhokhh(i,j  )              &
                                 + zpia(i,j+1,k) * zrhokhh(i,j+1) )            &
                                 / r_earth                                     &
                         * 0.5 * ( dzeta_dphi(i,j,k) + dzeta_dphi(i+1,j,k) )
            zth2(i,j,k)  = zth2(i,j,k) - zhfkh * zhfdy

            zhfkh        = 0.5_ireals * ( zrhokhh(i,j) + zrhokhh(i,j+1) )      &
                         / r_earth                                             &
                         * 0.5 * ( dzeta_dphi(i,j,k) + dzeta_dphi(i,j+1,k) )

            zhfdy        = 0.5 * ( ( qv(i,j  ,k+1,nnow) - qv(i,j  ,k,nnow) )   &
                                 + ( qv(i,j+1,k+1,nnow) - qv(i,j+1,k,nnow) ) )
            zqvh2(i,j,k) = zqvh2(i,j,k) - zhfkh * zhfdy

            zhfdy        = 0.5 * ( ( qc(i,j  ,k+1,nnow) - qc(i,j  ,k,nnow) )   &
                                 + ( qc(i,j+1,k+1,nnow) - qc(i,j+1,k,nnow) ) )
            zqch2(i,j,k) = zqch2(i,j,k) - zhfkh * zhfdy

!AK (30.03.2012)
!            DO iig = 1, ntracer
!              zhfdy        = 0.5 * ( ( tracer(i,j  ,k+1,nnow,iig) - tracer(i,j  ,k,nnow,iig) )   &
!                                   + ( tracer(i,j+1,k+1,nnow,iig) - tracer(i,j+1,k,nnow,iig) ) )
!              ztracerh2(i,j,k,iig) = ztracerh2(i,j,k,iig) - zhfkh * zhfdy
!            ENDDO
!AK (30.03.2012)


            IF ( lprog_qi ) THEN
              zhfdy = 0.5 * ( ( qi(i,j  ,k+1,nnow) - qi(i,j  ,k,nnow) )        &
                            + ( qi(i,j+1,k+1,nnow) - qi(i,j+1,k,nnow) ) )
              zqih2(i,j,k) = zqih2(i,j,k) - zhfkh * zhfdy
            END IF

! UB>> unnecessary ?
!!$            IF ( lkge2_prog_tke ) THEN
!!$              zhfkh         = ( tkhm(i,j,k) + tkhm(i+1,j,k) ) * acrlat(j,1)
!!$              zhfdy = 0.5 * ( ( tke(i,j  ,k+1,nnow) - tke(i,j  ,k,nnow) )      &
!!$                              + ( tke(i,j+1,k+1,nnow) - tke(i,j+1,k,nnow) ) )
!!$              ztkeh2(i,j,k) = ztkeh2(i,j,k) - zhfkh * zhfdy
!!$            END IF
! UB<<

          ENDDO
        ENDDO

      ELSE IF ( k==ke ) THEN

        DO j = jstart, jendv+1
          DO  i = istart, iendu+1
            ztau11(i,j,k) = ztau11(i,j,k)                                      &
                    - 2.0 * zrhokhm(i,j) * acrlat(j,1)                         &
                     *  dzeta_dlam(i,j,k)                                      &
                     * 0.5 * ( ( u(i-1,j,k,nnow) - u(i-1,j,k-1,nnow) )         &
                             + ( u(i  ,j,k,nnow) - u(i  ,j,k-1,nnow) ) )

            ztau12(i,j,k) =  ztau12(i,j,k)                                     &
                    -0.25_ireals *                                             &
                         ( zrhokhm(i,j  )  + zrhokhm(i+1,j  )                  &
                         + zrhokhm(i,j+1)  + zrhokhm(i+1,j+1) )                &
                     * acrlat(j,2)                                             &
                     * 0.25 * ( dzeta_dlam(i,j  ,k) + dzeta_dlam(i+1,j  ,k)    &
                              + dzeta_dlam(i,j+1,k) + dzeta_dlam(i+1,j+1,k) )  &
                     * 0.5  * ( ( v(i  ,j,k,nnow) - v(i  ,j,k-1,nnow) )        &
                              + ( v(i+1,j,k,nnow) - v(i+1,j,k-1,nnow) ) )

            ztau22(i,j,k) = ztau22(i,j,k)                                      &
                    -2.0 * zrhokhm(i,j) / r_earth                              &
                     *  dzeta_dphi(i,j,k)                                      &
                     *  0.5 *  ( ( v(i,j  ,k,nnow) - v(i,j  ,k-1,nnow) )       &
                               + ( v(i,j-1,k,nnow) - v(i,j-1,k-1,nnow) ) )

          ENDDO
        ENDDO

        ! scalar variables, 1. flux component
        DO j = jstart, jend
          DO i = istart-1, iend

            IF ( lmoist_turb ) THEN
              zhfdx = 0.5 * ( ( ztheta_l(i  ,j,k) - ztheta_l(i  ,j,k-1) )      &
                            + ( ztheta_l(i+1,j,k) - ztheta_l(i+1,j,k-1) ) )
            ELSE
              zhfdx = 0.5 * ( ( ztheta  (i  ,j,k) - ztheta  (i  ,j,k-1) )      &
                            + ( ztheta  (i+1,j,k) - ztheta  (i+1,j,k-1) ) )
            END IF
            zhfkh = 0.5_ireals * ( zpia(i  ,j,k) * zrhokhh(i  ,j)              &
                                 + zpia(i+1,j,k) * zrhokhh(i+1,j) )            &
                                 * acrlat(j,1)                                 &
                         * 0.5 * ( dzeta_dlam(i,j,k) + dzeta_dlam(i+1,j,k) )
            zth1(i,j,k)  = zth1(i,j,k) - zhfkh * zhfdx

            zhfkh = 0.5_ireals * ( zrhokhh(i,j) + zrhokhh(i+1,j) )             &
                                 * acrlat(j,1)                                 &
                         * 0.5 * ( dzeta_dlam(i,j,k) + dzeta_dlam(i+1,j,k) )

            zhfdx = 0.5 * ( ( qv(i  ,j,k,nnow) - qv(i  ,j,k-1,nnow) )          &
                          + ( qv(i+1,j,k,nnow) - qv(i+1,j,k-1,nnow) ) )
            zqvh1(i,j,k) = zqvh1(i,j,k) - zhfkh * zhfdx

            zhfdx = 0.5 * ( ( qc(i  ,j,k,nnow) - qc(i  ,j,k-1,nnow) )          &
                          + ( qc(i+1,j,k,nnow) - qc(i+1,j,k-1,nnow) ) )
            zqch1(i,j,k) = zqch1(i,j,k) - zhfkh * zhfdx

!AK (30.03.2012)
!            DO iig = 1, ntracer
!              zhfdx = 0.5 * ( ( tracer(i  ,j,k,nnow,iig) - tracer(i  ,j,k-1,nnow,iig) )          &
!                            + ( tracer(i+1,j,k,nnow,iig) - tracer(i+1,j,k-1,nnow,iig) ) )
!              ztracerh1(i,j,k,iig) = ztracerh1(i,j,k,iig) - zhfkh * zhfdx
!            ENDDO
!AK (30.03.2012)

            IF ( lprog_qi ) THEN
              zhfdx = 0.5 * ( ( qi(i  ,j,k,nnow) - qi(i  ,j,k-1,nnow) )        &
                            + ( qi(i+1,j,k,nnow) - qi(i+1,j,k-1,nnow) ) )
              zqih1(i,j,k) = zqih1(i,j,k) - zhfkh * zhfdx
            END IF

! UB>>
            IF ( lkge2_prog_tke ) THEN
!!$            IF ( lprog_tke ) THEN
              zhfkh        = ( tkhm(i,j,k) + tkhm(i+1,j,k) ) * acrlat(j,1)
              zhfdx = 0.5 * ( ( tke(i  ,j,k,nnow) - tke(i  ,j,k-1,nnow) )      &
                            + ( tke(i+1,j,k,nnow) - tke(i+1,j,k-1,nnow) ) )
              ztkeh1(i,j,k) = ztkeh1(i,j,k) - zhfkh * zhfdx
            END IF

          ENDDO
        ENDDO

        ! scalar variables, 2. flux component
! UB>>
!!$        DO j = jstart, jend
!!$          DO i = istart-1, iend
        DO j = jstart-1, jend
          DO i = istart, iend
! UB<<

            IF ( lmoist_turb ) THEN
              zhfdy = 0.5 * ( ( ztheta_l(i,j  ,k) - ztheta_l(i,j  ,k-1) )      &
                            + ( ztheta_l(i,j+1,k) - ztheta_l(i,j+1,k-1) ) )
            ELSE
              zhfdy = 0.5 * ( ( ztheta  (i,j  ,k) - ztheta  (i,j  ,k-1) )      &
                            + ( ztheta  (i,j+1,k) - ztheta  (i,j+1,k-1) ) )
            END IF
            zhfkh = 0.5_ireals * ( zpia(i,j  ,k) * zrhokhh(i,j  )              &
                                 + zpia(i,j+1,k) * zrhokhh(i,j+1) )            &
                                 / r_earth                                     &
                         * 0.5 * ( dzeta_dphi(i,j,k) + dzeta_dphi(i+1,j,k) )
            zth2(i,j,k)  = zth2(i,j,k) - zhfkh * zhfdy

            zhfkh = 0.5_ireals * ( zrhokhh(i,j) + zrhokhh(i,j+1) )             &
                                 / r_earth                                     &
                         * 0.5 * ( dzeta_dphi(i,j,k) + dzeta_dphi(i,j+1,k) )

            zhfdy = 0.5 * ( ( qv(i,j  ,k,nnow) - qv(i,j  ,k-1,nnow) )          &
                          + ( qv(i,j+1,k,nnow) - qv(i,j+1,k-1,nnow) ) )
            zqvh2(i,j,k) = zqvh2(i,j,k) - zhfkh * zhfdy

            zhfdy = 0.5 * ( ( qc(i,j  ,k,nnow) - qc(i,j  ,k-1,nnow) )          &
                          + ( qc(i,j+1,k,nnow) - qc(i,j+1,k-1,nnow) ) )
            zqch2(i,j,k) = zqch2(i,j,k) - zhfkh * zhfdy

!AK (30.03.2012)
!            DO iig = 1, ntracer
!              zhfdy = 0.5 * ( ( tracer(i,j  ,k,nnow,iig) - tracer(i,j  ,k-1,nnow,iig) )          &
!                            + ( tracer(i,j+1,k,nnow,iig) - tracer(i,j+1,k-1,nnow,iig) ) )
!              ztracerh2(i,j,k,iig) = ztracerh2(i,j,k,iig) - zhfkh * zhfdy
!            ENDDO
!AK (30.03.2012)

            IF ( lprog_qi ) THEN
              zhfdy = 0.5 * ( ( qi(i,j  ,k,nnow) - qi(i,j  ,k-1,nnow) )        &
                            + ( qi(i,j+1,k,nnow) - qi(i,j+1,k-1,nnow) ) )
              zqih2(i,j,k) = zqih2(i,j,k) - zhfkh * zhfdy
            END IF

! UB>>
            IF ( lkge2_prog_tke ) THEN
!!$            IF ( lprog_tke ) THEN
              zhfkh         = ( tkhm(i,j,k) + tkhm(i+1,j,k) ) * acrlat(j,1)
              zhfdy = 0.5 * ( ( tke(i,j  ,k,nnow) - tke(i,j  ,k-1,nnow) )      &
                            + ( tke(i,j+1,k,nnow) - tke(i,j+1,k-1,nnow) ) )
              ztkeh2(i,j,k) = ztkeh2(i,j,k) - zhfkh * zhfdy
            END IF

          ENDDO
        ENDDO

      END IF

    END IF
    ! end of metric terms for 3D turbulence

  ENDDO

  ! ----- (2) calculation of the tendencies: -------------------------

  DO k = 1, ke

    kup  = MAX( 2, k )
    klow = MIN( ke, k+1 )

! UB>> was missing!
    IF ( k >= 2 .AND. lprog_tke ) THEN
      lkge2_prog_tke = .TRUE.
    ELSE
      lkge2_prog_tke = .FALSE.
    END IF
! UB<<

    DO j = jstart, jend
      zcrlatr = 1.0_ireals / crlat(j,1)
      DO i = istart, iend

        zarhor = 2.0_ireals / ( rho(i,j,kup-1) + rho(i,j,klow-1) ) / r_earth

        wtens(i,j,k)   = wtens(i,j,k) - zarhor *                               &
             ( zcrlatr * ( ztau13(i,j,k) - ztau13(i-1,j  ,k) ) * eddlon        &
                       + ( ztau23(i,j,k) - ztau23(i  ,j-1,k) ) * eddlat )

        zarhor = 1.0_ireals / ( rho(i,j,k) * r_earth )

        ttens(i,j,k)   = ttens(i,j,k) - zarhor *                               &
             ( zcrlatr * ( zth1(i,j,k) - zth1(i-1,j  ,k) ) * eddlon            &
                       + ( zth2(i,j,k) - zth2(i  ,j-1,k) ) * eddlat )
        qvt_diff(i,j,k) = qvt_diff(i,j,k) - zarhor *                           &
             ( zcrlatr * ( zqvh1(i,j,k) - zqvh1(i-1,j  ,k) ) * eddlon          &
                       + ( zqvh2(i,j,k) - zqvh2(i  ,j-1,k) ) * eddlat )
        qctens(i,j,k)  = qctens(i,j,k) - zarhor *                              &
             ( zcrlatr * ( zqch1(i,j,k) - zqch1(i-1,j  ,k) ) * eddlon          &
                       + ( zqch2(i,j,k) - zqch2(i  ,j-1,k) ) * eddlat )
!AK (30.03.2012)
!        DO iig = 1, ntracer
!          tracertens(i,j,k,iig)  = tracertens(i,j,k,iig) - zarhor *                              &
!               ( zcrlatr * ( ztracerh1(i,j,k,iig) - ztracerh1(i-1,j  ,k,iig) ) * eddlon          &
!                         + ( ztracerh2(i,j,k,iig) - ztracerh2(i  ,j-1,k,iig) ) * eddlat )
!        ENDDO
!AK (30.03.2012)
        IF ( lprog_qi ) THEN
          qitens(i,j,k) = qitens(i,j,k) - zarhor *                             &
               ( zcrlatr * ( zqih1(i,j,k) - zqih1(i-1,j  ,k) ) * eddlon        &
                         + ( zqih2(i,j,k) - zqih2(i  ,j-1,k) ) * eddlat )
        END IF
        IF ( lkge2_prog_tke ) THEN
          tketens(i,j,k) = tketens(i,j,k) -                                    &
               ( zcrlatr * ( ztkeh1(i,j,k) - ztkeh1(i-1,j  ,k) ) * eddlon      &
                         + ( ztkeh2(i,j,k) - ztkeh2(i  ,j-1,k) ) * eddlat )    &
                                / r_earth
        END IF

      ENDDO
    ENDDO

    DO j = jstartu, jendu
      zcrlatr = 1.0_ireals / crlat(j,1)
      DO i = istartu, iendu

        zarhor = 2.0_ireals / ( rho(i,j,k) + rho(i+1,j,k) ) / r_earth

        utens(i,j,k)= utens(i,j,k) + zsqrtgrho_r_u(i,j,k) *                    &
                      ( ztaud13(i,j,klow) - ztaud13(i,j,k) )
        utens(i,j,k)= utens(i,j,k) - zarhor *                                  &
             ( zcrlatr * ( ztau11(i+1,j,k) - ztau11(i,j  ,k) ) * eddlon        &
                       + ( ztau12(i  ,j,k) - ztau12(i,j-1,k) ) * eddlat )

      ENDDO
    ENDDO

    DO j = jstartv, jendv
      zcrlatr = 1.0_ireals / crlat(j,2)
      DO i = istartv, iendv

        zarhor = 2.0_ireals / ( rho(i,j,k) + rho(i,j+1,k) ) / r_earth

        vtens(i,j,k)= vtens(i,j,k) + zsqrtgrho_r_v(i,j,k) *                    &
                   ( ztaud23(i,j,klow) - ztaud23(i,j,k) )
        vtens(i,j,k)= vtens(i,j,k) - zarhor *                                  &
             ( zcrlatr * ( ztau12(i,j  ,k) - ztau12(i-1,j,k) ) * eddlon        &
                       + ( ztau22(i,j+1,k) - ztau22(i  ,j,k) ) * eddlat )

      ENDDO
    ENDDO


    ! metric terms for 3D turbulence
    IF ( l3dturb_metr ) THEN

      IF ( ( k>=2 ) .AND. ( k<=ke-1 ) ) THEN

        DO j = jstart, jend
          zcrlatr = 1.0_ireals / crlat(j,1)
          DO i = istart, iend

            zarhor = 2.0_ireals / ( rho(i,j,kup-1)+rho(i,j,klow-1) ) / r_earth

            wtens(i,j,k) = wtens(i,j,k) - zarhor * (                           &
              zcrlatr  * 0.5 * ( dzeta_dlam(i,j,k) + dzeta_dlam(i,j,k-1) )     &
                      * 0.25 * ( ( ztau13(i  ,j,k+1) - ztau13(i  ,j,k-1) )     &
                               + ( ztau13(i-1,j,k+1) - ztau13(i-1,j,k-1) ) )   &
                    +   0.5 *  ( dzeta_dphi(i,j,k) + dzeta_dphi(i,j,k-1) )     &
                      * 0.25 * ( ( ztau23(i,j  ,k+1) - ztau23(i,j  ,k-1) )     &
                               + ( ztau23(i,j-1,k+1) - ztau23(i,j-1,k-1) ) ) )

            zarhor = 1.0_ireals / ( rho(i,j,k) * r_earth )

            ttens(i,j,k) = ttens(i,j,k) - zarhor * (                           &
              zcrlatr *           dzeta_dlam(i,j,k)                            &
                       * 0.25 * ( ( zth1(i  ,j,k+1) - zth1(i  ,j,k-1) )        &
                                + ( zth1(i-1,j,k+1) - zth1(i-1,j,k-1) ) )      &
                     +            dzeta_dphi(i,j,k)                            &
                       * 0.25 * ( ( zth2(i,j  ,k+1) - zth2(i,j  ,k-1) )        &
                                + ( zth2(i,j-1,k+1) - zth2(i,j-1,k-1) ) ) )

            qvt_diff(i,j,k) = qvt_diff(i,j,k) - zarhor * (                     &
              zcrlatr *           dzeta_dlam(i,j,k)                            &
                       * 0.25 * ( ( zqvh1(i  ,j,k+1) - zqvh1(i  ,j,k-1) )      &
                                + ( zqvh1(i-1,j,k+1) - zqvh1(i-1,j,k-1) ) )    &
                     +            dzeta_dphi(i,j,k)                            &
                       * 0.25 * ( ( zqvh2(i,j  ,k+1) - zqvh2(i,j  ,k-1) )      &
                                 + ( zqvh2(i,j-1,k+1) - zqvh2(i,j-1,k-1) ) ) )

            qctens(i,j,k) = qctens(i,j,k) - zarhor * (                         &
              zcrlatr *           dzeta_dlam(i,j,k)                            &
                       * 0.25 * ( ( zqch1(i  ,j,k+1) - zqch1(i  ,j,k-1) )      &
                                + ( zqch1(i-1,j,k+1) - zqch1(i-1,j,k-1) ) )    &
                     +            dzeta_dphi(i,j,k)                            &
                       * 0.25 * ( ( zqch2(i,j  ,k+1) - zqch2(i,j  ,k-1) )      &
                                + ( zqch2(i,j-1,k+1) - zqch2(i,j-1,k-1) ) ) )

!AK (30.03.2012)
!            DO iig = 1, ntracer
!              tracertens(i,j,k,iig) = tracertens(i,j,k,iig) - zarhor * (                         &
!                zcrlatr *          dzeta_dlam(i,j,k)                            &
!                         * 0.25 * ( ( ztracerh1(i  ,j,k+1,iig) - ztracerh1(i  ,j,k-1,iig) )      &
!                                  + ( ztracerh1(i-1,j,k+1,iig) - ztracerh1(i-1,j,k-1,iig) ) )    &
!                       +            dzeta_dphi(i,j,k)                               &
!                         * 0.25 * ( ( ztracerh2(i,j  ,k+1,iig) - ztracerh2(i,j  ,k-1,iig) )      &
!                                  + ( ztracerh2(i,j-1,k+1,iig) - ztracerh2(i,j-1,k-1,iig) ) ) )
!            ENDDO
!AK (30.03.2012)

            IF ( lprog_qi ) THEN
              qitens(i,j,k) = qitens(i,j,k) - zarhor * (                       &
                zcrlatr *        dzeta_dlam(i,j,k)                             &
                        * 0.25 * ( ( zqih1(i  ,j,k+1) - zqih1(i  ,j,k-1) )     &
                                 + ( zqih1(i-1,j,k+1) - zqih1(i-1,j,k-1) ) )   &
                      +          dzeta_dphi(i,j,k)                             &
                        * 0.25 * ( ( zqih2(i,j  ,k+1) - zqih2(i,j  ,k-1) )     &
                                 + ( zqih2(i,j-1,k+1) - zqih2(i,j-1,k-1) ) ) )
            END IF

            IF ( lkge2_prog_tke ) THEN
              tketens(i,j,k) = tketens(i,j,k) - zarhor * (                     &
                zcrlatr *      dzeta_dlam(i,j,k)                               &
                      * 0.25 * ( ( ztkeh1(i  ,j,k+1) - ztkeh1(i  ,j,k-1) )     &
                               + ( ztkeh1(i-1,j,k+1) - ztkeh1(i-1,j,k-1) ) )   &
                    +          dzeta_dphi(i,j,k)                               &
                      * 0.25 * ( ( ztkeh2(i,j  ,k+1) - ztkeh2(i,j  ,k-1) )     &
                               + ( ztkeh2(i,j-1,k+1) - ztkeh2(i,j-1,k-1) ) ) )
            END IF

          ENDDO
        ENDDO

        DO j = jstartu, jendu
          zcrlatr = 1.0_ireals / crlat(j,1)
          DO i = istartu, iendu

            zarhor = 2.0_ireals / ( rho(i,j,k) + rho(i+1,j,k) ) / r_earth

            utens(i,j,k)= utens(i,j,k) - zarhor * (                            &
              zcrlatr * 0.5  * ( dzeta_dlam(i,j,k) + dzeta_dlam(i+1,j,k) )     &
                    * 0.25 * ( ( ztau11(i  ,j,k+1) - ztau11(i  ,j,k-1) )       &
                             + ( ztau11(i+1,j,k+1) - ztau11(i+1,j,k-1) ) )     &
                  +     0.5  * ( dzeta_dphi(i,j,k) + dzeta_dphi(i+1,j,k) )     &
                    * 0.25 * ( ( ztau12(i,j  ,k+1) - ztau12(i,j  ,k-1) )       &
                             + ( ztau12(i,j-1,k+1) - ztau12(i,j-1,k-1) ) ) )

          ENDDO
        ENDDO

        DO j = jstartv, jendv
          zcrlatr = 1.0_ireals / crlat(j,2)
          DO i = istartv, iendv

            zarhor = 2.0_ireals / ( rho(i,j,k) + rho(i,j+1,k) ) / r_earth

            vtens(i,j,k)= vtens(i,j,k) - zarhor * (                            &
              zcrlatr * 0.5  * ( dzeta_dlam(i,j,k) + dzeta_dlam(i,j+1,k) )     &
                    * 0.25 * ( ( ztau12(i  ,j,k+1) - ztau12(i  ,j,k-1) )       &
                             + ( ztau12(i-1,j,k+1) - ztau12(i-1,j,k-1) ) )     &
                  +     0.5  * ( dzeta_dphi(i,j,k) + dzeta_dphi(i,j+1,k) )     &
                    * 0.25 * ( ( ztau22(i,j  ,k+1) - ztau22(i,j  ,k-1) )       &
                             + ( ztau22(i,j+1,k+1) - ztau22(i,j+1,k-1) ) ) )

          ENDDO
        ENDDO

      ELSE IF ( k==1 ) THEN

        DO j = jstart, jend
          zcrlatr = 1.0_ireals / crlat(j,1)
          DO i = istart, iend

            zarhor = 2.0_ireals / ( rho(i,j,kup-1)+rho(i,j,klow-1) ) / r_earth

            wtens(i,j,k) = wtens(i,j,k) - zarhor * (                           &
              zcrlatr          * dzeta_dlam(i,j,k)                             &
                    * 0.5  * ( ( ztau13(i  ,j,k+1) - ztau13(i  ,j,k) )         &
                             + ( ztau13(i-1,j,k+1) - ztau13(i-1,j,k) ) )       &
                  +              dzeta_dphi(i,j,k)                             &
                    * 0.5  * ( ( ztau23(i,j  ,k+1) - ztau23(i,j  ,k) )         &
                             + ( ztau23(i,j-1,k+1) - ztau23(i,j-1,k) ) ) )

            zarhor = 1.0_ireals / ( rho(i,j,k) * r_earth )

            ttens(i,j,k) = ttens(i,j,k) - zarhor * (                           &
              zcrlatr          * dzeta_dlam(i,j,k)                             &
                    * 0.5  * ( ( zth1(i  ,j,k+1) - zth1(i  ,j,k) )             &
                             + ( zth1(i-1,j,k+1) - zth1(i-1,j,k) ) )           &
                  +              dzeta_dphi(i,j,k)                             &
                    * 0.5  * ( ( zth2(i,j  ,k+1) - zth2(i,j  ,k) )             &
                             + ( zth2(i,j-1,k+1) - zth2(i,j-1,k) ) ) )

            qvt_diff(i,j,k) = qvt_diff(i,j,k) - zarhor * (                     &
              zcrlatr          * dzeta_dlam(i,j,k)                             &
                    * 0.5  * ( ( zqvh1(i  ,j,k+1) - zqvh1(i  ,j,k) )           &
                             + ( zqvh1(i-1,j,k+1) - zqvh1(i-1,j,k) ) )         &
                  +              dzeta_dphi(i,j,k)                             &
                    * 0.5  * ( ( zqvh2(i,j  ,k+1) - zqvh2(i,j  ,k) )           &
                             + ( zqvh2(i,j-1,k+1) - zqvh2(i,j-1,k) ) ) )

            qctens(i,j,k) = qctens(i,j,k) - zarhor * (                         &
              zcrlatr          * dzeta_dlam(i,j,k)                             &
                    * 0.5  * ( ( zqch1(i  ,j,k+1) - zqch1(i  ,j,k) )           &
                             + ( zqch1(i-1,j,k+1) - zqch1(i-1,j,k) ) )         &
                  +              dzeta_dphi(i,j,k)                             &
                    * 0.5  * ( ( zqch2(i,j  ,k+1) - zqch2(i,j  ,k) )           &
                             + ( zqch2(i,j-1,k+1) - zqch2(i,j-1,k) ) ) )

!AK (30.03.2012)
!            DO iig = 1, ntracer
!              tracertens(i,j,k,iig) = tracertens(i,j,k,iig) - zarhor * (                         &
!                zcrlatr          * dzeta_dlam(i,j,k)                             &
!                      * 0.5  * ( ( ztracerh1(i  ,j,k+1,iig) - ztracerh1(i  ,j,k,iig) )           &
!                               + ( ztracerh1(i-1,j,k+1,iig) - ztracerh1(i-1,j,k,iig) ) )         &
!                    +              dzeta_dphi(i,j,k)                                &
!                      * 0.5  * ( ( ztracerh2(i,j  ,k+1,iig) - ztracerh2(i,j  ,k,iig) )           &
!                               + ( ztracerh2(i,j-1,k+1,iig) - ztracerh2(i,j-1,k,iig) ) ) )
!            ENDDO
!AK (30.03.2012)

            IF ( lprog_qi ) THEN
              qitens(i,j,k) = qitens(i,j,k) - zarhor * (                       &
                zcrlatr        * dzeta_dlam(i,j,k)                             &
                    * 0.5  * ( ( zqih1(i  ,j,k+1) - zqih1(i  ,j,k) )           &
                             + ( zqih1(i-1,j,k+1) - zqih1(i-1,j,k) ) )         &
                  +              dzeta_dphi(i,j,k)                             &
                    * 0.5  * ( ( zqih2(i,j  ,k+1) - zqih2(i,j  ,k) )           &
                             + ( zqih2(i,j-1,k+1) - zqih2(i,j-1,k) ) ) )
            END IF

! UB>>  unnecessary ?
!!$            IF ( lkge2_prog_tke ) THEN
!!$              tketens(i,j,k) = tketens(i,j,k) - zarhor * (                     &
!!$                zcrlatr        * dzeta_dlam(i,j,k)                             &
!!$                    * 0.5  * ( ( ztkeh1(i  ,j,k+1) - ztkeh1(i  ,j,k) )         &
!!$                             + ( ztkeh1(i-1,j,k+1) - ztkeh1(i-1,j,k) ) )       &
!!$                  +              dzeta_dphi(i,j,k)                             &
!!$                    * 0.5 *  ( ( ztkeh2(i,j  ,k+1) - ztkeh2(i,j  ,k) )         &
!!$                             + ( ztkeh2(i,j-1,k+1) - ztkeh2(i,j-1,k) ) ) )
!!$            END IF
! UB<<
          ENDDO
        ENDDO

        DO j = jstartu, jendu
          zcrlatr = 1.0_ireals / crlat(j,1)
          DO i = istartu, iendu

            zarhor = 2.0_ireals / ( rho(i,j,k) + rho(i+1,j,k) ) / r_earth

            utens(i,j,k) = utens(i,j,k) - zarhor * (                           &
              zcrlatr * 0.5  * ( dzeta_dlam(i,j,k) + dzeta_dlam(i+1,j,k) )     &
                    * 0.5 *  ( ( ztau11(i  ,j,k+1) - ztau11(i  ,j,k) )         &
                             + ( ztau11(i+1,j,k+1) - ztau11(i+1,j,k) ) )       &
                  +     0.5  * ( dzeta_dphi(i,j,k) + dzeta_dphi(i+1,j,k) )     &
                    * 0.5 *  ( ( ztau12(i,j  ,k+1) - ztau12(i,j  ,k) )         &
                             + ( ztau12(i,j-1,k+1) - ztau12(i,j-1,k) ) ) )

          ENDDO
        ENDDO

        DO j = jstartv, jendv
          zcrlatr = 1.0_ireals / crlat(j,2)
          DO i = istartv, iendv

            zarhor = 2.0_ireals / ( rho(i,j,k) + rho(i,j+1,k) ) / r_earth

            vtens(i,j,k) = vtens(i,j,k) - zarhor * (                           &
              zcrlatr * 0.5 *  ( dzeta_dlam(i,j,k) + dzeta_dlam(i,j+1,k) )     &
                    * 0.5  * ( ( ztau12(i  ,j,k+1) - ztau12(i  ,j,k) )         &
                             + ( ztau12(i-1,j,k+1) - ztau12(i-1,j,k) ) )       &
                  +     0.5 *  ( dzeta_dphi(i,j,k) + dzeta_dphi(i,j+1,k) )     &
                    * 0.5  * ( ( ztau22(i,j  ,k+1) - ztau22(i,j  ,k) )         &
                             + ( ztau22(i,j+1,k+1) - ztau22(i,j+1,k) ) ) )

          ENDDO
        ENDDO

      ELSE IF ( k==ke ) THEN

        DO j = jstart, jend
          zcrlatr = 1.0_ireals / crlat(j,1)
          DO i = istart, iend

            zarhor = 2.0_ireals / ( rho(i,j,kup-1)+rho(i,j,klow-1) ) / r_earth

            wtens(i,j,k) = wtens(i,j,k) - zarhor * (                           &
              zcrlatr * 0.5 *  ( dzeta_dlam(i,j,k) + dzeta_dlam(i,j,k-1) )     &
                    * 0.5  * ( ( ztau13(i  ,j,k) - ztau13(i  ,j,k-1) )         &
                             + ( ztau13(i-1,j,k) - ztau13(i-1,j,k-1) ) )       &
                  +     0.5 *  ( dzeta_dphi(i,j,k) + dzeta_dphi(i,j,k-1) )     &
                    * 0.5 *  ( ( ztau23(i,j  ,k) - ztau23(i,j  ,k-1) )         &
                             + ( ztau23(i,j-1,k) - ztau23(i,j-1,k-1) ) ) )

            zarhor = 1.0_ireals / ( rho(i,j,k) * r_earth )

            ttens(i,j,k) = ttens(i,j,k) - zarhor * (                           &
              zcrlatr          * dzeta_dlam(i,j,k)                             &
                    * 0.5  * ( ( zth1(i  ,j,k) - zth1(i  ,j,k-1) )             &
                             + ( zth1(i-1,j,k) - zth1(i-1,j,k-1) ) )           &
                  +              dzeta_dphi(i,j,k)                             &
                    * 0.5  * ( ( zth2(i,j  ,k) - zth2(i,j  ,k-1) )             &
                             + ( zth2(i,j-1,k) - zth2(i,j-1,k-1) ) ) )

            qvt_diff(i,j,k) = qvt_diff(i,j,k) - zarhor * (                     &
              zcrlatr          * dzeta_dlam(i,j,k)                             &
                    * 0.5  * ( ( zqvh1(i  ,j,k) - zqvh1(i  ,j,k-1) )           &
                             + ( zqvh1(i-1,j,k) - zqvh1(i-1,j,k-1) ) )         &
                  +              dzeta_dphi(i,j,k)                             &
                    * 0.5  * ( ( zqvh2(i,j  ,k) - zqvh2(i,j  ,k-1) )           &
                             + ( zqvh2(i,j-1,k) - zqvh2(i,j-1,k-1) ) ) )

            qctens(i,j,k) = qctens(i,j,k) - zarhor * (                         &
              zcrlatr          * dzeta_dlam(i,j,k)                             &
                    * 0.5  * ( ( zqch1(i  ,j,k) - zqch1(i  ,j,k-1) )           &
                             + ( zqch1(i-1,j,k) - zqch1(i-1,j,k-1) ) )         &
                  +              dzeta_dphi(i,j,k)                             &
                    * 0.5  * ( ( zqch2(i,j  ,k) - zqch2(i,j  ,k-1) )           &
                             + ( zqch2(i,j-1,k) - zqch2(i,j-1,k-1) ) ) )

!AK (30.03.2012)
!            DO iig = 1, ntracer
!              tracertens(i,j,k,iig) = tracertens(i,j,k,iig) - zarhor * (                         &
!                zcrlatr          * dzeta_dlam(i,j,k)                             &
!                      * 0.5  * ( ( ztracerh1(i  ,j,k,iig) - ztracerh1(i  ,j,k-1,iig) )           &
!                               + ( ztracerh1(i-1,j,k,iig) - ztracerh1(i-1,j,k-1,iig) ) )         &
!                    +              dzeta_dphi(i,j,k)                                &
!                      * 0.5  * ( ( ztracerh2(i,j  ,k,iig) - ztracerh2(i,j  ,k-1,iig) )           &
!                               + ( ztracerh2(i,j-1,k,iig) - ztracerh2(i,j-1,k-1,iig) ) ) )
!            ENDDO
!AK (30.03.2012)

            IF ( lprog_qi ) THEN
              qitens(i,j,k) = qitens(i,j,k) - zarhor * (                       &
                zcrlatr        * dzeta_dlam(i,j,k)                             &
                    * 0.5  * ( ( zqih1(i  ,j,k) - zqih1(i  ,j,k-1) )           &
                             + ( zqih1(i-1,j,k) - zqih1(i-1,j,k-1) ) )         &
                  +              dzeta_dphi(i,j,k)                             &
                    * 0.5  * ( ( zqih2(i,j  ,k) - zqih2(i,j  ,k-1) )           &
                             + ( zqih2(i,j-1,k) - zqih2(i,j-1,k-1) ) ) )
            END IF

            IF ( lkge2_prog_tke ) THEN
              tketens(i,j,k) = tketens(i,j,k) - zarhor * (                     &
                zcrlatr        * dzeta_dlam(i,j,k)                             &
                    * 0.5  * ( ( ztkeh1(i  ,j,k) - ztkeh1(i  ,j,k-1) )         &
                             + ( ztkeh1(i-1,j,k) - ztkeh1(i-1,j,k-1) ) )       &
                  +              dzeta_dphi(i,j,k)                             &
                    * 0.5  * ( ( ztkeh2(i,j  ,k) - ztkeh2(i,j  ,k-1) )         &
                             + ( ztkeh2(i,j-1,k) - ztkeh2(i,j-1,k-1) ) ) )
            END IF

          ENDDO
        ENDDO

        DO j = jstartu, jendu
          zcrlatr = 1.0_ireals / crlat(j,1)
          DO i = istartu, iendu

            zarhor = 2.0_ireals / ( rho(i,j,k) + rho(i+1,j,k) ) / r_earth

            utens(i,j,k) = utens(i,j,k) - zarhor * (                           &
              zcrlatr * 0.5  * ( dzeta_dlam(i,j,k) + dzeta_dlam(i+1,j,k) )     &
                    * 0.5  * ( ( ztau11(i  ,j,k) - ztau11(i  ,j,k-1) )         &
                             + ( ztau11(i+1,j,k) - ztau11(i+1,j,k-1) ) )       &
                  +     0.5  * ( dzeta_dphi(i,j,k) + dzeta_dphi(i+1,j,k) )     &
                    * 0.5  * ( ( ztau12(i,j  ,k) - ztau12(i,j  ,k-1) )         &
                             + ( ztau12(i,j-1,k) - ztau12(i,j-1,k-1) ) ) )

          ENDDO
        ENDDO

        DO j = jstartv, jendv
          zcrlatr = 1.0_ireals / crlat(j,2)
          DO i = istartv, iendv

            zarhor = 2.0_ireals / ( rho(i,j,k) + rho(i,j+1,k) ) / r_earth

            vtens(i,j,k) = vtens(i,j,k) - zarhor * (                           &
              zcrlatr * 0.5 *  ( dzeta_dlam(i,j,k) + dzeta_dlam(i,j+1,k) )     &
                    * 0.5  * ( ( ztau12(i  ,j,k) - ztau12(i  ,j,k-1) )         &
                             + ( ztau12(i-1,j,k) - ztau12(i-1,j,k-1) ) )       &
                  +     0.5 *  ( dzeta_dphi(i,j,k) + dzeta_dphi(i,j+1,k) )     &
                    * 0.5  * ( ( ztau22(i,j  ,k) - ztau22(i,j  ,k-1) )         &
                             + ( ztau22(i,j+1,k) - ztau22(i,j+1,k-1) ) ) )

          ENDDO
        ENDDO

      END IF

    END IF
    ! end of metric terms for 3D turbulence

  ENDDO

END SUBROUTINE explicit_horizontal_diffusion

!==============================================================================

SUBROUTINE implicit_vert_diffusion_uvwt

  !----------------------------------------------------------------------------
  !
  ! Description:
  !   The vertical diffusion as a slow tendency 
  !   is computed here for all variables the vertical diffusion acts on. 
  !
  ! Method:
  !   The vertical diffusion is solved by a vertically implicit 
  !   scheme (modified Crank-Nicolson).
  !
  !----------------------------------------------------------------------------

  ! Declarations:

  ! Local scalars:
  ! -------------
  INTEGER (KIND=iintegers) ::  &
    i,  j,  k              !  Loop indices in lon., lat. and vert. direction

  INTEGER (KIND=iintegers) :: &
    km1, kp1
  
  REAL    (KIND=ireals   ) ::  &
    zgat, zgct,          & !
    zag, zas,            & !
    zcg, zcs,            & !
    zbg, zdg,            & !
    zd1g, zd2g,          & !
    znew, zz,            & !
    zzdtr, ztmcmq,       & !
    zlhvocp                !

  ! Local (automatic) arrays:
  ! ------------------------
  REAL    (KIND=ireals   ) ::  &
    zkm     (ie,je   ),    & !
    zgatz   (ie,je,ke),    & ! 
    zgctz   (ie,je,ke),    & ! 
    zc      (ie,je,ke),    & ! Upper main diagonal of ...
    zd1     (ie,je,ke),    & ! Right hand side of ...
    zd2     (ie,je,ke),    & ! Right hand side of ...
    ze      (ie,je,ke)       ! Soluton vector of ...

  ! End of header
  !============================================================================

  !----------------------------------------------------------------------------
  ! Begin Subroutine implicit_vert_diffusion_uvwt
  !----------------------------------------------------------------------------

  !----------------------------------------------------------------------------
  ! Section 1: Some preparations
  !----------------------------------------------------------------------------

  ! Setting of reciprocal time step
  zzdtr = 1.0_ireals / dt

  IF ( lvertdiff ) THEN
    
    !--------------------------------------------------------------------------
    ! Section 2: Temperature 
    !--------------------------------------------------------------------------

    IF ( lmoist_turb ) THEN
      
      zlhvocp = lh_v / cp_d
      
      ! Top layer
      DO j = jstart, jend
        DO i = istart, iend
          zgct       = - zkh(i,j,2)*zsqrtgrho_r_s(i,j,1)
          zcg        = zgct*za1t(2)*zpianf(i,j,2)/zpia(i,j,2)
          zbg        = zzdtr - zgct*za1t(2)*zpianf(i,j,2)/zpia(i,j,1)
          zdg        = zzdtr * t(i,j,1,nnow) + ttens(i,j,1)
          zdg        = zdg - zgct*zpianf(i,j,2) * (                          &
                         za2t(2)*  ( ztheta(i,j,2) - ztheta(i,j,1) )         &
                       - zlhvocp * ( qc(i,j,2,nnow)/zpia(i,j,2)              &
                                   - qc(i,j,1,nnow)/zpia(i,j,1) ) )
          zd1(i,j,1) = zdg/zbg
          zc(i,j,1)  = zcg/zbg
        ENDDO
      ENDDO

      ! The layers from k=2 to k=ke-1
      DO k = 2, ke-1
        DO j = jstart, jend
          DO i = istart, iend
            zgat       = - zkh(i,j,k  )*zsqrtgrho_r_s(i,j,k)
            zgct       = - zkh(i,j,k+1)*zsqrtgrho_r_s(i,j,k)
            zag        = zgat*za1t(k  )*zpianf(i,j,k  )/zpia(i,j,k-1)
            zcg        = zgct*za1t(k+1)*zpianf(i,j,k+1)/zpia(i,j,k+1)
            zbg        = zzdtr - zgat*za1t(k  )*zpianf(i,j,k  )/zpia(i,j,k)  &
                               - zgct*za1t(k+1)*zpianf(i,j,k+1)/zpia(i,j,k)
            zdg        = zzdtr * t(i,j,k,nnow) + ttens(i,j,k)
            zdg        = zdg - zgat * zpianf(i,j,k) *  (                     &
                  za2t(k)*  ( ztheta(i,j,k-1) - ztheta(i,j,k) )              &
                - zlhvocp * ( qc(i,j,k-1,nnow)/zpia(i,j,k-1)                 &
                            - qc(i,j,k  ,nnow)/zpia(i,j,k  ) ) )
            zdg        = zdg - zgct * zpianf(i,j,k+1)* (                     &
                  za2t(k+1)*( ztheta(i,j,k+1) - ztheta(i,j,k) )              &
                - zlhvocp * ( qc(i,j,k+1,nnow)/zpia(i,j,k+1)                 &
                            - qc(i,j,k  ,nnow)/zpia(i,j,k  ) ) )
            zz         = 1.0_ireals/( zbg - zag*zc(i,j,k-1) )
            zc (i,j,k) = zcg*zz
            zd1(i,j,k) = ( zdg -zag*zd1(i,j,k-1) )*zz
          ENDDO
        ENDDO
      ENDDO

      ! The bottom layer
      DO j = jstart, jend
        DO i = istart, iend
          zgat       = - zkh (i,j,ke)*zsqrtgrho_r_s(i,j,ke)
          zgct       = - ztch(i,j   )*zsqrtgrho_r_s(i,j,ke)
          zag        = zgat*za1t(ke)*zpianf(i,j,ke)/zpia(i,j,ke-1)
          zbg        = zzdtr - zgat*za1t(ke )*zpianf(i,j,ke )/zpia(i,j,ke)   &
                             - zgct*za1t_surf*zpianf(i,j,ke1)/zpia(i,j,ke)
          zdg        = zzdtr * t(i,j,ke,nnow) + ttens(i,j,ke)
          zdg        = zdg - zgat*zpianf(i,j,ke ) * (                        &
                za2t(ke)* ( ztheta(i,j,ke-1) - ztheta(i,j,ke) )              &
              - zlhvocp * ( qc(i,j,ke-1,nnow)/zpia(i,j,ke-1)                 &
                          - qc(i,j,ke  ,nnow)/zpia(i,j,ke  ) ) )
          zdg        = zdg + zgct*zpianf(i,j,ke1)* (                         &
                za2t_surf * ztheta(i,j,ke)                                   &
              - zlhvocp *   qc(i,j,ke  ,nnow)/zpia(i,j,ke  ) )               &
                           - zgct*t_g(i,j,nnow)
          znew       = ( zdg -zag*zd1(i,j,ke-1) ) / ( zbg - zag*zc(i,j,ke-1) )
          ttens(i,j,ke) = ( znew - t(i,j,ke,nnow) ) * zzdtr
          ze   (i,j,ke) = znew
        ENDDO
      ENDDO

    ELSE

      ! Top layer
      DO j = jstart, jend
        DO i = istart, iend
          zgct       = - zkh(i,j,2)*zsqrtgrho_r_s(i,j,1)
          zcg        = zgct*za1t(2)*zpianf(i,j,2)/zpia(i,j,2)
          zbg        = zzdtr - zgct*za1t(2)*zpianf(i,j,2)/zpia(i,j,1)
          zdg        = zzdtr * t(i,j,1,nnow) + ttens(i,j,1)            &
                     - zgct*za2t(2)*zpianf(i,j,2) *                    &
                       ( ztheta(i,j,2) - ztheta(i,j,1) )
          zd1(i,j,1) = zdg/zbg
          zc(i,j,1)  = zcg/zbg
        ENDDO
      ENDDO

      ! The layers from k=2 to k=ke-1
      DO k = 2, ke-1
        DO j = jstart, jend
          DO i = istart, iend
            zgat       = - zkh(i,j,k  )*zsqrtgrho_r_s(i,j,k)
            zgct       = - zkh(i,j,k+1)*zsqrtgrho_r_s(i,j,k)
            zag        = zgat*za1t(k  )*zpianf(i,j,k  )/zpia(i,j,k-1)
            zcg        = zgct*za1t(k+1)*zpianf(i,j,k+1)/zpia(i,j,k+1)
            zbg        = zzdtr - zgat*za1t(k  )*zpianf(i,j,k  )/zpia(i,j,k)  &
                       - zgct*za1t(k+1)*zpianf(i,j,k+1)/zpia(i,j,k)
            zdg        = zzdtr * t(i,j,k,nnow) + ttens(i,j,k)                &
              - zgat*za2t(k)*zpianf(i,j,k) *                                 &
                ( ztheta(i,j,k-1) - ztheta(i,j,k) )                          &
              - zgct*za2t(k+1)*zpianf(i,j,k+1)*                              &
                ( ztheta(i,j,k+1) - ztheta(i,j,k) )
            zz         = 1.0_ireals/( zbg - zag*zc(i,j,k-1) )
            zc (i,j,k) = zcg*zz
            zd1(i,j,k) = ( zdg -zag*zd1(i,j,k-1) )*zz
          ENDDO
        ENDDO
      ENDDO

      ! The bottom layer
      DO j = jstart, jend
        DO i = istart, iend
          zgat       = - zkh (i,j,ke)*zsqrtgrho_r_s(i,j,ke)
          zgct       = - ztch(i,j   )*zsqrtgrho_r_s(i,j,ke)
          zag        = zgat*za1t(ke)*zpianf(i,j,ke)/zpia(i,j,ke-1)
          zbg        = zzdtr - zgat*za1t(ke )*zpianf(i,j,ke )/zpia(i,j,ke)   &
            - zgct*za1t_surf*zpianf(i,j,ke1)/zpia(i,j,ke)
          zdg        = zzdtr * t(i,j,ke,nnow) + ttens(i,j,ke)                &
            - zgat*za2t(ke)*zpianf(i,j,ke ) *                                &
              ( ztheta(i,j,ke-1) - ztheta(i,j,ke) )                          &
            + zgct*za2t_surf *zpianf(i,j,ke1)*                               &
              ztheta(i,j,ke) - zgct*t_g(i,j,nnow)
          znew       = ( zdg -zag*zd1(i,j,ke-1) ) / ( zbg - zag*zc(i,j,ke-1) )
          ttens(i,j,ke) = ( znew - t(i,j,ke,nnow) ) * zzdtr
          ze   (i,j,ke) = znew
        ENDDO
      ENDDO

    END IF
  
    ! Backsubstitution and storage of the complete slow tendencies

    DO k = ke-1, 1, -1
      DO j = jstart, jend
        DO i = istart, iend
          ze   (i,j,k) = zd1(i,j,k) - zc(i,j,k)*ze(i,j,k+1)
          ttens(i,j,k) = ( ze(i,j,k) - t(i,j,k,nnow) ) * zzdtr
        ENDDO
      ENDDO
    ENDDO

    ! Calculation of the sensible heat flux at the surface
    ! This flux is integrated in time

    DO j = jstart, jend
      DO  i  = istart , iend
        shfl_s (i,j) = - ztch(i,j)*cp_d*                  &
          ( za2t_surf*(t_g(i,j,nnow) - zpianf(i,j,ke1)*   &
                       t  (i,j,ke,nnow)/zpia(i,j,ke) ) +  &
            za1t_surf*(t_g(i,j,nnew) - zpianf(i,j,ke1)*   &
                       ze(i,j,ke)/zpia(i,j,ke) ) )
      ENDDO
    ENDDO

    !--------------------------------------------------------------------------
    ! Section : Horizontal wind velocity u 
    !--------------------------------------------------------------------------

    ! Top layer:   k = 1
    DO j = jstartu, jendu
      DO i = istartu, iendu
        zkm(i,j)   = 0.5*( ztmkvm(i,j,2) + ztmkvm(i+1,j,2) )
        zgct       = - zkm(i,j)*zsqrtgrho_r_u(i,j,1)
        zcg        = zgct*za1t(2)
        zbg        = zzdtr - zcg
        zdg        = zzdtr * u(i,j,1,nnow) + utens(i,j,1)    &
          - za2t(2)*zgct * ( u(i,j,2,nnow) - u(i,j,1,nnow) )
        zd1(i,j,1) = zdg/zbg
        zc (i,j,1) = zcg/zbg
      ENDDO
    ENDDO

    ! The layers from k=2 to k=ke-1
    DO k = 2, ke-1
      DO j = jstartu, jendu
        DO i = istartu, iendu
          zgat       = - zkm(i,j)*zsqrtgrho_r_u(i,j,k)
          zkm(i,j)   = 0.5*( ztmkvm(i,j,k+1) + ztmkvm(i+1,j,k+1) )
          zgct       = - zkm(i,j)*zsqrtgrho_r_u(i,j,k)
          zag        = zgat*za1t(k)
          zas        = zgat*za2t(k)
          zcg        = zgct*za1t(k+1)
          zcs        = zgct*za2t(k+1)
          zbg        = zzdtr - zag - zcg
          zdg        = zzdtr * u(i,j,k,nnow) + utens(i,j,k)  &
            - zas * ( u(i,j,k-1,nnow) - u(i,j,k,nnow) )      &
            - zcs * ( u(i,j,k+1,nnow) - u(i,j,k,nnow) )
          zz         = 1.0_ireals/( zbg - zag*zc(i,j,k-1) )
          zc (i,j,k)  = zcg*zz
          zd1(i,j,k)  = ( zdg -zag*zd1(i,j,k-1) )*zz
        ENDDO
      ENDDO
    ENDDO

    ! The bottom layer:  k = ke
    ! Including the calculation of the u-momentum flux at the surface
    DO j = jstartu, jendu
      DO i = istartu, iendu
        zgat       = - zkm(i,j)*zsqrtgrho_r_u(i,j,ke)
        ztmcmq     = 0.5_ireals*( ztcm(i,j) + ztcm(i+1,j) )
        zgct       = - ztmcmq*zsqrtgrho_r_u(i,j,ke)
        zag        = zgat*za1t(ke)
        zas        = zgat*za2t(ke)
        zcg        = za1t(ke1)*zgct
        zcs        = za2t(ke1)*zgct
        zbg        = zzdtr - zag - zcg
        zdg        = zzdtr * u(i,j,ke,nnow) + utens(i,j,ke)  &
          - zas * ( u(i,j,ke-1,nnow) - u(i,j,ke,nnow) )      &
          + zcs * u(i,j,ke,nnow)
        znew       = ( zdg -zag*zd1(i,j,ke-1) ) / ( zbg - zag*zc(i,j,ke-1) )
        utens(i,j,ke) = ( znew - u(i,j,ke,nnow) ) * zzdtr
        umfl_s(i,j  ) = ztmcmq*( a2t(ke1)*u(i,j,ke,nnow) + a1t(ke1)*znew )
        ze   (i,j,ke) = znew
      ENDDO
    ENDDO

    ! Backsubstitution and storage of the complete slow tendencies

    DO k = ke-1, 1, -1
      DO j = jstartu, jendu
        DO i = istartu, iendu
          ze   (i,j,k) = zd1(i,j,k) - zc(i,j,k)*ze(i,j,k+1)
          utens(i,j,k) = ( ze(i,j,k) - u(i,j,k,nnow) ) * zzdtr
        ENDDO
      ENDDO
    ENDDO

    !--------------------------------------------------------------------------
    ! Section : Horizontal wind velocity v
    !--------------------------------------------------------------------------

    ! Top layer  k=1
    DO j = jstartv, jendv
      DO i = istartv, iendv
        zkm(i,j)   = 0.5*( ztmkvm(i,j,2) + ztmkvm(i,j+1,2) )
        zgct       = - zkm(i,j)*zsqrtgrho_r_v(i,j,1)
        zcg        = zgct*za1t(2)
        zbg        = zzdtr - zcg
        zdg        = zzdtr * v(i,j,1,nnow) + vtens(i,j,1)     &
          - za2t(2)*zgct * ( v(i,j,2,nnow) - v(i,j,1,nnow) )
        zd1(i,j,1) = zdg/zbg
        zc (i,j,1) = zcg/zbg
      ENDDO
    ENDDO

    ! The layers from k=2 to k=ke-1
    DO  k = 2, ke-1
      DO j = jstartv, jendv
        DO i = istartv, iendv
          zgat       = - zkm(i,j)*zsqrtgrho_r_v(i,j,k)
          zkm(i,j)   = 0.5*( ztmkvm(i,j,k+1) + ztmkvm(i,j+1,k+1) )
          zgct       = - zkm(i,j)*zsqrtgrho_r_v(i,j,k)
          zag        = zgat*za1t(k)
          zas        = zgat*za2t(k)
          zcg        = zgct*za1t(k+1)
          zcs        = zgct*za2t(k+1)
          zbg        = zzdtr - zag - zcg
          zdg        = zzdtr * v(i,j,k,nnow) + vtens(i,j,k)  &
            - zas * ( v(i,j,k-1,nnow) - v(i,j,k,nnow) )      &
            - zcs * ( v(i,j,k+1,nnow) - v(i,j,k,nnow) )
          zz         = 1.0_ireals/( zbg - zag*zc(i,j,k-1) )
          zc (i,j,k) = zcg*zz
          zd1(i,j,k) = ( zdg -zag*zd1(i,j,k-1) )*zz
        ENDDO
      ENDDO
    ENDDO

    ! The bottom layer k=ke
    ! Including the calculation of the v-momentum flux at the surface
    DO j = jstartv, jendv
      DO i = istartv, iendv
        zgat       = - zkm(i,j)*zsqrtgrho_r_v(i,j,ke)
        ztmcmq     = 0.5_ireals*( ztcm(i,j) + ztcm(i,j+1) )
        zgct       = - ztmcmq*zsqrtgrho_r_v(i,j,ke)
        zag        = zgat*za1t(ke)
        zas        = zgat*za2t(ke)
        zcg        = za1t(ke1)*zgct
        zcs        = za2t(ke1)*zgct
        zbg        = zzdtr - zag - zcg
        zdg        = zzdtr * v(i,j,ke,nnow) + vtens(i,j,ke)  &
          - zas * ( v(i,j,ke-1,nnow) - v(i,j,ke,nnow) )      &
          + zcs * v(i,j,ke,nnow)
        znew       = ( zdg -zag*zd1(i,j,ke-1) ) / ( zbg - zag*zc(i,j,ke-1) )
        vtens(i,j,ke) = ( znew - v(i,j,ke,nnow) ) * zzdtr
        vmfl_s(i,j  ) = ztmcmq*( a2t(ke1)*v(i,j,ke,nnow) + a1t(ke1)*znew )
        ze   (i,j,ke) = znew
      ENDDO
    ENDDO

    ! Backsubstitution and storage of the complete slow tendencies

    DO k = ke-1, 1, -1
      DO j = jstartv, jendv
        DO i = istartv, iendv
          ze   (i,j,k) = zd1(i,j,k) - zc(i,j,k)*ze(i,j,k+1)
          vtens(i,j,k) = ( ze(i,j,k) - v(i,j,k,nnow) ) * zzdtr
        ENDDO
      ENDDO
    ENDDO

  END IF
  
  !----------------------------------------------------------------------------
  ! Section : Vertical wind velocity w
  !----------------------------------------------------------------------------

  IF ( lvertdiff_w ) THEN
    
    DO j = jstart, jend

      DO  k = 2, ke
        km1 = MAX( 2, k-1 )
        kp1 = MIN( ke, k+1 )
      
        ! Precalculate some help variables to avoid divisions in subsequent code
        DO i = istart, iend
          zgatz(i,j,k) = - zsqrtgrho_r_w(i,j,k)   &
                         * ( ztmkvw(i,j,k)+ztmkvw(i,j,km1) )
          zgctz(i,j,k) = - zsqrtgrho_r_w(i,j,k)   &
                         * ( ztmkvw(i,j,kp1)+ztmkvw(i,j,k) )
        ENDDO
      ENDDO
      
      ! Top layer       
      DO i = istart, iend
        zag        = 0.5*(a1t(2)+a1t(1)) * zgatz(i,j,2)
        zas        = 0.5*(a2t(2)+a2t(1)) * zgatz(i,j,2)
        zcg        = 0.5*(a1t(3)+a1t(2)) * zgctz(i,j,2)
        zcs        = 0.5*(a2t(3)+a2t(2)) * zgctz(i,j,2)
        zbg        = zzdtr - zag - zcg
        zdg        = zzdtr * w(i,j,2,nnow) + wtens(i,j,2)                    &
                     + zas * w(i,j,2,nnow)                                   &
                     - zcs * ( w(i,j,3,nnow) - w(i,j,2,nnow) )
        zd1(i,j,2)  = zdg/zbg
        zc (i,j,2)  = zcg/zbg
      ENDDO

      ! The layers from k=3 to k=ke-1
      DO k = 3, ke-1
        DO i = istart, iend
          zag        = 0.5*(a1t(k)+a1t(k-1)) * zgatz(i,j,k)
          zas        = 0.5*(a2t(k)+a2t(k-1)) * zgatz(i,j,k)
          zcg        = 0.5*(a1t(k+1)+a1t(k)) * zgctz(i,j,k)
          zcs        = 0.5*(a2t(k+1)+a2t(k)) * zgctz(i,j,k)
          zbg        = zzdtr - zag - zcg
          zdg        = zzdtr * w(i,j,k,nnow) + wtens(i,j,k)                  &
                     - zas * ( w(i,j,k-1,nnow) - w(i,j,k,nnow) )             &
                     - zcs * ( w(i,j,k+1,nnow) - w(i,j,k,nnow) )
          zz         = 1.0_ireals/( zbg - zag*zc(i,j,k-1) )
          zc (i,j,k)  = zcg*zz
          zd1(i,j,k)  = ( zdg - zag*zd1(i,j,k-1) )*zz
        ENDDO
      ENDDO

      ! The bottom layer
      DO i = istart, iend
        zag        = 0.5*(a1t(ke)+a1t(ke-1)) * zgatz(i,j,ke)
        zas        = 0.5*(a2t(ke)+a2t(ke-1)) * zgatz(i,j,ke)
        zcg        = 0.5*(a1t(ke1)+a1t(ke))  * zgctz(i,j,ke)
        zcs        = 0.5*(a2t(ke1)+a2t(ke))  * zgctz(i,j,ke)
        zbg        = zzdtr - zag - zcg
        zdg        = zzdtr * w(i,j,ke,nnow) + wtens(i,j,ke)                 &
                     - zas * ( w(i,j,ke-1,nnow) - w(i,j,ke,nnow) )          &
                     - zcs * ( w(i,j,ke1,nnow) - w(i,j,ke,nnow))            &
                     - zcg * w(i,j,ke1,nnow) 
        znew       = ( zdg -zag*zd1(i,j,ke-1) ) / ( zbg - zag*zc(i,j,ke-1) )
        wtens(i,j,ke) = ( znew - w(i,j,ke,nnow) ) * zzdtr
        ze   (i,j,ke) = znew
      ENDDO

      ! Backsubstitution and storage of the complete slow tendencies

      DO k = ke-1, 2, -1
        DO i = istart, iend
          ze    (i,j,k) = zd1(i,j,k) - zc(i,j,k)*ze(i,j,k+1)
          wtens (i,j,k) = ( ze(i,j,k) - w(i,j,k,nnow) ) * zzdtr
        ENDDO
      ENDDO
    ENDDO

  END IF

  !----------------------------------------------------------------------------
  ! End of subroutine implicit_vert_diffusion_uvwt
  !----------------------------------------------------------------------------

END SUBROUTINE implicit_vert_diffusion_uvwt

!==============================================================================



SUBROUTINE complete_tend_uvwtpp_CN3Crow( dt_adv, nstar )

  !----------------------------------------------------------------------------
  !
  ! Description:
  !   This procedure calculates the tendencies of the vertical advection 
  !   (uadvt, vadvt, wadvt, ppadvt and tadvt)
  !   for the prognostic variables u,v,w,pp and T.
  !
  ! Method:
  !   the vertical advection is solved by a vertically implicit 
  !   scheme (modified Crank-Nicolson 3. order) starting at the states
  !   u(:,:,:,nstar), v(:,:,:,now), w(..,nstar), pp(..), T(..)
  !
  !----------------------------------------------------------------------------

  USE numeric_utilities, ONLY: solve_5banddiag, solve_5banddiag_vec

  USE data_parallel,      ONLY :  &
    my_cart_id     ! rank of this subdomain in the cartesian communicator


  IMPLICIT NONE

  ! Declarations:

  ! Subroutine arguments:
  ! ---------------------

  INTEGER (KIND=iintegers), INTENT(in) ::  &
    nstar

  REAL (KIND=ireals), INTENT(in) ::  &
    dt_adv


  ! Local scalars:
  ! -------------
  INTEGER (KIND=iintegers) ::  &
    i,  j,  k              !  Loop indices in lon., lat. and vert. direction

  REAL    (KIND=ireals   ) ::  &
    z_beta_v, z_beta_p, z_beta_m,     &  ! Crank-Nicholson-weights
    zC2, &           ! square of local Courant number
    zh1, zh2, &      ! help variable for storage of intermediate results
    z_dt_half, &     ! = dt_adv / 2
    z_dt_quart, &    ! = dt_adv / 4
    z_dt_recip, &    ! = 1 / dt
    z_beta_p_half, & ! = z_beta_p / 2
    z_beta_m_half    ! = z_beta_m / 2


  ! Local (automatic) arrays:
  ! ------------------------

  ! Arrays for the LGS with a  5-band-diagonal-matrix:
  REAL    (KIND=ireals   ), ALLOCATABLE ::  &
    z_Cour    (:,:,:),   &
    z_lgs     (:,:,:,:), &  ! the 5 banddiagonals of the lin. eq. system
    z_lgs_rhs (:,:,:),   &
    zh        (:,:,:)

  REAL    (KIND=ireals   ), ALLOCATABLE ::  &
    z_lgs_store (:,:,:,:) ! copy of the 5 banddiagonals of the lin. eq. system

  ! nicht sehr elegant, dass man hier extra Felder deklarieren muss:
  REAL    (KIND=ireals   ), ALLOCATABLE ::   &
    z_Cour_w    (:,:,:),   &
    z_lgs_w     (:,:,:,:), & ! the 5 banddiagonals of the lin. eq. system
    z_lgs_rhs_w (:,:,:),   &
    zh_w        (:,:,:)

  INTEGER (KIND=iintegers) ::  istat

  LOGICAL :: flag_vector_version  ! call Numerical Recipes routines, 
                                  ! which are optimized for vector computers

  ! End of header
  !============================================================================

  !----------------------------------------------------------------------------
  ! Begin Subroutine complete_tend_uvwtpp_CN3Crow
  !----------------------------------------------------------------------------

  !----------------------------------------------------------------------------
  ! Section 1: Some preparations
  !----------------------------------------------------------------------------

  flag_vector_version = .TRUE.

  ALLOCATE( z_Cour     (1:ie, 1:je, 1:ke), STAT=istat )
  ALLOCATE( z_lgs      (1:ie, 1:je, 1:ke, 1:5), STAT=istat )
  ALLOCATE( z_lgs_rhs  (1:ie, 1:je, 1:ke), STAT=istat )
  ALLOCATE( zh         (1:ie, 1:je, 1:ke), STAT=istat )

  ALLOCATE( z_Cour_w   (1:ie, 1:je, 1:ke1), STAT=istat )
  ALLOCATE( z_lgs_w    (1:ie, 1:je, 1:ke1, 1:5), STAT=istat )
  ALLOCATE( z_lgs_rhs_w(1:ie, 1:je, 1:ke1), STAT=istat )
  ALLOCATE( zh_w       (1:ie, 1:je, 1:ke1), STAT=istat )

  !IF ( my_cart_id == 0 ) THEN
  !  WRITE(*,*) "Subr. complete_tend_uvwtpp_CN3Crow ..."
  !  WRITE(*,*) "  flag_vector_version = ", flag_vector_version
  !END IF


  !MB: for debugging:
  !z_Cour     = 0.0
  !z_lgs      = 0.0
  !z_lgs_rhs  = 0.0
  !zh         = 0.0

  !z_Cour_w     = 0.0
  !z_lgs_w      = 0.0
  !z_lgs_rhs_w  = 0.0
  !zh_w         = 0.0

  ! precalculated factors for optimization 
  z_dt_half  = 0.5_ireals  * dt_adv 
  z_dt_quart = 0.25_ireals * dt_adv 
  z_dt_recip = 1.0_ireals / dt_adv

  ! Setting of parameters for implicit calculation of vertical advection       
  ! (z_beta_v: weight for the n+1 time level           
  !          z_beta_v=0: centered, =+1: full implicit,=-1: explicit)
  z_beta_v = 0.0
  z_beta_p = 0.5*( 1.0 + z_beta_v )
  z_beta_m = 0.5*( 1.0 - z_beta_v )

  z_beta_p_half = 0.5_ireals * z_beta_p
  z_beta_m_half = 0.5_ireals * z_beta_m

  !----------------------------------------------------------------------------
  ! Section 2: Setup of tridiagonal matrix systems resulting from the implicit
  !            numerical formulation of advection.
  !            - variables at scalar positions (T and pp)
  !----------------------------------------------------------------------------

  ! Band-diagonal-matrix

  ! Top layer k=1 (one-sided diff.)  
  DO j = jstart, jend
    DO i = istart, iend

      ! remark: wcon contains 'contravar. vertical velocity * sqrt(G)':
      z_Cour(i,j,1) = sqrtg_r_s(i,j,1) * ( wcon(i,j,1) + wcon(i,j,2) ) * z_dt_half 

      IF ( z_Cour(i,j,1) < 0.0_ireals ) THEN
      z_lgs(i,j,1,1) = 0.0_ireals
      z_lgs(i,j,1,2) = 0.0_ireals
      z_lgs(i,j,1,3) = 1.0_ireals - z_Cour(i,j,1) * z_beta_p
      z_lgs(i,j,1,4) =            + z_Cour(i,j,1) * z_beta_p
      z_lgs(i,j,1,5) = 0.0_ireals
      ELSE
        ! use explicit upwind:
        z_lgs(i,j,1,1) = 0.0_ireals
        z_lgs(i,j,1,2) = 0.0_ireals
        z_lgs(i,j,1,3) = 1.0_ireals
        z_lgs(i,j,1,4) = 0.0_ireals
        z_lgs(i,j,1,5) = 0.0_ireals
      END IF

    ENDDO
  ENDDO

  ! 2. Layer k=2 (centered diff. 2. order)
  DO j = jstart, jend
    DO i = istart, iend

      z_Cour(i,j,2) = sqrtg_r_s(i,j,2) * ( wcon(i,j,2) + wcon(i,j,3) ) * z_dt_half 

      z_lgs(i,j,2,1) = 0.0_ireals
      z_lgs(i,j,2,2) = - z_Cour(i,j,2) * z_beta_p_half
      z_lgs(i,j,2,3) = 1.0_ireals
      z_lgs(i,j,2,4) = - z_lgs(i,j,2,2)
      z_lgs(i,j,2,5) = 0.0_ireals

    ENDDO
  ENDDO

  ! The layers from k=3 to k=ke-2
  DO k = 3, ke-2
    DO j = jstart, jend
      DO i = istart, iend

        z_Cour(i,j,k) = sqrtg_r_s(i,j,k) * ( wcon(i,j,k) + wcon(i,j,k+1) ) * z_dt_half 
        zC2 = z_Cour(i,j,k) * z_Cour(i,j,k)
        zh2 = z_Cour(i,j,k) * z_beta_p

        IF ( z_Cour(i,j,k) > 0 ) THEN
          z_lgs(i,j,k,1) = zh2 * ( 1.0/6.0 + zC2/12.0 )
          z_lgs(i,j,k,2) = zh2 * ( - 1.0   - zC2/4.0  )
          z_lgs(i,j,k,3) = zh2 * ( 1.0/2.0 + zC2/4.0  ) + 1.0
          z_lgs(i,j,k,4) = zh2 * ( 1.0/3.0 - zC2/12.0 )
          z_lgs(i,j,k,5) = 0.0_ireals
        ELSE 
          z_lgs(i,j,k,1) = 0.0_ireals
          z_lgs(i,j,k,2) = - zh2 * ( 1.0/3.0 - zC2/12.0 )
          z_lgs(i,j,k,3) = - zh2 * ( 1.0/2.0 + zC2/4.0  ) + 1.0
          z_lgs(i,j,k,4) = - zh2 * ( - 1.0   - zC2/4.0  )
          z_lgs(i,j,k,5) = - zh2 * ( 1.0/6.0 + zC2/12.0 )
        END IF

      ENDDO
    ENDDO
  ENDDO

  ! Layer k=ke-1 (centered diff. 2. order)
  DO j = jstart, jend
    DO i = istart, iend

      z_Cour(i,j,ke-1) = sqrtg_r_s(i,j,ke-1) * ( wcon(i,j,ke-1) + wcon(i,j,ke) ) * z_dt_half 

      z_lgs(i,j,ke-1,1) = 0.0_ireals
      z_lgs(i,j,ke-1,2) = - z_Cour(i,j,ke-1) * z_beta_p_half
      z_lgs(i,j,ke-1,3) = 1.0_ireals
      z_lgs(i,j,ke-1,4) = - z_lgs(i,j,ke-1,2)
      z_lgs(i,j,ke-1,5) = 0.0_ireals

    ENDDO
  ENDDO

  ! The bottom layer k=ke
  DO j = jstart, jend
    DO i = istart, iend

      z_Cour(i,j,ke) = sqrtg_r_s(i,j,ke) * ( wcon(i,j,ke) + wcon(i,j,ke+1) ) * z_dt_half  

      IF ( z_Cour(i,j,k) > 0 ) THEN
        z_lgs(i,j,ke,1) = 0.0_ireals
        z_lgs(i,j,ke,2) =            - z_Cour(i,j,ke) * z_beta_p
        z_lgs(i,j,ke,3) = 1.0_ireals + z_Cour(i,j,ke) * z_beta_p
        z_lgs(i,j,ke,4) = 0.0_ireals
        z_lgs(i,j,ke,5) = 0.0_ireals
      ELSE
        ! use explicit upwind:
        z_lgs(i,j,ke,1) = 0.0_ireals
        z_lgs(i,j,ke,2) = 0.0_ireals
        z_lgs(i,j,ke,3) = 1.0_ireals
        z_lgs(i,j,ke,4) = 0.0_ireals
        z_lgs(i,j,ke,5) = 0.0_ireals
      END IF
    ENDDO
  ENDDO

  ! --- right hand side for T ---

  ! Top layer (k=1)
  DO j = jstart, jend
    DO i = istart, iend
      IF ( z_Cour(i,j,1) < 0.0_ireals ) THEN
      z_lgs_rhs(i,j,1) = T(i,j,1,nstar)                                   &
        &  - z_Cour(i,j,1) * z_beta_m * ( T(i,j,2,nstar)-T(i,j,1,nstar) ) 
      ELSE
        ! use explicit upwind:
        z_lgs_rhs(i,j,1) = T(i,j,1,nstar)                                   &
          &  - z_Cour(i,j,1) * 1.0_ireals * ( T(i,j,2,nstar)-T(i,j,1,nstar) ) 
      END IF
    ENDDO
  ENDDO

  ! 2. layer k=2
  DO j = jstart, jend
    DO i = istart, iend
      z_lgs_rhs(i,j,2) = T(i,j,2,nstar)                                          &
        &  - z_Cour(i,j,2) * ( T(i,j,3,nstar) - T(i,j,1,nstar) ) * z_beta_m_half 

    ENDDO
  ENDDO

  ! The layers from k=3 to k=ke-2
  DO k = 3, ke-2
    DO j = jstart, jend
      DO i = istart, iend

        zC2 = z_Cour(i,j,k) * z_Cour(i,j,k) 

        ! store advection-operator in help-variable zh1:
        IF ( z_Cour(i,j,k) > 0 ) THEN
          zh1 =   ( 1.0/6.0 + zC2/12.0 ) * T(i,j,k-2,nstar)  &
              & + ( - 1.0   - zC2/4.0  ) * T(i,j,k-1,nstar)  &
              & + ( 1.0/2.0 + zC2/4.0  ) * T(i,j,k  ,nstar)  &
              & + ( 1.0/3.0 - zC2/12.0 ) * T(i,j,k+1,nstar)
        ELSE 
          zh1 = - ( 1.0/6.0 + zC2/12.0 ) * T(i,j,k+2,nstar)  &
              & - ( - 1.0   - zC2/4.0  ) * T(i,j,k+1,nstar)  &
              & - ( 1.0/2.0 + zC2/4.0  ) * T(i,j,k  ,nstar)  &
              & - ( 1.0/3.0 - zC2/12.0 ) * T(i,j,k-1,nstar)
        END IF
        z_lgs_rhs(i,j,k) = T(i,j,k,nstar)                  &
          &   - z_Cour(i,j,k) * zh1 * z_beta_m  

      ENDDO
    ENDDO
  ENDDO

  ! layer k=ke-1
  DO j = jstart, jend
    DO i = istart, iend
      z_lgs_rhs(i,j,ke-1) = T(i,j,ke-1,nstar)                                          &
        &   - z_Cour(i,j,ke-1) * ( T(i,j,ke,nstar)-T(i,j,ke-2,nstar) ) * z_beta_m_half
    ENDDO
  ENDDO

  ! Bottom layer k=ke
  DO j = jstart, jend
    DO i = istart, iend
      IF ( z_Cour(i,j,k) > 0 ) THEN
        z_lgs_rhs(i,j,ke) = T(i,j,ke,nstar)              &
          &  - z_Cour(i,j,ke) * z_beta_m * ( T(i,j,ke,nstar)-T(i,j,ke-1,nstar) )
      ELSE
        ! use explicit upwind:
        z_lgs_rhs(i,j,ke) = T(i,j,ke,nstar)              &
          &  - z_Cour(i,j,ke) * 1.0_ireals * ( T(i,j,ke,nstar)-T(i,j,ke-1,nstar) )
      END IF
    ENDDO
  ENDDO

  IF ( flag_vector_version ) THEN 
    ALLOCATE( z_lgs_store( 1:ie, 1:je, 1:ke, 1:5 ), STAT=istat )

    z_lgs_store(:,:,:,:) = z_lgs(:,:,:,:)

    CALL solve_5banddiag_vec( z_lgs_store, z_lgs_rhs, zh, &
         ie, je, istart, iend, jstart, jend, ke ) 
    ! (remark: z_lgs_store now has changed)

    DEALLOCATE( z_lgs_store )
  ELSE

    CALL solve_5banddiag( z_lgs, z_lgs_rhs, zh, &
         ie, je, istart, iend, jstart, jend, ke ) 
  END IF

  DO k = 1, ke
    DO j = jstart, jend
      DO i = istart, iend
        tadvt(i,j,k) = ( zh(i,j,k) - T(i,j,k, nstar) ) * z_dt_recip  
      ENDDO
    ENDDO
  ENDDO

  ! --- right hand side for pp ---

  ! Top layer (k=1)      
  DO j = jstart, jend
    DO i = istart, iend
      IF ( z_Cour(i,j,1) < 0.0_ireals ) THEN
        z_lgs_rhs(i,j,1) = pp(i,j,1,nstar)                                   &
          &  - z_Cour(i,j,1) * z_beta_m * ( pp(i,j,2,nstar)-pp(i,j,1,nstar) )
      ELSE
        ! use explicit upwind:
        z_lgs_rhs(i,j,1) = pp(i,j,1,nstar)                                   &
          &  - z_Cour(i,j,1) * 1.0_ireals * ( pp(i,j,2,nstar)-pp(i,j,1,nstar) )
      END IF
    ENDDO
  ENDDO

  ! 2. layer k=2
  DO j = jstart, jend
    DO i = istart, iend
      z_lgs_rhs(i,j,2) = pp(i,j,2,nstar)                                           &
        &  - z_Cour(i,j,2) * ( pp(i,j,3,nstar) - pp(i,j,1,nstar) ) * z_beta_m_half

    ENDDO
  ENDDO

  ! The layers from k=3 to k=ke-2
  DO k = 3, ke-2
    DO j = jstart, jend
      DO i = istart, iend

        zC2 = z_Cour(i,j,k) * z_Cour(i,j,k) 

        ! store advection-operator in help-variable zh1:
        IF ( z_Cour(i,j,k) > 0 ) THEN
          zh1 =   ( 1.0/6.0 + zC2/12.0 ) * pp(i,j,k-2,nstar)  &
              & + ( - 1.0   - zC2/4.0  ) * pp(i,j,k-1,nstar)  &
              & + ( 1.0/2.0 + zC2/4.0  ) * pp(i,j,k  ,nstar)  &
              & + ( 1.0/3.0 - zC2/12.0 ) * pp(i,j,k+1,nstar)
        ELSE 
          zh1 = - ( 1.0/6.0 + zC2/12.0 ) * pp(i,j,k+2,nstar)  &
              & - ( - 1.0   - zC2/4.0  ) * pp(i,j,k+1,nstar)  &
              & - ( 1.0/2.0 + zC2/4.0  ) * pp(i,j,k  ,nstar)  &
              & - ( 1.0/3.0 - zC2/12.0 ) * pp(i,j,k-1,nstar)
        END IF
        z_lgs_rhs(i,j,k) = pp(i,j,k,nstar)                  &
          &   - z_Cour(i,j,k) * zh1 * z_beta_m 

      ENDDO
    ENDDO
  ENDDO

  ! layer k=ke-1
  DO j = jstart, jend
    DO i = istart, iend
      z_lgs_rhs(i,j,ke-1) = pp(i,j,ke-1,nstar)          &
        &   - z_Cour(i,j,ke-1) * ( pp(i,j,ke,nstar)-pp(i,j,ke-2,nstar) ) * z_beta_m_half
    ENDDO
  ENDDO

  ! Bottom layer k=ke
  DO j = jstart, jend
    DO i = istart, iend
      IF ( z_Cour(i,j,k) > 0 ) THEN
        z_lgs_rhs(i,j,ke) = pp(i,j,ke,nstar)              &
          &  - z_Cour(i,j,ke) * z_beta_m * ( pp(i,j,ke,nstar)-pp(i,j,ke-1,nstar) )
      ELSE
        ! use explicit upwind:
        z_lgs_rhs(i,j,ke) = pp(i,j,ke,nstar)              &
          &  - z_Cour(i,j,ke) * 1.0_ireals * ( pp(i,j,ke,nstar)-pp(i,j,ke-1,nstar) )
      END IF
    ENDDO
  ENDDO

  IF ( flag_vector_version ) THEN 
    CALL solve_5banddiag_vec( z_lgs, z_lgs_rhs, zh, &
        ie, je, istart, iend, jstart, jend, ke ) 
    ! attention: z_lgs has now changed

  ELSE
    CALL solve_5banddiag( z_lgs, z_lgs_rhs, zh, &
        ie, je, istart, iend, jstart, jend, ke ) 
  END IF

  DO k = 1, ke
    DO j = jstart, jend
      DO i = istart, iend
        ppadvt(i,j,k) = ( zh(i,j,k) - pp(i,j,k, nstar) ) * z_dt_recip  
      ENDDO
    ENDDO
  ENDDO

  !----------------------------------------------------------------------------
  ! Section 3: Setup of tridiagonal matrix systems resulting from the implicit
  !            numerical formulation of advection.
  !            -  variables at w-positions (w)
  !----------------------------------------------------------------------------


  ! Band-diagonal-matrix

  ! Top layer k=1 (one-sided diff.)  
  DO j = jstart, jend
    DO i = istart, iend

      z_Cour_w(i,j,1) = sqrtg_r_w(i,j,1) * wcon(i,j,1) * dt_adv 

      IF ( z_Cour(i,j,1) < 0.0_ireals ) THEN
        z_lgs_w(i,j,1,1) = 0.0_ireals
        z_lgs_w(i,j,1,2) = 0.0_ireals
        z_lgs_w(i,j,1,3) = 1.0_ireals - z_Cour_w(i,j,1) * z_beta_p
        z_lgs_w(i,j,1,4) =            + z_Cour_w(i,j,1) * z_beta_p
        z_lgs_w(i,j,1,5) = 0.0_ireals
      ELSE
        ! use explicit upwind
        z_lgs_w(i,j,1,1) = 0.0_ireals
        z_lgs_w(i,j,1,2) = 0.0_ireals
        z_lgs_w(i,j,1,3) = 1.0_ireals
        z_lgs_w(i,j,1,4) = 0.0_ireals
        z_lgs_w(i,j,1,5) = 0.0_ireals
      END IF

    ENDDO
  ENDDO

  ! 2. Layer k=2 (centered diff. 2. order)
  DO j = jstart, jend
    DO i = istart, iend

      z_Cour_w(i,j,2) = sqrtg_r_w(i,j,2) * wcon(i,j,2) * dt_adv 

      z_lgs_w(i,j,2,1) = 0.0_ireals
      z_lgs_w(i,j,2,2) = - z_Cour_w(i,j,2) * z_beta_p_half
      z_lgs_w(i,j,2,3) = 1.0_ireals
      z_lgs_w(i,j,2,4) = - z_lgs_w(i,j,2,2)
      z_lgs_w(i,j,2,5) = 0.0_ireals

    ENDDO
  ENDDO

  ! The layers from k=3 to k=ke-1
  DO k = 3, ke-1
    DO j = jstart, jend
      DO i = istart, iend

        z_Cour_w(i,j,k) = sqrtg_r_w(i,j,k) * wcon(i,j,k) * dt_adv 
        zC2 = z_Cour_w(i,j,k) * z_Cour_w(i,j,k)
        zh2 = z_Cour_w(i,j,k) * z_beta_p

        IF ( z_Cour_w(i,j,k) > 0 ) THEN
          z_lgs_w(i,j,k,1) = zh2 * ( 1.0/6.0 + zC2/12.0 )
          z_lgs_w(i,j,k,2) = zh2 * ( - 1.0   - zC2/4.0  )
          z_lgs_w(i,j,k,3) = zh2 * ( 1.0/2.0 + zC2/4.0  ) + 1.0
          z_lgs_w(i,j,k,4) = zh2 * ( 1.0/3.0 - zC2/12.0 )
          z_lgs_w(i,j,k,5) = 0.0_ireals
        ELSE 
          z_lgs_w(i,j,k,1) = 0.0_ireals
          z_lgs_w(i,j,k,2) = - zh2 * ( 1.0/3.0 - zC2/12.0 )
          z_lgs_w(i,j,k,3) = - zh2 * ( 1.0/2.0 + zC2/4.0  ) + 1.0
          z_lgs_w(i,j,k,4) = - zh2 * ( - 1.0   - zC2/4.0  )
          z_lgs_w(i,j,k,5) = - zh2 * ( 1.0/6.0 + zC2/12.0 ) 
        END IF

      ENDDO
    ENDDO
  ENDDO

  ! Layer k=ke (centered diff. 2. order)
  DO j = jstart, jend
    DO i = istart, iend

      z_Cour_w(i,j,ke) = sqrtg_r_w(i,j,ke) * wcon(i,j,ke) * dt_adv 

      z_lgs_w(i,j,ke,1) = 0.0_ireals
      z_lgs_w(i,j,ke,2) = - z_Cour_w(i,j,ke) * z_beta_p_half
      z_lgs_w(i,j,ke,3) = 1.0_ireals
      z_lgs_w(i,j,ke,4) = - z_lgs_w(i,j,ke,2)
      z_lgs_w(i,j,ke,5) = 0.0_ireals

    ENDDO
  ENDDO

  ! The bottom layer k=ke+1
  DO j = jstart, jend
    DO i = istart, iend

      z_Cour_w(i,j,ke1) = sqrtg_r_w(i,j,ke1) * wcon(i,j,ke1) * dt_adv  

      IF ( z_Cour(i,j,k) > 0 ) THEN
        z_lgs_w(i,j,ke1,1) = 0.0_ireals
        z_lgs_w(i,j,ke1,2) =            - z_Cour_w(i,j,ke1) * z_beta_p
        z_lgs_w(i,j,ke1,3) = 1.0_ireals + z_Cour_w(i,j,ke1) * z_beta_p
        z_lgs_w(i,j,ke1,4) = 0.0_ireals
        z_lgs_w(i,j,ke1,5) = 0.0_ireals
      ELSE
        z_lgs_w(i,j,ke1,1) = 0.0_ireals
        z_lgs_w(i,j,ke1,2) = 0.0_ireals
        z_lgs_w(i,j,ke1,3) = 1.0_ireals
        z_lgs_w(i,j,ke1,4) = 0.0_ireals
        z_lgs_w(i,j,ke1,5) = 0.0_ireals
      END IF
    ENDDO
  ENDDO

  ! --- right hand side for w ---

  ! Top layer (k=1)
  DO j = jstart, jend
    DO i = istart, iend
      IF ( z_Cour(i,j,1) < 0.0_ireals ) THEN
        z_lgs_rhs_w(i,j,1) = w(i,j,1,nstar)                                   &
          &  - z_Cour_w(i,j,1) * z_beta_m * ( w(i,j,2,nstar)-w(i,j,1,nstar) )
      ELSE
        ! use explicit upwind:
        z_lgs_rhs_w(i,j,1) = w(i,j,1,nstar)                                   &
          &  - z_Cour_w(i,j,1) * 1.0_ireals * ( w(i,j,2,nstar)-w(i,j,1,nstar) )
      END IF
    ENDDO
  ENDDO

  ! 2. layer k=2
  DO j = jstart, jend
    DO i = istart, iend
      z_lgs_rhs_w(i,j,2) = w(i,j,2,nstar)                                          &
        &  - z_Cour_w(i,j,2) * ( w(i,j,3,nstar) - w(i,j,1,nstar) ) * z_beta_m_half

    ENDDO
  ENDDO

  ! The layers from k=3 to k=ke-1
  DO k = 3, ke-1
    DO j = jstart, jend
      DO i = istart, iend

        zC2 = z_Cour_w(i,j,k) * z_Cour_w(i,j,k) 

        ! store advection-operator in help-variable zh1:
        IF ( z_Cour_w(i,j,k) > 0 ) THEN
          zh1 =   ( 1.0/6.0 + zC2/12.0 ) * w(i,j,k-2,nstar)  &
              & + ( - 1.0   - zC2/4.0  ) * w(i,j,k-1,nstar)  &
              & + ( 1.0/2.0 + zC2/4.0  ) * w(i,j,k  ,nstar)  &
              & + ( 1.0/3.0 - zC2/12.0 ) * w(i,j,k+1,nstar)
        ELSE 
          zh1 = - ( 1.0/6.0 + zC2/12.0 ) * w(i,j,k+2,nstar)  &
              & - ( - 1.0   - zC2/4.0  ) * w(i,j,k+1,nstar)  &
              & - ( 1.0/2.0 + zC2/4.0  ) * w(i,j,k  ,nstar)  &
              & - ( 1.0/3.0 - zC2/12.0 ) * w(i,j,k-1,nstar)
        END IF
        z_lgs_rhs_w(i,j,k) = w(i,j,k,nstar)                  &
          &   - z_Cour_w(i,j,k) * zh1 * z_beta_m

      ENDDO
    ENDDO
  ENDDO

  ! layer k=ke
  DO j = jstart, jend
    DO i = istart, iend
      z_lgs_rhs_w(i,j,ke) = w(i,j,ke,nstar)                                              &
        &   - z_Cour_w(i,j,ke) * ( w(i,j,ke+1,nstar)-w(i,j,ke-1,nstar) ) * z_beta_m_half
    ENDDO
  ENDDO

  ! Bottom layer k=ke+1
  DO j = jstart, jend
    DO i = istart, iend
      IF ( z_Cour(i,j,k) > 0 ) THEN
        z_lgs_rhs_w(i,j,ke1) = w(i,j,ke1,nstar)              &
          &  - z_Cour_w(i,j,ke1) * z_beta_m * ( w(i,j,ke1,nstar)-w(i,j,ke,nstar) )
      ELSE
        ! use explicit upwind:
        z_lgs_rhs_w(i,j,ke1) = w(i,j,ke1,nstar)              &
          &  - z_Cour_w(i,j,ke1) * 1.0_ireals * ( w(i,j,ke1,nstar)-w(i,j,ke,nstar) )
      END IF
    ENDDO
  ENDDO

  IF ( flag_vector_version ) THEN 
    CALL solve_5banddiag_vec( z_lgs_w, z_lgs_rhs_w, zh_w, &
        ie, je, istart, iend, jstart, jend, ke1 ) 
    ! attention: z_lgs_w has now changed
  ELSE 
    CALL solve_5banddiag( z_lgs_w, z_lgs_rhs_w, zh_w, &
        ie, je, istart, iend, jstart, jend, ke1 ) 
  END IF

  DO k = 1, ke1
    DO j = jstart, jend
      DO i = istart, iend
        wadvt(i,j,k) = ( zh_w(i,j,k) - w(i,j,k, nstar) ) * z_dt_recip  
      ENDDO
    ENDDO
  ENDDO

  !----------------------------------------------------------------------------
  ! Section 4: Setup of tridiagonal matrix systems resulting from the implicit
  !            numerical formulation of advection.
  !            - variables at u-positions (u)
  !----------------------------------------------------------------------------

  ! Band-diagonal-matrix

  ! Top layer k=1 (one-sided diff.)  
  DO j = jstartu, jendu
    DO i = istartu, iendu

      z_Cour(i,j,1) = sqrtg_r_u(i,j,1) *                         &
        &             ( wcon(i,  j,1) + wcon(i,  j,2)            &
        &             + wcon(i+1,j,1) + wcon(i+1,j,2) ) * z_dt_quart 

      IF ( z_Cour(i,j,1) < 0.0_ireals ) THEN
        z_lgs(i,j,1,1) = 0.0_ireals
        z_lgs(i,j,1,2) = 0.0_ireals
        z_lgs(i,j,1,3) = 1.0_ireals - z_Cour(i,j,1) * z_beta_p
        z_lgs(i,j,1,4) =            + z_Cour(i,j,1) * z_beta_p
        z_lgs(i,j,1,5) = 0.0_ireals
      ELSE
        ! use explicit upwind:
        z_lgs(i,j,1,1) = 0.0_ireals
        z_lgs(i,j,1,2) = 0.0_ireals
        z_lgs(i,j,1,3) = 1.0_ireals
        z_lgs(i,j,1,4) = 0.0_ireals
        z_lgs(i,j,1,5) = 0.0_ireals
      END IF
    ENDDO
  ENDDO

  ! 2. Layer k=2 (centered diff. 2. order)
  DO j = jstartu, jendu
    DO i = istartu, iendu

      z_Cour(i,j,2) =  sqrtg_r_u(i,j,2) *                        &
        &             ( wcon(i,  j,2) + wcon(i,  j,3)            &
        &             + wcon(i+1,j,2) + wcon(i+1,j,3) ) * z_dt_quart 

      z_lgs(i,j,2,1) = 0.0_ireals
      z_lgs(i,j,2,2) = - z_Cour(i,j,2) * z_beta_p_half
      z_lgs(i,j,2,3) = 1.0_ireals
      z_lgs(i,j,2,4) = - z_lgs(i,j,2,2)
      z_lgs(i,j,2,5) = 0.0_ireals

    ENDDO
  ENDDO

  ! The layers from k=3 to k=ke-2
  DO k = 3, ke-2
    DO j = jstartu, jendu
      DO i = istartu, iendu

        z_Cour(i,j,k) =  sqrtg_r_u(i,j,k) *                         &
          &             ( wcon(i,  j,k) + wcon(i,  j,k+1)           &
          &             + wcon(i+1,j,k) + wcon(i+1,j,k+1) ) * z_dt_quart 

        zC2 = z_Cour(i,j,k) * z_Cour(i,j,k)
        zh2 = z_Cour(i,j,k) * z_beta_p

        IF ( z_Cour(i,j,k) > 0 ) THEN
          z_lgs(i,j,k,1) = zh2 * ( 1.0/6.0 + zC2/12.0 )
          z_lgs(i,j,k,2) = zh2 * ( - 1.0   - zC2/4.0  )
          z_lgs(i,j,k,3) = zh2 * ( 1.0/2.0 + zC2/4.0  ) + 1.0
          z_lgs(i,j,k,4) = zh2 * ( 1.0/3.0 - zC2/12.0 )
          z_lgs(i,j,k,5) = 0.0_ireals
        ELSE 
          z_lgs(i,j,k,1) = 0.0_ireals
          z_lgs(i,j,k,2) = - zh2 * ( 1.0/3.0 - zC2/12.0 )
          z_lgs(i,j,k,3) = - zh2 * ( 1.0/2.0 + zC2/4.0  ) + 1.0
          z_lgs(i,j,k,4) = - zh2 * ( - 1.0   - zC2/4.0  )
          z_lgs(i,j,k,5) = - zh2 * ( 1.0/6.0 + zC2/12.0 )
        END IF

      ENDDO
    ENDDO
  ENDDO

  ! Layer k=ke-1 (centered diff. 2. order)
  DO j = jstartu, jendu
    DO i = istartu, iendu

      z_Cour(i,j,ke-1) =  sqrtg_r_u(i,j,ke-1) *                        &
        &                ( wcon(i,  j,ke-1)+ wcon(i,  j,ke)            &
        &                + wcon(i+1,j,ke-1)+ wcon(i+1,j,ke) ) * z_dt_quart  

      z_lgs(i,j,ke-1,1) = 0.0_ireals
      z_lgs(i,j,ke-1,2) = - z_Cour(i,j,ke-1) * z_beta_p_half
      z_lgs(i,j,ke-1,3) = 1.0_ireals
      z_lgs(i,j,ke-1,4) = - z_lgs(i,j,ke-1,2)
      z_lgs(i,j,ke-1,5) = 0.0_ireals

    ENDDO
  ENDDO

  ! The bottom layer k=ke
  DO j = jstartu, jendu
    DO i = istartu, iendu

      z_Cour(i,j,ke) =  sqrtg_r_u(i,j,ke) *                          &
        &              ( wcon(i,  j,ke)+ wcon(i,  j,ke+1)            &
        &              + wcon(i+1,j,ke)+ wcon(i+1,j,ke+1) ) * z_dt_quart  

      IF ( z_Cour(i,j,k) > 0 ) THEN
        z_lgs(i,j,ke,1) = 0.0_ireals
        z_lgs(i,j,ke,2) =            - z_Cour(i,j,ke) * z_beta_p
        z_lgs(i,j,ke,3) = 1.0_ireals + z_Cour(i,j,ke) * z_beta_p
        z_lgs(i,j,ke,4) = 0.0_ireals
        z_lgs(i,j,ke,5) = 0.0_ireals
      ELSE
        ! use explicit upwind:
        z_lgs(i,j,ke,1) = 0.0_ireals
        z_lgs(i,j,ke,2) = 0.0_ireals
        z_lgs(i,j,ke,3) = 1.0_ireals
        z_lgs(i,j,ke,4) = 0.0_ireals
        z_lgs(i,j,ke,5) = 0.0_ireals
      END IF
    ENDDO
  ENDDO

  ! --- right hand side for u ---

  ! Top layer (k=1)
  DO j = jstartu, jendu
    DO i = istartu, iendu
      IF ( z_Cour(i,j,1) < 0.0_ireals ) THEN
        z_lgs_rhs(i,j,1) = u(i,j,1,nstar)                                   &
          &  - z_Cour(i,j,1) * z_beta_m * ( u(i,j,2,nstar)-u(i,j,1,nstar) ) 
      ELSE
        ! use explicit upwind:
        z_lgs_rhs(i,j,1) = u(i,j,1,nstar)                                   &
          &  - z_Cour(i,j,1) * 1.0_ireals * ( u(i,j,2,nstar)-u(i,j,1,nstar) ) 
      END IF
    ENDDO
  ENDDO

  ! 2. layer k=2
  DO j = jstartu, jendu
    DO i = istartu, iendu
      z_lgs_rhs(i,j,2) = u(i,j,2,nstar)                                          &
        &  - z_Cour(i,j,2) * ( u(i,j,3,nstar) - u(i,j,1,nstar) ) * z_beta_m_half 

    ENDDO
  ENDDO

  ! The layers from k=3 to k=ke-2
  DO k = 3, ke-2
    DO j = jstartu, jendu
      DO i = istartu, iendu

        zC2 = z_Cour(i,j,k) * z_Cour(i,j,k) 

        ! store advection-operator in help-variable zh1:
        IF ( z_Cour(i,j,k) > 0 ) THEN
          zh1 =   ( 1.0/6.0 + zC2/12.0 ) * u(i,j,k-2,nstar)  &
              & + ( - 1.0   - zC2/4.0  ) * u(i,j,k-1,nstar)  &
              & + ( 1.0/2.0 + zC2/4.0  ) * u(i,j,k  ,nstar)  &
              & + ( 1.0/3.0 - zC2/12.0 ) * u(i,j,k+1,nstar)
        ELSE 
          zh1 = - ( 1.0/6.0 + zC2/12.0 ) * u(i,j,k+2,nstar)  &
              & - ( - 1.0   - zC2/4.0  ) * u(i,j,k+1,nstar)  &
              & - ( 1.0/2.0 + zC2/4.0  ) * u(i,j,k  ,nstar)  &
              & - ( 1.0/3.0 - zC2/12.0 ) * u(i,j,k-1,nstar)
        END IF
        z_lgs_rhs(i,j,k) = u(i,j,k,nstar)                  &
          &   - z_Cour(i,j,k) * zh1 * z_beta_m 

      ENDDO
    ENDDO
  ENDDO

  ! layer k=ke-1
  DO j = jstartu, jendu
    DO i = istartu, iendu
      z_lgs_rhs(i,j,ke-1) = u(i,j,ke-1,nstar)                                          &
        &   - z_Cour(i,j,ke-1) * ( u(i,j,ke,nstar)-u(i,j,ke-2,nstar) ) * z_beta_m_half
    ENDDO
  ENDDO

  ! Bottom layer k=ke
  DO j = jstartu, jendu
    DO i = istartu, iendu
      IF ( z_Cour(i,j,k) > 0 ) THEN
        z_lgs_rhs(i,j,ke) = u(i,j,ke,nstar)              &
          &  - z_Cour(i,j,ke) * z_beta_m * ( u(i,j,ke,nstar)-u(i,j,ke-1,nstar) )
      ELSE
        ! use explicit upwind:
        z_lgs_rhs(i,j,ke) = u(i,j,ke,nstar)              &
          &  - z_Cour(i,j,ke) * 1.0_ireals * ( u(i,j,ke,nstar)-u(i,j,ke-1,nstar) )
      END IF
    ENDDO
  ENDDO

  IF ( flag_vector_version ) THEN 
    CALL solve_5banddiag_vec( z_lgs, z_lgs_rhs, zh,  &
       ie, je, istartu, iendu, jstartu, jendu, ke ) 
    ! attention: z_lgs has now changed
  ELSE
    CALL solve_5banddiag( z_lgs, z_lgs_rhs, zh,  &
       ie, je, istartu, iendu, jstartu, jendu, ke ) 
  END IF

  DO k = 1, ke
    DO j = jstartu, jendu
      DO i = istartu, iendu
        uadvt(i,j,k) = ( zh(i,j,k) - u(i,j,k, nstar) ) * z_dt_recip  
      ENDDO
    ENDDO
  ENDDO

  !----------------------------------------------------------------------------
  ! Section 5: Setup of tridiagonal matrix systems resulting from the implicit
  !            numerical formulation of advection.
  !            - variables at v-positions (v)
  !----------------------------------------------------------------------------

  ! Band-diagonal-matrix

  ! Top layer k=1 (one-sided diff.)  
  DO j = jstartv, jendv
    DO i = istartv, iendv

      z_Cour(i,j,1) =  sqrtg_r_v(i,j,1) *                         &
        &             ( wcon(i,j,  1) + wcon(i,j,  2)             &
        &             + wcon(i,j+1,1) + wcon(i,j+1,2) ) * z_dt_quart

      IF ( z_Cour(i,j,1) < 0.0_ireals ) THEN
        z_lgs(i,j,1,1) = 0.0_ireals
        z_lgs(i,j,1,2) = 0.0_ireals
        z_lgs(i,j,1,3) = 1.0_ireals - z_Cour(i,j,1) * z_beta_p
        z_lgs(i,j,1,4) =            + z_Cour(i,j,1) * z_beta_p
        z_lgs(i,j,1,5) = 0.0_ireals
      ELSE
        ! use explicit upwind:
        z_lgs(i,j,1,1) = 0.0_ireals
        z_lgs(i,j,1,2) = 0.0_ireals
        z_lgs(i,j,1,3) = 1.0_ireals
        z_lgs(i,j,1,4) = 0.0_ireals
        z_lgs(i,j,1,5) = 0.0_ireals
      END IF

    ENDDO
  ENDDO

  ! 2. Layer k=2 (centered diff. 2. order)
  DO j = jstartv, jendv
    DO i = istartv, iendv

      z_Cour(i,j,2) =  sqrtg_r_v(i,j,2) *                         &
        &             ( wcon(i,j,  2) + wcon(i,j,  3)             &
        &             + wcon(i,j+1,2) + wcon(i,j+1,3) ) * z_dt_quart

      z_lgs(i,j,2,1) = 0.0_ireals
      z_lgs(i,j,2,2) = - z_Cour(i,j,2) * z_beta_p_half
      z_lgs(i,j,2,3) = 1.0_ireals
      z_lgs(i,j,2,4) = - z_lgs(i,j,2,2)
      z_lgs(i,j,2,5) = 0.0_ireals

    ENDDO
  ENDDO

  ! The layers from k=3 to k=ke-2
  DO k = 3, ke-2
    DO j = jstartv, jendv
      DO i = istartv, iendv

        z_Cour(i,j,k) =  sqrtg_r_v(i,j,k) *                         &
          &             ( wcon(i,j,  k) + wcon(i,j,  k+1)           &
          &             + wcon(i,j+1,k) + wcon(i,j+1,k+1) ) * z_dt_quart

        zC2 = z_Cour(i,j,k) * z_Cour(i,j,k)
        zh2 = z_Cour(i,j,k) * z_beta_p

        IF ( z_Cour(i,j,k) > 0 ) THEN
          z_lgs(i,j,k,1) = zh2 * ( 1.0/6.0 + zC2/12.0 )
          z_lgs(i,j,k,2) = zh2 * ( - 1.0   - zC2/4.0  )
          z_lgs(i,j,k,3) = zh2 * ( 1.0/2.0 + zC2/4.0  ) + 1.0
          z_lgs(i,j,k,4) = zh2 * ( 1.0/3.0 - zC2/12.0 )
          z_lgs(i,j,k,5) = 0.0_ireals
        ELSE 
          z_lgs(i,j,k,1) = 0.0_ireals
          z_lgs(i,j,k,2) = - zh2 * ( 1.0/3.0 - zC2/12.0 )
          z_lgs(i,j,k,3) = - zh2 * ( 1.0/2.0 + zC2/4.0  ) + 1.0
          z_lgs(i,j,k,4) = - zh2 * ( - 1.0   - zC2/4.0  )
          z_lgs(i,j,k,5) = - zh2 * ( 1.0/6.0 + zC2/12.0 ) 
        END IF

      ENDDO
    ENDDO
  ENDDO

  ! Layer k=ke-1 (centered diff. 2. order)
  DO j = jstartv, jendv
    DO i = istartv, iendv
      z_Cour(i,j,ke-1) =  sqrtg_r_v(i,j,ke-1) *                      &
        &                ( wcon(i,j,  ke-1)+ wcon(i,j,  ke)          &
        &                + wcon(i,j+1,ke-1)+ wcon(i,j+1,ke) ) * z_dt_quart
 
      z_lgs(i,j,ke-1,1) = 0.0_ireals
      z_lgs(i,j,ke-1,2) = - z_Cour(i,j,ke-1) * z_beta_p_half
      z_lgs(i,j,ke-1,3) = 1.0_ireals
      z_lgs(i,j,ke-1,4) = - z_lgs(i,j,ke-1,2)
      z_lgs(i,j,ke-1,5) = 0.0_ireals

    ENDDO
  ENDDO

  ! The bottom layer k=ke
  DO j = jstartv, jendv
    DO i = istartv, iendv

      z_Cour(i,j,ke) =  sqrtg_r_v(i,j,ke) *                        &
        &              ( wcon(i,j,  ke)+ wcon(i,j,  ke+1)          &
        &              + wcon(i,j+1,ke)+ wcon(i,j+1,ke+1) ) * z_dt_quart
 
      IF ( z_Cour(i,j,k) > 0 ) THEN
        z_lgs(i,j,ke,1) = 0.0_ireals
        z_lgs(i,j,ke,2) =            - z_Cour(i,j,ke) * z_beta_p
        z_lgs(i,j,ke,3) = 1.0_ireals + z_Cour(i,j,ke) * z_beta_p
        z_lgs(i,j,ke,4) = 0.0_ireals
        z_lgs(i,j,ke,5) = 0.0_ireals
      ELSE
        ! use explicit upwind:
        z_lgs(i,j,ke,1) = 0.0_ireals
        z_lgs(i,j,ke,2) = 0.0_ireals
        z_lgs(i,j,ke,3) = 1.0_ireals
        z_lgs(i,j,ke,4) = 0.0_ireals
        z_lgs(i,j,ke,5) = 0.0_ireals
      END IF
    ENDDO
  ENDDO

  ! --- right hand side for v ---

  ! Top layer (k=1)      
  DO j = jstartv, jendv
    DO i = istartv, iendv
      IF ( z_Cour(i,j,1) < 0.0_ireals ) THEN
        z_lgs_rhs(i,j,1) = v(i,j,1,nstar)                                   &
          &  - z_Cour(i,j,1) * z_beta_m * ( v(i,j,2,nstar)-v(i,j,1,nstar) ) 
      ELSE
        ! use explicit upwind:
        z_lgs_rhs(i,j,1) = v(i,j,1,nstar)                                   &
          &  - z_Cour(i,j,1) * 1.0_ireals * ( v(i,j,2,nstar)-v(i,j,1,nstar) ) 
      END IF
    ENDDO
  ENDDO

  ! 2. layer k=2
  DO j = jstartv, jendv
    DO i = istartv, iendv
      z_lgs_rhs(i,j,2) = v(i,j,2,nstar)                                          &
        &  - z_Cour(i,j,2) * ( v(i,j,3,nstar) - v(i,j,1,nstar) ) * z_beta_m_half 

    ENDDO
  ENDDO

  ! The layers from k=3 to k=ke-2
  DO k = 3, ke-2
    DO j = jstartv, jendv
      DO i = istartv, iendv

        zC2 = z_Cour(i,j,k) * z_Cour(i,j,k) 

        ! store advection-operator in help-variable zh1:
        IF ( z_Cour(i,j,k) > 0 ) THEN
          zh1 =   ( 1.0/6.0 + zC2/12.0 ) * v(i,j,k-2,nstar)  &
              & + ( - 1.0   - zC2/4.0  ) * v(i,j,k-1,nstar)  &
              & + ( 1.0/2.0 + zC2/4.0  ) * v(i,j,k  ,nstar)  &
              & + ( 1.0/3.0 - zC2/12.0 ) * v(i,j,k+1,nstar)
        ELSE 
          zh1 = - ( 1.0/6.0 + zC2/12.0 ) * v(i,j,k+2,nstar)  &
              & - ( - 1.0   - zC2/4.0  ) * v(i,j,k+1,nstar)  &
              & - ( 1.0/2.0 + zC2/4.0  ) * v(i,j,k  ,nstar)  &
              & - ( 1.0/3.0 - zC2/12.0 ) * v(i,j,k-1,nstar)
        END IF
        z_lgs_rhs(i,j,k) = v(i,j,k,nstar)                  &
          &   - z_Cour(i,j,k) * zh1 * z_beta_m 

      ENDDO
    ENDDO
  ENDDO

  ! layer k=ke-1
  DO j = jstartv, jendv
    DO i = istartv, iendv
      z_lgs_rhs(i,j,ke-1) = v(i,j,ke-1,nstar)                                          &
        &   - z_Cour(i,j,ke-1) * ( v(i,j,ke,nstar)-v(i,j,ke-2,nstar) ) * z_beta_m_half 
    ENDDO
  ENDDO

  ! Bottom layer k=ke
  DO j = jstartv, jendv
    DO i = istartv, iendv
      IF ( z_Cour(i,j,k) > 0 ) THEN
        z_lgs_rhs(i,j,ke) = v(i,j,ke,nstar)              &
          &  - z_Cour(i,j,ke) * z_beta_m * ( v(i,j,ke,nstar)-v(i,j,ke-1,nstar) ) 
      ELSE
        ! use explicit upwind:
        z_lgs_rhs(i,j,ke) = v(i,j,ke,nstar)              &
          &  - z_Cour(i,j,ke) * 1.0_ireals * ( v(i,j,ke,nstar)-v(i,j,ke-1,nstar) ) 
      END IF
    ENDDO
  ENDDO

  IF ( flag_vector_version ) THEN 
    CALL solve_5banddiag_vec( z_lgs, z_lgs_rhs, zh,  &
        ie, je, istartv, iendv, jstartv, jendv, ke ) 
    ! attention: z_lgs has now changed
  ELSE
    CALL solve_5banddiag( z_lgs, z_lgs_rhs, zh,  &
        ie, je, istartv, iendv, jstartv, jendv, ke ) 
  END IF

  DO k = 1, ke
    DO j = jstartv, jendv
      DO i = istartv, iendv
        vadvt(i,j,k) = ( zh(i,j,k) - v(i,j,k, nstar) ) * z_dt_recip  
      ENDDO
    ENDDO
  ENDDO

  DEALLOCATE( z_Cour, z_lgs, z_lgs_rhs, zh )

  DEALLOCATE( z_Cour_w, z_lgs_w, z_lgs_rhs_w, zh_w )

  !----------------------------------------------------------------------------
  ! End of subroutine complete_tend_uvwtpp_CN3Crow
  !----------------------------------------------------------------------------

END SUBROUTINE complete_tend_uvwtpp_CN3Crow


END MODULE src_slow_tendencies_rk
