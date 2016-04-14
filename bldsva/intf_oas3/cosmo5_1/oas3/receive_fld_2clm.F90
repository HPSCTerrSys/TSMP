SUBROUTINE receive_fld_2clm

!---------------------------------------------------------------------
! Description:
!  This routine receives fluxes from CLM3.5 and updates the mixing 
!  mixing coefficients in COSMO (Important surface layer physics involved)
!
! Current Code Owner: TR32, Z4: Prabhakar Shrestha
!    phone: 0228733453
!    email: pshrestha@uni-bonn.de
!
! History:
! Version    Date       Name
! ---------- ---------- ----
! 1.1.1        2011/11/28 Prabhakar Shrestha 
!   Modfied and Implemented in COSMO4.11, Initial release
! 1.2.1        2012/09/21 Markus Uebel 
!   Inclusion of CO2 fluxes (photosynthesis rate)
! 1.3.1        2013/01/31 Prabhakar Shrestha
!   Inclusion of aerodynamic resistance from CLM
!    ram1, rah1, raw1
!   Inclusion of surface temperature, humidity
!    t_sf, q_sf
!   This new option updates the surface temperature, surface humidity and
!   transfer coefficients directly without inversion. But now coupling frequency
!   between COSMO and CLM and timesteps should be the same.
! 1.3.2        2014/01/05 Prabhakar Shrestha
!   Inclusion of masked coupling between cosmo and clm (see zmask)
! 1.3.3        2015/08/26 Prabhakar Shrestha
!   Inclusion of full "frcv" including the outerbound halos
!
! @VERSION@    @DATE@     <Your name>
!  <Modification comments>         
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!==============================================================================

! Declarations:
!
! Modules used:

USE oas_cos_vardef
USE data_modelconfig, ONLY :  ke, dt, ie, je,                         &
                              istartpar, iendpar, jstartpar, jendpar, &
                              istart, iend, jstart, jend,             &
                              ie_tot, je_tot, ke1, idt_qv

USE data_parameters,  ONLY :  wp, iintegers
USE data_parallel,    ONLY :                &
                            isubpos,        &        ! positions of the subdomains in the total domain
                            nboundlines,    &    
                            my_cart_id

USE data_runcontrol , ONLY :                &
          lphys,        & ! forecast with physical parametrizations
          llake,        & ! forecst with lake model FLake
          lgsp,         & ! forecast with grid-scale precipitation
          lprog_qi,     & ! if .TRUE., running with cloud ice
          ltur,         & ! forecast with turbulent diffusion
          lcpfluc,      & ! consideration of fluctuations of the heat capacity of air
          hincrad,      & ! increment for running the radiation in hours 
          itype_turb,   & ! type of turbulent diffusion parametrization
          imode_turb,   & ! mode of turbulent diffusion parametrization
          l2tls     ,   & ! forecast with 1-TL integration scheme !CPS added
          nold  ,       & ! corresponds to ntstep - 1
          nnow  ,       & ! corresponds to ntstep
          nnew  ,       & ! corresponds to ntstep + 1
          nstop ,       & ! last time step of the forecast period
          ntstep

USE data_constants  , ONLY :                &
          r_d,          & ! gas constant for dry air
          rdv,          & ! r_d / r_v
          o_m_rdv,      & ! 1 - r_d/r_v
          rvd_m_o,      & ! r_v/r_d - 1
          cp_d,         & ! specific heat of dry air at constant pressure
          rdocp,        & ! r_d / cp_d
          cpdr,         & ! 1 / cp_d
          gamma,        & ! 1 / (1 - rdocp)   ( = cp_d/cv_d)
          lh_v,         & ! latent heat of vapourization
          g,            & ! acceleration due to gravity
          gq,           & ! g*g
          gr,           & ! 1 / g
          r_earth,      & ! radius of the earth/acrlat
          sigma         ! Boltzmann-constant

USE data_soil ,      ONLY :   ctalb       ! thermal albedo (1-emissivity) ?!

#ifdef COUP_OAS_COS
USE data_fields     , ONLY :         &
                     tcw,            & ! turbulent transfer coefficient for moisture   ( -- )  !CPS
                     trad_clm          ! radiative temperature for CLM as bottom BC     ( K )  !CPS
#endif

USE data_fields     , ONLY :                &
                     llandmask  ,    & ! landpoint mask (fr_land>=0.5)
                     fr_land,        & ! land fraction
                     alb_rad     ,   & ! albedo of the ground                            --
                     tch         ,   & ! turbulent transfer coefficient for heat       ( -- )
                     tcm         ,   & ! turbulent transfer coefficient for momentum   ( -- )
                     tfh         ,   & ! factor of laminar transfer of scalars         ( -- )
                     a1t, a2t  ,     & ! implicit weight of vertical diffusion
                     hhl       ,     & ! geometrical height of half levels             ( m   )
                     hsurf     ,     & ! height of surface topography                  ( m   )
    
                 ! 1. constant fields for the reference atmosphere                     (unit)
                 ! -----------------------------------------------
                     p0         ,    & ! base state pressure                           (Pa)
    
                 ! 3. prognostic variables                                             (unit)
                 ! -----------------------
                     u          ,    & ! zonal wind speed                              ( m/s )
                     v          ,    & ! meridional wind speed                         ( m/s )
                     t          ,    & ! temperature                                   (  k  )
!CPS                     qv         ,    & ! specific water vapor content                  (kg/kg)
                     pp         ,    & ! deviation from the reference pressure         ( pa  )
                     ps        ,     & ! surface pressure                              ( pa  )
!MU (19.09.2012)
                     rho        ,    & ! total density of air                          (kg/m3)
!MU (19.09.2012)

                 !    fields for surface values and soil model variables               (unit )
                 ! -----------------------------------------------------
                     t_g       ,     & ! weighted surface temperature                  (  K  )
                     t_s       ,     & ! temperature of the ground surface             (  k  )
                     qv_s      ,     & ! specific water vapor content on the surface   (kg/kg)

                     qvsflx     ,    & ! surface flux of water vapour                  (1/m2s)
                     umfl_s    ,     & ! u-momentum flux (surface)                     ( N/m2)
                     vmfl_s    ,     & ! v-momentum flux (surface)                     ( N/m2)
                     shfl_s    ,     & ! sensible heat flux (surface)                  ( W/m2)
                     lhfl_s            ! latent heat flux (surface)                    ( W/m2)

USE src_tracer,         ONLY: trcr_get, trcr_errorstr            !CPS
USE environment,        ONLY: model_abort                        !CPS
!MU (12.04.13)
!CPS
!USE data_tracer     , ONLY :                &
!CPS                     molmass_co2,&       ! molar mass of CO2                           (g/mol)
!CPS                     co2fl_s,&           ! CO2 tendency due to land-atmosphere exchange( 1/s )
!CPS                     psn_tens,&          ! CO2 tendency due to photosynthesis          ( 1/s )
!CPS                     plres_tens          ! CO2 tendency due to plant respiration       ( 1/s )
!MU (12.04.13)
USE data_turbulence , ONLY :   &
                     vel_min             ! minimal velocity scale [m/s]     !CPS
USE netcdf

!==============================================================================

IMPLICIT NONE

!==============================================================================

! Local Variable

INTEGER, PARAMETER         ::   jps_taux   =  1    !  zonal wind stress
INTEGER, PARAMETER         ::   jps_tauy   =  2    !  meridional wind stress
INTEGER, PARAMETER         ::   jps_lat    =  3    !  total latent heat flux (W/m**2)
INTEGER, PARAMETER         ::   jps_sens   =  4    !  total sensible heat flux (W/m**2)
INTEGER, PARAMETER         ::   jps_ir     =  5    !  emitted infrared (longwave) radiation (W/m**2)
INTEGER, PARAMETER         ::   jps_albd   =  6    !  direct albedo
INTEGER, PARAMETER         ::   jps_albi   =  7    !  diffuse albedo
!MU (18.09.2012)
INTEGER, PARAMETER         ::   jps_co2fl  =  8    !  net CO2 flux (now only photosynthesis) (u mol CO2 m-2s-1)
!MU (18.09.2012)
INTEGER, PARAMETER         ::   jps_ram1   =  9    !  Aerodynamic Resistance (s/m)
INTEGER, PARAMETER         ::   jps_rah1   =  10    !  Aerodynamic Resistance (s/m)
INTEGER, PARAMETER         ::   jps_raw1   =  11    !  Aerodynamic Resistance (s/m)
INTEGER, PARAMETER         ::   jps_tsf1   =  12    !  Surface Temperature (K) 
INTEGER, PARAMETER         ::   jps_qsf1   =  13    !  Surface Humidity (kg/kg)
!MU (12.04.13)
INTEGER, PARAMETER         ::   jps_fpsn   =  14    !  photosynthesis rate (u mol CO2 m-2s-1)
INTEGER, PARAMETER         ::   jps_fplres =  15    !  plant respiration (u mol CO2 m-2s-1)
!MU (12.04.13) 


INTEGER(KIND=iintegers)    ::   isec, info         ! temporary integer
INTEGER(KIND=iintegers)    ::   jn,ctr             ! index

INTEGER(KIND=iintegers) , DIMENSION(krcv) ::   nrcvinfo           ! OASIS info argument

INTEGER                    :: g1,i,j,n, k          ! do loop indices
INTEGER                    :: ier                  ! return error code
INTEGER                    :: nx            ! time index
INTEGER                    :: im1, jm1      ! i-1, j-1
INTEGER                    :: ip1, jp1      ! i+1, j+1

REAL (KIND=wp)         ::              & 
         ztvb  , zvbke,                    & !
         za1t_surf, za2t_surf,             & !
         zfpi,zppa,zpp,znew,ztmcmq,        &
         dzke, tcm_epsi, coef_lim, tch_epsi  ! by Oliver Fuhrer

! node-domain arrays
REAL (KIND=wp)         ::              &
         za1t    ( ke1),                   & !
         za2t    ( ke1),                   & !
         zpia    (ie,je),                  & !
         zpianf  (ie,je),                  & !
         ze   (ie,je,ke),                  & !
         ztcm    (ie,je),                  & !
         ztch    (ie,je),                  & !
         ztcw    (ie,je)                     !CPS


REAL(KIND=wp) ::   ztmp1(ie,je),       &      ! temporary array
                       umass(ie,je), vmass(ie,je)

REAL(kind=wp), SAVE ::   tic, toc

!CPS qx variables moved from data_fields to tracers
REAL (KIND=wp),     POINTER :: &
  qv  (:,:,:)    => NULL()             ! QV at given time-level 
CHARACTER (LEN=255)        :: yzerrmsg
CHARACTER (LEN=25)         :: yzroutine
!CPS 

INTEGER :: il_var_id(15),  dimids(3), il_file_id, mype, ib, npes, ierror, status        !CPS: dimids 2 to 3, MU: il_var_id 13 to 15

INTEGER :: cplstep, cplstop   !CPS cpl step

!------------------------------------------------------------------------------
!- End of header
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!- Begin Subroutine receive_fld_2clm
!------------------------------------------------------------------------------

 ! Statement function zfpi for calculation of the exner function
 ! where the dummy argument zppa is pressure
 zfpi(zppa) = (1.E-5_wp*zppa)**rdocp

 !  Coupling only with PE with at least one land point
 IF ( lpe_cpl ) THEN                  !CPS always true so redundant with OASIS3 ?

   nrcvinfo = OASIS_idle
   ztmp1=0._wp

   isec = ntstep * dt

   cplfreq = NINT ( hincrad * 3600.0_wp)        !CPS
!  CPS new values to remove kinks for startup from midnight
!CPS   IF (isec <= cplfreq ) THEN
!CPS     frcv(:,:,jps_taux) =   -0.00009_wp                     !CPS bandaid
!CPS     frcv(:,:,jps_tauy) =   -0.00009_wp                      !CPS bandaid 
!CPS     frcv(:,:,jps_lat ) =   4._wp                      !CPS bandaid
!CPS     frcv(:,:,jps_sens) =   -11.0_wp                      !CPS bandaid 
!CPS     frcv(:,:,jps_ir  ) = 386.079_wp                    !CPS bandaid 
!CPS     frcv(:,:,jps_albd) = 1._wp                       !CPS bandaid
!CPS     frcv(:,:,jps_albi) = 1._wp                       !CPS bandaid
!MU (18.09.12)
!CPS     frcv(:,:,jps_co2fl)=   0._wp
!MU (18.09.12)
!CPS     frcv(:,:,jps_ram1 )=   1089.77_wp
!CPS     frcv(:,:,jps_rah1 )=   1089.77_wp
!CPS     frcv(:,:,jps_raw1 )=   1089.77_wp
!CPS     frcv(:,:,jps_tsf1 )=   287.589_wp
!CPS     frcv(:,:,jps_qsf1 )=   0.010_wp
!MU (12.04.13)
!CPS     frcv(:,:,jps_fpsn)=   0._wp
!CPS     frcv(:,:,jps_fplres)=   0._wp
!MU (12.04.13)
!CPS   END IF

   cplstop = INT((nstop*dt+cplfreq)/cplfreq) - 1 !CPS

   IF (MOD(isec,cplfreq) == 0) THEN   !CPS
      cplstep = isec/cplfreq + 1      !CPS local variable so gets reset each time 
   ENDIF                              !CPS

   !  Performance measurement
   IF ( IOASISDEBUGLVL == 2 .and. my_cart_id == 0 ) THEN
      tic=MPI_Wtime()
   ENDIF

   ! frcv modified only at coupling time step
   ! frcv keep values from previous coupling (array from module) at other time step
   DO jn = 1, krcv
     IF( srcv(jn)%laction )   CALL oas_cos_rcv( jn, isec, ztmp1(:,:), nrcvinfo(jn) )   !CPS fix
!CPS_removebandaid     IF( nrcvinfo(jn) == OASIS_Rcv .and. isec>= cplfreq ) THEN
     IF( nrcvinfo(jn) == OASIS_Rcv) THEN
!CPS       frcv(nldi:nlei, nldj:nlej, jn)=ztmp1(nldi:nlei, nldj:nlej) 
        frcv(:,:,jn) = ztmp1(:,:)
     ENDIF
   ENDDO

   !  Performance measurement
   IF ( IOASISDEBUGLVL == 2 ) THEN
     CALL MPI_Barrier(kl_comm, ierror)
     IF( nrcvinfo(1) == OASIS_Rcv .and. my_cart_id == 0 ) THEN
       WRITE(6, *)  'COSMO calc_time ', tic-toc
       toc=MPI_Wtime()
       WRITE(6, *)  'COSMO After rcv ', toc
       WRITE(6, *)  'COSMO rcv_time ', toc-tic
       CALL flush(6)
     ENDIF
   ENDIF                 

 ENDIF                   !lpe_cpl


 ! CPS Write coupling fields for debug at coupling step only 
 IF ( IOASISDEBUGLVL == 1 .and. MOD(isec,cplfreq) == 0) THEN                 

   CALL MPI_Comm_rank(kl_comm, mype, ier)
   CALL MPI_Comm_size(kl_comm, npes, ier)

   IF (mype == 0 .and. isec == 0 ) then   !CPS open debug file at beginning only
     status = nf90_create("debugrcv_cos_clm.nc", NF90_CLOBBER, il_file_id)
     status = nf90_def_dim(il_file_id, "longitude", ie_tot - 2 * nboundlines, dimids(1))
     status = nf90_def_dim(il_file_id, "latitude", je_tot - 2 * nboundlines, dimids(2))
     status = nf90_def_dim(il_file_id, "time", cplstop, dimids(3))  !CPS Added time dimension
     ctr = 0
     DO jn = 1, krcv
       IF (srcv(jn)%laction) THEN
          ctr = ctr+1
          status = nf90_def_var(il_file_id, srcv(jn)%clname, NF90_DOUBLE, dimids, il_var_id(ctr))
       ENDIF
     ENDDO
     status = nf90_enddef(il_file_id)
     status = nf90_close(il_file_id)
   ENDIF

   DO ib = 0, npes - 1

     CALL MPI_Barrier(kl_comm, ierror)
  
     IF (mype == ib .AND. lpe_cpl ) THEN
       status = nf90_open("debugrcv_cos_clm.nc", NF90_WRITE, il_file_id)
       ctr = 0
       DO jn = 1, krcv
       IF (srcv(jn)%laction) THEN
          ctr = ctr+1
          status = nf90_inq_varid(il_file_id, srcv(jn)%clname , il_var_id(ctr)) 
          status = nf90_put_var(il_file_id, il_var_id(ctr),  frcv(nldi:nlei, nldj:nlej, jn), &
                                        start = (/ isubpos(my_cart_id,1)-nboundlines, &
                                                   isubpos(my_cart_id,2)-nboundlines, cplstep /), &
                                        count = (/ jih, jjh, 1 /) )
       ENDIF
       ENDDO
       status = nf90_close(il_file_id)

     ENDIF
    
   ENDDO

   ! End write coupling fields
 ENDIF

 CALL MPI_Barrier(kl_comm, ierror)


 nx = nnow          !CPS timelevel for coupling

 IF ( lpe_cpl ) THEN        

 ! Land points updated at each time step, frcv(:,:) updated only at coupling time step

   ! CESM / CLM3.5 fluxes (+up) ,taux and tauy are always -ve, COSMO fluxes (+down) 
   IF( srcv(jps_taux)%laction .and. nrcvinfo(jps_taux) == OASIS_Rcv ) frcv(:,:,jps_taux) = -frcv(:,:,jps_taux)
   IF( srcv(jps_tauy)%laction .and. nrcvinfo(jps_tauy) == OASIS_Rcv ) frcv(:,:,jps_tauy) = -frcv(:,:,jps_tauy)

   IF( srcv(jps_lat)%laction  .and. nrcvinfo(jps_lat)  == OASIS_Rcv ) frcv(:,:,jps_lat)  = -frcv(:,:,jps_lat)
   IF( srcv(jps_sens)%laction .and. nrcvinfo(jps_sens) == OASIS_Rcv ) frcv(:,:,jps_sens) = -frcv(:,:,jps_sens)
!MU (18.09.12)
   ! net CO2 tendency (now only photosynthesis)
   IF( srcv(jps_co2fl)%laction .and. nrcvinfo(jps_co2fl) == OASIS_Rcv ) frcv(:,:,jps_co2fl) = -frcv(:,:,jps_co2fl)
!MU (12.04.13)
   ! photosynthesis tendency
   IF( srcv(jps_fpsn)%laction .and. nrcvinfo(jps_fpsn) == OASIS_Rcv ) frcv(:,:,jps_fpsn) = -frcv(:,:,jps_fpsn)
   ! plant respiration tendency
   IF( srcv(jps_fplres)%laction .and. nrcvinfo(jps_fplres) == OASIS_Rcv ) frcv(:,:,jps_fplres) = frcv(:,:,jps_fplres)

!! CPS Update of COSMO variables consistent with masking, for masked points we
!! retain the COSMO values.

   DO j = jstartpar, jendpar
   jm1 = MAX(j-1,1)                !CPS
   DO i = istartpar, iendpar
   im1 = MAX(i-1,1)                !CPS

!CPS   IF (llandmask(i,j)) THEN        !land points only
   IF (zmask(i,j).eq.0) THEN       !CPS MASK


   IF( srcv(jps_co2fl)%laction ) THEN
!CPS     co2fl_s(i,j) = frcv(i,j,jps_co2fl) * 1.e-6_wp * (molmass_co2 * 0.001) / &
!CPS                                             ( rho(i,j,ke) * ( hhl(i,j,ke) - hhl(i,j,ke1) ) )
   ELSE
!CPS     co2fl_s(i,j) = 0.
   ENDIF
!MU (18.09.12)
   IF( srcv(jps_fpsn)%laction ) THEN
!CPS     psn_tens(i,j) = frcv(i,j,jps_fpsn) * 1.e-6_wp * (molmass_co2 * 0.001) / &
!CPS                                             ( rho(i,j,ke) * ( hhl(i,j,ke) - hhl(i,j,ke1) ) )
   ELSE
!CPS     psn_tens(i,j) = 0.
   ENDIF

   IF( srcv(jps_fplres)%laction ) THEN
!CPS     plres_tens(i,j) = frcv(i,j,jps_fplres) * 1.e-6_wp * (molmass_co2 * 0.001) / &
!CPS                                             ( rho(i,j,ke) * ( hhl(i,j,ke) - hhl(i,j,ke1) ) )
   ELSE
!CPS     plres_tens(:,:) = 0.
   ENDIF
!MU (12.04.13)

   IF( srcv(jps_taux)%laction ) THEN
     umfl_s(i,j) = frcv(i,j,jps_taux)
   ENDIF

   IF( srcv(jps_tauy)%laction ) THEN
     vmfl_s(i,j) = frcv(i,j,jps_tauy)
   ENDIF

   IF( srcv(jps_lat)%laction ) THEN
     lhfl_s(i,j) = frcv(i,j,jps_lat) 
   ENDIF

   IF( srcv(jps_sens)%laction ) THEN
     shfl_s(i,j) = frcv(i,j,jps_sens)
   ENDIF

#ifdef COUP_OAS_COS
!CPS alb_rad dimension changed for coupling
   IF( srcv(jps_albd)%laction ) THEN
    alb_rad(i,j,1) = frcv(i,j,jps_albd) 
   ENDIF

   IF( srcv(jps_albi)%laction ) THEN
    alb_rad(i,j,2) = frcv(i,j,jps_albi)
   ENDIF

   IF( srcv(jps_ir)%laction ) THEN
    trad_clm(i,j) =  (frcv(i,j,jps_ir) / sigma / (1._wp - ctalb))**0.25_wp
    t_g(i,j,nnew) = trad_clm(i,j)   !FOR cpl_scheme = False, else replaced below
   ENDIF
#endif
   
   IF( srcv(jps_ram1)%laction ) THEN
     tcm(i,j) = frcv(i,j,jps_ram1)   !CPS need to non-dimensionalize this later 
   !IF (my_cart_id <= 7) THEN
   ! PRINT*, "CPS DEBUG I, J ",my_cart_id, i, j, " TCM ", tcm(i,j)
   !ENDIF
   ENDIF

   IF( srcv(jps_rah1)%laction ) THEN
     tch(i,j) = frcv(i,j,jps_rah1)   !CPS need to non-dimensionalize this later 
   ENDIF

#ifdef COUP_OAS_COS
   IF( srcv(jps_raw1)%laction ) THEN
     tcw(i,j) = frcv(i,j,jps_raw1)   !CPS need to non-dimensionalize this later 
   ENDIF
#endif

   IF( srcv(jps_tsf1)%laction ) THEN
     t_g(i,j,nnew) = frcv(i,j,jps_tsf1)   !CPS Surface Temperature 
   ENDIF

   IF( srcv(jps_qsf1)%laction ) THEN
     qv_s(i,j,nnew) = frcv(i,j,jps_qsf1)   !CPS Surface Moisture 
   ENDIF

   t_g(i,j,nx)   = t_g(i,j,nnew)           !Surface temperture 
   t_s(i,j,nx)   = t_g(i,j,nnew)           !Surface temperature
   t_s(i,j,nnew) = t_g(i,j,nnew)
   qv_s(i,j,nx)  = qv_s(i,j,nnew)          !Surface humidity
!*****************************************************************************************************
!***  Compute COSMO transfer coefficients                                                 ************ 
!    itype_turb = 3 , imode_turb = 1 (Neumann boundary conditions for heat and moisture transport    *
!                       at the lower boundary (speciffied fluxes)                         ************
!*****************************************************************************************************

!CPS   DO j = jstartpar, jendpar
!CPS     jm1 = MAX(j-1,1)                !CPS
!CPS     DO i = istartpar, iendpar
!CPS       im1 = MAX(i-1,1)                !CPS
 
!CPS       IF (llandmask(i,j)) THEN     !CPS changed for land-points only

         zpp         = p0(i,j,ke) + pp(i,j,ke,nx)
         dzke        = hhl(i,j,ke) - hhl(i,j,ke1)
         zpia(i,j)   = zfpi( zpp )
         zpianf(i,j) = zfpi( ps(i,j,nx) )
         zvbke       = 0.5*SQRT ( (u(i,j,ke,nx) + u(im1,j,ke,nx))**2     &
                                +(v(i,j,ke,nx) + v(i,jm1,ke,nx))**2 )
         ztvb        = t_g (i,j,nx)*(1.0 + rvd_m_o*qv_s(i,j,nx))

!CPS
     IF (cpl_scheme) THEN
        !CPS non-dimensionlized conductance for air
         tcm(i,j) = 1._wp/(tcm(i,j)*MAX(zvbke,vel_min))   !Estimated based on CLM
         !IF (my_cart_id == 0) THEN
         !  PRINT*, "CPS DEBUG receive_fld_cos I, J ", i, j, tcm(i,j)
         !ENDIF
         tch(i,j) = 1._wp/(tch(i,j)*MAX(zvbke,vel_min))   !Estimated based on CLM
#ifdef COUP_OAS_COS
         tcw(i,j) = 1._wp/(tcw(i,j)*MAX(zvbke,vel_min))   !Estimated based on CLM
#endif

     ELSE

     ! Sensible heat flux inversion: derive new tch from clm flux
         tch(i,j) = -shfl_s (i,j) / (                                 &
                  zvbke*g*ps(i,j,nx)/(r_d*ztvb)*gr*cp_d*              &
                  ( t_g(i,j,nx) - zpianf(i,j)*t(i,j,ke,nx)/zpia(i,j) ) )
         IF (zvbke == 0) tch(i,j) = 0._wp

!CPS Some fixes for real case simulation 
        tch_epsi=1.0
        coef_lim=tch_epsi*dzke/(zvbke*dt)
        IF (ABS(tch(i,j))>coef_lim) THEN 
          tch(i,j)=SIGN(coef_lim,tch(i,j))
        END IF  
!CPS Some fixes for real case simulation
#ifdef COUP_OAS_COS
         tcw(i,j) = tch(i,j)
#endif

    ! Latent heat flux inversion: derive new qv_s from clm flux
        !CPS
        CALL trcr_get(ierror, idt_qv, ptr_tlev = nx, ptr = qv)
        IF (ierror /= 0) THEN
          yzerrmsg = trcr_errorstr(ierror)
          CALL model_abort(my_cart_id, ierror, yzerrmsg, yzroutine)
        ENDIF
        !CPS

        IF (tch(i,j) .NE. 0) THEN
          ztch(i,j)  = tch(i,j)*zvbke*g*ps(i,j,nx)/(r_d*ztvb)
          qvsflx(i,j) = lhfl_s(i,j) / lh_v
          qv_s(i,j,nx) = -qvsflx(i,j) / ( ztch(i,j)*gr)  + qv(i,j,ke)      
          qv_s(i,j,nnew) = qv_s(i,j,nx)
        ENDIF

! Momentum flux inversion: derive new 6tcm from clm3.5 momentum flux(use either u or v ) !CPS
!
        umass(i,j) = 0.5_wp *(u(i,j,ke,nx) + u(im1,j,ke,nx))
        vmass(i,j) = 0.5_wp *(v(i,j,ke,nx) + v(im1,j,ke,nx))
        IF (umass(i,j) < 0.) umfl_s(i,j) = -umfl_s(i,j)    !CPS correct sign, clm u always positive
        IF (vmass(i,j) < 0.) vmfl_s(i,j) = -vmfl_s(i,j)    !CPS correct sign, clm v always positive

! umfl , vmfl from CLM are in mass points and tcm, tch is also in mass point
        IF ( zvbke == 0 .OR. umass(i,j) == 0) THEN         !CPS nonzero denominator
          tcm(i,j) = 0._wp !CPS stable condition
        ELSE
          tcm(i,j) = (umfl_s(i,j)*r_d*ztvb)/(zvbke*g*ps(i,j,nx)*gr * umass(i,j))  
        ENDIF
 
! Some fixes for real case simulation 
        tcm(i,j) = ABS(tcm(i,j))
        tcm_epsi=0.1
        coef_lim=tcm_epsi*dzke/(zvbke*dt)
        IF (ABS(tcm(i,j))>coef_lim) THEN
           tcm(i,j)=SIGN(coef_lim,tcm(i,j))
        END IF    
        !CPS Some fixes for real case simulation
 
         ! lw radiation is coupled through updated t_g (no change in src_radiation.f90)
         ! sw radiation is coupled via albedo alb_rad (see change in src_radiation.f90)

       ENDIF                      ! Type of coupling scheme

       ENDIF                      ! CPS MASK for land-points only
     ENDDO                        ! i loop 
   ENDDO                          ! j loop
      !
 END IF                           ! coupling on valid PE only

!CPS WHERE ( fr_land < 0.5 )   t_g(:,:,nnew) = t_s(:,:,nx)

 CALL MPI_Barrier(kl_comm, ierror)

!------------------------------------------------------------------------------
!- End of the Subroutine
!------------------------------------------------------------------------------

END SUBROUTINE receive_fld_2clm
