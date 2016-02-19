SUBROUTINE send_fld_2clm

!---------------------------------------------------------------------
! Description:
!  This routine sends coupling fields to CLM3.5 
!
! References:
!  CEREFACS/ETH: E. Maisonnave, Edoward Davin
!
! Current Code Owner: TR32, Z4: Prabhakar Shrestha
!    phone: 0228733453
!    email: pshrestha@uni-bonn.de
!
! History:
! Version    Date       Name
! ---------- ---------- ----
! 1.1        2011/11/28 Prabhakar Shrestha 
!   Modfied and Implemented in COSMO4.11, Initial release
! 1.2        2013/01/15 Markus Uebel, Prabhakar Shrestha 
!   Adding CO2 fields and implemented in COSMO4.21
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
USE data_runcontrol,  ONLY :  ntstep, nnow, itype_gscp, nstop

USE data_modelconfig, ONLY :  ke, dt, ie, je, ie_tot, je_tot,     &
                              istartpar, iendpar, jstartpar, jendpar,idt_qv 
USE data_parameters,  ONLY :  wp, iintegers
USE data_parallel,    ONLY :                        &
                              isubpos,              &        ! positions of the subdomains in the total domain
                              nboundlines,          &    
                              my_cart_id

USE data_fields     , ONLY :                        &
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
!CPS                                    qv         ,    & ! specific water vapor content                  (kg/kg)
!MU (13.09.2012)
!CPS                                    qc         ,    & ! specific cloud water content                  (kg/kg)
!CPS                                    qi         ,    & ! specific cloud ice content                    (kg/kg)
!CPS                                    qr         ,    & ! specific rain content                         (kg/kg)
!CPS                                    qs         ,    & ! specific snow content                         (kg/kg)
!CPS                                    qg         ,    & ! specific graupel content                      (kg/kg)
!MU (13.09.2012)
                                    pp         ,    & ! deviation from the reference pressure         ( pa  )

                                ! 6. fields that are computed in the parametrization and dynamics     (unit )
                                ! ---------------------------------------------------------------
                                !   fields of convective and grid-scale precipitation
                                    prr_con     ,   & ! precipitation rate of rain, convective        (kg/m2*s)
                                    prs_con     ,   & ! precipitation rate of snow, convective        (kg/m2*s)
                                    prr_gsp     ,   & ! precipitation rate of rain, grid-scale        (kg/m2*s)
                                    prs_gsp     ,   & ! precipitation rate of snow, grid-scale        (kg/m2*s)
                                    prg_gsp     ,   & ! precipitation rate of graupel, grid-scale     (kg/m2*s)

                                !   fields of the radiation
                                    swdir_s     ,   & ! direct shortwave downward radiation           ( W/m2)
                                    swdifd_s    ,   & ! diffuse shortwave downward radiation          ( W/m2)
                                    lwd_s             ! thermal downward radiation at the ground      ( W/m2)

USE src_tracer,         ONLY: trcr_get, trcr_errorstr            !CPS
USE environment,        ONLY: model_abort                        !CPS

!MU (13.09.2012)
!CPSUSE data_tracer      , ONLY :   &
!CPS                                    tracer, &                 ! tracer concentration                  (kg/kg)
!CPS                                    ntracer, &                ! number of tracers
!CPS                                    molmass_da, &             ! molar mass of dry air                 (g/mol)
!CPS                                    molmass_w, &              ! molar mass of water                   (g/mol)
!CPS                                    molmass_co2, &            ! molar mass of CO2                     (g/mol)
!CPS                                    ltracer                   ! switch for tracer
!MU (13.09.2012)

USE netcdf


!==============================================================================

IMPLICIT NONE

!==============================================================================

! Local Variables

INTEGER, PARAMETER ::   jps_t   =  1            ! temperature
INTEGER, PARAMETER ::   jps_u   =  2            ! u wind
INTEGER, PARAMETER ::   jps_v   =  3            ! v wind
INTEGER, PARAMETER ::   jps_q   =  4            ! specific water vapor content
INTEGER, PARAMETER ::   jps_th  =  5            ! thickness of lowest level (m)
INTEGER, PARAMETER ::   jps_pr  =  6            ! surface pressure (Pa)
INTEGER, PARAMETER ::   jps_rs  =  7            ! direct shortwave downward radiation (W/m2)
INTEGER, PARAMETER ::   jps_fs  =  8            ! diffuse shortwave downward radiation (W/m2)
INTEGER, PARAMETER ::   jps_lw  =  9            ! longwave downward radiation (W/m2)
INTEGER, PARAMETER ::   jps_cr  = 10            ! convective rain precipitation      (kg/m2*s)
INTEGER, PARAMETER ::   jps_cs  = 11            ! convective snow precipitation      (kg/m2*s)
INTEGER, PARAMETER ::   jps_gr  = 12            ! gridscale rain precipitation
INTEGER, PARAMETER ::   jps_gs  = 13            ! gridscale snow precipitation
INTEGER, PARAMETER ::   jps_gg  = 14            ! gridscale graupel precipitation
INTEGER, PARAMETER ::   jps_cp  = 15            ! total convective precipitation
INTEGER, PARAMETER ::   jps_gp  = 16            ! total gridscale precipitation
!MU (13.09.2012)
INTEGER, PARAMETER ::   jps_co2 = 17            ! CO2 partial pressure (Pa)
!MU (13.09.2012)

INTEGER(KIND=iintegers) ::   isec, info         ! temporary integer
REAL(KIND=wp)       ::   ztmp1(ie,je),     &! temporary arrays
                             umass(ie,je),     &!
                             vmass(ie,je),     &!
                             winds(ie,je) 
REAL(KIND=wp)       :: toc

INTEGER                 :: i, j, im1, jm1, iig      ! i-1, j-1, !MU iig: loop index for tracer

!MU (13.09.2012)
REAL (KIND=wp)   ::          &
    molmass_ma(ie,je,ke),&                      ! molar mass of moist air         (g/mol)
    co2_s(ie,je)                                ! CO2 partial pressure at surface ( Pa  )
!MU (13.09.2012)

INTEGER :: il_var_id(17),  dimids(2), il_file_id, mype, ib, npes, ierror, status, ier, jn 
           !MU: il_var_id(16) changed to il_var_id(17) due to CO2 coupling

!CPS qx variables moved from data_fields to tracers
REAL (KIND=wp),     POINTER :: &
  qv  (:,:,:)    => NULL()             ! QV at given time-level 
CHARACTER (LEN=255)        :: yzerrmsg
CHARACTER (LEN=25)         :: yzroutine
!CPS 
!------------------------------------------------------------------------------
!- End of header
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!- Begin Subroutine  
!------------------------------------------------------------------------------

 !CPS QV is moved to tracer variables
 CALL trcr_get(ierror, idt_qv, ptr_tlev = nnow, ptr = qv)
   IF (ierror /= 0) THEN
   yzerrmsg = trcr_errorstr(ierror)
   CALL model_abort(my_cart_id, ierror, yzerrmsg, yzroutine)
 ENDIF
 !CPS

 IF ( IOASISDEBUGLVL == 2 ) THEN
   CALL MPI_Barrier(kl_comm, info)
   toc=MPI_Wtime()
 ENDIF

 ! Write coupling fields for restart file
 IF ( IOASISDEBUGLVL == 3 ) THEN

   CALL MPI_Comm_rank(kl_comm, mype, ier)
   CALL MPI_Comm_size(kl_comm, npes, ier)

   IF (mype == 0) then
     status = nf90_create("cosmo.nc", NF90_CLOBBER, il_file_id)
     status = nf90_def_dim(il_file_id, "x", ie_tot - 2 * nboundlines, dimids(1))
     status = nf90_def_dim(il_file_id, "y", je_tot - 2 * nboundlines, dimids(2))
     status = nf90_def_var(il_file_id, "COSTEMPE", NF90_DOUBLE, dimids, il_var_id(1))
     status = nf90_def_var(il_file_id, "COSUWIND", NF90_DOUBLE, dimids, il_var_id(2))
     status = nf90_def_var(il_file_id, "COSVWIND", NF90_DOUBLE, dimids, il_var_id(3))
     status = nf90_def_var(il_file_id, "COSSPWAT", NF90_DOUBLE, dimids, il_var_id(4))
     status = nf90_def_var(il_file_id, "COSTHICK", NF90_DOUBLE, dimids, il_var_id(5))
     status = nf90_def_var(il_file_id, "COSPRESS", NF90_DOUBLE, dimids, il_var_id(6))
     status = nf90_def_var(il_file_id, "COSDIRSW", NF90_DOUBLE, dimids, il_var_id(7))
     status = nf90_def_var(il_file_id, "COSDIFSW", NF90_DOUBLE, dimids, il_var_id(8))
     status = nf90_def_var(il_file_id, "COSLONGW", NF90_DOUBLE, dimids, il_var_id(9))
     ! Coupling separately snow, rain and graupel
     IF ( ssnd(jps_cr)%laction ) THEN
        status = nf90_def_var(il_file_id, "COSCVRAI", NF90_DOUBLE, dimids, il_var_id(10))
        status = nf90_def_var(il_file_id, "COSCVSNW", NF90_DOUBLE, dimids, il_var_id(11))
        status = nf90_def_var(il_file_id, "COSGSRAI", NF90_DOUBLE, dimids, il_var_id(12))
        status = nf90_def_var(il_file_id, "COSGSSNW", NF90_DOUBLE, dimids, il_var_id(13))
        status = nf90_def_var(il_file_id, "COSGRAUP", NF90_DOUBLE, dimids, il_var_id(14))
     ELSE
        status = nf90_def_var(il_file_id, "COSCVPRE", NF90_DOUBLE, dimids, il_var_id(15))
        status = nf90_def_var(il_file_id, "COSGSPRE", NF90_DOUBLE, dimids, il_var_id(16))
     ENDIF
!MU (13.09.2012)
     status = nf90_def_var(il_file_id, "COSCO2PP", NF90_DOUBLE, dimids, il_var_id(17))
!MU (13.09.2012)
        status = nf90_enddef(il_file_id)
        status = nf90_close(il_file_id)
   ENDIF               ! mype=0

   DO ib = 0, npes - 1                        !LOOP npes

     CALL MPI_Barrier(kl_comm, ierror)
 
     IF (mype == ib .AND. lpe_cpl ) THEN
        status = nf90_open("cosmo.nc", NF90_WRITE, il_file_id)
        status = nf90_inq_varid(il_file_id, "COSTEMPE" , il_var_id(1))
        status = nf90_inq_varid(il_file_id, "COSUWIND" , il_var_id(2))
        status = nf90_inq_varid(il_file_id, "COSVWIND" , il_var_id(3))
        status = nf90_inq_varid(il_file_id, "COSSPWAT" , il_var_id(4))
        status = nf90_inq_varid(il_file_id, "COSTHICK" , il_var_id(5))
        status = nf90_inq_varid(il_file_id, "COSPRESS" , il_var_id(6))
        status = nf90_inq_varid(il_file_id, "COSDIRSW" , il_var_id(7))
        status = nf90_inq_varid(il_file_id, "COSDIFSW" , il_var_id(8))
        status = nf90_inq_varid(il_file_id, "COSLONGW" , il_var_id(9))
        IF ( ssnd(jps_cr)%laction ) THEN
           status = nf90_inq_varid(il_file_id, "COSCVRAI" , il_var_id(10))
           status = nf90_inq_varid(il_file_id, "COSCVSNW" , il_var_id(11))
           status = nf90_inq_varid(il_file_id, "COSGSRAI" , il_var_id(12))
           status = nf90_inq_varid(il_file_id, "COSGSSNW" , il_var_id(13))
           status = nf90_inq_varid(il_file_id, "COSGRAUP" , il_var_id(14))
        ELSE
           status = nf90_inq_varid(il_file_id, "COSCVPRE" , il_var_id(15))
           status = nf90_inq_varid(il_file_id, "COSGSPRE" , il_var_id(16))
        ENDIF
!MU (13.09.2012)
        status = nf90_inq_varid(il_file_id, "COSCO2PP" , il_var_id(17))
!MU (13.09.2012)

        DO jn = 1,17       !MU: changed to 17 for CO2 coupling

          IF ( jn == 1 ) ztmp1(:,:) = t(:,:,ke,nnow)
          IF ( jn == 2 ) ztmp1(:,:) = u(:,:,ke,nnow)
          IF ( jn == 3 ) ztmp1(:,:) = v(:,:,ke,nnow)
          IF ( jn == 4 ) ztmp1(:,:) = qv(:,:,ke)
          IF ( jn == 5 ) ztmp1(:,:) = (hhl(:,:,ke) - hsurf(:,:)) * 0.5  !CPS
          IF ( jn == 6 ) ztmp1(:,:) = p0(:,:,ke) + pp(:,:,ke,nnow)
          IF ( jn == 7 ) ztmp1(:,:) = swdir_s(:,:)
          IF ( jn == 8 ) ztmp1(:,:) = swdifd_s(:,:)
          IF ( jn == 9 ) ztmp1(:,:) = lwd_s(:,:)
          IF ( ssnd(jps_cr)%laction ) THEN
            IF ( jn == 10 ) ztmp1(:,:) = prr_con(:,:)
            IF ( jn == 11 ) ztmp1(:,:) = prs_con(:,:)
            IF ( jn == 12 ) ztmp1(:,:) = prr_gsp(:,:)
            IF ( jn == 13 ) ztmp1(:,:) = prs_gsp(:,:)
            IF ( jn == 14 .and. itype_gscp == 4 ) ztmp1(:,:) = prg_gsp(:,:)
            IF ( jn == 14 .and. itype_gscp /= 4 ) ztmp1(:,:) = 0.
          ELSE
            IF ( jn == 15 ) ztmp1(:,:) = prr_con(:,:) + prs_con(:,:)
            IF ( jn == 16 ) ztmp1(:,:) = prr_gsp(:,:) + prs_gsp(:,:)
            IF ( jn == 16 .and. itype_gscp == 4 ) ztmp1(:,:) = ztmp1(:,:) + prg_gsp(:,:)
          ENDIF
!MU (13.09.2012)
!CPS          DO iig=1, ntracer
!CPS            IF ( jn == 17 .and. ltracer(7,iig) == 1 ) ztmp1(:,:) = tracer(:,:,ke,nnow,iig)
!CPS          ENDDO
!MU (13.09.2012)
!CPS                IF ( ssnd(jps_cr)%laction ) &
          status = nf90_put_var(il_file_id, il_var_id(jn),  ztmp1(nldi:nlei, nldj:nlej), &
                                start = (/ isubpos(my_cart_id,1)-nboundlines,            &
                                            isubpos(my_cart_id,2)-nboundlines /),        &
                                count = (/ jih, jjh /) )
        ENDDO

        status = nf90_close(il_file_id)

     ENDIF              !mype == ib .and. lpe_cpl 
    
   ENDDO                !LOOP npes


    ! End write coupling fields
   WRITE(6,*) ' cosmo.nc restart file written'
   WRITE(6,*) ' we stop simulation'
   STOP 

 ENDIF                  !IOASISDEBUGLVL=3

 ! Coupling only on PE with at least one land point
 IF ( lpe_cpl ) THEN

   isec = ntstep * dt

   IF( ssnd(jps_t)%laction )   CALL oas_cos_snd( jps_t, isec, t(:,:,ke,nnow), info )
   ! CPS Need to send u and v at mass point for CLM
   DO j = jstartpar, jendpar
     jm1 = MAX(j-1,1)                !CPS
     DO i = istartpar, iendpar
       im1 = MAX(i-1,1)                !CPS
       umass(i,j)  = 0.5_wp*(u(i,j,ke,nnow) + u(im1,j,ke,nnow))
       vmass(i,j)  = 0.5_wp*(v(i,j,ke,nnow) + v(i,jm1,ke,nnow))
       winds(i,j)  = SQRT(umass(i,j)**2 + vmass(i,j)**2)
   ! IF (my_cart_id == 0) THEN
   ! PRINT*, "CPS DEBUG I, J ",i, j, vmass(i,j),v(i,j,ke,nnow), v(i,jm1,ke,nnow), umass(i,j),0.0_wp
   !ENDIF
     ENDDO
   ENDDO

   IF( ssnd(jps_u)%laction )   CALL oas_cos_snd( jps_u, isec, umass(:,:), info )    !CPS
   IF( ssnd(jps_v)%laction )   CALL oas_cos_snd( jps_v, isec, vmass(:,:), info )    !CPS
   IF( ssnd(jps_q)%laction )   CALL oas_cos_snd( jps_q, isec, qv(:,:,ke), info )

   ! Send only total height
   ztmp1(:,:) = (hhl(:,:,ke) - hsurf(:,:)) * 0.5_wp   !CPS
   IF( ssnd(jps_th)%laction )   CALL oas_cos_snd( jps_th, isec, ztmp1(:,:), info )

   ! Send only total pressure
   ztmp1(:,:) = p0(:,:,ke) + pp(:,:,ke,nnow)
   IF( ssnd(jps_pr)%laction )   CALL oas_cos_snd( jps_pr, isec, ztmp1(:,:), info )

   IF( ssnd(jps_rs)%laction )   CALL oas_cos_snd( jps_rs, isec, swdir_s(:,:), info )
   IF( ssnd(jps_fs)%laction )   CALL oas_cos_snd( jps_fs, isec, swdifd_s(:,:), info )
   IF( ssnd(jps_lw)%laction )   CALL oas_cos_snd( jps_lw, isec, lwd_s(:,:), info )
   IF( ssnd(jps_cr)%laction )   CALL oas_cos_snd( jps_cr, isec, prr_con(:,:), info )
   IF( ssnd(jps_cs)%laction )   CALL oas_cos_snd( jps_cs, isec, prs_con(:,:), info )
   IF( ssnd(jps_gr)%laction )   CALL oas_cos_snd( jps_gr, isec, prr_gsp(:,:), info )
   IF( ssnd(jps_gs)%laction )   CALL oas_cos_snd( jps_gs, isec, prs_gsp(:,:), info )

   ! Send graupel if any, 0 otherwise
   ztmp1(:,:) = 0.
   IF ( itype_gscp == 4 )       ztmp1(:,:) = prg_gsp(:,:)
   IF( ssnd(jps_gg)%laction )   CALL oas_cos_snd( jps_gg, isec, ztmp1(:,:), info )

   ! Send convective and gridscale total precipitations only
   ztmp1(:,:) = prr_con(:,:) + prs_con(:,:)
   IF( ssnd(jps_cp)%laction )   CALL oas_cos_snd( jps_cp, isec, ztmp1(:,:), info )
   ztmp1(:,:) = prr_gsp(:,:) + prs_gsp(:,:)
   IF ( itype_gscp == 4 )       ztmp1(:,:) = ztmp1(:,:) + prg_gsp(:,:)
   IF( ssnd(jps_gp)%laction )   CALL oas_cos_snd( jps_gp, isec, ztmp1(:,:), info )

!MU (13.09.2012)
   ! Need to send CO2 partial pressure for CLM
     ! Calculation of the molar mass of moist air
!CPS     IF (itype_gscp == 4) THEN
!CPS       molmass_ma(:,:,ke) = ( 1 - qv(:,:,ke,nnow) - qc(:,:,ke,nnow) - qi(:,:,ke,nnow) &
!CPS                                - qr(:,:,ke,nnow) - qs(:,:,ke,nnow) - qg(:,:,ke,nnow) ) * molmass_da &
!CPS                            + (   qv(:,:,ke,nnow) + qc(:,:,ke,nnow) + qi(:,:,ke,nnow) &
!CPS                                + qr(:,:,ke,nnow) + qs(:,:,ke,nnow) + qg(:,:,ke,nnow) ) * molmass_w
!CPS     ELSE
!CPS       molmass_ma(:,:,ke) = ( 1 - qv(:,:,ke,nnow) - qc(:,:,ke,nnow) - qi(:,:,ke,nnow) &
!CPS                                - qr(:,:,ke,nnow) - qs(:,:,ke,nnow) ) * molmass_da &
!CPS                            + (   qv(:,:,ke,nnow) + qc(:,:,ke,nnow) + qi(:,:,ke,nnow) &
!CPS                                + qr(:,:,ke,nnow) + qs(:,:,ke,nnow) ) * molmass_w
!CPS     ENDIF

!CPS   DO iig=1, ntracer
!CPS     IF ( ltracer(7,iig) == 1) THEN
!CPS       co2_s(:,:) = tracer(:,:,ke,nnow,iig) * molmass_ma(:,:,ke) / molmass_co2 &
!CPS                                            * ( p0(:,:,ke) + pp(:,:,ke,nnow) )
!CPS     ENDIF
!CPS   ENDDO

  !MU (15.01.2013)
   ! Definition of a dummy if CO2 is not initialized in COSMO --> CO2 content from CLM is used
!CPS   IF ( sum(ltracer(7,:)) == 0) THEN
     co2_s(:,:) = 39.5_wp !CPS  -999.
!CPS   ENDIF
  !MU (15.01.2013)   

   IF( ssnd(jps_co2)%laction )  CALL oas_cos_snd( jps_co2, isec, co2_s(:,:), info )
!MU (13.09.2012)

   !  Performance measurement
   IF ( IOASISDEBUGLVL == 2 .AND. my_cart_id == 0 ) THEN
     IF ( info == PRISM_Sent     .OR. info == PRISM_ToRest .OR.   &
            & info == PRISM_SentOut  .OR. info == PRISM_ToRestOut ) &
            WRITE(6, *)  'oascos: send_fld_2clm: before send ', toc
   END IF

 END IF   ! coupling on valid PE only

!------------------------------------------------------------------------------
!- End of the Subroutine
!------------------------------------------------------------------------------

END SUBROUTINE send_fld_2clm
