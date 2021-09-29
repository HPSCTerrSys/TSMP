!>
!! @brief turbulent diagnosis for LES physics 
!!
!! <Describe the concepts of the procedures and algorithms used in the module.>
!! <Details of procedures are documented below with their definitions.>
!! <Include any applicable external references inline as module::procedure,>
!! <external_procedure(), or by using @see.>
!! <Don't forget references to literature.>
!!
!! @author <name, affiliation>
!! @author <name, affiliation>
!!
!!
!! @par Revision History
!! first implementation  by Anurag Dipankar, MPIM (2014-01)
!!
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!

!----------------------------
#include "omp_definitions.inc"
!----------------------------

MODULE mo_turbulent_diagnostic


  USE mo_kind,               ONLY: wp
  USE mo_nonhydro_state,     ONLY: p_nh_state ! SBr,CHa
  USE mo_impl_constants_grf, ONLY: grf_bdywidth_c
  USE mo_loopindices,        ONLY: get_indices_c
  USE mo_exception,          ONLY: message, message_text, finish
  USE mo_model_domain,       ONLY: t_patch
  USE mo_run_config,         ONLY: msg_level, iqv, iqc, iqi, iqr, iqs, iqg, iqh, dtime
  USE mo_nonhydro_types,     ONLY: t_nh_prog, t_nh_diag, t_nh_metrics
  USE mo_nwp_phy_types,      ONLY: t_nwp_phy_diag, t_nwp_phy_tend
  USE mo_nwp_lnd_types,      ONLY: t_lnd_prog, t_lnd_diag 
  USE mo_parallel_config,    ONLY: nproma
  USE mo_statistics,         ONLY: levels_horizontal_mean
  USE mo_les_nml,            ONLY: turb_profile_list, turb_tseries_list
  USE mo_les_config,         ONLY: les_config
  USE mo_mpi,                ONLY: my_process_is_stdio
  USE mo_write_netcdf      
  USE mo_impl_constants,     ONLY: min_rlcell_int
  USE mo_physical_constants, ONLY: cpd, grav, alv, vtmpc1
  USE mo_atm_phy_nwp_config, ONLY: atm_phy_nwp_config
  USE mo_satad,              ONLY: sat_pres_water, spec_humi
  USE mtime,                 ONLY: datetime, timeDelta, newTimedelta,             &
    &                              deallocateTimedelta, getTimedeltaFromDatetime, &
    &                              getTotalMillisecondsTimedelta
  USE mo_time_config,        ONLY: time_config
  USE mo_util_phys,          ONLY: cal_cape_cin
  USE mo_nwp_parameters,     ONLY: t_phy_params
  USE mo_ext_data_types,     ONLY: t_external_data

  IMPLICIT NONE

  LOGICAL, ALLOCATABLE :: is_at_full_level(:)

  INTEGER  :: nrec_tseries,   nrec_profile
  INTEGER  :: fileid_tseries, fileid_profile
  INTEGER  :: sampl_freq_step, avg_interval_step
  LOGICAL  :: is_sampling_time, is_writing_time, is_rh_out

  !Some indices: think of better way
  INTEGER  :: idx_sgs_th_flx, idx_sgs_qv_flx, idx_sgs_qc_flx
  INTEGER  :: idx_sgs_u_flx, idx_sgs_v_flx
  
  CHARACTER(20) :: tname     = 'time'
  CHARACTER(20) :: tlongname = 'Time'

  PRIVATE

  
  PUBLIC  :: les_cloud_diag
  PUBLIC  :: calculate_turbulent_diagnostics, write_vertical_profiles, write_time_series
  PUBLIC  :: init_les_turbulent_output, close_les_turbulent_output
  PUBLIC  :: sampl_freq_step, avg_interval_step, is_sampling_time, is_writing_time
  PUBLIC  :: idx_sgs_th_flx, idx_sgs_qv_flx, idx_sgs_qc_flx, idx_sgs_u_flx, idx_sgs_v_flx

CONTAINS

  !> AD: 28 July 2014- more diag yet to be added
  !!
  !! Calculates cloud diagnostics for realistic LES runs when convective 
  !! parameterization is off> !!  Very preliminary for now
  !!
  !! Most of the diagnostics are from mo_nwp_diagnosis/nwp_diag_for_output
  !! routine. Some of them which were very specific to NWP have been deleted.
  !!
  !! @par Revision History
  !!
  SUBROUTINE les_cloud_diag(  kstart_moist,               & !in
                            & ih_clch, ih_clcm,           & !in
                            & phy_params,                 & !in
                            & p_patch, p_metrics,         & !in
                            & p_prog,                     & !in
                            & p_prog_rcf,                 & !in
                            & p_diag,                     & !in
                            & p_diag_land,                & !in
                            & p_prog_land,                & !in
                            & prm_diag,                   & !inout
                            & ext_data                    ) !inout    

    !>
    ! !INPUT PARAMETERS:
    INTEGER                , INTENT(in)   :: kstart_moist
    INTEGER                , INTENT(IN)   :: ih_clch, ih_clcm
    TYPE(t_phy_params)     , INTENT(IN)   :: phy_params

    TYPE(t_patch),   TARGET, INTENT(in)   :: p_patch    !<grid/patch info.
    TYPE(t_lnd_prog),        INTENT(in)   :: p_prog_land
    TYPE(t_lnd_diag),        INTENT(in)   :: p_diag_land
    TYPE(t_nh_diag), TARGET, INTENT(in)   :: p_diag     !<the diagnostic variables
    TYPE(t_nh_prog), TARGET, INTENT(in)   :: p_prog_rcf !<the prognostic variables (with
    TYPE(t_nh_prog), TARGET, INTENT(in)   :: p_prog     !<the prognostic variables
    TYPE(t_nh_metrics), INTENT(in)        :: p_metrics
    TYPE(t_nwp_phy_diag)   , INTENT(inout):: prm_diag
    TYPE(t_external_data),INTENT(IN)  ::ext_data       !< external data

    REAL(wp), PARAMETER :: qc_min = 1.e-8_wp
    REAL(wp), PARAMETER :: grav_o_cpd = grav/cpd
    REAL(wp), PARAMETER :: zundef = -999._wp   ! undefined value for 0 deg C level

    REAL(wp):: zbuoy, zqsat, zcond
    REAL(wp):: ri_no
    REAL(wp):: ztp(nproma), zqp(nproma)

    INTEGER :: found_cltop, found_clbas
    INTEGER :: nlev, nlevp1
    INTEGER :: rl_start, rl_end
    INTEGER :: i_startblk, i_endblk    !> blocks
    INTEGER :: i_startidx, i_endidx    !< slices
    INTEGER :: i_nchdom, jg            !< domain index
    INTEGER :: jc,jk,jb                !block index
    INTEGER :: mtop_min
    INTEGER :: mlab(nproma)

    nlev      = p_patch%nlev 
    nlevp1    = p_patch%nlev+1 

    i_nchdom  = MAX(1,p_patch%n_childdom)
    jg        = p_patch%id

    rl_start   = grf_bdywidth_c+1
    rl_end     = min_rlcell_int
    i_startblk = p_patch%cells%start_blk(rl_start,1)
    i_endblk   = p_patch%cells%end_blk(rl_end,i_nchdom)

    ! minimum top index for dry convection
    mtop_min = (ih_clch+ih_clcm)/2    

   
    prm_diag%locum(:,:) = .FALSE.

!$OMP PARALLEL 
!$OMP DO PRIVATE(jb,jc,jk,i_startidx,i_endidx,found_cltop,found_clbas,ri_no,&
!$OMP            mlab,ztp,zqp,zbuoy,zqsat,zcond) ICON_OMP_DEFAULT_SCHEDULE 
    DO jb = i_startblk,i_endblk
       CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                          i_startidx, i_endidx, rl_start, rl_end)

       !Cloud base and cloud top model levels
       DO jc = i_startidx, i_endidx

         !cloud top
         found_cltop = -1
         DO jk = kstart_moist, nlev-1
           IF(jk.EQ.kstart_moist.AND.p_prog_rcf%tracer(jc,jk,jb,iqc)>qc_min.AND.found_cltop/=1)THEN
             prm_diag%mtop_con(jc,jb) = jk
             found_cltop = 1
           ELSEIF(p_prog_rcf%tracer(jc,jk,jb,iqc)<qc_min.AND. &
               p_prog_rcf%tracer(jc,jk+1,jb,iqc)>qc_min.AND.found_cltop/=1)THEN
             prm_diag%mtop_con(jc,jb) = jk
             found_cltop = 1
           END IF
         END DO

         !cloud base
         found_clbas = -1
         DO jk = nlev, kstart_moist+1, -1
           IF(jk.EQ.nlev.AND.p_prog_rcf%tracer(jc,jk,jb,iqc)>qc_min.AND.found_clbas/=1)THEN !Fog
             prm_diag%mbas_con(jc,jb) = jk    
             found_clbas = 1
           ELSEIF(p_prog_rcf%tracer(jc,jk,jb,iqc)<qc_min.AND. &
               p_prog_rcf%tracer(jc,jk-1,jb,iqc)>qc_min.AND.found_clbas/=1)THEN !otherwise
             prm_diag%mbas_con(jc,jb) = jk    
             found_clbas = 1
           END IF
         END DO

         !Accept only when both top and bottom exist and height of bottom is 
         !lower than that of top
         IF(found_clbas==1.AND.found_cltop==1)THEN
           IF(prm_diag%mtop_con(jc,jb) < prm_diag%mbas_con(jc,jb)) &
              prm_diag%locum(jc,jb) = .TRUE.
         ELSE
           prm_diag%mbas_con(jc,jb) = -1
           prm_diag%mtop_con(jc,jb) = -1
         END IF

       END DO!jc

!  -included calculation of boundary layer height (Anurag Dipankar, MPI Octo 2013).
!   using Bulk richardson number approach. 
       DO jc = i_startidx, i_endidx           
         DO jk = nlev-1, kstart_moist, -1

            ri_no = (grav/p_prog%theta_v(jc,nlev,jb)) * &
                    ( p_prog%theta_v(jc,jk,jb)-p_prog%theta_v(jc,nlev,jb) ) *  &
                    ( p_metrics%z_mc(jc,jk,jb)-p_metrics%z_mc(jc,nlev,jb) ) /  &
                    MAX( 1.e-6_wp,(p_diag%u(jc,jk,jb)**2+p_diag%v(jc,jk,jb)**2) )

            IF(ri_no > 0.28_wp)THEN
               prm_diag%z_pbl(jc,jb) = p_metrics%z_mc(jc,jk,jb)
               EXIT
            END IF
         END DO 
       ENDDO
       !
       ! height of ccloud base and top: hbas_con, htop_con
       ! 
       DO jc = i_startidx, i_endidx
         IF ( prm_diag%locum(jc,jb)) THEN
           prm_diag%hbas_con(jc,jb) = p_metrics%z_ifc( jc, prm_diag%mbas_con(jc,jb), jb)
           prm_diag%htop_con(jc,jb) = p_metrics%z_ifc( jc, prm_diag%mtop_con(jc,jb), jb)
         ELSE
           prm_diag%hbas_con(jc,jb) = zundef
           prm_diag%htop_con(jc,jb) = zundef
         END IF
       ENDDO  ! jc
       !
       ! height of the top of dry convection
       !
       DO jc = i_startidx, i_endidx 
         prm_diag%htop_dc(jc,jb) = zundef
         mlab(jc) = 1
         ztp (jc) = p_diag%temp(jc,nlev,jb) + 0.25_wp
         zqp (jc) = p_prog_rcf%tracer(jc,nlev,jb,iqv)
       ENDDO

       DO jk = nlev-1, mtop_min, -1
         DO jc = i_startidx, i_endidx 
           IF ( mlab(jc) == 1) THEN
             ztp(jc) = ztp(jc)  - grav_o_cpd*( p_metrics%z_mc(jc,jk,jb)    &
            &                                 -p_metrics%z_mc(jc,jk+1,jb) )
             zbuoy = ztp(jc)*( 1._wp + vtmpc1*zqp(jc) ) - p_diag%tempv(jc,jk,jb)
             zqsat = spec_humi( sat_pres_water(ztp(jc)), p_diag%pres(jc,jk,jb) )
             zcond = zqp(jc) - zqsat

             IF ( zcond < 0._wp .AND. zbuoy > 0._wp) THEN
               prm_diag%htop_dc(jc,jb) = p_metrics%z_ifc(jc,jk,jb)
             ELSE
               mlab(jc) = 0
             END IF
           END IF
         ENDDO
       ENDDO

       DO jc = i_startidx, i_endidx 
         IF ( prm_diag%htop_dc(jc,jb) > zundef) THEN
           prm_diag%htop_dc(jc,jb) = MIN( prm_diag%htop_dc(jc,jb),        &
          &                p_metrics%z_ifc(jc,nlevp1,jb) + 3000._wp )
           IF ( prm_diag%locum(jc,jb)) THEN
             prm_diag%htop_dc(jc,jb) = MIN( prm_diag%htop_dc(jc,jb),      &
            &                               prm_diag%hbas_con(jc,jb) )
           END IF
         ELSE
           prm_diag%htop_dc(jc,jb) = MIN( 0._wp, p_metrics%z_ifc(jc,nlevp1,jb) )
         END IF
       ENDDO
       ! 
       ! Compute wind speed in 10m
       ! 
       IF (atm_phy_nwp_config(jg)%inwp_turb > 0 ) THEN
         DO jc = i_startidx, i_endidx
           prm_diag%sp_10m(jc,jb) = SQRT(prm_diag%u_10m(jc,jb)**2 &
             &                    +      prm_diag%v_10m(jc,jb)**2 )
         ENDDO
       ENDIF
      
       !
       !Extended diagnostics for HDCP2
       !
      
       !Temperature and pressure at cloud base and top
       DO jc = i_startidx, i_endidx
         IF ( prm_diag%locum(jc,jb)) THEN
           prm_diag%t_cbase(jc,jb) = p_diag%temp(jc, prm_diag%mbas_con(jc,jb), jb)
           prm_diag%p_cbase(jc,jb) = p_diag%pres(jc, prm_diag%mbas_con(jc,jb), jb)

           prm_diag%t_ctop(jc,jb)  = p_diag%temp(jc, prm_diag%mtop_con(jc,jb), jb)
           prm_diag%p_ctop(jc,jb)  = p_diag%pres(jc, prm_diag%mtop_con(jc,jb), jb)
         ELSE
           prm_diag%t_cbase(jc,jb) = zundef
           prm_diag%p_cbase(jc,jb) = zundef
           prm_diag%t_ctop(jc,jb)  = zundef
           prm_diag%p_ctop(jc,jb)  = zundef
         END IF
       ENDDO  ! jc

      !
      !  CAPE and CIN of mean surface layer parcel
      !
      !  start level (kmoist) is limited to pressure heights above p=60hPa,
      !  in order to avoid unphysically low test parcel temperature.
      !  Otherwise computation crashes in sat_pres_water
      CALL cal_cape_cin( i_startidx, i_endidx,                     &
        &                kmoist  = MAX(kstart_moist,phy_params%k060), & !in
        &                te      = p_diag%temp(:,:,jb)          , & !in
        &                qve     = p_prog_rcf%tracer(:,:,jb,iqv), & !in
        &                prs     = p_diag%pres(:,:,jb)          , & !in
        &                hhl     = p_metrics%z_ifc(:,:,jb)       , & !in
        &                cape_ml = prm_diag%cape_ml(:,jb)        , & !in
        &                cin_ml  = prm_diag%cin_ml(:,jb) )


      ! Compute resolved surface drag: ps * del(orog)
 
      DO jc = i_startidx, i_endidx
         prm_diag%drag_u_grid(jc,jb) = p_diag%pres_ifc(jc,nlevp1,jb) * ext_data%atm%grad_topo(1,jc,jb)
         prm_diag%drag_v_grid(jc,jb) = p_diag%pres_ifc(jc,nlevp1,jb) * ext_data%atm%grad_topo(2,jc,jb)
      ENDDO

    ENDDO  ! jb
!$OMP END DO
!$OMP END PARALLEL  


  END SUBROUTINE les_cloud_diag


  !>
  !! <Calculates 1D and 0D turbulent diagnostics>
  !!
  !! <Describe the purpose of the subroutine and its algorithm(s).>
  !! <Include any applicable external references inline as module::procedure,>
  !! <external_procedure(), or by using @see.>
  !! <Don't forget references to literature.>
  !!
  !! @par Revision History
  !!
  SUBROUTINE calculate_turbulent_diagnostics(             &
                            & p_patch,                    & !in
                            & p_prog,   p_prog_rcf,       & !in
                            & p_diag,                     & !in
                            & p_prog_land, p_diag_land,   & !in
                            & phy_tend,                   & !in
                            & prm_diag                )     !inout
                            

    !>
    ! !INPUT PARAMETERS:

    TYPE(t_patch),   TARGET, INTENT(in)   :: p_patch    !<grid/patch info.
    TYPE(t_nh_diag), TARGET, INTENT(in)   :: p_diag     !<the diagnostic variables
    TYPE(t_nh_prog), TARGET, INTENT(in)   :: p_prog     !<the prognostic variables
    TYPE(t_nh_prog), TARGET, INTENT(in)   :: p_prog_rcf !<the prognostic variables (with
                                                        !< red. calling frequency for tracers!

    TYPE(t_lnd_prog),        INTENT(in)   :: p_prog_land
    TYPE(t_lnd_diag),        INTENT(in)   :: p_diag_land
    TYPE(t_nwp_phy_tend),TARGET,INTENT(in):: phy_tend    !< atm tend vars
    TYPE(t_nwp_phy_diag)   , INTENT(inout):: prm_diag

    ! Local
  
    REAL(wp), ALLOCATABLE, DIMENSION(:,:,:)  :: var3df, var3dh, theta, w_mc
    REAL(wp), ALLOCATABLE, DIMENSION(:)   :: &
              umean, vmean, thmean, qvmean, qcmean, wmean, outvar, thvmean
    REAL(wp) :: outvar0d, w_loc, th_loc, wth_loc, day_sec

    ! Local array bounds:

    INTEGER :: nlev, nlevp1            !< number of full levels
    INTEGER :: rl_start, rl_end
    INTEGER :: i_startblk, i_endblk    !> blocks
    INTEGER :: i_startidx, i_endidx    !< slices
    INTEGER :: i_nchdom                !< domain index
    INTEGER :: jc,jk,jb,jg             !block index
    INTEGER :: nvar, n, ilc1, ibc1, ilc2, ibc2, ilc3, ibc3
    CHARACTER(len=*), PARAMETER :: routine = 'mo_turbulent_diagnostic:calculate_turbulent_diagnostics'

    IF(msg_level>18) & 
      CALL message(routine,'Start!')

    day_sec = 86400._wp

    jg         = p_patch%id
    nlev       = p_patch%nlev
    nlevp1     = nlev + 1
    
    !allocation
    ALLOCATE( var3df(nproma,nlev,p_patch%nblks_c), var3dh(nproma,nlevp1,p_patch%nblks_c), &
              theta(nproma,nlev,p_patch%nblks_c),  w_mc(nproma,nlev,p_patch%nblks_c), &
              outvar(nlevp1) )

    i_nchdom  = MAX(1,p_patch%n_childdom)

    rl_start   = grf_bdywidth_c
    rl_end     = min_rlcell_int-1  !for wthsfs
    i_startblk = p_patch%cells%start_blk(rl_start,1)
    i_endblk   = p_patch%cells%end_blk(rl_end,i_nchdom)

    !Get w and theta at full levels
!$OMP PARALLEL 
!$OMP DO PRIVATE(jc,jb,jk,i_startidx,i_endidx)
    DO jb = i_startblk,i_endblk
       CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                          i_startidx, i_endidx, rl_start, rl_end)
       DO jk = 1 , nlev
         DO jc = i_startidx, i_endidx
           w_mc(jc,jk,jb) = ( p_prog%w(jc,jk,jb) + p_prog%w(jc,jk+1,jb) ) * 0.5_wp

           theta(jc,jk,jb)  = p_diag%temp(jc,jk,jb)/p_prog%exner(jc,jk,jb)
         END DO
       END DO
    END DO
!$OMP END DO
!$OMP END PARALLEL


    !For diagnostics
    rl_start   = grf_bdywidth_c+1
    rl_end     = min_rlcell_int
    i_startblk = p_patch%cells%start_blk(rl_start,1)
    i_endblk   = p_patch%cells%end_blk(rl_end,i_nchdom)

!======================================================================================
                 !Some vertical profiles
!======================================================================================
    
    nvar = SIZE(turb_profile_list,1)

    !Loop over all variables
    DO n = 1 , nvar

     outvar = 0_wp

     SELECT CASE (TRIM(turb_profile_list(n)))

     CASE('u')

       ALLOCATE(umean(1:nlev))
       CALL levels_horizontal_mean(p_diag%u, p_patch%cells%area, p_patch%cells%owned, outvar(1:nlev))
       umean = outvar(1:nlev)

     CASE('v')

       ALLOCATE(vmean(1:nlev))
       CALL levels_horizontal_mean(p_diag%v, p_patch%cells%area, p_patch%cells%owned, outvar(1:nlev))
       vmean = outvar(1:nlev)

     CASE('w')

       ALLOCATE(wmean(1:nlev))
       CALL levels_horizontal_mean(w_mc, p_patch%cells%area, p_patch%cells%owned, outvar(1:nlev))
       wmean = outvar(1:nlev)

     CASE('thv')

       ALLOCATE(thvmean(1:nlev))
       CALL levels_horizontal_mean(p_prog%theta_v, p_patch%cells%area, p_patch%cells%owned, outvar(1:nlev))
       thvmean = outvar(1:nlev)

     CASE('th') !theta mean

       ALLOCATE(thmean(1:nlev))
       CALL levels_horizontal_mean(theta, p_patch%cells%area, p_patch%cells%owned, outvar(1:nlev))
       thmean = outvar(1:nlev)

     CASE('exner')

       CALL levels_horizontal_mean(p_prog%exner, p_patch%cells%area, p_patch%cells%owned, outvar(1:nlev))

     CASE('rho')

       CALL levels_horizontal_mean(p_prog%rho, p_patch%cells%area, p_patch%cells%owned, outvar(1:nlev))

     CASE('qv')

       ALLOCATE(qvmean(1:nlev))
       CALL levels_horizontal_mean(p_prog_rcf%tracer(:,:,:,iqv),  &
             p_patch%cells%area, p_patch%cells%owned, outvar(1:nlev))
       qvmean = outvar(1:nlev) 

     CASE('qc')

       ALLOCATE(qcmean(1:nlev))
       CALL levels_horizontal_mean(p_prog_rcf%tracer(:,:,:,iqc),  &
            p_patch%cells%area, p_patch%cells%owned, outvar(1:nlev))
       qcmean = outvar(1:nlev)
       
     CASE('wu')

       IF(ALLOCATED(wmean).AND.ALLOCATED(umean))THEN
!$OMP PARALLEL 
!$OMP DO PRIVATE(jb,jc,jk,i_startidx,i_endidx)
        DO jb = i_startblk,i_endblk
          CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                             i_startidx, i_endidx, rl_start, rl_end)
          DO jk = 1 , nlev
            DO jc = i_startidx, i_endidx
             var3df(jc,jk,jb) = (w_mc(jc,jk,jb)-wmean(jk))*(p_diag%u(jc,jk,jb)-umean(jk))*p_prog%rho(jc,jk,jb) 
            END DO
          END DO
        END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL
       ELSE
         CALL finish(routine,'put <wu> after <w> and <u> in the namelist')
       END IF  

       CALL levels_horizontal_mean(var3df, p_patch%cells%area, p_patch%cells%owned, outvar(1:nlev))

     CASE('wv')!At full levels

       IF(ALLOCATED(wmean).AND.ALLOCATED(vmean))THEN
!$OMP PARALLEL 
!$OMP DO PRIVATE(jb,jc,jk,i_startidx,i_endidx)
        DO jb = i_startblk,i_endblk
          CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                             i_startidx, i_endidx, rl_start, rl_end)
          DO jk = 1 , nlev
            DO jc = i_startidx, i_endidx
             var3df(jc,jk,jb) = (w_mc(jc,jk,jb)-wmean(jk))*(p_diag%v(jc,jk,jb)-vmean(jk))*p_prog%rho(jc,jk,jb) 
            END DO
          END DO
        END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL
       ELSE
         CALL finish(routine,'put <wv> after <w> and <v> in the namelist')
       END IF  

       CALL levels_horizontal_mean(var3df, p_patch%cells%area, p_patch%cells%owned, outvar(1:nlev))

     CASE('wth')

       IF(ALLOCATED(wmean).AND.ALLOCATED(thmean))THEN
!$OMP PARALLEL 
!$OMP DO PRIVATE(jb,jc,i_startidx,i_endidx)
        DO jb = i_startblk,i_endblk
          CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                             i_startidx, i_endidx, rl_start, rl_end)
          DO jk = 1 , nlev
            DO jc = i_startidx, i_endidx
             var3df(jc,jk,jb) = (w_mc(jc,jk,jb)-wmean(jk))*(theta(jc,jk,jb)-thmean(jk))*p_prog%rho(jc,jk,jb) 
             p_nh_state(1)%metrics%wth(jc,jk,jb) = p_nh_state(1)%metrics%wth(jc,jk,jb) + var3df(jc,jk,jb) !SBr,CHa: output 3D wth
            END DO
          END DO
        END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL
       ELSE
         CALL finish(routine,'put <wth> after <w> and <th> in the namelist')
       END IF  

       CALL levels_horizontal_mean(var3df, p_patch%cells%area, p_patch%cells%owned, outvar(1:nlev))
       outvar(1:nlev) = outvar(1:nlev) * cpd           
  
     CASE('wthv')

       IF(ALLOCATED(wmean).AND.ALLOCATED(thvmean))THEN
!$OMP PARALLEL 
!$OMP DO PRIVATE(jb,jc,jk,i_startidx,i_endidx)
        DO jb = i_startblk,i_endblk
          CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                             i_startidx, i_endidx, rl_start, rl_end)
          DO jk = 1 , nlev
            DO jc = i_startidx, i_endidx
             var3df(jc,jk,jb) = (w_mc(jc,jk,jb)-wmean(jk))*(p_prog%theta_v(jc,jk,jb)-thvmean(jk))*p_prog%rho(jc,jk,jb) 
            END DO
          END DO
        END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL
       ELSE
         CALL finish(routine,'put <wthv> after <w> and <thv> in the namelist')
       END IF  

       CALL levels_horizontal_mean(var3df, p_patch%cells%area, p_patch%cells%owned, outvar(1:nlev))
       outvar(1:nlev) = outvar(1:nlev) * cpd           

     CASE('wqv')

       IF(ALLOCATED(wmean).AND.ALLOCATED(qvmean))THEN
!$OMP PARALLEL 
!$OMP DO PRIVATE(jb,jc,jk,i_startidx,i_endidx)
        DO jb = i_startblk,i_endblk
          CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                             i_startidx, i_endidx, rl_start, rl_end)
          DO jk = 1 , nlev
            DO jc = i_startidx, i_endidx
             var3df(jc,jk,jb) = (w_mc(jc,jk,jb)-wmean(jk))*(p_prog_rcf%tracer(jc,jk,jb,iqv)-qvmean(jk))*p_prog%rho(jc,jk,jb) 
             p_nh_state(1)%metrics%wqv(jc,jk,jb) = p_nh_state(1)%metrics%wqv(jc,jk,jb) + var3df(jc,jk,jb) !SBr,CHa: output 3D wqv
            END DO
          END DO
        END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL
       ELSE
         CALL finish(routine,'put <wqv> after <w> and <qv> in the namelist')
       END IF  

       CALL levels_horizontal_mean(var3df, p_patch%cells%area, p_patch%cells%owned, outvar(1:nlev))
       outvar(1:nlev) = outvar(1:nlev) * alv

     CASE('wqc')

       IF(ALLOCATED(wmean).AND.ALLOCATED(qcmean))THEN
!$OMP PARALLEL 
!$OMP DO PRIVATE(jb,jc,jk,i_startidx,i_endidx)
        DO jb = i_startblk,i_endblk
          CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                             i_startidx, i_endidx, rl_start, rl_end)
          DO jk = 1 , nlev
            DO jc = i_startidx, i_endidx
             var3df(jc,jk,jb) = (w_mc(jc,jk,jb)-wmean(jk))*(p_prog_rcf%tracer(jc,jk,jb,iqc)-qcmean(jk))*p_prog%rho(jc,jk,jb) 
            END DO
          END DO
        END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL
       ELSE
         CALL finish(routine,'put <wqc> after <w> and <qc> in the namelist')
       END IF  

       CALL levels_horizontal_mean(var3df, p_patch%cells%area, p_patch%cells%owned, outvar(1:nlev))
       outvar(1:nlev) = outvar(1:nlev) * alv

     CASE('ww')

       IF(ALLOCATED(wmean))THEN
!$OMP PARALLEL 
!$OMP DO PRIVATE(jb,jc,jk,i_startidx,i_endidx)
        DO jb = i_startblk,i_endblk
          CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                             i_startidx, i_endidx, rl_start, rl_end)
          DO jk = 1 , nlev
            DO jc = i_startidx, i_endidx
             var3df(jc,jk,jb) = (w_mc(jc,jk,jb)-wmean(jk))**2 
            END DO
          END DO
        END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL
       ELSE
         CALL finish(routine,'put <ww> after <w> in the namelist')
       END IF  

       CALL levels_horizontal_mean(var3df, p_patch%cells%area, p_patch%cells%owned, outvar(1:nlev))

     CASE('thth')

       IF(ALLOCATED(thmean))THEN
!$OMP PARALLEL 
!$OMP DO PRIVATE(jb,jc,jk,i_startidx,i_endidx)
        DO jb = i_startblk,i_endblk
          CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                             i_startidx, i_endidx, rl_start, rl_end)
          DO jk = 1 , nlev
            DO jc = i_startidx, i_endidx
             var3df(jc,jk,jb) = (theta(jc,jk,jb)-thmean(jk))**2 
            END DO
          END DO
        END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL
       ELSE
         CALL finish(routine,'put <thth> after <th> in the namelist')
       END IF  

       CALL levels_horizontal_mean(var3df, p_patch%cells%area, p_patch%cells%owned, outvar(1:nlev))

     CASE('qvqv')

       IF(ALLOCATED(qvmean))THEN
!$OMP PARALLEL 
!$OMP DO PRIVATE(jb,jc,jk,i_startidx,i_endidx)
        DO jb = i_startblk,i_endblk
          CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                             i_startidx, i_endidx, rl_start, rl_end)
          DO jk = 1 , nlev
            DO jc = i_startidx, i_endidx
             var3df(jc,jk,jb) = (p_prog_rcf%tracer(jc,jk,jb,iqv)-qvmean(jk))**2 
            END DO
          END DO
        END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL
       ELSE
         CALL finish(routine,'put <qvqv> after <qv> in the namelist')
       END IF  

       CALL levels_horizontal_mean(var3df, p_patch%cells%area, p_patch%cells%owned, outvar(1:nlev))

     CASE('qcqc')

       IF(ALLOCATED(qcmean))THEN
!$OMP PARALLEL 
!$OMP DO PRIVATE(jb,jc,jk,i_startidx,i_endidx)
        DO jb = i_startblk,i_endblk
          CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                             i_startidx, i_endidx, rl_start, rl_end)
          DO jk = 1 , nlev
            DO jc = i_startidx, i_endidx
             var3df(jc,jk,jb) = (p_prog_rcf%tracer(jc,jk,jb,iqc)-qcmean(jk))**2 
            END DO
          END DO
        END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL
       ELSE
         CALL finish(routine,'put <qcqc> after <qc> in the namelist')
       END IF  

       CALL levels_horizontal_mean(var3df, p_patch%cells%area, p_patch%cells%owned, outvar(1:nlev))

     CASE('uu')

       IF(ALLOCATED(umean))THEN
!$OMP PARALLEL 
!$OMP DO PRIVATE(jb,jc,jk,i_startidx,i_endidx)
        DO jb = i_startblk,i_endblk
          CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                             i_startidx, i_endidx, rl_start, rl_end)
          DO jk = 1 , nlev
            DO jc = i_startidx, i_endidx
             var3df(jc,jk,jb) = (p_diag%u(jc,jk,jb)-umean(jk))**2 
            END DO
          END DO
        END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL
       ELSE
         CALL finish(routine,'put <uu> after <u> in the namelist')
       END IF  

       CALL levels_horizontal_mean(var3df, p_patch%cells%area, p_patch%cells%owned, outvar(1:nlev))

     CASE('vv')

       IF(ALLOCATED(vmean))THEN
!$OMP PARALLEL 
!$OMP DO PRIVATE(jb,jc,jk,i_startidx,i_endidx)
        DO jb = i_startblk,i_endblk
          CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                             i_startidx, i_endidx, rl_start, rl_end)
          DO jk = 1 , nlev
            DO jc = i_startidx, i_endidx
             var3df(jc,jk,jb) = (p_diag%v(jc,jk,jb)-vmean(jk))**2 
            END DO
          END DO
        END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL
       ELSE
         CALL finish(routine,'put <vv> after <v> in the namelist')
       END IF  

       CALL levels_horizontal_mean(var3df, p_patch%cells%area, p_patch%cells%owned, outvar(1:nlev))

     CASE('kh')

       CALL levels_horizontal_mean(prm_diag%tkvh, p_patch%cells%area, p_patch%cells%owned, outvar(1:nlevp1))

     CASE('km')

       CALL levels_horizontal_mean(prm_diag%tkvm, p_patch%cells%area, p_patch%cells%owned, outvar(1:nlevp1))

     CASE('bruvais')

       CALL levels_horizontal_mean(prm_diag%bruvais,p_patch%cells%area,p_patch%cells%owned,outvar(1:nlevp1))
       outvar(1)      = outvar(2) 
       outvar(nlevp1) = outvar(nlev) 

     CASE('mechprd')
       !Mechanical production term: prm_diag%mech_prod / 2
       CALL levels_horizontal_mean(prm_diag%mech_prod, p_patch%cells%area,  &
                                   p_patch%cells%owned, outvar(1:nlevp1))
       outvar = outvar * 0.5_wp          
       outvar(1)      = outvar(2) 
       outvar(nlevp1) = outvar(nlev) 

     CASE('wthsfs')!subfilter scale flux: see Erlebacher et al. 1992

!$OMP PARALLEL 
!$OMP DO PRIVATE(jb,jc,jk,i_startidx,i_endidx,ilc1,ibc1,ilc2,ibc2,ilc3,ibc3, &
!$OMP            w_loc,th_loc,wth_loc)
      DO jb = i_startblk,i_endblk
         CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                            i_startidx, i_endidx, rl_start, rl_end)
         DO jk = 1 , nlev
           DO jc = i_startidx, i_endidx
             ilc1 = p_patch%cells%neighbor_idx(jc,jb,1)
             ibc1 = p_patch%cells%neighbor_blk(jc,jb,1)
             ilc2 = p_patch%cells%neighbor_idx(jc,jb,2)
             ibc2 = p_patch%cells%neighbor_blk(jc,jb,2)
             ilc3 = p_patch%cells%neighbor_idx(jc,jb,3)
             ibc3 = p_patch%cells%neighbor_blk(jc,jb,3)
           
             !Use averaging over neighboring cells to mimic test filter twice the grid size
             w_loc  = 0.25_wp*(w_mc(jc,jk,jb)+w_mc(ilc1,jk,ibc1)+w_mc(ilc2,jk,ibc2)+w_mc(ilc3,jk,ibc3))
             th_loc = 0.25_wp*(theta(jc,jk,jb)+theta(ilc1,jk,ibc1)+theta(ilc2,jk,ibc2)+theta(ilc3,jk,ibc3))
             wth_loc= 0.25_wp*(w_mc(jc,jk,jb)*theta(jc,jk,jb)+w_mc(ilc1,jk,ibc1)*theta(ilc1,jk,ibc1)+ &
                          w_mc(ilc2,jk,ibc2)*theta(ilc2,jk,ibc2)+w_mc(ilc3,jk,ibc3)*theta(ilc3,jk,ibc3))
        
             var3df(jc,jk,jb) = (wth_loc - w_loc*th_loc)*p_prog%rho(jc,jk,jb) 
           END DO
         END DO
      END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

       CALL levels_horizontal_mean(var3df, p_patch%cells%area, p_patch%cells%owned, outvar(1:nlev))
       outvar(1:nlev) = outvar(1:nlev) * cpd
 
     CASE('rh')
       IF(is_rh_out) &
         CALL levels_horizontal_mean(prm_diag%rh, p_patch%cells%area,  &
                                     p_patch%cells%owned, outvar(1:nlev))
     CASE('clc')
       CALL levels_horizontal_mean(prm_diag%clc, p_patch%cells%area,  &
                                   p_patch%cells%owned, outvar(1:nlev))
     CASE('qi')
       IF(atm_phy_nwp_config(jg)%inwp_gscp>0) &
       CALL levels_horizontal_mean(p_prog_rcf%tracer(:,:,:,iqi), p_patch%cells%area,  &
                                   p_patch%cells%owned, outvar(1:nlev))
     CASE('qs')
       IF(atm_phy_nwp_config(jg)%inwp_gscp>0) &
       CALL levels_horizontal_mean(p_prog_rcf%tracer(:,:,:,iqs), p_patch%cells%area,  &
                                   p_patch%cells%owned, outvar(1:nlev))
     CASE('qr')
       CALL levels_horizontal_mean(p_prog_rcf%tracer(:,:,:,iqr), p_patch%cells%area,  &
                                   p_patch%cells%owned, outvar(1:nlev))
     CASE('qg')
       IF(ANY((/4,5/) == atm_phy_nwp_config(jg)%inwp_gscp))&
         CALL levels_horizontal_mean(p_prog_rcf%tracer(:,:,:,iqg), p_patch%cells%area,  &
                                     p_patch%cells%owned, outvar(1:nlev))
     CASE('qh')
       IF(ANY((/4,5/) == atm_phy_nwp_config(jg)%inwp_gscp))&
         CALL levels_horizontal_mean(p_prog_rcf%tracer(:,:,:,iqh), p_patch%cells%area,  &
                                     p_patch%cells%owned, outvar(1:nlev))
     CASE('lwf')
       IF(atm_phy_nwp_config(jg)%inwp_radiation>0) &
       CALL levels_horizontal_mean(prm_diag%lwflxall, p_patch%cells%area,  &
                                   p_patch%cells%owned, outvar(1:nlevp1))
     CASE('swf')
       IF(atm_phy_nwp_config(jg)%inwp_radiation>0)THEN
       CALL levels_horizontal_mean(prm_diag%trsolall, p_patch%cells%area,  &
                                   p_patch%cells%owned, outvar(1:nlevp1))
       outvar0d = 0._wp
       CALL levels_horizontal_mean(prm_diag%flxdwswtoa, p_patch%cells%area,  &
                                   p_patch%cells%owned, outvar0d)
       outvar = outvar*outvar0d
       END IF
     CASE('dt_t_sw')
       IF(atm_phy_nwp_config(jg)%inwp_radiation>0) &
       CALL levels_horizontal_mean(phy_tend%ddt_temp_radsw, p_patch%cells%area,  &
                                   p_patch%cells%owned, outvar(1:nlev))
     CASE('dt_t_lw')
       IF(atm_phy_nwp_config(jg)%inwp_radiation>0) &
       CALL levels_horizontal_mean(phy_tend%ddt_temp_radlw, p_patch%cells%area,  &
                                   p_patch%cells%owned, outvar(1:nlev))
     CASE('dt_t_tb')
       IF(atm_phy_nwp_config(jg)%inwp_turb>0) &
       CALL levels_horizontal_mean(phy_tend%ddt_temp_turb, p_patch%cells%area,  &
                                   p_patch%cells%owned, outvar(1:nlev))
     CASE('dt_t_mc')
       IF(atm_phy_nwp_config(jg)%inwp_gscp>0) &
         CALL levels_horizontal_mean(phy_tend%ddt_temp_gscp, p_patch%cells%area,  &
                                     p_patch%cells%owned, outvar(1:nlev))
     CASE DEFAULT !In case calculations are performed somewhere else
      
       outvar = 0._wp
       
     END SELECT

     !Calculate time mean
     IF(is_at_full_level(n))THEN
       prm_diag%turb_diag_1dvar(1:nlev,n) = prm_diag%turb_diag_1dvar(1:nlev,n)+outvar(1:nlev)
     ELSE
       prm_diag%turb_diag_1dvar(1:nlevp1,n) = prm_diag%turb_diag_1dvar(1:nlevp1,n)+outvar(1:nlevp1)
     END IF

    END DO!nvar

!======================================================================================
       !Some time series
!======================================================================================

    nvar = SIZE(turb_tseries_list,1)

    !Loop over all variables
    DO n = 1 , nvar

     outvar0d = 0._wp
     SELECT CASE (TRIM(turb_tseries_list(n)))

     CASE('ccover')
       CALL levels_horizontal_mean(prm_diag%clct, p_patch%cells%area, p_patch%cells%owned, outvar0d)
     CASE('shflx')
       CALL levels_horizontal_mean(prm_diag%shfl_s, p_patch%cells%area, p_patch%cells%owned, outvar0d)
     CASE('lhflx')
       CALL levels_horizontal_mean(prm_diag%lhfl_s, p_patch%cells%area, p_patch%cells%owned, outvar0d)
     CASE('ustress')
       CALL levels_horizontal_mean(prm_diag%umfl_s, p_patch%cells%area, p_patch%cells%owned, outvar0d)
     CASE('vstress')
       CALL levels_horizontal_mean(prm_diag%vmfl_s, p_patch%cells%area, p_patch%cells%owned, outvar0d)
     CASE('tsfc')
       CALL levels_horizontal_mean(p_prog_land%t_g, p_patch%cells%area, p_patch%cells%owned, outvar0d)
     CASE('qsfc')
       CALL levels_horizontal_mean(p_diag_land%qv_s, p_patch%cells%area, p_patch%cells%owned, outvar0d)
     CASE('hbl')
       CALL levels_horizontal_mean(prm_diag%z_pbl, p_patch%cells%area, p_patch%cells%owned, outvar0d)
     CASE('psfc')
       CALL levels_horizontal_mean(p_diag%pres_sfc, p_patch%cells%area, p_patch%cells%owned, outvar0d)
     CASE('swf_tom')
       IF(atm_phy_nwp_config(jg)%inwp_radiation>0) &
       CALL levels_horizontal_mean(prm_diag%swflxtoa, p_patch%cells%area, p_patch%cells%owned, outvar0d)
     CASE('lwf_tom')
       IF(atm_phy_nwp_config(jg)%inwp_radiation>0) &
       CALL levels_horizontal_mean(prm_diag%lwflxall(:,1,:), p_patch%cells%area, p_patch%cells%owned, outvar0d)
     CASE('swf_sfc')
       IF(atm_phy_nwp_config(jg)%inwp_radiation>0) &
       CALL levels_horizontal_mean(prm_diag%swflxsfc, p_patch%cells%area, p_patch%cells%owned, outvar0d)
     CASE('lwf_sfc')
       IF(atm_phy_nwp_config(jg)%inwp_radiation>0) &
       CALL levels_horizontal_mean(prm_diag%lwflxsfc, p_patch%cells%area, p_patch%cells%owned, outvar0d)
     CASE('precp_t')
       CALL levels_horizontal_mean(prm_diag%tot_prec_rate_avg, p_patch%cells%area, p_patch%cells%owned, outvar0d)
       outvar0d = outvar0d * day_sec
     CASE('precp_r')
       IF(atm_phy_nwp_config(jg)%inwp_gscp>0) &
       CALL levels_horizontal_mean(prm_diag%rain_gsp_rate, p_patch%cells%area, p_patch%cells%owned, outvar0d)
       outvar0d = outvar0d * day_sec
     CASE('precp_s')
       IF(atm_phy_nwp_config(jg)%inwp_gscp>0) &
       CALL levels_horizontal_mean(prm_diag%snow_gsp_rate, p_patch%cells%area, p_patch%cells%owned, outvar0d)
       outvar0d = outvar0d * day_sec
     CASE('precp_g')
       IF(ANY((/4,5/) == atm_phy_nwp_config(jg)%inwp_gscp))THEN
         CALL levels_horizontal_mean(prm_diag%graupel_gsp_rate, p_patch%cells%area, p_patch%cells%owned, outvar0d)
         outvar0d = outvar0d * day_sec
       END IF
     CASE('precp_h')
       IF(ANY((/4,5/) == atm_phy_nwp_config(jg)%inwp_gscp))THEN
         CALL levels_horizontal_mean(prm_diag%hail_gsp_rate, p_patch%cells%area, p_patch%cells%owned, outvar0d)
         outvar0d = outvar0d * day_sec
       END IF
     CASE('precp_i')
       IF(ANY((/4,5/) == atm_phy_nwp_config(jg)%inwp_gscp))THEN
         CALL levels_horizontal_mean(prm_diag%ice_gsp_rate, p_patch%cells%area, p_patch%cells%owned, outvar0d)
         outvar0d = outvar0d * day_sec
       END IF
     END SELECT  

     prm_diag%turb_diag_0dvar(n) = outvar0d 
 
    END DO
 
    DEALLOCATE( umean, vmean, thmean, qvmean, qcmean, wmean, outvar, var3df, var3dh, theta, w_mc )

    IF(msg_level>18) & 
      CALL message(routine,'Over!')

  END SUBROUTINE calculate_turbulent_diagnostics

!===========================================================================
!>
  !! <write out profile>
  !!
  !! <Describe the purpose of the subroutine and its algorithm(s).>
  !! <Include any applicable external references inline as module::procedure,>
  !! <external_procedure(), or by using @see.>
  !! <Don't forget references to literature.>
  !!
  !! @par Revision History
  !!
  SUBROUTINE write_vertical_profiles(outvar, this_datetime, ncount)
    REAL(wp),                INTENT(IN)  :: outvar(:,:)
    TYPE(datetime), POINTER, INTENT(IN)  :: this_datetime
    INTEGER,                 INTENT(IN)  :: ncount

    INTEGER                  :: nvar, n
    REAL(wp)                 :: inv_ncount
    REAL(wp)                 :: sim_time     !< elapsed simulation time on this grid level
    TYPE(timeDelta), POINTER :: time_diff

    
    ! calculate elapsed simulation time in seconds
    time_diff  => newTimedelta("PT0S")
    time_diff  =  getTimeDeltaFromDateTime(this_datetime, time_config%tc_startdate)
    sim_time   =  getTotalMillisecondsTimedelta(time_diff, this_datetime)*1.e-3_wp
    CALL deallocateTimedelta(time_diff)
 
    !Write profiles
      
    inv_ncount = 1._wp / REAL(ncount,wp)

    nvar = SIZE(turb_profile_list,1)

    !Loop over all variables
    IF( my_process_is_stdio() )THEN

      !First write time
      CALL writevar_nc(fileid_profile, tname, sim_time, nrec_profile) 

      DO n = 1 , nvar       
       CALL writevar_nc(fileid_profile, TRIM(turb_profile_list(n)),  &
                        outvar(:,n)*inv_ncount, nrec_profile) 
      END DO

    END IF

  END SUBROUTINE write_vertical_profiles

!===========================================================================
!>
  !! <write out time series>
  !!
  !! <Describe the purpose of the subroutine and its algorithm(s).>
  !! <Include any applicable external references inline as module::procedure,>
  !! <external_procedure(), or by using @see.>
  !! <Don't forget references to literature.>
  !!
  !! @par Revision History
  !!
  SUBROUTINE write_time_series(outvar, this_datetime)
    REAL(wp),                INTENT(IN) :: outvar(:)
    TYPE(datetime), POINTER, INTENT(IN) :: this_datetime

    INTEGER                  :: nvar, n
    REAL(wp)                 :: sim_time     !< elapsed simulation time on this grid level
    TYPE(timeDelta), POINTER :: time_diff
    
    ! calculate elapsed simulation time in seconds
    time_diff  => newTimedelta("PT0S")
    time_diff  =  getTimeDeltaFromDateTime(this_datetime, time_config%tc_startdate)
    sim_time   =  getTotalMillisecondsTimedelta(time_diff, this_datetime)*1.e-3_wp
    CALL deallocateTimedelta(time_diff)

    !Write time series 
    nvar = SIZE(turb_tseries_list,1)

    !Loop over all variables
    IF( my_process_is_stdio() )THEN

      !First write time
      CALL writevar_nc(fileid_tseries, tname, sim_time, nrec_tseries) 

      DO n = 1 , nvar
        CALL writevar_nc(fileid_tseries, turb_tseries_list(n), outvar(n), nrec_tseries) 
      END DO

    END IF

  END SUBROUTINE write_time_series

!===========================================================================
!>
  !! <initialize turbulent output>
  !!
  !! <Describe the purpose of the subroutine and its algorithm(s).>
  !! <Include any applicable external references inline as module::procedure,>
  !! <external_procedure(), or by using @see.>
  !! <Don't forget references to literature.>
  !!
  !! @par Revision History
  !!
  SUBROUTINE init_les_turbulent_output(p_patch, p_metrics, this_datetime, l_rh, ldelete)
   TYPE(t_patch),   TARGET, INTENT(in)   :: p_patch    !<grid/patch info.
   TYPE(t_nh_metrics)     , INTENT(in)   :: p_metrics
   TYPE(datetime), POINTER, INTENT(IN)   :: this_datetime
   LOGICAL, INTENT(IN), OPTIONAL         :: ldelete
   LOGICAL, INTENT(IN), OPTIONAL         :: l_rh  !if rh to be output or not
  
   CHARACTER (40), ALLOCATABLE, DIMENSION(:) :: dimname, dimlongname, dimunit
   CHARACTER (LEN=80)                        :: longname, unit
   REAL(wp), ALLOCATABLE                     :: dimvalues(:,:)
   INTEGER,  ALLOCATABLE                     :: dimsize(:)
   INTEGER :: n, nlev, nlevp1, nvar, jg
   REAL(wp) :: z_mc_avg(p_patch%nlev), z_ic_avg(p_patch%nlev+1)
   CHARACTER(len=*), PARAMETER :: routine = 'mo_turbulent_diagnostic:init_les_turbulent_output'
   TYPE(timeDelta), POINTER            :: time_diff
   REAL(wp)                            :: p_sim_time     !< elapsed simulation time on this grid level
 
   ! calculate elapsed simulation time in seconds
   time_diff  => newTimedelta("PT0S")
   time_diff  =  getTimeDeltaFromDateTime(this_datetime, time_config%tc_startdate)
   p_sim_time =  getTotalMillisecondsTimedelta(time_diff, this_datetime)*1.e-3_wp
   CALL deallocateTimedelta(time_diff)

   jg = p_patch%id

   !Check if diagnostics are to be calculated or not
   IF(.NOT.les_config(jg)%ldiag_les_out)THEN
    sampl_freq_step   = 0
    avg_interval_step = 0
    RETURN
   END IF

   IF(msg_level>18)CALL message(routine,'Start!')

   !Sampling and output frequencies in terms of time steps
   sampl_freq_step   = NINT(les_config(jg)%sampl_freq_sec/dtime)
   avg_interval_step = NINT(les_config(jg)%avg_interval_sec/dtime)

   is_rh_out = l_rh

   nlev   = p_patch%nlev
   nlevp1 = nlev + 1

   !Dimensions
   ALLOCATE( dimname(2), dimlongname(2), dimunit(2), dimsize(2), dimvalues(nlevp1,2) )

   !Calculate average height
   CALL levels_horizontal_mean(p_metrics%z_mc, p_patch%cells%area, p_patch%cells%owned,  z_mc_avg)
   CALL levels_horizontal_mean(p_metrics%z_ifc, p_patch%cells%area, p_patch%cells%owned, z_ic_avg)

   !open profile file
   IF( my_process_is_stdio() ) &
     CALL open_nc(TRIM(les_config(jg)%expname)//'_profile.nc', fileid_profile, nrec_profile, p_sim_time, ldelete)
 
   !addvar
   nvar = SIZE(turb_profile_list,1)

   ALLOCATE(is_at_full_level(nvar))

   DO n = 1 , nvar

     is_at_full_level(n) = .TRUE.

     SELECT CASE (TRIM(turb_profile_list(n)))
    
     CASE('u')
      longname = 'zonal wind'
      unit     = 'm/s'
     CASE('v')
      longname = 'meridional wind'
      unit     = 'm/s'
     CASE('w')
      longname = 'vertical wind'
      unit     = 'm/s'
     CASE('th') !theta mean
      longname = 'potential temperature'
      unit     = 'K'
     CASE('thv') !thetav mean
      longname = 'virtual potential temperature'
      unit     = 'K'
     CASE('exner')
      longname = 'exner function'
      unit     = ''
     CASE('rho')
      longname = 'density'
      unit     = 'kg/m3'
     CASE('qv')
      longname = 'specific humidity'
      unit     = 'kg/kg'
     CASE('qc')
      longname = 'cloud water'
      unit     = 'kg/kg'
     CASE('wu')
       longname = 'resolved zonal wind flux'
       unit     = 'm2/s2'
     CASE('wv')
       longname = 'resolved meridional wind flux'
       unit     = 'm2/s2'
     CASE('wth')
       longname = 'resolved potential temperature flux'
       unit     = 'W/m2'
     CASE('wthv')
       longname = 'resolved virtual potential temperature flux'
       unit     = 'W/m2'
     CASE('wqv')
       longname = 'resolved specific humidity flux'
       unit     = 'W/m2'
     CASE('wqc')
       longname = 'resolved cloud water flux'
       unit     = 'W/m2'
     CASE('ww')
       longname = 'resolved vertical velocity variance'
       unit     = 'm2/s2'
     CASE('thth')
       longname = 'resolved potential temperature variance'
       unit     = 'K2'
     CASE('qvqv')
       longname = 'resolved specific humidity variance'
       unit     = 'kg2/kg2'
     CASE('qcqc')
       longname = 'resolved cloud water variance'
       unit     = 'kg2/kg2'
     CASE('uu')
       longname = 'resolved zonal wind variance'
       unit     = 'm2/s2'
     CASE('vv')
       longname = 'resolved meridional wind variance'
       unit     = 'm2/s2'
     CASE('kh')
       longname = 'mass weighted eddy diffusivity'
       unit     = 'kg/(ms)'
       is_at_full_level(n) = .FALSE.
     CASE('km')
       longname = 'mass weighted eddy viscosity'
       unit     = 'kg/(ms)'
       is_at_full_level(n) = .FALSE.
     CASE('wud') !diffuse u flux
       longname = '(mass) subgrid zonal wind flux'
       unit     = 'kg/ms2'
       is_at_full_level(n) = .FALSE.
       idx_sgs_u_flx = n
     CASE('wvd') !diffuse flux
       longname = '(mass) subgrid meridional wind flux'
       unit     = 'kg/ms2'
       is_at_full_level(n) = .FALSE.
       idx_sgs_v_flx = n
     CASE('wthd') !diffuse flux
       longname = 'subgrid potential temperature flux'
       unit     = 'W/m2'
       is_at_full_level(n) = .FALSE.
       idx_sgs_th_flx = n
     CASE('wqvd') !diffuse flux
       longname = 'subgrid specific humidity flux'
       unit     = 'W/m2'
       is_at_full_level(n) = .FALSE.
       idx_sgs_qv_flx = n
     CASE('wqcd') !diffuse flux
       longname = 'subgrid cloud water flux'
       unit     = 'W/m2'
       is_at_full_level(n) = .FALSE.
       idx_sgs_qc_flx = n
     CASE('bruvais') 
       longname = 'Brunt Vaisala Frequency'
       unit     = '1/s2'
       is_at_full_level(n) = .FALSE.
     CASE('mechprd') 
       longname = 'Mechanical production term in TKE'
       unit     = '1/s2'
       is_at_full_level(n) = .FALSE.
     CASE('wthsfs') 
       longname = 'sub-filter scale flux'
       unit     = 'W/m2'
     CASE('rh') 
       longname = 'relative humidity'
       unit     = '%'
     CASE('clc') 
       longname = 'cloud cover'
       unit     = '-'
     CASE('qi') 
       longname = 'specific ice content'
       unit     = 'kg/kg'
     CASE('qs') 
       longname = 'specific snow content'
       unit     = 'kg/kg'
     CASE('qr') 
       longname = 'specific rain content'
       unit     = 'kg/kg'
     CASE('qg') 
       longname = 'specific graupel content'
       unit     = 'kg/kg'
     CASE('qh') 
       longname = 'specific hail content'
       unit     = 'kg/kg'
     CASE('lwf') 
       longname = 'net longwave flux'
       unit     = 'W/m2'
       is_at_full_level(n) = .FALSE.
     CASE('swf') 
       longname = 'net shortwave flux'
       unit     = 'W/m2'
       is_at_full_level(n) = .FALSE.
     CASE('dt_t_sw') 
       longname = 'shortwave temp tendency'
       unit     = 'K/s'
     CASE('dt_t_lw') 
       longname = 'longwave temp tendency'
       unit     = 'K/s'
     CASE('dt_t_tb') 
       longname = 'turbulent temp tendency'
       unit     = 'K/s'
     CASE('dt_t_mc') 
       longname = 'microphysics temp tendency'
       unit     = 'K/s'
     CASE DEFAULT 
         WRITE(message_text,'(a)')TRIM(turb_profile_list(n))
         CALL finish(routine,'Variable '//TRIM(message_text)//' is not listed in les_nml')
     END SELECT

     dimname(2) = tname
     dimlongname(2) = tlongname
     dimunit(1) = 'm'
     dimunit(2) = 's'
     dimvalues = 0._wp
     IF(is_at_full_level(n))THEN
      dimname(1) = 'zf'
      dimlongname(1) = 'Full level height'
      dimsize = (/nlev,0/)
      dimvalues(1:nlev,1) = z_mc_avg(1:nlev)     
     ELSE
      dimname(1) = 'zh'
      dimlongname(1) = 'Half level height'
      dimsize = (/nlevp1,0/)
      dimvalues(1:nlevp1,1) = z_ic_avg(1:nlevp1)     
     END IF

     IF( my_process_is_stdio() ) &
       CALL addvar_nc(fileid_profile, TRIM(turb_profile_list(n)), TRIM(longname), TRIM(unit), &
                      dimname, dimlongname, dimunit, dimsize, dimvalues)

    END DO!nvar
    

    !deallocate
    DEALLOCATE( dimname, dimlongname, dimunit, dimsize, dimvalues )
    ALLOCATE( dimname(1), dimlongname(1), dimunit(1) )


   !open time series file
   IF( my_process_is_stdio() ) &
      CALL open_nc(TRIM(les_config(jg)%expname)//'_tseries.nc', fileid_tseries, nrec_tseries, p_sim_time, ldelete)
 
   !addvar
   nvar = SIZE(turb_tseries_list,1)

   DO n = 1 , nvar

     SELECT CASE (TRIM(turb_tseries_list(n)))
 
     CASE('ccover')
       longname = 'cloud cover'
       unit     = ' '
     CASE('shflx')
       longname = 'surface sensible heat flux'
       unit     = 'W/m2'
     CASE('lhflx')
       longname = 'surface latent heat flux'
       unit     = 'W/m2'
     CASE('ustress')
       longname = 'surface zonal stress'
       unit     = 'Kg/ms2'
     CASE('vstress')
       longname = 'surface meridional stress'
       unit     = 'Kg/ms2'
     CASE('tsfc')
       longname = 'surface temperature'
       unit     = 'K'
     CASE('qsfc')
       longname = 'surface humidity'
       unit     = 'kg/kg'
     CASE('hbl')
       longname = 'boundary layer height'
       unit     = 'm'
     CASE('psfc')
       longname = 'surface pressure'
       unit     = 'Pa'
     CASE('swf_tom')
       longname = 'net shortwave flux at model top'
       unit     = 'W/m2'
     CASE('lwf_tom')
       longname = 'net longwave flux at model top'
       unit     = 'W/m2'
     CASE('swf_sfc')
       longname = 'net shortwave flux at sfc'
       unit     = 'W/m2'
     CASE('lwf_sfc')
       longname = 'net longwave flux at sfc'
       unit     = 'W/m2'
     CASE('precp_t')
       longname = 'time avg total gridscale precipitation'
       unit     = 'mm/day'
     CASE('precp_r')
       longname = 'gridscale rain rate'
       unit     = 'mm/day'
     CASE('precp_s')
       longname = 'gridscale snow rate'
       unit     = 'mm/day'
     CASE('precp_g')
       longname = 'gridscale graupel rate'
       unit     = 'mm/day'
     CASE('precp_h')
       longname = 'gridscale hail rate'
       unit     = 'mm/day'
     CASE('precp_i')
       longname = 'gridscale ice rate'
       unit     = 'mm/day'
     CASE DEFAULT 
         WRITE(message_text,'(a)')TRIM(turb_tseries_list(n))
         CALL finish(routine,'Variable '//TRIM(message_text)//' is not listed in les_nml')
     END SELECT  
  
     dimname(1) = tname
     dimlongname(1) = tlongname
     dimunit(1) = 's'
     IF( my_process_is_stdio() ) &
        CALL addvar_nc(fileid_tseries, TRIM(turb_tseries_list(n)), TRIM(longname), TRIM(unit), &
                       dimname, dimlongname, dimunit)

   END DO!nvar

   DEALLOCATE( dimname, dimlongname, dimunit )

   IF(msg_level>18)CALL message(routine,'Over!')

  END SUBROUTINE init_les_turbulent_output

!===========================================================================
!>
  !! <close turbulent output>
  !!
  !! <Describe the purpose of the subroutine and its algorithm(s).>
  !! <Include any applicable external references inline as module::procedure,>
  !! <external_procedure(), or by using @see.>
  !! <Don't forget references to literature.>
  !!
  !! @par Revision History
  !!
  SUBROUTINE close_les_turbulent_output(jg)
   INTEGER, INTENT(IN) :: jg

   IF(.NOT.les_config(jg)%ldiag_les_out)THEN
    RETURN
   END IF

   IF( my_process_is_stdio() ) THEN
     CALL close_nc(fileid_profile) 
     CALL close_nc(fileid_tseries) 
   END IF
   DEALLOCATE(is_at_full_level)

  END SUBROUTINE close_les_turbulent_output

!===========================================================================

END MODULE mo_turbulent_diagnostic

