!>
!! This module is the interface between nwp_nh_interface to the 
!! surface parameterisations:
!! inwp_sfc  == 1 == surface scheme TERRA run in COSMO
!!
!! @author Kristina Froehlich, DWD, Offenbach (2010-01-25)
!!
!! @par Revision History
!! Initial Kristina Froehlich, DWD, Offenbach (2010-01-25)
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

MODULE mo_nwp_sfc_interface

  USE mo_kind,                ONLY: wp
  USE mo_exception,           ONLY: message, message_text, finish
  USE mo_model_domain,        ONLY: t_patch
  USE mo_impl_constants,      ONLY: min_rlcell_int, iedmf, icosmo
  USE mo_impl_constants_grf,  ONLY: grf_bdywidth_c
  USE mo_loopindices,         ONLY: get_indices_c
  USE mo_ext_data_types,      ONLY: t_external_data
  USE mo_nonhydro_types,      ONLY: t_nh_prog, t_nh_diag, t_nh_metrics
  USE mo_nwp_lnd_types,       ONLY: t_lnd_prog, t_wtr_prog, t_lnd_diag
  USE mo_nwp_phy_types,       ONLY: t_nwp_phy_diag
  USE mo_nwp_phy_state,       ONLY: phy_params
  USE mo_parallel_config,     ONLY: nproma
  USE mo_run_config,          ONLY: iqv, iqi, msg_level
  USE mo_atm_phy_nwp_config,  ONLY: atm_phy_nwp_config
  USE mo_lnd_nwp_config,      ONLY: nlev_soil, nlev_snow, ibot_w_so, ntiles_total,    &
    &                               ntiles_water, lseaice, llake, lmulti_snow,        &
    &                               ntiles_lnd, lsnowtile, isub_water, isub_seaice,   &
    &                               isub_lake, itype_interception, l2lay_rho_snow,    &
    &                               lprog_albsi, itype_trvg, itype_snowevap, zml_soil
  USe mo_extpar_config,       ONLY: itype_vegetation_cycle
  USE mo_ensemble_pert_config,ONLY: sst_pert_corrfac
  USE mo_satad,               ONLY: sat_pres_water, sat_pres_ice, spec_humi, dqsatdT_ice
  USE sfc_terra,              ONLY: terra
  USE mo_nwp_sfc_utils,       ONLY: diag_snowfrac_tg, update_idx_lists_lnd, update_idx_lists_sea
  USE sfc_flake,              ONLY: flake_interface
  USE sfc_flake_data,         ONLY: h_Ice_min_flk
  USE sfc_seaice,             ONLY: seaice_timestep_nwp
  USE sfc_terra_data                ! soil and vegetation parameters for TILES
  USE mo_physical_constants,  ONLY: tmelt, grav, salinity_fac, rhoh2o
  USE mo_nwp_gpu_util,        ONLY: gpu_d2h_nh_nwp, gpu_h2d_nh_nwp
#ifdef COUP_OAS_ICON
  USE oas_icon_define
  USE mo_physical_constants,  ONLY: lh_v => alv, lh_s => als, t0_melt =>tmelt,        &
                                    r_d   => rd, rvd_m_o=>vtmpc1 ! SPo add rd & r_v/r_d - 1 
  !USE mo_nwp_lnd_state,       ONLY: p_lnd_state
  !USE mo_dynamics_config,     ONLY: nnow
  USE oas_icon_define,        ONLY: oas_rcv_field_icon
#endif

  IMPLICIT NONE 

  PRIVATE



  PUBLIC  ::  nwp_surface


#ifdef __SX__
! parameter for loop unrolling
INTEGER, PARAMETER :: nlsoil= 8
#endif


CONTAINS
  !!
  !!-------------------------------------------------------------------------
  !!
  SUBROUTINE nwp_surface( tcall_sfc_jg,                   & !>in
                        & p_patch,                        & !>in
                        & ext_data,                       & !>in
                        & p_prog, p_prog_rcf,             & !>in/inout
                        & p_diag, p_metrics,              & !>inout
                        & prm_diag,                       & !>inout 
                        & lnd_prog_now, lnd_prog_new,     & !>inout
                        & p_prog_wtr_now, p_prog_wtr_new, & !>inout
                        & lnd_diag,                       & !>inout
                        & lacc                            ) !>in

    TYPE(t_patch),        TARGET,INTENT(in)   :: p_patch       !< grid/patch info
    TYPE(t_external_data),       INTENT(inout):: ext_data      !< external data
    TYPE(t_nh_prog),      TARGET,INTENT(inout):: p_prog        !< dynamic prognostic vars
    TYPE(t_nh_prog),      TARGET,INTENT(inout):: p_prog_rcf    !< call freq
    TYPE(t_nh_diag),      TARGET,INTENT(inout):: p_diag        !< diag vars
    TYPE(t_nh_metrics),   TARGET,INTENT(in)   :: p_metrics     !< metrics vars
    TYPE(t_nwp_phy_diag),        INTENT(inout):: prm_diag      !< atm phys vars
    TYPE(t_lnd_prog),            INTENT(inout):: lnd_prog_now  !< prog vars for sfc
    TYPE(t_lnd_prog),            INTENT(inout):: lnd_prog_new  !< prog vars for sfc
    TYPE(t_wtr_prog),            INTENT(inout):: p_prog_wtr_now !< prog vars for wtr
    TYPE(t_wtr_prog),            INTENT(inout):: p_prog_wtr_new !< prog vars for wtr
    TYPE(t_lnd_diag),            INTENT(inout):: lnd_diag      !< diag vars for sfc
    REAL(wp),                    INTENT(in)   :: tcall_sfc_jg  !< time interval for 
    LOGICAL, OPTIONAL,           INTENT(in)   :: lacc          !< GPU flag
    LOGICAL :: lzacc

    ! Local array bounds:
    !
    INTEGER :: rl_start, rl_end
    INTEGER :: i_startblk, i_endblk    !> blocks
    INTEGER :: i_startidx, i_endidx    !< slices
    INTEGER :: i_nchdom                !< domain index
    INTEGER :: nlev                    !< number of full levels
    INTEGER :: isubs, isubs_snow, icant

    ! Local scalars:
    !
    INTEGER  :: jc,jb,jg,jk     ! loop indices
    REAL(wp) :: area_frac       ! tile area fraction

    REAL(wp) :: ps_t        (nproma)
    REAL(wp) :: prr_con_t   (nproma)
    REAL(wp) :: prs_con_t   (nproma)
    REAL(wp) :: conv_frac   (nproma)
    REAL(wp) :: prr_gsp_t   (nproma)
    REAL(wp) :: prs_gsp_t   (nproma)
    REAL(wp) :: prg_gsp_t   (nproma)

    REAL(wp) :: u_t (nproma)
    REAL(wp) :: v_t (nproma)
    REAL(wp) :: t_t (nproma)
    REAL(wp) :: qv_t(nproma)
    REAL(wp) :: p0_t(nproma)

    REAL(wp) :: sso_sigma_t(nproma)
    INTEGER  :: lc_class_t (nproma)

    REAL(wp) :: t_snow_now_t (nproma)
    REAL(wp) :: t_snow_new_t (nproma)

    REAL(wp) :: t_s_now_t  (nproma)
    REAL(wp) :: t_s_new_t  (nproma)

    REAL(wp) :: t_sk_now_t (nproma)
    REAL(wp) :: t_sk_new_t (nproma)

    REAL(wp) :: t_g_t      (nproma)
    REAL(wp) :: qv_s_t     (nproma)

    REAL(wp) :: w_snow_now_t(nproma)
    REAL(wp) :: w_snow_new_t(nproma)
  
    REAL(wp) :: rho_snow_now_t (nproma)
    REAL(wp) :: rho_snow_new_t (nproma)

    REAL(wp) :: h_snow_t (nproma)
    REAL(wp) :: h_snow_gp_t (nproma)
    REAL(wp) :: meltrate (nproma)

    REAL(wp) :: w_i_now_t (nproma)
    REAL(wp) :: w_i_new_t (nproma)

    REAL(wp) :: w_p_now_t (nproma)
    REAL(wp) :: w_p_new_t (nproma)

    REAL(wp) :: w_s_now_t (nproma)
    REAL(wp) :: w_s_new_t (nproma)

    REAL(wp) :: u_10m_t    (nproma)
    REAL(wp) :: v_10m_t    (nproma)
    REAL(wp) :: freshsnow_t(nproma)
    REAL(wp) :: snowfrac_t (nproma)
    REAL(wp) :: snowfrac_lcu_t (nproma)

    REAL(wp) :: tch_t      (nproma)
    REAL(wp) :: tcm_t      (nproma)
    REAL(wp) :: tfv_t      (nproma)

    REAL(wp) :: sobs_t     (nproma)
    REAL(wp) :: thbs_t     (nproma)
    REAL(wp) :: pabs_t     (nproma)

    REAL(wp) :: runoff_s_inst_t (nproma)
    REAL(wp) :: runoff_g_inst_t (nproma)

    INTEGER  :: soiltyp_t (nproma)
    REAL(wp) :: plcov_t   (nproma)
    REAL(wp) :: rootdp_t  (nproma)
    REAL(wp) :: sai_t     (nproma)
    REAL(wp) :: tai_t     (nproma)
    REAL(wp) :: laifac_t  (nproma)
    REAL(wp) :: eai_t     (nproma)
    REAL(wp) :: skinc_t   (nproma)
    REAL(wp) :: rsmin2d_t (nproma)

    ! local dummy variable for precipitation rate of graupel, grid-scale
    REAL(wp), TARGET  :: dummy_graupel_gsp_rate(nproma,p_patch%nblks_c)
    ! pointer to actual or dummy precipitation rate of graupel
    REAL(wp), POINTER :: p_graupel_gsp_rate(:,:)

    REAL(wp) :: t_snow_mult_now_t(nproma, nlev_snow+1)
    REAL(wp) :: t_snow_mult_new_t(nproma, nlev_snow+1)

    REAL(wp) :: rho_snow_mult_now_t(nproma, nlev_snow)
    REAL(wp) :: rho_snow_mult_new_t(nproma, nlev_snow)

    REAL(wp) :: wliq_snow_now_t(nproma, nlev_snow)
    REAL(wp) :: wliq_snow_new_t(nproma, nlev_snow)

    REAL(wp) :: wtot_snow_now_t(nproma, nlev_snow)
    REAL(wp) :: wtot_snow_new_t(nproma, nlev_snow)

    REAL(wp) :: dzh_snow_now_t(nproma, nlev_snow)
    REAL(wp) :: dzh_snow_new_t(nproma, nlev_snow)

    REAL(wp) :: t_so_now_t(nproma, nlev_soil+1)
    REAL(wp) :: t_so_new_t(nproma, nlev_soil+1)

    REAL(wp) :: w_so_now_t(nproma, nlev_soil)
    REAL(wp) :: w_so_new_t(nproma, nlev_soil)

    REAL(wp) :: w_so_ice_now_t(nproma, nlev_soil)
    REAL(wp) :: w_so_ice_new_t(nproma, nlev_soil)

    INTEGER  :: i_count, i_count_snow, ic, icount_init, is1, is2, init_list(nproma), it1(nproma), it2(nproma)
    REAL(wp) :: tmp1, tmp2, tmp3, qsat1, dqsdt1, qsat2, dqsdt2
    REAL(wp) :: frac_sv(nproma), frac_snow_sv(nproma), fact1(nproma), fact2(nproma), tsnred(nproma), &
                sntunefac(nproma), sntunefac2(nproma, ntiles_total)
    REAL(wp) :: rain_gsp_rate(nproma, ntiles_total)
    REAL(wp) :: snow_gsp_rate(nproma, ntiles_total)
    REAL(wp) :: rain_con_rate(nproma, ntiles_total)
    REAL(wp) :: snow_con_rate(nproma, ntiles_total)
    REAL(wp) :: graupel_gsp_rate(nproma, ntiles_total)
    REAL(wp), PARAMETER :: small = 1.E-06_wp

    REAL(wp) :: t_g_s(nproma)
    REAL(wp) :: shfl_s_t    (nproma) ! sensible heat flux sfc
    REAL(wp) :: lhfl_s_t    (nproma) ! latent heat flux sfc
    REAL(wp) :: qhfl_s_t    (nproma) ! moisture flux sfc
    REAL(wp) :: shfl_soil_t (nproma) ! sensible heat flux sfc (snow free)
    REAL(wp) :: lhfl_soil_t (nproma) ! latent heat flux sfc   (snow free)
    REAL(wp) :: shfl_snow_t (nproma) ! sensible heat flux sfc (snow covered)
    REAL(wp) :: lhfl_snow_t (nproma) ! latent heat flux sfc   (snow covered)
    REAL(wp) :: lhfl_bs_t   (nproma)
    REAL(wp) :: lhfl_pl_t   (nproma, nlev_soil)
    REAL(wp) :: plevap_t    (nproma)
    REAL(wp) :: rstom_t     (nproma)
    REAL(wp) :: z0_t        (nproma)

#ifdef COUP_OAS_ICON
    REAL(wp) :: qhfl_help(nproma,p_patch%nblks_c)
    REAL(wp) :: qv_s_help(nproma,p_patch%nblks_c)
#endif

    CHARACTER(len=*), PARAMETER :: routine = 'mo_nwp_sfc_interface:nwp_surface'

!--------------------------------------------------------------
    IF(PRESENT(lacc)) THEN
        lzacc = lacc
    ELSE
        lzacc = .FALSE.
    ENDIF

    ! get patch ID
    jg = p_patch%id


    IF (atm_phy_nwp_config(jg)%lhave_graupel) THEN
      ! COSMO-DE (3-cat ice: snow, cloud ice, graupel)
      p_graupel_gsp_rate => prm_diag%graupel_gsp_rate(:,:)
    ELSE
      ! initialize dummy variable (precipitation rate of graupel, grid-scale)
      dummy_graupel_gsp_rate(:,:) = 0._wp
      p_graupel_gsp_rate => dummy_graupel_gsp_rate(:,:)
    ENDIF


    ! local variables related to the blocking

    i_nchdom  = MAX(1,p_patch%n_childdom)

    ! number of vertical levels
    nlev   = p_patch%nlev

    ! exclude nest boundary and halo points
    rl_start = grf_bdywidth_c+1
    rl_end   = min_rlcell_int

    i_startblk = p_patch%cells%start_blk(rl_start,1)
    i_endblk   = p_patch%cells%end_blk(rl_end,i_nchdom)

    ! canopy-type needed by TERRA:
    SELECT CASE (atm_phy_nwp_config(jg)%inwp_turb)
    CASE(icosmo,iedmf)
       icant=2 !canopy-treatment related to Raschendorfer-transfer-scheme
    CASE DEFAULT
       icant=1 !canopy-treatment related to Louis-transfer-scheme
    END SELECT

    IF (msg_level >= 15) THEN
      CALL message('mo_nwp_sfc_interface: ', 'call land-surface scheme')
    ENDIF

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,jk,i_startidx,i_endidx,isubs,i_count,ic,isubs_snow,i_count_snow,                     &
!$OMP   tmp1,tmp2,tmp3,fact1,fact2,frac_sv,frac_snow_sv,icount_init,init_list,it1,it2,is1,is2,              &
!$OMP   rain_gsp_rate,snow_gsp_rate,rain_con_rate,snow_con_rate,ps_t,prr_con_t,prs_con_t,                   &
!$OMP   prr_gsp_t,prs_gsp_t,u_t,v_t,t_t,qv_t,p0_t,sso_sigma_t,lc_class_t,t_snow_now_t,t_s_now_t,            &
!$OMP   t_g_t,qv_s_t,w_snow_now_t,rho_snow_now_t,w_i_now_t,w_p_now_t,w_s_now_t,freshsnow_t,                 &
!$OMP   snowfrac_t,runoff_s_inst_t,runoff_g_inst_t,u_10m_t,v_10m_t,tch_t,tcm_t,tfv_t,sobs_t,thbs_t,pabs_t,  &
!$OMP   soiltyp_t,plcov_t,rootdp_t,sai_t,tai_t,eai_t,rsmin2d_t,t_snow_mult_now_t,wliq_snow_now_t,           &
!$OMP   rho_snow_mult_now_t,wtot_snow_now_t,dzh_snow_now_t,t_so_now_t,w_so_now_t,w_so_ice_now_t,            &
!$OMP   t_s_new_t,w_snow_new_t,rho_snow_new_t,h_snow_t,w_i_new_t,w_p_new_t,w_s_new_t,t_so_new_t,            &
!$OMP   lhfl_bs_t,rstom_t,shfl_s_t,lhfl_s_t,qhfl_s_t,t_snow_mult_new_t,rho_snow_mult_new_t,                 &
!$OMP   wliq_snow_new_t,wtot_snow_new_t,dzh_snow_new_t,w_so_new_t,w_so_ice_new_t,lhfl_pl_t,                 &
!$OMP   shfl_soil_t,lhfl_soil_t,shfl_snow_t,lhfl_snow_t,t_snow_new_t,graupel_gsp_rate,prg_gsp_t,            &
!$OMP   meltrate,h_snow_gp_t,conv_frac,t_sk_now_t,t_sk_new_t,skinc_t,tsnred,plevap_t,z0_t,laifac_t,         &
!$OMP   qsat1,dqsdt1,qsat2,dqsdt2,sntunefac,sntunefac2,snowfrac_lcu_t) ICON_OMP_GUIDED_SCHEDULE

    !$acc enter data copyin (zml_soil)                                             &
    !$acc            create (sntunefac, sntunefac2, rain_con_rate, snow_con_rate)  &
    !$acc            create (rain_gsp_rate, snow_gsp_rate, graupel_gsp_rate)       &
    !$acc            create (init_list, it1, it2, fact1, fact2, frac_sv)           &
    !$acc            create (frac_snow_sv),                                        &
    !$acc if(lzacc)

    DO jb = i_startblk, i_endblk

      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
        & i_startidx, i_endidx, rl_start, rl_end)
      IF( atm_phy_nwp_config(jg)%inwp_surface == 0) THEN
        ! check dry case
        IF( atm_phy_nwp_config(jg)%inwp_satad == 0) THEN
          lnd_diag%qv_s (:,jb) = 0._wp
        ELSE
          ! 
          !> adjust humidity at water surface because of changed surface pressure
          !
          !$acc parallel if(lzacc)
          !$acc loop gang vector
          DO jc = i_startidx, i_endidx
            lnd_diag%qv_s (jc,jb) = &
              &            spec_humi(sat_pres_water(lnd_prog_now%t_g(jc,jb)),&
              &                                      p_diag%pres_sfc(jc,jb) )
          ENDDO
          !$acc end parallel
        ENDIF
      ELSE  ! inwp_surface/=0
         ! 
         !> adjust humidity at water surface because of changing surface pressure
         !
!$NEC ivdep
         !$acc parallel if(lzacc)
         !$acc loop gang vector private(jc)
         DO ic=1,ext_data%atm%list_seawtr%ncount(jb)
           jc = ext_data%atm%list_seawtr%idx(ic,jb)

           ! salinity_fac accounts for the average reduction of saturation pressure caused by the salt content of oceans
           ! sst_pert_corrfac is a tuning factor to compensate the increased evaporation due to SST ensemble perturbations
           lnd_diag%qv_s_t(jc,jb,isub_water) = salinity_fac * sst_pert_corrfac *      &
             &         spec_humi(sat_pres_water(lnd_prog_now%t_g_t(jc,jb,isub_water)),&
             &                                   p_diag%pres_sfc(jc,jb) )
         ENDDO
         !$acc end parallel
      ENDIF


      IF ( atm_phy_nwp_config(jg)%inwp_surface == 1 ) THEN

       IF (ext_data%atm%list_land%ncount(jb) == 0) CYCLE ! skip loop if there is no land point

       ! Copy precipitation fields for subsequent downscaling
       !$acc parallel if(lzacc)
       !$acc loop gang vector private(i_count)
       DO isubs = 1,ntiles_total
         i_count = ext_data%atm%gp_count_t(jb,isubs) 
         IF (i_count == 0) CYCLE ! skip loop if the index list for the given tile is empty
!$NEC ivdep
         !$acc loop private(jc)
         DO ic = 1, i_count
           jc = ext_data%atm%idx_lst_t(ic,jb,isubs)
           rain_gsp_rate(jc,isubs)    = prm_diag%rain_gsp_rate(jc,jb)
           snow_gsp_rate(jc,isubs)    = prm_diag%snow_gsp_rate(jc,jb)
           rain_con_rate(jc,isubs)    = prm_diag%rain_con_rate(jc,jb)
           snow_con_rate(jc,isubs)    = prm_diag%snow_con_rate(jc,jb)
           graupel_gsp_rate(jc,isubs) = p_graupel_gsp_rate    (jc,jb)
         END DO
         IF( atm_phy_nwp_config(jg)%l2moment) THEN
           DO ic = 1, i_count
             jc = ext_data%atm%idx_lst_t(ic,jb,isubs)
             ! for 1mom ice is included in snow_gsp, but for 2mom not
             snow_gsp_rate(jc,isubs)    = snow_gsp_rate(jc,isubs) + prm_diag%ice_gsp_rate(jc,jb)
             ! here we ignore the different densities of graupel and hail in TERRA (at least for now)
             graupel_gsp_rate(jc,isubs) = graupel_gsp_rate(jc,isubs) + prm_diag%hail_gsp_rate(jc,jb) 
           END DO
         ENDIF
       END DO
       !$acc end parallel


       IF (lsnowtile .AND. itype_snowevap == 3) THEN
         !$acc parallel if(lzacc)
         !$acc loop gang vector 
         DO jc = i_startidx, i_endidx
           IF (lnd_diag%h_snow(jc,jb) > 5.e-4_wp) THEN ! traces of snow are ignored
             lnd_diag%hsnow_max(jc,jb) = MAX(lnd_diag%hsnow_max(jc,jb), lnd_diag%h_snow(jc,jb))
             lnd_diag%snow_age(jc,jb)  = MIN(365._wp,lnd_diag%snow_age(jc,jb) + tcall_sfc_jg/86400._wp)
             ! Tuning factor for reduced snow evaporation (stronger reduction for long-lasting snow cover and during melting phase)
             IF (lnd_diag%snow_age(jc,jb) >= 30._wp) THEN
               sntunefac(jc) = 1._wp + MIN(2._wp,MAX(0._wp,(MIN(120._wp,lnd_diag%snow_age(jc,jb))-30._wp)/45._wp* &
                 (0.5_wp+1.5_wp*(lnd_diag%hsnow_max(jc,jb)-lnd_diag%h_snow(jc,jb))/lnd_diag%hsnow_max(jc,jb)) ))
             ELSE
               sntunefac(jc) = 0.75_wp + MAX(0._wp,0.25_wp*(lnd_diag%snow_age(jc,jb)-10._wp)/20._wp)
             ENDIF
           ELSE
             lnd_diag%hsnow_max(jc,jb) = 0._wp
             lnd_diag%snow_age(jc,jb)  = 0._wp
             sntunefac(jc) = 1._wp
           ENDIF
         ENDDO
         !$acc end parallel

         !$acc parallel if(lzacc)
         !$acc loop gang vector private(i_count)
         DO isubs = ntiles_lnd+1, ntiles_total
           i_count = ext_data%atm%gp_count_t(jb,isubs) 
!$NEC ivdep
           !$acc loop private(jc,tmp1)
           DO ic = 1, i_count
             jc = ext_data%atm%idx_lst_t(ic,jb,isubs)
             ! Another tuning factor in order to treat partial snow cover different for fresh snow and 'old' snow
             IF (sntunefac(jc) < 1._wp) THEN
               sntunefac2(jc,isubs) = 4._wp*(1._wp-sntunefac(jc))*lnd_diag%snowfrac_lc_t(jc,jb,isubs) + &
                                      4._wp*(sntunefac(jc)-0.75_wp)
             ELSE
               sntunefac2(jc,isubs) = 1._wp
             ENDIF
             !
             ! parameterization for snow drift, treated as a source term for cloud ice (restricted to glaciers in order
             ! to avoid erroneous side effects on snow density)
             ! Note that, consistent with the approximations made for evaporation and deposition of precipitation,
             ! the related change of total air mass is neglected here
             !
             IF (ext_data%atm%lc_class_t(jc,jb,isubs) == ext_data%atm%i_lc_snow_ice) THEN
               tmp1 = tcall_sfc_jg * 7.5e-10_wp * (600._wp-lnd_prog_now%rho_snow_t(jc,jb,isubs))* &
                 MAX(0._wp,prm_diag%gust10(jc,jb)-7.5_wp)**2

               p_prog_rcf%tracer(jc,nlev,jb,iqi) = p_prog_rcf%tracer(jc,nlev,jb,iqi) + tmp1 * &
                 ext_data%atm%frac_t(jc,jb,isubs) / (p_prog%rho(jc,nlev,jb)*p_metrics%ddqz_z_full(jc,nlev,jb))

               lnd_prog_now%w_snow_t(jc,jb,isubs) = lnd_prog_now%w_snow_t(jc,jb,isubs) - tmp1/rhoh2o
               lnd_diag%h_snow_t(jc,jb,isubs) = rhoh2o*lnd_prog_now%w_snow_t(jc,jb,isubs)/lnd_prog_now%rho_snow_t(jc,jb,isubs)
             ENDIF

           ENDDO
         ENDDO
         !$acc end parallel
       ELSE IF (lsnowtile) THEN
         !$acc parallel if(lzacc)
         !$acc loop gang vector
         DO jc = i_startidx, i_endidx
           sntunefac(jc) = 1._wp
         ENDDO
         !$acc end parallel
         !$acc parallel if (lzacc)
         !$acc loop gang vector private(i_count)
         DO isubs = ntiles_lnd+1, ntiles_total
           i_count = ext_data%atm%gp_count_t(jb,isubs) 
           !$acc loop private(jc)
           DO ic = 1, i_count
             jc = ext_data%atm%idx_lst_t(ic,jb,isubs)
             sntunefac2(jc,isubs) = lnd_diag%snowfrac_lc_t(jc,jb,isubs)
           ENDDO
         ENDDO
         !$acc end parallel
       ENDIF

!---------- Preparations for TERRA in the case if snow tiles are considered
       IF(lsnowtile) THEN      ! snow is considered as separate tiles
         !$acc parallel if(lzacc)
         !$acc loop gang vector private(isubs_snow,i_count_snow)
         DO isubs = 1, ntiles_lnd

           isubs_snow = isubs + ntiles_lnd
           i_count_snow = ext_data%atm%gp_count_t(jb,isubs_snow) 

!$NEC ivdep
           !$acc loop private(jc, tmp3, tmp1, tmp2)
           DO ic = 1, i_count_snow
             jc = ext_data%atm%idx_lst_t(ic,jb,isubs_snow)
  
             ! Snow and rain fall onto snow-covered tile surface only, 
             ! if 
             ! 1) the corresponding snow tile already exists and 
             ! 2) the temperature of snow-free tile is below freezing point (with transition zone between 0 and 1 deg C).
             ! If the temperature of snow-free tile is above freezing point,
             ! precipitation over it will be processed by this tile itself (no snow is created).
             ! If there is no snow tile so far at all, precipitation falls on the snow-free tile,
             ! and the snow tile will be created after TERRA.
             !
             IF (lnd_diag%snowfrac_lc_t(jc,jb,isubs) < 1._wp .AND. lnd_prog_now%t_snow_t(jc,jb,isubs) < tmelt+1._wp) THEN

               ! transition factor to avoid discontinuity at freezing point of soil in snow-free tile
               tmp3 = tmelt + 1._wp - MAX(tmelt,lnd_prog_now%t_snow_t(jc,jb,isubs))
               ! enhancement factor in snow tile
               tmp1 = MAX(1._wp,tmp3/MAX(lnd_diag%snowfrac_lc_t(jc,jb,isubs),0.01_wp))
               ! factor for snow-free tile to ensure that no water gets lost
               tmp2 = (1._wp-tmp1*lnd_diag%snowfrac_lc_t(jc,jb,isubs))/(1._wp-lnd_diag%snowfrac_lc_t(jc,jb,isubs))

               rain_gsp_rate(jc,isubs)    = rain_gsp_rate(jc,isubs)*tmp2
               snow_gsp_rate(jc,isubs)    = snow_gsp_rate(jc,isubs)*tmp2
               rain_con_rate(jc,isubs)    = rain_con_rate(jc,isubs)*tmp2
               snow_con_rate(jc,isubs)    = snow_con_rate(jc,isubs)*tmp2
               graupel_gsp_rate(jc,isubs) = graupel_gsp_rate(jc,isubs)*tmp2
               rain_gsp_rate(jc,isubs_snow)    = rain_gsp_rate(jc,isubs_snow)*tmp1
               snow_gsp_rate(jc,isubs_snow)    = snow_gsp_rate(jc,isubs_snow)*tmp1
               rain_con_rate(jc,isubs_snow)    = rain_con_rate(jc,isubs_snow)*tmp1
               snow_con_rate(jc,isubs_snow)    = snow_con_rate(jc,isubs_snow)*tmp1
               graupel_gsp_rate(jc,isubs_snow) = graupel_gsp_rate(jc,isubs_snow)*tmp1
             END IF
           END DO
         END DO
         !$acc end parallel
       END IF

!---------- Copy input fields for each tile

       !$acc enter data create (soiltyp_t, plcov_t, rootdp_t, sai_t, eai_t, tai_t, laifac_t,      &
       !$acc                    skinc_t, rsmin2d_t, u_t, v_t, t_t, qv_t, p0_t, ps_t, h_snow_gp_t, &
       !$acc                    u_10m_t, v_10m_t, prr_con_t, prs_con_t, conv_frac, prr_gsp_t,     &
       !$acc                    prs_gsp_t, prg_gsp_t, sobs_t, thbs_t, pabs_t, tsnred,             &
       !$acc                    t_snow_now_t, t_s_now_t, t_sk_now_t, t_g_t, qv_s_t, w_snow_now_t, &
       !$acc                    rho_snow_now_t, h_snow_t, w_i_now_t, w_p_now_t, w_s_now_t,        &
       !$acc                    freshsnow_t, snowfrac_t, tch_t, tcm_t, tfv_t, runoff_s_inst_t,    &
       !$acc                    runoff_g_inst_t, t_snow_mult_now_t, rho_snow_mult_now_t,          &
       !$acc                    wliq_snow_now_t, wtot_snow_now_t, dzh_snow_now_t, t_so_now_t,     &
       !$acc                    w_so_now_t, w_so_ice_now_t, t_snow_new_t, t_s_new_t, t_sk_new_t,  &
       !$acc                    w_snow_new_t, rho_snow_new_t, meltrate, w_i_new_t, w_p_new_t,     &
       !$acc                    w_s_new_t, shfl_soil_t, lhfl_soil_t, shfl_snow_t, lhfl_snow_t,    &
       !$acc                    rstom_t, lhfl_bs_t, t_snow_mult_new_t, rho_snow_mult_new_t,       &
       !$acc                    wliq_snow_new_t, wtot_snow_new_t, dzh_snow_new_t, t_so_new_t,     &
       !$acc                    w_so_new_t, w_so_ice_new_t, lhfl_pl_t, shfl_s_t, lhfl_s_t,        &
       !$acc                    qhfl_s_t, plevap_t, z0_t, sso_sigma_t,                            &
       !$acc                    snowfrac_lcu_t, lc_class_t),                                      &
       !$acc if(lzacc)

!----------------------------------
       DO isubs = 1,ntiles_total
!----------------------------------

        i_count = ext_data%atm%gp_count_t(jb,isubs) 

        IF (i_count == 0) CYCLE ! skip loop if the index list for the given tile is empty


!$NEC ivdep
        !$acc parallel if(lzacc)
        !$acc loop gang vector private(jc)
        DO ic = 1, i_count
          jc = ext_data%atm%idx_lst_t(ic,jb,isubs)

          ps_t(ic)      =  p_diag%pres_sfc(jc,jb)    
          prr_con_t(ic) =  rain_con_rate(jc,isubs)
          prs_con_t(ic) =  snow_con_rate(jc,isubs)
          conv_frac(ic) =  phy_params(jg)%rcucov*     (1._wp - prm_diag%tropics_mask(jc,jb)) + &
                           phy_params(jg)%rcucov_trop*         prm_diag%tropics_mask(jc,jb)
          prr_gsp_t(ic) =  rain_gsp_rate(jc,isubs)
          prs_gsp_t(ic) =  snow_gsp_rate(jc,isubs)
          prg_gsp_t(ic) =  graupel_gsp_rate(jc,isubs)

          u_t(ic)       =  p_diag%u         (jc,nlev,jb)
          v_t(ic)       =  p_diag%v         (jc,nlev,jb)
          t_t(ic)       =  p_diag%temp      (jc,nlev,jb)     
          qv_t(ic)      =  p_prog_rcf%tracer(jc,nlev,jb,iqv) 
          p0_t(ic)      =  p_diag%pres      (jc,nlev,jb) 
    
          sso_sigma_t(ic)       = ext_data%atm%sso_stdh(jc,jb)
          lc_class_t(ic)        = ext_data%atm%lc_class_t(jc,jb,isubs)

          t_snow_now_t(ic)          =  lnd_prog_now%t_snow_t(jc,jb,isubs) 
          t_s_now_t(ic)             =  lnd_prog_now%t_s_t(jc,jb,isubs)   
          t_sk_now_t(ic)            =  lnd_prog_now%t_sk_t(jc,jb,isubs)
          t_g_t (ic)                =  lnd_prog_now%t_g_t(jc,jb,isubs)
          qv_s_t(ic)                =  lnd_diag%qv_s_t(jc,jb,isubs)
          w_snow_now_t(ic)          =  lnd_prog_now%w_snow_t(jc,jb,isubs)
          rho_snow_now_t(ic)        =  lnd_prog_now%rho_snow_t(jc,jb,isubs)
          w_i_now_t(ic)             =  lnd_prog_now%w_i_t(jc,jb,isubs)
          h_snow_t(ic)              =  lnd_diag%h_snow_t(jc,jb,isubs)
          freshsnow_t(ic)           =  lnd_diag%freshsnow_t(jc,jb,isubs)
          snowfrac_t(ic)            =  lnd_diag%snowfrac_t(jc,jb,isubs)
          snowfrac_lcu_t(ic)        =  lnd_diag%snowfrac_lcu_t(jc,jb,isubs)

          IF (isubs > ntiles_lnd) THEN ! snowtiles
            ! grid-point averaged snow depth needed for snow aging parameterization
            h_snow_gp_t(ic)         =  MAX(lnd_diag%snowfrac_lc_t(jc,jb,isubs),0.01_wp)*h_snow_t(ic)
          ELSE
            h_snow_gp_t(ic)         =  h_snow_t(ic)
          ENDIF

          IF (itype_interception == 2) THEN
            w_p_now_t(ic)             =  lnd_prog_now%w_p_t(jc,jb,isubs)
            w_s_now_t(ic)             =  lnd_prog_now%w_s_t(jc,jb,isubs)
          ELSE
            w_p_now_t(ic)             =  0._wp
            w_s_now_t(ic)             =  0._wp
          END IF

          IF (itype_trvg == 3) THEN
            plevap_t(ic)            =  lnd_diag%plantevap_t(jc,jb,isubs)
          ELSE
            plevap_t(ic)            =  0._wp
          ENDIF

          IF (itype_vegetation_cycle >= 2) THEN
            laifac_t(ic)            =  ext_data%atm%laifac_t(jc,jb,isubs)
          ELSE
            laifac_t(ic)            =  1._wp
          ENDIF

          IF (isubs > ntiles_lnd) THEN
            z0_t(ic)                =  prm_diag%gz0_t(jc,jb,isubs-ntiles_lnd)/grav
          ELSE
            z0_t(ic)                =  prm_diag%gz0_t(jc,jb,isubs)/grav
          ENDIF

          ! note: we reset "runoff_s_inst_t", "runoff_g_inst_t" in
          ! order to obtain the instantaneous values (and not the sum
          ! over forecast) from terra:
          runoff_s_inst_t(ic)       =  0._wp 
          runoff_g_inst_t(ic)       =  0._wp

          u_10m_t(ic)               =  prm_diag%u_10m_t(jc,jb,isubs)
          v_10m_t(ic)               =  prm_diag%v_10m_t(jc,jb,isubs)  
          tch_t(ic)                 =  prm_diag%tch_t(jc,jb,isubs)
          tcm_t(ic)                 =  prm_diag%tcm_t(jc,jb,isubs)
          tfv_t(ic)                 =  prm_diag%tfv_t(jc,jb,isubs)
          sobs_t(ic)                =  prm_diag%swflxsfc_t(jc,jb,isubs) 
          thbs_t(ic)                =  prm_diag%lwflxsfc_t(jc,jb,isubs) 
          pabs_t(ic)                =  prm_diag%swflx_par_sfc(jc,jb) 

          soiltyp_t(ic)             =  ext_data%atm%soiltyp_t(jc,jb,isubs)
          plcov_t(ic)               =  ext_data%atm%plcov_t(jc,jb,isubs)
          rootdp_t(ic)              =  ext_data%atm%rootdp_t(jc,jb,isubs)
          sai_t(ic)                 =  ext_data%atm%sai_t(jc,jb,isubs)
          tai_t(ic)                 =  ext_data%atm%tai_t(jc,jb,isubs)
          eai_t(ic)                 =  ext_data%atm%eai_t(jc,jb,isubs)
          skinc_t(ic)               =  ext_data%atm%skinc_t(jc,jb,isubs)
          rsmin2d_t(ic)             =  ext_data%atm%rsmin2d_t(jc,jb,isubs)

          t_so_now_t(ic,nlev_soil+1)= lnd_prog_now%t_so_t(jc,nlev_soil+1,jb,isubs)

          IF(lmulti_snow) THEN
            t_snow_mult_now_t(ic,nlev_snow+1) = lnd_prog_now%t_snow_mult_t(jc,nlev_snow+1,jb,isubs)
          ENDIF

          IF (l2lay_rho_snow) THEN ! only level 1 is actually needed as input
            rho_snow_mult_now_t(ic,1) = lnd_prog_now%rho_snow_mult_t(jc,1,jb,isubs)
          ENDIF

        ENDDO
        !$acc end parallel

        IF (itype_snowevap == 1 .OR. .NOT. lsnowtile) THEN
          !$acc parallel if(lzacc)
          !$acc loop gang vector
          DO ic = 1, i_count
            tsnred(ic) = 0._wp
          ENDDO
          !$acc end parallel
        ELSE IF (isubs > ntiles_lnd) THEN
          ! compute temperature offset for reducing snow evaporation in vegetated areas,
          ! parameterizing the temperature difference between the snow and the snow-vegetation-mixture
          ! represented by the variable t_snow and the related snow albedo

          !$acc parallel if(lzacc)
          !$acc loop gang vector private(jc,tmp1,qsat1,dqsdt1,tmp2,qsat2,dqsdt2,tmp2)
          DO ic = 1, i_count
            jc = ext_data%atm%idx_lst_t(ic,jb,isubs)
            tmp1 = 0.06_wp * sntunefac(jc) * sobs_t(ic) * (1._wp - MIN(1._wp,                     &
              (1._wp - (csalb_snow_min + freshsnow_t(ic)*(csalb_snow_max-csalb_snow_min))) /      &
              (1._wp - prm_diag%albdif_t(jc,jb,isubs)) )) * sntunefac2(jc,isubs)
            qsat1 = spec_humi(sat_pres_ice(t_snow_now_t(ic)),ps_t(ic) )
            dqsdt1 = dqsatdT_ice(qsat1,t_snow_now_t(ic))
            tmp2 = tmp1 * (0.1_wp + 1000._wp*MAX(0._wp,qsat1-qv_t(ic))) / (1._wp + tmp1*1000._wp*dqsdt1)
            qsat2 = spec_humi(sat_pres_ice(t_snow_now_t(ic)-tmp2),ps_t(ic) )
            dqsdt2 = dqsatdT_ice(qsat2,t_snow_now_t(ic)-tmp2)
            tmp2 = tmp1 * (0.1_wp + 1000._wp*MAX(0._wp,qsat1-qv_t(ic))) / (1._wp + tmp1*500._wp*(dqsdt1+dqsdt2))
            tsnred(ic) = MIN(10._wp*SQRT(sntunefac(jc)),tmp2) / (1._wp + 4.e-5_wp*z0_t(ic)*sso_sigma_t(ic)**2)
          ENDDO
          !$acc end parallel
        ELSE
          ! If the snow-cover fraction is artificially reduced by the melting-rate parameterization, the bare soil evaporation
          ! in TERRA is turned off on the corresponding snow-free tile.
          ! This is controlled by negative values of tsnred

          !$acc parallel if(lzacc)
          !$acc loop gang vector private(jc)
          DO ic = 1, i_count
            jc = ext_data%atm%idx_lst_t(ic,jb,isubs)
            tsnred(ic) = MIN(2._wp,sntunefac(jc))*(lnd_diag%snowfrac_lc_t(jc,jb,isubs)-snowfrac_lcu_t(ic))/ &
                         (1._wp-lnd_diag%snowfrac_lc_t(jc,jb,isubs))
          ENDDO
          !$acc end parallel
        ENDIF

       MSNOWI: IF(lmulti_snow) THEN
        
       !$acc parallel if(lzacc)
#ifdef __LOOP_EXCHANGE
        !$acc loop gang vector private (jc)
        DO ic = 1, i_count
          jc = ext_data%atm%idx_lst_t(ic,jb,isubs)
          !$acc loop
          DO jk=1,nlev_snow
#else
        !$acc loop gang vector collapse(2) private(jc)  
        DO jk=1,nlev_snow
          DO ic = 1, i_count
            jc = ext_data%atm%idx_lst_t(ic,jb,isubs)
#endif
            t_snow_mult_now_t  (ic,jk) = lnd_prog_now%t_snow_mult_t  (jc,jk,jb,isubs) 
            rho_snow_mult_now_t(ic,jk) = lnd_prog_now%rho_snow_mult_t(jc,jk,jb,isubs)
            wliq_snow_now_t    (ic,jk) = lnd_prog_now%wliq_snow_t    (jc,jk,jb,isubs) 
            wtot_snow_now_t    (ic,jk) = lnd_prog_now%wtot_snow_t    (jc,jk,jb,isubs)
            dzh_snow_now_t     (ic,jk) = lnd_prog_now%dzh_snow_t     (jc,jk,jb,isubs) 
          ENDDO
        ENDDO
        !$acc end parallel
       END IF MSNOWI

        !$acc parallel if(lzacc)
#ifdef __LOOP_EXCHANGE
        !$acc loop gang vector private(jc)
        DO ic = 1, i_count   
          jc = ext_data%atm%idx_lst_t(ic,jb,isubs)
          !$acc loop 
          DO jk=1,nlev_soil
#else
!$NEC outerloop_unroll(nlsoil)
        !$acc loop gang vector collapse(2) private(jc)
        DO jk=1,nlev_soil
          DO ic = 1, i_count
            jc = ext_data%atm%idx_lst_t(ic,jb,isubs)
#endif
            t_so_now_t    (ic,jk) = lnd_prog_now%t_so_t    (jc,jk,jb,isubs) 
            w_so_now_t    (ic,jk) = lnd_prog_now%w_so_t    (jc,jk,jb,isubs) 
            w_so_ice_now_t(ic,jk) = lnd_prog_now%w_so_ice_t(jc,jk,jb,isubs)
          ENDDO
        ENDDO
       !$acc end parallel

!---------- END Copy index list fields
#ifndef COUP_OAS_ICON
        CALL terra (                                           &
        &  nvec         = nproma                             , & !IN array dimensions
        &  ivstart      = 1                                  , & !IN optional start/end indicies
        &  ivend        = i_count                            , & !IN optional start/end indicies
        &  iblock       = jb                                 , & !IN actual block number
        &  ke_soil      = nlev_soil-1                        , & !IN without lowermost (climat.) soil layer
        &  ke_snow      = nlev_snow                          , & !IN without lowermost (climat.) soil layer
        &  ke_soil_hy   = ibot_w_so                          , & !IN number of hydrological active soil layers
        &  zmls         = zml_soil                           , & !IN processing soil level structure 
        &  icant        = icant                              , & !IN canopy-type
        &  nclass_gscp  = atm_phy_nwp_config(jg)%nclass_gscp , & !IN number of hydrometeor classes
        &  dt           = tcall_sfc_jg                       , & !IN 
        &  soiltyp_subs = soiltyp_t                          , & !IN type of the soil (keys 0-9)         --    
        &  plcov        = plcov_t                            , & !IN fraction of plant cover             --
        &  rootdp       = rootdp_t                           , & !IN depth of the roots                ( m  )
        &  sai          = sai_t                              , & !IN surface area index                  --
        &  tai          = tai_t                              , & !IN surface area index                  --
        &  laifac       = laifac_t                           , & !IN ratio between current LAI and laimax                 --
        &  eai          = eai_t                              , & !IN surface area index                  --
        &  skinc        = skinc_t                            , & !IN skin conductivity                 ( W/m**2/K )
        &  rsmin2d      = rsmin2d_t                          , & !IN minimum stomata resistance        ( s/m )
        &  z0           = z0_t                               , & !IN vegetation roughness length        ( m )
!
        &  u            =  u_t                               , & !IN zonal wind speed
        &  v            =  v_t                               , & !IN meridional wind speed 
        &  t            =  t_t                               , & !IN temperature                       (  K  )
        &  qv           =  qv_t                              , & !IN specific water vapor content      (kg/kg)
        &  ptot         =  p0_t                              , & !IN base state pressure               ( Pa  ) 
        &  ps           =  ps_t                              , & !IN surface pressure                  ( Pa  )
!
        &  t_snow_now    = t_snow_now_t                      , & !INOUT temperature of the snow-surface (  K  )
        &  t_snow_new    = t_snow_new_t                      , & !OUT temperature of the snow-surface   (  K  )
!
        &  t_snow_mult_now = t_snow_mult_now_t               , & !INOUT temperature of the snow-surface (  K  )
        &  t_snow_mult_new = t_snow_mult_new_t               , & !OUT temperature of the snow-surface   (  K  )
!
        &  t_s_now       = t_s_now_t                         , & !INOUT temperature of the ground surface (  K  )
        &  t_s_new       = t_s_new_t                         , & !OUT temperature of the ground surface   (  K  )
!
        &  t_sk_now      = t_sk_now_t                        , & !INOUT skin temperature                  (  K  )
        &  t_sk_new      = t_sk_new_t                        , & !OUT skin temperature                    (  K  )
!
        &  t_g           = t_g_t                             , & !INOUT weighted surface temperature      (  K  )
        &  qv_s          = qv_s_t                            , & !INOUT specific humidity at the surface  (kg/kg)
!
        &  w_snow_now    = w_snow_now_t                      , & !INOUT water content of snow         (m H2O) 
        &  w_snow_new    = w_snow_new_t                      , & !OUT water content of snow           (m H2O) 
!
        &  rho_snow_now      = rho_snow_now_t                , & !IN  snow density                    (kg/m**3)
        &  rho_snow_new      = rho_snow_new_t                , & !OUT snow density                    (kg/m**3)
!
        &  rho_snow_mult_now = rho_snow_mult_now_t           , & !INOUT snow density               (kg/m**3) 
        &  rho_snow_mult_new = rho_snow_mult_new_t           , & !OUT snow density                 (kg/m**3) 
!
        &  h_snow        = h_snow_t                          , & !INOUT snow height
        &  h_snow_gp     = h_snow_gp_t                       , & !IN grid-point averaged snow height
        &  meltrate      = meltrate                          , & !OUT snow melting rate
        &  tsnred        = tsnred                            , & !IN temperature offset for computing snow evaporation
!
        &  w_i_now       = w_i_now_t                         , & !INOUT water content of interception water(m H2O)
        &  w_i_new       = w_i_new_t                         , & !OUT water content of interception water(m H2O)
!
        &  w_p_now       = w_p_now_t                         , & !INOUT water content of interception water(m H2O)
        &  w_p_new       = w_p_new_t                         , & !OUT water content of interception water(m H2O)
!
        &  w_s_now       = w_s_now_t                         , & !INOUT water content of interception water(m H2O)
        &  w_s_new       = w_s_new_t                         , & !OUT water content of interception water(m H2O)
!
        &  t_so_now      = t_so_now_t                        , & !INOUT soil temperature (main level)    (  K  )
        &  t_so_new      = t_so_new_t                        , & !OUT soil temperature (main level)      (  K  )
!
        &  w_so_now      = w_so_now_t                        , & !IN  total water content (ice + liquid water) (m H20)
        &  w_so_new      = w_so_new_t                        , & !OUT total water content (ice + liquid water) (m H20)
!
        &  w_so_ice_now  = w_so_ice_now_t                    , & !IN  ice content   (m H20)
        &  w_so_ice_new  = w_so_ice_new_t                    , & !OUT ice content   (m H20)
!
        &  u_10m         = u_10m_t                           , & !IN zonal wind in 10m                 ( m/s )
        &  v_10m         = v_10m_t                           , & !IN meridional wind in 10m            ( m/s )
        &  freshsnow     = freshsnow_t                       , & !INOUT indicator for age of snow in top of snow layer (  -  )
        &  zf_snow       = snowfrac_t                        , & !INOUT snow-cover fraction                            (  -  )
!
        &  wliq_snow_now = wliq_snow_now_t                   , & !INOUT liquid water content in the snow     (m H2O)
        &  wliq_snow_new = wliq_snow_new_t                   , & !OUT liquid water content in the snow       (m H2O)
!                                                            
        &  wtot_snow_now = wtot_snow_now_t                   , & !INOUT total (liquid + solid) water content of snow  (m H2O)
        &  wtot_snow_new = wtot_snow_new_t                   , & !OUT total (liquid + solid) water content of snow  (m H2O)
!
        &  dzh_snow_now  = dzh_snow_now_t                    , & !INOUT layer thickness between half levels in snow (  m  )
        &  dzh_snow_new  = dzh_snow_new_t                    , & !OUT layer thickness between half levels in snow   (  m  )
!
        &  prr_con       = prr_con_t                         , & !IN precipitation rate of rain, convective       (kg/m2*s)
        &  prs_con       = prs_con_t                         , & !IN precipitation rate of snow, convective       (kg/m2*s)
        &  conv_frac     = conv_frac                         , & !IN convective area fraction
        &  prr_gsp       = prr_gsp_t                         , & !IN precipitation rate of rain, grid-scale       (kg/m2*s)
        &  prs_gsp       = prs_gsp_t                         , & !IN precipitation rate of snow, grid-scale       (kg/m2*s)
        &  prg_gsp       = prg_gsp_t                         , & !IN precipitation rate of graupel, grid-scale    (kg/m2*s)
!
        &  tch           = tch_t                             , & !INOUT turbulent transfer coefficient for heat     ( -- )
        &  tcm           = tcm_t                             , & !INOUT turbulent transfer coefficient for momentum ( -- )
        &  tfv           = tfv_t                             , & !INOUT laminar reduction factor for evaporation    ( -- )
!
        &  sobs          = sobs_t                            , & !IN solar radiation at the ground               (W/m2)
        &  thbs          = thbs_t                            , & !IN thermal radiation at the ground             (W/m2)
        &  pabs          = pabs_t                            , & !IN photosynthetic active radiation             (W/m2)
!
        &  runoff_s      = runoff_s_inst_t                   , & !INOUT surface water runoff   (kg/m2)
        &  runoff_g      = runoff_g_inst_t                   , & !INOUT soil water runoff      (kg/m2)
!
        &  zshfl_s       = shfl_soil_t                       , & !OUT sensible heat flux soil/air interface    (W/m2) 
        &  zlhfl_s       = lhfl_soil_t                       , & !OUT latent   heat flux soil/air interface    (W/m2) 
        &  zshfl_snow    = shfl_snow_t                       , & !OUT sensible heat flux snow/air interface    (W/m2) 
        &  zlhfl_snow    = lhfl_snow_t                       , & !OUT latent   heat flux snow/air interface    (W/m2) 
        &  lhfl_bs       = lhfl_bs_t                         , & !OUT latent heat flux from bare soil evap.    (W/m2)
        &  lhfl_pl       = lhfl_pl_t                         , & !OUT latent heat flux from bare soil evap.    (W/m2)
        &  plevap        = plevap_t                          , & !INOUT accumulated plant evaporation          (kg/m2)
        &  rstom         = rstom_t                           , & !OUT stomatal resistance                      ( s/m )
        &  zshfl_sfc     = shfl_s_t                          , & !OUT sensible heat flux surface interface     (W/m2) 
        &  zlhfl_sfc     = lhfl_s_t                          , & !OUT latent   heat flux surface interface     (W/m2) 
        &  zqhfl_sfc     = qhfl_s_t                          , & !OUT moisture flux surface interface          (kg/m2/s)
        &  lacc          = lzacc                             )   !IN flag to run OpenACC code
#else
        prm_diag%umfl_s_t(:,:,1) = -oas_rcv_field_icon(:,:,4)
        prm_diag%vmfl_s_t(:,:,1) = -oas_rcv_field_icon(:,:,5)
        prm_diag%shfl_s_t(:,:,1) = -oas_rcv_field_icon(:,:,6) 
        prm_diag%lhfl_s_t(:,:,1) = -oas_rcv_field_icon(:,:,7)

        lnd_prog_new%t_s_t     (:,:,1) = oas_rcv_field_icon(:,:,8)
        lnd_prog_new%t_g_t     (:,:,1) = oas_rcv_field_icon(:,:,8)

        ! SBr, CHa: follow steps as in terra_multlay
         DO jc = i_startidx, i_endidx
           !qhfl_help(jc,jb) = 0.5_wp+SIGN(0.5_wp, p_lnd_state(1)%p_prog_lnd(nnow(1))%t_s_t(jc,jb,1))
           qhfl_help(jc,jb) = 0.5_wp+SIGN(0.5_wp, oas_rcv_field_icon(jc,jb,8) - t0_melt)
           prm_diag%qhfl_s_t(jc,jb,1) = prm_diag%lhfl_s_t(jc,jb,1) / &
                 (qhfl_help(jc,jb)*lh_v + (1._wp-qhfl_help(jc,jb))*lh_s) !SPo correct calculation
         END DO

        ! SPo calculate effective qv_s based on clm fluxes following
        ! terra_multlay
        DO jc = i_startidx, i_endidx
           qv_s_help(jc,jb) =  prm_diag%tvh_t(jc,jb,1) * &
               p_diag%pres_sfc(jc,jb)/(r_d*lnd_prog_new%t_g_t(jc,jb,1) * &
               (1._wp + rvd_m_o * lnd_diag%qv_s_t(jc,jb,1))) + 1.E-6_wp
           lnd_diag%qv_s_t(jc,jb,1) = p_prog_rcf%tracer(jc,nlev,jb,iqv) - &
                                      prm_diag%qhfl_s_t(jc,jb,1) / qv_s_help(jc,jb)
           ! SPo qv_s may violate the saturation constraint 
           lnd_diag%qv_s_t(jc,jb,1) = MIN(lnd_diag%qv_s_t(jc,jb,1), &
               spec_humi(sat_pres_water(lnd_prog_new%t_g_t(jc,jb,1)),p_diag%pres_sfc(jc,jb)))
!           !SPo limiter
!           lnd_diag%qv_s_t(jc,jb,1) = min(max(1.E-5_wp,lnd_diag%qv_s_t(jc,jb,1)),1E-1_wp)
        END DO


#endif

        ! Multiply w_snow with old snow fraction in order to obtain the area-average SWE needed for
        ! diagnosing the new snow fraction
        IF (lsnowtile .AND. isubs > ntiles_lnd) THEN
          !$acc parallel if(lzacc)
          !$acc loop gang vector private(jc)
          DO ic = 1, i_count
            jc = ext_data%atm%idx_lst_t(ic,jb,isubs)
            w_snow_now_t(ic) = w_snow_new_t(ic)*MAX(lnd_diag%snowfrac_lc_t(jc,jb,isubs),0.01_wp)
            meltrate(ic) = meltrate(ic)*MAX(lnd_diag%snowfrac_lc_t(jc,jb,isubs),0.01_wp)
          ENDDO
          !$acc end parallel
        ELSE
          !$acc parallel if(lzacc)
          !$acc loop gang vector
          DO ic = 1, i_count
            w_snow_now_t(ic) = w_snow_new_t(ic)
          ENDDO
          !$acc end parallel
        ENDIF


        CALL diag_snowfrac_tg(                      &
          &  istart     = 1, iend = i_count       , & ! start/end indices
          &  lc_class   = lc_class_t              , & ! land-cover class
          &  i_lc_urban = ext_data%atm%i_lc_urban , & ! land-cover class index for urban areas
          &  t_snow     = t_snow_new_t            , & ! snow temperature
          &  t_soiltop  = t_sk_new_t              , & ! soil top temperature or skin temperature
          &  w_snow     = w_snow_now_t            , & ! snow WE
          &  rho_snow   = rho_snow_new_t          , & ! snow density
          &  freshsnow  = freshsnow_t             , & ! fresh snow fraction
          &  meltrate   = meltrate                , & ! snow melting rate
          &  sso_sigma  = sso_sigma_t             , & ! sso stdev
          &  z0         = z0_t                    , & ! vegetation roughness length
          &  snowfrac   = snowfrac_t              , & ! OUT: snow cover fraction
          &  snowfrac_u = snowfrac_lcu_t          , & ! OUT: unmodified snow cover fraction
          &  t_g        = t_g_t                   , & ! OUT: averaged ground temperature
          &  lacc       = lzacc                     ) ! flag for OpenACC


!$NEC ivdep
        !$acc parallel if(lzacc)
        !$acc loop gang vector private(jc,tmp1,tmp2)
        DO ic = 1, i_count
          jc = ext_data%atm%idx_lst_t(ic,jb,isubs)

!---------- Further processing of snow-cover fraction in case of artificial reduction during melting phase

          ! Avoid spreading of melting snow on warm surface before sunset
          IF (isubs > ntiles_lnd .AND. snowfrac_t(ic) > lnd_diag%snowfrac_lc_t(jc,jb,isubs)) THEN 
            IF (meltrate(ic) > 0._wp) THEN
              tmp1 = MAX(0._wp,0.02_wp*(50._wp-prm_diag%swflxsfc_t(jc,jb,isubs-ntiles_lnd)))
              tmp2 = MIN(1._wp,MAX(0._wp,tmelt+1._wp-lnd_prog_new%t_s_t(jc,jb,isubs-ntiles_lnd)))
              snowfrac_t(ic) = MIN(snowfrac_t(ic),lnd_diag%snowfrac_lc_t(jc,jb,isubs)+MAX(tmp1,tmp2)*tcall_sfc_jg/10800._wp)
            ELSE IF (prs_gsp_t(ic) + prs_con_t(ic) + prg_gsp_t(ic) == 0._wp) THEN
              snowfrac_t(ic) = MIN(snowfrac_t(ic),lnd_diag%snowfrac_lc_t(jc,jb,isubs)+tcall_sfc_jg/7200._wp)
            ELSE
              snowfrac_t(ic) = MIN(snowfrac_t(ic),lnd_diag%snowfrac_lc_t(jc,jb,isubs)+tcall_sfc_jg/1800._wp)
            ENDIF
          ENDIF

          ! Remark: snowfrac_t and snowfrac_lc_t differ only if lsnowtile=true (see below)  
          lnd_diag%snowfrac_lc_t (jc,jb,isubs) = snowfrac_t    (ic) 
          lnd_diag%snowfrac_t    (jc,jb,isubs) = snowfrac_t    (ic)
          lnd_diag%snowfrac_lcu_t(jc,jb,isubs) = snowfrac_lcu_t(ic)

!---------- Copy remaining index list fields back to state fields

          lnd_prog_new%t_snow_t  (jc,jb,isubs) = t_snow_new_t  (ic)        
#ifndef COUP_OAS_ICON 
          lnd_prog_new%t_s_t     (jc,jb,isubs) = t_s_new_t     (ic)              
          lnd_prog_new%t_sk_t    (jc,jb,isubs) = t_sk_new_t    (ic)
          lnd_prog_new%t_g_t     (jc,jb,isubs) = t_g_t         (ic)
          ! qv_s may violate the saturation constraint in cases of numerical instability
          lnd_diag%qv_s_t        (jc,jb,isubs) = MIN(qv_s_t    (ic), &
            spec_humi(sat_pres_water(t_g_t(ic)),ps_t(ic)) )
#endif
          lnd_prog_new%w_snow_t  (jc,jb,isubs) = w_snow_new_t  (ic)          
          lnd_prog_new%rho_snow_t(jc,jb,isubs) = rho_snow_new_t(ic)        
          lnd_diag%h_snow_t      (jc,jb,isubs) = h_snow_t      (ic)
          lnd_prog_new%w_i_t     (jc,jb,isubs) = w_i_new_t     (ic)
          IF (itype_interception == 2) THEN
            lnd_prog_new%w_p_t     (jc,jb,isubs) = w_p_new_t     (ic)             
            lnd_prog_new%w_s_t     (jc,jb,isubs) = w_s_new_t     (ic)     
          END IF
          lnd_diag%freshsnow_t   (jc,jb,isubs) = freshsnow_t   (ic) 
          lnd_diag%runoff_s_inst_t    (jc,jb,isubs) = runoff_s_inst_t    (ic)  
          lnd_diag%runoff_g_inst_t    (jc,jb,isubs) = runoff_g_inst_t    (ic)  

          lnd_prog_new%t_so_t(jc,nlev_soil+1,jb,isubs) = t_so_new_t(ic,nlev_soil+1)

          prm_diag%lhfl_bs_t     (jc,jb,isubs) = lhfl_bs_t     (ic)
          lnd_diag%rstom_t       (jc,jb,isubs) = rstom_t       (ic)
#ifndef COUP_OAS_ICON
          prm_diag%shfl_s_t      (jc,jb,isubs) = shfl_s_t      (ic)
          prm_diag%lhfl_s_t      (jc,jb,isubs) = lhfl_s_t      (ic)
          prm_diag%qhfl_s_t      (jc,jb,isubs) = qhfl_s_t      (ic)
#endif

          IF (itype_trvg == 3) lnd_diag%plantevap_t(jc,jb,isubs) = plevap_t(ic)     

          IF(lmulti_snow) THEN
            lnd_prog_new%t_snow_mult_t(jc,nlev_snow+1,jb,isubs) = t_snow_mult_new_t(ic,nlev_snow+1)
          ENDIF

          IF (l2lay_rho_snow) THEN
            lnd_prog_new%rho_snow_mult_t(jc,1:2,jb,isubs) = rho_snow_mult_new_t(ic,1:2)
          ENDIF

        ENDDO
        !$acc end parallel 

        IF (lsnowtile .AND. isubs > ntiles_lnd) THEN ! copy snowfrac_t to snow-free tile
!$NEC ivdep                                          ! (needed for index list computation)
          !$acc parallel if(lzacc)
          !$acc loop gang vector private(jc)
          DO ic = 1, i_count
            jc = ext_data%atm%idx_lst_t(ic,jb,isubs)
            lnd_diag%snowfrac_lc_t(jc,jb,isubs-ntiles_lnd)  = lnd_diag%snowfrac_lc_t(jc,jb,isubs)
            lnd_diag%snowfrac_lcu_t(jc,jb,isubs-ntiles_lnd) = lnd_diag%snowfrac_lcu_t(jc,jb,isubs)
          ENDDO
          !$acc end parallel
        ENDIF


        MSNOWO: IF(lmulti_snow) THEN

        !$acc parallel if(lzacc)
#ifdef __LOOP_EXCHANGE
        !$acc loop gang vector private(jc)
        DO ic = 1, i_count
          jc = ext_data%atm%idx_lst_t(ic,jb,isubs)
          !$acc loop
          DO jk=1,nlev_snow
#else
        !$acc loop gang vector private(jc) collapse(2)   
        DO jk=1,nlev_snow
!$NEC ivdep
          DO ic = 1, i_count
            jc = ext_data%atm%idx_lst_t(ic,jb,isubs)
#endif
            lnd_prog_new%t_snow_mult_t  (jc,jk,jb,isubs) = t_snow_mult_new_t  (ic,jk)   
            lnd_prog_new%rho_snow_mult_t(jc,jk,jb,isubs) = rho_snow_mult_new_t(ic,jk) 
            lnd_prog_new%wliq_snow_t    (jc,jk,jb,isubs) = wliq_snow_new_t    (ic,jk)     
            lnd_prog_new%wtot_snow_t    (jc,jk,jb,isubs) = wtot_snow_new_t    (ic,jk)     
            lnd_prog_new%dzh_snow_t     (jc,jk,jb,isubs) = dzh_snow_new_t     (ic,jk)      
          ENDDO
        ENDDO
        !$acc end parallel
        END IF MSNOWO

        !$acc parallel if(lzacc)
#ifdef __LOOP_EXCHANGE
        !$acc loop gang vector private(jc)
        DO ic = 1, i_count
          jc = ext_data%atm%idx_lst_t(ic,jb,isubs)
          !$acc loop
          DO jk=1,nlev_soil
#else
!$NEC outerloop_unroll(nlsoil)
        !$acc loop gang vector private(jc) collapse(2)
        DO jk=1,nlev_soil
!$NEC ivdep
          DO ic = 1, i_count
            jc = ext_data%atm%idx_lst_t(ic,jb,isubs)
#endif
            lnd_prog_new%t_so_t    (jc,jk,jb,isubs) = t_so_new_t    (ic,jk)          
            lnd_prog_new%w_so_t    (jc,jk,jb,isubs) = w_so_new_t    (ic,jk)          
            lnd_prog_new%w_so_ice_t(jc,jk,jb,isubs) = w_so_ice_new_t(ic,jk)

            ! diagnostic field
            prm_diag%lhfl_pl_t     (jc,jk,jb,isubs) = lhfl_pl_t     (ic,jk)     
          ENDDO
        ENDDO
        !$acc end parallel


       END DO ! isubs - loop over tiles

       !$acc exit data delete (soiltyp_t, plcov_t, rootdp_t, sai_t, eai_t, tai_t, laifac_t,      &
       !$acc                   skinc_t, rsmin2d_t, u_t, v_t, t_t, qv_t, p0_t, ps_t, h_snow_gp_t, &
       !$acc                   u_10m_t, v_10m_t, prr_con_t, prs_con_t, conv_frac, prr_gsp_t,     &
       !$acc                   prs_gsp_t, prg_gsp_t, sobs_t, thbs_t, pabs_t, tsnred,             &
       !$acc                   t_snow_now_t, t_s_now_t, t_sk_now_t, t_g_t, qv_s_t, w_snow_now_t, &
       !$acc                   rho_snow_now_t, h_snow_t, w_i_now_t, w_p_now_t, w_s_now_t,        &
       !$acc                   freshsnow_t, snowfrac_t, tch_t, tcm_t, tfv_t, runoff_s_inst_t,    &
       !$acc                   runoff_g_inst_t, t_snow_mult_now_t, rho_snow_mult_now_t,          &
       !$acc                   wliq_snow_now_t, wtot_snow_now_t, dzh_snow_now_t, t_so_now_t,     &
       !$acc                   w_so_now_t, w_so_ice_now_t, t_snow_new_t, t_s_new_t, t_sk_new_t,  &
       !$acc                   w_snow_new_t, rho_snow_new_t, meltrate, w_i_new_t, w_p_new_t,     &
       !$acc                   w_s_new_t, shfl_soil_t, lhfl_soil_t, shfl_snow_t, lhfl_snow_t,    &
       !$acc                   rstom_t, lhfl_bs_t, t_snow_mult_new_t, rho_snow_mult_new_t,       &
       !$acc                   wliq_snow_new_t, wtot_snow_new_t, dzh_snow_new_t, t_so_new_t,     &
       !$acc                   w_so_new_t, w_so_ice_new_t, lhfl_pl_t, shfl_s_t, lhfl_s_t,        &
       !$acc                   qhfl_s_t, plevap_t, z0_t, sso_sigma_t,                            &
       !$acc                   snowfrac_lcu_t, lc_class_t),                                      &
       !$acc if(lzacc)


       IF(lsnowtile) THEN      ! snow is considered as separate tiles
         DO isubs = 1, ntiles_lnd

           isubs_snow = isubs + ntiles_lnd

           ! save previous area fractions for subsequent redistribution computations
           !$acc kernels if(lzacc)
           frac_sv(:)      = ext_data%atm%frac_t(:,jb,isubs)
           frac_snow_sv(:) = ext_data%atm%frac_t(:,jb,isubs_snow)
           !$acc end kernels

           ! update index lists for snow tiles
           CALL update_idx_lists_lnd (idx_lst_lp       = ext_data%atm%idx_lst_lp_t(:,jb,isubs),         &
                                    lp_count           = ext_data%atm%lp_count_t(jb,isubs),             &
                                    idx_lst            = ext_data%atm%idx_lst_t(:,jb,isubs),            &
                                    gp_count           = ext_data%atm%gp_count_t(jb,isubs),             &
                                    idx_lst_snow       = ext_data%atm%idx_lst_t(:,jb,isubs_snow),       &
                                    gp_count_snow      = ext_data%atm%gp_count_t(jb,isubs_snow),        &
                                    lc_frac            = ext_data%atm%lc_frac_t(:,jb,isubs),            &
                                    partial_frac       = ext_data%atm%frac_t(:,jb,isubs),               &
                                    partial_frac_snow  = ext_data%atm%frac_t(:,jb,isubs_snow),          &
                                    snowtile_flag      = ext_data%atm%snowtile_flag_t(:,jb,isubs),      &
                                    snowtile_flag_snow = ext_data%atm%snowtile_flag_t(:,jb,isubs_snow), &
                                    snowfrac           = lnd_diag%snowfrac_lc_t(:,jb,isubs), &
                                    lacc               = lzacc)
  
           i_count = ext_data%atm%gp_count_t(jb,isubs) 
           i_count_snow = ext_data%atm%gp_count_t(jb,isubs_snow)

           ! Check for newly activated grid points that need to be initialized
           icount_init = 0
!$NEC ivdep
           !$acc parallel if(lzacc)
           !$acc loop seq 
           DO ic = 1, i_count
             jc = ext_data%atm%idx_lst_t(ic,jb,isubs)
             IF (ext_data%atm%snowtile_flag_t(jc,jb,isubs) == 2) THEN
               icount_init = icount_init + 1
               init_list(icount_init) = jc
               it1(icount_init) = isubs      ! target of copy operation
               it2(icount_init) = isubs_snow ! source of copy operation
             ENDIF
           ENDDO
           !$acc end parallel
!$NEC ivdep
           !$acc parallel if(lzacc)
           !$acc loop seq
           DO ic = 1, i_count_snow
             jc = ext_data%atm%idx_lst_t(ic,jb,isubs_snow)
             IF (ext_data%atm%snowtile_flag_t(jc,jb,isubs_snow) == 2) THEN
               icount_init = icount_init + 1
               init_list(icount_init) = jc
               it1(icount_init) = isubs_snow ! target of copy operation
               it2(icount_init) = isubs      ! source of copy operation
             ENDIF
           ENDDO
           !$acc end parallel
!$NEC ivdep
           !$acc parallel if(lzacc)
           !$acc loop gang vector private(jc, is1, is2)
           DO ic = 1, icount_init
             jc = init_list(ic)
             is1 = it1(ic)
             is2 = it2(ic)
             lnd_prog_new%t_snow_t  (jc,jb,is1) = lnd_prog_new%t_snow_t  (jc,jb,is2)        
             lnd_prog_new%t_s_t     (jc,jb,is1) = lnd_prog_new%t_s_t     (jc,jb,is2)       
             lnd_prog_new%t_sk_t    (jc,jb,is1) = lnd_prog_new%t_sk_t    (jc,jb,is2)
             lnd_prog_new%t_g_t     (jc,jb,is1) = lnd_prog_new%t_g_t     (jc,jb,is2) 
             lnd_diag%qv_s_t        (jc,jb,is1) = lnd_diag%qv_s_t        (jc,jb,is2)             
             lnd_prog_new%w_snow_t  (jc,jb,is1) = lnd_prog_new%w_snow_t  (jc,jb,is2)     
             lnd_prog_new%rho_snow_t(jc,jb,is1) = lnd_prog_new%rho_snow_t(jc,jb,is2)
             lnd_diag%h_snow_t      (jc,jb,is1) = lnd_diag%h_snow_t      (jc,jb,is2)
             lnd_prog_new%w_i_t     (jc,jb,is1) = lnd_prog_new%w_i_t     (jc,jb,is2)        

             lnd_diag%freshsnow_t   (jc,jb,is1) = lnd_diag%freshsnow_t   (jc,jb,is2)
             lnd_diag%snowfrac_lc_t (jc,jb,is1) = lnd_diag%snowfrac_lc_t (jc,jb,is2) 
             lnd_diag%snowfrac_lcu_t(jc,jb,is1) = lnd_diag%snowfrac_lcu_t(jc,jb,is2) 
             lnd_diag%snowfrac_t    (jc,jb,is1) = lnd_diag%snowfrac_t    (jc,jb,is2) 
             lnd_diag%runoff_s_inst_t (jc,jb,is1) = lnd_diag%runoff_s_inst_t    (jc,jb,is2)
             lnd_diag%runoff_g_inst_t (jc,jb,is1) = lnd_diag%runoff_g_inst_t    (jc,jb,is2)

             prm_diag%lhfl_bs_t     (jc,jb,is1) = prm_diag%lhfl_bs_t     (jc,jb,is2)
             lnd_diag%rstom_t       (jc,jb,is1) = lnd_diag%rstom_t       (jc,jb,is2)
             prm_diag%shfl_s_t      (jc,jb,is1) = prm_diag%shfl_s_t      (jc,jb,is2)
             prm_diag%lhfl_s_t      (jc,jb,is1) = prm_diag%lhfl_s_t      (jc,jb,is2)
             prm_diag%qhfl_s_t      (jc,jb,is1) = prm_diag%qhfl_s_t      (jc,jb,is2)
             prm_diag%albdif_t      (jc,jb,is1) = prm_diag%albdif_t      (jc,jb,is2)
             !$acc loop
             DO jk= 1, nlev_soil+1
               lnd_prog_new%t_so_t    (jc,jk,jb,is1) = lnd_prog_new%t_so_t    (jc,jk,jb,is2)          
             ENDDO
             !$acc loop
             DO jk = 1, nlev_soil
               lnd_prog_new%w_so_t    (jc,jk,jb,is1) = lnd_prog_new%w_so_t    (jc,jk,jb,is2)        
               lnd_prog_new%w_so_ice_t(jc,jk,jb,is1) = lnd_prog_new%w_so_ice_t(jc,jk,jb,is2)
               prm_diag%lhfl_pl_t     (jc,jk,jb,is1) = prm_diag%lhfl_pl_t     (jc,jk,jb,is2)     
             ENDDO 
             IF (itype_trvg == 3) lnd_diag%plantevap_t(jc,jb,is1) = lnd_diag%plantevap_t(jc,jb,is2)     

             IF (l2lay_rho_snow .OR. lmulti_snow) THEN
               !$acc loop
               DO jk=1,nlev_snow
                 lnd_prog_new%rho_snow_mult_t(jc,jk,jb,is1) = lnd_prog_new%rho_snow_mult_t(jc,jk,jb,is2)
               ENDDO
             ENDIF

             IF (lmulti_snow) THEN
               !$acc loop
               DO jk=1,nlev_snow+1
                 lnd_prog_new%t_snow_mult_t(jc,jk,jb,is1) = lnd_prog_new%t_snow_mult_t  (jc,jk,jb,is2)
               ENDDO
               !$acc loop
               DO jk=1,nlev_snow
                 lnd_prog_new%wliq_snow_t  (jc,jk,jb,is1) = lnd_prog_new%wliq_snow_t    (jc,jk,jb,is2)
                 lnd_prog_new%wtot_snow_t  (jc,jk,jb,is1) = lnd_prog_new%wtot_snow_t    (jc,jk,jb,is2)
                 lnd_prog_new%dzh_snow_t   (jc,jk,jb,is1) = lnd_prog_new%dzh_snow_t     (jc,jk,jb,is2)
               ENDDO
             ENDIF

             IF (itype_interception == 2) THEN
               lnd_prog_new%w_p_t(jc,jb,is1) = lnd_prog_new%w_p_t(jc,jb,is2)        
               lnd_prog_new%w_s_t(jc,jb,is1) = lnd_prog_new%w_s_t(jc,jb,is2)        
             END IF

           ENDDO
           !$acc end parallel
!$NEC ivdep
           !$acc parallel if(lzacc)
           !$acc loop gang vector private(jc)
           DO ic = 1, i_count_snow
             jc = ext_data%atm%idx_lst_t(ic,jb,isubs_snow)

             IF (ext_data%atm%snowtile_flag_t(jc,jb,isubs_snow) == 1 .AND. &
                 ext_data%atm%snowtile_flag_t(jc,jb,isubs)      == 1) THEN

               ! compute factors for redistribution of heat and moisture
               fact1(jc) = MIN(1._wp,frac_sv(jc)/     MAX(small,ext_data%atm%frac_t(jc,jb,isubs)     ))
               fact2(jc) = MIN(1._wp,frac_snow_sv(jc)/MAX(small,ext_data%atm%frac_t(jc,jb,isubs_snow)))
             ENDIF

           END DO
           !$acc end parallel

           ! redistribution of heat and moisture between snow-covered and snow-free tiles 
           ! according to their new fractions, in order to keep heat and moisture balances
           !$acc parallel if(lzacc)
           !$acc loop gang vector private(jc, tmp1, tmp2, tmp3) collapse(2)
           DO jk = 1, nlev_soil
!$NEC ivdep
             DO ic = 1, i_count_snow
               jc = ext_data%atm%idx_lst_t(ic,jb,isubs_snow)

               IF (ext_data%atm%snowtile_flag_t(jc,jb,isubs_snow) == 1 .AND. &
                   ext_data%atm%snowtile_flag_t(jc,jb,isubs)      == 1) THEN

                 tmp1 = lnd_prog_new%t_so_t(jc,jk,jb,isubs) 
                 tmp2 = lnd_prog_new%w_so_t(jc,jk,jb,isubs)
                 tmp3 = lnd_prog_new%w_so_ice_t(jc,jk,jb,isubs)
  
                 lnd_prog_new%t_so_t    (jc,jk,jb,isubs) = lnd_prog_new%t_so_t    (jc,jk,jb,isubs)*fact1(jc) &
                   &                       + lnd_prog_new%t_so_t    (jc,jk,jb,isubs_snow)*(1._wp - fact1(jc))
                 lnd_prog_new%w_so_t    (jc,jk,jb,isubs) = lnd_prog_new%w_so_t    (jc,jk,jb,isubs)*fact1(jc) &
                   &                       + lnd_prog_new%w_so_t    (jc,jk,jb,isubs_snow)*(1._wp - fact1(jc))
                 lnd_prog_new%w_so_ice_t(jc,jk,jb,isubs) = lnd_prog_new%w_so_ice_t(jc,jk,jb,isubs)*fact1(jc) &
                   &                       + lnd_prog_new%w_so_ice_t(jc,jk,jb,isubs_snow)*(1._wp - fact1(jc))
 
                 lnd_prog_new%t_so_t    (jc,jk,jb,isubs_snow) = tmp1*(1._wp - fact2(jc)) &
                   &              + lnd_prog_new%t_so_t    (jc,jk,jb,isubs_snow)*fact2(jc)
                 lnd_prog_new%w_so_t    (jc,jk,jb,isubs_snow) = tmp2*(1._wp - fact2(jc)) &
                   &              + lnd_prog_new%w_so_t    (jc,jk,jb,isubs_snow)*fact2(jc)
                 lnd_prog_new%w_so_ice_t(jc,jk,jb,isubs_snow) = tmp3*(1._wp - fact2(jc)) &
                  &               + lnd_prog_new%w_so_ice_t(jc,jk,jb,isubs_snow)*fact2(jc)

                 IF (jk == 1) THEN
                   lnd_prog_new%t_s_t(jc,jb,isubs)       = lnd_prog_new%t_so_t(jc,jk,jb,isubs)
                   lnd_prog_new%t_s_t(jc,jb,isubs_snow)  = lnd_prog_new%t_so_t(jc,jk,jb,isubs_snow)
                   lnd_prog_new%t_sk_t(jc,jb,isubs)      = lnd_prog_new%t_so_t(jc,jk,jb,isubs)
                   lnd_prog_new%t_sk_t(jc,jb,isubs_snow) = lnd_prog_new%t_so_t(jc,jk,jb,isubs_snow)

                   tmp1 = lnd_prog_new%w_i_t(jc,jb,isubs) 
                   lnd_prog_new%w_i_t(jc,jb,isubs) = tmp1*fact1(jc)           &
                     + lnd_prog_new%w_i_t(jc,jb,isubs_snow)*(1._wp - fact1(jc))
                   lnd_prog_new%w_i_t(jc,jb,isubs_snow) = tmp1*(1._wp - fact2(jc)) &
                     + lnd_prog_new%w_i_t(jc,jb,isubs_snow)*fact2(jc)
                 ENDIF
               ENDIF

             END DO
           END DO        ! soil layers
           !$acc end parallel
!$NEC ivdep

           !$acc parallel if(lzacc)
           !$acc loop gang vector private(jc)
           DO ic = 1, i_count_snow
             jc = ext_data%atm%idx_lst_t(ic,jb,isubs_snow)

             IF (ext_data%atm%snowtile_flag_t(jc,jb,isubs_snow) == 2 .AND. .NOT. lmulti_snow) THEN ! new snow point
               ! in this case, the h_snow and w_snow are not yet rescaled according to the snow-cover fraction
               ! ** In principle, this rescaling also needs to be made for the multi-layer scheme, but this leads         **
               ! ** to a crash because of a division by zero. Adding the rescaling for wliq_snow, wtot_snow and dzh_snow  **
               ! ** (which is still missing here) does NOT cure this problem                                              **
               lnd_prog_new%w_snow_t(jc,jb,isubs_snow) = lnd_prog_new%w_snow_t(jc,jb,isubs_snow) / &
                 MAX(0.01_wp,lnd_diag%snowfrac_t(jc,jb,isubs_snow))
               lnd_diag%h_snow_t(jc,jb,isubs_snow)     = lnd_diag%h_snow_t(jc,jb,isubs_snow) / &
                  MAX(0.01_wp,lnd_diag%snowfrac_t(jc,jb,isubs_snow))
             ELSE
               ! Rescale SWE and snow depth according to changes in the snow cover fraction
               lnd_prog_new%w_snow_t(jc,jb,isubs_snow) = (lnd_prog_new%w_snow_t(jc,jb,isubs_snow) * &
                 frac_snow_sv(jc) + lnd_prog_new%w_snow_t(jc,jb,isubs) * &
                 frac_sv(jc) )/MAX(small,ext_data%atm%frac_t(jc,jb,isubs_snow))
               lnd_diag%h_snow_t(jc,jb,isubs_snow)     = (lnd_diag%h_snow_t(jc,jb,isubs_snow) * &
                 frac_snow_sv(jc) + lnd_diag%h_snow_t(jc,jb,isubs) * &
                 frac_sv(jc) ) /MAX(small,ext_data%atm%frac_t(jc,jb,isubs_snow))
             ENDIF

             ! reset field for actual snow-cover for grid points / land-cover classes for which there
             ! are seperate snow-free and snow-covered tiles 
             lnd_diag%snowfrac_t(jc,jb,isubs)      = 0._wp
             lnd_prog_new%w_snow_t(jc,jb,isubs)    = 0._wp
             lnd_diag%h_snow_t(jc,jb,isubs)        = 0._wp
             lnd_prog_new%t_snow_t(jc,jb,isubs)    = lnd_prog_new%t_s_t(jc,jb,isubs)
             lnd_prog_new%t_g_t(jc,jb,isubs)       = lnd_prog_new%t_sk_t(jc,jb,isubs)

             ! copy rho_snow and freshsnow in order to get the right tile-averaged values
             lnd_prog_new%rho_snow_t(jc,jb,isubs)  = lnd_prog_new%rho_snow_t(jc,jb,isubs_snow)
             lnd_diag%freshsnow_t(jc,jb,isubs)     = lnd_diag%freshsnow_t(jc,jb,isubs_snow)

             ! to prevent numerical stability problems, we require at least 1 cm of snow in order to
             ! have a snow-cover fraction of 1 on snow tiles (not critical for the single-layer
             ! snow scheme, but the multi-layer snow model becomes numerically unstable within a few
             ! time steps when associating traces of snow with a snow-cover fraction of 1)
             lnd_diag%snowfrac_t(jc,jb,isubs_snow) = MIN(1._wp,lnd_diag%h_snow_t(jc,jb,isubs_snow)*100._wp)

             ! Rediagnose t_g according to the modified snow-cover fraction
             lnd_prog_new%t_g_t(jc,jb,isubs_snow) =  &
               lnd_diag%snowfrac_t(jc,jb,isubs_snow) * lnd_prog_new%t_snow_t(jc,jb,isubs_snow) + &
               (1._wp-lnd_diag%snowfrac_t(jc,jb,isubs_snow))*lnd_prog_new%t_sk_t(jc,jb,isubs_snow)

             IF (lmulti_snow) THEN
               lnd_prog_new%t_snow_mult_t(jc,nlev_snow+1,jb,isubs) = lnd_prog_new%t_s_t(jc,jb,isubs)
             ENDIF

             IF (l2lay_rho_snow) THEN
               lnd_prog_new%rho_snow_mult_t(jc,1,jb,isubs) = lnd_prog_new%rho_snow_mult_t(jc,1,jb,isubs_snow)
               lnd_prog_new%rho_snow_mult_t(jc,2,jb,isubs) = lnd_prog_new%rho_snow_mult_t(jc,2,jb,isubs_snow)
             ENDIF

           END DO
           !$acc end parallel

           IF (lmulti_snow) THEN
             !$acc parallel if(lzacc)
             !$acc loop gang vector private(jc) collapse(2)
             DO jk=1,nlev_snow
!$NEC ivdep
               DO ic = 1, i_count_snow
                 jc = ext_data%atm%idx_lst_t(ic,jb,isubs_snow)
                 lnd_prog_new%t_snow_mult_t(jc,jk,jb,isubs) = lnd_prog_new%t_s_t(jc,jb,isubs)
                 lnd_prog_new%wliq_snow_t(jc,jk,jb,isubs) = 0._wp
                 lnd_prog_new%wtot_snow_t(jc,jk,jb,isubs) = 0._wp
                 lnd_prog_new%dzh_snow_t (jc,jk,jb,isubs) = 0._wp

                 ! Rescale mass-related variables according to changes in the snow cover fraction
                 lnd_prog_new%wliq_snow_t(jc,jk,jb,isubs_snow) = lnd_prog_new%wliq_snow_t(jc,jk,jb,isubs_snow) * &
                   frac_snow_sv(jc)/MAX(small,ext_data%atm%frac_t(jc,jb,isubs_snow))
                 lnd_prog_new%wtot_snow_t(jc,jk,jb,isubs_snow) = lnd_prog_new%wtot_snow_t(jc,jk,jb,isubs_snow) * &
                   frac_snow_sv(jc)/MAX(small,ext_data%atm%frac_t(jc,jb,isubs_snow))
                 lnd_prog_new%dzh_snow_t (jc,jk,jb,isubs_snow) = lnd_prog_new%dzh_snow_t (jc,jk,jb,isubs_snow) * &
                   frac_snow_sv(jc)/MAX(small,ext_data%atm%frac_t(jc,jb,isubs_snow))

                 ! copy rho_snow_mult in order to get the right tile-averaged values
                 lnd_prog_new%rho_snow_mult_t(jc,jk,jb,isubs)   = lnd_prog_new%rho_snow_mult_t(jc,jk,jb,isubs_snow)

               ENDDO
             ENDDO
             !$acc end parallel
           ENDIF

         END DO
       ENDIF  !snow tiles

    
      ELSE IF ( atm_phy_nwp_config(jg)%inwp_surface == 2 ) THEN 

          !-------------------------------------------------------------------------
          !> ECHAM version 
          !-------------------------------------------------------------------------
     

     
      ENDIF !inwp_sfc

    ENDDO  

!$OMP END DO
!$OMP END PARALLEL

    !$acc exit data delete (zml_soil)                                             &
    !$acc           delete (sntunefac, sntunefac2, rain_con_rate, snow_con_rate)  &
    !$acc           delete (rain_gsp_rate, snow_gsp_rate, graupel_gsp_rate)       &
    !$acc           delete (init_list, it1, it2, fact1, fact2, frac_sv)           &
    !$acc           delete (frac_snow_sv),                                        &
    !$acc if(lzacc)

    !
    ! Call seaice parameterization
    !
    IF ( (atm_phy_nwp_config(jg)%inwp_surface == 1) .AND. (lseaice) ) THEN
      
#ifdef _OPENACC
      CALL finish (routine, 'nwp_seaice:  OpenACC version currently not implemented')
#endif

      CALL nwp_seaice(p_patch, p_diag, prm_diag, p_prog_wtr_now, p_prog_wtr_new, &
        &             lnd_prog_now, lnd_prog_new, ext_data, lnd_diag, tcall_sfc_jg)
    ENDIF

    !
    ! Call fresh water lake model (Flake)
    !

    IF ( (atm_phy_nwp_config(jg)%inwp_surface == 1) .AND. (llake) ) THEN
      CALL nwp_lake(p_patch, p_diag, prm_diag, p_prog_wtr_now, p_prog_wtr_new, &
        &           lnd_prog_now, lnd_prog_new, ext_data, lnd_diag, tcall_sfc_jg, lacc=lzacc)
    ENDIF

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Final step: aggregate t_g, qv_s and surface fluxes !!
    !                                                    !!
    ! Loop over all points (land AND water points)       !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,jk,isubs,i_startidx,i_endidx,t_g_s,area_frac)
    !$acc enter data create(t_g_s)
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
        & i_startidx, i_endidx, rl_start, rl_end)

       IF (ntiles_total == 1) THEN 
         !$acc parallel if(lzacc)
         !$acc loop gang vector
         DO jc = i_startidx, i_endidx
           prm_diag%shfl_s (jc,jb)  = prm_diag%shfl_s_t (jc,jb,1) 
           prm_diag%lhfl_s (jc,jb)  = prm_diag%lhfl_s_t (jc,jb,1)
           prm_diag%qhfl_s (jc,jb)  = prm_diag%qhfl_s_t (jc,jb,1)
           prm_diag%lhfl_bs(jc,jb)  = prm_diag%lhfl_bs_t(jc,jb,1) 
           prm_diag%umfl_s (jc,jb)  = prm_diag%umfl_s_t (jc,jb,1)
           prm_diag%vmfl_s (jc,jb)  = prm_diag%vmfl_s_t (jc,jb,1) 
         ENDDO
         !$acc end parallel
         IF (atm_phy_nwp_config(jg)%inwp_surface > 0) THEN
           !$acc parallel if(lzacc)
           !$acc loop gang vector
           DO jc = i_startidx, i_endidx
             lnd_prog_new%t_g(jc,jb)  = lnd_prog_new%t_g_t(jc,jb,1)
             lnd_diag%qv_s   (jc,jb)  = lnd_diag%qv_s_t   (jc,jb,1)
             lnd_diag%h_snow (jc,jb)  = lnd_diag%h_snow_t (jc,jb,1)
           ENDDO
           !$acc loop gang vector collapse(2)
           DO jk=1,nlev_soil
             DO jc = i_startidx, i_endidx
               prm_diag%lhfl_pl(jc,jk,jb)= prm_diag%lhfl_pl_t(jc,jk,jb,1)
             ENDDO  ! jc
           ENDDO  ! jk
           !$acc end parallel
         ENDIF
       ELSE ! aggregate fields over tiles

         !$acc parallel if(lzacc)
         !$acc loop gang vector 
         DO jc = i_startidx, i_endidx
           t_g_s(jc)      = 0._wp
           lnd_diag%qv_s   (jc,jb) = 0._wp
           lnd_diag%h_snow (jc,jb) = 0._wp
           prm_diag%shfl_s (jc,jb) = 0._wp
           prm_diag%lhfl_s (jc,jb) = 0._wp
           prm_diag%qhfl_s (jc,jb) = 0._wp
           prm_diag%umfl_s (jc,jb) = 0._wp
           prm_diag%vmfl_s (jc,jb) = 0._wp
           prm_diag%lhfl_bs(jc,jb) = 0._wp
           DO jk = 1, nlev_soil
             prm_diag%lhfl_pl(jc,jk,jb) = 0._wp
           ENDDO
         ENDDO
         !$acc end parallel

         !$acc parallel if(lzacc)
         !$acc loop seq 
         DO isubs = 1,ntiles_total+ntiles_water
           !$acc loop gang vector private(area_frac)
           DO jc = i_startidx, i_endidx
             area_frac = ext_data%atm%frac_t(jc,jb,isubs)
             t_g_s (jc)           = t_g_s(jc)  + lnd_prog_new%t_g_t(jc,jb,isubs)**4 * area_frac 
             lnd_diag%qv_s(jc,jb) = lnd_diag%qv_s(jc,jb) + lnd_diag%qv_s_t(jc,jb,isubs) * area_frac
             prm_diag%shfl_s(jc,jb) = prm_diag%shfl_s(jc,jb)                    &
               &                    + prm_diag%shfl_s_t (jc,jb,isubs) * area_frac 
             prm_diag%lhfl_s(jc,jb) = prm_diag%lhfl_s(jc,jb)                    &
               &                    + prm_diag%lhfl_s_t (jc,jb,isubs) * area_frac 
             prm_diag%qhfl_s(jc,jb) = prm_diag%qhfl_s(jc,jb)                    &
               &                    + prm_diag%qhfl_s_t (jc,jb,isubs) * area_frac 
             prm_diag%umfl_s(jc,jb) = prm_diag%umfl_s(jc,jb)                    &
               &                    + prm_diag%umfl_s_t (jc,jb,isubs) * area_frac
             prm_diag%vmfl_s(jc,jb) = prm_diag%vmfl_s(jc,jb)                    &
               &                    + prm_diag%vmfl_s_t (jc,jb,isubs) * area_frac 
           ENDDO
         ENDDO
         !$acc end parallel

         !$acc parallel if(lzacc)
         !$acc loop seq  
         DO isubs = 1,ntiles_total
           !$acc loop gang vector private(area_frac)
           DO jc = i_startidx, i_endidx
             ! use rescaled area fraction on mixed land-water points in order to obtain aggregated values
             ! representative for the land part
             area_frac = ext_data%atm%frac_t(jc,jb,isubs)*ext_data%atm%inv_frland_from_tiles(jc,jb)
             prm_diag%lhfl_bs(jc,jb) = prm_diag%lhfl_bs(jc,jb) + prm_diag%lhfl_bs_t(jc,jb,isubs) * area_frac
             lnd_diag%h_snow(jc,jb)  = lnd_diag%h_snow(jc,jb) + lnd_diag%h_snow_t(jc,jb,isubs) * area_frac
           ENDDO  ! jc
           !$acc loop gang vector collapse(2)
           DO jk=1,nlev_soil
             DO jc = i_startidx, i_endidx
               prm_diag%lhfl_pl(jc,jk,jb) = prm_diag%lhfl_pl(jc,jk,jb) + ext_data%atm%frac_t(jc,jb,isubs) &
                 &      * ext_data%atm%inv_frland_from_tiles(jc,jb) * prm_diag%lhfl_pl_t(jc,jk,jb,isubs)
             ENDDO  ! jc
           ENDDO  ! jk
         ENDDO  ! isubs
         !$acc end parallel

         !$acc parallel if(lzacc)
         !$acc loop gang vector
         DO jc = i_startidx, i_endidx
           lnd_prog_new%t_g(jc,jb)  = SQRT(SQRT(t_g_s(jc)))
         ENDDO  ! jc
         !$acc end parallel
       ENDIF    ! with or without tiles

    ENDDO  ! jb
    !$acc exit data delete(t_g_s)
!$OMP END DO
!$OMP END PARALLEL
 
  END SUBROUTINE nwp_surface



  !>
  !! Interface for seaice parameterization
  !!
  !! Interface for seaice parameterization. Calls seaice time integration scheme 
  !! seaice_timestep_nwp and updates the dynamic seaice index lists.
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2012-08-31)
  !!
  !! Modifications by Dmitrii Mironov, DWD (2016-08-08)
  !! - Call to "seaice_timestep_nwp" is modified 
  !!   to allow prognostic treatment of the sea-ice albedo.
  !!
  SUBROUTINE nwp_seaice (p_patch, p_diag, prm_diag, p_prog_wtr_now,  &
    &                    p_prog_wtr_new, lnd_prog_now, lnd_prog_new, &
    &                    ext_data, p_lnd_diag, dtime)

    TYPE(t_patch),        TARGET,INTENT(in)   :: p_patch        !< grid/patch info
    TYPE(t_nh_diag),      TARGET,INTENT(in)   :: p_diag         !< diag vars
    TYPE(t_nwp_phy_diag),        INTENT(in)   :: prm_diag       !< atm phys vars
    TYPE(t_wtr_prog),            INTENT(inout):: p_prog_wtr_now !< prog vars for wtr
    TYPE(t_wtr_prog),            INTENT(inout):: p_prog_wtr_new !< prog vars for wtr
    TYPE(t_lnd_prog),            INTENT(inout):: lnd_prog_now   !< prog vars for sfc
    TYPE(t_lnd_prog),            INTENT(inout):: lnd_prog_new   !< prog vars for sfc
    TYPE(t_external_data),       INTENT(inout):: ext_data       !< external data
    TYPE(t_lnd_diag),            INTENT(inout):: p_lnd_diag     !< diag vars for sfc
    REAL(wp),                    INTENT(in)   :: dtime          !< time interval for 
                                                                !< surface

    ! Local arrays  (local copies)
    !
    REAL(wp) :: shfl_s   (nproma)   ! sensible heat flux at the surface               [W/m^2]
    REAL(wp) :: lhfl_s   (nproma)   ! latent heat flux at the surface                 [W/m^2]
    REAL(wp) :: lwflxsfc (nproma)   ! net long-wave radiation flux at the surface     [W/m^2] 
    REAL(wp) :: swflxsfc (nproma)   ! net solar radiation flux at the surface         [W/m^2]
    REAL(wp) :: snow_rate(nproma)   ! snow rate (convecive + grid-scale)              [kg/(m^2 s)]
    REAL(wp) :: rain_rate(nproma)   ! rain rate (convecive + grid-scale)              [kg/(m^2 s)]
    REAL(wp) :: tice_now (nproma)   ! temperature of ice upper surface at previous time  [K]
    REAL(wp) :: hice_now (nproma)   ! ice thickness at previous time level               [m]
    REAL(wp) :: tsnow_now(nproma)   ! temperature of snow upper surface at previous time [K]
    REAL(wp) :: hsnow_now(nproma)   ! snow thickness at previous time level              [m]
    REAL(wp) :: albsi_now(nproma)   ! sea-ice albedo at previous time level              [-]
    REAL(wp) :: tice_new (nproma)   ! temperature of ice upper surface at new time       [K]
    REAL(wp) :: hice_new (nproma)   ! ice thickness at new time level                    [m]
    REAL(wp) :: tsnow_new(nproma)   ! temperature of snow upper surface at new time      [K]
    REAL(wp) :: hsnow_new(nproma)   ! snow thickness at new time level                   [m]
    REAL(wp) :: albsi_new(nproma)   ! sea-ice albedo at new time level                   [-]

    ! Local array bounds:
    !
    INTEGER :: rl_start, rl_end
    INTEGER :: i_startblk, i_endblk    !> blocks
    INTEGER :: i_nchdom                !< domain index

    ! Local scalars:
    !
    INTEGER :: jc, jb, ic              !loop indices
    INTEGER :: i_count

    CHARACTER(len=*), PARAMETER :: routine = 'mo_nwp_sfc_interface:nwp_seaice'
    !-------------------------------------------------------------------------

    ! exclude nest boundary and halo points
    rl_start = grf_bdywidth_c+1
    rl_end   = min_rlcell_int

    i_nchdom  = MAX(1,p_patch%n_childdom)

    i_startblk = p_patch%cells%start_blk(rl_start,1)
    i_endblk   = p_patch%cells%end_blk(rl_end,i_nchdom)

    IF (msg_level >= 15) THEN
      CALL message(routine, 'call nwp_seaice scheme')
    ENDIF

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_count,ic,jc,shfl_s,lhfl_s,lwflxsfc,swflxsfc,snow_rate,rain_rate, &
!$OMP            tice_now, hice_now,tsnow_now,hsnow_now,tice_new,hice_new,tsnow_new,   &
!$OMP            hsnow_new,albsi_now,albsi_new) ICON_OMP_GUIDED_SCHEDULE
    DO jb = i_startblk, i_endblk

      !
      ! Copy input fields
      !
      i_count = ext_data%atm%list_seaice%ncount(jb)


      IF (i_count == 0) CYCLE ! skip loop if the index list for the given block is empty

      DO ic = 1, i_count
        jc = ext_data%atm%list_seaice%idx(ic,jb)

        shfl_s   (ic) = prm_diag%shfl_s_t  (jc,jb,isub_seaice)   ! sensible heat flux at sfc    [W/m^2]
        lhfl_s   (ic) = prm_diag%lhfl_s_t  (jc,jb,isub_seaice)   ! latent heat flux at sfc      [W/m^2]
        lwflxsfc (ic) = prm_diag%lwflxsfc_t(jc,jb,isub_seaice)   ! net lw radiation flux at sfc [W/m^2]
        swflxsfc (ic) = prm_diag%swflxsfc_t(jc,jb,isub_seaice)   ! net solar radiation flux at sfc [W/m^2]
        snow_rate(ic) = prm_diag%snow_gsp_rate(jc,jb)  &         ! snow rate (convecive + grid-scale) [kg/(m^2 s)]
          &           + prm_diag%snow_con_rate(jc,jb)
        rain_rate(ic) = prm_diag%rain_gsp_rate(jc,jb)  &         !  rain rate (convecive + grid-scale) [kg/(m^2 s)]
          &           + prm_diag%rain_con_rate(jc,jb)
        tice_now (ic) = p_prog_wtr_now%t_ice    (jc,jb)
        hice_now (ic) = p_prog_wtr_now%h_ice    (jc,jb)
        tsnow_now(ic) = p_prog_wtr_now%t_snow_si(jc,jb)
        hsnow_now(ic) = p_prog_wtr_now%h_snow_si(jc,jb)
        albsi_now(ic) = p_prog_wtr_now%alb_si(jc,jb)             ! sea-ice albedo [-]
      ENDDO  ! ic


      ! call seaice time integration scheme
      !
      CALL seaice_timestep_nwp (                               &
                            &   dtime   = dtime,               &
                            &   nsigb   = i_count,             & !in
                            &   qsen    = shfl_s(:),           & !in 
                            &   qlat    = lhfl_s(:),           & !in
                            &   qlwrnet = lwflxsfc(:),         & !in
                            &   qsolnet = swflxsfc(:),         & !in
                            &   snow_rate = snow_rate(:),      & !in
                            &   rain_rate = rain_rate(:),      & !in
                            &   tice_p  = tice_now(:),         & !in
                            &   hice_p  = hice_now(:),         & !in
                            &   tsnow_p = tsnow_now(:),        & !in    ! DUMMY: not used yet
                            &   hsnow_p = hsnow_now(:),        & !in    ! DUMMY: not used yet
                            &   albsi_p = albsi_now(:),        & !in 
                            &   tice_n  = tice_new(:),         & !out
                            &   hice_n  = hice_new(:),         & !out
                            &   tsnow_n = tsnow_new(:),        & !out   ! DUMMY: not used yet
                            &   hsnow_n = hsnow_new(:),        & !out   ! DUMMY: not used yet
                            &   albsi_n = albsi_new(:)         ) !out
! optional arguments dticedt, dhicedt, dtsnowdt, dhsnowdt (tendencies) are neglected



      !  Recover fields from index list
      !
!$NEC ivdep
      DO ic = 1, i_count
        jc = ext_data%atm%list_seaice%idx(ic,jb)

        p_prog_wtr_new%t_ice(jc,jb)     = tice_new(ic)
        p_prog_wtr_new%h_ice(jc,jb)     = hice_new(ic)
        p_prog_wtr_new%t_snow_si(jc,jb) = tsnow_new(ic)
        p_prog_wtr_new%h_snow_si(jc,jb) = hsnow_new(ic)
        IF (lprog_albsi) THEN
          p_prog_wtr_new%alb_si(jc,jb)  = albsi_new(ic)
        ENDIF

        lnd_prog_new%t_g_t(jc,jb,isub_seaice) = tice_new(ic)
        ! surface saturation specific humidity (uses saturation water vapor pressure 
        ! over ice)
        p_lnd_diag%qv_s_t(jc,jb,isub_seaice) = spec_humi(sat_pres_ice(tice_new(ic)),&
          &                                   p_diag%pres_sfc(jc,jb) )
      ENDDO  ! ic


      ! Update dynamic sea-ice index list
      !
      CALL update_idx_lists_sea (                                                 &
        &              hice_n           = p_prog_wtr_new%h_ice(:,jb),             &!in
        &              pres_sfc         = p_diag%pres_sfc(:,jb),                  &!in
        &              list_seawtr_idx  = ext_data%atm%list_seawtr%idx(:,jb),     &!inout
        &              list_seawtr_count= ext_data%atm%list_seawtr%ncount(jb),    &!inout
        &              list_seaice_idx  = ext_data%atm%list_seaice%idx(:,jb),     &!inout
        &              list_seaice_count= ext_data%atm%list_seaice%ncount(jb),    &!inout
        &              frac_t_ice       = ext_data%atm%frac_t(:,jb,isub_seaice),  &!inout
        &              frac_t_water     = ext_data%atm%frac_t(:,jb,isub_water),   &!inout
        &              lc_frac_t_water  = ext_data%atm%lc_frac_t(:,jb,isub_water),&!inout
        &              fr_seaice        = p_lnd_diag%fr_seaice(:,jb),             &!inout
        &              hice_old         = p_prog_wtr_now%h_ice(:,jb),             &!inout
        &              tice_old         = p_prog_wtr_now%t_ice(:,jb),             &!inout
        &              albsi_now        = p_prog_wtr_now%alb_si(:,jb),            &!inout
        &              albsi_new        = p_prog_wtr_new%alb_si(:,jb),            &!inout
        &              t_g_t_now        = lnd_prog_now%t_g_t(:,jb,isub_water),    &!inout
        &              t_g_t_new        = lnd_prog_new%t_g_t(:,jb,isub_water),    &!inout
        &              t_s_t_now        = lnd_prog_now%t_s_t(:,jb,isub_water),    &!inout
        &              t_s_t_new        = lnd_prog_new%t_s_t(:,jb,isub_water),    &!inout
        &              t_sk_t_now       = lnd_prog_now%t_sk_t(:,jb,isub_water),   &!inout
        &              t_sk_t_new       = lnd_prog_new%t_sk_t(:,jb,isub_water),   &!inout
        &              qv_s_t           = p_lnd_diag%qv_s_t(:,jb,isub_water),     &!inout
        &              t_seasfc         = p_lnd_diag%t_seasfc(:,jb)               )!inout

    ENDDO  ! jb
!$OMP END DO
!$OMP END PARALLEL

  END SUBROUTINE nwp_seaice



  !>
  !! Interface for fresh water lake (Flake) parameterization
  !!
  !! Interface for fresh water lake (Flake) parameterization. Calls time 
  !! integration routine flake_interface and updates the prognostic Flake variables 
  !! as well as t_g_t and qv_s_t.
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2013-06-26)
  !!
  SUBROUTINE nwp_lake (p_patch, p_diag, prm_diag, p_prog_wtr_now,    &
    &                    p_prog_wtr_new, lnd_prog_now, lnd_prog_new, &
    &                    ext_data, p_lnd_diag, dtime, lacc)

    TYPE(t_patch),        TARGET,INTENT(in)   :: p_patch        !< grid/patch info
    TYPE(t_nh_diag),      TARGET,INTENT(in)   :: p_diag         !< diag vars
    TYPE(t_nwp_phy_diag),        INTENT(in)   :: prm_diag       !< atm phys vars
    TYPE(t_wtr_prog),            INTENT(in)   :: p_prog_wtr_now !< prog vars for wtr
    TYPE(t_wtr_prog),            INTENT(inout):: p_prog_wtr_new !< prog vars for wtr
    TYPE(t_lnd_prog),            INTENT(in)   :: lnd_prog_now   !< prog vars for sfc
    TYPE(t_lnd_prog),            INTENT(inout):: lnd_prog_new   !< prog vars for sfc
    TYPE(t_external_data),       INTENT(in)   :: ext_data       !< external data
    TYPE(t_lnd_diag),            INTENT(inout):: p_lnd_diag     !< diag vars for sfc
    REAL(wp),                    INTENT(in)   :: dtime          !< time interval for 
    LOGICAL, OPTIONAL,           INTENT(IN)   :: lacc           !< openACC flag
                                                                !< surface

    ! Local arrays  (local copies)
    !
    REAL(wp) :: f_c      (nproma)
    REAL(wp) :: depth_lk (nproma)
    REAL(wp) :: fetch_lk (nproma)
    REAL(wp) :: dp_bs_lk (nproma)
    REAL(wp) :: t_bs_lk  (nproma)
    REAL(wp) :: gamso_lk (nproma)
    REAL(wp) :: qmom     (nproma)
    REAL(wp) :: shfl_s   (nproma)
    REAL(wp) :: lhfl_s   (nproma)
    REAL(wp) :: swflxsfc (nproma)
    REAL(wp) :: lwflxsfc (nproma)
    REAL(wp) :: t_snow_lk_now(nproma), t_snow_lk_new(nproma)
    REAL(wp) :: h_snow_lk_now(nproma), h_snow_lk_new(nproma)
    REAL(wp) :: t_ice_now(nproma), t_ice_new(nproma)
    REAL(wp) :: h_ice_now(nproma), h_ice_new(nproma)
    REAL(wp) :: t_mnw_lk_now(nproma), t_mnw_lk_new(nproma)
    REAL(wp) :: t_wml_lk_now(nproma), t_wml_lk_new(nproma)
    REAL(wp) :: t_bot_lk_now(nproma), t_bot_lk_new(nproma)
    REAL(wp) :: c_t_lk_now(nproma), c_t_lk_new(nproma)
    REAL(wp) :: h_ml_lk_now(nproma), h_ml_lk_new(nproma)
    REAL(wp) :: t_b1_lk_now(nproma), t_b1_lk_new(nproma)
    REAL(wp) :: h_b1_lk_now(nproma), h_b1_lk_new(nproma)
    REAL(wp) :: t_scf_lk_now(nproma), t_scf_lk_new(nproma)


    ! Local array bounds:
    !
    INTEGER :: rl_start, rl_end
    INTEGER :: i_startblk, i_endblk    !> blocks
    INTEGER :: i_nchdom                !< domain index

    ! Local scalars:
    !
    INTEGER :: jc, jb, ic              !loop indices
    INTEGER :: icount_flk

    ! openACC flag
    !
    LOGICAL :: lzacc

    ! routine name
    !
    CHARACTER(len=*), PARAMETER :: routine = 'mo_nwp_sfc_interface:nwp_lake'
    !-------------------------------------------------------------------------

    IF(PRESENT(lacc)) THEN
      lzacc = lacc
    ELSE
      lzacc = .FALSE.
    ENDIF

    ! put local variables on gpu
    !$acc enter data &
    !$acc create (f_c, depth_lk, fetch_lk, dp_bs_lk, t_bs_lk, gamso_lk, qmom, shfl_s, lhfl_s) &
    !$acc create (swflxsfc, lwflxsfc, t_snow_lk_now, h_snow_lk_now, t_ice_now, h_ice_now, t_mnw_lk_now) &
    !$acc create (t_wml_lk_now, t_bot_lk_now, c_t_lk_now, h_ml_lk_now, t_b1_lk_now, h_b1_lk_now, t_scf_lk_now) &
    !$acc create (t_snow_lk_new, h_snow_lk_new,t_ice_new,h_ice_new,t_mnw_lk_new,t_wml_lk_new) &
    !$acc create (t_bot_lk_new, c_t_lk_new, h_ml_lk_new, t_b1_lk_new, h_b1_lk_new, t_scf_lk_new) if (lzacc)

    ! exclude nest boundary and halo points
    rl_start = grf_bdywidth_c+1
    rl_end   = min_rlcell_int

    i_nchdom  = MAX(1,p_patch%n_childdom)

    i_startblk = p_patch%cells%start_blk(rl_start,1)
    i_endblk   = p_patch%cells%end_blk(rl_end,i_nchdom)

    IF (msg_level >= 15) THEN
      CALL message(routine, 'call nwp_lake scheme')
    ENDIF

    

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,ic,jc,icount_flk,f_c,depth_lk,fetch_lk,dp_bs_lk,t_bs_lk,  &
!$OMP            gamso_lk,qmom,shfl_s,lhfl_s,swflxsfc,lwflxsfc,t_snow_lk_now, &
!$OMP            h_snow_lk_now,t_ice_now,h_ice_now,t_mnw_lk_now,              &
!$OMP            t_wml_lk_now,t_bot_lk_now,c_t_lk_now,h_ml_lk_now,t_b1_lk_now,&
!$OMP            h_b1_lk_now,t_scf_lk_now,t_snow_lk_new,h_snow_lk_new,        &
!$OMP            t_ice_new,h_ice_new,t_mnw_lk_new,t_wml_lk_new,               &
!$OMP            t_bot_lk_new,c_t_lk_new,h_ml_lk_new,t_b1_lk_new,h_b1_lk_new, &
!$OMP            t_scf_lk_new) ICON_OMP_GUIDED_SCHEDULE
    DO jb = i_startblk, i_endblk

      !
      ! Copy input fields
      !
      icount_flk = ext_data%atm%list_lake%ncount(jb) 

      ! Collect data for lake points in 1D-arrays

      !$acc parallel default (present) if (lzacc)
      !$acc loop gang vector private (jc)
      DO ic=1,icount_flk

        jc = ext_data%atm%list_lake%idx(ic,jb)

        f_c      (ic) = p_patch%cells%f_c     (jc,jb)    ! Coriolis parameter   [s^-1]
 
        depth_lk (ic) = ext_data%atm%depth_lk (jc,jb)    ! lake depth           [m]
        fetch_lk (ic) = ext_data%atm%fetch_lk (jc,jb)    ! wind fetch over lake [m]
        dp_bs_lk (ic) = ext_data%atm%dp_bs_lk (jc,jb)
        t_bs_lk  (ic) = ext_data%atm%t_bs_lk  (jc,jb)
        gamso_lk (ic) = ext_data%atm%gamso_lk (jc,jb)
        
        ! absolute value of momentum flux at sfc
        qmom     (ic) = SQRT(prm_diag%umfl_s_t(jc,jb,isub_lake)**2  &
          &             +    prm_diag%vmfl_s_t(jc,jb,isub_lake)**2 )
        shfl_s   (ic) = prm_diag%shfl_s_t  (jc,jb,isub_lake)   ! sensible heat flux at sfc [W/m^2]
        lhfl_s   (ic) = prm_diag%lhfl_s_t  (jc,jb,isub_lake)   ! latent heat flux at sfc   [W/m^2]
        swflxsfc (ic) = prm_diag%swflxsfc_t(jc,jb,isub_lake)   ! net shortwave flux at sfc [W/m^2]
        lwflxsfc (ic) = prm_diag%lwflxsfc_t(jc,jb,isub_lake)   ! net longwave flux at sfc  [W/m^2]

        t_snow_lk_now(ic) = p_prog_wtr_now%t_snow_lk(jc,jb)
        h_snow_lk_now(ic) = p_prog_wtr_now%h_snow_lk(jc,jb)
        t_ice_now    (ic) = p_prog_wtr_now%t_ice    (jc,jb)    ! ice temperature
        h_ice_now    (ic) = p_prog_wtr_now%h_ice    (jc,jb)    ! ice depth
        t_mnw_lk_now (ic) = p_prog_wtr_now%t_mnw_lk (jc,jb)
        t_wml_lk_now (ic) = p_prog_wtr_now%t_wml_lk (jc,jb)
        t_bot_lk_now (ic) = p_prog_wtr_now%t_bot_lk (jc,jb)
        c_t_lk_now   (ic) = p_prog_wtr_now%c_t_lk   (jc,jb)
        h_ml_lk_now  (ic) = p_prog_wtr_now%h_ml_lk  (jc,jb)
        t_b1_lk_now  (ic) = p_prog_wtr_now%t_b1_lk  (jc,jb)
        h_b1_lk_now  (ic) = p_prog_wtr_now%h_b1_lk  (jc,jb)
        t_scf_lk_now (ic) = lnd_prog_now%t_g_t      (jc,jb,isub_lake) ! only required to compute the time
                                                                      ! tendency of the lake surface temperature,
                                                                      ! which is omitted so far 
                                                                      ! (optional arg of flake_interface) 
      ENDDO
      !$acc end parallel
      
      CALL flake_interface (                                  & !in
                     &  dtime       = dtime           ,       & !in
                     &  nflkgb      = icount_flk      ,       & !in
                     &  coriolispar = f_c          (:),       & !in
                     &  depth_lk    = depth_lk     (:),       & !in
                     &  fetch_lk    = fetch_lk     (:),       & !in
                     &  dp_bs_lk    = dp_bs_lk     (:),       & !in
                     &  t_bs_lk     = t_bs_lk      (:),       & !in
                     &  gamso_lk    = gamso_lk     (:),       & !in
                     &  qmom        = qmom         (:),       & !in
                     &  qsen        = shfl_s       (:),       & !in
                     &  qlat        = lhfl_s       (:),       & !in
                     &  qlwrnet     = lwflxsfc     (:),       & !in
                     &  qsolnet     = swflxsfc     (:),       & !in
                     &  t_snow_p    = t_snow_lk_now(:),       & !in
                     &  h_snow_p    = h_snow_lk_now(:),       & !in
                     &  t_ice_p     = t_ice_now    (:),       & !in
                     &  h_ice_p     = h_ice_now    (:),       & !in
                     &  t_mnw_lk_p  = t_mnw_lk_now (:),       & !in
                     &  t_wml_lk_p  = t_wml_lk_now (:),       & !in
                     &  t_bot_lk_p  = t_bot_lk_now (:),       & !in
                     &  c_t_lk_p    = c_t_lk_now   (:),       & !in
                     &  h_ml_lk_p   = h_ml_lk_now  (:),       & !in
                     &  t_b1_lk_p   = t_b1_lk_now  (:),       & !in
                     &  h_b1_lk_p   = h_b1_lk_now  (:),       & !in       
                     &  t_scf_lk_p  = t_scf_lk_now (:),       & !in 
                     &  t_snow_n    = t_snow_lk_new(:),       & !out
                     &  h_snow_n    = h_snow_lk_new(:),       & !out
                     &  t_ice_n     = t_ice_new    (:),       & !out
                     &  h_ice_n     = h_ice_new    (:),       & !out
                     &  t_mnw_lk_n  = t_mnw_lk_new (:),       & !out
                     &  t_wml_lk_n  = t_wml_lk_new (:),       & !out
                     &  t_bot_lk_n  = t_bot_lk_new (:),       & !out
                     &  c_t_lk_n    = c_t_lk_new   (:),       & !out
                     &  h_ml_lk_n   = h_ml_lk_new  (:),       & !out
                     &  t_b1_lk_n   = t_b1_lk_new  (:),       & !out
                     &  h_b1_lk_n   = h_b1_lk_new  (:),       & !out
                     &  t_scf_lk_n  = t_scf_lk_new (:),       & !out
                     &  lacc        = lzacc                   ) !in openACC flag
! optional arguments (tendencies) are omitted


      !  Recover fields from index list
      !
!$NEC ivdep

      !$acc parallel default (present) if (lzacc)
      !$acc loop gang vector
      DO ic = 1,icount_flk
        jc = ext_data%atm%list_lake%idx(ic,jb)

        p_prog_wtr_new%t_snow_lk(jc,jb)     = t_snow_lk_new(ic)
        p_prog_wtr_new%h_snow_lk(jc,jb)     = h_snow_lk_new(ic)
        p_prog_wtr_new%t_ice    (jc,jb)     = t_ice_new    (ic)
        p_prog_wtr_new%h_ice    (jc,jb)     = h_ice_new    (ic)
        p_prog_wtr_new%t_mnw_lk (jc,jb)     = t_mnw_lk_new (ic)
        p_prog_wtr_new%t_wml_lk (jc,jb)     = t_wml_lk_new (ic)
        p_prog_wtr_new%t_bot_lk (jc,jb)     = t_bot_lk_new (ic)
        p_prog_wtr_new%c_t_lk   (jc,jb)     = c_t_lk_new   (ic)
        p_prog_wtr_new%h_ml_lk  (jc,jb)     = h_ml_lk_new  (ic)
        p_prog_wtr_new%t_b1_lk  (jc,jb)     = t_b1_lk_new  (ic)
        p_prog_wtr_new%h_b1_lk  (jc,jb)     = h_b1_lk_new  (ic)

        lnd_prog_new%t_g_t(jc,jb,isub_lake) = t_scf_lk_new (ic)

        ! for consistency, set 
        ! t_so(0) = t_wml_lk       mixed-layer temperature (273.15K if the lake is frozen)
        lnd_prog_new%t_s_t (jc,jb,isub_lake) = p_prog_wtr_new%t_wml_lk (jc,jb)

        ! surface saturation specific humidity over water/ice 
        !
        IF ( h_ice_new (ic) > h_Ice_min_flk ) THEN
          p_lnd_diag%qv_s_t(jc,jb,isub_lake)  = spec_humi(sat_pres_ice(t_scf_lk_new(ic)),&
            &                                   p_diag%pres_sfc(jc,jb) )
          ! keep fr_seaice synchronized with h_ice
          p_lnd_diag%fr_seaice(jc,jb) = 1._wp
        ELSE
          p_lnd_diag%qv_s_t(jc,jb,isub_lake)  = spec_humi(sat_pres_water(t_scf_lk_new(ic)),&
            &                                   p_diag%pres_sfc(jc,jb) )
          ! keep fr_seaice synchronized with h_ice
          p_lnd_diag%fr_seaice(jc,jb) = 0._wp
        ENDIF

      ENDDO  ! ic
      !$acc end parallel

    ENDDO  !jb
!$OMP END DO
!$OMP END PARALLEL

! remove local variables from gpu
!$acc exit data &
!$acc delete (f_c, depth_lk, fetch_lk, dp_bs_lk, t_bs_lk, gamso_lk, qmom, shfl_s, lhfl_s) &
!$acc delete (swflxsfc, lwflxsfc, t_snow_lk_now, h_snow_lk_now, t_ice_now, h_ice_now, t_mnw_lk_now) &
!$acc delete (t_wml_lk_now, t_bot_lk_now, c_t_lk_now, h_ml_lk_now, t_b1_lk_now, h_b1_lk_now, t_scf_lk_now) &
!$acc delete (t_snow_lk_new, h_snow_lk_new,t_ice_new,h_ice_new,t_mnw_lk_new,t_wml_lk_new) &
!$acc delete (t_bot_lk_new, c_t_lk_new, h_ml_lk_new, t_b1_lk_new, h_b1_lk_new, t_scf_lk_new) if (lzacc)

  END SUBROUTINE nwp_lake

END MODULE mo_nwp_sfc_interface

