!OPTION! -cont -msg o
!! this command should fix the problem of copying arrays in a subroutine call
!>
!! This module is the interface between nwp_nh_interface to the
!! turbulence parameterisations:
!! inwp_turb == 1 == turbulence scheme by M. Raschendorfer run in COSMO
!! inwp_turb == 2 == turbulence scheme imported from the GME
!! This module handles the computation of surface transfer coefficients, only.
!!
!! @author Kristina Froehlich, DWD, Offenbach (2010-01-25)
!!
!! @par Revision History
!! Initial Kristina Froehlich, DWD, Offenbach (2010-01-25)
!! Added gust calculation by Helmut Frank, DWD, Offenbach (2013-03-13)
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
#if defined __xlC__
@PROCESS SPILL(988)
#endif
MODULE mo_nwp_turbtrans_interface

  USE mo_kind,                 ONLY: wp
  USE mo_exception,            ONLY: message, finish, warning
  USE mo_model_domain,         ONLY: t_patch
  USE mo_impl_constants,       ONLY: min_rlcell_int, icosmo, igme, ismag, iedmf, iprog
  USE mo_impl_constants_grf,   ONLY: grf_bdywidth_c
  USE mo_loopindices,          ONLY: get_indices_c
  USE mo_physical_constants,   ONLY: rd_o_cpd, grav, lh_v=>alv, lh_s=>als, rd, cpd
  USE mo_ext_data_types,       ONLY: t_external_data
  USE mo_nonhydro_types,       ONLY: t_nh_prog, t_nh_diag, t_nh_metrics
  USE mo_nwp_phy_types,        ONLY: t_nwp_phy_diag
  USE mo_nwp_phy_state,        ONLY: phy_params
  USE mo_nwp_lnd_types,        ONLY: t_lnd_prog, t_wtr_prog, t_lnd_diag
  USE mo_parallel_config,      ONLY: nproma
  USE mo_run_config,           ONLY: msg_level, iqv, iqc, iqtke
  USE mo_atm_phy_nwp_config,   ONLY: atm_phy_nwp_config
  USE mo_advection_config,     ONLY: advection_config
  USE turb_data,               ONLY: get_turbdiff_param
  USE sfc_flake_data,          ONLY: h_Ice_min_flk, tpl_T_f
  USE turb_transfer,           ONLY: turbtran
  USE mo_satad,                ONLY: sat_pres_water, spec_humi
  USE mo_gme_turbdiff,         ONLY: parturs, nearsfc
  USE mo_util_phys,            ONLY: nwp_dyn_gust
  USE mo_run_config,           ONLY: ltestcase
  USE mo_nh_testcases_nml,     ONLY: nh_test_name
  USE mo_lnd_nwp_config,       ONLY: ntiles_total, ntiles_lnd, ntiles_water, llake,  &
    &                                isub_lake, lseaice
  USE mo_vupdz0_tile,          ONLY: vupdz0_tile
  USE mo_vexcs,                ONLY: vexcs

!$ser verbatim USE mo_ser_nwp_tutra, ONLY: serialize_turbtrans_interface_input,&
!$ser verbatim                             serialize_turbtrans_interface_output

  IMPLICIT NONE

  PRIVATE



  PUBLIC  ::  nwp_turbtrans

CONTAINS
  !!
  !!-------------------------------------------------------------------------
  !!
SUBROUTINE nwp_turbtrans  ( tcall_turb_jg,                     & !>in
                          & p_patch, p_metrics,                & !>in
                          & ext_data,                          & !>in
                          & p_prog,                            & !>in
                          & p_prog_rcf,                        & !>inout
                          & p_diag ,                           & !>inout
                          & prm_diag,                          & !>inout
                          & wtr_prog_new,                      & !>in
                          & lnd_prog_new,                      & !>inout
                          & lnd_diag,                          & !>inout
                          & lacc                               ) !>in


  TYPE(t_patch),        TARGET,INTENT(in)   :: p_patch        !!<grid/patch info.
  TYPE(t_external_data),TARGET,INTENT(in)   :: ext_data        !< external data
  TYPE(t_nh_metrics)          ,INTENT(in)   :: p_metrics
  TYPE(t_nh_prog),      TARGET,INTENT(inout):: p_prog          !<the prog vars
  TYPE(t_nh_prog)             ,INTENT(inout):: p_prog_rcf      !< current time levels
  TYPE(t_nh_diag),      TARGET,INTENT(inout):: p_diag          !<the diag vars
  TYPE(t_nwp_phy_diag),        INTENT(inout):: prm_diag        !< atm phys vars
  TYPE(t_wtr_prog),            INTENT(in)   :: wtr_prog_new    !< prog vars for wtr
  TYPE(t_lnd_prog),            INTENT(inout):: lnd_prog_new    !< prog vars for sfc
  TYPE(t_lnd_diag),            INTENT(inout):: lnd_diag        !< diag vars for sfc
  REAL(wp),                    INTENT(in)   :: tcall_turb_jg   !< time interval for
                                                               !< turbulence
  LOGICAL, OPTIONAL,           INTENT(in)   :: lacc            !< GPU flag
  LOGICAL :: lzacc

  CHARACTER(len=*),PARAMETER :: routine = 'mo_nwp_turbtrans_interface:nwp_turbtrans'

  ! Local array bounds

  INTEGER :: rl_start, rl_end
  INTEGER :: i_startblk, i_endblk    !> blocks
  INTEGER :: i_startidx, i_endidx    !< slices
  INTEGER :: nzprv                   !< argument to turbtran

  ! Local scalars:

  INTEGER :: jc,jb,jk,jt,jg,ic,i_count      !loop indices
  INTEGER :: jk_gust(nproma)

  ! local variables for turbdiff

  INTEGER :: ierrstat
  CHARACTER (LEN=25) :: eroutine=''
  CHARACTER (LEN=80) :: errormsg=''

  INTEGER  :: nlev, nlevp1, nlevcm                  !< number of full, half and canopy levels
  INTEGER  :: lc_class                              !< land-cover class

  ! SQRT(2*TKE) turbulence velocity scale [m/s]
  REAL(wp) :: z_tvs(nproma,3,1), tvs_t(nproma,3,1,ntiles_total+ntiles_water)

  REAL(wp) :: fr_land_t(nproma),depth_lk_t(nproma),h_ice_t(nproma),area_frac,z0_mod,fact_z0rough

  ! Local fields needed to reorder turbtran input/output fields for tile approach

  ! 1D fields
  REAL(wp), DIMENSION(nproma)   :: pres_sfc_t, l_hori

  ! 2D half-level fields
  REAL(wp), DIMENSION(nproma,3) :: z_ifc_t
  REAL(wp), DIMENSION(nproma,3,ntiles_total+ntiles_water) :: rcld_t

  ! 2D full level fields (no tiles)
  REAL(wp), DIMENSION(nproma,2) :: u_t, v_t, temp_t, pres_t, qv_t, qc_t, epr_t

  ! 3D full-level fields (tiles)
  REAL(wp), DIMENSION(nproma,3,ntiles_total+ntiles_water) :: tkvm_t, tkvh_t

  ! 2D fields (tiles)
  REAL(wp), DIMENSION(nproma,ntiles_total+ntiles_water) :: &
   gz0_t, tcm_t, tch_t, tfm_t, tfh_t, tfv_t, tvm_t, tvh_t, tkr_t, &
   t_2m_t, qv_2m_t, td_2m_t, rh_2m_t, u_10m_t, v_10m_t, t_g_t, qv_s_t, sai_t, shfl_s_t,  &
   lhfl_s_t, qhfl_s_t, umfl_s_t, vmfl_s_t

  REAL(wp), DIMENSION(nproma) :: &
   PGEOMLEV,PCPTGZLEV,PCPTSTI,PUCURR,PVCURR,ZCFMTI,PCFHTI,PCFQTI,ZBUOMTI,ZZDLTI,  &
   ZZ0MTI,ZZ0HTI,ZZ0QTI
  REAL(wp) :: rho_s

  INTEGER,  POINTER :: ilist(:)       ! pointer to tile index list

!--------------------------------------------------------------

  IF(PRESENT(lacc)) THEN
      lzacc = lacc
  ELSE
      lzacc = .FALSE.
  ENDIF

  IF (msg_level >= 15) CALL message('mo_nwp_turbtrans_interface:', 'turbulence')

  ! number of vertical levels
  nlev   = p_patch%nlev
  nlevp1 = p_patch%nlevp1

  ! domain
  jg = p_patch%id

!$acc data create(gz0_t, tcm_t, tch_t, tfm_t, tfh_t, tfv_t, tvm_t, tvh_t, tkr_t, t_2m_t, qv_2m_t, td_2m_t, rh_2m_t, &
!$acc             u_10m_t, v_10m_t, t_g_t, qv_s_t, sai_t, shfl_s_t, lhfl_s_t, qhfl_s_t, umfl_s_t, vmfl_s_t, &
!$acc             tkvm_t, tkvh_t, u_t, v_t, temp_t, pres_t, qv_t, qc_t, epr_t, rcld_t, z_ifc_t, pres_sfc_t, l_hori, &
!$acc             z_tvs, tvs_t, fr_land_t, depth_lk_t, h_ice_t, jk_gust) &
!$acc      if(lzacc)

  ! exclude boundary interpolation zone of nested domains
  rl_start = grf_bdywidth_c+1
  rl_end   = min_rlcell_int

  i_startblk = p_patch%cells%start_block(rl_start)
  i_endblk   = p_patch%cells%end_block(rl_end)


  IF ( ANY( (/icosmo,iedmf/)==atm_phy_nwp_config(jg)%inwp_turb ) ) THEN
     CALL get_turbdiff_param(jg)
  ENDIF

  ! Scaling factor for SSO contribution to roughness length ("Erdmann Heise formula")
  fact_z0rough = 1.e-5_wp*ATAN(phy_params(jg)%mean_charlen/2250._wp)


  ! Serialbox2 input fields serialization
  !$ser verbatim IF(lzacc) THEN
  !$ser verbatim     CALL serialize_turbtrans_interface_input(jg, nproma, nlev,&
  !$ser verbatim                                              p_metrics, p_prog, p_prog_rcf, p_diag, prm_diag,&
  !$ser verbatim                                              lnd_prog_new, lnd_diag, wtr_prog_new, ext_data)
  !$ser verbatim ENDIF


!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jt,jc,jk,ic,ilist,i_startidx,i_endidx,i_count,ierrstat,errormsg,eroutine,   &
!$OMP nzprv,lc_class,z_tvs,z0_mod,gz0_t,tcm_t,tch_t,tfm_t,tfh_t,tfv_t,tvm_t,tvh_t,tkr_t,l_hori, &
!$OMP t_g_t,qv_s_t,t_2m_t,qv_2m_t,td_2m_t,rh_2m_t,u_10m_t,v_10m_t,tvs_t,pres_sfc_t,u_t,v_t,     &
!$OMP temp_t,pres_t,qv_t,qc_t,tkvm_t,tkvh_t,z_ifc_t,rcld_t,sai_t,fr_land_t,depth_lk_t,h_ice_t,  &
!$OMP area_frac,shfl_s_t,lhfl_s_t,qhfl_s_t,umfl_s_t,vmfl_s_t,nlevcm,jk_gust,epr_t,              &
!$OMP PGEOMLEV,PCPTGZLEV,PCPTSTI,PUCURR,PVCURR,ZCFMTI,PCFHTI,PCFQTI,ZBUOMTI,ZZDLTI,             &
!$OMP ZZ0MTI,ZZ0HTI,ZZ0QTI,rho_s ) ICON_OMP_GUIDED_SCHEDULE
!MR:>

  DO jb = i_startblk, i_endblk

    CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
      & i_startidx, i_endidx, rl_start, rl_end)

   !-------------------------------------------------------------------------
   !<  turbulent transfer
   !-------------------------------------------------------------------------

   !<  NOTE: since  turbulence is a fast process it is
   !!        allowed to do a sequential updating except for wind speed
   !!        because back-and-forth interpolation would cause too large errors
   !!  (GZ, 2011-08-29): Nevertheless, tendency fields are now passed to turbdiff
   !!        to have them available for extended diagnostic output


    IF( atm_phy_nwp_config(jg)%inwp_surface == 0) THEN
      ! check dry case
      IF( atm_phy_nwp_config(jg)%inwp_satad == 0) THEN
        lnd_diag%qv_s (:,jb) = 0._wp
      ELSE IF ( ANY( (/icosmo,igme,iedmf/)==atm_phy_nwp_config(jg)%inwp_turb ) ) THEN
        IF ( ltestcase .AND. nh_test_name == 'wk82') THEN

!DR Note that this must be re-checked, once turbtran is called at the very end
!DR of the fast physics part.
!DIR$ IVDEP
         !$acc parallel default(present) if(lzacc)
         !$acc loop gang vector
         DO jc = i_startidx, i_endidx
          lnd_prog_new%t_g(jc,jb) = p_diag%temp(jc,nlev,jb)*  &
                      ((p_diag%pres_sfc(jc,jb))/p_diag%pres(jc,nlev,jb))**rd_o_cpd
          lnd_diag%qv_s (jc,jb) = &
             &         spec_humi(sat_pres_water(lnd_prog_new%t_g(jc,jb)),&
             &                                   p_diag%pres_sfc(jc,jb) )
          lnd_diag%qv_s(jc,jb) = MIN (lnd_diag%qv_s(jc,jb) ,p_prog_rcf%tracer(jc,nlev,jb,iqv))
         END DO
         !$acc end parallel
        ELSE
         !
         !> adjust humidity at water surface because of changed surface pressure
         !
!DIR$ IVDEP
         !$acc parallel default(present) if(lzacc)
         !$acc loop gang vector
         DO jc = i_startidx, i_endidx
           lnd_diag%qv_s (jc,jb) = &
             &         spec_humi(sat_pres_water(lnd_prog_new%t_g(jc,jb)),&
             &                                   p_diag%pres_sfc(jc,jb) )
         ENDDO
         !$acc end parallel
        END IF
      ENDIF
    ELSE IF (atm_phy_nwp_config(jg)%itype_z0 >= 2) THEN
      ! specify land-cover-related roughness length over land points
      ! NOTE:  open water, lake and sea-ice points are set in turbtran
      DO jt = 1, ntiles_total
!$NEC ivdep
        !$acc parallel default(present) if(lzacc)
        !$acc loop gang vector private(jc,lc_class,z0_mod)
        DO ic = 1, ext_data%atm%gp_count_t(jb,jt)
          ! works for the following two cases
          ! 1. snow-covered and snow_free tiles treated separately
          ! 2. snow_covered and snow_free areas combined in one tile
          jc = ext_data%atm%idx_lst_t(ic,jb,jt)
          lc_class = MAX(1,ext_data%atm%lc_class_t(jc,jb,jt)) ! to avoid segfaults
          ! Reduction of land-cover related roughness length when the vegetation is below 75%
          IF (ext_data%atm%z0_lcc(lc_class) >= 0.5_wp) THEN ! high vegetation; maximum reduction to 40%
            z0_mod = ext_data%atm%z0_lcc(lc_class) * SQRT( MAX(0.16_wp, &
              MIN(1._wp,1.3333_wp*ext_data%atm%plcov_t(jc,jb,jt)) ))
          ELSE ! lower vegetation, maximum reduction to 70%
            z0_mod = ext_data%atm%z0_lcc(lc_class) * SQRT( MAX(0.5_wp, &
              MIN(1._wp,1.3333_wp*ext_data%atm%plcov_t(jc,jb,jt)) ))
          ENDIF
          ! ensure that z0 does not fall below the minimum allowed value
          z0_mod = MAX(z0_mod,ext_data%atm%z0_lcc_min(lc_class))
          ! Modify roughness length depending on snow cover
          prm_diag%gz0_t(jc,jb,jt) = grav *( (1._wp-lnd_diag%snowfrac_t(jc,jb,jt)**2)*z0_mod + &
            lnd_diag%snowfrac_t(jc,jb,jt)**2*ext_data%atm%z0_lcc_min(lc_class) )
          ! Set gz0 on empty snow-free tiles because this is used in the snow-cover fraction diagnosis
          IF (jt > ntiles_lnd .AND. lnd_diag%snowfrac_lc_t(jc,jb,jt) > 0.999_wp) THEN
            prm_diag%gz0_t(jc,jb,jt-ntiles_lnd) = grav*z0_mod
          ENDIF
        ENDDO
        !$acc end parallel
        IF (atm_phy_nwp_config(jg)%itype_z0 == 3) THEN ! Add SSO contribution to tile-specific roughness length
          !$acc parallel default(present) if(lzacc)
          !$acc loop gang vector private(jc)
!$NEC ivdep
          DO ic = 1, ext_data%atm%gp_count_t(jb,jt)
            jc = ext_data%atm%idx_lst_t(ic,jb,jt)
            prm_diag%gz0_t(jc,jb,jt) = prm_diag%gz0_t(jc,jb,jt) + grav * &
              MIN(fact_z0rough*ext_data%atm%sso_stdh_raw(jc,jb)**2,7.5_wp)
          ENDDO
          !$acc end parallel
        ENDIF
      ENDDO
    ELSE ! uniform tile-averaged roughness length if SSO contribution is to be included
      !$acc parallel default(present) if(lzacc)
      !$acc loop gang vector collapse(2)
      DO jt = 1, ntiles_total + ntiles_water
        DO jc = i_startidx, i_endidx
          prm_diag%gz0_t(jc,jb,jt) = prm_diag%gz0(jc,jb)
        ENDDO
      ENDDO
      !$acc end parallel
    ENDIF

    IF (atm_phy_nwp_config(jg)%inwp_sso > 0) THEN
      !$acc parallel default(present) if(lzacc)
      !$acc loop gang vector
      DO jc = i_startidx, i_endidx
        IF (prm_diag%ktop_envel(jc,jb) < nlev) THEN
          jk_gust(jc) = prm_diag%ktop_envel(jc,jb) - 1
        ELSE
          jk_gust(jc) = nlev
        ENDIF
      ENDDO
      !$acc end parallel
    ELSE
      !$acc kernels default(present) if(lzacc)
      jk_gust(:) = nlev
      !$acc end kernels
    ENDIF

    SELECT CASE (atm_phy_nwp_config(jg)%inwp_turb)

    CASE(icosmo)

!-------------------------------------------------------------------------
!< COSMO turbulence scheme by M. Raschendorfer
!-------------------------------------------------------------------------
 
      ierrstat = 0

      ! note that TKE must be converted to the turbulence velocity scale SQRT(2*TKE)
      ! for turbdiff
      ! INPUT to turbtran is timestep new
      !$acc parallel default(present) if(lzacc)
      !$acc loop gang vector collapse(2)
      DO jk = 1, 3
        DO jc = i_startidx, i_endidx
          z_tvs(jc,jk,1) = SQRT(2._wp * p_prog_rcf%tke(jc,nlev-2+jk,jb))
        ENDDO
      ENDDO
      !$acc end parallel
      ! First call of turbtran for all grid points (water points with > 50% water
      ! fraction and tile 1 of the land points)
      IF (ntiles_total == 1) THEN ! tile approach not used; use tile-averaged fields from extpar

        !should be dependent on location in future!
        !$acc kernels default(present) if(lzacc)
        l_hori(i_startidx:i_endidx)=phy_params(jg)%mean_charlen
        !$acc end kernels

        nzprv  = 1
        nlevcm = 3

        ! turbtran
        CALL turbtran (               & ! only surface-layer turbulence
          &  iini=0,                  & !
          &  ltkeinp=.FALSE.,         & !
          &  lgz0inp=.FALSE.,         & !
          &  lstfnct=.TRUE. ,         & ! with stability function
          &  lsrflux=.TRUE. ,         & !
          &  lnsfdia=.TRUE. ,         & ! including near-surface diagnostics
          &  lrunscm=.FALSE.,         & ! no single column run
          &  dt_tke=tcall_turb_jg,                                                     & !in
          &  nprv=nzprv, ntur=1, ntim=1,                                               & !in
          &  nvec=nproma, ke=2,    ke1=3,      kcm=nlevcm, iblock=jb,                  & !in
          &  ivstart=i_startidx, ivend=i_endidx,                                       & !in
          &  l_hori=l_hori,                                                            & !in
          &  hhl=p_metrics%z_ifc(:,nlev-1:nlevp1,jb),                                  & !in
          &  fr_land=ext_data%atm%fr_land(:,jb),                                       & !in
          &  depth_lk=ext_data%atm%depth_lk(:,jb),                                     & !in
          &  h_ice=wtr_prog_new%h_ice(:,jb),                                           & !in
          &  sai=ext_data%atm%sai_t(:,jb,1),                                           & !in
          &  gz0=prm_diag%gz0_t(:,jb,1),                                               & !inout
          &  t_g=lnd_prog_new%t_g(:,jb),                                               & !in
          &  qv_s=lnd_diag%qv_s(:,jb),                                                 & !in
          &  ps=p_diag%pres_sfc(:,jb),                                                 & !in
          &  u=p_diag%u(:,nlev-1:nlev,jb),                                             & !in
          &  v=p_diag%v(:,nlev-1:nlev,jb),                                             & !in
          &  t=p_diag%temp(:,nlev-1:nlev,jb),                                          & !in
          &  prs=p_diag%pres(:,nlev-1:nlev,jb),                                        & !in
          &  qv=p_prog_rcf%tracer(:,nlev-1:nlev,jb,iqv),                               & !in
          &  qc=p_prog_rcf%tracer(:,nlev-1:nlev,jb,iqc),                               & !in
          &  epr=p_prog%exner(:,nlev-1:nlev,jb),                                       & !in
          &  tcm=prm_diag%tcm_t(:,jb,1),                                               & !out
          &  tch=prm_diag%tch_t(:,jb,1),                                               & !out
          &  tvm=prm_diag%tvm_t(:,jb,1),                                               & !inout
          &  tvh=prm_diag%tvh_t(:,jb,1),                                               & !inout
          &  tfm=prm_diag%tfm(:,jb),                                                   & !inout
          &  tfh=prm_diag%tfh(:,jb),                                                   & !inout
          &  tfv=prm_diag%tfv_t(:,jb,1),                                               & !inout
          &  tkr=prm_diag%tkr_t(:,jb,1),                                               & !inout
          &  tke=z_tvs(:,:,:),                                                         & !inout
          &  tkvm=prm_diag%tkvm(:,nlev-1:nlevp1,jb),                                   & !inout
          &  tkvh=prm_diag%tkvh(:,nlev-1:nlevp1,jb),                                   & !inout
          &  rcld=prm_diag%rcld(:,nlev-1:nlevp1,jb),                                   & !inout
          &  t_2m=prm_diag%t_2m(:,jb),                                                 & !inout
          &  qv_2m=prm_diag%qv_2m(:,jb),                                               & !out
          &  td_2m=prm_diag%td_2m(:,jb),                                               & !out
          &  rh_2m=prm_diag%rh_2m(:,jb),                                               & !out
          &  u_10m=prm_diag%u_10m_t(:,jb,1),                                           & !out
          &  v_10m=prm_diag%v_10m_t(:,jb,1),                                           & !out
          &  shfl_s=prm_diag%shfl_s_t(:,jb,1),                                         & !out
          &  qvfl_s=prm_diag%qhfl_s_t(:,jb,1),                                         & !out
          &  umfl_s=prm_diag%umfl_s_t(:,jb,1),                                         & !out
          &  vmfl_s=prm_diag%vmfl_s_t(:,jb,1),                                         & !out
          &  ierrstat=ierrstat, yerrormsg=errormsg, yroutine=eroutine,                 & !inout
          &  lacc=lzacc                                                                 ) !in

        !$acc kernels default(present) if(lzacc)
        prm_diag%lhfl_s_t(i_startidx:i_endidx,jb,1) = &
          &  prm_diag%qhfl_s_t(i_startidx:i_endidx,jb,1) * lh_v
        !$acc end kernels

        ! fix latent heat flux over seaice
        !$acc parallel default(present) if(lzacc)
        !$acc loop gang vector
        DO jc = i_startidx, i_endidx
          IF (wtr_prog_new%h_ice(jc,jb) > 0._wp ) THEN
            prm_diag%lhfl_s_t(jc,jb,1) = (lh_s/lh_v) * prm_diag%lhfl_s_t(jc,jb,1)
          ENDIF
        ENDDO
        !$acc end parallel

        ! copy
        !$acc kernels default(present) if(lzacc)
        prm_diag%gz0(i_startidx:i_endidx,jb)   = prm_diag%gz0_t(i_startidx:i_endidx,jb,1)
        prm_diag%tcm(i_startidx:i_endidx,jb)   = prm_diag%tcm_t(i_startidx:i_endidx,jb,1)
        prm_diag%tch(i_startidx:i_endidx,jb)   = prm_diag%tch_t(i_startidx:i_endidx,jb,1)
        prm_diag%tfv(i_startidx:i_endidx,jb)   = prm_diag%tfv_t(i_startidx:i_endidx,jb,1)
        prm_diag%tvm(i_startidx:i_endidx,jb)   = prm_diag%tvm_t(i_startidx:i_endidx,jb,1)
        prm_diag%tvh(i_startidx:i_endidx,jb)   = prm_diag%tvh_t(i_startidx:i_endidx,jb,1)
        prm_diag%tkr(i_startidx:i_endidx,jb)   = prm_diag%tkr_t(i_startidx:i_endidx,jb,1)
        prm_diag%u_10m(i_startidx:i_endidx,jb) = prm_diag%u_10m_t(i_startidx:i_endidx,jb,1)
        prm_diag%v_10m(i_startidx:i_endidx,jb) = prm_diag%v_10m_t(i_startidx:i_endidx,jb,1)

        ! instantaneous max/min 2m temperature over tiles (trivial operation for 1 tile)
        prm_diag%t_tilemax_inst_2m(i_startidx:i_endidx,jb) = prm_diag%t_2m(i_startidx:i_endidx,jb)
        prm_diag%t_tilemin_inst_2m(i_startidx:i_endidx,jb) = prm_diag%t_2m(i_startidx:i_endidx,jb)
        prm_diag%tmax_2m(i_startidx:i_endidx,jb) = MAX(prm_diag%t_2m(i_startidx:i_endidx,jb), &
          &                                        prm_diag%tmax_2m(i_startidx:i_endidx,jb) )
        prm_diag%tmin_2m(i_startidx:i_endidx,jb) = MIN(prm_diag%t_2m(i_startidx:i_endidx,jb), &
          &                                        prm_diag%tmin_2m(i_startidx:i_endidx,jb) )
        !$acc end kernels

      ELSE ! tile approach used

        ! preset variables for land tile indices
        !$acc kernels default(present) if(lzacc)
        fr_land_t(:)  = 1._wp
        depth_lk_t(:) = 0._wp
        h_ice_t(:)    = 0._wp
        !$acc end kernels


        ! Loop over land tile points, sea, lake points and seaice points
        ! Each tile has a separate index list
        DO  jt = 1, ntiles_total + ntiles_water

          IF (jt <= ntiles_total) THEN ! land tile points
            i_count =  ext_data%atm%gp_count_t(jb,jt)
            ilist   => ext_data%atm%idx_lst_t(:,jb,jt)
          ELSE IF (jt == ntiles_total + 1) THEN ! sea points (open water)
            i_count =  ext_data%atm%list_seawtr%ncount(jb)
            ilist   => ext_data%atm%list_seawtr%idx(:,jb)
            !$acc kernels default(present) if(lzacc)
            fr_land_t(:) = 0._wp
            !$acc end kernels
          ELSE IF (jt == ntiles_total + 2) THEN ! lake points
            i_count =  ext_data%atm%list_lake%ncount(jb)
            ilist   => ext_data%atm%list_lake%idx(:,jb)
            !$acc kernels default(present) if(lzacc)
            fr_land_t (:) = 0._wp
            depth_lk_t(:) = 1._wp
            !$acc end kernels
            IF (llake) THEN
              !$acc parallel default(present) if(lzacc)
              !$acc loop gang vector private (jc)
              DO ic= 1, i_count
                jc = ilist(ic)
                h_ice_t(ic) = MERGE(1._wp, 0._wp, wtr_prog_new%h_ice(jc,jb)>=h_Ice_min_flk)
              ENDDO
              !$acc end parallel
            ELSE
              !$acc parallel default(present) if(lzacc)
              !$acc loop gang vector private (jc)
              DO ic= 1, i_count
                jc = ilist(ic)
                h_ice_t(ic) = MERGE(1._wp, 0._wp, lnd_prog_new%t_g_t(jc,jb,isub_lake)< tpl_T_f)
              ENDDO
              !$acc end parallel
            ENDIF
          ELSE IF (jt == ntiles_total + 3) THEN ! seaice points
            ! Note that if the sea-ice scheme is not used (lseaice=.FALSE.), list_seaice%ncount=0.
            i_count =  ext_data%atm%list_seaice%ncount(jb)
            ilist   => ext_data%atm%list_seaice%idx(:,jb)
            !$acc kernels default(present) if(lzacc)
            fr_land_t (:) = 0._wp
            depth_lk_t(:) = 0._wp
            h_ice_t   (:) = 1._wp  ! Only needed for checking whether ice is present or not
            !$acc end kernels
          ELSE
            ! Paranoia
            CALL finish( TRIM(routine),'wrong value of ntiles_total + ntiles_water')
          ENDIF

          IF (i_count == 0) CYCLE ! skip loop if the index list for the given tile is empty

          ! Copy input fields to the local re-indexed variables
          ! It remains to be determined which of the model levels are actually needed for non-init calls
          !
          !MR: Hauptflaechengroessen nur fuer level nlev
!$NEC ivdep
          !$acc kernels default(present) if(lzacc)
          !$acc loop
          DO ic = 1, i_count
            jc = ilist(ic)

            z_ifc_t (ic,1:3)    = p_metrics%z_ifc   (jc,nlev-1:nlevp1,jb)
            u_t    (ic,1:2)     = p_diag%u          (jc,nlev-1:nlev  ,jb)
            v_t    (ic,1:2)     = p_diag%v          (jc,nlev-1:nlev  ,jb)
            temp_t (ic,1:2)     = p_diag%temp       (jc,nlev-1:nlev  ,jb)
            pres_t (ic,1:2)     = p_diag%pres       (jc,nlev-1:nlev  ,jb)
            epr_t  (ic,1:2)     = p_prog%exner      (jc,nlev-1:nlev  ,jb)
            qv_t   (ic,1:2)     = p_prog_rcf%tracer (jc,nlev-1:nlev  ,jb,iqv)
            qc_t   (ic,1:2)     = p_prog_rcf%tracer (jc,nlev-1:nlev  ,jb,iqc)
            pres_sfc_t(ic)      = p_diag%pres_sfc   (jc,jb)

            gz0_t  (ic,jt)      = prm_diag%gz0_t    (jc,jb,jt)
            t_g_t  (ic,jt)      = lnd_prog_new%t_g_t(jc,jb,jt)
            qv_s_t (ic,jt)      = lnd_diag%qv_s_t   (jc,jb,jt)
            sai_t  (ic,jt)      = ext_data%atm%sai_t(jc,jb,jt)
!MR: rcld: benoetigt nur fuer level nlev (als Nebenflaechenvariable)
            rcld_t (ic,1:3,jt)  = prm_diag%rcld     (jc,nlev-1:nlevp1,jb)
            tvs_t  (ic,1:2,1,jt)= z_tvs             (jc,1:2,1)
            tvs_t  (ic,3,1,jt)  = prm_diag%tvs_s_t  (jc,jb,jt)     ! tile-specific for lowest level
            tkvm_t (ic,1:2,jt)  = prm_diag%tkvm     (jc,nlev-1:nlev,jb)
            tkvm_t (ic,3,jt)    = prm_diag%tkvm_s_t (jc,jb,jt)     ! tile-specific for lowest level
            tkvh_t (ic,1:2,jt)  = prm_diag%tkvh     (jc,nlev-1:nlev,jb)
            tkvh_t (ic,3,jt)    = prm_diag%tkvh_s_t (jc,jb,jt)     ! tile-specific for lowest level
            tkr_t  (ic,jt)      = prm_diag%tkr_t    (jc,jb,jt)

            !should be dependent on location in future!
            l_hori(ic)=phy_params(jg)%mean_charlen
            
          ENDDO
          !$acc end kernels

          nlevcm = 3
          nzprv  = 1

          ! turbtran
          CALL turbtran (               & ! only surface-layer turbulence
            &  iini=0,                  & !
            &  ltkeinp=.FALSE.,         & !
            &  lgz0inp=.FALSE.,         & !
            &  lstfnct=.TRUE. ,         & ! with stability function
            &  lsrflux=.TRUE. ,         & !
            &  lnsfdia=.TRUE. ,         & ! including near-surface diagnostics
            &  lrunscm=.FALSE.,         & ! no single column run
!
            &  dt_tke=tcall_turb_jg,                                                     & !in
            &  nprv=nzprv, ntur=1, ntim=1,                                               & !in
            &  nvec=nproma, ke=2, ke1=3, kcm=nlevcm, iblock=jb,             & !in
            &  ivstart=1, ivend=i_count,                                    & !in
!
            &  l_hori=l_hori,                                               & !in
            &  hhl=z_ifc_t(:,:),                                            & !in
            &  fr_land=fr_land_t(:),                                        & !in
            &  depth_lk=depth_lk_t(:),                                      & !in
            &  h_ice=h_ice_t(:),                                            & !in
            &  sai=sai_t(:,jt),                                             & !in
            &  gz0=gz0_t(:,jt),                                             & !inout
            &  t_g=t_g_t(:,jt),                                             & !in
            &  qv_s=qv_s_t(:,jt),                                           & !in
            &  ps=pres_sfc_t(:),                                            & !in
            &  u=u_t(:,:),                                                  & !in
            &  v=v_t(:,:),                                                  & !in
            &  t=temp_t(:,:),                                               & !in
            &  prs=pres_t(:,:),                                             & !in
            &  qv=qv_t(:,:),                                                & !in
            &  qc=qc_t(:,:),                                                & !in
            &  epr=epr_t(:,:),                                              & !in
            &  tcm=tcm_t(:,jt),                                             & !out
            &  tch=tch_t(:,jt),                                             & !out
            &  tvm=tvm_t(:,jt),                                             & !inout
            &  tvh=tvh_t(:,jt),                                             & !inout
            &  tkr=tkr_t(:,jt),                                             & !inout
            &  tfm=tfm_t(:,jt),                                             & !inout
            &  tfh=tfh_t(:,jt),                                             & !inout
            &  tfv=tfv_t(:,jt),                                             & !inout
            &  tke=tvs_t(:,:,:,jt),                                         & !inout
            &  tkvm=tkvm_t(:,:,jt),                                         & !inout
            &  tkvh=tkvh_t(:,:,jt),                                         & !inout
            &  rcld=rcld_t(:,:,jt),                                         & !inout
            &  t_2m=t_2m_t(:,jt),                                           & !out
            &  qv_2m=qv_2m_t(:,jt),                                         & !out
            &  td_2m=td_2m_t(:,jt),                                         & !out
            &  rh_2m=rh_2m_t(:,jt),                                         & !out
            &  u_10m=u_10m_t(:,jt),                                         & !out
            &  v_10m=v_10m_t(:,jt),                                         & !out
            &  shfl_s=shfl_s_t(:,jt),                                       & !out
            &  qvfl_s=qhfl_s_t(:,jt),                                       & !out
            &  umfl_s=umfl_s_t(:,jt),                                       & !out
            &  vmfl_s=vmfl_s_t(:,jt),                                       & !out
            &  ierrstat=ierrstat, yerrormsg=errormsg, yroutine=eroutine,    & !inout
            &  lacc=lzacc                                                   ) !in

          ! Decision as to "ice" vs. "no ice" is made on the basis of h_ice_t(:).
          !$acc parallel default(present) if(lzacc)
          !$acc loop gang vector
          DO ic= 1, i_count
            lhfl_s_t(ic,jt) = MERGE(qhfl_s_t(ic,jt)*lh_v, qhfl_s_t(ic,jt)*lh_s,  &
              &                     h_ice_t(ic)<h_Ice_min_flk)
          ENDDO
          !$acc end parallel

        ENDDO ! loop over tiles


        ! Aggregate tile-based output fields of turbtran over tiles
        ! i) initialize fields to zero before starting the summation
        !$acc kernels default(present) if(lzacc)
        prm_diag%gz0   (i_startidx:i_endidx,jb) = 0._wp
        prm_diag%tcm   (i_startidx:i_endidx,jb) = 0._wp
        prm_diag%tch   (i_startidx:i_endidx,jb) = 0._wp
        prm_diag%tfm   (i_startidx:i_endidx,jb) = 0._wp
        prm_diag%tfh   (i_startidx:i_endidx,jb) = 0._wp
        prm_diag%tfv   (i_startidx:i_endidx,jb) = 0._wp
        prm_diag%tvm   (i_startidx:i_endidx,jb) = 0._wp
        prm_diag%tvh   (i_startidx:i_endidx,jb) = 0._wp
        prm_diag%tkr   (i_startidx:i_endidx,jb) = 0._wp
        prm_diag%t_2m  (i_startidx:i_endidx,jb) = 0._wp
        prm_diag%qv_2m (i_startidx:i_endidx,jb) = 0._wp
        prm_diag%td_2m (i_startidx:i_endidx,jb) = 0._wp
        prm_diag%rh_2m (i_startidx:i_endidx,jb) = 0._wp
        prm_diag%u_10m (i_startidx:i_endidx,jb) = 0._wp
        prm_diag%v_10m (i_startidx:i_endidx,jb) = 0._wp

        prm_diag%t_2m_land (i_startidx:i_endidx,jb) = 0._wp
        prm_diag%td_2m_land(i_startidx:i_endidx,jb) = 0._wp
        prm_diag%rh_2m_land(i_startidx:i_endidx,jb) = 0._wp

        z_tvs        (i_startidx:i_endidx,3, 1)      = 0._wp
        prm_diag%tkvm(i_startidx:i_endidx,nlevp1,jb) = 0._wp
        prm_diag%tkvh(i_startidx:i_endidx,nlevp1,jb) = 0._wp
        prm_diag%rcld(i_startidx:i_endidx,nlevp1,jb) = 0._wp

        ! re-initialize dyn_gust for proper maximum computation over tiles
        ! note that dyn_gust is an instantaneous field
        prm_diag%dyn_gust(i_startidx:i_endidx,jb) = 0._wp

        prm_diag%t_tilemax_inst_2m(i_startidx:i_endidx,jb) = -999._wp
        prm_diag%t_tilemin_inst_2m(i_startidx:i_endidx,jb) = 999._wp
        !$acc end kernels

         ! ii) loop over index lists
        DO  jt = 1, ntiles_total + ntiles_water

          IF (jt <= ntiles_total) THEN ! land tile points
            i_count = ext_data%atm%gp_count_t(jb,jt)
            ilist => ext_data%atm%idx_lst_t(:,jb,jt)
          ELSE IF (jt == ntiles_total + 1) THEN ! sea points (seaice points excluded)
            i_count = ext_data%atm%list_seawtr%ncount(jb)
            ilist => ext_data%atm%list_seawtr%idx(:,jb)
          ELSE IF (jt == ntiles_total + 2) THEN ! lake points
            i_count = ext_data%atm%list_lake%ncount(jb)
            ilist => ext_data%atm%list_lake%idx(:,jb)
          ELSE ! IF (jt == ntiles_total + 3) THEN ! seaice points
            i_count = ext_data%atm%list_seaice%ncount(jb)
            ilist => ext_data%atm%list_seaice%idx(:,jb)
          ENDIF

          IF (i_count == 0) CYCLE ! skip loop if the index list for the given tile is empty

!$NEC ivdep
          !$acc parallel default(present) if(lzacc)
          !$acc loop gang vector private(jc,area_frac)
          DO ic = 1, i_count
            jc = ilist(ic)
            area_frac = ext_data%atm%frac_t(jc,jb,jt)

            ! Aggregate
            prm_diag%gz0(jc,jb) = prm_diag%gz0(jc,jb)+gz0_t(ic,jt) * area_frac
            prm_diag%tcm(jc,jb) = prm_diag%tcm(jc,jb)+tcm_t(ic,jt) * area_frac
            prm_diag%tch(jc,jb) = prm_diag%tch(jc,jb)+tch_t(ic,jt) * area_frac
            prm_diag%tfm(jc,jb) = prm_diag%tfm(jc,jb)+tfm_t(ic,jt) * area_frac
            prm_diag%tfh(jc,jb) = prm_diag%tfh(jc,jb)+tfh_t(ic,jt) * area_frac
            prm_diag%tfv(jc,jb) = prm_diag%tfv(jc,jb)+tfv_t(ic,jt) * area_frac
            prm_diag%tvm(jc,jb) = prm_diag%tvm(jc,jb)+tvm_t(ic,jt) * area_frac
            prm_diag%tvh(jc,jb) = prm_diag%tvh(jc,jb)+tvh_t(ic,jt) * area_frac
            prm_diag%tkr(jc,jb) = prm_diag%tkr(jc,jb)+tkr_t(ic,jt) * area_frac !not necessary

            z_tvs(jc,3,1)               = z_tvs(jc,3,1)+tvs_t(ic,3,1,jt)              * area_frac
            prm_diag%tkvm(jc,nlevp1,jb) = prm_diag%tkvm(jc,nlevp1,jb)+tkvm_t(ic,3,jt) * area_frac
            prm_diag%tkvh(jc,nlevp1,jb) = prm_diag%tkvh(jc,nlevp1,jb)+tkvh_t(ic,3,jt) * area_frac
            prm_diag%rcld(jc,nlevp1,jb) = prm_diag%rcld(jc,nlevp1,jb)+rcld_t(ic,3,jt) * area_frac

            prm_diag%t_2m  (jc,jb) = prm_diag%t_2m(jc,jb)  + t_2m_t(ic,jt)  * area_frac
            prm_diag%qv_2m (jc,jb) = prm_diag%qv_2m(jc,jb) + qv_2m_t(ic,jt) * area_frac
            prm_diag%td_2m (jc,jb) = prm_diag%td_2m(jc,jb) + td_2m_t(ic,jt) * area_frac
            prm_diag%rh_2m (jc,jb) = prm_diag%rh_2m(jc,jb) + rh_2m_t(ic,jt) * area_frac
            prm_diag%u_10m (jc,jb) = prm_diag%u_10m(jc,jb) + u_10m_t(ic,jt) * area_frac
            prm_diag%v_10m (jc,jb) = prm_diag%v_10m(jc,jb) + v_10m_t(ic,jt) * area_frac

            prm_diag%t_tilemax_inst_2m(jc,jb) = MAX(t_2m_t(ic,jt),prm_diag%t_tilemax_inst_2m(jc,jb))
            prm_diag%t_tilemin_inst_2m(jc,jb) = MIN(t_2m_t(ic,jt),prm_diag%t_tilemin_inst_2m(jc,jb))


            ! Store
            prm_diag%shfl_s_t(jc,jb,jt) = shfl_s_t(ic,jt)
            prm_diag%lhfl_s_t(jc,jb,jt) = lhfl_s_t(ic,jt)
            prm_diag%qhfl_s_t(jc,jb,jt) = qhfl_s_t(ic,jt)
            prm_diag%umfl_s_t(jc,jb,jt) = umfl_s_t(ic,jt)
            prm_diag%vmfl_s_t(jc,jb,jt) = vmfl_s_t(ic,jt)
            prm_diag%u_10m_t (jc,jb,jt) = u_10m_t(ic,jt)  ! needed by TERRA and turbtran
            prm_diag%v_10m_t (jc,jb,jt) = v_10m_t(ic,jt)  ! needed by TERRA and turbtran
            prm_diag%tch_t   (jc,jb,jt) = tch_t(ic,jt)    ! needed by TERRA
            prm_diag%tcm_t   (jc,jb,jt) = tcm_t(ic,jt)    ! needed by TERRA
            prm_diag%tfv_t   (jc,jb,jt) = tfv_t(ic,jt)    ! needed by TERRA

            prm_diag%tvm_t   (jc,jb,jt) = tvm_t(ic,jt)    ! to be used by TERRA instead of 'tcm'
            prm_diag%tvh_t   (jc,jb,jt) = tvh_t(ic,jt)    ! to be used by TERRA instead of 'tch'
            prm_diag%tkr_t   (jc,jb,jt) = tkr_t(ic,jt)    ! input for turbtran
            prm_diag%gz0_t   (jc,jb,jt) = gz0_t(ic,jt)    ! input for turbtran at n+1
                                                          ! over land
            prm_diag%tvs_s_t (jc,jb,jt) = tvs_t(ic,3,1,jt)! needed as input for turbtran
            prm_diag%tkvm_s_t(jc,jb,jt) = tkvm_t(ic,3,jt) ! needed as input for turbtran
            prm_diag%tkvh_s_t(jc,jb,jt) = tkvh_t(ic,3,jt) ! needed as input for turbtran

          ENDDO
          !$acc end parallel

          ! averages over land fraction of mixed land-water points
          IF (jt <= ntiles_total) THEN
            !$acc parallel default(present) if(lzacc)
            !$acc loop gang vector private(jc,area_frac)
!$NEC ivdep
            DO ic = 1, i_count
              jc = ilist(ic)
              area_frac = ext_data%atm%frac_t(jc,jb,jt)/ext_data%atm%fr_land(jc,jb)
              prm_diag%t_2m_land  (jc,jb) = prm_diag%t_2m_land(jc,jb)  + t_2m_t(ic,jt)  * area_frac
              prm_diag%td_2m_land (jc,jb) = prm_diag%td_2m_land(jc,jb) + td_2m_t(ic,jt) * area_frac
              prm_diag%rh_2m_land (jc,jb) = prm_diag%rh_2m_land(jc,jb) + rh_2m_t(ic,jt) * area_frac
            ENDDO
            !$acc end parallel
          ELSE
            !$acc parallel default(present) if(lzacc)
            !$acc loop gang vector private(jc,area_frac)
!$NEC ivdep
            DO ic = 1, i_count
              jc = ilist(ic)
              IF (ext_data%atm%fr_land(jc,jb) == 0._wp) THEN
                area_frac = ext_data%atm%frac_t(jc,jb,jt)
                prm_diag%t_2m_land  (jc,jb) = prm_diag%t_2m_land(jc,jb)  + t_2m_t(ic,jt)  * area_frac
                prm_diag%td_2m_land (jc,jb) = prm_diag%td_2m_land(jc,jb) + td_2m_t(ic,jt) * area_frac
                prm_diag%rh_2m_land (jc,jb) = prm_diag%rh_2m_land(jc,jb) + rh_2m_t(ic,jt) * area_frac
              ENDIF
            ENDDO
            !$acc end parallel
          ENDIF

        ENDDO  ! jt

      ENDIF ! tiles / no tiles


      ! Dynamic gusts are diagnosed from averaged values in order to avoid artifacts along coastlines
      !$acc parallel default(present) if(lzacc)
      !$acc loop gang vector
      DO jc = i_startidx, i_endidx
        prm_diag%dyn_gust(jc,jb) =  nwp_dyn_gust (prm_diag%u_10m(jc,jb),      &
          &                                       prm_diag%v_10m(jc,jb),      &
          &                                       prm_diag%tcm  (jc,jb),      &
          &                                       p_diag%u      (jc,nlev,jb), &
          &                                       p_diag%v      (jc,nlev,jb), &
          &                                       p_diag%u(jc,jk_gust(jc),jb),&
          &                                       p_diag%v(jc,jk_gust(jc),jb),&
          &                                  p_metrics%mask_mtnpoints_g(jc,jb))
      ENDDO
      !$acc end parallel

      ! transform updated turbulent velocity scale back to TKE
      !$acc kernels default(present) if(lzacc)
      p_prog_rcf%tke(i_startidx:i_endidx,nlevp1,jb)= 0.5_wp*(z_tvs(i_startidx:i_endidx,3,1))**2
      !$acc end kernels



      ! Re-compute TKE at lowest main levels. Note that slight temporal inconsistencies are
      ! ignored at this point.
      IF (advection_config(jg)%iadv_tke > 0) THEN
        !$acc parallel default(present) if(lzacc)
        !$acc loop gang vector
        DO jc=i_startidx, i_endidx
          p_prog_rcf%tracer(jc,nlev,jb,iqtke) = 0.5_wp* ( z_tvs(jc,2,1) + z_tvs(jc,3,1) )
        ENDDO
        !$acc end parallel
      ENDIF

    CASE(igme,ismag,iprog)

#ifdef _OPENACC
      CALL warning('GPU:mo_nwp_turbtrans_interface:nwp_turbtrans', 'igme/ismag/iprog unsupported. Only cosmo turbulence is supported on GPU!')
#endif
!-------------------------------------------------------------------------
!> GME turbulence scheme
!-------------------------------------------------------------------------

      ! turbulent diffusion coefficients at the surface
#ifndef COUP_OAS_ICON
      CALL parturs( zsurf=p_metrics%z_ifc(:,nlevp1,jb), z1=p_metrics%z_mc(:,nlev,jb),   & !in
        &           u1=p_diag%u(:,nlev,jb), v1=p_diag%v(:,nlev,jb),                     & !in
        &           t1=p_diag%temp(:,nlev,jb), qv1=p_prog_rcf%tracer(:,nlev,jb,iqv),    & !in
        &           t_g=lnd_prog_new%t_g(:,jb), qv_s=lnd_diag%qv_s(:,jb),               & !in
        &           ps=p_diag%pres_ifc(:,nlevp1,jb),                                    & !in
        &           fr_land=ext_data%atm%fr_land(:,jb), h_ice=wtr_prog_new%h_ice(:,jb), & !in
        &           ie=nproma, i_startidx=i_startidx, i_endidx=i_endidx,                & !in
        &           tcm=prm_diag%tcm(:,jb), tch=prm_diag%tch(:,jb),                     & !out
        &           gz0=prm_diag%gz0(:,jb),       shfl_s=prm_diag%shfl_s(:,jb),         & !inout, out
        &           lhfl_s=prm_diag%lhfl_s(:,jb), qhfl_s=prm_diag%qhfl_s(:,jb),         & !out, out
        &           umfl_s=prm_diag%umfl_s_t(:,jb,1), vmfl_s=prm_diag%vmfl_s_t(:,jb,1))   !out, out
#endif


      !DR inside "nearsfc", lhfl_s is converted to qhfl_s via
      !DR qhfl_s = lhfl_s/lh_v. This is incorrect over snow and ice.
      !DR Shouldn't we simply pass qhfl_s ?
      !
      ! diagnose 2 m temperature, humidity, 10 m wind
      CALL nearsfc( t=p_diag%temp(:,:,jb), qv=p_prog_rcf%tracer(:,:,jb,iqv),            & !in
        &           u=p_diag%u(:,:,jb),    v=p_diag%v(:,:,jb),                          & !in
        &           zf=p_metrics%z_mc(:,:,jb), ps=p_diag%pres_ifc(:,nlevp1,jb),         & !in
        &           t_g=lnd_prog_new%t_g(:,jb),                                         & !in
        &           tcm=prm_diag%tcm(:,jb), tch=prm_diag%tch(:,jb),                     & !in
        &           gz0=prm_diag%gz0(:,jb),                                             & !in
        &           shfl_s=prm_diag%shfl_s(:,jb), lhfl_s=prm_diag%lhfl_s(:,jb),         & !in
        &           umfl_s=prm_diag%umfl_s_t(:,jb,1), vmfl_s=prm_diag%vmfl_s_t(:,jb,1), & !in
        &           zsurf=p_metrics%z_ifc(:,nlevp1,jb),                                 & !in
        &           fr_land=ext_data%atm%fr_land(:,jb), pf1=p_diag%pres(:,nlev,jb),     & !in
        &           qv_s=lnd_diag%qv_s(:,jb), ie=nproma, ke=nlev,                       & !in
        &           i_startidx=i_startidx, i_endidx=i_endidx,                           & !in
        &           t_2m=prm_diag%t_2m(:,jb), qv_2m=prm_diag%qv_2m(:,jb),               & !out
        &           td_2m=prm_diag%td_2m(:,jb), rh_2m=prm_diag%rh_2m(:,jb),             & !out
        &           u_10m=prm_diag%u_10m(:,jb), v_10m=prm_diag%v_10m(:,jb)              ) !out


      ! dynamic gusts
      DO jc = i_startidx, i_endidx
        prm_diag%dyn_gust(jc,jb) = nwp_dyn_gust(prm_diag%u_10m(jc,jb), prm_diag%v_10m(jc,jb), &
          &  prm_diag%tcm(jc,jb), p_diag%u(jc,nlev,jb), p_diag%v(jc,nlev,jb),                 &
          &  p_diag%u(jc,jk_gust(jc),jb), p_diag%v(jc,jk_gust(jc),jb),p_metrics%mask_mtnpoints_g(jc,jb) )
      ENDDO

      ! instantaneous max/min 2m temperature over tiles (trivial operation for 1 tile)
      prm_diag%t_tilemax_inst_2m(i_startidx:i_endidx,jb) = prm_diag%t_2m(i_startidx:i_endidx,jb)
      prm_diag%t_tilemin_inst_2m(i_startidx:i_endidx,jb) = prm_diag%t_2m(i_startidx:i_endidx,jb)


      DO jt = 1, ntiles_total+ntiles_water
        DO jc = i_startidx, i_endidx

          ! Copy transfer coefficients to tile-based variables, which are used in TERRA
          prm_diag%tcm_t(jc,jb,jt) = prm_diag%tcm(jc,jb)
          prm_diag%tch_t(jc,jb,jt) = prm_diag%tch(jc,jb)
          ! the GME turbulence scheme does not have tfv. Set tfv=1
          prm_diag%tfv_t(jc,jb,jt) = 1._wp    !   prm_diag%tfv(jc,jb)

          ! Copy transfer u_10m/v_10m to tile-based variables, which are used in TERRA
          prm_diag%u_10m_t(jc,jb,jt) = prm_diag%u_10m(jc,jb)
          prm_diag%v_10m_t(jc,jb,jt) = prm_diag%v_10m(jc,jb)

          ! Copy sensible and latent heat fluxes to tile-based variables
          ! (needed by Flake, sea-ice model)
          prm_diag%shfl_s_t(jc,jb,jt) = prm_diag%shfl_s(jc,jb)
          prm_diag%lhfl_s_t(jc,jb,jt) = prm_diag%lhfl_s(jc,jb)
          prm_diag%qhfl_s_t(jc,jb,jt) = prm_diag%qhfl_s(jc,jb)
        ENDDO
      ENDDO


    CASE(iedmf)

#ifdef _OPENACC
      CALL warning('GPU:mo_nwp_turbtrans_interface:nwp_turbtrans', 'iedmf unsupported. Only cosmo turbulence is supported on GPU!')
#endif
!-------------------------------------------------------------------------
!> EDMF surface layer
!-------------------------------------------------------------------------

!---------------------------------------------------------------------------------------
! Notes on EDMF surface layer:
!   - z0 should be different between momentum and heat/moisture (see chapter 11 ECMWF docu)
!   - z0m for land points (not sea/seaice/lake) is calculated above: prm_diag%gz0_t(jc,jb,jt)
!     note: TURBTRAN calculates z0h internally in a physical way (Matthias)
!           GME Louis diagnoses z0h as z0h=min(z0m,10cm)
!   - surface fluxes calculated as over ocean, other fluxes will be recalculated by TERRA
!   - attention: SEAICE and LAKE need fluxes as input
!---------------------------------------------------------------------------------------

      DO jt = 1, ntiles_total+ntiles_water

        DO jc = i_startidx, i_endidx
          PGEOMLEV(jc)  = p_metrics%geopot_agl(jc,nlev,jb)   
          PCPTGZLEV(jc) = p_metrics%geopot_agl(jc,nlev,jb)       + p_diag%temp(jc,nlev,jb)      * cpd
          PCPTSTI(jc)   = p_metrics%geopot_agl_ifc(jc,nlevp1,jb) + lnd_prog_new%t_g_t(jc,jb,jt) * cpd
          PUCURR(jc)    = 0.0_wp
          PVCURR(jc)    = 0.0_wp
          ZZ0MTI(jc)    = prm_diag%gz0_t(jc,jb,jt) / grav
          ZZ0HTI(jc)    = prm_diag%gz0_t(jc,jb,jt) / grav            ! modified (/10) in vupdz0_tile
          ZZ0QTI(jc)    = prm_diag%gz0_t(jc,jb,jt) / grav            ! same
        ENDDO      
      
        CALL VUPDZ0_TILE(i_startidx, i_endidx, nproma, jt, 1, 1.5_wp, &    ! kstep=1 (no initialization)
         & p_diag%u(:,nlev,jb), p_diag%v(:,nlev,jb),                  &
         & p_diag%temp(:,nlev,jb), p_prog_rcf%tracer(:,nlev,jb,iqv),  &
         & p_diag%pres_ifc(:,nlevp1,jb), PGEOMLEV,                    &
         & prm_diag%umfl_s_t (:,jb,jt), prm_diag%vmfl_s_t(:,jb,jt),   &
         & prm_diag%shfl_s_t (:,jb,jt), prm_diag%qhfl_s_t(:,jb,jt),   &
         & wtr_prog_new%h_ice(:,jb),   lnd_prog_new%t_g_t(:,jb,jt),   &
         & PUCURR, PVCURR, ZZ0MTI, ZZ0HTI, ZZ0QTI, ZBUOMTI, ZZDLTI)

        CALL VEXCS(i_startidx, i_endidx, nproma, 3, 0, .true., 100.0_wp, 1.5_wp, &  ! tcall_turb_jgrvdifts
         & p_diag%u(:,nlev,jb), p_diag%v(:,nlev,jb), p_diag%temp(:,nlev,jb),     &
         & p_prog_rcf%tracer(:,nlev,jb,iqv), p_diag%pres_ifc(:,nlevp1,jb),       &
         & PGEOMLEV, PCPTGZLEV, PCPTSTI, lnd_diag%qv_s_t(:,jb,jt),               &
         & ZZ0MTI, ZZ0HTI, ZZ0QTI,                                               &
         & ZZDLTI, ZBUOMTI, PUCURR, PVCURR, ZCFMTI, PCFHTI, PCFQTI,              &
         & prm_diag%tkvm_s_t(:,jb,jt), prm_diag%tkvh_s_t(:,jb,jt),               &
         & prm_diag%tcm_t   (:,jb,jt), prm_diag%tch_t   (:,jb,jt))                  ! K used in TERRA

        DO jc = i_startidx, i_endidx
          rho_s     = p_diag%pres_ifc(jc,nlevp1,jb) / ( rd * lnd_prog_new%t_g_t(jc,jb,jt) )
  
          prm_diag%shfl_s_t(jc,jb,jt) = prm_diag%tkvh_s_t(jc,jb,jt) * rho_s              &
            & * ( cpd  * ( p_diag%temp(jc,nlev,jb)     - lnd_prog_new%t_g_t(jc,jb,jt)  ) &
                + grav * ( p_metrics%z_mc (jc,nlev,jb) - p_metrics%z_ifc(jc,nlevp1,jb) ) )
  
          prm_diag%qhfl_s_t(jc,jb,jt) = prm_diag%tkvh_s_t(jc,jb,jt) * rho_s              &
            & * ( p_prog_rcf%tracer(jc,nlev,jb,iqv) - lnd_diag%qv_s_t(jc,jb,jt) )
  
          IF ( lseaice .AND. ( wtr_prog_new%h_ice(jc,jb) > h_ice_min_flk ) ) THEN
            prm_diag%lhfl_s_t(jc,jb,jt) = lh_s*prm_diag%qhfl_s_t(jc,jb,jt)
          ELSE
            prm_diag%lhfl_s_t(jc,jb,jt) = lh_v*prm_diag%qhfl_s_t(jc,jb,jt)
          END IF
  
! --- momentum flux:
          prm_diag%umfl_s_t(jc,jb,jt) = prm_diag%tkvm_s_t(jc,jb,jt) * rho_s * p_diag%u(jc,nlev,jb)
          prm_diag%vmfl_s_t(jc,jb,jt) = prm_diag%tkvm_s_t(jc,jb,jt) * rho_s * p_diag%v(jc,nlev,jb)
 
          prm_diag%gz0_t(jc,jb,jt) = ZZ0MTI(jc) * grav    ! arbitrary saving z0m in gz0
        ENDDO

      ENDDO ! tiles 

      ! -- diagnostics taken from GME code

      ! diagnose 2 m temperature, humidity, 10 m wind
      CALL nearsfc( t=p_diag%temp(:,:,jb), qv=p_prog_rcf%tracer(:,:,jb,iqv),            & !in
        &           u=p_diag%u(:,:,jb),    v=p_diag%v(:,:,jb),                          & !in
        &           zf=p_metrics%z_mc(:,:,jb), ps=p_diag%pres_ifc(:,nlevp1,jb),         & !in
        &           t_g=lnd_prog_new%t_g(:,jb),                                         & !in
        &           tcm=prm_diag%tcm_t(:,jb,1), tch=prm_diag%tch_t(:,jb,1),             & !in
        &           gz0=prm_diag%gz0_t(:,jb,1),                                         & !in
        &           shfl_s=prm_diag%shfl_s_t(:,jb,1), lhfl_s=prm_diag%lhfl_s_t(:,jb,1), & !in
        &           umfl_s=prm_diag%umfl_s_t(:,jb,1), vmfl_s=prm_diag%vmfl_s_t(:,jb,1), & !in
        &           zsurf=p_metrics%z_ifc(:,nlevp1,jb),                                 & !in
        &           fr_land=ext_data%atm%fr_land(:,jb), pf1=p_diag%pres(:,nlev,jb),     & !in
        &           qv_s=lnd_diag%qv_s_t(:,jb,1), ie=nproma, ke=nlev,                   & !in
        &           i_startidx=i_startidx, i_endidx=i_endidx,                           & !in
        &           t_2m=prm_diag%t_2m(:,jb), qv_2m=prm_diag%qv_2m(:,jb),               & !out
        &           td_2m=prm_diag%td_2m(:,jb), rh_2m=prm_diag%rh_2m(:,jb),             & !out
        &           u_10m=prm_diag%u_10m(:,jb), v_10m=prm_diag%v_10m(:,jb)              ) !out

      ! dynamic gusts
      DO jc = i_startidx, i_endidx
        prm_diag%dyn_gust(jc,jb) = nwp_dyn_gust(prm_diag%u_10m(jc,jb), prm_diag%v_10m(jc,jb), &
          &  prm_diag%tcm(jc,jb), p_diag%u(jc,nlev,jb), p_diag%v(jc,nlev,jb),                 &
          &  p_diag%u(jc,jk_gust(jc),jb), p_diag%v(jc,jk_gust(jc),jb),p_metrics%mask_mtnpoints_g(jc,jb) )
      ENDDO

      ! instantaneous max/min 2m temperature over tiles (trivial operation for 1 tile)
      prm_diag%t_tilemax_inst_2m(i_startidx:i_endidx,jb) = prm_diag%t_2m(i_startidx:i_endidx,jb)
      prm_diag%t_tilemin_inst_2m(i_startidx:i_endidx,jb) = prm_diag%t_2m(i_startidx:i_endidx,jb)


      ! -- end GME code

      ! Aggregate tile-based output fields over tiles - initialize to zero before summation
      prm_diag%gz0 (i_startidx:i_endidx,jb)        = 0._wp
      prm_diag%tcm (i_startidx:i_endidx,jb)        = 0._wp
      prm_diag%tch (i_startidx:i_endidx,jb)        = 0._wp
      prm_diag%tkvm(i_startidx:i_endidx,nlevp1,jb) = 0._wp
      prm_diag%tkvh(i_startidx:i_endidx,nlevp1,jb) = 0._wp
        
      DO jt = 1, ntiles_total+ntiles_water
        DO jc = i_startidx, i_endidx

          ! Aggregate
          area_frac = ext_data%atm%frac_t(jc,jb,jt)

          prm_diag%gz0(jc,jb)         = prm_diag%gz0(jc,jb)        +prm_diag%gz0_t   (jc,jb,jt) * area_frac
          prm_diag%tcm(jc,jb)         = prm_diag%tcm(jc,jb)        +prm_diag%tcm_t   (jc,jb,jt) * area_frac
          prm_diag%tch(jc,jb)         = prm_diag%tch(jc,jb)        +prm_diag%tch_t   (jc,jb,jt) * area_frac
          prm_diag%tkvm(jc,nlevp1,jb) = prm_diag%tkvm(jc,nlevp1,jb)+prm_diag%tkvm_s_t(jc,jb,jt) * area_frac
          prm_diag%tkvh(jc,nlevp1,jb) = prm_diag%tkvh(jc,nlevp1,jb)+prm_diag%tkvh_s_t(jc,jb,jt) * area_frac

          ! the EDMF turbulence scheme does not have tfv. Set tfv=1
          prm_diag%tfv_t(jc,jb,jt)    = 1._wp    !   prm_diag%tfv(jc,jb)

          ! Copy transfer u_10m/v_10m to tile-based variables, which are used in TERRA
          prm_diag%u_10m_t(jc,jb,jt)  = prm_diag%u_10m(jc,jb)
          prm_diag%v_10m_t(jc,jb,jt)  = prm_diag%v_10m(jc,jb)

        ENDDO
      ENDDO


    END SELECT !inwp_turb


    ! Compute wind speed in 10m
    ! used by mo_albedo (albedo_whitecap=1)
    ! 
    !$acc parallel default(present) if(lzacc)
    !$acc loop gang vector
    DO jc = i_startidx, i_endidx
      prm_diag%sp_10m(jc,jb) = SQRT(prm_diag%u_10m(jc,jb)**2 + prm_diag%v_10m(jc,jb)**2 )
    ENDDO
    !$acc end parallel

  ENDDO ! jb
!$OMP END DO NOWAIT
!$OMP END PARALLEL

  ! Serialbox2 input fields serialization
  !$ser verbatim IF(lzacc) THEN
  !$ser verbatim     CALL serialize_turbtrans_interface_output(jg, nproma, nlev,&
  !$ser verbatim                                               p_prog, p_prog_rcf, p_diag, prm_diag,&
  !$ser verbatim                                               lnd_prog_new, lnd_diag)
  !$ser verbatim ENDIF

!$acc end data

END SUBROUTINE nwp_turbtrans

END MODULE mo_nwp_turbtrans_interface
