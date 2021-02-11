!>
!! mo_sgs_turbulence
!!
!! Calculates 3D subgrid-scale viscosity and diffusivity in the nonhydrostatic model
!!
!! @author Anurag Dipankar, MPI-M
!!
!!
!! @par Revision History
!! Initial release by Anurag Dipankar, MPI-M (2013-02-20)
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

MODULE mo_sgs_turbulence

  USE mo_kind,                ONLY: wp
  USE mo_exception,           ONLY: message, finish,message_text, debug_messages_on
  USE mo_nonhydro_types,      ONLY: t_nh_prog, t_nh_diag, t_nh_metrics
  USE mo_model_domain,        ONLY: t_patch
  USE mo_intp_data_strc,      ONLY: t_int_state
  USE mo_intp_rbf,            ONLY: rbf_vec_interpol_vertex, rbf_vec_interpol_edge
  USE mo_intp,                ONLY: verts2edges_scalar, edges2verts_scalar, &
                                    cells2verts_scalar, cells2edges_scalar, &
                                    edges2cells_scalar, verts2cells_scalar
  USE mo_intp_rbf,            ONLY: rbf_vec_interpol_cell
  USE mo_parallel_config,     ONLY: nproma, p_test_run
  USE mo_run_config,          ONLY: iqv, iqc, msg_level
  USE mo_loopindices,         ONLY: get_indices_e, get_indices_c, get_indices_v
  USE mo_impl_constants    ,  ONLY: min_rledge, min_rlcell, min_rlvert, &
                                    min_rledge_int, min_rlcell_int, min_rlvert_int
  USE mo_math_constants,      ONLY: dbl_eps, pi
  USE mo_math_utilities,      ONLY: tdma_solver
  USE mo_sync,                ONLY: SYNC_E, SYNC_C, SYNC_V, sync_patch_array, &
                                    sync_patch_array_mult
  USE mo_physical_constants,  ONLY: cpd, rcvd, p0ref, grav, rcpd, alv
  USE mo_nwp_lnd_types,       ONLY: t_lnd_prog, t_wtr_prog, t_lnd_diag
  USE mo_surface_les,         ONLY: surface_conditions
  USE mo_nwp_phy_types,       ONLY: t_nwp_phy_diag, t_nwp_phy_tend
  USE mo_les_config,          ONLY: les_config
  USE mo_impl_constants_grf,  ONLY: grf_bdywidth_c, grf_bdywidth_e
  USE mo_util_dbg_prnt,       ONLY: dbg_print
  USE mo_turbulent_diagnostic,ONLY: is_sampling_time, idx_sgs_th_flx, &
                                    idx_sgs_qv_flx, idx_sgs_qc_flx,   &
                                    idx_sgs_u_flx, idx_sgs_v_flx
  USE mo_statistics,          ONLY: levels_horizontal_mean
  USE mo_les_utilities,       ONLY: brunt_vaisala_freq, vert_intp_full2half_cell_3d
  USE mo_fortran_tools,       ONLY: init

  IMPLICIT NONE

  PRIVATE
  REAL(wp),         PARAMETER :: z_1by3  = 1._wp/3._wp

  !Parameter for vertical scheme type
  INTEGER, PARAMETER :: iexplicit = 1
  INTEGER, PARAMETER :: iimplicit = 2

  PUBLIC :: drive_subgrid_diffusion

  !Variables for the module
  REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: visc_smag_iv, visc_smag_ie, visc_smag_c, &
                                             rho_ic, DIV_c, u_vert, v_vert, &
                                             w_ie, w_vert

  CHARACTER(len=*), PARAMETER :: inmodule = 'mo_sgs_turbulence:'
  CHARACTER(len=12)  :: str_module = 'sgs_turb'  ! Output of module for 1 line debug
  INTEGER            :: idt_src    = 4           ! Determines level of detail for 1 line debug

  CONTAINS

  !-------------------------------------------------------------------------------------
  !>
  !! drive_subgrid_diffusion
  !!------------------------------------------------------------------------
  !! Driver for computing the sgs viscosity and diffusivity using Smagorisnky model
  !! and evaluating diffusion term on triangles
  !!------------------------------------------------------------------------
  !! @par Revision History
  !! Initial release by Anurag Dipankar, MPI-M (2013-03-05)
  SUBROUTINE drive_subgrid_diffusion(p_nh_prog, p_nh_prog_rcf, p_nh_diag, p_nh_metrics,&
                                     p_patch, p_int, p_prog_lnd_now, p_prog_lnd_new,   &
                                     p_diag_lnd, prm_diag, prm_nwp_tend, dt)

    TYPE(t_nh_prog),   INTENT(inout)     :: p_nh_prog     !< single nh prognostic state
    TYPE(t_nh_prog),   INTENT(in)        :: p_nh_prog_rcf !< rcf nh prognostic state
    TYPE(t_nh_diag),   INTENT(inout)     :: p_nh_diag     !< single nh diagnostic state
    TYPE(t_nh_metrics),INTENT(inout),TARGET :: p_nh_metrics  !< single nh metric state, SBr, CHa: in->inout
    TYPE(t_patch),     INTENT(in),TARGET :: p_patch       !< single patch
    TYPE(t_int_state), INTENT(in),TARGET :: p_int         !< single interpolation state
    TYPE(t_lnd_prog),  INTENT(in)        :: p_prog_lnd_now!<land prog state
    TYPE(t_lnd_prog),  INTENT(inout)     :: p_prog_lnd_new!<land prog state
    TYPE(t_lnd_diag),  INTENT(inout)     :: p_diag_lnd    !<land diag state
    TYPE(t_nwp_phy_diag),   INTENT(inout):: prm_diag      !< atm phys vars
    TYPE(t_nwp_phy_tend), TARGET,INTENT(inout):: prm_nwp_tend    !< atm tend vars
    REAL(wp),          INTENT(in)        :: dt

    REAL(wp), ALLOCATABLE :: theta(:,:,:), theta_v(:,:,:)

    INTEGER :: nlev, nlevp1
    INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx, i_nchdom
    INTEGER :: rl_start, rl_end
    INTEGER :: jk, jb, jc, jg


    !CALL debug_messages_on

    IF (msg_level >= 15) &
         CALL message(TRIM(inmodule), 'drive_subgrid_diffusion')

    jg = p_patch%id

    nlev   = p_patch%nlev
    nlevp1 = nlev+1
    i_nchdom   = MAX(1,p_patch%n_childdom)

    ALLOCATE( u_vert(nproma,nlev,p_patch%nblks_v),           &
              v_vert(nproma,nlev,p_patch%nblks_v),           &
              w_vert(nproma,nlevp1,p_patch%nblks_v),         &
              w_ie(nproma,nlevp1,p_patch%nblks_e),           &
              visc_smag_iv(nproma,nlevp1,p_patch%nblks_v),   &
              visc_smag_c(nproma,nlev,p_patch%nblks_c),      &
              visc_smag_ie(nproma,nlevp1,p_patch%nblks_e),   &
              theta(nproma,nlev,p_patch%nblks_c),            &
              theta_v(nproma,nlev,p_patch%nblks_c),          &
              DIV_c(nproma,nlev,p_patch%nblks_c),            &
              rho_ic(nproma,nlevp1,p_patch%nblks_c)          &
             )

    !Initialize

!$OMP PARALLEL
    CALL init(visc_smag_iv(:,:,:))
    CALL init(visc_smag_c(:,:,:))
    CALL init(visc_smag_ie(:,:,:))

    IF(p_test_run)THEN
      CALL init(u_vert(:,:,:))
      CALL init(v_vert(:,:,:))
      CALL init(w_vert(:,:,:))
    END IF
!$OMP END PARALLEL

    !Convert temperature to potential temperature: all routines within
    !use theta.

    rl_start   = 1
    rl_end     = min_rlcell
    i_startblk = p_patch%cells%start_blk(rl_start,1)
    i_endblk   = p_patch%cells%end_blk(rl_end,i_nchdom)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,i_startidx,i_endidx)
    DO jb = i_startblk,i_endblk
       CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                          i_startidx, i_endidx, rl_start, rl_end)
       DO jc = i_startidx, i_endidx
         theta(jc,1:nlev,jb)  = p_nh_diag%temp(jc,1:nlev,jb) / &
                                p_nh_prog%exner(jc,1:nlev,jb)

         theta_v(jc,1:nlev,jb)  = p_nh_diag%tempv(jc,1:nlev,jb) / &
                                  p_nh_prog%exner(jc,1:nlev,jb)
       END DO
    END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    !Get rho at interfaces to be used later
    CALL vert_intp_full2half_cell_3d(p_patch, p_nh_metrics, p_nh_prog%rho, rho_ic, &
                                     2, min_rlcell_int-2)

    CALL surface_conditions(p_nh_metrics, p_patch, p_nh_diag, p_int, p_prog_lnd_now, &
                            p_prog_lnd_new, p_diag_lnd, prm_diag, theta,             &
                            p_nh_prog%tracer(:,:,:,iqv))

    !Calculate Brunt Vaisala Frequency
    CALL brunt_vaisala_freq(p_patch, p_nh_metrics, theta_v, prm_diag%bruvais)

    CALL smagorinsky_model(p_nh_prog, p_nh_diag, p_nh_metrics, p_patch, p_int, prm_diag, theta)

    CALL diffuse_hori_velocity(p_nh_prog, p_nh_diag, p_nh_metrics, p_patch, p_int, prm_diag, &
                               prm_nwp_tend%ddt_u_turb, prm_nwp_tend%ddt_v_turb, dt)

    !Vertical velocity is updated here
    CALL diffuse_vert_velocity(p_nh_prog, p_nh_diag, p_nh_metrics, p_patch, p_int, &
                               prm_diag%tkvm, prm_nwp_tend%ddt_w_turb, dt)

    CALL diffuse_scalar(theta, p_nh_metrics, p_patch, p_int, p_nh_diag,  &
                        prm_nwp_tend%ddt_temp_turb, p_nh_prog%exner,     &
                        prm_diag, p_nh_prog%rho, dt, 'theta')

    !For qv and qc: implement for qr as well
    !SBr,CHa: let moisture go into atmosphere for dry CBL setup
    !IF(.NOT.les_config(jg)%is_dry_cbl)THEN
      CALL diffuse_scalar(p_nh_prog%tracer(:,:,:,iqv), p_nh_metrics, p_patch, p_int, &
                          p_nh_diag, prm_nwp_tend%ddt_tracer_turb(:,:,:,iqv),        &
                          p_nh_prog%exner, prm_diag, p_nh_prog%rho, dt, 'qv')

      CALL diffuse_scalar(p_nh_prog%tracer(:,:,:,iqc), p_nh_metrics, p_patch, p_int, &
                          p_nh_diag, prm_nwp_tend%ddt_tracer_turb(:,:,:,iqc),        &
                          p_nh_prog%exner, prm_diag, p_nh_prog%rho, dt, 'qc')
    !ELSE
!$OMP PARALLEL
    !  CALL init(prm_nwp_tend%ddt_tracer_turb(:,:,:,iqv))
    !  CALL init(prm_nwp_tend%ddt_tracer_turb(:,:,:,iqc))
!$OMP END PARALLEL
    !END IF


    DEALLOCATE(u_vert, v_vert, w_vert, w_ie, visc_smag_iv, visc_smag_ie,  &
               theta, visc_smag_c, rho_ic, DIV_c, theta_v)


  END SUBROUTINE drive_subgrid_diffusion
  !-------------------------------------------------------------------------------------


  !-------------------------------------------------------------------------------------
  !>
  !! smagorisnky_model
  !!------------------------------------------------------------------------
  !! Computes the sgs viscosity and diffusivity using Smagorisnky model
  !! \tau_ij = KD_ij where D_ij = du_i/dx_j+du_j/dx_i
  !!
  !! and  K = cs*\Delta**2*D/sqrt(2) where D = sqrt(D_ijD_ij)
  !!
  !! and D**2 = D_11**2 + D_22**2 + D_33**2 + 2D_12**2 + 2D_13**2 + 2D_23**2
  !!
  !! where, D_11 = 2du_1/dx_1
  !!        D_22 = 2d_u2/dx_2
  !!        D_33 = 2d_u3/dx_3
  !!        D_12 = du_1/dx_2 + du_2/dx_1
  !!        D_13 = du_1/dx_3 + du_3/dx_1
  !!        D_23 = du_2/dx_3 + du_3/dx_2
  !! For triangles: 1=normal, 2=tangential, and 3 = z directions
  !!------------------------------------------------------------------------
  !! @par Revision History
  !! Initial release by Anurag Dipankar, MPI-M (2013-02-20)
  SUBROUTINE smagorinsky_model(p_nh_prog, p_nh_diag, p_nh_metrics, p_patch, p_int, &
                               prm_diag, theta)

    TYPE(t_patch),     INTENT(in),TARGET :: p_patch    !< single patch
    TYPE(t_int_state), INTENT(in),TARGET :: p_int      !< single interpolation state
    TYPE(t_nh_prog),   INTENT(inout)     :: p_nh_prog  !< single nh prognostic state
    TYPE(t_nh_diag),   INTENT(in)     :: p_nh_diag  !< single nh diagnostic state
    TYPE(t_nh_metrics),INTENT(inout),TARGET :: p_nh_metrics  !< single nh metric state, SBr, CHa: in->inout
    TYPE(t_nwp_phy_diag),   INTENT(inout):: prm_diag      !< atm phys vars
    REAL(wp),          INTENT(in)        :: theta(:,:,:)   !< atm phys vars, SBr, CHa

    ! local variables
    REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: div_of_stress, mech_prod_e, &
                                               vn_ie, vt_ie
    REAL(wp) :: vn_vert1, vn_vert2, vn_vert3, vn_vert4, vt_vert1, vt_vert2, vt_vert3, &
                vt_vert4, w_full_c1
    REAL(wp) :: w_full_c2, w_full_v1, w_full_v2
    REAL(wp) :: D_11, D_12, D_13, D_22, D_23, D_33, z1, z2
    REAL(wp), POINTER :: diff_smag_ic(:,:,:)

    INTEGER  :: nlev, nlevp1             !< number of full levels
    INTEGER,  DIMENSION(:,:,:), POINTER :: ividx, ivblk, iecidx, iecblk, ieidx, ieblk
    INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx, i_nchdom
    INTEGER :: rl_start, rl_end, jg
    INTEGER :: jk, jb, jc, je, ic, jkp1, jv, jkm1, jcn, jbn, jvn

    !--------------------------------------------------------------------------

    IF (msg_level >= 18) &
         CALL message(TRIM(inmodule), 'smagorinsky_model')

    jg = p_patch%id

    ! number of vertical levels
    nlev   = p_patch%nlev
    nlevp1 = nlev+1
    i_nchdom   = MAX(1,p_patch%n_childdom)

    !Allocation
    ALLOCATE( vn_ie(nproma,nlevp1,p_patch%nblks_e),        &
              vt_ie(nproma,nlevp1,p_patch%nblks_e),        &
              mech_prod_e(nproma,nlev,p_patch%nblks_e),    &
              div_of_stress(nproma,nlev,p_patch%nblks_e)   &
            )

    diff_smag_ic => prm_diag%tkvh

    !Initialize
    IF(p_test_run)THEN
      diff_smag_ic(:,:,:) = 0._wp
    END IF


    !--------------------------------------------------------------------------
    !1) Interpolate velocities at desired locations- mostly around the quadrilateral
    !
    !<It assumes that prog values are all synced at this stage while diag values might not>
    !--------------------------------------------------------------------------

    CALL cells2verts_scalar(p_nh_prog%w, p_patch, p_int%cells_aw_verts, w_vert, opt_rlend=min_rlvert_int)
    CALL cells2edges_scalar(p_nh_prog%w, p_patch, p_int%c_lin_e, w_ie, opt_rlend=min_rledge_int-2)

    ! RBF reconstruction of velocity at vertices: include halos
    CALL rbf_vec_interpol_vertex( p_nh_prog%vn, p_patch, p_int, &
                                  u_vert, v_vert, opt_rlend=min_rlvert_int )

    !sync them
    CALL sync_patch_array_mult(SYNC_V, p_patch, 3, w_vert, u_vert, v_vert)

    !Get vn at interfaces and then get vt at interfaces
    !Boundary values are extrapolated like dynamics although
    !they are not required in current implementation

!$OMP PARALLEL PRIVATE (rl_start,rl_end,i_startblk,i_endblk)
    rl_start = 2
    rl_end   = min_rledge_int-3

    i_startblk = p_patch%edges%start_blk(rl_start,1)
    i_endblk   = p_patch%edges%end_blk(rl_end,i_nchdom)

!$OMP DO PRIVATE(jb,jk,je,i_startidx,i_endidx)
    DO jb = i_startblk,i_endblk
       CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
                          i_startidx, i_endidx, rl_start, rl_end)

#ifdef __LOOP_EXCHANGE
       DO je = i_startidx, i_endidx
         DO jk = 2 , nlev
#else
       DO jk = 2, nlev
         DO je = i_startidx, i_endidx
#endif
            vn_ie(je,jk,jb) = p_nh_metrics%wgtfac_e(je,jk,jb)*p_nh_prog%vn(je,jk,jb) &
                    + (1._wp - p_nh_metrics%wgtfac_e(je,jk,jb))*p_nh_prog%vn(je,jk-1,jb)
         END DO
       END DO
       DO je = i_startidx, i_endidx
          vn_ie(je,1,jb) =                              &
            p_nh_metrics%wgtfacq1_e(je,1,jb)*p_nh_prog%vn(je,1,jb) + &
            p_nh_metrics%wgtfacq1_e(je,2,jb)*p_nh_prog%vn(je,2,jb) + &
            p_nh_metrics%wgtfacq1_e(je,3,jb)*p_nh_prog%vn(je,3,jb)

          vn_ie(je,nlevp1,jb) =                           &
            p_nh_metrics%wgtfacq_e(je,1,jb)*p_nh_prog%vn(je,nlev,jb) +   &
            p_nh_metrics%wgtfacq_e(je,2,jb)*p_nh_prog%vn(je,nlev-1,jb) + &
            p_nh_metrics%wgtfacq_e(je,3,jb)*p_nh_prog%vn(je,nlev-2,jb)
       END DO
    END DO
!$OMP END DO

    CALL rbf_vec_interpol_edge(vn_ie, p_patch, p_int, vt_ie, opt_rlstart=3, &
                               opt_rlend=min_rledge_int-2)

    !--------------------------------------------------------------------------
    !2) Compute horizontal strain rate tensor at full levels
    !--------------------------------------------------------------------------
    ividx => p_patch%edges%vertex_idx
    ivblk => p_patch%edges%vertex_blk

    iecidx => p_patch%edges%cell_idx
    iecblk => p_patch%edges%cell_blk

    ieidx => p_patch%cells%edge_idx
    ieblk => p_patch%cells%edge_blk


    rl_start = 4
    rl_end   = min_rledge_int-2

    i_startblk = p_patch%edges%start_blk(rl_start,1)
    i_endblk   = p_patch%edges%end_blk(rl_end,i_nchdom)

!$OMP DO PRIVATE(jb,jk,je,i_startidx,i_endidx,vn_vert1,vn_vert2,vn_vert3,vn_vert4,  &
!$OMP            vt_vert1,vt_vert2,vt_vert3,vt_vert4,w_full_c1,w_full_c2,w_full_v1,  &
!$OMP            w_full_v2,D_11,D_12,D_13,D_22,D_23,D_33,jkp1)
    DO jb = i_startblk,i_endblk

       CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
                          i_startidx, i_endidx, rl_start, rl_end)
#ifdef __LOOP_EXCHANGE
       DO je = i_startidx, i_endidx
         DO jk = 1 , nlev
#else
       DO jk = 1, nlev
         DO je = i_startidx, i_endidx
#endif
            jkp1 = jk + 1

            vn_vert1 = u_vert(ividx(je,jb,1),jk,ivblk(je,jb,1)) * &
                       p_patch%edges%primal_normal_vert(je,jb,1)%v1 + &
                       v_vert(ividx(je,jb,1),jk,ivblk(je,jb,1)) * &
                       p_patch%edges%primal_normal_vert(je,jb,1)%v2

            vn_vert2 = u_vert(ividx(je,jb,2),jk,ivblk(je,jb,2)) * &
                       p_patch%edges%primal_normal_vert(je,jb,2)%v1 + &
                       v_vert(ividx(je,jb,2),jk,ivblk(je,jb,2)) * &
                       p_patch%edges%primal_normal_vert(je,jb,2)%v2

            vt_vert2 = u_vert(ividx(je,jb,2),jk,ivblk(je,jb,2)) * &
                       p_patch%edges%dual_normal_vert(je,jb,2)%v1 + &
                       v_vert(ividx(je,jb,2),jk,ivblk(je,jb,2)) * &
                       p_patch%edges%dual_normal_vert(je,jb,2)%v2

            vt_vert1 = u_vert(ividx(je,jb,1),jk,ivblk(je,jb,1)) * &
                       p_patch%edges%dual_normal_vert(je,jb,1)%v1 + &
                       v_vert(ividx(je,jb,1),jk,ivblk(je,jb,1)) * &
                       p_patch%edges%dual_normal_vert(je,jb,1)%v2

            vn_vert3 = u_vert(ividx(je,jb,3),jk,ivblk(je,jb,3)) * &
                       p_patch%edges%primal_normal_vert(je,jb,3)%v1 + &
                       v_vert(ividx(je,jb,3),jk,ivblk(je,jb,3)) * &
                       p_patch%edges%primal_normal_vert(je,jb,3)%v2

            vn_vert4 = u_vert(ividx(je,jb,4),jk,ivblk(je,jb,4)) * &
                       p_patch%edges%primal_normal_vert(je,jb,4)%v1 + &
                       v_vert(ividx(je,jb,4),jk,ivblk(je,jb,4)) * &
                       p_patch%edges%primal_normal_vert(je,jb,4)%v2

            vt_vert4 = u_vert(ividx(je,jb,4),jk,ivblk(je,jb,4)) * &
                       p_patch%edges%dual_normal_vert(je,jb,4)%v1 + &
                       v_vert(ividx(je,jb,4),jk,ivblk(je,jb,4)) * &
                       p_patch%edges%dual_normal_vert(je,jb,4)%v2

            vt_vert3 = u_vert(ividx(je,jb,3),jk,ivblk(je,jb,3)) * &
                       p_patch%edges%dual_normal_vert(je,jb,3)%v1 + &
                       v_vert(ividx(je,jb,3),jk,ivblk(je,jb,3)) * &
                       p_patch%edges%dual_normal_vert(je,jb,3)%v2

            !W at full levels
            w_full_c1  = 0.5_wp * (                                        &
                         p_nh_prog%w(iecidx(je,jb,1),jk,iecblk(je,jb,1)) + &
                         p_nh_prog%w(iecidx(je,jb,1),jkp1,iecblk(je,jb,1)) )

            w_full_c2  = 0.5_wp * (                                        &
                         p_nh_prog%w(iecidx(je,jb,2),jk,iecblk(je,jb,2)) + &
                         p_nh_prog%w(iecidx(je,jb,2),jkp1,iecblk(je,jb,2)) )

            !W at full levels vertices from w at vertices at interface levels
            w_full_v1  = 0.5_wp * (                                    &
                         w_vert(ividx(je,jb,1),jk,ivblk(je,jb,1)) +    &
                         w_vert(ividx(je,jb,1),jkp1,ivblk(je,jb,1)) )

            w_full_v2  = 0.5_wp * (                                    &
                         w_vert(ividx(je,jb,2),jk,ivblk(je,jb,2)) +    &
                         w_vert(ividx(je,jb,2),jkp1,ivblk(je,jb,2)) )

            !Strain rates at edge center
            D_11       =      2._wp * (vn_vert4-vn_vert3) * &
                              p_patch%edges%inv_vert_vert_length(je,jb)

            D_12       =     p_patch%edges%tangent_orientation(je,jb) *  &
                            (vn_vert2-vn_vert1) *                       &
                             p_patch%edges%inv_primal_edge_length(je,jb)&
                             + (vt_vert4-vt_vert3)*                     &
                             p_patch%edges%inv_vert_vert_length(je,jb)

            D_13       =     (vn_ie(je,jk,jb)-vn_ie(je,jkp1,jb))* &
                              p_nh_metrics%inv_ddqz_z_full_e(je,jk,jb)  +   &
                             (w_full_c2 - w_full_c1) *                      &
                              p_patch%edges%inv_dual_edge_length(je,jb)

            D_22       =     2._wp*(vt_vert2-vt_vert1)*p_patch%edges%tangent_orientation(je,jb) * &
                             p_patch%edges%inv_primal_edge_length(je,jb)

            D_23       =    (vt_ie(je,jk,jb) - vt_ie(je,jkp1,jb)) *          &
                             p_nh_metrics%inv_ddqz_z_full_e(je,jk,jb)  +     &
                             p_patch%edges%tangent_orientation(je,jb) *       &
                            (w_full_v2 - w_full_v1) *                        &
                             p_patch%edges%inv_primal_edge_length(je,jb)

            D_33       =     2._wp * (w_ie(je,jk,jb) - w_ie(je,jkp1,jb)) *   &
                             p_nh_metrics%inv_ddqz_z_full_e(je,jk,jb)

            !Mechanical prod is half of this value divided by km
            mech_prod_e(je,jk,jb) = D_11**2 + D_22**2 + D_33**2 + &
                                     2._wp * ( D_12**2 + D_13**2 + D_23**2 )

            !calculate divergence to get the deviatoric part of stress tensor in
            !diffusion: D_11-1/3*(D_11+D_22+D_33)
            div_of_stress(je,jk,jb)  = 0.5_wp * ( D_11 + D_22 + D_33 )

         ENDDO
       ENDDO
    ENDDO
!$OMP END DO

    !Interpolate mech production term from mid level edge to interface level cell
    !except top and bottom boundaries
    rl_start = 3
    rl_end   = min_rlcell_int-1
    i_startblk = p_patch%cells%start_blk(rl_start,1)
    i_endblk   = p_patch%cells%end_blk(rl_end,i_nchdom)

!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx,jkm1)
    DO jb = i_startblk,i_endblk
       CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                          i_startidx, i_endidx, rl_start, rl_end)
#ifdef __LOOP_EXCHANGE
       DO jc = i_startidx, i_endidx
         DO jk = 2 , nlev
#else
       DO jk = 2, nlev
         DO jc = i_startidx, i_endidx
#endif
            jkm1  = jk - 1
            prm_diag%mech_prod(jc,jk,jb) = p_nh_metrics%wgtfac_c(jc,jk,jb) * (          &
              mech_prod_e(ieidx(jc,jb,1),jk,ieblk(jc,jb,1))*p_int%e_bln_c_s(jc,1,jb)  + &
              mech_prod_e(ieidx(jc,jb,2),jk,ieblk(jc,jb,2))*p_int%e_bln_c_s(jc,2,jb)  + &
              mech_prod_e(ieidx(jc,jb,3),jk,ieblk(jc,jb,3))*p_int%e_bln_c_s(jc,3,jb) )+ &
                                 (1._wp - p_nh_metrics%wgtfac_c(jc,jk,jb)) * (          &
              mech_prod_e(ieidx(jc,jb,1),jkm1,ieblk(jc,jb,1))*p_int%e_bln_c_s(jc,1,jb)+ &
              mech_prod_e(ieidx(jc,jb,2),jkm1,ieblk(jc,jb,2))*p_int%e_bln_c_s(jc,2,jb)+ &
              mech_prod_e(ieidx(jc,jb,3),jkm1,ieblk(jc,jb,3))*p_int%e_bln_c_s(jc,3,jb) )
         END DO
       END DO
    END DO
!$OMP END DO

    !div_of_stress from edge to cell-scalar interpolation
    rl_start = grf_bdywidth_c+1
    rl_end   = min_rlcell_int-1 !-1 for its use in hor diffusion
    i_startblk = p_patch%cells%start_blk(rl_start,1)
    i_endblk   = p_patch%cells%end_blk(rl_end,i_nchdom)

!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx)
    DO jb = i_startblk,i_endblk
       CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                          i_startidx, i_endidx, rl_start, rl_end)
#ifdef __LOOP_EXCHANGE
       DO jc = i_startidx, i_endidx
         DO jk = 1 , nlev
#else
       DO jk = 1 , nlev
         DO jc = i_startidx, i_endidx
#endif
           DIV_c(jc,jk,jb) =                                                               &
              (div_of_stress(ieidx(jc,jb,1),jk,ieblk(jc,jb,1))*p_int%e_bln_c_s(jc,1,jb)  + &
               div_of_stress(ieidx(jc,jb,2),jk,ieblk(jc,jb,2))*p_int%e_bln_c_s(jc,2,jb)  + &
               div_of_stress(ieidx(jc,jb,3),jk,ieblk(jc,jb,3))*p_int%e_bln_c_s(jc,3,jb))
         END DO
       END DO
    END DO
!$OMP END DO

    !--------------------------------------------------------------------------
    !3) Classical Smagorinsky model with stability correction due to Lilly 1962
    !   at interface cell centers. At this point mech_prod is twice the actual
    !   mechanical production term.
    !--------------------------------------------------------------------------
    !MP = Mechanical prod term calculated above
    !visc = mixing_length_sq * SQRT(MP)/2 * SQRT(1-Ri/Pr) where
    !Ri = (g/theta)*d_theta_dz/(MP/2), where Brunt_vaisala_freq (byncy prod term/kh)
    !   = (g/theta)*d_theta_dz.
    !After simplification: visc = mixing_length_sq/SQRT(2) * SQRT[MP/2 - (Brunt_vaisala_frq/Pr)]
    !Note that the factor SQRT(2) with mixing_length_sq is considered into the Smag constant

    rl_start = 3
    rl_end   = min_rlcell_int
    i_startblk = p_patch%cells%start_blk(rl_start,1)
    i_endblk   = p_patch%cells%end_blk(rl_end,i_nchdom)

!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx)
    DO jb = i_startblk,i_endblk
       CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                          i_startidx, i_endidx, rl_start, rl_end)
#ifdef __LOOP_EXCHANGE
       DO jc = i_startidx, i_endidx
         DO jk = 2 , nlev
#else
       DO jk = 2 , nlev
         DO jc = i_startidx, i_endidx
#endif
           diff_smag_ic(jc,jk,jb) = rho_ic(jc,jk,jb) * les_config(jg)%rturb_prandtl * &
                     p_nh_metrics%mixing_length_sq(jc,jk,jb)             * &
                     SQRT(MAX(0._wp, prm_diag%mech_prod(jc,jk,jb)*0.5_wp - &
                     les_config(jg)%rturb_prandtl*prm_diag%bruvais(jc,jk,jb)))  
         END DO
       END DO
       DO jc = i_startidx, i_endidx
          diff_smag_ic(jc,1,jb)      = diff_smag_ic(jc,2,jb)
          diff_smag_ic(jc,nlevp1,jb) = diff_smag_ic(jc,nlev,jb)
       END DO
    END DO
!$OMP END DO
!$OMP END PARALLEL

    CALL sync_patch_array(SYNC_C, p_patch, diff_smag_ic)

    IF (is_sampling_time) THEN
!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx)
      DO jb = i_startblk,i_endblk
         CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                          i_startidx, i_endidx, rl_start, rl_end)
#ifdef __LOOP_EXCHANGE
        DO jc = i_startidx, i_endidx
          DO jk = 2 , nlev
#else 
        DO jk = 2 , nlev
          DO jc = i_startidx, i_endidx
#endif
            p_nh_metrics%wthd(jc,jk,jb) = -diff_smag_ic(jc,jk,jb) * (theta(jc,jk-1,jb) - theta(jc,jk,jb)) * &
                                   p_nh_metrics%inv_ddqz_z_half(jc,jk,jb)
          END DO
        END DO
        DO jc = i_startidx, i_endidx
            p_nh_metrics%wthd(jc,1,jb) = 0._wp
            !extrapolate for surface density
            p_nh_metrics%wthd(jc,nlevp1,jb) = prm_diag%shfl_s(jc,jb) * rcpd
        END DO
      END DO
    END IF

    !--------------------------------------------------------------------------
    !4) Interpolate difusivity (viscosity) to different locations: calculate them for
    !   halos also because they will be used later in diffusion
    !--------------------------------------------------------------------------

    !4a) visc at cell center
    rl_start = grf_bdywidth_c
    rl_end   = min_rlcell_int-1
    i_startblk = p_patch%cells%start_blk(rl_start,1)
    i_endblk   = p_patch%cells%end_blk(rl_end,i_nchdom)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx)
    DO jb = i_startblk,i_endblk
       CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                          i_startidx, i_endidx, rl_start, rl_end)
#ifdef __LOOP_EXCHANGE
       DO jc = i_startidx, i_endidx
         DO jk = 1 , nlev
#else
       DO jk = 1 , nlev
         DO jc = i_startidx, i_endidx
#endif
           visc_smag_c(jc,jk,jb) = MAX( les_config(jg)%km_min, &
                                  (diff_smag_ic(jc,jk,jb)+diff_smag_ic(jc,jk+1,jb)) * &
                                   0.5_wp * les_config(jg)%turb_prandtl )
         END DO
       END DO
    END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    !4b) visc at vertices
    CALL cells2verts_scalar(diff_smag_ic, p_patch, p_int%cells_aw_verts, visc_smag_iv, &
                            opt_rlstart=5, opt_rlend=min_rlvert_int-1)
    visc_smag_iv = MAX( les_config(jg)%km_min, visc_smag_iv * les_config(jg)%turb_prandtl )

    !4c) Now calculate visc_smag at half levels at edge
    CALL cells2edges_scalar(diff_smag_ic, p_patch, p_int%c_lin_e, visc_smag_ie, &
                            opt_rlstart=grf_bdywidth_e, opt_rlend=min_rledge_int-1)
    visc_smag_ie = MAX( les_config(jg)%km_min, visc_smag_ie * les_config(jg)%turb_prandtl )

    !4d)Get visc_smag_ic
    prm_diag%tkvm = MAX( les_config(jg)%km_min, prm_diag%tkvh * les_config(jg)%turb_prandtl )

    !DEALLOCATE variables
    DEALLOCATE( div_of_stress, mech_prod_e, vn_ie, vt_ie )

  END SUBROUTINE smagorinsky_model
  !-------------------------------------------------------------------------------------


  !-------------------------------------------------------------------------------------
  !>
  !! diffuse_hori_velocity
  !!------------------------------------------------------------------------
  !! Calculate the SGS diffusion term for normal velocity component
  !! - Uses the forward Euler time scheme in time split (sequential) manner adopted
  !!   in the NH version.
  !! - Option to switch on implicit scheme in vertical
  !!
  !! d_vn/d_t =  d_tau_11/d_x1 + d_tau_12/d_x2 + d_tau_13/d_x3
  !!------------------------------------------------------------------------
  !! @par Revision History
  !! Initial release by Anurag Dipankar, MPI-M (2013-02-05)
  SUBROUTINE diffuse_hori_velocity(p_nh_prog, p_nh_diag, p_nh_metrics, p_patch, p_int, &
                                   prm_diag, ddt_u, ddt_v, dt)

    TYPE(t_nh_prog),   INTENT(in)        :: p_nh_prog    !< single nh prognostic state
    TYPE(t_nh_diag),   INTENT(in)        :: p_nh_diag    !< single nh diagnostic state
    TYPE(t_nh_metrics),INTENT(in),TARGET :: p_nh_metrics !< single nh metric state
    TYPE(t_patch), TARGET, INTENT(in)    :: p_patch      !< single patch
    TYPE(t_int_state), INTENT(in),TARGET :: p_int        !< single interpolation state
    TYPE(t_nwp_phy_diag),INTENT(inout)   :: prm_diag     !< atm phys vars
    REAL(wp),   TARGET, INTENT(inout)    :: ddt_u(:,:,:) !< u tendency
    REAL(wp),   TARGET, INTENT(inout)    :: ddt_v(:,:,:) !< v tendency
    REAL(wp),           INTENT(in)       :: dt           !< dt turb

    REAL(wp) :: flux_up_e, flux_dn_e, flux_up_v, flux_dn_v, flux_up_c, flux_dn_c
    REAL(wp) :: stress_uc, stress_vc, stress_c1n, stress_c2n, inv_mwind
    REAL(wp) :: vn_vert1, vn_vert2, vn_vert3, vn_vert4, dvt, dwdn, inv_dt
    REAL(wp) :: inv_rhoe(nproma,p_patch%nlev,p_patch%nblks_e)
    REAL(wp) :: vn_new(nproma,p_patch%nlev,p_patch%nblks_e)
    REAL(wp) :: unew(nproma,p_patch%nlev,p_patch%nblks_c)
    REAL(wp) :: vnew(nproma,p_patch%nlev,p_patch%nblks_c)
    REAL(wp) :: tot_tend(nproma,p_patch%nlev,p_patch%nblks_e)

    REAL(wp), DIMENSION(nproma,p_patch%nlev) :: a, b, c, rhs
    REAL(wp), DIMENSION(p_patch%nlev) :: var_new, outvar

    INTEGER,  DIMENSION(:,:,:), POINTER :: ividx, ivblk, iecidx, iecblk
    INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx, i_nchdom
    INTEGER :: rl_start, rl_end
    INTEGER :: jk, jkp1, jb, je, jcn, jbn, jvn, jc
    INTEGER :: nlev, jg, nlevp1

    IF (msg_level >= 18) &
         CALL message(TRIM(inmodule), 'diffuse_hori_velocity')

    !patch id
    jg = p_patch%id

    ! number of vertical levels
    nlev     = p_patch%nlev
    nlevp1   = nlev+1
    i_nchdom = MAX(1,p_patch%n_childdom)

    inv_dt  = 1._wp / dt

    ividx  => p_patch%edges%vertex_idx
    ivblk  => p_patch%edges%vertex_blk

    iecidx => p_patch%edges%cell_idx
    iecblk => p_patch%edges%cell_blk

    !Some initializations
    a = 0._wp; c = 0._wp

    !total tendency
    tot_tend(:,:,:) = 0._wp

    !new vn
    vn_new(:,:,:) = p_nh_prog%vn(:,:,:)


    !density at edge
    CALL cells2edges_scalar(p_nh_prog%rho, p_patch, p_int%c_lin_e, inv_rhoe, &
                            opt_rlstart=grf_bdywidth_e+1, opt_rlend=min_rledge_int)

    rl_start = grf_bdywidth_e+1
    rl_end   = min_rledge_int

    i_startblk = p_patch%edges%start_blk(rl_start,1)
    i_endblk   = p_patch%edges%end_blk(rl_end,i_nchdom)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,je,i_startidx,i_endidx)
    DO jb = i_startblk,i_endblk
      CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, rl_start, rl_end)
#ifdef __LOOP_EXCHANGE
       DO je = i_startidx, i_endidx
         DO jk = 1 , nlev
#else
      DO jk = 1 , nlev
       DO je = i_startidx, i_endidx
#endif
         inv_rhoe(je,jk,jb) = 1._wp / inv_rhoe(je,jk,jb)
       END DO
      END DO
    END DO
!$OMP END DO

    ! 1) First get the horizontal tendencies

!$OMP DO PRIVATE(jb,jk,je,i_startidx,i_endidx,vn_vert1,vn_vert2,vn_vert3,vn_vert4,&
!$OMP   dvt,jcn,jbn,flux_up_c,flux_dn_c,jvn,flux_up_v,flux_dn_v)
    DO jb = i_startblk,i_endblk
      CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, rl_start, rl_end)
#ifdef __LOOP_EXCHANGE
       DO je = i_startidx, i_endidx
         DO jk = 1 , nlev
#else
      DO jk = 1 , nlev
       DO je = i_startidx, i_endidx
#endif
         vn_vert1 = u_vert(ividx(je,jb,1),jk,ivblk(je,jb,1)) * &
                    p_patch%edges%primal_normal_vert(je,jb,1)%v1 + &
                    v_vert(ividx(je,jb,1),jk,ivblk(je,jb,1)) * &
                    p_patch%edges%primal_normal_vert(je,jb,1)%v2

         vn_vert2 = u_vert(ividx(je,jb,2),jk,ivblk(je,jb,2)) * &
                    p_patch%edges%primal_normal_vert(je,jb,2)%v1 + &
                    v_vert(ividx(je,jb,2),jk,ivblk(je,jb,2)) * &
                    p_patch%edges%primal_normal_vert(je,jb,2)%v2

         vn_vert3 = u_vert(ividx(je,jb,3),jk,ivblk(je,jb,3)) * &
                    p_patch%edges%primal_normal_vert(je,jb,3)%v1 + &
                    v_vert(ividx(je,jb,3),jk,ivblk(je,jb,3)) * &
                    p_patch%edges%primal_normal_vert(je,jb,3)%v2

         vn_vert4 = u_vert(ividx(je,jb,4),jk,ivblk(je,jb,4)) * &
                    p_patch%edges%primal_normal_vert(je,jb,4)%v1 + &
                    v_vert(ividx(je,jb,4),jk,ivblk(je,jb,4)) * &
                    p_patch%edges%primal_normal_vert(je,jb,4)%v2

         dvt      = u_vert(ividx(je,jb,4),jk,ivblk(je,jb,4)) * &
                    p_patch%edges%dual_normal_vert(je,jb,4)%v1 + &
                    v_vert(ividx(je,jb,4),jk,ivblk(je,jb,4)) * &
                    p_patch%edges%dual_normal_vert(je,jb,4)%v2 - &
                    (u_vert(ividx(je,jb,3),jk,ivblk(je,jb,3)) * &
                    p_patch%edges%dual_normal_vert(je,jb,3)%v1 + &
                    v_vert(ividx(je,jb,3),jk,ivblk(je,jb,3)) * &
                    p_patch%edges%dual_normal_vert(je,jb,3)%v2)


         !tendency in normal direction:
         !flux = visc*(D_11-2/3DIV) = visc*(2*delta_v/(vert_vert_len/2)-2/3*div_of_stress)

         jcn     = iecidx(je,jb,2)
         jbn     = iecblk(je,jb,2)
         flux_up_c = visc_smag_c(jcn,jk,jbn) * &
                   ( 4._wp*(vn_vert4-p_nh_prog%vn(je,jk,jb))   * &
                     p_patch%edges%inv_vert_vert_length(je,jb) - &
                     2._wp*z_1by3*DIV_c(jcn,jk,jbn) )

         jcn     = iecidx(je,jb,1)
         jbn     = iecblk(je,jb,1)
         flux_dn_c = visc_smag_c(jcn,jk,jbn) * &
                   ( 4._wp*(p_nh_prog%vn(je,jk,jb)-vn_vert3)   * &
                     p_patch%edges%inv_vert_vert_length(je,jb) - &
                     2._wp*z_1by3*DIV_c(jcn,jk,jbn) )

         !tendency in tangential direction

         !D_12 between edge center and the vertex: delta_v/(primal_edge_len/2) +
         ! ((vt4+vt2)/2-(vt3+vt2)/2)/(distance_opp_edges)
         !flux = D_12*visc

         !Note that the tangential velocities at vertices are used in D_12 is an
         !approximation for speed. Better way is to use vt reconstructed from vn at
         !each edges. Also, visc at somewhere between edge mid point and the vertex
         !should be used but this is a good approximation

         jvn     = ividx(je,jb,2)
         jbn     = ivblk(je,jb,2)
         flux_up_v = 0.5_wp * (visc_smag_iv(jvn,jk,jbn)+visc_smag_iv(jvn,jk+1,jbn)) *     &
             ( p_patch%edges%tangent_orientation(je,jb)*(vn_vert2-p_nh_prog%vn(je,jk,jb))* &
              p_patch%edges%inv_primal_edge_length(je,jb)*2._wp +                         &
              dvt*p_patch%edges%inv_vert_vert_length(je,jb) )

         jvn     = ividx(je,jb,1)
         jbn     = ivblk(je,jb,1)
         flux_dn_v = 0.5_wp * (visc_smag_iv(jvn,jk,jbn)+visc_smag_iv(jvn,jk+1,jbn)) *     &
           ( p_patch%edges%tangent_orientation(je,jb)*(p_nh_prog%vn(je,jk,jb)-vn_vert1) *  &
             p_patch%edges%inv_primal_edge_length(je,jb)*2._wp +                          &
             dvt*p_patch%edges%inv_vert_vert_length(je,jb) )


         tot_tend(je,jk,jb) = ( (flux_up_c-flux_dn_c)*p_patch%edges%inv_dual_edge_length(je,jb) + &
                        p_patch%edges%tangent_orientation(je,jb) * (flux_up_v-flux_dn_v)  * &
                        p_patch%edges%inv_primal_edge_length(je,jb) * 2._wp ) * inv_rhoe(je,jk,jb)

       END DO
      END DO
    END DO
!$OMP END DO
!$OMP END PARALLEL

    !
    ! 2) Vertical tendency
    !
    SELECT CASE(les_config(jg)%vert_scheme_type)

    CASE(iexplicit)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jkp1,je,i_startidx,i_endidx,flux_up_e,flux_dn_e,stress_c1n,stress_c2n)
      DO jb = i_startblk,i_endblk
        CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
                           i_startidx, i_endidx, rl_start, rl_end)
#ifdef __LOOP_EXCHANGE
       DO je = i_startidx, i_endidx
         DO jk = 2 , nlev-1
#else
        DO jk = 2 , nlev-1
         DO je = i_startidx, i_endidx
#endif
           !tendency in vertical direction
           jkp1 = jk + 1

           flux_up_e = visc_smag_ie(je,jk,jb) *                             &
                  ( (p_nh_prog%vn(je,jk-1,jb) - p_nh_prog%vn(je,jk,jb)) * &
                     p_nh_metrics%inv_ddqz_z_half_e(je,jk,jb) +           &
                    (p_nh_prog%w(iecidx(je,jb,2),jk,iecblk(je,jb,2)) -    &
                     p_nh_prog%w(iecidx(je,jb,1),jk,iecblk(je,jb,1))) *   &
                     p_patch%edges%inv_dual_edge_length(je,jb) )

           flux_dn_e = visc_smag_ie(je,jkp1,jb) *                           &
                  ( (p_nh_prog%vn(je,jk,jb) - p_nh_prog%vn(je,jkp1,jb)) * &
                     p_nh_metrics%inv_ddqz_z_half_e(je,jkp1,jb) +         &
                    (p_nh_prog%w(iecidx(je,jb,2),jkp1,iecblk(je,jb,2)) -  &
                     p_nh_prog%w(iecidx(je,jb,1),jkp1,iecblk(je,jb,1))) * &
                     p_patch%edges%inv_dual_edge_length(je,jb) )

          tot_tend(je,jk,jb) = tot_tend(je,jk,jb) +  (flux_up_e - flux_dn_e) * &
                            p_nh_metrics%inv_ddqz_z_full_e(je,jk,jb) * inv_rhoe(je,jk,jb)

         END DO
        END DO

         !-----------------------------------------------------------------
         !3) Boundary treatment in vertical: use surface fluxes from
         !   surface_conditions
         !-----------------------------------------------------------------

         !-----------------------------------------------------------------
         !jk = 1
         !-----------------------------------------------------------------
         DO je = i_startidx, i_endidx
           flux_dn_e = visc_smag_ie(je,2,jb) *                         &
                  ( (p_nh_prog%vn(je,1,jb) - p_nh_prog%vn(je,2,jb)) *  &
                     p_nh_metrics%inv_ddqz_z_half_e(je,2,jb) +         &
                    (p_nh_prog%w(iecidx(je,jb,2),2,iecblk(je,jb,2)) -  &
                     p_nh_prog%w(iecidx(je,jb,1),2,iecblk(je,jb,1))) * &
                     p_patch%edges%inv_dual_edge_length(je,jb) )

           tot_tend(je,1,jb) = tot_tend(je,1,jb) - flux_dn_e *  &
                        p_nh_metrics%inv_ddqz_z_full_e(je,1,jb) * inv_rhoe(je,1,jb)
         END DO

         !-----------------------------------------------------------------
         ! jk = nlev
         !-----------------------------------------------------------------
         DO je = i_startidx, i_endidx
           flux_up_e = visc_smag_ie(je,nlev,jb) *                             &
                  ( (p_nh_prog%vn(je,nlev-1,jb) - p_nh_prog%vn(je,nlev,jb)) * &
                     p_nh_metrics%inv_ddqz_z_half_e(je,nlev,jb) +             &
                    (p_nh_prog%w(iecidx(je,jb,2),nlev,iecblk(je,jb,2)) -      &
                     p_nh_prog%w(iecidx(je,jb,1),nlev,iecblk(je,jb,1))) *     &
                     p_patch%edges%inv_dual_edge_length(je,jb) )

           !Get net shear stress in the direction of vn at surface

           !shear stress in normal direction from cell 1
           stress_c1n = prm_diag%umfl_s(iecidx(je,jb,1),iecblk(je,jb,1)) * &
                        p_patch%edges%primal_normal_cell(je,jb,1)%v1     + &
                        prm_diag%vmfl_s(iecidx(je,jb,1),iecblk(je,jb,1)) * &
                        p_patch%edges%primal_normal_cell(je,jb,1)%v2

           !shear stress in normal direction from cell 2
           stress_c2n = prm_diag%umfl_s(iecidx(je,jb,2),iecblk(je,jb,2)) * &
                        p_patch%edges%primal_normal_cell(je,jb,2)%v1     + &
                        prm_diag%vmfl_s(iecidx(je,jb,2),iecblk(je,jb,2)) * &
                        p_patch%edges%primal_normal_cell(je,jb,2)%v2

           !Net stress at the edge
           flux_dn_e    = stress_c1n * p_int%c_lin_e(je,1,jb) + &
                          stress_c2n * p_int%c_lin_e(je,2,jb)

           tot_tend(je,nlev,jb) = tot_tend(je,nlev,jb) + (flux_up_e - flux_dn_e) * &
                p_nh_metrics%inv_ddqz_z_full_e(je,nlev,jb) * inv_rhoe(je,nlev,jb)
         END DO
      END DO
!$OMP END DO
!$OMP END PARALLEL

    CASE(iimplicit)

      !vertical tendency- implicit solver
      !a*x(k-1)+b*x(k)+c*x(k+1)=rhs(k)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jkp1,je,i_startidx,i_endidx,dwdn,a,b,c,rhs,var_new,&
!$OMP            flux_dn_e,stress_c1n,stress_c2n)
      DO jb = i_startblk,i_endblk
        CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
                           i_startidx, i_endidx, rl_start, rl_end)
#ifdef __LOOP_EXCHANGE
       DO je = i_startidx, i_endidx
         DO jk = 2 , nlev-1
#else
         DO jk = 2, nlev-1
           DO je = i_startidx, i_endidx
#endif
             jkp1  = jk + 1

             a(je,jk)   = - visc_smag_ie(je,jk,jb) * p_nh_metrics%inv_ddqz_z_full_e(je,jk,jb) * &
                               p_nh_metrics%inv_ddqz_z_half_e(je,jk,jb) * inv_rhoe(je,jk,jb)

             c(je,jk)   = - visc_smag_ie(je,jkp1,jb) * p_nh_metrics%inv_ddqz_z_full_e(je,jk,jb) * &
                               p_nh_metrics%inv_ddqz_z_half_e(je,jkp1,jb) * inv_rhoe(je,jk,jb)

             b(je,jk)   =  inv_dt - a(je,jk) - c(je,jk)

             !term due to dwdn- goes to RHS
             dwdn          = ( visc_smag_ie(je,jk,jb) * p_patch%edges%inv_dual_edge_length(je,jb) * &
                              (p_nh_prog%w(iecidx(je,jb,2),jk,iecblk(je,jb,2))  -  &
                               p_nh_prog%w(iecidx(je,jb,1),jk,iecblk(je,jb,1))) -  &
                               visc_smag_ie(je,jkp1,jb) * p_patch%edges%inv_dual_edge_length(je,jb)* &
                              (p_nh_prog%w(iecidx(je,jb,2),jkp1,iecblk(je,jb,2))  -  &
                               p_nh_prog%w(iecidx(je,jb,1),jkp1,iecblk(je,jb,1)))    &
                             ) * p_nh_metrics%inv_ddqz_z_full_e(je,jk,jb) * inv_rhoe(je,jk,jb)

             rhs(je,jk) =  p_nh_prog%vn(je,jk,jb) * inv_dt + dwdn
           END DO
         END DO

         !TOP Boundary treatment
         ! jk = 1
         !
         DO je = i_startidx, i_endidx
           c(je,1)   = - visc_smag_ie(je,2,jb) * p_nh_metrics%inv_ddqz_z_full_e(je,1,jb) * &
                            p_nh_metrics%inv_ddqz_z_half_e(je,2,jb) * inv_rhoe(je,1,jb)

           b(je,1)   =  inv_dt - c(je,1)

           !term due to dwdn- goes to RHS
           dwdn          =  -visc_smag_ie(je,2,jb) * p_patch%edges%inv_dual_edge_length(je,jb) * &
                            (p_nh_prog%w(iecidx(je,jb,2),2,iecblk(je,jb,2))  -  &
                             p_nh_prog%w(iecidx(je,jb,1),2,iecblk(je,jb,1))) *  &
                             p_nh_metrics%inv_ddqz_z_full_e(je,1,jb) * inv_rhoe(je,1,jb)

           rhs(je,1) =  p_nh_prog%vn(je,1,jb) * inv_dt + dwdn
         END DO

         !SFC Boundary treatment
         !jk = nlev
         !
         DO je = i_startidx, i_endidx
           a(je,nlev)  = - visc_smag_ie(je,nlev,jb) * p_nh_metrics%inv_ddqz_z_full_e(je,nlev,jb) * &
                              p_nh_metrics%inv_ddqz_z_half_e(je,nlev,jb) * inv_rhoe(je,nlev,jb)

           b(je,nlev)  = inv_dt - a(je,nlev)

           !term due to dwdn- goes to RHS
           dwdn          =  visc_smag_ie(je,nlev,jb) * p_patch%edges%inv_dual_edge_length(je,jb) * &
                            (p_nh_prog%w(iecidx(je,jb,2),nlev,iecblk(je,jb,2))  -  &
                             p_nh_prog%w(iecidx(je,jb,1),nlev,iecblk(je,jb,1))) *  &
                             p_nh_metrics%inv_ddqz_z_full_e(je,nlev,jb) * inv_rhoe(je,nlev,jb)

           !Get net shear stress in the direction of vn at surface

           !shear stress in normal direction from cell 1
           stress_c1n = prm_diag%umfl_s(iecidx(je,jb,1),iecblk(je,jb,1)) * &
                        p_patch%edges%primal_normal_cell(je,jb,1)%v1     + &
                        prm_diag%vmfl_s(iecidx(je,jb,1),iecblk(je,jb,1)) * &
                        p_patch%edges%primal_normal_cell(je,jb,1)%v2

           !shear stress in normal direction from cell 2
           stress_c2n = prm_diag%umfl_s(iecidx(je,jb,2),iecblk(je,jb,2)) * &
                        p_patch%edges%primal_normal_cell(je,jb,2)%v1     + &
                        prm_diag%vmfl_s(iecidx(je,jb,2),iecblk(je,jb,2)) * &
                        p_patch%edges%primal_normal_cell(je,jb,2)%v2

           !Net stress at the edge
           flux_dn_e    = stress_c1n * p_int%c_lin_e(je,1,jb) + stress_c2n * p_int%c_lin_e(je,2,jb)

           rhs(je,nlev) = p_nh_prog%vn(je,nlev,jb) * inv_dt + dwdn - flux_dn_e * &
                          p_nh_metrics%inv_ddqz_z_full_e(je,nlev,jb) * inv_rhoe(je,nlev,jb)
         END DO

         !CALL TDMA
         DO je = i_startidx, i_endidx
           CALL tdma_solver(a(je,:),b(je,:),c(je,:),rhs(je,:),nlev,var_new(:))
           tot_tend(je,:,jb) = tot_tend(je,:,jb) + (var_new(:) - p_nh_prog%vn(je,:,jb)) * inv_dt
         END DO

      END DO!block loop
!$OMP END DO
!$OMP END PARALLEL

    END SELECT !vert_scheme


   ! 4) Update vn: it makes more sense to first apply diffusion on vn
   ! and then get ddt_u/v than to interpolating tot_tend directly to
   ! ddt_u/v. Proof: during the test phase it was found that the latter slowed
   ! down the computation by 15%, although the results look same.

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,je,i_startidx,i_endidx)
    DO jb = i_startblk,i_endblk
      CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, rl_start, rl_end)
#ifdef __LOOP_EXCHANGE
       DO je = i_startidx, i_endidx
         DO jk = 1 , nlev
#else
      DO jk = 1 , nlev
       DO je = i_startidx, i_endidx
#endif
         vn_new(je,jk,jb) = p_nh_prog%vn(je,jk,jb) + dt * tot_tend(je,jk,jb)
       END DO
      END DO
    END DO
!$OMP END DO
!$OMP END PARALLEL


    !5) Get turbulent tendency at cell center
    CALL sync_patch_array(SYNC_E, p_patch, vn_new)
    CALL rbf_vec_interpol_cell(vn_new, p_patch, p_int, unew, vnew, opt_rlend=min_rlcell_int)

    rl_start = grf_bdywidth_c+1
    rl_end   = min_rlcell_int

    i_startblk = p_patch%cells%start_blk(rl_start,1)
    i_endblk   = p_patch%cells%end_blk(rl_end,i_nchdom)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx)
    DO jb = i_startblk,i_endblk
      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, rl_start, rl_end)
#ifdef __LOOP_EXCHANGE
       DO jc = i_startidx, i_endidx
         DO jk = 1 , nlev
#else
      DO jk = 1 , nlev
       DO jc = i_startidx, i_endidx
#endif
         ddt_u(jc,jk,jb) = ( unew(jc,jk,jb) - p_nh_diag%u(jc,jk,jb) ) * inv_dt
         ddt_v(jc,jk,jb) = ( vnew(jc,jk,jb) - p_nh_diag%v(jc,jk,jb) ) * inv_dt
       END DO
      END DO
    END DO
!$OMP END DO
!$OMP END PARALLEL

!    CALL sync_patch_array(SYNC_E, p_patch, tot_tend)
!    CALL rbf_vec_interpol_cell(tot_tend, p_patch, p_int, ddt_u, ddt_v, opt_rlend=min_rlcell_int-1)

    !subgrid fluxes: using vn_new, unew, vnew to store sgs_fluxes

    IF(is_sampling_time)THEN

!$OMP PARALLEL PRIVATE(rl_start,rl_end,i_startblk,i_endblk)

      CALL init(unew(:,:,:))
      CALL init(vnew(:,:,:))
      CALL init(vn_new(:,:,:))
!$OMP BARRIER

      rl_start = grf_bdywidth_e+1
      rl_end   = min_rledge_int
      i_startblk = p_patch%edges%start_blk(rl_start,1)
      i_endblk   = p_patch%edges%end_blk(rl_end,i_nchdom)

!$OMP DO PRIVATE(jb,jk,je,i_startidx,i_endidx,stress_c1n,stress_c2n)
      DO jb = i_startblk,i_endblk
        CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
                           i_startidx, i_endidx, rl_start, rl_end)
#ifdef __LOOP_EXCHANGE
       DO je = i_startidx, i_endidx
         DO jk = 2 , nlev
#else
        DO jk = 2 , nlev
         DO je = i_startidx, i_endidx
#endif
           vn_new(je,jk-1,jb) = -visc_smag_ie(je,jk,jb) *                  &
                  ( (p_nh_prog%vn(je,jk-1,jb) - p_nh_prog%vn(je,jk,jb)) * &
                     p_nh_metrics%inv_ddqz_z_half_e(je,jk,jb) +           &
                    (p_nh_prog%w(iecidx(je,jb,2),jk,iecblk(je,jb,2)) -    &
                     p_nh_prog%w(iecidx(je,jb,1),jk,iecblk(je,jb,1))) *   &
                     p_patch%edges%inv_dual_edge_length(je,jb) )
         END DO
        END DO

        !-----------------------------------------------------------------
        ! jk = nlev+1
        !-----------------------------------------------------------------
        DO je = i_startidx, i_endidx

           !Get net shear stress in the direction of vn at surface

           !shear stress in normal direction from cell 1
           stress_c1n = prm_diag%umfl_s(iecidx(je,jb,1),iecblk(je,jb,1)) * &
                        p_patch%edges%primal_normal_cell(je,jb,1)%v1     + &
                        prm_diag%vmfl_s(iecidx(je,jb,1),iecblk(je,jb,1)) * &
                        p_patch%edges%primal_normal_cell(je,jb,1)%v2

           !shear stress in normal direction from cell 2
           stress_c2n = prm_diag%umfl_s(iecidx(je,jb,2),iecblk(je,jb,2)) * &
                        p_patch%edges%primal_normal_cell(je,jb,2)%v1     + &
                        prm_diag%vmfl_s(iecidx(je,jb,2),iecblk(je,jb,2)) * &
                        p_patch%edges%primal_normal_cell(je,jb,2)%v2

           !Net stress at the edge
           vn_new(je,nlev,jb)  = -stress_c1n * p_int%c_lin_e(je,1,jb) - &
                                  stress_c2n * p_int%c_lin_e(je,2,jb)
         END DO
      END DO
!$OMP END DO
!$OMP END PARALLEL

      !Get sgs flux at cell center
      CALL sync_patch_array(SYNC_E, p_patch, vn_new)
      CALL rbf_vec_interpol_cell(vn_new, p_patch, p_int, unew, vnew, &
                                 opt_rlend=min_rlcell_int)

      !u sgs flux
      CALL levels_horizontal_mean(unew, p_patch%cells%area, p_patch%cells%owned, outvar)
      prm_diag%turb_diag_1dvar(1,idx_sgs_u_flx) = 0._wp
      DO jk = 2 , nlevp1
        prm_diag%turb_diag_1dvar(jk,idx_sgs_u_flx) =  &
              prm_diag%turb_diag_1dvar(jk,idx_sgs_u_flx)+outvar(jk-1)
      END DO

      !v sgs flux
      CALL levels_horizontal_mean(vnew, p_patch%cells%area, p_patch%cells%owned, outvar)
      prm_diag%turb_diag_1dvar(1,idx_sgs_v_flx) = 0._wp
      DO jk = 2 , nlevp1
        prm_diag%turb_diag_1dvar(jk,idx_sgs_v_flx) =  &
              prm_diag%turb_diag_1dvar(jk,idx_sgs_v_flx)+outvar(jk-1)
      END DO

    END IF!is_sampling_time

  END SUBROUTINE diffuse_hori_velocity
  !-------------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------------
  !>
  !! diffuse_vert_velocity
  !!------------------------------------------------------------------------
  !! Calculate the SGS diffusion term for vertical velocity component
  !! - Uses the forward Euler time scheme in time split (sequential) manner adopted
  !!   in the NH version.
  !! - Option to switch on implicit scheme in vertical
  !! - only solves for jk=2 to nlev. The bottom and top boundaries are left untouched
  !!------------------------------------------------------------------------
  !! @par Revision History
  !! Initial release by Anurag Dipankar, MPI-M (2013-02-05)
  SUBROUTINE diffuse_vert_velocity(p_nh_prog, p_nh_diag, p_nh_metrics, &
                                   p_patch, p_int, visc_smag_ic, ddt_w, dt)

    TYPE(t_nh_prog),   INTENT(inout)     :: p_nh_prog    !< single nh prognostic state
    TYPE(t_nh_diag),   INTENT(in)        :: p_nh_diag    !< single nh diagnostic state
    TYPE(t_nh_metrics),INTENT(in),TARGET :: p_nh_metrics !< single nh metric state
    TYPE(t_patch), TARGET, INTENT(in)    :: p_patch      !< single patch
    TYPE(t_int_state), INTENT(in),TARGET :: p_int        !< single interpolation state
    REAL(wp),          INTENT(in)        :: visc_smag_ic(:,:,:)
    REAL(wp),   TARGET, INTENT(inout)    :: ddt_w(:,:,:) !< w tendency
    REAL(wp),          INTENT(in)        :: dt

    REAL(wp) :: flux_up_c, flux_dn_c, dvn1, dvn2, dvt1, dvt2, flux_up_v, flux_dn_v
    REAL(wp) :: vt_e(nproma,p_patch%nlev,p_patch%nblks_e), inv_dt
    REAL(wp), DIMENSION(nproma,p_patch%nlev) :: a, b, c, rhs
    REAL(wp)                                 :: var_new(p_patch%nlev)

    !interface level variables but only nlev quantities are needed
    REAL(wp) :: hor_tend(nproma,p_patch%nlev,p_patch%nblks_e)
    REAL(wp), POINTER :: tot_tend(:,:,:)
    REAL(wp) :: inv_rho_ic(nproma,p_patch%nlev,p_patch%nblks_c)!not necessary to allocate for nlev+1

    INTEGER,  DIMENSION(:,:,:), POINTER :: ividx, ivblk, iecidx, iecblk, ieidx, ieblk
    INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx, i_nchdom
    INTEGER :: rl_start, rl_end, jg
    INTEGER :: jk, jkp1, jb, je, jkm1, jc, jcn, jbn, jvn
    INTEGER  :: nlev

    IF (msg_level >= 18) &
         CALL message(TRIM(inmodule), 'diffuse_vert_velocity')

    !patch id
    jg = p_patch%id

    ! number of vertical levels
    nlev     = p_patch%nlev
    i_nchdom = MAX(1,p_patch%n_childdom)

    inv_dt  = 1._wp / dt

    ividx  => p_patch%edges%vertex_idx
    ivblk  => p_patch%edges%vertex_blk

    iecidx => p_patch%edges%cell_idx
    iecblk => p_patch%edges%cell_blk

    ieidx => p_patch%cells%edge_idx
    ieblk => p_patch%cells%edge_blk

    tot_tend => ddt_w

    !Some initializations
    a = 0._wp; c = 0._wp

    IF(p_test_run)THEN
      hor_tend(:,:,:) = 0._wp
    END IF

    CALL rbf_vec_interpol_edge(p_nh_prog%vn, p_patch, p_int, vt_e, opt_rlend=min_rledge_int-1)

    !Calculate rho at interface for vertical diffusion
    rl_start = grf_bdywidth_c+1
    rl_end   = min_rlcell_int

    i_startblk = p_patch%cells%start_blk(rl_start,1)
    i_endblk   = p_patch%cells%end_blk(rl_end,i_nchdom)

!$OMP PARALLEL
!$OMP DO PRIVATE(jc,jb,jk,i_startidx,i_endidx)
    DO jb = i_startblk,i_endblk
       CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                          i_startidx, i_endidx, rl_start, rl_end)
#ifdef __LOOP_EXCHANGE
       DO jc = i_startidx, i_endidx
         DO jk = 2 , nlev
#else
       DO jk = 2 , nlev
         DO jc = i_startidx, i_endidx
#endif
           inv_rho_ic(jc,jk,jb) = 1._wp / rho_ic(jc,jk,jb)
         END DO
       END DO
    END DO
!$OMP END DO
!$OMP END PARALLEL


    ! 1) First get the horizontal tendencies at the half level edges

    rl_start = grf_bdywidth_e
    rl_end   = min_rledge_int-1

    i_startblk = p_patch%edges%start_blk(rl_start,1)
    i_endblk   = p_patch%edges%end_blk(rl_end,i_nchdom)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,je,i_startidx,i_endidx,jkm1,jcn,jbn,dvn1,dvn2,flux_up_c,flux_dn_c,&
!$OMP            jvn,dvt1,dvt2,flux_up_v,flux_dn_v)
    DO jb = i_startblk,i_endblk

      CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, rl_start, rl_end)
#ifdef __LOOP_EXCHANGE
       DO je = i_startidx, i_endidx
         DO jk = 2 , nlev
#else
      DO jk = 2 , nlev
       DO je = i_startidx, i_endidx
#endif
         jkm1    = jk-1

         !tendency in normal direction
         !flux = visc_c * D_31_c where D_31(=D_13) is calculated at half level
         !cell center

         jcn     = iecidx(je,jb,2)
         jbn     = iecblk(je,jb,2)

         dvn2  = p_nh_diag%u(jcn,jkm1,jbn)*p_patch%edges%primal_normal_cell(je,jb,2)%v1 + &
                 p_nh_diag%v(jcn,jkm1,jbn)*p_patch%edges%primal_normal_cell(je,jb,2)%v2 - &
                 p_nh_diag%u(jcn,jk,jbn)*p_patch%edges%primal_normal_cell(je,jb,2)%v1   - &
                 p_nh_diag%v(jcn,jk,jbn)*p_patch%edges%primal_normal_cell(je,jb,2)%v2

         flux_up_c = visc_smag_ic(jcn,jk,jbn) * (                    &
                     dvn2*p_nh_metrics%inv_ddqz_z_half(jcn,jk,jbn) + &
                     (w_vert(ividx(je,jb,4),jk,ivblk(je,jb,4))-w_ie(je,jk,jb)) * &
                     p_patch%edges%inv_vert_vert_length(je,jb)*2.0_wp )

         jcn     = iecidx(je,jb,1)
         jbn     = iecblk(je,jb,1)

         dvn1  = p_nh_diag%u(jcn,jkm1,jbn)*p_patch%edges%primal_normal_cell(je,jb,1)%v1 + &
                 p_nh_diag%v(jcn,jkm1,jbn)*p_patch%edges%primal_normal_cell(je,jb,1)%v2 - &
                 p_nh_diag%u(jcn,jk,jbn)*p_patch%edges%primal_normal_cell(je,jb,1)%v1   - &
                 p_nh_diag%v(jcn,jk,jbn)*p_patch%edges%primal_normal_cell(je,jb,1)%v2


         flux_dn_c = visc_smag_ic(jcn,jk,jbn) * (                    &
                     dvn1*p_nh_metrics%inv_ddqz_z_half(jcn,jk,jbn) + &
                     (w_ie(je,jk,jb)-w_vert(ividx(je,jb,3),jk,ivblk(je,jb,3))) *  &
                     p_patch%edges%inv_vert_vert_length(je,jb)*2.0_wp )


         !tendency in tangential direction
         !flux = visc_v * D_32_v where D_32(=D_23) is calculated at half level
         !between vertex and edge center

         jvn     = ividx(je,jb,2)
         jbn     = ivblk(je,jb,2)

         dvt2    = ( u_vert(jvn,jkm1,jbn)*p_patch%edges%dual_normal_vert(je,jb,2)%v1 + &
                     v_vert(jvn,jkm1,jbn)*p_patch%edges%dual_normal_vert(je,jb,2)%v2 + &
                     vt_e(je,jkm1,jb) ) * 0.5_wp           - &
                   ( u_vert(jvn,jk,jbn)*p_patch%edges%dual_normal_vert(je,jb,2)%v1 + &
                     v_vert(jvn,jk,jbn)*p_patch%edges%dual_normal_vert(je,jb,2)%v2 + &
                     vt_e(je,jk,jb) ) * 0.5_wp

         flux_up_v = visc_smag_iv(jvn,jk,jbn) * ( &
                     dvt2*p_nh_metrics%inv_ddqz_z_half_v(jvn,jk,jbn) + &
                     p_patch%edges%tangent_orientation(je,jb)*(w_vert(jvn,jk,jbn)-w_ie(je,jk,jb)) / &
                     p_patch%edges%edge_cell_length(je,jb,2) )


         jvn     = ividx(je,jb,1)
         jbn     = ivblk(je,jb,1)

         dvt1    = ( u_vert(jvn,jkm1,jbn)*p_patch%edges%dual_normal_vert(je,jb,1)%v1 + &
                     v_vert(jvn,jkm1,jbn)*p_patch%edges%dual_normal_vert(je,jb,1)%v2 + &
                     vt_e(je,jkm1,jb) ) * 0.5_wp           - &
                   ( u_vert(jvn,jk,jbn)*p_patch%edges%dual_normal_vert(je,jb,1)%v1 + &
                     v_vert(jvn,jk,jbn)*p_patch%edges%dual_normal_vert(je,jb,1)%v2 + &
                     vt_e(je,jk,jb) ) * 0.5_wp

         flux_dn_v = visc_smag_iv(jvn,jk,jbn) * ( &
                     dvt1*p_nh_metrics%inv_ddqz_z_half_v(jvn,jk,jbn) + &
                     p_patch%edges%tangent_orientation(je,jb)*(w_ie(je,jk,jb)-w_vert(jvn,jk,jbn)) / &
                     p_patch%edges%edge_cell_length(je,jb,1) )

         hor_tend(je,jk,jb) = (flux_up_c - flux_dn_c) * p_patch%edges%inv_dual_edge_length(je,jb) + &
                               p_patch%edges%tangent_orientation(je,jb) * (flux_up_v - flux_dn_v)  * &
                               p_patch%edges%inv_primal_edge_length(je,jb) * 2._wp

       END DO
      END DO
    END DO
!$OMP END DO
!$OMP END PARALLEL

    !CALL sync_patch_array(SYNC_E, p_patch, hor_tend)

    !Interpolate horizontal tendencies to w point: except top and bottom boundaries
    !w==0 at these boundaries

    rl_start = grf_bdywidth_c+1
    rl_end   = min_rlcell_int

    i_startblk = p_patch%cells%start_blk(rl_start,1)
    i_endblk   = p_patch%cells%end_blk(rl_end,i_nchdom)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx)
    DO jb = i_startblk,i_endblk
       CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                          i_startidx, i_endidx, rl_start, rl_end)
#ifdef __LOOP_EXCHANGE
       DO jc = i_startidx, i_endidx
         DO jk = 2 , nlev
#else
       DO jk = 2 , nlev
         DO jc = i_startidx, i_endidx
#endif
           tot_tend(jc,jk,jb) =                  inv_rho_ic(jc,jk,jb)                 * &
               ( hor_tend(ieidx(jc,jb,1),jk,ieblk(jc,jb,1))*p_int%e_bln_c_s(jc,1,jb)  + &
                 hor_tend(ieidx(jc,jb,2),jk,ieblk(jc,jb,2))*p_int%e_bln_c_s(jc,2,jb)  + &
                 hor_tend(ieidx(jc,jb,3),jk,ieblk(jc,jb,3))*p_int%e_bln_c_s(jc,3,jb) )
         END DO
       END DO
    END DO
!$OMP END DO
!$OMP END PARALLEL


    !
    ! 2) Vertical tendency: evaluated at w point
    !

    SELECT CASE(les_config(jg)%vert_scheme_type)

    CASE(iexplicit)
    !Vertical tendency - explicit solver

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx,flux_up_c,flux_dn_c, &
!$OMP            jkm1)
      DO jb = i_startblk,i_endblk

        CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                           i_startidx, i_endidx, rl_start, rl_end)
#ifdef __LOOP_EXCHANGE
       DO jc = i_startidx, i_endidx
         DO jk = 2 , nlev
#else
        DO jk = 2 , nlev
         DO jc = i_startidx, i_endidx
#endif
           jkm1 = jk - 1

           flux_up_c = visc_smag_c(jc,jkm1,jb) * ( (p_nh_prog%w(jc,jkm1,jb)-p_nh_prog%w(jc,jk,jb))* &
                       p_nh_metrics%inv_ddqz_z_full(jc,jkm1,jb) - z_1by3*DIV_c(jc,jkm1,jb) )

           flux_dn_c = visc_smag_c(jc,jk,jb) * ( (p_nh_prog%w(jc,jk,jb)-p_nh_prog%w(jc,jk+1,jb)) * &
                       p_nh_metrics%inv_ddqz_z_full(jc,jk,jb) - z_1by3*DIV_c(jc,jk,jb) )

           tot_tend(jc,jk,jb) = tot_tend(jc,jk,jb) + 2._wp * (flux_up_c - flux_dn_c) *  &
                                p_nh_metrics%inv_ddqz_z_half(jc,jk,jb) * inv_rho_ic(jc,jk,jb)

         END DO
        END DO
      END DO
!$OMP END DO
!$OMP END PARALLEL

    CASE(iimplicit)
    !vertical tendency- implicit solver
    !a*x(k-1)+b*x(k)+c*x(k+1)=rhs(k)

!$OMP PARALLEL
!$OMP DO PRIVATE(jc,jb,jk,jkm1,i_startidx,i_endidx,a,b,c,rhs,var_new)
      DO jb = i_startblk,i_endblk
         CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                            i_startidx, i_endidx, rl_start, rl_end)
#ifdef __LOOP_EXCHANGE
       DO jc = i_startidx, i_endidx
         DO jk = 3 , nlev-1
#else
         DO jk = 3,nlev-1
          DO jc = i_startidx, i_endidx
#endif
             jkm1 = jk - 1

             a(jc,jk)   = - 2._wp * visc_smag_c(jc,jkm1,jb) * p_nh_metrics%inv_ddqz_z_full(jc,jkm1,jb) * &
                               p_nh_metrics%inv_ddqz_z_half(jc,jk,jb) * inv_rho_ic(jc,jk,jb)

             c(jc,jk)   = - 2._wp * visc_smag_c(jc,jk,jb) * p_nh_metrics%inv_ddqz_z_full(jc,jk,jb) * &
                               p_nh_metrics%inv_ddqz_z_half(jc,jk,jb) * inv_rho_ic(jc,jk,jb)

             b(jc,jk)   =  inv_dt - a(jc,jk) - c(jc,jk)

             rhs(jc,jk) =   p_nh_prog%w(jc,jk,jb) * inv_dt + &
                               2._wp*( visc_smag_c(jc,jk,jb)*z_1by3*DIV_c(jc,jk,jb) -       &
                                       visc_smag_c(jc,jkm1,jb)*z_1by3*DIV_c(jc,jkm1,jb) ) * &
                                       p_nh_metrics%inv_ddqz_z_half(jc,jk,jb) * inv_rho_ic(jc,jk,jb)
          END DO
         END DO

          !TOP Boundary treatment:w==0
          !jk = 2
          !
          DO jc = i_startidx, i_endidx
            c(jc,2)   = - 2._wp * visc_smag_c(jc,2,jb) * p_nh_metrics%inv_ddqz_z_full(jc,2,jb) * &
                             p_nh_metrics%inv_ddqz_z_half(jc,2,jb) * inv_rho_ic(jc,2,jb)

            b(jc,2)   = inv_dt - c(jc,2) + 2._wp * visc_smag_c(jc,1,jb) * &
                           p_nh_metrics%inv_ddqz_z_full(jc,1,jb) *              &
                           p_nh_metrics%inv_ddqz_z_half(jc,2,jb) * inv_rho_ic(jc,2,jb)

            rhs(jc,2) = p_nh_prog%w(jc,2,jb) * inv_dt +  &
                           2._wp * ( visc_smag_c(jc,2,jb)*z_1by3*DIV_c(jc,2,jb) -   &
                                     visc_smag_c(jc,1,jb)*z_1by3*DIV_c(jc,1,jb) ) * &
                            p_nh_metrics%inv_ddqz_z_half(jc,2,jb) * inv_rho_ic(jc,2,jb)
          END DO

          !SFC Boundary treatment:w==0
          !jk = nlev
          !
          DO jc = i_startidx, i_endidx
            a(jc,nlev)  = - visc_smag_c(jc,nlev-1,jb) * p_nh_metrics%inv_ddqz_z_full(jc,nlev-1,jb) * &
                               p_nh_metrics%inv_ddqz_z_half(jc,nlev,jb) * 2._wp * inv_rho_ic(jc,nlev,jb)

            b(jc,nlev)  = inv_dt - a(jc,nlev) + 2._wp * visc_smag_c(jc,nlev,jb) * &
                             p_nh_metrics%inv_ddqz_z_full(jc,nlev,jb) *                 &
                             p_nh_metrics%inv_ddqz_z_half(jc,nlev,jb) * inv_rho_ic(jc,nlev,jb)

            rhs(jc,nlev)= p_nh_prog%w(jc,nlev,jb) * inv_dt +  &
                             2._wp * ( visc_smag_c(jc,nlev,jb)*z_1by3*DIV_c(jc,nlev,jb) -       &
                                       visc_smag_c(jc,nlev-1,jb)*z_1by3*DIV_c(jc,nlev-1,jb) ) * &
                                       p_nh_metrics%inv_ddqz_z_half(jc,nlev,jb) * inv_rho_ic(jc,nlev,jb)
          END DO

          !CALL TDMA
          DO jc = i_startidx, i_endidx
            CALL tdma_solver(a(jc,2:nlev),b(jc,2:nlev),c(jc,2:nlev), &
                             rhs(jc,2:nlev),nlev-1,var_new(2:nlev))

            tot_tend(jc,2:nlev,jb) = tot_tend(jc,2:nlev,jb) +  &
                                     ( var_new(2:nlev)-p_nh_prog%w(jc,2:nlev,jb) ) * inv_dt
         END DO

      END DO !block
!$OMP END DO
!$OMP END PARALLEL

    END SELECT !vert_scheme

    !Now update w

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx)
      DO jb = i_startblk,i_endblk
        CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                           i_startidx, i_endidx, rl_start, rl_end)
#ifdef __LOOP_EXCHANGE
       DO jc = i_startidx, i_endidx
         DO jk = 2 , nlev
#else
        DO jk = 2 , nlev
         DO jc = i_startidx, i_endidx
#endif
           p_nh_prog%w(jc,jk,jb) = p_nh_prog%w(jc,jk,jb) + dt * tot_tend(jc,jk,jb)
         END DO
        END DO
      END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL


   CALL sync_patch_array(SYNC_C, p_patch, p_nh_prog%w)

  END SUBROUTINE diffuse_vert_velocity
  !-------------------------------------------------------------------------------------


  !-------------------------------------------------------------------------------------
  !>
  !! diffuse_scalar
  !!------------------------------------------------------------------------
  !! Calculate the SGS diffusion term for cell based scalars
  !! - Uses the forward Euler time scheme in time split (sequential) manner adopted
  !!   in the NH version.
  !! - Option to switch on implicit scheme in vertical
  !!------------------------------------------------------------------------
  !! @par Revision History
  !! Initial release by Anurag Dipankar, MPI-M (2013-02-05)
  SUBROUTINE diffuse_scalar(var, p_nh_metrics, p_patch, p_int, p_nh_diag, tot_tend,  &
                            exner, prm_diag, rho, dt, scalar_name)

    REAL(wp),          INTENT(in)        :: var(:,:,:)   !input scalar
    TYPE(t_nh_metrics),INTENT(in),TARGET :: p_nh_metrics !< single nh metric state
    TYPE(t_patch),     INTENT(in),TARGET :: p_patch      !< single patch
    TYPE(t_int_state), INTENT(in),TARGET :: p_int        !< single interpolation state
    TYPE(t_nh_diag),   INTENT(in)        :: p_nh_diag    !< single nh diagnostic state
    REAL(wp),        INTENT(inout),TARGET:: tot_tend(:,:,:)!<total tendency
    REAL(wp),          INTENT(in)        :: exner(:,:,:)   !
    REAL(wp),          INTENT(in)        :: rho(:,:,:)     !density at cell center
    TYPE(t_nwp_phy_diag),  INTENT(inout) :: prm_diag       !< atm phys vars
    REAL(wp),          INTENT(in)        :: dt
    CHARACTER(*), INTENT(in)             :: scalar_name

    !Local variables
    REAL(wp) :: flux_up, flux_dn, inv_dt
    REAL(wp) :: nabla2_e(nproma,p_patch%nlev,p_patch%nblks_e)
    REAL(wp) :: fac(nproma,p_patch%nlev,p_patch%nblks_c)
    REAL(wp) :: exner_ic(nproma,p_patch%nlev+1,p_patch%nblks_c)
    REAL(wp) :: sgs_flux(nproma,p_patch%nlev+1,p_patch%nblks_c)
    REAL(wp) :: exner_me(nproma,p_patch%nlev,p_patch%nblks_e)
    REAL(wp) :: sflux(nproma,p_patch%nblks_c), outvar(p_patch%nlev+1)

    REAL(wp), DIMENSION(nproma,p_patch%nlev) :: a, b, c, rhs
    REAL(wp)                                 :: var_new(p_patch%nlev)
    REAL(wp), POINTER :: diff_smag_ic(:,:,:)

    INTEGER,  DIMENSION(:,:,:), POINTER :: iecidx, iecblk, ieidx, ieblk
    INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx, i_nchdom
    INTEGER :: rl_start, rl_end
    INTEGER :: jk, jb, je, jc, jg
    INTEGER  :: nlev, nlevp1

    IF (msg_level >= 18) &
         CALL message(TRIM(inmodule), 'diffuse_scalar')

    !patch id
    jg = p_patch%id

    ! number of vertical levels
    nlev = p_patch%nlev
    nlevp1 = nlev+1
    i_nchdom   = MAX(1,p_patch%n_childdom)
    inv_dt  = 1._wp / dt

    iecidx => p_patch%edges%cell_idx
    iecblk => p_patch%edges%cell_blk

    ieidx => p_patch%cells%edge_idx
    ieblk => p_patch%cells%edge_blk

    diff_smag_ic => prm_diag%tkvh

    !multiply by exner to convert from theta tend to temp tend
    !assuming that exner perturbation are small compared to temp

    !1) First set exner local vars to 1 for other scalars
    !   Soon get different routines for different scalars

!$OMP PARALLEL
    CALL init(exner_me(:,:,:), 1._wp)
    CALL init(exner_ic(:,:,:), 1._wp)
    CALL init(a(:,:))
    CALL init(c(:,:))
!$OMP END PARALLEL

    !2) Calculate exner at edge for horizontal diffusion
     IF(TRIM(scalar_name)=='theta') &
        CALL cells2edges_scalar(exner, p_patch, p_int%c_lin_e, exner_me, opt_rlend=min_rledge_int-2)

    !3) Calculate exner at interface for vertical diffusion

!$OMP PARALLEL PRIVATE(rl_start, rl_end, i_startblk, i_endblk)
    rl_start = grf_bdywidth_c+1
    rl_end   = min_rlcell_int

    i_startblk = p_patch%cells%start_blk(rl_start,1)
    i_endblk   = p_patch%cells%end_blk(rl_end,i_nchdom)

    IF(TRIM(scalar_name)=='theta')THEN

!$OMP DO PRIVATE(jc,jb,jk,i_startidx,i_endidx)
      DO jb = i_startblk,i_endblk
         CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                            i_startidx, i_endidx, rl_start, rl_end)
#ifdef __LOOP_EXCHANGE
       DO jc = i_startidx, i_endidx
         DO jk = 2 , nlev
#else
         DO jk = 2, nlev
           DO jc = i_startidx, i_endidx
#endif
             exner_ic(jc,jk,jb) = ( p_nh_metrics%wgtfac_c(jc,jk,jb)*exner(jc,jk,jb) + &
                                  (1._wp-p_nh_metrics%wgtfac_c(jc,jk,jb))*exner(jc,jk-1,jb) )
           END DO
         END DO
         DO jc = i_startidx, i_endidx
            exner_ic(jc,nlev+1,jb) =                                &
              p_nh_metrics%wgtfacq_c(jc,1,jb)*exner(jc,nlev  ,jb) + &
              p_nh_metrics%wgtfacq_c(jc,2,jb)*exner(jc,nlev-1,jb) + &
              p_nh_metrics%wgtfacq_c(jc,3,jb)*exner(jc,nlev-2,jb)
         ENDDO
      END DO
!$OMP END DO

    END IF

    !Special boundary treatment for different scalars

    IF(TRIM(scalar_name)=='theta')THEN
!$OMP DO PRIVATE(jc,jb,jk,i_startidx,i_endidx)
        DO jb = i_startblk,i_endblk
          CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                             i_startidx, i_endidx, rl_start, rl_end)
#ifdef __LOOP_EXCHANGE
          DO jc = i_startidx, i_endidx
            DO jk = 1 , nlev
#else
          DO jk = 1, nlev
            DO jc = i_startidx, i_endidx
#endif
              fac(jc,jk,jb) = cpd * rcvd / rho(jc,jk,jb)
            END DO
          END DO
          DO jc = i_startidx, i_endidx
             sflux(jc,jb) = prm_diag%shfl_s(jc,jb) * rcpd
          END DO
        END DO
!$OMP END DO
    ELSEIF(TRIM(scalar_name)=='qv')THEN
!$OMP DO PRIVATE(jc,jb,jk,i_startidx,i_endidx)
        DO jb = i_startblk,i_endblk
          CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                             i_startidx, i_endidx, rl_start, rl_end)
#ifdef __LOOP_EXCHANGE
          DO jc = i_startidx, i_endidx
            DO jk = 1 , nlev
#else
          DO jk = 1, nlev
            DO jc = i_startidx, i_endidx
#endif
              fac(jc,jk,jb) = 1._wp / rho(jc,jk,jb)
            END DO
          END DO
          DO jc = i_startidx, i_endidx
             sflux(jc,jb) = prm_diag%lhfl_s(jc,jb) / alv
          END DO
        END DO
!$OMP END DO
    ELSEIF(TRIM(scalar_name)=='qc')THEN
!$OMP DO PRIVATE(jc,jb,jk,i_startidx,i_endidx)
        DO jb = i_startblk,i_endblk
          CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                             i_startidx, i_endidx, rl_start, rl_end)
#ifdef __LOOP_EXCHANGE
          DO jc = i_startidx, i_endidx
            DO jk = 1 , nlev
#else
          DO jk = 1, nlev
            DO jc = i_startidx, i_endidx
#endif
              fac(jc,jk,jb) = 1._wp / rho(jc,jk,jb)
            END DO
          END DO
          DO jc = i_startidx, i_endidx
             sflux(jc,jb) = 0._wp
          END DO
        END DO
!$OMP END DO
    END IF

    ! use conservative discretization div(k*grad(var))-horizontal part
    ! following mo_nh_diffusion

    !include halo points and boundary points because these values will be
    !used in next loop
    rl_start = grf_bdywidth_e
    rl_end   = min_rledge_int-1

    i_startblk = p_patch%edges%start_blk(rl_start,1)
    i_endblk   = p_patch%edges%end_blk(rl_end,i_nchdom)

!$OMP DO PRIVATE(jk,je,jb,i_startidx,i_endidx)
    DO jb = i_startblk,i_endblk
       CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
                          i_startidx, i_endidx, rl_start, rl_end)
#ifdef __LOOP_EXCHANGE
       DO je = i_startidx, i_endidx
         DO jk = 1 , nlev
#else
       ! compute kh_smag_e * grad(var)
       DO jk = 1, nlev
         DO je = i_startidx, i_endidx
#endif
             nabla2_e(je,jk,jb) = 0.5_wp*(visc_smag_ie(je,jk,jb)+visc_smag_ie(je,jk+1,jb))* &
                                  les_config(jg)%rturb_prandtl * exner_me(je,jk,jb) * &
                                  p_patch%edges%inv_dual_edge_length(je,jb)*   &
                                 (var(iecidx(je,jb,2),jk,iecblk(je,jb,2)) -    &
                                  var(iecidx(je,jb,1),jk,iecblk(je,jb,1)))
         ENDDO
       ENDDO
    ENDDO
!$OMP END DO

    ! now compute the divergence of the quantity above
    rl_start = grf_bdywidth_c+1
    rl_end   = min_rlcell_int

    i_startblk = p_patch%cells%start_blk(rl_start,1)
    i_endblk   = p_patch%cells%end_blk(rl_end,i_nchdom)

!$OMP DO PRIVATE(jc,jb,jk,i_startidx,i_endidx)
    DO jb = i_startblk,i_endblk
       CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                          i_startidx, i_endidx, rl_start, rl_end)
#ifdef __LOOP_EXCHANGE
       DO jc = i_startidx, i_endidx
         DO jk = 1 , nlev
#else
       DO jk = 1, nlev
         DO jc = i_startidx, i_endidx
#endif
           !horizontal tendency
           tot_tend(jc,jk,jb)  =  fac(jc,jk,jb) * (                                 &
             nabla2_e(ieidx(jc,jb,1),jk,ieblk(jc,jb,1))*p_int%geofac_div(jc,1,jb) + &
             nabla2_e(ieidx(jc,jb,2),jk,ieblk(jc,jb,2))*p_int%geofac_div(jc,2,jb) + &
             nabla2_e(ieidx(jc,jb,3),jk,ieblk(jc,jb,3))*p_int%geofac_div(jc,3,jb) )
         ENDDO
       ENDDO
    ENDDO
!$OMP END DO
!$OMP END PARALLEL


    !---------------------------------------------------------------
    !Vertical diffusion
    !---------------------------------------------------------------

    rl_start = grf_bdywidth_c+1
    rl_end   = min_rlcell_int

    i_startblk = p_patch%cells%start_blk(rl_start,1)
    i_endblk   = p_patch%cells%end_blk(rl_end,i_nchdom)


    SELECT CASE(les_config(jg)%vert_scheme_type)

    CASE(iexplicit)
    !Vertical tendency - explicit solver

!$OMP PARALLEL
!$OMP DO PRIVATE(jc,jb,jk,i_startidx,i_endidx,flux_up,flux_dn)
      DO jb = i_startblk,i_endblk
         CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                            i_startidx, i_endidx, rl_start, rl_end)
#ifdef __LOOP_EXCHANGE
         DO jc = i_startidx, i_endidx
           DO jk = 2 , nlev-1
#else
         DO jk = 2, nlev-1
           DO jc = i_startidx, i_endidx
#endif
             flux_up = diff_smag_ic(jc,jk,jb) * (var(jc,jk-1,jb) - var(jc,jk,jb)) * &
                       p_nh_metrics%inv_ddqz_z_half(jc,jk,jb) * exner_ic(jc,jk,jb)

             flux_dn = diff_smag_ic(jc,jk+1,jb) * (var(jc,jk,jb) - var(jc,jk+1,jb)) * &
                       p_nh_metrics%inv_ddqz_z_half(jc,jk+1,jb) * exner_ic(jc,jk+1,jb)

             tot_tend(jc,jk,jb) = tot_tend(jc,jk,jb) + (flux_up - flux_dn) *  &
                             p_nh_metrics%inv_ddqz_z_full(jc,jk,jb) * fac(jc,jk,jb)

           ENDDO
         ENDDO

         !Boundary treatment

         !--------------------------------------------------------
         !jk = 1
         !--------------------------------------------------------
         DO jc = i_startidx, i_endidx
            flux_up = 0._wp

            flux_dn = diff_smag_ic(jc,2,jb) * (var(jc,1,jb) - var(jc,2,jb)) * &
                      p_nh_metrics%inv_ddqz_z_half(jc,2,jb) * exner_ic(jc,2,jb)

            tot_tend(jc,1,jb) = tot_tend(jc,1,jb) + (flux_up- flux_dn) * &
                       p_nh_metrics%inv_ddqz_z_full(jc,1,jb) * fac(jc,1,jb)
         ENDDO
         !--------------------------------------------------------
         !jk = nlev
         !--------------------------------------------------------
         DO jc = i_startidx, i_endidx
           flux_up = diff_smag_ic(jc,nlev,jb) * (var(jc,nlev-1,jb) - var(jc,nlev,jb)) * &
                     p_nh_metrics%inv_ddqz_z_half(jc,nlev,jb) * exner_ic(jc,nlev,jb)

           flux_dn  = -sflux(jc,jb) * exner_ic(jc,nlevp1,jb)

           tot_tend(jc,nlev,jb) = tot_tend(jc,nlev,jb) + (flux_up - flux_dn) * &
                   p_nh_metrics%inv_ddqz_z_full(jc,nlev,jb) * fac(jc,nlev,jb)
         ENDDO
      END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    CASE(iimplicit)

    !vertical tendency- implicit solver
    !a*x(k-1)+b*x(k)+c*x(k+1)=rhs(k)

!$OMP PARALLEL
!$OMP DO PRIVATE(jc,jb,jk,i_startidx,i_endidx,a,b,c,rhs,var_new)
      DO jb = i_startblk,i_endblk
         CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                            i_startidx, i_endidx, rl_start, rl_end)
#ifdef __LOOP_EXCHANGE
         DO jc = i_startidx, i_endidx
           DO jk = 2 , nlev-1
#else
         DO jk = 2, nlev-1
           DO jc = i_startidx, i_endidx
#endif
              a(jc,jk)  = - diff_smag_ic(jc,jk,jb) * p_nh_metrics%inv_ddqz_z_full(jc,jk,jb) * &
                   p_nh_metrics%inv_ddqz_z_half(jc,jk,jb) * exner_ic(jc,jk,jb) * fac(jc,jk,jb)

              c(jc,jk)  = - diff_smag_ic(jc,jk+1,jb) * p_nh_metrics%inv_ddqz_z_full(jc,jk,jb) * &
               p_nh_metrics%inv_ddqz_z_half(jc,jk+1,jb) * exner_ic(jc,jk+1,jb) * fac(jc,jk,jb)

              b(jc,jk)   =  inv_dt - a(jc,jk) - c(jc,jk)

              rhs(jc,jk) =  var(jc,jk,jb) * inv_dt
           END DO
         END DO

         !TOP Boundary treatment
         ! jk = 1
         !
         DO jc = i_startidx, i_endidx
           c(jc,1)   = -diff_smag_ic(jc,2,jb) * p_nh_metrics%inv_ddqz_z_full(jc,1,jb) * &
                      p_nh_metrics%inv_ddqz_z_half(jc,2,jb) * exner_ic(jc,2,jb) * fac(jc,1,jb)

           b(jc,1)   = inv_dt - c(jc,1)

           rhs(jc,1) = var(jc,1,jb) * inv_dt
         END DO

         !SFC Boundary treatment
         !jk = nlev
         !
         DO jc = i_startidx, i_endidx
           a(jc,nlev) = -diff_smag_ic(jc,nlev,jb) * p_nh_metrics%inv_ddqz_z_full(jc,nlev,jb) * &
               p_nh_metrics%inv_ddqz_z_half(jc,nlev,jb) * exner_ic(jc,nlev,jb) * fac(jc,nlev,jb)

           b(jc,nlev)  = inv_dt - a(jc,nlev)

           rhs(jc,nlev)= var(jc,nlev,jb) * inv_dt + sflux(jc,jb) * exner_ic(jc,nlevp1,jb) * &
                          p_nh_metrics%inv_ddqz_z_full(jc,nlev,jb) * fac(jc,nlev,jb)
         END DO

         !CALL TDMA
         DO jc = i_startidx, i_endidx
            CALL tdma_solver(a(jc,:),b(jc,:),c(jc,:),rhs(jc,:),nlev,var_new)
            tot_tend(jc,:,jb) = tot_tend(jc,:,jb) + ( var_new(:) - var(jc,:,jb) ) * inv_dt
         END DO

         END DO!block
!$OMP END DO
!$OMP END PARALLEL

    END SELECT !vert_scheme


    !subgrid fluxes for different scalars

    IF(is_sampling_time)THEN
!$OMP PARALLEL
!$OMP DO PRIVATE(jc,jb,jk,i_startidx,i_endidx)
      DO jb = i_startblk,i_endblk
         CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                            i_startidx, i_endidx, rl_start, rl_end)
#ifdef __LOOP_EXCHANGE
         DO jc = i_startidx, i_endidx
           DO jk = 2 , nlev
#else
         DO jk = 2, nlev
           DO jc = i_startidx, i_endidx
#endif
              sgs_flux(jc,jk,jb) = -diff_smag_ic(jc,jk,jb) * (var(jc,jk-1,jb) - var(jc,jk,jb)) * &
                                   p_nh_metrics%inv_ddqz_z_half(jc,jk,jb)
           END DO
         END DO
         DO jc = i_startidx, i_endidx
            sgs_flux(jc,1,jb) = 0._wp
            !extrapolate for surface density
            sgs_flux(jc,nlevp1,jb) = sflux(jc,jb)
         END DO
      END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

      IF(TRIM(scalar_name)=='theta')THEN

        CALL levels_horizontal_mean(sgs_flux, p_patch%cells%area, p_patch%cells%owned, &
                                    outvar)

        outvar = outvar * cpd
        prm_diag%turb_diag_1dvar(1:nlevp1,idx_sgs_th_flx) =  &
               prm_diag%turb_diag_1dvar(1:nlevp1,idx_sgs_th_flx)+outvar(1:nlevp1)

      ELSEIF(TRIM(scalar_name)=='qv')THEN

        CALL levels_horizontal_mean(sgs_flux, p_patch%cells%area, p_patch%cells%owned, &
                                     outvar)

        outvar = outvar * alv
        prm_diag%turb_diag_1dvar(1:nlevp1,idx_sgs_qv_flx) =  &
               prm_diag%turb_diag_1dvar(1:nlevp1,idx_sgs_qv_flx)+outvar(1:nlevp1)

      ELSEIF(TRIM(scalar_name)=='qc')THEN

        CALL levels_horizontal_mean(sgs_flux, p_patch%cells%area, p_patch%cells%owned, &
                                    outvar)

        outvar = outvar * alv
        prm_diag%turb_diag_1dvar(1:nlevp1,idx_sgs_qc_flx) =  &
               prm_diag%turb_diag_1dvar(1:nlevp1,idx_sgs_qc_flx)+outvar(1:nlevp1)

      END IF

    END IF !is_sampling_time


  END SUBROUTINE diffuse_scalar

!---------------------------------------------------------------------------------

END MODULE mo_sgs_turbulence
