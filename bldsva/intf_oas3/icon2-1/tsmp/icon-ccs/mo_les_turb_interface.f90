!OPTION! -cont -msg o
!! this command should fix the problem of copying arrays in a subroutine call
!>
!!AD: Interface between the les_phy_interface and various turbulence schems suitable
!!    for LES. At present it only has classical Smagorinsky scheme. 
!!
!! @author Kristina Froehlich, DWD, Offenbach (2010-01-25)
!!
!! @par Revision History
!! Initial Kristina Froehlich, DWD, Offenbach (2010-01-25)
!! Modified by Anurag Dipankar, MPIM (2013-07-01) 
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

MODULE mo_les_turb_interface

  USE mo_kind,                 ONLY: wp
  USE mo_exception,            ONLY: message
  USE mo_model_domain,         ONLY: t_patch
  USE mo_intp_data_strc,       ONLY: t_int_state
  USE mo_impl_constants,       ONLY: min_rlcell_int, ismag
  USE mo_impl_constants_grf,   ONLY: grf_bdywidth_c
  USE mo_loopindices,          ONLY: get_indices_c
  USE mo_nonhydro_types,       ONLY: t_nh_prog, t_nh_diag, t_nh_metrics
  USE mo_nwp_phy_types,        ONLY: t_nwp_phy_diag, t_nwp_phy_tend
  USE mo_nwp_lnd_types,        ONLY: t_lnd_prog, t_lnd_diag
  USE mo_run_config,           ONLY: msg_level, iqv, iqc
  USE mo_atm_phy_nwp_config,   ONLY: atm_phy_nwp_config
  USE mo_nonhydrostatic_config,ONLY: kstart_moist
  USE mo_sgs_turbulence,       ONLY: drive_subgrid_diffusion
  USE mo_sgs_turbmetric,       ONLY: drive_subgrid_diffusion_m
  USE mo_les_config,           ONLY: les_config

  IMPLICIT NONE

  PRIVATE

  PUBLIC  ::  les_turbulence

CONTAINS
  !!
  !!-------------------------------------------------------------------------
  !!
SUBROUTINE les_turbulence  ( tcall_turb_jg,                   & !>in
                          & p_patch,                          & !>in
                          & p_metrics,                        & !>in
                          & p_int,                            & !>in
                          & p_prog,                           & !>in
                          & p_prog_rcf,                       & !>inout
                          & p_diag ,                          & !>inout
                          & prm_diag, prm_nwp_tend,           & !>inout
                          & lnd_prog_now,                     & !>in 
                          & lnd_prog_new,                     & !>inout only for idealized LES
                          & lnd_diag                          ) !>in


  TYPE(t_patch),        TARGET,INTENT(in)   :: p_patch        !!<grid/patch info.
  TYPE(t_int_state),    INTENT(in),TARGET   :: p_int          !< single interpolation state
  TYPE(t_nh_metrics)          ,INTENT(inout)   :: p_metrics  ! SBr, CHa:in->inout
  TYPE(t_nh_prog),      TARGET,INTENT(inout):: p_prog          !<the prog vars
  TYPE(t_nh_prog),      TARGET,INTENT(inout):: p_prog_rcf      !<call freq
  TYPE(t_nh_diag),      TARGET,INTENT(inout):: p_diag          !<the diag vars
  TYPE(t_nwp_phy_diag),        INTENT(inout):: prm_diag        !< atm phys vars
  TYPE(t_nwp_phy_tend), TARGET,INTENT(inout):: prm_nwp_tend    !< atm tend vars
  TYPE(t_lnd_prog),            INTENT(in)   :: lnd_prog_now    !< prog vars for sfc
  TYPE(t_lnd_prog),            INTENT(inout):: lnd_prog_new    !< prog vars for sfc
  TYPE(t_lnd_diag),            INTENT(inout):: lnd_diag        !< diag vars for sfc
  REAL(wp),                    INTENT(in)   :: tcall_turb_jg   !< time interval for turbulence

  ! Local array bounds

  INTEGER :: rl_start, rl_end
  INTEGER :: i_startblk, i_endblk    !> blocks
  INTEGER :: i_startidx, i_endidx    !< slices
  INTEGER :: i_nchdom                !< domain index

  ! Local scalars:
  INTEGER :: jc,jk,jb,jg      !loop indices

  INTEGER  :: nlev            !< number of full levels

!--------------------------------------------------------------

  ! number of vertical levels
  nlev   = p_patch%nlev

  ! local variables related to the blocking
  i_nchdom  = MAX(1,p_patch%n_childdom)
  jg        = p_patch%id

  IF (msg_level >= 15) CALL message('mo_les_turb_interface:', 'turbulence')
  
  !For 3D turbulence the whole patch needs to be passed. Therefore, this call
  !is made outside the block loop next. However, the tendencies it calculates
  !is then used inside the block loop (see at the end) to update u,v,t,qv,qc
  IF ( atm_phy_nwp_config(jg)%inwp_turb==ismag )THEN

    ! if les metrics is choosen, drive the subgrid diffusion from mo_sgs_turbmetric
    IF (les_config(jg)%les_metric) THEN
      CALL drive_subgrid_diffusion_m(p_prog,       & !inout for w (it is updated inside)
                                     p_prog_rcf,   & !in
                                     p_diag,       & !inout
                                     p_metrics,    & !in
                                     p_patch,      & !in
                                     p_int,        & !in
                                     lnd_prog_now, & !in
                                     lnd_prog_new, & !inout only for idealized cases
                                     lnd_diag,     & !inout
                                     prm_diag,     & !inout
                                     prm_nwp_tend, & !inout
                                     tcall_turb_jg & !in
                                     )


    ELSE
      CALL drive_subgrid_diffusion(p_prog,       & !inout for w (it is updated inside)
                                   p_prog_rcf,   & !in
                                   p_diag,       & !inout
                                   p_metrics,    & !inout, SBr, CHa:in->inout
                                   p_patch,      & !in
                                   p_int,        & !in
                                   lnd_prog_now, & !in
                                   lnd_prog_new, & !inout only for idealized cases
                                   lnd_diag,     & !inout
                                   prm_diag,     & !inout
                                   prm_nwp_tend, & !inout
                                   tcall_turb_jg & !in
                                   )
    END IF

  END IF

  
  ! exclude boundary interpolation zone of nested domains
  rl_start = grf_bdywidth_c+1
  rl_end   = min_rlcell_int

  i_startblk = p_patch%cells%start_blk(rl_start,1)
  i_endblk   = p_patch%cells%end_blk(rl_end,i_nchdom)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,jk,i_startidx,i_endidx) ICON_OMP_GUIDED_SCHEDULE
  DO jb = i_startblk, i_endblk

    CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
      & i_startidx, i_endidx, rl_start, rl_end)

    ! Update wind speed, QV and temperature with turbulence tendencies
    ! Note: the update of wind speed is done here in order to pass u and v at the correct time level
    ! to turbtran and the convection scheme. However, the update of the prognostic variable vn
    ! is done at the end of the NWP interface by first interpolating the u/v tendencies to the 
    ! velocity points (in order to minimize interpolation errors) and then adding the tendencies
    ! to vn
    DO jk = 1, nlev
      DO jc = i_startidx, i_endidx
        p_prog_rcf%tracer(jc,jk,jb,iqv) =MAX(0._wp, p_prog_rcf%tracer(jc,jk,jb,iqv) &
             &           + tcall_turb_jg*prm_nwp_tend%ddt_tracer_turb(jc,jk,jb,iqv))
        p_diag%temp(jc,jk,jb) = p_diag%temp(jc,jk,jb)  &
         &  + tcall_turb_jg*prm_nwp_tend%ddt_temp_turb(jc,jk,jb)
        p_diag%u(jc,jk,jb) = p_diag%u(jc,jk,jb) + tcall_turb_jg*prm_nwp_tend%ddt_u_turb(jc,jk,jb)
        p_diag%v(jc,jk,jb) = p_diag%v(jc,jk,jb) + tcall_turb_jg*prm_nwp_tend%ddt_v_turb(jc,jk,jb)
      ENDDO
    ENDDO
    ! QC is updated only in that part of the model domain where moisture physics is active
    DO jk = kstart_moist(jg), nlev
      DO jc = i_startidx, i_endidx
        p_prog_rcf%tracer(jc,jk,jb,iqc) =MAX(0._wp, p_prog_rcf%tracer(jc,jk,jb,iqc) &
             &           + tcall_turb_jg*prm_nwp_tend%ddt_tracer_turb(jc,jk,jb,iqc))
      ENDDO
    ENDDO

  ENDDO ! jb
!$OMP END DO NOWAIT
!$OMP END PARALLEL

END SUBROUTINE les_turbulence

END MODULE mo_les_turb_interface
