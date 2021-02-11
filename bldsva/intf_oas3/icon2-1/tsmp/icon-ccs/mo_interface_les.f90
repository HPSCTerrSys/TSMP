!>
!!
!! AD:This is the interface for LES (large eddy simulation) physics. At present
!! it follows the same structure as nwp_phy_interface and uses the same physics
!! except turbulence However, in future we will have different radiation and
!! microphysics routines and at that point of time it will be independent
!! of nwp physics (mostly to get rid of nwp_phy config states)
!!
!! @par Revision History
!!  first implementation by Kristina Froehlich, DWD (2009-06-12)
!!  Modified for les physics by Anurag Dipankar, MPIM (2013-07-01)
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

MODULE mo_interface_les

  USE mo_nonhydro_state,     ONLY: p_nh_state ! SBr,CHa
  USE mtime,                 ONLY: datetime, timeDelta, newTimedelta,     &
    &                              deallocateTimedelta, getTimedeltaFromDatetime, &
    &                              getTotalMillisecondsTimedelta
  USE mo_time_config,        ONLY: time_config
  USE mo_kind,               ONLY: wp
  USE mo_timer
  USE mo_exception,          ONLY: message, message_text
  USE mo_impl_constants,     ONLY: itccov, itrad, itgscp,         &
    &                              itsatad, itturb, itsfc, itradheat, &
    &                              itfastphy, max_char_length,    &
    &                              min_rlcell_int, min_rledge_int, min_rlcell
  USE mo_impl_constants_grf, ONLY: grf_bdywidth_c, grf_bdywidth_e
  USE mo_loopindices,        ONLY: get_indices_c, get_indices_e
  USE mo_intp_rbf,           ONLY: rbf_vec_interpol_cell
  USE mo_model_domain,       ONLY: t_patch
  USE mo_intp_data_strc,     ONLY: t_int_state
  USE mo_nonhydro_types,     ONLY: t_nh_prog, t_nh_diag, t_nh_metrics
  USE mo_nonhydrostatic_config, ONLY: kstart_moist, lhdiff_rcf, ih_clch, ih_clcm
  USE mo_nwp_lnd_types,      ONLY: t_lnd_prog, t_wtr_prog, t_lnd_diag
  USE mo_ext_data_types,     ONLY: t_external_data
  USE mo_nwp_phy_types,      ONLY: t_nwp_phy_diag, t_nwp_phy_tend
  USE mo_parallel_config,    ONLY: nproma, p_test_run, use_icon_comm, use_physics_barrier
  USE mo_diffusion_config,   ONLY: diffusion_config
  USE mo_run_config,         ONLY: ntracer, iqv, iqc, iqi, iqr, iqs, iqm_max,   &
    &                              msg_level, ltimer, timers_level, nqtendphy,  &
    &                              ltransport, lart, iqni, iqnc
  USE mo_physical_constants, ONLY: rd, rd_o_cpd, vtmpc1, p0ref, cvd, cvv

  USE mo_nh_diagnose_pres_temp,ONLY: diagnose_pres_temp

  USE mo_atm_phy_nwp_config, ONLY: atm_phy_nwp_config
  USE mo_lnd_nwp_config,     ONLY: ntiles_total, ntiles_water
  USE mo_cover_koe,          ONLY: cover_koe
  USE mo_satad,              ONLY: satad_v_3D
  USE mo_radiation,          ONLY: radheat, pre_radiation_nwp
  USE mo_radiation_config,   ONLY: irad_aero
  USE mo_nwp_gscp_interface, ONLY: nwp_microphysics
  USE mo_les_turb_interface, ONLY: les_turbulence
  USE mo_nwp_sfc_interface,  ONLY: nwp_surface
  USE mo_nwp_rad_interface,  ONLY: nwp_radiation
  USE mo_sync,               ONLY: sync_patch_array, sync_patch_array_mult, SYNC_E, &
                                   SYNC_C, SYNC_C1, global_sum_array
  USE mo_mpi,                ONLY: my_process_is_mpi_all_parallel, work_mpi_barrier
  USE mo_nwp_diagnosis,      ONLY: nwp_statistics, nwp_diag_output_1, nwp_diag_output_2
  USE mo_icon_comm_lib,      ONLY: new_icon_comm_variable, &
     & icon_comm_sync_all, is_ready, until_sync
  USE mo_art_washout_interface,  ONLY:art_washout_interface
  USE mo_linked_list,         ONLY: t_var_list
  USE mo_ls_forcing_nml,      ONLY: is_ls_forcing
  USE mo_ls_forcing,          ONLY: apply_ls_forcing
  USE mo_advection_config,    ONLY: advection_config
  USE mo_nwp_turbtrans_interface, ONLY: nwp_turbtrans
  USE mo_turbulent_diagnostic, ONLY: calculate_turbulent_diagnostics, &
                                     write_vertical_profiles, write_time_series, &
                                     avg_interval_step, sampl_freq_step,  &
                                     is_sampling_time, is_writing_time
  USE mo_les_utilities,       ONLY: init_vertical_grid_for_les
  USE mo_fortran_tools,       ONLY: copy

  IMPLICIT NONE

  PRIVATE


  REAL(wp), PARAMETER :: rd_o_p0ref = rd / p0ref
  REAL(wp), PARAMETER :: cpd_o_rd = 1._wp / rd_o_cpd
  INTEGER             :: ncount = 0

  PUBLIC :: les_phy_interface, init_les_phy_interface

  CHARACTER(len=12)  :: str_module = 'les_interface'  ! Output of module for 1 line debug

CONTAINS
  !
  !-----------------------------------------------------------------------
  !
  SUBROUTINE init_les_phy_interface(jg, p_patch, p_int_state, p_metrics)
    INTEGER,                   INTENT(in)     :: jg
    TYPE(t_patch),     TARGET, INTENT(in)     :: p_patch
    TYPE(t_int_state),         INTENT(in)     :: p_int_state
    TYPE(t_nh_metrics),        INTENT(inout)  :: p_metrics

    ! precompute ddz_z metrics
    ! could be done in init_vertical_grid later on
    CALL init_vertical_grid_for_les(jg, p_patch, p_int_state, p_metrics)
  END SUBROUTINE init_les_phy_interface

  !
  !-----------------------------------------------------------------------
  !
  SUBROUTINE les_phy_interface(lcall_phy_jg, linit, lredgrid,      & !input
                            & dt_loc, dt_phy_jg,                   & !input
                            & nstep, mtime_current,                & !input
                            & pt_patch, pt_int_state, p_metrics,   & !input
                            & pt_par_patch,                        & !input
                            & ext_data,                            & !input
                            & pt_prog,                             & !inout
                            & pt_prog_now_rcf, pt_prog_rcf,        & !in/inout
                            & pt_diag ,                            & !inout
                            & prm_diag, prm_nwp_tend, lnd_diag,    &
                            & lnd_prog_now, lnd_prog_new,          & !inout
                            & wtr_prog_now, wtr_prog_new,          & !inout
                            & p_prog_list                          ) !in

    !>
    ! !INPUT PARAMETERS:

    LOGICAL, INTENT(IN)          ::   &             !< physics package time control (switches)
         &                          lcall_phy_jg(:) !< for domain jg
    LOGICAL, INTENT(IN)          :: linit           !< .TRUE. if initialization call (this switch is currently used
                                                    !  to call turbtran in addition to the slow-physics routines
    LOGICAL, INTENT(IN)          :: lredgrid        !< use reduced grid for radiation
    REAL(wp),INTENT(in)          :: dt_loc          !< (advective) time step applicable to local grid level
    REAL(wp),INTENT(in)          :: dt_phy_jg(:)    !< time interval for all physics on jg
    INTEGER, INTENT(in)          :: nstep           !time step counter
    TYPE(datetime), POINTER      :: mtime_current
    TYPE(t_patch),        TARGET,INTENT(in):: pt_patch         !<grid/patch info.
    TYPE(t_patch),        TARGET,INTENT(in):: pt_par_patch     !<grid/patch info (parent grid)
    TYPE(t_int_state),    TARGET,INTENT(in)   :: pt_int_state  !< interpolation state
    TYPE(t_nh_metrics)   ,       INTENT(inout)   :: p_metrics ! SBr,CHa:in->inout
    TYPE(t_external_data),       INTENT(inout):: ext_data
    TYPE(t_nh_diag), TARGET, INTENT(inout)    :: pt_diag       !<the diagnostic variables
    TYPE(t_nh_prog), TARGET, INTENT(inout)    :: pt_prog       !<the prognostic variables
    TYPE(t_nh_prog), TARGET, INTENT(inout)    :: pt_prog_now_rcf !<old state for tke
    TYPE(t_nh_prog), TARGET, INTENT(inout)    :: pt_prog_rcf   !<the RCF prognostic variables
    TYPE(t_nwp_phy_diag),       INTENT(inout) :: prm_diag
    TYPE(t_nwp_phy_tend),TARGET,INTENT(inout) :: prm_nwp_tend
    TYPE(t_lnd_prog),           INTENT(inout) :: lnd_prog_now, lnd_prog_new
    TYPE(t_wtr_prog),           INTENT(inout) :: wtr_prog_now, wtr_prog_new
    TYPE(t_lnd_diag),           INTENT(inout) :: lnd_diag

    TYPE(t_var_list), INTENT(in) :: p_prog_list !current prognostic state list

    ! !OUTPUT PARAMETERS:            !<variables induced by the whole physics
    ! Local array bounds:

    INTEGER :: nlev, nlevp1            !< number of full levels
    INTEGER :: rl_start, rl_end
    INTEGER :: i_startblk, i_endblk    !> blocks
    INTEGER :: i_startidx, i_endidx    !< slices
    INTEGER :: i_nchdom                !< domain index

    ! Local scalars:

    INTEGER :: jc,jk,jb,jce      !loop indices
    INTEGER :: jg                !domain id

    LOGICAL :: ltemp, lpres, ltemp_ifc, l_any_fastphys, l_any_slowphys

    INTEGER,  POINTER ::  iidx(:,:,:), iblk(:,:,:)

    REAL(wp), TARGET :: &                                     !> temporal arrays for
      & z_ddt_temp  (nproma,pt_patch%nlev,pt_patch%nblks_c)   !< Temperature tendency

    REAL(wp) :: z_exner_sv(nproma,pt_patch%nlev,pt_patch%nblks_c)
    REAL(wp) :: z_temp_old(nproma,pt_patch%nlev,pt_patch%nblks_c)

    !< vertical interfaces

    REAL(wp) :: z_airmass (nproma,pt_patch%nlev) !< needed for radheat
    REAL(wp) :: zsct ! solar constant (at time of year)
    REAL(wp) :: zcosmu0 (nproma,pt_patch%nblks_c)

    REAL(wp) :: inv_dt_fastphy

    REAL(wp) :: z_qsum(nproma,pt_patch%nlev)  !< summand of virtual increment
    REAL(wp) :: z_ddt_qsum                    !< summand of tendency of virtual increment

    ! Variables for dpsdt diagnostic
    REAL(wp) :: dps_blk(pt_patch%nblks_c), dpsdt_avg
    INTEGER  :: npoints_blk(pt_patch%nblks_c), npoints

    ! Variables for EDMF DUALM
    REAL(wp) :: qtvar(nproma,pt_patch%nlev)

    ! communication ids, these do not need to be different variables,
    ! since they are not treated individualy
    INTEGER :: ddt_u_tot_comm, ddt_v_tot_comm, tracers_comm, tempv_comm, exner_pr_comm

    CHARACTER(len=max_char_length), PARAMETER :: routine = 'mo_interface_les:les_phy_interface:'

    ! Pointer to IDs of tracers which contain prognostic condensate.
    ! Required for computing the water loading term 
    INTEGER, POINTER :: condensate_list(:)

    TYPE(timeDelta), POINTER            :: time_diff
    REAL(wp)                            :: p_sim_time     !< elapsed simulation time on this grid level

    ! calculate elapsed simulation time in seconds (local time for
    ! this domain!)
    time_diff  => newTimedelta("PT0S")
    time_diff  =  getTimeDeltaFromDateTime(mtime_current, time_config%tc_exp_startdate)
    p_sim_time =  getTotalMillisecondsTimedelta(time_diff, mtime_current)*1.e-3_wp
    CALL deallocateTimedelta(time_diff)

    IF (ltimer) CALL timer_start(timer_physics)

    ! local variables related to the blocking

    i_nchdom  = MAX(1,pt_patch%n_childdom)
    jg        = pt_patch%id

    ! number of vertical levels
    nlev   = pt_patch%nlev
    nlevp1 = pt_patch%nlevp1

    !define pointers
    iidx  => pt_patch%edges%cell_idx
    iblk  => pt_patch%edges%cell_blk

    ! inverse of fast physics time step
    inv_dt_fastphy = 1._wp / dt_phy_jg(itfastphy)

    IF (lcall_phy_jg(itsatad) .OR. lcall_phy_jg(itgscp) .OR. &
        lcall_phy_jg(itturb)  .OR. lcall_phy_jg(itsfc)) THEN
      l_any_fastphys = .TRUE.
    ELSE
      l_any_fastphys = .FALSE.
    ENDIF


    IF (lcall_phy_jg(itrad) .OR. lcall_phy_jg(itccov)) THEN
      l_any_slowphys = .TRUE.
    ELSE
      l_any_slowphys = .FALSE.
    ENDIF

    ! condensate tracer IDs
    condensate_list => advection_config(jg)%ilist_hydroMass

    !Check if time to sample data
    IF(sampl_freq_step > 0)THEN
      IF( .NOT.linit .AND. MOD(nstep,sampl_freq_step)==0 )THEN
         is_sampling_time = .TRUE.
         ncount = ncount + 1
      ELSE
         is_sampling_time = .FALSE.
      END IF
    ELSE
      is_sampling_time = .FALSE.
    END IF

    !Check if time to write data
    IF(avg_interval_step > 0)THEN
      IF( .NOT.linit .AND. MOD(nstep,avg_interval_step)==0 )THEN
         is_writing_time = .TRUE.
      ELSE
         is_writing_time = .FALSE.
      END IF
    ELSE
      is_writing_time = .FALSE.
    END IF

    !-------------------------------------------------------------------------
    !>  Additional syns required for physics before 3D turbulence
    !-------------------------------------------------------------------------
    IF(diffusion_config(jg)%lhdiff_temp .AND. lhdiff_rcf) &
        CALL sync_patch_array_mult(SYNC_C, pt_patch, 2, pt_prog%theta_v, &
                                   pt_prog%exner)

    IF(diffusion_config(jg)%lhdiff_w .AND. lhdiff_rcf) &
       CALL sync_patch_array(SYNC_C, pt_patch, pt_prog%w)

    !add all tracers that are used in satad and turbulence
    IF(ltransport) &
      CALL sync_patch_array_mult(SYNC_C, pt_patch, iqm_max, &
           f4din=pt_prog_rcf%tracer(:,:,:,1:iqm_max)) 

    !-------------------------------------------------------------------------
    !>  Update the tracer for every advective timestep,
    !!  all other updates are done in dynamics
    !-------------------------------------------------------------------------

    IF (.NOT. linit) THEN

      IF (msg_level >= 15) &
           & CALL message(TRIM(routine), 'update_tracers')

      IF (timers_level > 2) CALL timer_start(timer_update_prog_phy)

      rl_start = 1
      rl_end   = min_rlcell

      CALL update_tracers_les ( pt_patch            ,& !in
           &                  dt_phy_jg(itfastphy)  ,& !in
           &                  prm_nwp_tend          ,& !in
           &                  pt_prog_rcf           ,& !inout tracer
           &                  rl_start, rl_end)

      IF (timers_level > 2) CALL timer_stop(timer_update_prog_phy)

    ENDIF

    IF ( lcall_phy_jg(itturb) .OR. linit ) THEN

      !-------------------------------------------------------------------------
      !>
      !!   Interpolation from v_n onto u,v =>  Reconstruct u and v
      !!   This is needed for turbulence, convection and SSO/GWdrag
      !!
      !-------------------------------------------------------------------------

      IF (msg_level >= 15) &
           & CALL message(TRIM(routine), 'reconstruct u/v')

      IF (timers_level > 3) CALL timer_start(timer_phys_u_v)

      CALL rbf_vec_interpol_cell(pt_prog%vn,            & !< normal wind comp.
        &                        pt_patch,              & !< patch
        &                        pt_int_state,          & !< interpolation state
        &                        pt_diag%u, pt_diag%v )   !<  reconstr. u,v wind

      IF (timers_level > 3) CALL timer_stop(timer_phys_u_v)

    ENDIF ! diagnose u/v

    IF (l_any_fastphys .OR. linit) THEN

      ! Diagnose temperature if any of the fast physics schemes is called
      CALL diagnose_pres_temp (p_metrics, pt_prog, pt_prog_rcf,    &
           &                              pt_diag, pt_patch,       &
           &                              opt_calc_temp=.TRUE.,    &
           &                              opt_calc_pres=.FALSE.,   &
           &                              opt_rlend=min_rlcell_int-2 )

      IF (msg_level >= 20) THEN ! Initial diagnostic output
        CALL nwp_diag_output_1(pt_patch, pt_diag, pt_prog_rcf)
      ENDIF

    ENDIF ! fast physics activated

    !!-------------------------------------------------------------------------
    !> Initial saturation adjustment (a second one follows at the end of the microphysics)
    !!-------------------------------------------------------------------------

    IF (lcall_phy_jg(itsatad)) THEN

      IF (msg_level >= 15) CALL message(TRIM(routine), 'satad')
      IF (timers_level > 2) CALL timer_start(timer_satad_v_3D)

      rl_start = 1
      rl_end   = min_rlcell

      i_startblk = pt_patch%cells%start_blk(rl_start,1)
      i_endblk   = pt_patch%cells%end_blk(rl_end,i_nchdom)


!$OMP PARALLEL

        ! Store exner function for sound-wave reduction and open upper boundary condition
        ! this needs to be done for all grid points (including halo points)
      CALL copy(pt_prog%exner(:,:,:), z_exner_sv(:,:,:))
!$OMP BARRIER

!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx,z_qsum) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = i_startblk, i_endblk

        CALL get_indices_c(pt_patch, jb, i_startblk, i_endblk, &
          & i_startidx, i_endidx, rl_start, rl_end)

        !-------------------------------------------------------------------------
        !   call the saturation adjustment
        !-------------------------------------------------------------------------

!#ifdef __BOUNDCHECK
          CALL satad_v_3D( &
               & maxiter  = 10                             ,& !> IN
               & tol      = 1.e-3_wp                       ,& !> IN
               & te       = pt_diag%temp       (:,:,jb)    ,& !> INOUT
               & qve      = pt_prog_rcf%tracer (:,:,jb,iqv),& !> INOUT
               & qce      = pt_prog_rcf%tracer (:,:,jb,iqc),& !> INOUT
               & rhotot   = pt_prog%rho        (:,:,jb)    ,& !> IN
               & idim     = nproma                         ,& !> IN
               & kdim     = nlev                           ,& !> IN
               & ilo      = i_startidx                     ,& !> IN
               & iup      = i_endidx                       ,& !> IN
               & klo      = kstart_moist(jg)               ,& !> IN
               & kup      = nlev                            & !> IN
              !& count, errstat,                              !> OUT
               )

        DO jk = kstart_moist(jg), nlev
          DO jc = i_startidx, i_endidx

            ! calculate virtual temperature from condens' output temperature
            ! taken from SUBROUTINE update_tempv_geopot in hydro_atmos/mo_ha_update_diag.f90
            z_qsum(jc,jk) = SUM(pt_prog_rcf%tracer (jc,jk,jb,condensate_list))

          ENDDO
        ENDDO


        DO jk = kstart_moist(jg), nlev
!DIR$ IVDEP
          DO jc = i_startidx, i_endidx

            pt_diag%tempv(jc,jk,jb) =  pt_diag%temp(jc,jk,jb)                        &
              &                   * ( 1._wp +  vtmpc1                                &
              &                   * pt_prog_rcf%tracer(jc,jk,jb,iqv) - z_qsum(jc,jk) )

            pt_prog%exner(jc,jk,jb) = EXP(rd_o_cpd*LOG(rd_o_p0ref                   &
              &                     * pt_prog%rho(jc,jk,jb)*pt_diag%tempv(jc,jk,jb)))

          ENDDO
        ENDDO
      ENDDO ! nblks

!$OMP END DO NOWAIT
!$OMP END PARALLEL

      IF (timers_level > 2) CALL timer_stop(timer_satad_v_3D)

    ELSE ! satad turned off

!$OMP PARALLEL
      ! Store exner function for sound-wave reduction and open upper boundary condition
      CALL copy(pt_prog%exner(:,:,:), z_exner_sv(:,:,:))
!$OMP END PARALLEL

    ENDIF ! satad

    !!-------------------------------------------------------------------------
    !>  turbulent transfer and diffusion and microphysics
    !!
    !!  Because we consider the following physical processes as fast ones
    !!  we allow here the update of prognostic variables inside the subroutines
    !!  This means that the conversion back to the ICON-prognostic variables
    !!  has to be done afterwards
    !!-------------------------------------------------------------------------

    IF (lcall_phy_jg(itgscp) .OR. lcall_phy_jg(itturb) .OR. lcall_phy_jg(itsfc)) THEN

      IF (msg_level >= 15) &
        & CALL message(TRIM(routine), 'diagnose pressure for fast physics')

      !-------------------------------------------------------------------------
      !> temperature and virtual temperature are already up to date:
      !! thus diagnose only pressure on main and interface levels
      !! =>  opt_calc_pres_nh=.TRUE.
      !-------------------------------------------------------------------------
      CALL diagnose_pres_temp (p_metrics, pt_prog, pt_prog_rcf, &
        & pt_diag, pt_patch,      &
        & opt_calc_temp =.FALSE., &
        & opt_calc_pres =.TRUE.,  &
        & opt_rlend=min_rlcell_int-2)

    ENDIF


    IF ( lcall_phy_jg(itsfc) ) THEN

         !> as pressure is needed only for an approximate adiabatic extrapolation
         !! of the temperature at the lowest model level towards ground level,
         !! a recalculation is not required

         CALL nwp_surface    (  dt_phy_jg(itfastphy),              & !>input
                               & pt_patch,                         & !>input
                               & ext_data,                         & !>input
                               & pt_prog_rcf,                      & !>in/inout rcf=reduced calling freq.
                               & pt_diag ,                         & !>inout
                               & prm_diag,                         & !>inout
                               & lnd_prog_now, lnd_prog_new,       & !>inout
                               & wtr_prog_now, wtr_prog_new,       & !>inout
                               & lnd_diag                          ) !>input

    END IF


    !Call to turbulent parameterization schemes
    IF (  lcall_phy_jg(itturb) ) THEN

      IF (timers_level > 1) CALL timer_start(timer_nwp_turbulence)

      CALL les_turbulence (  dt_phy_jg(itfastphy),              & !>in
                            & pt_patch, p_metrics,              & !>inout, SBr,CHa:in->inout
                            & pt_int_state,                     & !>in
                            & pt_prog,                          & !>in
                            & pt_prog_rcf,                      & !>inout
                            & pt_diag ,                         & !>inout
                            & prm_diag,prm_nwp_tend,            & !>inout
                            & lnd_prog_now,                     & !>in
                            & lnd_prog_new,                     & !>inout ONLY for idealized LES
                            & lnd_diag                          ) !>in

      IF (timers_level > 1) CALL timer_stop(timer_nwp_turbulence)

    END IF


    !-------------------------------------------------------------------------
    !  prognostic microphysic and precipitation scheme
    !-------------------------------------------------------------------------
    IF ( lcall_phy_jg(itgscp)) THEN

      IF (msg_level >= 15) &
        & CALL message(TRIM(routine), 'microphysics')

      !> temperature and tracers have been updated by turbulence;
      !! an update of the pressure field is not needed because pressure
      !! is not needed at high accuracy in the microphysics scheme

      IF (timers_level > 1) CALL timer_start(timer_nwp_microphysics)

      !Copy temp to calculate its tendency next
      CALL copy(pt_diag%temp(:,:,:), z_temp_old(:,:,:)) 

      CALL nwp_microphysics ( dt_phy_jg(itfastphy),             & !>input
                            & pt_patch, p_metrics,              & !>input
                            & pt_prog,                          & !>inout
                            & pt_prog_rcf,                      & !>inout
                            & pt_diag ,                         & !>inout
                            & prm_diag                          ) !>inout

      !Calculate temp tendency due to microphysics in interior points
      rl_start = grf_bdywidth_c+1
      rl_end   = min_rlcell_int

      i_startblk = pt_patch%cells%start_blk(rl_start,1)
      i_endblk   = pt_patch%cells%end_blk(rl_end,i_nchdom)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,i_startidx, i_endidx ) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = i_startblk, i_endblk
        CALL get_indices_c(pt_patch, jb, i_startblk, i_endblk, &
          & i_startidx, i_endidx, rl_start, rl_end )
         DO jk = kstart_moist(jg), nlev
          DO jc =  i_startidx, i_endidx
            prm_nwp_tend%ddt_temp_gscp(jc,jk,jb) =  &
                 ( pt_diag%temp(jc,jk,jb) - z_temp_old(jc,jk,jb) ) * inv_dt_fastphy
          END DO
         END DO
      END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

      IF (timers_level > 1) CALL timer_stop(timer_nwp_microphysics)

    ENDIF


    IF (lart) THEN

      CALL art_washout_interface(pt_prog, pt_diag ,             & !>in
                 &          dt_phy_jg(itfastphy),               & !>in
                 &          pt_patch,                           & !>in
                 &          prm_diag,                           & !>in
                 &          p_metrics,                          & !>in
                 &          pt_prog_rcf%tracer)                   !>inout
    ENDIF !lart



    IF (lcall_phy_jg(itsatad) .OR. lcall_phy_jg(itgscp) .OR. lcall_phy_jg(itturb)) THEN

      IF (timers_level > 1) CALL timer_start(timer_fast_phys)


      ! Remark: in the (unusual) case that satad is used without any other physics,
      ! recalculation of the thermodynamic variables is duplicated here. However,
      ! this is the easiest way to combine minimization of halo communications
      ! with a failsafe flow control

      IF (msg_level >= 15) &
        & CALL message(TRIM(routine), 'recalculate thermodynamic variables')


      ! exclude boundary interpolation zone of nested domains
      rl_start = grf_bdywidth_c+1
      rl_end   = min_rlcell_int

      i_startblk = pt_patch%cells%start_blk(rl_start,1)
      i_endblk   = pt_patch%cells%end_blk(rl_end,i_nchdom)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,i_startidx, i_endidx, z_qsum ) ICON_OMP_DEFAULT_SCHEDULE

      DO jb = i_startblk, i_endblk
        CALL get_indices_c(pt_patch, jb, i_startblk, i_endblk, &
          & i_startidx, i_endidx, rl_start, rl_end )

        !-------------------------------------------------------------------------
        !>
        !! re-calculate saturation adjustment to avoid needle clouds
        !!
        !-------------------------------------------------------------------------

        IF(ANY((/4,5/) == atm_phy_nwp_config(jg)%inwp_gscp))THEN
          CALL satad_v_3D( &
               & maxiter  = 10 ,& !> IN
               & tol      = 1.e-3_wp ,& !> IN
               & te       = pt_diag%temp       (:,:,jb) ,& !> INOUT
               & qve      = pt_prog_rcf%tracer (:,:,jb,iqv),& !> INOUT
               & qce      = pt_prog_rcf%tracer (:,:,jb,iqc),& !> INOUT
               & rhotot   = pt_prog%rho        (:,:,jb) ,& !> IN
               & idim     = nproma ,& !> IN
               & kdim     = nlev ,& !> IN
               & ilo      = i_startidx ,& !> IN
               & iup      = i_endidx ,& !> IN
               & klo      = kstart_moist(jg) ,& !> IN
               & kup      = nlev & !> IN
              !& count, errstat, !> OUT
               )
        ENDIF

        !-------------------------------------------------------------------------
        !>
        !! re-calculate scalar prognostic variables out of physics variables!
        !!
        !-------------------------------------------------------------------------

        IF (kstart_moist(jg) > 1) z_qsum(:,1:kstart_moist(jg)-1) = 0._wp

        DO jk = kstart_moist(jg), nlev
          DO jc = i_startidx, i_endidx

            z_qsum(jc,jk) = SUM(pt_prog_rcf%tracer (jc,jk,jb,condensate_list))

          ENDDO
        ENDDO


        DO jk = 1, nlev
!DIR$ IVDEP
          DO jc =  i_startidx, i_endidx

            pt_diag%tempv(jc,jk,jb) =  pt_diag%temp(jc,jk,jb)          &
&                                  * ( 1._wp +  vtmpc1                 &
&                                  *  pt_prog_rcf%tracer(jc,jk,jb,iqv) &
&                                   - z_qsum(jc,jk) )

            pt_prog%exner(jc,jk,jb) = EXP(rd_o_cpd*LOG(rd_o_p0ref                   &
              &                     * pt_prog%rho(jc,jk,jb)*pt_diag%tempv(jc,jk,jb)))

            pt_diag%exner_pr(jc,jk,jb) = pt_diag%exner_pr(jc,jk,jb) + &
              pt_prog%exner(jc,jk,jb) - z_exner_sv(jc,jk,jb)

            pt_prog%theta_v(jc,jk,  jb) = pt_diag%tempv(jc,jk,jb) &
&                                       / pt_prog%exner(jc,jk,jb)

          ENDDO
        ENDDO


      ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

     IF (timers_level > 1) CALL timer_stop(timer_fast_phys)
    ENDIF ! end of fast physics part

    IF (lcall_phy_jg(itturb) .OR. linit .OR. l_any_slowphys) THEN
      ! re-diagnose pressure
      CALL diagnose_pres_temp (p_metrics, pt_prog, pt_prog_rcf,   &
        &                      pt_diag, pt_patch,                 &
        &                      opt_calc_temp     = .FALSE.,       &
        &                      opt_calc_pres     = .TRUE.,        &
        &                      opt_rlend         = min_rlcell_int )
    ENDIF


    !Calculate turbulent surface exchange coefficient as in mo_nh_interface_nwp
    !but called only if nwp_surface is ON. For LES it uses the method used for GME
    !turbulence i.e NO tiles. For now check if NO TILES approach work. If it does then
    !remove the modifications in mo_surface for tiles approach, and do extra cleaning.
    !If it doesn't then try to make that one work- Status as on 11.09.2013 (AD)

    IF ( (lcall_phy_jg(itturb) .OR. linit) .AND. atm_phy_nwp_config(jg)%inwp_surface>0 ) THEN
!
      IF (timers_level > 1) CALL timer_start(timer_nwp_turbulence)

      ! compute turbulent transfer coefficients (atmosphere-surface interface)
      CALL nwp_turbtrans  ( dt_phy_jg(itfastphy),             & !>in
                          & pt_patch, p_metrics,              & !>in
                          & ext_data,                         & !>in
                          & pt_prog_rcf,                      & !>inout
                          & pt_diag,                          & !>inout
                          & prm_diag,                         & !>inout
                          & wtr_prog_new,                     & !>in
                          & lnd_prog_new,                     & !>inout
                          & lnd_diag                          ) !>inout

      IF (timers_level > 1) CALL timer_stop(timer_nwp_turbulence)
    ENDIF !lcall(itturb)


    !!-------------------------------------------------------------------------
    !!  slow physics part
    !!-------------------------------------------------------------------------

    IF (l_any_slowphys) THEN

      IF (msg_level >= 15) &
         CALL message(TRIM(routine), 'diagnose pres/temp for slow physics')

      ! If slow physics is called without fast physics (which should happen
      ! at the initial time step only), temperature needs to be calculated
      ! Otherwise, temperature is up to date
      IF ( .NOT. (lcall_phy_jg(itgscp) .OR. lcall_phy_jg(itturb))) THEN
        ltemp = .TRUE.
      ELSE
        ltemp = .FALSE.
      ENDIF


      ! Pressure has already been updated at the end of the fast physics part
      lpres = .FALSE.

      ! Temperature at interface levels is needed if irad_aero = 5 or 6
      ! or if Ritter-Geleyn radiation is called
      IF ( lcall_phy_jg(itrad) .AND. ( irad_aero == 5 .OR. irad_aero == 6 &
           .OR. irad_aero == 9 .OR. atm_phy_nwp_config(jg)%inwp_radiation == 2 ) ) THEN
        ltemp_ifc = .TRUE.
      ELSE
        ltemp_ifc = .FALSE.
      ENDIF

      CALL diagnose_pres_temp (p_metrics, pt_prog, pt_prog_rcf,   &
        &                      pt_diag, pt_patch,                 &
        &                      opt_calc_temp     = ltemp,         &
        &                      opt_calc_pres     = lpres,         &
        &                      lnd_prog          = lnd_prog_new,  &
        &                      opt_calc_temp_ifc = ltemp_ifc,     &
        &                      opt_rlend         = min_rlcell_int )

    ENDIF


    !-------------------------------------------------------------------------
    !> Cloud cover
    !-------------------------------------------------------------------------

    IF ( lcall_phy_jg(itccov) ) THEN

      ! When using a reduced grid for radiation, part of the boundary points need
      ! to be included in order to compute spatial gradients for back-interpolation
      IF (lredgrid) THEN
        rl_start = grf_bdywidth_c-1
      ELSE
        rl_start = grf_bdywidth_c+1
      ENDIF
      rl_end   = min_rlcell_int

      i_startblk = pt_patch%cells%start_blk(rl_start,1)
      i_endblk   = pt_patch%cells%end_blk(rl_end,i_nchdom)

      IF (msg_level >= 15) &
        &  CALL message(TRIM(routine), 'cloud cover')

      IF (timers_level > 2) CALL timer_start(timer_cover_koe)

      !-------------------------------------------------------------------------
      !> Cloud water distribution: cloud cover, cloud water, cloud ice
      !  inwp_cldcover =
      !  (0) no clouds
      !  (1) diagnostic cloud cover
      !  (2) prognostic total water variance (not yet started)
      !  (3) clouds as in COSMO
      !  (4) clouds as in turbulence
      !  (5) grid-scale cloud cover [1 or 0]
      !-------------------------------------------------------------------------

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx) ICON_OMP_GUIDED_SCHEDULE
      DO jb = i_startblk, i_endblk

        CALL get_indices_c(pt_patch, jb, i_startblk, i_endblk, &
&                       i_startidx, i_endidx, rl_start, rl_end)

        CALL cover_koe &
&             (kidia  = i_startidx ,   kfdia  = i_endidx  ,       & !! in:  horizonal begin, end indices
&              klon = nproma,  kstart = kstart_moist(jg)  ,       & !! in:  horiz. and vert. vector length
&              klev   = nlev,                                     &
&              icldscheme = atm_phy_nwp_config(jg)%inwp_cldcover, & !! in:  cloud cover option
&              inwp_turb  = atm_phy_nwp_config(jg)%inwp_turb,     & !! in:  turbulence scheme number
&              tt     = pt_diag%temp         (:,:,jb)     ,       & !! in:  temperature at full levels
&              pp     = pt_diag%pres         (:,:,jb)     ,       & !! in:  pressure at full levels
&              ps     = pt_diag%pres_sfc     (:,jb)       ,       & !! in:  surface pressure at full levels
&              t_g    = lnd_prog_new%t_g     (:,jb)       ,       & !! in:  surface temperature
&              pgeo   = p_metrics%geopot_agl (:,:,jb)     ,       & !! in:  geopotential height
&              rho    = pt_prog%rho          (:,:,jb  )   ,       & !! in:  density
&              rcld   = prm_diag%rcld        (:,:,jb)     ,       & !! in:  standard deviation of saturation deficit
&              ldland = ext_data%atm%llsm_atm_c (:,jb)    ,       & !! in:  land/sea mask
&              ldcum  = prm_diag%locum       (:,jb)       ,       & !! in:  convection on/off
&              kcbot  = prm_diag%mbas_con    (:,jb)       ,       & !! in:  convective cloud base
&              kctop  = prm_diag%mtop_con    (:,jb)       ,       & !! in:  convective cloud top
&              pmfude_rate = prm_diag%con_udd(:,:,jb,3)   ,       & !! in:  convective updraft detrainment rate
&              plu         = prm_diag%con_udd(:,:,jb,7)   ,       & !! in:  updraft condensate
&              qv     = pt_prog_rcf%tracer   (:,:,jb,iqv) ,       & !! in:  spec. humidity
&              qc     = pt_prog_rcf%tracer   (:,:,jb,iqc) ,       & !! in:  cloud water
&              qi     = pt_prog_rcf%tracer   (:,:,jb,iqi) ,       & !! in:  cloud ice
&              qs     = pt_prog_rcf%tracer   (:,:,jb,iqs) ,       & !! in:  snow
&              qtvar  = qtvar                             ,       & !! in:  qtvar !ONLY for inwp_turb==iedmf
&              cc_tot = prm_diag%clc         (:,:,jb)     ,       & !! out: cloud cover
&              qv_tot = prm_diag%tot_cld     (:,:,jb,iqv) ,       & !! out: qv       -"-
&              qc_tot = prm_diag%tot_cld     (:,:,jb,iqc) ,       & !! out: clw      -"-
&              qi_tot = prm_diag%tot_cld     (:,:,jb,iqi) )         !! out: ci       -"-

      ENDDO

!$OMP END DO NOWAIT
!$OMP END PARALLEL

      IF (timers_level > 2) CALL timer_stop(timer_cover_koe)

    ENDIF! cloud cover

    !-------------------------------------------------------------------------
    !> Radiation
    !-------------------------------------------------------------------------

    IF ( lcall_phy_jg(itrad) ) THEN

      IF (ltimer) CALL timer_start(timer_nwp_radiation)
      CALL nwp_radiation (lredgrid,              & ! in
           &              p_sim_time,            & ! in
           &              mtime_current,         & ! in
           &              pt_patch,              & ! in
           &              pt_par_patch,          & ! in
           &              ext_data,              & ! in
           &              lnd_diag,              & ! in
           &              pt_prog,               & ! inout
           &              pt_diag,               & ! inout
           &              prm_diag,              & ! inout
           &              lnd_prog_new,          & ! in
           &              wtr_prog_new           ) ! in
      IF (ltimer) CALL timer_stop(timer_nwp_radiation)

    ENDIF


    IF ( lcall_phy_jg(itradheat) ) THEN

      IF (msg_level >= 15) &
        & CALL message(TRIM(routine), 'radiative heating')


      IF (timers_level > 1) CALL timer_start(timer_pre_radiation_nwp)

      CALL pre_radiation_nwp (                      &
        & kbdim      = nproma,                      &
        & p_inc_rad  = dt_phy_jg(itfastphy),        &
        & p_sim_time = p_sim_time,                  &
        & pt_patch   = pt_patch,                    &
        & zsmu0      = zcosmu0,                     &
        & zsct       = zsct )
      IF (timers_level > 1) CALL timer_stop(timer_pre_radiation_nwp)

      ! exclude boundary interpolation zone of nested domains
      rl_start = grf_bdywidth_c+1
      rl_end   = min_rlcell_int

      i_startblk = pt_patch%cells%start_blk(rl_start,1)
      i_endblk   = pt_patch%cells%end_blk(rl_end,i_nchdom)

      IF (timers_level > 2) CALL timer_start(timer_radheat)
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,z_airmass) ICON_OMP_DEFAULT_SCHEDULE
!
      DO jb = i_startblk, i_endblk
        !
        CALL get_indices_c(pt_patch, jb, i_startblk, i_endblk, &
&                       i_startidx, i_endidx, rl_start, rl_end)

        zcosmu0 (i_startidx:i_endidx,jb) = MAX(zcosmu0(i_startidx:i_endidx,jb),0.0_wp)

        !calculate solar incoming flux at TOA
        prm_diag%flxdwswtoa(i_startidx:i_endidx,jb) = zcosmu0(i_startidx:i_endidx,jb) &
          &                                         * zsct                 !zsct by pre_radiation

        !DR computed twice - could be replaced by prep_adv(jg)%rhodz_mc_new
        z_airmass(i_startidx:i_endidx,:) = p_metrics%ddqz_z_full(i_startidx:i_endidx,:,jb) * &
                                           pt_prog%rho(i_startidx:i_endidx,:,jb)

        prm_diag%swflxsfc (:,jb)=0._wp
        prm_diag%lwflxsfc (:,jb)=0._wp
        prm_diag%swflxtoa (:,jb)=0._wp

        IF (atm_phy_nwp_config(jg)%inwp_surface >= 1) THEN

          prm_diag%swflxsfc_t (:,jb,:)=0._wp
          prm_diag%lwflxsfc_t (:,jb,:)=0._wp

          CALL radheat (                   &
          !
          ! input
          ! -----
          !
          & jcs=i_startidx                         ,&! in     start index of inner do loop
          & jce=i_endidx                           ,&! in     end index of inner do loop
          & kbdim=nproma                           ,&! in     loop length and dimension size
          & klev=nlev                              ,&! in     vertical dimension size
          & klevp1=nlevp1                          ,&! in     vertical dimension size
          & ntiles=ntiles_total                    ,&! in     number of tiles of sfc flux fields
          & ntiles_wtr=ntiles_water                ,&! in     number of extra tiles for ocean and lakes
          & pmair=z_airmass                        ,&! in     layer air mass             [kg/m2]
          & pqv=prm_diag%tot_cld(:,:,jb,iqv)       ,&! in     specific moisture           [kg/kg]
          & pcd=cvd                                ,&! in    specific heat of dry air  [J/kg/K]
          & pcv=cvv                                ,&! in    specific heat of vapor    [J/kg/K]
          & pi0=prm_diag%flxdwswtoa(:,jb)          ,&! in     solar incoming flux at TOA  [W/m2]
          & pemiss=ext_data%atm%emis_rad(:,jb)     ,&! in     lw sfc emissivity
          & pqc=prm_diag%tot_cld    (:,:,jb,iqc)   ,&! in     specific cloud water        [kg/kg]
          & pqi=prm_diag%tot_cld    (:,:,jb,iqi)   ,&! in     specific cloud ice          [kg/kg]
          & ppres_ifc=pt_diag%pres_ifc(:,:,jb)     ,&! in     pressure at layer boundaries [Pa]
          & albedo=prm_diag%albdif(:,jb)           ,&! in     grid-box average shortwave albedo
          & albedo_t=prm_diag%albdif_t(:,jb,:)     ,&! in     tile-specific shortwave albedo
          & lp_count=ext_data%atm%lp_count(jb)     ,&! in     number of land points
          & gp_count_t=ext_data%atm%gp_count_t(jb,:),&! in   number of land points per tile
          & spi_count =ext_data%atm%spi_count(jb)  ,&! in     number of seaice points
          & fp_count  =ext_data%atm%fp_count(jb)   ,&! in     number of (f)lake points
          & idx_lst_lp=ext_data%atm%idx_lst_lp(:,jb), &! in   index list of land points
          & idx_lst_t=ext_data%atm%idx_lst_t(:,jb,:), &! in   index list of land points per tile
          & idx_lst_spi=ext_data%atm%idx_lst_spi(:,jb),&! in  index list of seaice points
          & idx_lst_fp=ext_data%atm%idx_lst_fp(:,jb),&! in    index list of (f)lake points
          & cosmu0=zcosmu0(:,jb)                   ,&! in     cosine of solar zenith angle
          & opt_nh_corr=.TRUE.                     ,&! in     switch for NH mode
          & ptsfc=lnd_prog_new%t_g(:,jb)           ,&! in     surface temperature         [K]
          & ptsfc_t=lnd_prog_new%t_g_t(:,jb,:)     ,&! in     tile-specific surface temperature         [K]
          & ptsfctrad=prm_diag%tsfctrad(:,jb)      ,&! in     sfc temp. used for pflxlw   [K]
          & ptrmsw=prm_diag%trsolall (:,:,jb)      ,&! in     shortwave net tranmissivity []
          & pflxlw=prm_diag%lwflxall (:,:,jb)      ,&! in     longwave net flux           [W/m2]
          & lwflx_up_sfc_rs=prm_diag%lwflx_up_sfc_rs(:,jb), &! in longwave upward flux at surface [W/m2]
          & trsol_up_toa=prm_diag%trsol_up_toa(:,jb),   & ! in shortwave upward transm. at the top of the atmosphere
          & trsol_up_sfc=prm_diag%trsol_up_sfc(:,jb),   & ! in shortwave upward transm. at the surface
          & trsol_par_sfc=prm_diag%trsol_par_sfc(:,jb), & ! in photosynthetically active downward transm. at the surface
          & trsol_dn_sfc_diff=prm_diag%trsol_dn_sfc_diff(:,jb),&! in shortwave diffuse downward transm. at the surface
          !
          ! output
          ! ------
          !
          & pdtdtradsw=prm_nwp_tend%ddt_temp_radsw(:,:,jb),&! out    rad. heating by SW        [K/s]
          & pdtdtradlw=prm_nwp_tend%ddt_temp_radlw(:,:,jb),&! out    rad. heating by lw        [K/s]
          & pflxsfcsw =prm_diag%swflxsfc (:,jb)   ,&        ! out shortwave surface net flux [W/m2]
          & pflxsfclw =prm_diag%lwflxsfc (:,jb)   ,&        ! out longwave surface net flux  [W/m2]
          & pflxsfcsw_t=prm_diag%swflxsfc_t (:,jb,:)   ,&   ! out tile-specific shortwave surface net flux [W/m2]
          & pflxsfclw_t=prm_diag%lwflxsfc_t (:,jb,:)   ,&   ! out tile-specific longwave surface net flux  [W/m2]
          & pflxtoasw =prm_diag%swflxtoa (:,jb)        ,&   ! out shortwave toa net flux     [W/m2]
          & lwflx_up_sfc=prm_diag%lwflx_up_sfc(:,jb)   ,&   ! out longwave upward flux at surface [W/m2]
          & swflx_up_toa=prm_diag%swflx_up_toa(:,jb)   ,&   ! out shortwave upward flux at the TOA [W/m2]
          & swflx_up_sfc=prm_diag%swflx_up_sfc(:,jb)   ,&   ! out shortwave upward flux at the surface [W/m2]
          & swflx_par_sfc=prm_diag%swflx_par_sfc(:,jb) ,&   ! out shortwave upward flux at the surface [W/m2]
          & swflx_dn_sfc_diff=prm_diag%swflx_dn_sfc_diff(:,jb) ) ! out shortwave diffuse downward flux at the surface [W/m2]

        ELSE
          CALL radheat (                   &
          !
          ! input
          ! -----
          !
          & jcs=i_startidx                         ,&! in     start index of inner do loop
          & jce=i_endidx                           ,&! in     end index of inner do loop
          & kbdim=nproma                           ,&! in     loop length and dimension size
          & klev=nlev                              ,&! in     vertical dimension size
          & klevp1=nlevp1                          ,&! in     vertical dimension size
          & ntiles=1                               ,&! in     number of tiles of sfc flux fields
          & ntiles_wtr=0                           ,&! in     number of extra tiles for ocean and lakes
          & pmair=z_airmass                        ,&! in     layer air mass             [kg/m2]
          & pqv=prm_diag%tot_cld(:,:,jb,iqv)       ,&! in     specific moisture           [kg/kg]
          & pcd=cvd                                ,&! in     specific heat of dry air  [J/kg/K]
          & pcv=cvv                                ,&! in     specific heat of vapor    [J/kg/K]
          & pi0=prm_diag%flxdwswtoa(:,jb)          ,&! in     solar incoming flux at TOA  [W/m2]
          & pemiss=ext_data%atm%emis_rad(:,jb)     ,&! in     lw sfc emissivity
          & pqc=prm_diag%tot_cld    (:,:,jb,iqc)   ,&! in     specific cloud water        [kg/kg]
          & pqi=prm_diag%tot_cld    (:,:,jb,iqi)   ,&! in     specific cloud ice          [kg/kg]
          & ppres_ifc=pt_diag%pres_ifc(:,:,jb)     ,&! in     pressure at layer boundaries [Pa]
          & cosmu0=zcosmu0(:,jb)                   ,&! in     cosine of solar zenith angle
          & opt_nh_corr=.TRUE.                     ,&! in     switch for NH mode
          & ptsfc=lnd_prog_new%t_g(:,jb)           ,&! in     surface temperature         [K]
          & ptsfctrad=prm_diag%tsfctrad(:,jb)      ,&! in     sfc temp. used for pflxlw   [K]
          & ptrmsw=prm_diag%trsolall (:,:,jb)      ,&! in     shortwave net tranmissivity []
          & pflxlw=prm_diag%lwflxall (:,:,jb)      ,&! in     longwave net flux           [W/m2]
          & lwflx_up_sfc_rs=prm_diag%lwflx_up_sfc_rs(:,jb), &! in longwave upward flux at surface [W/m2]
          & trsol_up_toa=prm_diag%trsol_up_toa(:,jb),   & ! in shortwave upward transm. at the top of the atmosphere
          & trsol_up_sfc=prm_diag%trsol_up_sfc(:,jb),   & ! in shortwave upward transm. at the surface
          & trsol_par_sfc=prm_diag%trsol_par_sfc(:,jb), & ! in photosynthetically active downward transm. at the surface
          & trsol_dn_sfc_diff=prm_diag%trsol_dn_sfc_diff(:,jb),&! in shortwave diffuse downward transm. at the surface
          !
          ! output
          ! ------
          !
          & pdtdtradsw=prm_nwp_tend%ddt_temp_radsw(:,:,jb),&! out    rad. heating by SW        [K/s]
          & pdtdtradlw=prm_nwp_tend%ddt_temp_radlw(:,:,jb),&! out    rad. heating by lw        [K/s]
          & pflxsfcsw =prm_diag%swflxsfc (:,jb)   ,&        ! out shortwave surface net flux [W/m2]
          & pflxsfclw =prm_diag%lwflxsfc (:,jb)   ,&        ! out longwave surface net flux  [W/m2]
          & pflxtoasw =prm_diag%swflxtoa (:,jb)   ,&        ! out shortwave toa net flux     [W/m2]
          & lwflx_up_sfc=prm_diag%lwflx_up_sfc(:,jb)   ,&   ! out longwave upward flux at surface [W/m2]
          & swflx_up_toa=prm_diag%swflx_up_toa(:,jb)   ,&   ! out shortwave upward flux at the TOA [W/m2]
          & swflx_up_sfc=prm_diag%swflx_up_sfc(:,jb)   ,&   ! out shortwave upward flux at the surface [W/m2]
          & swflx_par_sfc=prm_diag%swflx_par_sfc(:,jb) ,&   ! out shortwave upward flux at the surface [W/m2]
          & swflx_dn_sfc_diff=prm_diag%swflx_dn_sfc_diff(:,jb) ) ! out shortwave diffuse downward flux at the surface [W/m2]

        ENDIF

      ENDDO ! blocks

!$OMP END DO NOWAIT
!$OMP END PARALLEL

      IF (timers_level > 2) CALL timer_stop(timer_radheat)

    ENDIF  ! inwp_radiation


    !-------------------------------------------------------------------------
    ! Anurag Dipankar MPIM (2013-May-29)
    ! Large-scale forcing is to be applied at the end of all physics so that
    ! the most updated variable is used. Ideally it should be "next" timestep
    ! variable. Also note that its not actually a part of physics (sub-grid
    ! activity). It is called here to take advantage of u,v.
    !
    ! These LS forcing act as slow process so the tendencies from them are
    ! accumulated with the slow physics tendencies next
    !
    !(2013-25-June) LS forcing is called every physics step
    !-------------------------------------------------------------------------
    IF(is_ls_forcing)THEN

      IF (msg_level >= 15) &
        &  CALL message(TRIM(routine), 'LS forcing')

      IF (timers_level > 3) CALL timer_start(timer_ls_forcing)

      ! exclude boundary interpolation zone of nested domains
      rl_start = grf_bdywidth_c+1
      rl_end   = min_rlcell_int

      CALL apply_ls_forcing ( pt_patch,          &  !>in
        &                     p_metrics,         &  !>in
        &                     pt_prog,           &  !>in
        &                     pt_diag,           &  !>in
        &                     pt_prog_rcf%tracer(:,:,:,iqv),  & !>in
        &                     rl_start,                       & !>in
        &                     rl_end,                         & !>in
        &                     prm_nwp_tend%ddt_u_ls,          & !>out
        &                     prm_nwp_tend%ddt_v_ls,          & !>out
        &                     prm_nwp_tend%ddt_temp_ls,       & !>out
        &                     prm_nwp_tend%ddt_tracer_ls(:,iqv) ) !>out

      IF (timers_level > 3) CALL timer_stop(timer_ls_forcing)

    END IF


    IF (timers_level > 2) CALL timer_start(timer_phys_acc)
    !-------------------------------------------------------------------------
    !>  accumulate tendencies of slow_physics:
    !-------------------------------------------------------------------------
    IF( (l_any_slowphys .OR. lcall_phy_jg(itradheat)) .OR. is_ls_forcing) THEN

      IF (timers_level > 3) CALL timer_start(timer_phys_acc_1)

      ! exclude boundary interpolation zone of nested domains
      rl_start = grf_bdywidth_c+1
      rl_end   = min_rlcell_int

      i_startblk = pt_patch%cells%start_blk(rl_start,1)
      i_endblk   = pt_patch%cells%end_blk(rl_end,i_nchdom)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,i_startidx, i_endidx , z_qsum, z_ddt_qsum &
!$OMP  ) ICON_OMP_DEFAULT_SCHEDULE
!
      DO jb = i_startblk, i_endblk
!
        CALL get_indices_c(pt_patch, jb, i_startblk, i_endblk, &
&                       i_startidx, i_endidx, rl_start, rl_end)

        z_ddt_temp(i_startidx:i_endidx,:,jb) =                                                   &
   &                                       prm_nwp_tend%ddt_temp_radsw(i_startidx:i_endidx,:,jb) &
   &                                    +  prm_nwp_tend%ddt_temp_radlw(i_startidx:i_endidx,:,jb)

        IF (kstart_moist(jg) > 1) z_qsum(:,1:kstart_moist(jg)-1) = 0._wp

        DO jk = kstart_moist(jg), nlev
          DO jc = i_startidx, i_endidx

            z_qsum(jc,jk) = SUM(pt_prog_rcf%tracer (jc,jk,jb,condensate_list))

          ENDDO
        ENDDO



        ! Convert temperature tendency into Exner function tendency
        DO jk = 1, nlev
!DIR$ IVDEP
          DO jc = i_startidx, i_endidx
            pt_diag%ddt_exner_phy(jc,jk,jb) = rd_o_cpd / pt_prog%theta_v(jc,jk,jb)           &
              &                             * z_ddt_temp(jc,jk,jb)                           &
              &                             *(1._wp + vtmpc1*pt_prog_rcf%tracer(jc,jk,jb,iqv)&
              &                             - z_qsum(jc,jk))
          ENDDO
        ENDDO


        !-------------------------------------------------------------------------
        !>  accumulate tendencies of slow_physics when LS forcing is ON
        !-------------------------------------------------------------------------
        IF (is_ls_forcing) THEN

          DO jk = 1, nlev
            DO jc = i_startidx, i_endidx
              z_ddt_temp(jc,jk,jb) = z_ddt_temp(jc,jk,jb)               &
                                +  prm_nwp_tend%ddt_temp_ls(jk)

              ! Convert temperature tendency into Exner function tendency
              z_ddt_qsum =   prm_nwp_tend%ddt_tracer_ls(jk,iqc)          &
                &          + prm_nwp_tend%ddt_tracer_ls(jk,iqi)

              pt_diag%ddt_exner_phy(jc,jk,jb) = rd_o_cpd / pt_prog%theta_v(jc,jk,jb)           &
                &                             * ( z_ddt_temp(jc,jk,jb)                         &
                &                             *(1._wp + vtmpc1*pt_prog_rcf%tracer(jc,jk,jb,iqv)&
                &                             - z_qsum(jc,jk)) + pt_diag%temp(jc,jk,jb)        &
                &           * (vtmpc1 * prm_nwp_tend%ddt_tracer_ls(jk,iqv) - z_ddt_qsum) )

            END DO  ! jc
          END DO  ! jk

        ENDIF ! END of LS forcing tendency accumulation

      ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL
      IF (timers_level > 3) CALL timer_stop(timer_phys_acc_1)

    END IF!END OF slow physics tendency accumulation



    !--------------------------------------------------------
    ! Final section: Synchronization of updated prognostic variables,
    !                interpolation of u/v tendendies to edge points,
    !                and diagnostic computations
    !--------------------------------------------------------

    ! Synchronize tracers if any of the updating (fast-physics) processes was active.
    ! In addition, tempv needs to be synchronized, and in case of lhdiff_rcf, also exner_pr
    IF (l_any_fastphys) THEN

      IF (timers_level > 3) CALL timer_start(timer_phys_sync_tracers)

      IF (use_icon_comm) THEN ! use communication library

        tracers_comm = new_icon_comm_variable(pt_prog_rcf%tracer, pt_patch%sync_cells_not_in_domain,  &
          & status=is_ready, scope=until_sync, name="pt_prog_rcf%tracer")
        tempv_comm = new_icon_comm_variable(pt_diag%tempv, pt_patch%sync_cells_not_in_domain, &
          & status=is_ready, scope=until_sync, name="pt_diag%tempv")

        IF (lhdiff_rcf) THEN
          exner_pr_comm = new_icon_comm_variable(pt_diag%exner_pr, &
            & pt_patch%sync_cells_not_in_domain, &
            & status=is_ready, scope=until_sync, name="pt_diag%exner_pr")
        ENDIF

      ELSE
        IF (lhdiff_rcf) THEN
          CALL sync_patch_array_mult(SYNC_C, pt_patch, ntracer+2, pt_diag%tempv, &
                                     pt_diag%exner_pr, f4din=pt_prog_rcf%tracer)
        ELSE
          CALL sync_patch_array_mult(SYNC_C, pt_patch, ntracer+1, pt_diag%tempv, f4din=pt_prog_rcf%tracer)
        ENDIF

      ENDIF

      IF (timers_level > 3) THEN
        CALL timer_stop(timer_phys_sync_tracers)
      ENDIF
    ENDIF

    !------------------------------------------------------------
    ! sync here the slowphys for aggregation
    !-------------------------------------------------------------------
    IF (use_physics_barrier) THEN
      CALL timer_start(timer_barrier)
      CALL work_mpi_barrier()
      CALL timer_stop(timer_barrier)
    ENDIF
    !-------------------------------------------------------------------
    IF (timers_level > 3) CALL timer_start(timer_phys_sync_ddt_u)
    IF (use_icon_comm) THEN

      IF (lcall_phy_jg(itturb) ) THEN
        ddt_u_tot_comm = new_icon_comm_variable(prm_nwp_tend%ddt_u_turb, &
          & pt_patch%sync_cells_one_edge_in_domain, status=is_ready, scope=until_sync, &
          & name="prm_nwp_tend%ddt_u_turb")
        ddt_v_tot_comm = new_icon_comm_variable(prm_nwp_tend%ddt_v_turb, &
          & pt_patch%sync_cells_one_edge_in_domain, status=is_ready, scope=until_sync, &
          & name="prm_nwp_tend%ddt_v_turb")
      ENDIF

       ! sync everything here
      CALL icon_comm_sync_all()

    ELSE

      IF ( lcall_phy_jg(itturb) ) THEN
        CALL sync_patch_array_mult(SYNC_C1, pt_patch, 2, prm_nwp_tend%ddt_u_turb, &
                                 prm_nwp_tend%ddt_v_turb)
      ENDIF
    ENDIF

    IF (timers_level > 3) CALL timer_stop(timer_phys_sync_ddt_u)
    !------------------------------------------------------------


    !------------------------------------------------------------
    ! compute on the halos
    IF (timers_level > 4) CALL timer_start(timer_phys_acc_par)
    IF (l_any_fastphys) THEN
      IF (my_process_is_mpi_all_parallel() ) THEN

        rl_start = min_rlcell_int-1
        rl_end   = min_rlcell

        i_startblk = pt_patch%cells%start_blk(rl_start,1)
        i_endblk   = pt_patch%cells%end_blk(rl_end,i_nchdom)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,i_startidx, i_endidx) ICON_OMP_DEFAULT_SCHEDULE

        DO jb = i_startblk, i_endblk
          CALL get_indices_c(pt_patch, jb, i_startblk, i_endblk, &
            & i_startidx, i_endidx, rl_start, rl_end )

          IF (lhdiff_rcf) THEN
            DO jk = 1, nlev
              DO jc =  i_startidx, i_endidx

                IF (p_metrics%mask_prog_halo_c(jc,jb)) THEN
                  pt_prog%exner(jc,jk,jb) = EXP(rd_o_cpd*LOG(rd_o_p0ref                   &
                    &                     * pt_prog%rho(jc,jk,jb)*pt_diag%tempv(jc,jk,jb)))

                  pt_prog%theta_v(jc,jk,  jb) = pt_diag%tempv(jc,jk,jb) &
&                                             / pt_prog%exner(jc,jk,jb)

                ENDIF

              ENDDO
            ENDDO
          ELSE
            DO jk = 1, nlev
              DO jc =  i_startidx, i_endidx

                IF (p_metrics%mask_prog_halo_c(jc,jb)) THEN
                  pt_prog%exner(jc,jk,jb) = EXP(rd_o_cpd*LOG(rd_o_p0ref                   &
                    &                     * pt_prog%rho(jc,jk,jb)*pt_diag%tempv(jc,jk,jb)))

                  pt_diag%exner_pr(jc,jk,jb) = pt_diag%exner_pr(jc,jk,jb) + &
                    pt_prog%exner(jc,jk,jb) - z_exner_sv(jc,jk,jb)

                  pt_prog%theta_v(jc,jk,  jb) = pt_diag%tempv(jc,jk,jb) &
&                                             / pt_prog%exner(jc,jk,jb)

                ENDIF

              ENDDO
            ENDDO
          ENDIF
        ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

      ENDIF ! my_process_is_mpi_all_parallel
    ENDIF ! fast-physics synchronization
    IF (timers_level > 4) CALL timer_stop(timer_phys_acc_par)


    ! Initialize fields for runtime diagnostics
    ! In case that average ABS(dpsdt) is diagnosed
    IF (msg_level >= 11) THEN
      dps_blk(:)   = 0._wp
      npoints_blk(:) = 0
    ENDIF

    !-------------------------------------------------------------------------
    !>
    !!    @par Interpolation from  u,v onto v_n
    !!      ddt_vn_phy  =interpol(ddt_u_tot)+interpol(ddt_v_tot)
    !!      Calculate normal velocity at edge midpoints
    !-------------------------------------------------------------------------

    IF (timers_level > 4)  CALL timer_start(timer_phys_acc_2)
!$OMP PARALLEL PRIVATE(rl_start,rl_end,i_startblk,i_endblk)

    ! exclude boundary interpolation zone of nested domains
    rl_start = grf_bdywidth_e+1
    rl_end   = min_rledge_int

    i_startblk = pt_patch%edges%start_blk(rl_start,1)
    i_endblk   = pt_patch%edges%end_blk(rl_end,i_nchdom)


!$OMP DO PRIVATE(jb,jk,jce,i_startidx,i_endidx) ICON_OMP_DEFAULT_SCHEDULE

    DO jb = i_startblk, i_endblk

      CALL get_indices_e(pt_patch, jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, rl_start, rl_end)

      IF ( is_ls_forcing .AND. lcall_phy_jg(itturb) ) THEN

#ifdef __LOOP_EXCHANGE
        DO jce = i_startidx, i_endidx
          DO jk = 1, nlev
#else
!CDIR UNROLL=5
        DO jk = 1, nlev
          DO jce = i_startidx, i_endidx
#endif

            pt_diag%ddt_vn_phy(jce,jk,jb) =   pt_int_state%c_lin_e(jce,1,jb)           &
&                                   * (       prm_nwp_tend%ddt_u_ls(jk)                &
&                                   *  pt_patch%edges%primal_normal_cell(jce,jb,1)%v1  &
&                                   +         prm_nwp_tend%ddt_v_ls(jk)                &
&                                   *  pt_patch%edges%primal_normal_cell(jce,jb,1)%v2 )&
&                                                 + pt_int_state%c_lin_e(jce,2,jb)     &
&                                   * (       prm_nwp_tend%ddt_u_ls(jk)                &
&                                   *  pt_patch%edges%primal_normal_cell(jce,jb,2)%v1  &
&                                   +         prm_nwp_tend%ddt_v_ls(jk)                &
&                                   *  pt_patch%edges%primal_normal_cell(jce,jb,2)%v2 )

            pt_prog%vn(jce,jk,jb) = pt_prog%vn(jce,jk,jb) + dt_loc * (                 &
                                              pt_int_state%c_lin_e(jce,1,jb)           &
&                     * ( prm_nwp_tend%ddt_u_turb(iidx(jce,jb,1),jk,iblk(jce,jb,1))    &
&                                   *  pt_patch%edges%primal_normal_cell(jce,jb,1)%v1  &
&                       + prm_nwp_tend%ddt_v_turb(iidx(jce,jb,1),jk,iblk(jce,jb,1))    &
&                                   *  pt_patch%edges%primal_normal_cell(jce,jb,1)%v2 )&
&                                                 + pt_int_state%c_lin_e(jce,2,jb)     &
&                     * ( prm_nwp_tend%ddt_u_turb(iidx(jce,jb,2),jk,iblk(jce,jb,2))    &
&                                    * pt_patch%edges%primal_normal_cell(jce,jb,2)%v1  &
&                      +  prm_nwp_tend%ddt_v_turb(iidx(jce,jb,2),jk,iblk(jce,jb,2))    &
&                                  *  pt_patch%edges%primal_normal_cell(jce,jb,2)%v2 ) )

          ENDDO
        ENDDO

      ELSE IF (lcall_phy_jg(itturb) ) THEN
#ifdef __LOOP_EXCHANGE
        DO jce = i_startidx, i_endidx
          DO jk = 1, nlev
#else
!CDIR UNROLL=8
        DO jk = 1, nlev
          DO jce = i_startidx, i_endidx
#endif

            pt_prog%vn(jce,jk,jb) = pt_prog%vn(jce,jk,jb) + dt_loc * (                 &
                                              pt_int_state%c_lin_e(jce,1,jb)           &
&                     * ( prm_nwp_tend%ddt_u_turb(iidx(jce,jb,1),jk,iblk(jce,jb,1))    &
&                                   *  pt_patch%edges%primal_normal_cell(jce,jb,1)%v1  &
&                       + prm_nwp_tend%ddt_v_turb(iidx(jce,jb,1),jk,iblk(jce,jb,1))    &
&                                   *  pt_patch%edges%primal_normal_cell(jce,jb,1)%v2 )&
&                                                 + pt_int_state%c_lin_e(jce,2,jb)     &
&                     * ( prm_nwp_tend%ddt_u_turb(iidx(jce,jb,2),jk,iblk(jce,jb,2))    &
&                                    * pt_patch%edges%primal_normal_cell(jce,jb,2)%v1  &
&                      +  prm_nwp_tend%ddt_v_turb(iidx(jce,jb,2),jk,iblk(jce,jb,2))    &
&                                  *  pt_patch%edges%primal_normal_cell(jce,jb,2)%v2 ) )

          ENDDO
        ENDDO

      ENDIF

    ENDDO
!$OMP END DO


    ! Diagnosis of ABS(dpsdt) if msg_level >= 11
    IF (msg_level >= 11) THEN

      rl_start = grf_bdywidth_c+1
      rl_end   = min_rlcell_int

      i_startblk = pt_patch%cells%start_blk(rl_start,1)
      i_endblk   = pt_patch%cells%end_blk(rl_end,i_nchdom)

!$OMP DO PRIVATE(jb,jc,i_startidx,i_endidx) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = i_startblk, i_endblk

        CALL get_indices_c(pt_patch, jb, i_startblk, i_endblk, &
                           i_startidx, i_endidx, rl_start, rl_end)

        DO jc = i_startidx, i_endidx
          ! Note: division by time step follows below
          dps_blk(jb) = dps_blk(jb) + &
            ABS(pt_diag%pres_sfc(jc,jb)-pt_diag%pres_sfc_old(jc,jb))
          npoints_blk(jb) = npoints_blk(jb) + 1
          pt_diag%pres_sfc_old(jc,jb) = pt_diag%pres_sfc(jc,jb)
        ENDDO
      ENDDO
!$OMP END DO NOWAIT
    ENDIF

!$OMP END PARALLEL

    IF (timers_level > 4) CALL timer_stop(timer_phys_acc_2)
    IF (timers_level > 3) CALL timer_start(timer_phys_sync_vn)
    IF (lcall_phy_jg(itturb)) CALL sync_patch_array(SYNC_E, pt_patch, pt_prog%vn)
    IF (timers_level > 3) CALL timer_stop(timer_phys_sync_vn)
    IF (timers_level > 2) CALL timer_stop(timer_phys_acc)


    ! dpsdt diagnostic - omitted in the case of a parallization test (p_test_run) because this
    ! is a purely diagnostic quantitiy, for which it does not make sense to implement an order-invariant
    ! summation
    IF (.NOT. p_test_run .AND. msg_level >= 11) THEN
      dpsdt_avg = SUM(dps_blk)
      npoints   = SUM(npoints_blk)
      dpsdt_avg = global_sum_array(dpsdt_avg)
      npoints   = global_sum_array(npoints)
      dpsdt_avg = dpsdt_avg/(REAL(npoints,wp)*dt_loc)
      ! Exclude initial time step where pres_sfc_old is zero
      IF (dpsdt_avg < 10000._wp/dt_loc) THEN
        WRITE(message_text,'(a,f12.6,a,i3)') 'average |dPS/dt| =',dpsdt_avg,' Pa/s in domain',jg
        CALL message(TRIM(routine), TRIM(message_text))
      ENDIF
    ENDIF

    IF (msg_level >= 20) THEN ! extended diagnostic
      CALL nwp_diag_output_2(pt_patch, pt_prog_rcf, prm_nwp_tend, lcall_phy_jg(itturb))
    ENDIF


    ! time averages, accumulations and vertical integrals
    CALL nwp_statistics(lcall_phy_jg,                    & !in
                        & dt_phy_jg,p_sim_time,          & !in
                        & kstart_moist(jg),              & !in
                        & ih_clch(jg), ih_clcm(jg),      & !in
                        & pt_patch, p_metrics,           & !in
                        & pt_prog, pt_prog_rcf,          & !in
                        & pt_diag,                       & !inout
                        & prm_diag                       ) !inout


    !Special diagnostics for LES runs- 1D, time series
    IF( is_sampling_time )THEN
      CALL calculate_turbulent_diagnostics(                 &
                              & pt_patch,                   & !in
                              & pt_prog,  pt_prog_rcf,      & !in
                              & pt_diag,                    & !in
                              & lnd_prog_new, lnd_diag,     & !in
                              & prm_nwp_tend,               & !in
                              & prm_diag                )     !inout

      !write out time series
      CALL write_time_series(prm_diag%turb_diag_0dvar, mtime_current)
    END IF

    IF( is_writing_time )THEN
      CALL write_vertical_profiles(prm_diag%turb_diag_1dvar, mtime_current, ncount)
      !SBr, CHa
      p_nh_state(1)%metrics%wth3d = p_nh_state(1)%metrics%wth / REAL(ncount,wp)
      p_nh_state(1)%metrics%wth = 0._wp
      p_nh_state(1)%metrics%wqv3d = p_nh_state(1)%metrics%wqv / REAL(ncount,wp)
      p_nh_state(1)%metrics%wqv = 0._wp
      ncount = 0
      prm_diag%turb_diag_1dvar = 0._wp
    END IF


    IF (ltimer) CALL timer_stop(timer_physics)


  END SUBROUTINE les_phy_interface

  !-------------------------------------------------------------------------
  !-------------------------------------------------------------------------


  SUBROUTINE update_tracers_les( pt_patch, pdtime, prm_nwp_tend, &
    &                            pt_prog_rcf, rl_start, rl_end )

    TYPE(t_patch),       INTENT(IN)   :: pt_patch     !!grid/patch info.
    TYPE(t_nwp_phy_tend),TARGET, INTENT(IN):: prm_nwp_tend   !< atm tend vars
    TYPE(t_nh_prog),     INTENT(INOUT):: pt_prog_rcf  !!the tracer field at
                                                      !!reduced calling frequency
    REAL(wp),INTENT(in)            :: pdtime
    INTEGER, INTENT(in)            :: rl_start, rl_end

    ! Local array bounds:
    INTEGER :: i_startblk, i_endblk    !! blocks
    INTEGER :: i_startidx, i_endidx    !! slices
    INTEGER :: i_nchdom                !! domain index

    ! Local scalars:
    INTEGER  :: nlev        !< number of full levels
    INTEGER  :: jb          !block index
    INTEGER  :: jt          !tracers
    INTEGER  :: jk,jc,jg
    REAL(wp) :: zqcn, zqin, zqc, zqi

    jg = pt_patch%id

    ! number of vertical levels
    nlev = pt_patch%nlev

    ! local variables related to the blocking

    i_nchdom  = MAX(1,pt_patch%n_childdom)

    i_startblk = pt_patch%cells%start_blk(rl_start,1)
    i_endblk   = pt_patch%cells%end_blk(rl_end,i_nchdom)


!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,jt,i_startidx,i_endidx,zqc,zqi,zqcn,zqin) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk
        CALL get_indices_c(pt_patch, jb, i_startblk, i_endblk, &
                           i_startidx, i_endidx, rl_start, rl_end)

      IF(is_ls_forcing)THEN
        DO jt=1, nqtendphy  ! qv,qc,qi
          DO jk = kstart_moist(jg), nlev
            DO jc = i_startidx, i_endidx
              pt_prog_rcf%tracer(jc,jk,jb,jt) =MAX(0._wp, pt_prog_rcf%tracer(jc,jk,jb,jt)    &
                &                       + pdtime*prm_nwp_tend%ddt_tracer_ls(jk,jt))
            ENDDO
          ENDDO
        ENDDO
      END IF

!Clipping tracers
      DO jk = kstart_moist(jg), nlev
        DO jc = i_startidx, i_endidx

          zqc = pt_prog_rcf%tracer(jc,jk,jb,iqc)
          zqi = pt_prog_rcf%tracer(jc,jk,jb,iqi)

          zqcn = MIN(0._wp, zqc)
          zqin = MIN(0._wp, zqi)

          pt_prog_rcf%tracer(jc,jk,jb,iqc) = MAX(0._wp, zqc)
          pt_prog_rcf%tracer(jc,jk,jb,iqi) = MAX(0._wp, zqi)

          ! Subtract moisture generated by artificial clipping of QC and QI from QV
          pt_prog_rcf%tracer(jc,jk,jb,iqv) = MAX(0._wp, pt_prog_rcf%tracer(jc,jk,jb,iqv)+zqcn+zqin)

        ENDDO
      ENDDO

      DO jt=iqr, iqm_max  ! qr,qs,etc.
        DO jk = kstart_moist(jg), nlev
          DO jc = i_startidx, i_endidx
            pt_prog_rcf%tracer(jc,jk,jb,jt) =MAX(0._wp, pt_prog_rcf%tracer(jc,jk,jb,jt))
          ENDDO
        ENDDO
      ENDDO

!Clipping for number concentrations
      IF(ANY((/4,5/) == atm_phy_nwp_config(jg)%inwp_gscp))THEN
         DO jt=iqni, iqnc  ! qni,qnr,qns,qng,qnh,qnc
            DO jk = kstart_moist(jg), nlev
               DO jc = i_startidx, i_endidx
                  pt_prog_rcf%tracer(jc,jk,jb,jt) =MAX(0._wp, pt_prog_rcf%tracer(jc,jk,jb,jt))
               ENDDO
            ENDDO
         ENDDO
      END IF

    END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL


  END SUBROUTINE update_tracers_les

END MODULE mo_interface_les

