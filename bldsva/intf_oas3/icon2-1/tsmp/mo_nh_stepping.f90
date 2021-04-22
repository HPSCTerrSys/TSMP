!>
!! Initializes and controls the time stepping in the nonhydrostatic model.
!!
!!
!! @par Revision History
!! Initial release by Almut Gassmann, MPI-M (2009-02-06)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
!! The time stepping does eventually perform an IAU with the follwing
!! characteristics:
!!
!! IAU iteration
!!
!!                     input
!!                       /
!!                      /
!!                     / 
!!          ........../
!!         /         
!!        /
!!       /
!!      /
!!     /
!!  -90min               0min              90min         
!! ---|------------------|------------------|------------->
!!    |//////////////////| - - - - - - - - - - - - - - - ->                  
!!                       /       free forecast (iteration = false)             
!!                      /
!!                     /
!!          ........../
!!         /   reset           
!!        /   
!!       /
!!      /
!!     /
!!  -90min               0min              90min         
!! ---|------------------|------------------|------------->
!!    |//////////////////|//////////////////| free forecast                
!!
!!    \_______IAU________/  
!!
!----------------------------
#include "omp_definitions.inc"
!----------------------------

MODULE mo_nh_stepping
!-------------------------------------------------------------------------
!
!    ProTeX FORTRAN source: Style 2
!    modified for ICON project, DWD/MPI-M 2006
!
!-------------------------------------------------------------------------
!
!

  USE mo_kind,                     ONLY: wp, vp
  USE mo_io_units
  USE mo_nonhydro_state,           ONLY: p_nh_state, p_nh_state_lists
  USE mo_nonhydrostatic_config,    ONLY: lhdiff_rcf, itime_scheme, divdamp_order,                     &
    &                                    divdamp_fac, divdamp_fac_o2, ih_clch, ih_clcm, kstart_moist, &
    &                                    ndyn_substeps, ndyn_substeps_var, ndyn_substeps_max
  USE mo_diffusion_config,         ONLY: diffusion_config
  USE mo_dynamics_config,          ONLY: nnow,nnew, nnow_rcf, nnew_rcf, nsav1, nsav2, idiv_method
  USE mo_io_config,                ONLY: is_totint_time, n_diag
  USE mo_parallel_config,          ONLY: nproma, itype_comm, iorder_sendrecv, use_async_restart_output, &
                                         num_prefetch_proc
  USE mo_run_config,               ONLY: ltestcase, dtime, nsteps, ldynamics, ltransport,   &
    &                                    ntracer, iforcing, msg_level, test_mode,           &
    &                                    output_mode, lart
  USE mo_echam_phy_config,         ONLY: echam_phy_config
  USE mo_advection_config,         ONLY: advection_config
  USE mo_radiation_config,         ONLY: albedo_type
  USE mo_timer,                    ONLY: ltimer, timers_level, timer_start, timer_stop,   &
    &                                    timer_total, timer_model_init, timer_nudging,    &
    &                                    timer_bdy_interp, timer_feedback, timer_nesting, &
    &                                    timer_integrate_nh, timer_nh_diagnostics,        &
    &                                    timer_iconam_echam
  USE mo_atm_phy_nwp_config,       ONLY: dt_phy, atm_phy_nwp_config, iprog_aero
  USE mo_nwp_phy_init,             ONLY: init_nwp_phy, init_cloud_aero_cpl
  USE mo_nwp_phy_state,            ONLY: prm_diag, prm_nwp_tend, phy_params
  USE mo_lnd_nwp_config,           ONLY: nlev_soil, nlev_snow, sstice_mode
  USE mo_nwp_lnd_state,            ONLY: p_lnd_state
  USE mo_ext_data_state,           ONLY: ext_data
  USE mo_ext_data_init,            ONLY: interpol_monthly_mean
  USE mo_extpar_config,            ONLY: itopo
  USE mo_limarea_config,           ONLY: latbc_config
  USE mo_model_domain,             ONLY: p_patch, t_patch
  USE mo_time_config,              ONLY: time_config
  USE mo_grid_config,              ONLY: n_dom, lfeedback, ifeedback_type, l_limited_area, &
    &                                    n_dom_start, lredgrid_phys, start_time, end_time, patch_weight
  USE mo_nh_testcases,             ONLY: init_nh_testcase
  USE mo_nh_testcases_nml,         ONLY: nh_test_name, rotate_axis_deg, lcoupled_rho, is_toy_chem
  USE mo_nh_pa_test,               ONLY: set_nh_w_rho
  USE mo_nh_df_test,               ONLY: get_nh_df_velocity
  USE mo_nh_dcmip_hadley,          ONLY: set_nh_velocity_hadley
  USE mo_nh_dcmip_terminator,      ONLY: dcmip_terminator_interface
  USE mo_nh_supervise,             ONLY: supervise_total_integrals_nh, print_maxwinds,  &
    &                                    init_supervise_nh, finalize_supervise_nh
  USE mo_intp_data_strc,           ONLY: p_int_state, t_int_state
  USE mo_intp_rbf,                 ONLY: rbf_vec_interpol_cell
  USE mo_intp,                     ONLY: verts2cells_scalar
  USE mo_grf_intp_data_strc,       ONLY: p_grf_state
  USE mo_gridref_config,           ONLY: l_density_nudging, grf_intmethod_e
  USE mo_grf_bdyintp,              ONLY: interpol_scal_grf
  USE mo_nh_nest_utilities,        ONLY: compute_tendencies, boundary_interpolation,    &
                                         complete_nesting_setup, prep_bdy_nudging,      &
                                         outer_boundary_nudging, nest_boundary_nudging, &
                                         prep_rho_bdy_nudging, density_boundary_nudging,&
                                         prep_outer_bdy_nudging, save_progvars
  USE mo_nh_feedback,              ONLY: feedback, relax_feedback
  USE mo_exception,                ONLY: message, message_text, finish
  USE mo_impl_constants,           ONLY: SUCCESS, MAX_CHAR_LENGTH, iphysproc, iphysproc_short,     &
    &                                    itconv, itccov, itrad, itradheat, itsso, itsatad, itgwd,  &
    &                                    inoforcing, iheldsuarez, inwp, iecham,                    &
    &                                    itturb, itgscp, itsfc,                                    &
    &                                    MODE_IAU, MODE_IAU_OLD, MODIS, SSTICE_ANA_CLINC,          &
    &                                    SSTICE_CLIM, SSTICE_AVG_MONTHLY, SSTICE_AVG_DAILY, max_dom
  USE mo_math_divrot,              ONLY: rot_vertex, div_avg !, div
  USE mo_solve_nonhydro,           ONLY: solve_nh
  USE mo_update_dyn,               ONLY: add_slowphys
  USE mo_advection_stepping,       ONLY: step_advection
  USE mo_advection_aerosols,       ONLY: aerosol_2D_advection, setup_aerosol_advection
  USE mo_integrate_density_pa,     ONLY: integrate_density_pa
  USE mo_nh_dtp_interface,         ONLY: prepare_tracer, compute_airmass
  USE mo_nh_diffusion,             ONLY: diffusion
  USE mo_mpi,                      ONLY: proc_split, push_glob_comm, pop_glob_comm, p_bcast, p_comm_work

#ifdef NOMPI
  USE mo_mpi,                      ONLY: my_process_is_mpi_all_seq
#endif

  USE mo_sync,                     ONLY: sync_patch_array_mult, sync_patch_array, SYNC_C, global_max
  USE mo_nh_interface_nwp,         ONLY: nwp_nh_interface
  USE mo_interface_iconam_echam,   ONLY: interface_iconam_echam
  USE mo_echam_phy_memory,         ONLY: prm_tend
  USE mo_phys_nest_utilities,      ONLY: interpol_phys_grf, feedback_phys_diag, interpol_rrg_grf, copy_rrg_ubc
  USE mo_vertical_grid,            ONLY: set_nh_metrics
  USE mo_nh_diagnose_pres_temp,    ONLY: diagnose_pres_temp
  USE mo_nh_held_suarez_interface, ONLY: held_suarez_nh_interface
  USE mo_master_config,            ONLY: isRestart
  USE mo_restart_attributes,       ONLY: t_RestartAttributeList, getAttributesForRestarting
  USE mo_meteogram_config,         ONLY: meteogram_output_config
  USE mo_meteogram_output,         ONLY: meteogram_sample_vars, meteogram_is_sample_step
  USE mo_name_list_output,         ONLY: write_name_list_output, istime4name_list_output
  USE mo_name_list_output_init,    ONLY: output_file
  USE mo_pp_scheduler,             ONLY: new_simulation_status, pp_scheduler_process
  USE mo_pp_tasks,                 ONLY: t_simulation_status
  USE mo_art_diagnostics_interface,ONLY: art_diagnostics_interface
  USE mo_art_emission_interface,   ONLY: art_emission_interface
  USE mo_art_sedi_interface,       ONLY: art_sedi_interface
  USE mo_art_tools_interface,      ONLY: art_tools_interface

  USE mo_nwp_sfc_utils,            ONLY: aggregate_landvars, update_sstice, update_ndvi
  USE mo_nh_init_nest_utils,       ONLY: initialize_nest
  USE mo_nh_init_utils,            ONLY: hydro_adjust_downward, compute_iau_wgt, save_initial_state, &
                                         restore_initial_state
  USE mo_td_ext_data,              ONLY: set_actual_td_ext_data
  USE mo_initicon_config,          ONLY: init_mode, timeshift, init_mode_soil, is_avgFG_time, &
                                         iterate_iau, dt_iau
  USE mo_initicon_utils,           ONLY: average_first_guess, reinit_average_first_guess
  USE mo_synsat_config,            ONLY: lsynsat
  USE mo_rttov_interface,          ONLY: rttov_driver, copy_rttov_ubc
  USE mo_ls_forcing_nml,           ONLY: is_ls_forcing
  USE mo_ls_forcing,               ONLY: init_ls_forcing
  USE mo_sync_latbc,               ONLY: prepare_latbc_data,                    &
    &                                    read_latbc_data_sync=>read_latbc_data, &
    &                                    deallocate_latbc_data, p_latbc_data,   &
    &                                    read_latbc_tlev, last_latbc_tlev,      &
    &                                    update_lin_interc
  USE mo_interface_les,            ONLY: les_phy_interface
  USE mo_restart,                  ONLY: t_RestartDescriptor, createRestartDescriptor, deleteRestartDescriptor
  USE mo_nh_prepadv_types,         ONLY: prep_adv, t_prepare_adv, jstep_adv
  USE mo_action,                   ONLY: reset_act
  USE mo_output_event_handler,     ONLY: get_current_jfile
  USE mo_nwp_diagnosis,            ONLY: nwp_diag_for_output
  USE mo_turbulent_diagnostic,     ONLY: calculate_turbulent_diagnostics, &
                                         write_vertical_profiles, write_time_series, &
                                         sampl_freq_step, les_cloud_diag
  USE mo_opt_diagnostics,          ONLY: update_opt_acc, reset_opt_acc, &
    &                                    calc_mean_opt_acc, p_nh_opt_diag
  USE mo_var_list,                 ONLY: nvar_lists, var_lists, print_var_list
  USE mo_async_latbc_utils,        ONLY: recv_latbc_data, update_lin_interpolation
  USE mo_async_latbc_types,        ONLY: t_latbc_data
  USE mo_nonhydro_types,           ONLY: t_nh_state
  USE mo_interface_les,            ONLY: init_les_phy_interface
  USE mo_fortran_tools,            ONLY: swap, copy, init
  USE mtime,                       ONLY: datetime, newDatetime, deallocateDatetime, datetimeToString,     &
       &                                 timedelta, newTimedelta, deallocateTimedelta, timedeltaToString, &
       &                                 MAX_DATETIME_STR_LEN, MAX_TIMEDELTA_STR_LEN, newDatetime,        &
       &                                 MAX_MTIME_ERROR_STR_LEN, no_error, mtime_strerror,               &
       &                                 OPERATOR(-), OPERATOR(+), OPERATOR(>), OPERATOR(*),              &
       &                                 ASSIGNMENT(=), OPERATOR(==), OPERATOR(>=), OPERATOR(/=),         &
       &                                 event, eventGroup, newEvent,                                     &
       &                                 addEventToEventGroup, isCurrentEventActive,                      &
       &                                 getTotalMillisecondsTimedelta, getTotalSecondsTimedelta,     &
       &                                 getTimedeltaFromDatetime
  USE mo_event_manager,            ONLY: initEventManager, addEventGroup, getEventGroup, printEventGroup
  USE mo_derived_variable_handling, ONLY: perform_accumulation, reset_accumulation

#if defined( _OPENACC )
  USE mo_nonhydro_gpu_types,       ONLY: save_convenience_pointers, refresh_convenience_pointers
  USE mo_mpi,                      ONLY: i_am_accel_node, my_process_is_work
#endif

#if defined( MESSY )
  USE mo_timer,                    ONLY: timer_extra1, timer_extra2, timer_extra3, timer_extra4    &
       &                               , timer_extra5, timer_extra6, timer_extra7
#endif

#ifdef COUP_OAS_ICON
  USE oas_icon_define
  USE mo_parallel_config,          ONLY: idx_1d
  USE mo_nonhydro_types,           ONLY: t_nh_metrics, t_nh_prog, t_nh_diag
  USE mo_run_config,               ONLY: iqv
  USE mo_impl_constants,           ONLY: min_rlcell_int, grf_bdywidth_c
  USE mo_loopindices,              ONLY: get_indices_c
#endif

  IMPLICIT NONE

  PRIVATE

  !> module name string
  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_nh_stepping'


  ! additional flow control variables that need to be dimensioned with the
  ! number of model domains
  LOGICAL, ALLOCATABLE :: lcall_phy(:,:) ! contains information which physics package
                                         ! must be called at the current timestep
                                         ! and on the current domain.

  REAL(wp), ALLOCATABLE :: t_elapsed_phy(:,:)  ! time (in s) since the last call of
                                               ! the corresponding physics package
                                               ! (fast physics packages are treated as one)

  LOGICAL, ALLOCATABLE :: linit_dyn(:)  ! determines whether dynamics must be initialized
                                        ! on given patch


  ! additional time control variables which are not dimensioned with the number
  ! of model domains
  LOGICAL :: map_phyproc(iphysproc,iphysproc_short) !< mapping matrix
  INTEGER :: iproclist(iphysproc)  !< x-axis of mapping matrix

  ! event handling manager, wrong place, have to move later

  TYPE(eventGroup), POINTER :: checkpointEventGroup => NULL()

  PUBLIC :: prepare_nh_integration
  PUBLIC :: perform_nh_stepping

  TYPE t_datetime_ptr
    TYPE(datetime), POINTER :: ptr => NULL()
  END TYPE t_datetime_ptr

  CONTAINS

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!
!
  !>
  !! Initialisation of the nonhydrostatic state and initial conditions.
  !!
  !!
  !! @par Revision History
  !! Initial release by Almut Gassmann, (2009-03-06)
  !!
  SUBROUTINE prepare_nh_integration
!
  INTEGER :: ntl, jg

!-----------------------------------------------------------------------

  ! for the split explict scheme, ntl is always 2
  ntl = 2


  CALL set_nh_metrics(p_patch(1:), p_nh_state, p_int_state(1:), ext_data)

  IF (n_dom > 1) THEN
    CALL complete_nesting_setup()
  ENDIF

  IF (ltestcase) THEN

    CALL init_nh_testcase(p_patch(1:), p_nh_state, p_int_state(1:), p_lnd_state(1:), &
                        & ext_data, ntl)

    IF(is_ls_forcing) &
       CALL init_ls_forcing(p_nh_state(1)%metrics)
  ENDIF

  IF (iforcing == inwp) THEN
    CALL setup_time_ctrl_physics( )
  END IF

  ! init LES
  DO jg = 1 , n_dom
   IF(atm_phy_nwp_config(jg)%is_les_phy) THEN
     CALL init_les_phy_interface(jg, p_patch(jg), p_int_state(jg), &
       p_nh_state(jg)%metrics)
   END IF
  END DO

  END SUBROUTINE prepare_nh_integration
  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  !>
  !! Organizes nonhydrostatic time stepping
  !! Currently we assume to have only one grid level.
  !!
  !! @par Revision History
  !! Initial release by Almut Gassmann, (2009-04-15)
  !!
  SUBROUTINE perform_nh_stepping (mtime_current, latbc)
    !
    TYPE(datetime),     POINTER       :: mtime_current     !< current datetime (mtime)
    TYPE(t_latbc_data), INTENT(INOUT) :: latbc             !< data structure for async latbc prefetching

  TYPE(t_simulation_status)            :: simulation_status

  CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
    &  routine = modname//':perform_nh_stepping'

  INTEGER                              :: jg, jgc, jn
  INTEGER                              :: ierr
#ifdef COUP_OAS_ICON
  INTEGER :: rl_start, rl_end, i_startblk, i_endblk, jb, jc, &
    i_startidx, i_endidx, ii, oas_i
  CHARACTER(len=1024) :: oas_message
  INTEGER :: oas_prt = 6
#endif

!!$  INTEGER omp_get_num_threads
!!$  INTEGER omp_get_max_threads
!!$  INTEGER omp_get_max_active_levels
!-----------------------------------------------------------------------

  IF (timers_level > 3) CALL timer_start(timer_model_init)

  CALL allocate_nh_stepping ()

  ! Compute diagnostic dynamics fields for initial output and physics initialization
  CALL diag_for_output_dyn ()

  IF (.NOT. isRestart()) THEN
    IF (timeshift%dt_shift < 0._wp) THEN
      ! add IAU time shift interval to current date (mtime):
      mtime_current = mtime_current + timeshift%mtime_shift
    ENDIF
  ENDIF

  ! diagnose airmass from \rho(now) for both restart and non-restart runs
  ! airmass_new required by initial physics call (init_slowphysics)
  ! airmass_now not needed, since ddt_temp_dyn is not computed during the
  ! initial slow physics call.
  DO jg=1, n_dom
    CALL compute_airmass(p_patch(jg),                  &
      &                  p_nh_state(jg)%metrics,       &
      &                  p_nh_state(jg)%prog(nnow(jg)),&
      &                  p_nh_state(jg)%diag, itlev = 2)

    ! initialize exner_pr if the model domain is active
    IF (p_patch(jg)%ldom_active .AND. .NOT. isRestart()) CALL init_exner_pr(jg, nnow(jg))
  ENDDO


  IF (ANY((/SSTICE_CLIM,SSTICE_AVG_MONTHLY,SSTICE_AVG_DAILY/) == sstice_mode) .AND. iforcing == inwp) THEN
    ! t_seasfc and fr_seaice have to be set again from the ext_td_data files;
    ! the values from the analysis have to be overwritten.
    ! In the case of a restart, the call is required to open the file and read the data
    CALL set_actual_td_ext_data (.TRUE., mtime_current, mtime_current, sstice_mode,  &
                                &  p_patch(1:), ext_data, p_lnd_state)
  END IF

  ! Save initial state if IAU iteration mode is chosen
  IF (iterate_iau .AND. .NOT. isRestart()) THEN
    CALL save_initial_state(p_patch(1:), p_nh_state, prm_diag, p_lnd_state, ext_data)
    WRITE(message_text,'(a)') 'IAU iteration is activated: Start of first cycle with halved IAU window'
    CALL message('',message_text)
  ENDIF

  SELECT CASE (iforcing)
  CASE (inwp)
    DO jg=1, n_dom
      IF (.NOT. p_patch(jg)%ldom_active) CYCLE
      CALL init_nwp_phy(                            &
           & p_patch(jg)                           ,&
           & p_nh_state(jg)%metrics                ,&
           & p_nh_state(jg)%prog(nnow(jg))         ,&
           & p_nh_state(jg)%diag                   ,&
           & prm_diag(jg)                          ,&
           & prm_nwp_tend(jg)                      ,&
           & p_lnd_state(jg)%prog_lnd(nnow_rcf(jg)),&
           & p_lnd_state(jg)%prog_lnd(nnew_rcf(jg)),&
           & p_lnd_state(jg)%prog_wtr(nnow_rcf(jg)),&
           & p_lnd_state(jg)%prog_wtr(nnew_rcf(jg)),&
           & p_lnd_state(jg)%diag_lnd              ,&
           & ext_data(jg)                          ,&
           & phy_params(jg), mtime_current )

      IF (.NOT.isRestart()) THEN
        CALL init_cloud_aero_cpl (mtime_current, p_patch(jg), p_nh_state(jg)%metrics, ext_data(jg), prm_diag(jg))
      ENDIF

      IF (iprog_aero == 1) CALL setup_aerosol_advection(p_patch(jg))

    ENDDO
    IF (.NOT.isRestart()) THEN

#ifdef COUP_OAS_ICON
    !CALL diagnose_pres_temp (p_metrics, pt_prog, pt_prog_rcf,   &
    !     &                      pt_diag, pt_patch,                 &
    !     &                      lnd_prog          = lnd_prog_new)
    
    CALL diag_for_output_dyn()

    IF(.NOT.atm_phy_nwp_config(1)%is_les_phy) THEN
      ! diagnostics which are only required for output
      CALL nwp_diag_for_output(mtime_current, kstart_moist(1), & !in
        &                      ih_clch(1), ih_clcm(1), & !in
        &                      phy_params(1), & !in
        &                      p_patch(1), & !in
        &                      p_nh_state(1)%metrics, & !in
        &                      p_nh_state(1)%prog(nnow(1)), & !in  !nnow or nnew?
        &                      p_nh_state(1)%prog(nnow_rcf(1)), & !in  !nnow or nnew?
        &                      p_nh_state(1)%diag, & !in
        &                      p_lnd_state(1)%diag_lnd, & !in
        & p_lnd_state(1)%prog_lnd(nnow_rcf(1)), & !in
        & p_lnd_state(1)%prog_wtr(nnow_rcf(1)), & !inout
        &                      ext_data(1), & !in
        &                      prm_diag(1) ) !inout
    END IF

    IF(msg_level >= 30 ) then !SPo
      WRITE(oas_prt,*) 'iconoas: ',TRIM(routine),' oasis send fields'
      FLUSH(oas_prt)
    END IF

    oas_snd_field(:,7) = 0._wp
    oas_snd_field(:,8) = 0._wp
    
    rl_start = grf_bdywidth_c+1
    rl_end   = min_rlcell_int
    i_startblk = p_patch(1)%cells%start_blk(rl_start, 1)
    i_endblk   = p_patch(1)%cells%end_blk(rl_end, MAX(1,p_patch(1)%n_childdom))
    DO jb = i_startblk, i_endblk
      CALL get_indices_c(p_patch(1), jb, i_startblk, i_endblk, i_startidx, &
        i_endidx, rl_start, rl_end)
      DO jc = i_startidx, i_endidx
        ii = idx_1d(jc,jb)
        oas_snd_field(ii,1)  = p_nh_state(1)%diag%temp(jc,p_patch(1)%nlev,jb)
        oas_snd_field(ii,2)  = p_nh_state(1)%diag%u(jc,p_patch(1)%nlev,jb)  ! Slavko: check on which level !!
        oas_snd_field(ii,3)  = p_nh_state(1)%diag%v(jc,p_patch(1)%nlev,jb)
        oas_snd_field(ii,4)  = p_nh_state(1)%prog(nnow(1))%tracer(jc,p_patch(1)%nlev,jb,iqv)
        oas_snd_field(ii,5)  = &
          p_nh_state(1)%metrics%z_mc(jc,p_patch(1)%nlev,jb) - p_nh_state(1)%metrics%z_ifc(jc,p_patch(1)%nlev+1,jb)
        oas_snd_field(ii,6)  = p_nh_state(1)%diag%pres_sfc(jc,jb)
        oas_snd_field(ii,7)  = MAX(0._wp, prm_diag(1)%swflxsfc(jc,jb) - prm_diag(1)%swflx_dn_sfc_diff(jc,jb))
        oas_snd_field(ii,8)  = prm_diag(1)%swflx_dn_sfc_diff(jc,jb)
        oas_snd_field(ii,9)  = -prm_diag(1)%lwflxsfc(jc,jb) - prm_diag(1)%lwflx_up_sfc(jc,jb)
        oas_snd_field(ii,10) = prm_diag(1)%rain_con_rate(jc,jb) + prm_diag(1)%snow_con_rate(jc,jb)
        oas_snd_field(ii,11) = prm_diag(1)%rain_gsp_rate(jc,jb) + prm_diag(1)%snow_gsp_rate(jc,jb)! + &
        !  prm_diag(1)%ice_gsp_rate(jc,jb) + prm_diag(1)%graupel_gsp_rate(jc,jb) + &
        ! prm_diag(1)%hail_gsp_rate(jc,jb)
      END DO
    END DO

    DO ii = 1, SIZE(oas_snd_meta)
      CALL oasis_put(oas_snd_meta(ii)%vid, 0, oas_snd_field(:,ii), oas_error)
      WRITE(oas_message,*) 'Init sending  ', oas_snd_meta(ii)%clpname, 'at 0: ', &
        MINVAL(oas_snd_field(:,ii)), MAXVAL(oas_snd_field(:,ii))
      CALL message(routine, oas_message)
      IF (oas_error .NE. OASIS_Ok .AND. oas_error .LT. OASIS_Sent) THEN
        WRITE(oas_message,*) 'Init failure in oasis_put of ', oas_snd_meta(ii)%clpname
        CALL oasis_abort(oas_comp_id, oas_comp_name, oas_message)
      END IF
    END DO

    IF(msg_level >= 30 ) then !SPo
      WRITE(oas_prt,*) 'iconoas: ',TRIM(routine),' oasis receive fields'
      FLUSH(oas_prt)
    END IF

    DO oas_i = 1, SIZE(oas_rcv_meta)
      WRITE(oas_message,*) 'Init receiving  ', oas_rcv_meta(oas_i)%clpname,&
        " at t=1"
      CALL message(routine, oas_message)
      CALL oasis_get(oas_rcv_meta(oas_i)%vid, 0, oas_rcv_field(:,oas_i), oas_error)
      WRITE(oas_message,*) 'Init received  ', oas_rcv_meta(oas_i)%clpname, 'at t=0'
      CALL message(routine, oas_message)
      IF (oas_error .NE. OASIS_Ok .AND. oas_error .LT. OASIS_Recvd) THEN
        WRITE(oas_message,*) 'Init failure in oasis_get of ', oas_rcv_meta(oas_i)%clpname
        CALL oasis_abort(oas_comp_id, oas_comp_name, oas_message)
      END IF

      ! transform thru oasis recieved variable into something icon can use
      !
      !DO jb = 1, p_patch(jg)%nblks_c
      !  DO jc = 
      !    glb_idx_1d = p_patch(jg)%cells%decomp_info%glb_index(idx_1d(jc,jb))
      !    oas_rcv_field(jc,jb) = oas_rcv_field(oas_i,glb_idx_1d)
      !  END DO
      !END DO
      
      IF(msg_level > 30 ) THEN
        WRITE(oas_message,*) "ICONOAS: Init icon transforming received vars"
        CALL message(routine, oas_message)
      END IF

      rl_start = grf_bdywidth_c+1
      rl_end   = min_rlcell_int
      i_startblk = p_patch(1)%cells%start_blk(rl_start, 1)
      i_endblk   = p_patch(1)%cells%end_blk(rl_end, MAX(1,p_patch(1)%n_childdom))
      DO jb = i_startblk, i_endblk
        CALL get_indices_c(p_patch(1), jb, i_startblk, i_endblk, i_startidx,&
          i_endidx, rl_start, rl_end)
        DO jc = i_startidx, i_endidx
          ii = idx_1d(jc,jb)
          oas_rcv_field_icon(jc,jb,oas_i) = oas_rcv_field(ii,oas_i)
        END DO
      END DO

      IF(msg_level > 30 ) THEN
        WRITE(oas_message,*) "ICONOAS: Init icon transformed received vars"
        CALL message(routine, oas_message)

      ! check:
        WRITE(oas_message,*) "ICONOAS: Init for ", oas_i, " got rank0-local min, max=", &
          MINVAL(oas_rcv_field(:,oas_i)), MAXVAL(oas_rcv_field(:,oas_i))
        CALL message(routine, oas_message)
      END IF
        END DO
              
#endif

      ! Compute diagnostic physics fields
      CALL aggr_landvars
      ! Initial call of (slow) physics schemes, including computation of transfer coefficients
      CALL init_slowphysics (mtime_current, 1, dtime)

      DO jg = 1, n_dom

        IF (.NOT. p_patch(jg)%ldom_active) CYCLE

        IF(.NOT.atm_phy_nwp_config(jg)%is_les_phy) THEN

          ! diagnostics which are only required for output
          CALL nwp_diag_for_output(mtime_current, kstart_moist(jg),           & !in
               &                      ih_clch(jg), ih_clcm(jg),               & !in
               &                      phy_params(jg),                         & !in
               &                      p_patch(jg),                            & !in
               &                      p_nh_state(jg)%metrics,                 & !in
               &                      p_nh_state(jg)%prog(nnow(jg)),          & !in  !nnow or nnew?
               &                      p_nh_state(jg)%prog(nnow_rcf(jg)),      & !in  !nnow or nnew?
               &                      p_nh_state(jg)%diag,                    & !in
               &                      p_lnd_state(jg)%diag_lnd,               & !in
               &                      p_lnd_state(jg)%prog_lnd(nnow_rcf(jg)), & !in
               &                      p_lnd_state(jg)%prog_wtr(nnow_rcf(jg)), & !inout
               &                      ext_data(jg),                           & !in
               &                      prm_diag(jg)                            ) !inout

          ! In case of vertical nesting, copy upper levels of synsat input fields to local parent grid
          DO jn = 1, p_patch(jg)%n_childdom
            jgc = p_patch(jg)%child_id(jn)
            IF (.NOT. p_patch(jgc)%ldom_active) CYCLE
            IF (lsynsat(jgc) .AND. p_patch(jgc)%nshift > 0) CALL copy_rttov_ubc (jg, jgc)
          ENDDO
          ! Compute synthetic sat images
          IF (lsynsat(jg)) CALL rttov_driver (jg, p_patch(jg)%parent_id, nnow_rcf(jg))

        ELSE !is_les_phy

           !LES specific diagnostics only for output
           CALL les_cloud_diag    ( kstart_moist(jg),                       & !in
             &                      ih_clch(jg), ih_clcm(jg),               & !in
             &                      phy_params(jg),                         & !in
             &                      p_patch(jg),                            & !in
             &                      p_nh_state(jg)%metrics,                 & !in
             &                      p_nh_state(jg)%prog(nnow(jg)),          & !in  !nnow or nnew?
             &                      p_nh_state(jg)%prog(nnow_rcf(jg)),      & !in  !nnow or nnew?
             &                      p_nh_state(jg)%diag,                    & !in
             &                      p_lnd_state(jg)%diag_lnd,               & !in
             &                      p_lnd_state(jg)%prog_lnd(nnow_rcf(jg)), & !in
             &                      prm_diag(jg),                           & !inout
             &                      ext_data(jg)                            ) !in

         END IF!is_les_phy
      ENDDO!jg

      CALL fill_nestlatbc_phys

    ENDIF!is_restart
  CASE (iecham)
    IF (.NOT.isRestart()) THEN
      CALL init_slowphysics (mtime_current, 1, dtime)
    END IF
  END SELECT ! iforcing

  !------------------------------------------------------------------
  !  get and write out some of the initial values
  !------------------------------------------------------------------
  IF (.NOT.isRestart() .AND. (mtime_current >= time_config%tc_exp_startdate)) THEN

    !--------------------------------------------------------------------------
    ! loop over the list of internal post-processing tasks, e.g.
    ! interpolate selected fields to p- and/or z-levels
    simulation_status = new_simulation_status(l_first_step   = .TRUE.,                  &
      &                                       l_output_step  = .TRUE.,                  &
      &                                       l_dom_active   = p_patch(1:)%ldom_active, &
      &                                       i_timelevel_dyn= nnow, i_timelevel_phy= nnow_rcf)
    CALL pp_scheduler_process(simulation_status)

    CALL perform_accumulation(nnow(1),nnow_rcf(1))
    IF (p_nh_opt_diag(1)%acc%l_any_m) THEN
      CALL update_opt_acc(p_nh_opt_diag(1)%acc,            &
        &                 p_nh_state(1)%prog(nnow_rcf(1)), &
        &                 p_nh_state(1)%prog(nnow(1))%rho, &
        &                 p_nh_state(1)%diag,              &
        &                 p_patch(1)%cells%owned,          &
        &                 p_patch(1)%nlev                  )
    END IF

    IF (output_mode%l_nml) THEN
      CALL write_name_list_output(jstep=0)
    END IF

    IF (p_nh_opt_diag(1)%acc%l_any_m) THEN
      CALL reset_opt_acc(p_nh_opt_diag(1)%acc)
    END IF
    CALL reset_accumulation

    ! sample meteogram output
    DO jg = 1, n_dom
      IF (.NOT. output_mode%l_none .AND. &    ! meteogram output is not initialized for output=none
        & p_patch(jg)%ldom_active  .AND. &
        & meteogram_is_sample_step( meteogram_output_config(jg), 0 ) ) THEN
        CALL meteogram_sample_vars(jg, 0, time_config%tc_startdate)
      END IF
    END DO

    !AD: Also output special diagnostics for LES on torus
    IF(atm_phy_nwp_config(1)%is_les_phy .AND. sampl_freq_step>0)THEN
      CALL calculate_turbulent_diagnostics(                      &
                             & p_patch(1),                       & !in
                             & p_nh_state(1)%prog(nnow(1)),      &
                             & p_nh_state(1)%prog(nnow_rcf(1)),  & !in
                             & p_nh_state(1)%diag,                   & !in
                             & p_lnd_state(1)%prog_lnd(nnow_rcf(1)), &
                             & p_lnd_state(1)%diag_lnd,              &
                             & prm_nwp_tend(1),                      &
                             & prm_diag(1)                )     !inout

      !write out time series
      CALL write_time_series(prm_diag(1)%turb_diag_0dvar, mtime_current)
      CALL write_vertical_profiles(prm_diag(1)%turb_diag_1dvar, mtime_current, 1)
      prm_diag(1)%turb_diag_1dvar = 0._wp
    END IF


#ifdef MESSY
    ! MESSy initial output
!    CALL messy_write_output
#endif

  END IF ! not isRestart()

  IF (timers_level > 3) CALL timer_stop(timer_model_init)

!   IF (parallel_radiation_omp) THEN
!
!     !---------------------------------------
!     CALL init_ompthread_radiation()
!
! !$    CALL omp_set_nested(.TRUE.)
! !$    CALL omp_set_num_threads(2)
! !$    write(0,*) 'omp_get_max_active_levels=',omp_get_max_active_levels
! !$    write(0,*) 'omp_get_max_threads=',omp_get_max_threads()
! !$OMP PARALLEL SECTIONS
! !$OMP SECTION
! !$  CALL omp_set_num_threads(nh_stepping_ompthreads)
! !$    write(0,*) 'This is the nh_timeloop, max threads=',omp_get_max_threads()
! !$    write(0,*) 'omp_get_num_threads=',omp_get_num_threads()
!
!     CALL perform_nh_timeloop (datetime_current, jfile, l_have_output )
!     CALL model_end_ompthread()
!
! !$OMP SECTION
! !$  write(0,*) 'This is the nwp_parallel_radiation_thread, max threads=',&
! !$    omp_get_max_threads()
!   CALL nwp_start_radiation_ompthread()
! !$OMP END PARALLEL SECTIONS
!
!   ELSE
    !---------------------------------------

    CALL perform_nh_timeloop (mtime_current, latbc)
!   ENDIF

  CALL deallocate_nh_stepping (latbc)


  END SUBROUTINE perform_nh_stepping
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !! Organizes nonhydrostatic time stepping
  !! Currently we assume to have only one grid level.
  !!
  !! @par Revision History
  !! Initial release by Almut Gassmann, (2009-04-15)
  !!
  SUBROUTINE perform_nh_timeloop (mtime_current, latbc)
    !
    CHARACTER(len=*), PARAMETER :: routine = modname//':perform_nh_timeloop'
    TYPE(t_latbc_data),     INTENT(INOUT)  :: latbc !< data structure for async latbc prefetching
    TYPE(datetime),         POINTER        :: mtime_current     ! current datetime (mtime)

  INTEGER                              :: jg, jn, jgc
  INTEGER                              :: ierr
  LOGICAL                              :: l_compute_diagnostic_quants,  &
    &                                     l_nml_output, lprint_timestep, &
    &                                     lwrite_checkpoint, lcfl_watch_mode
  TYPE(t_simulation_status)            :: simulation_status
  TYPE(datetime),   POINTER            :: mtime_old         ! copy of current datetime (mtime)

  INTEGER                              :: i, iau_iter
  REAL(wp)                             :: elapsed_time_global
  INTEGER                              :: jstep   ! step number
  INTEGER                              :: jstep0  ! step for which the restart file
                                                  ! was produced
  INTEGER                              :: kstep   ! step number relative to restart step
  INTEGER                              :: jstep_shift ! start counter for time loop
  INTEGER, ALLOCATABLE                 :: output_jfile(:)

  TYPE(timedelta), POINTER             :: forecast_delta
  CHARACTER(LEN=128)                   :: forecast_delta_str

  TYPE(timedelta), POINTER             :: model_time_step => NULL()

  TYPE(datetime), POINTER              :: eventRefDate      => NULL(), &
       &                                  eventStartDate    => NULL(), &
       &                                  eventEndDate      => NULL()
  TYPE(datetime), POINTER              :: checkpointRefDate => NULL(), &
       &                                  restartRefDate    => NULL()
  TYPE(timedelta), POINTER             :: eventInterval     => NULL()
  TYPE(event), POINTER                 :: checkpointEvent   => NULL()
  TYPE(event), POINTER                 :: restartEvent      => NULL()

  INTEGER                              :: checkpointEvents
  LOGICAL                              :: lret
  TYPE(t_datetime_ptr)                 :: datetime_current(max_dom) 
  TYPE(timeDelta), POINTER             :: time_diff
  TYPE(t_RestartAttributeList), POINTER :: restartAttributes
  CLASS(t_RestartDescriptor), POINTER  :: restartDescriptor

  CHARACTER(LEN=MAX_TIMEDELTA_STR_LEN)   :: td_string
  CHARACTER(LEN=MAX_DATETIME_STR_LEN)    :: dt_string, dstring
  CHARACTER(len=MAX_MTIME_ERROR_STR_LEN) :: errstring
  
  REAL(wp)                             :: sim_time         !< elapsed simulation time
  INTEGER                              :: sim_time_oas     !< elapsed simulation time
  CHARACTER(len=128)                   :: oas_message
#ifdef COUP_OAS_ICON
  INTEGER :: jc, jb, nlev, rl_start, rl_end, i_startblk, &
    i_endblk, i_startidx, i_endidx, ii
#endif

  LOGICAL :: l_isStartdate, l_isExpStopdate, l_isRestart, l_isCheckpoint, l_doWriteRestart
  
!!$  INTEGER omp_get_num_threads

!-----------------------------------------------------------------------

  IF (ltimer) CALL timer_start(timer_total)
  
  ! calculate elapsed simulation time in seconds
  time_diff  => newTimedelta("PT0S")
  time_diff  =  getTimeDeltaFromDateTime(mtime_current, time_config%tc_exp_startdate)
  sim_time   =  getTotalMillisecondsTimedelta(time_diff, mtime_current)*1.e-3_wp

  IF (iterate_iau .AND. .NOT. isRestart()) THEN
    iau_iter = 1
  ELSE
    iau_iter = 0
  ENDIF

  ! allocate temporary variable for restarting purposes
  IF (output_mode%l_nml) THEN
    ALLOCATE(output_jfile(SIZE(output_file)), STAT=ierr)
    IF (ierr /= SUCCESS)  CALL finish (routine, 'ALLOCATE failed!')
  ENDIF
  
  ! If the testbed mode is selected, reset iorder_sendrecv to 0 in order to suppress
  ! MPI communication from now on.
  IF (test_mode > 0) iorder_sendrecv = 0

  IF (timeshift%dt_shift < 0._wp  .AND. .NOT. isRestart()) THEN
    jstep_shift = NINT(timeshift%dt_shift/dtime)
    WRITE(message_text,'(a,i6,a)') 'Model start shifted backwards by ', ABS(jstep_shift),' time steps'
    CALL message(TRIM(routine),message_text)
    atm_phy_nwp_config(:)%lcalc_acc_avg = .FALSE.
  ELSE
    jstep_shift = 0
  ENDIF

  mtime_old => newDatetime(mtime_current)
  DO jg=1, n_dom
    datetime_current(jg)%ptr => newDatetime(mtime_current)
  END DO

  restartDescriptor => createRestartDescriptor("atm")

  jstep0 = 0

  restartAttributes => getAttributesForRestarting()
  IF (isRestart()) THEN
    ! get start counter for time loop from restart file:
    jstep0 = restartAttributes%getInteger("jstep")
  END IF

  ! for debug purposes print var lists: for msg_level >= 13 short and for >= 20 long format
  IF  (.NOT. ltestcase) THEN
    IF (msg_level >= 13) THEN
      DO i = 1, nvar_lists
        IF (msg_level >= 20) THEN
          CALL print_var_list(var_lists(i), lshort=.FALSE.)
        ELSE
          CALL print_var_list(var_lists(i), lshort=.TRUE.)
        ENDIF
      ENDDO
    ENDIF
  ENDIF

  ! Check if current number of dynamics substeps is larger than the default value
  ! (this can happen for restarted runs only at this point)
  IF (ANY(ndyn_substeps_var(1:n_dom) > ndyn_substeps)) THEN
    lcfl_watch_mode = .TRUE.
  ELSE
    lcfl_watch_mode = .FALSE.
  ENDIF

  ! init routine for mo_nh_supervise module (eg. opening of files)
  CALL init_supervise_nh()

  ! set events, group and the events

  CALL message('','')

  eventRefDate   => time_config%tc_exp_refdate
  eventStartDate => time_config%tc_exp_startdate
  eventEndDate   => time_config%tc_exp_stopdate

  ! for debugging purposes the referenece (anchor) date for checkpoint
  ! and restart may be switched to be relative to current jobs start
  ! date instead of the experiments start date.

  IF (time_config%is_relative_time) THEN
    checkpointRefDate => time_config%tc_startdate
    restartRefDate    => time_config%tc_startdate
  ELSE
    checkpointRefDate => time_config%tc_exp_startdate
    restartRefDate    => time_config%tc_exp_startdate
  ENDIF

  ! create an event manager, ie. a collection of different events
  CALL initEventManager(time_config%tc_exp_refdate)

  ! --- create an event group for checkpointing and restart
  checkpointEvents =  addEventGroup('checkpointEventGroup')
  checkpointEventGroup => getEventGroup(checkpointEvents)

  ! --- --- create checkpointing event
  eventInterval  => time_config%tc_dt_checkpoint
  checkpointEvent => newEvent('checkpoint', checkpointRefDate, eventStartDate, eventEndDate, eventInterval, errno=ierr)
  IF (ierr /= no_Error) THEN
    ! give an elaborate error message:
    CALL datetimeToString(checkpointRefDate, dt_string)
    WRITE (0,*) "event reference date: ",    dt_string
    CALL datetimeToString(eventStartDate,    dt_string)
    WRITE (0,*) "event start date    : ",    dt_string
    CALL datetimeToString(eventEndDate,      dt_string)
    WRITE (0,*) "event end date      : ",    dt_string
    CALL timedeltaToString(eventInterval,    td_string)
    WRITE (0,*) "event interval      : ",    td_string
    CALL mtime_strerror(ierr, errstring)
    CALL finish('perform_nh_timeloop', "event 'checkpoint': "//errstring)
  ENDIF
  lret = addEventToEventGroup(checkpointEvent, checkpointEventGroup)

  ! --- --- create restart event, ie. checkpoint + model stop
  eventInterval  => time_config%tc_dt_restart
  restartEvent => newEvent('restart', restartRefDate, eventStartDate, eventEndDate, eventInterval, errno=ierr)
  IF (ierr /= no_Error) THEN
    ! give an elaborate error message:
    CALL datetimeToString(restartRefDate, dt_string)
    WRITE (0,*) "event reference date: ", dt_string
    CALL datetimeToString(eventStartDate, dt_string)
    WRITE (0,*) "event start date    : ", dt_string
    CALL datetimeToString(eventEndDate,   dt_string)
    WRITE (0,*) "event end date      : ", dt_string
    CALL timedeltaToString(eventInterval, td_string)
    WRITE (0,*) "event interval      : ", td_string
    CALL mtime_strerror(ierr, errstring)
    CALL finish('perform_nh_timeloop', "event 'restart': "//errstring)
  ENDIF
  lret = addEventToEventGroup(restartEvent, checkpointEventGroup)

  CALL printEventGroup(checkpointEvents)

  ! set time loop properties
  model_time_step => time_config%tc_dt_model
  
  CALL message('','')
  CALL datetimeToString(mtime_current, dstring)
  WRITE(message_text,'(a,a)') 'Start date of this run: ', dstring
  CALL message('',message_text)
  CALL datetimeToString(time_config%tc_stopdate, dstring)
  WRITE(message_text,'(a,a)') 'Stop date of this run:  ', dstring
  CALL message('',message_text)
  CALL message('','')

  jstep = jstep0+jstep_shift+1

#if defined( _OPENACC )
  i_am_accel_node = my_process_is_work()    ! Activate GPUs

  CALL save_convenience_pointers( )

!$ACC DATA COPYIN( p_int_state, p_patch, p_nh_state, prep_adv, advection_config ), IF ( i_am_accel_node )

  CALL refresh_convenience_pointers( )
#endif

  TIME_LOOP: DO

    ! Check if a nested domain needs to be turned off
    DO jg=2, n_dom
      IF (p_patch(jg)%ldom_active .AND. (sim_time >= end_time(jg))) THEN
        p_patch(jg)%ldom_active = .FALSE.
        WRITE(message_text,'(a,i2,a,f12.2)') 'domain ',jg,' stopped at time ',sim_time
        CALL message('perform_nh_timeloop', TRIM(message_text))
      ENDIF
    ENDDO

#ifdef MESSYTIMER
    CALL messy_time(1)
#endif

    ! update model date and time mtime based
    mtime_current = mtime_current + model_time_step

    ! store state of output files for restarting purposes
    IF (output_mode%l_nml .AND. jstep>=0 ) THEN
      DO i=1,SIZE(output_file)
        output_jfile(i) = get_current_jfile(output_file(i)%out_event)
      END DO
    ENDIF

    ! turn on calculation of averaged and accumulated quantities at the first regular time step
    IF (jstep-jstep0 == 1) atm_phy_nwp_config(:)%lcalc_acc_avg = .TRUE.

    ! read boundary data if necessary
    IF (l_limited_area .AND. latbc_config%itype_latbc > 0 .AND. num_prefetch_proc == 0) &
      CALL read_latbc_data_sync(p_patch(1), p_nh_state(1), ext_data(1), p_int_state(1), mtime_current)

    IF (msg_level > 2) THEN
      lprint_timestep = .TRUE.
    ELSE
      lprint_timestep = MOD(jstep,25) == 0
    ENDIF

    ! always print the first and the last time step
    lprint_timestep = lprint_timestep .OR. (jstep == jstep0+1) .OR. (jstep == jstep0+nsteps)

    IF (lprint_timestep) THEN
      ! compute current forecast time (delta):
      forecast_delta => newTimedelta("P01D")
      forecast_delta = mtime_current - time_config%tc_exp_startdate
      ! we append the forecast time delta as an ISO 8601 conforming
      ! string (where, for convenience, the 'T' token has been
      ! replaced by a blank character)
      IF (forecast_delta%ms /= 0) THEN
        WRITE (forecast_delta_str,'(i0,a,3(i2.2,a),i3.3,a)') &
             &                                   forecast_delta%day, 'D ',   &
             &                                   forecast_delta%hour, 'H',   &
             &                                   forecast_delta%minute, 'M', &
             &                                   forecast_delta%second, '.', &
             &                                   forecast_delta%ms, 'S'
      ELSE
        WRITE (forecast_delta_str,'(i0,a,3(i2.2,a))') &
             &                                   forecast_delta%day, 'D ',   &
             &                                   forecast_delta%hour, 'H',   &
             &                                   forecast_delta%minute, 'M', &
             &                                   forecast_delta%second, 'S'
      ENDIF

      IF (iforcing == inwp) THEN
        WRITE(message_text,'(a,i8,a,i0,a,5(i2.2,a),i3.3,a,a)') &
             &             'Time step: ', jstep, ', model time: ',                              &
             &             mtime_current%date%year,   '-', mtime_current%date%month,    '-',    &
             &             mtime_current%date%day,    ' ', mtime_current%time%hour,     ':',    &
             &             mtime_current%time%minute, ':', mtime_current%time%second,   '.',    &
             &             mtime_current%time%ms, ' forecast time ', TRIM(forecast_delta_str)
      ELSE
        WRITE(message_text,'(a,i8,a,i0,a,4(i2.2,a),i2.2)') &
             &             'Time step: ', jstep, ' model time ',                                &
             &             mtime_current%date%year,   '-', mtime_current%date%month,    '-',    &
             &             mtime_current%date%day,    ' ', mtime_current%time%hour,     ':',    &
             &             mtime_current%time%minute, ':', mtime_current%time%second
      ENDIF

      CALL message('',message_text)

      CALL deallocateTimedelta(forecast_delta)
    ENDIF

    ! Update the following surface fields, if a new day is coming
    !
    ! - ndviratio, plcov_t, tai_t, sai_t
    ! - SST, fr_seaice (depending on sstice_mode)
    ! - MODIS albedo fields alb_dif, albuv_dif, albni_dif
    !

    ! The update is skipped in IAU iteration mode if the model is reset to the initial state at the
    ! end of the current time step
    IF ( (mtime_current%date%day /= mtime_old%date%day) .AND. .NOT. (jstep == 0 .AND. iau_iter == 1) ) THEN

      WRITE(message_text,'(a,i10,a,i10)') 'New day  day_old: ', mtime_old%date%day, &
                &                 ' ,  day: ', mtime_current%date%day
      CALL message(TRIM(routine),message_text)

      !Update ndvi normalized differential vegetation index
      IF (itopo == 1 .AND. iforcing == inwp .AND.                  &
        & ALL(atm_phy_nwp_config(1:n_dom)%inwp_surface >= 1)) THEN
        DO jg=1, n_dom
          ! Note: here we can assume that hour=minute=second=0, due to
          ! the surrounding if-clause.
          CALL interpol_monthly_mean(p_patch(jg), mtime_current,            &! in
            &                        ext_data(jg)%atm_td%ndvi_mrat,         &! in
            &                        ext_data(jg)%atm%ndviratio      )! out
        ENDDO

        ! after updating ndvi_mrat, probably plcov_t and tai_t have to be updated also.
        ! So it is better not to update ndvi_mrat till this is clarified
        CALL update_ndvi(p_patch(1:), ext_data)
      END IF

      ! Check if the SST and Sea ice fraction have to be updated (sstice_mode 3,4,5)
      IF ( ANY((/SSTICE_ANA_CLINC, SSTICE_CLIM,SSTICE_AVG_MONTHLY,SSTICE_AVG_DAILY/) == sstice_mode) &
        &   .AND. iforcing == inwp  ) THEN

        CALL set_actual_td_ext_data (.FALSE., mtime_current, mtime_old, sstice_mode,  &
                                  &  p_patch(1:), ext_data, p_lnd_state)

        CALL update_sstice( p_patch(1:), ext_data, p_lnd_state,                    &
          &                 p_nh_state, sstice_mode, time_config%tc_exp_startdate, &
          &                 mtime_current )

      END IF  !sstice_mode>1


      ! Check if MODIS albedo needs to be updated
      IF (iforcing == inwp .AND. albedo_type == MODIS) THEN
        ! Note that here only an update of the external parameter fields is
        ! performed. The actual update happens in mo_albedo.
        !
        ! Note: here we can assume that hour=minute=second=0, due to
        ! the surrounding if-clause.
        DO jg = 1, n_dom
          CALL interpol_monthly_mean(p_patch(jg), mtime_current,            &! in
            &                        ext_data(jg)%atm_td%alb_dif,           &! in
            &                        ext_data(jg)%atm%alb_dif          )! out

          CALL interpol_monthly_mean(p_patch(jg), mtime_current,            &! in
            &                        ext_data(jg)%atm_td%albuv_dif,         &! in
            &                        ext_data(jg)%atm%albuv_dif        )! out

          CALL interpol_monthly_mean(p_patch(jg), mtime_current,            &! in
            &                        ext_data(jg)%atm_td%albni_dif,         &! in
            &                        ext_data(jg)%atm%albni_dif        )! out
        ENDDO
      ENDIF

      mtime_old = mtime_current

    END IF ! end update of surface parameter fields


    !--------------------------------------------------------------------------
    ! Set output flags
    !--------------------------------------------------------------------------

    l_nml_output = output_mode%l_nml .AND. jstep >= 0 .AND. istime4name_list_output(jstep)

    ! In IAU iteration mode, output at the nominal initial date is written only at the
    ! end of the first cycle, providing an initialized analysis to which the analysis 
    ! increments have been completely added
    IF (jstep == 0 .AND. iau_iter == 2) l_nml_output = .FALSE.

    ! Computation of diagnostic quantities may also be necessary for
    ! meteogram sampling:
!DR Note that this may be incorrect for meteograms in case that
!DR meteogram_output_config is not the same for all domains.
    l_compute_diagnostic_quants = l_nml_output
    DO jg = 1, n_dom
      l_compute_diagnostic_quants = l_compute_diagnostic_quants .OR. &
        &          meteogram_is_sample_step(meteogram_output_config(jg), jstep )
    END DO
    l_compute_diagnostic_quants = jstep >= 0 .AND. l_compute_diagnostic_quants .AND. &
      &                           .NOT. output_mode%l_none


    ! Calculations for enhanced sound-wave and gravity-wave damping during the spinup phase
    ! if mixed second-order/fourth-order divergence damping (divdamp_order=24) is chosen.
    ! Includes increased vertical wind off-centering during the first 2 hours of integration.
    IF (divdamp_order==24 .AND. .NOT. isRestart()) THEN
      elapsed_time_global = (REAL(jstep,wp)-0.5_wp)*dtime
      IF (elapsed_time_global <= 7200._wp+0.5_wp*dtime .AND. .NOT. ltestcase) THEN
        CALL update_spinup_damping(elapsed_time_global)
      ENDIF
    ELSE IF (divdamp_order==24) THEN
      divdamp_fac_o2 = 0._wp
    ENDIF


    !--------------------------------------------------------------------------
    !
    ! dynamics stepping
    !
    CALL integrate_nh(datetime_current, 1, jstep-jstep_shift, iau_iter, dtime, model_time_step, 1, latbc)

    ! Compute diagnostics for output if necessary
    IF (l_compute_diagnostic_quants .OR. iforcing==iecham .OR. iforcing==inoforcing) THEN

      CALL diag_for_output_dyn ()

      IF (iforcing == inwp) THEN
        CALL aggr_landvars

        DO jg = 1, n_dom
          IF (.NOT. p_patch(jg)%ldom_active) CYCLE

          IF(.NOT.atm_phy_nwp_config(jg)%is_les_phy) THEN
            ! diagnostics which are only required for output
            CALL nwp_diag_for_output(mtime_current, kstart_moist(jg),           & !in
                 &                      ih_clch(jg), ih_clcm(jg),               & !in
                 &                      phy_params(jg),                         & !in
                 &                      p_patch(jg),                            & !in
                 &                      p_nh_state(jg)%metrics,                 & !in
                 &                      p_nh_state(jg)%prog(nnow(jg)),          & !in  !nnow or nnew?
                 &                      p_nh_state(jg)%prog(nnow_rcf(jg)),      & !in  !nnow or nnew?
                 &                      p_nh_state(jg)%diag,                    & !in
                 &                      p_lnd_state(jg)%diag_lnd,               & !in
                 &                      p_lnd_state(jg)%prog_lnd(nnow_rcf(jg)), & !in
                 &                      p_lnd_state(jg)%prog_wtr(nnow_rcf(jg)), & !inout
                 &                      ext_data(jg),                           & !in
                 &                      prm_diag(jg)                            ) !inout

            ! In case of vertical nesting, copy upper levels of synsat input fields to local parent grid
            DO jn = 1, p_patch(jg)%n_childdom
              jgc = p_patch(jg)%child_id(jn)
              IF (.NOT. p_patch(jgc)%ldom_active) CYCLE
              IF (lsynsat(jgc) .AND. p_patch(jgc)%nshift > 0) CALL copy_rttov_ubc (jg, jgc)
            ENDDO
            ! Compute synthetic sat images
            IF (lsynsat(jg)) CALL rttov_driver (jg, p_patch(jg)%parent_id, nnow_rcf(jg))

          ELSE !is_les_phy

            !LES specific diagnostics only for output
            CALL les_cloud_diag    ( kstart_moist(jg),                       & !in
              &                      ih_clch(jg), ih_clcm(jg),               & !in
              &                      phy_params(jg),                         & !in
              &                      p_patch(jg),                            & !in
              &                      p_nh_state(jg)%metrics,                 & !in
              &                      p_nh_state(jg)%prog(nnow(jg)),          & !in  !nnow or nnew?
              &                      p_nh_state(jg)%prog(nnow_rcf(jg)),      & !in  !nnow or nnew?
              &                      p_nh_state(jg)%diag,                    & !in
              &                      p_lnd_state(jg)%diag_lnd,               & !in
              &                      p_lnd_state(jg)%prog_lnd(nnow_rcf(jg)), & !in
              &                      prm_diag(jg),                           & !inout
              &                      ext_data(jg)                            ) !in

          END IF!is_les_phy

        ENDDO!jg

        CALL fill_nestlatbc_phys

      END IF !iforcing=inwp

      ! Unit conversion for output from mass mixing ratios to densities
      ! and calculation of ART diagnostics
      DO jg = 1, n_dom
        IF (.NOT. p_patch(jg)%ldom_active) CYCLE
        ! Call the ART diagnostics
        CALL art_diagnostics_interface(p_patch(jg),                              &
          &                            p_nh_state(jg)%prog(nnew(jg))%rho,        &
          &                            p_nh_state(jg)%diag%pres,                 &
          &                            p_nh_state(jg)%prog(nnow_rcf(jg))%tracer, &
          &                            p_nh_state(jg)%metrics%ddqz_z_full,       &
          &                            p_nh_state(jg)%metrics%z_mc, jg)
        ! Call the ART unit conversion 
        CALL art_tools_interface('unit_conversion',                            & !< in
          &                      p_nh_state_lists(jg)%prog_list(nnow_rcf(jg)), & !< in
          &                      p_nh_state(jg)%prog(nnow_rcf(jg))%tracer,     & !< in
          &                      p_nh_state(jg)%prog(nnew_rcf(jg))%tracer,     & !< out
          &                      p_nh_state(jg)%prog(nnew(jg))%rho)              !< in
      END DO

    ENDIF


    ! Adapt number of dynamics substeps if necessary
    !
    IF (lcfl_watch_mode .OR. MOD(jstep-jstep_shift,5) == 0) THEN
      CALL set_ndyn_substeps(lcfl_watch_mode)
    ENDIF

    !--------------------------------------------------------------------------
    ! loop over the list of internal post-processing tasks, e.g.
    ! interpolate selected fields to p- and/or z-levels
    !
    ! Mean sea level pressure needs to be computed also at
    ! no-output-steps for accumulation purposes; set by l_accumulation_step
    simulation_status = new_simulation_status(l_output_step  = l_nml_output,             &
      &                                       l_last_step    = (jstep==(nsteps+jstep0)), &
      &                                       l_accumulation_step = (iforcing == iecham),&
      &                                       l_dom_active   = p_patch(1:)%ldom_active,  &
      &                                       i_timelevel_dyn= nnow, i_timelevel_phy= nnow_rcf)
    CALL pp_scheduler_process(simulation_status)

#ifdef MESSY
    IF (timers_level > 3) CALL timer_start(timer_extra1)
    CALL messy_write_output
    IF (timers_level > 3) CALL timer_stop(timer_extra1)
#endif

    ! update accumlated values
    CALL perform_accumulation(nnow(1),nnow_rcf(1))
    IF (p_nh_opt_diag(1)%acc%l_any_m) THEN
      CALL update_opt_acc(p_nh_opt_diag(1)%acc,            &
        &                 p_nh_state(1)%prog(nnow_rcf(1)), &
        &                 p_nh_state(1)%prog(nnow(1))%rho, &
        &                 p_nh_state(1)%diag,              &
        &                 p_patch(1)%cells%owned,          &
        &                 p_patch(1)%nlev)
    IF (l_nml_output) CALL calc_mean_opt_acc(p_nh_opt_diag(1)%acc)
    END IF

    ! output of results
    ! note: nnew has been replaced by nnow here because the update
    IF (l_nml_output) THEN
      CALL write_name_list_output(jstep)
    ENDIF

    CALL reset_accumulation


    ! sample meteogram output
    DO jg = 1, n_dom
      IF (.NOT. output_mode%l_none .AND. &    ! meteogram output is not initialized for output=none
        & p_patch(jg)%ldom_active  .AND. .NOT. (jstep == 0 .AND. iau_iter == 2) .AND. &
        & meteogram_is_sample_step(meteogram_output_config(jg), jstep)) THEN
        CALL meteogram_sample_vars(jg, jstep, mtime_current)
      END IF
    END DO



    ! Diagnostics: computation of total integrals
    !
    ! Diagnostics computation is not yet properly MPI-parallelized
    !
    IF (output_mode%l_totint .AND. is_totint_time(current_step =jstep,   &
      &                                           restart_step = jstep0, &
      &                                           n_diag       = n_diag, &
      &                                           n_steps      = nsteps) ) THEN

      kstep = jstep-jstep0

#ifdef NOMPI
      IF (my_process_is_mpi_all_seq()) &
#endif
        CALL supervise_total_integrals_nh( kstep, p_patch(1:), p_nh_state, p_int_state(1:), &
        &                                  nnow(1:n_dom), nnow_rcf(1:n_dom), jstep == (nsteps+jstep0))
    ENDIF


    ! re-initialize MAX/MIN fields with 'resetval'
    ! must be done AFTER output

    CALL reset_act%execute(slack=dtime, mtime_date=mtime_current)

    ! re-initialization for FG-averaging. Ensures that average is centered in time.
    IF (is_avgFG_time(mtime_current)) THEN
      IF (p_nh_state(1)%diag%nsteps_avg(1) == 0) THEN
        CALL reinit_average_first_guess(p_patch(1), p_nh_state(1)%diag, p_nh_state(1)%prog(nnow_rcf(1)))
      END IF
    ENDIF


    !--------------------------------------------------------------------------
    ! Write restart file
    !--------------------------------------------------------------------------
    ! check whether time has come for writing restart file

    CALL message('','')
    !
    ! default is to assume we do not write a checkpoint/restart file
    lwrite_checkpoint = .FALSE.
    ! if thwe model is not supposed to write output, do not write checkpoints
    IF (.NOT. output_mode%l_none ) THEN
      ! to clarify the decision tree we use shorter and more expressive names:

      l_isStartdate    = (time_config%tc_startdate == mtime_current)
      l_isExpStopdate  = (time_config%tc_exp_stopdate == mtime_current)
      l_isRestart      = isCurrentEventActive(restartEvent, mtime_current)
      l_isCheckpoint   = isCurrentEventActive(checkpointEvent, mtime_current)
      l_doWriteRestart = time_config%tc_write_restart

      IF ( &
           !  if normal checkpoint or restart cycle has been reached, i.e. checkpoint+model stop
           &         (l_isRestart .OR. l_isCheckpoint)                     &
           &  .AND.                                                        &
           !  and the current date differs from the start date
           &        .NOT. l_isStartdate                                    &
           &  .AND.                                                        &
           !  and end of run has not been reached or restart writing has been disabled
           &        (.NOT. l_isExpStopdate .OR. l_doWriteRestart)          &
           & ) THEN
        lwrite_checkpoint = .TRUE.
      END IF
    END IF


    IF (lwrite_checkpoint) THEN

      ! apply nest boundary filling; this has no impact on the correctness (in the sense of reproducibility)
      ! of the restart but facilitates debugging
      CALL diag_for_output_dyn ()
      IF (iforcing == inwp) THEN
        CALL aggr_landvars
        CALL fill_nestlatbc_phys
      END IF

        DO jg = 1, n_dom
            CALL restartDescriptor%updatePatch(p_patch(jg), &
              & opt_t_elapsed_phy          = t_elapsed_phy(jg,:),        &
              & opt_lcall_phy              = lcall_phy(jg,:),            &
              & opt_ndyn_substeps          = ndyn_substeps_var(jg),      &
              & opt_jstep_adv_marchuk_order= jstep_adv(jg)%marchuk_order,&
              & opt_depth_lnd              = nlev_soil,                  &
              & opt_nlev_snow              = nlev_snow,                  &
              & opt_ndom                   = n_dom)
        ENDDO

        CALL restartDescriptor%writeRestart(mtime_current, jstep, opt_output_jfile = output_jfile)

#ifdef MESSY
        IF (timers_level > 3) CALL timer_start(timer_extra2)
        CALL messy_write_restart
        IF (timers_level > 3) CALL timer_stop(timer_extra2)
#endif
    END IF  ! lwrite_checkpoint

#ifdef MESSYTIMER
    ! model time management synchronisation
    CALL messy_time(2)
#endif

    ! prefetch boundary data if necessary
    IF(num_prefetch_proc >= 1 .AND. latbc_config%itype_latbc > 0) THEN
      CALL recv_latbc_data(latbc, p_patch(1), p_nh_state(1), p_int_state(1), &
        &                  cur_datetime=mtime_current, lcheck_read=.TRUE., ltime_incr=.TRUE.,     &
        &                  tlev=latbc%new_latbc_tlev)
    ENDIF

    IF (mtime_current >= time_config%tc_stopdate) THEN
       ! leave time loop
       EXIT TIME_LOOP
    END IF

    ! Reset model to initial state if IAU iteration is selected and the first iteration cycle has been completed

    IF (jstep == 0 .AND. iau_iter == 1) THEN
      jstep_adv(:)%marchuk_order = 0
      t_elapsed_phy(:,:)         = 0._wp
      linit_dyn(:)               = .TRUE.
      !time_config%sim_time(:)    = timeshift%dt_shift
      mtime_current = mtime_current + timeshift%mtime_shift
      mtime_old = mtime_current
      DO jg = 1, n_dom
        datetime_current(jg)%ptr = mtime_current
      END DO
      CALL reset_to_initial_state(mtime_current)
      iau_iter = 2
      jstep = jstep0+jstep_shift+1
    ELSE
      jstep = jstep + 1
    ENDIF

    time_diff  = getTimeDeltaFromDateTime(mtime_current, time_config%tc_exp_startdate)
    sim_time   = getTotalMillisecondsTimedelta(time_diff, mtime_current)*1.e-3_wp
     
  ENDDO TIME_LOOP

#if defined( _OPENACC )
  CALL save_convenience_pointers( )
!$ACC END DATA
  CALL refresh_convenience_pointers( )
  i_am_accel_node = .FALSE.                 ! Deactivate GPUs
#endif

  ! clean-up routine for mo_nh_supervise module (eg. closing of files)
  CALL finalize_supervise_nh()

  IF (use_async_restart_output) CALL deleteRestartDescriptor(restartDescriptor)

  IF (ltimer) CALL timer_stop(timer_total)

  ! clean up
  IF (output_mode%l_nml) THEN
    DEALLOCATE(output_jfile, STAT=ierr)
    IF (ierr /= SUCCESS)  CALL finish (routine, 'DEALLOCATE failed!')
  ENDIF

  CALL deallocateTimedelta(time_diff)  
  CALL deallocateDatetime(mtime_old)
  DO jg=1,n_dom
    IF (ASSOCIATED(datetime_current(jg)%ptr)) &
      &  CALL deallocateDatetime(datetime_current(jg)%ptr)
  END DO

  END SUBROUTINE perform_nh_timeloop


  !-------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  !>
  !! integrate_nh
  !!
  !! Performs dynamics time stepping:  Rotational modes (helicity bracket) and
  !! divergent modes (Poisson bracket) are split using Strang splitting.
  !!
  !!
  !! @par Revision History
  !! Initial release by Almut Gassmann, MPI-M (2009-08-25)
  !! Adaptation for grid refinement by Guenther Zaengl, DWD (2010-02-09)
  !! Modification by Daniel Reinert, DWD (2010-04-15)
  !!  - Implementation of tracer transport
  !! Modification by Daniel Reinert, DWD (2010-07-23)
  !!  - optional reduced calling frequency for transport and physics
  !!
  RECURSIVE SUBROUTINE integrate_nh (datetime_local, jg, nstep_global,   &
    &                                iau_iter, dt_loc, mtime_dt_loc, num_steps, latbc )

    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: routine = modname//':integrate_nh'

    TYPE(t_datetime_ptr)    :: datetime_local(:)     !< current datetime in mtime format (for each patch)

    INTEGER , INTENT(IN)    :: jg           !< current grid level
    INTEGER , INTENT(IN)    :: nstep_global !< counter of global time step
    INTEGER , INTENT(IN)    :: num_steps    !< number of time steps to be executed
    INTEGER , INTENT(IN)    :: iau_iter     !< counter for IAU iteration
    REAL(wp), INTENT(IN)    :: dt_loc       !< time step applicable to local grid level
    TYPE(timedelta), POINTER :: mtime_dt_loc !< time step applicable to local grid level (mtime format)
    TYPE(t_latbc_data), INTENT(INOUT) :: latbc

    ! Local variables

    ! Time levels
    INTEGER :: n_now_grf, n_now, n_save
    INTEGER :: n_now_rcf, n_new_rcf         ! accounts for reduced calling frequencies (rcf)

    INTEGER :: jstep, jgp, jgc, jn, ii, jc, jb

    REAL(wp):: dt_sub                ! (advective) timestep for next finer grid level
    TYPE(timedelta), POINTER :: mtime_dt_sub
    REAL(wp):: rdt_loc,  rdtmflx_loc ! inverse time step for local grid level

    LOGICAL :: lnest_active, lcall_rrg, lbdy_nudging

    INTEGER, PARAMETER :: nsteps_nest=2 ! number of time steps executed in nested domain

    TYPE(timeDelta), POINTER             :: time_diff
    REAL(wp)                             :: sim_time !< elapsed simulation time on this grid level
    INTEGER                              :: sim_time_oas, oas_i
    CHARACTER(len=128)                   :: oas_message

#ifdef COUP_OAS_ICON
    INTEGER :: rl_start, rl_end, i_startblk, i_endblk, i_startidx, i_endidx, nlev
#endif

    ! calculate elapsed simulation time in seconds (local time for
    ! this domain!)
    time_diff  => newTimedelta("PT0S")
    time_diff  =  getTimeDeltaFromDateTime(datetime_local(jg)%ptr, time_config%tc_exp_startdate)
    sim_time =  getTotalMillisecondsTimedelta(time_diff, datetime_local(jg)%ptr)*1.e-3_wp
    CALL deallocateTimedelta(time_diff)

    !--------------------------------------------------------------------------
    ! This timer must not be called in nested domain because the model crashes otherwise
    IF (jg == 1 .AND. ltimer) CALL timer_start(timer_integrate_nh)

    ! Determine parent domain ID
    IF ( jg > 1) THEN
      jgp = p_patch(jg)%parent_id
    ELSE IF (n_dom_start == 0) THEN
      jgp = 0
    ELSE
      jgp = 1
    ENDIF

    ! If the limited-area mode is used, save initial state in the coarse domain
    ! The save time level is later on used for boundary relaxation in the case of
    ! fixed boundary conditions.
    ! If time-dependent data from a driving model are provided,
    ! they should be written to the save time level, so that the relaxation routine
    ! automatically does the right thing

    IF (jg == 1 .AND. l_limited_area .AND. linit_dyn(jg)) THEN

      n_save = nsav2(jg)
      n_now = nnow(jg)
#ifndef _OPENACC
!$OMP PARALLEL
#endif
      CALL copy(p_nh_state(jg)%prog(n_now)%vn, &
           p_nh_state(jg)%prog(n_save)%vn)
      CALL copy(p_nh_state(jg)%prog(n_now)%w, &
           p_nh_state(jg)%prog(n_save)%w)
      CALL copy(p_nh_state(jg)%prog(n_now)%rho, &
           p_nh_state(jg)%prog(n_save)%rho)
      CALL copy(p_nh_state(jg)%prog(n_now)%theta_v, &
           p_nh_state(jg)%prog(n_save)%theta_v)
#ifndef _OPENACC
!$OMP END PARALLEL
#endif

    ENDIF

    ! This executes one time step for the global domain and two steps for nested domains
    DO jstep = 1, num_steps


      IF (ifeedback_type == 1 .AND. (jstep == 1) .AND. jg > 1 ) THEN
        ! Save prognostic variables at current timestep to compute
        ! feedback increments (not needed in global domain)
        n_now = nnow(jg)
        n_save = nsav2(jg)
#ifndef _OPENACC
!$OMP PARALLEL
#endif
        CALL copy(p_nh_state(jg)%prog(n_now)%vn, &
             p_nh_state(jg)%prog(n_save)%vn)
        CALL copy(p_nh_state(jg)%prog(n_now)%w, &
             p_nh_state(jg)%prog(n_save)%w)
        CALL copy(p_nh_state(jg)%prog(n_now)%rho, &
             p_nh_state(jg)%prog(n_save)%rho)
        CALL copy(p_nh_state(jg)%prog(n_now)%theta_v, &
             p_nh_state(jg)%prog(n_save)%theta_v)
#ifndef _OPENACC
!$OMP END PARALLEL
#endif
      ENDIF


      ! update several switches which decide upon
      ! - switching order of operators in case of Marchuk-splitting
      !
      ! simplified setting (may be removed lateron)
      jstep_adv(jg)%marchuk_order = jstep_adv(jg)%marchuk_order + 1



      IF ( p_patch(jg)%n_childdom > 0 .AND. ndyn_substeps_var(jg) > 1) THEN

        lbdy_nudging = .FALSE.
        lnest_active = .FALSE.
        DO jn = 1, p_patch(jg)%n_childdom
          jgc = p_patch(jg)%child_id(jn)
          IF (p_patch(jgc)%ldom_active) THEN
            lnest_active = .TRUE.
            IF (.NOT. lfeedback(jgc)) lbdy_nudging = .TRUE.
          ENDIF
        ENDDO

        ! Save prognostic variables at current timestep to compute
        ! interpolation tendencies
        n_now  = nnow(jg)
        n_save = nsav1(jg)
        IF (lbdy_nudging) THEN ! full copy needed
#ifndef _OPENACC
!$OMP PARALLEL
#endif
          CALL copy(p_nh_state(jg)%prog(n_now)%vn,p_nh_state(jg)%prog(n_save)%vn)
          CALL copy(p_nh_state(jg)%prog(n_now)%w,p_nh_state(jg)%prog(n_save)%w)
          CALL copy(p_nh_state(jg)%prog(n_now)%rho,p_nh_state(jg)%prog(n_save)%rho)
          CALL copy(p_nh_state(jg)%prog(n_now)%theta_v,p_nh_state(jg)%prog(n_save)%theta_v)
#ifndef _OPENACC
!$OMP END PARALLEL
#endif
        ELSE IF (lnest_active) THEN ! optimized copy restricted to nest boundary points
          CALL save_progvars(jg,p_nh_state(jg)%prog(n_now),p_nh_state(jg)%prog(n_save))
        ENDIF

      ENDIF


      ! Set local variable for rcf-time levels
      n_now_rcf = nnow_rcf(jg)
      n_new_rcf = nnew_rcf(jg)

#ifdef MESSY
      IF (timers_level > 3) CALL timer_start(timer_extra3)
      CALL messy_global_start(jg)
      CALL messy_local_start(jg)
      CALL messy_vdiff(jg)
      IF (timers_level > 3) CALL timer_stop(timer_extra3)
#endif

      !
      ! Update model date (for local patch!) - Note that for the
      ! top-level patch, this is omitted, since the update has already
      ! happened in the calling subroutine.
      datetime_local(jg)%ptr = datetime_local(jg)%ptr + mtime_dt_loc
      time_diff  => newTimedelta("PT0S")
      time_diff  =  getTimeDeltaFromDateTime(datetime_local(jg)%ptr, time_config%tc_exp_startdate)
      sim_time =  getTotalMillisecondsTimedelta(time_diff, datetime_local(jg)%ptr)*1.e-3_wp
      CALL deallocateTimedelta(time_diff)

      !time_diff  => newTimedelta("PT0S")
      !time_diff  =  getTimeDeltaFromDateTime(datetime_local(jg)%ptr, time_config%tc_exp_startdate)
      !sim_time =  getTotalMillisecondsTimedelta(time_diff, datetime_local(jg)%ptr)*1.e-3_wp
      !CALL deallocateTimedelta(time_diff)
#ifdef COUP_OAS_ICON
    ! compute oasis' send variables
    !
    nlev = p_patch(1)%nlev

    ! SBr, CHa: coupling counter
    sim_time_oas = FLOOR(sim_time * 10.)

    !CALL diagnose_pres_temp (p_metrics, pt_prog, pt_prog_rcf,   &
    !     &                      pt_diag, pt_patch,                 &
    !     &                      lnd_prog          = lnd_prog_new)
    
    CALL diag_for_output_dyn()

    IF(.NOT.atm_phy_nwp_config(1)%is_les_phy) THEN
      ! diagnostics which are only required for output
      CALL nwp_diag_for_output(datetime_local(jg)%ptr, kstart_moist(1), & !in
        &                      ih_clch(1), ih_clcm(1), & !in
        &                      phy_params(1), & !in
        &                      p_patch(1), & !in
        &                      p_nh_state(1)%metrics, & !in
        &                      p_nh_state(1)%prog(nnow(1)), & !in  !nnow or nnew?
        &                      p_nh_state(1)%prog(nnow_rcf(1)), & !in  !nnow or nnew?
        &                      p_nh_state(1)%diag, & !in
        &                      p_lnd_state(1)%diag_lnd, & !in
        & p_lnd_state(1)%prog_lnd(nnow_rcf(1)), & !in
        & p_lnd_state(1)%prog_wtr(nnow_rcf(1)), & !inout
        &                      ext_data(1), & !in
        &                      prm_diag(1) ) !inout
    END IF

    oas_snd_field(:,7) = 0._wp
    oas_snd_field(:,8) = 0._wp
    
    rl_start = grf_bdywidth_c+1
    rl_end   = min_rlcell_int
    i_startblk = p_patch(1)%cells%start_blk(rl_start, 1)
    i_endblk   = p_patch(1)%cells%end_blk(rl_end, MAX(1,p_patch(1)%n_childdom))
    DO jb = i_startblk, i_endblk
      CALL get_indices_c(p_patch(1), jb, i_startblk, i_endblk, i_startidx, &
        i_endidx, rl_start, rl_end)
      DO jc = i_startidx, i_endidx
        ii = idx_1d(jc,jb)
        oas_snd_field(ii,1)  = p_nh_state(1)%diag%temp(jc,nlev,jb)
        oas_snd_field(ii,2)  = p_nh_state(1)%diag%u(jc,nlev,jb)  ! Slavko: check on which level !!
        oas_snd_field(ii,3)  = p_nh_state(1)%diag%v(jc,nlev,jb)
        oas_snd_field(ii,4)  = p_nh_state(1)%prog(nnow(1))%tracer(jc,nlev,jb,iqv)
        oas_snd_field(ii,5)  = &
          p_nh_state(1)%metrics%z_mc(jc,nlev,jb) - p_nh_state(1)%metrics%z_ifc(jc,nlev+1,jb)
        oas_snd_field(ii,6)  = p_nh_state(1)%diag%pres_sfc(jc,jb)
        oas_snd_field(ii,7)  = MAX(0._wp, prm_diag(1)%swflxsfc(jc,jb) - prm_diag(1)%swflx_dn_sfc_diff(jc,jb))
        oas_snd_field(ii,8)  = prm_diag(1)%swflx_dn_sfc_diff(jc,jb)
        oas_snd_field(ii,9)  = -prm_diag(1)%lwflxsfc(jc,jb) - prm_diag(1)%lwflx_up_sfc(jc,jb)
        oas_snd_field(ii,10) = prm_diag(1)%rain_con_rate(jc,jb) + prm_diag(1)%snow_con_rate(jc,jb)
        oas_snd_field(ii,11) = prm_diag(1)%rain_gsp_rate(jc,jb) + prm_diag(1)%snow_gsp_rate(jc,jb)! + &
        !  prm_diag(1)%ice_gsp_rate(jc,jb) + prm_diag(1)%graupel_gsp_rate(jc,jb) + &
        ! prm_diag(1)%hail_gsp_rate(jc,jb)
      END DO
    END DO

    DO ii = 1, SIZE(oas_snd_meta)
      CALL oasis_put(oas_snd_meta(ii)%vid, sim_time_oas, oas_snd_field(:,ii), oas_error)
      IF(msg_level > 30 ) &
        WRITE(oas_message,*) 'Sending  ', oas_snd_meta(ii)%clpname, 'at 0'
      CALL message(routine, oas_message)
      IF (oas_error .NE. OASIS_Ok .AND. oas_error .LT. OASIS_Sent) THEN
        WRITE(oas_message,*) 'Failure in oasis_put of ', oas_snd_meta(ii)%clpname
        CALL oasis_abort(oas_comp_id, oas_comp_name, oas_message)
      END IF
    END DO

    IF(msg_level > 30 ) THEN
      WRITE(oas_message,*) 'ICONOAS: sim_time= ', sim_time
      CALL message(routine, oas_message)
    END IF

    DO oas_i = 1, SIZE(oas_rcv_meta)
      IF(msg_level > 30 ) THEN
        WRITE(oas_message,*) 'Receiving  ', oas_rcv_meta(oas_i)%clpname,&
          " at t=", sim_time_oas
        CALL message(routine, oas_message)
      END IF

      CALL oasis_get(oas_rcv_meta(oas_i)%vid, sim_time_oas, oas_rcv_field(:,oas_i), oas_error)

      IF(msg_level > 30 ) THEN
        WRITE(oas_message,*) 'Received  ', oas_rcv_meta(oas_i)%clpname, 'at t=', sim_time_oas
        CALL message(routine, oas_message)
      END IF

      IF (oas_error .NE. OASIS_Ok .AND. oas_error .LT. OASIS_Recvd) THEN
        WRITE(oas_message,*) 'Failure in oasis_get of ', oas_rcv_meta(oas_i)%clpname
        CALL oasis_abort(oas_comp_id, oas_comp_name, oas_message)
      END IF

      ! transform thru oasis recieved variable into something icon can use
      !
      !DO jb = 1, p_patch(jg)%nblks_c
      !  DO jc = 
      !    glb_idx_1d = p_patch(jg)%cells%decomp_info%glb_index(idx_1d(jc,jb))
      !    oas_rcv_field(jc,jb) = oas_rcv_field(oas_i,glb_idx_1d)
      !  END DO
      !END DO
      
      IF(msg_level > 30 ) THEN
        WRITE(oas_message,*) "ICONOAS: icon transforming received vars"
        CALL message(routine, oas_message)
      END IF

      rl_start = grf_bdywidth_c+1
      rl_end   = min_rlcell_int
      i_startblk = p_patch(1)%cells%start_blk(rl_start, 1)
      i_endblk   = p_patch(1)%cells%end_blk(rl_end, MAX(1,p_patch(1)%n_childdom))
      DO jb = i_startblk, i_endblk
        CALL get_indices_c(p_patch(1), jb, i_startblk, i_endblk, i_startidx,&
          i_endidx, rl_start, rl_end)
        DO jc = i_startidx, i_endidx
          ii = idx_1d(jc,jb)
          oas_rcv_field_icon(jc,jb,oas_i) = oas_rcv_field(ii,oas_i)
        END DO
      END DO

      IF(msg_level > 30 ) THEN
        WRITE(oas_message,*) "ICONOAS: icon transformed received vars"
        CALL message(routine, oas_message)

      ! check:
        WRITE(oas_message,*) "ICONOAS: for ", oas_i, " got rank0-local min, max=", &
          MINVAL(oas_rcv_field(:,oas_i)), MAXVAL(oas_rcv_field(:,oas_i))
        CALL message(routine, oas_message)
      END IF
    END DO
    
#endif

      IF (itime_scheme == 1) THEN
        !------------------
        ! Pure advection
        !------------------

        ! Print control output for maximum horizontal and vertical wind speed
        !
        ! 2 Cases:
        ! msg_level E [12, inf[: print max/min output for every domain and every transport step
        ! msg_level E [ 8,  11]: print max/min output for global domain and every transport step
        IF (msg_level >= 12) THEN
          CALL print_maxwinds(p_patch(jg), p_nh_state(jg)%prog(nnow(jg))%vn,   &
            p_nh_state(jg)%prog(nnow(jg))%w)
        ELSE IF (msg_level >= 8) THEN
          IF (jg == 1) THEN
            CALL print_maxwinds(p_patch(jg), p_nh_state(jg)%prog(nnow(jg))%vn, &
              p_nh_state(jg)%prog(nnow(jg))%w)
          ENDIF
        ENDIF

#ifdef MESSY
        IF (timers_level > 3) CALL timer_start(timer_extra4)
        CALL messy_beforeadv(jg)
        IF (timers_level > 3) CALL timer_stop(timer_extra4)
#endif

        SELECT CASE ( TRIM(nh_test_name) )

        CASE ('PA') ! solid body rotation

          ! set time-variant vertical velocity
          CALL set_nh_w_rho( p_patch(jg),p_nh_state(jg)%metrics,                    &! in
            & jstep_adv(jg)%marchuk_order, dt_loc, sim_time-dt_loc,                 &! in
            &               p_nh_state(jg)%prog(nnew(jg))%w,                        &! inout
            &               p_nh_state(jg)%diag%pres,                               &! inout
            &               p_nh_state(jg)%diag%rho_ic                              )! inout

        CASE ('DF1', 'DF2', 'DF3', 'DF4') ! deformational flow

          ! get velocity field
          CALL get_nh_df_velocity( p_patch(jg), p_nh_state(jg)%prog(nnew(jg)), &
            &                     nh_test_name, rotate_axis_deg,               &
            &                     sim_time-dt_loc+dt_loc )


          ! get mass flux and new \rho. The latter one is only computed,
          ! if the density equation is re-integrated.
          CALL integrate_density_pa(p_patch(jg), p_int_state(jg),  & !in
            &                     p_nh_state(jg)%prog(nnow(jg)),   & !in
            &                     p_nh_state(jg)%prog(nnew(jg)),   & !in
            &                     p_nh_state(jg)%metrics,          & !in
            &                     p_nh_state(jg)%diag, dt_loc,     & !inout,in
            &                     jstep_adv(jg)%marchuk_order,     & !in
            &                     lcoupled_rho                     )


        CASE ('DCMIP_PA_12', 'dcmip_pa_12')

          ! get velocity field for the DCMIP Hadley-like meridional circulation test
          !
          CALL set_nh_velocity_hadley( p_patch(jg), p_nh_state(jg)%prog(nnew(jg)), & !in,inout
            &                          p_nh_state(jg)%diag, p_int_state(jg),       & !in
            &                          p_nh_state(jg)%metrics,                     & !in
            &                          sim_time-dt_loc+dt_loc)                       !in

          ! get mass flux and updated density for the DCMIP Hadley-like
          ! meridional circulation test
          !
          CALL integrate_density_pa(p_patch(jg), p_int_state(jg),  & !in
            &                     p_nh_state(jg)%prog(nnow(jg)),   & !in
            &                     p_nh_state(jg)%prog(nnew(jg)),   & !in
            &                     p_nh_state(jg)%metrics,          & !in
            &                     p_nh_state(jg)%diag, dt_loc,     & !inout,in
            &                     jstep_adv(jg)%marchuk_order,     & !in
            &                     lcoupled_rho                     )
        END SELECT


        ! Diagnose some velocity-related quantities for the tracer
        ! transport scheme
        CALL prepare_tracer( p_patch(jg), p_nh_state(jg)%prog(nnow(jg)),  &! in
          &         p_nh_state(jg)%prog(nnew(jg)),                        &! in
          &         p_nh_state(jg)%metrics, p_int_state(jg),              &! in
          &         ndyn_substeps_var(jg), .TRUE., .TRUE.,                &! in
          &         advection_config(jg)%lfull_comp,                      &! in
          &         p_nh_state(jg)%diag,                                  &! inout
          &         prep_adv(jg)%vn_traj, prep_adv(jg)%mass_flx_me,       &! inout
          &         prep_adv(jg)%w_traj, prep_adv(jg)%mass_flx_ic,        &! inout
          &         prep_adv(jg)%topflx_tra                               )! out

        CALL compute_airmass(p_patch(jg),                   &
          &                  p_nh_state(jg)%metrics,        &
          &                  p_nh_state(jg)%prog(nnow(jg)), &
          &                  p_nh_state(jg)%diag, itlev = 1)


        ! Update air mass in layer.  Air mass is needed by both the transport and physics.
        CALL compute_airmass(p_patch(jg),                   &
          &                  p_nh_state(jg)%metrics,        &
          &                  p_nh_state(jg)%prog(nnew(jg)), &
          &                  p_nh_state(jg)%diag, itlev = 2)

        CALL step_advection( p_patch(jg), p_int_state(jg), dt_loc,       & !in
          &        jstep_adv(jg)%marchuk_order,                          & !in
          &        p_nh_state(jg)%prog(n_now_rcf)%tracer,                & !in
          &        prep_adv(jg)%mass_flx_me, prep_adv(jg)%vn_traj,       & !in
          &        prep_adv(jg)%mass_flx_ic, prep_adv(jg)%w_traj,        & !in
          &        p_nh_state(jg)%metrics%ddqz_z_full,                   & !in
          &        p_nh_state(jg)%diag%airmass_new,                      & !in
          &        p_nh_state(jg)%diag%airmass_now,                      & !in
          &        p_nh_state(jg)%diag%grf_tend_tracer,                  & !inout
          &        p_nh_state(jg)%prog(n_new_rcf)%tracer,                & !inout
          &        p_nh_state(jg)%diag%hfl_tracer,                       & !out
          &        p_nh_state(jg)%diag%vfl_tracer,                       & !out
          &        opt_topflx_tra=prep_adv(jg)%topflx_tra,               & !in
          &        opt_q_int=p_nh_state(jg)%diag%q_int,                  & !out
          &        opt_ddt_tracer_adv=p_nh_state(jg)%diag%ddt_tracer_adv ) !out

#ifdef MESSY
        IF (timers_level > 3) CALL timer_start(timer_extra5)
        CALL messy_afteradv(jg)
        IF (timers_level > 3) CALL timer_stop(timer_extra5)
#endif

      ELSE  ! itime_scheme /= 1


        ! artificial forcing (Held-Suarez test forcing)
        !!!!!!!!
        ! re-check: iadv_rcf -> ndynsubsteps
        !!!!!!!!
        IF ( iforcing == iheldsuarez) THEN
          CALL held_suarez_nh_interface (p_nh_state(jg)%prog(nnow(jg)), p_patch(jg), &
                                         p_int_state(jg),p_nh_state(jg)%metrics,  &
                                         p_nh_state(jg)%diag)
        ENDIF

        ! For real-data runs, perform an extra diffusion call before the first time
        ! step because no other filtering of the interpolated velocity field is done
        !
        ! For the time being, we hand over the dynamics time step and replace iadv_rcf by
        ! ndyn_substeps (for bit-reproducibility).
        IF (ldynamics .AND. .NOT.ltestcase .AND. linit_dyn(jg) .AND. diffusion_config(jg)%lhdiff_vn .AND. &
            init_mode /= MODE_IAU .AND. init_mode /= MODE_IAU_OLD) THEN
          CALL diffusion(p_nh_state(jg)%prog(nnow(jg)), p_nh_state(jg)%diag,       &
            p_nh_state(jg)%metrics, p_patch(jg), p_int_state(jg), dt_loc/ndyn_substeps, .TRUE.)
        ENDIF

        IF (itype_comm == 1) THEN

          IF (ldynamics) THEN

            ! dynamics integration with substepping
            !
            CALL perform_dyn_substepping (p_patch(jg), p_nh_state(jg), p_int_state(jg), &
              &                           prep_adv(jg), jstep, iau_iter, dt_loc, datetime_local(jg)%ptr)

            ! diffusion at physics time steps
            !
            IF (diffusion_config(jg)%lhdiff_vn .AND. lhdiff_rcf) THEN
              CALL diffusion(p_nh_state(jg)%prog(nnew(jg)), p_nh_state(jg)%diag,     &
                &            p_nh_state(jg)%metrics, p_patch(jg), p_int_state(jg),   &
                &            dt_loc/ndyn_substeps, .FALSE.)
            ENDIF

          ELSE IF (iforcing == inwp .OR. iforcing == iecham) THEN
            CALL add_slowphys(p_nh_state(jg), p_patch(jg), nnow(jg), nnew(jg), dt_loc)
          ENDIF
        ELSE
          CALL finish (routine, 'itype_comm /= 1 currently not implemented')
        ENDIF


#ifdef MESSY
        IF (timers_level > 3) CALL timer_start(timer_extra4)
        CALL messy_beforeadv(jg)
        IF (timers_level > 3) CALL timer_stop(timer_extra4)
#endif

        ! 5. tracer advection
        !-----------------------
        IF ( ltransport) THEN

          IF (lart) THEN
            CALL art_emission_interface(                       &
              &      ext_data(jg),                             &!in
              &      p_patch(jg),                              &!in
              &      dt_loc,                                   &!in
              &      p_nh_state(jg),                           &!in
              &      prm_diag(jg),                             &!in
              &      p_lnd_state(jg)%diag_lnd,                 &!in
              &      p_nh_state(jg)%prog(nnew(jg))%rho,        &!in
              &      datetime_local(jg)%ptr,                    &!in
              &      p_nh_state(jg)%prog(n_now_rcf)%tracer)     !inout
          ENDIF


          IF (msg_level >= 13) THEN
            WRITE(message_text,'(a,i2)') 'call advection  DOM:',jg
            CALL message('integrate_nh', TRIM(message_text))
          ENDIF

          CALL step_advection( p_patch(jg), p_int_state(jg), dt_loc,         & !in
            &          jstep_adv(jg)%marchuk_order,                          & !in
            &          p_nh_state(jg)%prog(n_now_rcf)%tracer,                & !in
            &          prep_adv(jg)%mass_flx_me, prep_adv(jg)%vn_traj,       & !in
            &          prep_adv(jg)%mass_flx_ic, prep_adv(jg)%w_traj,        & !in
            &          p_nh_state(jg)%metrics%ddqz_z_full,                   & !in
            &          p_nh_state(jg)%diag%airmass_new,                      & !in
            &          p_nh_state(jg)%diag%airmass_now,                      & !in
            &          p_nh_state(jg)%diag%grf_tend_tracer,                  & !inout
            &          p_nh_state(jg)%prog(n_new_rcf)%tracer,                & !inout
            &          p_nh_state(jg)%diag%hfl_tracer,                       & !out
            &          p_nh_state(jg)%diag%vfl_tracer,                       & !out
            &          opt_topflx_tra=prep_adv(jg)%topflx_tra,               & !in
            &          opt_q_int=p_nh_state(jg)%diag%q_int,                  & !out
            &          opt_ddt_tracer_adv=p_nh_state(jg)%diag%ddt_tracer_adv ) !out

          IF (iprog_aero == 1) THEN

            CALL aerosol_2D_advection( p_patch(jg), p_int_state(jg), dt_loc, & !in
            &          prm_diag(jg)%aerosol, prep_adv(jg)%vn_traj,           & !inout, in
            &          prep_adv(jg)%mass_flx_me, prep_adv(jg)%mass_flx_ic,   & !in
            &          p_nh_state(jg)%metrics%ddqz_z_full_e,                 & !in
            &          p_nh_state(jg)%diag%airmass_now,                      & !in
            &          p_nh_state(jg)%diag%airmass_new                       ) !in

          ENDIF

        ! ART tracer sedimentation:
        !     Internal substepping with ndyn_substeps_var(jg)
        !-----------------------
          IF (lart) THEN
            CALL art_sedi_interface( p_patch(jg),             &!in
               &      dt_loc,                                 &!in
               &      p_nh_state(jg)%prog(n_new_rcf),         &!in
               &      p_nh_state(jg)%metrics,                 &!in
               &      p_nh_state(jg)%prog(nnew(jg))%rho,      &!in
               &      p_nh_state(jg)%diag,                    &!in
               &      prm_diag(jg),                           &!in
               &      ndyn_substeps_var(jg),                  &!in
               &      p_nh_state(jg)%prog(n_new_rcf)%tracer,  &!inout
               &      .TRUE.)                                  !print CFL number
          ENDIF ! lart

        ENDIF !ltransport

#ifdef MESSY
        IF (timers_level > 3) CALL timer_start(timer_extra5)
        CALL messy_afteradv(jg)
        IF (timers_level > 3) CALL timer_stop(timer_extra5)
#endif



        ! Apply boundary nudging in case of one-way nesting
        IF (jg > 1 ) THEN

          IF (lfeedback(jg) .AND. l_density_nudging .AND. grf_intmethod_e <= 4) THEN
            IF (ltimer)            CALL timer_start(timer_nesting)
            IF (timers_level >= 2) CALL timer_start(timer_nudging)
            CALL density_boundary_nudging(jg,nnew(jg),REAL(ndyn_substeps,wp))
            IF (timers_level >= 2) CALL timer_stop(timer_nudging)
            IF (ltimer)            CALL timer_stop(timer_nesting)
          ELSE IF (.NOT. lfeedback(jg)) THEN
            IF (ltimer)            CALL timer_start(timer_nesting)
            IF (timers_level >= 2) CALL timer_start(timer_nudging)
            CALL nest_boundary_nudging(jg,nnew(jg),nnew_rcf(jg),REAL(ndyn_substeps,wp))
            IF (timers_level >= 2) CALL timer_stop(timer_nudging)
            IF (ltimer)            CALL timer_stop(timer_nesting)
          ENDIF

        ENDIF

        IF ( ( iforcing==inwp .OR. iforcing==iecham ) ) THEN

          ! Determine which physics packages must be called/not called at the current
          ! time step
          IF ( iforcing==inwp ) THEN
            CALL time_ctrl_physics ( dt_phy, dt_loc, jg,        &! in
              &                      .FALSE.,                   &! in
              &                      t_elapsed_phy,             &! inout
              &                      lcall_phy )                 ! out

            IF (msg_level >= 13) THEN
              WRITE(message_text,'(a,i2,a,5l2,a,6l2)') 'call phys. proc DOM:', &
                &  jg ,'   SP:', lcall_phy(jg,1:5), '   FP:',lcall_phy(jg,6:10)
              CALL message(TRIM(routine), TRIM(message_text))
            END IF
          END IF

          IF (atm_phy_nwp_config(jg)%is_les_phy) THEN


            ! les physics
            CALL les_phy_interface(lcall_phy(jg,:), .FALSE.,         & !in
              &                  lredgrid_phys(jg),                  & !in
              &                  dt_loc,                             & !in
              &                  t_elapsed_phy(jg,:),                & !in
              &                  nstep_global,                       & !in
              &                  datetime_local(jg)%ptr,              & !in
              &                  p_patch(jg)  ,                      & !in
              &                  p_int_state(jg),                    & !in
              &                  p_nh_state(jg)%metrics ,            & !in
              &                  p_patch(jgp),                       & !in
              &                  ext_data(jg)           ,            & !in
              &                  p_nh_state(jg)%prog(nnew(jg)) ,     & !inout
              &                  p_nh_state(jg)%prog(n_now_rcf),     & !in for tke
              &                  p_nh_state(jg)%prog(n_new_rcf) ,    & !inout
              &                  p_nh_state(jg)%diag ,               & !inout
              &                  prm_diag  (jg),                     & !inout
              &                  prm_nwp_tend(jg),                   &
              &                  p_lnd_state(jg)%diag_lnd,           &
              &                  p_lnd_state(jg)%prog_lnd(n_now_rcf),& !inout
              &                  p_lnd_state(jg)%prog_lnd(n_new_rcf),& !inout
              &                  p_lnd_state(jg)%prog_wtr(n_now_rcf),& !inout
              &                  p_lnd_state(jg)%prog_wtr(n_new_rcf),& !inout
              &                  p_nh_state_lists(jg)%prog_list(n_new_rcf) ) !in

          ELSE ! is_les_phy

            SELECT CASE (iforcing)

            CASE (inwp) ! iforcing

              ! nwp physics
              CALL nwp_nh_interface(lcall_phy(jg,:), .FALSE.,          & !in
                &                  lredgrid_phys(jg),                  & !in
                &                  dt_loc,                             & !in
                &                  t_elapsed_phy(jg,:),                & !in
                &                  datetime_local(jg)%ptr,              & !in
                &                  p_patch(jg)  ,                      & !in
                &                  p_int_state(jg),                    & !in
                &                  p_nh_state(jg)%metrics ,            & !in
                &                  p_patch(jgp),                       & !in
                &                  ext_data(jg)           ,            & !in
                &                  p_nh_state(jg)%prog(nnew(jg)) ,     & !inout
                &                  p_nh_state(jg)%prog(n_now_rcf),     & !in for tke
                &                  p_nh_state(jg)%prog(n_new_rcf) ,    & !inout
                &                  p_nh_state(jg)%diag ,               & !inout
                &                  prm_diag  (jg),                     & !inout
                &                  prm_nwp_tend(jg),                   &
                &                  p_lnd_state(jg)%diag_lnd,           &
                &                  p_lnd_state(jg)%prog_lnd(n_now_rcf),& !inout
                &                  p_lnd_state(jg)%prog_lnd(n_new_rcf),& !inout
                &                  p_lnd_state(jg)%prog_wtr(n_now_rcf),& !inout
                &                  p_lnd_state(jg)%prog_wtr(n_new_rcf),& !inout
                &                  p_nh_state_lists(jg)%prog_list(n_new_rcf) ) !in

            CASE (iecham) ! iforcing

              ! echam physics
              IF (ltimer) CALL timer_start(timer_iconam_echam)
              CALL interface_iconam_echam( dt_loc                         ,& !in
                &                          datetime_local(jg)%ptr          ,& !in
                &                          p_patch(jg)                    ,& !in
                &                          p_int_state(jg)                ,& !in
                &                          p_nh_state(jg)%metrics         ,& !in
                &                          p_nh_state(jg)%prog(nnew(jg))  ,& !inout
                &                          p_nh_state(jg)%prog(n_new_rcf) ,& !inout
                &                          p_nh_state(jg)%diag,            &            
                &                          p_nh_state_lists(jg)%prog_list(n_new_rcf))
            
              IF (ltimer) CALL timer_stop(timer_iconam_echam)

            END SELECT ! iforcing

          END IF ! is_les_phy

          ! Boundary interpolation of land state variables entering into radiation computation
          ! if a reduced grid is used in the child domain(s)
          IF (ltimer)            CALL timer_start(timer_nesting)
          IF (timers_level >= 2) CALL timer_start(timer_bdy_interp)
          DO jn = 1, p_patch(jg)%n_childdom

            jgc = p_patch(jg)%child_id(jn)
            IF (.NOT. p_patch(jgc)%ldom_active) CYCLE

            IF (patch_weight(jgc) > 0._wp) THEN
              CALL p_bcast(lcall_phy(jgc,itrad),      p_patch(jgc)%proc0, p_comm_work)
              CALL p_bcast(t_elapsed_phy(jgc,itrad),  p_patch(jgc)%proc0, p_comm_work)
            ENDIF

            ! Determine if radiation will be called in the nested domain during the subsequent two (small) time steps
            IF (lredgrid_phys(jgc) .AND. atm_phy_nwp_config(jgc)%lproc_on(itrad) .AND. .NOT. lcall_phy(jgc,itrad) &
              .AND. t_elapsed_phy(jgc,itrad) + dt_loc >= 0.99999999_wp*dt_phy(jgc,itrad) ) THEN
              lcall_rrg = .TRUE.
            ELSE
              lcall_rrg = .FALSE.
            ENDIF

            IF (lcall_rrg .AND. atm_phy_nwp_config(jgc)%inwp_surface >= 1) THEN
              CALL interpol_rrg_grf(jg, jgc, jn, nnew_rcf(jg))
            ENDIF
            IF (lcall_rrg .AND. atm_phy_nwp_config(jgc)%latm_above_top) THEN
              CALL copy_rrg_ubc(jg, jgc)
            ENDIF
          ENDDO
          IF (timers_level >= 2) CALL timer_stop(timer_bdy_interp)
          IF (ltimer)            CALL timer_stop(timer_nesting)

        ENDIF !iforcing

        ! Terminator toy chemistry
        !
        ! So far it can only be activated for testcases and not for real-cases, 
        ! since the initialization is done in init_nh_testcase. However, 
        ! nothing speaks against combining toy chemistry with real case runs.
        IF (ltestcase .AND. is_toy_chem) THEN
          CALL dcmip_terminator_interface (p_patch(jg),            & !in
            &                              p_nh_state(jg)%metrics, & !in
            &                              p_nh_state(jg)%prog,    & !inout
            &                              p_nh_state(jg)%diag,    & !inout
            &                              datetime_local(jg)%ptr, & !in
            &                              dt_loc                  ) !in
        ENDIF

#ifdef MESSY
        IF (timers_level > 3) CALL timer_start(timer_extra6)
        CALL messy_physc(jg)
        IF (timers_level > 3) CALL timer_stop(timer_extra6)
#endif

      ENDIF  ! itime_scheme

      ! Update nudging tendency fields for limited-area mode
      IF (jg == 1 .AND. l_limited_area) THEN

         IF (latbc_config%itype_latbc > 0) THEN ! use time-dependent boundary data

            IF (num_prefetch_proc >= 1) THEN

              ! Asynchronous LatBC read-in:
              ! update the coefficients for the linear interpolation
              CALL update_lin_interpolation(latbc, datetime_local(jg)%ptr)
              CALL prep_outer_bdy_nudging(p_patch(jg),p_nh_state(jg)%prog(nnew(jg)),  &
                &  p_nh_state(jg)%prog(n_new_rcf),                                    &
                &  p_nh_state(jg)%metrics,p_nh_state(jg)%diag,                        &
                &  p_latbc_old=latbc%latbc_data(latbc%prev_latbc_tlev())%atm,         &
                &  p_latbc_new=latbc%latbc_data(latbc%new_latbc_tlev)%atm)
            ELSE
              
              ! update the coefficients for the linear interpolation
              CALL update_lin_interc(datetime_local(jg)%ptr)
              CALL prep_outer_bdy_nudging(p_patch(jg),p_nh_state(jg)%prog(nnew(jg)),     &
                &                         p_nh_state(jg)%prog(n_new_rcf),                &  
                &                         p_nh_state(jg)%metrics,p_nh_state(jg)%diag,    &
                &                         p_latbc_old=p_latbc_data(last_latbc_tlev)%atm, &
                &                         p_latbc_new=p_latbc_data(read_latbc_tlev)%atm)

            ENDIF

         ELSE ! constant lateral boundary data

            CALL prep_outer_bdy_nudging(p_patch(jg),p_nh_state(jg)%prog(nnew(jg)),p_nh_state(jg)%prog(n_new_rcf), &
                 p_nh_state(jg)%metrics,p_nh_state(jg)%diag,p_latbc_const=p_nh_state(jg)%prog(nsav2(jg)))

         ENDIF

         ! Apply nudging at the lateral boundaries
         CALL outer_boundary_nudging (jg, nnew(jg), n_new_rcf, REAL(ndyn_substeps,wp))

      ENDIF



      ! Check if at least one of the nested domains is active
      !
      IF (p_patch(jg)%n_childdom > 0) THEN
        lnest_active = .FALSE.
        DO jn = 1, p_patch(jg)%n_childdom
          jgc = p_patch(jg)%child_id(jn)
          IF (p_patch(jgc)%ldom_active) lnest_active = .TRUE.
        ENDDO
      ENDIF

      ! If there are nested domains...
      IF (p_patch(jg)%n_childdom > 0 .AND. lnest_active ) THEN


        IF (ndyn_substeps_var(jg) == 1) THEN
          n_now_grf  = nnow(jg)
        ELSE
          n_now_grf  = nsav1(jg)
        ENDIF

        rdt_loc     = 1._wp/dt_loc
        dt_sub      = dt_loc/2._wp    ! (adv.) time step on next refinement level
        mtime_dt_sub => newTimedelta(mtime_dt_loc)
        mtime_dt_sub = mtime_dt_sub*0.5_wp
        rdtmflx_loc = 1._wp/(dt_loc*(REAL(MAX(1,ndyn_substeps_var(jg)-1),wp)/REAL(ndyn_substeps_var(jg),wp)))

        IF (ltimer)            CALL timer_start(timer_nesting)
        IF (timers_level >= 2) CALL timer_start(timer_bdy_interp)

        ! Compute time tendencies for interpolation to refined mesh boundaries
        CALL compute_tendencies (jg,nnew(jg),n_now_grf,n_new_rcf,n_now_rcf, &
          &                      rdt_loc,rdtmflx_loc)

        ! Loop over nested domains
        DO jn = 1, p_patch(jg)%n_childdom

          jgc = p_patch(jg)%child_id(jn)

          ! Interpolate tendencies to lateral boundaries of refined mesh (jgc)
          IF (p_patch(jgc)%ldom_active) THEN
            CALL boundary_interpolation(jg, jgc,                   &
              &  n_now_grf,nnow(jgc),n_now_rcf,nnow_rcf(jgc),      &
              &  prep_adv(jg)%mass_flx_me,prep_adv(jgc)%mass_flx_me)
          ENDIF

        ENDDO
        IF (timers_level >= 2) CALL timer_stop(timer_bdy_interp)

        ! prep_bdy_nudging can not be called using delayed requests!
        DO jn = 1, p_patch(jg)%n_childdom

          jgc = p_patch(jg)%child_id(jn)
          IF (.NOT. p_patch(jgc)%ldom_active) CYCLE
          ! If feedback is turned off for child domain, compute parent-child
          ! differences for boundary nudging
          ! *** prep_bdy_nudging adapted for reduced calling frequency of tracers ***
          IF (lfeedback(jgc) .AND. l_density_nudging .AND. grf_intmethod_e <= 4) THEN
            IF (timers_level >= 2) CALL timer_start(timer_nudging)
            CALL prep_rho_bdy_nudging(jg,jgc)
            IF (timers_level >= 2) CALL timer_stop(timer_nudging)
          ELSE IF (.NOT. lfeedback(jgc)) THEN
            IF (timers_level >= 2) CALL timer_start(timer_nudging)
            CALL prep_bdy_nudging(jg,jgc)
            IF (timers_level >= 2) CALL timer_stop(timer_nudging)
          ENDIF
        ENDDO
        IF (ltimer)            CALL timer_stop(timer_nesting)

        DO jn = 1, p_patch(jg)%n_childdom

          jgc = p_patch(jg)%child_id(jn)
          IF (.NOT. p_patch(jgc)%ldom_active) CYCLE

          IF(p_patch(jgc)%n_patch_cells > 0) THEN
            IF(proc_split) CALL push_glob_comm(p_patch(jgc)%comm, p_patch(jgc)%proc0)
            ! Recursive call to process_grid_level for child grid level
            CALL integrate_nh( datetime_local, jgc, nstep_global, iau_iter, &
              &                dt_sub, mtime_dt_sub, nsteps_nest, latbc )
            IF(proc_split) CALL pop_glob_comm()
          ENDIF

        ENDDO

        ! clean up
        CALL deallocateTimedelta(mtime_dt_sub)

        IF (ltimer)            CALL timer_start(timer_nesting)
        DO jn = 1, p_patch(jg)%n_childdom

          ! Call feedback to copy averaged prognostic variables from refined mesh back
          ! to the coarse mesh (i.e. from jgc to jg)
          jgc = p_patch(jg)%child_id(jn)
          IF (.NOT. p_patch(jgc)%ldom_active) CYCLE

          IF (lfeedback(jgc)) THEN
            IF (timers_level >= 2) CALL timer_start(timer_feedback)
            IF (ifeedback_type == 1) THEN
              CALL feedback(p_patch, p_nh_state, p_int_state, p_grf_state, p_lnd_state, &
                &           jgc, jg)
            ELSE
              CALL relax_feedback(  p_patch(n_dom_start:n_dom),            &
                & p_nh_state(1:n_dom), p_int_state(n_dom_start:n_dom),     &
                & p_grf_state(n_dom_start:n_dom), prm_diag, jgc, jg, dt_loc)
            ENDIF
            ! Note: the last argument of "feedback" ensures that tracer feedback is
            ! only done for those time steps in which transport and microphysics are called
            IF (timers_level >= 2) CALL timer_stop(timer_feedback)
          ENDIF
        ENDDO
        IF (ltimer)            CALL timer_stop(timer_nesting)

      ENDIF


      ! Average atmospheric variables needed as first guess for data assimilation
      !
      IF ( jg == 1 .AND. is_avgFG_time(datetime_local(jg)%ptr))  THEN
        CALL average_first_guess(p_patch(jg), p_int_state(jg), p_nh_state(jg)%diag, &
          p_nh_state(jg)%prog(nnew(jg)), p_nh_state(jg)%prog(nnew_rcf(jg)))
      ENDIF


      IF (test_mode <= 0) THEN ! ... normal execution of time stepping
        ! Finally, switch between time levels now and new for next time step
        CALL swap(nnow(jg), nnew(jg))

        ! Special treatment for processes (i.e. advection) which can be treated with
        ! reduced calling frequency. Switch between time levels now and new immediately
        ! AFTER the last transport timestep.
        CALL swap(nnow_rcf(jg), nnew_rcf(jg))

      ENDIF


      ! Check if nested domains have to be activated
      IF ( p_patch(jg)%n_childdom > 0 ) THEN

        ! Loop over nested domains
        DO jn = 1, p_patch(jg)%n_childdom
          jgc = p_patch(jg)%child_id(jn)

          IF ( .NOT. p_patch(jgc)%ldom_active .AND. &
            &  (sim_time >= start_time(jgc))  .AND. &
            &  (sim_time <  end_time(jgc))) THEN
            p_patch(jgc)%ldom_active = .TRUE.

            jstep_adv(jgc)%marchuk_order = 0
            datetime_local(jgc)%ptr      = datetime_local(jg)%ptr
            t_elapsed_phy(jgc,:)         = 0._wp
            linit_dyn(jgc)               = .TRUE.
            dt_sub                       = dt_loc/2._wp

            IF (  atm_phy_nwp_config(jgc)%inwp_surface == 1 ) THEN
              CALL aggregate_landvars(p_patch(jg), ext_data(jg),                &
                p_lnd_state(jg)%prog_lnd(nnow_rcf(jg)), p_lnd_state(jg)%diag_lnd)
            ENDIF

            CALL initialize_nest(jg, jgc)

            ! Apply hydrostatic adjustment, using downward integration
            CALL hydro_adjust_downward(p_patch(jgc), p_nh_state(jgc)%metrics,                     &
              p_nh_state(jgc)%prog(nnow(jgc))%rho, p_nh_state(jgc)%prog(nnow(jgc))%exner,         &
              p_nh_state(jgc)%prog(nnow(jgc))%theta_v )

            CALL init_exner_pr(jgc, nnow(jgc))

            ! Activate cold-start mode in TERRA-init routine irrespective of what has been used for the global domain
            init_mode_soil = 1

            IF (iforcing == inwp) THEN
              CALL init_nwp_phy(                           &
                & p_patch(jgc)                            ,&
                & p_nh_state(jgc)%metrics                 ,&
                & p_nh_state(jgc)%prog(nnow(jgc))         ,&
                & p_nh_state(jgc)%diag                    ,&
                & prm_diag(jgc)                           ,&
                & prm_nwp_tend(jgc)                       ,&
                & p_lnd_state(jgc)%prog_lnd(nnow_rcf(jgc)),&
                & p_lnd_state(jgc)%prog_lnd(nnew_rcf(jgc)),&
                & p_lnd_state(jgc)%prog_wtr(nnow_rcf(jgc)),&
                & p_lnd_state(jgc)%prog_wtr(nnew_rcf(jgc)),&
                & p_lnd_state(jgc)%diag_lnd               ,&
                & ext_data(jgc)                           ,&
                & phy_params(jgc), datetime_local(jgc)%ptr,&
                & lnest_start=.TRUE. )

              CALL init_cloud_aero_cpl (datetime_local(jgc)%ptr, p_patch(jgc), p_nh_state(jgc)%metrics, &
                &                       ext_data(jgc), prm_diag(jgc))

              IF (iprog_aero == 1) CALL setup_aerosol_advection(p_patch(jgc))
            ENDIF

            CALL compute_airmass(p_patch(jgc),                   &
              &                  p_nh_state(jgc)%metrics,        &
              &                  p_nh_state(jgc)%prog(nnow(jgc)),&
              &                  p_nh_state(jgc)%diag, itlev = 2 )

            IF ( lredgrid_phys(jgc) ) THEN
              CALL interpol_rrg_grf(jg, jgc, jn, nnow_rcf(jg))
              IF (atm_phy_nwp_config(jgc)%latm_above_top) THEN
                CALL copy_rrg_ubc(jg, jgc)
              ENDIF
            ENDIF

            CALL init_slowphysics (datetime_local(jgc)%ptr, jgc, dt_sub)

            WRITE(message_text,'(a,i2,a,f12.2)') 'domain ',jgc,' started at time ',sim_time
            CALL message('integrate_nh', TRIM(message_text))

          ENDIF
        ENDDO
      ENDIF

#ifdef MESSY
      IF (timers_level > 3) CALL timer_start(timer_extra7)
      CALL messy_local_end(jg)
      CALL messy_global_end(jg)
      IF (timers_level > 3) CALL timer_stop(timer_extra7)
#endif

    ENDDO

    IF (jg == 1 .AND. ltimer) CALL timer_stop(timer_integrate_nh)

  END SUBROUTINE integrate_nh


  !>
  !! Performs dynamical core substepping with respect to physics/transport.
  !!
  !! Perform dynamical core substepping with respect to physics/transport.
  !! Number of substeps is given by ndyn_substeps.
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2014-10-28)
  !!
  SUBROUTINE perform_dyn_substepping (p_patch, p_nh_state, p_int_state, prep_adv, &
    &                                 jstep, iau_iter, dt_phy, mtime_current)

    TYPE(t_patch)       ,INTENT(IN)    :: p_patch

    TYPE(t_nh_state)    ,INTENT(INOUT) :: p_nh_state

    TYPE(t_int_state)   ,INTENT(IN)    :: p_int_state

    TYPE(t_prepare_adv) ,INTENT(INOUT) :: prep_adv

    INTEGER             ,INTENT(IN)    :: jstep     ! number of current (large) time step
                                                    ! performed in current domain
    INTEGER             ,INTENT(IN)    :: iau_iter  ! counter for IAU iteration
    REAL(wp)            ,INTENT(IN)    :: dt_phy    ! physics time step for current patch

    TYPE(datetime)      ,INTENT(IN)    :: mtime_current

    ! local variables
    INTEGER                  :: jg                ! domain ID
    INTEGER                  :: nstep             ! timestep counter
    INTEGER                  :: ndyn_substeps_tot ! total number of dynamics substeps
                                                  ! since last boundary update
    REAL(wp)                 :: dt_dyn            ! dynamics time step
    REAL(wp)                 :: cur_time          ! current time (for IAU)

    LOGICAL                  :: lclean_mflx       ! .TRUE.: first substep
    LOGICAL                  :: l_recompute       ! .TRUE.: recompute velocity tendencies for predictor
                                                  ! (first substep)
    LOGICAL                  :: lsave_mflx
    LOGICAL                  :: lprep_adv         !.TRUE.: do computations for preparing tracer advection in solve_nh
    LOGICAL                  :: l_bdy_nudge
    LOGICAL                  :: llast             !.TRUE.: this is the last substep
    TYPE(timeDelta), POINTER :: time_diff
    !-------------------------------------------------------------------------

    ! get domain ID
    jg = p_patch%id

    ! compute dynamics timestep
    dt_dyn = dt_phy/ndyn_substeps_var(jg)


    IF (jg > 1 .AND. .NOT. lfeedback(jg) .OR. jg == 1 .AND. l_limited_area) THEN
      ! apply boundary nudging if feedback is turned off and in limited-area mode
      l_bdy_nudge = .TRUE.
    ELSE
      l_bdy_nudge = .FALSE.
    ENDIF

    IF ( idiv_method == 1 .AND. (ltransport .OR. p_patch%n_childdom > 0 .AND. grf_intmethod_e >= 5)) THEN
      lprep_adv = .TRUE. ! do computations for preparing tracer advection in solve_nh
    ELSE
      lprep_adv = .FALSE.
    ENDIF

    ! compute airmass \rho*\Delta z [kg m-2] for nnow
    CALL compute_airmass(p_patch,                   &
      &                  p_nh_state%metrics,        &
      &                  p_nh_state%prog(nnow(jg)), &
      &                  p_nh_state%diag, itlev = 1)



    ! perform dynamics substepping
    !
    SUBSTEPS: DO nstep = 1, ndyn_substeps_var(jg)
      
      ! Print control output for maximum horizontal and vertical wind speed
      !
      ! 3 Cases:
      ! msg_level E [12, inf[: print max/min output for every domain and every substep
      ! msg_level E [ 8,  11]: print max/min output for global domain and every substep
      ! msg_level E [ 5,   7]: print max/min output for global domain and first substep
      !
      IF (msg_level >= 12) THEN
        CALL print_maxwinds(p_patch, p_nh_state%prog(nnow(jg))%vn,   &
          p_nh_state%prog(nnow(jg))%w)
      ELSE IF (msg_level >= 8) THEN
        IF (jg == 1) THEN
          CALL print_maxwinds(p_patch, p_nh_state%prog(nnow(jg))%vn, &
            p_nh_state%prog(nnow(jg))%w)
        ENDIF
      ELSE IF (msg_level >= 5) THEN
        IF ( (jg == 1) .AND. (nstep == 1) ) THEN
          CALL print_maxwinds(p_patch, p_nh_state%prog(nnow(jg))%vn, &
            p_nh_state%prog(nnow(jg))%w)
        ENDIF
      ENDIF

      ! total number of dynamics substeps since last boundary update
      ! applicable to refined domains only
      ndyn_substeps_tot = (jstep-1)*ndyn_substeps_var(jg) + nstep

      ! nullify prep_adv fields at first substep
      lclean_mflx = MERGE(.TRUE.,.FALSE.,nstep==1)
      l_recompute = lclean_mflx

      ! logical checking for the last substep
      llast = MERGE(.TRUE.,.FALSE.,nstep==ndyn_substeps_var(jg))

      ! save massflux at first substep
      IF (p_patch%n_childdom > 0 .AND. nstep == 1 ) THEN
        lsave_mflx = .TRUE.
      ELSE
        lsave_mflx = .FALSE.
      ENDIF

      IF ( ANY((/MODE_IAU,MODE_IAU_OLD/)==init_mode) ) THEN ! incremental analysis mode
        time_diff  => newTimedelta("PT0S")
        time_diff  =  getTimeDeltaFromDateTime(mtime_current, time_config%tc_exp_startdate)
        cur_time = REAL(getTotalSecondsTimedelta(time_diff, mtime_current)                  &
             &         -getTotalSecondsTimedelta(timeshift%mtime_shift, mtime_current),wp)  &
             &    +(REAL(nstep-ndyn_substeps_var(jg),wp)-0.5_wp)*dt_dyn
        CALL deallocateTimedelta(time_diff)
        IF (iau_iter == 1) THEN
          CALL compute_iau_wgt(cur_time, dt_dyn, 0.5_wp*dt_iau, lclean_mflx)
        ELSE
          CALL compute_iau_wgt(cur_time, dt_dyn, dt_iau, lclean_mflx)
        ENDIF
      ENDIF

      ! integrate dynamical core
      CALL solve_nh(p_nh_state, p_patch, p_int_state, prep_adv,     &
        &           nnow(jg), nnew(jg), linit_dyn(jg), l_recompute, &
        &           lsave_mflx, lprep_adv, lclean_mflx,             &
        &           nstep, ndyn_substeps_tot-1, l_bdy_nudge, dt_dyn)

      ! now reset linit_dyn to .FALSE.
      linit_dyn(jg) = .FALSE.

      ! compute diffusion at every dynamics substep (.NOT. lhdiff_rcf)
      IF (diffusion_config(jg)%lhdiff_vn .AND. .NOT. lhdiff_rcf)   &
        CALL diffusion(p_nh_state%prog(nnew(jg)), p_nh_state%diag, &
          &            p_nh_state%metrics, p_patch, p_int_state,   &
          &            dt_dyn, .FALSE.)

      IF (llast .OR. advection_config(jg)%lfull_comp) &
        CALL prepare_tracer( p_patch, p_nh_state%prog(nnow(jg)),        &! in
          &                  p_nh_state%prog(nnew(jg)),                 &! in
          &                  p_nh_state%metrics, p_int_state,           &! in
          &                  ndyn_substeps_var(jg), llast, lclean_mflx, &! in
          &                  advection_config(jg)%lfull_comp,           &! in
          &                  p_nh_state%diag,                           &! inout
          &                  prep_adv%vn_traj, prep_adv%mass_flx_me,    &! inout
          &                  prep_adv%w_traj,  prep_adv%mass_flx_ic,    &! inout
          &                  prep_adv%topflx_tra                        )! out

      ! Finally, switch between time levels now and new for next iteration
      !
      ! Note, that we do not swap during the very last iteration.
      ! This final swap is postponed till the end of the integration step.
      IF ( .NOT. llast ) THEN
        CALL swap(nnow(jg), nnew(jg))
      ENDIF

    END DO SUBSTEPS


    ! compute airmass \rho*\Delta z [kg m-2] for nnew
    CALL compute_airmass(p_patch,                   &
      &                  p_nh_state%metrics,        &
      &                  p_nh_state%prog(nnew(jg)), &
      &                  p_nh_state%diag, itlev = 2)


  END SUBROUTINE perform_dyn_substepping


  !-------------------------------------------------------------------------
  !>
  !! Driver routine for initial call of physics routines.
  !! Apart from the full set of slow physics parameterizations, also turbulent transfer is
  !! called, in order to have proper transfer coefficients available at the initial time step.
  !!
  !! This had to be moved ahead of the initial output for the physics fields to be more complete
  !!
  !! @par Revision History
  !! Developed by Guenther Zaengl, DWD (2013-01-04)
  !!
  RECURSIVE SUBROUTINE init_slowphysics (mtime_current, jg, dt_loc)

    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: routine = modname//':init_slowphysics'

    TYPE(datetime), POINTER :: mtime_current
    INTEGER , INTENT(IN)    :: jg           !< current grid level
    REAL(wp), INTENT(IN)    :: dt_loc       !< time step applicable to local grid level

    ! Local variables
    INTEGER                             :: n_now_rcf, nstep
    INTEGER                             :: jgp, jgc, jn
    REAL(wp)                            :: dt_sub       !< (advective) timestep for next finer grid level
    INTEGER                             :: oas_i, rl_start, rl_end, jb, i_startidx, i_endidx, jc, ii, &
                                           i_startblk, i_endblk
    CHARACTER(len=2556)                 :: oas_message

    ! Determine parent domain ID
    IF ( jg > 1) THEN
      jgp = p_patch(jg)%parent_id
    ELSE IF (n_dom_start == 0) THEN
      jgp = 0
    ELSE
      jgp = 1
    ENDIF

    ! Set local variable for rcf-time levels
    n_now_rcf = nnow_rcf(jg)

    IF (iforcing == inwp) THEN
      CALL time_ctrl_physics ( dt_phy, dt_loc, jg,     &! in
        &                      .TRUE.,                 &! in
        &                      t_elapsed_phy,          &! inout
        &                      lcall_phy )              ! out
    END IF


    IF (msg_level >= 7) THEN
      WRITE(message_text,'(a,i2)') 'initial call of (slow) physics, domain ', jg
      CALL message(TRIM(routine), TRIM(message_text))
    ENDIF

    IF (atm_phy_nwp_config(jg)%is_les_phy) THEN

      nstep = 0
      CALL les_phy_interface(lcall_phy(jg,:), .TRUE.,          & !in
        &                  lredgrid_phys(jg),                  & !in
        &                  dt_loc,                             & !in
        &                  dt_phy(jg,:),                       & !in
        &                  nstep,                              & !in
        &                  mtime_current,                      & !in
        &                  p_patch(jg)  ,                      & !in
        &                  p_int_state(jg),                    & !in
        &                  p_nh_state(jg)%metrics ,            & !in
        &                  p_patch(jgp),                       & !in
        &                  ext_data(jg)           ,            & !in
        &                  p_nh_state(jg)%prog(nnow(jg)) ,     & !inout
        &                  p_nh_state(jg)%prog(n_now_rcf) ,    & !inout
        &                  p_nh_state(jg)%prog(n_now_rcf) ,    & !inout
        &                  p_nh_state(jg)%diag,                & !inout
        &                  prm_diag  (jg),                     & !inout
        &                  prm_nwp_tend(jg)                ,   &
        &                  p_lnd_state(jg)%diag_lnd,           &
        &                  p_lnd_state(jg)%prog_lnd(n_now_rcf),& !inout
        &                  p_lnd_state(jg)%prog_lnd(n_now_rcf),& !inout
        &                  p_lnd_state(jg)%prog_wtr(n_now_rcf),& !inout
        &                  p_lnd_state(jg)%prog_wtr(n_now_rcf),& !inout
        &                  p_nh_state_lists(jg)%prog_list(n_now_rcf) ) !in

    ELSE ! is_les_phy

      SELECT CASE (iforcing)

      CASE (inwp) ! iforcing

        ! nwp physics, slow physics forcing
        CALL nwp_nh_interface(lcall_phy(jg,:), .TRUE.,           & !in
          &                  lredgrid_phys(jg),                  & !in
          &                  dt_loc,                             & !in
          &                  dt_phy(jg,:),                       & !in
          &                  mtime_current,                      & !in
          &                  p_patch(jg)  ,                      & !in
          &                  p_int_state(jg),                    & !in
          &                  p_nh_state(jg)%metrics ,            & !in
          &                  p_patch(jgp),                       & !in
          &                  ext_data(jg)           ,            & !in
          &                  p_nh_state(jg)%prog(nnow(jg)) ,     & !inout
          &                  p_nh_state(jg)%prog(n_now_rcf) ,    & !inout
          &                  p_nh_state(jg)%prog(n_now_rcf) ,    & !inout
          &                  p_nh_state(jg)%diag,                & !inout
          &                  prm_diag  (jg),                     & !inout
          &                  prm_nwp_tend(jg)                ,   &
          &                  p_lnd_state(jg)%diag_lnd,           &
          &                  p_lnd_state(jg)%prog_lnd(n_now_rcf),& !inout
          &                  p_lnd_state(jg)%prog_lnd(n_now_rcf),& !inout
          &                  p_lnd_state(jg)%prog_wtr(n_now_rcf),& !inout
          &                  p_lnd_state(jg)%prog_wtr(n_now_rcf),& !inout
          &                  p_nh_state_lists(jg)%prog_list(n_now_rcf) ) !in

      ! This is to enforce another slow physics call at the end of the first time step
      lcall_phy(jg,:)     = .FALSE.
      t_elapsed_phy(jg,:) = dt_phy(jg,:) - dt_loc

      CASE (iecham) ! iforcing

        SELECT CASE (echam_phy_config%idcphycpl)

        CASE (1) ! idcphycpl

          ! echam physics, fast physics coupling
          ! the physics forcing in the dynamical core is zero
          p_nh_state(jg)%diag%ddt_exner_phy(:,:,:)   = 0._wp
          p_nh_state(jg)%diag%ddt_vn_phy(:,:,:)      = 0._wp
          prm_tend  (jg)%qtrc(:,:,:,:)               = 0._wp

        CASE (2) ! idcphycpl

          ! echam physics, slow physics coupling
          IF (ltimer) CALL timer_start(timer_iconam_echam)
          CALL interface_iconam_echam( dt_loc                         ,& !in
            &                          mtime_current                  ,& !in
            &                          p_patch(jg)                    ,& !in
            &                          p_int_state(jg)                ,& !in
            &                          p_nh_state(jg)%metrics         ,& !in
            &                          p_nh_state(jg)%prog(nnow(jg))  ,& !inout
            &                          p_nh_state(jg)%prog(n_now_rcf) ,& !inout
            &                          p_nh_state(jg)%diag            ,&
            &                          p_nh_state_lists(jg)%prog_list(n_now_rcf))
          IF (ltimer) CALL timer_stop(timer_iconam_echam)

        CASE DEFAULT ! idcphycpl

          CALL finish (routine, 'echam_phy_config%idcphycpl /= 1,2 currently not implemented')

        END SELECT ! idcphycpl

      END SELECT ! iforcing

    END IF ! is_les_phy

    ! Boundary interpolation of land state variables entering into radiation computation
    ! if a reduced grid is used in the child domain(s)
    DO jn = 1, p_patch(jg)%n_childdom

      jgc = p_patch(jg)%child_id(jn)
      IF (.NOT. p_patch(jgc)%ldom_active) CYCLE

      IF ( lredgrid_phys(jgc) ) THEN
        CALL interpol_rrg_grf(jg, jgc, jn, nnow_rcf(jg))
        IF (atm_phy_nwp_config(jgc)%latm_above_top) THEN
          CALL copy_rrg_ubc(jg, jgc)
        ENDIF
      ENDIF
    ENDDO

    IF (p_patch(jg)%n_childdom > 0) THEN

      dt_sub     = dt_loc/2._wp    ! dyn. time step on next refinement level

      DO jn = 1, p_patch(jg)%n_childdom

        jgc = p_patch(jg)%child_id(jn)
        IF (.NOT. p_patch(jgc)%ldom_active) CYCLE

        IF(p_patch(jgc)%n_patch_cells > 0) THEN
          IF(proc_split) CALL push_glob_comm(p_patch(jgc)%comm, p_patch(jgc)%proc0)
          CALL init_slowphysics( mtime_current, jgc, dt_sub )
          IF(proc_split) CALL pop_glob_comm()
        ENDIF

      ENDDO

    ENDIF

  END SUBROUTINE init_slowphysics

  !-------------------------------------------------------------------------
  !>
  !! Diagnostic computations for output - dynamics fields
  !!
  !! This routine encapsulates calls to diagnostic computations required at output
  !! times only
  !!
  !! @par Revision History
  !! Developed by Guenther Zaengl, DWD (2012-05-09)
  !!
  SUBROUTINE diag_for_output_dyn ()

!!$    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
!!$      &  routine = 'mo_nh_stepping:diag_for_output_dyn'

    ! Local variables
    INTEGER :: jg, jgc, jn ! loop indices

    REAL(wp), DIMENSION(:,:,:), POINTER  :: p_vn   => NULL()

    IF (ltimer) CALL timer_start(timer_nh_diagnostics)

    DO jg = 1, n_dom

      IF(p_patch(jg)%n_patch_cells == 0) CYCLE
      IF (.NOT. p_patch(jg)%ldom_active) CYCLE

      p_vn  => p_nh_state(jg)%prog(nnow(jg))%vn


      CALL rbf_vec_interpol_cell(p_vn,p_patch(jg),p_int_state(jg),&
                                 p_nh_state(jg)%diag%u,p_nh_state(jg)%diag%v)

      !CALL div(p_vn, p_patch(jg), p_int_state(jg), p_nh_state(jg)%diag%div)
      CALL div_avg(p_vn, p_patch(jg), p_int_state(jg),p_int_state(jg)%c_bln_avg,&
                                                          p_nh_state(jg)%diag%div)

      CALL rot_vertex (p_vn, p_patch(jg), p_int_state(jg), p_nh_state(jg)%diag%omega_z)

      ! Diagnose relative vorticity on cells
      CALL verts2cells_scalar(p_nh_state(jg)%diag%omega_z, p_patch(jg), &
        p_int_state(jg)%verts_aw_cells, p_nh_state(jg)%diag%vor)


      CALL diagnose_pres_temp (p_nh_state(jg)%metrics, p_nh_state(jg)%prog(nnow(jg)), &
        &                      p_nh_state(jg)%prog(nnow_rcf(jg)),                     &
        &                      p_nh_state(jg)%diag,p_patch(jg),                       &
        &                      opt_calc_temp=.TRUE.,                                  &
        &                      opt_calc_pres=.TRUE.                                   )

    ENDDO ! jg-loop

    ! Fill boundaries of nested domains
    DO jg = n_dom, 1, -1

      IF(p_patch(jg)%n_patch_cells == 0 .OR. p_patch(jg)%n_childdom == 0) CYCLE
      IF (.NOT. p_patch(jg)%ldom_active) CYCLE

      CALL sync_patch_array_mult(SYNC_C, p_patch(jg), 3, p_nh_state(jg)%diag%u,      &
        p_nh_state(jg)%diag%v, p_nh_state(jg)%diag%div)


      DO jn = 1, p_patch(jg)%n_childdom
        jgc = p_patch(jg)%child_id(jn)
        IF (.NOT. p_patch(jgc)%ldom_active) CYCLE

        CALL interpol_scal_grf (p_patch(jg), p_patch(jgc), p_grf_state(jg)%p_dom(jn), 3, &
             p_nh_state(jg)%diag%u, p_nh_state(jgc)%diag%u, p_nh_state(jg)%diag%v,       &
             p_nh_state(jgc)%diag%v, p_nh_state(jg)%diag%div, p_nh_state(jgc)%diag%div   )

      ENDDO

    ENDDO ! jg-loop

    IF (ltimer) CALL timer_stop(timer_nh_diagnostics)

  END SUBROUTINE diag_for_output_dyn



  !-------------------------------------------------------------------------
  !>
  !! Wrapper for computation of aggregated land variables
  !!
  !!
  !! @par Revision History
  !! Developed by Guenther Zaengl, DWD (2014-07-21)
  !!
  SUBROUTINE aggr_landvars

    ! Local variables
    INTEGER :: jg ! loop indices

    IF (ltimer) CALL timer_start(timer_nh_diagnostics)

    DO jg = 1, n_dom

      IF(p_patch(jg)%n_patch_cells == 0) CYCLE
      IF (.NOT. p_patch(jg)%ldom_active) CYCLE

      IF (  atm_phy_nwp_config(jg)%inwp_surface == 1 ) THEN
        CALL aggregate_landvars( p_patch(jg), ext_data(jg),                 &
             p_lnd_state(jg)%prog_lnd(nnow_rcf(jg)), p_lnd_state(jg)%diag_lnd)
      ENDIF

    ENDDO ! jg-loop

    IF (ltimer) CALL timer_stop(timer_nh_diagnostics)

  END SUBROUTINE aggr_landvars

  !-------------------------------------------------------------------------
  !>
  !! Fills nest boundary cells for physics fields
  !!
  !!
  !! @par Revision History
  !! Developed by Guenther Zaengl, DWD (2014-07-21)
  !!
  SUBROUTINE fill_nestlatbc_phys

    ! Local variables
    INTEGER :: jg, jgc, jn ! loop indices

    IF (ltimer) CALL timer_start(timer_nh_diagnostics)

    ! Fill boundaries of nested domains
    DO jg = n_dom, 1, -1

      IF(p_patch(jg)%n_patch_cells == 0 .OR. p_patch(jg)%n_childdom == 0) CYCLE
      IF (.NOT. p_patch(jg)%ldom_active) CYCLE

      CALL sync_patch_array(SYNC_C, p_patch(jg), p_nh_state(jg)%prog(nnow_rcf(jg))%tke)

      DO jn = 1, p_patch(jg)%n_childdom
        jgc = p_patch(jg)%child_id(jn)
        IF (.NOT. p_patch(jgc)%ldom_active) CYCLE

        CALL interpol_phys_grf(ext_data, jg, jgc, jn)

        IF (lfeedback(jgc) .AND. ifeedback_type==1) CALL feedback_phys_diag(jgc, jg)

        CALL interpol_scal_grf (p_patch(jg), p_patch(jgc), p_grf_state(jg)%p_dom(jn), 1, &
           p_nh_state(jg)%prog(nnow_rcf(jg))%tke, p_nh_state(jgc)%prog(nnow_rcf(jgc))%tke)

      ENDDO

    ENDDO ! jg-loop

    IF (ltimer) CALL timer_stop(timer_nh_diagnostics)

  END SUBROUTINE fill_nestlatbc_phys


  !-------------------------------------------------------------------------
  !>
  !! Update of vertical wind offcentering and divergence damping
  !!
  !! This routine handles the increased sound-wave damping (by increasing the vertical wind offcentering)
  !! and mixed second-order/fourth-order divergence damping during the initial spinup phase
  !!
  !! @par Revision History
  !! Developed by Guenther Zaengl, DWD (2013-06-04)
  !!
  SUBROUTINE update_spinup_damping(elapsed_time)

    REAL(wp), INTENT(IN) :: elapsed_time
    REAL(wp) :: time1, time2

    time1 = 1800._wp  ! enhanced damping during the first half hour of integration
    time2 = 7200._wp  ! linear decrease of enhanced damping until time2

    IF (elapsed_time <= time1) THEN ! apply slightly super-implicit weights
      divdamp_fac_o2 = 8._wp*divdamp_fac
    ELSE IF (elapsed_time <= time2) THEN ! linearly decrease minimum weights to 0.5
      divdamp_fac_o2 = 8._wp*divdamp_fac*(time2-elapsed_time)/(time2-time1)
    ELSE
      divdamp_fac_o2 = 0._wp
    ENDIF


  END SUBROUTINE update_spinup_damping


  !-------------------------------------------------------------------------
  !> Auxiliary routine to encapsulate initialization of exner_pr variable
  !!
  SUBROUTINE init_exner_pr(jg, nnow)

    INTEGER, INTENT(IN) :: jg   ! domain ID
    INTEGER, INTENT(IN) :: nnow ! time step indicator


#ifndef _OPENACC
!$OMP PARALLEL
#endif
    CALL copy(p_nh_state(jg)%prog(nnow)%exner-REAL(p_nh_state(jg)%metrics%exner_ref_mc,wp), &
         p_nh_state(jg)%diag%exner_pr)
#ifndef _OPENACC
!$OMP END PARALLEL
#endif

  END SUBROUTINE init_exner_pr

  !-------------------------------------------------------------------------
  !> Driver routine to reset the model to its initial state if IAU iteration is selected
  !!
  SUBROUTINE reset_to_initial_state(datetime_current)

    TYPE(datetime), POINTER :: datetime_current
    INTEGER :: jg

    nnow(:)     = 1
    nnow_rcf(:) = 1
    nnew(:)     = 2
    nnew_rcf(:) = 2

    WRITE(message_text,'(a)') 'Reset model to initial state, repeat IAU with full incrementation window'
    CALL message('',message_text)

    atm_phy_nwp_config(:)%lcalc_acc_avg = .FALSE.

    CALL restore_initial_state(p_patch(1:), p_nh_state, prm_diag, prm_nwp_tend, p_lnd_state, ext_data)

    DO jg=1, n_dom

      IF (.NOT. p_patch(jg)%ldom_active) CYCLE

      CALL diagnose_pres_temp (p_nh_state(jg)%metrics, p_nh_state(jg)%prog(nnow(jg)), &
        &                      p_nh_state(jg)%prog(nnow_rcf(jg)),                     &
        &                      p_nh_state(jg)%diag,p_patch(jg),                       &
        &                      opt_calc_temp=.TRUE.,                                  &
        &                      opt_calc_pres=.TRUE.                                   )

      CALL rbf_vec_interpol_cell(p_nh_state(jg)%prog(nnow(jg))%vn,p_patch(jg),p_int_state(jg),&
                                 p_nh_state(jg)%diag%u,p_nh_state(jg)%diag%v)

      CALL compute_airmass(p_patch(jg),                  &
        &                  p_nh_state(jg)%metrics,       &
        &                  p_nh_state(jg)%prog(nnow(jg)),&
        &                  p_nh_state(jg)%diag, itlev = 2)

      CALL init_exner_pr(jg, nnow(jg))

      CALL init_nwp_phy(                            &
           & p_patch(jg)                           ,&
           & p_nh_state(jg)%metrics                ,&
           & p_nh_state(jg)%prog(nnow(jg))         ,&
           & p_nh_state(jg)%diag                   ,&
           & prm_diag(jg)                          ,&
           & prm_nwp_tend(jg)                      ,&
           & p_lnd_state(jg)%prog_lnd(nnow_rcf(jg)),&
           & p_lnd_state(jg)%prog_lnd(nnew_rcf(jg)),&
           & p_lnd_state(jg)%prog_wtr(nnow_rcf(jg)),&
           & p_lnd_state(jg)%prog_wtr(nnew_rcf(jg)),&
           & p_lnd_state(jg)%diag_lnd              ,&
           & ext_data(jg)                          ,&
           & phy_params(jg)                        ,&
           & datetime_current)

    ENDDO

      CALL aggr_landvars

      CALL init_slowphysics (datetime_current, 1, dtime)

      CALL fill_nestlatbc_phys

  END SUBROUTINE reset_to_initial_state

  !-------------------------------------------------------------------------
  !> Control routine for adaptive number of dynamic substeps
  !!
  SUBROUTINE set_ndyn_substeps(lcfl_watch_mode)

    LOGICAL, INTENT(INOUT) :: lcfl_watch_mode

    INTEGER :: jg
    REAL(wp) :: mvcfl(n_dom)
    LOGICAL :: lskip

    lskip = .FALSE.

    mvcfl(1:n_dom) = p_nh_state(1:n_dom)%diag%max_vcfl_dyn

    p_nh_state(1:n_dom)%diag%max_vcfl_dyn = 0._vp

    mvcfl = global_max(mvcfl)
    IF (ANY(mvcfl(1:n_dom) > 0.85_wp) .AND. .NOT. lcfl_watch_mode) THEN
      WRITE(message_text,'(a)') 'High CFL number for vertical advection in dynamical core, entering watch mode'
      CALL message('',message_text)
      lcfl_watch_mode = .TRUE.
    ENDIF

    IF (lcfl_watch_mode) THEN
      DO jg = 1, n_dom
        IF (mvcfl(jg) > 0.95_wp .OR. ndyn_substeps_var(jg) > ndyn_substeps) THEN
          WRITE(message_text,'(a,i3,a,f7.4)') 'Maximum vertical CFL number in domain ', &
            jg,':', mvcfl(jg)
          CALL message('',message_text)
        ENDIF
        IF (mvcfl(jg) > 1.05_wp) THEN
          ndyn_substeps_var(jg) = MIN(ndyn_substeps_var(jg)+1,ndyn_substeps_max)
          advection_config(jg)%ivcfl_max = ndyn_substeps_var(jg)
          WRITE(message_text,'(a,i3,a,i3)') 'Number of dynamics substeps in domain ', &
            jg,' increased to ', ndyn_substeps_var(jg)
          CALL message('',message_text)
        ENDIF
        IF (ndyn_substeps_var(jg) > ndyn_substeps .AND.                                            &
            mvcfl(jg)*REAL(ndyn_substeps_var(jg),wp)/REAL(ndyn_substeps_var(jg)-1,wp) < 0.95_wp) THEN
          ndyn_substeps_var(jg) = ndyn_substeps_var(jg)-1
          advection_config(jg)%ivcfl_max = ndyn_substeps_var(jg)
          WRITE(message_text,'(a,i3,a,i3)') 'Number of dynamics substeps in domain ', &
            jg,' decreased to ', ndyn_substeps_var(jg)
          CALL message('',message_text)
          lskip = .TRUE.
        ENDIF
      ENDDO
    ENDIF

    IF (ALL(ndyn_substeps_var(1:n_dom) == ndyn_substeps) .AND. ALL(mvcfl(1:n_dom) < 0.8_wp) .AND. &
        lcfl_watch_mode .AND. .NOT. lskip) THEN
      WRITE(message_text,'(a)') 'CFL number for vertical advection has decreased, leaving watch mode'
      CALL message('',message_text)
      lcfl_watch_mode = .FALSE.
    ENDIF

  END SUBROUTINE set_ndyn_substeps

  !-------------------------------------------------------------------------
  !>
  !! Physics time control
  !!
  !! Time control for slow and fast physics. This function provides a 2D array
  !! of type LOGICAL. For each physical process there is one column of length
  !! n_dom. A physical process (iphys) on domain (jg) will be called, if
  !! lcall_phy(jg,iphys)=.TRUE.. Whether it is .TRUE. or .FALSE. depends
  !! on the current time and the prescribed calling period listed in
  !! dt_phy.
  !!
  !! @par Revision History
  !! Developed by Daniel Reinert, DWD (2010-09-23)
  !! Modification by Daniel Reinert (2011-12-14)
  !! - replaced zoo of fast physics time steps by a single time step, named
  !!   dt_fastphy.
  !!
  SUBROUTINE time_ctrl_physics ( dt_phy, dt_loc, jg, linit, t_elapsed_phy, lcall_phy )

    REAL(wp), INTENT(IN)    ::   &      !< Field of calling-time interval (seconds) for
      &  dt_phy(:,:)                    !< each domain and physical process

    LOGICAL, INTENT(IN)     ::   &      !< special initialization of lcall_phy and
      &  linit                          !< t_elapsed_phy for the first call of
                                        !< nwp_nh_interface before the first dynamcs step

    REAL(wp), INTENT(INOUT) ::   &      !< elapsed time after the last call of physics
      &  t_elapsed_phy(:,:)             !< packages

    LOGICAL, INTENT(OUT)    ::   &
      &  lcall_phy(:,:)

    REAL(wp), INTENT(IN) :: dt_loc      !< transport/physics time step

    INTEGER, INTENT(IN) :: jg           !< domain number

    INTEGER :: ip, ips                  !< loop index


  !-------------------------------------------------------------------------

    ! special treatment for the first physics call prior to the very first
    ! dynamics step
    IF (linit) THEN

      ! Initialize lcall_phy with .false. Only slow physics will be set to .TRUE. initially;
      ! turbulent transfer is called by specifying the initialization mode of the NWP interface
      lcall_phy(jg,:)  = .FALSE.

      ! slow physics
      IF ( atm_phy_nwp_config(jg)%lproc_on(itconv) ) lcall_phy(jg,itconv) = .TRUE.

      IF ( atm_phy_nwp_config(jg)%lproc_on(itccov) ) lcall_phy(jg,itccov) = .TRUE.

      IF ( atm_phy_nwp_config(jg)%lproc_on(itrad)  ) lcall_phy(jg,itrad)  = .TRUE.

      IF ( atm_phy_nwp_config(jg)%lproc_on(itradheat) ) lcall_phy(jg,itradheat) = .TRUE.

      IF ( atm_phy_nwp_config(jg)%lproc_on(itsso)  ) lcall_phy(jg,itsso)  = .TRUE.

      IF ( atm_phy_nwp_config(jg)%lproc_on(itgwd)  ) lcall_phy(jg,itgwd)  = .TRUE.

    ELSE
      !
      ! all physical processes are forced to run at a multiple of
      ! the advective time step. Note that fast physics are treated
      ! as a combined process in the case of t_elapsed_phy and dt_phy,
      ! but treated individually in the case of lcall_phy.
      !
      DO ips = 1, iphysproc_short

        ! If a physics package has been called at previous timestep,
        ! reset the time counter.
        !
        IF ( ANY(lcall_phy(jg,PACK(iproclist,map_phyproc(1:iphysproc,ips)))) ) THEN
           t_elapsed_phy(jg,ips)  = 0._wp
        ENDIF

        ! update time counter
        !
        t_elapsed_phy(jg,ips) = t_elapsed_phy(jg,ips) + dt_loc  ! transport/physics timestep !!


        ! The purpose of the 0.99999999 factor is to prevent pathological cases of
        ! truncation error accumulation; will become obsolete when changing to
        ! integer arithmetics for time control
        IF( t_elapsed_phy(jg,ips) >= 0.99999999_wp*dt_phy(jg,ips) ) THEN
          lcall_phy(jg,PACK(iproclist,map_phyproc(1:iphysproc,ips)))  = .TRUE.
        ELSE
          lcall_phy(jg,PACK(iproclist,map_phyproc(1:iphysproc,ips)))  = .FALSE.
        ENDIF

      ENDDO  ! ips


      ! In addition, it must be checked, whether the individual processes
      ! are switched on at all (lproc_on =.TRUE.). If not, lcall_phy is
      ! reset to false.
      !
      DO ip = 1, iphysproc   ! not that we have to loop over ALL processes
        IF (.NOT. atm_phy_nwp_config(jg)%lproc_on(ip) ) THEN
          lcall_phy(jg,ip)  = .FALSE.
        ENDIF
      ENDDO  ! ip

    ENDIF

  END SUBROUTINE time_ctrl_physics


  !-------------------------------------------------------------------------
  !>
  !! Setup of physics time control
  !!
  !! Setup of time control for slow and fast physics. The mapping matrix is
  !! initialized, which provides mapping rules required when mapping
  !! between variables of size iphysproc and iphysproc_short. Typical examples
  !! are lcall_phy and t_elapsed_phy.
  !!
  !! @par Revision History
  !! Developed by Daniel Reinert, DWD (2011-12-14)
  !!
  SUBROUTINE setup_time_ctrl_physics ( )

  !-------------------------------------------------------------------------

  ! list of physical processes (x-axis of mapping matrix)
    iproclist = (/ itconv,itccov,itrad,itsso,itgwd,itsatad,itturb,&
      &            itgscp,itsfc,itradheat /)

    map_phyproc(1:iphysproc,1:iphysproc_short) = .FALSE. ! initialization

    map_phyproc(1   ,1) = .TRUE.  ! simple one to one mapping
    map_phyproc(2   ,2) = .TRUE.  ! simple one to one mapping
    map_phyproc(3   ,3) = .TRUE.  ! simple one to one mapping
    map_phyproc(4   ,4) = .TRUE.  ! simple one to one mapping
    map_phyproc(5   ,5) = .TRUE.  ! simple one to one mapping
    map_phyproc(6:10,6) = .TRUE.  ! mapping of fast physics processes to single one

  END SUBROUTINE setup_time_ctrl_physics
  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  !>
  !! @par Revision History
  !!
  SUBROUTINE deallocate_nh_stepping(latbc)

    TYPE (t_latbc_data), INTENT(INOUT) :: latbc
  INTEGER                              ::  jg, ist

  !-----------------------------------------------------------------------
  !
  ! deallocate auxiliary fields for tracer transport and rcf
  !
  DO jg = 1, n_dom
    DEALLOCATE( prep_adv(jg)%mass_flx_me, prep_adv(jg)%mass_flx_ic,    &
      &         prep_adv(jg)%vn_traj, prep_adv(jg)%w_traj,             &
      &         prep_adv(jg)%topflx_tra, STAT=ist                      )
    IF (ist /= SUCCESS) THEN
      CALL finish ( modname//': perform_nh_stepping',            &
        &    'deallocation for mass_flx_me, mass_flx_ic, vn_traj,' // &
        &    'w_traj, topflx_tra failed' )
    ENDIF
  ENDDO

  DEALLOCATE( prep_adv, STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ( modname//': perform_nh_stepping',              &
      &    'deallocation for prep_adv failed' )
  ENDIF

  DEALLOCATE( jstep_adv, STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ( modname//': perform_nh_stepping',              &
      &    'deallocation for jstep_adv failed' )
  ENDIF

  !
  ! deallocate flow control variables
  !
  DEALLOCATE( lcall_phy, linit_dyn, t_elapsed_phy, STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ( modname//': perform_nh_stepping',          &
      &    'deallocation for lcall_phy, t_elapsed_phy ' //        &
      &    'failed' )
  ENDIF

  IF (l_limited_area .AND. latbc_config%itype_latbc > 0) THEN
    IF (num_prefetch_proc >= 1) THEN
      CALL latbc%finalize()
    ELSE
      CALL deallocate_latbc_data()
    ENDIF
  ENDIF

  END SUBROUTINE deallocate_nh_stepping
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !! @par Revision History
  !!
  SUBROUTINE allocate_nh_stepping
!
  INTEGER                              :: jg, jp !, nlen
  INTEGER                              :: ist
  CHARACTER(len=MAX_CHAR_LENGTH)       :: attname   ! attribute name
  TYPE(t_RestartAttributeList), POINTER :: restartAttributes

!-----------------------------------------------------------------------

  !
  ! allocate axiliary fields for transport
  !
  ALLOCATE(prep_adv(n_dom), STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ( modname//': perform_nh_stepping',           &
    &      'allocation for prep_adv failed' )
  ENDIF

  ALLOCATE(jstep_adv(n_dom), STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ( modname//': perform_nh_stepping',           &
    &      'allocation for jstep_adv failed' )
  ENDIF


  ! allocate flow control variables for transport and slow physics calls
  ALLOCATE(lcall_phy(n_dom,iphysproc), linit_dyn(n_dom), &
    &      t_elapsed_phy(n_dom,iphysproc_short),         &
    &      STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ( modname//': perform_nh_stepping',           &
    &      'allocation for flow control variables failed' )
  ENDIF
  !
  ! initialize
  restartAttributes => getAttributesForRestarting()
  IF (ASSOCIATED(restartAttributes)) THEN
    !
    ! Get t_elapsed_phy and lcall_phy from restart file
    DO jg = 1,n_dom
      WRITE(attname,'(a,i2.2)') 'ndyn_substeps_DOM',jg
      ndyn_substeps_var(jg) = restartAttributes%getInteger(TRIM(attname))
      WRITE(attname,'(a,i2.2)') 'jstep_adv_marchuk_order_DOM',jg
      jstep_adv(jg)%marchuk_order = restartAttributes%getInteger(TRIM(attname))

      DO jp = 1,iphysproc_short
        WRITE(attname,'(a,i2.2,a,i2.2)') 't_elapsed_phy_DOM',jg,'_PHY',jp
        t_elapsed_phy(jg,jp) = restartAttributes%getReal(TRIM(attname))
      ENDDO
      DO jp = 1,iphysproc
        WRITE(attname,'(a,i2.2,a,i2.2)') 'lcall_phy_DOM',jg,'_PHY',jp
        lcall_phy(jg,jp) = restartAttributes%getLogical(TRIM(attname))
      ENDDO
    ENDDO
    linit_dyn(:)      = .FALSE.
  ELSE
    jstep_adv(:)%marchuk_order = 0
    t_elapsed_phy(:,:)         = 0._wp
    linit_dyn(:)               = .TRUE.
  ENDIF

  DO jg=1, n_dom
    ALLOCATE(                                                                      &
      &  prep_adv(jg)%mass_flx_me (nproma,p_patch(jg)%nlev  ,p_patch(jg)%nblks_e), &
      &  prep_adv(jg)%mass_flx_ic (nproma,p_patch(jg)%nlevp1,p_patch(jg)%nblks_c), &
      &  prep_adv(jg)%vn_traj     (nproma,p_patch(jg)%nlev,  p_patch(jg)%nblks_e), &
      &  prep_adv(jg)%w_traj      (nproma,p_patch(jg)%nlevp1,p_patch(jg)%nblks_c), &
      &  prep_adv(jg)%topflx_tra  (nproma,p_patch(jg)%nblks_c,MAX(1,ntracer)),     &
      &       STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ( modname//': perform_nh_stepping',           &
      &      'allocation for mass_flx_me, mass_flx_ic, vn_traj, ' // &
      &      'w_traj, topflx_tra failed' )
    ENDIF
    !
    ! initialize (as long as restart output is synchroinzed with advection,
    ! these variables do not need to go into the restart file)
!$OMP PARALLEL
    CALL init(prep_adv(jg)%mass_flx_me)
    CALL init(prep_adv(jg)%mass_flx_ic)
    CALL init(prep_adv(jg)%vn_traj)
    CALL init(prep_adv(jg)%w_traj)
    CALL init(prep_adv(jg)%topflx_tra)
!$OMP END PARALLEL

  ENDDO

  IF (l_limited_area .AND. latbc_config%itype_latbc > 0 .AND. num_prefetch_proc == 0) THEN
    CALL prepare_latbc_data(p_patch(1), p_int_state(1), p_nh_state(1), ext_data(1))
  ENDIF

END SUBROUTINE allocate_nh_stepping
  !-----------------------------------------------------------------------------

END MODULE mo_nh_stepping

