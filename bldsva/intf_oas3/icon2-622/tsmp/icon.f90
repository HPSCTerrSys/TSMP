!>
!! This is the master program of the ICON model.
!!
!!
!! @par Revision History
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
PROGRAM icon

#if defined (__INTEL_COMPILER) || defined (__PGI) || defined (NAGFOR)
#ifdef VARLIST_INITIZIALIZE_WITH_NAN
  USE, INTRINSIC :: ieee_features
  USE, INTRINSIC :: ieee_arithmetic
  USE, INTRINSIC :: ieee_exceptions

  USE mo_kind, ONLY: wp
#endif
#endif
#if defined (__INTEL_COMPILER) && ! defined (VARLIST_INITIZIALIZE_WITH_NAN)
  USE, INTRINSIC :: ieee_arithmetic
#endif
  USE mo_exception,           ONLY: message_text, message, finish
  USE mo_io_units,            ONLY: filename_max
  USE mo_mpi,                 ONLY: start_mpi, stop_mpi, my_process_is_global_root
  USE mo_master_control,      ONLY: init_master_control,                                &
    &                               get_my_namelist_filename, get_my_process_type,      &
    &                               atmo_process, ocean_process, ps_radiation_process,  &
    &                               hamocc_process
#ifndef __NO_ICON_TESTBED__
  USE mo_master_control,      ONLY: testbed_process
#endif
  USE mo_time_config,         ONLY: time_config
  USE mtime,                  ONLY: OPERATOR(>)
  USE mo_util_signal
  USE mo_util_vcs,            ONLY: show_version

#ifndef __NO_ICON_OCEAN__
  USE mo_ocean_model,         ONLY: ocean_model
  USE mo_hamocc_model,        ONLY: hamocc_model  
#endif

#ifndef __NO_ICON_TESTBED__
  USE mo_icon_testbed,        ONLY: icon_testbed
#endif

#ifndef __NO_ICON_ATMO__
  USE mo_atmo_model,          ONLY: atmo_model
#ifdef COUP_OAS_ICON
  USE oas_icon_define
#endif
#endif

#ifndef __NO_ICON_PS_RAD__
  USE mo_ps_radiation_model, ONLY: ps_radiation_model
#endif

  IMPLICIT NONE

  INTEGER                     :: master_control_status, my_process_component
  CHARACTER(len=filename_max) :: my_namelist_filename
  CHARACTER(len=filename_max) :: master_namelist_filename="icon_master.namelist"

#if defined (__INTEL_COMPILER) || defined (__PGI) || defined (NAGFOR)
#ifdef VARLIST_INITIZIALIZE_WITH_NAN
  TYPE(ieee_status_type)      :: saved_fpscr
  LOGICAL                     :: halting_mode,  current_flags(size(ieee_all))
  REAL(wp)                    :: r
#endif
#endif

#if defined (__xlC__)
  INTEGER                     :: core_dump_flag
  INTEGER                     :: signals(1)
  INTEGER                     :: iret
#endif

  ! handling of comand-line arguments:
  TYPE t_cmdarg_option
    CHARACTER(len=1024) :: arg   !< (case-sensitive) option
    CHARACTER(len=1024) :: help  !< help string
  END TYPE t_cmdarg_option

  ENUM, BIND(C)
    ENUMERATOR :: ARG_UNKNOWN = 0, &
      &           ARG_HELP,        &
      &           ARG_VERSION
  END ENUM

  TYPE(t_cmdarg_option), PARAMETER :: cmdarg_options(2) =            &
    & [                                                              &
    &  t_cmdarg_option("--help",    "print this help message."),     &
    &  t_cmdarg_option("--version", "print version info and exit.")  &
    & ]

  INTEGER :: i,j
  LOGICAL :: lmatch, lmatch_ij, lcmdarg(0:SIZE(cmdarg_options))
  CHARACTER(len=1024) :: arg


!--------------------------------------------------------------------

#if defined (__INTEL_COMPILER) || defined (__PGI) || defined (NAGFOR)
#ifdef VARLIST_INITIZIALIZE_WITH_NAN
  CALL ieee_get_status(saved_fpscr)
  CALL ieee_set_halting_mode(ieee_all, .TRUE.)
#endif
#endif

#if defined (__INTEL_COMPILER) && ! defined (VARLIST_INITIZIALIZE_WITH_NAN)
  ! Important on Intel: disable underflow exceptions:
  CALL ieee_set_halting_mode(ieee_underflow, .FALSE.)
#endif

  !-------------------------------------------------------------------
  ! Initialize MPI, this should always be the first call
  CALL start_mpi('ICON')

  !-------------------------------------------------------------------
  !set up signal trapping on IBM: export USE_SIGNAL_HANDLING=yes

#if defined (__xlC__)
  core_dump_flag = 0
  signals(1)     = 0

  iret = signal_trap(core_dump_flag, signals)

  IF (iret == -2) THEN
    CALL message('', 'Signal trapping disabled by environment')
  ELSE IF (iret == -1) THEN
     WRITE(message_text,'(a,i0)') 'Error: ', iret
     CALL message('', message_text)
  ELSE IF (iret == 0) THEN
    CALL message('', 'FPE trapping is not set')
  ELSE
    WRITE(message_text,'(a,i0)') 'FPE trapping mode =', iret
    CALL message('', message_text)
  END IF
#endif
#if defined (__INTEL_COMPILER) || defined (__PGI) || defined (NAGFOR)
#ifdef VARLIST_INITIZIALIZE_WITH_NAN
  WRITE(message_text,'(a,l1)') ' IEEE standard supported: ', ieee_support_standard(r)
  CALL message('', message_text)
#endif
#endif

  ! print info on the current version:
  CALL show_version()

  ! When executing ICON, it is now possible to invoke command-line
  ! arguments (logical switches).
  i = 1
  lcmdarg(:) = .FALSE.
  DO
    CALL get_command_argument(i, arg)
    IF (LEN_TRIM(arg) == 0) EXIT

    lmatch = .FALSE.
    DO j=1,SIZE(cmdarg_options)
      lmatch_ij  = (TRIM(arg) == TRIM(cmdarg_options(j)%arg))
      lcmdarg(j) = lcmdarg(j) .OR. lmatch_ij
      lmatch     = lmatch     .OR. lmatch_ij
    END DO
    IF (.NOT. lmatch) THEN
      lcmdarg(ARG_UNKNOWN) = .TRUE.
      CALL message("", "command-line argument '"//TRIM(arg)//"' unknown!")      
    END IF
    i = i+1
  END DO

  ! print a list of available options and exit
  IF (ANY([lcmdarg(ARG_UNKNOWN), lcmdarg(ARG_HELP)])) THEN
    CALL message("", "")
    CALL message("", "list of available command-line options:")
    DO j=1,SIZE(cmdarg_options)
      WRITE(message_text,'(a,A20,a)') "    ", [CHARACTER(50) :: TRIM(cmdarg_options(j)%arg)], TRIM(cmdarg_options(j)%help)
      CALL message('', message_text)
    END DO
    CALL message("", "")
  END IF

  ! for '--version' ICON prints the same information as usually at the
  ! beginning of stdout and aborts then.
  IF (ANY([lcmdarg(ARG_UNKNOWN), lcmdarg(ARG_HELP), lcmdarg(ARG_VERSION)])) THEN
    CALL stop_mpi ! Shut down MPI
    STOP 'icon'
  END IF


  !-------------------------------------------------------------------
  ! Initialize the master control

  master_control_status = init_master_control(TRIM(master_namelist_filename))

  my_namelist_filename = get_my_namelist_filename()
  my_process_component = get_my_process_type()

  SELECT CASE (my_process_component)

#ifndef __NO_ICON_ATMO__
  CASE (atmo_process)
    CALL atmo_model  (my_namelist_filename, TRIM(master_namelist_filename))
#endif

#ifndef __NO_ICON_OCEAN__
  CASE (ocean_process)
    CALL ocean_model (my_namelist_filename, TRIM(master_namelist_filename))
#endif

#ifndef __NO_ICON_PS_RAD__
  CASE (ps_radiation_process)
    CALL ps_radiation_model  (my_namelist_filename, TRIM(master_namelist_filename))
#endif

#ifndef __NO_ICON_OCEAN__
  CASE (hamocc_process)
    CALL hamocc_model (my_namelist_filename, TRIM(master_namelist_filename))
#endif

#ifndef __NO_ICON_TESTBED__
  CASE (testbed_process)
    CALL icon_testbed(my_namelist_filename, TRIM(master_namelist_filename))
#endif


  CASE default
    CALL finish("icon","my_process_component is unknown")

  END SELECT

  IF (ASSOCIATED(time_config%tc_exp_stopdate) .AND. ASSOCIATED(time_config%tc_stopdate)) THEN
    ! write the control.status file
    IF (my_process_is_global_root()) THEN
      OPEN (500, FILE="finish.status")
      IF ((time_config%tc_exp_stopdate > time_config%tc_stopdate) .AND. time_config%tc_write_restart) THEN
        WRITE(500,*) "RESTART"
      ELSE
        WRITE(500,*) "OK"
      ENDIF
      CLOSE(500)
    END IF
  END IF

  ! Shut down MPI
  CALL stop_mpi

#if defined (__INTEL_COMPILER) || defined (__PGI) || defined (NAGFOR)
#ifdef VARLIST_INITIZIALIZE_WITH_NAN
  CALL ieee_set_status(saved_fpscr)
#endif
#endif

END PROGRAM icon
