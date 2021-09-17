!> Basic module initializing the MPI communication and handling most
!> of the MPI calls (wrapper functions).
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
!!
!!
!!  Processors are divided into
!!    1.    worker PEs    : majority of MPI tasks, doing the actual work
!!    2.    I/O PEs       : dedicated I/O server tasks        (only for parallel_nml::num_io_procs > 0)
!!    3.    one test PE   : for verification runs             (only for parallel_nml::p_test_run == .TRUE.)
!!    4.    restart PEs   : for asynchronous restart writing  (only for parallel_nml::num_restart_procs > 0)
!!    5.    prefetch PEs  : for prefetching of data           (only for parallel_nml::num_prefetch_proc > 0)
!!
!!
!!  The communicators are split like this:
!!
!!          0    p_work_pe0    p_io_pe0    p_restart_pe0    p_pref_pe0   process_mpi_all_size
!!
!!          |         |            |             |               |                |
!!          V         V            V             V               V                V
!!
!!          +---------------------------------------------------------------------+
!!          !                      process_mpi_all_comm                           !
!!          +---------+------------+-------------+---------------+----------------+
!!          | test PE | worker PEs |   I/O PEs   |  restart PEs  |  prefetch PEs  |
!!          +---------+------------+-------------+---------------+----------------+
!!          |    A    |     B      |     C       |      D        |       E        | p_comm_work
!!          |    A    |     A      |             |               |                | p_comm_work_test
!!          |    A    |     B      |     B       |               |                | p_comm_work_io
!!          |    A    |            |     B       |               |                | p_comm_io (B is worker PE 0 if num_io_procs == 0)
!!          |         |     A      |             |      A        |                | p_comm_work_restart
!!          |         |     A      |             |               |       A        | p_comm_work_pref
!!          +---------+------------+-------------+---------------+----------------+
!!
!!  Note that there are actually two different p_comm_work_io communicators:
!!  One that spans the worker AND I/O PEs as the NAME implies, AND one that IS ONLY defined on the test PE.
!!  Similarly, there IS another ghost communicator p_comm_io defined on the test PE.
!!  This has the consequence that `my_process_is_io()` IS NOT equivalent to `p_comm_io /= MPI_COMM_NULL`
!!
!!  Process groups with specific main procs, all of these are called from mo_atmo_model:
!!    * I/O (mo_name_list_output): name_list_io_main_proc()
!!    * restart (mo_async_restart): restart_main_proc()
!!    * prefetch (mo_async_latbc): prefetch_main_proc()
!!
!!
!!  List of MPI communicators:
!!  --------------------------
!!
!!       global_mpi_communicator
!!
!!         description  : MPI communicator spanning all PEs running  (= MPI_COMM_WORLD).
!!         size         : global_mpi_size
!!         this PE's ID : my_global_mpi_id (= get_my_global_mpi_id())
!!
!!
!!       process_mpi_all_comm  (= get_my_mpi_all_communicator())
!!
!!         description  : MPI communicator containing all PEs that are running this
!!                        model component. Different from global_mpi_communicator,
!!                        if master_nml::no_of_models > 1
!!         size         : process_mpi_all_size
!!         this PE's ID : my_process_mpi_all_id (= get_my_mpi_all_id())
!!
!!
!!       p_comm_work  (=get_my_mpi_work_communicator())
!!
!!         description  : MPI communicator for work group. On I/O and test PEs this
!!                        defaults to process_mpi_all_comm.
!!         size         : num_work_procs (on worker PEs: num_work_procs==p_n_work)
!!         this PE's ID : p_pe_work (= get_my_mpi_work_id())
!!
!!
!!       ... less important MPI communicators ...
!!
!!       p_comm_work_test (size = 0/1)
!!         description  : MPI communicator spanning work group and test PE
!!                        in verification mode parallel_nml::p_test_run == .TRUE.
!!       p_comm_work_io
!!         description  : MPI communicator spanning work group and I/O PEs
!!       p_comm_io
!!         description  : MPI communicator spanning I/O PEs
!!       p_comm_work_2_io
!!         description  : Inter(!)communicator work PEs - I/O PEs
!!       p_comm_input_bcast
!!         description  : MPI communicator for broadcasts in NetCDF input
!!       p_comm_work_pref
!!         description  : MPI Communicator spanning work group and prefetch PEs
!!       p_comm_work_2_pref
!!         description  : Inter(!)communicator work PEs - prefetching PEs
!!
!!
!!  Processor splitting
!!  -------------------
!!
!!       In order to improve the parallel load balancing, the domains
!!       of the first refinement level can be distributed to disjoint
!!       processor subsets. This distribution process is weighted by
!!       the namelist parameter grid_nml::patch_weight.  Since the
!!       processor sets are disjoint, MPI calls happen independently on
!!       each PE subset.
!!
!!       For global operations (sum, min, max in mo_sync), each patch
!!       contains additional data members which are initialized in
!!       mo_complete_subdivision::set_patch_communicators:
!!
!!       p_patch % comm
!!         description  : MPI communicator for this patch. If no processor
!!                        splitting enabled: p_patch%comm==p_comm_work.
!!         size         : p_patch%n_proc
!!         this PE's ID : p_patch%rank
!!
!!       p_patch % proc0
!!         description  : global id of processor with rank 0 (within the
!!                        working set p_comm_work)
!!
!!       The global communicator which is currently in use is stored by
!!       a call to mo_mpi::push_glob_comm in the data structure
!!       mo_mpi::glob_comm(0:max_lev), where the level number is equal
!!       to 0 for the global grid, 1 for the first generation of child
!!       domains etc.
!!
!!
!! Split of process_mpi_all_comm and stdio process
!! -----------------------------
!!    The process_mpi_all_comm is the whole model communicator and is split
!!    to test/work/io/restart, in this order
!!
!!    The process_mpi_stdio_id is always the 0 process of the process_mpi_all_comm
!!    If is p_test, then the testing process is also the stdio process
!!
!!    The my_process_is_mpi_workroot() is TRUE for the 0 process of the
!!    p_comm_work communicator, this communicator does not include the test process,
!!    the io and the restart processes. This communicator and and the mpi_workroot
!!    process should be used in the case of gather and scatter procedures that do
!!    not use the io process (the third part). This will be different from using
!!    the stdio process only in the case of p_test

MODULE mo_mpi

  ! Comment: Please use basic WRITE to nerr for messaging in the whole
  !          MPI package to achieve proper output.

  USE ISO_C_BINDING, ONLY: C_CHAR
  ! actual method (MPI-2)
#ifndef NOMPI
#if !defined (__SUNPRO_F95)
  USE mpi
#endif
#endif

#ifdef _OPENMP
  USE omp_lib, ONLY: omp_get_max_threads, omp_set_num_threads
#endif
#ifdef COUP_OAS_ICON
  USE oas_icon_define
#endif

  USE mo_kind
  USE mo_impl_constants, ONLY: SUCCESS
  USE mo_io_units,       ONLY: nerr

  IMPLICIT NONE

  PRIVATE                          ! all declarations are private

#ifndef NOMPI
#if defined (__SUNPRO_F95)
  INCLUDE "mpif.h"
#endif
#endif

  ! start/stop methods
  PUBLIC :: start_mpi, stop_mpi, abort_mpi

  ! split the global communicator to _process_mpi_communicator
  PUBLIC :: split_global_mpi_communicator
  !The given communicator will be the all communicator for this component
!   PUBLIC :: set_process_mpi_communicator
  ! Sets the test, work, i/o and prefetch communicators
  PUBLIC :: set_mpi_work_communicators
  ! Sets the p_comm_input_bcast
  PUBLIC :: set_comm_input_bcast
  ! set other parameters
  PUBLIC :: set_process_mpi_name

  PUBLIC :: push_glob_comm, pop_glob_comm, get_glob_proc0

  ! Logical functions
  PUBLIC :: run_is_global_mpi_parallel
  PUBLIC :: my_process_is_stdio, my_process_is_mpi_parallel, my_process_is_mpi_all_parallel
  PUBLIC :: my_process_is_mpi_seq, my_process_is_mpi_test, my_process_is_mpi_workroot
  PUBLIC :: my_process_is_mpi_ioroot
  PUBLIC :: my_process_is_mpi_all_seq, my_process_is_io
  PUBLIC :: my_process_is_global_root
  PUBLIC :: my_process_is_mpi_prefroot, my_process_is_pref
  PUBLIC :: my_process_is_restart, my_process_is_mpi_restartroot,my_process_is_work

  ! get parameters
  PUBLIC :: get_my_global_mpi_communicator   ! essentially a copy of MPI_COMM_WORLD   
  PUBLIC :: get_my_mpi_all_communicator   ! the communicator for the specific component, ie the process_mpi_all_comm
  PUBLIC :: get_my_mpi_all_comm_size   ! this is the the size of the communicator for the specific component
  PUBLIC :: get_my_mpi_work_communicator   ! the communicator for the workers of this component
  PUBLIC :: get_my_mpi_work_comm_size   ! this is the the size of the workers

  PUBLIC :: get_mpi_all_workroot_id, get_my_global_mpi_id, get_my_mpi_all_id
  PUBLIC :: get_mpi_all_ioroot_id, get_mpi_all_restartroot_id
  PUBLIC :: get_my_mpi_work_id, get_mpi_prefroot_id
  PUBLIC :: default_comm_type, null_comm_type


  ! some public communicators
  PUBLIC :: process_mpi_all_comm
  PUBLIC :: process_mpi_all_test_id, process_mpi_all_workroot_id, &
    &       process_mpi_all_ioroot_id, process_mpi_all_restartroot_id, &
    &       process_mpi_all_prefroot_id, p_comm_work_pref_compute_pe0

  PUBLIC :: p_comm_work, p_comm_work_test
  PUBLIC :: p_comm_work_2_io, p_comm_input_bcast, p_comm_work_io, &
    &       p_comm_io, p_comm_work_pref, p_comm_work_2_pref
  !restart communicators
  PUBLIC :: p_comm_work_2_restart, p_comm_work_restart
  PUBLIC :: p_communicator_a, p_communicator_b, p_communicator_d

  PUBLIC :: process_mpi_io_size, process_mpi_restart_size, process_mpi_pref_size

  PUBLIC :: process_mpi_stdio_id
  PUBLIC :: process_mpi_root_id
  PUBLIC :: process_work_io0

  ! Main communication methods
  PUBLIC :: global_mpi_barrier
  PUBLIC :: work_mpi_barrier

  PUBLIC :: p_comm_size
  PUBLIC :: p_comm_rank
  PUBLIC :: p_send, p_recv, p_sendrecv, p_bcast, p_barrier
! PUBLIC :: p_bcast_achar
  PUBLIC :: p_get_bcast_role
  PUBLIC :: p_isend, p_irecv, p_wait, p_wait_any,         &
    &       p_irecv_packed, p_send_packed, p_recv_packed, &
    &       p_bcast_packed,                               &
    &       p_pack_int, p_pack_bool, p_pack_real,         &
    &       p_pack_int_1d, p_pack_real_1d,                &
    &       p_pack_string, p_pack_real_2d,                &
    &       p_unpack_int, p_unpack_bool, p_unpack_real,   &
    &       p_unpack_int_1d, p_unpack_real_1d,            &
    &       p_unpack_string, p_unpack_real_2d, p_test
  PUBLIC :: p_scatter, p_gather, p_max, p_min, p_sum, p_global_sum, p_field_sum
  PUBLIC :: p_probe
  PUBLIC :: p_gatherv, p_allgather, p_allgatherv
  PUBLIC :: p_scatterv
  PUBLIC :: p_allreduce_max
  PUBLIC :: p_commit_type_struct
  PUBLIC :: p_alltoall
  PUBLIC :: p_alltoallv
  PUBLIC :: p_alltoallv_p2p
  PUBLIC :: p_clear_request
  PUBLIC :: p_mpi_wtime
  PUBLIC :: p_comm_is_intercomm
  PUBLIC :: p_comm_remote_size
  PUBLIC :: get_mpi_comm_world_ranks

  PUBLIC :: p_isEqual

  !----------- to be removed -----------------------------------------
  PUBLIC :: p_pe, p_io
  PUBLIC :: num_test_procs, num_work_procs,     &
       &    p_work_pe0, p_io_pe0,               &
       &    p_n_work, p_pe_work, p_restart_pe0, &
       &    p_pref_pe0
  !     p_test_pe,
  !--------------------------------------------------------------------

#ifndef NOMPI
#ifdef  __SUNPRO_F95
  INTEGER,PUBLIC :: MPI_INTEGER, MPI_STATUS_SIZE, MPI_SUCCESS, &
            MPI_INFO_NULL, MPI_ADDRESS_KIND, MPI_COMM_SELF, &
            MPI_UNDEFINED, mpi_max, mpi_in_place
#else
  PUBLIC :: MPI_INTEGER, MPI_STATUS_SIZE, MPI_SUCCESS, &
            MPI_INFO_NULL, MPI_ADDRESS_KIND, MPI_COMM_SELF, &
            MPI_UNDEFINED, mpi_max, mpi_in_place, mpi_op_null, &
            mpi_datatype_null, mpi_sum, mpi_lor
#endif
  PUBLIC :: MPI_2INTEGER
#endif
  PUBLIC :: MPI_ANY_SOURCE, MPI_COMM_NULL

  ! real data type matching real type of MPI implementation
  PUBLIC :: p_real_dp, p_real_sp
  PUBLIC :: p_int
  PUBLIC :: p_int_i8
  PUBLIC :: p_bool
  PUBLIC :: p_address_kind
  PUBLIC :: p_real_dp_byte, p_real_sp_byte, p_int_byte
  PUBLIC :: p_int_i4_byte, p_int_i8_byte

  PUBLIC ::get_mpi_time,ellapsed_mpi_time

  ! old fashioned method (MPI-1)

!!$#ifndef NOMPI
!!$  INCLUDE 'mpif.h'
!!$#endif

  ! general run time information

#ifndef NOMPI
  INTEGER :: version, subversion   ! MPI version
 ! MPI call inherent variables
  INTEGER :: p_error                     ! MPI error number
  INTEGER :: p_status(MPI_STATUS_SIZE)   ! standard information of MPI_RECV

  INTEGER, ALLOCATABLE, SAVE :: p_request(:) ! request values for non blocking calls
  INTEGER :: p_irequest ! the first p_irequest values of p_request are in use
  INTEGER :: p_mrequest ! actual size of p_request
  INTEGER, PARAMETER :: p_request_alloc_size = 4096
  INTEGER, PARAMETER :: p_address_kind = MPI_ADDRESS_KIND
#else
  INTEGER, PARAMETER :: p_address_kind = i8    ! should not get touched at all
  INTEGER, PARAMETER :: MPI_COMM_NULL  = 0
  ! dummy arguments for function calls:
  INTEGER, PARAMETER :: MPI_ANY_SOURCE = 0
#endif

  ! public parallel run information
  CHARACTER(len=64) :: global_mpi_name
  CHARACTER(len=64) :: process_mpi_name
  ! Do not change the following parameters
  INTEGER, PARAMETER :: process_mpi_stdio_id = 0
  INTEGER, PARAMETER :: process_mpi_root_id = 0
  INTEGER, PARAMETER :: default_comm_type = 1
  INTEGER, PARAMETER :: null_comm_type = 0

  ! communicator sets
  ! this is the global communicator
  INTEGER :: global_mpi_communicator  ! replaces MPI_COMM_WORLD
  INTEGER :: global_mpi_size          ! total number of processes in global world
  INTEGER :: my_global_mpi_id         ! process id in global world
  LOGICAL :: is_global_mpi_parallel

  ! this is the communicator for the whole component (atmo/ocean/etc)
  INTEGER :: process_mpi_all_comm     ! communicator in the whole model-component
  INTEGER :: process_mpi_all_size     ! total number of processes in the whole model-component
  INTEGER :: my_process_mpi_all_id        ! the id in the whole model communicator
  INTEGER :: process_mpi_all_workroot_id  ! the root process in component
  INTEGER :: process_mpi_all_ioroot_id    ! the first I/O process
  INTEGER :: process_mpi_all_test_id  ! the test process in component
  INTEGER :: process_mpi_all_restartroot_id ! the id of the first restart output process
  INTEGER :: process_work_io0
  INTEGER :: process_mpi_all_prefroot_id ! the id of the first prefetch process
  INTEGER :: p_comm_work_pref_compute_pe0 ! the ID for Communicator spanning work group and prefetch PEs

  LOGICAL :: process_is_mpi_parallel
  LOGICAL :: process_is_stdio
  LOGICAL :: is_mpi_test_run = .false.
  LOGICAL :: is_openmp_test_run = .false.

  ! this is the local work communicator (computation, i/o, etc)
!   INTEGER :: process_mpi_local_comm     ! communicator in the work group
!   INTEGER :: process_mpi_local_size     ! total number of processes in the whole model-component
!   INTEGER :: my_process_mpi_local_id

  INTEGER :: my_mpi_function  ! test, work, i/o, restart_output or prefetch
  INTEGER, PARAMETER :: test_mpi_process = 1
  INTEGER, PARAMETER :: work_mpi_process = 2
  INTEGER, PARAMETER :: io_mpi_process = 3
  INTEGER, PARAMETER :: restart_mpi_process = 4
  INTEGER, PARAMETER :: pref_mpi_process = 5

  !------------------------------------------------------------
  ! Processor distribution:
  ! num_test_procs:      0 or 1
  ! num_work_procs:      number of procs running in parallel on the model
  ! process_mpi_io_size: number of procs for I/O
  ! num_restart_procs:   number of procs used for writing restart files
  ! num_prefetch_proc:      number of procs used for prefetching of data
  ! num_test_procs + num_work_procs + process_mpi_io_size + num_restart_procs + num_prefetch_proc = process_mpi_all_size
  INTEGER :: num_test_procs
  INTEGER :: num_work_procs
  INTEGER :: process_mpi_io_size
  INTEGER :: process_mpi_restart_size
  INTEGER :: process_mpi_pref_size

  ! Note: p_test_pe, p_work_pe0, p_io_pe0 are identical on all PEs

  ! In a verification run, p_test_pe is the number of the PE running the complete model,
  ! otherwise it contains -1

  INTEGER :: p_test_pe     ! Number of test PE
  INTEGER :: p_work_pe0    ! Number of workgroup PE 0 within all PEs
  INTEGER :: p_io_pe0      ! Number of I/O PE 0 within all PEs (process_mpi_all_size if no I/O PEs)
  INTEGER :: p_restart_pe0 ! Number of Restart Output PE 0 within all PEs
                           ! (process_mpi_all_size if no restart PEs)
  INTEGER :: p_pref_pe0    ! Number of the prefetching PE 0 within all PEs (process_mpi_all_size
  !                          if no prefetching PEs)

  ! Note: p_n_work, p_pe_work are NOT identical on all PEs

  ! p_n_work: Number of PEs working together:
  ! - num_work_procs    for non-verification runs
  ! - num_work_procs    for verification runs on pes != p_test_pe
  ! - 1                 for verification runs on p_test_pe
  ! - num_io_procs      always on I/O pes
  ! - num_restart_procs on restart PEs

  INTEGER :: p_n_work
  INTEGER :: p_pe_work        ! PE number within work group

  ! MPI communicators
  INTEGER :: p_comm_work           ! Communicator for work group
  INTEGER :: p_comm_work_test      ! Communicator spanning work group and test PE
  INTEGER :: p_comm_work_io        ! Communicator spanning work group and I/O PEs
  INTEGER :: p_comm_io             ! Communicator spanning the I/O PEs
  INTEGER :: p_comm_work_2_io      ! Inter(!)communicator work PEs - I/O PEs
  INTEGER :: p_comm_input_bcast    ! Communicator for broadcasts in NetCDF input
  INTEGER :: p_comm_work_restart   ! Communicator spanning work group and Restart Output PEs
  INTEGER :: p_comm_work_2_restart ! Inter(!)communicator work PEs - Restart PEs
  INTEGER :: p_comm_work_pref           ! Communicator spanning work group and prefetch PEs
  INTEGER :: p_comm_work_2_pref    ! Inter(!)communicator work PEs - prefetching PEs

  INTEGER :: p_communicator_a ! for Set A
  INTEGER :: p_communicator_b ! for Set B
  INTEGER :: p_communicator_d ! for debug node

  INTEGER :: p_pe     = 0     ! this is the PE number of this task
  INTEGER :: p_io     = 0     ! PE number of PE handling IO

! non blocking calls

  ! module intrinsic names

!!#ifndef NOMPI
!!LK  INTEGER :: iope                  ! PE able to do IO
!!#endif

#ifdef DEBUG
  INTEGER :: nbcast                ! counter for broadcasts for debugging
#endif

  ! MPI native transfer types

  INTEGER :: p_int     = 0
  INTEGER :: p_real    = 0
  INTEGER :: p_bool    = 0
  INTEGER :: p_char    = 0

  ! MPI transfer types

  INTEGER :: p_real_dp = 0
  INTEGER :: p_real_sp = 0
  INTEGER :: p_int_i4  = 0
  INTEGER :: p_int_i8  = 0

  ! MPI native transfer types byte size

  INTEGER :: p_int_byte     = 0
  INTEGER :: p_real_byte    = 0
!  INTEGER :: p_bool_byte    = 0
!  INTEGER :: p_char_byte    = 0

  ! MPI transfer types byte size

  INTEGER :: p_real_dp_byte = 0
  INTEGER :: p_real_sp_byte = 0
  INTEGER :: p_int_i4_byte  = 0
  INTEGER :: p_int_i8_byte  = 0

  CHARACTER(len=256) :: message_text = ''

  ! Flag if processor splitting is active
  LOGICAL, PUBLIC :: proc_split = .FALSE.
#ifdef _OPENACC
  LOGICAL, PUBLIC :: i_am_accel_node = .FALSE.
#endif

  ! communicator stack for global sums
  INTEGER, PARAMETER :: max_lev = 10 ! 2 is sufficient
  INTEGER, PUBLIC :: comm_lev = 0, glob_comm(0:max_lev), comm_proc0(0:max_lev)

  ! Storage for buffered non-blocking point-to-point communication. This is
  ! especially needed for communication between I/O servers. These tasks may
  ! be completely asynchronous and therefore ISENDs may be launched before a
  ! corresponding IRECVs has been issued.
  INTEGER :: mpi_buffer(10000)


  ! define generic interfaces to allow proper compiling with picky compilers
  ! like NAG f95 for clean argument checking and shortening the call sequence.

  INTERFACE p_send
     MODULE PROCEDURE p_send_char
     MODULE PROCEDURE p_send_real
     MODULE PROCEDURE p_send_sreal
     MODULE PROCEDURE p_send_int
     MODULE PROCEDURE p_send_bool
     MODULE PROCEDURE p_send_real_1d
     MODULE PROCEDURE p_send_int_1d
     MODULE PROCEDURE p_send_bool_1d
     MODULE PROCEDURE p_send_char_1d
     MODULE PROCEDURE p_send_real_2d
     MODULE PROCEDURE p_send_int_2d
     MODULE PROCEDURE p_send_bool_2d
     MODULE PROCEDURE p_send_real_3d
     MODULE PROCEDURE p_send_int_3d
     MODULE PROCEDURE p_send_bool_3d
     MODULE PROCEDURE p_send_real_4d
     MODULE PROCEDURE p_send_int_4d
     MODULE PROCEDURE p_send_bool_4d
     MODULE PROCEDURE p_send_real_5d
  END INTERFACE

  INTERFACE p_isend
     MODULE PROCEDURE p_isend_char
     MODULE PROCEDURE p_isend_real
     MODULE PROCEDURE p_isend_sreal
     MODULE PROCEDURE p_isend_int
     MODULE PROCEDURE p_isend_bool
     MODULE PROCEDURE p_isend_real_1d
     MODULE PROCEDURE p_isend_int_1d
     MODULE PROCEDURE p_isend_bool_1d
     MODULE PROCEDURE p_isend_real_2d
     MODULE PROCEDURE p_isend_sreal_2d
     MODULE PROCEDURE p_isend_int_2d
     MODULE PROCEDURE p_isend_bool_2d
     MODULE PROCEDURE p_isend_real_3d
     MODULE PROCEDURE p_isend_int_3d
     MODULE PROCEDURE p_isend_bool_3d
     MODULE PROCEDURE p_isend_real_4d
     MODULE PROCEDURE p_isend_int_4d
     MODULE PROCEDURE p_isend_bool_4d
     MODULE PROCEDURE p_isend_real_5d
  END INTERFACE

  INTERFACE p_recv
     MODULE PROCEDURE p_recv_char
     MODULE PROCEDURE p_recv_real
     MODULE PROCEDURE p_recv_sreal
     MODULE PROCEDURE p_recv_int
     MODULE PROCEDURE p_recv_bool
     MODULE PROCEDURE p_recv_real_1d
     MODULE PROCEDURE p_recv_int_1d
     MODULE PROCEDURE p_recv_bool_1d
     MODULE PROCEDURE p_recv_char_1d
     MODULE PROCEDURE p_recv_real_2d
     MODULE PROCEDURE p_recv_int_2d
     MODULE PROCEDURE p_recv_bool_2d
     MODULE PROCEDURE p_recv_real_3d
     MODULE PROCEDURE p_recv_int_3d
     MODULE PROCEDURE p_recv_bool_3d
     MODULE PROCEDURE p_recv_real_4d
     MODULE PROCEDURE p_recv_int_4d
     MODULE PROCEDURE p_recv_bool_4d
     MODULE PROCEDURE p_recv_real_5d
  END INTERFACE

  INTERFACE p_irecv
     MODULE PROCEDURE p_irecv_char
     MODULE PROCEDURE p_irecv_real
     MODULE PROCEDURE p_irecv_sreal
     MODULE PROCEDURE p_irecv_int
     MODULE PROCEDURE p_irecv_bool
     MODULE PROCEDURE p_irecv_real_1d
     MODULE PROCEDURE p_irecv_int_1d
     MODULE PROCEDURE p_irecv_bool_1d
     MODULE PROCEDURE p_irecv_real_2d
     MODULE PROCEDURE p_irecv_sreal_2d
     MODULE PROCEDURE p_irecv_int_2d
     MODULE PROCEDURE p_irecv_bool_2d
     MODULE PROCEDURE p_irecv_real_3d
     MODULE PROCEDURE p_irecv_int_3d
     MODULE PROCEDURE p_irecv_bool_3d
     MODULE PROCEDURE p_irecv_real_4d
     MODULE PROCEDURE p_irecv_int_4d
     MODULE PROCEDURE p_irecv_bool_4d
  END INTERFACE

  INTERFACE p_sendrecv
     MODULE PROCEDURE p_sendrecv_real_1d
     MODULE PROCEDURE p_sendrecv_real_2d
     MODULE PROCEDURE p_sendrecv_real_3d
     MODULE PROCEDURE p_sendrecv_real_4d
     MODULE PROCEDURE p_sendrecv_char_array
  END INTERFACE

  INTERFACE p_bcast
     MODULE PROCEDURE p_bcast_real
     MODULE PROCEDURE p_bcast_real_single
     MODULE PROCEDURE p_bcast_int_i4
     MODULE PROCEDURE p_bcast_int_i8
     MODULE PROCEDURE p_bcast_bool
     MODULE PROCEDURE p_bcast_real_1d
     MODULE PROCEDURE p_bcast_real_1d_single
     MODULE PROCEDURE p_bcast_int_1d
     MODULE PROCEDURE p_bcast_int_i8_1d
     MODULE PROCEDURE p_bcast_bool_1d
     MODULE PROCEDURE p_bcast_real_2d
     MODULE PROCEDURE p_bcast_real_2d_single
     MODULE PROCEDURE p_bcast_int_2d
     MODULE PROCEDURE p_bcast_bool_2d
     MODULE PROCEDURE p_bcast_real_3d
     MODULE PROCEDURE p_bcast_int_3d
     MODULE PROCEDURE p_bcast_bool_3d
     MODULE PROCEDURE p_bcast_real_4d
     MODULE PROCEDURE p_bcast_int_4d
     MODULE PROCEDURE p_bcast_bool_4d
     MODULE PROCEDURE p_bcast_real_5d
     MODULE PROCEDURE p_bcast_int_7d
     MODULE PROCEDURE p_bcast_real_7d
     MODULE PROCEDURE p_bcast_char
     MODULE PROCEDURE p_bcast_cchar
     MODULE PROCEDURE p_bcast_char_1d
  END INTERFACE

  INTERFACE p_scatter
     MODULE PROCEDURE p_scatter_real_1d1d
     MODULE PROCEDURE p_scatter_single_1d1d
     MODULE PROCEDURE p_scatter_int_1d1d
  END INTERFACE

  INTERFACE p_gather
     MODULE PROCEDURE p_gather_real_0d1d
     MODULE PROCEDURE p_gather_real_1d2d
     MODULE PROCEDURE p_gather_real_5d6d
     MODULE PROCEDURE p_gather_real_1d1d
     MODULE PROCEDURE p_gather_int_0d1d
     MODULE PROCEDURE p_gather_int_1d1d
     MODULE PROCEDURE p_gather_char_0d1d
  END INTERFACE

  INTERFACE p_allgather
     MODULE PROCEDURE p_allgather_int_0d1d
     MODULE PROCEDURE p_allgather_int_1d2d
  END INTERFACE

  INTERFACE p_allgatherv
     MODULE PROCEDURE p_allgatherv_real_1d
     MODULE PROCEDURE p_allgatherv_int_1d
  END INTERFACE

  INTERFACE p_scatterv
    MODULE PROCEDURE p_scatterv_real1D2D
    MODULE PROCEDURE p_scatterv_real1D1D
    MODULE PROCEDURE p_scatterv_single1D1D
  END INTERFACE

  INTERFACE p_gatherv
    MODULE PROCEDURE p_gatherv_int
    MODULE PROCEDURE p_gatherv_real2D1D
    MODULE PROCEDURE p_gatherv_real3D1D
    MODULE PROCEDURE p_gatherv_int2D1D
    MODULE PROCEDURE p_gatherv_real2D2D
    MODULE PROCEDURE p_gatherv_sreal2D2D
    MODULE PROCEDURE p_gatherv_int2D2D
  END INTERFACE

  INTERFACE p_max
     MODULE PROCEDURE p_max_0d
     MODULE PROCEDURE p_max_int_0d
     MODULE PROCEDURE p_max_1d
     MODULE PROCEDURE p_max_int_1d
     MODULE PROCEDURE p_max_2d
     MODULE PROCEDURE p_max_3d
  END INTERFACE

  INTERFACE p_min
     MODULE PROCEDURE p_min_0d
     MODULE PROCEDURE p_min_int_0d
     MODULE PROCEDURE p_min_1d
     MODULE PROCEDURE p_min_int_1d
     MODULE PROCEDURE p_min_2d
     MODULE PROCEDURE p_min_3d
  END INTERFACE

  INTERFACE p_sum
     MODULE PROCEDURE p_sum_dp_0s
     MODULE PROCEDURE p_sum_dp_0d
     MODULE PROCEDURE p_sum_dp_1d
     MODULE PROCEDURE p_sum_i8_1d
     MODULE PROCEDURE p_sum_i_1d
     MODULE PROCEDURE p_sum_i_0d
  END INTERFACE

  INTERFACE p_global_sum
     MODULE PROCEDURE p_global_sum_1d
  END INTERFACE

  INTERFACE p_field_sum
     MODULE PROCEDURE p_field_sum_1d
  END INTERFACE

  !> generic interface for MPI communication calls
  INTERFACE p_allreduce_max
    MODULE PROCEDURE p_allreduce_max_int_1d
  END INTERFACE

  INTERFACE p_alltoall
    MODULE PROCEDURE p_alltoall_int
  END INTERFACE

  INTERFACE p_alltoallv
    MODULE PROCEDURE p_alltoallv_int
    MODULE PROCEDURE p_alltoallv_real_2d
    MODULE PROCEDURE p_alltoallv_sreal_2d
    MODULE PROCEDURE p_alltoallv_int_2d
  END INTERFACE

  INTERFACE p_alltoallv_p2p
    MODULE PROCEDURE p_alltoallv_p2p_real_2d
    MODULE PROCEDURE p_alltoallv_p2p_int_2d
  END INTERFACE

  INTERFACE p_isEqual
    MODULE PROCEDURE p_isEqual_int
    MODULE PROCEDURE p_isEqual_charArray
  END INTERFACE

  CHARACTER(*), PARAMETER :: modname = "mo_mpi"

#if defined( _OPENACC )
#define ACC_DEBUG NOACC
#if defined(__MPI_NOACC)
  LOGICAL, PARAMETER ::  acc_on = .FALSE.
#else
  LOGICAL, PARAMETER ::  acc_on = .TRUE.
#endif
#endif

CONTAINS

  !------------------------------------------------------------------------------
  REAL(dp) FUNCTION get_mpi_time()
#ifdef NOMPI
    get_mpi_time = 0.0_dp
#else
    get_mpi_time = MPI_Wtime()
#endif
  END FUNCTION get_mpi_time
  !------------------------------------------------------------------------------

  !------------------------------------------------------------------------------
  REAL(dp) FUNCTION ellapsed_mpi_time(start_time)
    REAL(dp), INTENT(in) :: start_time
#ifdef NOMPI
    ellapsed_mpi_time = 0.0_dp
#else
    ellapsed_mpi_time = MPI_Wtime() - start_time
#endif
  END FUNCTION ellapsed_mpi_time
  !------------------------------------------------------------------------------

  !------------------------------------------------------------------------------
  SUBROUTINE set_process_mpi_name(name)
    CHARACTER(len=*), INTENT(in) ::name

    process_mpi_name = TRIM(name)

  END SUBROUTINE set_process_mpi_name
  !------------------------------------------------------------------------------

  !------------------------------------------------------------------------------
  ! Rene: this has become obsolete and should be erased.
  ! It is used to idendify globally the proccess that caused an error
  INTEGER FUNCTION get_my_global_mpi_id()
    get_my_global_mpi_id = my_global_mpi_id
  END FUNCTION get_my_global_mpi_id
  !------------------------------------------------------------------------------

  !------------------------------------------------------------------------------
  INTEGER FUNCTION get_my_global_mpi_communicator()
    get_my_global_mpi_communicator = global_mpi_communicator
  END FUNCTION get_my_global_mpi_communicator
  !------------------------------------------------------------------------------

  !------------------------------------------------------------------------------
  INTEGER FUNCTION get_my_mpi_all_communicator()
    get_my_mpi_all_communicator = process_mpi_all_comm
  END FUNCTION get_my_mpi_all_communicator
  !------------------------------------------------------------------------------

  !------------------------------------------------------------------------------
  INTEGER FUNCTION get_my_mpi_all_comm_size()
    get_my_mpi_all_comm_size = process_mpi_all_size
  END FUNCTION get_my_mpi_all_comm_size
  !------------------------------------------------------------------------------

  !------------------------------------------------------------------------------
  INTEGER FUNCTION get_my_mpi_all_id()
    get_my_mpi_all_id = my_process_mpi_all_id
  END FUNCTION get_my_mpi_all_id
  !------------------------------------------------------------------------------

  !------------------------------------------------------------------------------
  INTEGER FUNCTION get_mpi_all_workroot_id()
    get_mpi_all_workroot_id = process_mpi_all_workroot_id
  END FUNCTION get_mpi_all_workroot_id
  !------------------------------------------------------------------------------

  !------------------------------------------------------------------------------
  INTEGER FUNCTION get_mpi_all_ioroot_id()
    get_mpi_all_ioroot_id = process_mpi_all_ioroot_id
  END FUNCTION get_mpi_all_ioroot_id
  !------------------------------------------------------------------------------

  !------------------------------------------------------------------------------
  INTEGER FUNCTION get_mpi_all_restartroot_id()
    get_mpi_all_restartroot_id = process_mpi_all_restartroot_id
  END FUNCTION get_mpi_all_restartroot_id
  !------------------------------------------------------------------------------

  !------------------------------------------------------------------------------
  INTEGER FUNCTION get_mpi_prefroot_id()
    get_mpi_prefroot_id = process_mpi_all_prefroot_id
  END FUNCTION get_mpi_prefroot_id
  !------------------------------------------------------------------------------

  !------------------------------------------------------------------------------
  INTEGER FUNCTION get_my_mpi_work_communicator()
    get_my_mpi_work_communicator = p_comm_work
  END FUNCTION get_my_mpi_work_communicator
  !------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  INTEGER FUNCTION get_my_mpi_work_comm_size()
    get_my_mpi_work_comm_size = p_n_work
  END FUNCTION get_my_mpi_work_comm_size
  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  INTEGER FUNCTION get_my_mpi_work_id()
    get_my_mpi_work_id = p_pe_work
  END FUNCTION get_my_mpi_work_id
  !------------------------------------------------------------------------------



  !------------------------------------------------------------------------------
  LOGICAL FUNCTION my_process_is_stdio()
#ifdef NOMPI
    my_process_is_stdio = .TRUE.
#else
    my_process_is_stdio = process_is_stdio
#endif
  END FUNCTION my_process_is_stdio
  !------------------------------------------------------------------------------

  !------------------------------------------------------------------------------
  LOGICAL FUNCTION my_process_is_io()
     my_process_is_io = (my_mpi_function == io_mpi_process)
  END FUNCTION my_process_is_io
  !------------------------------------------------------------------------------

  !------------------------------------------------------------------------------
  !>
  LOGICAL FUNCTION my_process_is_mpi_ioroot()
    my_process_is_mpi_ioroot = (my_process_mpi_all_id == process_mpi_all_ioroot_id)
  END FUNCTION my_process_is_mpi_ioroot
  !------------------------------------------------------------------------------

  !------------------------------------------------------------------------------
  !>
  LOGICAL FUNCTION my_process_is_global_root()
    my_process_is_global_root = (my_global_mpi_id == process_mpi_stdio_id) ! == 0
  END FUNCTION my_process_is_global_root
  !------------------------------------------------------------------------------

  !------------------------------------------------------------------------------
  !>
  LOGICAL FUNCTION my_process_is_restart()
     my_process_is_restart = (my_mpi_function == restart_mpi_process)
  END FUNCTION my_process_is_restart
  !------------------------------------------------------------------------------

  !------------------------------------------------------------------------------
  !>
  LOGICAL FUNCTION my_process_is_work()
     my_process_is_work = (my_mpi_function == work_mpi_process)
  END FUNCTION my_process_is_work
  !------------------------------------------------------------------------------

  !------------------------------------------------------------------------------
  !>
  LOGICAL FUNCTION my_process_is_mpi_restartroot()
    my_process_is_mpi_restartroot = (my_process_is_restart() .AND. &
     &                   (my_process_mpi_all_id == process_mpi_all_restartroot_id))
  END FUNCTION my_process_is_mpi_restartroot
  !------------------------------------------------------------------------------

  !------------------------------------------------------------------------------
  LOGICAL FUNCTION my_process_is_pref()
    my_process_is_pref = (my_mpi_function == pref_mpi_process)
  END FUNCTION my_process_is_pref
  !------------------------------------------------------------------------------

  !------------------------------------------------------------------------------
  !>
  LOGICAL FUNCTION my_process_is_mpi_prefroot()
    my_process_is_mpi_prefroot = (my_process_is_pref() .AND. &
     &                   (my_process_mpi_all_id == process_mpi_all_prefroot_id))
  END FUNCTION my_process_is_mpi_prefroot
  !------------------------------------------------------------------------------

  !------------------------------------------------------------------------------
  LOGICAL FUNCTION my_process_is_mpi_test()
    my_process_is_mpi_test = (my_mpi_function == test_mpi_process)
  END FUNCTION my_process_is_mpi_test
  !------------------------------------------------------------------------------

  !------------------------------------------------------------------------------
  !>
  ! If is mpi parallel and not a test process
  !! Note: mpi i/o processes do not count in mpi parallel work
  !!        Only computational processes are checked for running in parallel
  LOGICAL FUNCTION my_process_is_mpi_parallel()
    my_process_is_mpi_parallel = process_is_mpi_parallel
  END FUNCTION my_process_is_mpi_parallel
  !------------------------------------------------------------------------------


  !------------------------------------------------------------------------------
  !>
  ! If is mpi all parallel
  LOGICAL FUNCTION my_process_is_mpi_all_parallel()
    my_process_is_mpi_all_parallel = (process_mpi_all_size > 1)
  END FUNCTION my_process_is_mpi_all_parallel
  !------------------------------------------------------------------------------

  !------------------------------------------------------------------------------
  !>
  ! If is mpi all parallel
  LOGICAL FUNCTION my_process_is_mpi_all_seq()
    my_process_is_mpi_all_seq = (process_mpi_all_size <= 1)
  END FUNCTION my_process_is_mpi_all_seq
  !------------------------------------------------------------------------------

  !------------------------------------------------------------------------------
  !>
  !! If is not mpi work parallel or this is a test process
  !! returns TRUE
  !! Note: mpi i/o processes do not count is mpi parallel work
  !!        Only computational processes are checked for running in parallel
  LOGICAL FUNCTION my_process_is_mpi_seq()
    my_process_is_mpi_seq = .NOT. process_is_mpi_parallel
  END FUNCTION my_process_is_mpi_seq
  !------------------------------------------------------------------------------

  !------------------------------------------------------------------------------
  !>
  LOGICAL FUNCTION my_process_is_mpi_workroot()
!     my_process_is_mpi_workroot = (my_process_mpi_all_id == process_mpi_all_workroot_id)
    my_process_is_mpi_workroot = (p_pe_work == process_mpi_root_id)
  END FUNCTION my_process_is_mpi_workroot
  !------------------------------------------------------------------------------


  !------------------------------------------------------------------------------
  !>
  ! If is not mpi parallel or is a test process
  LOGICAL FUNCTION run_is_global_mpi_parallel()
    run_is_global_mpi_parallel = is_global_mpi_parallel
  END FUNCTION run_is_global_mpi_parallel
  !------------------------------------------------------------------------------

  !------------------------------------------------------------------------------
  SUBROUTINE print_info_stderr (name, text)
    CHARACTER (len=*), INTENT(in) :: name
    CHARACTER (len=*), INTENT(in) :: text

    IF (my_process_is_stdio() ) THEN
      WRITE (nerr,'(a,a,a,a)') " ", TRIM(name), ": ", TRIM(text)
    ENDIF

  END SUBROUTINE print_info_stderr
  !------------------------------------------------------------------------------


  !------------------------------------------------------------------------------
  SUBROUTINE finish (name, text)
    CHARACTER (len=*), INTENT(in) :: name
    CHARACTER (len=*), INTENT(in) :: text

    WRITE (nerr,'(i4.4,a,a,a,a)') my_global_mpi_id, ": ", TRIM(name), ": ", TRIM(text)
    CALL abort_mpi

  END SUBROUTINE finish
  !------------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  !
  !> Pushes the communicator and proc0 onto the communicator stack.
  !! The communicator stack is needed for global sums if the processor
  !! set is split among different 1st level patches.

  SUBROUTINE push_glob_comm(comm, proc0)
    INTEGER, INTENT(IN) :: comm, proc0

    ! Safety check
    IF(comm_lev>=max_lev) &
      CALL finish('push_glob_comm','max_lev exceeded')

    comm_lev = comm_lev+1
    glob_comm(comm_lev) = comm
    comm_proc0(comm_lev) = proc0

  END SUBROUTINE push_glob_comm

  !-------------------------------------------------------------------------
  !
  !> Pops one level of the communicator stack

  SUBROUTINE pop_glob_comm()

    ! Safety check
    IF(comm_lev<=0) &
      CALL finish('pop_glob_comm','stack empty')

    comm_lev = comm_lev-1

  END SUBROUTINE pop_glob_comm


  !-------------------------------------------------------------------------
  !
  !> @return current top of communicator stack

  FUNCTION get_glob_proc0()
    INTEGER :: get_glob_proc0
    get_glob_proc0 = comm_proc0(comm_lev)
  END FUNCTION get_glob_proc0


  !------------------------------------------------------------------------------
  !>
  !! Sets the p_comm_input_bcast
  !! If comm_flag == null_comm_type then
  !!    only the test process reads and
  !!    no broadcast takes place
  !! Otherwise
  !!    the test or the root process reads
  !!    and broadcasts to the rest
  SUBROUTINE set_comm_input_bcast (comm_flag)
    INTEGER, INTENT(in), OPTIONAL:: comm_flag

    INTEGER :: comm_type

    comm_type = default_comm_type

    IF (PRESENT(comm_flag)) THEN
       comm_type = comm_flag
    ENDIF

#ifndef NOMPI
    SELECT CASE(comm_type)

    CASE(null_comm_type)
      IF(my_process_is_mpi_test()) THEN
        p_comm_input_bcast = MPI_COMM_SELF ! i.e. effectively no broadcast
      ELSE
        p_comm_input_bcast = MPI_COMM_NULL ! Must not be used!
      ENDIF

    CASE default

      IF (my_process_is_io() .OR.(my_process_is_restart().OR. my_process_is_pref())) THEN
        ! I/O PEs and Restart PEs and prefetching PEs never participate in reading
        p_comm_input_bcast = MPI_COMM_NULL
      ELSE
        IF(is_mpi_test_run) THEN
          ! Test PE reads and broadcasts to workers
          p_comm_input_bcast = p_comm_work_test
        ELSE
          ! PE 0 reads and broadcasts
          p_comm_input_bcast = p_comm_work
        ENDIF
      ENDIF

    END SELECT
#endif

  END SUBROUTINE set_comm_input_bcast
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  SUBROUTINE set_mpi_work_communicators(p_test_run, l_test_openmp, num_io_procs, &
    &                                   num_restart_procs, num_prefetch_proc)
    LOGICAL,INTENT(INOUT) :: p_test_run, l_test_openmp
    INTEGER,INTENT(INOUT) :: num_io_procs
    INTEGER,INTENT(INOUT) :: num_restart_procs
    INTEGER,INTENT(INOUT), OPTIONAL :: num_prefetch_proc

!   !local variables
    INTEGER :: my_color, peer_comm, peer_comm_restart, peer_comm_pref, p_error
    CHARACTER(*), PARAMETER :: method_name = "set_mpi_work_communicators"
    INTEGER :: grp_process_mpi_all_comm, grp_comm_work_io, input_ranks(1), &
               translated_ranks(1)
    INTEGER :: sizeof_prefetch_processes

    IF (PRESENT(num_prefetch_proc)) THEN
      sizeof_prefetch_processes = num_prefetch_proc
    ELSE
      sizeof_prefetch_processes = 0
    ENDIF


! check l_test_openmp
#ifndef _OPENMP
    IF (l_test_openmp) THEN
      CALL print_info_stderr(method_name, &
        & 'l_test_openmp has no effect if the model is compiled without OpenMP support')
      CALL print_info_stderr(method_name, &
        & '--> l_test_openmp set to .FALSE.')
      l_test_openmp = .FALSE.
    END IF
#endif

    ! check p_test_run, num_io_procs and num_restart_procs
#ifdef NOMPI
    ! Unconditionally set p_test_run to .FALSE. and num_io_procs to 0,
    ! all other variables are already set correctly
    IF (p_test_run) THEN
      CALL print_info_stderr(method_name, &
       & 'p_test_run has no effect if the model is compiled with the NOMPI compiler directive')
      CALL print_info_stderr(method_name, &
       & '--> p_test_run set to .FALSE.')
      p_test_run = .FALSE.
    END IF
    IF (num_io_procs /= 0) THEN
      CALL print_info_stderr(method_name, &
       & 'num_io_procs has no effect if the model is compiled with the NOMPI compiler directive')
      CALL print_info_stderr(method_name, &
       & '--> num_io_procs set to 0')
      num_io_procs = 0
    END IF
    IF (num_restart_procs /= 0) THEN
      CALL print_info_stderr(method_name, &
       & 'num_restart_procs has no effect if the model is compiled with the NOMPI compiler directive')
      CALL print_info_stderr(method_name, &
       & '--> num_restart_procs set to 0')
      num_restart_procs = 0
    END IF
    IF (sizeof_prefetch_processes /= 0) THEN
      CALL print_info_stderr(method_name, &
      & 'sizeof_prefetch_processes has no effect if the model is compiled with the NOMPI compiler directive')
      CALL print_info_stderr(method_name, &
      & '--> sizeof_prefetch_processes set to 0')
      sizeof_prefetch_processes = 0
    END IF


    ! set the sequential values
    p_work_pe0 = 0
    num_io_procs = 0
    sizeof_prefetch_processes = 0
    num_work_procs = 1

#else

    ! A run on 1 PE is never a verification run,
    ! correct this if the user should set it differently
    IF (process_mpi_all_size < 2) THEN
      IF (p_test_run) THEN
        CALL print_info_stderr(method_name, &
            & 'p_test_run has no effect in seq run')
        CALL print_info_stderr(method_name, &
            & '--> p_test_run set to .FALSE.')
        p_test_run = .FALSE.
      ENDIF
      IF (num_io_procs > 0) THEN
        CALL print_info_stderr(method_name, &
            & 'num_io_procs cannot be > 0 in seq run')
        CALL print_info_stderr(method_name, &
            & '--> num_io_procs set to 0')
        num_io_procs = 0
      ENDIF
      IF (num_restart_procs > 0) THEN
        CALL print_info_stderr(method_name, &
            & 'num_restart_procs cannot be > 0 in seq run')
        CALL print_info_stderr(method_name, &
            & '--> num_restart_procs set to 0')
        num_restart_procs = 0
      ENDIF
      IF (sizeof_prefetch_processes > 0) THEN
        CALL print_info_stderr(method_name, &
            & 'sizeof_prefetch_processes cannot be > 0 in seq run')
        CALL print_info_stderr(method_name, &
            & '--> sizeof_prefetch_processes set to 0')
        sizeof_prefetch_processes = 0
      ENDIF
    ENDIF
    IF(num_io_procs < 0) num_io_procs = 0           ! for safety only
    IF(num_restart_procs < 0) num_restart_procs = 0 ! for safety only
    IF(sizeof_prefetch_processes < 0) sizeof_prefetch_processes = 0 ! for safety only
    ! -----------------------------------------
    ! Set if test
    IF(p_test_run) THEN
      num_test_procs = 1
      p_test_pe = 0
    ELSE
      num_test_procs = 0
      p_test_pe = -1
    ENDIF

    ! -----------------------------------------
    ! how many work processors?
    num_work_procs = process_mpi_all_size - num_test_procs - num_io_procs - num_restart_procs - sizeof_prefetch_processes

    ! Check if there are sufficient PEs at all
    IF(num_work_procs < 1) THEN
      CALL finish(method_name, &
      & 'not enough processors for given values of p_test_run/num_io_procs/num_restart_procs/sizeof_prefetch_processes')
    ELSE IF (p_test_run .AND. num_work_procs == 1) THEN
#ifdef _OPENACC
      PRINT *, "Testing 1 CPU against 1 GPU"
#else
      CALL finish(method_name, &
      & 'running p_test_run with only 1 work processor does not make sense')
#endif
    ENDIF

    WRITE(message_text,'(5(a,i0))') 'Number of procs for test: ',num_test_procs, &
      & ', work: ',num_work_procs, &
      & ', I/O: ',num_io_procs, &
      & ', Restart: ',num_restart_procs, &
      & ', Prefetching: ',sizeof_prefetch_processes
    CALL print_info_stderr(method_name, message_text)

    ! Everything seems ok. Proceed to setup the communicators and ids
    ! Set up p_test_pe, p_work_pe0, p_io_pe0, p_restart_pe0, p_pref_pe0
    ! which are identical on all PEs
    p_work_pe0    = num_test_procs
    p_io_pe0      = num_test_procs + num_work_procs
    p_restart_pe0 = num_test_procs + num_work_procs + num_io_procs
    p_pref_pe0    = num_test_procs + num_work_procs + num_io_procs + num_restart_procs

    ! Set up p_n_work and p_pe_work which are NOT identical on all PEs
    IF(p_pe < p_work_pe0) THEN
      ! Test PE (if present)
      p_n_work  = 1          ! 1 PE in verification work group
      p_pe_work = 0          ! PE number within work group
    ELSE IF(p_pe < p_io_pe0) THEN
      ! Work PE
      p_n_work  = num_work_procs
      p_pe_work = p_pe - num_test_procs
    ELSE IF(p_pe < p_restart_pe0) THEN
      ! I/O PE (if present)
      p_n_work  = num_io_procs
      p_pe_work = p_pe - num_test_procs - num_work_procs
    ELSE IF(p_pe < p_pref_pe0) THEN
      ! Restart PE (if present)
      p_n_work  = num_restart_procs
      p_pe_work = p_pe - num_test_procs - num_work_procs - num_io_procs
    ELSE
      p_n_work  = sizeof_prefetch_processes
      p_pe_work = p_pe - num_test_procs - num_work_procs - num_io_procs - num_restart_procs
    ENDIF


    ! Set communicators
    ! =================

    ! Split communicator process_mpi_all_comm between test/work/io/restart/prefetching
    ! to get  which is the communicator for
    ! usage WITHIN every group of the 5 different types
    IF(p_pe < p_work_pe0) THEN
      my_mpi_function = test_mpi_process
    ELSE IF(p_pe < p_io_pe0) THEN
      my_mpi_function = work_mpi_process
    ELSE IF(p_pe < p_restart_pe0) THEN
      my_mpi_function = io_mpi_process
    ELSE IF(p_pe < p_pref_pe0) THEN
      my_mpi_function = restart_mpi_process
    ELSE
      my_mpi_function = pref_mpi_process
    ENDIF

    CALL MPI_Comm_split(process_mpi_all_comm, my_mpi_function, p_pe, p_comm_work, p_error)

    ! Set p_comm_work_test, the communicator spanning work group and test PE
    IF(p_test_run) THEN
      IF(p_pe < p_io_pe0) THEN
        my_color = 1
      ELSE
        my_color = MPI_UNDEFINED ! p_comm_work_test must never be used on I/O PEs
      ENDIF

      CALL MPI_Comm_split(process_mpi_all_comm, my_color, p_pe, p_comm_work_test, p_error)
    ELSE
      ! If not a test run, p_comm_work_test must not be used at all
      p_comm_work_test = MPI_COMM_NULL
    ENDIF

    ! Set p_comm_work_io, the communicator spanning work group and I/O PEs
    IF (p_test_run .AND. (p_pe < p_work_pe0)) THEN
       my_color = 1
    ELSE IF ((num_io_procs > 0) .AND. (p_pe >= p_work_pe0) .AND. (p_pe < p_restart_pe0)) THEN
       my_color = 2
    ELSE IF ((num_io_procs == 0) .AND. (p_pe < p_restart_pe0)) THEN
       my_color = 3
    ELSE
       my_color = MPI_UNDEFINED
    END IF
    CALL MPI_Comm_split(process_mpi_all_comm, my_color, p_pe, p_comm_work_io, p_error)
    if (my_color == MPI_UNDEFINED)  p_comm_work_io = MPI_COMM_NULL

    ! Set p_comm_io, the communicator spanning the I/O PEs
    IF (p_test_run .AND. (p_pe < p_work_pe0)) THEN
       my_color = 1
    ELSE IF (((num_io_procs > 0) .AND. ((p_pe >= p_io_pe0) .AND. (p_pe < p_restart_pe0))) .OR.  &
         &   ((num_io_procs == 0) .AND. (p_pe == p_work_pe0))) THEN
       my_color = 2
    ELSE
       my_color = MPI_UNDEFINED ! p_comm_work_io must only be used on test, work, I/O PEs
    END IF
    CALL MPI_Comm_split(process_mpi_all_comm, my_color, p_pe, p_comm_io, p_error)
    if (my_color == MPI_UNDEFINED)  p_comm_io = MPI_COMM_NULL

    ! translate the rank "p_io_pe0" to the communicator "p_comm_work_io":
    IF (p_comm_work_io /= MPI_COMM_NULL) THEN
      CALL MPI_Comm_group(process_mpi_all_comm, grp_process_mpi_all_comm, p_error)
      CALL MPI_Comm_group(p_comm_work_io,       grp_comm_work_io,         p_error)
      IF (num_io_procs > 0) THEN
        input_ranks(1) = p_io_pe0
      ELSE IF (p_test_run .AND. (p_pe < p_work_pe0)) THEN
        input_ranks(1) = 0
      ELSE
        input_ranks(1) = p_work_pe0
      END IF
      CALL MPI_group_translate_ranks(grp_process_mpi_all_comm, 1, input_ranks, &
        &                            grp_comm_work_io, translated_ranks, p_error)
      process_work_io0 = translated_ranks(1)
      CALL MPI_group_free(grp_process_mpi_all_comm, p_error)
      CALL MPI_group_free(grp_comm_work_io, p_error)
    END IF

    ! Set p_comm_work_restart, the communicator spanning work group and Restart Ouput PEs
    IF(num_restart_procs > 0) THEN
      IF(p_pe < p_work_pe0 .OR. &
      &  ((num_io_procs > 0 .AND. (p_pe >= p_io_pe0 .AND. p_pe < p_restart_pe0)) .OR. &
      &  (sizeof_prefetch_processes > 0 .AND. p_pe >= p_pref_pe0))) THEN
        my_color = MPI_UNDEFINED ! p_comm_work_restart must never be used on test and IO PE
      ELSE
        my_color = 1    ! This is set only for all workers and for all restart PEs
      ENDIF

      CALL MPI_Comm_split(process_mpi_all_comm, my_color, p_pe, p_comm_work_restart, p_error)
    ELSE
      ! If no I/O PEs are present, p_comm_work_io must not be used at all
      p_comm_work_restart = MPI_COMM_NULL
    ENDIF

    ! Set p_comm_work_pref, the communicator spanning work group and prefetching PEs
    IF(sizeof_prefetch_processes > 0) THEN
      IF(p_pe < p_work_pe0 .OR. &
      &  ((num_io_procs > 0 .AND. (p_pe >= p_io_pe0 .AND. p_pe < p_restart_pe0)) .OR. &
      &  ( num_restart_procs > 0 .AND. (p_pe >= p_restart_pe0 .AND. p_pe < p_pref_pe0)))) THEN
        my_color = MPI_UNDEFINED ! p_comm_work_restart must never be used on test and IO PE
      ELSE
        my_color = 1    ! This is set only for all workers and for all prefetching PEs
      ENDIF

      CALL MPI_Comm_split(process_mpi_all_comm, my_color, p_pe, p_comm_work_pref, p_error)
    ELSE
      ! If no prefetching PEs are present, p_comm_inp_pref must not be used at all
      p_comm_work_pref = MPI_COMM_NULL
    ENDIF

!     The following is moved to set_comm_input_bcast
!     ! Set p_comm_input_bcast, the communicator for broadcasting the NetCDF input
!     IF(lrestore_states) THEN
!       ! NetCDF input is only read by the test pe and MUST NOT be broadcast
!       IF(p_pe == p_test_pe) THEN
!         p_comm_input_bcast = MPI_COMM_SELF ! i.e. effectively no broadcast
!       ELSE
!         p_comm_input_bcast = MPI_COMM_NULL ! Must not be used!
!       ENDIF
!     ELSE
!       IF(p_pe < p_io_pe0) THEN
!         IF(p_test_run) THEN
!           ! Test PE reads and broadcasts to workers
!           p_comm_input_bcast = p_comm_work_test
!         ELSE
!           ! PE 0 reads and broadcasts
!           p_comm_input_bcast = p_comm_work
!         ENDIF
!       ELSE
!         ! I/O PEs never participate in reading
!         p_comm_input_bcast = MPI_COMM_NULL
!       ENDIF
!     ENDIF

    ! Create Intercommunicator work PEs - I/O PEs

    ! From MPI-Report:
    ! Advice to users: We recommend using a dedicated peer communicator, such as a
    ! duplicate of MPI_COMM_WORLD, to avoid trouble with peer communicators.

    ! No idea what these troubles may be in reality, but let us follow the advice

    CALL MPI_Comm_dup(process_mpi_all_comm, peer_comm, p_error)

    IF(p_pe /= p_test_pe .AND. num_io_procs>0) THEN

      IF(p_pe < p_io_pe0) THEN
        CALL MPI_Intercomm_create(p_comm_work, 0, peer_comm, p_io_pe0, &
          & 1, p_comm_work_2_io, p_error)
      ELSE IF(p_pe < p_restart_pe0) THEN
        CALL MPI_Intercomm_create(p_comm_work, 0, peer_comm, p_work_pe0,&
          & 1, p_comm_work_2_io, p_error)
      ELSE
        ! No Intercommunicator for restart PEs
        p_comm_work_2_io = MPI_COMM_NULL
      ENDIF
    ELSE
      ! No Intercommunicator for test PE or for all, when no IO PEs are defined
      p_comm_work_2_io = MPI_COMM_NULL
    ENDIF

    ! Perform the same as above, but create the intra-communicators between
    ! the worker PEs and the restart PEs.
    CALL MPI_Comm_dup(process_mpi_all_comm, peer_comm_restart, p_error)

    IF(p_pe /= p_test_pe .AND. num_restart_procs>0) THEN

      IF(p_pe < p_io_pe0) THEN ! All workers
        CALL MPI_Intercomm_create(p_comm_work, 0, peer_comm_restart, p_restart_pe0, &
          & 2, p_comm_work_2_restart, p_error)
      ELSE IF((p_pe >= p_restart_pe0).AND. (p_pref_pe0 > p_pe)) THEN ! All restart PEs
        CALL MPI_Intercomm_create(p_comm_work, 0, peer_comm_restart, p_work_pe0,&
          & 2, p_comm_work_2_restart, p_error)
      ELSE
        ! No Intercommunicator for IO PEs
        p_comm_work_2_restart = MPI_COMM_NULL
      ENDIF
    ELSE
      ! No Intercommunicator for Test PE or for all, when no restart PEs
      p_comm_work_2_restart = MPI_COMM_NULL
    ENDIF

    ! Perform the same as above, but create the intra-communicators between
    ! the worker PEs and the prefetching PEs.
    CALL MPI_Comm_dup(process_mpi_all_comm, peer_comm_pref, p_error)

    IF(p_pe /= p_test_pe .AND. sizeof_prefetch_processes>0) THEN

      IF(p_pe < p_io_pe0) THEN ! All workers
        CALL MPI_Intercomm_create(p_comm_work, 0, peer_comm_pref, p_pref_pe0, &
          & 3, p_comm_work_2_pref, p_error)
      ELSE IF(p_pe >= p_pref_pe0) THEN ! All prefetching PEs
        CALL MPI_Intercomm_create(p_comm_work, 0, peer_comm_pref, p_work_pe0,&
          & 3, p_comm_work_2_pref, p_error)
      ELSE
        ! No Intercommunicator for prefetching PEs
        p_comm_work_2_pref = MPI_COMM_NULL
      ENDIF
    ELSE
      ! No Intercommunicator for Test PE or for all, when no restart PEs
      p_comm_work_2_pref = MPI_COMM_NULL
    ENDIF

    ! if OpenMP is used, the test PE uses only 1 thread in order to check
    ! the correctness of the OpenMP implementation
    ! Currently the I/O PEs are also single threaded!
#ifdef _OPENMP
    IF (l_test_openmp .AND. p_pe == p_test_pe) CALL OMP_SET_NUM_THREADS(1)
    IF (p_pe >= p_io_pe0) CALL OMP_SET_NUM_THREADS(1)
#endif

#endif

    ! fill some derived variables
    process_mpi_all_workroot_id     = p_work_pe0
    process_mpi_all_ioroot_id       = p_io_pe0
    process_mpi_all_test_id         = p_work_pe0-1
    process_mpi_all_restartroot_id  = p_restart_pe0
    process_mpi_restart_size        = num_restart_procs
    process_mpi_io_size             = num_io_procs
    process_mpi_all_prefroot_id     = p_pref_pe0
    process_mpi_pref_size           = sizeof_prefetch_processes
    p_comm_work_pref_compute_pe0    = p_work_pe0

    ! In case of test run, only the test process is stdio
    process_is_stdio = (my_process_mpi_all_id == process_mpi_stdio_id)
    process_is_mpi_parallel = (num_work_procs > 1)
    IF (my_process_is_mpi_test()) process_is_mpi_parallel = .false.

    ! still to be filled
!     process_mpi_local_comm  = process_mpi_all_comm
!     process_mpi_local_size  = process_mpi_all_size
!     my_process_mpi_local_id = my_process_mpi_all_id

#ifdef DEBUG
      WRITE (nerr,'(a,a,i5)') method_name, ' p_pe=',            p_pe
      WRITE (nerr,'(a,a,i5)') method_name, ' p_work_pe0=',      p_work_pe0
      WRITE (nerr,'(a,a,i5)') method_name, ' p_io_pe0=',        p_io_pe0
      WRITE (nerr,'(a,a,i5)') method_name, ' p_restart_pe0=',   p_restart_pe0
      WRITE (nerr,'(a,a,i5)') method_name, ' p_comm_work=',     p_comm_work
      WRITE (nerr,'(a,a,i5)') method_name, ' p_comm_work_2_restart=', &
        & p_comm_work_2_restart
      WRITE (nerr,'(a,a,i5)') method_name, ' my_mpi_function=', my_mpi_function
#endif

    ! fill my  parameters
    is_mpi_test_run = p_test_run
    is_openmp_test_run = l_test_openmp
    ! fill other default
    CALL set_comm_input_bcast()

    IF (PRESENT(num_prefetch_proc)) THEN
      num_prefetch_proc = sizeof_prefetch_processes
    ENDIF

  END SUBROUTINE set_mpi_work_communicators
  !-------------------------------------------------------------------------

  !------------------------------------------------------------------------------
  !>
  SUBROUTINE set_default_mpi_work_variables()

    ! fill some derived variables
    process_is_stdio = (my_process_mpi_all_id == process_mpi_stdio_id)
    process_is_mpi_parallel = (process_mpi_all_size > 1)
    my_mpi_function = work_mpi_process

    process_mpi_all_test_id     = -1
    process_mpi_all_workroot_id = 0
    process_mpi_io_size         = 0
    process_mpi_restart_size    = 0
    process_mpi_pref_size       = 0
    p_comm_work_pref_compute_pe0 = 0
    is_mpi_test_run = .false.
    is_openmp_test_run = .false.

!     process_mpi_local_comm  = process_mpi_all_comm
!     process_mpi_local_size  = process_mpi_all_size
!     my_process_mpi_local_id = my_process_mpi_all_id

    ! set some of the old variables
    ! should be removed once the old variables are cleaned
    p_pe           = my_process_mpi_all_id
    p_io           = 0
    num_test_procs = 0
    num_work_procs = process_mpi_all_size
    p_test_pe      = -1
    p_work_pe0     = 0
    p_io_pe0       = process_mpi_all_size    ! Number of I/O PE 0 within all PEs (process_mpi_all_size if no I/O PEs)
    ! Number of restart PE 0 within all PEs (process_mpi_all_size if no restart PEs)
    p_restart_pe0  = process_mpi_all_size
    ! Number of prefetching PE 0 within all PEs (process_mpi_all_size if no prefetching PEs)
    p_pref_pe0     = process_mpi_all_size
    p_n_work       = process_mpi_all_size
    p_pe_work      = my_process_mpi_all_id

    p_comm_work             = process_mpi_all_comm
    p_comm_input_bcast      = process_mpi_all_comm
    p_comm_work_io          = MPI_COMM_NULL
    p_comm_work_test        = MPI_COMM_NULL
    p_comm_work_2_io        = MPI_COMM_NULL
    p_comm_work_restart     = MPI_COMM_NULL
    p_comm_work_2_restart   = MPI_COMM_NULL
    p_comm_io               = MPI_COMM_NULL
    P_comm_work_pref        = MPI_COMM_NULL
    P_comm_work_2_pref      = MPI_COMM_NULL

    ! print some info
    IF ( .NOT. process_is_mpi_parallel) THEN
      WRITE (nerr,'(/,a,/)')  ' Single processor run.'
    ELSEIF (process_is_stdio) THEN
      WRITE (nerr,'(/,a,a,a,i0,a,/)') ' ', TRIM(process_mpi_name), &
        ' runs on ', process_mpi_all_size, ' mpi processes.'
    END IF

  END SUBROUTINE set_default_mpi_work_variables
  !------------------------------------------------------------------------------

  !------------------------------------------------------------------------------
  !>
  !! Splits the global communicator into this component's communicator
  !! Should be called before the component configuration
  SUBROUTINE split_global_mpi_communicator(component_no)
    INTEGER, INTENT(in) :: component_no

    INTEGER :: new_communicator
    LOGICAL             :: l_mpi_is_initialised
    CHARACTER(len=*), PARAMETER :: method_name = 'split_process_mpi_communicator'
#ifdef NOMPI
    RETURN
#else
    !--------------------------------------------
    ! check if mpi is initialized
    CALL MPI_INITIALIZED(l_mpi_is_initialised, p_error)
    IF (p_error /= MPI_SUCCESS) THEN
      WRITE (nerr,'(a,a)') method_name, ' MPI_INITITIALIZED failed.'
      WRITE (nerr,'(a,i4)') ' Error =  ', p_error
      STOP
    END IF
    !--------------------------------------------
    ! split global_mpi_communicator
    CALL MPI_Comm_split(global_mpi_communicator, component_no, my_global_mpi_id, &
      & new_communicator, p_error)
    IF (p_error /= MPI_SUCCESS) THEN
      WRITE (nerr,'(a,a)') method_name, ' MPI_Comm_split failed.'
      WRITE (nerr,'(a,i4)') ' Error =  ', p_error
      STOP
    END IF
    CALL set_process_mpi_communicator(new_communicator)
#endif

  END SUBROUTINE split_global_mpi_communicator
  !------------------------------------------------------------------------------

  !------------------------------------------------------------------------------
  !>
  !! Set this component's communicator
  !! The new communicator id is duplicated from the input communicator
  !! Should be called before the component configuration
  SUBROUTINE set_process_mpi_communicator(new_communicator)
    INTEGER, INTENT(in) :: new_communicator

    LOGICAL             :: l_mpi_is_initialised
    CHARACTER(len=*), PARAMETER :: method_name = 'set_process_mpi_communicator'

#ifdef NOMPI
    process_mpi_all_comm    = new_communicator
    process_mpi_all_size    = 1
    my_process_mpi_all_id   = 0

#else
    ! since here we define the process all communicator
    ! the work communicator is identical to the all communicator
    ! and the no test or i/o processes are present
    CALL MPI_INITIALIZED(l_mpi_is_initialised, p_error)

    IF (p_error /= MPI_SUCCESS) THEN
      WRITE (nerr,'(a,a)') method_name, ' MPI_INITITIALIZED failed.'
      WRITE (nerr,'(a,i4)') ' Error =  ', p_error
      STOP
    END IF

    IF ( .NOT. l_mpi_is_initialised ) THEN
       WRITE (nerr,'(a,a)') method_name, &
         & ' MPI_Init or start_mpi needs to be called first.'
       STOP
    ENDIF

    IF ( process_mpi_all_comm /= MPI_COMM_NULL) THEN
      ! free original communicator
      CALL MPI_COMM_FREE(process_mpi_all_comm, p_error)
      IF (p_error /= MPI_SUCCESS) THEN
        WRITE (nerr,'(a,a)') method_name, &
          & ' MPI_COMM_FREE failed. start_mpi needs to be called before.'
        WRITE (nerr,'(a,i4)') ' Error =  ', p_error
        CALL abort_mpi
      END IF
    ENDIF
    ! assign MPI communicator generated elsewhere to process_mpi_all_comm

    CALL MPI_COMM_DUP(new_communicator, process_mpi_all_comm, p_error)
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,a)') method_name,&
         & ' MPI_COMM_DUP failed for process_mpi_all_comm.'
       WRITE (nerr,'(a,i4)') ' Error =  ', p_error
       CALL abort_mpi
    END IF

    ! get local PE identification
    CALL MPI_COMM_RANK (process_mpi_all_comm, my_process_mpi_all_id, p_error)

    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,a)') method_name, ' MPI_COMM_RANK failed.'
       WRITE (nerr,'(a,i4)') ' Error =  ', p_error
       CALL abort_mpi
    ELSE
#ifdef DEBUG
       WRITE (nerr,'(a,a,i4,a)') method_name, ' my_process_mpi_all_id ', &
         & my_process_mpi_all_id, ' started.'
#endif
    END IF


    CALL MPI_COMM_SIZE (process_mpi_all_comm, process_mpi_all_size, p_error)

    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,a,i4,a)') method_name, ' PE: ',&
         & my_process_mpi_all_id, ' MPI_COMM_SIZE failed.'
       WRITE (nerr,'(a,i4)') ' Error =  ', p_error
       CALL abort_mpi
    END IF


!     CALL MPI_COMM_DUP(process_mpi_all_comm,p_comm_work,p_error)
!     IF (p_error /= MPI_SUCCESS) THEN
!        WRITE (nerr,'(a,a)') method_name, ' MPI_COMM_DUP failed for p_comm_work.'
!        WRITE (nerr,'(a,i4)') ' Error =  ', p_error
!        CALL abort_mpi
!     END IF
!
!     CALL MPI_COMM_DUP(process_mpi_all_comm,p_comm_work_test,p_error)
!     IF (p_error /= MPI_SUCCESS) THEN
!        WRITE (nerr,'(a,a)') method_name, ' MPI_COMM_DUP failed for p_comm_work_test.'
!        WRITE (nerr,'(a,i4)') ' Error =  ', p_error
!        CALL abort_mpi
!     END IF

#endif

    CALL set_default_mpi_work_variables()

  END SUBROUTINE set_process_mpi_communicator
  !------------------------------------------------------------------------------

  !------------------------------------------------------------------------------
  SUBROUTINE start_mpi(global_name)

!#ifdef _OPENMP
!    USE mo_util_string, ONLY: toupper
!#endif

#ifndef NOMPI
#if defined (__prism) && defined (use_comm_MPI1)
    USE mod_prism_proto, ONLY: prism_ok
#endif
#endif

    CHARACTER(len=*), INTENT(in), OPTIONAL :: global_name

    ! variables are required for determing I/O size in bytes of the defined
    ! KIND types for assigning the right MPI data types with the used kinds
    INTEGER     :: iig = 0
    INTEGER(i4) :: ii4 = 0_i4
    INTEGER(i8) :: ii8 = 0_i8
    REAL        :: rrg = 0.0
    REAL(sp)    :: rsp = 0.0_sp
    REAL(dp)    :: rdp = 0.0_dp

    CHARACTER(len=132) :: yname

    ! variables used for determing the OpenMP threads
    ! suitable as well for coupled models
#if (defined _OPENMP)
!    CHARACTER(len=32) :: env_name
!    CHARACTER(len=32) :: thread_num
!    INTEGER :: threads
    INTEGER :: global_no_of_threads
#ifndef NOMPI
    INTEGER :: provided
#endif
#ifndef __SX__
    ! status
!    INTEGER :: istat
#else
    EXTERNAL :: getenv
#endif
#endif

#ifndef NOMPI
#if defined (__prism) && defined (use_comm_MPI1)
    INTEGER :: prism_model_number
    CHARACTER(len=132) :: prism_model_name

    EXTERNAL :: prism_abort_proto
    EXTERNAL :: prism_init_comp_proto
    EXTERNAL :: prism_get_localcomm_proto
#endif
#endif
    CHARACTER(len=*), PARAMETER :: method_name = 'start_mpi'


! #ifdef _OPENMP
!     global_no_of_threads = 1
! #endif


    ! start MPI
#ifndef NOMPI
#ifdef _OPENMP
#ifdef __MULTIPLE_MPI_THREADS
    CALL MPI_INIT_THREAD(MPI_THREAD_MULTIPLE,provided,p_error)
#else
    CALL MPI_INIT_THREAD(MPI_THREAD_FUNNELED,provided,p_error)
#endif
#elif (defined COUP_OAS_ICON)
    ! init OASIS coupler
    CALL oasis_init_comp(oas_comp_id, oas_comp_name, oas_error)
    IF (oas_error /= 0) &
      CALL oasis_abort(oas_comp_id, 'oas_icon_init', 'Failure in oasis_init')
    CALL oasis_get_localcomm(kl_comm, oas_error)
    IF (oas_error /= 0) &
      CALL oasis_abort(oas_comp_id, 'oas_icon_init', 'Failure in oasis_get_localcomm')
#else
    CALL MPI_INIT (p_error)
#endif
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,a)') method_name, ' MPI_INIT failed.'
       WRITE (nerr,'(a,i4)') ' Error =  ', p_error
       STOP
    END IF

#ifdef _OPENMP
    ! Check if MPI_INIT_THREAD returned at least MPI_THREAD_FUNNELED in "provided"
#ifdef __MULTIPLE_MPI_THREADS
    IF (provided < MPI_THREAD_MULTIPLE) THEN
       WRITE (nerr,'(a,a)') method_name, &
         & ' MPI_INIT_THREAD did not return desired level of thread support'
       WRITE (nerr,'(a,i0)') " provided: ", provided
       WRITE (nerr,'(a,i0)') " required: ", MPI_THREAD_MULTIPLE
       CALL MPI_Finalize(p_error)
       STOP
    END IF
#else
    IF (provided < MPI_THREAD_FUNNELED) THEN
       WRITE (nerr,'(a,a)') method_name, &
         & ' MPI_INIT_THREAD did not return desired level of thread support'
       WRITE (nerr,'(a,i0)') " provided: ", provided
       WRITE (nerr,'(a,i0)') " required: ", MPI_THREAD_FUNNELED
       ! CALL MPI_Finalize(p_error)
       ! STOP
    END IF
#endif
#endif

    CALL MPI_BUFFER_ATTACH(mpi_buffer, SIZE(mpi_buffer), p_error)
    IF (p_error /= 0) THEN
       WRITE (0,*) "Error in MPI_BUFFER_ATTACH."
       STOP
    END IF

    process_mpi_all_comm = MPI_COMM_NULL
    IF (PRESENT(global_name)) THEN
      yname = TRIM(global_name)
    ELSE
      yname = '(unnamed)'
    END IF
#ifdef COUP_OAS_ICON
    yname = TRIM("icon_COUPLED")
#endif
    global_mpi_name = TRIM(yname)

    ! create communicator for this process alone before
    ! potentially joining MPI2
#if defined (__prism) && defined (use_comm_MPI1)

    prism_model_name = TRIM(yname)

    CALL prism_init_comp_proto (prism_model_number, TRIM(prism_model_name), &
         p_error)

    IF (p_error /= prism_ok) THEN
      WRITE (nerr,*) method_name, ' prism_init_comp_proto failed'
      CALL prism_abort_proto(prism_model_number, TRIM(yname),'abort1')
    ENDIF

    CALL prism_get_localcomm_proto(global_mpi_communicator, p_error)

    IF (p_error /= prism_ok) THEN
      WRITE (nerr,*) method_name, ' prism_get_localcomm_proto failed'
      CALL prism_abort_proto(prism_model_number, TRIM(yname),'abort2')
    ENDIF

#else

#ifdef COUP_OAS_ICON
    CALL MPI_COMM_DUP (kl_comm, global_mpi_communicator, p_error)
#else
    CALL MPI_COMM_DUP (MPI_COMM_WORLD, global_mpi_communicator, p_error)
#endif

    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,a)') method_name, ' MPI_COMM_DUP failed.'
       WRITE (nerr,'(a,i4)') ' Error =  ', p_error
       CALL abort_mpi
    END IF

#endif

    ! get local PE identification
    CALL MPI_COMM_RANK (global_mpi_communicator, my_global_mpi_id, p_error)

    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,a)') method_name, ' MPI_COMM_RANK failed.'
       WRITE (nerr,'(a,i4)') ' Error =  ', p_error
       CALL abort_mpi
    END IF

    ! Information ...
    IF (my_global_mpi_id == 0) THEN
       WRITE (nerr,'(/,a,a,a)') ' ', &
            TRIM(yname), ' MPI interface runtime information:'
    END IF

#ifdef DEBUG
    WRITE (nerr,'(a,a,i4,a)') method_name, ' PE ', my_global_mpi_id, ' started.'
#endif

    ! get number of available PEs
    CALL MPI_COMM_SIZE (global_mpi_communicator, global_mpi_size, p_error)

    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,a,i4,a)') method_name ,' PE: ', my_global_mpi_id, &
       &  ' MPI_COMM_SIZE failed.'
       WRITE (nerr,'(a,i4)') ' Error =  ', p_error
       CALL abort_mpi
    END IF

    ! for non blocking calls
    p_mrequest = p_request_alloc_size
    ALLOCATE(p_request(p_mrequest))
    p_irequest = 0

    ! lets check the available MPI version
    CALL MPI_GET_VERSION (version, subversion, p_error)

    IF (p_error /= MPI_SUCCESS) THEN
      WRITE (nerr,'(a,a)') method_name , ' MPI_GET_VERSION failed.'
      WRITE (nerr,'(a,i4)') ' Error =  ', p_error
      CALL abort_mpi
    END IF

    IF (my_global_mpi_id == 0) THEN
      WRITE (nerr,'(a,a,a,i0,a1,i0)') " ", method_name, &
           ' Used MPI version: ', version, '.', subversion
    END IF

    IF (my_global_mpi_id == 0) THEN
      WRITE (nerr,'(a,a,a,a,a,i0,a)') " ", method_name, " ", &
        & TRIM(yname), ': Globally run on ',&
        & global_mpi_size, ' mpi processes.'
    END IF

   ! due to a possible circular dependency with mo_machine and other
    ! modules, we determine here locally the I/O size of the different
    ! kind types (assume 8 bit/byte. This is than used for determing
    ! the right MPI send/receive type parameters.
    CALL MPI_SIZEOF(iig, p_int_byte, p_error)
    CALL MPI_SIZEOF(ii4, p_int_i4_byte, p_error)
    CALL MPI_SIZEOF(ii8, p_int_i8_byte, p_error)
    CALL MPI_SIZEOF(rrg, p_real_byte, p_error)
    CALL MPI_SIZEOF(rsp, p_real_sp_byte, p_error)
    CALL MPI_SIZEOF(rdp, p_real_dp_byte, p_error)

    p_int     = MPI_INTEGER
    p_real    = MPI_REAL
    p_bool    = MPI_LOGICAL
    p_char    = MPI_CHARACTER

    CALL MPI_TYPE_MATCH_SIZE(MPI_TYPECLASS_REAL, p_real_sp_byte, p_real_sp, p_error)
    CALL MPI_TYPE_MATCH_SIZE(MPI_TYPECLASS_REAL, p_real_dp_byte, p_real_dp, p_error)
    CALL MPI_TYPE_MATCH_SIZE(MPI_TYPECLASS_INTEGER, p_int_i4_byte, p_int_i4, p_error)
    CALL MPI_TYPE_MATCH_SIZE(MPI_TYPECLASS_INTEGER, p_int_i8_byte, p_int_i8, p_error)


#ifdef DEBUG
    WRITE (nerr,'(/,a)')    ' MPI transfer sizes [bytes]:'
    WRITE (nerr,'(a,i4)') '  INTEGER generic:', p_int_byte
    WRITE (nerr,'(a,i4)') '  INTEGER 4 byte :', p_int_i4_byte
    WRITE (nerr,'(a,i4)') '  INTEGER 8 byte :', p_int_i8_byte
    WRITE (nerr,'(a,i4)') '  REAL generic   :', p_real_byte
    WRITE (nerr,'(a,i4)') '  REAL single    :', p_real_sp_byte
    WRITE (nerr,'(a,i4)') '  REAL double    :', p_real_dp_byte
#endif

! MPI ends here
#else

    WRITE (nerr,'(a,a)')  method_name, ' No MPI: Single processor run.'
    ! set defaults for sequential run
    global_mpi_communicator = MPI_COMM_NULL
    global_mpi_size  = 1        ! total number of processes in global world
    my_global_mpi_id = 0        ! process id in global world
    is_global_mpi_parallel = .false.
    process_mpi_name = 'uknown'
    process_mpi_all_comm = MPI_COMM_NULL

#endif


    IF (global_mpi_size > 1) THEN
      is_global_mpi_parallel = .TRUE.
    ELSE
      global_mpi_size = 1 ! just in case we got wrong size
    ENDIF

#ifdef _OPENMP
    ! The number of threads, if varying, will be defined via
    ! namelists
    global_no_of_threads = omp_get_max_threads()
#ifndef NOMPI
    ! Make number of threads from environment available to all model PEs
!     CALL MPI_BCAST (global_no_of_threads, 1, MPI_INTEGER, 0, global_mpi_communicator, p_error)
#endif

    IF (my_global_mpi_id == 0) THEN

      IF (is_global_mpi_parallel) THEN
        WRITE (nerr,'(a,a,a)') " ", method_name, &
          & ': Running globally hybrid OpenMP-MPI mode.'
      ELSE
        WRITE (nerr,'(a,a,a)') " ", method_name,': Running globally OpenMP mode.'
      ENDIF
      WRITE (nerr,'(a,a, a, i0)') " ", method_name, &
        & ' global_no_of_threads is ', global_no_of_threads
    ENDIF

#endif


    ! The number of threads, if varying, will be defined via
    ! namelists
! #ifdef _OPENMP
!     ! Expect that PE 0 did got the information of OMP_NUM_THREADS.
!     ! That might be wrong in the coupled case when the model is
!     ! started via MPI dynamic process creation. So we have to check
!     ! the environment variable too.
!     IF (my_global_mpi_id == 0) THEN
!
!       IF (is_global_mpi_parallel) THEN
!         WRITE (nerr,'(/,a)') ' Running globally hybrid OpenMP-MPI mode.'
!       ELSE
!         WRITE (nerr,'(/,a)') ' Running globally OpenMP mode.'
!       ENDIF
!
!       env_name = toupper(TRIM(yname)) // '_THREADS'
! #ifdef __SX__
!       CALL getenv(TRIM(env_name), thread_num)
!
!       IF (thread_num /= ' ') THEN
! #else
!       CALL get_environment_variable(name=TRIM(env_name), value=thread_num, &
!           status=istat)
!       IF (istat == 0) THEN
! #endif
!          READ(thread_num,*) global_no_of_threads
!        ELSE
!          WRITE (nerr,'(1x,a,/)') ' Global number of OpenMP threads not given!'
!          WRITE (nerr,'(1x,a,a,a,/,1x,a,a,a)') &
!              ' Environment variable ', TRIM(env_name), ' either not set,', &
!              ' or not available to ', TRIM(yname), ' root PE.'
!          global_no_of_threads = omp_get_max_threads()
!        ENDIF
!     ENDIF ! (my_global_mpi_id == 0)
!
! #ifndef NOMPI
!     ! Make number of threads from environment available to all model PEs
!     CALL MPI_BCAST (global_no_of_threads, 1, MPI_INTEGER, 0, global_mpi_communicator, p_error)
! #endif
!     ! Inform on OpenMP thread usage
! !     CALL OMP_SET_NUM_THREADS(threads)
! !     threads = OMP_GET_MAX_THREADS()
!
!     IF (my_global_mpi_id == 0) THEN
!        WRITE (nerr,*)
!        WRITE (nerr,'(1x,a,i3)') ' global_no_of_threads is ', global_no_of_threads
!     ENDIF
! #endif


    ! by default, the global communicator is the process communicator
    CALL set_process_mpi_name(global_mpi_name)
    CALL set_process_mpi_communicator(global_mpi_communicator)

  END SUBROUTINE start_mpi
  !------------------------------------------------------------------------------

  !------------------------------------------------------------------------------
  SUBROUTINE stop_mpi

    ! finish MPI and clean up all PEs

#ifndef NOMPI
    ! to prevent abort due to unfinished communication
    CALL p_barrier(process_mpi_all_comm)

#ifdef COUP_OAS_ICON
    CALL oasis_terminate(oas_error)
    IF (oas_error /= 0) CALL oasis_abort(oas_comp_id, oas_comp_name, &
      'Failure in oasis_terminate')
#else
    CALL MPI_FINALIZE (p_error)
#endif

    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a)') ' MPI_FINALIZE failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL abort_mpi
    END IF
    process_is_mpi_parallel = .FALSE.
    DEALLOCATE(p_request)
#endif

  END SUBROUTINE stop_mpi
  !------------------------------------------------------------------------------

  !------------------------------------------------------------------------------
  SUBROUTINE abort_mpi

    ! this routine should be used instead of abort, util_abort() or STOP
    ! in all routines for proper clean up of all PEs

#ifndef __STANDALONE
    EXTERNAL :: util_exit
#endif

#ifndef NOMPI
#ifdef COUP_OAS_ICON
    CALL oasis_abort(oas_comp_id, oas_comp_name, 'Abort mpi through oasis')
#else
    CALL MPI_ABORT (MPI_COMM_WORLD, 0, p_error)
#endif

    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a)') ' MPI_ABORT failed.'
       WRITE (nerr,'(a,i4)') ' Error =  ', p_error
       STOP
    END IF
#else
#ifndef __STANDALONE
    CALL util_exit(1)
#else
    STOP 'mo_mpi: abort_mpi ..'
#endif
#endif

  END SUBROUTINE abort_mpi
  !------------------------------------------------------------------------------

  ! communicator set up
  SUBROUTINE p_set_communicator (nproca, nprocb, mapmesh, debug_parallel)

    INTEGER, INTENT(in) :: nproca, nprocb
    INTEGER, INTENT(in) :: mapmesh(0:,0:)
    INTEGER, INTENT(in) :: debug_parallel

#ifndef NOMPI
    INTEGER :: all_debug_pes(SIZE(mapmesh))

    INTEGER :: group_world, group_a, group_b, group_d
    INTEGER :: p_communicator_tmp

    INTEGER :: n, members

    INTEGER :: ranks(1) = 0

    ! first set global group

    CALL MPI_COMM_GROUP (process_mpi_all_comm, group_world, p_error)

    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a)') ' PE: ', my_process_mpi_all_id, &
         & ' MPI_COMM_GROUP failed.'
       WRITE (nerr,'(a,i4)') ' Error =  ', p_error
       CALL abort_mpi
    END IF

    ! communicator is process_mpi_all_comm

    IF (debug_parallel >= 0 ) THEN

       CALL MPI_GROUP_INCL (group_world, 1, ranks, group_d, p_error)

       IF (p_error /= MPI_SUCCESS) THEN
          WRITE (nerr,'(a,i4,a)') ' PE: ', my_process_mpi_all_id, &
            &' MPI_GROUP_INCL failed.'
          WRITE (nerr,'(a,i4)') ' Error =  ', p_error
          CALL abort_mpi
       END IF

       CALL MPI_COMM_CREATE (process_mpi_all_comm, group_d, p_communicator_tmp, &
            p_error)

       IF (p_error /= MPI_SUCCESS) THEN
          WRITE (nerr,'(a,i4,a)') ' PE: ', my_process_mpi_all_id, &
            &' MPI_COMM_CREATE failed.'
          WRITE (nerr,'(a,i4)') ' Error =  ', p_error
          CALL abort_mpi
       END IF

       IF (my_process_mpi_all_id == 0) p_communicator_d = p_communicator_tmp

       DO n = 1, SIZE(mapmesh)
          all_debug_pes(n) = n
       END DO

       CALL MPI_GROUP_INCL (group_world, SIZE(mapmesh), all_debug_pes, &
            group_d, p_error)

       IF (p_error /= MPI_SUCCESS) THEN
          WRITE (nerr,'(a,i4,a)') ' PE: ', my_process_mpi_all_id, &
            & ' MPI_GROUP_INCL failed.'
          WRITE (nerr,'(a,i4)') ' Error =  ', p_error
          CALL abort_mpi
       END IF

       CALL MPI_COMM_CREATE (process_mpi_all_comm, group_d, p_communicator_tmp, &
            p_error)

       IF (p_error /= MPI_SUCCESS) THEN
          WRITE (nerr,'(a,i4,a)') ' PE: ', my_process_mpi_all_id, &
            & ' MPI_COMM_CREATE failed.'
          WRITE (nerr,'(a,i4)') ' Error =  ', p_error
          CALL abort_mpi
       END IF

       IF (my_process_mpi_all_id /= 0) p_communicator_d = p_communicator_tmp

    ELSE
       p_communicator_d = process_mpi_all_comm
    END IF

    DO n = 0, nproca-1
       members = nprocb
       CALL MPI_GROUP_INCL (group_world, members, mapmesh(:,n), group_a, &
            p_error)

       IF (p_error /= MPI_SUCCESS) THEN
          WRITE (nerr,'(a,i4,a)') ' PE: ', my_process_mpi_all_id, &
            & ' MPI_GROUP_INCL failed.'
          WRITE (nerr,'(a,i4)') ' Error =  ', p_error
          CALL abort_mpi
       END IF

       CALL MPI_COMM_CREATE (process_mpi_all_comm, group_a, p_communicator_tmp, &
            p_error)

       IF (p_error /= MPI_SUCCESS) THEN
          WRITE (nerr,'(a,i4,a)') ' PE: ', my_process_mpi_all_id, &
            & ' MPI_COMM_CREATE failed.'
          WRITE (nerr,'(a,i4)') ' Error =  ', p_error
          CALL abort_mpi
       END IF
       IF(p_communicator_tmp/=MPI_COMM_NULL) &
         p_communicator_a = p_communicator_tmp

    END DO

    ! create groups for set Bs

    DO n = 0, nprocb-1
       members = nproca
       CALL MPI_GROUP_INCL (group_world, members, mapmesh(n,:), group_b, &
            p_error)

       IF (p_error /= MPI_SUCCESS) THEN
          WRITE (nerr,'(a,i4,a)') ' PE: ', my_process_mpi_all_id, &
            & ' MPI_GROUP_INCL failed.'
          WRITE (nerr,'(a,i4)') ' Error =  ', p_error
          CALL abort_mpi
       END IF

       CALL MPI_COMM_CREATE (process_mpi_all_comm, group_b, p_communicator_tmp, &
            p_error)

       IF (p_error /= MPI_SUCCESS) THEN
          WRITE (nerr,'(a,i4,a)') ' PE: ', my_process_mpi_all_id, &
            & ' MPI_COMM_CREATE failed.'
          WRITE (nerr,'(a,i4)') ' Error =  ', p_error
          CALL abort_mpi
       END IF
       IF(p_communicator_tmp/=MPI_COMM_NULL) &
         p_communicator_b = p_communicator_tmp

    END DO

    CALL MPI_BARRIER (process_mpi_all_comm, p_error)

    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a)') ' PE: ', my_process_mpi_all_id, '&
         & MPI_BARRIER failed.'
       WRITE (nerr,'(a,i4)') ' Error =  ', p_error
       CALL abort_mpi
    END IF

    IF (debug_parallel >= 0 .AND. my_process_mpi_all_id == 0) THEN
      p_communicator_a = p_communicator_d
      p_communicator_b = p_communicator_d
    ENDIF

#ifdef DEBUG
    WRITE (nerr,'(a,i4,a,3i8)') &
         'p_set_communicator on PE ', my_process_mpi_all_id, ': ', &
         p_communicator_d, &
         p_communicator_a, &
         p_communicator_b
#endif
#endif
  END SUBROUTINE p_set_communicator

!=========================================================================

  !---------------------------------------------------------------------------------------------------------------------------------
  !> wrapper for MPI_Comm_size()
  !---------------------------------------------------------------------------------------------------------------------------------
  FUNCTION p_comm_size(communicator)
    INTEGER :: p_comm_size
    INTEGER, INTENT(IN) :: communicator

    CHARACTER(LEN=*), PARAMETER :: routine = 'p_comm_size'
    INTEGER :: ierr

#ifndef NOMPI
    CALL MPI_Comm_size(communicator, p_comm_size, ierr)
    IF(ierr /= MPI_SUCCESS) CALL finish(routine, 'Error in MPI_COMM_SIZE operation!')
#else
    p_comm_size = 1
#endif
  END FUNCTION p_comm_size


  !---------------------------------------------------------------------------------------------------------------------------------
  !> wrapper for MPI_Comm_rank()
  !---------------------------------------------------------------------------------------------------------------------------------
  FUNCTION p_comm_rank(communicator)
    INTEGER :: p_comm_rank
    INTEGER, INTENT(IN) :: communicator

    CHARACTER(LEN=*), PARAMETER :: routine = 'p_comm_rank'
    INTEGER :: ierr

#ifndef NOMPI
    CALL MPI_COMM_RANK(communicator, p_comm_rank, ierr)
    IF(ierr /= MPI_SUCCESS) CALL finish(routine, 'Error in MPI_COMM_RANK operation!')
#else
    p_comm_rank = 0
#endif
  END FUNCTION p_comm_rank


!=========================================================================

  ! send implementation

  SUBROUTINE p_send_real (t_buffer, p_destination, p_tag, p_count, comm)

    REAL (dp), INTENT(in) :: t_buffer
    INTEGER,   INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

#ifdef __USE_G2G
!$ACC DATA PRESENT( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
!$ACC HOST_DATA USE_DEVICE( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
#endif

    IF (PRESENT(p_count)) THEN
       CALL MPI_SEND (t_buffer, p_count, p_real_dp, p_destination, p_tag, &
            p_comm, p_error)
    ELSE
       CALL MPI_SEND (t_buffer, 1, p_real_dp, p_destination, p_tag, &
            p_comm, p_error)
    END IF

#ifdef __USE_G2G
!$ACC END HOST_DATA
!$ACC END DATA
#endif

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_SEND from ', my_process_mpi_all_id, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL abort_mpi
    END IF
#endif
#endif

  END SUBROUTINE p_send_real


  ! send implementation

  SUBROUTINE p_send_sreal (t_buffer, p_destination, p_tag, p_count, comm)

    REAL (sp), INTENT(in) :: t_buffer
    INTEGER,   INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

#ifdef __USE_G2G
!$ACC DATA PRESENT( t_buffer ), IF ( i_am_accel_node )
!$ACC HOST_DATA USE_DEVICE( t_buffer ), IF ( i_am_accel_node )
#endif

    IF (PRESENT(p_count)) THEN
       CALL MPI_SEND (t_buffer, p_count, p_real_sp, p_destination, p_tag, &
            p_comm, p_error)
    ELSE
       CALL MPI_SEND (t_buffer, 1, p_real_sp, p_destination, p_tag, &
            p_comm, p_error)
    END IF

#ifdef __USE_G2G
!$ACC END HOST_DATA
!$ACC END DATA
#endif

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_SEND from ', my_process_mpi_all_id, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL abort_mpi
    END IF
#endif
#endif

  END SUBROUTINE p_send_sreal


  SUBROUTINE p_send_real_1d (t_buffer, p_destination, p_tag, p_count, comm)

    REAL (dp), INTENT(in) :: t_buffer(:)
    INTEGER,   INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

#ifdef __USE_G2G
!$ACC DATA PRESENT( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
!$ACC HOST_DATA USE_DEVICE( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
#endif

    IF (PRESENT(p_count)) THEN
       CALL MPI_SEND (t_buffer, p_count, p_real_dp, p_destination, p_tag, &
            p_comm, p_error)
    ELSE
       CALL MPI_SEND (t_buffer, SIZE(t_buffer), p_real_dp, p_destination, p_tag, &
            p_comm, p_error)
    END IF

#ifdef __USE_G2G
!$ACC END HOST_DATA
!$ACC END DATA
#endif

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_SEND from ', my_process_mpi_all_id, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL abort_mpi
    END IF
#endif
#endif

  END SUBROUTINE p_send_real_1d

  SUBROUTINE p_send_real_2d (t_buffer, p_destination, p_tag, p_count, comm)

    REAL (dp), INTENT(in) :: t_buffer(:,:)
    INTEGER,   INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

#ifdef __USE_G2G
!$ACC DATA PRESENT( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
!$ACC HOST_DATA USE_DEVICE( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
#endif

    IF (PRESENT(p_count)) THEN
       CALL MPI_SEND (t_buffer, p_count, p_real_dp, p_destination, p_tag, &
            p_comm, p_error)
    ELSE
       CALL MPI_SEND (t_buffer, SIZE(t_buffer), p_real_dp, p_destination, p_tag, &
            p_comm, p_error)
    END IF

#ifdef __USE_G2G
!$ACC END HOST_DATA
!$ACC END DATA
#endif

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_SEND from ', my_process_mpi_all_id, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL abort_mpi
    END IF
#endif
#endif

  END SUBROUTINE p_send_real_2d

  SUBROUTINE p_send_real_3d (t_buffer, p_destination, p_tag, p_count, comm)

    REAL (dp), INTENT(in) :: t_buffer(:,:,:)
    INTEGER,   INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

#ifdef __USE_G2G
!$ACC DATA PRESENT( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
!$ACC HOST_DATA USE_DEVICE( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
#endif

    IF (PRESENT(p_count)) THEN
       CALL MPI_SEND (t_buffer, p_count, p_real_dp, p_destination, p_tag, &
            p_comm, p_error)
    ELSE
       CALL MPI_SEND (t_buffer, SIZE(t_buffer), p_real_dp, p_destination, p_tag, &
            p_comm, p_error)
    END IF

#ifdef __USE_G2G
!$ACC END HOST_DATA
!$ACC END DATA
#endif

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_SEND from ', my_process_mpi_all_id, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL abort_mpi
    END IF
#endif
#endif

  END SUBROUTINE p_send_real_3d

  SUBROUTINE p_send_real_4d (t_buffer, p_destination, p_tag, p_count, comm)

    REAL (dp), INTENT(in) :: t_buffer(:,:,:,:)
    INTEGER,   INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

#ifdef __USE_G2G
!$ACC DATA PRESENT( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
!$ACC HOST_DATA USE_DEVICE( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
#endif

    IF (PRESENT(p_count)) THEN
       CALL MPI_SEND (t_buffer, p_count, p_real_dp, p_destination, p_tag, &
            p_comm, p_error)
    ELSE
       CALL MPI_SEND (t_buffer, SIZE(t_buffer), p_real_dp, p_destination, p_tag, &
            p_comm, p_error)
    END IF

#ifdef __USE_G2G
!$ACC END HOST_DATA
!$ACC END DATA
#endif

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_SEND from ', my_process_mpi_all_id, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL abort_mpi
    END IF
#endif
#endif

  END SUBROUTINE p_send_real_4d

  SUBROUTINE p_send_real_5d (t_buffer, p_destination, p_tag, p_count, comm)

    REAL (dp), INTENT(in) :: t_buffer(:,:,:,:,:)
    INTEGER,   INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

#ifdef __USE_G2G
!$ACC DATA PRESENT( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
!$ACC HOST_DATA USE_DEVICE( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
#endif

    IF (PRESENT(p_count)) THEN
       CALL MPI_SEND (t_buffer, p_count, p_real_dp, p_destination, p_tag, &
            p_comm, p_error)
    ELSE
       CALL MPI_SEND (t_buffer, SIZE(t_buffer), p_real_dp, p_destination, p_tag, &
            p_comm, p_error)
    END IF

#ifdef __USE_G2G
!$ACC END HOST_DATA
!$ACC END DATA
#endif

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_SEND from ', my_process_mpi_all_id, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL abort_mpi
    END IF
#endif
#endif

  END SUBROUTINE p_send_real_5d

  SUBROUTINE p_send_int (t_buffer, p_destination, p_tag, p_count, comm)

    INTEGER, INTENT(in) :: t_buffer
    INTEGER, INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

#ifdef __USE_G2G
!$ACC DATA PRESENT( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
!$ACC HOST_DATA USE_DEVICE( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
#endif

    IF (PRESENT(p_count)) THEN
       CALL MPI_SEND (t_buffer, p_count, p_int, p_destination, p_tag, &
            p_comm, p_error)
    ELSE
       CALL MPI_SEND (t_buffer, 1, p_int, p_destination, p_tag, &
            p_comm, p_error)
    END IF

#ifdef __USE_G2G
!$ACC END HOST_DATA
!$ACC END DATA
#endif

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_SEND from ', my_process_mpi_all_id, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL abort_mpi
    END IF
#endif
#endif

  END SUBROUTINE p_send_int

  SUBROUTINE p_send_int_1d (t_buffer, p_destination, p_tag, p_count, comm)

    INTEGER, INTENT(in) :: t_buffer(:)
    INTEGER, INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

#ifdef __USE_G2G
!$ACC DATA PRESENT( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
!$ACC HOST_DATA USE_DEVICE( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
#endif

    IF (PRESENT(p_count)) THEN
       CALL MPI_SEND (t_buffer, p_count, p_int, p_destination, p_tag, &
            p_comm, p_error)
    ELSE
       CALL MPI_SEND (t_buffer, SIZE(t_buffer), p_int, p_destination, p_tag, &
            p_comm, p_error)
    END IF

#ifdef __USE_G2G
!$ACC END HOST_DATA
!$ACC END DATA
#endif

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_SEND from ', my_process_mpi_all_id, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL abort_mpi
    END IF
#endif
#endif

  END SUBROUTINE p_send_int_1d

  SUBROUTINE p_send_int_2d (t_buffer, p_destination, p_tag, p_count, comm)

    INTEGER, INTENT(in) :: t_buffer(:,:)
    INTEGER, INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

#ifdef __USE_G2G
!$ACC DATA PRESENT( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
!$ACC HOST_DATA USE_DEVICE( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
#endif

    IF (PRESENT(p_count)) THEN
       CALL MPI_SEND (t_buffer, p_count, p_int, p_destination, p_tag, &
            p_comm, p_error)
    ELSE
       CALL MPI_SEND (t_buffer, SIZE(t_buffer), p_int, p_destination, p_tag, &
            p_comm, p_error)
    END IF

#ifdef __USE_G2G
!$ACC END HOST_DATA
!$ACC END DATA
#endif

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_SEND from ', my_process_mpi_all_id, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL abort_mpi
    END IF
#endif
#endif

  END SUBROUTINE p_send_int_2d

  SUBROUTINE p_send_int_3d (t_buffer, p_destination, p_tag, p_count, comm)

    INTEGER, INTENT(in) :: t_buffer(:,:,:)
    INTEGER, INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

#ifdef __USE_G2G
!$ACC DATA PRESENT( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
!$ACC HOST_DATA USE_DEVICE( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
#endif

    IF (PRESENT(p_count)) THEN
       CALL MPI_SEND (t_buffer, p_count, p_int, p_destination, p_tag, &
            p_comm, p_error)
    ELSE
       CALL MPI_SEND (t_buffer, SIZE(t_buffer), p_int, p_destination, p_tag, &
            p_comm, p_error)
    END IF

#ifdef __USE_G2G
!$ACC END HOST_DATA
!$ACC END DATA
#endif

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_SEND from ', my_process_mpi_all_id, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL abort_mpi
    END IF
#endif
#endif

  END SUBROUTINE p_send_int_3d

  SUBROUTINE p_send_int_4d (t_buffer, p_destination, p_tag, p_count, comm)

    INTEGER, INTENT(in) :: t_buffer(:,:,:,:)
    INTEGER, INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

#ifdef __USE_G2G
!$ACC DATA PRESENT( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
!$ACC HOST_DATA USE_DEVICE( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
#endif

    IF (PRESENT(p_count)) THEN
       CALL MPI_SEND (t_buffer, p_count, p_int, p_destination, p_tag, &
            p_comm, p_error)
    ELSE
       CALL MPI_SEND (t_buffer, SIZE(t_buffer), p_int, p_destination, p_tag, &
            p_comm, p_error)
    END IF

#ifdef __USE_G2G
!$ACC END HOST_DATA
!$ACC END DATA
#endif

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_SEND from ', my_process_mpi_all_id, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL abort_mpi
    END IF
#endif
#endif

  END SUBROUTINE p_send_int_4d


  SUBROUTINE p_send_bool (t_buffer, p_destination, p_tag, p_count, comm)

    LOGICAL, INTENT(in) :: t_buffer
    INTEGER, INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

#ifdef __USE_G2G
!$ACC DATA PRESENT( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
!$ACC HOST_DATA USE_DEVICE( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
#endif

    IF (PRESENT(p_count)) THEN
       CALL MPI_SEND (t_buffer, p_count, p_bool, p_destination, p_tag, &
            p_comm, p_error)
    ELSE
       CALL MPI_SEND (t_buffer, 1, p_bool, p_destination, p_tag, &
            p_comm, p_error)
    END IF

#ifdef __USE_G2G
!$ACC END HOST_DATA
!$ACC END DATA
#endif

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_SEND from ', my_process_mpi_all_id, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL abort_mpi
    END IF
#endif
#endif

  END SUBROUTINE p_send_bool

  SUBROUTINE p_send_bool_1d (t_buffer, p_destination, p_tag, p_count, comm)

    LOGICAL, INTENT(in) :: t_buffer(:)
    INTEGER, INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

#ifdef __USE_G2G
!$ACC DATA PRESENT( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
!$ACC HOST_DATA USE_DEVICE( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
#endif

    IF (PRESENT(p_count)) THEN
       CALL MPI_SEND (t_buffer, p_count, p_bool, p_destination, p_tag, &
            p_comm, p_error)
    ELSE
       CALL MPI_SEND (t_buffer, SIZE(t_buffer), p_bool, p_destination, p_tag, &
            p_comm, p_error)
    END IF

#ifdef __USE_G2G
!$ACC END HOST_DATA
!$ACC END DATA
#endif

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_SEND from ', my_process_mpi_all_id, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL abort_mpi
    END IF
#endif
#endif

  END SUBROUTINE p_send_bool_1d


  SUBROUTINE p_send_char_1d (t_buffer, p_destination, p_tag, p_count, comm)

    CHARACTER(LEN=*), INTENT(in) :: t_buffer(:)
    INTEGER, INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

#ifdef __USE_G2G
!$ACC DATA PRESENT( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
!$ACC HOST_DATA USE_DEVICE( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
#endif

    IF (PRESENT(p_count)) THEN
       CALL MPI_SEND (t_buffer, LEN(t_buffer(1))*p_count, p_char, p_destination, p_tag, &
            p_comm, p_error)
    ELSE
       CALL MPI_SEND (t_buffer, LEN(t_buffer(1))*SIZE(t_buffer), p_char, p_destination, p_tag, &
            p_comm, p_error)
    END IF

#ifdef __USE_G2G
!$ACC END HOST_DATA
!$ACC END DATA
#endif

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_SEND from ', my_process_mpi_all_id, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL abort_mpi
    END IF
#endif
#endif

  END SUBROUTINE p_send_char_1d


  SUBROUTINE p_send_bool_2d (t_buffer, p_destination, p_tag, p_count, comm)

    LOGICAL, INTENT(in) :: t_buffer(:,:)
    INTEGER, INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

#ifdef __USE_G2G
!$ACC DATA PRESENT( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
!$ACC HOST_DATA USE_DEVICE( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
#endif

    IF (PRESENT(p_count)) THEN
       CALL MPI_SEND (t_buffer, p_count, p_bool, p_destination, p_tag, &
            p_comm, p_error)
    ELSE
       CALL MPI_SEND (t_buffer, SIZE(t_buffer), p_bool, p_destination, p_tag, &
            p_comm, p_error)
    END IF

#ifdef __USE_G2G
!$ACC END HOST_DATA
!$ACC END DATA
#endif

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_SEND from ', my_process_mpi_all_id, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL abort_mpi
    END IF
#endif
#endif

  END SUBROUTINE p_send_bool_2d

  SUBROUTINE p_send_bool_3d (t_buffer, p_destination, p_tag, p_count, comm)

    LOGICAL, INTENT(in) :: t_buffer(:,:,:)
    INTEGER, INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

#ifdef __USE_G2G
!$ACC DATA PRESENT( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
!$ACC HOST_DATA USE_DEVICE( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
#endif

    IF (PRESENT(p_count)) THEN
       CALL MPI_SEND (t_buffer, p_count, p_bool, p_destination, p_tag, &
            p_comm, p_error)
    ELSE
       CALL MPI_SEND (t_buffer, SIZE(t_buffer), p_bool, p_destination, p_tag, &
            p_comm, p_error)
    END IF

#ifdef __USE_G2G
!$ACC END HOST_DATA
!$ACC END DATA
#endif

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_SEND from ', my_process_mpi_all_id, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL abort_mpi
    END IF
#endif
#endif

  END SUBROUTINE p_send_bool_3d

  SUBROUTINE p_send_bool_4d (t_buffer, p_destination, p_tag, p_count, comm)

    LOGICAL, INTENT(in) :: t_buffer(:,:,:,:)
    INTEGER, INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

#ifdef __USE_G2G
!$ACC DATA PRESENT( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
!$ACC HOST_DATA USE_DEVICE( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
#endif

    IF (PRESENT(p_count)) THEN
       CALL MPI_SEND (t_buffer, p_count, p_bool, p_destination, p_tag, &
            p_comm, p_error)
    ELSE
       CALL MPI_SEND (t_buffer, SIZE(t_buffer), p_bool, p_destination, p_tag, &
            p_comm, p_error)
    END IF

#ifdef __USE_G2G
!$ACC END HOST_DATA
!$ACC END DATA
#endif

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_SEND from ', my_process_mpi_all_id, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL abort_mpi
    END IF
#endif
#endif

  END SUBROUTINE p_send_bool_4d

  SUBROUTINE p_send_char (t_buffer, p_destination, p_tag, p_count, comm)

    CHARACTER(len=*),  INTENT(in) :: t_buffer
    INTEGER,           INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

#ifdef __USE_G2G
!$ACC DATA PRESENT( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
!$ACC HOST_DATA USE_DEVICE( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
#endif

    IF (PRESENT(p_count)) THEN
       CALL MPI_SEND (t_buffer, p_count, p_char, p_destination, p_tag, &
            p_comm, p_error)
    ELSE
       CALL MPI_SEND (t_buffer, LEN(t_buffer), p_char, p_destination, p_tag, &
            p_comm, p_error)
    END IF

#ifdef __USE_G2G
!$ACC END HOST_DATA
!$ACC END DATA
#endif

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_SEND from ', my_process_mpi_all_id, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL abort_mpi
    END IF
#endif
#endif

  END SUBROUTINE p_send_char

! non-blocking sends

  SUBROUTINE p_inc_request
    INTEGER, ALLOCATABLE :: tmp(:)

#ifndef NOMPI
    p_irequest = p_irequest + 1
    IF (p_irequest > p_mrequest) THEN
      ALLOCATE(tmp(p_mrequest))
      tmp(:) = p_request(:)
      DEALLOCATE(p_request)
      ALLOCATE(p_request(p_mrequest+p_request_alloc_size))
      p_request(1:p_mrequest) = tmp(:)
      p_mrequest = p_mrequest+p_request_alloc_size
      DEALLOCATE(tmp)
    ENDIF
#endif

  END SUBROUTINE p_inc_request

  SUBROUTINE p_isend_real (t_buffer, p_destination, p_tag, p_count, comm)

    REAL (dp), INTENT(inout) :: t_buffer
    INTEGER,   INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

#ifdef __USE_G2G
!$ACC DATA PRESENT( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
!$ACC HOST_DATA USE_DEVICE( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
#endif

    IF (PRESENT(p_count)) THEN
       CALL p_inc_request
       CALL MPI_ISEND (t_buffer, p_count, p_real_dp, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    ELSE
       CALL p_inc_request
       CALL MPI_ISEND (t_buffer, 1, p_real_dp, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    END IF

#ifdef __USE_G2G
!$ACC END HOST_DATA
!$ACC END DATA
#endif

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_ISEND from ', my_process_mpi_all_id, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL abort_mpi
    END IF
#endif
#endif

  END SUBROUTINE p_isend_real


  SUBROUTINE p_isend_sreal (t_buffer, p_destination, p_tag, p_count, comm)

    REAL (sp), INTENT(inout) :: t_buffer
    INTEGER,   INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

#ifdef __USE_G2G
!$ACC DATA PRESENT( t_buffer ), IF ( i_am_accel_node )
!$ACC HOST_DATA USE_DEVICE( t_buffer ), IF ( i_am_accel_node )
#endif

    IF (PRESENT(p_count)) THEN
       CALL p_inc_request
       CALL MPI_ISEND (t_buffer, p_count, p_real_sp, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    ELSE
       CALL p_inc_request
       CALL MPI_ISEND (t_buffer, 1, p_real_sp, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    END IF

#ifdef __USE_G2G
!$ACC END HOST_DATA
!$ACC END DATA
#endif

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_ISEND from ', my_process_mpi_all_id, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL abort_mpi
    END IF
#endif
#endif

  END SUBROUTINE p_isend_sreal


  SUBROUTINE p_isend_real_1d (t_buffer, p_destination, p_tag, p_count, comm, request)

    REAL (dp), INTENT(inout) :: t_buffer(:)
    INTEGER,   INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
    INTEGER,   INTENT(INOUT), OPTIONAL :: request
#ifndef NOMPI
    INTEGER :: p_comm, icount, out_request

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

    icount = SIZE(t_buffer)
    IF (PRESENT(p_count))  icount = p_count

#ifdef __USE_G2G
!$ACC DATA PRESENT( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
!$ACC HOST_DATA USE_DEVICE( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
#endif

    CALL p_inc_request
    CALL MPI_ISEND (t_buffer, icount, p_real_dp, p_destination, p_tag, &
      &             p_comm, out_request, p_error)

#ifdef __USE_G2G
!$ACC END HOST_DATA
!$ACC END DATA
#endif

    IF (PRESENT(request)) THEN
      request               = out_request
    ELSE
      p_request(p_irequest) = out_request
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_ISEND from ', my_process_mpi_all_id, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL abort_mpi
    END IF
#endif
#endif

  END SUBROUTINE p_isend_real_1d

  SUBROUTINE p_isend_real_2d (t_buffer, p_destination, p_tag, p_count, comm)

    REAL (dp), INTENT(inout) :: t_buffer(:,:)
    INTEGER,   INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

#ifdef __USE_G2G
!$ACC DATA PRESENT( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
!$ACC HOST_DATA USE_DEVICE( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
#endif

    IF (PRESENT(p_count)) THEN
       CALL p_inc_request
       CALL MPI_ISEND (t_buffer, p_count, p_real_dp, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    ELSE
       CALL p_inc_request
       CALL MPI_ISEND (t_buffer, SIZE(t_buffer), p_real_dp, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    END IF

#ifdef __USE_G2G
!$ACC END HOST_DATA
!$ACC END DATA
#endif

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_ISEND from ', my_process_mpi_all_id, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL abort_mpi
    END IF
#endif
#endif

  END SUBROUTINE p_isend_real_2d


  SUBROUTINE p_isend_sreal_2d (t_buffer, p_destination, p_tag, p_count, comm)

    REAL (sp), INTENT(inout) :: t_buffer(:,:)
    INTEGER,   INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

#ifdef __USE_G2G
!$ACC DATA PRESENT( t_buffer ), IF ( i_am_accel_node )
!$ACC HOST_DATA USE_DEVICE( t_buffer ), IF ( i_am_accel_node )
#endif

    IF (PRESENT(p_count)) THEN
       CALL p_inc_request
       CALL MPI_ISEND (t_buffer, p_count, p_real_sp, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    ELSE
       CALL p_inc_request
       CALL MPI_ISEND (t_buffer, SIZE(t_buffer), p_real_sp, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    END IF

#ifdef __USE_G2G
!$ACC END HOST_DATA
!$ACC END DATA
#endif

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_ISEND from ', my_process_mpi_all_id, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL abort_mpi
    END IF
#endif
#endif

  END SUBROUTINE p_isend_sreal_2d


  SUBROUTINE p_isend_real_3d (t_buffer, p_destination, p_tag, p_count, comm)

    REAL (dp), INTENT(inout) :: t_buffer(:,:,:)
    INTEGER,   INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

#ifdef __USE_G2G
!$ACC DATA PRESENT( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
!$ACC HOST_DATA USE_DEVICE( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
#endif

    IF (PRESENT(p_count)) THEN
       CALL p_inc_request
       CALL MPI_ISEND (t_buffer, p_count, p_real_dp, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    ELSE
       CALL p_inc_request
       CALL MPI_ISEND (t_buffer, SIZE(t_buffer), p_real_dp, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    END IF

#ifdef __USE_G2G
!$ACC END HOST_DATA
!$ACC END DATA
#endif

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_ISEND from ', my_process_mpi_all_id, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL abort_mpi
    END IF
#endif
#endif

  END SUBROUTINE p_isend_real_3d

  SUBROUTINE p_isend_real_4d (t_buffer, p_destination, p_tag, p_count, comm)

    REAL (dp), INTENT(inout) :: t_buffer(:,:,:,:)
    INTEGER,   INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

#ifdef __USE_G2G
!$ACC DATA PRESENT( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
!$ACC HOST_DATA USE_DEVICE( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
#endif

    IF (PRESENT(p_count)) THEN
       CALL p_inc_request
       CALL MPI_ISEND (t_buffer, p_count, p_real_dp, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    ELSE
       CALL p_inc_request
       CALL MPI_ISEND (t_buffer, SIZE(t_buffer), p_real_dp, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    END IF

#ifdef __USE_G2G
!$ACC END HOST_DATA
!$ACC END DATA
#endif

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_ISEND from ', my_process_mpi_all_id, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL abort_mpi
    END IF
#endif
#endif

  END SUBROUTINE p_isend_real_4d

  SUBROUTINE p_isend_real_5d (t_buffer, p_destination, p_tag, p_count, comm)

    REAL (dp), INTENT(inout) :: t_buffer(:,:,:,:,:)
    INTEGER,   INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

#ifdef __USE_G2G
!$ACC DATA PRESENT( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
!$ACC HOST_DATA USE_DEVICE( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
#endif

    IF (PRESENT(p_count)) THEN
       CALL p_inc_request
       CALL MPI_ISEND (t_buffer, p_count, p_real_dp, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    ELSE
       CALL p_inc_request
       CALL MPI_ISEND (t_buffer, SIZE(t_buffer), p_real_dp, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    END IF

#ifdef __USE_G2G
!$ACC END HOST_DATA
!$ACC END DATA
#endif

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_ISEND from ', my_process_mpi_all_id, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL abort_mpi
    END IF
#endif
#endif

  END SUBROUTINE p_isend_real_5d

  SUBROUTINE p_isend_int (t_buffer, p_destination, p_tag, p_count, comm)

    INTEGER, INTENT(inout) :: t_buffer
    INTEGER, INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

#ifdef __USE_G2G
!$ACC DATA PRESENT( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
!$ACC HOST_DATA USE_DEVICE( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
#endif

    IF (PRESENT(p_count)) THEN
       CALL p_inc_request
       CALL MPI_ISEND (t_buffer, p_count, p_int, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    ELSE
       CALL p_inc_request
       CALL MPI_ISEND (t_buffer, 1, p_int, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    END IF

#ifdef __USE_G2G
!$ACC END HOST_DATA
!$ACC END DATA
#endif

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_ISEND from ', my_process_mpi_all_id, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL abort_mpi
    END IF
#endif
#endif

  END SUBROUTINE p_isend_int

  SUBROUTINE p_isend_int_1d (t_buffer, p_destination, p_tag, p_count, comm)

    INTEGER, INTENT(inout) :: t_buffer(:)
    INTEGER, INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

#ifdef __USE_G2G
!$ACC DATA PRESENT( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
!$ACC HOST_DATA USE_DEVICE( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
#endif

    IF (PRESENT(p_count)) THEN
       CALL p_inc_request
       CALL MPI_ISEND (t_buffer, p_count, p_int, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    ELSE
       CALL p_inc_request
       CALL MPI_ISEND (t_buffer, SIZE(t_buffer), p_int, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    END IF

#ifdef __USE_G2G
!$ACC END HOST_DATA
!$ACC END DATA
#endif

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_ISEND from ', my_process_mpi_all_id, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL abort_mpi
    END IF
#endif
#endif

  END SUBROUTINE p_isend_int_1d

  SUBROUTINE p_isend_int_2d (t_buffer, p_destination, p_tag, p_count, comm)

    INTEGER, INTENT(inout) :: t_buffer(:,:)
    INTEGER, INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

#ifdef __USE_G2G
!$ACC DATA PRESENT( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
!$ACC HOST_DATA USE_DEVICE( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
#endif

    IF (PRESENT(p_count)) THEN
       CALL p_inc_request
       CALL MPI_ISEND (t_buffer, p_count, p_int, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    ELSE
       CALL p_inc_request
       CALL MPI_ISEND (t_buffer, SIZE(t_buffer), p_int, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    END IF

#ifdef __USE_G2G
!$ACC END HOST_DATA
!$ACC END DATA
#endif

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_ISEND from ', my_process_mpi_all_id, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL abort_mpi
    END IF
#endif
#endif

  END SUBROUTINE p_isend_int_2d

  SUBROUTINE p_isend_int_3d (t_buffer, p_destination, p_tag, p_count, comm)

    INTEGER, INTENT(inout) :: t_buffer(:,:,:)
    INTEGER, INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

#ifdef __USE_G2G
!$ACC DATA PRESENT( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
!$ACC HOST_DATA USE_DEVICE( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
#endif

    IF (PRESENT(p_count)) THEN
       CALL p_inc_request
       CALL MPI_ISEND (t_buffer, p_count, p_int, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    ELSE
       CALL p_inc_request
       CALL MPI_ISEND (t_buffer, SIZE(t_buffer), p_int, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    END IF

#ifdef __USE_G2G
!$ACC END HOST_DATA
!$ACC END DATA
#endif

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_ISEND from ', my_process_mpi_all_id, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL abort_mpi
    END IF
#endif
#endif

  END SUBROUTINE p_isend_int_3d

  SUBROUTINE p_isend_int_4d (t_buffer, p_destination, p_tag, p_count, comm)

    INTEGER, INTENT(inout) :: t_buffer(:,:,:,:)
    INTEGER, INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

#ifdef __USE_G2G
!$ACC DATA PRESENT( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
!$ACC HOST_DATA USE_DEVICE( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
#endif

    IF (PRESENT(p_count)) THEN
       CALL p_inc_request
       CALL MPI_ISEND (t_buffer, p_count, p_int, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    ELSE
       CALL p_inc_request
       CALL MPI_ISEND (t_buffer, SIZE(t_buffer), p_int, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    END IF

#ifdef __USE_G2G
!$ACC END HOST_DATA
!$ACC END DATA
#endif

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_ISEND from ', my_process_mpi_all_id, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL abort_mpi
    END IF
#endif
#endif

  END SUBROUTINE p_isend_int_4d


  SUBROUTINE p_isend_bool (t_buffer, p_destination, p_tag, p_count, comm)

    LOGICAL, INTENT(inout) :: t_buffer
    INTEGER, INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

#ifdef __USE_G2G
!$ACC DATA PRESENT( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
!$ACC HOST_DATA USE_DEVICE( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
#endif

    IF (PRESENT(p_count)) THEN
       CALL p_inc_request
       CALL MPI_ISEND (t_buffer, p_count, p_bool, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    ELSE
       CALL p_inc_request
       CALL MPI_ISEND (t_buffer, 1, p_bool, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    END IF

#ifdef __USE_G2G
!$ACC END HOST_DATA
!$ACC END DATA
#endif

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_ISEND from ', my_process_mpi_all_id, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL abort_mpi
    END IF
#endif
#endif

  END SUBROUTINE p_isend_bool

  SUBROUTINE p_isend_bool_1d (t_buffer, p_destination, p_tag, p_count, comm)

    LOGICAL, INTENT(inout) :: t_buffer(:)
    INTEGER, INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

#ifdef __USE_G2G
!$ACC DATA PRESENT( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
!$ACC HOST_DATA USE_DEVICE( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
#endif

    IF (PRESENT(p_count)) THEN
       CALL p_inc_request
       CALL MPI_ISEND (t_buffer, p_count, p_bool, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    ELSE
       CALL p_inc_request
       CALL MPI_ISEND (t_buffer, SIZE(t_buffer), p_bool, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    END IF

#ifdef __USE_G2G
!$ACC END HOST_DATA
!$ACC END DATA
#endif

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_ISEND from ', my_process_mpi_all_id, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL abort_mpi
    END IF
#endif
#endif

  END SUBROUTINE p_isend_bool_1d

  SUBROUTINE p_isend_bool_2d (t_buffer, p_destination, p_tag, p_count, comm)

    LOGICAL, INTENT(inout) :: t_buffer(:,:)
    INTEGER, INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

#ifdef __USE_G2G
!$ACC DATA PRESENT( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
!$ACC HOST_DATA USE_DEVICE( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
#endif

    IF (PRESENT(p_count)) THEN
       CALL p_inc_request
       CALL MPI_ISEND (t_buffer, p_count, p_bool, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    ELSE
       CALL p_inc_request
       CALL MPI_ISEND (t_buffer, SIZE(t_buffer), p_bool, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    END IF

#ifdef __USE_G2G
!$ACC END HOST_DATA
!$ACC END DATA
#endif

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_ISEND from ', my_process_mpi_all_id, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL abort_mpi
    END IF
#endif
#endif

  END SUBROUTINE p_isend_bool_2d

  SUBROUTINE p_isend_bool_3d (t_buffer, p_destination, p_tag, p_count, comm)

    LOGICAL, INTENT(inout) :: t_buffer(:,:,:)
    INTEGER, INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

#ifdef __USE_G2G
!$ACC DATA PRESENT( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
!$ACC HOST_DATA USE_DEVICE( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
#endif

    IF (PRESENT(p_count)) THEN
       CALL p_inc_request
       CALL MPI_ISEND (t_buffer, p_count, p_bool, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    ELSE
       CALL p_inc_request
       CALL MPI_ISEND (t_buffer, SIZE(t_buffer), p_bool, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    END IF

#ifdef __USE_G2G
!$ACC END HOST_DATA
!$ACC END DATA
#endif

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_ISEND from ', my_process_mpi_all_id, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL abort_mpi
    END IF
#endif
#endif

  END SUBROUTINE p_isend_bool_3d

  SUBROUTINE p_isend_bool_4d (t_buffer, p_destination, p_tag, p_count, comm)

    LOGICAL, INTENT(inout) :: t_buffer(:,:,:,:)
    INTEGER, INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

#ifdef __USE_G2G
!$ACC DATA PRESENT( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
!$ACC HOST_DATA USE_DEVICE( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
#endif

    IF (PRESENT(p_count)) THEN
       CALL p_inc_request
       CALL MPI_ISEND (t_buffer, p_count, p_bool, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    ELSE
       CALL p_inc_request
       CALL MPI_ISEND (t_buffer, SIZE(t_buffer), p_bool, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    END IF

#ifdef __USE_G2G
!$ACC END HOST_DATA
!$ACC END DATA
#endif

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_ISEND from ', my_process_mpi_all_id, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL abort_mpi
    END IF
#endif
#endif

  END SUBROUTINE p_isend_bool_4d

  SUBROUTINE p_isend_char (t_buffer, p_destination, p_tag, p_count, comm)

    CHARACTER(len=*),  INTENT(inout) :: t_buffer
    INTEGER,           INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

#ifdef __USE_G2G
!$ACC DATA PRESENT( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
!$ACC HOST_DATA USE_DEVICE( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
#endif

    IF (PRESENT(p_count)) THEN
       CALL p_inc_request
       CALL MPI_ISEND (t_buffer, p_count, p_char, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    ELSE
       CALL p_inc_request
       CALL MPI_ISEND (t_buffer, LEN(t_buffer), p_char, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    END IF

#ifdef __USE_G2G
!$ACC END HOST_DATA
!$ACC END DATA
#endif

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_ISEND from ', my_process_mpi_all_id, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL abort_mpi
    END IF
#endif
#endif

  END SUBROUTINE p_isend_char

  ! recv implementation

  SUBROUTINE p_recv_real (t_buffer, p_source, p_tag, p_count, comm)

    REAL (dp), INTENT(out) :: t_buffer
    INTEGER,   INTENT(in)  :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

#ifdef __USE_G2G
!$ACC DATA PRESENT( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
!$ACC HOST_DATA USE_DEVICE( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
#endif

    IF (PRESENT(p_count)) THEN
       CALL MPI_RECV (t_buffer, p_count, p_real_dp, p_source, p_tag, &
            p_comm, p_status, p_error)
    ELSE
       CALL MPI_RECV (t_buffer, 1, p_real_dp, p_source, p_tag, &
            p_comm, p_status, p_error)
    END IF

#ifdef __USE_G2G
!$ACC END HOST_DATA
!$ACC END DATA
#endif

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_RECV on ', my_process_mpi_all_id, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL abort_mpi
    END IF
#endif
#endif

  END SUBROUTINE p_recv_real

  ! recv implementation

  SUBROUTINE p_recv_sreal (t_buffer, p_source, p_tag, p_count, comm)

    REAL (sp), INTENT(out) :: t_buffer
    INTEGER,   INTENT(in)  :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

#ifdef __USE_G2G
!$ACC DATA PRESENT( t_buffer ), IF ( i_am_accel_node )
!$ACC HOST_DATA USE_DEVICE( t_buffer ), IF ( i_am_accel_node )
#endif

    IF (PRESENT(p_count)) THEN
       CALL MPI_RECV (t_buffer, p_count, p_real_sp, p_source, p_tag, &
            p_comm, p_status, p_error)
    ELSE
       CALL MPI_RECV (t_buffer, 1, p_real_sp, p_source, p_tag, &
            p_comm, p_status, p_error)
    END IF

#ifdef __USE_G2G
!$ACC END HOST_DATA
!$ACC END DATA
#endif

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_RECV on ', my_process_mpi_all_id, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL abort_mpi
    END IF
#endif
#endif

  END SUBROUTINE p_recv_sreal

  SUBROUTINE p_recv_real_1d (t_buffer, p_source, p_tag, p_count, comm, displs)

    REAL (dp), INTENT(out) :: t_buffer(:)
    INTEGER,   INTENT(in)  :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm, displs
#ifndef NOMPI
    INTEGER :: p_comm, icount, idispls

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

    icount = SIZE(t_buffer)
    IF (PRESENT(p_count))  icount = p_count

    idispls = 1
    IF (PRESENT(displs)) idispls = displs

#ifdef __USE_G2G
!$ACC DATA PRESENT( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
!$ACC HOST_DATA USE_DEVICE( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
#endif

    CALL MPI_RECV (t_buffer(idispls), icount, p_real_dp, p_source, p_tag, &
      &            p_comm, p_status, p_error)

#ifdef __USE_G2G
!$ACC END HOST_DATA
!$ACC END DATA
#endif

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_RECV on ', my_process_mpi_all_id, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL abort_mpi
    END IF
#endif
#endif

  END SUBROUTINE p_recv_real_1d

  SUBROUTINE p_recv_real_2d (t_buffer, p_source, p_tag, p_count, comm)

    REAL (dp), INTENT(out) :: t_buffer(:,:)
    INTEGER,   INTENT(in)  :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

#ifdef __USE_G2G
!$ACC DATA PRESENT( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
!$ACC HOST_DATA USE_DEVICE( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
#endif

    IF (PRESENT(p_count)) THEN
       CALL MPI_RECV (t_buffer, p_count, p_real_dp, p_source, p_tag, &
            p_comm, p_status, p_error)
    ELSE
       CALL MPI_RECV (t_buffer, SIZE(t_buffer), p_real_dp, p_source, p_tag, &
            p_comm, p_status, p_error)
    END IF

#ifdef __USE_G2G
!$ACC END HOST_DATA
!$ACC END DATA
#endif

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_RECV on ', my_process_mpi_all_id, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL abort_mpi
    END IF
#endif
#endif

  END SUBROUTINE p_recv_real_2d

  SUBROUTINE p_recv_real_3d (t_buffer, p_source, p_tag, p_count, comm)

    REAL (dp), INTENT(out) :: t_buffer(:,:,:)
    INTEGER,   INTENT(in)  :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

#ifdef __USE_G2G
!$ACC DATA PRESENT( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
!$ACC HOST_DATA USE_DEVICE( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
#endif

    IF (PRESENT(p_count)) THEN
       CALL MPI_RECV (t_buffer, p_count, p_real_dp, p_source, p_tag, &
            p_comm, p_status, p_error)
    ELSE
       CALL MPI_RECV (t_buffer, SIZE(t_buffer), p_real_dp, p_source, p_tag, &
            p_comm, p_status, p_error)
    END IF

#ifdef __USE_G2G
!$ACC END HOST_DATA
!$ACC END DATA
#endif

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_RECV on ', my_process_mpi_all_id, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL abort_mpi
    END IF
#endif
#endif

  END SUBROUTINE p_recv_real_3d

  SUBROUTINE p_recv_real_4d (t_buffer, p_source, p_tag, p_count, comm)

    REAL (dp), INTENT(out) :: t_buffer(:,:,:,:)
    INTEGER,   INTENT(in)  :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

#ifdef __USE_G2G
!$ACC DATA PRESENT( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
!$ACC HOST_DATA USE_DEVICE( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
#endif

    IF (PRESENT(p_count)) THEN
       CALL MPI_RECV (t_buffer, p_count, p_real_dp, p_source, p_tag, &
            p_comm, p_status, p_error)
    ELSE
       CALL MPI_RECV (t_buffer, SIZE(t_buffer), p_real_dp, p_source, p_tag, &
            p_comm, p_status, p_error)
    END IF

#ifdef __USE_G2G
!$ACC END HOST_DATA
!$ACC END DATA
#endif

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_RECV on ', my_process_mpi_all_id, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL abort_mpi
    END IF
#endif
#endif

  END SUBROUTINE p_recv_real_4d

  SUBROUTINE p_recv_real_5d (t_buffer, p_source, p_tag, p_count, comm)

    REAL (dp), INTENT(out) :: t_buffer(:,:,:,:,:)
    INTEGER,   INTENT(in)  :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

#ifdef __USE_G2G
!$ACC DATA PRESENT( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
!$ACC HOST_DATA USE_DEVICE( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
#endif

    IF (PRESENT(p_count)) THEN
       CALL MPI_RECV (t_buffer, p_count, p_real_dp, p_source, p_tag, &
            p_comm, p_status, p_error)
    ELSE
       CALL MPI_RECV (t_buffer, SIZE(t_buffer), p_real_dp, p_source, p_tag, &
            p_comm, p_status, p_error)
    END IF

#ifdef __USE_G2G
!$ACC END HOST_DATA
!$ACC END DATA
#endif

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_RECV on ', my_process_mpi_all_id, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL abort_mpi
    END IF
#endif
#endif

  END SUBROUTINE p_recv_real_5d

  SUBROUTINE p_recv_int (t_buffer, p_source, p_tag, p_count, comm)

    INTEGER, INTENT(out) :: t_buffer
    INTEGER, INTENT(in)  :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

#ifdef __USE_G2G
!$ACC DATA PRESENT( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
!$ACC HOST_DATA USE_DEVICE( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
#endif

    IF (PRESENT(p_count)) THEN
       CALL MPI_RECV (t_buffer, p_count, p_int, p_source, p_tag, &
            p_comm, p_status, p_error)
    ELSE
       CALL MPI_RECV (t_buffer, 1, p_int, p_source, p_tag, &
            p_comm, p_status, p_error)
    END IF

#ifdef __USE_G2G
!$ACC END HOST_DATA
!$ACC END DATA
#endif

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_RECV on ', my_process_mpi_all_id, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL abort_mpi
    END IF
#endif
#endif

  END SUBROUTINE p_recv_int

  SUBROUTINE p_recv_int_1d (t_buffer, p_source, p_tag, p_count, comm)

    INTEGER, INTENT(out) :: t_buffer(:)
    INTEGER, INTENT(in)  :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

#ifdef __USE_G2G
!$ACC DATA PRESENT( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
!$ACC HOST_DATA USE_DEVICE( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
#endif

    IF (PRESENT(p_count)) THEN
       CALL MPI_RECV (t_buffer, p_count, p_int, p_source, p_tag, &
            p_comm, p_status, p_error)
    ELSE
       CALL MPI_RECV (t_buffer, SIZE(t_buffer), p_int, p_source, p_tag, &
            p_comm, p_status, p_error)
    END IF

#ifdef __USE_G2G
!$ACC END HOST_DATA
!$ACC END DATA
#endif

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_RECV on ', my_process_mpi_all_id, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL abort_mpi
    END IF
#endif
#endif

  END SUBROUTINE p_recv_int_1d

  SUBROUTINE p_recv_int_2d (t_buffer, p_source, p_tag, p_count, comm)

    INTEGER, INTENT(out) :: t_buffer(:,:)
    INTEGER, INTENT(in)  :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

#ifdef __USE_G2G
!$ACC DATA PRESENT( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
!$ACC HOST_DATA USE_DEVICE( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
#endif

    IF (PRESENT(p_count)) THEN
       CALL MPI_RECV (t_buffer, p_count, p_int, p_source, p_tag, &
            p_comm, p_status, p_error)
    ELSE
       CALL MPI_RECV (t_buffer, SIZE(t_buffer), p_int, p_source, p_tag, &
            p_comm, p_status, p_error)
    END IF

#ifdef __USE_G2G
!$ACC END HOST_DATA
!$ACC END DATA
#endif

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_RECV on ', my_process_mpi_all_id, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL abort_mpi
    END IF
#endif
#endif

  END SUBROUTINE p_recv_int_2d

  SUBROUTINE p_recv_int_3d (t_buffer, p_source, p_tag, p_count, comm)

    INTEGER, INTENT(out) :: t_buffer(:,:,:)
    INTEGER, INTENT(in)  :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

#ifdef __USE_G2G
!$ACC DATA PRESENT( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
!$ACC HOST_DATA USE_DEVICE( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
#endif

    IF (PRESENT(p_count)) THEN
       CALL MPI_RECV (t_buffer, p_count, p_int, p_source, p_tag, &
            p_comm, p_status, p_error)
    ELSE
       CALL MPI_RECV (t_buffer, SIZE(t_buffer), p_int, p_source, p_tag, &
            p_comm, p_status, p_error)
    END IF

#ifdef __USE_G2G
!$ACC END HOST_DATA
!$ACC END DATA
#endif

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_RECV on ', my_process_mpi_all_id, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL abort_mpi
    END IF
#endif
#endif

  END SUBROUTINE p_recv_int_3d

  SUBROUTINE p_recv_int_4d (t_buffer, p_source, p_tag, p_count, comm)

    INTEGER, INTENT(out) :: t_buffer(:,:,:,:)
    INTEGER, INTENT(in)  :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

#ifdef __USE_G2G
!$ACC DATA PRESENT( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
!$ACC HOST_DATA USE_DEVICE( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
#endif

    IF (PRESENT(p_count)) THEN
       CALL MPI_RECV (t_buffer, p_count, p_int, p_source, p_tag, &
            p_comm, p_status, p_error)
    ELSE
       CALL MPI_RECV (t_buffer, SIZE(t_buffer), p_int, p_source, p_tag, &
            p_comm, p_status, p_error)
    END IF

#ifdef __USE_G2G
!$ACC END HOST_DATA
!$ACC END DATA
#endif

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_RECV on ', my_process_mpi_all_id, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL abort_mpi
    END IF
#endif
#endif

  END SUBROUTINE p_recv_int_4d


  SUBROUTINE p_recv_bool (t_buffer, p_source, p_tag, p_count, comm)

    LOGICAL, INTENT(out) :: t_buffer
    INTEGER, INTENT(in)  :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

#ifdef __USE_G2G
!$ACC DATA PRESENT( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
!$ACC HOST_DATA USE_DEVICE( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
#endif

    IF (PRESENT(p_count)) THEN
       CALL MPI_RECV (t_buffer, p_count, p_bool, p_source, p_tag, &
            p_comm, p_status, p_error)
    ELSE
       CALL MPI_RECV (t_buffer, 1, p_bool, p_source, p_tag, &
            p_comm, p_status, p_error)
    END IF

#ifdef __USE_G2G
!$ACC END HOST_DATA
!$ACC END DATA
#endif

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_RECV on ', my_process_mpi_all_id, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL abort_mpi
    END IF
#endif
#endif

  END SUBROUTINE p_recv_bool

  SUBROUTINE p_recv_bool_1d (t_buffer, p_source, p_tag, p_count, comm)

    LOGICAL, INTENT(out) :: t_buffer(:)
    INTEGER, INTENT(in)  :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

#ifdef __USE_G2G
!$ACC DATA PRESENT( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
!$ACC HOST_DATA USE_DEVICE( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
#endif

    IF (PRESENT(p_count)) THEN
       CALL MPI_RECV (t_buffer, p_count, p_bool, p_source, p_tag, &
            p_comm, p_status, p_error)
    ELSE
       CALL MPI_RECV (t_buffer, SIZE(t_buffer), p_bool, p_source, p_tag, &
            p_comm, p_status, p_error)
    END IF

#ifdef __USE_G2G
!$ACC END HOST_DATA
!$ACC END DATA
#endif

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_RECV on ', my_process_mpi_all_id, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL abort_mpi
    END IF
#endif
#endif

  END SUBROUTINE p_recv_bool_1d


  SUBROUTINE p_recv_char_1d (t_buffer, p_source, p_tag, p_count, comm)

    CHARACTER(LEN=*), INTENT(out) :: t_buffer(:)
    INTEGER, INTENT(in)  :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

#ifdef __USE_G2G
!$ACC DATA PRESENT( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
!$ACC HOST_DATA USE_DEVICE( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
#endif

    IF (PRESENT(p_count)) THEN
       CALL MPI_RECV (t_buffer, LEN(t_buffer(1))*p_count, p_char, p_source, p_tag, &
            p_comm, p_status, p_error)
    ELSE
       CALL MPI_RECV (t_buffer, LEN(t_buffer(1))*SIZE(t_buffer), p_char, p_source, p_tag, &
            p_comm, p_status, p_error)
    END IF

#ifdef __USE_G2G
!$ACC END HOST_DATA
!$ACC END DATA
#endif

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_RECV on ', my_process_mpi_all_id, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL abort_mpi
    END IF
#endif
#endif

  END SUBROUTINE p_recv_char_1d


  SUBROUTINE p_recv_bool_2d (t_buffer, p_source, p_tag, p_count, comm)

    LOGICAL, INTENT(out) :: t_buffer(:,:)
    INTEGER, INTENT(in)  :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

#ifdef __USE_G2G
!$ACC DATA PRESENT( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
!$ACC HOST_DATA USE_DEVICE( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
#endif

    IF (PRESENT(p_count)) THEN
       CALL MPI_RECV (t_buffer, p_count, p_bool, p_source, p_tag, &
            p_comm, p_status, p_error)
    ELSE
       CALL MPI_RECV (t_buffer, SIZE(t_buffer), p_bool, p_source, p_tag, &
            p_comm, p_status, p_error)
    END IF

#ifdef __USE_G2G
!$ACC END HOST_DATA
!$ACC END DATA
#endif

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_RECV on ', my_process_mpi_all_id, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL abort_mpi
    END IF
#endif
#endif

  END SUBROUTINE p_recv_bool_2d

  SUBROUTINE p_recv_bool_3d (t_buffer, p_source, p_tag, p_count, comm)

    LOGICAL, INTENT(out) :: t_buffer(:,:,:)
    INTEGER, INTENT(in)  :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

#ifdef __USE_G2G
!$ACC DATA PRESENT( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
!$ACC HOST_DATA USE_DEVICE( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
#endif

    IF (PRESENT(p_count)) THEN
       CALL MPI_RECV (t_buffer, p_count, p_bool, p_source, p_tag, &
            p_comm, p_status, p_error)
    ELSE
       CALL MPI_RECV (t_buffer, SIZE(t_buffer), p_bool, p_source, p_tag, &
            p_comm, p_status, p_error)
    END IF

#ifdef __USE_G2G
!$ACC END HOST_DATA
!$ACC END DATA
#endif

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_RECV on ', my_process_mpi_all_id, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL abort_mpi
    END IF
#endif
#endif

  END SUBROUTINE p_recv_bool_3d

  SUBROUTINE p_recv_bool_4d (t_buffer, p_source, p_tag, p_count, comm)

    LOGICAL, INTENT(out) :: t_buffer(:,:,:,:)
    INTEGER, INTENT(in)  :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

#ifdef __USE_G2G
!$ACC DATA PRESENT( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
!$ACC HOST_DATA USE_DEVICE( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
#endif

    IF (PRESENT(p_count)) THEN
       CALL MPI_RECV (t_buffer, p_count, p_bool, p_source, p_tag, &
            p_comm, p_status, p_error)
    ELSE
       CALL MPI_RECV (t_buffer, SIZE(t_buffer), p_bool, p_source, p_tag, &
            p_comm, p_status, p_error)
    END IF

#ifdef __USE_G2G
!$ACC END HOST_DATA
!$ACC END DATA
#endif

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_RECV on ', my_process_mpi_all_id, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL abort_mpi
    END IF
#endif
#endif

  END SUBROUTINE p_recv_bool_4d

  SUBROUTINE p_recv_char (t_buffer, p_source, p_tag, p_count, comm)

    CHARACTER(len=*),  INTENT(out) :: t_buffer
    INTEGER,           INTENT(in)  :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
#ifdef DEBUG
    CHARACTER(LEN=*), PARAMETER :: routine = 'p_recv_char'
#endif
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

#ifdef __USE_G2G
!$ACC DATA PRESENT( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
!$ACC HOST_DATA USE_DEVICE( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
#endif

    IF (PRESENT(p_count)) THEN
       CALL MPI_RECV (t_buffer, p_count, p_char, p_source, p_tag, &
            p_comm, p_status, p_error)
    ELSE
       CALL MPI_RECV (t_buffer, LEN(t_buffer), p_char, p_source, p_tag, &
            p_comm, p_status, p_error)
    END IF

#ifdef __USE_G2G
!$ACC END HOST_DATA
!$ACC END DATA
#endif

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) CALL finish (routine, 'MPI function call failed')
#endif
#endif

  END SUBROUTINE p_recv_char

  ! non-blocking receives

  !================================================================================================
  ! CHARACTER SECTION -----------------------------------------------------------------------------
  !
  SUBROUTINE p_irecv_char (t_buffer, p_source, p_tag, p_count, comm)

    CHARACTER(len=*), INTENT(inout) :: t_buffer
    INTEGER,   INTENT(in) :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

#ifdef __USE_G2G
!$ACC DATA PRESENT( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
!$ACC HOST_DATA USE_DEVICE( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
#endif

    IF (PRESENT(p_count)) THEN
       CALL p_inc_request
       CALL MPI_IRECV (t_buffer, p_count, p_char, p_source, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    ELSE
       CALL p_inc_request
       CALL MPI_IRECV (t_buffer, 1, p_char, p_source, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    END IF

#ifdef __USE_G2G
!$ACC END HOST_DATA
!$ACC END DATA
#endif

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_IRECV on ', my_process_mpi_all_id, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL abort_mpi
    END IF
#endif
#endif

  END SUBROUTINE p_irecv_char
  !================================================================================================
  ! REAL SECTION ----------------------------------------------------------------------------------
  !
  SUBROUTINE p_irecv_real (t_buffer, p_source, p_tag, p_count, comm)

    REAL(dp),  INTENT(inout) :: t_buffer
    INTEGER,   INTENT(in) :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

#ifdef __USE_G2G
!$ACC DATA PRESENT( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
!$ACC HOST_DATA USE_DEVICE( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
#endif

    IF (PRESENT(p_count)) THEN
       CALL p_inc_request
       CALL MPI_IRECV (t_buffer, p_count, p_real_dp, p_source, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    ELSE
       CALL p_inc_request
       CALL MPI_IRECV (t_buffer, 1, p_real_dp, p_source, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    END IF

#ifdef __USE_G2G
!$ACC END HOST_DATA
!$ACC END DATA
#endif

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_IRECV on ', my_process_mpi_all_id, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL abort_mpi
    END IF
#endif
#endif

  END SUBROUTINE p_irecv_real


  SUBROUTINE p_irecv_sreal (t_buffer, p_source, p_tag, p_count, comm)

    REAL(sp),  INTENT(inout) :: t_buffer
    INTEGER,   INTENT(in) :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

#ifdef __USE_G2G
!$ACC DATA PRESENT( t_buffer ), IF ( i_am_accel_node )
!$ACC HOST_DATA USE_DEVICE( t_buffer ), IF ( i_am_accel_node )
#endif

    IF (PRESENT(p_count)) THEN
       CALL p_inc_request
       CALL MPI_IRECV (t_buffer, p_count, p_real_sp, p_source, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    ELSE
       CALL p_inc_request
       CALL MPI_IRECV (t_buffer, 1, p_real_sp, p_source, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    END IF

#ifdef __USE_G2G
!$ACC END HOST_DATA
!$ACC END DATA
#endif

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_IRECV on ', my_process_mpi_all_id, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL abort_mpi
    END IF
#endif
#endif

  END SUBROUTINE p_irecv_sreal


  SUBROUTINE p_irecv_real_1d (t_buffer, p_source, p_tag, p_count, comm)

    REAL(dp),  INTENT(inout) :: t_buffer(:)
    INTEGER,   INTENT(in) :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

#ifdef __USE_G2G
!$ACC DATA PRESENT( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
!$ACC HOST_DATA USE_DEVICE( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
#endif

    IF (PRESENT(p_count)) THEN
       CALL p_inc_request
       CALL MPI_IRECV (t_buffer, p_count, p_real_dp, p_source, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    ELSE
       CALL p_inc_request
       CALL MPI_IRECV (t_buffer, SIZE(t_buffer), p_real_dp, p_source, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    END IF

#ifdef __USE_G2G
!$ACC END HOST_DATA
!$ACC END DATA
#endif

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_IRECV on ', my_process_mpi_all_id, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL abort_mpi
    END IF
#endif
#endif

  END SUBROUTINE p_irecv_real_1d

  SUBROUTINE p_irecv_real_2d (t_buffer, p_source, p_tag, p_count, comm)

    REAL(dp),  INTENT(inout) :: t_buffer(:,:)
    INTEGER,   INTENT(in) :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

#ifdef __USE_G2G
!$ACC DATA PRESENT( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
!$ACC HOST_DATA USE_DEVICE( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
#endif

    IF (PRESENT(p_count)) THEN
       CALL p_inc_request
       CALL MPI_IRECV (t_buffer, p_count, p_real_dp, p_source, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    ELSE
       CALL p_inc_request
       CALL MPI_IRECV (t_buffer, SIZE(t_buffer), p_real_dp, p_source, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    END IF

#ifdef __USE_G2G
!$ACC END HOST_DATA
!$ACC END DATA
#endif

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_IRECV on ', my_process_mpi_all_id, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL abort_mpi
    END IF
#endif
#endif

  END SUBROUTINE p_irecv_real_2d


  SUBROUTINE p_irecv_sreal_2d (t_buffer, p_source, p_tag, p_count, comm)

    REAL(sp),  INTENT(inout) :: t_buffer(:,:)
    INTEGER,   INTENT(in) :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

#ifdef __USE_G2G
!$ACC DATA PRESENT( t_buffer ), IF ( i_am_accel_node )
!$ACC HOST_DATA USE_DEVICE( t_buffer ), IF ( i_am_accel_node )
#endif

    IF (PRESENT(p_count)) THEN
       CALL p_inc_request
       CALL MPI_IRECV (t_buffer, p_count, p_real_sp, p_source, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    ELSE
       CALL p_inc_request
       CALL MPI_IRECV (t_buffer, SIZE(t_buffer), p_real_sp, p_source, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    END IF

#ifdef __USE_G2G
!$ACC END HOST_DATA
!$ACC END DATA
#endif

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_IRECV on ', my_process_mpi_all_id, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL abort_mpi
    END IF
#endif
#endif

  END SUBROUTINE p_irecv_sreal_2d


  SUBROUTINE p_irecv_real_3d (t_buffer, p_source, p_tag, p_count, comm)

    REAL(dp),  INTENT(inout) :: t_buffer(:,:,:)
    INTEGER,   INTENT(in) :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

#ifdef __USE_G2G
!$ACC DATA PRESENT( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
!$ACC HOST_DATA USE_DEVICE( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
#endif

    IF (PRESENT(p_count)) THEN
       CALL p_inc_request
       CALL MPI_IRECV (t_buffer, p_count, p_real_dp, p_source, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    ELSE
       CALL p_inc_request
       CALL MPI_IRECV (t_buffer, SIZE(t_buffer), p_real_dp, p_source, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    END IF

#ifdef __USE_G2G
!$ACC END HOST_DATA
!$ACC END DATA
#endif

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_IRECV on ', my_process_mpi_all_id, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL abort_mpi
    END IF
#endif
#endif

  END SUBROUTINE p_irecv_real_3d

  SUBROUTINE p_irecv_real_4d (t_buffer, p_source, p_tag, p_count, comm)

    REAL(dp),  INTENT(inout) :: t_buffer(:,:,:,:)
    INTEGER,   INTENT(in) :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

#ifdef __USE_G2G
!$ACC DATA PRESENT( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
!$ACC HOST_DATA USE_DEVICE( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
#endif

    IF (PRESENT(p_count)) THEN
       CALL p_inc_request
       CALL MPI_IRECV (t_buffer, p_count, p_real_dp, p_source, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    ELSE
       CALL p_inc_request
       CALL MPI_IRECV (t_buffer, SIZE(t_buffer), p_real_dp, p_source, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    END IF

#ifdef __USE_G2G
!$ACC END HOST_DATA
!$ACC END DATA
#endif

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_IRECV on ', my_process_mpi_all_id, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL abort_mpi
    END IF
#endif
#endif

  END SUBROUTINE p_irecv_real_4d
  !================================================================================================
  ! INTEGER SECTION -------------------------------------------------------------------------------
  !
  SUBROUTINE p_irecv_int (t_buffer, p_source, p_tag, p_count, comm)

    INTEGER, INTENT(inout) :: t_buffer
    INTEGER,   INTENT(in) :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

#ifdef __USE_G2G
!$ACC DATA PRESENT( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
!$ACC HOST_DATA USE_DEVICE( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
#endif

    IF (PRESENT(p_count)) THEN
       CALL p_inc_request
       CALL MPI_IRECV (t_buffer, p_count, p_int, p_source, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    ELSE
       CALL p_inc_request
       CALL MPI_IRECV (t_buffer, 1, p_int, p_source, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    END IF

#ifdef __USE_G2G
!$ACC END HOST_DATA
!$ACC END DATA
#endif

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_IRECV on ', my_process_mpi_all_id, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL abort_mpi
    END IF
#endif
#endif

  END SUBROUTINE p_irecv_int

  SUBROUTINE p_irecv_int_1d (t_buffer, p_source, p_tag, p_count, comm)

    INTEGER, INTENT(inout) :: t_buffer(:)
    INTEGER,   INTENT(in) :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

#ifdef __USE_G2G
!$ACC DATA PRESENT( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
!$ACC HOST_DATA USE_DEVICE( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
#endif

    IF (PRESENT(p_count)) THEN
       CALL p_inc_request
       CALL MPI_IRECV (t_buffer, p_count, p_int, p_source, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    ELSE
       CALL p_inc_request
       CALL MPI_IRECV (t_buffer, SIZE(t_buffer), p_int, p_source, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    END IF

#ifdef __USE_G2G
!$ACC END HOST_DATA
!$ACC END DATA
#endif

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_IRECV on ', my_process_mpi_all_id, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL abort_mpi
    END IF
#endif
#endif

  END SUBROUTINE p_irecv_int_1d

  SUBROUTINE p_irecv_int_2d (t_buffer, p_source, p_tag, p_count, comm)

    INTEGER, INTENT(inout) :: t_buffer(:,:)
    INTEGER,   INTENT(in) :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

#ifdef __USE_G2G
!$ACC DATA PRESENT( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
!$ACC HOST_DATA USE_DEVICE( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
#endif

    IF (PRESENT(p_count)) THEN
       CALL p_inc_request
       CALL MPI_IRECV (t_buffer, p_count, p_int, p_source, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    ELSE
       CALL p_inc_request
       CALL MPI_IRECV (t_buffer, SIZE(t_buffer), p_int, p_source, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    END IF

#ifdef __USE_G2G
!$ACC END HOST_DATA
!$ACC END DATA
#endif

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_IRECV on ', my_process_mpi_all_id, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL abort_mpi
    END IF
#endif
#endif

  END SUBROUTINE p_irecv_int_2d

  SUBROUTINE p_irecv_int_3d (t_buffer, p_source, p_tag, p_count, comm)

    INTEGER, INTENT(inout) :: t_buffer(:,:,:)
    INTEGER,   INTENT(in) :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

#ifdef __USE_G2G
!$ACC DATA PRESENT( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
!$ACC HOST_DATA USE_DEVICE( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
#endif

    IF (PRESENT(p_count)) THEN
       CALL p_inc_request
       CALL MPI_IRECV (t_buffer, p_count, p_int, p_source, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    ELSE
       CALL p_inc_request
       CALL MPI_IRECV (t_buffer, SIZE(t_buffer), p_int, p_source, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    END IF

#ifdef __USE_G2G
!$ACC END HOST_DATA
!$ACC END DATA
#endif

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_IRECV on ', my_process_mpi_all_id, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL abort_mpi
    END IF
#endif
#endif

  END SUBROUTINE p_irecv_int_3d

  SUBROUTINE p_irecv_int_4d (t_buffer, p_source, p_tag, p_count, comm)

    INTEGER, INTENT(inout) :: t_buffer(:,:,:,:)
    INTEGER,   INTENT(in) :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

#ifdef __USE_G2G
!$ACC DATA PRESENT( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
!$ACC HOST_DATA USE_DEVICE( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
#endif

    IF (PRESENT(p_count)) THEN
       CALL p_inc_request
       CALL MPI_IRECV (t_buffer, p_count, p_int, p_source, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    ELSE
       CALL p_inc_request
       CALL MPI_IRECV (t_buffer, SIZE(t_buffer), p_int, p_source, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    END IF

#ifdef __USE_G2G
!$ACC END HOST_DATA
!$ACC END DATA
#endif

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_IRECV on ', my_process_mpi_all_id, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL abort_mpi
    END IF
#endif
#endif

  END SUBROUTINE p_irecv_int_4d

  !================================================================================================
  ! LOGICAL SECTION -------------------------------------------------------------------------------
  !
  SUBROUTINE p_irecv_bool (t_buffer, p_source, p_tag, p_count, comm)

    LOGICAL, INTENT(inout) :: t_buffer
    INTEGER,   INTENT(in) :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL p_inc_request
       CALL MPI_IRECV (t_buffer, p_count, p_bool, p_source, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    ELSE
       CALL p_inc_request
       CALL MPI_IRECV (t_buffer, 1, p_bool, p_source, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    END IF

#ifdef __USE_G2G
!$ACC END HOST_DATA
!$ACC END DATA
#endif

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_IRECV on ', my_process_mpi_all_id, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL abort_mpi
    END IF
#endif
#endif

  END SUBROUTINE p_irecv_bool

  SUBROUTINE p_irecv_bool_1d (t_buffer, p_source, p_tag, p_count, comm)

    LOGICAL, INTENT(inout) :: t_buffer(:)
    INTEGER,   INTENT(in) :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

#ifdef __USE_G2G
!$ACC DATA PRESENT( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
!$ACC HOST_DATA USE_DEVICE( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
#endif

    IF (PRESENT(p_count)) THEN
       CALL p_inc_request
       CALL MPI_IRECV (t_buffer, p_count, p_bool, p_source, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    ELSE
       CALL p_inc_request
       CALL MPI_IRECV (t_buffer, SIZE(t_buffer), p_bool, p_source, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    END IF

#ifdef __USE_G2G
!$ACC END HOST_DATA
!$ACC END DATA
#endif

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_IRECV on ', my_process_mpi_all_id, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL abort_mpi
    END IF
#endif
#endif

  END SUBROUTINE p_irecv_bool_1d

  SUBROUTINE p_irecv_bool_2d (t_buffer, p_source, p_tag, p_count, comm)

    LOGICAL, INTENT(inout) :: t_buffer(:,:)
    INTEGER,   INTENT(in) :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

#ifdef __USE_G2G
!$ACC DATA PRESENT( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
!$ACC HOST_DATA USE_DEVICE( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
#endif

    IF (PRESENT(p_count)) THEN
       CALL p_inc_request
       CALL MPI_IRECV (t_buffer, p_count, p_bool, p_source, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    ELSE
       CALL p_inc_request
       CALL MPI_IRECV (t_buffer, SIZE(t_buffer), p_bool, p_source, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    END IF

#ifdef __USE_G2G
!$ACC END HOST_DATA
!$ACC END DATA
#endif

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_IRECV on ', my_process_mpi_all_id, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL abort_mpi
    END IF
#endif
#endif

  END SUBROUTINE p_irecv_bool_2d

  SUBROUTINE p_irecv_bool_3d (t_buffer, p_source, p_tag, p_count, comm)

    LOGICAL, INTENT(inout) :: t_buffer(:,:,:)
    INTEGER,   INTENT(in) :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

#ifdef __USE_G2G
!$ACC DATA PRESENT( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
!$ACC HOST_DATA USE_DEVICE( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
#endif

    IF (PRESENT(p_count)) THEN
       CALL p_inc_request
       CALL MPI_IRECV (t_buffer, p_count, p_bool, p_source, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    ELSE
       CALL p_inc_request
       CALL MPI_IRECV (t_buffer, SIZE(t_buffer), p_bool, p_source, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    END IF

#ifdef __USE_G2G
!$ACC END HOST_DATA
!$ACC END DATA
#endif

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_IRECV on ', my_process_mpi_all_id, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL abort_mpi
    END IF
#endif
#endif

  END SUBROUTINE p_irecv_bool_3d

  SUBROUTINE p_irecv_bool_4d (t_buffer, p_source, p_tag, p_count, comm)

    LOGICAL, INTENT(inout) :: t_buffer(:,:,:,:)
    INTEGER,   INTENT(in) :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

#ifdef __USE_G2G
!$ACC DATA PRESENT( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
!$ACC HOST_DATA USE_DEVICE( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
#endif

    IF (PRESENT(p_count)) THEN
       CALL p_inc_request
       CALL MPI_IRECV (t_buffer, p_count, p_bool, p_source, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    ELSE
       CALL p_inc_request
       CALL MPI_IRECV (t_buffer, SIZE(t_buffer), p_bool, p_source, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    END IF

#ifdef __USE_G2G
!$ACC END HOST_DATA
!$ACC END DATA
#endif

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_IRECV on ', my_process_mpi_all_id, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL abort_mpi
    END IF
#endif
#endif

  END SUBROUTINE p_irecv_bool_4d

  SUBROUTINE p_pack_int (t_var, t_buffer, p_pos, comm)

    INTEGER,   INTENT(IN)    :: t_var
    CHARACTER, INTENT(INOUT) :: t_buffer(:)
    INTEGER,   INTENT(INOUT) :: p_pos
    INTEGER, OPTIONAL, INTENT(IN) :: comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

    CALL MPI_PACK(t_var, 1, p_int, t_buffer, SIZE(t_buffer), p_pos, p_comm, p_error)
#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) CALL finish ("p_pack_int", 'MPI call failed')
#endif
#endif
  END SUBROUTINE p_pack_int

  SUBROUTINE p_pack_bool (t_var, t_buffer, p_pos, comm)

    LOGICAL,   INTENT(IN)    :: t_var
    CHARACTER, INTENT(INOUT) :: t_buffer(:)
    INTEGER,   INTENT(INOUT) :: p_pos
    INTEGER, OPTIONAL, INTENT(IN) :: comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

    CALL MPI_PACK(t_var, 1, p_bool, t_buffer, SIZE(t_buffer), p_pos, p_comm, p_error)
#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) CALL finish ("p_pack_int", 'MPI call failed')
#endif
#endif
  END SUBROUTINE p_pack_bool

  SUBROUTINE p_pack_real (t_var, t_buffer, p_pos, comm)

    REAL(wp),  INTENT(IN)    :: t_var
    CHARACTER, INTENT(INOUT) :: t_buffer(:)
    INTEGER,   INTENT(INOUT) :: p_pos
    INTEGER, OPTIONAL, INTENT(IN) :: comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

    CALL MPI_PACK(t_var, 1, p_real_dp, t_buffer, SIZE(t_buffer), p_pos, p_comm, p_error)
#ifdef DEBUG
   IF (p_error /= MPI_SUCCESS) CALL finish ("p_pack_real", 'MPI call failed')
#endif
#endif
  END SUBROUTINE p_pack_real

  SUBROUTINE p_pack_int_1d (t_var, p_count, t_buffer, p_pos, comm)

    INTEGER,   INTENT(IN)    :: t_var(:)
    INTEGER,   INTENT(IN)    :: p_count
    CHARACTER, INTENT(INOUT) :: t_buffer(:)
    INTEGER,   INTENT(INOUT) :: p_pos
    INTEGER, OPTIONAL, INTENT(IN) :: comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

    CALL MPI_PACK(t_var, p_count, p_int, t_buffer, SIZE(t_buffer), p_pos, p_comm, p_error)
#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) CALL finish ("p_pack_int_1d", 'MPI call failed')
#endif
#endif
  END SUBROUTINE p_pack_int_1d

  SUBROUTINE p_pack_real_1d (t_var, p_count, t_buffer, p_pos, comm)

    REAL(wp),  INTENT(IN)    :: t_var(:)
    INTEGER,   INTENT(IN)    :: p_count
    CHARACTER, INTENT(INOUT) :: t_buffer(:)
    INTEGER,   INTENT(INOUT) :: p_pos
    INTEGER, OPTIONAL, INTENT(IN) :: comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

    CALL MPI_PACK(t_var, p_count, p_real_dp, t_buffer, SIZE(t_buffer), p_pos, p_comm, p_error)
#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) CALL finish ("p_pack_real_1d", 'MPI call failed')
#endif
#endif
  END SUBROUTINE p_pack_real_1d


  SUBROUTINE p_pack_string (t_var, t_buffer, p_pos, comm)

    CHARACTER(LEN=*), INTENT(IN) :: t_var
    CHARACTER, INTENT(INOUT) :: t_buffer(:)
    INTEGER,   INTENT(INOUT) :: p_pos
    INTEGER, OPTIONAL, INTENT(IN) :: comm
#ifndef NOMPI
    INTEGER :: p_comm, ilength
    CHARACTER(LEN=128) :: tmp_var

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

    tmp_var =     TRIM(t_var)
    ilength = LEN_TRIM(t_var)
    ! first pack string length, then the character sequence
    CALL MPI_PACK(ilength,       1,  p_int, t_buffer, SIZE(t_buffer), p_pos, p_comm, p_error)
#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) CALL finish ("p_pack_char_1d", 'MPI call failed')
#endif
    CALL MPI_PACK(tmp_var, ilength, p_char, t_buffer, SIZE(t_buffer), p_pos, p_comm, p_error)
#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) CALL finish ("p_pack_char_1d", 'MPI call failed')
#endif
#endif
  END SUBROUTINE p_pack_string

  SUBROUTINE p_pack_real_2d (t_var, p_count, t_buffer, p_pos, comm)

    REAL(wp),  INTENT(IN)    :: t_var(:,:)
    INTEGER,   INTENT(IN)    :: p_count
    CHARACTER, INTENT(INOUT) :: t_buffer(:)
    INTEGER,   INTENT(INOUT) :: p_pos
    INTEGER, OPTIONAL, INTENT(IN) :: comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

    CALL MPI_PACK(t_var, p_count, p_real_dp, t_buffer, SIZE(t_buffer), p_pos, p_comm, p_error)
#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) CALL finish ("p_pack_real_2d", 'MPI call failed')
#endif
#endif
  END SUBROUTINE p_pack_real_2d

  SUBROUTINE p_unpack_int (t_buffer, p_pos, t_var, comm)

    CHARACTER, INTENT(IN)    :: t_buffer(:)
    INTEGER,   INTENT(INOUT) :: p_pos
    INTEGER,   INTENT(OUT)   :: t_var
    INTEGER, OPTIONAL, INTENT(IN) :: comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

    CALL MPI_UNPACK(t_buffer, SIZE(t_buffer), p_pos, t_var, 1, p_int, p_comm, p_error)
#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) CALL finish ("p_unpack_int", 'MPI call failed')
#endif
#endif
  END SUBROUTINE p_unpack_int

  SUBROUTINE p_unpack_bool (t_buffer, p_pos, t_var, comm)

    CHARACTER, INTENT(IN)    :: t_buffer(:)
    INTEGER,   INTENT(INOUT) :: p_pos
    LOGICAL,   INTENT(OUT)   :: t_var
    INTEGER, OPTIONAL, INTENT(IN) :: comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

    CALL MPI_UNPACK(t_buffer, SIZE(t_buffer), p_pos, t_var, 1, p_bool, p_comm, p_error)
#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) CALL finish ("p_unpack_int", 'MPI call failed')
#endif
#endif
  END SUBROUTINE p_unpack_bool

  SUBROUTINE p_unpack_real (t_buffer, p_pos, t_var, comm)

    CHARACTER, INTENT(IN)    :: t_buffer(:)
    INTEGER,   INTENT(INOUT) :: p_pos
    REAL(wp),  INTENT(OUT)   :: t_var
    INTEGER, OPTIONAL, INTENT(IN) :: comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

    CALL MPI_UNPACK(t_buffer, SIZE(t_buffer), p_pos, t_var, 1, p_real_dp, p_comm, p_error)
#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) CALL finish ("p_unpack_real", 'MPI call failed')
#endif
#endif
  END SUBROUTINE p_unpack_real

  SUBROUTINE p_unpack_int_1d (t_buffer, p_pos, t_var, p_count, comm)

    CHARACTER, INTENT(IN)    :: t_buffer(:)
    INTEGER,   INTENT(INOUT) :: p_pos
    INTEGER,   INTENT(INOUT) :: t_var(:)
    INTEGER,   INTENT(IN)    :: p_count
    INTEGER, OPTIONAL, INTENT(IN) :: comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

    CALL MPI_UNPACK(t_buffer, SIZE(t_buffer), p_pos, t_var, p_count, p_int, p_comm, p_error)
#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) CALL finish ("p_unpack_int_1d", 'MPI call failed')
#endif
#endif
  END SUBROUTINE p_unpack_int_1d

  SUBROUTINE p_unpack_real_1d (t_buffer, p_pos, t_var, p_count, comm)

    CHARACTER, INTENT(IN)    :: t_buffer(:)
    INTEGER,   INTENT(INOUT) :: p_pos
    REAL(wp),  INTENT(INOUT) :: t_var(:)
    INTEGER,   INTENT(IN)    :: p_count
    INTEGER, OPTIONAL, INTENT(IN) :: comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

    CALL MPI_UNPACK(t_buffer, SIZE(t_buffer), p_pos, t_var, p_count, p_real_dp, p_comm, p_error)
#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) CALL finish ("p_unpack_real_1d", 'MPI call failed')
#endif
#endif
  END SUBROUTINE p_unpack_real_1d


  SUBROUTINE p_unpack_string (t_buffer, p_pos, t_var, comm)

    CHARACTER, INTENT(IN) :: t_buffer(:)
    INTEGER,   INTENT(INOUT) :: p_pos
    CHARACTER(LEN=*), INTENT(INOUT) :: t_var
    INTEGER, OPTIONAL, INTENT(IN) :: comm
#ifndef NOMPI
    INTEGER :: p_comm, ilength

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF
    ! first unpack string length, then the character sequence
    CALL MPI_UNPACK(t_buffer, SIZE(t_buffer), p_pos, ilength,       1,  p_int, p_comm, p_error)
#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) CALL finish ("p_unpack_char_1d", 'MPI call failed')
#endif
    t_var = " "
    CALL MPI_UNPACK(t_buffer, SIZE(t_buffer), p_pos,   t_var, ilength, p_char, p_comm, p_error)
#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) CALL finish ("p_unpack_char_1d", 'MPI call failed')
#endif
#endif
  END SUBROUTINE p_unpack_string


  SUBROUTINE p_unpack_real_2d (t_buffer, p_pos, t_var, p_count, comm)

    CHARACTER, INTENT(IN)    :: t_buffer(:)
    INTEGER,   INTENT(INOUT) :: p_pos
    REAL(wp),  INTENT(INOUT) :: t_var(:,:)
    INTEGER,   INTENT(IN)    :: p_count
    INTEGER, OPTIONAL, INTENT(IN) :: comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

    CALL MPI_UNPACK(t_buffer, SIZE(t_buffer), p_pos, t_var, p_count, p_real_dp, p_comm, p_error)
#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) CALL finish ("p_unpack_real_2d", 'MPI call failed')
#endif
#endif
  END SUBROUTINE p_unpack_real_2d

  SUBROUTINE p_recv_packed (t_buffer, p_source, p_tag, p_count, comm)

    CHARACTER, INTENT(INOUT) :: t_buffer(:)
    INTEGER,   INTENT(IN)    :: p_source, p_tag
    INTEGER,   INTENT(IN)    :: p_count
    INTEGER, OPTIONAL, INTENT(IN) :: comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

    CALL MPI_RECV (t_buffer, p_count, MPI_PACKED, p_source, p_tag, &
      p_comm, p_status, p_error)

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_RECV on ', my_process_mpi_all_id, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_recv_packed

  SUBROUTINE p_irecv_packed (t_buffer, p_source, p_tag, p_count, comm)

    CHARACTER, INTENT(INOUT) :: t_buffer(:)
    INTEGER,   INTENT(IN)    :: p_source, p_tag
    INTEGER,   INTENT(IN)    :: p_count
    INTEGER, OPTIONAL, INTENT(IN) :: comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

    CALL p_inc_request
    CALL MPI_IRECV (t_buffer, p_count, MPI_PACKED, p_source, p_tag, &
      p_comm, p_request(p_irequest), p_error)

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_IRECV on ', my_process_mpi_all_id, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL abort_mpi
    END IF
#endif
#endif

  END SUBROUTINE p_irecv_packed

  SUBROUTINE p_send_packed (t_buffer, p_destination, p_tag, p_count, comm)

    CHARACTER, INTENT(in) :: t_buffer(:)
    INTEGER,   INTENT(in) :: p_destination, p_tag
    INTEGER,   INTENT(in) :: p_count
    INTEGER, OPTIONAL, INTENT(in) :: comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF
    CALL MPI_SEND (t_buffer, p_count, MPI_PACKED, p_destination, p_tag, &
      p_comm, p_error)

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_SEND from ', my_process_mpi_all_id, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL abort_mpi
    END IF
#endif
#endif

  END SUBROUTINE p_send_packed

  SUBROUTINE p_bcast_packed (t_buffer, p_source, p_count, comm)
    CHARACTER, INTENT(inout) :: t_buffer(:)
    INTEGER,   INTENT(in)    :: p_source
    INTEGER,   INTENT(in)    :: p_count
    INTEGER, OPTIONAL, INTENT(in) :: comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

    CALL MPI_BCAST (t_buffer, p_count, MPI_PACKED, p_source, &
            p_comm, p_error)

#ifdef DEBUG
    WRITE (nerr,'(a,i4,a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' with broadcast successful.'

    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif

#endif
  END SUBROUTINE p_bcast_packed

  !
  !================================================================================================

  ! sendrecv implementation

  SUBROUTINE p_sendrecv_real_1d (sendbuf, p_dest, recvbuf, p_source, &
                                  p_tag, comm)

    REAL(dp), INTENT(in)           :: sendbuf (:)
    INTEGER,  INTENT(in)           :: p_dest
    REAL(dp), INTENT(out)          :: recvbuf (:)
    INTEGER,  INTENT(in)           :: p_source
    INTEGER,  INTENT(in)           :: p_tag
    INTEGER,  INTENT(in) ,OPTIONAL :: comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

#ifdef __USE_G2G
!$ACC DATA PRESENT( sendbuf, recvbuf ), IF ( i_am_accel_node .AND. acc_on )
!$ACC HOST_DATA USE_DEVICE( sendbuf, recvbuf ), IF ( i_am_accel_node .AND. acc_on )
#endif

    CALL MPI_SENDRECV (sendbuf, SIZE(sendbuf), p_real_dp, p_dest,   p_tag, &
                        recvbuf, SIZE(recvbuf), p_real_dp, p_source, p_tag, &
                        p_comm, p_status, p_error)

#ifdef __USE_G2G
!$ACC END HOST_DATA
!$ACC END DATA
#endif

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i4,a,i6,a)') ' MPI_SENDRECV by ', my_process_mpi_all_id, &
            ' to ', p_dest, ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL abort_mpi
    END IF
#endif
#endif

  END SUBROUTINE p_sendrecv_real_1d

  SUBROUTINE p_sendrecv_real_2d (sendbuf, p_dest, recvbuf, p_source, &
                                  p_tag, comm)

    REAL(dp), INTENT(in)           :: sendbuf (:,:)
    INTEGER,  INTENT(in)           :: p_dest
    REAL(dp), INTENT(out)          :: recvbuf (:,:)
    INTEGER,  INTENT(in)           :: p_source
    INTEGER,  INTENT(in)           :: p_tag
    INTEGER,  INTENT(in) ,OPTIONAL :: comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

#ifdef __USE_G2G
!$ACC DATA PRESENT( sendbuf, recvbuf ), IF ( i_am_accel_node .AND. acc_on )
!$ACC HOST_DATA USE_DEVICE( sendbuf, recvbuf ), IF ( i_am_accel_node .AND. acc_on )
#endif

    CALL MPI_SENDRECV (sendbuf, SIZE(sendbuf), p_real_dp, p_dest,   p_tag, &
                        recvbuf, SIZE(recvbuf), p_real_dp, p_source, p_tag, &
                        p_comm, p_status, p_error)

#ifdef __USE_G2G
!$ACC END HOST_DATA
!$ACC END DATA
#endif

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i4,a,i6,a)') ' MPI_SENDRECV by ', my_process_mpi_all_id, &
            ' to ', p_dest, ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL abort_mpi
    END IF
#endif
#endif

  END SUBROUTINE p_sendrecv_real_2d

  SUBROUTINE p_sendrecv_real_3d (sendbuf, p_dest, recvbuf, p_source, &
                                  p_tag, comm)

    REAL(dp), INTENT(in)           :: sendbuf (:,:,:)
    INTEGER,  INTENT(in)           :: p_dest
    REAL(dp), INTENT(out)          :: recvbuf (:,:,:)
    INTEGER,  INTENT(in)           :: p_source
    INTEGER,  INTENT(in)           :: p_tag
    INTEGER,  INTENT(in) ,OPTIONAL :: comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

#ifdef __USE_G2G
!$ACC DATA PRESENT( sendbuf, recvbuf ), IF ( i_am_accel_node .AND. acc_on )
!$ACC HOST_DATA USE_DEVICE( sendbuf, recvbuf ), IF ( i_am_accel_node .AND. acc_on )
#endif

    CALL MPI_SENDRECV (sendbuf, SIZE(sendbuf), p_real_dp, p_dest,   p_tag, &
                        recvbuf, SIZE(recvbuf), p_real_dp, p_source, p_tag, &
                        p_comm, p_status, p_error)

#ifdef __USE_G2G
!$ACC END HOST_DATA
!$ACC END DATA
#endif

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i4,a,i6,a)') ' MPI_SENDRECV by ', my_process_mpi_all_id, &
            ' to ', p_dest, ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL abort_mpi
    END IF
#endif
#endif

  END SUBROUTINE p_sendrecv_real_3d

  SUBROUTINE p_sendrecv_real_4d (sendbuf, p_dest, recvbuf, p_source, &
                                  p_tag, comm)

    REAL(dp), INTENT(in)           :: sendbuf (:,:,:,:)
    INTEGER,  INTENT(in)           :: p_dest
    REAL(dp), INTENT(out)          :: recvbuf (:,:,:,:)
    INTEGER,  INTENT(in)           :: p_source
    INTEGER,  INTENT(in)           :: p_tag
    INTEGER,  INTENT(in) ,OPTIONAL :: comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

#ifdef __USE_G2G
!$ACC DATA PRESENT( sendbuf, recvbuf ), IF ( i_am_accel_node .AND. acc_on )
!$ACC HOST_DATA USE_DEVICE( sendbuf, recvbuf ), IF ( i_am_accel_node .AND. acc_on )
#endif

    CALL MPI_SENDRECV (sendbuf, SIZE(sendbuf), p_real_dp, p_dest,   p_tag, &
                        recvbuf, SIZE(recvbuf), p_real_dp, p_source, p_tag, &
                        p_comm, p_status, p_error)

#ifdef __USE_G2G
!$ACC END HOST_DATA
!$ACC END DATA
#endif

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i4,a,i6,a)') ' MPI_SENDRECV by ', my_process_mpi_all_id, &
            ' to ', p_dest, ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL abort_mpi
    END IF
#endif
#endif

  END SUBROUTINE p_sendrecv_real_4d

  SUBROUTINE p_sendrecv_char_array (sendbuf, p_dest, recvbuf, p_source, p_tag, comm)
    CHARACTER(KIND = C_CHAR), INTENT(in) :: sendbuf (:)
    CHARACTER(KIND = C_CHAR), INTENT(out) :: recvbuf (:)
    INTEGER, VALUE :: p_dest, p_source, p_tag
    INTEGER, INTENT(in), OPTIONAL :: comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

#ifdef __USE_G2G
!$ACC DATA PRESENT( sendbuf, recvbuf ), IF ( i_am_accel_node .AND. acc_on )
!$ACC HOST_DATA USE_DEVICE( sendbuf, recvbuf ), IF ( i_am_accel_node .AND. acc_on )
#endif

    CALL MPI_Sendrecv(sendbuf(1), SIZE(sendbuf, 1), MPI_CHARACTER, p_dest,   p_tag, &
    &                 recvbuf(1), SIZE(recvbuf, 1), MPI_CHARACTER, p_source, p_tag, &
    &                 p_comm, p_status, p_error)

#ifdef __USE_G2G
!$ACC END HOST_DATA
!$ACC END DATA
#endif

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i4,a,i6,a)') ' MPI_SENDRECV by ', my_process_mpi_all_id, &
            ' to ', p_dest, ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL abort_mpi
    END IF
#endif
#endif

  END SUBROUTINE p_sendrecv_char_array

  ! bcast implementation

  SUBROUTINE p_bcast_real (t_buffer, p_source, comm)

    REAL (dp), INTENT(inout) :: t_buffer
    INTEGER,   INTENT(in)    :: p_source
    INTEGER, OPTIONAL, INTENT(in) :: comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

#ifdef DEBUG
    nbcast = nbcast+1
#endif

    IF (process_mpi_all_size == 1) THEN
       RETURN
    ELSE

#ifdef __USE_G2G
!$ACC DATA PRESENT( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
!$ACC HOST_DATA USE_DEVICE( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
#endif

       CALL MPI_BCAST (t_buffer, 1, p_real_dp, p_source, &
            p_comm, p_error)

#ifdef __USE_G2G
!$ACC END HOST_DATA
!$ACC END DATA
#endif

    ENDIF

#ifdef DEBUG
    WRITE (nerr,'(a,i4,a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' with broadcast number ', nbcast, ' successful.'

    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL abort_mpi
    END IF
#endif
#endif

  END SUBROUTINE p_bcast_real

  !---------------------------------------------------------------------------------------------------------------------------------
  !> wrapper for MPI_Bcast
  !---------------------------------------------------------------------------------------------------------------------------------
  SUBROUTINE p_bcast_real_single (t_buffer, p_source, comm)

    REAL (sp), INTENT(inout) :: t_buffer
    INTEGER,   INTENT(in)    :: p_source
    INTEGER, OPTIONAL, INTENT(in) :: comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

#ifdef DEBUG
    nbcast = nbcast+1
#endif

    IF (process_mpi_all_size == 1) THEN
       RETURN
    ELSE

#ifdef __USE_G2G
!$ACC DATA PRESENT( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
!$ACC HOST_DATA USE_DEVICE( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
#endif

       CALL MPI_BCAST (t_buffer, 1, p_real_sp, p_source, &
            p_comm, p_error)

#ifdef __USE_G2G
!$ACC END HOST_DATA
!$ACC END DATA
#endif

    ENDIF

#ifdef DEBUG
    WRITE (nerr,'(a,i4,a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' with broadcast number ', nbcast, ' successful.'

    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL abort_mpi
    END IF
#endif
#endif

  END SUBROUTINE p_bcast_real_single

  SUBROUTINE p_bcast_real_1d (t_buffer, p_source, comm)

    REAL (dp), INTENT(inout) :: t_buffer(:)
    INTEGER,   INTENT(in)    :: p_source
    INTEGER, OPTIONAL, INTENT(in) :: comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

#ifdef DEBUG
    nbcast = nbcast+1
#endif

    IF (process_mpi_all_size == 1) THEN
       RETURN
    ELSE

#ifdef __USE_G2G
!$ACC DATA PRESENT( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
!$ACC HOST_DATA USE_DEVICE( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
#endif

       CALL MPI_BCAST (t_buffer, SIZE(t_buffer), p_real_dp, p_source, &
            p_comm, p_error)

#ifdef __USE_G2G
!$ACC END HOST_DATA
!$ACC END DATA
#endif

    ENDIF

#ifdef DEBUG
    WRITE (nerr,'(a,i4,a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' with broadcast number ', nbcast, ' successful.'

     IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL abort_mpi
    END IF
#endif
#endif

  END SUBROUTINE p_bcast_real_1d

  !---------------------------------------------------------------------------------------------------------------------------------
  !> wrapper for MPI_Bcast
  !---------------------------------------------------------------------------------------------------------------------------------
  SUBROUTINE p_bcast_real_1d_single (t_buffer, p_source, comm)

    REAL (sp), INTENT(inout) :: t_buffer(:)
    INTEGER,   INTENT(in)    :: p_source
    INTEGER, OPTIONAL, INTENT(in) :: comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

#ifdef DEBUG
    nbcast = nbcast+1
#endif

    IF (process_mpi_all_size == 1) THEN
       RETURN
    ELSE

#ifdef __USE_G2G
!$ACC DATA PRESENT( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
!$ACC HOST_DATA USE_DEVICE( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
#endif

       CALL MPI_BCAST (t_buffer, SIZE(t_buffer), p_real_sp, p_source, &
            p_comm, p_error)

#ifdef __USE_G2G
!$ACC END HOST_DATA
!$ACC END DATA
#endif

    ENDIF

#ifdef DEBUG
    WRITE (nerr,'(a,i4,a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' with broadcast number ', nbcast, ' successful.'

     IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL abort_mpi
    END IF
#endif
#endif

  END SUBROUTINE p_bcast_real_1d_single

  SUBROUTINE p_bcast_real_2d (t_buffer, p_source, comm)

    REAL (dp), INTENT(inout) :: t_buffer(:,:)
    INTEGER,   INTENT(in)    :: p_source
    INTEGER, OPTIONAL, INTENT(in) :: comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

#ifdef DEBUG
    nbcast = nbcast+1
#endif

    IF (process_mpi_all_size == 1) THEN
       RETURN
    ELSE

#ifdef __USE_G2G
!$ACC DATA PRESENT( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
!$ACC HOST_DATA USE_DEVICE( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
#endif

       CALL MPI_BCAST (t_buffer, SIZE(t_buffer), p_real_dp, p_source, &
            p_comm, p_error)

#ifdef __USE_G2G
!$ACC END HOST_DATA
!$ACC END DATA
#endif

    ENDIF

#ifdef DEBUG
    WRITE (nerr,'(a,i4,a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' with broadcast number ', nbcast, ' successful.'

     IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL abort_mpi
    END IF
#endif
#endif

  END SUBROUTINE p_bcast_real_2d


  SUBROUTINE p_bcast_real_2d_single (t_buffer, p_source, comm)

    REAL (sp), INTENT(inout) :: t_buffer(:,:) ! SINGLE PRECISION
    INTEGER,   INTENT(in)    :: p_source
    INTEGER, OPTIONAL, INTENT(in) :: comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

#ifdef DEBUG
    nbcast = nbcast+1
#endif

    IF (process_mpi_all_size == 1) THEN
       RETURN
    ELSE

#ifdef __USE_G2G
!$ACC DATA PRESENT( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
!$ACC HOST_DATA USE_DEVICE( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
#endif

       CALL MPI_BCAST (t_buffer, SIZE(t_buffer), p_real_sp, p_source, &
            p_comm, p_error)

#ifdef __USE_G2G
!$ACC END HOST_DATA
!$ACC END DATA
#endif

    ENDIF

#ifdef DEBUG
    WRITE (nerr,'(a,i4,a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' with broadcast number ', nbcast, ' successful.'

     IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL abort_mpi
    END IF
#endif
#endif

  END SUBROUTINE p_bcast_real_2d_single


  SUBROUTINE p_bcast_real_3d (t_buffer, p_source, comm)

    REAL (dp), INTENT(inout) :: t_buffer(:,:,:)
    INTEGER,   INTENT(in)    :: p_source
    INTEGER, OPTIONAL, INTENT(in) :: comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

#ifdef DEBUG
    nbcast = nbcast+1
#endif

    IF (process_mpi_all_size == 1) THEN
       RETURN
    ELSE

#ifdef __USE_G2G
!$ACC DATA PRESENT( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
!$ACC HOST_DATA USE_DEVICE( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
#endif

       CALL MPI_BCAST (t_buffer, SIZE(t_buffer), p_real_dp, p_source, &
            p_comm, p_error)

#ifdef __USE_G2G
!$ACC END HOST_DATA
!$ACC END DATA
#endif

    ENDIF

#ifdef DEBUG
    WRITE (nerr,'(a,i4,a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' with broadcast number ', nbcast, ' successful.'

    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL abort_mpi
    END IF
#endif
#endif

  END SUBROUTINE p_bcast_real_3d

  SUBROUTINE p_bcast_real_4d (t_buffer, p_source, comm)

    REAL (dp), INTENT(inout) :: t_buffer(:,:,:,:)
    INTEGER,   INTENT(in)    :: p_source
    INTEGER, OPTIONAL, INTENT(in) :: comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

#ifdef DEBUG
    nbcast = nbcast+1
#endif

    IF (process_mpi_all_size == 1) THEN
       RETURN
    ELSE

#ifdef __USE_G2G
!$ACC DATA PRESENT( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
!$ACC HOST_DATA USE_DEVICE( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
#endif

       CALL MPI_BCAST (t_buffer, SIZE(t_buffer), p_real_dp, p_source, &
            p_comm, p_error)

#ifdef __USE_G2G
!$ACC END HOST_DATA
!$ACC END DATA
#endif

    ENDIF

#ifdef DEBUG
    WRITE (nerr,'(a,i4,a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' with broadcast number ', nbcast, ' successful.'

     IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL abort_mpi
    END IF
#endif
#endif

  END SUBROUTINE p_bcast_real_4d

  SUBROUTINE p_bcast_real_5d (t_buffer, p_source, comm)

    REAL (dp), INTENT(inout) :: t_buffer(:,:,:,:,:)
    INTEGER,   INTENT(in)    :: p_source
    INTEGER, OPTIONAL, INTENT(in) :: comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

#ifdef DEBUG
    nbcast = nbcast+1
#endif

    IF (process_mpi_all_size == 1) THEN
       RETURN
    ELSE

#ifdef __USE_G2G
!$ACC DATA PRESENT( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
!$ACC HOST_DATA USE_DEVICE( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
#endif

       CALL MPI_BCAST (t_buffer, SIZE(t_buffer), p_real_dp, p_source, &
            p_comm, p_error)

#ifdef __USE_G2G
!$ACC END HOST_DATA
!$ACC END DATA
#endif

    ENDIF

#ifdef DEBUG
    WRITE (nerr,'(a,i4,a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' with broadcast number ', nbcast, ' successful.'

     IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL abort_mpi
    END IF
#endif
#endif

  END SUBROUTINE p_bcast_real_5d

  SUBROUTINE p_bcast_real_7d (t_buffer, p_source, comm)

    REAL (dp), INTENT(inout) :: t_buffer(:,:,:,:,:,:,:)
    INTEGER,   INTENT(in)    :: p_source
    INTEGER, OPTIONAL, INTENT(in) :: comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

#ifdef DEBUG
    nbcast = nbcast+1
#endif

    IF (process_mpi_all_size == 1) THEN
       RETURN
    ELSE

#ifdef __USE_G2G
!$ACC DATA PRESENT( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
!$ACC HOST_DATA USE_DEVICE( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
#endif

       CALL MPI_BCAST (t_buffer, SIZE(t_buffer), p_real_dp, p_source, &
            p_comm, p_error)

#ifdef __USE_G2G
!$ACC END HOST_DATA
!$ACC END DATA
#endif

    ENDIF

#ifdef DEBUG
    WRITE (nerr,'(a,i4,a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' with broadcast number ', nbcast, ' successful.'

     IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL abort_mpi
    END IF
#endif
#endif

  END SUBROUTINE p_bcast_real_7d

  SUBROUTINE p_bcast_int_i4 (t_buffer, p_source, comm)

    INTEGER (i4), INTENT(inout) :: t_buffer
    INTEGER, INTENT(in)    :: p_source
    INTEGER, OPTIONAL, INTENT(in) :: comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

#ifdef DEBUG
    nbcast = nbcast+1
#endif

    IF (process_mpi_all_size == 1) THEN
       RETURN
    ELSE

#ifdef __USE_G2G
!$ACC DATA PRESENT( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
!$ACC HOST_DATA USE_DEVICE( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
#endif

       CALL MPI_BCAST (t_buffer, 1, p_int_i4, p_source, &
            p_comm, p_error)

#ifdef __USE_G2G
!$ACC END HOST_DATA
!$ACC END DATA
#endif

    ENDIF

#ifdef DEBUG
    WRITE (nerr,'(a,i4,a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' with broadcast number ', nbcast, ' successful.'

    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL abort_mpi
    END IF
#endif
#endif

  END SUBROUTINE p_bcast_int_i4

  SUBROUTINE p_bcast_int_i8 (t_buffer, p_source, comm)

    INTEGER (i8), INTENT(inout) :: t_buffer
    INTEGER, INTENT(in)    :: p_source
    INTEGER, OPTIONAL, INTENT(in) :: comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

#ifdef DEBUG
    nbcast = nbcast+1
#endif

    IF (process_mpi_all_size == 1) THEN
       RETURN
    ELSE

#ifdef __USE_G2G
!$ACC DATA PRESENT( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
!$ACC HOST_DATA USE_DEVICE( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
#endif

       CALL MPI_BCAST (t_buffer, 1, p_int_i8, p_source, &
            p_comm, p_error)

#ifdef __USE_G2G
!$ACC END HOST_DATA
!$ACC END DATA
#endif

    ENDIF

#ifdef DEBUG
    WRITE (nerr,'(a,i4,a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' with broadcast number ', nbcast, ' successful.'

     IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL abort_mpi
    END IF
#endif
#endif

  END SUBROUTINE p_bcast_int_i8

  SUBROUTINE p_bcast_int_1d (t_buffer, p_source, comm)

    INTEGER, INTENT(inout) :: t_buffer(:)
    INTEGER, INTENT(in)    :: p_source
    INTEGER, OPTIONAL, INTENT(in) :: comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

#ifdef DEBUG
    nbcast = nbcast+1
#endif

    IF (process_mpi_all_size == 1) THEN
       RETURN
    ELSE

#ifdef __USE_G2G
!$ACC DATA PRESENT( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
!$ACC HOST_DATA USE_DEVICE( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
#endif

       CALL MPI_BCAST (t_buffer, SIZE(t_buffer), p_int, p_source, &
            p_comm, p_error)

#ifdef __USE_G2G
!$ACC END HOST_DATA
!$ACC END DATA
#endif

    ENDIF

#ifdef DEBUG
    WRITE (nerr,'(a,i4,a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' with broadcast number ', nbcast, ' successful.'

     IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL abort_mpi
    END IF
#endif
#endif

  END SUBROUTINE p_bcast_int_1d

  SUBROUTINE p_bcast_int_i8_1d (t_buffer, p_source, comm)

    INTEGER(i8), INTENT(inout) :: t_buffer(:)
    INTEGER, INTENT(in)    :: p_source
    INTEGER, OPTIONAL, INTENT(in) :: comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

#ifdef DEBUG
    nbcast = nbcast+1
#endif

    IF (process_mpi_all_size == 1) THEN
       RETURN
    ELSE

#ifdef __USE_G2G
!$ACC DATA PRESENT( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
!$ACC HOST_DATA USE_DEVICE( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
#endif

       CALL MPI_BCAST (t_buffer, SIZE(t_buffer), p_int_i8, p_source, &
            p_comm, p_error)

#ifdef __USE_G2G
!$ACC END HOST_DATA
!$ACC END DATA
#endif

    ENDIF

#ifdef DEBUG
    WRITE (nerr,'(a,i4,a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' with broadcast number ', nbcast, ' successful.'

     IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL abort_mpi
    END IF
#endif
#endif

  END SUBROUTINE p_bcast_int_i8_1d

  SUBROUTINE p_bcast_int_2d (t_buffer, p_source, comm)

    INTEGER, INTENT(inout) :: t_buffer(:,:)
    INTEGER, INTENT(in)    :: p_source
    INTEGER, OPTIONAL, INTENT(in) :: comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

#ifdef DEBUG
    nbcast = nbcast+1
#endif

    IF (process_mpi_all_size == 1) THEN
       RETURN
    ELSE

#ifdef __USE_G2G
!$ACC DATA PRESENT( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
!$ACC HOST_DATA USE_DEVICE( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
#endif

       CALL MPI_BCAST (t_buffer, SIZE(t_buffer), p_int, p_source, &
            p_comm, p_error)

#ifdef __USE_G2G
!$ACC END HOST_DATA
!$ACC END DATA
#endif

    ENDIF

#ifdef DEBUG
    WRITE (nerr,'(a,i4,a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' with broadcast number ', nbcast, ' successful.'

     IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL abort_mpi
    END IF
#endif
#endif

  END SUBROUTINE p_bcast_int_2d

  SUBROUTINE p_bcast_int_3d (t_buffer, p_source, comm)

    INTEGER, INTENT(inout) :: t_buffer(:,:,:)
    INTEGER, INTENT(in)    :: p_source
    INTEGER, OPTIONAL, INTENT(in) :: comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

#ifdef DEBUG
    nbcast = nbcast+1
#endif

    IF (process_mpi_all_size == 1) THEN
       RETURN
    ELSE

#ifdef __USE_G2G
!$ACC DATA PRESENT( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
!$ACC HOST_DATA USE_DEVICE( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
#endif

       CALL MPI_BCAST (t_buffer, SIZE(t_buffer), p_int, p_source, &
            p_comm, p_error)

#ifdef __USE_G2G
!$ACC END HOST_DATA
!$ACC END DATA
#endif

    ENDIF

#ifdef DEBUG
    WRITE (nerr,'(a,i4,a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' with broadcast number ', nbcast, ' successful.'

     IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL abort_mpi
    END IF
#endif
#endif

  END SUBROUTINE p_bcast_int_3d

  SUBROUTINE p_bcast_int_4d (t_buffer, p_source, comm)

    INTEGER, INTENT(inout) :: t_buffer(:,:,:,:)
    INTEGER, INTENT(in)    :: p_source
    INTEGER, OPTIONAL, INTENT(in) :: comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

#ifdef DEBUG
    nbcast = nbcast+1
#endif

    IF (process_mpi_all_size == 1) THEN
       RETURN
    ELSE

#ifdef __USE_G2G
!$ACC DATA PRESENT( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
!$ACC HOST_DATA USE_DEVICE( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
#endif

       CALL MPI_BCAST (t_buffer, SIZE(t_buffer), p_int, p_source, &
            p_comm, p_error)

#ifdef __USE_G2G
!$ACC END HOST_DATA
!$ACC END DATA
#endif

    ENDIF

#ifdef DEBUG
    WRITE (nerr,'(a,i4,a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' with broadcast number ', nbcast, ' successful.'

     IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL abort_mpi
    END IF
#endif
#endif

  END SUBROUTINE p_bcast_int_4d

  SUBROUTINE p_bcast_int_7d (t_buffer, p_source, comm)

    INTEGER, INTENT(inout) :: t_buffer(:,:,:,:,:,:,:)
    INTEGER, INTENT(in)    :: p_source
    INTEGER, OPTIONAL, INTENT(in) :: comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

#ifdef DEBUG
    nbcast = nbcast+1
#endif

    IF (process_mpi_all_size == 1) THEN
       RETURN
    ELSE

#ifdef __USE_G2G
!$ACC DATA PRESENT( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
!$ACC HOST_DATA USE_DEVICE( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
#endif

       CALL MPI_BCAST (t_buffer, SIZE(t_buffer), p_int, p_source, &
            p_comm, p_error)

#ifdef __USE_G2G
!$ACC END HOST_DATA
!$ACC END DATA
#endif

    ENDIF

#ifdef DEBUG
    WRITE (nerr,'(a,i4,a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' with broadcast number ', nbcast, ' successful.'

     IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL abort_mpi
    END IF
#endif
#endif

  END SUBROUTINE p_bcast_int_7d

  SUBROUTINE p_bcast_bool (t_buffer, p_source, comm)

    LOGICAL, INTENT(inout) :: t_buffer
    INTEGER, INTENT(in)    :: p_source
    INTEGER, OPTIONAL, INTENT(in) :: comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

#ifdef DEBUG
    nbcast = nbcast+1
#endif

    IF (process_mpi_all_size == 1) THEN
       RETURN
    ELSE

#ifdef __USE_G2G
!$ACC DATA PRESENT( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
!$ACC HOST_DATA USE_DEVICE( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
#endif

       CALL MPI_BCAST (t_buffer, 1, p_bool, p_source, &
            p_comm, p_error)

#ifdef __USE_G2G
!$ACC END HOST_DATA
!$ACC END DATA
#endif

    ENDIF

#ifdef DEBUG
    WRITE (nerr,'(a,i4,a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' with broadcast number ', nbcast, ' successful.'

     IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL abort_mpi
    END IF
#endif
#endif

  END SUBROUTINE p_bcast_bool

  SUBROUTINE p_bcast_bool_1d (t_buffer, p_source, comm)

    LOGICAL, INTENT(inout) :: t_buffer(:)
    INTEGER, INTENT(in) :: p_source
    INTEGER, OPTIONAL, INTENT(in) :: comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

#ifdef DEBUG
    nbcast = nbcast+1
#endif

    IF (process_mpi_all_size == 1) THEN
       RETURN
    ELSE

#ifdef __USE_G2G
!$ACC DATA PRESENT( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
!$ACC HOST_DATA USE_DEVICE( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
#endif

       CALL MPI_BCAST (t_buffer, SIZE(t_buffer), p_bool, p_source, &
            p_comm, p_error)

#ifdef __USE_G2G
!$ACC END HOST_DATA
!$ACC END DATA
#endif

    ENDIF

#ifdef DEBUG
    WRITE (nerr,'(a,i4,a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' with broadcast number ', nbcast, ' successful.'

    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL abort_mpi
    END IF
#endif
#endif

  END SUBROUTINE p_bcast_bool_1d

  SUBROUTINE p_bcast_bool_2d (t_buffer, p_source, comm)

    LOGICAL, INTENT(inout) :: t_buffer(:,:)
    INTEGER, INTENT(in)    :: p_source
    INTEGER, OPTIONAL, INTENT(in) :: comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

#ifdef DEBUG
    nbcast = nbcast+1
#endif

    IF (process_mpi_all_size == 1) THEN
       RETURN
    ELSE

#ifdef __USE_G2G
!$ACC DATA PRESENT( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
!$ACC HOST_DATA USE_DEVICE( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
#endif

       CALL MPI_BCAST (t_buffer, SIZE(t_buffer), p_bool, p_source, &
            p_comm, p_error)

#ifdef __USE_G2G
!$ACC END HOST_DATA
!$ACC END DATA
#endif

    ENDIF

#ifdef DEBUG
    WRITE (nerr,'(a,i4,a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' with broadcast number ', nbcast, ' successful.'

     IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL abort_mpi
    END IF
#endif
#endif

  END SUBROUTINE p_bcast_bool_2d

  SUBROUTINE p_bcast_bool_3d (t_buffer, p_source, comm)

    LOGICAL, INTENT(inout) :: t_buffer(:,:,:)
    INTEGER, INTENT(in)    :: p_source
    INTEGER, OPTIONAL, INTENT(in) :: comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

#ifdef DEBUG
    nbcast = nbcast+1
#endif

    IF (process_mpi_all_size == 1) THEN
       RETURN
    ELSE

#ifdef __USE_G2G
!$ACC DATA PRESENT( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
!$ACC HOST_DATA USE_DEVICE( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
#endif

       CALL MPI_BCAST (t_buffer, SIZE(t_buffer), p_bool, p_source, &
            p_comm, p_error)

#ifdef __USE_G2G
!$ACC END HOST_DATA
!$ACC END DATA
#endif

    ENDIF

#ifdef DEBUG
    WRITE (nerr,'(a,i4,a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' with broadcast number ', nbcast, ' successful.'

    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL abort_mpi
    END IF
#endif
#endif

  END SUBROUTINE p_bcast_bool_3d

  SUBROUTINE p_bcast_bool_4d (t_buffer, p_source, comm)

    LOGICAL, INTENT(inout) :: t_buffer(:,:,:,:)
    INTEGER, INTENT(in)    :: p_source
    INTEGER, OPTIONAL, INTENT(in) :: comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

#ifdef DEBUG
    nbcast = nbcast+1
#endif

    IF (process_mpi_all_size == 1) THEN
       RETURN
    ELSE

#ifdef __USE_G2G
!$ACC DATA PRESENT( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
!$ACC HOST_DATA USE_DEVICE( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
#endif

       CALL MPI_BCAST (t_buffer, SIZE(t_buffer), p_bool, p_source, &
            p_comm, p_error)

#ifdef __USE_G2G
!$ACC END HOST_DATA
!$ACC END DATA
#endif

    ENDIF

#ifdef DEBUG
    WRITE (nerr,'(a,i4,a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' with broadcast number ', nbcast, ' successful.'

     IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL abort_mpi
    END IF
#endif
#endif

  END SUBROUTINE p_bcast_bool_4d

  SUBROUTINE p_bcast_char (t_buffer, p_source, comm)

    CHARACTER(len=*),  INTENT(inout) :: t_buffer
    INTEGER,           INTENT(in)    :: p_source
    INTEGER, OPTIONAL, INTENT(in) :: comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

#ifdef DEBUG
    nbcast = nbcast+1
#endif

    IF (process_mpi_all_size == 1) THEN
       RETURN
    ELSE

#ifdef __USE_G2G
!$ACC DATA PRESENT( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
!$ACC HOST_DATA USE_DEVICE( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
#endif

       CALL MPI_BCAST (t_buffer, LEN(t_buffer), p_char, p_source, &
            p_comm, p_error)

#ifdef __USE_G2G
!$ACC END HOST_DATA
!$ACC END DATA
#endif

    ENDIF

#ifdef DEBUG
    WRITE (nerr,'(a,i4,a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' with broadcast number ', nbcast, ' successful.'

     IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL abort_mpi
    END IF
#endif
#endif

  END SUBROUTINE p_bcast_char

  SUBROUTINE p_bcast_char_1d (t_buffer, p_source, comm)

    CHARACTER (*), INTENT(inout) :: t_buffer(:)
    INTEGER,       INTENT(in)    :: p_source
    INTEGER, OPTIONAL, INTENT(in) :: comm
    INTEGER                      :: lexlength, flength

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

#ifdef DEBUG
    nbcast = nbcast+1
#endif

    IF (process_mpi_all_size == 1) THEN
       RETURN
    ELSE
       lexlength=LEN(t_buffer(1))
       flength=SIZE(t_buffer)
       lexlength=lexlength*flength

#ifdef __USE_G2G
!$ACC DATA PRESENT( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
!$ACC HOST_DATA USE_DEVICE( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
#endif

       CALL MPI_BCAST (t_buffer, lexlength, p_char, p_source, p_comm, p_error)

#ifdef __USE_G2G
!$ACC END HOST_DATA
!$ACC END DATA
#endif

#ifdef DEBUG
       WRITE (nerr,'(a,i4,a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' with broadcast number ', nbcast, ' successful.'

       IF (p_error /= MPI_SUCCESS) THEN
          WRITE (nerr,'(a,i4,a)') ' MPI_BCAST from ', p_source, &
               ' failed.'
          WRITE (nerr,'(a,i4)') ' Error = ', p_error
          CALL abort_mpi
       END IF
#endif
    ENDIF
#endif

  END SUBROUTINE p_bcast_char_1d



  SUBROUTINE p_bcast_cchar (t_buffer, buflen, p_source, p_comm)

    USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_SIGNED_CHAR

    INTEGER,           INTENT(IN)    :: buflen
    INTEGER(C_SIGNED_CHAR), INTENT(INOUT) :: t_buffer(buflen)

    INTEGER,           INTENT(IN)    :: p_source
    INTEGER,           INTENT(IN)    :: p_comm

#ifndef NOMPI

#ifdef DEBUG
    nbcast = nbcast+1
#endif

    IF (process_mpi_all_size == 1) THEN
       RETURN
    ELSE

#ifdef __USE_G2G
!$ACC DATA PRESENT( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
!$ACC HOST_DATA USE_DEVICE( t_buffer ), IF ( i_am_accel_node .AND. acc_on )
#endif

       CALL MPI_BCAST (t_buffer, buflen, p_char, p_source, p_comm, p_error)

#ifdef __USE_G2G
!$ACC END HOST_DATA
!$ACC END DATA
#endif

#ifdef DEBUG
       WRITE (nerr,'(a,i4,a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' with broadcast number ', nbcast, ' successful.'

       IF (p_error /= MPI_SUCCESS) THEN
          WRITE (nerr,'(a,i4,a)') ' MPI_BCAST from ', p_source, &
               ' failed.'
          WRITE (nerr,'(a,i4)') ' Error = ', p_error
          CALL abort_mpi
       END IF
#endif
    ENDIF
#endif

  END SUBROUTINE p_bcast_cchar

  ! A bcast variant that handles ALLOCATABLE strings. Cannot be overloaded with p_bcast() because that would make the CALL ambigous.
  ! After this CALL, the string will always be ALLOCATED, even IF its length IS zero
  !
  !XXX: Deactivated due to a gfortran bug that looses the contents of allocatable strings on return.
!  SUBROUTINE p_bcast_achar(string, source, comm)
!    CHARACTER(:), ALLOCATABLE, INTENT(INOUT) :: string
!    INTEGER, VALUE :: source, comm
!
!    INTEGER :: length, error
!    CHARACTER(*), PARAMETER :: routine = modname//":p_bcast_achar"
!
!#ifndef NOMPI
!    ! inform the receivers about the length of the string
!    length = 0
!    IF(ALLOCATED(string)) length = LEN(string)
!    CALL p_bcast(length, source, comm)
!
!    ! ensure that the string IS ALLOCATED to the correct length
!    IF(ALLOCATED(string)) THEN
!        IF(length /= LEN(string) .OR. length == 0) DEALLOCATE(string)
!    END IF
!    IF(.NOT.ALLOCATED(string)) THEN
!        ALLOCATE(CHARACTER(LEN = length) :: string, STAT = error)
!        IF(error /= SUCCESS) CALL finish(routine, "memory allocation failed")
!    END IF
!
!    ! actually TRANSFER the string
!    CALL p_bcast(string, source, comm)
!#endif
!    END SUBROUTINE p_bcast_achar

  ! Collective CALL to determine whether this process IS a sender/receiver IN a broadcast operation.
  ! This routine IS robust IN the presence of inter-communicators (which IS its reason d'etre).
  SUBROUTINE p_get_bcast_role(root, communicator, lIsSender, lIsReceiver)
    INTEGER, VALUE :: root, communicator
    LOGICAL, INTENT(OUT) :: lIsSender, lIsReceiver

#ifndef NOMPI
    INTEGER :: error
    LOGICAL :: lIsInterCommunicator
    CHARACTER(LEN = *), PARAMETER :: routine = modname//":p_get_bcast_role"

    ! IS this an inter-communicator?
    CALL MPI_Comm_test_inter(communicator, lIsInterCommunicator, error)
    IF(error /= MPI_SUCCESS) CALL finish(routine, "MPI_Comm_test_inter() returned an error")
    IF(lIsInterCommunicator) THEN
        ! inter-communicator, we have to check for the special values MPI_ROOT AND MPI_PROC_NULL IN the root variable
        lIsSender = .FALSE.
        IF(root == MPI_ROOT) lIsSender = .TRUE.

        lIsReceiver = .FALSE.
        IF(root /= MPI_ROOT .AND. root /= MPI_PROC_NULL) lIsReceiver = .TRUE.
    ELSE
        ! intra-communicator, we have to compare the root variable against our own rank
        lIsSender = root == p_comm_rank(communicator)
        lIsReceiver = .NOT.lIsSender
    END IF
#else
    ! same behavior as IF there IS ONLY a single MPI process, which consequently must be the root of the bcast() operation
    lIsSender = .TRUE.
    lIsReceiver = .FALSE.
#endif
  END SUBROUTINE p_get_bcast_role


  ! probe implementation

  SUBROUTINE p_probe (p_tagcount, p_tagtable, p_source, &
       p_tag, p_count, comm)

    INTEGER,   INTENT(in)  :: p_tagcount, p_tagtable(:)
    INTEGER,   INTENT(out) :: p_source, p_tag, p_count
    INTEGER, OPTIONAL, INTENT(in) :: comm

#ifndef NOMPI
    INTEGER :: p_comm, i
    LOGICAL :: flag

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

    p_tag = -1
    DO WHILE (p_tag == -1)
       DO i = 1, p_tagcount
          CALL MPI_IPROBE (MPI_ANY_SOURCE, p_tagtable(i), p_comm, &
               flag, p_status, p_error)
#ifdef DEBUG
          IF (p_error /= MPI_SUCCESS) THEN
             WRITE (nerr,'(a,i4,a,i4,a)') ' MPI_IPROBE on ', my_process_mpi_all_id, &
                  ' for tag ', p_tagtable(i), ' failed.'
             WRITE (nerr,'(a,i4)') ' Error = ', p_error
             CALL abort_mpi
          END IF
#endif
          IF (flag) THEN
             p_source = p_status(MPI_SOURCE)
             p_tag = p_status(MPI_TAG)
             CALL MPI_GET_COUNT(p_status, p_real_dp, p_count, p_error)
#ifdef DEBUG
             IF (p_error /= MPI_SUCCESS) THEN
                WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_GET_COUNT on ', &
                  & my_process_mpi_all_id, &
                     ' for tag ', p_tag, ' from ' , p_source, ' failed.'
                WRITE (nerr,'(a,i4)') ' Error = ', p_error
                CALL abort_mpi
             END IF
#endif
             EXIT
          ELSE
             p_tag = -1
          END IF
       END DO
    END DO

#endif

  END SUBROUTINE p_probe
  !------------------------------------------------------

  !------------------------------------------------------
  SUBROUTINE p_wait(request)
    INTEGER, INTENT(INOUT), OPTIONAL :: request
#ifndef NOMPI
    INTEGER :: p_status_wait(MPI_STATUS_SIZE,p_irequest)

    IF (PRESENT(request)) THEN
      CALL MPI_WAIT(request,p_status_wait(:,1), p_error)
      CALL p_clear_request(request)
    ELSE
      CALL MPI_WAITALL(p_irequest, p_request, p_status_wait, p_error)
      p_irequest = 0
    END IF
#endif
  END SUBROUTINE p_wait

  SUBROUTINE p_wait_any(return_pe)

    INTEGER, INTENT(out) :: return_pe
#ifndef NOMPI
    INTEGER :: i

    CALL MPI_WAITANY(p_irequest, p_request, i, p_status, p_error)
    IF (i == MPI_UNDEFINED) THEN
      p_irequest = 0
      return_pe = -1
    ELSE
      return_pe = p_status(MPI_SOURCE)
    ENDIF
#else
    return_pe = 0
#endif
  END SUBROUTINE p_wait_any
  !------------------------------------------------------


  !------------------------------------------------------
  FUNCTION p_test() RESULT(ret)
    LOGICAL :: ret
#ifndef NOMPI
    INTEGER :: p_status_wait(MPI_STATUS_SIZE,p_irequest)
    CALL MPI_TESTALL(p_irequest, p_request, ret, p_status_wait, p_error)
#else
    ret = .TRUE.
#endif
  END FUNCTION p_test


  !------------------------------------------------------
  SUBROUTINE p_barrier (comm)
  INTEGER ,INTENT(IN) ,OPTIONAL :: comm
#ifndef NOMPI
    INTEGER :: com
    com = MPI_COMM_WORLD; IF(PRESENT(comm)) com = comm
    CALL MPI_BARRIER (com, p_error)

!#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a)') ' MPI_BARRIER on ', my_process_mpi_all_id, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL abort_mpi
    END IF
!#endif
#endif

  END SUBROUTINE p_barrier
  !------------------------------------------------------

  !------------------------------------------------------
  SUBROUTINE global_mpi_barrier
#ifndef NOMPI
    CALL MPI_BARRIER (global_mpi_communicator, p_error)

!#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a)') ' global_mpi_barrier on ', get_my_global_mpi_id(), ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL abort_mpi
    END IF
!#endif
#endif

  END SUBROUTINE global_mpi_barrier
  !------------------------------------------------------

  !------------------------------------------------------
  SUBROUTINE work_mpi_barrier
#ifndef NOMPI
    CALL MPI_BARRIER (p_comm_work, p_error)

!#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a)') ' global_mpi_barrier on ', get_my_global_mpi_id(), ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL abort_mpi
    END IF
!#endif
#endif

  END SUBROUTINE work_mpi_barrier
  !------------------------------------------------------

  !------------------------------------------------------
  !> computes a global sum of real single precision numbers
  !
  !  perform an ALLREDUCE operation
  FUNCTION p_sum_dp_0s (zfield, comm) RESULT (p_sum)

    REAL(sp)                        :: p_sum
    REAL(sp),  INTENT(in)           :: zfield
    INTEGER,   INTENT(in)           :: comm
#ifndef NOMPI
    INTEGER :: p_comm

    p_comm = comm

    IF (my_process_is_mpi_parallel()) THEN
        CALL MPI_ALLREDUCE (zfield, p_sum, 1, p_real_sp, &
          MPI_SUM, p_comm, p_error)
    ELSE
       p_sum = zfield
    END IF
#else
    p_sum = zfield
#endif
  END FUNCTION p_sum_dp_0s
  !------------------------------------------------------

  !------------------------------------------------------
  !> computes a global sum of real numbers
  !
  ! @param[in]    root     (Optional:) root PE, otherwise we perform an
  !                                    ALLREDUCE operation
  FUNCTION p_sum_dp_0d (zfield, comm, root) RESULT (p_sum)

    REAL(dp)                        :: p_sum
    REAL(dp),  INTENT(in)           :: zfield
    INTEGER,   INTENT(in)           :: comm
    INTEGER,   INTENT(in), OPTIONAL :: root
#ifndef NOMPI
    INTEGER :: p_comm, my_rank

    p_comm = comm

    IF (my_process_is_mpi_parallel()) THEN
      IF (PRESENT(root)) THEN
        CALL MPI_REDUCE (zfield, p_sum, 1, p_real_dp, &
          MPI_SUM, root, p_comm, p_error)
        ! get local PE identification
        CALL MPI_COMM_RANK (p_comm, my_rank, p_error)
        ! do not use the result on all the other ranks:
        IF (root /= my_rank)  p_sum = zfield
      ELSE
        CALL MPI_ALLREDUCE (zfield, p_sum, 1, p_real_dp, &
          MPI_SUM, p_comm, p_error)
      END IF
    ELSE
       p_sum = zfield
    END IF
#else
    p_sum = zfield
#endif
  END FUNCTION p_sum_dp_0d
  !------------------------------------------------------

  !------------------------------------------------------
  FUNCTION p_sum_dp_1d (zfield, comm, root) RESULT (p_sum)

    REAL(dp),          INTENT(in) :: zfield(:)
    INTEGER, OPTIONAL, INTENT(in) :: comm, root
    REAL(dp)                      :: p_sum (SIZE(zfield))

#ifndef NOMPI
    INTEGER :: p_comm, my_rank

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

    IF (my_process_is_mpi_all_parallel()) THEN
      IF (PRESENT(root)) THEN
        CALL mpi_reduce(zfield, p_sum, SIZE(zfield), p_real_dp, &
             mpi_sum, root, p_comm, p_error)
        ! get local PE identification
        CALL mpi_comm_rank(p_comm, my_rank, p_error)
        ! do not use the result on all the other ranks:
        IF (root /= my_rank) p_sum = zfield
      ELSE
        CALL mpi_allreduce (zfield, p_sum, SIZE(zfield), p_real_dp, &
             mpi_sum, p_comm, p_error)
      END IF
    ELSE
       p_sum = zfield
    END IF
#else
    p_sum = zfield
#endif

  END FUNCTION p_sum_dp_1d

  FUNCTION p_sum_i8_1d (kfield, comm) RESULT (p_sum)

    INTEGER(i8),       INTENT(in) :: kfield(:)
    INTEGER, OPTIONAL, INTENT(in) :: comm
    INTEGER(i8)                   :: p_sum (SIZE(kfield))

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

    IF (my_process_is_mpi_all_parallel()) THEN
       CALL MPI_ALLREDUCE (kfield, p_sum, SIZE(kfield), p_int_i8, &
            MPI_SUM, p_comm, p_error)
    ELSE
       p_sum = kfield
    END IF
#else
    p_sum = kfield
#endif

  END FUNCTION p_sum_i8_1d

  FUNCTION p_sum_i_1d (kfield, comm) RESULT (p_sum)

    INTEGER,       INTENT(in) :: kfield(:)
    INTEGER, OPTIONAL, INTENT(in) :: comm
    INTEGER                   :: p_sum (SIZE(kfield))

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

    IF (my_process_is_mpi_all_parallel()) THEN
       CALL MPI_ALLREDUCE (kfield, p_sum, SIZE(kfield), p_int, &
            MPI_SUM, p_comm, p_error)
    ELSE
       p_sum = kfield
    END IF
#else
    p_sum = kfield
#endif

  END FUNCTION p_sum_i_1d


  FUNCTION p_sum_i_0d (kfield, comm) RESULT (p_sum)

    INTEGER,       INTENT(in) :: kfield
    INTEGER, OPTIONAL, INTENT(in) :: comm
    INTEGER                   :: p_sum

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

    IF (my_process_is_mpi_all_parallel()) THEN
       CALL MPI_ALLREDUCE (kfield, p_sum, 1, p_int, &
            MPI_SUM, p_comm, p_error)
    ELSE
       p_sum = kfield
    END IF
#else
    p_sum = kfield
#endif

  END FUNCTION p_sum_i_0d

  FUNCTION p_global_sum_1d (zfield, comm) RESULT (p_sum)

    REAL(dp),          INTENT(in) :: zfield(:)
    INTEGER, OPTIONAL, INTENT(in) :: comm
    REAL(dp)                      :: p_sum
    REAL(dp)                      :: pe_sums(SIZE(zfield))

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

    IF (my_process_is_mpi_all_parallel()) THEN
       CALL MPI_REDUCE (zfield, pe_sums, SIZE(zfield), p_real_dp, &
            MPI_SUM, process_mpi_root_id, p_comm, p_error)
       p_sum = SUM(pe_sums)
    ELSE
       p_sum = SUM(zfield)
    END IF
#else
    p_sum = SUM(zfield)
#endif

  END FUNCTION p_global_sum_1d

  FUNCTION p_field_sum_1d (zfield, comm) RESULT (p_sum)

    REAL(dp),          INTENT(in) :: zfield(:)
    INTEGER, OPTIONAL, INTENT(in) :: comm
    REAL(dp)                      :: p_sum (SIZE(zfield))

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

    IF (my_process_is_mpi_all_parallel()) THEN
       CALL MPI_REDUCE (zfield, p_sum, SIZE(zfield), p_real_dp, &
            MPI_SUM, process_mpi_root_id, p_comm, p_error)
       IF (.NOT. my_process_is_stdio()) p_sum = 0.0_dp
    ELSE
       p_sum = zfield
    END IF
#else
    p_sum = zfield
#endif

  END FUNCTION p_field_sum_1d

  !> common code of min/max reductions
  SUBROUTINE p_minmax_common(in_field, out_field, n, op, loc_op, &
       proc_id, keyval, comm, root)
    INTEGER, INTENT(in) :: n, op, loc_op
    REAL(dp), INTENT(in) :: in_field(n)
    REAL(dp), INTENT(out) :: out_field(n)

    INTEGER, OPTIONAL, INTENT(inout) :: proc_id(n)
    INTEGER, OPTIONAL, INTENT(inout) :: keyval(n)
    INTEGER, OPTIONAL, INTENT(in)    :: root
    INTEGER, OPTIONAL, INTENT(in)    :: comm

#ifndef NOMPI
    INTEGER  :: p_comm, rank, comm_size
    LOGICAL :: compute_ikey
    INTEGER, ALLOCATABLE  :: meta_info(:), ikey(:)
#ifndef SLOW_MPI_MAXMINLOC
    DOUBLE PRECISION, ALLOCATABLE :: val_loc(:,:,:)
#endif

    IF (PRESENT(comm)) THEN
       p_comm = comm
       comm_size = p_comm_size(comm)
    ELSE
       p_comm = process_mpi_all_comm
       comm_size = process_mpi_all_size
    ENDIF

    IF (comm_size > 1) THEN

      IF (PRESENT(proc_id) .OR. PRESENT(keyval)) THEN
        ! encode meta information
        ALLOCATE(meta_info(n), ikey(n))
        meta_info = 0
        IF (PRESENT(keyval))  meta_info = meta_info + keyval*comm_size
        IF (PRESENT(proc_id)) meta_info = meta_info + proc_id
        ! use mpi_minloc to transfer additional data
#ifndef SLOW_MPI_MAXMINLOC
        ALLOCATE(val_loc(2, n, 2))
#endif
        IF (PRESENT(root)) THEN
          CALL MPI_COMM_RANK(p_comm, rank, p_error)
#ifdef SLOW_MPI_MAXMINLOC
          ! on BG/Q, {max|min}loc is slow
          CALL mpi_allreduce(in_field, out_field, n, mpi_double_precision, &
               op, p_comm, p_error)
          ikey = MERGE(meta_info, HUGE(1), in_field == out_field)
          CALL mpi_reduce(ikey, meta_info, n, mpi_integer, &
               mpi_min, root, p_comm, p_error)
#else
          val_loc(1, :, 1) = DBLE(in_field)
          val_loc(2, :, 1) = DBLE(meta_info)
          CALL mpi_reduce(val_loc(:, :, 1), val_loc(:, :, 2), &
               n, mpi_2double_precision, loc_op, root, p_comm, p_error)
          IF (rank == root) THEN
             out_field = val_loc(1, :, 2)
          ELSE
             out_field = 0.
          END IF
#endif
          compute_ikey = rank == root
        ELSE
#ifdef SLOW_MPI_MAXMINLOC
          CALL mpi_allreduce(in_field, out_field, n, mpi_double_precision, &
               op, p_comm, p_error)
          ikey = MERGE(meta_info, HUGE(1), in_field == out_field)
          CALL mpi_allreduce(ikey, meta_info, n, mpi_integer, &
               mpi_min, p_comm, p_error)
#else
          val_loc(1, :, 1) = DBLE(in_field)
          val_loc(2, :, 1) = DBLE(meta_info)
          CALL mpi_allreduce(val_loc(:, :, 1), val_loc(:, :, 2), &
               n, mpi_2double_precision, loc_op, p_comm, p_error)
          out_field = val_loc(1, :, 2)
#endif
          compute_ikey = .TRUE.
        END IF
        ! decode meta info:
        IF (compute_ikey) THEN
#ifdef SLOW_MPI_MAXMINLOC
          ikey = meta_info / comm_size
          IF (PRESENT(proc_id)) proc_id = mod(meta_info,comm_size)
#else
          ikey = NINT(val_loc(2, :, 2)) / comm_size
          IF (PRESENT(proc_id)) proc_id = mod(nint(val_loc(2, :, 2)),comm_size)
#endif
          IF (PRESENT(keyval)) keyval = ikey
        END IF
      ELSE
        ! compute simple (standard) minimum
        IF (PRESENT(root)) THEN
          CALL mpi_reduce(in_field, out_field, n, p_real_dp, &
               op, root, p_comm, p_error)
        ELSE
          CALL mpi_allreduce(in_field, out_field, n, p_real_dp, &
               op, p_comm, p_error)
        END IF
     END IF
    ELSE
      out_field = in_field
    END IF
#else
    out_field = in_field
#endif
  END SUBROUTINE p_minmax_common

  !> computes a global maximum of real numbers
  !
  ! @param[out]   proc_id  (Optional:) PE number of maximum value
  ! @param[inout] keyval   (Optional:) additional meta information
  ! @param[in]    root     (Optional:) root PE, otherwise we perform an
  !                                    ALL-TO-ALL operation
  !
  ! The parameter @p keyval can be used to communicate
  ! additional data on the maximum value, e.g., the level
  ! index where the maximum occurred.
  !
  FUNCTION p_max_0d (zfield, proc_id, keyval, comm, root) RESULT (p_max)

    REAL(dp)                         :: p_max
    REAL(dp),          INTENT(in)    :: zfield
    INTEGER, OPTIONAL, INTENT(inout) :: proc_id
    INTEGER, OPTIONAL, INTENT(inout) :: keyval
    INTEGER, OPTIONAL, INTENT(in)    :: root
    INTEGER, OPTIONAL, INTENT(in)    :: comm

    REAL(dp) :: temp_in(1), temp_out(1)
    INTEGER :: temp_keyval(1), temp_proc_id(1)
#ifdef NOMPI
    INTEGER, PARAMETER :: mpi_max = 3333, mpi_maxloc = 3333
#endif
    temp_in(1) = zfield
    IF (PRESENT(proc_id) .AND. PRESENT(keyval)) THEN
      temp_keyval(1) = keyval; temp_proc_id(1) = proc_id
      CALL p_minmax_common(temp_in, temp_out, 1, mpi_max, mpi_maxloc, &
           proc_id=temp_proc_id, keyval=temp_keyval, comm=comm, root=root)
      keyval = temp_keyval(1); proc_id = temp_proc_id(1)
    ELSE IF (PRESENT(proc_id)) THEN
      temp_proc_id(1) = proc_id
      CALL p_minmax_common(temp_in, temp_out, 1, mpi_max, mpi_maxloc, &
           proc_id=temp_proc_id, comm=comm, root=root)
      proc_id = temp_proc_id(1)
    ELSE IF (PRESENT(keyval)) THEN
      temp_keyval(1) = keyval
      CALL p_minmax_common(temp_in, temp_out, 1, mpi_max, mpi_maxloc, &
         keyval=temp_keyval, comm=comm, root=root)
      keyval = temp_keyval(1)
    ELSE ! .not. present(keyval) .and. .not. present(proc_id)
      CALL p_minmax_common(temp_in, temp_out, 1, mpi_max, mpi_maxloc, &
           comm=comm, root=root)
    END IF

    p_max = temp_out(1)

  END FUNCTION p_max_0d

  FUNCTION p_max_int_0d (zfield, comm) RESULT (p_max)

    INTEGER                          :: p_max
    INTEGER,           INTENT(in)    :: zfield
    INTEGER, OPTIONAL, INTENT(in)    :: comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

    IF (my_process_is_mpi_all_parallel()) THEN
       CALL MPI_ALLREDUCE (zfield, p_max, 1, p_int, &
            MPI_MAX, p_comm, p_error)
    ELSE
       p_max = zfield
    END IF
#else
    p_max = zfield
#endif

  END FUNCTION p_max_int_0d

  !> computes a global maximum of a real number array
  !
  ! @param[out]   proc_id  (Optional:) PE number of maximum value
  ! @param[inout] keyval   (Optional:) additional meta information
  ! @param[in]    root     (Optional:) root PE, otherwise we perform an
  !                                    ALL-TO-ALL operation
  !
  ! The parameter @p keyval can be used to communicate
  ! additional data on the maximum value, e.g., the level
  ! index where the maximum occurred.
  !
  FUNCTION p_max_1d (zfield, proc_id, keyval, comm, root) RESULT (p_max)

    REAL(dp),          INTENT(in)    :: zfield(:)
    INTEGER, OPTIONAL, INTENT(inout) :: proc_id(SIZE(zfield))
    INTEGER, OPTIONAL, INTENT(inout) :: keyval(SIZE(zfield))
    INTEGER, OPTIONAL, INTENT(in)    :: root
    INTEGER, OPTIONAL, INTENT(in)    :: comm
    REAL(dp)                         :: p_max (SIZE(zfield))
#ifdef NOMPI
    INTEGER, PARAMETER :: mpi_max = 3333, mpi_maxloc = 3333
#endif

    CALL p_minmax_common(zfield, p_max, SIZE(zfield), mpi_max, mpi_maxloc, &
           proc_id=proc_id, keyval=keyval, comm=comm, root=root)

  END FUNCTION p_max_1d

  ! Computes maximum of a 1D field of integers.
  !
  ! @param[in] root Optional root PE, otherwise we perform an
  !                 ALL-TO-ALL operation
  FUNCTION p_max_int_1d (zfield, comm, root) RESULT (p_max)

    INTEGER,           INTENT(in) :: zfield(:)
    INTEGER, OPTIONAL, INTENT(in) :: comm
    INTEGER, OPTIONAL, INTENT(in) :: root
    INTEGER                       :: p_max (SIZE(zfield))

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

    IF (my_process_is_mpi_all_parallel()) THEN
      IF (PRESENT(root)) THEN
        CALL MPI_REDUCE (zfield, p_max, SIZE(zfield), p_int, &
          MPI_MAX, root, p_comm, p_error)
      ELSE
        CALL MPI_ALLREDUCE (zfield, p_max, SIZE(zfield), p_int, &
          MPI_MAX, p_comm, p_error)
      END IF ! if present(root)
    ELSE
       p_max = zfield
    END IF
#else
    p_max = zfield
#endif

  END FUNCTION p_max_int_1d

  FUNCTION p_max_2d (zfield, comm) RESULT (p_max)

    REAL(dp),          INTENT(in) :: zfield(:,:)
    INTEGER, OPTIONAL, INTENT(in) :: comm
    REAL(dp)                      :: p_max (SIZE(zfield,1),SIZE(zfield,2))

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

    IF (my_process_is_mpi_all_parallel()) THEN
       CALL MPI_ALLREDUCE (zfield, p_max, SIZE(zfield), p_real_dp, &
            MPI_MAX, p_comm, p_error)
    ELSE
       p_max = zfield
    END IF
#else
    p_max = zfield
#endif

  END FUNCTION p_max_2d

  FUNCTION p_max_3d (zfield, comm) RESULT (p_max)

    REAL(dp),          INTENT(in) :: zfield(:,:,:)
    INTEGER, OPTIONAL, INTENT(in) :: comm
    REAL(dp)                      :: p_max (SIZE(zfield,1),SIZE(zfield,2)&
                                           ,SIZE(zfield,3))

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

    IF (my_process_is_mpi_all_parallel()) THEN
       CALL MPI_ALLREDUCE (zfield, p_max, SIZE(zfield), p_real_dp, &
            MPI_MAX, p_comm, p_error)
    ELSE
       p_max = zfield
    END IF
#else
    p_max = zfield
#endif

  END FUNCTION p_max_3d


  !> computes a global minimum of real numbers
  !
  ! @param[out]   proc_id  (Optional:) PE number of maximum value
  ! @param[inout] keyval   (Optional:) additional meta information
  ! @param[in]    root     (Optional:) root PE, otherwise we perform an
  !                                    ALL-TO-ALL operation
  !
  ! The parameter @p keyval can be used to communicate
  ! additional data on the maximum value, e.g., the level
  ! index where the maximum occurred.
  !
  FUNCTION p_min_0d (zfield, proc_id, keyval, comm, root) RESULT (p_min)

    REAL(dp)                         :: p_min
    REAL(dp),          INTENT(in)    :: zfield
    INTEGER, OPTIONAL, INTENT(inout) :: proc_id
    INTEGER, OPTIONAL, INTENT(inout) :: keyval
    INTEGER, OPTIONAL, INTENT(in)    :: root
    INTEGER, OPTIONAL, INTENT(in)    :: comm

    REAL(dp) :: temp_in(1), temp_out(1)
    INTEGER :: temp_keyval(1), temp_proc_id(1)
#ifdef NOMPI
    INTEGER, PARAMETER :: mpi_min = 4444, mpi_minloc = 4444
#endif
    temp_in(1) = zfield
    IF (PRESENT(proc_id) .AND. PRESENT(keyval)) THEN
      temp_keyval(1) = keyval; temp_proc_id(1) = proc_id
      CALL p_minmax_common(temp_in, temp_out, 1, mpi_min, mpi_minloc, &
           proc_id=temp_proc_id, keyval=temp_keyval, comm=comm, root=root)
      keyval = temp_keyval(1); proc_id = temp_proc_id(1)
    ELSE IF (PRESENT(proc_id)) THEN
      temp_proc_id(1) = proc_id
      CALL p_minmax_common(temp_in, temp_out, 1, mpi_min, mpi_minloc, &
           proc_id=temp_proc_id, comm=comm, root=root)
      proc_id = temp_proc_id(1)
    ELSE IF (PRESENT(keyval)) THEN
      temp_keyval(1) = keyval
      CALL p_minmax_common(temp_in, temp_out, 1, mpi_min, mpi_minloc, &
         keyval=temp_keyval, comm=comm, root=root)
      keyval = temp_keyval(1)
    ELSE ! .not. present(keyval) .and. .not. present(proc_id)
      CALL p_minmax_common(temp_in, temp_out, 1, mpi_min, mpi_minloc, &
           comm=comm, root=root)
    END IF

    p_min = temp_out(1)

  END FUNCTION p_min_0d

  FUNCTION p_min_int_0d (zfield, comm) RESULT (p_min)

    INTEGER,           INTENT(in) :: zfield
    INTEGER, OPTIONAL, INTENT(in) :: comm
    INTEGER                       :: p_min
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

    IF (my_process_is_mpi_all_parallel()) THEN
       CALL MPI_ALLREDUCE (zfield, p_min, 1, p_int, &
            MPI_MIN, p_comm, p_error)
    ELSE
       p_min = zfield
    END IF
#else
    p_min = zfield
#endif

  END FUNCTION p_min_int_0d

  FUNCTION p_min_1d (zfield, proc_id, keyval, comm, root) RESULT (p_min)

    REAL(dp),          INTENT(in) :: zfield(:)
    INTEGER, OPTIONAL, INTENT(inout) :: proc_id(SIZE(zfield))
    INTEGER, OPTIONAL, INTENT(inout) :: keyval(SIZE(zfield))
    INTEGER, OPTIONAL, INTENT(in) :: root
    INTEGER, OPTIONAL, INTENT(in) :: comm
    REAL(dp)                      :: p_min (SIZE(zfield))
#ifdef NOMPI
    INTEGER, PARAMETER :: mpi_min = 4444, mpi_minloc = 4444
#endif

    CALL p_minmax_common(zfield, p_min, SIZE(zfield), mpi_min, mpi_minloc, &
           proc_id=proc_id, keyval=keyval, comm=comm, root=root)

  END FUNCTION p_min_1d

  FUNCTION p_min_int_1d (zfield, comm) RESULT (p_min)

    INTEGER,           INTENT(in) :: zfield(:)
    INTEGER, OPTIONAL, INTENT(in) :: comm
    INTEGER                       :: p_min (SIZE(zfield))

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

    IF (my_process_is_mpi_all_parallel()) THEN
       CALL MPI_ALLREDUCE (zfield, p_min, SIZE(zfield), p_int, &
            MPI_MIN, p_comm, p_error)
    ELSE
       p_min = zfield
    END IF
#else
    p_min = zfield
#endif

  END FUNCTION p_min_int_1d

  FUNCTION p_min_2d (zfield, comm) RESULT (p_min)

    REAL(dp),          INTENT(in) :: zfield(:,:)
    INTEGER, OPTIONAL, INTENT(in) :: comm
    REAL(dp)                      :: p_min (SIZE(zfield,1),SIZE(zfield,2))

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

    IF (my_process_is_mpi_all_parallel()) THEN
       CALL MPI_ALLREDUCE (zfield, p_min, SIZE(zfield), p_real_dp, &
            MPI_MIN, p_comm, p_error)
    ELSE
       p_min = zfield
    END IF
#else
    p_min = zfield
#endif

  END FUNCTION p_min_2d

  FUNCTION p_min_3d (zfield, comm) RESULT (p_min)

    REAL(dp),          INTENT(in) :: zfield(:,:,:)
    INTEGER, OPTIONAL, INTENT(in) :: comm
    REAL(dp)                      :: p_min (SIZE(zfield,1),SIZE(zfield,2)&
                                           ,SIZE(zfield,3))
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

    IF (my_process_is_mpi_all_parallel()) THEN
       CALL MPI_ALLREDUCE (zfield, p_min, SIZE(zfield), p_real_dp, &
            MPI_MIN, p_comm, p_error)
    ELSE
       p_min = zfield
    END IF
#else
    p_min = zfield
#endif

  END FUNCTION p_min_3d


  ! Computes maximum of a 1D field of integers.
  !
  ! @param[in] root Optional root PE, otherwise we perform an
  !                 ALL-TO-ALL operation
  SUBROUTINE p_allreduce_max_int_1d (zfield, comm, root)
    INTEGER, INTENT(inout) :: zfield(:)
    INTEGER, INTENT(in) :: comm
    INTEGER, INTENT(in), OPTIONAL :: root

#if !defined(NOMPI)
    INTEGER :: p_comm, p_error

    p_comm = comm
    IF (PRESENT(root)) THEN
      CALL MPI_REDUCE (MPI_IN_PLACE, zfield, SIZE(zfield), p_int, &
        MPI_MAX, root, p_comm, p_error)
    ELSE
      CALL MPI_ALLREDUCE (MPI_IN_PLACE, zfield, SIZE(zfield), p_int, &
        MPI_MAX, p_comm, p_error)
    END IF ! if present(root)
#endif
  END SUBROUTINE p_allreduce_max_int_1d


  !---------------------------------------------------------------------------------------------------------------------------------
  !> wrapper for MPI_Scatter
  !---------------------------------------------------------------------------------------------------------------------------------
  SUBROUTINE p_scatter_real_1d1d(sendbuf, recvbuf, p_src, comm)
    REAL(dp),          INTENT(inout) :: sendbuf(:), recvbuf(:)
    INTEGER,           INTENT(in) :: p_src
    INTEGER, OPTIONAL, INTENT(in) :: comm

#ifndef NOMPI
    CHARACTER(*), PARAMETER :: routine = TRIM("mo_mpi:p_scatter_real_1d1d")
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

    CALL MPI_Scatter(sendbuf, SIZE(recvbuf), p_real_dp, &
    &                recvbuf, SIZE(recvbuf), p_real_dp, &
    &                p_src, p_comm, p_error)
    IF(p_error /= MPI_SUCCESS) CALL finish(routine, 'Error in MPI_Scatter operation!')
#else
     recvbuf = sendbuf
#endif
   END SUBROUTINE p_scatter_real_1d1d

  !---------------------------------------------------------------------------------------------------------------------------------
  !> wrapper for MPI_Scatter
  !---------------------------------------------------------------------------------------------------------------------------------
  SUBROUTINE p_scatter_single_1d1d(sendbuf, recvbuf, p_src, comm)
    REAL(sp),          INTENT(inout) :: sendbuf(:), recvbuf(:)
    INTEGER,           INTENT(in) :: p_src
    INTEGER, OPTIONAL, INTENT(in) :: comm

#ifndef NOMPI
    CHARACTER(*), PARAMETER :: routine = TRIM("mo_mpi:p_scatter_single_1d1d")
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

    CALL MPI_Scatter(sendbuf, SIZE(recvbuf), p_real_sp, &
    &                recvbuf, SIZE(recvbuf), p_real_sp, &
    &                p_src, p_comm, p_error)
    IF(p_error /= MPI_SUCCESS) CALL finish(routine, 'Error in MPI_Scatter operation!')
#else
     recvbuf = sendbuf
#endif
   END SUBROUTINE p_scatter_single_1d1d

  !---------------------------------------------------------------------------------------------------------------------------------
  !> wrapper for MPI_Scatter
  !---------------------------------------------------------------------------------------------------------------------------------
  SUBROUTINE p_scatter_int_1d1d(sendbuf, recvbuf, p_src, comm)
    INTEGER,           INTENT(inout) :: sendbuf(:), recvbuf(:)
    INTEGER,           INTENT(in) :: p_src
    INTEGER, OPTIONAL, INTENT(in) :: comm

#ifndef NOMPI
    CHARACTER(*), PARAMETER :: routine = TRIM("mo_mpi:p_scatter_int_1d1d")
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

    CALL MPI_Scatter(sendbuf, SIZE(recvbuf), p_int, &
    &                recvbuf, SIZE(recvbuf), p_int, &
    &                p_src, p_comm, p_error)
    IF(p_error /= MPI_SUCCESS) CALL finish(routine, 'Error in MPI_Scatter operation!')
#else
     recvbuf = sendbuf
#endif
   END SUBROUTINE p_scatter_int_1d1d

  SUBROUTINE p_gather_real_0d1d (sendbuf, recvbuf, p_dest, comm)
    REAL(dp),          INTENT(inout) :: sendbuf, recvbuf(:)
    INTEGER,           INTENT(in) :: p_dest
    INTEGER, OPTIONAL, INTENT(in) :: comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

     CALL MPI_GATHER(sendbuf, 1, p_real_dp, &
                     recvbuf, 1, p_real_dp, &
                     p_dest, p_comm, p_error)
#else
     recvbuf(:) = sendbuf
#endif
   END SUBROUTINE p_gather_real_0d1d

  SUBROUTINE p_gather_real_1d2d (sendbuf, recvbuf, p_dest, comm)
    REAL(dp),          INTENT(inout) :: sendbuf(:), recvbuf(:,:)
    INTEGER,           INTENT(in) :: p_dest
    INTEGER, OPTIONAL, INTENT(in) :: comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

     CALL MPI_GATHER(sendbuf, SIZE(sendbuf), p_real_dp, &
                     recvbuf, SIZE(sendbuf), p_real_dp, &
                     p_dest, p_comm, p_error)
#else
     recvbuf(:,LBOUND(recvbuf,2)) = sendbuf(:)
#endif
   END SUBROUTINE p_gather_real_1d2d

   SUBROUTINE p_gather_real_5d6d (sendbuf, recvbuf, p_dest, comm)

     REAL(dp)                      :: sendbuf(:,:,:,:,:), recvbuf(:,:,:,:,:,:)
     INTEGER,           INTENT(in) :: p_dest
     INTEGER, OPTIONAL, INTENT(in) :: comm

#ifndef NOMPI
     INTEGER :: p_comm

     IF (PRESENT(comm)) THEN
       p_comm = comm
     ELSE
       p_comm = process_mpi_all_comm
     ENDIF

     CALL MPI_GATHER(sendbuf, SIZE(sendbuf), p_real_dp, &
                     recvbuf, SIZE(sendbuf), p_real_dp, &
                     p_dest, p_comm, p_error)

     IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a)') ' p_gather_real_5d6d failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       STOP
     END IF

#else
     recvbuf(:,:,:,:,:,LBOUND(recvbuf,6)) = sendbuf(:,:,:,:,:)
#endif
   END SUBROUTINE p_gather_real_5d6d


  SUBROUTINE p_gather_real_1d1d (sendbuf, recvbuf, p_dest, comm)
    REAL(dp),          INTENT(inout) :: sendbuf(:), recvbuf(:)
    INTEGER,           INTENT(in) :: p_dest
    INTEGER, OPTIONAL, INTENT(in) :: comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

     CALL MPI_GATHER(sendbuf, SIZE(sendbuf), p_real_dp, &
                     recvbuf, SIZE(sendbuf), p_real_dp, &
                     p_dest, p_comm, p_error)
#else
     recvbuf(:) = sendbuf(:)
#endif
   END SUBROUTINE p_gather_real_1d1d


  !---------------------------------------------------------------------------------------------------------------------------------
  !> wrapper for MPI_Gather()
  !---------------------------------------------------------------------------------------------------------------------------------
   SUBROUTINE p_gather_int_0d1d (sendbuf, recvbuf, p_dest, comm)
     INTEGER,           INTENT(inout) :: sendbuf, recvbuf(:)
     INTEGER,           INTENT(in) :: p_dest
     INTEGER, OPTIONAL, INTENT(in) :: comm

#ifndef NOMPI
     CHARACTER(*), PARAMETER :: routine = TRIM("mo_mpi:p_gather_int_0d1d")
     INTEGER :: p_comm

     IF (PRESENT(comm)) THEN
       p_comm = comm
     ELSE
       p_comm = process_mpi_all_comm
     ENDIF

     CALL MPI_GATHER(sendbuf, 1, MPI_INTEGER, &
       &             recvbuf, 1, MPI_INTEGER, &
       &             p_dest, p_comm, p_error)
     IF (p_error /=  MPI_SUCCESS) CALL finish (routine, 'Error in MPI_GATHER operation!')
#else
     recvbuf = sendbuf
#endif
   END SUBROUTINE p_gather_int_0d1d


  !---------------------------------------------------------------------------------------------------------------------------------
  !> wrapper for MPI_Gather()
  !---------------------------------------------------------------------------------------------------------------------------------
   SUBROUTINE p_gather_int_1d1d (sendbuf, recvbuf, p_dest, comm)
     INTEGER,           INTENT(inout) :: sendbuf(:), recvbuf(:)
     INTEGER,           INTENT(in) :: p_dest
     INTEGER, OPTIONAL, INTENT(in) :: comm

#ifndef NOMPI
     CHARACTER(*), PARAMETER :: routine = TRIM("mo_mpi:p_gather_int_1d1d")
     INTEGER :: p_comm

     IF (PRESENT(comm)) THEN
       p_comm = comm
     ELSE
       p_comm = process_mpi_all_comm
     ENDIF

     CALL MPI_GATHER(sendbuf, SIZE(sendbuf), MPI_INTEGER, &
       &             recvbuf, SIZE(sendbuf), MPI_INTEGER, &
       &             p_dest, p_comm, p_error)
     IF (p_error /=  MPI_SUCCESS) CALL finish (routine, 'Error in MPI_GATHER operation!')
#else
     recvbuf = sendbuf
#endif
   END SUBROUTINE p_gather_int_1d1d


  !---------------------------------------------------------------------------------------------------------------------------------
  !> wrapper for MPI_Gather()
  !---------------------------------------------------------------------------------------------------------------------------------
   SUBROUTINE p_gather_char_0d1d (sbuf, recvbuf, p_dest, comm)
     CHARACTER(len=*),  INTENT(in)    ::  sbuf
     CHARACTER(len=*),  INTENT(inout) ::  recvbuf(:)
     INTEGER,           INTENT(in) :: p_dest
     INTEGER, OPTIONAL, INTENT(in) :: comm
     CHARACTER(*), PARAMETER :: routine = "mo_mpi:p_gather_char_0d1d"

#ifndef NOMPI
     INTEGER :: p_comm

     IF (PRESENT(comm)) THEN
       p_comm = comm
     ELSE
       p_comm = process_mpi_all_comm
     ENDIF

     ! recvbuf argument is only significant on root
     IF (p_comm_rank(p_comm) == p_dest) THEN
       IF (LEN(sbuf) /= LEN(recvbuf(1))) THEN
         CALL finish (routine, 'Internal error: String lengths do not match!')
       END IF
     END IF

     CALL MPI_GATHER(sbuf, LEN(sbuf), p_char,    &
       &             recvbuf, LEN(sbuf), p_char, &
       &             p_dest, p_comm, p_error)
     IF (p_error /=  MPI_SUCCESS) CALL finish (routine, 'Error in MPI_GATHER operation!')
#else
     IF (LEN(sbuf) /= LEN(recvbuf(1))) THEN
       CALL finish (routine, 'Internal error: String lengths do not match!')
     END IF
     recvbuf = sbuf
#endif
   END SUBROUTINE p_gather_char_0d1d


   SUBROUTINE p_gatherv_int (sendbuf, sendcount, recvbuf, recvcounts, &
     &                       displs, p_dest, comm)
     INTEGER,           INTENT(in)    :: sendbuf(:), sendcount
     INTEGER,           INTENT(inout) :: recvbuf(:), recvcounts(:), displs(:)
     INTEGER,           INTENT(in)    :: p_dest
     INTEGER, OPTIONAL, INTENT(in)    :: comm

#ifndef NOMPI
     CHARACTER(*), PARAMETER :: routine = TRIM("mo_mpi:p_gatherv_int")
     INTEGER :: p_comm

     IF (PRESENT(comm)) THEN
       p_comm = comm
     ELSE
       p_comm = process_mpi_all_comm
     ENDIF

     CALL MPI_GATHERV(sendbuf, sendcount,  MPI_INTEGER, &
       &              recvbuf, recvcounts, displs, MPI_INTEGER, &
       &              p_dest, p_comm, p_error)
     IF (p_error /=  MPI_SUCCESS) CALL finish (routine, 'Error in MPI_GATHERV operation!')
#else
     recvbuf((displs(1)+1):(displs(1)+sendcount)) = sendbuf(1:sendcount)
#endif
   END SUBROUTINE p_gatherv_int


   SUBROUTINE p_gatherv_real2D2D (sendbuf, sendcount, recvbuf, recvcounts, &
     &                            displs, p_dest, comm)
     REAL(DP), INTENT(IN) :: sendbuf(:,:)
     INTEGER, INTENT(IN)  :: sendcount
     REAL(DP), INTENT(OUT) :: recvbuf(:,:)
     INTEGER, INTENT(IN)  :: recvcounts(:)
     INTEGER, INTENT(IN)  :: displs(:)
     INTEGER, INTENT(IN)  :: p_dest
     INTEGER, INTENT(IN)  :: comm

#ifndef NOMPI
     CHARACTER(*), PARAMETER :: routine = TRIM("mo_mpi:p_gatherv_real2D2D")

     INTEGER :: dim1_size

     dim1_size = SIZE(sendbuf, 1)

     CALL MPI_GATHERV(sendbuf, sendcount*dim1_size,  p_real_dp, &
       &              recvbuf, recvcounts(:)*dim1_size, displs*dim1_size, &
       &              p_real_dp, p_dest, comm, p_error)
     IF (p_error /=  MPI_SUCCESS) &
       CALL finish (routine, 'Error in MPI_GATHERV operation!')
#else
     recvbuf(:, (displs(1)+1):(displs(1)+sendcount)) = sendbuf(:, 1:sendcount)
#endif
   END SUBROUTINE p_gatherv_real2D2D


   SUBROUTINE p_gatherv_sreal2D2D (sendbuf, sendcount, recvbuf, recvcounts, &
     &                            displs, p_dest, comm)
     REAL(SP), INTENT(IN) :: sendbuf(:,:)
     INTEGER, INTENT(IN)  :: sendcount
     REAL(SP), INTENT(OUT) :: recvbuf(:,:)
     INTEGER, INTENT(IN)  :: recvcounts(:)
     INTEGER, INTENT(IN)  :: displs(:)
     INTEGER, INTENT(IN)  :: p_dest
     INTEGER, INTENT(IN)  :: comm

#ifndef NOMPI
     CHARACTER(*), PARAMETER :: routine = TRIM("mo_mpi:p_gatherv_sreal2D2D")

     INTEGER :: dim1_size

     dim1_size = SIZE(sendbuf, 1)

     CALL MPI_GATHERV(sendbuf, sendcount*dim1_size,  p_real_sp, &
       &              recvbuf, recvcounts(:)*dim1_size, displs*dim1_size, &
       &              p_real_sp, p_dest, comm, p_error)
     IF (p_error /=  MPI_SUCCESS) &
       CALL finish (routine, 'Error in MPI_GATHERV operation!')
#else
     recvbuf(:, (displs(1)+1):(displs(1)+sendcount)) = sendbuf(:, 1:sendcount)
#endif
   END SUBROUTINE p_gatherv_sreal2D2D


   SUBROUTINE p_gatherv_int2D2D (sendbuf, sendcount, recvbuf, recvcounts, &
     &                            displs, p_dest, comm)
     INTEGER, INTENT(IN)  :: sendbuf(:,:)
     INTEGER, INTENT(IN)  :: sendcount
     INTEGER, INTENT(OUT)  :: recvbuf(:,:)
     INTEGER, INTENT(IN)  :: recvcounts(:)
     INTEGER, INTENT(IN)  :: displs(:)
     INTEGER, INTENT(IN)  :: p_dest
     INTEGER, INTENT(IN)  :: comm

#ifndef NOMPI
     CHARACTER(*), PARAMETER :: routine = TRIM("mo_mpi:p_gatherv_int2D2D")

     INTEGER :: dim1_size

     dim1_size = SIZE(sendbuf, 1)

     CALL MPI_GATHERV(sendbuf, sendcount*dim1_size,  p_int, &
       &              recvbuf, recvcounts(:)*dim1_size, displs*dim1_size, &
       &              p_int, p_dest, comm, p_error)
     IF (p_error /=  MPI_SUCCESS) &
       CALL finish (routine, 'Error in MPI_GATHERV operation!')
#else
     recvbuf(:, (displs(1)+1):(displs(1)+sendcount)) = sendbuf(:, 1:sendcount)
#endif
   END SUBROUTINE p_gatherv_int2D2D


   SUBROUTINE p_gatherv_real2D1D (sendbuf, sendcount, recvbuf, recvcounts, displs, p_dest, comm)
     REAL(dp),          INTENT(IN)    :: sendbuf(:,:)
     INTEGER,           INTENT(IN)    :: sendcount
     REAL(dp),          INTENT(INOUT) :: recvbuf(:)
     INTEGER,           intent(IN)    :: recvcounts(:), displs(:)
     INTEGER,           INTENT(in)    :: p_dest
     INTEGER,           INTENT(in)    :: comm

#if !defined(NOMPI)
     CHARACTER(*), PARAMETER :: routine = TRIM("mo_mpi:p_gatherv_real2D1D")
     INTEGER :: p_error

     CALL MPI_GATHERV(sendbuf, sendcount, p_real_dp,   &    ! sendbuf, sendcount, sendtype
       &              recvbuf, recvcounts, displs,     &    ! recvbuf, recvcounts, displs
       &              p_real_dp, p_dest, p_comm_work, p_error)  ! recvtype, root, comm, error
     IF (p_error /=  MPI_SUCCESS) CALL finish (routine, 'Error in MPI_GATHERV operation!')
#else
     recvbuf(:) = RESHAPE(sendbuf, (/ SIZE(recvbuf) /) )
#endif
   END SUBROUTINE p_gatherv_real2D1D


  SUBROUTINE p_gatherv_int2D1D (sendbuf, sendcount, recvbuf, recvcounts, displs, p_dest, comm)
    INTEGER,           INTENT(IN)    :: sendbuf(:,:)
    INTEGER,           INTENT(IN)    :: sendcount
    INTEGER,           INTENT(INOUT) :: recvbuf(:)
    INTEGER,           INTENT(IN)    :: recvcounts(:), displs(:)
    INTEGER,           INTENT(in)    :: p_dest
    INTEGER,           INTENT(in)    :: comm

#if !defined(NOMPI)
    CHARACTER(*), PARAMETER :: routine = TRIM("mo_mpi:p_gatherv_int2D1D")
    INTEGER :: p_error

    CALL MPI_GATHERV(sendbuf, sendcount, p_int,       &    ! sendbuf, sendcount, sendtype
      &              recvbuf, recvcounts, displs,     &    ! recvbuf, recvcounts, displs
      &              p_int, p_dest, p_comm_work, p_error)  ! recvtype, root, comm, error
    IF (p_error /=  MPI_SUCCESS) CALL finish (routine, 'Error in MPI_GATHERV operation!')
#else
    recvbuf(:) = RESHAPE(sendbuf, (/ SIZE(recvbuf) /) )
#endif
  END SUBROUTINE p_gatherv_int2D1D


   SUBROUTINE p_gatherv_real3D1D (sendbuf, sendcount, recvbuf, recvcounts, displs, p_dest, comm)
     REAL(dp),          INTENT(IN)    :: sendbuf(:,:,:)
     INTEGER,           INTENT(IN)    :: sendcount
     REAL(dp),          INTENT(INOUT) :: recvbuf(:)
     INTEGER,           intent(IN)    :: recvcounts(:), displs(:)
     INTEGER,           INTENT(in)    :: p_dest
     INTEGER,           INTENT(in)    :: comm

#if !defined(NOMPI)
     CHARACTER(*), PARAMETER :: routine = TRIM("mo_mpi:p_gatherv_real2D1D")
     INTEGER :: p_error

     CALL MPI_GATHERV(sendbuf, sendcount, p_real_dp,   &    ! sendbuf, sendcount, sendtype
       &              recvbuf, recvcounts, displs,     &    ! recvbuf, recvcounts, displs
       &              p_real_dp, p_dest, p_comm_work, p_error)  ! recvtype, root, comm, error
     IF (p_error /=  MPI_SUCCESS) CALL finish (routine, 'Error in MPI_GATHERV operation!')
#else
     recvbuf(:) = RESHAPE(sendbuf, (/ SIZE(recvbuf) /) )
#endif
   END SUBROUTINE p_gatherv_real3D1D


   SUBROUTINE p_scatterv_real1D2D (sendbuf, sendcounts, displs, recvbuf, recvcount, p_dest, comm)
     REAL(dp),          INTENT(IN)    :: sendbuf(:)
     INTEGER,           INTENT(IN)    :: sendcounts(:), displs(:)
     REAL(dp),          INTENT(INOUT) :: recvbuf(:,:)
     INTEGER,           INTENT(IN)    :: recvcount
     INTEGER,           INTENT(in)    :: p_dest
     INTEGER,           INTENT(in)    :: comm

#if !defined(NOMPI)
     CHARACTER(*), PARAMETER :: routine = TRIM("mo_mpi:p_scatterv_real1D2D")
     INTEGER :: p_error

     ! FIXME: I may be wrong, but this looks like a bug to me:
     !        This call *ignores* the communicator that is passed in and uses p_comm_work instead.
     CALL MPI_SCATTERV(sendbuf, sendcounts, displs,   &           ! sendbuf, sendcount, displs
       &               p_real_dp, recvbuf, recvcount, &           ! sendtype, recvbuf, recvcounts,
       &               p_real_dp, p_dest, p_comm_work, p_error)   ! recvtype, root, comm, error
     IF (p_error /=  MPI_SUCCESS) CALL finish (routine, 'Error in MPI_SCATTERV operation!')
#else
     recvbuf(:,:) = RESHAPE(sendbuf, (/ SIZE(recvbuf,1), SIZE(recvbuf,2) /))
#endif
   END SUBROUTINE p_scatterv_real1D2D


  !---------------------------------------------------------------------------------------------------------------------------------
  !> wrapper for MPI_Scatterv()
  !---------------------------------------------------------------------------------------------------------------------------------
   SUBROUTINE p_scatterv_real1D1D (sendbuf, sendcounts, displs, recvbuf, recvcount, p_src, comm)
        implicit none
        REAL(wp), INTENT(IN) :: sendbuf(:)
        INTEGER, INTENT(IN)  :: sendcounts(:)
        REAL(wp), INTENT(INOUT) :: recvbuf(:)
        INTEGER, INTENT(IN)  :: recvcount
        INTEGER, INTENT(IN)  :: displs(:)
        INTEGER, INTENT(IN)  :: p_src
        INTEGER, INTENT(IN)  :: comm

#ifndef NOMPI
        CHARACTER(*), PARAMETER :: routine = TRIM("mo_mpi:p_scatterv_real1D1D")
        INTEGER :: ierr

        CALL MPI_Scatterv(sendbuf, sendcounts, displs, p_real_dp, &
        &                 recvbuf, recvcount, p_real_dp, &
        &                 p_src, comm, ierr)
        IF (ierr /=  MPI_SUCCESS) CALL finish (routine, 'Error in MPI_Scatterv operation!')
#else
        recvbuf(1:recvcount) = sendbuf((displs(1)+1):(displs(1)+recvcount))
#endif
   END SUBROUTINE p_scatterv_real1D1D


  !---------------------------------------------------------------------------------------------------------------------------------
  !> wrapper for MPI_Scatterv()
  !---------------------------------------------------------------------------------------------------------------------------------
   SUBROUTINE p_scatterv_single1D1D (sendbuf, sendcounts, displs, recvbuf, recvcount, p_src, comm)
        implicit none
        REAL(sp), INTENT(IN) :: sendbuf(:)
        INTEGER, INTENT(IN)  :: sendcounts(:)
        REAL(sp), INTENT(INOUT) :: recvbuf(:)
        INTEGER, INTENT(IN)  :: recvcount
        INTEGER, INTENT(IN)  :: displs(:)
        INTEGER, INTENT(IN)  :: p_src
        INTEGER, INTENT(IN)  :: comm

#ifndef NOMPI
        CHARACTER(*), PARAMETER :: routine = TRIM("mo_mpi:p_scatterv_single1D1D")
        INTEGER :: ierr

        CALL MPI_Scatterv(sendbuf, sendcounts, displs, p_real_sp, &
        &                 recvbuf, recvcount, p_real_sp, &
        &                 p_src, comm, ierr)
        IF (ierr /=  MPI_SUCCESS) CALL finish (routine, 'Error in MPI_Scatterv operation!')
#else
        recvbuf(1:recvcount) = sendbuf((displs(1)+1):(displs(1)+recvcount))
#endif
   END SUBROUTINE p_scatterv_single1D1D

   SUBROUTINE p_allgather_int_0d1d(sendbuf, recvbuf, comm)
     INTEGER,           INTENT(inout) :: recvbuf(:)
     INTEGER,           INTENT(in) :: sendbuf
     INTEGER, OPTIONAL, INTENT(in) :: comm

#ifndef NOMPI
     CHARACTER(*), PARAMETER :: routine = TRIM("mo_mpi:p_allgather_int_0d1d")
     INTEGER :: p_comm

     IF (PRESENT(comm)) THEN
       p_comm = comm
     ELSE
       p_comm = process_mpi_all_comm
     ENDIF

     CALL mpi_allgather(sendbuf, 1, mpi_integer, &
          &             recvbuf, 1, mpi_integer, &
          &             p_comm, p_error)
     IF (p_error /=  MPI_SUCCESS) CALL finish (routine, 'Error in mpi_allgather operation!')
#else
     recvbuf = sendbuf
#endif
   END SUBROUTINE p_allgather_int_0d1d

   SUBROUTINE p_allgather_int_1d2d(sendbuf, recvbuf, comm)
     INTEGER,           INTENT(inout) :: recvbuf(:,:)
     INTEGER,           INTENT(in) :: sendbuf(:)
     INTEGER, OPTIONAL, INTENT(in) :: comm

#ifndef NOMPI
     CHARACTER(*), PARAMETER :: routine = TRIM("mo_mpi:p_allgather_int_0d1d")
     INTEGER :: p_comm, n

     IF (PRESENT(comm)) THEN
       p_comm = comm
     ELSE
       p_comm = process_mpi_all_comm
     ENDIF
     n = SIZE(sendbuf)
     CALL mpi_allgather(sendbuf, n, mpi_integer, &
          &             recvbuf, n, mpi_integer, &
          &             p_comm, p_error)
     IF (p_error /=  MPI_SUCCESS) CALL finish (routine, 'Error in mpi_allgather operation!')
#else
     recvbuf(:, 1) = sendbuf
#endif
   END SUBROUTINE p_allgather_int_1d2d

   SUBROUTINE p_allgatherv_real_1d(sendbuf, recvbuf, recvcounts, comm)
     REAL(dp),          INTENT(in)    :: sendbuf(:)
     REAL(dp),          INTENT(inout) :: recvbuf(:)
     INTEGER,           INTENT(in)    :: recvcounts(:)
     INTEGER, OPTIONAL, INTENT(in)    :: comm

#ifndef NOMPI
     CHARACTER(*), PARAMETER :: routine = TRIM("mo_mpi:p_allgatherv_real_1d")
     INTEGER :: p_comm, sendcount, comm_size, i
     INTEGER, ALLOCATABLE :: displs(:)

     IF (PRESENT(comm)) THEN
       p_comm = comm
     ELSE
       p_comm = process_mpi_all_comm
     ENDIF

     IF (p_comm_is_intercomm(p_comm)) THEN
      comm_size = p_comm_remote_size(p_comm)
     ELSE
      comm_size = p_comm_size(p_comm)
     END IF

     IF ((comm_size > SIZE(recvcounts, 1)) .OR. &
      &  (SUM(recvcounts) > SIZE(recvbuf, 1))) &
       CALL finish(routine, "invalid recvcounts")

     ALLOCATE(displs(comm_size))
     displs(1) = 0
     DO i = 2, comm_size
       displs(i) = displs(i-1) + recvcounts(i-1)
     END DO

     sendcount = SIZE(sendbuf)
     CALL mpi_allgatherv(sendbuf, sendcount, p_real_dp, &
          &              recvbuf, recvcounts, displs, p_real_dp, &
          &              p_comm, p_error)
     IF (p_error /=  MPI_SUCCESS) &
       CALL finish (routine, 'Error in mpi_allgatherv operation!')

     DEALLOCATE(displs)
#else
     recvbuf = sendbuf
#endif
   END SUBROUTINE p_allgatherv_real_1d

   SUBROUTINE p_allgatherv_int_1d(sendbuf, recvbuf, recvcounts, comm)
     INTEGER,           INTENT(in)    :: sendbuf(:)
     INTEGER,           INTENT(inout) :: recvbuf(:)
     INTEGER,           INTENT(in)    :: recvcounts(:)
     INTEGER, OPTIONAL, INTENT(in)    :: comm

#ifndef NOMPI
     CHARACTER(*), PARAMETER :: routine = TRIM("mo_mpi:p_allgatherv_int_1d")
     INTEGER :: p_comm, sendcount, comm_size, i
     INTEGER, ALLOCATABLE :: displs(:)

     IF (PRESENT(comm)) THEN
       p_comm = comm
     ELSE
       p_comm = process_mpi_all_comm
     ENDIF

     IF (p_comm_is_intercomm(p_comm)) THEN
      comm_size = p_comm_remote_size(p_comm)
     ELSE
      comm_size = p_comm_size(p_comm)
     END IF

     IF ((comm_size > SIZE(recvcounts, 1)) .OR. &
      &  (SUM(recvcounts) > SIZE(recvbuf, 1))) &
       CALL finish(routine, "invalid recvcounts")

     ALLOCATE(displs(comm_size))
     displs(1) = 0
     DO i = 2, comm_size
       displs(i) = displs(i-1) + recvcounts(i-1)
     END DO

     sendcount = SIZE(sendbuf)
     CALL mpi_allgatherv(sendbuf, sendcount, mpi_integer, &
          &              recvbuf, recvcounts, displs, mpi_integer, &
          &              p_comm, p_error)
     IF (p_error /=  MPI_SUCCESS) &
       CALL finish (routine, 'Error in mpi_allgatherv operation!')

     DEALLOCATE(displs)
#else
     recvbuf = sendbuf
#endif
   END SUBROUTINE p_allgatherv_int_1d


   ! Commits a user-defined MPI type
   !
   FUNCTION p_commit_type_struct(oldtypes, blockcounts) RESULT(newtype)
     INTEGER :: newtype
     INTEGER, INTENT(IN) :: oldtypes(2), blockcounts(2)
#if !defined(NOMPI)
     INTEGER :: ierr, offsets(2), extent
     ! define structured type and commit it
     CALL MPI_TYPE_EXTENT(oldtypes(1), extent, ierr)
     offsets(:) = (/ 0, extent /)
     CALL MPI_TYPE_STRUCT(2, blockcounts, offsets, oldtypes, newtype, ierr)
     CALL MPI_TYPE_COMMIT(newtype, ierr)
#else
     newtype = 0
#endif
   END FUNCTION p_commit_type_struct


   SUBROUTINE p_alltoall_int (sendbuf, recvbuf, comm)
     INTEGER,           INTENT(inout) :: sendbuf(:), recvbuf(:)
     INTEGER,           INTENT(in) :: comm
#if !defined(NOMPI)
     CHARACTER(*), PARAMETER :: routine = TRIM("mo_mpi:p_alltoall_int")
     INTEGER :: p_comm, p_error

     p_comm = comm
     CALL MPI_ALLTOALL(sendbuf, 1, p_int, recvbuf, 1, p_int, p_comm, p_error)
     IF (p_error /=  MPI_SUCCESS) CALL finish (routine, 'Error in MPI_ALLTOALL operation!')
#else
     recvbuf(:) = sendbuf(:)
#endif
   END SUBROUTINE p_alltoall_int


   SUBROUTINE p_alltoallv_real_2d (sendbuf, sendcounts, sdispls, &
     &                             recvbuf, recvcounts, rdispls, comm)
     REAL(dp),          INTENT(in) :: sendbuf(:,:)
     INTEGER,           INTENT(in) :: sendcounts(:), sdispls(:)
     REAL(dp),          INTENT(inout) :: recvbuf(:,:)
     INTEGER,           INTENT(in) :: recvcounts(:), rdispls(:)
     INTEGER,           INTENT(in) :: comm
#if !defined(NOMPI)
     CHARACTER(*), PARAMETER :: routine = TRIM("mo_mpi:p_alltoallv_real_2d")
     INTEGER :: p_comm, p_error, dim1_size

     p_comm = comm
     dim1_size = SIZE(sendbuf, 1)
     CALL MPI_ALLTOALLV(sendbuf, sendcounts(:)*dim1_size, &
       &                sdispls(:)*dim1_size, p_real_dp, recvbuf, &
       &                recvcounts(:)*dim1_size, rdispls(:)*dim1_size, &
       &                p_real_dp, p_comm, p_error)
     IF (p_error /=  MPI_SUCCESS) &
       CALL finish (routine, 'Error in MPI_ALLTOALLV operation!')
#else
     ! displs are zero based -> have to add 1
     recvbuf(:,rdispls(1)+1:rdispls(1)+recvcounts(1)) = &
       sendbuf(:,sdispls(1)+1:sdispls(1)+sendcounts(1))
#endif
   END SUBROUTINE p_alltoallv_real_2d


   SUBROUTINE p_alltoallv_sreal_2d (sendbuf, sendcounts, sdispls, &
     &                              recvbuf, recvcounts, rdispls, comm)
     REAL(sp),          INTENT(in) :: sendbuf(:,:)
     INTEGER,           INTENT(in) :: sendcounts(:), sdispls(:)
     REAL(sp),          INTENT(inout) :: recvbuf(:,:)
     INTEGER,           INTENT(in) :: recvcounts(:), rdispls(:)
     INTEGER,           INTENT(in) :: comm
#if !defined(NOMPI)
     CHARACTER(*), PARAMETER :: routine = TRIM("mo_mpi:p_alltoallv_real_2d")
     INTEGER :: p_comm, p_error, dim1_size

     p_comm = comm
     dim1_size = SIZE(sendbuf, 1)
     CALL MPI_ALLTOALLV(sendbuf, sendcounts(:)*dim1_size, &
       &                sdispls(:)*dim1_size, p_real_sp, recvbuf, &
       &                recvcounts(:)*dim1_size, rdispls(:)*dim1_size, &
       &                p_real_sp, p_comm, p_error)
     IF (p_error /=  MPI_SUCCESS) &
       CALL finish (routine, 'Error in MPI_ALLTOALLV operation!')
#else
     ! displs are zero based -> have to add 1
     recvbuf(:,rdispls(1)+1:rdispls(1)+recvcounts(1)) = &
       sendbuf(:,sdispls(1)+1:sdispls(1)+sendcounts(1))
#endif
   END SUBROUTINE p_alltoallv_sreal_2d


   SUBROUTINE p_alltoallv_int_2d (sendbuf, sendcounts, sdispls, &
     &                            recvbuf, recvcounts, rdispls, comm)
     INTEGER,           INTENT(in) :: sendbuf(:,:)
     INTEGER,           INTENT(in) :: sendcounts(:), sdispls(:)
     INTEGER,           INTENT(inout) :: recvbuf(:,:)
     INTEGER,           INTENT(in) :: recvcounts(:), rdispls(:)
     INTEGER,           INTENT(in) :: comm
#if !defined(NOMPI)
     CHARACTER(*), PARAMETER :: routine = TRIM("mo_mpi:p_alltoallv_int_2d")
     INTEGER :: p_comm, p_error, dim1_size

     p_comm = comm
     dim1_size = SIZE(sendbuf, 1)
     CALL MPI_ALLTOALLV(sendbuf, sendcounts(:)*dim1_size, &
       &                sdispls(:)*dim1_size, p_int, recvbuf, &
       &                recvcounts(:)*dim1_size, rdispls(:)*dim1_size, &
       &                p_int, p_comm, p_error)
     IF (p_error /=  MPI_SUCCESS) &
       CALL finish (routine, 'Error in MPI_ALLTOALLV operation!')
#else
     ! displs are zero based -> have to add 1
     recvbuf(:,rdispls(1)+1:rdispls(1)+recvcounts(1)) = &
       sendbuf(:,sdispls(1)+1:sdispls(1)+sendcounts(1))
#endif
   END SUBROUTINE p_alltoallv_int_2d


   SUBROUTINE p_alltoallv_int (sendbuf, sendcounts, sdispls, &
     &                         recvbuf, recvcounts, rdispls, comm)
     INTEGER,           INTENT(in) :: sendbuf(:), sendcounts(:), sdispls(:)
     INTEGER,           INTENT(inout) :: recvbuf(:)
     INTEGER,           INTENT(in) :: recvcounts(:), rdispls(:)
     INTEGER,           INTENT(in) :: comm
#if !defined(NOMPI)
     CHARACTER(*), PARAMETER :: routine = TRIM("mo_mpi:p_alltoallv_int")
     INTEGER :: p_comm, p_error

     p_comm = comm
     CALL MPI_ALLTOALLV(sendbuf, sendcounts, sdispls, p_int, &
       &                recvbuf, recvcounts, rdispls, p_int, p_comm, p_error)
     IF (p_error /=  MPI_SUCCESS) &
       CALL finish (routine, 'Error in MPI_ALLTOALLV operation!')
#else
     ! displs are zero based -> have to add 1
     recvbuf(rdispls(1)+1:rdispls(1)+recvcounts(1)) = &
       sendbuf(sdispls(1)+1:sdispls(1)+sendcounts(1))
#endif
   END SUBROUTINE p_alltoallv_int

#if !defined(NOMPI)
   SUBROUTINE p_alltoallv_p2p_real_2d_core(dim1_size, sendbuf, sendcounts, sdispls, &
        &                             recvbuf, recvcounts, rdispls, comm)
     INTEGER, INTENT(in) :: dim1_size
     REAL(wp),          INTENT(in) :: sendbuf(dim1_size,*)
     INTEGER,           INTENT(in) :: sendcounts(:), sdispls(:)
     REAL(wp),          INTENT(inout) :: recvbuf(dim1_size,*)
     INTEGER,           INTENT(in) :: recvcounts(:), rdispls(:)
     INTEGER,           INTENT(in) :: comm
     CHARACTER(*), PARAMETER :: routine = "mo_mpi:p_alltoallv_p2p_real_2d"
     INTEGER :: i, comm_size, tag, ofs, datatype

     CALL p_wait

     comm_size = p_comm_size(comm)
     tag = 1
     SELECT CASE (wp)
       CASE(sp)
         datatype = p_real_sp
       CASE(dp)
         datatype = p_real_dp
       CASE DEFAULT
         datatype = -1
         CALL finish (routine, 'invalid read type')
     END SELECT
     DO i = 1, comm_size
       IF (recvcounts(i) > 0) THEN
         ofs = 1 + rdispls(i)
         CALL p_inc_request
         CALL mpi_irecv(recvbuf(:, ofs:ofs+recvcounts(i)-1), &
              &         recvcounts(i) * dim1_size, datatype, i - 1, tag, &
              &         comm, p_request(p_irequest), p_error)
#ifdef DEBUG
         IF (p_error /= MPI_SUCCESS) THEN
           WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_IRECV on ', &
                my_process_mpi_all_id, &
                ' from ', i - 1, ' for tag ', tag, ' failed.'
           WRITE (nerr,'(a,i4)') ' Error = ', p_error
           CALL abort_mpi
         END IF
#endif
       END IF
       IF (sendcounts(i) > 0) THEN
         ofs = 1 + sdispls(i)
         CALL p_inc_request
         CALL mpi_isend(sendbuf(:, ofs:ofs+sendcounts(i)-1), &
              &         sendcounts(i) * dim1_size, datatype, i - 1, tag, &
              &         comm, p_request(p_irequest), p_error)
#ifdef DEBUG
         IF (p_error /= MPI_SUCCESS) THEN
           WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_ISEND on ', &
                my_process_mpi_all_id, &
                ' from ', i - 1, ' for tag ', tag, ' failed.'
           WRITE (nerr,'(a,i4)') ' Error = ', p_error
           CALL abort_mpi
         END IF
#endif
       END IF
     END DO
     CALL p_wait
   END SUBROUTINE p_alltoallv_p2p_real_2d_core
#endif

   SUBROUTINE p_alltoallv_p2p_real_2d (sendbuf, sendcounts, sdispls, &
     &                             recvbuf, recvcounts, rdispls, comm)
     REAL(wp),          INTENT(in) :: sendbuf(:,:)
     INTEGER,           INTENT(in) :: sendcounts(:), sdispls(:)
     REAL(wp),          INTENT(inout) :: recvbuf(:,:)
     INTEGER,           INTENT(in) :: recvcounts(:), rdispls(:)
     INTEGER,           INTENT(in) :: comm
#if !defined(NOMPI)
     CHARACTER(*), PARAMETER :: routine = "mo_mpi:p_alltoallv_p2p_real_2d"
     INTEGER :: dim1_size

     dim1_size = SIZE(sendbuf, 1)
     CALL p_alltoallv_p2p_real_2d_core(dim1_size, &
          &                            sendbuf, sendcounts, sdispls, &
          &                            recvbuf, recvcounts, rdispls, comm)
#else
     ! displs are zero based -> have to add 1
     recvbuf(:,rdispls(1)+1:rdispls(1)+recvcounts(1)) = &
       sendbuf(:,sdispls(1)+1:sdispls(1)+sendcounts(1))
#endif
   END SUBROUTINE p_alltoallv_p2p_real_2d

#if !defined(NOMPI)
   SUBROUTINE p_alltoallv_p2p_int_2d_core(dim1_size, sendbuf, sendcounts, sdispls, &
        &                                 recvbuf, recvcounts, rdispls, comm)
     INTEGER, INTENT(in) :: dim1_size
     INTEGER,           INTENT(in) :: sendbuf(dim1_size,*)
     INTEGER,           INTENT(in) :: sendcounts(:), sdispls(:)
     INTEGER,           INTENT(inout) :: recvbuf(dim1_size,*)
     INTEGER,           INTENT(in) :: recvcounts(:), rdispls(:)
     INTEGER,           INTENT(in) :: comm
     CHARACTER(*), PARAMETER :: routine = "mo_mpi:p_alltoallv_p2p_int_2d"
     INTEGER :: i, comm_size, tag, ofs

     CALL p_wait

     comm_size = p_comm_size(comm)
     tag = 1
     DO i = 1, comm_size
       IF (recvcounts(i) > 0) THEN
         ofs = 1 + rdispls(i)
         CALL p_inc_request
         CALL mpi_irecv(recvbuf(:, ofs:ofs+recvcounts(i)-1), &
              &         recvcounts(i) * dim1_size, p_int, i - 1, tag, &
              &         comm, p_request(p_irequest), p_error)
#ifdef DEBUG
         IF (p_error /= MPI_SUCCESS) THEN
           WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_IRECV on ', &
                my_process_mpi_all_id, &
                ' from ', i - 1, ' for tag ', tag, ' failed.'
           WRITE (nerr,'(a,i4)') ' Error = ', p_error
           CALL abort_mpi
         END IF
#endif
       END IF
       IF (sendcounts(i) > 0) THEN
         ofs = 1 + sdispls(i)
         CALL p_inc_request
         CALL mpi_isend(sendbuf(:, ofs:ofs+sendcounts(i)-1), &
              &         sendcounts(i) * dim1_size, p_int, i - 1, tag, &
              &         comm, p_request(p_irequest), p_error)
#ifdef DEBUG
         IF (p_error /= MPI_SUCCESS) THEN
           WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_ISEND on ', &
                my_process_mpi_all_id, &
                ' from ', i - 1, ' for tag ', tag, ' failed.'
           WRITE (nerr,'(a,i4)') ' Error = ', p_error
           CALL abort_mpi
         END IF
#endif
       END IF
     END DO
     CALL p_wait
   END SUBROUTINE p_alltoallv_p2p_int_2d_core
#endif

   SUBROUTINE p_alltoallv_p2p_int_2d (sendbuf, sendcounts, sdispls, &
     &                                recvbuf, recvcounts, rdispls, comm)
     INTEGER,           INTENT(in) :: sendbuf(:,:)
     INTEGER,           INTENT(in) :: sendcounts(:), sdispls(:)
     INTEGER,           INTENT(inout) :: recvbuf(:,:)
     INTEGER,           INTENT(in) :: recvcounts(:), rdispls(:)
     INTEGER,           INTENT(in) :: comm
#if !defined(NOMPI)
     CHARACTER(*), PARAMETER :: routine = "mo_mpi:p_alltoallv_p2p_int_2d"
     INTEGER :: dim1_size

     dim1_size = SIZE(sendbuf, 1)
     CALL p_alltoallv_p2p_int_2d_core(dim1_size, &
          &                           sendbuf, sendcounts, sdispls, &
          &                           recvbuf, recvcounts, rdispls, comm)
#else
     ! displs are zero based -> have to add 1
     recvbuf(:,rdispls(1)+1:rdispls(1)+recvcounts(1)) = &
       sendbuf(:,sdispls(1)+1:sdispls(1)+sendcounts(1))
#endif
   END SUBROUTINE p_alltoallv_p2p_int_2d

  SUBROUTINE p_clear_request(request)
    INTEGER, INTENT(INOUT) :: request
#if !defined(NOMPI)
    request = MPI_REQUEST_NULL
#endif
  END SUBROUTINE p_clear_request


  FUNCTION p_mpi_wtime()
    REAL(dp) :: p_mpi_wtime
#ifndef NOMPI
    p_mpi_wtime = MPI_Wtime()
#else
    p_mpi_wtime = 0d0
#endif
  END FUNCTION p_mpi_wtime


  !--------------------------------------------------------------------


  !> @return Global MPI ranks within communicator "comm"
  !
  SUBROUTINE get_mpi_comm_world_ranks(comm, global_ranks, nranks)
    INTEGER, INTENT(IN)  :: comm               !< MPI communicator
    INTEGER, INTENT(OUT) :: global_ranks(:)    !< Output: list of global MPI ranks in communicator "comm"
    INTEGER, INTENT(OUT) :: nranks             !< Output: number of entries in rank list
    ! local variables
    CHARACTER(*), PARAMETER :: routine = TRIM("mo_mpi:get_mpi_comm_world_ranks")
    INTEGER              :: p_error, grp_comm, grp_comm_world, i
    INTEGER, ALLOCATABLE :: comm_ranks(:)

#if !defined(NOMPI)
    nranks = 0
    IF (comm /= MPI_COMM_NULL) THEN
      nranks = p_comm_size(comm)    ! inquire communicator size

      IF (nranks > SIZE(global_ranks))  CALL finish (routine, 'Input array too small!')
      ALLOCATE(comm_ranks(nranks))
      comm_ranks(1:nranks) = (/ (i, i=0,(nranks-1)) /)

      CALL MPI_COMM_GROUP(comm, grp_comm, p_error)
      IF (p_error /= MPI_SUCCESS)  CALL finish (routine, 'Error in MPI_COMM_GROUP operation!')
      CALL MPI_COMM_GROUP(MPI_COMM_WORLD, grp_comm_world,  p_error)
      IF (p_error /= MPI_SUCCESS)  CALL finish (routine, 'Error in MPI_COMM_GROUP operation!')

      global_ranks(:) = 0
      CALL MPI_GROUP_TRANSLATE_RANKS(grp_comm, nranks, comm_ranks, &
        &                            grp_comm_world, global_ranks, p_error)
      IF (p_error /= MPI_SUCCESS)  CALL finish (routine, 'Error in MPI_GROUP_TRANSLATE_RANKS operation!')

      CALL MPI_GROUP_FREE(grp_comm,       p_error)
      IF (p_error /= MPI_SUCCESS)  CALL finish (routine, 'Error in MPI_GROUP_FREE operation!')
      CALL MPI_GROUP_FREE(grp_comm_world, p_error)
      IF (p_error /= MPI_SUCCESS)  CALL finish (routine, 'Error in MPI_GROUP_FREE operation!')

      DEALLOCATE(comm_ranks)
    END IF
#else
    nranks = 1
    global_ranks(1) = 0
#endif
  END SUBROUTINE get_mpi_comm_world_ranks

  !--------------------------------------------------------------------

  LOGICAL FUNCTION p_comm_is_intercomm(intercomm)

    INTEGER, INTENT(IN) :: intercomm
    LOGICAL :: flag
    INTEGER :: p_error

#ifndef NOMPI
    CALL MPI_COMM_TEST_INTER(intercomm, flag, p_error)
    IF (p_error /= MPI_SUCCESS) &
      CALL finish ("p_comm_is_intercomm", &
        &          'Error in MPI_Comm_test_inter operation!')

    p_comm_is_intercomm = flag
#else
    p_comm_is_intercomm = .FALSE.
#endif
  END FUNCTION p_comm_is_intercomm

  !--------------------------------------------------------------------

  INTEGER FUNCTION p_comm_remote_size(intercomm)

    INTEGER, INTENT(IN) :: intercomm
    INTEGER :: remote_size, p_error

#ifndef NOMPI
    CALL MPI_COMM_REMOTE_SIZE(intercomm, remote_size, p_error)
    IF (p_error /= MPI_SUCCESS) &
      CALL finish ("p_comm_remote_size", &
        &          'Error in MPI_Comm_remote_size operation!')

    p_comm_remote_size = remote_size
#else
    p_comm_remote_size = 0
#endif

  END FUNCTION p_comm_remote_size


  LOGICAL FUNCTION p_isEqual_int(value, comm) RESULT(resultVar)
    INTEGER, VALUE :: value
    INTEGER, OPTIONAL, INTENT(IN) :: comm

#ifndef NOMPI
    CHARACTER(*), PARAMETER :: routine = modname//":p_isEqual_int"
    INTEGER :: sendBuffer(2), minmax(2)

    !Compute the MIN AND MAX of the values using the fact that MAX_i(x_i) = -min_i(-x_i)
    !i. e. we compute both with a single MPI MIN reduction.
    sendBuffer(1) = VALUE
    sendBuffer(2) = -VALUE
    minmax = p_min(sendBuffer, comm = comm)
    resultVar = minmax(1) == -minmax(2)
#else
    resultVar = .TRUE.
#endif
  END FUNCTION p_isEqual_int

  LOGICAL FUNCTION p_isEqual_charArray(charArray, comm) RESULT(resultVar)
    CHARACTER(KIND = C_CHAR), INTENT(IN) :: charArray(:)
    INTEGER, OPTIONAL, INTENT(IN) :: comm

#ifndef NOMPI
    INTEGER :: myProc, nextProc, prevProc, p_comm, commSize, i, error
    CHARACTER(KIND = C_CHAR) :: prevArray(SIZE(charArray, 1))
    LOGICAL :: stringsEqual

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = process_mpi_all_comm
    ENDIF

    !Check whether all processes have the same string SIZE.
    resultVar = p_isEqual(SIZE(charArray, 1), comm = p_comm)
    IF(.NOT. resultVar) RETURN

    !Get the neighbor ranks for a cyclic DATA exchange.
    myProc = p_comm_rank(p_comm)
    commSize = p_comm_size(p_comm)
    nextProc = myProc + 1
    IF(nextProc == commSize) nextProc = 0
    prevProc = myProc - 1
    IF(prevProc == -1) prevProc = commSize - 1

    !Do a cyclic exchange of the strings, compare the local string to the one passed IN from the prevProc, AND reduce whether all strings compared equal.
    CALL p_sendrecv_char_array(charArray, nextProc, prevArray, prevProc, 0, comm = p_comm)
    stringsEqual = .TRUE.
    DO i = 1, SIZE(charArray, 1)
        IF(charArray(i) /= prevArray(i)) stringsEqual = .FALSE.
    END DO

    CALL MPI_Allreduce(stringsEqual, resultVar, 1, MPI_LOGICAL, MPI_LAND, p_comm, error)
#else
    resultVar = .TRUE.
#endif
  END FUNCTION p_isEqual_charArray

END MODULE mo_mpi
