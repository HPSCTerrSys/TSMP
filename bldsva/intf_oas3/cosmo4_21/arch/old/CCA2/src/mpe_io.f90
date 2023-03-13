!+ Module for parallel I/O in meteorological codes

MODULE mpe_io

! Description:
!   This module provides routines and data structures for 
!   conducting parallel I/O in several modes and on different
!   devices (DWD file I/O, DWD database I/O). The design is
!   kept in layers to support the inclusion of new devices.
!
! Current Code Owner: Pallas GmbH, Heinrich Bockhorst 
!   phone : +49-(0)2232-1896-0
!   direct: +49-(0)2232-1896-27 
!   email : bockhorst@pallas.com
!
! History:          Date        Name
! ---------         ----        ----
! 1.00              1999/01/21  Heinrich Bockhorst
!  Initial release
! 1.10              1999/02/26  H. J. Plum
!                               Integration of namelist variables
!                               database/idbg_level, database/iretry
!                               stdout tidied up
!
!  
!
!
!----------------------------------------------------------------------------

! include statements
INCLUDE "mpif.h"

! public interfaces of mpe_io

PUBLIC  mpe_io_init          & ! initialization
       ,mpe_io_reconfig      & ! reconfiguration
       ,mpe_io_open          & ! opens a file
       ,mpe_io_write         & ! writes to file
       ,mpe_io_complete      & ! complete a (collective) mpe_io_write
       ,mpe_io_read          & ! read from file
       ,mpe_io_close         & ! close file on I/O PEs
       ,mpe_is_compute       & ! identifies compute PEs
       ,mpe_is_io            & ! identifies io      PEs
       ,mpe_io_node          & ! program for "dumb" io node
       ,mpe_io_shutdown      & ! send shutdown signal to io_nodes
       ,mpe_io_ready         & ! actions at the end of a timestep
       ,mpe_io_ready_wait    & ! wait if io accumulates 
       ,mpe_readnldb         & ! NAMELIST /database/ input
       ,mpe_db_tab           & ! edits 'ty=<table>' database order 
       ,mpe_get_db_tab         ! retrieves current db order



! private routines 

PRIVATE mpe_db_init          & ! database initialization 
      , char_pack              ! character to integer conversion
            
PRIVATE

! KIND -parameters for declaring the variables that are passed to the
! grib library:  values for the DWDLIB
INTEGER, parameter :: gribc_kind = kind(1)
INTEGER, parameter :: gribf_kind = kind(1)
!INTEGER, parameter :: gribf_kind = 4    ! only for using libgrib1 on T3E

! Communicator and process identification

INTEGER :: io_config          ! identify io configuration

INTEGER :: MPE_COMM_IO       &! io PE communicator
         , num_io            &! number of PEs in io communicator
         , me_io             &! rank in io communicator
         , MPE_COMM_COMPUTE  &! compute communicator
         , num_compute       &! number of compute PEs
         , me_compute        &! rank in io communicator
         , MPE_COMM_CO2IO     ! communicator between COMPUTE 
                              ! and IO PEs.  

! special tags used inside the mpe_io  system

INTEGER ,parameter :: inter_tag     =  7777   ! tag for 
                                              ! inter_comm_create
INTEGER ,parameter :: ready_tag     = 13472
INTEGER ,parameter :: tag_open      =  6666   ! tag for open 
INTEGER ,parameter :: tag_write_end = 22222   ! this tag is reserved 
                                              ! for the "end-message"
                                              ! of write 
INTEGER , parameter :: tag_read     = 12345   ! tag for read

! status flags for the mpe_io system

LOGICAL :: got_end_signal                    ! end_signal has been recv.

! Field description 

INTEGER :: field_snd      &! index of the current field to be sent
         , field_rcv       ! index of the current field to be received  

LOGICAL :: is_compute_pe     &! self explaining
         , is_io_pe

!--- Send buffer management -------------------------------------------

INTEGER, parameter                  :: send_max_request                = 3
INTEGER                             :: send_requests(send_max_request) = MPI_REQUEST_NULL
INTEGER                             :: send_current_request            = 1

INTEGER(KIND=gribf_kind), allocatable :: send_buffer(:,:)


!--- Buffer environment ------------------------------------------------

INTEGER , parameter    :: BUFF_length = 500000  ! value for dwdlib...

INTEGER (KIND = gribf_kind) &
        , allocatable  :: BUFF_space( :,:)     ! buffer for write


INTEGER                :: BUFF_max_requests        &! number of requests
                         ,BUFF_request              ! current request

INTEGER, allocatable   :: BUFF_len_curr(:)          ! current buffer_length

!--- Device specific stuff: ----------------------------------------------

INTEGER, parameter   :: dev_name_len = 20  ! length of device identifier

CHARACTER (LEN=dev_name_len) :: out_device ! output system == 'file'
                                           !                , 'bank', ...
CHARACTER (LEN=dev_name_len) ::  in_device ! input  sysyem == ...

CHARACTER (LEN=dev_name_len) ::  ori_out_device ! original device 
                                                ! in case of error

LOGICAL                      :: error_handling  ! flag to identify whether
                                                ! error handling is on

CHARACTER(LEN=1) :: current_mode  ! = 'r', 'w' : sysytem is in read/write
                               !              mode

INTEGER, parameter :: max_name_len  = 250     ! maximal length for file names

CHARACTER (LEN=max_name_len) :: current_datname ! backup of filename for 
                                              ! error handling

CHARACTER (LEN=3)          :: mode_backup     ! backup of exact mode
                                              ! is not equal to current_mode
                                              ! for mode = 'a'


LOGICAL     :: end_read  


!------- Database parameters ---------------------------------------------


CHARACTER(LEN=300) :: DB_order                        ! order to DB
CHARACTER(LEN=300) :: DB_init_order='ak=nix'          ! initial order to DB
CHARACTER(LEN=300) :: DB_output_order=' '             ! output order to DB
CHARACTER(LEN=300) :: DB_input_order=' '              ! input order to DB
CHARACTER(LEN=300) :: DB_ana_tab='*****'              ! DB table of analysis 
                                                      ! input data
CHARACTER(LEN=300) :: DB_bd_tab='*****'               ! DB table of boundary 
                                                      ! input data
INTEGER, parameter :: DB_MAX_SOCKETS = 100            ! max number of Sock.
INTEGER            :: DB_OUT_SOCKETS = 0              ! used sockets
INTEGER            :: DB_IN_SOCKETS = 0               ! used sockets
INTEGER            :: DB_retry = 300                 ! #seconds to re-try
                                                      ! on DB-failure
INTEGER            :: DB_request(DB_MAX_SOCKETS)      ! async. requests
INTEGER            :: DB_hits                         ! sucessfull == 1
INTEGER            :: DB_err                          ! error code

INTEGER            :: MPE_dbg_level = 0

LOGICAL            :: DB_read_done=.false.

!------- end Database para.  ---------------------------------------------


!------- Backup system ( GRIB_SAVE)  -------------------------------------

LOGICAL                       :: DB_do_backup=.false. ! backup en/disabled
INTEGER                       :: GRIB_SAVE_SIZE=-1    ! size of incore buffer
CHARACTER (LEN=max_name_len)  :: GRIB_SAVE_NAME='.'   ! backup directory

!------- status parameters for diagnosis ---------------------------------

INTEGER                       :: open_count = 0       ! counter for unclosed
                                                      ! calls to mpe_io_open 

!--------------- error testing --------------------------------------------

INTEGER ::  db_init_count = 0, db_req_count   = 0,&
            db_end_count = 0 , db_close_count = 0

INTEGER ::  db_init_err = -1, db_req_err   = -1,&
            db_end_err  = -1, db_close_err = -1 


CONTAINS

!--------------------------------------------------------------------------
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!--------------------------------------------------------------------------

SUBROUTINE mpe_readnldb (knuin, knuout, iostat)
!
!=======================================================================
!
!     *mpe_readnldb* reads the NAMELIST(s) /database/ from INPUT_GME
!      
!
!=======================================================================
!
!     Author: Pallas GmbH,  Jan. 1999
!
!
!
!=======================================================================
!
  IMPLICIT NONE

!
!=======================================================================
!
!     Input
!
  INTEGER, INTENT(IN)  ::  knuin       ! unit number for input  file 
  INTEGER, INTENT(IN)  ::  knuout      ! unit number for output file
  INTEGER, INTENT(OUT) ::  iostat      ! io status
!
!=======================================================================
! 
!-------------------
! NAMELIST variables
!-------------------
!
CHARACTER(LEN=300) :: yinit_order                  ! initial order to DB
CHARACTER(LEN=300) :: yana_tab                     ! DB table of analysis
                                                   ! input data
CHARACTER(LEN=300) :: ybd_tab                      ! DB table of boundary
                                                   ! input data
INTEGER            :: nout_sockets                 ! used sockets
INTEGER            :: nin_sockets                  ! used sockets
INTEGER            :: iretry                       ! #seconds to re-try
                                                   ! on DB-failure
CHARACTER(LEN=100) :: ybackup_dir                  ! Directory for out-core
INTEGER            :: ibackup_size                 ! Size of in-core backup 
                                                   ! buffer
INTEGER            :: idbg_level                   ! Controls debug 

!
!
!=======================================================================
! 
  NAMELIST /database/         &
       yinit_order,           &
       yana_tab,              &
       ybd_tab,               &
       nout_sockets,          &
       nin_sockets,           &
       iretry,                &
       ybackup_dir,           &
       ibackup_size,          &
       idbg_level             
!
!=======================================================================
!
!--------------------------------------------------------
!
! Local variables
!
  logical:: do_out

  do_out = .false.

  IF( knuout > 0 ) THEN
   INQUIRE(UNIT=knuout,opened=do_out)
  ENDIF

  iostat = 0

  yinit_order     = DB_init_order(1:LEN_TRIM(DB_init_order))
  yana_tab        = DB_ana_tab(1:LEN_TRIM(DB_ana_tab))
  ybd_tab         = DB_bd_tab(1:LEN_TRIM(DB_bd_tab))
  nout_sockets    = DB_OUT_SOCKETS
  nin_sockets     = DB_IN_SOCKETS 
  iretry          = DB_retry
  ybackup_dir     = GRIB_SAVE_NAME(1:LEN_TRIM(GRIB_SAVE_NAME))
  ibackup_size    = GRIB_SAVE_SIZE
  idbg_level      = MPE_dbg_level

  read(knuin,database,end=1110,err=1111)
  GOTO 1110
1111 continue

  WRITE(0,*) 'READNLDB (database NAMELIST) '
  WRITE(0,*) 'WARNING!! S.TH. WRONG WITH NAMELIST STRUCTURE'
  WRITE(0,*) 'DATABASE DISABLED '
  iostat = 1002

1110 continue

  DB_init_order  = yinit_order(1:LEN_TRIM(yinit_order))
  DB_ana_tab     = yana_tab(1:LEN_TRIM(yana_tab))
  DB_bd_tab      = ybd_tab(1:LEN_TRIM(ybd_tab))
  DB_OUT_SOCKETS = nout_sockets
  DB_IN_SOCKETS  = nin_sockets
  DB_retry       = iretry
  GRIB_SAVE_NAME = ybackup_dir(1:LEN_TRIM(ybackup_dir))
  GRIB_SAVE_SIZE = ibackup_size
  MPE_dbg_level  = idbg_level

  DB_do_backup = (GRIB_SAVE_SIZE>=0 .and. DB_OUT_SOCKETS > 0 )

  DB_read_done = .true.

! print content of namelists
!
  IF (do_out) THEN
    WRITE(knuout,'(//)')
    WRITE(knuout,'(A,/)') ' Namelist /database/'
 
    WRITE(knuout,*) ' init_order   ',DB_init_order(1:LEN_TRIM(DB_init_order))
    WRITE(knuout,*) ' ana_tab      ',DB_ana_tab(1:LEN_TRIM(DB_ana_tab))
    WRITE(knuout,*) ' bd_tab       ',DB_bd_tab(1:LEN_TRIM(DB_bd_tab))
    WRITE(knuout,*) ' OUT_SOCKETS  ',DB_OUT_SOCKETS
    WRITE(knuout,*) ' IN_SOCKETS   ',DB_IN_SOCKETS
    WRITE(knuout,*) ' Retry-time   ',DB_retry
    IF( DB_do_backup ) THEN
      WRITE(knuout,*) ' Backup mechanism, out-core directory: '
      WRITE(knuout,*) GRIB_SAVE_NAME(1:LEN_TRIM(GRIB_SAVE_NAME))
      WRITE(knuout,*) ' In-core size: '
      WRITE(knuout,*) GRIB_SAVE_SIZE
    ELSE
      WRITE(knuout,*) ' Backup mechanism disabled '
    ENDIF

  ENDIF

!=======================================================================
!
END SUBROUTINE mpe_readnldb


SUBROUTINE char_pack(ystring, istring, islen, idir, ierr)

IMPLICIT NONE

CHARACTER *(*), INTENT(INOUT):: ystring
INTEGER, INTENT(IN)          :: idir
INTEGER, INTENT(INOUT)       :: islen
INTEGER, INTENT(INOUT)       :: istring(islen), ierr

INTEGER:: slen,i

IF( idir == 1 ) THEN

 slen = LEN_TRIM(ystring)

 IF( slen > islen ) THEN
  WRITE(0,*) ' Erroneous use of routine char_pack, string truncated '
 ELSE
  islen = slen
 ENDIF

 DO i=1,islen
  istring(i) = ichar(ystring(i:i))
 ENDDO

ELSE

 slen = LEN(ystring)

 IF( slen < islen ) THEN
  WRITE(0,*) ' Erroneous use of routine char_pack, string truncated '
 ELSE
  slen = islen
 ENDIF

 DO i=1,slen
  ystring(i:i) = achar(istring(i))
 ENDDO

ENDIF

END SUBROUTINE char_pack

!USUS------------------------------------------------------------------------

SUBROUTINE mpe_db_init( )

!--------------------------------------------------------------------------
!
! $Revision: 1.1.1.6 $
! $DATE$
!
!
! Purpose      : Set the basic parameters
!
!                in_device, out_device
!
!                DB_init_order
!                DB_ana_tab/DB_bd_tab
!                DB_OUT_SOCKETS 
!                DB_IN_SOCKETS 
!                DB_retry
!                DB_do_backup
!
!                GRIB_SAVE_SIZE
!                GRIB_SAVE_NAME
!
!                This subroutine has to be called after setup of
!                the parameters
!                         num_io, is_io_pe
!                in -> mpe_io_init
!
!--------------------------------------------------------------------------

  IMPLICIT NONE

!--------------------------------------------------------------------------
!
INTEGER:: max_data
INTEGER, DIMENSION(:), ALLOCATABLE:: buf 

#ifdef __MPICH2
INTEGER, DIMENSION(:), ALLOCATABLE :: bsend_buf
INTEGER :: bsend_size
#endif

INTEGER:: iroot,root,pe,itag,ilen,ierr

INTEGER:: ndata, status(MPI_STATUS_SIZE)

!--------------------------------------------------------------------------

max_data=max(LEN(DB_init_order),LEN(GRIB_SAVE_NAME),LEN(db_ana_tab), &
             LEN(DB_bd_tab), 5)

ALLOCATE (buf(max_data), STAT=ierr)

#ifdef __MPICH2
  IF( me_compute == 0 ) THEN
    bsend_size = max_data*5*num_io
    ALLOCATE (bsend_buf(bsend_size))
    CALL MPI_Buffer_attach(bsend_buf, bsend_size, ierr)
  ENDIF
#endif

IF( ierr /= 0 ) THEN
  WRITE(0,*) me_compute,' : SERIOUS MEMORY LACK IN mpe_io_init '
  RETURN
ENDIF

IF( is_compute_pe ) THEN
  ! First determine root PE which has read namelist

  IF( DB_read_done ) THEN
    iroot = me_compute
  ELSE
    iroot = -1
  ENDIF

  CALL MPI_ALLREDUCE(iroot, root, 1, MPI_INTEGER, MPI_MAX, mpe_comm_compute, ierr)

  IF( root < 0 ) THEN

    IF( me_compute == 0 ) THEN
      WRITE(0,*) 'Warning mpe_io_init: readnldb NOT called !! '
      WRITE(0,*) 'DATABASE disabled '
    ENDIF
    root = 0
    DB_OUT_SOCKETS = 0
    DB_IN_SOCKETS  = 0
    DB_retry       = 0
    GRIB_SAVE_SIZE = -1
    MPE_dbg_level  = 0

  ENDIF

  buf(1) = DB_OUT_SOCKETS
  buf(2) = DB_IN_SOCKETS
  buf(3) = DB_retry
  buf(4) = GRIB_SAVE_SIZE
  buf(5) = MPE_dbg_level
  ndata  = 5

  CALL MPI_BCAST (buf,ndata,MPI_INTEGER,root,mpe_comm_compute, ierr)

  DB_OUT_SOCKETS = buf(1)
  DB_IN_SOCKETS  = buf(2)
  DB_retry       = buf(3)
  GRIB_SAVE_SIZE = buf(4)
  MPE_dbg_level  = buf(5)

  IF( me_compute == root ) THEN

    itag = tag_open-1
    IF (io_config == 2) THEN
      ! for io_config==1, the broadcast above did the job already
      DO pe = 0,num_io-1
#ifdef __MPICH2
        CALL MPI_BSEND(buf, ndata, MPI_INTEGER, pe, itag, MPE_COMM_CO2IO, ierr)
#else
        CALL MPI_SEND(buf, ndata, MPI_INTEGER, pe, itag, MPE_COMM_CO2IO, ierr)
#endif
      ENDDO
    ENDIF

    DO pe = 0,num_io-1

      IF ( (io_config == 2) .OR. (pe > 0) ) THEN
        ! io_config==1 and pe==0: I do not need to send to myself

        IF (MAX(DB_OUT_SOCKETS,DB_IN_SOCKETS) > 0 ) THEN
          ndata=max_data
          CALL char_pack(DB_init_order,buf,ndata,1,ierr)
#ifdef __MPICH2
          CALL MPI_BSEND(buf, ndata, MPI_INTEGER, pe, itag, MPE_COMM_CO2IO, ierr)
#else
          CALL MPI_SEND(buf, ndata, MPI_INTEGER, pe, itag, MPE_COMM_CO2IO, ierr)
#endif

          ndata=max_data
          CALL char_pack(DB_ana_tab,buf,ndata,1,ierr)
#ifdef __MPICH2  
          CALL MPI_BSEND(buf, ndata, MPI_INTEGER, pe, itag, MPE_COMM_CO2IO, ierr)
#else
          CALL MPI_SEND(buf, ndata, MPI_INTEGER, pe, itag, MPE_COMM_CO2IO, ierr)
#endif
     
          ndata=max_data
          CALL char_pack(DB_bd_tab,buf,ndata,1,ierr)
#ifdef __MPICH2
          CALL MPI_BSEND(buf, ndata, MPI_INTEGER, pe, itag, MPE_COMM_CO2IO, ierr)
#else
          CALL MPI_SEND(buf, ndata, MPI_INTEGER, pe, itag, MPE_COMM_CO2IO, ierr)
#endif
        ENDIF

        IF( GRIB_SAVE_SIZE >= 0 ) THEN
          CALL char_pack(GRIB_SAVE_NAME,buf,ndata,1,ierr)
#ifdef __MPICH2
          CALL MPI_BSEND(buf, ndata, MPI_INTEGER, pe, itag, MPE_COMM_CO2IO, ierr)
#else
          CALL MPI_SEND(buf, ndata, MPI_INTEGER, pe, itag, MPE_COMM_CO2IO, ierr)
#endif
        ENDIF
      ENDIF

    ENDDO

  ENDIF

ENDIF   ! is_compute_pe

IF( is_io_pe ) THEN
  itag = tag_open-1

  IF (io_config == 2) THEN
    ! for io_config==1, the broadcast above did the job already
    CALL MPI_RECV(buf,max_data,MPI_INTEGER, MPI_ANY_SOURCE, &
                  itag, MPE_COMM_CO2IO, status, ierr)

    DB_OUT_SOCKETS = buf(1)
    DB_IN_SOCKETS = buf(2)
    DB_retry = buf(3)
    GRIB_SAVE_SIZE = buf(4)
    MPE_dbg_level = buf(5)

    !write(0,*) me_io,' Got OUT_SOCKETs   =',DB_OUT_SOCKETS
    !write(0,*) me_io,' Got  IN_SOCKETs   =',DB_IN_SOCKETS
    !write(0,*) me_io,' Got DB_retry      =',DB_retry
    !write(0,*) me_io,' Got DB_bcksize    =',GRIB_SAVE_SIZE
    !write(0,*) me_io,' Got MPE_dbg_level =',MPE_dbg_level
  ENDIF


  IF (io_config == 2) THEN
    IF( max(DB_OUT_SOCKETS,DB_IN_SOCKETS) > 0 ) THEN
      CALL MPI_RECV(buf,max_data,MPI_INTEGER, MPI_ANY_SOURCE, &
                    itag, MPE_COMM_CO2IO, status, ierr)
      CALL MPI_GET_COUNT(status, MPI_INTEGER, ndata, ierr)
      CALL char_pack(DB_init_order,buf,ndata,-1,ierr)
      !write(0,*) ' Received order: ',DB_init_order

      CALL MPI_RECV(buf,max_data,MPI_INTEGER, MPI_ANY_SOURCE, &
                    itag, MPE_COMM_CO2IO, status, ierr)
      CALL MPI_GET_COUNT(status, MPI_INTEGER, ndata, ierr)
      CALL char_pack(DB_ana_tab,buf,ndata,-1,ierr)

      CALL MPI_RECV(buf,max_data,MPI_INTEGER, MPI_ANY_SOURCE, &
                    itag, MPE_COMM_CO2IO, status, ierr)
      CALL MPI_GET_COUNT(status, MPI_INTEGER, ndata, ierr)
      CALL char_pack(DB_bd_tab,buf,ndata,-1,ierr)
      !write(0,*) ' Received tables ',DB_ana_tab, ' ', DB_bd_tab

      IF( GRIB_SAVE_SIZE >= 0 ) THEN
        CALL MPI_RECV(buf,max_data,MPI_INTEGER, MPI_ANY_SOURCE, &
                      itag, MPE_COMM_CO2IO, status, ierr)
        CALL MPI_GET_COUNT(status, MPI_INTEGER, ndata, ierr)
        CALL char_pack(GRIB_SAVE_NAME,buf,ndata,-1,ierr)
        !write(0,*) ' Received backup dir: ',GRIB_SAVE_NAME
      ENDIF
    ENDIF
  ENDIF

ENDIF  ! is_io_pe

DEALLOCATE(buf)

#ifdef __MPICH2
IF( me_compute == 0 ) THEN
  CALL MPI_Buffer_detach(bsend_buf, bsend_size, ierr)
  DEALLOCATE(bsend_buf)
ENDIF
#endif

DB_do_backup = (GRIB_SAVE_SIZE >= 0 .and. DB_OUT_SOCKETS > 0 )

IF( DB_IN_SOCKETS > 0 ) THEN
  in_device  = 'bank'
ELSE
  in_device = 'file'
ENDIF

IF( DB_OUT_SOCKETS > 0 ) THEN
  out_device  = 'bank'
ELSE
  out_device = 'file'
ENDIF

END SUBROUTINE mpe_db_init

! -------------------------------------------------------------

SUBROUTINE mpe_db_tab(ytabtyp,ydate,yorder,iolen)

IMPLICIT NONE

CHARACTER *(*), INTENT(IN) :: ytabtyp, ydate
CHARACTER *(*), INTENT(OUT):: yorder
INTEGER, INTENT(IN):: iolen

INTEGER:: curr_len

!print*, ' new db_tab ' 


curr_len = LEN(yorder)
yorder(iolen+1:curr_len) = ' '

IF( ytabtyp == 'ana' ) THEN

IF( iolen > 0 ) THEN
yorder = yorder(1:iolen)//',d='//ydate//','//DB_ana_tab(1:LEN_TRIM(DB_ana_tab))
ELSE
yorder = 'd='//ydate//','//DB_ana_tab(1:LEN_TRIM(DB_ana_tab))
ENDIF


ELSE

IF( iolen > 0 ) THEN
yorder = yorder(1:iolen)//',d='//ydate//','//DB_bd_tab(1:LEN_TRIM(DB_bd_tab))
ELSE
yorder = 'd='//ydate//','//DB_bd_tab(1:LEN_TRIM(DB_bd_tab))
ENDIF

ENDIF

DB_input_order = yorder

END SUBROUTINE mpe_db_tab

SUBROUTINE mpe_get_db_tab(ymode,ytabtyp)

IMPLICIT NONE

CHARACTER *(*), INTENT(IN ):: ymode
CHARACTER *(*), INTENT(OUT):: ytabtyp

IF ( ymode == 'out' ) THEN

  ytabtyp = TRIM(DB_output_order)

ELSE

  ytabtyp = TRIM(DB_input_order)

ENDIF

END SUBROUTINE mpe_get_db_tab



!+ Initialization of the io module (Constructor)
!
SUBROUTINE mpe_io_init( icomm_ori, icompute_pes, io_pes, icomm_compute, ierror )
!
!--------------------------------------------------------------------------
!
! Definitions  : compute PE : PE that does computations and calls
!                             mpe_io_write to store arrays 
!
!                io      PE : PE that receives data from compute 
!                             PEs and dumps it directly to a 
!                             device
!                             PE can be io- and compute PE 
!
!
! Description  : Build the communicator for this 
!                io system and a new communicator for compute-PEs.
!                Initialize the chosen io devices.     
!                          
!                This subroutine has to be called after MPI_INIT! 
!
! Method       : 
!
!                There are 2 IO Configurations:
!
!                1. Every io PE is a compute PE ( traditional method )
!                2. Every io PE is not a compute PE ( separat io PEs )
!
!                The two configurations are determined by the number 
!                of compute PEs ( icompute_pes ), the number 
!                of io PEs ( io_pes ) together with the number 
!                of PEs inside icomm_ori (ori_pes).
!                
!                1: icompute_pes          = ori_pes
!                2: icompute_pes + io_pes = ori_pes
!
!                Currently io_pes <= icompute_pes is necessary.
!
!                Communication is made more safe by using different
!                communicators for different (sub)systems. 
! 
!                Communicator:       Purpose:
!                --------------------------------------------
!
!                MPI_COMM_COMPUTE    communicator for computation PEs
!                MPI_COMM_IO         communicator for io PEs
!                MPI_COMM_CO2IO      communicator for communication
!                                    between io and compute PEs
!
!                In case of separat io and compute group 
!                MPI_COMM_CO2IO is an intercommunicator.
!
! error treatment:
!
!             -  MPI is not initialized 
!             -  Use of a nonsense configuration                            
!             -  Device failed ( no Database response )    
!             -  allocation error   
!
! debug output   : 
!             
!--------------------------------------------------------------------------

  IMPLICIT NONE

!--------------------------------------------------------------------------
!
! Parameters       :

INTEGER, INTENT(IN)  ::  icomm_ori           &! original communicator
                       , icompute_pes        &! number of PEs for computation
                       , io_pes              ! number of PEs for io
             
INTEGER, INTENT(OUT) ::  icomm_compute        ! new communicator for 
                                             ! compute PEs

INTEGER, INTENT(OUT) ::  ierror               ! error return value
!
!--------------------------------------------------------------------------
!
! Local parameters :

INTEGER           :: num_ori     &! number of PEs in ori. communicator
                   , me_ori      &! rank in original comm 
                   , key         &! key for MPI_COMM_SPLIT
                   , color_io    &! color   MPI_COMM_SPLIT (I/O)
                   , color_comp  &! color   MPI_COMM_SPLIT (compute)
                   , ierr        &! error parameter
                   , r_leader    &! remote leeder PE for intra comm..
                   , comm_local   ! local communicator

!---------------- DB parameter -------------------------------------------

INTEGER          :: dummy_len, dummy_data

!-------------------------------------------------------------------------

! error handling status

  error_handling = .FALSE.
  ori_out_device = ''

! define error variable 

  ierror = 0

! default for output requests

  BUFF_max_requests = 1

   num_io      = io_pes           ! copy dummy parameter to module variable
   num_compute = icompute_pes

   ! determine number of PEs

   CALL MPI_COMM_SIZE( icomm_ori, num_ori, ierr)
   CALL MPI_COMM_RANK( icomm_ori, me_ori , ierr)

   ! correct nonsense configurations

   IF(num_io > num_ori) num_io = num_ori

   IF( num_ori == 1     ) THEN   ! no sep io PEs

      num_compute = num_ori
      num_io      = 1

   ENDIF


!--------------------------------------------------------------------------

   ! Initialisation for the first io-configuration with 
   ! "included" io PEs

   IF( num_compute == num_ori ) THEN    

      io_config = 1

      IF( me_ori < num_io ) THEN 
         is_io_pe    = .TRUE.
         color_io    =  0
      ELSE
         is_io_pe    = .FALSE.
         color_io    =  MPI_UNDEFINED 
      ENDIF

      ! define I/O communicator
      
      key = 0
      CALL MPI_COMM_SPLIT( icomm_ori, color_io, key, MPE_COMM_IO, ierr)

      IF ( is_io_pe ) THEN
         CALL MPI_COMM_RANK( MPE_COMM_IO, me_io, ierr)
      ELSE
         me_io = -1
      ENDIF

      ! compute communicator is unchanged

      icomm_compute = icomm_ori
      is_compute_pe = .TRUE.

      ! "inter" communicator MPE_COMM_CO2IO is just icomm_compute

      CALL MPI_COMM_DUP( icomm_compute, MPE_COMM_CO2IO, ierr)

!-----------------------------------------------------------------------

   ELSE IF ( (num_compute + num_io) == num_ori ) THEN 
      
      ! configuration with separat io PEs

      io_config = 2

      r_leader  = num_ori - num_io ! remote leader PE

      IF( me_ori >= r_leader ) THEN     ! last PEs = I/O PE
         is_io_pe      = .TRUE.
         color_io      =  0
         is_compute_pe = .FALSE.
         color_comp    =  MPI_UNDEFINED
      ELSE
         is_io_pe      = .FALSE.
         color_io      =  MPI_UNDEFINED
         is_compute_pe = .TRUE.
         color_comp    =  0 
      ENDIF

      ! define I/O communicator

      key = 0
      CALL MPI_COMM_SPLIT( icomm_ori, color_io, key, MPE_COMM_IO, ierr)

      IF ( is_io_pe ) THEN
         CALL MPI_COMM_RANK( MPE_COMM_IO, me_io, ierr)
         me_compute = -1
      ELSE
         me_io = -1
      ENDIF

      ! PE 0,..., num_pe-io_pes-1 are compute PEs

      CALL MPI_COMM_SPLIT( icomm_ori, color_comp, key, icomm_compute, ierr )

      ! io group and comm group are disjunct the will communicate
      ! via the true inter communicator MPE_COMM_CO2IO
      ! first we have to build the inter communicator

      if( is_compute_pe ) THEN

         comm_local = icomm_compute

      ELSE

         comm_local = MPE_COMM_IO
         r_leader   = 0

      ENDIF

      CALL MPI_INTERCOMM_CREATE( comm_local          &! local communicator
                                  , 0                &! local leader PE
                                  , icomm_ori         &! "peer" communicator
                                  , r_leader         &! remote leader 
                                  , inter_tag        &! safe tag
                                  , MPE_COMM_CO2IO   &! inter communicator
                                  , ierr )



   ELSE !------------------------------------------------------------------

      print*, " non implemented config: !"
      print*, " num_compute =            ", icompute_pes
      print*, " num_io      =            ", num_io
      print*, " num_ori     =            ", num_ori

   ENDIF !-----------------------------------------------------------------

   IF( is_compute_pe )  THEN

      ! duplicate this communicator for internal use

      CALL MPI_COMM_DUP( icomm_compute, MPE_COMM_COMPUTE, ierr)

      CALL MPI_COMM_SIZE( MPE_COMM_COMPUTE , num_compute, ierr)
      CALL MPI_COMM_RANK( MPE_COMM_COMPUTE , me_compute , ierr)
      
   ENDIF

!  Initialize the device ( bank,...)

   CALL DEV_Init( ierror )

   CALL DB_err_handler( ierror, .TRUE., DB_err)

!  allocate space for the Buffer on io_pes

   IF( is_io_pe ) THEN

      allocate( BUFF_space(BUFF_length, BUFF_max_requests), STAT= ierror)

      allocate( BUFF_len_curr( BUFF_max_requests))

      IF( ierror /= 0 ) write(0,*) 'allocation error in init'

   ENDIF

!  allocate space for send buffer

   IF (is_compute_pe) THEN

      allocate( send_buffer(BUFF_length, send_max_request), STAT = ierror)

      IF( ierror /= 0 ) write(0,*) 'allocation error in init'

   ENDIF

END SUBROUTINE mpe_io_init

!--------------------------------------------------------------------------

SUBROUTINE DEV_Init( error )

IMPLICIT NONE

INTEGER, INTENT(OUT) :: error

LOGICAL::first=.TRUE.

INTEGER              :: back_err, cso_dbg_level

   error = 0

!------------- DATABASE code ----------------------------------------------

   IF( first) THEN
     CALL mpe_db_init()
     first=.FALSE.
   ENDIF

   IF ( (out_device == 'bank' .OR. in_device == 'bank') .AND. is_io_pe ) THEN

      ! initialize bank

      IF( MPE_dbg_level > 1 ) THEN
      print*, 'USING CSOBANK with ', DB_OUT_SOCKETS, ' OUT_SOCKETS '&
                                   , DB_IN_SOCKETS , ' IN_SOCKETS' 
      ENDIF

      IF( MPE_dbg_level > 0 ) THEN
      cso_dbg_level = MPE_dbg_level + me_io
      ELSE
      cso_dbg_level = 0
      ENDIF
      CALL csodb_init( DB_init_order, DB_OUT_SOCKETS, DB_IN_SOCKETS, &
                       DB_retry, cso_dbg_level, error)


      db_init_count = db_init_count + 1


      IF( db_init_count == db_init_err ) error = 6

      IF( error /= 0 ) THEN

        write(0,*) ' open bank with DB_err = ', error,' on ', me_io

      ENDIF

      IF( out_device == 'bank' .AND. error == 0) THEN

        ! initialize backup system

        IF( DB_do_backup ) THEN

          CALL GRIB_SAVE_INIT( trim(GRIB_SAVE_NAME) &! Backup directory
                          , me_io               &! Backup index
                          , GRIB_SAVE_SIZE      &! Backup incore size
                          , back_err   )         ! Backup system error

          IF( back_err /= 0 ) THEN
            write(0,*) 'error in GRIB_SAVE_INIT on ', me_io
          ENDIF

        ENDIF
    
        BUFF_max_requests = DB_OUT_SOCKETS

      ENDIF 

   ENDIF

END SUBROUTINE DEV_Init

!--------------------------------------------------------------------------

SUBROUTINE mpe_io_reconfig( icomm_compute, ierr )

!--------------------------------------------------------------------------
!
! $Revision: 1.1.1.6 $
! $DATE$
!
! Purpose          :  reconfigures this io module by using 
!                     a different communicator for the compute 
!                     PEs. IO PEs will stay the same.
!
!                     It is collective and blocking for  
!                     current icompute_pes
!
!--------------------------------------------------------------------------

  IMPLICIT NONE

!--------------------------------------------------------------------------
!
! Parameters       :
!
                    
  INTEGER, INTENT(IN)  :: icomm_compute   ! compute PE communicator
  INTEGER, INTENT(OUT) :: ierr            ! error status
                                        

!--------------------------------------------------------------------------
!
! Local parameters :
!

  INTEGER          :: result      ! result of comm_compare  

!--------------------------------------------------------------------------
!
! Error checks     :
!
!--------------------------------------------------------------------------

  ierr = 0

  CALL MPI_COMM_COMPARE( icomm_compute, MPE_COMM_COMPUTE, result, ierr)
  IF( result /= MPI_CONGRUENT ) THEN 
     ierr = 1 
     write(0,*)  'needs reconfiguration'
  ENDIF 

 END SUBROUTINE mpe_io_reconfig

 !---------------------------------------------------------------------

 LOGICAL FUNCTION mpe_is_compute()

 mpe_is_compute = is_compute_pe

 END FUNCTION

 !---------------------------------------------------------------------

 LOGICAL FUNCTION mpe_is_io()

 mpe_is_io = is_io_pe

 END FUNCTION

 !---------------------------------------------------------------------

 SUBROUTINE mpe_io_open( nudat, datname, db_in_order, mode, lstop_io, ierror )

 !--------------------------------------------------------------------------
 !
 ! $Revision: 1.1.1.6 $
 ! $DATE$
 !
 ! Purpose          : opens io file only on I/O PEs
 !
 ! Description      : filename has to be provided on the 
 !                    first num_io compute nodes. This information
 !                    is broadcasted to the io nodes. 
 !
 !                    the current method will not work if 
 !                    num_io > num_compute
 !--------------------------------------------------------------------------

   IMPLICIT NONE

 !--------------------------------------------------------------------------
 !
 ! Parameters       :
 !

   INTEGER             , INTENT(OUT) :: nudat     ! file descriptor

   CHARACTER (LEN = *) , INTENT(IN)  :: datname   ! name of datafile
   CHARACTER (LEN = *) , INTENT(IN)  :: db_in_order &
                                                  ! database input order
                                                  ! in case mode = 'r'
                                      , mode      ! file mode

   LOGICAL             , INTENT(OUT) :: lstop_io   ! tells pure io 
                                                  ! nodes to stop

   INTEGER, INTENT(OUT), OPTIONAL    :: ierror
 
 !--------------------------------------------------------------------------
 !
 ! Local parameters :
 !

  INTEGER                           :: i,i0        ! loop_counter

  INTEGER (KIND=gribc_kind)         :: nudatc     &! local file descriptor
                                     , ierrc       ! error for file open

  INTEGER, parameter      :: buff_len = 256

  INTEGER                 :: dat_info(0:buff_len)  ! integer array for
                                                   ! datname and mode

  INTEGER                 :: dat_len, mode_len    &! length for datname, moder
                           , dbord_len, info_len   ! and dborder

  INTEGER                 :: tag                   ! tag for sending dat_info

  CHARACTER  (LEN = 250)  :: datname_new           ! new name and mode for 
  CHARACTER (LEN = 3)     :: mode_new              ! datname (parallel io) 

  INTEGER                 :: ierr                 &! error parameter
                        , status(MPI_STATUS_SIZE)  

  INTEGER                 :: dummy, dummy_len = 0

  INTEGER                 :: error

 !--------------------------------------------------------------------------
 !
 ! Error checks     :
 !
 !--------------------------------------------------------------------------

 !  write(0,*) 'step into mpe_io_open, me = ', me_compute, me_io

  lstop_io = .FALSE.
  got_end_signal = .FALSE.
  mode_new(1:3)  = '   '

  ! datname and mode is copied to temporary variables on compute
  ! nodes. 
   
  IF (PRESENT(ierror)) THEN
    ierror = 0
  ENDIF

  dat_len   = LEN_TRIM(datname)
  mode_len  = LEN_TRIM(mode)
  dbord_len = LEN_TRIM(db_in_order)
  DB_input_order = db_in_order(1:dbord_len)

  IF ( is_compute_pe .AND. is_io_pe ) THEN

     datname_new(1:dat_len) = datname(1:dat_len)

     mode_new(1:3)  = '   '
     mode_new(1:mode_len)   = mode(1:mode_len)

  ENDIF

  IF ( io_config == 2 ) THEN

     ! datname and mode have to communicated to io nodes
     ! code datname and mode on dat_info

     end_read = .FALSE.

     IF ( is_compute_pe ) THEN


        IF( mode(1:1) == 'r' .and. in_device == 'bank' ) THEN

            dat_len = dbord_len
            dat_info(0) = dbord_len
            info_len = dbord_len

            CALL char_pack(DB_input_order(1:dbord_len),dat_info(1),info_len, &
                           1,error)

            !WRITE(0,*) me_compute,' Send db_order ',DB_input_order(1:dbord_len)

          ELSE

            dat_info(0) = dat_len  
            info_len = dat_len
            datname_new(1:dat_len) = datname(1:dat_len)

            CALL char_pack(datname_new(1:dat_len),dat_info(1),info_len, &
                           1,error)

          ENDIF

          dat_info(dat_len+1) = mode_len

          i0 = dat_len + 2
          DO i = i0, dat_len + mode_len + 1
               dat_info(i) = ichar(mode(i-i0+1:i-i0+1))
               !print*, 'dat_inf, mode =', dat_info(i), mode(i-i0+1:i-i0+1)
        ENDDO

        info_len = dat_len + mode_len + 2

        ! send dat_info to io PEs

        DO i = 0, num_io - 1
           IF( me_compute == i ) THEN

             CALL MPI_SEND( dat_info, info_len, MPI_INTEGER, i, tag_open, &
                            MPE_COMM_CO2IO, ierr ) 
           ENDIF
        ENDDO

     ELSE


       DO

          CALL MPI_RECV( dat_info, buff_len, MPI_INTEGER, me_io, tag_open, &
                             MPE_COMM_CO2IO, status, ierr ) 


          dat_len = dat_info(0)


          IF( MPE_dbg_level > 1 ) THEN
          write(0,*) 'received open message with dat_len = ', dat_len
          ENDIF

          IF ( dat_len == -1) THEN

             lstop_io = .TRUE.
             return

          ! ready message: compute PEs have finished current timestep
          ! and send a message to the waiting IO-PEs

          ELSE IF ( dat_len == -2 ) THEN 

             CALL MPI_SEND( dummy, dummy_len, MPI_INTEGER, me_io, ready_tag,&
                           MPE_COMM_CO2IO, ierr)


          ELSE IF ( dat_len == -3 ) THEN

             CALL DEV_ready( dat_info )

          ELSE

            EXIT

          ENDIF

        ENDDO

        ! decode datname (or input_order) and mode


        mode_len = dat_info(dat_len+1)
        mode_new(1:3)  = '   '     
        i0 = dat_len + 2
        DO i = dat_len+2, dat_len + mode_len + 1
           mode_new(i-i0+1:i-i0+1) = achar( dat_info(i))
        ENDDO

        ! WRITE(0,*) ' mode = ',mode_new,' in_device = ',in_device

        IF(  mode_new(1:1) == 'r' .and. in_device == 'bank' ) THEN

          do i = 1, dat_len
           DB_input_order(i:i) = achar( dat_info(i))
          enddo

          ! WRITE(0,*) ' unpacked order: ',DB_input_order

        ELSE

          do i = 1, dat_len
           datname_new(i:i) = achar( dat_info(i))
          enddo

        ENDIF

      ENDIF

  ENDIF


  IF( mode_new(1:1) == 'r' ) THEN
     current_mode = 'r'
  ELSE
     current_mode = 'w'
  ENDIF

  IF ( is_io_pe ) THEN

 !write(0,*) 'out_device = ', out_device, ' mode_new = ', mode_new(1:1), num_io

     CALL DEV_Open( nudat, datname_new, dat_len,  mode_new, error )

     IF( in_device == 'bank' .AND. mode_new(1:1) == 'r' ) THEN


     DB_input_order = 'ak=re,arr=j,'//DB_input_order(1:LEN_TRIM(DB_input_order))
     DB_order = TRIM(DB_input_order)
     
     DB_hits = 2

     IF( MPE_dbg_level > 1 ) THEN
     write(0,*) 'DB_input_order = ', DB_order(1:LEN_TRIM(DB_order))
     ENDIF

     ENDIF

     field_rcv = me_io - num_io  

  ENDIF

  field_snd = me_compute - num_compute 

  !write(0,*) 'step off mpe_io_open, me = ', me_compute, me_io, ' ',&
  !           field_rcv, field_snd

 END SUBROUTINE mpe_io_open

!--------------------------------------------------------------------------

SUBROUTINE DEV_Open( nudat, datname, dat_len,  mode, ierror )

IMPLICIT NONE


INTEGER           , INTENT(OUT)   :: nudat
CHARACTER  (LEN=*), INTENT(IN)    :: datname
INTEGER           , INTENT(IN)    :: dat_len
CHARACTER  (LEN=3), INTENT(IN)    :: mode
INTEGER           , INTENT(OUT)   :: ierror

!--------------------------------------------------------------------------

INTEGER                             :: new_len
CHARACTER    (LEN= dat_len+3)       :: datname_new
INTEGER      (KIND=gribc_kind)      :: nudatc     &! local file descriptor
                                      ,ierrc       ! local file error
new_len     = dat_len
datname_new = datname(1:dat_len)

!nudatc = nudat
ierror = 0

IF( out_device == 'file' .AND. current_mode == 'w' ) THEN

   IF( num_io > 1 ) THEN

      !write(0,*) 'num_io = ', num_io

!<Change 1.2
!     write( datname_new(dat_len+1:dat_len+3),'(A1,i1)') '_',me_io
      if( me_io < 10 ) then
      write( datname_new(dat_len+1:dat_len+2),'(A1,I1)') '_',me_io
      else if( me_io < 100 ) then
      write( datname_new(dat_len+1:dat_len+3),'(A1,I2)') '_',me_io
      else if( me_io < 1000 ) then
      write( datname_new(dat_len+1:dat_len+4),'(A1,I3)') '_',me_io
      endif
!>

      new_len = dat_len+3

   ENDIF

   !write(0,*) 'open ', datname_new(1:new_len), ' and mode', mode(1:1)

#ifdef GRIBDWD
   CALL copen(nudatc,datname_new(1:new_len), mode(1:3), ierrc) 
   ierror = ierrc
#endif

ENDIF

IF( in_device == 'file' .AND. current_mode == 'r' .AND. me_io == 0) THEN

#ifdef GRIBDWD
   CALL copen(nudatc,datname_new(1:new_len), mode(1:3), ierrc)
   ierror = ierrc
#endif

ENDIF

IF( out_device == 'bank' .and. DB_do_backup .AND. current_mode == 'w') THEN

   !write(0,*) 'before gribsave open: ', datname_new(1:new_len)

   CALL GRIB_SAVE_OPEN( datname_new(1:new_len), me_io, ierror) 

   IF( ierror /= 0 ) THEN
     write(0,*) 'error in GRIB_SAVE_OPEN : ', ierror 
   ENDIF
   
   
ENDIF

IF( ierror == 0 ) THEN
  open_count = open_count + 1
ENDIF


current_datname = ' '
IF( current_mode == 'w' ) current_datname = datname(1:dat_len)
mode_backup     = mode(1:3) 
nudat    = nudatc


END SUBROUTINE DEV_Open

!--------------------------------------------------------------------------

!+ checks whether a postprocessing step has finished

SUBROUTINE mpe_io_ready(file_name)

IMPLICIT NONE

CHARACTER (LEN=*) :: file_name      ! filename of the ready file



INTEGER           :: file_len
INTEGER           :: dummy(0:max_name_len), len, ierr,i


!write(0,*) 'jump into  ready file', me_io, me_compute

                        
! Encode filename on integer, dummy(0) = -3 is identifier in open

IF( is_compute_pe ) THEN

file_len = len_trim( file_name )

dummy(0) = -3
len      =  max_name_len + 1

DO i = 1, file_len
   dummy(i) = ichar(file_name(i:i))
ENDDO  

DO i = file_len+1, max_name_len
   dummy(i) = ichar(' ')
ENDDO


IF( is_io_pe ) THEN

  CALL DEV_ready( dummy )

ELSE

  !send ready message to io nodes 

  IF( me_compute < num_io ) THEN

     CALL MPI_SEND( dummy, len, MPI_INTEGER, me_compute, tag_open, &
          MPE_COMM_CO2IO, ierr )

  ENDIF

ENDIF

ENDIF

END SUBROUTINE mpe_io_ready

!--------------------------------------------------------------------------

SUBROUTINE DEV_ready( dummy ) 

IMPLICIT NONE


INTEGER                      :: dummy(0:*)
INTEGER, parameter           :: max_files = 100
INTEGER                      :: i, error,  err_close, err_init, err_in, unit
INTEGER                      :: status, nudat, err
CHARACTER (LEN=max_name_len) :: file_name
LOGICAL                      :: dummy_pend(BUFF_max_requests)

!< Change 1.2
INTEGER                      :: ihits
!>

IF( is_io_pe ) THEN

CALL MPI_BARRIER( MPE_COMM_IO, error)

IF( out_device == 'bank' .OR. ori_out_device == 'bank' ) THEN

! close database


  !write(0,*) 'before csodb_close', me_io
!< Change 1.2
! CALL csodb_close( err_close )
  CALL csodb_arr_wait_out(ihits, err_close)
!>

  ! write(0,*) 'after csodb_close me_io , err', me_io, err_close
  ! in case of error write down buffers in files

  IF ( err_close /= 0 .AND. .NOT.error_handling ) THEN

      CALL DB_err_handler_write( err_close, nudat, BUFF_length*gribf_kind&
            ,dummy_pend, error ) 

      IF( open_count > 0 ) CALL mpe_io_close( nudat, error )

  ENDIF

  ! delete GRIB_SAVE buffer

   IF( DB_do_backup ) THEN

     CALL GRIB_SAVE_DELETE( me_io, error)

     IF( error /= 0 ) THEN
       write(0,*) 'error in GRIB_SAVE_DELETE', me_io
     ENDIF
    

     CALL GRIB_SAVE_INIT( trim(GRIB_SAVE_NAME) &! Backup directory
                         , me_io               &! Backup index
                         , GRIB_SAVE_SIZE      &! Backup incore size
                         , error   )            ! Backup system error

     IF( error /= 0 ) THEN
       write(0,*) 'error in GRIB_SAVE_DELETE', me_io
     ENDIF


   ENDIF

ENDIF

! communicate error handling

err = 0
IF( error_handling ) err = 100

CALL MPI_ALLREDUCE( err, status, 1, MPI_INTEGER, MPI_MAX&
     , MPE_COMM_IO, error) 

  IF( status /= 0 ) error_handling = .TRUE.

  IF( me_io == 0 ) THEN

    ! write ready file on first io pe

    DO i = 1, max_name_len
       file_name(i:i)  = achar( dummy(i) )
    ENDDO 

    unit = 13
    open (unit, file=trim(file_name), form='formatted')

    IF (error_handling ) THEN
      write(unit, '(a)') 'database failed'
    ELSE
      write(unit, '(a)') 'ready'
    ENDIF

    close(unit)

  ENDIF

  IF( out_device == 'bank' .OR. ori_out_device == 'bank' ) THEN

    ! Try to reopen DATABASE

    out_device = 'bank'

    CALL DEV_Init( err_init )


    IF ( err_init == 0 .AND. error_handling ) THEN

      ! change from error handling to normal mode

      write(0,*) ' END ERROR HANDLING, go back to out_device ', out_device&
                ,' on PE ', me_io 

      ori_out_device = ''
      error_handling = .FALSE.

    ENDIF

    if( error_handling ) out_device = 'file'

    IF ( err_init /= 0 .AND. .NOT. error_handling ) THEN

      ! error handling for DEV_Init

      CALL DB_err_handler( err_init, .TRUE., error)

    ENDIF

  ENDIF

ENDIF



END SUBROUTINE DEV_ready

!---------------------------------------------------------------------------

!+ waits for completion of an postprocessing step

SUBROUTINE mpe_io_ready_wait( nstep )

IMPLICIT none

INTEGER, INTENT(IN) :: nstep


LOGICAL           :: lreceived
INTEGER           :: status(MPI_STATUS_SIZE)
INTEGER           :: dummy, len, ierr, ierr1
REAL              :: time_current = 0.0, time_init
REAL              :: time_delta = 0.0, time_stamp, tprev, tintrpt
REAL, PARAMETER   :: x_delta = 10.0 ! Maximum assumed idle time
DOUBLE PRECISION MPI_WTIME          ! in waiting for IO PEs

len = 1
dummy = -2


IF( io_config == 2 ) THEN

  IF( is_compute_pe ) THEN

     IF ( nstep > 0 .AND. me_compute < num_io ) THEN

         CALL MPI_SEND( dummy, len, MPI_INTEGER, me_compute, tag_open, &
                     MPE_COMM_CO2IO, ierr )

         time_init = MPI_WTIME()
         time_current = 0.0
         tprev = 0.0
         tintrpt = 0.0

         lreceived = .FALSE.

         write(0,'(A,I4,A,I4,A)')  &
                 ' Idle waiting for IO-PE on compute PE ',me_compute, &
                 '; maximum time will be ',DB_retry,'s'

         DO WHILE( time_current < DB_retry )

           CALL MPI_IPROBE( me_compute, ready_tag, &
                       MPE_COMM_CO2IO,lreceived, status, ierr) 
        
           IF( lreceived ) EXIT

           time_stamp = MPI_WTIME()-time_init
           time_delta = time_stamp-tprev
           tprev = time_stamp
  
           IF( time_delta > x_delta ) THEN
! time_delta >= x_delta: interpreted as  irregular interrupt
! (suspension) and thus ignored
             tintrpt = tintrpt+time_delta
           ENDIF

           time_current = time_stamp-tintrpt

         ENDDO

         IF( lreceived ) THEN

           write(0,'(A,f4.1,A,I4,A)')  &
                 ' ... took ',time_current,'s on ',me_compute

         
           CALL MPI_RECV( dummy,len, MPI_INTEGER, me_compute, ready_tag, &
                      MPE_COMM_CO2IO, status, ierr)

         ELSE

           write(0,*) me_compute,' WAITED ', time_current, &
                      ' seconds for IO PE to respond; now ABORTING!! '      

           ierr1 = 1
           CALL MPI_ABORT( MPI_COMM_WORLD, ierr1, ierr)

         ENDIF


       ENDIF ! nstep

       CALL MPI_BARRIER( MPE_COMM_COMPUTE, ierr)

    ENDIF    ! is compute_pe

ENDIF        ! io_config == 2

END SUBROUTINE mpe_io_ready_wait

!--------------------------------------------------------------------------

 SUBROUTINE mpe_io_next_snd( pe_target, tag_target, field_idx)

!--------------------------------------------------------------------------
!
! $Revision: 1.1.1.6 $
! $DATE$
!
! Purpose          : determine next target PE for I/O
!
! Modules          :
! 
! Module parameters:
!
!--------------------------------------------------------------------------

  IMPLICIT NONE

!--------------------------------------------------------------------------
!
! Parameters       :
!

INTEGER, INTENT(OUT) :: pe_target
INTEGER, INTENT(OUT) :: tag_target
INTEGER, INTENT(IN)  :: field_idx

!
!--------------------------------------------------------------------------
!
! Local parameters :
!
!--------------------------------------------------------------------------
!
! Error checks     :
!
!--------------------------------------------------------------------------


IF( is_compute_pe .and. mod(field_idx, num_compute) == me_compute) THEN
   field_snd = field_snd + num_compute
   pe_target = mod(field_snd, num_io)

   !write(0,*) 'Send data #',field_snd,' from ', me_compute, ' to ', pe_target


!   if( io_config == 2 ) pe_target = pe_target + num_compute 

ELSE
   pe_target = MPI_PROC_NULL
ENDIF

tag_target = field_snd


END SUBROUTINE mpe_io_next_snd


!--------------------------------------------------------------------------

 SUBROUTINE mpe_io_next_recv( pe_source, tag_source, field_idx)

!--------------------------------------------------------------------------
!
! $Revision: 1.1.1.6 $
! $DATE$
!
! Purpose          : determine next source PE for I/O 
!
! Modules          :
! 
! Module parameters:
!
!--------------------------------------------------------------------------

  IMPLICIT NONE

!--------------------------------------------------------------------------
!
! Parameters       :
!

INTEGER, INTENT(OUT) :: pe_source
INTEGER, INTENT(OUT) :: tag_source
INTEGER, INTENT(IN)  :: field_idx

!
!--------------------------------------------------------------------------
!
! Local parameters :
!
!--------------------------------------------------------------------------
!
! Error checks     :
!
!--------------------------------------------------------------------------


IF( is_io_pe .AND. MOD(field_idx, num_io) == me_io ) THEN
   field_rcv = field_rcv + num_io
   pe_source = mod(field_rcv, num_compute)
ELSE
   pe_source = MPI_PROC_NULL
ENDIF

tag_source = field_rcv

END SUBROUTINE mpe_io_next_recv

!--------------------------------------------------------------------------

SUBROUTINE mpe_io_send( data, len, ierr)

!--------------------------------------------------------------------------
!
! $Revision: 1.1.1.6 $
! $DATE$
!
! Purpose          :
!
! Modules          :
! 
! Module parameters:
!
!--------------------------------------------------------------------------

  IMPLICIT NONE

!--------------------------------------------------------------------------
!
! Parameters       :
!

INTEGER  (KIND = gribf_kind)         &
        , INTENT(IN)      :: data(*)   ! array to be written  

INTEGER, INTENT(IN)  :: len
INTEGER, INTENT(OUT) :: ierr


!--------------------------------------------------------------------------
!
! Local parameters :
!
INTEGER :: isend_status(MPI_STATUS_SIZE)

INTEGER :: pe_target, tag_target, izlen

INTEGER :: current_buffer_size

!--------------------------------------------------------------------------
!
! Error checks     :
!
!--------------------------------------------------------------------------


IS_COMP: IF( is_compute_pe ) THEN

     CALL mpe_io_next_snd( pe_target, tag_target, me_compute)

     IF ( pe_target /= me_io ) THEN

             !write(0,*) 'send from ', me_compute, ' field #', tag_target &
             !          ,' to PE#', pe_target,' bytes = ', len

        ! On little endian machines, the usage of MPI_BYTE is problematic,
        ! because in the last word of data (which is actually an integer array),
        ! the last bytes could get lost in the byte-stream, if not the full
        ! integer word is filled up:
        ! the easiest solution is: just send 4 additional bytes
        izlen = len + 4

        CALL MPI_Wait( send_requests(send_current_request), isend_status, ierr)

        current_buffer_size = izlen/gribf_kind + 1

        send_buffer(1:current_buffer_size,send_current_request) = data(1:current_buffer_size)

        CALL MPI_Isend( send_buffer(:,send_current_request), izlen, MPI_BYTE, pe_target, tag_target &
             , MPE_COMM_CO2IO, send_requests(send_current_request), ierr)


        send_current_request = mod( send_current_request, send_max_request) + 1

     ENDIF


ENDIF IS_COMP



END SUBROUTINE mpe_io_send


!--------------------------------------------------------------------------

SUBROUTINE mpe_io_recv( data, recv_len, stop_flag, ierr)

!--------------------------------------------------------------------------
!
! $Revision: 1.1.1.6 $
! $DATE$
!
! Purpose          :
!
! Modules          :
! 
! Module parameters:
!
!--------------------------------------------------------------------------

  IMPLICIT NONE

!--------------------------------------------------------------------------
!
! Parameters       :
!


INTEGER  (KIND = gribf_kind)         &
        , INTENT(INOUT)   :: data(*)   ! array to be written  

INTEGER,  INTENT(INOUT)   :: recv_len
LOGICAl,  INTENT(OUT)     :: stop_flag
INTEGER,  INTENT(OUT)     :: ierr


!--------------------------------------------------------------------------
!
! Local parameters :
!
INTEGER  :: pe_source, tag_source  &! pe and tag for source PE
           ,pe_target, tag_target  &! pe and tag for target PE
           ,dummy_buf

INTEGER  :: local_field, ilr, status(MPI_STATUS_SIZE)

INTEGER  :: min_tag, max_tag

LOGICAL  :: flag

!--------------------------------------------------------------------------
!
! Error checks     :
!
!--------------------------------------------------------------------------

           stop_flag = .FALSE.

           CALL mpe_io_next_recv( pe_source, tag_source, me_io)

           IF( pe_source /= me_compute .OR. pe_source == MPI_PROC_NULL ) THEN

             recv_len = 0                

           ENDIF

           !write(0,*) 'try to receive from Pe: ', pe_source

           IF( pe_source /= MPI_PROC_NULL .AND. pe_source /= me_compute ) THEN

              ! check next incomming message

              flag = .false.

              !HB311001 Do while ( .not. flag )

                 !write(0,*) ' waiting for message from pe:',pe_source

                 !HB311001 CALL MPI_IPROBE( pe_source, MPI_ANY_TAG, MPE_COMM_CO2IO &
                 !HB311001      , flag, status, ierr )

              !HB311001 ENDDO

              CALL MPI_PROBE( pe_source, MPI_ANY_TAG, MPE_COMM_CO2IO &
                      , status, ierr )

              flag = .TRUE.

              IF( status(MPI_TAG) == tag_source ) THEN

                 CALL MPI_GET_COUNT( status, MPI_BYTE, recv_len, ierr) 

                 CALL MPI_RECV( data, recv_len, MPI_BYTE, pe_source,&
                      tag_source &
                      , MPE_COMM_CO2IO, status, ierr)

                 ! Because of the problems on little endian machines (see above)
                 ! 4 additional bytes are send, which are now subtracted again:
                 recv_len = recv_len - 4
                 !write(0,*) 'receive #', tag_source,' on ',&
                 !            me_io, ' from ', pe_source, 'len ',recv_len

                 
              ELSE IF ( status(MPI_TAG) /= tag_write_end ) THEN
                 
                 write(0,*) 'wrong tag ', status(MPI_TAG)

                 write(0,*) 'should recv from ', pe_source, ' field #' &
                        , tag_source &
                        ,' to PE#', me_compute


              ENDIF

              IF ( status(MPI_TAG) == tag_write_end ) THEN

                 !write(0,*) 'end io-step handling'


                 ! end the mpe_io_write call on a io PE.
                 ! determine minimal tag (field index) for the last sent 
                 ! field. Receive "end messages from compute PEs"

                 CALL MPI_ALLREDUCE( tag_source, min_tag, 1, MPI_INTEGER &
                      , MPI_MIN                 &
                      , MPE_COMM_IO, ierr)

                 max_tag = min_tag + num_compute 

                 DO WHILE ( tag_source .lt. max_tag ) 

                    !write(0,*) max_tag, tag_source, pe_source

                    if( io_config == 1 .AND. pe_source == me_compute ) &
                         pe_source = MPI_PROC_NULL


                    CALL MPI_RECV( dummy_buf, 1, MPI_BYTE   &
                         , pe_source, tag_write_end       &
                         , MPE_COMM_CO2IO, status, ierr)

                    !write(0,*) 'receive done'

                    CALL mpe_io_next_recv( pe_source, tag_source, me_io )

                 ENDDO

                 got_end_signal = .TRUE.              

                 stop_flag = .TRUE.

             ENDIF

           ENDIF

           IF ( is_compute_pe .AND. pe_source >= (num_compute-num_io))  &
                                                  stop_flag = .TRUE.


END SUBROUTINE mpe_io_recv

!--------------------------------------------------------------------------

SUBROUTINE mpe_io_write( nudat, data, ilen, ilfd, ydborder, ierror)

!--------------------------------------------------------------------------
!
! $Revision: 1.1.1.6 $
! $DATE$
!
! Purpose          : Send Data from compute PEs to I/O PEs and 
!                    write it on I/O PEs
!
! Modules          :
! 
! Module parameters:
!
!--------------------------------------------------------------------------

  IMPLICIT NONE

!--------------------------------------------------------------------------
!
! Parameters       :
!

INTEGER, INTENT(INOUT)    :: nudat    &! file descriptor
                           , ilen      ! number of Bytes to be written
                           
INTEGER, INTENT(IN)       :: ilfd      ! buffer size in intgribf
                     
INTEGER  (KIND = gribf_kind)         &
        , INTENT(INOUT)   :: data(*)   ! array to be written  

INTEGER, INTENT(OUT), OPTIONAL :: ierror  ! for error handling

!------------------------- DB-Parameters ----------------------------------

CHARACTER (LEN=5), INTENT(IN) :: ydborder ! type of DB

!--------------------------------------------------------------------------
!
! Local parameters :
!

CHARACTER (LEN=5) :: DB_type ! type of DB 
INTEGER           :: DB_out_req, DB_finish_req 


INTEGER  (KIND= gribc_kind ) :: nudatc, ilenc, ierrc       
                           ! corresponding variables for C-routines

INTEGER  :: pe_source, tag_source  &! pe and tag for source PE
           ,pe_target, tag_target   ! pe and tag for target PE

INTEGER  :: ierr, local_field, ilr, ilen_new, ifirst, i, status(MPI_STATUS_SIZE)

INTEGER, save  :: write_count = 0

INTEGER  :: min_tag, max_tag

LOGICAL  :: stop_flag, pending_req( DB_MAX_SOCKETS )

INTEGER  :: data_int_len,error, err_in

!--------------------------------------------------------------------------
!
! Error checks     :
!
!--------------------------------------------------------------------------

IF (PRESENT(ierror)) THEN
  ierror = 0
ENDIF

 got_end_signal = .FALSE.
 DB_type        = ydborder(1:5)
 ilen_new       = ilen

! "Translate" C parameters

  ilenc  = ilen

! SEND data to I/O PE

IF ( (out_device == 'bank' .OR.  ori_out_device == 'bank')&
                           .AND. is_compute_pe .AND. ilen > 0 ) THEN

ilen_new = ilen_new + 5*gribf_kind
ifirst  = ilen/gribf_kind + 1
data( ifirst: ifirst+4 ) = (/(ichar(DB_type(i:i)),i=1,5)/)


ENDIF

CALL mpe_io_send( data, ilen_new, ierr)


! receive buffer

IS_IO : IF( is_io_pe ) THEN


   IF( is_compute_pe ) THEN  ! copy buffer from compute PE

      BUFF_request = mod( write_count, BUFF_max_requests ) + 1
      data_int_len = ilen_new/gribf_kind + 5 
      BUFF_space(1:data_int_len,BUFF_request) = data(1:data_int_len)


   ENDIF

   stop_flag = .FALSE.
   pending_req = .false.

   DO WHILE ( .not.stop_flag ) 

      BUFF_request = mod( write_count, BUFF_max_requests ) + 1

      IF ( out_device == 'bank' ) THEN


         IF ( pending_req( BUFF_request ) ) THEN

           CALL DEV_Wait( DB_request( BUFF_request ), DB_err )

           !write(0,*) 'write to DB, DB_hits = ', DB_hits,' ', write_count,&
           !           ' request = ', BUFF_request

           CALL DB_err_handler_write( DB_err, nudat, ilfd*gribf_kind&
                                 ,pending_req, error )

           if( DB_err /= 0 ) THEN

              write(0,*) 'Database failed on mpe_io = ', me_io &
                        ,' gribs written = ', write_count - 1
              !CALL MPI_Abort( MPI_COMM_WORLD, DB_err, ierr)

           ENDIF
           pending_req( BUFF_request) = .false.
        ENDIF

      ENDIF

      ilr     = ilen_new

      !write(0,*) 'before mpe_io_recv BUFF_request, BUFF_max_req = '&
      !     , BUFF_request, BUFF_max_requests

      call mpe_io_recv( BUFF_space(1,BUFF_request),&
                                              ilr, stop_flag, ierr)

      !write(0,*) ' after mpe_io_recv', BUFF_request

      ! in case of error handling the compute PEs don't know that
      ! we are now in filemode and they don't need the last 5
      ! ints which contain the database type

      IF( error_handling ) ilr = ilr - 5*gribf_kind

      !write(0,*) 'before  write g#',write_count,' requestno, DB_req = ',&
      !           BUFF_request,DB_request(BUFF_request) 
    
      CALL DEV_Write( nudat, BUFF_space(1,BUFF_request)     &
               , BUFF_length*gribf_kind, ilr, DB_request(BUFF_request), err_in)

      if( ilr > 0 ) THEN
         write_count                  = write_count + 1
         pending_req(BUFF_request)    = .true.
         BUFF_len_curr( BUFF_request) = ilr

         IF( ilr < 1000 ) write(0,*) 'zero byte ?', me_io, ilr


         !write(0,*) 'wrote grib #',write_count,' error = ', err_in

         CALL DB_err_handler_write( err_in, nudat, BUFF_length*gribf_kind&
                                 ,pending_req, error ) 

      ENDIF


   ENDDO


   IF( out_device == 'bank' ) THEN

     DO i = 1, BUFF_max_requests

        DB_finish_req = mod( BUFF_request + i - 1, BUFF_max_requests) + 1

        IF( pending_req(DB_finish_req) ) THEN

           CALL DEV_Wait(DB_request(DB_finish_req), DB_err )

           CALL DB_err_handler_write( DB_err, nudat, ilfd*gribf_kind&
                                    ,pending_req, error ) 

           IF( DB_hits == 0 ) THEN

              write(0,*) 'Database failed on mpe_io = ', me_io
              !CALL MPI_Abort( MPI_COMM_WORLD, DB_err, ierr)

           ENDIF

        ENDIF
     ENDDO

   ENDIF

ENDIF IS_IO


END SUBROUTINE mpe_io_write

!-------------------------------------------------------------------

SUBROUTINE DEV_Write( nudat, data, buff_len, ilen, request, error)

IMPLICIT NONE

INTEGER                  ,INTENT(INOUT)  :: nudat

INTEGER (KIND = gribf_kind), INTENT(IN)  :: data(*)

INTEGER                    , INTENT(IN)  :: buff_len

INTEGER                    , INTENT(IN)  :: ilen

INTEGER                    , INTENT(INOUT)  :: request

INTEGER                    , INTENT(OUT) :: error


INTEGER  (KIND= gribc_kind ) :: nudatc, ilenc, ierrc       
                           ! corresponding variables for C-routines

INTEGER                                 :: ifirst,err,i

CHARACTER (LEN=5) :: DB_type ! type of DB


error = 0 
nudatc = nudat

IF( ilen > 0 ) THEN
   ilenc = ilen 

   IF ( out_device == 'file' ) THEN

      ierrc = 0 

#ifdef GRIBDWD
      !write(0,*) 'before cuegex ', ierrc, ilenc

      call cuegex( nudatc, data, ilenc, ierrc)

      !write(0,*) 'after cuegex ', ierrc, ilenc
#endif

      error = ierrc

   ELSE

      ilenc = ilenc - 5*gribf_kind
      
      ifirst = ilenc/gribf_kind + 1
      write(DB_type,'(5A)')&
           char(data(ifirst)),char(data(ifirst+1)),char(data(ifirst+2)),&
           char(data(ifirst+3)),char(data(ifirst+4))
! Note: This cannot be compiled by the GNU fortran compiler, since it has not been implemented
! yet. But code is not active for normal applications.
!      write(DB_type,'(5A)')&
!           (/(char(data(i)),i=ifirst,ifirst+4)/)   

      !WRITE(0,*) ' Received and read oder from: ',ifirst,'<',DB_type(1:5),'>'

      DB_order = ' '
!< Change 1.2
!     DB_order = 'arr=j,ak=st,ty='//DB_type(1:5)
      DB_order = 'arr=j,ak=st,test=n,unload=n,ty='//DB_type(1:5)
!>

      IF( MPE_dbg_level > 1 ) THEN
      write(0,*) me_io,' send to database: '
      write(0,*) DB_order(1:len_trim(DB_order))
      ENDIF

      call csodb_arrout_req( DB_order, request,&
           DB_hits,error, buff_len, ilenc, data)

      db_req_count = db_req_count + 1 

      IF( db_req_count == db_req_err) error = 6 

      !write(0,*) 'arrout req on : ', me_io, ' with err: ', error 

      IF(DB_do_backup) THEN
        !write(0,*) 'save: ', ilenc, ' in GRIB_SAVE'
        CALL GRIB_SAVE( trim(current_datname), me_io, data, ilenc, err)

        IF( err /= 0 ) THEN
          write(0,*) 'error in GRIB_SAVE', me_io
        ENDIF

      ENDIF

   ENDIF

ENDIF

nudat = nudatc

END SUBROUTINE DEV_Write

!--------------------------------------------------------------------

SUBROUTINE DEV_Wait( request, error)

IMPLICIT NONE

INTEGER, INTENT(INOUT)  :: request
INTEGER, INTENT(OUT)    :: error

error = 0

SELECT CASE ( out_device )

   CASE ( 'bank' ) 
     CALL csodb_arrout_end( request, DB_hits, error )

END SELECT

END SUBROUTINE DEV_Wait

!--------------------------------------------------------------------

SUBROUTINE mpe_io_complete( )

!--------------------------------------------------------------------------
!
! $Revision: 1.1.1.6 $
! $DATE$
!
! Purpose          : completes a mpe_io_write call
!
! Modules          :
! 
! Module parameters:
!
!--------------------------------------------------------------------------

  IMPLICIT NONE

!--------------------------------------------------------------------------
!
! Parameters       :
!
!--------------------------------------------------------------------------
!
! Local parameters :

INTEGER :: ilen, n_io, pe_target, tag_target, ierr, data   

INTEGER :: end_fld, rest_fld

INTEGER :: max_tag, min_tag, status(MPI_STATUS_SIZE)

INTEGER :: pe_source, tag_source

INTEGER :: isend_status(MPI_STATUS_SIZE)
!
!--------------------------------------------------------------------------
!
! Error checks     :



!
!--------------------------------------------------------------------------

ilen  =  0  ! send "end" messages to end receive loop

!write(0,*) 'jump into  complete ', me_io, me_compute

IF( mpe_is_compute()) THEN

   CALL mpe_io_next_snd( pe_target, tag_target, me_compute)

   !write(0,*) 'sending end message ', tag_target, me_compute

   !don't send to myself

   if( io_config == 1 .AND. pe_target == me_compute ) &
                          pe_target = MPI_PROC_NULL

   CALL MPI_SEND( data, ilen, MPI_BYTE, pe_target, tag_write_end &
        , MPE_COMM_CO2IO, ierr) 

   !CALL MPI_Wait( isend_request, isend_status, ierr)

ENDIF

IF( mpe_is_io() .AND. .NOT.got_end_signal ) THEN

   CALL mpe_io_next_recv( pe_source, tag_source, me_io )

   CALL MPI_ALLREDUCE( tag_source, min_tag, 1, MPI_INTEGER &
        , MPI_MIN                 &
        , MPE_COMM_IO, ierr)

   max_tag = min_tag + num_compute 
   
   DO WHILE ( tag_source .lt. max_tag ) 

      !write(0,*) max_tag, tag_source, pe_source

      ! don't receive from myself

      if( io_config == 1 .AND. pe_source == me_compute ) &
           pe_source = MPI_PROC_NULL

      
      CALL MPI_RECV( data, 1, MPI_BYTE   &
           , pe_source, tag_write_end       &
           , MPE_COMM_CO2IO, status, ierr)

      !write(0,*) 'receive done'

      CALL mpe_io_next_recv( pe_source, tag_source, me_io )

   ENDDO
   


ENDIF


END SUBROUTINE mpe_io_complete

!-------------------------------------------------------------------------

!+ reads from device 

SUBROUTINE mpe_io_read( nudat, data, ilen, ilfd, error )

!--------------------------------------------------------------------------
!
! $Revision: 1.1.1.6 $
! $DATE$
!
! Purpose          : Read data on I/O PEs and 
!                    send it to compute PEs
!
! Modules          :
! 
! Module parameters:
!
!--------------------------------------------------------------------------

  IMPLICIT NONE

!--------------------------------------------------------------------------
!
! Parameters       :
!

INTEGER, INTENT(IN)       :: nudat    &! file descriptor
                           , ilfd      ! buffer size in intgribf

INTEGER, INTENT(OUT)      :: ilen      ! number of Bytes read
                     
INTEGER  (KIND = gribf_kind)         &
        , INTENT(INOUT)   :: data(*)   ! array for read data 

INTEGER , INTENT(OUT), OPTIONAL     :: error

!--------------------------------------------------------------------------
!
!
! Local Parameters

INTEGER  (KIND= gribc_kind ) :: nudatc, ilfdc, ilenc, ierrc       
                           ! corresponding variables for C-routines


INTEGER ierr, status(MPI_STATUS_SIZE)

INTEGER, save  :: totalbytes = 0


!--------------------------------------------------------------------------

IF (PRESENT(error)) THEN
  error = 0
ENDIF

!if( is_io_pe ) then
!WRITE(0,*) me_io,' into mpe_io_read , curr. hist /totalbytes = '&
!                ,DB_hits,totalbytes,' in_device = ',TRIM(in_device)
!endif


DO
   ilen = 0  

   IF( in_device == 'file' .AND. me_io == 0 ) THEN

      ilfdc = ilfd

#ifdef GRIBDWD
      call cuegin( nudat, ilfdc, data, ilenc, ierrc )

      ! write(0,*) 'cuegin done on, len,err = ', ilenc,' ',ierrc,ilfdc
#endif

      ilen = ilenc

   ENDIF

   IF( in_device == 'bank' .AND. me_io == 0 ) THEN

      ilen = 0

      IF( totalbytes == 0 .OR. DB_hits > 1 ) THEN


!     write(0,*) 'Read from database ',DB_order

      call csodb_arrin_req(DB_order, DB_request, DB_err )


      call csodb_arrin_end(DB_request, DB_hits, DB_err, ilfd, &
                                     ilen, data)

      totalbytes = totalbytes + ilen

      if( MPE_dbg_level > 1 ) THEN
       write(0,*) 'Read from database with DB_hits, ilen , totalbytes ',&
                                           DB_hits,' ', ilen ,' ', totalbytes 
      ENDIF

      IF( DB_hits < 1 ) THEN

        write(0,*) 'Read from Database failed' 

      ENDIF 

      ENDIF


   ENDIF

   IF( io_config == 2 ) THEN

      IF( me_io == 0 ) THEN  ! send data to me_compute = 0

         CALL MPI_SEND( data, ilen, MPI_BYTE, 0, tag_read &
              , MPE_COMM_CO2IO, ierr)


      ENDIF

      IF ( me_compute == 0 ) THEN


         IF ( .NOT. end_read ) THEN

           CALL MPI_RECV( data, ilfd, MPI_BYTE, 0, tag_read &
              , MPE_COMM_CO2IO, status, ierr)
         
           CALL MPI_GET_COUNT( status, MPI_BYTE, ilen, ierr) 

           !write(0, *) ilen , ' Bytes received on compute node '

           if( ilen == 0 ) end_read = .TRUE.

        ELSE
           ilen = 0
        ENDIF 

      ENDIF

  ENDIF 

  IF( is_compute_pe .OR.  ilen == 0  ) EXIT
  
ENDDO


END SUBROUTINE mpe_io_read

!-------------------------------------------------------------------------

SUBROUTINE mpe_io_close( nudat, ierror )

!--------------------------------------------------------------------------
!
! $Revision: 1.1.1.6 $
! $DATE$
!
! Purpose          :
!
! Modules          :
! 
! Module parameters:
!
!--------------------------------------------------------------------------

  IMPLICIT NONE

!--------------------------------------------------------------------------
!
! Parameters       :
!

INTEGER, INTENT(IN)  :: nudat  ! filedescriptor
INTEGER, INTENT(OUT) :: ierror ! error parameter

!--------------------------------------------------------------------------
!
! Local parameters :
!

INTEGER  (KIND= gribc_kind ) :: nudatc &! file descriptor for C-routine      
                               ,ierrc   ! error parameter for C-routine

!--------------------------------------------------------------------------
!
! Error checks     :
!
!--------------------------------------------------------------------------

ierror = 0
ierrc  = 0

IF ( is_io_pe ) THEN

#ifdef GRIBDWD
   nudatc = nudat
   IF( out_device == 'file' .AND. current_mode == 'w' ) &
        CALL cclose(nudatc,'exi',ierrc)

   IF( in_device  == 'file' .AND. current_mode == 'r' .AND. me_io == 0) &
        CALL cclose(nudatc,'exi',ierrc)

   ierror = ierrc
#endif

ENDIF

IF( ierror == 0 ) open_count = open_count - 1

current_datname =  ''

END SUBROUTINE mpe_io_close

!-------------------------------------------------------------------------

SUBROUTINE mpe_io_node()

!--------------------------------------------------------------------------
!
! $Revision: 1.1.1.6 $
! $DATE$
!
! Purpose          : program for "dumb" io nodes
!
! Modules          :
! 
! Module parameters:
!
!--------------------------------------------------------------------------

  IMPLICIT NONE

!--------------------------------------------------------------------------
!
! Parameters       :
!
!--------------------------------------------------------------------------
!
! Local parameters :
!

   LOGICAL                      :: stop_io
   INTEGER                      :: nudat, len_dummy
   CHARACTER(LEN=1)             :: datname_dummy, mode_dummy, db_order_dummy
   INTEGER, parameter           :: buff_len = BUFF_length*gribf_kind
   INTEGER  (KIND= gribf_kind ) :: data( buff_len )   

   CHARACTER (LEN=5)            :: dummy_type

!--------------------------------------------------------------------------
!
! Error checks     :

   INTEGER                      :: error

!
!--------------------------------------------------------------------------

   stop_io = .FALSE.

   DO

      IF( MPE_dbg_level > 1 ) THEN
      write(0,*) 'before open in mpe_io_node #', me_io
      ENDIF

      CALL mpe_io_open(  nudat, datname_dummy&
                              , db_order_dummy, mode_dummy, stop_io )

      IF( MPE_dbg_level > 1 ) THEN
      write(0,*) 'open on io node done ', me_io
      ENDIF

      if( stop_io ) EXIT

      IF (current_mode == 'w' ) THEN

         IF( MPE_dbg_level > 1 ) THEN
         write(0,*) 'before write on node done', me_io
         ENDIF

         CALL mpe_io_write( nudat, data, len_dummy, buff_len, dummy_type)    

         IF( MPE_dbg_level > 1 ) THEN
         write(0,*) 'write on node done', me_io
         ENDIF

      ENDIF

      IF( current_mode == 'r' ) THEN

         !write(0,*) 'before read on io-node '

         CALL mpe_io_read( nudat, data, len_dummy, buff_len)
 
         !write(0,*) 'read on node done'

      ENDIF 

      IF( MPE_dbg_level > 1 ) THEN
      write(0,*) 'before close on io node done ', me_io
      ENDIF

      CALL mpe_io_close( nudat, error)

      IF( MPE_dbg_level > 1 ) THEN
      write(0,*) 'close on io node done ', me_io
      ENDIF

   ENDDO

  IF( MPE_dbg_level > 1 ) THEN
   write(0,*) 'io-node #', me_io, ' has been shut down'
  ENDIF

END SUBROUTINE mpe_io_node

!-------------------------------------------------------------------------

    SUBROUTINE mpe_io_shutdown()

!--------------------------------------------------------------------------
!
! $Revision: 1.1.1.6 $
! $DATE$
!
! Purpose          : Shutdown io nodes
!
! Modules          :
! 
! Module parameters:
!
!--------------------------------------------------------------------------

  IMPLICIT NONE

!--------------------------------------------------------------------------
!
! Parameters       :
!
!--------------------------------------------------------------------------
!
! Local parameters :
!
   INTEGER          :: tag, dummy(1), ierr

   INTEGER          :: i                   ! loop index


   INTEGER          :: isend_status(MPI_STATUS_SIZE, send_max_request)
!
!--------------------------------------------------------------------------
!
! Error checks     :
!
!--------------------------------------------------------------------------

   dummy(1) = -1

   IF ( io_config == 2 ) THEN

      IF ( is_compute_pe ) THEN

        DO i = 0, num_io - 1
           IF( me_compute == i ) THEN

             CALL MPI_SEND( dummy, 1, MPI_INTEGER, i, tag_open, &
                            MPE_COMM_CO2IO, ierr ) 
           ENDIF
        ENDDO

         
      ENDIF

   ENDIF

   IF( is_compute_pe ) THEN

     CALL MPI_Waitall(send_max_request, send_requests, isend_status, ierr)

     DEALLOCATE( send_buffer )

   ENDIF

   IF( is_io_pe ) THEN

     deallocate( BUFF_space , STAT = ierr )
     deallocate( BUFF_len_curr )   

   ENDIF

   IF( out_device == 'bank' .or. in_device == 'bank' ) THEN
! Shutdown database
     IF( is_io_pe ) THEN

       CALL csodb_close(DB_err)
!       CALL csodb_arr_close()

      print*, 'Bank closed with error ', DB_err 

       CALL DB_statistics( MPE_COMM_IO)


     ENDIF
   ENDIF


END SUBROUTINE mpe_io_shutdown

!--------------------------------------------------------------------------

SUBROUTINE DB_statistics(DB_Comm)

IMPLICIT NONE

INTEGER:: DB_Comm

REAL:: output_time, input_time, env_time
REAL:: OMbytes, IMbytes, msg(6)
INTEGER:: inp_calls, inp_succ, outp_calls, outp_succ, env_calls, env_succ
INTEGER:: rank, size, mpi_err, imsg(12)

      call DB_stat_out   &
      (input_time, IMBytes, inp_calls, inp_succ, &
       output_time, OMBytes, outp_calls, outp_succ, &
       env_time, env_calls, env_succ)


If( DB_Comm /= MPI_COMM_NULL ) THEN

      CALL MPI_COMM_SIZE( DB_Comm, size, mpi_err )
      CALL MPI_COMM_RANK( DB_Comm, rank, mpi_err )

      msg(1) = input_time
      msg(2) = output_time
      msg(3) = env_time

      CALL MPI_REDUCE(msg(1), msg(4), 3, MPI_REAL,MPI_MAX, 0, DB_Comm, mpi_err)

      input_time = msg(4)
      output_time = msg(5)
      env_time = msg(6)

      imsg(1) = inp_calls
      imsg(2) = inp_succ
      imsg(3) = outp_calls
      imsg(4) = outp_succ
      imsg(5) = env_calls
      imsg(6) = env_succ


      CALL MPI_REDUCE(imsg(1), imsg(7), 6, MPI_INTEGER&
                             ,MPI_SUM, 0, DB_Comm, mpi_err)

      inp_calls = imsg(7)
      inp_succ  = imsg(8)
      outp_calls = imsg(9)
      outp_succ  = imsg(10)
      env_calls = imsg(11)
      env_succ  = imsg(12)

      msg(1) = IMBytes
      msg(2) = OMbytes

      CALL MPI_REDUCE(msg(1), msg(3), 2, MPI_REAL,MPI_SUM, 0, DB_Comm, mpi_err)

      IMBytes = msg(3)
      OMBytes = msg(4)

      IF( rank == 0 .and. max(DB_OUT_SOCKETS,DB_IN_SOCKETS)>0 ) THEN

      PRINT 1,' '
      PRINT 1,'========================================================'
      PRINT 1,' '
      PRINT 1,' OVERALL DATA BASE STATISTICS: '
      PRINT 1,' '
      PRINT 1,' #IO-PE:                   ',num_io
      PRINT 1,' '
      IF( DB_OUT_SOCKETS > 0 ) THEN

      PRINT 1,' #OUTPUT SOCKETS           ',DB_OUT_SOCKETS*num_io
      PRINT 1,' '
      PRINT 1,' Data base output requests: '
      PRINT 1,' '

      PRINT 1,' Number:                   ',outp_calls
      PRINT 1,' Unsuccessful:             ',outp_calls-outp_succ
      PRINT 2,' Time:                     ',output_time
      PRINT 2,' MBytes:                   ',OMbytes
      IF( output_time > 0. ) THEN
      PRINT 2,' OUTPUT MBytes/s:          ',OMbytes/output_time
      ENDIF

      ENDIF
      IF( DB_IN_SOCKETS > 0 ) THEN

      PRINT 1,' '
      PRINT 1,' Data base input requests: '
      PRINT 1,' '

      PRINT 1,' Number:                   ',inp_calls
      PRINT 1,' Unsuccessful:             ',inp_calls-inp_succ
      PRINT 2,' Time:                     ',input_time
      PRINT 2,' MBytes:                   ',IMbytes
      IF( input_time > 0. ) THEN
      PRINT 2,' INPUT MBytes/s:           ',IMbytes/input_time
      ENDIF

      ENDIF

      PRINT 1,' '
      PRINT 1,' Data base open/close/commit requests: '
      PRINT 1,' '

      PRINT 1,' Number:                   ',env_calls
      PRINT 1,' Unsuccessful:             ',env_calls-env_succ
      PRINT 2,' Time:                     ',env_time


      ENDIF

! Statistics from gribsave library:
      IF( is_io_pe ) THEN

        CALL IO_STAT(imsg(1), msg(1))
        CALL MPI_REDUCE(imsg(1), size, 1, &
              MPI_INTEGER,MPI_SUM, 0, DB_Comm, mpi_err)
        CALL MPI_REDUCE(msg(1), output_time, 1, &
              MPI_INTEGER,MPI_MAX, 0, DB_Comm, mpi_err)

        IF( me_io == 0 ) THEN

        PRINT 1,' '
        PRINT 2,' Backup mechanism  time:   ',output_time
        PRINT 2,' Backup MBytes:            ',1.e-6*size
        IF( output_time > 0. ) THEN

         PRINT 2,' Backup MBytes/s:          ',1.e-6*size/output_time

        ENDIF
        PRINT 1,' '
        PRINT 1,'========================================================'
        PRINT 1,' '

        ENDIF

      ENDIF


ENDIF

1 FORMAT(A,I12)
2 FORMAT(A,f12.2)

END SUBROUTINE DB_statistics


!------------------------------------------------------------------

! Error handling routines:

!------------------------------------------------------------------

SUBROUTINE DB_err_handler( error_in, collective, error_out)

IMPLICIT NONE

INTEGER, INTENT(IN)    :: error_in
LOGICAL, INTENT(IN)    :: collective
INTEGER, INTENT(OUT)   :: error_out

INTEGER                :: err_snd, err_max

error_out = 0


IF( is_io_pe ) THEN

  err_snd = abs( error_in)
  err_max = abs( error_in)

  IF ( collective ) THEN

    CALL MPI_ALLREDUCE( err_snd, err_max, 1, MPI_INTEGER, MPI_MAX&
                  , MPE_COMM_IO, error_out) 

  ENDIF

  IF( out_device == 'bank' .AND. err_max /= 0 ) THEN

    write(0,*) ' BEGIN ERROR HANDLING on I/O PE ', me_io
    write(0,*) ' SWITCHING TO FILE-IO '

    ori_out_device = 'bank'
    error_handling = .TRUE.
    out_device     = 'file'
    error_out      =  100

  ENDIF

ENDIF

END SUBROUTINE DB_err_handler

!--------------------------------------------------------------

SUBROUTINE DB_err_handler_write( error_in, nudat, buff_len&
                                         , pending, error_out )
IMPLICIT NONE

INTEGER, INTENT(IN)    :: error_in
INTEGER, INTENT(INOUT) :: nudat
INTEGER, INTENT(IN)    :: buff_len
LOGICAL, INTENT(INOUT) :: pending(*)
INTEGER, INTENT(OUT)   :: error_out

INTEGER                :: ilen, dat_len,i,request, nudat_tmp
INTEGER                :: finish_req, grib_length, status
CHARACTER (LEN=max_name_len) :: name_tmp
INTEGER, parameter           :: max_files = 128

nudat_tmp = 0

IF( error_in /= 0 ) THEN

   CALL DB_err_handler( error_in,.FALSE., error_out)

   IF( DB_do_backup ) THEN

      CALL GRIB_SAVE_CLOSE( me_io, error_out)

      IF( error_out /= 0 ) THEN
        write(0,*) 'error in GRIB_SAVE_CLOSE on me_io', me_io
      ENDIF

      DO i = 1, max_files

        ! write filenames to "current_datname" this string is 
        ! a module variable and is used by DB_err_handler_write

         name_tmp = ' '

        CALL GRIB_SAVE_FILENAME( name_tmp, me_io, i, status)

        !get rid of trailing ^ 
        name_tmp = name_tmp(1:len_trim(name_tmp)-1)

        IF ( status /= 0 ) EXIT

        IF( MPE_dbg_level > 1 ) THEN
        write(0,*) 'GRIB_SAVE name = ', name_tmp
        ENDIF

        dat_len        = len_trim( name_tmp )

        CALL DEV_Open( nudat_tmp, trim(name_tmp), dat_len&
           , 'w  ', error_out)

        request = 1

        DO

         CALL GRIB_RETRIEVE(trim(name_tmp)   &! backup name
                         ,me_io                   &! index
                         ,BUFF_space(1,1)         &! buffer
                         ,buff_len                &! buffer dimension
                         ,grib_length             &! length to be written
                         ,error_out)

         IF( error_out /= 0 ) THEN
           write(0,*) 'error in retrieve on me_io', me_io
         ENDIF

         IF( MPE_dbg_level > 1 ) THEN
         write(0,*) 'retrieve ', grib_length,' bytes'  
         ENDIF

 
         IF( grib_length <= 0 ) EXIT

         CALL DEV_Write( nudat_tmp, BUFF_space(1,request)&
                                  , buff_len, grib_length&
                                  , request, error_out)

        ENDDO

        IF( trim(name_tmp) /= trim(current_datname) ) THEN
          CALL mpe_io_close( nudat_tmp, error_out )
        ELSE
          write(0,*) 'going on with file I/O ',trim(current_datname) 
          nudat     = nudat_tmp
          nudat_tmp = 0
        ENDIF

      ENDDO 
    ENDIF

    DO i = 1, BUFF_max_requests
        pending(i) = .FALSE.
    ENDDO

ENDIF

END SUBROUTINE DB_err_handler_write 

END MODULE mpe_io


