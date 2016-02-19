!+ Utilitiy module for the parallel environment
!==============================================================================

MODULE  environment

!==============================================================================
!
! Description:
!   This module provides routines that have to deal with defining the 
!   environment of the program. In case of running on parallel computers
!   with distributed memory, the MPI library is used as message passing
!   library. If running on a sequential machine without MPI, dummies are
!   provided for linking the program.
!
!   Routines (module procedures) currently contained:
!
!     - collapse:
!         changes the computational indices to allow for a loop collapsing
!         on NEC-SX architectures
!     - extend_field:
!         extends a field with additional rows and columns
!     - init_environment:
!         initializes the MP-library, number of processors present and 
!         rank of own process.
!     - final_environment:
!         cleanup of the message passing library at the end of the program.
!     - get_free_unit:
!         for ascii file handling: this routine determines the next free unit
!     - release_unit:
!         for ascii file handling: a unit number can be set free again
!     - model_abort:   
!         one process stops the whole parallel program
!     - init_procgrid:
!         creates the virtual topology of the processor grid.
!     - exchg_boundaries:
!         performs the data exchange between neighbouring processors
!     - putbuf:
!         puts values for data exchange into the sending buffer
!     - getbuf:
!         gets values for data exchange from the sending buffer
!     - comm_barrier:
!         sets a synchronization point
!     - setup_data_type
!         defines and allocates MPI-Datatypes for the boundary exchange
!       
! Method: 
!   Calls to the message-passing library MPI
!
!
! Current Code Owner: DWD, Ulrich Schaettler
!  phone:  +49  69  8062 2739
!  fax:    +49  69  8062 3721
!  email:  ulrich.schaettler@dwd.de
!
! History:
! Version    Date       Name
! ---------- ---------- ----
! 1.1        1998/03/11 Ulrich Schaettler
!  Initial release
! 1.8        1998/08/03 Ulrich Schaettler
!  Eliminated dependencies on grib parameters.
! 1.10       1998/09/29 Ulrich Schaettler
!  Eliminated Parameters nmax* from Use lists.
! 1.29       1999/05/11 Ulrich Schaettler
!  Adaptations to use this module also in GME2LM
! 1.34       1999/12/10 Ulrich Schaettler
!  Changes in putbuf and getbuf for vectorization and changed model_abort
! 2.8        2001/07/06 Ulrich Schaettler
!  Introduced lreorder as parameter for subroutine init_procgrid
! 2.13       2002/01/18 Ulrich Schaettler
!  Correction of error code messaging.
! 2.17       2002/05/08 Ulrich Schaettler
!  Additions to perform the communications for I/O in irealgrib-format
! 3.5        2003/09/02 Ulrich Schaettler
!  Optimizations in the routine for the boundary exchange
! 3.7        2004/02/18 Ulrich Schaettler
!  Increased number of optional variables for exchg_boundaries.
!  Introduced new routines get_free_unit and release_unit for ascii-file
!  handling; added closing of ascii-files to model_abort
! 3.13       2004/12/03 Ulrich Schaettler
!  Adaptations to version in new INT2LM.
!  Determine domain boundaries for ncomm_type=3 in exchg_boundaries
! 3.14       2005/01/25 Ulrich Schaettler
!  Bug fix in the routine extend_field
! 3.21       2006/12/04 Ulrich Schaettler
!  Eliminated some unused variables
! V4_9         2009/07/16 Ulrich Schaettler
!  Implemented SR collapse to change loop indices for collapsing on NEC-SX
! V4_11        2009/11/30 Ulrich Schaettler
!  Adaptation to run the MPI data types also for 64-bit environments
!  (Based on the work by C. Pospiech, IBM)
! V4_13        2010/05/11 Michael Gertz
!  Adaptions to SVN
! V4_17        2011/02/24 Ulrich Blahak
!  Adapted interface of exchg_boundaries; corrected kzdims(1:20) -> kzdims(1:24).
!  Eliminated my_peri_neigh. 
!  Included 1-processor periodic/2dim exchanges
!    into exchg_boundaries, which now has to be called also by 1-processor runs.
!    This greatly simplifies the implementation of periodic and other BCs.
!  Added lperi_x, lperi_y, l2dim to calling list of extend_field.
! V4_19        2011/08/01 Ulrich Schaettler
!  Introduced conditional compilation for GRIBDWD: necessary here for fsleep
! V4_20        2011/08/31 Ulrich Schaettler
!  Implemented interface to OASIS coupler using conditional compilation with -DCOUP_OAS
!   (by CLM Community)
! V4_23        2012/05/10 Ulrich Schaettler
!  Unification with latest INT2LM Version 1.19
! V4_25        2012/09/28 Carlos Osuna, Ulrich Blahak
!  Introduce asynchronous netcdf in the distribution of PE into
!  compute PE and I/O PE.
! Ulrich Blahak: some bug fixes
!  SR init_procgrid: added correct setting of my_cart_neigh in case of serial
!                    runs with periodic BCs.
!  SR exchg_boundaries:
!  - in error-case nlines > nboundlines, changed RETURN to MODEL_ABORT() for user safety.
!  - for periodic BCs, exchange all nboundlines regardless of the parameter nlines,
!    to get "beautiful" periodic output. In this case, performance is usually not
!    so important to the users.
! V4_27        2013/03/19 Ulrich Blahak, Astrid Kerkweg
!  bugfix of bugfix concerning serial periodic BCs: corrected lperi_x to lperi_y (UB)
!  MESSy interface introduced: here interface for MMD (Multi-Model-Driver) (AK)
! V4_28        2013/07/12 Ulrich Schaettler
!  Extended interface for init_environment, to pass integer type imp_integ_ga
!    for grib_api integer (for length of message in bytes)
!  Boundary exchange in extend_field only for num_compute > 1
! V4_30        2013/11/08 Ulrich Schaettler
!  Corrected usage of uninitialized variable yzerrmsg in SR exchg_boundaries
! V5_1         2014-11-28 Ulrich Schaettler, Oliver Fuhrer
!  Adaptations to latest INT2LM version
!  Replaced ireals by wp (working precision) (OF)
!  Added parameters imp_single and imp_double (SP and DP REAL kind for MPI) (OF)
!  Replaced call to fsleep with dosleep from utilities.f90
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!==============================================================================
!
! Declarations:
!
! Modules used:
!------------------------------------------------------------------------------

USE data_parameters , ONLY :   &
  wp,        & ! KIND-type parameters for real variables
  iintegers, & ! kind-type parameter for "normal" integer variables
  int_ga       ! integer precision for grib_api: length of message in bytes

USE utilities,        ONLY: dosleep

#ifdef MESSYMMD
  USE mmd_handle_communicator, ONLY: MMD_get_model_communicator &
                                   , MMD_Print_Error_Message    &
                                   , MMD_FreeMem_communicator   &
                                   , MMD_STATUS_OK
#endif

!------------------------------------------------------------------------------

IMPLICIT NONE

!------------------------------------------------------------------------------

! include statements
INCLUDE "mpif.h"

#ifdef MESSYMMD
  INTEGER(KIND=iintegers),  SAVE  :: MMD_comm_world  ! mz_kk_20081107
#endif

!------------------------------------------------------------------------------

! Global variables
  INTEGER(KIND=iintegers), PARAMETER       ::                         &
    iexchg_MPI_type_len    = 200
       ! length of global vector iexchg_MPI_type used in environment.f90


  INTEGER   (KIND=iintegers   ) ::  &
    iexchg_MPI_types(iexchg_MPI_type_len) = MPI_DATATYPE_NULL,        &
                  ! List of MPI data types used in exchg_datatypes
                  ! routine. If set to MPI_DATATYPE_NULL, the
                  ! corresponding entry has not been properly set up.
    iexchg_counts(iexchg_MPI_type_len)
                  ! List of counts of these data types
                  ! used in exchg_datatypes routine.
                  ! Has meaningful value only if the corresponding
                  ! vector element of iexchg_MPI_types is not
                  ! MPI_DATATYPE_NULL.

! Variables for ASCII file handling
  INTEGER(KIND=iintegers)                  ::                         &
    iunit_start = 21,                     &
    iunit_end   = 100,                    &
    iunit_table(21:100) = 0_iintegers

!------------------------------------------------------------------------------

CONTAINS
! Define procedures contained in this module.

!==============================================================================
!+ Subroutine that extends a given field with additional rows and columns
!------------------------------------------------------------------------------

SUBROUTINE collapse (onoff, ie, istart, iend, jstart, jend)

!------------------------------------------------------------------------------
!
!  This subroutine changes the computational indices to 1, ie, 1, je instead
!  of using istartpar, iendpar, jstartpar, jendpar. This enables a loop
!  collapsing on NEC-SX machines.
!
!------------------------------------------------------------------------------

LOGICAL,                  INTENT(IN)    :: onoff
INTEGER (KIND=iintegers), INTENT(IN)    :: ie

INTEGER (KIND=iintegers), INTENT(INOUT) :: istart,iend,jstart,jend

! Local variable that have to be saved for the next call, to reset the indices
! again
INTEGER (KIND=iintegers),SAVE :: istartmem,iendmem,jstartmem,jendmem

IF (onoff) THEN
  istartmem=istart
  iendmem=iend
  jstartmem=jstart
  jendmem=jend
  IF ((istart == 1) .AND. (iend == ie)) THEN
    jend=jstart
    iend=ie*(jendmem+1-jstartmem)
  END IF
ELSE
  istart=istartmem
  iend=iendmem
  jstart=jstartmem
  jend=jendmem 
END IF

END SUBROUTINE collapse

!==============================================================================
!+ Subroutine that extends a given field with additional rows and columns
!------------------------------------------------------------------------------

SUBROUTINE extend_field ( field_in,  ie_in,  je_in,                          &
                          field_out, ie_out, je_out, kedim,                  &
                          nextlines, nboundlines, jstartpar_ext, jendpar_ext,&
                          sendbuf, isblen, imp_reals, icomm, neighbors,      &
                          num_compute, lperi_x, lperi_y, l2dim)

!------------------------------------------------------------------------------
!
! Description:
!   This routine extends field_in with nextlines additional rows and columns
!   around. At the borders of the total domain, all values are set to the
!   values of the neighboring row/column of field_in. If field_in has extra
!   boundary lines (nboundlines in LM/INT2LM), values for these boundary
!   lines must be given at the borders of the total domain.
!
!   At the borders to the other subdomains, values are determined through
!   a boundary exchange.
!
! Method:
!
!------------------------------------------------------------------------------

! Subroutine arguments:
! ---------------------
  INTEGER (KIND=iintegers), INTENT (IN) ::  &
    ie_in,  je_in,         & ! horizontal dimensions of field_in
    ie_out, je_out, kedim, & ! dimensions of field_out
    nextlines,             & ! number of additional rows/columns at each side
                             ! for the extension
    nboundlines,           & ! number of boundary lines for field_in from
                             ! the decomposition
    isblen,                & ! length of sendbuffer
    imp_reals,             & ! determines the REAL type used
    icomm,                 & ! communicator for virtual cartesian topology
    neighbors(4),          & ! process-id's of the neighbors in the grid
    num_compute              ! number of compute processors

  INTEGER (KIND=iintegers), INTENT(OUT) ::  &
    jstartpar_ext, jendpar_ext ! start- and end-indices for the boundary
                               ! exchange for the extended field

  REAL (KIND=wp),           INTENT (INOUT) ::    &
    sendbuf (isblen, 8)        ! send buffer

  REAL    (KIND=wp),        INTENT (IN) ::  &
    field_in (ie_in, je_in, kedim)

  REAL    (KIND=wp),        INTENT(OUT) ::  &
    field_out(ie_out,je_out,kedim)

  LOGICAL,                  INTENT(IN)  ::  &
    lperi_x, lperi_y, l2dim

! Local scalars:
! -------------
INTEGER (KIND=iintegers) ::  &
     izerror, i, j, k, nbdext, kzdims(24)

CHARACTER (LEN=200)        ::  &
  yzerrmsg    ! error message for error handling

!------------------------------------------------------------------------------

  izerror   = 0
  yzerrmsg  = '   '

  nbdext    = nboundlines + nextlines

  field_out(:,:,:) = 0.0_wp

  DO k = 1, kedim
    DO j = 1, je_in
      DO i = 1, ie_in
        field_out(i+nextlines,j+nextlines,k) = field_in(i,j,k)
      END DO
    END DO
  END DO

  ! set values for the boundary lines in the extended frame
  ! western border
  IF (neighbors(1) == -1) THEN
    DO k = 1, kedim
      DO j = 1, je_in
        DO i = 1, nextlines
          field_out(i,j+nextlines,k) = field_in(1,j,k)
        ENDDO
      ENDDO
    ENDDO
  ENDIF

  ! eastern border
  IF (neighbors(3) == -1) THEN
    DO k = 1, kedim
      DO j = 1, je_in
        DO i = ie_in+1, ie_in+nextlines
          field_out(i+nextlines,j+nextlines,k) = field_in(ie_in,j,k)
        ENDDO
      ENDDO
    ENDDO
  ENDIF

  ! southern border
  IF (neighbors(4) == -1) THEN
    jstartpar_ext = 1
    DO k = 1, kedim
      DO j = 1, nextlines
        DO i = 1, ie_in
          field_out(i+nextlines,j,k) = field_in(i,1,k)
        ENDDO
      ENDDO
    ENDDO
  ELSE
    jstartpar_ext = 1 + nboundlines + nextlines
  ENDIF

  ! northern border
  IF (neighbors(2) == -1) THEN
    jendpar_ext = je_out
    DO k = 1, kedim
      DO j = je_in+1, je_in+nextlines
        DO i = 1, ie_in
          field_out(i+nextlines,j+nextlines,k) = field_in(i,je_in,k)
        ENDDO
      ENDDO
    ENDDO
  ELSE
    jendpar_ext = je_in - nboundlines + nextlines
  ENDIF


  ! southwest corner
  IF ( (neighbors(4) == -1) .AND. (neighbors(1) == -1) ) THEN
    DO k = 1, kedim
      DO j = 1, nextlines
        DO i = 1, nextlines
          field_out(i,j,k) = field_in(1,1,k)
        ENDDO
      ENDDO
    ENDDO
  ENDIF

  ! northwest corner
  IF ( (neighbors(2) == -1) .AND. (neighbors(1) == -1) ) THEN
    DO k = 1, kedim
      DO j = je_in+1, je_in+nextlines
        DO i = 1, nextlines
          field_out(i,j+nextlines,k) = field_in(1,je_in,k)
        ENDDO
      ENDDO
    ENDDO
  ENDIF

  ! southeast corner
  IF ( (neighbors(4) == -1) .AND. (neighbors(3) == -1) ) THEN
    DO k = 1, kedim
      DO j = 1, nextlines
        DO i = ie_in+1, ie_in+nextlines
          field_out(i+nextlines,j,k) = field_in(ie_in,1,k)
        ENDDO
      ENDDO
    ENDDO
  ENDIF

  ! northeast corner
  IF ( (neighbors(2) == -1) .AND. (neighbors(3) == -1) ) THEN
    DO k = 1, kedim
      DO j = je_in+1, je_in+nextlines
        DO i = ie_in+1, ie_in+nextlines
          field_out(i+nextlines,j+nextlines,k) = field_in(ie_in,je_in,k)
        ENDDO
      ENDDO
    ENDDO
  ENDIF

  IF (num_compute > 1) THEN
    kzdims(1:24) = (/ kedim,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /)
    CALL exchg_boundaries                                                      &
      ( 0, sendbuf, isblen, imp_reals, icomm, num_compute, ie_out,             &
        je_out, kzdims, jstartpar_ext, jendpar_ext, nbdext, nbdext, neighbors, &
        lperi_x, lperi_y, l2dim, 2000,  .FALSE., 1, izerror, yzerrmsg,         &
        field_out )
  ENDIF

END SUBROUTINE extend_field

!==============================================================================
!+ Subroutine that initializes variables for the parallel environment
!------------------------------------------------------------------------------

SUBROUTINE init_environment (nproc, my_world_id, icomm_world, igroup_world, &
                             imp_int, imp_real, imp_single, imp_double,     &
                             imp_grib, imp_byte, imp_char,                  &
                             imp_logical, imp_intga, iexch_req,             &
                             irealgrib, yerrmsg, ierror                     )

!------------------------------------------------------------------------------
!
! Description:
!   This subroutine initializes the environment of the program.
!   MPI, the nodes and several organizational variables for the 
!   parallelization are also initialized. These are
!     - nproc:        the number of processors present for running the program
!     - igroup_world: identifier for the global group
!     - my_world_id:  the rank of this processor in the global group
!     - iexch_req:    identifier for send- and recv-requests in the data
!                     exchange between subdomains.
!   In the sequential program, these variables are set accordingly.
!
!   The following variables exist also as MPI variables. But because not every
!   routine should include the header file mpif.h, they are also declared as
!   model variables:
!     - imp_integers: integer data type used
!     - imp_reals:    real data type used
!     - icomm_world:  identifier for the communicator of the global group
!
! Method:
!   MPI-calls and error handling.
!
!------------------------------------------------------------------------------

#ifdef COUP_OAS_COS
    USE oas_cos_vardef,  only: kl_comm
#endif

! Subroutine arguments
  INTEGER (KIND=iintegers), INTENT(IN)   ::       &
    irealgrib          ! Kind parameter for the grib library

  INTEGER (KIND=iintegers), INTENT(OUT)  ::       &
    nproc,           & ! total number of processors: nprocx * nprocy + nprocio
    my_world_id,     & ! rank of this subdomain in the global communicator
    icomm_world,     & ! communicator that belongs to igroup_world
                       ! (is MPI_COMM_WORLD, if no coupling is used)
    igroup_world,    & ! group that belongs to MPI_COMM_WORLD, i.e. all processors
    imp_real,        & ! determines the correct REAL type used in the model for MPI
    imp_single,      & ! single precision REAL type for MPI
    imp_double,      & ! double precision REAL type for MPI
    imp_grib,        & ! determines the REAL type for the GRIB library
    imp_int,         & ! determines the correct INTEGER type used in the model for MPI
    imp_byte,        & ! determines the correct BYTE type used in the model for MPI
    imp_char,        & ! determines the correct CHARACTER type used in the model for MPI
    imp_logical,     & ! determines the correct LOGICAL   type used in the model for MPI
    imp_intga,       & ! determines the correct INTEGER type used for grib_api

    iexch_req(4),    & ! stores the sends requests for the neighbor-exchange
                       ! that can be used by MPI_WAIT to identify the send
    ierror             ! error status variable

  CHARACTER (LEN=100),      INTENT(OUT)  ::       &
    yerrmsg            ! for MPI error message

! Local variables:

INTEGER (KIND=iintegers)           :: izmplcode

!------------------------------------------------------------------------------

! Begin subroutine init_environment

  ierror     = 0
  izmplcode  = 0
  yerrmsg(:) = ' '

#if defined(MESSYMMD) && defined(COUP_OAS_COS)
  WRITE (*,'(a)') 'ERROR: MESSy-MMD (-DMESSYMMD) and OASIS (-DCOUP_OAS) cannot run at the same time!'
  ierror = 33
  yerrmsg = 'ERROR: MESSy-MMD (-DMESSYMMD) and OASIS (-DCOUP_OAS) cannot run at the same time!'
  RETURN
#endif

! If the coupler OASIS is active, the COSMO-Model gets its communicator
! from OASIS
#ifdef COUP_OAS_COS
  CALL oas_cos_init
  icomm_world = kl_comm
#else
  ! Initialize MPI
  CALL MPI_INIT (izmplcode)
  IF (izmplcode /= 0) THEN
    ierror  = izmplcode
    yerrmsg = 'MPI_INIT'
    RETURN
  ENDIF
#ifdef MESSYMMD
  call MMD_get_model_communicator (MMD_comm_world, izmplcode)
  if(izmplcode /= MMD_STATUS_OK)   then
    ierror  = izmplcode
    yerrmsg = 'MMD_model_communicator'
    RETURN

  end if
  icomm_world = MMD_comm_world ! mz_kk_20081107
#else
  icomm_world = MPI_COMM_WORLD
#endif
#endif

  ! Determination of number of present processors
  CALL MPI_COMM_SIZE (icomm_world, nproc, izmplcode)
  IF (izmplcode /= 0) THEN
    ierror  = izmplcode
    yerrmsg = 'MPI_COMM_SIZE'
    RETURN
  ENDIF

  ! Determination of own rank (= process_id) in the overall communicator 
  CALL MPI_COMM_RANK (icomm_world, my_world_id, izmplcode)
  IF (izmplcode /= 0) THEN
    ierror  = izmplcode
    yerrmsg = 'MPI_COMM_RANK'
    RETURN
  ENDIF

  ! Determination of the global group and the corresponding communicator
  CALL MPI_COMM_GROUP (icomm_world, igroup_world, izmplcode)
  IF (izmplcode /= 0) THEN
    ierror  = izmplcode
    yerrmsg = 'MPI_COMM_GROUP'
    RETURN
  ENDIF

  ! Determine the type of REALs and INTEGERs for MPI and other variables
  ! If the KIND-type parameters in data_parameters are changed, the
  ! variables here have to be changed accordingly.
  ! Model Real variables
  imp_single = MPI_REAL
  imp_double = MPI_DOUBLE_PRECISION
  IF     (KIND (1.0)   == wp) THEN
    imp_real = imp_single
  ELSEIF (KIND (1.0D0) == wp) THEN
    imp_real = imp_double
  ELSEIF (KIND (1.0) == 8 .AND. wp == 4) THEN
    ! it seems that this is a T3E where 4 Byte REALs are used
    imp_real = MPI_REAL4
  ELSE
    ierror = 1
    yerrmsg = 'cannot find MPI REAL Type for working precision wp'
    RETURN
  ENDIF
  ! GRIB Real variables
  IF     (KIND (1.0)   == irealgrib) THEN
    imp_grib = MPI_REAL
  ELSEIF (KIND (1.0D0) == irealgrib) THEN
    imp_grib = MPI_DOUBLE_PRECISION
  ELSEIF (KIND (1.0) == 8 .AND. irealgrib == 4) THEN
    ! it seems that this is a T3E where 4 Byte REALs are used
    imp_real = MPI_REAL4
  ELSE
    ierror = 1
    yerrmsg = 'cannot find MPI REAL Type for irealgrib'
    RETURN
  ENDIF

  IF     (KIND (1)   == iintegers) THEN
    imp_int     = MPI_INTEGER
  ELSE
    ierror = 2
    yerrmsg = 'cannot find MPI INTEGER Type for iintegers'
    RETURN
  ENDIF

! This could be used with MPI 3.0
! CALL MPI_SIZEOF (1_int_ga, izsize, izmplcode)
! CALL MPI_TYPE_MATCH_SIZE (MPI_TYPECLASS_INTEGER, izsize, imp_intga, izmplcode)

! and this is the old style
  IF     (HUGE(1_int_ga) == INT(HUGE(1_iintegers), int_ga)) THEN
    imp_intga   = MPI_INTEGER
  ELSEIF (HUGE(1_int_ga) >  INT(HUGE(1_iintegers), int_ga)) THEN
    imp_intga   = MPI_INTEGER8
  ELSE
    ierror = 2
    yerrmsg = 'cannot find MPI INTEGER Type for int_ga'
    RETURN
  ENDIF

  imp_byte    = MPI_BYTE
  imp_char    = MPI_CHARACTER
  imp_logical = MPI_LOGICAL

  ! Set all exchange requests to null 
  iexch_req(:) = MPI_REQUEST_NULL

END SUBROUTINE init_environment

!==============================================================================
!+ Subroutine that finalizes the parallel environment
!------------------------------------------------------------------------------

SUBROUTINE final_environment (ierror, yerrmsg)

!------------------------------------------------------------------------------
!
! Description:
!   This subroutine finalizes the environment of the parallel program.
!   No more MPI-routines can be called after this operation.
!
! Method:
!   MPI-FINALIZE
!
!------------------------------------------------------------------------------

! Subroutine arguments
  INTEGER (KIND=iintegers), INTENT(OUT)  ::       &
    ierror             ! error status variable

  CHARACTER (LEN=25),       INTENT(OUT)  ::       &
    yerrmsg            ! for MPI error message

! Local variables:
INTEGER (KIND=iintegers)           :: izmplcode, istat

!------------------------------------------------------------------------------

! Begin subroutine final_environment

  ierror     = 0
  yerrmsg    = '    '
  izmplcode  = 0
  istat      = 0

#ifdef COUP_OAS_COS
  CALL oas_cos_finalize
#else
  CALL MPI_FINALIZE (izmplcode)
  IF (izmplcode /= 0) THEN
    ierror  = izmplcode
    yerrmsg = 'MPI_FINALIZE'
  ENDIF
#endif

END SUBROUTINE final_environment

!==============================================================================
!+ This routine determines the next free unit number
!------------------------------------------------------------------------------

SUBROUTINE get_free_unit (iunit)

!------------------------------------------------------------------------------
!
! Description:
!   A table with the unit numbers from iunit_start to iunit_end (values see
!   above) is defined. This routine gets the next free unit number for the
!   calling routine. All unit numbers for the program shall be determined 
!   with this routine in order to avoid problems with the same unit number for
!   different files
!
! Method:
!   The entries in the table iunit_table are 0, if the unit is not used and
!   they are 1, if the unit is used already.
!
!------------------------------------------------------------------------------

! Subroutine arguments
  INTEGER (KIND=iintegers), INTENT(OUT)   ::       &
    iunit              ! is set to the next free unit number

  INTEGER (KIND=iintegers) :: i

! Begin subroutine get_free_unit

  iunit = -1
  search_loop: DO i = iunit_start, iunit_end
    IF (iunit_table(i) == 0) THEN
      iunit_table(i) = 1
      iunit          = i
      EXIT search_loop
    ENDIF
  ENDDO search_loop

  IF (iunit == -1) THEN
    PRINT *, ' *** WARNING:  No more free unit numbers available *** '
  ENDIF

END SUBROUTINE get_free_unit

!==============================================================================
!+ This routine defines a unit number as free again
!------------------------------------------------------------------------------

SUBROUTINE release_unit (iunit)

!------------------------------------------------------------------------------
!
! Description:
!   A table with the unit numbers from iunit_start to iunit_end (values see
!   above) is defined. This routine sets the number iunit as free again.
!
! Method:
!   The entry "iunit" in the table iunit_table is set to 0 again.
!
!------------------------------------------------------------------------------

! Subroutine arguments
  INTEGER (KIND=iintegers), INTENT(IN)    ::       &
    iunit              ! is set to the next free unit number

! Begin subroutine release_unit

  IF (iunit_table(iunit) == 1) THEN
    iunit_table(iunit) = 0
  ELSE
    PRINT *, ' *** WARNING:  The unit ', iunit, ' is already free *** '
  ENDIF

END SUBROUTINE release_unit

!==============================================================================
!+ Subroutine for stopping the program in case of errors
!------------------------------------------------------------------------------

SUBROUTINE model_abort (my_id,ierrorcode,yerrorstring,yroutine, implerrorcode)

!------------------------------------------------------------------------------
!
! Description:
!   This subroutine stops the program in case of errors. Model_abort prints 
!   a self-defined error message, the errorcode and the calling routine.
!   If the error occurs in the routine of the message passing library used, 
!   the optional parameter implerrorcode must be present. In this case
!   also an appropriate error message is printed. This is determined via 
!   MPI_ERROR_STRING and implerrorcode.
!   If the error occurs in several nodes, this message will also be printed
!   several times.
!
! Method:
!   Printing an error message and then MPI_ABORT.
!
!------------------------------------------------------------------------------

! Parameter list:

#ifdef COUP_OAS_COS
USE data_parallel, ONLY : icomm_world
#endif

INTEGER (KIND=iintegers),INTENT(IN)   ::                                     &
  my_id,        & ! id of this processor
  ierrorcode      ! self-defined integer code of the error detected

CHARACTER (LEN=*), INTENT(IN) ::     &
   yerrorstring   ! self-defined error message
CHARACTER (LEN=*), INTENT(IN) ::     &
   yroutine       ! calling routine

INTEGER (KIND=iintegers), OPTIONAL, INTENT(IN)   ::                          &
  implerrorcode   ! error-code of the message passing library

! Local variables:
INTEGER (KIND=iintegers)   :: i, nzerrcode, nzlen
LOGICAL                    :: lzopen
CHARACTER (LEN=100) ymplmsg   ! for MPI error message

!------------------------------------------------------------------------------

  IF (my_id == 0) THEN

    ! print the error message
    IF (PRESENT (implerrorcode)) THEN
      ! this is parallel mode
      CALL MPI_ERROR_STRING (implerrorcode, ymplmsg, nzlen, nzerrcode)

      PRINT *,'*------------------------------------------------------------*'
      PRINT *,'*    PROGRAM TERMINATED BECAUSE OF MPI ERRORS DETECTED'
      PRINT *,'*              IN ROUTINE:   ',yroutine
      PRINT *,'*'
      PRINT *,'*    ERROR CODE is ',ierrorcode,': '
      PRINT *,'*    MPI ROUTINE:  ',yerrorstring(1:LEN_TRIM(yerrorstring))
      PRINT *,'*    MPI ERROR CODE is ', implerrorcode,':  ',ymplmsg
      PRINT *,'*------------------------------------------------------------*'
    ELSE
      PRINT *,'*------------------------------------------------------------*'
      PRINT *,'*    PROGRAM TERMINATED BECAUSE OF ERRORS DETECTED'
      PRINT *,'*              IN ROUTINE:   ',yroutine
      PRINT *,'*'
      PRINT *,'*    ERROR CODE is ',ierrorcode
      PRINT *,'*    ', yerrorstring(1:LEN_TRIM(yerrorstring))
      PRINT *,'*------------------------------------------------------------*'
    ENDIF

    ! Check, whether there are open files and close them
    DO i = iunit_start, iunit_end
      IF (iunit_table(i) == 1) THEN
        ! inquire, whether unit i is open
        lzopen = .FALSE.
        INQUIRE (UNIT=i, OPENED=lzopen)
        IF (lzopen) THEN
          CLOSE (i)
        ENDIF
      ENDIF
    ENDDO

  ELSE

    ! all other processors first sleep for a while, to give my_id == 0 the 
    ! opportunity to print the error message and close all files
    ! If my_id == 0 does not call model_abort, then one of the other PEs
    ! will print the message and abort

    i = dosleep(30)

!US#ifdef GRIBDWD
!US    CALL fsleep (30)
!US#endif

    IF (PRESENT (implerrorcode)) THEN
      ! this is parallel mode
      CALL MPI_ERROR_STRING (implerrorcode, ymplmsg, nzlen, nzerrcode)

      PRINT *,'*------------------------------------------------------------*'
      PRINT *,'*    PROGRAM TERMINATED BECAUSE OF MPI ERRORS DETECTED'
      PRINT *,'*              IN ROUTINE:   ',yroutine
      PRINT *,'*'
      PRINT *,'*    ERROR CODE is ',ierrorcode,': '
      PRINT *,'*    MPI ROUTINE:  ',yerrorstring(1:LEN_TRIM(yerrorstring))
      PRINT *,'*    MPI ERROR CODE is ', implerrorcode,':  ',ymplmsg
      PRINT *,'*------------------------------------------------------------*'
    ELSE
      PRINT *,'*------------------------------------------------------------*'
      PRINT *,'*    PROGRAM TERMINATED BECAUSE OF ERRORS DETECTED'
      PRINT *,'*              IN ROUTINE:   ',yroutine
      PRINT *,'*'
      PRINT *,'*    ERROR CODE is ',ierrorcode
      PRINT *,'*    ', yerrorstring(1:LEN_TRIM(yerrorstring))
      PRINT *,'*------------------------------------------------------------*'
    ENDIF

  ENDIF

!US
! to avoid dependency on icomm_world and also to avoid adding icomm_world
! to the argument list (which would change hundreds of calls to model_abort)
! we implemented the call to MPI_ABORT also with the ifdef COUP_OAS

! Question: wouldn't it be useful, also in a coupled environment, to abort the
!           whole application, if the model fails?

#ifdef COUP_OAS_COS
  CALL MPI_ABORT (icomm_world, ierrorcode, nzerrcode)
#else
  CALL MPI_ABORT (MPI_COMM_WORLD, ierrorcode, nzerrcode)
#endif

END SUBROUTINE model_abort

!==============================================================================
!+ Subroutine that creates the virtual processor topology.
!------------------------------------------------------------------------------

SUBROUTINE init_procgrid                                                        &
     (nproc, nprocx, nprocy, nprocio, nc_asyn_io, lperi_x, lperi_y, lreproduce, &
      lreorder, icomm_world, my_world_id, icomm_compute,                        &
      icomm_asynio, icomm_cart, igroup_cart, my_cart_id, my_cart_pos,           &
      my_cart_neigh, icomm_row, lcompute_pe, yerrmsg, ierror )

!------------------------------------------------------------------------------
!
! Description:
!   This subroutine creates the communicator for the compute PEs and the 
!   virtual topology of the twodimensional cartesian processor grid. Both
!   communicators are related to the same group (the compute PEs) but in a 
!   different context (unordered vs. cartesian grid).
!   The following organizational variables are determined using MPI-calls.
!    - icomm_compute:    MPI-communicator the (unordered) group of compute PEs
!    - lcompute_pe:      .TRUE., if this PE belongs to icomm_compute
!    - icomm_cart:       MPI-communicator for the cartesian grid
!    - igroup_cart:      MPI-group for this communicator icomm_cart
!    - my_cart_id:       rank and id of this processor in the virtual topology
!    - my_cart_pos(2):   position in the cartesian processor grid in 
!                           x- and y-direction
!    - my_cart_neigh(4): neighbours of this processor in the order west, north,
!                        east, south
!    - icomm_row:        MPI-communicator for east-west processor rows of the
!                        cartesian grid in case of desired reproducibility.
!
!   For the sequential program the variables are set accordingly.
!
! Method:
!   MPI-calls and computations.
!
!------------------------------------------------------------------------------

! Subroutine arguments
  INTEGER (KIND=iintegers), INTENT(IN)   ::       &
    icomm_world,       & ! communicator for all processors
    my_world_id,       & ! ID of this PE in icomm_world
    nproc,             & ! number of total processes
    nprocx,            & ! number of processes in east-west-direction
    nprocy,            & ! number of processes in north-south-direction
    nprocio,           & ! number of extra processes for IO (GRIB)
    nc_asyn_io           ! number of extra processes for IO (NetCDF)

  LOGICAL,                  INTENT(IN)   ::       &
    lperi_x,           & ! special treatment for periodic boundary conditions in x-dir.
    lperi_y,           & ! special treatment for periodic boundary conditions in y-dir.
    lreproduce,        & ! additional communication necessary for creating
                         ! reproducible results
    lreorder             ! during the creation of the virtual topology the
                         ! ranking of the processors may be reordered

  INTEGER (KIND=iintegers), INTENT(OUT)  ::       &
    icomm_compute,     & ! communicator for compute PEs
    icomm_asynio,      & ! communicator for netcdf asynchronous PEs
    icomm_cart,        & ! communicator for the virtual cartesian topology
    igroup_cart,       & ! MPI-group for icomm_cart
    icomm_row,         & ! communicator for the group of a east-west processor
                         ! row (needed in the LM for reproducible results)
    my_cart_id,        & ! rank of this subdomain in the cartesian communicator
    my_cart_pos(2),    & ! position of this subdomain in the cartesian grid
                         ! in x- and y-direction
    my_cart_neigh(4),  & ! neighbors of this subdomain in the cartesian grid
    ierror               ! error status variable

  LOGICAL,                  INTENT(OUT)  ::       &
    lcompute_pe          ! TRUE, if this PE belongs to icomm_compute

  CHARACTER (LEN=50),       INTENT(OUT)  ::       &
    yerrmsg              ! for MPI error message

!------------------------------------------------------------------------------

! Local variables
  INTEGER (KIND=iintegers)   ::       &
    nznumdims,         & ! number of dimensions of the virtual cartesian 
                         ! topology of the processor grid
    nzsizdims(2),      & ! number of processes in each dimension 
    izranks(0:nproc-1),& ! for building groups and communicators
    izgroup,           & ! group for processors rows
    izcomm,            & ! communicator for processors rows
    izmplcode,         & ! local error status variable for MPI-calls
    nzneigh, ij(2), ix, iy, nziope0, color_comp   ! dummy variables

  LOGICAL                    ::       &
    lzperiods(2)     ! specifies whether the processor grid should be periodic
                     ! (.TRUE.) in each dimension

!------------------------------------------------------------------------------
!- End of header -
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!- Section 1: Initializations
!------------------------------------------------------------------------------

  ierror     = 0
  izmplcode  = 0 
  yerrmsg    = '    '
 
  ! Set defaults as if this is a sequential program
  icomm_compute    = icomm_world
  icomm_asynio     = icomm_world
  lcompute_pe      = .TRUE.
  my_cart_id       = 0
  icomm_cart       = icomm_world
  icomm_row        = icomm_world
  my_cart_pos(:)   = 0
  my_cart_neigh(:) = -1

  parallel: IF (nproc > 1) THEN

!------------------------------------------------------------------------------
!- Section 2: Create the communicator for the compute PEs
!------------------------------------------------------------------------------

    nziope0 = 0

    IF (nprocio > 0) THEN
      nziope0 = nproc - nprocio ! first PE of IO-group
    ELSE IF (nc_asyn_io > 0) THEN
      nziope0 = nproc - nc_asyn_io ! first PE of IO-group
    ENDIF

    ! nziope0 has to be greater than 0, restricted by namelists
    IF( nziope0 /= 0 ) THEN

      IF( my_world_id >= nziope0) THEN     ! last PEs = I/O PE
        lcompute_pe   = .FALSE.
        color_comp    =  MPI_UNDEFINED
      ELSE
        lcompute_pe   = .TRUE.
        color_comp    =  0
      ENDIF

      ! PE 0,..., nproc-nprocio-1 are compute PEs
      CALL MPI_COMM_SPLIT( icomm_world, color_comp, 0, icomm_compute, &
                           izmplcode)
      IF (izmplcode /= 0) THEN
        ierror  = izmplcode
        yerrmsg = 'MPI_COMM_SPLIT'
        RETURN
      ENDIF

      ! in the else-case the initial settings above are valid
    ENDIF

!------------------------------------------------------------------------------
!- Section 3: Create the virtual topology of the processor grid
!             Determine the exact shape of the processor grid
!------------------------------------------------------------------------------

    compute: IF (lcompute_pe .AND. (nprocx*nprocy > 1)) THEN

      !------------------------------------------------------------------------
      !- Section 3.1: Determine the exact shape of the processor grid
      !------------------------------------------------------------------------

      ! The following initializations determine the exact shape of the virtual
      ! topology for the cartesian processor grid. The dimensions and the 
      ! size of each dimension is specified.
      ! The creation of the grid may reorder the ranks (i.e. id's) of the 
      ! processors, if lreorder = .TRUE.. Therefore the rank 
      ! has to be specified again afterwards
      nznumdims    = 2
      nzsizdims(1) = nprocx
      nzsizdims(2) = nprocy
      lzperiods(1) = lperi_x
      lzperiods(2) = lperi_y

      ! The following call creates the virtual processor topology for the 
      ! cartesian processor grid according to the specifications in above

      CALL MPI_CART_CREATE (    &
               icomm_compute,   & ! old communicator
               nznumdims,       & ! number of dimensions
               nzsizdims,       & ! number of processes in each dimension
               lzperiods,       & ! periodic grid or not
               lreorder,        & ! reordering of the ranks
               icomm_cart,      & ! new communicator with cartesian topology
               izmplcode)         ! error status variable
 
      IF (izmplcode /= 0) THEN
        ierror  = izmplcode
        yerrmsg = 'MPI_CART_CREATE'
        RETURN
      ENDIF

      CALL MPI_COMM_GROUP (icomm_cart, igroup_cart, izmplcode)
      IF (izmplcode /= 0) THEN
        ierror  = izmplcode
        yerrmsg = 'MPI_COMM_GROUP'
        RETURN
      ENDIF

      !------------------------------------------------------------------------
      !- Section 3.2: Get grid characteristics (rank, coordinates)
      !------------------------------------------------------------------------

      ! Get the rank of this processor in the new cartesian communicator 
      CALL MPI_COMM_RANK (icomm_cart, my_cart_id, izmplcode)
      IF (izmplcode /= 0) THEN
        ierror  = izmplcode
        yerrmsg = 'MPI_COMM_RANK'
        RETURN
      ENDIF

      ! Get the position of this processor in the 2D cartesian gid
      CALL MPI_CART_COORDS (icomm_cart, my_cart_id , nznumdims,           &
                                        my_cart_pos, izmplcode)
      IF (izmplcode /= 0) THEN
        ierror  = izmplcode
        yerrmsg = 'MPI_CART_COORDS'
        RETURN
      ENDIF

      !------------------------------------------------------------------------
      !- Section 3.3: Determine the neighbours of this processor
      !------------------------------------------------------------------------

      IF (my_cart_pos(1) == 0) THEN
        ! no neighbour to the west (left)
        my_cart_neigh(1) = -1
      ELSE
        ij(1) = my_cart_pos(1) - 1
        ij(2) = my_cart_pos(2)
        CALL MPI_CART_RANK (icomm_cart, ij, nzneigh, izmplcode)

        IF (izmplcode /= 0) THEN
          ierror  = izmplcode
          yerrmsg = 'MPI_CART_RANK'
          RETURN
        ENDIF
        my_cart_neigh(1) = nzneigh
      ENDIF

      IF (my_cart_pos(2) == nprocy-1) THEN
        ! no neighbour to the north (up)
        my_cart_neigh(2) = -1
      ELSE
        ij(1) = my_cart_pos(1)
        ij(2) = my_cart_pos(2)+1
        CALL MPI_CART_RANK (icomm_cart, ij, nzneigh, izmplcode)
  
        IF (izmplcode /= 0) THEN
          ierror  = izmplcode
          yerrmsg = 'MPI_CART_RANK'
          RETURN
        ENDIF
  
        my_cart_neigh(2) = nzneigh
      ENDIF
  
      IF (my_cart_pos(1) == nprocx-1) THEN
        ! no neighbour to the east (right)
        my_cart_neigh(3) = -1
      ELSE
        ij(1) = my_cart_pos(1)+1
        ij(2) = my_cart_pos(2)
        CALL MPI_CART_RANK (icomm_cart, ij, nzneigh, izmplcode)
  
        IF (izmplcode /= 0) THEN
          ierror  = izmplcode
          yerrmsg = 'MPI_CART_RANK'
          RETURN
        ENDIF
  
        my_cart_neigh(3) = nzneigh
      ENDIF
  
      IF (my_cart_pos(2) == 0) THEN
        ! no neighbour to the south (down)
        my_cart_neigh(4) = -1
      ELSE
        ij(1) = my_cart_pos(1)
        ij(2) = my_cart_pos(2)-1
        CALL MPI_CART_RANK (icomm_cart, ij, nzneigh, izmplcode)
  
        IF (izmplcode /= 0) THEN
          ierror  = izmplcode
          yerrmsg = 'MPI_CART_RANK'
          RETURN
        ENDIF
  
        my_cart_neigh(4) = nzneigh
      ENDIF

      !------------------------------------------------------------------------
      !- Section 3.4: Determine the neighbours in case of periodic boundaries
      !------------------------------------------------------------------------

      IF ( lperi_x ) THEN
  
        IF (my_cart_pos(1) == 0) THEN
          ! periodic neighbour to the west (left)
          ij(1) = nprocx-1
          ij(2) = my_cart_pos(2)
          CALL MPI_CART_RANK (icomm_cart, ij, nzneigh, izmplcode)

          IF (izmplcode /= 0) THEN
            ierror  = izmplcode
            yerrmsg = 'MPI_CART_RANK'
            RETURN
          ENDIF

          my_cart_neigh(1) = nzneigh
        ENDIF
  
        IF (my_cart_pos(1) == nprocx-1) THEN
          ! periodic neighbour to the east (right)
          ij(1) = 0
          ij(2) = my_cart_pos(2)
          CALL MPI_CART_RANK (icomm_cart, ij, nzneigh, izmplcode)

          IF (izmplcode /= 0) THEN
            ierror  = izmplcode
            yerrmsg = 'MPI_CART_RANK'
            RETURN
          ENDIF

          my_cart_neigh(3) = nzneigh
        ENDIF

      ENDIF   ! periodic boundaries in x-dir.

      IF ( lperi_y) THEN

        IF (my_cart_pos(2) == nprocy-1) THEN
          ! periodic neighbour to the north (up)
          ij(1) = my_cart_pos(1)
          ij(2) = 0
          CALL MPI_CART_RANK (icomm_cart, ij, nzneigh, izmplcode)

          IF (izmplcode /= 0) THEN
            ierror  = izmplcode
            yerrmsg = 'MPI_CART_RANK'
            RETURN
          ENDIF

          my_cart_neigh(2) = nzneigh
        END IF

        IF (my_cart_pos(2) == 0) THEN
          ! periodic neighbour to the south (down)
          ij(1) = my_cart_pos(1)
          ij(2) = nprocy-1
          CALL MPI_CART_RANK (icomm_cart, ij, nzneigh, izmplcode)

          IF (izmplcode /= 0) THEN
            ierror  = izmplcode
            yerrmsg = 'MPI_CART_RANK'
            RETURN
          ENDIF

          my_cart_neigh(4) = nzneigh
        ENDIF

      ENDIF   ! periodic boundaries in y-dir.

      !------------------------------------------------------------------------
      !- Section 3.5: Get communicator for east-west processor row
      !------------------------------------------------------------------------
 
      IF (lreproduce) THEN
        DO iy = 0, nprocy-1
          DO ix = 0, nprocx-1
            ! determine ranks of the group iy
            izranks(ix) = iy + ix * nprocy
          ENDDO
    
          ! Determine group and communicator for row iy
          CALL MPI_GROUP_INCL                                                   &
                  ( igroup_cart, nprocx, izranks(0), izgroup, izmplcode)
          IF (izmplcode /= 0) THEN
            ierror  = izmplcode
            yerrmsg = 'MPI_GROUP_INCL'
            RETURN
          ENDIF
    
          CALL MPI_BARRIER ( icomm_cart, izmplcode)
    
          CALL MPI_COMM_CREATE                                                  &
                  ( icomm_cart, izgroup, izcomm, izmplcode)
          IF (izmplcode /= 0) THEN
            ierror  = izmplcode
            yerrmsg = 'MPI_COMM_CREATE'
            RETURN
          ENDIF
  
          ! Store the own communicator
          IF (my_cart_pos(2) == iy) THEN
            icomm_row = izcomm
          ENDIF
        ENDDO
      ENDIF
    ENDIF  compute

!------------------------------------------------------------------------------
!- Section 4: Set some of these values for IO-PEs
!------------------------------------------------------------------------------

    IF (.NOT. lcompute_pe) THEN
      my_cart_id = MPI_UNDEFINED
      icomm_cart = MPI_UNDEFINED
      icomm_row  = MPI_UNDEFINED
    ENDIF

  ELSE   ! serial run

    IF (lperi_x) THEN
      my_cart_neigh(1) = 0
      my_cart_neigh(3) = 0
    END IF

    IF (lperi_y) THEN
      my_cart_neigh(2) = 0
      my_cart_neigh(4) = 0
    END IF

  ENDIF  parallel

!------------------------------------------------------------------------------
!- End of the Subroutine
!------------------------------------------------------------------------------

END SUBROUTINE init_procgrid

!==============================================================================
!==============================================================================
!+ This subroutine performs the data exchange between boundaries
!------------------------------------------------------------------------------


!!! No extra treatment of staggered u/v gridpoints in periodic exchange
!!! necessary!!!

SUBROUTINE exchg_boundaries                                                  &
               ( icase, sendbuf, isendbuflen, imp_type, icomm, num_compute,  &
                 idim, jdim, kdim, jstartpar, jendpar, nlines, nboundlines,  &
                 neighbors, lperi_x, lperi_y, l2dim, ntag, lmpi_types, ntype,&
                 ierror, yerrmsg,                                            &
                 var01, var02, var03, var04, var05, var06, var07, var08,     &
                 var09, var10, var11, var12, var13, var14, var15, var16,     &
                 var17, var18, var19, var20, var21, var22, var23, var24 )

!------------------------------------------------------------------------------
!
! Description:
!   This subroutine performs the boundary exchange of up to 20 variables. Only
!   one variable has to occur, the others are optional. 
!
! Method:
!   At the moment there are 3 different MPI-communications implemented:
!     1) immediate send, blocking receive and wait
!     2) immediate receive, blocking send and wait
!     3) MPI_SendRecv
!   Also there is the choice of an explicit buffering (putbuf, getbuf) or 
!   implicit buffering (MPI-Datatypes) of the data to be send.
!
!------------------------------------------------------------------------------

! Subroutine arguments
  INTEGER (KIND=iintegers), INTENT (IN)         ::    &
    icase,              & ! tag for exchange scenario
    isendbuflen,        & ! length of sendbuffer
    imp_type,           & ! determines the REAL type used
    icomm,              & ! communicator for virtual cartesian topology
    idim, jdim,         & ! horizontal dimensions of the fields
    kdim(24),           & ! array for the vertical dimensions of var1..var20
    jstartpar,          & ! start index in j-direction
    jendpar,            & ! end index in j-direction
    nlines,             & ! number of lines that have to be exchanged
                          ! (<= nboundlines)
    nboundlines,        & ! number of overlapping boundary lines
    neighbors(4),       & ! process-id's of the neighbors in the grid
    ntag,               & ! tag of the message
    ntype,              & ! indicates how the communication should be done
    num_compute           ! number of compute processors

  LOGICAL, INTENT(IN)                           ::    &
    lmpi_types,         & ! whether implicit (with MPI-Datatypes) or explicit
                          ! (putbuf, getbuf) buffering of data is used
    lperi_x, lperi_y, l2dim

  INTEGER (KIND=iintegers), INTENT (OUT)        ::    &
    ierror                ! error status variable

  CHARACTER (LEN=*),        INTENT(OUT)  ::       &
    yerrmsg               ! for MPI error message

  REAL (KIND=wp),           INTENT (INOUT)      ::    &
    sendbuf (isendbuflen, 8)   ! send buffer

  REAL (KIND=wp),           INTENT (INOUT), TARGET      ::    &
    var01 (idim, jdim, kdim( 1))   ! first field that has to occur

  REAL (KIND=wp),     OPTIONAL, INTENT (INOUT), TARGET  ::    &
    var02 (idim, jdim, kdim( 2)),& ! additional optional fields
    var03 (idim, jdim, kdim( 3)),& ! additional optional fields
    var04 (idim, jdim, kdim( 4)),& ! additional optional fields
    var05 (idim, jdim, kdim( 5)),& ! additional optional fields
    var06 (idim, jdim, kdim( 6)),& ! additional optional fields
    var07 (idim, jdim, kdim( 7)),& ! additional optional fields
    var08 (idim, jdim, kdim( 8)),& ! additional optional fields
    var09 (idim, jdim, kdim( 9)),& ! additional optional fields
    var10 (idim, jdim, kdim(10)),& ! additional optional fields
    var11 (idim, jdim, kdim(11)),& ! additional optional fields
    var12 (idim, jdim, kdim(12)),& ! additional optional fields
    var13 (idim, jdim, kdim(13)),& ! additional optional fields
    var14 (idim, jdim, kdim(14)),& ! additional optional fields
    var15 (idim, jdim, kdim(15)),& ! additional optional fields
    var16 (idim, jdim, kdim(16)),& ! additional optional fields
    var17 (idim, jdim, kdim(17)),& ! additional optional fields
    var18 (idim, jdim, kdim(18)),& ! additional optional fields
    var19 (idim, jdim, kdim(19)),& ! additional optional fields
    var20 (idim, jdim, kdim(20)),& ! additional optional fields
    var21 (idim, jdim, kdim(21)),& ! additional optional fields
    var22 (idim, jdim, kdim(22)),& ! additional optional fields
    var23 (idim, jdim, kdim(23)),& ! additional optional fields
    var24 (idim, jdim, kdim(24))   ! additional optional fields
! Local variables

  INTEGER (KIND=iintegers)   ::       &
    ! the following numbers are for filling/emptying the buffers for each
    ! neighbor
    izlo_lr, izup_lr, jzlo_lr, jzup_lr,     & ! left , receive
    izlo_rr, izup_rr, jzlo_rr, jzup_rr,     & ! right, receive
    izlo_ur, izup_ur, jzlo_ur, jzup_ur,     & ! upper, receive
    izlo_dr, izup_dr, jzlo_dr, jzup_dr,     & ! down , receive
    izlo_ls, izup_ls, jzlo_ls, jzup_ls,     & ! left , send
    izlo_rs, izup_rs, jzlo_rs, jzup_rs,     & ! right, send
    izlo_us, izup_us, jzlo_us, jzup_us,     & ! upper, send
    izlo_ds, izup_ds, jzlo_ds, jzup_ds,     & ! down , send
    nzcount_ls, nzcount_rs,     & ! counting the values
    nzcount_us, nzcount_ds,     & ! counting the values
    nzcount_lr, nzcount_rr,     & ! counting the values
    nzcount_ur, nzcount_dr,     & ! counting the values
    nzrequest(MPI_STATUS_SIZE), & ! for MPI-receive
    nzstatus (MPI_STATUS_SIZE), & ! for MPI-WAIT
    ncount, type_handle,        & ! return values from setup_data_type
    MPI_neighbors(4), i, j, k,  & ! same as neighbors, if neighbor exists
    ilocalreq(4),               & ! the local requests for the ISEND and IRECV
    nlines_x, nlines_y            ! to adapt nlines=nboundlines in case of parallel periodic exchange

  INTEGER (KIND=iintegers)   ::       &
    izmplcode                   ! for MPI error code

  LOGICAL :: zlvarpresent(24), ldebugflag
  CHARACTER  (LEN=200)          :: yzerrmsg

  TYPE :: pointerto3d
    REAL(KIND=wp),     POINTER, DIMENSION(:,:,:) :: p
  END TYPE pointerto3d

  TYPE(pointerto3d) :: varxxp(24)


!------------------------------------------------------------------------------
!- End of header -
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!- Section 1: Initializations
!------------------------------------------------------------------------------

  DO i=1,24
    NULLIFY(varxxp(i)%p)
  END DO

  !.. Enable/disable debugging mode:
  ldebugflag = .TRUE.

  ierror     = 0
  izmplcode  = 0
  yerrmsg    = '    '

  ! check whether nlines <= nboundlines
  IF (nlines > nboundlines) THEN
    ierror  = 9011
    yerrmsg     = ' *** nlines > nboundlines *** '
    CALL model_abort (0, ierror, TRIM(yerrmsg), 'EXCHG_BOUNDARIES()')
  ENDIF

  IF (lperi_x) THEN
    nlines_x = nboundlines
  ELSE
    nlines_x = nlines
  END IF

  IF (lperi_y) THEN
    nlines_y = nboundlines
  ELSE
    nlines_y = nlines
  END IF

  ! Determine the start- and end-indices for routines putbuf, getbuf
  izlo_ls = nboundlines + 1
  izup_ls = nboundlines + nlines_x
  jzlo_ls = jstartpar
  jzup_ls = jendpar

  izlo_lr = nboundlines + 1 - nlines_x
  izup_lr = nboundlines
  jzlo_lr = jstartpar
  jzup_lr = jendpar

  izlo_us = 1
  izup_us = idim
  jzlo_us = jdim - nboundlines - nlines_y + 1
  jzup_us = jdim - nboundlines

  izlo_ur = 1
  izup_ur = idim
  jzlo_ur = jdim - nboundlines + 1
  jzup_ur = jdim - nboundlines + nlines_y

  izlo_rs = idim - nboundlines - nlines_x + 1
  izup_rs = idim - nboundlines
  jzlo_rs = jstartpar
  jzup_rs = jendpar

  izlo_rr = idim - nboundlines + 1
  izup_rr = idim - nboundlines + nlines_x
  jzlo_rr = jstartpar
  jzup_rr = jendpar

  izlo_ds = 1
  izup_ds = idim
  jzlo_ds = nboundlines + 1
  jzup_ds = nboundlines + nlines_y

  izlo_dr = 1
  izup_dr = idim
  jzlo_dr = nboundlines + 1 - nlines_y
  jzup_dr = nboundlines

  nzcount_lr = 0
  nzcount_rr = 0
  nzcount_ur = 0
  nzcount_dr = 0
  nzcount_ls = 0
  nzcount_rs = 0
  nzcount_us = 0
  nzcount_ds = 0

  !.. Check whether optional exchange-variables are present but are 
  !   erroneously assiged a Z-dimension of 0:
  !
  !.. 1) Assign pointers to the input fields if they are present.

  zlvarpresent( 1) = .TRUE.
  zlvarpresent( 2) = PRESENT(var02)
  zlvarpresent( 3) = PRESENT(var03)
  zlvarpresent( 4) = PRESENT(var04)
  zlvarpresent( 5) = PRESENT(var05)
  zlvarpresent( 6) = PRESENT(var06)
  zlvarpresent( 7) = PRESENT(var07)
  zlvarpresent( 8) = PRESENT(var08)
  zlvarpresent( 9) = PRESENT(var09)
  zlvarpresent(10) = PRESENT(var10)
  zlvarpresent(11) = PRESENT(var11)
  zlvarpresent(12) = PRESENT(var12)
  zlvarpresent(13) = PRESENT(var13)
  zlvarpresent(14) = PRESENT(var14)
  zlvarpresent(15) = PRESENT(var15)
  zlvarpresent(16) = PRESENT(var16)
  zlvarpresent(17) = PRESENT(var17)
  zlvarpresent(18) = PRESENT(var18)
  zlvarpresent(19) = PRESENT(var19)
  zlvarpresent(20) = PRESENT(var20)
  zlvarpresent(21) = PRESENT(var21)
  zlvarpresent(22) = PRESENT(var22)
  zlvarpresent(23) = PRESENT(var23)
  zlvarpresent(24) = PRESENT(var24)

  !..  2) check kdim against the present variables:
  
  IF (ldebugflag) THEN
    DO k=1,24
      IF (zlvarpresent(k) .AND. kdim(k) <= 0) THEN
        yzerrmsg = REPEAT(' ', LEN(yzerrmsg))
        WRITE(yzerrmsg,'(a,i2.2,a,i2,a)') ' *** Error: var',k,&
             ' is exchanged but kdim(',k,') is 0! This may lead to erroneous exchange! *** '
        ierror = 54321
        CALL model_abort (0, ierror, TRIM(yzerrmsg), 'EXCHG_BOUNDARIES()')
      END IF
    END DO
  END IF

  !.. Assign pointers to the present exchange fields:

  IF (zlvarpresent( 1)) varxxp( 1)%p => var01
  IF (zlvarpresent( 2)) varxxp( 2)%p => var02
  IF (zlvarpresent( 3)) varxxp( 3)%p => var03
  IF (zlvarpresent( 4)) varxxp( 4)%p => var04
  IF (zlvarpresent( 5)) varxxp( 5)%p => var05
  IF (zlvarpresent( 6)) varxxp( 6)%p => var06
  IF (zlvarpresent( 7)) varxxp( 7)%p => var07
  IF (zlvarpresent( 8)) varxxp( 8)%p => var08
  IF (zlvarpresent( 9)) varxxp( 9)%p => var09
  IF (zlvarpresent(10)) varxxp(10)%p => var10
  IF (zlvarpresent(11)) varxxp(11)%p => var11
  IF (zlvarpresent(12)) varxxp(12)%p => var12
  IF (zlvarpresent(13)) varxxp(13)%p => var13
  IF (zlvarpresent(14)) varxxp(14)%p => var14
  IF (zlvarpresent(15)) varxxp(15)%p => var15
  IF (zlvarpresent(16)) varxxp(16)%p => var16
  IF (zlvarpresent(17)) varxxp(17)%p => var17
  IF (zlvarpresent(18)) varxxp(18)%p => var18
  IF (zlvarpresent(19)) varxxp(19)%p => var19
  IF (zlvarpresent(20)) varxxp(20)%p => var20
  IF (zlvarpresent(21)) varxxp(21)%p => var21
  IF (zlvarpresent(22)) varxxp(22)%p => var22
  IF (zlvarpresent(23)) varxxp(23)%p => var23
  IF (zlvarpresent(24)) varxxp(24)%p => var24


  IF (l2dim .OR. (num_compute == 1 .AND. lperi_y)) THEN
    ! Do serial exchange(s): periodic and/or 2-dim; presently all nboundlines are
    ! exchanged, not only nlines:

    ! Exchange needed in y-direction:
    DO k = 1, 24
      IF (zlvarpresent(k)) THEN
        DO j=1, nboundlines
          varxxp(k)%p(:,nboundlines+1-j   ,:) = varxxp(k)%p(:,jdim-nboundlines+1-j,:)
          varxxp(k)%p(:,jdim-nboundlines+j,:) = varxxp(k)%p(:,nboundlines+j       ,:)
        ENDDO
      ENDIF
    ENDDO

  ENDIF

  IF (num_compute == 1) THEN

    IF (lperi_x) THEN
      ! periodic BCs needed in x-direction:
      DO k = 1, 24
        IF (zlvarpresent(k)) THEN
          DO i=1,nboundlines
            varxxp(k)%p(nboundlines+1-i   ,:,:) = varxxp(k)%p(idim-nboundlines+1-i,:,:)
            varxxp(k)%p(idim-nboundlines+i,:,:) = varxxp(k)%p(nboundlines+i       ,:,:)
          END DO
        ENDIF
      ENDDO

    ENDIF

  ELSE  ! num_compute > 1

    ! Do parallel exchange, inner boundaries and periodic exchange
    ! together:

    ! Fix list of neighbors (use MPI_PROC_NULL rather than -1 to indicate
    ! missing neighbor).
    DO i= 1, 4
      IF ( neighbors(i) /= -1 ) THEN
        MPI_neighbors(i) = neighbors(i)
      ELSE
        MPI_neighbors(i) = MPI_PROC_NULL
      ENDIF
    ENDDO

  !------------------------------------------------------------------------------
  !- Section 2: Determine the necessary datatypes
  !------------------------------------------------------------------------------
  
    IF (lmpi_types) THEN
      ! Exchange with left and right neighbor
      IF ( iexchg_MPI_types(2*icase-1) == MPI_DATATYPE_NULL ) THEN
        IF ( MPI_neighbors(1) /= MPI_PROC_NULL ) THEN
          CALL setup_data_type                                               &
           ( var01,var02,var03,var04,var05,var06,var07,var08,var09,var10,    &
             var11,var12,var13,var14,var15,var16,var17,var18,var19,var20,    &
             var21,var22,var23,var24,                                        &
             idim, jdim, kdim, izlo_ls, izup_ls, jzlo_ls, jzup_ls,           &
             ierror, yerrmsg, imp_type, ncount, type_handle )
          iexchg_MPI_types(2*icase-1) = type_handle
          iexchg_counts   (2*icase-1) = ncount
        ELSE
          CALL setup_data_type                                               &
           ( var01,var02,var03,var04,var05,var06,var07,var08,var09,var10,    &
             var11,var12,var13,var14,var15,var16,var17,var18,var19,var20,    &
             var21,var22,var23,var24,                                        &
             idim, jdim, kdim, izlo_rr, izup_rr, jzlo_rr, jzup_rr,           &
             ierror, yerrmsg, imp_type, ncount, type_handle )
          iexchg_MPI_types(2*icase-1) = type_handle
          iexchg_counts   (2*icase-1) = ncount
        ENDIF
      ENDIF
    
      ! Exchange with upper and lower neighbor
      IF ( iexchg_MPI_types(2*icase) == MPI_DATATYPE_NULL ) THEN
        IF ( MPI_neighbors(2) /= MPI_PROC_NULL ) THEN
          CALL setup_data_type                                               &
           ( var01,var02,var03,var04,var05,var06,var07,var08,var09,var10,    &
             var11,var12,var13,var14,var15,var16,var17,var18,var19,var20,    &
             var21,var22,var23,var24,                                        &
             idim, jdim, kdim, izlo_us, izup_us, jzlo_us, jzup_us,           &
             ierror, yerrmsg, imp_type, ncount, type_handle )
          iexchg_MPI_types(2*icase) = type_handle
          iexchg_counts   (2*icase) = ncount
        ELSE
          CALL setup_data_type                                               &
           ( var01,var02,var03,var04,var05,var06,var07,var08,var09,var10,    &
             var11,var12,var13,var14,var15,var16,var17,var18,var19,var20,    &
             var21,var22,var23,var24,                                        &
             idim, jdim, kdim, izlo_dr, izup_dr, jzlo_dr, jzup_dr,           &
             ierror, yerrmsg, imp_type, ncount, type_handle )
          iexchg_MPI_types(2*icase) = type_handle
          iexchg_counts   (2*icase) = ncount
        ENDIF
      ENDIF
    ENDIF
  
  !------------------------------------------------------------------------------
  !- Section 3: Exchange with immediate Send and blocking Recv
  !------------------------------------------------------------------------------
  
    IF   (ntype == 1) THEN
  
      IF (lmpi_types) THEN
  
        !------------------------------------------------------------------------
        !- Section 3.1: exchange with left and right neighbors using datatypes
        !------------------------------------------------------------------------
  
        IF (neighbors(1) /= -1) THEN
          ! left neighbor is present
          CALL MPI_ISEND ( var01(izlo_ls,jzlo_ls,1), iexchg_counts(2*icase-1), &
                           iexchg_MPI_types(2*icase-1), MPI_neighbors(1),      &
                           ntag, icomm, ilocalreq(1), izmplcode)
          IF (izmplcode /= 0) THEN
            ierror  = izmplcode
            yerrmsg = 'MPI_ISEND'
            RETURN
          ENDIF
        ENDIF
  
        IF (neighbors(3) /= -1) THEN
          ! right neighbor is present
          ! send the data
          CALL MPI_ISEND ( var01(izlo_rs,jzlo_rs,1), iexchg_counts(2*icase-1), &
                           iexchg_MPI_types(2*icase-1), MPI_neighbors(3),      &
                           ntag, icomm, ilocalreq(3), izmplcode)
          IF (izmplcode /= 0) THEN
            ierror  = izmplcode
            yerrmsg = 'MPI_ISEND'
            RETURN
          ENDIF
        ENDIF
  
        IF (neighbors(3) /= -1) THEN
          ! right neighbor is present
          ! receive the data
          CALL MPI_RECV (var01(izlo_rr,jzlo_rr,1), iexchg_counts(2*icase-1),   &
                         iexchg_MPI_types(2*icase-1), MPI_neighbors(3),        &
                         ntag, icomm, nzrequest, izmplcode)
          IF (izmplcode /= 0) THEN
            ierror  = izmplcode
            yerrmsg = 'MPI_RECV'
            RETURN
          ENDIF
        ENDIF
  
        IF (neighbors(1) /= -1) THEN
          ! left neighbor is present
          ! receive the data
          CALL MPI_RECV ( var01(izlo_lr,jzlo_lr,1), iexchg_counts(2*icase-1),  &
                          iexchg_MPI_types(2*icase-1), MPI_neighbors(1),       &
                          ntag, icomm, nzrequest, izmplcode)
          IF (izmplcode /= 0) THEN
            ierror  = izmplcode
            yerrmsg = 'MPI_RECV'
            RETURN
          ENDIF
        ENDIF
  
        IF (neighbors(1) /= -1) THEN
          ! wait for the completion of the last send to neighbors(1)
          CALL MPI_WAIT (ilocalreq(1), nzstatus, izmplcode)
          IF (izmplcode /= 0) THEN
            ierror  = izmplcode
            yerrmsg = 'MPI_WAIT'
            RETURN
          ENDIF
        ENDIF
  
        IF (neighbors(3) /= -1) THEN
          ! wait for the completion of the last send to neighbors(3)
          CALL MPI_WAIT (ilocalreq(3), nzstatus, izmplcode)
          IF (izmplcode /= 0) THEN
            ierror  = izmplcode
            yerrmsg = 'MPI_WAIT'
            RETURN
          ENDIF
        ENDIF
  
        !------------------------------------------------------------------------
        !- Section 3.2: exchange with upper and lower neighbors using datatypes
        !------------------------------------------------------------------------
  
        IF (neighbors(2) /= -1) THEN
          ! upper neighbor is present
          ! send the data
          CALL MPI_ISEND (var01(izlo_us,jzlo_us,1), iexchg_counts(2*icase),    &
                          iexchg_MPI_types(2*icase), MPI_neighbors(2),         &
                          ntag, icomm, ilocalreq(2), izmplcode)
          IF (izmplcode /= 0) THEN
            ierror  = izmplcode
            yerrmsg = 'MPI_ISEND'
            RETURN
          ENDIF
        ENDIF
  
        IF (neighbors(4) /= -1) THEN
          ! send the data
          CALL MPI_ISEND (var01(izlo_ds,jzlo_ds,1), iexchg_counts(2*icase),    &
                          iexchg_MPI_types(2*icase), MPI_neighbors(4),         &
                          ntag, icomm, ilocalreq(4), izmplcode)
          IF (izmplcode /= 0) THEN
            ierror  = izmplcode
            yerrmsg = 'MPI_ISEND'
            RETURN
          ENDIF
        ENDIF
  
        IF (neighbors(4) /= -1) THEN
          ! lower neighbor is present
          CALL MPI_RECV (var01(izlo_dr,jzlo_dr,1), iexchg_counts(2*icase),    &
                         iexchg_MPI_types(2*icase), MPI_neighbors(4),         &
                         ntag, icomm, nzrequest, izmplcode)
          IF (izmplcode /= 0) THEN
            ierror  = izmplcode
            yerrmsg = 'MPI_RECV'
            RETURN
          ENDIF
        ENDIF
  
        IF (neighbors(2) /= -1) THEN
          ! upper neighbor is present
          CALL MPI_RECV (var01(izlo_ur,jzlo_ur,1), iexchg_counts(2*icase),     &
                          iexchg_MPI_types(2*icase), MPI_neighbors(2),         &
                          ntag, icomm, nzrequest, izmplcode)
          IF (izmplcode /= 0) THEN
            ierror  = izmplcode
            yerrmsg = 'MPI_RECV'
            RETURN
          ENDIF
        ENDIF
    
        IF (neighbors(2) /= -1) THEN
          ! wait for the completion of the last send to neighbors(2)
          CALL MPI_WAIT (ilocalreq(2), nzstatus, izmplcode)
          IF (izmplcode /= 0) THEN
            ierror  = izmplcode
            yerrmsg = 'MPI_WAIT'
            RETURN
          ENDIF
        ENDIF
    
        IF (neighbors(4) /= -1) THEN
          ! wait for the completion of the last send to neighbors(4)
          CALL MPI_WAIT (ilocalreq(4), nzstatus, izmplcode)
          IF (izmplcode /= 0) THEN
            ierror  = izmplcode
            yerrmsg = 'MPI_WAIT'
            RETURN
          ENDIF
        ENDIF
  
      ELSE
  
        !------------------------------------------------------------------------
        !- Section 3.3: exchange with left and right neigh. using explict buff.
        !------------------------------------------------------------------------
  
        IF (neighbors(1) /= -1) THEN
          ! left neighbor is present
    
          ! determine start- and end-indices for routine putbuf
          nzcount_ls = 0
          CALL putbuf ( var01, var02, var03, var04, var05, var06, var07, var08,&
                        var09, var10, var11, var12, var13, var14, var15 ,var16,&
                        var17, var18, var19, var20, var21, var22, var23 ,var24,&
                        sendbuf, isendbuflen, idim, jdim, kdim,                &
                        izlo_ls, izup_ls, jzlo_ls, jzup_ls, nzcount_ls, 1 )
  
          ! send the data
          CALL MPI_ISEND ( sendbuf(1,1), nzcount_ls, imp_type, neighbors(1),   &
                           ntag, icomm, ilocalreq(1), izmplcode)
          IF (izmplcode /= 0) THEN
            ierror  = izmplcode
            yerrmsg = 'MPI_ISEND'
            RETURN
          ENDIF
        ENDIF
  
        IF (neighbors(3) /= -1) THEN
          ! right neighbor is present
  
          ! determine start- and end-indices for routine putbuf
          nzcount_rs = 0
          CALL putbuf ( var01, var02, var03, var04, var05, var06, var07, var08,&
                        var09, var10, var11, var12, var13, var14, var15 ,var16,&
                        var17, var18, var19, var20, var21, var22, var23 ,var24,&
                        sendbuf, isendbuflen, idim, jdim, kdim,                &
                        izlo_rs, izup_rs, jzlo_rs, jzup_rs, nzcount_rs, 3 )
  
          ! send the data
          CALL MPI_ISEND ( sendbuf(1,3), nzcount_rs, imp_type, neighbors(3),   &
                           ntag, icomm, ilocalreq(3), izmplcode)
          IF (izmplcode /= 0) THEN
            ierror  = izmplcode
            yerrmsg = 'MPI_ISEND'
            RETURN
          ENDIF
        ENDIF
  
        IF (neighbors(3) /= -1) THEN
          ! receive the data
          CALL MPI_RECV ( sendbuf(1,7), isendbuflen, imp_type, neighbors(3),   &
                          ntag, icomm, nzrequest, izmplcode)
          IF (izmplcode /= 0) THEN
            ierror  = izmplcode
            yerrmsg = 'MPI_RECV'
            RETURN
          ENDIF
  
          CALL getbuf ( var01, var02, var03, var04, var05, var06, var07, var08,&
                        var09, var10, var11, var12, var13, var14, var15 ,var16,&
                        var17, var18, var19, var20, var21, var22, var23 ,var24,&
                        sendbuf, isendbuflen, idim, jdim, kdim,                &
                        izlo_rr, izup_rr, jzlo_rr, jzup_rr, nzcount_rr, 7 )
        ENDIF
  
        IF (neighbors(1) /= -1) THEN
          ! left neighbor is present
          ! receive the data
  
          CALL MPI_RECV ( sendbuf(1,5), isendbuflen, imp_type, neighbors(1),   &
                          ntag, icomm, nzrequest, izmplcode)
          IF (izmplcode /= 0) THEN
            ierror  = izmplcode
            yerrmsg = 'MPI_RECV'
            RETURN
          ENDIF
  
          CALL getbuf ( var01, var02, var03, var04, var05, var06, var07, var08,&
                        var09, var10, var11, var12, var13, var14, var15 ,var16,&
                        var17, var18, var19, var20, var21, var22, var23 ,var24,&
                        sendbuf, isendbuflen, idim, jdim, kdim,                &
                        izlo_lr, izup_lr, jzlo_lr, jzup_lr, nzcount_lr, 5 )
        ENDIF
  
        IF (neighbors(1) /= -1) THEN
          ! wait for the completion of the last send to neighbors(1)
          ! to safely reuse sendbuf(1,1)
          CALL MPI_WAIT (ilocalreq(1), nzstatus, izmplcode)
          IF (izmplcode /= 0) THEN
            ierror  = izmplcode
            yerrmsg = 'MPI_WAIT'
            RETURN
          ENDIF
        ENDIF
  
        IF (neighbors(3) /= -1) THEN
          ! wait for the completion of the last send to neighbors(3)
          ! to safely reuse sendbuf(1,3)
          CALL MPI_WAIT (ilocalreq(3), nzstatus, izmplcode)
            IF (izmplcode /= 0) THEN
              ierror  = izmplcode
              yerrmsg = 'MPI_WAIT'
              RETURN
            ENDIF
        ENDIF
  
        !------------------------------------------------------------------------
        !- Section 3.4: exchange with lower and upper neigh. using explict buff.
        !------------------------------------------------------------------------
  
        IF (neighbors(2) /= -1) THEN
          nzcount_us = 0
          CALL putbuf ( var01, var02, var03, var04, var05, var06, var07, var08,&
                        var09, var10, var11, var12, var13, var14, var15 ,var16,&
                        var17, var18, var19, var20, var21, var22, var23 ,var24,&
                        sendbuf, isendbuflen, idim, jdim, kdim,                &
                        izlo_us, izup_us, jzlo_us, jzup_us, nzcount_us, 2 )
  
          ! send the data
          CALL MPI_ISEND ( sendbuf(1,2), nzcount_us, imp_type, neighbors(2),   &
                           ntag, icomm, ilocalreq(2), izmplcode)
          IF (izmplcode /= 0) THEN
            ierror  = izmplcode
            yerrmsg = 'MPI_ISEND'
            RETURN
          ENDIF
        ENDIF
  
        IF (neighbors(4) /= -1) THEN
          ! lower neighbor is present
  
          nzcount_ds = 0
          CALL putbuf ( var01, var02, var03, var04, var05, var06, var07, var08,&
                        var09, var10, var11, var12, var13, var14, var15 ,var16,&
                        var17, var18, var19, var20, var21, var22, var23 ,var24,&
                        sendbuf, isendbuflen, idim, jdim, kdim,                &
                        izlo_ds, izup_ds, jzlo_ds, jzup_ds, nzcount_ds, 4 )
  
          ! send the data
          CALL MPI_ISEND ( sendbuf(1,4), nzcount_ds, imp_type, neighbors(4),   &
                           ntag, icomm, ilocalreq(4), izmplcode)
          IF (izmplcode /= 0) THEN
            ierror  = izmplcode
            yerrmsg = 'MPI_ISEND'
            RETURN
          ENDIF
        ENDIF
  
        IF (neighbors(4) /= -1) THEN
          ! lower neighbor is present
          ! receive the data
    
          CALL MPI_RECV ( sendbuf(1,8), isendbuflen, imp_type, neighbors(4),   &
                          ntag, icomm, nzrequest, izmplcode)
          IF (izmplcode /= 0) THEN
            ierror  = izmplcode
            yerrmsg = 'MPI_RECV'
            RETURN
          ENDIF
    
          CALL getbuf ( var01, var02, var03, var04, var05, var06, var07, var08,&
                        var09, var10, var11, var12, var13, var14, var15 ,var16,&
                        var17, var18, var19, var20, var21, var22, var23 ,var24,&
                        sendbuf, isendbuflen, idim, jdim, kdim,                &
                        izlo_dr, izup_dr, jzlo_dr, jzup_dr, nzcount_dr, 8 )
        ENDIF
    
        IF (neighbors(2) /= -1) THEN
          ! upper neighbor is present
          ! receive the data
    
          CALL MPI_RECV ( sendbuf(1,6), isendbuflen, imp_type, neighbors(2),   &
                          ntag, icomm, nzrequest, izmplcode)
          IF (izmplcode /= 0) THEN
            ierror  = izmplcode
            yerrmsg = 'MPI_RECV'
            RETURN
          ENDIF
    
          CALL getbuf ( var01, var02, var03, var04, var05, var06, var07, var08,&
                        var09, var10, var11, var12, var13, var14, var15 ,var16,&
                        var17, var18, var19, var20, var21, var22, var23 ,var24,&
                        sendbuf, isendbuflen, idim, jdim, kdim,                &
                        izlo_ur, izup_ur, jzlo_ur, jzup_ur, nzcount_ur, 6 )
        ENDIF
    
        IF (neighbors(2) /= -1) THEN
          ! wait for the completion of the last send to neighbors(2)
          ! to safely reuse sendbuf(1,2)
          CALL MPI_WAIT (ilocalreq(2), nzstatus, izmplcode)
          IF (izmplcode /= 0) THEN
            ierror  = izmplcode
            yerrmsg = 'MPI_WAIT'
            RETURN
          ENDIF
        ENDIF
    
        IF (neighbors(4) /= -1) THEN
          ! wait for the completion of the last send to neighbors(4)
          ! to safely reuse sendbuf(1,4)
          CALL MPI_WAIT (ilocalreq(4), nzstatus, izmplcode)
          IF (izmplcode /= 0) THEN
            ierror  = izmplcode
            yerrmsg = 'MPI_WAIT'
            RETURN
          ENDIF
        ENDIF
  
      ENDIF
  
  !------------------------------------------------------------------------------
  !- Section 4: Exchange with immediate Recv and blocking Send
  !------------------------------------------------------------------------------
  
    ELSEIF (ntype == 2) THEN
  
      IF (lmpi_types) THEN
  
        !------------------------------------------------------------------------
        !- Section 4.1: exchange with left and right neighbors
        !------------------------------------------------------------------------
  
        IF (neighbors(3) /= -1) THEN
          CALL MPI_IRECV (var01(izlo_rr,jzlo_rr,1), iexchg_counts(2*icase-1), &
                          iexchg_MPI_types(2*icase-1), MPI_neighbors(3),      &
                          ntag, icomm, ilocalreq(3), izmplcode)
          IF (izmplcode /= 0) THEN
            ierror  = izmplcode
            yerrmsg = 'MPI_IRECV'
            RETURN
          ENDIF
        ENDIF
  
        IF (neighbors(1) /= -1) THEN
          CALL MPI_IRECV (var01(izlo_lr,jzlo_lr,1), iexchg_counts(2*icase-1),&
                          iexchg_MPI_types(2*icase-1), MPI_neighbors(1),     &
                          ntag, icomm, ilocalreq(1), izmplcode)
          IF (izmplcode /= 0) THEN
            ierror  = izmplcode
            yerrmsg = 'MPI_IRECV'
            RETURN
          ENDIF
        ENDIF
  
        IF (neighbors(1) /= -1) THEN
          ! now send the data to the left neighbor
          CALL MPI_SEND ( var01(izlo_ls,jzlo_ls,1), iexchg_counts(2*icase-1),&
                          iexchg_MPI_types(2*icase-1), MPI_neighbors(1),     &
                          ntag, icomm, izmplcode)
          IF (izmplcode /= 0) THEN
            ierror  = izmplcode
            yerrmsg = 'MPI_SEND'
            RETURN
          ENDIF
        ENDIF
  
        IF (neighbors(3) /= -1) THEN
          ! now send the data to the right neighbor
          CALL MPI_SEND ( var01(izlo_rs,jzlo_rs,1), iexchg_counts(2*icase-1),&
                          iexchg_MPI_types(2*icase-1), MPI_neighbors(3),     &
                          ntag, icomm, izmplcode)
          IF (izmplcode /= 0) THEN
            ierror  = izmplcode
            yerrmsg = 'MPI_SEND'
            RETURN
          ENDIF
        ENDIF
  
        IF (neighbors(3) /= -1) THEN
          CALL MPI_WAIT (ilocalreq(3), nzstatus, izmplcode)
          IF (izmplcode /= 0) THEN
            ierror  = izmplcode
            yerrmsg = 'MPI_WAIT'
            RETURN
          ENDIF
        ENDIF
  
        IF (neighbors(1) /= -1) THEN
          CALL MPI_WAIT (ilocalreq(1), nzstatus, izmplcode)
          IF (izmplcode /= 0) THEN
            ierror  = izmplcode
            yerrmsg = 'MPI_WAIT'
            RETURN
          ENDIF
        ENDIF
  
        !------------------------------------------------------------------------
        !- Section 4.2: exchange with upper and lower neighbors
        !------------------------------------------------------------------------
  
        IF (neighbors(2) /= -1) THEN
          CALL MPI_IRECV (var01(izlo_ur,jzlo_ur,1), iexchg_counts(2*icase),  &
                          iexchg_MPI_types(2*icase), MPI_neighbors(2),       &
                          ntag, icomm, ilocalreq(2), izmplcode)
          IF (izmplcode /= 0) THEN
            ierror  = izmplcode
            yerrmsg = 'MPI_IRECV'
            RETURN
          ENDIF
        ENDIF
  
        IF (neighbors(4) /= -1) THEN
          CALL MPI_IRECV (var01(izlo_dr,jzlo_dr,1), iexchg_counts(2*icase),  &
                          iexchg_MPI_types(2*icase), MPI_neighbors(4),       &
                          ntag, icomm, ilocalreq(4), izmplcode)
          IF (izmplcode /= 0) THEN
            ierror  = izmplcode
            yerrmsg = 'MPI_IRECV'
            RETURN
          ENDIF
        ENDIF
  
        IF (neighbors(4) /= -1) THEN
          CALL MPI_SEND (var01(izlo_ds,jzlo_ds,1), iexchg_counts(2*icase),  &
                         iexchg_MPI_types(2*icase), MPI_neighbors(4),       &
                         ntag, icomm, izmplcode)
          IF (izmplcode /= 0) THEN
            ierror  = izmplcode
            yerrmsg = 'MPI_SEND'
            RETURN
          ENDIF
        ENDIF
  
        IF (neighbors(2) /= -1) THEN
          CALL MPI_SEND (var01(izlo_us,jzlo_us,1), iexchg_counts(2*icase),  &
                         iexchg_MPI_types(2*icase), MPI_neighbors(2),       &
                         ntag, icomm, izmplcode)
          IF (izmplcode /= 0) THEN
            ierror  = izmplcode
            yerrmsg = 'MPI_SEND'
            RETURN
          ENDIF
        ENDIF
  
        IF (neighbors(2) /= -1) THEN
          CALL MPI_WAIT (ilocalreq(2), nzstatus, izmplcode)
          IF (izmplcode /= 0) THEN
            ierror  = izmplcode
            yerrmsg = 'MPI_WAIT'
            RETURN
          ENDIF
        ENDIF
  
        IF (neighbors(4) /= -1) THEN
          CALL MPI_WAIT (ilocalreq(4), nzstatus, izmplcode)
          IF (izmplcode /= 0) THEN
            ierror  = izmplcode
            yerrmsg = 'MPI_WAIT'
            RETURN
          ENDIF
        ENDIF
  
      ELSE
  
        !------------------------------------------------------------------------
        !- Section 4.3: exchange with left and right neighbors
        !------------------------------------------------------------------------
  
        IF (neighbors(3) /= -1) THEN
          CALL MPI_IRECV ( sendbuf(1,7), isendbuflen, imp_type, neighbors(3),  &
                          ntag, icomm, ilocalreq(3), izmplcode)
          IF (izmplcode /= 0) THEN
            ierror  = izmplcode
            yerrmsg = 'MPI_IRECV'
            RETURN
          ENDIF
        ENDIF
  
        IF (neighbors(1) /= -1) THEN
          CALL MPI_IRECV ( sendbuf(1,5), isendbuflen, imp_type, neighbors(1),  &
                          ntag, icomm, ilocalreq(1), izmplcode)
          IF (izmplcode /= 0) THEN
            ierror  = izmplcode
            yerrmsg = 'MPI_IRECV'
            RETURN
          ENDIF
        ENDIF
  
        IF (neighbors(1) /= -1) THEN
          nzcount_ls = 0
          CALL putbuf ( var01, var02, var03, var04, var05, var06, var07, var08,&
                        var09, var10, var11, var12, var13, var14, var15 ,var16,&
                        var17, var18, var19, var20, var21, var22, var23 ,var24,&
                        sendbuf, isendbuflen, idim, jdim, kdim,                &
                        izlo_ls, izup_ls, jzlo_ls, jzup_ls, nzcount_ls, 1 )
  
          CALL MPI_SEND ( sendbuf(1,1), nzcount_ls, imp_type, neighbors(1),    &
                          ntag, icomm, izmplcode)
          IF (izmplcode /= 0) THEN
            ierror  = izmplcode
            yerrmsg = 'MPI_SEND'
            RETURN
          ENDIF
        ENDIF
  
        IF (neighbors(3) /= -1) THEN
          nzcount_rs = 0
          CALL putbuf ( var01, var02, var03, var04, var05, var06, var07, var08,&
                        var09, var10, var11, var12, var13, var14, var15 ,var16,&
                        var17, var18, var19, var20, var21, var22, var23 ,var24,&
                        sendbuf, isendbuflen, idim, jdim, kdim,                &
                        izlo_rs, izup_rs, jzlo_rs, jzup_rs, nzcount_rs, 3 )
  
          CALL MPI_SEND ( sendbuf(1,3), nzcount_rs, imp_type, neighbors(3),    &
                          ntag, icomm, izmplcode)
          IF (izmplcode /= 0) THEN
            ierror  = izmplcode
            yerrmsg = 'MPI_SEND'
            RETURN
          ENDIF
        ENDIF
  
        IF (neighbors(3) /= -1) THEN
          CALL MPI_WAIT (ilocalreq(3), nzstatus, izmplcode)
          IF (izmplcode /= 0) THEN
            ierror  = izmplcode
            yerrmsg = 'MPI_WAIT'
            RETURN
          ENDIF
  
          CALL getbuf ( var01, var02, var03, var04, var05, var06, var07, var08,&
                        var09, var10, var11, var12, var13, var14, var15 ,var16,&
                        var17, var18, var19, var20, var21, var22, var23 ,var24,&
                        sendbuf, isendbuflen, idim, jdim, kdim,                &
                        izlo_rr, izup_rr, jzlo_rr, jzup_rr, nzcount_rr, 7 )
        ENDIF
  
        IF (neighbors(1) /= -1) THEN
          CALL MPI_WAIT (ilocalreq(1), nzstatus, izmplcode)
          IF (izmplcode /= 0) THEN
            ierror  = izmplcode
            yerrmsg = 'MPI_WAIT'
            RETURN
          ENDIF
  
          CALL getbuf ( var01, var02, var03, var04, var05, var06, var07, var08,&
                        var09, var10, var11, var12, var13, var14, var15 ,var16,&
                        var17, var18, var19, var20, var21, var22, var23 ,var24,&
                        sendbuf, isendbuflen, idim, jdim, kdim,                &
                        izlo_lr, izup_lr, jzlo_lr, jzup_lr, nzcount_lr, 5 )
        ENDIF
  
        !------------------------------------------------------------------------
        !- Section 4.4: exchange with upper and lower neighbors
        !------------------------------------------------------------------------
  
        IF (neighbors(2) /= -1) THEN
          CALL MPI_IRECV ( sendbuf(1,6), isendbuflen, imp_type, neighbors(2),  &
                          ntag, icomm, ilocalreq(2), izmplcode)
          IF (izmplcode /= 0) THEN
            ierror  = izmplcode
            yerrmsg = 'MPI_IRECV'
            RETURN
          ENDIF
        ENDIF
  
        IF (neighbors(4) /= -1) THEN
          CALL MPI_IRECV ( sendbuf(1,8), isendbuflen, imp_type, neighbors(4),  &
                          ntag, icomm, ilocalreq(4), izmplcode)
          IF (izmplcode /= 0) THEN
            ierror  = izmplcode
            yerrmsg = 'MPI_IRECV'
            RETURN
          ENDIF
        ENDIF
  
        IF (neighbors(4) /= -1) THEN
          nzcount_ds = 0
          CALL putbuf ( var01, var02, var03, var04, var05, var06, var07, var08,&
                        var09, var10, var11, var12, var13, var14, var15 ,var16,&
                        var17, var18, var19, var20, var21, var22, var23 ,var24,&
                        sendbuf, isendbuflen, idim, jdim, kdim,                &
                        izlo_ds, izup_ds, jzlo_ds, jzup_ds, nzcount_ds, 4 )
  
          CALL MPI_SEND ( sendbuf(1,4), nzcount_ds, imp_type, neighbors(4),    &
                          ntag, icomm, izmplcode)
          IF (izmplcode /= 0) THEN
            ierror  = izmplcode
            yerrmsg = 'MPI_SEND'
            RETURN
          ENDIF
        ENDIF
  
        IF (neighbors(2) /= -1) THEN
          nzcount_us = 0
          CALL putbuf ( var01, var02, var03, var04, var05, var06, var07, var08,&
                        var09, var10, var11, var12, var13, var14, var15 ,var16,&
                        var17, var18, var19, var20, var21, var22, var23 ,var24,&
                        sendbuf, isendbuflen, idim, jdim, kdim,                &
                        izlo_us, izup_us, jzlo_us, jzup_us, nzcount_us, 2 )
   
          CALL MPI_SEND ( sendbuf(1,2), nzcount_us, imp_type, neighbors(2),    &
                          ntag, icomm, izmplcode)
          IF (izmplcode /= 0) THEN
            ierror  = izmplcode
            yerrmsg = 'MPI_SEND'
            RETURN
          ENDIF
        ENDIF
  
        IF (neighbors(2) /= -1) THEN
          CALL MPI_WAIT (ilocalreq(2), nzstatus, izmplcode)
          IF (izmplcode /= 0) THEN
            ierror  = izmplcode
            yerrmsg = 'MPI_WAIT'
            RETURN
          ENDIF
  
          CALL getbuf ( var01, var02, var03, var04, var05, var06, var07, var08,&
                        var09, var10, var11, var12, var13, var14, var15 ,var16,&
                        var17, var18, var19, var20, var21, var22, var23 ,var24,&
                        sendbuf, isendbuflen, idim, jdim, kdim,                &
                        izlo_ur, izup_ur, jzlo_ur, jzup_ur, nzcount_ur, 6 )
        ENDIF
  
   
        IF (neighbors(4) /= -1) THEN
          ! Now wait until the data have arrived
          CALL MPI_WAIT (ilocalreq(4), nzstatus, izmplcode)
          IF (izmplcode /= 0) THEN
            ierror  = izmplcode
            yerrmsg = 'MPI_WAIT'
            RETURN
          ENDIF
  
          CALL getbuf ( var01, var02, var03, var04, var05, var06, var07, var08,&
                        var09, var10, var11, var12, var13, var14, var15 ,var16,&
                        var17, var18, var19, var20, var21, var22, var23 ,var24,&
                        sendbuf, isendbuflen, idim, jdim, kdim,                &
                        izlo_dr, izup_dr, jzlo_dr, jzup_dr, nzcount_dr, 8 )
        ENDIF
  
      ENDIF
  
  !------------------------------------------------------------------------------
  !- Section 5: Exchange with SendRecv
  !------------------------------------------------------------------------------
  
    ELSEIF (ntype == 3) THEN
  
      IF (lmpi_types) THEN
  
        !--------------------------------------------------------------------------
        !- Section 5.1: Send data to the left and receive from the right neighbor
        !--------------------------------------------------------------------------
   
        CALL MPI_SENDRECV ( var01(izlo_ls,jzlo_ls,1), iexchg_counts(2*icase-1),&
                            iexchg_MPI_types(2*icase-1), MPI_neighbors(1),     &
                            ntag,                                              &
                            var01(izlo_rr,jzlo_rr,1), iexchg_counts(2*icase-1),&
                            iexchg_MPI_types(2*icase-1), MPI_neighbors(3),     &
                            ntag,                                              &
                            icomm, nzstatus, izmplcode)
        IF (izmplcode /= 0) THEN
          ierror  = izmplcode
          yerrmsg = 'MPI_SENDRECV-1'
          RETURN
        ENDIF
  
        !--------------------------------------------------------------------------
        !- Section 5.2: Receive data from the left and send to the right neighbor
        !--------------------------------------------------------------------------
  
        CALL MPI_SENDRECV ( var01(izlo_rs,jzlo_rs,1), iexchg_counts(2*icase-1),&
                            iexchg_MPI_types(2*icase-1), MPI_neighbors(3),     &
                            ntag+1,                                            &
                            var01(izlo_lr,jzlo_lr,1), iexchg_counts(2*icase-1),&
                            iexchg_MPI_types(2*icase-1), MPI_neighbors(1),     &
                            ntag+1,                                            &
                            icomm, nzstatus, izmplcode)
        IF (izmplcode /= 0) THEN
          ierror  = izmplcode
          yerrmsg = 'MPI_SENDRECV-2'
          RETURN
        ENDIF
  
        !--------------------------------------------------------------------------
        !- Section 5.3: Send data to the upper and receive from the lower neighbor
        !--------------------------------------------------------------------------
  
        CALL MPI_SENDRECV ( var01(izlo_us,jzlo_us,1), iexchg_counts(2*icase),  &
                            iexchg_MPI_types(2*icase), MPI_neighbors(2),       &
                            ntag+2,                                            &
                            var01(izlo_dr,jzlo_dr,1), iexchg_counts(2*icase),  &
                            iexchg_MPI_types(2*icase), MPI_neighbors(4),       &
                            ntag+2,                                            &
                            icomm, nzstatus, izmplcode)
        IF (izmplcode /= 0) THEN
          ierror  = izmplcode
          yerrmsg = 'MPI_SENDRECV-3'
          RETURN
        ENDIF
   
  
        !--------------------------------------------------------------------------
        !- Section 5.4: Receive data from the upper and send to the lower neighbor
        !--------------------------------------------------------------------------
   
        CALL MPI_SENDRECV ( var01(izlo_ds,jzlo_ds,1), iexchg_counts(2*icase),  &
                            iexchg_MPI_types(2*icase), MPI_neighbors(4),       &
                            ntag+3,                                            &
                            var01(izlo_ur,jzlo_ur,1), iexchg_counts(2*icase),  &
                            iexchg_MPI_types(2*icase), MPI_neighbors(2),       &
                            ntag+3,                                            &
                            icomm, nzstatus, izmplcode)
        IF (izmplcode /= 0) THEN
          ierror  = izmplcode
          yerrmsg = 'MPI_SENDRECV-4'
          RETURN
        ENDIF
  
      ELSE
  
        !--------------------------------------------------------------------------
        !- Section 5.5: Send data to the left and receive from the right neighbor
        !--------------------------------------------------------------------------
  
        nzcount_ls = 0
        IF (MPI_neighbors(1) /= MPI_PROC_NULL) THEN
          CALL putbuf ( var01, var02, var03, var04, var05, var06, var07, var08,&
                        var09, var10, var11, var12, var13, var14, var15 ,var16,&
                        var17, var18, var19, var20, var21, var22, var23 ,var24,&
                        sendbuf, isendbuflen, idim, jdim, kdim,                &
                        izlo_ls, izup_ls, jzlo_ls, jzup_ls, nzcount_ls, 1 )
        ENDIF
  
        CALL MPI_SENDRECV                                                      &
             ( sendbuf(1,1), nzcount_ls,  imp_type, MPI_neighbors(1), ntag,    &
               sendbuf(1,7), isendbuflen, imp_type, MPI_neighbors(3), ntag,    &
               icomm, nzstatus, izmplcode)
        IF (izmplcode /= 0) THEN
          ierror  = izmplcode
          yerrmsg = 'MPI_SENDRECV'
          RETURN
        ENDIF
  
        IF (MPI_neighbors(3) /= MPI_PROC_NULL) THEN
          CALL getbuf ( var01, var02, var03, var04, var05, var06, var07, var08,&
                        var09, var10, var11, var12, var13, var14, var15 ,var16,&
                        var17, var18, var19, var20, var21, var22, var23 ,var24,&
                        sendbuf, isendbuflen, idim, jdim, kdim,                &
                        izlo_rr, izup_rr, jzlo_rr, jzup_rr, nzcount_rr, 7 )
        ENDIF
  
        !--------------------------------------------------------------------------
        !- Section 5.6: Send data to the right and receive from the left neighbor
        !--------------------------------------------------------------------------
  
        nzcount_rs = 0
        IF (MPI_neighbors(3) /= MPI_PROC_NULL) THEN
          CALL putbuf ( var01, var02, var03, var04, var05, var06, var07, var08,&
                        var09, var10, var11, var12, var13, var14, var15 ,var16,&
                        var17, var18, var19, var20, var21, var22, var23 ,var24,&
                        sendbuf, isendbuflen, idim, jdim, kdim,                &
                        izlo_rs, izup_rs, jzlo_rs, jzup_rs, nzcount_rs, 3 )
        ENDIF
  
        CALL MPI_SENDRECV                                                      &
             ( sendbuf(1,3), nzcount_rs,  imp_type, MPI_neighbors(3), ntag,    &
               sendbuf(1,5), isendbuflen, imp_type, MPI_neighbors(1), ntag,    &
               icomm, nzstatus, izmplcode)
        IF (izmplcode /= 0) THEN
          ierror  = izmplcode
          yerrmsg = 'MPI_SENDRECV'
          RETURN
        ENDIF
  
        IF (MPI_neighbors(1) /= MPI_PROC_NULL) THEN
          CALL getbuf ( var01, var02, var03, var04, var05, var06, var07, var08,&
                        var09, var10, var11, var12, var13, var14, var15 ,var16,&
                        var17, var18, var19, var20, var21, var22, var23 ,var24,&
                        sendbuf, isendbuflen, idim, jdim, kdim,                &
                        izlo_lr, izup_lr, jzlo_lr, jzup_lr, nzcount_lr, 5 )
        ENDIF
  
        !--------------------------------------------------------------------------
        !- Section 5.7: Send data to the upper and receive from the lower neighbor
        !--------------------------------------------------------------------------
   
        nzcount_us = 0
        IF (MPI_neighbors(2) /= MPI_PROC_NULL) THEN
          CALL putbuf ( var01, var02, var03, var04, var05, var06, var07, var08,&
                        var09, var10, var11, var12, var13, var14, var15 ,var16,&
                        var17, var18, var19, var20, var21, var22, var23 ,var24,&
                        sendbuf, isendbuflen, idim, jdim, kdim,                &
                        izlo_us, izup_us, jzlo_us, jzup_us, nzcount_us, 2 )
        ENDIF
  
        CALL MPI_SENDRECV                                                      &
             ( sendbuf(1,2), nzcount_us,  imp_type, MPI_neighbors(2), ntag,    &
               sendbuf(1,8), isendbuflen, imp_type, MPI_neighbors(4), ntag,    &
               icomm, nzstatus, izmplcode)
        IF (izmplcode /= 0) THEN
          ierror  = izmplcode
          yerrmsg = 'MPI_SENDRECV'
          RETURN
        ENDIF
  
        IF (MPI_neighbors(4) /= MPI_PROC_NULL) THEN
          CALL getbuf ( var01, var02, var03, var04, var05, var06, var07, var08,&
                        var09, var10, var11, var12, var13, var14, var15 ,var16,&
                        var17, var18, var19, var20, var21, var22, var23 ,var24,&
                        sendbuf, isendbuflen, idim, jdim, kdim,                &
                        izlo_dr, izup_dr, jzlo_dr, jzup_dr, nzcount_dr, 8 )
        ENDIF
  
        !--------------------------------------------------------------------------
        !- Section 5.8: Send data to the lower and receive from the upper neighbor
        !--------------------------------------------------------------------------
   
        nzcount_ds = 0
        IF (MPI_neighbors(4) /= MPI_PROC_NULL) THEN
          CALL putbuf ( var01, var02, var03, var04, var05, var06, var07, var08,&
                        var09, var10, var11, var12, var13, var14, var15 ,var16,&
                        var17, var18, var19, var20, var21, var22, var23 ,var24,&
                        sendbuf, isendbuflen, idim, jdim, kdim,                &
                        izlo_ds, izup_ds, jzlo_ds, jzup_ds, nzcount_ds, 4 )
        ENDIF
  
        CALL MPI_SENDRECV                                                      &
             ( sendbuf(1,4), nzcount_ds,  imp_type, MPI_neighbors(4), ntag,    &
               sendbuf(1,6), isendbuflen, imp_type, MPI_neighbors(2), ntag,    &
               icomm, nzstatus, izmplcode)
        IF (izmplcode /= 0) THEN
          ierror  = izmplcode
          yerrmsg = 'MPI_SENDRECV'
          RETURN
        ENDIF
  
        IF (MPI_neighbors(2) /= MPI_PROC_NULL) THEN
          CALL getbuf ( var01, var02, var03, var04, var05, var06, var07, var08,&
                        var09, var10, var11, var12, var13, var14, var15 ,var16,&
                        var17, var18, var19, var20, var21, var22, var23 ,var24,&
                        sendbuf, isendbuflen, idim, jdim, kdim,                &
                        izlo_ur, izup_ur, jzlo_ur, jzup_ur, nzcount_ur, 6 )
        ENDIF
  
      ENDIF
  
    ENDIF

  END IF  ! num_compute > 1

!------------------------------------------------------------------------------
!- End of the Subroutine
!------------------------------------------------------------------------------

END SUBROUTINE exchg_boundaries

!==============================================================================
!==============================================================================
!+ This subroutine puts all necessary values into sendbuf
!------------------------------------------------------------------------------

SUBROUTINE putbuf  (var01, var02, var03, var04, var05, var06, var07, var08, &
                    var09, var10, var11, var12, var13, var14, var15, var16, &
                    var17, var18, var19, var20, var21, var22, var23, var24, &
                    sendbuf, isendbuflen, idim, jdim, kdim,                 &
                    ilo, iup, jlo, jup, ncount, nentry )

!------------------------------------------------------------------------------
!
! Description:
!   This subroutine puts the necessary values from the present variables
!   (determined by ilo, iup, jlo, jup) into sendbuf.
!
! Method:
!   Check which variables are present.
!
!------------------------------------------------------------------------------

! Subroutine arguments
  INTEGER (KIND=iintegers), INTENT (IN)         ::    &
    isendbuflen,                  & ! length of sendbuffer
    idim, jdim, kdim(24),         & ! dimensions of the fields
    ilo, iup, jlo, jup,           & ! start- and end-indices
    nentry                          ! specifies the row of the sendbuf

  INTEGER (KIND=iintegers), INTENT (INOUT)      ::    &
    ncount                          ! counts the variables

  REAL (KIND=wp),     INTENT (INOUT)            ::    &
    sendbuf (isendbuflen, 8),     & ! send buffer
    var01  (idim, jdim, kdim( 1))   ! first field that has to occur

  REAL (KIND=wp),     OPTIONAL, INTENT (INOUT)  ::    &
    var02  (idim, jdim, kdim( 2)),& ! additional optional fields
    var03  (idim, jdim, kdim( 3)),& ! additional optional fields
    var04  (idim, jdim, kdim( 4)),& ! additional optional fields
    var05  (idim, jdim, kdim( 5)),& ! additional optional fields
    var06  (idim, jdim, kdim( 6)),& ! additional optional fields
    var07  (idim, jdim, kdim( 7)),& ! additional optional fields
    var08  (idim, jdim, kdim( 8)),& ! additional optional fields
    var09  (idim, jdim, kdim( 9)),& ! additional optional fields
    var10  (idim, jdim, kdim(10)),& ! additional optional fields
    var11  (idim, jdim, kdim(11)),& ! additional optional fields
    var12  (idim, jdim, kdim(12)),& ! additional optional fields
    var13  (idim, jdim, kdim(13)),& ! additional optional fields
    var14  (idim, jdim, kdim(14)),& ! additional optional fields
    var15  (idim, jdim, kdim(15)),& ! additional optional fields
    var16  (idim, jdim, kdim(16)),& ! additional optional fields
    var17  (idim, jdim, kdim(17)),& ! additional optional fields
    var18  (idim, jdim, kdim(18)),& ! additional optional fields
    var19  (idim, jdim, kdim(19)),& ! additional optional fields
    var20  (idim, jdim, kdim(20)),& ! additional optional fields
    var21  (idim, jdim, kdim(21)),& ! additional optional fields
    var22  (idim, jdim, kdim(22)),& ! additional optional fields
    var23  (idim, jdim, kdim(23)),& ! additional optional fields
    var24  (idim, jdim, kdim(24))   ! additional optional fields

! Local variables

  INTEGER (KIND=iintegers)   ::       &
    i, j, k, nzc

  LOGICAL                    ::       &
    lpres02, lpres03, lpres04, lpres05, lpres06, lpres07, lpres08, lpres09, &
    lpres10, lpres11, lpres12, lpres13, lpres14, lpres15, lpres16, lpres17, &
    lpres18, lpres19, lpres20, lpres21, lpres22, lpres23, lpres24

!------------------------------------------------------------------------------
!- End of header -
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!- Section 1: Initializations
!------------------------------------------------------------------------------

  ! check which variables are present
  lpres02 = PRESENT (var02)
  lpres03 = PRESENT (var03)
  lpres04 = PRESENT (var04)
  lpres05 = PRESENT (var05)
  lpres06 = PRESENT (var06)
  lpres07 = PRESENT (var07)
  lpres08 = PRESENT (var08)
  lpres09 = PRESENT (var09)
  lpres10 = PRESENT (var10)
  lpres11 = PRESENT (var11)
  lpres12 = PRESENT (var12)
  lpres13 = PRESENT (var13)
  lpres14 = PRESENT (var14)
  lpres15 = PRESENT (var15)
  lpres16 = PRESENT (var16)
  lpres17 = PRESENT (var17)
  lpres18 = PRESENT (var18)
  lpres19 = PRESENT (var19)
  lpres20 = PRESENT (var20)
  lpres21 = PRESENT (var21)
  lpres22 = PRESENT (var22)
  lpres23 = PRESENT (var23)
  lpres24 = PRESENT (var24)

!------------------------------------------------------------------------------
!- Section 2: Put data into the buffer
!------------------------------------------------------------------------------

  ! use nzc as a local counter  (based on a work from Mike O'Neill to 
  ! improve vectorization of putbuf and getbuf)
  nzc = ncount

  ! first variable that has to be present
  DO k = 1, kdim( 1)
    DO j = jlo, jup
      DO i = ilo, iup
        nzc = nzc + 1
        sendbuf (nzc,nentry) = var01(i,j,k)
      ENDDO
    ENDDO
  ENDDO

  ! optional variables that are present
  IF (lpres02 .EQV. .TRUE.) THEN
    DO k = 1, kdim( 2)
      DO j = jlo, jup
        DO i = ilo, iup
          nzc = nzc + 1
          sendbuf (nzc,nentry) = var02(i,j,k)
        ENDDO
      ENDDO
    ENDDO
  ENDIF
  
  ! optional variables that are present
  IF (lpres03 .EQV. .TRUE.) THEN
    DO k = 1, kdim( 3)
      DO j = jlo, jup
        DO i = ilo, iup
          nzc = nzc + 1
          sendbuf (nzc,nentry) = var03(i,j,k)
        ENDDO
      ENDDO
    ENDDO
  ENDIF
  
  ! optional variables that are present
  IF (lpres04 .EQV. .TRUE.) THEN
    DO k = 1, kdim( 4)
      DO j = jlo, jup
        DO i = ilo, iup
          nzc = nzc + 1
          sendbuf (nzc,nentry) = var04(i,j,k)
        ENDDO
      ENDDO
    ENDDO
  ENDIF
  
  ! optional variables that are present
  IF (lpres05 .EQV. .TRUE.) THEN
    DO k = 1, kdim( 5)
      DO j = jlo, jup
        DO i = ilo, iup
          nzc = nzc + 1
          sendbuf (nzc,nentry) = var05(i,j,k)
        ENDDO
      ENDDO
    ENDDO
  ENDIF
  
  ! optional variables that are present
  IF (lpres06 .EQV. .TRUE.) THEN
    DO k = 1, kdim( 6)
      DO j = jlo, jup
        DO i = ilo, iup
          nzc = nzc + 1
          sendbuf (nzc,nentry) = var06(i,j,k)
        ENDDO
      ENDDO
    ENDDO
  ENDIF
  
  ! optional variables that are present
  IF (lpres07 .EQV. .TRUE.) THEN
    DO k = 1, kdim( 7)
      DO j = jlo, jup
        DO i = ilo, iup
          nzc = nzc + 1
          sendbuf (nzc,nentry) = var07(i,j,k)
        ENDDO
      ENDDO
    ENDDO
  ENDIF
  
  ! optional variables that are present
  IF (lpres08 .EQV. .TRUE.) THEN
    DO k = 1, kdim( 8)
      DO j = jlo, jup
        DO i = ilo, iup
          nzc = nzc + 1
          sendbuf (nzc,nentry) = var08(i,j,k)
        ENDDO
      ENDDO
    ENDDO
  ENDIF
  
  ! optional variables that are present
  IF (lpres09 .EQV. .TRUE.) THEN
    DO k = 1, kdim( 9)
      DO j = jlo, jup
        DO i = ilo, iup
          nzc = nzc + 1
          sendbuf (nzc,nentry) = var09(i,j,k)
        ENDDO
      ENDDO
    ENDDO
  ENDIF
  
  ! optional variables that are present
  IF (lpres10 .EQV. .TRUE.) THEN
    DO k = 1, kdim(10)
      DO j = jlo, jup
        DO i = ilo, iup
          nzc = nzc + 1
          sendbuf (nzc,nentry) = var10(i,j,k)
        ENDDO
      ENDDO
    ENDDO
  ENDIF
  
  ! optional variables that are present
  IF (lpres11 .EQV. .TRUE.) THEN
    DO k = 1, kdim(11)
      DO j = jlo, jup
        DO i = ilo, iup
          nzc = nzc + 1
          sendbuf (nzc,nentry) = var11(i,j,k)
        ENDDO
      ENDDO
    ENDDO
  ENDIF
  
  ! optional variables that are present
  IF (lpres12 .EQV. .TRUE.) THEN
    DO k = 1, kdim(12)
      DO j = jlo, jup
        DO i = ilo, iup
          nzc = nzc + 1
          sendbuf (nzc,nentry) = var12(i,j,k)
        ENDDO
      ENDDO
    ENDDO
  ENDIF
  
  ! optional variables that are present
  IF (lpres13 .EQV. .TRUE.) THEN
    DO k = 1, kdim(13)
      DO j = jlo, jup
        DO i = ilo, iup
          nzc = nzc + 1
          sendbuf (nzc,nentry) = var13(i,j,k)
        ENDDO
      ENDDO
    ENDDO
  ENDIF
  
  ! optional variables that are present
  IF (lpres14 .EQV. .TRUE.) THEN
    DO k = 1, kdim(14)
      DO j = jlo, jup
        DO i = ilo, iup
          nzc = nzc + 1
          sendbuf (nzc,nentry) = var14(i,j,k)
        ENDDO
      ENDDO
    ENDDO
  ENDIF
  
  ! optional variables that are present
  IF (lpres15 .EQV. .TRUE.) THEN
    DO k = 1, kdim(15)
      DO j = jlo, jup
        DO i = ilo, iup
          nzc = nzc + 1
          sendbuf (nzc,nentry) = var15(i,j,k)
        ENDDO
      ENDDO
    ENDDO
  ENDIF
  
  ! optional variables that are present
  IF (lpres16 .EQV. .TRUE.) THEN
    DO k = 1, kdim(16)
      DO j = jlo, jup
        DO i = ilo, iup
          nzc = nzc + 1
          sendbuf (nzc,nentry) = var16(i,j,k)
        ENDDO
      ENDDO
    ENDDO
  ENDIF
  
  ! optional variables that are present
  IF (lpres17 .EQV. .TRUE.) THEN
    DO k = 1, kdim(17)
      DO j = jlo, jup
        DO i = ilo, iup
          nzc = nzc + 1
          sendbuf (nzc,nentry) = var17(i,j,k)
        ENDDO
      ENDDO
    ENDDO
  ENDIF
  
  ! optional variables that are present
  IF (lpres18 .EQV. .TRUE.) THEN
    DO k = 1, kdim(18)
      DO j = jlo, jup
        DO i = ilo, iup
          nzc = nzc + 1
          sendbuf (nzc,nentry) = var18(i,j,k)
        ENDDO
      ENDDO
    ENDDO
  ENDIF
  
  ! optional variables that are present
  IF (lpres19 .EQV. .TRUE.) THEN
    DO k = 1, kdim(19)
      DO j = jlo, jup
        DO i = ilo, iup
          nzc = nzc + 1
          sendbuf (nzc,nentry) = var19(i,j,k)
        ENDDO
      ENDDO
    ENDDO
  ENDIF
  
  ! optional variables that are present
  IF (lpres20 .EQV. .TRUE.) THEN
    DO k = 1, kdim(20)
      DO j = jlo, jup
        DO i = ilo, iup
          nzc = nzc + 1
          sendbuf (nzc,nentry) = var20(i,j,k)
        ENDDO
      ENDDO
    ENDDO
  ENDIF
  
  ! optional variables that are present
  IF (lpres21 .EQV. .TRUE.) THEN
    DO k = 1, kdim(21)
      DO j = jlo, jup
        DO i = ilo, iup
          nzc = nzc + 1
          sendbuf (nzc,nentry) = var21(i,j,k)
        ENDDO
      ENDDO
    ENDDO
  ENDIF
  
  ! optional variables that are present
  IF (lpres22 .EQV. .TRUE.) THEN
    DO k = 1, kdim(22)
      DO j = jlo, jup
        DO i = ilo, iup
          nzc = nzc + 1
          sendbuf (nzc,nentry) = var22(i,j,k)
        ENDDO
      ENDDO
    ENDDO
  ENDIF
  
  ! optional variables that are present
  IF (lpres23 .EQV. .TRUE.) THEN
    DO k = 1, kdim(23)
      DO j = jlo, jup
        DO i = ilo, iup
          nzc = nzc + 1
          sendbuf (nzc,nentry) = var23(i,j,k)
        ENDDO
      ENDDO
    ENDDO
  ENDIF
  
  ! optional variables that are present
  IF (lpres24 .EQV. .TRUE.) THEN
    DO k = 1, kdim(24)
      DO j = jlo, jup
        DO i = ilo, iup
          nzc = nzc + 1
          sendbuf (nzc,nentry) = var24(i,j,k)
        ENDDO
      ENDDO
    ENDDO
  ENDIF
  
  ! put nzc to global counter
  ncount = nzc

!------------------------------------------------------------------------------
!- End of the Subroutine
!------------------------------------------------------------------------------
  
END SUBROUTINE putbuf 

!==============================================================================
!==============================================================================
!+ This subroutine gets all necessary values from sendbuf
!------------------------------------------------------------------------------

SUBROUTINE getbuf  (var01, var02, var03, var04, var05, var06, var07, var08, &
                    var09, var10, var11, var12, var13, var14, var15, var16, &
                    var17, var18, var19, var20, var21, var22, var23, var24, &
                    sendbuf, isendbuflen, idim, jdim, kdim,                 &
                    ilo, iup, jlo, jup, ncount, nentry )

!------------------------------------------------------------------------------
!
! Description:
!   This subroutine gets the necessary values for the present variables
!   (determined by ilo, iup, jlo, jup) from sendbuf.
!
! Method:
!   Check which variables are present.
!
!------------------------------------------------------------------------------

! Subroutine arguments
  INTEGER (KIND=iintegers), INTENT (IN)         ::    &
    isendbuflen,                  & ! length of sendbuffer
    idim, jdim, kdim(24),         & ! dimensions of the fields
    ilo, iup, jlo, jup,           & ! start- and end-indices
    nentry                          ! specifies the row of sendbuf to be used

  INTEGER (KIND=iintegers), INTENT (INOUT)      ::    &
    ncount                          ! counts the variables

  REAL (KIND=wp),     INTENT (INOUT)            ::    &
    sendbuf (isendbuflen, 8),     & ! send buffer
    var01  (idim, jdim, kdim( 1))   ! first field that has to occur

  REAL (KIND=wp),     OPTIONAL, INTENT (INOUT)  ::    &
    var02  (idim, jdim, kdim( 2)),& ! additional optional fields
    var03  (idim, jdim, kdim( 3)),& ! additional optional fields
    var04  (idim, jdim, kdim( 4)),& ! additional optional fields
    var05  (idim, jdim, kdim( 5)),& ! additional optional fields
    var06  (idim, jdim, kdim( 6)),& ! additional optional fields
    var07  (idim, jdim, kdim( 7)),& ! additional optional fields
    var08  (idim, jdim, kdim( 8)),& ! additional optional fields
    var09  (idim, jdim, kdim( 9)),& ! additional optional fields
    var10  (idim, jdim, kdim(10)),& ! additional optional fields
    var11  (idim, jdim, kdim(11)),& ! additional optional fields
    var12  (idim, jdim, kdim(12)),& ! additional optional fields
    var13  (idim, jdim, kdim(13)),& ! additional optional fields
    var14  (idim, jdim, kdim(14)),& ! additional optional fields
    var15  (idim, jdim, kdim(15)),& ! additional optional fields
    var16  (idim, jdim, kdim(16)),& ! additional optional fields
    var17  (idim, jdim, kdim(17)),& ! additional optional fields
    var18  (idim, jdim, kdim(18)),& ! additional optional fields
    var19  (idim, jdim, kdim(19)),& ! additional optional fields
    var20  (idim, jdim, kdim(20)),& ! additional optional fields
    var21  (idim, jdim, kdim(21)),& ! additional optional fields
    var22  (idim, jdim, kdim(22)),& ! additional optional fields
    var23  (idim, jdim, kdim(23)),& ! additional optional fields
    var24  (idim, jdim, kdim(24))   ! additional optional fields

! Local variables

  INTEGER (KIND=iintegers)   ::       &
    i, j, k, nzc

  LOGICAL                    ::       &
    lpres02, lpres03, lpres04, lpres05, lpres06, lpres07, lpres08, lpres09, &
    lpres10, lpres11, lpres12, lpres13, lpres14, lpres15, lpres16, lpres17, &
    lpres18, lpres19, lpres20, lpres21, lpres22, lpres23, lpres24

!------------------------------------------------------------------------------
!- End of header -
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!- Section 1: Initializations
!------------------------------------------------------------------------------

  ! check which variables are present
  lpres02 = PRESENT (var02)
  lpres03 = PRESENT (var03)
  lpres04 = PRESENT (var04)
  lpres05 = PRESENT (var05)
  lpres06 = PRESENT (var06)
  lpres07 = PRESENT (var07)
  lpres08 = PRESENT (var08)
  lpres09 = PRESENT (var09)
  lpres10 = PRESENT (var10)
  lpres11 = PRESENT (var11)
  lpres12 = PRESENT (var12)
  lpres13 = PRESENT (var13)
  lpres14 = PRESENT (var14)
  lpres15 = PRESENT (var15)
  lpres16 = PRESENT (var16)
  lpres17 = PRESENT (var17)
  lpres18 = PRESENT (var18)
  lpres19 = PRESENT (var19)
  lpres20 = PRESENT (var20)
  lpres21 = PRESENT (var21)
  lpres22 = PRESENT (var22)
  lpres23 = PRESENT (var23)
  lpres24 = PRESENT (var24)

!------------------------------------------------------------------------------
!- Section 2: Get data from the buffer
!------------------------------------------------------------------------------

  ! use nzc as a local counter  (based on a work from Mike O'Neill to 
  ! improve vectorization of putbuf and getbuf)
  nzc = ncount

  ! first variable that has to be present
  DO k = 1, kdim( 1)
    DO j = jlo, jup
      DO i = ilo, iup
        nzc = nzc + 1
        var01(i,j,k) = sendbuf (nzc,nentry)
      ENDDO
    ENDDO
  ENDDO

  ! optional variables that are present
  IF (lpres02.EQV. .TRUE.) THEN
    DO k = 1, kdim( 2)
      DO j = jlo, jup
        DO i = ilo, iup
          nzc = nzc + 1
          var02(i,j,k) = sendbuf (nzc,nentry)
        ENDDO
      ENDDO
    ENDDO
  ENDIF
  
  ! optional variables that are present
  IF (lpres03.EQV. .TRUE.) THEN
    DO k = 1, kdim( 3)
      DO j = jlo, jup
        DO i = ilo, iup
          nzc = nzc + 1
          var03(i,j,k) = sendbuf (nzc,nentry)
        ENDDO
      ENDDO
    ENDDO
  ENDIF
  
  ! optional variables that are present
  IF (lpres04.EQV. .TRUE.) THEN
    DO k = 1, kdim( 4)
      DO j = jlo, jup
        DO i = ilo, iup
          nzc = nzc + 1
          var04(i,j,k) = sendbuf (nzc,nentry)
        ENDDO
      ENDDO
    ENDDO
  ENDIF
  
  ! optional variables that are present
  IF (lpres05.EQV. .TRUE.) THEN
    DO k = 1, kdim( 5)
      DO j = jlo, jup
        DO i = ilo, iup
          nzc = nzc + 1
          var05(i,j,k) = sendbuf (nzc,nentry)
        ENDDO
      ENDDO
    ENDDO
  ENDIF
  
  ! optional variables that are present
  IF (lpres06.EQV. .TRUE.) THEN
    DO k = 1, kdim( 6)
      DO j = jlo, jup
        DO i = ilo, iup
          nzc = nzc + 1
          var06(i,j,k) = sendbuf (nzc,nentry)
        ENDDO
      ENDDO
    ENDDO
  ENDIF
  
  ! optional variables that are present
  IF (lpres07.EQV. .TRUE.) THEN
    DO k = 1, kdim( 7)
      DO j = jlo, jup
        DO i = ilo, iup
          nzc = nzc + 1
          var07(i,j,k) = sendbuf (nzc,nentry)
        ENDDO
      ENDDO
    ENDDO
  ENDIF
  
  ! optional variables that are present
  IF (lpres08.EQV. .TRUE.) THEN
    DO k = 1, kdim( 8)
      DO j = jlo, jup
        DO i = ilo, iup
          nzc = nzc + 1
          var08(i,j,k) = sendbuf (nzc,nentry)
        ENDDO
      ENDDO
    ENDDO
  ENDIF
  
  ! optional variables that are present
  IF (lpres09 .EQV. .TRUE.) THEN
    DO k = 1, kdim( 9)
      DO j = jlo, jup
        DO i = ilo, iup
          nzc = nzc + 1
          var09(i,j,k) = sendbuf (nzc,nentry)
        ENDDO
      ENDDO
    ENDDO
  ENDIF
  
  ! optional variables that are present
  IF (lpres10 .EQV. .TRUE.) THEN
    DO k = 1, kdim(10)
      DO j = jlo, jup
        DO i = ilo, iup
          nzc = nzc + 1
          var10(i,j,k) = sendbuf (nzc,nentry)
        ENDDO
      ENDDO
    ENDDO
  ENDIF
  
  ! optional variables that are present
  IF (lpres11 .EQV. .TRUE.) THEN
    DO k = 1, kdim(11)
      DO j = jlo, jup
        DO i = ilo, iup
          nzc = nzc + 1
          var11(i,j,k) = sendbuf (nzc,nentry)
        ENDDO
      ENDDO
    ENDDO
  ENDIF
  
  ! optional variables that are present
  IF (lpres12 .EQV. .TRUE.) THEN
    DO k = 1, kdim(12)
      DO j = jlo, jup
        DO i = ilo, iup
          nzc = nzc + 1
          var12(i,j,k) = sendbuf (nzc,nentry)
        ENDDO
      ENDDO
    ENDDO
  ENDIF
  
  ! optional variables that are present
  IF (lpres13 .EQV. .TRUE.) THEN
    DO k = 1, kdim(13)
      DO j = jlo, jup
        DO i = ilo, iup
          nzc = nzc + 1
          var13(i,j,k) = sendbuf (nzc,nentry)
        ENDDO
      ENDDO
    ENDDO
  ENDIF
  
  ! optional variables that are present
  IF (lpres14 .EQV. .TRUE.) THEN
    DO k = 1, kdim(14)
      DO j = jlo, jup
        DO i = ilo, iup
          nzc = nzc + 1
          var14(i,j,k) = sendbuf (nzc,nentry)
        ENDDO
      ENDDO
    ENDDO
  ENDIF
  
  ! optional variables that are present
  IF (lpres15 .EQV. .TRUE.) THEN
    DO k = 1, kdim(15)
      DO j = jlo, jup
        DO i = ilo, iup
          nzc = nzc + 1
          var15(i,j,k) = sendbuf (nzc,nentry)
        ENDDO
      ENDDO
    ENDDO
  ENDIF
  
  ! optional variables that are present
  IF (lpres16 .EQV. .TRUE.) THEN
    DO k = 1, kdim(16)
      DO j = jlo, jup
        DO i = ilo, iup
          nzc = nzc + 1
          var16(i,j,k) = sendbuf (nzc,nentry)
        ENDDO
      ENDDO
    ENDDO
  ENDIF
  
  ! optional variables that are present
  IF (lpres17 .EQV. .TRUE.) THEN
    DO k = 1, kdim(17)
      DO j = jlo, jup
        DO i = ilo, iup
          nzc = nzc + 1
          var17(i,j,k) = sendbuf (nzc,nentry)
        ENDDO
      ENDDO
    ENDDO
  ENDIF
  
  ! optional variables that are present
  IF (lpres18 .EQV. .TRUE.) THEN
    DO k = 1, kdim(18)
      DO j = jlo, jup
        DO i = ilo, iup
          nzc = nzc + 1
          var18(i,j,k) = sendbuf (nzc,nentry)
        ENDDO
      ENDDO
    ENDDO
  ENDIF
  
  ! optional variables that are present
  IF (lpres19 .EQV. .TRUE.) THEN
    DO k = 1, kdim(19)
      DO j = jlo, jup
        DO i = ilo, iup
          nzc = nzc + 1
          var19(i,j,k) = sendbuf (nzc,nentry)
        ENDDO
      ENDDO
    ENDDO
  ENDIF
  
  ! optional variables that are present
  IF (lpres20 .EQV. .TRUE.) THEN
    DO k = 1, kdim(20)
      DO j = jlo, jup
        DO i = ilo, iup
          nzc = nzc + 1
          var20(i,j,k) = sendbuf (nzc,nentry)
        ENDDO
      ENDDO
    ENDDO
  ENDIF
  
  ! optional variables that are present
  IF (lpres21 .EQV. .TRUE.) THEN
    DO k = 1, kdim(21)
      DO j = jlo, jup
        DO i = ilo, iup
          nzc = nzc + 1
          var21(i,j,k) = sendbuf (nzc,nentry)
        ENDDO
      ENDDO
    ENDDO
  ENDIF
  
  ! optional variables that are present
  IF (lpres22 .EQV. .TRUE.) THEN
    DO k = 1, kdim(22)
      DO j = jlo, jup
        DO i = ilo, iup
          nzc = nzc + 1
          var22(i,j,k) = sendbuf (nzc,nentry)
        ENDDO
      ENDDO
    ENDDO
  ENDIF
  
  ! optional variables that are present
  IF (lpres23 .EQV. .TRUE.) THEN
    DO k = 1, kdim(23)
      DO j = jlo, jup
        DO i = ilo, iup
          nzc = nzc + 1
          var23(i,j,k) = sendbuf (nzc,nentry)
        ENDDO
      ENDDO
    ENDDO
  ENDIF
  
  ! optional variables that are present
  IF (lpres24 .EQV. .TRUE.) THEN
    DO k = 1, kdim(24)
      DO j = jlo, jup
        DO i = ilo, iup
          nzc = nzc + 1
          var24(i,j,k) = sendbuf (nzc,nentry)
        ENDDO
      ENDDO
    ENDDO
  ENDIF
  
  ! put nzc to global counter
  ncount = nzc

!------------------------------------------------------------------------------
!- End of the Subroutine
!------------------------------------------------------------------------------
  
END SUBROUTINE getbuf 

!==============================================================================
!==============================================================================
!+ Calls MPI_BARRIER
!------------------------------------------------------------------------------

SUBROUTINE comm_barrier (icomm, ierror, yerrmsg)

!------------------------------------------------------------------------------
!
! Description:
!
! Method:
!   MPI-routine MPI_BARRIER
!
!------------------------------------------------------------------------------

! Subroutine arguments
  INTEGER (KIND=iintegers), INTENT (IN)      ::       &
    icomm            ! communicator to be used

  INTEGER (KIND=iintegers), INTENT (OUT)     ::       &
    ierror           ! error-status variable

  CHARACTER (LEN=*),        INTENT (OUT)     ::       &
    yerrmsg            ! for MPI error message

! Local variables

  INTEGER (KIND=iintegers)   ::       &
    izmplcode                           ! for MPI error code

!------------------------------------------------------------------------------
!- End of header -
!------------------------------------------------------------------------------

ierror    = 0
yerrmsg   = '   '
izmplcode = 0

CALL MPI_BARRIER (icomm, izmplcode)

IF (izmplcode /= 0) THEN
  ierror = izmplcode
  yerrmsg = 'MPI_BARRIER'
ENDIF

!------------------------------------------------------------------------------
!- End of the Subroutine
!------------------------------------------------------------------------------

END SUBROUTINE comm_barrier

!==============================================================================
!==============================================================================
!+ This subroutine defines and allocates MPI data types for exchg_boundaries
!------------------------------------------------------------------------------

SUBROUTINE setup_data_type                                                   &
     ( var01, var02, var03, var04, var05, var06, var07, var08, var09, var10, &
       var11, var12, var13, var14, var15, var16, var17, var18, var19, var20, &
       var21, var22, var23, var24,                                           &
       idim, jdim, kdim, ilo, iup, jlo, jup,                                 &
       ierror, yerrmsg, imp_type, ncount, type_handle)

!------------------------------------------------------------------------------
!
! Description:
!   This subroutine allocates and commits a data type consisting of a
!   certain subgrid of up to 14 variables. The subgrid is specified
!   via  idim, jdim, kdim, ilo, iup, jlo, jup, klo, kup . Only
!   one variable has to occur, the other are optional.
!
!   The values of "var1(ilo, jlo, klo), ncount, datatype" 
!   should constitute the right entries for the dummy variables 
!   "BUF, COUNT, DATATYPE" in a call to
!   MPI_SEND/MPI_RECV. This should describe the subgrid for all
!   up to 14 variables.
!
!   As a consequence, the same data type can only be used, if the same
!   variables (same position in memory !) are used. 
!
!  Author: C. Pospiech, IBM
!
!------------------------------------------------------------------------------

! Subroutine arguments
  INTEGER (KIND=iintegers), INTENT (IN)         ::    &
    imp_type,                     & ! determines the REAL type used
    idim, jdim,                   & ! horizontal dimensions of the fields
    kdim(24),                     & ! vertical dimensions of var01..var20
    ilo, iup, jlo, jup              ! start- and end-indices

  INTEGER (KIND=iintegers), INTENT (OUT)      ::    &
    ncount,                       & ! how many copies of type_handle
    type_handle                     ! handle for MPI data type

  INTEGER (KIND=iintegers), INTENT (OUT)        ::    &
    ierror                ! error status variable

  CHARACTER (LEN=*),        INTENT(OUT)  ::       &
    yerrmsg               ! for MPI error message

  REAL (KIND=wp),     INTENT (INOUT)            ::    &
    var01 (idim,jdim,kdim( 1))         ! first field that has to occur

  REAL (KIND=wp),     OPTIONAL, INTENT (INOUT)  ::    &
    var02 (idim,jdim,kdim( 2)),& ! additional optional fields
    var03 (idim,jdim,kdim( 3)),& ! additional optional fields
    var04 (idim,jdim,kdim( 4)),& ! additional optional fields
    var05 (idim,jdim,kdim( 5)),& ! additional optional fields
    var06 (idim,jdim,kdim( 6)),& ! additional optional fields
    var07 (idim,jdim,kdim( 7)),& ! additional optional fields
    var08 (idim,jdim,kdim( 8)),& ! additional optional fields
    var09 (idim,jdim,kdim( 9)),& ! additional optional fields
    var10 (idim,jdim,kdim(10)),& ! additional optional fields
    var11 (idim,jdim,kdim(11)),& ! additional optional fields
    var12 (idim,jdim,kdim(12)),& ! additional optional fields
    var13 (idim,jdim,kdim(13)),& ! additional optional fields
    var14 (idim,jdim,kdim(14)),& ! additional optional fields
    var15 (idim,jdim,kdim(15)),& ! additional optional fields
    var16 (idim,jdim,kdim(16)),& ! additional optional fields
    var17 (idim,jdim,kdim(17)),& ! additional optional fields
    var18 (idim,jdim,kdim(18)),& ! additional optional fields
    var19 (idim,jdim,kdim(19)),& ! additional optional fields
    var20 (idim,jdim,kdim(20)),& ! additional optional fields
    var21 (idim,jdim,kdim(21)),& ! additional optional fields
    var22 (idim,jdim,kdim(22)),& ! additional optional fields
    var23 (idim,jdim,kdim(23)),& ! additional optional fields
    var24 (idim,jdim,kdim(24))   ! additional optional fields

! Local variables

  INTEGER (KIND=MPI_ADDRESS_KIND) ::  &
    meta_addr(24),         &   ! Vector of addresses of the varxx
    meta_disp(24)              ! Displacements of varxx in memory

  INTEGER (KIND=iintegers)   ::       &
    nzc,                   &
    sect1d, sect2d,        &   ! Variables to hold intermediate 
    sect3d(24), meta_vect, &   ! MPI data types
    meta_blklen(24),       &   ! some intermediate variable
    num_meta_entries,      &   ! how many varxx are present
    disp(2), blocklen(2),  &   ! Variables needed to define
    vartype(2),            &   ! MPI data type of a certain extent
    sizeofreal                 ! size of data type in byte

  INTEGER (KIND=iintegers)   ::       &
    izmplcode                   ! for MPI error code

  LOGICAL                    ::       &
    lpres02, lpres03, lpres04, lpres05, lpres06, lpres07, lpres08, lpres09,  &
    lpres10, lpres11, lpres12, lpres13, lpres14, lpres15, lpres16, lpres17,  &
    lpres18, lpres19, lpres20, lpres21, lpres22, lpres23, lpres24

!------------------------------------------------------------------------------
!- End of header -
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!- Section 1: Initializations
!------------------------------------------------------------------------------

  ! check which variables are present
  lpres02  = PRESENT (var02)
  lpres03  = PRESENT (var03)
  lpres04  = PRESENT (var04)
  lpres05  = PRESENT (var05)
  lpres06  = PRESENT (var06)
  lpres07  = PRESENT (var07)
  lpres08  = PRESENT (var08)
  lpres09  = PRESENT (var09)
  lpres10  = PRESENT (var10)
  lpres11  = PRESENT (var11)
  lpres12  = PRESENT (var12)
  lpres13  = PRESENT (var13)
  lpres14  = PRESENT (var14)
  lpres15  = PRESENT (var15)
  lpres16  = PRESENT (var16)
  lpres17  = PRESENT (var17)
  lpres18  = PRESENT (var18)
  lpres19  = PRESENT (var19)
  lpres20  = PRESENT (var20)
  lpres21  = PRESENT (var21)
  lpres22  = PRESENT (var22)
  lpres23  = PRESENT (var23)
  lpres24  = PRESENT (var24)
  sect3d(:) = MPI_DATATYPE_NULL

!------------------------------------------------------------------------------
!- Section 2: Set up of MPI data types *** subarrays
!------------------------------------------------------------------------------

  ! set up 1-dimensional section
  nzc = iup - ilo + 1
  CALL MPI_TYPE_CONTIGUOUS  (nzc, imp_type, sect1d, izmplcode)
  IF (izmplcode /= 0) THEN
    ierror  = izmplcode
    yerrmsg = 'MPI_TYPE_CONTIGUOUS'
    RETURN
  ENDIF

  ! set up 2-dimensional section
  nzc = jup - jlo +1
  CALL MPI_TYPE_EXTENT(imp_type, sizeofreal, izmplcode)
  IF (izmplcode /= 0) THEN
    ierror  = izmplcode
    yerrmsg = 'MPI_TYPE_EXTENT'
    RETURN
  ENDIF
  CALL MPI_TYPE_HVECTOR    (nzc, 1, idim*sizeofreal,         &
                            sect1d, sect2d, izmplcode)
  IF (izmplcode /= 0) THEN
    ierror  = izmplcode
    yerrmsg = 'MPI_TYPE_HVECTOR-2'
    RETURN
  ENDIF 

!US: this must be done for every optional entry
  ! set up 3-dimensional section
  CALL MPI_TYPE_HVECTOR      (kdim( 1), 1, idim*jdim*sizeofreal,    &
                              sect2d, sect3d( 1), izmplcode)
  IF (lpres02) THEN
    CALL MPI_TYPE_HVECTOR    (kdim( 2), 1, idim*jdim*sizeofreal,    &
                              sect2d, sect3d( 2), izmplcode)
  ENDIF
  IF (lpres03) THEN
    CALL MPI_TYPE_HVECTOR    (kdim( 3), 1, idim*jdim*sizeofreal,    &
                              sect2d, sect3d( 3), izmplcode)
  ENDIF
  IF (lpres04) THEN
    CALL MPI_TYPE_HVECTOR    (kdim( 4), 1, idim*jdim*sizeofreal,    &
                              sect2d, sect3d( 4), izmplcode)
  ENDIF
  IF (lpres05) THEN
    CALL MPI_TYPE_HVECTOR    (kdim( 5), 1, idim*jdim*sizeofreal,    &
                              sect2d, sect3d( 5), izmplcode)
  ENDIF
  IF (lpres06) THEN
    CALL MPI_TYPE_HVECTOR    (kdim( 6), 1, idim*jdim*sizeofreal,    &
                              sect2d, sect3d( 6), izmplcode)
  ENDIF
  IF (lpres07) THEN
    CALL MPI_TYPE_HVECTOR    (kdim( 7), 1, idim*jdim*sizeofreal,    &
                              sect2d, sect3d( 7), izmplcode)
  ENDIF
  IF (lpres08) THEN
    CALL MPI_TYPE_HVECTOR    (kdim( 8), 1, idim*jdim*sizeofreal,    &
                              sect2d, sect3d( 8), izmplcode)
  ENDIF
  IF (lpres09) THEN
    CALL MPI_TYPE_HVECTOR    (kdim( 9), 1, idim*jdim*sizeofreal,    &
                              sect2d, sect3d( 9), izmplcode)
  ENDIF
  IF (lpres10) THEN
    CALL MPI_TYPE_HVECTOR    (kdim(10), 1, idim*jdim*sizeofreal,    &
                              sect2d, sect3d(10), izmplcode)
  ENDIF
  IF (lpres11) THEN
    CALL MPI_TYPE_HVECTOR    (kdim(11), 1, idim*jdim*sizeofreal,    &
                              sect2d, sect3d(11), izmplcode)
  ENDIF
  IF (lpres12) THEN
    CALL MPI_TYPE_HVECTOR    (kdim(12), 1, idim*jdim*sizeofreal,    &
                              sect2d, sect3d(12), izmplcode)
  ENDIF
  IF (lpres13) THEN
    CALL MPI_TYPE_HVECTOR    (kdim(13), 1, idim*jdim*sizeofreal,    &
                              sect2d, sect3d(13), izmplcode)
  ENDIF
  IF (lpres14) THEN
    CALL MPI_TYPE_HVECTOR    (kdim(14), 1, idim*jdim*sizeofreal,    &
                              sect2d, sect3d(14), izmplcode)
  ENDIF
  IF (lpres15) THEN
    CALL MPI_TYPE_HVECTOR    (kdim(15), 1, idim*jdim*sizeofreal,    &
                              sect2d, sect3d(15), izmplcode)
  ENDIF
  IF (lpres16) THEN
    CALL MPI_TYPE_HVECTOR    (kdim(16), 1, idim*jdim*sizeofreal,    &
                              sect2d, sect3d(16), izmplcode)
  ENDIF
  IF (lpres17) THEN
    CALL MPI_TYPE_HVECTOR    (kdim(17), 1, idim*jdim*sizeofreal,    &
                              sect2d, sect3d(17), izmplcode)
  ENDIF
  IF (lpres18) THEN
    CALL MPI_TYPE_HVECTOR    (kdim(18), 1, idim*jdim*sizeofreal,    &
                              sect2d, sect3d(18), izmplcode)
  ENDIF
  IF (lpres19) THEN
    CALL MPI_TYPE_HVECTOR    (kdim(19), 1, idim*jdim*sizeofreal,    &
                              sect2d, sect3d(19), izmplcode)
  ENDIF
  IF (lpres20) THEN
    CALL MPI_TYPE_HVECTOR    (kdim(20), 1, idim*jdim*sizeofreal,    &
                              sect2d, sect3d(20), izmplcode)
  ENDIF
  IF (lpres21) THEN
    CALL MPI_TYPE_HVECTOR    (kdim(21), 1, idim*jdim*sizeofreal,    &
                              sect2d, sect3d(21), izmplcode)
  ENDIF
  IF (lpres22) THEN
    CALL MPI_TYPE_HVECTOR    (kdim(22), 1, idim*jdim*sizeofreal,    &
                              sect2d, sect3d(22), izmplcode)
  ENDIF
  IF (lpres23) THEN
    CALL MPI_TYPE_HVECTOR    (kdim(23), 1, idim*jdim*sizeofreal,    &
                              sect2d, sect3d(23), izmplcode)
  ENDIF
  IF (lpres24) THEN
    CALL MPI_TYPE_HVECTOR    (kdim(24), 1, idim*jdim*sizeofreal,    &
                              sect2d, sect3d(24), izmplcode)
  ENDIF

  IF (izmplcode /= 0) THEN
    ierror  = izmplcode
    yerrmsg = 'MPI_TYPE_HVECTOR-3'
    RETURN
  ENDIF 

!------------------------------------------------------------------------------
!- Section 3: Set up of MPI data types *** meta structure from all varxx
!------------------------------------------------------------------------------

  num_meta_entries = 1
  CALL MPI_GET_ADDRESS (var01, meta_addr(num_meta_entries), izmplcode)
  IF (izmplcode /= 0) THEN
    ierror  = izmplcode
    yerrmsg = 'MPI_GET_ADDRESS-01'
    RETURN
  ENDIF

  IF ( lpres02 .EQV. .TRUE.) THEN
     num_meta_entries = num_meta_entries + 1
     CALL MPI_GET_ADDRESS (var02, meta_addr(num_meta_entries), izmplcode)
     IF (izmplcode /= 0) THEN
       ierror  = izmplcode
       yerrmsg = 'MPI_GET_ADDRESS-02'
       RETURN
     ENDIF
  ENDIF

  IF ( lpres03 .EQV. .TRUE.) THEN
     num_meta_entries = num_meta_entries + 1
     CALL MPI_GET_ADDRESS (var03, meta_addr(num_meta_entries), izmplcode)
     IF (izmplcode /= 0) THEN
       ierror  = izmplcode
       yerrmsg = 'MPI_GET_ADDRESS-03'
       RETURN
     ENDIF
  ENDIF

  IF ( lpres04 .EQV. .TRUE.) THEN
     num_meta_entries = num_meta_entries + 1
     CALL MPI_GET_ADDRESS (var04, meta_addr(num_meta_entries), izmplcode)
     IF (izmplcode /= 0) THEN
       ierror  = izmplcode
       yerrmsg = 'MPI_GET_ADDRESS-04'
       RETURN
     ENDIF
  ENDIF

  IF ( lpres05 .EQV. .TRUE.) THEN
     num_meta_entries = num_meta_entries + 1
     CALL MPI_GET_ADDRESS (var05, meta_addr(num_meta_entries), izmplcode)
     IF (izmplcode /= 0) THEN
       ierror  = izmplcode
       yerrmsg = 'MPI_GET_ADDRESS-05'
       RETURN
     ENDIF
  ENDIF

  IF ( lpres06 .EQV. .TRUE.) THEN
     num_meta_entries = num_meta_entries + 1
     CALL MPI_GET_ADDRESS (var06, meta_addr(num_meta_entries), izmplcode)
     IF (izmplcode /= 0) THEN
       ierror  = izmplcode
       yerrmsg = 'MPI_GET_ADDRESS-06'
       RETURN
     ENDIF
  ENDIF

  IF ( lpres07 .EQV. .TRUE.) THEN
     num_meta_entries = num_meta_entries + 1
     CALL MPI_GET_ADDRESS (var07, meta_addr(num_meta_entries), izmplcode)
     IF (izmplcode /= 0) THEN
       ierror  = izmplcode
       yerrmsg = 'MPI_GET_ADDRESS-07'
       RETURN
     ENDIF
  ENDIF

  IF ( lpres08 .EQV. .TRUE.) THEN
     num_meta_entries = num_meta_entries + 1
     CALL MPI_GET_ADDRESS (var08, meta_addr(num_meta_entries), izmplcode)
     IF (izmplcode /= 0) THEN
       ierror  = izmplcode
       yerrmsg = 'MPI_GET_ADDRESS-08'
       RETURN
     ENDIF
  ENDIF

  IF ( lpres09 .EQV. .TRUE.) THEN
     num_meta_entries = num_meta_entries + 1
     CALL MPI_GET_ADDRESS (var09, meta_addr(num_meta_entries), izmplcode)
     IF (izmplcode /= 0) THEN
       ierror  = izmplcode
       yerrmsg = 'MPI_GET_ADDRESS-09'
       RETURN
     ENDIF
  ENDIF

  IF ( lpres10 .EQV. .TRUE.) THEN
     num_meta_entries = num_meta_entries + 1
     CALL MPI_GET_ADDRESS (var10, meta_addr(num_meta_entries), izmplcode)
     IF (izmplcode /= 0) THEN
       ierror  = izmplcode
       yerrmsg = 'MPI_GET_ADDRESS-10'
       RETURN
     ENDIF
  ENDIF

  IF ( lpres11 .EQV. .TRUE.) THEN
     num_meta_entries = num_meta_entries + 1
     CALL MPI_GET_ADDRESS (var11, meta_addr(num_meta_entries), izmplcode)
     IF (izmplcode /= 0) THEN
       ierror  = izmplcode
       yerrmsg = 'MPI_GET_ADDRESS-11'
       RETURN
     ENDIF
  ENDIF

  IF ( lpres12 .EQV. .TRUE.) THEN
     num_meta_entries = num_meta_entries + 1
     CALL MPI_GET_ADDRESS (var12, meta_addr(num_meta_entries), izmplcode)
     IF (izmplcode /= 0) THEN
       ierror  = izmplcode
       yerrmsg = 'MPI_GET_ADDRESS-12'
       RETURN
     ENDIF
  ENDIF

  IF ( lpres13 .EQV. .TRUE.) THEN
     num_meta_entries = num_meta_entries + 1
     CALL MPI_GET_ADDRESS (var13, meta_addr(num_meta_entries), izmplcode)
     IF (izmplcode /= 0) THEN
       ierror  = izmplcode
       yerrmsg = 'MPI_GET_ADDRESS-13'
       RETURN
     ENDIF
  ENDIF

  IF ( lpres14 .EQV. .TRUE.) THEN
     num_meta_entries = num_meta_entries + 1
     CALL MPI_GET_ADDRESS (var14, meta_addr(num_meta_entries), izmplcode)
     IF (izmplcode /= 0) THEN
       ierror  = izmplcode
       yerrmsg = 'MPI_GET_ADDRESS-14'
       RETURN
     ENDIF
  ENDIF

  IF ( lpres15 .EQV. .TRUE.) THEN
     num_meta_entries = num_meta_entries + 1
     CALL MPI_GET_ADDRESS (var15, meta_addr(num_meta_entries), izmplcode)
     IF (izmplcode /= 0) THEN
       ierror  = izmplcode
       yerrmsg = 'MPI_GET_ADDRESS-15'
       RETURN
     ENDIF
  ENDIF

  IF ( lpres16 .EQV. .TRUE.) THEN
     num_meta_entries = num_meta_entries + 1
     CALL MPI_GET_ADDRESS (var16, meta_addr(num_meta_entries), izmplcode)
     IF (izmplcode /= 0) THEN
       ierror  = izmplcode
       yerrmsg = 'MPI_GET_ADDRESS-16'
       RETURN
     ENDIF
  ENDIF

  IF ( lpres17 .EQV. .TRUE.) THEN
     num_meta_entries = num_meta_entries + 1
     CALL MPI_GET_ADDRESS (var17, meta_addr(num_meta_entries), izmplcode)
     IF (izmplcode /= 0) THEN
       ierror  = izmplcode
       yerrmsg = 'MPI_GET_ADDRESS-17'
       RETURN
     ENDIF
  ENDIF

  IF ( lpres18 .EQV. .TRUE.) THEN
     num_meta_entries = num_meta_entries + 1
     CALL MPI_GET_ADDRESS (var18, meta_addr(num_meta_entries), izmplcode)
     IF (izmplcode /= 0) THEN
       ierror  = izmplcode
       yerrmsg = 'MPI_GET_ADDRESS-18'
       RETURN
     ENDIF
  ENDIF

  IF ( lpres19 .EQV. .TRUE.) THEN
     num_meta_entries = num_meta_entries + 1
     CALL MPI_GET_ADDRESS (var19, meta_addr(num_meta_entries), izmplcode)
     IF (izmplcode /= 0) THEN
       ierror  = izmplcode
       yerrmsg = 'MPI_GET_ADDRESS-19'
       RETURN
     ENDIF
  ENDIF

  IF ( lpres20 .EQV. .TRUE.) THEN
     num_meta_entries = num_meta_entries + 1
     CALL MPI_GET_ADDRESS (var20, meta_addr(num_meta_entries), izmplcode)
     IF (izmplcode /= 0) THEN
       ierror  = izmplcode
       yerrmsg = 'MPI_GET_ADDRESS-20'
       RETURN
     ENDIF
  ENDIF

  IF ( lpres21 .EQV. .TRUE.) THEN
     num_meta_entries = num_meta_entries + 1
     CALL MPI_GET_ADDRESS (var21, meta_addr(num_meta_entries), izmplcode)
     IF (izmplcode /= 0) THEN
       ierror  = izmplcode
       yerrmsg = 'MPI_GET_ADDRESS-21'
       RETURN
     ENDIF
  ENDIF

  IF ( lpres22 .EQV. .TRUE.) THEN
     num_meta_entries = num_meta_entries + 1
     CALL MPI_GET_ADDRESS (var22, meta_addr(num_meta_entries), izmplcode)
     IF (izmplcode /= 0) THEN
       ierror  = izmplcode
       yerrmsg = 'MPI_GET_ADDRESS-22'
       RETURN
     ENDIF
  ENDIF

  IF ( lpres23 .EQV. .TRUE.) THEN
     num_meta_entries = num_meta_entries + 1
     CALL MPI_GET_ADDRESS (var23, meta_addr(num_meta_entries), izmplcode)
     IF (izmplcode /= 0) THEN
       ierror  = izmplcode
       yerrmsg = 'MPI_GET_ADDRESS-23'
       RETURN
     ENDIF
  ENDIF

  IF ( lpres24 .EQV. .TRUE.) THEN
     num_meta_entries = num_meta_entries + 1
     CALL MPI_GET_ADDRESS (var24, meta_addr(num_meta_entries), izmplcode)
     IF (izmplcode /= 0) THEN
       ierror  = izmplcode
       yerrmsg = 'MPI_GET_ADDRESS-24'
       RETURN
     ENDIF
  ENDIF

  meta_disp(:)   = meta_addr(:) - meta_addr(1)
  meta_blklen(:) = 1
  CALL MPI_TYPE_CREATE_STRUCT  (num_meta_entries, meta_blklen, meta_disp,   &
                                sect3d, meta_vect, izmplcode)
  IF (izmplcode /= 0) THEN
    ierror  = izmplcode
    yerrmsg = 'MPI_TYPE_CREATE_STRUCT'
    RETURN
  ENDIF
 
!------------------------------------------------------------------------------
!- Section 4: Reset extent of this new data type by defining upper bound
!------------------------------------------------------------------------------

  blocklen(:)   = 1
  disp(1)     = 0
  disp(2)     = sizeofreal*(iup - ilo + 2)
  vartype(1)  = meta_vect
  vartype(2)  = MPI_UB
  CALL MPI_TYPE_STRUCT     (2, blocklen, disp, vartype, type_handle, izmplcode)
  IF (izmplcode /= 0) THEN
    ierror  = izmplcode
    yerrmsg = 'MPI_TYPE_STRUCT'
    RETURN
  ENDIF

!------------------------------------------------------------------------------
!- Section 5: Commit the data type
!------------------------------------------------------------------------------

  CALL MPI_TYPE_COMMIT       (type_handle,izmplcode)
  IF (izmplcode /= 0) THEN
    ierror  = izmplcode
    yerrmsg = 'MPI_TYPE_COMMIT'
    RETURN
  ENDIF
  ncount = 1

!------------------------------------------------------------------------------
!- End of the Subroutine
!------------------------------------------------------------------------------

END SUBROUTINE setup_data_type

!==============================================================================

END MODULE environment
