!+ Source Code for reading cosmo binary restart files
!------------------------------------------------------------------------------

PROGRAM readCOSMObinr

!------------------------------------------------------------------------------
!
! Description:
!  This code is based on COSMO5.1 routines for reading COSMO binary restart files
!  Current Code Owner: P. Shrestha 
!    email: pshrestha@uni-bonn.de
!  History:
!  Version    Date       Name
!  ---------- ---------- ----
!  1.1        2016/06/03 P. Shrestha 
!  Initial release
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!==============================================================================
 

IMPLICIT NONE
 INTEGER, PARAMETER       :: intgribf  = KIND(1)
 INTEGER, PARAMETER       :: iintegers = KIND(1)
 INTEGER, PARAMETER       :: sp = SELECTED_REAL_KIND( 6, 37) 
 INTEGER, PARAMETER       :: dp = SELECTED_REAL_KIND(12,307) 
 INTEGER, PARAMETER       :: wp = dp

 INTEGER (KIND=intgribf),  PARAMETER      :: &
               npds  = 321_intgribf, & ! product definition section 
               ngds  = 626_intgribf    ! grid description section 

 INTEGER (KIND=iintegers), PARAMETER       ::  &
               idim  = 30, &
               jdim  = 20, &
               kdim  = 50
 INTEGER (KIND=iintegers) :: nudat, izerr,iz_countl 

 REAL(KIND=wp)            :: rbuf(idim*jdim)
 INTEGER(KIND=intgribf)   :: ipdsbuf(npds), igdsbuf(ngds)
 
 REAL(KIND=wp)            :: psm0, dsem0, msem0, kem0, qcm0
 INTEGER (KIND=iintegers) :: ntke
 INTEGER (KIND=iintegers) :: izvctype_read
 REAL(KIND=wp)            :: refatm_p0sl,   refatm_t0sl,   &
                             refatm_dt0lp, refatm_delta_t, &
                             refatm_h_scal, refatm_bvref,  &
                             vcoord_vcflat, zvc_params(kdim+1)
 CHARACTER (LEN=*), PARAMETER       ::  &
      datname = "/Users/pshrestha/work/ncl_current/DAtool/lrff00030000o"

!------------------------------------------------------------------------------
!Variable Defininiton Complete, program starts here

!------------------------------------------------------------------------------
!Section 1: INITIALIZE
!------------------------------------------------------------------------------
 
 iz_countl = 0
 nudat  = 100

!------------------------------------------------------------------------------
!Section 2: OPEN BINARY RESTART FILE
!------------------------------------------------------------------------------
 REWIND(nudat)
 OPEN (nudat, FILE=TRIM(datname), FORM='UNFORMATTED', STATUS='OLD', &
                   ACTION='READ', IOSTAT=izerr)
 IF (izerr /= 0) STOP

 READ (nudat, IOSTAT=izerr) psm0, dsem0, msem0, kem0, qcm0, ntke
 IF (izerr /= 0) STOP

 READ (nudat,IOSTAT=izerr) izvctype_read, refatm_p0sl,   refatm_t0sl,   &
                                 refatm_dt0lp, vcoord_vcflat, zvc_params
 IF (izerr /= 0) STOP 

 IF     ( (izvctype_read >   0) .AND. (izvctype_read <= 100) ) THEN
   PRINT*, "izvctype_read = ", izvctype_read 
 ELSEIF ( (izvctype_read > 100) .AND. (izvctype_read <= 200) ) THEN
   READ (nudat,IOSTAT=izerr) refatm_delta_t, refatm_h_scal
   IF (izerr /= 0) STOP 
 ELSEIF ( (izvctype_read > 200) .AND. (izvctype_read <= 300) ) THEN
   READ (nudat,IOSTAT=izerr) refatm_bvref
   IF (izerr /= 0) STOP 
 ENDIF 

!------------------------------------------------------------------------------
!Section 3: READ ALL RECORDS
!------------------------------------------------------------------------------

 read_loop: DO

   READ (nudat, IOSTAT=izerr) ipdsbuf, igdsbuf, rbuf

   WRITE (*,*) iz_countl, igdsbuf(4), ipdsbuf(8), MINVAL(rbuf), MAXVAL(rbuf)

   IF (izerr < 0) THEN
     WRITE (*,*) '       EOF REACHED'
     EXIT read_loop
   ELSEIF (izerr > 0) THEN
     WRITE (*,*) '    ERROR READING RESTART FILE'
     EXIT read_loop
   ENDIF
   iz_countl = iz_countl + 1

 ENDDO read_loop

 CLOSE (nudat, IOSTAT=izerr)
END
