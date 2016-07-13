!+ Source Code for reading cosmo binary restart files
!------------------------------------------------------------------------------

program readCOSMObinr

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
!
!> @TODO ... what is the table and code for pressure or pressure perturbations
!> @TODO ... what is the table and code for cloud ice
!> @TODO ... what temperature is the prognostic variable
!> @TODO ... vertical ... everything ...
!
!==============================================================================

IMPLICIT NONE

INTEGER, PARAMETER :: intgribf  = KIND(1)
INTEGER, PARAMETER :: iintegers = KIND(1)
INTEGER, PARAMETER :: sp = SELECTED_REAL_KIND( 6, 37) 
INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(12,307) 
INTEGER, PARAMETER :: wp = dp

! Table to decode the record contents of igdsbuf
integer, parameter :: indx_numEW    =  5, &
                      indx_numNS    =  6, &
                      indx_startlat =  7, &
                      indx_startlon =  8, &
                      indx_endlat   = 10, &
                      indx_endlon   = 11

! Table to decode the record contents of ipdsbuf
!> @TODO no seconds?
integer, parameter :: indx_gribver   =   2, &
                      indx_var       =   7, &
                      indx_zlevtyp   =   8, &
                      indx_zlevtop   =   9, &
                      indx_zlevbot   =  10, &
                      indx_year      =  11, &
                      indx_mm        =  12, &
                      indx_dd        =  13, &
                      indx_hh        =  14, &
                      indx_min       =  15, &
                      indx_startstep =  17, &
                      indx_endstep   =  18, &
                      indx_nztri     =  19, &
                      indx_cc        =  22

INTEGER (KIND=intgribf),  PARAMETER ::  &
              npds  = 321_intgribf, & ! dimension of product definition section 
              ngds  = 626_intgribf    ! dimension of grid description section 

INTEGER (KIND=iintegers), PARAMETER ::  &
              idim  = 30, &           ! dimesion of longitude
              jdim  = 20, &           ! dimension of latitutde
              kdim  = 249              ! dimension of vertical co-ordinate

INTEGER (KIND=iintegers) ::  &
              nudat,      &           ! unit number
              izerr,      &           ! error status
              ivar,       &           ! variable reference number based on iver
              iver,       &           ! version number of GRIB1 indicator table
              iz_countl               ! counter for binary data

integer :: i, nx, ny, ilev, ilevp1, ilevtyp
integer :: icc, iyy, imm, idd, ihh, imin, iccyy
integer :: istartstep, iendstep, nztri
real(kind=wp) lat1, latN, lon1, lonN

REAL(KIND=wp) ::  &
              rbuf(idim*jdim),  &     ! data to be read
              psm0,             &     ! initial value for mean surface pressure ps
              dsem0,            &     ! initial value for mean dry static energy
              msem0,            &     ! initial value for mean moist static energy
              kem0,             &     ! initial value for mean kinetic energy
              qcm0,             &     ! initial value for mean cloudwater content
              refatm_p0sl,      &     ! constant reference pressure on sea-level
              refatm_t0sl,      &     ! constant reference temperature on sea-level
              refatm_dt0lp,     &     ! d (t0) / d (ln p0)
              refatm_delta_t,   &     ! temperature difference between sea level and stratosphere (for irefatm=2)
              refatm_h_scal,    &     ! scale height (for irefatm=2)
              refatm_bvref,     &     ! constant Brund-Vaisala-frequency for irefatm=3
              vcoord_vcflat,    &     ! coordinate where levels become flat
              zvc_params(kdim+1)      ! vertical coordinate - heights?

INTEGER(KIND=intgribf) :: &
             ipdsbuf(npds), &         ! pds: product definition section
             igdsbuf(ngds)            ! gds: grid definition section

INTEGER(KIND=iintegers) :: &
             ntke,         &          ! time level for TKE
             izvctype_read            ! check vertical coordinate type in restarts

CHARACTER(LEN=*), PARAMETER :: datname = "lrff00030000o"

!------------------------------------------------------------------------------
!Variable Defininiton Complete, program starts here

!------------------------------------------------------------------------------
!Section 1: INITIALIZE
!------------------------------------------------------------------------------

write(*,*)'iintegers is ',iintegers
 
nudat = 100

!------------------------------------------------------------------------------
!Section 2: OPEN BINARY RESTART FILE
!------------------------------------------------------------------------------

OPEN (nudat, FILE=TRIM(datname), FORM='UNFORMATTED', STATUS='OLD', &
                  ACTION='READ', IOSTAT=izerr)
IF (izerr /= 0) then
   write(*,*)'ERROR: unable to open file "'//trim(datname)//'".'
   STOP
endif

REWIND(nudat)

READ (nudat, IOSTAT=izerr) psm0, dsem0, msem0, kem0, qcm0, ntke
IF (izerr /= 0) then
   write(*,*)'ERROR: unable to read first record.'
   STOP
endif

write(*,*)'first record had ', psm0, dsem0, msem0, kem0, qcm0, ntke

READ (nudat,IOSTAT=izerr) izvctype_read, refatm_p0sl,   refatm_t0sl,   &
                                refatm_dt0lp, vcoord_vcflat, zvc_params
IF (izerr /= 0) then
   write(*,*)'ERROR: unable to read second record.'
   STOP
endif

write(*,*)'second record had', izvctype_read, refatm_p0sl,   refatm_t0sl,   &
                                refatm_dt0lp, vcoord_vcflat, zvc_params

IF     ( (izvctype_read >   0) .AND. (izvctype_read <= 100) ) THEN
  PRINT*, "izvctype_read = ", izvctype_read 
ELSEIF ( (izvctype_read > 100) .AND. (izvctype_read <= 200) ) THEN
  READ (nudat,IOSTAT=izerr) refatm_delta_t, refatm_h_scal
  write(*,*)'aux record was ',refatm_delta_t, refatm_h_scal
  IF (izerr /= 0) STOP 
ELSEIF ( (izvctype_read > 200) .AND. (izvctype_read <= 300) ) THEN
  READ (nudat,IOSTAT=izerr) refatm_bvref
  write(*,*)'aux record was ', refatm_bvref
  IF (izerr /= 0) STOP 
ELSEIF ( izvctype_read > 300 ) THEN
  write(*,*)'ERROR: izvctype_read is ',izvctype_read,' which is greater than max expected (300)'
  stop
ENDIF 

!------------------------------------------------------------------------------
!Section 3: READ ALL RECORDS
!------------------------------------------------------------------------------

iz_countl = 0 

read_loop: DO

  READ (nudat, IOSTAT=izerr) ipdsbuf, igdsbuf, rbuf

  IF (izerr < 0) THEN
    WRITE (*,*) ' EOF reached after ',iz_countl,' data records.'
    exit read_loop
  ELSEIF (izerr > 0) THEN
    WRITE (*,*) 'ERROR READING RESTART FILE around data record ',iz_countl
    exit read_loop
  ENDIF

  iz_countl = iz_countl + 1

  ! Decode ipdsbuf

  iver    = ipdsbuf(indx_gribver)
  ivar    = ipdsbuf(indx_var)
  ilevtyp = ipdsbuf(indx_zlevtyp)
  ilev    = ipdsbuf(indx_zlevtop)
  ilevp1  = ipdsbuf(indx_zlevbot)

  icc     = ipdsbuf(indx_cc)-1
  iyy     = ipdsbuf(indx_year)
  iccyy   = iyy + icc*100
  imm     = ipdsbuf(indx_mm)
  idd     = ipdsbuf(indx_dd)
  ihh     = ipdsbuf(indx_hh)
  imin    = ipdsbuf(indx_min)

  istartstep = ipdsbuf(indx_startstep)
  iendstep   = ipdsbuf(indx_endstep)
  nztri      = ipdsbuf(indx_nztri)

! write(*,'(A,3x,i8,5(1x,i2),3(1x,i4))') 'time for iz_countl = ', &
!        iz_countl,icc,iyy,imm,idd,ihh,imin,nztri,istartstep,iendstep 

  ! Decode igdsbuf

  nx   = igdsbuf(indx_numEW)
  ny   = igdsbuf(indx_numNS)
  lat1 = real(igdsbuf(indx_startlat),wp) * 0.001_wp
  lon1 = real(igdsbuf(indx_startlon),wp) * 0.001_wp
  latN = real(igdsbuf(indx_endlat  ),wp) * 0.001_wp
  lonN = real(igdsbuf(indx_endlon  ),wp) * 0.001_wp

  ! Since Fortran binary reads have no way to know if they fail or not
  ! we are going to compare our input with the sizes encoded in the file.

  if (nx .ne. idim .or. ny .ne. jdim) then
     write(*,*)'ERROR: (file) nx /= idim ',nx,idim, ' or '
     write(*,*)'ERROR: (file) ny /= jdim ',ny,jdim
     stop
  endif

  ! debug code ... 3300 does not exist
  if (iver == 2 .and. ivar == 3300) then
     do i = 1,npds
        write(*,*)iver, ivar, 'ipdsbuf(',i,') = ',ipdsbuf(i)
     enddo

     do i = 1,ngds
        write(*,*)iver, ivar, 'igdsbuf(',i,') = ',igdsbuf(i)
     enddo
  endif

  write (*,'( 7(1x,i4),2(1xf15.8))') iver, ivar, nx, ny, ilev, ilevp1, ilevtyp, MINVAL(rbuf), MAXVAL(rbuf)

ENDDO read_loop

CLOSE (nudat, IOSTAT=izerr)

end program readCOSMObinr

