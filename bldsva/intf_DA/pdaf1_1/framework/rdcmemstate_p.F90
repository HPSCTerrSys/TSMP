! Copyright 2006-2014 ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities granted to it by
! virtue of its status as an intergovernmental organisation nor does it submit to any jurisdiction.
!
! This subroutine is created by Shaoning Lv
! it imitate CMEM ASCii reading subroutine to get data from state_p
! state_p is a vector in length of 2*nlay_soil_ls+18
 
SUBROUTINE RDCMEMSTATE_P(state_p) 

! Purpose :
! -------
!     Read input data for state_p vector for PDAF
   
! Externals :
! ---------

! Internal datafields :
! -------------------
  
!     Patricia de Rosnay  01-10-2007
!     Patricia de Rosnay  02-11-2009 SSS
!     Patricia de Rosnay  February 2011 fix in multi-layer atmospheric model
!     Patricia de Rosnay, 19-01-2012 add HTessel confi with variable LAI
!------------------------------------------------------------------------------

USE PARKIND1, ONLY : JPIM, JPRM
USE YOMLUN, ONLY : NULOUT, NULTMP
USE YOMCMEMPAR, ONLY : CIDVEG, CITVEG, CIATM,  LGPRINT,CNAMEID, LOFIELDEXP, PPG, &
                    & LOMASK_OCEAN, LOMASK_AUTO,nlay_soil_ls
USE YOMCMEMFIELDS, ONLY : N, JJ, doy,ndate,ntime,nbpt &
                    &  ,fTVL, fTVH,  fs_laiL, ftfrac, fsnowd, frsnow,ftl_lsm, ftair,ftveg &
                    &  ,ftskin,fwc_lsm,fsand,fclay,fwater,mask 
USE YOMCMEMATM, ONLY : fZ
USE YOMCMEMSOIL, ONLY : sal_sea 

IMPLICIT NONE

INTEGER(KIND=JPIM) :: i
REAL(KIND=JPRM) :: rvcov(0:20)
REAL(KIND=JPRM) :: rvlai(0:20)
REAL(KIND=JPRM) :: b(0:7)
REAL(KIND=JPRM) :: b1(0:7)
REAL(KIND=JPRM) :: b2(0:7)
REAL(KIND=JPRM) :: b3(0:7)
REAL(KIND=JPRM) :: VWC(0:7)
REAL(KIND=JPRM) :: ZNrh(0:7)
REAL(KIND=JPRM) :: ZNrv(0:7)
REAL(KIND=JPRM) :: Ztth(0:7)
REAL(KIND=JPRM) :: Zttv(0:7)
REAL(KIND=JPRM) :: Zhr(0:7)
REAL(KIND=JPRM) :: Zw_eff(0:7,2)
INTEGER(KIND=JPIM) :: RVTV(0:20)

REAL(KIND=JPRM),ALLOCATABLE :: cvegl(:), cvegh(:)
REAL(KIND=JPRM),ALLOCATABLE :: fbare(:)  
REAL(KIND=JPRM),ALLOCATABLE :: fvegl(:)  
REAL(KIND=JPRM),ALLOCATABLE :: fvegh(:)    
REAL(KIND=JPRM),ALLOCATABLE :: sncov(:)

REAL(KIND=JPRM) ::  zdoy,zsand,zclay, zfZ,zfTVL,zfTVH,zfvegl,zfvegh,zfwater,zlaiH

!LSN: ARGUMENTS:
REAL, INTENT(in)    :: state_p(dim_p)

!
!------------------------------------------------------------------------------
!
!
!Define default mask
!
mask(:) = 1_JPIM
!

!----------------------------------------------------------
! 1.  Read main input file 
!----------------------------------------------------------

! LSN:  OPEN(NULTMP,FILE='forcing_cmem_main.asc',status='old')
! LSN:  DO JJ = 1,N
! LSN:  READ (NULTMP,*) nbpt(JJ), ndate(JJ),ntime(JJ),(fwc_lsm(JJ,i),i=1,nlay_soil_ls),ftskin(JJ)&
! LSN:               & ,(ftl_lsm(JJ,i),i=1,nlay_soil_ls),fsnowd(JJ),frsnow(JJ),doy(JJ)
! LSN:  ENDDO 
! LSN:  CLOSE(NULTMP)
!-----------------------------! LSN 
! reading input from state_p vector, no loop
JJ = 1
nbpt(JJ)   = 1
ndate(JJ)  = 1
ntime(JJ)  = 1 
fwc_lsm(JJ,1:nlay_soil_ls) = state_p(2:nlay_soil_ls+1)
ftl_lsm(JJ,1:nlay_soil_ls) = state_p(nlay_soil_ls+2:2*nlay_soil_ls+1)
fsnowd(JJ) = state_p(2*nlay_soil_ls+2)
frsnow(JJ) = state_p(2*nlay_soil_ls+3)
ftskin(JJ) = state_p(2*nlay_soil_ls+4)
doy(JJ)    = 1
!-----------------------------! LSN end

fsnowd(:) = MAX(1.e-2_JPRM,fsnowd(:))
WHERE(fsnowd == 1.e-2_JPRM) fsnowd(:) = 0.0

! Update the mask to eliminate unrealistic points
   SELECT CASE (LOMASK_AUTO)
   CASE (.True.)
     do i = 1,nlay_soil_ls
      WHERE (ftl_lsm(:,i) <100.0_JPRM .or. ftskin(:) < 100.0_JPRM) mask(:) = 0_JPIM  ! unrealistic temperature
      WHERE ( fwc_lsm(:,i) < 0.0_JPRM ) mask(:) = 0_JPIM  ! Unrealistic SM
     enddo
   END SELECT

     ! in any case ensure sm is positive
     fwc_lsm(:,:) = MAX(0.0_JPRM, fwc_lsm(:,:))


!----------------------------------------------------------
! 2.  Read soil and vegetation input files
!----------------------------------------------------------

! 2.1 Read data texture and vegetation types
!---------------------------------------------

  ALLOCATE (fbare(N))
  ALLOCATE (fvegl(N))
  ALLOCATE (fvegh(N))
  ALLOCATE (cvegl(N))
  ALLOCATE (cvegh(N))
  ALLOCATE (sncov(N))

  ! default vegetation cover is 100 % on vegetated tiles
  cvegh(:) = 1.0_JPRM
  cvegl(:) = 1.0_JPRM



SELECT CASE (LOFIELDEXP)
       
   CASE (.TRUE.)  ! Field experiment

       ! constant parameters to be read (%, %, m2/s2)

       ! LSN:  OPEN(NULTMP,FILE='forcing_cmem_soil_constant.asc',status='old')
       ! LSN:  READ (NULTMP,*) zsand,zclay,zfZ
       ! LSN:  CLOSE(NULTMP)
       !-------------------------! LSN
       fsand(JJ) = state_p(2*nlay_soil_ls+5)
       fclay(JJ) = state_p(2*nlay_soil_ls+6)
       zfZ(JJ)    = state_p(2*nlay_soil_ls+7)
       !-------------------------! LSN end

       ! Convert geopotential at surface from (m2/s2) to height (km):    
       zfZ = zfZ / 1000._JPRM / PPG

       ! Update the mask to eliminate unrealistic points
       SELECT CASE (LOMASK_AUTO)
         CASE (.True.)
         WHERE (fZ > 9.0_JPRM .or. fZ < -1.0_JPRM ) mask = 0_JPIM  ! unexpected elevation on the Earth
       END SELECT

       SELECT CASE (CIDVEG)

          CASE ('Tessel','HTessel')

           ! fraction in (-)
           ! LSN: OPEN(NULTMP,FILE='forcing_cmem_veg-tessel_constant.asc',status='old')
           ! LSN: READ (NULTMP,*) zfTVL,zfTVH,zfvegl,zfvegh,zfwater
           ! LSN: CLOSE(NULTMP)
           !-------------------------! LSN
           zfTVL     = state_p(2*nlay_soil_ls+11)
           zfTVH     = state_p(2*nlay_soil_ls+10)
           zfvegl    = state_p(2*nlay_soil_ls+9)
           zfvegh    = state_p(2*nlay_soil_ls+8)
           zfwater   = state_p(2*nlay_soil_ls+12)
           !-------------------------! LSN end

          CASE ('Ecoclimap')

           ! fraction in (-)
           ! LSN: OPEN(NULTMP,FILE='forcing_cmem_veg-ecoclimap_constant.asc',status='old')
           ! LSN: WRITE(NULOUT,*) 'Read forcing_cmem_veg'
           ! LSN: READ (NULTMP,*) zfTVL,zfTVH,zfvegl,zfvegh,zfwater
           ! LSN: CLOSE(NULTMP)
           !-------------------------! LSN
           zfTVL     = state_p(2*nlay_soil_ls+11)
           zfTVH     = state_p(2*nlay_soil_ls+10)
           zfvegl    = state_p(2*nlay_soil_ls+9)
           zfvegh    = state_p(2*nlay_soil_ls+8)
           zfwater   = state_p(2*nlay_soil_ls+12)
           !-------------------------! LSN end


       END SELECT

       ! Apply constant values to 1D file
       fsand = zsand
       fclay = zclay
       fZ = zfZ
       fTVL = zfTVL
       fTVH = zfTVH
       fvegl = zfvegl
       fvegh = zfvegh
       fwater = zfwater

  CASE (.FALSE.)

       ! Surface conditions vary
        
       ! LSN: OPEN(NULTMP,FILE='forcing_cmem_soil.asc',status='old')
       ! LSN:  DO JJ = 1,N
       ! LSN:   READ (NULTMP,*) nbpt(JJ),ndate(JJ),ntime(JJ),fsand(JJ),fclay(JJ),fZ(JJ),doy(JJ)
       ! LSN:  ENDDO
       !-------------------------! LSN
       fsand(JJ) = state_p(2*nlay_soil_ls+5)
       fclay(JJ) = state_p(2*nlay_soil_ls+6)
       fz(JJ)    = state_p(2*nlay_soil_ls+7)
       !-------------------------! LSN end
       
       ! Convert geopotential at surface from (m2/s2) to height (km):    
       fZ = fZ / 1000._JPRM / PPG
       ! CLOSE(NULTMP) !LSN

       SELECT CASE (CIDVEG)

         CASE ('Tessel','HTessel')

           ! LSN: OPEN(NULTMP,FILE='forcing_cmem_veg-tessel.asc',status='old') 
           ! LSN: WRITE(NULOUT,*) 'Read forcing_cmem_veg-tessel' 
           ! LSN: DO JJ = 1,N ! LSN
           ! LSN:   READ (NULTMP,*) nbpt(JJ),ndate(JJ),ntime(JJ),fTVL(JJ),fTVH(JJ),fvegl(JJ),fvegh(JJ),fwater(JJ),doy(JJ) 
           ! LSN: ENDDO ! LSN
           ! LSN: CLOSE(NULTMP) 

           !-------------------------! LSN
           fTVL     = state_p(2*nlay_soil_ls+11)
           fTVH     = state_p(2*nlay_soil_ls+10)
           fvegl    = state_p(2*nlay_soil_ls+9)
           fvegh    = state_p(2*nlay_soil_ls+8)
           fwater   = state_p(2*nlay_soil_ls+12)
           !-------------------------! LSN end

         CASE ('Ecoclimap')

           ! LSN: OPEN(NULTMP,FILE='forcing_cmem_veg-ecoclimap.asc',status='old')
           ! LSN: WRITE(NULOUT,*) 'Read forcing_cmem_veg-ecoclimap'
           ! LSN: 
           ! LSN: DO JJ = 1,N
           ! LSN:   READ (NULTMP,*) nbpt(JJ),ndate(JJ),ntime(JJ),fTVL(JJ),fTVH(JJ),fvegl(JJ),fvegh(JJ),fwater(JJ),doy(JJ)
           ! LSN: ENDDO
           ! LSN: CLOSE(NULTMP)
           
           !-------------------------! LSN
           fTVL     = state_p(2*nlay_soil_ls+11)
           fTVH     = state_p(2*nlay_soil_ls+10)
           fvegl    = state_p(2*nlay_soil_ls+9)
           fvegh    = state_p(2*nlay_soil_ls+8)
           fwater   = state_p(2*nlay_soil_ls+12)
           !-------------------------! LSN end
         
           WHERE ( fTVL > 7_JPIM .and. fwater(:) <  0.5_JPRM ) fTVL = 4_JPIM  ! Defaults low vegetation is grass
           WHERE ( fTVL > 7_JPIM .and. fwater(:) >= 0.5_JPRM ) fTVL = 0_JPIM  ! No vegetation on  water pixel
           WHERE ( fTVH > 7_JPIM .and. fwater(:) <  0.5_JPRM ) fTVH = 1_JPIM  ! Defaults low vegetation is grass
           WHERE ( fTVH > 7_JPIM .and. fwater(:) >= 0.5_JPRM ) fTVH = 0_JPIM  ! No vegetation on  water pixel

       END SELECT


  

END SELECT  ! LOFIELDEXP


! 2.2  Read lai 
!--------------


SELECT CASE (CIDVEG)

   CASE ('Ecoclimap')

    !OPEN(NULTMP,FILE='forcing_cmem_lai.asc',status='old')

    ! LSN:  DO JJ = 1,N
    ! LSN: READ (NULTMP,*) nbpt(JJ),ndate(JJ),ntime(JJ),fs_laiL(JJ),zlaiH,doy(JJ)
    ! LSN: ENDDO
    ! LSN: CLOSE(NULTMP)
    !-----------------------------! LSN
    fs_laiL(JJ) = state_p(2*nlay_soil_ls+13)
    zlaiH       = state_p(2*nlay_soil_ls+13)
    !-----------------------------! LSN end

    CASE ('Tessel')

    CALL VEGTABLE (RVLAI,RVCOV,RVTV,b,VWC,b1,b2,b3,ZNrh,ZNrv,Ztth,Zttv,Zhr,Zw_eff)
    fs_laiL(:) = rvlai(fTVL(:))
    ! TESSEL vegetation map needs to account for bare soil
    ! Once LAI is computed, convert to ECOCLIMAP
    cvegl(:) = rvcov(fTVL(:)) 
    cvegh(:) = rvcov(fTVH(:))  
    ! Convert model vegetation types to CMEM vegetation requirement
    fTVL(:) = RVTV(fTVL(:))
    fTVH(:) = RVTV(fTVH(:))

    CASE ('HTessel')
    
    ! read LAI
    ! LSN:  OPEN(NULTMP,FILE='forcing_cmem_lai.asc',status='old') !LSN
    ! LSN:  DO JJ = 1,N !LSN
    ! LSN:  READ (NULTMP,*) nbpt(JJ),ndate(JJ),ntime(JJ),fs_laiL(JJ),zlaiH,doy(JJ) !LSN
    ! LSN:  ENDDO !LSN
    ! LSN:  CLOSE(NULTMP) !LSN
    !-----------------------------! LSN
    fs_laiL(JJ) = state_p(2*nlay_soil_ls+13)
    zlaiH       = state_p(2*nlay_soil_ls+13)
    !-----------------------------! LSN end

    ! HTESSEL vegetation map needs to account for bare soil
    ! Once LAI is computed, convert to ECOCLIMAP
    CALL VEGTABLE (RVLAI,RVCOV,RVTV,b,VWC,b1,b2,b3,ZNrh,ZNrv,Ztth,Zttv,Zhr,Zw_eff)
    cvegl(:) = rvcov(fTVL(:)) 
    cvegh(:) = rvcov(fTVH(:))  
    ! Convert model vegetation types to CMEM vegetation requirement
    fTVL(:) = RVTV(fTVL(:))
    fTVH(:) = RVTV(fTVH(:))



END SELECT


! 2.3  Calculate tile fractions
!------------------------------

    fvegl(:) =  fvegl(:) * cvegl(:)      
    fvegh(:) =  fvegh(:) * cvegh(:)
    fbare(:) = 1. - fvegl(:) - fvegh(:)  

    ! Snow cover on land tiles

    WHERE (  fsnowd(:) > 0.0_JPRM) sncov(:) = 1.0_JPRM  
    WHERE (  fsnowd(:) <= 0.0_JPRM) sncov(:) = 0.0_JPRM  


    ftfrac(:,1) = fbare(:) * (1. - fwater(:)) * (1.-sncov(:)) ! bare soil 
    ftfrac(:,2) = fbare(:) * (1. - fwater(:)) * sncov(:)      ! bare soil with snow   
    ftfrac(:,3) = fvegl(:) * (1. - fwater(:)) * (1.-sncov(:)) ! low veg
    ftfrac(:,4) = fvegl(:) * (1. - fwater(:)) * sncov(:)      ! low veg with snow
    ftfrac(:,5) = fvegh(:) * (1. - fwater(:)) * (1.-sncov(:)) ! high veg
    ftfrac(:,6) = fvegh(:) * (1. - fwater(:)) * sncov(:)      ! high veg with snow
    ftfrac(:,7) = fwater(:)                                ! water


   SELECT CASE (LOMASK_OCEAN)
       CASE (.True.)
        WHERE ( fwater(:) >= 0.5 )  mask(:) = 0.0_JPRM
    END SELECT



  DEALLOCATE (fbare)
  DEALLOCATE (fvegl)
  DEALLOCATE (fvegh)
  DEALLOCATE (cvegl)
  DEALLOCATE (cvegh)
  DEALLOCATE (sncov)
  

!-------------------------------------------------------------------
! 3.  Read atmospheric temperature input  (Only CASE CITVEG = 'Tair' or CASE CIATM = 'Pellarin')
!-------------------------------------------------------------------


SELECT CASE (CITVEG)

   CASE ('Tair')

     ! LSN:  OPEN(NULTMP,FILE='forcing_cmem_tair.asc',status='old') ! LSN
     ! LSN: DO JJ = 1, N ! LSN
     ! LSN:  READ (NULTMP,*) nbpt(JJ),ndate(JJ),ntime(JJ),ftair(JJ),doy(JJ) ! LSN
     ! LSN: 
     ! LSN: ENDDO ! LSN
     ! LSN: CLOSE(NULTMP) ! LSN     
     !-----------------------------! LSN
     ftair(JJ) = state_p(2*nlay_soil_ls+14)
     !-----------------------------! LSN end

     ftveg(:) = ftair(:)

!    Update the mask to eliminate unrealistic points
     SELECT CASE (LOMASK_AUTO)
       CASE (.True.)
       WHERE (ftveg < 100.0_JPRM ) mask = 0_JPIM  ! Unrealistic Temperature
     END SELECT
!

   CASE ('Tsurf')

      !-----------------------------! LSN
     ftveg(:) = state_p(1*nlay_soil_ls+1)
     !-----------------------------! LSN end
!     ftveg(:) = ftl_lsm(:,1)

END SELECT

SELECT CASE (CIATM)

   CASE ('Pellarin', 'Ulaby')

     ! LSN:  OPEN(NULTMP,FILE='forcing_cmem_tair.asc',status='old') 
     ! LSN: DO JJ = 1, N
     ! LSN:  READ (NULTMP,*) nbpt(JJ),ndate(JJ),ntime(JJ),ftair(JJ),doy(JJ)
     ! LSN: 
     ! LSN:  ENDDO
     ! LSN: CLOSE(NULTMP)
     !-----------------------------! LSN
     ftair(JJ) = state_p(2*nlay_soil_ls+14)
     !-----------------------------! LSN end

!    Update the mask to eliminate unrealistic points
      SELECT CASE (LOMASK_AUTO)
       CASE (.True.)
       WHERE (ftair < 100.0_JPRM ) mask = 0_JPIM  ! Unrealistic Temperature
      END SELECT
!


END SELECT

!-------------------------------------------------------------------
! 4.  Read SSS
!-------------------------------------------------------------------

sal_sea(:) = 32.5_JPRM

         
END SUBROUTINE RDCMEMSTATE_P
