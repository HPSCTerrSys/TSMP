! Copyright 2006-2014 ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities granted to it by
! virtue of its status as an intergovernmental organisation nor does it submit to any jurisdiction.
!
SUBROUTINE CMEM_SETUP 

! Purpose
! -------
!   Initialize CMEM model for computing microwave emission of the surface 

! Interface
! ---------
!   none

! Method
! ------
! First check if used online or offline
! Define defalut values
! Read Namelists
! Check option consitency
! Define runing options and output filename CNAMEID
!
! Externals
! ---------

! Reference
! ---------

! Author
! ------
!  November 2006 T. Holmes
!  September 2007 P. de Rosnay 
!  January 2008 P. de Rosnay
!  January 2009 P. de Rosnay (Grib API as default option)
!  November 2009 P. de Rosnay (V3.0 Salinity moved out from setup to RD files)
!  February 2012 P. de Rosnay  (Tessel option, bug fix in TEFF and VEGTABLE)

!---------------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRM
USE YOMLUN   , ONLY : NULOUT, NULTMP,NULNAM

USE YOMCMEMPAR
USE YOMCMEMATM, ONLY : t_sky
USE YOMCMEMVEG, ONLY : sal_vw, a_geoL, a_geoH
USE YOMCMEMSOIL
USE YOMCMEMFIELDS, ONLY :  NFCASCMAX

IMPLICIT NONE
INTEGER(KIND=JPIM) :: nlay_soil_ls_default 
INTEGER(KIND=JPIM) :: i

#include "namopt.h"
#include "namrad.h"
#include "namlev.h"
#include "namdef.h"

!---------------------------------------------------------------------------
! 1.0 Offline mode: if CMEM is not used in the ECMWF IFS, then initialize unit numbers:
!---------------------------------------------------------------------------

!Default run definition:
LOFFCMEM = .True.       ! Offline
LOFIELDEXP = .False.     ! Grid simulation
LOMASK_AUTO = .True.    ! Mask points with non physical inputs
LOMASK_OCEAN = .False.  ! Mask Ocean
LGPRINT = .False.       ! No debug info
JPHISTLEV = 1_JPIM      ! Basic output (TBH, TBH, TEFF)
CFINOUT = 'clm'        ! Default Input/output is GribAPI
INPUTSATINFO = '/p/project/chbn29/hbn29q/command/SATOPE_CLM4CMEM-master/forcing/SMOS_L1orbit_neckar.nc' ! PSG: default SAT

IF (LOFFCMEM) NULTMP = 77
IF (LOFFCMEM) NULNAM = 14 

! Open namelist:
OPEN(NULNAM,FILE='/p/project/chbn29/hbn29q/terrsysmp/bldsva/intf_DA/pdaf1_1/framework/input',STATUS='OLD')

! Read Definition namelist
READ (NULNAM,NML=namdef)

WRITE(NULOUT,*) 'Run definition :'
WRITE(NULOUT,*) 'OFFLINE = :',LOFFCMEM
WRITE(NULOUT,'(a22,a7)') ' Input/output format = ',CFINOUT
WRITE(NULOUT,*) 'Field experiment (1 grid point) = :',LOFIELDEXP
WRITE(NULOUT,*) 'Automatic Masking of wrong points = ',LOMASK_AUTO
WRITE(NULOUT,*) 'Ocean Masking = ',LOMASK_OCEAN
WRITE(NULOUT,*) 'SAELLITE INFO FILE = ',INPUTSATINFO ! PSG

!----------------------------
! 2.0 Default *SETUP* parameters
!----------------------------

! 2.1 Default model parameterization options (NAMOPT):

CIDIEL = 'Wang'
CITEFF = 'Wigneron'
CISMR = 'Fresnel'
CIRGHR = 'Wsimple' 
CIVEG = 'Wigneron'
CIATM = 'Pellarin' 
CITVEG = 'Tsurf'
CIDVEG = 'HTessel'
CITDIEL = 'Teff'

! Summarizes default options: 

CNAMEID = CIDIEL(1:2)//CITEFF(1:2)//CISMR(1:2)//CIRGHR(1:2)//CIVEG(1:2)//CIATM(1:2)//CITVEG(1:2)//CIDVEG(1:2)
WRITE(NULOUT,*)  'Defaults options: ', CNAMEID


! 2.2 Default radiometrique configuration (NAMRAD):

fghz = 1.4_JPRM 
theta = 50._JPRM 

! 2.3 Define the soil input soil layer depth (NAMLEV):

nlay_soil_mw = 1.0_JPRM  ! number of soil layer for radiometric model
nlay_soil_ls_default = 3.0_JPRM  ! number of soil layer for radiometric model
nlay_soil_ls   = nlay_soil_ls_default  ! default number of soil layer for radiometric model
ALLOCATE (z_lsm(nlay_soil_ls_default)) ! Land Surface model layers
! HTESSEL
z_lsm(1) = 0.07_JPRM        ! input soil layers depths
z_lsm(2) = 0.28_JPRM
z_lsm(3) = 1.0_JPRM
WRITE(NULOUT,*) 'LSM layers depth: ',z_lsm(:)

! ORCHIDEE
!z1 = 0.02151_JPRM
!z2 = 0.1642_JPRM
!z3 = 1.314_JPRM
NFCASCMAX = 2000000 ! max size of ascii input can be changed


! ----------------
! 3.0 read namelist
! ----------------


REWIND(NULNAM)
READ (NULNAM,NML=namopt)
REWIND(NULNAM)
READ (NULNAM,NML=namrad)
REWIND(NULNAM)
READ (NULNAM,NML=namlev)


WRITE(CFREQ,'(I3.3)') NINT(fghz*10._JPRM)
WRITE(CANGLE,'(I2)') INT(THETA) 
WRITE(NULOUT,*) 'Number of Layer in the MW model: ',nlay_soil_mw
WRITE(NULOUT,*) 'Number of Layer in the LSM: ',nlay_soil_ls

IF ( nlay_soil_ls .ne. nlay_soil_ls_default ) Then

 SELECT CASE (CFINOUT)
 !
  CASE ('netcdf') ! netcdfcase
 !   
  DEALLOCATE(z_lsm) ! netcdfcase
  ALLOCATE (z_lsm(nlay_soil_ls)) ! netcdfcase
  OPEN(NULTMP,FILE='LSM_VERTICAL_RESOL.asc',status='old') ! netcdfcase
  READ (NULTMP,*) (z_lsm(i),i=1,nlay_soil_ls)  ! netcdfcase
  WRITE(NULOUT,*) ' NETCDF IO: LSM layers depth: ',z_lsm(:) ! netcdfcase
  CLOSE(NULTMP)    ! netcdfcase
 !
 !-------------------------------------! PSG
  CASE('clm') ! PSG: for clm case
   DEALLOCATE(z_lsm) ! PSG: CLMcase
   ALLOCATE (z_lsm(nlay_soil_ls)) ! PSG: CLMcase    
   ! PSG: the elements for z_lsm will be passed by memory from
   ! PSG: the subroutine: memory_cmem_forcing()
 !-------------------------------------! PSG end 
 !
 !CASE ('ifs')
 !
 !
 END SELECT


ELSE

  WRITE(NULOUT,*) 'DEFAULT vertical resolution:'
  WRITE(NULOUT,*) 'HTESSEL 3-layer input considered'

ENDIF

! 3.1 check namelist option exist
! -------------------------------

IF (CIDIEL/= 'Wang'.AND.CIDIEL/= 'Dobson'.AND.CIDIEL/= 'Mironov') CALL ABOR1('Wrong CIDIEL choice. Choose Wang, Dobson or Mironov.')
IF (CITEFF/= 'Tsoil'.AND.CITEFF/= 'Choudhury'.AND.CITEFF/= 'Wigneron'.AND.CITEFF/= 'Holmes' &
     & .AND.CITEFF/= 'Lv') CALL ABOR1('Wrong CITEFF choice.')
IF (CISMR/= 'Fresnel'.AND.CISMR/= 'Wilheit') CALL ABOR1('Wrong CISMR choice. Choose Fresnel or Wilheit.')
IF(CIRGHR/= 'No'.AND.CIRGHR/= 'Choudhury'.AND.CIRGHR/= 'Wsimple'.AND.CIRGHR/= 'Wegmueller' & 
     & .AND. CIRGHR/= 'Wtexture' .AND. CIRGHR/= 'Wigneron') &
     &   CALL ABOR1('Wrong CIRGHR choice.')
IF(CIVEG/= 'No'.AND.CIVEG/= 'Kirdyashev'.AND.CIVEG/= 'Wegmueller'.AND.CIVEG/= 'Wigneron'&
     & .AND.CIVEG/= 'Jackson') CALL ABOR1('Wrong CIVEG choice.')
IF(CIATM/= 'No'.AND.CIATM/= 'Pellarin'.AND.CIATM/= 'Ulaby') CALL ABOR1('Wrong CIATM choice.')
IF(CITDIEL/= 'Teff'.AND.CITDIEL/= 'Tsurf') CALL ABOR1('Wrong CITDIEL choice. Choose Teff or Tsurf')
IF (CIDVEG /= 'Tessel'.AND.CIDVEG /= 'HTessel'.AND.CIDVEG /= 'Ecoclimap') &
     & CALL ABOR1('Wrong CIDVEG choice. Choose Tessel or HTessel or Ecoclimap.')
IF(CITVEG/= 'Tsurf'.AND.CITVEG/= 'Tair'.AND.CITVEG/= 'Tir') CALL ABOR1('Wrong CITVEG choice. Choose Tsurf, Tair or Tr')

!-------------------------------------! PSG: following line includes 'clm' as option
IF(CFINOUT/= 'netcdf'.and.CFINOUT/= 'clm'.and.CFINOUT/= 'pdaf' )&
& CALL ABOR1('Wrong CFINOUT choice.  Choose netcdf or ascii')
!-------------------------------------
! 4 Set Constants, and make conversions
!-------------------------------------

 PPG = 9.81_JPRM                                ! Gravity constant m/s2
 C = 2.998E8_JPRM                               ! Light speed
 tfreeze = 273.15_JPRM                          !
 rhowat = 1000._JPRM  
 pi = acos(-1._JPRM)
 f = fghz * 1E9_JPRM
 omega = 2._JPRM * pi * f
 lam = c / f           
 k = 2._JPRM * pi / lam    
 lamcm = lam * 100._JPRM 
 kcm = k / 100._JPRM 
 costheta = cos(theta*pi/180._JPRM)
 sintheta = sin(theta*pi/180._JPRM)

! Set atmospheric constants 
t_sky = 2.7_JPRM

! Set vegetation constants 
sal_vw = 6.0_JPRM

a_geoL(:) = (/0.33,0.33/)
a_geoH(:) = (/0.66,0.66/)

! Set soil constants
! CMEM conv   rho_s = 2.664_JPRM
rho_s = 2.66_JPRM

sal_soil = 0._JPRM
!sal_sea = 32.5_JPRM
ip_rgh_surf = 2.2_JPRM 

w0 = 0.41_JPRM  
bw = 0.35_JPRM  


eps0 = 0.08_JPRM
bh = 0.87_JPRM
eps_winf  = 4.9_JPRM
eps_0 = 8.854e-12_JPRM
ea = (1.0, 0.0)
er = (5.5, 0.2)  
ef = (5.0, 0.5)

! Derived soil parameters

  ip_Q = 0.35_JPRM * (1.0_JPRM - exp(-0.6_JPRM * ip_rgh_surf**2._JPRM * fghz))
  IF ( fghz < 2._JPRM ) ip_Q=0._JPRM   ! Q is assumed zero at low frequency



!-----------------------------------
! 5.0 Check consistency between options
!-----------------------------------

SELECT CASE (CITEFF)
  CASE ( 'Choudhury' ) 
    IF ( fghz<2.5 ) THEN
      WRITE(NULOUT,*) '*WARNING* Choudhury model neglects the strong soil moisture'
      WRITE(NULOUT,*) '          effect on sensing depth at Frequencies below 2.5 GHZ'
    ENDIF
  CASE ( 'Wigneron','Holmes' ) 
    IF ( fghz>2.5 ) THEN
      WRITE(NULOUT,*) '*WARNING* Effective temperature model calibrated for L-band Frequencies'
      WRITE(NULOUT,*) '          CITEFF set to Choudhury'
      CITEFF = 'Choudhury'
    ENDIF
ENDSELECT

IF ((CIDIEL == 'Wang' .AND. fghz > 11.) .OR. (CIDIEL == 'Dobson' .AND. fghz > 20.)) THEN
  WRITE(NULOUT,*) '*WARNING* Dielectric model not validated'
ENDIF

SELECT CASE (CISMR)
  CASE ( 'Fresnel' ) 
    IF ( nlay_soil_mw /= 1_JPIM) THEN
      nlay_soil_mw = 1_JPIM
      WRITE(NULOUT,*) '*WARNING* Fresnel model uses only one soil layer :'
      WRITE(NULOUT,*) '          nlay_soil_mw set to 1'
    ENDIF
  CASE ( 'Wilheit' ) 
    IF ( nlay_soil_mw < 2_JPIM ) THEN
      nlay_soil_mw = 10_JPIM 
      WRITE(NULOUT,*) '*WARNING* Wilheit model requires more than one soil layer :'
      WRITE(NULOUT,*) '          nlay_soil_mw set to default value of 10'
    ENDIF
    WRITE(NULOUT,*) 'Wilheit model, number of layers is:' , nlay_soil_mw 
    CITEFF = 'Tsoil'  ! wilhiet model calculats own Tef
    CITDIEL = 'Tsurf' ! wilhiet model calculats dielectric constant per soillayer
ENDSELECT

IF ((CIRGHR == 'Wsimple' .OR. CIRGHR == 'Wtexture' .OR. CIRGHR == 'Wigneron' ) .AND. fghz >= 2. ) THEN
  WRITE(NULOUT,*) '*WARNING* Roughness model not valid'
  IF (fghz>10.) THEN 
    WRITE(NULOUT,*) '          changed to Wegmueller'
    CIRGHR = 'Wegmueller'
  ELSE
    WRITE(NULOUT,*) '          changed to Choudhury'
    CIRGHR = 'Choudhury'
  ENDIF    
ELSEIF (CIRGHR == 'Wegmueller' .AND. theta >= 70. ) THEN
  CALL ABOR1('CMEM_SETUP: CIRGHR == Wegmueller AND theta >= 70 ')
ENDIF

IF ((CIVEG == 'Kirdyashev' .AND. fghz > 7.5) .OR. (CIVEG == 'Wigneron' .AND. fghz > 11.)) THEN
  WRITE(NULOUT,*) '*WARNING* Vegetation model not valid'
  WRITE(NULOUT,*) '          changed to Wegmueller'
  CIVEG = 'Wegmueller'
ELSEIF (CIVEG == 'Wigneron' .AND. fghz > 2.5 ) THEN
  WRITE(NULOUT,*) '*WARNING* Vegetation b-parameter is for L-band'
ENDIF

IF (CIATM == 'Pellarin' .AND. fghz >= 10. ) THEN
  WRITE(NULOUT,*) '*WARNING* Atmospheric model is only valid for freqeuncies below 10GHz'
  WRITE(NULOUT,*) '          CIATM set to Ulaby'
  CIATM = 'Ulaby'
ENDIF

IF ( LOFIELDEXP .eqv. .True. .AND. CFINOUT /= 'ascii') CALL ABOR1('Stop: ascii format required for field experiments')


!---------------------------
! 6 Retained run configuration 
!---------------------------


CNAMEID = CIDIEL(1:2)//CITEFF(1:2)//CISMR(1:2)//CIRGHR(1:2)//CIVEG(1:2)//CIATM(1:2)//CITVEG(1:2)//CIDVEG(1:2)

WRITE(NULOUT,*)  'Runing options: ', CNAMEID



END SUBROUTINE CMEM_SETUP
