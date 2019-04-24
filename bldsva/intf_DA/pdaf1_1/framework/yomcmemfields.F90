! Copyright 2006-2014 ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities granted to it by
! virtue of its status as an intergovernmental organisation nor does it submit to any jurisdiction.
!
MODULE YOMCMEMFIELDS

! Module containing the I/O datafields

! ftair : Air temperature at screen level (2m) (K)
!   ftskin skin temperature
! ftveg : Canopy temperature (K)
! ftskin : Skin Temperature (K)
! ftl_lsm : Soil Temp on land surface model (lsm) layers (K)
! fwc_lsm : water content of soil on the lsm layers [m3/m3]

! ftb_soil  surface brightness temperature (K)
! ftb_toa  TOA brightness temperature (K)
! fsurf_emis : global surface emissivity [-]

! fsand [%]
! fclay [%]
! frho_b : bulk density  [kg/m3]
! fp : porosity [%]
! fWP : Wilting point [%]
! falpha : alpha for WANG dielectric model

! ftau_atm
! ftb_ad [K]
! ftb_au [K]

! fsnowd : snow depth, water equivalent (m)
! frsnow : snow density (kg/m3)
! fLSM(:)   
! fTVL(:)    
! fTVH(:)    
! ftfrac(:,:)  
! fwater : water fraction (0-1)
! mask(:) : point on which TB is going to be computed (1), or masked (0)
!  fs_laiL : global LAI for Low vegetation of ecoclimap
!  fs_laiH : global LAI for High vegetation of ecoclimap, NOT USED

! fwc_veg : vegetation water content for low (1) and high vegetation (2) [kg/m3]
! fh : effective roughness height [cm]
! fb : b parameter (Jackson)
! fsal_wat : salinity of water  (psu = ppt(weight) )
! ftauN: vegetation optical thickness at Nadir
! fteffC: diagnostif teff averaged for the grid box: C (1), TeffH(2), TeffV(3)

! Control variables :
!  N : number of data cells in field
!  JJ : current field cell of computation
!  JCMEMTILE : current tile
! for ascci input/outputs : doy,ndate,ntime,nbpt,NFCASCMAX

!---------------------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM,JPRM
IMPLICIT NONE

REAL(KIND = JPRM),ALLOCATABLE :: doy(:) 
INTEGER(KIND = JPIM),SAVE :: NFCASCMAX ! max size of ascii input 
INTEGER(KIND = JPIM),ALLOCATABLE :: NBPT(:) ! pixel id
INTEGER(KIND = JPIM),ALLOCATABLE :: ndate(:) 
INTEGER(KIND = JPIM),ALLOCATABLE :: ntime(:) 

REAL(KIND = JPRM),ALLOCATABLE :: ftair(:) 
REAL(KIND = JPRM),ALLOCATABLE :: ftveg(:) 
REAL(KIND = JPRM),ALLOCATABLE :: ftskin(:) 
REAL(KIND = JPRM),ALLOCATABLE :: ftl_lsm(:,:) 

REAL(KIND = JPRM),ALLOCATABLE :: fwc_lsm(:,:)    

REAL(KIND = JPRM),ALLOCATABLE :: ftb_soil(:,:)      
REAL(KIND = JPRM),ALLOCATABLE :: ftb_toa(:,:)       
REAL(KIND = JPRM),ALLOCATABLE :: fsurf_emis(:,:)    

REAL(KIND = JPRM),ALLOCATABLE :: fsand(:)    
REAL(KIND = JPRM),ALLOCATABLE :: fclay(:)    
REAL(KIND = JPRM),ALLOCATABLE :: frho_b(:)   
REAL(KIND = JPRM),ALLOCATABLE :: fp(:)   
REAL(KIND = JPRM),ALLOCATABLE :: fWP(:)  
REAL(KIND = JPRM),ALLOCATABLE :: falpha(:)  
REAL(KIND = JPRM),ALLOCATABLE :: fteffC(:,:)  

REAL(KIND = JPRM),ALLOCATABLE :: ftau_atm(:)  
REAL(KIND = JPRM),ALLOCATABLE :: ftb_au(:)  
REAL(KIND = JPRM),ALLOCATABLE :: ftb_ad(:)  

REAL(KIND = JPRM),ALLOCATABLE :: fsnowd(:) 
REAL(KIND = JPRM),ALLOCATABLE :: frsnow(:) 
REAL(KIND = JPRM),ALLOCATABLE :: fLSM(:)   
INTEGER(KIND=JPIM),ALLOCATABLE :: fTVL(:)    
INTEGER(KIND=JPIM),ALLOCATABLE :: fTVH(:)    
INTEGER(KIND=JPIM),ALLOCATABLE :: mask(:)    
REAL(KIND = JPRM),ALLOCATABLE :: ftfrac(:,:)  
REAL(KIND = JPRM),ALLOCATABLE :: fwater(:)    
REAL(KIND = JPRM),ALLOCATABLE :: fs_laiL(:)  
!REAL(KIND = JPRM),ALLOCATABLE :: fs_laiH(:)  

REAL(KIND = JPRM),ALLOCATABLE :: fNrh_L(:)
REAL(KIND = JPRM),ALLOCATABLE :: fNrh_H(:)
REAL(KIND = JPRM),ALLOCATABLE :: fNrv_L(:)
REAL(KIND = JPRM),ALLOCATABLE :: fNrv_H(:)

REAL(KIND = JPRM),ALLOCATABLE :: ftau_veg(:,:)  
REAL(KIND = JPRM),ALLOCATABLE :: ftth(:,:)
REAL(KIND = JPRM),ALLOCATABLE :: fttv(:,:)
REAL(KIND = JPRM),ALLOCATABLE :: fw_effL(:,:)
REAL(KIND = JPRM),ALLOCATABLE :: fw_effH(:,:)
REAL(KIND = JPRM),ALLOCATABLE :: fhrmodel(:,:)


REAL(KIND = JPRM),ALLOCATABLE :: fwc_veg(:,:) 
REAL(KIND = JPRM),ALLOCATABLE :: fh(:) 
REAL(KIND = JPRM),ALLOCATABLE :: fb(:,:) 
REAL(KIND = JPRM),ALLOCATABLE :: ftauN(:,:) 
REAL(KIND = JPRM),ALLOCATABLE :: fsal_wat(:) 
REAL(KIND = JPRM),ALLOCATABLE :: ftheta_inc(:,:)
! PSG: pixel-based incidence angle
real, allocatable, dimension(:,:,:,:,:) :: TB_HV 
! real(KIND = JPRM), allocatable, dimension(:,:,:,:,:) :: TB_HV
! PSG: TB with dims as [lons,lats,pol,inc,time]

INTEGER(KIND=JPIM) :: N
INTEGER(KIND=JPIM) :: JJ
INTEGER(KIND=JPIM) :: JCMEMTILE
INTEGER(KIND=JPIM) :: JJINC  ! PSG: index for incidence angle
!
CHARACTER(len=80) :: clname
!---------------------------------------------------------------------------

END MODULE YOMCMEMFIELDS
