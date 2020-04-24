! Copyright 2006-2014 ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities granted to it by
! virtue of its status as an intergovernmental organisation nor does it submit to any jurisdiction.
!
MODULE YOMCMEMPAR

! Module containing the model *SETUP* parameters,
!  derived parameters and universal constants

! Model Parameter options :

!   Physical Parameters :
!     fghz : microwave frequency (GHz)
!     theta : incidence angle (degrees)

!   Model Configurations: See the NOTICE file for references
!     CIDIEL    dielectric mixing model : 
!               1 : 'Wang' Wang and Schmugge model
!               2 : 'Dobson' Dobson model
!               3 : 'Mironov' Mironov model
!     CITEFF    effective temperature parametrization  
!               0 : 'Tsoil' teff = tsoil1
!               1 : 'Choudhury' 
!               2 : 'Wigneron' 
!               3 : 'Holmes' 
!               4 : 'Lv' et al., 2014
!     CISMR     model for Smooth Surface Emissivity
!               1 : 'Fresnel' (as described by Njoku) 
!               2 : 'Wilheit'
!     CIRGHR    surface roughness model
!               0 : 'No' Use Smooth Surface Emissivity
!               1 : 'Choudhury'  
!               2 : 'Wsimple' 
!               3 : 'Wegmueller' 
!               4 : 'Wtexture' 
!               5 : 'Wigneron'  
!     CIVEG     vegetation opacity model
!               0 : 'No' no vegetation
!               1 : 'Kirdyashev' 
!               2 : 'Wegmueller' 
!               3 : 'Wigneron'   
!               4 : 'Jackson'    
!     CIATM     atmospheric radiative transfer model
!               0 : 'No' no atmosphere and no cosmic background radiation
!               1 : 'Pellarin' 
!               2 : 'Ulaby'
!    CITDIEL : Temperature for dielectric model
!               1 : 'Teff'   tdiel=teff   
!               2 : 'Tsurf'  tdiel=tsurf  
!    CITVEG : Temperature of vegetation
!               1 : 'Tsurf'   tveg=tsurf   
!               2 : 'Tair'    tveg=tair    
!               3 : 'Tir'     tveg=tskin           
!    CIDVEG : vegetation cover input data
!               1 : 'Ecoclimap'  ecoclimap
!               2 : 'Tessel' ECMWF constant LAI
!               3 : 'HTessel' ECMWF varying LAI

!   Aditional soil model parameters :
!    nlay_soil_mw : number of soil layers in the MW emission model
!    nlay_soil_ls : number of soil layers in the Land surface model
!   ip_rgh_surf : soil surface roughness [cm]
!              (Attention: surface roughness in Choudhury (~0.35)
!              is not the same as in Wegmueller (~0-5cm))
!  ip_Q : Cross polarization parameter


! Derived Parameters :
!   f : microwave frequency (Hz)
!   costheta : cosine of incedence angle theta
!   sintheta : sine of incedence angle theta
!   lam : wavelength (m)
!   k : wavenumber (1/m)
!   lamcm : wavelength (cm)
!   kcm : wavenumber (1/cm)
!   omega : radian frequency (omega = 2. * pi * f), with f in Herz 

! Universal constants
!   c : speed of light (m/s)
!   tfreeze : freezing point of water (K)
!   rhowat : water density (kg/m3)
!   pi
!---------------------------------------------------------------------------
!    NAME       TYPE      PURPOSE
!    ----    :  ----   : ----------------------------------------------

! NAMOPT:
! CIDIEL   C : CHARACTER: IDENTIFIES THE DIELCTRIC MODEL
! CITEFF   C : CHARACTER: IDENTIFIES THE EFFECTIVE TEMPRATURE MODEL 
! CISMR    C : CHARACTER: IDENTIFIES THE SOIL MODEL 
! CIRGHR   C : CHARACTER: IDENTIFIES THE ROUGHNESS MODEL
! CIVEG    C : CHARACTER: IDENTIFIES THE VEGETATION OPACITY MODEL
! CIATM    C : CHARACTER: IDENTIFIES THE ATMOSPHERIC OPACITY MODEL
! CITVEG   C : CHARACTER: IDENTIFIES THE VEGETAION TEMPERATURE
! CIDVEG   C : CHARACTER: IDENTIFIES THE INPUT VEGETATION DATA SET
! CITDIEL  C : CHARACTER: IDENTIFIES THE TEMPERATURE FOR DIELCTRIC MODEL 
! CNAMEID  C : CHARACTER: NAME OF SIMUALTION, SUMMARIZES THE 8 FIRST OPTIONS CHOOSEN
! CFREQ    C : CHARACTER: Frequency of observation (GHz)
! CANGLE   C : CHARACTER: Angle of observation (degrees)

! NAMDEF:
! LOFFCMEM L : LOGICAL: False when offline
! LOMASK_OCEAN L : LOGICAL: False when global 
! LOMASK_AUTO  L : LOGICAL: True to remove wrong points (example TSKIN <=0K) 
! LOFIELDEXP L : LOGICAL: False when grid, true for 1pt -  only with ascii forcing 

! LGPRINT  L : LOGICAL: TRUE IF DEBUG INFO IS REQUIRED
! JPHISTLEV I : Number of output files
!             : Level1 (TB, Teff)
!             : Level2 (TB, Teff, VWC, frbare)
!             : Level3 (TB, Teff, VWC, frbare,salwat,Bpa, WP, falpha, fh, ffraci)
! CFINOUT   C : File format for input/output files: 'netcdf','grib,'ascii'

!NAMERAD:
! THETA    R : REAL: Incidence angle (in degres)
! FGHZ     R : REAL: Microwave frequency in Ghz

! NAMELEV:
! NLAY_SOIL_MW I: INTEGER: Number of soil layer in the MW model in input data
! NLAY_SOIL_LS I: INTEGER: Number of soil layer in the Land Surface Model




USE PARKIND1  ,ONLY : JPIM,JPRM 
IMPLICIT NONE


INTEGER(KIND=JPIM), PARAMETER::JPCHARLEN=200             !! Length of CNAMEID
CHARACTER(LEN=JPCHARLEN) ::  CIDIEL, CITEFF, CISMR, CIRGHR, CIVEG, CIATM, CITVEG, CIDVEG, CITDIEL
INTEGER(KIND=JPIM), PARAMETER::JPNAMEIDLEN=16             !! Length of CNAMEID
CHARACTER(LEN=JPNAMEIDLEN) :: CNAMEID
CHARACTER(LEN=3) CFREQ
CHARACTER(LEN=2) CANGLE

INTEGER(KIND=JPIM), PARAMETER::JPCMEMTILE=7             !! number of tiles in CMEM




INTEGER(KIND=JPIM)             :: JPHISTLEV      !! Defines levels of outputs
LOGICAL                        :: LGPRINT        !! Print debug info
LOGICAL                   :: LOFFCMEM       !! Offline mode
LOGICAL                   :: LOFIELDEXP       !! For field experiment application (1pixel)
LOGICAL                   :: LOMASK_OCEAN       !! For ocean masking in case no input data 
LOGICAL                   :: LOMASK_AUTO       !! For masking in case of wrong input data 
CHARACTER(LEN=JPCHARLEN) ::  CFINOUT        !! Input/output file format

!!!type(SATELLITE) :: SAT    ! PSG: to hold Sensor's all information
CHARACTER(LEN=JPCHARLEN) ::  INPUTNAMLST, INPUTSATINFO,CLMNAME,SURFNAME ! PSG: input files
!CHARACTER(LEN=JPCHARLEN) :: SURFNAME,CLMNAME,CLM_fname,surf_fname,inparam_fname

REAL(KIND = JPRM) :: fghz        
REAL(KIND = JPRM) :: theta      

INTEGER(KIND=JPIM) :: nlay_soil_mw
INTEGER(KIND=JPIM) :: nlay_soil_ls



REAL(KIND = JPRM) :: ip_sal_soil
REAL(KIND = JPRM) :: ip_rgh_surf
REAL(KIND = JPRM) :: ip_wcveg(2)
REAL(KIND = JPRM) :: ip_Q

REAL(KIND = JPRM) :: f
REAL(KIND = JPRM) :: costheta
REAL(KIND = JPRM) :: sintheta
REAL(KIND = JPRM) :: lam
REAL(KIND = JPRM) :: k
REAL(KIND = JPRM) :: lamcm
REAL(KIND = JPRM) :: kcm
REAL(KIND = JPRM) :: omega

REAL(KIND = JPRM) :: c
REAL(KIND = JPRM) :: tfreeze
REAL(KIND = JPRM) :: rhowat
REAL(KIND = JPRM) :: pi
REAL(KIND = JPRM) :: PPG 


 
END MODULE YOMCMEMPAR
