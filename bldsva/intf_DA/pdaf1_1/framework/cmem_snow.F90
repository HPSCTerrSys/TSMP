! Copyright 2006-2014 ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities granted to it by
! virtue of its status as an intergovernmental organisation nor does it submit to any jurisdiction.
!
SUBROUTINE CMEM_SNOW(rsn,PDSW,DLUMI,esn)

! Purpose :
! -------
!   HUT-SNOW EMISSION MODEL FOR SINGLE SNOW LAYER
   


! Author :
! ------
!   2-Oct-2006 Thomas Holmes   *ECMWF*
!   January 2008 Patricia de Rosnay, ECMWF coding standards
!   November 2013 Patricia de Rosnay, include snow grain size function 

! Modifications :
! -------------
! End Modifications

! OUTPUT PARAMETERS:
! tb_tov = SNOW COVERED TERRAIN BRIGHTNESS TEMPERATURE AT H and V-POLARIZATION
! esn = EMISSIVITY OF SNOW COVERED TERRAIN AT 1: H-POLARIZATION 2: V-POLARIZATION

! INPUT:
! rsn   = REFLECTIVITY BETWEEN THE SNOW AND GROUND AT (1, H-POL. 2, V.)
! TLUMI = TEMPERATURE OF SNOW (DEGREES C)
! PDSW  = SNOW DENSITY (kg/M^3)
! DLUMI= THICKNESS OF SNOW LAYER (M)
! PPDIAM    = SNOW GRAIN SIZE (DIAMETER) (MM) from field
! PPDIAM computed using snow grain size function parameter ZGRAINA & ZGRAINB following:
! ZGRAINA = snow grain size function parameter a 
! ZGRAINB = snow grain size function parameter b 

! not used :  TSKY  = SKY BRIGHTNESS TEMPERATURE (K)
! PPQ2MOD : EMPIRICAL PARAMETER OF MODIFIED HUT MODEL
! eps_ice : dielectric constant of ice
! eir : REAL PART OF ICE EPSILON 
! eii : imaginary PART OF ICE EPSILON 
! sal_snow : SALINITY OF SNOW (PPM)
!------------------------------------------------------------------------------

USE PARKIND1, ONLY : JPIM, JPRM
USE YOMCMEMPAR, ONLY : fghz, omega, k, sintheta, costheta, tfreeze, rhowat, pi, LGPRINT
USE YOMCMEMATM, ONLY : tb_tov
USE YOMCMEMSOIL, ONLY : eps_winf, eps_0, XMV, tsoil
USE YOMLUN   , ONLY : NULOUT
IMPLICIT NONE

COMPLEX :: CIM, EPSA, EPSB, EPSC, EPSMARKA, CN2, RF

INTEGER(KIND=JPIM) :: i

REAL(KIND=JPRM) :: rsn(2)
REAL(KIND=JPRM) :: esn(2)
REAL(KIND=JPRM) :: PICE, Y0
!REAL(KIND=JPRM),INTENT(IN) :: PDSW,DLUMI
REAL(KIND=JPRM) :: TLUMI,PDSW,DLUMI
REAL(KIND=JPRM) :: PPDIAM 
REAL(KIND=JPRM), PARAMETER :: PPQ2MOD = 0.96
REAL(KIND=JPRM), PARAMETER :: PPsal_snow = 0.
REAL(KIND=JPRM) :: EE(2), XM, YSUORA(2)
REAL(KIND=JPRM) :: PDS, WEQ, REDS, XN1
REAL(KIND=JPRM) :: A, B, C, AP, BP, CP, AA, AB, AC
REAL(KIND=JPRM) :: EPSSW, FF0

REAL(KIND=JPRM) :: eir, eii, eiis, eiip, deltaeii
REAL(KIND=JPRM) :: VI, XIEDS, F0A, F0B, F0C, EPSINFA, EPSINFB, EPSINFC
REAL(KIND=JPRM) :: EPSSA, EPSSB, EPSSC, REWS, XIEWS, ALFA, BETA
REAL(KIND=JPRM) :: K2X, PPP, QQ, XJAKAJA, O2, TCOH, TSA, TRANS, SEC,KED,KE
REAL(KIND=JPRM) :: KABSD, KABS, XINTA, XINTE, SUHDE, L2, LA, XR1
REAL(KIND=JPRM) :: KS, L2APU, XML, EMISS1, EMISS2, SULAPU, SUL

REAL(KIND=JPRM), PARAMETER      :: ZGRAINA = 1.6e-4   
REAL(KIND=JPRM), PARAMETER      :: ZGRAINB = 1.1e-13  

COMPLEX :: eps_ice
!------------------------------------------------------------------------------
! LSN: SNOW DENSITY transfer from (kg/m^3) to (G/CM^3) 
PDSW = PDSW/1000.

IF (LGPRINT) WRITE(*,*) 'cmem_snow:',rsn,PDSW,DLUMI

! snow temperature [C] set equal to the surface soil temperature
TLUMI = tsoil(1)-tfreeze
! CONSTANTS:
 CIM=(0.,1.)
 PICE = 0.916 ! DENSITY OF ICE
 Y0 = 4.*pi*1.e-7! MYY NOLLA

PDS = (PDSW - XMV)/(1.0 - XMV) ! DENSITY OF DRY SNOW
WEQ = rhowat * DLUMI * PDSW   ! SNOW WATER EQUIVALENT
! TRANSMISSIVITY FROM GROUND TO SNOW (E = 1 - [S0V S0H])
EE(1) = 1. - rsn(1) ! H
EE(2) = 1. - rsn(2) ! V

! -------- IMAGINARY PART OF ICE EPSILON
CALL DIEL_ICE (TLUMI+tfreeze,eps_ice)
eir = REAL(eps_ice)
eii = AIMAG(eps_ice)
 ! -------- IMPURE ICE -5 C (MATZLER)
  A=0.0026
  B=0.00023
  C=0.87
  eiis=A/fghz+B*(fghz**C)
 ! -------- PURE ICE -5 C (MATZLER)
  AP=6.e-4               
  BP=6.5e-5
  CP=1.07
  eiip=AP/fghz+BP*(fghz**CP)
 ! -------- DIFFERENCE BETWEEN PURE AND IMPURE ICE AT -5 C
  deltaeii = eiis - eiip 
 ! -------- EFFECT OF SALINITY (USING OF THE DIFFERENCE AT -5 C):
 ! -------- IMAGINARY PART OF ICE EPSILON
  eii = eii + deltaeii*PPsal_snow/13. 

! -------- REAL PART OF DRY SNOW EPSILON
REDS = 1.+ 1.58*PDS/(1.-0.365*PDS)
XN1 = SQRT(Y0/eps_0)

! -------- IMAGINARY PART OF DRY SNOW EPSILON (PVS MODEL)
VI=PDS/PICE
XIEDS = 3.*VI*eii*(REDS**2.)*(2.*REDS+1.) /((eir+2.*REDS)*(eir+2.*REDS**2.))

! -------- DIELECTRIC CONSTANT OF WET SNOW (MATZLER)
IF (XMV > 0.) THEN
    AA = 0.005
    AB = 0.4975
    AC = 0.4975
    EPSSW = 88.
    FF0= 9.

    F0A = FF0 *(1.+ AA*(EPSSW-eps_winf)/(REDS+AA*(eps_winf-REDS)) )
    F0B = FF0 *(1.+ AB*(EPSSW-eps_winf)/(REDS+AB*(eps_winf-REDS)) )
    F0C = FF0 *(1.+ AC*(EPSSW-eps_winf)/(REDS+AC*(eps_winf-REDS)) )

    EPSINFA = (XMV/3.)*(eps_winf-REDS)/(1.+AA*((eps_winf/REDS)-1.))
    EPSINFB = (XMV/3.)*(eps_winf-REDS)/(1.+AB*((eps_winf/REDS)-1.))
    EPSINFC = (XMV/3.)*(eps_winf-REDS)/(1.+AC*((eps_winf/REDS)-1.))

    EPSSA = (XMV/3.) * (EPSSW-REDS)/(1.+AA*((EPSSW/REDS)-1.))
    EPSSB = (XMV/3.) * (EPSSW-REDS)/(1.+AB*((EPSSW/REDS)-1.))
    EPSSC = (XMV/3.) * (EPSSW-REDS)/(1.+AC*((EPSSW/REDS)-1.))

    EPSA = EPSINFA + (EPSSA-EPSINFA)/(1.+CIM*fghz/F0A)
    EPSB = EPSINFB + (EPSSB-EPSINFB)/(1.+CIM*fghz/F0B)
    EPSC = EPSINFC + (EPSSC-EPSINFC)/(1.+CIM*fghz/F0C)
       
    EPSMARKA = EPSA + EPSB + EPSC + (REDS-CIM*XIEDS)
    REWS = REAL(EPSMARKA)
    XIEWS = -1. * AIMAG(EPSMARKA)
ELSE
    REWS = REDS
    XIEWS = XIEDS
ENDIF
  
ALFA = k*ABS(AIMAG(SQRT(REWS-CIM*XIEWS)))
BETA = k*REAL(SQRT(REWS-CIM*XIEWS))

K2X = k * sintheta
PPP = 2. * ALFA * BETA
QQ = BETA**2. - ALFA**2. - (k**2.)*(sintheta)**2.
XJAKAJA = (1./(SQRT(2.)))*SQRT(SQRT(PPP**2.+QQ**2.)+QQ)
! -------- PROPAGATION ANGLE IN SNOW
O2 = ATAN(K2X/XJAKAJA)                             
!ETKULLUMI=O2*180./PI;      % NON USEFUL

! -------- WAVE IMBEDANCE IN SNOW
 CN2 = SQRT((Y0/eps_0)/(REWS-CIM*XIEWS))                  

! polarization H, V
DO i=1,2 
    
  ! -------- FRESNEL REFLECTION COEFFICIENTS BETWEEN SNOW AND AIR
  IF (i == 1) THEN
    RF = (CN2*costheta-XN1*COS(O2))/(CN2*costheta+XN1*COS(O2))
  ELSE
    RF = (XN1*costheta-CN2*COS(O2))/(XN1*costheta+CN2*COS(O2))
  ENDIF

  ! -------- POWER TRANSMISSION COEFFICIENT BETWEEN AIR AND SNOW:
  TCOH = (1. - (ABS(RF))**2.)

  TSA=TCOH
  TRANS=TSA

  ! -------- DRY SNOW EXTINCTION COEFFICIENT (DB/M => NP/M)
  SEC = 1./COS(O2)
  ! -------- IF fghz<60  KED=F(fghz,PPDIAM)
  PPDIAM = MIN(1000*(ZGRAINA + ZGRAINB*((PDSW*1000.0_JPRM)**4)),3.0_JPRM) 
  KED = 0.0018*(fghz**2.8)*(PPDIAM**2.0)
  KED = KED/4.3429  ! -------- (NP/M)

  ! -------- DRY SNOW ABSORPTION COEFFICIENT (NP/M)
  KABSD = 2.*omega*SQRT(Y0*eps_0*REDS) * SQRT(0.5*(SQRT(1.+(XIEDS/REDS)**2.)-1.))
  ! ----- CORRECTION TO EXTINCTION COEFFICIENT FOR SMALL FREQUENCIES (F<18 GHZ):
  IF (KED < KABSD)  KED = KABSD

  ! -------- ABSORPTION COEFFICIENT OF WET SNOW:
  KABS = 2.*omega* SQRT(Y0*eps_0*REWS) * SQRT(0.5*(SQRT(1.+(XIEWS/REWS)**2.)-1.))

  ! -------- ATTENUATION DUE TO ABSORPTION IN WET SNOW (NP/M):
  XINTA=0.001*KABS/PDSW*WEQ*SEC 

  ! -------- TOTAL EXTINCTION:             ! (ASSUMING THAT SCATTERING 
  !        COEFF. (KAPPAS) IS THE SAME AS IN THE CASE OF DRY SNOW!!!
   KE = (KED - KABSD) + KABS
  ! -------- TOTAL ATTENUATION (NP/M)
  XINTE = 0.001*KE/PDSW*WEQ*SEC           
  SUHDE=XINTE/XINTA

  ! -------- TOTAL ATTENUATION
  L2=EXP(XINTE)       
  ! -------- ATTENUATION DUE TO ABSORPTION
  LA=EXP(XINTA)       
  
  ! -------- REFLECTIVITY OF SNOW-AIR BOUNDARY:
  XR1 = 1. - TRANS

  ! -------- GROUND EMISSION CONTRIBUTION
  KS = KE - KABS
  L2APU = EXP((KE-PPQ2MOD*KS)*SEC*DLUMI)
  XML = EE(i) *tsoil(1)/L2APU * (1.-XR1)/(1.-(1.-rsn(i) )*(XR1)/L2APU**2.)
  YSUORA(1)=XML
  EMISS1 = YSUORA(1)/tsoil(1)

  ! -------- SNOW EMISSION CONTRIBUTION
  SULAPU = (1.-XR1) * (TLUMI+tfreeze) * (KABS/(KE-PPQ2MOD*KS)) * (1. - 1./L2APU)
  SUL = (1. + rsn(i)/L2APU) * SULAPU/(1.-rsn(i)*(XR1)/L2APU**2.)
  YSUORA(2)=SUL 
  EMISS2 = YSUORA(2)/(TLUMI+tfreeze)

  ! -------- BRIGHTNESS TEMPERATURE OF SNOW COVERED TERRAIN:
  tb_tov(i) = YSUORA(1)+YSUORA(2)

  ! -------- EMISSIVITY OF SNOW COVERED TERRAIN
  esn(i) = EMISS1 + EMISS2
  
ENDDO  
 
END SUBROUTINE CMEM_SNOW
