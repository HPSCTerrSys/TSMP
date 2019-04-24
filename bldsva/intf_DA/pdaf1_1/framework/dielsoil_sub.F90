! Copyright 2006-2014 ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities granted to it by
! virtue of its status as an intergovernmental organisation nor does it submit to any jurisdiction.
!
!---------------------------------------------------------------------------
! Models to calculate the dielectric constant of a soil medium following:
!   1. Wang and Schmugge
!   2. Dobson 
!   3. Mironov
!---------------------------------------------------------------------------

SUBROUTINE DIELWANG (ew)

! Purpose :
!   Calculate the dielectric constant of a wet soil 
!   Developed and validated for 1.4 and 5 GHz.

! External :
!   dielsal_sub.F90

! Input/Output Variables:
!   wt : transition moisture point (cm3/cm3)
!   gamma : fitting parameter
!   ecl : conductivity loss
!   ei : dielectric constant of ice, for initially absorbed water
!   ex : dielectric constant of the initially absorbed water
!---------------------------------------------------------------------------

USE PARKIND1, ONLY : JPIM, JPRM
USE YOMCMEMSOIL, ONLY : wc, eps, p, wp, alpha, wp, ea, er

IMPLICIT NONE

REAL(KIND = JPRM) :: wt
REAL(KIND = JPRM) :: gamma
REAL(KIND = JPRM) :: ecl
COMPLEX :: ew, ex, ei
!---------------------------------------------------------------------------

ei = (3.2, 0.1)

! 1. Calculate dielectric constant of soil-water mixture

gamma = -0.57_JPRM * wp + 0.481_JPRM
wt = 0.49_JPRM * wp + 0.165_JPRM
IF (wc <= wt) THEN
  ex = ei + (ew-ei) * (wc/wt) * gamma
  eps = wc*ex + (p-wc) * ea + (1.-p) * er
ELSE
  ex = ei + (ew-ei) * gamma
  eps = wt*ex + (wc-wt) * ew + (p-wc) * ea + (1.-p) * er
ENDIF

! 2. add conductivity loss

ecl = alpha * wc**2._JPRM
eps = eps + (0.,1.) * ecl

END SUBROUTINE DIELWANG

!===========================================================================

SUBROUTINE DIELDOBSON (ew)

! Purpose :
!   Calculate the dielectric constant of a wet soil 
!   Developed and validated for 1.4 and 18 GHz.
!   alphas : constant shape factor
!---------------------------------------------------------------------------

USE YOMLUN   , ONLY : NULOUT
USE PARKIND1, ONLY : JPIM, JPRM

USE YOMCMEMSOIL, ONLY : wc, eps, sand, clay, rho_s, rho_b
USE YOMCMEMPAR, ONLY :  LGPRINT

IMPLICIT NONE

REAL(KIND=JPRM) :: eps_s, beta, epsr, epsi, eaa
COMPLEX :: ew
REAL(KIND=JPRM), PARAMETER :: alphas = 0.65_JPRM
!---------------------------------------------------------------------------

! IF (wc < 0.001_JPRM) wc = 0.001_JPRM     ! to avoid dividing by zero
wc = MAX(wc,0.001_JPRM) ! to avoid dividing by zero

eps_s = (1.01_JPRM + 0.44_JPRM * rho_s)**2._JPRM - 0.062_JPRM  ! compare to ATBD: 4.7

beta = (127.48_JPRM - 0.519_JPRM * sand - 0.152_JPRM * clay) / 100._JPRM
eaa = 1.0_JPRM + (rho_b / rho_s) * (eps_s ** alphas - 1.0_JPRM)   &
  & + (wc ** beta) * (real(ew) ** alphas) - wc
epsr = eaa ** (1._JPRM/alphas)

beta = (133.797_JPRM - 0.603_JPRM * sand - 0.166_JPRM * clay) / 100._JPRM
eaa= (wc ** beta) * (abs(aimag(ew)) ** alphas)
epsi = eaa ** (1._JPRM/alphas)

eps = cmplx(epsr,epsi)

IF (LGPRINT) WRITE(NULOUT,*) '--- DOBSON:',eps,ew

END SUBROUTINE DIELDOBSON


!=====================================================================================

SUBROUTINE DIELMIRONOV 

! Purpose :
!   Calculate the dielectric constant of a wet soil
!   Developed and validated from 1 to 10 GHz.
!   adapted for a large range of soil moisture


USE PARKIND1,ONLY : JPIM,JPRM
USE YOMCMEMSOIL, ONLY : wc, eps, sand, clay, rho_s, rho_b, eps_0,eps_winf
USE YOMCMEMPAR, ONLY : f,pi

IMPLICIT NONE

REAL(KIND = JPRM) :: znd,zkd,zxmvt,zep0b,ztaub,zsigmab,zep0u,ztauu,zsigmau,zcxb,zepwux,zepwuy &
                    &, zcxu,zepwbx,zepwby,znb,zkb,znu,zku,znm,zkm,zflag,zxmvt2,zepmx,zepmy 

complex :: c  = ( 5.5,0.2) 

!INTEGER(KIND=JPIM) :: 


! Initializing the GRMDM spectroscopic parameters with clay (fraction)
!-----------------------------------------------------------------------

!  RI & NAC of dry soils
!-----------------------
znd = 1.634_JPRM - 0.539_JPRM * (clay/100._JPRM) + 0.2748_JPRM * (clay/100._JPRM)**2._JPRM
zkd = 0.03952_JPRM - 0.04038_JPRM * (clay / 100._JPRM)

! Maximum bound water fraction
!-----------------------------
zxmvt = 0.02863_JPRM + 0.30673_JPRM * clay / 100._JPRM

! Bound water parameters
!-----------------------
zep0b = 79.8_JPRM - 85.4_JPRM  * (clay / 100._JPRM) + 32.7_JPRM  * (clay / 100._JPRM)*(clay / 100._JPRM)
ztaub = 1.062e-11_JPRM + 3.450e-12_JPRM * (clay / 100._JPRM)
zsigmab = 0.3112_JPRM + 0.467_JPRM * (clay / 100._JPRM)

! Unbound (free) water parameters
!-------------------------------
zep0u = 100._JPRM
ztauu = 8.5e-12_JPRM
zsigmau = 0.3631_JPRM + 1.217_JPRM * (clay / 100._JPRM)


! Computation of epsilon water (bound & unbound)
!----------------------------------------------

zcxb = (zep0b - eps_winf) / (1._JPRM + (2._JPRM*pi*f*ztaub)**2._JPRM)
zepwbx = eps_winf + zcxb
zepwby =  zcxb * (2._JPRM*pi*f*ztaub) + zsigmab / (2._JPRM*pi*eps_0*f)

zcxu = (zep0u - eps_winf) / (1._JPRM + (2._JPRM*pi*f*ztauu)**2._JPRM)
zepwux = eps_winf + zcxu
zepwuy =  zcxu * (2._JPRM*pi*f*ztauu) + zsigmau/(2._JPRM*pi*eps_0*f)

! Computation of refractive index of water (bound & unbound)
! -----------------------------------------------------------

znb= sqrt( sqrt( zepwbx**2._JPRM + zepwby**2._JPRM) + zepwbx ) / sqrt(2._JPRM)
zkb= sqrt( sqrt( zepwbx**2._JPRM + zepwby**2._JPRM) - zepwbx ) / sqrt(2._JPRM)

znu= sqrt( sqrt( zepwux**2._JPRM + zepwuy**2._JPRM) + zepwux ) / sqrt(2._JPRM)
zku= sqrt( sqrt( zepwux**2._JPRM + zepwuy**2._JPRM) - zepwux ) / sqrt(2._JPRM)

! Computation of soil refractive index (nm & km):
! xmv can be a vector
! -----------------------------------------------------------

zxmvt2= min (wc, zxmvt)
zflag = 0._JPRM
IF ( wc >= zxmvt ) zflag = 1._JPRM

znm = znd + (znb - 1._JPRM) * zxmvt2 + (znu - 1._JPRM) * (wc-zxmvt) * zflag
zkm = zkd + zkb * zxmvt2 + zku * (wc-zxmvt) * zflag

! computation of soil dielectric constant:
! -----------------------------------------------------------

zepmx= znm ** 2._JPRM - zkm ** 2._JPRM     
zepmy= znm * zkm * 2._JPRM      

! As Prof. Su suggested: organic-rich
! The organic(matrix) part is according to Zheng et al 2017
! -----------------------------------------------------------
!zepmx= (wc*zepmx** 0.5+ (1-0.563)*5.5** 0.5+ (0.563-wc)*1** 0.5)**2._JPRM
!zepmy= (wc*zepmy** 0.5+ (1-0.563)*0.2** 0.5+ (0.563-wc)*1** 0.5)**2._JPRM

! The organic(matrix) part is according to Bircher et al 2016
! -----------------------------------------------------------
!zepmx= 50.69*wc ** 3._JPRM + 18.81*wc ** 2._JPRM + 25.0*wc + 1.636
!zepmy= 10.61*wc ** 3._JPRM - 11.08*wc ** 2._JPRM + 9.613*wc + 0.1211

eps = cmplx(zepmx,zepmy)
eps = (wc*sqrt(eps)+(0.563-wc)+(1-0.563)*sqrt(c))**2._JPRM


END SUBROUTINE DIELMIRONOV


