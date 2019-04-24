! Copyright 2006-2014 ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities granted to it by
! virtue of its status as an intergovernmental organisation nor does it submit to any jurisdiction.
!
! vegetation opacity models (iveg)
!               0 = no vegetation
!               1 = Kirdyashev 
!               2 = Wegmueller 
!               3 = Wigneron
!---------------------------------------------------------------------------
!  About vegetation structure coefficient a_geo:
!
!    a_geo(pol) = u(pol) / 3
!
!    u takes vegetation geometry into account
!
!    Theoretical examples:
!       horizontal leaves:  a_geo(h) = 1     a_geo(v) = (cos(theta))**2
!       isotropic leaves:   a_geo(h) = 2/3   a_geo(v) = 2/3
!       vertical stalks:    a_geo(h) = 0     a_geo(v) = (sin(theta))**2
!       isotropic stalks:   a_geo(h) = 1/3   a_geo(v) = 1/3
!---------------------------------------------------------------------------

SUBROUTINE VEGKIRD 

! Purpose :
!  Calculate vegetation opacity using 'Effective Medium theory'
!  For low frequencies (<7.5GHz)
!---------------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM
USE YOMLUN   , ONLY : NULOUT
USE YOMCMEMPAR, ONLY : costheta, k, rhowat, LGPRINT
USE YOMCMEMVEG, ONLY : wc_veg, tau_veg, eps_vw, a_geo

IMPLICIT NONE

INTEGER(KIND=JPIM) :: i
!---------------------------------------------------------------------------

DO i = 1, 2
  tau_veg(i) = a_geo(i) * k * (wc_veg / rhowat) * aimag(eps_vw) / costheta
ENDDO
!WRITE(6,*) eps_vw, (wc_veg / rhowat)

END SUBROUTINE VEGKIRD

!===========================================================================

SUBROUTINE VEGWEGM

! Purpose :
!  Calculate vegetation opacity using 'Geometrical Optics theory'
!  For low to high  frequencies (<100GHz)

! Internal Variables:
!  rho_veg : vegetation density (kg/m3)
!  m_d     : dry mass fraction of vegetation
!  d_leaf  : leaf thickness (m)
!  eps_veg : Dielectric constant of leaves
!  BB      : wet biomass
!---------------------------------------------------------------------------

USE YOMLUN   , ONLY : NULOUT
USE PARKIND1  ,ONLY : JPIM,JPRM
USE YOMCMEMPAR, ONLY : theta, pi, k, costheta, sintheta, LGPRINT
USE YOMCMEMVEG, ONLY : wc_veg, tau_veg, eps_vw, a_geo

IMPLICIT NONE

REAL(KIND = JPRM), PARAMETER :: rho_veg = 950._JPRM
REAL(KIND = JPRM), PARAMETER :: m_d = 0.3_JPRM
REAL(KIND = JPRM), PARAMETER :: d_leaf = 0.2e-3_JPRM
REAL(KIND = JPRM) :: kz0
REAL(KIND = JPRM) :: th, tv
REAL(KIND = JPRM) :: BB

COMPLEX(KIND = JPRM) :: eps_veg
COMPLEX(KIND = JPRM) :: kz1
COMPLEX(KIND = JPRM) :: Rh, Rv, j, hfrac, vfrac, exp1, exp2
!---------------------------------------------------------------------------

eps_veg = (0.51_JPRM - 0.64_JPRM * m_d) * eps_vw  +  3.2_JPRM * m_d  +  0.49_JPRM

! 1. Calculate transmissivity of a single leaf
!    -----------------------------------------
kz0 = k * costheta
kz1 = k * sqrt(eps_veg - sintheta * sintheta)

Rh = (kz0 - kz1) / (kz0 + kz1)
Rv = (eps_veg * kz0 - kz1) / (eps_veg * kz0 + kz1)

j = (0._JPRM,-1._JPRM)
exp1 = exp(j * (kz0 - kz1) * d_leaf)
exp2 = exp(-j * 2._JPRM * kz1 * d_leaf)

hfrac = (4._JPRM * kz0 * kz1 * exp1) / ((kz0 + kz1)**2._JPRM * (1. - Rh**2._JPRM * exp2))
vfrac = (4._JPRM * eps_veg * kz0 * kz1 * exp1) /   &
        ((eps_veg * kz0 + kz1)**2._JPRM * (1._JPRM - Rv**2._JPRM * exp2))

th = abs(hfrac) ** 2._JPRM
tv = abs(vfrac) ** 2._JPRM

! 2. Calculate vegetation opacity
!    ----------------------------
BB = wc_veg / (1._JPRM - m_d)   
tau_veg(1) = a_geo(1) * k * ( BB / rho_veg ) * aimag(eps_veg) * th / costheta
tau_veg(2) = a_geo(2) * k * ( BB / rho_veg )  * aimag(eps_veg) * tv / costheta
!WRITE(6,*) wc_veg, m_d, eps_vw, eps_veg , ( B / rho_veg )

END SUBROUTINE VEGWEGM

!===========================================================================

SUBROUTINE VEGWIGN

! Purpose :
!  CALCULATE VEGETATION OPACITY for frequencies < 10GHz
! tauN related to wc_veg 
! tauV and tauH both related to tauN (Wigneron)
!    generalized case; ttH=ttV=1 --> tauH=tauV=tauN
!      tau assumed to be independent of both polarization and incedence angle

! Internal Variables:
!  tauN : nadir value of vegetation optical depth (index 1.h-pol, 2.v-pol)
!  ttH, ttV : polarization correction parameters
!---------------------------------------------------------------------------

USE YOMLUN   , ONLY : NULOUT
USE PARKIND1  ,ONLY : JPIM,JPRM
USE YOMCMEMPAR, ONLY : costheta , sintheta, LGPRINT
USE YOMCMEMVEG, ONLY : wc_veg, tauN, tau_veg, bj,tth,ttv

IMPLICIT NONE


! 1. Polarized vegetation optical depth
!    ----------------------------------
tau_veg(1) = tauN * (costheta**2 + ttH * sintheta**2)
tau_veg(2) = tauN * (costheta**2 + ttV * sintheta**2)


IF (LGPRINT) WRITE(NULOUT,*) '--- Tauh, tauv:',tau_veg
IF (LGPRINT) WRITE(NULOUT,*) '--- Tau nadir:',tauN

! 2. Calculate vegetation opacity
!    ----------------------------
tau_veg(1) = tau_veg(1) / costheta
tau_veg(2) = tau_veg(2) / costheta


END SUBROUTINE VEGWIGN
!===========================================================================

SUBROUTINE VEGJACK

! Purpose :
!  CALCULATE VEGETATION OPACITY for frequencies < 10GHz
! tauN related to wc_veg (Jackson)
! tauV and tauH both related to tauN (Wigneron, )


! Internal Variables:
!  tauN : nadir value of vegetation optical depth (index 1.h-pol, 2.v-pol)
!  ttH, ttV : polarization correction parameters
!---------------------------------------------------------------------------

USE YOMLUN   , ONLY : NULOUT
USE PARKIND1  ,ONLY : JPIM,JPRM
USE YOMCMEMPAR, ONLY : costheta , sintheta, LGPRINT
USE YOMCMEMVEG, ONLY : wc_veg, tauN, tau_veg, bj,tth,ttv

IMPLICIT NONE

!REAL(KIND = JPRM) :: tauN
!REAL(KIND = JPRM),parameter :: ttH=1.
!REAL(KIND = JPRM),parameter :: ttV=1.
!---------------------------------------------------------------------------

tauN = bj * wc_veg   


! 1. Polarized vegetation optical depth
!    ----------------------------------
tau_veg(1) = tauN 
tau_veg(2) = tauN 

IF (LGPRINT) WRITE(NULOUT,*) '--- Tauh, tauv:',tau_veg
IF (LGPRINT) WRITE(NULOUT,*) '--- Tau nadir:',tauN

! 2. Calculate vegetation opacity
!    ----------------------------
tau_veg(1) = tau_veg(1) / costheta
tau_veg(2) = tau_veg(2) / costheta


END SUBROUTINE VEGJACK
