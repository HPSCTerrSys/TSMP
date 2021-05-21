! Copyright 2006-2014 ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities granted to it by
! virtue of its status as an intergovernmental organisation nor does it submit to any jurisdiction.
!
SUBROUTINE DIEL_WAT (medium,isal,T,sal,ew)

! Purpose : 
!   Calculate dielectric constant of water in three different media : 
!   pure water, sea water, soil water


! Interface :
!  medium = pure water(0) sea water(1) soil water(2)
!  isal = Stogryn (1) Klein and Swift (2)

! local variables :
!  N : normality from salinity (Stogryn, modified by Klein and Swift)
!  T : temperature of water (C)
!   ew  : dielectric constant of water
!  sal : water salinity (psu = ppt(weight) )
!  eps_w0 : static dielectric constant of pure water (Klein and Swift) 
!  eps_sw0 : Static dielectric constant of soil water
!  tau_w : relaxation time of pure water (Stogryn)
!  tau_sw : relaxation time of saline water
!  sigma : ionic conductivity
!  sigma_eff : effective conductivity of water (S/m)
!---------------------------------------------------------------------------

USE PARKIND1, ONLY : JPIM, JPRM
USE YOMCMEMPAR, ONLY : omega, LGPRINT
USE YOMCMEMSOIL, ONLY : eps_winf, eps_0, rho_b , rho_s, sand, clay, wc

IMPLICIT NONE

INTEGER(KIND=JPIM) :: medium, isal
REAL(KIND=JPRM) :: T
REAL(KIND=JPRM) :: sal
REAL(KIND=JPRM) :: N
REAL(KIND=JPRM) :: sigma, sigma_eff
REAL(KIND=JPRM) :: tau_w, tau_sw
REAL(KIND=JPRM) :: eps_w0, eps_sw0, a, bb
!COMPLEX(KIND=JPRM) :: ew, j
COMPLEX :: ew, j
!---------------------------------------------------------------------------

j = (0.,1.)

tau_w = 1.768e-11_JPRM  - 6.068e-13_JPRM * T  + 1.104e-14_JPRM * T**2_JPRM  - 8.111e-17_JPRM * T**3_JPRM
! same as:
!tau_w = 1./(2.*pi) * (1.1109e-10_JPRM - 3.824e-12_JPRM * T + 6.938e-14_JPRM * T**2_JPRM  &
!     - 5.096e-16_JPRM * T**3_JPRM)

SELECT CASE (isal)
  
  CASE ( 1 )
    N = 0.9141_JPRM * sal * (1.707e-2_JPRM + 1.205e-5_JPRM * sal + 4.058e-9_JPRM * sal**2_JPRM)
    
    eps_sw0 = 87.74_JPRM  - 0.4008_JPRM * T  + 9.398e-4_JPRM * T**2_JPRM  + 1.410e-6_JPRM * T**3_JPRM
    a = 1.0_JPRM  - 0.2551_JPRM * N  + 5.151e-2_JPRM * N**2_JPRM  - 6.889e-3_JPRM * N**3_JPRM
    eps_sw0 = eps_sw0 * a

    bb = 1.0_JPRM - 0.04896_JPRM * N - 0.02967_JPRM * N**2_JPRM + 5.644e-3_JPRM &
       & * N**3_JPRM + 0.1463e-2_JPRM * N * T
    tau_sw  = tau_w * bb
   
  CASE ( 2 )
  
    eps_sw0 = 87.134_JPRM  - 1.949e-1_JPRM * T  - 1.276e-2_JPRM * T**2_JPRM  + 2.491e-4_JPRM * T**3_JPRM
    a = 1.000_JPRM  + 1.613e-5_JPRM * sal * T  - 3.656e-3_JPRM * sal  + 3.210e-5_JPRM * sal**2_JPRM  &
             -  4.232e-7_JPRM * sal**3_JPRM
    eps_sw0 = eps_sw0 * a
    
    bb = 1.000_JPRM  + 2.282e-5_JPRM * sal * T  - 7.638e-4_JPRM * sal  - 7.760e-6_JPRM * sal**2_JPRM  &
               + 1.105e-8_JPRM * sal**3_JPRM
    tau_sw  = tau_w * bb
   
END SELECT

SELECT CASE (medium)

  CASE ( 0 ) ! pure water
    eps_w0 = 88.045_JPRM - 0.4147_JPRM * T + 6.295e-4_JPRM * T**2_JPRM + 1.075e-5_JPRM * T**3_JPRM 
    ew = eps_winf + (eps_w0 - eps_winf) / (1._JPRM - j * omega * tau_w)
  CASE ( 1 )
    CALL ION_CONDUCT (T,sal,sigma)
    ! Debye expression '41
    ew = eps_winf + (eps_sw0 - eps_winf) / (1._JPRM - j * omega * tau_sw)  &
         + j * sigma / (omega * eps_0)
  
  CASE ( 2 )
    sigma_eff = -1.645_JPRM  + 1.939_JPRM * rho_b  - 0.02256_JPRM * sand  + 0.01594_JPRM * clay
    !  sigma_eff gets negative for very sandy soils
    !   with low bulk densities. If this happens set sigma_eff to zero.
    IF (sigma_eff < 0.) sigma_eff = 0._JPRM
    
    ! Modified Debye expression, Dobson '85
    !IF (wc < 0.001_JPRM) wc = 0.001_JPRM     ! to avoid dividing by zero
     wc = MAX(0.001_JPRM, wc)  ! to avoid dividing by zero
    ew = eps_winf + (eps_sw0 - eps_winf) / (1._JPRM - j * omega * tau_sw)  &
         + j * sigma_eff / (omega * eps_0) * (rho_s - rho_b) / (rho_s * wc)

END SELECT

END SUBROUTINE DIEL_WAT

!===========================================================================
! Internal subroutines
!===========================================================================

SUBROUTINE ION_CONDUCT (T,sal,sigma)

! Purpose :
!  Calculate ionic conductivity for saline water
!  following Stogryn

! local variables :
!  T : temperature of water (C)
!  sigma25 : ionic conductivity for sea water at 25 degree C (S/m)
!  sigma : ionic conductivity for saline water (S/m)
!---------------------------------------------------------------------------

USE PARKIND1, ONLY : JPIM, JPRM
IMPLICIT NONE

REAL(KIND = JPRM) :: T
REAL(KIND = JPRM) :: sal
REAL(KIND = JPRM) :: alphac
REAL(KIND = JPRM) :: delta
REAL(KIND = JPRM) :: sigma25
REAL(KIND = JPRM) :: sigma
!---------------------------------------------------------------------------

delta  =  25._JPRM - T

alphac  =  2.033e-2_JPRM  +  1.266e-4_JPRM * delta  +  2.464e-6_JPRM * delta ** 2.0_JPRM   &
          - sal * (1.849e-5_JPRM  -  2.551e-7_JPRM * delta  +  2.551e-8_JPRM * delta ** 2.0_JPRM)

sigma25  =  sal * (0.182521_JPRM  -  1.46192e-3_JPRM * sal        &
                             +  2.09324e-5_JPRM * sal ** 2.0_JPRM   &
                             -  1.28205e-7_JPRM * sal ** 3.0_JPRM)

sigma   = sigma25 * exp(-1._JPRM * delta * alphac)

END SUBROUTINE ION_CONDUCT
