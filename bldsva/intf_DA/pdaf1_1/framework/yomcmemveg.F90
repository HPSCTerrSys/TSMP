! Copyright 2006-2014 ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities granted to it by
! virtue of its status as an intergovernmental organisation nor does it submit to any jurisdiction.
!
MODULE YOMCMEMVEG

! Module containing the variables for the vegetation model for a single cell

!   t_veg   : vegetation temperature (K)
!   wc_veg  : vegetation water content (kg/m2)
!   sal_vw  : salinity of vegetation water (psu = ppt (weight))
!   bj       : b-parameter Jackson (tauN=bj*wc_veg)
!   w_eff   : single scattering albedo of vegetation (index 1.h-pol, 2.v-pol)
!             ATBD (0,0) ,  C-band (0.06)
!   w_effL   : single scattering albedo of LOW vegetation 
!   w_effH   : single scattering albedo of HIGH vegetation 
!   a_geo   : vegetation structure coefficient (index 1.h-pol, 2.v-pol)
!   a_geoL   : LOW vegetation structure coefficient (index 1.h-pol, 2.v-pol)
!   a_geoH   : HIGH vegetation structure coefficient (index 1.h-pol, 2.v-pol)
!   tau_veg : effective vegetation opacity (index 1.h-pol, 2.v-pol)
!   tb_veg  : microwave emission from vegetation (K) (index 1.h-pol, 2.v-pol)
! tth, ttv  : empirical parameter to account for incidence angle in tau at horiz. and vertical pol.

!   eps_vw  : dielectric constant of vegetation water
!---------------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM,JPRM
IMPLICIT NONE

REAL(KIND = JPRM) :: t_veg
REAL(KIND = JPRM) :: wc_veg
REAL(KIND = JPRM) :: sal_vw
REAL(KIND = JPRM) :: bj
REAL(KIND = JPRM) :: w_eff(2)
REAL(KIND = JPRM) :: w_effL(2)
REAL(KIND = JPRM) :: w_effH(2)
REAL(KIND = JPRM) :: a_geo(2)
REAL(KIND = JPRM) :: a_geoL(2)
REAL(KIND = JPRM) :: a_geoH(2)
REAL(KIND = JPRM) :: tau_veg(2)
REAL(KIND = JPRM) :: tauN,ttv,tth
REAL(KIND = JPRM) :: tb_veg(2)

COMPLEX(KIND = JPRM) :: eps_vw

END MODULE YOMCMEMVEG
