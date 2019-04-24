! Copyright 2006-2014 ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities granted to it by
! virtue of its status as an intergovernmental organisation nor does it submit to any jurisdiction.
!
MODULE YOMCMEMATM

! Module containing the variables for the atmosphere and vegetation for a single cell

! Atmosphere variables :
!   tb_toa :: top-of-atmosphere brightness temperature (K)
!   tb_tov :: top-of-atmosphere brightness temperature (K)
!   tb_au :: upward atmospheric radiation
!   tb_ad :: downward atmospheric radiation
!   tau_atm :: optical thickness of atmosphere (zenith opacity / costheta)
!   tair : 2m air temperature (K)
!   Z :  Average altitude above sea level for the pixel (km)

! t_sky : cosmic background radiation (K)
!---------------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM,JPRM

IMPLICIT NONE

REAL(KIND = JPRM) :: tb_toa(2)
REAL(KIND = JPRM) :: tb_tov(2)
REAL(KIND = JPRM) :: tb_au
REAL(KIND = JPRM) :: tb_ad
REAL(KIND = JPRM) :: tau_atm
REAL(KIND = JPRM) :: tair
REAL(KIND = JPRM) :: Z

REAL(KIND = JPRM), allocatable, dimension (:) :: fZ 

REAL(KIND = JPRM) :: t_sky
!---------------------------------------------------------------------------

END MODULE YOMCMEMATM
