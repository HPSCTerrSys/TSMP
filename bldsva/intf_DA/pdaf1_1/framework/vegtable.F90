! Copyright 2006-2014 ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities granted to it by
! virtue of its status as an intergovernmental organisation nor does it submit to any jurisdiction.
!
SUBROUTINE VEGTABLE (ZRVLAI,ZRVCOV,ZRVTV,Zb,ZVWC,Zb1,Zb2,Zb3,ZNrh,ZNrv,Ztth,Zttv,Zhr,Zw_eff)

! Purpose :
! -------
!     Build tables of vegetation characteristics
! Convert TESSEL (ECMWF) Vegetation biome classes to ECOCLIMAP 7 biomes

! Drusch, M. et al 2009 JHM DOI: 10.1175/2008JHM964.1
! de Rosnay P. et al 2009 JGR doi:10.1029/2008JD010724
!------------------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRM
USE YOMLUN   , ONLY : NULOUT

IMPLICIT NONE

REAL(KIND=JPRM) :: Zw_eff(0:7,2)
REAL(KIND=JPRM) :: Zhr(0:7)
REAL(KIND=JPRM) :: Ztth(0:7)
REAL(KIND=JPRM) :: Zttv(0:7)
REAL(KIND=JPRM) :: ZNrh(0:7)
REAL(KIND=JPRM) :: ZNrv(0:7)
REAL(KIND=JPRM) :: Zb(0:7)
REAL(KIND=JPRM) :: Zb1(0:7)
REAL(KIND=JPRM) :: Zb2(0:7)
REAL(KIND=JPRM) :: Zb3(0:7)
REAL(KIND=JPRM) :: ZVWC(0:7)
REAL(KIND=JPRM) :: ZRVLAI(0:20)
REAL(KIND=JPRM) :: ZRVCOV(0:20)
INTEGER(KIND=JPIM) :: ZRVTV(0:20)
!------------------------------------------------------------------------------

! 1.  LEAF AREA INDEX
!     ---------------
ZRVLAI(1)=3.     ! Crops, Mixed Farming
ZRVLAI(2)=2.     ! Short Grass
ZRVLAI(3)=5.     ! Evergreen Needleleaf Trees
ZRVLAI(4)=5.     ! Deciduous Needleleaf Trees
ZRVLAI(5)=5.     ! Deciduous Broadleaf Trees
ZRVLAI(6)=6.     ! Evergreen Broadleaf Trees
ZRVLAI(7)=2.     ! Tall Grass
ZRVLAI(8)=0.5      ! Desert
ZRVLAI(9)=1.0       ! Tundra
ZRVLAI(10)=3.    ! Irrigated Crops
ZRVLAI(11)=0.5     ! Semidesert
ZRVLAI(12)=0.0     ! Ice Caps and Glaciers
ZRVLAI(13)=4.    ! Bogs and Marshes
ZRVLAI(14)=0.0     ! Inland Water
ZRVLAI(15)=0.0     ! Ocean
ZRVLAI(16)=3.    ! Evergreen Shrubs
ZRVLAI(17)=1.5   ! Deciduous Shrubs
ZRVLAI(18)=5.    ! Mixed Forest/woodland
ZRVLAI(19)=2.5   ! Interrupted Forest
ZRVLAI(20)=4.    ! Water and Land Mixtures
ZRVLAI(0)=0.0      ! no continent

! 2.  VEGETATION COVER
!     ----------------


ZRVCOV(1)=0.9     ! Crops, Mixed Farming
ZRVCOV(2)=0.85    ! Short Grass
ZRVCOV(3)=0.9     ! Evergreen Needleleaf Trees
ZRVCOV(4)=0.9     ! Deciduous Needleleaf Trees
ZRVCOV(5)=0.9     ! Deciduous Broadleaf Trees
ZRVCOV(6)=0.99    ! Evergreen Broadleaf Trees
ZRVCOV(7)=0.7     ! Tall Grass
ZRVCOV(8)=0.       ! Desert
ZRVCOV(9)=0.5       ! Tundra
ZRVCOV(10)=0.9    ! Irrigated Crops
ZRVCOV(11)=0.1    ! Semidesert
ZRVCOV(12)=0.      ! Ice Caps and Glaciers
ZRVCOV(13)=0.6    ! Bogs and Marshes
ZRVCOV(14)=0.      ! Inland Water
ZRVCOV(15)=0.      ! Ocean
ZRVCOV(16)=0.5      ! Evergreen Shrubs
ZRVCOV(17)=0.5      ! Deciduous Shrubs
ZRVCOV(18)=0.9    ! Mixed Forest/woodland
ZRVCOV(19)=0.9    ! Interrupted Forest
ZRVCOV(20)=0.6    ! Water and Land Mixtures
ZRVCOV(0)=0.0

! 3.  Corresponding biome in 7 class ECOCLIMAP + 1
!     ----------------------------------------
ZRVTV(1)=6    ! Crops, Mixed Farming
ZRVTV(2)=4    ! Short Grass
ZRVTV(3)=2    ! Evergreen Needleleaf Trees
ZRVTV(4)=2    ! Deciduous Needleleaf Trees
ZRVTV(5)=1    ! Deciduous Broadleaf Trees
ZRVTV(6)=3    ! Evergreen Broadleaf Trees
ZRVTV(7)=4    ! Tall Grass
ZRVTV(8)=0    ! Desert
ZRVTV(9)=4    ! Tundra
ZRVTV(10)=6   ! Irrigated Crops
ZRVTV(11)=4   ! Semidesert
ZRVTV(12)=0   ! Ice Caps and Glaciers
ZRVTV(13)=4   ! Bogs and Marshes
ZRVTV(14)=0   ! Inland Water
ZRVTV(15)=0   ! Ocean
ZRVTV(16)=4   ! Evergreen Shrubs
ZRVTV(17)=4   ! Deciduous Shrubs
ZRVTV(18)=1   ! Mixed Forest/woodland
ZRVTV(19)=1   ! Interrupted Forest
ZRVTV(20)=0   ! Water and Land Mixtures
ZRVTV(0)=0    ! Empty

! 4  vegetation b parameter for L-band (used by Jackson vegetation model)
!    -----------------------------------------
Zb(0)=0.0      ! no veg
Zb(1)=0.33     ! Decidious forests
Zb(2)=0.33     ! Coniferous forests
Zb(3)=0.33     ! Rain forests
Zb(4)=0.20     ! C3 Grasslands
Zb(5)=0.20     ! C4 Grasslands
Zb(6)=0.15     ! C3 Crops
Zb(7)=0.15     ! C4 Crops

! 5  Vegetation water content (Jackson)
!    ------------------------
ZVWC(0)=0.         ! No vegetation
ZVWC(1)=4.         ! Decidious forests (4)
ZVWC(2)=3.         ! Coniferous forests (3)
ZVWC(3)=10.        ! Rain forests (6)
ZVWC(4)=0.! f(LAI) ! Grasslands
ZVWC(5)=0.! f(LAI) ! Grasslands
ZVWC(6)=0.! f(LAI) ! Crops 
ZVWC(7)=0.! f(LAI) ! Crops 


! 6 Vegetation b1 and b2 parameter (Wigneron- replace b and VWC in Jackson)
!
Zb1(0)=0.0       ! no veg
Zb1(1)=0.226     ! Decidious forests
Zb1(2)=0.260     ! Coniferous forests
Zb1(3)=0.226     ! Rain forests
Zb1(4)=0.0375    ! C3 Grasslands
Zb1(5)=0.0375    ! C4 Grasslands
Zb1(6)=0.05      ! C3 Crops
Zb1(7)=0.05      ! C4 Crops

Zb2(0)=0.0       ! no veg
Zb2(1)=0.001     ! Decidious forests
Zb2(2)=0.006     ! Coniferous forests
Zb2(3)=0.001     ! Rain forests
Zb2(4)=0.05      ! C3 Grasslands
Zb2(5)=0.05      ! C4 Grasslands
Zb2(6)=0.0       ! C3 Crops
Zb2(7)=0.0       ! C4 Crops

Zb3(0)=0.0       ! no veg
Zb3(1)=0.7       ! Decidious forests
Zb3(2)=0.69      ! Coniferous forests
Zb3(3)=0.7       ! Rain forests
Zb3(4)=0.0       ! C3 Grasslands
Zb3(5)=0.0       ! C4 Grasslands
Zb3(6)=0.0       ! C3 Crops
Zb3(7)=0.0       ! C4 Crops


! 7 Soil roughness parameter Nr as function of vegetation type 
! in case Wigneron

ZNrh(0)=0.0       ! Bare soil
ZNrh(1)=1.0       ! Decidious forests
ZNrh(2)=1.75      ! Coniferous forests
ZNrh(3)=1.0       ! Rain forests
ZNrh(4)=1.0       ! C3 Grasslands
ZNrh(5)=1.0       ! C4 Grasslands
ZNrh(6)=0.0       ! C3 Crops
ZNrh(7)=0.0       ! C4 Crops

ZNrv(0)=-1.0       ! Bare soil
ZNrv(1)=2.0       ! Decidious forests
ZNrv(2)=0.0      ! Coniferous forests
ZNrv(3)=0.0       ! Rain forests
ZNrv(4)=0.0       ! C3 Grasslands
ZNrv(5)=0.0       ! C4 Grasslands
ZNrv(6)=-1.0       ! C3 Crops
ZNrv(7)=-1.0       ! C4 Crops

! 8 Empirical parameter to account for incidence angle in tau
!   Wigneron

Ztth(0)=1.0        ! Bare soil
Ztth(1)=0.49       ! Decidious forests
Ztth(2)=0.8        ! Coniferous forests
Ztth(3)=1.0        ! Rain forests
Ztth(4)=1.0        ! C3 Grasslands
Ztth(5)=1.0        ! C4 Grasslands
Ztth(6)=1.0       ! C3 Crops
Ztth(7)=2.0       ! C4 Crops

Zttv(0)=1.0       ! Bare soil
Zttv(1)=0.46       ! Decidious forests
Zttv(2)=0.8        ! Coniferous forests
Zttv(3)=1.0        ! Rain forests
Zttv(4)=1.0        ! C3 Grasslands
Zttv(5)=1.0        ! C4 Grasslands
Zttv(6)=2.0       ! C3 Crops
Zttv(7)=1.0       ! C4 Crops

! 9 Empirical roughness parameter for Wigneron

Zhr(0)=0.1        ! Bare soil
Zhr(1)=1.0        ! Decidious forests
Zhr(2)=1.2        ! Coniferous forests
Zhr(3)=1.3        ! Rain forests
Zhr(4)=0.0        ! C3 Grasslands computed in cmem_init as 1.3 - 1.13 sm
Zhr(5)=0.0        ! C4 Grasslands computed in cmem_init as 1.3 - 1.13 sm
Zhr(6)=0.1        ! C3 Crops
Zhr(7)=0.6        ! C4 Crops

! 10 Omega, horizontal and vertical polarisation

Zw_eff(0,:)=(0.0,0.0)          ! Bare soil
Zw_eff(1,:)=(0.07,0.07)        ! Decidious forests
Zw_eff(2,:)=(0.08,0.08)        ! Coniferous forests
Zw_eff(3,:)=(0.095,0.095)      ! Rain forests
Zw_eff(4,:)=(0.05,0.05)        ! C3 Grasslands 
Zw_eff(5,:)=(0.05,0.05)        ! C4 Grasslands
Zw_eff(6,:)=(0.0,0.0)          ! C3 Crops
Zw_eff(7,:)=(0.05,0.05)        ! C4 Crops



END SUBROUTINE VEGTABLE
