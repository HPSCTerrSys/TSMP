! Copyright 2006-2014 ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities granted to it by
! virtue of its status as an intergovernmental organisation nor does it submit to any jurisdiction.
!
MODULE YOMCMEMSOIL

! Module containing the variables for the soil emissivity model for a single cell

! Soil variables :
!   tb_soil : soil brightness temperature (K), Index 1 = hor. pol., 2 = ver. pol.
!   tsoil : temperature of each soil layer (K)
!   tsoildeep : temperature of deep soil layer (10-35cm?) (K)
!   tc : soil temperature in Celsius (C)
!   zsoil :  thickness of the soil layers (cm)
!   z_lsm: depth of the land surface model layers (m)
!   wc1 :  volumetric soil water content of surface (cm3/cm3)
!   wcsoil :  volumetric soil water content of profile (cm3/cm3)
!   wc :  volumetric soil water content of layer (cm3/cm3)
!   t_eff : effective soil temperature (K), Index 1 for C param, 2 or 3 for H & V teff
!    (put in fteffC for mean grib box average diagnostic)
!   surf_emis : surface emissivity  (-), Index 1 = hor. pol., 2 = ver. pol.
!   r_s : reflectivity of smooth soil surface (-), index 1 = h-pol., 2 = v-pol.
!   r_r : reflectivity of rough soil surface (-),  index 1 = h-pol., 2 = v-pol.

! Snow variables :
!   XMV : snow moisture [cm3/cm3]

! Soil Texture Parameters:
!   p : porosity of dry soil (cm3/cm3)
!   rho_b : soil bulk density (g/cm3) ~1.3,  f(soiltexture) 
!   rho_s : soil specific density (g/cm3), ATBD 2.664
!   sand  : sand fraction (% of dry weight of soil)
!   clay  : clay fraction (% of dry weight of soil)
!   wp : wilting point (cm3/cm3)
!   alpha : fitting parameter (Wang and Schmugge)

!   sal_soil : soil water salinity (psu = ppt(weight) ) , ATBD 0, (~0.65)
!   sal_sea : sea water salinity (psu = ppt(weight) ) ,  32.5, range 6-60psu
!   sal_wat : water salinity (psu = ppt(weight) )   
!   Q : Cross polarization parameter [-]
!   Nrv : exponent in surface roughness expression
!   Nrh

! parameters for effective temperature :
!   w0 : teff wign parameter (f(texture)) avignon (0.3) smosrex (0.33) ATBD (0.3)
!   bw : teff wign parameter (f(texture)) avignon (0.3) smosrex (0.63) ATBD (0.3)
!   eps0 : teff holm parameter (f(texture))  smosrex (0.08) 
!   bh : teff holm parameter (f(texture)) smosrex (0.87) 
            
! Dielectric Properties :
!   eps : Dielectric constant of the soil layer
!   eps_soil : Dielectric constant for each soil layer
!   eps_winf : dielectric constant at infinite frequency (Stogryn) 
!       also called: high frequency dielectric constant (4.9 Lane and Saxton )
!   eps_0 : permittivity of free space (Klein and Swift) [Farads/meter] 
!   ea : dielectric constant of air
!   er : dielectric constant of rock
!  ef :  dielectric constant of frozen soil 

! Control variables :
!   layer : current soil layer of computation
!---------------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM,JPRM
IMPLICIT NONE

REAL(KIND = JPRM) :: tb_soil(2) 
REAL(KIND = JPRM), allocatable :: tsoil(:) 
REAL(KIND = JPRM) :: tsoildeep
REAL(KIND = JPRM) :: tc 
REAL(KIND = JPRM), allocatable :: zsoil(:)
REAL(KIND = JPRM), allocatable :: z_lsm(:) 
REAL(KIND = JPRM) :: wc1 
REAL(KIND = JPRM), allocatable  :: wcsoil(:) 
REAL(KIND = JPRM) :: wc 
REAL(KIND = JPRM) :: t_eff(3)  
REAL(KIND = JPRM) :: surf_emis(2)
REAL(KIND = JPRM) :: r_r(2)
REAL(KIND = JPRM) :: r_s(2)

REAL(KIND = JPRM) :: XMV

REAL(KIND = JPRM) :: p
REAL(KIND = JPRM) :: rho_b
REAL(KIND = JPRM) :: rho_s
REAL(KIND = JPRM) :: sand
REAL(KIND = JPRM) :: clay
REAL(KIND = JPRM) :: wp
REAL(KIND = JPRM) :: alpha

REAL(KIND = JPRM) :: sal_soil
REAL(KIND = JPRM),allocatable :: sal_sea(:)
REAL(KIND = JPRM) :: sal_wat
REAL(KIND = JPRM) :: Nrv
REAL(KIND = JPRM) :: Nrh
REAL(KIND = JPRM) :: hrmodel

REAL(KIND = JPRM) :: w0
REAL(KIND = JPRM) :: bw 
REAL(KIND = JPRM) :: eps0
REAL(KIND = JPRM) :: bh

REAL(KIND = JPRM) :: eps_winf
REAL(KIND = JPRM) :: eps_0
COMPLEX :: eps
COMPLEX, allocatable :: eps_soil(:)
COMPLEX :: ea
COMPLEX :: er
COMPLEX :: ef

INTEGER :: JLAYER
!---------------------------------------------------------------------------

END MODULE YOMCMEMSOIL
