! Copyright 2006-2014 ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities granted to it by
! virtue of its status as an intergovernmental organisation nor does it submit to any jurisdiction.
!
!---------------------------------------------------------------------------
! Model to calculate the dielectric constant ice
!---------------------------------------------------------------------------

SUBROUTINE DIEL_ICE (T,eps_ice)

! Purpose :
!   Calculate the relative permittivty of pure ice in the microwave region

!       EPX=3.15,  EPY=-1e-3,  EPSIC=CMPLX(EPX,EPY) (ULABY)    

! Input/Output Variables:
! T (K)
! eps_ice : dielectric constant of ice
! rei : real part of dielectric constant of ice 
!---------------------------------------------------------------------------

USE YOMCMEMPAR, ONLY : fghz, tfreeze, LGPRINT

IMPLICIT NONE

REAL :: T
REAL :: B1, B2, BB, DELTABETA, BETAM, BETA,THETA,ALFA
REAL :: eir
COMPLEX :: eps_ice, E
!---------------------------------------------------------------------------

E=(0.,1.)

! -------- IMPROVED VERSION 
B1 = 0.0207
B2 = 1.16e-11
BB = 335.
DELTABETA = EXP(-10.02 + 0.0364*(T-tfreeze))
BETAM = (B1/T) * ( EXP(BB/T) / ((EXP(BB/T)-1.)**2.) ) + B2*(fghz**2.)
BETA = BETAM + DELTABETA
! -------- PURE ICE AS A FUNCTION OF TEMPERATURE 
THETA = 300. / T - 1.
ALFA = (0.00504 + 0.0062 * THETA) * EXP(-22.1 * THETA)

! -------- REAL PART OF ICE EPSILON 
eir = 3.1884 + 9.1e-4 * (T-tfreeze)

eps_ice = eir - E*(ALFA/fghz + BETA *fghz)
   
END SUBROUTINE DIEL_ICE
