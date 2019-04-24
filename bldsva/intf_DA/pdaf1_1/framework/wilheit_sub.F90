! Copyright 2006-2014 ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities granted to it by
! virtue of its status as an intergovernmental organisation nor does it submit to any jurisdiction.
!
SUBROUTINE WILHEIT

! Purpose :
!  CALCULATE REFLECTION AND EMISSION OF SOIL 
!  CMEM v4.0 Change in t_eff dimensions (1 for C, 2 for TEFF H, 3 for TEFF V) 
!            to ensure a correct output of teff
!  CMEM v4.1 P de Rosnay May 2012: Improved Wilheit: flexible LSM and MW model vertical grid.
!                              interpolation from LSM to MW model vbertical grid
!---------------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM,JPRM
USE YOMCMEMPAR, ONLY : theta , nlay_soil_mw, lamcm, c, pi, LGPRINT
USE YOMCMEMSOIL, ONLY : eps_soil , r_s , tsoil , zsoil , t_eff 

IMPLICIT NONE

REAL(KIND = JPRM) :: thrad
REAL(KIND = JPRM):: zsoil_wil(nlay_soil_mw+1), fa(nlay_soil_mw)
COMPLEX(KIND = JPRM):: nref(nlay_soil_mw+1), cp(nlay_soil_mw+1) 
INTEGER(KIND=JPIM) :: i
REAL(KIND = JPRM):: z
REAL(KIND = JPRM):: sumtfa, sumfa, sumdz, sumzfa
!---------------------------------------------------------------------------

! Conversions

thrad = theta * pi / 180.   ! incidence angle (degrees to radians)

! Wilheit-model works with atmosphere as an additional layer.
! That means that the first soil layer is layer 2 in the model, ...

zsoil_wil(1) = 1.         ! thickness of atmosphere (arbitrary value!)
do i = 1, nlay_soil_mw
  zsoil_wil(i+1) = zsoil(i)
enddo

! Calculate index of refraction

nref(1) = (1.0,0.0)   ! atmosphere
do i = 1, nlay_soil_mw
  nref(i+1) = sqrt(eps_soil(i))
enddo

! Calculate reflection and effective soil temperature

call hpld (thrad, nlay_soil_mw+1, nref, zsoil_wil, cp)

r_s(1) = real(cp(1))

sumtfa = 0.
sumdz = 0.
sumfa = 0.
sumzfa = 0.
do i = 1, nlay_soil_mw
  fa(i) = real(cp(i+1))
  sumtfa = sumtfa + fa(i) * tsoil(i)
  sumfa = sumfa + fa(i)
  sumdz = sumdz + zsoil(i)
  z = sumdz - zsoil(i)/2.              ! mid-depth of actual layer
  sumzfa = sumzfa + z * fa(i)
enddo
t_eff(2) = sumtfa / sumfa

! Calculate reflection and emission for vertical polarization

if (theta /= 0.0) then

  call vpld (thrad, nlay_soil_mw+1, nref, zsoil_wil, cp)

  r_s(2) = real(cp(1))

  sumtfa = 0.
  sumdz = 0.
  sumfa = 0.
  sumzfa = 0.
  do i = 1, nlay_soil_mw
    fa(i) = real(cp(i+1))
    sumtfa = sumtfa + fa(i) * tsoil(i)
    sumfa = sumfa + fa(i)
    sumdz = sumdz + zsoil(i)
    z = sumdz - zsoil(i)/2.              ! mid-depth of actual layer
    sumzfa = sumzfa + z * fa(i)
  enddo
  t_eff(3) = sumtfa / sumfa

else

  r_s(2) = r_s(1)

endif

END SUBROUTINE WILHEIT

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------

SUBROUTINE HPLD (THRAD, N, CN, DEL, CP)

! Purpose : 
!-- CALCULATE SURFACE REFLECTIVITY AND LAYER ABSORPTION FRACTIONS ----------
!----------- AT HORIZONTAL POLARIZATION ------------------------------------
 
! Input/Output Variables:
!   Name    Type     I/O  Description
!   THRAD   real      I   incidence angle [rad]
!   N       integer   I   number of layers (= number of soil layers plus 1,
!                         layer 1 = atmosphere)
!   CN      complex   I   index of refraction for atmosphere (CN(1)) 
!                         and each soil layer
!   DEL     real      I   layer thicknesses [cm]
!                         (the first and last layer are assumed to be 
!                         semi-infinite, their thickness value is ignored)
!   CP      complex   O   At the end of this subroutine, the real part
!                         of CP(1) contains the soil-surface reflectivity
!                         and the real parts of CP(2)..CP(N) contain the
!                         absorption fraction of each layer.
!---------------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM,JPRM
USE YOMCMEMPAR ,ONLY : lamcm , pi, LGPRINT

IMPLICIT NONE

INTEGER(KIND=JPIM) :: N, NL, NMAX
REAL(KIND = JPRM) :: THRAD
REAL(KIND = JPRM) :: DEL(N)

INTEGER(KIND=JPIM) :: I, J, JJ, LL
REAL(KIND = JPRM) :: R, S, ARG, E2, DP, X
COMPLEX(KIND = JPRM*2) :: CS, CC, CARG
COMPLEX(KIND = JPRM*2) :: CSJ, CCJ, CSJP1, CCJP1, CA, CB, CX, CXP
COMPLEX(KIND = JPRM*2) :: CN(N), CP(N), CEP(N), CEM(N)
!---------------------------------------------------------------------------

! Calculate squares of interface propagators.
! CP(I) here is equal to p(j)**2 in section II of the paper.

S = SIN(THRAD)
CP(1) = (1.,0.)
NL = N - 1
DO I = 2, NL
  NMAX = I + 1   ! paper says NMAX=I, but this seems to be wrong
  CS = CN(1) * S / CN(I)
  CC = CSQRT( (1.,0.) - CS*CS )  ! in paper CC is cos(theta)(j)
  ARG = DEL(I) * 2. * pi / lamcm
  CARG = 2. * ARG * CN(I) * CC * (0.,1.)
  CP(I) = CEXP(CARG) * CP(I-1)
  IF (CABS(CP(I)) .LT. 0.0001) EXIT
ENDDO

! Calculate electric fields within each layer.
! CEP(J) and CEM(J) here are equal to E+(j) and E-(j) in the paper.

CEP(NMAX) = (1.,0.)
CEM(NMAX) = (0.,0.)
DO JJ = 2, NMAX
  J = NMAX - JJ + 1
  CSJ = CN(1) * S / CN(J)
  CCJ = CSQRT( (1.,0.) - CSJ*CSJ )
  CSJP1 = CN(1) * S / CN(J+1)
  CCJP1 = CSQRT( (1.,0.) - CSJP1*CSJP1 )
  CA = 2. * CN(J) * CCJ / ( CN(J)*CCJ + CN(J+1)*CCJP1 )
  CB = ( CN(J)*CCJ - CN(J+1)*CCJP1 ) / ( (CN(J)*CCJ + CN(J+1)*CCJP1) * CP(J) )
  CEP(J) = CEP(J+1)/CA + CB*CEM(J+1)/CA
  CEM(J) = CEM(J+1) + (CEP(J+1)-CEP(J))*CP(J)
ENDDO
CX = CEP(1)
DO J = 1, NMAX
  CEP(J) = CEP(J) / CX
  CEM(J) = CEM(J) / CX
ENDDO

! Calculate absorption fraction of each layer (f(i) in paper).

DO J = NMAX, N
  CP(J) = (0.,1.) / (1.E15)
ENDDO
LL = NMAX - 1
DO JJ = 1, LL
  J = NMAX - JJ + 1
  CS = SIN(THRAD) / CN(J)
  CC = CSQRT((1.,0.)-CS*CS)
  R = CABS(CP(J))
  S = CABS(CP(J-1))
  E2 = (S-R) * CABS(CEP(J))**2. + (1./R-1./S) * CABS(CEM(J))**2.
  DP = E2 * REAL(CN(J)*CC) / COS(THRAD)
  CXP = CEP(J) * CONJG(CEM(J))
  X = 2. * AIMAG(CN(J)*CC/COS(THRAD))  &
         * ( AIMAG(CXP*CP(J-1)/CABS(CP(J-1))) - AIMAG(CXP*CP(J)/CABS(CP(J))) )
  DP = DP - X
  CP(J) = CMPLX(DP,0.)
ENDDO
R = CABS(CEM(1))**2. * REAL(CN(1))
CP(1) = CMPLX(R,0.)

END SUBROUTINE HPLD

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------

SUBROUTINE VPLD (THRAD, N, CN, DEL, CP)

! Purpose : 
!-- CALCULATE SURFACE REFLECTIVITY AND LAYER ABSORPTION FRACTIONS ----------
!----------- AT VERTICAL POLARIZATION ------------------------------------

! Input/Output Variables:
!   Name    Type     I/O  Description
!   THRAD   real      I   incidence angle [rad]
!   N       integer   I   number of layers (= number of soil layers plus 1,
!                         layer 1 = atmosphere)
!   CN      complex   I   index of refraction for atmosphere (CN(1)) 
!                         and each soil layer
!   DEL     real      I   layer thicknesses [cm]
!                         (the first and last layer are assumed to be 
!                         semi-infinite, their thickness value is ignored)
!   CP      complex   O   At the end of this subroutine, the real part
!                         of CP(1) contains the soil-surface reflectivity
!                         and the real parts of CP(2)..CP(N) contain the
!                         absorption fraction of each layer.
!---------------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM,JPRM
USE YOMCMEMPAR, ONLY : lamcm , pi, LGPRINT

IMPLICIT NONE

INTEGER(KIND=JPIM) :: N, NL, NMAX
REAL(KIND = JPRM) :: THRAD
REAL(KIND = JPRM) :: DEL(N)

INTEGER(KIND=JPIM) :: I, J, JJ, LL
REAL(KIND = JPRM) :: R, S, ARG, E2, DP, X
COMPLEX(KIND = JPRM*2) :: CS, CC, CARG
COMPLEX(KIND = JPRM*2) :: CSJ, CCJ, CSJP1, CCJP1, CD, CA, CB, CR, CX, CXP
COMPLEX(KIND = JPRM*2) :: CN(N), CP(N), CEP(N), CEM(N)
!---------------------------------------------------------------------------

! Calculate squares of interface propagators.
! CP(I) here is equal to p(j)**2 in section II of the paper.

S = SIN(THRAD)
CP(1) = (1.,0.)
NL = N - 1
DO I = 2, NL
  NMAX = I + 1   ! paper says NMAX=I, but this seems to be wrong
  CS = CN(1) * S / CN(I)
  CC = CSQRT( (1.,0.) - CS*CS )  ! in paper CC is cos(theta)(j)
  ARG = DEL(I) * 2. * pi / lamcm
  CARG = 2. * ARG * CN(I) * CC * (0.,1.)
  CP(I) = CEXP(CARG) * CP(I-1)
  IF (CABS(CP(I)) .LT. 0.0001) EXIT
ENDDO

! Calculate electric fields within each layer.
! CEP(J) and CEM(J) here are equal to E+(j) and E-(j) in the paper.

CEP(NMAX) = (1.,0.)
CEM(NMAX) = (0.,0.)
DO JJ = 2, NMAX
  J = NMAX - JJ + 1
  CSJ = CN(1) * S / CN(J)
  CCJ = CSQRT( (1.,0.) - CSJ*CSJ )
  CSJP1 = CN(1) * S / CN(J+1)
  CCJP1 = CSQRT( (1.,0.) - CSJP1*CSJP1 )
  CD = 2. * CN(J) * CCJ
  CA = CN(J)*CCJP1 + CN(J+1)*CCJ
  CB = CN(J)*CCJP1 - CN(J+1)*CCJ
  CEP(J) = CA*CEP(J+1)/CD + CB*CEM(J+1)/(CD*CP(J))
  CR = CN(J+1) / CN(J)
  CEM(J) = CR*CEM(J+1) + (CEP(J)-CEP(J+1)*CR)*CP(J)
ENDDO
CX = CEP(1)
DO J = 1, NMAX
  CEP(J) = CEP(J) / CX
  CEM(J) = CEM(J) / CX
ENDDO

! Calculate absorption fraction of each layer (f(i) in paper).

DO J = NMAX, N
  CP(J) = (0.,1.) / (1.E15)
ENDDO
LL = NMAX - 1
DO JJ = 1, LL
  J = NMAX - JJ + 1
  CS = SIN(THRAD) / CN(J)
  CC = CSQRT((1.,0.)-CS*CS)
  R = CABS(CP(J))
  S = CABS(CP(J-1))
  E2 = (S-R) * CABS(CEP(J))**2. + (1./R-1./S) * CABS(CEM(J))**2.
  DP = E2 * REAL(CN(J)*CC) / COS(THRAD)
  CXP = CEP(J) * CONJG(CEM(J))
  X = 2. * AIMAG(CN(J)*CC/COS(THRAD))  &
          * ( AIMAG(CXP*CP(J-1)/CABS(CP(J-1))) - AIMAG(CXP*CP(J)/CABS(CP(J))) )
  DP = DP - X
  CP(J) = CMPLX(DP,0.)
ENDDO
R = CABS(CEM(1))**2. * REAL(CN(1))
CP(1) = CMPLX(R,0.)

END SUBROUTINE VPLD
