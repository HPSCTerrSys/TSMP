! Copyright 2006-2014 ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities granted to it by
! virtue of its status as an intergovernmental organisation nor does it submit to any jurisdiction.
!
MODULE YOMLUN


USE PARKIND1  ,ONLY : JPIM

USE YOMLUN_IFSAUX, ONLY : NULOUT, NULERR


IMPLICIT NONE


SAVE


!     ------------------------------------------------------------------


!*    Logical units used by code


!     NULOUT :   output unit (now from YOMLUN_IFSAUX)
!     NULNAM :   unit number for namelist
!     NCMAFL :   UNIT NUMBERS FOR OCEAN COUPLER FILES
!     NULCL1 :   unit number for climatological fields (month before)
!     NULCL2 :   unit number for climatological fields (month after)
!     NTRJSH :   unit number for trajectory spectral data          WRTRA
!     NINMSH :   unit number for initial point of the minimization SUVAZX
!     NINISH :   unit number for initial spectral data             SUSPEC
!     NINIGG :   unit number for initial grid-point data           SUSPEC
!     NFGISH :   unit number for first-guess spectral data
!     NFGIGG :   unit number for first-guess grid-point data
!     NPOSSH :   output unit number (spectral fields)              CPREP1
!     NTIDE  :   unit number for the LFI file containing the total tendencies
!     NPDIRL :   unit number for post-processing directory listing
!     NPPPSH :   unit number for post-processed spherical harmonics WRPLPP
!     NULDILA:   unit number for dilatation matrix (SUDIL,DILAT,SPDILA)
!     NULCONT:   unit number for contraction matrix (SUDIL,DILAT,SPDILA)
!     NULROTC:   unit number for upper troncature rotation matrix (SUROT,SPORTS)
!     NULCO  :   unit number for coupled fields (ICMCO)
!     NPODDH :   unit number for mask diagnostic files (DDH)
!     NULRCF :   unit number for restart control file
!     NULHWF :   unit number for history witness file
!     NBIAS  :   unit number for bias (dig. filt. guess - guess)
!     NEFLS  :   unit number for coupling ALADIN file
!     NEFLSS :   unit number for coupling ALADIN file (initialisation)
!     NULUSR1:   unit numbers for user defined files
!     NULSTAT:   unit number for status  file
!     NULASE :   unit number for CANARI statistics (forecast error s.d.)
!     NULASS :   unit number for CANARI statistics (analysis error s.d.)
!     NULUSR2
!     NULUSR3
!     NULUSR4
!     NULUSR5
!     NULTMP :   unit numbers for file opened and closed in the same routine
!     NULFPxx    unit numbers for Full-POS output files
!     NSCRTCH:   unit number for Full-POS scratch file (for in-line post-proc.)
!     NULFPOS:   unit number for Full-POS control file (end of post-processing
!                in conf. 001 ; auxilary namelist file in conf. 927)
!     NULERR :   unit number for comparison with reference run (now from YOMLUN_IFSAUX)
!     NULREF :   unit number for storing reference run
!     NULRAD :   unit number for writing radiation diagnostics
!     NULRTL :   unit number for reading RTLIMB coefficient files
!     NUO3CH1:   unit number for reading ozone chemistry file 1
!     NUO3CH2:   unit number for reading ozone chemistry file 2
!     NTCSR  :   unit number for fields of radiation coefficients
!     NSCATAB    SCAT. SIGMA0 TABLE
!     NSCASIG    SCAT. SIGMA0 BIAS CORRECTION
!     NSCASPE    SCAT. SPEED BIAS CORRECTION
!     NEGASH     UNIT NUMBER FOR JK INPUT
!     NULTRAJHR: unit number for high resolution trajectory (option LTRAJHR)
!     NULTRAJBG: unit number for background (option LBACKGR)
!!INTEGER(KIND=JPIM) :: NULOUT
INTEGER(KIND=JPIM) :: NULNAM
INTEGER(KIND=JPIM) :: NCMAFL(10)
INTEGER(KIND=JPIM) :: NPOSSH
INTEGER(KIND=JPIM) :: NTIDE
INTEGER(KIND=JPIM) :: NTRJSH
INTEGER(KIND=JPIM) :: NINMSH
INTEGER(KIND=JPIM) :: NINISH
INTEGER(KIND=JPIM) :: NINIGG
INTEGER(KIND=JPIM) :: NFGISH
INTEGER(KIND=JPIM) :: NFGIGG
INTEGER(KIND=JPIM) :: NPPPSH
INTEGER(KIND=JPIM) :: NULTMP
INTEGER(KIND=JPIM) :: NPODDH
INTEGER(KIND=JPIM) :: NULCL1
INTEGER(KIND=JPIM) :: NULCL2
INTEGER(KIND=JPIM) :: NULASE
INTEGER(KIND=JPIM) :: NULASS
INTEGER(KIND=JPIM) :: NULDILA
INTEGER(KIND=JPIM) :: NULCONT
INTEGER(KIND=JPIM) :: NULROTC
INTEGER(KIND=JPIM) :: NULRCF
INTEGER(KIND=JPIM) :: NULHWF
INTEGER(KIND=JPIM) :: NULUSR1
INTEGER(KIND=JPIM) :: NULUSR2
INTEGER(KIND=JPIM) :: NULUSR3
INTEGER(KIND=JPIM) :: NULUSR4
INTEGER(KIND=JPIM) :: NULUSR5
INTEGER(KIND=JPIM) :: NULCO
INTEGER(KIND=JPIM) :: NEFLS
INTEGER(KIND=JPIM) :: NEFLSS
INTEGER(KIND=JPIM) :: NBIAS
INTEGER(KIND=JPIM) :: NPDIRL
INTEGER(KIND=JPIM) :: NULSTAT
INTEGER(KIND=JPIM) :: NULFP01
INTEGER(KIND=JPIM) :: NULFP02
INTEGER(KIND=JPIM) :: NULFP03
INTEGER(KIND=JPIM) :: NULFP04
INTEGER(KIND=JPIM) :: NULFP05
INTEGER(KIND=JPIM) :: NULFP06
INTEGER(KIND=JPIM) :: NULFP07
INTEGER(KIND=JPIM) :: NULFP08
INTEGER(KIND=JPIM) :: NULFP09
INTEGER(KIND=JPIM) :: NULFP10
INTEGER(KIND=JPIM) :: NULFP11
INTEGER(KIND=JPIM) :: NULFP12
INTEGER(KIND=JPIM) :: NULFP13
INTEGER(KIND=JPIM) :: NULFP14
INTEGER(KIND=JPIM) :: NULFP15
INTEGER(KIND=JPIM) :: NULFPOS
INTEGER(KIND=JPIM) :: NSCRTCH
!!INTEGER(KIND=JPIM) :: NULERR
INTEGER(KIND=JPIM) :: NULREF
INTEGER(KIND=JPIM) :: NULRAD
INTEGER(KIND=JPIM) :: NULRTL
INTEGER(KIND=JPIM) :: NUO3CH1
INTEGER(KIND=JPIM) :: NUO3CH2
INTEGER(KIND=JPIM) :: NTCSR
INTEGER(KIND=JPIM) :: NSCATAB
INTEGER(KIND=JPIM) :: NSCASIG
INTEGER(KIND=JPIM) :: NSCASPE
INTEGER(KIND=JPIM) :: NEGASH
INTEGER(KIND=JPIM) :: NULTRAJHR
INTEGER(KIND=JPIM) :: NULTRAJBG


!     ------------------------------------------------------------------


END MODULE YOMLUN

