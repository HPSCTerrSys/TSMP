! Copyright 2006-2014 ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities granted to it by
! virtue of its status as an intergovernmental organisation nor does it submit to any jurisdiction.
!
! This subroutine is created by Shaoning Lv
! it imitate CMEM ASCii reading subroutine to get data from state_p
! state_p is a vector in length of 2*nlay_soil_ls+18
SUBROUTINE WRCMEMSTATE_P

! Purpose :
! -------
!  WRITE CMEM OUTPUT for state_p vector for PDAF
   
! Interface :
! ---------

! Method :
! ------

! Externals :
! ---------

! Internal variables
! ------------------

! Author :
!     Patricia de Rosnay  October 2007

!------------------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM 
USE YOMLUN   , ONLY : NULOUT,NULTMP

USE YOMCMEMPAR, ONLY : JPHISTLEV,CNAMEID,CFREQ,CANGLE
USE YOMCMEMFIELDS, ONLY: N,JJ,doy,ndate,ntime, nbpt,CLNAME &
                     & ,fwc_veg,fb,ftfrac,fh,fWP,falpha,fsal_wat,ftau_atm,ftb_au &
                     & ,ftb_toa,fteffC,fsurf_emis,ftau_veg


IMPLICIT NONE

!------------------------------------------------------------------------------
! LSN: write output in memory
REAL, INTENT(out) :: m_state_p ! PE-local observed state

! LEVEL 1 outputs
!----------------

!  CLNAME='out_tbtoat_'//CNAMEID//'_'//cfreq//'_'//cangle//'.asc'

!  OPEN(NULTMP,FILE=CLNAME,status='replace')
!  WRITE(NULTMP,'(a72)') '#   nbpt   date        hour     doy    TBH(K)  TBV(K)   TEFF_H(K)   TEFF_V(k)'
!  DO JJ = 1, N
!  WRITE(NULTMP,'(3I9,f9.2,5(f9.3))') nbpt(JJ), ndate(JJ),ntime(JJ), doy(JJ), ftb_toa(JJ,1),ftb_toa(JJ,2),fteffC(JJ,2),fteffC(JJ,3)
!  ENDDO
!  CLOSE(NULTMP)

! LSN: h-polarization or v-polarization
m_state_p = ftb_toa(JJ,1)
! m_state_p = ftb_toa(JJ,2)


       
END SUBROUTINE WRCMEMSTATE_P

