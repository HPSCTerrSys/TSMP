! Copyright 2006-2014 ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities granted to it by
! virtue of its status as an intergovernmental organisation nor does it submit to any jurisdiction.
!
SUBROUTINE CMEM_INIT

! Purpose :
! -------
! Compute Optical depths and soil roughness 
   
! Interface :
! ---------

! Method :
! ------
! Externals :
! ---------

! Internal variables
! ------------------
! Authors :
!     Thomas Holmes and Patricia de Rosnay  August 2007
!     Modifications: 
!     Patricia de Rosnay ECMWF  v1.1 Dec  2007 (Grib and ASCII options)
!                               v1.2 February 2008 (Netcdf option added)
!                               v1.3 April 2008
!                               v1.4 November 2008
!                               v2.0 December 2008 (Grib API option added: default option instead of grib)
!                               v3.0 November 2009 (bugfix Wsimple)
!                               v4.0 February 2012  HTessel option, bug fix in TEFF and VEGTABLE
!    Hans Lievens Ghent Univ.   v4.2 August   2012  Bug fix (patch from v4.1) in the porosity computation 

!------------------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM, JPRM
USE YOMLUN   , ONLY : NULOUT

USE YOMCMEMPAR, ONLY : CIDVEG,CIRGHR, CIATM, CIVEG, LGPRINT,kcm, ip_rgh_surf, tfreeze,  fghz 
USE YOMCMEMFIELDS, ONLY: N, JJ, ftl_lsm,ftair, fTVL, fTVH, fwc_veg, fb, fs_laiL, ftfrac &
          & , fwc_lsm,fsand,fclay,frho_b,fp,falpha,fwp,fsal_wat,fwater, fh,ftau_atm, ftb_au &
          & ,ftb_ad,ftauN,fNrh_L,fNrh_H,fNrv_L,fNrv_H,ftth,fttv,fhrmodel,fw_effL,fw_effH& 
          & ,fteffC
USE YOMCMEMSOIL, ONLY : rho_s, sal_sea, sal_soil


IMPLICIT NONE

INTEGER(KIND=JPIM) :: i
REAL(KIND=JPRM) :: rvcov(0:20)
REAL(KIND=JPRM) :: rvlai(0:20)
INTEGER(KIND=JPIM) :: RVTV(0:20)
REAL(KIND=JPRM) :: b(0:7)
REAL(KIND=JPRM) :: b1(0:7)
REAL(KIND=JPRM) :: b2(0:7)
REAL(KIND=JPRM) :: b3(0:7)
REAL(KIND=JPRM) :: VWC(0:7)

REAL(KIND=JPRM) :: ZNrh(0:7)
REAL(KIND=JPRM) :: ZNrv(0:7)
REAL(KIND=JPRM) :: Ztth(0:7)
REAL(KIND=JPRM) :: Zttv(0:7)
REAL(KIND=JPRM) :: Zhr(0:7)
REAL(KIND=JPRM) :: Zw_eff(0:7,2)

REAL(KIND=JPRM) :: hr_min
REAL(KIND=JPRM) :: hr
REAL(KIND=JPRM) :: fc(N)
REAL(KIND=JPRM) :: wt(N)
REAL(KIND=JPRM) :: xlc

!------------------------------------------------------------------------------
! 1. ATMOSPHERIC MODULE
! ---------------------------------------------------------------------------


CALL CMEM_ATM


! ---------------------------------------------------------------------------
! 2. Initialize vegetation parameters (tauN/wc_veg relation)
! ---------------------------------------------------------------------------


CALL VEGTABLE (RVLAI,RVCOV,RVTV,b,VWC,b1,b2,b3,ZNrh,ZNrv,Ztth,Zttv,Zhr,Zw_eff) 
    

  
! 2.3 Low and High Vegetation water content (Kirdyashev, Wegmuller, Jackson)

fwc_veg(:,1)= 0.5 * fs_laiL(:)
fwc_veg(:,2)= VWC(fTVH(:))

fw_effL(:,1) = 0.05
fw_effL(:,2) = 0.05
fw_effH(:,1) = 0.05
fw_effH(:,2) = 0.05
 
SELECT CASE (CIVEG)

   CASE ('Jackson')
   ! 2.4 b parameter for vegetation model (Jackson)
      
      fb(:,1)= b(fTVL(:))
      fb(:,2)= b(fTVH(:))       
  
   CASE ('Wigneron')
   ! 2.5 b1 and b2 parameters used to compute taunadir 
 
      ftauN(:,1) = b1(fTVL(:)) * fs_laiL(:) + b2(fTVL(:)) 
      ftauN(:,2) = b3(fTVH(:))  ! for high veg cste value assigned

   ! 2.6 tth and ttv function of vegetation type 

     ftth(:,1) = Ztth(fTVL(:))
     ftth(:,2) = Ztth(fTVH(:))

     fttv(:,1) = Zttv(fTVL(:))
     fttv(:,2) = Zttv(fTVH(:))
    
   ! 2.7  Wigneron single diffusion albedo
    
     fw_effL(:,1) = Zw_eff(fTVL(:),1)
     fw_effL(:,2) = Zw_eff(fTVL(:),2)
     fw_effH(:,1) = Zw_eff(fTVH(:),1)
     fw_effH(:,2) = Zw_eff(fTVH(:),2)

END SELECT 


! ---------------------------------------------------------------------------
! 3.  Initialize soil characteristics
! ---------------------------------------------------------------------------

! 3.0 rougness parameters (Wigneron)

! fNrh and fNrv read from vegtable and then transmitted to CMEM_main 

fNrh_L(:) = ZNrh(fTVL(:)) 
fNrh_H(:) = ZNrh(fTVH(:)) 

fNrv_L(:) = ZNrv(fTVL(:)) 
fNrv_H(:) = ZNrv(fTVH(:)) 


! 3.1 bulk density

frho_b(:)= (fsand(:)*1.6_JPRM +  fclay(:)*1.1_JPRM + (100._JPRM-fsand(:)-fclay(:))*1.2_JPRM)/100._JPRM
! for LMEB convergence
! frho_b(:)= 1.4


! 3.2 porosity

fp(:)= 1.0_JPRM - frho_b(:)/rho_s


! 3.3 Wilting point, for Wang dielectric model

fwp(:) = 0.06774_JPRM - 0.00064_JPRM * fsand(:) + 0.00478_JPRM * fclay(:)

! 3.4 alpha, for conductivity loss in Wang dielectric model

SELECT CASE ( fghz > 2.5_JPRM )
      CASE ( .TRUE. )
        falpha(:)=0._JPRM
      CASE ( .FALSE. )
        falpha(:)=100._JPRM * fWP(:)
        WHERE ( falpha(:) > 26._JPRM) falpha(:)=26._JPRM
END SELECT

! 3.5 initialize TEFF to zero
fteffC(:,:) = 0._JPRM


! ---------------------------------------------------------------------------
! 4. Water salinity
! ---------------------------------------------------------------------------

WHERE ( fwater(:) < 0.5 )  fsal_wat(:) = sal_soil
WHERE ( fwater(:) >= 0.5 )  fsal_wat(:) = sal_sea
!
!
! ---------------------------------------------------------------------------
! 5. Roughness factor H [cm] global variable, now constant
! ---------------------------------------------------------------------------
!
SELECT CASE (CIRGHR)
  CASE ( 'No' )
    fh(:) = 0._JPRM
    fhrmodel(:,:) = 0._JPRM
  CASE ( 'Choudhury' )
    fh(:) = (2._JPRM * kcm * ip_rgh_surf) ** 2._JPRM
    fhrmodel(:,1) = fh(:)  ! on low veg
    fhrmodel(:,2) = fh(:)  ! on high veg
     fNrh_L(:) = 0.0_JPRM 
     fNrh_H(:) = 0.0_JPRM 
     fNrv_L(:) = 0.0_JPRM 
     fNrv_H(:) = 0.0_JPRM 
  CASE ( 'Wsimple' )   !    Wigneron model (see NOTICE file)
    xlc = 6._JPRM
    fh(:) = 1.3972_JPRM* ( ip_rgh_surf / xlc ) ** 0.5879_JPRM
    fhrmodel(:,1) = fh(:)  ! on low veg
    fhrmodel(:,2) = fh(:)  ! on high veg
    fNrh_L(:) = 0.0_JPRM 
    fNrh_H(:) = 0.0_JPRM 
    fNrv_L(:) = 0.0_JPRM 
    fNrv_H(:) = 0.0_JPRM 
  CASE ( 'Wegmueller' )
    fh(:) = kcm * ip_rgh_surf
    fhrmodel(:,1) = fh(:)  ! on low veg
    fhrmodel(:,2) = fh(:)  ! on high veg
  CASE ( 'Wtexture' )  !    Wigneron model (see NOTICE file)
    !   hr : standard deviation of surface height (cm)
    hr_min = 0.05_JPRM
    hr = 0.10_JPRM
    wt(:) = 0.49_JPRM * fwp(:) + 0.165_JPRM
    fc(:) = wt + 0.10_JPRM * fclay(:) / 100._JPRM
    fh(:) = hr_min
    WHERE ( fwc_lsm(:,1) < fc(:) )  fh(:) = hr - (hr-hr_min)*((fwc_lsm(:,1)-wt(:))/(fc(:)-wt(:)))
    WHERE ( fwc_lsm(:,1) < wt(:) )  fh(:) = hr
    fhrmodel(:,1) = fh(:)  ! on low veg
    fhrmodel(:,2) = fh(:)  ! on high veg 
  CASE ( 'Wigneron' )   !    Wigneron (see NOTICE file)
     fhrmodel(:,1) = Zhr(fTVL(:))
     fhrmodel(:,2) = Zhr(fTVH(:))
     WHERE ( fTVL(:) == 4_JPIM )  fhrmodel(:,1) = 1.3_JPRM - 1.13_JPRM *  fwc_lsm(:,1)
     WHERE ( fTVL(:) == 5_JPIM )  fhrmodel(:,1) = 1.3_JPRM - 1.13_JPRM *  fwc_lsm(:,1)

END SELECT
        
       
WHERE ( ftl_lsm(:,1)-tfreeze < -5._JPRM ) fhrmodel(:,1) = 0._JPRM
WHERE ( ftl_lsm(:,1)-tfreeze < -5._JPRM ) fhrmodel(:,2) = 0._JPRM
WHERE ( ftfrac(:,7) >  0.99999_JPRM )  fhrmodel(:,1) = 0._JPRM 
WHERE ( ftfrac(:,7) >  0.99999_JPRM )  fhrmodel(:,2) = 0._JPRM 
 

IF (LGPRINT) WRITE(NULOUT,*) 'Roughness parameter h', fhrmodel(1,:)

       
END SUBROUTINE CMEM_INIT

