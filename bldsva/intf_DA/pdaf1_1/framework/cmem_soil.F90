! Copyright 2006-2014 ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities granted to it by
! virtue of its status as an intergovernmental organisation nor does it submit to any jurisdiction.
!
SUBROUTINE CMEM_SOIL

! Purpose :
! -------
!   CMEM SOIL MODULE
!   Calculate surface emissivity 
   
! Interface :
! ---------

! Externals :
! ---------
!   dielco_sub.F90 : contains all dielectric models
!   smref_sub.F90 : contains models for smooth surface reflectivity
!   rghref_sub.F90 : contains models for rough surface reflectivity

! Method :
!-------

! Reference :
! ---------

! Author :
! ------
! 23-Aug-2006 Thomas Holmes   *ECMWF*
! Revised:      
!     Patricia de Rosnay  August 2007
!                         Nov 2008 (Tskin for water bodies: lakes and sea)
!  February 2012 P. de Rosnay  (Tessel option, bug fix in TEFF and VEGTABLE)
!  May 2012      P. de Rosnay ECMWF: Wilheit model flexible interpolation
!                                    of soil moisture and temperature profiles
!                                    from the LSM vertical layers to the MW model layers 
!                                    Bug fixed on 29 May 2012 (line 277 hrmodel required to test roughness)

! Modifications :
! -------------
! End Modifications

!   alpha : fitting parameter (Wang and Schmugge)
! z_mw : depth of center of layer [m]
! dz : layer thickness [m]
!   medium : choose pure water (0), sea water (1), soil water (2)
!   isal : choose model for diel_wat, stogryn (1), Klein&Swift (2)
!------------------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRM
USE YOMLUN   , ONLY : NULOUT

USE YOMCMEMPAR ! ip_rgh_surf
USE YOMCMEMFIELDS
USE YOMCMEMSOIL ! alpha, ef

IMPLICIT NONE

INTEGER(KIND=JPIM) :: i, jmw,jls
INTEGER(KIND=JPIM) :: medium, isal
REAL(KIND=JPRM) :: frostfrac, sal
REAL(KIND=JPRM) :: z_mw, dz_mw, z_mw_max, prev_zmw, lev_zmw,prev_zls,lev_zls
REAL(KIND=JPRM)  :: interp_coef(nlay_soil_mw, nlay_soil_ls)

REAL,ALLOCATABLE   ::tlv(:)
REAL,ALLOCATABLE   ::zlv(:)
REAL,ALLOCATABLE   ::wclv(:)

COMPLEX(KIND=JPRM) :: ew, j

!------------------------------------------------------------------------------
! 0.0 Initialize
!------------------------------------------------------------------------------

surf_emis(:) = (/0.,0./)
tb_soil(:) = (/0.,0./)
    
! 1. Input from datafields
!    ---------------------
 wc1 = fwc_lsm(JJ,1)
 sand = fsand(JJ)
 clay = fclay(JJ)
 rho_b = frho_b(JJ)
 p = fp(JJ)
 sal_wat = fsal_wat(JJ)

ALLOCATE (eps_soil(nlay_soil_mw))
! ALLOCATE (tsoil(nlay_soil_mw))
ALLOCATE (zsoil(nlay_soil_mw))
ALLOCATE (wcsoil(nlay_soil_mw))

ALLOCATE (tlv(nlay_soil_ls))
ALLOCATE (zlv(nlay_soil_ls))
ALLOCATE (wclv(nlay_soil_ls))

! Interpolate the land surface model input on the MW model vertical resolution
! LSM variables: ftl_lsm,fwc_lsm,z_lsm,nlay_soil_ls
! MW  variables: tsoil,wcsoil,z,nlay_soil_mw
zsoil(1) = z_lsm(1)
tsoil(1) = ftl_lsm(JJ,1) 
wcsoil(1) = fwc_lsm(JJ,1)

IF (nlay_soil_mw > 1) THEN
! OLD interpol adapted to arrays for lsm and mw, but still not flexible (and still wrong...)
!  z_mw_max = z_lsm(nlay_soil_ls) 
!  dz_mw = z_mw_max/(nlay_soil_mw)
!  z_mw = -dz_mw/2.
!  DO JLAYER = 1,nlay_soil_mw
!   z_mw = z_mw + dz_mw   
!   zsoil(JLAYER) = dz_mw*100.
!
!   IF (z_mw <= z_lsm(1)) THEN
!     tsoil(JLAYER) = ftl_lsm(JJ,1)
!     wcsoil(JLAYER) = fwc_lsm(JJ,1)
!   ELSEIF (z_mw <= z_lsm(2)) THEN
!     tsoil(JLAYER) = ftl_lsm(JJ,1) + &
!                &   (ftl_lsm(JJ,2)-ftl_lsm(JJ,1)) * ((z_mw-z_lsm(1))/(z_lsm(2)-z_lsm(1)))
!     wcsoil(JLAYER) = fwc_lsm(JJ,1) + (fwc_lsm(JJ,2)-fwc_lsm(JJ,1)) * ((z_mw-z_lsm(1))/(z_lsm(2)-z_lsm(1)))
!   ELSEIF (z_mw <= z_lsm(3)) THEN
!     tsoil(JLAYER) = ftl_lsm(JJ,2) + &
!                &   (ftl_lsm(JJ,3)-ftl_lsm(JJ,2)) * ((z_mw-z_lsm(2))/(z_lsm(3)-z_lsm(2)))
!     wcsoil(JLAYER) = fwc_lsm(JJ,2) + (fwc_lsm(JJ,3)-fwc_lsm(JJ,2)) * ((z_mw-z_lsm(2))/(z_lsm(3)-z_lsm(2)))
!   ELSE
!     tsoil(JLAYER) = ftl_lsm(JJ,3)
!     wcsoil(JLAYER) = fwc_lsm(JJ,3)
!   ENDIF
!  ENDDO


! Interpolation of soil moisture and temperature profiles from nlay_soil_ls to nlay_soil_mw
 
        z_mw_max = z_lsm(nlay_soil_ls) 
        dz_mw = z_mw_max/(nlay_soil_mw)

        interp_coef(:,:) = 0.0_JPRM
        prev_zmw = 0.0_JPRM             ! mw model: depth of above layer 
        DO jmw = 1, nlay_soil_mw        ! loop on mw model layers
          zsoil(jmw) = dz_mw*100.       ! convert into cm, to be used in wilheit subroutine
          lev_zmw = prev_zmw + dz_mw    ! depth of current layer for mw model vertical grid
          prev_zls = 0.0_JPRM           ! depth of above layer for the ls model vertical grid

          DO jls = 1, nlay_soil_ls      ! loop on lsm vertical layers
             IF ( jls == nlay_soil_ls .AND. z_lsm(jls) < lev_zmw ) THEN 
                ! should not happen because z_mw_max = z_lsm(nlay_soil_ls)
                lev_zls = lev_zmw
             ELSE
                lev_zls = z_lsm(jls)
             ENDIF
             interp_coef(jmw,jls) = MAX(MIN(lev_zmw,lev_zls)-MAX(prev_zmw, prev_zls), 0.0)/(lev_zmw-prev_zmw)
             prev_zls = lev_zls
          ENDDO
           prev_zmw = lev_zmw
        ENDDO
        !
        IF ( LGPRINT ) THEN
           WRITE(*,*) 'cmem_soil interpolation of LSM vertical to MW grid for Wilheit'
           DO jmw = 1,  nlay_soil_mw
              WRITE(NULOUT,*) jmw, '-', interp_coef(jmw,1:nlay_soil_ls)
           ENDDO
           WRITE(*,*) "SUM of interpolation coefficient"
           DO jmw = 1, nlay_soil_mw 
              WRITE(NULOUT,*) jmw, '-', SUM(interp_coef(jmw,1:nlay_soil_ls))
           ENDDO
        ENDIF
        !
    tsoil(:) = 0._JPRM
    wcsoil(:) = 0._JPRM 

    DO jls = 1, nlay_soil_ls
      DO jmw = 1, nlay_soil_mw
          wcsoil(jmw) = wcsoil(jmw) + fwc_lsm(JJ,jls)*interp_coef(jmw,jls)
          tsoil(jmw)  = tsoil(jmw) + ftl_lsm(JJ,jls)*interp_coef(jmw,jls)
          zsoil(jmw)  = z_lsm(jls)
      ENDDO          
    ENDDO

ENDIF

DO jls = 1, nlay_soil_ls
      wclv(jls) = fwc_lsm(JJ,jls)
      tlv(jls)  = ftl_lsm(JJ,jls)
      zlv(jls)  = z_lsm(jls)
ENDDO


! 2. Compute Smooth Surface Reflectivity of surface
!    -----------------------------------------------
TILE: SELECT CASE (JCMEMTILE)  

! 2A. Water tile
!------------------------------------------------------------------------------
  CASE ( 7 ) TILE
    fh(JJ) = 0.
    
    ! 2A.1 Effective temperature : tsea = tsoil1 (sea+ice)  
    ! t_eff(:) = (C,tsoil(1)/)
    ! t_eff(:) = (/1.0,ftskin(JJ),ftskin(JJ)/)
    t_eff(1) = tsoil(1)
    t_eff(2) = tsoil(1)
    t_eff(3) = tsoil(1)
    
    ! 2A.2 Compute Dielectric Constant of surface medium
    SELECT CASE ( ((t_eff(2)-tfreeze) < -0.5) )
      CASE ( .TRUE. )  ! water is pure ice
        CALL DIEL_ICE (t_eff(2),ew)
      CASE DEFAULT
        medium = 1
        isal = 2
        CALL DIEL_WAT (medium,isal,t_eff(2)-tfreeze,sal_wat,ew)
    END SELECT
      
    ! 2A.3. Compute Smooth Surface Reflectivity
    CALL FRESNEL (ew)           
IF (LGPRINT) WRITE(NULOUT,*) '--- cmem soil 2A.3:', ew,t_eff(:)
    
! 2B. Soil tiles
!------------------------------------------------------------------------------
  CASE ( :6 ) TILE
  
    ! 2B.1 Effective temperature of soil medium
    SELECT CASE (CITEFF)
      CASE ( 'Tsoil' )
        t_eff(2) = tsoil(1)
        t_eff(3) = tsoil(1)
      CASE ( 'Choudhury', 'Wigneron','Holmes','Lv' )  
          ! IF (LGPRINT) WRITE(NULOUT,*) ,wclv,tlv,zlv
        CALL TEFF_SUB(tsoil(1),tsoildeep,wclv,tlv,zlv)
    END SELECT

  ! 2B.2 Compute Dielectric Constant of soil medium
  DO JLAYER = 1, nlay_soil_mw   
    ! Temperature of soil JLAYER
    SELECT CASE (CITDIEL)
      CASE ( 'Teff' )  
        tc = t_eff(2) - tfreeze
      CASE ( 'Tsurf' )  
        tc = tsoil(JLAYER) - tfreeze
    END SELECT    
    ! Water content of soil JLAYER
    ! wc = wcsoil(JLAYER)
    wc = wclv(1)
    ! dielectric constant of surface medium  
    DRY: SELECT CASE ((wc < 0.02) .AND. (sand > 90.) .AND. fghz<10.)
    
      CASE ( .TRUE. ) DRY ! Microwave 1-10GHz permittivity of dry sand (matzler '98)    
          j = (0.,1.)     !   relaxation frequency of 0.27 GHz 
          eps = 2.53 + (2.79 - 2.53) / (1. - j * (fghz/0.27)) + j * 0.002
  
      CASE DEFAULT DRY
        ! Calculate dielectric constant of soil water
        ICE: SELECT CASE ( (tc < -0.5) )
          CASE ( .TRUE. ) ICE
            CALL DIEL_ICE (tc+tfreeze,ew)
          CASE DEFAULT ICE
            isal = 2
            medium = 2
            CALL DIEL_WAT (medium,isal,tc,sal_soil,ew)
        END SELECT ICE
        ! dielectric of surface medium    
        SELECT CASE (CIDIEL)
          CASE ( 'Wang' )  
            wp = fWP(JJ)
            alpha = falpha(JJ)
            CALL DIELWANG (ew)
          CASE ( 'Dobson' )            
            CALL DIELDOBSON (ew)
          CASE ( 'Mironov')
            CALL DIELMIRONOV 
        END SELECT    

     END SELECT DRY
     
     ! mix dielectric constant of frozen and non-frozen soil
     IF (tc < -5.) THEN ! frozen soil
       frostfrac = 1.
     ELSEIF (tc < -0.5) THEN
       frostfrac = 1.0-2.6945*ABS(tc-0.05)**(-0.6104)/100
     ELSE 
       frostfrac = 0.
     ENDIF
     eps = eps * (1.-frostfrac) + ef * frostfrac 
     eps_soil(JLAYER) = eps
     
  ENDDO

  ! 2B.3. Compute Smooth Surface Reflectivity
  !       -----------------------------------
  SELECT CASE (CISMR)
    CASE ( 'Fresnel' )  
      CALL FRESNEL (eps_soil(1))            
    CASE ( 'Wilheit' )  
      CALL WILHEIT
  END SELECT
   

IF (LGPRINT) WRITE(NULOUT,*) '--- REFLSOL :',r_s(:)
IF (LGPRINT) WRITE(NULOUT,*) '--- TEFF cmem soil 2B.3 :',t_eff(:)
END SELECT TILE 
!------------------------------------------------------------------------------

! 4. Compute Rough Surface Reflectivity
!   ----------------------------------
IF ( hrmodel == 0. .AND. ip_Q == 0. ) THEN
  r_r(:) = r_s(:)
ELSE
  SELECT CASE (CIRGHR)
    CASE ( 'No' )  
      r_r(:) = r_s(:)
    CASE ( 'Choudhury','Wsimple','Wtexture','Wigneron' ) 
      CALL RGHCHOU
    CASE ( 'Wegmueller' )
      CALL RGHWEGM
  END SELECT
ENDIF
IF (LGPRINT) WRITE(NULOUT,*) '--- REFL ROUGH:', r_r

! 5. Compute Surface Emissivity 
!   ------------------------------
DO i = 1,2    ! h- and v-polarization
  surf_emis(i) = 1.0 - r_r(i)
ENDDO

! 6. Compute surface brightness temperatures
!     --------------------------------------
DO i = 1,2    ! h- and v-polarization; t_eff(i+1) is teff for  H (2) & V(3)
  tb_soil(i) = t_eff(i+1) * surf_emis(i)
ENDDO
IF (LGPRINT) WRITE(NULOUT,*) '--- Tb soil :', tb_soil

! 7. Clean
!    -----
DEALLOCATE (eps_soil)
DEALLOCATE (zsoil)
!DEALLOCATE (tsoil)
DEALLOCATE (wcsoil) 
DEALLOCATE (tlv)
DEALLOCATE (zlv)
DEALLOCATE (wclv)
 
IF (LGPRINT) WRITE(NULOUT,*) '--- teff cmem soil end:', t_eff(:) 
END SUBROUTINE CMEM_SOIL
