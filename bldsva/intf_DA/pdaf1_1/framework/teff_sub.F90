SUBROUTINE TEFF_SUB (TSURF, TDEEP,wclv,tlv,zlv)

!---------------------------------------------------------------------------
! Models to parameterize the effective temperature for microwave observations, 
! using a surface and a deep soil temperature.
! CITEFF  =    Tsoil: teff = tsoil1

!              Choudhury:  Choudhury et al., 1982
!   Reference :
!   Choudhury, B., T. Schmugge, and T. Mo (1982), A parameterization of effective 
!    soil temperature for microwave emission, J. Geophys. Res., 87, 1301–1304

!              2) Wigneron et al., 2001
!   Reference :
!   Wigneron, J.-P., L. Laguerre, and Y. Kerr (2001), A simple parmeterization 
!    of the L-band microwave emission from rough agricultural soils, IEEE Trans. 
!    Geosci. Remote. Sens., 39, 1697–1707.

!              3) Holmes et al., 2006
! Reference :
!   Holmes, T. R. H., P. de Rosnay, R. de Jeu, R. J.-P. Wigneron, Y. Kerr, J.-C. 
!    Calvet, M. J. Escorihuela, K. Saleh, and F. Lemaître (2006), A new 
!    parameterization of the effective !temperature for L band radiometry, 
!    Geophys. Res. Lett., 33, L07405, doi:10.1029/2006GL025724.

! Input Variables :
!   Tsurf: Surface Temperature, e.g. at 1-5cm (K)
!   Tdeep : Deep Temperature, e.g. at 50-100cm (K)
!   C : fitting parameter (m/m)

!   P. de Rosnay CMWEM v4.0 Fix in TEFF computation and archiving when Holmes is used - Feb 2012
!---------------------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM,JPRM
USE YOMLUN   , ONLY : NULOUT
USE YOMCMEMPAR, ONLY : CITEFF, lamcm, tfreeze, LGPRINT, nlay_soil_ls, f,pi
USE YOMCMEMSOIL, ONLY : T_EFF, sal_soil,tc, wc1, wc, w0, bw, eps, eps0, bh, alpha, wp, ef
USE YOMCMEMFIELDS, ONLY : JJ, fWP, falpha

IMPLICIT NONE

INTEGER(KIND=JPIM) :: medium,isal,jls
REAL(KIND = JPRM) :: tsurf
REAL(KIND = JPRM) :: tdeep
REAL(KIND = JPRM) :: Depth
REAL(KIND = JPRM),intent(in)   :: wclv(nlay_soil_ls)
REAL(KIND = JPRM),intent(in)   :: tlv(nlay_soil_ls)
REAL(KIND = JPRM),intent(in)   :: zlv(nlay_soil_ls)
REAL(KIND = JPRM) :: Bsoil
REAL(KIND = JPRM) :: C 
REAL(KIND = JPRM) :: Bsoil_sum 

REAL(KIND = JPRM) :: frostfrac
COMPLEX(KIND=JPRM) :: ew
!---------------------------------------------------------------------------
  
SELECT CASE (CITEFF)

  CASE ( 'Choudhury' )
    IF ( lamcm < 4.4 )  THEN
      C = 0.802
    ELSEIF ( lamcm < 8.5 ) THEN
      C = 0.667
    ELSEIF ( lamcm < 16. ) THEN
      C = 0.480
    ELSEIF ( lamcm < 35. ) THEN
      C = 0.246
    ELSE 
      C = 0.084
    ENDIF
  
  CASE ( 'Wigneron' )
    C = MAX(0.001,(wc1/w0)**bw)
  
  CASE ( 'Holmes' )
     tc = tsurf - tfreeze
     wc = wc1
     wp = fWP(JJ)
     alpha = falpha(JJ)
        ! Calculate dielectric constant of soil water
        ICE: SELECT CASE ( (tc < -0.5) )
          CASE ( .TRUE. ) ICE
            CALL DIEL_ICE (tc+tfreeze,ew)
          CASE DEFAULT ICE
            isal = 2
            medium = 0 ! teff calibrated with permittivity of saline water
            CALL DIEL_WAT (medium,isal,tc,sal_soil,ew)
        END SELECT ICE
     CALL dielwang (ew)
     ! parameterization does not take conductivity loss into account 
     ! remove loss factor from imaginary part
     eps=eps - (0.,1.) * alpha * wc**2
       ! frozen soil influence
     IF (tc < -5.) THEN ! frozen soil
       frostfrac = 1.
     ELSEIF (tc < -0.5) THEN
       frostfrac = 0.5
     ELSE 
       frostfrac = 0.
     ENDIF
     eps = eps * (1.-frostfrac) + ef * frostfrac 
     C = MAX(0.001_JPRM,MIN(1.0_JPRM,( ( aimag(eps) / real(eps) ) /eps0 )**bh))
     ! re-initialize eps and ew
     eps =(0.,0.)
     ew  =(0.,0.)

  CASE ( 'Lv' )     

     ! 6 layers as standard configuration: 5/10/20/40/80/160 cm
     ! (1) for the layer 0-7.5 cm
     ! (2) for the layer 7.5-15 cm
     ! (3) for the layer 15-30 cm
     ! (4) for the layer 30-60 cm
     ! (5) for the layer 60-120 cm

     DO jls = 1, nlay_soil_ls-1
       tc = tlv(jls) - tfreeze
       wc = wclv(jls)
       IF (jls == 1) THEN
            Depth = zlv(1)
       ELSE
            Depth = zlv(jls)
     !       Depth = zlv(jls+1)-zlv(jls)
       ENDIF

       wp = fWP(JJ)
       alpha = falpha(JJ)
          ! Calculate dielectric constant of soil water
          SELECT CASE ( (tc < -0.5) )
            CASE ( .TRUE. )
              CALL DIEL_ICE (tc+tfreeze,ew)
            CASE DEFAULT
              isal = 2
              medium = 0 ! teff calibrated with permittivity of saline water
              CALL DIEL_WAT (medium,isal,tc,sal_soil,ew)
          END SELECT
       CALL dielwang (ew)
       ! parameterization does not take conductivity loss into account 
       ! remove loss factor from imaginary part
       eps=eps - (0.,1.) * alpha * wc**2
           ! frozen soil influence
       IF (tc < -5.) THEN ! frozen soil
         frostfrac = 1.
       ELSEIF (tc < -0.5) THEN
         frostfrac = 1.0-2.6945*ABS(tc-0.05)**(-0.6104)/100
       ELSE 
         frostfrac = 0.
       ENDIF
       eps = eps * (1.-frostfrac) + ef * frostfrac 
       Bsoil = Depth* pi* &
                    2.0 / (3.0e8/f) * aimag(eps) / real(eps) **0.5
       ! re-initialize eps and ew
       eps =(0.,0.)
       ew =(0.,0.)

       ! loop on lsm vertical layers  
       IF (jls == 1) THEN
            tdeep = tsurf*(1-exp(-Bsoil))
            Bsoil_sum = Bsoil
       ELSE
            tdeep = tdeep+tlv(jls)*(1-exp(-Bsoil))*exp(-Bsoil_sum)  
            Bsoil_sum = Bsoil_sum + Bsoil
       ENDIF

     ENDDO
              
     C        = 0

ENDSELECT

t_eff(1) = C
t_eff(2) = tdeep + (tsurf-tdeep)*C
t_eff(3) = t_eff(2)

IF (LGPRINT) WRITE(NULOUT,*) '--- TEFF (teff_sub.F90):',C, tdeep,tsurf,t_eff(:)


END SUBROUTINE TEFF_SUB
