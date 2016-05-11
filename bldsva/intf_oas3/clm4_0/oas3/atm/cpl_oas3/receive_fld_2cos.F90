SUBROUTINE receive_fld_2cos(nstep, dtime, a2x, lcoupled)

!---------------------------------------------------------------------
! Description:
!  This routine receives the atmospheric fields from COSMO and updates
!  the DATM a2x arrays 
!
! Current Code Owner: TR32, Z4: Prabhakar Shrestha
!    phone: 0228733453
!    email: pshrestha@uni-bonn.de
!
! History:
! Version    Date       Name
! ---------- ---------- ----
! 2.1.0        2016/02/29 Prabhakar Shrestha
! Implementation for CESM 1.2.1
! @VERSION@    @DATE@     <Your name>
!  <Modification comments>         
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!==============================================================================

! Declarations:
!
! Modules used:

USE oas_clm_vardef
USE shr_kind_mod ,            ONLY : R8 => SHR_KIND_R8
USE mct_mod
use shr_const_mod,            ONLY : SHR_CONST_TKFRZ,  SHR_CONST_RDAIR
!==============================================================================

IMPLICIT NONE

!==============================================================================

! Parameter List
INTEGER, INTENT(IN)                    ::   nstep           ! time step
INTEGER, INTENT(IN)                    ::   dtime           ! dt
LOGICAL, INTENT(OUT)                   ::   lcoupled        ! CPS flag for coupling 
TYPE(mct_aVect) ,INTENT(INOUT)         ::   a2x             !MCT Vector a to x
!
INTEGER, PARAMETER ::   jps_t   =  1            ! temperature
INTEGER, PARAMETER ::   jps_u   =  2            ! u wind
INTEGER, PARAMETER ::   jps_v   =  3            ! v wind
INTEGER, PARAMETER ::   jps_q   =  4            ! specific water vapor content
INTEGER, PARAMETER ::   jps_th  =  5            ! thickness of lowest level (m)
INTEGER, PARAMETER ::   jps_pr  =  6            ! surface pressure (Pa)
INTEGER, PARAMETER ::   jps_rs  =  7            ! direct shortwave downward radiation (W/m2)
INTEGER, PARAMETER ::   jps_fs  =  8            ! diffuse shortwave downward radiation (W/m2)
INTEGER, PARAMETER ::   jps_lw  =  9            ! longwave downward radiation (W/m2) 
INTEGER, PARAMETER ::   jps_cr  = 10            ! convective rain precipitation      (kg/m2*s)
INTEGER, PARAMETER ::   jps_cs  = 11            ! convective snow precipitation      (kg/m2*s)
INTEGER, PARAMETER ::   jps_gr  = 12            ! gridscale rain precipitation
INTEGER, PARAMETER ::   jps_gs  = 13            ! gridscale snow precipitation
INTEGER, PARAMETER ::   jps_gg  = 14            ! gridscale graupel precipitation
INTEGER, PARAMETER ::   jps_cp  = 15            ! total convective precipitation
INTEGER, PARAMETER ::   jps_gp  = 16            ! total gridscale precipitation
INTEGER, PARAMETER ::   jps_co2 = 17            ! CO2 partial pressure (Pa)

INTEGER                                :: n, k, k1, k2   !INDICES
REAL(KIND=r8),PARAMETER                :: tkFrz = SHR_CONST_TKFRZ ! freezing T of fresh water ~ K 
REAL(KIND=r8),PARAMETER                :: rdair  = SHR_CONST_RDAIR! dry air gas constant   ~ J/K/kg
REAL(KIND=r8)                          :: vp, frac, qbot

INTEGER                                :: jn, isec, ier , begg,endg
INTEGER, DIMENSION(krcv)               :: nrcvinfo           ! OASIS info argument
INTEGER                                :: rank
REAL(KIND=r8), ALLOCATABLE             :: ztmp1(:)
!------------------------------------------------------------------------------
!- End of header
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!- Begin Subroutine receive_fld2cos 
!------------------------------------------------------------------------------

 CALL MPI_Comm_Rank(kl_comm, rank, nerror) 
 IF (nerror /= 0) CALL prism_abort_proto(ncomp_id, 'MPI_Comm_Rank', 'Failure in receive_fld_2cos') 

 isec = dtime * nstep 

 ! Default value (changed when coupling)
 lcoupled = .FALSE.
 nrcvinfo = OASIS_idle
 begg     = start1d
 endg     = start1d+length1d-1

 ALLOCATE ( ztmp1(begg:endg), stat=nerror)
 IF (nerror /= 0) THEN
   CALL prism_abort_proto( ncomp_id, 'receive_fld_2cos', 'Failure in allocating ztmp1' )
   RETURN
 ENDIF
!
! Update DATM variable with the received input from COSMO
 DO jn = 1, krcv
   ztmp1 = -9999._r8
   nrcvinfo = OASIS_idle
   IF(srcv(jn)%laction) CALL oas_clm_rcv( jn, isec, ztmp1,begg,endg, nrcvinfo(jn) )
   IF( nrcvinfo(jn)==OASIS_Rcv) THEN   
      exfld(:,jn)=ztmp1(:)
   ENDIF
 ENDDO

 ! temperature at the lowest model level (K)
 IF( srcv(jps_t)%laction ) THEN
   k              = mct_aVect_indexRA(a2x,'Sa_tbot')
   a2x%rAttr(k,:) = exfld(:,jps_t)
   k              = mct_aVect_indexRA(a2x,'Sa_ptem')
   a2x%rAttr(k,:) = exfld(:,jps_t)
 ENDIF
 ! zonal wind at the lowest model level (m s-1)
 IF( srcv(jps_u)%laction ) THEN
   k              = mct_aVect_indexRA(a2x,'Sa_u')
   a2x%rAttr(k,:) = exfld(:,jps_u)
 ENDIF
 ! meridional wind at the lowest model level (m s-1)
 IF( srcv(jps_v)%laction ) THEN
   k              = mct_aVect_indexRA(a2x,'Sa_v')
   a2x%rAttr(k,:) = exfld(:,jps_v)
 ENDIF
 ! Specific humidity at the lowest model level (kg kg-1)
 IF( srcv(jps_q)%laction ) THEN
   k              = mct_aVect_indexRA(a2x,'Sa_shum')
   DO n = 1, length1d
      qbot = exfld(n,jps_q)
      a2x%rAttr(k,n) = qbot/(qbot + 1._r8)   !CPS Mixing Ratio to Sp. Humidity 
   ENDDO
 ENDIF
 ! height at the lowest model level (m)
 IF( srcv(jps_th)%laction ) THEN
   k              = mct_aVect_indexRA(a2x,'Sa_z')
   a2x%rAttr(k,:) = exfld(:,jps_th)
 ENDIF
 ! pressure at the lowest model level (Pa)
 IF( srcv(jps_pr)%laction ) THEN
   k              = mct_aVect_indexRA(a2x,'Sa_pbot')
   a2x%rAttr(k,:) = exfld(:,jps_pr)
   k              = mct_aVect_indexRA(a2x,'Sa_pslv')
   a2x%rAttr(k,:) = exfld(:,jps_pr)
 ENDIF
 ! direct near-infrared/visible incident solar radiation (W m-2) 
 IF( srcv(jps_rs)%laction ) THEN
   k              = mct_aVect_indexRA(a2x,'Faxa_swndr')
   a2x%rAttr(k,:) = exfld(:,jps_rs) * 0.50_r8
   k              = mct_aVect_indexRA(a2x,'Faxa_swvdr')
   a2x%rAttr(k,:) = exfld(:,jps_rs) * 0.50_r8
 ENDIF
 ! diffuse near-infrared/visible incident solar radiation (W m-2)
 IF( srcv(jps_fs)%laction ) THEN
   k              = mct_aVect_indexRA(a2x,'Faxa_swndf')
   a2x%rAttr(k,:) = exfld(:,jps_fs) * 0.50_r8
   k              = mct_aVect_indexRA(a2x,'Faxa_swvdf')
   a2x%rAttr(k,:) = exfld(:,jps_fs) * 0.50_r8
 ENDIF
 ! downward longwave heat flux (W m-2)
 IF( srcv(jps_lw)%laction ) THEN
   k              = mct_aVect_indexRA(a2x,'Faxa_lwdn')
   a2x%rAttr(k,:) = exfld(:,jps_lw)
 ENDIF 
 ! convective precipitaiton rate, rain/snow  (kg m-2 s-1)
 IF( srcv(jps_cp)%laction ) THEN
   k1             = mct_aVect_indexRA(a2x,'Faxa_rainc')
   k2             = mct_aVect_indexRA(a2x,'Faxa_snowc')
   DO n = 1, length1d
     frac            = (exfld(n,jps_t) - tkFrz) * 0.5_r8   !ramp near freezing,, NCEPCLM
     frac            = MIN(1.0_r8, MAX(0.0_r8,frac))
     a2x%rAttr(k1,n) = MAX(0.0_r8, exfld(n,jps_cp) * (         frac) ) 
     a2x%rAttr(k2,n) = MAX(0.0_r8, exfld(n,jps_cp) * (1.0_r8 - frac) ) 
   ENDDO
 ENDIF
 ! large-scale (stable) precipitation rate (kg m-2 s-1)
 IF( srcv(jps_gp)%laction ) THEN
   k1             = mct_aVect_indexRA(a2x,'Faxa_rainl')
   k2             = mct_aVect_indexRA(a2x,'Faxa_snowl')
   DO n = 1, length1d
     frac            = (exfld(n,jps_t) - tkFrz) * 0.5_r8   !ramp near freezing,, NCEPCLM
     frac            = MIN(1.0_r8, MAX(0.0_r8,frac))
     a2x%rAttr(k1,n) = MAX(0.0_r8, exfld(n,jps_gp) * (         frac) )
     a2x%rAttr(k2,n) = MAX(0.0_r8, exfld(n,jps_gp) * (1.0_r8 - frac) )        
   ENDDO 
 ENDIF
 ! prognostic CO2 at the lowest model level (1e-6 mol/mol)
 IF( srcv(jps_co2)%laction ) THEN
   k              = mct_aVect_indexRA(a2x,'Sa_co2prog',perrWith='quiet')
   a2x%rAttr(k,:) = exfld(:,jps_co2)
 ENDIF

 k  = mct_aVect_indexRA(a2x,'Sa_dens')
 k1 = mct_aVect_indexRA(a2x,'Sa_shum')
 !--- density ---
 DO n = 1, length1d
   vp = (a2x%rAttr(k1,n)*exfld(n,jps_pr)) / (0.622_R8 + 0.378_R8 * a2x%rAttr(k1,n))
   a2x%rAttr(k,n) = (exfld(n,jps_pr) - 0.378_R8 * vp) / (exfld(n,jps_t)*rdair)
 ENDDO
! IF (rank == 0) PRINT*, "CPS DENSITY", MINVAL(a2x%rAttr(k,:)), MAXVAL(a2x%rAttr(k,:))

 DEALLOCATE(ztmp1)

 !  Wait coupling result (involved PE and non involved PE)
 !CALL MPI_Bcast( lcoupled, 1, dummy, 0, kl_comm, ier )

!------------------------------------------------------------------------------
!- End of the Subroutine
!------------------------------------------------------------------------------

END SUBROUTINE receive_fld_2cos
