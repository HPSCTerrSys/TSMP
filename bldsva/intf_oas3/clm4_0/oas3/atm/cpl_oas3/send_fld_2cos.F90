SUBROUTINE send_fld_2cos(nstep, dtime, x2a)

!---------------------------------------------------------------------
! Description:
!  This routine sends the fluxes from CESM to COSMO, via x2a vectors
!  CESM sign convention is that fluxes are positive downward

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
USE mct_mod
!

!==============================================================================

IMPLICIT NONE

!==============================================================================

INTEGER,  INTENT(IN)             :: nstep, dtime
TYPE(mct_aVect) ,INTENT(IN)      :: x2a            !MCT vector x to a, activates with prognostic

! Local Variables:
INTEGER                          :: isec
INTEGER                          :: k, k1, k2      ! mct_vectors ID
! processor bounds indices
INTEGER, PARAMETER ::   jps_taux   =  1     !  zonal wind stress
INTEGER, PARAMETER ::   jps_tauy   =  2     !  meridional wind stress
INTEGER, PARAMETER ::   jps_lat    =  3     !  total latent heat flux (W/m**2)
INTEGER, PARAMETER ::   jps_sens   =  4     !  total sensible heat flux (W/m**2)
INTEGER, PARAMETER ::   jps_ir     =  5     !  emitted infrared (longwave) radiation (W/m**2)
INTEGER, PARAMETER ::   jps_albd   =  6     !  direct albedo
INTEGER, PARAMETER ::   jps_albi   =  7     !  diffuse albedo
INTEGER, PARAMETER ::   jps_co2fl  =  8     !  net CO2 flux (now only photosynthesis rate) (umol CO2 m-2s-1)
INTEGER, PARAMETER ::   jps_ram1   =  9     !  Aerodynamic Resistance (s/m). !CPS
INTEGER, PARAMETER ::   jps_rah1   =  10    !  Aerodynamic Resistance (s/m). !CPS
INTEGER, PARAMETER ::   jps_raw1   =  11    !  Aerodynamic Resistance (s/m). !CPS
INTEGER, PARAMETER ::   jps_tsf1   =  12    !  Surface Temperature (K)  !CPS
INTEGER, PARAMETER ::   jps_qsf1   =  13    !  Surface Humidity (kg/kg) !CPS
INTEGER, PARAMETER ::   jps_fpsn   =  14    !  photosynthesis rate (umol CO2 m-2s-1)
INTEGER, PARAMETER ::   jps_fplres =  15    !  plant respiration (umol CO2 m-2s-1)
!
INTEGER                          :: rank,info,n
REAL(KIND=r8), ALLOCATABLE       :: fsnd(:,:)      ! temporary arrays
!------------------------------------------------------------------------------
!- End of header
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!- Begin Subroutine send_fld2cos 
!------------------------------------------------------------------------------

   CALL MPI_Comm_Rank(kl_comm, rank, nerror) 
   IF (nerror /= 0) CALL prism_abort_proto(ncomp_id, 'MPI_Comm_Rank', 'Failure in send_fld_2cos') 

   ALLOCATE ( fsnd(start1d:start1d+length1d-1, ksnd), stat=nerror)  
   IF (nerror /= 0) THEN 
     CALL prism_abort_proto( ncomp_id, 'send_fld_2cos', 'Failure in allocating fsnd' )
     RETURN
   ENDIF

 ! zero on unmasked points
 fsnd = -999999._r8
  
 isec = nstep*dtime

 ! Retrieve x2a Fields to send to COSMO
 ! flux prefix, Faxx , between a and x, computed by x
 ! Sign change done in COSMO, except for IR, to be consistent with CLM3.5
 k                = mct_aVect_indexRA(x2a,'Faxx_taux')
 fsnd(:,jps_taux) = x2a%rAttr(k,:)
 k                = mct_aVect_indexRA(x2a,'Faxx_tauy') 
 fsnd(:,jps_tauy) = x2a%rAttr(k,:)
 k                = mct_aVect_indexRA(x2a,'Faxx_lat')
 fsnd(:,jps_lat)  = x2a%rAttr(k,:)
 k                = mct_aVect_indexRA(x2a,'Faxx_sen')
 fsnd(:,jps_sens) = x2a%rAttr(k,:)
 k                = mct_aVect_indexRA(x2a,'Faxx_lwup')
 fsnd(:,jps_ir)   = -1._r8 * x2a%rAttr(k,:)
 ! state prefix, Sx_, from coupler 
 k1               = mct_aVect_indexRA(x2a,'Sx_avsdr')
 k2               = mct_aVect_indexRA(x2a,'Sx_anidr')
 fsnd(:,jps_albd) = 0.5_r8 * x2a%rAttr(k1,:) + 0.5_r8 * x2a%rAttr(k2,:)
 k1               = mct_aVect_indexRA(x2a,'Sx_avsdf')
 k2               = mct_aVect_indexRA(x2a,'Sx_anidf')
 fsnd(:,jps_albi) = 0.5_r8 * x2a%rAttr(k1,:) + 0.5_r8 * x2a%rAttr(k2,:)
 ! 
 ! Missing Co2 from land via coupler
 fsnd(:,jps_co2fl)= 300._r8  !CPS

 ! Zonal surface stress (N m-2) 
 IF( ssnd(jps_taux)%laction )  &
   CALL oas_clm_snd( jps_taux, isec, fsnd(:,jps_taux),start1d,start1d+length1d-1,info )
 ! Meridional surface stress (N m-2)
 IF( ssnd(jps_tauy)%laction )  &
   CALL oas_clm_snd( jps_tauy, isec, fsnd(:,jps_tauy),start1d,start1d+length1d-1,info )
 ! Surface latent heat flux (W m-2)
 IF( ssnd(jps_lat)%laction )  &
   CALL oas_clm_snd( jps_lat, isec, fsnd(:,jps_lat),start1d,start1d+length1d-1,info )
 ! Surface sensible heat flux (W m-2)
 IF( ssnd(jps_sens)%laction )  &
   CALL oas_clm_snd( jps_sens, isec, fsnd(:,jps_sens),start1d,start1d+length1d-1,info )
 ! Surface upward longwave heat flux (W m-2)
 IF( ssnd(jps_ir)%laction )  &
   CALL oas_clm_snd( jps_ir, isec, fsnd(:,jps_ir),start1d,start1d+length1d-1,info )
 ! Direct albedo (visible radiation and near-infrared radiation)
 IF( ssnd(jps_albd)%laction )  &
   CALL oas_clm_snd( jps_albd, isec, fsnd(:,jps_albd),start1d,start1d+length1d-1,info )
 ! Direct albedo (visible radiation and near-infrared radiation)
 IF( ssnd(jps_albi)%laction )  &
   CALL oas_clm_snd( jps_albi, isec, fsnd(:,jps_albi),start1d,start1d+length1d-1,info )
 ! surface_upward_flux_of_carbon_dioxide_where_land (moles m-2 s-1) 
 IF( ssnd(jps_co2fl)%laction )  &
   CALL oas_clm_snd( jps_co2fl, isec, fsnd(:,jps_co2fl),start1d,start1d+length1d-1,info )
 
 !CALL MPI_Barrier(kl_comm, nerror)

 DEALLOCATE(fsnd)
!------------------------------------------------------------------------------
!- End of the Subroutine
!------------------------------------------------------------------------------


END SUBROUTINE send_fld_2cos
