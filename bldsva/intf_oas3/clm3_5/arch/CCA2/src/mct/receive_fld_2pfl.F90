SUBROUTINE receive_fld_2pfl(kt)

!---------------------------------------------------------------------
! Description:
!  This routine receives the subsurface fields (soil moisture / pressure head)
!  from ParFlow and converts into volumetric soil moisture / pressure head in 
!  CLM3.5 column level. A new routine g2p has been added to downscale data from
!  grid level to column level in the module subgridAveMod 
!
! Current Code Owner: TR32, Z4: Prabhakar Shrestha
!    phone: 0228733453
!    email: pshrestha@uni-bonn.de
!
! History:
! Version    Date       Name
! ---------- ---------- ----
! 1.1        2012/01/20 Prabhakar Shrestha 
!   Initial release for CLM3.5 oasis coupler
! @VERSION@    @DATE@     <Your name>
! 1.1        2012/01/30 Mauro Sulis
! Soil matrix potential received from Parflow points to pfl_psi variable
! defined in clmtype.F90 and allocated in clmtypeInitMod.F90
! 1.1        2012/02/27  Prabhakar Shrestha
! Bug fix for receive broadcast logical flags (pfl_coupled_sat & pfl_coupled_psi) 
! 1.1        2012/05/02  Mauro Sulis
! The volumetric moisture content received from Parflow points to the new variable pfl_h2osoi_vol
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
USE shr_kind_mod ,            ONLY : r8 => shr_kind_r8
USE clm_time_manager ,        ONLY : dtime,                                      &! timestep in second
                                     nelapse,                                    &! total simulation step
                                     get_nstep                                    ! return timestep number
USE clm_varpar,               ONLY : nlevsoi
USE clm_varcon,               ONLY : denh2o                                       ! density of liquid water [kg/m3]
USE clmtype,                  ONLY : clm3                                         ! input fields
USE spmdMod ,                 ONLY : masterproc,iam,npes,mpicom, dummy => MPI_LOGICAL,           &! CPS
                                     dummy_int => MPI_INTEGER                     ! CMS
USE decompMod ,               ONLY : get_proc_clumps, get_clump_bounds,          &
                                     get_proc_bounds,                            &
                                     get_proc_global, ldecomp
USE spmdGathScatMod,          ONLY : scatter_data_from_master
USE filterMod,                ONLY : filter, setFilters
USE netcdf

!==============================================================================

IMPLICIT NONE

!==============================================================================

! Parameter List
INTEGER, INTENT(IN)                    :: kt             ! time step
LOGICAL, DIMENSION(vrcv)               :: pfl_coupled_sat! CMS flag for coupling
LOGICAL, DIMENSION(vrcv)               :: pfl_coupled_psi! CMS flag for coupling 
!
! Local Variables
INTEGER                                :: g,c,fc,nc      ! indices
! processor bounds indices
INTEGER                                :: numg           ! total number of gridcells across all processors
INTEGER                                :: numl           ! total number of landunits across all processors
INTEGER                                :: numc           ! total number of columns across all processors
INTEGER                                :: nump           ! total number of pfts across all processors
INTEGER                                :: begg,endg      ! local beg/end gridcells gdc
INTEGER                                :: begl,endl      ! local beg/end landunits
INTEGER                                :: begc,endc      ! local beg/end columns
INTEGER                                :: begp,endp      ! local beg/end pfts
INTEGER                                :: nclumps        !
INTEGER                                :: num_soilc      ! number of column soil points in column filter
INTEGER, ALLOCATABLE                   :: filter_soilc(:)! column filter for soil points

! Pointers
REAL(KIND=r8), POINTER                 :: pfl_vol(:,:)   ! volumetric soil water (nlevsoi) 
REAL(KIND=r8), POINTER                 :: pfl_soi(:,:)   ! liquid water
REAL(KIND=r8), POINTER                 :: pfl_ice(:,:)   ! ice water   !CPS
REAL(KIND=r8), POINTER                 :: pfl_psi(:,:)   ! soil matrix potential received from ParFlow (mm)
REAL(KIND=r8), POINTER                 :: watsat(:,:)    ! volumetric soil water at saturation (porosity)
REAL(KIND=r8), POINTER                 :: dz(:,:)        ! layer depth (m)

REAL(KIND=r8), allocatable, DIMENSION(:)   :: rcv_field

INTEGER                                :: jn, jnid, isec, ier,it1,it2
INTEGER,  DIMENSION(vrcv)              :: k
INTEGER , DIMENSION(vrcv)              :: nrcvinfo           ! OASIS info argument

INTEGER ,DIMENSION(4)                  :: dimids                ! CPS 3 to 4
INTEGER ,DIMENSION(vrcv)               :: il_var_id
INTEGER                                :: il_file_id, status
REAL(KIND=r8), SAVE                    :: tic, toc

!------------------------------------------------------------------------------
!- End of header
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!- Begin Subroutine receive_fld2pfl 
!------------------------------------------------------------------------------

 isec = dtime * kt

 pfl_soi    => clm3%g%l%c%cws%h2osoi_liq 
 pfl_ice    => clm3%g%l%c%cws%h2osoi_ice  !CPS
 pfl_vol    => clm3%g%l%c%cws%pfl_h2osoi_vol
 pfl_psi    => clm3%g%l%c%cps%pfl_psi
 watsat     => clm3%g%l%c%cps%watsat
 dz         => clm3%g%l%c%cps%dz 

! Get global information 
 CALL get_proc_global(numg,numl,numc,nump)
 CALL get_proc_bounds(begg,endg,begl,endl,begc,endc,begp,endp)

     ALLOCATE(rcv_field(begg:endg), stat=ier)
     IF (ier /= 0) THEN
        CALL prism_abort_proto( ncomp_id, 'receive_fld_2pfl', 'Failure in allocating rcv_field' )
        RETURN
    ENDIF

 ! Default value (changed when coupling)
 pfl_coupled_sat(:) = .FALSE.
 pfl_coupled_psi(:) = .FALSE.
 nrcvinfo = OASIS_idle

CALL MPI_Barrier(kl_comm, ier )


 ! Scatter PFL soil moisture and pressure head to all nodes
 DO jn = 1, vrcv            !jn loop
   jnid = jn + 100 


   rcv_field(:) = 0.        
   IF( srcv(jnid)%laction ) THEN
      CALL oas_clm_rcv( jnid, isec, rcv_field(:),begg,endg, nrcvinfo(jn) )
      IF (srcv(jnid)%ref .eq. 'SAT') THEN
        pfl_coupled_sat(jn) = .TRUE.
        k(jn) = srcv(jnid)%level
      ELSE IF (srcv(jnid)%ref .eq. 'PSI') THEN
        pfl_coupled_psi(jn) = .TRUE.
        k(jn) = srcv(jnid)%level
      ENDIF
   ENDIF



   ! Get clumps
   nclumps = get_proc_clumps()
   DO nc = 1, nclumps     !Loop over clumps
     CALL get_clump_bounds(nc, begg, endg, begl, endl, begc, endc, begp, endp)
     
     ALLOCATE(filter_soilc(endc-begc+1), stat=ier)
     IF (ier /= 0) THEN
       CALL prism_abort_proto( ncomp_id, 'receive_fld_2pfl', 'Failure in allocating filter_soilc' )
       RETURN
     ENDIF
     

   num_soilc    = filter(nc)%num_soilc
   filter_soilc = filter(nc)%soilc
    
   IF (pfl_coupled_sat(jn)) THEN       !CMS From 1-10, values associated to soil water saturation
     DO fc =  1, num_soilc !begc, endc
       c = filter_soilc(fc)       
       IF (pfl_ice(c,k(jn)) <= 0.1_r8)  THEN   !DECOUPLING MOISTURE UPDATE FOR FROZEN SOILS CPS
         pfl_vol(c,k(jn)) = rcv_field(c)*watsat(c,k(jn))
         pfl_soi(c,k(jn)) = pfl_vol(c,k(jn))*dz(c,k(jn))*denh2o
       ENDIF
     ENDDO
   ELSE IF (pfl_coupled_psi(jn)) THEN  !CMS From 11-20, values associated to soil matrix potential
     DO fc =  1, num_soilc
       c = filter_soilc(fc)       
       IF (pfl_ice(c,k(jn)) <= 0.1_r8)  THEN   !DECOUPLING MOISTURE UPDATE FOR FROZEN SOILS CPS
         pfl_psi(c,k(jn)) = rcv_field(c)
       ENDIF
     ENDDO
   ENDIF
   ENDDO                 !Loop over clumps
   DEALLOCATE(filter_soilc) 
 END DO                  ! jn loop


!CPS WRITE(6, *)  "oasclm: receive_fld_2pfl: Receive subsurface field complete" 

deallocate(rcv_field)

 !  Wait coupling result (involved PE and non involved PE)
  CALL MPI_Barrier(kl_comm, ier )
 !CALL MPI_Bcast( lcoupled, 1, dummy, 0, kl_comm, ier )


!------------------------------------------------------------------------------
!- End of the Subroutine
!------------------------------------------------------------------------------

END SUBROUTINE receive_fld_2pfl
