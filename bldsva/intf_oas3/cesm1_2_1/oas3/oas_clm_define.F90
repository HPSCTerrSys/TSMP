SUBROUTINE oas_clm_define(SDATM)

!---------------------------------------------------------------------
! Description:
!  This routine sends invariant CLM3.5 fields to OASIS3 coupler:
!      1) grid information
!      2) domain decomposition
!      3) coupling variable definition
!
! Current Code Owner: TR32, Z4: Prabhakar Shrestha
!    phone: 0228733453
!    email: pshrestha@uni-bonn.de
!
! History:
! Version    Date       Name
! ---------- ---------- ----
! 1.1.1        2011/11/28 Prabhakar Shrestha 
!   Modfied and Implemented in CLM3.5, Initial release
! @VERSION@    @DATE@     <Your name>
! 1.1.1        2012/01/30 Mauro Sulis
! Definition of 10 sending fields CLM2PFL (source/sink term) and 
! 20 receiving PFL2CLM (water saturation and soil matrix potential)
! 1.2.1        2013/01/17 Markus Uebel, Prabhakar Shrestha
! Implementation of CO2 coupling
! 1.3.1        2013/01/31 Prabhakar Shrestha
! Implementation of aerodynamic resistance exchange
! Implementation of surface temperature and moisture exchange
! 2.1.0        2016/02/29 Prabhakar Shrestha
! Implementation for CESM 1.2.1
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!==============================================================================

! Declarations:
!
! Modules used:

USE oas_clm_vardef
USE mod_prism_grids_writing
USE mod_prism_def_partition_proto
USE shr_strdata_mod
USE mct_mod

!==============================================================================

IMPLICIT NONE

!==============================================================================

! Local Parmaters 

TYPE(shr_strdata_type), INTENT(IN)      :: SDATM

! Local Variables
INTEGER                    :: igrid                         ! ids returned by prism_def_grid

INTEGER                    :: var_nodims(2) ! 
INTEGER                    :: ipshape(2)    ! 
INTEGER, ALLOCATABLE       :: igparal(:)                    ! shape of arrays passed to PSMILe
INTEGER                    :: NULOUT=6

INTEGER                    :: k, ji, jj, jg, jgm1              ! local loop indicees
INTEGER ,POINTER           :: dmask(:)
INTEGER ,ALLOCATABLE       :: oas_mask(:)

CHARACTER(len=4)           :: clgrd='gclm'                  ! CPS

REAL(KIND=r8)              :: dx,dy
REAL(KIND=r8), ALLOCATABLE :: tmp_2D(:,:)                   ! Store Corners
REAL(KIND=r8), ALLOCATABLE :: zlon(:), zlat(:)
INTEGER                    :: write_aux_files
INTEGER                    :: rank, nprocs                    !CPScesm
!------------------------------------------------------------------------------
!- End of header
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!- Begin Subroutine oas_clm_define 
!------------------------------------------------------------------------------

 CALL MPI_Barrier(kl_comm, nerror)
 CALL MPI_Comm_Rank(kl_comm, rank, nerror)  !CPScesm 
 IF (nerror /= 0) CALL prism_abort_proto(ncomp_id, 'MPI_Comm_Rank', 'Failure in oas_clm_define') !CPScesm
 CALL MPI_Comm_Size ( kl_comm, nprocs, nerror )
 IF (nerror /= 0) CALL prism_abort_proto(ncomp_id, 'MPI_Comm_Size', 'Failure in oas_clm_define')

 IF (rank == 0) THEN
   ndlon = SDATM%nxg
   ndlat = SDATM%nyg
   PRINT*, "CPS DEBUG", ndlon, ndlat

   ALLOCATE( dmask(ndlon*ndlat), stat = nerror )
   IF ( nerror > 0 ) THEN
     CALL prism_abort_proto( ncomp_id, 'oas_clm_define', 'Failure in allocating dmask' )
     RETURN
   ENDIF

   k     = mct_aVect_indexRA(SDATM%grid%data,'mask')
   dmask = SDATM%grid%data%rAttr(k,:)
   PRINT*, "CPS DEBUG", MINVAL(dmask), MAXVAL(dmask)
   !
   !
   ALLOCATE( zlon(ndlon*ndlat), stat = nerror )
   IF ( nerror > 0 ) THEN
     CALL prism_abort_proto( ncomp_id, 'oas_clm_define', 'Failure in allocating zlon' )
     RETURN
   ENDIF
   ALLOCATE( zlat(ndlat*ndlon), stat = nerror )
   IF ( nerror > 0 ) THEN
     CALL prism_abort_proto( ncomp_id, 'oas_clm_define', 'Failure in allocating zlat' )
     RETURN
   ENDIF
   ALLOCATE( llmask(ndlon,ndlat,1), stat = nerror )
   IF ( nerror > 0 ) THEN
      CALL prism_abort_proto( ncomp_id, 'oas_clm_define', 'Failure in allocating llmask' )
      RETURN
   ENDIF
   ALLOCATE( tmp_2D(ndlon*ndlat,8), stat = nerror )
   IF ( nerror > 0 ) THEN
     CALL prism_abort_proto( ncomp_id, 'oas_clm_define', 'Failure in allocating tmp_2D' )
     RETURN
   ENDIF
   ALLOCATE(oas_mask(ndlon*ndlat))
   IF ( nerror > 0 ) THEN
      CALL prism_abort_proto( ncomp_id, 'oas_clm_define', 'Failure in allocating oas_mask' )
      RETURN
   ENDIF

      
   ! -----------------------------------------------------------------
   ! ... Define the elements, i.e. specify the corner points for each
   !     volume element. 
   !     We only need to give the 4 horizontal corners
   !     for a volume element plus the vertical position of the upper
   !     and lower face. Nevertheless the volume element has 8 corners.
   ! -----------------------------------------------------------------

   ! -----------------------------------------------------------------
   ! ... Define centers and corners 
   ! -----------------------------------------------------------------
   !  1: lower left corner. 2,: lower right corner.
   !  3: upper right corner. 4,: upper left corner.
   !
   ! Longitudes
   k     = mct_aVect_indexRA(SDATM%grid%data,'lon')
   zlon  = SDATM%grid%data%rAttr(k,ji)

   ! Latitudes
   k     = mct_aVect_indexRA(SDATM%grid%data,'lat')
   zlat  = SDATM%grid%data%rAttr(k,:)
    
   ! Grid Resolution 
   dx    = ABS(zlon(2)-zlon(1))
   dy    = ABS(zlat(ndlon+1)-zlat(1))  

   tmp_2D(:,1) = zlon - dx*0.5_r8
   tmp_2D(:,2) = zlon + dx*0.5_r8
   tmp_2D(:,3) = tmp_2D(:,2)
   tmp_2D(:,4) = tmp_2D(:,1) 

   tmp_2D(:,5) = zlat - dy*0.5_r8
   tmp_2D(:,6) = tmp_2D(:,5) 
   tmp_2D(:,7) = zlat + dy*0.5_r8
   tmp_2D(:,8) = tmp_2D(:,7) 

   ! -----------------------------------------------------------------
   ! ... Define the mask
   ! -----------------------------------------------------------------
   DO jj = 1, ndlat
   DO ji = 1, ndlon
     jg    = (jj-1)*ndlon + ji
     IF ( dmask(jg) == 1 ) THEN        
        llmask(ji,jj,1) = .TRUE.
        oas_mask(jg) = 0
     ELSE
        llmask(ji,jj,1) = .FALSE.
        oas_mask(jg) = 1
     ENDIF
   ENDDO
   ENDDO

   WRITE(nulout,*) 'oasclm: oas_clm_define: Land point number / total number', COUNT(llmask), ndlon*ndlat
   CALL flush(nulout)

   ! -----------------------------------------------------------------
   ! ... Write info on OASIS auxillary files (if needed)
   ! ----------------------------------------------------------------
   CALL prism_start_grids_writing (write_aux_files)
    
  IF ( write_aux_files == 1 ) THEN
 
 
     CALL prism_write_grid (clgrd, ndlon*ndlat, 1, zlon, zlat)

     CALL prism_write_corner (clgrd, ndlon*ndlat, 1, 4, tmp_2d(:,1:4), tmp_2d(:,5:8))

     CALL prism_write_mask (clgrd, ndlon*ndlat, 1, oas_mask)

     CALL prism_terminate_grids_writing()
        
  ENDIF                      ! write_aux_files = 1

  WRITE(nulout,*) ' oasclm: oas_clm_define: prism_terminate_grids'
  CALL flush(nulout)

 ENDIF                       ! masterproc


 IF (rank == 0) THEN !CPScesm
 ! -----------------------------------------------------------------
 ! ... Define the partition 
 ! -----------------------------------------------------------------
     
  ALLOCATE(igparal(4))
  ! Compute global offsets and local extents
  igparal(1) = 3              ! ORANGE style partition
  igparal(2) = 1              ! partitions number
  igparal(3) = 0              ! Global offset
  igparal(4) = ndlon*ndlat    ! Local extent

  CALL prism_def_partition_proto( igrid, igparal, nerror )
  IF( nerror /= PRISM_Success )   CALL prism_abort_proto (ncomp_id, 'oas_clm_define',   &
            &                                                        'Failure in prism_def_partition' )

  WRITE(nulout,*) ' oasclm: oas_clm_define: prism_def_partition' 
  CALL flush(nulout)

  ! -----------------------------------------------------------------
  ! ... Variable definition
  ! ----------------------------------------------------------------

  ! Default values
  ssnd(1:nmaxfld)%laction=.FALSE.  ; srcv(1:nmaxfld)%laction=.FALSE.

  ssnd(1)%clname='FSENDMD1'      !CPScesm dummy
!CPScesm  ssnd(1)%clname='CLM_TAUX'      !  zonal wind stress
  ssnd(2)%clname='CLM_TAUY'      !  meridional wind stress
  ssnd(3)%clname='CLMLATEN'      !  total latent heat flux (W/m**2)
  ssnd(4)%clname='CLMSENSI'      !  total sensible heat flux (W/m**2)
  ssnd(5)%clname='CLMINFRA'      ! emitted infrared (longwave) radiation (W/m**2)
  ssnd(6)%clname='CLMALBED'      ! direct albedo
  ssnd(7)%clname='CLMALBEI'      ! diffuse albedo
!MU (17.01.13)
  ssnd(8)%clname='CLMCO2FL'      ! net CO2 flux (now only photosynthesis rate) (umol CO2 m-2s-1)
!MU (17.01.13)
  ssnd(9)%clname='CLM_RAM1'      ! Aerodynamic resistance (s/m)   !CPS
  ssnd(10)%clname='CLM_RAH1'      ! Aerodynamic resistance (s/m)   !CPS
  ssnd(11)%clname='CLM_RAW1'      ! Aerodynamic resistance (s/m)   !CPS
  ssnd(12)%clname='CLM_TSF1'      ! Surface Temperature (K)   !CPS
  ssnd(13)%clname='CLM_QSF1'      ! Surface Humidity (kg/kg)   !CPS
!MU (12.04.13)
  ssnd(14)%clname='CLMPHOTO'      ! photosynthesis rate (umol CO2 m-2s-1)
  ssnd(15)%clname='CLMPLRES'      ! plant respiration (umol CO2 m-2s-1)
!MU (12.04.13)

  !CMS: from 101 to 200 are the sending fields from CLM to PFL
  ssnd(101)%clname='CLMFLX01'    !  evapotranspiration fluxes sent to PFL for each soil layer  
  ssnd(102)%clname='CLMFLX02'
  ssnd(103)%clname='CLMFLX03'
  ssnd(104)%clname='CLMFLX04'
  ssnd(105)%clname='CLMFLX05'
  ssnd(106)%clname='CLMFLX06'
  ssnd(107)%clname='CLMFLX07'
  ssnd(108)%clname='CLMFLX08'
  ssnd(109)%clname='CLMFLX09'
  ssnd(110)%clname='CLMFLX10'

  srcv(1)%clname='CLMTEMPE'
  srcv(2)%clname='CLMUWIND'
  srcv(3)%clname='CLMVWIND'
  srcv(4)%clname='CLMSPWAT'   ! specific water vapor content
  srcv(5)%clname='CLMTHICK'   ! thickness of lowest level (m)
  srcv(6)%clname='CLMPRESS'   ! surface pressure (Pa)
  srcv(7)%clname='CLMDIRSW'   ! direct shortwave downward radiation (W/m2)
  srcv(8)%clname='CLMDIFSW'   ! diffuse shortwave downward radiation (W/m2)
  srcv(9)%clname='CLMLONGW'   ! longwave downward radiation (W/m2)
  srcv(10)%clname='CLMCVRAI'  ! convective rain precipitation      (kg/m2*s)
  srcv(11)%clname='CLMCVSNW'  ! convective snow precipitation      (kg/m2*s)
  srcv(12)%clname='CLMGSRAI'  ! gridscale rain precipitation
  srcv(13)%clname='CLMGSSNW'  ! gridscale snow precipitation
  srcv(14)%clname='CLMGRAUP'  ! gridscale graupel precipitation
  srcv(15)%clname='CLMCVPRE'  ! total convective precipitation
  srcv(16)%clname='CLMGSPRE'  ! total gridscale precipitation
  srcv(17)%clname='CLMCO2PP'  ! CO2 partial pressure (Pa)  !CMU

  !CMS: from 101 to 200 are the receiving fields from PFL to CLM
  srcv(101)%clname='CLMSAT01' ! water saturation received from PFL for each soil layer
  srcv(102)%clname='CLMSAT02'
  srcv(103)%clname='CLMSAT03'
  srcv(104)%clname='CLMSAT04'   
  srcv(105)%clname='CLMSAT05'  
  srcv(106)%clname='CLMSAT06'
  srcv(107)%clname='CLMSAT07' 
  srcv(108)%clname='CLMSAT08'  
  srcv(109)%clname='CLMSAT09'   
  srcv(110)%clname='CLMSAT10'  

  srcv(101)%level= 1 ! # of soil layer
  srcv(102)%level= 2
  srcv(103)%level= 3
  srcv(104)%level= 4   
  srcv(105)%level= 5  
  srcv(106)%level= 6
  srcv(107)%level= 7 
  srcv(108)%level= 8  
  srcv(109)%level= 9   
  srcv(110)%level= 10  

  srcv(101:110)%ref='SAT'  

  srcv(111)%clname='CLMPSI01' ! pressure head received from PFL for each soil layer
  srcv(112)%clname='CLMPSI02'
  srcv(113)%clname='CLMPSI03'
  srcv(114)%clname='CLMPSI04'   
  srcv(115)%clname='CLMPSI05'  
  srcv(116)%clname='CLMPSI06'
  srcv(117)%clname='CLMPSI07' 
  srcv(118)%clname='CLMPSI08'  
  srcv(119)%clname='CLMPSI09'   
  srcv(120)%clname='CLMPSI10'  
  
  srcv(111)%level= 1 ! # of soil layer
  srcv(112)%level= 2
  srcv(113)%level= 3
  srcv(114)%level= 4   
  srcv(115)%level= 5  
  srcv(116)%level= 6
  srcv(117)%level= 7 
  srcv(118)%level= 8  
  srcv(119)%level= 9   
  srcv(120)%level= 10  

  srcv(111:120)%ref='PSI'
 
! Send/Receive Variable Selection
#ifdef COUP_OAS_COS

ssnd(1)%laction=.TRUE. !CPScesm
!CPScesm  IF (cpl_scheme) THEN         !CPS
!CPScesm     ssnd(5:7)%laction=.TRUE.
    !ssnd(6:7)%laction=.TRUE.     !CPS
!MU (12.04.13)
!CPScesm    ssnd(8)%laction=.TRUE.
!CPScesm    ssnd(14)%laction=.FALSE.
!CPScesm    ssnd(15)%laction=.FALSE.
!MU (12.04.13)
!CPScesm    ssnd(9:13)%laction=.TRUE.    !CPS
!CPScesm  ELSE
!CPScesm    ssnd(1:7)%laction=.TRUE.
!CPScesm    ssnd(8)%laction=.TRUE.
!CPScesm    ssnd(14)%laction=.FALSE.
!CPScesm    ssnd(15)%laction=.FALSE. 
!CPScesm  ENDIF                         !CPS

!CPScesm  srcv(1:9)%laction=.TRUE.
!CPScesm  srcv(15:16)%laction=.TRUE. ! Coupling only total convective and gridscale precipitations 
!MU (17.01.13)
!CPScesm  srcv(17)%laction=.TRUE.    ! always true
!  IF (srcv(17)%laction==.FALSE.) THEN
!    PRINT*, 'ERROR (oas_clm_define): srcv(17)%laction has to be .TRUE.'
!    PRINT*, '----- If CO2 is not initialized in COSMO a dummy is sent to CLM and CO2 content from CLM is used.'
!  ENDIF
!MU (17.01.13)
#endif

#ifdef COUP_OAS_PFL
! CMS From 101 to 110 the ssnd variables are defined for ParFlow
  ssnd(101:110)%laction=.TRUE.
! CMS From 101 to 120 the srcv variables are defined for parflow
  srcv(101:110)%laction=.TRUE.
  srcv(111:120)%laction=.TRUE.
#endif

  var_nodims(1) = 1           ! Dimension number of exchanged arrays
  var_nodims(2) = 1           ! number of bundles (always 1 for OASIS3)

  ipshape(1) = 1             ! minimum index for each dimension of the coupling field array
  ipshape(2) = ndlon*ndlat   ! maximum index for each dimension of the coupling field array

  ! ... Announce send variables. 
  !
!      ksnd=0                           !CPS
  DO ji = 1, nmaxfld
    IF ( ssnd(ji)%laction ) THEN 
            
      CALL prism_def_var_proto(ssnd(ji)%nid, ssnd(ji)%clname, igrid, &
                                     var_nodims, PRISM_Out, ipshape, PRISM_Real, nerror )
      IF ( nerror /= PRISM_Success )   CALL prism_abort_proto( ssnd(ji)%nid, 'oas_clm_define',   &
               &                                               'Failure in prism_def_var for '//TRIM(ssnd(ji)%clname))
 !           ksnd = ksnd + 1             !CPS now defined in oas_clm_vardef
    ENDIF
  END DO
  !
  ! ... Announce received variables. 
  !
  DO ji = 1, nmaxfld
    IF ( srcv(ji)%laction ) THEN 

    CALL prism_def_var_proto(srcv(ji)%nid, srcv(ji)%clname, igrid, &
                         var_nodims, PRISM_In, ipshape, PRISM_Real, nerror )
    IF ( nerror /= PRISM_Success )   CALL prism_abort_proto( srcv(ji)%nid, 'oas_clm_define',   &
               &                                               'Failure in prism_def_var for '//TRIM(srcv(ji)%clname))
    ENDIF
  END DO
      
  !
  ! ... Allocate memory for data exchange and initilize it
  !
  ALLOCATE( exfld(ndlon, ndlat, krcv), stat = nerror )
  IF ( nerror > 0 ) THEN
     CALL prism_abort_proto( ncomp_id, 'oas_cos_define', 'Failure in allocating exfld' )
     RETURN
  ENDIF
     
  exfld = -99999._r8

  !------------------------------------------------------------------
  ! End of definition phase
  !------------------------------------------------------------------
  DEALLOCATE( tmp_2D, zlon, zlat, igparal)

  ! must be done by coupled processors only (OASIS3)
  CALL prism_enddef_proto( nerror )

  IF ( nerror /= PRISM_Success )   CALL prism_abort_proto ( ncomp_id, 'oas_clm_define', 'Failure in prism_enddef')

 ENDIF                 !masterproc
      
 CALL MPI_Barrier(kl_comm, nerror)
 WRITE(nulout,*) 'oasclm: oas_clm_define: prism enddef' 
 CALL flush(nulout)

!------------------------------------------------------------------------------
!- End of the Subroutine
!------------------------------------------------------------------------------

END SUBROUTINE oas_clm_define
