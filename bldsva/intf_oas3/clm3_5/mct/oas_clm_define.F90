SUBROUTINE oas_clm_define(filenam)

!---------------------------------------------------------------------
! Description:
!  This routine sends invariant CLM3.5 fields to OASIS3 coupler:
!      1) grid information
!      2) domain decomposition
!      3) coupling variable definition
!
! References:
!  CEREFACS/ETH: E. Maisonnave, Edoward Davin
!
! Current Code Owner: TR32, Z4: Prabhakar Shrestha
!    phone: 0228733453
!    email: pshrestha@uni-bonn.de
!
! History:
! Version    Date       Name
! ---------- ---------- ----
! 1.1        2011/11/28 Prabhakar Shrestha 
!   Modfied and Implemented in CLM3.5, Initial release
! @VERSION@    @DATE@     <Your name>
! 1.1        2012/01/30 Mauro Sulis
! Definition of 10 sending fields CLM2PFL (source/sink term) and 
! 20 receiving PFL2CLM (water saturation and soil matrix potential)
! 2.1        2013/01/17 Markus Uebel, Prabhakar Shrestha
! Implementation of CO2 coupling
! 3.1        2013/01/31 Prabhakar Shrestha
! Implementation of aerodynamic resistance exchange
! Implementation of surface temperature and moisture exchange
! 4.1        2013/10/31 Fabian Gasper
! Implementation of a parallel coupling   
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!==============================================================================

! Declarations:
!
! Modules used:

USE oas_clm_vardef

!USE mod_prism_grid
!USE mod_prism_def_partition_proto

USE domainMod   , ONLY : latlon_type, latlon_check, alatlon,amask , adomain
USE spmdMod     , ONLY : masterproc , mpicom, iam
USE surfrdMod   , ONLY : surfrd_get_latlon
USE decompMod    , only : get_proc_global, get_proc_bounds, adecomp
!==============================================================================

IMPLICIT NONE

!==============================================================================

! Local Parmaters 

CHARACTER(len=256), INTENT(IN) :: filenam             !CPS clm grid data

! Local Variables

Integer :: begg,endg, begl,endl,begc,endc,begp,endp,c,off,leng,it1,i,j,im1,jm1
INTEGER              :: igrid                         ! ids returned by prism_def_grid
INTEGER              :: iptid                         ! ids returned by prism_set_points

INTEGER              :: imskid                        ! ids returned by prism_set_mask

INTEGER              :: iextent(1,3)  ! 
INTEGER              :: ioffset(1,3)  ! 
INTEGER              :: var_nodims(2) ! 
INTEGER              :: ipshape(2)    ! 

INTEGER, ALLOCATABLE :: igparal(:)                    ! shape of arrays passed to PSMILe
INTEGER              :: NULOUT=6

INTEGER              :: ji, jj, jg, jgm1              ! local loop indicees
INTEGER ,ALLOCATABLE :: oas_mask(:)

CHARACTER(len=4)     :: clgrd='gclm'                  ! CPS

REAL(KIND=r8)        :: dx
REAL(KIND=r8), ALLOCATABLE :: zclo(:,:), zcla(:,:), tmp_2D(:,:)
REAL(KIND=r8), ALLOCATABLE :: zlon(:), zlat(:)
INTEGER                    :: write_aux_files

Integer :: c_comm



integer :: ai, aj ,ani , anj , an , owner, last_owner



!------------------------------------------------------------------------------
!- End of header
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!- Begin Subroutine oas_clm_define 
!------------------------------------------------------------------------------

! Define coupling scheme between COSMO and CLM
#ifdef CPL_SCHEME_F
 cpl_scheme = .True. !TRN Scheme
#else
 cpl_scheme = .False. !INV Scheme
#endif


 ndlon = alatlon%ni
 ndlat = alatlon%nj

 

 ! For the moment OASIS4 not able to process 2D subdomain with non contiguous unmasked points (land points) 
 ! (12/2010)
 ! To compare OASIS3/OASIS4 behaviours, we do the same with OASIS3
 !
 CALL MPI_Barrier(kl_comm, nerror)
#define WITHAUXFILES
#ifdef WITHAUXFILES
 IF ( masterproc ) THEN
   !
   CALL prism_start_grids_writing (write_aux_files)
   !
   ALLOCATE( zclo(ndlon, 4), stat = nerror )
   IF ( nerror > 0 ) THEN
     CALL prism_abort_proto( ncomp_id, 'oas_clm_define', 'Failure in allocating zclo' )
     RETURN
   ENDIF
   ALLOCATE( zcla(ndlat, 4), stat = nerror )
   IF ( nerror > 0 ) THEN
     CALL prism_abort_proto( ncomp_id, 'oas_clm_define', 'Failure in allocating zcla' )
     RETURN
   ENDIF
   ALLOCATE( zlon(ndlon), stat = nerror )
   IF ( nerror > 0 ) THEN
     CALL prism_abort_proto( ncomp_id, 'oas_clm_define', 'Failure in allocating zlon' )
     RETURN
   ENDIF
   ALLOCATE( zlat(ndlat), stat = nerror )
   IF ( nerror > 0 ) THEN
     CALL prism_abort_proto( ncomp_id, 'oas_clm_define', 'Failure in allocating zlat' )
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
   !  using latlon%edges :  global edges (N,E,S,W)
   !
   ! Latitudes
   ! Assumes lats are constant on an i line
   !
    zlat(:) = alatlon%latc(1:ndlat)
   !
   IF (alatlon%latc(ndlat) > alatlon%latc(1)) THEN  ! South to North grid
      zcla(1,1) = alatlon%edges(3)
      zcla(ndlat,2) = alatlon%edges(1)
    !
      DO jj = 2, ndlat
         jgm1 = jj - 1
         zcla(jj,1) = ( alatlon%latc(jj) + alatlon%latc(jgm1) ) * 0.5_r8
         zcla(jgm1,2) = zcla(jj,1)
      ENDDO 

   ELSE                                      ! North to South grid
      zcla(ndlat,1) = alatlon%edges(3)
      zcla(1,2) = alatlon%edges(1)
      !
      DO jj = 2, ndlat
        jgm1 = jj - 1
        zcla(jj,2) = ( alatlon%latc(jj) + alatlon%latc(jgm1) ) * 0.5_r8
        zcla(jgm1,1) = zcla(jj,2)
      ENDDO

   ENDIF

   ! Longitudes
   !
   zlon(:) = alatlon%lonc(1:ndlon)
   !
   zclo(:,1) = alatlon%edges(4)
   zclo(:,2) = alatlon%edges(2)

   dx = zlon(2)-zlon(1)

   DO ji = 2, ndlon
      zclo(ji,1) = zclo(ji,1) + ( ji -1 ) * dx
      zclo(ji-1,2) = zclo(ji,1)
   ENDDO 

   DO ji = 1, ndlon  ! make sure that all centers are nearby corners regarding the longitude...
     IF( zlon(ji) - zclo(ji,1) < -300._r8 )   zlon(ji) = zlon(ji) + 360._r8
     IF( zlon(ji) - zclo(ji,1) >  300._r8 )   zlon(ji) = zlon(ji) - 360._r8
   END DO


   ! -----------------------------------------------------------------
   ! ... Define the mask
   ! -----------------------------------------------------------------
   DO jj = 1, ndlat
   DO ji = 1, ndlon
     jg    = (jj-1)*ndlon + ji
     IF ( amask(jg) == 1 ) THEN
        oas_mask(jg) = 0
     ELSE
        oas_mask(jg) = 1
     ENDIF
   ENDDO
   ENDDO

   CALL flush(nulout)

   ! -----------------------------------------------------------------
   ! ... Write info on OASIS auxillary files (if needed)
   ! ----------------------------------------------------------------

   IF ( write_aux_files == 1 ) THEN
 
    DO jj = 1, ndlat
    DO ji = 1, ndlon
       tmp_2D(ji+(jj-1)*ndlon,1) = zlon(ji)
       tmp_2D(ji+(jj-1)*ndlon,2) = zlat(jj)
    ENDDO
    ENDDO

   CALL prism_write_grid (clgrd, ndlon*ndlat, 1, tmp_2D(:,1), tmp_2D(:,2))

    DO jj = 1, ndlat
    DO ji = 1, ndlon
     tmp_2D(ji+(jj-1)*ndlon,1) = zclo(ji,1)
     tmp_2D(ji+(jj-1)*ndlon,2) = zclo(ji,2)
     tmp_2D(ji+(jj-1)*ndlon,5) = zcla(jj,1)
     tmp_2D(ji+(jj-1)*ndlon,7) = zcla(jj,2)
    ENDDO
    ENDDO

    ! Fill missing corners for longitude
    tmp_2D(:,3) = tmp_2D(:,2)
    tmp_2D(:,4) = tmp_2D(:,1)

    ! Fill missing corners for latitude
    tmp_2D(:,6) = tmp_2D(:,5)
    tmp_2D(:,8) = tmp_2D(:,7)

    CALL prism_write_corner (clgrd, ndlon*ndlat, 1, 4, tmp_2D(:,1:4), tmp_2D(:,5:8))

    ! tbd CALL prism_write_angle (clgrd, ndlon*ndlat, 1, angle)

    CALL prism_write_mask (clgrd, ndlon*ndlat, 1, oas_mask)

    ! EM tbd Fill areas with true values
    tmp_2D(:,1) = 1 
!CPS        CALL prism_write_area (clgrd, ndlon*ndlat, 1, tmp_2D(:,1))

    CALL prism_terminate_grids_writing()
        
  ENDIF                      ! write_aux_files = 1

  WRITE(nulout,*) ' oasclm: oas_clm_define: prism_terminate_grids'
  CALL flush(nulout)
  DEALLOCATE(tmp_2D, oas_mask, zclo, zcla, zlon, zlat)
 ENDIF                       ! masterproc
#endif

!FG  The following piece of code querries the number, length and offset of 
!    contiguous parts in the lat/lon representation. This is done in a 
!    brutforce way by travercing through begg-endg and determining if
!    predecessor in begg-endg is also predecessor in lat/lon.
!    I couldn't find a CLM intrinsic way to get this information, if there is
!    one, please let me know. 
!    At the moment the orange partition of oasis supports only 200 parts,
!    luckily CLM seem to use always 20 parts regardless of the domainsize.

ALLOCATE(igparal(200))


igparal(1) = 3  !orange partition

#if (1 == 0)
CALL get_proc_bounds(begg,endg,begl,endl,begc,endc,begp,endp)

c=1
off=(adecomp%gdc2j(begg)-1)*ndlon +adecomp%gdc2i(begg) -1
leng=1

do it1=begg+1,endg

  i=adecomp%gdc2i(it1)
  im1=adecomp%gdc2i(it1-1)
  j=adecomp%gdc2j(it1)
  if(j.eq.1) then
        jm1=-1
  else
    jm1=adecomp%gdc2j(it1-1)
  endif

  if( (i-1.eq.im1) .or. (i.eq.1 .and. jm1.eq.j-1 .and. im1.eq.ndlon)  ) then
    leng=leng+1
  else
    igparal(c*2+1) = off
    igparal(c*2+2) = leng
    igparal(2) = c
    off=(adecomp%gdc2j(it1)-1)*ndlon +adecomp%gdc2i(it1) -1
    c=c+1
    leng=1
  endif

enddo
igparal(c*2+1) = off
igparal(c*2+2) = leng
igparal(2) = c

#endif

ani = adomain%ni
anj = adomain%nj

total_part_len=0
c=0
leng=0
last_owner=-1
do aj = 1,anj
  do ai = 1,ani
    an = (aj-1)*ani + ai
    owner = adecomp%glo2owner(an)
    if(owner==iam ) then
       if(last_owner/=iam) then       !start of new partiontion
          if(last_owner/=-1) leng=0   !no leading masked cells have to be added to the first partition
          c=c+1
          if(last_owner==-1) then
             igparal(c*2+1) = 0
          else
             igparal(c*2+1) = an-1
          endif
       endif
       leng=leng+1
    else  
       if(owner==-1) then
          if(last_owner== iam .or. last_owner== -1) leng=leng+1 !adding also masked cells to the partition
       else                                                     !other procs partiton start
          if(last_owner==iam) then                               ! if this cell ends ongoing partition, length is written out
             igparal(c*2+2) = leng            
             total_part_len=total_part_len+leng
          endif
       endif 
    endif    
    if(owner/=-1) last_owner=owner
  enddo
enddo
if(last_owner==iam)then   !write length if this proc had last partition
  igparal(c*2+2) = leng  
  total_part_len=total_part_len+leng
endif
igparal(2)=c

CALL MPI_Barrier(kl_comm, nerror)

  ! -----------------------------------------------------------------
  ! ... Define the partition 
  ! -----------------------------------------------------------------
     
        
  CALL prism_def_partition_proto( igrid, igparal, nerror )
  IF( nerror /= PRISM_Success )   CALL prism_abort_proto (ncomp_id, 'oas_clm_define',   &
            &                                                        'Failure in prism_def_partition' )




  ! -----------------------------------------------------------------
  ! ... Variable definition
  ! ----------------------------------------------------------------

  ! Default values
  ssnd(1:nmaxfld)%laction=.FALSE.  ; srcv(1:nmaxfld)%laction=.FALSE.

  ssnd(1)%clname='CLM_TAUX'      !  zonal wind stress
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

IF (cpl_scheme) THEN         !CPS
  ssnd(5:7)%laction=.TRUE.     !CPS
  ssnd(8)%laction=.TRUE.
  ssnd(14)%laction=.FALSE.
  ssnd(15)%laction=.FALSE.
  ssnd(9:13)%laction=.TRUE.    !CPS
ELSE
  ssnd(1:7)%laction=.TRUE.
  ssnd(8)%laction=.TRUE.
  ssnd(14)%laction=.FALSE.
  ssnd(15)%laction=.FALSE.
ENDIF 

  srcv(1:9)%laction=.TRUE.
  srcv(15:16)%laction=.TRUE. ! Coupling only total convective and gridscale precipitations 
!MU (17.01.13)
  srcv(17)%laction=.TRUE.    ! always true
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


  !------------------------------------------------------------------
  ! End of definition phase
  !------------------------------------------------------------------

DEALLOCATE( igparal)


  CALL prism_enddef_proto( nerror )
  IF ( nerror /= PRISM_Success )   CALL prism_abort_proto ( ncomp_id, 'oas_clm_define', 'Failure in prism_enddef')


 CALL MPI_Barrier(kl_comm, nerror)
 WRITE(nulout,*) 'oasclm: oas_clm_define: prism enddef' 
 CALL flush(nulout)

!------------------------------------------------------------------------------
!- End of the Subroutine
!------------------------------------------------------------------------------

END SUBROUTINE oas_clm_define
