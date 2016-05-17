SUBROUTINE oas_cos_define

!---------------------------------------------------------------------
! Description:
!  This routine defines grids, domain decomposition and coupling variables
!  for the OASIS3 coupler
!
! Current Code Owner: TR32, Z4: Prabhakar Shrestha
!    phone: 0228733453
!    email: pshrestha@uni-bonn.de
!
! History:
! Version    Date       Name
! ---------- ---------- ----
! 1.1.1        2011/11/28 Prabhakar Shrestha 
!   Modfied and Implemented in COSMO4.11, Initial release
! 1.2.1        2013/01/17 Markus Uebel 
!   CO2 coupling (photosynthesis rate) included and implemented in COSMO4.21
! 1.2.2        2013/12/30 Prabhakar Shrestha
!   Add readclm for masked simulations
! @VERSION@    @DATE@     <Your name>
!  <Modification comments>         
! 1.3.1        2013/01/31 Prabhakar Shrestha
!   Added new receive variables from CLM
!     aerodynamic resistance for heat, moisture, momentum
!     surface temperature and humidity
! 1.3.2        2014/01/05 Prabhakar Shrestha
!   Added readclm flag for idealized masked simulation
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!==============================================================================

! Declarations:
!
! Modules used:

USE data_parameters,  ONLY : wp, iintegers
USE oas_cos_vardef
USE mod_prism_grids_writing
USE mod_prism_def_partition_proto
USE data_modelconfig, ONLY : ie, je, dlon, dlat, startlon_tot, startlat_tot, raddeg,  &
                                   ie_tot, je_tot
USE data_fields,      ONLY : rlat, rlon, fr_land
USE data_parallel,    ONLY :   &
                            my_cart_id,      & ! rank of this subdomain in the cartesian communicator
                            isubpos,         & ! positions of the subdomains in the total domain. Given
                                               ! are the i- and the j-indices of the lower left and the
                                               ! upper right grid point in the order
                                               !                  i_ll, j_ll, i_ur, j_ur.
                                               ! Only the interior of the domains are considered, not
                                               ! the boundary lines.
                            nboundlines,     & ! number of boundary lines of the domain for which
                                               ! no forecast is computed = overlapping boundary
                                               ! lines of the subdomains
                            num_compute        ! Total number of PE's CPS
USE parallel_utilities,       ONLY :         &
                            gather_field,    & ! gathers the parts of a total field from all subdomains
                           distribute_field    ! distribute total field to PEs
USE netcdf
!==============================================================================

IMPLICIT NONE

!==============================================================================

! Local Variable

INTEGER(KIND=iintegers)         :: igrid     ! ids returned by prism_def_grid
INTEGER(KIND=iintegers)         :: iptid     ! ids returned by prism_set_points

INTEGER(KIND=iintegers)         :: imskid    ! ids returned by prism_set_mask

INTEGER(KIND=iintegers)         :: iextent(1,3)    ! 
INTEGER(KIND=iintegers)         :: ioffset(1,3)    ! 

INTEGER(KIND=iintegers)         :: ji, jj, jg
INTEGER(KIND=iintegers)         :: jg_end, jh_beg     ! index
INTEGER(KIND=iintegers)         :: NULOUT=6
INTEGER(KIND=iintegers)         :: id_part    !
INTEGER(KIND=iintegers)         :: write_aux_files    !
INTEGER(KIND=iintegers)         :: paral(5)       ! OASIS3 box partition
INTEGER(KIND=iintegers)         :: var_nodims(2)       ! OASIS3 box partition
INTEGER(KIND=iintegers)         :: ishape(2,2) ! shape of arrays passed to PSMILe


CHARACTER(len=4)                :: clgrd

REAL(KIND=wp), DIMENSION(ie,je,4)             :: zclo, zcla
REAL(KIND=wp), DIMENSION(ie,je)               :: zlon, zlat
!CPS REAL(KIND=wp), DIMENSION(ie,je)               :: zmask
INTEGER(KIND=iintegers), DIMENSION(ie_tot,je_tot) :: imask_tot
REAL(KIND=wp), DIMENSION(ie_tot,je_tot)       :: zmask_tot
REAL(KIND=wp), DIMENSION(ie_tot,je_tot,4)     :: zclo_tot, zcla_tot
REAL(KIND=wp), DIMENSION(ie_tot,je_tot)       :: zlon_tot, zlat_tot

INTEGER(KIND=iintegers)         :: nlei_tot, nlej_tot  ! upper halo limits on a global subdomain
INTEGER(KIND=iintegers)         :: jih_tot, jjh_tot    ! global subdomain size without halo
REAL(KIND=wp)               :: start_lonc, start_latc

INTEGER                         :: readclm = 0         ! 1 or 0 to read clm mask
INTEGER                         :: status, cosncid, cosvarid(7)

!------------------------------------------------------------------------------
!- End of header
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!- Begin Subroutine oas_cos_define 
!------------------------------------------------------------------------------

 ! Define coupling scheme between COSMO and CLM
#ifdef CPL_SCHEME_F
 cpl_scheme = .True. !TRN Scheme
#else
 cpl_scheme = .False. !INV Scheme
#endif

 !  array size without halo
 !
 ! Local
 jih = isubpos(my_cart_id,3) - isubpos(my_cart_id,1) + 1
 jjh = isubpos(my_cart_id,4) - isubpos(my_cart_id,2) + 1
 !  
 ! Global
 jih_tot = ie_tot - 2 * nboundlines
 jjh_tot = je_tot - 2 * nboundlines

 ! halo limits
 !
 ! Local
 nldi = 1 + nboundlines
 nldj = 1 + nboundlines
 nlei = ie - nboundlines
 nlej = je - nboundlines
 !
 ! Global
 nlei_tot = ie_tot - nboundlines
 nlej_tot = je_tot - nboundlines

 ! First determine if any land point in the domain
 ! If not, process should not be involved in the coupling
 !
 ! CPS Allocate zmask (zmask =0 for llandmask = True, fr_land>=0.5) 
 ALLOCATE( zmask(ie, je), stat = nerror )
 IF ( nerror > 0 ) THEN
   CALL prism_abort_proto(ncomp_id,'oas_cos_define','Failure in allocating zmask')
   RETURN
 ENDIF

 zmask(:,:) = 0.

 WHERE (fr_land(nldi:nlei, nldj:nlej) < 0.5_wp) zmask(nldi:nlei, nldj:nlej) = 1.
 !
 IF ( COUNT( zmask(nldi:nlei, nldj:nlej) < 0.5) < 1 ) THEN
   WRITE(NULOUT,*) ' oas_cos_define : All PE grid points masked '
   CALL FLUSH(NULOUT)
   lpe_cpl = .FALSE.
 ELSE
   lpe_cpl = .TRUE.
 ENDIF
 !
 ! With OASIS3, all domains involved
  lpe_cpl = .TRUE.

 ! -----------------------------------------------------------------
 ! ... Define the elements, i.e. specify the corner points for each
 !     volume element. 
 ! -----------------------------------------------------------------
      
 ! -----------------------------------------------------------------
 ! ... Define the center points
 ! -----------------------------------------------------------------

 ! conversion to degrees
 zlon = rlon * raddeg
 zlat = rlat * raddeg

 ! -----------------------------------------------------------------
 ! ... Define the corners
 ! -----------------------------------------------------------------
 !  1: lower left corner. 2: lower right corner 
 !  3: upper right corner 4: upper left corner
 !  (anti-clockwise, starting from the bottom left corner)
 !
 start_lonc = startlon_tot - dlon * 0.5_wp
 start_latc = startlat_tot - dlat * 0.5_wp

 CALL constant_corners(zclo(:,:,1), zcla(:,:,1), start_lonc, start_latc)

 start_lonc = startlon_tot + dlon * 0.5_wp
 start_latc = startlat_tot - dlat * 0.5_wp

 CALL constant_corners(zclo(:,:,2), zcla(:,:,2), start_lonc, start_latc)

 start_lonc = startlon_tot + dlon * 0.5_wp
 start_latc = startlat_tot + dlat * 0.5_wp

 CALL constant_corners(zclo(:,:,3), zcla(:,:,3), start_lonc, start_latc)

 start_lonc = startlon_tot - dlon * 0.5_wp
 start_latc = startlat_tot + dlat * 0.5_wp

 CALL constant_corners(zclo(:,:,4), zcla(:,:,4), start_lonc, start_latc)

 zclo = zclo * raddeg
 zcla = zcla * raddeg

 ! Gather information
 CALL gather_field(zlon,ie,je,zlon_tot,ie_tot,je_tot,0,nerror) 
 CALL gather_field(zlat,ie,je,zlat_tot,ie_tot,je_tot,0,nerror) 
 CALL gather_field(zmask,ie,je,zmask_tot,ie_tot,je_tot,0,nerror) 
 DO jg =1, 4
    CALL gather_field(zclo(:,:,jg),ie,je,zclo_tot(:,:,jg),ie_tot,je_tot,0,nerror) 
    CALL gather_field(zcla(:,:,jg),ie,je,zcla_tot(:,:,jg),ie_tot,je_tot,0,nerror) 
 ENDDO

 ! -----------------------------------------------------------------
 ! ... Master process writes info on OASIS3 auxillary files (if needed) 
 ! ----------------------------------------------------------------

!CPS      CALL prism_start_grids_writing (write_aux_files)

!! CPS READ CLM GRID AND MASK PARAMETERS
!! CPS This will overwride COSMO grids and mask, these grids will be sent to
!! OASIS. This is done to make consistent variable exchange without
!! interpolation between the COSMO grid and CLM pseduo-atmosphere grid.

  IF (my_cart_id == 0) THEN
     imask_tot=INT(zmask_tot)
  ENDIF

  IF (readclm == 1) THEN
  IF (my_cart_id == 0) THEN
     
     WRITE(nulout,*) "oascos: oas_cos_define: Reading clmgrid ....."
     status = nf90_open("clmgrid.nc", NF90_NOWRITE, cosncid)
     status = nf90_inq_varid(cosncid, "LONGXY" , cosvarid(1))
     status = nf90_inq_varid(cosncid, "LATIXY" , cosvarid(2))
     status = nf90_inq_varid(cosncid, "LONW" , cosvarid(3))
     status = nf90_inq_varid(cosncid, "LONE" , cosvarid(4))
     status = nf90_inq_varid(cosncid, "LATS" , cosvarid(5))
     status = nf90_inq_varid(cosncid, "LATN" , cosvarid(6))
     status = nf90_inq_varid(cosncid, "LANDMASK" , cosvarid(7))
     status = nf90_get_var(cosncid, cosvarid(1), zlon_tot(nldi:nlei_tot,nldj:nlej_tot))
     status = nf90_get_var(cosncid, cosvarid(2), zlat_tot(nldi:nlei_tot,nldj:nlej_tot))

     status = nf90_get_var(cosncid, cosvarid(3), zclo_tot(nldi:nlei_tot,nldj:nlej_tot,1))
     status = nf90_get_var(cosncid, cosvarid(4), zclo_tot(nldi:nlei_tot,nldj:nlej_tot,2))
     status = nf90_get_var(cosncid, cosvarid(4), zclo_tot(nldi:nlei_tot,nldj:nlej_tot,3))
     status = nf90_get_var(cosncid, cosvarid(3), zclo_tot(nldi:nlei_tot,nldj:nlej_tot,4))
       
     status = nf90_get_var(cosncid, cosvarid(5), zcla_tot(nldi:nlei_tot,nldj:nlej_tot,1))
     status = nf90_get_var(cosncid, cosvarid(5), zcla_tot(nldi:nlei_tot,nldj:nlej_tot,2))
     status = nf90_get_var(cosncid, cosvarid(6), zcla_tot(nldi:nlei_tot,nldj:nlej_tot,3))
     status = nf90_get_var(cosncid, cosvarid(6), zcla_tot(nldi:nlei_tot,nldj:nlej_tot,4)) 

     status = nf90_get_var(cosncid, cosvarid(7), imask_tot(nldi:nlei_tot,nldj:nlej_tot))
     status = nf90_close(cosncid)
     ! maskland = 0  
     zmask_tot(nldi:nlei_tot,nldj:nlej_tot) = ABS(imask_tot(nldi:nlei_tot,nldj:nlej_tot)-1)
     imask_tot(nldi:nlei_tot,nldj:nlej_tot) = ABS(imask_tot(nldi:nlei_tot,nldj:nlej_tot)-1)

  ENDIF           ! my_cart_id .eq. 0

    ! CPS Distribute to all PE's, this should be called by all PE's 
    IF (num_compute == 1) THEN
      zmask = zmask_tot
    ELSE
      CALL distribute_field (zmask_tot(:,:), ie_tot, je_tot,           &
                             zmask(:,:), ie,     je, 0,    nerror)
    ENDIF
    ! CPS Distribute

 ENDIF            ! readclm == 1


 IF ( my_cart_id == 0 ) THEN


  IF ( IOASISDEBUGLVL == 0 ) THEN
   WRITE(nulout,*) ' dx dy  ', zlon_tot(2,1)-zlon_tot(1,1), zlat_tot(1,2)-zlat_tot(1,1) 
   WRITE(nulout,*) ' zlon   ', MINVAL(zlon_tot(nldi:nlei_tot, nldj:nlej_tot)), MAXVAL(zlon_tot(nldi:nlei_tot, nldj:nlej_tot))
   WRITE(nulout,*) ' zlat   ', MINVAL(zlat_tot(nldi:nlei_tot, nldj:nlej_tot)), MAXVAL(zlat_tot(nldi:nlei_tot, nldj:nlej_tot))
   WRITE(nulout,*) ' zclo ll', MINVAL(zclo_tot(nldi:nlei_tot, nldj:nlej_tot,1)), MAXVAL(zclo_tot(nldi:nlei_tot, nldj:nlej_tot,1))
   WRITE(nulout,*) ' zclo ur', MINVAL(zclo_tot(nldi:nlei_tot, nldj:nlej_tot,3)), MAXVAL(zclo_tot(nldi:nlei_tot, nldj:nlej_tot,3))
   WRITE(nulout,*) ' zcla ll', MINVAL(zcla_tot(nldi:nlei_tot, nldj:nlej_tot,1)), MAXVAL(zcla_tot(nldi:nlei_tot, nldj:nlej_tot,1))
   WRITE(nulout,*) ' zcla ur', MINVAL(zcla_tot(nldi:nlei_tot, nldj:nlej_tot,3)), MAXVAL(zcla_tot(nldi:nlei_tot, nldj:nlej_tot,3))
   WRITE(nulout,*) ' zmask   ', MINVAL(zmask_tot(nldi:nlei_tot, nldj:nlej_tot)), MAXVAL(zmask_tot(nldi:nlei_tot, nldj:nlej_tot))
   CALL flush(nulout)
 ENDIF



   CALL prism_start_grids_writing (write_aux_files)               !CPS
   clgrd = 'gcos'

   IF ( write_aux_files == 1 ) THEN

     CALL prism_write_grid (clgrd, jih_tot, jjh_tot,            &
                        zlon_tot(nldi:nlei_tot, nldj:nlej_tot), &
                                 zlat_tot(nldi:nlei_tot, nldj:nlej_tot))

     CALL prism_write_corner (clgrd, jih_tot, jjh_tot, 4,       &
                      zclo_tot(nldi:nlei_tot, nldj:nlej_tot,:), &
                      zcla_tot(nldi:nlei_tot, nldj:nlej_tot,:))

! tbd CALL prism_write_angle (clgrd, jih_tot, jjh_tot, angle)

     CALL prism_write_mask (clgrd, jih_tot, jjh_tot,            &
                      imask_tot(nldi:nlei_tot, nldj:nlej_tot))

!CPS          CALL prism_write_area (clgrd, jih_tot, jjh_tot, &
!CPS                                 zlat_tot(nldi:nlei_tot, nldj:nlej_tot)*0.+1.)

     CALL prism_terminate_grids_writing()

     IF ( IOASISDEBUGLVL > 0 ) THEN
        WRITE(nulout,*) " oascos: oas_cos_define: prism_terminate_grids_writing"
        CALL flush(nulout)
     ENDIF

   ENDIF                 ! write_aux_files == 1
 ENDIF                   ! my_cart_id == 0

 ! Wait OASIS signal to continue simulation
!CPS      IF ( write_aux_files == 1 ) call MPI_Barrier(mpi_comm, nerror)
 CALL MPI_Barrier(kl_comm, nerror)

 IF ( lpe_cpl ) THEN

   ! ... Allocate memory for data exchange and partition defintion
   !
   ALLOCATE( exfld(nldi:nlei, nldj:nlej), stat = nerror )
   IF ( nerror > 0 ) THEN
     CALL prism_abort_proto( ncomp_id, 'oas_cos_define', 'Failure in allocating exfld' )
     RETURN
   ENDIF

   ! -----------------------------------------------------------------
   ! ... Define the partition
   ! -----------------------------------------------------------------

   ! whole domain decomposition can be represented by a rectangle
   paral(1) = 2  ! 2 means box partitioning

   ! Global extent in x
   paral(5) = jih_tot
   ! Upper left corner global offset
   paral(2) = isubpos(my_cart_id,1)-1-nboundlines + ( isubpos(my_cart_id,2)-1-nboundlines) * paral(5)
   ! Local extent in x
   paral(3) = jih
   ! Local extent in y
   paral(4) = jjh

   ! -----------------------------------------------------------------
   ! ... Define the shape of the valid region without the halo and overlaps between cpus
   ! -----------------------------------------------------------------

   ishape(1,1) = 1
   ishape(2,1) = jih
   ishape(1,2) = 1
   ishape(2,2) = jjh

   ! -----------------------------------------------------------------
   ! ... Define the partition 
   ! -----------------------------------------------------------------

   ! Compute global offsets and local extents
 
   IF ( IOASISDEBUGLVL > 1 ) THEN
     WRITE(NULOUT,*) 'local ', paral(3:4), my_cart_id
     WRITE(NULOUT,*) 'global ', paral(2), paral(5), my_cart_id
     WRITE(NULOUT,*) 'isubpos ', my_cart_id,isubpos(my_cart_id,1),isubpos(my_cart_id,2)
   CALL FLUSH(NULOUT)
   ENDIF
  
   CALL prism_def_partition_proto( id_part, paral, nerror )
   IF( nerror /= PRISM_Success )   CALL prism_abort_proto (ncomp_id, 'oas_cos_define',   &
        &                                           'Failure in prism_def_partition' )

   IF ( IOASISDEBUGLVL > 0 ) THEN
     WRITE(nulout,*) ' oascos: oas_cos_define: prism_def_partition '
   ENDIF
   ! -----------------------------------------------------------------
   ! ... Variable definition
   ! ----------------------------------------------------------------

   ! Default values
   ssnd(:)%laction=.FALSE.  ; srcv(:)%laction=.FALSE.
   ssnd(:)%clgrid=clgrd     ; srcv(:)%clgrid=clgrd

   ssnd(1)%clname='COSTEMPE'
   ssnd(2)%clname='COSUWIND'
   ssnd(3)%clname='COSVWIND'
   ssnd(4)%clname='COSSPWAT'   ! specific water vapor content
   ssnd(5)%clname='COSTHICK'   ! thickness of lowest level (m)
   ssnd(6)%clname='COSPRESS'   ! surface pressure (Pa)
   ssnd(7)%clname='COSDIRSW'   ! direct shortwave downward radiation (W/m2)
   ssnd(8)%clname='COSDIFSW'   ! diffuse shortwave downward radiation (W/m2)
   ssnd(9)%clname='COSLONGW'   ! longwave downward radiation (W/m2)
   ssnd(10)%clname='COSCVRAI'  ! convective rain precipitation      (kg/m2*s)
   ssnd(11)%clname='COSCVSNW'  ! convective snow precipitation      (kg/m2*s)
   ssnd(12)%clname='COSGSRAI'  ! gridscale rain precipitation
   ssnd(13)%clname='COSGSSNW'  ! gridscale snow precipitation
   ssnd(14)%clname='COSGRAUP'  ! gridscale graupel precipitation
   ssnd(15)%clname='COSCVPRE'  ! total convective precipitation
   ssnd(16)%clname='COSGSPRE'  ! total gridscale precipitation
!MU (13.09.2012)
   ssnd(17)%clname='COSCO2PP'  ! CO2 partial pressure
!MU (13.09.2012)

!    ksnd = 16                 !CPS now defined in oas_cos_vardef
 
   srcv(1)%clname='COS_TAUX'      !  zonal wind stress
   srcv(2)%clname='COS_TAUY'      !  meridional wind stress
   srcv(3)%clname='COSLATEN'      !  total latent heat flux (W/m**2)
   srcv(4)%clname='COSSENSI'      !  total sensible heat flux (W/m**2)
   srcv(5)%clname='COSINFRA'      ! emitted infrared (longwave) radiation (W/m**2)
   srcv(6)%clname='COSALBED'      ! direct albedo
   srcv(7)%clname='COSALBEI'      ! diffuse albedo
!MU (13.09.12)
   srcv(8)%clname='COSCO2FL'      ! net CO2 flux (now only photosynthesis rate)
!MU (13.09.12)
   srcv(9)%clname='COS_RAM1'      ! Aerodynamic Resistance (s/m)
   srcv(10)%clname='COS_RAH1'      ! Aerodynamic Resistance (s/m)
   srcv(11)%clname='COS_RAW1'      ! Aerodynamic Resistance (s/m)
   srcv(12)%clname='COS_TSF1'      ! Surface Temperature (K)   
   srcv(13)%clname='COS_QSF1'      ! Surface Humidity (kg/kg)
!MU (12.04.13)
   srcv(14)%clname='COSPHOTO'      ! photosynthesis rate (umol CO2 m-2s-1)
   srcv(15)%clname='COSPLRES'      ! plant respiration (umol CO2 m-2s-1)
!MU (12.04.13)

! Forcing Variable selection
   ssnd(1:9)%laction=.TRUE.  
   ssnd(15:16)%laction=.TRUE.   !Coupling only total convective and gridscale precipitations
!MU (17.01.2013)
   ! Coupling CO2 partial pressure
     !If CO2 is not initialized in COSMO a dummy is sent to CLM and CO2 content from CLM is used
   ssnd(17)%laction=.TRUE. !always true for this version
!C   IF (ssnd(17)%laction==.FALSE.) THEN
!     PRINT*, 'ERROR (oas_cos_define): ssnd(17)%laction has to be .TRUE.'
!     PRINT*, '----- If CO2 is not initialized in COSMO a dummy is sent to CLM and CO2 content from CLM is used.'
!     STOP
!   ENDIF
!MU (17.01.2013)

 IF (cpl_scheme) THEN    !CPS
! Lower Boundary variable selection
   srcv(5:7)%laction=.TRUE.
   srcv(8)%laction=.TRUE.
   srcv(14)%laction=.FALSE.
   srcv(15)%laction=.FALSE.
   srcv(9:13)%laction=.TRUE.!CPS
 ELSE
   srcv(1:7)%laction=.TRUE.
   srcv(8)%laction=.TRUE.
   srcv(14)%laction=.FALSE.
   srcv(15)%laction=.FALSE.
 END IF


   var_nodims(1) = 2           ! Dimension number of exchanged arrays
   var_nodims(2) = 1           ! number of bundles (always 1 for OASIS3)

   ! ... Announce send variables. 
   !
   DO ji = 1, nmaxfld
     IF ( ssnd(ji)%laction ) THEN 
            
        CALL prism_def_var_proto( ssnd(ji)%nid, ssnd(ji)%clname, id_part, var_nodims, &
                                      PRISM_Out, ishape, PRISM_REAL, nerror )
        IF ( nerror /= PRISM_Success )   CALL prism_abort_proto( ssnd(ji)%nid, 'oas_cos_define',   &
               &                             'Failure in prism_def_var for '//TRIM(ssnd(ji)%clname))
      ENDIF
   END DO
   !
   ! ... Announce received variables. 
   !
!CPS    krcv = 0
   DO ji = 1, nmaxfld
     IF ( srcv(ji)%laction ) THEN 

        CALL prism_def_var_proto( srcv(ji)%nid, srcv(ji)%clname, id_part, var_nodims, &
                                      PRISM_In, ishape, PRISM_REAL, nerror )
        IF ( nerror /= PRISM_Success )   CALL prism_abort_proto( srcv(ji)%nid, 'oas_cos_define',   &
               &                                               'Failure in prism_def_var for '//TRIM(srcv(ji)%clname))
!CPS            krcv = krcv + 1
     ENDIF
   END DO
   !
   ! Allocate array to store received fields between two coupling steps
   ALLOCATE( frcv(ie, je, krcv), stat = nerror )
   IF ( nerror > 0 ) THEN
     CALL prism_abort_proto( ncomp_id, 'oas_cos_define', 'Failure in allocating exfld' )
     RETURN
   ENDIF
      
   !------------------------------------------------------------------
   ! End of definition phase (must be call by all processes , including PE not involved in the coupling
   !------------------------------------------------------------------
   !
      
   CALL prism_enddef_proto( nerror )
   IF ( nerror /= PRISM_Success )   CALL prism_abort_proto ( ncomp_id, 'oas_cos_define', 'Failure in prism_enddef')

 ENDIF        ! lpe_cpl 
     
 IF ( IOASISDEBUGLVL > 0 ) THEN 
    WRITE(nulout,*) "oascos: prism_enddef"
 ENDIF
  
 ! Get the coupling frequency of COSMO with CLM, only TEMP variable is choosen because cplfreq is same for all variables 
 CALL prism_get_freq(ssnd(1)%nid, cplfreq, nerror)
 IF ( IOASISDEBUGLVL >= 0 ) THEN
    WRITE(nulout,*) "oascos: coupling frequency :", cplfreq
 ENDIF

 CALL MPI_Barrier(kl_comm, nerror)

!------------------------------------------------------------------------------
!- End of the Subroutine
!------------------------------------------------------------------------------

END SUBROUTINE oas_cos_define

!------------------------------------------------------------------------------
!- External sub-subroutine
!------------------------------------------------------------------------------

SUBROUTINE constant_corners(longc,latc,startlonc,startlatc)

!---------------------------------------------------------------------------------
! Obtain geographical corners of the COSMO grid using the rotated grid as the input
!--------------------------------------------------------------------------------

USE data_parameters,  ONLY : wp, iintegers
USE data_parallel,    ONLY :   &
              my_cart_id,      & ! rank of this subdomain in the cartesian communicator
              isubpos,         & ! positions of the subdomains in the total domain. Given
                                 ! are the i- and the j-indices of the lower left and the
                                 ! upper right grid point in the order
                                 !                  i_ll, j_ll, i_ur, j_ur.
                                 ! Only the interior of the domains are considered, not
                                               ! the boundary lines.
              nboundlines        ! number of boundary lines of the domain for which
                                 ! no forecast is computed = overlapping boundary
                                 ! lines of the subdomains
USE data_modelconfig, ONLY : ie, je, dlon, dlat, pollon, pollat, polgam, degrad
USE utilities,        ONLY :    &
      phirot2phi,               &!
      rlarot2rla                 !


IMPLICIT NONE

REAL (KIND=wp),   INTENT (IN)    ::  startlonc,startlatc
REAL (KIND=wp),   INTENT (OUT)   ::  longc(ie,je),latc(ie,je)

! Local variable:
REAL (KIND=wp)        ::  zlats, zlons

INTEGER (KIND=iintegers)  ::  i, j, i_td, j_td


 ! constant fields related to the grid (rlat, rlon)
 ! ---------------------------------------------------------------------

 j_td = isubpos(my_cart_id,2) - nboundlines - 1

 DO j = 1 , je
   j_td = j_td + 1
   ! cos (lat) and 1 / cos (lat)
   zlats        = startlatc + (j_td-1) * dlat

   i_td = isubpos(my_cart_id,1) - nboundlines - 1

   DO i = 1 , ie
       i_td = i_td + 1
       ! geographical latitude and longitude
       zlons  = startlonc + (i_td-1) * dlon

       IF (zlons  > 180.0) THEN
          zlons  = zlons  - 360.0_wp
       ENDIF

       latc(i,j) = phirot2phi ( zlats , zlons , pollat, pollon, polgam) * degrad
       longc(i,j) = rlarot2rla ( zlats , zlons , pollat, pollon, polgam) * degrad

   ENDDO
 ENDDO

END SUBROUTINE constant_corners

