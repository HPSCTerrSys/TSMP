SUBROUTINE receive_fld_2cos(kt, lcoupled)

!---------------------------------------------------------------------
! Description:
!  This routine receives the atmospheric fields from COSMO and updates
!  the CLM3.53.53.5 variable x(:,:,:)
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
! 2.1        2012/09/26 Markus Uebel, Prabhakar Shrestha 
!   CO2 coupling included
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
USE shr_kind_mod ,            ONLY : r8 => shr_kind_r8
USE clm_time_manager ,        ONLY : dtime,                            &! timestep in second
                                     nelapse,                          &!
                                     get_nstep                          ! return timestep number

USE atmdrvMod,                ONLY : x                                  ! input fields
USE spmdMod ,                 ONLY : masterproc, dummy => MPI_LOGICAL   ! CPS
USE netcdf

!==============================================================================

IMPLICIT NONE

!==============================================================================

! Parameter List
INTEGER, INTENT(IN)                    ::   kt                           ! time step

LOGICAL, INTENT(OUT)                   ::  lcoupled                      ! CPS flag for coupling 
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
!MU (13.09.2012)
INTEGER, PARAMETER ::   jps_co2 = 17            ! CO2 partial pressure (Pa)
!MU (13.09.2012)

! Local Variables

INTEGER                                :: jn, isec, ier, ctr   ! CPS added ctr
INTEGER , DIMENSION(krcv)              ::   nrcvinfo           ! OASIS info argument
INTEGER                                :: ji, jj

REAL(KIND=r8) , DIMENSION(ndlon,ndlat) ::   ztmp1   

CHARACTER(8)                           :: dateoas
CHARACTER(10)                          :: timeoas
CHARACTER(5)                           :: zoneoas
INTEGER ,DIMENSION(8)                  :: valuesoas
REAL(KIND=r8)                          :: millisec
REAL(KIND=r8)                          :: sec_date
INTEGER ,DIMENSION(3)                  :: dimids                ! CPS 2 to 3
INTEGER ,DIMENSION(krcv)               :: il_var_id
INTEGER                                :: il_file_id, status
REAL(KIND=r8), SAVE                    :: tic, toc

!------------------------------------------------------------------------------
!- End of header
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!- Begin Subroutine receive_fld2cos 
!------------------------------------------------------------------------------

 isec = dtime * kt

 IF ( masterproc ) THEN                   !masterporc

 ! Default value (changed when coupling)
 lcoupled = .FALSE.
 nrcvinfo = OASIS_idle

 ! Performances measurement
 IF ( IOASISDEBUGLVL == 2 ) THEN
   tic=MPI_Wtime()
 ENDIF

 DO jn = 1, krcv
   IF( srcv(jn)%laction ) THEN
     CALL oas_clm_rcv( jn, isec, ztmp1(:,:), nrcvinfo(jn) )
     ! change masked point value to avoid warning message
       IF ( jn == jps_t )               &
         WHERE ( .not. llmask(:,:,1) ) ztmp1(:,:) = 60._r8 !CPS??

       ! Update arrays which contains coupling fields only at coupling time step
       IF( nrcvinfo(jn) == OASIS_Rcv ) THEN
          exfld(:,:,jn)=ztmp1(:,:)
          lcoupled = .TRUE.
       ENDIF
   ENDIF
 ENDDO

 ! Performances measurement
 IF ( IOASISDEBUGLVL == 2 ) THEN
   IF( nrcvinfo(1) == OASIS_Rcv ) THEN
     WRITE(6, *)  'oasclm: receive_fld_2cos: calc_time ', tic-toc
     toc=MPI_Wtime()
     WRITE(6, *)  'oasclm: receive_fld_2cos: After rcv ', toc
     WRITE(6, *)  'oasclm: receive_fld_2cos: rcv_time ', toc-tic
     CALL flush(6)
   ENDIF
 ENDIF

 ! Write coupling fields for debug or new mask setting
 IF ( IOASISDEBUGLVL == 1 ) THEN
   IF (isec == 0 ) THEN                            !CPS create file at beginning 
     status =  nf90_create("debugrcv_clm_cos.nc", NF90_CLOBBER, il_file_id)
     status =  nf90_def_dim(il_file_id, "longitude", ndlon, dimids(1))
     status =  nf90_def_dim(il_file_id, "latitude", ndlat, dimids(2))
     status =  nf90_def_dim(il_file_id, "time", nelapse, dimids(3))    !CPS 
     ctr = 0
     DO jn = 1, krcv
       IF (srcv(jn)%laction) THEN
         ctr = ctr+1
         status = nf90_def_var(il_file_id, srcv(jn)%clname, NF90_DOUBLE, dimids, il_var_id(ctr))
       ENDIF
     ENDDO
     status =  nf90_enddef(il_file_id)
     status =  nf90_close(il_file_id)
   !debug file end definition
   ENDIF                                  !CPS isec == 0

   !debug file output start
   status = nf90_open("debugrcv_clm_cos.nc", NF90_WRITE, il_file_id)
   ctr = 0
   DO jn = 1,krcv
     IF (srcv(jn)%laction) THEN
       ctr = ctr+1
       status = nf90_inq_varid(il_file_id, srcv(jn)%clname , il_var_id(ctr))
!CPS       WRITE(6, *) "CPS CLM " ,jn, kt, MINVAL(exfld(:,:,jn)),MAXVAL(exfld(:,:,jn))
       status =  nf90_put_var(il_file_id, il_var_id(ctr), exfld(:,:,jn),       &
                                         start = (/ 1, 1, kt+1/),              &
                                         count = (/ ndlon, ndlat, 1 /) )
     ENDIF
   ENDDO

   status =  nf90_close(il_file_id)

 ENDIF            !IOASISDEBUGLVL=1
 
! Update CLM variable with the received input from COSMO

 IF( srcv(jps_t)%laction )  x(:,:,1)   = exfld(:,:,jps_t)     ! air temperature (K)
 IF( srcv(jps_u)%laction )  x(:,:,15)  = exfld(:,:,jps_u)     
 IF( srcv(jps_v)%laction )  x(:,:,16)  = exfld(:,:,jps_v)     
 IF( srcv(jps_u)%laction .AND. srcv(jps_v)%laction)          &
         x(:,:,2) = sqrt(x(:,:,15)**2 + x(:,:,16)**2)         ! wind (m/s)
 IF( srcv(jps_q)%laction )  x(:,:,3)   = exfld(:,:,jps_q)     ! specific humidity (kg/kg)

 IF( srcv(jps_th)%laction )  x(:,:,6)  = exfld(:,:,jps_th)    ! thickness of lowest level (m)

 IF( srcv(jps_pr)%laction )  x(:,:,7)  = exfld(:,:,jps_pr)    ! surface pressure (Pa)

 IF( srcv(jps_rs)%laction )  x(:,:,9)  = exfld(:,:,jps_rs)    ! direct shortwave downward radiation (W/m2)
 IF( srcv(jps_fs)%laction )  x(:,:,10) = exfld(:,:,jps_fs)    ! diffuse shortwave downward radiation (W/m2)
 IF( srcv(jps_lw)%laction )  x(:,:,11) = exfld(:,:,jps_lw)    ! longwave downward radiation (W/m2)

 ! Coupling separately snow, rain and graupel  
 IF( srcv(jps_cr)%laction .AND. srcv(jps_cs)%laction )       &
         x(:,:,13) = exfld(:,:,jps_cr) + exfld(:,:,jps_cs)    ! convective precip (mm/s)
 IF( srcv(jps_gr)%laction .AND. srcv(jps_gs)%laction ) THEN
         x(:,:,14) = exfld(:,:,jps_gr) + exfld(:,:,jps_gs)    ! gridscale precip (mm/s) 
  IF( srcv(jps_gg)%laction )                                 &
         x(:,:,14) = x(:,:,14) +  exfld(:,:,jps_gg)

 ENDIF     ! if coupled graupel
!MU (13.09.2012)
 IF( srcv(jps_co2)%laction )  x(:,:,17) = exfld(:,:,jps_co2)    ! CO2 partial pressure (Pa)
!MU (13.09.2012)

 ! Coupling only total convective and gridscale precipitations
 IF( srcv(jps_cp)%laction )  x(:,:,13) = exfld(:,:,jps_cp)    ! total convective precip (mm/s)
 IF( srcv(jps_gp)%laction )  x(:,:,14) = exfld(:,:,jps_gp)    ! total gridscale precip (mm/s)

! Uncomment in case of emergency only
!         x(:,:,1) = 288.15
!         x(:,:,15) = 1.
!         x(:,:,16) = 1.
!         x(:,:,2) = sqrt(2.)
!         x(:,:,3) = 0.009
!         x(:,:,6) = 135.
!         x(:,:,7) = 102360.
!         x(:,:,9) = 0.
!         x(:,:,10) = 0.
!         x(:,:,11) = 365.
!         x(:,:,13) = 0.
!         x(:,:,14) = 0.
!MU (26.09.2012)
!         x(:,:,17) = 39.5
!MU (26.09.2012)

 ENDIF       ! masterproc

!CPS WRITE(6, *)  "oasclm: receive_fld_2cos: Receive atmospheric field complete" 

 !  Wait coupling result (involved PE and non involved PE)
 CALL MPI_Bcast( lcoupled, 1, dummy, 0, kl_comm, ier )

!------------------------------------------------------------------------------
!- End of the Subroutine
!------------------------------------------------------------------------------

END SUBROUTINE receive_fld_2cos
