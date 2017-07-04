SUBROUTINE oas_cos_rcv( kid, kstep, pdata, kinfo )

!---------------------------------------------------------------------
! Description:
!  This routine call fields from OASIS3 coupler at each coupling time step.
!
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
! 1.2.1        2015/08/24 Prabhakar Shrestha
!   Added update of outerbound halos using inner domain data
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

USE data_parameters,     ONLY:  wp, iintegers
USE data_modelconfig,    ONLY:  ie, je, jstartpar, jendpar
USE oas_cos_vardef
USE mod_prism_get_proto
!CPS
USE environment,        ONLY : exchg_boundaries
USE data_parallel,      ONLY :  &
    nprocx,          & ! number of processors in x-direction
    nprocy,          & ! number of processors in y-direction
    my_cart_id,      & ! rank of this subdomain in the cartesian communicator
    my_cart_pos,     & ! position of this subdomain in the cartesian grid
                       ! in x- and y-direction
    my_cart_neigh,   & ! neighbors of this subdomain in the cartesian grid
    nboundlines,     & !
    ncomm_type,      & !
    sendbuf,         & !
    isendbuflen,     & !
    imp_reals,       & !
    icomm_cart,      & !
    num_compute        !
USE data_runcontrol , ONLY :   &
    lperi_x,      & ! lartif_data=.TRUE.:  periodic boundary conditions
                    !            =.FALSE.: with Davies conditions
    lperi_y,      & ! lartif_data=.TRUE.:  periodic boundary conditions
                    !            =.FALSE.: with Davies conditions
    l2dim,        & !
    nbl_exchg       !
!CPS
!==============================================================================

IMPLICIT NONE

!==============================================================================

! Arguments

INTEGER(KIND=iintegers), INTENT(IN)   :: kid    ! variable intex in the array
INTEGER(KIND=iintegers), INTENT(IN)   :: kstep  ! ocean time-step in seconds
REAL(KIND=wp), DIMENSION(ie,je), INTENT(OUT)   :: pdata
INTEGER(KIND=iintegers), INTENT(OUT)  :: kinfo  ! OASIS4 info argument
!
LOGICAL                   :: llaction
INTEGER(KIND=iintegers)   :: NULOUT=6
!CPS
INTEGER(KIND=iintegers)   :: ii, jj, izerror, kzdims(24)
CHARACTER (LEN=250)       :: yzerrmsg
!------------------------------------------------------------------------------
!- End of header
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!- Begin Subroutine oas_cos_rcv 
!------------------------------------------------------------------------------

 !  Masked point are not modified by OASIS4
 !  = buffer set to 0 before calling oasis
 exfld=0._wp
 pdata=0._wp   !CPS Infinity values if not initialized
 !
 IF (lpe_cpl) THEN
   !
   ! receive local data from OASIS3 on every process involved in the coupling
   !
   CALL prism_get_proto( srcv(kid)%nid, kstep, exfld(nldi:nlei, nldj:nlej), kinfo )         

   llaction = .false.
   IF( kinfo == PRISM_Recvd   .OR. kinfo == PRISM_FromRest .OR.   &
      kinfo == PRISM_RecvOut .OR. kinfo == PRISM_FromRestOut )   llaction = .TRUE.

   IF ( IOASISDEBUGLVL == 1 )  WRITE(NULOUT,*) "oascos: oas_cos_rcv:  &
     llaction, kinfo, kstep, clname: " , llaction, kinfo, kstep, srcv(kid)%clname

   ! If coupling time step
   IF ( llaction ) THEN

     ! Declare to calling routine that OASIS provided coupling field
     kinfo = OASIS_Rcv

     ! Update array which contains coupling field (only on valid shape)
     pdata(nldi:nlei, nldj:nlej) = exfld(nldi:nlei, nldj:nlej)

     !CPS BANDAID FOR OUTERBOUND HALOs
     !Removes Infinity with inversion of tcm/tch 1/0.
     IF (fhalo) THEN
     !West
     IF (my_cart_neigh(1) == -1 .OR. (lperi_x .AND. my_cart_pos(1) == 0) ) THEN
       DO ii = 1, nldi-1 
         pdata(ii,:) = pdata(nldi,:)
       ENDDO
     ENDIF
     !South
     IF (my_cart_neigh(4) == -1 .OR.(lperi_y .AND. my_cart_pos(2) == 0) ) THEN
       DO jj = 1, nldj-1
         pdata(:,jj) = pdata(:,nldj)
       ENDDO
     ENDIF
     !East
     IF (my_cart_neigh(3) == -1 .OR. (lperi_x .AND. my_cart_pos(1) == nprocx-1) ) THEN
       DO ii = nlei+1,ie
         pdata(ii,:) = pdata(nlei,:)
       ENDDO
     ENDIF
     !North
     IF (my_cart_neigh(2) == -1 .OR. (lperi_y .AND. my_cart_pos(2) == nprocy-1) ) THEN
       DO jj = nlej+1,je
         pdata(:,jj) = pdata(:,nlej)
       ENDDO
     ENDIF
     
     ENDIF !fhalo
!     ! exchange CLM data (for safety and for periodic BCs):
     kzdims(1:24)=(/1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
     CALL exchg_boundaries                                                  &
           ( 0,  sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute, ie, je, &
           kzdims, jstartpar, jendpar, nbl_exchg, nboundlines, my_cart_neigh,      &
           lperi_x, lperi_y, l2dim, &
           10000, .FALSE., ncomm_type, izerror, yzerrmsg,                  &
           pdata )

   IF ( IOASISDEBUGLVL == 1 ) THEN
     WRITE(NULOUT,*) '****************'
     WRITE(NULOUT,*) 'oascos: prism_get: Incoming ', srcv(kid)%clname
     WRITE(NULOUT,*) 'oascos: prism_get: ivarid '  , srcv(kid)%nid
     WRITE(NULOUT,*) 'oascos: prism_get:   kstep', kstep
     WRITE(NULOUT,*) 'oascos: prism_get:   info ', kinfo
     WRITE(NULOUT,*) '     - Minimum value is ', MINVAL(pdata(nldi:nlei, nldj:nlej))
     WRITE(NULOUT,*) '     - Maximum value is ', MAXVAL(pdata(nldi:nlei, nldj:nlej))
     WRITE(NULOUT,*) '     -     Sum value is ', SUM(pdata(nldi:nlei, nldj:nlej))
     WRITE(NULOUT,*) '****************'
   ENDIF
         
   ELSE
     ! Declare to calling routine that OASIS did not provide coupling field
     kinfo = OASIS_idle     
   ENDIF

 ENDIF             ! lpe_cpl

!------------------------------------------------------------------------------
!- End of the Subroutine
!------------------------------------------------------------------------------

END SUBROUTINE oas_cos_rcv
