SUBROUTINE oas_clm_rcv( kid, kstep, pdata, begg,endg, kinfo )

!---------------------------------------------------------------------
! Description:
!  This routine receiveds atmospheric field from OASIS3 coupler at
!  each timestep
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

USE shr_kind_mod , ONLY : r8 => shr_kind_r8
USE spmdMod ,ONLY : iam,npes,mpicom
USE oas_clm_vardef

 USE decompMod ,               ONLY :  adecomp
USE domainMod   , ONLY : adomain
!==============================================================================

IMPLICIT NONE

!==============================================================================

! Arguments
INTEGER,                          INTENT(IN)        :: kid    ! variable intex in the array
INTEGER,                          INTENT(IN)        :: kstep  ! time-step in seconds
INTEGER, INTENT(IN)                               :: begg,endg
REAL(KIND=r8), DIMENSION(begg:endg), INTENT(OUT)  :: pdata
INTEGER,                          INTENT(OUT)       :: kinfo  ! OASIS info argument

! Local Variables
LOGICAL                                           :: llaction
INTEGER                                           :: NULOUT=6,ier
real(kind=r8), dimension(:) , allocatable :: buffer_array
integer :: c,cl,c1,ani,anj,an,owner,last_owner,ai,aj

!------------------------------------------------------------------------------
!- End of header
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!- Begin Subroutine oas_clm_rcv 
!------------------------------------------------------------------------------


pdata=0.


allocate(buffer_array(begg:(begg+total_part_len-1)))


CALL prism_get_proto( srcv(kid)%nid, kstep, buffer_array , kinfo )


ani = adomain%ni
anj = adomain%nj
last_owner=-1
c=0
c1=0
cl=0
do aj = 1,anj
   do ai = 1,ani
     an = (aj-1)*ani + ai
     owner = adecomp%glo2owner(an)
     if(owner == iam) then
        pdata(begg+c1)=buffer_array(begg+c+cl)
        c=c+1
        c1=c1+1 
     else
        if(owner == -1) then
           if(last_owner==iam) c=c+1
           if(last_owner==-1)   cl=cl+1
        else
           if(last_owner==-1)   cl=0
        endif
     endif
     if(owner/=-1) last_owner=owner
   enddo
enddo

deallocate(buffer_array)

 llaction = .false.
 IF( kinfo == PRISM_Recvd   .OR. kinfo == PRISM_FromRest .OR.   &
     kinfo == PRISM_RecvOut .OR. kinfo == PRISM_FromRestOut )   llaction = .TRUE.


 IF ( llaction ) THEN

  ! Declare to calling routine that OASIS provided coupling field
  kinfo = OASIS_Rcv
         
 ELSE
  ! Declare to calling routine that OASIS did not provide coupling field
  kinfo = OASIS_idle     
 ENDIF

!------------------------------------------------------------------------------
!- End of the Subroutine
!------------------------------------------------------------------------------


END SUBROUTINE oas_clm_rcv
