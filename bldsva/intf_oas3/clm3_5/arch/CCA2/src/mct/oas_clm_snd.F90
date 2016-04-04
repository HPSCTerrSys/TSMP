SUBROUTINE oas_clm_snd( kid, kstep, pdata,begg,endg, kinfo )

!---------------------------------------------------------------------
! Description:
!  This routine sends CLM3.5 fields to OASIS3 coupler at each coupling
!  time step defined in namcouple
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

USE shr_kind_mod ,       ONLY : r8 => shr_kind_r8
USE oas_clm_vardef
USE spmdMod     , ONLY : iam, masterproc
USE domainMod   , ONLY : adomain
USE decompMod    , only : adecomp

!==============================================================================

IMPLICIT NONE

!==============================================================================

! * Arguments
!
INTEGER,                 INTENT(IN)    :: kid    ! variable index in the array
INTEGER,                 INTENT(OUT)   :: kinfo  ! OASIS4 info argument
INTEGER,                 INTENT(IN)    :: kstep  ! CLM time-step in seconds
INTEGER, INTENT(IN)                               :: begg,endg
REAL(KIND=r8), DIMENSION(begg:endg), INTENT(IN)   :: pdata
!

real(kind=r8), dimension(:), allocatable :: buffer_array
integer :: c,cl,c1,ani,anj,an,owner,last_owner,ai,aj,ier


!------------------------------------------------------------------------------
!- End of header
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!- Begin Subroutine oas_clm_snd 
!------------------------------------------------------------------------------


! do ai = 0 , 127
!   if(ai==iam)then
     
    !   write(6,*)"###FG: iam: ",iam," begg: ",begg," endg: ",endg," total_part_len: ",total_part_len
     
!   endif
!   CALL MPI_Barrier(kl_comm, nerror)
! enddo


allocate(buffer_array(begg:(begg+total_part_len-1)))


ani = adomain%ni
anj = adomain%nj
buffer_array = -999999._r8
last_owner=-1
cl=0
c=0
c1=0
do aj = 1,anj
   do ai = 1,ani
     an = (aj-1)*ani + ai
     owner = adecomp%glo2owner(an)
     if(owner == iam) then
        !if(masterproc)write(6,*)"###FG: ",(begg+c+cl),(begg+c1),an
        buffer_array(begg+c+cl)= pdata(begg+c1)
        c=c+1
        c1=c1+1
     else 
        if(owner == -1)then 
           if(last_owner==iam) then 
              c=c+1
             ! if(masterproc)write(6,*)"###FG: -99999",an
           endif
           if(last_owner==-1) cl=cl+1  
         else
           if(last_owner==-1)   cl=0
        endif
     endif
     if(owner/=-1) last_owner=owner
   enddo
enddo


CALL MPI_Barrier(kl_comm, ier)

 !
 ! Call OASIS at each time step but field sent to other model only at coupling time step
 ! (accumulation otherwise, if asked in the SMIOC configuration file)
 !
 CALL prism_put_proto( ssnd(kid)%nid, kstep, buffer_array, kinfo )

deallocate(buffer_array)

!------------------------------------------------------------------------------
!- End of the Subroutine
!------------------------------------------------------------------------------

END SUBROUTINE oas_clm_snd
