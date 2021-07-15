!-------------------------------------------------------------------------------------------
!Copyright (c) 2013-2016 by Wolfgang Kurtz and Guowei He (Forschungszentrum Juelich GmbH)
!
!This file is part of TerrSysMP-PDAF
!
!TerrSysMP-PDAF is free software: you can redistribute it and/or modify
!it under the terms of the GNU Lesser General Public License as published by
!the Free Software Foundation, either version 3 of the License, or
!(at your option) any later version.
!
!TerrSysMP-PDAF is distributed in the hope that it will be useful,
!but WITHOUT ANY WARRANTY; without even the implied warranty of
!MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!GNU LesserGeneral Public License for more details.
!
!You should have received a copy of the GNU Lesser General Public License
!along with TerrSysMP-PDAF.  If not, see <http://www.gnu.org/licenses/>.
!-------------------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------------------
!obs_op_pdaf.F90: TerrSysMP-PDAF implementation of routine
!                 'obs_op_pdaf' (PDAF online coupling)
!-------------------------------------------------------------------------------------------

!$Id: obs_op_pdaf.F90 1441 2013-10-04 10:33:42Z lnerger $
!BOP
!
! !ROUTINE: obs_op_pdaf --- Implementation of observation operator
!
! !INTERFACE:
SUBROUTINE obs_op_pdaf(step, dim_p, dim_obs_p, state_p, m_state_p)

! !DESCRIPTION:
! User-supplied routine for PDAF.
! Used in the filters: SEEK/SEIK/EnKF/ETKF/ESTKF
!
! The routine is called during the analysis step.
! It has to perform the operation of the
! observation operator acting on a state vector.
! For domain decomposition, the action is on the
! PE-local sub-domain of the state and has to
! provide the observed sub-state for the PE-local
! domain.
!
! !REVISION HISTORY:
! 2013-02 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
   USE mod_assimilation, &
        ONLY: obs_index_p, dim_obs,obs_filename,obs_index_p_TB,&
              ens_TB,member_TB, sc_p, &
        depth_obs_p !, obs_p
   USE mod_read_obs, ONLY: crns_flag !clm_obs
   use mod_parallel_model, &
        only: model,mype_model,npes_model,mype_world,npes_world,COMM_model
   use mod_tsmp, &
        only: tag_model_clm,tag_model_parflow, &
              nprocpf,nprocclm,lcmem,nx_local,ny_local,nz_local, & !SPo add lcmem
              soilay, soilay_fortran, nz_glob 
   USE, INTRINSIC :: iso_c_binding
!   USE mod_parallel_model, ONLY: tcycle 

   ! LSN: module load for the implementation of CMEM model
   USE mod_parallel_pdaf, &     ! Parallelization variables fro assimilation
        ONLY:n_modeltasks,task_id,COMM_filter,mype_filter,npes_filter,&
            MPI_DOUBLE_PRECISION,MPI_INTEGER,filterpe,&
            MPI_SUCCESS,mype_couple         
   USE spmdMod      , only: masterproc,iam,mpicom,npes

#if defined CLMSA
   USE enkf_clm_mod, & 
        ONLY : clm_varsize, clm_paramarr, clmupdate_swc, clmupdate_T
#endif

   IMPLICIT NONE
! !ARGUMENTS:
  INTEGER, INTENT(in)   :: step               ! Currrent time step
  INTEGER, INTENT(in)   :: dim_p              ! PE-local dimension of state
  INTEGER, INTENT(in)   :: dim_obs_p          ! Dimension of observed state
  REAL, INTENT(in)      :: state_p(dim_p)     ! PE-local model state
  REAL, INTENT(out)     :: m_state_p(dim_obs_p) ! PE-local observed state

  character*200         :: inparam_fname 
  REAL,DIMENSION(1)     :: TB(dim_obs)
  REAL,ALLOCATABLE      :: ens_TB_tmp(:,:)
  INTEGER               :: i,j,k,nerror,nproc
  CHARACTER (len = 110) :: OBSFILE
  !REAL,ALLOCATABLE      :: ens_TB(:)
  !common /coeff/ ens_TB
 
! CALLING SEQUENCE:
! Called by: PDAF_seek_analysis   (as U_obs_op)
! Called by: PDAF_seik_analysis, PDAF_seik_analysis_newT
! Called by: PDAF_enkf_analysis_rlm, PDAF_enkf_analysis_rsm
! EOP

! *** local variables ***

! hcp test with hardcoding variable declaration
real(8), dimension(:), allocatable :: soide !soil depth
!real(8), dimension(0:12), parameter :: &
! soide=(/0.d0,  0.02d0,  0.05d0,  0.1d0,  0.17d0, 0.3d0,  0.5d0, &
!                0.8d0,   1.3d0,   2.d0,  3.d0, 5.d0,  12.d0/) !soil depth

real(8) :: tot, avesm
integer :: nsc
! end of hcp 


! *********************************************
! *** Perform application of measurement    ***
! *** operator H on vector or matrix column ***
! *********************************************
 !SPo add switch for cmem 
 IF (lcmem) THEN 

   CALL MPI_Barrier(COMM_filter, nerror)
   ! if (model == tag_model_clm) then

      if(.not.allocated(ens_TB_tmp)) ALLOCATE(ens_TB_tmp(dim_obs,n_modeltasks))
      if(.not.allocated(ens_TB)) ALLOCATE(ens_TB(dim_obs,n_modeltasks))
      if(.not.allocated(member_TB)) then
         ALLOCATE(member_TB(nprocpf+nprocclm))
         do i = 1, nprocpf+nprocclm
            member_TB(i) = 1
         end do
      end if 
     
      ! LSN: Broadcasts brightness temperature ensemble (ens_TB) from 
      ! the master process to all other processes of COMM_filter   
      CALL MPI_Bcast(ens_TB, dim_obs*n_modeltasks, MPI_DOUBLE_PRECISION, nprocpf, COMM_filter, nerror)
      IF (mype_filter == nprocpf.AND.member_TB(mype_filter+1)==1) THEN 
         write(*,*) 'obs_op_pdaf: ens_TB is', ens_TB
      END IF 

      DO i = 1, dim_obs
         DO j = 1, n_modeltasks
            ens_TB_tmp(i,j) = ens_TB(i,j)
         END DO
      END DO      
    
      DO i = 1, dim_obs
         TB(i) = ens_TB_tmp(i,member_TB(mype_filter+1))
      END DO   
      
      ! LSN: member_TB is a set of counter for member loop
      IF (member_TB(mype_filter+1) == n_modeltasks) THEN
         member_TB(mype_filter+1) = 1
      ELSE 
         member_TB(mype_filter+1) = member_TB(mype_filter+1)+1
      END IF

      ! LSN: assign ens_TB instead of state_p to m_state_p
      DO i = 1, dim_obs_p
         m_state_p(i) = TB(obs_index_p_TB(i))
         ! WRITE(*,*) 'obs_op_pdaf: member_TB is', member_TB(mype_filter+1)-1
         WRITE(*,*) 'obs_op_pdaf: obs_index_p_TB(i) is ',obs_index_p_TB(i)
         WRITE(*,*) 'obs_op_pdaf: TB(obs_index_p_TB(i)) is ',TB(obs_index_p_TB(i))
      END DO
   
  
      IF (ALLOCATED(ens_TB_tmp)) DEALLOCATE(ens_TB_tmp)
      !LSN: the end of the implementation of CMEM model here
   CALL MPI_Barrier(COMM_filter, nerror)
   ! end if  
 ELSE
   ! point observations
   ! LSN: the original observation operator of soil moisture in Parflow
   DO i = 1, dim_obs_p
      m_state_p(i) = state_p(obs_index_p(i))
   END DO
 END IF

#if defined CLMSA
if (clmupdate_swc.NE.0) then
  DO i = 1, dim_obs_p
     m_state_p(i) = state_p(obs_index_p(i))
  END DO
endif
if (clmupdate_T.EQ.1) then
  DO i = 1, dim_obs_p
     m_state_p(i) &
    = (exp(-0.5*clm_paramarr(obs_index_p(i))) &
                     *state_p(obs_index_p(i))**4 & 
       +(1.-exp(-0.5*clm_paramarr(obs_index_p(i)))) &
                     *state_p(clm_varsize+obs_index_p(i))**4)**0.25
  END DO
!  write(*,*) 'now is cycle ', tcycle
!  write(*,*) 'obs_p =', obs_p(:)
!  write(*,*) 'model LST', m_state_p(:)
!  write(*,*) 'TG', state_p(obs_index_p(:))
!  write(*,*) 'TV', state_p(clm_varsize+obs_index_p(:))
endif
#else
 if (crns_flag.EQ.1) then
     call C_F_POINTER(soilay,soilay_fortran,[nz_glob])
     Allocate(soide(0:nz_glob))
     soide(0)=0.d0
     do i=1,nz_glob
       soide(i)=soide(i-1)+soilay_fortran(nz_glob-i+1) 
     enddo
     do i = 1, dim_obs_p
       nsc= size(sc_p(i)%scol_obs_in(:))
       avesm=0.d0
       do j=1, nsc-1
           avesm=avesm+(1.d0-0.5d0*(soide(j)+soide(j-1))/depth_obs_p(i))*(soide(j)-soide(j-1)) &
                 *state_p(sc_p(i)%scol_obs_in(j))/depth_obs_p(i)
       enddo
       avesm=avesm+(1.d0-0.5d0*(depth_obs_p(i)+soide(nsc-1))/depth_obs_p(i))*(depth_obs_p(i)-soide(nsc-1)) &
             *state_p(sc_p(i)%scol_obs_in(nsc))/depth_obs_p(i)
       tot=0.d0
       do j=1, nsc-1
           tot=tot+(1.d0-0.5d0*(soide(j)+soide(j-1))/depth_obs_p(i))*(soide(j)-soide(j-1)) &
               /depth_obs_p(i)
       enddo
       tot=tot+(1.d0-0.5d0*(depth_obs_p(i)+soide(nsc-1))/depth_obs_p(i))*(depth_obs_p(i)-soide(nsc-1)) &
          /depth_obs_p(i)

       avesm=avesm/tot
       m_state_p(i)=avesm
     enddo
     deallocate(soide)
 else 
  DO i = 1, dim_obs_p
     m_state_p(i) = state_p(obs_index_p(i))
  END DO
 endif
#endif

END SUBROUTINE obs_op_pdaf
