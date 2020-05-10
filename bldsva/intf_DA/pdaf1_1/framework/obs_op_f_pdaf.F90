!-------------------------------------------------------------------------------------------
!Copyright (c) 2013-2016 by Wolfgang Kurtz, Guowei He and Mukund Pondkule (Forschungszentrum Juelich GmbH)
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
!obs_op_f_pdaf.F90: TerrSysMP-PDAF implementation of routine
!                   'obs_op_f_pdaf' (PDAF online coupling)
!-------------------------------------------------------------------------------------------

!$Id: obs_op_f_pdaf.F90 1441 2013-10-04 10:33:42Z lnerger $
!BOP
!
! !ROUTINE: obs_op_f_pdaf --- Implementation of observation operator
!
! !INTERFACE:
SUBROUTINE obs_op_f_pdaf(step, dim_p, dim_obs_f, state_p, m_state_f)

  ! !DESCRIPTION:
  ! User-supplied routine for PDAF.
  ! Used in the filters: LSEIK/LETKF/LESTKF
  !
  ! The routine is called in PDAF\_X\_update
  ! before the loop over all local analysis domains
  ! is entered.  The routine has to perform the 
  ! operation of the observation operator acting on 
  ! a state vector.  The full vector of all 
  ! observations required for the localized analysis
  ! on the PE-local domain has to be initialized.
  ! This is usually data on the PE-local domain plus 
  ! some region surrounding the PE-local domain. 
  ! This data is gathered by MPI operations. The 
  ! gathering has to be done here, since in the loop 
  ! through all local analysis domains, no global
  ! MPI operations can be performed, because the 
  ! number of local analysis domains can vary from 
  ! PE to PE.
  !
  ! !REVISION HISTORY:
  ! 2013-09 - Lars Nerger - Initial code
  ! Later revisions - see svn log
  !
  ! !USES:
  USE mod_assimilation, &
       ONLY: obs_index_p, local_dims_obs, obs_id_p, m_id_f, &
       var_id_obs, dim_obs_p, dim_obs,obs_filename,obs_index_p_TB
  USE mod_parallel_pdaf, &
       ONLY: mype_filter, npes_filter, comm_filter, &
             n_modeltasks,task_id,COMM_filter,mype_filter,npes_filter, &
             MPI_DOUBLE, MPI_DOUBLE_PRECISION, MPI_INT,MPI_SUM
  
  ! LSN: module load for the implementation of CMEM model 
  use mod_parallel_model, &
        only: model,mype_model,npes_model,mype_world,npes_world
  use mod_tsmp, &
        only: tag_model_clm,nprocpf
  USE spmdMod      , only: masterproc,iam,mpicom,npes
  use clm4cmem     , only: SATELLITE,CLM_DATA
  USE rdclm4pdaf   , only: read_CLM_pdaf
  USE get_tb_cmem  , only: cmem_main
  USE rdclm_wrcmem , only: read_satellite_info
  USE YOMCMEMPAR   , only: INPUTNAMLST,LGPRINT

  IMPLICIT NONE

  ! !ARGUMENTS:
  INTEGER, INTENT(in) :: step                 ! Current time step
  INTEGER, INTENT(in) :: dim_p                ! PE-local dimension of state
  INTEGER, INTENT(in) :: dim_obs_f            ! Dimension of observed state
  REAL, INTENT(in)    :: state_p(dim_p)       ! PE-local model state
  REAL, INTENT(inout) :: m_state_f(dim_obs_f) ! PE-local observed state

  ! !CALLING SEQUENCE:
  ! Called by: PDAF_lseik_update   (as U_obs_op)
  ! Called by: PDAF_lestkf_update  (as U_obs_op)
  ! Called by: PDAF_letkf_update   (as U_obs_op)
  !EOP

  ! local variables
  INTEGER :: ierror, max_var_id
  INTEGER :: i,nerror                  ! Counter
  REAL, ALLOCATABLE :: m_state_tmp(:)  ! Temporary process-local state vector
  INTEGER, ALLOCATABLE :: local_dis(:) ! Displacement array for gathering
  INTEGER, ALLOCATABLE :: m_id_p_tmp(:)
  INTEGER, ALLOCATABLE :: m_id_f_tmp(:)
  
  CHARACTER*200:: inparam_fname
  TYPE(SATELLITE),ALLOCATABLE :: SAT
  TYPE(CLM_DATA),ALLOCATABLE :: CLMVARS
  REAL,DIMENSION(1) :: TB(dim_obs)
  CHARACTER (len = 110) :: current_observation_filename
  ! *********************************************
  ! *** Perform application of measurement    ***
  ! *** operator H on vector or matrix column ***
  ! *********************************************

  ! Initialize process-local observed state
  ALLOCATE(m_state_tmp(local_dims_obs(mype_filter+1)))
  allocate(m_id_p_tmp (local_dims_obs(mype_filter+1)))

  if(allocated(m_id_f)) deallocate(m_id_f)
  allocate(m_id_f(dim_obs_f))
  if(allocated(m_id_f_tmp)) deallocate(m_id_f_tmp)
  allocate(m_id_f_tmp(dim_obs_f))

  !DO i = 1, local_dims_obs(mype_filter+1)
  !   m_state_tmp(i) = state_p(obs_index_p(i))
  !   m_id_p_tmp(i)  = obs_id_p(obs_index_p(i))
  !END DO
  
  ! LSN: the implementation of CMEM model here
   IF (model == tag_model_clm) THEN
     
     allocate(CLMVARS)  ! assign CLMVARS array
     inparam_fname = './obs/input'
     INPUTNAMLST   = trim(inparam_fname)
     
     LGPRINT = .False.
       
     write(current_observation_filename, '(a, i5.5)') trim(obs_filename)//'.',step
     write(*,*) 'The current observation file is ', current_observation_filename
     allocate(SAT)  !assign SAT array
     call read_satellite_info(current_observation_filename,SAT)

     call read_CLM_pdaf(LS=CLMVARS,SAT=SAT)
     WRITE(*,*) 'Read CLM memory is done'
     
     IF (masterproc) THEN      
        call cmem_main(LS=CLMVARS,SATinfo=SAT,TB=TB,step=step)    
     END IF
     
  END IF
  
  CALL MPI_Barrier(COMM_filter, nerror) 

  !Broadcasts brightness temperature (TB) from the master process to all
  !other processes of COMM_filter
  CALL MPI_Bcast(TB, dim_obs, MPI_DOUBLE_PRECISION, nprocpf, COMM_filter,nerror)
  
  DO i = 1, dim_obs_p
     m_state_tmp(i)  = TB(obs_index_p_TB(i))
     m_id_p_tmp(i)   = obs_id_p(obs_index_p_TB(i))
  END DO
 
  IF (ALLOCATED(SAT)) DEALLOCATE(SAT)
  IF (ALLOCATED(CLMVARS)) DEALLOCATE(CLMVARS)
  ! LSN: the end of the implementation of CMEM model here
 
  ! Gather full observed state
  ALLOCATE(local_dis(npes_filter))

  !print *,'local_dims_obs(mype_filter+1) ', local_dims_obs(mype_filter+1)
  !print *,'dim_obs_p ', dim_obs_p

  local_dis(1) = 0
  DO i = 2, npes_filter
     local_dis(i) = local_dis(i-1) + local_dims_obs(i-1)
  END DO

  ! gather local observed states of different sizes in a vector 
  CALL mpi_allgatherv(m_state_tmp, local_dims_obs(mype_filter+1), &
       MPI_DOUBLE_PRECISION, m_state_f, local_dims_obs, local_dis, &  
       MPI_DOUBLE_PRECISION, comm_filter, ierror)

  ! gather m_id_p 
  CALL mpi_allgatherv(m_id_p_tmp, local_dims_obs(mype_filter+1), &
       MPI_INT, m_id_f, local_dims_obs, local_dis, &  
       MPI_INT, comm_filter, ierror)

  ! resort m_id_p
  do i=1,dim_obs_f
     m_id_f_tmp(i) = m_id_f(i)     
  enddo
  do i=1,dim_obs_f
    ! print *,'m_id_f_tmp(i) ', m_id_f_tmp(i)
     m_id_f(m_id_f_tmp(i)) = i
  enddo

  ! Clean up
  DEALLOCATE(m_state_tmp, local_dis,m_id_p_tmp,m_id_f_tmp)

END SUBROUTINE obs_op_f_pdaf
