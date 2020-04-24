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
        ONLY: obs_index_p, dim_obs,obs_filename,obs_index_p_TB
   use mod_parallel_model, &
        only: model,mype_model,npes_model,mype_world,npes_world
   use mod_tsmp, &
        only: tag_model_clm,nprocpf 

   ! LSN: module load for the implementation of CMEM model
   USE mod_parallel_pdaf, &     ! Parallelization variables fro assimilation
        ONLY:n_modeltasks,task_id,COMM_filter,mype_filter,npes_filter,&
            MPI_DOUBLE_PRECISION         
   USE spmdMod      , only : masterproc,iam,mpicom,npes
   use clm4cmem     , only: SATELLITE,CLM_DATA
   USE rdclm4pdaf   , only: read_CLM_pdaf
   USE get_tb_cmem  , only: cmem_main
   USE rdclm_wrcmem , only: read_satellite_info
   USE YOMCMEMPAR   , only: INPUTNAMLST, CLMNAME, SURFNAME,LGPRINT
   IMPLICIT NONE
! !ARGUMENTS:
  INTEGER, INTENT(in) :: step               ! Currrent time step
  INTEGER, INTENT(in) :: dim_p              ! PE-local dimension of state
  INTEGER, INTENT(in) :: dim_obs_p          ! Dimension of observed state
  REAL, INTENT(in)    :: state_p(dim_p)     ! PE-local model state
  REAL, INTENT(out) :: m_state_p(dim_obs_p) ! PE-local observed state

  character*200:: CLM_fname, inparam_fname, surf_fname
  TYPE(SATELLITE),ALLOCATABLE :: SAT
  TYPE(CLM_DATA), ALLOCATABLE :: CLMVARS             
  REAL,DIMENSION(1) :: TB(dim_obs)  
  INTEGER :: i, nerror
  CHARACTER (len = 110) :: current_observation_filename
! type(SATELLITE), intent(out) :: SAT 
! !CALLING SEQUENCE:
! Called by: PDAF_seek_analysis   (as U_obs_op)
! Called by: PDAF_seik_analysis, PDAF_seik_analysis_newT
! Called by: PDAF_enkf_analysis_rlm, PDAF_enkf_analysis_rsm
!EOP

! *** local variables ***

! *********************************************
! *** Perform application of measurement    ***
! *** operator H on vector or matrix column ***
! *********************************************
  ! LSN: the original observation operator of soil moisture in Parflow
  !DO i = 1, dim_obs_p
  !   m_state_p(i) = state_p(obs_index_p(i))
  !END DO
  
  ! LSN: the implementation of CMEM model here
  IF (model == tag_model_clm) THEN
     
     allocate(CLMVARS)  ! assign CLMVARS array
     surf_fname    = './obs/surf_fname.nc'
     SURFNAME      = trim(surf_fname)
     CLM_fname     = './obs/CLM_fname.nc'
     CLMNAME       = trim(CLM_fname)     
     inparam_fname = './obs/input'
     INPUTNAMLST   = trim(inparam_fname)

     LGPRINT = .True.
  
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
  CALL MPI_Bcast(TB, dim_obs, MPI_DOUBLE_PRECISION, nprocpf, COMM_filter, nerror)
  
  DO i = 1, dim_obs_p
     m_state_p(i) = TB(obs_index_p_TB(i))
  END DO
 
  IF (ALLOCATED(SAT)) DEALLOCATE(SAT)
  IF (ALLOCATED(CLMVARS)) DEALLOCATE(CLMVARS)
  
  !LSN: the end of the implementation of CMEM model here
END SUBROUTINE obs_op_pdaf
