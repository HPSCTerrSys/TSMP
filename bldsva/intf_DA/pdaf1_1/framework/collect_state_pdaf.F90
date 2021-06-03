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
!collect_state_pdaf.F90: TerrSysMP-PDAF implementation of routine
!                        'collect_state_pdaf' (PDAF online coupling)
!-------------------------------------------------------------------------------------------

!$Id: collect_state_pdaf.F90 1441 2013-10-04 10:33:42Z lnerger $
!BOP
!
! !ROUTINE: collect_state_pdaf --- Initialize state vector from model fields
!
! !INTERFACE:
SUBROUTINE collect_state_pdaf(dim_p, state_p)

! !DESCRIPTION:
! User-supplied routine for PDAF.
! Used in the filters: SEEK/EnKF/SEIK/LSEIK/ETKF/LETKF/ESTKF/LESTKF
!
! This subroutine is called during the forecast 
! phase from PDAF\_put\_state\_X after the 
! propagation of each ensemble member. 
! The supplied state vector has to be initialized
! from the model fields (typically via a module). 
! With parallelization, MPI communication might be 
! required to initialize state vectors for all 
! subdomains on the model PEs. 
!
! The routine is executed by each process that is
! participating in the model integrations.
!
! !REVISION HISTORY:
! 2004-11 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
!    use mod_tsmp, &
!        only: pf_statevec_fortran, tag_model_parflow, tag_model_clm
!    use mod_parallel_model, &
!        only: model,mype_world,mype_model
!    USE spmdMod      , only: masterproc
   USE mod_assimilation, &
        ONLY: obs_index_p, dim_obs,obs_filename,obs_index_p_TB,&
              dim_state,dim_state_p_count, dim_state_p_stride, &
              ens_TB, step_TB
   use mod_parallel_model, &
        only: model,mype_model,npes_model,mype_world,npes_world,&
              COMM_model,mpi_comm_world
   use mod_tsmp, &
        only: tag_model_clm,tag_model_parflow,tag_model_cosmo,&
              pf_statevec_fortran, &
              nprocpf,nprocclm,lcmem !SPo add lcmem

   ! LSN: module load for the implementation of CMEM model
   USE mod_parallel_pdaf, &     ! Parallelization variables fro assimilation
        ONLY:n_modeltasks,task_id,COMM_filter,mype_filter,npes_filter,&
            MPI_DOUBLE_PRECISION,MPI_INTEGER,filterpe,&
            MPI_SUCCESS,COMM_couple,mype_couple
   USE spmdMod      , only: masterproc,iam,mpicom,npes
   use clm4cmem     , only: SATELLITE,CLM_DATA
   USE rdclm4pdaf   , only: read_CLM_pdaf
   USE get_tb_cmem  , only: cmem_main
   USE rdclm_wrcmem , only: read_satellite_info
   USE YOMCMEMPAR   , only: INPUTNAMLST,LGPRINT
   use clm_time_manager, only : get_nstep   

#if defined CLMSA
    !kuw: get access to clm variables
    USE clmtype     , only : clm3
    USE clm_varpar  , only : nlevsoi
    use shr_kind_mod, only: r8 => shr_kind_r8
    use enkf_clm_mod, only: clm_statevec
    !kuw end
#endif
    ! To set the PDAF state to COSMO state
#if (defined COUP_OAS_COS) && (!defined COUP_OAS_PFL)
    USE enkf_cosmo_mod, ONLY: cos_statevec
#endif

  IMPLICIT NONE
  
! !ARGUMENTS:
  INTEGER, INTENT(in) :: dim_p           ! PE-local state dimension
  REAL, INTENT(inout) :: state_p(dim_p)  ! local state vector

  ! local variables
 ! integer i
  INTEGER  :: step               ! Currrent time step
!  INTEGER, INTENT(in)   :: dim_p              ! PE-local dimension of state
!  INTEGER, INTENT(in)   :: dim_obs_p          ! Dimension of observed state
!  REAL, INTENT(in)      :: state_p(dim_p)     ! PE-local model state
!  REAL, INTENT(out)     :: m_state_p(dim_obs_p) ! PE-local observed state

  character*200         :: inparam_fname
  TYPE(SATELLITE),ALLOCATABLE :: SAT
  TYPE(CLM_DATA), ALLOCATABLE :: CLMVARS
  REAL,DIMENSION(1),ALLOCATABLE :: TB(:)
  REAL,DIMENSION(1),ALLOCATABLE :: ens_TB_tmp(:)
  !REAL,ALLOCATABLE      :: ens_TB(dim_obs,n_modeltask)
  INTEGER               :: i,j, nerror,nproc,n_obs
  CHARACTER (len = 110) :: OBSFILE

#if defined CLMSA
  real(r8), pointer :: sw_c(:,:)
#endif
! !CALLING SEQUENCE:
! Called by: PDAF_put_state_X   (as U_coll_state)
!EOP
  

! *************************************************
! *** Initialize state vector from model fields ***
! *************************************************
 !call mpi_bcast(dim_obs, 1, MPI_INTEGER, 3, mpi_comm_world, nerror)

 if (model == tag_model_parflow) then
     !print *, "Parflow: collect_state_pdaf, from subvecs to state_p"
     state_p = pf_statevec_fortran
 end if

 IF (lcmem) THEN
    
    write(*,*) 'collect_state_pdaf:',masterproc,mype_world,mype_model
    ! LSN: the implementation of CMEM model here
    IF (model == tag_model_clm) THEN

     LGPRINT = .False.

     ! 1. read 'input' for CMEM namelist
     inparam_fname = './obs/input'
     INPUTNAMLST   = trim(inparam_fname)

     ! 2. assign SAT array from observation     
     !step = get_nstep()-1
     WRITE(OBSFILE,'(a, i5.5)') trim(obs_filename)//'.',step_TB
     WRITE(*,*) 'The currentate_pdaf: observation file is ',OBSFILE     
     allocate(SAT)
     call read_satellite_info(OBSFILE,SAT)

     ! 3. assign CLMVARS array and write forcings as .nc file
     allocate(CLMVARS)
     call read_CLM_pdaf(LS=CLMVARS,SAT=SAT)
     WRITE(*,*) 'collect_state_pdaf: Read CLM memory is done'

     ! 4. the CMEM forwarding part
     IF (masterproc) THEN
   
         n_obs  = size(SAT%lon_foprt)
         write(*,*) 'collect_state_pdaf: n_obs is',n_obs
         ALLOCATE(ens_TB_tmp(n_obs*n_modeltasks))
         ALLOCATE(TB(n_obs))
         call cmem_main(LS=CLMVARS,SATinfo=SAT,TB=TB,n_obs=n_obs,step=step_TB)
         write(*,*) 'collect_state_pdaf: TB are',TB
 
         CALL MPI_Barrier(COMM_couple, nerror)
         CALL MPI_Gather(TB, n_obs,MPI_DOUBLE_PRECISION,ens_TB_tmp,n_obs,&
                   MPI_DOUBLE_PRECISION,0,COMM_couple,nerror)
         write(*,*) 'collect_state_pdaf: Gather parameters of',n_obs,n_modeltasks
         write(*,*) 'collect_state_pdaf: mype_filter is',mype_filter
         write(*,*) 'collect_state_pdaf: mype_couple is',mype_couple
         write(*,*) 'collect_state_pdaf: mype_world is',mype_world

      END IF
   
      !write(*,*) 'collect_state_pdaf: Gather parameters of',dim_obs,n_modeltasks
      !write(*,*) 'collect_state_pdaf: mype_filter is',mype_filter
      !write(*,*) 'collect_state_pdaf: mype_couple is',mype_couple
      !write(*,*) 'collect_state_pdaf: mype_world is',mype_world

      if (mype_world == nprocpf) then
     
         write(*,*) 'collect_state_pdaf: ens_TB_tmp is',ens_TB_tmp  
         ! IF(.NOT.ALLOCATED(ens_TB)) ALLOCATE(ens_TB(n_obs,n_modeltasks))
         IF (ALLOCATED(ens_TB)) DEALLOCATE(ens_TB)
         ALLOCATE(ens_TB(n_obs,n_modeltasks))
         do i = 1,n_obs
            do j = 1,n_modeltasks
               ens_TB(i,j) = ens_TB_tmp((j-1)*n_obs+i)
            end do
         end do
      end if

   IF (ALLOCATED(SAT)) DEALLOCATE(SAT)
   IF (ALLOCATED(CLMVARS)) DEALLOCATE(CLMVARS)
   IF (ALLOCATED(TB)) DEALLOCATE(TB)
   IF (ALLOCATED(ens_TB_tmp)) DEALLOCATE(ens_TB_tmp) 

   END IF
 END IF

 CALL MPI_Barrier(mpi_comm_world, nerror)

#if defined CLMSA
 !kuw: define state vector for clm
 if (model == tag_model_clm) then
     !sw_c  => clm3%g%l%c%cws%h2osoi_vol
     !state_p = reshape(sw_c,shape(state_p))
!     print *, "lsn: check here"     
!print *, clm3%g%l%c%cws%h2osoi_vol
     state_p = clm_statevec
 end if
 !kuw end
#endif
    ! To set the PDAF state to COSMO state
#if (defined COUP_OAS_COS) && (!defined COUP_OAS_PFL)
    if (model == tag_model_cosmo) then
        state_p = cos_statevec
    end if
#endif

  
END SUBROUTINE collect_state_pdaf
