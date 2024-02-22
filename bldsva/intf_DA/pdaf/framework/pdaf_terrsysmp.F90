!-------------------------------------------------------------------------------------------
!Copyright (c) 2013-2016 by Wolfgang Kurtz and Guowei He (Forschungszentrum Juelich GmbH)
!
!This file is part of TSMP-PDAF
!
!TSMP-PDAF is free software: you can redistribute it and/or modify
!it under the terms of the GNU Lesser General Public License as published by
!the Free Software Foundation, either version 3 of the License, or
!(at your option) any later version.
!
!TSMP-PDAF is distributed in the hope that it will be useful,
!but WITHOUT ANY WARRANTY; without even the implied warranty of
!MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!GNU LesserGeneral Public License for more details.
!
!You should have received a copy of the GNU Lesser General Public License
!along with TSMP-PDAF.  If not, see <http://www.gnu.org/licenses/>.
!-------------------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------------------
!pdaf_terrsysmp.F90: Driver program for TSMP-PDAF
!-------------------------------------------------------------------------------------------

!> @author Wolfgang Kurtz, Guowei He
!> @date 28.02.2022
!> @brief Main program for TSMP-PDAF
!> @details
!> Main TSMP-PDAF program.
program pdaf_terrsysmp

    use mpi, &
      only: MPI_INIT, MPI_FINALIZE
    !, MPI_BARRIER, MPI_COMM_WORLD, MPI_SUCCESS
  
    use mod_parallel_pdaf, &
        only : mype_world, MPIerr

    use mod_tsmp, &
      only: initialize_tsmp, integrate_tsmp, update_tsmp, finalize_tsmp, &
      tcycle, total_steps

    use mod_assimilation, &
      only: screen

    implicit none

    ! initialize mpi
    call MPI_INIT(MPIerr)

    ! intitialize parallel pdaf (communicators et al.)
    call init_parallel_pdaf(0, 1)

    ! Read TSMP-PDAF input from "enkfpf.par"
    ! initialize TSMP instances
    call initialize_tsmp()

    ! initialize pdaf variables
    call init_pdaf()

    ! barrier before model integration starts
    ! call MPI_BARRIER(MPI_COMM_WORLD, MPIerr)
    ! if (MPIerr .ne. MPI_SUCCESS) then
    !     print *, "barrier failed"
    ! end if

    ! time loop
    do tcycle = 1, total_steps
        if (mype_world > -1 .and. screen > 2) then
            print *, "TSMP-PDAF mype(w)=", mype_world, ": time loop", tcycle
        endif

        ! forward simulation of component models
        call integrate_tsmp()

        ! assimilation step
        call assimilate_pdaf()

        !call MPI_BARRIER(MPI_COMM_WORLD, MPIerr)
        !print *,"Finished assimilation", tcycle

        !call print_update_pfb()
        call update_tsmp()

        !call MPI_BARRIER(MPI_COMM_WORLD, MPIerr)
        !print *,"Finished complete assimilation cycle", tcycle
    enddo

    ! barrier after model integrations
    !call MPI_BARRIER(MPI_COMM_WORLD, MPIerr)
    !print *, "pdaf: advancing finished, rank ", mype_world

    ! close pdaf
    call finalize_pdaf()
    !print *, "pdaf: finalized, rank ", mype_world

    ! close model instances
    call finalize_tsmp()
    !call MPI_BARRIER(MPI_COMM_WORLD, MPIerr)
    !print *, "model: finalized, rank ", mype_world

    ! close mpi
    call MPI_FINALIZE(MPIerr)

end program pdaf_terrsysmp
