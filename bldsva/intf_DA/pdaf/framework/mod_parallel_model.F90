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
!mod_parallel_model.F90: Module for control of parallel models of TerrSysMP
!-------------------------------------------------------------------------------------------

MODULE mod_parallel_model

    USE mpi
    USE iso_c_binding, ONLY: c_int, c_double

    SAVE

    ! mpi related
    INTEGER :: npes_parflow
    INTEGER :: coupcol
#if defined PARFLOW_STAND_ALONE
! Parflow stand alone directly uses binded communicator
    INTEGER(c_int),bind(c,name='comm_model_pdaf') :: comm_model
#else
! CLM stand alone uses comm_model directly, while TerrSysMP use this and da_comm(_clm)
    INTEGER :: comm_model
#endif
    INTEGER(c_int), BIND(c) :: mype_model
    INTEGER(c_int), BIND(c) :: npes_model
    INTEGER(c_int), BIND(c) :: mype_world
    INTEGER(c_int), BIND(c) :: npes_world
    ! model input parameters
    REAL(c_double), BIND(c) :: t_start
    INTEGER(c_int), BIND(c) ::  model
    INTEGER(c_int), BIND(c) :: tcycle
    INTEGER(c_int), BIND(c) :: tstartcycle
    INTEGER(c_int), BIND(c) :: total_steps

    INTERFACE
        SUBROUTINE read_enkfpar(parname) BIND(C, name='read_enkfpar')
            USE iso_c_binding
            IMPLICIT NONE
            CHARACTER, DIMENSION(*), INTENT(in) :: parname
        END SUBROUTINE read_enkfpar
    END INTERFACE
END MODULE mod_parallel_model
