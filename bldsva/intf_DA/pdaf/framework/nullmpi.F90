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
!nullmpi.F90: TSMP-PDAF implementation of routine
!             'nullmpi' (PDAF online coupling)
!-------------------------------------------------------------------------------------------

!$Id: nullmpi.F90 1383 2013-05-03 12:26:53Z lnerger $
!BOP
!
! !ROUTINE: mpi_init() --- Pseudo-implementation of MPI_init
!
! !INTERFACE:
SUBROUTINE mpi_init(i)

! !DESCRIPTION:
! This routine simulates MPI functionality for
! a program running on a single processor. Its
! purpose is to avoid the need of a real MPI
! library when running serial jobs. 
!
! !REVISION HISTORY:
! 2004-10 - Lars Nerger - Initial code
! Later revisions - see svn log
!EOP

  IMPLICIT NONE

  INTEGER :: i
  
  i = 0

END SUBROUTINE mpi_init

! ------------------------------------------------------------------------------
SUBROUTINE mpi_finalize(i)

  IMPLICIT NONE

  INTEGER :: i
  
  i=0

END SUBROUTINE mpi_finalize

! ------------------------------------------------------------------------------
SUBROUTINE MPI_Comm_Size(comm, npes_world, i)

  IMPLICIT NONE

  INTEGER :: comm
  INTEGER :: npes_world
  INTEGER :: i

  npes_world = 1
  i = 0

END SUBROUTINE MPI_Comm_Size

! ------------------------------------------------------------------------------
SUBROUTINE MPI_Comm_Rank(comm, mype_world, i)

  IMPLICIT NONE

  INTEGER :: comm
  INTEGER :: mype_world
  INTEGER :: i

  mype_world = 0
  i = 0

END SUBROUTINE MPI_Comm_Rank

! ------------------------------------------------------------------------------
SUBROUTINE MPI_Comm_split(comm_a, my_comm, mype_a, comm_b, mpierr)

  IMPLICIT NONE

  INTEGER :: comm_a
  INTEGER :: my_comm
  INTEGER :: mype_a
  INTEGER :: comm_b
  INTEGER :: MPIerr

  comm_b = 1
  my_comm = 1
  mype_a = 0
  mpierr = 0

END SUBROUTINE MPI_Comm_split

! ------------------------------------------------------------------------------
SUBROUTINE MPI_Barrier(comm, mpierr)

  IMPLICIT NONE

  INTEGER :: comm
  INTEGER :: mpierr

  mpierr = 0

END SUBROUTINE MPI_Barrier

! ------------------------------------------------------------------------------
SUBROUTINE MPI_Send(field, dim, datatype, pe_source, pe_target, &
     comm, mpierr)

  IMPLICIT NONE

  INTEGER :: dim
  INTEGER :: field(dim)
  INTEGER :: datatype
  INTEGER :: pe_source
  INTEGER :: pe_target
  INTEGER :: comm
  INTEGER :: mpierr

  mpierr = 0

END SUBROUTINE MPI_Send

! ------------------------------------------------------------------------------
SUBROUTINE MPI_Recv(field, dim, datatype, pe_source, pe_target, &
     comm, flag, mpierr)

  IMPLICIT NONE

  INTEGER :: dim
  INTEGER :: field(dim)
  INTEGER :: datatype
  INTEGER :: pe_source
  INTEGER :: pe_target
  INTEGER :: comm
  INTEGER :: flag
  INTEGER :: mpierr

  mpierr = 0

END SUBROUTINE MPI_Recv

! ------------------------------------------------------------------------------
SUBROUTINE MPI_BCast(field, dim, datatype, pe_source, comm, mpierr)

  IMPLICIT NONE

  INTEGER :: dim
  INTEGER :: field(dim)
  INTEGER :: datatype
  INTEGER :: pe_source
  INTEGER :: comm
  INTEGER :: mpierr

  mpierr = 0

END SUBROUTINE MPI_BCast

! ------------------------------------------------------------------------------
SUBROUTINE MPI_Allreduce(field_in, field_out, dim, fieldtype, operation, &
     comm, mpierr)

  IMPLICIT NONE

  INTEGER :: dim
  REAL    :: field_in(dim)
  REAL    :: field_out(dim)
  INTEGER :: fieldtype
  INTEGER :: operation
  INTEGER :: comm
  INTEGER :: mpierr

  field_out = field_in
  mpierr = 0

END SUBROUTINE MPI_Allreduce

! ------------------------------------------------------------------------------
SUBROUTINE MPI_Reduce(field_in, field_out, dim, fieldtype, operation, &
     pe_root, comm, mpierr)

  IMPLICIT NONE

  INTEGER :: dim
  REAL    :: field_in(dim)
  REAL    :: field_out(dim)
  INTEGER :: fieldtype
  INTEGER :: operation
  INTEGER :: pe_root
  INTEGER :: comm
  INTEGER :: mpierr

  field_out = field_in
  mpierr = 0

END SUBROUTINE MPI_REDUCE

! ------------------------------------------------------------------------------
SUBROUTINE MPI_Allgather(field_in, dim_in, type_in, field_out, dim_out, &
     type_out, comm, mpierr)

  IMPLICIT NONE

  INTEGER :: dim_in
  REAL    :: field_in(dim_in)
  INTEGER :: type_in
  INTEGER :: dim_out
  REAL    :: field_out(dim_out)
  INTEGER :: dis
  INTEGER :: type_out
  INTEGER :: comm
  INTEGER :: mpierr

  field_out = field_in
  mpierr = 0

END SUBROUTINE MPI_ALLGATHER

! ------------------------------------------------------------------------------
SUBROUTINE MPI_AllGatherV(field_in, dim_in, type_in, field_out, dim_out, &
     dis, type_out, comm, mpierr)

  IMPLICIT NONE

  INTEGER :: dim_in
  REAL    :: field_in(dim_in)
  INTEGER :: type_in
  INTEGER :: dim_out
  REAL    :: field_out(dim_out)
  INTEGER :: dis
  INTEGER :: type_out
  INTEGER :: comm
  INTEGER :: mpierr

  field_out = field_in
  mpierr = 0

END SUBROUTINE MPI_AllGatherV
