! Copyright 2006-2014 ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities granted to it by
! virtue of its status as an intergovernmental organisation nor does it submit to any jurisdiction.
!

SUBROUTINE ABOR1(CDTEXT)


! SCM dummy: exit with error


USE YOMLUN   , ONLY : NULOUT


IMPLICIT NONE


CHARACTER(LEN=*)                 :: CDTEXT ! UNDETERMINED INTENT


WRITE(NULOUT,'(1X,A)') CDTEXT


CLOSE(NULOUT)


STOP 1

END SUBROUTINE ABOR1

