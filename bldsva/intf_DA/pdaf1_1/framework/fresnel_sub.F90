! Copyright 2006-2014 ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities granted to it by
! virtue of its status as an intergovernmental organisation nor does it submit to any jurisdiction.
!
SUBROUTINE FRESNEL (eps_surf)

! Purpose :
!  Fresnel Law to computes the reflectivities of flat surfaces
! eps_surf : dielectric constant of the surface
!---------------------------------------------------------------------------

USE YOMCMEMPAR ,ONLY : sintheta, costheta, LGPRINT
USE YOMCMEMSOIL, ONLY : r_s

IMPLICIT NONE

COMPLEX :: g
COMPLEX :: eps_surf
!---------------------------------------------------------------------------

g = sqrt( eps_surf - sintheta*sintheta )

r_s(1) = abs((costheta-g)/(costheta+g)) ** 2.
r_s(2) = abs((costheta*eps_surf-g)/(costheta*eps_surf+g)) ** 2.

END SUBROUTINE FRESNEL
