! Copyright 2006-2014 ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities granted to it by
! virtue of its status as an intergovernmental organisation nor does it submit to any jurisdiction.
!
SUBROUTINE IO_CMEMNETCDF(CLNAME,LD_firstcall,LD_infocoord)

! PURPOSE
! -------
!  To read data from / to netcdf files
   
! INTERFACE
! ---------
!  CALL IO_CMEMNETCDF

! METHOD
! ------
!
! EXTERNALS
! ---------
! netcdf: nf90_inq_varid, nf90_inq_varid, nf90_get_var
!
! INTERNAL ROUTINES
! -----------------
!  Patricia de Rosnay  December 2007
!------------------------------------------------------------------------------

USE YOMCMEMNETCDF

IMPLICIT NONE

CHARACTER(len=80) :: clname
LOGICAL :: LD_firstcall, LD_infocoord
!------------------------------------------------------------------------------


    CALL INNETCDF(CLNAME,LD_firstcall,LD_infocoord)

RETURN

CONTAINS

!------------------------------------------------------------------------------

SUBROUTINE INNETCDF(INFILE,LD_firstcall, LD_infocoord)

USE netcdf
USE YOMLUN, ONLY : NULOUT
USE YOMCMEMNETCDF
USE YOMCMEMPAR, ONLY : LGPRINT
!
IMPLICIT NONE

INTEGER(KIND=JPIM) :: nlon_varid, nlat_varid,ntime_varid,nlvl_varid
INTEGER(KIND=JPIM) :: field_varid
INTEGER(KIND=JPIM) :: ncid,status
INTEGER(KIND=JPIM) :: ilvl, ilat, ilon, irec, i
INTEGER(KIND=JPIM), DIMENSION(nf90_max_var_dims) :: DimIds

CHARACTER(len=80) infile

LOGICAL :: LD_firstcall, LD_infocoord
REAL, ALLOCATABLE :: var_in_tr(:,:,:,:)



! 1.  Open the input file
!------------------------

 WRITE(NULOUT,*)'INPUT_NC, read file  ',trim(infile)
 status=nf90_open(INFILE, nf90_nowrite, ncid)
 IF(status /= nf90_noerr) CALL ABOR1(trim(nf90_strerror(status)))
!
!
! 2. Get the varids and dimension of the field netCDF variables
!----------------------------------------------------------------
!
! Default dimension:
  NLATS = 1_JPIM
  NLONS = 1_JPIM
  NTIMES = 1_JPIM
  NLVLS = 1_JPIM
  NDIMS = 1_JPIM
  DimIds(:) = 0_JPIM
!
  status =  nf90_inq_varid(ncid, CCVARNC_NAME, field_varid)
  IF(status /= nf90_noerr) CALL ABOR1(trim(nf90_strerror(status)))
!
  status = nf90_inquire_variable(ncid, field_varid, ndims = nDims)
  IF(status /= nf90_noerr) CALL ABOR1(trim(nf90_strerror(status)))
!
  IF (LGPRINT) WRITE(NULOUT,*) 'ndims',ndims
!
  status = nf90_inquire_variable(ncid, field_varid,dimids = DimIds)
  IF(status /= nf90_noerr) CALL ABOR1(trim(nf90_strerror(status)))
!
  IF (LGPRINT) WRITE(NULOUT,*) 'DimIds',DimIds

SELECT CASE (DimIds(1) /= 0)
  CASE (.TRUE.)
   status = nf90_inquire_dimension(ncid,dimids(1),len=nLons)
   if (status /= nf90_noerr) CALL ABOR1(trim(nf90_strerror(status))) 
   IF (LGPRINT) WRITE(NULOUT,*) 'nLons',nLons
END SELECT
!
SELECT CASE (DimIds(2) /= 0)
  CASE (.TRUE.)
   status = nf90_inquire_dimension(ncid,dimids(2),len=nlats)
   if (status /= nf90_noerr) CALL ABOR1(trim(nf90_strerror(status))) 
   IF (LGPRINT) WRITE(NULOUT,*) 'nlats',nlats
END SELECT
!
SELECT CASE (DimIds(4) /= 0)
  CASE (.TRUE.)
   status = nf90_inquire_dimension(ncid,dimids(4),len=ntimes)
   if (status /= nf90_noerr) CALL ABOR1(trim(nf90_strerror(status))) 
   IF (LGPRINT) WRITE(NULOUT,*) 'ntimes',ntimes
END SELECT
!
SELECT CASE (DimIds(3) /= 0)
  CASE (.TRUE.)
   status = nf90_inquire_dimension(ncid,dimids(3),len=nlvls)
   if (status /= nf90_noerr) CALL ABOR1(trim(nf90_strerror(status))) 
   IF (LGPRINT) WRITE(NULOUT,*)'nlvls',nlvls
END SELECT
!
!
! 3. Allocate and Read input variable according to its dimensions 
!-----------------------------------------------------------------
!
SELECT CASE (LD_firstcall)
!
  CASE(.False.)
!
   IF( ALLOCATED (fncfield)) DEALLOCATE (fncfield)
   ALLOCATE(fncfield(NLONS*NLATS*NTIMES,NLVLS))
!
   IF( ALLOCATED (var_in)) DEALLOCATE (var_in)
   ALLOCATE(var_in(NLONS,NLATS,NLVLS,NTIMES))
   IF( ALLOCATED (var_in_tr)) DEALLOCATE (var_in_tr)
   ALLOCATE(var_in_tr(NLONS,NLATS,NTIMES,NLVLS))
!
   status=nf90_get_var(ncid, field_varid, var_in)
   if (status /= nf90_noerr) CALL ABOR1(trim(nf90_strerror(status)))
!
   do ilvl = 1, nlvls
   do irec = 1,ntimes
    var_in_tr (:,:,irec,ilvl) = var_in(:,:,ilvl,irec)
   enddo
   enddo
!
! reshape the input variable before it goes to CMEM
!
  fncfield = reshape(var_in_tr, (/NLONS*NLATS*NTIMES,NLVLS/))
!
END SELECT
!
!
! 4.  Get latitude, longitude, time, levels coordinate variables
!---------------------------------------------------------------
!
SELECT CASE (LD_infocoord)
!
  CASE (.True.)
  IF (LGPRINT) WRITE(NULOUT,*) 'Read coordinates informations'

   NLONS_SM = NLONS
   NLATS_SM = NLATS
   NLVLS_SM = NLVLS
   NTIMES_SM = NTIMES
   NDIMS_SM = NDIMS
!
   SELECT CASE (DimIds(1) /= 0)
     CASE (.TRUE.)
!    4.1 Get longitude info
!
       status=nf90_inq_varid(ncid, CCLON_NAME, nlon_varid)
       IF(status /= nf90_noerr)  CALL ABOR1(trim(nf90_strerror(status)))
!
       IF( ALLOCATED (xlons)) DEALLOCATE (xlons)
       ALLOCATE(xlons(NLONS))
       status=nf90_get_var(ncid, nlon_varid, xlons)
       IF(status /= nf90_noerr) CALL ABOR1(trim(nf90_strerror(status)))
!
   END SELECT
!
   SELECT CASE (DimIds(2) /= 0)
     CASE (.TRUE.)
!    4.2 Get latitude info 
!
       status=nf90_inq_varid(ncid, CCLAT_NAME, nlat_varid)
       IF(status /= nf90_noerr)  CALL ABOR1(trim(nf90_strerror(status)))

       IF( ALLOCATED (xlats)) DEALLOCATE (xlats)
       ALLOCATE(xlats(NLATS))
       status= nf90_get_var(ncid, nlat_varid, xlats)
       IF(status /= nf90_noerr) CALL ABOR1(trim(nf90_strerror(status)))
!
   END SELECT
!
   SELECT CASE (DimIds(3) /= 0)
     CASE (.TRUE.)
!    4.3 Get Level info
!
         status=nf90_inq_varid(ncid, CCLVL_NAME, nlvl_varid)
         IF(status /= nf90_noerr)  CALL ABOR1(trim(nf90_strerror(status)))
!
         IF( ALLOCATED (xlvls)) DEALLOCATE (xlvls)
         ALLOCATE(xlvls(NLVLS))
         status= nf90_get_var(ncid, nlvl_varid, xlvls)
         IF(status /= nf90_noerr) CALL ABOR1(trim(nf90_strerror(status)))
!
   END SELECT
!
   SELECT CASE (DimIds(4) /= 0)
     CASE (.TRUE.)
!    4.4 Get time info
!
       status=nf90_inq_varid(ncid, CCTIME_NAME, ntime_varid)
       IF(status /= nf90_noerr)  CALL ABOR1(trim(nf90_strerror(status)))
!
       IF( ALLOCATED (xtimes)) DEALLOCATE (xtimes)
       ALLOCATE(xtimes(NTIMES))
       status= nf90_get_var(ncid, ntime_varid, xtimes)
       IF(status /= nf90_noerr) CALL ABOR1(trim(nf90_strerror(status)))
!
   END SELECT
!
!
END SELECT


! 5. Deallocate local variables and close the file. This frees up any internal netCDF resources
!----------------------------------------------------------------------------------------------

   IF( ALLOCATED (var_in)) DEALLOCATE (var_in)
   IF( ALLOCATED (var_in_tr)) DEALLOCATE (var_in_tr)

    status=nf90_close(ncid)
    IF(status /= nf90_noerr) CALL ABOR1(trim(nf90_strerror(status)))


RETURN


END SUBROUTINE INNETCDF

!------------------------------------------------------------------------------


END SUBROUTINE IO_CMEMNETCDF


