! Copyright 2006-2014 ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities granted to it by
! virtue of its status as an intergovernmental organisation nor does it submit to any jurisdiction.
!
SUBROUTINE WRCMEMNETCDF

! Purpose :
! -------
!  WRITE CMEM OUTPUT DATA IN NETCDF FORMAT
   
! Interface :
! ---------

! Method :
! ------

! Externals :
! ---------

! Internal variables
! ------------------

! Author :
!     Patricia de Rosnay  January 2008

!------------------------------------------------------------------------------
!
USE netcdf
!
USE PARKIND1, ONLY : JPIM, JPRM 
USE YOMLUN, ONLY : NULOUT,NULTMP
!
USE YOMCMEMPAR, ONLY : JPHISTLEV,CNAMEID,CFREQ,CANGLE
USE YOMCMEMFIELDS, ONLY: N,JJ,CLNAME,fwc_veg,fb,ftfrac,fh,fWP,ftau_atm,ftb_au &
                     & ,ftb_toa,fteffC,fsurf_emis,ftau_veg
USE YOMCMEMNETCDF, ONLY : NTIMES_SM,NLATS_SM, NLONS_SM, NLVLS_SM,NTIMES,NLATS, NLONS, NLVLS &
                     & , CCLVL_NAME,CCLAT_NAME,CCLON_NAME,CCTIME_NAME,CCLAT_UNITS,CCLON_UNITS &
                     & ,CCTIME_UNITS,CCTEMP_UNITS,CCVWC_UNITS,CCUNITS,CCPFILED_UNITS &
                     & ,xlats,xlons,xlvls,xtimes
!
USE YOMCMEMATM, ONLY : fZ
!
!
IMPLICIT NONE
!
INTEGER(KIND=JPIM) :: ncid,status
INTEGER(KIND=JPIM) :: nlon_varid, nlat_varid,ntime_varid,nlvl_varid
INTEGER(KIND=JPIM) :: londimid,latdimid,timedimid,lvldimid
INTEGER(KIND=JPIM) :: itime
INTEGER(KIND=JPIM), DIMENSION(4) :: DimIds
INTEGER(KIND=JPIM) :: itbv_varid, itbh_varid, iteff_varid &
                     & ,iftauvegh_varid,iftauvegv_varid &
                     & ,ifbar_varid, ivwc_varid, itauatm_varid,itbatmup_varid &
                     & ,ibarenosnow_varid,ilownosnow_varid, ivwcl_varid, ivwch_varid &
                     & ,ibparaml_varid, ibparamh_varid, irugoh_varid, iwiltsm_varid &
                     & , ieh_varid, iev_varid, icparam_varid
!
! LEVEL1 OUTPUTS
CHARACTER (LEN = *), PARAMETER :: CLTBH_NAME="TBH"
CHARACTER (LEN = *), PARAMETER :: CLTBV_NAME="TBV"
CHARACTER (LEN = *), PARAMETER :: CLTEFF_NAME="EFFECTIVE_TEMP"
!
! LEVEL2 OUTPUTS
CHARACTER (LEN = *), PARAMETER :: CLFTAUVEGH_NAME="TAU_VEG_H"
CHARACTER (LEN = *), PARAMETER :: CLFTAUVEGV_NAME="TAU_VEG_V"
CHARACTER (LEN = *), PARAMETER :: CLFBAR_NAME="BARE_FRACT"
CHARACTER (LEN = *), PARAMETER :: CLVWC_NAME="VEG_WATER_CONTENT"
CHARACTER (LEN = *), PARAMETER :: CLTAUATM_NAME="TAU_ATM"
CHARACTER (LEN = *), PARAMETER :: CLTBATMUP_NAME="TB_ATM_UP"
!
! LEVEL3 OUTPUTS
CHARACTER (LEN = *), PARAMETER :: CLBARENOSNOW_NAME="BARE_FRACT_NOSNOW"
CHARACTER (LEN = *), PARAMETER :: CLLOWNOSNOW_NAME="LOW_VEG_FRACT_NOSNOW"
CHARACTER (LEN = *), PARAMETER :: CLVWCL_NAME="LOW_VEG_WATER_CONTENT"
CHARACTER (LEN = *), PARAMETER :: CLVWCH_NAME="HIGH_VEG_WATER_CONTENT"
CHARACTER (LEN = *), PARAMETER :: CLBPARAML_NAME="LOW_VEG_B_PARAM"
CHARACTER (LEN = *), PARAMETER :: CLBPARAMH_NAME="HIGH_VEG_B_PARAM"
CHARACTER (LEN = *), PARAMETER :: CLRUGOH_NAME="RUGO_H"
CHARACTER (LEN = *), PARAMETER :: CLWILTSM_NAME="WILTING_PT_SM"
CHARACTER (LEN = *), PARAMETER :: CLEH_NAME="EMIS_H"
CHARACTER (LEN = *), PARAMETER :: CLEV_NAME="EMIS_V"
CHARACTER (LEN = *), PARAMETER :: CLCPARAM_NAME="C_PARAM_TEFF"


REAL(KIND=JPRM), DIMENSION(NLONS_SM,NLATS_SM,NLVLS_SM,NTIMES_SM) :: tbv_out,tbh_out,teff_out &
                & , ftauvegh_out,ftauvegv_out,fbar_out,vwc_out,tauatm_out,tbatmup_out &
                & , barenosnow_out, lownosnow_out, vwcl_out, vwch_out &
                & , bparaml_out, bparamh_out, rugoh_out, wiltsm_out &
                & , eh_out, ev_out, cparam_out
 
REAL(KIND=JPRM), DIMENSION(NLONS_SM,NLATS_SM,NTIMES_SM) :: var_in_tr
REAL(KIND=JPRM), DIMENSION(N) :: fvar
!
!------------------------------------------------------------------------------
!
!
!
NLONS = NLONS_SM
NLATS = NLATS_SM
NLVLS = 1_JPIM
NTIMES = NTIMES_SM
!
!
! 1. LEVEL 1 outputs: first netcdf file
!---------------------------------------
!
! 1.0 put variable in right format 
!---------------------------------
!
!  TBH
var_in_tr = reshape(ftb_toa(:,1), (/NLONS,NLATS,NTIMES/))
tbh_out(:,:,1,:) = var_in_tr(:,:,:) 
!
!  TBV
var_in_tr = reshape(ftb_toa(:,2), (/NLONS,NLATS,NTIMES/))
tbv_out(:,:,1,:) = var_in_tr(:,:,:) 
!
!  TEFF
var_in_tr = reshape(fteffC(:,2), (/NLONS,NLATS,NTIMES/))
teff_out(:,:,1,:) = var_in_tr(:,:,:) 
!
!
! 1.1 Open the output file
!-------------------------
!
CLNAME='out_level1_'//CNAMEID//'_'//cfreq//'_'//cangle//'.nc'
status = NF90_CREATE(CLNAME,NF90_CLOBBER,ncid)
IF(status /= nf90_noerr) CALL ABOR1(trim(nf90_strerror(status)))
!
!
! 1.2 Define file dimensions
!---------------------------
!
status = NF90_DEF_DIM(ncid,CCLON_NAME,NLONS,londimid)
IF(status /= nf90_noerr) CALL ABOR1(trim(nf90_strerror(status)))
!
status = NF90_DEF_DIM(ncid,CCLAT_NAME,NLATS,latdimid)
IF(status /= nf90_noerr) CALL ABOR1(trim(nf90_strerror(status)))
!
status = NF90_DEF_DIM(ncid,CCLVL_NAME,NLVLS,lvldimid)
IF(status /= nf90_noerr) CALL ABOR1(trim(nf90_strerror(status)))
!
status = NF90_DEF_DIM(ncid, CCTIME_NAME, NF90_UNLIMITED, timedimid)
IF(status /= nf90_noerr) CALL ABOR1(trim(nf90_strerror(status)))
!
!status = NF90_DEF_DIM(ncid, CCTIME_NAME, NTIMES, timedimid)
IF(status /= nf90_noerr) CALL ABOR1(trim(nf90_strerror(status)))
!
!
status = NF90_DEF_VAR(ncid, CCLON_NAME, NF90_REAL, londimid, nlon_varid)
IF(status /= nf90_noerr) CALL ABOR1(trim(nf90_strerror(status)))
!
status = NF90_DEF_VAR(ncid, CCLAT_NAME, NF90_REAL, latdimid, nlat_varid)
IF(status /= nf90_noerr) CALL ABOR1(trim(nf90_strerror(status)))
!
status = NF90_DEF_VAR(ncid, CCLVL_NAME, NF90_REAL, lvldimid, nlvl_varid)
IF(status /= nf90_noerr) CALL ABOR1(trim(nf90_strerror(status)))
!
status = NF90_DEF_VAR(ncid, CCTIME_NAME, NF90_REAL, timedimid, ntime_varid)
IF(status /= nf90_noerr) CALL ABOR1(trim(nf90_strerror(status)))
!
!
! Assign units attributes to coordinate variables.
!
status = nf90_put_att(ncid, nlat_varid, CCUNITS, CCLAT_UNITS)
IF(status /= nf90_noerr) CALL ABOR1(trim(nf90_strerror(status)))
!
status = nf90_put_att(ncid, nlon_varid, CCUNITS, CCLON_UNITS)
IF(status /= nf90_noerr) CALL ABOR1(trim(nf90_strerror(status)))
!
status = nf90_put_att(ncid, ntime_varid, CCUNITS, CCTIME_UNITS)
IF(status /= nf90_noerr) CALL ABOR1(trim(nf90_strerror(status)))
!
dimids = (/ londimid, latdimid, lvldimid, timedimid /)
!
!
!
! 1.3 Define the netCDF variables for the output data
!----------------------------------------------------
!
status = nf90_def_var(ncid, CLTBH_NAME, NF90_REAL, dimids, itbh_varid) 
IF(status /= nf90_noerr) CALL ABOR1(trim(nf90_strerror(status)))
!
status = nf90_def_var(ncid, CLTBV_NAME, NF90_REAL, dimids, itbv_varid) 
IF(status /= nf90_noerr) CALL ABOR1(trim(nf90_strerror(status)))
!
status = nf90_def_var(ncid, CLTEFF_NAME, NF90_REAL, dimids, iteff_varid) 
IF(status /= nf90_noerr) CALL ABOR1(trim(nf90_strerror(status)))
!
!
! Assign units attributes to the netCDF variables.
!
status = nf90_put_att(ncid, itbh_varid, CCUNITS, CCTEMP_UNITS)
IF(status /= nf90_noerr) CALL ABOR1(trim(nf90_strerror(status)))
!
status = nf90_put_att(ncid, itbv_varid, CCUNITS, CCTEMP_UNITS)
IF(status /= nf90_noerr) CALL ABOR1(trim(nf90_strerror(status)))
!
status = nf90_put_att(ncid, iteff_varid, CCUNITS, CCTEMP_UNITS)
IF(status /= nf90_noerr) CALL ABOR1(trim(nf90_strerror(status)))
!
!
! 1.4 End define mode
!--------------------


status =  NF90_ENDDEF(ncid)
IF(status /= nf90_noerr) CALL ABOR1(trim(nf90_strerror(status)))

! 1.5 Write the coordinate variable data (latitudes and longitudes) into the netCDF file
!---------------------------------------------------------------------------------------

status =  NF90_PUT_VAR(ncid, nlat_varid, xlats) 
IF(status /= nf90_noerr) CALL ABOR1(trim(nf90_strerror(status)))

status =  NF90_PUT_VAR(ncid, nlon_varid, xlons) 
IF(status /= nf90_noerr) CALL ABOR1(trim(nf90_strerror(status)))
  
status =  NF90_PUT_VAR(ncid, ntime_varid, xtimes) 
IF(status /= nf90_noerr) CALL ABOR1(trim(nf90_strerror(status)))
  

! 1.6 Write the data 
!-------------------
!
status =  NF90_PUT_VAR(ncid, itbh_varid, tbh_out) 
IF(status /= nf90_noerr) CALL ABOR1(trim(nf90_strerror(status)))
!
status =  NF90_PUT_VAR(ncid, itbv_varid, tbv_out) 
IF(status /= nf90_noerr) CALL ABOR1(trim(nf90_strerror(status)))
!
status =  NF90_PUT_VAR(ncid, iteff_varid, teff_out) 
IF(status /= nf90_noerr) CALL ABOR1(trim(nf90_strerror(status)))
!
! 1.7 Close the file
!-------------------
!
status = nf90_close(ncid) 
IF(status /= nf90_noerr) CALL ABOR1(trim(nf90_strerror(status)))
!
!
!
!
! 2.0 LEVEL 2 Outputs: second netcdf file
!----------------------------------------
!
SELECT CASE (JPHISTLEV)
!
  CASE ( 2_JPIM, 3_JPIM )
!
! 2.0 put variable in right format 
!---------------------------------
!
!  FTAUVEGH 
   fvar(:) = ftau_veg(:,1)
   var_in_tr = reshape(fvar(:), (/NLONS,NLATS,NTIMES/))
   ftauvegh_out(:,:,1,:) = var_in_tr(:,:,:) 
!  FTAUVEGV 
   fvar(:) = ftau_veg(:,2)
   var_in_tr = reshape(fvar(:), (/NLONS,NLATS,NTIMES/))
   ftauvegv_out(:,:,1,:) = var_in_tr(:,:,:) 
!  FBAR 
   fvar(:) = ftfrac(:,1)+ftfrac(:,2) 
   var_in_tr = reshape(fvar(:), (/NLONS,NLATS,NTIMES/))
   fbar_out(:,:,1,:) = var_in_tr(:,:,:) 
!
!  VWC 
   fvar(:) = fwc_veg(:,1)*(ftfrac(:,3)+ftfrac(:,4)) +  fwc_veg(:,2)*(ftfrac(:,5)+ftfrac(:,6))
   var_in_tr = reshape(fvar(:), (/NLONS,NLATS,NTIMES/))
   vwc_out(:,:,1,:) = var_in_tr(:,:,:) 
!
!  TAU_ATM 
   var_in_tr = reshape(ftau_atm(:), (/NLONS,NLATS,NTIMES/))
   tauatm_out(:,:,1,:) = var_in_tr(:,:,:) 
!
!  TB_ATM_UP 
   var_in_tr = reshape(ftb_au(:), (/NLONS,NLATS,NTIMES/))
   tbatmup_out(:,:,1,:) = var_in_tr(:,:,:) 
!
!
! 2.1 Open the output file
!-------------------------
!
   CLNAME='out_level2_'//CNAMEID//'_'//cfreq//'_'//cangle//'.nc'
   status = NF90_CREATE(CLNAME,NF90_CLOBBER,ncid)
   IF(status /= nf90_noerr) CALL ABOR1(trim(nf90_strerror(status)))
!
!
! 2.2 Define file dimensions
!---------------------------
!
   status = NF90_DEF_DIM(ncid,CCLON_NAME,NLONS,londimid)
   IF(status /= nf90_noerr) CALL ABOR1(trim(nf90_strerror(status)))
!
   status = NF90_DEF_DIM(ncid,CCLAT_NAME,NLATS,latdimid)
   IF(status /= nf90_noerr) CALL ABOR1(trim(nf90_strerror(status)))
!
   status = NF90_DEF_DIM(ncid,CCLVL_NAME,NLVLS,lvldimid)
   IF(status /= nf90_noerr) CALL ABOR1(trim(nf90_strerror(status)))
!
!status = NF90_DEF_DIM(ncid, CCTIME_NAME, NF90_UNLIMITED, timedimid)
   IF(status /= nf90_noerr) CALL ABOR1(trim(nf90_strerror(status)))
!
   status = NF90_DEF_DIM(ncid, CCTIME_NAME, NTIMES, timedimid)
   IF(status /= nf90_noerr) CALL ABOR1(trim(nf90_strerror(status)))
!
!
   status = NF90_DEF_VAR(ncid, CCLON_NAME, NF90_REAL, londimid, nlon_varid)
   IF(status /= nf90_noerr) CALL ABOR1(trim(nf90_strerror(status)))
!
   status = NF90_DEF_VAR(ncid, CCLAT_NAME, NF90_REAL, latdimid, nlat_varid)
   IF(status /= nf90_noerr) CALL ABOR1(trim(nf90_strerror(status)))
!
   status = NF90_DEF_VAR(ncid, CCLVL_NAME, NF90_REAL, lvldimid, nlvl_varid)
   IF(status /= nf90_noerr) CALL ABOR1(trim(nf90_strerror(status)))
!
   status = NF90_DEF_VAR(ncid, CCTIME_NAME, NF90_REAL, timedimid, ntime_varid)
   IF(status /= nf90_noerr) CALL ABOR1(trim(nf90_strerror(status)))
!
!
! Assign units attributes to coordinate variables.
!
   status = nf90_put_att(ncid, nlat_varid, CCUNITS, CCLAT_UNITS)
   IF(status /= nf90_noerr) CALL ABOR1(trim(nf90_strerror(status)))
!
   status = nf90_put_att(ncid, nlon_varid, CCUNITS, CCLON_UNITS)
   IF(status /= nf90_noerr) CALL ABOR1(trim(nf90_strerror(status)))
!
   status = nf90_put_att(ncid, ntime_varid, CCUNITS, CCTIME_UNITS)
   IF(status /= nf90_noerr) CALL ABOR1(trim(nf90_strerror(status)))
!
   dimids = (/ londimid, latdimid, lvldimid, timedimid /)
!
!
!
! 2.3 Define the netCDF variables for the output data
!----------------------------------------------------
!
   status = nf90_def_var(ncid, CLFTAUVEGH_NAME, NF90_REAL, dimids, iftauvegh_varid) 
   IF(status /= nf90_noerr) CALL ABOR1(trim(nf90_strerror(status)))
!
   status = nf90_def_var(ncid, CLFTAUVEGV_NAME, NF90_REAL, dimids, iftauvegv_varid) 
   IF(status /= nf90_noerr) CALL ABOR1(trim(nf90_strerror(status)))
!
   status = nf90_def_var(ncid, CLFBAR_NAME, NF90_REAL, dimids, ifbar_varid) 
   IF(status /= nf90_noerr) CALL ABOR1(trim(nf90_strerror(status)))
!
   status = nf90_def_var(ncid, CLVWC_NAME, NF90_REAL, dimids, ivwc_varid) 
   IF(status /= nf90_noerr) CALL ABOR1(trim(nf90_strerror(status)))
!
   status = nf90_def_var(ncid, CLTAUATM_NAME, NF90_REAL, dimids, itauatm_varid) 
   IF(status /= nf90_noerr) CALL ABOR1(trim(nf90_strerror(status)))
!
   status = nf90_def_var(ncid, CLTBATMUP_NAME, NF90_REAL, dimids, itbatmup_varid) 
   IF(status /= nf90_noerr) CALL ABOR1(trim(nf90_strerror(status)))
!
! Assign units attributes to the netCDF variables.
!
   status = nf90_put_att(ncid, iftauvegv_varid, CCUNITS, CCPFILED_UNITS)
   IF(status /= nf90_noerr) CALL ABOR1(trim(nf90_strerror(status)))
!
   status = nf90_put_att(ncid, iftauvegh_varid, CCUNITS, CCPFILED_UNITS)
   IF(status /= nf90_noerr) CALL ABOR1(trim(nf90_strerror(status)))
!
   status = nf90_put_att(ncid, ifbar_varid, CCUNITS, CCPFILED_UNITS)
   IF(status /= nf90_noerr) CALL ABOR1(trim(nf90_strerror(status)))
!
   status = nf90_put_att(ncid, ivwc_varid, CCUNITS, CCVWC_UNITS)
   IF(status /= nf90_noerr) CALL ABOR1(trim(nf90_strerror(status)))
!
   status = nf90_put_att(ncid, itauatm_varid, CCUNITS, CCPFILED_UNITS)
   IF(status /= nf90_noerr) CALL ABOR1(trim(nf90_strerror(status)))
!
   status = nf90_put_att(ncid, itbatmup_varid, CCUNITS, CCTEMP_UNITS)
   IF(status /= nf90_noerr) CALL ABOR1(trim(nf90_strerror(status)))
!
   status =  NF90_ENDDEF(ncid)
   IF(status /= nf90_noerr) CALL ABOR1(trim(nf90_strerror(status)))

! 2.5 Write the coordinate variable data (latitudes and longitudes) into the netCDF file
!---------------------------------------------------------------------------------------

   status =  NF90_PUT_VAR(ncid, nlat_varid, xlats) 
   IF(status /= nf90_noerr) CALL ABOR1(trim(nf90_strerror(status)))
   
   status =  NF90_PUT_VAR(ncid, nlon_varid, xlons) 
   IF(status /= nf90_noerr) CALL ABOR1(trim(nf90_strerror(status)))
  
   status =  NF90_PUT_VAR(ncid, ntime_varid, xtimes) 
   IF(status /= nf90_noerr) CALL ABOR1(trim(nf90_strerror(status)))
  

! 2.6 Write the data 
!-------------------
!
   status =  NF90_PUT_VAR(ncid, iftauvegh_varid, ftauvegh_out) 
   IF(status /= nf90_noerr) CALL ABOR1(trim(nf90_strerror(status)))
!
   status =  NF90_PUT_VAR(ncid, iftauvegv_varid, ftauvegv_out) 
   IF(status /= nf90_noerr) CALL ABOR1(trim(nf90_strerror(status)))
!
   status =  NF90_PUT_VAR(ncid, ifbar_varid, fbar_out) 
   IF(status /= nf90_noerr) CALL ABOR1(trim(nf90_strerror(status)))
!
   status =  NF90_PUT_VAR(ncid, ivwc_varid, vwc_out) 
   IF(status /= nf90_noerr) CALL ABOR1(trim(nf90_strerror(status)))
!
   status =  NF90_PUT_VAR(ncid, itauatm_varid, tauatm_out) 
   IF(status /= nf90_noerr) CALL ABOR1(trim(nf90_strerror(status)))
!
   status =  NF90_PUT_VAR(ncid, itbatmup_varid, tbatmup_out) 
   IF(status /= nf90_noerr) CALL ABOR1(trim(nf90_strerror(status)))
!
! 2.7 Close the file
!-------------------
!
   status = nf90_close(ncid) 
   IF(status /= nf90_noerr) CALL ABOR1(trim(nf90_strerror(status)))
!
!
!

 
END SELECT


! LEVEL 3 outputs: third netcdf file
!-----------------------------------

SELECT CASE (JPHISTLEV)

 CASE ( 3_JPIM )
 
!
! 3.0 put variable in right format 
!---------------------------------
!
!  BARNOSNOW 
   var_in_tr = reshape(ftfrac(:,1), (/NLONS,NLATS,NTIMES/))
   barenosnow_out(:,:,1,:) = var_in_tr(:,:,:) 
!
!  LOWNOSNOW 
   var_in_tr = reshape(ftfrac(:,3), (/NLONS,NLATS,NTIMES/))
   lownosnow_out(:,:,1,:) = var_in_tr(:,:,:) 
!
!  VWCL
   var_in_tr = reshape(fwc_veg(:,1), (/NLONS,NLATS,NTIMES/))
   vwcl_out(:,:,1,:) = var_in_tr(:,:,:) 
!
!  VWCH
   var_in_tr = reshape(fwc_veg(:,2), (/NLONS,NLATS,NTIMES/))
   vwch_out(:,:,1,:) = var_in_tr(:,:,:) 
!
!  BPARAML 
   var_in_tr = reshape(fb(:,1), (/NLONS,NLATS,NTIMES/))
   bparaml_out(:,:,1,:) = var_in_tr(:,:,:) 
!
!  BPARAMH 
   var_in_tr = reshape(fb(:,2), (/NLONS,NLATS,NTIMES/))
   bparamh_out(:,:,1,:) = var_in_tr(:,:,:) 
!
!  RUGOH 
   var_in_tr = reshape(fh(:), (/NLONS,NLATS,NTIMES/))
   rugoh_out(:,:,1,:) = var_in_tr(:,:,:) 
!
!  WILTSM 
   var_in_tr = reshape(fWP(:), (/NLONS,NLATS,NTIMES/))
   wiltsm_out(:,:,1,:) = var_in_tr(:,:,:) 
!
!  EH
   var_in_tr = reshape(fsurf_emis(:,1), (/NLONS,NLATS,NTIMES/))
   eh_out(:,:,1,:) = var_in_tr(:,:,:) 
!
!  EV
   var_in_tr = reshape(fsurf_emis(:,2), (/NLONS,NLATS,NTIMES/))
   ev_out(:,:,1,:) = var_in_tr(:,:,:) 
!
!  CPARAM 
   var_in_tr = reshape(fteffC(:,1), (/NLONS,NLATS,NTIMES/))
   cparam_out(:,:,1,:) = var_in_tr(:,:,:) 
!
!
! 3.1 Open the output file
!-------------------------
!
   CLNAME='out_level3_'//CNAMEID//'_'//cfreq//'_'//cangle//'.nc'
   status = NF90_CREATE(CLNAME,NF90_CLOBBER,ncid)
   IF(status /= nf90_noerr) CALL ABOR1(trim(nf90_strerror(status)))
!
!
! 3.2 Define file dimensions
!---------------------------
!
   status = NF90_DEF_DIM(ncid,CCLON_NAME,NLONS,londimid)
   IF(status /= nf90_noerr) CALL ABOR1(trim(nf90_strerror(status)))
!
   status = NF90_DEF_DIM(ncid,CCLAT_NAME,NLATS,latdimid)
   IF(status /= nf90_noerr) CALL ABOR1(trim(nf90_strerror(status)))
!
   status = NF90_DEF_DIM(ncid,CCLVL_NAME,NLVLS,lvldimid)
   IF(status /= nf90_noerr) CALL ABOR1(trim(nf90_strerror(status)))
!
!status = NF90_DEF_DIM(ncid, CCTIME_NAME, NF90_UNLIMITED, timedimid)
   IF(status /= nf90_noerr) CALL ABOR1(trim(nf90_strerror(status)))
!
   status = NF90_DEF_DIM(ncid, CCTIME_NAME, NTIMES, timedimid)
   IF(status /= nf90_noerr) CALL ABOR1(trim(nf90_strerror(status)))
!
!
   status = NF90_DEF_VAR(ncid, CCLON_NAME, NF90_REAL, londimid, nlon_varid)
   IF(status /= nf90_noerr) CALL ABOR1(trim(nf90_strerror(status)))
!
   status = NF90_DEF_VAR(ncid, CCLAT_NAME, NF90_REAL, latdimid, nlat_varid)
   IF(status /= nf90_noerr) CALL ABOR1(trim(nf90_strerror(status)))
!
   status = NF90_DEF_VAR(ncid, CCLVL_NAME, NF90_REAL, lvldimid, nlvl_varid)
   IF(status /= nf90_noerr) CALL ABOR1(trim(nf90_strerror(status)))
!
   status = NF90_DEF_VAR(ncid, CCTIME_NAME, NF90_REAL, timedimid, ntime_varid)
   IF(status /= nf90_noerr) CALL ABOR1(trim(nf90_strerror(status)))
!
!
! Assign units attributes to coordinate variables.
!
   status = nf90_put_att(ncid, nlat_varid, CCUNITS, CCLAT_UNITS)
   IF(status /= nf90_noerr) CALL ABOR1(trim(nf90_strerror(status)))
!
   status = nf90_put_att(ncid, nlon_varid, CCUNITS, CCLON_UNITS)
   IF(status /= nf90_noerr) CALL ABOR1(trim(nf90_strerror(status)))
!
   status = nf90_put_att(ncid, ntime_varid, CCUNITS, CCTIME_UNITS)
   IF(status /= nf90_noerr) CALL ABOR1(trim(nf90_strerror(status)))
!
   dimids = (/ londimid, latdimid, lvldimid, timedimid /)
!
!
!
! 3.3 Define the netCDF variables for the output data
!----------------------------------------------------
!
   status = nf90_def_var(ncid, CLBARENOSNOW_NAME, NF90_REAL, dimids, ibarenosnow_varid) 
   IF(status /= nf90_noerr) CALL ABOR1(trim(nf90_strerror(status)))
!
   status = nf90_def_var(ncid, CLLOWNOSNOW_NAME, NF90_REAL, dimids, ilownosnow_varid) 
   IF(status /= nf90_noerr) CALL ABOR1(trim(nf90_strerror(status)))
!
   status = nf90_def_var(ncid, CLVWCL_NAME, NF90_REAL, dimids, ivwcl_varid) 
   IF(status /= nf90_noerr) CALL ABOR1(trim(nf90_strerror(status)))
!
   status = nf90_def_var(ncid, CLVWCH_NAME, NF90_REAL, dimids, ivwch_varid) 
   IF(status /= nf90_noerr) CALL ABOR1(trim(nf90_strerror(status)))
!
   status = nf90_def_var(ncid, CLBPARAML_NAME, NF90_REAL, dimids, ibparaml_varid) 
   IF(status /= nf90_noerr) CALL ABOR1(trim(nf90_strerror(status)))
!
   status = nf90_def_var(ncid, CLBPARAMH_NAME, NF90_REAL, dimids, ibparamh_varid) 
   IF(status /= nf90_noerr) CALL ABOR1(trim(nf90_strerror(status)))
!
   status = nf90_def_var(ncid, CLRUGOH_NAME, NF90_REAL, dimids, irugoh_varid) 
   IF(status /= nf90_noerr) CALL ABOR1(trim(nf90_strerror(status)))
!
   status = nf90_def_var(ncid, CLWILTSM_NAME, NF90_REAL, dimids, iwiltsm_varid) 
   IF(status /= nf90_noerr) CALL ABOR1(trim(nf90_strerror(status)))
!
   status = nf90_def_var(ncid, CLEH_NAME, NF90_REAL, dimids, ieh_varid) 
   IF(status /= nf90_noerr) CALL ABOR1(trim(nf90_strerror(status)))
!
   status = nf90_def_var(ncid, CLEV_NAME, NF90_REAL, dimids, iev_varid) 
   IF(status /= nf90_noerr) CALL ABOR1(trim(nf90_strerror(status)))
!
   status = nf90_def_var(ncid, CLCPARAM_NAME, NF90_REAL, dimids, icparam_varid) 
   IF(status /= nf90_noerr) CALL ABOR1(trim(nf90_strerror(status)))
!
! Assign units attributes to the netCDF variables.
!
   status = nf90_put_att(ncid, ibarenosnow_varid, CCUNITS, CCPFILED_UNITS)
   IF(status /= nf90_noerr) CALL ABOR1(trim(nf90_strerror(status)))
!
   status = nf90_put_att(ncid, ilownosnow_varid, CCUNITS, CCPFILED_UNITS)
   IF(status /= nf90_noerr) CALL ABOR1(trim(nf90_strerror(status)))
!
   status = nf90_put_att(ncid, ivwcl_varid, CCUNITS, CCVWC_UNITS)
   IF(status /= nf90_noerr) CALL ABOR1(trim(nf90_strerror(status)))
!
   status = nf90_put_att(ncid, ivwch_varid, CCUNITS, CCVWC_UNITS)
   IF(status /= nf90_noerr) CALL ABOR1(trim(nf90_strerror(status)))
!
   status = nf90_put_att(ncid, ibparaml_varid, CCUNITS, CCPFILED_UNITS)
   IF(status /= nf90_noerr) CALL ABOR1(trim(nf90_strerror(status)))
!
   status = nf90_put_att(ncid, ibparamh_varid, CCUNITS, CCPFILED_UNITS)
   IF(status /= nf90_noerr) CALL ABOR1(trim(nf90_strerror(status)))
!
   status = nf90_put_att(ncid, irugoh_varid, CCUNITS, CCPFILED_UNITS)
   IF(status /= nf90_noerr) CALL ABOR1(trim(nf90_strerror(status)))
!
   status = nf90_put_att(ncid, iwiltsm_varid, CCUNITS, CCPFILED_UNITS)
   IF(status /= nf90_noerr) CALL ABOR1(trim(nf90_strerror(status)))
!
   status = nf90_put_att(ncid, ieh_varid, CCUNITS, CCPFILED_UNITS)
   IF(status /= nf90_noerr) CALL ABOR1(trim(nf90_strerror(status)))
!
   status = nf90_put_att(ncid, iev_varid, CCUNITS, CCPFILED_UNITS)
   IF(status /= nf90_noerr) CALL ABOR1(trim(nf90_strerror(status)))
!
   status = nf90_put_att(ncid, icparam_varid, CCUNITS, CCTEMP_UNITS)
   IF(status /= nf90_noerr) CALL ABOR1(trim(nf90_strerror(status)))

!
!
! 3.4 End define mode
!--------------------


   status =  NF90_ENDDEF(ncid)
   IF(status /= nf90_noerr) CALL ABOR1(trim(nf90_strerror(status)))

! 3.5 Write the coordinate variable data (latitudes and longitudes) into the netCDF file
!---------------------------------------------------------------------------------------

   status =  NF90_PUT_VAR(ncid, nlat_varid, xlats) 
   IF(status /= nf90_noerr) CALL ABOR1(trim(nf90_strerror(status)))
   
   status =  NF90_PUT_VAR(ncid, nlon_varid, xlons) 
   IF(status /= nf90_noerr) CALL ABOR1(trim(nf90_strerror(status)))
  
   status =  NF90_PUT_VAR(ncid, ntime_varid, xtimes) 
   IF(status /= nf90_noerr) CALL ABOR1(trim(nf90_strerror(status)))
  

! 3.6 Write the data 
!-------------------
!
   status =  NF90_PUT_VAR(ncid, ibarenosnow_varid, barenosnow_out) 
   IF(status /= nf90_noerr) CALL ABOR1(trim(nf90_strerror(status)))
!
   status =  NF90_PUT_VAR(ncid, ilownosnow_varid, lownosnow_out) 
   IF(status /= nf90_noerr) CALL ABOR1(trim(nf90_strerror(status)))
!
   status =  NF90_PUT_VAR(ncid, ivwcl_varid, vwcl_out) 
   IF(status /= nf90_noerr) CALL ABOR1(trim(nf90_strerror(status)))
!
   status =  NF90_PUT_VAR(ncid, ivwch_varid, vwch_out) 
   IF(status /= nf90_noerr) CALL ABOR1(trim(nf90_strerror(status)))
!
   status =  NF90_PUT_VAR(ncid, ibparaml_varid, bparaml_out) 
   IF(status /= nf90_noerr) CALL ABOR1(trim(nf90_strerror(status)))
!
   status =  NF90_PUT_VAR(ncid, ibparamh_varid, bparamh_out) 
   IF(status /= nf90_noerr) CALL ABOR1(trim(nf90_strerror(status)))
!
   status =  NF90_PUT_VAR(ncid, irugoh_varid, rugoh_out) 
   IF(status /= nf90_noerr) CALL ABOR1(trim(nf90_strerror(status)))
!
   status =  NF90_PUT_VAR(ncid, iwiltsm_varid, wiltsm_out) 
   IF(status /= nf90_noerr) CALL ABOR1(trim(nf90_strerror(status)))
!
   status =  NF90_PUT_VAR(ncid, ieh_varid, eh_out) 
   IF(status /= nf90_noerr) CALL ABOR1(trim(nf90_strerror(status)))
!
   status =  NF90_PUT_VAR(ncid, iev_varid, ev_out) 
   IF(status /= nf90_noerr) CALL ABOR1(trim(nf90_strerror(status)))
!
   status =  NF90_PUT_VAR(ncid, icparam_varid, cparam_out) 
   IF(status /= nf90_noerr) CALL ABOR1(trim(nf90_strerror(status)))
!
! 3.7 Close the file
!-------------------
!
   status = nf90_close(ncid) 
   IF(status /= nf90_noerr) CALL ABOR1(trim(nf90_strerror(status)))
!
!
!

 
END SELECT

       
END SUBROUTINE WRCMEMNETCDF

