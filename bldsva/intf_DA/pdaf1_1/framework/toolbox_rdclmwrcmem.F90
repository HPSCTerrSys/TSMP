module rdclm_wrcmem 

  use netcdf
  use clm4cmem, only: CLM_DATA, SHOWINFO
  use parkind1, only: JPIM, JPRM
  USE YOMCMEMNETCDF, ONLY : NTIMES_SM,NLATS_SM, NLONS_SM, NLVLS_SM,NTIMES, xlats,xlons,xlvls,xtimes
contains

  ! ===============================================================
  ! SUBROUTINE to storage the output SATELLITE synthetic data
  ! as NetCDF files in the cmem directory.
  ! ---------------------------------------------------------------
  ! (c) 2017 Pablo Saavedra G. (email: pablosaa@uni-bonn.de)
  ! UNIVERSITY OF BONN, GERMANY
  ! See LICENCE
  ! ---------------------------------------------------------------
  subroutine write_satellite_operator(SAT,CLMVARS,OUT_fname,nc_out)
    use clm4cmem, only : SATELLITE,CLM_DATA
    implicit none
    ! INPUT VARIABLES:
    type(SATELLITE), intent(in)  :: SAT
    type(CLM_DATA),  intent(in)  :: CLMVARS
    character(len=*), intent(in) :: OUT_fname
    character(len=*), optional, intent(in) :: nc_out

    ! Local VARIABLES
    character (len=:), allocatable :: nc_path_out
    integer :: I, J, idxstr ! JJ, N
    integer :: pix_dimid,lon_varid,lat_varid,time_dimid,inc_dimid,lev_dimid
    integer :: glon_varid,glat_varid,lat_dimid,lon_dimid
    integer :: wrncid, inc_varid, t_varid, tbh_varid,&
                  tbv_varid,gtbh_varid,gtbv_varid
    integer :: gsm_varid,gst_varid,gsd_varid,gRSN_varid,gtskin_varid, &
                  ginci_varid,gz_varid, gsand_varid,gclay_varid,gcvh_varid,&
                  gcvl_varid,gtair_varid,gtvh_varid,gtvl_varid,glsm_varid,&
                  gwater_varid,glai_varid,gslope_varid,gaspect_varid
    integer :: NPIX, NLONS, NLATS, NTIME, NINC,NLEV
    integer :: tb_dimid(4),SAT_tb_dimid(3),&
                  f4D_dimid(4),f3D_dimid(3),f2D_dimid(2)
    integer :: datum_uhrzeit(8)
    character (len=21) ::DATESTRING
    character (len = *), parameter :: UNITS = "units"
    character (len = *), parameter :: LONG_NAME = "long_name"
    character (len = *), parameter :: SHORT_NAME = "short_name"
    character (len = *), parameter :: MISS_VAL = "missing_value"
    character (len = *), parameter :: FILL_VAL = "_FillValue"
    character (len = *), parameter :: VALID_RA = "valid_range"
    character (len =100) :: FILETMP

    ! Populating variables:
    NPIX  = size(SAT%lon_foprt) ! should be same with lat_foprt
    NINC  = size(SAT%theta)
    NTIME = size(SAT%time)
    NLONS = size(SAT%lon)
    NLATS = size(SAT%lat)
    NLEV  = size(CLMVARS%levsoi)
    
    IF(.not.present(nc_out)) THEN
       nc_path_out = './'
    ELSE
       nc_path_out = trim(nc_out)
    END IF
    call date_and_time(VALUES=datum_uhrzeit)
    write(DATESTRING,'(I2.2"."I2.2"."I4 "T" I2.2":"I2.2":"I2.2)') datum_uhrzeit((/3,2,1,5,6,7/))
    ! 1-RD: Geopotential at surface [km]
    ! -----------------------------------------------------------------
    ! Create the file.
    if (SHOWINFO) print*, "Creating satellite_output NetCDF"
      
    call check( nf90_create(nc_path_out//trim(OUT_fname), nf90_clobber, wrncid))
    ! Define the dimensions.
    call check( nf90_def_dim(wrncid, "NPIXEL", NPIX, pix_dimid))
    call check( nf90_def_dim(wrncid, "GRIDLONG", NLONS, lon_dimid))
    call check( nf90_def_dim(wrncid, "GRIDLAT", NLATS, lat_dimid))
    call check( nf90_def_dim(wrncid, "INC_ANGLE", NINC, inc_dimid))
    call check( nf90_def_dim(wrncid, "TIME", NTIME, time_dimid))
    call check( nf90_def_dim(wrncid, "LEVEL", NLEV, lev_dimid))

    ! Define the variable
    call check( nf90_def_var(wrncid, "LONGITUDE", NF90_REAL,pix_dimid,lon_varid))
    call check( nf90_put_att(wrncid, lon_varid, LONG_NAME, "Longitude"))
    call check( nf90_put_att(wrncid, lon_varid, SHORT_NAME, "lon"))
    call check( nf90_put_att(wrncid, lon_varid, UNITS, "degrees_east"))

    call check( nf90_def_var(wrncid, "LATITUDE", NF90_REAL,pix_dimid,lat_varid))
    call check( nf90_put_att(wrncid, lat_varid, LONG_NAME, "Latitude"))
    call check( nf90_put_att(wrncid, lat_varid, SHORT_NAME, "lat"))
    call check( nf90_put_att(wrncid, lat_varid, UNITS, "degrees_north"))

    call check( nf90_def_var(wrncid, "INC_ANGLE", NF90_REAL,inc_dimid,inc_varid))
    call check( nf90_put_att(wrncid, inc_varid, LONG_NAME, "Incidence angle"))
    call check( nf90_put_att(wrncid, inc_varid, SHORT_NAME, "INC"))
    call check( nf90_put_att(wrncid, inc_varid, UNITS, "degrees_nadir"))

    call check( nf90_def_var(wrncid, "TIME", NF90_REAL,time_dimid,t_varid))
    call check( nf90_put_att(wrncid, t_varid, LONG_NAME, "time counter"))
    call check( nf90_put_att(wrncid, t_varid, SHORT_NAME, "time"))
    call check( nf90_put_att(wrncid, t_varid, UNITS, "hour"))

    ! * Satellite Brightness temperatures H-pol
    SAT_tb_dimid = (/ pix_dimid, inc_dimid, time_dimid /)
    call check( nf90_def_var(wrncid, "TBSAT_H", NF90_REAL,SAT_tb_dimid,tbh_varid))
    call check( nf90_put_att(wrncid, tbh_varid, LONG_NAME, "Brightness Temperature H-pol"))
    call check( nf90_put_att(wrncid, tbh_varid, SHORT_NAME, "TB_H"))
    call check( nf90_put_att(wrncid, tbh_varid, UNITS, "K"))
    call check( nf90_put_att(wrncid, tbh_varid, MISS_VAL, -99.))
    !call check( nf90_put_att(wrncid, tbh_varid, VALID_RA, (/0.,400./)))
    !call check( nf90_put_att(wrncid, tbh_varid, FILL_VAL, NF90_FILL_REAL))
    ! * Satellite Brightness temperatures V-pol
    call check( nf90_def_var(wrncid, "TBSAT_V", NF90_REAL,SAT_tb_dimid,tbv_varid))
    call check( nf90_put_att(wrncid, tbv_varid, LONG_NAME, "Brightness Temperature V-pol"))
    call check( nf90_put_att(wrncid, tbv_varid, SHORT_NAME, "TB_V"))
    call check( nf90_put_att(wrncid, tbv_varid, UNITS, "K"))
    call check( nf90_put_att(wrncid, tbv_varid, MISS_VAL, -99.))
    !call check( nf90_put_att(wrncid, tbv_varid, VALID_RA, (/0.,400./)))
    !call check( nf90_put_att(wrncid, tbv_varid, FILL_VAL, NF90_FILL_REAL))

    ! * Grid Satellite Brightness temperatures H-pol
    tb_dimid = (/ lon_dimid,lat_dimid,inc_dimid, time_dimid /)
    call check( nf90_def_var(wrncid, "TB_H",NF90_REAL,tb_dimid,gtbh_varid))
    call check( nf90_put_att(wrncid, gtbh_varid, LONG_NAME, "Grid Brightness Temperature H-pol"))
    call check( nf90_put_att(wrncid, gtbh_varid, SHORT_NAME, "gTB_H"))
    call check( nf90_put_att(wrncid, gtbh_varid, UNITS, "K"))
    call check( nf90_put_att(wrncid, gtbh_varid, MISS_VAL, -99.))
    ! * Grid Satellite Brightness temperatures V-pol
    call check( nf90_def_var(wrncid, "TB_V",NF90_REAL,tb_dimid,gtbv_varid))
    call check( nf90_put_att(wrncid, gtbv_varid, LONG_NAME, "Grid Brightness Temperature V-pol"))
    call check( nf90_put_att(wrncid, gtbv_varid, SHORT_NAME, "gTB_V"))
    call check( nf90_put_att(wrncid, gtbv_varid, UNITS, "K"))
    call check( nf90_put_att(wrncid, gtbv_varid, MISS_VAL, -99.))
    
    f4D_dimid = (/ lon_dimid,lat_dimid,lev_dimid,time_dimid /)
    f3D_dimid = (/ lon_dimid,lat_dimid,time_dimid /)
    f2D_dimid = (/ lon_dimid,lat_dimid /)
    ! * Grid longitude
    call check( nf90_def_var(wrncid, "GRIDSLONG",NF90_REAL,f2D_dimid,glon_varid))
    call check( nf90_put_att(wrncid, glon_varid, LONG_NAME, "Grids Longitude"))
    call check( nf90_put_att(wrncid, glon_varid, SHORT_NAME, "glon"))
    call check( nf90_put_att(wrncid, glon_varid, UNITS, "degrees_east"))
    ! * Grid latitude
    call check( nf90_def_var(wrncid, "GRIDSLAT", NF90_REAL,f2D_dimid,glat_varid))
    call check( nf90_put_att(wrncid, glat_varid, LONG_NAME, "Grids Latitude"))
    call check( nf90_put_att(wrncid, glat_varid, SHORT_NAME, "glat"))
    call check( nf90_put_att(wrncid, glat_varid, UNITS, "degrees_north"))
    ! * Grid soil moisture
    call check( nf90_def_var(wrncid, "H2OSOI",NF90_REAL,f4D_dimid,gsm_varid))
    call check( nf90_put_att(wrncid, gsm_varid, LONG_NAME, "soil moisture"))
    call check( nf90_put_att(wrncid, gsm_varid, SHORT_NAME, "SWVL"))
    call check( nf90_put_att(wrncid, gsm_varid, UNITS, "%/%"))
    call check( nf90_put_att(wrncid, gsm_varid, MISS_VAL, -99.))
    ! * Grid soil temperature
    call check( nf90_def_var(wrncid, "TSOI",NF90_REAL,f4D_dimid,gst_varid))
    call check( nf90_put_att(wrncid, gst_varid, LONG_NAME, "soil temperature"))
    call check( nf90_put_att(wrncid, gst_varid, SHORT_NAME, "STL"))
    call check( nf90_put_att(wrncid, gst_varid, UNITS, "K"))
    call check( nf90_put_att(wrncid, gst_varid, MISS_VAL, -99.))
    ! * Grid snow depth
    call check( nf90_def_var(wrncid, "SNOWDEPTH",NF90_REAL,f3D_dimid,gsd_varid))
    call check( nf90_put_att(wrncid, gsd_varid, LONG_NAME, "snow depth"))
    call check( nf90_put_att(wrncid, gsd_varid, SHORT_NAME, "SD"))
    call check( nf90_put_att(wrncid, gsd_varid, UNITS, "m"))
    call check( nf90_put_att(wrncid, gsd_varid, MISS_VAL, -99.))
    ! * Grid RSN
    call check( nf90_def_var(wrncid, "RSN",NF90_REAL,f3D_dimid,gRSN_varid))
    call check( nf90_put_att(wrncid, gRSN_varid, LONG_NAME, "RSN"))
    call check( nf90_put_att(wrncid, gRSN_varid, SHORT_NAME, "RSN"))
    call check( nf90_put_att(wrncid, gRSN_varid, UNITS, "unitless"))
    call check( nf90_put_att(wrncid, gRSN_varid, MISS_VAL, -99.))
    ! * Grid skin temperature
    call check( nf90_def_var(wrncid, "TSKIN",NF90_REAL,f3D_dimid,gtskin_varid))
    call check( nf90_put_att(wrncid, gtskin_varid, LONG_NAME, "soil skin temperature"))
    call check( nf90_put_att(wrncid, gtskin_varid, SHORT_NAME, "TSKIN"))
    call check( nf90_put_att(wrncid, gtskin_varid, UNITS, "K"))
    call check( nf90_put_att(wrncid, gtskin_varid, MISS_VAL, -99.))
    ! * Grid clay fraction
    call check( nf90_def_var(wrncid, "CLAY",NF90_REAL,f3D_dimid,gclay_varid))
    call check( nf90_put_att(wrncid, gclay_varid, LONG_NAME, "clay fraction"))
    call check( nf90_put_att(wrncid, gclay_varid, SHORT_NAME, "CLAY"))
    call check( nf90_put_att(wrncid, gclay_varid, UNITS, "%"))
    call check( nf90_put_att(wrncid, gclay_varid, MISS_VAL, -99.))
    ! * Grid sand fraction
    call check( nf90_def_var(wrncid, "SAND",NF90_REAL,f3D_dimid,gsand_varid))
    call check( nf90_put_att(wrncid, gsand_varid, LONG_NAME, "sand fraction"))
    call check( nf90_put_att(wrncid, gsand_varid, SHORT_NAME, "SAND"))
    call check( nf90_put_att(wrncid, gsand_varid, UNITS, "%"))
    call check( nf90_put_att(wrncid, gsand_varid, MISS_VAL, -99.))
    ! * Grid coverage of high vegetation
    call check( nf90_def_var(wrncid, "CVH",NF90_REAL,f3D_dimid,gcvh_varid))
    call check( nf90_put_att(wrncid, gcvh_varid, LONG_NAME, "coverage of high vegetation"))
    call check( nf90_put_att(wrncid, gcvh_varid, SHORT_NAME, "CVH"))
    call check( nf90_put_att(wrncid, gcvh_varid, UNITS, "binary"))
    call check( nf90_put_att(wrncid, gcvh_varid, MISS_VAL, -99.))
    ! * Grid coverage of low vegetation
    call check( nf90_def_var(wrncid, "CVL",NF90_REAL,f3D_dimid,gcvl_varid))
    call check( nf90_put_att(wrncid, gcvl_varid, LONG_NAME, "coverage of low vegetation"))
    call check( nf90_put_att(wrncid, gcvl_varid, SHORT_NAME, "CVL"))
    call check( nf90_put_att(wrncid, gcvl_varid, UNITS, "binary"))
    call check( nf90_put_att(wrncid, gcvl_varid, MISS_VAL, -99.))
    ! * Grid type of high vegetation
    call check( nf90_def_var(wrncid, "TVH",NF90_REAL,f3D_dimid,gtvh_varid))
    call check( nf90_put_att(wrncid, gtvh_varid, LONG_NAME, "type of high vegetation"))
    call check( nf90_put_att(wrncid, gtvh_varid, SHORT_NAME, "TVH"))
    call check( nf90_put_att(wrncid, gtvh_varid, UNITS, "index"))
    call check( nf90_put_att(wrncid, gtvh_varid, MISS_VAL, -99.))
    ! * Grid type of low vegetation
    call check( nf90_def_var(wrncid, "TVL",NF90_REAL,f3D_dimid,gtvl_varid))
    call check( nf90_put_att(wrncid, gtvl_varid, LONG_NAME, "type of low vegetation"))
    call check( nf90_put_att(wrncid, gtvl_varid, SHORT_NAME, "TVL"))
    call check( nf90_put_att(wrncid, gtvl_varid, UNITS, "index"))
    call check( nf90_put_att(wrncid, gtvl_varid, MISS_VAL, -99.))
    ! * Grid land coverage
    call check( nf90_def_var(wrncid, "LSM",NF90_REAL,f3D_dimid,glsm_varid))
    call check( nf90_put_att(wrncid, glsm_varid, LONG_NAME, "land coverage"))
    call check( nf90_put_att(wrncid, glsm_varid, SHORT_NAME, "LSM"))
    call check( nf90_put_att(wrncid, glsm_varid, UNITS, "binary"))
    call check( nf90_put_att(wrncid, glsm_varid, MISS_VAL, -99.))
    ! * Grid water coverage
    call check( nf90_def_var(wrncid, "WATER",NF90_REAL,f3D_dimid,gwater_varid))
    call check( nf90_put_att(wrncid, gwater_varid, LONG_NAME, "water coverage"))
    call check( nf90_put_att(wrncid, gwater_varid, SHORT_NAME, "WATER"))
    call check( nf90_put_att(wrncid, gwater_varid, UNITS, "binary"))
    call check( nf90_put_att(wrncid, gwater_varid, MISS_VAL, -99.))
    ! * Grid leaf area index
    call check( nf90_def_var(wrncid, "LAI",NF90_REAL,f3D_dimid,glai_varid))
    call check( nf90_put_att(wrncid, glai_varid, LONG_NAME, "leaf area index"))
    call check( nf90_put_att(wrncid, glai_varid, SHORT_NAME, "LAI"))
    call check( nf90_put_att(wrncid, glai_varid, UNITS, "index"))
    call check( nf90_put_att(wrncid, glai_varid, MISS_VAL, -99.))
    ! * Grid 2 meter air temperature
    call check( nf90_def_var(wrncid, "TAIR",NF90_REAL,f3D_dimid,gtair_varid))
    call check( nf90_put_att(wrncid, gtair_varid, LONG_NAME, "2 meter air temperature"))
    call check( nf90_put_att(wrncid, gtair_varid, SHORT_NAME, "TAIR"))
    call check( nf90_put_att(wrncid, gtair_varid, UNITS, "K"))
    call check( nf90_put_att(wrncid, gtair_varid, MISS_VAL, -99.))
    ! * Grid incidence angle
    call check( nf90_def_var(wrncid, "INCIDENCE",NF90_REAL,f3D_dimid,ginci_varid))
    call check( nf90_put_att(wrncid, ginci_varid, LONG_NAME, "incidence angel"))
    call check( nf90_put_att(wrncid, ginci_varid, SHORT_NAME, "INCI"))
    call check( nf90_put_att(wrncid, ginci_varid, UNITS, "degree"))
    call check( nf90_put_att(wrncid, ginci_varid, MISS_VAL, -99.))
     ! * Grid geopotential height
    call check( nf90_def_var(wrncid,"Geopotential height",NF90_REAL,f3D_dimid,gz_varid))
    call check( nf90_put_att(wrncid, gz_varid, LONG_NAME, "Geopotential height"))
    call check( nf90_put_att(wrncid, gz_varid, SHORT_NAME, "Z"))
    call check( nf90_put_att(wrncid, gz_varid, UNITS, "km"))
    call check( nf90_put_att(wrncid, gz_varid, MISS_VAL, -99.))
    ! * Grid slope
    call check( nf90_def_var(wrncid,"SLOPE",NF90_REAL,f3D_dimid,gslope_varid))
    call check( nf90_put_att(wrncid, gslope_varid, LONG_NAME, "slope"))
    call check( nf90_put_att(wrncid, gslope_varid, SHORT_NAME, "SLOPE"))
    call check( nf90_put_att(wrncid, gslope_varid, UNITS, "Hellingsgraad"))
    call check( nf90_put_att(wrncid, gslope_varid, MISS_VAL, -99.))
     ! * Grid aspect
    call check( nf90_def_var(wrncid,"ASPECT",NF90_REAL,f3D_dimid,gaspect_varid))
    call check( nf90_put_att(wrncid, gaspect_varid, LONG_NAME, "aspect"))
    call check( nf90_put_att(wrncid, gaspect_varid, SHORT_NAME, "ASPECT"))
    call check( nf90_put_att(wrncid, gaspect_varid, UNITS, "NaN"))
    call check( nf90_put_att(wrncid, gaspect_varid, MISS_VAL, -99.))


    ! Putting global attributes
    call check( nf90_put_att(wrncid, NF90_GLOBAL, "title", "Satellite TB_HV"))
    call check( nf90_put_att(wrncid, NF90_GLOBAL, "orbit", trim(SAT%OrbitFileName)))
    call check( nf90_put_att(wrncid, NF90_GLOBAL, "sensor",trim(SAT%name)))
    call check( nf90_put_att(wrncid, NF90_GLOBAL, "creation",DATESTRING))
    call check( nf90_put_att(wrncid, NF90_GLOBAL, "contact","slvv@uni-bonn.de"))

    call check(nf90_enddef(wrncid))  ! End of definition

    !call check( nf90_put_var(wrncid, lon_dimid, SAT%lon))
    !call check( nf90_put_var(wrncid, lat_dimid, SAT%lat))
    !call check( nf90_put_var(wrncid, lev_dimid, CLMVARS%levsoi))
    !call check( nf90_put_var(wrncid, time_dimid, SAT%time))

    ! Putting satellite values in variables
    if(allocated(SAT%theta).and.allocated(SAT%lon_foprt).and.&
         allocated(SAT%lat_foprt).and.allocated(SAT%time).and.&
         allocated(SAT%TB_H).and.allocated(SAT%TB_V).and.&
         allocated(SAT%TBSAT_HV)) then

       call check(nf90_put_var(wrncid, lon_varid, SAT%lon_foprt))
       call check(nf90_put_var(wrncid, lat_varid, SAT%lat_foprt))
     !  call check(nf90_put_var(wrncid, glon_varid, SAT%lon))
     !  call check(nf90_put_var(wrncid, glat_varid, SAT%lat))
       call check(nf90_put_var(wrncid, inc_varid, SAT%theta))
       call check(nf90_put_var(wrncid, t_varid, SAT%time))

       call check(nf90_put_var(wrncid, tbh_varid, SAT%TBSAT_HV(:,1,:,:)))
       call check(nf90_put_var(wrncid, tbv_varid, SAT%TBSAT_HV(:,2,:,:)))
       call check(nf90_put_var(wrncid, gtbh_varid, SAT%TB_H(:,:,:,:)))
       call check(nf90_put_var(wrncid, gtbv_varid, SAT%TB_V(:,:,:,:)))
    else
       print*, 'Some SATELLITE variables are not allocated for NetCDF'
    end if

    ! Putting input values in variables
    if(allocated(CLMVARS%SWVL).and.allocated(CLMVARS%STL).and.&
         allocated(CLMVARS%SD).and.allocated(CLMVARS%RSN).and.&
         allocated(CLMVARS%TSKIN).and.allocated(CLMVARS%longxy).and.&
         allocated(CLMVARS%latixy)) then
       call check(nf90_put_var(wrncid, glon_varid, CLMVARS%longxy))
       call check(nf90_put_var(wrncid, glat_varid, CLMVARS%latixy))
       call check(nf90_put_var(wrncid, gsm_varid, CLMVARS%SWVL))
       call check(nf90_put_var(wrncid, gst_varid, CLMVARS%STL))
       call check(nf90_put_var(wrncid, gsd_varid, CLMVARS%SD))
       call check(nf90_put_var(wrncid, gRSN_varid, CLMVARS%RSN))
       call check(nf90_put_var(wrncid, gtskin_varid, CLMVARS%TSKIN))
       call check(nf90_put_var(wrncid, gsand_varid, CLMVARS%SAND))
       call check(nf90_put_var(wrncid, gclay_varid, CLMVARS%CLAY))
       call check(nf90_put_var(wrncid, gcvh_varid, CLMVARS%CVH))
       call check(nf90_put_var(wrncid, gcvl_varid, CLMVARS%CVL))
       call check(nf90_put_var(wrncid, gtvh_varid, CLMVARS%TVH))
       call check(nf90_put_var(wrncid, gtvl_varid, CLMVARS%TVL))
       call check(nf90_put_var(wrncid, glsm_varid, CLMVARS%LSM))
       call check(nf90_put_var(wrncid, gwater_varid, CLMVARS%WATER))
       call check(nf90_put_var(wrncid, glai_varid, CLMVARS%LAIL))
       call check(nf90_put_var(wrncid, gtair_varid, CLMVARS%TAIR))
       call check(nf90_put_var(wrncid, ginci_varid, CLMVARS%theta_inc))
       call check(nf90_put_var(wrncid, gz_varid, CLMVARS%Z))
       call check(nf90_put_var(wrncid, gslope_varid, CLMVARS%slope))
       call check(nf90_put_var(wrncid, gaspect_varid, CLMVARS%aspect))
    else
       print*, 'Some CLM variables are not allocated for NetCDF'
    end if
    ! c) Closing NC file
    call check(nf90_close(wrncid))

    return
  end subroutine write_satellite_operator
  ! --------------end of write satellite operator -----------------
  !:)

  ! ===============================================================
  ! SUBROUTINE to read the information for SATELLITE to be used
  ! The information is storaged as NetCDF file.
  ! ---------------------------------------------------------------
  ! (c) 2017 Pablo Saavedra G. (email: pablosaa@uni-bonn.de)
  ! UNIVERSITY OF BONN, GERMANY
  ! See LICENCE.TXT
  ! ---------------------------------------------------------------
  subroutine read_satellite_info(IN_fname,SAT)
    use clm4cmem, only : SATELLITE
    USE YOMCMEMPAR,    ONLY: LGPRINT  
    implicit none

    ! INPUT VARIABLES:
    character(len=*), intent(in) :: IN_fname
    type(SATELLITE), intent(out) :: SAT

    ! -- local variables
    integer :: rdncid, pix_dimid, inc_dimid
    integer :: varid, NPIXELS, NINCS

    call check( nf90_open(IN_fname, nf90_nowrite, rdncid))
    ! call check( nf90_open(current_observation_filename, nf90_nowrite, rdncid))
    ! Define the dimensions.
  
    call check( nf90_inq_dimid(rdncid, "NINC", inc_dimid))
    call check( nf90_inq_dimid(rdncid, "NPIXEL", pix_dimid))

    call check( nf90_inquire_dimension(rdncid,pix_dimid,len=NPIXELS))
    call check( nf90_inquire_dimension(rdncid,inc_dimid,len=NINCS))
 
    ! Allocating variables in SATELLITE:
    ! NOTE: SAT%TBSAT_HV is allocated by TBSAT_OPERATOR()
    if(allocated(SAT%theta)) deallocate(SAT%theta)
    allocate(SAT%theta(NINCS))
     
    if(allocated(SAT%lon_foprt)) deallocate(SAT%lon_foprt)
    allocate(SAT%lon_foprt(NPIXELS))
     
    if(allocated(SAT%lat_foprt)) deallocate(SAT%lat_foprt)
    allocate(SAT%lat_foprt(NPIXELS))

    if(allocated(SAT%incl_foprt)) deallocate(SAT%incl_foprt)
    allocate(SAT%incl_foprt(NPIXELS))

    ! ASSIGNING values to SAT members:
    SAT%OrbitFileName = trim(IN_fname)
    call check( nf90_inq_varid(rdncid,"THETA_INC", varid))
    call check( nf90_get_var(rdncid,varid,SAT%theta))
    IF (LGPRINT) WRITE(*,*) 'SAT%theta', SAT%theta
    call check( nf90_inq_varid(rdncid,"lon", varid))
    call check( nf90_get_var(rdncid,varid,SAT%lon_foprt))
    IF (LGPRINT) WRITE(*,*) 'SAT%lon_foprt', SAT%lon_foprt
    call check( nf90_inq_varid(rdncid,"lat", varid))
    call check( nf90_get_var(rdncid,varid,SAT%lat_foprt))
     
    call check( nf90_inq_varid(rdncid,"INCLI", varid))
    call check( nf90_get_var(rdncid,varid,SAT%incl_foprt))
     
    call check( nf90_get_att(rdncid,nf90_global,"SATELLITE_name",SAT%name))
    call check( nf90_get_att(rdncid,nf90_global,"Orbit_altitude_km",SAT%orbit))
    call check( nf90_get_att(rdncid,nf90_global,"Orbit_azimuth_deg",SAT%azimuth))
    call check( nf90_get_att(rdncid,nf90_global,"SENSOR_antenna_m",SAT%antenna))
    call check( nf90_get_att(rdncid,nf90_global,"SENSOR_wavelength_m",SAT%wavelength))
     
    call check( nf90_close(rdncid))
   
    return
  end subroutine read_satellite_info
  ! ------------- end of read satellite information ---------------
  !:)

  ! ===============================================================
  ! ---------------------------------------------------------------
  ! SUBROUTINE for creating the forcing for CMEM
  ! as NetCDF files in the specifiec output directory or ./
  ! (c) 2014 Pablo Saavedra G. (email: pablosaa@uni-bonn.de)
  ! UNIVERSITY OF BONN, GERMANY
  ! See LICENCE
  ! ---------------------------------------------------------------
  subroutine write_cmem_forcing(LS,nc_out)

    implicit none
    ! INPUT VARIABLES:
    type(CLM_DATA), intent(in) :: LS
    character(len=*), optional, intent(in) :: nc_out

    ! Local VARIABLES
    character (len=:), allocatable :: nc_path_out
    integer :: I, J, idxstr ! JJ, N
    integer :: lon_dimid, lat_dimid, time_dimid, lev_dimid
    integer :: wrncid, varid, z_varid
    integer :: NLEV, NLONS, NLATS, NTIME, NINC
    integer :: z_dimid(4)
    character (len = *), parameter :: UNITS = "units"
    character (len = *), parameter :: LONG_NAME = "long_name"
    character (len = *), parameter :: SHORT_NAME = "short_name"
    character (len = *), parameter :: MISS_VAL = "missing_value"
    character (len = *), parameter :: FILL_VAL = "_FillValue"
    character (len = *), parameter :: AXIS_NAME = "axis"
    character (len =100) :: FILETMP
    character (len = :), allocatable :: CLM_fname

    ! Populating variables:
    NLEV  = size(LS%levsoi)
    NLONS = size(LS%lons)
    NLATS = size(LS%lats)
    NTIME = size(LS%time)
    NINC  = size(LS%theta_inc, 3)
    CLM_fname = trim(LS%CLM_fname)
    IF(.not.present(nc_out)) THEN
       nc_path_out = './'
    ELSE
       nc_path_out = trim(nc_out)
    END IF

    if(SHOWINFO) print*, 'OUTPUTS at: -->'//trim(nc_path_out)
    ! -----------------------------------------------------------------
    ! 0-LSM: ASCII file with soil level profile
    !  -----------------------------------------------------------------
    open(unit=9,file=nc_path_out//'LSM_VERTICAL_RESOL.asc',status='REPLACE')
    write(9,*) (/ (LS%levsoi(I),I=1,NLEV) /)
    close(unit=9)

    ! 1-RD: Geopotential at surface [km]
    ! -----------------------------------------------------------------
    ! Create the file.
    if (SHOWINFO) print*, "Creating Z.nc"
    call check( nf90_create(nc_path_out//"Z.nc", nf90_clobber, wrncid))
    ! Define the dimensions.
    call check( nf90_def_dim(wrncid, "LONGITUDE", NLONS, lon_dimid))
    call check( nf90_def_dim(wrncid, "LATITUDE", NLATS, lat_dimid))
    call check( nf90_def_dim(wrncid, "TIME", 1, time_dimid))
    call check( nf90_def_dim(wrncid, "LEV", 1, lev_dimid))

    ! Define the variable
    call check( nf90_def_var(wrncid, "LONGITUDE", NF90_REAL,lon_dimid,varid))
    call check( nf90_put_att(wrncid, varid, LONG_NAME, "Longitud"))
    call check( nf90_put_att(wrncid, varid, SHORT_NAME, "x"))
    call check( nf90_put_att(wrncid, varid, UNITS, "degrees_east"))

    call check( nf90_def_var(wrncid, "LATITUDE", NF90_REAL,lat_dimid,varid))
    call check( nf90_put_att(wrncid, varid, LONG_NAME, "Latitude"))
    call check( nf90_put_att(wrncid, varid, SHORT_NAME, "y"))
    call check( nf90_put_att(wrncid, varid, UNITS, "degrees_north"))

    call check( nf90_def_var(wrncid, "LEV", NF90_REAL,lev_dimid,varid))
    call check( nf90_put_att(wrncid, varid, "axis", "1D"))
    call check( nf90_put_att(wrncid, varid, LONG_NAME, "Tile index"))
    call check( nf90_put_att(wrncid, varid, "point_spacing", "even"))
    call check( nf90_put_att(wrncid, varid, UNITS, "-"))

    call check( nf90_def_var(wrncid, "TIME", NF90_REAL,time_dimid,varid))
    call check( nf90_put_att(wrncid, varid, LONG_NAME, "time counter"))
    call check( nf90_put_att(wrncid, varid, SHORT_NAME, "time"))
    call check( nf90_put_att(wrncid, varid, UNITS, "doy"))

    ! * Geopotential
    z_dimid = (/ lon_dimid, lat_dimid, lev_dimid, time_dimid /)
    call check( nf90_def_var(wrncid, "Z", NF90_REAL,z_dimid,z_varid))
    call check( nf90_put_att(wrncid, z_varid, LONG_NAME, "Geopotential"))
    call check( nf90_put_att(wrncid, z_varid, SHORT_NAME, "Z"))
    call check( nf90_put_att(wrncid, z_varid, UNITS, "KM"))
    call check( nf90_put_att(wrncid, z_varid, MISS_VAL, 999.9))
    ! Putting global attributes
    call check( nf90_put_att(wrncid, NF90_GLOBAL, "title", "Geopotential"))
    call check( nf90_put_att(wrncid, NF90_GLOBAL, "model", LS%CLMver))
    call check( nf90_put_att(wrncid, NF90_GLOBAL, "origin",trim(LS%CLM_fname)))
    call check( nf90_put_att(wrncid, NF90_GLOBAL, "creation",LS%DATESTRING))

    call check(nf90_enddef(wrncid))  ! End of definition
    ! Putting the values in variables
    call check(nf90_put_var(wrncid, lon_dimid, LS%lons))
    call check(nf90_put_var(wrncid, lat_dimid, LS%lats))
    call check(nf90_put_var(wrncid, lev_dimid, 1))
    call check(nf90_put_var(wrncid, time_dimid, 1))
    call check(nf90_put_var(wrncid, z_varid, LS%Z))

    ! c) Closing NC file
    call check(nf90_close(wrncid))
    ! ------------------ End of Z.NC file --------------------------

    ! 2-RD: Snow Depth [m]
    ! -----------------------------------------------------------------
    ! Create the file.
    if (SHOWINFO) print*, "Creating SD.nc"
    call check( nf90_create(nc_path_out//"SD.nc", nf90_clobber, wrncid))
    ! Define the dimensions.
    call check( nf90_def_dim(wrncid, "LONGITUDE", NLONS, lon_dimid))
    call check( nf90_def_dim(wrncid, "LATITUDE", NLATS, lat_dimid))
    call check( nf90_def_dim(wrncid, "LEV", 1, lev_dimid))
    call check( nf90_def_dim(wrncid, "TIME", NTIME, time_dimid))
    ! Define the variable
    call check( nf90_def_var(wrncid, "LONGITUDE", NF90_REAL,lon_dimid,varid))
    call check( nf90_put_att(wrncid, varid, LONG_NAME, "Longitud"))
    call check( nf90_put_att(wrncid, varid, SHORT_NAME, "x"))
    call check( nf90_put_att(wrncid, varid, UNITS, "degrees_east"))

    call check( nf90_def_var(wrncid, "LATITUDE", NF90_REAL,lat_dimid,varid))
    call check( nf90_put_att(wrncid, varid, LONG_NAME, "Latitude"))
    call check( nf90_put_att(wrncid, varid, SHORT_NAME, "y"))
    call check( nf90_put_att(wrncid, varid, UNITS, "degrees_north"))

    call check( nf90_def_var(wrncid, "LEV", NF90_REAL,lev_dimid,varid))
    call check( nf90_put_att(wrncid, varid, "axis", "1D"))
    call check( nf90_put_att(wrncid, varid, LONG_NAME, "Tile index"))
    call check( nf90_put_att(wrncid, varid, "point_spacing", "even"))
    call check( nf90_put_att(wrncid, varid, UNITS, "-"))

    call check( nf90_def_var(wrncid, "TIME", NF90_REAL,time_dimid,varid))
    call check( nf90_put_att(wrncid, varid, AXIS_NAME, "T"))
    call check( nf90_put_att(wrncid, varid, LONG_NAME, "time counter"))
    call check( nf90_put_att(wrncid, varid, SHORT_NAME, "time"))
    call check( nf90_put_att(wrncid, varid, UNITS, "doy"))

    ! * Snow Depth
    z_dimid = (/ lon_dimid, lat_dimid,lev_dimid, time_dimid /)
    call check( nf90_def_var(wrncid, "SD", NF90_REAL,z_dimid,z_varid))
    call check( nf90_put_att(wrncid, z_varid, LONG_NAME, "Snow Depth"))
    call check( nf90_put_att(wrncid, z_varid, SHORT_NAME, "SD"))
    call check( nf90_put_att(wrncid, z_varid, UNITS, "m"))
    call check( nf90_put_att(wrncid, z_varid, MISS_VAL, -1E+34))
    call check( nf90_put_att(wrncid, z_varid, FILL_VAL, -1E+34))
    ! Putting global attributes
    call check( nf90_put_att(wrncid, NF90_GLOBAL, "title", "SNOWDP, snow height, time instanteneous"))
    call check( nf90_put_att(wrncid, NF90_GLOBAL, "model", LS%CLMver))
    call check( nf90_put_att(wrncid, NF90_GLOBAL, "origin",trim(LS%CLM_fname)))
    call check( nf90_put_att(wrncid, NF90_GLOBAL, "creation",LS%DATESTRING))

    call check(nf90_enddef(wrncid))

    ! Putting the values in variables
    call check(nf90_put_var(wrncid, lon_dimid, LS%lons))
    call check(nf90_put_var(wrncid, lat_dimid, LS%lats))
    call check(nf90_put_var(wrncid, time_dimid, LS%time))
    call check(nf90_put_var(wrncid, lev_dimid, 1))

    call check(nf90_put_var(wrncid, z_varid, LS%SD))

    ! c) Closing NC file
    call check(nf90_close(wrncid))
    ! ------------------ End of SD.NC file --------------------------

    ! 3-RD: Snow Density [kg/m^3]
    ! -----------------------------------------------------------------
    ! Create the file.
    if (SHOWINFO) print*, "Creating RSN.nc"
    call check( nf90_create(nc_path_out//"RSN.nc", nf90_clobber, wrncid))
    ! Define the dimensions.
    call check( nf90_def_dim(wrncid, "LONGITUDE", NLONS, lon_dimid))
    call check( nf90_def_dim(wrncid, "LATITUDE", NLATS, lat_dimid))
    call check( nf90_def_dim(wrncid, "LEV", 1, lev_dimid))
    call check( nf90_def_dim(wrncid, "TIME", NTIME, time_dimid))
    ! Define the variable
    call check( nf90_def_var(wrncid, "LONGITUDE", NF90_REAL,lon_dimid,varid))
    call check( nf90_put_att(wrncid, varid, LONG_NAME, "Longitud"))
    call check( nf90_put_att(wrncid, varid, SHORT_NAME, "x"))
    call check( nf90_put_att(wrncid, varid, UNITS, "degrees_east"))

    call check( nf90_def_var(wrncid, "LATITUDE", NF90_REAL,lat_dimid,varid))
    call check( nf90_put_att(wrncid, varid, LONG_NAME, "Latitude"))
    call check( nf90_put_att(wrncid, varid, SHORT_NAME, "y"))
    call check( nf90_put_att(wrncid, varid, UNITS, "degrees_north"))

    call check( nf90_def_var(wrncid, "LEV", NF90_REAL,lev_dimid,varid))
    call check( nf90_put_att(wrncid, varid, "axis", "1D"))
    call check( nf90_put_att(wrncid, varid, LONG_NAME, "Tile index"))
    call check( nf90_put_att(wrncid, varid, "point_spacing", "even"))
    call check( nf90_put_att(wrncid, varid, UNITS, "-"))

    call check( nf90_def_var(wrncid, "TIME", NF90_REAL,time_dimid,varid))
    call check( nf90_put_att(wrncid, varid, AXIS_NAME, "T"))
    call check( nf90_put_att(wrncid, varid, LONG_NAME, "time counter"))
    call check( nf90_put_att(wrncid, varid, SHORT_NAME, "time"))
    call check( nf90_put_att(wrncid, varid, UNITS, "doy"))

    ! * Snow Depth
    z_dimid = (/ lon_dimid, lat_dimid,lev_dimid, time_dimid /)
    call check( nf90_def_var(wrncid, "RSN", NF90_REAL,z_dimid,z_varid))
    call check( nf90_put_att(wrncid, z_varid, LONG_NAME, "Snow Density"))
    call check( nf90_put_att(wrncid, z_varid, SHORT_NAME, "SDen"))
    call check( nf90_put_att(wrncid, z_varid, UNITS, "kg/m^3"))
    call check( nf90_put_att(wrncid, z_varid, MISS_VAL, -1E+34))
    call check( nf90_put_att(wrncid, z_varid, FILL_VAL, -1E+34))
    ! Putting global attributes
    call check( nf90_put_att(wrncid, NF90_GLOBAL, "title", "SNOWLIQ, snow density, time instanteneous"))
    call check( nf90_put_att(wrncid, NF90_GLOBAL, "model", LS%CLMver))
    call check( nf90_put_att(wrncid, NF90_GLOBAL, "origin",trim(CLM_fname)))
    call check( nf90_put_att(wrncid, NF90_GLOBAL, "creation",LS%DATESTRING))

    call check(nf90_enddef(wrncid))

    ! Putting the values in variables
    call check(nf90_put_var(wrncid, lon_dimid, LS%lons))
    call check(nf90_put_var(wrncid, lat_dimid, LS%lats))
    call check(nf90_put_var(wrncid, time_dimid, LS%time))
    call check(nf90_put_var(wrncid, lev_dimid, 1))

    call check(nf90_put_var(wrncid, z_varid, LS%RSN))

    ! c) Closing NC file
    call check(nf90_close(wrncid))
    ! ------------------ End of RSN.NC file --------------------------

    ! 4-RD: Soil Temperature [K] (three layers for TESSEL)
    ! -----------------------------------------------------------------
    do I=1,nlev
       ! Create the file.
       if(I.lt.10) then
          write(FILETMP,'(A,I1,A)') "STL",I,".nc"
       else
          write(FILETMP,'(A,I2,A)') "STL",I,".nc"
       endif
       idxstr = index(FILETMP,'.',back=.true.)-1
       if (SHOWINFO) print*, "Creating ", FILETMP
       call check( nf90_create(nc_path_out//FILETMP, nf90_clobber, wrncid))
       ! Define the dimensions.
       call check( nf90_def_dim(wrncid, "LONGITUDE", NLONS, lon_dimid))
       call check( nf90_def_dim(wrncid, "LATITUDE", NLATS, lat_dimid))
       call check( nf90_def_dim(wrncid, "LEV", 1, lev_dimid))
       call check( nf90_def_dim(wrncid, "TIME", NTIME, time_dimid))
       ! Define the variable
       call check( nf90_def_var(wrncid, "LONGITUDE", NF90_REAL,lon_dimid,varid))
       call check( nf90_put_att(wrncid, varid, LONG_NAME, "Longitud"))
       call check( nf90_put_att(wrncid, varid, SHORT_NAME, "x"))
       call check( nf90_put_att(wrncid, varid, UNITS, "degrees_east"))

       call check( nf90_def_var(wrncid, "LATITUDE", NF90_REAL,lat_dimid,varid))
       call check( nf90_put_att(wrncid, varid, LONG_NAME, "Latitude"))
       call check( nf90_put_att(wrncid, varid, SHORT_NAME, "y"))
       call check( nf90_put_att(wrncid, varid, UNITS, "degrees_north"))

       call check( nf90_def_var(wrncid, "LEV", NF90_REAL,lev_dimid,varid))
       call check( nf90_put_att(wrncid, varid, "axis", "1D"))
       call check( nf90_put_att(wrncid, varid, LONG_NAME, "Tile index"))
       call check( nf90_put_att(wrncid, varid, "point_spacing", "even"))
       call check( nf90_put_att(wrncid, varid, UNITS, "-"))

       call check( nf90_def_var(wrncid, "TIME", NF90_REAL,time_dimid,varid))
       call check( nf90_put_att(wrncid, varid, AXIS_NAME, "T"))
       call check( nf90_put_att(wrncid, varid, LONG_NAME, "time counter"))
       call check( nf90_put_att(wrncid, varid, SHORT_NAME, "time"))
       call check( nf90_put_att(wrncid, varid, UNITS, "doy"))

       ! * Soil temperature
       z_dimid = (/ lon_dimid, lat_dimid,lev_dimid, time_dimid /)
       call check( nf90_def_var(wrncid, FILETMP(1:idxstr), NF90_REAL,z_dimid,z_varid))
       call check( nf90_put_att(wrncid, z_varid, LONG_NAME, "Soil Temp"))
       call check( nf90_put_att(wrncid, z_varid, SHORT_NAME, "ST"))
       call check( nf90_put_att(wrncid, z_varid, UNITS, "K"))
       call check( nf90_put_att(wrncid, z_varid, MISS_VAL, -1E+34))
       call check( nf90_put_att(wrncid, z_varid, FILL_VAL, -1E+34))
       ! Putting global attributes
       call check( nf90_put_att(wrncid, NF90_GLOBAL, "title", "TSOI, Soil Temperature, time instanteneous"))
       call check( nf90_put_att(wrncid, NF90_GLOBAL, "model", LS%CLMver))
       call check( nf90_put_att(wrncid, NF90_GLOBAL, "origin",trim(CLM_fname)))
       call check( nf90_put_att(wrncid, NF90_GLOBAL, "creation",LS%DATESTRING))

       call check(nf90_enddef(wrncid))

       ! Putting the values in variables
       call check(nf90_put_var(wrncid, lon_dimid, LS%lons))
       call check(nf90_put_var(wrncid, lat_dimid, LS%lats))
       call check(nf90_put_var(wrncid, time_dimid, LS%time))
       call check(nf90_put_var(wrncid, lev_dimid, 1))
       call check(nf90_put_var(wrncid, z_varid, LS%STL(:,:,I:I,:)))

       ! c) Closing NC file
       call check(nf90_close(wrncid))
    end do    ! end for I loop over files for soil temperature
    ! ------------------ End of STL/1,2,3/.NC file --------------------------

    ! 4.1-RD: Skin Temperature [K]
    ! -----------------------------------------------------------------
    ! Create the file.
    if (SHOWINFO) print*, "Creating TSKIN.nc"
    call check( nf90_create(nc_path_out//"TSKIN.nc", nf90_clobber, wrncid))
    ! Define the dimensions.
    call check( nf90_def_dim(wrncid, "LONGITUDE", NLONS, lon_dimid))
    call check( nf90_def_dim(wrncid, "LATITUDE", NLATS, lat_dimid))
    call check( nf90_def_dim(wrncid, "LEV", 1, lev_dimid))
    call check( nf90_def_dim(wrncid, "TIME", NTIME, time_dimid))

    ! Define the variable
    call check( nf90_def_var(wrncid, "LONGITUDE", NF90_REAL,lon_dimid,varid))
    call check( nf90_put_att(wrncid, varid, LONG_NAME, "Longitud"))
    call check( nf90_put_att(wrncid, varid, SHORT_NAME, "x"))
    call check( nf90_put_att(wrncid, varid, UNITS, "degrees_east"))

    call check( nf90_def_var(wrncid, "LATITUDE", NF90_REAL,lat_dimid,varid))
    call check( nf90_put_att(wrncid, varid, LONG_NAME, "Latitude"))
    call check( nf90_put_att(wrncid, varid, SHORT_NAME, "y"))
    call check( nf90_put_att(wrncid, varid, UNITS, "degrees_north"))

    call check( nf90_def_var(wrncid, "LEV", NF90_REAL,lev_dimid,varid))
    call check( nf90_put_att(wrncid, varid, "axis", "1D"))
    call check( nf90_put_att(wrncid, varid, LONG_NAME, "Tile index"))
    call check( nf90_put_att(wrncid, varid, "point_spacing", "even"))
    call check( nf90_put_att(wrncid, varid, UNITS, "-"))

    call check( nf90_def_var(wrncid, "TIME", NF90_REAL,time_dimid,varid))
    call check( nf90_put_att(wrncid, varid, LONG_NAME, "time counter"))
    call check( nf90_put_att(wrncid, varid, SHORT_NAME, "time"))
    call check( nf90_put_att(wrncid, varid, UNITS, "doy"))

    ! * Ground temperature
    z_dimid = (/ lon_dimid, lat_dimid, lev_dimid, time_dimid /)
    call check( nf90_def_var(wrncid, "TSKIN", NF90_REAL,z_dimid,z_varid))
    call check( nf90_put_att(wrncid, z_varid, LONG_NAME, "Skin Temperature"))
    call check( nf90_put_att(wrncid, z_varid, SHORT_NAME, "Ts"))
    call check( nf90_put_att(wrncid, z_varid, UNITS, "K"))
    call check( nf90_put_att(wrncid, z_varid, MISS_VAL, -1E+34))
    call check( nf90_put_att(wrncid, z_varid, FILL_VAL, -1E+34))
    ! Putting global attributes
    call check( nf90_put_att(wrncid, NF90_GLOBAL, "title", "Ground Temperature, time: instantaneous"))
    call check( nf90_put_att(wrncid, NF90_GLOBAL, "model", LS%CLMver))
    call check( nf90_put_att(wrncid, NF90_GLOBAL, "origin",trim(CLM_fname)))
    call check( nf90_put_att(wrncid, NF90_GLOBAL, "creation",LS%DATESTRING))

    call check(nf90_enddef(wrncid))  ! End of definition

    ! Putting the values in variables
    call check(nf90_put_var(wrncid, lon_dimid, LS%lons))
    call check(nf90_put_var(wrncid, lat_dimid, LS%lats))
    call check(nf90_put_var(wrncid, lev_dimid, 1))
    call check(nf90_put_var(wrncid, time_dimid, LS%time))

    call check(nf90_put_var(wrncid, z_varid, LS%TSKIN))

    ! c) Closing NC file
    call check(nf90_close(wrncid))
    ! ------------------ End of TSKIN.NC file --------------------------

    ! 5-RD: Soil water [m^3/m^3] (three layers...)
    ! -----------------------------------------------------------------
    do I=1,nlev
       ! Create the file.
       if(I.lt.10) then
          write(FILETMP,'(A,I1,A)') "SWVL",I,".nc"
       else
          write(FILETMP,'(A,I2,A)') "SWVL",I,".nc"
       endif
       idxstr = index(FILETMP,'.',back=.true.)-1
       if (SHOWINFO) print*, "Creating ", FILETMP
       call check( nf90_create(nc_path_out//FILETMP, nf90_clobber, wrncid))
       ! Define the dimensions.
       call check( nf90_def_dim(wrncid, "LONGITUDE", NLONS, lon_dimid))
       call check( nf90_def_dim(wrncid, "LATITUDE", NLATS, lat_dimid))
       call check( nf90_def_dim(wrncid, "LEV", 1, lev_dimid))
       call check( nf90_def_dim(wrncid, "TIME", NTIME, time_dimid))
       ! Define the variable
       call check( nf90_def_var(wrncid, "LONGITUDE", NF90_REAL,lon_dimid,varid))
       call check( nf90_put_att(wrncid, varid, LONG_NAME, "Longitud"))
       call check( nf90_put_att(wrncid, varid, SHORT_NAME, "x"))
       call check( nf90_put_att(wrncid, varid, UNITS, "degrees_east"))

       call check( nf90_def_var(wrncid, "LATITUDE", NF90_REAL,lat_dimid,varid))
       call check( nf90_put_att(wrncid, varid, LONG_NAME, "Latitude"))
       call check( nf90_put_att(wrncid, varid, SHORT_NAME, "y"))
       call check( nf90_put_att(wrncid, varid, UNITS, "degrees_north"))

       call check( nf90_def_var(wrncid, "LEV", NF90_REAL,lev_dimid,varid))
       call check( nf90_put_att(wrncid, varid, "axis", "1D"))
       call check( nf90_put_att(wrncid, varid, LONG_NAME, "Tile index"))
       call check( nf90_put_att(wrncid, varid, "point_spacing", "even"))
       call check( nf90_put_att(wrncid, varid, UNITS, "-"))

       call check( nf90_def_var(wrncid, "TIME", NF90_REAL,time_dimid,varid))
       call check( nf90_put_att(wrncid, varid, AXIS_NAME, "T"))
       call check( nf90_put_att(wrncid, varid, LONG_NAME, "time counter"))
       call check( nf90_put_att(wrncid, varid, SHORT_NAME, "time"))
       call check( nf90_put_att(wrncid, varid, UNITS, "doy"))

       ! * Volumetric Soil water
       z_dimid = (/ lon_dimid, lat_dimid,lev_dimid, time_dimid /)
       call check( nf90_def_var(wrncid, FILETMP(1:idxstr), NF90_REAL,z_dimid,z_varid))
       call check( nf90_put_att(wrncid, z_varid, LONG_NAME, "Soil Moisture"))
       call check( nf90_put_att(wrncid, z_varid, SHORT_NAME, "SWVL"))
       call check( nf90_put_att(wrncid, z_varid, UNITS, "m^3/m^3"))
       call check( nf90_put_att(wrncid, z_varid, MISS_VAL, -1E+34))
       call check( nf90_put_att(wrncid, z_varid, FILL_VAL, -1E+34))
       ! Putting global attributes
       call check( nf90_put_att(wrncid, NF90_GLOBAL, "title", "H2OSOI, Volumetric Soil Water, time instanteneous"))
       call check( nf90_put_att(wrncid, NF90_GLOBAL, "model", LS%CLMver))
       call check( nf90_put_att(wrncid, NF90_GLOBAL, "origin",trim(CLM_fname)))
       call check( nf90_put_att(wrncid, NF90_GLOBAL, "creation",LS%DATESTRING))

       call check(nf90_enddef(wrncid))

       ! Putting the values in variables
       call check(nf90_put_var(wrncid, lon_dimid, LS%lons))
       call check(nf90_put_var(wrncid, lat_dimid, LS%lats))
       call check(nf90_put_var(wrncid, time_dimid, LS%time))
       call check(nf90_put_var(wrncid, lev_dimid, 1))
       call check(nf90_put_var(wrncid, z_varid, LS%SWVL(:,:,I:I,:)))

       ! c) Closing NC file
       call check(nf90_close(wrncid))
    end do    ! end for I loop over files for soil temperature
    ! ------------------ End of SWVL/1,2,3,.../.NC file --------------------------

    ! ------------------------------------------------------------------------
    ! 6-RD: Texture
    ! Create the file.
    if (SHOWINFO) print*, "Creating SAND.nc"
    call check( nf90_create(nc_path_out//"sand.nc", nf90_clobber, wrncid))
    ! Define the dimensions.
    call check( nf90_def_dim(wrncid, "LONGITUDE", NLONS, lon_dimid))
    call check( nf90_def_dim(wrncid, "LATITUDE", NLATS, lat_dimid))

    ! Define the variable
    call check( nf90_def_var(wrncid, "LONGITUDE", NF90_REAL,lon_dimid,varid))
    call check( nf90_put_att(wrncid, varid, LONG_NAME, "Longitud"))
    call check( nf90_put_att(wrncid, varid, SHORT_NAME, "x"))
    call check( nf90_put_att(wrncid, varid, UNITS, "degrees_east"))

    call check( nf90_def_var(wrncid, "LATITUDE", NF90_REAL,lat_dimid,varid))
    call check( nf90_put_att(wrncid, varid, LONG_NAME, "Latitude"))
    call check( nf90_put_att(wrncid, varid, SHORT_NAME, "y"))
    call check( nf90_put_att(wrncid, varid, UNITS, "degrees_north"))

    ! * water fraction
    !z_dimid = (/ lon_dimid, lat_dimid /)
    call check( nf90_def_var(wrncid, "SAND", NF90_REAL,(/ lon_dimid, lat_dimid /),z_varid))
    call check( nf90_put_att(wrncid, z_varid, LONG_NAME, "Sand fraction"))
    call check( nf90_put_att(wrncid, z_varid, SHORT_NAME, "SAND"))
    call check( nf90_put_att(wrncid, z_varid, UNITS, "-"))
    call check( nf90_put_att(wrncid, z_varid, MISS_VAL, -1E+34))
    ! Putting global attributes
    call check( nf90_put_att(wrncid, NF90_GLOBAL, "title", "Sand % (surfdata)"))
    call check( nf90_put_att(wrncid, NF90_GLOBAL, "model", LS%CLMver))
    call check( nf90_put_att(wrncid, NF90_GLOBAL, "origin",trim(CLM_fname)))
    call check( nf90_put_att(wrncid, NF90_GLOBAL, "creation",LS%DATESTRING))

    call check(nf90_enddef(wrncid))  ! End of definition

    ! Putting the values in variables
    call check(nf90_put_var(wrncid, lon_dimid, LS%lons))
    call check(nf90_put_var(wrncid, lat_dimid, LS%lats))
    call check(nf90_put_var(wrncid, z_varid, LS%SAND))

    ! c) Closing NC file
    call check(nf90_close(wrncid))
    ! ------------------ End of SAND.NC file --------------------------

    ! Create the file.
    if (SHOWINFO) print*, "Creating CLAY.nc"
    call check( nf90_create(nc_path_out//"clay.nc", nf90_clobber, wrncid))
    ! Define the dimensions.
    call check( nf90_def_dim(wrncid, "LONGITUDE", NLONS, lon_dimid))
    call check( nf90_def_dim(wrncid, "LATITUDE", NLATS, lat_dimid))

    ! Define the variable
    call check( nf90_def_var(wrncid, "LONGITUDE", NF90_REAL,lon_dimid,varid))
    call check( nf90_put_att(wrncid, varid, LONG_NAME, "Longitud"))
    call check( nf90_put_att(wrncid, varid, SHORT_NAME, "x"))
    call check( nf90_put_att(wrncid, varid, UNITS, "degrees_east"))

    call check( nf90_def_var(wrncid, "LATITUDE", NF90_REAL,lat_dimid,varid))
    call check( nf90_put_att(wrncid, varid, LONG_NAME, "Latitude"))
    call check( nf90_put_att(wrncid, varid, SHORT_NAME, "y"))
    call check( nf90_put_att(wrncid, varid, UNITS, "degrees_north"))

    ! * water fraction
    !z_dimid = (/ lon_dimid, lat_dimid /)
    call check( nf90_def_var(wrncid, "CLAY", NF90_REAL,(/ lon_dimid, lat_dimid /),z_varid))
    call check( nf90_put_att(wrncid, z_varid, LONG_NAME, "Clay fraction"))
    call check( nf90_put_att(wrncid, z_varid, SHORT_NAME, "CLAY"))
    call check( nf90_put_att(wrncid, z_varid, UNITS, "-"))
    call check( nf90_put_att(wrncid, z_varid, MISS_VAL, -1E+34))
    ! Putting global attributes
    call check( nf90_put_att(wrncid, NF90_GLOBAL, "title", "clay % (surfdata)"))
    call check( nf90_put_att(wrncid, NF90_GLOBAL, "model", LS%CLMver))
    call check( nf90_put_att(wrncid, NF90_GLOBAL, "origin",trim(CLM_fname)))
    call check( nf90_put_att(wrncid, NF90_GLOBAL, "creation",LS%DATESTRING))

    call check(nf90_enddef(wrncid))  ! End of definition

    ! Putting the values in variables
    call check(nf90_put_var(wrncid, lon_dimid, LS%lons))
    call check(nf90_put_var(wrncid, lat_dimid, LS%lats))
    call check(nf90_put_var(wrncid, z_varid, LS%CLAY))

    ! c) Closing NC file
    call check(nf90_close(wrncid))
    ! ------------------ End of CLAY.NC file --------------------------

    ! -----------------------------------------------------------------
    ! 7-RD: Vegetation
    ! 7.1 : High/Low vegetation fraction
    ! Create the file.
    if (SHOWINFO) print*, "Creating ECOCVH.nc"
    call check( nf90_create(nc_path_out//"ECOCVH.nc", nf90_clobber, wrncid))
    ! Define the dimensions.
    call check( nf90_def_dim(wrncid, "LONGITUDE", NLONS, lon_dimid))
    call check( nf90_def_dim(wrncid, "LATITUDE", NLATS, lat_dimid))

    ! Define the variable
    call check( nf90_def_var(wrncid, "LONGITUDE", NF90_REAL,lon_dimid,varid))
    call check( nf90_put_att(wrncid, varid, LONG_NAME, "Longitud"))
    call check( nf90_put_att(wrncid, varid, SHORT_NAME, "x"))
    call check( nf90_put_att(wrncid, varid, UNITS, "degrees_east"))

    call check( nf90_def_var(wrncid, "LATITUDE", NF90_REAL,lat_dimid,varid))
    call check( nf90_put_att(wrncid, varid, LONG_NAME, "Latitude"))
    call check( nf90_put_att(wrncid, varid, SHORT_NAME, "y"))
    call check( nf90_put_att(wrncid, varid, UNITS, "degrees_north"))

    ! * fraction of high veg
    !z_dimid = (/ lon_dimid, lat_dimid /)
    call check( nf90_def_var(wrncid, "CVH", NF90_REAL,(/ lon_dimid, lat_dimid /),z_varid))
    call check( nf90_put_att(wrncid, z_varid, LONG_NAME, "High veg factor"))
    call check( nf90_put_att(wrncid, z_varid, SHORT_NAME, "cvh"))
    call check( nf90_put_att(wrncid, z_varid, UNITS, "-"))
    call check( nf90_put_att(wrncid, z_varid, MISS_VAL, -1E+34))
    ! Putting global attributes
    call check( nf90_put_att(wrncid, NF90_GLOBAL, "title", "=(2 & 1 ECOCLIM)"))
    call check( nf90_put_att(wrncid, NF90_GLOBAL, "model", LS%CLMver))
    call check( nf90_put_att(wrncid, NF90_GLOBAL, "origin",trim(CLM_fname)))
    call check( nf90_put_att(wrncid, NF90_GLOBAL, "creation",LS%DATESTRING))

    call check(nf90_enddef(wrncid))  ! End of definition

    ! Putting the values in variables
    call check(nf90_put_var(wrncid, lon_dimid, LS%lons))
    call check(nf90_put_var(wrncid, lat_dimid, LS%lats))
    call check(nf90_put_var(wrncid, z_varid, LS%CVH))

    ! c) Closing NC file
    call check(nf90_close(wrncid))
    ! ===========================================================
    ! Create the file.
    if (SHOWINFO) print*, "Creating ECOCVL.nc"
    call check( nf90_create(nc_path_out//"ECOCVL.nc", nf90_clobber, wrncid))
    ! Define the dimensions.
    call check( nf90_def_dim(wrncid, "LONGITUDE", NLONS, lon_dimid))
    call check( nf90_def_dim(wrncid, "LATITUDE", NLATS, lat_dimid))

    ! Define the variable
    call check( nf90_def_var(wrncid, "LONGITUDE", NF90_REAL,lon_dimid,varid))
    call check( nf90_put_att(wrncid, varid, LONG_NAME, "Longitud"))
    call check( nf90_put_att(wrncid, varid, SHORT_NAME, "x"))
    call check( nf90_put_att(wrncid, varid, UNITS, "degrees_east"))

    call check( nf90_def_var(wrncid, "LATITUDE", NF90_REAL,lat_dimid,varid))
    call check( nf90_put_att(wrncid, varid, LONG_NAME, "Latitude"))
    call check( nf90_put_att(wrncid, varid, SHORT_NAME, "y"))
    call check( nf90_put_att(wrncid, varid, UNITS, "degrees_north"))

    ! * fraction of high veg
    !z_dimid = (/ lon_dimid, lat_dimid /)
    call check( nf90_def_var(wrncid, "CVL", NF90_REAL,(/ lon_dimid, lat_dimid /),z_varid))
    call check( nf90_put_att(wrncid, z_varid, LONG_NAME, "Low veg factor"))
    call check( nf90_put_att(wrncid, z_varid, SHORT_NAME, "cvl"))
    call check( nf90_put_att(wrncid, z_varid, UNITS, "-"))
    call check( nf90_put_att(wrncid, z_varid, MISS_VAL, -1E+34))
    ! Putting global attributes
    call check( nf90_put_att(wrncid, NF90_GLOBAL, "title", "=(7 & 5)"))
    call check( nf90_put_att(wrncid, NF90_GLOBAL, "model", LS%CLMver))
    call check( nf90_put_att(wrncid, NF90_GLOBAL, "origin",trim(CLM_fname)))
    call check( nf90_put_att(wrncid, NF90_GLOBAL, "creation",LS%DATESTRING))

    call check(nf90_enddef(wrncid))  ! End of definition

    ! Putting the values in variables
    call check(nf90_put_var(wrncid, lon_dimid, LS%lons))
    call check(nf90_put_var(wrncid, lat_dimid, LS%lats))
    call check(nf90_put_var(wrncid, z_varid, LS%CVL))

    ! c) Closing NC file
    call check(nf90_close(wrncid))
    ! ------------------ End of ECOCVL.NC  and ECOCVH.NC file --------------------------

    ! 7.2 : High/Low vegetation type
    ! Create the file.
    if (SHOWINFO) print*, "Creating ECOTVL.nc"
    call check( nf90_create(nc_path_out//"ECOTVL.nc", nf90_clobber, wrncid))
    ! Define the dimensions.
    call check( nf90_def_dim(wrncid, "LONGITUDE", NLONS, lon_dimid))
    call check( nf90_def_dim(wrncid, "LATITUDE", NLATS, lat_dimid))

    ! Define the variable
    call check( nf90_def_var(wrncid, "LONGITUDE", NF90_REAL,lon_dimid,varid))
    call check( nf90_put_att(wrncid, varid, LONG_NAME, "Longitud"))
    call check( nf90_put_att(wrncid, varid, SHORT_NAME, "x"))
    call check( nf90_put_att(wrncid, varid, UNITS, "degrees_east"))

    call check( nf90_def_var(wrncid, "LATITUDE", NF90_REAL,lat_dimid,varid))
    call check( nf90_put_att(wrncid, varid, LONG_NAME, "Latitude"))
    call check( nf90_put_att(wrncid, varid, SHORT_NAME, "y"))
    call check( nf90_put_att(wrncid, varid, UNITS, "degrees_north"))

    ! * fraction of high veg
    !z_dimid = (/ lon_dimid, lat_dimid /)
    call check( nf90_def_var(wrncid, "TVL", NF90_REAL,(/ lon_dimid, lat_dimid /),z_varid))
    call check( nf90_put_att(wrncid, z_varid, LONG_NAME, "Low veg type"))
    call check( nf90_put_att(wrncid, z_varid, SHORT_NAME, "tvl"))
    call check( nf90_put_att(wrncid, z_varid, UNITS, "-"))
    call check( nf90_put_att(wrncid, z_varid, MISS_VAL, -1E+34))
    ! Putting global attributes
    call check( nf90_put_att(wrncid, NF90_GLOBAL, "title", "=4 C3 Grassland &7"))
    call check( nf90_put_att(wrncid, NF90_GLOBAL, "model", LS%CLMver))
    call check( nf90_put_att(wrncid, NF90_GLOBAL, "origin",trim(CLM_fname)))
    call check( nf90_put_att(wrncid, NF90_GLOBAL, "creation",LS%DATESTRING))

    call check(nf90_enddef(wrncid))  ! End of definition

    ! Putting the values in variables
    call check(nf90_put_var(wrncid, lon_dimid, LS%lons))
    call check(nf90_put_var(wrncid, lat_dimid, LS%lats))
    call check(nf90_put_var(wrncid, z_varid, LS%TVL))

    ! c) Closing NC file
    call check(nf90_close(wrncid))
    ! =========================================================================
    ! Create the file.
    if (SHOWINFO) print*, "Creating ECOTVH.nc"
    call check( nf90_create(nc_path_out//"ECOTVH.nc", nf90_clobber, wrncid))
    ! Define the dimensions.
    call check( nf90_def_dim(wrncid, "LONGITUDE", NLONS, lon_dimid))
    call check( nf90_def_dim(wrncid, "LATITUDE", NLATS, lat_dimid))

    ! Define the variable
    call check( nf90_def_var(wrncid, "LONGITUDE", NF90_REAL,lon_dimid,varid))
    call check( nf90_put_att(wrncid, varid, LONG_NAME, "Longitud"))
    call check( nf90_put_att(wrncid, varid, SHORT_NAME, "x"))
    call check( nf90_put_att(wrncid, varid, UNITS, "degrees_east"))

    call check( nf90_def_var(wrncid, "LATITUDE", NF90_REAL,lat_dimid,varid))
    call check( nf90_put_att(wrncid, varid, LONG_NAME, "Latitude"))
    call check( nf90_put_att(wrncid, varid, SHORT_NAME, "y"))
    call check( nf90_put_att(wrncid, varid, UNITS, "degrees_north"))

    ! * fraction of high veg
    !z_dimid = (/ lon_dimid, lat_dimid /)
    call check( nf90_def_var(wrncid, "TVH", NF90_REAL,(/ lon_dimid, lat_dimid /),z_varid))
    call check( nf90_put_att(wrncid, z_varid, LONG_NAME, "High veg type"))
    call check( nf90_put_att(wrncid, z_varid, SHORT_NAME, "tvh"))
    call check( nf90_put_att(wrncid, z_varid, UNITS, "-"))
    call check( nf90_put_att(wrncid, z_varid, MISS_VAL, -1E+34))
    ! Putting global attributes
    call check( nf90_put_att(wrncid, NF90_GLOBAL, "title", "=1&2 Decidious forests"))
    call check( nf90_put_att(wrncid, NF90_GLOBAL, "model", LS%CLMver))
    call check( nf90_put_att(wrncid, NF90_GLOBAL, "origin",trim(CLM_fname)))
    call check( nf90_put_att(wrncid, NF90_GLOBAL, "creation",LS%DATESTRING))

    call check(nf90_enddef(wrncid))  ! End of definition

    ! Putting the values in variables
    call check(nf90_put_var(wrncid, lon_dimid, LS%lons))
    call check(nf90_put_var(wrncid, lat_dimid, LS%lats))
    call check(nf90_put_var(wrncid, z_varid, LS%TVH))

    ! c) Closing NC file
    call check(nf90_close(wrncid))
    ! ------------------ End of ECOTVL.NC  and ECOTVH.NC file --------------------------

    ! 7.3 : Water fraction/Land fraction
    ! -----------------------------------------------------------------
    ! Create the file.
    if (SHOWINFO) print*, "Creating ECOWAT.nc"
    call check( nf90_create(nc_path_out//"ECOWAT.nc", nf90_clobber, wrncid))
    ! Define the dimensions.
    call check( nf90_def_dim(wrncid, "LONGITUDE", NLONS, lon_dimid))
    call check( nf90_def_dim(wrncid, "LATITUDE", NLATS, lat_dimid))

    ! Define the variable
    call check( nf90_def_var(wrncid, "LONGITUDE", NF90_REAL,lon_dimid,varid))
    call check( nf90_put_att(wrncid, varid, LONG_NAME, "Longitud"))
    call check( nf90_put_att(wrncid, varid, SHORT_NAME, "x"))
    call check( nf90_put_att(wrncid, varid, UNITS, "degrees_east"))

    call check( nf90_def_var(wrncid, "LATITUDE", NF90_REAL,lat_dimid,varid))
    call check( nf90_put_att(wrncid, varid, LONG_NAME, "Latitude"))
    call check( nf90_put_att(wrncid, varid, SHORT_NAME, "y"))
    call check( nf90_put_att(wrncid, varid, UNITS, "degrees_north"))

    ! * water fraction
    !z_dimid = (/ lon_dimid, lat_dimid /)
    call check( nf90_def_var(wrncid, "WATER", NF90_REAL,(/ lon_dimid, lat_dimid /),z_varid))
    call check( nf90_put_att(wrncid, z_varid, LONG_NAME, "Water fraction"))
    call check( nf90_put_att(wrncid, z_varid, SHORT_NAME, "WF"))
    call check( nf90_put_att(wrncid, z_varid, UNITS, "-"))
    call check( nf90_put_att(wrncid, z_varid, MISS_VAL, -1E+34))
    ! Putting global attributes
    call check( nf90_put_att(wrncid, NF90_GLOBAL, "title", "1-land_fraction"))
    call check( nf90_put_att(wrncid, NF90_GLOBAL, "model", LS%CLMver))
    call check( nf90_put_att(wrncid, NF90_GLOBAL, "origin",trim(CLM_fname)))
    call check( nf90_put_att(wrncid, NF90_GLOBAL, "creation",LS%DATESTRING))

    call check(nf90_enddef(wrncid))  ! End of definition

    ! Putting the values in variables
    call check(nf90_put_var(wrncid, lon_dimid, LS%lons))
    call check(nf90_put_var(wrncid, lat_dimid, LS%lats))
    call check(nf90_put_var(wrncid, z_varid, LS%WATER))

    ! c) Closing NC file
    call check(nf90_close(wrncid))

    ! ************* LSM file **************
    if (SHOWINFO) print*, "Creating LSM.nc"
    call check( nf90_create(nc_path_out//"LSM.nc", nf90_clobber, wrncid))
    ! Define the dimensions.
    call check( nf90_def_dim(wrncid, "LONGITUDE", NLONS, lon_dimid))
    call check( nf90_def_dim(wrncid, "LATITUDE", NLATS, lat_dimid))

    ! Define the variable
    call check( nf90_def_var(wrncid, "LONGITUDE", NF90_REAL,lon_dimid,varid))
    call check( nf90_put_att(wrncid, varid, LONG_NAME, "Longitud"))
    call check( nf90_put_att(wrncid, varid, SHORT_NAME, "x"))
    call check( nf90_put_att(wrncid, varid, UNITS, "degrees_east"))

    call check( nf90_def_var(wrncid, "LATITUDE", NF90_REAL,lat_dimid,varid))
    call check( nf90_put_att(wrncid, varid, LONG_NAME, "Latitude"))
    call check( nf90_put_att(wrncid, varid, SHORT_NAME, "y"))
    call check( nf90_put_att(wrncid, varid, UNITS, "degrees_north"))

    ! * land fraction
    !z_dimid = (/ lon_dimid, lat_dimid /)
    call check( nf90_def_var(wrncid, "LSM", NF90_REAL,(/ lon_dimid, lat_dimid /),z_varid))
    call check( nf90_put_att(wrncid, z_varid, LONG_NAME, "Land fraction"))
    call check( nf90_put_att(wrncid, z_varid, SHORT_NAME, "LF"))
    call check( nf90_put_att(wrncid, z_varid, UNITS, "-"))
    call check( nf90_put_att(wrncid, z_varid, MISS_VAL, -1E+34))
    ! Putting global attributes
    call check( nf90_put_att(wrncid, NF90_GLOBAL, "title", "LANDFRAC_PFT (surfdata)"))
    call check( nf90_put_att(wrncid, NF90_GLOBAL, "model", LS%CLMver))
    call check( nf90_put_att(wrncid, NF90_GLOBAL, "origin",trim(CLM_fname)))
    call check( nf90_put_att(wrncid, NF90_GLOBAL, "creation",LS%DATESTRING))

    call check(nf90_enddef(wrncid))  ! End of definition

    ! Putting the values in variables
    call check(nf90_put_var(wrncid, lon_dimid, LS%lons))
    call check(nf90_put_var(wrncid, lat_dimid, LS%lats))
    call check(nf90_put_var(wrncid, z_varid, LS%LSM))

    ! c) Closing NC file
    call check(nf90_close(wrncid))
    ! ------------------ End of ECOWAT.NC / LSM.nc file --------------------------
    ! 7.4-RD: low vegetation LAI
    ! Create the file.
    if (SHOWINFO) print*, "Creating ECOLAIL.nc"
    call check( nf90_create(nc_path_out//"ECOLAIL.nc", nf90_clobber, wrncid))
    ! Define the dimensions.
    call check( nf90_def_dim(wrncid, "LONGITUDE", NLONS, lon_dimid))
    call check( nf90_def_dim(wrncid, "LATITUDE", NLATS, lat_dimid))
    call check( nf90_def_dim(wrncid, "LEV", 1, lev_dimid))
    call check( nf90_def_dim(wrncid, "TIME", NTIME, time_dimid))

    ! Define the variable
    call check( nf90_def_var(wrncid, "LONGITUDE", NF90_REAL,lon_dimid,varid))
    call check( nf90_put_att(wrncid, varid, LONG_NAME, "Longitud"))
    call check( nf90_put_att(wrncid, varid, SHORT_NAME, "x"))
    call check( nf90_put_att(wrncid, varid, UNITS, "degrees_east"))

    call check( nf90_def_var(wrncid, "LATITUDE", NF90_REAL,lat_dimid,varid))
    call check( nf90_put_att(wrncid, varid, LONG_NAME, "Latitude"))
    call check( nf90_put_att(wrncid, varid, SHORT_NAME, "y"))
    call check( nf90_put_att(wrncid, varid, UNITS, "degrees_north"))

    call check( nf90_def_var(wrncid, "LEV", NF90_REAL,lev_dimid,varid))
    call check( nf90_put_att(wrncid, varid, "axis", "1D"))
    call check( nf90_put_att(wrncid, varid, LONG_NAME, "Tile index"))
    call check( nf90_put_att(wrncid, varid, "point_spacing", "even"))
    call check( nf90_put_att(wrncid, varid, UNITS, "-"))

    call check( nf90_def_var(wrncid, "TIME", NF90_REAL,time_dimid,varid))
    call check( nf90_put_att(wrncid, varid, LONG_NAME, "time counter"))
    call check( nf90_put_att(wrncid, varid, SHORT_NAME, "time"))
    call check( nf90_put_att(wrncid, varid, UNITS, "doy"))

    ! * LAI low vegetation
    z_dimid = (/ lon_dimid, lat_dimid, lev_dimid, time_dimid /)
    call check( nf90_def_var(wrncid, "LAI", NF90_REAL,z_dimid,z_varid))
    call check( nf90_put_att(wrncid, z_varid, LONG_NAME, "LAI low veg"))
    call check( nf90_put_att(wrncid, z_varid, SHORT_NAME, "LAI"))
    call check( nf90_put_att(wrncid, z_varid, UNITS, "-"))
    call check( nf90_put_att(wrncid, z_varid, MISS_VAL, -1E+34))
    call check( nf90_put_att(wrncid, z_varid, FILL_VAL, -1E+34))
    ! Putting global attributes
    call check( nf90_put_att(wrncid, NF90_GLOBAL, "title", "LAIL=TLAI*CVL (??)"))
    call check( nf90_put_att(wrncid, NF90_GLOBAL, "model", LS%CLMver))
    call check( nf90_put_att(wrncid, NF90_GLOBAL, "origin",trim(CLM_fname)))
    call check( nf90_put_att(wrncid, NF90_GLOBAL, "creation",LS%DATESTRING))

    call check(nf90_enddef(wrncid))  ! End of definition

    ! Putting the values in variables
    call check(nf90_put_var(wrncid, lon_dimid, LS%lons))
    call check(nf90_put_var(wrncid, lat_dimid, LS%lats))
    call check(nf90_put_var(wrncid, lev_dimid, 1))
    call check(nf90_put_var(wrncid, time_dimid, LS%time))

    call check(nf90_put_var(wrncid, z_varid, LS%LAIL))

    ! c) Closing NC file
    call check(nf90_close(wrncid))
    ! ------------------ End of ECOLAIL.NC file --------------------------

    ! -----------------------------------------------------------------
    ! 8.1-RD: 2m Air Temperature [K]
    ! -----------------------------------------------------------------
    ! Create the file.
    if (SHOWINFO) print*, "Creating 2T.nc"
    call check( nf90_create(nc_path_out//"2T.nc", nf90_clobber, wrncid))
    ! Define the dimensions.
    call check( nf90_def_dim(wrncid, "LONGITUDE", NLONS, lon_dimid))
    call check( nf90_def_dim(wrncid, "LATITUDE", NLATS, lat_dimid))
    call check( nf90_def_dim(wrncid, "LEV", 1, lev_dimid))
    call check( nf90_def_dim(wrncid, "TIME", NTIME, time_dimid))

    ! Define the variable
    call check( nf90_def_var(wrncid, "LONGITUDE", NF90_REAL,lon_dimid,varid))
    call check( nf90_put_att(wrncid, varid, LONG_NAME, "Longitud"))
    call check( nf90_put_att(wrncid, varid, SHORT_NAME, "x"))
    call check( nf90_put_att(wrncid, varid, UNITS, "degrees_east"))

    call check( nf90_def_var(wrncid, "LATITUDE", NF90_REAL,lat_dimid,varid))
    call check( nf90_put_att(wrncid, varid, LONG_NAME, "Latitude"))
    call check( nf90_put_att(wrncid, varid, SHORT_NAME, "y"))
    call check( nf90_put_att(wrncid, varid, UNITS, "degrees_north"))

    call check( nf90_def_var(wrncid, "LEV", NF90_REAL,lev_dimid,varid))
    call check( nf90_put_att(wrncid, varid, "axis", "1D"))
    call check( nf90_put_att(wrncid, varid, LONG_NAME, "Tile index"))
    call check( nf90_put_att(wrncid, varid, "point_spacing", "even"))
    call check( nf90_put_att(wrncid, varid, UNITS, "-"))

    call check( nf90_def_var(wrncid, "TIME", NF90_REAL,time_dimid,varid))
    call check( nf90_put_att(wrncid, varid, LONG_NAME, "time counter"))
    call check( nf90_put_att(wrncid, varid, SHORT_NAME, "time"))
    call check( nf90_put_att(wrncid, varid, UNITS, "doy"))

    ! * 2m air temperature
    z_dimid = (/ lon_dimid, lat_dimid, lev_dimid, time_dimid /)
    call check( nf90_def_var(wrncid, "TAIR", NF90_REAL,z_dimid,z_varid))
    call check( nf90_put_att(wrncid, z_varid, LONG_NAME, "Air Temperature"))
    call check( nf90_put_att(wrncid, z_varid, SHORT_NAME, "TAIR"))
    call check( nf90_put_att(wrncid, z_varid, UNITS, "K"))
    call check( nf90_put_att(wrncid, z_varid, MISS_VAL, -1E+34))
    call check( nf90_put_att(wrncid, z_varid, FILL_VAL, -1E+34))
    ! Putting global attributes
    call check( nf90_put_att(wrncid, NF90_GLOBAL, "title", "2m Air Temperature, time: instantaneous"))
    call check( nf90_put_att(wrncid, NF90_GLOBAL, "model", LS%CLMver))
    call check( nf90_put_att(wrncid, NF90_GLOBAL, "origin",trim(CLM_fname)))
    call check( nf90_put_att(wrncid, NF90_GLOBAL, "creation",LS%DATESTRING))

    call check(nf90_enddef(wrncid))  ! End of definition

    ! Putting the values in variables
    call check(nf90_put_var(wrncid, lon_dimid, LS%lons))
    call check(nf90_put_var(wrncid, lat_dimid, LS%lats))
    call check(nf90_put_var(wrncid, lev_dimid, 1))
    call check(nf90_put_var(wrncid, time_dimid, LS%time))

    call check(nf90_put_var(wrncid, z_varid, LS%TAIR))

    ! c) Closing NC file
    call check(nf90_close(wrncid))
    ! ------------------ End of 2T.NC file --------------------------

    ! -----------------------------------------------------------------
    ! 9.1-THETA_INC: SLOPE, ASPECT and THETA_inc
    ! -----------------------------------------------------------------
    ! Create the file.
    if (SHOWINFO) print*, "Creating THETA_INC.nc"
    call check( nf90_create(nc_path_out//"THETA_INC.nc", nf90_clobber, wrncid))
    ! Define the dimensions.
    call check( nf90_def_dim(wrncid, "LONGITUDE", NLONS, lon_dimid))
    call check( nf90_def_dim(wrncid, "LATITUDE", NLATS, lat_dimid))
    call check( nf90_def_dim(wrncid, "INCIDENCE", NINC, lev_dimid))

    ! Define the variable
    call check( nf90_def_var(wrncid, "LONGITUDE", NF90_REAL,lon_dimid,varid))
    call check( nf90_put_att(wrncid, varid, LONG_NAME, "Longitud"))
    call check( nf90_put_att(wrncid, varid, SHORT_NAME, "x"))
    call check( nf90_put_att(wrncid, varid, UNITS, "degrees_east"))

    call check( nf90_def_var(wrncid, "LATITUDE", NF90_REAL,lat_dimid,varid))
    call check( nf90_put_att(wrncid, varid, LONG_NAME, "Latitude"))
    call check( nf90_put_att(wrncid, varid, SHORT_NAME, "y"))
    call check( nf90_put_att(wrncid, varid, UNITS, "degrees_north"))

    call check( nf90_def_var(wrncid, "INCIDENCE", NF90_REAL,lev_dimid,varid))
    call check( nf90_put_att(wrncid, varid, LONG_NAME, "Incidence Angle"))
    call check( nf90_put_att(wrncid, varid, SHORT_NAME, "theta"))
    call check( nf90_put_att(wrncid, varid, UNITS, "degrees_zenith"))

    ! * Pixel-wise Incidence Angle
    call check( nf90_def_var(wrncid, "THETA_INC", NF90_REAL,(/ lon_dimid, lat_dimid, lev_dimid/),z_varid))
    call check( nf90_put_att(wrncid, z_varid, LONG_NAME, "Incidence Angle"))
    call check( nf90_put_att(wrncid, z_varid, SHORT_NAME, "TINC"))
    call check( nf90_put_att(wrncid, z_varid, UNITS, "DEG"))
    call check( nf90_put_att(wrncid, z_varid, MISS_VAL, -1E+34))
    call check( nf90_put_att(wrncid, z_varid, FILL_VAL, -1E+34))
    ! Putting global attributes
    call check( nf90_put_att(wrncid, NF90_GLOBAL, "title", "Pixel-wise Incidence Angle, 3rd-dim: global incidene angle"))
    call check( nf90_put_att(wrncid, NF90_GLOBAL, "model", LS%CLMver))
    call check( nf90_put_att(wrncid, NF90_GLOBAL, "origin",trim(CLM_fname)))
    call check( nf90_put_att(wrncid, NF90_GLOBAL, "creation",LS%DATESTRING))

    call check(nf90_enddef(wrncid))  ! End of definition

    ! Putting the values in variables
    call check(nf90_put_var(wrncid, lon_dimid, LS%lons))
    call check(nf90_put_var(wrncid, lat_dimid, LS%lats))
    call check(nf90_put_var(wrncid, lev_dimid, LS%theta))

    call check(nf90_put_var(wrncid, z_varid, LS%theta_inc))

    ! c) Closing NC file
    call check(nf90_close(wrncid))
    ! ------------------ End of THETA_INC.NC file -----------------
  end subroutine write_cmem_forcing
  ! ========= END SUBROUTINE write_cmem_forcing ===================
  !:)

  ! ===============================================================
  ! ANCILLARY SUBROUTINE FOR read_CLM_file
  ! ---------------------------------------------------------------
  ! Get dimensions of CLM NetCDF input file and check consistency
  ! (not very useful but needed to match CMEM structure)
  !
  ! (c) 2016 P. Saavedra Garfias (pablosaa@uni-bonn.de) UNI BONN
  ! see: LICENSE.TXT
  subroutine info_CLM_file(CLM_fname,Ntot,SURF,inhr,ilev)

    !use clm4cmem, only: NLONS, NLATS, NTIMES
    implicit none
    ! Info: INPUT arguments:
    character (len=*), intent(in) :: CLM_fname
    integer(kind=JPIM), intent(out) :: Ntot(3)
    character (len=*), intent(in), optional :: SURF
    integer, intent(in), optional :: inhr, ilev

    ! Info: Local variables
    integer :: NLEVLAK, NLEVSOI, NLEVPFT, NLATS, NLONS, NTIMES
    integer :: varid, lat_varid, lon_varid, doy_varid,levlak_varid, levsoi_varid, levpft_varid
    integer :: dimid, ncid, status
    character (len=10) :: dimname
    ! Info: following variables should be standard CLM dim names
    character (len = *), parameter :: LAT_NAME = "lat"
    character (len = *), parameter :: LON_NAME = "lon"
    character (len = *), parameter :: TIME_NAME = "time"
    character (len = *), parameter :: LEVLAK_NAME = "levlak"
    character (len = *), parameter :: LEVSOI_NAME = "levsoi"
    character (len = *), parameter :: PFT_NAME = "lsmpft"

    ! Info: *Open CLM netcdf input file:
    status = nf90_open(CLM_fname,nf90_nowrite, ncid)
    call check(status, trim(CLM_fname)//' CLM file.')
    ! Info: *Getting Latitude and Longitude dims and values
    call check( nf90_inq_dimid(ncid,LAT_NAME,dimid))
    call check( nf90_Inquire_Dimension(ncid,dimid,dimname,NLATS))
    call check( nf90_inq_dimid(ncid,LON_NAME,dimid))
    call check( nf90_Inquire_Dimension(ncid,dimid,dimname,NLONS))
    call check( nf90_inq_dimid(ncid,TIME_NAME,dimid))
    call check( nf90_Inquire_Dimension(ncid,dimid,dimname,NTIMES))

    ! Info: *Read dimension and value for soil levels
    call check( nf90_inq_dimid(ncid,LEVSOI_NAME,dimid))
    call check( nf90_Inquire_Dimension(ncid,dimid,dimname,NLEVSOI))

    ! Info: make checking for SURF too?
    ! ...

    call check( nf90_close(ncid))
    ! Info: checking consistency of dimensions with input arguments
    if(present(inhr)) then
       if(inhr.lt.1.or.inhr.gt.NTIMES) stop 'Non consistent TIME dim'
       NTIMES = 1
    else
       if(SHOWINFO) print*, 'Number of time dim ', NTIMES
    end if
    if(present(ilev)) then
       if(ilev.lt.1.or.ilev.gt.NLEVSOI) stop 'Non consistent SOIL LEVEL dim'
       NLEVSOI = ilev
    else
       if(SHOWINFO) print*, 'Number of levels ', NLEVSOI
    end if
    print*,'going out of info'
    !Ntot = NLONS*NLATS*NTIMES
    Ntot = (/NLONS,NLATS,NTIMES/)  ! output variable
    return
  end subroutine info_CLM_file
  ! ================ END INFO CLM NETCDF =============================
  ! ------------------------------------------------------------------
  ! SUBROUTINE TO READ CLM INPUT FILE AND PROCESS THE DATA FOR CMEM
  ! ------------------------------------------------------------------
  ! INPUTS:
  ! - clm_fname : NetCDF input CLM file
  ! - sm_factor : (OPTIONAL) [factor, bias] to convert soil moisture
  ! - nlev      : (OPTIONAL) number of SM level to consider
  ! - idxtime   : (OPTIONAL) index of the time dimension to consider
  ! OUTPUTS:
  ! *
  !
  ! ------------------------------------------------------------------

  subroutine read_CLM_file(CLM_fname,SAT,SURF, SMf,inhr,ilev, LS)

    !use constants, only : RHOw
    use clm4cmem, only: SATELLITE, local_incidence, gradientm
    implicit none

    ! INPUT arguments:
    character (len=*), intent(in) :: CLM_fname
    type(SATELLITE), intent(inout) :: SAT
    character (len=*), intent(in), optional :: SURF
    real, dimension(2), intent(in), optional :: SMf
    integer, intent(in), optional :: inhr, ilev
    type(CLM_DATA), intent(out) :: LS

    ! SUBROUTINE arguments:
    real, dimension(2) :: SM_faktor
    character(len=:), allocatable :: SURF_fname, nc_path_out
    CHARACTER(LEN=100) ::  CFINOUT='clm'        !! Input/output file format
    integer :: I, idxstr, idxtime, nlev
    integer :: ncid, suid
    integer :: dimid, ndims_in, nvars_in, ngatts_in, unlimdimid_in, status
    integer :: lon_dimid, lat_dimid, lev_dimid, time_dimid
    integer :: start3(3), count3(3), start4(4), count4(4)
    character (len=10) :: dimname, varname, varunits
    integer, dimension(8) :: datum_uhrzeit

    character (len=90) :: FILETMP
    character (len = *), parameter :: LAT_NAME = "lat"
    character (len = *), parameter :: LON_NAME = "lon"
    character (len = *), parameter :: DOY_NAME = "mcsec"
    character (len = *), parameter :: LEVLAK_NAME = "levlak"
    character (len = *), parameter :: LEVSOI_NAME = "levsoi"
    character (len = *), parameter :: PFT_NAME = "lsmpft"

    ! For the lat lon coordinate netCDF variables.
    integer :: NINC, NLATS, NLONS, NTIME, NLEVLAK, NLEVSOI, NLEVPFT

    real, dimension(:), allocatable :: time, levsoi
    !real, dimension(:), allocatable :: lats, lons, time, levlak, levsoi
    !real, dimension(:,:), allocatable :: Z, WATER, LSM, CVL, CVH, TVL, TVH, CLAY, SAND   ! 2D var
    !real, dimension(:,:), allocatable :: latixy, longxy, slope, aspect
    real, dimension(:,:,:), allocatable :: TMP
    real, dimension(:,:,:,:), allocatable :: TMP4D
    !real, dimension(:,:,:), allocatable :: TMP, LAI, theta_inc  ! 3D var
    !real, dimension(:,:,:,:), allocatable :: TMP4D, LAIL     ! 4D var
    !real, dimension(:,:,:,:), allocatable :: SD, RSN, STL  ! 4D var
    !real, allocatable :: TSKIN(:,:,:,:),SWVL(:,:,:,:), TAIR(:,:,:,:)  ! 4D var
    integer :: varid, lat_varid, lon_varid, doy_varid,levlak_varid, levsoi_varid, levpft_varid
    integer, parameter :: stime=1  ! which time-step (30=7.5hr)
    !  level indexes e.g. idxlev(nlev) = (/1, 2, 3/) 1st, 2nd, 3rd levels
    integer, dimension(:), allocatable :: idxlev
    !logical ::
    SHOWINFO = .true.

    ! * Open CLM netcdf input file:
    status = nf90_open(CLM_fname,nf90_nowrite, ncid)
    call check(status, trim(CLM_fname)//' CLM file.')

    IF(.not.present(SURF)) THEN
       ! Temporal for retrieving the surf_fname from global attribute in CLM file:
       idxstr = index(CLM_fname,'/',back=.true.)
       status = nf90_get_att(ncid, nf90_global, "Surface_dataset", FILETMP)
       call check(status, trim(FILETMP)//' problem?.')
       SURF_fname = CLM_fname(1:idxstr)//trim(FILETMP)
       if(SHOWINFO) print*, 'Attribute Surface_dataset = '//SURF_fname
    ELSE
       SURF_fname = trim(SURF)
    END IF

    ! * Open SURF netcdf input file:
    status = nf90_open(SURF_fname,nf90_nowrite, suid)
    call check(status, trim(SURF_fname)//' SURFACE file.')

    status = nf90_inquire(ncid,ndims_in, nvars_in, ngatts_in, unlimdimid_in)
    call check(status)
    if(SHOWINFO) print*, ndims_in, nvars_in, ngatts_in, unlimdimid_in
    ! Getting Latitude and Longitude dims and values
    status = nf90_inq_dimid(ncid,LAT_NAME,dimid)
    status = nf90_Inquire_Dimension(ncid,dimid,dimname,NLATS)

    status = nf90_inq_dimid(ncid,LON_NAME,dimid)
    status = nf90_Inquire_Dimension(ncid,dimid,dimname,NLONS)

    ! Get the varids of the latitude and longitude coordinate variables.
    call check( nf90_inq_varid(ncid, LAT_NAME, lat_varid) )
    call check( nf90_inq_varid(ncid, LON_NAME, lon_varid) )
    ! Read the latitude and longitude data.
    allocate(LS%lats(NLATS))
    allocate(LS%lons(NLONS))
    call check( nf90_get_var(ncid, lat_varid, LS%lats) )
    call check( nf90_get_var(ncid, lon_varid, LS%lons) )

    ! Read the time data (and convert to DoY?)
!!!    if(.not.present(inhr)) then !.AND.idxtime.eq.0) then
       call check( nf90_inq_dimid(ncid,'time',dimid))
       call check( nf90_Inquire_Dimension(ncid,dimid,dimname,NTIME))
    allocate(time(NTIME))
    call check( nf90_inq_varid(ncid,DOY_NAME,doy_varid))
    call check( nf90_get_var(ncid, doy_varid, time) )

       idxtime= 1  ! read all time-steps in CLM [1...NTIME]
!!!    else
       if(present(inhr)) then
          idxtime = inhr
          NTIME = 1
       end if

!!!    call check( nf90_get_var(ncid, doy_varid, time, start=(/idxtime/), count=(/NTIME/)) )

    ! Read dimension and value for lake levels
    call check( nf90_inq_dimid(ncid,LEVLAK_NAME,dimid))
    call check( nf90_Inquire_Dimension(ncid,dimid,dimname,NLEVLAK))
    call check( nf90_inq_varid(ncid,LEVLAK_NAME,levlak_varid))
    allocate(LS%levlak(NLEVLAK))
    call check( nf90_get_var(ncid,levlak_varid,LS%levlak))

    ! Read dimension and value for soil levels
    call check( nf90_inq_dimid(ncid,LEVSOI_NAME,dimid))
    call check( nf90_Inquire_Dimension(ncid,dimid,dimname,NLEVSOI))
    call check( nf90_inq_varid(ncid,LEVSOI_NAME,levsoi_varid))
    allocate(levsoi(NLEVSOI))
    call check( nf90_get_var(ncid,levsoi_varid,levsoi))

    ! Read dimension and value for PFT levels
    call check( nf90_inq_dimid(suid,PFT_NAME,dimid))
    call check( nf90_Inquire_Dimension(suid,dimid,dimname,NLEVPFT))

    ! Defining the range of data to read: case for the time with only 1 value
    start3 = (/1,1,idxtime/)
    count3 = (/NLONS,NLATS,NTIME/)
    start4 = (/1,1,1,idxtime/)
    count4 = (/NLONS,NLATS,NLEVSOI,NTIME/)

    ! Setting default values for input variables in case not passed in call
    if(.not.present(ilev)) then !nlev.eq.0) then
       nlev = NLEVSOI
    else
       nlev = ilev
    end if
    allocate(idxlev(nlev))
    idxlev = (/(I,I=1,nlev)/)

    if(.not.allocated(LS%levsoi)) allocate(LS%levsoi(nlev))
    LS%levsoi = levsoi(idxlev)

    if(.not.allocated(LS%time)) allocate(LS%time(NTIME))
    LS%time = time(idxtime)

    if(.not.present(SMf)) then
       SM_faktor=(/1,0/)
    else
       SM_faktor=SMf
    end if

    NINC = size(SAT%theta)
    if(allocated(SAT%time)) deallocate(SAT%time)
    allocate(SAT%time(NTIME))
    SAT%time = time(idxtime)

    ! ---------------------------------------------------
    ! Summary of parameters:
    if (SHOWINFO) then
       print*, 'CLM input: -->'//trim(CLM_fname)
       print*, 'SURF input:-->'//trim(SURF_fname)
       print*, 'Soil Moisture reduction factor: -->',SM_faktor(1),', Bias -->',SM_faktor(2)
       print*, 'Num. of soil levels: --> ', nlev
       print*, 'UTC Time span: --> ',time(idxtime)/3600.0,' to ', time(idxtime+NTIME-1)/3600.0, 'Hrs'
       print*, 'Number of incidences --> ', NINC
    end if
    ! ---------------------------------------------------

    ! READING VARIABLES:
    ! 1 -> Reading TOPO:
    allocate(LS%Z(NLONS,NLATS))
    call check( nf90_inq_varid(ncid,"topo",varid))
    call check( nf90_get_var(ncid,varid,LS%Z))
    LS%Z = LS%Z*1.0E-3   ! Geopotential convertion to KM

    ! 1.1 -> reading the longitude and latitude grid-boxes:
    allocate(LS%latixy(NLONS,NLATS), LS%longxy(NLONS,NLATS))
    call check( nf90_inq_varid(ncid,"latixy",varid))
    call check( nf90_get_var(ncid,varid,LS%latixy))
    call check( nf90_inq_varid(ncid,"longxy",varid))
    call check( nf90_get_var(ncid,varid,LS%longxy))

    ! 1.2 -> calculating slope and aspect from TOPO Z:
    allocate(LS%slope(NLONS,NLATS), LS%aspect(NLONS,NLATS))
    allocate(LS%theta_inc(NLONS,NLATS,NINC))

    call gradientm(LS%Z, LS%longxy, LS%latixy, LS%slope,&
         LS%aspect, LS%resol_km)
    !!write(*,*) (slope(I,1),I=1,10)
    !!print*, ':D'
    !!write(*,*) (aspect(I,1),I=1,10)
    ! 1.3 -> calculating pixel-wise incidence angle theta_inc:
    LS%theta = SAT%theta
    LS%theta_inc = local_incidence(LS%slope, LS%aspect, SAT)

    ! 2 -> Making up snow depth [m]
    if (allocated(TMP)) deallocate(TMP)
    allocate(TMP(NLONS,NLATS,NTIME))
    allocate(LS%SD(NLONS,NLATS,1,NTIME))
    call check( nf90_inq_varid(ncid,"SNOWDP",varid))  ! [m]
    call check( nf90_get_var(ncid,varid,TMP,start=start3, count=count3))
    where(TMP.EQ.1E+36) TMP = -1E+34
    LS%SD(:,:,1,:) = TMP

    ! 3 -> Making up snow density [kg/m^2]
    if (allocated(TMP)) deallocate(TMP)
    allocate(TMP(NLONS,NLATS,NTIME))
    allocate(LS%RSN(NLONS,NLATS,1,NTIME))
    call check( nf90_inq_varid(ncid,"SNOWLIQ",varid))   ! [kg/m^2]
    call check( nf90_get_var(ncid,varid,TMP,start=start3, count=count3))
    LS%RSN(:,:,1,:) = TMP/LS%SD(:,:,1,:)    ! Snow density [kg/m^3]
    where(TMP.EQ.1E+36) LS%RSN(:,:,1,:) = -1E+34
    where(TMP.LE.0.0.OR.LS%SD(:,:,1,:).EQ.-1E+34) LS%RSN(:,:,1,:) = -1E+34
    !LS%SD(:,:,1,:) = TMP/RHOw  ! SD/RHOw  !  Snow depth water equivalent [m]

    ! 4 -> Reading Soil Temperature [K]
    if (allocated(TMP4D)) deallocate(TMP4D)
    allocate(TMP4D(NLONS,NLATS,NLEVSOI,NTIME))
    allocate(LS%STL(NLONS,NLATS,nlev,NTIME))
    call check( nf90_inq_varid(ncid,"TSOI",varid))  ! Soil Temp [K]
    call check( nf90_get_var(ncid,varid,TMP4D,start=start4,count=count4))
    where(TMP4D.EQ.1E+36) TMP4D = -1E+34
    LS%STL=TMP4D     ! e.g. depths: 0.07, 0.21, 0.72m for TESSEL

    ! --> Making up the Skin Temperature [K]
    if (allocated(TMP)) deallocate(TMP)
    allocate(TMP(NLONS,NLATS,NTIME))
    call check( nf90_inq_varid(ncid,"TG",varid))  ! Ground Temp [K]
    call check( nf90_get_var(ncid,varid,TMP,start=start3,count=count3))
    where(TMP.EQ.1E+36) TMP = -1E+34
    allocate(LS%TSKIN(NLONS,NLATS,1,NTIME))
    LS%TSKIN(:,:,1,:) = TMP

    ! 5 -> Reading Volumetric Soil water [mm^3/mm^3]
    if (allocated(TMP4D)) deallocate(TMP4D)
    allocate(TMP4D(NLONS,NLATS,NLEVSOI,NTIME))
    allocate(LS%SWVL(NLONS,NLATS,nlev,NTIME))
    call check( nf90_inq_varid(ncid,"H2OSOI",varid))  ! Soil water [mm^3/mm^3]
    call check( nf90_get_var(ncid,varid,TMP4D,start=start4,count=count4))
    where(TMP4D.EQ.1E+36) TMP4D = -1E+34
    LS%SWVL=SM_faktor(1)*TMP4D !TEMPORAL (:,:,idxlev,:) ! depths: 0.07, 0.21, 0.72m TESSEL


    ! 6 ->  Texture
    if (allocated(LS%SAND)) deallocate(LS%SAND)
    allocate(LS%SAND(NLONS,NLATS))
    if (allocated(TMP)) deallocate(TMP)
    allocate(TMP(NLONS,NLATS,NLEVSOI))
    call check( nf90_inq_varid(suid,"PCT_SAND",varid)) ! new from surfdata
    call check( nf90_get_var(suid,varid,TMP))
    LS%SAND = TMP(:,:,1) ! new sand [%]
    if (allocated(LS%CLAY)) deallocate(LS%CLAY)
    allocate(LS%CLAY(NLONS,NLATS))
    if (allocated(TMP)) deallocate(TMP)
    allocate(TMP(NLONS,NLATS,NLEVSOI))
    call check( nf90_inq_varid(suid,"PCT_CLAY",varid)) ! new from surfdata
    call check( nf90_get_var(suid,varid,TMP))
    LS%CLAY = TMP(:,:,1) ! new clay [%]

    ! 7 -> Vegetation
    ! Vegetation types    | CLM VR0  | ECOCLIMAP | TESSEL :
    ! ---------------------------------------------------------------
    ! bare_soil           | 0        |           |
    ! needle_leaf_forest  | 1        |  2        | 3      (High Veg)
    ! broad_leaf_forest   | 7        |  1        | 4      (High Veg)
    ! grassland           | 14       |  5        | 2      (Low Veg)
    ! crops/agricul_land  | 15       |  6        | 1      (Low Veg)
    ! ---------------------------------------------------------------
    ! 7.1 -> high/low vegetation fraction
    if(allocated(LS%CVL)) deallocate(LS%CVL)
    allocate(LS%LAI(NLONS,NLATS,NTIME))
    if(allocated(TMP)) deallocate(TMP)
    allocate(TMP(NLONS,NLATS,NLEVPFT))
    call check( nf90_inq_varid(ncid,"TLAI",varid))  ! total projected LAI
    call check( nf90_get_var(ncid,varid,LS%LAI,start=start3,count=count3))
    call check( nf90_inq_varid(suid,"PCT_PFT",varid)) ! percent of PFT
    call check( nf90_get_var(suid,varid,TMP))
    TMP = TMP/100  ! change from percent to franciton

    allocate(LS%CVL(NLONS,NLATS))   ! TESSEL fraction of low veg
    allocate(LS%CVH(NLONS,NLATS))   ! TESSEL fraction of high veg
    LS%CVH = TMP(:,:,2)+TMP(:,:,8)  ! index CLM 1 and 7
    LS%CVL = TMP(:,:,15)+TMP(:,:,16)! index CLM 14 and 15

    ! 7.2 -> high/low vegetation type
    allocate(LS%TVL(NLONS,NLATS))   ! type of low veg
    allocate(LS%TVH(NLONS,NLATS))   ! type of high veg
    LS%TVL = TMP(:,:,15)*5+TMP(:,:,16)*6 ! ECOCLIM: 5=short_grass, 6=Crops
    LS%TVH = TMP(:,:,2)*2+TMP(:,:,8)*1   ! ECOCLIM: 2=Needleleaf, 1=Deciduous
    !TVL = TMP(:,:,15)*2+TMP(:,:,16)*1 ! TESSEL: 2=short_grass, 1=Crops
    !TVH = TMP(:,:,2)*3+TMP(:,:,8)*4   ! TESSEL: 3=Needleleaf, 4=Deciduous

    ! 7.3 -> water fraction [-] / land franction [-]
    if(allocated(LS%WATER)) deallocate(LS%WATER)
    allocate(LS%WATER(NLONS,NLATS))
    if(allocated(LS%LSM)) deallocate(LS%LSM)
    allocate(LS%LSM(NLONS,NLATS))
    call check( nf90_inq_varid(suid,"LANDFRAC_PFT",varid)) ! land fraction
    call check( nf90_get_var(suid,varid,LS%LSM))
    where(LS%LSM.EQ.1E+35)
       LS%WATER=-1E+34
    elsewhere
       LS%WATER = 1.0 - LS%LSM  ! converting land to water fraction
    endwhere
    ! 7.4 -> LAI of low vegetation
    if(allocated(LS%LAIL)) deallocate(LS%LAIL)
    allocate(LS%LAIL(NLONS,NLATS,1,NTIME))
    if(allocated(TMP)) deallocate(TMP)
    allocate(TMP(NLONS,NLATS,NTIME))
    call check( nf90_inq_varid(ncid,"TLAI",varid))  ! total LAI
    call check( nf90_get_var(ncid,varid,TMP, start=start3, count=count3)) ! total LAI
    LS%LAIL(:,:,1,:) = spread(LS%CVL,3,NTIME)*TMP ! low vegetation LAI

    ! 8 -> 2m Air Temperature [K]
    if (allocated(LS%TAIR)) deallocate(LS%TAIR)
    allocate(LS%TAIR(NLONS,NLATS,1,NTIME))
    call check( nf90_inq_varid(ncid,"TSA",varid))  ! 2m air temperature [K]
    call check( nf90_get_var(ncid,varid,LS%TAIR, start=start3, count=count3))
    where(LS%TAIR.EQ.1E+36) LS%TAIR = -1E+34

    !
    ! * Generating input data for CMEM
    call date_and_time(VALUES=datum_uhrzeit)
    write(LS%DATESTRING,'(I2.2"."I2.2"."I4 "T" I2.2":"I2.2":"I2.2)') datum_uhrzeit((/3,2,1,5,6,7/))


    status = nf90_close(ncid)
    IF(status /= nf90_noerr) call check(status)
    status = nf90_close(suid)
    IF(status /= nf90_noerr) call check(status)

    if (SHOWINFO) print*, 'NC Successfully close!'

    ! Deallocating variables
    if (allocated(time)) deallocate(time)
    if (allocated(levsoi)) deallocate(levsoi)
    if (allocated(TMP)) deallocate(TMP)
    if (allocated(TMP4D)) deallocate(TMP4D)

    ! --------------------------------------------
    ! Declaring global variables for CMEM which
    ! are needed for storing final variables as
    ! NetCDF ouput:
    NLONS_SM = NLONS
    NLATS_SM = NLATS
    NLVLS_SM = 1_JPIM
    NTIMES_SM = NTIME
    IF( ALLOCATED (xlons)) DEALLOCATE (xlons)
    ALLOCATE(xlons(NLONS))
    xlons = LS%lons
    IF( ALLOCATED (xlats)) DEALLOCATE (xlats)
    ALLOCATE(xlats(NLATS))
    xlats = LS%lats
    IF( ALLOCATED (xlvls)) DEALLOCATE (xlvls)
    ALLOCATE(xlvls(nlev))
    xlvls = LS%levsoi
    IF( ALLOCATED (xtimes)) DEALLOCATE (xtimes)
    ALLOCATE(xtimes(NTIME))
    xtimes = LS%time

    return

  end subroutine read_CLM_file
  ! ===================== END OF READING CLM FILE SUBROUTINE =========



  ! ---------------------------------------------------------------
  ! ---------------------------------------------------------------
  ! SUBROUTINE FOR CHECKING STATUS OF NETCDF INTRINSIC FUNCTIONS
  !----------------------------------------------------------------
  subroutine check(status, messages)
    use netcdf, only : nf90_strerror, nf90_noerr
    integer, intent ( in) :: status
    character (len=*), optional, intent(in) :: messages

    if(status /= nf90_noerr) then
       print *, 'Error: ',trim(nf90_strerror(status)),'->'//messages
       stop
    end if
    RETURN
  end subroutine check
  ! ----------------------------------------------------------------

end module rdclm_wrcmem
