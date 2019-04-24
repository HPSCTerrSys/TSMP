! *************************************************************************
! SET OF SUBROUTINES AND CLASSES FOR THE COUPLING OF LAND-SURFACE-MODEL (CLM)
! TO THE COMUNITY MICROWAVE EMISION MODEL (CMEM) AND THE SATELLITE SIMULATOR
! FOR PRODUCING SYNTHETIC DATA ALIKE NASA'S SMAP AND ESA'S SMOS MISSIONS.
! 
! Part of the DFG project FOR2131, sub-project P2.
! 
! (c) 2015-2017 Pablo Saavedra G. (email: pablosaa@uni-bonn.de)
! UNIVERSITY OF BONN, GERMANY
! *************************************************************************

module constants
  use parkind1, only: JPIM, JPRM
  implicit none

!  integer, parameter :: JPIM = SELECTED_INT_KIND(9)
!  integer, parameter :: JPRM = SELECTED_REAL_KIND(6,37)
  real, parameter :: Ro = 6371.0_JPRM  ! approx. Earth's radious [km]
  real, parameter :: pPI = 4._JPRM*atan(1._JPRM)
  real, parameter :: RAD2deg = 180._JPRM/pPI
  real, parameter :: DEG2rad = pPI/180._JPRM
  real, parameter :: RHOw = 1.0E3 ! Water density [kg/m^3]

!!$  &SMOS
!!$  name = 'SMOS',
!!$  orbit = 755.0_JPRM,
!!$  FOV  = 3.0_JPRM,    ! [deg]
!!$  antenna = 5.0_JPRM, ! [m]
!!$  azimuth = 312.18_JPRM-180._JPRM, ! [deg]
!!$  theta = 30 40 50  ! [deg]
!!$  /
end module constants

! =========================================================================
! MODULE INCLUDING ALL SUBROUTINES AND FUNCTIONS PART OF CLM4CMEM TOOLBOX
! LIST OF CLASSES:
! * SATELLITE: information regarding sensor and satellite
! * CLM_DATA : list of variables to extract from CLM and auxillary data
! LIST OF SUBROUTINES:
! * local_incidence() : calculate the local incidence depending on topography
! * gradientm() : to calculate the slope and aspect from given topography
! * memory_cmem_forcing() : put CLM_DATA class to memory to be used by CMEM
! *
! -------------------------------------------------------------------------
module clm4cmem
  use constants
  implicit none

  logical :: SHOWINFO = .true.
  logical :: DO_NETCDF = .true.

  type SATELLITE
     ! name: string with Sensor's name
     ! orbit: Orbit altitude of Satellite [km]
     ! FOV: field-of-view [deg]
     ! antenna: Sensor's antenna diameter [m]
     ! wavelength: Sensor's wavelength [m]
     ! theta: Vector with incidence angles [deg] (can be several)
     ! lon_foprt, lat_foprt: longitude and latitude of footprints [deg]
     ! azimuth: Along-track orbit azimuth [deg] from North CW
     ! incl_foprt: Footprint inclination angle [deg] from E CCW 
     ! TBSAT_HV: Bightness Temperature [K] (pixel,pol,inc,time)
     ! OrbitFileName: string with full-path of Satellite coordinates.
     sequence
     character(len=10) :: name
     real(kind=JPRM) :: orbit
     real(kind=JPRM) :: antenna, azimuth, wavelength
     real(kind=JPRM), allocatable, dimension(:) :: theta, incl_foprt, time
     real(kind=JPRM), allocatable, dimension(:) :: lon_foprt, lat_foprt
     real(kind=JPRM), allocatable, dimension(:,:,:,:) :: TBSAT_HV
     character(len=300) :: OrbitFileName
  end type SATELLITE

  type CLM_DATA
     sequence
     character (len=200) :: CLM_fname  ! CLM file name
     character (len=21) :: DATESTRING  ! date and time of creation
     character :: CLMver = "CLM3.5"

     ! 1D variables:
     real, dimension(2) :: resol_km
     real, dimension(:), allocatable :: theta
     real, dimension(:), allocatable :: lats, lons, time, levlak, levsoi
     ! 2D variables:
     real, dimension(:,:), allocatable :: Z, WATER, LSM, CVL, CVH
     real, dimension(:,:), allocatable :: TVL, TVH, CLAY, SAND
     real, dimension(:,:), allocatable :: latixy, longxy, slope, aspect
     ! 3D variables:
     real, dimension(:,:,:), allocatable :: LAI, theta_inc
     ! 4D variables:
     real, allocatable, dimension(:,:,:,:) :: TMP4D, LAIL
     real, allocatable, dimension(:,:,:,:) :: SD, RSN, STL 
     real, allocatable, dimension(:,:,:,:) :: TSKIN, SWVL, TAIR
  end type CLM_DATA
contains
  ! -------------------------------------------------------
  ! SUBROUTINE to compute the SLOPE and ASPECT for a given TOPO
  ! INPUTS:
  ! - Z(longitude, latitude) : TOPO data in [km]
  ! - longxy(longitude, latitude) : LONGITUDE grid in [deg] west to east
  ! - latixy(longitude, latitude) : LATITUDE  grid in [deg] south to north
  ! OUTPUTS:
  ! * slope(longitude, latitude) : SLOPE of terrain in [deg] from horizontal
  ! * aspect(longitude, latitude) : ASPECT for slope in [deg] from North
  ! * reslon_km and reslat_km     : grid resolution for lon and lat in [km]
  ! --
  ! (c) 2016, Pablo Saavedra (pablosaa@uni-bonn.de), Uni BONN
  ! ----------------------------------------------------------
  subroutine gradientm(Z, longxy, latixy, slope, aspect, LonLatresol_km)
    !use constants
    implicit none
    real, dimension(:,:), intent(in) :: Z, longxy, latixy
    real, dimension(:,:), intent(out) :: slope, aspect
    real, dimension(2), optional, intent(out) :: LonLatresol_km
    ! internal variables:
    integer(kind=JPIM) :: Nx, Ny
    !real, parameter :: Ro = 6371.0  ! approx. Earth's radious [km]
    !real, parameter :: PI = 4.*atan(1.)
    !real, parameter :: RAD2deg = 180./PI
    !real, parameter :: DEG2rad = PI/180.
    real(kind=JPRM), dimension(size(Z,1),size(Z,2)) :: delWE, delSN, faktor
    real(kind=JPRM), dimension(size(Z,1),size(Z,2)) :: delX, delY, Zwe, Zsn

    Nx = size(Z,1)
    Ny = size(Z,2)
    ! Calculating grid West-East angular resolution:
    delWE = 0.0_JPRM
    delWE(1,:) = longxy(2,:)-longxy(1,:)
    delWE(Nx,:) = longxy(Nx,:)-longxy(Nx-1,:)
    delWE(2:Nx-1,:) = (longxy(3:Nx,:)-longxy(1:Nx-2,:))/2.
    delWE = delWE*DEG2rad   ! [rad]
    ! Calculating grid South-North angular resolution:
    delSN = 0.0_JPRM
    delSN(:,1) = latixy(:,2)-latixy(:,1)
    delSN(:,Ny) = latixy(:,Ny)-latixy(:,Ny-1)
    delSN(:,2:Ny-1) = (latixy(:,3:Ny)-latixy(:,1:Ny-2))/2.
    delSN = delSN*DEG2rad   ! [rad]
    ! latitudinal adjustment factor for distance over meridians:
    faktor = cos(latixy*DEG2rad)
    delX = Ro*faktor*delWE   ! grid-res over longitude [km]
    delY = Ro*delSN          ! grid-res over latitude [km]
    ! Calculating grid West-East grid gradient:
    Zwe = 0.0_JPRM
    Zwe(1,:) = Z(2,:)-Z(1,:)
    Zwe(Nx,:) = Z(Nx,:)-Z(Nx-1,:)
    Zwe(2:Nx,:) = (Z(3:Nx,:)-Z(1:Nx-2,:))/2.  ! [km]
    ! Calculating grid South-North grid gradient:
    Zsn = 0.0_JPRM
    Zsn(:,1) = Z(:,2)-Z(:,1)
    Zsn(:,Ny) = Z(:,Ny)-Z(:,Ny-1)
    Zsn(:,2:Ny-1) = (Z(:,3:Ny)-Z(:,1:Ny-2))/2.  ! [km]

    ! Calculating the slope in [deg]
    slope = sqrt((Zwe/delX)**2 + (Zsn/delY)**2)*RAD2deg
    ! Calculating the aspect in [deg] from North Clockwise
    aspect = 270.0_JPRM - RAD2deg*atan2(Zsn/delY,Zwe/delX)
    where(aspect.ge.360.) aspect = aspect - 360.0_JPRM
    
    if(present(LonLatresol_km)) LonLatresol_km(1) = sum(delX(:,1))/Nx
    if(present(LonLatresol_km)) LonLatresol_km(2) = sum(delY(1,:))/Ny

    return
  end subroutine gradientm

  ! -----------------------------------------------------------------
  ! FUNCTION to calculate the local incidence angle taking into
  ! account the slope and aspect derived from the TOPO grid
  ! INPUTS:
  ! - slope : matrix (longitude,latitude) in [deg] from horizontal
  ! - aspect: matrix (longitude,latitude) in [deg] from North CW
  ! - SAT: SATELLITE structure including satellite information
  ! OUTPUTS:
  ! * theta_hat: matrix (longitude,latitude,Ninc) in [deg] from Normal
  ! Where Ninc is the number of incidence angles to consider
  !
  ! (c) 2017 Pablo Saavedra G. (pablosaa@uni-bonn.de)
  !-------------------------------------------------------------------
  function local_incidence(slope,aspect,SAT) result(theta_hat)
    implicit none

    real, dimension(:,:), intent(in) :: slope, aspect
    type(SATELLITE), intent(in) :: SAT
    real(kind=JPRM), dimension(size(slope,1),size(slope,2),size(SAT%theta,1)) :: theta_hat
    real(kind=JPRM), dimension(size(SAT%theta,1)) :: theta_i
    real(kind=JPRM), dimension(size(slope,1),size(slope,2)) :: Axy
    integer(kind=JPIM) :: k, Ninc

    Ninc = size(SAT%theta)

    theta_i = asin((1.0_JPRM+(SAT%orbit/Ro))*sin(DEG2rad*SAT%theta))  ! [rad]
    print*,'Number of angles: ',size(SAT%theta), shape(theta_hat)
    Axy = sin(aspect*DEG2rad)*sin(SAT%azimuth*DEG2rad) + &
         cos(aspect*DEG2rad)*cos(SAT%azimuth*DEG2rad)
    do k=1,Ninc
       theta_hat(:,:,k) = acos(sin(slope*DEG2rad)*sin(theta_i(k))*Axy + &
            cos(slope*DEG2rad)*cos(theta_i(k)))*RAD2deg
    end do

    !!print*, theta_hat(7,1,1:Ninc)
    return
  end function local_incidence
  ! ============= END FUNCTION FOR THETA_INC ======================

  ! ---------------------------------------------------------------
  ! SUBROUTINE for creating the forcing for CMEM
  ! Either in memory or as NetCDF files
  ! ---------------------------------------------------------------
  subroutine memory_cmem_forcing(LS)
    use yomcmempar, only: CIDVEG, CITVEG, CIATM, LGPRINT,CNAMEID &
         & , LOMASK_OCEAN,LOMASK_AUTO, nlay_soil_ls
    use yomcmematm, only: fZ
    use yomcmemfields, only: N, JJ, CLNAME, fTVL, fTVH,  fs_laiL, ftfrac, fs_laiL &
         &  , fsnowd, frsnow, fwater, ftl_lsm, ftair,ftveg &
         & , ftskin,fwc_lsm,fsand,fclay,fwater, mask, ftheta_inc ! PSG: insert theta_inc
    use yomcmemsoil, only: sal_sea, z_lsm
    implicit none

    ! INPUT VARIABLES:
    type(CLM_DATA), intent(in) :: LS

    ! LOCAL VARIABLES:
    integer :: I, J !, idxstr, JJ, N
    integer :: NLEV , NLONS, NLATS, NTIMES, NINC
    real(kind=JPRM), allocatable, dimension(:,:) :: fncfield
    real(kind=JPRM), allocatable, dimension(:,:,:,:) :: var_in_tr

    ! TEMPORAL VARIABLES FROM CMEM for VEGTABLE subroutine
    real(KIND=JPRM) :: rvcov(0:20)
    real(KIND=JPRM) :: rvlai(0:20)
    real(KIND=JPRM) :: b(0:7)
    real(KIND=JPRM) :: b1(0:7)
    real(KIND=JPRM) :: b2(0:7)
    real(KIND=JPRM) :: b3(0:7)
    real(KIND=JPRM) :: VWC(0:7)
    real(KIND=JPRM) :: ZNrh(0:7)
    real(KIND=JPRM) :: ZNrv(0:7)
    real(KIND=JPRM) :: Ztth(0:7)
    real(KIND=JPRM) :: Zttv(0:7)
    real(KIND=JPRM) :: Zhr(0:7)
    real(KIND=JPRM) :: Zw_eff(0:7,2)
    integer(KIND=JPIM) :: RVTV(0:20)
    ! -- 
    REAL(KIND=JPRM),ALLOCATABLE :: cvegl(:), cvegh(:)
    REAL(KIND=JPRM),ALLOCATABLE :: fbare(:)
    REAL(KIND=JPRM),ALLOCATABLE :: fvegl(:)
    REAL(KIND=JPRM),ALLOCATABLE :: fvegh(:)
    REAL(KIND=JPRM),ALLOCATABLE :: sncov(:)
    ! END OF VARIABLE DEFINITION

    NLEV  = size(LS%levsoi)
    NLONS = size(LS%lons)
    NLATS = size(LS%lats)
    NTIMES = size(LS%time)
    NINC  = size(LS%theta_inc, 3)

    ! temporal OPT variables assigned (should be done by cmem_init)
    !N=NLONS*NLATS*NTIME
    !!! if(.not.allocated(fTVL)) allocate(fTVL(N)) ! must be in cmem_main
    !!! if(.not.allocated(fTVH)) allocate(fTVH(N)) ! must be in cmem_main
    !!! if(.not.allocated(fwater)) allocate(fwater(N)) ! must be in cmem_main
    !!! if(.not.allocated(ftfrac)) allocate(ftfrac(N,7))  ! must be in cmem_main
    !!! CIATM='Prellarin'  ! should be in cmem_init or input
    !!! CIDVEG='Ecoclimap'  ! should be in cmem_init or input
    !!! CITVEG='Tsurf'  ! should be in cmem_init or input
    ! end of OPT variables assigned
    print*, N, NLONS, NLATS, NLEV, NTIMES, NINC !, CLM_fname
    if(SHOWINFO) print*,LS%DATESTRING

    ! ===============================================================
    ! Setting CMEM variables to memory (no NetCDF files are created)
    ! ---------------------------------------------------------------
    MASK = 1_JPIM   ! Initializaton for mask
    ! 0- z_lsm var: depth of the land surface model layers (m)
    z_lsm = LS%levsoi(1:nlay_soil_ls)
    WRITE(*,*) ' CLM IO: LSM layers depth: ',z_lsm(:) ! CLMcase
    ! 1- fZ variable: -- topography Z [km] ------
    JJ = NLONS*NLATS
    if(.not.allocated(fZ)) allocate(fZ(N))
    if(allocated(fncfield)) deallocate(fncfield)
    allocate(fncfield(JJ,1))
    fncfield = reshape(LS%Z,(/JJ,1/))
    print*, 'fcnfield : ', shape(fncfield)
    do I=1,N-JJ+1,JJ ! filling up fZ with NTIME times
       fZ(I:) = fncfield(:,1)
    end do
    if(SHOWINFO) print*, 'out of fZ!'
    ! ===== end variable fZ ====================

    ! 2- fsnowd variable: -- snow depth SD [m] --------
    JJ = NLONS*NLATS
    if(.not.allocated(fsnowd)) allocate(fsnowd(N))
    if(allocated(fncfield)) deallocate(fncfield)
    allocate(fncfield(JJ,1))
    fncfield = reshape(LS%SD,(/JJ,1/))
    do I=1,N-JJ+1,JJ ! filling up fZ with NTIME times
       fsnowd(I:) = max(1.E-2_JPRM,fncfield(:,1))
    end do
    where(fsnowd.eq.1.E-2_JPRM) fsnowd = 0.0_JPRM
    if(SHOWINFO) print*, 'out of fsnowd!'
    ! ===== end variable fsnowd ====================

    ! 3- frsnow variable: -- snow density [kg/m^3] ---
    JJ = NLONS*NLATS
    if(.not.allocated(frsnow)) allocate(frsnow(N))
    if(allocated(fncfield)) deallocate(fncfield)
    allocate(fncfield(JJ,1))
    fncfield = reshape(LS%RSN,(/JJ,1/))
    do I=1,N-JJ+1,JJ ! filling up fZ with NTIME times
       frsnow(I:) = fncfield(:,1)
    end do
    if(SHOWINFO) print*, 'out of frsnow!'
    ! ===== end variable frsnow ====================

    ! 4- ftl_lsm variable: -- Soil Temperature [K] ---
    if(.not.allocated(ftl_lsm)) allocate(ftl_lsm(N,NLEV))
    if(allocated(var_in_tr)) deallocate(var_in_tr)
    allocate(var_in_tr(NLONS,NLATS,NTIMES,NLEV))
    var_in_tr = reshape(LS%STL,(/NLONS,NLATS,NTIMES,NLEV/),ORDER=(/1,2,4,3/))
    ftl_lsm = reshape(var_in_tr,(/N,NLEV/))
    if(SHOWINFO) print*, 'out of ftl_lsm!'
    ! ===== end variable ftl_lsm ======================

    ! 4.1- ftskin variable: -- Skin Temperature [k] ---
    JJ = NLONS*NLATS
    if(.not.allocated(ftskin)) allocate(ftskin(N))
    if(allocated(fncfield)) deallocate(fncfield)
    allocate(fncfield(JJ,1))
    fncfield = reshape(LS%TSKIN,(/JJ,1/))
    do I=1,N-JJ+1,JJ ! filling up fZ with NTIME times
       ftskin(I:) = fncfield(:,1)
    end do
    ! for LOMAST_AUTO to eliminate unrealistic points:
    do I=1,NLEV
       where(ftl_lsm(:,I).lt.100.0_JPRM.or.ftskin.lt.100.0_JPRM) mask(:) = 0_JPIM
    end do
    if(SHOWINFO) print*, 'out of ftskin!'
    ! ===== end variable ftskin =======================

    ! 5- fwc_lsm variable: -- Soil moisture [m^3/m^3]
    if(.not.allocated(fwc_lsm)) allocate(fwc_lsm(N,NLEV))
    if(allocated(var_in_tr)) deallocate(var_in_tr)
    allocate(var_in_tr(NLONS,NLATS,NTIMES,NLEV))
    var_in_tr = reshape(LS%SWVL,(/NLONS,NLATS,NTIMES,NLEV/),ORDER=(/1,2,4,3/))
    fwc_lsm = reshape(var_in_tr,(/N,NLEV/))
    ! update the mask to eliminate unrealistic points
    do I=1,NLEV
       where(fwc_lsm(:,I).lt.0.0_JPRM) mask(:) = 0_JPIM
    end do
    fwc_lsm = max(0.0_JPRM,fwc_lsm)  ! to ensure only positive values
    if(SHOWINFO) print*, 'out of fwc_lsm!'
    ! ===== end variable fwc_lsm =======================

    ! 6.1- fsand variable: -- Soil Texture, SAND [%] --
    JJ = NLONS*NLATS
    if(.not.allocated(fsand)) allocate(fsand(N))
    if(allocated(fncfield)) deallocate(fncfield)
    allocate(fncfield(JJ,1))
    fncfield = reshape(LS%SAND,(/JJ,1/))
    do I=1,N-JJ+1,JJ ! filling up fZ with NTIME times
       fsand(I:) = fncfield(:,1)
    end do
    if(SHOWINFO) print*, 'out of fsand!'
    ! ===== end variable fsand =======================
    ! 6.2- fclay variable: -- Soil Texture, CLAY [%] --
    JJ = NLONS*NLATS
    if(.not.allocated(fclay)) allocate(fclay(N))
    if(allocated(fncfield)) deallocate(fncfield)
    allocate(fncfield(JJ,1))
    fncfield = reshape(LS%CLAY,(/JJ,1/))
    do I=1,N-JJ+1,JJ ! filling up fZ with NTIME times
       fclay(I:) = fncfield(:,1)
    end do
    if(SHOWINFO) print*, 'out of fclay!'
    ! ===== end variable fclay =======================

    ! 7- variables: -- Vegetation --
    if(.not.allocated(fbare)) allocate(fbare(N))
    if(.not.allocated(fvegl)) allocate(fvegl(N))
    if(.not.allocated(fvegh)) allocate(fvegh(N))
    if(.not.allocated(cvegl)) allocate(cvegl(N))
    if(.not.allocated(cvegh)) allocate(cvegh(N))
    if(.not.allocated(sncov)) allocate(sncov(N))

    cvegh = 1.0_JPRM
    cvegl = 1.0_JPRM
    select case(CIDVEG)
    case ('Ecoclimap')
       ! low vegetation fraction [--]
       if(allocated(fncfield)) deallocate(fncfield)
       allocate(fncfield(JJ,1))
       fncfield = reshape(LS%CVL,(/JJ,1/))
       do I=1,N-JJ+1,JJ ! filling up fZ with NTIME times
          fvegl(I:) = fncfield(:,1)
       end do
       if(SHOWINFO) print*,' sub out of low veg'
       ! high vegetation fraction [--]
       if(allocated(fncfield)) deallocate(fncfield)
       allocate(fncfield(JJ,1))
       fncfield = reshape(LS%CVH,(/JJ,1/))
       do I=1,N-JJ+1,JJ ! filling up fZ with NTIME times
          fvegh(I:) = fncfield(:,1)
       end do
       if(SHOWINFO) print*,' sub out of hig veg'
       ! low vegetation types
       if(allocated(fncfield)) deallocate(fncfield)
       allocate(fncfield(JJ,1))
       fncfield = reshape(LS%TVL,(/JJ,1/))
       do I=1,N-JJ+1,JJ ! filling up fZ with NTIME times
          fTVL(I:) = fncfield(:,1)
       end do
       if(SHOWINFO) print*,' sub out of low veg type'
       ! high vegetation types
       if(allocated(fncfield)) deallocate(fncfield)
       allocate(fncfield(JJ,1))
       fncfield = reshape(LS%TVH,(/JJ,1/))
       do I=1,N-JJ+1,JJ ! filling up fZ with NTIME times
          fTVH(I:) = fncfield(:,1)
       end do
       if(SHOWINFO) print*,' sub out of high veg type'
       ! water fraction [--]
       if(allocated(fncfield)) deallocate(fncfield)
       allocate(fncfield(JJ,1))
       fncfield = reshape(LS%WATER,(/JJ,1/))
       do I=1,N-JJ+1,JJ ! filling up fZ with NTIME times
          fwater(I:) = fncfield(:,1)
       end do
       if(SHOWINFO) print*,' sub out of water frac'
       ! LAI of low vegetation
       JJ = NLONS*NLATS*NTIMES
       if(allocated(fncfield)) deallocate(fncfield)
       allocate(fncfield(JJ,1))
       fncfield = reshape(LS%LAIL,(/JJ,1/))
       fs_laiL = fncfield(:,1)
       if(SHOWINFO) print*,' sub out of laiL'
       ! clean up vegetation tables
       where(fTVL.gt.7_JPIM.and.fwater.lt.0.5_JPRM) fTVL = 4_JPIM ! default low vegetation is grass
       where(fTVL.gt.7_JPIM.and.fwater.ge.0.5_JPRM) fTVL = 0_JPIM ! no vegetation on water pixel
       where(fTVH.gt.7_JPIM.and.fwater.lt.0.5_JPRM) fTVH = 1_JPIM ! default high vegetation in grass
       where(fTVH.gt.7_JPIM.and.fwater.ge.0.5_JPRM) fTVH = 0_JPIM ! no vegetation on water pixel

    case ('Tessel','HTessel')
       ! low vegetation fraction [--] (Tessel)
       if(allocated(fncfield)) deallocate(fncfield)
       allocate(fncfield(JJ,1))
       fncfield = reshape(LS%CVH,(/JJ,1/))
       do I=1,N-JJ+1,JJ ! filling up fZ with NTIME times
          fvegl(I:) = fncfield(:,1)
       end do
       ! high vegetation fraction [--] (Tessel)
       if(allocated(fncfield)) deallocate(fncfield)
       allocate(fncfield(JJ,1))
       fncfield = reshape(LS%CVL,(/JJ,1/))
       do I=1,N-JJ+1,JJ ! filling up fZ with NTIME times
          fvegh(I:) = fncfield(:,1)
       end do
       ! low vegetation type (Tessel)
       if(allocated(fncfield)) deallocate(fncfield)
       allocate(fncfield(JJ,1))
       fncfield = reshape(LS%TVL,(/JJ,1/))
       do I=1,N-JJ+1,JJ ! filling up fZ with NTIME times
          fTVL(I:) = fncfield(:,1)
       end do
       ! high vegetation types (Tessel)
       if(allocated(fncfield)) deallocate(fncfield)
       allocate(fncfield(JJ,1))
       fncfield = reshape(LS%TVH,(/JJ,1/))
       do I=1,N-JJ+1,JJ ! filling up fZ with NTIME times
          fTVH(I:) = fncfield(:,1)
       end do
       ! water fraction (from LSM in TESSEL) [--]
       if(allocated(fncfield)) deallocate(fncfield)
       allocate(fncfield(JJ,1))
       fncfield = reshape(LS%LSM,(/JJ,1/))
       do I=1,N-JJ+1,JJ ! filling up fZ with NTIME times
          fwater(I:) = 1.0_JPIM - fncfield(:,1)
       end do
       ! LAI for low vegetation (Tessel)
       if(CIDVEG.eq.'Tessel') then
          call vegtable(RVLAI,RVCOV,RVTV,b,VWC,b1,b2,b3,ZNrh,ZNrv,Ztth,Zttv,Zhr,Zw_eff)
          fs_laiL = rvlai(fTVL) ! ???
       else
          if(allocated(fncfield)) deallocate(fncfield)
          allocate(fncfield(JJ,1))
          fncfield = reshape(LS%LAIL,(/JJ,1/))  !?? check correct var
          do I=1,N-JJ+1,JJ ! filling up fZ with NTIME times
             fs_laiL(I:) = fncfield(:,1)
          end do
       end if
    case default
       print*,'WARNING: vegetation specification not valid!!! :', CIDVEG
    end select
    if(SHOWINFO) print*, 'out of fl_laiL'
    ! -- end for VEGETATION variables ----------------------

    ! 7.1 -- Calculate tile fraction -----
    !     -- TESSEL vegetation map needs to account for bare soil

    if(CIDVEG.eq.'Tessel') then
       cvegl = rvcov(fTVL)
       cvegh = rvcov(fTVH)
       ! convert model vegetation types to CMEM vegetation requirements
       fTVL = RVTV(fTVL)
       fTVH = RVTV(fTVH)
    end if
    fvegl = fvegl*cvegl
    fvegh = fvegh*cvegh
    fbare = 1.0 - fvegl - fvegh
    !    -- snow cover on land tiles
    where(fsnowd.gt.0.0_JPRM) sncov = 1.0_JPRM
    where(fsnowd.le.0.0_JPRM) sncov = 0.0_JPRM

    ftfrac(:,1) = fbare * (1.0-fwater) * (1.0-sncov)  ! bare soil
    ftfrac(:,2) = fbare * (1.0-fwater) * (sncov)      ! bare soil with snow
    ftfrac(:,3) = fvegl * (1.0-fwater) * (1.0-sncov)  ! low vegetation
    ftfrac(:,4) = fvegl * (1.0-fwater) * (sncov)      ! low vegetation with snow
    ftfrac(:,5) = fvegh * (1.0-fwater) * (1.0-sncov)  ! high vegetation
    ftfrac(:,6) = fvegh * (1.0-fwater) * (sncov)      ! high vegetation with snow
    ftfrac(:,7) = fwater                              ! water

    if(LOMASK_OCEAN) where(fwater.ge.0.5) mask = 0.0_JPRM
    ! DEALLOCATE TEMPORARY VARIABES
    if(allocated(fbare)) DEALLOCATE (fbare)
    if(allocated(fvegl)) DEALLOCATE (fvegl)
    if(allocated(fvegh)) DEALLOCATE (fvegh)
    if(allocated(cvegl)) DEALLOCATE (cvegl)
    if(allocated(cvegh)) DEALLOCATE (cvegh)
    if(allocated(sncov)) DEALLOCATE (sncov)
    if(SHOWINFO) print*,' out of vegi'
    ! ======= end VEGETATION tails ===================

    ! 8- 2T variable: -- Vegetation temperature [K] --
    if(.not.allocated(ftveg)) allocate(ftveg(N))
    select case (CITVEG)
    case ('Tsurf')
       ftveg = ftl_lsm(:,1)   ! (ftl_lsm(N,nlay_soil_ls))
    case ('Tair')
       ! RD: Temperature air 2m [K]
       JJ = NLONS*NLATS*NTIMES
       if(.not.allocated(fncfield)) allocate(fncfield(JJ,1))
       fncfield = reshape(LS%TAIR,(/JJ,1/))
       ftveg = fncfield(:,1)
       ! Update the mask to eliminate unrealistic Temperatures
       if(LOMASK_AUTO) where(ftveg.lt.100.0_JPRM) mask = 0_JPIM
    end select
    if(SHOWINFO) print*,'out of veg temperature'
    ! 8.1- Air temperature [K] --
    select case(CIATM)
    case('Pellarin','Ulaby')
       if(.not.allocated(ftair)) allocate(ftair(N))
       ! RD: T air 2m [K]
       JJ = NLONS*NLATS*NTIMES
       if(.not.allocated(fncfield)) allocate(fncfield(JJ,1))
       fncfield = reshape(LS%TAIR,(/JJ,1/))
       ftair = fncfield(:,1)
       ! Update the mask to eliminate unrealistic Temperatures
       if(LOMASK_AUTO) where(ftair.lt.100.0_JPRM) mask = 0_JPIM
    end select
    if(SHOWINFO) print*,'out of air temperature'
    ! ===== end TAIR variable ==================================

    ! 9- sal_sea variable: -- Read SSS
    if(.not.allocated(sal_sea)) allocate(sal_sea(N))
    sal_sea = 32.5_JPRM

    ! 10- theta_inc variable: -- Local Incidence angle for every pixel
    JJ = NLONS*NLATS
    if(.not.allocated(ftheta_inc)) allocate(ftheta_inc(N,NINC))
    if(allocated(fncfield)) deallocate(fncfield)
    allocate(fncfield(JJ,NINC))
    print*, 'memory: shape of theta_inc: ', shape(LS%theta_inc)
    fncfield = reshape(LS%theta_inc,(/JJ,NINC/))
    do I=1,N-JJ+1,JJ ! filling up fZ with NTIME times
       ftheta_inc(I:,:) = fncfield !(:,1)
    end do
    if(SHOWINFO) print*, 'out of ftheta_inc!'
    ! ===== end variable fZ ====================


    deallocate(fncfield, var_in_tr)

    return
  end subroutine memory_cmem_forcing
  ! ---- END SUBROUTINE TO SET CMEM VARIABLES TO MEMORY ----------
  ! ==============================================================


end module clm4cmem
! ---------------------------------------------------------------
! ************* End of Module Definition ************************
! ---------------------------------------------------------------
