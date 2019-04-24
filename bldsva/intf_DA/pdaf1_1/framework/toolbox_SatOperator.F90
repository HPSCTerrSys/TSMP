module SatOperator
  use clm4cmem, only : SATELLITE, SHOWINFO,CLM_DATA
  use constants

  contains
  ! ===============================================================
  ! SUBROUTINE TO APPLY THE SATELLITE OPERATOR on Hi-res TB_HV
  ! ---------------------------------------------------------------
  ! Inputs:
  ! * SAT       : SATELLITE type with sensor information
  ! * lons, lats: 1D CLM longitude and latitude
  ! * resol_km  : 1D CLM grid resolution [km] (lon_resol, lat_resol)
  ! * TB        : 4D matrix with brightness temperature (lon,lat,pol,inc)
  ! Outputs:
  ! * SAT%lon_foprt : 1D longitude for satellite Footprint [deg] 
  ! * SAT%lat_foprt : 1D latitude for satellite Footprint  [deg]
  ! * SAT%incl_foprt: 1D inclination angle for Footprint   [deg]
  ! * SAT%TBSAT_HV  : 3D matrix with TBs (N_lonlat,pol,inc)
  !
  ! (c) 2017 P. Saavedra Garfias (pablosaa@uni-bonn.de) UNI BONN
  ! see: LICENCE.TXT
  ! --------------------------------------------------------------
  subroutine TBSAT_OPERATOR(SAT,CLMVARS, TB)

    implicit none

    type(SATELLITE), intent(inout) :: SAT
    type(CLM_DATA),  intent(in) :: CLMVARS
    real, allocatable, dimension(:) :: lons, lats, resol_km
    real, TARGET, intent(in), dimension(:,:,:,:,:) :: TB
    ! -- internal variables 
    integer :: NLONS, NLATS, NINC, NTIME
    integer :: Ntot, UFILE, I, J, IINC, IT, nx, ny, ij_org(2)
    integer, allocatable, dimension(:) :: xmin, ymin, idxlon, idxlat
    integer :: lon1, lon2, lat1, lat2
    real, target, allocatable, dimension(:,:) :: Fnn, theta_nn
    real, allocatable, dimension(:,:,:) :: Fxy, theta_xy
    real, allocatable, dimension(:) :: diff
    real :: eps_LON, eps_LAT, zetha
    real, pointer, dimension(:,:) :: Wnn, Axy
    logical :: ROTATE = .false.
    
    lons  = CLMVARS%lons
    lats  = CLMVARS%lats
    resol_km = CLMVARS%resol_km
    NINC  = size(SAT%theta)
    NLONS = size(lons)
    NLATS = size(lats)
    NTIME = size(SAT%time)
    Ntot  = size(SAT%lon_foprt)
    ! average angular resolutions for longitude and latitude:
    ! (assuming that CLM has always equally spaced grid)
    eps_LON = sum(lons(2:)-lons(1:NLONS-1))/NLONS
    eps_LAT = sum(lats(2:)-lats(1:NLATS-1))/NLATS

    ! Loading the (lon,lat) data from the satellite footprints:
!!$    open(unit=UFILE,file=SAT%OrbitFileName,status='old')
!!$    read(UFILE,'(I5)') Ntot
!!$
!!$    if(.not.allocated(SAT%lon_foprt)) allocate(SAT%lon_foprt(Ntot))
!!$    if(.not.allocated(SAT%lat_foprt)) allocate(SAT%lat_foprt(Ntot))
!!$    if(.not.allocated(SAT%incl_foprt)) allocate(SAT%incl_foprt(Ntot))

!!$
!!$    read(UFILE,*) (SAT%lon_foprt(I), SAT%lat_foprt(I),&
!!$         SAT%incl_foprt(I), I=1,Ntot)
!!$    close(UFILE)

    if(.not.allocated(SAT%TBSAT_HV)) allocate(SAT%TBSAT_HV(Ntot,2,NINC,NTIME))
    call AntennaPattern_2D(SAT,resol_km,Fxy,theta_xy)

    nx = size(Fxy,1)
    ny = size(Fxy,2)
    ij_org(:) = int((/nx/2., ny/2./))   ! Origin for rotation

    ! Allocating local variables:
    allocate(xmin(Ntot))
    allocate(ymin(Ntot))
    allocate(diff(Ntot))
    allocate(idxlon(nx),idxlat(ny))
    allocate(Fnn(nx,ny), theta_nn(nx,ny))

    if(SHOWINFO) print*, '% reading orbit file ',Ntot, ' EPS_lon=',eps_LON, ' EPS_lat=',eps_LAT

    !!print*, 'minloc=['
    if(associated(Axy)) nullify(Wnn)
    if(associated(Wnn)) nullify(Wnn)
    do I=1,Ntot
       diff = abs(SAT%lon_foprt(I) - lons)
       xmin(I) = minloc(diff, dim=1, mask=diff.lt.eps_LON)
       diff = abs(SAT%lat_foprt(I) - lats)
       ymin(I) = minloc(diff, dim=1, mask=diff.lt.eps_LAT)
       ! Applying the Antenna Filter to satellite pixels:
       if(xmin(I).eq.0.or.ymin(I).eq.0) then
          ! no nearby pixel found, assigning NaN
          SAT%TBSAT_HV(I,:,:,:) = -99.  ! should be NaN 
          cycle
       end if

       ! creating index vectors to extract TB sub-array to apply
       ! antenna pattern filter:
       idxlon = 0
       idxlat = 0
       idxlon = (/(int(xmin(I)-(floor(nx/2.)-J)),J=0,nx-1)/)
       idxlat = (/(int(ymin(I)-(floor(ny/2.)-J)),J=0,ny-1)/)
       ! selecting only indeces inside the grid:
       lon1 = minloc(idxlon, dim=1, mask=idxlon.gt.0)
       lon2 = maxloc(idxlon, dim=1, mask=idxlon.le.NLONS)
       lat1 = minloc(idxlat, dim=1, mask=idxlat.gt.0)
       lat2 = maxloc(idxlat, dim=1, mask=idxlat.le.NLATS)

       ! Rotation for F00 and theta_00:
       ! Rotational Matrix:
       if(SAT%incl_foprt(I).gt.0.0.and.SAT%incl_foprt(I).lt.180.0) ROTATE=.true.

       if(I.gt.1.and.abs(SAT%incl_foprt(I)-SAT%incl_foprt(I-1)).le.5) ROTATE = .false.
       do IINC=1, NINC
          if(ROTATE) then
             zetha = SAT%incl_foprt(I)-90.0  ! [deg] rotation matrix
             Fnn = rotate2D_ccw(Fxy(:,:,IINC),zetha,pivot=ij_org)
             theta_nn = rotate2D_ccw(theta_xy(:,:,IINC),zetha,pivot=ij_org)
          else
             Fnn = Fxy(:,:,IINC)
             theta_nn = theta_xy(:,:,IINC)
          end if

          ! -------------------------------
          ! WARNING: this is temporal....
          if(I==-99.and.IINC.eq.1) then
          write(*,'("Fnn",I1,"=[")') IINC
          do J=1,nx
             print*, Fnn(J,:),';'
          end do
          print*, '];'
          write(*,'("THETA",I1,"=[")') IINC
          do J=1,nx
             print*, theta_nn(J,:),';'
          end do
          print*, '];'
          end if
          ! WARNING here ends the temporal....
          ! ---------------------------------

          ! Convolution between TB-hires and Antenna Pattern:
          Wnn => Fnn(lon1:lon2,lat1:lat2)
          do IT=1,NTIME
             ! for H-pol:
             Axy => TB(idxlon(lon1):idxlon(lon2),&
                  idxlat(lat1):idxlat(lat2),1,IINC,IT)

             SAT%TBSAT_HV(I,1,IINC,IT) = &
                  sum(sum(Axy*Wnn,dim=2),dim=1)/sum(sum(Wnn,dim=2),dim=1)
             ! for V-pol:
             nullify(Axy)
             Axy => TB(idxlon(lon1):idxlon(lon2),&
                  idxlat(lat1):idxlat(lat2),2,IINC,IT)
             SAT%TBSAT_HV(I,2,IINC,IT) = &
                  sum(sum(Axy*Wnn,dim=2),dim=1)/sum(sum(Wnn,dim=2),dim=1)
             nullify(Axy)
          end do  ! end over NTIME (time domain from CLM)
       end do  ! end over NINC (number of incidence angles)
    end do  ! end over Ntot (nmber of satellite lon,lat pixels)

    if(associated(Axy).or.associated(Wnn)) nullify(Axy,Wnn)
    ! deallocating local variables:
    deallocate(xmin,ymin,diff,idxlon,idxlat,Fnn,theta_nn,Fxy,theta_xy)
    
    return
  end subroutine TBSAT_OPERATOR
  ! ---------------- end of Satellite Operator -------------------

  ! ================================================================
  ! SUBROUTINE TO ESTIMATE THE ANTENNA PATTERN FOR THE
  ! SPECIFIC SATELLITE'S SENSOR CHARACTERISTICS
  ! ---------------------------------------------------------------
  ! INPUT:
  ! * SAT -> SATELLITE structure with sensor information
  ! * resol_km -> grid resolution in [km] (longiture, latiture)
  ! OUTPUT:
  ! * Fnn  -> Antenna pattern (NxM) for every incidence angle
  ! * theta_nn -> (optional) distance from origin for Fnn grid [km]
  ! Reference:
  !     Microwave Remote Sensing, Active & Passive.
  !     Vol. I, Microwave Remote Sensing Fundamentals and Radiometry
  !
  ! (c) 2016 P. Saavedra Garfias (pablosaa@uni-bonn.de) UNI BONN
  ! see: LICENCE.TXT
  ! ----------------------------------------------------------------
  subroutine AntennaPattern_2D(SAT, resol_km, Fnn, theta_nn)

    implicit none
    type(SATELLITE), intent(in) :: SAT
    real, dimension(2), intent(in) :: resol_km
    real, allocatable, dimension(:,:,:), intent(out) :: Fnn
    real, allocatable, dimension(:,:,:), intent(out), optional :: theta_nn
    !-- Internal variables:
    integer :: NXpix, NYpix, Npix, Nang
    integer :: I, J, K
    integer :: ij_org(2) ! temporal
    real :: Horb, theta_FWHM, lambda_SAT
    real :: D_antenna
    real, dimension(size(SAT%theta)) :: theta_inc, theta_hat,&
         r_orb, FoPrint_a, FoPrint_b
    real :: d_l, alpha_FWHM
    real, allocatable, dimension(:) :: rr_a, rr_b
    real, allocatable, dimension(:,:) :: theta_a, theta_b, Rxy
    real, allocatable, dimension(:,:) :: Fij, thetaij


    ! Defining satellite's sensor features:
    Horb = SAT%orbit   ! [km]
    D_antenna = SAT%antenna   ! [m]
    lambda_SAT = SAT%wavelength ! [m] e.g. 0.21 for L-band
    theta_inc = DEG2rad*SAT%theta  ! satellite incidence angles [rad]
    alpha_FWHM = DEG2rad*(70.0*lambda_SAT/D_antenna) ! [rad] FOV at -3dB
    d_l = 1./alpha_FWHM     ! 1/[rad], altenative d_l=1/FOV
    Nang = size(theta_inc)
    theta_hat = asin((1+Horb/Ro)*sin(theta_inc)) ! local inc angle [rad]
    r_orb = Ro*(sin(theta_hat)/tan(theta_inc) - cos(theta_hat))  ! [km]
    
    ! Calculating FootPrint over Earth surface
    ! with axis cross-track (_a) and along-track (_b):
    FoPrint_a = 2.0*r_orb*tan(alpha_FWHM/2.0)   ! [km]
    FoPrint_b = 2.0*r_orb*cos(theta_hat)/tan(alpha_FWHM/2.0)/&
         ((cos(theta_hat)/tan(alpha_FWHM/2.0))**2-sin(theta_hat)**2)

    ! Number of pixels to cover the FOV_FWHM dimensions over Earth:
    NXpix = int(maxval(FoPrint_a)/resol_km(1))
    NYpix = int(maxval(FoPrint_b)/resol_km(2))
    Npix  = 150 !2*max(NXpix,NYpix)  ! for a Fnn square matrix
    ij_org(:) = int((/Npix/2., Npix/2./))   ! Temporal
    
    if(.not.allocated(rr_a)) allocate(rr_a(Npix))
    if(.not.allocated(rr_b)) allocate(rr_b(Npix))
    if(.not.allocated(Rxy)) allocate(Rxy(Npix,Npix))
    if(.not.allocated(theta_a)) allocate(theta_a(Npix,Npix))
    if(.not.allocated(theta_b)) allocate(theta_b(Npix,Npix))
    if(.not.allocated(Fij)) allocate(Fij(Npix,Npix))
    if(.not.allocated(thetaij)) allocate(thetaij(Npix,Npix))

    ! following variables will be deallocated out-side
    if(.not.allocated(Fnn)) allocate(Fnn(Npix,Npix,Nang))
    if(.not.allocated(theta_nn)) allocate(theta_nn(Npix,Npix,Nang))
    
    do I=1,Nang
       ! Eliptic Footprint grid-bins for cross-track semi-axes a [rad]
       rr_a = linspace(-resol_km(1)*Npix/2,resol_km(1)*Npix/2,Npix)
       rr_a = atan(rr_a/r_orb(I))
       ! Eliptic Footprint grid-bins for along-track semi-axes b [rad]
       rr_b = linspace(-resol_km(2)*Npix/2,resol_km(2)*Npix/2,Npix)
       where(rr_b.ge.0.0) 
          rr_b = atan(rr_b*cos(theta_hat(I))/&
               (r_orb(I)+rr_b*sin(theta_hat(I))))
       elsewhere
          rr_b = atan(rr_b*cos(theta_hat(I))/&
               (r_orb(I)-rr_b*sin(theta_hat(I))))
       end where
       ! Creating grid-cells:
       call meshgrid(rr_a, rr_b,theta_a, theta_b)

       Rxy = sqrt(sin(theta_a)**2 + sin(theta_b)**2)+1.E-20
       ! Antenna power
       Fnn(:,:,I) = abs(sin(pPI*d_l*Rxy)/(pPI*d_l*Rxy))
       if(present(theta_nn)) theta_nn(:,:,I) = asin(Rxy)*r_orb(I)/cos(theta_hat(I))
       !*Rxy/cos(asin(Rxy)+theta_hat(I))
       
       Fnn(ij_org(1),ij_org(1):,I) = 1.E-5 ! ONLY testing, to take out

    end do
    
    ! deallocating local variables:
    deallocate(rr_a, rr_b, theta_a, theta_b, Rxy) 
    return
  end subroutine AntennaPattern_2D
  ! ---------- END OF ANTENNA PATTERN SUBTOUTINE -----------------

  ! =======================================================
  ! FUNCTION TO ROTATE 2D MATRIX COUNTER-CLOCKWISE
  ! -------------------------------------------------------
  ! Where:
  ! * A        : NxM original matrix to rotate
  ! * angle_deg: angle from abcise to rotate counter-clockwise
  ! * pivot    : (optional) indexes [i0,j0] as origin to rotate around
  ! * FillVal  : (optional) value to fill X where elements is assigned
  ! Output:
  ! * X        : NxM rotated matrix X = rotational(A)
  !
  ! (c) 2016 P. Saavedra Garfias (pablosaa@uni-bonn.de) UNI BONN
  ! see: LICENCE.TXT
  ! -------------------------------------------------------------

  function rotate2D_ccw(A,angle_deg,pivot,FillVal) result(X)
    implicit none
    real, intent(in), dimension(:,:) :: A
    real, intent(in) :: angle_deg
    integer, intent(in), optional, dimension(2)   :: pivot
    real, intent(in), optional :: FillVal
    real, dimension(size(A,1),size(A,2)) :: X
    ! -- internal variables
    integer :: K, J, Nx, Ny
    integer :: nm(2), kj_org(2)
    real, dimension(2,2) :: ROTA
    real, parameter :: deg2rad = 4.*atan(1.)/180.

    ! assigning variables:
    X = -99.  ! default fill Value
    Nx=size(A,1)
    Ny=size(A,2)
    kj_org = 0
    if(present(pivot)) kj_org = pivot
    if(present(FillVal)) X = FillVal
    ROTA(1,:) = (/ cos(angle_deg*deg2rad), -sin(angle_deg*deg2rad) /)
    ROTA(2,:) = (/ sin(angle_deg*deg2rad),  cos(angle_deg*deg2rad) /)
    do K=1,Nx
       do J=1,Ny
          nm = nint(matmul(((/K,J/) - kj_org),ROTA)) + kj_org
          if(nm(1)>0.and.nm(1)<=Nx.and.nm(2)>0.and.nm(2)<=Ny) then
             X(K,J) = A(nm(1),nm(2))
          end if
       end do
    end do

    return
  end function rotate2D_ccw
  ! --------------- end of rotate matric function ----------------

  ! =============================================================
  ! FUNCTION TO MAKE A VECTOR LINEARLY SPACED BETWEEN TWO POINTS
  ! INPUT:
  ! * x and y: scalars with x minimum and y maximum interval
  ! * n : number of elements for result vector from x to y
  ! OUTPUT:
  ! * VAR: 1D vector with n-elements from [x...y]
  !
  ! (c) 2017 P. Saavedra Garfias (pablosaa@uni-bonn.de) UNI BONN
  ! see: LICENCE.TXT
  ! -------------------------------------------------------------  
  function linspace(x,y,n) result(var)
    implicit none
    real, intent(in) :: x,y
    integer, intent(in) :: n
    real, dimension(n) :: var
    ! -- internal variables                                              
    integer :: I
    real :: delVAR
    
    delVAR = (y-x)/(n-1)
    var(1) = x
    var(n) = y
    var(2:n-1) = var(1) + (/((I-1)*delVAR, I = 2,n-1)/)
    return                                                               
  end function linspace
  ! ------------ end of linspace() function ---------------------

  ! =============================================================
  ! SUBROUTINE TO MAKE A 2d MESH-GRID FROM TWO VECTORS
  ! INPUT:
  ! * xx, yy: vectors with M and N elements respectively
  ! OUTPUT:
  ! * X, Y: 2D matrices with MxN elements both
  !
  ! (c) 2017 P. Saavedra Garfias (pablosaa@uni-bonn.de) UNI BONN
  ! see: LICENCE.TXT
  ! -------------------------------------------------------------
  subroutine meshgrid(xx,yy,X,Y)
    implicit none
    real, dimension(:), intent(in) :: xx, yy
    real, dimension(:,:),intent(out) :: X, Y
    ! -- internal variables                                               
    integer :: nx, ny
    nx = size(xx)
    ny = size(yy)
    X = reshape(spread(xx,2,ny),(/nx, ny/))
    Y = reshape(spread(yy,1,nx),(/nx, ny/))
    return
  end subroutine meshgrid
  ! ----------- end of function meshgrid -------------------------

end module SatOperator
