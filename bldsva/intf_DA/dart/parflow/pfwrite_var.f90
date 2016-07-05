  SUBROUTINE pfwrite_var(fname,nx,ny,nz,dx,dy,dz,xd,yd,zd,nxs,nys,pfvar)

   IMPLICIT NONE
!
   CHARACTER*(*),   INTENT(IN)        :: fname
   INTEGER(KIND=4), INTENT(IN)        :: nx,ny,nz,nxs,nys
   REAL(KIND=8),    INTENT(IN)        :: pfvar(nx,ny,nz)
   REAL(KIND=8),    INTENT(IN)        :: dx, dy, dz
   REAL(KIND=8),    INTENT(IN)        :: xd, yd, zd
   INTEGER(KIND=4)                    :: i,j,k, ix, iy, iz, is, ns,       &
                                         rx, ry, rz, nnx, nny, nnz,       &
                                         ixs, iys
! 
! 
  OPEN(100,file=trim(fname),status='new',access='stream',convert='BIG_ENDIAN',form='unformatted') 
         !binary outputfile of Parflow, gfortran
!        recordtype='stream',convert='BIG_ENDIAN',status='old')     !binary outputfile of Parflow, ifort

! Start: Writing of domain spatial information
  WRITE(100) xd !X
  WRITE(100) yd !Y 
  WRITE(100) zd !Z

  WRITE(100) nx !NX
  WRITE(100) ny !NY
  WRITE(100) nz !NZ
 
  WRITE(100) dx !DX
  WRITE(100) dy !DY
  WRITE(100) dz !DZ

  ns = INT(nxs*nys)
  WRITE(100) ns !num_subgrids
! End: Writing of domain spatial information

! Start: loop over number of sub grids
  nnx = INT(nx/nxs)
  nny = INT(ny/nys)
  nnz = nz
  iz = 0
!
  do iys = 0, nys-1 
  do ixs = 0, nxs-1
  print *,"Doing subgrid",ixs,iys

! Start: Writing of sub-grid spatial information

  ix = INT(nnx*ixs)
  iy = INT(nny*iys)
 
  WRITE(100) ix
  WRITE(100) iy
  WRITE(100) iz

  WRITE(100) nnx
  WRITE(100) nny
  WRITE(100) nnz

  rx = 0; ry = 0; rz=0
  WRITE(100) rx
  WRITE(100) ry
  WRITE(100) rz

! End: Writing of sub-grid spatial information

! Start: Write in data from each individual subgrid
  DO  k=iz +1 , iz + nnz
   DO  j=iy +1 , iy + nny
    DO  i=ix +1 , ix + nnx
     WRITE(100) pfvar(i,j,k) 
    END DO
   END DO
  END DO

! End: Write in data from each individual subgrid

  END DO
  END DO
! End: loop over number of sub grids

  CLOSE(100)
  RETURN
 END SUBROUTINE pfwrite_var 
