  SUBROUTINE pfread_var(fname,nx,ny,nz,pfvar)

     IMPLICIT NONE
     !
     CHARACTER*(*),   INTENT(IN)        :: fname
     INTEGER(KIND=4), INTENT(IN)        :: nx,ny,nz
     REAL(KIND=8),    INTENT(OUT)       :: pfvar(nx,ny,nz)
     REAL(KIND=8)                       :: dx, dy, dz, x1, y1, z1
     INTEGER(KIND=4)                    :: cps
     INTEGER(KIND=4)                    :: i,j,k, ix, iy, iz, is,       &
                                           ns,  rx, ry, rz,nnx, nny, nnz
! 
! 
  OPEN(100,file=trim(fname),form='unformatted' , & 
         access='stream',convert='BIG_ENDIAN',status='old', &
         position='REWIND')         !binary outputfile of Parflow, gfortran
!        recordtype='stream',convert='BIG_ENDIAN',status='old')     !binary outputfile of Parflow, ifort
   
! Read in header info
! Start: READing of domain spatial information
  READ(100) x1 !X
  READ(100) y1 !Y 
  READ(100) z1 !Z

  READ(100) cps !NX
  READ(100) cps !NY
  READ(100) cps !NZ
 
  READ(100) dx !DX
  READ(100) dy !DY
  READ(100) dz !DZ

  READ(100) ns !num_subgrids
! End: READing of domain spatial information

! Start: loop over number of sub grids
  do is = 0, ns-1 

! Start: READing of sub-grid spatial information
   READ(100) ix
   READ(100) iy
   READ(100) iz

   READ(100) nnx
   READ(100) nny
   READ(100) nnz

   READ(100) rx
   READ(100) ry
   READ(100) rz

! End: Reading of sub-grid spatial information

! Start: Read in data from each individual subgrid
  DO  k=iz +1 , iz + nnz
   DO  j=iy +1 , iy + nny
    DO  i=ix +1 , ix + nnx
     READ(100) pfvar(i,j,k)
    END DO
   END DO
  END DO
! End: Read in data from each individual subgrid

  END DO
! End: loop over number of sub grids

  CLOSE(100)
  RETURN
 END SUBROUTINE pfread_var 
