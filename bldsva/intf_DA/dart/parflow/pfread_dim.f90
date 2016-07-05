  SUBROUTINE pfread_dim(fname,nx,ny,nz,dx,dy,dz)

   IMPLICIT NONE
!
   CHARACTER*(*),   INTENT(IN)        :: fname
   INTEGER(KIND=4), INTENT(OUT)       :: nx,ny,nz
   REAL(KIND=8)                       :: dx, dy, dz
   REAL(KIND=8)                       :: x1, y1, z1 
!  
  OPEN(100,file=trim(fname),form='unformatted' , & 
         access='stream',convert='BIG_ENDIAN',status='old')         !binary outputfile of Parflow, gfortran
!        recordtype='stream',convert='BIG_ENDIAN',status='old')     !binary outputfile of Parflow, ifort
   
  ! Read in header info

! Start: READing of domain spatial information
  READ(100) x1 !X
  READ(100) y1 !Y 
  READ(100) z1 !Z

  READ(100) nx !NX
  READ(100) ny !NY
  READ(100) nz !NZ

  READ(100) dx !dX
  READ(100) dy !dY 
  READ(100) dz !dZ

  RETURN
 END SUBROUTINE pfread_dim 
