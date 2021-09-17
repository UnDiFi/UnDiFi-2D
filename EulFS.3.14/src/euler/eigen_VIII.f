      SUBROUTINE Eigen_VIII(Matrix,LDA,DVDZ,DUDV,NDIM,LDB)
C
C     Conserved (pressure,velocity) variables ..
C
      IMPLICIT NONE
C
      INCLUDE 'paramt.h'
      INCLUDE 'constants.h'
C
C
      INCLUDE 'three.com'
      INCLUDE 'chorin.com'
C
C
      INTEGER LDA,NDIM,LDB
C
C
      DOUBLE PRECISION Matrix(LDA,LDA,*),DVDZ(*),DUDV(*)
C
C
      INTEGER IDIM
C
      DOUBLE PRECISION DDOT
      EXTERNAL         DDOT
C
C     Jacobian matrix of the inviscid fluxes
C
      IDIM = 1
C
      Matrix(1,1,IDIM) = ZERO
      Matrix(1,2,IDIM) = BETA
      Matrix(1,3,IDIM) = ZERO
      Matrix(1,4,IDIM) = ZERO
C
      Matrix(2,1,IDIM) = ONE
      Matrix(2,2,IDIM) = TWO*ZAVG(2)
      Matrix(2,3,IDIM) = ZERO
      Matrix(2,4,IDIM) = ZERO
C
      Matrix(3,1,IDIM) = ZERO
      Matrix(3,2,IDIM) = ZAVG(3)
      Matrix(3,3,IDIM) = ZAVG(2)
      Matrix(3,4,IDIM) = ZERO
C
      Matrix(4,1,IDIM) = ZERO
      Matrix(4,2,IDIM) = ZAVG(4)
      Matrix(4,3,IDIM) = ZERO
      Matrix(4,4,IDIM) = ZAVG(2)
*
      IDIM = 2
C
      Matrix(1,1,IDIM) = ZERO
      Matrix(1,2,IDIM) = ZERO
      Matrix(1,3,IDIM) = BETA
      Matrix(1,4,IDIM) = ZERO
C
      Matrix(2,1,IDIM) = ZERO
      Matrix(2,2,IDIM) = ZAVG(3)
      Matrix(2,3,IDIM) = ZAVG(2)
      Matrix(2,4,IDIM) = ZERO
C
      Matrix(3,1,IDIM) = ONE
      Matrix(3,2,IDIM) = ZERO
      Matrix(3,3,IDIM) = TWO*ZAVG(3)
      Matrix(3,4,IDIM) = ZERO
C
      Matrix(4,1,IDIM) = ZERO
      Matrix(4,2,IDIM) = ZERO
      Matrix(4,3,IDIM) = ZAVG(4)
      Matrix(4,4,IDIM) = ZAVG(3)
*
      IDIM = 3
C
      Matrix(1,1,IDIM) = ZERO
      Matrix(1,2,IDIM) = ZERO
      Matrix(1,3,IDIM) = ZERO
      Matrix(1,4,IDIM) = BETA
C
      Matrix(2,1,IDIM) = ZERO
      Matrix(2,2,IDIM) = ZAVG(4)
      Matrix(2,3,IDIM) = ZERO
      Matrix(2,4,IDIM) = ZAVG(2)
C
      Matrix(3,1,IDIM) = ZERO
      Matrix(3,2,IDIM) = ZERO
      Matrix(3,3,IDIM) = ZAVG(4)
      Matrix(3,4,IDIM) = ZAVG(3)
C
      Matrix(4,1,IDIM) = ONE
      Matrix(4,2,IDIM) = ZERO
      Matrix(4,3,IDIM) = ZERO
      Matrix(4,4,IDIM) = TWO*ZAVG(4)
C
      RETURN
      END
