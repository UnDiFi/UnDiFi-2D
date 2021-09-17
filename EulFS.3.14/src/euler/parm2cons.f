!> \par Purpose
!>
!> Compute the transformation matrix from parameter vector to conserved variables
!>
!> \f[
!> \frac{\partial U}{\partial Z} = \left[
!> \begin{array}{ccccc}
!> 2 z_1 & 0 & 0 & 0 & 0 \\
!> \frac{z_2}{\gamma} & \frac{z_1}{\gamma} &
!> \frac{\gamma-1}{\gamma} z_3 & \frac{\gamma-1}{\gamma} z_4 &
!> \frac{\gamma-1}{\gamma} z_5 \\
!> z_3 & 0 & z_1 & 0 & 0 \\
!> z_4 & 0 & 0 & z_1 & 0 \\
!> z_5 & 0 & 0 & 0 & z_1
!> \end{array} \right]
!> \f]
!
!> @param[in]  ZROE Roe's parameter vector
!> @param[out]  DUDZ is the matrix \f$ \frac{\partial U}{\partial Z} \f$
!> @param[in]  NOFVAR leading dimension of DUDZ
!> @param[in]  NDIM dimension of the space
!
!> \author $Author: abonfi $
!> \version $Revision: 1.5 $
!> \date $Date: 2013/08/21 10:30:45 $
!
      SUBROUTINE PARM2CONS(ZROE,DUDZ,NOFVAR,NDIM)
C
C     $Id: parm2cons.f,v 1.5 2013/08/21 10:30:45 abonfi Exp $
C
      IMPLICIT NONE
C
C
C
C     Assembles the dUdZ matrix ...
C
C     .. Parameters ..
      INCLUDE 'constants.h'
      INCLUDE 'pfcgas.com'
C     ..
C     .. Scalar Arguments ..
      INTEGER NDIM,NOFVAR
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION DUDZ(NOFVAR,NOFVAR),ZROE(NOFVAR)
C     ..
C     .. Local Scalars ..
C     ..
      DUDZ(1,1) = TWO*ZROE(1)
      DUDZ(1,2) = ZERO
      DUDZ(1,3) = ZERO
      DUDZ(1,4) = ZERO
C
      DUDZ(2,1) = GINV*ZROE(2)
      DUDZ(2,2) = GINV*ZROE(1)
      DUDZ(2,3) = GM1OG*ZROE(3)
      DUDZ(2,4) = GM1OG*ZROE(4)
C
      DUDZ(3,1) = ZROE(3)
      DUDZ(3,2) = ZERO
      DUDZ(3,3) = ZROE(1)
      DUDZ(3,4) = ZERO
C
      DUDZ(4,1) = ZROE(4)
      DUDZ(4,2) = ZERO
      DUDZ(4,3) = ZERO
      DUDZ(4,4) = ZROE(1)
C
      IF (NDIM.EQ.2) RETURN
C
      DUDZ(1,5) = ZERO
      DUDZ(2,5) = GM1OG*ZROE(5)
      DUDZ(3,5) = ZERO
      DUDZ(4,5) = ZERO
C
      DUDZ(5,1) = ZROE(5)
      DUDZ(5,2) = ZERO
      DUDZ(5,3) = ZERO
      DUDZ(5,4) = ZERO
      DUDZ(5,5) = ZROE(1)
C
C
      RETURN

      END
