!> \par Purpose
!>
!> Compute the trasformation matrix \c DZDV \f$ = \frac{\partial Z}{\partial \tilde{U}} \f$
!> from parameter vector to the set of Euler symmetrizing variables
!> \f[
!>   \frac{\partial Z}{\partial \tilde{U}} = \frac{1}{\sqrt{\rho}} \left[
!>  \begin{array}{ccc}
!> \frac{1}{2} & \frac{1}{2} \rho / a & \mathbf{0} \\\
!>\frac{1}{2} \left( 2 k - H \right) & \frac{\rho}{a}
!> \left[ \left(\gamma-\frac{1}{2}\right) H
!> -\left(\gamma-1\right) k \right] & \rho \mathbf{u}^t \\\
!> \frac{1}{2} \mathbf{u} & \frac{1}{2} \rho \mathbf{u} / a  & \rho \, I_{d\times d}
!> \end{array} \right] \\\ 
!> = \left[
!> \begin{array}{ccccc}
!> \frac{1}{2 z_1} & \frac{z_1}{2 \,a} & 0 & 0 & 0 \\
!> \frac{1}{2 z_1^2} \left( \left( z_3^2+z_4^2+z_5^2 \right)/z_1 - z_2
!> \right) & \frac{1}{a} \left[ \left(\gamma-\frac{1}{2}\right) z_2
!> -\frac{\gamma-1 }{2} \left( z_3^2+z_4^2+z_5^2 \right)/z_1 \right] &
!> z_3 & z_4 & z_5 \\
!> \frac{z_3}{2 \, z_1^2} & \frac{1}{2a} z_3 & z_1 & 0 & 0 \\
!> \frac{z_4}{2 \, z_1^2} & \frac{1}{2a} z_4 & 0 & z_1 & 0 \\
!> \frac{z_5}{2 \, z_1^2} & \frac{1}{2a} z_5 & 0 & 0 & z_1
!> \end{array} \right]
!> \f]
!
!> @param[in]  ZROE Roe's parameter vector
!> @param[out] DZDV is the matrix \f$ \frac{\partial Z}{\partial \tilde{U}} \f$
!> @param[in]  NOFVAR leading dimension of DZDV
!> @param[in]  NDIM dimension of the space
!
!> \author $Author: abonfi $
!> \version $Revision: 1.5 $
!> \date $Date: 2013/08/22 11:15:08 $
!
      SUBROUTINE SYMM2PARM(ZROE,DZDV,NOFVAR,NDIM)
C
C     $Id: symm2parm.f,v 1.5 2013/08/22 11:15:08 abonfi Exp $
C
      IMPLICIT NONE
C
C     transformation matrix from symmetrizing to
C     parameter vector
C
C     .. Parameters ..
      INCLUDE 'constants.h'
      INCLUDE 'pfcgas.com'
C     ..
C     .. Scalar Arguments ..
      INTEGER NDIM,NOFVAR
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION DZDV(NOFVAR,NOFVAR),ZROE(NOFVAR)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION AINV,ASQR,TEMP,ZINV,ZINVSQR
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC SQRT
C     ..
      TEMP = ZROE(3)*ZROE(3) + ZROE(4)*ZROE(4)
      IF (NDIM.EQ.3) TEMP = TEMP + ZROE(5)*ZROE(5)
      TEMP = HALF*TEMP
C
      ZINV = ONE/ZROE(1)
      ZINVSQR = ZINV*ZINV
      ASQR = GM1* (ZROE(2)-TEMP*ZINV)*ZINV
      AINV = ONE/SQRT(ASQR)
C
C     Transformation matrix from parameter vector to symmetrizing variables
C
      DZDV(1,1) = HALF*ZINV
      DZDV(2,1) = HALF*ZINVSQR* (TWO*TEMP*ZINV-ZROE(2))
      DZDV(3,1) = HALF*ZROE(3)*ZINVSQR
      DZDV(4,1) = HALF*ZROE(4)*ZINVSQR
C
      DZDV(1,2) = HALF*ZROE(1)*AINV
      DZDV(2,2) = ((GAM-HALF)*ZROE(2)-GM1*TEMP*ZINV)*AINV
      DZDV(3,2) = HALF*ZROE(3)*AINV
      DZDV(4,2) = HALF*ZROE(4)*AINV
C
      DZDV(1,3) = ZERO
      DZDV(2,3) = ZROE(3)
      DZDV(3,3) = ZROE(1)
      DZDV(4,3) = ZERO
C
      DZDV(1,4) = ZERO
      DZDV(2,4) = ZROE(4)
      DZDV(3,4) = ZERO
      DZDV(4,4) = ZROE(1)
C
      IF (NDIM.EQ.2) RETURN
C
      DZDV(5,1) = HALF*ZROE(5)*ZINVSQR
      DZDV(5,2) = HALF*ZROE(5)*AINV
      DZDV(5,3) = ZERO
      DZDV(5,4) = ZERO
C
      DZDV(1,5) = ZERO
      DZDV(2,5) = ZROE(5)
      DZDV(3,5) = ZERO
      DZDV(4,5) = ZERO
      DZDV(5,5) = ZROE(1)
C
      RETURN

      END
