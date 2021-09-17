!> \par Purpose
!>
!> Compute the transformation matrix \c DVDZ = \f$ \frac{\partial \tilde{U}}{\partial Z} \f$
!> from parameter vector to the set of the Euler symmetrizing variables
!> \f[
!> \frac{\partial \tilde{U}}{\partial Z} = \sqrt{\rho} \left[
!> \begin{array}{ccccc}
!> 2 - \frac{\gamma-1 }{\gamma} \frac{H}{a^2} & -\frac{\gamma-1}{\gamma} \frac{1}{a^2} &
!> \frac{\gamma-1}{\gamma} \frac{u}{a^2} & \frac{\gamma-1}{\gamma} \frac{v}{a^2} &
!> \frac{\gamma-1}{\gamma} \frac{w}{a^2} \\\
!> \frac{\gamma-1}{\gamma} \frac{H}{\rho a} & \frac{\gamma-1}{\gamma} \frac{1}{\rho a} &
!> -\frac{\gamma-1}{\gamma} \frac{u}{\rho a} & -\frac{\gamma-1}{\gamma} \frac{v}{\rho a} &
!> -\frac{\gamma-1}{\gamma} \frac{w}{\rho a} \\\
!> -\frac{u}{\rho} & 0 & \frac{1}{\rho} & 0 & 0 \\\
!> -\frac{v}{\rho} & 0 & 0 & \frac{1}{\rho} & 0 \\\
!> -\frac{w}{\rho} & 0 & 0 & 0 & \frac{1}{\rho}
!> \end{array} \right]
!> =
!> \left[
!> \begin{array}{ccccc}
!> 2 z_1 - \frac{\gamma-1 }{\gamma} \frac{z_2}{a^2} & -\frac{\gamma-1}{\gamma} \frac{z_1}{ a^2} &
!> \frac{\gamma-1}{\gamma} \frac{z_3}{a^2} & \frac{\gamma-1}{\gamma} \frac{z_4}{a^2} &
!> \frac{\gamma-1}{\gamma} \frac{z_5}{a^2} \\\
!> \frac{\gamma-1}{\gamma} \frac{z_2}{\rho a} & \frac{\gamma-1}{\gamma} \frac{z_1}{\rho a} &
!> -\frac{\gamma-1}{\gamma} \frac{z_3}{\rho a} & -\frac{\gamma-1}{\gamma} \frac{z_4}{\rho a} &
!> -\frac{\gamma-1}{\gamma} \frac{z_5}{\rho a} \\\
!> -\frac{z_3}{z_1^2} & 0 & \frac{1}{z_1} & 0 & 0 \\\
!> -\frac{z_4}{z_1^2} & 0 & 0 & \frac{1}{z_1} & 0 \\\
!> -\frac{z_5}{z_1^2} & 0 & 0 & 0 & \frac{1}{z_1}
!> \end{array} \right]
!> \f]
!
!> @param[in]  ZROE Roe's parameter vector
!> @param[out] DVDZ is the matrix \f$ \frac{\partial \tilde{U}}{\partial Z} \f$
!> @param[in]  NOFVAR leading dimension of DUDZ
!> @param[in]  NDIM dimension of the space
!
!> \author $Author: abonfi $
!> \version $Revision: 1.6 $
!> \date $Date: 2013/08/22 11:57:10 $
      SUBROUTINE PARM2SYMM(ZROE,DVDZ,NOFVAR,NDIM)
C
C     $Id: parm2symm.f,v 1.6 2013/08/22 11:57:10 abonfi Exp $
C
C
      IMPLICIT NONE
C
C     Transformation matrix from parameter vector to
C     symmetrizing variables
C
C     .. Parameters ..
      INCLUDE 'constants.h'
      INCLUDE 'pfcgas.com'
C     ..
C     .. Scalar Arguments ..
      INTEGER NDIM,NOFVAR
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION DVDZ(NOFVAR,NOFVAR),ZROE(NOFVAR)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION AINV,ASQR,RAINV,TEMP,ZINV,ZINVSQR
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC DSQRT
C     ..
      TEMP = ZROE(3)*ZROE(3) + ZROE(4)*ZROE(4)
      IF (NDIM.EQ.3) TEMP = TEMP + ZROE(5)*ZROE(5)
      TEMP = HALF*TEMP
C
      ZINV = ONE/ZROE(1)
      ZINVSQR = ZINV*ZINV
      ASQR = GM1* (ZROE(2)-TEMP*ZINV)*ZINV
      AINV = ONE/DSQRT(ASQR)
      RAINV = ZINVSQR*AINV
C
C
      DVDZ(1,1) = TWO*ZROE(1) - GM1OG*ZROE(2)/ASQR
      DVDZ(2,1) = GM1OG*RAINV*ZROE(2)
      DVDZ(3,1) = -ZROE(3)*ZINVSQR
      DVDZ(4,1) = -ZROE(4)*ZINVSQR
C
      DVDZ(2,2) = GM1OG*ZROE(1)*RAINV
      DVDZ(1,2) = -GM1OG*ZROE(1)/ASQR
      DVDZ(3,2) = ZERO
      DVDZ(4,2) = ZERO
C
      DVDZ(2,3) = -GM1OG*ZROE(3)*RAINV
      DVDZ(1,3) = GM1OG*ZROE(3)/ASQR
      DVDZ(3,3) = ONE/ZROE(1)
      DVDZ(4,3) = ZERO
C
      DVDZ(2,4) = -GM1OG*ZROE(4)*RAINV
      DVDZ(1,4) = GM1OG*ZROE(4)/ASQR
      DVDZ(3,4) = ZERO
      DVDZ(4,4) = ONE/ZROE(1)
C
      IF (NDIM.EQ.2) RETURN
C
      DVDZ(5,1) = -ZROE(5)*ZINVSQR
      DVDZ(5,2) = ZERO
      DVDZ(5,3) = ZERO
      DVDZ(5,4) = ZERO
C
      DVDZ(1,5) = GM1OG*ZROE(5)/ASQR
      DVDZ(2,5) = -GM1OG*ZROE(5)*RAINV
      DVDZ(3,5) = ZERO
      DVDZ(4,5) = ZERO
      DVDZ(5,5) = ONE/ZROE(1)
C
      RETURN

      END
