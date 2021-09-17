!> \par Purpose
!>
!> Source term for the convection equation
!> \f[
!>  u_t +   \lambda_x u_x+\lambda_y u_y = f
!> \f]
!>
!>
!> the source term is:
!> \f[ 
!> f = -\left(x\cos\delta+x\sin\delta\right)
!> \f]
!> and the exact steady solution is:
!> \f[
!>    u = \sin\left(2\pi\left(y\cos\delta-x\sin\delta\right)\right)+ x'y' -xy
!> \f]
!> where
!> \f{eqnarray*}{
!>    x' &=& x - \cos(\delta) (x\cos(\delta)+y\sin(\delta)) \\
!>    y' &=& y - \sin(\delta) (x\cos(\delta)+y\sin(\delta))
!> \f}
!>
!> @param[in] X Cartesian x coordinate
!> @param[in] Y Cartesian y coordinate
!> @param[in] DELTA is the angle \f$\delta\f$
!> @return The pointwise value of the source term \c f.
!> \author $Author: abonfi $
!> \version $Revision: 1.8 $
!> \date $Date: 2013/08/20 14:48:46 $
!>
!>
      DOUBLE PRECISION FUNCTION FUNSOU2(X,Y,DELTA)
C
C     $Id: fun2.f,v 1.2 2013/01/24 08:08:43 abonfi Exp abonfi $
C
      IMPLICIT NONE
C     .. Scalar Arguments ..
      DOUBLE PRECISION DELTA,X,Y
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC COS,SIN
C     ..
      FUNSOU2 = - (Y*COS(DELTA)+X*SIN(DELTA))

      RETURN

      END
