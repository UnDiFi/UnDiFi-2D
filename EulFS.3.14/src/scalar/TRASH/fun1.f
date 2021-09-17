!> \par Purpose
!>
!> Source term for the convection-diffusion equation
!> \f[
!>  u_t +   \lambda_xu_x+\lambda_yu_y - \varepsilon(u_{xx}+u_{yy}) = f
!> \f]
!>
!>  The convection speed is \c (A,B) = \f$\lambda_x\mathbf{e}_x+\lambda_y\mathbf{e}_y\f$
!>
!>  \c EPS is the diffusion coefficient \f$\varepsilon\f$ which is set using the runtime option:
!>  \verbatim
!>  -Reynolds [value]
!>  \endverbatim
!>
!> where \c value = \f$1/\varepsilon\f$.
!>
!> the source term is:
!> \f[ 
!> f = [(\lambda_x-2\varepsilon)y+(\lambda_y-2\varepsilon)x]e^{(x+y)} [(\lambda_x-2\varepsilon)+(\lambda_y-2\varepsilon)]xye^{(x+y)}
!> \f]
!> and the exact steady solution is:
!> \f[
!>    u = xy\exp(x+y)
!> \f]
!>
!> @param[in] X Cartesian x coordinate
!> @param[in] Y Cartesian y coordinate
!> @param[in] A Cartesian x component of the convection speed, i.e. \f$\lambda_x\f$
!> @param[in] B Cartesian y component of the convection speed, i.e. \f$\lambda_y\f$
!> @param[in] EPS is the diffusion coefficient \f$\varepsilon\f$.
!> @return The pointwise value of the source term \c f.
!> \author $Author: abonfi $
!> \version $Revision: 1.8 $
!> \date $Date: 2013/08/20 14:48:46 $
!>
!>
      DOUBLE PRECISION FUNCTION FUNSOU1(X,Y,A,B,EPS)
C
C     $Id: fun1.f,v 1.2 2013/01/24 08:08:43 abonfi Exp abonfi $
C
      IMPLICIT NONE
C
C
      INCLUDE 'constants.h'

C     .. Scalar Arguments ..
      DOUBLE PRECISION A,B,EPS,X,Y
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION TEMPA,TEMPB
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC EXP
C     ..
      TEMPA = A - TWO*EPS
      TEMPB = B - TWO*EPS

      FUNSOU1 = ((TEMPA*Y+TEMPB*X)+ (TEMPA+TEMPB)*X*Y)*EXP(X+Y)

      RETURN

      END
