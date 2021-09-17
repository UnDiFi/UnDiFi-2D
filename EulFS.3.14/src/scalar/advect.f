!> \par Purpose
!>
!> This subroutine returns
!> the cell averaged advection speed \f$\hat{\lambda}\f$
!> and the integral over the cell of the source term, i.e. \f$\int_{T_e}f\,\mathrm{d}V\f$
!>
!> different cases are available, depending on \c NDIM and the value of \c ICASE which is set
!> through the runtime option:
!> \verbatim
!> -testcase ICASE
!> \endverbatim
!>
!> 1. \c NDIM=2 and \c ICASE=1 : \f$ \lambda = 2 \mathbf{e}_x + \mathbf{e}_y \quad\quad f = 0\f$.
!> 2. \c NDIM=2 and \c ICASE=2 : \f$ \lambda = \frac{1}{d+1}\left( \sum_{j=1}^{d+1} u_j \right) \mathbf{e}_x + \mathbf{e}_y \quad\quad f = 0\f$.
!> 3. \c NDIM=2 and \c ICASE=3 : \f$ \lambda = 2 \mathbf{e}_x + \mathbf{e}_y \quad\quad f = 0\f$.
!> 4. \c NDIM=2 and \c ICASE=4 : \f$ \lambda = 2 \mathbf{e}_x + \mathbf{e}_y \quad\quad f = 4x -2y\f$.
!> 5. \c NDIM=2 and \c ICASE=5 : \f$ \lambda = 2 \mathbf{e}_x + \mathbf{e}_y \quad\quad f = xye^{\left(x+y\right)}\f$.
!> 6. \c NDIM=2 and \c ICASE=7 : \f$ \lambda = \cos\delta \mathbf{e}_x + \sin\delta \mathbf{e}_y \quad\quad f(x,y) = -\left(x\cos\delta+x\sin\delta\right)\f$.
!> 7. \c NDIM=2 and \c ICASE=8 : \f$ \lambda = \mathbf{e}_x + 2 \mathbf{e}_y \quad\quad f = 0\f$.
!> 8. \c NDIM=2 and \c ICASE=6 : \f$ \lambda = 2\pi y \mathbf{e}_x - 2\pi x \mathbf{e}_y \quad\quad f = 0\f$. This is for the rotating hump testcase, see [De Palma et al. Journal of Computational Physics 208 (2005) 1–3](http://dx.doi.org/doi:10.1016/j.jcp.2004.11.023)
!> 9. \c NDIM=2 and \c ICASE=9 : \f$ \lambda = \mathbf{e}_x \quad\quad f = 0\f$.
!> 10. \c NDIM=3 and \c ICASE=1 : \f$ \lambda = 0.75 \mathbf{e}_x + 0.875 \mathbf{e}_y + 1. \mathbf{e}_z \quad\quad f = 0\f$.
!> 11. \c NDIM=3 and \c ICASE=2 : \f$ \lambda = x \mathbf{e}_x + 0.2 \mathbf{e}_y + y \mathbf{e}_z \quad\quad f = 0\f$.
!>
!> @param[in] IELEM is the current simplicial cell
!> @param[out] X is the cell averaged convection speed \f$\hat{\lambda}\f$ is the current simplicial cell
!> @param[in] XYZ are the \c NDIM Cartesian coordinates of the \c NOFVERT vertices of cell \c IELEM
!> @param[in] U are the \c NOFVERT values of the dependent variable within the \c NOFVERT vertices of cell \c IELEM
!> @param[in] NDIM nof dimensions of the space
!> @param[in] NOFVERT \c =NDIM nof vertices of the current simplex
!> @param[out] FSOU is the integral of the source term, i.e. \f$\int_{T_e}f\,\mathrm{d}V\f$
!>
!> \author $Author: abonfi $
!> \version $Revision: 1.9 $
!> \date $Date: 2013/09/23 11:23:35 $
!>
!>
!>
!>
!>
      SUBROUTINE ADVECT(IELEM,X,XYZ,U,NDIM,NOFVERT,FSOU)
c
c     $Id: advect.f,v 1.9 2013/09/23 11:23:35 abonfi Exp $
c
      IMPLICIT NONE
c
c
      INCLUDE 'constants.h'
      INCLUDE 'scalar.com'
      INCLUDE 'visco.com'
C
C
C
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION FSOU
      INTEGER IELEM,NDIM,NOFVERT
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION U(NOFVERT),X(NDIM),XYZ(NDIM,NOFVERT)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION DELTA,SUMY,SUMX
      INTEGER IVERT,IERR,IOPT
      CHARACTER*30 ERRMSG
C     ..
C     .. External Functions ..
      DOUBLE PRECISION FUNSOU1,FUNSOU2
      EXTERNAL FUNSOU1,FUNSOU2
C     ..
C     .. External Subroutines ..
      EXTERNAL SETERR
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC COS,SIN
C     ..
      DATA ERRMSG(1:25)/'ADVECT - INVALID ICASE = '/
C
      FSOU = ZERO
      DELTA = 22.5d0*PI/180.d0
C
      GOTO (200,300) NDIM - 1
c
c       2D scalar testcases
c
  200 CONTINUE
c
c Linear advection
c
      IF (ICASE.EQ.1) THEN
          X(1) = TWO
          X(2) = ONE

      ELSEIF (ICASE.EQ.2) THEN
          SUMX = 0.00
          DO 3 IVERT = 1,NOFVERT
    3     SUMX = SUMX + U(IVERT)
          X(1) = SUMX/NOFVERT
          X(2) = ONE

      ELSEIF (ICASE.EQ.3) THEN
          X(1) = TWO
          X(2) = ONE

      ELSEIF (ICASE.EQ.4) THEN
          X(1) = TWO
          X(2) = ONE
C
C     source term 4*x - 2*y
C
          SUMX = 0.d0
          SUMY = 0.d0
          DO 5 IVERT = 1,NOFVERT
              SUMX = SUMX + XYZ(1,IVERT)
              SUMY = SUMY + XYZ(2,IVERT)
    5     CONTINUE
          SUMX = 4.d0*SUMX/NOFVERT
          SUMY = -2.d0*SUMY/NOFVERT
          FSOU = SUMX + SUMY

      ELSEIF (ICASE.EQ.5) THEN
          X(1) = TWO
          X(2) = ONE
C
C     source term for U(X,Y) = X*Y*EXP(X+Y)
C
C     first option: evaluate the source term in the
C     centroid of the triangle
C
CXX       SUMX = 0.d0
CXX       SUMY = 0.d0
CXX       DO 7 IVERT = 1,NOFVERT
CXX           SUMX = SUMX + XYZ(1,IVERT)
CXX           SUMY = SUMY + XYZ(2,IVERT)
CXX 7     CONTINUE
CXX       SUMX = SUMX/NOFVERT
CXX       SUMY = SUMY/NOFVERT
CXX       FSOU = FUNSOU1(SUMX,SUMY,X(1),X(2),REINV)
C
C     second option: take the arithmetic mean of the
C     values in the nodes
C
          SUMX=0.d0
          DO 8 IVERT = 1 , NOFVERT
             SUMX = SUMX +
     +       FUNSOU1(XYZ(1,IVERT),XYZ(2,IVERT),X(1),X(2),REINV)
    8     CONTINUE
          FSOU = SUMX/NOFVERT
      ELSEIF (ICASE.EQ.6) THEN
          X(1) = COS(DELTA)
          X(2) = SIN(DELTA)
C
C     source term for U(X,Y) =
C
C     first option: evaluate the source term in the
C     centroid of the triangle
C
          SUMX = 0.d0
          SUMY = 0.d0
          DO 9 IVERT = 1,NOFVERT
              SUMX = SUMX + XYZ(1,IVERT)
              SUMY = SUMY + XYZ(2,IVERT)
    9     CONTINUE
          SUMX = SUMX/NOFVERT
          SUMY = SUMY/NOFVERT
          FSOU = FUNSOU2(SUMX,SUMY,DELTA)
C
C     second option: take the arithmetic mean of the
C     values in the nodes
C
C         SUMX=0.d0
C         DO 8 IVERT = 1 , NOFVERT
C            SUMX = SUMX +
C    +       FUNSOU2(XYZ(1,IVERT),XYZ(2,IVERT),DELTA)
C   8     CONTINUE
C         FSOU = SUMX/NOFVERT
      ELSEIF (ICASE.EQ.7) THEN
          X(1) = ONE
          X(2) = TWO
      ELSEIF (ICASE.EQ.8) THEN
          SUMX = 0.d0
          SUMY = 0.d0
          DO IVERT = 1,NOFVERT
              SUMX = SUMX + XYZ(1,IVERT)
              SUMY = SUMY + XYZ(2,IVERT)
          ENDDO
          SUMX = SUMX/NOFVERT
          SUMY = SUMY/NOFVERT
          X(1) = TWO*PI*SUMY 
          X(2) =-TWO*PI*SUMX
      ELSEIF (ICASE.EQ.9) THEN
          X(1) = ONE
          X(2) = ZERO
      ELSE
          GOTO 666

      ENDIF

      RETURN
c
c       3D scalar testcases
c
  300 CONTINUE
c
c Linear advection
c
      IF (ICASE.EQ.1) THEN
          X(1) = 0.75d0
          X(2) = 0.875d0
          X(3) = ONE
c
c Spiral testcase
c
      ELSEIF (ICASE.EQ.2) THEN
          SUMX = ZERO
          SUMY = ZERO
          DO 1 IVERT = 1,NOFVERT
              SUMX = SUMX + XYZ(3,IVERT)
    1     SUMY = SUMY - XYZ(1,IVERT)
          X(1) = SUMX/NOFVERT
          X(2) = 0.2d0
          X(3) = SUMY/NOFVERT
c
      ELSE
          GOTO 666

      ENDIF

C
      RETURN

  666 CONTINUE
      WRITE(ERRMSG(26:30),FMT=100)ICASE
  100 FORMAT(I5.5)
      IERR=666
      IOPT=1
      CALL SETERR(ERRMSG,30,666,1)
C
      END
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
!> \version $Revision: 1.9 $
!> \date $Date: 2013/09/23 11:23:35 $
!>
!>
      DOUBLE PRECISION FUNCTION FUNSOU1(X,Y,A,B,EPS)
C
C     $Id: advect.f,v 1.9 2013/09/23 11:23:35 abonfi Exp $
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
!> \version $Revision: 1.9 $
!> \date $Date: 2013/09/23 11:23:35 $
!>
!>
      DOUBLE PRECISION FUNCTION FUNSOU2(X,Y,DELTA)
C
C     $Id: advect.f,v 1.9 2013/09/23 11:23:35 abonfi Exp $
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
