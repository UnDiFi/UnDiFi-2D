!> \brief \b TSC
!> \par Purpose
!>
!> (Spalart & Allmaras model: version: SA-fv3)
!> vorticity module computation
!>
!> The equations are the same as for the "standard" version (SA), with the following exceptions: 
!>
!> \f[
!> \tilde{S} = f_{v_3} \, \Omega + \frac{1}{Re} \frac{\tilde{\nu}}{\kappa^2d^2} f_{v_2}
!> \f]
!> \f[
!> f_{v_2} = \frac{1}{\left(1+\chi/c_{v_2}\right)^3}
!> \f]
!> \f[
!> f_{v_3} = \frac{\left(1+\chi \,f_{v_1}\right)\left(1-f_{v_2}\right)}{\chi}
!> \f]
!> \f[
!> c_{v_3} = 5
!> \f]
!>
!>
!> @param[in] OMEGA is \f$ \Omega \f$
!> @param[in] TD is the distance from the nearest wall
!> @param[in] TVI is the turbulent working variable \f$ \tilde{\nu} \f$
!> @param[in] VI is laminar kinematic viscosity
!> \author $Author: abonfi $
!> \version $Revision: 1.11 $
!> \date $Date: 2020/02/05 15:05:33 $
      DOUBLE PRECISION FUNCTION TSC(OMEGA,TD,TVI,VI)
c2345678
C
      IMPLICIT NONE
C
      INCLUDE 'paramt.h'
      INCLUDE 'bnd.h'
      INCLUDE 'turb.com'
      INCLUDE 'stream.com'
      INCLUDE 'visco.com'
C

C
C     .. Scalar Arguments ..
      DOUBLE PRECISION OMEGA,TD,TVI,VI
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION TCHI,FV2,FV3
C     ..
C     .. External Functions ..
      DOUBLE PRECISION TFV2,TFV3
      EXTERNAL TFV2,TFV3
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC SQRT
C     ..
C
      TCHI = TVI/VI
      FV2 = TFV2(TCHI)
      FV3 = TFV3(TCHI,FV2)
      IF (DABS(TD).LT.1.D-10) TD=1.D-10
      TSC = OMEGA*FV3 + TVI/ (TK*TD)**2.*FV2*REINV

      RETURN

      END
!> \brief \b TSCSA
!> \par Purpose
!>
!> vorticity module computation
!> (Spalart & Allmaras model: version: SA)
!>
!> The equations are those of the "standard" version (SA):
!>
!> \f[
!> \tilde{S} = \Omega + \frac{1}{Re} \frac{\tilde{\nu}}{\kappa^2d^2} f_{v_2}
!> \f]
!> \f[
!> f_{v_1} = \frac{\chi^3}{\left(\chi^3+c_{v_1}^3\right)}
!> \f]
!> \f[
!> f_{v_2} = 1-\frac{\chi}{\left(1+\chi\,f_{v_1}\right)}
!> \f]
!>
!>
!> @param[in] OMEGA is \f$ \Omega \f$
!> @param[in] TD is the distance from the nearest wall
!> @param[in] TVI is the turbulent working variable \f$ \tilde{\nu} \f$
!> @param[in] VI is laminar kinematic viscosity
!> \author $Author: abonfi $
!> \version $Revision: 1.11 $
!> \date $Date: 2020/02/05 15:05:33 $
      DOUBLE PRECISION FUNCTION TSCSA(OMEGA,TD,TVI,VI)
c2345678
C
      IMPLICIT NONE
C
      INCLUDE 'paramt.h'
      INCLUDE 'bnd.h'
      INCLUDE 'turb.com'
      INCLUDE 'stream.com'
      INCLUDE 'visco.com'
C
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION OMEGA,TD,TVI,VI
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION TCHI,FV2
C     ..
C     .. External Functions ..
      DOUBLE PRECISION TFV2SA
      EXTERNAL TFV2SA
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC DABS
C     ..
C
      TCHI = TVI/VI
      FV2 = TFV2SA(TCHI)
      IF (DABS(TD).LT.1.D-10) TD=1.D-10
      TSCSA = OMEGA + TVI/ (TK*TD)**2*FV2*REINV
      RETURN
      END
!> \brief \b TSCSAc
!> \par Purpose
!>
!> vorticity module computation with limiting (c)
!> (Spalart & Allmaras model: version: SA)
!>
!> The equations are those of the "standard" version (SA):
!>
!> \f[
!> \tilde{S} = \Omega + \frac{1}{Re} \frac{\tilde{\nu}}{\kappa^2d^2} f_{v_2}
!> \f]
!> \f[
!> f_{v_1} = \frac{\chi^3}{\left(\chi^3+c_{v_1}^3\right)}
!> \f]
!> \f[
!> f_{v_2} = 1-\frac{\chi}{\left(1+\chi\,f_{v_1}\right)}
!> \f]
!>
!>
!> @param[in] OMEGA is \f$ \Omega \f$
!> @param[in] TD is the distance from the nearest wall
!> @param[in] TVI is the turbulent working variable \f$ \tilde{\nu} \f$
!> @param[in] VI is laminar kinematic viscosity
!> \author $Author: abonfi $
!> \version $Revision: 1.11 $
!> \date $Date: 2020/02/05 15:05:33 $
      DOUBLE PRECISION FUNCTION TSCSAc(OMEGA,TD,TVI,VI)
c2345678
C
      IMPLICIT NONE
C
      INCLUDE 'paramt.h'
      INCLUDE 'constants.h'
      INCLUDE 'bnd.h'
      INCLUDE 'turb.com'
      INCLUDE 'stream.com'
      INCLUDE 'visco.com'
C
C     ..
C     .. Parameters ..
      DOUBLE PRECISION C2,C3
      PARAMETER(C2=0.7d0,C3=0.9d0)
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION OMEGA,TD,TVI,VI
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION TCHI,FV2,HATS
C     ..
C     .. External Functions ..
      DOUBLE PRECISION TFV2SA
      EXTERNAL TFV2SA
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC DABS
C     ..
C
      TCHI = TVI/VI
      FV2 = TFV2SA(TCHI)
      IF (DABS(TD).LT.1.D-10) TD=1.D-10
      HATS = TVI/ (TK*TD)**2*FV2*REINV
      IF(HATS/OMEGA.LT.-C2)THEN
         TSCSAc = OMEGA *
     &  (ONE + ((C2*C2*OMEGA+C3*HATS)/((C3-TWO*C2)*OMEGA-HATS)))
      ELSE 
         TSCSAc = OMEGA + HATS
      ENDIF 
      RETURN
      END
