!> \brief \b TFV3
!> \par Purpose
!>
!> (Spalart & Allmaras model: version: SA-fv3)
!>
!> Compute the function
!>
!> \f[
!> f_{v_3} = \frac{\left(1+\chi \,f_{v_1}\right)\left(1-f_{v_2}\right)}{\chi}
!> \f]
!>
!> In the denominator we take \f$ \max\left(\chi,0.001\right) \f$ following a suggestions by G.A. Ashford
!>
!> @param[in] TCHI is \f$ \frac{\tilde{\nu}}{\nu} \f$
!> @param[in] FV2 is \f$ f_{v_2}\left(\chi\right) \f$
!> \author $Author: abonfi $
!> \version $Revision: 1.3 $
!> \date $Date: 2015/05/11 07:23:09 $
!> \warning Is is maybe not such a good idea to pass FV2: it saves a calculation, but maybe prone to errors
      DOUBLE PRECISION FUNCTION TFV3(TCHI,FV2)

      implicit none
      include 'constants.h'

C     .. Parameters ..
      DOUBLE PRECISION EPS
      PARAMETER (EPS=0.001d0)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION FV2,TCHI
C     ..
C     .. External Functions ..
      DOUBLE PRECISION TFV1
      EXTERNAL TFV1
C     ..
      TFV3 = (ONE+TCHI*TFV1(TCHI))* (ONE-FV2)/MAX(TCHI,EPS)
      RETURN
      END
