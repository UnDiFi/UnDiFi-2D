!> \brief \b TBDEST
!> \par Purpose
!>
!> (Spalart & Allmaras model: versions: all, except Bassi)
!>
!> Compute the destruction term
!>
!> \f[
!> s_{D} = -\left[c_w f_w -\frac{c_{b1}}{\kappa^2} f_{t2} \right]\left(\frac{\nu}{d}\right)^2
!> \f]
!>
!>
!> @param[in] TD is the distance \f$ d \f$ from the nearest wall
!> @param[in] TS is the modified strain \f$ \tilde{S} \f$
!> @param[in] TV is turbulent laminar viscosity \f$ \tilde{\nu} \f$
!> @param[in] VI is kinematic laminar viscosity \f$ \nu \f$
!> \author $Author: abonfi $
!> \version $Revision: 1.3 $
!> \date $Date: 2016/11/10 09:48:44 $
!> \warning are there any warnings?
      DOUBLE PRECISION FUNCTION TBDEST(TD,TS,TV,VI)
c2345678
      IMPLICIT NONE
      INCLUDE 'turb.com'
      INCLUDE 'visco.com'

C     .. Scalar Arguments ..
      DOUBLE PRECISION TD,TS,TV,VI
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION DUM1,TCHI,TFT2,TR,DENOM
C     ..
C     .. External Functions ..
      DOUBLE PRECISION TFW
      EXTERNAL TFW
C     ..
C     IF(TS.eq.0.0)TS=1.0d-22
      DENOM=TS* (TK*TD)**2.
      IF(DABS(DENOM).lt.1.d-22) DENOM=1.d-22
Caldo TR = TV/ (TS* (TK*TD)**2.)*REINV
      TR = TV/ DENOM*REINV
      IF (TR.GT.10.) THEN
C         WRITE (*,FMT=*) 'tr>10',TR
C         PAUSE

          TR = 10.
      ENDIF

      DUM1 = 0.
      IF(TTFLAG.EQ.1)THEN
        TCHI=TV/VI
        DUM1=TCB1/TK**2*TFT2(TCHI)
      ENDIF
      TBDEST = - (TCW1*TFW(TR)-DUM1)* (TV/TD)**2. *REINV
      RETURN

      END
