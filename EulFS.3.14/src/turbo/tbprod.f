!> \brief \b TBPROD
!> \par Purpose
!>
!> (Spalart & Allmaras model: versions: all, except Bassi)
!>
!> Compute the production term
!>
!> \f[
!> s_{P} = c_{b1}\left(1-f_{t2}\right) \tilde{S} \tilde{\nu}
!> \f]
!>
!>
!> @param[in] TV is turbulent laminar viscosity \f$ \tilde{\nu} \f$
!> @param[in] VI is kinematic laminar viscosity \f$ \nu \f$
!> @param[in] TS is the modified strain \f$ \tilde{S} \f$
!> \author $Author: abonfi $
!> \version $Revision: 1.2 $
!> \date $Date: 2015/05/22 08:11:10 $
!> \warning are there any warnings?
      DOUBLE PRECISION FUNCTION TBPROD(TV,VI,TS)
c2345678
c versione semplificata del termine di produzione
      IMPLICIT NONE
      INCLUDE 'constants.h'
      INCLUDE 'turb.com'

C     .. Scalar Arguments ..
      DOUBLE PRECISION TS,TV,VI
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION DUM1,TCHI,TFT2
C     ..
      DUM1 = ONE
      IF(TTFLAG.EQ.1)THEN
         TCHI = TV/VI
         DUM1=(ONE-TFT2(TCHI))
      ENDIF 
      TBPROD = TCB1*TS*TV*DUM1
      RETURN

      END
