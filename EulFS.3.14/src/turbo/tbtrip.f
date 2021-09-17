C///////////////////////////////////////////////////////////////////////

C       Identification Keywords:

C       @(#)tbtrip.f      1.0.0.1   95/08/22        12:00:00
C       /alpha1/users/mubadmin/scalar/sccs/s.tbtrip.f
C

C///////////////////////////////////////////////////////////////////////

C//////////////////////////////////////////////////////////////////////
C
      DOUBLE PRECISION FUNCTION TBTRIP(TD,TTD,UX,UY,UZ)
*  ==============================================================
*
*
*      FILENAME: tbtrip.f
*
*
*
*      DESCRIPTION:
*      -----------
*
*      Trip term term computation
*      (Spalart & Allmort?? model: easy version)
*
*      FILE HISTORY:
*      ------------
*
*      LANGUAGE : Fortran 77
*      LAST CHANGE : August 95
*
*  ==============================================================

c2345678
c termine di trip
      IMPLICIT NONE
      INCLUDE 'visco.com'

C     .. Scalar Arguments ..
      DOUBLE PRECISION TD,TTD,UX,UY,UZ
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION DU
C     ..
C     .. External Functions ..
      DOUBLE PRECISION TDU,TFT1
      EXTERNAL TDU,TFT1
C     ..
      DU = TDU(UX,UY,UZ)
      TBTRIP = TFT1(TD,TTD,DU)*DU*DU*RE
      RETURN

      END
