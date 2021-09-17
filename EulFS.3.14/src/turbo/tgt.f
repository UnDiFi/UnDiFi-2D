C
      DOUBLE PRECISION FUNCTION TGT(DU)
*  ==============================================================
*
*
*      FILENAME: tft2.f
*
*
*
*      DESCRIPTION:
*      -----------
*
*      turbulence model function
*      (Spalart & Allmaras model: version)
*
*
*
*  ==============================================================

c2345678
      IMPLICIT NONE 
      INCLUDE 'turb.com'

C     .. Scalar Arguments ..
      DOUBLE PRECISION DU
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC MIN
C     ..
      TGT = MIN(0.1D0,DU/TST/TDXT)
C     TGT=0.1D0
      RETURN

      END
