C
C     $Id: limiters.f,v 1.1 2001/09/21 07:13:11 abonfi Exp $
C
      DOUBLE PRECISION FUNCTION MINMOD(X,Y)
C
      IMPLICIT NONE
C
C
C
C
C
C
C
C     .. Parameters ..
      DOUBLE PRECISION HALF,ONE
      PARAMETER (HALF=0.5D0,ONE=1.D0)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION X,Y
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,DMIN1,SIGN
C     ..
      MINMOD = HALF* (ONE+SIGN(ONE,X*Y))*HALF*
     +         (SIGN(ONE,X)+SIGN(ONE,Y))*DMIN1(ABS(X),ABS(Y))
      RETURN

      END
C
C
      DOUBLE PRECISION FUNCTION HARMONIC(X,Y)
C
      IMPLICIT NONE
C
C
C
C
C     .. Parameters ..
      DOUBLE PRECISION HALF,ONE,TWO
      PARAMETER (HALF=0.5D0,ONE=1.D0,TWO=2.D0)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION X,Y
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC SIGN
C     ..
      HARMONIC = HALF* (ONE+SIGN(ONE,X*Y))*TWO*X*Y/ (X+Y)
      RETURN

      END
C
C
      DOUBLE PRECISION FUNCTION VANALBADA(X,Y)
C
      IMPLICIT NONE
C
C
C
C
C     .. Parameters ..
      DOUBLE PRECISION HALF,ONE
      PARAMETER (HALF=0.5D0,ONE=1.D0)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION X,Y
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC SIGN
C     ..
      VANALBADA = HALF* (ONE+SIGN(ONE,X*Y))*X*Y* (X+Y)/ (X*X+Y*Y)
      RETURN

      END
C
C
      DOUBLE PRECISION FUNCTION SUPERBEE(X,Y)
C
      IMPLICIT NONE
C
C
C
C
C
C     .. Parameters ..
      DOUBLE PRECISION ZERO,HALF,ONE,TWO
      PARAMETER (ZERO=0.D0,HALF=0.5D0,ONE=1.D0,TWO=2.D0)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION X,Y
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC DMAX1,DMIN1,SIGN
C     ..
      IF (Y.LE.1.D-15) THEN
          SUPERBEE = ZERO

      ELSE
          SUPERBEE = HALF* (ONE+SIGN(ONE,X*Y))*Y*
     +               DMAX1(DMIN1(TWO*X/Y,ONE),DMIN1(X/Y,TWO))
      ENDIF

      RETURN

      END
