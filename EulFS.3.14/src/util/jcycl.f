      INTEGER FUNCTION JCYCL (I)
C
      IMPLICIT NONE
C
C     //////////////////////////////////////////////////////////////
C
C     This function brings the value of I back into the interval
C     [1,3] in a cyclic way.
C
C     .. Scalar Arguments ..
C
      INTEGER I
C
C     .. Local Scalars ..
C
      INTEGER IM
C
C     .. Intrinsic Functions ..
C
      INTRINSIC MOD,ISIGN
C
C     .. Executable Statements ..
C
      IM    = MOD(I,3)
      JCYCL = IM + 3*((1-ISIGN(1,(IM-1)))/2)
C
      RETURN
      END
