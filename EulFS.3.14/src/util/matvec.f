C
C
c...
c...
c...
      SUBROUTINE LUDECO(A,ORDER)
c...
c...
C     .. Scalar Arguments ..
      INTEGER ORDER
C     ..
C     .. Array Arguments ..
      REAL*8 A(ORDER,1)
C     ..
C     .. Local Scalars ..
      REAL*8 SUM
      INTEGER JC,JM,JR,JRJC,JRJCM1,JRJCP1
C     ..
      DO 8 JC = 2,ORDER
    8 A(1,JC) = A(1,JC)/A(1,1)
      JRJC = 1
   10 CONTINUE
      JRJC = JRJC + 1
      JRJCM1 = JRJC - 1
      JRJCP1 = JRJC + 1
      DO 14 JR = JRJC,ORDER
          SUM = A(JR,JRJC)
          DO 12 JM = 1,JRJCM1
   12     SUM = SUM - A(JR,JM)*A(JM,JRJC)
   14 A(JR,JRJC) = SUM
      IF (JRJC.EQ.ORDER) RETURN
      DO 18 JC = JRJCP1,ORDER
          SUM = A(JRJC,JC)
          DO 16 JM = 1,JRJCM1
   16     SUM = SUM - A(JRJC,JM)*A(JM,JC)
   18 A(JRJC,JC) = SUM/A(JRJC,JRJC)
      GOTO 10

      END
C
c...
c...SUBROUTINE TO SOLVE LINEAR ALGEBRAIC SYSTEM OF
c...EQUATIONS A*C=B AND STORE RESULTS IN VECTOR C.
c...
c...
c...
c...
      SUBROUTINE LUSOLV(A,B,C,ORDER)
c...
c...
c...FIRST L(INV)*B
c...
C     .. Scalar Arguments ..
      INTEGER ORDER
C     ..
C     .. Array Arguments ..
      REAL*8 A(ORDER,1),B(1),C(1)
C     ..
C     .. Local Scalars ..
      REAL*8 SUM
      INTEGER JM,JMJM,JR,JRJR,JRM1,JRP1
C     ..
      C(1) = B(1)/A(1,1)
      DO 14 JR = 2,ORDER
          JRM1 = JR - 1
          SUM = B(JR)
          DO 12 JM = 1,JRM1
   12     SUM = SUM - A(JR,JM)*C(JM)
   14 C(JR) = SUM/A(JR,JR)
c...
c...NEXT U(INV) OF L(INV)*B
c...
      DO 18 JRJR = 2,ORDER
          JR = ORDER - JRJR + 1
          JRP1 = JR + 1
          SUM = C(JR)
          DO 16 JMJM = JRP1,ORDER
              JM = ORDER - JMJM + JRP1
   16     SUM = SUM - A(JR,JM)*C(JM)
   18 C(JR) = SUM
c...
      RETURN

      END
