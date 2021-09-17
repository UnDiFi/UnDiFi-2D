c...
c...
c...
      SUBROUTINE MLUSOL(A,B,C,ORDER)
c...
c...
c...FIRST L(INV)*B
c...
C     .. Scalar Arguments ..
      INTEGER ORDER
C     ..
C     .. Array Arguments ..
      REAL*8 A(ORDER,ORDER),B(ORDER,ORDER),C(ORDER,ORDER)
      REAL*8 SUM(5)
C     ..
C     .. Local Scalars ..
      INTEGER JM,JMJM,JR,JRJR,JRM1,JRP1
C     ..
      TEMP=1.d0/A(1,1)
      DO 10 JR=1,ORDER
         C(1,JR) = B(1,JR)*TEMP
   10 CONTINUE
      DO 14 JR = 2,ORDER
          JRM1 = JR - 1
          DO 13 IR=1,ORDER
          SUM(IR) = B(JR,IR)
   13     CONTINUE
          DO 14 IR = 1,ORDER
          DO 12 JM = 1,JRM1
   12     SUM(IR) = SUM(IR) - A(JR,JM)*C(JM,IR)
C         DO 14 IR = 1,ORDER
   14 C(JR,IR) = SUM(IR)/A(JR,JR)
c...
c...NEXT U(INV) OF L(INV)*B
c...
      DO 18 JRJR = 2,ORDER
          JR = ORDER - JRJR + 1
          JRP1 = JR + 1
           DO 19 IR=1,ORDER
              SUM(IR) = C(JR,IR)
   19   CONTINUE
          DO 16 JMJM = JRP1,ORDER
              JM = ORDER - JMJM + JRP1
              DO 16 IR = 1,ORDER
   16     SUM(IR) = SUM(IR) - A(JR,JM)*C(JM,IR)
              DO 18 IR = 1,ORDER
   18 C(JR,IR) = SUM(IR)
c...
      RETURN

      END
