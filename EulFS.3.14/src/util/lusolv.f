c...
c...
      SUBROUTINE xxSOLV(A,B,C,ORDER)
c...
c...
C     .. Scalar Arguments ..
      INTEGER ORDER
C     ..
C     .. Array Arguments ..
      REAL*8 A(ORDER,*),B(*),C(*)
C     ..
C     .. Local Scalars ..
      REAL*8 SUM
      GOTO(200,300,400,500)ORDER-1
  200 CONTINUE
       C(1) = B(1)/A(1,1)
       SUM = B(2)
       SUM = SUM - A(2,1)*C(1)
       C(2) = SUM/A(2,2)
       SUM = C(1)
       SUM = SUM - A(1,2)*C(2)
       C(1) = SUM
       RETURN
  300 CONTINUE
       C(1) = B(1)/A(1,1)
       SUM = B(2)
       SUM = SUM - A(2,1)*C(1)
       C(2) = SUM/A(2,2)
       SUM = B(3)
       SUM = SUM - A(3,1)*C(1)
       SUM = SUM - A(3,2)*C(2)
       C(3) = SUM/A(3,3)
       SUM = C(2)
       SUM = SUM - A(2,3)*C(3)
       C(2) = SUM
       SUM = C(1)
       SUM = SUM - A(1,3)*C(3)
       SUM = SUM - A(1,2)*C(2)
       C(1) = SUM
       RETURN
  400 CONTINUE
       C(1) = B(1)/A(1,1)
       SUM = B(2)
       SUM = SUM - A(2,1)*C(1)
       C(2) = SUM/A(2,2)
       SUM = B(3)
       SUM = SUM - A(3,1)*C(1)
       SUM = SUM - A(3,2)*C(2)
       C(3) = SUM/A(3,3)
       SUM = B(4)
       SUM = SUM - A(4,1)*C(1)
       SUM = SUM - A(4,2)*C(2)
       SUM = SUM - A(4,3)*C(3)
       C(4) = SUM/A(4,4)
       SUM = C(3)
       SUM = SUM - A(3,4)*C(4)
       C(3) = SUM
       SUM = C(2)
       SUM = SUM - A(2,4)*C(4)
       SUM = SUM - A(2,3)*C(3)
       C(2) = SUM
       SUM = C(1)
       SUM = SUM - A(1,4)*C(4)
       SUM = SUM - A(1,3)*C(3)
       SUM = SUM - A(1,2)*C(2)
       C(1) = SUM
       RETURN
  500 CONTINUE
       C(1) = B(1)/A(1,1)
       SUM = B(2)
       SUM = SUM - A(2,1)*C(1)
       C(2) = SUM/A(2,2)
       SUM = B(3)
       SUM = SUM - A(3,1)*C(1)
       SUM = SUM - A(3,2)*C(2)
       C(3) = SUM/A(3,3)
       SUM = B(4)
       SUM = SUM - A(4,1)*C(1)
       SUM = SUM - A(4,2)*C(2)
       SUM = SUM - A(4,3)*C(3)
       C(4) = SUM/A(4,4)
       SUM = B(5)
       SUM = SUM - A(5,1)*C(1)
       SUM = SUM - A(5,2)*C(2)
       SUM = SUM - A(5,3)*C(3)
       SUM = SUM - A(5,4)*C(4)
       C(5) = SUM/A(5,5)
       SUM = C(4)
       SUM = SUM - A(4,5)*C(5)
       C(4) = SUM
       SUM = C(3)
       SUM = SUM - A(3,5)*C(5)
       SUM = SUM - A(3,4)*C(4)
       C(3) = SUM
       SUM = C(2)
       SUM = SUM - A(2,5)*C(5)
       SUM = SUM - A(2,4)*C(4)
       SUM = SUM - A(2,3)*C(3)
       C(2) = SUM
       SUM = C(1)
       SUM = SUM - A(1,5)*C(5)
       SUM = SUM - A(1,4)*C(4)
       SUM = SUM - A(1,3)*C(3)
       SUM = SUM - A(1,2)*C(2)
       C(1) = SUM
       RETURN
       END
