c...
c...
c...
      SUBROUTINE xxDECO(A,ORDER)
c...
c...
C     .. Scalar Arguments ..
      INTEGER ORDER
C     ..
C     .. Array Arguments ..
      REAL*8 A(ORDER,*)
C     ..
C     .. Local Scalars ..
      REAL*8 SUM
      GOTO(200,300,400,500)ORDER-1
  200 CONTINUE
       A(1,2) = A(1,2)/A(1,1)
       SUM = A(2,2)
       SUM = SUM - A(2,1)*A(1,2)
       A(2,2) = SUM
       RETURN
  300 CONTINUE
       A(1,2) = A(1,2)/A(1,1)
       A(1,3) = A(1,3)/A(1,1)
       SUM = A(2,2)
       SUM = SUM - A(2,1)*A(1,2)
       A(2,2) = SUM
       SUM = A(3,2)
       SUM = SUM - A(3,1)*A(1,2)
       A(3,2) = SUM
       SUM = A(2,3)
       SUM = SUM - A(2,1)*A(1,3)
       A(2,3) = SUM/A(2,2)
       SUM = A(3,3)
       SUM = SUM - A(3,1)*A(1,3)
       SUM = SUM - A(3,2)*A(2,3)
       A(3,3) = SUM
       RETURN
  400 CONTINUE
       A(1,2) = A(1,2)/A(1,1)
       A(1,3) = A(1,3)/A(1,1)
       A(1,4) = A(1,4)/A(1,1)
       SUM = A(2,2)
       SUM = SUM - A(2,1)*A(1,2)
       A(2,2) = SUM
       SUM = A(3,2)
       SUM = SUM - A(3,1)*A(1,2)
       A(3,2) = SUM
       SUM = A(4,2)
       SUM = SUM - A(4,1)*A(1,2)
       A(4,2) = SUM
       SUM = A(2,3)
       SUM = SUM - A(2,1)*A(1,3)
       A(2,3) = SUM/A(2,2)
       SUM = A(2,4)
       SUM = SUM - A(2,1)*A(1,4)
       A(2,4) = SUM/A(2,2)
       SUM = A(3,3)
       SUM = SUM - A(3,1)*A(1,3)
       SUM = SUM - A(3,2)*A(2,3)
       A(3,3) = SUM
       SUM = A(4,3)
       SUM = SUM - A(4,1)*A(1,3)
       SUM = SUM - A(4,2)*A(2,3)
       A(4,3) = SUM
       SUM = A(3,4)
       SUM = SUM - A(3,1)*A(1,4)
       SUM = SUM - A(3,2)*A(2,4)
       A(3,4) = SUM/A(3,3)
       SUM = A(4,4)
       SUM = SUM - A(4,1)*A(1,4)
       SUM = SUM - A(4,2)*A(2,4)
       SUM = SUM - A(4,3)*A(3,4)
       A(4,4) = SUM
       RETURN
  500 CONTINUE
       A(1,2) = A(1,2)/A(1,1)
       A(1,3) = A(1,3)/A(1,1)
       A(1,4) = A(1,4)/A(1,1)
       A(1,5) = A(1,5)/A(1,1)
       SUM = A(2,2)
       SUM = SUM - A(2,1)*A(1,2)
       A(2,2) = SUM
       SUM = A(3,2)
       SUM = SUM - A(3,1)*A(1,2)
       A(3,2) = SUM
       SUM = A(4,2)
       SUM = SUM - A(4,1)*A(1,2)
       A(4,2) = SUM
       SUM = A(5,2)
       SUM = SUM - A(5,1)*A(1,2)
       A(5,2) = SUM
       SUM = A(2,3)
       SUM = SUM - A(2,1)*A(1,3)
       A(2,3) = SUM/A(2,2)
       SUM = A(2,4)
       SUM = SUM - A(2,1)*A(1,4)
       A(2,4) = SUM/A(2,2)
       SUM = A(2,5)
       SUM = SUM - A(2,1)*A(1,5)
       A(2,5) = SUM/A(2,2)
       SUM = A(3,3)
       SUM = SUM - A(3,1)*A(1,3)
       SUM = SUM - A(3,2)*A(2,3)
       A(3,3) = SUM
       SUM = A(4,3)
       SUM = SUM - A(4,1)*A(1,3)
       SUM = SUM - A(4,2)*A(2,3)
       A(4,3) = SUM
       SUM = A(5,3)
       SUM = SUM - A(5,1)*A(1,3)
       SUM = SUM - A(5,2)*A(2,3)
       A(5,3) = SUM
       SUM = A(3,4)
       SUM = SUM - A(3,1)*A(1,4)
       SUM = SUM - A(3,2)*A(2,4)
       A(3,4) = SUM/A(3,3)
       SUM = A(3,5)
       SUM = SUM - A(3,1)*A(1,5)
       SUM = SUM - A(3,2)*A(2,5)
       A(3,5) = SUM/A(3,3)
       SUM = A(4,4)
       SUM = SUM - A(4,1)*A(1,4)
       SUM = SUM - A(4,2)*A(2,4)
       SUM = SUM - A(4,3)*A(3,4)
       A(4,4) = SUM
       SUM = A(5,4)
       SUM = SUM - A(5,1)*A(1,4)
       SUM = SUM - A(5,2)*A(2,4)
       SUM = SUM - A(5,3)*A(3,4)
       A(5,4) = SUM
       SUM = A(4,5)
       SUM = SUM - A(4,1)*A(1,5)
       SUM = SUM - A(4,2)*A(2,5)
       SUM = SUM - A(4,3)*A(3,5)
       A(4,5) = SUM/A(4,4)
       SUM = A(5,5)
       SUM = SUM - A(5,1)*A(1,5)
       SUM = SUM - A(5,2)*A(2,5)
       SUM = SUM - A(5,3)*A(3,5)
       SUM = SUM - A(5,4)*A(4,5)
       A(5,5) = SUM
       RETURN
       END
