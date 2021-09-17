      SUBROUTINE PLANE(X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3,A,B,C,D)
C
C     calcola i coefficienti dell'eqn. del piano
C
      IMPLICIT NONE
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION A,B,C,D,X1,X2,X3,Y1,Y2,Y3,Z1,Z2,Z3
      DOUBLE PRECISION X,Y,Z
C     ..
      A = y1*(z2 - z3) + y2*(z3 - z1) + y3*(z1 - z2)
      B = z1*(x2 - x3) + z2*(x3 - x1) + z3*(x1 - x2)
      C = x1*(y2 - y3) + x2*(y3 - y1) + x3*(y1 - y2)
      D = x1*(y2*z3 - y3*z2) + x2*(y3*z1 - y1*z3) + 
     &    x3*(y1*z2 - y2*z1)
      D = -D
C
      RETURN

      END
