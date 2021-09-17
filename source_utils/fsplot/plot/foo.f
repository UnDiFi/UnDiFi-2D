      program foo
      real*8 XA,YA,ZA,XB,YB,ZB,XC,YC,ZC,AA,BB,CC,DD
      PIANO(X,Y,Z) = AA*X+BB*Y+CC*Z+DD
      read(5,*)XA,YA,ZA,XB,YB,ZB,XC,YC,ZC
      write(6,*)xa,ya,za,xb,yb,zb,xc,yc,zc
      call PLANE(XA,YA,ZA,XB,YB,ZB,XC,YC,ZC,AA,BB,CC,DD)
      write(6,*)aa,bb,cc,dd
      write(6,*)piano(xa,ya,za) 
      write(6,*)piano(xb,yb,zb) 
      write(6,*)piano(xc,yc,zc) 
      end

      SUBROUTINE PLANE(X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3,A,B,C,D)
C
C     calcola i coefficienti dell'eqn. del piano
C
      IMPLICIT NONE
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION A,B,C,D,X1,X2,X3,Y1,Y2,Y3,Z1,Z2,Z3
      DOUBLE PRECISION X,Y,Z
      DOUBLE PRECISION PIANO
C
      PIANO(X,Y,Z) = A*X+B*Y+C*Z+D
C     ..
      A = y1*(z2 - z3) + y2*(z3 - z1) + y3*(z1 - z2)
      B = z1*(x2 - x3) + z2*(x3 - x1) + z3*(x1 - x2)
      C = x1*(y2 - y3) + x2*(y3 - y1) + x3*(y1 - y2)
      D = x1*(y2*z3 - y3*z2) + x2*(y3*z1 - y1*z3) + 
     &    x3*(y1*z2 - y2*z1)
      D = -D
C
C
      return
c
c     there's a bug, maybe !?!?!?!?!?!
c
c
c     if(
c    +abs(piano(xa,ya,za)) .gt. 1.d-8 .or.
c    +abs(piano(xb,yb,zb)) .gt. 1.d-8 .or.
c    +abs(piano(xc,yc,zc)) .gt. 1.d-8 )then
c     write(6,*)xa,ya,za,xb,yb,zb,xc,yc,zc
c     write(6,*)piano(xa,ya,za) 
c     write(6,*)piano(xb,yb,zb) 
c     write(6,*)piano(xc,yc,zc) 
c     pause
c     endif

      RETURN

      END
