      SUBROUTINE XAA23S(X,Y,Z,rbgn,rend,POUT,a,INFO)
C
      IMPLICIT NONE
C
C     This routine finds the values of the
C     dependent variables along the line
C     of equation AA*x+BB*y+CC=0.
C     The line is traced starting from its intersection
C     with the boundary edges colored ICLR
C
C     IPNTR(1:2,1:NWFAC) gives the element and vertex
C                      of the no-slip faces
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION RBGN(*),REND(*),a(*)
      INTEGER InFO
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION X(3),Y(3),Z(3),POUT(3)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION AA,BB,CC,DD
      DOUBLE PRECISION DENOM,S,T,X0,Y0,Z0
      DOUBLE PRECISION X1,Y1,Z1,X2,Y2,Z2
      INTEGER I,J,K,L
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION TMP(5),DR(3,4),QMAT(3,3),DR3(3,4)
      DOUBLE PRECISION XP(4),YP(4),TOLER
      INTEGER IDXS(4)
      PARAMETER(TOLER=1.d-10)
C     ..
C     .. External Functions ..
      INTEGER ICYCL
      DOUBLE PRECISION PLPT,AREA
      EXTERNAL ICYCL,PLPT,AREA
C     ..
C     .. External Subroutines ..
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,MAX,MIN,SQRT
C     ..
C     .. Data statements ..

      DATA TMP/5*0.D0/
C     ..

      x1 = rbgn(1) 
      y1 = rbgn(2) 
      z1 = rbgn(3) 
      x2 = rend(1) 
      y2 = rend(2) 
      z2 = rend(3) 
      pout(1) = 1.d+38
      pout(2) = 1.d+38
      pout(3) = 1.d+38
C
C  calcola l'eq del piano
C
      CALL plane(x(1),y(1),z(1),x(2),y(2),z(2),x(3),y(3),z(3),
     +aa,bb,cc,dd)
C
C ... checks where the plane intersect with the line
C
      t = PLPT(aa,bb,cc,dd,x1,x2,y1,y2,z1,z2,info)
      if(info .NE. 0)then
         return
      endif
C
      X0 = x1 + (x2 - x1)*t
      Y0 = y1 + (y2 - y1)*t
      Z0 = z1 + (z2 - z1)*t
!     write(6,*)'a,b,c,d, ',aa,bb,cc,dd
!     write(6,*)x1,y1,z1
!     write(6,*)x2,y2,z2
!     write(6,*)t,x0,y0,z0
!     write(6,*)aa*x0+bb*y0+cc*z0+dd
!     write(6,*)aa*x(1)+bb*y(1)+cc*z(1)+dd
!     write(6,*)aa*x(2)+bb*y(2)+cc*z(2)+dd
!     write(6,*)aa*x(3)+bb*y(3)+cc*z(3)+dd
C
C ... checks if the intersection belongs to the
C     triangular face
C
C     we transform the coords of the vertices of the
C     face into a reference frame (xi,eta,zeta)
C     where zeta is parallel to the normal to
C     the plane and xi is aligned with the
C     shortest edge
C     we form the rotation matrix Q = (xi,eta,zeta)
C     then we tranform any vector a as follows
C
C     a_(xi,eta,zeta) = Q^t a_(x,y,z)
C
C     find the edge of min length
C
      t = 1.d38
      l = 0
      do i = 1,3
         j = icycl(i+1,3)
         k = icycl(i+2,3)
         dr(1,i) = (x(k)-x(j))
         dr(2,i) = (y(k)-y(j))
         dr(3,i) = (z(k)-z(j))
         s = sqrt( dr(1,i)*dr(1,i) + dr(2,i)*dr(2,i) + dr(3,i)*dr(3,i) )
         if( s .LT. t )then
             t = s
             l = i
         endif
      enddo
c
c     the smallest edge is the one facing vertex l
c
      s = sqrt( dr(1,l)*dr(1,l) + dr(2,l)*dr(2,l) + dr(3,l)*dr(3,l) )
      s = 1.d0/s
c
c     transformation matrix
c
      qmat(1,1) = dr(1,l) * s
      qmat(2,1) = dr(2,l) * s
      qmat(3,1) = dr(3,l) * s
c
c
c
      t = sqrt(aa*aa+bb*bb+cc*cc)
      t = 1.d0/t
      qmat(1,3) = aa*t
      qmat(2,3) = bb*t
      qmat(3,3) = cc*t
      call cross_prod( qmat(1,3), qmat(1,1), qmat(1,2) )
c
c     compute P(i)-P(l) in the orig ref frame
c
      do j = 1,3
         dr3(1,j) = x(j) - x(l)
         dr3(2,j) = y(j) - y(l)
         dr3(3,j) = z(j) - z(l)
      enddo
c
c     compute P(0)-P(l) in the orig ref frame
c
      dr3(1,4) = x0 - x(l)
      dr3(2,4) = y0 - y(l)
      dr3(3,4) = z0 - z(l)
!     do i = 1,3
!     write(6,*)l,(dr3(i,j),j=1,4)
!     enddo
c
c     transform into local coordinates
c
      call dgemm('Trans','No',3,4,3,1.d0,qmat,3,dr3,3,0.d0,dr,3)
!     do i = 1,3
!     write(6,*)(dr3(i,j),j=1,4)
!     enddo
c     do i = 1,3
c     write(6,*)(dr(i,j),j=1,3)
c     enddo
c     pause
      do i = 1,4
         xp(i) = dr(1,i)
         yp(i) = dr(2,i)
      enddo
!     write(6,*)'coords ',(xp(j),j=1,4)
!     write(6,*)'coords ',(yp(j),j=1,4)
c
c     the third component of dr2 should be zero
c
      idxs(1) = 1
      idxs(2) = 2
      idxs(3) = 3
      denom = 1.d0/area(xp,yp,3,idxs)
      idxs(3) = 4
      do i = 1,3
         idxs(1) = icycl(1+i,3)
         idxs(2) = icycl(2+i,3)
         a(i) = area(xp,yp,3,idxs)*denom
         if(ABS(a(i)).LE.TOLER)a(i) = 0.d0
      enddo
      s = min( a(1), a(2), a(3) ) 
      t = max( a(1), a(2), a(3) ) 
      if( ( s .GE. 0.d0 .AND. s .LE. 1.d0 ) .AND.
     +    ( t .GE. 0.d0 .AND. t .LE. 1.d0 ) )then
          info = 0
          pout(1) = x0
          pout(2) = y0
          pout(3) = z0
          write(16,*)'info = ',info,' Area coords: ',(a(j),j=1,3)
!         write(6,*)'info = ',info,' Area coords: ',(a(j),j=1,3)
      else
          info = 1
          write(16,*)'info = ',info,' Area coords: ',(a(j),j=1,3)
          write(6,*)'info = ',info,' Area coords: ',(a(j),j=1,3)
      endif
!     write(16,*)'Area coords ',(a(j),j=1,3),1.d0/denom
!     pause
      RETURN

!  80 FORMAT (8 (E16.9,1X))

      END
