      DOUBLE PRECISION FUNCTION YDIST(ICELNOD,NOFVERT,CORG,NDIM,
     &IELEM,IVERT)
      IMPLICIT NONE
      INTEGER NOFVERT,NDIM,IELEM,IVERT
      INTEGER  ICELNOD(NOFVERT,*)
      DOUBLE PRECISION CORG(NDIM,*)
      DOUBLE PRECISION XA,XB,XC,YA,YB,YC,ZA,ZB,ZC,
     &XI,YI,ZI,X0,Y0,Z0,  HELP
      DOUBLE PRECISION AA,BB,CC,dd,DA,DB,DC,DP,R
      INTEGER IA,IB,IC
      INTEGER ICYCL
c
      IA = ICELNOD(IVERT,IELEM)

      X0 = CORG(1,IA)
      Y0 = CORG(2,IA)
      Z0 = CORG(3,IA)

      IA = ICELNOD(ICYCL(IVERT+1,NOFVERT),IELEM)
      IB = ICELNOD(ICYCL(IVERT+2,NOFVERT),IELEM)
      IC = ICELNOD(ICYCL(IVERT+3,NOFVERT),IELEM)

      XA = CORG(1,IA)
      XB = CORG(1,IB)
      XC = CORG(1,IC)
      YA = CORG(2,IA)
      YB = CORG(2,IB)
      YC = CORG(2,IC)
      ZA = CORG(3,IA)
      ZB = CORG(3,IB)
      ZC = CORG(3,IC)
C
C Find the eqn. of the plane
C
      CALL PLANE(XA,YA,ZA,XB,YB,ZB,XC,YC,ZC,AA,BB,CC,DD)
C
C     A point belonging to the line passing through
C     (X0,Y0,Z0) and perpendicular to the plane is
C     X = X0 + R*AA
C     Y = Y0 + R*BB
C     Z = Z0 + R*CC
C     the plane is described by eqn.:
C     AA*X+BB*Y+CC*Z+DD = 0
C
C Subject 5.05: How do I find the intersection of a line and a plane?
C
C    If the plane is defined as:
C
C        a*x + b*y + c*z + d = 0
C
C    and the line is defined as:
C
C        x = x1 + (x2 - x1)*t = x1 + i*t
C        y = y1 + (y2 - y1)*t = y1 + j*t
C        z = z1 + (z2 - z1)*t = z1 + k*t
C
C    Then just substitute these into the plane equation. You end up
C    with:
C
C        t = - (a*x1 + b*y1 + c*z1 + d)/(a*i + b*j + c*k)
C
C    When the denominator is zero, the line is contained in the plane
C    if the numerator is also zero (the point at t=0 satisfies the
C    plane equation), otherwise the line is parallel to the plane.
C
C
      R = - (AA*X0+BB*Y0+CC*Z0+DD)/ (AA*AA+BB*BB+CC*CC)

C   Let I be the point of perpendicular projection of C onto the plane
C
      XI = X0 + R*AA
      YI = Y0 + R*BB
      ZI = Z0 + R*CC
C
      HELP = AA*XI+BB*YI+CC*ZI+DD
      IF(ABS(HELP).GE.1.e-7)THEN
         write(6,*)' check is = ',help
      ENDIF
C
caldo         write(6,*)xa,ya,za
caldo         write(6,*)xb,yb,zb
caldo         write(6,*)xc,yc,zc
caldo         write(6,*)xi,yi,zi
C
C    if the point (XI,YI,ZI) belongs to the face defined
C    by vertices A,B,C that's the shortest distance
C
C
      DA = (X0-XA)**2 + (Y0-YA)**2 + (Z0-ZA)**2
      DB = (X0-XB)**2 + (Y0-YB)**2 + (Z0-ZB)**2
      DC = (X0-XC)**2 + (Y0-YC)**2 + (Z0-ZC)**2
      DP = (X0-XI)**2 + (Y0-YI)**2 + (Z0-ZI)**2

      DA=SQRT(DA)
      DB=SQRT(DB)
      DC=SQRT(DC)
      DP=SQRT(DP)
C
C     distanza dal baricentro
C
caldo XI = (XA+XB+XC)/3.d0
caldo YI = (YA+YB+YC)/3.d0
caldo ZI = (ZA+ZB+ZC)/3.d0
caldo DP = (X0-XI)**2 + (Y0-YI)**2 + (Z0-ZI)**2
caldo DP=SQRT(DP)
C
      YDIST = DP
      RETURN
      END
C
