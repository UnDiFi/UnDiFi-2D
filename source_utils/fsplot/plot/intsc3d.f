      SUBROUTINE INTSC3D(ICLR,XA,XB,YA,YB,ZA,ZB,Z,NOFVAR,COOR,NDIM,
     +           NPOIN,ICELNOD,ICELCEL,NOFVERT,NELEM,IBNDFAC,NBFAC,
     +           IPNTR,SKINF,NWFAC,FNAME)
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
      REAL*8 AA,BB,CC,DD,REYNO
      INTEGER ICLR,NBFAC,NDIM,NELEM,NOFVAR,NOFVERT,NPOIN,NWFAC
      CHARACTER FNAME* (*)
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION COOR(NDIM,NPOIN),SKINF(NWFAC),Z(NOFVAR,NPOIN)
      INTEGER IBNDFAC(3,NBFAC),ICELCEL(NOFVERT,NELEM),
     +        ICELNOD(NOFVERT,NELEM),IPNTR(2,NWFAC)
C     ..
C     .. Local Scalars ..
      REAL*8 DIST,FUNY,S,T,TAUW,UTAU,X0,X1,X2,XA,XB,XC,XMAX,XMIN,Y0,Y1,
     +       Y2,YA,YB,YC,YMAX,YMIN,U,V,R,H,PRES,DX,DY,ANX,ANY,RR,YY,XX,
     +       TEMP,UX,ZMIN,ZMAX,Z0,Z1,Z2,ZA,ZB,ZC,ZZ
      INTEGER I,IDUMMY,IE,IFACE,INTERSECT,IPNT,IUNIT,IV,IVAR,IXDRS,J,JV,
     +        N1,N2,KUNIT,IFAIL
      LOGICAL COMPRESSIBLE,NORMAL_TO_THE_WALL
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION TMP(5)
C     ..
C     .. External Functions ..
      INTEGER ICYCL,INITXDR
      EXTERNAL ICYCL,INITXDR
C     ..
C     .. External Subroutines ..
      INTEGER IXDRIMAT,IXDRINT
      EXTERNAL IXDRIMAT,IXDRINT
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,MAX,MIN,SQRT
C     ..
C     .. Data statements ..

      DATA TMP/5*0.D0/
C     ..

C
C  calcola l'eq del piano
C
      IF (ABS(BB).GT.ABS(AA)) THEN
          XA = XMIN
          YA = - (AA*XA+CC)/BB
          XB = XMAX
          YB = - (AA*XB+CC)/BB

      ELSE
          YA = YMIN
          XA = - (BB*YA+CC)/AA
          YB = YMAX
          XB = - (BB*YB+CC)/AA
      ENDIF
C
      RR = 2.d0 *MAX(ABS(XMAX-XMIN),ABS(YMAX-YMIN),ABS(ZMAX-ZMIN)) 
C
C     write(6,*)xa,ya,xb,yb
C
      IUNIT = 20
      KUNIT = 50
C
      DO 3 IV = 1,(NOFVERT-1)
          JVERT = ICELNOD(ICYCL(IVERT+IV,NOFVERT),IE)
          
C ... find the intersection of the rake with the boundary
          IF (IBNDFAC(3,I).NE.ICLR) GOTO 3
          IE = IBNDFAC(1,I)
          IV = IBNDFAC(2,I)
C ... nodes on the boundary
          N1 = ICELNOD(ICYCL(IV+1,NOFVERT),IE)
          N2 = ICELNOD(ICYCL(IV+2,NOFVERT),IE)
          N3 = ICELNOD(ICYCL(IV+3,NOFVERT),IE)
          X1 = COOR(1,N1)
          Y1 = COOR(2,N1)
          Z1 = COOR(3,N1)
          X2 = COOR(1,N2)
          Y2 = COOR(2,N2)
          Z2 = COOR(3,N2)
          X3 = COOR(1,N3)
          Y3 = COOR(2,N3)
          Z3 = COOR(3,N3)
C
C     equazione del piano
C
          call piano(X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3,AA,BB,CC,DD)
C
C ... check whether these intersect with the line
          T = PLPT(AA,BB,CC,DD,XA,XB,YA,YB,ZA,ZB)
      t = plpt(aa,bb,cc,dd,x1,x2,y1,y2,z1,z2)
      X0 = x1 + (x2 - x1)*t
      Y0 = y1 + (y2 - y1)*t
      Z0 = z1 + (z2 - z1)*t
C ... if not go to the next boundary edge
          IF (INTERSECT.NE.1) GOTO 3
C
          DX = X2-X1 
          DY = Y2-Y1 
          TEMP = 1.d0/SQRT(DX*DX+DY*DY)
          DX = DX*TEMP
          DY = DY*TEMP
C
C     normal to the wall
C
          ANX = -DY
          ANY =  DX
C
C     compute the intersection
C
          X0 = (1.D0-S)*X1 + S*X2
          Y0 = (1.D0-S)*Y1 + S*Y2
C
C     coordinates of a far away point, on the normal to the surface.
C
          XX = X0 + RR*ANX
          YY = Y0 + RR*ANY
C
C         write(6,*)x0,y0,xx,yy
C         pause
C
          WRITE (6,FMT=*) 'Found an intersection at ',X0,Y0
          WRITE (6,FMT=*) 'Neighbouring meshpoints are ',n1,n2,n3
C
C     compute values at X0,Y0 by linear interpolation
C
          DO 4 IVAR = 1,NOFVAR
              TMP(IVAR) = (1.D0-S)*Z(IVAR,N1) + S*Z(IVAR,N2)
C             write(6,*)Z(IVAR,N1),Z(IVAR,N2),TMP(IVAR)
    4     CONTINUE
C         stop
          IF(COMPRESSIBLE)THEN
                 TMP(4) = TMP(4)/TMP(1)
                 TMP(3) = TMP(3)/TMP(1)
                 TMP(2) = TMP(2)/TMP(1)
                 TMP(1) = TMP(1)*TMP(1)
          ENDIF
          IF(NORMAL_TO_THE_WALL)THEN
              IF(COMPRESSIBLE)THEN
                  TMP(3) = TMP(3)*DX+TMP(4)*DY
                  TMP(4) = TMP(3)*ANX+TMP(4)*ANY
              ELSE
                  TMP(2) = TMP(2)*DX+TMP(3)*DY
                  TMP(3) = TMP(2)*ANX+TMP(3)*ANY
              ENDIF
          ENDIF
C
C    find the friction velocity at the wall
C    (at present only for incompressible flows)
C
          TAUW = 1.D0
          FUNY = 1.D0
          IF (NOFVAR.EQ.1) THEN
              UTAU = 1.d0
              GOTO 15
          ENDIF
          DO 13 IFACE = 1,NWFAC
C    find the wall shear stress for the current edge
C    remind that shear stress are associated with edges
              IF (IE.NE.IPNTR(1,IFACE)) GOTO 13
              IF (IV.NE.IPNTR(2,IFACE)) GOTO 13
              TAUW = ABS(SKINF(IFACE))
              UTAU = SQRT(TAUW)
              FUNY = SQRT(TAUW)*REYNO
              GOTO 15

   13     CONTINUE
          WRITE(6,FMT=*)'Could not find the corresponding shear stress'
          TAUW = 1.d0
          UTAU = 1.d0
          FUNY = 1.d0

   15     CONTINUE
C
          WRITE (6,FMT=180)
  180 FORMAT(5X,'wall distance, y+, flow variables, u+')
          IUNIT = IUNIT + 1
          KUNIT = KUNIT + 1
          WRITE (6,FMT=*) 'Writing flow variables to UNIT ',IUNIT
          WRITE (6,FMT=*) 'Writing b.l. variables to UNIT ',KUNIT
C
C    writes wall distance, y+, variables, u+ for the wall node
C
          IF( COMPRESSIBLE )THEN
              TMP(5) = PRES(TMP(1),TMP(2),TMP(3),TMP(4))
          ELSE
              TMP(5) = 0.D0
          ENDIF 
          WRITE (IUNIT,FMT=80) X0,Y0,0.D0,
     +(TMP(IVAR),IVAR=1,NOFVAR),TMP(5)
          WRITE (KUNIT,FMT=80) X0,Y0,0.D0,0.D0
C
C    loop over the other edges of the current triangle
C    and check which of the two intersects the rake
C
    2     DO 5 J = 1,2
              JV = ICYCL(IV+J,NOFVERT)
              N1 = ICELNOD(ICYCL(JV+1,NOFVERT),IE)
              N2 = ICELNOD(ICYCL(JV+2,NOFVERT),IE)
              X1 = COOR(1,N1)
              Y1 = COOR(2,N1)
              X2 = COOR(1,N2)
              Y2 = COOR(2,N2)
              IF(NORMAL_TO_THE_WALL)THEN
                  INTERSECT = GR_LINE_CROSS(X0,Y0,XX,YY,X1,Y1,X2,Y2,S,T)
              ELSE
                  INTERSECT = GR_LINE_CROSS(XA,YA,XB,YB,X1,Y1,X2,Y2,S,T)
              ENDIF
              IF (INTERSECT.EQ.1) GOTO 7
    5     CONTINUE
          WRITE(6,*) 'Uh! Oh! Smthg. went wrong in subroutine CUT'
          CALL EXIT(1)

    7     CONTINUE
C    compute the intersection point
          XC = (1.D0-S)*X1 + S*X2
          YC = (1.D0-S)*Y1 + S*Y2
C
C     write(6,*)xc,yc,utau
C
C     compute values at XC,YC by linear interpolation
C     variables are:
C     density,total enthalpy,u,v (compressible)
C     pressure,u,v (incompressible)
C
          DO 9 IVAR = 1,NOFVAR
              TMP(IVAR) = (1.D0-S)*Z(IVAR,N1) + S*Z(IVAR,N2)
    9     CONTINUE
          IF(COMPRESSIBLE)THEN
                 TMP(4) = TMP(4)/TMP(1)
                 TMP(3) = TMP(3)/TMP(1)
                 TMP(2) = TMP(2)/TMP(1)
                 TMP(1) = TMP(1)*TMP(1)
C
          ENDIF
          IF(NORMAL_TO_THE_WALL)THEN
              IF(COMPRESSIBLE)THEN
                  TMP(3) = TMP(3)*DX+TMP(4)*DY
                  TMP(4) = TMP(3)*ANX+TMP(4)*ANY
              ELSE
                  TMP(2) = TMP(2)*DX+TMP(3)*DY
                  TMP(3) = TMP(2)*ANX+TMP(3)*ANY
              ENDIF
          ENDIF
          DIST = SQRT((XC-X0)**2+ (YC-Y0)**2)
C
C    writes wall distance, y+, variables, u+
C    >>> incompressible <<<
C
          IF( COMPRESSIBLE )THEN
C    compute pressure
              TMP(5) = PRES(TMP(1),TMP(2),TMP(3),TMP(4))
              UX = TMP(3)
          ELSE
C    compute velocity magnitude
              TMP(5) = SQRT(TMP(2)*TMP(2)+TMP(3)*TMP(3))
              UX = TMP(2)
          ENDIF 
          WRITE (IUNIT,FMT=80) XC,YC,DIST,
     +      (TMP(IVAR),IVAR=1,NOFVAR),TMP(5)
C
          WRITE (KUNIT,FMT=80) XC,YC,DIST*FUNY,UX/UTAU
C
          IE = ICELCEL(JV,IE)
C
C                  +
C                /   \
C               /     \
C              /       \
C             /         \
C            /           \
C          + ------------ +
C            \           /
C             \         /
C              \       /
C               \     /
C                \   /
C                  +
C
C
          IF (IE.EQ.0 .OR. IE.GT.NELEM) GOTO 3
          DO 11 IV = 1,NOFVERT
              IF (ICELNOD(IV,IE).NE.N1 .AND.
     +            ICELNOD(IV,IE).NE.N2) GOTO 2
   11     CONTINUE
          WRITE(6,*) 'Uh! Oh! Smthg. went wrong in subroutine CUT'
          CALL EXIT(1)

    3 CONTINUE
      RETURN

   80 FORMAT (8 (E16.9,1X))
  100 FORMAT (/10("*"),'Probing normal to the wall',/)
  200 FORMAT(5X,'DATA FOR ',A,' FLOWS')

      END
