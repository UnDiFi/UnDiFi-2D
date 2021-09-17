      SUBROUTINE BL(ICLR,AA,BB,CC,Z,NOFVAR,COOR,NDIM,NPOIN,ICELNOD,
     +               ICELCEL,NOFVERT,NELEM,IBNDFAC,NBFAC,IPNTR,SKINF,
     +               NWFAC,REYNO,M_INFTY,NORMAL_TO_THE_WALL,FILE)
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
      REAL*8 AA,BB,CC,REYNO,M_INFTY,ETA,ETABAR,YBAR,DELX,DELY,
     &XSAVE,YSAVE,RA,RB
      INTEGER ICLR,NBFAC,NDIM,NELEM,NOFVAR,NOFVERT,NPOIN,NWFAC
      CHARACTER FILE* (*)
C
      DOUBLE PRECISION GA
      PARAMETER(GA=1.4d0) 
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION COOR(NDIM,NPOIN),SKINF(NWFAC),Z(NOFVAR,NPOIN)
      INTEGER IBNDFAC(3,NBFAC),ICELCEL(NOFVERT,NELEM),
     +        ICELNOD(NOFVERT,NELEM),IPNTR(2,NWFAC)
C     ..
C     .. Local Scalars ..
      REAL*8 DIST,FUNY,S,T,TAUW,UTAU,X0,X1,X2,XA,XB,XC,XMAX,XMIN,Y0,Y1,
     +       Y2,YA,YB,YC,YMAX,YMIN,U,V,R,H,PRES,DX,DY,ANX,ANY,RR,YY,XX,
     +       TEMP,UX,UY
      INTEGER I,IDUMMY,IE,IFACE,INTERSECT,IPNT,IUNIT,IV,IVAR,IXDRS,J,JV,
     +        N1,N2,KUNIT,IFAIL
      LOGICAL COMPRESSIBLE,NORMAL_TO_THE_WALL
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION TMP(6)
C     ..
C     .. External Functions ..
      INTEGER GR_LINE_CROSS,ICYCL,INITXDR
      EXTERNAL GR_LINE_CROSS,ICYCL,INITXDR
C     ..
C     .. External Subroutines ..
      INTEGER IXDRIMAT,IXDRINT
      EXTERNAL IXDRIMAT,IXDRINT
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,MAX,MIN,SQRT
C     ..
C     .. Data statements ..

      DATA TMP/6*0.D0/
C     ..
C
C     compute pressure out of R,H,U,V
C
      PRES(R,H,U,V) = R*(GA-1.)/GA*(H-0.5*(U*U+V*V))


      IF (NDIM.NE.2) THEN
         WRITE(6,FMT=*) 'Routine CUT works only in 2D'
         CALL EXIT(1)
      ENDIF
      COMPRESSIBLE = (NOFVAR .EQ. 4)
      IF(COMPRESSIBLE)THEN
         WRITE(6,200)"Compressible"
      ELSE
         WRITE(6,200)"InCompressible"
      ENDIF
      WRITE(6,140)AA,"X",BB,"Y",CC,"="
      IF( NORMAL_TO_THE_WALL)WRITE(6,100)
C
C     IF( COMPRESSIBLE )THEN
C         U = 3
C     ELSE
C         U = 2
C     ENDIF
C     read neighbours from file FILE
C
      IXDRS = INITXDR(FILE,'r',.FALSE.)

      IFAIL = IXDRINT(IXDRS,IDUMMY)
      IFAIL = IXDRINT(IXDRS,IDUMMY)
      IFAIL = IXDRIMAT(IXDRS,NOFVERT*NELEM,ICELCEL)

      IPNT = 1
      XMIN = COOR(1,IPNT)
      XMAX = COOR(1,IPNT)
      YMIN = COOR(2,IPNT)
      YMAX = COOR(2,IPNT)
      DO 10 IPNT = 2,NPOIN
          XMIN = MIN(XMIN,COOR(1,IPNT))
          XMAX = MAX(XMAX,COOR(1,IPNT))
          YMIN = MIN(YMIN,COOR(2,IPNT))
          YMAX = MAX(YMAX,COOR(2,IPNT))
   10 CONTINUE
      XMIN = 1.05*XMIN
      XMAX = 1.05*XMAX
      YMIN = 1.05*YMIN
      YMAX = 1.05*YMAX
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
      RR = 2.d0 *MAX(ABS(XMAX-XMIN),ABS(YMAX-YMIN)) 
C
C     write(6,*)xa,ya,xb,yb
C
      IUNIT = 20
      KUNIT = 50
C
C     loop over all faces of the boundary to see which
C     one (even more than one) intersects the rake
C     X1,X2 are the coordinates of the vertices of the side
C
      DO 3 I = 1,NBFAC
C ... find the intersection of the rake with the boundary
          IF (IBNDFAC(3,I).NE.ICLR) GOTO 3
          IE = IBNDFAC(1,I)
          IV = IBNDFAC(2,I)
C ... nodes on the boundary
          N1 = ICELNOD(ICYCL(IV+1,NOFVERT),IE)
          N2 = ICELNOD(ICYCL(IV+2,NOFVERT),IE)
          X1 = COOR(1,N1)
          Y1 = COOR(2,N1)
          X2 = COOR(1,N2)
          Y2 = COOR(2,N2)
C ... check whether these intersect with the line
          INTERSECT = GR_LINE_CROSS(XA,YA,XB,YB,X1,Y1,X2,Y2,S,T)
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
          XSAVE = X0
          YSAVE = Y0
          WRITE (6,FMT=*) 'Neighbouring meshpoints are ',n1,n2
C
C     compute values at X0,Y0 by linear interpolation
C
          DO 4 IVAR = 1,NOFVAR
              TMP(IVAR) = (1.D0-S)*Z(IVAR,N1) + S*Z(IVAR,N2)
C             write(6,*)Z(IVAR,N1),Z(IVAR,N2),TMP(IVAR)
    4     CONTINUE
          IF(COMPRESSIBLE)THEN
C
C     transform from parameter vector into (r,H,u,v) 
C
                 TMP(4) = TMP(4)/TMP(1)
                 TMP(3) = TMP(3)/TMP(1)
                 TMP(2) = TMP(2)/TMP(1)
                 TMP(1) = TMP(1)*TMP(1)
C     save density ratio for Howarth Dorodnitsyn coordinate scaling
                 RA = TMP(1)
          ELSE
                 RA = 1.d0
          ENDIF
          IF(NORMAL_TO_THE_WALL)THEN
C     rotate velocity components
              IF(COMPRESSIBLE)THEN
                  UX = TMP(3)
                  UY = TMP(4)
                  TMP(3) = UX*DX +UY*DY
                  TMP(4) = UX*ANX+UY*ANY
              ELSE
                  UX = TMP(2)
                  UY = TMP(3)
                  TMP(2) = UX*DX +UY*DY
                  TMP(3) = UX*ANX+UY*ANY
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
              WRITE(6,FMT=*)'Wall friction velocity is ',utau
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
          IUNIT = IUNIT + 1
          KUNIT = KUNIT + 1
          WRITE (6,FMT=*) 'Writing flow variables to UNIT ',IUNIT
          WRITE (6,FMT=*) 'Writing b.l. variables to UNIT ',KUNIT
C
C    writes wall distance, y+, variables, u+ for the wall node
C
          IF( COMPRESSIBLE )THEN
              TMP(5) = PRES(TMP(1),TMP(2),TMP(3),TMP(4))
              TMP(6) = GA*M_INFTY*M_INFTY*TMP(5)/TMP(1)
          ELSE
              TMP(5) = 0.D0
              TMP(6) = 0.D0
          ENDIF 
          YBAR = 0.d0
          WRITE (IUNIT,FMT=80) X0,Y0,0.D0,YBAR,
     +(TMP(IVAR),IVAR=1,NOFVAR),TMP(5),TMP(6)
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
          WRITE(6,FMT=*) 'Uh! Oh! Smthg. went wrong in subroutine CUT'
          CALL EXIT(1)

    7     CONTINUE
C    compute the intersection point
          XC = (1.D0-S)*X1 + S*X2
          YC = (1.D0-S)*Y1 + S*Y2
          DELX = XC-XSAVE
          DELY = YC-YSAVE
          XSAVE = XC
          YSAVE = YC
C
C     write(6,*)xc,yc,utau
C
C     compute values at XC,YC by linear interpolation;
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
                 RB = TMP(1)
          ELSE
                 RB = 1.d0
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
C
C    YBAR = \int_{0}{y} \rho/\rho_{\infty} dy
C
          WRITE(6,*)0.5d0*(RA+RB)
          YBAR = YBAR + 0.5d0*(RA+RB)*DSQRT(DELX*DELX+DELY*DELY)
          DIST = SQRT((XC-X0)**2+ (YC-Y0)**2)
          RA = RB
C
C    really for b.l 
C
          ETA = DIST/XC*SQRT(0.5d0*REYNO*XC)
          ETABAR = YBAR/XC*SQRT(0.5d0*REYNO*XC)
C
C    writes wall distance, y+, variables, u+
C    >>> incompressible <<<
C
          IF( COMPRESSIBLE )THEN
C    compute pressure
              TMP(5) = PRES(TMP(1),TMP(2),TMP(3),TMP(4))
              TMP(6) = GA*M_INFTY*M_INFTY*TMP(5)/TMP(1)
              UX = TMP(3)
          ELSE
C    compute velocity magnitude
              TMP(5) = SQRT(TMP(2)*TMP(2)+TMP(3)*TMP(3))
              UX = TMP(2)
              TMP(6) = 0.d0
          ENDIF 
          WRITE (IUNIT,FMT=80) XC,YC,ETA,ETABAR,
     +      (TMP(IVAR),IVAR=1,NOFVAR),TMP(5),TMP(6)
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
          WRITE(6,FMT=*) 'Uh! Oh! Smthg. went wrong in subroutine CUT'
          CALL EXIT(1)

    3 CONTINUE
      RETURN

   80 FORMAT (10 (E16.9,1X))
  100 FORMAT (/10("*"),'Probing normal to the wall',/)
  140 FORMAT (/,'Probing along the straight line of equation ',/,
     &        3(E16.9,1X,1A,1X),'0',/)
  180 FORMAT(5X,'wall distance, y+, flow variables, u+')
  200 FORMAT(5X,'DATA FOR ',A,' FLOWS')

      END
