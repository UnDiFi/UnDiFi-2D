      SUBROUTINE CF(SKINF,QFLUX,TINDX,EDGNOD,IWORK,IRANK,LWORK,ICELNOD,
     +              COOR,ZROE,NDIM,NOFVERT,NOFVAR,REYNO,M_INFTY,
     +              FILENAME)
C
      IMPLICIT NONE
      INCLUDE 'constants'
C
C     .. Parameters ..
      INTEGER MBODIES
      PARAMETER (MBODIES=10)
      DOUBLE PRECISION TK
      PARAMETER (TK=0.41D0)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION M_INFTY,REYNO
      INTEGER LWORK,NDIM,NOFVAR,NOFVERT
      CHARACTER FILENAME* (*)
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION COOR(NDIM,*),SKINF(LWORK),TINDX(*),
     +ZROE(NOFVAR,*),QFLUX(LWORK)
      INTEGER EDGNOD(LWORK,2),ICELNOD(NOFVERT,*),IRANK(LWORK),
     +        IWORK(2,LWORK)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION A,ASQR,DENS,DNX,DNY,DX,DY,TAUW,UTAU,VISC,
     +                 XG,YDIST,YG,YPLUS,Z1,Z2,Z3,Z4,HEAT,X1,Y1,X2,Y2
      INTEGER I,IBODY,IC,IEDGE,IELEM,IFLAG,IFRST,ILAST,IPOIN,IUNIT,
     +        IVERT,IXDRS,IXDRS2,J,N1,N2,NBODY6,NEDGE,NWFAC,IFAIL
      LOGICAL COMPRESSIBLE,TURBO
C     ..
C     .. Local Arrays ..
      INTEGER IBGN(MBODIES+1),NFACE(MBODIES)
C     ..
C     .. External Functions ..
      DOUBLE PRECISION SUTHERLAW
      INTEGER INITXDR,JCYCL
      EXTERNAL SUTHERLAW,INITXDR,JCYCL
C     ..
C     .. External Subroutines ..
      INTEGER IXDRDMAT,IXDRIMAT,IXDRINT,IXDRCLOSE
      EXTERNAL IXDRDMAT,IXDRIMAT,IXDRINT,IXDRCLOSE
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,SQRT
C     ..
C
      IF (NDIM.NE.2) THEN
         WRITE(6,FMT=*) 'SUBROUTINE CF only works with NDIM=2'
         CALL EXIT(1)
      ENDIF
C
      COMPRESSIBLE = (NOFVAR.EQ.4)
C
      IUNIT = 30
C
      IXDRS = INITXDR(FILENAME,'r',.FALSE.)
C
C     file017.dat is the file containing the turbulent index
C
      INQUIRE(FILE='file017.dat',EXIST=TURBO)
      IF(TURBO)THEN
         IXDRS2 = INITXDR('file017.dat','r',.FALSE.)
      ELSE
         WRITE(6,*)'No turbulent data'
      ENDIF
C
      IFAIL = IXDRINT(IXDRS,NBODY6)
      IFAIL = IXDRINT(IXDRS,NWFAC)
C
      IBGN(1) = 1
      NWFAC = 0
C
      DO 10 IBODY = 1,NBODY6
C
C     number of faces for the current body
C
          IFAIL = IXDRINT(IXDRS,NFACE(IBODY))
C
          IBGN(IBODY+1) = IBGN(IBODY) + NFACE(IBODY)
      WRITE (6,FMT=*) 'No-slip body # ',IBODY,NFACE(IBODY),' faces'

          NWFAC = NWFAC + NFACE(IBODY)
   10 CONTINUE

      IFAIL = IXDRIMAT(IXDRS,2*NWFAC,IWORK)
      IFAIL = IXDRDMAT(IXDRS,NWFAC,SKINF)
      IFAIL = IXDRDMAT(IXDRS,NWFAC,QFLUX)
      IF(TURBO)THEN
          IFAIL = IXDRDMAT(IXDRS2,NWFAC,TINDX)
      ELSE
          CALL DINIT(NWFAC,ZERO,TINDX,1)
      ENDIF
         
C
C     loop over bodies
C
      WRITE (6,FMT=*) 'Ordering edges .......' 
      DO 15 IBODY = 1,NBODY6
          WRITE (6,FMT=*) 'No-slip body # ',IBODY,NFACE(IBODY),' faces'
          IUNIT = IUNIT + 1
          WRITE (6,FMT=*) 'Writing skin friction coefficient to unit =
     &',IUNIT
          WRITE (6,FMT=*) 'Reynolds number is ',REYNO
          WRITE (6,FMT=*) 'Freestream Mach is ',M_INFTY
          WRITE (6,FMT=*) 'x   y   Cf   Ch   y+   a'
C
          WRITE (IUNIT,FMT=*) '# Reynolds number is ',REYNO
          WRITE (IUNIT,FMT=*) '# Freestream Mach is ',M_INFTY
          WRITE (IUNIT,FMT=*) '# x   y   Cf   Ch   y+   a   yn'
          IFRST = IBGN(IBODY)
          ILAST = IBGN(IBODY+1) - 1
          IC = 0
          IFLAG = 0
          DO 17 J = IFRST,ILAST
              IELEM = IWORK(1,J)
              IVERT = IWORK(2,J)
              IC = IC + 1
C
C
              EDGNOD(J,1) = ICELNOD(JCYCL(IVERT+1),IELEM)
              EDGNOD(J,2) = ICELNOD(JCYCL(IVERT+2),IELEM)
   17     CONTINUE
C
C     first of all we assume that the body is closed
C     so we try to order edges in sequence;
C
C
C     NEDGE is the number of edges for the current body
C
          NEDGE = ILAST - IFRST + 1
C
C     start from the first edge (no matter which one it is)
C
          IC = IFRST
          IEDGE = 0
C
    4     IEDGE = IEDGE + 1
          IF (IEDGE.GT.NEDGE) GOTO 16
          IRANK(IEDGE) = IC
C
          IELEM = IWORK(1,IC)
          IVERT = IWORK(2,IC)
C
          N1 = EDGNOD(IC,1)
          N2 = EDGNOD(IC,2)
C
          IF (IFLAG.EQ.0) THEN
C
C     here we use the fact that, if all edges of a given body
C     are oriented in a consistent manner, i.e. if (N1,N2) are
C     the first and second node of the edge, then the vector
C     [-(Y(N2)-Y(N1)),(X(N2)-X(N1))] normal to the edge
C     points inside (or outside) the domain for ALL edges, then
C     node N1 will appear once as first and once as second node
C     in the edge to node list EDGNOD
C
              IPOIN = N2
C
C     find the edge having IPOIN in its first column
C
              DO 3 I = IFRST,ILAST
                  IF (IPOIN.EQ.EDGNOD(I,1)) THEN
                      IC = I
                      GOTO 4

                  ENDIF

    3         CONTINUE
C
C    If we get here it (probably) means that we have found
C    an open boundary (such a flat plate): therefore we
C    give up the idea of ordering the edges and proceed
C    following their natural order
C
              WRITE (6,FMT=*) IEDGE,
     +' apparently belongs to an open boundary '
              IFLAG = 1
              IC = IFRST
              IEDGE = 0
              GOTO 4

          ELSE
              IC = IC + 1
              GOTO 4

          ENDIF

   16     CONTINUE
C
C     now write according to the rank
C
          DO 19 J = 1,NEDGE
              IEDGE = IRANK(J)
              IELEM = IWORK(1,IEDGE)
              IVERT = IWORK(2,IEDGE)
C
              N1 = EDGNOD(IEDGE,1)
              N2 = EDGNOD(IEDGE,2)
              X1 = COOR(1,N1)
              Y1 = COOR(2,N1)
              X2 = COOR(1,N2)
              Y2 = COOR(2,N2)
C
C     baricentro del lato
C
              XG = (COOR(1,N1)+COOR(1,N2))/2.d0
              YG = (COOR(2,N1)+COOR(2,N2))/2.d0
C
C     normal to the face
C
              DNX = -COOR(2,N2) + COOR(2,N1)
              DNY = COOR(1,N2) - COOR(1,N1)
              A = ONE/SQRT(DNX*DNX+DNY*DNY)
              DNX = A*DNX
              DNY = A*DNY
C
C     distance from G to the vertex in front of the face
C
              IPOIN = ICELNOD(IVERT,IELEM)
              DX = COOR(1,IPOIN) - XG
              DY = COOR(2,IPOIN) - YG
C
              YDIST = ABS((DX*DNX + DY*DNY))
C
              TAUW = SKINF(IEDGE)
              HEAT = QFLUX(IEDGE)
C
              IF (COMPRESSIBLE) THEN
                  Z1 = HALF * (ZROE(1,N1)+ZROE(1,N2))
                  Z2 = HALF * (ZROE(2,N1)+ZROE(2,N2))
                  Z3 = HALF * (ZROE(3,N1)+ZROE(3,N2))
                  Z4 = HALF * (ZROE(4,N1)+ZROE(4,N2))
                  DENS = Z1*Z1
                  ASQR = GM1/DENS* (Z1*Z2-HALF* (Z3*Z3+Z4*Z4))
                  IF (ASQR.LT.ZERO) THEN
                      WRITE (6,FMT=*) 'Negative a^2 in edge ',IEDGE,ASQR
                  ENDIF

                  A = SQRT(ASQR)
                  VISC = SUTHERLAW(M_INFTY,A,ASQR)/DENS

              ELSE
                  VISC = ONE
                  DENS = ONE
              ENDIF
C
              UTAU = SQRT(ABS(TAUW)/DENS)
              YPLUS = REYNO*UTAU/VISC*YDIST
C
              A = TINDX(IEDGE)/ (REYNO*TK*UTAU)
C
              WRITE (IUNIT,FMT=100) XG,YG,2.d0*TAUW,HEAT,YPLUS,A,YDIST
!             WRITE (IUNIT,FMT=100) X1,Y1,2.d0*TAUW,HEAT,YPLUS,A,YDIST
!             WRITE (IUNIT,FMT=100) X2,Y2,2.d0*TAUW,HEAT,YPLUS,A,YDIST
   19     CONTINUE
   15 CONTINUE
C
      CLOSE (IUNIT)

      LWORK = NWFAC
      IFAIL = IXDRCLOSE( IXDRS )
      IF(TURBO)IFAIL = IXDRCLOSE( IXDRS2 )

! 100 FORMAT (F14.9,2X,7 (E12.6,1X))
  100 FORMAT (7(E12.6,1X))

      RETURN

      END
