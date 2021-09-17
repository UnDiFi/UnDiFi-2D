      SUBROUTINE LAPLACE4MOTION(IELEM,VCN,NDIM,NOFVERT,
     +                          VISCL,STIFEL,VOLUME)
C
      IMPLICIT NONE
C
C     $Id: laplace4motion.f,v 1.1 2013/06/28 09:21:44 abonfi Exp $
C
      INCLUDE 'paramt.h'
      INCLUDE 'constants.h'
C
C     ..
C     .. Common blocks ..
C     ..
C
      INCLUDE 'flags.com'
      INCLUDE 'io.com'
      INCLUDE 'three.com'
C
C
C     .. Parameters ..
C     ..
C     .. Scalar Arguments ..
      INTEGER NDIM,NOFVAR,NOFVERT
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION VISCL
      DOUBLE PRECISION STIFEL(NOFVERT,NOFVERT),VCN(NDIM,*)
C     ..
C     .. Subroutine Arguments ..
C     ..
C     .. Arrays in Common ..
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION S,VOLUME,TEMPB
      INTEGER IELEM,IFAIL,IOFF,I,J,IPOIN,IVERT
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION TMPIJ(MAXNOFVERT,MAXNOFVERT)
C     ..
C     .. External Functions ..
      DOUBLE PRECISION DDOT
      EXTERNAL DDOT
C     ..
C     .. External Subroutines ..
      EXTERNAL DINIT
C     ..
C     .. Intrinsic Functions ..
C
C     ..
      TEMPB = VISCL/ (NDIM*NDIM*VOLUME)
C
      DO 1 J = 1,NOFVERT
          DO 1 I = 1,NOFVERT
              IF (J.LE.I) THEN
                  STIFEL(I,J) =-TEMPB*DDOT(NDIM,VCN(1,I),1,VCN(1,J),1)

              ELSE
                  STIFEL(I,J) = STIFEL(J,I)
              ENDIF

    1 CONTINUE
C
      RETURN
C
      END
