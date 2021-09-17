      SUBROUTINE UNSTEADY(BETA,Z,NODRES,STIFC,NDIM,NOFVERT)
C
      IMPLICIT NONE
C
C     $Id:$
C
C
      INCLUDE 'paramt.h'
      INCLUDE 'constants.h'
      INCLUDE 'time.h'
      INCLUDE 'time.com'
C
C     .. Parameters ..
C     ..
C     .. Scalar Arguments ..
      INTEGER NDIM,NOFVERT
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION BETA(NOFVERT),Z(NOFVERT,*),NODRES(NOFVERT),
     +                 STIFC(NOFVERT,NOFVERT)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION RD,S,HELP
      INTEGER I,J,IADDR
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION MMAT(MAX_NOFVERT_SQR)
C     ..
C     .. External Functions ..
      DOUBLE PRECISION DDOT
      EXTERNAL DDOT
C     ..
C
      HELP = (ONE+HALF*GAMT)/(NOFVERT*DTVOL)
      RD = ONE/REAL(NOFVERT)
      GOTO (10,20,30,40,50) MMTYPE
   10 CONTINUE ! lumped mass matrix
      IADDR = 0
      DO 5 J = 1,NOFVERT
         DO 5 I = 1,NOFVERT
            IADDR = IADDR + 1 
            IF(J.EQ.I)THEN
               MMAT(IADDR) = ONE
            ELSE
               MMAT(IADDR) = ZERO
            ENDIF 
    5 CONTINUE
      GOTO 100
   20 CONTINUE ! Petrov-Galerkin
      IADDR = 0
      DO 3 J = 1,NOFVERT
         DO 3 I = 1,NOFVERT
            IADDR = IADDR + 1 
            IF(J.EQ.I)THEN
               MMAT(IADDR) = BETA(I)+ONE/6.d0 ! change for 3D
            ELSE
               MMAT(IADDR) = BETA(I)-ONE/12.d0 ! change for 3D
            ENDIF 
    3 CONTINUE
      GOTO 100
   30 CONTINUE ! Consistent Upwind
      IADDR = 0
      DO 1 J = 1,NOFVERT
         DO 1 I = 1,NOFVERT
            IADDR = IADDR + 1 
            IF(J.EQ.I)THEN
               MMAT(IADDR) = BETA(I)*(TWO-BETA(J))
            ELSE
               MMAT(IADDR) = BETA(I)*(ONE-BETA(J))
            ENDIF 
    1 CONTINUE
      GOTO 100
   40 CONTINUE ! Simple Upwind
      IADDR = 0
      DO 7 J = 1,NOFVERT
         DO 7 I = 1,NOFVERT
            IADDR = IADDR + 1 
            MMAT(IADDR) = BETA(I)
    7 CONTINUE
      GOTO 100
   50 CONTINUE ! Centred
      IADDR = 0
      DO 9 J = 1,NOFVERT
         DO 9 I = 1,NOFVERT
            IADDR = IADDR + 1 
            MMAT(IADDR) = RD
    9 CONTINUE
      GOTO 100
!
!         CALL X04CAF('G',' ',NOFVERT,NOFVERT,MMAT(1),
!    +                NOFVERT,'Mass matrix ',INFO)
!         CALL X04CAF('G',' ',NOFVERT,NOFVERT,STIFC(1,1),
!    +                NOFVERT,'C_ij matrix ',INFO)
  100 CONTINUE
C update the matrix by adding the mass matrix
      IADDR = 0
      DO J = 1,NOFVERT
         DO I = 1,NOFVERT
            IADDR = IADDR + 1 
            STIFC(I,J) = THETAT*STIFC(I,J)-HELP*MMAT(IADDR)
         ENDDO
      ENDDO
C
C update the rhs by adding the contribution
C from the previous time levels
C
      IF(GAMT.NE.ZERO)THEN
         HELP = (HALF*GAMT)/(NOFVERT*DTVOL)
C
C        compute u^n-u^{n-1}
C
         DO J = 1, NOFVERT
            Z(J,2) = Z(J,2) - Z(J,3) 
         ENDDO
         DO I = 1,NOFVERT
            S = ZERO 
            DO J = 1,NOFVERT
               IADDR = (J-1)*NOFVERT+I
               S = S + MMAT(IADDR)*Z(J,2)
            ENDDO
            NODRES(I) = NODRES(I) + HELP*S
         ENDDO
      ENDIF 
      RETURN
      END 
