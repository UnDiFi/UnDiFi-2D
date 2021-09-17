      SUBROUTINE TEST( X , Y , TOLER, IELEM , NVAR )
C
C     $Id: test.f,v 1.1 2002/10/17 02:37:02 abonfi Exp $
C
      IMPLICIT NONE
C
      DOUBLE PRECISION ONE
      PARAMETER(ONE=1.D0)
C
C     .. Scalar Arguments ..
C
      INTEGER IELEM,NVAR
      DOUBLE PRECISION TOLER
C
C     .. Array Arguments ..
C
      DOUBLE PRECISION X(*),Y(*)
C
C     .. Local Scalars ..
C
      LOGICAL WARN
      DOUBLE PRECISION S
      INTEGER IVAR
C
C     .. Intrinsic Functions ..
C
      INTRINSIC DLOG10
C
C     .. Executable Statements ..
C
      WARN = .FALSE.
C
      DO 30 IVAR = 1 , NVAR
         S = DABS( X(IVAR) - Y(IVAR) )
         IF( S .GT. TOLER )WARN = .TRUE.
   30 CONTINUE
      IF( WARN )THEN
      WRITE(6,*)
      WRITE(6,"('Flux divergence test in elem #.',I6)")IELEM
        DO 32 IVAR = 1 , NVAR
caldo   IF( DABS(Y(IVAR)) .GT. 1.E-10 )THEN
caldo      S = DABS(ONE-X(IVAR)/Y(IVAR))
caldo      IF( DABS(S) .LE. 1.D-15 )THEN
caldo         S = 0.D0
caldo      ELSE
caldo         S = DLOG10(S)
caldo      ENDIF
caldo   ELSE
caldo      S = X(IVAR)-Y(IVAR)
caldo   ENDIF
           S = X(IVAR)-Y(IVAR)
        WRITE(6,'(i5,3(e20.12))')IVAR,X(IVAR),Y(IVAR),S
   32   CONTINUE
C       PAUSE
      ENDIF
      RETURN
      END
C
