      SUBROUTINE PARM2CONS(ZROE,DUDZ,NOFVAR,NDIM)
C
C     $Id: parm2consSA.f,v 1.2 2013/01/26 12:33:05 abonfi Exp $
C
      IMPLICIT NONE
C
C     transformation matrix from conserved variables to
C     parameter vector; version for the SA model 
C
C
C
C
C     Assembles the dUdZ matrix ...
C
C     .. Parameters ..
      INCLUDE 'constants.h'
C     ..
C     .. Scalar Arguments ..
      INTEGER NDIM,NOFVAR
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION DUDZ(NOFVAR,NOFVAR),ZROE(NOFVAR)
C     ..
C     .. Local Scalars ..
C     ..
      DUDZ(1,1) = TWO*ZROE(1)
      DUDZ(1,2) = ZERO
      DUDZ(1,3) = ZERO
      DUDZ(1,4) = ZERO
C
      DUDZ(2,1) = GINV*ZROE(2)
      DUDZ(2,2) = GINV*ZROE(1)
      DUDZ(2,3) = GM1OG*ZROE(3)
      DUDZ(2,4) = GM1OG*ZROE(4)
C
      DUDZ(3,1) = ZROE(3)
      DUDZ(3,2) = ZERO
      DUDZ(3,3) = ZROE(1)
      DUDZ(3,4) = ZERO
C
      DUDZ(4,1) = ZROE(4)
      DUDZ(4,2) = ZERO
      DUDZ(4,3) = ZERO
      DUDZ(4,4) = ZROE(1)
C
C     turbulent variable
C
      DUDZ(1,NOFVAR) = ZERO
      DUDZ(2,NOFVAR) = ZERO
      DUDZ(3,NOFVAR) = ZERO
      DUDZ(4,NOFVAR) = ZROE(1)
C
      DUDZ(NOFVAR,1) = ZAVG(NOFVAR)
      DUDZ(NOFVAR,2) = ZERO
      DUDZ(NOFVAR,3) = ZERO
      DUDZ(NOFVAR,4) = ZERO
      DUDZ(NOFVAR,NOFVAR) = ZAVG(1)
C
      IF (NDIM.EQ.2) RETURN
C
      DUDZ(1,5) = ZERO
      DUDZ(2,5) = GM1OG*ZROE(5)
      DUDZ(3,5) = ZERO
      DUDZ(4,5) = ZERO
C
      DUDZ(5,1) = ZROE(5)
      DUDZ(5,2) = ZERO
      DUDZ(5,3) = ZERO
      DUDZ(5,4) = ZERO
      DUDZ(5,5) = ZROE(1)
C
C     turbulent var
C
      DUDZ(5,NOFVAR) = ZERO
      DUDZ(NOFVAR,5) = ZERO
      RETURN

      END
