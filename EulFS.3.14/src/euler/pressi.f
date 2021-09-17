      DOUBLE PRECISION FUNCTION PRESSI(NDIM,ZROE)
C
C    .. Incompressible PRESSURE = ZROE(1,*)
C
C
C
C     .. Scalar Arguments ..
      INTEGER NDIM
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION ZROE(*)
C     ..
      PRESSI = ZROE(1)
C
      RETURN

      END
