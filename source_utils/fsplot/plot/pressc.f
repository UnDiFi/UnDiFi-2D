      DOUBLE PRECISION FUNCTION PRESSC( NDIM , ZROE )
C
C    .. This function computes PRESSURE from Roe's
C       parameter vector ..
C
      IMPLICIT NONE
C
      INCLUDE 'constants'
C
      INTEGER NDIM
      DOUBLE PRECISION ZROE(*)
      DOUBLE PRECISION TEMP
C
      TEMP       = ZROE(3)*ZROE(3) + ZROE(4)*ZROE(4)
      IF( NDIM .EQ. 3 )TEMP = TEMP + ZROE(5)*ZROE(5)
      TEMP = HALF * TEMP
      PRESSC = GM1OG * ( ZROE(1)*ZROE(2) - TEMP )
C
      RETURN
      END
C
