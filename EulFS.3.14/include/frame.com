      DOUBLE PRECISION RotationMatrix(3,3)
      COMMON / STRFRM / RotationMatrix
C
C     RotationMatrix relates a stream-aligned ref. frame
C     to a cartesian one
C
