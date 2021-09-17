      DOUBLE PRECISION TD,TTD
      COMMON/SACOM/TD,TTD
C
C     used in setupRHS and SA turbulence models
C
