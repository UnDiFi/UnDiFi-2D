      DOUBLE PRECISION C1,C2,C3,C4
      COMMON /VISCOM/ C1,C2,C3,C4
C
C     Constants for the non-dimensional form
C     of sutherland's law:
C
C     viscl(a) = C1 * a^3 * C4 / (C2 * a^2 + C3 )
C
C     where a is the sound speed and the constant
C     depend upon the selected linearisation
