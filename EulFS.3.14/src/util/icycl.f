      INTEGER FUNCTION ICYCL (I,NCYC)
C
C     //////////////////////////////////////////////////////////////
C
C     This function brings the value of I back into the interval
C     [1,NCYC] in a cyclic way.
C
C     For instance, if NCYC = 5
C
C      -7 -> 3
C      -6 -> 4
C      -5 -> 5
C      -4 -> 1
C      -3 -> 2
C      -2 -> 3
C      -1 -> 4
C       0 -> 5
C       1 -> 1
C       2 -> 2
C       3 -> 3
C       4 -> 4
C       5 -> 5
C       6 -> 1
C       7 -> 2
C       8 -> 3
C       9 -> 4
C
C     ..............................................................
C
C
      INTEGER I,NCYC
C
      INTEGER IM
C
      INTRINSIC MOD,ISIGN
C
      IM    = MOD(I,NCYC)
      ICYCL = IM + NCYC*((1-ISIGN(1,(IM-1)))/2)
C
      RETURN
      END
