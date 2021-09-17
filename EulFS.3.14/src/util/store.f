      DOUBLE PRECISION FUNCTION STORE(X)
C
C***********************************************************
C
C                                               From TRIPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                             (817) 565-2767
C                                                   03/18/90
C
C   This function forces its argument X to be stored in a
C memory location, thus providing a means of determining
C floating point number characteristics (such as the machine
C precision) when it is necessary to avoid computation in
C high precision registers.
C
C
C On input:
C
C       X = Value to be stored.
C
C X is not altered by this function.
C
C On output:
C
C       STORE = Value of X after it has been stored and
C               possibly truncated or rounded to the single
C               precision word length.
C
C Modules required by STORE:  None
C
C***********************************************************
C
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION X
C     ..
C     .. Scalars in Common ..
      DOUBLE PRECISION Y
C     ..
C     .. Common blocks ..
      COMMON /STCOM/Y
C     ..
      Y = X
      STORE = Y
      RETURN

      END
