      SUBROUTINE   IINIT   ( N, A, X, INCX )
C
C     ==================================================================
C     ==================================================================
C     ====  IINIT -- INITIALIZE INTEGER VECTOR TO CONSTANT          ====
C     ==================================================================
C     ==================================================================
C
C     PURPOSE ... INITIALIZES INTEGER VECTOR TO A CONSTANT VALUE 'A'
C
C     CREATED       ... MAR. 8, 1985
C     LAST MODIFIED ... APR. 19, 1985
C
C     ==================================================================
C
C     -------------
C     ... ARGUMENTS
C     -------------
C
      INTEGER             N, INCX
C
      INTEGER             A, X (*)
C
C     -------------------
C     ... LOCAL VARIABLES
C     -------------------
C
      INTEGER             XADDR, I
C
C     ==================================================================
C
      IF  ( INCX .EQ. 1 )  THEN
C
C         ----------------------------------
C         ... UNIT INCREMENT (STANDARD CASE)
C         ----------------------------------
C
          DO 100 I = 1, N
              X(I) = A
  100     CONTINUE
C
      ELSE
C
C         ----------------------
C         ... NON-UNIT INCREMENT
C         ----------------------
C
          XADDR = 1
          IF  ( INCX .LT. 0 )  THEN
              XADDR = (-N+1)*INCX + 1
          ENDIF
C
          DO 200 I = 1, N
              X (XADDR) = A
              XADDR     = XADDR + INCX
  200     CONTINUE
C
      ENDIF
C
      RETURN
C
      END
