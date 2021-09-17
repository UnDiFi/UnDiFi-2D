      SUBROUTINE   ICOPY   ( N, X, INCX, Y, INCY )
C
C     ==================================================================
C     ==================================================================
C     ====  ICOPY -- COPY ONE INTEGER VECTOR TO ANOTHER             ====
C     ==================================================================
C     ==================================================================
C
C     PURPOSE ... (VARIANT OF 'SCOPY')
C                 COPY ONE INTEGER VECTOR TO ANOTHER.
C                 STANDARD INCREMENT OF 1 SHOULD BE USED FOR FORWARD
C                 COPY WITHIN SAME VECTOR.
C
C     CREATED       ... MAR. 12, 1985
C     LAST MODIFIED ... APR. 19, 1985
C
C     ==================================================================
C
C     -------------
C     ... ARGUMENTS
C     -------------
C
      INTEGER             N, INCX, INCY
C
      INTEGER             X (*), Y (*)
C
C     -------------------
C     ... LOCAL VARIABLES
C     -------------------
C
      INTEGER             XADDR, YADDR, I
C
C     ==================================================================
C
      IF  ( INCX .EQ. 1  .AND.  INCY .EQ. 1 )  THEN
C
C         -----------------------------------
C         ... UNIT INCREMENTS (STANDARD CASE)
C         -----------------------------------
C
          DO 100 I = 1, N
              Y (I) = X (I)
  100     CONTINUE
C
      ELSE
C
C         -------------------------
C         ... NON-UNIT INCREMENTS
C             (-1) USED FOR REVERSE
C             COPYING IN SAME ARRAY
C         -------------------------
C
          XADDR = 1
          YADDR = 1
C
          IF  ( INCX .LT. 0 )  THEN
              XADDR = (-N+1)*INCX + 1
          ENDIF
C
          IF  ( INCY .LT. 0 )  THEN
              YADDR = (-N+1)*INCY + 1
          ENDIF
C
          DO 200 I = 1, N
              Y (YADDR) = X (XADDR)
              XADDR     = XADDR + INCX
              YADDR     = YADDR + INCY
  200     CONTINUE
C
      ENDIF
C
      RETURN
C
      END
