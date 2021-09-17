      SUBROUTINE   ZCOPY   ( N, X, INCX, Y, INCY )
C
C     ==================================================================
C     ==================================================================
C     ====  ZCOPY -- COPY ONE DOUBLE COMPLEX VECTOR TO ANOTHER      ====
C     ==================================================================
C     ==================================================================
C
C     PURPOSE ... STANDARD  BLAS  
C                 COPY ONE DOUBLE COMPLEX VECTOR TO ANOTHER.
C                 STANDARD INCREMENT OF 1 SHOULD BE USED FOR FORWARD
C                 COPY WITHIN SAME VECTOR.
C
C     ==================================================================
C
C     -------------
C     ... ARGUMENTS
C     -------------
C
      INTEGER             N, INCX, INCY
C
      COMPLEX * 16        X (*), Y (*)
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
C     ==================================================================
C     ==================================================================
C     ====  ZINIT -- INITIALIZE DOUBLE COMPLEX VECTOR TO CONSTANT   ====
C     ==================================================================
C     ==================================================================
C
      SUBROUTINE   ZINIT   ( N, A, X, INCX )
C
C     ==================================================================
C
C     PURPOSE ... INITIALIZES DOUBLE COMPLEX VECTOR TO 
C                 A CONSTANT VALUE 'A'
C
C     CREATED ... APR. 14, 1987
C
C     ==================================================================
C
C     -------------
C     ... ARGUMENTS
C     -------------
C
      INTEGER             N, INCX
C
      COMPLEX * 16        A, X (*)
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

