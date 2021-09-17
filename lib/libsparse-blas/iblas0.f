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
      SUBROUTINE   GNINDX   ( NZ, N, ICLOBR, KINDX, INDX )
C
C     ==================================================================
C     ==================================================================
C     ====  GNINDX -- GENERATE INDEX ARRAY PATTERNS                 ====
C     ==================================================================
C     ==================================================================
C
C     GNINDX GENERATES VARIOUS PATTERNS FOR THE ARRAY INDX BASED
C     ON THE KEY KINDX.  THE GENERATED INDX ARRAY HAS NZ SIGNIFICANT
C     COMPONENTS.  THE REMAINING N-NZ COMPONENTS ARE SET TO 
C     ICLOBR.
C
C     ==================================================================
C
C     -------------
C     ... ARGUMENTS
C     -------------
C
      INTEGER             NZ, N, ICLOBR, KINDX, INDX (*)
C
C     -------------------
C     ... LOCAL VARIABLES
C     -------------------
C
      INTEGER             I,  L
C
C     --------------------
C     ... SUBPROGRAMS USED
C     --------------------
C
      EXTERNAL            IINIT
C
C     ==================================================================
C
      IF  ( N .LE. 0 )  RETURN
C
      L = MAX ( N, N-NZ )
      CALL IINIT ( L, ICLOBR, INDX, 1 )
C
      IF  ( NZ .LE. 0 )  RETURN
C
      KINDX = MAX ( KINDX, 1 )
      KINDX = MIN ( KINDX, 5 )
C
C     -------------------
C     ... BRANCH ON KINDX
C     -------------------
C
      GO TO ( 100, 200, 300, 400, 500 ), KINDX
C
C     -----------------------------------
C     ... ASCENDING ORDER - 1, 2, ..., NZ
C     -----------------------------------
C
  100 DO 110 I = 1, NZ
          INDX(I) = I
  110 CONTINUE
      GO TO 900
C
C     ------------------------------------------
C     ... ASCENDING ORDER - N-NZ+1, N-NZ, ..., N
C     ------------------------------------------
C
  200 L = N - NZ
      DO 210 I = 1, NZ
          INDX(I) = L + I
  210 CONTINUE
      GO TO 900
C
C     ---------------------------------------
C     ... DESCENDING ORDER - NZ, NZ-1, ..., 1
C     ---------------------------------------
C
  300 L = NZ
      DO 310 I = 1, NZ
          INDX(I) = L
          L       = L -1
  310 CONTINUE
      GO TO 900
C
C     ------------------------------------------
C     ... DESCENDING ORDER - N, N-1, ..., N-NZ+1
C     ------------------------------------------
C
  400 L = N
      DO 410 I = 1, NZ
          INDX(I) = L
          L       = L - 1
  410 CONTINUE
      GO TO 900
C
C     --------------------------------------------------------
C     ... ALTERNATING ORDER WITH EVEN NUMBERS IN REVERSE ORDER
C     --------------------------------------------------------
C
  500 DO 510 I = 1, NZ, 2
          INDX(I) = I
  510 CONTINUE
C
      L = N
      DO 520 I = 2, NZ, 2
          INDX(I) = L
          L       = L - 2
  520 CONTINUE
      GO TO 900
C
C     ==================================================================
C
  900 RETURN
      END
      LOGICAL FUNCTION  IVSAME   ( N, IX, IY )
C
C     ==================================================================
C
C     LOGICAL FUNCTION  IVSAME  DETERMINES IF THE VECTORS  IX  AND  IY
C     AGREE EXACTLY WITH EACH OTHER.
C
C     ==================================================================
C
C     ------------------------
C     ... VARIABLE DECLARATION
C     ------------------------
C
      INTEGER             I, N, IX (*), IY (*)
C
C     ==================================================================
C
      IVSAME = .TRUE.
C
      IF  ( N .LE. 0 )  RETURN
C
      DO 10 I = 1, N
          IF  ( IX(I) .NE. IY(I) )  THEN
              IVSAME = .FALSE.
              GO TO 20
          ENDIF
   10 CONTINUE
C
   20 RETURN
C
      END

