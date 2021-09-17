      SUBROUTINE   TSXPYI   ( NOUT,   EPSILN, THRESH, NZMAX2, 
     1                        NUMNZ,  NZVALU, NUMA,   AVALUE,
     2                        X,      XSAVE,  XTRUE,  Y,      YSAVE, 
     3                        YTRUE , INDX,   INDXT,  LIST,   ERRCNT, 
     4                        ERRMAX )
C
C     ==================================================================
C     ==================================================================
C     ====  TSXPYI  -- CERTIFY  SAXPYI                              ====
C     ==================================================================
C     ==================================================================
C
C     SUBROUTINE  TSXPYI  IS THE CERTIFICATION MODULE FOR THE SPARSE
C     BASIC LINEAR ALGEBRA SUBROUTINE MODULE  SAXPYI.
C
C     WRITTEN BY      ROGER G GRIMES
C                     APRIL 1987
C
C     ==================================================================
C
C     -------------
C     ... ARGUMENTS
C     -------------
C
      INTEGER             NOUT,   NZMAX2, NUMNZ,  NUMA,   ERRCNT,
     1                    ERRMAX
C
      INTEGER             NZVALU (*),  INDX (*),    INDXT (*),
     1                    LIST (*)
C
      REAL                EPSILN, THRESH
C
      REAL                AVALUE (*),  
     1                    X (*),       XSAVE (*),   XTRUE (*),
     2                    Y (*),       YSAVE (*),   YTRUE (*)
C
C     -------------------
C     ... LOCAL VARIABLES
C     -------------------
C
      REAL                A,      ATRUE,  CLOBBR
C
      INTEGER             COUNT,  I,      ICLOBR, J,      KA,
     1                    KINDX,  KNZ,    N,      NZ,     NZTRUE
C
      REAL                ERR,    S,      T
C
C     --------------------
C     ... SUBPROGRAMS USED
C     --------------------
C
      LOGICAL             IVSAME, SVSAME
C
      EXTERNAL            ICOPY,  SCOPY,  IINIT,  SINIT,  GNINDX, 
     1                    IVSAME, SVSAME, SAXPYI
C
C     ==================================================================
C
C     ------------------
C     ... INITIALIZATION
C     ------------------
C
      COUNT     =   0
C
      CLOBBR    =   -1.0E10
      ICLOBR    =   -10000000
C
C     ------------------------------------
C     ... GENERATE SOME VALUES FOR X AND Y
C     ------------------------------------
C
      DO 100 I = 1, NZMAX2
         XSAVE(I) = COS ( .6*FLOAT(I) )
         YSAVE(I) = SIN ( .7*FLOAT(I) )
  100 CONTINUE
C
C     ------------------------
C     ... FOR EACH VALUE OF NZ
C     ------------------------
C
      DO 700 KNZ = 1, NUMNZ
C
          NZTRUE = NZVALU(KNZ)
          N      = 2 * MAX ( NZTRUE, 1 )
C
C         -----------------------
C         ... FOR EACH VALUE OF A
C         -----------------------
C
          DO 600 KA = 1, NUMA
C
              ATRUE = AVALUE(KA)
C
C             -------------------------------
C             ... FOR EACH KIND OF INDX ARRAY
C             -------------------------------
C
              DO 500 KINDX = 1, 5
C
                  CALL GNINDX ( NZTRUE, N, ICLOBR, KINDX, INDXT )
C
                  CALL IINIT ( N, -1, LIST, 1 )
C
                  DO 150 I = 1, NZTRUE
                      LIST (INDXT(I)) = I
  150             CONTINUE
C
C                 -----------------------
C                 ... GENERATE INPUT DATA
C                 -----------------------
C
                  I = MIN ( N, N-NZTRUE )
                  J = N - I + 1
                  CALL SCOPY ( NZTRUE, XSAVE,  1, XTRUE, 1 )
                  CALL SINIT ( I,      CLOBBR, XTRUE(J), 1 )
                  CALL SINIT ( N,      CLOBBR, YTRUE, 1 )
C
                  DO 200 I = 1, NZTRUE
                      YTRUE (INDXT(I)) = YSAVE (INDXT(I))
  200             CONTINUE
C
C                 -------------------
C                 ... COPY TRUE INPUT
C                 -------------------
C
                  A  = ATRUE
                  NZ = NZTRUE
C
                  CALL SCOPY ( N, YTRUE, 1, Y, 1 )
                  CALL SCOPY ( N, XTRUE, 1, X, 1 )
                  CALL ICOPY ( N, INDXT, 1, INDX, 1 )
C
C                 --------------------------
C                 ... COMPUTE IN-LINE RESULT
C                 --------------------------
C
                  DO 300 I = 1, NZTRUE
                      YTRUE (INDXT(I)) = YTRUE (INDXT(I))  + 
     1                                   ATRUE * XTRUE(I)
  300             CONTINUE
C
C                 ---------------
C                 ... CALL SAXPYI
C                 ---------------
C
                  CALL SAXPYI ( NZ, A, X, INDX, Y )
C
C                 -----------------------------------------
C                 ... TEST ARGUMENTS OF SAXPYI THAT ARE NOT
C                     SUPPOSE TO CHANGE.
C                 -----------------------------------------
C
                  IF  ( NZ .NE. NZTRUE )  THEN
                      COUNT = COUNT + 1
                      IF  ( COUNT .LE. ERRMAX )  THEN 
                          WRITE ( NOUT, 1000 ) NZTRUE, ATRUE, KINDX,
     1                                         NZ
                      END IF
                  END IF
C
                  IF  ( A .NE. ATRUE )  THEN
                      COUNT = COUNT + 1
                      IF  ( COUNT .LE. ERRMAX )  THEN 
                          WRITE ( NOUT, 1100 ) NZTRUE, ATRUE, KINDX,
     1                                         A
                      END IF
                  END IF
C
                  IF  ( .NOT. SVSAME ( N, X, XTRUE ) )  THEN
                      COUNT = COUNT + 1
                      IF  ( COUNT .LE. ERRMAX )  THEN 
                          WRITE ( NOUT, 1200 ) NZTRUE, ATRUE, KINDX
                      END IF
                  END IF
C
                  IF  ( .NOT. IVSAME ( N, INDX, INDXT ) )  THEN
                      COUNT = COUNT + 1
                      IF  ( COUNT .LE. ERRMAX )  THEN 
                          WRITE ( NOUT, 1300 ) NZTRUE, ATRUE, KINDX
                      END IF
                  END IF
C
C                 ---------------------------
C                 ... TEST OUTPUT FROM SAXPYI
C                 ---------------------------
C
                  DO 400 J = 1, N
                      IF  ( LIST(J) .EQ. -1 )  THEN 
                          IF  ( Y(J) .NE. YTRUE(J) )  THEN
                              COUNT = COUNT + 1
                              IF  ( COUNT .LE. ERRMAX )  THEN 
                                  WRITE ( NOUT, 1400 ) NZTRUE, ATRUE, 
     1                                                 KINDX, J, 
     2                                                 Y(J), YTRUE(J)
                              END IF
                          END IF
C
                      ELSE
C
                          S   = ABS ( Y(J) - YTRUE(J) )
                          T   = ABS ( ATRUE) * ABS ( XTRUE (LIST(J)))  + 
     1                          ABS ( YSAVE(J))
                          ERR = S / ( EPSILN * T )
                          IF  ( ERR .GT. THRESH )  THEN
                              COUNT = COUNT + 1
                              IF  ( COUNT .LE. ERRMAX )  THEN 
                                  WRITE ( NOUT, 1500 ) NZTRUE, ATRUE, 
     1                                                 KINDX, J, Y(J),
     2                                                 YTRUE(J), ERR
                              END IF
                          END IF
C
                      END IF
C
  400             CONTINUE
C
  500         CONTINUE
C
  600     CONTINUE
C
  700 CONTINUE
C
C     ==================================================================
C
C     ------------------
C     ... END OF TESTING
C     ------------------
C
      ERRCNT = ERRCNT + COUNT
      IF  ( COUNT .NE. 0 )  GO TO 800
C
C     -----------------------------------
C     ... WRITE PASSED MESSAGE AND RETURN
C     -----------------------------------
C
      WRITE ( NOUT, 2700 )
      GO TO 900
C
C     -----------------------------------
C     ... WRITE FAILED MESSAGE AND RETURN
C     -----------------------------------
C
  800 WRITE ( NOUT, 2800 ) COUNT
C
C     ------------------------
C     ... END OF MODULE TSXPYI
C     ------------------------
C
  900 CONTINUE
      RETURN
C
C     ==================================================================
C
C     -----------
C     ... FORMATS
C     -----------
C
 1000 FORMAT ( 5X, 'SAXPYI ALTERED NZ FOR TEST WITH NZ = ', I5,
     1             ' A =', 1PE15.5,
     2             ' AND THE INDX TYPE NO. ', I5,
     3             '.  ALTERED VALUE OF NZ = ', I5 )
C
 1100 FORMAT ( 5X, 'SAXPYI ALTERED A FOR TEST WITH NZ = ', I5,
     1             ' A =', 1PE15.5,
     2             ' AND THE INDX TYPE NO. ', I5,
     3             '.  ALTERED VALUE OF A =', 1PE15.5 )
C
 1200 FORMAT ( 5X, 'SAXPYI ALTERED ARRAY X FOR TEST WITH NZ = ', I5,
     1             ' A =', 1PE15.5,
     2             ' AND THE INDX TYPE NO. ', I5 )
C
 1300 FORMAT ( 5X, 'SAXPYI ALTERED ARRAY INDX FOR TEST WITH NZ = ', I5,
     1             ' A =', 1PE15.5,
     2             ' AND THE INDX TYPE NO. ', I5 )
C
 1400 FORMAT ( 5X, 'SAXPYI OUTPUT ARRAY Y IS INCORRECT FOR TEST WITH ',
     1             'NZ = ', I5, ' A =', 1PE15.5,
     2             ' AND THE INDX TYPE NO. ', I5
     3        /5X, 'INCORRECT COMPONENT NO. ', I5, ' HAS VALUE =', 
     4             1PE15.5,
     5             ' TRUE VALUE =', 1PE15.5 )
C
 1500 FORMAT ( 5X, 'SAXPYI OUTPUT ARRAY Y IS INACCURATE FOR TEST WITH ',
     1             'NZ = ', I5, ' A =', 1PE15.5,
     2             ' AND THE INDX TYPE NO. ', I5
     3        /5X, 'INACCURATE COMPONENT NO. ', I5, ' HAS VALUE =', 
     4             1PE15.5, ' TRUE VALUE =',
     5             1PE15.5, ' ERROR = ', 1PE12.1 )
C
 2700 FORMAT ( /5X, 'SAXPYI PASSED ALL TESTS.' ) 
C
 2800 FORMAT ( /5X, 'SAXPYI FAILED', I10, ' TESTS.'  )
C
C     ==================================================================
C
      END
      SUBROUTINE   TSDOTI   ( NOUT,   EPSILN, THRESH, NZMAX2, 
     1                        NUMNZ,  NZVALU, 
     2                        X,      XSAVE,  XTRUE,  Y,      YSAVE, 
     3                        YTRUE , INDX,   INDXT,  ERRCNT, ERRMAX )
C
C     ==================================================================
C     ==================================================================
C     ====  TSDOTI  --  CERTIFY  SDOTI                             ====
C     ==================================================================
C     ==================================================================
C
C     SUBROUTINE  TSDOTI  IS THE CERTIFICATION MODULE FOR THE SPARSE
C     BASIC LINEAR ALGEBRA SUBROUTINE MODULE  SDOTI.
C
C     WRITTEN BY      ROGER G GRIMES
C                     APRIL 1987
C
C     ==================================================================
C
C     -------------
C     ... ARGUMENTS
C     -------------
C
      INTEGER             NOUT,   NZMAX2, NUMNZ,  ERRCNT,
     1                    ERRMAX
C
      INTEGER             NZVALU (*),  INDX (*),    INDXT (*)
C
      REAL                EPSILN, THRESH
C
      REAL                X (*),       XSAVE (*),   XTRUE (*),
     1                    Y (*),       YSAVE (*),   YTRUE (*)
C
C     -------------------
C     ... LOCAL VARIABLES
C     -------------------
C
      INTEGER             COUNT,  I,      ICLOBR, J,      KINDX,
     1                    KNZ,    N,      NZ,     NZTRUE
C
      REAL                ERR,    S,      T
C
      REAL                CLOBBR, V,      W
C
C     --------------------
C     ... SUBPROGRAMS USED
C     --------------------
C
      LOGICAL             IVSAME, SVSAME
C
      REAL                SDOTI
C
      EXTERNAL            ICOPY,  SCOPY,  SINIT,  GNINDX,
     1                    IVSAME, SVSAME, SDOTI
C
C     ==================================================================
C
C     ------------------
C     ... INITIALIZATION
C     ------------------
C
      COUNT     =   0
C
      CLOBBR    =   -1.0E10
      ICLOBR    =   -10000000
C
C     ------------------------------------
C     ... GENERATE SOME VALUES FOR X AND Y
C     ------------------------------------
C
      DO 100 I = 1, NZMAX2
         XSAVE(I) = COS ( .6*FLOAT(I) )
         YSAVE(I) = SIN ( .7*FLOAT(I) )
  100 CONTINUE
C
C     ------------------------
C     ... FOR EACH VALUE OF NZ
C     ------------------------
C
      DO 600 KNZ = 1, NUMNZ
C
          NZTRUE = NZVALU(KNZ)
          N      = 2 * MAX ( NZTRUE, 1 )
C
C         -------------------------------
C         ... FOR EACH KIND OF INDX ARRAY
C         -------------------------------
C
          DO 500 KINDX = 1, 5
C
              CALL GNINDX ( NZTRUE, N, ICLOBR, KINDX, INDXT )
C
C             -----------------------
C             ... GENERATE INPUT DATA
C             -----------------------
C
              I = MIN ( N, N-NZTRUE )
              J = N - I + 1
              CALL SCOPY ( NZTRUE, XSAVE,  1, XTRUE, 1 )
              CALL SINIT ( I,      CLOBBR, XTRUE(J), 1 )
              CALL SINIT ( N,      CLOBBR, YTRUE, 1 )
C
              DO 200 I = 1, NZTRUE
                  YTRUE (INDXT(I)) = YSAVE (INDXT(I))
  200         CONTINUE
C
C             -------------------
C             ... COPY TRUE INPUT
C             -------------------
C
              NZ = NZTRUE
C
              CALL SCOPY ( N, YTRUE, 1, Y, 1 )
              CALL SCOPY ( N, XTRUE, 1, X, 1 )
              CALL ICOPY ( N, INDXT, 1, INDX, 1 )
C
C             --------------------------
C             ... COMPUTE IN-LINE RESULT
C             --------------------------
C
              V = 0.0E0
C
              DO 300 I = 1, NZTRUE
                  V = V +  XTRUE(I) * YTRUE (INDXT(I))
  300         CONTINUE
C
C             --------------
C             ... CALL SDOTI
C             --------------
C
              W = SDOTI ( NZ, X, INDX, Y )
C
C             ----------------------------------------
C             ... TEST ARGUMENTS OF SDOTI THAT ARE NOT
C                     SUPPOSE TO CHANGE.
C             ----------------------------------------
C
              IF  ( NZ .NE. NZTRUE )  THEN
                  COUNT = COUNT + 1
                  IF  ( COUNT .LE. ERRMAX )  THEN 
                      WRITE ( NOUT, 1000 ) NZTRUE, KINDX, NZ
                  END IF
              END IF
C
              IF  ( .NOT. SVSAME ( N, X, XTRUE ) )  THEN
                  COUNT = COUNT + 1
                  IF  ( COUNT .LE. ERRMAX )  THEN 
                      WRITE ( NOUT, 1100 ) NZTRUE, KINDX
                  END IF
              END IF
C
              IF  ( .NOT. IVSAME ( N, INDX, INDXT ) )  THEN
                  COUNT = COUNT + 1
                  IF  ( COUNT .LE. ERRMAX )  THEN 
                      WRITE ( NOUT, 1200 ) NZTRUE, KINDX
                  END IF
              END IF
C
              IF  ( .NOT. SVSAME ( N, Y, YTRUE ) )  THEN
                  COUNT = COUNT + 1
                  IF  ( COUNT .LE. ERRMAX )  THEN 
                      WRITE ( NOUT, 1300 ) NZTRUE, KINDX
                  END IF
              END IF
C
C             --------------------------
C             ... TEST OUTPUT FROM SDOTI
C             --------------------------
C
              S = ABS ( V - W )
C
              T = 0.0E0
              DO 400 I = 1, NZTRUE
                  T = T + ABS ( XTRUE(I) * YTRUE (INDXT(I)) )
  400         CONTINUE
C
              IF  ( T .EQ. 0.0E0 )  T = 1.0E0
C
              ERR = S / ( EPSILN * T )
C
              IF  ( ERR .GT. THRESH )  THEN
                  COUNT = COUNT + 1
                  IF  ( COUNT .LE. ERRMAX )  THEN 
                      WRITE ( NOUT, 1400 ) NZTRUE, KINDX, 
     1                                     W, V, ERR
                  END IF
              END IF
C
  500     CONTINUE
C     
  600 CONTINUE
C
C     ==================================================================
C
C     ------------------
C     ... END OF TESTING
C     ------------------
C
      ERRCNT = ERRCNT + COUNT
      IF  ( COUNT .NE. 0 )  GO TO 800
C
C     -----------------------------------
C     ... WRITE PASSED MESSAGE AND RETURN
C     -----------------------------------
C
      WRITE ( NOUT, 2700 )
      GO TO 900
C
C     -----------------------------------
C     ... WRITE FAILED MESSAGE AND RETURN
C     -----------------------------------
C
  800 WRITE ( NOUT, 2800 ) COUNT
C
C     ------------------------
C     ... END OF MODULE TSDOTI
C     ------------------------
C
  900 CONTINUE
      RETURN
C
C     ==================================================================
C
C     -----------
C     ... FORMATS
C     -----------
C
 1000 FORMAT ( 5X, 'SDOTI ALTERED NZ FOR TEST WITH NZ = ', I5,
     1             ' AND THE INDX TYPE NO. ', I5,
     2             '.  ALTERED VALUE OF NZ = ', I5 )
C
 1100 FORMAT ( 5X, 'SDOTI ALTERED ARRAY X FOR TEST WITH NZ = ', I5,
     1             ' AND THE INDX TYPE NO. ', I5 )
C
 1200 FORMAT ( 5X, 'SDOTI ALTERED ARRAY INDX FOR TEST WITH NZ = ', I5,
     1             ' AND THE INDX TYPE NO. ', I5 )
C
 1300 FORMAT ( 5X, 'SDOTI ALTERED ARRAY Y FOR TEST WITH NZ = ', I5,
     1             ' AND THE INDX TYPE NO. ', I5 )
C
 1400 FORMAT ( 5X, 'SDOTI OUTPUT W IS INACCURATE FOR TEST WITH ',
     1             'NZ = ', I5, ' AND THE INDX TYPE NO. ', I5
     2        /5X, 'SDOTI HAS VALUE =', 1PE15.5,
     3             ' TRUE VALUE =', 1PE15.5,
     4             ' ERROR = ', 1PE12.1 )
C
 2700 FORMAT ( /5X, 'SDOTI  PASSED ALL TESTS.' ) 
C
 2800 FORMAT ( /5X, 'SDOTI  FAILED', I10, ' TESTS.'  )
C
C     ==================================================================
C
      END
      SUBROUTINE   TSGTHR   ( NOUT,   NZMAX2, NUMNZ,  NZVALU, 
     1                        X,      XSAVE,  XTRUE,  Y,      YSAVE, 
     2                        YTRUE , INDX,   INDXT,  ERRCNT, ERRMAX )
C
C     ==================================================================
C     ==================================================================
C     ====  TSGTHR  --  CERTIFY  SGTHR                              ====
C     ==================================================================
C     ==================================================================
C
C     SUBROUTINE  TSGTHR  IS THE CERTIFICATION MODULE FOR THE SPARSE
C     BASIC LINEAR ALGEBRA SUBROUTINE MODULE  SGTHR.
C
C     WRITTEN BY      ROGER G GRIMES
C                     APRIL 1987
C
C     ==================================================================
C
C     -------------
C     ... ARGUMENTS
C     -------------
C
      INTEGER             NOUT,   NZMAX2, NUMNZ,  ERRCNT,
     1                    ERRMAX
C
      INTEGER             NZVALU (*),  INDX (*),    INDXT (*)
C
      REAL                X (*),       XSAVE (*),   XTRUE (*),
     1                    Y (*),       YSAVE (*),   YTRUE (*)
C
C     -------------------
C     ... LOCAL VARIABLES
C     -------------------
C
      INTEGER             COUNT,  I,      ICLOBR, KINDX,
     1                    KNZ,    N,      NZ,     NZTRUE
C
      REAL                CLOBBR   
C
C     --------------------
C     ... SUBPROGRAMS USED
C     --------------------
C
      LOGICAL             IVSAME, SVSAME
C
      EXTERNAL            ICOPY,  SCOPY,  SINIT,  GNINDX,
     1                    IVSAME, SVSAME, SGTHR
C
C     ==================================================================
C
C     ------------------
C     ... INITIALIZATION
C     ------------------
C
      COUNT     =   0
C
      CLOBBR    =   -1.0E10
      ICLOBR    =   -10000000
C
C     ------------------------------------
C     ... GENERATE SOME VALUES FOR X AND Y
C     ------------------------------------
C
      DO 100 I = 1, NZMAX2
         XSAVE(I) = COS ( .6*FLOAT(I) )
         YSAVE(I) = SIN ( .7*FLOAT(I) )
  100 CONTINUE
C
C     ------------------------
C     ... FOR EACH VALUE OF NZ
C     ------------------------
C
      DO 600 KNZ = 1, NUMNZ
C
          NZTRUE = NZVALU(KNZ)
          N      = 2 * MAX ( NZTRUE, 1 )
C
C         -------------------------------
C         ... FOR EACH KIND OF INDX ARRAY
C         -------------------------------
C
          DO 500 KINDX = 1, 5
C
              CALL GNINDX ( NZTRUE, N, ICLOBR, KINDX, INDXT )
C
C             -----------------------
C             ... GENERATE INPUT DATA
C             -----------------------
C
              CALL SINIT ( N, CLOBBR, XTRUE, 1 )
              CALL SINIT ( N, CLOBBR, YTRUE, 1 )
C
              DO 200 I = 1, NZTRUE
                  YTRUE (INDXT(I)) = YSAVE (INDXT(I))
  200         CONTINUE
C
C             -------------------
C             ... COPY TRUE INPUT
C             -------------------
C
              NZ = NZTRUE
C
              CALL SCOPY ( N, YTRUE, 1, Y, 1 )
              CALL SCOPY ( N, XTRUE, 1, X, 1 )
              CALL ICOPY ( N, INDXT, 1, INDX, 1 )
C         
C             --------------------------
C             ... COMPUTE IN-LINE RESULT
C             --------------------------
C
              DO 300 I = 1, NZTRUE
                  XTRUE (I) = YTRUE (INDXT(I))
  300         CONTINUE
C
C             --------------
C             ... CALL SGTHR
C             --------------
C
              CALL SGTHR ( NZ, Y, X, INDX )
C
C             ----------------------------------------
C             ... TEST ARGUMENTS OF SGTHR THAT ARE NOT
C                 SUPPOSE TO CHANGE.
C             ----------------------------------------
C
              IF  ( NZ .NE. NZTRUE )  THEN
                  COUNT = COUNT + 1
                  IF  ( COUNT .LE. ERRMAX )  THEN 
                      WRITE ( NOUT, 1000 ) NZTRUE, KINDX, NZ
                  END IF
              END IF
C
              IF  ( .NOT. SVSAME ( N, Y, YTRUE ) )  THEN
                  COUNT = COUNT + 1
                  IF  ( COUNT .LE. ERRMAX )  THEN 
                      WRITE ( NOUT, 1100 ) NZTRUE, KINDX
                  END IF
              END IF
C
              IF  ( .NOT. IVSAME ( N, INDX, INDXT ) )  THEN
                  COUNT = COUNT + 1
                  IF  ( COUNT .LE. ERRMAX )  THEN 
                      WRITE ( NOUT, 1200 ) NZTRUE, KINDX
                  END IF
              END IF
C
C             --------------------------
C             ... TEST OUTPUT FROM SGTHR
C             --------------------------
C
              DO 400 I = 1, N
                  IF  ( X(I) .NE. XTRUE(I) )  THEN
                      COUNT = COUNT + 1
                      IF  ( COUNT .LE. ERRMAX )  THEN 
                          WRITE ( NOUT, 1300 ) NZTRUE, KINDX, I, 
     1                                         X(I), XTRUE(I)
                      END IF
                  END IF
  400         CONTINUE
C
  500     CONTINUE
C     
  600 CONTINUE
C
C     ==================================================================
C
C     ------------------
C     ... END OF TESTING
C     ------------------
C
      ERRCNT = ERRCNT + COUNT
      IF  ( COUNT .NE. 0 )  GO TO 800
C
C     -----------------------------------
C     ... WRITE PASSED MESSAGE AND RETURN
C     -----------------------------------
C
      WRITE ( NOUT, 2700 )
      GO TO 900
C
C     -----------------------------------
C     ... WRITE FAILED MESSAGE AND RETURN
C     -----------------------------------
C
  800 WRITE ( NOUT, 2800 ) COUNT
C
C     ------------------------
C     ... END OF MODULE TSGTHR
C     ------------------------
C
  900 CONTINUE
      RETURN
C
C     ==================================================================
C
C     -----------
C     ... FORMATS
C     -----------
C
 1000 FORMAT ( 5X, 'SGTHR ALTERED NZ FOR TEST WITH NZ = ', I5,
     1             ' AND THE INDX TYPE NO. ', I5,
     2             '.  ALTERED VALUE OF NZ = ', I5 )
C
 1100 FORMAT ( 5X, 'SGTHR ALTERED ARRAY Y FOR TEST WITH NZ = ', I5,
     1             ' AND THE INDX TYPE NO. ', I5 )
C
 1200 FORMAT ( 5X, 'SGTHR ALTERED ARRAY INDX FOR TEST WITH NZ = ', I5,
     1             ' AND THE INDX TYPE NO. ', I5 )
C
 1300 FORMAT ( 5X, 'SGTHR OUTPUT ARRAY X IS INCORRECT FOR TEST WITH ',
     1             'NZ = ', I5, ' AND THE INDX TYPE NO. ', I5
     2        /5X, 'INACCURATE COMPONENT NO. ', I5, ' HAS VALUE =', 
     3             1PE15.5, ' TRUE VALUE = ', 1PE15.5 )
C
 2700 FORMAT ( /5X, 'SGTHR  PASSED ALL TESTS.' ) 
C
 2800 FORMAT ( /5X, 'SGTHR  FAILED', I10, ' TESTS.'  )
C
C     ==================================================================
C
      END
      SUBROUTINE   TSGTHZ   ( NOUT,   NZMAX2, NUMNZ,  NZVALU, 
     1                        X,      XSAVE,  XTRUE,  Y,      YSAVE, 
     2                        YTRUE , INDX,   INDXT,  ERRCNT, ERRMAX )
C
C     ==================================================================
C     ==================================================================
C     ====  TSGTHZ  --  CERTIFY  SGTHRZ                             ====
C     ==================================================================
C     ==================================================================
C
C     SUBROUTINE  TSGTHZ  IS THE CERTIFICATION MODULE FOR THE SPARSE
C     BASIC LINEAR ALGEBRA SUBROUTINE MODULE  SGTHRZ.
C
C     WRITTEN BY      ROGER G GRIMES
C                     APRIL 1987
C
C     ==================================================================
C
C     -------------
C     ... ARGUMENTS
C     -------------
C
      INTEGER             NOUT,   NZMAX2, NUMNZ,  ERRCNT,
     1                    ERRMAX
C
      INTEGER             NZVALU (*),  INDX (*),    INDXT (*)
C
      REAL                X (*),       XSAVE (*),   XTRUE (*),
     1                    Y (*),       YSAVE (*),   YTRUE (*)
C
C     -------------------
C     ... LOCAL VARIABLES
C     -------------------
C
      INTEGER             COUNT,  I,      ICLOBR, KINDX,
     1                    KNZ,    N,      NZ,     NZTRUE
C
      REAL                CLOBBR   
C
C     --------------------
C     ... SUBPROGRAMS USED
C     --------------------
C
      LOGICAL             IVSAME, SVSAME
C
      EXTERNAL            ICOPY,  SCOPY,  SINIT,  GNINDX,
     1                    IVSAME, SVSAME, SGTHRZ
C
C     ==================================================================
C
C     ------------------
C     ... INITIALIZATION
C     ------------------
C
      COUNT     =   0
C
      CLOBBR    =   -1.0E10
      ICLOBR    =   -10000000
C
C     ------------------------------------
C     ... GENERATE SOME VALUES FOR X AND Y
C     ------------------------------------
C
      DO 100 I = 1, NZMAX2
         XSAVE(I) = COS ( .6*FLOAT(I) )
         YSAVE(I) = SIN ( .7*FLOAT(I) )
  100 CONTINUE
C
C     ------------------------
C     ... FOR EACH VALUE OF NZ
C     ------------------------
C
      DO 600 KNZ = 1, NUMNZ
C
          NZTRUE = NZVALU(KNZ)
          N      = 2 * MAX ( NZTRUE, 1 )
C
C         -------------------------------
C         ... FOR EACH KIND OF INDX ARRAY
C         -------------------------------
C
          DO 500 KINDX = 1, 5
C
              CALL GNINDX ( NZTRUE, N, ICLOBR, KINDX, INDXT )
C
C             -----------------------
C             ... GENERATE INPUT DATA
C             -----------------------
C
              CALL SINIT ( N, CLOBBR, XTRUE, 1 )
              CALL SINIT ( N, CLOBBR, YTRUE, 1 )
C
              DO 200 I = 1, NZTRUE
                  YTRUE (INDXT(I)) = YSAVE (INDXT(I))
  200         CONTINUE
C
C             -------------------
C             ... COPY TRUE INPUT
C             -------------------
C
              NZ = NZTRUE
C
              CALL SCOPY ( N, YTRUE, 1, Y, 1 )
              CALL SCOPY ( N, XTRUE, 1, X, 1 )
              CALL ICOPY ( N, INDXT, 1, INDX, 1 )
C         
C             --------------------------
C             ... COMPUTE IN-LINE RESULT
C             --------------------------
C
              DO 300 I = 1, NZTRUE
                  XTRUE (I) = YTRUE (INDXT(I))
                  YTRUE(INDXT(I)) = 0.0E0
  300         CONTINUE
C
C             ---------------
C             ... CALL SGTHRZ
C             ---------------
C
              CALL SGTHRZ ( NZ, Y, X, INDX )
C
C             -----------------------------------------
C             ... TEST ARGUMENTS OF SGTHRZ THAT ARE NOT
C                 SUPPOSE TO CHANGE.
C             -----------------------------------------
C
              IF  ( NZ .NE. NZTRUE )  THEN
                  COUNT = COUNT + 1
                  IF  ( COUNT .LE. ERRMAX )  THEN 
                      WRITE ( NOUT, 1000 ) NZTRUE, KINDX, NZ
                  END IF
              END IF
C
              IF  ( .NOT. IVSAME ( N, INDX, INDXT ) )  THEN
                  COUNT = COUNT + 1
                  IF  ( COUNT .LE. ERRMAX )  THEN 
                      WRITE ( NOUT, 1100 ) NZTRUE, KINDX
                  END IF
              END IF
C
C             ---------------------------
C             ... TEST OUTPUT FROM SGTHRZ
C             ---------------------------
C
              DO 400 I = 1, N
C
                  IF  ( X(I) .NE. XTRUE(I) )  THEN
                      COUNT = COUNT + 1
                      IF  ( COUNT .LE. ERRMAX )  THEN 
                          WRITE ( NOUT, 1200 ) NZTRUE, KINDX, I, 
     1                                         X(I), XTRUE(I)
                      END IF
                  END IF
C
                  IF  ( Y(I) .NE. YTRUE(I) )  THEN
                      COUNT = COUNT + 1
                      IF  ( COUNT .LE. ERRMAX )  THEN 
                          WRITE ( NOUT, 1300 ) NZTRUE, KINDX, I, 
     1                                         Y(I), YTRUE(I)
                      END IF
                  END IF
C
  400         CONTINUE
C
  500     CONTINUE
C     
  600 CONTINUE
C
C     ==================================================================
C
C     ------------------
C     ... END OF TESTING
C     ------------------
C
      ERRCNT = ERRCNT + COUNT
      IF  ( COUNT .NE. 0 )  GO TO 800
C
C     -----------------------------------
C     ... WRITE PASSED MESSAGE AND RETURN
C     -----------------------------------
C
      WRITE ( NOUT, 2700 )
      GO TO 900
C
C     -----------------------------------
C     ... WRITE FAILED MESSAGE AND RETURN
C     -----------------------------------
C
  800 WRITE ( NOUT, 2800 ) COUNT
C
C     ------------------------
C     ... END OF MODULE TSGTHZ
C     ------------------------
C
  900 CONTINUE
      RETURN
C
C     ==================================================================
C
C     -----------
C     ... FORMATS
C     -----------
C
 1000 FORMAT ( 5X, 'SGTHRZ ALTERED NZ FOR TEST WITH NZ = ', I5,
     1             ' AND THE INDX TYPE NO. ', I5,
     2             '.  ALTERED VALUE OF NZ = ', I5 )
C
 1100 FORMAT ( 5X, 'SGTHRZ ALTERED ARRAY INDX FOR TEST WITH NZ = ', I5,
     1             ' AND THE INDX TYPE NO. ', I5 )
C
 1200 FORMAT ( 5X, 'SGTHRZ OUTPUT ARRAY X IS INCORRECT FOR TEST WITH ',
     1             'NZ = ', I5, ' AND THE INDX TYPE NO. ', I5
     2        /5X, 'INACCURATE COMPONENT NO. ', I5, ' HAS VALUE =', 
     3             1PE15.5, ' TRUE VALUE =', 1PE15.5 )
C
 1300 FORMAT ( 5X, 'SGTHRZ OUTPUT ARRAY Y IS INCORRECT FOR TEST WITH ',
     1             'NZ = ', I5, ' AND THE INDX TYPE NO. ', I5
     2        /5X, 'INACCURATE COMPONENT NO. ', I5, ' HAS VALUE =', 
     3             1PE15.5, ' TRUE VALUE =', 1PE15.5 )
C
 2700 FORMAT ( /5X, 'SGTHRZ PASSED ALL TESTS.' ) 
C
 2800 FORMAT ( /5X, 'SGTHRZ FAILED', I10, ' TESTS.'  )
C
C     ==================================================================
C
      END
      SUBROUTINE   TSROTI   ( NOUT,   EPSILN, THRESH, NZMAX2, 
     1                        NUMNZ,  NZVALU, NUMG,   CVALUE, SVALUE,
     2                        X,      XSAVE,  XTRUE,  Y,      YSAVE, 
     3                        YTRUE , INDX,   INDXT,  LIST,   ERRCNT, 
     4                        ERRMAX )
C
C     ==================================================================
C     ==================================================================
C     ====  TSROTI  --  CERTIFY  SROTI                              ====
C     ==================================================================
C     ==================================================================
C
C     SUBROUTINE  TSROTI  IS THE CERTIFICATION MODULE FOR THE SPARSE
C     BASIC LINEAR ALGEBRA SUBROUTINE MODULE  SROTI.
C
C     WRITTEN BY      ROGER G GRIMES
C                     APRIL 1987
C
C     ==================================================================
C
C     -------------
C     ... ARGUMENTS
C     -------------
C
      INTEGER             NOUT,   NZMAX2, NUMNZ,  NUMG,   ERRCNT,
     1                    ERRMAX
C
      INTEGER             NZVALU (*),  INDX (*),    INDXT (*),
     1                    LIST (*)
C
      REAL                EPSILN, THRESH
C
      REAL                CVALUE (*),  SVALUE (*),
     1                    X (*),       XSAVE (*),   XTRUE (*),
     2                    Y (*),       YSAVE (*),   YTRUE (*)
C
C     -------------------
C     ... LOCAL VARIABLES
C     -------------------
C
      INTEGER             COUNT,  I,      ICLOBR, J,      KG,
     1                    KINDX,  KNZ,    N,      NZ,     NZTRUE
C
      REAL                C,      CLOBBR, CTRUE,  ERR,    S,     
     1                    STRUE,  V,      W
C
C     --------------------
C     ... SUBPROGRAMS USED
C     --------------------
C
      LOGICAL             IVSAME
C
      EXTERNAL            SCOPY,  SINIT,  ICOPY,  IINIT,  GNINDX,
     1                    IVSAME, SROTI
C
C     ==================================================================
C
C     ------------------
C     ... INITIALIZATION
C     ------------------
C
      COUNT     =   0
C
      CLOBBR    =   -1.0E10
      ICLOBR    =   -10000000
C
C     ------------------------------------
C     ... GENERATE SOME VALUES FOR X AND Y
C     ------------------------------------
C
      DO 100 I = 1, NZMAX2
         XSAVE(I) = COS ( .6E0 * FLOAT(I) )
         YSAVE(I) = SIN ( .7E0 * FLOAT(I) )
  100 CONTINUE
C
C     ------------------------
C     ... FOR EACH VALUE OF NZ
C     ------------------------
C
      DO 700 KNZ = 1, NUMNZ
C
          NZTRUE = NZVALU(KNZ)
          N      = 2 * MAX ( NZTRUE, 1 )
C
C         -----------------------------
C         ... FOR EACH VALUE OF C AND S
C         -----------------------------
C
          DO 600 KG = 1, NUMG
C
              CTRUE = CVALUE(KG)
              STRUE = SVALUE(KG)
C
C             -------------------------------
C             ... FOR EACH KIND OF INDX ARRAY
C             -------------------------------
C
              DO 500 KINDX = 1, 5
C
                  CALL GNINDX ( NZTRUE, N, ICLOBR, KINDX, INDXT )
C
                  CALL IINIT ( N, -1, LIST, 1 )
C
                  DO 150 I = 1, NZTRUE
                      LIST (INDXT(I)) = I
  150             CONTINUE
C
C                 -----------------------
C                 ... GENERATE INPUT DATA
C                 -----------------------
C
                  I = MIN ( N, N-NZTRUE )
                  J = N - I + 1
                  CALL SCOPY ( NZTRUE, XSAVE,  1, XTRUE, 1 )
                  CALL SINIT ( I,      CLOBBR, XTRUE(J), 1 )
                  CALL SINIT ( N,      CLOBBR, YTRUE   , 1 )
C
                  DO 200 I = 1, NZTRUE
                      YTRUE (INDXT(I)) = YSAVE (INDXT(I))
  200             CONTINUE
C
C                 -------------------
C                 ... COPY TRUE INPUT
C                 -------------------
C
                  C  = CTRUE
                  S  = STRUE
                  NZ = NZTRUE
C
                  CALL SCOPY ( N, YTRUE, 1, Y, 1 )
                  CALL SCOPY ( N, XTRUE, 1, X, 1 )
                  CALL ICOPY ( N, INDXT, 1, INDX, 1  )
C
C                 --------------------------
C                 ... COMPUTE IN-LINE RESULT
C                 --------------------------
C
                  DO 300 I = 1, NZTRUE
                      V                = XTRUE(I)
                      XTRUE(I)         =  CTRUE * V  + 
     1                                    STRUE * YTRUE (INDXT(I))
                      YTRUE (INDXT(I)) = -STRUE * V  + 
     1                                    CTRUE * YTRUE (INDXT(I))
  300             CONTINUE
C
C                 --------------
C                 ... CALL SROTI
C                 --------------
C
                  CALL SROTI ( NZ, X, INDX, Y, C, S )
C
C                 ----------------------------------------
C                 ... TEST ARGUMENTS OF SROTI THAT ARE NOT
C                     SUPPOSE TO CHANGE.
C                 ----------------------------------------
C
                  IF  ( NZ .NE. NZTRUE )  THEN
                      COUNT = COUNT + 1
                      IF  ( COUNT .LE. ERRMAX )  THEN 
                          WRITE ( NOUT, 1000 ) NZTRUE, CTRUE, STRUE, 
     1                                         KINDX,  NZ
                      END IF
                  END IF
C
                  IF  ( C .NE. CTRUE )  THEN
                      COUNT = COUNT + 1
                      IF  ( COUNT .LE. ERRMAX )  THEN 
                          WRITE ( NOUT, 1100 ) NZTRUE, CTRUE, STRUE,
     1                                         KINDX,  C,     S
                      END IF
                  END IF
C
                  IF  ( S .NE. STRUE )  THEN
                      COUNT = COUNT + 1
                      IF  ( COUNT .LE. ERRMAX )  THEN 
                          WRITE ( NOUT, 1200 ) NZTRUE, CTRUE, STRUE,
     1                                         KINDX,  C,     S
                      END IF
                  END IF
C
                  IF  ( .NOT. IVSAME ( N, INDX, INDXT ) )  THEN
                      COUNT = COUNT + 1
                      IF  ( COUNT .LE. ERRMAX )  THEN 
                          WRITE ( NOUT, 1300 ) NZTRUE, CTRUE, STRUE, 
     1                                         KINDX
                      END IF
                  END IF
C
C                 --------------------------
C                 ... TEST OUTPUT FROM SROTI
C                 --------------------------
C
                  DO 400 J = 1, N
C
                      IF  ( LIST(J) .EQ. -1 )  THEN 
C
                          IF  ( X(J) .NE. XTRUE(J) )  THEN
                              COUNT = COUNT + 1
                              IF  ( COUNT .LE. ERRMAX )  THEN 
                                  WRITE ( NOUT, 1400 ) NZTRUE, CTRUE, 
     1                                                 STRUE, KINDX, J, 
     2                                                 X(J), XTRUE(J)
                              END IF
                          END IF
C
                          IF  ( Y(J) .NE. YTRUE(J) )  THEN
                              COUNT = COUNT + 1
                              IF  ( COUNT .LE. ERRMAX )  THEN 
                                  WRITE ( NOUT, 1500 ) NZTRUE, CTRUE, 
     1                                                 STRUE, KINDX, J, 
     2                                                 Y(J), YTRUE(J)
                              END IF
                          END IF
C
                      ELSE
C
                          V = ABS ( X (LIST(J)) - XTRUE (LIST(J)) )
                          W = ABS ( CTRUE ) * ABS ( XSAVE (LIST(J)) )  +
     1                        ABS ( STRUE ) * ABS ( YSAVE(J) )
                          IF  ( W .EQ. 0.0E0 )  W = 1.0E0
                          ERR = V / ( EPSILN * W )
                          IF  ( ERR .GT. THRESH )  THEN
                              COUNT = COUNT + 1
                              IF  ( COUNT .LE. ERRMAX )  THEN 
                                  WRITE ( NOUT, 1600 ) NZTRUE, CTRUE, 
     1                                                 STRUE, KINDX, I,
     2                                                 X (LIST(J)), 
     3                                                 XTRUE (LIST(J)), 
     4                                                 ERR
                              END IF
                          END IF
C
                          V = ABS ( Y(J) - YTRUE(J) )
                          W = ABS ( STRUE ) * ABS ( XSAVE (LIST(J)) )  +
     1                        ABS ( CTRUE ) * ABS ( YSAVE(J) )
                          IF  ( W .EQ. 0.0E0 )  W = 1.0E0
                          ERR = V / ( EPSILN * W )
                          IF  ( ERR .GT. THRESH )  THEN
                              COUNT = COUNT + 1
                              IF  ( COUNT .LE. ERRMAX )  THEN 
                                  WRITE ( NOUT, 1700 ) NZTRUE, CTRUE, 
     1                                                 STRUE, KINDX, J,
     2                                                 Y(J), YTRUE(J), 
     3                                                 ERR
                              END IF
                          END IF
C
                      END IF
C
  400             CONTINUE
C
  500         CONTINUE
C
  600     CONTINUE
C
  700 CONTINUE
C
C     ==================================================================
C
C     ------------------
C     ... END OF TESTING
C     ------------------
C
      ERRCNT = ERRCNT + COUNT
      IF  ( COUNT .NE. 0 )  GO TO 800
C
C     -----------------------------------
C     ... WRITE PASSED MESSAGE AND RETURN
C     -----------------------------------
C
      WRITE ( NOUT, 2700 )
      GO TO 900
C
C     -----------------------------------
C     ... WRITE FAILED MESSAGE AND RETURN
C     -----------------------------------
C
  800 WRITE ( NOUT, 2800 ) COUNT
C
C     ------------------------
C     ... END OF MODULE TSROTI
C     ------------------------
C
  900 CONTINUE
      RETURN
C
C     ==================================================================
C
C     -----------
C     ... FORMATS
C     -----------
C
 1000 FORMAT ( 5X, 'SROTI ALTERED NZ FOR TEST WITH NZ = ', I5,
     1             ' C, S = ', 1P, 2E15.5, ' AND THE INDX TYPE NO. ', I5
     2        /5X, 'ALTERED VALUE OF NZ = ', I5 )
C
 1100 FORMAT ( 5X, 'SROTI ALTERED C FOR TEST WITH NZ = ', I5,
     1             ' C, S = ', 1P, 2E15.5, ' AND THE INDX TYPE NO. ', I5
     2        /5X, 'ALTERED VALUE OF C = ', 1PE15.5 )
C
 1200 FORMAT ( 5X, 'SROTI ALTERED S FOR TEST WITH NZ = ', I5,
     1             ' C, S = ', 1P, 2E15.5, ' AND THE INDX TYPE NO. ', I5
     2        /5X, 'ALTERED VALUE OF S = ', 1PE15.5 )
C
 1300 FORMAT ( 5X, 'SROTI ALTERED ARRAY INDX FOR TEST WITH NZ = ', I5,
     1             ' C, S = ', 1P, 2E15.5, ' AND THE INDX TYPE NO. ', 
     2             I5 )
C
 1400 FORMAT ( 5X, 'SROTI OUTPUT ARRAY X IS INCORRECT FOR TEST WITH ',
     1             'NZ = ', I5, ' C, S = ', 1P, 2E15.5,
     2             ' AND THE INDX TYPE NO. ', I5
     3        /5X, 'INCORRECT COMPONENT NO. ', I5, ' HAS VALUE = ', 
     4             1PE15.5, ' TRUE VALUE = ', 1PE15.5 )
C
 1500 FORMAT ( 5X, 'SROTI OUTPUT ARRAY Y IS INCORRECT FOR TEST WITH ',
     1             'NZ = ', I5, ' C, S = ', 1P, 2E15.5,
     2             ' AND THE INDX TYPE NO. ', I5
     3        /5X, 'INCORRECT COMPONENT NO. ', I5, ' HAS VALUE = ', 
     4             1PE15.5, ' TRUE VALUE = ', 1PE15.5 )
C
 1600 FORMAT ( 5X, 'SROTI OUTPUT ARRAY X IS INACCURATE FOR TEST WITH ',
     1             'NZ = ', I5, ' C, S = ', 1P, 2E15.5,
     2             ' AND THE INDX TYPE NO. ', I5
     3        /5X, 'INACCURATE COMPONENT NO. ', I5, ' HAS VALUE = ', 
     4             1PE15.5, ' TRUE VALUE = ', 1PE15.5, ' ERROR = ', 
     5             1PE12.1 )
C
 1700 FORMAT ( 5X, 'SROTI OUTPUT ARRAY Y IS INACCURATE FOR TEST WITH ',
     1             'NZ = ', I5, ' C, S = ', 1P, 2E15.5,
     2             ' AND THE INDX TYPE NO. ', I5
     3        /5X, 'INACCURATE COMPONENT NO. ', I5, ' HAS VALUE = ', 
     4             1PE15.5, ' TRUE VALUE = ', 1PE15.5, ' ERROR = ', 
     5             1PE12.1 )
C
 2700 FORMAT ( /5X, 'SROTI  PASSED ALL TESTS.' ) 
C
 2800 FORMAT ( /5X, 'SROTI  FAILED', I10, ' TESTS.'  )
C
C     ==================================================================
C
      END
      SUBROUTINE   TSSCTR   ( NOUT,   NZMAX2, NUMNZ,  NZVALU, 
     1                        X,      XSAVE,  XTRUE,  Y,      YSAVE, 
     2                        YTRUE , INDX,   INDXT,  ERRCNT, ERRMAX )
C
C     ==================================================================
C     ==================================================================
C     ====  TSSCTR  --  CERTIFY  SSCTR                              ====
C     ==================================================================
C     ==================================================================
C
C     SUBROUTINE  TSSCTR  IS THE CERTIFICATION MODULE FOR THE SPARSE
C     BASIC LINEAR ALGEBRA SUBROUTINE MODULE  SSCTR.
C
C     WRITTEN BY      ROGER G GRIMES
C                     APRIL 1987
C
C     ==================================================================
C
C     -------------
C     ... ARGUMENTS
C     -------------
C
      INTEGER             NOUT,   NZMAX2, NUMNZ,  ERRCNT,
     1                    ERRMAX
C
      INTEGER             NZVALU (*),  INDX (*),    INDXT (*)
C
      REAL                X (*),       XSAVE (*),   XTRUE (*),
     1                    Y (*),       YSAVE (*),   YTRUE (*)
C
C     -------------------
C     ... LOCAL VARIABLES
C     -------------------
C
      INTEGER             COUNT,  I,      ICLOBR, J,      KINDX,
     1                    KNZ,    N,      NZ,     NZTRUE
C
      REAL                CLOBBR   
C
C     --------------------
C     ... SUBPROGRAMS USED
C     --------------------
C
      LOGICAL             IVSAME, SVSAME
C
      EXTERNAL            ICOPY,  SCOPY,  SINIT,  GNINDX,
     1                    IVSAME, SVSAME, SSCTR
C
C     ==================================================================
C
C     ------------------
C     ... INITIALIZATION
C     ------------------
C
      COUNT     =   0
C
      CLOBBR    =   -1.0E10
      ICLOBR    =   -10000000
C
C     ------------------------------------
C     ... GENERATE SOME VALUES FOR X AND Y
C     ------------------------------------
C
      DO 100 I = 1, NZMAX2
         XSAVE(I) = COS ( .6*FLOAT(I) )
         YSAVE(I) = SIN ( .7*FLOAT(I) )
  100 CONTINUE
C
C     ------------------------
C     ... FOR EACH VALUE OF NZ
C     ------------------------
C
      DO 600 KNZ = 1, NUMNZ
C
          NZTRUE = NZVALU(KNZ)
          N      = 2 * MAX ( NZTRUE, 1 )
C
C         -------------------------------
C         ... FOR EACH KIND OF INDX ARRAY
C         -------------------------------
C
          DO 500 KINDX = 1, 5
C
              CALL GNINDX ( NZTRUE, N, ICLOBR, KINDX, INDXT )
C
C             -----------------------
C             ... GENERATE INPUT DATA
C             -----------------------
C
              I = MIN ( N, N-NZTRUE )
              J = N - I + 1
              CALL SCOPY ( NZTRUE, XSAVE,  1, XTRUE, 1 )
              CALL SINIT ( I,      CLOBBR, XTRUE(J), 1 )
              CALL SINIT ( N,      CLOBBR, YTRUE, 1 )
C
C             -------------------
C             ... COPY TRUE INPUT
C             -------------------
C
              NZ = NZTRUE
C
              CALL SCOPY ( N, YTRUE, 1, Y, 1 )
              CALL SCOPY ( N, XTRUE, 1, X, 1 )
              CALL ICOPY ( N, INDXT, 1, INDX, 1 )
C         
C             --------------------------
C             ... COMPUTE IN-LINE RESULT
C             --------------------------
C
              DO 300 I = 1, NZTRUE
                  YTRUE (INDXT(I)) = XTRUE (I)
  300         CONTINUE
C
C             --------------
C             ... CALL SSCTR
C             --------------
C
              CALL SSCTR ( NZ, X, INDX, Y )
C
C             ----------------------------------------
C             ... TEST ARGUMENTS OF SSCTR THAT ARE NOT
C                 SUPPOSE TO CHANGE.
C             ----------------------------------------
C
              IF  ( NZ .NE. NZTRUE )  THEN
                  COUNT = COUNT + 1
                  IF  ( COUNT .LE. ERRMAX )  THEN 
                      WRITE ( NOUT, 1000 ) NZTRUE, KINDX, NZ
                  END IF
              END IF
C                                          
              IF  ( .NOT. SVSAME ( N, X, XTRUE ) )  THEN
                  COUNT = COUNT + 1
                  IF  ( COUNT .LE. ERRMAX )  THEN 
                      WRITE ( NOUT, 1100 ) NZTRUE, KINDX
                  END IF
              END IF
C
              IF  ( .NOT. IVSAME ( N, INDX, INDXT ) )  THEN
                  COUNT = COUNT + 1
                  IF  ( COUNT .LE. ERRMAX )  THEN 
                      WRITE ( NOUT, 1200 ) NZTRUE, KINDX
                  END IF
              END IF
C
C             --------------------------
C             ... TEST OUTPUT FROM SSCTR
C             --------------------------
C
              DO 400 I = 1, N
                  IF  ( Y(I) .NE. YTRUE(I) )  THEN
                      COUNT = COUNT + 1
                      IF  ( COUNT .LE. ERRMAX )  THEN 
                          WRITE ( NOUT, 1300 ) NZTRUE, KINDX, I, 
     1                                         Y(I), YTRUE(I)
                      END IF
                  END IF
  400         CONTINUE
C
  500     CONTINUE
C     
  600 CONTINUE
C
C     ==================================================================
C
C     ------------------
C     ... END OF TESTING
C     ------------------
C
      ERRCNT = ERRCNT + COUNT
      IF  ( COUNT .NE. 0 )  GO TO 800
C
C     -----------------------------------
C     ... WRITE PASSED MESSAGE AND RETURN
C     -----------------------------------
C
      WRITE ( NOUT, 2700 )
      GO TO 900
C
C     -----------------------------------
C     ... WRITE FAILED MESSAGE AND RETURN
C     -----------------------------------
C
  800 WRITE ( NOUT, 2800 ) COUNT
C
C     ------------------------
C     ... END OF MODULE TSSCTR
C     ------------------------
C
  900 CONTINUE
      RETURN
C
C     ==================================================================
C
C     -----------
C     ... FORMATS
C     -----------
C
 1000 FORMAT ( 5X, 'SSCTR ALTERED NZ FOR TEST WITH NZ = ', I5,
     1             ' AND THE INDX TYPE NO. ', I5,
     2             '.  ALTERED VALUE OF NZ = ', I5 )
C
 1100 FORMAT ( 5X, 'SSCTR ALTERED ARRAY X FOR TEST WITH NZ = ', I5,
     1             ' AND THE INDX TYPE NO. ', I5 )
C
 1200 FORMAT ( 5X, 'SSCTR ALTERED ARRAY INDX FOR TEST WITH NZ = ', I5,
     1             ' AND THE INDX TYPE NO. ', I5 )
C
 1300 FORMAT ( 5X, 'SSCTR OUTPUT ARRAY Y IS INCORRECT FOR TEST WITH ',
     1             'NZ = ', I5, ' AND THE INDX TYPE NO. ', I5
     2        /5X, 'INACCURATE COMPONENT NO. ', I5, ' HAS VALUE =', 
     3             1PE15.5, ' TRUE VALUE =', 1PE15.5 )
C
 2700 FORMAT ( /5X, 'SSCTR  PASSED ALL TESTS.' ) 
C
 2800 FORMAT ( /5X, 'SSCTR  FAILED', I10, ' TESTS.'  )
C
C     ==================================================================
C
      END
