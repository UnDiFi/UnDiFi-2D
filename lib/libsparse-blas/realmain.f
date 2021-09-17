      PROGRAM TSSPBL
C
C
C     ==================================================================
C     ==================================================================
C     ====  TSSPBL -- CERTIFY SINGLE PRECISION SPARSE BLAS          ====
C     ==================================================================
C     ==================================================================
C
C     TSSPBL IS THE CERTIFICATION PROGRAM FOR THE SINGLE PRECISION
C     SPARSE BLAS.  THE APPROACH USED TO CERTIFY THE SPARSE BLAS
C     IS AS FOLLOWS:
C
C     1.  READ IN USER SPECIFIED INPUT ON OUTPUT UNIT, THRESHOLD VALUE
C         FOR TEST RATIO, AND THE SPECIFICATIONS FOR NZ, A, C AND S.
C     2.  VERIFY THE CORRECTNESS OF THE USER SPECIFIED INPUT AND
C         ECHO TO THE OUTPUT UNIT.
C     3.  FOR EACH SUBPROGRAM IN THE REAL PRECISION SPARSE BLAS
C         PERFORM ALL THE USER SPECIFIED TESTS AND PRINT A PASS/FAIL
C         MESSAGE.  TESTS WHICH FAIL GENERATE ADDITIONAL OUTPUT.
C
C     SPARSE BLAS SUBPROGRAMS WHICH ARE CERTIFIED BY THIS PROGRAM ARE
C
C         SAXPYI          SGTHR           SROTI
C         SDOTI           SGTHRZ          SSCTR
C
C     THIS PROGRAM REQUIRES AN INPUT FILE ASSIGNED TO UNIT NIN 
C     (CURRENTLY SET TO 5 BY A PARAMETER STATEMENT).  THE DATA ON
C     THIS INPUT FILE CONTROLS THE OUTPUT UNIT, THE THRESHOLD VALUE
C     FOR THE NUMERICAL TESTING, AND THE SPECIFICATIONS FOR THE
C     TEST VALUES FOR THE LENGTH OF THE SPARSE VECTORS AND THE SCALARS 
C     USED BY THE VARIOUS SUBPROGRAMS.  AN EXAMPLE OF THE INPUT FILE 
C     FOLLOWS
C
C LINE  1     'SBLATS.SUMM'           NAME OF OUTPUT FILE
C LINE  2     6                       UNIT NUMBER OF OUTPUT FILE
C LINE  3     100                     MAX. NO. OF PRINTED ERROR MESSAGES
C LINE  4     5.0                     THRESHOLD VALUE OF TEST RATIO
C LINE  5     6                       NUMBER OF VALUES OF NZ
C LINE  6     -1 0 1 2 5 9            VALUES OF NZ
C LINE  7     3                       NUMBER OF VALUES OF A FOR -AXPYI
C LINE  8     0.0 1.0 0.7             VALUES OF A
C LINE  9     5                       NUMBER OF VALUES OF C,S FOR -ROTI
C LINE 10     0. 1. 0. -.6  .8        VALUES OF C
C LINE 11     0. 0. 1   .8 -.6        VALUES OF S
C
C
C     THIS INPUT FILE IS READ USING FORTRAN-77 STANDARD LIST DIRECTED
C     INPUT.  SINGLE QUOTES ARE REQUIRED AROUND THE NAME OF THE OUTPUT
C     FILE ON LINE 1.  THE NUMBERS ON LINES 6, 8, 10, AND 11 CAN BE
C     DELIMITED BY BLANKS OR COMMAS.
C
C     THIS PROGRAM WAS WRITTEN BY ROGER G. GRIMES, BOEING
C     COMPUTER SERVICES, BELLEVUE, WA. DURING APRIL, 1987.
C
C     ==================================================================
C
C     ------------------------------------
C     ... PROBLEM SPECIFICATION PARAMETERS
C     ------------------------------------
C
C         NIN         INPUT UNIT
C         NZMAX       MAXIMUM VALUE OF ANY SINGLE NZ
C         NNZMAX      MAXIMUM NUMBER OF VALUES OF NZ 
C         NAMAX       MAXIMUM NUMBER OF VALUES OF A (-AXPYI
C                     SCALAR)
C         NGMAX       MAXIMUM NUMBER OF VALUES OF C AND S
C                     (-ROTI SCALARS FOR GIVENS ROTATION)
C
C     ==================================================================
C
      INTEGER             NIN,    NZMAX,  NNZMAX, NAMAX,  NGMAX
C
      PARAMETER         ( NIN = 5,        NZMAX = 320,
     1                    NNZMAX = 24,    NAMAX = 7,      NGMAX = 7 )
C
C     -----------------------
C     ... COMPUTED PARAMETERS
C     -----------------------
C
      INTEGER             NZMAX2
C
      PARAMETER         ( NZMAX2 = 2 * NZMAX )
C
C     ==================================================================
C
C     ------------------------
C     ... VARIABLE DECLARATION
C     ------------------------
C
      CHARACTER*32        NAMOUT
C
      INTEGER             ERRCNT, ERRMAX, I,      NOUT,   NUMA,   
     1                    NUMG,   NUMNZ
C
      INTEGER             INDX  (NZMAX2),     INDXT (NZMAX2),
     1                    LIST  (NZMAX2),     NZVALU(NNZMAX)
C
      REAL                EPSILN, EPSSAV, THRESH
C
      REAL                X     (NZMAX2),     Y     (NZMAX2),
     1                    XTRUE (NZMAX2),     YTRUE (NZMAX2),
     2                    XSAVE (NZMAX2),     YSAVE (NZMAX2),
     3                    AVALUE(NAMAX),      CVALUE(NGMAX),
     4                    SVALUE(NGMAX)
C
C     --------------------
C     ... SUBPROGRAMS USED
C     --------------------
C
      REAL                SDIFF
C
      EXTERNAL            TSXPYI, TSDOTI, TSGTHR, TSGTHZ, TSROTI,
     1                    TSSCTR, SDIFF
C
C     ==================================================================
C
      ERRCNT = 0
C
C     ------------------------------------------------
C     ... READ IN USER SPECIFIED INPUT FOR OUTPUT UNIT
C     ------------------------------------------------
C
      READ ( NIN, * ) NAMOUT
      READ ( NIN, * ) NOUT
C
C     --------------------
C     ... OPEN OUTPUT UNIT
C     --------------------
C
      OPEN ( UNIT = NOUT, FILE = NAMOUT, STATUS = 'NEW' )
C
C     ------------------------------
C     ... READ IN REMAINDER OF INPUT
C     ------------------------------
C
      READ ( NIN, * ) ERRMAX
      READ ( NIN, * ) THRESH
      READ ( NIN, * ) NUMNZ
C
      IF ( NUMNZ .GT. NNZMAX ) THEN
          ERRCNT = 1
          WRITE ( NOUT, 1100 ) NUMNZ, NNZMAX
          GO TO 900
      END IF
C
      READ ( NIN, * ) ( NZVALU(I), I = 1, NUMNZ )
C
      READ ( NIN, * ) NUMA
C
      IF ( NUMA .GT. NAMAX ) THEN
          ERRCNT = 1
          WRITE ( NOUT, 1110 ) NUMA, NAMAX
          GO TO 900
      END IF
C
      READ ( NIN, * ) ( AVALUE(I), I = 1, NUMA  )
C
      READ ( NIN, * ) NUMG
C
      IF ( NUMG .GT. NGMAX ) THEN
          ERRCNT = 1
          WRITE ( NOUT, 1120 ) NUMG, NGMAX
          GO TO 900
      END IF
C
      READ ( NIN, * ) ( CVALUE(I), I = 1, NUMG  )
      READ ( NIN, * ) ( SVALUE(I), I = 1, NUMG  )
C
C     ------------------------------
C     ... PRINT USER SPECIFIED INPUT
C     ------------------------------
C
      WRITE ( NOUT, 1000 ) NAMOUT, NOUT, ERRMAX, THRESH
      WRITE ( NOUT, 1010 ) NUMNZ
      WRITE ( NOUT, 1020 ) ( NZVALU(I), I = 1, NUMNZ )
      WRITE ( NOUT, 1030 ) NUMA
      WRITE ( NOUT, 1040 ) ( AVALUE(I), I = 1, NUMA  )
      WRITE ( NOUT, 1050 ) NUMG
      WRITE ( NOUT, 1060 ) ( CVALUE(I), I = 1, NUMG  )
      WRITE ( NOUT, 1070 ) ( SVALUE(I), I = 1, NUMG  )
C
C     -------------------------------
C     ... VERIFY USER SPECIFIED INPUT
C     -------------------------------
C
      IF  ( THRESH .LE. 0.0E0 )  THEN
          WRITE ( NOUT, 1130 ) THRESH
          THRESH = 10.0E0
      END IF
C
      IF  ( NUMNZ .LE. 0 )  THEN
          WRITE ( NOUT, 1140 ) NUMNZ
          ERRCNT = 1
      END IF
C
      DO 100 I = 1, NUMNZ
          IF  ( NZVALU(I) .GT. NZMAX )  THEN
              WRITE ( NOUT, 1150 ) I, NZVALU(I), NZMAX
              NZVALU(I) = NZMAX
          END IF
  100 CONTINUE
C
      IF  ( ERRCNT .NE. 0 )  GO TO 900
C
C     ---------------------------
C     ... COMPUTE MACHINE EPSILON
C     ---------------------------
C
      EPSILN = 1.0E0
      EPSSAV = 1.0E0
C
  200 IF  ( SDIFF ( 1.0E0 + EPSILN, 1.0E0 ) .EQ. 0.0E0 )  GO TO 210
C
          EPSSAV = EPSILN
          EPSILN = EPSILN * .5E0
          GO TO 200
C
  210 EPSILN = EPSSAV
C
C     ==================================================================
C
C     -----------------------------
C     ... TEST THE REAL SPARSE BLAS
C     -----------------------------
C
C     ------------------
C     ... CERTIFY SAXPYI
C     ------------------
C
      CALL TSXPYI (   NOUT,   EPSILN, THRESH, NZMAX2,
     1                NUMNZ,  NZVALU, NUMA,   AVALUE,
     2                X,      XSAVE,  XTRUE,  Y,      YSAVE,  YTRUE,
     3                INDX,   INDXT,  LIST,   ERRCNT, ERRMAX )
C
C     -----------------
C     ... CERTIFY SDOTI
C     -----------------
C
      CALL TSDOTI (   NOUT,   EPSILN, THRESH, NZMAX2,
     1                NUMNZ,  NZVALU, 
     2                X,      XSAVE,  XTRUE,  Y,      YSAVE,  YTRUE,
     3                INDX,   INDXT,  ERRCNT, ERRMAX )
C
C     -----------------
C     ... CERTIFY SGTHR
C     -----------------
C
      CALL TSGTHR (   NOUT,   NZMAX2, NUMNZ,  NZVALU, 
     1                X,      XSAVE,  XTRUE,  Y,      YSAVE,  YTRUE,
     2                INDX,   INDXT,  ERRCNT, ERRMAX )
C
C     ------------------
C     ... CERTIFY SGTHRZ
C     ------------------
C
      CALL TSGTHZ (   NOUT,   NZMAX2, NUMNZ,  NZVALU, 
     1                X,      XSAVE,  XTRUE,  Y,      YSAVE,  YTRUE,
     2                INDX,   INDXT,  ERRCNT, ERRMAX )
C
C     -----------------
C     ... CERTIFY SROTI
C     -----------------
C
      CALL TSROTI (   NOUT,   EPSILN, THRESH, NZMAX2,
     1                NUMNZ,  NZVALU, NUMG,   CVALUE, SVALUE,
     2                X,      XSAVE,  XTRUE,  Y,      YSAVE,  YTRUE,
     3                INDX,   INDXT,  LIST,   ERRCNT, ERRMAX )
C
C     -----------------
C     ... CERTIFY SSCTR
C     -----------------
C
      CALL TSSCTR (   NOUT,   NZMAX2, NUMNZ,  NZVALU, 
     1                X,      XSAVE,  XTRUE,  Y,      YSAVE,  YTRUE,
     2                INDX,   INDXT,  ERRCNT, ERRMAX )
C
C     ==================================================================
C
C     -------------------------------------
C     ... PRINT GLOBAL PASS OR FAIL MESSAGE
C     -------------------------------------
C
  900 IF  ( ERRCNT .EQ. 0 )  THEN
          WRITE ( NOUT, 2000 )
      ELSE
          WRITE ( NOUT, 2100 ) ERRCNT
      END IF
C
C     -----------------------------------------
C     ... END OF CERTIFICATION PROGRAM FOR REAL
C         SPARSE BLAS
C     -----------------------------------------
C
      STOP
C
C     ==================================================================
C
C     -----------
C     ... FORMATS
C     -----------
C
 1000 FORMAT( '1' ///
     1          5X, 'START OF CERTIFICATION PROGRAM FOR THE REAL ',
     2              'PRECISION SPARSE BLAS'
     3         /5X, '--------------------------------------------',
     4              '---------------------'
     5        //5X, 'NAME   OF OUTPUT UNIT              = ', A
     6         /5X, 'NUMBER OF OUTPUT UNIT              = ', I10
     7         /5X, 'MAX. NO. OF PRINTED ERROR MESSAGES = ', I10    
     8         /5X, 'THRESHOLD VALUE OF TEST RATIO      = ', F10.1 )
C 
 1010 FORMAT ( /5X, 'NUMBER OF VALUES OF NZ        = ', I10 )
C
 1020 FORMAT ( /5X, 'VALUES OF NZ = ', 10I5 )
C
 1030 FORMAT ( /5X, 'NUMBER OF VALUES OF A         = ', I10 )
C
 1040 FORMAT ( /5X, 'VALUES OF A = ', 1P, 5E13.4 )
C
 1050 FORMAT ( /5X, 'NUMBER OF VALUES OF C AND S   = ', I10 )
C
 1060 FORMAT ( /5X, 'VALUES OF C = ', 1P, 5E13.4 )
C
 1070 FORMAT ( /5X, 'VALUES OF S = ', 1P, 5E13.4 )
C
 1100 FORMAT ( /5X, 'USER SPECIFIED NUMBER OF TEST CASES FOR THE ',
     1              'NUMBER OF NONZEROES EXCEEDS PROGRAM LIMIT.'
     2         /5X, 'NUMBER SPECIFIED = ', I10, 2X, 'PROGRAM LIMIT =',
     3              I10 )
C
 1110 FORMAT ( /5X, 'USER SPECIFIED NUMBER OF TEST CASES FOR THE ',
     1              'SCALAR A EXCEEDS PROGRAM LIMIT.'
     2         /5X, 'NUMBER SPECIFIED = ', I10, 2X, 'PROGRAM LIMIT =',
     3              I10 )
C
 1120 FORMAT ( /5X, 'USER SPECIFIED NUMBER OF TEST CASES FOR THE ',
     1              'SCALARS C AND S EXCEEDS PROGRAM LIMIT.'
     2         /5X, 'NUMBER SPECIFIED = ', I10, 2X, 'PROGRAM LIMIT =',
     3              I10 )
C
 1130 FORMAT ( /5X, 'USER SPECIFIED VALUE FOR THRESHOLD IS ', 1PE15.5,
     1              ' WHICH IS NONPOSITIVE.  IT HAS BEEN RESET TO 10.')
C
 1140 FORMAT ( /5X, 'USER SPECIFIED NUMBER OF VALUES OF NZ IS ', I5,
     1              ' WHICH IS NONPOSITIVE.  NO TESTING WILL OCCUR.' )
C
 1150 FORMAT ( /5X, 'THE ', I3, '-TH USER SPECIFIED VALUE OF NZ IS ', 
     1              I8, ' IS LARGER THAN THE MAXIMUM ALLOWABLE ',
     2              'VALUE OF NZ.  IT HAS BEEN RESET TO ', I5 )
C
 2000 FORMAT ( /5X, 'REAL PRECISION SPARSE BLAS HAVE PASSED ALL ',
     1              'TESTS.' )
C
 2100 FORMAT ( /5X, 'REAL SPARSE BLAS HAVE FAILED ', I10,
     1              ' TESTS.  SEE ABOVE PRINTED ERROR MESSAGES.' )
C
C     ==================================================================
C
      END
      REAL FUNCTION SDIFF ( X, Y )
C
C     ==================================================================
C
C     SDIFF IS USED BY THE MAIN PROGRAM TO COMPARE 1.0 + EPSILN WITH
C     1.0.  ITS SOLE USE IS TO FOOL AN OPTIMIZING COMPILER.
C
C     ==================================================================
C
C     ------------------------
C     ... VARIABLE DECLARATION
C     ------------------------
C
      REAL                X, Y
C
C     ==================================================================
C
      SDIFF = X - Y
C
C     ==================================================================
C
      RETURN
      END
      LOGICAL FUNCTION   SVSAME   ( N, SX, SY )
C
C     ==================================================================
C
C     LOGICAL FUNCTION  SVSAME  DETERMINES IF THE VECTORS  SX  AND  SY
C     AGREE EXACTLY WITH EACH OTHER.
C
C     ==================================================================
C
C     ------------------------
C     ... VARIABLE DECLARATION
C     ------------------------
C
      INTEGER             I, N 
C
      REAL                SX (*), SY (*)
C
C     ==================================================================
C
      SVSAME = .TRUE.
C
      DO 10 I = 1, N
          IF  ( SX(I) .NE. SY(I) )  THEN
              SVSAME = .FALSE.
              GO TO 20
          ENDIF
   10 CONTINUE
C
   20 RETURN
      END  
