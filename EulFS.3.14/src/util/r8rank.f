      SUBROUTINE R8RANK(RV,M1,M2,IRANK,IFAIL)
      IMPLICIT NONE
C
C     R8RANK RE-ARRANGES A VECTOR OF REAL NUMBERS INTO
C     THE ORDER SPECIFIED BY A VECTOR OF RANKS.
C
C     .. Parameters ..
      CHARACTER*6       SUBNAME
      PARAMETER         (SUBNAME='R8RANK')
      INTEGER           IFAIL, M1, M2
      DOUBLE PRECISION  RV(M2)
      INTEGER           IRANK(M2)
      DOUBLE PRECISION  A, B
      INTEGER           I, IERR, J, K
      CHARACTER*80      ERRMSG(2)
      INTRINSIC         ABS
C
C
      IF (M2.LT.1 .OR. M1.LT.1 .OR. M1.GT.M2) THEN
         IERR = 1
         WRITE (ERRMSG,FMT=99999) M1, M2
      ELSE
         IERR = 0
         DO 20 I = M1, M2
            J = IRANK(I)
            IF ((J.LT.M1) .OR. (J.GT.M2)) GO TO 100
            IF (I.NE.J) IRANK(I) = -J
   20    CONTINUE
C
C
         DO 60 I = M1, M2
            K = -IRANK(I)
            IF (K.GE.0) THEN
               J = I
               A = RV(I)
   40          IRANK(J) = K
               B = RV(K)
               RV(K) = A
               J = K
               A = B
               K = -IRANK(J)
               IF (K.GT.0) GO TO 40
               IF (J.NE.I) GO TO 120
            END IF
   60    CONTINUE
      END IF
C
C       RETURN
C
   80 IF (IERR.NE.0) THEN
         IF (IERR.EQ.1)THEN
            CALL SETERR(SUBNAME//ERRMSG(1),72,IERR,1)
         ELSE
            CALL SETERR(SUBNAME//ERRMSG(2),72,IERR,1)
         ENDIF
         IFAIL = IERR
      ELSE
         IFAIL = 0
      END IF
      RETURN
  100 IERR = 2
      WRITE (ERRMSG(2),FMT=99997) I, J
      GO TO 140
  120 IERR = 3
      WRITE (ERRMSG(2),FMT=99996) J
  140 WRITE (ERRMSG(1),FMT=99998)
C
C     RESTORE IRANK
C
      DO 160 J = M1, M2
         IRANK(J) = ABS(IRANK(J))
  160 CONTINUE
      GO TO 80
C
99999 FORMAT (' ** On entry: M1 =',I16,'  M2 =',I16)
99998 FORMAT (' ** IRANK(M1:M2) does not contain a permutation of the ',
     *  'integers M1 to M2')
99997 FORMAT (' IRANK(',I6,') contains an out-of-range value',I16)
99996 FORMAT (' IRANK contains a repeated value',I16)
      END
