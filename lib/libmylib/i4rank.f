      SUBROUTINE I4RANK(IV,M1,M2,IRANK,IFAIL)
      IMPLICIT NONE
C
C     I4RANK RE-ARRANGES A VECTOR OF INTEGER NUMBERS INTO
C     THE ORDER SPECIFIED BY A VECTOR OF RANKS.
C
C
      CHARACTER*6       SUBNAME
      PARAMETER         (SUBNAME='I4RANK')
      INTEGER           IFAIL, M1, M2
      INTEGER           IRANK(M2), IV(M2)
      INTEGER           A, B, I, IERR, J, K
      CHARACTER*80      ERRSTR(2)
      INTRINSIC         ABS
C
C
      IF (M2.LT.1 .OR. M1.LT.1 .OR. M1.GT.M2) THEN
         IERR = 1
         WRITE (ERRSTR,FMT=999) M1, M2
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
               A = IV(I)
   40          IRANK(J) = K
               B = IV(K)
               IV(K) = A
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
            CALL SETERR(SUBNAME//ERRSTR(1),72,IERR,1)
         ELSE
            CALL SETERR(SUBNAME//ERRSTR(2),72,IERR,1)
         ENDIF
         IFAIL = IERR
      ELSE
         IFAIL = 0
      END IF
      RETURN
  100 IERR = 2
      WRITE (ERRSTR(2),FMT=997) I, J
      GO TO 140
  120 IERR = 3
      WRITE (ERRSTR(2),FMT=996) J
  140 WRITE (ERRSTR(1),FMT=998)
C
C     RESTORE IRANK
C
      DO 160 J = M1, M2
         IRANK(J) = ABS(IRANK(J))
  160 CONTINUE
      GO TO 80
C
  999 FORMAT (' ** On entry: M1 =',I16,'  M2 =',I16)
  998 FORMAT (' ** IRANK(M1:M2) does not contain a permutation of the ',
     *  'integers M1 to M2')
  997 FORMAT (' IRANK(',I6,') contains an out-of-range value',I16)
  996 FORMAT (' IRANK contains a repeated value',I16)
      END
