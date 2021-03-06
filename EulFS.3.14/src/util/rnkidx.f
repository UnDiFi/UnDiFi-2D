      SUBROUTINE RNKIDX(IPERM,M1,M2,IFAIL)
      IMPLICIT NONE
C
C     RNKIDX INVERTS A PERMUTATION, AND HENCE CONVERTS A RANK VECTOR TO
C     AN INDEX VECTOR, OR VICE VERSA.
C
      CHARACTER*6       SUBNAME
      PARAMETER         (SUBNAME='RNKIDX')
      INTEGER           IFAIL, M1, M2
      INTEGER           IPERM(M2)
      INTEGER           I, IERR, J, K, L
      CHARACTER*80      ERRMSG(2)
C
C
      IF (M2.LT.1 .OR. M1.LT.1 .OR. M1.GT.M2) THEN
         IERR = 1
         WRITE (ERRMSG,FMT=99999) M1, M2
      ELSE
         IERR = 0
         DO 20 I = M1, M2
            J = IPERM(I)
            IF ((J.LT.M1) .OR. (J.GT.M2)) GO TO 100
            IF (I.NE.J) IPERM(I) = -J
   20    CONTINUE
C
C
         DO 60 I = M1, M2
            K = -IPERM(I)
            IF (K.GE.0) THEN
               J = I
   40          L = -IPERM(K)
               IF (L.GT.0) THEN
                  IPERM(K) = J
                  J = K
                  K = L
                  GO TO 40
               END IF
               IF (IPERM(I).LT.0) GO TO 120
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
      WRITE (ERRMSG(2),FMT=99996) K
  140 WRITE (ERRMSG(1),FMT=99998)
      GO TO 80
C
99999 FORMAT (' ** On entry: M1 =',I16,'  M2 =',I16)
99998 FORMAT (' ** IPERM(M1:M2) does not contain a permutation of the ',
     *  'integers M1 to M2')
99997 FORMAT (' IPERM(',I6,') contains an out-of-range value',I16)
99996 FORMAT (' IPERM contains a repeated value',I16)
      END
