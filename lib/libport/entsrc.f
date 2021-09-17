      SUBROUTINE ENTSRC(IROLD,IRNEW)
C
C  THIS ROUTINE RETURNS IROLD = LRECOV AND SETS LRECOV = IRNEW.
C
C  IF THERE IS AN ACTIVE ERROR STATE, THE MESSAGE IS PRINTED
C  AND EXECUTION STOPS.
C
C  IRNEW = 0 LEAVES LRECOV UNCHANGED, WHILE
C  IRNEW = 1 GIVES RECOVERY AND
C  IRNEW = 2 TURNS RECOVERY OFF.
C
C  ERROR STATES -
C
C    1 - ILLEGAL VALUE OF IRNEW.
C    2 - CALLED WHILE IN AN ERROR STATE.
C
      IF (IRNEW.LT.0 .OR. IRNEW.GT.2)
     1   CALL SETERR(31HENTSRC - ILLEGAL VALUE OF IRNEW,31,1,2)
C
      IROLD=I8SAVE(2,IRNEW,IRNEW.NE.0)
C
C  IF HAVE AN ERROR STATE, STOP EXECUTION.
C
      IF (I8SAVE(1,0,.FALSE.) .NE. 0) CALL SETERR
     1   (39HENTSRC - CALLED WHILE IN AN ERROR STATE,39,2,2)
C
      RETURN
C
      END
