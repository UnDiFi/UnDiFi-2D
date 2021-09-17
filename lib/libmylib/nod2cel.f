      SUBROUTINE NOD2CEL(NPOIN,NELEM,NOFVERT,ICELNOD,JA,LENJA,IA,LENIA,
     +                   IERR)
C
      IMPLICIT NONE
C
C     This routine returns the linKed list of the elements each node
C     is surrounded by
C     IA(K) points to the position of JA() where the list of node K starts
C
C
C     .. Scalar Arguments ..
      INTEGER IERR,LENIA,LENJA,NELEM,NOFVERT,NPOIN
C     ..
C     .. Array Arguments ..
      INTEGER IA(*),ICELNOD(NOFVERT,*),JA(*)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION A
      INTEGER IELEM,II,J,JBGN,JEND,K,KNOD,KSAV,KSAVN
C     ..
C     .. External Subroutines ..
      EXTERNAL IINIT
C     ..
      IF (LENIA.LT. (NPOIN+1)) THEN
          IERR = 1
          RETURN

      ENDIF
C
      CALL IINIT(NPOIN+1,0,IA(1),1)
C
C     Count the number of elements surrounding each node
C
      DO 3 IELEM = 1,NELEM
          DO 4 J = 1,NOFVERT
              KNOD = ICELNOD(J,IELEM)
              IF((KNOD.LT.1).OR.(KNOD.GT.NPOIN))THEN
                  IERR = 3
                  RETURN
              ENDIF
              IA(KNOD) = IA(KNOD) + 1
    4     CONTINUE
    3 CONTINUE
C
C     ... construct the connectivity list ...
C
      KSAV = IA(1)
      IA(1) = 1
      DO 6 J = 2,NPOIN + 1
          KSAVN = IA(J)
          IA(J) = IA(J-1) + KSAV
          KSAV = KSAVN
    6 CONTINUE
C
C     DO 8 K = 1, IA(NPOIN+1)-1
C        JA(K) = 0
C   8 CONTINUE
      IF (LENJA.LT. (IA(NPOIN+1)-1)) THEN
          LENJA = IA(NPOIN+1)-1
          IERR = 2
          RETURN

      ENDIF

      CALL IINIT(IA(NPOIN+1)-1,0,JA(1),1)
C
C     Main loop
C
      DO 12 IELEM = 1,NELEM
c
c     get nodal points
c
          DO 14 J = 1,NOFVERT
c
              KNOD = ICELNOD(J,IELEM)
c
              JBGN = IA(KNOD)
              JEND = IA(KNOD+1) - 1
              DO 13 K = JBGN,JEND
                  IF (JA(K).NE.0) GOTO 13
                  JA(K) = IELEM
                  GOTO 14

   13         CONTINUE
   14     CONTINUE
   12 CONTINUE
C
C
      LENJA = IA(NPOIN+1)-IA(1)
C
C
      IERR = 0
      RETURN

      END
