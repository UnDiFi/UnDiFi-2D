      SUBROUTINE SORTN2C(NPOIN,NELEM,NOFVERT,ICELCEL,JA,LENJA,IA,LENIA,
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
      INTEGER IA(*),ICELCEL(NOFVERT,*),JA(*)
C     ..
C     .. Local Scalars ..
      INTEGER IPOIN,NN,J,JBGN,JEND,K,JELEM,IPOS,ILAST,IELEM,NITEMS,IFAIL
C     ..
C     .. Local Arrays ..
      INTEGER IWKL
      PARAMETER (IWKL=10)
      INTEGER IWKSP(IWKL)
C     ..
C     .. External Subroutines ..
      INTEGER BISRCH
      EXTERNAL IINIT
C     ..
      write(6,*)'This is where I sort'
      IF (LENIA.LT. (NPOIN+1)) THEN
          IERR = 1
          RETURN

      ENDIF
C
C     Count the number of elements surrounding each node
C
      DO 3 IPOIN = 1,NPOIN
          JBGN = IA(IPOIN)
          JEND = IA(IPOIN+1)-1
          NITEMS = JEND-JBGN+1
          CALL M01CBF(JA,JBGN,JEND,"A",IFAIL)
          IF(NITEMS.GT.IWKL)stop 'IWKSP too short in sortn2c'
          CALL IINIT(NITEMS,0,IWKSP,1) 
c         loop over triangles neighbouring a given gridpoint
c         pick up the first cell among the boundary cells,
c         otherwise (if the node is NOT on the bounadry) will pick up the last one
          DO J = JBGN,JEND
             IELEM = JA(J)
             DO K = 1,NOFVERT
                IF( (ICELCEL(K,IELEM) .EQ. 0) .OR.
     &              (ICELCEL(K,IELEM) .GT. NELEM) )GOTO 9
             ENDDO
          ENDDO 
    9     CONTINUE
          NN = 1
          IWKSP(NN) = IELEM
    8     DO 4 J = 1,NOFVERT
               JELEM = ICELCEL(J,IELEM)
C         try to locate element JELEM
               CALL BINSRC(JELEM,JA(JBGN),NITEMS,IPOS,ILAST)
               IF(IPOS.EQ.0)THEN
               ELSE
C         see if it is already in the list
                  DO 7 K = 1,NN
                     IF(IWKSP(K).EQ.JELEM)THEN
                        GOTO 4
                     ENDIF 
    7             CONTINUE
C              JELEM is NOT already in the list: add 
                  NN = NN+1
                  IWKSP(NN) = JELEM
                  IELEM = JELEM
                  GOTO 8
               ENDIF
    4     CONTINUE
          IF(NN.NE.NITEMS)THEN
             WRITE(6,*)(JA(J),J=JBGN,JEND)
             WRITE(6,*)(IWKSP(J),J=1,NITEMS)
             DO J = JBGN,JEND
                WRITE(6,*)JA(J),(ICELCEL(K,JA(J)),K=1,NOFVERT)
             ENDDO 
             STOP
          ELSE 
!            WRITE(6,*)(JA(J),J=JBGN,JEND)
!            WRITE(6,*)(IWKSP(J),J=1,NITEMS)
             K = 0
             DO J = JBGN,JEND
                K = K+1
                JA(J) = IWKSP(K)
             ENDDO
          ENDIF
    3 CONTINUE
C
      IERR = 0
      RETURN

      END
