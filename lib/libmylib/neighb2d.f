      SUBROUTINE NEIGHB2D(ICELNOD,ICELCEL,NELEM,DEGREE,IA,JA,LENJA,
     +NPOIN)

      IMPLICIT NONE

      INTEGER ICELNOD(3,*),
     +        ICELCEL(3,*),DEGREE(*),JA(*),IA(*)

      INTEGER NELEM,LENJA,NPOIN,NBELEM,MINDEG,MAXDEG
      INTEGER IELEM,INODE,I,II,JJ,KK,LL,ISUM
      INTEGER  ICYCL,ISAME2D
      EXTERNAL ICYCL,ISAME2D
C
      DO I = 1 , NPOIN
       DEGREE(I) = 0
      ENDDO
C
C      counts the number of elements surrounding a given node
C      and store them in DEGREE
C
      DO 100 IELEM = 1 , NELEM
         DO 100 I = 1 , 3
            INODE = ICELNOD(I,IELEM)
            DEGREE(INODE) = DEGREE(INODE) + 1
  100 CONTINUE

      write(6,*)npoin
      ISUM = 0      
      MINDEG = 10000
      MAXDEG = 0
      DO 90 I = 1 , NPOIN
       ISUM = ISUM + DEGREE(I)
         MINDEG = MIN(MINDEG,DEGREE(I))
         MAXDEG = MAX(MAXDEG,DEGREE(I))
       if(DEGREE(I).LE.0)WRITE(6,*)'Non positive degree ',DEGREE(I),
     + ' vertex ',i
   90 CONTINUE
      WRITE(6,*)'Min/Avg/Max Vertex degree ',MINDEG,ISUM/NPOIN,MAXDEG
      write(16,*)(DEGREE(I),i=1,NPOIN)

      IF( ISUM .GT. LENJA )THEN
       WRITE(6,*)'ISUM,LENJA are : ',ISUM,LENJA
       STOP
      ENDIF
C
C     IA(I) points to the position in JA where
C     data of node I begin
C
      IA(1) = 1
      DO 80 I = 2 , NPOIN
       IA(I) = IA(I-1) + DEGREE(I-1)
* write(16,*)i,ia(i)
   80 CONTINUE
C
      DO 85 I = 1 , NPOIN
       DEGREE(I) = 0
   85 CONTINUE

      DO 70 IELEM = 1 , NELEM
         DO 70 I = 1 , 3
            INODE = ICELNOD(I,IELEM)
            JA(IA(INODE)+DEGREE(INODE)) = IELEM
            DEGREE(INODE) = DEGREE(INODE) + 1
   70 CONTINUE
*      write(18,*)(DEGREE(I),i=1,NPOIN)
C
C     The neighbouring element opposite the i-th vertex
C     in the IELEM-th cell is given by the intersection
C     among the three lists of neighbouring elements
C     of the three nodes forming the face opposite node ICELNOD(i,IELEM)
C     The intersection will give the element IELEM and the one
C     we are looking for. This check is done in the function ISAME
C
      NBELEM = 0
C
      DO 60 IELEM = 1 , NELEM
          DO 60 I = 1 , 3
             II = ICELNOD(I,IELEM)
             JJ = ICELNOD(ICYCL(I+1,3),IELEM)
             KK = ICELNOD(ICYCL(I+2,3),IELEM)
             ICELCEL(I,IELEM) = 
     +ISAME2D( IELEM, JA(IA(JJ)), DEGREE(JJ), JA(IA(KK)), DEGREE(KK))
      IF ( ICELCEL(I,IELEM) .EQ. 0 )NBELEM=NBELEM+1
   60 CONTINUE 

      WRITE(6,*)'There seem to be ',NBELEM,' boundary elements'
      RETURN
      END


      INTEGER FUNCTION ISAME2D( IELEM, IA, LIA, IB, LIB ) 

      IMPLICIT NONE

      INTEGER IA(*),IB(*),LIA,LIB,I,J,IELEM

      DO 3 I = 1, LIA
          DO 2 J = 1, LIB
             IF( (IA(I) .EQ. IB(J)) .AND. IA(I) .NE. IELEM )THEN
                ISAME2D = IA(I)
                RETURN
             ENDIF
    2     CONTINUE
    3 CONTINUE

      ISAME2D = 0

      RETURN
      END
