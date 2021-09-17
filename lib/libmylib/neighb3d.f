
      SUBROUTINE NEIGHB3D(ICELNOD,ICELCEL,NELEM,DEGREE,IA,JA,LENJA,
     +NPOIN)

      IMPLICIT NONE

      INTEGER ICELNOD(4,*),
     +        ICELCEL(4,*),DEGREE(*),JA(*),IA(*)

      INTEGER NELEM,LENJA,NPOIN,NBELEM,MINDEG,MAXDEG
      INTEGER IELEM,INODE,I,II,JJ,KK,LL,ISUM,IFAIL
      INTEGER  ICYCL,ISAME
      EXTERNAL ICYCL,ISAME
C
      DO I = 1 , NPOIN
       DEGREE(I) = 0
      ENDDO
C
C      counts the number of elements surrounding a given node
C      and store them in DEGREE
C
      DO 100 IELEM = 1 , NELEM
         DO 100 I = 1 , 4
            INODE = ICELNOD(I,IELEM)
            DEGREE(INODE) = DEGREE(INODE) + 1
  100 CONTINUE

      ISUM = 0      
      MINDEG = 10000
      MAXDEG = 0
      IFAIL = 0
      DO 90 I = 1 , NPOIN
         ISUM = ISUM + DEGREE(I)
         MINDEG = MIN(MINDEG,DEGREE(I))
         MAXDEG = MAX(MAXDEG,DEGREE(I))
         if(DEGREE(I).LE.0)THEN
            WRITE(6,*)'Non positive degree ',DEGREE(I),
     + ' vertex ',i
            IFAIL = IFAIL+1
         ENDIF
   90 CONTINUE
      WRITE(6,*)'Min/Avg/Max Vertex degree ',MINDEG,ISUM/NPOIN,MAXDEG
*     write(16,*)(DEGREE(I),i=1,NPOIN)

      IF(IFAIL.NE.0)THEN
         WRITE(6,*)IFAIL,' vertices have non positive degree'
         STOP
      ENDIF
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
         DO 70 I = 1 , 4
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
          DO 60 I = 1 , 4
             II = ICELNOD(I,IELEM)
             JJ = ICELNOD(ICYCL(I+1,4),IELEM)
             KK = ICELNOD(ICYCL(I+2,4),IELEM)
             LL = ICELNOD(ICYCL(I+3,4),IELEM)
             ICELCEL(I,IELEM) = 
     +ISAME( IELEM, JA(IA(JJ)), DEGREE(JJ), JA(IA(KK)), DEGREE(KK),
     +              JA(IA(LL)), DEGREE(LL))
      IF ( ICELCEL(I,IELEM) .EQ. 0 )NBELEM=NBELEM+1
   60 CONTINUE 

      WRITE(6,*)'There seem to be ',NBELEM,' boundary elements'
      RETURN
      END
      INTEGER FUNCTION ISAME( IELEM, IA, LIA, IB, LIB, IC, LIC ) 

      IMPLICIT NONE

      INTEGER IA(*),IB(*),IC(*),LIA,LIB,LIC,I,J,K,IELEM

      DO 3 I = 1, LIA
          DO 2 J = 1, LIB
             IF( (IA(I) .EQ. IB(J)) .AND. IA(I) .NE. IELEM )THEN
                ISAME = IA(I)
C
C     check if ISAME is present also in the third list
C
                DO 5 K = 1, LIC
                   IF( IC(K) .EQ. ISAME )RETURN
    5           CONTINUE
             ENDIF
    2     CONTINUE
    3 CONTINUE

      ISAME = 0

      RETURN
      END
