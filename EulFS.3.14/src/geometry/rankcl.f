      SUBROUTINE RANKCL(INDX,A,JA,IA,LDA,NCL)
C
C     $Id: rankcl.f,v 1.1 2005/12/23 09:45:40 abonfi Exp $
C
C     A subroutine that re-orders the data structure
C     of the c-lines
C
      IMPLICIT NONE
      INTEGER LDA ,NCL
      DOUBLE PRECISION A(LDA,*)
      INTEGER INDX(*),JA(LDA,*),IA(NCL+1)
      INTEGER NNZR ,I
      NNZR = IA(NCL+1) -IA(1)
      DO I = 1,NNZR
           JA(4,I) = INDX(JA(4,I))
      ENDDO
      RETURN
      END
