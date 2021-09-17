      subroutine testlhs(A,ZRoe,RHS,Y,V,NODCODE,NPOIN,NOFVAR)
C
C
C     this rouine was written to check the correctness of
C     the Picard linearization, i.e. that
C     rhs = - A*z
C     this should be true up to machine accuracy
C     It is no longer true on far field boundaries
C     due to the presence of a "source" term in the r.h.s 
C
C     $Id: test.F,v 1.15 2020/04/23 09:37:34 abonfi Exp $
C
C
#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscmat.h"
      use petscvec
      use petscmat
C
      IMPLICIT NONE
C
      INCLUDE 'constants.h'
C
      INTEGER NOFVAR,NPOIN
      INTEGER I,J,IFAIL
      INTEGER NODCODE(NPOIN)
      Mat A
      Vec RHS,VecY,VecV,ZRoe
      DOUBLE PRECISION Y(NOFVAR,NPOIN),V(NOFVAR,NPOIN)
      DOUBLE PRECISION ERR,dnrm2
      integer istkgt
      integer istak(1)
      double precision dstak(1)
      common/cstak/dstak
      equivalence(istak(1),dstak(1))
      character*5 file(4)
      data file /'ftn22','ftn23','ftn24','ftn25'/
C
C     Create VecV to store a cop of RHS
C
      CALL VecCreateSeqWithArray(MPI_COMM_SELF,NOFVAR,NPOIN*NOFVAR,
     +                           V,VecV,IFAIL)
      WRITE(6,*)'VecCreateSeqWithArray has ret IFAIL = ',IFAIL
C
C     copy rhs into V
C
      CALL VecCopy(RHS,VecV,IFAIL)
      WRITE(6,*)'VecCopy has ret IFAIL = ',IFAIL
C
C     create a vector y to store A*Zroe
C
      CALL VecCreateSeqWithArray(MPI_COMM_SELF,NOFVAR,NPOIN*NOFVAR,
     +                           Y,VecY,IFAIL)
      WRITE(6,*)'VecCreateSeqWithArray has ret IFAIL = ',IFAIL
C
C     compute y = A*v, should be equal to -RHS
C
      CALL MatMult(A,ZRoe,VecY,IFAIL)
      WRITE(6,*)'MatMult has ret IFAIL = ',IFAIL
C     CALL VecNorm(VecY,NORM_2,ERR,IFAIL)
C     write(6,*)' y = A * z ',ERR 

C
      WRITE(6,FMT=*)'Writing test output to ftn22'
C
      DO 1 I = 1, NPOIN
         WRITE(22,*)'Node # ',I,' Code # ',NODCODE(I)
C        WRITE(22,100)(Y(J,I),J=1,NOFVAR)
         DO 1 J = 1, NOFVAR
         WRITE(22,100)Y(J,I),V(J,I),Y(J,I)+V(J,I)
    1 CONTINUE
C
C     Do 2 I = 1, NOFVAR 
C        CALL DCOPY(NPOIN,Y(I,1),NOFVAR,Z,1)
C     CALL SOLZNE(file(I),Z,1,NPOIN,'w')
C   2 CONTINUE
C
C     Uncomment the following if you want to have
C     y = rhs - A*z
C
      CALL VecAXPY(VecY,ONE,VecV,IFAIL)
      WRITE(6,*)'VecAXPY has ret IFAIL = ',IFAIL
      CALL VecNorm(VecY,NORM_2,ERR,IFAIL)
      WRITE(6,*)'VecNorm has ret IFAIL = ',IFAIL
      write(6,*)'norm of RHS - C*z ',ERR 
C
      CALL VecDestroy(VecV,IFAIL)
      CALL VecDestroy(VecY,IFAIL)
      STOP 
C
      RETURN
  100 FORMAT(4(E14.8,1X))
      END
