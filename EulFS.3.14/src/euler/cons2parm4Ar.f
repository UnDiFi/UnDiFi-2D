      SUBROUTINE CONS2PARM4Ar(ZROE,dZdU,NDIM,NOFVAR)
      ENTRY MatdZdU4Ar(ZROE,dZdU,NDIM,NOFVAR)
C
      IMPLICIT NONE
      INCLUDE 'paramt.h'
C
      INCLUDE 'constants.h' 
      INCLUDE 'plasma.h'
C      
      INCLUDE 'dofs.com'
      INCLUDE 'three.com'
      INCLUDE 'four.com'
C
      INTEGER NDIM,NOFVAR
      INTEGER I,J
      DOUBLE PRECISION ZROE(NOFVAR),dZdU(NOFVAR,*)
      DOUBLE PRECISION ZRINV,ZRSQRINV,HELP
C
C
C     Assembles the dZdU matrix ...
C      WRITE(*,*)'SQRTR ZAVG',SQRTR
C
      SQRTR = ZERO
      DO I = 1,NSP
         SQRTR = SQRTR + ZROE(I)
      ENDDO      
C
C      write(*,*)'SQRTR ZROE',SQRTR
C
      ZRINV = ONE/SQRTR
      ZRSQRINV = ZRINV*ZRINV
C
      KINETIC = ZROE(IX)*ZROE(IX) + ZROE(IY)*ZROE(IY)
      IF (NDIM.EQ.3) KINETIC = KINETIC + ZROE(IZ)*ZROE(IZ)
      KINETIC = HALF*KINETIC*ZRSQRINV
!      write(*,*) KINETIC,ZROE(IX),ZROE(IY)
C
      DO I = 1 , NSP
!         write(*,*)I
!         write(*,*) DR(I)
         DR(I) = CHI(I) + KAPPA * KINETIC
!         write(*,*) DR(I)
      ENDDO
C
      DO I = 1 , NSP
        DO J = 1 , NSP 
            dZdU(I,J) = - ZROE(I)*HALF*ZRSQRINV
            IF(I.EQ.J) THEN
                dZdU(I,J) = dZdU(I,J) + ZRINV
            ENDIF
        ENDDO
      ENDDO
C
      DO J = 1 , NSP
        dZdU(IE,J) = DR(j) * ZRINV - HALF * ZROE(IE) * ZRSQRINV 
        dZdU(IX,J) = -HALF * ZROE(IX) * ZRSQRINV 
        dZdU(IY,J) = -HALF * ZROE(IY) * ZRSQRINV 
      ENDDO
      dZdU(IE,IE) = (DE+ONE) * ZRINV
      dZdU(IE,IX) =  -DE * ZROE(IX) * ZRSQRINV
      dZdU(IE,IY) =  -DE * ZROE(IY) * ZRSQRINV
C   
      dZdU(IX,IX) = ZRINV
C
      dZdU(IY,IY) = ZRINV
C
      IF (NDIM.EQ.2) RETURN
C
      dZdU(IE,IZ) = -DE * ZROE(IZ) * ZRSQRINV

      DO J = 1 , NSP
        dZdU(IZ,J) = -HALF * ZROE(IZ) * ZRSQRINV 
      ENDDO

      dZdU(IZ,IZ) = ZRINV
C
      RETURN

      END
