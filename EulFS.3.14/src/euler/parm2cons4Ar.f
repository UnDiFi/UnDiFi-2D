      SUBROUTINE PARM2CONS4Ar(ZROE,DUDZ,NOFVAR,NDIM)
C
C     $Id: parm2cons4Ar.f,v 1.3 2013/04/30 07:42:54 abonfi Exp $
C
      IMPLICIT NONE
C
C     transformation matrix from
C     parameter vector to conserved variables
C
C
C
C
C     Assembles the dUdZ matrix ...
C
C     .. Parameters ..
      INCLUDE 'paramt.h'
      INCLUDE 'constants.h'
      INCLUDE 'plasma.h'
      INCLUDE 'dofs.com'
      INCLUDE 'three.com'
      INCLUDE 'four.com'
C     ..
C     .. Scalar Arguments ..
      INTEGER NDIM,NOFVAR
      INTEGER I,J
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION DUDZ(NOFVAR,NOFVAR),ZROE(NOFVAR)
      DOUBLE PRECISION DPDZS(NSP),DPDZH,DPDZU(3),SUMDPZS ,HELP
C
C
      SUMDPZS = ZERO
      DO I = 1 , NSP
        SUMDPZS = SUMDPZS + DR(I) * ZAVG(I)
      ENDDO
C
C     .. Local Scalars ..
C     ..
C
C     PRESSURE JACOBIAN  DP/DZ 
C
C     singles species
C
      DO I = 1 , NSP
         DPDZS(I) = ONE/(ONE+DE) * (SUMDPZS + DR(I)*SQRTR + DE*ZAVG(IE)
     &   + DM(1)*ZAVG(IX) + DM(2)*ZAVG(IY))
         IF (NDIM.EQ.3) THEN
           DPDZS(I) = DPDZS(I) +  DM(3)*ZAVG(IZ)/(ONE+DE)   
         ENDIF  
      ENDDO
      HELP = SQRTR/(ONE+DE)
C
C     energy	
C
      DPDZH = DE*HELP
C
C     x-momentum
C
      DPDZU(1) = DM(1)*HELP
C
C     y-momentum
C      
      DPDZU(2) = DM(2)*HELP
C
C     z-momentum
C
      IF (NDIM.EQ.3) THEN
         DPDZU(3) = DM(3)*HELP
      ENDIF 
C
C
C     JACOBIAN DU/DZ
C
C     single species 	
C
      DO I = 1 , NSP
        DO J = 1 , NSP  
          DUDZ(I,J) = ZROE(I)
          IF (I.EQ.J) THEN
            DUDZ(I,J) = DUDZ(I,J) + SQRTR
          ENDIF
        ENDDO
      ENDDO

      DO I = 1 , NSP
        DUDZ(I,IE) = ZERO
        DUDZ(I,IX) = ZERO
        DUDZ(I,IY) = ZERO
      ENDDO
C
C     energy
C
      DO I = 1 , NSP
        DUDZ(IE,I) = ZAVG(IE) - DPDZS(I)
      ENDDO
      DUDZ(IE,IE) = SQRTR - DPDZH
      DUDZ(IE,IX) = - DPDZU(1)
      DUDZ(IE,IY) = - DPDZU(2)
C
C     x-momentum
C
      DO I = 1 , NSP
        DUDZ(IX,I) = ZAVG(IX)               
      ENDDO    
      DUDZ(IX,IE) = ZERO
      DUDZ(IX,IX) = SQRTR
      DUDZ(IX,IY) = ZERO
C
C     y-momentum
C
      DO I = 1 , NSP
        DUDZ(IY,I) = ZAVG(IY)
      ENDDO
      DUDZ(IY,IE) = ZERO
      DUDZ(IY,IX) = ZERO
      DUDZ(IY,IY) = SQRTR
C
      IF (NDIM.EQ.2) RETURN
C     
C     
      DO I = 1 , NSP
        DUDZ(I,IZ) = ZERO
      ENDDO
      DUDZ(IE,IZ) = -DPDZU(3)
      DUDZ(IX,IZ) = ZERO
      DUDZ(IY,IZ) = ZERO
C
C     z-momentum
C
      DO I = 1 , NSP
        DUDZ(IZ,I) = ZAVG(IZ)
      ENDDO
      DUDZ(IZ,IE) = ZERO
      DUDZ(IZ,IX) = ZERO
      DUDZ(IZ,IY) = ZERO
      DUDZ(IZ,IZ) = SQRTR
C
C
      RETURN
C
      END
