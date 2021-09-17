      SUBROUTINE ROTATE(ZROE,NOFVAR,NPNOD)
      IMPLICIT NONE
C
C     $Id: rotate.f,v 1.4 2013/04/30 07:23:51 abonfi Exp $
C
C Subroutine for rotating velocities
C in the periodic AND annular case v,w
C velocities are stored in Z(IY,*) Z(IZ,*)
C it is assumed that x coincides with the
C axis of the turbomachine
C
      INCLUDE 'paramt.h'
      INCLUDE 'periodic.com'
      INCLUDE 'dofs.com'
C
C     .. Scalar Arguments ..
      INTEGER NOFVAR,NPNOD
C     .. Array Arguments ..
      DOUBLE PRECISION ZROE(NOFVAR,*)
C     ..
C     .. Local Scalars ..
      INTEGER I
C     ..
C
      DO 3 I = 1, NPNOD
         CYY = ZROE(IY,I)
         CZZ = ZROE(IZ,I)
         ZROE(IY,I) = COSALPHA*CYY-SINALPHA*CZZ
         ZROE(IZ,I) = SINALPHA*CYY+COSALPHA*CZZ
    3 CONTINUE
C
      RETURN
 
      END
      SUBROUTINE INIQMAT(AMAT,NOFVAR)
C
C     Initialize the rotation matrix
C     to be used for annular cascade flows
C     this routine should be called only once
C
C     AMAT is nothing but QMAT
C     the reason for passing QMAT as AMAT
C     is that QMAT is stored as a 1D array
C     while here we want to access AMAT as
C     a 2D array
C
      IMPLICIT NONE
      INCLUDE 'paramt.h'
      INCLUDE 'constants.h'
      INCLUDE 'periodic.com'
      INCLUDE 'dofs.com'
C
      INTEGER NOFVAR
      DOUBLE PRECISION AMAT(NOFVAR,*)
      INTEGER I,J
C
      DO J = 1,NOFVAR
         DO I = 1,NOFVAR
            AMAT(I,J) = ZERO
         ENDDO
      ENDDO
      DO I = 1,NOFVAR
         AMAT(I,I) = ONE
      ENDDO
      AMAT(IY,IY) = COSALPHA
      AMAT(IY,IZ) =-SINALPHA
      AMAT(IZ,IY) = SINALPHA
      AMAT(IZ,IZ) = COSALPHA
      RETURN
      END
C
      SUBROUTINE ROTATECIJ(A,NOFVAR,NOFVERT,IMAX,JMAX)
C
C     Rotate the Jacobian due to periodic
C     bcs in an annular cascade
C
      IMPLICIT NONE
      INCLUDE 'paramt.h'
      INCLUDE 'constants.h'
      INCLUDE 'periodic.com'
C
      INTEGER NOFVAR,NOFVERT,IMAX,JMAX
      DOUBLE PRECISION A(NOFVAR,NOFVAR,NOFVERT,NOFVERT)
      DOUBLE PRECISION WKSP(MAX_NOFVAR_SQR)
      INTEGER IVERT,JVERT
C
C     loop over vertices
C
      DO IVERT = 1,IMAX
         IF( PFLAG(IVERT) )THEN 
            DO JVERT = 1,JMAX
C
C     does C_{il} := Q^t C_{il} forall l
C
               CALL DCOPY(NOFVAR*NOFVAR,A(1,1,IVERT,JVERT),1,WKSP,1)
               CALL DGEMM('Transpose','No Transpose',NOFVAR,NOFVAR,
     &         NOFVAR,ONE,QMAT,NOFVAR,WKSP,NOFVAR,ZERO,
     &         A(1,1,IVERT,JVERT),NOFVAR)
            ENDDO
C
C     does C_{ii} := C_{ii} Q
C
            CALL DCOPY(NOFVAR*NOFVAR,A(1,1,IVERT,IVERT),1,WKSP,1)
            CALL DGEMM('No Transpose','No Transpose',NOFVAR,NOFVAR,
     &      NOFVAR,ONE,WKSP,NOFVAR,QMAT,NOFVAR,ZERO,
     &      A(1,1,IVERT,IVERT),NOFVAR)
         ENDIF
      ENDDO
      RETURN
      END
