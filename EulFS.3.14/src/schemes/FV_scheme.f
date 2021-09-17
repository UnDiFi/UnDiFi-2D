!> \par Purpose
!>
!> This routine computes the scalar version of the Finite Volume (FV) scheme on one 
!> triangle (2D) or tetrahedron (3D)
!>
!> @param[in] IELEM the current simplicial element
!> @param[in] VCN the NDIM cartesian component of the inward face normal to all NOFVERT vertices, scaled by its measure
!> @param[in] ADVECTION the NDIM cartesian component of the advection speed
!> @param[out] CELRES the elemental residual
!> @param[in] SOURCE the volume integral of the elemental source term
!> @param[in] Q the NOFVERT values of the dependent variable
!> @param[in,out] DT elemental contribution to the (inverse of the) time step
!> @param[out] NODRES the signals sent to each of the NOFVERT vertices
!> @param[out] BETA the distribution matrices
!> @param[out] STIFC the elemental contribution to the implicit matrix (computed only if MATRIX_ASSEMBLY is true)
!> @param[in] NDIM is the dimension of the space
!> @param[out] NOFVERT is number of vertices of the simplicial elements (= NDIM+1)
!> @param[in] MATRIX_ASSEMBLY when set == .TRUE. the STIFC matrix will be assembled
!>
      SUBROUTINE FV_scheme(IELEM,VCN,ADVECTION,CELRES,SOURCE,Q,DT,
     &                     NODRES,BETA,STIFC,NDIM,NOFVERT,
     &                     MATRIX_ASSEMBLY)
C
C     $Id: FV_scheme.f,v 1.6 2013/08/22 15:10:57 abonfi Exp $
C
      IMPLICIT NONE
C
C
      include 'paramt.h'
      include 'constants.h'
C
C     .. Scalar Arguments ..
C
      INTEGER IELEM,NDIM,NOFVERT
      DOUBLE PRECISION CELRES,SOURCE
      LOGICAL MATRIX_ASSEMBLY
C
C     .. Array Arguments ..
C
      DOUBLE PRECISION VCN(NDIM,NOFVERT),ADVECTION(NDIM),DT(NOFVERT),
     &NODRES(NOFVERT),BETA(NOFVERT),STIFC(NOFVERT,NOFVERT),Q(NOFVERT)
C
C     .. Local Scalars ..
C
      INTEGER INODE,IVERT,IVAR,I,J,Ni,Nj
      DOUBLE PRECISION EdgeRes,dK,DENOM
C
C     .. Local Arrays ..
C
      DOUBLE PRECISION K(MAXNOFVERT)
C
C     .. External Functions ..
C
      DOUBLE PRECISION DDOT
      EXTERNAL         DDOT
C
C     .. Intrinsic Functions ..
C
      INTRINSIC DBLE
C
C     .. Executable Statements ..
C
      CELRES = ZERO  ! residual = - fluctuation
      DENOM = ONE/REAL(NOFVERT)
C
      IF(MATRIX_ASSEMBLY)CALL DINIT(NOFVERT*NOFVERT,ZERO,STIFC,1)
C
      DO 10 IVERT = 1 , NOFVERT
c
c Dotting advection speed with normal
c
        K(IVERT) = DDOT( NDIM , VCN(1,IVERT) , 1 , ADVECTION , 1 ) 
     &  / REAL(NDIM)
C
   10 CONTINUE
C
C     .. Loop over the edges of the element ..
C
      DO 15 I = 1 , NOFVERT
        DO 15 J = I+1 , NOFVERT
C
c Dotting advection speed with edge normal
c
          dK = K(J) - K(I)
C
C here we compute \phi_{ji} = (k_j-k_i)^+ (u_j-u_i)
C
          EdgeRes = dK * ( Q(J) - Q(I) ) * DENOM
c
          CELRES = CELRES + EdgeRes
c
          IF( dK .GT. ZERO )THEN
            DT(J) = DT(J) + dK*DENOM
            NODRES(J) = NODRES(J) - EdgeRes
            IF(MATRIX_ASSEMBLY)STIFC(J,I) = STIFC(J,I) + dK*DENOM
          ELSE
            DT(I) = DT(I) - dK*DENOM
            NODRES(I) = NODRES(I) - EdgeRes
            IF(MATRIX_ASSEMBLY)STIFC(I,J) = STIFC(I,J) - dK*DENOM
          ENDIF
C
C
   15 CONTINUE ! End loop over the edges
C
      IF(MATRIX_ASSEMBLY)THEN
        DO 7 I = 1, NOFVERT
           DO 9 J = 1,NOFVERT
              IF(J.EQ.I)GOTO 9
              STIFC(I,I) = STIFC(I,I) - STIFC(I,J)
    9      CONTINUE
    7   CONTINUE
      ENDIF
C
      RETURN
      END
C
