!> \par Purpose
!>
!> This routine computes the time-accurate LW system scheme on one tetrahedron;
!>
!> it differs from the LDA version 1 system scheme since here the matrix \f$ K_i^+ \f$
!> multiplies \f$ \left( \sum_{j=1}^{d+1} K_j^+ \right)^{-1} \f$ from the left, i.e.
!>
!> the distribution matrix for the time-accurate LW scheme reads:
!> \f[
!> B_i = 
!> \f]
!>
!>
!> @param[in] MATRIXSPLITTER is the subroutine use to split the Jacobian matrix
!> @param[in] W is the set of NORDER dependent variables in the NOFVERT vertices of the current cell
!> @param[out] NODRES are the signals sent to each of the NOFVERT vertices
!> @param[in,out] TSTEP is the elemental contribution to the (inverse of the) time step
!> @param[out] BETA the NOFVERT distribution matrices
!> @param[out] STIFC the elemental contribution to the implicit matrix (computed only if MATRIX_ASSEMBLY is true)
!> @param[in] NORDER is the order of the hyperbolic system being solved, it could be < than NDOF
!> @param[in] NDOF is the number of dofs in the gridpoints, also the leading dimension of W, STIFC, NODRES, TSTEP
!> @param[out] NOFVERT is number of vertices of the simplicial elements (= NDIM+1)
!> @param[in] VCN the NDIM cartesian component of the inward face normal to all NOFVERT vertices, scaled by its measure
!> @param[in] NDIM is the dimension of the space
!> @param[in] DFGHDU the NDIM cartesian components of the Jacobian matrix
!> @param[in] LDJ is the leading dimension of DFGHDU
!> @param[out] CELRES the elemental residual
!> @param[in] SOURCE the elemental source term
!> @param[in] IELEM is the current simplicial element
!> @param[in] MATRIX_ASSEMBLY when set == .TRUE. the STIFC matrix will be assembled
!> \author $Author: abonfi $
!> \version $Revision: 1.9 $
!> \date $Date: 2020/03/28 09:49:28 $
!> \bug Not sure this is ok for 3D
      SUBROUTINE LWTSYS_SCHEME(MATRIXSPLITTER,W,NODRES,TSTEP,BETA,STIFC,
     +                         NORDER,NDOF,NOFVERT,VCN,NDIM,DFGHDU,
     +                         LDJ,CELRES,SOURCE,IELEM,
     &                         MATRIX_ASSEMBLY)
C
C     $Id: LWaSys_scheme.F,v 1.9 2020/03/28 09:49:28 abonfi Exp $
C
C
C
C
C
      IMPLICIT NONE
      INCLUDE 'paramt.h'
      INCLUDE 'constants.h'
      INCLUDE 'time.h'
      INCLUDE 'flags.com'
      INCLUDE 'time.com'
C
C
cnag  EXTERNAL F07ADF,F07AEF,DSCAL
C
C
C
C     .. Parameters ..
C
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION VOLUME
      INTEGER IELEM,LDJ,NDIM,NOFVERT,NDOF,NORDER
      LOGICAL MATRIX_ASSEMBLY
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION DFGHDU(LDJ,*),NODRES(NDOF,NOFVERT),
     +                 TSTEP(NDOF,NOFVERT),SOURCE(*),CELRES(*),
     3                 BETA(*),
     4                 STIFC(NDOF,NDOF,NOFVERT,NOFVERT),
     +                 VCN(NDIM,NOFVERT),W(NDOF,NOFVERT)
C     ..
C     .. Subroutine Arguments ..
      EXTERNAL MATRIXSPLITTER
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION GALERKIN,TEMP,HELP
      DOUBLE PRECISION CFLELEM
      INTEGER I,INFO,IROW,IVAR,IVERT,JVERT,J,N,M,ORDSQR,IADD,JADD
      LOGICAL LFLAG
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION 
     + KNEG(MAX_NOFVAR_SQR*MAXNOFVERT),KPOS(MAX_NOFVAR_SQR*MAXNOFVERT),
     +    K(MAX_NOFVAR_SQR*MAXNOFVERT),
     + VRIGHT(MAX_NOFVAR_SQR),VLEFT(MAX_NOFVAR_SQR),
     & WKSP2(MAXNOFVAR,MAXNOFVAR),
     & WNEG(MAXNOFVAR),WPOS(MAXNOFVAR),WR(MAXNOFVAR)
C     ..
C     .. External Functions ..
      LOGICAL NULLMAT,UNITMAT
      EXTERNAL NULLMAT,UNITMAT
C     ..
C     .. External Subroutines ..
      EXTERNAL DGEMM,DGEMV,DGETRF,DGETRS,DINIT,DSCAL,MATSUM,R8Mat_Print
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC DABS
C     ..
C     .. Statement Functions ..
      INTEGER IDX
      IDX(I,J,N,M) = (((J-1)*M+I-1)*N*N)+1
C     ..
C
      IF(.NOT.LTIME)THEN
         WRITE(*,*) 'This version of the LW scheme should be used for ti
     &me-accurate only'
         CALL EXIT(1)
      ELSE
         IF(DUALTS)THEN
            WRITE(*,*) 'The LW scheme should be used WITHOUT dual time-s
     &tepping'
            CALL EXIT(1)
         ENDIF
      ENDIF
C
      ORDSQR = NORDER*NORDER
      DO 3 I = 1,2*NORDER
          CELRES(I) = ZERO
    3 CONTINUE
C
C --------------- Debugging code starts here ---------------
C
      IF (ICHECK.NE.0) THEN
          CALL DINIT(MAX_NOFVAR_SQR,ZERO,WKSP2,1)
      ENDIF
C
C --------------- Debugging code ends here ---------------
C
      GALERKIN = ONE / NOFVERT
      CFLELEM = HALF*DTVOL ! not sure this is ok for 3D
!     IF(NDIM.EQ.3)STOP 'LW scheme has not been checked in 3D'
C
C     Loop over the vertices of the element ..
C
      DO 10 IVERT = 1,NOFVERT
C
C     The matrix IS split into its positive and negative parts
C     (this could be avoided)
C
C     probabilmente posso usare .FALSE. 
C
          IADD = ORDSQR*(IVERT-1)+1
          CALL MATRIXSPLITTER(IELEM,NDIM,NORDER,VCN(1,IVERT),DFGHDU,LDJ,
     +                        K(IADD),KPOS(IADD),KNEG(IADD),
     &                        VLEFT,VRIGHT,NORDER,
!    +                        WR,WPOS,WNEG,.TRUE.)
     +                        WR,WPOS,WNEG,.FALSE.)
C
C       Timestep ..
C
          DO 32 IVAR = 1,NORDER
              TSTEP(IVAR,IVERT) = TSTEP(IVAR,IVERT) + WPOS(IVAR)
   32     CONTINUE
C
C       the residual is computed as Sum_j K(j) * U(j) ..
C
         CALL DGEMV('NoTranspose',NORDER,NORDER,ONE,K(IADD),NORDER,
     +   W(1,IVERT),1,ONE,CELRES,1)
C
   10 CONTINUE ! end loop over vertices
C
C       .. Loop again over all vertices to compute the distribution
C               matrices (overwritten on KPOS)  ..
C
      DO 30 IVERT = 1,NOFVERT
C
C     compute the distribution matrices
C
         IADD = ORDSQR*(IVERT-1)+1
         CALL DCOPY(ORDSQR,K(IADD),1,BETA(IADD),1)
C
C        .. add the Galerkin term
C
         CALL DSCAL(ORDSQR,CFLELEM,BETA(IADD),1)
         DO 33 IROW = 1, NORDER
               JADD = IADD + (IROW-1)*NORDER + IROW -1
               BETA(JADD) = GALERKIN + BETA(JADD)
   33    CONTINUE
!         CALL R8Mat_Print('General',' ',NORDER,NORDER,BETA(IADD),
!    +                 NORDER,'  |B| matrix ',INFO)
C
C       NODRES(IVERT) = -BETA(IVERT)*CELRES
C
         CALL DGEMV('NoTranspose',NORDER,NORDER,MONE,BETA(IADD),NORDER,
     +   CELRES,1,ZERO,NODRES(1,IVERT),1)
C
C
C --------------- Debugging code starts here ---------------
C
          IF (ICHECK.NE.0) THEN
C
C       .. Sum up the distribution matrices ..
C
              CALL MATSUM(WKSP2,MAXNOFVAR,BETA(IADD),NORDER,
     +        NORDER,NORDER)
C
          ENDIF
C
C --------------- Debugging code ends here ---------------
C
   30 CONTINUE ! End loop over vertices
C
C       the source term is distributed using Galerkin
C
      DO IROW = 1,NORDER
         HELP = SOURCE(IROW)/REAL(NOFVERT)
         DO IVERT = 1, NOFVERT
            NODRES(IROW,IVERT) = NODRES(IROW,IVERT) -HELP
         ENDDO  
      ENDDO
C
      DO IROW = 1, NORDER
         CELRES(IROW) = CELRES(IROW) + SOURCE(IROW)
      ENDDO 
C
C --------------- Debugging code starts here ---------------
C
      include 'test1.inc'
C
C --------------- Debugging code ends here ---------------
C
      IF (MATRIX_ASSEMBLY) THEN
         STOP 'LW scheme should be used with explicit time-stepping'
C
C     .. Assembling the element stiffness matrix for the LW scheme ..
C
         DO 31 JVERT = 1,NOFVERT
C
             JADD = IDX(JVERT,1,NORDER,NOFVERT)
             DO 31 IVERT = 1,NOFVERT
C
C     ... C_{ij}^{LW} = - Beta_i^{+} K_j
C         Beta_i^{+} had been written in BETA
C
               IADD = IDX(IVERT,1,NORDER,NOFVERT)
               CALL DGEMM('N','N',NORDER,NORDER,NORDER,MONE,BETA(IADD),
     +         NORDER,K(JADD),NORDER,ZERO,STIFC(1,1,IVERT,JVERT),NDOF)
C
   31    CONTINUE
C
C --------------- Debugging code starts here ---------------
C
      include 'test2.inc'
      ENDIF
C
C --------------- Debugging code ends here ---------------
C
      RETURN
      END
