!> \brief \b LDASYS_SCHEME
!> \par Purpose
!>
!> This routine computes the LDA version 1 system scheme on one tetrahedron;
!>
!> it differs from the LDA version 2 system scheme since here the matrix \f$ K_i^+ \f$
!> multiplies \f$ \left( \sum_{j=1}^{d+1} K_j^+ \right)^{-1} \f$ from the right, i.e.
!>
!> the distribution matrix for the LDA verion 1 scheme reads:
!> \f[
!> B_i = \left( \sum_{j=1}^{d+1} K_j^+ \right)^{-1} K_i^+
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
!> \version $Revision: 1.15 $
!> \date $Date: 2020/03/28 09:49:28 $
!> \warning Define \c FSPL_USE_PETSC_KERNELS if you want to use PETSc's kernels to invert matrix \f$ \left( \sum_{j=1}^{d+1} K_j^+\right) \f$ (the preferred choice)
!> \warning Define \c FSPL_USE_LAPACK if you want to use LAPACK to invert matrix \f$ \left( \sum_{j=1}^{d+1} K_j^+\right) \f$ (presumably much slower than the former)
      SUBROUTINE LDASYS_SCHEME(MATRIXSPLITTER,W,NODRES,TSTEP,BETA,STIFC,
     +           NORDER,NDOF,NOFVERT,VCN,NDIM,DFGHDU,LDJ,CELRES,
     +           SOURCE,IELEM,MATRIX_ASSEMBLY)
C
C     $Id: LDASys_scheme.F,v 1.15 2020/03/28 09:49:28 abonfi Exp $
C
#define FSPL_USE_PETSC_KERNELS
CCCC#define FSPL_USE_LAPACK
C#define DEBUG
C
C
      IMPLICIT NONE
      INCLUDE 'paramt.h'
      INCLUDE 'constants.h'
      INCLUDE 'flags.com'
C
C
C     .. Parameters ..
C     ..
C     .. Scalar Arguments ..
      INTEGER IELEM,LDJ,NDIM,NDOF,NOFVERT,NORDER
      DOUBLE PRECISION VOLUME
      LOGICAL MATRIX_ASSEMBLY
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION DFGHDU(LDJ,*),NODRES(NDOF,NOFVERT),
     +                 TSTEP(NDOF,NOFVERT),CELRES(*),SOURCE(*),
     +                 BETA(*),
     4                 STIFC(NDOF,NDOF,NOFVERT,NOFVERT),
     +                 VCN(NDIM,NOFVERT),W(NDOF,NOFVERT)
C     ..
C     .. Subroutine Arguments ..
      EXTERNAL MatrixSplitter
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION TEMP,HELP
      INTEGER INFO,IVAR,IVERT,I,J,IADD,JADD,ORDSQR,N,M,JVERT
      LOGICAL LFLAG
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION K(MAX_NOFVAR_SQR*MAXNOFVERT),
     2KNEG(MAX_NOFVAR_SQR*MAXNOFVERT),KPOS(MAX_NOFVAR_SQR*MAXNOFVERT),
     3SUM_K_POS(MAX_NOFVAR_SQR),
     4VLEFT(MAX_NOFVAR_SQR),VRIGHT(MAX_NOFVAR_SQR),
     5WNEG(MAXNOFVAR),WPOS(MAXNOFVAR),WR(MAXNOFVAR),WKSP(MAXNOFVAR)
      EQUIVALENCE(WKSP(1),WPOS(1))
C
#ifdef FSPL_USE_LAPACK
C
      INTEGER LWORK
      PARAMETER(LWORK=2*MAXNOFVAR)
      INTEGER IPIV(MAXNOFVAR)
      DOUBLE PRECISION WORK(LWORK)
#endif
C     ..
C     .. External Functions ..
      LOGICAL NULLMAT,UNITMAT
      EXTERNAL NULLMAT,UNITMAT
      DOUBLE PRECISION WKSP2(MAXNOFVAR,MAXNOFVAR)
C     DOUBLE PRECISION WKSP2(MAX_NOFVAR_SQR)
C     ..
C     .. External Subroutines ..
      EXTERNAL DGEMM,DGEMV,DGETRF,DGETRS,DINIT,DSCAL,MATSUM,R8Mat_Print
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC DABS
C     ..
C     .. Statement Functions ..
cold  INTEGER IDX
cold  IDX(I,J,N,M) = (((J-1)*M+I-1)*N*N)+1
C     ..
C     ..
      ORDSQR = NORDER*NORDER
      CALL DINIT(ORDSQR,ZERO,SUM_K_POS,1)
      DO 3 I = 1,2*NDOF
          CELRES(I) = ZERO
    3 CONTINUE
C
C --------------- Debugging code starts here ---------------
C
      IF (ICHECK.NE.0) THEN
          CALL DINIT(MAXNOFVAR*MAXNOFVAR,ZERO,WKSP2,1)
      ENDIF
C
C --------------- Debugging code ends here ---------------
C
C
C       .. Loop over the vertices of the element ..
C
      DO 10 IVERT = 1,NOFVERT
C
C       .. The matrix is split into its positive and negative parts
C
C       Numerically or Analitically ..
C
          IADD = (IVERT-1)*ORDSQR+1
          CALL MATRIXSPLITTER(IELEM,NDIM,NORDER,VCN(1,IVERT),DFGHDU,LDJ,
     +                       K(IADD),KPOS(IADD),BETA(IADD),VLEFT,VRIGHT,
     +                       NORDER,WR,WPOS,WNEG,.TRUE.)
C
C       The positive jacobians are accumulated in SUM_K_POS ..
C
          CALL DAXPY(ORDSQR,ONE,KPOS(IADD),1,SUM_K_POS,1)
C
C       Timestep ..
C
          DO 32 IVAR = 1,NORDER
              TSTEP(IVAR,IVERT) = TSTEP(IVAR,IVERT) + WPOS(IVAR)
   32     CONTINUE
C
C       the residual is computed as Sum_j K(j) * U(j) ..
C
         CALL DGEMV('N',NORDER,NORDER,ONE,K(IADD),NORDER,
     +              W(1,IVERT),1,ONE,CELRES,1)
C
   10 CONTINUE
C
C       .. Invert SUM_K_POS ..
C
#ifdef FSPL_USE_LAPACK
      CALL DGETRF(NORDER,NORDER,SUM_K_POS,NORDER,IPIV,INFO)
      CALL DGETRI(NORDER,SUM_K_POS,NORDER,IPIV,WORK,LWORK,INFO)
#else
      CALL SOLVEN( NORDER, SUM_K_POS )
#endif
C
      DO 14 IVAR = 1,NORDER
         WKSP(IVAR) = CELRES(IVAR)+SOURCE(IVAR)
   14 CONTINUE
C
C       .. Loop again over all vertices to compute the distribution
C               matrices (put into BETA)  ..
C
      DO 30 IVERT = 1,NOFVERT
C
C       .. solution ..
C
          IADD = (IVERT-1)*ORDSQR+1
          CALL DGEMM('NoTranspose','NoTranspose',NORDER,NORDER,NORDER,
     +             ONE,SUM_K_POS,NORDER,KPOS(IADD),NORDER,ZERO,
     &             BETA(IADD),NORDER)
C
C         NODRES(IVERT) = -BETA(IVERT)*CELRES
C
          CALL DGEMV('N',NORDER,NORDER,MONE,BETA(IADD),NORDER,WKSP,
     +                1,ZERO,NODRES(1,IVERT),1)
C
C
C --------------- Debugging code starts here ---------------
C
      include 'test0.inc'
C
C --------------- Debugging code ends here ---------------
C
C
   30 CONTINUE
C
C --------------- Debugging code starts here ---------------
C
      include 'test1.inc'
C
C --------------- Debugging code ends here ---------------
C
      IF (MATRIX_ASSEMBLY) THEN
C
C     .. Assembling the element stiffness matrix for the LDA scheme ..
C
         DO 31 J = 1,NOFVERT
C
            JADD = (J-1)*ORDSQR+1
            DO 31 I = 1,NOFVERT
C
C     ... C_{ij}^{LDA} = - Delta_i^{+} K_j
C         Delta_i^{+} had been overwritten onto KPOS
C
              IADD = (I-1)*ORDSQR+1
              CALL DGEMM('N','N',NORDER,NORDER,NORDER,MONE,BETA(IADD),
     +             NORDER,K(JADD),NORDER,ZERO,STIFC(1,1,I,J),NDOF)
C
   31    CONTINUE
      ENDIF
C
C --------------- Debugging code starts here ---------------
C
      include 'test2.inc'
C
C --------------- Debugging code ends here ---------------
C
      RETURN

      END
