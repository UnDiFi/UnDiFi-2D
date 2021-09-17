      SUBROUTINE LDAS_FAST(MATRIXSPLITTER,W,NODRES,TSTEP,BETA,STIFC,
     +                     NORDER,NDOF,NOFVERT,VCN,NDIM,DFGHDU,LDJ,
     +                     CELRES,SOURCE,IELEM,MATRIX_ASSEMBLY)
C
C     $Id: LDAS_fast.F,v 1.9 2020/03/28 09:49:28 abonfi Exp $
C
C#define BLAS_LAPACK
C#define DEBUG
#define UNROLL
C
C
C NORDER                is the order of the system to solve for, i.e.
C                       the order of the matrix dFGHdU.
C FrstEq                is the first among the NDOF variable which
C                       appears into the system of order NORDER.
C IELEM                 is the current element.
C dFGHdU(LDJ,*) is the Jacobian Matrix of the system.
C LDJ                   is the leading dimension of dFGHdU.
C W                  stores by columns the NORDER variables of the
C                       NOFVERT vertices.
C LNODRES                   is the leading dimension of W.
C CELRES[1:2*NDOF]  OUT Real
C CELRES[1:NDOF]    stores the residual computed by the Matrix scheme
C                       as \sum K_j U_j (explicit part of the scheme)
C CELRES[NDOF+1:2*NDOF]
C                       stores the residual computed by the Matrix scheme
C                       as \sum C_{ij} U_j (implicit part of the scheme)
C MatrixSplitter        is the procedure used to compute the eigenvector
C                       decomposition of the matrix dFGHdU.
C TSTEP                 is the nodal timestep.
C DU(NDOF,*)          is the nodal residual.
C
C
C
C This routine computes the LDA system scheme on one tetrahedron
C
C
      IMPLICIT NONE
      INCLUDE 'paramt.h'
      INCLUDE 'constants.h'
      INCLUDE 'flags.com'
C
C
cnag  EXTERNAL F07ADF,F07AEF,DSCAL
C
C
C
C     .. Parameters ..
C     ..
C     .. Scalar Arguments ..
      INTEGER IELEM,LDJ,NDIM,NOFVERT,NDOF,NORDER
      LOGICAL MATRIX_ASSEMBLY
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION DFGHDU(LDJ,*),NODRES(NDOF,NOFVERT),
     +                 TSTEP(NDOF,NOFVERT),SOURCE(*),
     +                 CELRES(*),BETA(*),
     4                 STIFC(NDOF,NDOF,NOFVERT,NOFVERT),
     +                 VCN(NDIM,NOFVERT),W(NDOF,NOFVERT)
C     ..
C     .. Subroutine Arguments ..
      EXTERNAL MATRIXSPLITTER
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION SUM,TEMP
      INTEGER I,INFO,IROW,IVAR,IVERT,J,JCOL,IADD,JADD,ORDSQR
      LOGICAL LFLAG
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION K(MAXNOFVAR*MAXNOFVAR*VMAX),
     +KNEG(MAXNOFVAR*MAXNOFVAR*VMAX),
     +SUM_K_POS(MAXNOFVAR*MAXNOFVAR),VLEFT(MAXNOFVAR*MAXNOFVAR),
     +VRIGHT(MAXNOFVAR*MAXNOFVAR),WKSP1(MAXNOFVAR),
     +WKSP2(MAXNOFVAR,MAXNOFVAR),WNEG(MAXNOFVAR),
     +WPOS(MAXNOFVAR),WR(MAXNOFVAR),WORK(MAXNOFVAR)
      EQUIVALENCE(WORK(1),WPOS(1))
      INTEGER IPIV(MAXNOFVAR)
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
      ORDSQR = NORDER*NORDER
      CALL DINIT(ORDSQR,ZERO,SUM_K_POS,1)
      DO 3 I = 1,2*NDOF
          CELRES(I) = ZERO
    3 CONTINUE
C
C --------------- Debugging code starts here ---------------
C
#ifdef DEBUG
      IF (ICHECK.NE.0) THEN
          CALL DINIT(MAXNOFVAR,ZERO,WKSP1,1)
          CALL DINIT(MAXNOFVAR*MAXNOFVAR,ZERO,WKSP2,1)
      ENDIF
#endif
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
     +                       K(IADD),BETA(IADD),KNEG(IADD),VLEFT,VRIGHT,
     +                       NORDER,WR,WPOS,WNEG,.TRUE.)
C
C       The positive jacobians are accumulated in SUM_K_POS ..
C
          CALL DAXPY(ORDSQR,ONE,BETA(IADD),1,SUM_K_POS,1)
C
C       Timestep ..
C
          DO 32 IVAR = 1,NORDER
              TSTEP(IVAR,IVERT) = TSTEP(IVAR,IVERT) + WPOS(IVAR)
   32     CONTINUE
C
C       the residual is computes as Sum_j K(j) * U(j) ..
C
#ifdef BLAS_LAPACK
         CALL DGEMV('N',NORDER,NORDER,ONE,K(IADD),NORDER,
     +              W(1,IVERT),1,ONE,CELRES,1)
#else
         CALL MY_GEMV('N',NORDER,NORDER,ONE,K(IADD),NORDER,
     +              W(1,IVERT),1,ONE,CELRES,1)
#endif
C
C  add the source term
C
   10 CONTINUE
      DO 14 IVAR = 1,NORDER
         WORK(IVAR) = CELRES(IVAR) + SOURCE(IVAR)
   14 CONTINUE
C
C       .. LU factorization of SUM_K_POS ..
C
#ifdef BLAS_LAPACK
      CALL DGETRF(NORDER,NORDER,SUM_K_POS,NORDER,IPIV,INFO)
#else
#ifdef UNROLL
      CALL xxDECO(SUM_K_POS,NORDER)
#else
      CALL LUDECO(SUM_K_POS,NORDER)
#endif
#endif
C
C       .. Loop again over all vertices to compute the distribution
C               matrices (overwritten on BETA)  ..
C
      DO 30 IVERT = 1,NOFVERT
C
C       .. solution ..
C
          IADD = (IVERT-1)*ORDSQR+1
#ifdef BLAS_LAPACK
          CALL DGETRS('N',NORDER,NORDER,SUM_K_POS,NORDER,IPIV,
     +                BETA(IADD),NORDER,INFO)
#else
CCC       DO 7 J = 1,NORDER
CCC          JADD = IADD + (J-1)* NORDER
CCC          CALL xxSOLV(SUM_K_POS,BETA(JADD),BETA(JADD),NORDER)
CCC 7     CONTINUE
#ifdef UNROLL
          CALL MLUSOLX(SUM_K_POS,BETA(IADD),BETA(IADD),NORDER)
#else
          CALL MLUSOL(SUM_K_POS,BETA(IADD),BETA(IADD),NORDER)
#endif
#endif
C
C         NODRES(IVERT) = -BETA(IVERT)*CELRES
C
#ifdef BLAS_LAPACK
          CALL DGEMV('N',NORDER,NORDER,-ONE,BETA(IADD),NORDER,WORK,
     +                1,ZERO,NODRES(1,IVERT),1)
#else
          CALL MY_GEMV('N',NORDER,NORDER,-ONE,BETA(IADD),NORDER,
     +                WORK,1,ZERO,NODRES(1,IVERT),1)
#endif
C
C
C --------------- Debugging code starts here ---------------
C
#ifdef DEBUG
          IF (ICHECK.NE.0) THEN
C
C       .. Sum up the distribution matrices ..
C
          CALL MATSUM(WKSP2,MAXNOFVAR,BETA(IADD),NORDER,NORDER,NORDER)
C
          ENDIF
#endif
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
C         Delta_i^{+} had been overwritten onto BETA
C
              IADD = (I-1)*ORDSQR+1
              CALL DGEMM('N','N',NORDER,NORDER,NORDER,MONE,BETA(IADD),
     +             NORDER,K(JADD),NORDER,ZERO,STIFC(1,1,I,J),NDOF)
C
   31 CONTINUE
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
