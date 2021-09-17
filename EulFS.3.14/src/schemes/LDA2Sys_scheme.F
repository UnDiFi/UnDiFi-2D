!> \brief \b LDA2SYS_SCHEME
!> \par Purpose
!>
!> This routine computes the LDA version 2 system scheme on one tetrahedron;
!>
!> it differs from the LDA version 1 system scheme since here the matrix \f$ K_i^+ \f$
!> multiplies \f$ \left( \sum_{j=1}^{d+1} K_j^+ \right)^{-1} \f$ from the left, i.e.
!>
!> the distribution matrix for the LDA verion 2 scheme reads:
!> \f[
!> B_i = K_i^+ \left( \sum_{j=1}^{d+1} K_j^+ \right)^{-1}
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
!> \version $Revision: 1.14 $
!> \date $Date: 2020/03/28 09:49:28 $
!> \warning Define \c FSPL_USE_PETSC_KERNELS if you want to use PETSc's kernels to invert matrix \f$ \left( \sum_{j=1}^{d+1} K_j^+ \right) \f$ (the preferred choice)
!> \warning Define \c FSPL_USE_LAPACK if you want to use LAPACK to invert matrix \f$ \left( \sum_{j=1}^{d+1} K_j^+ \right) \f$ (presumably much slower than the former)
!> \bug The source term is NOT distributed in the current release
      SUBROUTINE LDA2SYS_SCHEME(MATRIXSPLITTER,W,NODRES,TSTEP,BETA,
     &                         STIFC,NORDER,NDOF,NOFVERT,VCN,NDIM,
     &                         DFGHDU,LDJ,CELRES,SOURCE,IELEM,
     &                         MATRIX_ASSEMBLY)
C
C     $Id: LDA2Sys_scheme.F,v 1.14 2020/03/28 09:49:28 abonfi Exp $
C     
C     Only one of the following can be defined
C
#define FSPL_USE_PETSC_KERNELS
CCCC#define FSPL_USE_LAPACK
C
C
      IMPLICIT NONE
      INCLUDE 'paramt.h'
      INCLUDE 'constants.h'
      INCLUDE 'flags.com'
C
cnag  EXTERNAL F07ADF,F07AEF,DSCAL
C
C     .. Parameters ..
C     ..
C     .. Scalar Arguments ..
      INTEGER IELEM,LDJ,NDIM,NOFVERT,NDOF,NORDER
      DOUBLE PRECISION VOLUME
      LOGICAL MATRIX_ASSEMBLY
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION DFGHDU(LDJ,*),NODRES(NDOF,NOFVERT),
     &                 TSTEP(NDOF,NOFVERT),SOURCE(*),
     &                 CELRES(*),BETA(*),
     &                 STIFC(NDOF,NDOF,NOFVERT,NOFVERT),
     &                 VCN(NDIM,NOFVERT),W(NDOF,NOFVERT)
C     ..
C     .. Subroutine Arguments ..
      EXTERNAL MATRIXSPLITTER
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION TEMP,HELP
      INTEGER I,J,INFO,IVAR,IVERT,JVERT,IADD,JADD,ORDSQR,N,M
      LOGICAL LFLAG
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION K(MAX_NOFVAR_SQR*VMAX),SUM_K_POS(MAX_NOFVAR_SQR),
     &KNEG(MAX_NOFVAR_SQR*MAXNOFVERT),KPOS(MAX_NOFVAR_SQR*MAXNOFVERT),
     &VLEFT(MAX_NOFVAR_SQR),VRIGHT(MAX_NOFVAR_SQR),
     &WNEG(MAXNOFVAR),WPOS(MAXNOFVAR),WR(MAXNOFVAR),
     &WKSP2(MAXNOFVAR,MAXNOFVAR)
C
#ifdef FSPL_USE_LAPACK
      INTEGER LWORK
      PARAMETER(LWORK=2*MAXNOFVAR)
      INTEGER IPIV(MAXNOFVAR)
      DOUBLE PRECISION WORK(LWORK)
#endif
C
C     ..
C     .. External Functions ..
      LOGICAL  NULLMAT,UNITMAT
      EXTERNAL NULLMAT,UNITMAT
C     ..
C     .. External Subroutines ..
      EXTERNAL DGEMM,DGEMV,DGETRF,DGETRS,DINIT,DSCAL,MATSUM,R8Mat_Print
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC DABS
C     ..
C     .. Statement Functions ..
C     ..
      ORDSQR = NORDER*NORDER
      CALL DINIT(MAX_NOFVAR_SQR,ZERO,SUM_K_POS,1)
C
C residual = - fluctuation
      DO 3 I = 1,2*NDOF
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
C       .. Loop over the vertices of the element ..
C
      DO 10 IVERT = 1,NOFVERT
C
          IADD = (IVERT-1)*ORDSQR+1
C
C       .. The matrix is split into its positive and negative parts
C
C       Numerically or Analitically ..
C
          CALL MATRIXSPLITTER(IELEM,NDIM,NORDER,VCN(1,IVERT),DFGHDU,LDJ,
     +                        K(IADD),KPOS(IADD),
     +                        KNEG(IADD),VLEFT,VRIGHT,NORDER,WR,WPOS,
     +                        WNEG,.TRUE.)
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
C       the residual is computes as Sum_j K(j) * U(j) ..
C
          CALL DGEMV('N',NORDER,NORDER,ONE,K(IADD),NORDER,W(1,IVERT),1,
     +    ONE,CELRES,1) ! should be + source !!
C
   10 CONTINUE
C
C       .. LU factorization of SUM_K_POS ..
C
#ifdef FSPL_USE_LAPACK
C
cnag  CALL F07ADF(NORDER,NORDER,SUM_K_POS,MAXNOFVAR,IPIV,INFO)
      CALL DGETRF(NORDER,NORDER,SUM_K_POS,NORDER,IPIV,INFO)
C
      IF (INFO.GT.0) THEN
          WRITE (6,FMT=99999) IELEM
99999     FORMAT (5X,'Matrix SUM_K_POS is singular in IELEM = ',I6)
          CALL EXIT(INFO)
      ENDIF
C
C     compute the inverse of \sum K^+
C
      CALL DGETRI(NORDER,SUM_K_POS,NORDER,IPIV,WORK,LWORK,INFO)
C
#endif
#ifdef FSPL_USE_PETSC_KERNELS
      CALL SOLVEN( NORDER, SUM_K_POS )
#endif
C
C       .. Loop again over all vertices to compute the distribution
C               matrices (overwritten on KPOS)  ..
C
      DO 30 IVERT = 1,NOFVERT
          IADD = (IVERT-1)*ORDSQR+1
C
C     compute the distribution matrices as B_i =  K^+ (\sum K^+)^{-1}
C
          CALL DGEMM('No Transpose','No Transpose',NORDER,NORDER,NORDER,
     +              ONE,KPOS(IADD),NORDER,SUM_K_POS,NORDER,ZERO,
     +              BETA(IADD),NORDER)
caldo
!             CALL R8Mat_Print('G',' ',NORDER,NORDER,beta(IADD),
!    +                    NORDER,'Beta matrix',INFO)
caldo
C
C       NODRES(IVERT) = -BETA(IVERT)*CELRES
C
         CALL DGEMV('N',NORDER,NORDER,MONE,BETA(IADD),NORDER,
     +   CELRES,1,ZERO,NODRES(1,IVERT),1)
C
C --------------- Debugging code starts here ---------------
C
          IF (ICHECK.NE.0) THEN
C
C       .. Sum up the distribution matrices ..
C
              CALL MATSUM(WKSP2,MAXNOFVAR,BETA(IADD),
     &                    NORDER,NORDER,NORDER)
C
          ENDIF
C
C --------------- Debugging code ends here ---------------
C
C
   30 CONTINUE ! End loop over vertices
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
          JADD = (J-1)*ORDSQR+1
C
          DO 31 I = 1,NOFVERT
          IADD = (I-1)*ORDSQR+1
C
C     ... C_{ij}^{LDA} = - Delta_i^{+} K_j
C         Delta_i^{+} had been overwritten onto KPOS
C
              CALL DGEMM('N','N',NORDER,NORDER,NORDER,MONE,BETA(IADD),
     +             NORDER,K(JADD),NORDER,ZERO,STIFC(1,1,I,J),NDOF)
C
   31 CONTINUE
C
C --------------- Debugging code starts here ---------------
C
      include 'test2.inc'
C
C --------------- Debugging code ends here ---------------
C
      ENDIF
      RETURN

      END
