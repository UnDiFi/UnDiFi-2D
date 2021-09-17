!> \brief \b UCVSYS_SCHEME
!> \par Purpose
!>
!> This routine computes the UCV system scheme on one tetrahedron
!> the UCV scheme (Upwind Control Volume) is credited to Mike Giles and described in \c doi:10.1016/j.jcp.2004.11.023
!>
!> The distribution matrix for the UCV scheme reads
!> \f[
!> B_i = \frac{1}{d+1} \left[ I_{m\times m} + d \, \left( \sum_{j=1}^{d+1} |K_j|\right)^{-1} K_i \right]
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
!> \version $Revision: 1.4 $
!> \date $Date: 2020/03/28 09:49:28 $
!> \warning Define \c FSPL_USE_PETSC_KERNELS if you want to use PETSc's kernels to invert matrix \f$ \left( \sum_{j=1}^{d+1} |K_j|\right) \f$ (the preferred choice)
!> \warning Define \c FSPL_USE_LAPACK if you want to use LAPACK to invert matrix \f$ \left( \sum_{j=1}^{d+1} |K_j|\right) \f$ (presumably much slower than the former)
!> \bug The source term is NOT distributed in the current release
      SUBROUTINE UCVSYS_SCHEME(MATRIXSPLITTER,W,NODRES,TSTEP,BETA,STIFC,
     +                         NORDER,NDOF,NOFVERT,VCN,NDIM,DFGHDU,
     +                         LDJ,CELRES,SOURCE,IELEM,
     &                         MATRIX_ASSEMBLY)
C
C     $Id: UCVSys_scheme.F,v 1.4 2020/03/28 09:49:28 abonfi Exp $
C
C
C     Only one of the following can be defined:
C
#define FSPL_USE_PETSC_KERNELS
CCCC#define FSPL_USE_LAPACK
C
C
C
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
      DOUBLE PRECISION CFLELEM
C
#ifdef FSPL_USE_LAPACK
      INTEGER LWORK
      PARAMETER(LWORK=2*MAXNOFVAR)
      INTEGER IPIV(LWORK)
      DOUBLE PRECISION WORK(LWORK)
#endif
C
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION VOLUME
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
      DOUBLE PRECISION GALERKIN,TEMP,HELP
      INTEGER I,INFO,IROW,IVAR,IVERT,JVERT,J,N,M,ORDSQR,IADD,JADD
      LOGICAL LFLAG
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION SUM_K_ABS(MAX_NOFVAR_SQR),
     + KNEG(MAX_NOFVAR_SQR*MAXNOFVERT),KPOS(MAX_NOFVAR_SQR*MAXNOFVERT),
     + K(MAX_NOFVAR_SQR*MAXNOFVERT),
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
      CFLELEM = REAL(NDIM)/(NDIM+ONE)
C
      CALL DINIT(ORDSQR,ZERO,SUM_K_ABS,1)
C
C     Loop over the vertices of the element ..
C
      DO 10 IVERT = 1,NOFVERT
C
C     The matrix IS split into its positive and negative parts
C
          IADD = ORDSQR*(IVERT-1)+1
          CALL MATRIXSPLITTER(IELEM,NDIM,NORDER,VCN(1,IVERT),DFGHDU,LDJ,
     +                        K(IADD),KPOS(IADD),KNEG(IADD),
     &                        VLEFT,VRIGHT,NORDER,
     +                        WR,WPOS,WNEG,.TRUE.)
C
C       The absolute value jacobians are accumulated in SUM_K_ABS ..
C
          CALL DAXPY(ORDSQR, ONE,KPOS(IADD),1,SUM_K_ABS,1)
          CALL DAXPY(ORDSQR,MONE,KNEG(IADD),1,SUM_K_ABS,1)
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
#ifdef FSPL_USE_LAPACK
C
C       .. LU factorization of SUM_K_ABS ..
C
cnag  CALL F07ADF(NORDER,NORDER,SUM_K_ABS,NORDER,IPIV,INFO)
      CALL DGETRF(NORDER,NORDER,SUM_K_ABS,NORDER,IPIV,INFO)
!         CALL R8Mat_Print('General',' ',NORDER,NORDER,SUM_K_ABS,
!    +                 NORDER,' sum |K| matrix ',INFO)
!         pause
C
      IF (INFO.GT.0) THEN
          WRITE (6,FMT=99999) IELEM
99999     FORMAT (5X,'Matrix SUM_K_ABS is singular in IELEM = ',I6)
          CALL R8Mat_Print('General',' ',NORDER,NORDER,SUM_K_ABS,
     +                 NORDER,' sum |K| matrix ',INFO)
      DO 19 IVERT = 1,NOFVERT
          IADD = ORDSQR*(IVERT-1)+1
          CALL R8Mat_Print('General',' ',NORDER,NORDER,K(IADD),
     +                 NORDER,' sum |K| matrix ',INFO)
   19 CONTINUE
          STOP
      ENDIF
C
C     compute the inverse of (\sum K^+)
C
      CALL DGETRI(NORDER,SUM_K_ABS,NORDER,IPIV,WORK,LWORK,INFO)
C
      IF (INFO.GT.0) THEN
          WRITE (6,FMT=99998) INFO
99998     FORMAT (5X,'DGETRI HAS RETURNED INFO = ',I6)
          STOP
      ENDIF
C
#endif
#ifdef FSPL_USE_PETSC_KERNELS
      CALL SOLVEN( NORDER, SUM_K_ABS )
#endif
C       .. Loop again over all vertices to compute the distribution
C               matrices (overwritten on KPOS)  ..
C
      DO 30 IVERT = 1,NOFVERT
C
C     compute the distribution matrices as B_i =  (\sum K^+)^{-1} K_i^+ 
C
         IADD = ORDSQR*(IVERT-1)+1
         CALL DGEMM('No Transpose','No Transpose',NORDER,NORDER,NORDER,
     +              CFLELEM,SUM_K_ABS,NORDER,K(IADD),NORDER,
     +              ZERO,BETA(IADD),NORDER)
C
C        .. add the Galerkin term
C
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
     +   CELRES,1,ZERO,NODRES(1,IVERT),1) ! should be + SOURCE
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
C
C --------------- Time accurate stuff starts here ---------------
C
      RETURN

      END
