      SUBROUTINE CDSSYS_SCHEME(MATRIXSPLITTER,W,NodRes,TSTEP,BETA,STIFC,
     +                       NORDER,NDOF,NOFVERT,VCN,NDIM,DFGHDU,LDJ,
     +                       CELRES,SOURCE,IELEM,MATRIX_ASSEMBLY)
C
      IMPLICIT NONE
C
C
C
C NORDER                IN Integer
C                       is the order of the system to solve for, i.e.
C                       the order of the matrix DFGHDU.
C IELEM                 IN Integer
C                       is the current element.
C DFGHDU(LDJ,*) is the Jacobian Matrix of the system.
C LDJ                   IN Integer
C                       is the leading dimension of DFGHDU.
C W(NORDER:NOFVERT)     stores by columns the NORDER variables of the
C                       NOFVERT vertices.
C LNodRes                   IN Integer
C                       is the leading dimension of Q.
C CELRES[1:2*NORDER]  OUT Real
C CELRES[1:NORDER]    stores the residual computed by the Matrix scheme
C                       as \sum K_j U_j (explicit part of the scheme)
C CELRES[NORDER+1:2*NORDER]
C                       stores the residual computed by the Matrix scheme
C                       as \sum C_{ij} U_j (implicit part of the scheme)
C MatrixSplitter        is the procedure used to compute the eigenvector
C                       decomposition of the matrix DFGHDU.
C TSTEP                 is the nodal timestep.
C NodRes(NORDER,NOFVERT)    is the nodal residual.
C
C
C
C
C This routine computes the system N scheme on one tetrahedron
C
C
C
C
C
C
C
C
C
C
C
C
cnag  EXTERNAL F07ADF,F07AEF,DSCAL
C
C
C
C     .. Parameters ..
C
      include 'paramt.h'
      include 'constants.h'
      include 'flags.com'
C
C     ..
C     .. Scalar Arguments ..
      INTEGER IELEM,LDJ,NDIM,NOFVERT,NORDER
      LOGICAL MATRIX_ASSEMBLY
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION DFGHDU(LDJ,*),NodRes(NDOF,NOFVERT),
     +                 TSTEP(NDOF,NOFVERT),CELRES(*),SOURCE(*),
     +                 STIFC(NDOF,NDOF,NOFVERT,NOFVERT),
     +                 VCN(NDIM,NOFVERT),W(NDOF,NOFVERT),BETA(*)
C     ..
C     .. Subroutine Arguments ..
      EXTERNAL MATRIXSPLITTER
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION DUM,TEMP,S
      INTEGER I,INFO,IROW,IVAR,IVERT,J,JCOL,NDOF,IADD,ORDSQR
      LOGICAL LFLAG
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION KMAT(MAXNOFVAR,MAXNOFVAR),
     +KNEG(MAXNOFVAR,MAXNOFVAR,VMAX),KPOS(MAXNOFVAR,MAXNOFVAR,VMAX),
     +SUM_K_NEG(MAXNOFVAR,MAXNOFVAR),UNEG(MAXNOFVAR),
     +VLEFT(MAXNOFVAR,MAXNOFVAR),VRIGHT(MAXNOFVAR,MAXNOFVAR),
     +WKSP1(MAXNOFVAR),WKSP2(MAXNOFVAR),WKSP3(MAXNOFVAR,MAXNOFVAR),
     +WNEG(MAXNOFVAR),WPOS(MAXNOFVAR),WR(MAXNOFVAR)
      INTEGER IPIV(MAXNOFVAR)
C     ..
C     .. External Functions ..
      DOUBLE PRECISION DNRM2
      LOGICAL NULLMAT,UNITMAT
      EXTERNAL DNRM2,NULLMAT,UNITMAT
C     ..
C     .. External Subroutines ..
      EXTERNAL DAXPY,DGEMM,DGEMV,DGETRF,DGETRS,DINIT,DSCAL,MATSUB,
     +         MATSUM,R8Mat_Print
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC DABS
C     ..
      DO 1 J = 1,MAXNOFVAR
          UNEG(J) = ZERO
          DO 1 I = 1,MAXNOFVAR
              SUM_K_NEG(I,J) = ZERO
    1 CONTINUE
      DO 3 I = 1,2*NORDER
          CELRES(I) = ZERO
    3 CONTINUE
C
C --------------- Debugging code starts here ---------------
C
      IF (ICHECK.NE.0) THEN
          CALL DINIT(MAXNOFVAR,ZERO,WKSP1,1)
          CALL DINIT(MAXNOFVAR,ZERO,WKSP2,1)
          CALL DINIT(MAXNOFVAR**2,ZERO,WKSP3,1)
      ENDIF
C
C --------------- Debugging code ends here ---------------
C
      DO 10 IVERT = 1,NOFVERT
C
C
C    The matrix is split into its positive and negative parts
C
C    Numerically or Analitically ..
C
          CALL MATRIXSPLITTER(IELEM,NDIM,NORDER,VCN(1,IVERT),DFGHDU,LDJ,
     +         KMAT,KPOS(1,1,IVERT),KNEG(1,1,IVERT),
     +         VLEFT,VRIGHT,MAXNOFVAR,WR,WPOS,WNEG,.TRUE.)
C
C    The negative jacobians are accumulated in SUM_K_NEG ..
C
          DO 15 JCOL = 1,NORDER
              DO 15 IROW = 1,NORDER
                  SUM_K_NEG(IROW,JCOL) = SUM_K_NEG(IROW,JCOL) +
     +                                   KNEG(IROW,JCOL,IVERT)
   15     CONTINUE
C
C    Computes KNEG(IVERT)*U(IVERT) and adds it to UNEG ..
C
cblas    CALL DGEMV('N',NORDER,NORDER,ONE,KNEG(1,1,IVERT),MAXNOFVAR,W
cblas+   (1,IVERT),1,ONE,UNEG,1)
C
          DO 11 I = 1,NORDER
              DUM = ZERO
              DO 13 J = 1,NORDER
                  DUM = DUM + KNEG(I,J,IVERT)*W(J,IVERT)
   13         CONTINUE
              UNEG(I) = UNEG(I) + DUM
   11     CONTINUE
C
C
C --------------- Debugging code starts here ---------------
C
caldo     IF (ICHECK.NE.0) THEN
C
C    the residual is computes as Sum_j K(j) * U(j) ..
C
              CALL DGEMV('N',NORDER,NORDER,ONE,KMAT,MAXNOFVAR,
     +                   W(1,IVERT),1,ONE,CELRES,1)
caldo     ENDIF
C
C --------------- Debugging code ends here ---------------
C
C    Timestep ..
C
          DO 34 I = 1,NORDER
              TSTEP(1,IVERT) = TSTEP(1,IVERT) + WPOS(I)
   34     CONTINUE
C
   10 CONTINUE
C
C    Finds the generalized inflow point
C       solving [\sum_j K_j^-] U_{-} = \sum_j K_j^- U_j
C
C    LU factorization ..
C
cnag  CALL F07ADF(NORDER,NORDER,SUM_K_NEG,MAXNOFVAR,IPIV,INFO)
      CALL DGETRF(NORDER,NORDER,SUM_K_NEG,MAXNOFVAR,IPIV,INFO)
C
      IF (INFO.GT.0) THEN
          WRITE (6,FMT=99999) IELEM

99999     FORMAT (5X,'Matrix SUM_K_NEG is singular in IELEM = ',I6)
*
          DO 9 IVERT = 1,NOFVERT
              WRITE (6,FMT=*) '    Vertice # ',IVERT
              CALL R8Mat_Print('G',' ',NORDER,NORDER,KPOS(1,1,IVERT),
     +                    MAXNOFVAR,'K(+) ',INFO)
              CALL R8Mat_Print('G',' ',NORDER,NORDER,KNEG(1,1,IVERT),
     +                    MAXNOFVAR,'K(-) ',INFO)
    9     CONTINUE
          CALL R8Mat_Print('G',' ',NORDER,NORDER,SUM_K_NEG,MAXNOFVAR,
     +                'SUM_j K(-) ',INFO)
          STOP

      ENDIF
C
C    solution ..
C
cnag  CALL F07AEF('N',NORDER,1,SUM_K_NEG,MAXNOFVAR,IPIV,UNEG,MAXNOFVAR,INFO)
      CALL DGETRS('N',NORDER,1,SUM_K_NEG,MAXNOFVAR,IPIV,UNEG,
     +            MAXNOFVAR,INFO)
C
C Loops over nodes
C
      S = ONE/REAL(NOFVERT)
C
      DO 30 IVERT = 1,NOFVERT
C
          DO 35 IVAR = 1,NORDER
              W(IVAR,IVERT) = W(IVAR,IVERT) - UNEG(IVAR)
   35     CONTINUE
C
C    NodRes[INODE] = - KPOS(IVERT)*[U(IVERT)-UNEG]
C
cblas    CALL DGEMV('N',NORDER,NORDER,ONE,KPOS(1,1,IVERT),MAXNOFVAR, W
cblas+   (1,IVERT),1,ZERO,WKSP,1)
C
          DO 23 IVAR = 1,NORDER
              NodRes(IVAR,IVERT) = -(CelRes(ivar)+SOURCE(ivar))*S
   23     CONTINUE
C
C --------------- Debugging code starts here ---------------
C
          IF (ICHECK.NE.0) THEN
C
C    The residual is computed as Sum_j KPOS(j) * [U(j)-Uin]
C       and stored in WKSP1
C
              CALL DAXPY(NORDER,-ONE,NodRes(1,IVERT),1,WKSP1,1)
          ENDIF
C
C --------------- Debugging code ends here ---------------
C
   30 CONTINUE
C
C --------------- Debugging code starts here ---------------
C
      include 'test1.inc'
C
C --------------- Debugging code ends here ---------------
C
      IF (.NOT.MATRIX_ASSEMBLY) RETURN
C
C     Assembling the element stiffness matrix for the N scheme ..
C
      DO 31 J = 1,NOFVERT
C
C     Solve [\sum_K^{-}] Delta_j^{-} = K_j^{-}
C         Delta_j^{-} is overwritten onto K_j^{-}
C
          CALL DGETRS('N',NORDER,NORDER,SUM_K_NEG,MAXNOFVAR,IPIV,
     +                KNEG(1,1,J),MAXNOFVAR,INFO)
C
          DO 31 I = 1,NOFVERT
C
              IF (I.EQ.J) GOTO 31
C
C     C_{ij}^N = K_i^{+} Delta_j^{-} i neq j
C
              CALL DGEMM('N','N',NORDER,NORDER,NORDER,ONE,KPOS(1,1,I),
     +             MAXNOFVAR,KNEG(1,1,J),MAXNOFVAR,ZERO,STIFC(1,1,I,J),
     +             NDOF)
C
C     C_{ii} = - \sum_{j \neq i} C_{ij}
C                = - K_i^{+} [ \sum_{j \neq i} \Delta_j^- ] ...
C
              CALL MATSUB(STIFC(1,1,I,I),NDOF,STIFC(1,1,I,J),NDOF,
     +                    NORDER,NORDER)
C
   31 CONTINUE
C
C --------------- Debugging code starts here ---------------
C
      IF (ICHECK.NE.0) THEN
C
C        check that \sum_j \Delta_j^- = Identity matrix
C
          CALL DINIT(MAXNOFVAR**2,ZERO,WKSP3,1)
          DO 14 J = 1,NOFVERT
              CALL MATSUM(WKSP3,MAXNOFVAR,KNEG(1,1,J),MAXNOFVAR,
     +        NORDER,NORDER)
   14     CONTINUE
          LFLAG = UNITMAT(WKSP3,NORDER,NORDER,MAXNOFVAR,1.D-14)
          IF (.NOT.LFLAG) THEN
              WRITE (6,FMT=*) IELEM,I
              CALL R8Mat_Print('General',' ',NORDER,NORDER,WKSP3,
     +             MAXNOFVAR,' Sum Delta_j^- equals identity ?',INFO)
              PAUSE

          ENDIF
      ENDIF
C
      include 'test2.inc'
C
C --------------- Debugging code ends here ---------------
C
      RETURN

      END
