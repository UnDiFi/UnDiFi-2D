!> \copydetails LDASYS_SCHEME()
      SUBROUTINE FVSys_scheme(MATRIXSPLITTER,W,NODRES,TSTEP,BETA,STIFC,
     &NORDER,NDOF,NOFVERT,VCN,NDIM,DFGHDU,LDJ,CELRES,SOURCE,IELEM,
     &MATRIX_ASSEMBLY)
C
      IMPLICIT NONE
C
C
!> This routine computes the FV scheme for systems on one tetrahedron
!! or triangle
C
      include 'paramt.h'
      include 'constants.h'
      include 'flags.com'
      include 'time.com'
C
C     .. External Arguments ..
C
      EXTERNAL  MATRIXSPLITTER
C
C     .. Scalar Arguments ..
C
      INTEGER NORDER,FrstEq,IELEM,LDJ,LDQ,NDOF,NDIM,NOFVERT
C
C     .. Array Arguments ..
C
      DOUBLE PRECISION DFGHDU(LDJ,*),W(NDOF,NOFVERT),CELRES(*),
     &TSTEP(NDOF,NOFVERT),NODRES(NDOF,NOFVERT),SOURCE(*),
     &VCN(NDIM,NOFVERT),BETA(*),STIFC(NDOF,NDOF,NOFVERT,NOFVERT)
C
C     .. Local Scalars ..
C
      INTEGER INODE,IVERT,IVAR,IDIM,IROW,I,J
      DOUBLE PRECISION TEMP1,DENOM
      LOGICAL MATRIX_ASSEMBLY
C
C     .. Local Arrays ..
C
      DOUBLE PRECISION KPOS(MAXNOFEQN,MAXNOFEQN),
     &KNEG(MAXNOFEQN,MAXNOFEQN),KMAT(MAXNOFEQN,MAXNOFEQN,VMAX),
     &J_dot_N(MAXNOFEQN,MAXNOFEQN),
     &DQ(MAXNOFEQN),WKSP4(MAXNOFEQN),WKSP5(MAXNOFEQN),
     &VLEFT(MAXNOFEQN,MAXNOFEQN),VRIGHT(MAXNOFEQN,MAXNOFEQN),
     &VN(3),WR(MAXNOFEQN),WPOS(MAXNOFEQN),WNEG(MAXNOFEQN)
C
C     .. External Functions ..
C
C
C     .. External Subroutines ..
C
      EXTERNAL DSCAL
C
C     .. Intrinsic Functions ..
C
      INTRINSIC DBLE
C
C     .. Data Statements ..
C
      DATA VN /3*ZERO/
C
C     .. Executable Statements ..
C
      DENOM = ONE/REAL(NOFVERT)
      CALL DSCAL(NORDER,ZERO,CELRES,1)! residual = - fluctuation
C
      IF(ICHECK.NE.0)THEN
C
         CALL DSCAL(MAXNOFEQN,ZERO,WKSP5,1)
         CALL DINIT(MAXNOFEQN*MAXNOFEQN*VMAX,ZERO,KMAT,1)
C
C Compute the K's
C
         DO 10 IVERT =  1, NOFVERT
C
C.. The jacobian matrix is dotted with the face normal
C   and saved in KMAT; the fluctuation is computed
C
C
cold        CALL MATRIXSPLITTER(IELEM,NORDER,VN,DFGHDU,LDJ,
cold +      K(1,1,IVERT),KPOS(1,1),KNEG(1,1),VLEFT,VRIGHT,MAXNOFEQN,WR,
cold +      WPOS,WNEG,.FALSE.)
            CALL MATRIXSPLITTER(IELEM,NDIM,NORDER,VCN(1,IVERT),DFGHDU,
     +      LDJ,KMAT(1,1,IVERT),KPOS(1,1),KNEG(1,1),VLEFT,VRIGHT,
     +      MAXNOFEQN,WR,WPOS,WNEG,.FALSE.)
C
            CALL DGEMV('N',NORDER,NORDER,ONE,KMAT(1,1,IVERT),MAXNOFEQN,
     +      W(1,IVERT),1,ONE,CELRES,1)
C
   10    CONTINUE ! End loop over vertices
      ENDIF ! ICHECK
C
C.. Loop over the edges ..
C
      DO 15 I = 1 , NOFVERT
         DO 14 J = 1, NOFVERT
            IF(J.EQ.I)GOTO 14
c
c.. Compute [n(j)-n(i)]/d
c
            DO 17 IDIM =  1, NDIM
   17       VN(IDIM) = VCN(IDIM,J) - VCN(IDIM,I)
C
            DO 19 IROW =  1, NORDER
   19       DQ(IROW) = W(IROW,J) - W(IROW,I)
C
C
C       .. The matrix is split into its positive and negative parts
C
C       .. Numerically or Analitically ..
C
            CALL MATRIXSPLITTER(IELEM,NDIM,NORDER,VN,DFGHDU,LDJ,
     +      J_dot_N,KPOS(1,1),KNEG(1,1),VLEFT,VRIGHT,MAXNOFEQN,WR,WPOS,
     +      WNEG,.TRUE.)
 
C
C *************** Contribution to node J ***************
C
C.. NODRES[J] := NODRES[J] - [K(JI)]^+*[U(J)-U(I)]/(d+1)
C
C.. WKSP4 := [KPOS(J)-KPOS(I)]*[U(J)-U(I)]/(d+1)
C
            CALL DGEMV('N',NORDER,NORDER,-DENOM,KPOS,MAXNOFEQN,
     +      DQ, 1,ZERO,NODRES(1,J),1)
C
            IF(ICHECK.NE.0)CALL DAXPY(NORDER,-ONE,NODRES(1,J),1,WKSP5,1)
C
            DO 7 IVAR =  1, NORDER
    7       TSTEP(IVAR,J) = TSTEP(IVAR,J) + WPOS(IVAR)
C
C
C
   14 CONTINUE ! End loop over edges
   15 CONTINUE ! End loop over edges
C
      IF(ICHECK.NE.0)THEN
C
         DO IVAR =  1, NORDER
            TEMP1 = WKSP5(IVAR) - CELRES(IVAR)
            IF(DABS(TEMP1).GT.1.D-15)THEN
               WRITE(6,*)'FV system scheme, elem ',IELEM,' var # ',IVAR,
     +         ' computed ',WKSP5(IVAR),' "true" ',CELRES(IVAR)
            ENDIF
         ENDDO
      ENDIF
      IF(MATRIX_ASSEMBLY)STOP ' Picard iteration NOT implemented in FV'
      IF(LTIME)STOP ' Mass-matrix NOT implemented in FV'
C
      RETURN
      END
