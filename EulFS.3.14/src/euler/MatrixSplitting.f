C
      SUBROUTINE MatSplitNum(IELEM,NDIM,NORDER,VNOR,DUMMYA,DUMMYB,
     +JacobianMatrix,LDJ,Cmat,Cpos,Cneg,VLEFT,VRIGHT,LD,
     +WR,LPOS,LNEG,EIGENDECO)
C
C************************************************************
C
C     IELEM:    is the current element (INPUT)
C     NORDER:   is the order of the system to be solved (INPUT)
C     VNOR:     is the face normal (INPUT)
C JacobianMatrix:  is the jacobian matrix times the face normal (OUTPUT)
C     Cpos:     is the "positive" part of the jacobian matrix (OUTPUT)
C     Cneg:     is the "negative" part of the jacobian matrix (OUTPUT)
C     VLEFT:    is left eigenvector matrix of the jacobian matrix (OUTPUT)
C     VRIGHT:   is right eigenvector matrix of the jacobian matrix (OUTPUT)
C     LD:       is the leading dimension of the previous matrices
C     WR:       eigenvalues of the jacobian matrix (OUTPUT)
C     LPOS:     positive eigenvalues of the jacobian matrix (OUTPUT)
C     LNEG:     negative eigenvalues of the jacobian matrix (OUTPUT)
C
C     EIGENDECO .TRUE. if the eigenvector decomposition of the matrix
C               is needed, if .FALSE. only the K matrix is computed (INPUT)
C
C************************************************************
C
      IMPLICIT NONE
C
C     .. Parameters ..
C
      INCLUDE 'paramt.h'
C
C     .. Commons ..
C
      INCLUDE 'flags.com'
C
C     .. Parameters ..
C
      INTEGER    LWORK,LDH,LDQ,LDA
      PARAMETER(LWORK=64*NMAX,LDH=NMAX,LDQ=NMAX,LDA=NMAX)
      DOUBLE PRECISION ZERO,ONE
      PARAMETER(ZERO=0.D0,ONE=1.D0)
C
C     .. Scalar Arguments ..
C
      INTEGER IELEM,NDIM,NORDER,LD,LDJ
      LOGICAL EIGENDECO
C
C     .. Array Arguments ..
C
      DOUBLE PRECISION VNOR(1),JacobianMatrix(LDJ,1), Cmat(LD,1), Cpos
     +(LD,1),Cneg(LD,1), VLEFT(LD,1),VRIGHT(LD,1),WR(1),LPOS(1),LNEG(1),
     +DUMMYA(*),DUMMYB(*)
C
C     .. Local Scalars ..
C
      INTEGER INFO,I,J,M,IFAIL,IROW,JCOL,IDIM,JDIM
      DOUBLE PRECISION TEMP
      LOGICAL LFLAG
C
C     .. Local Arrays ..
C
      DOUBLE PRECISION A
      DIMENSION    A(LDA,NMAX)
      INTEGER   IPIV
      DIMENSION IPIV(NMAX)
      INTEGER   IWORK(LWORK)
C
      DOUBLE PRECISION WORK,WI,WORK1,WORK2
      DIMENSION WORK(LWORK),WI(NMAX), WORK1(NMAX,NMAX),WORK2(NMAX,NMAX)
C
      DOUBLE PRECISION WKSP2
      DIMENSION WKSP2(NMAX,NMAX)
C
      DATA LFLAG / .FALSE. /
C
C     .. External Subroutines ..
C
cnag  EXTERNAL    F07ADF,F07AJF,R8Mat_Print
      EXTERNAL    DGETRF,DGETRI,R8Mat_Print
C
C     .. External Functions ..
C
C
C     .. Intrinsic Functions ..
C
      INTRINSIC    DABS
C
C     .. Data Statements ..
C
C     .. Executable Statements ..
C
C     .. Initializing the matrix K ..
C
      CALL DSCAL(LD*LD,ZERO,Cmat,1)
C
C     .. Assembling the matrix K = ( JacobianMatrix . Normal ) / d ..
C
      DO 32 IDIM =  1, NDIM
         JDIM = (IDIM-1)*LDJ
         TEMP = VNOR(IDIM)/NDIM
         DO 22 IROW =  1, NORDER
            DO 22 JCOL = 1 , NORDER
               Cmat(IROW,JCOL) = Cmat(IROW,JCOL) +
     &         JacobianMatrix(IROW,JDIM+JCOL)*TEMP
C
C    .. The original matrix is saved in A ..
C
               A(IROW,JCOL) = Cmat(IROW,JCOL)
   22    CONTINUE
   32 CONTINUE
C
      IF( EIGENDECO .EQV. .FALSE. )RETURN
C
      INFO = 0
*
*    Calculate Right eigenvectors of A
*
      CALL rg(LDA,NORDER,A,wr,wi,1,VRIGHT,IWORK,WORK,INFO)
c
c     .. The eigenvectors are divided by STAGFIX ..
c
*     CALL DSCAL(NORDER,ONE/STAGFIX,WR,1)
C
C    .. VRIGHT is copied into VLEFT
C
      CALL DCOPY(LD*LD,VRIGHT,1,VLEFT,1)
cnag  CALL F06QFF('General',NORDER,NORDER,VRIGHT,LD,VLEFT,LD)
C
C    .. VLEFT is factorized ..
C
cnag  CALL F07ADF(NORDER,NORDER,VLEFT,LD,IPIV,INFO)
      CALL DGETRF(NORDER,NORDER,VLEFT,LD,IPIV,INFO)
C
C    .. VLEFT is inverted ..
C
cnag  CALL F07AJF(NORDER,VLEFT,LD,IPIV,WORK,LWORK,INFO)
      CALL DGETRI(NORDER,VLEFT,LD,IPIV,WORK,LWORK,INFO)
C
      DO 10 I = 1 , NORDER
         LPOS(I) = 0.5D0 * (WR(I) + DABS(WR(I)) )
         LNEG(I) = 0.5D0 * (WR(I) - DABS(WR(I)) )
   10 CONTINUE
C
C    .. Computes Lambda(+/-) * VLEFT
C
      DO 20 I = 1 , NORDER
         DO 20 J =  1, NORDER
            WORK1(I,J) = LPOS(I) * VLEFT(I,J)
            WORK2(I,J) = LNEG(I) * VLEFT(I,J)
   20 CONTINUE
C
C    .. Computes VRIGHT * Lambda(+/-) * VLEFT
C
      CALL DGEMM('N','N',NORDER,NORDER,NORDER,1.D0,VRIGHT,LD, WORK1,
     +NMAX,0.D0,Cpos,LD)
      CALL DGEMM('N','N',NORDER,NORDER,NORDER,1.D0,VRIGHT,LD, WORK2,
     +NMAX,0.D0,Cneg,LD)
C
      IF(ICHECK.EQ.2)THEN
         DO 30 I =  1, NORDER
            DO 30 J =  1, NORDER
               WKSP2(I,J) = Cpos(I,J)+Cneg(I,J)
               WORK2(I,J) = WKSP2(I,J)-Cmat(I,J)
               IF( DABS(WORK2(I,J)) .GT. 1.D- 14)LFLAG = .TRUE.
   30    CONTINUE
C
C        lflag = .true. 
         IF( LFLAG .EQV. .TRUE. )THEN
            LFLAG = .FALSE.
C
            IF(ICHECK.EQ.2)THEN
               WRITE(6,99998)IELEM
               WRITE(6,99999)(' (',WR(I),',',WI(I),')',I=1,NORDER)
99998 FORMAT(5X,'Element # ',I6)
99999 FORMAT(1X,A,F8.4,A,F8.4,A)
C
               CALL R8Mat_Print('General',' ',NORDER,NORDER,VRIGHT,LD,
     +         'Right eigenvector matrix',IFAIL)
               CALL R8Mat_Print('General',' ',NORDER,NORDER,VLEFT,LD,
     +         'Left eigenvector matrix',IFAIL)
            ENDIF
C
            CALL R8Mat_Print('General',' ',NORDER,NORDER,Cmat,LD,
     +      'Original matrix',IFAIL)
            CALL R8Mat_Print('General',' ',NORDER,NORDER,WKSP2,NMAX,
     +      'Reassembled matrix',IFAIL)
            CALL R8Mat_Print('General',' ',NORDER,NORDER,WORK2,NMAX,
     +      'Error matrix',IFAIL)
            PAUSE
         ENDIF
      ENDIF
      RETURN
      END
