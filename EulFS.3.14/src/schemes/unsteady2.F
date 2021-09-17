      SUBROUTINE UNSTEADY2(WORK,BETA,Z,NOFVAR,NODRES,STIFC,NDOF,NDIM,
     &                     NOFVERT,MATRIX_ASSEMBLY)
C
C     $Id: unsteady2.F,v 1.2 2020/03/28 09:49:28 abonfi Exp $
C
      IMPLICIT NONE
C
C     compute the mass matrix for the incompressible equations
C
      INCLUDE 'paramt.h'
      INCLUDE 'constants.h'
      INCLUDE 'time.h'
      INCLUDE 'time.com'
C
C     .. Parameters ..
C     ..
C     .. Scalar Arguments ..
      INTEGER NOFVAR,NDOF,NDIM,NOFVERT
      LOGICAL MATRIX_ASSEMBLY
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION BETA(NDOF,NDOF,NOFVERT),Z(NOFVAR,NOFVERT,*),
     +       NODRES(NOFVAR,NOFVERT),WORK(*),
     2       STIFC(NDOF,NDOF,NOFVERT,NOFVERT)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION RD,S,HELP,ALPHA,OFFD,DIAG
      INTEGER I,J,K,L,IADDR,ORDSQR,INFO,JADDR
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION MMAT(MAX_NOFVERT_SQR*MAX_NOFVAR_SQR),
     &                 WKSP(MAXNOFVAR),
     2                 MATHELP(MAX_NOFVAR_SQR*MAXTIMLEVS)
#ifdef DEBUG
      DOUBLE PRECISION MATXXX(MAX_NOFVAR_SQR),MATYYY(MAX_NOFVAR_SQR),
     2       ROWSUM(MAX_NOFVAR_SQR),ROWTRU(MAX_NOFVAR_SQR*MAXNOFVERT)
#endif
C     ..
C     .. External Functions ..
C     ..
C     .. Statement Functions ..
      INTEGER IDX
      IDX(I,J) = ((J-1)*NOFVERT+I-1)*ORDSQR+1
C
      ORDSQR = NDOF*NDOF
!     HELP = (ONE+HALF*GAMT)/(NOFVERT*DTVOL)
      RD = ONE/REAL(NOFVERT)
      GOTO (10,20,30,40,50) MMTYPE
   10 CONTINUE ! lumped mass matrix
      STOP 'MM-LM yet to be implemented'
      GOTO 100
   20 CONTINUE ! Petrov-Galerkin
      CALL DINIT(ORDSQR*NOFVERT*NOFVERT,ZERO,MMAT,1)! set mass matrix to zero
      DIAG = (12.d0-NDIM)/60.d0
      OFFD = (TWO*NDIM-9.d0)/60.d0
      IADDR = 0
      DO 3 J = 1,NOFVERT
         DO 3 I = 1,NOFVERT
            IADDR = IDX(I,J)
            CALL DCOPY(ORDSQR,BETA(1,1,I),1,MMAT(IADDR),1)
            IADDR = IADDR -1
            IF(J.EQ.I)THEN
                DO K = 1,ORDSQR,(NDOF+1)
                    MMAT(IADDR+K) = MMAT(IADDR+K)+DIAG
                ENDDO
            ELSE
                DO K = 1,ORDSQR,(NDOF+1)
                    MMAT(IADDR+K) = MMAT(IADDR+K)+OFFD
                ENDDO
            ENDIF
    3 CONTINUE
      GOTO 100
   30 CONTINUE ! Consistent Upwind
      STOP 'MM-CU yet to be implemented'
      DO 1 J = 1,NOFVERT
         DO 1 I = 1,NOFVERT
            IADDR = IDX(I,J)
C     compute -B_i B_j
            CALL DGEMM('No','No',NDOF,NDOF,NDOF,MONE,BETA(1,1,I),NDOF,
     &                 BETA(1,1,J),NDOF,ZERO,MMAT(IADDR),NDOF)
    1 CONTINUE
      DO 2 J = 1,NOFVERT
         DO 2 I = 1,NOFVERT
            IADDR = IDX(I,J)
            IF(J.EQ.I)THEN
               ALPHA = REAL(NDIM) ! not sure it is ok for 3D 
            ELSE
               ALPHA = ONE 
            ENDIF
            CALL DAXPY(ORDSQR,ALPHA,BETA(1,1,I),1,MMAT(IADDR),1)
    2 CONTINUE
      GOTO 100
   40 CONTINUE ! Simple Upwind
      STOP 'MM-SU yet to be implemented'
      IADDR = 0
      DO 7 J = 1,NOFVERT
         DO 7 I = 1,NOFVERT
            IADDR = IDX(I,J)
            CALL DCOPY(ORDSQR,BETA(1,1,I),1,MMAT(IADDR),1)
    7 CONTINUE
      GOTO 100
   50 CONTINUE ! Centred
      STOP 'MM-CC yet to be implemented'
      CALL DINIT(NOFVERT*NOFVERT*ORDSQR,ZERO,MMAT,1)
      IADDR = 0
      DO 8 J = 1,NOFVERT
         DO 8 I = 1,NOFVERT
            IADDR = IDX(I,J) -1
            DO K = 1,ORDSQR,(NDOF+1)
                MMAT(IADDR+K) = RD
         ENDDO
    8 CONTINUE
      GOTO 100
  100 CONTINUE
!        DO 9 J = 1,NOFVERT
!           DO 9 I = 1,NOFVERT
!              IADDR = ((J-1)*NOFVERT+I)*ORDSQR+1
!              IADDR = IDX(I,J)
!              CALL DINIT(NDOF,ZERO,MMAT(IADDR),NDOF)
!   9    CONTINUE
!        CALL DINIT(NDOF*NOFVERT*NOFVERT,ZERO,MMAT(1),NDOF)
      call zerofirstrowcol(mmat,ndof,nofvert)
c
!     DO J = 1,NOFVERT
!     DO I = 1,NOFVERT
!        IADDR = IDX(I,J)
!        WRITE(6,*)'I = ',I,' J = ',J 
!        CALL R8Mat_Print('G',' ',NDOF,NDOF,MMAT(IADDR),
!    +                NDOF,'Mass matrix ',INFO)
!         CALL R8Mat_Print('G',' ',NOFVERT,NOFVERT,STIFC(1,1),
!    +                NOFVERT,'C_ij matrix ',INFO)
!     ENDDO
!     ENDDO
C
      IF(DUALTS)THEN ! dual time stepping
C
C update the rhs by adding the contribution
C of the time derivative term
C
C        compute (u^n-u^{n-1}) .AND.  (u^{n+1}-u^{n})
C
!     CALL R8Mat_Print('G',' ',NOFVERT,3,Z(1,1),
!    +                NOFVERT,'Z matrix ',INFO)
C
caldo
!     CALL R8Mat_Print('G',' ',ndof,NOFVERT*NTIMLEVS,Z,
!    +                ndof,'Z matrix ',INFO)
!     CALL R8Mat_Print('G',' ',ndof,NOFVERT*NTIMLEVS,Z,
!    +                ndof,'Z matrix before  ',INFO)
!     CALL R8Mat_Print('G',' ',ndof,NOFVERT*NTIMLEVS,NODRES,
!    +                ndof,'R matrix before  ',INFO)
caldo
!     write(6,*)'coeff = ',(a(k),k=-1,1)
!        CALL DINIT(NDOF,ZERO,WORK,1)
!        CALL DINIT(NDOF,ZERO,WKSP,1)
!        CALL DINIT(ORDSQR*MAXNOFVERT,ZERO,MATHELP,1)
C
C     set to zero R^{n} and R^{n-1}
C
!        CALL DINIT(NDOF*NOFVERT*(NTIMLEVS-1),ZERO,NODRES(1,1,2),1)
C
         HELP = ONE/(NOFVERT*DTVOL)
!        HELP = ONE/DTVOL
         DO K = 1,NTIMLEVS
            ALPHA =-HELP*TCOEF(2-K)
            DO I = 1,NOFVERT
!              CALL DAXPY(NDOF,ALPHA,Z(1,I,K),1,WKSP,1)
               DO J = 1,NOFVERT
                  IADDR = IDX(I,J)
                  JADDR = ((K-1)*NOFVERT+(J-1))*NOFVAR+1
                  CALL DGEMV('No Transpose',NDOF,NDOF,ALPHA,MMAT(IADDR),
     2                    NDOF,Z(1,J,K),1,ONE,NODRES(1,I),1)
!                 CALL DGEMV('No Transpose',NDOF,NDOF,ALPHA,MMAT(IADDR),
!    2                    NDOF,Z(JADDR),1,ONE,NODRES(1,I,K),1)
               ENDDO ! j
            ENDDO ! i
         ENDDO ! k
!     CALL R8Mat_Print('G',' ',ndof,NOFVERT*NTIMLEVS,NODRES,
!    +                ndof,'R matrix  ',INFO)
!     write(6,*)'ielem = ',-1,MATRIX_ASSEMBLY
!        do j = 1,nofvert 
!           k = (j-1)*ordsqr+1 
!           do i = 1,nofvert 
!              IADDR = IDX(I,J)
!              CALL daxpy(NDOF*NDOF,ONE,MMAT(IADDR),1,MATHELP(k),1)
!           enddo
!        enddo
!     CALL R8Mat_Print('G',' ',ndof,NOFVERT*Ndof,bighelp,
!    +                ndof,'I matrix ????  ',INFO)
  110 FORMAT(4(E12.6,1X))
!     write(6,*)'ielem = ',ielem
C
C update the matrix by adding the mass matrix
C
         IF(MATRIX_ASSEMBLY)THEN
            HELP =-(ONE+HALF*GAMT)/(NOFVERT*DTVOL)
            DO J = 1,NOFVERT
               DO I = 1,NOFVERT
                  IADDR = IDX(I,J)
                  CALL DAXPY(ORDSQR,HELP,MMAT(IADDR),1,STIFC(1,1,I,J),1)
               ENDDO
            ENDDO
         ENDIF
      ELSE ! no inner iterations (DUALTS == .FALSE.)
         STOP 'Unimplemented feature: Must use dual ts !'
      ENDIF 
      RETURN
      END
      subroutine zerofirstrowcol(a,ndof,nofvert)
      implicit none
      integer ndof,nofvert
      double precision a(ndof,ndof,nofvert,nofvert)
      integer i,j,k,l
      do k = 1,nofvert
         do l = 1,nofvert
c set to 0.d0 rows and cols corresponding to pressure
            do j = 1,ndof
               a(1,j,k,l) = 0.d0
            enddo
            do i = 2,ndof
               a(i,1,k,l) = 0.d0
            enddo
         enddo
      enddo
      return
      end
 
