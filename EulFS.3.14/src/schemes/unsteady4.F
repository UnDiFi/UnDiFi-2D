!> \par Purpose
!>
!> Add the time-dependent term for unsteady Euler problems
!>
!> If we deal with a non-linear PDE, a dual time-stepping strategy is adopted.
!> We can define a modified nodal residual (\f$R\f$) augmented with the time derivative term:
!> \f[
!> R_i\left(u^{n+1,k}\right) := R_i\left(u^{n+1,k}\right) -
!> \sum_{j} m_{ij} \frac{\left(1+\gamma_t/2\right)u_j^{n+1}-\left(1+\gamma_t\right)u_j^n
!> +\left(\gamma_t/2\right)u_j^{n-1}}{\Delta t}
!> \f]
!>
!> When using an ALE formulation the nodal residual should be updated as follows, see J. Dobes phD thesis on page 62, Eq.(3.119):
!>
!> \f[
!> R_i\left(u^{n+1,k}\right) := R_i\left(u^{n+1,k}\right) -
!> \sum_{j} m_{ij} \frac{\left(1+\gamma_t/2\right)|T^e|^{n+1}u_j^{n+1}-\left(1+\gamma_t\right)|T^e|^{n}u_j^n
!> +\left(\gamma_t/2\right)|T^e|^{n-1}u_j^{n-1}}{\Delta t} +
!> \frac{\left(1+\gamma_t/2\right)|T^e|^{n+1}-\left(1+\gamma_t\right)|T^e|^{n} +\left(\gamma_t/2\right)|T^e|^{n-1}}{\Delta t} 
!> \sum_{j} m_{ij} u_j^{n+1,k}
!> \f]
!>
!>
!> various kinds of mass-matrices could be assembled, depending on the value of \c MYTYPE and according to the table below:
!>
!> at present, however, only the mass-matrix described in [here] is available.
!> [here]: http://www.tandfonline.com/doi/abs/10.1080/10618562.2013.813491#.UuuQcf6nSnw "A. Bonfiglioli, R. Paciorri. A mass-matrix formulation of unsteady fluctuation splitting schemes consistent with Roe’s parameter vector, 2014"
!>
!
!> @param[in] dUdZ is \f$\partial U/\partial Z\f$
!> @param[in] BETA the distribution matrices \f$ \beta_i \f$
!> @param[in] ZROE nodal values of the dependent variable at the various time levels: \f$ Z(*,1) = Z^{n+1,k}; Z(*,2) = Z^{n}; Z(*,3) = Z^{n-1} \f$ 
!> @param[in] NOFVAR nof degrees of freedom (=1) used for compatibility with the calling sequence of similar routines; note that NOFVAR is the leading dimension of both Z and NODRES; when solving the RANS eqns NOFVAR > NDOF
!> @param[in,out] NODRES nodal residual updated with the contribution of the time-derivative term
!> @param[in,out] STIFC the Jacobian matrix updated with the time derivative term; only if MATRIX_ASSEMBLY == .TRUE.
!> @param[in] VOL the volumes at: current time level (e.g. \c n+a), \c n+1 \c n \c n-1
!> @param[in] NDOF leading dimension of dUdZ used for compatibility with the calling sequence of similar routines
!> @param[in] NDIM dimension of the space
!> @param[in] NOFVERT (\c =NDIM+1) nof of vertices of the current simplex
!> @param[in] MATRIX_ASSEMBLY when .TRUE. entries of STIFC are updated
!>
!> \author $Author: abonfi $
!> \version $Revision: 1.11 $
!> \date $Date: 2020/03/28 09:49:28 $
!> \warning Almost un-tested with \c DUALTS = \c .FALSE.
!>
!>
      SUBROUTINE UNSTEADY4(dUdZ,BETA,ZROE,NOFVAR,NODRES,STIFC,VOL,NDOF,
     &                     NDIM,NOFVERT,MATRIX_ASSEMBLY)
C
C     $Id: unsteady4.F,v 1.11 2020/03/28 09:49:28 abonfi Exp $
C
C
CCC#define DEBUG
CCC#define XPRINT
CCCC#define PRINT
CCCCCCC#define USE_BLAS
C
      IMPLICIT NONE
C
C     unsteady update for compressible flows
C
      INCLUDE 'paramt.h'
      INCLUDE 'constants.h'
      INCLUDE 'time.h'
      INCLUDE 'time.com'
      INCLUDE 'stream.com'
C
C     .. Parameters ..
      INTEGER LWORK
      PARAMETER (LWORK=MAXNOFVERT*MAX_NOFVAR_SQR)
C     ..
C     .. Scalar Arguments ..
      INTEGER NOFVAR,NDOF,NDIM,NOFVERT
      LOGICAL MATRIX_ASSEMBLY
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION dUdZ(NDOF,NDOF,NOFVERT),BETA(NDOF,NDOF,NOFVERT),
     1       ZROE(NOFVAR,NOFVERT,*),NODRES(NOFVAR,*),
     2       STIFC(NDOF,NDOF,NOFVERT,NOFVERT),VOL(*)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION S,HELP,ALPHA,DIVB,CNST
      INTEGER I,J,K,L,IR,IC,IADDR,ORDSQR,INFO,JADDR,IVERT,IVAR,IJADDR,
     &JIADDR,IBGN,IEND
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION MMAT(MAX_NOFVERT_SQR*MAX_NOFVAR_SQR),
     1                 WORKB(LWORK),WKSPC(LWORK),
     3                 dZdt(MAXNOFVAR*MAXNOFVERT),dZdU(MAX_NOFVAR_SQR)
      INTEGER IPIV(MAXNOFVAR)
#ifdef DEBUG
      DOUBLE PRECISION WORKA(LWORK), WORKC(LWORK)
#endif
C     ..
C     .. External Functions ..
      LOGICAL NULLMAT
C     ..
C     .. Statement Functions ..
      INTEGER IDX,JDX
      IDX(I,J) = ((J-1)*NOFVERT+I-1)*ORDSQR+1
      JDX(I,J) =  (J-1)*NDOF+I-1
C
      ORDSQR = NDOF*NDOF
#ifdef PRINT
      CALL R8Mat_Print('G',' ',NDOF,NDOF*NOFVERT,dudZ,
     +                NDOF,'nodal dudZ (<<<) ',INFO)
#endif
      GOTO (10,20,30,20,50) MMTYPE
   10 CONTINUE ! lumped mass matrix
      STOP 'MM-LM yet to be implemented'
      GOTO 100
   20 CONTINUE ! Petrov-Galerkin or Simple-Upwind (up to a certain stage they're the same)
C
C     compute (d+1)<(dU/dZ)> and store in MMAT(1,1)
C
#ifdef USE_BLAS
      CALL DCOPY(ORDSQR,dUdZ(1,1,1),1,MMAT,1)
      DO IVERT = 2,NOFVERT
         CALL DAXPY(ORDSQR,ONE,dUdZ(1,1,IVERT),1,MMAT,1) 
      ENDDO
#else
      DO I = 1,ORDSQR
         MMAT(I) = ZERO
      ENDDO
      DO IVERT = 1,NOFVERT
         DO IC = 1,NDOF
            DO IR = 1,NDOF
               HELP = MMAT(JDX(IR,IC)+1)
               MMAT(JDX(IR,IC)+1) = HELP + dUdZ(IR,IC,IVERT)
            ENDDO
         ENDDO
      ENDDO
#endif
#ifdef PRINT
      CALL R8Mat_Print('G',' ',NDOF,NDOF,MMAT(1),
     +                NDOF,'average dudZ (<<<) ',INFO)
#endif
      IF(MATRIX_ASSEMBLY)THEN
          CALL DCOPY(ORDSQR,MMAT,1,dZdU,1) 
          CALL DSCAL(ORDSQR,ONE/REAL(NOFVERT),dZdU,1)
          CALL DGETRF( NDOF,NDOF, dZdU, NDOF, IPIV, INFO )
          IF(INFO.NE.0)THEN
              CALL EXIT(INFO)
          ENDIF
          CALL DGETRI(NDOF,dZdU,NDOF,IPIV,WORKB,LWORK,INFO)
          IF(INFO.NE.0)THEN
              CALL EXIT(INFO)
          ENDIF
      ENDIF
C
C     0. copy (d+1)<(dU/dZ)> into M_{ij}
C
      DO 3 J = 1,NOFVERT
         DO 3 I = 1,NOFVERT
            IF((I.EQ.1).AND.(J.EQ.1))GOTO 3
            IJADDR = IDX(I,J)
#ifdef USE_BLAS
            CALL DCOPY(ORDSQR,MMAT(1),1,MMAT(IJADDR),1)
#else
            DO IC = 1,NDOF
               DO IR = 1,NDOF
                  MMAT(IJADDR+JDX(IR,IC)) = MMAT(JDX(IR,IC)+1)
               ENDDO
            ENDDO
#endif
#ifdef PRINT
      CALL R8Mat_Print('G',' ',NDOF,NDOF,MMAT(IJADDR),
     +                NDOF,'M_ij (0) ',INFO)
#endif
    3 CONTINUE ! end loop over vertices i,j
C
#ifdef DEBUG
C
C     compute the exact column sum (just for de-bugging)
C
      DO IVERT = 1,NOFVERT
         IADDR = (IVERT-1)*ORDSQR+1
         CALL DCOPY(ORDSQR,dUdZ(1,1,IVERT),1,WORKA(IADDR),1)
         CALL DCOPY(ORDSQR,BETA(1,1,IVERT),1,WORKB(IADDR),1)
      ENDDO
C
C     B_i := B_i-I/[(d+1)(d+2)]
C
      HELP = ONE/(NOFVERT*(NDIM+2))
      DO IVERT = 1,NOFVERT
         IBGN = (IVERT-1)*ORDSQR+1
         IEND = IBGN+ORDSQR-1
         DO I = IBGN,IEND,(NDOF+1) 
            WORKB(I) = WORKB(I) - HELP
         ENDDO
      ENDDO
      DO IVERT = 1,NOFVERT
         IADDR = (IVERT-1)*ORDSQR+1
         CALL DGEMM('No','No',NDOF,NDOF,NDOF,ONE/NOFVERT,WORKB(IADDR),
     &              NDOF,MMAT(1),NDOF,HELP,WORKA(IADDR),NDOF)
      ENDDO
#endif
C
C     compute C_j = (1/[(d+1)(d+2)])[(dU/dZ)_j + (d+1)<(dU/dZ)>]
C
C     1. copy dUdZ_k into C for all k=1,NOFVERT
C     2. copy B_i    into B for all k=1,NOFVERT
C
      CALL DCOPY(ORDSQR*NOFVERT,dUdZ(1,1,1),1,WKSPC,1)
      CALL DCOPY(ORDSQR*NOFVERT,BETA(1,1,1),1,WORKB,1)
C
#ifdef PRINT
      CALL R8Mat_Print('G',' ',NDOF,NDOF*NOFVERT,WKSPC,
     +                NDOF,'C_j (1) ',INFO)
      CALL R8Mat_Print('G',' ',NDOF,NDOF*NOFVERT,WORKB,
     +                NDOF,'B_i (2) ',INFO)
#endif
C
C     3. C_j := [(dU/dZ)_j + (d+1)<(dU/dZ)>]
C
      DO IVERT =1, NOFVERT
          IADDR = (IVERT-1)*ORDSQR+1
          CALL DAXPY(ORDSQR,ONE,MMAT(1),1,WKSPC(IADDR),1)
      ENDDO 
C
#ifdef PRINT
      CALL R8Mat_Print('G',' ',NDOF,NDOF*NOFVERT,WKSPC,
     +                NDOF,'C_j (3) ',INFO)
#endif
C
C     4. C_j := 1/[(d+1)(d+2)] * C_j
C
      HELP = ONE/(NOFVERT*(NDIM+2))
#ifdef USE_BLAS
      CALL DSCAL(ORDSQR*NOFVERT,HELP,WKSPC,1)
#else
      DO I = 1,ORDSQR*NOFVERT
           WKSPC(I) = HELP * WKSPC(I)
      ENDDO
#endif
#ifdef PRINT
      CALL R8Mat_Print('G',' ',NDOF,NDOF*NOFVERT,WKSPC,
     +                NDOF,'C_j (4) ',INFO)
#endif
C
C     Up to this point MM-PG and MM-SU are the same;
C     for MM-SU we can however skip certain operations
C
      IF(MMTYPE.EQ.4)GOTO 40
C
C     Here we do MM-PG
C
C     5. B_i := B_i - I/(d+1)
C
      HELP = ONE/NOFVERT
      DO IVERT = 1,NOFVERT
         IBGN = (IVERT-1)*ORDSQR+1
         IEND = IBGN+ORDSQR-1
         DO I = IBGN,IEND,(NDOF+1) 
            WORKB(I) = WORKB(I) - HELP
         ENDDO
      ENDDO
#ifdef PRINT
      CALL R8Mat_Print('G',' ',NDOF,NDOF*NOFVERT,WORKB,
     +                NDOF,'B_j (5) ',INFO)
#endif
C
C     the cost should be: (d+1)**2 GEMM + [(d+1)**3 +(d+1)] AXPY
C
C
      DO J = 1,NOFVERT
         JADDR = ORDSQR*(J-1)+1
         DO I = 1,J
            IADDR = ORDSQR*(I-1)+1
            IF(I.EQ.J)THEN
               ALPHA = TWO/(NDIM+3)
               HELP = -TWO/(10.d0*NDIM)
            ELSE
               ALPHA = ONE/(NDIM+3)
               HELP = MONE/(10.d0*NDIM)
            ENDIF
            HELP = HELP/NOFVERT
            IJADDR = IDX(I,J)
#ifdef USE_BLAS
            CALL DSCAL(ORDSQR,HELP,MMAT(IJADDR),1)
            CALL DAXPY(ORDSQR,ALPHA,WKSPC(IADDR),1,MMAT(IJADDR),1)
            CALL DAXPY(ORDSQR,ALPHA,WKSPC(JADDR),1,MMAT(IJADDR),1)
#else
            DO IC = 1,NDOF
               DO IR = 1,NDOF
                  S = HELP*MMAT(IJADDR+JDX(IR,IC))
                  S = S+(WKSPC(IADDR+JDX(IR,IC))+
     &                   WKSPC(JADDR+JDX(IR,IC)))*ALPHA
                  MMAT(IJADDR+JDX(IR,IC)) = S
               ENDDO 
            ENDDO
#endif
#ifdef PRINT
         WRITE(6,*)'I = ',I,' J = ',J 
         CALL R8Mat_Print('G',' ',NDOF,NDOF,MMAT(IJADDR),
     +                NDOF,'Mass matrix ',INFO)
#endif
         ENDDO ! loop over I
      ENDDO ! loop over J
C
C     Here we take advantage of symmetry
C
      DO J = 1,NOFVERT
         JADDR = ORDSQR*(J-1)+1
         DO I = J+1,NOFVERT
            IJADDR = IDX(I,J)
            JIADDR = IDX(J,I)
#ifdef USE_BLAS
            CALL DCOPY(ORDSQR,MMAT(JIADDR),1,MMAT(IJADDR),1)
#else
            DO IC = 1,NDOF
               DO IR = 1,NDOF
                  MMAT(IJADDR+JDX(IR,IC)) = MMAT(JIADDR+JDX(IR,IC))
               ENDDO 
            ENDDO
#endif
#ifdef PRINT
            WRITE(6,*)'I = ',I,' J = ',J 
            CALL R8Mat_Print('G',' ',NDOF,NDOF,MMAT(IJADDR),
     +                NDOF,'Mass matrix ',INFO)
            CALL R8Mat_Print('G',' ',NDOF,NDOF,MMAT(JIADDR),
     +                NDOF,'Mass matrix ',INFO)
#endif
          ENDDO ! end loop over i
       ENDDO ! end loop over j
C
   40 CONTINUE ! this is for both MM-SU and MM-PG
C
C     performs (d+1)**2 GEMMs
C     M_{ij} = M_{ij} + B_i x C_j for MM-PG
C     M_{ij} =          B_i x C_j for MM-SU
C
      IF(MMTYPE.EQ.2)THEN
         CNST = ONE 
      ELSE
         CNST = ZERO 
      ENDIF
      DO J = 1,NOFVERT 
         JADDR = ORDSQR*(J-1)+1
         DO I = 1,NOFVERT 
            IADDR = ORDSQR*(I-1)+1
            IJADDR = IDX(I,J)
#ifdef USE_BLAS
            CALL DGEMM('N','N',NDOF,NDOF,NDOF,ONE,WORKB(IADDR),
     &      NDOF,WKSPC(JADDR),NDOF,CNST,MMAT(IJADDR),NDOF)
#else
            DO IC = 1,NDOF
               DO IR = 1,NDOF
                  IF(MMTYPE.EQ.2)THEN
                     HELP = MMAT(IJADDR+JDX(IR,IC)) ! MM-PG
                  ELSE
                     HELP = ZERO ! MM-SU
                  ENDIF
                  DO L = 1,NDOF
                     HELP = HELP+WORKB(IADDR+JDX(IR,L))*
     &                           WKSPC(JADDR+JDX(L,IC))
                  ENDDO 
                  MMAT(IJADDR+JDX(IR,IC)) = HELP
               ENDDO 
            ENDDO 
#endif
#ifdef PRINT
         WRITE(6,*)'I = ',I,' J = ',J 
         CALL R8Mat_Print('G',' ',NDOF,NDOF,MMAT(IADDR),
     +                NDOF,'Mass matrix ',INFO)
#endif
         ENDDO
      ENDDO
      GOTO 45
   45 CONTINUE ! this is for MM-SU 
C
#ifdef DEBUG
C
C     compute the column sum
C
      CALL DINIT(ORDSQR*NOFVERT,ZERO,WORKB,1)
      DO I = 1,NOFVERT
         IADDR = (I-1)*ORDSQR+1
         DO J = 1,NOFVERT
            IJADDR = IDX(I,J)
            CALL DAXPY(ORDSQR,ONE,MMAT(IJADDR),1,WORKB(IADDR),1)
         ENDDO
      ENDDO
C
      CALL DCOPY(ORDSQR*NOFVERT,WORKA,1,WORKC,1)
      CALL DAXPY(ORDSQR*NOFVERT,MONE,WORKB,1,WORKC,1)
      IF( .NOT. NULLMAT(WORKC,NDOF,NDOF,NDOF,1.D-13) )THEN
C
      CALL R8Mat_Print('G',' ',NDOF,NDOF*NOFVERT,WORKB,
     +                NDOF,'ROWSUM matrix ',INFO)
      CALL R8Mat_Print('G',' ',NDOF,NDOF*NOFVERT,WORKA,
     +                NDOF,'exact ROWSUM matrix ',INFO)
      CALL R8Mat_Print('G',' ',NDOF,NDOF*NOFVERT,WORKC,
     +                NDOF,'error matrix ',INFO)
          pause
      ENDIF
C
#endif
      GOTO 100
   30 CONTINUE ! Consistent Upwind
      STOP 'MM-CU yet to be implemented'
      GOTO 100
   50 CONTINUE ! Centred
      STOP 'MM-CC yet to be implemented'
  100 CONTINUE
c
#ifdef PRINT
      DO J = 1,NOFVERT
      DO I = 1,NOFVERT
         IJADDR = IDX(I,J)
         WRITE(6,*)'I = ',I,' J = ',J 
         CALL R8Mat_Print('G',' ',NDOF,NDOF,MMAT(IJADDR),
     +                NDOF,'Mass matrix ',INFO)
!         CALL R8Mat_Print('G',' ',NOFVERT,NOFVERT,STIFC(1,1),
!    +                NOFVERT,'C_ij matrix ',INFO)
      ENDDO
      ENDDO
#endif
C
C
C
      IF(DUALTS)THEN ! dual time stepping
C
C update the rhs by adding the contribution
C of the time derivative term
C
caldo
!     CALL R8Mat_Print('G',' ',ndof,NOFVERT*NTIMLEVS,ZROE,
!    +                ndof,'Z matrix ',INFO)
!     CALL R8Mat_Print('G',' ',ndof,NOFVERT*NTIMLEVS,ZROE,
!    +                ndof,'ZROE matrix before  ',INFO)
!     CALL R8Mat_Print('G',' ',ndof,NOFVERT*NTIMLEVS,NODRES,
!    +                NOFVAR,'R matrix before  ',INFO)
caldo
C
C     dZdt(1:NDOF,1:NOFVERT) is the time derivative of the parameter vector
C
         IF(LALE)THEN 
            CALL DINIT(NDOF*NOFVERT,ZERO,dZdt,1)
            DIVB = ZERO
            DO K = 1,NTIMLEVS
               DIVB = DIVB + TCOEF(2-K)*VOL(K+1)
            ENDDO
            DIVB = DIVB/DELT
            DO IVERT = 1,NOFVERT
               DO K = 1,NTIMLEVS
!                 ALPHA = -TCOEF(2-K)/DTVOL ! *VOL(K+1)/DELT
                  ALPHA = -TCOEF(2-K)*VOL(K+1)/DELT! /DTVOL ! *VOL(K+1)/DELT
                  DO IVAR = 1,NDOF 
                     IADDR = (IVERT-1)*NDOF+IVAR
                     dZdt(IADDR) = dZdt(IADDR) + ALPHA * 
     &                             ZROE(IVAR,IVERT,K)
!    &                             U_infty(IVAR)
                  ENDDO
               ENDDO
            ENDDO
!     WRITE(6,*)'Vol(b) = ',(VOL(K),K=1,NTIMLEVS+1) 
!     WRITE(6,*)'Div(b) = ',divb 
!     CALL R8Mat_Print('G',' ',ndof,NTIMLEVS,dZdt,
!    +                NDOF,'matrix dZdt  ',INFO)
         ELSE ! fixed grid
C
            HELP = ONE/DTVOL
            CALL DINIT(NDOF*NOFVERT,ZERO,dZdt,1)
            DO IVERT = 1,NOFVERT
               DO K = 1,NTIMLEVS
                  ALPHA =-HELP*TCOEF(2-K)
                  DO IVAR = 1,NDOF 
                     IADDR = (IVERT-1)*NDOF+IVAR
                     dZdt(IADDR) = dZdt(IADDR) + ALPHA * 
     &                             ZROE(IVAR,IVERT,K)
                  ENDDO
               ENDDO
            ENDDO
#ifdef PRINT
         CALL R8Mat_Print('G',' ',ndof,NOFVERT,dzdt,ndof,'dzdt matrix',INFO)
#endif
         ENDIF
C             
C     unsteady contribution: applies to both fixed and moving meshes
C             
#ifdef XPRINT
         CALL DINIT(NDOF*NOFVERT,ZERO,WORKA,1)
         CALL DINIT(NDOF*NOFVERT,ZERO,WORKB,1)
#endif
         DO I = 1,NOFVERT ! loop over vertices
            DO J = 1,NOFVERT ! loop over vertices
               IADDR = IDX(I,J)
               JADDR = (J-1)*NDOF+1
               CALL DGEMV('No Transpose',NDOF,NDOF,ONE,MMAT(IADDR), ! notice that we account in ALPHA for the (-) sign
     2                    NDOF,dZdt(JADDR),1,ONE,NODRES(1,I),1)
            ENDDO ! j
         ENDDO ! i
#ifdef XPRINT
      CALL R8Mat_Print('G',' ',ndof,NOFVERT,dzdt,ndof,'dzdt matrix',INFO)
      CALL R8Mat_Print('G',' ',ndof,Nofvert,NODRES,
     +                nofvar,'R vector before ',INFO)
      write(6,*)'Div(b) = ',divb
#endif
C
C     this is the term: div (b) \sum m_{ij} Z_j^{n+1,k}
C
         IF(LALE)THEN
            DO I = 1,NOFVERT ! loop over vertices
               DO J = 1,NOFVERT ! loop over vertices
                  IADDR = IDX(I,J)
                  CALL DGEMV('No Transpose',NDOF,NDOF,DIVB,
     2                       MMAT(IADDR),NDOF,ZROE(1,J,1),1,ONE,
     3                       NODRES(1,I),1)
               ENDDO ! j
            ENDDO ! i
#ifdef XPRINT
      CALL R8Mat_Print('G',' ',ndof,Nofvert,NODRES,
     +                nofvar,'R vector after ',INFO)
      CALL DAXPY(NDOF*NOFVERT,ONE,WORKB,1,WORKA,1) 
      IF( .NOT. NULLMAT(WORKA,NDOF,NOFVERT,NDOF,1.D-13) )THEN
      WRITE(6,*)'Uh! Oh! *****'
      CALL R8Mat_Print('G',' ',ndof,Nofvert,WORKA,
     +                nofvar,'WORKA vector ',INFO)
!     CALL R8Mat_Print('G',' ',ndof,Nofvert,WORKB,
!    +                nofvar,'WORKB vector ',INFO)
      ENDIF
#endif
         ENDIF
C
C update the matrix by adding the mass matrix
C
         IF(MATRIX_ASSEMBLY)THEN
!           HELP =-(ONE+HALF*GAMT)/(NOFVERT*DTVOL)
!           HELP = +TCOEF(1)/DTVOL
            HELP = -TCOEF(1)/DTVOL
!           HELP = HELP/NOFVERT
            HELP = 0.d0
            DO J = 1,NOFVERT
               DO I = 1,NOFVERT
                  IADDR = IDX(I,J)
!     CALL R8Mat_Print('G',' ',ndof,Ndof,STIFC(1,1,I,J),
!    +                ndof,'C matrix before ',INFO)
cold              CALL DAXPY(ORDSQR,HELP,MMAT(IADDR),1,STIFC(1,1,I,J),1)
            CALL DGEMM('N','N',NDOF,NDOF,NDOF,HELP,MMAT(IADDR),
     &      NDOF,dZdU,NDOF,ONE,STIFC(1,1,I,J),NDOF)
!                                             ^
!                                             | use NOFVAR for turbulent flow cases?
!     CALL R8Mat_Print('G',' ',ndof,Ndof,STIFC(1,1,I,J),
!    +                ndof,'C matrix after ',INFO)
               ENDDO ! loop over I
            ENDDO ! loop over J
!           pause
         ENDIF
      ELSE ! no inner iterations (DUALTS == .FALSE.)
         STOP 'Unimplemented feature: Must use dual ts !'
      ENDIF 
      RETURN
      END
