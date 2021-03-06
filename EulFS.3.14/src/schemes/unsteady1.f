!> \par Purpose
!>
!> Add the time-dependent term for unsteady scalar problems
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
!> +\left(\gamma_t/2\right)|T^e|^{n-1}u_j^{n-1}}{\Delta t} -
!> \frac{\left(1+\gamma_t/2\right)|T^e|^{n+1}-\left(1+\gamma_t\right)|T^e|^{n} +\left(\gamma_t/2\right)|T^e|^{n-1}}{\Delta t} 
!> \sum_{j} m_{ij} u_j^{n+1,k}
!> \f]
!>
!>
!> various kinds of mass-matrices are assembled, depending on the value of \c MYTYPE and according to the table below:
!>
!> \f[
!> \begin{array}{c|c|c}
!> {\tt MYTYPE} & \omega^e_i\left(\mathbf{x}\right) & m^e_{ij} \\\\\hline
!> {\tt MM\_LUMPED} & & \frac{\delta_{ij}}{d+1} |T_e| \\\
!> {\tt MM\_PETROV\_GALERKIN} & N^e_i\left(\mathbf{x}\right) + \beta_i^e - \frac{1}{d+1} &  \frac{|T_e|}{d+1} \left( \frac{1+\delta_{ij}}{d+2} + \beta_i^e - \frac{1}{d+1} \right) \\\
!> {\tt MM\_CONS\_UPWIND} & & \frac{1}{d+1} \beta_i \left(1+\delta_{ij}-\beta_j\right) |T_e| \\\
!> {\tt MM\_SIMPLE\_UPWIND} & \beta_i & \frac{1}{d+1} \beta_i |T_e| \\\
!> {\tt MM\_CENTRED} & & \frac{|T_e|}{d+1} \\\
!> {\tt MM\_NEW} & &
!> \end{array}
!> \f]
!>
!
!> @param[in] dUdZ dummy argument used for compatibility with the calling sequence of similar routines
!> @param[in] BETA the distribution matrices \f$ \beta_i \f$
!> @param[in] Z nodal values of the dependent variable at the various time levels: \f$ Z(*,1) = u^{n+1,k}; Z(*,2) = u^{n}; Z(*,3) = u^{n-1} \f$ 
!> @param[in] NOFVAR nof degrees of freedom (=1) used for compatibility with the calling sequence of similar routines
!> @param[in,out] NODRES nodal residual updated with the contribution of the time-derivative term
!> @param[in,out] STIFC the Jacobian matrix updated with the time derivative term; only if MATRIX_ASSEMBLY == .TRUE.
!> @param[in] VOL the volumes at: current time level (e.g. \c n+a), \c n+1 \c n \c n-1
!> @param[in] NDOF leading dimension of dUdZ used for compatibility with the calling sequence of similar routines
!> @param[in] NDIM dimension of the space
!> @param[in] NOFVERT (\c =NDIM+1) nof of vertices of the current simplex
!> @param[in] MYTYPE type of mass-matrix to be computed, see time.h
!> @param[in] MATRIX_ASSEMBLY when .TRUE. entries of STIFC are updated
!>
!> \author $Author: abonfi $
!> \version $Revision: 1.9 $
!> \date $Date: 2020/03/28 09:49:28 $
!> \warning Almost un-tested with \c DUALTS = \c .FALSE.
!>
!>
      SUBROUTINE UNSTEADY1(dUdZ,BETA,Z,NOFVAR,NODRES,STIFC,VOL,
     2                     NDOF,NDIM,NOFVERT,MYTYPE,MATRIX_ASSEMBLY)
C
      IMPLICIT NONE
C
C     $Id: unsteady1.f,v 1.9 2020/03/28 09:49:28 abonfi Exp $
C
C
      INCLUDE 'paramt.h'
      INCLUDE 'constants.h'
      INCLUDE 'time.h'
      INCLUDE 'time.com'
C
C     .. Parameters ..
C     ..
C     .. Scalar Arguments ..
      INTEGER NDIM,NOFVERT,NOFVAR,NDOF,MYTYPE
      LOGICAL MATRIX_ASSEMBLY
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION dUdZ(NDOF,*),BETA(NOFVERT),Z(NOFVERT,*),
     +                 NODRES(NOFVERT),STIFC(NOFVERT,NOFVERT),VOL(*)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION RD,S,HELP,DIAG,OFFD,DIVB
      INTEGER I,J,IADDR,INFO,ILEV
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION MMAT(MAX_NOFVERT_SQR),DUDT(MAXNOFVERT)
C     ..
C     .. External Functions ..
      DOUBLE PRECISION DDOT
      EXTERNAL DDOT
C     ..
C
      IF(LALE)THEN
      ENDIF
      HELP = (ONE+HALF*GAMT)/(NOFVERT*DTVOL)
      RD = ONE/REAL(NOFVERT)
      DIAG = (12.d0-NDIM)/60.d0
      OFFD = (TWO*NDIM-9.d0)/60.d0
      GOTO (10,20,30,40,50,60) MYTYPE
C
C     MMAT is missing the term VOL/(d+1) which is add only later
C
   10 CONTINUE ! lumped mass matrix
      IADDR = 0
      DO 5 J = 1,NOFVERT
         DO 5 I = 1,NOFVERT
            IADDR = IADDR + 1 
            IF(J.EQ.I)THEN
               MMAT(IADDR) = ONE
            ELSE
               MMAT(IADDR) = ZERO
            ENDIF 
    5 CONTINUE
      GOTO 100
   20 CONTINUE ! Petrov-Galerkin
      IADDR = 0
      DO 3 J = 1,NOFVERT
         DO 3 I = 1,NOFVERT
            IADDR = IADDR + 1 
            IF(J.EQ.I)THEN
               MMAT(IADDR) = BETA(I)+DIAG
            ELSE
               MMAT(IADDR) = BETA(I)+OFFD
            ENDIF 
    3 CONTINUE
      GOTO 100
   30 CONTINUE ! Consistent Upwind
      IADDR = 0
      DO 1 J = 1,NOFVERT
         DO 1 I = 1,NOFVERT
            IADDR = IADDR + 1 
            IF(J.EQ.I)THEN
               MMAT(IADDR) = BETA(I)*(TWO-BETA(J))
            ELSE
               MMAT(IADDR) = BETA(I)*(ONE-BETA(J))
            ENDIF 
    1 CONTINUE
      GOTO 100
   40 CONTINUE ! Simple Upwind
      IADDR = 0
      DO 7 J = 1,NOFVERT
         DO 7 I = 1,NOFVERT
            IADDR = IADDR + 1 
            MMAT(IADDR) = BETA(I)
    7 CONTINUE
      GOTO 100
   50 CONTINUE ! Centred
      IADDR = 0
      DO 9 J = 1,NOFVERT
         DO 9 I = 1,NOFVERT
            IADDR = IADDR + 1 
            MMAT(IADDR) = RD
    9 CONTINUE
      GOTO 100
   60 CONTINUE ! Aldo Upwind
      IADDR = 0
      DO 19 J = 1,NOFVERT
         DO 19 I = 1,NOFVERT
            IADDR = IADDR + 1 
            IF(J.EQ.I)THEN
               MMAT(IADDR) = BETA(I)*BETA(J)
            ELSE
               MMAT(IADDR) = BETA(I)*(ONE+BETA(J))
            ENDIF 
   19 CONTINUE
      GOTO 100
!
!         CALL R8Mat_Print('G',' ',NOFVERT,NOFVERT,MMAT(1),
!    +                NOFVERT,'Mass matrix ',INFO)
!         CALL R8Mat_Print('G',' ',NOFVERT,NOFVERT,STIFC(1,1),
!    +                NOFVERT,'C_ij matrix ',INFO)
  100 CONTINUE
      IF(DUALTS)THEN ! dual time stepping
C
C update the rhs by adding the contribution of the time derivative term
C
!     CALL R8Mat_Print('G',' ',NOFVERT,NTIMLEVS,Z(1,1),
!    +                NOFVERT,'Z matrix ',INFO)
C
         DO I = 1,NOFVERT
            DUDT(I) = ZERO
         ENDDO
C
         IF(LALE)THEN
C
            HELP = TCOEF(1)*VOL(2) ! (1.+\gamma_t/2)T_e^{n+1}
            DIVB = HELP
            CALL DAXPY(NOFVERT,HELP,Z(1,1),1,DUDT,1) ! add time level n+1,k
C
            HELP = TCOEF(0)*VOL(3) ! -(1.+\gamma_t)T_e^{n}
            DIVB = DIVB + HELP
            CALL DAXPY(NOFVERT,HELP,Z(1,2),1,DUDT,1) ! add time level n
            IF(NTIMLEVS.EQ.3)THEN
               HELP = TCOEF(-1)*VOL(4) ! (\gamma_t/2)T_e^{n-1}
               CALL DAXPY(NOFVERT,HELP,Z(1,3),1,DUDT,1) ! add time level n-1
               DIVB = DIVB + HELP
            ENDIF
            DIVB = DIVB/DELT
!        WRITE(6,*)'vol = ',(VOL(i),i=1,4),' div(b) = ',DIVB
            CALL DGEMV('No',NOFVERT,NOFVERT,RD/DELT,MMAT,
     &              NOFVERT,DUDT,1,MONE,NODRES,1) ! update nodres with the term: \sum_j m_{ij} \frac{d}{dt} |T_e| u_j
!           divb = zero
            CALL DGEMV('No',NOFVERT,NOFVERT,DIVB*RD,MMAT,
     &              NOFVERT,Z(1,1),1,MONE,NODRES,1) ! update nodres with the term: div(b) \sum_j m_{ij} u_j^{n+1,k}
         ELSE
            HELP = ONE/(NOFVERT*DTVOL)
            DO ILEV = 1, NTIMLEVS ! loop over the time-levels
               CALL DAXPY(NOFVERT,TCOEF(2-ILEV),Z(1,ILEV),1,DUDT,1)
            ENDDO
               CALL DGEMV('No',NOFVERT,NOFVERT,-HELP,MMAT,
     &              NOFVERT,DUDT,1,ONE,NODRES,1) ! update nodres
         ENDIF
!     CALL R8Mat_Print('G',' ',NOFVERT,NTIMLEVS,NODRES,
!    +                NOFVERT,'update matrix ',INFO)
C
C update the matrix by adding the mass matrix
C
         IF(MATRIX_ASSEMBLY)THEN
            HELP = (ONE+HALF*GAMT)/(NOFVERT*DTVOL)
            IADDR = 0
            DO J = 1,NOFVERT
               DO I = 1,NOFVERT
                  IADDR = IADDR + 1 
                  STIFC(I,J) = THETAT*STIFC(I,J)-HELP*MMAT(IADDR)
               ENDDO
            ENDDO
         ENDIF
C
      ELSE ! no inner iterations (DUALTS == .FALSE.)
C
C update the matrix by adding the mass matrix
C
         IADDR = 0
         DO J = 1,NOFVERT
            DO I = 1,NOFVERT
               IADDR = IADDR + 1 
               STIFC(I,J) = THETAT*STIFC(I,J)-HELP*MMAT(IADDR)
            ENDDO
         ENDDO
C
C update the rhs by adding the contribution
C from the previous time levels
C
         IF(GAMT.NE.ZERO)THEN
            HELP = (HALF*GAMT)/(NOFVERT*DTVOL)
C
C        compute u^n-u^{n-1}
C
            DO J = 1, NOFVERT
               Z(J,2) = Z(J,2) - Z(J,3) 
            ENDDO
            DO I = 1,NOFVERT
               S = ZERO 
               DO J = 1,NOFVERT
                  IADDR = (J-1)*NOFVERT+I
                  S = S + MMAT(IADDR)*Z(J,2)
               ENDDO
               NODRES(I) = NODRES(I) + HELP*S
            ENDDO
         ENDIF 
      ENDIF ! DUALTS
      RETURN
      END
