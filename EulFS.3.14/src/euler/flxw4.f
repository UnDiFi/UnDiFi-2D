!> \brief \b FLXW4
!
!> \par Purpose
!>
!> compute inviscid wall b.c.'s for compressible flows: weak approach
!>
!> @param[in] NDIM the dimension of the space
!> @param[in] NOFVAR nof dofs
!> @param[in] NOFVERT nof vertices of a cell (=NDIM+1)
!> @param[in] STIFC the jacobian matrix
!> @param[in] WORK work array
!> @param[in] WORK2 work array
!> @param[in] VCZ nodal values in the NOFVERT vertices
!> @param[in] VCB grid velocity at all NOFVERT vertices
!> @param[in] VCN NDIM components of the normal to the boundary
!> @param[out] NODRES nodal residual
!> @param[in] PICARD .TRUE. if the jacobian matrix has to be assembled
!     
!>
!> In two-dimensional flows we use:
!>
!> \f{eqnarray*}{
!> R_1 &:=& R_1 + \frac{\alpha}{2} F_1 + \frac{1-\alpha}{2} F_2 \\
!> R_2 &:=& R_2 + \frac{\alpha}{2} F_2 + \frac{1-\alpha}{2} F_1
!>   \f} 
!     
!>
!> In three-dimensional flows we use:
!>
!> \f{eqnarray*}{
!> R_1 &:=& R_1 + \frac{\alpha}{3} F_1 + \frac{1-\alpha}{6} F_2 + \frac{1-\alpha}{6} F_3 \\
!> R_2 &:=& R_2 + \frac{\alpha}{3} F_2 + \frac{1-\alpha}{6} F_3 + \frac{1-\alpha}{6} F_1 \\
!> R_3 &:=& R_3 + \frac{\alpha}{3} F_3 + \frac{1-\alpha}{6} F_1 + \frac{1-\alpha}{6} F_2
!>   \f} 
!>
!> \f$\alpha\f$ is hardwired to 3/4
!     
!> \author $Author: abonfi $
!> \version $Revision: 1.10 $
!> \date $Date: 2020/03/28 09:51:15 $
!
      SUBROUTINE FLXW4(NDIM,NOFVAR,NOFVERT,STIFC,WORK,WORK2,VCZ,VCB,VCN,
     +                 NODRES,PICARD)
C
C     $Id: flxw4.f,v 1.10 2020/03/28 09:51:15 abonfi Exp $
C
      IMPLICIT NONE
C
C
C
      include 'paramt.h'
      include 'constants.h'
C
C     .. Parameters ..
      DOUBLE PRECISION ALPHA
      PARAMETER (ALPHA=0.75d0)
C     ..
C     .. Scalar Arguments ..
      INTEGER NDIM,NOFVAR,NOFVERT
      LOGICAL PICARD
C     NDIM   dimension of the space (2 or 3)
C     NOFVAR number of variables (degrees of freedom)
C            in each meshpoint
C     NOFVAR number of variables (degrees of freedom)
C            in each meshpoint
C     PICARD .TRUE. for Picard linearisation
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION NODRES(NOFVAR,NOFVERT),
     +                 STIFC(NOFVAR,NOFVAR,NOFVERT,NOFVERT),VCN(NDIM),
     +                 VCZ(NOFVAR,NOFVERT),WORK(NDIM+2,NDIM+2,NDIM),
     +                 WORK2(NDIM+2,NDIM+2,NDIM,NDIM),VCB(NDIM,NOFVERT)
C
C     On entry:
C     --------
C     VCN(1:NDIM) cartesian components of the normal
C                 to the boundary face
C     VCZ(1:NOFVAR,1:NOFVERT) dependent variables in the vertices
C                           of the current (IELEM) element
C                           the freestream values need
C                           to be stored in VCZ(1:NOFVAR,NOFVERT)
C     Upon return:
C     -----------
C     STIFC(1:NOFVAR,1:NOFVAR,1:NOFVERT,1:NOFVERT) 
C                   convection matrix
C                   in the NOFVERT-1 vertices of the boundary face
C     NODRES(1:NOFVAR,1:NOFVERT-1) nodal residual due to the incoming
C                   characteristics in the NOFVERT-1 vertices 
C                   of the boundary face
C
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION BETA,CNST,TMP
      INTEGER I,IADD,IFAIL,IVERT,J,JVERT,K,L,NOFEQN
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION DZDU(75),FLUX(MAXNOFVAR*MAXNOFVERT)
C     ..
C     .. External Subroutines ..
      EXTERNAL DAXPY,DGEMM,DINIT,GETDF4CORRDU,INVWLL,MATDZDU
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC REAL
C     ..
C     .. Data statements ..

      DATA DZDU/75*0.d0/
      INTEGER IPATCH
      COMMON/MYBUG/IPATCH
C     ..
C
C     Compute correction flux
C
      NOFEQN = NDIM+2
C
      DO 3 IVERT = 1,NOFVERT
          IADD = (IVERT-1)*NOFVAR + 1
          CALL INVWLL(NDIM,VCN,VCB(1,IVERT),VCZ(1,IVERT),FLUX(IADD))
    3 CONTINUE
      IF(IPATCH.EQ.-2)THEN
      WRITE(6,*)'Current patch is ',IPATCH
      CALL R8Mat_Print('G',' ',Nofvar,Nofvert,VCZ(1,1),Nofvar,
     +      'Nodal values of Z ',IFAIL)
      CALL R8Mat_Print('G',' ',NDIM,NDIM,VCB(1,1),NDIM,
     +      'Nodal values of B ',IFAIL)
      CALL R8Mat_Print('G',' ',NOFVAR,Nofvert,FLUX,NOFVAR,
     +      'Nodal values of the flux ',IFAIL)
      WRITE(6,*)
      ENDIF
C
      CALL DINIT(NOFVAR* (NOFVERT-1),ZERO,NODRES,1)
      BETA = (ONE-ALPHA)/REAL(NDIM-1)/REAL(NDIM)
      DO 5 IVERT = 1,NOFVERT - 1
          DO 7 JVERT = 1,NOFVERT - 1
              IADD = (JVERT-1)*NOFVAR + 1
              IF (JVERT.EQ.IVERT) THEN
                  CNST = ALPHA/NDIM

              ELSE
                  CNST = BETA
              ENDIF

              CALL DAXPY(NOFEQN,CNST,FLUX(IADD),1,NODRES(1,IVERT),1)
    7     CONTINUE
    5 CONTINUE
C
      IF (.NOT.PICARD) RETURN
C
C     Compute matrices dFcorr/dZ for all NOFVERT-1 vertices
C     of the boundary face
C
      DO 1 IVERT = 1,NDIM
          CALL GETDF4CORRDU(VCZ(1,IVERT),VCB(1,IVERT),VCN,NDIM,NOFEQN,
     +                      WORK(1,1,IVERT))
    1 CONTINUE

      BETA = (ONE-ALPHA)/REAL(NDIM-1)
      DO 8 I = 1,NDIM
          DO 8 J = 1,NDIM
              IF (J.EQ.I) THEN
                  TMP = ALPHA/REAL(NDIM)

              ELSE
                  TMP = BETA/REAL(NDIM)
              ENDIF

              DO 8 L = 1,NOFEQN
                  DO 8 K = 1,NOFEQN
                      WORK2(K,L,I,J) = TMP*WORK(K,L,J)
    8 CONTINUE
C
C
C     transform the convection stiffness matrix into
C     conserved variables as C_{ij} := 2 C_{ij} dZdU(j)
C
      DO 9 IVERT = 1,NDIM
          IADD = (IVERT-1)*NOFEQN*NOFEQN + 1
          CALL MATDZDU(VCZ(1,IVERT),DZDU(IADD),NDIM,NOFEQN)
          DO 9 I = 1,NDIM
              CALL DGEMM('No transpose','No transpose',NOFEQN,NOFEQN,
     +                   NOFEQN,ONE,WORK2(1,1,I,IVERT),NOFEQN,
     +                   DZDU(IADD),NOFEQN,ZERO,STIFC(1,1,I,IVERT),
     +                   NOFVAR)
C
!     write(6,*)ivert,i
!     CALL R8Mat_Print('G',' ',Nofvar,Nofvar,STIFC(1,1,I,IVERT),Nofvar,
!    +      'C(i,j) ',IFAIL)
!     pause
C
    9 CONTINUE
C
      RETURN

  564 FORMAT ((E12.6,1X))

      END
!> \brief \b FLXW44Ar
!
!> \par Purpose
!>
!> compute inviscid wall b.c.'s for an Argon plasma: weak approach
!>
!> @param[in] NDIM the dimension of the space
!> @param[in] NOFVAR nof dofs
!> @param[in] NOFVERT nof vertices of a cell (=NDIM+1)
!> @param[in] STIFC the jacobian matrix
!> @param[in] WORK work array
!> @param[in] WORK2 work array
!> @param[in] VCZ nodal values in the NOFVERT vertices
!> @param[in] VCN NDIM components of the normal to the boundary
!> @param[out] NODRES nodal residual
!> @param[in] PICARD .TRUE. if the jacobian matrix has to be assembled
!     
!> \author $Author: abonfi $
!> \version $Revision: 1.10 $
!> \date $Date: 2020/03/28 09:51:15 $
!
      SUBROUTINE FLXW44Ar(NDIM,NOFVAR,NOFVERT,STIFC,WORK,WORK2,VCZ,VCN,
     +                 NODRES,PICARD)
C
      IMPLICIT NONE
C
C
C     Purpose: 
C     ------
C     compute inviscid wall b.c.'s for plasma flows
C
      include 'paramt.h'
      include 'constants.h'
      include 'plasma.h'
C
C     .. Parameters ..
      INTEGER LENB
      DOUBLE PRECISION ALPHA
      PARAMETER (LENB=MAXNOFVERT*MAX_NOFVAR_SQR,ALPHA=0.75d0)
C     ..
C     .. Scalar Arguments ..
      INTEGER NDIM,NOFVAR,NOFVERT
      LOGICAL PICARD
C     NDIM   dimension of the space (2 or 3)
C     NOFVAR number of variables (degrees of freedom)
C            in each meshpoint
C     NOFVAR number of variables (degrees of freedom)
C            in each meshpoint
C     PICARD .TRUE. for Picard linearisation
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION NODRES(NOFVAR,NOFVERT),
     +                 STIFC(NOFVAR,NOFVAR,NOFVERT,NOFVERT),VCN(NDIM),
     +                 VCZ(NOFVAR,NOFVERT),WORK(NOFVAR,NOFVAR,NDIM),
     +                 WORK2(NOFVAR,NOFVAR,NDIM,NDIM)
C
C     On entry:
C     --------
C     VCN(1:NDIM) cartesian components of the normal
C                 to the boundary face
C     VCZ(1:NOFVAR,1:NOFVERT) dependent variables in the vertices
C                           of the current (IELEM) element
C                           the freestream values need
C                           to be stored in VCZ(1:NOFVAR,NOFVERT)
C     Upon return:
C     -----------
C     STIFC(1:NOFVAR,1:NOFVAR,1:NOFVERT,1:NOFVERT) 
C                   convection matrix
C                   in the NOFVERT-1 vertices of the boundary face
C     NODRES(1:NOFVAR,1:NOFVERT-1) nodal residual due to the incoming
C                   characteristics in the NOFVERT-1 vertices 
C                   of the boundary face
C
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION BETA,CNST,TMP
      INTEGER I,IADD,IFAIL,IVERT,J,JVERT,K,L,NOFEQN
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION DZDU(LENB),FLUX(NMAX*VMAX)
      DOUBLE PRECISION WKSP(MAXNOFVAR) 
C     ..
C     .. External Subroutines ..
      EXTERNAL DAXPY,DGEMM,DINIT,GETDF4CORRDU,INVWLL,MATDZDU
      EXTERNAL GETDF4CORRDU4Ar,INVWLL4Ar,MATDZDU4Ar 
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC REAL
C     ..
C     .. Data statements ..

      DATA DZDU/LENB*ZERO/
C
C     ..
C     Compute correction flux
C
      NOFEQN = NDIM+NSP+1
C
      DO 3 IVERT = 1,NOFVERT
          IADD = (IVERT-1)*NOFVAR + 1
          CALL INVWLL4Ar(NDIM,VCN,VCZ(1,IVERT),FLUX(IADD))
    3 CONTINUE
!     CALL SIMPSON(NDIM,NOFVERT,NOFVAR,VCZ,VCN,FLUX,INVWLL4Ar) !< da cambiare INVWLL
C
      CALL DINIT(NOFVAR* (NOFVERT-1),ZERO,NODRES,1)
      BETA = (1.-ALPHA)/REAL(NDIM-1)/REAL(NDIM)
      DO 5 IVERT = 1,NOFVERT - 1
          DO 7 JVERT = 1,NOFVERT - 1
              IADD = (JVERT-1)*NOFVAR + 1
              IF (JVERT.EQ.IVERT) THEN
                  CNST = ALPHA/NDIM

              ELSE
                  CNST = BETA
              ENDIF
!             IADD = (NOFVERT-1)*NOFVAR + 1
!             CNST = ONE/REAL(NOFVERT-1)

              CALL DAXPY(NOFEQN,CNST,FLUX(IADD),1,NODRES(1,IVERT),1)
    7     CONTINUE
    5 CONTINUE
C
      IF (.NOT.PICARD) RETURN
C
C     Compute matrices dFcorr/dZ for all NOFVERT-1 vertices
C     of the boundary face
C
      DO 1 IVERT = 1,NDIM
          CALL GETDF4CORRDU4Ar(VCZ(1,IVERT),VCN,NDIM,NOFEQN, !< da cambiare INVWLL
     +                      WORK(1,1,IVERT))
    1 CONTINUE      
      
      BETA = (1.d0-ALPHA)/REAL(NDIM-1)
      DO 8 I = 1,NDIM
          DO 8 J = 1,NDIM
              IF (J.EQ.I) THEN
                  TMP = ALPHA/REAL(NDIM)

              ELSE
                  TMP = BETA/REAL(NDIM)
              ENDIF

              DO 8 L = 1,NOFEQN
                  DO 8 K = 1,NOFEQN
                     WORK2(K,L,I,J) = TMP*WORK(K,L,J)
    8 CONTINUE
!     do j = 1,ndim 
!     write(6,*)'vertice ',J
!      CALL R8Mat_Print('G',' ',Nofvar,Nofvar,WORK(1,1,J),Nofvar,
!    +      'WORK ',IFAIL)
!     enddo
!     pause
C
C
C     transform the convection stiffness matrix into
C     conserved variables as C_{ij} := 2 C_{ij} dZdU(j)
C
      DO 9 IVERT = 1,NDIM
          IADD = (IVERT-1)*NOFEQN*NOFEQN + 1
          CALL MATDZDU4Ar(VCZ(1,IVERT),DZDU(IADD),NDIM,NOFEQN) ! <-- da cambiare
!      CALL R8Mat_Print('G',' ',Nofvar,Nofvar,DZDU(IADD),Nofvar,
!    +      'dZdU ',IFAIL)
          DO 9 I = 1,NDIM
              CALL DGEMM('No transpose','No transpose',NOFEQN,NOFEQN,
     +                   NOFEQN,ONE,WORK2(1,1,I,IVERT),NOFEQN,
     +                   DZDU(IADD),NOFEQN,ZERO,STIFC(1,1,I,IVERT),
     +                   NOFVAR)
C
!      write(6,*)ivert,i
!      CALL R8Mat_Print('G',' ',Nofvar,Nofvar,STIFC(1,1,I,IVERT),Nofvar,
!    +      'C(i,j) ',IFAIL)
!      pause
C
    9 CONTINUE
C
!      write(6,*) DZDU
!      pause
C
      RETURN

  564 FORMAT ((E12.6,1X))

      END
