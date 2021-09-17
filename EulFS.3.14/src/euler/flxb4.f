
      SUBROUTINE FLXB4(NDIM,NOFVAR,NOFVERT,STIFC,BFLX,WORK,VCZ,VCB,
     +                 VCN,NODRES,PICARD)
C
      IMPLICIT NONE
C
C
C     Purpose: 
C     ------
C     compute prescribed flux b.c.'s for compressible flows
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
      DOUBLE PRECISION STIFC(NOFVAR,NOFVAR,NOFVERT,NOFVERT),
     &                 BFLX(NOFVAR),WORK(*),VCZ(NOFVAR,NOFVERT),
     &                 VCB(NDIM),VCN(NDIM),NODRES(NOFVAR,NOFVERT)
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
      INTEGER I,IADDR,IFAIL,IVERT,J,JVERT,K,L,NOFEQN
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION FLUX(MAXNOFVAR*MAXNOFVERT),WKSP(5)
C     ..
C     .. External Subroutines ..
      EXTERNAL DAXPY,DGEMM,DINIT,GETDF4CORRDU,INVWLL,FLUX4
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC REAL
C     ..
C     .. Data statements ..

C     ..
C
C     Compute correction flux
C
      NOFEQN = NDIM+2
C
caldo write(6,*)'bndry flux ',(bflx(i),i=1,nofeqn)
C
      CALL SIMPSON(NDIM,NOFVERT,NOFVAR,VCZ,VCB,VCN,FLUX,FLUX4)
C
C     the flux through the face is returned in FLUX(1,NOFVERT)
C
      IADDR = (NOFVERT-1)*NOFVAR+1
      CALL DSCAL(NOFEQN,MONE,FLUX(IADDR),1)
!     CALL R8Mat_Print('G',' ',Nofvar,Nofvert,flux,nofvar,
!    +      'Flux before subtracting ',IFAIL)
      CALL DAXPY(NOFEQN,ONE,BFLX,1,FLUX(IADDR),1)
!     CALL R8Mat_Print('G',' ',Nofvar,Nofvert,flux,nofvar,
!    +      'Flux  after subtracting',IFAIL)
C
      CALL DINIT(NOFVAR* (NOFVERT-1),ZERO,NODRES,1)
!     BETA = (1.-ALPHA)/REAL(NDIM-1)/REAL(NDIM)
      BETA = ONE/REAL(NOFVERT-1)
      DO 5 IVERT = 1,NOFVERT - 1
           CALL DAXPY(NOFEQN,BETA,FLUX(IADDR),1,NODRES(1,IVERT),1)
    5 CONTINUE
C
!     write(6,*)
!     write(6,*)(flux(i),i=1,nofeqn)
!     write(6,*)(bflx(i),i=1,nofeqn)
!     CALL R8Mat_Print('G',' ',Nofvar,Nofvar,STIFC(1,1,I,IVERT),Nofvar,
!    +      'C(i,j) ',IFAIL)
!     pause
C
C
      IF (.FALSE.) THEN
!     IF (PICARD) THEN
          WRITE(6,*) 'Unimplemented feature in subr. flxb4'
          CALL EXIT(3)
      ENDIF
      RETURN

  564 FORMAT ((E12.6,1X))

      END
