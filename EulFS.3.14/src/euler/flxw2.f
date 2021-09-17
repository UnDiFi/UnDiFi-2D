!> \copydetails FLXW4()
      SUBROUTINE FLXW2(NDIM,NOFVAR,NOFVERT,STIFC,WORK,WORK2,VCZ,VCB,VCN,
     +                 NODRES,PICARD)
C
C     $Id: flxw2.f,v 1.4 2012/12/20 10:39:57 abonfi Exp $
C
      IMPLICIT NONE
C
C     Purpose: 
C     ------
C     compute inviscid wall b.c.'s for INcompressible flows
C
      include 'paramt.h'
      include 'constants.h'
C
C
C     .. Parameters ..
      DOUBLE PRECISION ALPHA
      PARAMETER (ALPHA=0.75d0)
C     ..
C     .. Scalar Arguments ..
      INTEGER NDIM,NOFVAR,NOFVERT
      LOGICAL PICARD
C
C     NDIM   dimension of the space (2 or 3)
C     NOFVAR number of variables (degrees of freedom)
C            in each meshpoint
C     NOFVAR number of variables (degrees of freedom)
C            in each meshpoint
C     PICARD .TRUE. for Picard linearisation
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION NODRES(NOFVAR,NDIM),
     +                 STIFC(NOFVAR,NOFVAR,NOFVERT,NOFVERT),VCN(NDIM),
     +                 VCZ(NOFVAR,NDIM),WORK(NDIM+1,NDIM+1,NOFVERT-1),
     +                 WORK2(NDIM+1,NDIM+1,NOFVERT-1,NOFVERT-1),
     &                 VCB(NDIM,NOFVERT)
C
C
C     On entry:
C     --------
C     VCN(1:NDIM,1:NOFVERT) cartesian components of the normals 
C                           to the element sides/faces
C                           the normal to the boundary face need
C                           to be stored in VCN(1:NDIM,NOFVERT)
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
      DOUBLE PRECISION BETA,CNST
      INTEGER I,IADD,IFAIL,IVERT,J,JVERT,K,L,NORDER
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION FLUX(MAXNOFVERT*MAXNOFVAR)
C     ..
C     .. External Subroutines ..
      EXTERNAL DAXPY,DINIT,GETDF2CORRDU,INVWLLI
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC REAL
C     ..
      NORDER = NDIM + 1
      DO 2 IVERT = 1,NOFVERT
          IADD = (IVERT-1)*NOFVAR + 1
          CALL INVWLLI(NDIM,VCN,VCB(1,IVERT),VCZ(1,IVERT),FLUX(IADD))
    2 CONTINUE
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

              CALL DAXPY(NORDER,CNST,FLUX(IADD),1,NODRES(1,IVERT),1)
    7     CONTINUE
    5 CONTINUE
C
      IF (.NOT.PICARD) RETURN

      DO 3 IVERT = 1,NOFVERT - 1
          CALL GETDF2CORRDU(VCZ(2,IVERT),VCB(1,IVERT),VCN,NDIM,NORDER,
C                               ^
C                               | 
C                               | 
C            adresses the location of the x-component of the velocity vector
C
     +                      WORK(1,1,IVERT))
    3 CONTINUE

      DO 8 I = 1,NOFVERT - 1
          DO 8 J = 1,NOFVERT - 1
              IF (J.EQ.I) THEN
                  CNST = ALPHA/REAL(NDIM)

              ELSE
                  CNST = BETA
              ENDIF

              DO 8 L = 1,NORDER
                  DO 8 K = 1,NORDER
caldo             STIFC(k,l,i,j) = 0.5d0*CNST*work(k,l,j)
                      STIFC(K,L,I,J) = CNST*WORK(K,L,J)
    8 CONTINUE

      RETURN

  564 FORMAT ((E12.6,1X))

      END
