!> \copydetails LDA_SCHEME()
      SUBROUTINE PSI_SCHEME(IELEM,VCN,ADVECTION,CELRES,SOURCE,Q,DT,
     &                      NODRES,BETA,STIFC,NDIM,NOFVERT,
     &                      MATRIX_ASSEMBLY)
C
C     $Id: PSI_scheme.f,v 1.4 2013/01/24 07:46:33 abonfi Exp $
C
C This routine computes the PSI scheme on one tetrahedron
C
C this is a FORTRAN implementation of the original
C C version by G. Bourgois
C
c Only linear timestepping has been implemented
c
C
      IMPLICIT NONE
C
C
C
C
C
C     .. Parameters ..
      INCLUDE "paramt.h"
      INCLUDE "constants.h"
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION CELRES,SOURCE
      INTEGER IELEM,NDIM,NOFVERT
      LOGICAL MATRIX_ASSEMBLY
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION ADVECTION(NDIM),BETA(NOFVERT),
     2                 STIFC(NOFVERT,NOFVERT),
     &                 DT(NOFVERT),NODRES(NOFVERT),Q(NOFVERT),
     &                 VCN(NDIM,NOFVERT)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION KPOS,S,HELP,UIN
      INTEGER I,IVERT,J,POSI,TAGI
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION DQ(MAXNOFVERT),K(MAXNOFVERT)
      INTEGER POS(MAXNOFVERT),TAG(MAXNOFVERT)
C     ..
C     .. External Functions ..
      DOUBLE PRECISION DDOT
      EXTERNAL DDOT
C     ..
C     .. External Subroutines ..
      EXTERNAL SETERR
C     ..
C
      IF (SOURCE.NE.ZERO) CALL SETERR(
     +           33HPSI SCHEME - NON ZERO SOURCE TERM,33,995,2)
C
      POSI = 0
      TAGI = 0
      CELRES = ZERO
      KPOS = ZERO
      UIN = ZERO
C
C Loops over all vertices
C
      DO 10 IVERT = 1,NOFVERT
c
c Dotting advection speed with normal
c
          HELP = DDOT(NDIM,VCN(1,IVERT),1,ADVECTION,1)
          K(IVERT) = HELP/NDIM
c
          CELRES = CELRES + Q(IVERT)*K(IVERT)
c
          IF (K(IVERT).GT.ZERO) THEN
              POSI = POSI + 1
              POS(POSI) = IVERT
              KPOS = KPOS + K(IVERT)

          ELSE
              UIN = UIN - K(IVERT)*Q(IVERT)
              NODRES(IVERT) = ZERO
          ENDIF
c
   10 CONTINUE
C
      IF (POSI.EQ.0) RETURN
C
      UIN = UIN/KPOS
C
C Target tracking
C
      HELP = ZERO
c ... Looping over downstream nodes
      DO 30 I = 1,POSI
          J = POS(I)
          DQ(J) = Q(J) - UIN
c
          IF (DQ(J)*CELRES.GT.ZERO) THEN
              TAGI = TAGI + 1
              TAG(TAGI) = J
              HELP = HELP + DQ(J)*K(J)
          ENDIF
c
c Linear time-step for all the downstream nodes
c
          DT(J) = DT(J) + K(J)
   30 CONTINUE
C
C Loops over downstream nodes (TARGET UPDATING)
C
      DO 20 I = 1,TAGI
          J = TAG(I)
          S = CELRES/HELP*DQ(J)*K(J)
c
          NODRES(J) = -S
c
   20 CONTINUE
C
C ---------- Distribution of the source term using LDA ----------
C
C     CELRES = CELRES + SOURCE
C ...  Looping over downstream nodes
C     DO 35 I = 1,POSI
C         J = POS(I)
C
C         S = K(J)*SOURCE/KPOS
C
C         NODRES(J) = -S
C
C  35 CONTINUE
C  37 CONTINUE
      IF (MATRIX_ASSEMBLY) CALL SETERR(
     +           52HPSI SCHEME IMPLICIT TIME INTEGRATION NOT IMPLEMENTED
     +                          ,52,999,2)
C
      RETURN

      END
