!> \copydoc LDA_SCHEME()
      SUBROUTINE LDASqr_SCHEME(IELEM,VCN,ADVECTION,CELRES,SOURCE,Q,DT,
     +                      NODRES,BETA,STIFC,NDIM,NOFVERT,
     3                      MATRIX_ASSEMBLY)
C     $Id: LDASqr_scheme.f,v 1.2 2013/01/24 07:46:33 abonfi Exp $
C
C this is a FORTRAN implementation of the original
C C version by G. Bourgois
C
C SOURCE IS THE VOLUME INTEGRAL OF THE SOURCE TERM
C
C
      INCLUDE 'paramt.h'
      INCLUDE 'constants.h'
C
C     .. Parameters ..
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION CELRES,SOURCE
      INTEGER IELEM,NDIM,NOFVERT
      LOGICAL MATRIX_ASSEMBLY
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION ADVECTION(NDIM),DT(NOFVERT),NODRES(NOFVERT),
     +                 Q(NOFVERT),BETA(NOFVERT),STIFC(NOFVERT,NOFVERT),
     +                 VCN(NDIM,NOFVERT)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION KNEGINV,KPOS,S,SUM
      INTEGER I,IROW,IVERT,J,JCOL,NEGI,POSI
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION K(MAXNOFVERT)
      INTEGER POS(MAXNOFVERT)
C     ..
C     .. External Functions ..
      DOUBLE PRECISION DDOT
      EXTERNAL DDOT
C     ..
      POSI = 0
      NEGI = NOFVERT + 1
      CELRES = ZERO
      KPOS = ZERO
C
      DO 10 IVERT = 1,NOFVERT
c
c Dotting advection speed with the face normal
c
          SUM = DDOT(NDIM,VCN(1,IVERT),1,ADVECTION,1)
          K(IVERT) = SUM/NDIM
c
          CELRES = CELRES + Q(IVERT)*K(IVERT)
c
          IF (K(IVERT).GT.ZERO) THEN
              POSI = POSI + 1
              POS(POSI) = IVERT
              KPOS = KPOS + K(IVERT)**2
              BETA(IVERT) = K(IVERT)**2
          ELSE
              NEGI = NEGI - 1
              POS(NEGI) = IVERT
              NODRES(IVERT)=ZERO
              BETA(IVERT) = ZERO
          ENDIF

   10 CONTINUE
C
      CELRES = CELRES + SOURCE
C
      IF(POSI.EQ.0)RETURN
C
C Loops over downstream nodes
C

      DO 20 I = 1,POSI
          J = POS(I)
          DT(J) = DT(J) + BETA(J)
          BETA(J) = BETA(J)/KPOS
          NODRES(J) = -BETA(J)*CELRES
   20 CONTINUE
C
C
C
      IF (MATRIX_ASSEMBLY) THEN
          KNEGINV = -ONE/KPOS
          DO 40 J = 1,NOFVERT
              DO 40 I = 1,NOFVERT
                  STIFC(I,J) = ZERO
   40     CONTINUE
          DO 30 I = 1,POSI
              IROW = POS(I)
              S = K(IROW)**2
              DO 28 JCOL = 1,NOFVERT
                  STIFC(IROW,JCOL) = S*K(JCOL)*KNEGINV
   28         CONTINUE
   30     CONTINUE
      ENDIF
C
C
      RETURN

      END
