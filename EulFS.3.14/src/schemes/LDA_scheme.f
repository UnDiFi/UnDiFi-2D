!> \brief Computes the LDA scheme for scalar problems on a triangle/tetrahedron
!> \f[
!> \beta_i^{LDA} = \delta_i^+ = \frac{ k_{i}^{+} }{ \sum_{\ell=1}^{d+1}  k_{\ell}^{+} }
!> \f]
!> 
!>
!> @param[in] IELEM the current simplicial element
!> @param[in] VCN the NDIM cartesian component of the inward face normal to all NOFVERT vertices, scaled by its measure
!> @param[in] ADVECTION the NDIM cartesian component of the advection speed
!> @param[out] CELRES the elemental residual
!> @param[in] SOURCE the volume integral of the elemental source term
!> @param[in] Q the NOFVERT values of the dependent variable
!> @param[in,out] DT elemental contribution to the (inverse of the) time step
!> @param[out] NODRES the signals sent to each of the NOFVERT vertices
!> @param[out] BETA the distribution matrices
!> @param[out] STIFC the elemental contribution to the implicit matrix (computed only if MATRIX_ASSEMBLY is true)
!> @param[in] NDIM is the dimension of the space
!> @param[out] NOFVERT is number of vertices of the simplicial elements (= NDIM+1)
!> @param[in] MATRIX_ASSEMBLY when set == .TRUE. the STIFC matrix will be assembled
!>
      SUBROUTINE LDA_SCHEME(IELEM,VCN,ADVECTION,CELRES,SOURCE,Q,DT,
     +                      NODRES,BETA,STIFC,NDIM,NOFVERT,
     3                      MATRIX_ASSEMBLY)
C
C     $Id: LDA_scheme.f,v 1.6 2013/08/22 15:10:57 abonfi Exp $
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
              KPOS = KPOS + K(IVERT)
              BETA(IVERT) = K(IVERT)
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
      IF(POSI.EQ.0)THEN
          WRITE(I1MACH(4),*)'At least one downstream vertex expected ',P
     &OSI,IELEM
          WRITE(I1MACH(4),*)'n = ',((VCN(J,IVERT),J=1,NDIM),IVERT=1,NOFV
     &ERT)
          WRITE(I1MACH(4),*)'u = ',(ADVECTION(IVERT),IVERT=1,NDIM)
          WRITE(I1MACH(4),*)'k = ',(K(IVERT),IVERT=1,NOFVERT)
          RETURN
      ENDIF
C
C Loops over downstream nodes
C

      DO 20 I = 1,POSI
          J = POS(I)
          S = K(J)/KPOS*CELRES
          NODRES(J) = -S
          DT(J) = DT(J) + K(J)
          BETA(J) = BETA(J)/KPOS
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
              S = K(IROW)
              DO 28 JCOL = 1,NOFVERT
                  STIFC(IROW,JCOL) = S*K(JCOL)*KNEGINV
   28         CONTINUE
   30     CONTINUE
      ENDIF
C
      RETURN

      END
