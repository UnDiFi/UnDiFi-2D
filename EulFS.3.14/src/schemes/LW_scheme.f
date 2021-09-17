!> \par Purpose
!>
!> Computes the UCV scheme for scalar problems on a triangle/tetrahedron. 
!>
!> The UCV scheme is attributed to [Giles et al. Upwind control volumes - A new upwind approach, AIAA 90-0104] (http://arc.aiaa.org/doi/abs/10.2514/6.1990-104)
!>
!> Its Fluctuation Splitting version is attributed to Paillere and described in [De Palma et al. Journal of Computational Physics 208 (2005) 1–3](http://dx.doi.org/doi:10.1016/j.jcp.2004.11.023)
!>
!> The UCV scheme is in fact the Lax Wendroff scheme with a particular choice
!> of the elemental time step; the distribution weights are:
!>
!> \f[
!> \beta_i^{UCV} = \left( \frac{1}{d+1} + \frac{2}{3} \frac{ k_{i} }{ \sum_{\ell=1}^{d+1}  |k_{\ell}| }  \right) =
!> \left( \frac{1}{d+1} + \frac{2}{3} \frac{ k_{i} }{ \sum_{\ell=1}^{d+1}  \left( k_{\ell}^{+} - k_{\ell}^{-} \right) } \right)
!> \f]
!> 
!>
!> @param[in] IELEM the current simplicial element
!> @param[in] VCN the \c NDIM cartesian component of the inward face normal to all \c NOFVERT vertices, scaled by its measure
!> @param[in] ADVECTION the \c NDIM cartesian component of the advection speed
!> @param[out] CELRES the elemental residual
!> @param[in] SOURCE the volume integral of the elemental source term
!> @param[in] Q the \c NOFVERT values of the dependent variable
!> @param[in,out] DT elemental contribution to the (inverse of the) time step
!> @param[in,out] NODRES is updated by summing the signals sent to each of the \c NOFVERT vertices of cell \c IELEM
!> @param[out] BETA the distribution matrices
!> @param[out] STIFC the elemental contribution to the implicit matrix (computed only if MATRIX_ASSEMBLY is true): \f$C_{ij} = -\beta_i k_j\f$.
!> @param[in] NDIM is the dimension of the space
!> @param[out] NOFVERT is number of vertices of the simplicial elements (= NDIM+1)
!> @param[in] MATRIX_ASSEMBLY when set == .TRUE. the \c STIFC matrix will be assembled
!> \warning Why 2/3 here? Check with Gee's report: it might not be ok in 3D
!> \author $Author: abonfi $
!> \version $Revision: 1.5 $
!> \date $Date: 2013/09/20 11:12:46 $
!>
      SUBROUTINE LW_SCHEME(IELEM,VCN,ADVECTION,CELRES,SOURCE,Q,DT,
     +                     NODRES,BETA,STIFC,NDIM,NOFVERT,
     3                     MATRIX_ASSEMBLY)
C
C     $Id: LW_scheme.f,v 1.5 2013/09/20 11:12:46 abonfi Exp $
C
C     .. Parameters ..
      INCLUDE 'paramt.h'
      INCLUDE 'constants.h'
      DOUBLE PRECISION CFLELEM
      PARAMETER (CFLELEM=TWO/3.d0)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION CELRES,SOURCE
      INTEGER IELEM,NDIM,NOFVERT
      LOGICAL MATRIX_ASSEMBLY
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION ADVECTION(NDIM),DT(NOFVERT),NODRES(NOFVERT),
     +                 Q(*),BETA(NOFVERT),STIFC(NOFVERT,NOFVERT),
     3                 VCN(NDIM,NOFVERT)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION KNEGINV,S,HELP,KSUM
      INTEGER I,IROW,IVERT,J,JCOL,NEGI,POSI
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION K(MAXNOFVERT)
C     ..
C     .. External Functions ..
      DOUBLE PRECISION DDOT
      EXTERNAL DDOT
C     ..
      CELRES = ZERO
      KSUM = ZERO
C
      DO 10 IVERT = 1,NOFVERT
C
C Dotting advection speed with the face normal
C
          HELP = DDOT(NDIM,VCN(1,IVERT),1,ADVECTION,1)
          K(IVERT) = HELP/NDIM
C
          CELRES = CELRES + Q(IVERT)*K(IVERT)
C
          KSUM = KSUM + ABS(K(IVERT))
C
   10 CONTINUE
C
      CELRES = CELRES + SOURCE
C
C Loops over all nodes
C
      HELP = ONE/REAL(NOFVERT)
      KSUM = CFLELEM/KSUM
      DO 20 I = 1,NOFVERT
          IF(K(I).GT.ZERO)DT(I) = DT(I) + K(I)
          BETA(I) = HELP + KSUM*K(I)
          NODRES(I) = -BETA(I)*CELRES
   20 CONTINUE
C
C
      IF (MATRIX_ASSEMBLY) THEN
          DO 40 J = 1,NOFVERT
              DO 40 I = 1,NOFVERT
                  STIFC(I,J) = -BETA(I)*K(J)
   40     CONTINUE
      ENDIF
C
      RETURN

      END
