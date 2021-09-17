!> \copydetails LDA_SCHEME()
      SUBROUTINE SUPG_SCHEME(IELEM,VCN,ADVECTION,CELRES,SOURCE,Q,DT,
     2                       NODRES,BETA,STIFC,NDIM,NOFVERT,
     3                       MATRIX_ASSEMBLY)
C
C     $Id: SUPG_scheme.f,v 1.7 2013/01/24 07:46:33 abonfi Exp $
C
      IMPLICIT NONE 
C
C     ..
C     .. Parameters ..
C     ..
      INCLUDE 'paramt.h'
      INCLUDE 'constants.h'
C
C     SUPG type scheme with distribution coefficients:
C
C          beta_i = 1/(d+1) + TAU * k_i / V
C          TAU = C_1 * h / |a|
C          a is the convection speed
C          take h = |x_+ - x_-| = (|a| * V) / \sum k^+
C          C_1 = d/(d+1) or 2.
C
C
C
C ugly,ugly,ugly
C
      INCLUDE 'visco.com'
C     ..
C     .. Scalars in Common .. common ABC is also in scalar/scalar.f
      DOUBLE PRECISION VOLUME
      COMMON /ABC/VOLUME
C
C ugly,ugly,ugly
C
C
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION CELRES,SOURCE
      INTEGER IELEM,NDIM,NOFVERT
      LOGICAL MATRIX_ASSEMBLY
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION ADVECTION(NDIM),STIFC(NOFVERT,NOFVERT),
     +                 DT(NOFVERT),NODRES(NOFVERT),Q(NOFVERT),
     +                 VCN(NDIM,NOFVERT),BETA(NOFVERT)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION C_1,GALERKIN,KPOS,S,TAU,Z,PEH,XI
      PARAMETER(C_1=HALF)
      INTEGER I,IVERT,J
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION K(MAXNOFVERT)
C     ..
C     .. External Functions ..
      DOUBLE PRECISION DDOT
      EXTERNAL DDOT
C     ..
C     .. External Subroutines ..
C     EXTERNAL MYFUN
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC MAX,REAL,DTANH
C     ..
C     .. Statement Functions ..
      DOUBLE PRECISION COTH
C     COTH(Z)=(EXP(2.d0*Z)+1.d0)/(EXP(2.d0*Z)-1.d0)
      COTH(Z)=ONE/DTANH(Z)
C
      CELRES = ZERO
      KPOS = ZERO
C
      DO 10 IVERT = 1,NOFVERT
c
c Dotting advection speed with the face normal
c
          K(IVERT) = DDOT(NDIM,VCN(1,IVERT),1,ADVECTION,1)/NDIM
c
          CELRES = CELRES + Q(IVERT)*K(IVERT)
c
          IF (K(IVERT).GT.ZERO) KPOS = KPOS + K(IVERT)
   10 CONTINUE
C
      CELRES = CELRES
C
      GALERKIN = ONE/REAL(NDIM+1)
C
!     PEH = DDOT(NDIM,ADVECTION,1,ADVECTION,1)*VOLUME/(2.0d0*REINV*KPOS)
      PEH = DDOT(NDIM,ADVECTION,1,ADVECTION,1)*VOLUME/(REINV*KPOS)
!     XI = COTH(PEH)-1.d0/PEH
c
c     this is the optimal design, see e.g. Carette PhD (2.30) 
c
      XI = COTH(HALF*PEH)-TWO/PEH
C     XI = 1.d0 ! advection dominated limit
      TAU = C_1*XI/KPOS
C
C Loops over downstream nodes
C

      DO 20 I = 1,NOFVERT
          BETA(I) = GALERKIN + TAU*K(I)
          S = BETA(I)*CELRES
   22     NODRES(I) = -S
          DT(I) = DT(I) + MAX(K(I),ZERO)
   20 CONTINUE
C
C     CALL FSOU(NDIM,NOFVERT,TAU,K,NODRES,VCP,myfun)
C
      IF (MATRIX_ASSEMBLY) THEN
          DO 30 I = 1,NOFVERT
              S = -BETA(I)
              DO 28 J = 1,NOFVERT
                  STIFC(I,J) = S*K(J)
   28         CONTINUE
   30     CONTINUE
      ENDIF
C
      RETURN

      END
