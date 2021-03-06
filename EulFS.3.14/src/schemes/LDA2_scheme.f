!> \copydoc LDA_SCHEME()
      SUBROUTINE LDA2_SCHEME(IELEM,VCN,ADVECTION,CELRES,SOURCE,Q,DT,
     +                       NODRES,BETA,STIFC,NDIM,NOFVERT,
     3                       MATRIX_ASSEMBLY)
C
C     $Id: LDA2_scheme.f,v 1.2 2013/01/24 07:46:33 abonfi Exp $
C
C SOURCE IS THE VOLUME INTEGRAL OF THE SOURCE TERM
C Hiro's versione of the LDA scheme, with a Peh dependence
C
C
      INCLUDE 'paramt.h'
      INCLUDE 'constants.h'
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
C
C     .. Parameters ..
      DOUBLE PRECISION C_1
      PARAMETER (C_1=ONE)
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
      DOUBLE PRECISION KNEGINV,KPOS,S,SM,PEH,AJ
      INTEGER I,IROW,IVERT,J,JCOL,NEGI,POSI
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION K(VMAX)
      INTEGER POS(VMAX)
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
          K(IVERT) = DDOT(NDIM,VCN(1,IVERT),1,ADVECTION,1)/NDIM
c
          CELRES = CELRES + Q(IVERT)*K(IVERT)
c
          IF (K(IVERT).GT.ZERO) THEN
              POSI = POSI + 1
              POS(POSI) = IVERT
              KPOS = KPOS + K(IVERT)

          ELSE
              NEGI = NEGI - 1
              POS(NEGI) = IVERT
              NODRES(IVERT)=ZERO
          ENDIF

   10 CONTINUE
C
      CELRES = CELRES + SOURCE
C
      IF(POSI.EQ.0)RETURN
C
!     PEH = DDOT(NDIM,ADVECTION,1,ADVECTION,1)*VOLUME/(2.0d0*REINV*KPOS)
      PEH = DDOT(NDIM,ADVECTION,1,ADVECTION,1)*VOLUME/(REINV*KPOS)
      DENO = 1.d0+C_1/PEH
C
C Loops over ALL nodes
C
      SM = 0.d0
      DO 20 J = 1,NOFVERT
          BETA(J) = ( MAX(0.d0,K(J)/KPOS) + (C_1/(NOFVERT*PEH)) )/DENO
!         sm = sm + aj
          S = BETA(J)*CELRES
          NODRES(J) = -S
          DT(J) = DT(J) + MAX(0.d0,K(J))
   20 CONTINUE
C
C         write(6,*)ielem,peh,sm
C
C
      IF (MATRIX_ASSEMBLY) THEN
          KNEGINV = -ONE/KPOS
          DO 40 J = 1,NOFVERT
              DO 40 I = 1,NOFVERT
                  STIFC(I,J) = ZERO
   40     CONTINUE
          DO 30 J = 1,NOFVERT
          DO 30 I = 1,NOFVERT
                  STIFC(I,J) = -BETA(I)*K(J)
   30     CONTINUE
      ENDIF
C
!     IF (LTIME) THEN
!         CALL UNSTEADY(BETA,Q,NODRES,STIFC,NDIM,NOFVERT,
!    &                  MATRIX_ASSEMBLY)
!     ENDIF ! LTIME
C
C
      RETURN

      END
