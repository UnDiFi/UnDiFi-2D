!> \brief Computes the N scheme for scalar problems on a triangle/tetrahedron
!> \f[
!> \phi_i^{N} = - \left( \sum_{\ell=1}^{d+1} k_{\ell}^+ \right)
!! \delta_i^+ \left( u_i - u_- \right) = - k_i^+ \left( u_i - u_- \right)
!> \f]
!> 
!> \copydetails LDA_SCHEME()
      SUBROUTINE N_SCHEME(IELEM,VCN,ADVECTION,CELRES,SOURCE,Q,DT,NODRES,
     &                    BETA,STIFC,NDIM,NOFVERT,MATRIX_ASSEMBLY)
C
C
C This routine computes the N scheme on one tetrahedron
C
C this is a FORTRAN implementation of the original
C C version by G. Bourgois
C
C
      IMPLICIT NONE
      INCLUDE 'paramt.h'
      INCLUDE 'constants.h'
C
C
C
C     .. Parameters ..
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION CELRES,SOURCE
      INTEGER IELEM,NDIM,NOFVERT
      LOGICAL MATRIX_ASSEMBLY
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION ADVECTION(NDIM),BETA(NOFVERT),STIFC(NOFVERT,NOFVE
     &RT),DT(NOFVERT),NODRES(NOFVERT),Q(NOFVERT),VCN(NDIM,NOFVERT)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION KNEGINV,KPOS,S,UIN
      INTEGER I,IFAIL,IROW,J,JCOL,NEGI,POSI
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION K(MAXNOFVERT)
      INTEGER POS(MAXNOFVERT)
C     ..
C     .. External Functions ..
      DOUBLE PRECISION DDOT
      EXTERNAL DDOT
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC DABS
C     ..
C
      IF (SOURCE.NE.0.D0) then
                 write(6,*)source
                 CALL SETERR(
     +           31HN SCHEME - NON ZERO SOURCE TERM 
     +                          ,31,999,2)
      endif
C
      POSI = 0
      NEGI = NOFVERT + 1
      CELRES = ZERO
      KPOS = ZERO
      UIN = ZERO
C
C Compute the advection vector
C
      DO 10 I = 1,NOFVERT
C
C Dotting the advection speed with normal
C
          K(I) = DDOT(NDIM,VCN(1,I),1,ADVECTION,1)/NDIM
c
          CELRES = CELRES + Q(I)*K(I)
c
          NODRES(I) = ZERO
c
          IF (K(I).GT.ZERO) THEN
              POSI = POSI + 1
              POS(POSI) = I
              KPOS = KPOS + K(I)

          ELSE
              NEGI = NEGI - 1
              POS(NEGI) = I
              UIN = UIN - K(I)*Q(I)
          ENDIF

   10 CONTINUE
C
      IF (DABS(KPOS).LE.1.D-15) RETURN
      UIN = UIN/KPOS
C
C Loops over downstream nodes
C
      DO 20 I = 1,POSI
          J = POS(I)
          S = K(J)* (Q(J)-UIN)
          NODRES(J) = -S
          DT(J) = DT(J) + K(J)
   20 CONTINUE
C
!     if( posi .EQ. 2 )then
!         s = - nodres(pos(2))/(nodres(pos(1)) +1.d-15)
!     else 
!         s = 0.d0
!     endif
!     write(12,*)ielem,s
C
      IF (.NOT.MATRIX_ASSEMBLY) RETURN
      KNEGINV = -ONE/KPOS
C
C     The convection matrix has to be zeroth since in the
C     subsequent loops (28,30) not all vertices are touched
C
      DO 40 J = 1,NOFVERT
          DO 40 I = 1,NOFVERT
              STIFC(I,J) = ZERO
   40 CONTINUE
      DO 30 I = 1,POSI
          IROW = POS(I)
          S = K(IROW)
          DO 28 J = NOFVERT,NEGI,-1
              JCOL = POS(J)
              STIFC(IROW,JCOL) = S*K(JCOL)*KNEGINV
   28     CONTINUE
          STIFC(IROW,IROW) = -S
   30 CONTINUE
C
      RETURN

      END
