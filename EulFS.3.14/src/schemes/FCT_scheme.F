C
!> This routine computes the FCT scheme on one tetrahedron
C
!> The blending is the one suggested by
!! Sermeus, Van Ransbeek, Deconinck AMIF 98
!>
!> Only linear timestepping has been implemented
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
      SUBROUTINE FCT_SCHEME(IELEM,VCN,ADVECTION,CELRES,SOURCE,Q,DT,
     &                      NODRES,BETA,STIFC,NDIM,NOFVERT,
     &                      MATRIX_ASSEMBLY)
C
      implicit none
c
C
C
C
C
C
C     .. Parameters ..
      INCLUDE "constants.h"
      INCLUDE "paramt.h"
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION CELRES,SOURCE
      INTEGER IELEM,NDIM,NOFVERT
      LOGICAL MATRIX_ASSEMBLY
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION ADVECTION(NDIM),BETA(NOFVERT),
     &                 STIFC(NOFVERT,NOFVERT),DT(NOFVERT),
     &                 NODRES(NOFVERT),Q(NOFVERT),
     &                 VCN(NDIM,NOFVERT)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION ALPH,DENOM,KPOS,S,DUM,TEMPA,TEMPB,UIN,UOUT
      DOUBLE PRECISION TEMP,TOLER
      PARAMETER(TOLER=-5.D-14)
      INTEGER I,IVAR,IVERT,J,JVERT,POSI,IROW,JCOL
      logical lflag
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION K(MAXNOFVERT),DELPOS(MAXNOFVERT),
     2                               DELNEG(MAXNOFVERT)
      DOUBLE PRECISION work(MAXNOFVERT),cijn,cijlda
      INTEGER POS(MAXNOFVERT)
C     ..
C     .. External Functions ..
      DOUBLE PRECISION DDOT,DNRM2
      EXTERNAL DDOT,DNRM2
C     ..
C     .. External Subroutines ..
      EXTERNAL SETERR
C     ..
C
C     .. Intrinsic Functions ..
      INTRINSIC ABS
C     ..
      IF (SOURCE.NE.ZERO) CALL SETERR(
     +                              33HFCT SCHEME - NON ZERO SOURCE TERM
     +                                ,33,995,2)
C
      POSI = 0
      CELRES = ZERO
      KPOS = ZERO
      UIN = ZERO
      UOUT = ZERO
C
C Loops over all vertices
C
      DO 10 IVERT = 1,NOFVERT
c
c Dotting advection speed with normal
c
          DUM = DDOT(NDIM,VCN(1,IVERT),1,ADVECTION,1)
          K(IVERT) = DUM/NDIM
c
          TEMP = Q(IVERT)*K(IVERT)
          CELRES = CELRES + TEMP
c
          IF (K(IVERT).GT.ZERO) THEN
              POSI = POSI + 1
              POS(POSI) = IVERT
              KPOS = KPOS + K(IVERT)
              UOUT = UOUT + TEMP

          ELSE
              UIN = UIN - TEMP
              NODRES(IVERT) = ZERO
          ENDIF
c
   10 CONTINUE
C
      IF (POSI.EQ.0) RETURN
C
      UIN = UIN/KPOS
      UOUT = UOUT/KPOS
C
C computation of alpha
C
c ... Looping over downstream nodes
C
      DENOM = ZERO
      DO 32 I = 1,POSI
          J = POS(I)
          DENOM = DENOM + ABS( K(J) * (Q(J)-UIN) )
c
c
   32 CONTINUE
      IF(POSI.EQ.1)THEN
         ALPH = ONE
      ELSE ! deals with more than one target
         DENOM = MAX(DENOM,1.D-14) ! this might be useless
         ALPH = ONE - ABS(CELRES)/DENOM
         IF( (ALPH .LT. ZERO) )THEN
!        IF( (ALPH .LT. ZERO) .AND. (ALPH .GT. TOLER) )THEN
             ALPH = ZERO
         ENDIF
      ENDIF
      write(12,*)ielem,posi,ALPH
!     write(12,*)ielem,posi,ALPH,celres,(q(j),j=1,nofvert),denom,uin,
!    &uout
      if(alph.lt.0.d0.or.alph.gt.1.d0)then
            write(6,*)'problem with alpha ',ielem,alph,celres,denom
c     write(6,*)ielem,alph
c           pause
      endif
C
C Target tracking
C
      DUM = ZERO
c ... Looping over downstream nodes
      DO 30 I = 1,POSI
          J = POS(I)
          TEMPA = -K(J)* (Q(J)-UIN)
          TEMPB =  K(J)* (Q(J)-UOUT)
c
          NODRES(J) = TEMPA + ALPH*TEMPB
c
c Linear time-step for all the downstream nodes
c
          DT(J) = DT(J) + K(J)
   30 CONTINUE
C
C ---------- 
C
      IF (.NOT.MATRIX_ASSEMBLY) RETURN
C
          DO 40 J = 1,NOFVERT
              DELPOS(J) = MAX(ZERO,K(J))/KPOS
              DELNEG(J) =-MIN(ZERO,K(J))/KPOS
              DO 40 I = 1,NOFVERT
                  STIFC(I,J) = ZERO
   40     CONTINUE
          DO 31 I = 1,POSI
              IROW = POS(I)
              S = K(IROW)
              DO 28 JCOL = 1,NOFVERT
              if(jcol.ne.irow)then
caldo         cijn = s * delneg(jcol)
caldo         cijlda = - s / kpos * k(jcol)
                  STIFC(IROW,JCOL)= S*(DELNEG(JCOL)-ALPH*DELPOS(JCOL))
caldo             STIFC(IROW,JCOL)= cijn * (one-alph) + alph * cijlda
              else
                  STIFC(IROW,IROW)=-S*(ONE+ALPH*(DELPOS(IROW)-ONE))
              endif
   28         CONTINUE
   31     CONTINUE
C
C
#if 0
C
C     This is used to check that the convection matrix is all right
C
      write(6,*)'hello world!'
      lflag=.false.
      do 34 i = 1,nofvert
      DUM=zero 
      do 33 j = 1,nofvert
         DUM = DUM+stifc(i,j)*q(j)
   33 continue
      work(i) = DUM
      if(abs(DUM-nodres(i)).gt.1.d-14)lflag=.true.
   34 continue
      if(lflag)then
           write(6,*)ielem
        do 35 i = 1,nofvert
           write(6,*)nodres(i),work(i)
   35   continue
      pause
!     else
!       do i = 1,nofvert
!          write(6,*)ielem,nodres(i),work(i)
!       enddo
      endif
#endif

C
      RETURN

      END
