      SUBROUTINE CHECK(IELEM,NDIM,NVAR)
c
      IMPLICIT NONE
c
c This routine computes the flux integral over the cell # IELEM
c as : div F = dF/dZ*dZ/dx + dG/dZ*dZ/dx + dH/dZ*dZ/dx
c where Z is the parameter vector (for compressible)
C or the primitive variables (p,u,v,w) for incompressible
C
C it is meant to be used as a debugging tool
C
C vec_$dot(U,V,Count)			result = SUM(U(i) * V(i))
C vec_$mult_constant(U,Count,a,R)	R(i) = a * U(i)
C vec_$mult_add(U,V,Count,a,R)		R(i) = (a * V(i)) + U(i)
C vec_$mult_sub(U,V,Count,a,R)		R(i) = (a * V(i)) - U(i)
C
C     .. Parameters ..
C
      INCLUDE 'paramt.h'
      INCLUDE 'constants.h'
      INCLUDE 'plasma.h'

      INTEGER LDA
      PARAMETER(LDA=MAXNOFVAR)
C
C     .. Commons ..
C
      INCLUDE 'dofs.com'
      INCLUDE 'three.com'
      INCLUDE 'four.com'
      INCLUDE 'chorin.com'
      INCLUDE 'flags.com'
      INCLUDE 'pfcgas.com'
C
C     .. Scalar Arguments ..
C
      INTEGER IELEM,NDIM,NVAR
C
C     .. Local Scalars ..
C
      INTEGER IDIM,IFAIL
C
C     .. Local Scalars ..
C
       INTEGER I,J
C
C     .. Local Arrays ..
C
C
C     .. Local Arrays ..
C
      DOUBLE PRECISION Jacobian(LDA,LDA,3)
      DOUBLE PRECISION DPDZS(NSP),DPDZH,DPDZU(3),SUMDPZS,HELP
C
C     .. External Subroutines ..
C
      EXTERNAL DGEMV
C
C     .. Executable Statements ..
C
      IF( IABS(KAN) .EQ. 2 )THEN
        GOTO 100
      ELSEIF( IABS(KAN) .EQ. 4 )THEN
        GOTO 200
      ELSEIF( IABS(KAN) .EQ. 3 )THEN
        GOTO 400
      ENDIF
c
  100 CONTINUE
c
c     JACOBIAN IN Parameter Vector (Incompressible)
c
c	Matrix dF/dZ
c
      Jacobian(1,1,1) = ZERO
      Jacobian(1,2,1) = BETA
      Jacobian(1,3,1) = ZERO
      Jacobian(1,4,1) = ZERO
c
      Jacobian(2,1,1) = ONE
      Jacobian(2,2,1) = 2*ZAVG(2)
      Jacobian(2,3,1) = ZERO
      Jacobian(2,3,1) = ZERO
c
      Jacobian(3,1,1) = ZERO
      Jacobian(3,2,1) = ZAVG(3)
      Jacobian(3,3,1) = ZAVG(2)
      Jacobian(3,4,1) = ZERO
c
      Jacobian(4,1,1) = ZERO
      Jacobian(4,2,1) = ZAVG(4)
      Jacobian(4,3,1) = ZERO
      Jacobian(4,4,1) = ZAVG(2)
c
c	MATRIX dG/dZ
c
      Jacobian(1,1,2) = ZERO
      Jacobian(1,2,2) = ZERO
      Jacobian(1,3,2) = BETA
      Jacobian(1,4,2) = ZERO
c
      Jacobian(2,1,2) = ZERO
      Jacobian(2,2,2) = ZAVG(3)
      Jacobian(2,3,2) = ZAVG(2)
      Jacobian(2,4,2) = ZERO
c
      Jacobian(3,1,2) = ONE
      Jacobian(3,2,2) = ZERO
      Jacobian(3,3,2) = TWO*ZAVG(3)
      Jacobian(3,4,2) = ZERO
c
      Jacobian(4,1,2) = ZERO
      Jacobian(4,2,2) = ZERO
      Jacobian(4,3,2) = ZAVG(4)
      Jacobian(4,4,2) = ZAVG(3)
c
c	Matrix dH/dZ
c
c
      Jacobian(1,1,3) = ZERO
      Jacobian(1,2,3) = ZERO
      Jacobian(1,3,3) = ZERO
      Jacobian(1,4,3) = BETA
c
      Jacobian(2,1,3) = ZERO
      Jacobian(2,2,3) = ZAVG(4)
      Jacobian(2,3,3) = ZERO
      Jacobian(2,4,3) = ZAVG(2)
c
      Jacobian(3,1,3) = ZERO
      Jacobian(3,2,3) = ZERO
      Jacobian(3,3,3) = ZAVG(4)
      Jacobian(3,4,3) = ZAVG(3)
c
      Jacobian(4,1,3) = ONE
      Jacobian(4,2,3) = ZERO
      Jacobian(4,3,3) = ZERO
      Jacobian(4,4,3) = TWO*ZAVG(4)
c
      GOTO 300
  200 CONTINUE
c
c     JACOBIAN IN Parameter Vector (Compressible)
c
c	Matrix dF/dZ
c
      Jacobian(1,1,1) = ZAVG(3)
      Jacobian(1,2,1) = ZERO
      Jacobian(1,3,1) = ZAVG(1)
      Jacobian(1,4,1) = ZERO
      Jacobian(1,5,1) = ZERO
c
      Jacobian(2,1,1) = ZERO
      Jacobian(2,2,1) = ZAVG(3)
      Jacobian(2,3,1) = ZAVG(2)
      Jacobian(2,4,1) = ZERO
      Jacobian(2,5,1) = ZERO
c
      Jacobian(3,1,1) = GM1OG * ZAVG(2)
      Jacobian(3,2,1) = GM1OG * ZAVG(1)
      Jacobian(3,3,1) = GP1OG * ZAVG(3)
      Jacobian(3,4,1) =-GM1OG * ZAVG(4)
      Jacobian(3,5,1) =-GM1OG * ZAVG(5)
c
      Jacobian(4,1,1) = ZERO
      Jacobian(4,2,1) = ZERO
      Jacobian(4,3,1) = ZAVG(4)
      Jacobian(4,4,1) = ZAVG(3)
      Jacobian(4,5,1) = ZERO
c
      Jacobian(5,1,1) = ZERO
      Jacobian(5,2,1) = ZERO
      Jacobian(5,3,1) = ZAVG(5)
      Jacobian(5,4,1) = ZERO
      Jacobian(5,5,1) = ZAVG(3)
c
c	MATRIX dG/dZ
c
      Jacobian(1,1,2) = ZAVG(4)
      Jacobian(1,2,2) = ZERO
      Jacobian(1,3,2) = ZERO
      Jacobian(1,4,2) = ZAVG(1)
      Jacobian(1,5,2) = ZERO
c
      Jacobian(2,1,2) = ZERO
      Jacobian(2,2,2) = ZAVG(4)
      Jacobian(2,3,2) = ZERO
      Jacobian(2,4,2) = ZAVG(2)
      Jacobian(2,5,2) = ZERO
c
      Jacobian(3,1,2) = ZERO
      Jacobian(3,2,2) = ZERO
      Jacobian(3,3,2) = ZAVG(4)
      Jacobian(3,4,2) = ZAVG(3)
      Jacobian(3,5,2) = ZERO
c
      Jacobian(4,1,2) = GM1OG * ZAVG(2)
      Jacobian(4,2,2) = GM1OG * ZAVG(1)
      Jacobian(4,3,2) =-GM1OG * ZAVG(3)
      Jacobian(4,4,2) = GP1OG * ZAVG(4)
      Jacobian(4,5,2) =-GM1OG * ZAVG(5)
c
      Jacobian(5,1,2) = ZERO
      Jacobian(5,2,2) = ZERO
      Jacobian(5,3,2) = ZERO
      Jacobian(5,4,2) = ZAVG(5)
      Jacobian(5,5,2) = ZAVG(4)
c
c	Matrix dH/dZ
c
      Jacobian(1,1,3) = ZAVG(5)
      Jacobian(1,2,3) = ZERO
      Jacobian(1,3,3) = ZERO
      Jacobian(1,4,3) = ZERO
      Jacobian(1,5,3) = ZAVG(1)
c
      Jacobian(2,1,3) = ZERO
      Jacobian(2,2,3) = ZAVG(5)
      Jacobian(2,3,3) = ZERO
      Jacobian(2,4,3) = ZERO
      Jacobian(2,5,3) = ZAVG(2)
c
      Jacobian(3,1,3) = ZERO
      Jacobian(3,2,3) = ZERO
      Jacobian(3,3,3) = ZAVG(5)
      Jacobian(3,4,3) = ZERO
      Jacobian(3,5,3) = ZAVG(3)
c
      Jacobian(4,1,3) = ZERO
      Jacobian(4,2,3) = ZERO
      Jacobian(4,3,3) = ZERO
      Jacobian(4,4,3) = ZAVG(5)
      Jacobian(4,5,3) = ZAVG(4)
c
      Jacobian(5,1,3) = GM1OG * ZAVG(2)
      Jacobian(5,2,3) = GM1OG * ZAVG(1)
      Jacobian(5,3,3) =-GM1OG * ZAVG(3)
      Jacobian(5,4,3) =-GM1OG * ZAVG(4)
      Jacobian(5,5,3) = GP1OG * ZAVG(5)
C
      GOTO 300
  400 CONTINUE
C
C     Sum of (dp/d(rhos)) * Zs          
C
      SUMDPZS = ZERO
      DO I = 1 , NSP
        SUMDPZS = SUMDPZS + DR(I) * ZAVG(I)
      ENDDO
C
      HELP = ONE / (ONE + DE) 
C
C     PRESSURE JACOBIAN  DP/DZ 
C
C     singles species
C
      DO I = 1 , NSP
         DPDZS(I) = HELP * (SUMDPZS + DR(I)*SQRTR + DE*ZAVG(IE)
     &   + DM(1)*ZAVG(IX) + DM(2)*ZAVG(IY))
         IF (NDIM.EQ.3) THEN
           DPDZS(I) = DPDZS(I) +  DM(3)*ZAVG(IZ) * HELP    
         ENDIF  
      ENDDO
!      HELP = SQRTR/(ONE+DE)
C
C     energy	
C
      DPDZH = DE*HELP*SQRTR
C
C     x-momentum
C
      DPDZU(1) = DM(1)*HELP*SQRTR
C
C     y-momentum
C      
      DPDZU(2) = DM(2)*HELP*SQRTR
C
C     z-momentum
C
      IF (NDIM.EQ.3) THEN
         DPDZU(3) = DM(3)*HELP*SQRTR
      ENDIF 
c
c     JACOBIAN IN Parameter Vector (Plasma)
c
c	Matrix dF/dZ
c
      DO I = 1 , NSP
        DO J = 1 , NSP  
          Jacobian(I,J,1) = ZERO
          IF (I.EQ.J) THEN
            Jacobian(I,J,1) = Jacobian(I,J,1) + ZAVG(IX)
          ENDIF
        ENDDO
      ENDDO
c      
      DO I = 1 , NSP
        Jacobian(I,IE,1) = ZERO
        Jacobian(I,IX,1) = ZAVG(I)
        Jacobian(I,IY,1) = ZERO
        Jacobian(I,IZ,1) = ZERO
      ENDDO
c     
      DO J = 1 , NSP
        Jacobian(IE,J,1) = ZERO
        Jacobian(IX,J,1) = DPDZS(J)
        Jacobian(IY,J,1) = ZERO
        Jacobian(IZ,J,1) = ZERO
      ENDDO
c
      Jacobian(IE,IE,1) = ZAVG(IX)
      Jacobian(IE,IX,1) = ZAVG(IE)
      Jacobian(IE,IY,1) = ZERO
      Jacobian(IE,IZ,1) = ZERO
c
      Jacobian(IX,IE,1) = DPDZH
      Jacobian(IX,IX,1) = TWO*ZAVG(IX) + DPDZU(1)
      Jacobian(IX,IY,1) = DPDZU(2)
      Jacobian(IX,IZ,1) = DPDZU(3)
c
      Jacobian(IY,IE,1) = ZERO
      Jacobian(IY,IX,1) = ZAVG(IY)
      Jacobian(IY,IY,1) = ZAVG(IX)
      Jacobian(IY,IZ,1) = ZERO
c
      Jacobian(IZ,IE,1) = ZERO
      Jacobian(IZ,IX,1) = ZAVG(IZ)
      Jacobian(IZ,IY,1) = ZERO
      Jacobian(IZ,IZ,1) = ZAVG(IX)
c
c	MATRIX dG/dZ
c
      DO I = 1 , NSP
        DO J = 1 , NSP  
          Jacobian(I,J,2) = ZERO
          IF (I.EQ.J) THEN
            Jacobian(I,J,2) = Jacobian(I,J,2) + ZAVG(IY)
          ENDIF
        ENDDO
      ENDDO
c
      DO I = 1 , NSP
        Jacobian(I,IE,2) = ZERO
        Jacobian(I,IX,2) = ZERO
        Jacobian(I,IY,2) = ZAVG(I)
        Jacobian(I,IZ,2) = ZERO
      ENDDO
c
      DO J = 1 , NSP
        Jacobian(IE,J,2) = ZERO
        Jacobian(IX,J,2) = ZERO
        Jacobian(IY,J,2) = DPDZS(J)
        Jacobian(IZ,J,2) = ZERO
      ENDDO
c
      Jacobian(IE,IE,2) = ZAVG(IY)
      Jacobian(IE,IX,2) = ZERO
      Jacobian(IE,IY,2) = ZAVG(IE)
      Jacobian(IE,IZ,2) = ZERO
c
      Jacobian(IX,IE,2) = ZERO
      Jacobian(IX,IX,2) = ZAVG(IY)
      Jacobian(IX,IY,2) = ZAVG(IX)
      Jacobian(IX,IZ,2) = ZERO
c
      Jacobian(IY,IE,2) = DPDZH
      Jacobian(IY,IX,2) = DPDZU(1)
      Jacobian(IY,IY,2) = TWO*ZAVG(IY) + DPDZU(2)
      Jacobian(IY,IZ,2) = DPDZU(3)
c
      Jacobian(IZ,IE,2) = ZERO
      Jacobian(IZ,IX,2) = ZERO
      Jacobian(IZ,IY,2) = ZAVG(IZ)
      Jacobian(IZ,IZ,2) = ZAVG(IY)
c
c	Matrix dH/dZ
c
      DO I = 1 , NSP
        DO J = 1 , NSP  
          Jacobian(I,J,3) = ZERO
          IF (I.EQ.J) THEN
            Jacobian(I,J,3) = Jacobian(I,J,3) + ZAVG(IZ)
          ENDIF
        ENDDO
      ENDDO
c
      DO I = 1 , NSP
        Jacobian(I,IE,3) = ZERO
        Jacobian(I,IX,3) = ZERO
        Jacobian(I,IY,3) = ZERO
        Jacobian(I,IZ,3) = ZAVG(I)
      ENDDO
c
      DO J = 1 , NSP
        Jacobian(IE,J,3) = ZERO
        Jacobian(IX,J,3) = ZERO     
        Jacobian(IY,J,3) = ZERO
        Jacobian(IZ,J,3) = DPDZS(J)
      ENDDO
c
      Jacobian(IE,IE,3) = ZAVG(IZ)
      Jacobian(IE,IX,3) = ZERO
      Jacobian(IE,IY,3) = ZERO
      Jacobian(IE,IZ,3) = ZAVG(IE)
c
      Jacobian(IX,IE,3) = ZERO
      Jacobian(IX,IX,3) = ZAVG(IZ)
      Jacobian(IX,IY,3) = ZERO
      Jacobian(IX,IZ,3) = ZAVG(IX)
c
      Jacobian(IY,IE,3) = ZERO
      Jacobian(IY,IX,3) = ZERO
      Jacobian(IY,IY,3) = ZAVG(IZ)
      Jacobian(IY,IZ,3) = ZAVG(IY)
c
      Jacobian(IZ,IE,3) = DPDZH
      Jacobian(IZ,IX,3) = DPDZU(1)
      Jacobian(IZ,IY,3) = DPDZU(2)
      Jacobian(IZ,IZ,3) = TWO*ZAVG(IZ) + DPDZU(3)
cC
c
c	Computes the conserved flux
c
  300 IDIM = 1
      CALL DGEMV( 'N', NVAR, NVAR, ONE, Jacobian(1,1,IDIM), LDA,
     .           GRAD_PARM(1,IDIM), 1, ZERO, DivFlux, 1 )
c
      DO 1 IDIM = 2 , NDIM
      CALL DGEMV( 'N', NVAR, NVAR, ONE, Jacobian(1,1,IDIM), LDA,
     .           GRAD_PARM(1,IDIM), 1, ONE, DivFlux, 1 )
   1  CONTINUE
C
C     .. The following call computes the flux integeral over
C        the cell by integrating the flux function over
C        the faces of the cell; it is intended to be used
C        to check that this coincides with the quasi linear
C        form computed above and stored in DivFlux ..
C
caldo CALL FLXBNC( -IELEM )
c
      RETURN
      ENTRY CHECK2(IELEM,NDIM,NVAR)
c
c     JACOBIAN IN Parameter Vector (Compressible)
c
c	Matrix dF/dZ
c
      Jacobian(1,1,1) = ZAVG(3)
      Jacobian(1,2,1) = ZERO
      Jacobian(1,3,1) = ZAVG(1)
      Jacobian(1,4,1) = ZERO
c
      Jacobian(2,1,1) = ZERO
      Jacobian(2,2,1) = ZAVG(3)
      Jacobian(2,3,1) = ZAVG(2)
      Jacobian(2,4,1) = ZERO
c
      Jacobian(3,1,1) = GM1OG * ZAVG(2)
      Jacobian(3,2,1) = GM1OG * ZAVG(1)
      Jacobian(3,3,1) = GP1OG * ZAVG(3)
      Jacobian(3,4,1) =-GM1OG * ZAVG(4)
c
      Jacobian(4,1,1) = ZERO
      Jacobian(4,2,1) = ZERO
      Jacobian(4,3,1) = ZAVG(4)
      Jacobian(4,4,1) = ZAVG(3)
c
      IF(NDIM.EQ.3)THEN
         Jacobian(1,5,1) = ZERO
         Jacobian(2,5,1) = ZERO
         Jacobian(3,5,1) =-GM1OG * ZAVG(5)
         Jacobian(4,5,1) = ZERO
c
         Jacobian(5,1,1) = ZERO
         Jacobian(5,2,1) = ZERO
         Jacobian(5,3,1) = ZAVG(5)
         Jacobian(5,4,1) = ZERO
         Jacobian(5,5,1) = ZAVG(3)
      ENDIF
c
      Jacobian(1,NVAR,1) = ZERO
      Jacobian(2,NVAR,1) = ZERO
      Jacobian(3,NVAR,1) = ZERO
      Jacobian(4,NVAR,1) = ZERO
      IF(NDIM.EQ.3)Jacobian(5,NVAR,1) = ZERO
c
      Jacobian(NVAR,1,1) = ZERO
      Jacobian(NVAR,2,1) = ZERO
      Jacobian(NVAR,3,1) = ZAVG(NVAR)
      Jacobian(NVAR,4,1) = ZERO
      IF(NDIM.EQ.3)Jacobian(NVAR,5,1) = ZERO
      Jacobian(NVAR,NVAR,1) = ZAVG(3)
c
c	MATRIX dG/dZ
c
      Jacobian(1,1,2) = ZAVG(4)
      Jacobian(1,2,2) = ZERO
      Jacobian(1,3,2) = ZERO
      Jacobian(1,4,2) = ZAVG(1)
c
      Jacobian(2,1,2) = ZERO
      Jacobian(2,2,2) = ZAVG(4)
      Jacobian(2,3,2) = ZERO
      Jacobian(2,4,2) = ZAVG(2)
c
      Jacobian(3,1,2) = ZERO
      Jacobian(3,2,2) = ZERO
      Jacobian(3,3,2) = ZAVG(4)
      Jacobian(3,4,2) = ZAVG(3)
c
      Jacobian(4,1,2) = GM1OG * ZAVG(2)
      Jacobian(4,2,2) = GM1OG * ZAVG(1)
      Jacobian(4,3,2) =-GM1OG * ZAVG(3)
      Jacobian(4,4,2) = GP1OG * ZAVG(4)
c
      IF(NDIM.EQ.3)THEN
         Jacobian(1,5,2) = ZERO
         Jacobian(2,5,2) = ZERO
         Jacobian(3,5,2) = ZERO
         Jacobian(4,5,2) =-GM1OG * ZAVG(5)
         Jacobian(5,1,2) = ZERO
         Jacobian(5,2,2) = ZERO
         Jacobian(5,3,2) = ZERO
         Jacobian(5,4,2) = ZAVG(5)
         Jacobian(5,5,2) = ZAVG(4)
      ENDIF
c
      Jacobian(1,NVAR,2) = ZERO
      Jacobian(2,NVAR,2) = ZERO
      Jacobian(3,NVAR,2) = ZERO
      Jacobian(4,NVAR,2) = ZERO
      IF(NDIM.EQ.3)Jacobian(5,NVAR,2) = ZERO
c
      Jacobian(NVAR,1,2) = ZERO
      Jacobian(NVAR,2,2) = ZERO
      Jacobian(NVAR,3,2) = ZERO
      Jacobian(NVAR,4,2) = ZAVG(NVAR)
      IF(NDIM.EQ.3)Jacobian(NVAR,5,2) = ZERO
      Jacobian(NVAR,NVAR,2) = ZAVG(4)
c
c	Matrix dH/dZ
c
      IF(NDIM.EQ.3)THEN
         Jacobian(1,1,3) = ZAVG(5)
         Jacobian(1,2,3) = ZERO
         Jacobian(1,3,3) = ZERO
         Jacobian(1,4,3) = ZERO
         Jacobian(1,5,3) = ZAVG(1)
c
         Jacobian(2,1,3) = ZERO
         Jacobian(2,2,3) = ZAVG(5)
         Jacobian(2,3,3) = ZERO
         Jacobian(2,4,3) = ZERO
         Jacobian(2,5,3) = ZAVG(2)
c
         Jacobian(3,1,3) = ZERO
         Jacobian(3,2,3) = ZERO
         Jacobian(3,3,3) = ZAVG(5)
         Jacobian(3,4,3) = ZERO
         Jacobian(3,5,3) = ZAVG(3)
c
         Jacobian(4,1,3) = ZERO
         Jacobian(4,2,3) = ZERO
         Jacobian(4,3,3) = ZERO
         Jacobian(4,4,3) = ZAVG(5)
         Jacobian(4,5,3) = ZAVG(4)
c
         Jacobian(5,1,3) = GM1OG * ZAVG(2)
         Jacobian(5,2,3) = GM1OG * ZAVG(1)
         Jacobian(5,3,3) =-GM1OG * ZAVG(3)
         Jacobian(5,4,3) =-GM1OG * ZAVG(4)
         Jacobian(5,5,3) = GP1OG * ZAVG(5)
c
         Jacobian(1,NVAR,3) = ZERO
         Jacobian(2,NVAR,3) = ZERO
         Jacobian(3,NVAR,3) = ZERO
         Jacobian(4,NVAR,3) = ZERO
         Jacobian(5,NVAR,3) = ZERO
c
         Jacobian(NVAR,1,3) = ZERO
         Jacobian(NVAR,2,3) = ZERO
         Jacobian(NVAR,3,3) = ZERO
         Jacobian(NVAR,4,3) = ZERO
         Jacobian(NVAR,5,3) = ZAVG(6)
         Jacobian(NVAR,NVAR,3) = ZAVG(5)
      ENDIF
cxxx  CALL R8Mat_Print('G',' ',NVAR,NDIM,GRAD_PARM,NMAX,
cxxx +      'grad(Z) ',IFAIL)
c
c	Computes the conserved flux
c
      IDIM = 1
      CALL DGEMV( 'N', NVAR, NVAR, ONE, Jacobian(1,1,IDIM), LDA,
     .           GRAD_PARM(1,IDIM), 1, ZERO, DivFlux, 1 )
cxxx  CALL R8Mat_Print('G',' ',NVAR,NVAR,Jacobian(1,1,IDIM),LDA,
cxxx +      'dFdZ(1) ',IFAIL)
c
      DO 7 IDIM = 2 , NDIM
cxxx  CALL R8Mat_Print('G',' ',NVAR,NVAR,Jacobian(1,1,IDIM),LDA,
cxxx +      'dFdZ(2,3) ',IFAIL)
      CALL DGEMV( 'N', NVAR, NVAR, ONE, Jacobian(1,1,IDIM), LDA,
     .           GRAD_PARM(1,IDIM), 1, ONE, DivFlux, 1 )
   7  CONTINUE
C
      RETURN
      END
C
