      SUBROUTINE VSFLX2(IELEM,ZROE,NodRes,TSTEP,NOFVAR,VCN,NDIM,NOFVERT,
     +                  VOLUME,STIFD,VISCL,VISCT,MATRIX_ASSEMBLY)
C
C     $Id: vsflx2.f,v 1.11 2020/03/28 09:52:52 abonfi Exp $
C
      IMPLICIT NONE 
C
C     Purpose:
C     -------
C
C     This routine computes the viscous fluxes
C     for the INCOMPRESSIBLE Navier Stokes eqns. 
C     for INTERIOR elements
C
C     input:
C     -----
C     IELEM         is the current element
C     NDIM          is the dimension of the space (=2 or 3)
C     NOFVERT       is the no of vertices of the cell (=NDIM+1)
C     NOFVAR        is the no of degrees of freedom in the vertices
C
C     ZROE(1:NOFVAR,1:NOFVERT) is the parameter vector
C     VCN(1:NDIM,1:NOFVERT) keeps the NDIM cartesian components of
C                           the NOFVERT face normals
C     VISCL is the non-dimensional laminar viscosity
C     VISCT is the non-dimensional turbulent viscosity
C
C
C     output:
C     ------
C     NodRes(1:NOFVAR,1:NOFVERT) is the nodal residual
C                             updated with the addition of the viscous
C                             fluxes
C     TSTEP(NOFVAR,1:NPOIN) is the nodal timestep divided by the median
C                 dual cell updated with the addition of the viscous
C                 contribution
C
C     STIFD(1:NOFVAR,1:NOFVAR,1:NOFVERT,1:NOFVERT) is the approximate
C            jacobian updated with the addition of the viscous terms
C
C
      INCLUDE 'paramt.h'
      INCLUDE 'constants.h'
C
C
      INCLUDE 'dofs.com'
      INCLUDE 'visco.com'
      INCLUDE 'three.com'
C
C
C
C
C
C
C
C     .. Parameters ..
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION VISCL,VISCT,VOLUME
      INTEGER IELEM,NDIM,NOFVAR,NOFVERT
      LOGICAL MATRIX_ASSEMBLY
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION TSTEP(NOFVAR,NOFVERT),NodRes(NOFVAR,NOFVERT),
     +                 STIFD(NOFVAR,NOFVAR,NOFVERT,NOFVERT),
     +                 VCN(NDIM,NOFVERT),ZROE(NOFVAR,NOFVERT)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION CNST,LI,LJ,MI,MJ,MU,NI,NI_DOT_NJ,NJ,DTV
      DOUBLE PRECISION A22,A23,A32,A33,A24,A42,A43,A34,A44
      INTEGER I,IFAIL,IV,IVERT,J,JV,N4,NOFEQN
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION GFLUX(NMAX,VMAX),TAU(3,3),VSFLX(5,VMAX)
C     ..
C     .. External Functions ..
      DOUBLE PRECISION DDOT
      EXTERNAL DDOT
C     ..
C     .. External Subroutines ..
      EXTERNAL DGEMV,DINIT,R8Mat_Print
C     ..
C     .. Data statements ..
      DATA TAU,NI,NJ/9*ZERO,2*ZERO/
C     ..
C
      NOFEQN = NDIM + 1
C
C     Compute viscosity (VISCL must be 1)
C
      MU = (VISCL+VISCT)*REINV
C
C     ... Compute the stress tensor ...
C     \mu ( grad u + grad^T u )
C
      TAU(1,1) = MU*TWO* (GRAD_PARM(2,1))
      TAU(1,2) = MU* (GRAD_PARM(2,2)+GRAD_PARM(3,1))
      TAU(1,3) = MU* (GRAD_PARM(2,3)+GRAD_PARM(4,1))
      TAU(2,1) = TAU(1,2)
      TAU(2,2) = MU*TWO* (GRAD_PARM(3,2))
      TAU(2,3) = MU* (GRAD_PARM(3,3)+GRAD_PARM(4,2))
      TAU(3,1) = TAU(1,3)
      TAU(3,2) = TAU(2,3)
      TAU(3,3) = MU*TWO* (GRAD_PARM(4,3))
C
      CNST = -REINV/ (NDIM*NDIM*VOLUME)
C
C     Compute the viscous fluxes for each of the NOFVERT vertices
C
      DO 17 IVERT = 1,NOFVERT
C
          DTV = -CNST * DDOT(NDIM,VCN(1,IVERT),1,VCN(1,IVERT),1)
          TSTEP(1,IVERT) = TSTEP(1,IVERT) + DTV
          TSTEP(IX,IVERT) = TSTEP(IX,IVERT) + DTV
          TSTEP(IY,IVERT) = TSTEP(IY,IVERT) + DTV
          IF(NDIM.EQ.3)TSTEP(IZ,IVERT) = TSTEP(IZ,IVERT) + DTV
C
C     Viscous flux ..
C
          GFLUX(1,IVERT) = ZERO
          GFLUX(2,IVERT) = -DDOT(NDIM,VCN(1,IVERT),1,TAU(1,1),1)/
     &                     NDIM
          GFLUX(3,IVERT) = -DDOT(NDIM,VCN(1,IVERT),1,TAU(1,2),1)/
     &                     NDIM
          IF (NDIM.EQ.3) GFLUX(4,IVERT) = -DDOT(NDIM,VCN(1,IVERT),1,
     +        TAU(1,3),1)/NDIM
C
          NodRes(2,IVERT) = NodRes(2,IVERT) + GFLUX(2,IVERT)
          NodRes(3,IVERT) = NodRes(3,IVERT) + GFLUX(3,IVERT)
          IF (NDIM.EQ.3) NodRes(4,IVERT) = NodRes(4,IVERT) + 
     &                   GFLUX(4,IVERT)
C
   17 CONTINUE
C
C    write(6,*)'ielem =',ielem
C     CALL R8Mat_Print('General',' ',NOFVAR,NOFVERT,GFLUX,
C    +NMAX,'Diffusion element matrix',IFAIL)
C
C     Compute the stress vectors
C
      IF (.NOT.MATRIX_ASSEMBLY) RETURN
C
C     contruct a matrix D_{i,j} such that the viscous terms
C     can be written as \sum_j D_{i,j} U_j
C
      MU = MU/ (NDIM*NDIM*VOLUME)
C
C     ... Loop over vertices ...
C
      DO 10 JV = 1,NOFVERT
C
          LJ = VCN(1,JV)
          MJ = VCN(2,JV)
          IF (NDIM.EQ.3) NJ = VCN(3,JV)
C
          DO 10 IV = 1,JV
C
              LI = VCN(1,IV)
              MI = VCN(2,IV)
              IF (NDIM.EQ.3) NI = VCN(3,IV)
C
              NI_DOT_NJ = LI*LJ + MI*MJ + NI*NJ
C
C     ... Momentum flux
C
              A22 = MU* (LI*LJ+NI_DOT_NJ)
              A33 = MU* (MI*MJ+NI_DOT_NJ)
              A23 = MU*MI*LJ
              A32 = MU*LI*MJ
              STIFD(2,2,IV,JV) = STIFD(2,2,IV,JV) + A22
              STIFD(3,3,IV,JV) = STIFD(3,3,IV,JV) + A33
              STIFD(2,3,IV,JV) = STIFD(2,3,IV,JV) + A23 
              STIFD(3,2,IV,JV) = STIFD(3,2,IV,JV) + A32
C
              IF (NDIM.EQ.3) THEN
                  A24 = MU*NI*LJ
                  A34 = MU*NI*MJ
                  A42 = MU*LI*NJ
                  A43 = MU*MI*NJ
                  A44 = MU* (NI*NJ+NI_DOT_NJ)
                  STIFD(2,4,IV,JV) = STIFD(2,4,IV,JV) + A24
                  STIFD(3,4,IV,JV) = STIFD(3,4,IV,JV) + A34
                  STIFD(4,2,IV,JV) = STIFD(4,2,IV,JV) + A42
                  STIFD(4,3,IV,JV) = STIFD(4,3,IV,JV) + A43
                  STIFD(4,4,IV,JV) = STIFD(4,4,IV,JV) + A44
              ENDIF
C
              IF (IV.EQ.JV) GOTO 10
C
C        ... D_{ji} = D_{ij}^T ...
C
              STIFD(2,2,JV,IV) = STIFD(2,2,JV,IV) + A22
              STIFD(3,3,JV,IV) = STIFD(3,3,JV,IV) + A33
              STIFD(2,3,JV,IV) = STIFD(2,3,JV,IV) + A32 
              STIFD(3,2,JV,IV) = STIFD(3,2,JV,IV) + A23
              IF (NDIM.EQ.3) THEN
                  STIFD(2,4,JV,IV) = STIFD(2,4,JV,IV) + A42
                  STIFD(3,4,JV,IV) = STIFD(3,4,JV,IV) + A43
                  STIFD(4,2,JV,IV) = STIFD(4,2,JV,IV) + A24
                  STIFD(4,3,JV,IV) = STIFD(4,3,JV,IV) + A34
                  STIFD(4,4,JV,IV) = STIFD(4,4,JV,IV) + A44
              ENDIF
C
caldo         DO 8 J = 2,NOFVAR
caldo             DO 8 I = 2,NOFVAR
caldo                 STIFD(I,J,JV,IV) = STIFD(J,I,IV,JV)
cald8         CONTINUE
C
   10 CONTINUE
C
C     comment the following RETURN if you wish to debug the
C     viscous flux calculation
C
      RETURN
C
C     ... Debugging stuff ...
C     the following will work ONLY if STIFD is 0.d0 upon entry
C
      CALL DINIT(5*VMAX,ZERO,VSFLX,1)
C
C     .. The "explicitly" calculated viscous flux is compared
C        with the "implicit" one ...
C
      DO 12 IV = 1,NOFVERT
          DO 14 JV = 1,NOFVERT
C
          CALL DGEMV('N',NOFEQN,NOFEQN,ONE,STIFD(1,1,IV,JV),NOFVAR,
     +                   ZROE(1,JV),1,ONE,VSFLX(1,IV),1)
C
          CALL R8Mat_Print('General',' ',NOFVAR,NOFVAR,STIFD(1,1,IV,JV),
     +                    NOFVAR,'Diffusion element matrix',IFAIL)
C
   14     CONTINUE
C
          WRITE (6,FMT=*) 'Element # ',IELEM,' vertex # ',IV
C
          DO 15 I = 2,NOFEQN
              WRITE (6,FMT=100) I,GFLUX(I,IV),VSFLX(I,IV),
     +          VSFLX(I,IV)/MAX(1.D-20,GFLUX(I,IV))
   15     CONTINUE
C
   12 CONTINUE
C
      PAUSE
C
      RETURN

  100 FORMAT (I1,3 (5X,E18.8))

      END
