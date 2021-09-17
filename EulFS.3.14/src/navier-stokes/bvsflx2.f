      SUBROUTINE BVSFLX2(ICLR,IELEM,KVERT,ZROE,NodRes,TSTEP,NOFVAR,
     +                   NOFVERT,VCN,NDIM,VOLUME,STIFC,VISCL,VISCT,
     +                   TAUW,QFLUX,WALL,MATRIX_ASSEMBLY)
C
C     $Id: bvsflx2.f,v 1.15 2020/03/28 09:52:52 abonfi Exp $
C
      IMPLICIT NONE 
C
C     This routine computes the viscous fluxes which appear
C     in the INCOMPRESSIBLE Navier Stokes eqns. on boundary cells
C
C
C
      INCLUDE 'paramt.h'
      INCLUDE 'bnd.h'
      INCLUDE 'constants.h'
C
C
      INCLUDE 'dofs.com'
      INCLUDE 'visco.com'
      INCLUDE 'bnd.com'
      INCLUDE 'three.com'
C
C
C
C     output:
C     ------
C     NodRes(1:NOFVAR,1:NOFVERT) is the nodal residual
C                             updated with the addition of the viscous
C                             fluxes
C     TSTEP(1:NPOIN) is the nodal timestep divided by the median
C                 dual cell updated with the addition of the viscous
C                 contribution
C
C     STIFC(1:NOFVAR,1:NOFVAR,1:NOFVERT,1:NOFVERT) is the approximate
C            jacobian updated with the addition of the viscous terms
C
C
C
C
C
C     .. Parameters ..
      INTEGER NMAX2
      PARAMETER (NMAX2=NMAX*NMAX)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION QFLUX,TAUW,VISCL,VISCT,VOLUME
      INTEGER ICLR,IELEM,KVERT,NDIM,NOFVAR,NOFVERT
      LOGICAL MATRIX_ASSEMBLY,WALL
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION TSTEP(NOFVAR,NOFVERT),NodRes(NOFVAR,NOFVERT),
     +                 STIFC(NOFVAR,NOFVAR,NOFVERT,NOFVERT),
     +                 VCN(NDIM,NOFVERT),ZROE(NOFVAR,NOFVERT)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION CNST,LI,LJ,MI,MJ,MU,NI,NI_DOT_NJ,NJ,SURFSQR,DTV
      INTEGER I,IFAIL,IV,IVAR,IVERT,JV,K,NONAME
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION GFLUX(NMAX,VMAX),SIGMA(3),TAU(3,3),
     +                 VSFLX(NMAX,VMAX)
C     ..
C     .. External Functions ..
      DOUBLE PRECISION DDOT,DNRM2
      INTEGER ICYCL
      EXTERNAL DDOT,DNRM2
C     ..
C     .. External Subroutines ..
      EXTERNAL DGEMV,DINIT
C     ..
C     .. Data statements ..
C
C
      DATA TAU,NI,NJ/9*ZERO,2*ZERO/
C     ..
C
      NONAME = NDIM + 1
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
C     Compute viscous stress only if on a solid surface
C
      IF (WALL) THEN
C
C     Compute the stress vector
C
          SIGMA(1) = DDOT(NDIM,TAU(1,1),1,VCN(1,KVERT),1)
          SIGMA(2) = DDOT(NDIM,TAU(1,2),1,VCN(1,KVERT),1)
          IF (NDIM.EQ.3) SIGMA(3) = DDOT(NDIM,TAU(1,3),1,VCN(1,KVERT),1)
C
          SURFSQR = ONE/DDOT(NDIM,VCN(1,KVERT),1,VCN(1,KVERT),1)
C
          IF (NDIM.EQ.2) THEN
              TAUW = - (-VCN(2,KVERT)*SIGMA(1)+VCN(1,KVERT)*SIGMA(2))*
     +               SURFSQR

          ELSE
              CNST = DDOT(NDIM,SIGMA(1),1,VCN(1,KVERT),1)*SURFSQR
              LI = SIGMA(1) - CNST*VCN(1,KVERT)
              MI = SIGMA(2) - CNST*VCN(2,KVERT)
              NI = SIGMA(3) - CNST*VCN(3,KVERT)
              TAUW = SQRT(LI*LI+MI*MI+NI*NI)
          ENDIF
C
          VISCF(1,ICLR) = VISCF(1,ICLR) + SIGMA(1)
          VISCF(2,ICLR) = VISCF(2,ICLR) + SIGMA(2)
          IF (NDIM.EQ.3) VISCF(3,ICLR) = VISCF(3,ICLR) + SIGMA(3)
C
          QFLUX = 0.d0
C
      ENDIF
C
C     CNST = 1.d0/(d * Re)
C
      CNST = -REINV/NDIM
C
C     Compute the viscous fluxes for each of the (NOFVERT-1) vertices
C     on the boundary
C
C     the viscous term looks like:
C     -\frac{1}{d} \n_{boundary} \cdot [ 0, \tau ]
C
      DO 17 I = 1, (NOFVERT-1)
C
          IVERT = ICYCL(KVERT+I,NOFVERT)
          DTV = -CNST*DDOT(NDIM,VCN(1,IVERT),1,VCN(1,IVERT),1)/
     &(NDIM*VOLUME)
C
          TSTEP(1,IVERT) = TSTEP(1,IVERT) + DTV
          TSTEP(IX,IVERT) = TSTEP(IX,IVERT) + DTV
          TSTEP(IY,IVERT) = TSTEP(IY,IVERT) + DTV
          IF(NDIM.EQ.3)TSTEP(IZ,IVERT) = TSTEP(IZ,IVERT) + DTV
C
          GFLUX(1,IVERT) = ZERO
          GFLUX(2,IVERT) = -DDOT(NDIM,VCN(1,KVERT),1,TAU(1,1),1)/NDIM
          GFLUX(3,IVERT) = -DDOT(NDIM,VCN(1,KVERT),1,TAU(1,2),1)/NDIM
          IF (NDIM.EQ.3) GFLUX(4,IVERT) = -DDOT(NDIM,VCN(1,KVERT),1,
     +        TAU(1,3),1)/NDIM
C
          NodRes(2,IVERT) = NodRes(2,IVERT) + GFLUX(2,IVERT)
          NodRes(3,IVERT) = NodRes(3,IVERT) + GFLUX(3,IVERT)
          IF (NDIM.EQ.3) NodRes(4,IVERT) = NodRes(4,IVERT) +
     & GFLUX(4,IVERT)
C
   17 CONTINUE
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
      DO 10 I = 1, (NOFVERT-1)
C
          IV = ICYCL(KVERT+I,NOFVERT)
C
          LI = VCN(1,KVERT)
          MI = VCN(2,KVERT)
          IF (NDIM.EQ.3) NI = VCN(3,KVERT)
C
          DO 10 JV = 1,NOFVERT
C
              LJ = VCN(1,JV)
              MJ = VCN(2,JV)
              IF (NDIM.EQ.3) NJ = VCN(3,JV)
C
              NI_DOT_NJ = LI*LJ + MI*MJ + NI*NJ
C
C     ... Momentum flux
C
        STIFC(2,2,IV,JV) = MU* (LI*LJ+NI_DOT_NJ)
        STIFC(2,3,IV,JV) = MU*MI*LJ
        STIFC(3,2,IV,JV) = MU*LI*MJ
        STIFC(3,3,IV,JV) = MU* (MI*MJ+NI_DOT_NJ)
C
              IF (NDIM.EQ.3) THEN
C
           STIFC(2,4,IV,JV) = MU*NI*LJ
           STIFC(3,4,IV,JV) = MU*NI*MJ
C
           STIFC(4,2,IV,JV) = MU*LI*NJ
           STIFC(4,3,IV,JV) = MU*MI*NJ
           STIFC(4,4,IV,JV) = MU* (NI*NJ+NI_DOT_NJ)
              ENDIF
C
   10 CONTINUE
C
C comment out the following RETURN if you wish to debug the
C viscous flux calculation
C
      RETURN
C
C     Debugging stuff ...
C
      CALL DINIT(NMAX*VMAX,ZERO,VSFLX,1)
C
C     The "explicitely" calculated viscous flux is compared
C        with the "implicit" one ...
C
      DO 12 K = 1, (NOFVERT-1)
          IV = ICYCL(KVERT+K,NOFVERT)
          DO 14 JV = 1,NOFVERT
C
              CALL DGEMV('N',NONAME,NONAME,ONE,STIFC(1,1,IV,JV),NOFVAR,
     +                   ZROE(1,JV),1,ONE,VSFLX(1,IV),1)

C        CALL R8Mat_Print('General',' ',NOFVAR,NOFVAR,STIFC(1,1,IV,JV),
C    +   NOFVAR,'Diffusion element matrix',IFAIL)
C
   14     CONTINUE

          WRITE (6,FMT=*) 'Element # ',IELEM,
     +      ' (on the boundary) vertex # ',IV
C
          DO 15 IVAR = 2,NONAME
              WRITE (6,FMT=100) IVAR,GFLUX(IVAR,IV),VSFLX(IVAR,IV),
     +          VSFLX(IVAR,IV)/GFLUX(IVAR,IV)
   15     CONTINUE
C
   12 CONTINUE
C
C     PAUSE
C
      RETURN

  100 FORMAT (I1,3 (5X,E18.8))

      END
