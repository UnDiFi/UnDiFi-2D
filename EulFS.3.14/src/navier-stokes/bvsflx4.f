!> \par Purpose
!>
!> Compute viscous fluxes on boundary cells (compressible case)
!>
!>
!> Whenever the node lies on the boundary, additional terms due to the
!> boundary fluxes have to be accounted for; triangle \c IELEM contributes with
!> a term:
!> \f[
!> - \frac{1}{d} \mathbf{n}_{k}^{ie} \cdot \mathbf{G}
!> \f]
!> where \c k is the vertex facing the boundary face and \c ie = \c IELEM
!>
!>
!> @param[in] ICLR colour of the current boundary patch
!> @param[in] IELEM counter of the current triangle/tetrahedron
!> @param[in] KVERT index of the vertex facing the boundary face
!> @param[in] ZROE is parameter vector in the NOFVERT vertices of cell \c IELEM
!> @param[in,out] NODRES nodal residual
!> @param[in,out] TSTEP nodal timestep
!> @param[in] NOFVAR is the nof dofs
!> @param[in] NOFVERT nof vertices (= NDIM+1)
!> @param[in] VCN Cartesian components of the normals to a face, multiplied by the face area
!> @param[in] NDIM dimension of the space
!> @param[in] VOLUME area/volume of the current element (triangle,tetrahedron)
!> @param[in,out] ELTMAT implicit matrix, built only when \c MATRIX_ASSEMBLY .EQV. .TRUE.
!> @param[in] VISCL laminar viscosity
!> @param[in] VISCT turbulent viscosity
!> @param[out] TAUW dimensionless shear stress
!> @param[out] HFLUX nof ghost nodes in the mesh
!> @param[in] WALL \c .EQV. \c .TRUE. if the current face is a boundary face
!> @param[in] MATRIX_ASSEMBLY \c .EQV. \c .TRUE. if the implicit matrix has to be built, i.e. Picard linearization
!>
!> \author $Author: abonfi $
!> \version $Revision: 1.22 $
!> \date $Date: 2020/02/27 10:46:41 $
!> \bug 3D wall shear stress is neither computed nor stored
!> \warning
!
      SUBROUTINE BVSFLX4(ICLR,IELEM,KVERT,ZROE,NodRes,TSTEP,NOFVAR,
     +                   NOFVERT,VCN,NDIM,VOLUME,ELTMAT,VISCL,VISCT,
     +                   TAUW,HFLUX,WALL,MATRIX_ASSEMBLY)
C
C     $Id: bvsflx4.f,v 1.22 2020/02/27 10:46:41 abonfi Exp abonfi $
C
      IMPLICIT NONE 
C
C     Compute viscous fluxes on boundary cells (compressible case)
C
C
C
      INCLUDE 'paramt.h'
      INCLUDE 'bnd.h'
      INCLUDE 'constants.h'
C
      INCLUDE 'bnd.com'
      INCLUDE 'dofs.com'
      INCLUDE 'pfcgas.com'
      INCLUDE 'three.com'
      INCLUDE 'transf.com'
      INCLUDE 'turb.com'
      INCLUDE 'visco.com'
C
C
C     .. Parameters ..
      DOUBLE PRECISION TWOTHIRD
      INTEGER NMAX2
      PARAMETER (TWOTHIRD=TWO/3.D0,NMAX2=NMAX*NMAX)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION HFLUX,TAUW,VISCL,VISCT,VOLUME,DUMMY
      INTEGER ICLR,IELEM,KVERT,NDIM,NOFVAR,NOFVERT
      LOGICAL MATRIX_ASSEMBLY,WALL
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION TSTEP(NOFVAR,NOFVERT),
     +       NodRes(NOFVAR,NOFVERT),ELTMAT(NOFVAR,NOFVAR,NOFVERT,
     +       NOFVERT),VCN(NDIM,NOFVERT),ZROE(NOFVAR,NOFVERT)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION CNST,DIVV,LD,LI,LJ,MI,MJ,MPLUSL,MU,NI,NI_DOT_NJ,
     +                 NJ,SURFSQR,TEMPA,TEMPB,TEMPC,TEMPD,
     +                 U_DOT_NI,U_DOT_NJ,DSUM,DTV
      INTEGER I,IFAIL,IV,IVAR,IVERT,JV,K,IADDR,NOFEQN
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION GFLUX(NMAX,VMAX),SIGMA(3),TAU(3,4),
     +                 VSFLX(NMAX,VMAX),WKSP1(NMAX*VMAX),WORK(NMAX,NMAX)
C     ..
C     .. External Functions ..
      DOUBLE PRECISION DDOT,DNRM2,sutherlaw
      INTEGER ICYCL
      EXTERNAL DDOT,DNRM2,ICYCL
C     ..
C     .. External Subroutines ..
      EXTERNAL DCOPY,DGEMM,DGEMV,DINIT,PARM_TO_CONS
C     ..
C     .. Data statements ..
C
      DATA WORK,NI,NJ/NMAX2*ZERO,2*ZERO/
      DATA TAU/12*ZERO/
C     ..
C
C     call dinit((NOFVAR*NOFVERT)**2,ZERO,ELTMAT,1)
      NOFEQN = NDIM + 2
C
      MU = (VISCL+VISCT)*REINV
      LD = -TWOTHIRD*MU
C
C     if on a solid wall (which we assume is NOT moving)
C     set the cell averaged velocity to 0.d0
C
      IF(WALL)THEN
         UAVG(3) = ZERO
         UAVG(4) = ZERO
         IF(NDIM.EQ.3)UAVG(5) = ZERO
      ENDIF
C
      CALL GRADPRIM(IELEM,NDIM,NDIM+2)
C
C     Compute the divergence of the velocity field
C
      DIVV = GRAD_PRIM(3,1) + GRAD_PRIM(4,2)
      IF(NDIM.EQ.3)DIVV = DIVV + GRAD_PRIM(5,3)
C
C     Compute the stress tensor ...
C     \mu ( grad u + grad^T u ) + \lambda div u ]
C
      TAU(1,1) = TWO*MU*GRAD_PRIM(3,1) + LD*DIVV
      TAU(2,2) = TWO*MU*GRAD_PRIM(4,2) + LD*DIVV
      TAU(1,2) = MU* (GRAD_PRIM(3,2)+GRAD_PRIM(4,1))
      TAU(2,1) = TAU(1,2)
      IF(NDIM.EQ.3)THEN
          TAU(3,3) = TWO*MU*GRAD_PRIM(5,3) + LD*DIVV
          TAU(1,3) = MU* (GRAD_PRIM(3,3)+GRAD_PRIM(5,1))
          TAU(2,3) = MU* (GRAD_PRIM(4,3)+GRAD_PRIM(5,2))
          TAU(3,1) = TAU(1,3)
          TAU(3,2) = TAU(2,3)
      ENDIF
C
C     Compute the temperature gradient and energy term
C     (the latter is stored in the last column of TAU)
C
      CNST = (VISCL/PRANDTL+VISCT/TPR)/(GM1*UAVG(1))*REINV
C
C     \nabla Q = viscl / (gamma - 1.d0) / Pr * \nabla (a^2)
C
C     \nabla (a^2) = ( \gamma \nabla p - a^2 \nabla \rho )/ \rho
C
C
      TAU(1,4) = CNST* (GAM*GRAD_PRIM(2,1)-ASQR*GRAD_PRIM(1,1)) +
     +           UAVG(3)*TAU(1,1) + UAVG(4)*TAU(2,1)
      TAU(2,4) = CNST* (GAM*GRAD_PRIM(2,2)-ASQR*GRAD_PRIM(1,2)) +
     +           UAVG(3)*TAU(1,2) + UAVG(4)*TAU(2,2)
      IF(NDIM.EQ.3)THEN
          TAU(1,4) = TAU(1,4) + UAVG(5)*TAU(3,1)
          TAU(2,4) = TAU(2,4) + UAVG(5)*TAU(3,2)
          TAU(3,4) = CNST* (GAM*GRAD_PRIM(2,3)-ASQR*GRAD_PRIM(1,3)) +
     +           UAVG(3)*TAU(1,3) +
     2           UAVG(4)*TAU(2,3) +
     3           UAVG(5)*TAU(3,3)
      ENDIF
c
c     Here we enforce the term u \cdot stress tensor
c     as WELL as the temperature gradient to be zero
c     this applies only if adiabatic (IADIA=0) bcs
c     are applied
c
c
C     IF (WALL .AND. (IADIA.EQ.0) ) THEN
C        TAU(1,4) = 0.d0
C        TAU(2,4) = 0.d0
C        IF(NDIM.EQ.3) TAU(3,4) = 0.d0
C     ENDIF
C
C     Compute viscous stress only if on a solid surface
C
      IF (WALL) THEN
C
C     Compute the stress vector
C
          SIGMA(1) = VCN(1,KVERT)*TAU(1,1)+VCN(2,KVERT)*TAU(2,1)
          SIGMA(2) = VCN(1,KVERT)*TAU(1,2)+VCN(2,KVERT)*TAU(2,2)
          TEMPC = VCN(1,KVERT)*TAU(1,4)+VCN(2,KVERT)*TAU(2,4)
          SURFSQR = VCN(1,KVERT)*VCN(1,KVERT)+
     1              VCN(2,KVERT)*VCN(2,KVERT) 
          
          IF (NDIM.EQ.3) THEN
              SIGMA(1) = SIGMA(1) + VCN(3,KVERT)*TAU(3,1)
              SIGMA(2) = SIGMA(2) + VCN(3,KVERT)*TAU(3,2)
              TEMPC = TEMPC + VCN(3,KVERT)*TAU(3,4)
              SIGMA(3) = VCN(1,KVERT)*TAU(1,3)+
     2                   VCN(2,KVERT)*TAU(2,3)+
     3                   VCN(3,KVERT)*TAU(3,3)
              SURFSQR = SURFSQR + VCN(3,KVERT)*VCN(3,KVERT)
          ENDIF
C
          SURFSQR = ONE/SURFSQR
C
          HFLUX = TEMPC*DSQRT(SURFSQR)
C         write(6,*)ielem,hflux
C
          IF (NDIM.EQ.2) THEN
              TAUW = (-VCN(2,KVERT)*SIGMA(1)+VCN(1,KVERT)*SIGMA(2))*
     +               SURFSQR

          ELSE
              TAUW = 0.d0
          ENDIF
C
C
          VISCF(1,ICLR) = VISCF(1,ICLR) + SIGMA(1)
          VISCF(2,ICLR) = VISCF(2,ICLR) + SIGMA(2)
          IF (NDIM.EQ.3) VISCF(3,ICLR) = VISCF(3,ICLR) + SIGMA(3)
C
      ENDIF
C
      CNST = -REINV/NDIM
      DTV = REINV/(NDIM*NDIM*VOLUME)
C
C     Compute the viscous fluxes for each of the (NOFVERT-1) vertices
C     on the boundary
C
C     the viscous term looks like:
C     -\frac{1}{d} \n_{boundary} \cdot [ 0, \nabla q, \tau ]
C
      TEMPA = VCN(1,KVERT)*TAU(1,4)+VCN(2,KVERT)*TAU(2,4)
      TEMPB = VCN(1,KVERT)*TAU(1,1)+VCN(2,KVERT)*TAU(2,1)
      TEMPC = VCN(1,KVERT)*TAU(1,2)+VCN(2,KVERT)*TAU(2,2)
      IF(NDIM.EQ.3)THEN
             TEMPA = TEMPA + VCN(3,KVERT)*TAU(3,4)
             TEMPB = TEMPB + VCN(3,KVERT)*TAU(3,1)
             TEMPC = TEMPC + VCN(3,KVERT)*TAU(3,2)
             TEMPD = (
     1       VCN(1,KVERT)*TAU(1,3)+VCN(2,KVERT)*TAU(2,3)+
     2       VCN(3,KVERT)*TAU(3,3) )/NDIM
      ENDIF
C
C     enforce n \cdot \nabla T = 0 regardless of whether
C     it is a wall or not;
C
C
C     IF (WALL .AND. (IADIA.EQ.0) ) TEMPA = 0.d0
      IF ((IADIA.EQ.0) ) TEMPA = 0.d0
C
      DO 17 I = 1, (NOFVERT-1)
C
          IVERT = ICYCL(KVERT+I,NOFVERT)
C
          DSUM = VCN(1,IVERT)*VCN(1,IVERT)+
     1           VCN(2,IVERT)*VCN(2,IVERT)
          IF(NDIM.EQ.3)DSUM = DSUM + VCN(3,IVERT)*VCN(3,IVERT)
C                                        ^
C    should probably be KVERT but it does|not really matter 
C
          TSTEP(1,IVERT) = TSTEP(1,IVERT) + DTV*DSUM
          TSTEP(IE,IVERT) = TSTEP(IE,IVERT) + DTV*DSUM
          TSTEP(IX,IVERT) = TSTEP(IX,IVERT) + DTV*DSUM
          TSTEP(IY,IVERT) = TSTEP(IY,IVERT) + DTV*DSUM
          IF(NDIM.EQ.3)TSTEP(IZ,IVERT) = TSTEP(IZ,IVERT) + DTV*DSUM
C
C
          GFLUX(1,IVERT) = ZERO
          GFLUX(2,IVERT) = -TEMPA/NDIM
          GFLUX(3,IVERT) = -TEMPB/NDIM
          GFLUX(4,IVERT) = -TEMPC/NDIM
          IF(NDIM.EQ.3)   GFLUX(5,IVERT) = -TEMPD/NDIM
C
          NodRes(2,IVERT) = NodRes(2,IVERT) + GFLUX(2,IVERT)
          NodRes(3,IVERT) = NodRes(3,IVERT) + GFLUX(3,IVERT)
          NodRes(4,IVERT) = NodRes(4,IVERT) + GFLUX(4,IVERT)
          IF (NDIM.EQ.3) NodRes(5,IVERT) = NodRes(5,IVERT) + 
     &    GFLUX(5,IVERT)
C
   17 CONTINUE
C
C     All the implicit stuff begins here
C
      IF (.NOT.MATRIX_ASSEMBLY) RETURN
C
      DO 33 IV = 1, NOFVERT
         IADDR = (IV-1)*NOFEQN*NOFEQN+1
         CALL MatdZdU(ZROE(1,IV),DZDU(IADDR),NDIM,NOFEQN)
   33 CONTINUE
C
C
C     contruct a matrix D_{i,j} such that the viscous terms
C     can be written as \sum_j D_{i,j} U_j
C
      MPLUSL = MU + LD
      CNST = -ASQR/GM1 + KINETIC
      TEMPA = MU/PRANDTL
      TEMPA = (VISCL/PRANDTL+VISCT/TPR)*REINV
      TEMPC = ONE/ (NDIM*NDIM*VOLUME*ZAVG(1))
C
C     Loop over vertices ...
C
      DO 10 I = 1, (NOFVERT-1)
C
          IV = ICYCL(KVERT+I,NOFVERT)
C
          LI = VCN(1,KVERT)
          MI = VCN(2,KVERT)
          IF (NDIM.EQ.3) NI = VCN(3,KVERT)
C
          U_DOT_NI = UAVG(3)*LI + UAVG(4)*MI + UAVG(5)*NI
caldo     U_DOT_NI = ZERO
C
          DO 10 JV = 1,NOFVERT
C
              LJ = VCN(1,JV)
              MJ = VCN(2,JV)
              IF (NDIM.EQ.3) NJ = VCN(3,JV)
C
              U_DOT_NJ = UAVG(3)*LJ + UAVG(4)*MJ + UAVG(5)*NJ
caldo         U_DOT_NI = ZERO
C
              NI_DOT_NJ = LI*LJ + MI*MJ + NI*NJ
C
C     Be VERY careful with WORK: make sure that non zero elements
C     are cleared (overwritten) when going from one vertex to the next
C
C     Energy flux
C
              WORK(2,3) = MU* (U_DOT_NJ*LI+NI_DOT_NJ*UAVG(3)) +
     +                    LD*U_DOT_NI*LJ
              WORK(2,4) = MU* (U_DOT_NJ*MI+NI_DOT_NJ*UAVG(4)) +
     +                    LD*U_DOT_NI*MJ
              WORK(2,5) = MU* (U_DOT_NJ*NI+NI_DOT_NJ*UAVG(5)) +
     +                    LD*U_DOT_NI*NJ
C
caldo       WORK(2,3) = ZERO
caldo       WORK(2,4) = ZERO
caldo       WORK(2,5) = ZERO
C
C     Momentum flux
C
              WORK(3,3) = MPLUSL*LI*LJ + MU*NI_DOT_NJ
              WORK(4,4) = MPLUSL*MI*MJ + MU*NI_DOT_NJ
              WORK(5,5) = MPLUSL*NI*NJ + MU*NI_DOT_NJ
C
              WORK(3,4) = MU*MI*LJ + LD*LI*MJ
              WORK(3,5) = MU*NI*LJ + LD*LI*NJ
C
              WORK(4,3) = MU*LI*MJ + LD*MI*LJ
              WORK(4,5) = MU*NI*MJ + LD*MI*NJ
C
              WORK(5,3) = MU*LI*NJ + LD*NI*LJ
              WORK(5,4) = MU*MI*NJ + LD*NI*MJ
C
C     Transform THE DIFFUSION matrix from primitive
C     into parameter vector; note that the density^{-1/2} term is omitted
C
              WORK(2,1) = -UAVG(3)*WORK(2,3) - UAVG(4)*WORK(2,4) -
     +                     UAVG(5)*WORK(2,5)
              WORK(3,1) = -UAVG(3)*WORK(3,3) - UAVG(4)*WORK(3,4) -
     +                     UAVG(5)*WORK(3,5)
              WORK(4,1) = -UAVG(3)*WORK(4,3) - UAVG(4)*WORK(4,4) -
     +                     UAVG(5)*WORK(4,5)
              WORK(5,1) = -UAVG(3)*WORK(5,3) - UAVG(4)*WORK(5,4) -
     +                     UAVG(5)*WORK(5,5)
C
!             WORK(2,1) = ZERO
!             WORK(3,1) = ZERO
!             WORK(4,1) = ZERO
!             WORK(5,1) = ZERO
C
C
C        Add now the heat flux term (which is already in
C            parameter vector) ...
C            NOTE that again the density^{-1/2} term is missing
C
C
              TEMPB = NI_DOT_NJ*TEMPA
C
C
              WORK(2,1) = WORK(2,1) + TEMPB*CNST
              WORK(2,2) = TEMPB
              WORK(2,3) = WORK(2,3) - TEMPB*UAVG(3)
              WORK(2,4) = WORK(2,4) - TEMPB*UAVG(4)
              WORK(2,5) = WORK(2,5) - TEMPB*UAVG(5)
C
C        Transforms the diffusion matrix from parameter vector
C            into conserved variables ...
C
         IADDR = (JV-1)*NOFEQN*NOFEQN+1
C
C  N.B. the diffusion stiffness matrix is over-written
C       (BETA = 0.d0) in DGEMM
C
   20         CALL DGEMM('N','N',NOFEQN,NOFEQN,NOFEQN,TWO*TEMPC,WORK,
     +                   NMAX,DZDU(IADDR),NOFEQN,ZERO,ELTMAT(1,1,IV,JV),
     +                   NOFVAR)
C
C
   10 CONTINUE
C
      RETURN
C
C     Debugging stuff ...
C
      CALL DCOPY(NOFVAR*NOFVERT,ZROE,1,WKSP1,1)
C     Compute U_i i=1,..,NOFVERT
      CALL PARM_TO_CONS(WKSP1,NDIM,NOFVAR,NOFVERT,.FALSE.,IFAIL)
C
      CALL DINIT(NMAX*VMAX,ZERO,VSFLX,1)
C
C     The "explicitely" calculated viscous flux is compared
C     with the "implicit" one ...
C
      DO 12 K = 1, (NOFVERT-1)
          IV = ICYCL(KVERT+K,NOFVERT)
          DO 14 JV = 1,NOFVERT
C
              CALL DGEMV('N',NOFEQN,NOFEQN,ONE,ELTMAT(1,1,IV,JV),NOFVAR,
     +                   WKSP1((JV-1)*NOFVAR+1),1,ONE,VSFLX(1,IV),1)
C
   14     CONTINUE

          WRITE (6,FMT=*) 'Element # ',-IELEM,
     +      ' (on the boundary) vertex # ',IV,WALL
C
          DO 15 IVAR = 2,NOFEQN
              WRITE (6,FMT=100) IVAR,GFLUX(IVAR,IV),-VSFLX(IVAR,IV),
     +          (ONE+VSFLX(IVAR,IV)/GFLUX(IVAR,IV))*100.
   15     CONTINUE
C
   12 CONTINUE
C
C     PAUSE
C
      RETURN

  100 FORMAT (I1,3 (5X,E18.8))

      END
