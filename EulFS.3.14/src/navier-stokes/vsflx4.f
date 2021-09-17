!>
!> Compute viscous fluxes in the interior cells (compressible case)
!>
!> Calculation of the nodal residual
!> =================================
!>
!> Each vertex \c i of cell \c T receives a contribution due to the
!> viscous fluxes which amounts to:
!> \f[
!> {\tt GFLUX(:,IVERT)} =
!>   -\frac{1}{d} \mathbf{n}_i^T \cdot \mathbf{G}^T
!> \f]
!>
!> where:
!> \f[
!> \mathbf{G}^T = \left(
!>    \begin{array}{c}
!>  0 \\ \mathbf{u}  \cdot \mathbf{\tau} + \nabla q \\ \mathbf{\tau}
!>    \end{array}
!> \right)
!> \f]
!> is cell-wise constant (hence the superscript  <tt>T</tt>).
!>
!> In order to be consistent with the use of parameter vector as the set of dependent variables,
!> the gradient of the primitive variables \f$\left(\rho,p,u,v,w\right)\f$
!> is computed from the gradient of the parameter vector and stored in: 
!> \f[
!> {\tt GRAD\_PRIM} = \left(
!>    \begin{array}{ccc}
!>        \rho_x & \rho_x & \rho_x \\
!>        p_x & p_x & p_x  \\
!>        u_x & u_x & u_x  \\
!>        v_x & v_x & v_x  \\
!>        w_x & w_x & w_x 
!>    \end{array} \right)
!> \f]
!>
!> The stress tensor, which is also cell-wise constant, is stored in \f${\tt TAU}\f$
!> and it is computed using the entries of \f${\tt GRAD\_PRIM}\f$
!>
!> \f{eqnarray*}{ 
!> {\tt TAU} &=& \mu \left[ \nabla \mathbf{u} + \left(\nabla\mathbf{u}\right)^t \right] + \lambda \left(\nabla \cdot \mathbf{u} \right)\\ 
!> &=&
!> \left(
!>    \begin{array}{ccc}
!>    2\,\mu\,u_x +\lambda\,{\tt DIVV} & \mu\left(u_y+v_x\right) & \mu\left(u_z+w_x\right) \\
!>     \mu\left(v_x+u_y\right) & 2\,\mu\,v_y +\lambda\,{\tt DIVV} & \mu\left(v_z+w_y\right) \\
!>     \mu\left(w_x+u_z\right) & \mu\left(w_y+v_z\right) & 2\,\mu\,w_z +\lambda\,{\tt DIVV}
!>    \end{array}
!> \right)
!> \f}
!>
!> where \f${\tt DIVV} = u_x + v_y + w_z\f$, \f$\mu = {\tt VISCL} + {\tt VISCT}\f$ and \f$\lambda = -2/3\,\mu\f$.
!>
!>
!> The term that accounts for thermal conductivity is:
!>
!> \f[ 
!> \nabla q = \frac{1}{\gamma-1} \; k \; \nabla a^2
!> \f]
!>
!> where:
!>
!> \f[ 
!> k = \left( \frac{\mu}{\mbox{Pr}} + \frac{\mu_T}{\mbox{Pr}_T} \right) = {\tt VISCL/PRANDTL+VISCT/TPR}
!> \f]
!>
!> \f[ 
!> \nabla a^2 = \gamma \nabla \left( \frac{p}{\rho} \right) =
!> \frac{1}{\hat{\rho}} \left[\gamma\,\nabla p - \hat{a}^2\,\nabla\rho\right]
!> \f]
!>
!> Having set:
!>
!> \f[ 
!> {\tt CNST} = \frac{k}{\left(\gamma-1\right)\hat{\rho}}
!> \f]
!>
!> the viscous contribution to the energy equation is stored in the fourth column of the array \f${\tt TAU}\f$:
!>
!> \f[ 
!> {\tt TAU(:,4)} = {\tt CNST} \left[\gamma\,\nabla p - \hat{a}^2\,\nabla\rho\right]
!> + \hat{\mathbf{u}} \cdot \mathbf{\tau}
!> \f]
!>
!> The cell-averaged flow velocity \f$\hat{\mathbf{u}}\f$ is stored in \f${UAVG(3:5)}\f$;
!> the entries of \f${UAVG}\f$ are \f$\left(\hat{\rho},\hat{H},\hat{u},\hat{v},\hat{w}\right)\f$, where:
!>
!> \f[
!> \hat{\rho} = \widehat{Z}_1^2, \quad
!> \hat{H} = \widehat{Z}_2/\widehat{Z}_1, \quad
!> \hat{u} = \widehat{Z}_3/\widehat{Z}_1, \quad
!> \hat{v} = \widehat{Z}_4/\widehat{Z}_1, \quad
!> \hat{w} = \widehat{Z}_5/\widehat{Z}_1
!> \f]
!> and the cell-averaged value of the parameter vector is:
!> \f[ 
!> \hat{Z} = \frac{1}{d+1} \sum_{j=1}^{d+1} Z_j
!> \f]
!>
!> where the subscript \f$j\f$ refers to the vertex, not to the individual components of \f$Z\f$.
!>
!> Calculation of the time-step
!> ============================
!>
!> In order to mimick the scalar case, the time-step restriction should read:
!>
!> \f[
!> \frac{\Delta t}{V_i} \le \frac{1}{\sum d_{ij}} \quad\quad\mbox{where}\quad\quad d_{ij} = \frac{\nu}{d^2 \,V^T} \mathbf{n}_i \cdot \mathbf{n}_j
!> \f]
!>
!> \f$\nu = \mu/\hat{\rho}\f$ being the cell-averaged kinematic viscosity.
!>
!> Instead, the dimensionless version we use is:
!>
!> \f[
!> d_{ij} = \frac{1}{d^2 \, \mbox{Re} \, V^T} \mathbf{n}_i \cdot \mathbf{n}_j
!> \f]
!>
!> which amounts to assume that the dimensionless kinematic viscosity equals 1.
!>
!> Matrix assembly (only when \c MATRIX_ASSEMBLY \c .EQV. \c TRUE )
!> ============================
!>
!> the \c DMAT(:,:,i,j) matrix is such that
!>
!>
!> \f[
!> -\frac{1}{d} \mathbf{n}_i^T \cdot \mathbf{G}^T = \sum_{j=1}^{d+1} D_{ij} U_j
!> \f]
!>
!> where the summation ranges over all \c d+1 vertices of the cell and \f$U_j\f$ is the value of the conservative variable in vertex <tt>j</tt>.
!>
!> It can be shown that:
!>
!> \f[
!> D_{ij}^T = \frac{1}{d^2 \, V^T \, \mbox{Re}}
!> \left[ D_{ij}^{\tau} + D_{ij}^{\nabla q} \right] 2 \left( \frac{\partial U}{\partial Z} \right)_{Z_j}.
!> \f]
!> where:
!> \f[
!> D_{ij}^{\tau} =
!> \frac{1}{\widehat{\sqrt{\rho}}} \left(
!> \begin{array}{ccc}  0 & 0 & {\bf 0} \\
!> - {\bf d}_{ij} \cdot \hat{\mathbf{u}} & 0 & {\bf d}_{ij}^t \\
!> -\tilde{D}_{ij} \cdot \hat{\mathbf{u}} & 0 & \tilde{D}_{ij} \end{array} \right),
!> \f]
!> \f[
!> D_{ij}^{\nabla q} =
!> \frac{k}{\widehat{\sqrt{\rho}}}\,
!> \left( \begin{array}{ccr}
!> 0 & 0 & {\bf 0}^t \\
!> \hat{H} - \hat{a}^2/\delta &
!> 1 & -\hat{\mathbf{u}} \\
!> {\bf 0} & {\bf 0} & {\bf 0} \end{array}
!> \right) \left( \mathbf{n}_i \cdot \mathbf{n}_j \right)
!>\f]
!>\f[
!> {\bf d}_{ij} =
!> \mu \left( \hat{\mathbf{u}} \cdot \mathbf{n}_j \right) \mathbf{n}_i +
!>        \mu \left( \mathbf{n}_j \cdot \mathbf{n}_i \right) \hat{\mathbf{u}} +
!>        \lambda \left( \hat{\mathbf{u}} \cdot \mathbf{n}_i \right) \mathbf{n}_j
!>\f]
!>\f[
!>   \tilde{D}_{ij} =
!>   \left[ \mu \left( \mathbf{n}_j   \mathbf{n}_i \right) + \mu
!>   \left( \mathbf{n}_i \cdot \mathbf{n}_j \right) I +
!>   \lambda \left( \mathbf{n}_i   \mathbf{n}_j \right) \right]
!>\f]
!>observe that \f$\tilde{D}_{ij} = \tilde{D}^{t}_{ji}\f$ and that matrix \f$D_{ij}^{\nabla q}\f$, apart from the term \f$\mathbf{n}_i \cdot \mathbf{n}_j\f$, is cell-wise constant. 
!>
!> This is to say that we might save some calculations, which we are not doing right now.
!>
!>the vector \f${\bf d}_{ij}\f$ is stored in <tt>WORK(2,3:5)</tt>.
!>
!>the matrix \f$\tilde{D}_{ij}\f$ is stored in <tt>WORK(3:5,3:5)</tt>.
!>
!>the scalar \f$-{\bf d}_{ij}\cdot\hat{\mathbf{u}}\f$ is stored in <tt>WORK(2,1)</tt>.
!>
!>the vector \f$-\tilde{D}_{ij}\cdot\hat{\mathbf{u}}\f$ is stored in <tt>WORK(3:5,1)</tt>.
!>
!>so that, at this stage, apart from the term \f$\widehat{\sqrt{\rho}}\f$, matrix \f$D_{ij}^{\tau}\f$ has been built and stored in <tt>WORK</tt>.
!>
!>we set \f${\tt TEMPA} = \frac{\mu}{\mbox{Pr}} + \frac{\mu_T}{\mbox{Pr}_T}\f$
!>
!>we set \f${\tt TEMPB} = {\tt TEMPA} \left( \mathbf{n}_i \cdot \mathbf{n}_j\right)\f$
!>
!>we set:
!>\f[
!> {\tt TEMPC} = \frac{1}{d^2\,\widehat{\sqrt{\rho}}\,\mbox{Re}\,V^T}
!>\f]
!>
!> then add matrix, apart from the term \f$\widehat{\sqrt{\rho}}\f$, \f$D_{ij}^{\nabla q}\f$ to <tt>WORK</tt>.
!>
!> Finally, we compute matrix \f$D_{ij}\f$ by calling the Blas routine <tt>DGEMM</tt> of each pair of vertices.
!>
!>
!> @param[in] IELEM counter of the current triangle/tetrahedron
!> @param[in] ZROE is parameter vector in the NOFVERT vertices of cell \c IELEM
!> @param[in,out] NodRes nodal residual
!> @param[in,out] TSTEP nodal timestep
!> @param[in] NOFVAR is the nof dofs
!> @param[in] VCN Cartesian components of the normals to a face, multiplied by the face area
!> @param[in] NDIM dimension of the space
!> @param[in] NOFVERT nof vertices (= NDIM+1)
!> @param[in] VOLUME area/volume of the current element (triangle,tetrahedron)
!> @param[in,out] STIFD implicit matrix, built only when \c MATRIX_ASSEMBLY .EQV. .TRUE.
!> @param[in] VISCL laminar viscosity
!> @param[in] VISCT turbulent viscosity
!> @param[in] MATRIX_ASSEMBLY \c .EQV. \c .TRUE. if the implicit matrix \c STIFD has to be built, i.e. we are using Picard linearization
!>
!>
!>
!> \author $Author: abonfi $
!> \version $Revision: 1.19 $
!> \date $Date: 2020/03/28 09:52:52 $
!> \bug The calculation of the time-step restriction assumes that the dimensionless kinematic viscosity equals 1
!> \warning The calculation of \c DMAT might be speeded up
!>
      SUBROUTINE VSFLX4(IELEM,ZROE,NodRes,TSTEP,NOFVAR,VCN,NDIM,NOFVERT,
     +                  VOLUME,STIFD,VISCL,VISCT,MATRIX_ASSEMBLY)
C
C     $Id: vsflx4.f,v 1.19 2020/03/28 09:52:52 abonfi Exp $
C
      IMPLICIT NONE 
C
C
C
C     This routine computes the viscous fluxes which appear
C        in the Navier Stokes eqns. in the INTERIOR elements ..
C
C
      INCLUDE 'paramt.h'
      INCLUDE 'bnd.h'
      INCLUDE 'constants.h'
C
      INCLUDE 'dofs.com'
      INCLUDE 'pfcgas.com'
      INCLUDE 'three.com'
      INCLUDE 'transf.com'
      INCLUDE 'turb.com'
      INCLUDE 'visco.com'
C
C
C
C
C     .. Parameters ..
      DOUBLE PRECISION TWOTHIRD
      INTEGER NMAX2
      PARAMETER (TWOTHIRD=TWO/3.D0,NMAX2=NMAX*NMAX)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION VISCL,VISCT,VOLUME
      INTEGER IELEM,NDIM,NOFVAR,NOFVERT
      LOGICAL MATRIX_ASSEMBLY
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION TSTEP(NOFVAR,NOFVERT),
     +       NodRes(NOFVAR,NOFVERT),STIFD(NOFVAR,NOFVAR,NOFVERT,
     +       NOFVERT),VCN(NDIM,NOFVERT),ZROE(NOFVAR,NOFVERT)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION CNST,DIVV,LD,LI,LJ,MI,MJ,MPLUSL,MU,NI,NI_DOT_NJ,
     +                 NJ,TEMPA,TEMPB,TEMPC,U_DOT_NI,U_DOT_NJ,DSUM,DTV
      INTEGER I,IFAIL,IV,IVERT,JV,IADDR,NOFEQN
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION GFLUX(5,VMAX),TAU(3,4),VSFLX(5,VMAX),
     +                 WKSP1(NMAX*VMAX),WORK(NMAX,NMAX)
C     ..
C     .. External Functions ..
      DOUBLE PRECISION DDOT
      EXTERNAL DDOT
C     ..
C     .. External Subroutines ..
      EXTERNAL DCOPY,DGEMM,DGEMV,DINIT
C     ..
C     .. Data statements ..
C
C
      DATA WORK/NMAX2*ZERO/
      DATA TAU,NI,NJ/12*ZERO,2*ZERO/
C     ..
C
      NOFEQN = NDIM + 2
C
C     could compute it directly so as to
C     avoid call overhead, but for the time being .....
C
      CALL GRADPRIM(IELEM,NDIM,NOFEQN)
C
      MU = VISCL + VISCT
      LD = -TWOTHIRD*MU
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
      TAU(1,3) = MU* (GRAD_PRIM(3,3)+GRAD_PRIM(5,1))
      TAU(2,3) = MU* (GRAD_PRIM(4,3)+GRAD_PRIM(5,2))
      TAU(2,1) = TAU(1,2)
      IF(NDIM.EQ.3)THEN
          TAU(3,3) = TWO*MU*GRAD_PRIM(5,3) + LD*DIVV
          TAU(3,1) = TAU(1,3)
          TAU(3,2) = TAU(2,3)
      ENDIF
C
C     Compute the temperature gradient and energy term
C     (the latter is stored in the last column of TAU)
C
      CNST = (VISCL/PRANDTL+VISCT/TPR)/(GM1*UAVG(1))
C
C     \nabla Q = viscl / (gamma - 1.d0) / Pr * \nabla (a^2)
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
C
C     CNST = 1.d0/(d * Re)
C
      CNST = -REINV/NDIM
      DTV = REINV/(NDIM*NDIM*VOLUME)
C
C     Compute the viscous fluxes for each of the NOFVERT vertices
C
      DO 17 IVERT = 1,NOFVERT
C
          DSUM = VCN(1,IVERT)*VCN(1,IVERT)+
     1           VCN(2,IVERT)*VCN(2,IVERT)
          IF(NDIM.EQ.3)DSUM = DSUM + VCN(3,IVERT)*VCN(3,IVERT)
C
C         time step
C
          TSTEP( 1,IVERT) = TSTEP( 1,IVERT) + DTV*DSUM
          TSTEP(IE,IVERT) = TSTEP(IE,IVERT) + DTV*DSUM
          TSTEP(IX,IVERT) = TSTEP(IX,IVERT) + DTV*DSUM
          TSTEP(IY,IVERT) = TSTEP(IY,IVERT) + DTV*DSUM
          IF(NDIM.EQ.3)TSTEP(IZ,IVERT) = TSTEP(IZ,IVERT) + DTV*DSUM
C
          TEMPA = VCN(1,IVERT)*TAU(1,4)+VCN(2,IVERT)*TAU(2,4)
          TEMPB = VCN(1,IVERT)*TAU(1,1)+VCN(2,IVERT)*TAU(2,1)
          TEMPC = VCN(1,IVERT)*TAU(1,2)+VCN(2,IVERT)*TAU(2,2)
          IF(NDIM.EQ.3)THEN
             TEMPA = TEMPA + VCN(3,IVERT)*TAU(3,4)
             TEMPB = TEMPB + VCN(3,IVERT)*TAU(3,1)
             TEMPC = TEMPC + VCN(3,IVERT)*TAU(3,2)
             GFLUX(5,IVERT) = CNST*(
     1       VCN(1,IVERT)*TAU(1,3)+VCN(2,IVERT)*TAU(2,3)+
     2       VCN(3,IVERT)*TAU(3,3) )
          ENDIF
          GFLUX(1,IVERT) = ZERO
          GFLUX(2,IVERT) = CNST*TEMPA
          GFLUX(3,IVERT) = CNST*TEMPB
          GFLUX(4,IVERT) = CNST*TEMPC
C
          NodRes(2,IVERT) = NodRes(2,IVERT) + GFLUX(2,IVERT)
          NodRes(3,IVERT) = NodRes(3,IVERT) + GFLUX(3,IVERT)
          NodRes(4,IVERT) = NodRes(4,IVERT) + GFLUX(4,IVERT)
          IF (NDIM.EQ.3) NodRes(5,IVERT) = NodRes(5,IVERT) + 
     &    GFLUX(5,IVERT)
C
   17 CONTINUE
C         write(6,*)ielem
C     CALL R8Mat_Print('General',' ',NOFVAR,NOFVERT,Nodres,
C    +NOFVAR,'Nodres ',IFAIL)
C
C     Compute the stress vectors
C
      IF (.NOT.MATRIX_ASSEMBLY) RETURN
C
C     contruct a matrix D_{i,j} such that the viscous terms
C     can be written as \sum_j D_{i,j} U_j
C
      MPLUSL = MU + LD
      CNST = -ASQR/GM1 + KINETIC
      TEMPA = (VISCL/PRANDTL+VISCT/TPR)
      TEMPC = REINV/ (NDIM*NDIM*VOLUME*ZAVG(1))
C
C     Loop over vertices ...
C
      DO 10 JV = 1,NOFVERT
C
          LJ = VCN(1,JV)
          MJ = VCN(2,JV)
          IF (NDIM.EQ.3) NJ = VCN(3,JV)
C
          U_DOT_NJ = UAVG(3)*LJ + UAVG(4)*MJ + UAVG(5)*NJ
C
          DO 10 IV = NOFVERT,1,-1
c        DO 10 IV = 1, NOFVERT
C
              LI = VCN(1,IV)
              MI = VCN(2,IV)
              IF (NDIM.EQ.3) NI = VCN(3,IV)
C
              U_DOT_NI = UAVG(3)*LI + UAVG(4)*MI + UAVG(5)*NI
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
C  N.B. the diffusion stiffness matrix is add to the convective one
C       (BETA = 1.d0) in DGEMM
C
              IADDR = (JV-1)*NOFEQN*NOFEQN+1
   20         CALL DGEMM('N','N',NOFEQN,NOFEQN,NOFEQN,TWO*TEMPC,WORK,
     +                   NMAX,DZDU(IADDR),NOFEQN,ONE,STIFD(1,1,IV,JV),
     +                   NOFVAR)
C
C
   10 CONTINUE
C
C     DO IV = 1,NOFVERT
C        DO JV = IV+1,NOFVERT
C        write(6,*)iv,jv
C        CALL R8Mat_Print('General',' ',NOFVAR,NOFVAR,STIFD(1,1,IV,JV),
C    +               NOFVAR,'Diffusion element matrix',IFAIL)
C        CALL R8Mat_Print('General',' ',NOFVAR,NOFVAR,STIFD(1,1,JV,IV),
C    +               NOFVAR,'Diffusion element matrix',IFAIL)
C        pause
C     ENDDO
C     ENDDO
C
      RETURN
C
C     Debugging stuff ...
C     the following will work ONLY if STIFD is 0.d0 upon entry
C
      CALL DCOPY(NOFVAR*NOFVERT,ZROE,1,WKSP1,1)
C     Compute U_i i=1,..,NOFVERT
      CALL PARM_TO_CONS(WKSP1,NDIM,NOFVAR,NOFVERT,.FALSE.,IFAIL)
C
      CALL DINIT(5*VMAX,ZERO,VSFLX,1)
C
C     CALL R8Mat_Print('General',' ',NOFVAR,NOFVERT,WKSP1,
C    +NOFVAR,'Diffusion flux ',IFAIL)
C
C     The "explicitely" calculated viscous flux is compared
C     with the "implicit" one ...
C
      DO 12 IV = 1,NOFVERT
          DO 14 JV = 1,NOFVERT
C
              CALL DGEMV('N',NOFEQN,NOFEQN,ONE,STIFD(1,1,IV,JV),NOFVAR,
     +                   WKSP1((JV-1)*NOFVAR+1),1,ONE,VSFLX(1,IV),1)
C
C        CALL R8Mat_Print('General',' ',NOFVAR,NOFVAR,STIFD(1,1,IV,JV),
C    +   NOFVAR,'Diffusion element matrix',IFAIL)

   14     CONTINUE

          WRITE (6,FMT=*) 'Element # ',IELEM,' vertex # ',IV,'  ',
     +      TEMPC
C
          DO 15 I = 2,NOFEQN
              WRITE (6,FMT=100) GFLUX(I,IV),VSFLX(I,IV),
     +          VSFLX(I,IV)/GFLUX(I,IV)
   15     CONTINUE
C
      pause
   12 CONTINUE
C
C     PAUSE
C
      RETURN

  100 FORMAT (3 (5X,E18.8))

      END
