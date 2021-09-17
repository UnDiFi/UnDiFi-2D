!> \brief \b SA
!> \par Purpose
!>
!> This routine discretizes the SA (standard) version of the Spalart-Allmaras turbulence model
!> transport equation for the KINEMATIC viscosity 
!> to be used in the coupled solution approach for compressible flows
!>
!> \f[
!> \frac{D\tilde{\nu}}{Dt} = c_{b_1} \left[ 1 - f_{t_2} \right] \tilde{S} \tilde{\nu}
!> - \frac{1}{Re} \left[ c_{w_1} f_{w} - \frac{c_{b_1}}{\kappa^2} f_{t_2} \right] \left[ \frac{\tilde{\nu}}{d} \right]^2
!> + \frac{1}{\sigma Re} \left[ \nabla \cdot \left( \left( \nu + \tilde{\nu} \right) \nabla \tilde{\nu} \right) + c_{b_2} \left( \nabla \tilde{\nu} \right)^2 \right]
!> \f]
!>
!> Nondimensionalization also leads to the redefinition of the \f$ \tilde{S} \f$ and \f$r\f$ variables:
!> \f[
!> \tilde{S} \equiv \frac{U_{0}}{L_{0}} \left( \Omega + \frac{1}{Re_{0}} \frac{\tilde{\nu} f_{v2}} {\kappa^2 d^2 } \right)
!> \f]
!> \f[
!> r \equiv \frac{1}{Re_{0}} \, \frac{\tilde{\nu}}{\tilde{S} \kappa^2 d^2 }
!> \f]
!>
!> As far as the limiting of \f$ \tilde{S} \f$ is concerned, we use:
!> \f[
!> \tilde{S} = \max \left( \tilde{S}, 0.3\,\Omega \right)
!> \f]
!>
!> which is referred to as "limiting (b)" in http://turbmodels.larc.nasa.gov/spalart.html
!>
!> @param[in] IELEM is the current simplicial element
!> @param[in] VCN the NDIM cartesian component of the inward face normal to all NOFVERT vertices, scaled by its measure
!> @param[in] VCZ is the set of NOFVAR dependent variables in the NOFVERT vertices of the current cell
!> @param[in] NDIM is the dimension of the space
!> @param[in] NOFVERT is number of vertices of the simplicial elements (= NDIM+1)
!> @param[in] NOFVAR is the number of dofs in the gridpoints, also the leading dimension of VCZ, STIFEL, NODRES, TSTEP
!> @param[in] NTURB is the number of turbulent transport equations being solved
!> @param[out] NODRES are the signals sent to each of the NOFVERT vertices
!> @param[in,out] TSTEP is the elemental contribution to the (inverse of the) time step
!> @param[out] STIFEL the elemental contribution to the implicit matrix (computed only if MATRIX_ASSEMBLY is true)
!> @param[in] VOLUME is the area/volume of the current cell
!> @param[in] MATRIX_ASSEMBLY when set == .TRUE. the STIFEL matrix will be assembled
!> @param[in] SCALARSCHEME is the subroutine used to discretize the turbulent transport equation
!> @param[in] MATRIXSCHEME is the subroutine used to discretize the turbulent transport coupled equations
!> @param[in] RWORK is a work array
!> @param[in] COMPRESSIBLE is logical flag, currently unused
!> @param[out] VISCL is the elemental (cell averaged) molecular laminar viscosity
!> @param[out] VISCT is the elemental (cell averaged) molecular turbulent viscosity
!> \author $Author: abonfi $
!> \version $Revision: 1.2 $
!> \date $Date: 2020/04/23 09:55:35 $
!> \warning MATRIXSCHEME is currently unused
!> \warning MATRIX_ASSEMBLY == .TRUE. is NOT supported; it is among the arguments for compatibility with similar calls
!> \warning STIFEL is currently unused; it is among the arguments for compatibility with similar calls
!> \warning COMPRESSIBLE is currently unused; it is among the arguments for compatibility with similar calls
!> \warning RWORK is currently unused; it is among the arguments for compatibility with similar calls
      SUBROUTINE SA(IELEM,VCN,VCZ,NDIM,NOFVERT,
     +         NOFVAR,NTURB,NODRES,
     +         TSTEP,STIFEL,VOLUME,MATRIX_ASSEMBLY,
     +         SCALARSCHEME,MATRIXSCHEME,RWORK,COMPRESSIBLE,
     4         VISCL,VISCT)
C
C     $Id: sa.F,v 1.2 2020/04/23 09:55:35 abonfi Exp $
C
      IMPLICIT NONE
C
C
CCC#define PRINT
C
C
C     .. Parameters ..
C     ..
      INCLUDE 'paramt.h'
      INCLUDE 'constants.h'
      INCLUDE 'time.h'
C
C     .. Common blocks ..
C     ..
C
      INCLUDE 'three.com'
      INCLUDE 'flags.com'
      INCLUDE 'turb.com'
      INCLUDE 'time.com'
      INCLUDE 'visco.com'
      INCLUDE 'nboun.com'
      INCLUDE 'sa.com'
C
C     Sets residual and local timestep to zero
C
C     .. Parameters ..
C     ..
C     .. Scalar Arguments ..
      INTEGER NDIM,NELEM,NOFVAR,NTURB,NOFVERT
      DOUBLE PRECISION VISCL,VISCT
      LOGICAL COMPRESSIBLE
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION TSTEP(NOFVAR,*),NODRES(NOFVAR,NOFVERT),
     +                 VCN(NDIM,NOFVERT),VCZ(NOFVAR,NOFVERT,*),RWORK(*),
     +                 STIFEL(NOFVAR,NOFVAR,NOFVERT,NOFVERT),VOLUME(*)
C     ..
C     .. Subroutine Arguments ..
      EXTERNAL SCALARSCHEME,MATRIXSCHEME
C     ..
C     .. Arrays in Common ..
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION CELRES,SOURCE,OMEGA,OMEX,OMEY,OMEZ,
     +                 TS,TBDIFFX,GRADNRM,NULAM,NUTUR
      INTEGER I,J,IELEM,IFAIL,JADD,IVERT
      LOGICAL MATRIX_ASSEMBLY
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION STIFC(MAX_NOFVERT_SQR),STIFD(MAX_NOFVERT_SQR),
     +                 WORK(MAXNOFVERT*MAXTIMLEVS),DWORK(MAXNOFVERT),
     2                 TAUX(MAXNOFVERT),BETA(MAXNOFVERT)
C     ..
C     .. External Functions ..
      DOUBLE PRECISION TBDEST,TBPROD,TBTRIP,TSCSA,
     1                 TFV1,SUTHERLAW
      EXTERNAL TBDEST,TBPROD,TBTRIP,TSCSA,TFV1,SUTHERLAW
C     ..
C     .. External Subroutines ..
      EXTERNAL DINIT,TDIFF,VISCO
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,MAX,SIGN,SQRT
C     ..
C
C
      CALL DINIT(NOFVERT,ZERO,TAUX,1)
C
C     this will be called again within the viscous routines
C
      CALL GRADPRIM(IELEM,NDIM,NDIM+2)
C
C     molecular viscosity
C
      VISCL = SUTHERLAW(ZERO,ABAR,ASQR)
C
C     kinematic viscosity
C
      NULAM = VISCL / UAVG(1)
C
C Cell averaged vorticity
C
      IF( NDIM .EQ. 3 )THEN
          OMEX = GRAD_PRIM(5,2) - GRAD_PRIM(4,3)
          OMEY = GRAD_PRIM(3,3) - GRAD_PRIM(5,1)
      ELSE
          OMEX = ZERO
          OMEY = ZERO
      ENDIF
      OMEZ = GRAD_PRIM(4,1) - GRAD_PRIM(3,2)
C
      OMEGA = SQRT(OMEX*OMEX+OMEY*OMEY+OMEZ*OMEZ)
C     write(6,*)'w = ',omex,omey,omez
C
C
C     taking the max bw. 0 and NUTUR seems to have
C     a de-stabilizing effect on convergence
C
Caldo NUTUR = MAX(0.d0,ZAVG(NOFVAR))
      NUTUR = ZAVG(NOFVAR)
C
C     CALL R8Mat_Print('General',' ',NOFVAR,3,GRAD_PRIM,NMAX,
C    +         'Gradiente ',IFAIL)
C     CALL R8Mat_Print('General',' ',NOFVAR,NOFVERT*NTIMLEVS,VCZ,NOFVAR,
C    +         'Nodal values ',IFAIL)
C
      DO 22 J = 1, NTIMLEVS
         JADD = (J-1)*NOFVERT
         DO 22 IVERT = 1, NOFVERT
             WORK(JADD+IVERT) = VCZ(NOFVAR,IVERT,J)
   22 CONTINUE
C
C
C ****************************************************
C          computation of the source term
C ****************************************************
C
C
C     compute "modified" vorticity
C
      TS = TSCSA(OMEGA,TD,NUTUR,NULAM)
      TS = MAX(TS,0.3d0*OMEGA)
#ifdef PRINT
!     write(6,*)'w = ',omega,' ts ',ts,' nut ',nutur
#endif
C
C   it would be nice to enforce the compiler to inline all these
C   functions
C
      GRADNRM = GRAD_PARM(NOFVAR,1)**2 + 
     >          GRAD_PARM(NOFVAR,2)**2
      IF(NDIM.EQ.3)GRADNRM = GRADNRM + GRAD_PARM(NOFVAR,3)**2
C
      TBDIFFX = TCB2/TPR1*REINV * GRADNRM
C
      SOURCE = TBPROD(NUTUR,NULAM,TS) + 
     2         TBDEST(TD,TS,NUTUR,NULAM) +
     3         TBDIFFX
C
      SOURCE = -SOURCE*VOLUME(1)
C
          IF (TTFLAG.EQ.1) THEN
              STOP 'trip term NOT implemented '
          ENDIF
C
C        the convection speed is the velocity vector
C
      CALL SCALARSCHEME(IELEM,VCN,UAVG(3),CELRES,SOURCE,WORK,
     +                  TAUX,DWORK,BETA,STIFC,NDIM,NOFVERT,
     +                  MATRIX_ASSEMBLY)
C
C     Diffusion term ...
C
          CALL VISCO(IELEM,WORK,DWORK,TAUX,NOFVAR,VCN,NDIM,
     +               NOFVERT,VOLUME(1),STIFD,(NUTUR+NULAM)*REINV/TPR1,
     +               ZERO,MATRIX_ASSEMBLY)
C
      IF(LTIME)THEN
          DTVOL = DELT/VOLUME(1)
          CALL UNSTEADY1(WORK,BETA,WORK,NTURB,DWORK,STIFC,VOLUME,
     1                   NTURB,NDIM,NOFVERT,MMTYPE,.FALSE.)
!    1                   NOFVERT,MM_LUMPED,.FALSE.)
      ENDIF 
C
C we could transpose NODRES: don't know what's better.....
C
         DO 33 IVERT = 1, NOFVERT
            NODRES(NOFVAR,IVERT) = DWORK(IVERT)
            TSTEP(NOFVAR,IVERT) = TSTEP(NOFVAR,IVERT) + TAUX(IVERT)
   33    CONTINUE
C
#ifdef PRINT
      write(6,*)'IELEM = ',ielem
      write(6,*)(work(ivert),ivert=1,nofvert)
      write(6,*)(dwork(ivert),ivert=1,nofvert)
#if 0
      CALL R8Mat_Print('General',' ',NOFVAR,NOFVERT,NODRES,NOFVAR,
     +         'Nodal update ',IFAIL)
!     write(6,*)'TD,TS,NT,SOU = ',td,ts,nutur,source
!     write(6,*)'OME, U,V ',omega,uavg(3),uavg(4)
      write(6,*)td,ts,nutur,source
      write(6,*)omega,uavg(3),uavg(4)
      write(6,*)'prod = ',TBPROD(NUTUR,NULAM,TS)
      write(6,*)'dest = ',TBDEST(TD,TS,NUTUR,NULAM)
      write(6,*)'diff = ',TBDIFFX
      write(6,*)'celr = ',CELRES
#endif
#endif
C
      IF (MATRIX_ASSEMBLY) THEN
             WRITE(6,*) 'Picard NOT implemented'
cxxx         STOP 'Picard NOT implemented'
      ENDIF
C
C     compute "real" turbulent viscosity
C
      VISCT = NUTUR*TFV1(NUTUR/NULAM)*UAVG(1)
C
      RETURN
      END
