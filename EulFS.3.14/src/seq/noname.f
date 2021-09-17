      SUBROUTINE NONAME(IELEM,VCN,VCB,VCZ,ZTURB,NDIM,NOFVERT,NOFVAR,
     +                 NTURB,VOLUME,PICARD,VISCOUS,COMPRESSIBLE,
     +                 EulerModel,TurbulenceModel,NSModel,
     +                 ScalarScheme,MatrixScheme,TModelScheme,
     +                 NodRes,TSTEP,ELTMAT)
C
C     $Id: noname.f,v 1.10 2020/03/28 09:42:33 abonfi Exp $
C
C     this subroutine tries to provide a unified interface
C     for the inviscid,turbulent and viscous components of the
C     eqns. while at the same time keeping them in different routines
C
C     given the nodal values of the dependent vars.
C     for the current cell
C     the current subroutine will update
C     -)NodRes
C     -)ELTMAT
C     -)TSTEP
C
      IMPLICIT NONE
C
      INCLUDE "time.h"
      INCLUDE "time.com"
      INCLUDE "flags.com"
C
C     input:
C     -----
      INTEGER IELEM,NDIM,NOFVERT,NOFVAR,NTURB
C
      LOGICAL PICARD,VISCOUS,COMPRESSIBLE
C
      DOUBLE PRECISION ZTURB(*),VOLUME(*),
     +VCZ(NOFVAR,NOFVERT),VCN(NDIM,NOFVERT),VCB(NDIM,NOFVERT)
C
      EXTERNAL EulerModel,TurbulenceModel,NSModel,
     2ScalarScheme,MatrixScheme,TModelScheme
C
C     output:
C     ------
      DOUBLE PRECISION TSTEP(NOFVAR*NOFVERT),
     +ELTMAT(NOFVAR,NOFVAR,NOFVERT,NOFVERT), NodRes(NOFVAR,NOFVERT)
C
C     local:
C
      DOUBLE PRECISION RWORK(10)
      DOUBLE PRECISION viscl,visct
      integer ifail
C
      CALL LINEARIZE(IELEM,LALE,VCN,VCB,NDIM,NOFVERT,
     +               VCZ,NOFVAR,VOLUME(1))
C
C     we would like to remove PARM2PRIM
C
      IF(COMPRESSIBLE)CALL PARM2PRIM(NDIM,IELEM)
C
C     EulerModel() will compute the nodal residual due
C     to the inviscid terms of the equations, updating:
C     -)NodRes
C     -)ELTMAT
C     -)TSTEP
C     valid routines are:
C     ------------------
C     SCALAR     scalar convection
C     EulerII    compressible Euler, Hyperbolic-Elliptic Splitting, supersonic 2D
C     EulerIIbis compressible Euler, Hyperbolic-Elliptic Splitting, general
C     EulerVII   compressible Euler, un-preconditioned, symmetrysing
C                variables
C     EulerVIII  incompressible Euler, un-preconditioned
C     EulerIX    incompressible Euler, Hyperbolic-Elliptic Splitting, general
C     EulerXI    compressible Euler, un-preconditioned, conserved variables
C
      CALL EulerModel(IELEM,VCN,VCB,VCZ,NDIM,NOFVERT,NOFVAR,
     +                NTURB,NodRes,TSTEP,ELTMAT,VOLUME,PICARD,
     +                ScalarScheme,MatrixScheme)
C
C
caldo
cxxxx call SCALAR2(IELEM,VCN,VCZ,NDIM,NOFVERT,NOFVAR,-1,
cxxxx+                  NodRes,TSTEP,ELTMAT,VOLUME,PICARD,
cxxxx+                  SCALARSCHEME,MATRIXSCHEME)
caldo
      IF(VISCOUS)THEN 
C
C     TurbulenceModel() returns both laminar and turbulent 
C     viscosity; in the coupled approach for RANS will
C     also solve the turbulence transport equation
C     valid routines are:
C
C     TDUMMY  scalar diffusion problems
C     LAMINAR laminar flows
C     EVALTSA Spalart-Allmaras model, segregated approach
C     SCALAR2 solves an uncoupled advection equation
C     SA7     Spalart-Allmaras model, compressible flow
C             equations, coupled approach
C     SA8     Spalart-Allmaras model, incompressible flow
C             equations, coupled approach
C
C
C     CALL TurbulenceModel(ZTURB,NTURB,NOFVERT,arg4,arg5,COMPRESSIBLE,
C    1                     VISCL,VISCT)
C
      call TurbulenceModel(IELEM,VCN,VCZ,NDIM,NOFVERT,NOFVAR,NTURB,
     &                  NodRes,TSTEP,ELTMAT,VOLUME,PICARD,
     &                  TModelScheme,MATRIXSCHEME,ZTURB,COMPRESSIBLE,
     &                  VISCL,VISCT)
C
C     I wanted to test the FD jacobian for the viscous
C     terms of the incompressible flow eqns.
C     the set the flag TEST_JACOBIAN in setupRHS()
C     in incompressible flow the diffusive terms are linear
C     so PICARD and NEWTON should give the same results
C     and they do!
C     I need to set eltmat and nodres to 0.d0
C
!         write(6,*)ielem,viscl,visct
!         CALL R8Mat_Print('General',' ',Nofvar,Nofvert,VCZ,
!    +    Nofvar,'Nodal values Matrix inside noname ',IFAIL)
!         CALL R8Mat_Print('General',' ',1,Nofvert,ZTURB,
!    +    1,'v_t values Matrix inside noname ',IFAIL)
!         CALL R8Mat_Print('General',' ',Nofvar,Nofvert,NodRes,
!    +    Nofvar,'REsidual Matrix inside noname ',IFAIL)
!         pause
C
cxxxx
caldo call dinit ((nofvar*nofvert)**2,0.d0,eltmat,1)
caldo call dinit ((nofvar*nofvert),0.d0,nodres,1)
cxxxx
C
C     compute viscous fluxes
C
C     valid routines are:
C
C     VSFLX2 incompressible flows
C     VSFLX4 compressible flows
C     DIFF   scalar diffusion
C
C     write(6,*)ielem,viscl,visct
      CALL NSModel(IELEM,VCZ,NodRes,TSTEP,NOFVAR,VCN,NDIM,
     +             NOFVERT,VOLUME,ELTMAT,VISCL,VISCT,PICARD)
      ENDIF
C
      IF( LAPLACE )THEN
         CALL POISSON(IELEM,VCN,VCZ,NDIM,NOFVERT,NOFVAR,
     +                NTURB,NodRes,TSTEP,ELTMAT,VOLUME,.FALSE.)
      ENDIF
C
C
      RETURN
      END
