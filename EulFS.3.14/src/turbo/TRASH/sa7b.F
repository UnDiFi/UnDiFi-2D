      SUBROUTINE SA7b(IELEM,VCN,VCZ,NDIM,NOFVERT,
     +         NOFVAR,NTURB,NODRES,
     +         TSTEP,STIFEL,VOLUME,MATRIX_ASSEMBLY,
     +         SCALARSCHEME,MATRIXSCHEME,RWORK,COMPRESSIBLE,
     4         VISCL,VISCT)
C
C     $Id: sa7b.F,v 1.2 2013/01/26 12:01:29 abonfi Exp $
C
      IMPLICIT NONE
C
C     Spalart-Allmaras model; 
C     transport equation for the MOLECULAR viscosity 
C     to be used in the coupled solution approach
C     for compressible flows
C     symmetrising variables, i.e.
C     
C     \sqrt(r) (1,H,u,nu_t) 
C     
C
#define DEBUG
C
C
      INCLUDE 'paramt.h'
C
C
      INCLUDE 'constants.h'
      INCLUDE 'three.com'
      INCLUDE 'flags.com'
      INCLUDE 'turb.com'
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
      DOUBLE PRECISION VOLUME,VISCL,VISCT
      LOGICAL COMPRESSIBLE
C     ..
C     .. Subroutine Arguments ..
      EXTERNAL SCALARSCHEME,MATRIXSCHEME
C     ..
C     .. Arrays in Common ..
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION CELRES,SOURCE,OMEGA,OMEX,OMEY,OMEZ,
     +                 TS,TCHI,TBDIFFX,GRADNRM,NUTUR,NULAM,HELP
      INTEGER I,IELEM,IFAIL,INODE,IVERT,J
      LOGICAL MATRIX_ASSEMBLY
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION TSTEP(NOFVAR,*),NODRES(NOFVAR,NOFVERT),
     +                 STIFC(VMAX*VMAX),STIFD(VMAX*VMAX),VCN(3*VMAX),
     +                 VCZ(NOFVAR,NOFVERT),RWORK(*),
     +                 STIFEL(nofvar,nofvar,nofvert,nofvert),
     +                 WORK(4),DWORK(4),TAUX(VMAX),GRAD_NUTUR(3)
C     ..
C     .. External Functions ..
      DOUBLE PRECISION TBDEST,TBPROD,TBTRIP,TSC,
     1                 TFV1,SUTHERLAW
      EXTERNAL TBDEST,TBPROD,TBTRIP,TSC,TFV1,SUTHERLAW
C     ..
C     .. External Subroutines ..
      EXTERNAL DINIT,TDIFF,VISCO
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,MAX,SIGN,SQRT
C     ..
C     .. Common blocks ..
C     ..
C     .. Equivalences ..
C     ..
C
C
C
      CALL DINIT(NOFVERT,ZERO,TAUX,1)
C
C     NUTUR is nu_t
C
      NUTUR = ZAVG(NOFVAR)/ZAVG(1)
C
C     this will be called again within the viscous routines
C
      CALL GRADPRIM(IELEM,NDIM,NDIM+2)
C
C     kinematic viscosity
C
      VISCL = SUTHERLAW(0.d0,ABAR,ASQR)
      NULAM = VISCL/UAVG(1)
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
C
C     write(6,*)'w = ',omex,omey,omez
C
C     compute the gradient of the kinematic turb. viscosity
C     according to:
C     d/dx (nu_t) = (dz_1/dx - nu_t dz_6/dx ) / z_1
C
C
      GRADNRM = ZERO
      DO 22 I= 1, NDIM
             HELP = ( GRAD_PARM(1,I) -
     &       ZAVG(NOFVAR)*GRAD_PARM(1,I) )/ZAVG(1)
             GRADNRM = GRADNRM + HELP*HELP
             GRAD_NUTUR(I) = HELP
   22 CONTINUE
C
C     CALL X04CAF('General',' ',NOFVAR,3,GRAD_PRIM,NMAX,
C    +         'Gradiente ',IFAIL)
C     CALL X04CAF('General',' ',NOFVAR,NOFVERT,VCZ,NOFVAR,
C    +         'Nodal values ',IFAIL)
C
C     compute working variable in the nodes
C
      DO 23 IVERT = 1, NOFVERT
         WORK(IVERT) = (VCZ(1,IVERT)-NUTUR*VCZ(NOFVAR,IVERT))/
     &   ZAVG(1)
   23 CONTINUE
C
C
C ****************************************************
C          computation of the source term
C ****************************************************
C
C
C     compute "modified" vorticity
C
      TS = TSC(OMEGA,TD,NUTUR,NULAM)
C
C     write(6,*)'w = ',omega,' ts ',ts
C
C   it would be nice to enforce the compiler to inline all these
C   functions
C
C
      TBDIFFX = TCB2/TPR1*REINV * GRADNRM
C
      SOURCE = TBPROD(NUTUR,NULAM,TS) + 
     2         TBDEST(TD,TS,NUTUR,NULAM) +
     3         TBDIFFX
C
C    REM: density multiplied source term
C
      SOURCE = -UAVG(1)*SOURCE*VOLUME
C
          IF (TTFLAG.EQ.1) THEN
              STOP 'trip term NOT implemented '
              SOURCE=SOURCE+TBTRIP(TD,RWORK(7),RWORK(4),
     >RWORK(5),RWORK(6))
C
          ENDIF
C
CXXXX     pause
CXXXX     write(6,*)source,(VISCT+VISCL)*REINV/TPR1
C
C
C        the convection speed is the velocity vector
C
C        Distribute the source term 
C
      CALL NS_SCHEME(IELEM,VCN,UAVG(3),CELRES,SOURCE,WORK,
     +                      TAUX,DWORK,STIFC,NDIM,NOFVERT,
     +                      MATRIX_ASSEMBLY)
C
C        Diffusion term ...
C
#ifdef DEBUG
          CALL VISCO(IELEM,WORK,DWORK,TAUX,NOFVAR,VCN,NDIM,
     +               NOFVERT,VOLUME,STIFD,0.01,
     +               ZERO,MATRIX_ASSEMBLY)
#else
          CALL VISCO(IELEM,WORK,DWORK,TAUX,NOFVAR,VCN,NDIM,
     +               NOFVERT,VOLUME,STIFD,
     +               UAVG(1)*(NUTUR+NULAM)*REINV/TPR1,
     +               ZERO,MATRIX_ASSEMBLY)
#endif
C
C we could transpose NODRES: don't know what's better.....
C
         DO 33 IVERT = 1, NOFVERT
            NODRES(NOFVAR,IVERT) = NODRES(NOFVAR,IVERT) + DWORK(IVERT)
            TSTEP(NOFVAR,IVERT) = TSTEP(NOFVAR,IVERT) + TAUX(IVERT)
   33    CONTINUE
C
!     write(6,*)ielem
!     write(6,*)(dwork(ivert),ivert=1,nofvert)
!     CALL X04CAF('General',' ',NOFVAR,NOFVERT,NODRES,NOFVAR,
!    +         'Nodal update ',IFAIL)
C     write(6,*)ielem,td,ts,visct,source
C     write(6,*)omega,uavg(3),uavg(4)
C     write(6,*)'prod = ',TBPROD(VISCT,VISCL,TS)
C     write(6,*)'dest = ',TBDEST(TD,TS,VISCT,VISCL)
C     write(6,*)'diff = ',TBDIFFX
C     write(6,*)'celr = ',CELRES
C     pause
C
      IF (MATRIX_ASSEMBLY) THEN
          CALL DSCAL(NOFVERT*NOFVERT,-ONE,STIFC,1) 
          CALL DCOPY(NOFVERT*NOFVERT,STIFC,1,STIFEL(NOFVAR,NOFVAR,1,1),
     >NOFVAR*NOFVAR) 
      ENDIF
C
C     compute "real" turbulent viscosity
C
C     write(6,*)ielem,visct
      VISCT = UAVG(1)*NUTUR
      TCHI = VISCT/VISCL
      VISCT = VISCT*TFV1(TCHI)
cxxx
#ifdef DEBUG
      VISCT = 0.d0
c
c     prova con VISCL=0.d0 ed il modello "intero"
c
#endif
cxxx
C
      RETURN
C
      END
