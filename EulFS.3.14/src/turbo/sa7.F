      SUBROUTINE SA7(IELEM,VCN,VCZ,NDIM,NOFVERT,
     +         NOFVAR,NTURB,NODRES,
     +         TSTEP,STIFEL,VOLUME,MATRIX_ASSEMBLY,
     +         SCALARSCHEME,MATRIXSCHEME,RWORK,COMPRESSIBLE,
     4         VISCL,VISCT)
#undef _FUNCT_
#define _FUNCT_ "SA7"
C
C     $Id: sa7.F,v 1.8 2020/04/23 09:56:14 abonfi Exp $
C
      IMPLICIT NONE
C
C     Spalart-Allmaras model; 
C     transport equation for the KINEMATIC viscosity 
C     to be used in the coupled solution approach
C     for compressible flows
C     
C     define DEBUG if you wish to solve a simple convection
C     diffusion equation
C
CCCCC#define DEBUG
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
     +                 VCN(3*MAXNOFVERT),VCZ(NOFVAR,NOFVERT,*),RWORK(*),
     +                 STIFEL(nofvar,nofvar,nofvert,nofvert),VOLUME(*)
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
cxxx  visct = 0.d0
cxxx  return
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
Caldo NUTUR = MAX(0.d0,ZAVG(NOFVAR))
C
C     taking the max bw. 0 and NUTUR seems to have
C     a de-stabilizing effect on convergence
C
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
      TS = TSC(OMEGA,TD,NUTUR,NULAM)
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
              SOURCE=SOURCE+TBTRIP(TD,RWORK(7),RWORK(4),
     >RWORK(5),RWORK(6))
C
          ENDIF
C
CXXXX     pause
CXXXX     write(6,*)source,(NUTUR+NULAM)*REINV/TPR1
C
C
C        the convection speed is the velocity vector
C
#ifdef DEBUG
      CALL NS_SCHEME(IELEM,VCN,UAVG(3),CELRES,ZERO,WORK,
     +                      TAUX,DWORK,BETA,STIFC,NDIM,NOFVERT,
     +                      MATRIX_ASSEMBLY)
          CALL VISCO(IELEM,WORK,DWORK,TAUX,NOFVAR,VCN,NDIM,
     +               NOFVERT,VOLUME(1),STIFD,0.01,
     +               ZERO,MATRIX_ASSEMBLY)
#else
      CALL SCALARSCHEME(IELEM,VCN,UAVG(3),CELRES,SOURCE,WORK,
     +                  TAUX,DWORK,BETA,STIFC,NDIM,NOFVERT,
     +                  MATRIX_ASSEMBLY)
C
C     Diffusion term ...
C
          CALL VISCO(IELEM,WORK,DWORK,TAUX,NOFVAR,VCN,NDIM,
     +               NOFVERT,VOLUME(1),STIFD,(NUTUR+NULAM)*REINV/TPR1,
     +               ZERO,MATRIX_ASSEMBLY)
#endif
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
!     pause
      CALL FLUSH(6)
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
#ifdef PRINT
!     write(6,FMT="(I6,4(1X,F12.6))")ielem,visct,viscl,uavg(1)
!     CALL FLUSH(6)
#endif
c
#ifdef DEBUG
      VISCT = 0.d0
#endif
C
      RETURN
      END
