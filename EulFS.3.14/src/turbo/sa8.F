      SUBROUTINE SA8(IELEM,VCN,VCZ,NDIM,NOFVERT,
     +         NOFVAR,NTURB,NODRES,
     +         TSTEP,STIFEL,VOLUME,MATRIX_ASSEMBLY,
     +         SCALARSCHEME,MATRIXSCHEME,RWORK,COMPRESSIBLE,
     4         VISCL,VISCT)
C
C     $Id: sa8.F,v 1.9 2020/04/23 09:56:14 abonfi Exp $
C
      IMPLICIT NONE
C
C     Spalart-Allmaras model; 
C     transport equation for the KINEMATIC viscosity 
C     to be used in the coupled solution approach
C     for incompressible flows
C     
C
C     .. Parameters ..
      INCLUDE 'paramt.h'
      INCLUDE 'time.h'
      INCLUDE 'constants.h'
C
C
      INCLUDE 'three.com'
      INCLUDE 'flags.com'
      INCLUDE 'turb.com'
      INCLUDE 'visco.com'
      INCLUDE 'nboun.com'
      INCLUDE 'time.com'
      INCLUDE 'sa.com'
C
C     ..
C     .. Scalar Arguments ..
      INTEGER NDIM,NELEM,NOFVAR,NTURB,NOFVERT
      DOUBLE PRECISION VISCL,VISCT
      LOGICAL COMPRESSIBLE
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION TSTEP(NOFVAR,*),NODRES(NOFVAR,NOFVERT),
     +                 VCN(*),VCZ(NOFVAR,NOFVERT,*),RWORK(*), 
     3                 STIFEL(NOFVAR,NOFVAR,NOFVERT,NOFVERT),
     4                 VOLUME(*)
C     ..
C     .. Subroutine Arguments ..
      EXTERNAL SCALARSCHEME,MATRIXSCHEME
C     ..
C     .. Arrays in Common ..
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION CELRES,SOURCE,OMEGA,OMEX,OMEY,OMEZ,
     +                 TS,TCHI,TBDIFFX,GRADNRM
      INTEGER I,IELEM,IFAIL,IVERT,J,JADD,NN,MM
      LOGICAL MATRIX_ASSEMBLY
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION STIFC(MAX_NOFVERT_SQR),STIFD(MAX_NOFVERT_SQR),
     +                 WORK(MAXNOFVERT*MAXTIMLEVS),DWORK(MAXNOFVERT),
     3                 TAUX(MAXNOFVERT),BETA(MAXNOFVERT)
C     ..
C     .. External Functions ..
      DOUBLE PRECISION TBDEST,TBPROD,TBTRIP,TSC,
     1                 TFV1
      EXTERNAL TBDEST,TBPROD,TBTRIP,TSC,TFV1
C     ..
C     .. External Subroutines ..
      EXTERNAL DINIT,TDIFF,VISCO
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC MAX,SQRT
C     ..
C     .. Common blocks ..
C     ..
C     .. Equivalences ..
C     ..
C
C
C
CCC#define DEBUG
C
C
C
      CALL DINIT(NOFVERT,ZERO,TAUX,1)
      IF(MATRIX_ASSEMBLY)THEN 
         MM = NOFVERT*NOFVERT
         CALL DINIT(MM,ZERO,STIFC,1) 
         CALL DINIT(MM,ZERO,STIFD,1) 
      ENDIF
C
C
      VISCL = ONE
C
C Cell averaged vorticity
C
      IF( NDIM .EQ. 3 )THEN
          OMEX = GRAD_PARM(4,2) - GRAD_PARM(3,3)
          OMEY = GRAD_PARM(2,3) - GRAD_PARM(4,1)
      ELSE
          OMEX = ZERO
          OMEY = ZERO
      ENDIF
      OMEZ = GRAD_PARM(3,1) - GRAD_PARM(2,2)
C
      OMEGA = SQRT(OMEX*OMEX+OMEY*OMEY+OMEZ*OMEZ)
C     write(6,*)'w = ',omex,omey,omez
!     wkarea(1) = TD
!     wkarea(2) = VISCL
!     wkarea(3) = OMEGA
!     wkarea(4) = ZAVG(2)
!     wkarea(5) = ZAVG(3)
!     wkarea(6) = ZAVG(4)
!     wkarea(7) = TTD
!     if(myflag.NE.0)write(6,*)ielem,(wkarea(j),j=1,7)
C
C     VISCT = MAX(0.d0,ZAVG(NOFVAR))
      VISCT = ZAVG(NOFVAR)
C
C     CALL R8Mat_Print('General',' ',NOFVAR,3,GRAD_PARM,NMAX,
C    +         'Gradiente ',IFAIL)
C     CALL R8Mat_Print('General',' ',NOFVAR,NOFVERT,VCZ,NOFVAR,
C    +         'Nodal values ',IFAIL)
C
C
      DO 22 J = 1, NTIMLEVS
         JADD = (J-1)*NOFVERT
         DO 22 IVERT = 1, NOFVERT
             WORK(JADD+IVERT) = VCZ(NOFVAR,IVERT,J)
   22 CONTINUE
C
C
C
C ****************************************************
C          computation of the source term
C ****************************************************
C
C
C     compute "modified" vorticity
C
      TS = TSC(OMEGA,TD,VISCT,VISCL)
C     write(6,*)'w = ',omega,' ts ',ts
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
      SOURCE = TBPROD(VISCT,VISCL,TS) + 
     2         TBDEST(TD,TS,VISCT,VISCL) +
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
CXXXX     write(6,*)source,(VISCT+VISCL)*REINV/TPR1
C
#ifdef DEBUG
          source = 0.d0
#endif
C
C        the convection speed is the velocity vector
C
      CALL SCALARSCHEME(IELEM,VCN,ZAVG(2),CELRES,SOURCE,WORK,
     +                      TAUX,DWORK,BETA,STIFC,NDIM,NOFVERT,
     +                      MATRIX_ASSEMBLY)
C
C     Diffusion term ...
C
#ifdef DEBUG
          CALL VISCO(IELEM,WORK,DWORK,TAUX,NOFVAR,VCN,NDIM,
     +               NOFVERT,VOLUME(1),STIFD,REINV,
     +               ZERO,MATRIX_ASSEMBLY)
#else
          CALL VISCO(IELEM,WORK,DWORK,TAUX,NOFVAR,VCN,NDIM,
     +               NOFVERT,VOLUME(1),STIFD,(VISCT+VISCL)*REINV/TPR1,
     +               ZERO,MATRIX_ASSEMBLY)
#endif
C
      IF(LTIME)THEN
          DTVOL = DELT/VOLUME(1)
          CALL UNSTEADY1(WORK,BETA,WORK,NTURB,DWORK,STIFC,VOLUME,
     1                   NTURB,NDIM,NOFVERT,MMTYPE,MATRIX_ASSEMBLY)
!    1                   NOFVERT,1,MATRIX_ASSEMBLY)
      ENDIF 
C
C
C we could transpose NODRES: don't know what's better.....
C
         DO 33 IVERT = 1, NOFVERT
            NODRES(NOFVAR,IVERT) = DWORK(IVERT)
            TSTEP(NOFVAR,IVERT) = TSTEP(NOFVAR,IVERT) + TAUX(IVERT)
   33    CONTINUE
C
!     write(6,*)ielem
!     write(6,*)(dwork(ivert),ivert=1,nofvert)
!     CALL R8Mat_Print('General',' ',NOFVAR,NOFVERT,NODRES,NOFVAR,
!    +         'Nodal update ',IFAIL)
!     write(6,*)ielem,td,ts,visct,source
C     write(6,*)zavg(2),zavg(3),omega
!     write(6,*)'prod = ',TBPROD(VISCT,VISCL,TS)
!     write(6,*)'dest = ',TBDEST(TD,TS,VISCT,VISCL)
!     write(6,*)'diff = ',TBDIFFX
!     write(6,*)'celr = ',CELRES
!     pause
C
      IF (MATRIX_ASSEMBLY) THEN
C
         NN = NOFVAR*NOFVAR
         CALL DAXPY(MM,MONE,STIFC,1,STIFD,1)
!        CALL DSCAL(MM,MONE,STIFC,1) 
         CALL DCOPY(MM,STIFC,1,STIFEL(NOFVAR,NOFVAR,1,1),NN)
!        CALL DINIT(MM,ONE,STIFEL(NOFVAR,NOFVAR,1,1),NN)
!     CALL R8Mat_Print('General',' ',NOFVAR,NOFVAR*NOFVERT*NOFVERT,STIFEL,
!    +         NOFVAR,'Matrice J ',IFAIL)
      ENDIF
C
C     compute "real" turbulent viscosity
C
C     write(6,*)ielem,visct
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
