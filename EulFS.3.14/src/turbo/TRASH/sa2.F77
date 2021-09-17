      SUBROUTINE SA2(IELEM,VCN,VCZ,NDIM,NOFVERT,
     +         NOFVAR,NTURB,NODRES,
     +         TSTEP,STIFEL,VOLUME,MATRIX_ASSEMBLY,
     +         SCALARSCHEME,MATRIXSCHEME,RWORK,COMPRESSIBLE,
     4         VISCL,VISCT)
C
      IMPLICIT NONE
C
C     Spalart-Allmaras model; transport equation for the
C     MOLECULAR viscosity (rather than kinematic, as in the
C     original version)
C
C
C     $Id: turbcomp.F,v 1.14 2000/06/07 08:12:02 aldo Exp aldo $
C
C
C     define BUDGET to write to file the budget of the
C     turbulent quantities
C
C#define DEBUG 
C#define BUDGET 
C
      INCLUDE 'paramt.h'
C
C
      INCLUDE 'constants'
      INCLUDE 'nloc'
      INCLUDE 'three'
      INCLUDE 'flags.com'
      INCLUDE 'turb.com'
C     INCLUDE 'trip.com'
      INCLUDE 'visco.com'
      INCLUDE 'nboun.com'
C
      double precision s3min,s3max,s3
      double precision s1min,s1max,s1
      double precision s2min,s2max,s2
      common/trashcom/s3min,s3max,s2min,s2max,s1min,s1max
#ifdef BUDGET
      double precision s1,s2,s3,d1(4),d2
      common /sorgente/s1,s2,s3
      common /diffune/d1,d2
      double precision tgt,tdu,tft1
#endif
#ifdef DEBUG
      double precision tmpconv(*),tmpdif1(*),tmpdif2(*),tmpsou1(*),
     &tmpdiff(*),tmpsum(*),tmpsou2(*),tmpsou3(*)
      INTEGER ICN(4)
      DOUBLE PRECISION vort(33000)
      DOUBLE PRECISION stilde(33000)
      DOUBLE PRECISION tfw,tr
      EXTERNAL tfw
#endif
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
      DOUBLE PRECISION RESIDUAL,SOURCE,OMEGA,OMEX,OMEY,OMEZ,
     +                 TS,TTD,NULAM,NUTUR,VI,TD,TCHI,DENSIT
#ifdef BUDGET
      DOUBLE PRECISION BETA
#endif
      INTEGER I,IELEM,IFAIL,INODE,IVERT,J
      LOGICAL MATRIX_ASSEMBLY
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION TSTEP(*),NODRES(NOFVAR,NOFVERT),ADVCT(3),
     +                 STIFC(VMAX*VMAX),STIFD(VMAX*VMAX),VCN(3*VMAX),
     +                 VCP(3*VMAX),VCZ(NOFVAR,NOFVERT),RWORK(*),
     +                 STIFEL(*),WORK(4),DWORK(4)
C     ..
C     .. External Functions ..
      DOUBLE PRECISION DDOT,DNRM2,TBDEST,TBPROD,TBTRIP,TSC,TBDIFF,
     1                 SUTHERLAW,TFV1
      EXTERNAL DDOT,DNRM2,TBDEST,TBPROD,TBTRIP,TSC,TBDIFF,SUTHERLAW,TFV1
C     ..
C     .. External Subroutines ..
      EXTERNAL DINIT,LINEARIZE,TDIFF,VISCO
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
#ifdef BUDGET
#endif
#ifdef DEBUG
      CALL DINIT(NELEM,ZERO,vort,1)
      CALL DINIT(NELEM,ZERO,stilde,1)
#endif
C
C     remark that the working variable Z(NOFVAR) is different from the
C     physical kinatic turbulent viscosity 
C     we shall call:
C     NULAM the laminar kinematic viscosity
C     NUTUR the kinematic eddy viscosity (working variable)
C     the physical laminar viscosity is simply:
C     VISCL = DENSIT * NULAM
C
C
      IF(COMPRESSIBLE)THEN
          VISCL = SUTHERLAW(0.d0,ABAR,ASQR)
          DENSIT = UAVG(1)
C
C Cell averaged vorticity
C
          OMEX = GRAD_PRIM(5,2) - GRAD_PRIM(4,3)
          OMEY = GRAD_PRIM(3,3) - GRAD_PRIM(5,1)
          OMEZ = GRAD_PRIM(4,1) - GRAD_PRIM(3,2)
C
          NUTUR = MAX(0.d0,UAVG(NOFVAR))
C
C tranform from parameter vector into \rho \tilde{\nu}_t}
C we need to copy VCZ(NOFVAR,*) into WORK because we need the
C NOFVERT values of the eddy viscosity in
C contiguous locations, but also because for the compressible
C case we need to transform from parameter to "conserved". 
C
C
          DO 12 IVERT = 1, NOFVERT
             WORK(IVERT) = ZAVG(NOFVAR)*VCZ(1,IVERT)+
     &                     ZAVG(1)*VCZ(NOFVAR,IVERT)
   12     CONTINUE
C
          ADVCT(1) = ZAVG(2)
          ADVCT(2) = ZAVG(3)
          ADVCT(3) = ZAVG(4)
C

      ELSE
          VISCL = ONE
          DENSIT = ONE
C
C Cell averaged vorticity
C
          OMEX = GRAD_PRIM(4,2) - GRAD_PRIM(3,3)
          OMEY = GRAD_PRIM(2,3) - GRAD_PRIM(4,1)
          OMEZ = GRAD_PRIM(3,1) - GRAD_PRIM(2,2)
C
          NUTUR = MAX(0.d0,ZAVG(NOFVAR))
C
          DO 22 IVERT = 1, NOFVERT
             WORK(IVERT) = VCZ(NOFVAR,IVERT)
   22     CONTINUE
C
          ADVCT(1) = DENSIT*UAVG(2)
          ADVCT(2) = DENSIT*UAVG(3)
          ADVCT(3) = DENSIT*UAVG(4)
C
      ENDIF
C
      OMEGA = SQRT(OMEX*OMEX+OMEY*OMEY+OMEZ*OMEZ)
C
C     compute turbulent viscosity
C
      NULAM = VISCL/DENSIT
      TCHI = NUTUR/NULAM
      VISCT = DENSIT*NUTUR*TFV1(TCHI)
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
C   it would be nice to enforce the compiler to inline all these
C   functions
C
C
          SOURCE = TBPROD(NUTUR,NULAM,TS) + 
     2             TBDEST(TD,TS,NUTUR,NULAM) +
     3             TBDIFF(GRAD_PARM,NDIM)
C
          SOURCE = -DENSIT*SOURCE*VOLUME
C
          IF (TTFLAG.EQ.1) THEN
              SOURCE=SOURCE+TBTRIP(TD,RWORK(7),RWORK(4),
     >RWORK(5),RWORK(6))
C
          ENDIF
C
C        the convection speed is the velocity vector
C
          CALL DINIT(NOFVERT,ZERO,NODRES,1)
          CALL DINIT(NOFVERT,ZERO,TSTEP,1)
C
          CALL SCALARSCHEME(IELEM,VCN,RWORK(4),RESIDUAL,SOURCE,WORK,
     +                      TSTEP,DWORK,STIFC,NDIM,NOFVERT,
     &                      MATRIX_ASSEMBLY)
C
C     Diffusion term ...
C
          CALL VISCO(IELEM,DENSIT*(NUTUR+NULAM)*REINV/TPR1,WORK,
     +               DWORK,TSTEP,VCN,STIFD,NDIM,NOFVERT,VOLUME)
C
C we could transpose NODRES: don't know what's better.....
C
         DO 33 IVERT = 1, NOFVERT
            NODRES(NOFVAR,IVERT) = DWORK(IVERT)
   33    CONTINUE
C
          IF (MATRIX_ASSEMBLY) THEN
             STOP 'Picard NOT implemented'
          ENDIF
C
      RETURN
C
#ifdef DEBUG
      call SOLZNE('vorticity.dat',vort,1,NELEM,'w')
      call SOLZNE('stilde.dat',stilde,1,NELEM,'w')
      STOP
#endif


#ifdef BUDGET
      do 43 i = 1,nofvert
         inode = icn(i)+1
         beta =-nodres(i)/(source+residual)
         beta =1./3.
         tmpconv(inode) = tmpconv(inode) + beta * residual
         tmpsou1(inode) = tmpsou1(inode) + beta * s1
         tmpsou2(inode) = tmpsou2(inode) + beta * s2
         tmpsou3(inode) = tmpsou3(inode) + beta * s3
   43 continue
#endif
#ifdef BUDGET
         s1=TBPROD(NUTUR,NULAM,TS)*volume
         s2=TBDEST(TD,TS,NUTUR,NULAM)*volume
c        s3=TBDIFF(GRAD_PARM,ndim)*VOLUME  
         s3=TBTRIP(TD,rwork(7),rwork(4),rwork(5),rwork(6))*VOLUME  
c        write(6,*)td,rwork(7),s3,tgt(tdu(rwork(4),rwork(5),rwork(6)))
c        write(6,*)ielem,rwork(7),tgt(tdu(rwork(4),rwork(5),rwork(6)))
c        write(6,*)ielem,rwork(7),tft1(td,rwork(7),
c    &   tdu(rwork(4),rwork(5),rwork(6)))
#endif
c        s3=TBTRIP(TD,rwork(7),rwork(4),rwork(5),rwork(6))*VOLUME  
c        s1=TBPROD(NUTUR,NULAM,TS)*volume
c        s2=TBDEST(TD,TS,NUTUR,NULAM)*volume
c        s3min=min(s3,s3min) 
c        s3max=max(s3,s3max) 
c        s1min=min(s1,s1min) 
c        s1max=max(s1,s1max) 
c        s2min=min(s2,s2min) 
c        s2max=max(s2,s2max) 
#ifdef DEBUG
         write(6,*)'TD,TS,NUTUR,NULAM',TD,TS,NUTUR,NULAM
         write(6,*)'prod.',NUTUR,NULAM,ts,TBPROD(NUTUR,NULAM,TS)
         write(6,*)'dest.',TD,TS,NUTUR,NULAM,TBDEST(TD,TS,NUTUR,NULAM)
         tr = NUTUR/ (TS* (TK*TD)**2.)*REINV
         TR = (NUTUR/(TS*RE)) / ((TK*TD)**2)
         vort(ielem) = TBTRIP(TD,TTD,ZAVG(2),ZAVG(3),ZAVG(4))
         stilde(ielem) = tfw(tr)
#endif
C
C     Non conservative diffusion term
C
caldo     CALL TDIFF(IELEM,TCB2*REINV/TPR1,WORK,DWORK,TSTEP,VCN,STIFD,
caldo+               NDIM,NOFVERT,VOLUME)
C
#ifdef BUDGET
      do 35 i = 1,nofvert
         inode = icn(i)+1
         tmpdif1(inode) = tmpdif1(inode) + d1(i)
         tmpdif2(inode) = tmpdif2(inode) + d2
   35 continue
#endif
C
C     **** The computation of the residual (r.h.s.) ends HERE ****
C
caldo     write(6,*)'after TDIFF',(nodres(i),i=1,nofvert)
caldo     pause
caldo
      END
