      SUBROUTINE TURBO(IELEM,VCN,VCZ,VCB,NDIM,NOFVERT,
     +         NTURB,NODRES,
     +         TSTEP,STIFEL,VOLUME,MATRIX_ASSEMBLY,
     +         SCALARSCHEME,MATRIXSCHEME,RWORK)
C
      IMPLICIT NONE
C
C     Driver routine for the Spalart-Allmaras model
C     un-coupled approach
C     this routine is used for both compressible and INcompressible flows
C
C     $Id: turbcomp.F,v 1.25 2020/04/23 09:55:35 abonfi Exp $
C
C
C     define BUDGET to write to file the budget of the
C     turbulent quantities
C
C#define DEBUG 
C#define BUDGET 
C#define PRINT 
C
C
C     .. Parameters ..
      INCLUDE 'paramt.h'
      INCLUDE 'constants.h'
C
C     ..
C     .. Common blocks ..
C
      INCLUDE 'nloc.com'
      INCLUDE 'three.com'
      INCLUDE 'flags.com'
      INCLUDE 'nboun.com'
      INCLUDE 'time.com'
      INCLUDE 'turb.com'
      INCLUDE 'visco.com'
C
#ifdef BUDGET
      double precision s1,s2,s3,d1(4),d2,s4
      common /diffune/d1,d2
      double precision tgt,tdu,tft1
      include 'tdebug.com'
#endif
#ifdef DEBUG
      DOUBLE PRECISION tfw,tr
      EXTERNAL tfw
#endif
C
C     ..
C     .. Scalar Arguments ..
      INTEGER NDIM,NTURB,NOFVERT
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION TSTEP(*),NODRES(*),
     +                 VCZ(NTURB,NOFVERT),RWORK(*),STIFEL(*),
     3                 VCB(NDIM,*),VOLUME(*)
C     ..
C     .. Subroutine Arguments ..
      EXTERNAL SCALARSCHEME,MATRIXSCHEME
C     ..
C     .. Arrays in Common ..
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION RESIDUAL,SOURCE,
     +                 TS,TV,VI,TD
#ifdef BUDGET
      DOUBLE PRECISION BETA
      INTEGER I,INODE
#endif
      INTEGER IELEM,IVERT
      LOGICAL MATRIX_ASSEMBLY
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION STIFC(MAX_NOFVERT_SQR),STIFD(MAX_NOFVERT_SQR),
     +                 VCN(3*MAXNOFVERT)
C     ..
C     .. External Functions ..
      DOUBLE PRECISION TBDEST,TBPROD,TBTRIP,TSC,TBDIFF
      EXTERNAL         TBDEST,TBPROD,TBTRIP,TSC,TBDIFF
C     ..
C     .. External Subroutines ..
      EXTERNAL DINIT,LINEARIZE,TDIFF,VISCO
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC MAX,SQRT
C     ..
C     .. Equivalences ..
C     ..
C
C
C
C
#ifdef BUDGET
#endif
#ifdef DEBUG
C     CALL DINIT(NELEM,ZERO,vort,1)
C     CALL DINIT(NELEM,ZERO,stilde,1)
#endif
C
C
C     COMPUTES THE GRADIENT OF THE turbulent viscosity
C
      CALL LINEARIZE(IELEM,LALE,VCN,VCB,NDIM,NOFVERT,VCZ,NTURB,
     2               VOLUME(1))
C     write(6,*)'elmt. ',ielem
C     write(6,*)GRAD_PARM(1,1),GRAD_PARM(1,2),GRAD_PARM(1,3)
C     CALL R8Mat_Print('General',' ',NTURB,NOFVERT,VCZ,NTURB,
C    +         'Nodal values ',IFAIL)
C
      TD = RWORK(1)
      VI = RWORK(2)
C
C ****************************************************
C          computation of the source term
C ****************************************************
C
C     the cell averaged turbulent viscosity is the aritmetic mean
C     of the nodal turbulent viscosity
C
          TV = ZAVG(1)
C
C     compute vorticity
C
          TS = TSC(RWORK(3),TD,TV,VI)
C
!     write(6,*)td,RWORK(3),ts,tv
!     if(myflag.NE.0)write(6,*)ielem,(RWORK(j),j=1,7)
!     pause
          SOURCE = TBPROD(TV,VI,TS) + TBDEST(TD,TS,TV,VI)
C
#ifdef BUDGET
c        write(6,*)td,rwork(7),s3,tgt(tdu(rwork(4),rwork(5),rwork(6)))
c        write(6,*)ielem,rwork(7),tgt(tdu(rwork(4),rwork(5),rwork(6)))
c        write(6,*)ielem,rwork(7),tft1(td,rwork(7),
c    &   tdu(rwork(4),rwork(5),rwork(6)))
#endif
c        s3=TBTRIP(TD,rwork(7),rwork(4),rwork(5),rwork(6))*VOLUME  
c        s1=TBPROD(TV,VI,TS)*volume
c        s2=TBDEST(TD,TS,TV,VI)*volume
c        s3min=min(s3,s3min) 
c        s3max=max(s3,s3max) 
c        s1min=min(s1,s1min) 
c        s1max=max(s1,s1max) 
c        s2min=min(s2,s2min) 
c        s2max=max(s2,s2max) 
#ifdef DEBUG
         write(6,*)'TD,TS,TV,VI',TD,TS,TV,VI
         write(6,*)'prod.',TV,VI,ts,TBPROD(TV,VI,TS)
         write(6,*)'dest.',TD,TS,TV,VI,TBDEST(TD,TS,TV,VI)
         tr = TV/ (TS* (TK*TD)**2.)*REINV
         TR = (TV/(TS*RE)) / ((TK*TD)**2)
C        vort(ielem) = TBTRIP(TD,TTD,ZAVG(2),ZAVG(3),ZAVG(4))
C        stilde(ielem) = tfw(tr)
#endif
C
          IF (TTFLAG.EQ.1) THEN
              SOURCE=SOURCE+TBTRIP(TD,RWORK(7),RWORK(4),
     >RWORK(5),RWORK(6))
C
          ENDIF
C
caldo
          SOURCE = SOURCE + TBDIFF(GRAD_PARM,NDIM)
caldo
          SOURCE = -SOURCE*VOLUME(1)
C
          CALL DINIT(NOFVERT,ZERO,NODRES,1)
          CALL DINIT(NOFVERT,ZERO,TSTEP,1)
C
C        the convection speed is the velocity vector
C
C        STIFD is being used as a dummy for BETA
C
C
          CALL SCALARSCHEME(IELEM,VCN,RWORK(4),RESIDUAL,SOURCE,VCZ,
     +                      TSTEP,NODRES,STIFD,STIFC,NDIM,NOFVERT,
     &                      MATRIX_ASSEMBLY)
C
#ifdef BUDGET
         s1=TBPROD(TV,VI,TS)*volume
         s2=TBDEST(TD,TS,TV,VI)*volume
         IF(TTFLAG.EQ.1)THEN
            s3=TBTRIP(TD,rwork(7),rwork(4),rwork(5),rwork(6))*VOLUME(1)
         ELSE
            s3=0.d0
         ENDIF
         s4=TBDIFF(GRAD_PARM,ndim)*VOLUME(1)
      do 43 i = 1,nofvert
         inode = kcn(i)+1
         beta =-nodres(i)/(source+residual)
! write(6,*)ielem,kcn(i),beta
         beta =1./3.
         tmpconv(inode) = tmpconv(inode) + beta * residual
         tmpsou1(inode) = tmpsou1(inode) + beta * s1
         tmpsou2(inode) = tmpsou2(inode) + beta * s2
         tmpsou3(inode) = tmpsou3(inode) + beta * s3
         tmpdif2(inode) = tmpdif2(inode) + beta * s4
   43 continue
#endif
C
C     Diffusion term ...
C
          CALL VISCO(IELEM,VCZ,NODRES,TSTEP,NTURB,VCN,NDIM,NOFVERT,
     1               VOLUME(1),STIFD,(TV+VI)*REINV/TPR1,ZERO,
     2               MATRIX_ASSEMBLY)
C
#ifdef PRINT
      write(6,*)'IELEM = ',ielem
      write(6,*)(nodres(i),i=1,nofvert)
#if 1
!     CALL R8Mat_Print('General',' ',NTURB,NOFVERT,NODRES,NTURB,
!    +         'Nodal update (segregated) ',IFAIL)
!     write(6,*)'TD,TS,TV,SOU ',td,ts,tv,source
!     write(6,*)'RWK = ',(rwork(i),i=3,5)
      write(6,*)td,ts,tv,source
      write(6,*)(rwork(i),i=3,5)
      write(6,*)'prod = ',TBPROD(TV,VI,TS)
      write(6,*)'dest = ',TBDEST(TD,TS,TV,VI)
      write(6,*)'diff = ',TBDIFF(GRAD_PARM,NDIM)
      write(6,*)'celr = ',RESIDUAL
#endif
!     pause
#endif
C
C
C     Non conservative diffusion term
C
caldo     CALL TDIFF(IELEM,TCB2*REINV/TPR1,VCZ,NODRES,TSTEP,VCN,STIFD,
caldo+               NDIM,NOFVERT,VOLUME(1))
C
#ifdef BUDGET
!     do 35 i = 1,nofvert
!        inode = kcn(i)+1
!        tmpdif1(inode) = tmpdif1(inode) + d1(i)
!  35 continue
#endif
C
C     **** The computation of the residual (r.h.s.) ends HERE ****
C
caldo     write(6,*)'after TDIFF',(nodres(i),i=1,nofvert)
caldo     pause
C
          IF (MATRIX_ASSEMBLY) THEN
          DO 34 IVERT = 1, NOFVERT*NOFVERT
              STIFEL(IVERT) = STIFD(IVERT)-STIFC(IVERT) 
   34     CONTINUE
          ENDIF
C
#ifdef DEBUG
C     call SOLZNE('vorticity.dat',vort,1,NELEM,'w')
C     call SOLZNE('stilde.dat',stilde,1,NELEM,'w')
C     STOP
#endif
C
!     write(6,FMT="(I6,2(1X,F12.6))")ielem,zavg(1),TBPROD(TV,VI,TS),
!    &TBDEST(TD,TS,TV,VI),TBDIFF(GRAD_PARM,NDIM)
      RETURN
      END
