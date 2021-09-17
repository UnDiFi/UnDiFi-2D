      SUBROUTINE GHOST(BC5,ICELNOD,ICELFAC,VFACNOR,VOL,ZROE,
     +                 STIFEL,NODRES,NEWRES,VCZ,VCDZ,
     +                 NELEM,NDIM,NOFVERT,NOFVAR,
     +                 IBNDFAC,NBFAC,DT,RHS,A)
C
C     Purpose: inviscid wall boundary conditions;
C              this routine creates a dummy element for each
C              vertex of a boundary face
C
#define NEW_STUFF
C
C     $Id: ghost.F,v 1.16 2000/11/14 16:47:34 aldo Exp abonfi $
C     $Header: /home1/navier/abonfi/EulFS.0.10.12/src/seq/RCS/ghost.F,v 1.16 2000/11/14 16:47:34 aldo Exp abonfi $
C
C
      IMPLICIT NONE
C
C
#include "include/finclude/petsc.h"
#include "include/finclude/petscvec.h"
#include "include/finclude/petscmat.h"
C
      Vec DT,RHS
      Mat A
C
      INCLUDE 'paramt.h'
      INCLUDE 'constants'
      INCLUDE 'implicit.h'
C
      INCLUDE 'nloc'
      INCLUDE 'bnd.h'
      INCLUDE 'bnd'
      INCLUDE 'stream.com'
      INCLUDE 'flags.com'
C
C#define PRINT
C
C
      INCLUDE 'three'
C
C
C     .. Scalar Arguments ..
      INTEGER NBFAC,NDIM,NELEM,NOFVAR,NOFVERT,NPOIN
C     ..
C     .. Subroutine Arguments ..
C
      INTEGER ICELNOD(NOFVERT,NELEM),ICELFAC(NOFVERT,NELEM)
      DOUBLE PRECISION VFACNOR(NDIM,*),VOL(NELEM),
     +                 ZROE(NOFVAR,*),
     +                 STIFEL(NOFVAR,NOFVAR,NOFVERT,NOFVERT),
     +                 NODRES(NOFVAR,NOFVERT),NEWRES(NOFVAR,NOFVERT),
     +                 VCZ(NOFVAR,NOFVERT),VCDZ(NOFVAR,NOFVERT)
      INTEGER IBNDFAC(3,NBFAC) 
C     ..
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION VOLUME,EPS,ROOT_MACHINE_EPS
      INTEGER I,IBC,ICLOR,IELEM,IFACE,N,
     +        IFAIL,IVERT,J,LOCA,LOCB,IVAR,JVAR,JVERT
      LOGICAL COMPRESSIBLE
C
      EXTERNAL BC5
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION VCN(12),TSTEP(4),DUMMY(400)
C
C     VCN(1:NDIM,1:NOFVERT)
C     VCZ(1:NDIM,1:NOFVERT)
C     TSTEP(1:NOFVERT-1)
C     NODRES(1:NOFVAR,1:NOFVERT-1)
C     STIFEL(1:NOFVAR,1:NOFVAR,1:NOFVERT-1)
C
Caldo
      double precision zin(5),zout(5),tmp_array(5)
      common/aaa/zin,zout
Caldo
      logical swapped
C
      INTEGER ICN(VMAX)
C     ..
C     .. External Functions ..
      DOUBLE PRECISION DDOT,DNRM2
      INTEGER ICYCL
      EXTERNAL DDOT,DNRM2,ICYCL
C     ..
C     .. External Subroutines ..
      EXTERNAL CELPTR,DAXPY,DCOPY,DGEMV,DSWAP
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,IABS,REAL
C     ..
C     ..
C     .. Statement Functions ..
      INTEGER IDX
C     ..
C     ..
C     .. Statement Function definitions ..
C
C
      IDX(I,J,N) = (I-1)*N + J
C     ..
C
C
      COMPRESSIBLE = ABS(KAN) .EQ. 4
      ROOT_MACHINE_EPS=1.d-07
C
      DO 100 IFACE = 1,NBFAC
C
C
          IELEM = IBNDFAC(1,IFACE)
          IVERT = IBNDFAC(2,IFACE)
          ICLOR = IBNDFAC(3,IFACE)
C
          IBC = ICOLOR(ICLOR)
C
C     Do NOT apply (inviscid) boundary conditions if:
C
C     IBC = 1 supersonic inlet
C     IBC = 3 supersonic outlet
C     IBC = 6 viscous wall
C     IBC = 7 prescribed velocity profile (maybe I should)
C
#ifdef NEW_STUFF
          IF (.NOT.(IBC.EQ.2.OR.IBC.EQ.8.OR.IBC.EQ.5))GOTO 100
#else
          IF (IBC.NE.5)GOTO 100
#endif
C
      CALL CELPTR(IELEM, ICELNOD, ICELFAC, VOL,
     +   ZROE, VFACNOR, NDIM, NOFVERT, NOFVAR, ICN,
     +   VCZ, VCN, VOLUME)
C
      IF(IVERT.NE.NOFVERT)THEN
         CALL DSWAP(NOFVAR,VCZ(1,IVERT),1,VCZ(1,NOFVERT),1)
         LOCA = IDX(IVERT,1,NDIM)
         LOCB = IDX(NOFVERT,1,NDIM)
         CALL DSWAP(NDIM,VCN(LOCA),1,VCN(LOCB),1)
         LOCA = ICN(NOFVERT)
         ICN(NOFVERT) = ICN(IVERT)
         ICN(IVERT) = LOCA
      ENDIF
C
Caldo
#ifdef NEW_STUFF
      if(ibc.eq.5)then
C    CALL DCOPY(NOFVAR,zin,1,VCZ(1,NOFVERT),1)
            CALL ghost8(VCZ(1,NOFVERT),tmp_array,flowdir,ndim,nofvar)
C           CALL DCOPY(NOFVAR,zin,1,tmp_array,1)
      endif
      if(ibc.eq.2)then
            CALL DCOPY(NOFVAR,zout,1,VCZ(1,NOFVERT),1)
	    CALL DCOPY(NOFVAR,zout,1,tmp_array,1)
      endif
C     if(ibc.eq.5)then
C          CALL DCOPY(NOFVAR,U_INFTY,1,VCZ(1,NOFVERT),1)
C          CALL DCOPY(NOFVAR,U_INFTY,1,tmp_array,1)
C     endif
#else
      if(ibc.eq.5)CALL DCOPY(NOFVAR,U_INFTY,1,VCZ(1,NOFVERT),1)
           CALL DCOPY(NOFVAR,U_INFTY,1,tmp_array,1)
#endif
Caldo
C
          CALL BC5(IELEM,VCN,VCZ,STIFEL,NODRES,TSTEP,
     +NDIM,NOFVERT,NOFVAR,tmp_array)
C
      CALL VecSetValues(DT,NOFVERT-1,ICN,TSTEP,ADD_VALUES,IFAIL)
      CALL VecSetValuesBlocked(RHS,NOFVERT-1,ICN,NODRES,
     +                      ADD_VALUES,IFAIL)
C
      IF(.NOT.TIMEIMPL)GOTO 100
#if 0
      do ivert=1,nofvert-1
      write(6,*)ielem,ivert,ivert
      CALL X04CAF('General',' ',Nofvar,Nofvar,STIFEL(1,1,IVERT,1),
     +            Nofvar,'Picard matrix ',IFAIL)
      enddo
#endif
      IF(PICARD)GOTO 50
C
C     Compute Jacobian by finite differences
C
      DO 15 IVERT = 1, NOFVERT-1
          DO 15 IVAR = 1, NOFVAR
             CALL DCOPY(NOFVAR*NOFVERT,VCZ,1,VCDZ,1)
             IF(COMPRESSIBLE)THEN
                CALL PARM_TO_CONS(VCZ(1,IVERT),NDIM,NOFVAR,1)
                CALL PARM_TO_CONS(VCDZ(1,IVERT),NDIM,NOFVAR,1)
             ENDIF
             EPS = ROOT_MACHINE_EPS*MAX(ABS(VCZ(IVAR,IVERT)),1.d0)*
     +             SIGN(1.d0,VCZ(IVAR,IVERT))
             VCDZ(IVAR,IVERT) = VCZ(IVAR,IVERT) + EPS
             IF(COMPRESSIBLE)THEN
                CALL CONS_TO_PARM(VCZ(1,IVERT),NDIM,NOFVAR,1)
                CALL CONS_TO_PARM(VCDZ(1,IVERT),NDIM,NOFVAR,1)
             ENDIF
C
#if 0
      CALL X04CAF('General',' ',Nofvar,Nofvert,vcz(1,IVERT),
     +            Nofvar,'Z ',IFAIL)
      CALL X04CAF('General',' ',Nofvar,Nofvert,vcdz(1,IVERT),
     +            Nofvar,'DZ ',IFAIL)
#endif
C
             CALL BC5(IELEM,VCN,VCDZ,DUMMY,NEWRES,TSTEP,
     +                   NDIM,NOFVERT,NOFVAR,tmp_array)
C
C     The residual sent to vertex IVERT only depends
C     on the node itself
C
             DO 15 JVAR = 1, NOFVAR
                 STIFEL(JVAR,IVAR,IVERT,1) =
     +           -(NEWRES(JVAR,IVERT)-NODRES(JVAR,IVERT))/EPS
C
   15 CONTINUE
C
   50 CONTINUE
      DO 25 IVERT = 1,NOFVERT-1 
C
#ifdef MPI
          CALL MatSetValuesBlockedLocal(A,1,ICN(IVERT),1,ICN(IVERT),
     +                                  STIFEL(1,1,IVERT,1),
     +                                  ADD_VALUES,IFAIL)
#else
          CALL MatSetValuesBlocked(A,1,ICN(IVERT),1,ICN(IVERT),
     +                             STIFEL(1,1,IVERT,1),
     +                             ADD_VALUES,IFAIL)
#endif
#if 0
      write(6,*)ielem,ivert,ivert
      CALL X04CAF('General',' ',Nofvar,Nofvar,STIFEL(1,1,IVERT,1),
     +            Nofvar,'Jacobian ',IFAIL)
#endif
C
   25 CONTINUE
C
  100 CONTINUE
      RETURN
      END

      subroutine ghost8(zin,zstar,flowdir,ndim,nofvar)
      IMPLICIT NONE
      integer ndim,nofvar
      double precision zin(nofvar),zstar(nofvar),flowdir(ndim)
      double precision uavg(5)
C
      INCLUDE'constants'
      integer nerr,iopt
      double precision asqr,abar,mach,machsqr,kinetic,pres,dens,
     +temp
      character*67 errmsg
C
C
      IF(NDIM.EQ.3)
     &UAVG(5) = ZIN(5)/ZIN(1) ! z componenet of the velocity vector
      UAVG(4) = ZIN(4)/ZIN(1) ! y componenet of the velocity vector
      UAVG(3) = ZIN(3)/ZIN(1) ! x componenet of the velocity vector
      UAVG(2) = ZIN(2)/ZIN(1) ! Total Enthalpy
      UAVG(1) = ZIN(1)*ZIN(1) ! Density
C
      KINETIC = UAVG(3)*UAVG(3) + UAVG(4)*UAVG(4)
      IF( NDIM .EQ. 3 )KINETIC = KINETIC + UAVG(5)*UAVG(5) 
      KINETIC = HALF * KINETIC
      ASQR = GM1 * ( UAVG(2) - KINETIC )
      IF( ASQR .LT. ZERO )THEN
         WRITE(ERRMSG(1:67),FMT=333)ASQR,-10
  333 FORMAT('PARM2PRIM',1X,'Negative averaged sound speed ',F7.3,
     &       ' in element ',I8)
         NERR = 6
         IOPT = 1
         CALL SETERR(ERRMSG(1:67),67,NERR,IOPT)
      ENDIF
      ABAR      = SQRT(ASQR)
c     MACHSQR   = DMAX1( TWO * KINETIC / ASQR, EPSQR_STAG) ! Barth & Deconinck
      MACHSQR   = TWO * KINETIC / ASQR
      MACH      = SQRT( MACHSQR )

      TEMP = ONE+HALF*GM1*MACHSQR
      DENS = TEMP**(-ONE/GM1)
      PRES = ONE/TEMP**GOGM1
      ABAR = SQRT(GAM*PRES/DENS)
      ZSTAR(1) = SQRT(DENS)
      ZSTAR(2) = ZSTAR(1) * GOGM1
      ZSTAR(3) = ZSTAR(1) * ABAR * MACH * FLOWDIR(1)
      ZSTAR(4) = ZSTAR(1) * ABAR * MACH * FLOWDIR(2)
C
      RETURN
      END
C
