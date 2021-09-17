      SUBROUTINE TSETUP(IELEM,NELEM,ZROE,NDIM,NOFVERT,NOFVAR,NP,
     >                  COMPRESSIBLE,RWORK)
C
C     $Id: tsetup.f,v 1.8 2020/04/23 09:55:35 abonfi Exp $
C
      IMPLICIT NONE
C
C     ..
C     .. Arrays in Common ..
      DOUBLE PRECISION DSTAK(1)
      INTEGER ISTAK(1)
C     .. Common blocks ..
      COMMON /CSTAK/DSTAK
C     .. Equivalences ..
      EQUIVALENCE (DSTAK(1),ISTAK(1))
C
      include 'paramt.h'
      include 'constants.h'
      include 'bnd.h'
      include 'three.com'
      include 'stream.com'
      include 'time.com'
      include 'turb.com'
      include 'nloc.com'
C
      INTEGER IELEM,NELEM,NDIM,NOFVERT,NOFVAR,NP
      INTEGER IFAIL,ICN(MAXNOFVERT)
      LOGICAL COMPRESSIBLE
      DOUBLE PRECISION ZROE(*),RWORK(*)
      DOUBLE PRECISION VCZ(MAXNOFVAR*MAXNOFVERT),VCN(3*MAXNOFVERT),
     2                 VCB(3*MAXNOFVERT)
      DOUBLE PRECISION UX,UY,UZ,OMEX,OMEY,OMEZ,OME
      DOUBLE PRECISION VISCL,TD,TTD,VOLUME(MAXTIMLEVS+1)
      DOUBLE PRECISION SUTHERLAW
      INTEGER IVERT
      EXTERNAL SUTHERLAW
C
C**********************************************
C
C     flow variables
C
C**********************************************
C
          CALL CELPTR(IELEM,NELEM,ISTAK(LCELNOD),ISTAK(LCELFAC),
     +    DSTAK(LVOL),ZROE,DSTAK(LFACNOR),DSTAK(LXYZDOT),NDIM,
     +    NOFVERT, NOFVAR, NP, ICN, VCZ, VCN, VCB, VOLUME)
C
C     COMPUTES THE GRADIENT OF THE flow VARIABLES
C
          CALL LINEARIZE(IELEM,LALE,VCN,VCB,NDIM,NOFVERT,
     +                   VCZ,NOFVAR,VOLUME(1))
C
          IF(COMPRESSIBLE)THEN
              CALL PARM2PRIM(NDIM,IELEM)
              CALL GRADPRIM(IELEM,NDIM,NDIM+2)
C
C     kinematic viscosity
C
              VISCL = SUTHERLAW(M_INFTY,ABAR,ASQR)/UAVG(1)
              UX = UAVG(3) 
              UY = UAVG(4) 
              UZ = UAVG(5) 
C
C Cell averaged vorticity
C
              IF(NDIM.EQ.3)THEN
                OMEX = GRAD_PRIM(5,2) - GRAD_PRIM(4,3)
                OMEY = GRAD_PRIM(3,3) - GRAD_PRIM(5,1)
              ELSE
                OMEX = ZERO
                OMEY = ZERO
              ENDIF
              OMEZ = GRAD_PRIM(4,1) - GRAD_PRIM(3,2)
          ELSE
              VISCL = 1.d0
              UX = ZAVG(2)
              UY = ZAVG(3) 
              UZ = ZAVG(4) 
C
C Cell averaged vorticity
C
              IF(NDIM.EQ.3)THEN
                OMEX = GRAD_PARM(4,2) - GRAD_PARM(3,3)
                OMEY = GRAD_PARM(2,3) - GRAD_PARM(4,1)
              ELSE
                OMEX = ZERO
                OMEY = ZERO
              ENDIF
              OMEZ = GRAD_PARM(3,1) - GRAD_PARM(2,2)
          ENDIF
C
          OME = SQRT(OMEX*OMEX+OMEY*OMEY+OMEZ*OMEZ)
C
C     compute a cell averaged wall distance
C
          TD = ZERO
          DO 10 IVERT = 1,NOFVERT
              TD = TD + DSTAK(LTD+ICN(IVERT))
   10     CONTINUE
          TD=TD/NOFVERT
          IF(TTFLAG.EQ.1)THEN 
C
C     compute a cell averaged wall distance
C     on periodic grids, we should rather address
C     the original cell to vertex pointer
C     which is stored in ICELNOD(*,IELEM+NELEM)
C
              TTD = ZERO
              DO 20 IVERT = 1,NOFVERT
                  TTD = TTD + DSTAK(LTTD+ICN(IVERT))
   20         CONTINUE
              TTD=TTD/NOFVERT
          ELSE
              TTD=ZERO
          ENDIF
C
          RWORK(1) = TD
          RWORK(2) = VISCL
          RWORK(3) = OME
          RWORK(4) = UX
          RWORK(5) = UY
          RWORK(6) = UZ
          RWORK(7) = TTD
C
C     CALL R8Mat_Print('General',' ',NOFVAR,3,GRAD_PARM,NMAX,
C    +         'gradZ ',IFAIL)
C         write(*,*)(rwork(ivert),ivert=1,7),omex,omey,omez
C
          RETURN
          END 
