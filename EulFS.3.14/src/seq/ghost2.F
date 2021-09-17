      SUBROUTINE GHOST2(BC28,ICELNOD,ICELFAC,VFACNOR,XYZDOT,VOL,ZROE, 
     &   VCORG,RCLR,ZCLB,STIFEL,NODRES,NEWRES,VCZ,VCDZ,VCB,
     &   NELEM,NPOIN,NGHOST,NPNOD,NDIM,NOFVERT,NOFVAR,
     &   IBNDFAC,NBFAC,KLIST,VLIST,PetscVecDT,PetscVecRHS,PetscMatA)
 
C
C     Purpose: inflow/outflow boundary conditions;
C              this routine creates a dummy element for each
C              vertex of a boundary face
C
C
C     $Id: ghost2.F,v 1.30 2020/04/23 09:51:45 abonfi Exp $
C
C
C
C
#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscmat.h"
      use petscvec
      use petscmat
      IMPLICIT NONE
C
      Vec PetscVecDT,PetscVecRHS
      Mat PetscMatA
C
      INCLUDE 'paramt.h'
      INCLUDE 'constants.h'
      INCLUDE 'implicit.h'
      INCLUDE 'bnd.h'
      INCLUDE 'bnd.com'
      INCLUDE 'ibc2.com'
      INCLUDE 'ibc8.com'
      INCLUDE 'dofs.com'
      INCLUDE 'stream.com'
      INCLUDE 'flags.com'
      INCLUDE 'periodic.com'
      INCLUDE 'newtoneps.com'
      INCLUDE 'time.com'
      INTEGER NITEMS
      PARAMETER(NITEMS=nVarsInlet*MAXNOFVERT)
C
      INTEGER MY_PE
      COMMON/MPICOM/MY_PE
C
C#define PRINT
C
C     .. Scalar Arguments ..
      INTEGER NBFAC,NDIM,NELEM,NPOIN,NGHOST,NPNOD,NOFVAR,NOFVERT
C     ..
C     .. Subroutine Arguments ..
      EXTERNAL BC28
C     ..
C     .. Array Arguments ..
      INTEGER ICELNOD(NOFVERT,NELEM),ICELFAC(NOFVERT,NELEM),
     &KLIST(NLIST),IBNDFAC(3,NBFAC)
      DOUBLE PRECISION VFACNOR(NDIM,*),XYZDOT(NDIM,*),VOL(NELEM),
     &ZROE(NOFVAR,*),VCORG(NDIM,NPOIN),RCLR(*),ZCLB(NCL,*), ! here ZCLB is transposed
     &VLIST(nVarsInlet,NLIST),STIFEL(NOFVAR,NOFVAR,NOFVERT,NOFVERT), 
     &NODRES(NOFVAR,NOFVERT),NEWRES(NOFVAR,NOFVERT), 
     &VCZ(NOFVAR,NOFVERT),VCDZ(NOFVAR,NOFVERT),VCB(NDIM,NOFVERT)
C     ..
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION VOLUME(MAXTIMLEVS+1),EPS
      INTEGER IVERT,IBC,ICLOR,IELEM,IFACE,N,IFAIL,I,J,LOCA,LOCB,
     +IVAR,JVAR,JVERT,IV,JV,NP
      LOGICAL COMPRESSIBLE,LWK,PLASMA
      INTEGER NERR,IOPT
      CHARACTER*72 ERRMSG
      PARAMETER(NERR = 5, IOPT = 1)
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION VCN(3*MAXNOFVERT),TSTEP(MAXNOFVAR*MAXNOFVERT),
     +DUMMY(MAX_NOFVERT_SQR*MAX_NOFVAR_SQR),PCN(nVarsInlet,MAXNOFVERT),
     &WORK(MAXNOFVAR*MAXNOFVERT)
      INTEGER ICN(MAXNOFVERT),ICN2(MAXNOFVERT)
C
C     VCN(1:NDIM,1:NOFVERT)
C     VCZ(1:NDIM,1:NOFVERT)
C     TSTEP(NOFVAR,1:NOFVERT-1)
C     NODRES(1:NOFVAR,1:NOFVERT-1)
C     STIFEL(1:NOFVAR,1:NOFVAR,1:NOFVERT-1)
C
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
      DATA PCN/NITEMS*ZERO/
C     ..
C     .. Statement Function definitions ..
C
C
      IDX(I,J,N) = (I-1)*N + J
C     ..
C
C
      COMPRESSIBLE = ABS(KAN) .EQ. 4
      PLASMA = ABS(KAN) .EQ. 3
      NP = NPOIN + NGHOST + NPNOD
caldo
      call dinit(nofvar*nofvar*nofvert*nofvert,ZERO,stifel,1)
      call dinit(nofvar*nofvert,ZERO,nodres,1)
      call dinit(nofvar*nofvert,ZERO,newres,1)
C
C     fill in PCN with default values whenever there is no
C     file with an inflow profile
C
C     IF(.NOT.LREAD)THEN
         DO 3 I = 1, NOFVERT
            PCN(1,I) = ONE ! total pressure/ref pressure
            PCN(2,I) = ONE ! total temperature/ref temperature
            PCN(3,I) = POUTLET
            PCN(4,I) = FLOWDIR(1)
            PCN(5,I) = FLOWDIR(2)
            PCN(6,I) = FLOWDIR(3)
    3    CONTINUE
C     ENDIF
C
      DO 100 IFACE = 1,NBFAC ! loop over ALL bndry faces
C
          IELEM = IBNDFAC(1,IFACE)
          IVERT = IBNDFAC(2,IFACE)
          ICLOR = IBNDFAC(3,IFACE)
C
          IBC = ICOLOR(ICLOR,1) ! colours for the mean flow equations
C
C     Apply (inflow/outflow) boundary conditions iff:
C
C     IBC = 0 subsonic inlet
C     IBC = 2 subsonic outlet
C
          IF (.NOT.(IBC.EQ.BC_TYPE_SUBS_OUTLET.OR.
     &              IBC.EQ.BC_TYPE_SUBS_INLET))GOTO 100
C
          CALL CELPTR(IELEM, NELEM, ICELNOD, ICELFAC, VOL, ZROE,
     +                VFACNOR, XYZDOT, NDIM, NOFVERT, NOFVAR, NP, ICN,
     3                VCZ, VCN, VCB, VOLUME)
C
C      ugly!ugly!ugly!ugly!ugly!ugly!ugly!ugly!
C
C      ICN2() keeps the "original" cell to node pointer:
C             NO re-mapping due to periodic bcs.
C
          DO 13 I = 1, NOFVERT
              ICN2(I) = ICELNOD(I,IELEM+NELEM)
   13     CONTINUE
C      ugly!ugly!ugly!ugly!ugly!ugly!ugly!ugly!
C
      INCLUDE 'ptest.f'
C
C     put the vertex facing the boundary face in the
C     NOFVERT-th position
C
          IF(IVERT.NE.NOFVERT)THEN
              CALL DSWAP(NOFVAR,VCZ(1,IVERT),1,VCZ(1,NOFVERT),1)
              LOCA = IDX(IVERT,1,NDIM)
              LOCB = IDX(NOFVERT,1,NDIM)
              CALL DSWAP(NDIM,VCN(LOCA),1,VCN(LOCB),1)
              IF(LALE)CALL DSWAP(NDIM,VCB(1,IVERT),1,VCB(1,NOFVERT),1)
              LOCA = ICN(NOFVERT)
              LOCB = ICN2(NOFVERT)
              ICN(NOFVERT) = ICN(IVERT)
              ICN2(NOFVERT) = ICN2(IVERT)
              ICN(IVERT) = LOCA
              ICN2(IVERT) = LOCB
              LWK = PFLAG(NOFVERT)
              PFLAG(NOFVERT) = PFLAG(IVERT)
              PFLAG(IVERT) = LWK
          ENDIF
!         write(6,*)(icn(i),i=1,nofvert)
!         write(6,*)(icn2(i),i=1,nofvert)
C
C     find inlet profiles for the vertices of the
C     current boundary face
C
#ifdef PRINT
      CALL R8Mat_Print('General',' ',nVarsInlet,NLIST,VLIST,
     +nVarsInlet,'VLIST matrix within ghost2 ',IFAIL)
#endif
          IF(LREAD(1).AND.IBC.EQ.BC_TYPE_SUBS_INLET)
     &    CALL SUBYY(ICN2,NOFVERT,KLIST,
     +               VLIST,nVarsInlet,NLIST,PCN)
!         write(6,*)'IBC = ',iface,ibc
!     CALL R8Mat_Print('General',' ',nVarsInlet,NOFVert,PCN,
!    +nVarsInlet,'PCN before SUBYZ in ghost2 ',IFAIL)
          IF(RAD_EQUI.AND.IBC.EQ.BC_TYPE_SUBS_OUTLET)THEN
              CALL SUBYZ(ICN2,NOFVERT,VCORG,NDIM,
     +               PCN,nVarsInlet,RCLR,ZCLB,NCL)
          ELSE
              DO 23 I = 1, NOFVERT-1
                  PCN(3,I) = POUTLET
   23         CONTINUE
          ENDIF
 
!     CALL R8Mat_Print('General',' ',nVarsInlet,NOFVert,PCN,
!    +nVarsInlet,'PCN after SUBYZ in ghost2 ',IFAIL)
C
#ifdef PRINT
      if(ibc.EQ.BC_TYPE_SUBS_OUTLET)then
      write(6,*)' ibc is ',ibc
      CALL R8Mat_Print('General',' ',nVarsInlet,NOFVERT,PCN,
     +nVarsInlet,'PCN matrix within ghost2 ',IFAIL)
      endif
#endif
C
C     compute residual
C
          CALL DINIT(NOFVAR*(NOFVERT-1),ZERO,TSTEP,1) 
          CALL BC28(IELEM,VCN,VCZ,VCB,PCN,STIFEL,NODRES,TSTEP, NDIM,
     +    NOFVERT,NOFVAR,IBC)
c
#ifdef PRINT
          CALL R8Mat_Print('General',' ',Nofvar,Nofvert,nodres, Nofvar,
     +    'nodres ',IFAIL)
#endif
C
          IF(.NOT.TIMEIMPL)GOTO 50
#if 0
          do ivert=1,nofvert-1
              write(6,*)ielem,ivert,ivert
              CALL R8Mat_Print('General',' ',Nofvar,Nofvar,STIFEL(1,1,
     +        IVERT,1), Nofvar,'Picard matrix ',IFAIL)
          enddo
#endif
          IF(NEWTON)THEN
C
C     Compute Jacobian by finite differences
C
          DO 85 IVERT = 1, NOFVERT-1
              DO 15 IVAR = 1, NOFVAR
                  CALL DCOPY(NOFVAR*NOFVERT,VCZ,1,VCDZ,1)
                  IF(COMPRESSIBLE)THEN
                     CALL PARM_TO_CONS(VCZ (1,IVERT),NDIM,NOFVAR,1,
     &                                 .FALSE.,IFAIL)
                     CALL PARM_TO_CONS(VCDZ(1,IVERT),NDIM,NOFVAR,1,
     &                                 .FALSE.,IFAIL)     
                  ELSEIF(PLASMA)THEN
                     CALL PARM_TO_CONS4Ar(VCZ (1,IVERT),NDIM,NOFVAR,1,
     &          .FALSE.,IFAIL)
                     CALL PARM_TO_CONS4Ar(VCDZ(1,IVERT),NDIM,NOFVAR,1,
     &          .FALSE.,IFAIL)
                  ENDIF
                  EPS = ROOT_MACHINE_EPS*MAX(ABS(VCZ(IVAR,IVERT)),ONE)*
     +            SIGN(ONE,VCZ(IVAR,IVERT))
                  VCDZ(IVAR,IVERT) = VCZ(IVAR,IVERT) + EPS
                  IF(COMPRESSIBLE)THEN
                     CALL CONS_TO_PARM(VCZ (1,IVERT),NDIM,NOFVAR,1,
     &                                 .FALSE.,IFAIL)
                     IF(IFAIL.NE.0)THEN
                        WRITE(ERRMSG(1:64),FMT=105)'GHOST2____',MY_PE,
     &                  "(1)"
                        CALL SETERR(ERRMSG(1:64),64,NERR,IOPT)
                     ENDIF  
                     CALL CONS_TO_PARM(VCDZ(1,IVERT),NDIM,NOFVAR,1,
     &                                 .FALSE.,IFAIL)
                     IF(IFAIL.NE.0)THEN
                        WRITE(ERRMSG(1:64),FMT=105)'GHOST2____',MY_PE,
     &                  "(2)"
                        CALL SETERR(ERRMSG(1:64),64,NERR,IOPT)
                     ENDIF   ! ifail
                  ELSEIF(PLASMA)THEN
                     CALL CONS_TO_PARM4Ar(VCZ (1,IVERT),NDIM,NOFVAR,1,
     &                                    .FALSE.,IFAIL)
                     IF(IFAIL.NE.0)THEN
                        WRITE(ERRMSG(1:64),FMT=105)'GHOST2____',MY_PE,
     &                  "(1)"
                        CALL SETERR(ERRMSG(1:64),64,NERR,IOPT)
                     ENDIF
                     CALL CONS_TO_PARM4Ar(VCDZ(1,IVERT),NDIM,NOFVAR,1,
     &                                    .FALSE.,IFAIL)
                     IF(IFAIL.NE.0)THEN
                        WRITE(ERRMSG(1:64),FMT=105)'GHOST2____',MY_PE,
     &                  "(2)"
                        CALL SETERR(ERRMSG(1:64),64,NERR,IOPT)
                     ENDIF   ! ifail
              ENDIF ! compressible or plasma
C
#if 0
                  write(6,*)'ivert,ivar ibc ',ivert,ivar,ibc
                  CALL R8Mat_Print('General',' ',Nofvar,Nofvert,
     +            vcz(1,1),Nofvar,'Z ',IFAIL)
                  CALL R8Mat_Print('General',' ',Nofvar,Nofvert,
     +            vcdz(1,1),Nofvar,'DZ ',IFAIL)
c
#endif
C
                  CALL BC28(IELEM,VCN,VCDZ,VCB,PCN,DUMMY,NEWRES,WORK,
     +            NDIM,NOFVERT,NOFVAR,IBC)
C
#if 0
                  CALL R8Mat_Print('General',' ',Nofvar,Nofvert,nodres,
     +                             Nofvar,'nodres ',IFAIL)
                  CALL R8Mat_Print('General',' ',Nofvar,Nofvert,newres,
     +                             Nofvar,'newres ',IFAIL)
                  call dcopy(Nofvar*Nofvert,nodres,1,dummy,1)
                  call matsub(dummy,nofvar,newres,nofvar,nofvar,nofvert)
                  CALL R8Mat_Print('General',' ',Nofvar,Nofvert,dummy,
     +                             Nofvar,'newres -nodres ',IFAIL)
#endif
C
                  DO 95 JVAR = 1, NOFVAR
                      DO 95 JVERT = 1, NOFVERT-1
                          STIFEL(JVAR,IVAR,JVERT,IVERT) = -(NEWRES
     +                    (JVAR,JVERT)-NODRES(JVAR,JVERT))/EPS
C
   95             CONTINUE
   15         CONTINUE ! end loop over vars
#if 0
      do 85 jvert=1,nofvert-1
           write(6,*)'I= ',ivert,' J = ',jvert
           CALL R8Mat_Print('General',' ',Nofvar,Nofvar,stifel(1,1,
     &     ivert,jvert), Nofvar,'Jacobian Matrix in Ghost2 (3) ',IFAIL)
#endif
C
   85     CONTINUE ! loop over vertices
C
C     Insert FD jacobian matrix into the global matrix
C
          IF( PERIODIC_MESH .AND. ANNULAR )
     &    CALL ROTATECIJ(STIFEL,NOFVAR,NOFVERT,NOFVERT-1,NOFVERT-1)
C
          DO IVERT = 1,NOFVERT-1
              DO JVERT = 1,NOFVERT-1
                 CALL DCOPY(NOFVAR*NOFVAR,STIFEL(1,1,IVERT,JVERT),1,
     &                      DUMMY,1)
#ifdef MPI
                  CALL MatSetValuesBlockedLocal(PetscMatA,1,ICN(IVERT),
     +            1,ICN(JVERT),DUMMY,ADD_VALUES,IFAIL)
 
#else
                  CALL MatSetValuesBlocked(PetscMatA,1,ICN(IVERT),1,
     +            ICN(JVERT),DUMMY,ADD_VALUES,IFAIL)
#endif
#if 0
                  if(ivert.EQ.jvert)then
                      write(6,*)ivert,jvert
                      CALL R8Mat_Print('General',' ',Nofvar,Nofvar,
     +                stifel(1,1,ivert,jvert), Nofvar,
     +                'Jacobian Matrix in Ghost2 (1)',IFAIL)
                  endif
#endif
            ENDDO
         ENDDO
C
      ELSEIF(PICARD)THEN
C
C     Insert Picard matrix into the global matrix
C
          IF( PERIODIC_MESH .AND. ANNULAR )
     &    CALL ROTATECIJ(STIFEL,NOFVAR,NOFVERT,NOFVERT-1,1)
C
C
          DO 25 IVERT = 1,NOFVERT-1
             CALL DCOPY(NOFVAR*NOFVAR,STIFEL(1,1,IVERT,1),1,
     &                      DUMMY,1)
C
#ifdef MPI
              CALL MatSetValuesBlockedLocal(PetscMatA,1,ICN(IVERT),1,
     +        ICN(IVERT),DUMMY,ADD_VALUES,IFAIL)
 
#else
              CALL MatSetValuesBlocked(PetscMatA,1,ICN(IVERT),1,
     +        ICN(IVERT),DUMMY,ADD_VALUES,IFAIL)
 
#endif
#if 0
              write(6,*)ielem,ivert,ivert
              CALL R8Mat_Print('General',' ',Nofvar,Nofvar,STIFEL(1,1,
     +        IVERT,1), Nofvar,'Jacobian in Ghost2 (2) ',IFAIL)
#endif
C
   25     CONTINUE
      ENDIF ! Picard
C
   50 CONTINUE
C
      INCLUDE 'rotaterhs2.f'
C
         CALL VecSetValuesBlocked(PetscVecDT,NOFVERT-1,ICN,TSTEP,
     +                             ADD_VALUES,IFAIL)
         CALL DCOPY(NOFVAR*(NOFVERT-1),NODRES,1,DUMMY,1)
         CALL VecSetValuesBlocked(PetscVecRHS,NOFVERT-1,ICN,DUMMY,
     +                             ADD_VALUES,IFAIL)
C
  100 CONTINUE
      RETURN
  105 FORMAT('CONS_TO_PARM failed when called from ',A10,' on PE # ',I4,
     &1X,A3)
      END
