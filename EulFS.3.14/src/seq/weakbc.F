!> @param[in] ICELNOD Cell to node pointer: \c ICELNOD(i,j) gives the global node number of the i-th vertex of the j-th element
!> @param[in] ICELFAC Cell to face pointer: \c ICELFAC(i,j} gives the global face number of the face opposite the i-th vertex of the j-th element
!> @param[in] VFACNOR Cartesian components of the normals to a face, multiplied by the face area
!> @param[in] XYZDOT the cartesian components of the grid velocity
!> @param[in] VOL area/volume of the simplicial elements (triangles,tetrahedra)
!> @param[in] ZROE Nodal values of the dependent variable;  compressible equations: \f$ Z = \sqrt{\rho} \left( 1, H, \mathbf{u} \right) \f$; incompressible equations: \f$ Z = \left( p, \mathbf{u} \right) \f$
!> @param[in,out] STIFEL elemental stiffness matrix
!> @param[in,out] NODRES 
!> @param[in,out] NEWRES 
!> @param[out] VCZ returns the \c NOFVAR dofs of the NOFVERT vertices of cell \c IELEM; VCZ(1:NOFVAR,*,1) stores the values at time level \c n+1,k ; \c VCZ(1:NOFVAR,*,2) those at tme level n and \c VCZ(1:NOFVAR,*,3) those at time level \c n-1
!> @param[out] VCDZ returns the \c NOFVAR dofs of the NOFVERT vertices of cell \c IELEM; VCZ(1:NOFVAR,*,1) stores the values at time level \c n+1,k ; \c VCZ(1:NOFVAR,*,2) those at tme level n and \c VCZ(1:NOFVAR,*,3) those at time level \c n-1
!> @param[out] VCB returns the \c NDIM Cartesian components of the nodal grid velocities at time \c n+1/2 of the NOFVERT vertices of cell \c IELEM; only if \c LALE is \c .TRUE.
!> @param[in] NELEM is the no. of processor owned elements (triangles/tetrahedra); global number of elements in the uni-processor case
!> @param[in] NPOIN nof interior nodes in the mesh
!> @param[in] NGHOST nof ghost nodes in the mesh
!> @param[in] NPNOD nof periodic nodes in the mesh
!> @param[in] NDIM dimension of the space
!> @param[in] NOFVERT number of vertices per element (=NDIM+1, since only triangles or tetrahedra are allowed)
!> @param[in] NOFVAR number of variables (degrees of freedom) in each meshpoint
!> @param[in] IBNDFAC
!> @param[in] NBFAC number of bndry faces
!> @param[in,out] PetscVecDT Petsc Vec where time-step is stored
!> @param[in,out] PetscVecRHS Petsc Vec where the nodal residual is stored
!> @param[in,out] PetscMatA Petsc Mat where the implicit matrix is stored
!> @param[in] LFLAG
      SUBROUTINE WEAKBC(ICELNOD,ICELFAC,VFACNOR,XYZDOT,VOL,ZROE,STIFEL,
     +                  NODRES,NEWRES,VCZ,VCDZ,VCB,NELEM,NPOIN,NGHOST,
     +                  NPNOD,NDIM,NOFVERT,NOFVAR,IBNDFAC,NBFAC,
     4                  PetscVecDT,PetscVecRHS,PetscMatA,LFLAG)
 
C
C#define PRINT
C#define DEBUG
C     weakly imposes boundary conditions:
C     1) slip-free (inviscid) wall boundary condition
C     2) prescribed flux
C
C     $Id: weakbc.F,v 1.38 2020/04/23 17:28:17 abonfi Exp $
C
#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscmat.h"
      use petscvec
      use petscmat
C
      IMPLICIT NONE
C
      INCLUDE 'paramt.h'
      INCLUDE 'implicit.h'
C
C     PETSc stuff
C
      Vec PetscVecRHS,PetscVecDT
      Mat PetscMatA
C
      INCLUDE 'constants.h'
      INCLUDE 'bnd.h'
C
      INCLUDE 'io.com'
      INCLUDE 'dofs.com'
      INCLUDE 'bnd.com'
      INCLUDE 'chorin.com'
      INCLUDE 'nloc.com'
      INCLUDE 'newtoneps.com'
      INCLUDE 'periodic.com'
      INCLUDE 'time.com'
      INCLUDE 'three.com'
C
      INTEGER MY_PE
      COMMON/MPICOM/MY_PE
#ifdef DEBUG
      COMMON/MYBUG/ICLR
#endif
C
C     FORTRAN stack
C
      DOUBLE PRECISION DSTAK(1)
      COMMON /CSTAK/ DSTAK
      INTEGER ISTAK(1)
      EQUIVALENCE(DSTAK(1),ISTAK(1))
C
      INTEGER NELEM,NPOIN,NGHOST,NPNOD,NBFAC,NDIM,NOFVERT,NOFVAR
      INTEGER I,J,IELEM,LOC,IFAIL,IPOS,LAST,IADDR,NM1,
     +IBC,IBFAC,LOCA,LOCB,K,L,N,LWORK1,LWORK2,ICLR
C
C
      LOGICAL LFLAG(*)
      LOGICAL COMPRESSIBLE,PLASMA
      INTEGER IDX
C
      INTEGER ICELNOD(NOFVERT,NELEM),ICELFAC(NOFVERT,NELEM),
     2IBNDFAC(3,NBFAC)
      DOUBLE PRECISION VFACNOR(NDIM,*),XYZDOT(NDIM,*),VOL(NELEM),
     2ZROE(NOFVAR,*),STIFEL(NOFVAR,NOFVAR,NOFVERT,NOFVERT), 
     3NODRES(NOFVAR,NOFVERT),NEWRES(NOFVAR,NOFVERT),VCZ(NOFVAR,NOFVERT),
     4VCDZ(NOFVAR,NOFVERT),VCB(NDIM,NOFVERT)
C
C
C
C
C     NODRES(1:NOFVAR,1:NOFVERT) is used to accumulate the
C         nodal residual in conserved variables and scatter
C         it to the PetscVecRHS PETSc vector
C     NEWRES(1:NOFVAR,1:NOFVERT) same as NODRES, but for the
C         perturbed dependent variable (only Newton iteration)
C     TSTEP(1:NOFVERT) is used to accumulate the timestep
C         and then scatter it to the PETSc vector PetscVecDT
C
C
      DOUBLE PRECISION VCN(3*VMAX),TSTEP,
     &DUMMY(MAX_NOFVAR_SQR*MAX_NOFVERT_SQR)
      INTEGER ICN(VMAX)
C
      INTEGER NERR,IOPT
      CHARACTER*72 ERRMSG
      PARAMETER(NERR = 5, IOPT = 1)
C
      DOUBLE PRECISION VOLUME(MAXTIMLEVS+1),EPS
      LOGICAL LWK
C
      INTEGER  IVAR,JVAR,IVERT,JVERT,IV,JV,NP
      INTEGER  ICYCL,ISTKGT
      EXTERNAL ICYCL,ISTKGT
      DOUBLE PRECISION DNRM2
      EXTERNAL DNRM2,PRESSI,PRESSC,PRESS4Ar
C
      IDX(I,J,N) = (I-1)*N + J
C
      NP = NPOIN + NGHOST + NPNOD
      NM1 = NOFVERT-1
C
      COMPRESSIBLE = LFLAG(1)
      PLASMA = LFLAG(2)
C
      LWORK1 = ISTKGT((NOFVERT*NOFVAR)**2,4)
      LWORK2 = ISTKGT((NOFVERT*NOFVAR)**2,4)
C
!      CALL R8Mat_Print('General',' ',NOFVAR,NBFLX(2),DSTAK(LBNDFLX(2)),
!     +            NOFVAR,'bndry pointer ',IFAIL)
C
      DO 100 IBFAC = 1,NBFAC
          ICLR = IBNDFAC(3,IBFAC)
          IBC = ICOLOR(ICLR,1) ! colours for the mean flow equations
          IF( .NOT.((IBC .EQ.  BC_TYPE_SLIP_FREE) .OR.
     1              (IBC .EQ.  BC_TYPE_PRESCRIBED_FLUX)))GOTO 100
          IVERT = IBNDFAC(2,IBFAC)
          IELEM = IBNDFAC(1,IBFAC)
C
C
          CALL CELPTR(IELEM, NELEM, ICELNOD, ICELFAC, VOL, ZROE,
     +                VFACNOR, XYZDOT, NDIM, NOFVERT, NOFVAR, NP, ICN,
     3                VCZ, VCN, VCB, VOLUME)
C
C 
      INCLUDE 'ptest.f'
C
          CALL LINEARIZE(-IELEM,LALE,VCN,VCB,NDIM,NOFVERT,VCZ,
     +    NOFVAR,VOLUME(1))
C
          IF(COMPRESSIBLE)CALL PARM2PRIM(NDIM,IELEM)
          IF(PLASMA)CALL PARM2PRIM4Ar(NDIM,IELEM) !!< da cambiare
C
C     Put the vertex opposite the boundary face in the last position
C
          IF(IVERT.NE.NOFVERT)THEN
              CALL DSWAP(NOFVAR,VCZ(1,IVERT),1,VCZ(1,NOFVERT),1)
              LOCA = IDX(IVERT,1,NDIM)
              LOCB = IDX(NOFVERT,1,NDIM)
              CALL DSWAP(NDIM,VCN(LOCA),1,VCN(LOCB),1)
              IF(LALE)CALL DSWAP(NDIM,VCB(1,IVERT),1,VCB(1,NOFVERT),1)
              LOCA = ICN(NOFVERT)
              ICN(NOFVERT) = ICN(IVERT)
              ICN(IVERT) = LOCA
              LWK = PFLAG(NOFVERT)
              PFLAG(NOFVERT) = PFLAG(IVERT)
              PFLAG(IVERT) = LWK
          ENDIF
C
C     normal to the boundary
C
          LOC = (NM1)*NDIM+1
C
          CALL DINIT((NOFVAR*NOFVERT)**2,ZERO,STIFEL,1)
C
          IF( IBC. EQ. BC_TYPE_PRESCRIBED_FLUX)THEN 
C
C         retrieve the prescribed flux
C
              CALL BINSRC(IBFAC,ISTAK(LBNDFLX(2)),NBFLX(2),IPOS,LAST)
              IF(IPOS.EQ.0)THEN
                  WRITE(NOUT,*)'Cannot find matching bndry face in weakb
     &c; iface= ',IBFAC
                  CALL I4Mat_Print('General',' ',3,NBFAC,ISTAK(LBNDFAC),
     +            3,'bndry pointer ',IFAIL)
                  WRITE(6,*)(ISTAK(I),I=LBNDFLX(2),
     &                                  LBNDFLX(2)+NBFLX(2)-1)
                  CALL EXIT(IPOS)
              ELSE
                  IF(.NOT.(IPOS.GT.0.AND.IPOS.LE.NBFLX(2)))THEN
                      write(6,*)ipos
                  ENDIF
              ENDIF ! ipos
              IADDR = (IPOS-1)*NOFVAR + LBNDFLX(1) ! should be nofeqn ?!?!
!             write(6,*)ibfac,' bndry face is ',ipos
!             write(6,*)(DSTAK(ipos),ipos=iaddr,iaddr+nofvar-1)
              IF(COMPRESSIBLE)THEN
                  CALL FLXB4(NDIM,NOFVAR,NOFVERT,STIFEL,DSTAK(IADDR),
     +                       DSTAK(LWORK2),VCZ,VCB,VCN(LOC),NODRES,
     &                       PICARD)
              ELSE
                  STOP 'Unimplemented feature BC_TYPE_PRESCRIBED_FLUX in
     &subr. weakbc'
              ENDIF ! on COMPRESSIBLE
          ELSE ! slip-free boundary condition
              IF(COMPRESSIBLE)THEN
                  CALL FLXW4(NDIM,NOFVAR,NOFVERT,STIFEL,DSTAK(LWORK1),
     +                       DSTAK(LWORK2),VCZ,VCB,VCN(LOC),NODRES,
     &                       PICARD)
                  IF(NDIM.EQ.2)THEN
                      CALL PFORCESimpson(ICLR,NOFVERT,VCN,NDIM,VCZ,
     +                                   NOFVAR,NOFVERT,PRESSC)
                  ELSE
                      CALL PFORCE(ICLR,NOFVERT,VCN,NDIM,VCZ,NOFVAR,
     +                            NOFVERT,PRESSC)
                  ENDIF 
              ELSEIF(PLASMA)THEN
                  CALL FLXW44Ar(NDIM,NOFVAR,NOFVERT,STIFEL,
     +                          DSTAK(LWORK1),DSTAK(LWORK2),VCZ,
     +                          VCN(LOC),NODRES,PICARD)
                  IF(NDIM.EQ.2)THEN
                      CALL PFORCESimpson(ICLR,NOFVERT,VCN,NDIM,VCZ,
     +                                   NOFVAR,NOFVERT,PRESS4Ar)
                  ELSE
                      CALL PFORCE(ICLR,NOFVERT,VCN,NDIM,VCZ,NOFVAR,
     +                            NOFVERT,PRESS4Ar)
                  ENDIF 
              ELSE ! Incomprimibile
                  CALL FLXW2(NDIM,NOFVAR,NOFVERT,STIFEL,DSTAK(LWORK1),
     +                       DSTAK(LWORK2),VCZ,VCB,VCN(LOC),NODRES,
     &                       PICARD)
                  CALL PFORCE(ICLR,NOFVERT,VCN,NDIM,VCZ,NOFVAR,NOFVERT,
     +                        PRESSI)
              ENDIF ! on COMPRESSIBLE
          ENDIF ! on IBC
C
#if 0
          do k=1,nofvert-1
              do l=1,nofvert-1
                  CALL R8Mat_Print('General',' ',NOFVAR,NOFVAR,STIFEL
     +            (1,1,k,l),NOFVAR,'Picard matrix for wall b.c. ',IFAIL)
              enddo
          enddo
#endif
C
          IF(COMPRESSIBLE.OR.PLASMA)THEN
              TSTEP = ABAR * DNRM2(NDIM,VCN(LOC),1) ! should be sqrt(dnrm2)
          ELSE
              TSTEP = BETA * DNRM2(NDIM,VCN(LOC),1)
          ENDIF
C
          IF(.NOT.TIMEIMPL)GOTO 50 ! Explicit timestepping
          IF(NEWTON)THEN
C 
C     Compute Jacobian by finite differences
C
          DO 15 IVERT = 1, NM1
              DO 25 IVAR = 1, NOFVAR
                  CALL DCOPY(NOFVAR*NOFVERT,VCZ,1,VCDZ,1)
C
                  IF(COMPRESSIBLE)THEN
                      CALL PARM_TO_CONS( VCZ(1,IVERT),NDIM,NOFVAR,1,
     &                                   .FALSE.,IFAIL)
                      CALL PARM_TO_CONS(VCDZ(1,IVERT),NDIM,NOFVAR,1,
     &                                   .FALSE.,IFAIL)
                  ELSEIF(PLASMA)THEN ! modificare qui
                      CALL PARM_TO_CONS4Ar( VCZ(1,IVERT),NDIM,NOFVAR,1,
     &                                   .FALSE.,IFAIL)
                      CALL PARM_TO_CONS4Ar(VCDZ(1,IVERT),NDIM,NOFVAR,1,
     &                                   .FALSE.,IFAIL)
                  ENDIF
C
                  EPS = ROOT_MACHINE_EPS*MAX(ABS(VCZ(IVAR,IVERT)),ONE)*
     +            SIGN(ONE,VCZ(IVAR,IVERT))
                  VCDZ(IVAR,IVERT) = VCZ(IVAR,IVERT) + EPS
C
                  IF(COMPRESSIBLE)THEN
                      CALL CONS_TO_PARM(VCZ(1,IVERT),NDIM,NOFVAR,1,
     &                                  .FALSE.,IFAIL)
                      IF(IFAIL.NE.0)THEN
                         WRITE(ERRMSG(1:64),FMT=105)'WEAKBC____',
     &                                              MY_PE,"(1)"
                         CALL SETERR(ERRMSG(1:64),64,NERR,IOPT)
                      ENDIF   ! ifail
                      CALL CONS_TO_PARM(VCDZ(1,IVERT),NDIM,NOFVAR,1,
     &                                  .FALSE.,IFAIL)
                      IF(IFAIL.NE.0)THEN
                         WRITE(ERRMSG(1:64),FMT=105)'WEAKBC____',
     &                                              MY_PE,"(2)"
                         CALL SETERR(ERRMSG(1:64),64,NERR,IOPT)
                      ENDIF   ! ifail
                  ELSEIF(PLASMA)THEN ! modificare qui
                      CALL CONS_TO_PARM4Ar(VCZ(1,IVERT),NDIM,NOFVAR,1,
     &                                     .FALSE.,IFAIL)
                      IF(IFAIL.NE.0)THEN
                         WRITE(ERRMSG(1:64),FMT=105)'WEAKBC____',
     &                                              MY_PE,"(1)"
                         CALL SETERR(ERRMSG(1:64),64,NERR,IOPT)
                      ENDIF   ! ifail
                      CALL CONS_TO_PARM4Ar(VCDZ(1,IVERT),NDIM,NOFVAR,1,
     &                                     .FALSE.,IFAIL)
                      IF(IFAIL.NE.0)THEN
                         WRITE(ERRMSG(1:64),FMT=105)'WEAKBC____',
     &                                              MY_PE,"(2)"
                         CALL SETERR(ERRMSG(1:64),64,NERR,IOPT)
 
                      ENDIF
                  ENDIF
C
                  IF( IBC. EQ. BC_TYPE_PRESCRIBED_FLUX)THEN 
                      IF(COMPRESSIBLE)THEN
                          CALL FLXB4(NDIM,NOFVAR,NOFVERT,DUMMY,
     2                               DSTAK(IADDR),DSTAK(LWORK2),VCDZ,
     3                               VCN(LOC),NEWRES,.FALSE.)
                      ELSE
                         STOP ' Unimplemented feature'
                      ENDIF
                  ELSE ! inviscid wall bndry conditions
                      IF(COMPRESSIBLE)THEN
                          CALL FLXW4(NDIM,NOFVAR,NOFVERT,DUMMY,
     2                               DSTAK(LWORK1),DSTAK(LWORK2),VCDZ,
     3                               VCB,VCN(LOC),NEWRES,.FALSE.)
                      ELSEIF(PLASMA)THEN
                          CALL FLXW44Ar(NDIM,NOFVAR,NOFVERT,DUMMY, 
     2                               DSTAK(LWORK1),DSTAK(LWORK2),VCDZ, 
     3                               VCN(LOC),NEWRES,.FALSE.) ! modifica
                      ELSE ! incomprimibile
                          CALL FLXW2(NDIM,NOFVAR,NOFVERT,DUMMY,
     2                               DSTAK(LWORK1),DSTAK(LWORK2),VCDZ,
     3                               VCB,VCN(LOC),NEWRES,.FALSE.)
                      ENDIF
                  ENDIF
#if 0
                  CALL R8Mat_Print('General',' ',NOFVAR,NOFVert,nodres,
     +            NOFVAR,'nodres for wall b.c. ',IFAIL)
                  CALL R8Mat_Print('General',' ',NOFVAR,NOFVert,newres,
     +            NOFVAR,'newres for wall b.c. ',IFAIL)
#endif
C
                  DO JVAR = 1, NOFVAR
                      DO JVERT = 1, NM1
                          STIFEL(JVAR,IVAR,JVERT,IVERT) = -(NEWRES
     +                    (JVAR,JVERT)-NODRES(JVAR,JVERT))/EPS
                      ENDDO
                  ENDDO
C
   25        CONTINUE
   15     CONTINUE
         ENDIF ! Newton
C
          IF( PERIODIC_MESH .AND. ANNULAR )
     &    CALL ROTATECIJ(STIFEL,NOFVAR,NOFVERT,NM1,NM1)
C
C
         N = NM1*NOFVAR
         DO IV = 1, NM1
            DO JV = 1, NM1
               DO JVAR = 1, NOFVAR
                  DO IVAR = 1, NOFVAR
                     K = (IV-1)*NOFVAR+IVAR
                     L = (JV-1)*NOFVAR+JVAR
                     IADDR = (L-1)*N+K
                     DUMMY(IADDR) = STIFEL(IVAR,JVAR,IV,JV)
                  ENDDO ! IVAR
               ENDDO ! JVAR
#if 0
               IF(IBC.EQ.IBC)THEN
                  WRITE(6,*)IV,JV
                  CALL R8Mat_Print('General',' ',NOFVAR,NOFVAR,
     +            STIFEL(1,1,IV,JV),NOFVAR,'STIFEL ',ifail)
               ENDIF
#endif
            ENDDO ! JV
         ENDDO ! I
C
C     Add the element diffusion stiffness matrix to the global one
C
C     call MatSetValuesBlocked only once per boundary face
C
#ifdef MPI
         CALL MatSetValuesBlockedLocal(PetscMatA,NM1,ICN,NM1,
     +                     ICN,DUMMY,ADD_VALUES,IFAIL)
#else
         CALL MatSetValuesBlocked(PetscMatA,NM1,ICN,NM1,ICN,
     +                     DUMMY,ADD_VALUES,IFAIL)
#endif

C
   50     CONTINUE
C
      INCLUDE 'rotaterhs2.f'
C
C
C     This ugly copy is needed because VecSetValuesBlocked only
C     allows 1D arrays
C
          CALL DCOPY(NM1*NOFVAR,NODRES,1,DUMMY,1)
          CALL VecSetValuesBlocked(PetscVecRHS,NM1,ICN,DUMMY,
     +    ADD_VALUES, IFAIL)
C
C         insert the scalar time step into the PETSc vec
C                                 |
C                                 |
C                                \|/
          CALL DINIT(NM1*NOFVAR,TSTEP,DUMMY,1)
          CALL VecSetValuesBlocked(PetscVecDT,NM1,ICN,DUMMY,
     +       ADD_VALUES, IFAIL)
C
 100  CONTINUE ! end loop over the boundary faces
C
C     release work arrays
C
      CALL ISTKRL(2)
C
      RETURN
  105 FORMAT('CONS_TO_PARM failed when called from ',A10,' on PE # ',I4,
     &1X,A3)
      END
C
