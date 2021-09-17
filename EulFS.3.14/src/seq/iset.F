      SUBROUTINE ISET(NodalBcs, IDX, IBNDFAC, NBFAC, ICELNOD, NOFVERT,
     &NODCODE, NDIM, NOFVAR, NPOIN, NGHOST, NPNOD, NELEM,
     &ROWBGN, COMPRESSIBLE )
C
C     $Id: iset.F,v 1.36 2020/04/23 09:46:19 abonfi Exp $
C
C     This routine creates PETSc Index Sets addressing
C          boundary meshpoints with various types of
C          boundary conditions.
C     Warnings:
C     --------
C     In the parallel case it is important that
C          the routine is called with NPOIN=# of processor
C          owned nodes, since the index sets are local
C          and shoud not address ghost nodes
C
C     NODCODE will address interior nodes as well ghost ones
C     so it must be allocated for NPOIN+NGHOST in the calling routine
C
C     Note that Index sets are created even if empty
C
C     One of the limitations in the current implementation is that
C     if a boundary vertex belongs to two (or more) boundary faces
C     where different types of b.c. are specified the routine won't
C     detect it (though it shouldn't be too difficult to modify it
C     so that it will) but, more importantly, conflicts may arise
C     (or may not) depending on whether the two b.c. affect different
C     variables or not.
C     For instance: take the trailing edge of a flat plate
C     and suppose the t.e. is also located on the outflow boundary
C     where a cnst. pressure boundary condition is specified.
C     This won't conflict with the no-slip boundary
C     condition but WILL conflict if you specify an iso-thermal b.c.
C     on the whole plate (you can give it a try).
C
#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscis.h"
      use petscvec
      use petscis
C
      IMPLICIT NONE
C
      INCLUDE 'bnd.h'
      INCLUDE 'bctype.h'
      INCLUDE 'iset.h'
C
      INCLUDE 'bnd.com'
      INCLUDE 'io.com'
      INCLUDE 'flags.com'
      INCLUDE 'visco.com'
      INCLUDE 'ibc8.com'
      INCLUDE 'time.com'
C
#ifdef MPI
      integer istatus(MPI_STATUS_SIZE)
      integer isum,source,dest,total,p,tag
#endif
C
      INTEGER MY_PE
      COMMON/MPICOM/MY_PE
C
C     Scalar arguments:
C
      INTEGER NPOIN,NDIM,NOFVAR,NOFVERT,NGHOST,NPNOD,NBFAC,NELEM,ROWBGN
      LOGICAL COMPRESSIBLE
C
C     Array arguments:
C
      INTEGER NODCODE(NPOIN+NGHOST),IDX(NPOIN),IBNDFAC(3,NBFAC)
      INTEGER ICELNOD(NOFVERT,*)
C
C     Local scalars:
C
      INTEGER IPOIN,NI,IFAIL,IOFF,IFRST,J,IVERT,I,ICLR,IFACE,IELEM
      INTEGER ICYCL
C
      IS, dimension(0:*) :: NodalBcs 
C
      WRITE(NOUT,1000)MY_PE
#ifdef MPI
      IF(MY_PE.EQ.0)WRITE(IWUNIT,1000)MY_PE
#endif
C
C
C     +-----------------------------------------------------+
C     Test for no-slip b.c. (IBC=6)
C     and prescribed velocity-profile (IBC=7)
C     +-----------------------------------------------------+
C
C
      CALL IINIT(NPOIN+NGHOST,0,NODCODE,1)
C
      DO 1 IFACE = 1 , NBFAC
C
          ICLR = ICOLOR(IBNDFAC(3,IFACE),1) ! colours for the mean flow equations
C
C        faces with no-slip b.c. or prescribed velocity profile
C
          IF( .NOT.(ICLR .EQ. BC_TYPE_NO_SLIP .OR. 
     &              ICLR .EQ. BC_TYPE_PROFILE) ) GOTO 1
          IELEM = IBNDFAC(1,IFACE)
          IVERT = IBNDFAC(2,IFACE)
C
          DO 2 I =  1, NOFVERT - 1
C
              IPOIN = ICELNOD( ICYCL( IVERT+I , NOFVERT ) , IELEM )
C
              NODCODE(IPOIN) = 1000
C
    2     CONTINUE
    1 CONTINUE
C
C     Create an Index Set of No-slip velocity components
C     (to be used to set b.c. for N-S computations)
C     IOFF is the # of variables stored BEFORE the velocity
C     components
C
      IF(COMPRESSIBLE)THEN
          IOFF = 2
      ELSE
          IOFF = 1
      ENDIF
C
      NI=0
      I=0
      DO 230 IPOIN = 1, NPOIN
          IF(NODCODE(IPOIN).EQ.1000)THEN
              I=I+1
              IFRST = (IPOIN-1)*NOFVAR+IOFF -  1+ ROWBGN
              DO 220 J = 1,NDIM
                  NI=NI+1
                  IDX(NI) = IFRST + J
  220         CONTINUE
          ENDIF
  230 CONTINUE
C
      CALL ISCreateGeneral(PETSC_COMM_SELF,NI,IDX,
     +PETSC_COPY_VALUES,NodalBcs(NoSlipVelocities),IFAIL)
      WRITE(NOUT,2000)NI/NDIM,CBTYPE(BC_TYPE_NO_SLIP)
C
#ifdef MPI
      CALL MPI_COMM_SIZE(PETSC_COMM_WORLD, p, ifail)
      TAG = 50
      DEST = 0
      if (MY_PE .EQ. 0) then
          total = NI
          do source = 1, p-1
              call MPI_RECV(isum, 1, MPI_INTEGER, source, tag, 
     +              PETSC_COMM_WORLD, istatus, ifail)
              total = total + isum
          enddo 
          WRITE(IWUNIT,2000)total/NDIM,CBTYPE(BC_TYPE_NO_SLIP)
      else
          call MPI_SEND(NI, 1, MPI_INTEGER, dest, 
     +          tag, PETSC_COMM_WORLD, ifail)
      endif
#endif
C
C     +-----------------------------------------------------+
C     Create an index set containing the row numbers
C     (0-based indexing) where \rhoE is stored for
C     meshpoints belonging to
C     a solid wall no-slip boundary (IBC=6)
C    (to be used to set b.c. for isothermal b.c.)
C     +-----------------------------------------------------+
C
      IF(IADIA.NE.0)THEN
C
!        CALL SETFLAG(IBNDFAC,NBFAC,NODCODE,NPOIN,NGHOST,
!    &                ICELNOD,NOFVERT,ICOLOR,BC_TYPE_NO_SLIP)
         CALL IINIT(NPOIN+NGHOST,0,NODCODE,1)
C
         DO 41 IFACE = 1 , NBFAC
C
             ICLR = ICOLOR(IBNDFAC(3,IFACE),1) ! colours for the mean flow equations
C
C        faces with no-slip b.c.
C
             IF( ICLR .NE. BC_TYPE_NO_SLIP) GOTO 41
             IELEM = IBNDFAC(1,IFACE)
             IVERT = IBNDFAC(2,IFACE)
C
             DO 42 I =  1, NOFVERT - 1
C
                 IPOIN = ICELNOD( ICYCL( IVERT+I , NOFVERT ) , IELEM )
C
                 NODCODE(IPOIN) = 1000
C
   42        CONTINUE
   41    CONTINUE
C
C     IOFF-1=+1 since the storage is: \rho,\rho E, etc.
C
         NI=0
         IOFF=2
         DO 330 IPOIN = 1, NPOIN
             IF(NODCODE(IPOIN).EQ.1000)THEN
                 NI=NI+1
                 IDX(NI) = (IPOIN-1)*NOFVAR+IOFF -  1+ ROWBGN
             ENDIF
  330    CONTINUE
C
      ELSE
          NI=0
      ENDIF
      CALL ISCreateGeneral(PETSC_COMM_SELF,NI,IDX,
     &                     PETSC_COPY_VALUES,NodalBcs(Isothermal),
     &                     IFAIL)
      WRITE(NOUT,2000)NI,"ISOTHERMAL BCS"
C
C
C     +-----------------------------------------------------+
C     Test for supersonic b.c. (IBC=1)
C     +-----------------------------------------------------+
C
!        CALL SETFLAG(IBNDFAC,NBFAC,NODCODE,NPOIN,NGHOST,
!    &                ICELNOD,NOFVERT,ICOLOR,BC_TYPE_SUPS_INLET)
C
      CALL IINIT(NPOIN+NGHOST,0,NODCODE,1)
C
      DO 3 IFACE = 1 , NBFAC
C
          ICLR = ICOLOR(IBNDFAC(3,IFACE),1) ! colours for the mean flow equations
C
C        faces with dirichlet b.c.
C
          IF( ICLR .NE. BC_TYPE_SUPS_INLET )GOTO 3
          IELEM = IBNDFAC(1,IFACE)
          IVERT = IBNDFAC(2,IFACE)
C
          DO 4 I =  1, NOFVERT - 1
C
              IPOIN = ICELNOD( ICYCL( IVERT+I , NOFVERT ) , IELEM )
C
              NODCODE(IPOIN) = BC_TYPE_SUPS_INLET
C
    4     CONTINUE
    3 CONTINUE
C
C
C     Create an index set containing the indices of the variables
C     (0-based indexing) of meshpoints belonging to
C     a supersonic (or Dirichlet) boundary
C     to be used in Euler/N-S computations
C
      NI=0
      I=0
      DO 130 IPOIN = 1, NPOIN
          IF(NODCODE(IPOIN).EQ.BC_TYPE_SUPS_INLET)THEN
              I=I+1
              NI=NI+1
              IDX(NI) = (IPOIN-1) + ROWBGN
!             IDX(NI) = (IPOIN-1)*NOFVAR + ROWBGN ! changed to comply to petsc-3.2
          ENDIF
  130 CONTINUE
      WRITE(NOUT,2000)I,CBTYPE(BC_TYPE_SUPS_INLET)
      IF(NOFVAR.EQ.1)THEN
          CALL ISCreateGeneral(PETSC_COMM_SELF,NI,IDX,PETSC_COPY_VALUES,
     +                         NodalBcs(SupersonicNodes),IFAIL)
      ELSE
          CALL ISCreateBlock(PETSC_COMM_SELF,NOFVAR,NI,IDX,
     +    PETSC_COPY_VALUES,NodalBcs(SupersonicVariables),IFAIL)
      ENDIF
C
#ifdef MPI
      call MPI_COMM_SIZE(PETSC_COMM_WORLD, p, ifail)
      tag = 50
      dest = 0
      if (MY_PE .EQ. 0) then
          total = NI
          do source = 1, p-1
              call MPI_RECV(isum, 1, MPI_INTEGER, source, tag, 
     +              PETSC_COMM_WORLD, istatus, ifail)
              total = total + isum
          enddo 
          WRITE(IWUNIT,2000)total,CBTYPE(BC_TYPE_SUPS_INLET)
      else
          call MPI_SEND(NI, 1, MPI_INTEGER, dest, 
     +          tag, PETSC_COMM_WORLD, ifail)
      endif
#endif
C
      IF( LALE .AND. ALE_LAPLACE_GRID )THEN 
C
C     +-----------------------------------------------------+
C     Test for Dirchlet boundary conditions for the Motion solver
C     all boundaries are of Dirchlet type
C     +-----------------------------------------------------+
C
         CALL IINIT(NPOIN+NGHOST+NPNOD,0,NODCODE,1)
C
         DO 17 IFACE = 1 , NBFAC
C
            IELEM = IBNDFAC(1,IFACE)
            IVERT = IBNDFAC(2,IFACE)
C
            DO 19 I =  1, NOFVERT - 1
C
                IPOIN = ICELNOD(ICYCL(IVERT+I,NOFVERT),IELEM)
C
                NODCODE(IPOIN) = NODCODE(IPOIN)+1
C
   19       CONTINUE
   17    CONTINUE
C
         NI=0
         DO 739 IPOIN = 1, NPOIN ! ONLY interior nodes should be addressed
            IF(NODCODE(IPOIN).NE.0)THEN
                IFRST = (IPOIN-1)+ ROWBGN ! 0-based indexing for Petsc
                NI=NI+1
                IDX(NI) = IFRST
            ENDIF
  739 CONTINUE
C
      CALL ISCreateGeneral(PETSC_COMM_SELF,NI,IDX,
     &                     PETSC_COPY_VALUES,NodalBcs(MotionSolverBCS),
     &                     IFAIL)
      WRITE(NOUT,2000)NI,'MotionSolverBCS '
C
      ENDIF ! LapLAcE's equation for the motion solver
C
C
      IF( LAPLACE )THEN 
C
C     +-----------------------------------------------------+
C     Test for Dirchlet boundary conditions for Poisson's eqn.
C     +-----------------------------------------------------+
C
         CALL IINIT(NPOIN+NGHOST+NPNOD,0,NODCODE,1)
C
         DO 7 IFACE = 1 , NBFAC
C
            ICLR = ICOLOR(IBNDFAC(3,IFACE),3) ! colours for the mean flow equations
C
C        faces with Dirichlet bndry conditions
C
          IF( ICLR .NE. BC_TYPE_SUPS_INLET )GOTO 7 ! pick up all boundary nodes
             IELEM = IBNDFAC(1,IFACE)
             IVERT = IBNDFAC(2,IFACE)
C
            DO 9 I =  1, NOFVERT - 1
C
C
C
                IPOIN = ICELNOD(ICYCL(IVERT+I,NOFVERT),IELEM)
C
                NODCODE(IPOIN) = NODCODE(IPOIN)+1
C
    9       CONTINUE
    7    CONTINUE
C
         NI=0
         IOFF = NOFVAR ! the potential is stored last
         DO 734 IPOIN = 1, NPOIN
            IF(NODCODE(IPOIN).NE.0)THEN
                IFRST = (IPOIN-1)*NOFVAR+IOFF -  1+ ROWBGN
                NI=NI+1
                IDX(NI) = IFRST
            ENDIF
  734 CONTINUE
C
      CALL ISCreateGeneral(PETSC_COMM_SELF,NI,IDX,
     &                     PETSC_COPY_VALUES,NodalBcs(Dirichlet4Poisson)
     &                    ,IFAIL)
      WRITE(NOUT,2000)NI,'Dirichlet 4 Poisson'
C
      ENDIF ! LapLAcE
C
C     +-----------------------------------------------------+
C     +-----------------------------------------------------+
C
C
C     +-----------------------------------------------------+
C     Test for hanging ndoes
C     +-----------------------------------------------------+
C
      CALL IINIT(NPOIN+NGHOST,-1,NODCODE,1)
C
      DO 31 IELEM = 1 , NELEM
C
          DO 32 IVERT =  1, NOFVERT
C
              IPOIN = ICELNOD( IVERT , IELEM )
              NODCODE(IPOIN) = 0
C
   32     CONTINUE
   31 CONTINUE
C
C     Create an Index Set of hanging nodes
C
      NI=0
      DO 732 IPOIN = 1, NPOIN
          IF(NODCODE(IPOIN).EQ.-1)THEN
              IFRST = (IPOIN-1)*NOFVAR-  1+ ROWBGN
              DO I= 1,NOFVAR
                 NI=NI+1
                 IDX(NI) = IFRST + I
              ENDDO
          ENDIF
  732 CONTINUE
C
      CALL ISCreateGeneral(PETSC_COMM_SELF,NI,IDX,
     &                     PETSC_COPY_VALUES,NodalBcs(HangingNodes),
     &                     IFAIL)
      WRITE(NOUT,2000)NI/NOFVAR,"HANGING NODES"
C
#ifdef MPI
      call MPI_COMM_SIZE(PETSC_COMM_WORLD, p, ifail)
      tag = 50
      dest = 0
      if (MY_PE .EQ. 0) then
          total = NI
          do source = 1, p-1
              call MPI_RECV(isum, 1, MPI_INTEGER, source, tag, 
     +              PETSC_COMM_WORLD, istatus, ifail)
              total = total + isum
          enddo 
          WRITE(IWUNIT,2000)total/NOFVAR,"HANGING NODES"
      else
          call MPI_SEND(NI, 1, MPI_INTEGER, dest, 
     +          tag, PETSC_COMM_WORLD, ifail)
      endif
#endif
C
      RETURN
 1000 FORMAT(/,5X,'NODAL BOUNDARY CONDITIONS ON PE # ',I4,/ 5X,38("="))
 
 2000 FORMAT(I7,' NODES with B.C. ',A)
      END
C
      SUBROUTINE SETFLAG(IBNDFAC,NBFAC,NODCODE,NP,
     &ICELNOD,NOFVERT,ICOLOR,BC_TYPE)
C
C
      INTEGER IBNDFAC(3,NBFAC),NODCODE(*),ICELNOD(NOFVERT,*),ICOLOR(*)
      INTEGER NBFAC,NP,BC_TYPE
      INTEGER IFACE,ICLR,IELEM,IVERT,I
      INTEGER ICYCL
C
      CALL IINIT(NP,0,NODCODE,1)
C
      DO 1 IFACE = 1 , NBFAC
C
          ICLR = ICOLOR(IBNDFAC(3,IFACE))
C
C        faces with no-slip b.c. or prescribed velocity profile
C
          IF( .NOT.(ICLR .EQ. BC_TYPE)) GOTO 1
C
          IELEM = IBNDFAC(1,IFACE)
          IVERT = IBNDFAC(2,IFACE)
C
          DO 2 I =  1, NOFVERT - 1
C
              IPOIN = ICELNOD( ICYCL( IVERT+I , NOFVERT ) , IELEM )
C
              NODCODE(IPOIN) = BC_TYPE
C
    2     CONTINUE
    1 CONTINUE
      RETURN
      END
