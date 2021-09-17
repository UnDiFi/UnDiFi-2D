      SUBROUTINE ISET_T(NodalBcs, IDX, IBNDFAC, NBFAC, ICELNOD, NOFVERT,
     &NODCODE, NDIM, NOFVAR, NPOIN, NGHOST, ROWBGN, COMPRESSIBLE )
C
C     $Id: iset_t.F,v 1.16 2020/04/23 09:46:19 abonfi Exp $
C
C     This routine creates PETSc Index Sets addressing
C          boundary meshpoints with various types of
C          boundary conditions pertaining to the turbulence modeling.
C     In the parallel case it is important that
C          the routine is called with NPOIN=# of processor
C          owned nodes, since the index sets are local
C          and shoud not address ghost nodes 
C
#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscis.h"
      use petscvec 
      use petscis 
C
      IMPLICIT NONE
C
      INCLUDE 'bnd.h'
      INCLUDE 'iset.h'
C
      INCLUDE 'bnd.com'
      INCLUDE 'flags.com'
      INCLUDE 'io.com'
      INCLUDE 'visco.com'
C
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
      INTEGER NPOIN,NDIM,NOFVAR,NOFVERT,NGHOST,NBFAC,ROWBGN
      LOGICAL COMPRESSIBLE
C
C     Array arguments:
C
      INTEGER NODCODE(NPOIN+NGHOST),IDX(NPOIN),IBNDFAC(3,NBFAC)
      INTEGER ICELNOD(NOFVERT,*)
C
C     Local scalars:
C
      INTEGER IPOIN,NI,IFAIL,IVERT,I,ICLR,IFACE,IELEM,
     &        NOFEQN,ICLRT
      INTEGER ICYCL
C
      IS, dimension(0:*) :: NodalBcs
C
      WRITE(NOUT,FMT=1000) MY_PE
#ifdef MPI
      IF(MY_PE.EQ.0)WRITE(IWUNIT,FMT=1000) MY_PE
#endif
C
C     set the number of mean flow eqns.
C
      IF( COMPRESSIBLE )THEN
          NOFEQN = NDIM + 2
      ELSE
          NOFEQN = NDIM + 1 
      ENDIF 
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
         ICLR  = ICOLOR(IBNDFAC(3,IFACE),1) ! colours for the mean flow equations
         ICLRT = ICOLOR(IBNDFAC(3,IFACE),2) ! colours for the turbulence transport equations
C
C        faces with no-slip b.c.
C        we do not consider faces with prescribed velocity profile
C        this should be handled with the b.c. freestream_turbulence
C        however we use -turbulence_colours to set 
C        the turbulent variable to zero on boundaries were
C        a bc other than no-slip is used for the mean vars
C
         IF( .NOT.(ICLR  .EQ. BC_TYPE_NO_SLIP
     &       .OR.  ICLRT .EQ. BC_TYPE_NO_SLIP) )   GOTO 1
         IELEM = IBNDFAC(1,IFACE)
         IVERT = IBNDFAC(2,IFACE)
C
         DO 2 I = 1 , NOFVERT - 1
C
            IPOIN = ICELNOD( ICYCL( IVERT+I , NOFVERT ) , IELEM )
C
            NODCODE(IPOIN) = 1000
C
    2    CONTINUE
    1 CONTINUE ! loop over faces
C
C     +-----------------------------------------------------+
C     Create an index set containing the node numbers
C     (0-based indexing) of meshpoints belonging to
C     a solid wall no-slip boundary (IBC=6) 
C     (to be used to set b.c. for turb model in N-S computations)
C     +-----------------------------------------------------+
C
C
      NI=0
      DO 200 IPOIN = 1, NPOIN
         IF(NODCODE(IPOIN).EQ.1000)THEN 
             NI=NI+1
             IF( SEGREGATED ) THEN
                 IDX(NI) = IPOIN-1+ROWBGN
             ELSE
                 IDX(NI) = (IPOIN-1)*NOFVAR+NOFEQN+ROWBGN
             ENDIF
         ENDIF 
  200 CONTINUE
C
      WRITE(NOUT,2000)NI,CBTYPE(BC_TYPE_NO_SLIP)
      CALL ISCreateGeneral(PETSC_COMM_SELF,NI,IDX,
     +                     PETSC_COPY_VALUES,NodalBcs(NoSlipNodes),
     &                     IFAIL)
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
          WRITE(IWUNIT,2000)total,CBTYPE(BC_TYPE_NO_SLIP)
      else
          call MPI_SEND(NI, 1, MPI_INTEGER, dest, 
     +          tag, PETSC_COMM_WORLD, ifail)
      endif
#endif
C
      CALL IINIT(NPOIN+NGHOST,0,NODCODE,1) 
C
      DO 3 IFACE = 1 , NBFAC 
C
         ICLR = ICOLOR(IBNDFAC(3,IFACE),2) ! col 2 means colour for the turbulent quantity
C
C        ICLR .EQ. 1 means prescribed turbulent quantity
C
         IF( ICLR .NE. 1 )GOTO 3
         IELEM = IBNDFAC(1,IFACE)
         IVERT = IBNDFAC(2,IFACE)
C
         DO 4 I = 1 , NOFVERT - 1
C
            IPOIN = ICELNOD( ICYCL( IVERT+I , NOFVERT ) , IELEM )
C
            NODCODE(IPOIN) = 1
C
    4 CONTINUE
C
    3 CONTINUE
C
C     Create an index set containing the node numbers
C     (0-based indexing) of meshpoints belonging to
C     a boundary where the turbulent variable should be held fixed.
C
      NI=0
      DO 100 IPOIN = 1, NPOIN
         IF(NODCODE(IPOIN).EQ.1)THEN 
             NI=NI+1
C
C     0-based indexing
C
             IF( SEGREGATED )THEN
                 IDX(NI) = (IPOIN-1)+ROWBGN
             ELSE ! coupled
                 IDX(NI) = (IPOIN-1)*NOFVAR+NOFEQN+ROWBGN
             ENDIF 
         ENDIF 
  100 CONTINUE
C
      CALL ISCreateGeneral(PETSC_COMM_SELF,NI,IDX,PETSC_COPY_VALUES,
     +                     NodalBcs(FreestreamTurbulence),IFAIL)
      WRITE(NOUT,2000)NI,'FREESTREAM TURBULENCE'
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
          WRITE(IWUNIT,2000)total,'FREESTREAM TURBULENCE'
      else
          call MPI_SEND(NI, 1, MPI_INTEGER, dest, 
     +          tag, PETSC_COMM_WORLD, ifail)
      endif
#endif
C
      RETURN
 1000 FORMAT(/,5X,'NODAL BOUNDARY CONDITIONS (TURBULENCE) ON PE # ',I4,/
     &       5X,47("="))
 2000 FORMAT(I7,' NODES with B.C. ',A)
      END 
C
