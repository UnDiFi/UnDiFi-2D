      SUBROUTINE ISET(IDX, IBNDFAC, NBFAC, ICELNOD, NOFVERT,
     &NODCODE, NDIM, NOFVAR, NPOIN, NGHOST, ROWBGN, COMPRESSIBLE )
C
C     $Id: iset.F,v 1.15 2001/01/20 13:50:34 aldo Exp abonfi $
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
      IMPLICIT NONE
C
#include "include/finclude/petsc.h"
#include "include/finclude/petscvec.h"
#include "include/finclude/petscis.h"
C
      INCLUDE 'bnd.h'
      INCLUDE 'bnd'
      INCLUDE 'io.com'
      INCLUDE 'visco.com'
C
#include "iset.com"
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
      INTEGER IPOIN,NI,IFAIL,IOFF,IFRST,J,IVERT,I,ICLR,IFACE,IELEM
      INTEGER ICYCL
C
      WRITE(NOUT,1000)MY_PE
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
          ICLR = ICOLOR(IBNDFAC(3,IFACE))
C
C        faces with no-slip b.c. or prescribed velocity profile
C
          IF( .NOT.(ICLR .EQ. BC_NO_SLIP .OR. ICLR .EQ. BC_PROFILE) )
     +    GOTO 1
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
      CALL ISCreateGeneral(PETSC_COMM_SELF,NI,IDX,NoSlipVelocities,
     +IFAIL)
      WRITE(NOUT,2000)NI/NDIM,CBTYPE(BC_NO_SLIP)
C
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
          DO 300 I= 1, NI
C
C     +1 since the storage is: \rho,\rho E, etc.
C
              IDX(I) = IDX(I)*NOFVAR+1 + ROWBGN
  300     CONTINUE
      ELSE
          NI=0
      ENDIF
      CALL ISCreateGeneral(PETSC_COMM_SELF,NI,IDX,Isothermal, IFAIL)
 
C
C
C     +-----------------------------------------------------+
C     Test for supersonic b.c. (IBC=1)
C     +-----------------------------------------------------+
C
C
      CALL IINIT(NPOIN+NGHOST,0,NODCODE,1)
C
      DO 3 IFACE = 1 , NBFAC
C
          ICLR = ICOLOR(IBNDFAC(3,IFACE))
C
C        faces with dirichlet b.c.
C
          IF( ICLR .NE. BC_SUPS_INLET )GOTO 3
          IELEM = IBNDFAC(1,IFACE)
          IVERT = IBNDFAC(2,IFACE)
C
          DO 4 I =  1, NOFVERT - 1
C
              IPOIN = ICELNOD( ICYCL( IVERT+I , NOFVERT ) , IELEM )
C
              NODCODE(IPOIN) = BC_SUPS_INLET
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
          IF(NODCODE(IPOIN).EQ.BC_SUPS_INLET)THEN
              I=I+1
              NI=NI+1
              IDX(NI) = (IPOIN-1)*NOFVAR + ROWBGN
          ENDIF
  130 CONTINUE
      WRITE(NOUT,2000)I,CBTYPE(BC_SUPS_INLET)
      IF(NOFVAR.EQ.1)THEN
          CALL ISCreateGeneral(PETSC_COMM_SELF,NI,IDX, SupersonicNodes,
     +    IFAIL)
      ELSE
          CALL ISCreateBlock(PETSC_COMM_SELF,NOFVAR,NI,IDX,
     +    SupersonicVariables,IFAIL)
      ENDIF
C
C     +-----------------------------------------------------+
C     Test for constant pressure b.c. (IBC=2)
C     +-----------------------------------------------------+
C
      CALL IINIT(NPOIN+NGHOST,0,NODCODE,1)
C
      DO 5 IFACE = 1 , NBFAC
C
          ICLR = ICOLOR(IBNDFAC(3,IFACE))
C
C        faces with dirichlet b.c.
C
          IF( ICLR .NE. BC_CNST_PRES )GOTO 5
          IELEM = IBNDFAC(1,IFACE)
          IVERT = IBNDFAC(2,IFACE)
C
          DO 6 I =  1, NOFVERT - 1
C
              IPOIN = ICELNOD( ICYCL( IVERT+I , NOFVERT ) , IELEM )
C
              NODCODE(IPOIN) = 3000
C
    6     CONTINUE
    5 CONTINUE
C
C     for compressible flows the Index Set addresses the energy (2nd) eqn.
C     while for incompressible it addresses the mass (1st) eqn.
C     this is the reason for using IOFF
C
      IF(COMPRESSIBLE)THEN
          IOFF=1
      ELSE
          IOFF=0
      ENDIF
      NI=0
      DO 100 IPOIN = 1, NPOIN
          IF(NODCODE(IPOIN).EQ.3000)THEN
              NI=NI+1
C
C     0-based indexing
C
              IDX(NI) = NOFVAR*(IPOIN-1)+IOFF + ROWBGN
          ENDIF
  100 CONTINUE
      CALL ISCreateGeneral(PETSC_COMM_SELF,NI,IDX,CnstPressure, IFAIL)
 
      WRITE(NOUT,2000)NI,CBTYPE(BC_CNST_PRES)
C
C     +-----------------------------------------------------+
C     Test for inflow b.c. (IBC=8) for an internal flow
C     +-----------------------------------------------------+
C
      CALL IINIT(NPOIN+NGHOST,0,NODCODE,1)
C
      DO 7 IFACE = 1 , NBFAC
C
          ICLR = ICOLOR(IBNDFAC(3,IFACE))
C
C        faces with dirichlet b.c.
C
          IF( ICLR .NE. BC_INTERNAL_INFLOW )GOTO 7
          IELEM = IBNDFAC(1,IFACE)
          IVERT = IBNDFAC(2,IFACE)
C
          DO 9 I =  1, NOFVERT - 1
C
              IPOIN = ICELNOD( ICYCL( IVERT+I , NOFVERT ) , IELEM )
C
              NODCODE(IPOIN) = BC_INTERNAL_INFLOW
C
    9     CONTINUE
    7 CONTINUE
      NI=0
      DO 135 IPOIN = 1, NPOIN
          IF(NODCODE(IPOIN).EQ.BC_INTERNAL_INFLOW)THEN
              IF(COMPRESSIBLE)THEN
C
C     Compressible case
C     0-based indexing, addresses the equations for:
C     \rho, \rhoE, \rho v, \rho w
C     (should be made more general concerning the flow
C      directions to be imposed)
C
                  NI=NI+1
                  IDX(NI) = NOFVAR*(IPOIN-1)+ ROWBGN
                  NI=NI+1
                  IDX(NI) = NOFVAR*(IPOIN-1)+ ROWBGN +1
                  NI=NI+1
                  IDX(NI) = NOFVAR*(IPOIN-1)+ ROWBGN +3
                  IF(NDIM.EQ.3)THEN
                      NI=NI+1
                      IDX(NI) = NOFVAR*(IPOIN-1)+ ROWBGN +4
                  ENDIF
C
              ELSE
C
C     Incompressible case
C     0-based indexing, addresses the equations for:
C     p, v, w
C     (should be made more general)
C
                  NI=NI+1
                  IDX(NI) = NOFVAR*(IPOIN-1)+ ROWBGN
                  NI=NI+1
                  IDX(NI) = NOFVAR*(IPOIN-1)+ ROWBGN +2
                  IF(NDIM.EQ.3)THEN
                      NI=NI+1
                      IDX(NI) = NOFVAR*(IPOIN-1)+ ROWBGN +3
                  ENDIF
              ENDIF
          ENDIF
  135 CONTINUE
      CALL ISCreateGeneral(PETSC_COMM_SELF,NI,IDX,Internal_inflow,
     +IFAIL)
      IF(COMPRESSIBLE)THEN
          IOFF=NDIM+1
      ELSE
          IOFF=NDIM
      ENDIF
      WRITE(NOUT,2000)NI/IOFF,CBTYPE(BC_INTERNAL_INFLOW)
#if 1
C
C     +-----------------------------------------------------+
C     Test for inflow b.c. (IBC=8) for an internal flow
C     +-----------------------------------------------------+
C
      CALL IINIT(NPOIN+NGHOST,0,NODCODE,1)
C
      DO 27 IFACE = 1 , NBFAC
C
          ICLR = ICOLOR(IBNDFAC(3,IFACE))
C
C        faces with dirichlet b.c.
C
          IF( ICLR .NE. BC_INTERNAL_INFLOW )GOTO 27
          IELEM = IBNDFAC(1,IFACE)
          IVERT = IBNDFAC(2,IFACE)
C
          DO 29 I =  1, NOFVERT - 1
C
              IPOIN = ICELNOD( ICYCL( IVERT+I , NOFVERT ) , IELEM )
C
              NODCODE(IPOIN) = BC_INTERNAL_INFLOW
C
   29     CONTINUE
   27 CONTINUE
      NI=0
      DO 145 IPOIN = 1, NPOIN
          IF(NODCODE(IPOIN).EQ.BC_INTERNAL_INFLOW)THEN
              IF(COMPRESSIBLE)THEN
C
C     Compressible case
C     0-based indexing, addresses the equations for:
C     \rho, \rhoE, \rho v, \rho w
C     (should be made more general concerning the flow
C      directions to be imposed)
C
                  NI=NI+1
                  IDX(NI) = NOFVAR*(IPOIN-1)+ ROWBGN +2
C
              ELSE
C
C     Incompressible case
C     0-based indexing, addresses the equations for:
C     p, v, w
C     (should be made more general)
C
                  NI=NI+1
                  IDX(NI) = NOFVAR*(IPOIN-1)+ ROWBGN +1
              ENDIF
          ENDIF
  145 CONTINUE
      CALL ISCreateGeneral(PETSC_COMM_SELF,NI,IDX,IsDebug,
     +IFAIL)
      WRITE(NOUT,2000)NI,"debug"
C
#endif
C
      RETURN
 1000 FORMAT(/,5X,'NODAL BOUNDARY CONDITIONS ON PE # ',I4,/ 5X,38("="))
 
 2000 FORMAT(I7,' NODES with B.C. ',A)
      END
