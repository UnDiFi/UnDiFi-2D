      SUBROUTINE INIT( NodalBcs, IBNDFAC, NBFAC, ICELNOD, NOFVERT, X,
     &                 VCORG, ROWBGN, NDIM, NOFVAR, NPOIN, COMPRESSIBLE,
     &                 RESTART )
C
C     $Id: init.F,v 1.21 2020/04/23 09:46:19 abonfi Exp $
C
CCCCC#define COUETTE
C
#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscis.h"
      use petscvec
      use petscis

      IMPLICIT NONE
C
C    the BIG warning is: GHOST nodes will NOT be
C    initialised, Index Sets ONLY address interior nodes
C
C    NOTE that this routine must be called AFTER
C    the flow variables have been initialized,
C    so that the velocity components can be set to 0.0 
C    on viscous walls
C
C    Some of the Index sets use global numbering to address
C    specific rows of Vecs and Mats. This is because
C    use of local numbering would require a NON-blocked
C    mapping for the matrices.
C    
C
C     Parameters ..
C
      INCLUDE 'constants.h'
      INCLUDE 'paramt.h'
      INCLUDE 'bnd.h'
      INCLUDE 'iset.h'
C
C     Commons ..
C
      INCLUDE 'bnd.com'
      INCLUDE 'io.com'
      INCLUDE 'visco.com'
      INCLUDE 'flags.com'
      INCLUDE 'stream.com'
      INCLUDE 'pfcgas.com'
C
      INTEGER MY_PE
      COMMON/MPICOM/MY_PE
C
      PetscOffset IDX_I
      INTEGER IDX_V(1)
C
C     Scalar Arguments ..
C
      INTEGER NDIM,NOFVAR,NPOIN,ROWBGN,NOFVERT,NBFAC
      LOGICAL COMPRESSIBLE,RESTART
C
C     Array Arguments ..
C
      DOUBLE PRECISION X(*)
      DOUBLE PRECISION VCORG(NDIM,NPOIN)
      INTEGER IBNDFAC(3,NBFAC),ICELNOD(NOFVERT,*)
C
C     Local Scalars ..
C
      INTEGER IPOIN,NI,IFAIL,IOFF,IFACE,IRE,IBC,IELEM,I,J,K,IVERT,IADD
      DOUBLE PRECISION ALPHA
      INTEGER NERR,IOPT
      CHARACTER*72 ERRMSG
      PARAMETER(NERR = 5, IOPT = 1)
C
      IS, dimension(0:*) :: NodalBcs 
C
C     External Functions ..
C
      INTEGER  ICYCL
      EXTERNAL ICYCL
C
C     Intrinsic Functions ..
C
C
C     Executable Statements ..
C
C     no-slip nodes are identyfied by looping over the faces
C     -NOT- as nodes with a prescribed velocity profile
C     as these also includes nodes where a velocity profile
C     is specified
C
      IF(COMPRESSIBLE)THEN
          IOFF =  2
      ELSE
          IOFF =  1
      ENDIF
C
C     IOFF is the offset between the location where 
C     the first variable is stored for a given node and the one where
C     the first velocity component is
C
      DO 3 IFACE = 1, NBFAC
         IBC=IBNDFAC(3,IFACE)
         IF(ICOLOR(IBC,1).NE.BC_TYPE_NO_SLIP)GOTO 3 ! colours for the mean flow equations
         IELEM=IBNDFAC(1,IFACE)
         IVERT=IBNDFAC(2,IFACE)
         DO 4 J = 1,NOFVERT-1
            IPOIN = ICELNOD(ICYCL(IVERT+J,NOFVERT),IELEM)
            IADD = (IPOIN-1)*NOFVAR + IOFF
C
C         set u=v=w=0. ...
C
            DO 4 K = 1, NDIM
                X(IADD+K) = ZERO
    4    CONTINUE
    3 CONTINUE
C
#ifdef COUETTE
      write(6,*)'Ad hoc for Couette!!!'
C
      DO 13 IFACE = 1, NBFAC
         IBC=IBNDFAC(3,IFACE)
         IF(IBC.NE.3)GOTO 13
         IELEM=IBNDFAC(1,IFACE)
         IVERT=IBNDFAC(2,IFACE)
         DO 14 J = 1,NOFVERT-1
            IPOIN = ICELNOD(ICYCL(IVERT+J,NOFVERT),IELEM)
            IADD = (IPOIN-1)*NOFVAR + IOFF
C
C         set u=v=w=0. ...
C
                X(IADD+1) = ONE
   14 CONTINUE
   13 CONTINUE
#endif
C
      IF(IADIA.EQ.0)RETURN
C
C     Isothermal bcs on all solid surfaces
C
      IF(RESTART)THEN
         IF(MY_PE.EQ.0)WRITE(IWUNIT,FMT=*)'Isothermal walls are not any 
     &longer initialized to wall temperature'
         RETURN
      ENDIF
      
C
      CALL ISGetSize(NodalBcs(Isothermal),NI,IFAIL)
      IF(NI.NE.0)THEN
         CALL ISGetIndices(NodalBcs(Isothermal),IDX_V,IDX_I,IFAIL)
         ALPHA = RSTAR * TWALL /GM1
         DO 120 I = 1,NI
C
C        set 1-based indexing
C
             IRE = IDX_V(IDX_I+I)+1-ROWBGN
C
C where total enthalpy is stored
C
C
C     .. Isothermal wall .. 
C
            CALL PARM_TO_CONS(X(IRE-1),NDIM,NOFVAR,1,.FALSE.,IFAIL)
            X(IRE) = ALPHA*X(IRE-1)
            CALL CONS_TO_PARM(X(IRE-1),NDIM,NOFVAR,1,.FALSE.,IFAIL)
            IF(IFAIL.NE.0)THEN
               WRITE(ERRMSG(1:64),FMT=105)'INIT______',MY_PE,"(0)"
               CALL SETERR(ERRMSG(1:64),64,NERR,IOPT)
            ENDIF   ! ifail
  120    CONTINUE
         CALL ISRestoreIndices(NodalBcs(Isothermal),IDX_V,IDX_I,IFAIL)
      ENDIF
C
      RETURN
  105 FORMAT('CONS_TO_PARM failed when called from ',A10,' on PE # ',I4,
     &1X,A3)
      END
