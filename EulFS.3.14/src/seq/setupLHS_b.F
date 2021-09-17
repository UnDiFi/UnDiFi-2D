      SUBROUTINE SetupLHS_b(ICELNOD,ICELFAC,VFACNOR,CORG,XYZDOT,VOL,
     2                      CDIST,
     3                      NELEM,NPOIN,NGHOST,NPNOD,NDIM,NOFVERT,
     4                      C)
C
C     $Id: setupLHS_b.F,v 1.6 2020/04/23 09:37:34 abonfi Exp $
C
C
#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscmat.h"
      use petscvec
      use petscmat
C
      IMPLICIT NONE
C
      INCLUDE 'paramt.h'
      INCLUDE 'bnd.h'
      INCLUDE 'implicit.h'
C
      Mat C
C
      INCLUDE 'constants.h'
      INCLUDE 'three.com'
      INCLUDE 'flags.com'
      INCLUDE 'time.com'
C
      INTEGER NDIM,NOFVERT,NELEM,NPOIN,NGHOST,NPNOD,NOFVAR
C
C     On entry:
C     --------
C
C     NDIM    dimension of the space (2 or 3)
C     NOFVERT number of vertices per element (=NDIM+1, since
C             only triangles or tetrahedra are allowed)
C     NELEM   no. of processor owned elements (triangles/tetrahedra);
C             global number of elements in the uni-processor case
C
      INTEGER ICELNOD(NOFVERT,NELEM),ICELFAC(NOFVERT,NELEM)
C
C     ICELNOD(1:NOFVERT,1:NELEM)
C            Cell to Node pointer : ICELNOD(i,ielem) gives the
C            global node number of the i-th vertex of the ielem-th cell
C
      DOUBLE PRECISION VFACNOR(NDIM,*),XYZDOT(NDIM,*),CORG(NDIM,*),
     2                 VOL(NELEM),CDIST(NELEM)

      DOUBLE PRECISION VCZ(MAXNOFVAR*MAXNOFVERT*MAXTIMLEVS),
     &VCB(3*MAXNOFVERT)
      DOUBLE PRECISION VCN(3*MAXNOFVERT),STIFEL(MAX_NOFVERT_SQR)
C
C     FACNOR(1:NDIM,1:NFACE)  cartesian components
C                             of the NFACE faces/edges
C
      INTEGER JVAR,IVAR,IVERT,JVERT,IELEM,NP,IADDR
      INTEGER IFAIL,N4
      DOUBLE PRECISION EPS,VOLUME(MAXTIMLEVS+1)
C
C     ICN stores the vertices of the current element (0-based indexing)
C
      INTEGER ICN(MAXNOFVERT*MAXTIMLEVS)
      LOGICAL VERBOSE
      PARAMETER( VERBOSE = .FALSE. )
C
C
C
C     Some initializations ....
C
      N4 = NOFVERT*NOFVERT
      NP = NPOIN + NGHOST + NPNOD 
      IF( VERBOSE) THEN
          CALL I4Mat_Print('General',' ',Nofvert,nelem,icelnod,
     +                Nofvert, 'cell connectivity ',IFAIL)
          CALL R8Mat_Print('General',' ',NDIM,NP,XYZDOT,
     +                NDIM, 'Grid velocity ',IFAIL)
      ENDIF
C
      DO 2000 IELEM = 1,NELEM
C
C     The local element stiffness matrix is initialized to 0.d0
C
          CALL DINIT(N4,ZERO,STIFEL,1)
C
C
C
!                                                           
          CALL CELPTR(IELEM, NELEM, ICELNOD, ICELFAC, VOL, CORG,
!                                                          ^
!                                                          |
!                                                          |
!          we pass the grid velocity because we have there are 3 time levels in the un-steady case
!          otherwise there will be an access out-of-bounds in celptr
     +                VFACNOR, XYZDOT, NDIM, NOFVERT, NDIM, NP, ICN, 
     3                VCZ, VCN, VCB, VOLUME)
C
C**********************************************
C**********************************************
C
C
C
      IF( VERBOSE )THEN
          CALL R8Mat_Print('General',' ',NDIM,NofvERT,VCB,
     +                NDIM, 'VCB ',IFAIL)
          CALL R8Mat_Print('General',' ',NDIM,NofvERT*NTIMLEVS,VCZ,
     +                NDIM,'VCZ ',IFAIL)
          CALL R8Mat_Print('General',' ',NDIM,NofvERT,VCN,
     +                NDIM, 'VCN ',IFAIL)
          write(6,*)(icn(jvert),jvert=1,nofvert) 
      ENDIF
C
C
C**********************************************
C
C
C**********************************************
C
          CALL LAPLACE4MOTION(IELEM,VCN,NDIM,NOFVERT,
     &                        CDIST(IELEM),STIFEL,VOLUME(1))
      IF( VERBOSE )THEN
          CALL R8Mat_Print('General',' ',Nofvert,NofvERT,STIFEL,
     +                Nofvert,
     &                'Stiffness matrix velocity (2)',IFAIL)
      pause
      ENDIF
C
C     Insert elemental matrices into the global one
C
          DO 35 IVERT = 1,NOFVERT
             DO 35 JVERT = 1,NOFVERT
                IADDR = (JVERT-1)*NOFVERT + IVERT
#ifdef MPI
                CALL MatSetValuesLocal(C,1,ICN(IVERT),1, ICN
     +                (JVERT),STIFEL(IADDR), ADD_VALUES,IFAIL)
 
#else
                CALL MatSetValues(C,1,ICN(IVERT),1,ICN(JVERT),
     +                STIFEL(IADDR),ADD_VALUES,IFAIL)
#endif
   35         CONTINUE
C
 2000 CONTINUE ! end loop over elements
caldo
      CALL MatAssemblyBegin(C,MAT_FINAL_ASSEMBLY,IFAIL)
      CALL MatAssemblyEnd(C,MAT_FINAL_ASSEMBLY,IFAIL)
caldo
      RETURN
      END
