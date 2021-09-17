      SUBROUTINE CHKFLX2(ICELNOD,ICELFAC,
     1                   VFACNOR,XYZDOT,VOL,ZROE,
     3                   NELEM,NPOIN,NGHOST,NPNOD,NDIM,NOFVERT,
     4                   NOFVAR,NTURB)
C
      IMPLICIT NONE
C
C     $Id: chkflx2.f,v 1.5 2020/02/18 08:47:40 abonfi Exp $
C
C     Here we verify the telescoping property of the residuals, i.e.
C     we sum up the flux balances of all cells and verify that these
C     equal the integral through all boundaries, which is computed
C     in subroutine chkflx
C
      INCLUDE 'paramt.h'
      INCLUDE 'constants.h'
      INCLUDE 'flags.com'
      INCLUDE 'nloc.com'
      INCLUDE 'time.com'
      INCLUDE 'three.com'
C
C
      INTEGER NPOIN,NGHOST,NPNOD,NDIM,NOFVERT,NOFVAR,NELEM,NTURB
      INTEGER ICELNOD(NOFVERT,NELEM),ICELFAC(NOFVERT,NELEM)
      DOUBLE PRECISION VFACNOR(NDIM,*),VOL(NELEM), ZROE(NOFVAR,*),
     2                 XYZDOT(NDIM,*)
C
C     On entry:
C     --------
C
C     NDIM    dimension of the space (2 or 3)
C     NOFVERT number of vertices per element (=NDIM+1, since
C             only triangles or tetrahedra are allowed)
C     NOFVAR  number of variables (degrees of freedom)
C             in each meshpoint
C     NELEM   no. of processor owned elements (triangles/tetrahedra);
C             global number of elements in the uni-processor case
C
C
C     ICELNOD(1:NOFVERT,1:NELEM)
C            Cell to Node pointer : ICELNOD(i,ielem) gives the
C            global node number of the i-th vertex of the ielem-th cell
C
C
C     FACNOR(1:NDIM,1:NFACE)  cartesian components
C                             of the NFACE faces/edges
C
      INTEGER MY_PE
      COMMON/MPICOM/MY_PE
C
C
      INTEGER IELEM,NP,J
      INTEGER NERR,IOPT
      PARAMETER(NERR = 5, IOPT = 1)
C
C     ICN stores the vertices of the current element (0-based indexing)
C
C     ..
C     .. Local Arrays ..
      INTEGER ICN(MAXNOFVERT)
      DOUBLE PRECISION VCN(3*MAXNOFVERT),VOLUME(4),WKSP(MAXNOFVAR)
      DOUBLE PRECISION VCZ(MAXNOFVAR,MAXNOFVERT),VCB(3*MAXNOFVERT)
C
C
C     Some initializations ....
C
C
      NP = NPOIN + NGHOST + NPNOD
C
C     set work array equal to zero
C
      CALL DINIT(NOFVAR,ZERO,WKSP,1)
C
      DO 2000 IELEM = 1,NELEM
C
C     The element stiffness matrix is initialized to 0.d0
C
          CALL CELPTR(IELEM, NELEM, ICELNOD, ICELFAC, VOL, ZROE,
     +                VFACNOR, XYZDOT, NDIM, NOFVERT, NOFVAR, NP, ICN,
     3                VCZ, VCN, VCB, VOLUME)
C
          CALL LINEARIZE(IELEM,LALE,VCN,VCB,NDIM,NOFVERT,
     +               VCZ,NOFVAR,VOLUME(1))
C
          CALL CHECK(IELEM,NDIM,NOFVAR)
C
          CALL DSCAL(NOFVAR,VOLUME(1),DivFlux,1)
          CALL DAXPY(NOFVAR,ONE,DivFlux,1,WKSP,1)
C
C
 2000 CONTINUE ! end loop over elements
      WRITE(6,*)
      WRITE(6,335)(WKSP(J),J=1,NOFVAR)
      RETURN
  335 FORMAT(10X,'CHECKING THE TELESCOPING PROPERTY',/,21X,
     &'OF THE INVISCID FLUXES',6(1X,E12.6))
      END
