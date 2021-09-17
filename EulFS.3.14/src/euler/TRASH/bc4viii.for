      SUBROUTINE BC5VIII(IELEM,VCN,VCZ,KPOS,NODRES,TSTEP,NDIM,NOFVERT,
     +                   NOFVAR)
C
      IMPLICIT NONE
C
C     Purpose:
C     -------
C     compute inviscid wall boundary conditions for INcompressible
C     flows (ghost cell approach)
C
      include 'paramt.h'
      include 'constants'
      include 'bnd.h'
      include 'stream.com'
C
C     .. Scalar Arguments ..
      INTEGER IELEM,NDIM,NOFVAR,NOFVERT
C
C     IELEM  current boundary element
C     NDIM   dimension of the space (2 or 3)
C     NOFVAR number of variables (degrees of freedom)
C            in each meshpoint
C     NOFVAR number of variables (degrees of freedom)
C            in each meshpoint
C
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION KPOS(NOFVAR,NOFVAR,NOFVERT-1),
     +                 NODRES(NOFVAR,NOFVERT-1),TSTEP(NOFVAR,NOFVERT-1),
     +                 VCN(NDIM,NOFVERT),VCZ(NOFVAR,NOFVERT)
C
C     On entry:
C     --------
C     VCN(1:NDIM,1:NOFVERT) cartesian components of the normals 
C                           to the element sides/faces
C                           the normal to the boundary face need
C                           to be stored in VCN(1:NDIM,NOFVERT)
C     VCZ(1:NOFVAR,1:NOFVERT) dependent variables in the vertices
C                           of the current (IELEM) element
C                           the freestream values need
C                           to be stored in VCZ(1:NOFVAR,NOFVERT)
C     Upon return:
C     -----------
C     KPOS(1:NOFVAR,1:NOFVAR,1:NOFVERT-1) positive inflow parameters
C                   in the NOFVERT-1 vertices of the boundary face
C     NODRES(1:NOFVAR,1:NOFVERT-1) nodal residual due to the incoming
C                   characteristics in the NOFVERT-1 vertices 
C                   of the boundary face
C     TSTEP(NOFVAR,1:NOFVERT) contribution to the time-step in the
C                   vertices of the current element
C
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION SUM,VOLUME
      INTEGER IFAIL,IVAR,IVERT,NOFEQN
C     ..
C     .. Local Arrays ..
C
C     matrices KMAT,DUMMY,etc. must be dimensioned at least as
C     KPOS(NOFVAR,NOFVAR,*) because of the call to MatSplit
C
      DOUBLE PRECISION DUMMY(25),KMAT(25),KNEG(25),VLEFT(25),VRIGHT(25),
     +                 WNEG(4),WPOS(4),WR(4)
C     ..
C     .. External Subroutines ..
      EXTERNAL DAXPY,DCOPY,DGEMM,LINEARIZE,MATSPLITVIII
C     ..
      VOLUME = ONE
      NOFEQN = NDIM+1
      CALL LINEARIZE(IELEM,.FALSE.,VCN,NDIM,NOFVERT,VCZ,NOFVAR,VOLUME)
C
      CALL MATSPLITVIII(IELEM,NDIM,NOFVAR,VCN(1,NOFVERT),DUMMY,NOFVAR,
     +                  KMAT,KPOS,KNEG,VLEFT,VRIGHT,NOFVAR,WR,WPOS,WNEG,
     +                  .TRUE.)
Cifdef PRINT
C     CALL X04CAF('General',' ',NOFVAR,NOFVAR,KPOS,
C    +NOFVAR,'KPOSITIVE matrix ',IFAIL)
Cendif
C
      SUM = ZERO
      DO 3 IVAR = 1,NOFEQN
          SUM = SUM + WPOS(IVAR)
    3 CONTINUE
C
C     Looping over the vertices of the boundary face;
C     notice that the vertex facing the boundary face is NOFVERT
C     and is skipped in the following loop
C     Compute U(i) := [U(i) - U(\infty)] and store in DUMMY
C
      CALL DCOPY(NOFVAR*(NOFVERT-1),VCZ,1,DUMMY,1)
      DO 50 IVERT = 1,NOFVERT - 1
          CALL DAXPY(NOFVAR,-ONE,U_INFTY,1,DUMMY((IVERT-1)*NOFVAR+1),1)
          TSTEP(1,IVERT) = SUM
          IF (IVERT.GT.1) CALL DCOPY(NOFVAR*NOFVAR,KPOS(1,1,1),1,
     +                               KPOS(1,1,IVERT),1)
   50 CONTINUE
C
C     Compute NODRES(i) = -K(+)*[U(i) - U(\infty)]
C
      CALL DGEMM('No','No',NOFEQN,NOFVERT-1,NOFEQN,-ONE,KPOS,NOFVAR,
     +           DUMMY,NOFVAR,ZERO,NODRES,NOFVAR)
C
C
      RETURN

      END
