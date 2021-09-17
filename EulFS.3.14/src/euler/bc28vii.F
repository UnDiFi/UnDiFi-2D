      SUBROUTINE BC28VII(IELEM,VCN,VCZ,VCB,PCN,KPOS,NODRES,TSTEP,NDIM,
     +                  NOFVERT,NOFVAR,IBC)
C
C     $Id: bc28vii.F,v 1.14 2020/03/28 09:51:14 abonfi Exp $ 
C
C#define PRINT
C
      IMPLICIT NONE
C
C     Purpose:
C     -------
C     compute in-outflow boundary conditions for compressible
C     flows (ghost cell approach)
C
      include 'paramt.h'
      include 'constants.h'
      include 'implicit.h'
      include 'bnd.h'
      include 'ibc8.com'
      include 'stream.com'
      include 'time.com'
      include 'three.com'
C
C     .. Scalar Arguments ..
      INTEGER IELEM,NDIM,NOFVAR,NOFVERT,IBC
C
C     On entry:
C     --------
C     IELEM  current boundary element
C     NDIM   dimension of the space (2 or 3)
C     NOFVAR number of variables (degrees of freedom)
C            in each meshpoint
C     NOFVERT number of vertices per element (=NDIM+1, since
C             only triangles or tetrahedra are allowed)
C
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION KPOS(NOFVAR,NOFVAR,*),
     +       NODRES(NOFVAR,*),TSTEP(NOFVAR,*),PCN(nVarsInlet,NOFVERT),
     +       VCN(NDIM,NOFVERT),VCZ(NOFVAR,NOFVERT),VCB(NDIM,NOFVERT)
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
      DOUBLE PRECISION HELP,VOLUME
      DOUBLE PRECISION ddot
      INTEGER IVAR,IVERT,ILOCA,NORDER,NOFEQN
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION DUMMY(25),KMAT(16),KNEG(16),VLEFT(16),VRIGHT(16),
     +                 WNEG(4),WPOS(4),WR(4),DVDZ(25),DUDS(25),
     +                 DZDU(75),TEMPA(25),
     +                 TEMPB(NMAX*NMAX),WORK(6),wksp(NMAX*NMAX)
C
      DATA DVDZ,DUDS,DZDU,WORK/125*ZERO,6*ONE/
C     ..
C     .. External Subroutines ..
      EXTERNAL DAXPY,DCOPY,DGEMM,LINEARIZE,MATSPLITVII,DDOT
C     ..
#ifdef PRINT
      CALL R8Mat_Print('General',' ',nVarsInlet,NOFVERT,PCN,
     +nVarsInlet,'PCN matrix ',IFAIL)
#endif
C
      VOLUME = ONE
      CALL DINIT((NOFVAR*NOFVERT)**2,ZERO,KPOS,1)
C
      DO 55 IVERT = 1,NOFVERT - 1 ! outer loop over the vertices of the bndry face
C
C     static pressure (outlet bcs)
C
         WORK(1) = POUTLET
C
C     resp. total pressure, total temperature, n_x, n_y, n_z
C
         WORK(2) = PCN(1,IVERT)
         WORK(3) = PCN(2,IVERT)
         WORK(4) = PCN(4,IVERT)
         WORK(5) = PCN(5,IVERT)
         WORK(6) = PCN(6,IVERT)
C
         IF(    IBC.EQ.BC_TYPE_SUBS_OUTLET)THEN
            CALL GHOST2VII(VCZ(1,IVERT),VCZ(1,NOFVERT),WORK,NDIM)
         ELSEIF(IBC.EQ.BC_TYPE_SUBS_INLET)THEN
            CALL GHOST8VII(VCZ(1,IVERT),VCZ(1,NOFVERT),WORK,NDIM)
         ENDIF
C
#ifdef PRINT
      write(6,*)' ibc is ',ibc,ivert,ielem
      CALL R8Mat_Print('General',' ',NOFVAR,NOFVERT,VCZ,
     +NOFVAR,'VCZ matrix ',IFAIL)
      do i = 1,nofvar
      write(6,FMT="(4(F20.12,1X))")(VCZ(i,j),j=1,nofvert)
      enddo
#endif
C
         CALL LINEARIZE(IELEM,LALE,VCN,VCB,NDIM,NOFVERT,VCZ,NOFVAR,
     2                  VOLUME)
C
         CALL PARM2PRIM(NDIM,IELEM)
C
         NORDER=NDIM+1
         NOFEQN=NDIM+2
         CALL MATSPLITVII(IELEM,NDIM,NORDER,VCN(1,NOFVERT),DUMMY,NORDER,
     +                 KMAT,TEMPA,KNEG,VLEFT,VRIGHT,NORDER,WR,WPOS,WNEG,
     +                 .TRUE.)
C
         HELP = DDOT(NDIM,UAVG(3),1,VCN(1,NOFVERT),1)/REAL(NDIM)
         HELP = HALF*(HELP+ABS(HELP))
         KPOS(1,1,IVERT) = HELP
C
C     timestep
C
      IF( CHAR_TIMESTEPPING )THEN
         HELP = MAX(ZERO,HELP) ! be careful: help must be KPOS(1,1,1)
         DO IVAR = 1,NORDER 
            HELP = MAX(HELP,WPOS(IVAR))
         ENDDO 
      ELSE ! "standard" implementation
         HELP = ZERO
         DO IVAR = 1, NORDER
            HELP = HELP + WPOS(IVAR) 
         ENDDO
      ENDIF ! test characteristic timestepping
      DO IVAR = 1,NOFEQN  ! =NORDER+1
         TSTEP(IVAR,IVERT) = HELP
      ENDDO 
C
C     Copy the matrix in symmetrizing variables into TEMPA
C
         CALL MATINS(KPOS(1,1,IVERT),NOFVAR,TEMPA,NORDER,1,1,1)
C
#ifdef PRINT
      CALL R8Mat_Print('General',' ',NORDER,NORDER,TEMPA,
     +NORDER,'Reduced K matrix ',IFAIL)
      CALL R8Mat_Print('General',' ',NOFVAR,NOFVAR,KPOS(1,1,IVERT),
     +NOFVAR,'KPOS(+) matrix ',IFAIL)
#endif
C
C     Compute transformation matrices
C
         CALL SYMM2CONS(ZAVG,DUDS,NOFEQN,NDIM)
C
         CALL PARM2SYMM(ZAVG,DVDZ,NOFEQN,NDIM)
C
#ifdef PRINT
      CALL R8Mat_Print('General',' ',NOFEQN,NOFEQN,DUDS,
     +NOFEQN,'DUDS matrix ',IFAIL)
      CALL R8Mat_Print('General',' ',NOFEQN,NOFEQN,DVDZ,
     +NOFEQN,'DVDZ matrix ',IFAIL)
#endif
C
C     K(+) := dU/dV * K(+)
C
         CALL DGEMM('No','No',NOFEQN,NOFEQN,NOFEQN,ONE,DUDS,NOFEQN,
     +           KPOS(1,1,IVERT),NOFVAR,ZERO,TEMPB,NOFEQN)
C
#ifdef PRINT
      CALL R8Mat_Print('General',' ',NOFEQN,NOFEQN,TEMPB,
     +NOFEQN,'TEMPB matrix ',IFAIL)
#endif
C
C     K(+) := dU/dV * K(+) * dV/dZ
C
         CALL DGEMM('No','No',NOFEQN,NOFEQN,NOFEQN,ONE,TEMPB,NOFEQN,
     +           DVDZ,NOFEQN,ZERO,TEMPA,NOFEQN)
#ifdef PRINT
      CALL R8Mat_Print('General',' ',NOFEQN,NOFEQN,TEMPA,
     +NOFEQN,'TEMPA matrix ',IFAIL)
#endif
C
C     Looping over the vertices of the boundary face;
C     notice that the vertex facing the boundary face is NOFVERT
C     and is skipped in the following loop
C     Compute U(i) := [U(i) - U(*)] and store in TEMPB
C
         DO 245 IVAR = 1, NOFEQN
            wksp((IVERT-1)*NOFVAR+IVAR) = 
     +      VCZ(IVAR,IVERT) - VCZ(IVAR,NOFVERT) 
  245    CONTINUE
C
C
C     Compute -K(+)*[U(i) - U(*)]
C
      CALL DGEMV('No',NOFEQN,NOFEQN,MONE,TEMPA,NOFEQN,
     +     wksp((IVERT-1)*NOFVAR+1),1,ZERO,NODRES(1,IVERT),1)
C
C     K(+) := K(+) * 2. * (dZ/dU)(ivert)
C
          ILOCA = (IVERT-1)*NOFEQN*NOFEQN+1
          IF(TIMEIMPL) 
     +    CALL MatdZdU(VCZ(1,IVERT),DZDU(ILOCA),NDIM,NOFEQN)
          CALL DGEMM('No','No',NOFEQN,NOFEQN,NOFEQN,TWO,TEMPA,NOFEQN,
     +               DZDU(ILOCA),NOFEQN,ZERO,KPOS(1,1,IVERT),
     +               NOFVAR)
#ifdef PRINT
      write(6,*)' ibc is ',ibc
      CALL R8Mat_Print('General',' ',NOFVAR,NOFVAR,KPOS(1,1,IVERT),
     +NOFVAR,'K(+) matrix ',IFAIL)
      CALL R8Mat_Print('General',' ',NOFVAR,NOFVert,VCZ,
     +NOFVAR,'VCZ matrix ',IFAIL)
!     pause
#endif
   55 CONTINUE ! End loop over the vertices
#ifdef PRINT
      CALL R8Mat_Print('General',' ',NOFVAR,NOFVERT-1,wksp,
     +NOFVAR,'delta matrix ',IFAIL)
      CALL R8Mat_Print('General',' ',NOFVAR,NOFVERT,NODRES,
     +NOFVAR,'Residual matrix ',IFAIL)
#endif
C
      RETURN

      END
