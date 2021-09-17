      SUBROUTINE BC28VIII(IELEM,VCN,VCZ,VCB,PCN,KPOS,NODRES,TSTEP,NDIM,
     +                   NOFVERT,NOFVAR,IBC)
C
C     $Id: bc28viii.F,v 1.14 2020/03/28 09:51:14 abonfi Exp $
C
      IMPLICIT NONE
C
C     Purpose:
C     -------
C     compute inflow/outflow boundary conditions for INcompressible
C     flows (ghost cell approach)
C
      include 'paramt.h'
      include 'constants.h'
      include 'bnd.h'
      include 'ibc8.com'
      include 'stream.com'
      include 'time.com'
C
C     .. Scalar Arguments ..
      INTEGER IELEM,NDIM,NOFVAR,NOFVERT,IBC
C
C     IELEM  current boundary element
C     NDIM   dimension of the space (2 or 3)
C     NOFVAR number of variables (degrees of freedom)
C            in each meshpoint
C     NOFVAR number of variables (degrees of freedom)
C            in each meshpoint
C
CCC#define PRINT
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION KPOS(NOFVAR,NOFVAR,NOFVERT-1),
     +       NODRES(NOFVAR,NOFVERT-1),TSTEP(NOFVAR,NOFVERT-1),
     +       VCN(NDIM,NOFVERT),VCZ(NOFVAR,NOFVERT),VCB(NDIM,NOFVERT),
     &       PCN(nVarsInlet,NOFVERT)
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
      DOUBLE PRECISION WKSP,VOLUME,HELP
      INTEGER IFAIL,IVAR,IVERT,NOFEQN,j
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION DUMMY(MAXNOFVAR*MAXNOFVERT),KMAT(MAX_NOFVAR_SQR),
     &KNEG(MAX_NOFVAR_SQR),VLEFT(MAX_NOFVAR_SQR),VRIGHT(MAX_NOFVAR_SQR),
     &WNEG(MAXNOFVAR),WPOS(MAXNOFVAR),WR(MAXNOFVAR),WORK(9)
C     ..
      DOUBLE PRECISION DNRM2
C     ..
C     .. External Subroutines ..
      EXTERNAL DAXPY,DCOPY,DGEMV,LINEARIZE,MATSPLITVIII
C     ..
      DATA WORK/9*ONE/
C
      VOLUME = ONE
      call DINIT(MAXNOFVAR*MAXNOFVERT,MONE,DUMMY,1)
C
#ifdef PRINT
      write(6,*)' ibc is ',ibc
      CALL R8Mat_Print('General',' ',nVarsInlet,NOFVERT,PCN,
     +nVarsInlet,'PCN matrix ',IFAIL)
#endif
C     unit face normal
      WORK(7) = VCN(1,NOFVERT)
      WORK(8) = VCN(2,NOFVERT)
      WORK(9) = VCN(3,NOFVERT)
      WKSP = ONE/DNRM2(NDIM,WORK(7),1)
      WORK(7) = WORK(7)*WKSP
      WORK(8) = WORK(8)*WKSP
      WORK(9) = WORK(9)*WKSP
C
C     loop over the vertices of the boundary face
C
      DO 1000 IVERT = 1, NOFVERT-1
C
         WORK(1) = PCN(3,IVERT) ! POUTLET !
         WORK(2) = PCN(1,IVERT)
         WORK(3) = PCN(2,IVERT)
         WORK(4) = PCN(4,IVERT)
         WORK(5) = PCN(5,IVERT)
         WORK(6) = PCN(6,IVERT)
c
c     compute the starred state and store in Z(1:NOFVAR,NOFVERT) 
c
         IF(IBC.EQ.BC_TYPE_SUBS_OUTLET)THEN
             CALL GHOST2VIII(VCZ(1,IVERT),VCZ(1,NOFVERT),WORK,NDIM)
         ELSEIF(IBC.EQ.BC_TYPE_SUBS_INLET)THEN
             CALL GHOST8VIII(VCZ(1,IVERT),VCZ(1,NOFVERT),WORK,NDIM)
         ENDIF
#ifdef PRINT
      write(6,*)'IELEM = ',ielem,' IVERT = ',ivert,' IBC = ',ibc
      CALL R8Mat_Print('General',' ',NOFVAR,NOFVERT,VCZ,
     +NOFVAR,'VCZ matrix ',IFAIL)
#endif
c
c     compute the cell averaged state based upon
c     the face values and the starred (ghost) state
c
         CALL LINEARIZE(IELEM,LALE,VCN,VCB,NDIM,NOFVERT,VCZ,NOFVAR,
     &VOLUME)
C
C     compute K_i(+)
C
         NOFEQN = NDIM +1
         CALL MATSPLITVIII(IELEM,NDIM,NOFEQN,VCN(1,NOFVERT),DUMMY,
     +                  NOFVAR,KMAT,KPOS(1,1,IVERT),KNEG,VLEFT,VRIGHT,
     +                  NOFVAR,WR,WPOS,WNEG,.TRUE.)
C
#ifdef PRINT
      CALL R8Mat_Print('General',' ',NOFVAR,NOFVAR,KPOS(1,1,IVERT),
     +NOFVAR,'KPOSITIVE matrix ',IFAIL)
#endif
C
         HELP = ZERO
         DO IVAR = 1,NOFEQN
            IF( CHAR_TIMESTEPPING )THEN
                HELP = MAX(HELP,WPOS(IVAR))
            ELSE
                HELP = HELP + WPOS(IVAR)
            ENDIF
         ENDDO
         DO 3 IVAR = 1,NOFEQN
            TSTEP(IVAR,IVERT) = HELP
    3    CONTINUE
C
C     Compute U(i) := [U(i) - U(*)] and store in DUMMY
C
         CALL DCOPY(NOFVAR,VCZ(1,IVERT),1,DUMMY,1)
         CALL DAXPY(NOFVAR,MONE,VCZ(1,NOFVERT),1,DUMMY,1)
C
C     Compute NODRES(i) = -K(+)*[U(i) - U(*)]
C
         CALL DGEMV('No',NOFEQN,NOFEQN,MONE,KPOS(1,1,IVERT),
     +           NOFVAR,DUMMY,1,ZERO,NODRES(1,IVERT),1)
C
C
#ifdef PRINT
         CALL R8Mat_Print('General',' ',NOFVAR,NOFVert,DUMMY,
     +NOFVAR,'DUMMY matrix ',IFAIL)
         CALL R8Mat_Print('General',' ',NOFVAR,NOFVert,nodres,
     +NOFVAR,'nodres matrix ',IFAIL)
#endif
 1000 CONTINUE ! End loop over the d vertices of the bndry face
#ifdef PRINT
         CALL R8Mat_Print('General',' ',NOFVAR,NOFVert,nodres,
     +NOFVAR,'nodres matrix ',IFAIL)
#endif
C
      RETURN

      END
