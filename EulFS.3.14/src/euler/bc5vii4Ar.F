      SUBROUTINE BC5VII4Ar(IELEM,VCN,VCZ,VCB,KPOS,NODRES,TSTEP,NDIM,
     +                  NOFVERT,NOFVAR)
C
C     $Id: bc5vii4Ar.F,v 1.3 2020/03/28 09:51:14 abonfi Exp $
C
C#define PRINT
C
      IMPLICIT NONE
C
C     Purpose:
C     -------
C     compute far-field boundary conditions for plasma
C     flows (ghost cell approach)
C     here we use conservative variables
C
      include 'paramt.h'
      include 'constants.h'
      include 'implicit.h'
      include 'bnd.h'
      include 'plasma.h'
      include 'stream.com'
      include 'time.com'
      include 'three.com'
C
      INTEGER MAXNORDER,MAXORDSQR,LENA,LENB
      PARAMETER(MAXNORDER=NSP+4,MAXORDSQR=MAXNORDER*MAXNORDER)
      PARAMETER(LENA=MAXNOFEQN*MAXNOFEQN,LENB=LENA*(MAXNOFVERT-1))
C
C     .. Scalar Arguments ..
      INTEGER IELEM,NDIM,NOFVAR,NOFVERT
C
C     On entry:
C     --------
C     IELEM  current boundary element
C     NDIM   dimension of the space (2 or 3)
C     NOFVAR number of variables (degrees of freedom)
C            in each meshpoint
C     NOFVERT number of vertices per element (=NDIM+1, since
C            only triangles or tetrahedra are allowed)
C
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION KPOS(NOFVAR,NOFVAR,*),
     +                 NODRES(NOFVAR,*),TSTEP(NOFVAR,*),
     +                 VCN(NDIM,NOFVERT),VCZ(NOFVAR,*),VCB(NDIM,*)
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
      INTEGER IFAIL,IVAR,IVERT,I,J,ILOCA,NORDER,NOFEQN,IADD
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION KMAT(MAXORDSQR),KNEG(MAXORDSQR),
     >                 VLEFT(MAXORDSQR),VRIGHT(MAXORDSQR),
     +                 WNEG(MAXNORDER),WPOS(MAXNORDER),WR(MAXNORDER),
     >                 DVDZ(LENA),DUDS(LENA),DUDZ(LENA),
     4                 TEMPB(MAXNOFVAR*MAXNOFVAR), ! a bit of extra room
     +                 DZDU(LENB), TEMPA(MAXNOFVAR*MAXNOFVAR)          
C
C
C#define PRINT
C
      DATA DVDZ/LENA*ZERO/
      DATA DUDZ/LENA*ZERO/
      DATA DUDS/LENA*ZERO/
      DATA DZDU/LENB*ZERO/
C     ..
C     .. External Subroutines ..
      EXTERNAL DAXPY,DCOPY,DGEMM,LINEARIZE,MATSPLITXII,DDOT
C     ...
      VOLUME = ONE
      CALL DINIT(NOFVERT*(NOFVAR**2),ZERO,KPOS,1)
      CALL LINEARIZE(IELEM,LALE,VCN,VCB,NDIM,NOFVERT,VCZ,NOFVAR,
     2VOLUME)
C
      CALL PARM2PRIM4Ar(NDIM,IELEM)              
C
      NOFEQN=NDIM+NSP+1
      NORDER=NOFEQN
      CALL MATSPLITXII(IELEM,NDIM,NORDER,VCN(1,NOFVERT),TEMPB,NORDER,
     +                 KMAT,KPOS,KNEG,VLEFT,VRIGHT,NORDER,WR,WPOS,WNEG,
     +                 .TRUE.)     
C
C
#ifdef PRINT
      CALL R8Mat_Print('General',' ',nofvar,nofvert,vcz,
     +nofvar,'vcz matrix ',IFAIL)
      CALL R8Mat_Print('General',' ',ndim,nofvert,vcn,
     +ndim,'vcn matrix ',IFAIL)
      CALL R8Mat_Print('General',' ',nofvar,nofvar,kPOS,
     +nofvar,'K(+) matrix ',IFAIL)
      CALL R8Mat_Print('General',' ',norder,norder,kNEG,
     +norder,'K(-) matrix ',IFAIL)
#endif
C
C     timestep
C
C
C     Compute transformation matrices
C
      CALL PARM2CONS4Ar(ZAVG,DUDZ,NOFEQN,NDIM)
!     write(13,*)ielem,(DVDZ(iadd),iadd=1,NOFEQN*NOFEQN)
C
C     K(+) := dU/dV * K(+) is stored in TEMPB(1:NOFEQN,1:NOFEQN)
C
C
C
C     K(+) := K(+) * dU/dZ is stored in TEMPA(1:NOFEQN,1:NOFEQN)
C
      CALL DGEMM('No','No',NOFEQN,NOFEQN,NOFEQN,ONE,KPOS(1,1,1),NOFEQN,
     +           DUDZ,NOFEQN,ZERO,TEMPA,NOFEQN)
#ifdef PRINT
      CALL R8Mat_Print('General',' ',NOFEQN,NOFEQN,dUdZ,
     +NOFEQN,'dUdZ matrix ',IFAIL)
      CALL R8Mat_Print('General',' ',NOFEQN,NOFEQN,TEMPA,
     +NOFEQN,' K(+) * dU/dZ matrix ',IFAIL)
#endif
!     write(14,*)ielem,(TEMPA(j),j=1,NOFEQN*NOFEQN)
C
C     REM: TEMPB is dimensioned NOFVAR*NOFVERT
C
      CALL DCOPY(NOFVAR*(NOFVERT-1),VCZ,1,TEMPB,1)
C
C     Looping over the vertices of the boundary face;
C     notice that the vertex facing the boundary face is NOFVERT
C     and is skipped in the following loop
C     Compute U(i) := [U(i) - U(\infty)] \forall i=1,..,NOFVERT-1
C     and store in TEMPB
C
      DO 50 IVERT = 1,NOFVERT - 1
          IADD = (IVERT-1)*NOFVAR+1
          CALL DAXPY(NOFVAR,MONE,U_INFTY,1,TEMPB(IADD),1)
   50 CONTINUE
C
      IF( CHAR_TIMESTEPPING )THEN
         DO IVERT = 1,NOFVERT-1
            HELP = ZERO
            DO IVAR = 1,NOFEQN 
               HELP = MAX(HELP,WPOS(IVAR))
            ENDDO 
            DO IVAR = 1,NOFEQN 
               TSTEP(IVAR,IVERT) = HELP
            ENDDO 
         ENDDO
      ELSE ! "standard" implementation
         HELP = ZERO
         DO IVAR = 1, NOFEQN
            HELP = HELP + WPOS(IVAR) 
         ENDDO
         DO IVERT = 1,NOFVERT-1 
            DO IVAR = 1,NOFEQN 
               TSTEP(IVAR,IVERT) = HELP
            ENDDO 
         ENDDO
      ENDIF
      
C      write(6,*)'TEMPB=',TEMPB
C
C     Compute -K(+)*[U(i) - U(\infty)] as
C             -K(+) * dU/dZ * [Z(i) - Z(\infty)]
C            ( -dU/dV * K(+) * dV/dZ * [Z(i) - Z(\infty)])
C
      CALL DGEMM('No','No',NOFEQN,NOFVERT-1,NOFEQN,MONE,TEMPA,NOFEQN,
     +           TEMPB,NOFVAR,ZERO,NODRES,NOFVAR)
#ifdef PRINT
      CALL R8Mat_Print('General',' ',nofvar,nofvert,TEMPB,
     +nofvar,'TEMPB matrix ',IFAIL)
      CALL R8Mat_Print('General',' ',nofvar,nofvert,NODRES,
     +nofvar,'nodres matrix ',IFAIL)
#endif
!     write(18,*)ielem,((nodres(iadd,ivert),iadd=1,NOFVAR),
!    &ivert=1,NOFVERT)
C
      IF(.NOT.TIMEIMPL)RETURN
C
C     K(+) := K(+) * 2. * (dZ/dU)(ivert)
C
      DO 55 IVERT = 1,NOFVERT - 1
          ILOCA = (IVERT-1)*NOFEQN*NOFEQN+1
          IF(TIMEIMPL) 
C
C     compute dZdU(ivert)
     +    CALL MatdZdU4Ar(VCZ(1,IVERT),DZDU(ILOCA),NDIM,NOFEQN)
          CALL DGEMM('No','No',NOFEQN,NOFEQN,NOFEQN,TWO,TEMPA,NOFEQN,
     +               DZDU(ILOCA),NOFEQN,ZERO,KPOS(1,1,IVERT),
     +               NOFVAR)
#ifdef PRINT
      CALL R8Mat_Print('General',' ',nofvar,nofvar,kpos(1,1,ivert),
     +nofvar,'kpos matrix ',IFAIL)
       write(6,*)ivert,ielem
      pause 
#endif
   55 CONTINUE
C
      RETURN

      END
