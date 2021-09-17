      SUBROUTINE BC28VII4Ar(IELEM,VCN,VCZ,VCB,PCN,KPOS,NODRES,TSTEP,
     +                   NDIM,NOFVERT,NOFVAR,IBC)
C
C     $Id: bc28vii4Ar.F,v 1.3 2020/03/28 09:51:14 abonfi Exp $ 
C
C#define PRINT
C
      IMPLICIT NONE
C
C     Purpose:
C     -------
C     compute in-outflow boundary conditions for plasma
C     flows (ghost cell approach)
C     here we use conservative variables
C
      include 'paramt.h'
      include 'constants.h'
      include 'implicit.h'
      include 'bnd.h'
      include 'ibc8.com'
      include 'plasma.h'
      include 'stream.com'
      include 'time.com'
      include 'three.com'
C
      INTEGER MAXNORDER,MAXORDSQR,LENA,LENB,LENC
      PARAMETER(MAXNORDER=NSP+4,MAXORDSQR=MAXNORDER*MAXNORDER)
      PARAMETER(LENA=MAXNOFEQN*MAXNOFEQN,LENB=LENA*(MAXNOFVERT-1))
      PARAMETER(LENC=LENA+LENA+LENB)
C
C     .. Scalar Arguments ..
      INTEGER IELEM,NDIM,NOFVAR,NOFVERT,IBC,i,j,IFAIL
C
C     On entry:
C     --------
C     IELEM  current boundary element
C     NDIM   dimension of the space (2 or 3)
C     NOFVAR number of variables (degrees of freedom)
C            in each meshpoint
C     NOFVERT number of vertices per element (=NDIM+1, since
C             only triangles or tetrahedra are t) = value(rr,rcline,zb,allowed)
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
      DOUBLE PRECISION VOLUME,HELP
      DOUBLE PRECISION ddot
      INTEGER IVAR,IVERT,ILOCA,NORDER,NOFEQN
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION DUMMY(LENA),KMAT(MAXORDSQR),KNEG(MAXORDSQR),
     +                 VLEFT(MAXORDSQR),VRIGHT(MAXORDSQR),
     +                 WNEG(MAXNORDER),WPOS(MAXNORDER),WR(MAXNORDER),
     +                 DVDZ(LENA),DUDS(LENA),DUDZ(LENA),
     +                 TEMPA(MAXNOFVAR*MAXNOFVAR),
     +                 TEMPB(MAXNOFVAR*MAXNOFVAR),
     +                 DZDU(LENB),WORK(6),wksp(NMAX*NMAX)
C
      DATA DVDZ,DUDZ,DZDU,WORK/LENC*ZERO,6*ONE/
C     ..
C     .. External Subroutines ..
      EXTERNAL DAXPY,DCOPY,DGEMM,LINEARIZE,MATSPLITXII,DDOT
C     ..
#ifdef PRINT
      CALL R8Mat_Print('General',' ',nVarsInlet,NOFVERT,PCN,
     +nVarsInlet,'PCN matrix ',IFAIL)
#endif
C
      VOLUME = ONE
      CALL DINIT((NOFVAR*NOFVERT)**2,ZERO,KPOS,1)
C
      DO 55 IVERT = 1,NOFVERT - 1 ! loop over the d vertices of the bndry face
C
C     static pressure (outlet bcs)
C
         WORK(1) = POUTLET
C
C     resp. total pressure, total temperature, n_x, n_y, n_z
C
C     Verificare le dimensioni di PCN e WORK per il plasma 
         WORK(2) = PCN(1,IVERT)
         WORK(3) = PCN(2,IVERT)
         WORK(4) = PCN(4,IVERT)
         WORK(5) = PCN(5,IVERT)
         WORK(6) = PCN(6,IVERT)
C
         IF(IBC.EQ.BC_TYPE_SUBS_OUTLET)THEN
            CALL GHOST2VII4Ar(VCZ(1,IVERT),VCZ(1,NOFVERT),WORK,NDIM)
         ELSEIF(IBC.EQ.BC_TYPE_SUBS_INLET)THEN
            CALL GHOST8VII4Ar(VCZ(1,IVERT),VCZ(1,NOFVERT),WORK,NDIM)
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
         CALL PARM2PRIM4Ar(NDIM,IELEM)
C
         NOFEQN=NDIM+1+NSP
         NORDER=NOFEQN
         CALL MATSPLITXII(IELEM,NDIM,NORDER,VCN(1,NOFVERT),TEMPB,NORDER,
     +                 KMAT,KPOS,KNEG,VLEFT,VRIGHT,NORDER,WR,WPOS,WNEG,
     +                 .TRUE.)
C
C
#ifdef PRINT
      CALL R8Mat_Print('General',' ',NOFVAR,NOFVAR,KPOS(1,1,IVERT),
     +NOFVAR,'KPOS(+) matrix ',IFAIL)
#endif
C
C     timestep
C
         HELP = ZERO
         DO IVAR = 1,NORDER
            IF( CHAR_TIMESTEPPING )THEN 
               HELP = MAX(HELP,WPOS(IVAR))
            ELSE
               HELP = HELP + WPOS(IVAR)
            ENDIF
         ENDDO
         DO IVAR = 1,NORDER
            TSTEP(IVAR,IVERT) = HELP
         ENDDO
C
C     Compute transformation matrices
C
         CALL PARM2CONS4Ar(ZAVG,DUDZ,NOFEQN,NDIM)
C
C
C     K(+) := K(+) * dU/dZ is stored in TEMPA(1:NOFEQN,1:NOFEQN)r.F
C
         CALL DGEMM('No','No',NOFEQN,NOFEQN,NOFEQN,ONE,KPOS(1,1,1),
     +           NOFEQN,DUDZ,NOFEQN,ZERO,TEMPA,NOFEQN)
C
#ifdef PRINT
      CALL R8Mat_Print('General',' ',NOFEQN,NOFEQN,DUDZ,
     +NOFEQN,'dUdZ matrix ',IFAIL)
      CALL R8Mat_Print('General',' ',NOFEQN,NOFEQN,TEMPA,
     +NOFEQN,'K(+) * dU/dZ  matrix ',IFAIL)
#endif
C
C     K(+) := dU/dV * K(+) * dV/dZ
C
C      CALL DGEMM('No','No',NOFEQN,NOFEQN,NOFEQN,ONE,TEMPB,NOFEQN,
C     +           DVDZ,NOFEQN,ZERO,TEMPA,NOFEQN)
C#ifdef PRINT
C      CALL R8Mat_Print('General',' ',NOFEQN,NOFEQN,TEMPA,
C     +            NOFEQN,'TEMPA matrix ',IFAIL)
C#endif
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
         CALL DGEMV('No',NOFEQN,NOFEQN,-ONE,TEMPA,NOFEQN,
     +     wksp((IVERT-1)*NOFVAR+1),1,ZERO,NODRES(1,IVERT),1)
C
C     K(+) := K(+) * 2. * (dZ/dU)(ivert)
C
         ILOCA = (IVERT-1)*NOFEQN*NOFEQN+1
         IF(TIMEIMPL) 
     +   CALL MatdZdU4Ar(VCZ(1,IVERT),DZDU(ILOCA),NDIM,NOFEQN)
         CALL DGEMM('No','No',NOFEQN,NOFEQN,NOFEQN,TWO,TEMPA,NOFEQN,
     +               DZDU(ILOCA),NOFEQN,ZERO,KPOS(1,1,IVERT),
     +               NOFVAR)
#ifdef PRINT
         write(6,*)' ibc is ',ibc
         CALL R8Mat_Print('General',' ',NOFVAR,NOFVAR,KPOS(1,1,IVERT),
     +NOFVAR,'K(+) matrix ',IFAIL)
       pause    
#endif
C
   55 CONTINUE ! End loop over the vertices
C
#ifdef PRINT
      CALL R8Mat_Print('General',' ',NOFVAR,NOFVERT-1,wksp,
     +NOFVAR,'delta matrix ',IFAIL)
      CALL R8Mat_Print('General',' ',NOFVAR,NOFVERT,NODRES,
     +NOFVAR,'Residual matrix ',IFAIL)
       write(6,*)ivert,ielem
      pause
#endif
C
      RETURN

      END
