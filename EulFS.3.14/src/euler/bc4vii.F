      SUBROUTINE BC4VII(IELEM,VCN,VCZ,VCB,KPOS,NODRES,TSTEP,NDIM,
     +                  NOFVERT,NOFVAR)
C
C     $Id: bc4vii.F,v 1.13 2020/03/28 09:51:14 abonfi Exp $
C
CCCCC#define PRINT
C
      IMPLICIT NONE
C
C     Purpose:
C     -------
C     compute inviscid wall boundary conditions for compressible
C     flows (ghost cell approach)
C     here we use conserved variables
C
      include 'paramt.h'
      include 'constants.h'
      include 'implicit.h'
      include 'bnd.h'
      include 'dofs.com'
      include 'stream.com'
      include 'time.com'
      include 'three.com'
C
      INTEGER MAXNORDER,MAXORDSQR,LENA,LENB
      PARAMETER(MAXNORDER=5,MAXORDSQR=MAXNORDER*MAXNORDER)
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
C             only triangles or tetrahedra are allowed)
C
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION KPOS(NOFVAR,NOFVAR,*),
     +                 NODRES(NOFVAR,*),TSTEP(NOFVAR,*),
     +                 VCN(NDIM,NOFVERT),VCZ(NOFVAR,*),VCB(NDIM,NOFVERT)
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
      DOUBLE PRECISION HELP,VOLUME,ALPHA
      DOUBLE PRECISION ddot,DNRM2
      INTEGER IFAIL,IVAR,IVERT,I,J,ILOCA,NORDER,NOFEQN,IADD
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION KMAT(MAXORDSQR),KNEG(MAXORDSQR),
     >                 VLEFT(MAXORDSQR),VRIGHT(MAXORDSQR),
     +                 WNEG(MAXNORDER),WPOS(MAXNORDER),WR(MAXNORDER),
     >                 DUDZ(LENA),
     4                 TEMPB(MAXNOFVAR*MAXNOFVAR), ! a bit of extra room
     +                 DZDU(LENB),VARRAY(MAXNOFVAR),
     7                 CG7(MAXNOFEQN,MAXNOFEQN)
C
C
C#define PRINT
C
      DATA DUDZ/LENA*ZERO/
      DATA DZDU/LENB*ZERO/
C     ..
C     .. External Subroutines ..
      EXTERNAL DAXPY,DCOPY,DGEMM,LINEARIZE,MATSPLITVII,DDOT
C     ..
C     here the averaged state is computed using the states of the boundary
C     cell including the vertex which is opposite to the bndry face
C
      NOFEQN = NDIM+2
      NORDER = NOFEQN
      VOLUME = ONE
C
      CALL DINIT(NOFVERT*(NOFVAR**2),ZERO,KPOS,1)
C
      DO 50 IVERT = 1, (NOFVERT-1) ! loop over the d vertices of the bndry face
#ifdef PRINT
         write(6,*)'ivert = ',ivert
         CALL R8Mat_Print('General',' ',nofvar,nofvert,vcz,
     +nofvar,'vcz matrix (before)',IFAIL)
         CALL R8Mat_Print('General',' ',ndim,nofvert,vcn,
     +ndim,'vcn matrix ',IFAIL)
#endif
C
C     compute the starred state and put it into VCZ(*,NOFVERT): this is needed to compute the averaged state 
C
         CALL CALCPZ2(VCZ(1,IVERT),NDIM,CG7,MAXNOFEQN,VCN(1,NOFVERT))
         CALL DGEMV('No',NOFEQN,NOFEQN,ONE,CG7,MAXNOFEQN,VCZ(1,IVERT),1,
     &              ZERO,VCZ(1,NOFVERT),1)
C
         CALL LINEARIZE(IELEM,LALE,VCN,VCB,NDIM,NOFVERT,VCZ,
     2                  NOFVAR,VOLUME)
C
         CALL PARM2PRIM(NDIM,IELEM)
C
         CALL MATSPLITXI(IELEM,NDIM,NORDER,VCN(1,NOFVERT),TEMPB,NORDER,
     +                KMAT,KPOS(1,1,IVERT),KNEG,VLEFT,VRIGHT,NORDER,WR,
     +                WPOS,WNEG,.TRUE.)
C
#ifdef PRINT
      CALL R8Mat_Print('General',' ',NOfeqn,NOfeqn,CG7,
     +MAXNOFEQN,'CG7 matrix ',IFAIL)
      CALL R8Mat_Print('General',' ',nofvar,nofvert,vcz,
     +nofvar,'vcz matrix ',IFAIL)
      CALL R8Mat_Print('General',' ',nofvar,nofvar,KPOS(1,1,IVERT),
     +nofvar,'K(+) matrix ',IFAIL)
      CALL R8Mat_Print('General',' ',norder,norder,kneg,
     +norder,'K(-) matrix ',IFAIL)
      pause 
#endif
C
C     timestep
C
         IF( CHAR_TIMESTEPPING )THEN
            HELP = ZERO
            DO IVAR = 1,NORDER 
               HELP = MAX(HELP,WPOS(IVAR))
            ENDDO 
         ELSE ! "standard" implementation
            HELP = ZERO
            DO IVAR = 1, NORDER
               HELP = HELP + WPOS(IVAR) 
            ENDDO
         ENDIF
         DO 4 IVAR = 1,NORDER
             TSTEP(IVAR,IVERT) = HELP
    4    CONTINUE
C
C
C     Compute transformation matrices
C
         CALL PARM2CONS(ZAVG,DUDZ,NOFEQN,NDIM)
!     write(12,*)ielem,(DUDZ(iadd),iadd=1,NOFEQN*NOFEQN)
C
C     compute tempA : = (dU/dZ) x (I-CG7)
C
         CALL DINIT(NOFEQN,ZERO,VARRAY,1)
         CALL DCOPY(NDIM,VCN(1,NOFVERT),1,VARRAY(IX),1)
         HELP = DNRM2(NDIM,VARRAY(IX),1)
         CALL DSCAL(NDIM,ONE/HELP,VARRAY(IX),1) ! make n a unit vector and store in varray
         ALPHA = TWO * DDOT(NDIM,VCZ(IX,IVERT),1,VARRAY(IX),1)
         IF(LALE)THEN
            ALPHA = ALPHA - VCZ(1,IVERT) * 
     &              DDOT(NDIM,VCB(1,IVERT),1,VARRAY(IX),1)
         ENDIF
C
C     compute tempB := - Kpos x (dU/dZ) 
C
         CALL DGEMM('No','No',NOFEQN,NOFEQN,NOFEQN,MONE, 
     +           KPOS(1,1,IVERT),NOFEQN,DUDZ,NOFEQN,ZERO,TEMPB,NOFEQN)
C
C     compute
C
         CALL DGEMV('No',NOFEQN,NOFEQN,ALPHA,TEMPB,NOFEQN,VARRAY,1,
     &           ZERO,NODRES(1,IVERT),1)
C
C
!     write(11,*)ielem,((KPOS(i,j,1),j=1,NOFVAR),i=1,NOFVAR)
!     write(16,*)ielem,(tempb(iadd),iadd=1,NOFEQN*NOFEQN)
#ifdef PRINT
      CALL R8Mat_Print('General',' ',NOFEQN,NOFEQN,DUDZ,
     +NOFEQN,'dU/dZ ',IFAIL)
      CALL R8Mat_Print('General',' ',NOFEQN,NOFEQN,TEMPB,
     +NOFEQN,'- K x dudZ ',IFAIL)
      CALL R8Mat_Print('General',' ',nofvar,nofvert,NODRES,
     +nofvar,'nodres matrix ',IFAIL)
#endif
C
      IF(.NOT.TIMEIMPL)GOTO 50
C
C     K(+) := K(+) * 2. * (dZ/dU)(ivert)
C
!     DO 55 IVERT = 1,NOFVERT - 1
          ILOCA = (IVERT-1)*NOFEQN*NOFEQN+1
          IF(PICARD) THEN
C
C     compute dZdU(ivert)
          CALL MatdZdU(VCZ(1,IVERT),DZDU(ILOCA),NDIM,NOFEQN)
          CALL DGEMM('No','No',NOFEQN,NOFEQN,NOFEQN,-TWO,TEMPB,NOFEQN,
     +               DZDU(ILOCA),NOFEQN,ZERO,KPOS(1,1,IVERT),
     +               NOFVAR)
          ENDIF
#ifdef PRINT
       write(6,*)ivert,TIMEIMPL,PICARD,NEWTON
      CALL R8Mat_Print('General',' ',nofvar,nofvar,kpos(1,1,ivert),
     +nofvar,'kpos matrix ',IFAIL)
      write(6,*)'IE = ',ielem,' iv = ',ivert
      pause 
#endif
   50 CONTINUE ! end loop over vertices
C
      RETURN

      END
!     cg5(3,1) = -2 * l * (l * u + m * v + w * n)
!     cg5(3,3) = 2 * l ** 2
!     cg5(3,4) = 2 * l * m
!     cg5(3,5) = 2 * l * n
!     cg5(4,1) = -2 * m * (l * u + m * v + w * n)
!     cg5(4,3) = 2 * l * m
!     cg5(4,4) = 2 * m ** 2
!     cg5(4,5) = 2 * n * m
!     cg5(5,1) = -2 * n * (l * u + m * v + w * n)
!     cg5(5,3) = 2 * l * n
!     cg5(5,4) = 2 * n * m
!     cg5(5,5) = 2 * n ** 2
C
      SUBROUTINE CALCPZ2(Z,NDIM,CG7,LDA,AN)
C
C
C#define PRINT
C
      IMPLICIT NONE
C
C     Purpose:
C     -------
C     compute inviscid wall boundary conditions for compressible
C     flows (ghost cell approach)
C     here we use conserved variables
C
      include 'paramt.h'
      include 'constants.h'
C
C
C     .. Scalar Arguments ..
      INTEGER NDIM,LDA
      DOUBLE PRECISION CG7(LDA,*),AN(*),Z(*)
      DOUBLE PRECISION L,M,N,ANN
C
      L = AN(1)
      M = AN(2)
      IF(NDIM.EQ.3)THEN
         N = AN(3)
      ELSE
         N = ZERO
      ENDIF
      ANN = ONE/SQRT(L*L+M*M+N*N)
      L = L*ANN
      M = M*ANN
      N = N*ANN
C
      cg7(1,1) = ONE
      cg7(1,2) = ZERO
      cg7(1,3) = ZERO
      cg7(1,4) = ZERO
C
      cg7(2,1) = ZERO
      cg7(2,2) = ONE
      cg7(2,3) = ZERO
      cg7(2,4) = ZERO
C
      cg7(3,1) = ZERO
      cg7(3,2) = ZERO
      cg7(3,3) = ONE - TWO * l ** 2
      cg7(3,4) = -TWO * l * m
C
      cg7(4,1) = ZERO
      cg7(4,2) = ZERO
      cg7(4,3) = cg7(3,4)
      cg7(4,4) = ONE - TWO * m ** 2
C
      IF(NDIM.EQ.2)RETURN
C
      cg7(1,5) = ZERO
      cg7(2,5) = ZERO
      cg7(3,5) = -TWO * l * n
      cg7(4,5) = -TWO * n * m
C
      cg7(5,1) = ZERO
      cg7(5,2) = ZERO
      cg7(5,3) = -TWO * l * n
      cg7(5,4) = -TWO * n * m
      cg7(5,5) = ONE - TWO * n ** 2
C
      RETURN
      END
C
      SUBROUTINE CALCPZ(Z,NDIM,CG7,LDA,AN)
C
C
C#define PRINT
C
      IMPLICIT NONE
C
C     Purpose:
C     -------
C     compute inviscid wall boundary conditions for compressible
C     flows (ghost cell approach)
C     here we use conserved variables
C
      include 'paramt.h'
      include 'constants.h'
      INCLUDE 'pfcgas.com'
C
C
C     .. Scalar Arguments ..
      INTEGER NDIM,LDA
      DOUBLE PRECISION CG7(LDA,*),AN(*),Z(*)
      DOUBLE PRECISION L,M,N,ANN,HELP,CNST
C
      L = AN(1)
      M = AN(2)
      IF(NDIM.EQ.3)THEN
         N = AN(3)
      ELSE
         N = ZERO
      ENDIF
      ANN = ONE/SQRT(L*L+M*M+N*N)
      L = L*ANN
      M = M*ANN
      N = N*ANN
C
      HELP = Z(3)*L+Z(4)*M
      IF(NDIM.EQ.3)HELP = HELP + Z(5)*N
C
      CNST = TWO * GM1 * HELP
C
      cg7(1,1) = ONE
      cg7(1,2) = ZERO
      cg7(1,3) = ZERO
      cg7(1,4) = ZERO
C
      cg7(2,1) = -CNST * HELP / (Z(1)*Z(1))
      cg7(2,2) = ONE
      cg7(2,3) = CNST * L / z(1)
      cg7(2,4) = CNST * M / z(1)
C
      cg7(3,1) = TWO / z(1) * l * HELP
      cg7(3,2) = ZERO
      cg7(3,3) = ONE - TWO * l ** 2
      cg7(3,4) = -TWO * l * m
C
      cg7(4,1) = TWO / z(1) * m * HELP
      cg7(4,2) = ZERO
      cg7(4,3) = -TWO * l * m
      cg7(4,4) = ONE - TWO * m ** 2
C
      IF(NDIM.EQ.2)RETURN
C
      cg7(1,5) = ZERO
      cg7(2,5) = CNST * N / z(1)
      cg7(3,5) = -TWO * l * n
      cg7(4,5) = -TWO * n * m
C
      cg7(5,1) = TWO / z(1) * n * HELP
      cg7(5,2) = ZERO
      cg7(5,3) = -TWO * l * n
      cg7(5,4) = -TWO * n * m
      cg7(5,5) = ONE - TWO * n ** 2
C
      RETURN
      END
