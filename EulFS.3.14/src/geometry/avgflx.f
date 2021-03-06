      SUBROUTINE AVGFLX(IBNDFAC,NBFAC,ICELFAC,ICELNOD,NOFVERT,NELEM,
     &VCN,ZROE,NOFVAR,NDIM,NFACE,COMPRESSIBLE,ITER)
C
C     $Id: avgflx.f,v 1.3 2020/03/28 09:45:40 abonfi Exp $
C
C     This routine checks the mass flux through the boundaries
C     of the computational domain
C
C
      IMPLICIT NONE
C
      INCLUDE 'paramt.h'
      INCLUDE 'constants.h'
      INCLUDE 'bnd.h'
      INCLUDE 'bnd.com'
      INCLUDE 'periodic.com'
      INCLUDE 'stream.com'
      INCLUDE 'io.com'
C
      INTEGER NBFAC,NOFVERT,NELEM,NDIM,NFACE,NOFVAR,ITER
      INTEGER NERR,IOPT,IDOF
      CHARACTER ERRMSG*72
      LOGICAL COMPRESSIBLE
C
      INTEGER IBNDFAC(3,*),ICELFAC(NOFVERT,*),ICELNOD(NOFVERT,*)
      DOUBLE PRECISION VCN(NDIM,*),ZROE(NOFVAR,*)
C
      INTEGER J,IBC,IFAIL,JVERT
      INTEGER IELEM,IBFAC,IFREQ,IVERT,IFACE,IPOIN
      DOUBLE PRECISION DOTM,RAVG,TEMP
C
      DOUBLE PRECISION WKSP(MAXNOFVAR),ZG(MAXNOFVAR),
     2ALPHA(MAXNOFVAR),FLXN(MAXNOFVAR,MAXNOFVERT),
     3VCZ(MAXNOFVAR*MAXNOFVERT)
C
C
      INTEGER MY_PE
      COMMON/MPICOM/MY_PE
C
      DOUBLE PRECISION DDOT
      INTEGER ICYCL
      EXTERNAL         DDOT,ICYCL
C
      INTRINSIC MAX0
C
      DATA ERRMSG(1:7)/'BNDCHK '/
C
      IF(.NOT.COMPRESSIBLE)RETURN
      IOPT = 1
C
      WRITE(NOUT,2000)MY_PE
2000  FORMAT(//' CHECKING BOUNDARY FLUXES ON PE #',I4,/' ',19('=')/)
C
C
C     each boundary face has a color, which is stored in IBNDFAC(3,*);
C     the option "-color" assigns a boundary type (see bnd.h) to each color
C
      IFREQ = MAX0( 1 , NBFAC / 20 )
      DO 16 IBFAC = 1 , NBFAC
         IBC = IBNDFAC(3,IBFAC)
C
         IELEM = IBNDFAC(1,IBFAC)
         IVERT = IBNDFAC(2,IBFAC)
C
         CALL DINIT(NOFVAR,ZERO,ZG,1) ! state at mid edge/face
         CALL DINIT(NOFVAR,ZERO,WKSP,1)
         IFACE = ICELFAC(IVERT,IELEM)
C
         DO 12 JVERT = 1,(NOFVERT-1) ! loop over the vertices of the face
            IPOIN = ICELNOD(ICYCL(IVERT+JVERT,NOFVERT),IELEM+NELEM)
            CALL FLUX4(NDIM,ZROE(1,IPOIN),VCN(1,IFACE),FLXN(1,JVERT)) ! compute flux in gridpoint IPOIN
            CALL DAXPY(NOFVAR,ONE,ZROE(1,IPOIN),1,ZG,1)
            CALL DAXPY(NOFVAR,ONE/6.d0,FLXN(1,JVERT),1,WKSP,1) ! only for 2D
            CALL DCOPY(NOFVAR,ZROE(1,IPOIN),1,VCZ((JVERT-1)*NOFVAR+1),1)
   12    CONTINUE ! end loop over the vertices
C
!        CALL R8Mat_Print('General',' ',NOFVAR,NOFVERT,FLXN,
!    +            MAXNOFVAR,'Nodal fluxes ',IFAIL)
C
         CALL DSCAL(NOFVAR,ONE/NDIM,ZG,1)
         CALL DCOPY(NOFVAR,ZG,1,VCZ(NDIM*NOFVAR+1),1) 
C
C     Apply Simpson which IS o.k. in 2D
C
         CALL FLUX4(NDIM,ZG,VCN(1,IFACE),FLXN(1,NOFVERT)) ! compute flux in the center of gravity of the face
         CALL R8Mat_Print('General',' ',NOFVAR,NOFVERT,FLXN,
     +            MAXNOFVAR,'Nodal fluxes ',IFAIL)
         CALL DAXPY(NOFVAR,4.d0/6.d0,FLXN(1,NOFVERT),1,WKSP,1) ! only for 2D
         CALL DCOPY(NOFVAR,WKSP,1,FLXN(1,NOFVERT),1)
         CALL R8Mat_Print('General',' ',NOFVAR,NOFVERT,FLXN,
     +            MAXNOFVAR,'Nodal fluxes ',IFAIL)
         CALL R8Mat_Print('General',' ',NOFVAR,NOFVERT,VCZ,
     +            NOFVAR,'Nodal values ',IFAIL)
C
         CALL GETAVGZ(NDIM,FLXN(1,NOFVERT),VCN(1,IFACE),ZG,ALPHA,IFAIL)
         IF(IFAIL.NE.0)THEN
            WRITE(6,*)'Rootfinding has failed ',IFAIL,' re-trying'
            DO JVERT = 1,(NOFVERT-1) ! loop over the vertices of the face
               IPOIN = ICELNOD(ICYCL(IVERT+JVERT,NOFVERT),IELEM+NELEM)
               CALL DCOPY(NDIM+2,ZROE(1,IPOIN),1,ZG,1)
               CALL GETAVGZ(NDIM,FLXN(1,NOFVERT),VCN(1,IFACE),
     &                      ZG,ALPHA,IFAIL)
               IF(IFAIL.EQ.0)THEN
                  GOTO 10
               ELSE
                  WRITE(6,*)'Changing initial guess did not work ',JVERT
               ENDIF
            ENDDO
         ENDIF
   10 CONTINUE
         CALL FLUX4(NDIM,ZG,VCN(1,IFACE),WKSP) ! compute flux in gridpoint IPOIN
C
      WRITE(6,FMT=*)IELEM
      WRITE(6,FMT=335)(WKSP(IDOF),IDOF=1,NOFVAR)
      WRITE(6,FMT=336)(FLXN(IDOF,NOFVERT),IDOF=1,NOFVAR)
      WRITE(6,FMT=*)(ZG(IDOF),IDOF=1,NOFVAR)
         CALL R8Mat_Print('General',' ',NOFVAR,NOFVERT,VCZ,
     +            NOFVAR,'Nodal values ',IFAIL)
   16 CONTINUE
C
      RETURN
C
C     I/O FORMATS
C
  330 FORMAT(10X,'MASS FLUX THROUGH BOUNDARY ',I2,6(1X,E12.6))
  335 FORMAT(10X,'MASS FLUX AT THE AVERAGED STATE ',6(1X,E12.6))
  336 FORMAT(10X,'MASS FLUX THROUGH THE BOUNDARY  ',6(1X,E12.6))
C
      END
C
      SUBROUTINE GETAVGZ(NDIM,FLXN,VCN,ZAVG,ALPHA,IFAIL)
      IMPLICIT NONE
      INCLUDE 'paramt.h'
      INCLUDE 'constants.h'
      INTEGER NDIM,IFAIL
      DOUBLE PRECISION VCN(NDIM),ZAVG(*),FLXN(*),ALPHA(*)
      INTEGER NITER,ITER
C
C     Input:
C     NDIM is the space dimensionality; unchanged on exit
C     FLXN is the flux through the face computed using Simpson's rule; unchanged on exit
C     VCN  stores the NDIM components of the normal to the face; unchanged on exit
C     ZAVG stores the NDIM+2 components of the initial guess for the averaged state
C
C     Output
C     ZAVG return the averaged state
C
C
      DOUBLE PRECISION DFDZ(MAX_NOFVAR_SQR) ,Z(MAXNOFVAR),DZ(MAXNOFVAR)
      DOUBLE PRECISION RHS(MAXNOFVAR)
      DOUBLE PRECISION S,T
      DOUBLE PRECISION DNRM2
      NITER = 200
C
C     F* = FLXN is the flux through the boundary
C     we whant to find the root of G(z) = F(x) - F*
C     where F(x) is calculated by calling subroutine FLUX4  
C     G(z) = G(z0) + dGdz(z0) * (z-z0)
C     setting G(z) = 0 we get:
C     dGdz(z0) * (z0-z) = G(z0) i.e. Ax = b
C     solving for x
C     we obtain z = z0 - x
C
      DO ITER = 1,NITER
         CALL GETDFDZ(NDIM,ZAVG,VCN,DFDZ,NDIM+2)
         CALL FLUX4(NDIM,ZAVG,VCN,RHS) ! compute flux in z0 and store in rhs
         CALL DAXPY(NDIM+2,MONE,FLXN,1,RHS,1) ! compute rhs = F(z0) - F* = G(z0)
         CALL LUDECO(DFDZ,NDIM+2) ! LU factor
         CALL LUSOLV(DFDZ,RHS,DZ,NDIM+2) ! solve Ax=b and store x in dz
         S = DNRM2(NDIM+2,DZ,1) 
         T = DNRM2(NDIM+2,RHS,1) 
         CALL DAXPY(NDIM+2,MONE,DZ,1,ZAVG,1) ! z0 := z0 - x
         WRITE(6,*)ITER,S,T
         IF(T.LE.1.D-10)THEN
            IFAIL = 0
            RETURN
         ENDIF
      ENDDO
      IFAIL = 1
      RETURN
      END
C
      SUBROUTINE GETDFDZ(NDIM,ZROE,VCN,DFDZ,NOFVAR)
      IMPLICIT NONE
      INCLUDE 'paramt.h'
      INCLUDE 'constants.h'
      INCLUDE 'pfcgas.com'
      INTEGER NDIM,NOFVAR
      DOUBLE PRECISION VCN(NDIM),ZROE(*),DFDZ(NOFVAR,*)
      DOUBLE PRECISION DOTM
      DOTM = ZROE(3)*VCN(1)+ZROE(4)*VCN(2)
      IF(NDIM.EQ.3)DOTM = DOTM + ZROE(5)*VCN(3)
C
      DFDZ(1,1) = DOTM
      DFDZ(1,2) = ZERO
      DFDZ(1,3) = ZROE(1)*VCN(1)
      DFDZ(1,4) = ZROE(1)*VCN(2)
      IF(NDIM.EQ.3)DFDZ(1,5) = ZROE(1)*VCN(3)
C
      DFDZ(2,1) = ZERO
      DFDZ(2,2) = DOTM
      DFDZ(2,3) = ZROE(2)*VCN(1)
      DFDZ(2,4) = ZROE(2)*VCN(2)
      IF(NDIM.EQ.3)DFDZ(2,5) = ZROE(2)*VCN(3)
C
      DFDZ(3,1) = GM1OG*DFDZ(2,3)
      DFDZ(4,1) = GM1OG*DFDZ(2,4)
      IF(NDIM.EQ.3)DFDZ(5,1) = GM1OG*DFDZ(2,5)
C
      DFDZ(3,2) = GM1OG*DFDZ(1,3)
      DFDZ(4,2) = GM1OG*DFDZ(1,4)
      IF(NDIM.EQ.3)DFDZ(5,2) = GM1OG*DFDZ(1,5)
C
      DFDZ(3,3) =-GM1OG*ZROE(3)*VCN(1)+TWO*ZROE(3)*VCN(1)
      DFDZ(3,4) =-GM1OG*ZROE(4)*VCN(1)+    ZROE(3)*VCN(2)
      IF(NDIM.EQ.3)
     &DFDZ(3,5) =-GM1OG*ZROE(5)*VCN(1)+    ZROE(3)*VCN(3)
C
      DFDZ(4,3) =-GM1OG*ZROE(3)*VCN(2)+    ZROE(4)*VCN(1)
      DFDZ(4,4) =-GM1OG*ZROE(4)*VCN(2)+TWO*ZROE(4)*VCN(2)
C
      IF(NDIM.EQ.2)RETURN
C
      DFDZ(4,5) =-GM1OG*ZROE(5)*VCN(2)+    ZROE(4)*VCN(3)
C
      DFDZ(5,3) =-GM1OG*ZROE(3)*VCN(3)+    ZROE(5)*VCN(1)
      DFDZ(5,4) =-GM1OG*ZROE(4)*VCN(3)+    ZROE(5)*VCN(2)
      DFDZ(5,5) =-GM1OG*ZROE(5)*VCN(3)+TWO*ZROE(5)*VCN(3)
C
      RETURN
      END
