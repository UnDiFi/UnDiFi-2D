      SUBROUTINE WTINDX(TAUW,TINDX,VOL,NWFAC,IBNDFAC,NBODY6,
     +NDIM,NOFVERT,NOFVAR,COMPRESSIBLE,FILENAME)
C
C     Write turbulent index data to a file
C
      IMPLICIT NONE
C
      INCLUDE 'paramt.h'
      INCLUDE 'bnd.h'
      INCLUDE 'bnd.com'
      INCLUDE 'nloc.com'
      INCLUDE 'three.com'
      INCLUDE 'stream.com'
      INCLUDE 'turb.com'
      INCLUDE 'visco.com'
C     .. Parameters ..
      DOUBLE PRECISION ZERO,ONE
      PARAMETER (ZERO=0.D0,ONE=1.D0)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION VOLUME,RHOW
      LOGICAL COMPRESSIBLE
      INTEGER NDIM,NOFVAR,NOFVERT,IOFF
C     ..
      DOUBLE PRECISION VCN(12),VCZ(20),VCP(12)
      INTEGER ICF(4),ICN(4),IDX_F(16),IDX_P(12),IDX_V(20)
C     ..
C     .. Arrays in Common ..
      DOUBLE PRECISION DSTAK(1)
C     ..
C     .. Local Scalars ..
C     ..
C     .. Local Arrays ..
      INTEGER ISTAK(1),IARRAY(4)
      DOUBLE PRECISION DENS(3)
C     ..
C     .. External Subroutines ..
      EXTERNAL DGTHR,DSCAL,R8Mat_Print
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,SIGN
C     ..
C     .. Common blocks ..
      COMMON /CSTAK/DSTAK
C     ..
C     .. Statement Functions ..
      INTEGER IDX
C     ..
C     .. Equivalences ..
      EQUIVALENCE (ISTAK(1),DSTAK(1))
C
C
C     .. Scalar Arguments ..
      INTEGER NBODY6,NWFAC
      CHARACTER FILENAME* (*)
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION TAUW(NWFAC),TINDX(NWFAC),VOL(*)
      INTEGER IBNDFAC(3,*)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION ALPHA,DVTDN,UTAU
      INTEGER I,IC,IELEM,IFACE,JFACE,IFRST,ILAST,IVERT,IXDRS
      INTEGER IDIM,IFAIL,IPOIN,J,IV,N,NZD
C     ..
C     .. External Functions ..
      INTEGER INITXDR,ICYCL
      DOUBLE PRECISION DNRM2,DDOT
      EXTERNAL INITXDR,DNRM2,DDOT,ICYCL
      INTEGER IXDRCLOSE,IXDRIMAT,IXDRINT
      EXTERNAL IXDRCLOSE,IXDRIMAT,IXDRINT
#ifdef PARCH_t3d
      INTEGER IXDRRMAT64
      EXTERNAL IXDRRMAT64
#else
      INTEGER IXDRDMAT
      EXTERNAL IXDRDMAT
#endif
C     ..
      DATA IARRAY/4*0/
C     ..
C     ..
C     .. Statement Function definitions ..
C
      IDX(I,J,N) = (I-1)*N + J
C
C
      IXDRS = INITXDR(FILENAME,'w',.FALSE.)
C
C
      IC = 0
      DO 10 I = 1,NBODY6
C
          IFRST = IBGN(I)
          ILAST = IEND(I)
C
          DO 10 IFACE = IFRST,ILAST
              IELEM = IBNDFAC(1,IFACE)
              IVERT = IBNDFAC(2,IFACE)
              IC = IC + 1
C
C
C     cell to node and cell to faces pointer
C
      IOFF = IDX(IELEM,1,NOFVERT)-1
      CALL ICOPY(NOFVERT,ISTAK(LCELNOD+IOFF),1,ICN,1)
      CALL ICOPY(NOFVERT,ISTAK(LCELFAC+IOFF),1,ICF,1)
C
      NZD = 0
      DO 3 J = 1,NOFVERT
          IPOIN = ICN(J )
          JFACE = ABS( ICF(J) )
          DO 3 IDIM = 1,NDIM
              NZD = NZD + 1
              IDX_P(NZD) = IDX(IPOIN,IDIM,NDIM)
              IDX_F(NZD) = IDX(JFACE,IDIM,NDIM)
    3     CONTINUE
C
C
C     retrieve
C     1) nodal values (turbulence)
C     3) normal components
C     from global arrays
C
      CALL DGTHR(NOFVERT,DSTAK(LTURB),VCZ,ICN)
      CALL DGTHR(NZD,DSTAK(LCORG),VCP,IDX_P)
      CALL DGTHR(NZD,DSTAK(LFACNOR),VCN,IDX_F)
C
      VOLUME = VOL(IELEM)
C
      DO 12 IV = 1,NOFVERT
          ALPHA = SIGN(1,ICF(IV))
          CALL DSCAL(NDIM,ALPHA,VCN(IDX(IV,1,NDIM)),1)
   12 CONTINUE
C
      ALPHA = DNRM2(NDIM,VCN(IDX(IVERT,1,NDIM)),1)
C
C     compute the gradient of the eddy viscosity
C
      CALL LINEARIZE(IELEM,.FALSE.,VCN,NDIM,NOFVERT,VCZ,1,VOLUME)
C
C     write(6,*)iface,ibndfac(3,iface),(icn(iv),iv=1,nofvert) ,
C    +zavg(1)
C     CALL R8Mat_Print('General',' ',NDIM,NOFVERT,VCP,NDIM,
C    +         'Nodal coord. ',IFAIL)
C     CALL R8Mat_Print('General',' ',NOFVERT,1,VCZ,NOFVERT,
C    +         'Nodal values ',IFAIL)
C     CALL R8Mat_Print('General',' ',5,NDIM,GRAD_PARM,NMAX,
C    +         'Gradient ',IFAIL)
C     pause
C
C     compute the component normal to the wall of the eddy viscosity
C
      DVTDN = DDOT(NDIM,VCN(IDX(IVERT,1,NDIM)),1,GRAD_PARM(1,1),NMAX) 
     +/ ALPHA
C
C     alpha = dnrm2(ndim,grad_parm,NMAX)
C
C     compute an averaged face density 
C
C     IF( COMPRESSIBLE )THEN
C     DO 9 IV = 1,NOFVERT-1
C  IPOIN = ICN( ICYCL(IVERT+IV,NOFVERT) )
C         IDX_V(IV) = IDX(IPOIN,1,NOFVAR)
C   9 CONTINUE
C     CALL DGTHR(NOFVERT-1,DSTAK(LZROE),DENS,IDX_V)
C     RHOW = ZERO
C     DO 19 IV = 1,NOFVERT-1
C         RHOW = RHOW + DENS(IV)
C  19 CONTINUE
C     RHOW = RHOW / (NOFVERT-1)
C     ELSE
C     RHOW = ONE
C     ENDIF
C
C     UTAU = SQRT( ABS(TAUW(IC)) / RHOW )
C
C     TINDX(IC) = DVTDN * SQRT(REINV) / (UTAU * TK)
      TINDX(IC) = DVTDN
C     TINDX(IC) = alpha
C
   10 CONTINUE

      IF (IC.NE.NWFAC) THEN
          WRITE (6,FMT=*) 'Ic NFAC',IC,NWFAC
          STOP

      ENDIF
#ifdef PARCH_t3d
      IFAIL = IXDRDMAT(IXDRS,NWFAC,TINDX)
#else
      IFAIL = IXDRDMAT(IXDRS,NWFAC,TINDX)
#endif
C
      IFAIL = IXDRCLOSE(IXDRS)
      WRITE(6,*)'       Turbulent index written to ',FILENAME 
C
      RETURN

      END
