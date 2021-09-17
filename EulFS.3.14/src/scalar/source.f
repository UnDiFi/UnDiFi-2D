      SUBROUTINE FSOU(NDIM,NOFVERT,TAU,K,NODRES,VCN,VCP,VOLUME,FNCTN)
C
      IMPLICIT NONE
C
      INCLUDE 'paramt.h'
      INCLUDE 'constants.h'
C
C     Integral of the source term using an SUPG approach
C     \int f ( w + TAU a \cdot \nabla w )
C
C     ( w + TAU a \cdot \nabla w ) is the SUPG test function
C       w is the shape function i.e. the tent function
C
C     Two integrals have to be evaluated:
C
C     1) \int f w_i dV
C     2) TAU k_i / V \int f dV
C
C
C     .. Parameters ..
      INTEGER IPMAX
      PARAMETER (IPMAX=7)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION TAU
      INTEGER NDIM,NOFVERT
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION K(NOFVERT),NODRES(NOFVERT),VCN(NDIM,NOFVERT),
     +VCP(NDIM,NOFVERT)
C     ..
C     .. Function Arguments ..
      DOUBLE PRECISION FNCTN
      EXTERNAL FNCTN
C     ..
C     .. Arrays in Common ..
      DOUBLE PRECISION GCOOR(IPMAX,3),GWGHT(IPMAX)
C     ..
C     .. Scalars in Common ..
      INTEGER NPTS
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION SHAPE,UPWD,V2,W,exact,VOLUME
      DOUBLE PRECISION CNST,VALUE,XGAUSSP,YGAUSSP,ZGAUSSP
      INTEGER INODE,IPT,IVERT
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION V1(VMAX)
C     ..
C     .. Common blocks ..
      COMMON /GAUSSBLOCK/GCOOR,GWGHT,NPTS
C     ..
C
      CNST = ONE/ (NDIM*VOLUME)
C
C     Compute the coordinates of each Gauss point
C
      DO 2 IVERT = 1,NOFVERT
          V1(IVERT) = ZERO
    2 CONTINUE
      V2 = ZERO
C
      DO 1 IPT = 1,NPTS
C
          XGAUSSP = ZERO
          YGAUSSP = ZERO
          ZGAUSSP = ZERO
C
          DO 3 IVERT = 1,NOFVERT
C
              W = GCOOR(IPT,IVERT)
              XGAUSSP = XGAUSSP + W*VCP(1,IVERT)
              YGAUSSP = YGAUSSP + W*VCP(2,IVERT)
              IF (NDIM.EQ.3) ZGAUSSP = ZGAUSSP + W*VCP(3,IVERT)
C
    3     CONTINUE
C
C     Evaluate the function in the Gauss point
C
          VALUE = FNCTN(XGAUSSP,YGAUSSP,ZGAUSSP)
C
C     Evaluate the integral \int f dV
C
          V2 = V2 + TWO*GWGHT(IPT)*VALUE
C
C     Evaluate the shape function of all vertices in the Gauss point
C
          DO 5 IVERT = 1,NOFVERT
C
              SHAPE = ONE + CNST* (VCN(1,IVERT)* (XGAUSSP-VCP(1,IVERT))+
     +                VCN(2,IVERT)* (YGAUSSP-VCP(2,IVERT))+
     +                VCN(3,IVERT)* (ZGAUSSP-VCP(3,IVERT)))
C
              V1(IVERT) = V1(IVERT) + TWO*GWGHT(IPT)*VALUE*SHAPE
C
    5     CONTINUE
C
    1 CONTINUE
C
C     exact = 1./3.*(fnctn( vcp(1,1),vcp(2,1),1. ) +
C    +             + fnctn( vcp(1,2),vcp(2,2),1. )
C    +             + fnctn( vcp(1,3),vcp(2,3),1. ) )
C     write(6,*)exact,v2
C     pause
C
      DO 7 IVERT = 1,NOFVERT
C
          V1(IVERT) = V1(IVERT) * VOLUME 
C
          UPWD = K(IVERT)*TAU*V2
C         write(6,*)ivert,k(ivert),v2,v1(ivert)
C         write(6,*)ivert,k(ivert),upwd/v2/volume,v1(ivert)/v2/volume 
C         if(v1(ivert).lt.0)pause
C
          NODRES(IVERT) = NODRES(IVERT) + V1(IVERT) + UPWD
C
    7 CONTINUE
      RETURN

      END

      double precision function myfun(x,y,z)

      double precision x,y,z

      myfun = 4.d0*x-2.d0*y

      return

      end  
