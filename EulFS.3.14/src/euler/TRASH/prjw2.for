C
      SUBROUTINE PRJW2(ZROE,VN,PPOSZ,PNEGZ,PPOSU,PNEGU,NDIM,NOFVAR)
C
C     Compute projectors for inviscid wall boundary conditions
C
      IMPLICIT NONE
C
      INCLUDE 'paramt.h'
      INCLUDE 'constants'
      INCLUDE 'chorin.com'
      INCLUDE 'three'
C
C#define DEBUG
C
#ifdef DEBUG
      INTEGER ifail,lwork,ipiv(25)
      LOGICAL WARNA
      parameter(lwork=5)
      DOUBLE PRECISION temp(25),work(lwork)
#endif
C
C     .. Scalar Arguments ..
C
      INTEGER NOFVAR,NDIM
C
C     .. Array Arguments ..
C
C     VNOR stores by columns an orthonormal base (n,s,t) where n
C          is multiplied by the surface of the face
C     ZROE stores the parameter vector of the face vertices
C
      DOUBLE PRECISION VN(NDIM),PPOSU(NOFVAR,NOFVAR),
     +PNEGU(NOFVAR,NOFVAR),PPOSZ(NOFVAR,NOFVAR),PNEGZ(NOFVAR,NOFVAR),
     +ZROE(NOFVAR)
C
C     .. Local Scalars ..
C
      INTEGER i,j,ielem
      DOUBLE PRECISION COEFF1,COEFF2,COEFF3,AREA,AREA2
      DOUBLE PRECISION COEFF4,COEFF5,COEFF6,DENOM
      DOUBLE PRECISION C,UDOTN,UNPOS,AA,BB
C
C     .. Local Arrays ..
C
      DOUBLE PRECISION WORK1(25)
C
      DOUBLE PRECISION DNRM2
      EXTERNAL         DNRM2
C
C     .. Data Statements ..
C
C     .. Executable Statements ..
C
C     .. Solid Wall Boundary Conditions
C
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION DUMMY(25),KMAT(25),VLEFT(5,5),VRIGHT(5,5),
     +                 WNEG(5),WPOS(5),WR(5)
      INTEGER IPIV(5)
C     ..
C     .. External Subroutines ..
      EXTERNAL DGETRF,DGETRS,MATSPLITVIII
C     ..
      IELEM = -1
C     write(6,*)'n = ',(vecn(i),i=1,ndim)
C     write(6,*)'Z = ',(ZROE(i),i=1,nofvar)
C
C
C
      AREA = VN(1)*VN(1)+ VN(2)*VN(2)
      IF(NDIM.EQ.3)AREA=AREA+VN(3)*VN(3)
      AREA2 = AREA
      AREA = SQRT(AREA2)
      DENOM=ONE/AREA
C
      UDOTN = ZAVG(2)*VN(1)+ZAVG(3)*VN(2)
      IF(NDIM.EQ.3)UDOTN=UDOTN+ZAVG(4)*VN(3)
      UDOTN=UDOTN/AREA
C
      C = SQRT(BETA+UDOTN*UDOTN)
C
      UNPOS = 0.5d0*(ONE+ABS(UDOTN))
C
      AA = 2.d0*UNPOS*UDOTN*UDOTN + BETA * (UDOTN+C)
      DENOM = UNPOS*(UDOTN+C)*C
C
C     bb gia` diviso U.n 
C
      BB = BETA * (2.d0*UNPOS-(UDOTN+C))
      AA = AA/DENOM
      BB = BB/DENOM
C
C     .. Z^{*} = sqrt{rho} [ 1,H,-u_n,u_s,u_t]
C
C     COEFF1 = - VN(1)*VN(1) + VS(1)*VS(1) + VT(1)*VT(1)
C     COEFF2 = - VN(2)*VN(2) + VS(2)*VS(2) + VT(2)*VT(2)
C     COEFF3 = - VN(3)*VN(3) + VS(3)*VS(3) + VT(3)*VT(3)
C     COEFF4 = - VN(1)*VN(2) + VS(1)*VS(2) + VT(1)*VT(2)
C     COEFF5 = - VN(1)*VN(3) + VS(1)*VS(3) + VT(1)*VT(3)
C     COEFF6 = - VN(2)*VN(3) + VS(2)*VS(3) + VT(2)*VT(3)
C
      COEFF1 = (ONE - AA)- BB*VN(1)*VN(1)/AREA2
      COEFF2 = (ONE - AA)- BB*VN(2)*VN(2)/AREA2
      COEFF4 = - BB*VN(1)*VN(2)/AREA2
      COEFF5 = - BB*VN(1)*VN(3)/AREA2
      IF(NDIM.EQ.3)THEN
      COEFF3 = (ONE - AA)- BB*VN(3)*VN(3)/AREA2
      COEFF6 = - BB*VN(2)*VN(3)/AREA2
      ENDIF
C
C     An other approach
C
C     we take Z^{*} = sqrt{rho}[1,H,0,u_t,u_s]
C
C     TEMP = UAVG(3)*VN(1)+UAVG(4)*VN(2)+UAVG(5)*VN(3)
C     TEMP = ONE/(QINV*QINV) - TEMP*TEMP
C     TEMP = ONE/(QINV*DSQRT(TEMP))
C
C     COEFF1 = VS(1)*VS(1) + VT(1)*VT(1)
C     COEFF2 = VS(2)*VS(2) + VT(2)*VT(2)
C     COEFF3 = VS(3)*VS(3) + VT(3)*VT(3)
C     COEFF4 = VS(1)*VS(2) + VT(1)*VT(2)
C     COEFF5 = VS(1)*VS(3) + VT(1)*VT(3)
C     COEFF6 = VS(2)*VS(3) + VT(2)*VT(3)
C
C     COEFF1 = TEMP*COEFF1
C     COEFF2 = TEMP*COEFF2
C     COEFF3 = TEMP*COEFF3
C     COEFF4 = TEMP*COEFF4
C     COEFF5 = TEMP*COEFF5
C     COEFF6 = TEMP*COEFF6
C
C
      PNEGZ(1,1) = ONE
      PNEGZ(1,2) = ZERO
      PNEGZ(1,3) = ZERO
C
      PNEGZ(2,1) = ZERO
      PNEGZ(2,2) = COEFF1
      PNEGZ(2,3) = COEFF4
C
      PNEGZ(3,1) = ZERO
      PNEGZ(3,2) = COEFF4
      PNEGZ(3,3) = COEFF2
C
      IF(NDIM.EQ.3)THEN
C
      PNEGZ(1,4) = ZERO
      PNEGZ(2,4) = COEFF5
      PNEGZ(3,4) = COEFF6
C
      PNEGZ(4,1) = ZERO
      PNEGZ(4,2) = COEFF5
      PNEGZ(4,3) = COEFF6
      PNEGZ(4,4) = COEFF3
C
      ENDIF
C
C
      DO 3 J = 1,NOFVAR
         DO 3 I = 1,NOFVAR
            IF( I .EQ. J )THEN
               PPOSZ(I,J) = ONE - PNEGZ(I,J)
            ELSE
               PPOSZ(I,J) = - PNEGZ(I,J)
            ENDIF
    3 CONTINUE
C
C     AUX(1) = AA * ZROE(2) + BB * VN(1) * AREAI
C     AUX(2) = AA * ZROE(3) + BB * VN(2) * AREAI
C     AUX(1) = AA * ZROE(4) + BB * VN(3) * AREAI
C
C     CALL MATSPLITVIII(IELEM,NDIM,NOFVAR,VN,DUMMY,NOFVAR,KMAT,PPOSU,
C    +                  PNEGU,VLEFT,VRIGHT,NOFVAR,WR,WPOS,WNEG,.TRUE.)
C     CALL DGEMM
C
      RETURN
      END
