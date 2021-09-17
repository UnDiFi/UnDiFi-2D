      SUBROUTINE PARM2MERKLE(ZROE,DUDZ,NOFVAR,NDIM)
C
C     $Id: parm2merkle.f,v 1.3 2010/11/13 11:06:31 abonfi Exp abonfi $
C     $Header: /msa20/home/abonfi/CFD_codes/EulFS.0.15.2/src/euler/RCS/parm2merkle.f,v 1.3 2010/11/13 11:06:31 abonfi Exp abonfi $
C
      IMPLICIT NONE
C
C     transformation matrix from
C     parameter vector to (p,u,v,w,T)
C
C
C
C
C     Assembles the dUdZ matrix ...
C
C     .. Parameters ..
      INCLUDE 'paramt.h'
      INCLUDE 'constants.h'
C
C     .. Common ..
      INCLUDE 'merkle.com'
      INCLUDE 'stream.com'
C     ..
C     .. Scalar Arguments ..
      INTEGER NDIM,NOFVAR
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION DUDZ(NOFVAR,NOFVAR),ZROE(NOFVAR)
C     ..
C     .. Local Scalars ..
C     ..
      DOUBLE PRECISION RHOINV,HELP,XI
C
      RHOINV = ONE/(ZROE(1)*ZROE(1))
C
      DUDZ(1,1) = GM1OG*ZROE(2) ! dp/dz1
      DUDZ(1,2) = GM1OG*ZROE(1) ! dp/dz2
      DUDZ(1,3) =-GM1OG*ZROE(3) ! dp/dz3
      DUDZ(1,4) =-GM1OG*ZROE(4) ! dp/dz4
C
      DUDZ(2,1) =-ZROE(3)*RHOINV  ! du/dz1
      DUDZ(2,2) = ZERO            ! du/dz2
      DUDZ(2,3) = ONE/ZROE(1)     ! du/dz3
      DUDZ(2,4) = ZERO            ! du/dz4
C
      DUDZ(3,1) =-ZROE(3)*RHOINV  ! dv/dz1
      DUDZ(3,2) = ZERO
      DUDZ(3,3) = ZERO
      DUDZ(3,4) = ZROE(1)
C
      HELP = ZROE(3)*ZROE(3)+ZROE(4)*ZROE(4) 
      IF(NDIM.EQ.3)HELP = HELP + ZROE(5)*ZROE(5)
      HELP = HELP/ZROE(1)
!     XI = GM1*M_INFTY*M_INFTY 
      XI = GM1OG/RSTAR
C
      DUDZ(4,1) = XI*RHOINV*(-ZROE(2)+HELP)        ! dT/dz1
      DUDZ(4,2) =-XI*ZROE(3)*RHOINV
      DUDZ(4,3) =-TWO*XI*ZROE(3)*ZROE(3)*RHOINV/ZROE(1)
      DUDZ(4,4) =-TWO*XI*ZROE(4)*ZROE(4)*RHOINV/ZROE(1)
C
      IF (NDIM.EQ.2) RETURN
C
      STOP 'Transformation un-implemented in 3D'
C
      DUDZ(1,5) =-GM1OG*ZROE(5)
      DUDZ(2,5) =-XI*ZROE(5)*RHOINV
      DUDZ(3,5) = ZERO
      DUDZ(4,5) = ZERO
C
      DUDZ(5,1) = ZROE(5)
      DUDZ(5,2) = ZERO
      DUDZ(5,3) = ZERO
      DUDZ(5,4) = ZERO
      DUDZ(5,5) = ONE/ZROE(1)
C
C
      RETURN

      END
