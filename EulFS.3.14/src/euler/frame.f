C
      SUBROUTINE StreamAlignedFrame(NDIM)
C
      IMPLICIT NONE
C
      INCLUDE 'paramt.h'
      INCLUDE 'constants.h'
C
      INCLUDE 'three.com'
      INCLUDE 'frame.com'
C
      INTEGER I,NDIM

      DOUBLE PRECISION SING,COSG,TEMP2
C
      EXTERNAL CROSS_PROD,DCOPY
C
      INTEGER  ISDMIN,JCYCL
      EXTERNAL ISDMIN,JCYCL
      DOUBLE PRECISION DNRM2
      EXTERNAL         DNRM2
C
      GOTO(20,30),NDIM-1
C
   20 CONTINUE
C
      RotationMatrix(1,1) = UAVG(3)*QINV
      RotationMatrix(2,1) = UAVG(4)*QINV
      RotationMatrix(3,1) = ZERO
C
      RotationMatrix(1,2) =-UAVG(4)*QINV
      RotationMatrix(2,2) = UAVG(3)*QINV
      RotationMatrix(3,2) = ZERO
C
      CALL CROSS_PROD( RotationMatrix(1,1) , RotationMatrix(1,2) ,
     &                 RotationMatrix(1,3) )
C
      RETURN
C
   30 CONTINUE

      RotationMatrix(1,1) = UAVG(3) * QINV
      RotationMatrix(2,1) = UAVG(4) * QINV
      RotationMatrix(3,1) = UAVG(5) * QINV
C
      TEMP2 = ONE / SQRT( UAVG(3)*UAVG(3)+UAVG(4)*UAVG(4) )
C
      COSG = ONE
      SING = ZERO
C
C
      RotationMatrix(1,2) =(-UAVG(4)*COSG-UAVG(3)*UAVG(5)*SING*
     1				QINV)*TEMP2
      RotationMatrix(2,2) =( UAVG(3)*COSG-UAVG(4)*UAVG(5)*SING*
     1				QINV)*TEMP2
      RotationMatrix(3,2) =  SING*QINV/TEMP2
C
      RotationMatrix(1,3) =( UAVG(4)*SING-UAVG(3)*UAVG(5)*COSG*
     1				QINV)*TEMP2
      RotationMatrix(2,3) =(-UAVG(3)*SING-UAVG(4)*UAVG(5)*COSG*
     1				QINV)*TEMP2
      RotationMatrix(3,3) =  COSG*QINV/TEMP2
C
      RETURN
      END
