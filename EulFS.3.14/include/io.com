      INTEGER      NOUT,IWUNIT,ITIM1,IHST1,IHST2,IHST3,IPROBE,
     &IOALE,IHST4,IWEFLX
      COMMON/IODEV/NOUT,IWUNIT,ITIM1,IHST1,IHST2,IHST3,IPROBE,
     &IOALE,IHST4,IWEFLX
C
C	NOUT   is the OUTPUT device number where processor related
C              info are written
C	IWUNIT is the OUTPUT device number where global info are
C              written
C
C       REM: the following units are set in blockdata
C
C       ITIM1 = 4  is the unit for timing
C       IHST1 = 1  is the unit for convergence history
C       IHST2 = 2  is the unit for convergence history
C       IHST3 = 3  is the unit for convergence history
C       IHST4 = 116  is the unit for the volume integral 
C       IOALE  = 115 is dumps info on grid velocity
C       IWEFLX  = 117 is dumps info on current through patches
C
