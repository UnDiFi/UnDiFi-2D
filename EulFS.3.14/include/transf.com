      DOUBLE PRECISION DUDV(MAX_NOFVAR_SQR*MAXTIMLEVS),DZDU(192),
     1                 DVDZ(MAX_NOFVAR_SQR*MAXTIMLEVS)
      COMMON/COMTRX/DUDV,DVDZ,DZDU
C
C     transformation matrices:
C     DUDV from characteristic/symmetrizing to conserved
C     DZDU from conserved to parameter vector
C     DVDZ from parameter vector to characteristic/symmetrizing
C
