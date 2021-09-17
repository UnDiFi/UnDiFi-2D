      CHARACTER BAKFILE*512,VISCTFILE*255,FOLDFILE*255
      INTEGER NWFAC,NBODY6
C
C     A nasty way to pass these variables to the subroutine
C     myTS()
C     rem that on a Compaq I had to split character and integer
C     data on two different common blocks
C
C
      COMMON/COMZZC/BAKFILE,VISCTFILE,FOLDFILE
      COMMON/COMZZI/NWFAC,NBODY6
