      DOUBLE PRECISION RE,REINV,PRANDTL,TWALL
      INTEGER IADIA,IADIA_REPLACE_EQN
      COMMON /VISCO_R/ RE,REINV,PRANDTL,TWALL
      COMMON /VISCO_I/ IADIA,IADIA_REPLACE_EQN
C
C     Applicability KAN = +1,+2,+4
C
C     RE      is the Reynolds number
C     REINV   is the inverse of the Reynolds number
C     PRANDTL is the Prandtl number
C     TWALL   is the non-dimensional wall temperature 
C             (meaningful for IADIA <> 0) 
C     IADIA   set to 0 means adiabatic wall else isothermal
C     IADIA_REPLACE_EQN  set to IADIA_REPLACE_CONT or IADIA_REPLACE_ENER
C                         these are defined in visco.h
