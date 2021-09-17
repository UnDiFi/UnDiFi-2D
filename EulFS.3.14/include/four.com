C     NSP: SPECIES NUMBER
C     NR: REACTIONS NUMBER
C     IE: ENERGY INDEX
C     IX,IY,IZ: MOMENTUM INDEX
      DOUBLE PRECISION DR(NSP), DE, DM(3), CHI(NSP), ALPHA(NSP)
      DOUBLE PRECISION KAPPA, DENS, DENSINV, SQRTR
C
      COMMON/FOUR/ DR, DE, DM, CHI, ALPHA, KAPPA, DENS, DENSINV, SQRTR
C
C ZAVG  is the averaged state in the cartesian ref. frame
C UAVG  is the averaged state in the cartesian ref. frame
C DivFlux is the flux divergence (used when ICHECK <> 0)
C R_SPEED(*,KWAVE) is the advection speed of the KWAVE-th wave
C GRAD_PRIM gradient of the primitive variables
C GRAD_PARM gradient of Roe's parameter vector
C GRAD_CHAR gradient of the characteristic variables
C ABAR averaged sound speed
C ASQR averaged squared sound speed
C KINETIC averaged kinetic energy
C MACH  averaged Mach number
C MACHSQR averaged squared Mach number
C QINV  1./ averaged velocity magnitude
C BAVG is the cell averaged grid velocity
C
