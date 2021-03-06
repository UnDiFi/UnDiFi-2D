      INTEGER MAX_ANGLES,LDW
      PARAMETER (MAX_ANGLES = 3, LDW = MAX_ANGLES*3+2)

      REAL*8 ZAVG(MAXNOFVAR*MAXTIMLEVS), UAVG(MAXNOFVAR),
     1  DivFlux(MAXNOFVAR) ,R_SPEED(3,LDW), 
     2  GRAD_PRIM(MAXNOFEQN,3), GRAD_PARM(MAXNOFVAR,3),
     3  GRAD_CHAR(LDW,3),
     4  ABAR, ASQR , KINETIC , MACH, MACHSQR, QINV, BAVG(3)
C
      COMMON/THREE/ ZAVG, UAVG, DivFlux , R_SPEED ,
     .  GRAD_PRIM, GRAD_PARM, GRAD_CHAR, 
     .  ABAR , ASQR , KINETIC, MACH, MACHSQR, QINV, BAVG
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
