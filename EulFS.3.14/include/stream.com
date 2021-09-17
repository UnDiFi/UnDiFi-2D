C
C     This common block contains data relevant to
C     the non-dimensionalisation used; the latter is mainly 
C     reflected in the non-dimensional form of
C     Sutherland's law 
C
C
      DOUBLE PRECISION M_infty,U_infty(MAXNOFVAR),POUTLET,TREF,
     &UREF,PREF,RREF,LREF,FLOWDIR(3),RSTAR
      COMMON/STREAM/   M_infty,U_infty,POUTLET,TREF,
     &UREF,PREF,RREF,LREF,FLOWDIR,RSTAR
C
C     M_infty is the freestream Mach number
C     U_infty stores the NOFVAR freestream values (parameter vector)
C     VISCF   is the viscous force acting on the body
C     PRESF   is the pressure force acting on the body
C
C     TREF,UREF,PREF,RREF,LREF are the REFERENCE
C
C     temperature, velocity, pressure, density, lenght scale
C
C     see setibc.F to understand how these are computed
C
