      DOUBLE PRECISION AMPSQR,MERKLE_CUTOFF
      COMMON/MERKLE/AMPSQR,MERKLE_CUTOFF
C
C     This is Merkle's parameter M_p for preconditiong, see e.g.
C     the VKI LS 1999-03
C
C     RSTAR is a constant in the non-dimensional equation of state:
C     p = \rho RSTAR T
C
C     T changes depending on the choice of reference variables,
C     i.e. depending on the flag -nondimensionalisation
C     for external flows RSTAR = 1./(\gamma * M_{\infty}**2)
C     for internal flows RSTAR = 1.
C
