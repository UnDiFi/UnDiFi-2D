!
!     This common block contains data relevant to
!     the internal non-dimensionalisation (plasma flow);
!
!
      DOUBLE PRECISION ALPHA1(NSP),RSSTAR(NSP),RMIXSTAR 
      COMMON/IOPLASMA/ALPHA1,RSSTAR,RMIXSTAR   
!
!     ALPHA1: vector of the inlet chemical species concentrations (\rho_i1/rho_1)
!     RSSTAR: adimensional gas constant for the single species
!     RMIXSTAR: global adimensional gas constant 
!
!     see setibc.F to understand how these are computed
!
