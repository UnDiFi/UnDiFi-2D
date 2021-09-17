!
!     This common block contains data relevant to
!     the non-dimensionalisation used; the latter is mainly 
!     reflected in the non-dimensional form of
!     Sutherland's law 
!     Copy of common stream.com modified to deal with .f90 routines
!
!
      DOUBLE PRECISION TREFP,UREFP,PREFP,RREFP,LREFP,RSTARP(NSP),HREFP
      COMMON/STREAMPLASMA/   TREFP,UREFP,PREFP,RREFP,LREFP,RSTARP,HREFP
!
!
!     TREFP,UREFP,PREFP,RREFP,LREFP,HREFP are the REFERENCE
!
!     temperature, velocity, pressure, density, lenght scale, enthalpy
!
!     RSTARP(i) is the gas constant for the i chemical species
!
!     see setibc.F to understand how these are computed
!
