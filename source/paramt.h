!     EPS         : displacement between the two faces of one internal boundary
!     GA          : scecific heat ratio
!     GM1         : GM1=GA-1
!     SNDMIN      : Max normalized distance of a phantom point
!     DXCELL      : length of a cell side
!     CFL         : stability coefficient for the shock motion
!     IBAK        : save solution every IBAK iterations
!     NADDHOLES   : Number of additional hole points
!     CADDhole    : Coordinates of additional hole point
!     NprdBnd     : Number of periodic boundaries
!     prdBndclr   : array of the pair of colors of each periodic boundaries
!     FLT_Dspeed  : filter on discontinuity speeds
!     IMTF        : iteration of mesh topology freezing
!
!     NSHMAX, NPSHMAX, NESHMAX, NDIM, NDOF, NSPMAX, NADDHOLESMAX,
!     NprdBndMAX and the ZERO/HALF/ONE/TWO/PI literals now live in
!     mod_constants (ROADMAP.md Phase 2.1, issue #13).
      integer(i4) IBAK, IMTF

!     .. Common area  ..

      real(wp)  EPS,SHRELAX,GA,GM1,SNDMIN,DXCELL,                    &
              CADDhole(NDIM, NADDHOLESMAX),FLT_Dspeed
      integer(i4) NADDHOLES,NprdBnd,prdBndclr(3,NprdBndMAX)
      COMMON/PARAMT/CADDhole,EPS,SNDMIN,DXCELL,GA,SHRELAX,           &
                    GM1,IBAK,NADDHOLES,NprdBnd,prdBndclr,            &
                    FLT_Dspeed,IMTF
