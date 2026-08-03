!     NSHMAX      . max number od shocks
!     NPSHMAX     : max number of shock points for each shock.
!     NESHMAX     : max number of shock element for each shocks
!     EPS         : displacement between the two faces of one internal boundary
!     GA          : scecific heat ratio
!     GM1         : GM1=GA-1
!     SNDMIN      : Max normalized distance of a phantom point
!     DXCELL      : length of a cell side
!     CFL         : stability coefficient for the shock motion
!     IBAK        : save solution every IBAK iterations
!     NADDHOLESMAX: MAX of additional hole points
!     NADDHOLES   : Number of additional hole points
!     CADDhole    : Coordinates of additional hole point
!     NprdBndMAX  : MAX of periodic boundaries
!     NprdBnd     : Number of periodic boundaries
!     prdBndclr   : array of the pair of colors of each periodic boundaries
!     FLT_Dspeed  : filter on discontinuity speeds
!     IMTF        : iteration of mesh topology freezing

!     .. Parameters ..

      REAL*8  ZERO,HALF,ONE,TWO,PI
      INTEGER*4 NSHMAX, NPSHMAX,NESHMAX,NDIM,NDOF,NSPMAX,IBAK
      INTEGER*4 NADDHOLESMAX,NprdBndMAX,IMTF
      PARAMETER (ZERO=0.00d0,                                        &
                 HALF=0.5d0,                                         &
                 ONE=1.00d0,                                         &
                 TWO=2.00d0,                                         &
                 NDIM=2,                                             &
                 NDOF=4,                                             &
                 NSHMAX=10,                                          &
                 NPSHMAX=500,                                        &
                 NADDHOLESMAX=10,                                    &
                 NprdBndMAX=1,                                       &
                 NESHMAX=NPSHMAX-1,                                  &
                 NSPMAX=12,                                          &
                 PI=3.141593d0)

!     .. Common area  ..

      REAL*8  EPS,SHRELAX,GA,GM1,SNDMIN,DXCELL,                      &
              CADDhole(NDIM, NADDHOLESMAX),FLT_Dspeed
      INTEGER*4 NADDHOLES,NprdBnd,prdBndclr(3,NprdBndMAX)
      COMMON/PARAMT/CADDhole,EPS,SNDMIN,DXCELL,GA,SHRELAX,           &
                    GM1,IBAK,NADDHOLES,NprdBnd,prdBndclr,            &
                    FLT_Dspeed,IMTF
