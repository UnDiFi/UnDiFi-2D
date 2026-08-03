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
!    +           EPS=0.20d-4,     ! EPS=0.20d-4,
!    +           SNDMIN=0.30,
!    +           SNDMIN=0.20,     ! circular cylinder - Type IV
!    +           DXCELL=0.10,     ! circular cylinder
!    +           DXCELL=0.01,     ! regular reflection
!    +           DXCELL=0.015,    ! mach reflection
!    +           DXCELL=0.050,    ! circular cylinder - Type IV
!    +           DXCELL=0.01,     ! S-S interaction
!
! S-S interaction
!    +           EPS=0.20d-4,     ! EPS=0.10d-4,
!    +           SNDMIN=0.20,
!    +           DXCELL=0.01,
!    +           SHRELAX=0.80,
!    +           IBAK=50,
!    +           GA=1.40d+0,

! mach reflection (test case Nasuti 1)
!    +           EPS=0.20d-4,     ! EPS=0.20d-4,
!    +           DXCELL=0.015,
!    +           SNDMIN=0.30,
!    +           SHRELAX=0.90,
!    +           IBAK=100,
!    +           GA=1.40d+0,
! mach reflection
!    +           EPS=0.20d-4,     ! EPS=0.20d-4,
!    +           DXCELL=0.015,
!    +           SNDMIN=0.30,
!    +           SHRELAX=0.91,
!    +           IBAK=50,
!    +           GA=1.40d+0,
!
! circular cylinder
!    +           EPS=0.20d-9,     ! EPS=0.20d-4,
!    +           SNDMIN=0.20,
!    +           DXCELL=0.00625, !l0 DXCELL=0.1 !l1 DXCELL=0.05 !l2 DXCELL=0.025 !l3 DXCELL=0.0125 !l4 DXCELL=0.00625
!    +           SHRELAX=0.5,
!    +           IBAK=1000,         ! IBAK=50
!    +           GA=1.40d+0,
! circular cylinder - Type IV
!    +           EPS=0.2d-4,     ! EPS=0.10d-4,
!    +           DXCELL=0.030,   ! DXCELL=0.050
!    +           SNDMIN=0.27,
!    +           SHRELAX=0.3,
!    +           IBAK=10,
!    +           GA=1.40d+0,
!
! regular reflection
!    +           EPS=0.10d-4,     ! EPS=0.20d-4,
!    +           DXCELL=0.01,
!    +           SNDMIN=0.30,
!    +           SHRELAX=0.91,
!    +           IBAK=50,
!    +           GA=1.40d+0,
!
! mach reflection (test case Ivanov-1)
!    +           EPS=0.20d-4,     ! EPS=0.20d-4,
!    +           DXCELL=0.00500,    ! DXCELL=0.010
!    +           SNDMIN=0.30,
!    +           SHRELAX=0.70,      ! SHRELAX=0.70
!    +           IBAK=50,
!    +           GA=1.67D+0,

! mach reflection (test case Ivanov-4)
!    +           EPS=0.20d-4,     ! EPS=0.20d-4,
!    +           DXCELL=0.00800,    ! DXCELL=0.010
!    +           SNDMIN=0.30,
!    +           SHRELAX=0.70,
!    +           IBAK=50,
!    +           GA=1.67D+0,

! Q1D
!    +           EPS=0.20d-9,     ! EPS=0.20d-4,
!    +           SNDMIN=0.20,
!    +           DXCELL=0.04,
!    +           SHRELAX=0.5,
!    +           IBAK=50,         ! IBAK=50
!    +           GA=1.40d+0,

!    +           GM1=GA-1.0d+0,
                 PI=3.141593d0)

!     .. Common area  ..

      REAL*8  EPS,SHRELAX,GA,GM1,SNDMIN,DXCELL,                      &
              CADDhole(NDIM, NADDHOLESMAX),FLT_Dspeed
      INTEGER*4 NADDHOLES,NprdBnd,prdBndclr(3,NprdBndMAX)
      COMMON/PARAMT/CADDhole,EPS,SNDMIN,DXCELL,GA,SHRELAX,           &
                    GM1,IBAK,NADDHOLES,NprdBnd,prdBndclr,            &
                    FLT_Dspeed,IMTF
