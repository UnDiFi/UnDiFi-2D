c     NSHMAX      . max number od shocks
c     NPSHMAX     : max number of shock points for each shock. 
c     NESHMAX     : max number of shock element for each shocks
c     EPS         : displacement between the two faces of one internal boundary
c     GA          : scecific heat ratio
c     GM1         : GM1=GA-1
c     SNDMIN      : Max normalized distance of a phantom point
c     DXCELL      : length of a cell side
c     CFL         : stability coefficient for the shock motion
c     IBAK        : save solution every IBAK iterations
C     NADDHOLESMAX: MAX of additional hole points
c     NADDHOLES   : Number of additional hole points
C     CADDhole    : Coordinates of additional hole point 
c     NprdBndMAX  : MAX of periodic boundaries
c     NprdBnd     : Number of periodic boundaries 
c     prdBndclr   : array of the pair of colors of each periodic boundaries
c     FLT_Dspeed  : filter on discontinuity speeds  
c     IMTF        : iteration of mesh topology freezing 
     
C     .. Parameters ..

      REAL*8  ZERO,HALF,ONE,TWO,PI
      INTEGER*4 NSHMAX, NPSHMAX,NESHMAX,NDIM,NDOF,NSPMAX,IBAK
      INTEGER*4 NADDHOLESMAX,NprdBndMAX,IMTF
      PARAMETER (ZERO=0.00d0,
     +           HALF=0.5d0,
     +           ONE=1.00d0,
     +           TWO=2.00d0,
     +           NDIM=2,
     +           NDOF=4,
     +           NSHMAX=10,
     +           NPSHMAX=500,
     +           NADDHOLESMAX=10,
     +           NprdBndMAX=1,
     +           NESHMAX=NPSHMAX-1,
     +           NSPMAX=12,
c    +           EPS=0.20d-4,     ! EPS=0.20d-4,
c    +           SNDMIN=0.30,
c    +           SNDMIN=0.20,     ! circular cylinder - Type IV
c    +           DXCELL=0.10,     ! circular cylinder
c    +           DXCELL=0.01,     ! regular reflection
c    +           DXCELL=0.015,    ! mach reflection
c    +           DXCELL=0.050,    ! circular cylinder - Type IV
c    +           DXCELL=0.01,     ! S-S interaction 
c
c S-S interaction
c    +           EPS=0.20d-4,     ! EPS=0.10d-4,
c    +           SNDMIN=0.20,     
c    +           DXCELL=0.01,      
c    +           SHRELAX=0.80,
c    +           IBAK=50,
c    +           GA=1.40d+0,
 
c mach reflection (test case Nasuti 1) 
c    +           EPS=0.20d-4,     ! EPS=0.20d-4,
c    +           DXCELL=0.015,
c    +           SNDMIN=0.30,
c    +           SHRELAX=0.90,
c    +           IBAK=100,
c    +           GA=1.40d+0,
c mach reflection
c    +           EPS=0.20d-4,     ! EPS=0.20d-4,
c    +           DXCELL=0.015,    
c    +           SNDMIN=0.30,
c    +           SHRELAX=0.91,
c    +           IBAK=50,
c    +           GA=1.40d+0,
c
c circular cylinder
c    +           EPS=0.20d-9,     ! EPS=0.20d-4,
c    +           SNDMIN=0.20,      
c    +           DXCELL=0.00625, !l0 DXCELL=0.1 !l1 DXCELL=0.05 !l2 DXCELL=0.025 !l3 DXCELL=0.0125 !l4 DXCELL=0.00625
c    +           SHRELAX=0.5,
c    +           IBAK=1000,         ! IBAK=50
c    +           GA=1.40d+0,
c circular cylinder - Type IV
c    +           EPS=0.2d-4,     ! EPS=0.10d-4,
c    +           DXCELL=0.030,   ! DXCELL=0.050
c    +           SNDMIN=0.27,
c    +           SHRELAX=0.3,
c    +           IBAK=10,
c    +           GA=1.40d+0,
c
c regular reflection
c    +           EPS=0.10d-4,     ! EPS=0.20d-4,
c    +           DXCELL=0.01,     
c    +           SNDMIN=0.30,
c    +           SHRELAX=0.91,
c    +           IBAK=50,
c    +           GA=1.40d+0,
c
c mach reflection (test case Ivanov-1)
c    +           EPS=0.20d-4,     ! EPS=0.20d-4,
c    +           DXCELL=0.00500,    ! DXCELL=0.010
c    +           SNDMIN=0.30,
c    +           SHRELAX=0.70,      ! SHRELAX=0.70
c    +           IBAK=50,
c    +           GA=1.67D+0,

c mach reflection (test case Ivanov-4)
c    +           EPS=0.20d-4,     ! EPS=0.20d-4,
c    +           DXCELL=0.00800,    ! DXCELL=0.010
c    +           SNDMIN=0.30,
c    +           SHRELAX=0.70,
c    +           IBAK=50,
c    +           GA=1.67D+0,

c Q1D
c    +           EPS=0.20d-9,     ! EPS=0.20d-4,
c    +           SNDMIN=0.20,
c    +           DXCELL=0.04,
c    +           SHRELAX=0.5,
c    +           IBAK=50,         ! IBAK=50
c    +           GA=1.40d+0,

c    +           GM1=GA-1.0d+0,
     +           PI=3.141593d0)

C     .. Common area  ..

      REAL*8  EPS,SHRELAX,GA,GM1,SNDMIN,DXCELL,
     +        CADDhole(NDIM, NADDHOLESMAX),FLT_Dspeed
      INTEGER*4 NADDHOLES,NprdBnd,prdBndclr(3,NprdBndMAX)
      COMMON/PARAMT/CADDhole,EPS,SNDMIN,DXCELL,GA,SHRELAX,
     +              GM1,IBAK,NADDHOLES,NprdBnd,prdBndclr,
     +              FLT_Dspeed,IMTF
