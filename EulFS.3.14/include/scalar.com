      INTEGER ICASE
      COMMON /CSCALAR/ ICASE
C
C     Common block related to scalar 
C     advection-diffusion problems (|KAN|=1)
C               
C     ICASE is used to select different convection speeds depending
C           on its value and the space dimension NDIM
C
C     ICASE    |    u   |   v   |   w   |
C ---------------------2D------------------
C       1      |    2.  |   1.  |   -   |   linear convection
C       2      |    u   |   1.  |   -   |   Burger's eqn.
C       3      |    1.  |   1.  |   -   |   linear convection
C ---------------------3D------------------
C       1      |   0.75 | 0.875 |  1.0  |   linear convection
C       2      |    z   | 0.200 |  -x   |   spiral convection
C
