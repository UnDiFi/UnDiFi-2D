!> \details
!> @param[in] XYZ0 the coordinates of the fixed grid
!> @param[out] XYZ0 the coordinates of grid at time TIME
!> @param[in] NDIM dimension of the space
!> @param[in] NP number of gridpoints = NPOIN+NGHOST+NPNOD
!> @param[in] TIME the at the end of the current timestep
      SUBROUTINE MOVEXYZ(XYZ0,XYZNEW,NDIM,NP,TIME)
      IMPLICIT NONE
      INTEGER NDIM,NP
      DOUBLE PRECISION TIME
      DOUBLE PRECISION XYZ0(NDIM,*),XYZNEW(NDIM,*) 
      INCLUDE "constants.h"
      INCLUDE "time.com"
      INTEGER I,IPOIN,IFAIL
      DOUBLE PRECISION HELP,dxdt,dydt,OMEGA,dx,dy,xp,yp

c
      double precision cost,sint,xsi,eta,pi,a0,da
      double precision kappa,h0,th0,phi,deg2rad,th
      parameter(aa=2.05d0,bb=1.505d0)
      DOUBLE PRECISION ellipse
      ellispe(x,y) = ((x-0.5d0)/aa)**2+((y/bb)**2)-1.d0
c
      deg2rad = pi/180.d0
      a0=0.016d0*deg2rad
      da=2.51d0*deg2rad
c
      phi = 75.*deg2rad
      th0 = da
      kappa = 0.0814d0
      h0 = 0.0d0
c
c     pitch and plunge
c
      xy0(1) = 0.5d0
      xy0(2) = 0.d0
      th = a0 + th0 * sin(2.d0*kappa*time)
c
      cost=cos(th)
      sint=sin(th)
C
      WRITE(6,*)
      WRITE(6,*)' Moving gridpoints at t=',TIME
      WRITE(6,*)' dx,dy =',dx,dy,' omega = ',omega
c
      do i = 1,NP
         if(ellipse(xyz0(1,i),xyz0(2,i)).LT.0.d0)then
            xsi = xyz0(1,i)-0.5d0
            eta = xyz0(2,i)
            xyznew(1,i) =  xy0(1) + xsi*cost-eta*sint
            xyznew(2,i) =  xy0(2) + xsi*sint+eta*cost
         else
            xyznew(1,i) =  xyz0(1,i)
            xyznew(2,i) =  xyz0(2,i)
         endif
      enddo

      WRITE(6,*)
!
!     CALL X04CAF('General',' ',NDIM,NP,XYZ0,
!    +            NDIM,'Fixed grid nodal coordinates ',IFAIL)
!     CALL X04CAF('General',' ',NDIM,NP,XYZNEW,
!    +            NDIM,'New nodal coordinates ',IFAIL)

      RETURN
      END
