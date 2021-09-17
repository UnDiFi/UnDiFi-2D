!> \par Purpose
!>
!> Driver routine for moving the grid:
!> the behaviour depends on the integer \c IALE
!> which is stored in time.com
!>
!> @param[in] XYZ0 the coordinates of the fixed grid
!> @param[out] XYZNEW the coordinates of grid at time TIME
!> @param[in] NDIM dimension of the space
!> @param[in] NP number of gridpoints = NPOIN+NGHOST+NPNOD
!> @param[in] TIME the at the end of the current timestep
      SUBROUTINE MOVEXYZ(XYZ0,XYZNEW,NDIM,NP,TIME)
C
C
      IMPLICIT NONE
      INTEGER NDIM,NP
      DOUBLE PRECISION TIME
      DOUBLE PRECISION XYZ0(NDIM,*),XYZNEW(NDIM,*) 
      INCLUDE "time.com"
      GOTO(1,2,3,4,5,6,7,8,9)IALE
      STOP 'IALE has a wrong value' 
    1 CALL MOVEXYZ_1(XYZ0,XYZNEW,NDIM,NP,TIME)
      RETURN
    2 CALL MOVEXYZ_2(XYZ0,XYZNEW,NDIM,NP,TIME)
      RETURN
    3 CALL MOVEXYZ_3(XYZ0,XYZNEW,NDIM,NP,TIME)
      RETURN
    4 CALL MOVEXYZ_4(XYZ0,XYZNEW,NDIM,NP,TIME)
      RETURN
    5 CALL MOVEXYZ_5(XYZ0,XYZNEW,NDIM,NP,TIME)
      RETURN
    6 CALL MOVEXYZ_6(XYZ0,XYZNEW,NDIM,NP,TIME)
      RETURN
    7 CALL MOVEXYZ_7(XYZ0,XYZNEW,NDIM,NP,TIME)
      RETURN
    8 CALL MOVEXYZ_8(XYZ0,XYZNEW,NDIM,NP,TIME)
      RETURN
    9 CALL MOVEXYZ_9(XYZ0,XYZNEW,NDIM,NP,TIME)
      RETURN
      END
!> \par Purpose
!>
!> Rigidly oscillates the grid horizontally (along x)
!>
!> \f{eqnarray*}{
!> x^{n+1} &=& x^0 + a\, \int_{0}^{n\,\Delta t}\sin\left(\omega t\right)  \; \mathrm{d}t  \\\
!> &=& x^0 + \frac{a}{\omega} \left[1-\cos\left(\omega \, n\,\Delta t\right)\right] \quad n = 0,1,...N
!> \f}
!> The period of the oscillation coincides with the total running time \f$T\f$ so that
!> at time \f$ N\Delta t = T \f$ the grid will return to its initial position; it follows that
!> \f$
!> \omega = \frac{2\pi}{T} = \frac{2\pi}{N\Delta t} 
!> \f$
!>
!> the velocity \f$ a\f$ is set to 1
!>
!> Use:
!> \verbatim 
!> -ale_grid_motion_type 1
!> \endverbatim
!> to activate this kind of grid motion.
!>
!>
!> @param[in] XYZ0 the coordinates of the fixed grid
!> @param[out] XYZNEW the coordinates of grid at time TIME
!> @param[in] NDIM dimension of the space
!> @param[in] NP number of gridpoints = NPOIN+NGHOST+NPNOD
!> @param[in] TIME the at the end of the current timestep
      SUBROUTINE MOVEXYZ_1(XYZ0,XYZNEW,NDIM,NP,TIME)
C
C     the vortex: rigid motion
C
      IMPLICIT NONE
      INTEGER NDIM,NP
      DOUBLE PRECISION TIME
      DOUBLE PRECISION XYZ0(NDIM,*),XYZNEW(NDIM,*) 
      INCLUDE "constants.h"
      INCLUDE "time.com"
      INTEGER I,IPOIN,IFAIL
      DOUBLE PRECISION dxyz(3),OMEGA
      DOUBLE PRECISION b0,b1,a1
      DOUBLE PRECISION U0
      LOGICAL VERBOSE
      PARAMETER(VERBOSE=.FALSE.)
!
      OMEGA = TWO*PI/(ITSTEP*DELT)
      U0 = ONE
c
      IF(NDIM.EQ.2)THEN
         dxyz(1) = U0/OMEGA*(ONE-COS(OMEGA*TIME))
         dxyz(2) = ZERO
      ELSE
         dxyz(1) = ZERO
         dxyz(2) = ZERO
         dxyz(3) = U0/OMEGA*(ONE-COS(OMEGA*TIME))
      ENDIF
      WRITE(6,*)
      WRITE(6,*)' Rigidly moving gridpoints at t=',TIME
      WRITE(6,*)' dx,dy,dz =',(dxyz(I),I=1,NDIM),' omega = ',omega
      DO IPOIN= 1,NP 
         DO I = 1,NDIM 
            XYZNEW(I,IPOIN) = XYZ0(I,IPOIN) + DXYZ(I)
         ENDDO
      ENDDO
      WRITE(6,*)
!
      IF(VERBOSE)THEN
         CALL R8Mat_Print('General',' ',NDIM,NP,XYZ0,
     +            NDIM,'Fixed grid nodal coordinates ',IFAIL)
         CALL R8Mat_Print('General',' ',NDIM,NP,XYZNEW,
     +            NDIM,'New nodal coordinates ',IFAIL)
      ENDIF

      RETURN
      END
!> \par Purpose
!>
!> Deforms the grid for the vortex: deforming grid
!>
!> In particular, we consider a mesh which deforms according to the following mapping:
!> \f{eqnarray*}{
!> \left(x^{n+1}-x_0\right) &=& \left(x^0-x_0\right) \left[B_0+A_1\sin\left(\omega t\right)+B_1\cos\left(\omega t\right)\right] \\\
!> \left(y^{n+1}-y_0\right) &=& \left(y^0-y_0\right) \left[B_0+A_1\sin\left(\omega t\right)+B_1\cos\left(\omega t\right)\right]
!> \f}
!>
!> Point \f$x_0,y_0\f$ is set equal to the initial location of the vortex centre, which can be set at runtime using the option:
!> \verbatim
!> -ale_motion_origin [x0,y0,z0]
!> \endverbatim
!> One should use:
!> 1. (x0,y0) = (0.5d0,0.5d0) for the scalar case
!> 2. (x0,y0) = (0.d0,0.d0) for the compressible vortex case
!>
!> Moreover, we choose \f$B_0+B_1 = 1\f$ so that the grid is un-deformed at \f$t = 0\f$ and \f$t = T\f$ and \f$A_1 = 0\f$, for simplicity.
!> Finally, we set: \f$B_1 = 1/(3\pi)\f$.
!>
!> Use:
!> \verbatim 
!> -ale_grid_motion_type 2
!> \endverbatim
!> to activate this kind of grid motion.
!>
!>
!> @param[in] XYZ0 the coordinates of the fixed grid
!> @param[out] XYZNEW the coordinates of grid at time TIME
!> @param[in] NDIM dimension of the space
!> @param[in] NP number of gridpoints = NPOIN+NGHOST+NPNOD
!> @param[in] TIME the at the end of the current timestep
      SUBROUTINE MOVEXYZ_2(XYZ0,XYZNEW,NDIM,NP,TIME)
C
C
      IMPLICIT NONE
      INTEGER NDIM,NP
      DOUBLE PRECISION TIME
      DOUBLE PRECISION XYZ0(NDIM,*),XYZNEW(NDIM,*) 
      INCLUDE "constants.h"
      INCLUDE "time.com"
      INTEGER I,IPOIN,IFAIL
      DOUBLE PRECISION HELP,OMEGA
      DOUBLE PRECISION b0,b1,a1,x0
!
      HELP = ONE/DELT
      OMEGA = TWO*PI/(ITSTEP*DELT)
c
      a1 = zero
      b1 = one/(1.5d0*omega)
      b0 = one-b1
      help = (b0+a1*sin(omega*time)+b1*cos(omega*time))
      WRITE(6,*)
      WRITE(6,*)' Moving gridpoints at t=',TIME
      WRITE(6,*)' omega = ',omega
      DO IPOIN= 1,NP 
         DO I = 1,NDIM
            x0 = XYZ_C(I)
            XYZNEW(I,IPOIN) = x0 + (XYZ0(I,IPOIN)-x0) * help
         ENDDO
      ENDDO
      WRITE(6,*)
!
!     CALL R8Mat_Print('General',' ',NDIM,NP,XYZ0,
!    +            NDIM,'Fixed grid nodal coordinates ',IFAIL)
!     CALL R8Mat_Print('General',' ',NDIM,NP,XYZNEW,
!    +            NDIM,'New nodal coordinates ',IFAIL)

      RETURN
      END
!> \par Purpose
!>
!> Move the grid for the oscillating NACA.
!> When using:
!> \verbatim
!> -ale_grid move
!> \endverbatim
!> only the grid within an ellipse will be rigidly moved and it will only work with the grid stored in:
!> \verbatim
!> $HOME/grids/2D/aerofoils/naca0012/euler/NP2355_within_ellipse
!> \endverbatim
!> Otherwise, when using:
!> \verbatim 
!> -ale_grid laplace 
!> \endverbatim
!> it will (presumably) work with all NACA0012 grids. 
!>
!> Use 
!> \verbatim 
!> -ale_grid_motion_type 3
!> \endverbatim
!> to activate this kind of grid motion.
!>
!> @param[in] XYZ0 the coordinates of the fixed grid
!> @param[out] XYZNEW the coordinates of the grid at time TIME
!> @param[in] NDIM dimension of the space
!> @param[in] NP number of gridpoints = NPOIN+NGHOST+NPNOD
!> @param[in] TIME the at the end of the current timestep
      SUBROUTINE MOVEXYZ_3(XYZ0,XYZNEW,NDIM,NP,TIME)
C
      IMPLICIT NONE
      INTEGER NDIM,NP
      DOUBLE PRECISION TIME
      DOUBLE PRECISION XYZ0(NDIM,*),XYZNEW(NDIM,*) 
      INCLUDE "constants.h"
      INCLUDE "time.com"
      INTEGER I,IPOIN,IFAIL
      DOUBLE PRECISION HELP,dxdt,dydt,OMEGA,x,y,aa,bb

c
      double precision cost,sint,xsi,eta
      double precision kappa,th0,deg2rad,th
      parameter(aa=2.275d0,bb=1.7850d0)
      DOUBLE PRECISION ellipse
      ellipse(x,y) = (((x-0.5d0)/aa)**2)+((y/bb)**2)-1.d0
!     ellipse(x,y) = (x**2+y**2)-6.25d0
c
      deg2rad = acos(-1.d0)/180.d0
c
c     pitch
c
      th = APITCH(3) * sin(OPITCH(3)*time)
c
      cost=cos(-th)
      sint=sin(-th)
C
      WRITE(6,*)
      WRITE(6,*)' Moving gridpoints at t=',TIME
      WRITE(6,*)' theta =',th/deg2rad,' omega = ',OPITCH(3)
c
      do i = 1,NP
         if(ellipse(xyz0(1,i),xyz0(2,i)).LT.ZERO)then
!        if(.FALSE.)then
            xsi = xyz0(1,i)-XYZ_C(1)
            eta = xyz0(2,i)-XYZ_C(2)
            xyznew(1,i) =  XYZ_C(1) + xsi*cost-eta*sint
            xyznew(2,i) =  XYZ_C(2) + xsi*sint+eta*cost
         else
            xyznew(1,i) =  xyz0(1,i)
            xyznew(2,i) =  xyz0(2,i)
         endif
      enddo

      WRITE(6,*)
!
!     CALL R8Mat_Print('General',' ',NDIM,NP,XYZ0,
!    +            NDIM,'Fixed grid nodal coordinates ',IFAIL)
!     CALL R8Mat_Print('General',' ',NDIM,NP,XYZNEW,
!    +            NDIM,'New nodal coordinates ',IFAIL)

      RETURN
      END
!> \par Purpose
!>
!> Move the grid for the piston in cylinder, using the testcase published here [FLD870]
!>
!> [FLD870]: http://dx.doi.org/doi:10.1002/(SICI)1097-0363(19990815)30:7<865::AID-FLD870>3.0.CO;2-5 "E. Lefrançois, G. Dhatt, D. Vandromme, International Journal for Numerical Methods in Fluids Volume 30, (7), 1999"
!>
!>
!> \f[
!> x_{\mbox{piston}}/L = c + \frac{R}{L}\left[1-\cos\left(\omega^*t^*\right)\right] + 
!> \left(\frac{R}{2L}\right)^2\left[1-\cos\left(2\omega^*t^*\right)\right]
!> \f]
!> where
!> \f[
!> t^* = t \frac{\sqrt{ R T^0 }}{L}
!> \f]
!> is the dimensionless time and the dimensionless frequency \f$\omega^*\f$ is:
!> \f[
!> \omega^* = \Omega \frac{L}{\sqrt{ R T^0 }}
!> \f]
!> the dimensional frequency is \f$\Omega\f$ = 100 Hz.
!>
!> grid files are in:
!>
!> \verbatim
!> $GRIDS/2D/misc/piston
!> \endverbatim
!>
!>
!> The remeshing process is made by linear interpolation of node displacement between the
!> fixed wall on the left and the moving wall on the right, i.e.
!>
!> \f[
!> x_i\left(t\right) = \frac{x_i\left(t=0\right)}{x_{\mbox{piston}}\left(t=0\right)} \; x_{\mbox{piston}} \left(t\right)
!> \f]
!>
!> Use:
!> \verbatim 
!> -ale_grid_motion_type 4
!> \endverbatim
!> to activate this kind of grid motion.
!>
!>
!>
!> @param[in] XYZ0 the coordinates of the fixed grid
!> @param[out] XYZNEW the coordinates of the grid at time TIME
!> @param[in] NDIM dimension of the space
!> @param[in] NP number of gridpoints = NPOIN+NGHOST+NPNOD
!> @param[in] TIME the at the end of the current timestep
      SUBROUTINE MOVEXYZ_4(XYZ0,XYZNEW,NDIM,NP,TIME)
C
      IMPLICIT NONE
      INTEGER NDIM,NP
      DOUBLE PRECISION TIME
      DOUBLE PRECISION XYZ0(NDIM,*),XYZNEW(NDIM,*) 
      INCLUDE "paramt.h"
      INCLUDE "constants.h"
      INCLUDE "time.com"
      INCLUDE "stream.com"
      INTEGER I,IPOIN,IFAIL
      DOUBLE PRECISION HELP,alpha,OMEGA,r,c,h,l,xp,c0,c1,c2
      parameter(r=30.d0,c=r,h=37.5d0,l=100.d0,c0=c/l,c1=r/l,
     &c2=(HALF*r/l)**2)
c
      omega = 100.d0 * lref/uref
      help = omega*TIME
c
C
      WRITE(6,*)
      WRITE(6,*)' Moving gridpoints at t=',TIME
      WRITE(6,*)' omega = ',omega
c
      xp = c0 + c1*(ONE-cos(HELP))+ c2*(ONE-cos(TWO*HELP))
      do i = 1,NP
         alpha = (xyz0(1,i)*l)/c
         xyznew(1,i) = alpha * xp
         xyznew(2,i) = xyz0(2,i)
      enddo

      WRITE(6,*)
!
!     CALL R8Mat_Print('General',' ',NDIM,NP,XYZ0,
!    +            NDIM,'Fixed grid nodal coordinates ',IFAIL)
!     CALL R8Mat_Print('General',' ',NDIM,NP,XYZNEW,
!    +            NDIM,'New nodal coordinates ',IFAIL)

      RETURN
      END
!> \par Purpose
!>
!> Move the grid for the Smooth inviscid flow in a piston (2D) testcase from J. Dobes' phD thesis
!>
!> grid files are in:
!>
!> \verbatim
!> $GRIDS/2D/misc/piston2 
!> $GRIDS/3D/misc/piston2 
!> \endverbatim
!>
!> The piston starts to accelerate with derivative of acceleration \f$ \left(\frac{\mathrm{d}a}{\mathrm{d}t}\right)_{t=0} = 0.2 \f$
!> so that the trajectory of the piston is:
!>
!> \f[
!> x_{\mbox{piston}} = \left(\frac{\mathrm{d}a}{\mathrm{d}t}\right)_{t=0} \frac{t^3}{6}
!> \f]
!>
!> The numerical solution should be plotted at time t = 4, when the piston has reached position x = 2.133.
!>
!>
!> The remeshing process is made by linear interpolation of node displacement between the
!> fixed wall on the right and the moving wall on the left, i.e.
!>
!> \f[
!> x_i\left(t\right) = \left(1-\alpha\right) x_{\mbox{piston}} \left(t\right) + \alpha l
!> \f]
!> where:
!> \f[
!> \alpha = x_i\left(t=0\right)/l
!> \f]
!>
!> in the 3D case it is the \f$ z \f$ axis that is aligned with the cylinder's axis
!>
!> Use:
!> \verbatim 
!> -ale_grid_motion_type 5
!> \endverbatim
!> to activate this kind of grid motion.
!>
!> @param[in] XYZ0 the coordinates of the fixed grid
!> @param[out] XYZNEW the coordinates of the grid at time TIME
!> @param[in] NDIM dimension of the space
!> @param[in] NP number of gridpoints = NPOIN+NGHOST+NPNOD
!> @param[in] TIME the at the end of the current timestep
      SUBROUTINE MOVEXYZ_5(XYZ0,XYZNEW,NDIM,NP,TIME)
C
      IMPLICIT NONE
      INTEGER NDIM,NP
      DOUBLE PRECISION TIME
      DOUBLE PRECISION XYZ0(NDIM,*),XYZNEW(NDIM,*) 
      INCLUDE "paramt.h"
      INCLUDE "constants.h"
      INCLUDE "time.com"
      INCLUDE "stream.com"
      INTEGER I,IPOIN,IFAIL
      DOUBLE PRECISION alpha,dadt,l,xp
      parameter(l=5.d0,dadt=0.2d0)
c
      xp = (dadt/6.d0) * (TIME*TIME*TIME)
C
      WRITE(6,*)
      WRITE(6,*)' Moving gridpoints at t=',TIME
      WRITE(6,*)' the piston is at x = ',xp
c
      IF(NDIM.EQ.2)THEN
         do i = 1,NP
            alpha = (xyz0(1,i)/l)
            xyznew(1,i) = (ONE-alpha) * xp + alpha * l
            xyznew(2,i) = xyz0(2,i)
         enddo
      ELSEIF(NDIM.EQ.3)THEN
         do i = 1,NP
            alpha = (xyz0(3,i)/l)
            xyznew(3,i) = (ONE-alpha) * xp + alpha * l
            xyznew(1,i) = xyz0(1,i)
            xyznew(2,i) = xyz0(2,i)
         enddo
      ENDIF

      WRITE(6,*)
!
!     CALL R8Mat_Print('General',' ',NDIM,NP,XYZ0,
!    +            NDIM,'Fixed grid nodal coordinates ',IFAIL)
!     CALL R8Mat_Print('General',' ',NDIM,NP,XYZNEW,
!    +            NDIM,'New nodal coordinates ',IFAIL)

      RETURN
      END
!> \par Purpose
!>
!> Move the grid for the expantion due to a smoothly moving piston (2D)
!>
!> grid files are in:
!>
!> \verbatim
!> $GRIDS/2D/misc/piston3
!> \endverbatim
!>
!> We consider a piston moving (to the left of the origin) with the following velocity:
!> \f[
!> v_p\left(t\right) = v_1 \left( 1-\exp\left(-\alpha t\right)\right) \quad\quad\mbox{with}\quad\quad \alpha = 1 \; \mbox{s}^{-1}
!> \f]
!> so that it will approach a constant velocity \f$v_1\f$ = -75 m/s after a sufficiently long time.
!>
!> The trajectory of the piston, having set \f$t_0 = 0\f$ as the initial time, is:
!> \f[
!> x_p\left(t\right) = v_1 \left[ t + \frac{1}{\alpha} \left( \exp\left(-\alpha t\right)-1\right)\right]
!> \f]
!>
!> The following settings should be used in the .petscrc file:
!> \verbatim
!> -nondimensionalisation internal
!> -inlet_total_temperature  223.99203583872573419611
!> -inlet_total_pressure  2.e5
!> -reference_length 1000.
!> -inlet_mach_number 0.0
!> \endverbatim
!>
!> \f$T^0\f$ has been chosen in such a way that the sound speed in the cylinder is initially \f$a_0\f$ = 300 m/s;
!> Total pressure is irrelevant.
!>
!> Reference quantities are:
!>
!> 1. velocity is \f$a_T = \sqrt{R T^0} = a_0/\sqrt{\gamma} \f$ = 253.54627641855498 m/s.
!> 2. time-scale is \f$T^* = L/a_T \f$ = 3.9440531887330774
!>
!> The numerical solution should be plotted at time \f$t_1^* \f$ = 0.75 which corresponds to a physical time \f$t_1 \f$ = 2.958039891549808050 s
!>
!> The various non-dimensional quantities are therefore:
!> \f{eqnarray*}{
!>  v_1^* &=& v_1/a_T = -0.29580398915498080 \\
!>  T^* &=& T/a_T = 3.9440531887330774 \\
!>  \alpha^* &=& \alpha/a_T = 3.9440531887330774 \\
!> \f}
!>
!>
!> The remeshing process is made by linear interpolation of node displacement between the
!> fixed wall on the right and the moving wall on the left, i.e.
!>
!> \f[
!> x_i\left(t\right) = \left(1-\beta\right) x_{\mbox{piston}} \left(t\right) + \beta l
!> \f]
!> where:
!> \f[
!> \beta = x_i\left(t=0\right)/l
!> \f]
!>
!> Use:
!> \verbatim 
!> -ale_grid_motion_type 6
!> \endverbatim
!> to activate this kind of grid motion/deformation.
!>
!> @param[in] XYZ0 the coordinates of the fixed grid
!> @param[out] XYZNEW the coordinates of the grid at time TIME
!> @param[in] NDIM dimension of the space
!> @param[in] NP number of gridpoints = NPOIN+NGHOST+NPNOD
!> @param[in] TIME at the end of the current timestep
      SUBROUTINE MOVEXYZ_6(XYZ0,XYZNEW,NDIM,NP,TIME)
C
      IMPLICIT NONE
      INTEGER NDIM,NP
      DOUBLE PRECISION TIME
      DOUBLE PRECISION XYZ0(NDIM,*),XYZNEW(NDIM,*) 
      INCLUDE "paramt.h"
      INCLUDE "constants.h"
      INCLUDE "time.com"
      INCLUDE "stream.com"
      INTEGER I,IPOIN,IFAIL
!
!     l (=1) gives the initial position of the closed end of the cylinder
!     xp gives the trajectory of the piston
!     
!
      DOUBLE PRECISION v1,alpha,v1star,alphastar,reftime,t,xstar,l,beta
      parameter(v1=-75.d0,alpha=ONE,l=ONE)
      DOUBLE PRECISION xp
      xp(t) = v1star*(t-(ONE-EXP(-ALPHASTAR*t))/ALPHASTAR)
c
      v1star = v1/uref
      reftime = lref/uref
      alphastar = alpha*reftime
      xstar = xp(TIME) ! the trajectory
C
      WRITE(6,*)
      WRITE(6,*)' Moving gridpoints at t=',TIME
      WRITE(6,*)' the piston is at x = ',xstar
      WRITE(6,*)' dimensionless alpha is = ',alphastar
      WRITE(6,*)' reference length scale = ',lref
      WRITE(6,*)' reference velocity is = ',uref
      WRITE(6,*)' reference time is = ',reftime
      WRITE(6,*)' dimensionless piston velocity is = ',v1star
      WRITE(6,*)
c
      do i = 1,NP
         beta = (xyz0(1,i)/l)
         xyznew(1,i) = (ONE-beta) * xstar + beta * l
         xyznew(2,i) = xyz0(2,i)
      enddo

      WRITE(6,*)
!
!     CALL R8Mat_Print('General',' ',NDIM,NP,XYZ0,
!    +            NDIM,'Fixed grid nodal coordinates ',IFAIL)
!     CALL R8Mat_Print('General',' ',NDIM,NP,XYZNEW,
!    +            NDIM,'New nodal coordinates ',IFAIL)

      RETURN
      END
!>
!> @param[in] XYZ0 the coordinates of the fixed grid
!> @param[out] XYZNEW the coordinates of the grid at time TIME
!> @param[in] NDIM dimension of the space
!> @param[in] NP number of gridpoints = NPOIN+NGHOST+NPNOD
!> @param[in] TIME at the end of the current timestep
!> \par Purpose
!>
!> Move the grid for the expansion interacting with a moving shock
!>
!> grid files are in:
!>
!> \verbatim
!> $GRIDS/2D/misc/piston4
!> \endverbatim
!>
!> A planar piston is impulsively set into motion with dimensionless velocity \f$v_p^*\f$, so that a shock
!> is created which moves to the right at constant speed \f$W\f$ = 450 m/s.
!>
!> The temperature in the undisturbed part of the cylinder is such that the
!> shock Mach number is \f$M_{r_1} = W/\sqrt{\gamma\,R\,T_1}\f$ = 1.534; it follows that: \f$T_1\f$ = 214,172830291 K.
!> 
!> We decide to solve the equations in non-dimensional form using the set of reference variables for internal flows:
!> The following settings should be used in the .petscrc file:
!> \verbatim
!> -nondimensionalisation internal
!> -inlet_total_temperature  247.92660034794133
!> -inlet_total_pressure  [whatever]
!> -reference_length 1.
!> -inlet_mach_number [whatever]
!> \endverbatim
!>
!> By doing so, we have to choose a (dimensional) reference temperature and pressure.
!> In the present case, these are set equal to the corresponding  (total) quantities in the gas at rest:
!> \f[
!> p_{\mbox{ref}} = \mbox{whatever} \quad\mbox{and}\quad T_{\mbox{ref}} = T_1 
!> \f]
!> The reference velocity (which is the isothermal sound speed) follows:
!> \f[
!> u_{\mbox{ref}} = \sqrt{R\,T_{\mbox{ref}}} = 247.92660034794133 \; \mbox{m/s}
!> \f]
!> Having set the upstream uniform state as the reference one, the non-dimensional upstream state becomes:
!> \f[
!> p_1^* = 1 \quad\quad \rho_1^* = 1 \quad\quad T_1^* = 1
!> \f]
!> using the jump relations in a reference frame moving with the shock, the (non-dimensional) downstream values are found to be:
!> \f[
!> p_2^* = 2.5786821287522521 \quad\quad \rho_2^* = 1.9201192819084407 \quad\quad T_2^* = 1.3429801747468804.
!> \f]
!> Concerning the shock-downstream velocity, which is the same as the piston velocity, we have:
!> \f[
!> v_p = u_2 = W - M_{r_2} a_2
!> \f]
!> with \f$M_{r_2}\f$ = 0.68938584364157141.
!> Therefore, the dimensionless piston velocity is:
!> \f[
!> v_p^* = \frac{u_2}{\sqrt{R T_1}} = \frac{W-M_{r_2} a_2}{\sqrt{R T_1}} = 
!> \sqrt{\gamma} \left( M_{r_1} - \sqrt{\frac{T_2}{T_1}} M_{r_2} \right) = 0.86977177750456247
!> \f]
!> The initial dimensions of the computational domain are: \f$[0,5L] \times [0,L]\f$.
!>
!> At time \f$t = 0\f$, when the shock is at \f$x/L\f$ = 0.3, the piston's speed starts decreasing
!> with the following law:
!> \f[
!> v_p^*\left(t\right) = v_p^* \, e^{-t}
!> \f]
!> so that the trajectory of the piston is:
!> \f[
!> x_{\mbox{piston}}/L = v_p^* \, \left(1-e^{-t} \right)
!> \f]
!>
!>
!>
!>
!> The remeshing process is made by linear interpolation of node displacement between the
!> fixed wall on the right and the moving wall on the left, i.e.
!>
!> \f[
!> x_i\left(t\right) = \left(1-\beta\right) x_{\mbox{piston}} \left(t\right) + \beta l
!> \f]
!> where:
!> \f[
!> \beta = x_i\left(t=0\right)/l
!> \f]
!>
!> Use:
!> \verbatim 
!> -ale_grid_motion_type 7
!> \endverbatim
!> to activate this kind of grid motion/deformation.
      SUBROUTINE MOVEXYZ_7(XYZ0,XYZNEW,NDIM,NP,TIME)
C
      IMPLICIT NONE
      INTEGER NDIM,NP
      DOUBLE PRECISION TIME
      DOUBLE PRECISION XYZ0(NDIM,*),XYZNEW(NDIM,*) 
      INCLUDE "paramt.h"
      INCLUDE "constants.h"
      INCLUDE "dofs.com"
      INCLUDE "stream.com"
      INCLUDE "time.com"
      INTEGER I,IPOIN,IFAIL
!
!     l (=1) gives the initial position of the closed end of the cylinder
!     xp gives the trajectory of the piston
!     
!
      DOUBLE PRECISION v1,v1star,reftime,t,xstar,l,beta
      parameter(l=5.d0,v1star=0.86977177750456247d0)
      DOUBLE PRECISION xp
      xp(t) = v1star*(ONE-EXP(-t))
c
      xstar = xp(TIME) ! the trajectory
      reftime = lref/uref
C
      WRITE(6,*)
      WRITE(6,*)' Moving gridpoints at t=',TIME
      WRITE(6,*)' the piston is at x = ',xstar
      WRITE(6,*)' reference length scale = ',lref
      WRITE(6,*)' reference velocity is = ',uref
      WRITE(6,*)' reference time is = ',reftime
      WRITE(6,*)' dimensionless piston velocity is = ',v1star
      WRITE(6,*)
c
      do i = 1,NP
         beta = (xyz0(1,i)/l)
         xyznew(1,i) = (ONE-beta) * xstar + beta * l
         xyznew(2,i) = xyz0(2,i)
      enddo

      WRITE(6,*)
!
!     CALL R8Mat_Print('General',' ',NDIM,NP,XYZ0,
!    +            NDIM,'Fixed grid nodal coordinates ',IFAIL)
!     CALL R8Mat_Print('General',' ',NDIM,NP,XYZNEW,
!    +            NDIM,'New nodal coordinates ',IFAIL)

      RETURN
      END
!>
!> @param[in] XYZ0 the coordinates of the fixed grid
!> @param[out] XYZNEW the coordinates of the grid at time TIME
!> @param[in] NDIM dimension of the space
!> @param[in] NP number of gridpoints = NPOIN+NGHOST+NPNOD
!> @param[in] TIME at the end of the current timestep
!> \par Purpose
!>
!> Move the grid for the centred compression: all \f$ \mathcal{C}^+\f$ characteristic lines meet
!> in the focal point F; the focal point has dimensionless coordinates: \f$ x^*_F, t^*_F \f$.
!>
!> Grid files are in:
!>
!> \verbatim
!> $GRIDS/2D/misc/piston4/gridlev?
!> \endverbatim
!>
!> The computational domain is the rectangle \f$ \left[0,5 \right] \times \left[0,1\right]\f$.
!>
!> We choose \f$ x^*_F = 5 \f$ so that \f$ t^*_F = x^*_F / a_0^* = 4.22577127364258288756 \f$.
!>
!> The piston's trajectory, velocity and acceleration are given below:
!>
!> \f{eqnarray}{
!>  x^*_p\left(t\right) &=& x^*_F + \frac{2 \, a_0^*}{\gamma-1}
!> \left(t^*_F-t^*\right) \left\{ 1 - \frac{\gamma+1}{2} \left[ \left(1-\frac{t^*}{t_f}\right)^{-\frac{\gamma-1}{\gamma+1}} \right] \right\} \\
!>  \dot{x}^*_p\left(t\right) &=& \frac{2 a_0^*}{\gamma-1}
!> \left\{ \left(1-\frac{t^*}{t_f}\right)^{-\frac{\gamma-1}{\gamma+1}} -1 \right\} \\
!>  \ddot{x}^*_p\left(t\right) &=& \frac{2 a_0^*}{\gamma+1}
!> \left(1-\frac{t^*}{t_f}\right)^{-\frac{2\,\gamma}{\gamma+1}}
!> \f}
!>
!> Before the piston is set into motion, the (dimensional) temperature in the cylinder is \f$T_0\f$ = 273.15 K.
!> 
!> We solve the equations in non-dimensional form using the set of reference variables for internal flows;
!> the following settings should be used in the .petscrc file:
!> \verbatim
!> -nondimensionalisation internal
!> -inlet_total_temperature  273.15
!> -inlet_total_pressure  [whatever]
!> -reference_length 1.
!> -inlet_mach_number 0.d0
!> \endverbatim
!>
!> By doing so, we have to choose a (dimensional) reference temperature and pressure.
!> In the present case, these are set equal to the corresponding  (total) quantities in the gas at rest:
!> \f[
!> p_{\mbox{ref}} = \mbox{whatever} \quad\mbox{and}\quad T_{\mbox{ref}} = T_1 
!> \f]
!> The reference velocity (which is the isothermal sound speed) follows:
!> \f[
!> u_{\mbox{ref}} = \sqrt{R\,T_{\mbox{ref}}} = 247.92660034794133 \; \mbox{m/s}
!> \f]
!> Having set the upstream uniform state as the reference one, the non-dimensional upstream state becomes:
!> \f[
!> p_1^* = 1 \quad\quad \rho_1^* = 1 \quad\quad T_1^* = 1
!> \f]
!>
!> Using characteristic theory, it is not difficult to show that for a given point \f$x^*,t^*\f$ in the space-time plane
!>
!> \f{eqnarray}{
!>  u^*\left(x^*,t^*\right) &=& \dot{x}^*_p\left(\tau\right) \\
!>  a^*\left(x^*,t^*\right) &=& a_p\left(\tau\right) = a^*_0 + \delta \dot{x}^*_p\left(\tau\right)
!> \f}
!> 
!> \f$\tau < t\f$ being the time instant when the straight characteristic line passing through \f$x^*,t^*\f$ meets the piston's trajectory.
!>
!> In order to find \f$\tau^*\f$, one has to solve the following non-linear equation (which represents the equation
!> of the straight characteristic line \f$\mathcal{C}^+\f$ in the unknown \f$\tau^*\f$
!>
!>
!>
!> The remeshing process is made by linear interpolation of node displacement between the
!> fixed wall on the right and the moving wall on the left, i.e.
!>
!> \f[
!> x_i\left(t\right) = \left(1-\beta\right) x_{\mbox{piston}} \left(t\right) + \beta l
!> \f]
!> where:
!> \f[
!> \beta = x_i\left(t=0\right)/l
!> \f]
!>
!> Use:
!> \verbatim 
!> -ale_grid_motion_type 8
!> \endverbatim
!> to activate this kind of grid motion/deformation.
      SUBROUTINE MOVEXYZ_8(XYZ0,XYZNEW,NDIM,NP,TIME)
C
      IMPLICIT NONE
      INTEGER NDIM,NP
      DOUBLE PRECISION TIME
      DOUBLE PRECISION XYZ0(NDIM,*),XYZNEW(NDIM,*) 
      INCLUDE "paramt.h"
      INCLUDE "constants.h"
      INCLUDE "dofs.com"
      INCLUDE "stream.com"
      INCLUDE "pfcgas.com"
      INCLUDE "time.com"
      INTEGER I,IPOIN,IFAIL
!
!     l (=1) gives the initial position of the closed end of the cylinder
!     xp gives the trajectory of the piston
!     
!
      DOUBLE PRECISION v1,xf,t,tf,reftime,xstar,l,beta,GP1
      parameter(l=5.d0,xf=l)
!
      DOUBLE PRECISION xp
      xp(t) = xf + TWO*SQRT(GAM)/GM1*(tf-t)*
     &(ONE-HALF*(GAM+ONE)*((ONE-t/tf)**(-GM1/(GAM+ONE))))
c
      tf = xf/sqrt(GAM)
      xstar = xp(TIME) ! the trajectory
      reftime = lref/uref
C
      WRITE(6,*)
      WRITE(6,*)' Moving gridpoints at t=',TIME
      WRITE(6,*)' the piston is at x/z = ',xstar
      WRITE(6,*)' reference length scale = ',lref
      WRITE(6,*)' reference velocity is = ',uref
      WRITE(6,*)' reference time is = ',reftime
      WRITE(6,*)
c
      IF(NDIM.EQ.2)THEN 
         do i = 1,NP
            beta = (xyz0(1,i)/l)
            xyznew(1,i) = (ONE-beta) * xstar + beta * l
            xyznew(2,i) = xyz0(2,i)
         enddo
      ELSE
         do i = 1,NP
            beta = (xyz0(3,i)/l)
            xyznew(3,i) = (ONE-beta) * xstar + beta * l
            xyznew(1,i) = xyz0(1,i)
            xyznew(2,i) = xyz0(2,i)
         enddo
      ENDIF

      WRITE(6,*)
!
!     CALL R8Mat_Print('General',' ',NDIM,NP,XYZ0,
!    +            NDIM,'Fixed grid nodal coordinates ',IFAIL)
!     CALL R8Mat_Print('General',' ',NDIM,NP,XYZNEW,
!    +            NDIM,'New nodal coordinates ',IFAIL)

      RETURN
      END
!> \par Purpose
!>
!> Move the grid for the pitching and plunging NACA.
!>
!> The plunging motion is given by:
!> \f[
!> h\left(t\right) = h_0 \sin \left(k\,t\right)
!> \f]
!> where:
!> \f[
!> h_0 = 0.75
!> \f]
!> is currently hardwired into the code.
!>
!> The pitching motion is given by:
!> \f[
!> \theta\left(t\right) = \theta_0 \sin \left(k\,t+\phi\right)
!> \f]
!> where \f$ \theta_0 \f$ can be set through the variable \c APITCH(3) using the runtime option
!> \verbatim
!> -ale_pitching_amplitude APITCH(1),APITCH(2),APITCH(3)
!> \endverbatim
!> and \f$ k\f$ can be set through the variable \c OPITCH(3) using the runtime option:
!> \verbatim
!> -ale_pitching_omega OPITCH(1),OPITCH(2),OPITCH(3)
!> \endverbatim
!>
!> The phase shift \f$ \phi \f$ is also currently hardwired into the code and set equal to 0.
!>
!> Use:
!> \verbatim
!> -ale_grid laplace
!> \endverbatim
!> to move the grid. 
!>
!> One suitable grid can be found in:
!> \verbatim
!> $HOME/grids/2D/aerofoils/naca0012/euler/pitching/nohanging/
!> \endverbatim
!>
!> Use 
!> \verbatim 
!> -ale_grid_motion_type 9
!> \endverbatim
!> to activate this kind of grid motion.
!>
!> @param[in] XYZ0 the coordinates of the fixed grid
!> @param[out] XYZNEW the coordinates of the grid at time TIME
!> @param[in] NDIM dimension of the space
!> @param[in] NP number of gridpoints = NPOIN+NGHOST+NPNOD
!> @param[in] TIME the at the end of the current timestep
      SUBROUTINE MOVEXYZ_9(XYZ0,XYZNEW,NDIM,NP,TIME)
C
      IMPLICIT NONE
      INTEGER NDIM,NP
      DOUBLE PRECISION TIME
      DOUBLE PRECISION XYZ0(NDIM,*),XYZNEW(NDIM,*) 
      INCLUDE "constants.h"
      INCLUDE "time.com"
      INTEGER I,IPOIN,IFAIL
      DOUBLE PRECISION x,y,hh,offset
c
      double precision cost,sint,xsi,eta
      double precision deg2rad,th
      double precision circle
      circle(x,y) = x*x+y*y-1.d0
c
      deg2rad = pi/180.d0
c
c     pitch
c
      offset = 75.d0*deg2rad
      offset = ZERO
      th = APITCH(3) * sin(OPITCH(3)*time+offset)
      hh = 0.75d0 * sin(OPITCH(3)*time)
c
      cost=cos(-th)
      sint=sin(-th)
C
      WRITE(6,*)
      WRITE(6,*)' Moving gridpoints at t=',TIME
      WRITE(6,*)' theta =',th/deg2rad,' omega = ',OPITCH(3),' plunge = '
     &,hh
      WRITE(6,*)
c
      do i = 1,NP
         if(circle(xyz0(1,i),xyz0(2,i)).LT.ZERO)then
            xsi = xyz0(1,i)-XYZ_C(1)
            eta = xyz0(2,i)-XYZ_C(2)
            xyznew(1,i) =  XYZ_C(1) + xsi*cost-eta*sint
            xyznew(2,i) =  XYZ_C(2) + xsi*sint+eta*cost + hh
         else
            xyznew(1,i) =  xyz0(1,i)
            xyznew(2,i) =  xyz0(2,i)
         endif
      enddo
!
!
!     CALL R8Mat_Print('General',' ',NDIM,NP,XYZ0,
!    +            NDIM,'Fixed grid nodal coordinates ',IFAIL)
!     CALL R8Mat_Print('General',' ',NDIM,NP,XYZNEW,
!    +            NDIM,'New nodal coordinates ',IFAIL)
      RETURN
      END
