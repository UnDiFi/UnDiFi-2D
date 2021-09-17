subroutine physical_values(Z1D,Z2D,Z3D,Z4D,Z1U,Z2U,Z3U,Z4U)

  implicit none

  real*8,parameter:: gm = 1.4                  ! Gamma

  real*8::rho1, u1, v1, p1, h1, tH1, M1, a1    ! Left state
  real*8::rho2, u2, v2, p2, h2, tH2, M2, a2    ! Right state
  real*8::a                                    ! Mach number
  real*8::us,Ms

  real*8,intent(out)::Z1D, Z2D, Z3D, Z4D       ! Downstream Roe's variables
  real*8,intent(out)::Z1U, Z2U, Z3U, Z4U       ! Upstream Roe's variables

  integer::option

  print*,' Choose an option '
  print*,' 1 / Steady planar shock '
  print*,' 2 / Moving Planar shock (case 6) '
  print*,' 3 / Moving Planar shock (case 15) '
  print*,' 4 / Moving Planar shock (physical values of simple duct)'
  print*,' 5 / Piston'
  read*,option

  select case (option)

  case(1)
     ! -----------
     ! Left state
     ! -----------
     rho1 = 1.4 
     p1   = 1.0 
     u1   = 2.0 
     v1   = 0.0 
     h1   = gm*p1/( ( gm - 1.0 )*rho1 ) 
     a    = sqrt(gm*p1/rho1) 
     M1   = u1/a 
     tH1  = (gm*p1)/((gm-1)*rho1) + (u1*u1 + v1*v1)/2. ! Total enthalpy

     ! ------------
     ! Right state
     ! ------------
     rho2 = 6.0*M1*M1*rho1/(M1*M1+5.0)
     u2   = (rho1*u1)/rho2
     p2   = (rho1*u1*u1+p1)-(rho2*u2*u2)
     v2   = 0.0
     h2   = h1+(u1*u1-u2*u2)/2.0 
     tH2  = (gm*p2)/((gm-1)*rho2) + (u2*u2 + v2*v2)/2. ! Total enthapy

     ! ==============================
     ! Conversion in Roe's variables 
     ! ==============================

     ! -----------
     ! Left state
     ! -----------
     Z1D = sqrt(rho1)
     Z2D = sqrt(rho1)*tH1
     Z3D = sqrt(rho1)*u1
     Z4D = sqrt(rho1)*v1

     ! ------------
     ! Right state
     ! ------------
     Z1U = sqrt(rho2)
     Z2U = sqrt(rho2)*tH2
     Z3U = sqrt(rho2)*u2
     Z4U = sqrt(rho2)*v2

  case(2)
     ! -----------
     ! Left state
     ! -----------
     rho1 = 1.4 
     p1   = 1.0 
     u1   = 0.0 
     v1   = 0.0 
     h1   = gm*p1/( ( gm - 1.0 )*rho1 ) 
     !a    = sqrt(gm*p1/rho1) 
     M1   = 2. 
     tH1 = (gm*p1)/((gm-1)*rho1) + (u1*u1 + v1*v1)/2. ! Total enthalpy

     ! ------------
     ! Right state
     ! ------------
     rho2 = rho1*( gm + 1.0 )*M1*M1/( 2.0 + ( gm - 1.0 )*M1*M1 )
     u2   = u1 + 2.0/( gm + 1.0 )*( M1 - 1.0/M1 )*sqrt( gm*p1/rho1 )
     p2   = p1*( 2.0*gm*M1*M1/( gm + 1.0 ) - ( gm -1.0 )/( gm + 1.0 ) )
     v2   = 0.0
     h2   = gm*p2/( ( gm - 1.0 )*rho2 )
     tH2 = (gm*p2)/((gm-1)*rho2) + (u2*u2 + v2*v2)/2. ! Total enthapy

     ! ==============================
     ! Conversion in Roe's variables 
     ! ==============================

     ! -----------
     ! Left state
     ! -----------
     Z1U = sqrt(rho1)
     Z2U = sqrt(rho1)*tH1
     Z3U = sqrt(rho1)*u1
     Z4U = sqrt(rho1)*v1

     ! ------------
     ! Right state
     ! ------------
     Z1D = sqrt(rho2)
     Z2D = sqrt(rho2)*tH2
     Z3D = sqrt(rho2)*u2
     Z4D = sqrt(rho2)*v2 

  case(3)
     ! Right State
     rho2 = 1.39996 
     p2   = 0.999973
     u2   = 0.0
     v2   = 0.0
     h2   = gm*p2/( ( gm - 1.0 )*rho2 )
     tH2 = (gm*p2)/((gm-1)*rho2) + (u2*u2 + v2*v2)/2. ! Total enthapy

     a2  = sqrt(gm*p2/rho2) ;
     M2   = u2/a2 ;

     !Position and velocity of the shock 
     Ms = 6.0
     us = a2 * Ms

     ! Left State 
     M1   = Ms
     p1   = p2*( 1.0 + 2.0*gm*(M1*M1 - 1.0)/(gm + 1.0) )
     rho1 = rho2*(1.0+(gm+1.0)*p1/((gm-1.0)*p2))/((gm+1.0)/(gm-1.0)+p1/p2)
     u1   = us*(1.0-rho2/rho1) + rho2*u2/rho1
     v1   = 0
     h1   = gm*p1/( ( gm - 1.0 )*rho1 )
     tH1 = (gm*p1)/((gm-1)*rho1) + (u1*u1 + v1*v1)/2. ! Total enthalpy

     ! ==============================
     ! Conversion in Roe's variables 
     ! ==============================

     ! -----------
     ! Left state
     ! -----------
     Z1U = sqrt(rho2)
     Z2U = sqrt(rho2)*tH2
     Z3U = sqrt(rho2)*u2
     Z4U = sqrt(rho2)*v2

     ! ------------
     ! Right state
     ! ------------
     Z1D = sqrt(rho1)
     Z2D = sqrt(rho1)*tH1
     Z3D = sqrt(rho1)*u1
     Z4D = sqrt(rho1)*v1 

  case(4)
     ! Right State
     rho2 = 1.39996
     p2   = 0.999973
     u2   = 0.0
     v2   = 0.0
     h2   = 2.5
     tH2 = (gm*p2)/((gm-1)*rho2) + (u2*u2 + v2*v2)/2. ! Total enthapy

     ! Left State 
     p1   = 2.57687
     rho1 = 2.68632
     u1   = 0.734594
     v1   = 0
     h1   = 3.62721
     tH1 = (gm*p1)/((gm-1)*rho1) + (u1*u1 + v1*v1)/2. ! Total enthalpy

     ! ==============================
     ! Conversion in Roe's variables 
     ! ==============================

     ! -----------
     ! Left state
     ! -----------
     Z1U = sqrt(rho2)
     Z2U = sqrt(rho2)*tH2
     Z3U = sqrt(rho2)*u2
     Z4U = sqrt(rho2)*v2

     ! ------------
     ! Right state
     ! ------------
     Z1D = sqrt(rho1)
     Z2D = sqrt(rho1)*tH1
     Z3D = sqrt(rho1)*u1
     Z4D = sqrt(rho1)*v1 

   case(5)
     ! -----------
     ! Left state
     ! -----------
     Z1D = 0.16392692858310E+001
     Z2D = 0.59449690418206E+001
     Z3D = 1.20422585275454E+000
     Z4D = 0.00000000000000E+000

     ! ------------
     ! Right state
     ! ------------
     Z1U = 0.11832159566199E+001
     Z2U = 0.29580398915498E+001
     Z3U = 0.00000000000000E+000
     Z4U = 0.00000000000000E+000 

  end select

end subroutine physical_values
