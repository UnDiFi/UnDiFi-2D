/***************************************************************************
                              initial_solution.c
                              ------------------
This is initial_solution: depending on the value of initial_state, it reads
          or just sets up the initial solution for the computation
                             --------------------
    begin                : Mon May 6 2002
    copyright            : (C) 2002 by Mario Ricchiuto
    email                : ricchiut@vki.ac.be,


***************************************************************************/
#include "common.h"

extern int initial_state ;

extern int NN ;

extern struct node_struct *node ;

extern double gm ;
extern double Mach_inf, alpha_inf, rho_inf , p_inf ; 

extern void read_solution() ;
extern void P_2_Z( double *, double * ) ;
extern void Z_2_P( double *, double * ) ;

void initial_solution()
{
  int n ;
  double Ms, Xs, rb, xc, yc, r ;
  double rho1, p1, u1, h1, v1 ;
  double rho2, p2, u2, h2, v2, omega, dp, pm, p4 ;
  double U0, U1, U2, U3 ;
  double U0b, U1b, U2b, U3b ;
  double V0, V1, V2, V3  ;
  double M1,a;
  double a_inf = sqrt( gm * p_inf / rho_inf ) ;


  switch ( initial_state )
    {
    case 0:
      read_solution() ;
      for ( n = 0 ; n < NN ; n ++ ) P_2_Z( node[n].P[2], node[n].Z[2] ) ;
      break ;
    case 1:/* Riemann Problem */
      for ( n = 0 ; n < NN ; n ++ )
	{
	  node[n].P[2][0] = 1. ;
	  node[n].P[2][1] = 0. ;
	  node[n].P[2][2] = 0. ;
	  node[n].P[2][3] = 0.001 ;
				 
	  if ( node[n].coordinate[0] == 0. ) 
	    if ( node[n].coordinate[1] == 0. )
	      node[n].P[2][3] = 2500. ;

	  // 10358.83/2.

	  P_2_Z( node[n].P[2], node[n].Z[2] ) ;
	}
      break ;

    case 2:/* Q-1D Oscillating Riemann Problem (Shu & Osher) */
      for ( n = 0 ; n < NN ; n ++ )
	{
	  if ( node[n].coordinate[0] > -4.0 )
	    {
	      p1 = 1.0 ;
	      u1 = 0.0 ;
	      v1 = 0.0 ;
	      h1 = 3.5/( 1.0 + 0.2*sin( 5.0*node[n].coordinate[0] ) ) ;
		
	      U0 = p1*gm/( h1*( gm - 1.0 ) ) ;

	      node[n].P[2][0] = U0  ;
	      node[n].P[2][1] = U0*u1 ;
	      node[n].P[2][2] = U0*v1 ;
	      node[n].P[2][3] = p1/( gm - 1.0 ) + 0.5*( u1*u1 + v1*v1 )*U0 ;
	    }
	  else
	    {
	      p1 = 10.333330 ;
	      u1 = 2.629367 ;
	      v1 = 0.0 ;
	      h1 = 9.3765398 ;
		
	      U0 = gm*p1/( ( gm - 1.0 )*h1 ) ;
	      U3 = h1/gm ;

	      node[n].P[2][0] = U0 ;
	      node[n].P[2][1] = U0*u1 ;
	      node[n].P[2][2] = U0*v1 ;
	      node[n].P[2][3] = U3*U0 + 0.5*( u1*u1 + v1*v1 )*U0 ;
	    }

	  P_2_Z( node[n].P[2], node[n].Z[2] ) ;
	}
      break ;
    case 3: /* Mach 3 Wind Tunnel with a Step */
      for ( n = 0 ; n < NN ; n ++ )
	{
	  p1 = 1.0 ;
	  u1 = 3.0 ;
	  v1 = 0.0 ;
	  h1 = 2.5 ;

	  if ( node[n].coordinate[0] == 0.6 )
	    if ( node[n].coordinate[1] <= 0.2 ) u1 = 0.0 ;
		
	  U0 = gm*p1/( ( gm - 1.0 )*h1 ) ;
	  U3 = h1/gm ;
		
	  node[n].P[2][0] = U0 ;
	  node[n].P[2][1] = U0*u1 ;
	  node[n].P[2][2] = U0*v1 ;
	  node[n].P[2][3] = U3*U0 + 0.5*( u1*u1 + v1*v1 )*U0 ;

	  P_2_Z( node[n].P[2], node[n].Z[2] ) ;
	}
      break ;
    case 4: /* 2D Riemann Problem on a Regular Mesh */
      for ( n = 0 ; n < NN ; n ++ )
	{
	  if ( node[n].coordinate[1] >= 0.8 )
	    {
	      if ( node[n].coordinate[0] >= 0.8 )
		{
		  p1 = 1.5 ;
		  u1 = 0.0 ;
		  v1 = 0.0 ;
		  h1 = 3.5 ;
		
		  U0 = gm*p1/( ( gm - 1.0 )*h1 ) ;
		  U3 = h1/gm ;
 
		  node[n].P[2][0] = U0 ;
		  node[n].P[2][1] = U0*u1 ;
		  node[n].P[2][2] = U0*v1 ;
		  node[n].P[2][3] = U3*U0 + 0.5*( u1*u1 + v1*v1 )*U0 ;
		}
	      else
		{
		  p1 = 0.3 ;
		  u1 = 1.2060454 ;
		  v1 = 0.0 ;
		  h1 = 1.9727271 ;
		
		  U0 = gm*p1/( ( gm - 1.0 )*h1 ) ;
		  U3 = h1/gm ;
 
		  node[n].P[2][0] = U0 ;
		  node[n].P[2][1] = U0*u1 ;
		  node[n].P[2][2] = U0*v1 ;
		  node[n].P[2][3] = U3*U0 + 0.5*( u1*u1 + v1*v1 )*U0 ;
		}
	    }
	  else
	    {
	      if ( node[n].coordinate[0] >= 0.8 )
		{
		  p1 = 0.3 ;
		  v1 = 1.2060454 ;
		  u1 = 0.0 ;
		  h1 = 1.9727271 ;
		
		  U0 = gm*p1/( ( gm - 1.0 )*h1 ) ;
		  U3 = h1/gm ;
 
		  node[n].P[2][0] = U0 ;
		  node[n].P[2][1] = U0*u1 ;
		  node[n].P[2][2] = U0*v1 ;
		  node[n].P[2][3] = U3*U0 + 0.5*( u1*u1 + v1*v1 )*U0 ;
		}
	      else
		{
		  p1 = 0.0290323 ;
		  u1 = 1.2060454 ;
		  v1 = 1.2060454 ;
		  h1 = 0.73636487 ;
		
		  U0 = gm*p1/( ( gm - 1.0 )*h1 ) ;
		  U3 = h1/gm ;
 
		  node[n].P[2][0] = U0 ;
		  node[n].P[2][1] = U0*u1 ;
		  node[n].P[2][2] = U0*v1 ;
		  node[n].P[2][3] = U3*U0 + 0.5*( u1*u1 + v1*v1 )*U0 ;
		}
	    }

	  P_2_Z( node[n].P[2], node[n].Z[2] ) ;
	}
      break ;
    case 5: /*Double Mach Reflection*/

      p1 = 116.5 ;
      u1 = 8.25*sqrt( 3.0 )/2.0  ;
      v1 = -8.25/2.0 ;
      h1 = 50.96875 ;
		
      U0 = gm*p1/( ( gm - 1.0 )*h1 ) ;
      U3 = h1/gm ;

      U1 = U0*u1 ;
      U2 = U0*v1 ;
      U3 = U0*U3 + 0.5*( u1*u1 + v1*v1 )*U0 ;
		     
      p2 = 1.0 ;
      u2 = 0.0 ;
      v2 = 0.0 ;
      h2 = 2.5 ;
		
      V0 = gm*p2/( ( gm - 1.0 )*h2 ) ;
      V3 = h2/gm ;

      V1 = V0*u2 ;
      V2 = V0*v2 ;
      V3 = V0*V3 + 0.5*( u2*u2 + v2*v2 )*V0 ;
		     
      for ( n = 0 ; n < NN ; n ++ )
	{
	  if ( node[n].coordinate[0] < 1.0/6.0 + node[n].coordinate[1]/( sqrt( 3.0 ) ) )
	    {
	      node[n].P[2][0] = U0 ; 	
	      node[n].P[2][1] = U1 ;	
	      node[n].P[2][2] = U2 ; 	
	      node[n].P[2][3] = U3 ;	
	    }
	  else			
	    {
	      node[n].P[2][0] = V0 ; 	
	      node[n].P[2][1] = V1 ;	
	      node[n].P[2][2] = V2 ; 	
	      node[n].P[2][3] = V3 ;	
	    }
		
	  P_2_Z( node[n].P[2], node[n].Z[2] ) ;
	}   	
      break ;
    case 6:/*Right Moving Planar Shock: Defined by Ms and Xs*/

      // Right state

      rho1 = 1.4 ;
      p1   = 1.0 ;
      u1   = 0.0 ;
      v1   = 0.0 ;
      h1   = gm*p1/( ( gm - 1.0 )*rho1 ) ;
		
      U0 = gm*p1/( ( gm - 1.0 )*h1 ) ;
      U3 = h1/gm ;

      U1 = U0*u1 ;
      U2 = U0*v1 ;
      U3 = U0*U3 + 0.5*( u1*u1 + v1*v1 )*U0 ;

      // Initial position and shock Mach number		    
 
      Xs   = .50 ;
      Ms   = 10. ; 

      // Left state from analytical jump conditions
      // for right-moving shocks

      rho2 = rho1*( gm + 1.0 )*Ms*Ms/( 2.0 + ( gm - 1.0 )*Ms*Ms ) ;
      p2   = p1*( 2.0*gm*Ms*Ms/( gm + 1.0 ) - ( gm -1.0 )/( gm + 1.0 ) ) ;
      u2   = u1 + 2.0/( gm + 1.0 )*( Ms - 1.0/Ms )*sqrt( gm*p1/rho1 ) ;
      v2   = 0.0 ;
      h2   = gm*p2/( ( gm - 1.0 )*rho2 ) ;

      V0 = gm*p2/( ( gm - 1.0 )*h2 ) ;
      V3 = h2/gm ;

      V1 = V0*u2 ;
      V2 = V0*v2 ;
      V3 = V0*V3 + 0.5*( u2*u2 + v2*v2 )*V0 ;

				
      for ( n = 0 ; n < NN ; n ++ )
	{
	  if ( node[n].coordinate[0] < Xs )
	    {
	      node[n].P[2][0] = V0 ; 	
	      node[n].P[2][1] = V1 ;	
	      node[n].P[2][2] = V2 ; 	
	      node[n].P[2][3] = V3 ;	
	    }
	  else			  
	    {
  
	      node[n].P[2][0] = U0 ; 	
	      node[n].P[2][1] = U1 ;	
	      node[n].P[2][2] = U2 ; 	
	      node[n].P[2][3] = U3 ;	
	    }
		
	  P_2_Z( node[n].P[2], node[n].Z[2] ) ;
	}                   	
      break ;
    case 7: /* 2D Riemann Problem on a Isotropic Mesh (half rotated domain) */
      for ( n = 0 ; n < NN ; n ++ )
	{
	  if ( node[n].coordinate[0] >= 0.8*sqrt(2.0) + node[n].coordinate[1] )
	    {
	      p1   = 1.5 ;
	      u1   = 0.0 ;
	      v1   = 0.0 ;
	      h1   = 3.5 ;
		
	      U0 = gm*p1/( ( gm - 1.0 )*h1 ) ;
	      U3 = h1/gm ;

	      U1 = U0*u1 ;
	      U2 = U0*v1 ;
	      U3 = U0*U3 + 0.5*( u1*u1 + v1*v1 )*U0 ;
		     
	      node[n].P[2][0] = U0 ; 	
	      node[n].P[2][1] = U1 ;	
	      node[n].P[2][2] = U2 ; 	
	      node[n].P[2][3] = U3 ;	
	    }
	  else
	    {
	      if ( node[n].coordinate[0] >= 0.8*sqrt(2.0) - node[n].coordinate[1] )
		{
		  p1 = 0.3 ;
		  u1 = 1.2060454*0.5*sqrt( 2.0 ) ;
		  v1 = -1.2060454*0.5*sqrt( 2.0 ) ;
		  h1 = 1.9727271 ;

		  U0 = gm*p1/( ( gm - 1.0 )*h1 ) ;
		  U3 = h1/gm ;
 
		  U1 = U0*u1 ;
		  U2 = U0*v1 ;
		  U3 = U0*U3 + 0.5*( u1*u1 + v1*v1 )*U0 ;
		     
		  node[n].P[2][0] = U0 ; 	
		  node[n].P[2][1] = U1 ;	
		  node[n].P[2][2] = U2 ; 	
		  node[n].P[2][3] = U3 ;	
		}
	      else
		{
		  p1 = 0.0290323 ;
		  u1 = 1.2060454*sqrt( 2.0 ) ;
		  v1 = 0.0 ;
		  h1 = 0.73636487 ;

		  U0 = gm*p1/( ( gm - 1.0 )*h1 ) ;
		  U3 = h1/gm ;
 
		  U1 = U0*u1 ;
		  U2 = U0*v1 ;
		  U3 = U0*U3 + 0.5*( u1*u1 + v1*v1 )*U0 ;
		     
		  node[n].P[2][0] = U0 ; 	
		  node[n].P[2][1] = U1 ;	
		  node[n].P[2][2] = U2 ; 	
		  node[n].P[2][3] = U3 ;	
		}
	    }

	  P_2_Z( node[n].P[2], node[n].Z[2] ) ;
	}
      break ; 	
    case 8: /* Abgrall's 2D Riemann Problem */
      p1 = 0.1 ;
      u1 = 0.0 ;
      v1 = 0.0 ;
      h1 = 1.4*p1/( 0.4*0.1 ) ;
		
      U0 = gm*p1/( ( gm - 1.0 )*h1 ) ;
      U3 = h1/gm ;

      U1 = U0*u1 ;
      U2 = U0*v1 ;
      U3 = U0*U3 + 0.5*( u1*u1 + v1*v1 )*U0 ;
		     
      p2 = 1. ;
      u2 = 0.0 ;
      v2 = 0.0 ;
      h2 = 1.4*p2/( 0.4*1.0 ) ;
		
      V0 = gm*p2/( ( gm - 1.0 )*h2 ) ;
      V3 = h2/gm ;

      V1 = V0*u2 ;
      V2 = V0*v2 ;
      V3 = V3*V0 + 0.5*( u2*u2 + v2*v2 )*V0 ;
		     
      for ( n = 0 ; n < NN ; n ++ )
	{
	  if ( node[n].coordinate[0]*node[n].coordinate[1] < 0. )
	    {
	      node[n].P[2][0] = U0 ;
	      node[n].P[2][1] = U1 ;
	      node[n].P[2][2] = U2 ;
	      node[n].P[2][3] = U3 ;
	    }
	  else
	    {
	      node[n].P[2][0] = V0 ;
	      node[n].P[2][1] = V1 ;
	      node[n].P[2][2] = V2 ;
	      node[n].P[2][3] = V3 ;
	    }

	  P_2_Z( node[n].P[2], node[n].Z[2] ) ;
	}
      break ; 		   
    case 9://Jirka's vortex
      for ( n = 0 ; n < NN ; n ++ )
	{
	  xc    = node[n].coordinate[0] - 0.5 ;	 
	  yc    = node[n].coordinate[1] - 0.5 ;	 
	  r     = sqrt( xc*xc + yc*yc ) ;
	  omega = 4.*pi*r ;
	  omega = 15.*( 1. + cos( omega ) ) ;			   
	  pm    = 100. ; 
	  r     = 4.*pi*r ;
	  dp    = ( 225.*1.4/( 16.*pi*pi ) )*( 2.*cos( r ) + 2.*r*sin( r )
					       + cos( 2.*r )/8. + r*sin( 2.*r )/4. + r*r*12./16. ) ;
	  dp   -= ( 225.*1.4/( 16.*pi*pi ) )*( 2.*cos( pi ) + 2.*pi*sin( pi )
					       + cos( 2.*pi )/8. + pi*sin( 2.*pi )/4. + pi*pi*12./16. ) ;          
			   
	  if ( r/( 4.*pi ) >= 0.25 ){ omega = 0. ; dp = 0. ;}

	  node[n].Z[2][0] = 1.4 ; 
	  node[n].Z[2][1] = 6. - omega*yc ;
	  node[n].Z[2][2] =  omega*xc ;
	  node[n].Z[2][3] = pm + dp ;

	  Z_2_P( node[n].Z[2], node[n].P[2] ) ;
	}
      break;
    case 10:/* Shock bubble interaction*/
      //http://www.math.ntnu.no/~holden/fronttrack/gas/sb/index.html
      // Right state outside the bubble
      //
      rho1 = 1. ;
      p1   = 1.0 ;
      u1   = 0.0 ;
      v1   = 0.0 ;
      h1   = gm*p1/( ( gm - 1.0 )*rho1 ) ;
		
      U0 = gm*p1/( ( gm - 1.0 )*h1 ) ;
      U3 = h1/gm ;

      U1 = U0*u1 ;
      U2 = U0*v1 ;
      U3 = U0*U3 + 0.5*( u1*u1 + v1*v1 )*U0 ;
      //
      // Right state inside the bubble
      //
      rho1 = 0.1 ;
      p1   = 1.0 ;
      u1   = 0.0 ;
      v1   = 0.0 ;
      h1   = gm*p1/( ( gm - 1.0 )*rho1 ) ;
		
      U0b = gm*p1/( ( gm - 1.0 )*h1 ) ;
      U3b = h1/gm ;

      U1b = U0b*u1 ;
      U2b = U0b*v1 ;
      U3b = U0b*U3b + 0.5*( u1*u1 + v1*v1 )*U0b ;

      // Right-moving shock
      // Ms = 2.95
      // Xs = 0.
 
      Xs   = 0.0 ;
      Ms   = 2.95 ;
      //
      // Bubble radius and center 
      //
      rb   = 0.2 ;
      xc   = 0.3 ;
      yc   = 0.0 ;
      //
      // right state from analytical jump conditions
      // for right-moving shocks
      //

      rho2 = rho1*( gm + 1.0 )*Ms*Ms/( 2.0 + ( gm - 1.0 )*Ms*Ms ) ;
      p2   = p1*( 2.0*gm*Ms*Ms/( gm + 1.0 ) - ( gm -1.0 )/( gm + 1.0 ) ) ;
      u2   = u1 + 2.0/( gm + 1.0 )*( Ms - 1.0/Ms )*sqrt( gm*p1/rho1 ) ;
      v2   = 0.0 ;
      h2   = gm*p2/( ( gm - 1.0 )*rho2 ) ;

      V0 = gm*p2/( ( gm - 1.0 )*h2 ) ;
      V3 = h2/gm ;

      V1 = V0*u2 ;
      V2 = V0*v2 ;
      V3 = V0*V3 + 0.5*( u2*u2 + v2*v2 )*V0 ;

				
      for ( n = 0 ; n < NN ; n ++ )
	{
	  if ( node[n].coordinate[0] < Xs )
	    {
	      node[n].P[2][0] = V0 ; 	
	      node[n].P[2][1] = V1 ;	
	      node[n].P[2][2] = V2 ; 	
	      node[n].P[2][3] = V3 ;	
	    }
	  else			  
	    {
	      r = ( node[n].coordinate[0] - xc )*( node[n].coordinate[0] - xc ) +
		( node[n].coordinate[1] - yc )*( node[n].coordinate[1] - yc ) ;  
	      r = sqrt( r ) ;

	      if ( r > rb ) {  
		node[n].P[2][0] = U0 ; 	
		node[n].P[2][1] = U1 ;	
		node[n].P[2][2] = U2 ; 	
		node[n].P[2][3] = U3 ;}
	      else {
		node[n].P[2][0] = U0b ;
		node[n].P[2][1] = U1b ;	
		node[n].P[2][2] = U2b ; 	
		node[n].P[2][3] = U3b ;
                                                          
	      }	
	    }
		
	  P_2_Z( node[n].P[2], node[n].Z[2] ) ;
	}                   	
      break ;
    case 11:/*double expansion*/

      // Right state

      rho1 = 1.4 ;
      p1   = 1.0 ;
      u1   = 1.0 ;
      v1   = 0.0 ;
      h1   = gm*p1/( ( gm - 1.0 )*rho1 ) ;
		
      U0 = gm*p1/( ( gm - 1.0 )*h1 ) ;
      U3 = h1/gm ;

      U1 = U0*u1 ;
      U2 = U0*v1 ;
      U3 = U0*U3 + 0.5*( u1*u1 + v1*v1 )*U0 ;

      // Initial position and shock Mach number		    
 
      Xs   = 1.0 ;
      Ms   = 10. ; 

      // Left state from analytical jump conditions
      // for right-moving shocks

      rho2 = 1.4 ;
      p2   = 1. ;
      u2   = -1. ;
      v2   = 0.0 ;
      h2   = gm*p2/( ( gm - 1.0 )*rho2 ) ;

      V0 = gm*p2/( ( gm - 1.0 )*h2 ) ;
      V3 = h2/gm ;

      V1 = V0*u2 ;
      V2 = V0*v2 ;
      V3 = V0*V3 + 0.5*( u2*u2 + v2*v2 )*V0 ;

				
      for ( n = 0 ; n < NN ; n ++ )
	{
	  if ( node[n].coordinate[0] < Xs )
	    {
	      node[n].P[2][0] = V0 ; 	
	      node[n].P[2][1] = V1 ;	
	      node[n].P[2][2] = V2 ; 	
	      node[n].P[2][3] = V3 ;	
	    }
	  else			  
	    {
  
	      node[n].P[2][0] = U0 ; 	
	      node[n].P[2][1] = U1 ;	
	      node[n].P[2][2] = U2 ; 	
	      node[n].P[2][3] = U3 ;	
	    }
		
	  P_2_Z( node[n].P[2], node[n].Z[2] ) ;
	}                   	
      break ;
    case 12:/*Colella&Woodward's blast*/
      for ( n = 0 ; n < NN ; n ++ )
	{
	  if ( node[n].coordinate[0] < 0.1 )
	    {
	      node[n].Z[2][0] = 1. ; 	
	      node[n].Z[2][1] = 0. ;	
	      node[n].Z[2][2] = 0. ; 	
	      node[n].Z[2][3] = 1000. ;	
	    }
	  else {
	    if ( node[n].coordinate[0] < 0.9 )			  
	      {
		node[n].Z[2][0] = 1. ; 	
		node[n].Z[2][1] = 0. ;	
		node[n].Z[2][2] = 0. ; 	
		node[n].Z[2][3] = 0.01 ;	
	      }
	    else
	      {
		node[n].Z[2][0] = 1. ; 	
		node[n].Z[2][1] = 0. ;	
		node[n].Z[2][2] = 0. ; 	
		node[n].Z[2][3] = 100. ;	
	      }
	  }
		
	  Z_2_P( node[n].Z[2], node[n].P[2] ) ;
	}                   	
      break ;
    case 13:/*Flap Deflection*/

      p1 = 1.0 ;
      u1 = 3.0 ;
      v1 = 0.0 ;
      h1 = 2.5 ;
		
      U0 = gm*p1/( ( gm - 1.0 )*h1 ) ;
      U3 = h1/gm ;

      for ( n = 0 ; n < NN ; n ++ )
			                  
	{
	  node[n].P[2][0] = U0 ;
	  node[n].P[2][1] = U0*u1 ;
	  node[n].P[2][2] = U0*v1 ;
	  node[n].P[2][3] = U3*U0 + 0.5*( u1*u1 + v1*v1 )*U0 ;

	  P_2_Z( node[n].P[2], node[n].Z[2] ) ;
	}                   	
      break ;
// ------------------------------------
// ------------------------------------      
    case 14:/* Centered compression */
      
      // Left state
      
      rho1 = 2.6872;
      p1   = 2.57724;
      u1   = 0.734611;
      v1   = 0.0;
      h1   = gm*p1/( ( gm - 1.0 )*rho1 ) ;
      a    = sqrt(gm*p1/rho1) ;
      M1   = u1/a ;
      
      U0 = gm*p1/( ( gm - 1.0 )*h1 ) ;
      U3 = h1/gm ;
      
      U1 = U0*u1 ;
      U2 = U0*v1 ;
      U3 = U0*U3 + 0.5*( u1*u1 + v1*v1 )*U0 ;
      
      // Initial position		    
      
      Xs   = .3 ;
      
      // Right state from analytical jump conditions for steady shock
      
      /* rho2 = 6.0*M1*M1*rho1/(M1*M1+5.0); */
      /* u2   = (rho1*u1)/rho2; */
      /* p2   = (rho1*u1*u1+p1)-(rho2*u2*u2); */
      /* v2   = 0.0; */
      /* h2   = h1+(u1*u1-u2*u2)/2.0 ; */
      rho2 = 1.4;
      u2 = 0.0;
      v2 = 0.0;
      p2 = 1.0;
      h2 = gm*p2/((gm-1.0)*rho2);
      
      V0 = gm*p2/( ( gm - 1.0 )*h2 ) ;
      V3 = h2/gm ;
      
      V1 = V0*u2 ;
      V2 = V0*v2 ;
      V3 = V0*V3 + 0.5*( u2*u2 + v2*v2 )*V0 ;
      
      
      for ( n = 0 ; n < NN ; n ++ )
	{
	  if ( node[n].coordinate[0] > Xs )
	    {
	      node[n].P[2][0] = V0 ; 	
	      node[n].P[2][1] = V1 ;	
	      node[n].P[2][2] = V2 ; 	
	      node[n].P[2][3] = V3 ;	
	    }
	  else			  
	    {
	      
	      node[n].P[2][0] = U0 ; 	
	      node[n].P[2][1] = U1 ;	
	      node[n].P[2][2] = U2 ; 	
	      node[n].P[2][3] = U3 ;	
	    }
	  
	  P_2_Z( node[n].P[2], node[n].Z[2] ) ;
	}                   	
      break ;
 
      case 15: /*NACA test-cases*/ 
     
      for ( n = 0 ; n < NN ; n ++ )
	{
	  p1 = p_inf ;
	  u1 = Mach_inf*a_inf * cos(alpha_inf) ;
	  v1 = Mach_inf*a_inf * sin(alpha_inf) ;
		
	  node[n].P[2][0] = rho_inf ;
	  node[n].P[2][1] = rho_inf * u1 ; 
	  node[n].P[2][2] = rho_inf * v1 ;
	  node[n].P[2][3] = 1.0/(gm-1.0)*p1 + 0.5*( u1*u1 + v1*v1 )*rho_inf ;

	  P_2_Z( node[n].P[2], node[n].Z[2] ) ;
	}
      break ;
    }
}

/***************************************************************************
 *                                                                         *
 *   This program has been developed at the von Karman Institute for Fluid *
 *   Dynamics by Mario Ricchiuto. Any change or update must be done in a   *
 *   clean way and it must be clearly commented.                           *
 *   Mario Ricchiuto                                                       *
 *                                              23/04/2002                 *
 ***************************************************************************/

