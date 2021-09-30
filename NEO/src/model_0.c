/***************************************************************************
                                  model_0.c
                                  ---------
  This is model_0: it contains all the model dependent functions used for
        the solution of the 2D EULER equations with PERFECT GAS EOS
                             -------------------
    Code developed by M. Ricchiuto, Inria - mario.ricchiuto@inria.fr
***************************************************************************/
#include "common.h"

extern struct node_struct *node ;
extern struct element_struct *element ;
extern struct boundary_struct  *boundary ;
extern struct numerical_int_struct numerical_int ;
extern struct boundary_nodes_struct *b_nodes ; 

extern int initial_state, size, face_q_pts, volume_q_pts, NBN, NE, NN ;

extern double sigma_barx, sigma_bary ;
extern double *vel, **W, time, dt, **Right, **Left, *Lambda, ***K_i_p1 ;
extern double *phi_a_w0, *phi_a_w1, *phi_node, gm, *phi_w, *temp_vector, *work_vector ;
extern double *work_vector0, *work_vector1, *P_bar, **temp_mat, *FLUX, *temp_normal ;
extern double pressure0, speed_of_sound0, total_enthalpy0, *FLUX ;
extern double pressure1, zero_r, speed_of_sound1, total_enthalpy1 ;
extern double Mach_in, rho_in, alpha_in ;
extern double p_out ;
extern double Mach_inf, rho_inf, alpha_inf, p_inf;

extern void residual_update( int ) ;
/* node global index, local residual distributed to node, time level to update (1/2) */

extern void A_times_v( double **, double *, double * ) ;
extern void A_times_B( double **, double **, double ** ) ;
extern void initvec( double * ) ;
extern void initmat( double ** ) ;
extern void invertmat( double **, double ** ) ;
extern void Eigenvectors( double * ) ;
extern void Waves( double * ) ;
extern void P_2_Z( double *, double * ) ;
extern void Z_2_P( double *, double * ) ;
extern void first_stremline( int ) ;


void ad_flux( double *norm, double *P )
{
  FLUX[0] = P[0]*( P[1]*norm[0] + P[2]*norm[1] ) ;
  FLUX[1] = FLUX[0]*P[1] + P[3]*norm[0] ;
  FLUX[2] = FLUX[0]*P[2] + P[3]*norm[1] ;
  FLUX[3] = FLUX[0]*( gm*P[3]/( P[0]*( gm - 1.0 ) ) + 0.5*( P[1]*P[1] + P[2]*P[2] ) ) ;
}


void P_2_Z( double *P, double *Z )
{
  if ( P[0] < zero_r )
    Z[0] = zero_r ;
  else
    Z[0] = P[0] ;
  if ( P[0] < zero_r )
    {
      Z[1] = 0. ;
      Z[2] = 0. ;
    }
  else {
    Z[1] = P[1]/P[0] ;/* x - velocity */
    Z[2] = P[2]/P[0] ; /* y - velocity */}
  Z[3] = ( gm - 1.0 )*( P[3] - 0.5*P[0]*( Z[1]*Z[1] + Z[2]*Z[2] ) ) ;/*pressure */
  if ( Z[3] < zero_r ) Z[3] = zero_r ;

}

void Z_2_P( double *Z, double *P )
{
  if ( Z[0] < zero_r )
    P[0] = zero_r ;
  else 
    P[0] = Z[0] ;
  if ( Z[0] < zero_r )
    {
      P[1] = 0. ;
      P[2] = 0. ;
    }
  else {
    P[1] = Z[0]*Z[1] ;
    P[2] = Z[0]*Z[2] ; }
  
  if ( Z[3] > zero_r )
    P[3] = Z[3]/( gm - 1.0 ) + 0.5*Z[0]*( Z[1]*Z[1] + Z[2]*Z[2] ) ;
  else 
    P[3] = zero_r/( gm - 1.0 ) + 0.5*Z[0]*( Z[1]*Z[1] + Z[2]*Z[2] ) ;
}

/****************************************/
/**   Eigenvalues for 2D Euler (ALE)   **/
/****************************************/

void Waves( double *n )
{
  double edge ;

  edge = sqrt( n[0]*n[0] + n[1]*n[1] ) ;

  Lambda[0] = ( P_bar[1]*n[0]+P_bar[2]*n[1] )/P_bar[0] - sigma_barx*n[0]-sigma_bary*n[1] ;
  Lambda[1] = Lambda[0] ;
  Lambda[2] = Lambda[0] + speed_of_sound0*edge ;
  Lambda[3] = Lambda[0] - speed_of_sound0*edge ;
}


/****************************************/
/**     Eigenvectors for 2D Euler      **/
/****************************************/

void  Eigenvectors( double *norm )
{
  int m,l ;
  double edge, nx, ny ;
  double rho, u, v, a, p, H, k ;

  edge = sqrt( norm[0]*norm[0] + norm[1]*norm[1] ) ;

  if ( edge < epsilon ) edge = 1.0 ;

  nx = norm[0]/edge ;
  ny = norm[1]/edge ;

  p   = pressure0 ;
  rho = P_bar[0] ;
  a   = speed_of_sound0 ;
  u   = P_bar[1]/rho ;
  v   = P_bar[2]/rho ;
  k   = ( u*u + v*v )*0.5 ;
  H   = total_enthalpy0 ;

  Right[0][0] = 1.0 ;
  Right[1][0] = u ;
  Right[2][0] = v ;
  Right[3][0] = k ;

  Right[0][1] = 0.0 ;
  Right[1][1] = -rho*ny ;
  Right[2][1] = rho*nx ;
  Right[3][1] = rho*( v*nx - u*ny ) ;

  Right[0][2] = rho/a ;
  Right[1][2] = rho*( nx + u/a ) ;
  Right[2][2] = rho*( ny + v/a ) ;
  Right[3][2] = rho*H/a + rho*( u*nx + v*ny )  ;

  Right[0][3] = rho/a ;
  Right[1][3] = rho*( -nx + u/a ) ;
  Right[2][3] = rho*( -ny + v/a ) ;
  Right[3][3] = rho*H/a - rho*( u*nx + v*ny )  ;

  Left[0][0] = 1.0 - 0.4*k/( a*a ) ;
  Left[0][1] = 0.4*u/( a*a ) ;
  Left[0][2] = 0.4*v/( a*a ) ;
  Left[0][3] = -0.4/( a*a ) ;

  Left[1][0] = ( u*ny - v*nx )/rho ;
  Left[1][1] = -ny/rho ;
  Left[1][2] = nx/rho ;
  Left[1][3] = 0.0 ;

  Left[2][0] = 0.5*0.4*k/( rho*a ) - 0.5*( u*nx + v*ny )/rho ;
  Left[2][1] = -0.4*0.5*u/( rho*a ) + 0.5*nx/rho ;
  Left[2][2] = -0.4*0.5*v/( rho*a ) + 0.5*ny/rho ;
  Left[2][3] = 0.5*0.4/( rho*a ) ;

  Left[3][0] = 0.5*0.4*k/( rho*a ) + 0.5*( u*nx + v*ny )/rho ;
  Left[3][1] = -0.4*0.5*u/( rho*a ) - 0.5*nx/rho ;
  Left[3][2] = -0.4*0.5*v/( rho*a ) - 0.5*ny/rho ;
  Left[3][3] = 0.5*0.4/( rho*a ) ;
}


/*****************************************************/
/**          local velocity unit vector             **/
/**       for residual scalar decomposition         **/
/*****************************************************/

void velocity()
{
  double norm, Ma ;

  norm = sqrt( P_bar[1]*P_bar[1] + P_bar[2]*P_bar[2] )/P_bar[0] ;
  Ma = norm/speed_of_sound0 ;

  if ( Ma > 0.001 )
    {
      vel[0] = P_bar[1]/( P_bar[0]*norm ) ;
      vel[1] = P_bar[2]/( P_bar[0]*norm ) ;
    }
  else
    {
      vel[0] = 1.0 ;
      vel[1] = 0.0 ;
    }
}



void outlet_2D_Euler( int m )
{
    int n = b_nodes[m].node ;
    
    double drho, rr  ;
    drho  = node[n].Res[0] ;
    rr = node[n].P[2][0] - drho/node[n].vol_mod ;
    double drhou, rru ;
    drhou = node[n].Res[1] ;
    rru = node[n].P[2][1] - drhou/node[n].vol_mod ;
    double drhov, rrv ;
    drhov = node[n].Res[2] ;
    rrv = node[n].P[2][2] - drhou/node[n].vol_mod ;
    double drhoe, rre, re ;
    drhoe = node[n].Res[3] ;
    rre = node[n].P[2][3] - drhoe/node[n].vol_mod ;
    
    double pp = (gm-1.)*( rre - 0.5*(rru*rru + rrv*rrv)/rr ) ;
    double Ma2 = ( rru*rru + rrv*rrv )/( rr*gm*pp ) ;
    
    double rho = node[n].P[2][0]  ;
    double rhou = node[n].P[2][1]  ;
    double rhov = node[n].P[2][2]  ;
    double rhoe = node[n].P[2][3]  ;
    double ppp = (gm-1.)*( rhoe - 0.5*(rhou*rhou + rhov*rhov)/rho ) ;
    double Ma22 = ( rhou*rhou + rhov*rhov )/( rho*gm*ppp ) ;
    
    
    if (Ma22 <= 1. )
    {
    double Minf = 20. ;
    double Cinf2 = 1./(Minf*Minf) ;
    double htot = Cinf2*( 1. + 0.5*(gm - 1. )*Minf*Minf ) ;
    
    re = rr*( htot + 0.5*( gm - 1. )*( rru*rru + rrv*rrv )/(rr*rr) )/gm ;
    
    node[n].Res[3] = ( node[n].P[2][3] - re )*node[n].vol_mod ;
    }
    
    
}


/*****************************************************/
/**                  Streamline BC                  **/
/*****************************************************/

void first_streamline( int  m )
{
  int i, j, k, f, n, n1, n2 ;
  double nx, ny, length ;
  double rho, u, v, p, rhoun, H ;

  n1 = boundary[m].node[0] ;
  n2 = boundary[m].node[1] ;

  nx = boundary[m].normal[0] ;
  ny = boundary[m].normal[1] ;

  length = sqrt( nx*nx + ny*ny ) ;  
  
// new BC
    
  for (k = 0; k < face_q_pts; k++ )
    {
      rho = numerical_int.face_coordinate[k][0]*node[n1].Z[0][0] +
            numerical_int.face_coordinate[k][1]*node[n2].Z[0][0] ;

      u   = numerical_int.face_coordinate[k][0]*node[n1].Z[0][1] +
            numerical_int.face_coordinate[k][1]*node[n2].Z[0][1] ;

      v   = numerical_int.face_coordinate[k][0]*node[n1].Z[0][2] +
            numerical_int.face_coordinate[k][1]*node[n2].Z[0][2] ;

      p   = numerical_int.face_coordinate[k][0]*node[n1].Z[0][3] +
            numerical_int.face_coordinate[k][1]*node[n2].Z[0][3] ;
      
      rhoun = rho * ( u * nx + v * ny ) ;
      H     = gm/(gm-1.0) * p/rho + 0.5* ( u*u +  v*v ) ;
/*
   if (n1 == 58 || n2 == 58)  
     {
	printf("rho,u,v,p = %le,%le,%le,%le\n",rho,u,v,p) ;
	printf("rhoun = %le, H = %le\n",rhoun,H) ;
     }
*/
      node[n1].Res[0] += dt * numerical_int.face_weight[k]*numerical_int.face_coordinate[k][0] * rhoun ; 
      node[n1].Res[1] += dt * numerical_int.face_weight[k]*numerical_int.face_coordinate[k][0] * rhoun*u ;
      node[n1].Res[2] += dt * numerical_int.face_weight[k]*numerical_int.face_coordinate[k][0] * rhoun*v ;
      node[n1].Res[3] += dt * numerical_int.face_weight[k]*numerical_int.face_coordinate[k][0] * rhoun*H ;

      node[n2].Res[0] += dt * numerical_int.face_weight[k]*numerical_int.face_coordinate[k][1] * rhoun ; 
      node[n2].Res[1] += dt * numerical_int.face_weight[k]*numerical_int.face_coordinate[k][1] * rhoun*u ;
      node[n2].Res[2] += dt * numerical_int.face_weight[k]*numerical_int.face_coordinate[k][1] * rhoun*v ;
      node[n2].Res[3] += dt * numerical_int.face_weight[k]*numerical_int.face_coordinate[k][1] * rhoun*H ;
    }
	
//   if (n1 == 58 && n2 == 59)  
//     {
//	printf("rho,u,v,p = %le,%le,%le,%le\n",rho,u,v,p) ;
//	printf("rhoun = %le, H = %le\n",rhoun,H) ;
//     }
//
//  
//   if (n1 == 58 && n2 == 59)  
//     {
//      printf("nx = %f, ny = %f\n",nx/length,ny/length) ;
//      printf("n1 = %d, Res = %le,%le,%le,%le\n",n1,node[n1].Res[0],node[n1].Res[1],node[n1].Res[2],node[n1].Res[3]) ; 
//      printf("n2 = %d, Res = %le,%le,%le,%le\n",n2,node[n2].Res[0],node[n2].Res[1],node[n2].Res[2],node[n2].Res[3]) ; 
//     }
//   if (n1 == 76 && n2 == 77) 
//     {
//      printf("nx = %f, ny = %f\n",nx/length,ny/length) ;
//      printf("n1 = %d, Res = %le,%le,%le,%le\n",n1,node[n1].Res[0],node[n1].Res[1],node[n1].Res[2],node[n1].Res[3]) ; 
//      printf("n2 = %d, Res = %le,%le,%le,%le\n",n2,node[n2].Res[0],node[n2].Res[1],node[n2].Res[2],node[n2].Res[3]) ; 
//     }
  

// old BC
/*
  double phin, un ;

  phin = node[n1].Res[1]*nx + node[n1].Res[2]*ny ;

  un = node[n1].P[0][1]*nx + node[n1].P[0][2]*ny ;

  node[n1].Res[1] -= (phin-un*node[n1].vol_mod)*nx ;
  node[n1].Res[2] -= (phin-un*node[n1].vol_mod)*ny ;

  phin = node[n2].Res[1]*nx + node[n2].Res[2]*ny ;

  un = node[n2].P[0][1]*nx + node[n2].P[0][2]*ny ;

  node[n2].Res[1] -= (phin-un*node[n2].vol_mod)*nx ;
  node[n2].Res[2] -= (phin-un*node[n2].vol_mod)*ny ;
*/
}


void inlet_2D_Euler( int m )
{ 
  int i, j, k, f, n, n1, n2 ;
  double length, nx, ny, u, v, H, rho, p ;

  n1 = boundary[m].node[0] ;
  n2 = boundary[m].node[1] ;

  nx = boundary[m].normal[0] ;
  ny = boundary[m].normal[1] ;
  
  length = sqrt(nx*nx + ny*ny) ;
 
  temp_normal[0] = nx / sqrt(nx*nx + ny*ny) ;
  temp_normal[1] = ny / sqrt(nx*nx + ny*ny) ;

  for (k = 0; k < face_q_pts; k++ )
{
      rho = numerical_int.face_coordinate[k][0]*node[n1].Z[0][0] + 
            numerical_int.face_coordinate[k][1]*node[n2].Z[0][0] ;
      
      u   = numerical_int.face_coordinate[k][0]*node[n1].Z[0][1] + 
            numerical_int.face_coordinate[k][1]*node[n2].Z[0][1] ;
      
      v   = numerical_int.face_coordinate[k][0]*node[n1].Z[0][2] + 
            numerical_int.face_coordinate[k][1]*node[n2].Z[0][2] ;
      
      p   = numerical_int.face_coordinate[k][0]*node[n1].Z[0][3] + 
            numerical_int.face_coordinate[k][1]*node[n2].Z[0][3] ;
      
      temp_vector[0] = rho ;  
      temp_vector[1] = u ;  
      temp_vector[2] = v ;  
      temp_vector[3] = p ;
       
      Z_2_P( temp_vector, work_vector0 ) ;

      rho = rho_in ;
      u   = Mach_in * sqrt(gm*p/rho) * cos(alpha_in) ;     
      v   = Mach_in * sqrt(gm*p/rho) * sin(alpha_in) ;     
      
      temp_vector[0] = rho ;  
      temp_vector[1] = u ;  
      temp_vector[2] = v ;  
      temp_vector[3] = p ;
       
      pressure0       = p ;
      speed_of_sound0 = sqrt(gm * p/rho) ;
      total_enthalpy0 = gm/(gm-1.0) * p/rho + 0.5*(u*u + v*v) ;
 
      Z_2_P( temp_vector, P_bar ) ;

      for (i = 0; i < 4; i++)
	work_vector0[i] -= P_bar[i] ;
       
      Eigenvectors( temp_normal ) ;
      Waves( temp_normal ) ;
      
      for (i = 0; i < 4; i++)
        {
	work_vector1[i] = 0.0 ;
      
        if (Lambda[i] < 0.0)
	   for (j = 0; j < 4; j++)
	       work_vector1[i] += Lambda[i] * Left[i][j] * work_vector0[j] ;
	  
	}

      for (i = 0; i < 4; i++)
      { 
          work_vector0[i] = 0.0 ;
          for (j = 0; j < 4; j++)
              work_vector0[i] += Right[i][j] * work_vector1[j] ;
      }	
     
      node[n1].Res[0] += dt * length * numerical_int.face_weight[k]*numerical_int.face_coordinate[k][0] * work_vector0[0] ;
      node[n1].Res[1] += dt * length * numerical_int.face_weight[k]*numerical_int.face_coordinate[k][0] * work_vector0[1] ;
      node[n1].Res[2] += dt * length * numerical_int.face_weight[k]*numerical_int.face_coordinate[k][0] * work_vector0[2] ;
      node[n1].Res[3] += dt * length * numerical_int.face_weight[k]*numerical_int.face_coordinate[k][0] * work_vector0[3] ;

      node[n2].Res[0] += dt * length * numerical_int.face_weight[k]*numerical_int.face_coordinate[k][1] * work_vector0[0] ;
      node[n2].Res[1] += dt * length * numerical_int.face_weight[k]*numerical_int.face_coordinate[k][1] * work_vector0[1] ;
      node[n2].Res[2] += dt * length * numerical_int.face_weight[k]*numerical_int.face_coordinate[k][1] * work_vector0[2] ;
      node[n2].Res[3] += dt * length * numerical_int.face_weight[k]*numerical_int.face_coordinate[k][1] * work_vector0[3] ;
}  

}

/*****************************************************/
/**              BCs for the 2D RP .....            **/
/*****************************************************/

void supersonic_residual_2DRP_Euler( int m )
{
  int  i, n ;
  double speedx_b, speedx_t ;
  double speedy_l, speedy_r ;
  double x, y, Xo, Yo ;
  double p1, u1, v1, h1 ;
  double U0, U3 ;

  n = b_nodes[m].node ;

  x = node[n].coordinate[0] ;
  y = node[n].coordinate[1] ;

  speedx_t = -0.5322581*1.2060454/( 1.5 - 0.5322581 ) ;
  speedx_b = 0.1379928*1.2060454/( 0.1379928 - 0.5322581 ) ;
  speedy_l = -0.1379928*1.2060454/( 0.5322581 - 0.1379928 ) ;
  speedy_r = -0.5322581*1.2060454/( 1.5 - 0.5322581 ) ;

  if ( y > 0.999 )
    {
      node[n].Res[2] = 0. ;
      node[n].P[0][2] = 0. ;
    }

  if ( y < 0.001 )
    {
      for ( i = 0 ; i < size ; i++ ) node[n].Res[i] = 0 ;

      Xo = 0.8 - fabs( speedx_b )*time ;

      if ( x >= Xo )
	{
	  p1 = 0.3 ;
	  u1 = 0.0 ;
	  v1 = 1.2060454 ;
	  h1 = 1.9727271 ;

	  U0 = gm*p1/( ( gm - 1.0 )*h1 ) ;
	  U3 = U0*h1/gm ;

	  node[n].P[0][0] = U0 ;
	  node[n].P[0][1] = U0*u1 ;
	  node[n].P[0][2] = U0*v1 ;
	  node[n].P[0][3] = U3 + U0*0.5*( u1*u1 + v1*v1 ) ;
	}
      else
	{
	  p1 = 0.0290323 ;
	  u1 = 1.2060454 ;
	  v1 = 1.2060454 ;
	  h1 = 0.73636487 ;

	  U0 = gm*p1/( ( gm - 1.0 )*h1 ) ;
	  U3 = U0*h1/gm ;
		
	  node[n].P[0][0] = U0 ;
	  node[n].P[0][1] = U0*u1 ;
	  node[n].P[0][2] = U0*v1 ;
	  node[n].P[0][3] = U3 + U0*0.5*( u1*u1 + v1*v1 ) ;
	}
    }

  if ( x >  0.999 )
    {
      node[n].Res[1] = 0. ;
      node[n].P[0][1] = 0. ;
    } 

  if ( x < 0.001 )
    {
      for ( i = 0 ; i < size ; i++ ) node[n].Res[i] = 0 ;
      
      Yo = 0.8 - fabs( speedy_l )*time ;

      if ( y >= Yo )
	{
	  p1 = 0.3 ;
	  u1 = 1.2060454 ;
	  v1 = 0.0 ;
	  h1 = 1.9727271 ;

	  U0 = gm*p1/( ( gm - 1.0 )*h1 ) ;
	  U3 = U0*h1/gm ;

	  node[n].P[0][0] = U0 ;
	  node[n].P[0][1] = U0*u1 ;
	  node[n].P[0][2] = U0*v1 ;
	  node[n].P[0][3] = U3 + U0*0.5*( u1*u1 + v1*v1 ) ;
	}
      else
	{
	  p1 = 0.0290323 ;
	  u1 = 1.2060454 ;
	  v1 = 1.2060454 ;
	  h1 = 0.73636487 ;

	  U0 = gm*p1/( ( gm - 1.0 )*h1 ) ;
	  U3 = U0*h1/gm ;
	
	  node[n].P[0][0] = U0 ;
	  node[n].P[0][1] = U0*u1 ;
	  node[n].P[0][2] = U0*v1 ;
	  node[n].P[0][3] = U3 + U0*0.5*( u1*u1 + v1*v1 ) ;
	}
    }

  P_2_Z( node[n].P[2], node[n].Z[2] ) ;
}

void supersonic_residual_2DRPiso_Euler( int m )
{
  int i, n ;
  double speedL, speedR ;
  double X2, X1 ;
  double p1, u1, v1, h1, x, y ;
  double U0, U3, phin, nx, ny, length ;


  n = b_nodes[m].node ;

  x = node[n].coordinate[0] ;
  y = node[n].coordinate[1] ;

  speedL = 0.1379928*1.2060454/( 0.5322581 - 0.1379928 ) ;
  speedR = 0.5322581*1.2060454/( 1.5 - 0.5322581 ) ;

  X2 = sqrt( 2.0 )*( 0.9 - 0.5*speedR*time ) ;
  X1 = sqrt( 2.0 )*( 0.4 - 0.5*speedL*time ) ;

  nx = b_nodes[m].normal[0] ;
  ny = b_nodes[m].normal[1] ;
 
  length = sqrt( nx*nx + ny*ny ) ;  

  nx /= length ;
  ny /= length ;

  if ( fabs( x ) < 0.0001 ){ nx = 0. ; ny = 1. ; }
  if ( fabs( x - sqrt( 2. ) ) < 0.0001 ){ nx = 0. ; ny = 1. ; }

  if ( y < 0.0000000001 ) {
    node[n].Res[2] = 0. ;
  }
  else
    {
      if ( x > sqrt( 2. )/2. )
        {
          nx = -0.5*sqrt( 2. ) ;
          ny = -0.5*sqrt( 2. ) ;

          phin = node[n].Res[1]*nx + node[n].Res[2]*ny ;

          node[n].Res[1] -= phin*nx ;
          node[n].Res[2] -= phin*ny ;
        }
      else
        {
	  for ( i = 0 ; i < size ; i++ ) 
	    node[n].Res[i] = 0 ;
          if ( x >= X1 )
	    {
	      p1 = 0.3 ;
	      u1 = 1.2060454*0.5*sqrt( 2.0 ) ;
	      v1 = -1.2060454*0.5*sqrt( 2.0 ) ;
	      h1 = 1.9727271 ;

	      U0 = gm*p1/( ( gm - 1.0 )*h1 ) ;
	      U3 = U0*h1/gm ;

	      node[n].P[2][0] = U0 ;
	      node[n].P[2][1] = U0*u1 ;
	      node[n].P[2][2] = U0*v1 ;
	      node[n].P[2][3] = U3 + U0*0.5*( u1*u1 + v1*v1 ) ;
	    }
          else
	    {
	      p1 = 0.0290323 ;
	      u1 = 1.2060454*sqrt( 2.0 ) ;
	      v1 = 0.0 ;
	      h1 = 0.73636487 ;

	      U0 = gm*p1/( ( gm - 1.0 )*h1 ) ;
	      U3 = U0*h1/gm ;

	      node[n].P[2][0] = U0 ;
	      node[n].P[2][1] = U0*u1 ;
	      node[n].P[2][2] = U0*v1 ;
	      node[n].P[2][3] = U3 + U0*0.5*( u1*u1 + v1*v1 ) ;
            }
	  for ( i = 0 ; i < size ; i ++ ) 
	    node[n].Res[i] = -node[n].vol*( node[n].P[2][i] - node[n].P[0][i] ) ;
        }
    }

}

/*****************************************************/
/**        BCs for the Double Mach Reflection       **/
/*****************************************************/

void supersonic_residual_dMach_Euler( int m )
{
  int i, n ;
  double x, y, Yo ;
  double p1, u1, v1, h1 ;
  double p2, u2, v2, h2 ;
  double U0, U1, U2, U3 ;
  double V0, V1, V2, V3 ;

  n = b_nodes[m].node ; 

  x = node[n].coordinate[0] ;
  y = node[n].coordinate[1] ;

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

  Yo = 1.0/6.0  + ( y + 20.0* ( time ) )/sqrt(3.0) ;
  /*MARIO: RESTARTING FROM THE LAST SOLUTION!!!!!!!!*/
  if ( y < 0.00001 )
    {
      if ( x >= 1.0/6.0 )
	node[n].Res[2] = 0.0 ;
      else for ( i = 0 ; i < size ; i++ ) node[n].Res[i] = 0 ;
    }
  else
    {
      for ( i = 0 ; i < size ; i++ ) node[n].Res[i] = 0 ;

      if ( x < Yo )
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
      for ( i = 0 ; i < size ; i ++ ) 
	node[n].Res[i] = -node[n].vol*( node[n].P[2][i] - node[n].P[0][i] ) ;
    }


  P_2_Z( node[n].P[2], node[n].Z[2] ) ;
}


void far_2D_Euler( int m )
{
  int i, j, k, f, n, n1, n2 ;
  double length, nx, ny, u, v, H, rho, p ;

  n1 = boundary[m].node[0] ;
  n2 = boundary[m].node[1] ;

  nx = boundary[m].normal[0] ;
  ny = boundary[m].normal[1] ;
 
  length = sqrt(nx*nx + ny*ny) ;
 
  temp_normal[0] = nx / length ;
  temp_normal[1] = ny / length ;

  for (k = 0; k < face_q_pts; k++ )
{
      rho = numerical_int.face_coordinate[k][0]*node[n1].Z[0][0] + 
            numerical_int.face_coordinate[k][1]*node[n2].Z[0][0] ;
      
      u   = numerical_int.face_coordinate[k][0]*node[n1].Z[0][1] + 
            numerical_int.face_coordinate[k][1]*node[n2].Z[0][1] ;
      
      v   = numerical_int.face_coordinate[k][0]*node[n1].Z[0][2] + 
            numerical_int.face_coordinate[k][1]*node[n2].Z[0][2] ;
      
      p   = numerical_int.face_coordinate[k][0]*node[n1].Z[0][3] + 
            numerical_int.face_coordinate[k][1]*node[n2].Z[0][3] ;
      
      temp_vector[0] = rho ;  
      temp_vector[1] = u ;  
      temp_vector[2] = v ;  
      temp_vector[3] = p ;
       
      Z_2_P( temp_vector, work_vector0 ) ;

      rho = rho_inf ;
      p   = p_inf ;
      u   = Mach_inf * sqrt(gm*p/rho) * cos(alpha_inf) ;     
      v   = Mach_inf * sqrt(gm*p/rho) * sin(alpha_inf) ;     
 
      temp_vector[0] = rho ;  
      temp_vector[1] = u ;  
      temp_vector[2] = v ;  
      temp_vector[3] = p ;
       
      pressure0       = p ;
      speed_of_sound0 = sqrt(gm * p/rho) ;
      total_enthalpy0 = gm/(gm-1.0) * p/rho + 0.5*(u*u + v*v) ;
 
      Z_2_P( temp_vector, P_bar ) ;

      for (i = 0; i < 4; i++)
	work_vector0[i] -= P_bar[i] ; /* Ub - Uinf */
       
      Eigenvectors( temp_normal ) ;
      Waves( temp_normal ) ;
     
      for (i = 0; i < 4; i++)
        {
	work_vector1[i] = 0.0 ;
      
        if (Lambda[i] > 0.0)
	   for (j = 0; j < 4; j++)
	       work_vector1[i] += Lambda[i] * Left[i][j] * work_vector0[j] ;
	  
	}

      for (i = 0; i < 4; i++)
      { 
          work_vector0[i] = 0.0 ;
          for (j = 0; j < 4; j++)
              work_vector0[i] += Right[i][j] * work_vector1[j] ;
      }	
     
      node[n1].Res[0] += dt * length * numerical_int.face_weight[k]*numerical_int.face_coordinate[k][0] * work_vector0[0] ;
      node[n1].Res[1] += dt * length * numerical_int.face_weight[k]*numerical_int.face_coordinate[k][0] * work_vector0[1] ;
      node[n1].Res[2] += dt * length * numerical_int.face_weight[k]*numerical_int.face_coordinate[k][0] * work_vector0[2] ;
      node[n1].Res[3] += dt * length * numerical_int.face_weight[k]*numerical_int.face_coordinate[k][0] * work_vector0[3] ;

      node[n2].Res[0] += dt * length * numerical_int.face_weight[k]*numerical_int.face_coordinate[k][1] * work_vector0[0] ;
      node[n2].Res[1] += dt * length * numerical_int.face_weight[k]*numerical_int.face_coordinate[k][1] * work_vector0[1] ;
      node[n2].Res[2] += dt * length * numerical_int.face_weight[k]*numerical_int.face_coordinate[k][1] * work_vector0[2] ;
      node[n2].Res[3] += dt * length * numerical_int.face_weight[k]*numerical_int.face_coordinate[k][1] * work_vector0[3] ;
}  
}


void out_2D_Euler( int m )
{
  int i, j, k, f, n, n1, n2 ;
  double length, nx, ny, rhoun, u, v, H, rho, p ;

  n1 = boundary[m].node[0] ;
  n2 = boundary[m].node[1] ;

  nx = boundary[m].normal[0] ;
  ny = boundary[m].normal[1] ;
 
  length = sqrt(nx*nx + ny*ny) ;
 
  temp_normal[0] = nx / length ;
  temp_normal[1] = ny / length ;

  for (k = 0; k < face_q_pts; k++ )
{
      rho = numerical_int.face_coordinate[k][0]*node[n1].Z[0][0] + 
            numerical_int.face_coordinate[k][1]*node[n2].Z[0][0] ;
      
      u   = numerical_int.face_coordinate[k][0]*node[n1].Z[0][1] + 
            numerical_int.face_coordinate[k][1]*node[n2].Z[0][1] ;
      
      v   = numerical_int.face_coordinate[k][0]*node[n1].Z[0][2] + 
            numerical_int.face_coordinate[k][1]*node[n2].Z[0][2] ;
      
      p   = numerical_int.face_coordinate[k][0]*node[n1].Z[0][3] + 
            numerical_int.face_coordinate[k][1]*node[n2].Z[0][3] ;
      
      temp_vector[0] = rho ;  
      temp_vector[1] = u ;  
      temp_vector[2] = v ;  
      temp_vector[3] = p ;
       
      Z_2_P( temp_vector, work_vector0 ) ;

      p   = p_out ;
 
      temp_vector[0] = rho ;  
      temp_vector[1] = u ;  
      temp_vector[2] = v ;  
      temp_vector[3] = p ;
       
      pressure0       = p ;
      speed_of_sound0 = sqrt(gm * p/rho) ;
      total_enthalpy0 = gm/(gm-1.0) * p/rho + 0.5*(u*u + v*v) ;
 
      Z_2_P( temp_vector, P_bar ) ;

      for (i = 0; i < 4; i++)
	work_vector0[i] -= P_bar[i] ;
       
      Eigenvectors( temp_normal ) ;
      Waves( temp_normal ) ;
      
      for (i = 0; i < 4; i++)
        {
	work_vector1[i] = 0.0 ;
      
        if (Lambda[i] > 0.0)
	   for (j = 0; j < 4; j++)
	       work_vector1[i] += Lambda[i] * Left[i][j] * work_vector0[j] ;
	  
	}

      for (i = 0; i < 4; i++)
      { 
          work_vector0[i] = 0.0 ;
          for (j = 0; j < 4; j++)
              work_vector0[i] += Right[i][j] * work_vector1[j] ;
      }	
     
      node[n1].Res[0] += dt * length * numerical_int.face_weight[k]*numerical_int.face_coordinate[k][0] * work_vector0[0] ;
      node[n1].Res[1] += dt * length * numerical_int.face_weight[k]*numerical_int.face_coordinate[k][0] * work_vector0[1] ;
      node[n1].Res[2] += dt * length * numerical_int.face_weight[k]*numerical_int.face_coordinate[k][0] * work_vector0[2] ;
      node[n1].Res[3] += dt * length * numerical_int.face_weight[k]*numerical_int.face_coordinate[k][0] * work_vector0[3] ;

      node[n2].Res[0] += dt * length * numerical_int.face_weight[k]*numerical_int.face_coordinate[k][1] * work_vector0[0] ;
      node[n2].Res[1] += dt * length * numerical_int.face_weight[k]*numerical_int.face_coordinate[k][1] * work_vector0[1] ;
      node[n2].Res[2] += dt * length * numerical_int.face_weight[k]*numerical_int.face_coordinate[k][1] * work_vector0[2] ;
      node[n2].Res[3] += dt * length * numerical_int.face_weight[k]*numerical_int.face_coordinate[k][1] * work_vector0[3] ;
}  
}


