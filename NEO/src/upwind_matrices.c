/***************************************************************************
                              upwind_matrices.c
                              -----------------
 This is upwind_matrices: the model-dependent computation of the positive
   upwind perameters K_PLUS fully using the eigenspectrum of the system
                              is performed here
                              -------------------
    Code developed by M. Ricchiuto, Inria - mario.ricchiuto@inria.fr
 ***************************************************************************/
#include "common.h"

extern struct element_struct *element ;
extern struct node_struct *node ;

extern double sigma_barx, sigma_bary ;
extern double zero_r, dt, **Right, **Left, *Lambda, **temp_mat, alpha, speed_of_sound0, total_enthalpy0 ;
extern double ***K_i_p0, ***K_i_p1, *temp_normal, gm, c_tau, ***K1, ref_vel, ubar0, vbar0 ;

extern int size, scheme ;

extern void Eigenvectors( double * ) ;       /* normals, variables, RIGHT EIGENVECTORS */
extern void Waves( double * ) ;        /* normals, variables, WAVES */
extern void A_times_B( double **, double **, double ** ) ;

void compute_upwind_Keys( int e )
{
     int i, j, vv, kk, n ;
     double  abss, volume, max_v, max_c, max_h, hh ;
     double H, un, u, v, k, nx, ny, gm1, twogm ;
     double sigman ;

    hh = 1.e-3*speed_of_sound0 ; 
    
    H     = total_enthalpy0 ;
    u     = ubar0 ;
    v     = vbar0 ;
    k     = 0.5*( u*u + v*v ) ;        
    gm1   = gm - 1. ;
    twogm = 2. - gm ;

     max_h = max_v = max_c = 0.0 ;

     for ( vv = 0 ; vv < 3 ; vv ++ )
         {
           temp_normal[0] = element[e].normal[vv][0]/2.0 ;
           temp_normal[1] = element[e].normal[vv][1]/2.0 ;

           Eigenvectors( temp_normal ) ;
           Waves( temp_normal ) ;
          
           sigman = temp_normal[0]*sigma_barx + temp_normal[1]*sigma_bary ;

           nx = temp_normal[0] ;
           ny = temp_normal[1] ;
           un = u*nx + v*ny ;

/*********************************************************/
/**           Assembling K1 : Advective part            **/
/*********************************************************/

           K1[vv][0][0] = 0. ;
           K1[vv][0][1] = nx ;
           K1[vv][0][2] = ny ;
           K1[vv][0][3] = 0. ;
     
           K1[vv][1][0] = gm1*k*nx - u*un ;
           K1[vv][1][1] = un + twogm*u*nx ;
           K1[vv][1][2] = u*ny - gm1*v*nx ;
           K1[vv][1][3] = gm1*nx ;
     
           K1[vv][2][0] = gm1*k*ny - v*un ;
           K1[vv][2][1] = v*nx - gm1*u*ny ;
           K1[vv][2][2] = un + twogm*v*ny ;
           K1[vv][2][3] = gm1*ny ;
     
           K1[vv][3][0] = ( gm1*k - H )*un ;
           K1[vv][3][1] = H*nx - gm1*u*un ;
           K1[vv][3][2] = H*ny - gm1*v*un ;
           K1[vv][3][3] = gm*un ;

           if ( scheme != 10 ) 
              {
                for ( i = 0 ; i < size ; i ++)
                    {

/*********************************************************/
/**              Assembling K1 : ALE part               **/
/*********************************************************/

		      K1[vv][i][i] -=  sigman ;

/*********************************************************/
/**                   Building K1 plus                  **/
/*********************************************************/

	              for ( j = 0 ; j < size ; j ++ )
                          {
			    K_i_p1[vv][i][j] = 0.5*K1[vv][i][j] ; 

                            for ( kk = 0 ; kk < size ; kk ++ )
                                {
			          abss = fabs( Lambda[kk] ) ;
			          if( abss < 2.*hh ) 
                                      abss = hh + 0.25*( abss*abss )/hh ;

                                  K_i_p1[vv][i][j] += 0.5*Right[i][kk]*abss*Left[kk][j] ;
			        }
                           }
                      }
               }
           }

// DA MODIFICARE: QUESTO è IL PARAMETRO DI LxF
      for ( vv = 0 ; vv < 3 ; vv ++ )
          {
            for ( j = 0 ; j < 2 ; j ++ )
                  temp_normal[j] = element[e].normal[vv][0] ;
 
            volume = sqrt( temp_normal[0]*temp_normal[0] + temp_normal[1]*temp_normal[1] ) ;
            if ( volume > max_h ) max_h = volume ;

            n = element[e].node[vv] ;

            volume = sqrt( (node[n].Z[2][1]-sigma_barx)*(node[n].Z[2][1]-sigma_barx) + 
                                                         + (node[n].Z[2][2]-sigma_bary)*(node[n].Z[2][2]-sigma_bary) ) ;
            if ( volume > max_v ) max_v = volume ;

            volume = node[n].Z[2][0] ; 
            if ( volume < zero_r ) volume = zero_r ;
 
            volume = sqrt( gm*node[n].Z[2][3]/volume ) ; 
            if ( volume > max_c ) max_c = volume ;
          }


     alpha = 0.5*( max_v + max_c )*max_h ;
}
