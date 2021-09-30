/***************************************************************************
                             LWLF_distribution.c
                               ----------------
This is LWLF_distribution: the distribution of the fluctuation is performed via
the blending  phi_T/3 + theta alpha (U_i-bar U) + (1-theta) Dt K phi_T/2|T|
                             -------------------
    Code developed by M. Ricchiuto, Inria - mario.ricchiuto@inria.fr   
 ***************************************************************************/
#include "common.h"

extern struct element_struct *element ;
extern struct node_struct *node ;

extern double *work_vector, *work_vector0, *work_vector1, *work_vector2 ;
extern double **PHI_d, **PHIN_d, *phi_a_w0, *phi_a_w1, *phi_w, **PHIN_scal ;
extern double ***K_i_p1, ***K_i_p0, **sum_K, **sum_K_1, ***K1 ;
extern double *temp_vector, *normal, *vel, *phi_node, **Left, **Right, **temp_mat ;
extern double **W,  *U_C, dt, **phi_t_w, alpha, **dU ;

extern double diag0, diag1, q_ratio, ubar0, vbar0, speed_of_sound0, shield_factor ;

extern int size, lump_type ;

extern void velocity() ;
extern void initmat( double ** ) ;
extern void add2mat( double **, double ** ) ;
extern void initvec( double * ) ;
extern void add2vec( double *, double * ) ;
extern void A_times_v( double **, double *, double * ) ;
extern void residual_update( int ) ; /* node global index, local residual distributed to node, time level to update (1/2) */
extern void invertmat( double **, double ** ) ;
extern void decompose( double *, double *, double * ) ;
/* direction for the decomposition, nodal N scheme residual, scalar residuals */
extern void limiter( double *, double **, double *, int, int ) ;
/* scalar components of the cell fluctuation, N scheme scalar residuals, limited scalar nodal residuals, target node, nodes involved in the limiting */
extern void recompose( double *, double *, double * ) ;
/* direction for the decomposition, scalar P scheme residuals, P scheme residual */

extern void Eigenvectors( double * ) ;

void LW_LF_unsteady( int e )
{
     int j ,v, n1, k, l ;
     double length, theta, a, b, c, tau ;
     double DDV[3][4], diss[3][4], lump ;


     lump = 1.0*lump_type ;     

//     initmat( sum_K ) ;
//
//     for( v = 0 ; v < 3 ; v++ ) 
//          add2mat( sum_K, K_i_p1[v] ) ;
//
//     invertmat( sum_K, sum_K_1 ) ;
		

     
     for ( j = 0 ; j < size ; j ++ )
         {
// LF dissipation
          diss[0][j]  = 2.*W[0][j] - W[1][j] - W[2][j] ; diss[0][j] *= alpha/3. ;
          diss[1][j]  = 2.*W[1][j] - W[2][j] - W[0][j] ; diss[1][j] *= alpha/3. ; 
          diss[2][j]  = 2.*W[2][j] - W[0][j] - W[1][j] ; diss[2][j] *= alpha/3. ;
//
// Mass matrix corrections
//
          DDV[0][j] = -( 1. - lump )*( 2.*dU[0][j] + dU[1][j] + dU[2][j] )/4. - lump*dU[0][j] ;
          DDV[1][j] = -( 1. - lump )*( 2.*dU[1][j] + dU[0][j] + dU[2][j] )/4. - lump*dU[1][j] ;
          DDV[2][j] = -( 1. - lump )*( 2.*dU[2][j] + dU[1][j] + dU[0][j] )/4. - lump*dU[2][j] ;
	}

     for ( j = 0 ; j < size ; j ++ )
         {
//LF scheme
           PHIN_d[0][j] = dU[0][j] + dt*phi_a_w1[j]/3. + dt*diss[0][j] ;
           PHIN_d[1][j] = dU[1][j] + dt*phi_a_w1[j]/3. + dt*diss[1][j] ;
           PHIN_d[2][j] = dU[2][j] + dt*phi_a_w1[j]/3. + dt*diss[2][j] ;
// Tau matrix
    	   for( l = 0 ; l < size ; l++ )
              {
                for( v = 0 ; v < 3 ; v++ )
                   {
                       K_i_p0[v][j][l] = 0.5*K1[v][j][l]/alpha ;
//                     K_i_p0[v][j][l] = 0.0 ;
//                     for ( k = 0 ; k < size ; k++ ) 
//                           K_i_p0[v][j][l] += 0.5*sum_K_1[j][k]*K1[v][k][l] ;
                   }
              }
         }


/* LW scheme*/


      for ( k = 0 ; k < size ; k ++ )
          {
            PHIN_scal[0][k] = phi_w[k]/3. ;
            PHIN_scal[1][k] = phi_w[k]/3. ;
            PHIN_scal[2][k] = phi_w[k]/3. ;

            for ( l = 0; l < size ; l ++ )
                {
                  PHIN_scal[0][k] += K_i_p0[0][k][l]*phi_w[l] ;
                  PHIN_scal[1][k] += K_i_p0[1][k][l]*phi_w[l] ;
                  PHIN_scal[2][k] += K_i_p0[2][k][l]*phi_w[l] ;
                }
          }


       for ( v = 0 ; v < 3 ; v++ )  
           {
             for ( k = 0 ; k < size ; k ++ )
                 {
                   theta = fabs( PHIN_d[0][k] ) +  fabs( PHIN_d[1][k] ) + fabs( PHIN_d[2][k] ) + 1.e-20 ;
                   theta = fabs( phi_w[k] )/theta ;

                   phi_node[k] = DDV[v][k] + theta*PHIN_d[v][k] + ( 1. - theta )*PHIN_scal[v][k] ;
                 }

              residual_update( element[e].node[v] ) ;
          }
}


