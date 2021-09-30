/***************************************************************************
                               N_distribution.c
                               ----------------
  This is N_distribution: here the distribution of the cell fluctuation is
     performed through the N scheme distribution function KP(U_i-U_c)
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

void LDA_unsteady( int e )
{
     int j ,v, n1, k, l ;
     double length, theta, a, b, c, tau ;
     double DDV[3][4], diss[3][4], lump ;
     double div, vel_x, vel_y, ale_corr_fact ;

     lump = 1.0*lump_type ;     

     initmat( sum_K ) ;

//
// computing the divergence div(sigma)
//
     div = 0. ;
     for ( v = 0 ; v < 3 ; v ++ )
        { 
          vel_x = node[element[e].node[v]].vel[0] ;
          vel_y = node[element[e].node[v]].vel[1] ; 

          div += 0.5* ( vel_x*element[e].normal[v][0] + vel_y*element[e].normal[v][1] ) ;								     
        } 

     for( v = 0 ; v < 3 ; v++ ) 
          add2mat( sum_K, K_i_p1[v] ) ;

     invertmat( sum_K, sum_K_1 ) ;
		
     for ( j = 0 ; j < size ; j ++ )
         {
//
// ALE Mass matrix corrections
//

          ale_corr_fact = ( 1. + dt/(2.*element[e].volume)*div ) ; 

          DDV[0][j] = ale_corr_fact* ( -( 1. - lump )*( 2.*dU[0][j] + dU[1][j] + dU[2][j] )/4. - lump*dU[0][j] ) ;
          DDV[1][j] = ale_corr_fact* ( -( 1. - lump )*( 2.*dU[1][j] + dU[0][j] + dU[2][j] )/4. - lump*dU[1][j] ) ;
          DDV[2][j] = ale_corr_fact* ( -( 1. - lump )*( 2.*dU[2][j] + dU[1][j] + dU[0][j] )/4. - lump*dU[2][j] ) ;

	}

     for ( j = 0 ; j < size ; j ++ )
         {
           work_vector0[j] = 0. ;
    	   for( l = 0 ; l < size ; l++ )
                work_vector0[j] += sum_K_1[j][l]*phi_w[l] ;
         }


/* LDA scheme*/



       for ( v = 0 ; v < 3 ; v++ )  
           {
             for ( k = 0 ; k < size ; k ++ )
                 {
                   phi_node[k] = DDV[v][k] ;

                   for( l = 0 ; l < size ; l++ )
                        phi_node[k] += K_i_p1[v][k][l]*work_vector0[l] ;
                 }

              residual_update( element[e].node[v] ) ;
          }
}


