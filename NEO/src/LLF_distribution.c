/***************************************************************************
                              LLF_distribution.c
                               ----------------
This is LLF_distribution:  distribution performed via the limited nonlinear
function  Psi(B_LF) phi_T with B_LF phi_T = phi_T/3 + alpha (U_i - bar U)
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
extern double **W,  **Z, *U_C, dt, **phi_t_w, alpha, **dU  ;

extern double diag0, diag1, q_ratio, ubar0, vbar0, speed_of_sound0, shield_factor ;

extern int size, lump_type , iteration ;

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

void LLF_unsteady( int e )
{
     int j ,v, n1, k, l ;
     double volume, betaP, arg, r_max, u_max ;
     double length, theta, a, b, c, tau ;
     double DDV[3][4], diss[3][4], lump, mariuz ;
     double p_max, p_min, coeff ;
     double div, vel_x, vel_y, ale_corr_fact, length2 ;

     lump = 1.0*lump_type ;     

     volume = element[e].volume/3.0 ;

     if ( iteration > 1 )
     {
       length = dt*element[e].volume ;
       length = pow( length, 1./3. ) ;
     }
     else 
    length = sqrt(element[e].volume) ;
    
     initmat( sum_K ) ; // for Bc

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
          add2mat( sum_K, K_i_p1[v] ) ; // for Bc

     invertmat( sum_K, sum_K_1 ) ; //for Bc
		

     
     for ( j = 0 ; j < size ; j ++ )
         {
// LF dissipation
          diss[0][j]  = 2.*W[0][j] - W[1][j] - W[2][j] ; diss[0][j] *= alpha/3. ;
          diss[1][j]  = 2.*W[1][j] - W[2][j] - W[0][j] ; diss[1][j] *= alpha/3. ; 
          diss[2][j]  = 2.*W[2][j] - W[0][j] - W[1][j] ; diss[2][j] *= alpha/3. ;
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
//LF scheme
           PHIN_d[0][j] = dU[0][j] + dt*phi_a_w1[j]/3. + dt*diss[0][j] ;
           PHIN_d[1][j] = dU[1][j] + dt*phi_a_w1[j]/3. + dt*diss[1][j] ;
           PHIN_d[2][j] = dU[2][j] + dt*phi_a_w1[j]/3. + dt*diss[2][j] ;
// Tau matrix
   	   for( l = 0 ; l < size ; l++ )
              {
                for( v = 0 ; v < 3 ; v++ )
                   {
                     K_i_p0[v][j][l] = 0.0 ;
//                     K_i_p0[v][j][l] = 0.5*K1[v][j][l]/alpha ;
                     for ( k = 0 ; k < size ; k++ )
                           K_i_p0[v][j][l] += 0.5*sum_K_1[j][k]*K1[v][k][l] ;
                   }
              }
         }


              

     velocity() ;

     Eigenvectors( vel ) ;

/* limiting and redistribution of the scalar residuals */


      for ( k = 0 ; k < size ; k ++ )
          {
            work_vector0[k] = 0. ;
            PHIN_scal[0][k] = 0. ;
            PHIN_scal[1][k] = 0. ;
            PHIN_scal[2][k] = 0. ;

            for ( l = 0; l < size ; l ++ )
                {
                  work_vector0[k] += Left[k][l]*phi_w[l] ;
                  PHIN_scal[0][k] += Left[k][l]*PHIN_d[0][l] ;
                  PHIN_scal[1][k] += Left[k][l]*PHIN_d[1][l] ;
                  PHIN_scal[2][k] += Left[k][l]*PHIN_d[2][l] ;
                }
          }

     if (  iteration > 1 )
          theta = length*length*length/( fabs( work_vector0[0] )  + 1.e-20 ) ;
     else 
          theta = length*length/( fabs( work_vector0[0] )  + 1.e-20 ) ;

     if ( theta >= 1. ) 
          theta = 1.0 ;

        mariuz = 1.e-15 ;

     for ( v = 0 ; v < 3 ; v ++ ) 
          node[element[e].node[v]].tt += (1./3.)*theta*element[e].volume/node[element[e].node[v]].vol ;

       for ( v = 0 ; v < 3 ; v ++ )
           {
             for ( j = 0 ; j < size ; j ++ )
                 {
                   if ( fabs( work_vector0[j] ) > 1.e-20 )
                      {
                        arg = PHIN_scal[0][j]/work_vector0[j] ;
                        if ( arg > 0. ) betaP  = arg + mariuz ;
                        else betaP = mariuz ;

                       arg = PHIN_scal[1][j]/work_vector0[j] ;
                       if ( arg > 0. ) betaP  += arg + mariuz ;
                       else betaP += mariuz ;

                       arg = PHIN_scal[2][j]/work_vector0[j] ;
                       if ( arg > 0. ) betaP  += arg + mariuz ;
                       else betaP += mariuz ;

                       arg    = PHIN_scal[v][j]/work_vector0[j] ;
                       if ( arg > 0. ) betaP  = ( arg + mariuz )/betaP ;
                       else betaP = mariuz/betaP ;

                       work_vector1[j] = betaP*work_vector0[j] ;
                     }
                  else work_vector1[j] = 0. ;
                }

            for ( k = 0 ; k < size ; k ++ )
                {
                  phi_node[k] = DDV[v][k]  + theta*phi_w[k]/3. ; // for Bc
                //  phi_node[k] = DDV[v][k] ; // for LLF

		 // Projecting back + dissipation......
                  for ( l = 0 ; l < size ; l ++ )
                       phi_node[k] += ( 1. - theta )*Right[k][l]*work_vector1[l] + theta*K_i_p0[v][k][l]*phi_w[l] ; // for Bc
//                        phi_node[k] += Right[k][l]*work_vector1[l] + theta*K_i_p0[v][k][l]*phi_w[l] ; // for LLf
                }

            residual_update( element[e].node[v] ) ;
         }
}

