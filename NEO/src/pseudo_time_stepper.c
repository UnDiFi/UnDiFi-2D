/***************************************************************************
                            pseudo_time_stepper.c
                            ---------------------
    This is pseudo_time_stepper : it drives an explicit pseudo time loop
     to convergence through the computation of the nodal residuals and
                              the nodal updates
                             -------------------
    begin                : Tue May 7 2002
    copyright            : (C) 2002 by Mario Ricchiuto
    email                : ricchiut@vki.ac.be
 ***************************************************************************/
#include "common.h"

extern struct node_struct *node ;

extern int iteration, iteration_max, save_solution, steady_check, counter_info ;
extern int time_step, time_step_max, info, NN, size, iteration_max  ;

extern double residual_treshold, residual_norm, steady_norm, ref_norm ;
extern double zero_r, *temp_vector, lim_norm, CFL , dt ;

extern void compute_norm() ;
extern void nodal_residuals() ;
extern void boundary_conditions() ;
extern void write_solution() ;
extern void P_2_Z( double *, double * ) ;
extern void Z_2_P( double *, double * ) ;


void pseudo_time_stepper()
{
     int n, i ;
     double k_dt ;

     iteration = 0 ;

/******************************************/
/**    starting the pseudo time loop     **/
/******************************************/

     do
       {
         iteration += 1 ;

         nodal_residuals() ;

         boundary_conditions() ;


/************************************/
/********    Nodal Update   *********/
/************************************/
         
         for ( n = 0 ; n < NN ; n ++ )
             {  
               for ( i = 0 ; i < size ; i++ )

                     node[n].P[2][i] = node[n].P[0][i] - node[n].Res[i]/node[n].vol_mod ;
                      
                     P_2_Z( node[n].P[2], node[n].Z[2] ) ;
                     Z_2_P( node[n].Z[2], node[n].P[2] ) ;

             }

         if ( iteration == 1 )
            {
              for ( n = 0 ; n < NN ; n ++ )
                    for ( i = 0 ; i < size ; i++ )
                        {  
                          node[n].P[1][i] = node[n].P[2][i];
                          node[n].Z[1][i] = node[n].Z[2][i];
                        }
            }

         compute_norm() ;

         /* write convergence file */
         FILE *fp = fopen("residual_norm.dat", "a");
         fprintf(fp, "%f\n", residual_norm );


        } while (  iteration < iteration_max  ) ;


    if ( counter_info == info )
       {
         printf( "\n" ) ;
         printf( "          *************************************\n" ) ;
         printf( "          *************************************\n" ) ;
         printf( "          **      timestep information:      **\n" ) ;
         printf( "          ** MAX number of iteration reached **\n" ) ;
         printf( "          *************************************\n" ) ;
         printf( "               LOG10( norm1_res ) = %le           \n", residual_norm ) ;
         printf( "               Initial norm       = %le           \n", ref_norm ) ;
         printf( "               Time step          = %le           \n", dt ) ;
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
