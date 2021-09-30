/***************************************************************************
                              nodal_residuals.c
                              -----------------
This is nodal_residuals: it performs the computation of the nodal resusuals
                       within the main element loop 
                             -------------------
    Code developed by M. Ricchiuto, Inria - mario.ricchiuto@inria.fr
 ***************************************************************************/
#include "common.h"

extern int NE, NN ;

extern struct element_struct *element ;
extern struct node_struct *node ;

extern void STfluctuation2( int )  ;
extern void ( *distribute_fluctuation )( int ) ;
extern void initvec( double * ) ;
extern void consistent_vars( int ) ;
extern void compute_upwind_Keys( int ) ;
extern void local_average( int ) ;

void nodal_residuals()
{
     int e, n ;

/********************************************/
/** initialization of flags and residuals  **/
/********************************************/

     for  ( n = 0 ; n < NN ; n ++ )
          {
            node[n].flag = 0 ;

            initvec( node[n].Res ) ;
 
            node[n].tt = 0. ;
          }

/********************************************/
/**       main loop over the elements      **/
/********************************************/

     for ( e = 0 ; e < NE ; e ++ )
         {
          consistent_vars( e ) ;

          local_average( e ) ; 

          compute_upwind_Keys( e ) ;

          STfluctuation2( e ) ;

          distribute_fluctuation( e ) ;
        }
}


