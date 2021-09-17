/***************************************************************************
                              nodal_residuals.c
                              -----------------
This is nodal_residuals: it performs the computation of the nodal resusuals
    in pseudo time by looping over all the elements of the spatial mesh
                             -------------------
    begin                : Fri May 10 2002
    copyright            : (C) 2002 by Mario Ricchiuto
    email                : ricchiut@vki.ac.be
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


/***************************************************************************
 *                                                                         *
 *   This program has been developed at the von Karman Institute for Fluid *
 *   Dynamics by Mario Ricchiuto. Any change or update must be done in a   *
 *   clean way and it must be clearly commented.                           *
 *   Mario Ricchiuto                                                       *
 *                                              23/04/2002                 *
 ***************************************************************************/
