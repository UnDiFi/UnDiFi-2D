/***************************************************************************
                             boundary_residual.c
                             -------------------
  This is boundary_residual: it drives the computation of the additional
            residual coming from the boundary conditions.
   Colored boundaries correspond to
    case 0: Fitted shocks - the soution is set already form the RH conditions
    case 1: Supersonic inlet - the solution is fully (strongly) imposted
    case 2: Supersonic outlet - nothing happens (as in case 0)
    case 3: Wall - the normal component of the momentum residual is set to zero 
    case 4: Sub-sonic outlet using weak conditions based on a Steger-Warming flux
    case 5: Far field/inled - Full state weakly imposed depending on the charateristics sign
                               (weak conditions based on a Steger-Warming flux)
                             -------------------
    Code developed by M. Ricchiuto, Inria - mario.ricchiuto@inria.fr
 ***************************************************************************/
#include "common.h"

extern struct node_struct *node ;
extern struct boundary_struct  *boundary ;
extern struct boundary_nodes_struct *b_nodes ;
extern char grid_file[MAX_CHAR] ;

extern double *temp_normal ;

extern int NN, NBN, initial_state, size, NBF ;

extern void periodicity( int ) ;/* Face global index */
extern void ( *supersonic_residual )( int ) ;/* Face global index */
extern void first_streamline( int ) ;
extern void outlet_2D_Euler( int ) ;
extern void inlet_2D_Euler( int ) ;
extern void far_2D_Euler( int ) ;
extern void out_2D_Euler( int ) ;

void boundary_conditions()
{
     int f, n ;
      int e1, e2, e3, e4;

//    printf("\n") ;printf("BC.s\n") ;printf("\n") ;

     for ( n = 0 ; n < NN ; n ++ ) node[n].flag = 0 ;

     for ( f = 0 ; f < NBF ; f ++ )
         {
           switch ( boundary[f].type )
                  {
                    case 10:/*Periodic BCs*/
                            periodicity( f ) ;
                            break ;
                    case 0:/*Fitted shocks*/
		            /* nothing happens */
                            break ;
                    case 1:/*Supersonic inlet*/
                            supersonic_residual( f ) ; /*strong BC*/
//                            inlet_2D_Euler( f ) ; /*weak BC*/
                            break ;
                    case 2:/*Supersonic outlet*/
		            /* nothing happens */
                            break ;
                    case 3:/*Weak Wall*/
                            first_streamline( f ) ;
                            break ;
                    case 4:/*outlet*/
                            out_2D_Euler( f ) ;
                            break ;
                    case 5:/*far field*/
                            far_2D_Euler( f ) ;
                            break ;
                  }
         }
}
