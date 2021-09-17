/***************************************************************************
                              select_function.c
                              -----------------
          This is select_function: function pointers are assigned here
                             -------------------
    begin                : Thu Apr 25 2002
    copyright            : (C) 2002 by Mario Ricchiuto
    email                : ricchiut@vki.ac.be
***************************************************************************/
#include "common.h"

extern int problem ;
extern int scheme, iteration_max ;

extern void ( *supersonic_residual )( int ) ; /* node global index  */

extern void supersonic_residual_std( int ) ;
extern void supersonic_residual_2DRP_Euler( int ) ;
extern void supersonic_residual_2DRPiso_Euler( int ) ;
extern void supersonic_residual_dMach_Euler( int ) ;


/**************************************************************************/
/**                   SCHEME dependent functions                         **/
/**************************************************************************/

extern void ( *distribute_fluctuation )( int ) ;

extern void LDA_unsteady( int ) ;
extern void LLF_unsteady( int ) ;
extern void LDA_N_unsteady( int ) ;
extern void LW_LF_unsteady( int ) ;
extern void LW_unsteady( int ) ;
extern void N_unsteady( int ) ;


extern void ( *fluctuation2 )( int ) ;

extern void fluctuationRK1( int ) ;
extern void fluctuationRK2( int ) ;
extern void fluctuationRK3( int ) ;

extern void ( *move_the_mesh )( int ) ;

extern void sinusoidal( int ) ;
extern void flap_deflection( int ) ;
extern void centered_compression( int) ;


/**************************************************************************/
/**************************************************************************/
/**                  here begins SELECT_FUNCTION                         **/
/**************************************************************************/
/**************************************************************************/

void select_function()
{
  /**************************************************************************/
  /**                  PROBLEM DEPENDENT SETTINGS !!!                      **/
  /**************************************************************************/

  switch ( problem )
    {
    case 0:
      supersonic_residual = supersonic_residual_std ;
      break ;
    case 1:/* 2D Riemann Problem on a regular mesh- Euler*/
      supersonic_residual = supersonic_residual_2DRP_Euler ;
      break ;
    case 2:/* Double Mach Reflection - Euler */
      supersonic_residual = supersonic_residual_dMach_Euler ;
      break ;
    case 3:/* 2D Riemann Problem on a isotropic mesh (half rotated domain) - Euler*/
      supersonic_residual = supersonic_residual_2DRPiso_Euler  ;
      break ;
    case 4:/* Channel with Deflection of a Flap  */
      supersonic_residual = supersonic_residual_std  ;
      break ;
    case 6:/* NACA0012_M080_A0  */
      supersonic_residual = supersonic_residual_std  ;
      break ;
    case 7:/* NACA0012_M095_A0_FISHTAIL-1 */
      supersonic_residual = supersonic_residual_std  ;
      break ;
    case 8:/* Q1D */
      supersonic_residual = supersonic_residual_std  ;
      break ;
    }

  switch ( problem )
    {
    case 0:
      /* Vortex Advection:  sinusoidal law for a square mesh VORT_T */	
      move_the_mesh = sinusoidal ;
      break ;
    case 1:/* Riemann Problem: sinusoidal law for a square mesh VORT_T */
      move_the_mesh = sinusoidal ;
      break ;
    case 2:/* Double Mach Reflection - Euler */
      printf( "          ALE NOT READY\n" ) ;
      break ;
    case 3:/* 2D Riemann Problem on a isotropic mesh (half rotated domain) - Euler*/
      printf( "          ALE NOT READY\n" ) ;
      break ;
    case 4:/* Channel with Deflection of a Flap: following the boundaries for mesh LINEAR20 */
      move_the_mesh = flap_deflection ;
      break ;
    case 5:/* Centered compression*/
      move_the_mesh = centered_compression;
      break;

    case 6:/* NACA0012_M080_A0  */
      printf( "          ALE NOT READY\n" ) ;
      break ;
    case 7:/* NACA0012_M095_A0_FISHTAIL-1 */
      printf( "          ALE NOT READY\n" ) ;
      break ;
    case 8:/* Q1D */
      printf( "          ALE NOT READY\n" ) ;
      break ;
		      

    }

  /**************************************************************************/
  /**                  MODEL DEPENDENT SETTINGS !!!                        **/
  /**************************************************************************/

  /**************************************************************************/
  /**              STEADY/UNSTEADY dependent settings                      **/
  /**************************************************************************/

  switch ( scheme )
    {
    case 0:
      distribute_fluctuation = LDA_unsteady ;
      break ;
    case 1:
      distribute_fluctuation = LLF_unsteady ;
      break ;
    case 2:
      distribute_fluctuation = LDA_N_unsteady ;
      break ;	     
    case 3:
      distribute_fluctuation = LW_LF_unsteady ;
      break ;	     
    case 4:
      distribute_fluctuation = LW_unsteady ;
      break ;	     
    case 5:
      distribute_fluctuation = N_unsteady ;
      break ;
    }

  switch ( iteration_max )
    {
    case 1 :
      fluctuation2 = fluctuationRK1 ;
      break ; 
    case 2 :
      fluctuation2 = fluctuationRK2 ;
      break ; 
    case 3 :
      fluctuation2 = fluctuationRK3 ;
      break ; 
    }

  /**************************************************************************/
  /**************************************************************************/
  /**                    here ends SELECT_FUNCTION                         **/
  /**************************************************************************/
  /**************************************************************************/
}

/***************************************************************************
 *                                                                         *
 *   This program has been developed at the von Karman Institute for Fluid *
 *   Dynamics by Mario Ricchiuto. Any change or update must be done in a   *
 *   clean way and it must be clearly commented.                           *
 *   Mario Ricchiuto                                                       *
 *                                              23/04/2002                 *
 ***************************************************************************/

