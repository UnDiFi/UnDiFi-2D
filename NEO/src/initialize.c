/***************************************************************************
                                 initialize.c
                                 ------------
This is initialize: it initializes some model, geometry, scheme and case
                      dependent constant and parameters
                             -------------------
    begin                : Wed Apr 24 2002
    copyright            : (C) 2002 by Mario Ricchiuto
    email                : ricchiut@vki.ac.be
 ***************************************************************************/
#include "common.h"

extern int model, size, time_step_max, scheme, problem ;
extern int vertex ;
extern int NBN, NN, NE, face_q_pts, volume_q_pts ;

extern double p_out, Mach_inf, alpha_inf, rho_inf , p_inf ; 
extern double Mach_in, alpha_in, rho_in , gm ; 
extern double c_tau, zero_r ;

extern void read_spatial_grid_header() ;

void initialize()
{
/**************************************************************************/
/**                SIZE of the system of equations                       **/
/**************************************************************************/

     size = 4 ;

/**************************************************************************/
/**         Number of VERTICES of each (spatial) sub-element             **/
/**************************************************************************/

      vertex = 3 ;

/**************************************************************************/
/**       READING the GRID_FILE HEADER to determine the number of        **/
/**   nodes (NN),  elements (NE), faces (NF) and boundary faces (NBF)    **/
/**************************************************************************/

      read_spatial_grid_header() ;

/*************************************************************************/
/**            Interpolation-Order Dependent Constants                  **/
/*************************************************************************/

      NBN = 2 ;

/*    face_q_pts = 4 ;*/
      face_q_pts = 3 ;
      volume_q_pts = 4 ;
      c_tau = 0.5 ;

/*************************************************************************/
/**            Boundary Conditions - outlet, far field                  **/
/*************************************************************************/

        switch ( problem ) 
        {                                                               
              case 6: /*NACA0012_M080_A0*/
                rho_inf   = 1.0  ;                                      
                alpha_inf = atan2(0.017452406,0.999847695)   ;
                Mach_inf  = 0.8 ;                                       
                p_inf     = rho_inf / (gm * Mach_inf*Mach_inf) ;        
              break ;
              case 7: /*NACA0012_M095_A0_FISHTAIL-1*/                   
                rho_inf   = 1.0  ;
                alpha_inf = 0.0  ;
                Mach_inf  = 0.95  ;
                p_inf     = rho_inf / (gm * Mach_inf*Mach_inf) ;
              break ;
              case 8: /*Q1D*/
                p_out = 0.66085998521698897 ;
              break ;
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
