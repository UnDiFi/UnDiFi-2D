/***************************************************************************
                         set_numerical_integration.c
                         ---------------------------
This is set_numerical_integration: the weights for the numerical evaluation
           of line, surface and volume integrals are set here
                             -------------------
    begin                : Thu May 2 2002
    copyright            : (C) 2002 by Mario Ricchiuto
    email                : ricchiut@vki.ac.be
 ***************************************************************************/
#include "common.h"

extern struct numerical_int_struct numerical_int ;

extern int iteration_max ;

extern double a1[3], a2[3], a3[3] ;
extern double b1[3], b2[3], b3[3] ;

void set_numerical_integration()
{
     double aa0, aa1, s ;

     /* high order numerical integration */ 

     s = sqrt( 0.6 ) ;

     aa0 = ( 1.0 - s )*0.5 ;
     aa1 = ( 1.0 + s )*0.5 ;

     numerical_int.face_coordinate[0][0] = aa0 ;
     numerical_int.face_coordinate[0][1] = aa1 ;

     numerical_int.face_weight[0] = 5.0/18.0 ;


     s = -sqrt( 0.6 ) ;

     aa0 = ( 1.0 - s )*0.5 ;
     aa1 = ( 1.0 + s )*0.5 ;


     numerical_int.face_coordinate[1][0] = aa0 ;
     numerical_int.face_coordinate[1][1] = aa1 ;

     numerical_int.face_weight[1] = 5.0/18.0 ;


     numerical_int.face_coordinate[2][0] = 0.5 ;
     numerical_int.face_coordinate[2][1] = 0.5 ;

     numerical_int.face_weight[2] = 8.0/18.0 ;


/* low order numerical integration */
/*
     aa0 = 1.0 ;
     aa1 = 0.0 ;

     numerical_int.face_coordinate[0][0] = aa0 ;
     numerical_int.face_coordinate[0][1] = aa1 ;

     numerical_int.face_weight[0] = 0.5 ;

     aa0 = 0.0 ;
     aa1 = 1.0 ;


     numerical_int.face_coordinate[1][0] = aa0 ;
     numerical_int.face_coordinate[1][1] = aa1 ;

     numerical_int.face_weight[1] = 0.5 ;


     numerical_int.face_coordinate[2][0] = 0.5 ;
     numerical_int.face_coordinate[2][1] = 0.5 ;

     numerical_int.face_weight[2] = 0.0 ;
*/


     /* Volume integration */
     numerical_int.volume_coordinate[0][0] = 1.0/3.0 ;
     numerical_int.volume_coordinate[0][1] = 1.0/3.0 ;
     numerical_int.volume_coordinate[0][2] = 1.0/3.0 ;

     numerical_int.volume_weight[0] = -27.0/48.0 ;


     numerical_int.volume_coordinate[1][0] = 0.6 ;
     numerical_int.volume_coordinate[1][1] = 0.2 ;
     numerical_int.volume_coordinate[1][2] = 0.2 ;

     numerical_int.volume_weight[1] = 25.0/48.0 ;


     numerical_int.volume_coordinate[2][0] = 0.2 ;
     numerical_int.volume_coordinate[2][1] = 0.2 ;
     numerical_int.volume_coordinate[2][2] = 0.6 ;

     numerical_int.volume_weight[2] = 25.0/48.0 ;


     numerical_int.volume_coordinate[3][0] = 0.2 ;
     numerical_int.volume_coordinate[3][1] = 0.6 ;
     numerical_int.volume_coordinate[3][2] = 0.2 ;

     numerical_int.volume_weight[3] = 25.0/48.0 ;

//
// Initializing the Runge-Kutta weights
//

     switch ( iteration_max )
               {
                 case 1:
                        a1[0] = 0. ;
                        a1[1] = 0. ;
                        a1[2] = 0. ;

                        a2[0] = 0. ;
                        a2[1] = 0. ;
                        a2[2] = 0. ;

                        a3[0] = 1. ;
                        a3[1] = 0. ;
                        a3[2] = 0. ;

                        b1[0] = 0. ;
                        b1[1] = 0. ;
                        b1[2] = 0. ;

                        b2[0] = 0. ;
                        b2[1] = 0. ;
                        b2[2] = 0. ;

                        b3[0] = 0. ;
                        b3[1] = 0. ;
                        b3[2] = 0. ;
                        break ;
                  case 2:
                        a1[0] = 0. ;
                        a1[1] = 0.5 ;
                        a1[2] = 0. ;

                        a2[0] = 0. ;
                        a2[1] = 0. ;
                        a2[2] = 0. ;

                        a3[0] = 1. ;
                        a3[1] = 0.5 ;
                        a3[2] = 0. ;

                        b1[0] = 0. ;
                        b1[1] = -1. ;
                        b1[2] = 0. ;

                        b2[0] = 0. ;
                        b2[1] = 0. ;
                        b2[2] = 0. ;

                        b3[0] = 0. ;
                        b3[1] = 1. ;
                        b3[2] = 0. ;
                        break ;
                   case 3:
                        a1[0] = 0. ;
                        a1[1] = 0.25 ;
                        a1[2] = 1./6. ;

                        a2[0] = 0. ;
                        a2[1] = 0. ;
                        a2[2] = 1./6. ;

                        a3[0] = 1. ;
                        a3[1] = 0.25 ;
                        a3[2] = 2./3. ;

                        b1[0] = 0.  ;
                        b1[1] = -0.5 ;
                        b1[2] = -2. ;

                        b2[0] = 0. ;
                        b2[1] = 0. ;
                        b2[2] = 0. ;

                        b3[0] = 0. ;
                        b3[1] = 0.5 ;
                        b3[2] = 2. ;
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
