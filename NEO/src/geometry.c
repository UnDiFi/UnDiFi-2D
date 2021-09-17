/***************************************************************************
                                  geometry.c
                                  ----------
 This is geometry: here the external grid is read and preprocessed to store
       and compute all the usefull informations for the computation
                             -------------------
    begin                : Thu May 2 2002
    copyright            : (C) 2002 by Mario Ricchiuto
    email                : ricchiut@vki.ac.be
 ***************************************************************************/
#include "common.h"

extern void read_spatial_grid() ;
extern void element_preprocessing() ;
extern void set_numerical_integration() ;

 void geometry()
{
/***********************************/
/**   reading the external grid   **/
/***********************************/

		printf( "          *************************************\n" ) ;
    printf( "          **     Reading the gridfile.....   **\n" ) ;
    printf( "          *************************************\n" ) ;
		printf( "\n"                                       ) ;

		read_spatial_grid() ;

/***********************************/
/**    geometry pre-processing    **/
/***********************************/		
			
	  printf( "          *************************************\n" ) ;
    printf( "          **  Processing element geometry... **\n" ) ;
    printf( "          *************************************\n" ) ;
	  printf( "\n"                                       ) ;

    element_preprocessing() ;

    printf( "          *************************************\n" ) ;
    printf( "          **    Numerical integration...     **\n" ) ;
    printf( "          *************************************\n" ) ;
          printf( "\n"                                       ) ;

          set_numerical_integration() ;
}

/***************************************************************************
 *                                                                         *
 *   This program has been developed at the von Karman Institute for Fluid *
 *   Dynamics by Mario Ricchiuto. Any change or update must be done in a   *
 *   clean way and it must be clearly commented.                           *
 *   Mario Ricchiuto                                                       *
 *                                              23/04/2002                 *
 ***************************************************************************/
