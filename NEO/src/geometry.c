/***************************************************************************
                                  geometry.c
                                  ----------
 This is geometry: here the external grid is read and preprocessing invoked
   to store and compute all the usefull informations for the computation
                              -------------------
    Code developed by M. Ricchiuto, Inria - mario.ricchiuto@inria.fr   
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

