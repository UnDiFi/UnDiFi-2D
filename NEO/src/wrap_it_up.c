/***************************************************************************
                                 wrap_it_up.c
                                 ------------
                  This is wrap_it_up: the run ends here 
                             -------------------
    Code developed by M. Ricchiuto, Inria - mario.ricchiuto@inria.fr
 ***************************************************************************/
#include "common.h"

extern void write_solution() ;

void wrap_it_up()
{
/**********************************/
/**   write the final solution   **/
/**********************************/

     write_solution() ;

     printf( "\n" ) ;
     printf( "          *************************************\n" ) ;
     printf( "          *************************************\n" ) ;
     printf( "          **                                 **\n" ) ;
     printf( "          **      NEO has finished.....      **\n" ) ;
     printf( "          **                                 **\n" ) ;
     printf( "          **       I hope you'll like        **\n" ) ;
     printf( "          **          the results....        **\n" ) ;
     printf( "          **                                 **\n" ) ;
     printf( "          *************************************\n" ) ;
     printf( "          **                                 **\n" ) ;
     printf( "          **              BYE.....           **\n" ) ;
     printf( "          **                                 **\n" ) ;
     printf( "          *************************************\n" ) ;
     printf( "          *************************************\n" ) ;
     printf( "\n" ) ;
}

