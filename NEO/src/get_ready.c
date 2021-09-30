/***************************************************************************
                                 get_ready.c
                                 ----------
        This is get_ready: anything that needs to be done to start up
                  the computation is done through get_ready.
                            -------------------
    Code developed by M. Ricchiuto, Inria - mario.ricchiuto@inria.fr 
 ***************************************************************************/
#include "common.h"

extern void input_parameters() ;
extern void initialize() ;
extern void select_function() ;
extern void memory_allocation() ;
extern void geometry() ;
extern void initial_solution() ;
extern void write_solution() ;

void get_ready()
{
         printf( "\n" ) ;
         printf( "          *************************************\n" ) ;
         printf( "          *************************************\n" ) ;
         printf( "          **                                 **\n" ) ;
         printf( "          **      NEO has started .....      **\n" ) ;
         printf( "          **                                 **\n" ) ;
         printf( "          **       Sit down and relax        **\n" ) ;
         printf( "          **                                 **\n" ) ;
         printf( "          *************************************\n" ) ;
         printf( "          **                                 **\n" ) ;
         printf( "          **    Getting ready to compute...  **\n" ) ;
         printf( "          **                                 **\n" ) ;
         printf( "          *************************************\n" ) ;
         printf( "          *************************************\n" ) ;
         printf( "\n" ) ;

         input_parameters() ;
         initialize() ;
         select_function() ;
         memory_allocation() ;
         geometry() ;
         initial_solution() ;
         write_solution() ;
}

