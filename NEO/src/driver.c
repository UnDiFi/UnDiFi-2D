/***************************************************************************
                                  driver.c
                                   -----
   This is the main: it starts the computation and drives it to the end.
                             -------------------
    Code developed by M. Ricchiuto, Inria - mario.ricchiuto@inria.fr
 ***************************************************************************/
#include "common.h"
#include "common_variables.h"

extern void get_ready() ;
extern void run() ;
extern void wrap_it_up() ;

int main()
{
     get_ready() ;

     run() ;

     wrap_it_up() ;

     return 0 ;
}

