/***************************************************************************
                                  driver.c
                                   -----
   This is the driver: it starts the computation and drives it to the end.
                             -------------------
    begin                : Tue Apr 23 2002
    copyright            : (C) 2002 by Mario Ricchiuto
    email                : ricchiut@vki.ac.be
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

/***************************************************************************
 *                                                                         *
 *   This program has been developed at the von Karman Institute for Fluid *
 *   Dynamics by Mario Ricchiuto. Any change or update must be done in a   *
 *   clean way and it must be clearly commented.                           *
 *   Mario Ricchiuto                                                       *
 *                                              23/04/2002                 *
 ***************************************************************************/
