/***************************************************************************
                                 wrap_it_up.c
                                 ------------
  This is wrap_it_up: it ends the computation and gives some...final infos
                             -------------------
    begin                : Wed May 15 2002
    copyright            : (C) 2002 by Mario Ricchiuto
    email                : ricchiut@vki.ac.be
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

/***************************************************************************
 *                                                                         *
 *   This program has been developed at the von Karman Institute for Fluid *
 *   Dynamics by Mario Ricchiuto. Any change or update must be done in a   *
 *   clean way and it must be clearly commented.                           *
 *   Mario Ricchiuto                                                       *
 *                                              23/04/2002                 *
 ***************************************************************************/
