/***************************************************************************
                                  common.h
                                  -------
This is the common: it defines the data structures used in the computations
                             and some constants.
                              -------------------
    Code developed by M. Ricchiuto, Inria - mario.ricchiuto@inria.fr
 ***************************************************************************/
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <malloc.h>
#include "f2c.h"
#include "clapack.h"

#define MAX(a,b) ((a)>(b) ? (a) : (b))
#define MIN(a,b) ((a)<(b) ? (a) : (b))

#define pi  3.14159265358979
#define MAX_CHAR 150
#define infinity ( double )( 1.e+50 )
#define epsilon ( double )( 1.e-12 )
#define small ( double )( 1.e-10 )
#define SIGN(a)  ((a)<=0 ? -1 : 1)
#define MINMOD(a) ( 0.5*( 1.0 + SIGN((a)) )*MIN(1,(a)) )



struct element_struct
{
       int *node ;
       double  **normal ;
       double volume ;
} ;

struct node_struct
{
       int mate ;
       int flag ;
       double **P ;
       double **Z ;
       double *Res ;
       double *coordinate_new ;
       double *coordinate_old ;
       double *coordinate ;
       double *coordinate_0 ;
       double dtau ;
       double vol ;
       double vol_mod ;
       double tt ;
       double *vel ;
} ;

struct boundary_struct
{
       int *node ;
       int *types ;
       int type ;
       double *normal ;
} ;


struct numerical_int_struct
{
       double **face_coordinate ;
       double *face_weight ;
       double **volume_coordinate ;
       double *volume_weight ;
} ;


struct boundary_nodes_struct
{
      int node ;
      int type ;
      int mate ;
      int flag ;
      double *normal ;
};


