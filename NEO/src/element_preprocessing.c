/***************************************************************************
                           element_preprocessing.c
                           -----------------------
This is element_preprocessing: the geometry of each element is computed and
 processed here to have all the necessary informations for the computation
                             -------------------
    begin                : Thu May 2 2002
    copyright            : (C) 2002 by Mario Ricchiuto
    email                : ricchiut@vki.ac.be
***************************************************************************/
#include "common.h"

extern int NN, NE, NBF, NBN, initial_state ;

extern struct element_struct *element ;
extern struct node_struct    *node ;
extern struct boundary_struct *boundary ;
extern struct boundary_nodes_struct   *b_nodes ; 

extern double dt ;

extern char grid_file[MAX_CHAR] ;

/************************************************/
/**  Preprocessing for linear P1 interpolation **/
/************************************************/

void element_preprocessing()
{
  int e, n, n1, n2, f, m, m1, m2, ff ;
  double x1, y1, nx1, ny1, nx2, ny2, nx, ny  ;
  double x2, y2, dummy;
  double div, vel_x, vel_y ;
  int  ver, e1, e2, e3, e4, dum ;
  FILE *out;

  /*******************************************/
  /**       Inizialization of Volumes       **/
  /*******************************************/

  for  ( n = 0 ; n < NN ; n ++ )
    {
      node[n].vol = 0. ;
      node[n].vol_mod = 0. ; 
    }

  /*******************************************/
  /**    Computation of the nodal nomals    **/
  /*******************************************/
     
  for ( e = 0 ; e < NE ; e++ )
    {
      n1 = element[e].node[1] ;
      n2 = element[e].node[2] ;
	
      x1 = node[n1].coordinate[0] ;
      y1 = node[n1].coordinate[1] ;
	
      x2 = node[n2].coordinate[0] ;
      y2 = node[n2].coordinate[1] ;
	
      element[e].normal[0][0] = y1 - y2 ;
      element[e].normal[0][1] = x2 - x1 ;
	
      n1 = element[e].node[2] ;
      n2 = element[e].node[0] ;
	
      x1 = node[n1].coordinate[0] ;
      y1 = node[n1].coordinate[1] ;
	
      x2 = node[n2].coordinate[0] ;
      y2 = node[n2].coordinate[1] ;
	
      element[e].normal[1][0] = y1 - y2 ;
      element[e].normal[1][1] = x2 - x1 ;
	
      n1 = element[e].node[0] ;
      n2 = element[e].node[1] ;
	
      x1 = node[n1].coordinate[0] ;
      y1 = node[n1].coordinate[1] ;
	
      x2 = node[n2].coordinate[0] ;
      y2 = node[n2].coordinate[1] ;
	
      element[e].normal[2][0] = y1 - y2 ;
      element[e].normal[2][1] = x2 - x1 ;
	
      /**************************************************/
      /**  Computation of the volumes of the element   **/
      /**************************************************/

      element[e].volume = 0.5*fabs( element[e].normal[2][0]*element[e].normal[1][1] -
				    element[e].normal[1][0]*element[e].normal[2][1] ) ;

      div = 0. ;
      for ( ver = 0 ; ver < 3 ; ver ++ )
	{ 
	  vel_x = node[element[e].node[ver]].vel[0] ;
	  vel_y = node[element[e].node[ver]].vel[1] ; 

	  div += 0.5* ( vel_x*element[e].normal[ver][0] + vel_y*element[e].normal[ver][1] ) ;								     
	} 

      for ( ver = 0 ; ver < 3 ; ver ++ )
	{  
	  /**   Computation of the volume of the element (Dual Cell) **/   
	  node[element[e].node[ver]].vol += element[e].volume/3. ;
	  /**   Computation of the modified volume of the element (ALE Dual Cell)  **/
	  node[element[e].node[ver]].vol_mod += element[e].volume/3. + dt/6.*div ; 												     
	}   
    }

// ADD -- ADD -- ADD -- ADD -- ADD -- ADD -- ADD -- ADD -- 1  
  for (n = 0; n<NN; n++)
    {
      if (node[n].vol == 0.00000000000000)
	{
	  node[n].vol     = 1.;
          node[n].vol_mod = 1.;
	}
     }
// ------------------------------------------------------- 1

  /**********************************************/
  /*******      Boundary Normals          *******/
  /**********************************************/

  m = 0 ;

  for ( f = 0 ; f < NBF ; f ++ )  
    {
      n1 = boundary[f].node[0] ;
      n2 = boundary[f].node[1] ;

      x1 = node[n1].coordinate[0] ;
      y1 = node[n1].coordinate[1] ;
	
      x2 = node[n2].coordinate[0] ;
      y2 = node[n2].coordinate[1] ;
	
      boundary[f].normal[0] = y1 - y2 ;
      boundary[f].normal[1] = x2 - x1 ;
      
//      if ( (fabs(x1 - 10.0) <= 1.e-7) && (fabs(x2 - 10.0) <= 1.e-7) )
//           boundary[f].type = 4 ;   
    
    }


  for ( n = 0 ; n < NN ; n++ ) node[n].flag = 0 ;


}
/***************************************************************************
 *                                                                         *
 *   This program has been developed at the von Karman Institute for Fluid *
 *   Dynamics by Mario Ricchiuto. Any change or update must be done in a   *
 *   clean way and it must be clearly commented.                           *
 *   Mario Ricchiuto                                                       *
 *                                              23/04/2002                 *
 ***************************************************************************/
