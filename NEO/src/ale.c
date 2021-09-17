/***************************************************************************
                                     ale.c
                           -----------------------
				   ale stuff
                             -------------------
    begin                : Thu May 2 2002
    copyright            : (C) 2002 by Mario Ricchiuto
    email                : ricchiut@vki.ac.be
***************************************************************************/

#include "common.h"

extern struct node_struct *node ;

extern int NN ;
extern double time, dt, t_max ;

extern void ( *move_the_mesh )( int n ) ;

//**************************************************************************

void move_compute_save_the_mesh()
{
  int n ;
  FILE *vel_file;

  /*********************************************/
  /***      Reading of the velocity file     ***/
  /*********************************************/
  vel_file = fopen("./NEO_data/input/vel.dat","r");
  fscanf(vel_file,"        u                         v\n");
  for (n = 0; n < NN; n++)
    {
      fscanf(vel_file,"%le %le\n",&node[n].vel[0],&node[n].vel[1]);
    }
  fclose(vel_file);
//printf("About to move the mesh");

  for ( n = 0 ; n < NN ; n++ )
    {
      
      //printf(" %le   %le \n",node[n].vel[0],node[n].vel[1]);
      
      /**********************************************/
      /*******        Moving the mesh         *******/
      /**********************************************/
      
      //move_the_mesh( n ) ;
      
      /**********************************************/
      /****  Computing Velocity according Farhat  ***/
      /**********************************************/
      
      //node[n].vel[0] = ( node[n].coordinate_new[0] - node[n].coordinate_old[0] )/dt ;
      //node[n].vel[1] = ( node[n].coordinate_new[1] - node[n].coordinate_old[1] )/dt ;

      /**********************************************/
      /****  Computing new_node according Farhat  ***/
      /**********************************************/
      node[n].coordinate_new[0] = node[n].vel[0]*dt + node[n].coordinate_old[0] ;
      node[n].coordinate_new[1] = node[n].vel[1]*dt + node[n].coordinate_old[1] ;

      /**********************************************/
      /****   Computing Midpoint configuration    ***/
      /**********************************************/

      node[n].coordinate[0] = 0.5*( node[n].coordinate_new[0] + node[n].coordinate_old[0] ) ;
      node[n].coordinate[1] = 0.5*( node[n].coordinate_new[1] + node[n].coordinate_old[1] ) ;

      /**********************************************/
      /*******         Saving the mesh        *******/
      /**********************************************/

      node[n].coordinate_old[0] = node[n].coordinate_new[0] ;
      node[n].coordinate_old[1] = node[n].coordinate_new[1] ;
      
    }
}

void sinusoidal( int n )  
{
  double csi;
  double eta;    

  csi = node[n].coordinate_0[0] ;
  eta = node[n].coordinate_0[1] ;        

  node[n].coordinate_new[0] = csi + 0.1*sin(2.*pi*csi)*sin(2.*pi*eta)*sin(2.*pi*time/t_max) ; 
  node[n].coordinate_new[1] = eta + 0.1*sin(2.*pi*csi)*sin(2.*pi*eta)*sin(2.*pi*time/t_max) ;
   
}

void flap_deflection( int n )  
{
  double csi, eta;
  double u0, tau ;
  /*-----Deflection: step data-----*/
  double switch_time1 = 1.25 ;
  double switch_time2 = 2.50 ;
  double switch_time3 = 3.75 ;
  double alpha_max = 20. ; 		// [deg°]
  double transient_time = 0.25 ;    

  csi = node[n].coordinate_0[0] ;
  eta = node[n].coordinate_0[1] ; 

  tau = transient_time/5. ;

  if ( csi > 0.25 )
    {
      if ( time <= switch_time1 ) 
	u0 = tan ( (alpha_max*pi/180.) * (1.-exp(-time/tau) ) ) *(csi-.25) ;
      if ( time > switch_time1 && time <= switch_time2 )
	u0 = tan ( alpha_max*pi/180. - 2*(alpha_max*pi/180.)*(1.-exp(-(time-switch_time1)/tau) ) ) *(csi-.25) ;
      if ( time > switch_time2 && time <= switch_time3)
	u0 = tan ( -alpha_max*pi/180. + 2*(alpha_max*pi/180.)*(1.-exp(-(time-switch_time2)/tau) ) ) *(csi-.25) ;
      if ( time > switch_time3 )
	u0 = tan ( alpha_max*pi/180. - 2*(alpha_max*pi/180.)*(1.-exp(-(time-switch_time3)/tau) ) ) *(csi-.25) ;
       
      node[n].coordinate_new[0] = csi ;
      node[n].coordinate_new[1] = eta + (1.-eta)*u0 ;
    }
  else
    {
      node[n].coordinate_new[0] = csi ;
      node[n].coordinate_new[1] = eta ;
    }
}

void centered_compression( int n )
{
  double csi, eta;
  double x0, xp, L, vp;
  double L0 = 5.0;
  double alpha;

  vp = 0.86977177750456247;

  csi = node[n].coordinate_0[0];
  eta = node[n].coordinate_0[1]; 

  // Position of the piston
  x0 = node[0].coordinate_0[0];

  // Distance between the piston and the wall
  L = L0 - x0;

  // Computation of alpha
  alpha = csi/L;

  // Piston's law
  //xp = L*vp*(1-exp(-time));
  xp = vp*(1-exp(-time));

  // Compute the position of new node
  node[n].coordinate_new[0] = (1-alpha)*xp + alpha*L  ;
  node[n].coordinate_new[1] = eta;
}
