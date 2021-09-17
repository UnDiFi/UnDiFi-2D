/***************************************************************************
                                 read_write.c
                                 ------------
            This is read_write: it contains every function reading
                           or writing the solution
                             -------------------
    begin                : Mon May 6 2002
    copyright            : (C) 2002 by Mario Ricchiuto
    email                : ricchiut@vki.ac.be
***************************************************************************/
#include "common.h"

extern int NN, NE, movie, time_step ;
extern int NBF, initial_state ;

extern double gm, time ;

extern struct element_struct   *element ;
extern struct node_struct      *node ;
extern struct boundary_struct *boundary ;

extern char out_file[MAX_CHAR] ;
extern char initial_solution_file[MAX_CHAR] ;

FILE *out, *checkout ;
FILE *initial ;

/*********************************************************/
/*********************************************************/
/**         2D Euler for perfect gas.........           **/
/*********************************************************/
/*********************************************************/

void write_solution()
{
  char tecplot_file_name[MAX_CHAR] ;
  char check_tecplot_file_name[MAX_CHAR] ;

  int e, n, TOT_N, TOT_E, lev, NNN ;
  
  double rho, u, v, p, Ma, s, T, en, h1, h2, xc, yc, omega, mm, r, pm, dp, h0 , H;
	
  lev   = 2 ;
  TOT_N = NN ;
  TOT_E = NE ;
  
  printf( "          *************************************\n" ) ;
  printf( "          **     Writing tecplot file.....   **\n" ) ;
  printf( "          *************************************\n" ) ;
  printf( "\n"                                                ) ;


  if ( movie )	
    {
      if (time_step < 100000)
	sprintf( tecplot_file_name, "./NEO_data/output/%s0%d.dat", out_file, time_step ) ;
      if (time_step < 10000)
	sprintf( tecplot_file_name, "./NEO_data/output/%s00%d.dat", out_file, time_step ) ;
      if (time_step < 1000)
	sprintf( tecplot_file_name, "./NEO_data/output/%s000%d.dat", out_file, time_step ) ;
      if (time_step < 100)
	sprintf( tecplot_file_name, "./NEO_data/output/%s0000%d.dat", out_file, time_step ) ;
      if (time_step < 10)
	sprintf( tecplot_file_name, "./NEO_data/output/%s00000%d.dat", out_file, time_step ) ;    
    }
  else          sprintf( tecplot_file_name, "./NEO_data/output/%s.dat", out_file ) ;
  

  FILE *fp;
  fp=fopen("check.dat", "w");
  fprintf( fp, "TITLE      =  Unstructured grid data \n"                                 ) ;
  fprintf( fp, "VARIABLES  =  x  y  a b c d \n"                              ) ;
  fprintf( fp, "ZONE    N  =  %d    E  =  %d    F = FEPOINT    ET = TRIANGLE \n", TOT_N, TOT_E ) ;
  fprintf( fp, "\n"            ) ;
  for ( n = 0 ; n < TOT_N ; n++ )
    {
     fprintf( fp, "%.20f %.20f %.20f %.20f %.20f %.20f\n",node[n].coordinate[0], node[n].coordinate[1], node[n].P[2][0],node[n].P[2][1],node[n].P[2][2],node[n].P[2][3] ) ;
    }
    for ( e = 0 ; e < NE ; e++ )
    fprintf( fp, "%d %d %d \n", element[e].node[0] + 1,
             element[e].node[1] + 1,
             element[e].node[2] + 1 ) ; 
   fclose(fp) ;

  out = fopen( tecplot_file_name, "w" ) ;

  fprintf( out, "TITLE      =  Unstructured grid data \n"                                 ) ;
  fprintf( out, "VARIABLES  =  x  y  rho  u  v p H  Ma  s T time  sensor\n"                              ) ;
  fprintf( out, "ZONE    N  =  %d    E  =  %d    F = FEPOINT    ET = TRIANGLE \n", TOT_N, TOT_E ) ;
  fprintf( out, "\n"                                                                      ) ;
  
  if  ( initial_state == 9 ) h1 = h2 = 0. ;

  for ( n = 0 ; n < TOT_N ; n++ )
    {
      rho = node[n].P[lev][0] ; rho = MAX( rho, 0.00000000000000000001 ) ;
      u   = node[n].P[lev][1]/rho ;
      v   = node[n].P[lev][2]/rho ;
      en  = node[n].P[lev][3]/rho - 0.50000000000000000000*( u*u + v*v ) ;
      
      p = ( gm - 1.0000000000000000000 )*rho*en ; p = MAX( p, 0.00000000000000000001 ) ;
      T = gm*en ;
      
      T *= gm/( gm - 1.00000000000000000000 ) ;
      s   = p/pow( rho, gm ) ;
      Ma  = sqrt( u*u + v*v )/sqrt( gm*p/rho ) ;
      H = (gm*p)/((gm-1.00000000000000000000)*rho) + (u*u + v*v)*0.50000000000000000000  ;
      
//      fprintf( out, "%le %le %le %le %le %le %le %le %le %le %le %le\n",node[n].coordinate[0], node[n].coordinate[1], rho, u, v, p, H, Ma, s, T, time, node[n].tt ) ;
//      fprintf( out, "%Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf\n",node[n].coordinate[0], node[n].coordinate[1], rho, u, v, p, H, Ma, s, T, time, node[n].tt ) ;
//      fprintf( out, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",fprintf( out, "%Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf\n",node[n].coordinate[0], node[n].coordinate[1], rho, u, v, p, H, Ma, s, T, time, node[n].tt ) ;
        fprintf( out, "%.20f %.20f %.20f %.20f %.20f %.20f %.20f %.20f %.20f %.20f %.20f %.20f\n",node[n].coordinate[0], node[n].coordinate[1], rho, u, v, p, H, Ma, s, T, time, node[n].tt ) ;

//////
//        printf( "%.20f %.20f %.20f %.20f %.20f %.20f\n",  node[n].P[2][0],node[n].P[2][1],node[n].P[2][2],node[n].P[2][3] ) ;
//        printf(checkout,"check me out\n") ;
/////

      if ( initial_state == 9 )
	{
	  xc    = node[n].coordinate[0] - 0.5 ;	 
	  yc    = node[n].coordinate[1] - 0.5 ;	 
	  r     = sqrt( xc*xc + yc*yc ) ;
          
	  pm    = 100. ; 
	  mm = p - pm ;
	  
	  dp = 0. ;
	  
	  if ( r <= 0.25 )
	    {
	      omega = 4.*pi*r ;
	      omega = 15.*( 1. + cos( omega ) ) ;			   
	      r     = 4.*pi*r ;
	      dp    = ( 225.*1.4/( 16.*pi*pi ) )*( 2.*cos( r ) + 2.*r*sin( r ) + cos( 2.*r )/8. + r*sin( 2.*r )/4. + r*r*12./16. ) ;
	      dp   -= ( 225.*1.4/( 16.*pi*pi ) )*( 2.*cos( pi ) + 2.*pi*sin( pi ) + cos( 2.*pi )/8. + pi*sin( 2.*pi )/4. + pi*pi*12./16. ) ;
	    }
	  
	  mm -= dp ;		   			   
	  
	  mm = fabs( mm )/pm ;
	  
	  if ( mm > h0 ) h0 = mm  ;
	  h1 += mm/(1.*NN) ;
	  h2 += mm*mm/(1.*NN) ;
	}
    }	
  if (initial_state == 9 ) {
    h2 = sqrt( h2 ) ;
    
    printf("Errors\n") ;
    printf("Infty: %le\n",log10(h0)) ;
    printf("one: %le\n",log10(h1)) ;
    printf("two: %le\n",log10(h2)) ;
    printf("1/40: %le\n",log(1./40.)) ;
    printf("1/80: %le\n",log(1./80.)) ;
    printf("1/160: %le\n",log(1./160.)) ;
    printf("1/200: %le\n",log(1./200.)) ;
    printf("Ntot:%d\n",TOT_N) ;
  }

  
  for ( e = 0 ; e < NE ; e++ ) 
    fprintf( out, "%d %d %d \n", element[e].node[0] + 1,
	     element[e].node[1] + 1,
	     element[e].node[2] + 1 ) ;
  
  fclose( out ) ;
}

//+++++++++++++++++++++++++++++++++++++++++++++++

void read_solution()
{
  char init_file_name[MAX_CHAR] ;
  
  int n, dummy1, dummy2, NTOT ;
  
  double rho, u, v, p, Ma, s, T, dummy3, dummy4, dummy5, dummy6, H;
  // 11111111111111111111111111111111111111111111111111111111111
  int e, TOT_N, TOT_E, NNN ;
  
  double en, h1, h2, xc, yc, omega, mm, r, pm, dp, h0;
  
  TOT_N = NN;
  // 1111111111111111111111111111111111111111111111111111111111
  
  
  printf( "          *************************************\n" ) ;
  printf( "          **     Reading tecplot file.....   **\n" ) ;
  printf( "          *************************************\n" ) ;
  printf( "\n"                                                ) ;
  
  sprintf( init_file_name, "./NEO_data/output/%s.dat", initial_solution_file ) ;

  initial = fopen( init_file_name, "r" ) ;
  
  if (initial == NULL)
    {
      printf("Pb ouverture fichier" );
    }
  
  NTOT = NN ;
  
  fscanf(initial,"TITLE      =  Unstructured grid data \n");
  fscanf(initial,"VARIABLES  =  x  y  rho  u  v p H  Ma  s T time  sensor\n");
  fscanf(initial,"ZONE    N  =  %d    E  =  %d     F = FEPOINT    ET = TRIANGLE \n",&dummy1, &dummy2);
  fscanf(initial, "\n") ;
  
  for ( n = 0 ; n < NTOT ; n++ )
    {

//    fscanf(initial, "%le %le %le %le %le %le %le %le %le %le %le %le\n", &dummy3, &dummy4, &rho, &u, &v, &p, &H, &Ma, &s, &T, &dummy5, &dummy6);
//    fscanf(initial, "%Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf\n", &dummy3, &dummy4, &rho, &u, &v, &p, &H, &Ma, &s, &T, &dummy5, &dummy6);
      fscanf(initial, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", &dummy3, &dummy4, &rho, &u, &v, &p, &H, &Ma, &s, &T, &dummy5, &dummy6);
    
      node[n].P[2][0] = rho ;
      node[n].P[2][1] = u*rho ;
      node[n].P[2][2] = v*rho ;
      node[n].P[2][3] = p/( gm - 1.0000000000000000000) + 0.50000000000000000000*rho*( u*u + v*v ) ;
    }	
  fclose(initial); 
}


/*
void write_initial_solution()
{
  char tecplot_file_name[MAX_CHAR] ;

  int e, n, TOT_N, TOT_E ;
  
  double rho, u, v, p, Ma, s, T, en ;
  
  TOT_N = NN ;
  TOT_E = NE ;
  
  printf( "          *************************************\n" ) ;
  printf( "          **     Writing tecplot file.....   **\n" ) ;
  printf( "          *************************************\n" ) ;
  printf( "\n"                                                ) ;
  
  sprintf( tecplot_file_name, "./NEO_data/output/initial_solution.dat" ) ;
  out = fopen( tecplot_file_name, "w" ) ;

  fprintf( out, "TITLE      =  Unstructured grid data \n"                                 ) ;
  fprintf( out, "VARIABLES  =  x  y  rho  u  v  p  Ma  s  T \n"                              ) ;
  fprintf( out, "ZONE    N  =  %d    E  =  %d    F = FEPOINT    ET = TRIANGLE \n", TOT_N, TOT_E ) ;
  fprintf( out, "\n"                                                                      ) ;
  
  for ( n = 0 ; n < TOT_N ; n++ )
    {
      rho = node[n].P[2][0] ; rho = MAX( rho, 0.00000000000000000001 ) ;
      u   = node[n].P[2][1]/rho ;
      v   = node[n].P[2][2]/rho ;
      en  = node[n].P[2][3]/rho - 0.5000000000000000000*( u*u + v*v ) ;
      
      p = ( gm - 1.00000000000000000000 )*rho*en ; p= MAX( p , 0.00000000000000000001 ) ;
      T = gm*en ;
      
      T *= gm/( gm - 1.00000000000000000000 ) ;
      s   = p/pow( rho, gm ) ;
      Ma  = sqrt( u*u + v*v )/sqrt( gm*p/rho ) ;
      
//    fprintf( out, "%le %le %le %le %le %le %le %le %le\n", node[n].coordinate[0], node[n].coordinate[1], rho, u, v, p, Ma, s, T ) ;
//    fprintf( out, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", node[n].coordinate[0], node[n].coordinate[1], rho, u, v, p, Ma, s, T ) ;
//    fprintf( out, "%Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf\n" ,node[n].coordinate[0], node[n].coordinate[1], rho, u, v, p, Ma, s, T ) ;
      fprintf( out, "%.20f %.20f %.20f %.20f %.20f %.20f %.20f %.20f %.20f %.20f %.20f %.20f\n",node[n].coordinate[0], node[n].coordinate[1], rho, u, v, p, Ma, s, T ) ;
    }
  
  for ( e = 0 ; e < NE ; e++ ) fprintf( out, "%d %d %d \n", element[e].node[0] + 1,
					element[e].node[1] + 1,
					element[e].node[2] + 1 ) ;
  fclose( out ) ;

}
*/



/***************************************************************************
 *                                                                         *
 *   This program has been developed at the von Karman Institute for Fluid *
 *   Dynamics by Mario Ricchiuto. Any change or update must be done in a   *
 *   clean way and it must be clearly commented.                           *
 *   Mario Ricchiuto                                                       *
 *                                              23/04/2002                 *
 ***************************************************************************/
