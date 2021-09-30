/***************************************************************************
                             read_spatial_grid.c
                             -------------------
          This is read_spatial_grid: every function reading from the
                   (externally provided) grid_file is here.
MIND: The .grd grid-file format is an in-house format structured as follows:
                               - header      : dimension NE NN NF NBF
                               - blank line
                               - e:=1 -> NE  : e_node_1 ... e_node_vertex
                                               e_neigh_1 ... e_neigh_vertex
                                               e_neigh_face_1 ... e_neigh_face_vertex
                               - blank line
                               - n:=1 -> NN  : n_x_1 ... n_x_dimension n_mate
                               - blank line
                               - f:=1 -> NBF : f_n1... f_ndimension
                                               type_n1... type_ndimension
                                               f_type
                          -------------------
    Code developed by M. Ricchiuto, Inria - mario.ricchiuto@inria.fr
 ***************************************************************************/
#include "common.h"

extern int NE, NN, NBF, NBN, initial_state ;

extern struct element_struct   *element ;
extern struct node_struct      *node ;
extern struct boundary_struct  *boundary ;

extern char grid_file[MAX_CHAR] ;

extern char grid_file[MAX_CHAR] ;

/*********************************************************/
/**          READ the external .grd GRID_FILE           **/
/*********************************************************/

void read_spatial_grid()
{
	FILE *grid ;

	int e, n, f, dime ;

	char grid_file_name[MAX_CHAR] ;

        double dudu ;

	sprintf( grid_file_name, "./NEO_data/input/%s.grd", grid_file ) ;

	grid = fopen( grid_file_name, "r" ) ;
	if ( !grid ) printf( "ERROR: External spatial grid file was not found !!\n" ) ;
	
	fscanf( grid, "%d %d %d %d\n", &dime, &NE, &NN, &NBF ) ;
	fscanf( grid, "\n" ) ;
	
	for ( e = 0 ; e < NE ; e++ )
	      fscanf( grid, "%d %d %d\n", &element[e].node[0],&element[e].node[1],&element[e].node[2] ) ;
	
	 fscanf( grid, "\n" ) ;

	 for ( n = 0 ; n < NN ; n++ ){
	       fscanf( grid, "%le %le %d\n", &node[n].coordinate_0[0],&node[n].coordinate_0[1],&node[n].mate ) ;
       
         node[n].coordinate_old[0] = node[n].coordinate_0[0] ;
         node[n].coordinate_old[1] = node[n].coordinate_0[1] ;


         node[n].coordinate_new[0] = node[n].coordinate_0[0] ;
         node[n].coordinate_new[1] = node[n].coordinate_0[1] ;

         node[n].coordinate[0] = node[n].coordinate_0[0] ;
         node[n].coordinate[1] = node[n].coordinate_0[1] ;
        
         node[n].flag = 0 ;

                  	
		if ( initial_state == 9 )
      {
	// vortex advection: mates for mesh 1
	if(!strcmp(grid_file, "vort_t1")){
	  if ( n > 39 ) if ( n < 79 ) node[n].mate =  198 - n ;
	  if ( n > 119 ) if ( n < 159 ) node[n].mate =  198 - n ;
	}
	// vortex advection: mates for mesh 2
	else if(!strcmp(grid_file, "vort_t2")){
	  if ( n > 79 ) if ( n < 159 ) node[n].mate =  398 - n ;
	  if ( n > 239 ) if ( n < 319 ) node[n].mate =  398 - n ;
	}
	// vortex advection: mates for mesh 3
	else if(!strcmp(grid_file, "vort_t3")){
	  if ( n > 159 ) if ( n < 319 ) node[n].mate =  798 - n ;
	  if ( n > 479 ) if ( n < 639 ) node[n].mate =  798 - n ;
	}
	// vortex advection: mates for mesh 4
	else if(!strcmp(grid_file, "vort_t4")){
	  if ( n > 199 ) if ( n < 399 ) node[n].mate =  998 - n ;
	  if ( n > 599 ) if ( n < 799 ) node[n].mate =  998 - n ;
          node[n].coordinate[0] = 0.5*( node[n].coordinate[0] + 1. ) ; 
          node[n].coordinate[1] = 0.5*( node[n].coordinate[1] + 1. ) ; 
	}
	else if(!strcmp(grid_file, "vort_t5")){
	  if ( n > 266 ) if ( n < 533 ) node[n].mate =  1333 - n ;
	  if ( n > 800 ) if ( n < 1067 ) node[n].mate =  1333 - n ;
          node[n].coordinate[0] = 0.5*( node[n].coordinate[0] + 1. ) ; 
          node[n].coordinate[1] = 0.5*( node[n].coordinate[1] + 1. ) ; 
	} 
      }
      }
		
	 fscanf( grid, "\n" ) ;

	 for ( f = 0 ; f < NBF ; f++ ){
         fscanf( grid, "%d %d %d %d %d\n", &boundary[f].node[0], &boundary[f].node[NBN-1],
	                                         &boundary[f].types[0], &boundary[f].types[NBN-1],
	                                         &boundary[f].type ) ;

          if ( initial_state == 9 )
      {

               if(!strcmp(grid_file, "linear20-1"))
                 {
               boundary[f].types[0] = 1 ;
               boundary[f].types[1] = 1 ;
               boundary[f].type = 1 ;
                 }
              if(!strcmp(grid_file, "linear20-2"))
                 {
               boundary[f].types[0] = 1 ;
               boundary[f].types[1] = 1 ;
               boundary[f].type = 1 ;
                 }
		if(!strcmp(grid_file, "linear20-3"))
                 {
               boundary[f].types[0] = 1 ;
               boundary[f].types[1] = 1 ;
               boundary[f].type = 1 ;
                 }
		if(!strcmp(grid_file, "linear20-4"))
                 {
               boundary[f].types[0] = 1 ;
               boundary[f].types[1] = 1 ;
               boundary[f].type = 1 ;
                 }
      }

	}
	
	fclose ( grid ) ;
}

/****************************************************************/
/**      Read the HEADER of the external .grd GRID_FILE        **/
/****************************************************************/

void read_spatial_grid_header()
{
	FILE *grid_header ;

	char grid_file_name[MAX_CHAR] ;
	
	int dummy ;
	
	sprintf( grid_file_name, "./NEO_data/input/%s.grd", grid_file ) ;

	grid_header = fopen( grid_file_name, "r" ) ;
	if ( !grid_header ) printf( "ERROR: External spatial grid file was not found !! \n" ) ;

	fscanf( grid_header, "%d %d %d %d \n", &dummy, &NE, &NN, &NBF ) ;
	fclose ( grid_header ) ;
}
