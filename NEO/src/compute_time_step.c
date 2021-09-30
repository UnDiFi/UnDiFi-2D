/***************************************************************************
                             compute_time_step.c
                             -------------------
    This is compute_time_step:  local and global time steps computed
      from a CFL condition based on the wost case contraint (LF scheme)
                              -------------------
    Code developed by M. Ricchiuto, Inria - mario.ricchiuto@inria.fr
 ***************************************************************************/
#include "common.h"

extern int NE, size, model, NN , NBF ;

extern struct element_struct *element ;
extern struct node_struct *node ;

extern struct boundary_struct *boundary ;

extern int initial_state, time_step ;

extern double *temp_vector, *normal, *Lambda ;

extern double dt, shield_factor, time, t_max, q_ratio, gm, zero_r ;

extern void Waves( double * ) ; /* normals, variables, WAVES */
extern void local_average( int ) ;

void compute_time_step()
{
  int e, v, j, w, n, n1, n2, n3 , f ;
  double ratio, delta_t0, volume, max_v, max_c, max_h, speed ;
  double sigma_barx, sigma_bary ;
  FILE *out;

  for ( n = 0 ; n < NN ; n ++ ) node[n].dtau = 0. ;

  q_ratio = dt ;

  if ( time_step == 1 )  q_ratio = 1. ;
     
  delta_t0 = infinity ;

  for ( e = 0 ; e < NE ; e ++ )
    {

      n1 = element[e].node[0] ;
      n2 = element[e].node[1] ;
      n3 = element[e].node[2] ;

      sigma_barx = ( node[n1].vel[0] + node[n2].vel[0] + node[n3].vel[0] )/3. ;
      sigma_bary = ( node[n1].vel[1] + node[n2].vel[1] + node[n3].vel[1] )/3. ;

      volume = element[e].volume ;

      max_v = max_c = max_h = 0. ;

      for ( v = 0 ; v < 3 ; v ++ )
	{
	  for ( j = 0 ; j < 2 ; j ++ )
	    normal[j] = element[e].normal[v][j] ;
 
	  volume = sqrt( normal[0]*normal[0] + normal[1]*normal[1] ) ;

	  if ( volume > max_h ) max_h = volume ;

	  n = element[e].node[v] ;

	  volume = sqrt( (node[n].Z[2][1]-sigma_barx)*(node[n].Z[2][1]-sigma_barx) + 
			 + (node[n].Z[2][2]-sigma_bary)*(node[n].Z[2][2]-sigma_bary) ) ;

	  if ( volume > max_v ) max_v = volume ;

	  volume = node[n].Z[2][0] ; if ( volume < zero_r ) volume = zero_r ;
 
	  volume = sqrt( gm*node[n].Z[2][3]/volume ) ; 
 
	  if ( volume > max_c ) max_c = volume ;
	}

      speed = 0.5*( max_v + max_c )*max_h ;    
 
      node[element[e].node[0]].dtau += speed ;
      node[element[e].node[1]].dtau += speed ;
      node[element[e].node[2]].dtau += speed ;
    }

  for ( f = 0 ; f < NBF ; f ++ )
    {

      n1 = boundary[f].node[0] ;
      n2 = boundary[f].node[1] ;

      sigma_barx = ( node[n1].vel[0] + node[n2].vel[0] )/2. ;
      sigma_bary = ( node[n1].vel[1] + node[n2].vel[1] )/2. ;

      volume = sqrt(boundary[f].normal[0]*boundary[f].normal[0] + boundary[f].normal[1]*boundary[f].normal[1]) ;

      max_h = volume ;
     
      max_v = max_c =  0. ;

      for ( v = 0 ; v < 2 ; v ++ )
	{

	  n = boundary[f].node[v] ;

	  volume = sqrt( (node[n].Z[2][1]-sigma_barx)*(node[n].Z[2][1]-sigma_barx) + 
			 + (node[n].Z[2][2]-sigma_bary)*(node[n].Z[2][2]-sigma_bary) ) ;

	  if ( volume > max_v ) max_v = volume ;

	  volume = node[n].Z[2][0] ; if ( volume < zero_r ) volume = zero_r ;
 
	  volume = sqrt( gm*node[n].Z[2][3]/volume ) ; 
 
	  if ( volume > max_c ) max_c = volume ;
	}

      speed = 0.5*( max_v + max_c )*max_h ;    
 
      node[boundary[f].node[0]].dtau += speed ;
      node[boundary[f].node[1]].dtau += speed ;
    }

  for ( n = 0 ; n < NN ; n ++ )
    {
      ratio = node[n].vol/( node[n].dtau + 1.e-20 ) ;
      if ( delta_t0 > ratio ) delta_t0 = ratio ;
    }

  delta_t0 *= shield_factor ;

  if ( time + delta_t0 > t_max ) delta_t0 = t_max - time ;

  /****************************************************************************/
  /************ below: sampling for the Mach 3 wind tunnel ********************/
  /****************************************************************************/
  if ( initial_state ==3 )
    {
      if ( time < 0.5 )
	{
	  if ( time + delta_t0 >= 0.5 ) delta_t0 = 0.5 - time ;
	}
      
      if ( time < 1. )
	{
	  if ( time + delta_t0 >= 1. ) delta_t0 = 1. - time ;
	}
      
      if ( time < 1.5 )
	{
	  if ( time + delta_t0 >= 1.5 ) delta_t0 = 1.5 - time ;
	}
      
      if ( time < 2. )
	{
	  if ( time + delta_t0 >= 2. ) delta_t0 = 2. - time ;
	}
      
      if ( time < 2.5 )
	{
	  if ( time + delta_t0 >= 2.5 ) delta_t0 = 2.5 - time ;
	}
      
      if ( time < 3. )
	{
	  if ( time + delta_t0 >= 3. ) delta_t0 = 3. - time ;
	}
      
      if ( time < 3.5 )
	{
	  if ( time + delta_t0 >= 3.5 ) delta_t0 = 3.5 - time ;
	}
    }
  /****************************************************************************/
  /************ above: sampling for the MAch 3 wind tunnel ********************/
  /****************************************************************************/
/* Campoli
  out = fopen("timesteps.dat","r");
  fscanf(out,"%le ", &dt);
  fclose(out);
*/
    dt = delta_t0 ;

  time += dt ;

  q_ratio /= dt ;
/* Campoli
  out = fopen("timesteps.dat","w");
  fprintf(out,"%le ", dt);
  fclose(out);
*/
/*printf("Timestep is %le ", dt);*/
}

