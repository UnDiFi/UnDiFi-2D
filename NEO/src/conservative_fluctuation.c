/***************************************************************************
                          conservative_fluctuation.c
                          --------------------------
       This is conservative_fluctuation: here the contour integral of
                    the (space-time) fluxes is performed
                             -------------------
    begin                : Fri May 10 2002
    copyright            : (C) 2002 by Mario Ricchiuto
    email                : ricchiut@vki.ac.be
 ***************************************************************************/
#include "common.h"

extern struct element_struct *element ;
extern struct node_struct *node ;
extern struct numerical_int_struct numerical_int ;

extern int size, face_q_pts, volume_q_pts, iteration ;
extern double sigma_barx, sigma_bary ;
extern double *temp_vector, *temp_normal, *F_loc, *PHI_loc, *X_temp, ***K0, ***K1, *work_vector ;
extern double **dU, **W, dt, **phi_t_w, *phi_w, *phi_a_w0, *phi_a_w1, *FLUX  ;
extern double  *work_vector1, *work_vector0, q_ratio ;
extern double a1[3], a2[3], a3[3] ;
extern double b1[3], b2[3], b3[3] ;


extern void ( *fluctuation2 )( int ) ;
extern void initvec( double * ) ;
extern void A_times_v( double **, double *, double * ) ;
extern void add2vec( double *, double * ) ;
extern void ad_flux( double *, double * ) ;
extern void Z_2_P( double *, double * ) ;
extern void P_2_Z( double *, double * ) ;

/*********************************************/
/*********************************************/
/**      unsteady homogeneous systems       **/
/*********************************************/
/*********************************************/

void STfluctuation2( int e )
{
     int i, v, itr, i1, i2, i3, vv ;
     double V_T ;

     itr = iteration - 1 ;

     V_T = element[e].volume/3.0 ;

     for ( v = 0 ; v < 3 ; v ++ )
         {
           vv = element[e].node[v] ;
           for ( i = 0 ; i < size ; i ++ )
                 dU[v][i] = V_T*( b3[itr]*node[vv].P[2][i] + b2[itr]*node[vv].P[1][i] + b1[itr]*node[vv].P[0][i] ) ;
         }

     fluctuation2( e ) ;

     for ( i = 0 ; i < size ; i ++ )
          phi_w[i] = dU[0][i] + dU[1][i] + dU[2][i] + dt*phi_a_w1[i] ;
}

void fluctuationRK1( int e )
{
     int f, i, i1, i2, i0, itr, v, vvv ;
     double sigman ;  

     itr = iteration - 1 ;

     initvec( work_vector ) ;

     i0 = element[e].node[0] ;
     i1 = element[e].node[1] ;
     i2 = element[e].node[2] ;

/**********************************************/
/**    		ADVECTIVE PART               **/
/**********************************************/

/**********************************************/
/**        edge 0: nodes 2 and 3             **/
/**********************************************/

     temp_normal[0] = -element[e].normal[0][0] ;
     temp_normal[1] = -element[e].normal[0][1] ;

// time tn

     for ( f = 0 ; f < face_q_pts; f ++ )
         {
           for ( i = 0 ; i < size ; i ++ )
                 temp_vector[i] = numerical_int.face_coordinate[f][0]*node[i1].P[0][i] +
                                  numerical_int.face_coordinate[f][1]*node[i2].P[0][i] ;

	   P_2_Z( temp_vector, work_vector0 ) ;

           ad_flux( temp_normal, work_vector0 ) ;

           for( i = 0 ; i < size ; i ++ )
                work_vector[i] += numerical_int.face_weight[f]*FLUX[i] ;
         }

/**********************************************/
/**           edge 1: nodes 3 and 1          **/
/**********************************************/

     temp_normal[0] = -element[e].normal[1][0] ;
     temp_normal[1] = -element[e].normal[1][1] ;

     for ( f = 0 ; f < face_q_pts; f ++ )
         {
           for ( i = 0 ; i < size ; i ++ )
                 temp_vector[i] = numerical_int.face_coordinate[f][0]*node[i2].P[0][i] +
                                  numerical_int.face_coordinate[f][1]*node[i0].P[0][i] ;

	   P_2_Z( temp_vector, work_vector0 ) ;

           ad_flux( temp_normal, work_vector0 ) ;

           for( i = 0 ; i < size ; i ++ )
                work_vector[i] += numerical_int.face_weight[f]*FLUX[i] ;
         }

/**********************************************/
/**         edge 2:  nodes 1 and 2           **/
/**********************************************/

     temp_normal[0] = -element[e].normal[2][0] ;
     temp_normal[1] = -element[e].normal[2][1] ;

     for ( f = 0 ; f < face_q_pts; f ++ )
         {
           for ( i = 0 ; i < size ; i ++ )
                 temp_vector[i] = numerical_int.face_coordinate[f][0]*node[i0].P[0][i] +
                                  numerical_int.face_coordinate[f][1]*node[i1].P[0][i] ;

	    P_2_Z( temp_vector, work_vector0 ) ;

           ad_flux( temp_normal, work_vector0 ) ;

           for( i = 0 ; i < size ; i ++ )
                work_vector[i] += numerical_int.face_weight[f]*FLUX[i] ;
         }

    for ( i = 0 ; i < size ; i ++ ) 
         phi_a_w1[i] = work_vector[i] ;

/**********************************************/
/**         	   ALE PART	             **/
/**********************************************/

    for ( v = 0 ; v < 3 ; v ++ )
      {

       vvv = element[e].node[v] ;
       sigman = 0.5*( sigma_barx*element[e].normal[v][0] + sigma_bary*element[e].normal[v][1] ) ;

       for ( i = 0; i < size ; i++ ) phi_a_w1[i] -= sigman* node[vvv].P[0][i] ;

      }
      

}


void fluctuationRK2( int e )
{
     int f, i, i1, i2, i0, itr, v, vvv ;
     double sigman ;


     itr = iteration - 1 ;

     switch ( itr )
        {    
         case 0 : 
         fluctuationRK1( e ) ;
         break ;
         case 1 :
         initvec( work_vector ) ;

         i0 = element[e].node[0] ;
         i1 = element[e].node[1] ;
         i2 = element[e].node[2] ;

/**********************************************/
/**    		ADVECTIVE PART               **/
/**********************************************/

/**********************************************/
/**        edge 0: nodes 2 and 3             **/
/**********************************************/

         temp_normal[0] = -element[e].normal[0][0] ;
         temp_normal[1] = -element[e].normal[0][1] ;

// time tn

         for ( f = 0 ; f < face_q_pts; f ++ )
             {
               for ( i = 0 ; i < size ; i ++ )
                     temp_vector[i] = numerical_int.face_coordinate[f][0]*node[i1].P[0][i] +
                                      numerical_int.face_coordinate[f][1]*node[i2].P[0][i] ;

	        P_2_Z( temp_vector, work_vector0 ) ;

                ad_flux( temp_normal, work_vector0 ) ;

                for( i = 0 ; i < size ; i ++ )
                     work_vector[i] += 0.5*numerical_int.face_weight[f]*FLUX[i] ;
             }

// time t1

         for ( f = 0 ; f < face_q_pts; f ++ )
             {
               for ( i = 0 ; i < size ; i ++ )
                     temp_vector[i] = numerical_int.face_coordinate[f][0]*node[i1].P[2][i] +
                                      numerical_int.face_coordinate[f][1]*node[i2].P[2][i] ;

	       P_2_Z( temp_vector, work_vector0 ) ;

               ad_flux( temp_normal, work_vector0 ) ;

               for( i = 0 ; i < size ; i ++ )
                    work_vector[i] += 0.5*numerical_int.face_weight[f]*FLUX[i] ;
              }

/**********************************************/
/**           edge 1: nodes 3 and 1          **/
/**********************************************/

         temp_normal[0] = -element[e].normal[1][0] ;
         temp_normal[1] = -element[e].normal[1][1] ;

// time tn

         for ( f = 0 ; f < face_q_pts; f ++ )
             {
               for ( i = 0 ; i < size ; i ++ )
                     temp_vector[i] = numerical_int.face_coordinate[f][0]*node[i2].P[0][i] +
                                      numerical_int.face_coordinate[f][1]*node[i0].P[0][i] ;

	       P_2_Z( temp_vector, work_vector0 ) ;

               ad_flux( temp_normal, work_vector0 ) ;

               for( i = 0 ; i < size ; i ++ )
                    work_vector[i] += 0.5*numerical_int.face_weight[f]*FLUX[i] ;
             }

// time t1

         for ( f = 0 ; f < face_q_pts; f ++ )
             {
               for ( i = 0 ; i < size ; i ++ )
                     temp_vector[i] = numerical_int.face_coordinate[f][0]*node[i2].P[2][i] +
                                      numerical_int.face_coordinate[f][1]*node[i0].P[2][i] ;

	        P_2_Z( temp_vector, work_vector0 ) ;

                ad_flux( temp_normal, work_vector0 ) ;

                for( i = 0 ; i < size ; i ++ )
                     work_vector[i] += 0.5*numerical_int.face_weight[f]*FLUX[i] ;
             }
 
/**********************************************/
/**         edge 2:  nodes 1 and 2           **/
/**********************************************/

     temp_normal[0] = -element[e].normal[2][0] ;
     temp_normal[1] = -element[e].normal[2][1] ;

// time tn

         for ( f = 0 ; f < face_q_pts; f ++ )
             {
               for ( i = 0 ; i < size ; i ++ )
                     temp_vector[i] = numerical_int.face_coordinate[f][0]*node[i0].P[0][i] +
                                      numerical_int.face_coordinate[f][1]*node[i1].P[0][i] ;

	        P_2_Z( temp_vector, work_vector0 ) ;

                ad_flux( temp_normal, work_vector0 ) ;

                for( i = 0 ; i < size ; i ++ )
                     work_vector[i] += 0.5*numerical_int.face_weight[f]*FLUX[i] ;
              }

// time t1

         for ( f = 0 ; f < face_q_pts; f ++ )
             {
               for ( i = 0 ; i < size ; i ++ )
                     temp_vector[i] = numerical_int.face_coordinate[f][0]*node[i0].P[2][i] +
                                      numerical_int.face_coordinate[f][1]*node[i1].P[2][i] ;

	        P_2_Z( temp_vector, work_vector0 ) ;

                ad_flux( temp_normal, work_vector0 ) ;

                for( i = 0 ; i < size ; i ++ )
                     work_vector[i] += 0.5*numerical_int.face_weight[f]*FLUX[i] ;
              }

// Final update

          for ( i = 0 ; i < size ; i ++ ) 
                phi_a_w1[i] = work_vector[i] ;

/**********************************************/
/**         	   ALE PART	             **/
/**********************************************/

    for ( v = 0 ; v < 3 ; v ++ )
      {

       vvv = element[e].node[v] ;
       sigman = 0.5*( sigma_barx*element[e].normal[v][0] + sigma_bary*element[e].normal[v][1] ) ;

       for ( i = 0; i < size ; i++ ) phi_a_w1[i] -= 0.5*sigman*( node[vvv].P[0][i] + node[vvv].P[2][i] ) ;

      }


// end switch  
 
       }
}

void fluctuationRK3( int e )
{
     int f, i, i1, i2, i0, itr ;


     itr = iteration - 1 ;

     switch ( itr )
        {    
         case 0 : 
         fluctuationRK1( e ) ;
         break ;
         case 1 :
         fluctuationRK2( e ) ;
         for ( i = 0 ; i < size ; i ++ ) 
               phi_a_w1[i] *= 0.5 ;
         break ;
         case 2 :
         initvec( work_vector ) ;

         i0 = element[e].node[0] ;
         i1 = element[e].node[1] ;
         i2 = element[e].node[2] ;

/**********************************************/
/**        edge 0: nodes 2 and 3             **/
/**********************************************/

         temp_normal[0] = -element[e].normal[0][0] ;
         temp_normal[1] = -element[e].normal[0][1] ;

// time tn

         for ( f = 0 ; f < face_q_pts; f ++ )
             {
               for ( i = 0 ; i < size ; i ++ )
                     temp_vector[i] = numerical_int.face_coordinate[f][0]*node[i1].P[0][i] +
                                      numerical_int.face_coordinate[f][1]*node[i2].P[0][i] ;

	        P_2_Z( temp_vector, work_vector0 ) ;

                ad_flux( temp_normal, work_vector0 ) ;

                for( i = 0 ; i < size ; i ++ )
                     work_vector[i] += numerical_int.face_weight[f]*FLUX[i]/6. ;
             }

// time t2

         for ( f = 0 ; f < face_q_pts; f ++ )
             {
               for ( i = 0 ; i < size ; i ++ )
                     temp_vector[i] = numerical_int.face_coordinate[f][0]*node[i1].P[2][i] +
                                      numerical_int.face_coordinate[f][1]*node[i2].P[2][i] ;

	       P_2_Z( temp_vector, work_vector0 ) ;

               ad_flux( temp_normal, work_vector0 ) ;

               for( i = 0 ; i < size ; i ++ )
                    work_vector[i] += 2.*numerical_int.face_weight[f]*FLUX[i]/3. ;
              }
// time t1

         for ( f = 0 ; f < face_q_pts; f ++ )
             {
               for ( i = 0 ; i < size ; i ++ )
                     temp_vector[i] = numerical_int.face_coordinate[f][0]*node[i1].P[1][i] +
                                      numerical_int.face_coordinate[f][1]*node[i2].P[1][i] ;

	       P_2_Z( temp_vector, work_vector0 ) ;

               ad_flux( temp_normal, work_vector0 ) ;

               for( i = 0 ; i < size ; i ++ )
                    work_vector[i] += numerical_int.face_weight[f]*FLUX[i]/6. ;
              }



/**********************************************/
/**           edge 1: nodes 3 and 1          **/
/**********************************************/

         temp_normal[0] = -element[e].normal[1][0] ;
         temp_normal[1] = -element[e].normal[1][1] ;

// time tn

         for ( f = 0 ; f < face_q_pts; f ++ )
             {
               for ( i = 0 ; i < size ; i ++ )
                     temp_vector[i] = numerical_int.face_coordinate[f][0]*node[i2].P[0][i] +
                                      numerical_int.face_coordinate[f][1]*node[i0].P[0][i] ;

	       P_2_Z( temp_vector, work_vector0 ) ;

               ad_flux( temp_normal, work_vector0 ) ;

               for( i = 0 ; i < size ; i ++ )
                    work_vector[i] += numerical_int.face_weight[f]*FLUX[i]/6. ;
             }

// time t2

         for ( f = 0 ; f < face_q_pts; f ++ )
             {
               for ( i = 0 ; i < size ; i ++ )
                     temp_vector[i] = numerical_int.face_coordinate[f][0]*node[i2].P[2][i] +
                                      numerical_int.face_coordinate[f][1]*node[i0].P[2][i] ;

	        P_2_Z( temp_vector, work_vector0 ) ;

                ad_flux( temp_normal, work_vector0 ) ;

                for( i = 0 ; i < size ; i ++ )
                     work_vector[i] += 2.*numerical_int.face_weight[f]*FLUX[i]/3. ;
             }
 
// time t1

         for ( f = 0 ; f < face_q_pts; f ++ )
             {
               for ( i = 0 ; i < size ; i ++ )
                     temp_vector[i] = numerical_int.face_coordinate[f][0]*node[i2].P[1][i] +
                                      numerical_int.face_coordinate[f][1]*node[i0].P[1][i] ;

	       P_2_Z( temp_vector, work_vector0 ) ;

               ad_flux( temp_normal, work_vector0 ) ;

               for( i = 0 ; i < size ; i ++ )
                    work_vector[i] += numerical_int.face_weight[f]*FLUX[i]/6. ;
             }

/**********************************************/
/**         edge 2:  nodes 1 and 2           **/
/**********************************************/

     temp_normal[0] = -element[e].normal[2][0] ;
     temp_normal[1] = -element[e].normal[2][1] ;

// time tn

         for ( f = 0 ; f < face_q_pts; f ++ )
             {
               for ( i = 0 ; i < size ; i ++ )
                     temp_vector[i] = numerical_int.face_coordinate[f][0]*node[i0].P[0][i] +
                                      numerical_int.face_coordinate[f][1]*node[i1].P[0][i] ;

	        P_2_Z( temp_vector, work_vector0 ) ;

                ad_flux( temp_normal, work_vector0 ) ;

                for( i = 0 ; i < size ; i ++ )
                     work_vector[i] += numerical_int.face_weight[f]*FLUX[i]/6. ;
              }

// time t2

         for ( f = 0 ; f < face_q_pts; f ++ )
             {
               for ( i = 0 ; i < size ; i ++ )
                     temp_vector[i] = numerical_int.face_coordinate[f][0]*node[i0].P[2][i] +
                                      numerical_int.face_coordinate[f][1]*node[i1].P[2][i] ;

	        P_2_Z( temp_vector, work_vector0 ) ;

                ad_flux( temp_normal, work_vector0 ) ;

                for( i = 0 ; i < size ; i ++ )
                     work_vector[i] += 2.*numerical_int.face_weight[f]*FLUX[i]/3. ;
              }

// time t1

         for ( f = 0 ; f < face_q_pts; f ++ )
             {
               for ( i = 0 ; i < size ; i ++ )
                     temp_vector[i] = numerical_int.face_coordinate[f][0]*node[i0].P[1][i] +
                                      numerical_int.face_coordinate[f][1]*node[i1].P[1][i] ;

	        P_2_Z( temp_vector, work_vector0 ) ;

                ad_flux( temp_normal, work_vector0 ) ;

                for( i = 0 ; i < size ; i ++ )
                     work_vector[i] += numerical_int.face_weight[f]*FLUX[i]/6. ;
              }

/**************************************************/
/**************************************************/
// Final update
/**************************************************/
/**************************************************/

          for ( i = 0 ; i < size ; i ++ ) 
                phi_a_w1[i] = work_vector[i] ;
          break ;
// end switch  
 
       }
}

/***************************************************************************
 *                                                                         *
 *   This program has been developed at the von Karman Institute for Fluid *
 *   Dynamics by Mario Ricchiuto. Any change or update must be done in a   *
 *   clean way and it must be clearly commented.                           *
 *   Mario Ricchiuto                                                       *
 *                                              23/04/2002                 *
 ***************************************************************************/
