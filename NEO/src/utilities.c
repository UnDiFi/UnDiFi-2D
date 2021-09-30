/***************************************************************************
                                 utilities.c
                                 -----------
  This is utilities: here are programmed all the utility functions used in
   the code, including small algebraic functions for matrix-matrix and
 matrix vector multiplications, matrix-matrix and vector-vector summation
                             -------------------
    Code developed by M. Ricchiuto, Inria - mario.ricchiuto@inria.fr
 ***************************************************************************/
#include "common.h"

extern struct node_struct *node ;
extern struct element_struct *element ;
extern struct boundary_struct  *boundary ;
extern struct boundary_nodes_struct   *b_nodes ;       

extern int NN, NE, size, iteration, initial_state, iteration ;
extern int convergence_variable, steady_check, NN, time_step, q_factor ;

extern long int *ipiv ;

extern double *temp_vector, *P_bar, **Right, **Left, *work, *phi_node, gm, *work_vector ;

extern double  sigma_barx, sigma_bary ;
extern double residual_norm, steady_norm, ref_norm, q_ratio ;
extern double pressure0, speed_of_sound0, total_enthalpy0, **W, **Z ;
extern double pressure1, speed_of_sound1, total_enthalpy1, ubar0, vbar0 ;
extern double a1[3], a2[3] , a3[3] ;
extern double b1[3], b2[3] , b3[3] ;

extern void A_times_v( double **, double *, double * ) ;
extern void Eigenvectors( double * ) ;       /* normals, variables, RIGHT EIGENVECTORS */
extern void initvec( double * ) ;

extern int dgetrf_(integer *, integer *, doublereal *, integer *, integer *, integer * ) ;
extern int dgetri_(integer *, doublereal *, integer *, integer *, doublereal *, integer *, integer * ) ;

extern void Z_2_P( double *, double * ) ;
extern void P_2_Z( double *, double * ) ;



/*****************************************/
/**            void function            **/
/*****************************************/

void void_function(){}

/*****************************************/
/**   Initialization of the solution    **/
/*****************************************/

void initialize_solution()
{
     int n, l ;

     for ( n = 0 ; n < NN ; n ++ )
           for ( l = 0 ; l < size ; l ++ )
               {
                 node[n].P[0][l] = node[n].P[2][l] ;
                 node[n].P[1][l] = node[n].P[2][l] ;

                 node[n].Z[0][l] = node[n].Z[2][l] ;
                 node[n].Z[1][l] = node[n].Z[2][l] ;
               }
}
/*****************************************/
/**      consistent variables           **/
/****************************************/

void consistent_vars( int e )
{
     int i, j, vv, itr ;

     itr = iteration - 1 ; 
    
     for ( j = 0 ; j < 3 ; j ++ )
         {
           vv = element[e].node[j] ;
           for ( i = 0 ; i < size ; i++ )
                 W[j][i] = a1[itr]*node[vv].P[0][i] + 
                           a2[itr]*node[vv].P[1][i] + 
                           a3[itr]*node[vv].P[2][i] ;

          P_2_Z( W[j], Z[j] ) ;
         }
}

/*****************************************/
/**    locally averaged variables       **/
/*****************************************/

void local_average( int e )
{
     int i, n1, n2, n3 ;
     double hh ;

     n1 = element[e].node[0] ;
     n2 = element[e].node[1] ;
     n3 = element[e].node[2] ;

     for ( i = 0 ; i < size ; i ++ )
           work_vector[i] = ( Z[0][i] + Z[1][i] + Z[2][i] )/3.0 ;

     pressure0 = work_vector[3] ;
     ubar0     = work_vector[1] ;
     vbar0     = work_vector[2] ;

     Z_2_P( work_vector, P_bar ) ;

     total_enthalpy0 = P_bar[3]/P_bar[0] + pressure0/P_bar[0] ;

     speed_of_sound0 = sqrt( gm*pressure0/P_bar[0] ) ;

     sigma_barx = ( node[n1].vel[0] + node[n2].vel[0] + node[n3].vel[0] )/3. ;
     sigma_bary = ( node[n1].vel[1] + node[n2].vel[1] + node[n3].vel[1] )/3. ; 

//     hh = fabs( ubar0 ) + fabs( vbar0 ) ;
//     if ( hh < 2.*1.e-3*speed_of_sound0 ) 
//          ubar0 = 1.e-3*speed_of_sound0 + 0.25*hh*hh/( 1.e-3*speed_of_sound0 ) ;
}




/*****************************************/
/**           residual norm             **/
/*****************************************/

void compute_norm()
{
	   int n, m, s ;
     double residual, norm, a, b, c ;

	   norm     = 0.0 ;	
	   residual = 0.0 ;

	    a = ( 2. + q_ratio )/( 1. + q_ratio ) ;
     b = ( 1. + q_ratio )/q_ratio ;
     c = 1./( q_ratio*( q_ratio + 1. ) ) ;

     steady_norm = 0.0 ;
	
	   m = convergence_variable ;

     if ( m == -1 )
	      {
          for ( n = 0 ; n < NN; n++ )
              {
	              for ( s = 0 ; s < size ; s ++ )
                      residual +=  fabs ( node[n].Res[s] )*1.0/( size*NN ) ;

                if ( steady_check )
                     for ( s = 0 ; s < size ; s ++ )
                           steady_norm += fabs( a*node[n].P[2][s] - b*node[n].P[1][s] 
					      + c*node[n].P[0][s] )*1.0/( size*NN ) ;
             }
        }
     else
        {
	        for ( n = 0 ; n < NN ; n++ )
	            {
 	              residual +=  fabs ( node[n].Res[m] )*1.0/NN ;

                if ( steady_check ) steady_norm += fabs( a*node[n].P[2][m] - 
				b*node[n].P[1][m] + c*node[n].P[0][m] )*1.0/NN ;
              }
        }

     residual_norm = log10( MAX( residual , 1.e-17 ) ) ;

     if ( steady_check )
          steady_norm   = log10( MAX( steady_norm , 1.e-17 ) ) ;

     if ( iteration == 1 ) ref_norm = residual_norm ;

}

/*****************************************/
/**          initialize matrix          **/
/*****************************************/

void initmat( double **A )
{
     int i, j ;

     for ( i = 0 ; i < size ; i ++ )
           for ( j = 0 ; j < size ; j ++ ) A[i][j] = 0.0 ;
}


/*****************************************/
/**          initialize vector          **/
/*****************************************/

void initvec( double *v )
{
     int i ;

     for ( i = 0 ; i < size ; i ++ ) v[i] = 0.0 ;
}

/*****************************************/
/**   add a matrix to an existing one   **/
/*****************************************/

void add2mat( double **A, double **B )
{
     int i, j ;

     for ( i = 0 ; i < size ; i ++ )
           for ( j = 0 ; j < size ; j ++ ) A[i][j] += B[i][j] ;
}

/*****************************************/
/**   add a vector to an existing one   **/
/*****************************************/

void add2vec( double *v, double *w )
{
     int i ;

     for ( i = 0 ; i < size ; i ++ ) v[i] += w[i] ;
}

/*****************************************/
/**        matrix-vector product        **/
/*****************************************/

void A_times_v( double **A, double *v, double *w )
{
     int i, j ;

     for ( i = 0 ; i < size ; i ++ )
         {
           w[i] = 0.0 ;

           for ( j = 0 ; j < size ; j ++ ) w[i] += A[i][j]*v[j] ;
         }
}

/*****************************************/
/**        matrix-matrix product        **/
/*****************************************/

void A_times_B( double **A, double **B, double **C )
{
     int i, j, k ;

     for ( i = 0 ; i < size ; i ++ )
           for ( j = 0 ; j < size ; j ++ )
               {
                 C[i][j] = 0.0 ;
                 for ( k = 0 ; k < size ; k ++ ) C[i][j] += A[i][k]*B[k][j] ;
               }
}

/****************************************************/
/**  update the local residual: no transformation  **/
/****************************************************/

void residual_update( int n )
{
     int i ;

     for ( i = 0 ; i < size ; i ++ )
           node[n].Res[i] += phi_node[i] ;
}

/***************************************************/
/**     decompose a residual in scalar waves      **/
/***************************************************/

void decompose( double *v, double *phiN_i, double *phi_i )
{
     Eigenvectors( v ) ;

     A_times_v( Left, phiN_i, phi_i ) ;
}

/********************************************************/
/** limit scalar residuals to get a HO positive scheme **/
/********************************************************/
/*
void limiter( double *phi, double **phiN, double *phiP, int node, int nodes )
{
     int  j, v ;
     double beta, arg ;

     for ( j = 0 ; j < size ; j ++ )
         {
           beta = 0.0 ;

           for ( v = 0 ; v < nodes ; v ++ )
               {
                 arg = phiN[v][j]*phi[j] ;

                 beta += ( 1.0 + SIGN( arg ) )*phiN[v][j] + small  ;
               }
           arg = phiN[node][j]*phi[j] ;

           beta = ( ( 1.0 + SIGN( arg ) )*phiN[node][j]  + small )/beta ;

           phiP[j] = beta*phi[j] ;
         }
} 
*/
void limiter( double *phi, double **phiN, double *phiP, int node, int nodes )
{
     int  j, v ;
     double beta, arg, beta1 ;

     for ( j = 0 ; j < size ; j ++ )
         {
           beta  = 0.0 ;
           beta1 = 0.0 ;

           phiP[j] = 0.0 ;

           if ( fabs( phi[j] ) > epsilon )
              {
                for ( v = 0 ; v < nodes ; v ++ )
                    {
                      arg = phiN[v][j]/phi[j] ;

                      beta1 += MAX( 0.0, arg ) ;
                    }
                 if ( beta1 > epsilon )
                    {
                      arg = phiN[node][j]/phi[j] ;

                      beta = MAX( 0.0, arg )/beta1 ;
                    }
                }

           phiP[j] = beta*phi[j] ;
         }
}

/*********************************************************/
/** reassemble a system residual from scalar components **/
/*********************************************************/

void recompose( double *v, double *phi_i, double *phiP_i )
{
     A_times_v( Right, phi_i, phiP_i ) ;
}





/***************************************************/
/**       matrix inversion using CLAPACK          **/
/***************************************************/

int inversion3( double **A )
{

	unsigned ret = 0 ;

	const double a = A[ 0 ][ 0 ] ;
	const double b = A[ 0 ][ 1 ] ;
	const double c = A[ 0 ][ 2 ] ;
	const double d = A[ 1 ][ 0 ] ;
	const double e = A[ 1 ][ 1 ] ;
	const double f = A[ 1 ][ 2 ] ;
	const double g = A[ 2 ][ 0 ] ;
	const double h = A[ 2 ][ 1 ] ;
	const double i = A[ 2 ][ 2 ] ;

        double D = - a * e * i 
		+ a * f * h 
		+ b * d * i 
		- b * f * g 
		- c * d * h 
		+ c * e * g ;

	if( fabs( D ) < 1.e-30 ) {
		ret = 1 ;
		D = 1 ;
	}

	A[ 0 ][ 0 ] = ( - e * i + f * h ) / D ;
	A[ 0 ][ 1 ] = (   b * i - c * h ) / D ;
	A[ 0 ][ 2 ] = ( - b * f + c * e ) / D ;  
	A[ 1 ][ 0 ] = (   d * i - f * g ) / D ;
	A[ 1 ][ 1 ] = ( - a * i + c * g ) / D ;
	A[ 1 ][ 2 ] = (   a * f - c * d ) / D ;
	A[ 2 ][ 0 ] = ( - d * h + e * g ) / D ;
	A[ 2 ][ 1 ] = (   a * h - b * g ) / D ;
	A[ 2 ][ 2 ] = ( - a * e + b * d ) / D ; 

	return ret ;

}

int inversion4( double **A )
{
	int info = 0 ;

	const double a = A[ 0 ][ 0 ] ;
	const double b = A[ 0 ][ 1 ] ;
	const double c = A[ 0 ][ 2 ] ;
	const double d = A[ 0 ][ 3 ] ;
	const double e = A[ 1 ][ 0 ] ;
	const double f = A[ 1 ][ 1 ] ;
	const double g = A[ 1 ][ 2 ] ;
	const double h = A[ 1 ][ 3 ] ;
	const double i = A[ 2 ][ 0 ] ;
	const double j = A[ 2 ][ 1 ] ;
	const double k = A[ 2 ][ 2 ] ;
	const double l = A[ 2 ][ 3 ] ;
	const double m = A[ 3 ][ 0 ] ;
	const double n = A[ 3 ][ 1 ] ;
	const double o = A[ 3 ][ 2 ] ;
	const double p = A[ 3 ][ 3 ] ;

// Hope to improve a bit Total Instruction Fetch Misses
	const double ah = a * h ;
	const double bi = b * i ;
	const double cj = c * j ;
	const double dk = d * k ;
	const double eo = e * o ;
	const double fm = f * m ;
	const double gn = g * n ;

	const double D = - a*f*k*p + a*f*l*o + a*g*j*p - a*gn*l 
		- ah*j*o + ah*k*n + b*e*k*p - b*eo*l 
		- bi*g*p + b*g*l*m + bi*h*o - b*h*k*m 
		- cj*e*p + c*e*l*n + c*f*i*p - c*fm*l 
		- c*h*i*n + cj*h*m + d*eo*j - dk*e*n 
		- d*f*i*o + dk*fm + d*gn*i - d*g*j*m ;

	if( fabs( D ) < 1.e-30 ) 
		return 1 ; 

	A[ 0 ][ 0 ]  = (- f*k*p + f*l*o + g*j*p - gn*l - h*j*o + h*k*n) / D ;
	A[ 0 ][ 1 ]  = (b*k*p - b*l*o - cj*p + c*l*n + d*j*o - dk*n)	 / D ;
	A[ 0 ][ 2 ]  = (- b*g*p + b*h*o + c*f*p - c*h*n - d*f*o + d*gn) / D ;
	A[ 0 ][ 3 ]  = (b*g*l - b*h*k - c*f*l + cj*h + dk*f - d*g*j)	 / D ;
	A[ 1 ][ 0 ]  = (e*k*p - eo*l - g*i*p + g*l*m + h*i*o - h*k*m)	 / D ;
	A[ 1 ][ 1 ]  = (- a*k*p + a*l*o + c*i*p - c*l*m - d*i*o + dk*m) / D ;
	A[ 1 ][ 2 ]  = (a*g*p - ah*o - c*e*p + c*h*m + d*eo - d*g*m)	 / D ;
	A[ 1 ][ 3 ]  = (- a*g*l + ah*k + c*e*l - c*h*i - dk*e + d*g*i) / D ;
	A[ 2 ][ 0 ]  = (- e*j*p + e*l*n + f*i*p - fm*l - h*i*n + h*j*m) / D ;
	A[ 2 ][ 1 ]  = (a*j*p - a*l*n - bi*p + b*l*m + d*i*n - d*j*m)	 / D ;
	A[ 2 ][ 2 ]  = (- a*f*p + ah*n + b*e*p - b*h*m - d*e*n + d*fm) / D ;
	A[ 2 ][ 3 ]  = (a*f*l - ah*j - b*e*l + bi*h + d*e*j - d*f*i)	 / D ;
	A[ 3 ][ 0 ]  = (eo*j - e*k*n - f*i*o + fm*k + gn*i - g*j*m)	 / D ;
	A[ 3 ][ 1 ]  = (- a*j*o + a*k*n + bi*o - b*k*m - c*i*n + cj*m) / D ;
	A[ 3 ][ 2 ]  = (a*f*o - a*gn - b*eo + b*g*m + c*e*n - c*fm)	 / D ;
	A[ 3 ][ 3 ]  = (- a*f*k + a*g*j + b*e*k - bi*g - cj*e + c*f*i) / D ; 


	return info ;

}

void invertmat( double **A, double **A_1 )
{

	int i, j ;

	for ( i = 0 ; i < size ; i ++ ) for ( j = 0 ; j < size ; j ++ ) A_1[i][j]=A[i][j] ;

	inversion4( A_1 ) ;


/***************************************************/
/**       matrix inversion using CLAPACK          **/
/***************************************************/
/*
 long int n,lda,lwork,info ;
      int i, j ;

      for ( i = 0 ; i < size ; i ++ ) for ( j = 0 ; j < size ; j ++ ) A_1[i][j]=A[i][j] ;

      n = size ;

      lda = size ;

      lwork = size ;

      dgetrf_(&n,&n,A_1[0],&lda,ipiv,&info) ;

      dgetri_(&n,A_1[0],&lda,ipiv,work,&lwork,&info) ; */
/**************************************/      
}

/***************************************************/
/**                Periodic BCs                   **/
/***************************************************/

void periodicity( int m )
{
     int i, mate, n ;
     double dummy ;

     n = b_nodes[m].node ;

     mate = node[n].mate ;

     if ( node[n].flag == 0 )
    {
     for ( i = 0 ; i < size ; i ++ )
         {
           dummy = node[n].Res[i] ;
           node[n].Res[i] += node[mate].Res[i] ;
           node[mate].Res[i] += dummy ;
         }
     node[n].flag = 1 ;
     node[mate].flag = 1 ;
	}
}


/***************************************************/
/**      Standard supeersonic Inlet BCs           **/
/***************************************************/

void supersonic_residual_std( int v )
{
     int i, n, n1, n2 ;
     
     n1 = boundary[v].node[0] ;    
     n2 = boundary[v].node[1] ;    

    if ( node[n1].flag == 0 )
{
       for ( i = 0 ; i < size ; i ++ )
           node[n1].Res[i] = 0.0 ;

       node[n1].flag = 1 ;
}

    if ( node[n2].flag == 0 )
{
       for ( i = 0 ; i < size ; i ++ )
           node[n2].Res[i] = 0.0 ;

       node[n2].flag = 1 ;
}


}

