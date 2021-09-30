/***************************************************************************
                          common_variables.h
                               --------
This is common_varibles: it contains and defines all the variables (dangerously)
                 accessible from anywhere in the program.
                              -------------------
    Code developed by M. Ricchiuto, Inria - mario.ricchiuto@inria.fr
 ***************************************************************************/

struct element_struct   *element ;
struct node_struct      *node ;
struct boundary_struct  *boundary ;
struct numerical_int_struct numerical_int ;
struct boundary_nodes_struct  *b_nodes ;  


int model, initial_state, scheme, gauss,  face_q_pts, volume_q_pts ;
int iteration_max, convergence_variable, time_step_max, stop, out_format, start_up ;
int movie, save_solution, size, NN, NE, NBF, NBN ;
int vertex, thermal_vars, problem, iteration, counter_info ;
int time_step, iteration, steady_check, info, lump_type  ;

long int *ipiv ;

char initial_solution_file[MAX_CHAR], grid_file[MAX_CHAR] ;
char out_file[MAX_CHAR] ;
char oneD_file[MAX_CHAR] ;

double zero_r, pressure0, speed_of_sound0, total_enthalpy0, c_tau, ref_vel ;
double pressure1, speed_of_sound1, total_enthalpy1, alpha  ;
double p_out, Mach_in, alpha_in, rho_in  ;
double Mach_inf, alpha_inf, rho_inf , p_inf ;
 
double diag0, diag1, ubar0, vbar0 ;
double a1[3], a2[3], a3[3] ;
double b1[3], b2[3], b3[3] ;

double gm, delta_t, shield_factor, CFL, t_max, residual_treshold, residual_norm, steady_norm ;
double time, q_ratio ;
double ref_norm, lim_norm ;

double sigma_barx, sigma_bary ;
double *vel, *work_vector, *work_vector0, *work_vector1, *work_vector2, dt, **temp_mat1 ;
double *temp_vector, *temp_normal, **P_bar, *normal, *Lambda,  **temp_mat ;
double *phi_a_w0, *phi_a_w1, **phi_t_w, *F_loc, *U_C,  *work, *PHI_loc, *phi_w, *phi1_w ;
double **W, **Z, ***U, **PHI_d, **PHIN_d, **PHIN_scal, *phi_node, *FLUX, **dU ;
double ***K0, ***K1, ***K_i_n0, ***K_i_n1 ;
double **Right, **Left, ***K_i, ***K_i_p0, ***K_i_p1, **sum_K, **sum_K_1, *X_temp ;


void ( *distribute_fluctuation )( int ) ;
void ( *fluctuation2 )( int ) ;
void ( *move_the_mesh )( int ) ;
void ( *supersonic_residual )( int ) ; /* Face global index  */



/***************************************************************************
 *                                                                         *
 *   This program has been developed at the von Karman Institute for Fluid *
 *   Dynamics by Mario Ricchiuto. Any change or update must be done in a   *
 *   clean way and it must be clearly commented.                           *
 *   Mario Ricchiuto                                                       *
 *                                              23/04/2002                 *
 ***************************************************************************/
