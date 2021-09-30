/***************************************************************************
                              input_parameters.c
                               ----------------
This is input_parameters: it parses the inputfile and and saves all the necessary 
         parameters  needed to set up and run the computation   
                              -------------------
    Code developed by M. Ricchiuto, Inria - mario.ricchiuto@inria.fr
 ***************************************************************************/
#include "common.h"

extern int model, initial_state, scheme, zero_r ;
extern int iteration_max, convergence_variable, time_step_max, stop, movie, info, lump_type ;
extern int save_solution, out_format, problem, start_up, steady_check, out_lev ;

extern double shield_factor, CFL, residual_treshold, t_max, gm, lim_norm, ref_vel ;

extern char initial_solution_file[MAX_CHAR] ;
extern char grid_file[MAX_CHAR] ;
extern char out_file[MAX_CHAR] ;

void input_parameters()
{
      int done ;
      FILE *inputfile ;

      inputfile = fopen( "./NEO_data/textinput/inputfile-exp.txt", "r" ) ;
      if ( !inputfile ) printf( "ERROR: inputfile-exp.txt not found!\n" ) ;

      printf( "          *************************************\n" ) ;
      printf( "          **    Reading the inputfile.....   **\n" ) ;
      printf( "          *************************************\n" ) ;
      printf( "\n" ) ;

      fscanf( inputfile, "***************************************************\n" ) ;
      fscanf( inputfile, "***************************************************\n" ) ;
      fscanf( inputfile, "**                                               **\n" ) ;
      fscanf( inputfile, "**           NEO's basic parameters              **\n" ) ;
      fscanf( inputfile, "**                                               **\n" ) ;
      fscanf( inputfile, "***************************************************\n" ) ;
      fscanf( inputfile, "***************************************************\n" ) ;
      fscanf( inputfile, "***************************************************\n" ) ;
      fscanf( inputfile, "**                  PHYSICS                      **\n" ) ;
      fscanf( inputfile, "***************************************************\n" ) ;
      fscanf( inputfile, "   Model                                     %d    \n", &model ) ;
      fscanf( inputfile, "   Initial state                             %d    \n", &initial_state ) ;
      fscanf( inputfile, "   Initial solution file                     %s    \n", initial_solution_file ) ;
      fscanf( inputfile, "   Problem                                   %d    \n", &problem ) ;
      fscanf( inputfile, "   Gamma                                     %le    \n", &gm ) ;
      fscanf( inputfile, "   Zero density                              %le    \n", &zero_r ) ;
      fscanf( inputfile, "***************************************************\n" ) ;
      fscanf( inputfile, "**                  GEOMETRY                     **\n" ) ;
      fscanf( inputfile, "***************************************************\n" ) ;
      fscanf( inputfile, "   External grid file                        %s    \n", grid_file ) ;
      fscanf( inputfile, "***************************************************\n" ) ;
      fscanf( inputfile, "**                  NUMERICS                     **\n" ) ;
      fscanf( inputfile, "***************************************************\n" ) ;
      fscanf( inputfile, "   0(LDA),1(LLF),2(LDAN),3(LWLF),4(LW),5(N)     %d \n", &scheme ) ;
      fscanf( inputfile, "   Lumping type: 0 (selective),1 (global)       %d \n", &lump_type ) ;
      fscanf( inputfile, "   Shield factor                               %le \n", &shield_factor );
      fscanf( inputfile, "***************************************************\n" ) ;
      fscanf( inputfile, "**                  ALGEBRA                      **\n" ) ;
      fscanf( inputfile, "***************************************************\n" ) ;
      fscanf( inputfile, "   Convergence treshold                     %le    \n", &residual_treshold ) ;
      fscanf( inputfile, "   Convergence limit                        %le    \n", &lim_norm ) ;
      fscanf( inputfile, "   Check Steady State                        %d    \n", &steady_check ) ;
      fscanf( inputfile, "   Maximum number of iterations              %d    \n", &iteration_max ) ;
      fscanf( inputfile, "   Convergence variable (0->Neq.s-1)         %d    \n", &convergence_variable ) ;
      fscanf( inputfile, "***************************************************\n" ) ;
      fscanf( inputfile, "**                   STOP                        **\n" ) ;
      fscanf( inputfile, "***************************************************\n" ) ;
      fscanf( inputfile, "   Final time                               %le    \n", &t_max ) ;
      fscanf( inputfile, "   Maximum number of time steps              %d    \n", &time_step_max ) ;
      fscanf( inputfile, "   Stop                                      %d    \n", &stop ) ;
      fscanf( inputfile, "***************************************************\n" ) ;
      fscanf( inputfile, "**                  OUTPUT                       **\n" ) ;
      fscanf( inputfile, "***************************************************\n" ) ;
      fscanf( inputfile, "   Movie (0/1)                               %d    \n", &movie ) ;
      fscanf( inputfile, "   Steps to save                             %d    \n", &save_solution ) ;
      fscanf( inputfile, "   Information Frequency                     %d    \n", &info ) ;
      fscanf( inputfile, "   Output format  (0/1)                      %d    \n", &out_format ) ;
      fscanf( inputfile, "   Output file                               %s    \n", out_file ) ;
      fscanf( inputfile, "***************************************************\n" ) ;
      fscanf( inputfile, "***************************************************\n" ) ;
      fscanf( inputfile, "                      %d                           \n", &done ) ;
      fscanf( inputfile, "***************************************************\n" ) ;
      fscanf( inputfile, "***************************************************\n" ) ;
      fscanf( inputfile, "***************************************************\n" ) ;

      if ( done != 8576 ) printf( "ERROR: inputfile was not read correctly!!" ) ;
      fclose( inputfile ) ;

      if ( initial_state == 9 ) t_max = 1./6. ;
}

