      DOUBLE PRECISION RESMAX(MAXNOFVAR,2),RESL2(MAXNOFVAR,2),
     &DELMAX(MAXNOFVAR,2),DELL2(MAXNOFVAR,2),RESL20(2),RESMAX0(2),
     +CFL(2),CFLMAX(2),CFLRATIO,TOLER,OMEGA,DTSMALL
      INTEGER INMAX(MAXNOFVAR,2),INDEL(MAXNOFVAR,2),ITER,NITER,NSUBIT,
     &ITMAX,IVCNVG,ISTMP,ISTART,IBAK,CFL_RAMP
C
      COMMON/CONV_R/RESMAX,RESL2,DELMAX,DELL2,RESL20,RESMAX0,CFL,
     &CFLMAX,CFLRATIO,TOLER,OMEGA,DTSMALL
      COMMON/CONV_I/INMAX,INDEL,ITER,NITER,NSUBIT,ITMAX,
     &IVCNVG,ISTMP,ISTART,IBAK,CFL_RAMP
C
C       Second column is for the turbulence model
C
C
C	RESMAX(1:NOFVAR)stores the L_infty norm of the residual
C	RESL2(1:NOFVAR) stores the L_2 norm of the residual for each
C                       variable
C	DELMAX(1:NOFVAR)stores the L_infty norm of the update
C	DELL2(1:NOFVAR) stores the L_2 norm of the update for each
C                       variable
C	RESL20 	stores the L_2(IVCNVG) norm of the residual after 
C               the 1st iter.
C	RESMAX0 stores the L_infty(IVCNVG) norm of the residual 
C               after the 1st iter.
C       CFL     is the CFL number
C       CFLMAX  is the maximum timestep allowed for implicit timestepping
C       CFL_RAMP  is the timestep ramping strategy
C               look in implicit.h for possible values
C       TOLER   is the convergence threshold, i.e. 
C               the iteration is halted when RESL2(IVCNVG) <= TOLER
C	INMAX,2(1:NOFVAR) stores the node where RESMAX is found
C	INDEL(1:NOFVAR) stores the node where DELMAX is found
C       ITER    is the iteration counter for the current run
C       NITER   is the "global" iteration counter (not implemented, though)
C       NSUBIT  is the number of sub iteration on the turbulence eqn.
C       ITMAX   maximum number of non linear iterations allowed
C               i.e. ITER=1,2,...,ITMAX
C       IVCNVG  is the variable on which the convergence test is based 
C               (see TOLER above)
C       ISTMP   solution is printed to NOUT each ISTMP ITERations
C       ISTART  restarts from a previous solution if ISTART <> 0
C       IBAK    solution is saved to disk each IBAK ITERations 
