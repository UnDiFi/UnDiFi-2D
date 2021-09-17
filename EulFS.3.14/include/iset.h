C
C    $Id: iset.h,v 1.1 2020/04/23 09:27:12 abonfi Exp $
C
C    These are integers used to address PETSc index sets stored in an array of index-sets
C
      integer 
     1SupersonicNodes , SupersonicVariables , NoSlipNodes      ,
     2NoSlipVelocities, FreestreamTurbulence, Isothermal       ,
     4MotionSolverBCS , HangingNodes        , Dirichlet4Poisson
c
c     starting with Petsc 3.8.* the index sets are stored in an array of derived types 
c     the first NCOLOR+1 locations are used to store the addresses of the gridpoints
c     located on the boundary patches from 0 to NCOLOR
c
      parameter(SupersonicNodes     =NCOLOR+1,
     &	        SupersonicVariables =NCOLOR+2,
     &	        NoSlipNodes         =NCOLOR+3,
     &          NoSlipVelocities    =NCOLOR+4,
     &          FreestreamTurbulence=NCOLOR+5,
     &          Isothermal          =NCOLOR+6,
     &          MotionSolverBCS     =NCOLOR+7,
     &          HangingNodes        =NCOLOR+8,
     &          Dirichlet4Poisson   =NCOLOR+9)
C
