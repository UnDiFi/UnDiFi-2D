C
C    $Id: iset.com,v 1.8 2013/07/17 10:35:37 abonfi Exp $
C
C    These are the index sets in PETSc style that address "nodal" boundary conditions
C
      COMMON/COMISET/
     1SupersonicNodes     , SupersonicVariables , NoSlipNodes         ,
     2NoSlipVelocities    , FreestreamTurbulence, Isothermal          ,
     4MotionSolverBCS     , HangingNodes        , Dirichlet4Poisson   ,
     5bndrynodes
C
C     SupersonicNodes
C     SupersonicVariables
C     NoSlipNodes
C     FreestreamTurbulence
C     Isothermal
C     MotionSolverBCS  
C     Dirichlet4Poisson
C
