C///////////////////////////////////////////////////////////////////////
C     $Id: turb.com,v 1.5 2013/01/25 08:15:43 abonfi Exp $
C///////////////////////////////////////////////////////////////////////

      DOUBLE PRECISION TCB1,TCB2,TPR,TPR1,TK,TCW1,TCW2,TCW3,TCV1,TCV2,
     &        TCT1,TCT2,TCT3,TCT4,TDXT,TST
      INTEGER TTFLAG,TTELEM
C
C**************************************************************
C
C Set Turbulence model parameters
C
C
C*************************************************************
      PARAMETER( TCB1 = 0.1355d0, TCB2 = 0.622d0, TPR = 0.9d0,
     1TPR1 = 2.d0/3.d0, TK = 0.41d0, 
     2TCW1 = TCB1/TK**2 + (1.0+TCB2)/TPR1, TCW2 = 0.3d0,
     3TCW3 = 2.d0, TCV1 = 7.1d0, TCV2 = 5.d0, TCT1 = 10.,
     4TCT2 = 2.0, TCT3 = 1.2, TCT4 = 0.5)

      COMMON  / TURBPAR / TDXT,TST,TTFLAG,TTELEM

C     Variable name            Comment
C     TCB1                     Calibration constants (Spalart & Allmaras model)
C     TCB2                     Calibration constants (Spalart & Allmaras model)
C     TCW1                     Calibration constants (Spalart & Allmaras model)
C     TCW2                     Calibration constants (Spalart & Allmaras model)
C     TCW3                     Calibration constants (Spalart & Allmaras model)
C     TCV1                     Calibration constants (Spalart & Allmaras model)
C     TCT1                     Calibration constants (Spalart & Allmaras model)
C     TCT2                     Calibration constants (Spalart & Allmaras model)
C     TCT3                     Calibration constants (Spalart & Allmaras model)
C     TCT4                     Calibration constants (Spalart & Allmaras model)
C     TPR                      Turbulent Prandtl number (conducibility comp.)
C     TPR1                     Turbulent Prandtl number (diffusion term comp.)
C     TK                       Von Karman constant
C     TTFLAG                   Flag for trip term activation
C     TELEM                    cell used to compute the vorticity TST
C
C
C     Variable name            Comment
C     TTD                      TRIP POINT DISTANCE
C     TST                      VORTICITY MODULE IN TRIP POINT
C     TDXT                     Avg. mesh spacing at trip point
