      SUBROUTINE SETIBC(NDIM,NOFVAR)
C
C     This routine SETs Initial and Boundary Conditions
C
C     $Id: setibc.F,v 1.32 2020/04/23 09:37:34 abonfi Exp $
C
C
#include "petsc/finclude/petscsys.h"
      use petscsys
C
      IMPLICIT NONE
C
C
      PetscBool IFLAG
C
      INCLUDE 'paramt.h'
      INCLUDE 'constants.h'
      INCLUDE 'plasma.h'
      include 'bnd.h'
      INCLUDE 'implicit.h'
      INCLUDE 'visco.h'
      INCLUDE 'chem.h'
      INCLUDE 'paramchem.h'
C
      INCLUDE 'bnd.com'
      INCLUDE 'conv.com'
      INCLUDE 'dofs.com'
      INCLUDE 'flags.com'
      INCLUDE 'io.com'
      INCLUDE 'merkle.com'
      INCLUDE 'periodic.com'
      INCLUDE 'pfcgas.com'
      INCLUDE 'stream.com'
      INCLUDE 'suther.com'
      INCLUDE 'verbose.com'
      INCLUDE 'visco.com'
      INCLUDE 'streamplasma.com'
      INCLUDE 'ioplasma.com' 
      INCLUDE 'commonchem.inc' 
      INCLUDE 'electric.com' 
C
      INTEGER NDIM,NOFVAR
      DOUBLE PRECISION SDIM,VISC0,T0DIM,R
      PARAMETER(SDIM=110.5d0,VISC0=0.1716d-4,T0DIM=273.1d0,R=SDIM/T0DIM)
      INTEGER INTERNAL_FLOW,EXTERNAL_FLOW
      PARAMETER(INTERNAL_FLOW=0,EXTERNAL_FLOW=1)
C
      INTEGER NONDIM
C
      INTEGER MY_PE
      COMMON/MPICOM/MY_PE
C
      CHARACTER*80 STRINGA
      CHARACTER*72 ERRMSG
      INTEGER IERR,NVALS,ISP,IR,IFIT,IMDA 
      INTEGER NERR,IOPT
      DOUBLE PRECISION RTOT,PTOT,TTOT,TEMP,HREF,T
      DOUBLE PRECISION A1,R1,T1,P1,H1,U1,V1,W1,VISCL,Z1,MACH1,RF
      DOUBLE PRECISION IONDEGREE,RHOS(NSP),ALPHAREF(NSP),RMIX,HFTOT,
     +                 CONC_M,SOURCE(NSP),DAS(NSP),HELP,
     +                 DAMAX,RGASS(NSP),VGEN
      LOGICAL COMPRESSIBLE,INCOMPRESSIBLE,VISCOUS,LFLAG,PLASMA
C
      DOUBLE PRECISION PRESSC,PRESS4AR
 
      EXTERNAL PRESSC,PRESS4AR
C
C     Sutherland's law in dimensional form
C
      VISCL(T) = VISC0 * ((ONE+R)/(T/T0DIM+R))* ((T/T0DIM)**1.5d0)
C
C     Executable Statements ..
C
      COMPRESSIBLE = (ABS(KAN).EQ.4)
      INCOMPRESSIBLE = (.NOT.COMPRESSIBLE)
      VISCOUS = (KAN.EQ.2).OR.(KAN.EQ.4)
      PLASMA = (KAN.EQ.-3)
C
C
      IF( PERIODIC_MESH .AND. ANNULAR )THEN
         CALL INIQMAT(QMAT,NOFVAR)
!        CALL R8Mat_Print('General',' ',Nofvar,Nofvar,QMAT,
!    +   Nofvar,'Rotation Matrix inside setibc ',IERR)
!        pause
      ENDIF
C
C        
C
      IF(MY_PE.EQ.0)WRITE (IWUNIT,FMT=2000) 
 2000 FORMAT (/,/,' INITIAL-BOUNDARY CONDITIONS ',/,1X,28('='),/)
C
C
C     *********************************************
C     Read number of chemical reactions
C     *********************************************  
C         
      CALL PetscOptionsGetInt(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,
     +      '-chem',CHEM,IFLAG,IERR)
      IF( IFLAG .EQV. PETSC_FALSE )THEN
                CHEM = 0
      ENDIF
      WRITE(6,*)'Chemical source ','CHEM ',CHEM
!
      CALL PetscOptionsGetInt(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,
     +      '-reaction_number',NREAC,IFLAG,IERR)
      IF( IFLAG .EQV. PETSC_FALSE )THEN
                NREAC = 3
      ENDIF
      WRITE(6,*)'Chemical reactions number ','NREAC  ',NREAC
C
      IF( LAPLACE )THEN
C
C     *********************************************
C     Read electromagnetic quantities
C     *********************************************  
C         
            CALL PetscOptionsGetInt(PETSC_NULL_OPTIONS,
     +           PETSC_NULL_CHARACTER,'-ohm_source',OHM,IFLAG,IERR)
            IF( IFLAG .EQV. PETSC_FALSE )THEN
                OHM = 0
            ENDIF
           WRITE(6,*)' ','Ohmic heat source ',OHM

C
C         Read generator electrical potential 
C
caldo
              CALL PetscOptionsGetReal(PETSC_NULL_OPTIONS,
     +             PETSC_NULL_CHARACTER,'-generator_potential',VGEN,
     &             IFLAG,IERR)
              IF( IFLAG .EQV. PETSC_FALSE )THEN                 
                  VGEN = ONE
              ENDIF
              IF(MY_PE.EQ.0)
     &        WRITE(IWUNIT,320)'Generator Potential','VGEN',VGEN
              PHIREF = VGEN
caldo
C
C         Read circuit resistence 
C
              CALL PetscOptionsGetReal(PETSC_NULL_OPTIONS,
     +             PETSC_NULL_CHARACTER,'-circuit_resistence',RC,IFLAG,
     &             IERR)
              IF( IFLAG .EQV. PETSC_FALSE )THEN
                  RC = 0.d0
              ENDIF
              IF(MY_PE.EQ.0)
     &        WRITE(IWUNIT,320)'Circuit Resistence','RC',RC
!
      ENDIF ! test on Laplace
C
C     *********************************************
C     Check non-dimensionalisation of the equations
C     *********************************************  
C
      CALL PetscOptionsGetString(PETSC_NULL_OPTIONS,
     +PETSC_NULL_CHARACTER,'-nondimensionalisation',STRINGA,IFLAG,IERR)
      IF( IFLAG .EQV. PETSC_FALSE )THEN
          WRITE(ERRMSG(1:56),300)'-nondimensionalisation'
          NERR = 4
          IOPT = 1
          CALL SETERR(ERRMSG,56,NERR,IOPT)
      ELSE
          IF(STRINGA(1:8).EQ.'external') THEN
              NONDIM = EXTERNAL_FLOW
          ELSEIF(STRINGA(1:8).EQ.'internal') THEN
              NONDIM = INTERNAL_FLOW
          ELSE
              WRITE(ERRMSG,310)STRINGA,'-nondimensionalisation'
              NERR = 4
              IOPT = 1
              CALL SETERR(ERRMSG,72,NERR,IOPT)
          ENDIF
      ENDIF
C
C         read adiabatic index and gas constant
C
              CALL PetscOptionsGetReal(PETSC_NULL_OPTIONS,
     +             PETSC_NULL_CHARACTER,'-adiabatic_index',GAM,IFLAG,
     &             IERR)
C
C         If not explicitely set, we set the adiabatic index equal to 1.4
C
              IF( IFLAG .EQV. PETSC_FALSE )THEN
                  GAM = 1.4d0 ! diatomic gas
              ENDIF
              GM1 = GAM-ONE
              GM1OG = GM1/GAM
              GP1OG = (GAM+ONE)/GAM
              GOGM1 = GAM/GM1
              GINV = ONE/GAM
C
              CALL PetscOptionsGetReal(PETSC_NULL_OPTIONS,
     +             PETSC_NULL_CHARACTER,'-specific_gas_constant',RGAS,
     &             IFLAG,IERR)
C
C         If not explicitely set, we set the specific gas constant equal to that of air = 287
C
              IF( IFLAG .EQV. PETSC_FALSE )THEN
                  RGAS = 287.d0 ! specific gas constant for air
              ENDIF
C
C
      IF(MY_PE.EQ.0)THEN
         WRITE(IWUNIT,320)'Specific heats ratio ','GAMMA',GAM
         WRITE(IWUNIT,320)'Specific gas constant ','RGAS',RGAS
      ENDIF
C
      IF( NONDIM .EQ. EXTERNAL_FLOW )THEN
C
          IF(COMPRESSIBLE.OR.PLASMA)THEN
C
C         read freestream Mach number
C
              CALL PetscOptionsGetReal(PETSC_NULL_OPTIONS,
     +PETSC_NULL_CHARACTER,'-freestream_Mach_number',M_INFTY,IFLAG,IERR)
              IF( IFLAG .EQV. PETSC_FALSE )THEN
                  WRITE(ERRMSG(1:56),300)'-freestream_Mach_number'
                  NERR = 4
                  IOPT = 1
                  CALL SETERR(ERRMSG,56,NERR,IOPT)
              ENDIF
              IF(MY_PE.EQ.0)
     &        WRITE(IWUNIT,320)'Freestream Mach number','M_INFTY',
     +M_INFTY
C
C         set constant for Merkle's preconditioner and isothermal bcs
C
              RSTAR = ONE/(GAM*M_INFTY*M_INFTY)
          ENDIF
C
          IF( (COMPRESSIBLE.AND.VISCOUS) .OR. PLASMA )THEN
C
C         read Freestream Static Temperature
C
              CALL PetscOptionsGetReal(PETSC_NULL_OPTIONS,
     +             PETSC_NULL_CHARACTER,'-freestream_temperature',TREF,
     &             IFLAG,IERR)
              IF( IFLAG .EQV. PETSC_FALSE )THEN
                  TREF = 273.d0
              ENDIF
              IF(MY_PE.EQ.0)
     &        WRITE(IWUNIT,320)'Freestream Temperature','TREF',TREF
          ENDIF
C
         IF(PLASMA) THEN
C
!          CALL PetscOptionsGetInt(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,
!     +        '-reactions_number',NR,IFLAG,IERR)
!              IF( IFLAG .EQV. PETSC_FALSE )THEN
!                  NR = 3
!              ENDIF

C         read Freestream Static Pressure
C
              CALL PetscOptionsGetReal(PETSC_NULL_OPTIONS,
     +             PETSC_NULL_CHARACTER,'-freestream_pressure',PREF,
     &             IFLAG,IERR)
              IF( IFLAG .EQV. PETSC_FALSE )THEN
                  PREF = 1.d5
              ENDIF
              IF(MY_PE.EQ.0)
     &        WRITE(IWUNIT,320)'Freestream Pressure','PREF',PREF
C
C         read Freestream Ionization degree
C
              CALL PetscOptionsGetReal(PETSC_NULL_OPTIONS,
     +             PETSC_NULL_CHARACTER,'-ionization_degree',IONDEGREE,
     &             IFLAG,IERR)
              IF( IFLAG .EQV. PETSC_FALSE )THEN
C                  CALL IONDEGREE4Ar(PREF,TREF,IONDEGREE)
C                  IONDEGREE = 0.d0
              ENDIF
              IF(MY_PE.EQ.0)
     &        WRITE(IWUNIT,320)'Ionization Degree','IONDEGREE',IONDEGREE
C
C         read initial molar concentration of metastable Argon     
C
              CALL PetscOptionsGetReal(PETSC_NULL_OPTIONS,
     +             PETSC_NULL_CHARACTER,'-conc_metastable',CONC_M,
     &             IFLAG,IERR)
              IF( IFLAG .EQV. PETSC_FALSE )THEN
!                  CONC_M = 0.d0
                   CALL BOLTZEQ4Ar(PREF,TREF,IONDEGREE,CONC_M)
              ENDIF
C              
              IF( CONC_M.GT.(1.0d0-2.0d0*IONDEGREE) ) THEN
                 CALL BOLTZEQ4Ar(PREF,TREF,IONDEGREE,CONC_M)
              ENDIF
C
              IF(MY_PE.EQ.0)
     &        WRITE(IWUNIT,320)'Ar* concentration','CONC_M',CONC_M
C                      
C          compute Freestream Density  
C
             CALL PLASFLOW_INF(PREF,TREF,IONDEGREE,RHOS,RREF,RMIX,
     &       RGASS,HFTOT,CONC_M)
C              
C          compute Freestream Speed
C         
              UREF = M_INFTY * SQRT(GAM*RMIX*TREF)
C
C          compute Freestream Enthalpy           
C
              HREF = UREF**2
C
              WRITE(IWUNIT,320)'Freestream Speed','UREF',UREF
C
              DO ISP=1,NSP
                 ALPHAREF(ISP) = RHOS(ISP) / RREF
              ENDDO         

C         compute Freestream chemical source term
C
!             CALL CHEMSOURCEII(TREF,PREF,RHOS,SOURCE)              
C
C         read Reference Length
C
              CALL PetscOptionsGetReal(PETSC_NULL_OPTIONS,
     +             PETSC_NULL_CHARACTER,'-reference_length',LREF,IFLAG,
     &             IERR)
              IF( IFLAG .EQV. PETSC_FALSE )THEN
                  LREF = 1.d0
              ENDIF
              IF(MY_PE.EQ.0)
     &        WRITE(IWUNIT,320)'Reference Length','LREF',LREF
C
C         compute Damkohler numbers 
C
!           DO ISP = 1,NSP
!              DAS(ISP) = SOURCE(ISP)*LREF/(RREF*UREF)
!           ENDDO
C
C         compute reference Damkohler number
C
!           DAMAX=0.0d0
!           DO ISP = 1,NSP
!           WRITE(IWUNIT,320)'Damkohler number','DAS',DAS(ISP)
!           WRITE(IWUNIT,320)'Source term','SOURCE',SOURCE(ISP)
!              IF(ABS(DAS(ISP)).GT.DAMAX)THEN
!                 DAMAX=ABS(DAS(ISP))
!                 IMDA=ISP
!              ENDIF                     
!           ENDDO                       
C
!           WRITE(IWUNIT,320)'Max Damkohler number','DAMAX',DAMAX
!           WRITE(IWUNIT,320)'Species Index Da Max','ISP',IMDA       
C          
C         common in streamplasma.com
C          
          PREFP = RREF*UREF**2
          RREFP = RREF
          UREFP = UREF
          HREFP = HREF
          TREFP = TREF
          LREFP = LREF
          DO ISP=1,NSP
             RSTARP(ISP)=RGASS(ISP)*TREF/UREF**2
          ENDDO
C
          WRITE(IWUNIT,320)'Freestream Pressure plasma','PREFP',PREFP  
C        
          ENDIF
C
C         read flow angles
C
          NVALS = NDIM
          CALL PetscOptionsGetRealArray(PETSC_NULL_OPTIONS,
     +         PETSC_NULL_CHARACTER,'-flow_angles',FLOWDIR,NVALS,IFLAG,
     &         IERR)
          IF( IFLAG .EQV. PETSC_FALSE )THEN
              WRITE(ERRMSG(1:56),300)'-flow_angles'
              NERR = 4
              IOPT = 1
              CALL SETERR(ERRMSG,56,NERR,IOPT)
          ENDIF
C
          IF(VISCOUS)THEN
C
C    read Reynolds number
C
              CALL PetscOptionsGetReal(PETSC_NULL_OPTIONS,
     +             PETSC_NULL_CHARACTER,'-Reynolds',RE,IFLAG,IERR)
              IF( IFLAG .EQV. PETSC_FALSE )THEN
                  WRITE(ERRMSG,FMT=300)'-Reynolds'
                  NERR = 4
                  IOPT = 1
                  CALL SETERR(ERRMSG,72,NERR,IOPT)
              ENDIF
              IF(MY_PE.EQ.0)WRITE(IWUNIT,320)'Reynolds Number ','RE',RE
              REINV=1.d0/RE
          ELSE
              RE = 1.E+38
              REINV = ZERO
              PRANDTL = 1.E+38
          ENDIF
C
C     Set freestream values (to be used for far-field
C         boundary conditions)
C
          IF(COMPRESSIBLE)THEN
C
C     We take freestream density and velocity equal to 1
C     and compute enthalpy accordingly (based on Mach number)
C
C
              HREF = (TWO+GM1*M_INFTY*M_INFTY)/(TWO*GM1*M_INFTY*M_INFTY)
              U_INFTY(1) = ONE
              U_INFTY(2) = HREF
              U_INFTY(3) = FLOWDIR(1)
              U_INFTY(4) = FLOWDIR(2)
              U_INFTY(5) = FLOWDIR(3)
          ELSEIF(PLASMA)THEN
              HREF = (TWO+GM1*M_INFTY*M_INFTY)/(TWO*GM1*M_INFTY*M_INFTY)
              HREF = HREF + HFTOT/(UREF**2)
              U_INFTY(1) = ALPHAREF(1)
              U_INFTY(2) = ALPHAREF(2)
              U_INFTY(3) = ALPHAREF(3)
              U_INFTY(4) = ALPHAREF(4)
              U_INFTY(5) = HREF
              U_INFTY(6) = FLOWDIR(1)
              U_INFTY(7) = FLOWDIR(2)
              U_INFTY(8) = FLOWDIR(3)
          ELSE
              U_INFTY(1) = ONE
              U_INFTY(2) = FLOWDIR(1)
              U_INFTY(3) = FLOWDIR(2)
              U_INFTY(4) = FLOWDIR(3)
          ENDIF
C
C     outlet static pressure (NON dimensional)
C     this is only required if any of the boundaries
C     is of type BC_CNST_PRESSURE
C
          CALL PetscOptionsGetReal(PETSC_NULL_OPTIONS,
     +         PETSC_NULL_CHARACTER,'-outlet_static_pressure',POUTLET,
     &         IFLAG,IERR)
          IF( IFLAG .EQV. PETSC_FALSE )THEN
C
C     If one does not specify the value of the outlet
C     static pressure, this is
C     computed using the freestream conditions
C
              IF(COMPRESSIBLE)THEN
                  POUTLET = PRESSC(NDIM,U_INFTY(1))
              ELSEIF(PLASMA)THEN
                  POUTLET = PRESS4Ar( NDIM,U_INFTY(1))                  
              ELSE
                  POUTLET = ONE
              ENDIF
          ENDIF
C
C
      ELSEIF( NONDIM .EQ. INTERNAL_FLOW )THEN
C
C         set constant for Merkle's preconditioner and isothermal bcs
C
          RSTAR = ONE
C
          CALL PetscOptionsGetReal(PETSC_NULL_OPTIONS,
     +         PETSC_NULL_CHARACTER,'-inlet_total_temperature',TTOT,
     &         IFLAG,IERR)
          IF( IFLAG .EQV. PETSC_FALSE )THEN
              WRITE(ERRMSG(1:56),300)'-inlet_total_temperature'
              NERR = 4
              IOPT = 1
              CALL SETERR(ERRMSG,56,NERR,IOPT)
          ENDIF
C
          TREF=TTOT
C
          IF(MY_PE.EQ.0)
     &    WRITE(IWUNIT,320)'Inlet Total Temperature','TREF',TREF
C
          CALL PetscOptionsGetReal(PETSC_NULL_OPTIONS,
     +         PETSC_NULL_CHARACTER,'-inlet_total_pressure',PTOT,
     &         IFLAG,IERR)
          IF( IFLAG .EQV. PETSC_FALSE )THEN
              WRITE(ERRMSG(1:56),300)'-inlet_total_pressure'
              NERR = 4
              IOPT = 1
              CALL SETERR(ERRMSG,56,NERR,IOPT)
          ENDIF
          IF(MY_PE.EQ.0)
     &    WRITE(IWUNIT,320)'Inlet Total Pressure','PTOT',PTOT
C
          PREF=PTOT
C
          CALL PetscOptionsGetReal(PETSC_NULL_OPTIONS,
     +         PETSC_NULL_CHARACTER,'-outlet_static_pressure',POUTLET,
     &         IFLAG,IERR)
 
          IF( IFLAG .EQV. PETSC_FALSE )THEN
              WRITE(ERRMSG(1:56),300)'-outlet_static_pressure'
              NERR = 4
              IOPT = 1
              CALL SETERR(ERRMSG,56,NERR,IOPT)
          ENDIF
          IF(MY_PE.EQ.0)
     &    WRITE(IWUNIT,320)'Outlet Static Pressure','POUTLET',POUTLET
          POUTLET=POUTLET/PTOT
C
          CALL PetscOptionsGetReal(PETSC_NULL_OPTIONS,
     +         PETSC_NULL_CHARACTER,'-reference_length',LREF,IFLAG,IERR)
          IF( IFLAG .EQV. PETSC_FALSE )THEN
              WRITE(ERRMSG(1:56),300)'-reference_length'
              NERR = 4
              IOPT = 1
              CALL SETERR(ERRMSG,56,NERR,IOPT)
          ENDIF
          IF(MY_PE.EQ.0)WRITE(IWUNIT,320)'Reference Length ','LREF',
     +  LREF
C
C         read flow angles
C
          NVALS = NDIM
          CALL PetscOptionsGetRealArray(PETSC_NULL_OPTIONS,
     +         PETSC_NULL_CHARACTER,'-flow_angles',FLOWDIR,NVALS,IFLAG,
     &         IERR)
          IF( IFLAG .EQV. PETSC_FALSE )THEN
              WRITE(ERRMSG(1:56),300)'-flow_angles'
              NERR = 4
              IOPT = 1
              CALL SETERR(ERRMSG,56,NERR,IOPT)
          ENDIF
C
        IF ( PLASMA )THEN
C       Setting of initial and bc for internal adimensionalization and plasma flow           
C            STOP 'Internal flows NOT yet available for plasma'
C
C         read inlet Ionization degree
C
              CALL PetscOptionsGetReal(PETSC_NULL_OPTIONS,
     +             PETSC_NULL_CHARACTER,'-ionization_degree',IONDEGREE,
     &             IFLAG,IERR)
              IF( IFLAG .EQV. PETSC_FALSE )THEN
                  WRITE(ERRMSG(1:56),300)'-ionization_degree'
                  CALL IONDEGREE4Ar(PREF,TREF,IONDEGREE)
!                  IONDEGREE = 0.d0
                  write(6,*)'iondegree not declared'
                  pause
              ENDIF
!              CALL IONDEGREE4Ar(PREF,TREF,IONDEGREE)
              IF(MY_PE.EQ.0)                               
     &        WRITE(IWUNIT,320)'Ionization Degree','IONDEGREE',IONDEGREE
C
C         read inlet molar concentration of metastable Argon     
C
              CALL PetscOptionsGetReal(PETSC_NULL_OPTIONS,
     +             PETSC_NULL_CHARACTER,'-conc_metastable',CONC_M,
     &             IFLAG,IERR)
              IF( IFLAG .EQV. PETSC_FALSE )THEN
                  WRITE(ERRMSG(1:56),300)'-conc_metastable'
!                  CALL BOLTZEQ4Ar(PREF,TREF,IONDEGREE,CONC_M)
                  write(6,*)'Conc_M not declared'
                  pause
              ENDIF
C              
              IF( CONC_M.GT.(1.0d0-2.0d0*IONDEGREE) ) THEN
!                  CALL BOLTZEQ4Ar(PREF,TREF,IONDEGREE,CONC_M)
              ENDIF
!              CALL BOLTZEQ4Ar(PREF,TREF,IONDEGREE,CONC_M)
C
              IF(MY_PE.EQ.0)              
     &        WRITE(IWUNIT,320)'Ar* concentration','CONC_M',CONC_M
C                      
C          compute inlet density  
C
             CALL PLASFLOW_INF(PTOT,TTOT,IONDEGREE,RHOS,RTOT,RMIX,
     &       RGASS,HFTOT,CONC_M)
C
C          compute inlet adimensional gas constant (R_i1^*=\rho^0 T^0/p^0 *Ri1 )
C
             DO ISP = 1,NSP
                 RSSTAR(ISP) = RTOT*TTOT/PTOT * RGASS(ISP)
             ENDDO
C
C          compute species concentrations (\rho_i1^0/rho_1^0) 
C
              DO ISP=1,NSP
                 ALPHAREF(ISP) = RHOS(ISP) / RTOT
              ENDDO              
C
C          compute dimensional reference values 
C
             UREF = SQRT(RMIX*TTOT)
             RREF = RTOT
             HREF = RMIX*TTOT
C
C         common in streamplasma.com
C          
            PREFP = PREF
            RREFP = RREF
            UREFP = UREF
            HREFP = HREF
            TREFP = TREF
            LREFP = LREF   
!            RSTARP = RGASS
!            RSTARP = RSSTAR
            RMIXSTAR = 0.d0
            DO ISP = 1,NSP
               RMIXSTAR = RMIXSTAR + ALPHAREF(ISP)*RSSTAR(ISP)
               RSTARP(ISP) = RSSTAR(ISP)
            ENDDO
           WRITE(IWUNIT,320)'Reference Density','RREF',RREF
!           WRITE(IWUNIT,320)'Dimensionles Mixture Constant ',
!     +                      'RMIXSTAR', RMIXSTAR
C           WRITE(IWUNIT,320)'Dimensional Mixture Constant','RMIX',RMIX
C           WRITE(IWUNIT,320)'Reference Temperature','TREF',TREF
C           WRITE(IWUNIT,320)'Reference Pressure','PREF',PREF
C
C         common in ioplasma.com
C
            DO ISP = 1,NSP
               ALPHA1(ISP) = ALPHAREF(ISP)
            ENDDO      
      
C     the inlet mach number is only required to initialise
C     the flowfield
C
          CALL PetscOptionsGetReal(PETSC_NULL_OPTIONS,
     +    PETSC_NULL_CHARACTER,'-inlet_mach_number',MACH1,IFLAG,IERR)
          IF( IFLAG .EQV. PETSC_FALSE )THEN
              WRITE(ERRMSG(1:56),300)'-inlet_mach_number'
              NERR = 4
              IOPT = 1
              CALL SETERR(ERRMSG,56,NERR,IOPT)
          ENDIF
C
C     what follows is only required to compute
C     the initial conditions 
C
          TEMP = ONE+HALF*GM1*MACH1*MACH1
          t1=ttot/temp
          p1=ptot/(temp**GOGM1)
          r1=p1/(RMIX*t1)
          CALL PLASFLOW_INF(p1,t1,IONDEGREE,RHOS,r1,RMIX,
     &                  RGASS,HFTOT,CONC_M)
         
!          Do isp=1,NSP
!              write(6,*)'alpha',isp,RHOS(isp)/r1
!              write(6,*)'alpharef',isp,ALPHAREF(isp)
!          enddo  
!          pause  
          
          a1=sqrt(GAM*RMIX*t1)
          u1=mach1*a1*flowdir(1)
          v1=mach1*a1*flowdir(2)
          w1=mach1*a1*flowdir(3)
          temp=HALF/GM1*(TWO+GM1*mach1*mach1)
          h1 = a1*a1*temp + HFTOT/HREF            
          r1=r1/rtot
          h1=h1/href
          u1=u1/uref
          v1=v1/uref
          w1=w1/uref 
          z1=sqrt(r1)
C
          DO ISP = 1,NSP 
              U_infty(ISP) = z1*ALPHAREF(ISP)
          ENDDO
          U_infty(NSP+1) = z1 * h1
          U_infty(NSP+2) = z1 * u1
          U_infty(NSP+3) = z1 * v1
          U_infty(NSP+4) = z1 * w1          
C
!       Help = LREF/RREF/UREF 
!       write(6,*)'L/(rho*u)=',Help

C         compute Freestream chemical source term
!
!              CALL CHEMSOURCEII(TREF,PREF,RHOS,SOURCE)   
C
C         compute Damkohler numbers 
C
!           DO ISP = 1,NSP
!              DAS(ISP) = SOURCE(ISP)*LREF/(RREF*UREF)
!           ENDDO
C
C         compute reference Damkohler number
C
!           DAMAX=0.0d0
!           DO ISP = 1,NSP
!           WRITE(IWUNIT,320)'Damkohler number','DAS',DAS(ISP)
!           WRITE(IWUNIT,320)'Source term','SOURCE',SOURCE(ISP)
!              IF(ABS(DAS(ISP)).GT.DAMAX)THEN
!                 DAMAX=ABS(DAS(ISP))
!                 IMDA=ISP
!              ENDIF
!           ENDDO
C
!           WRITE(IWUNIT,320)'Max Damkohler number','DAMAX',DAMAX
!           WRITE(6,*)'Species Index Da Max','ISP',IMDA

           WRITE(IWUNIT,320)'Reference Speed','UREF',UREF      
C
C
       ELSEIF ( COMPRESSIBLE )THEN
C
C     dimensional reference values: to be used to compute the
C     Reynolds number
C
          UREF = SQRT(RGAS*TTOT)
          RTOT = PTOT/(RGAS*TTOT)
          RREF = RTOT
          HREF = RGAS*TTOT
          RE = RTOT*UREF*LREF/VISCL(TTOT)
          REINV = ONE/RE
          IF(MY_PE.EQ.0)WRITE(IWUNIT,320)'Reynolds Number ','RE',RE
C
C     the inlet mach number is only required to initialise
C     the flowfield
C
          CALL PetscOptionsGetReal(PETSC_NULL_OPTIONS,
     +    PETSC_NULL_CHARACTER,'-inlet_mach_number',MACH1,IFLAG,IERR)
          IF( IFLAG .EQV. PETSC_FALSE )THEN
              WRITE(ERRMSG(1:56),300)'-inlet_mach_number'
              NERR = 4
              IOPT = 1
              CALL SETERR(ERRMSG,56,NERR,IOPT)
          ENDIF
C
C     what follows is only required to compute
C     the initial conditions
C
          TEMP = ONE+HALF*GM1*MACH1*MACH1
          t1=ttot/temp
          p1=ptot/(temp**GOGM1)
          r1=p1/(rgas*t1)
          a1=sqrt(GAM*rgas*t1)
          u1=mach1*a1*flowdir(1)
          v1=mach1*a1*flowdir(2)
          w1=mach1*a1*flowdir(3)
          temp=HALF/GM1*(TWO+GM1*mach1*mach1)
          h1 = a1*a1*temp
C
          r1=r1/rtot
          h1=h1/href
          u1=u1/uref
          v1=v1/uref
          w1=w1/uref
          z1=sqrt(r1)
C
          U_infty(1) = z1
          U_infty(2) = z1 * h1
          U_infty(3) = z1 * u1
          U_infty(4) = z1 * v1
          U_infty(5) = z1 * w1
C
C        write(6,*)'sono in compressible internal'
C
          ELSEIF( INCOMPRESSIBLE )THEN
C
C     dimensional reference values: to be used to compute the
C     Reynolds number
C
C     set constant density 
C
          CALL PetscOptionsGetReal(PETSC_NULL_OPTIONS,
     +    PETSC_NULL_CHARACTER,'-constant_density',RTOT,IFLAG,IERR)
          IF( IFLAG .EQV. PETSC_FALSE )THEN
              WRITE(ERRMSG(1:56),300)'-constant_density'
              NERR = 4
              IOPT = 1
              CALL SETERR(ERRMSG,56,NERR,IOPT)
          ENDIF
          IF(MY_PE.EQ.0)WRITE(IWUNIT,320)'Constant density',
     +    'kg/m**3',RTOT
          RREF=RTOT 
C
          IF(MY_PE.EQ.0)WRITE(IWUNIT,320)'Total temperature',
     +    'K',TTOT
          TREF=TTOT
C
C         reference velocity
C
          UREF = SQRT(PTOT/RTOT)
          IF(MY_PE.EQ.0)WRITE(IWUNIT,320)'Reference velocity',
     +    'm/s',UREF
          RE = RTOT*UREF*LREF/VISCL(TTOT)
          REINV = ONE/RE
          IF(MY_PE.EQ.0)THEN
             WRITE(IWUNIT,320)'Viscosity ','kg/(m*s)',VISCL(TTOT)
             WRITE(IWUNIT,320)'Kinematic Viscosity ',
     +                        'm**2/s',VISCL(TTOT)/RTOT
             WRITE(IWUNIT,320)'Reynolds Number ','RE',RE
          ENDIF
C
          a1=sqrt(GAM*rgas*TTOT)
          u1=0.1d0*a1*flowdir(1)
          v1=0.1d0*a1*flowdir(2)
          w1=0.1d0*a1*flowdir(3)
          temp = half * ( u1*u1 + v1*v1 + w1*w1 )
          p1 = PTOT - temp
          U_infty(1) = p1/PTOT
          U_infty(2) = u1/UREF
          U_infty(3) = v1/UREF
          U_infty(4) = w1/UREF
C
C        write(6,*) 'sono in incopressible internal'
C
          ENDIF
C
      ENDIF ! kind of non-dimensionalization
      IF(MY_PE.EQ.0)
     &WRITE(IWUNIT,320)"Dimensionless gas constant ",'RSTAR',RSTAR
C
C     A few options follow that are only meaningful for
C     compressible flows
C
      IF(COMPRESSIBLE.AND.VISCOUS)THEN
          CALL PetscOptionsGetReal(PETSC_NULL_OPTIONS,
     +         PETSC_NULL_CHARACTER,'-Prandtl',PRANDTL, IFLAG,IERR)
C
C     Prandtl number Defaults to the value for air, if not
C             explicitely given
C
          IF( IFLAG .EQV. PETSC_FALSE )PRANDTL = 0.72d0
          IF(MY_PE.EQ.0)
     &    WRITE(IWUNIT,320)"Prandtl Number ",'PRANDTL',PRANDTL
C
C     test for isothermal walls
C
          CALL PetscOptionsHasName(PETSC_NULL_OPTIONS,
     +         PETSC_NULL_CHARACTER,'-isothermal_wall',IFLAG,IERR)
          IF( IFLAG .EQV. PETSC_TRUE )THEN
              IADIA = 1
              CALL PetscOptionsGetReal(PETSC_NULL_OPTIONS,
     +             PETSC_NULL_CHARACTER,'-wall_temperature',TWALL,
     &             IFLAG,IERR)
C
C     uses recovery temperature, if left unspecified in the
C     control file
C
              IF( IFLAG .EQV. PETSC_FALSE )THEN
                  IF(TURBULENT)THEN
                     RF = PRANDTL**(ONE/3.d0)
                  ELSE
                     RF = SQRT(PRANDTL)
                  ENDIF
                  TWALL = RF*HALF*GM1*M_INFTY*M_INFTY+ONE
                  TWALL = TWALL*TREF
                  IF(MY_PE.EQ.0)WRITE(IWUNIT,FMT=315)"Using recovery tem
     &perature"
              ENDIF
              IF(MY_PE.EQ.0)WRITE(IWUNIT,320)"Wall temperature ",'TWALL'
     &,TWALL
              CALL PetscOptionsGetInt(PETSC_NULL_OPTIONS,
     +             PETSC_NULL_CHARACTER,'-isothermal_bc_replace_eqn',
     &             IADIA_REPLACE_EQN,IFLAG,IERR)
              IF( IFLAG .EQV. PETSC_FALSE )THEN
                  IADIA_REPLACE_EQN = IADIA_REPLACE_ENER
              ELSE
                  LFLAG = (IADIA_REPLACE_EQN .EQ. IADIA_REPLACE_CONT)
     &               .OR. (IADIA_REPLACE_EQN .EQ. IADIA_REPLACE_ENER)
                  IF(.NOT.LFLAG)THEN
                      STOP 'Invalid IADIA_REPLACE_EQN'
                  ENDIF
              ENDIF
          ELSE
              IADIA = 0
              IADIA_REPLACE_EQN = -12438341
          ENDIF
          TWALL = TWALL/TREF
C
C     CONSTANTS FOR SUTHERLAND's LAW
C
          IF( NONDIM .EQ. EXTERNAL_FLOW )THEN
              C2 = M_INFTY*M_INFTY
              C1 = C2*M_INFTY
          ELSEIF( NONDIM .EQ. INTERNAL_FLOW )THEN
              C1 = GAM**(-1.5d0)
              C2 = GINV
          ENDIF
          C3 = SDIM/TREF
          C4 = ONE+C3
      ENDIF
C
  300 FORMAT('SETIBC Need to specify option: ',A25)
  310 FORMAT('Invalid argument: ',A22,' for option ',A18)
  315 FORMAT(5X,A30)
  320 FORMAT(5X,A30,2X,'(',A10,') =',E12.6)
C
      RETURN
      END
