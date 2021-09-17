      SUBROUTINE SETIBC(NDIM,NOFVAR)
C
C     This routine SETs Initial and Boundary Conditions
C
      IMPLICIT NONE
C
C     $Id: setibc.F,v 1.7 2000/12/14 12:17:51 aldo Exp abonfi $
C
C
#include "include/finclude/petsc.h"
#include "include/finclude/petscviewer.h"
C
      PetscTruth IFLAG
C
      INCLUDE 'paramt.h'
      include 'bnd.h'
C
      INCLUDE 'implicit.h'
      INCLUDE 'visco.com'
      INCLUDE 'constants'
      INCLUDE 'conv.com'
      INCLUDE 'stream.com'
      INCLUDE 'flags.com'
      INCLUDE 'bnd'
      INCLUDE 'io.com'
      INCLUDE 'suther.com'
      INCLUDE 'verbose.com'
Caldo
      double precision zin(5),zout(5)
      common/aaa/zin,zout
Caldo
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
      INTEGER I,IERR,NVALS
      INTEGER NERR,IOPT
      DOUBLE PRECISION RTOT,PTOT,TTOT,TEMP,U,UREF,CHORD,HREF,T
      DOUBLE PRECISION A1,R1,T1,P1,H1,U1,V1,W1,VISCL,Z1,MACH1
      LOGICAL COMPRESSIBLE,INCOMPRESSIBLE,VISCOUS
C
      INTEGER I1MACH
      DOUBLE PRECISION DNRM2,PRESSC
 
      EXTERNAL PRESSC
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
C
      IF(MY_PE.EQ.0)WRITE (IWUNIT,FMT=2000) 
 2000 FORMAT (/,/,' INITIAL-BOUNDARY CONDITIONS ',/,1X,28('='),/)
C
C     *********************************************
C     Check non-dimensionalisation of the equations
C     *********************************************
C
      CALL OptionsGetString(PETSC_NULL_CHARACTER,
     +'-nondimensionalisation',STRINGA,IFLAG,IERR)
      IF( IFLAG .EQ. PETSC_FALSE )THEN
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
      IF( NONDIM .EQ. EXTERNAL_FLOW )THEN
C
          IF(COMPRESSIBLE)THEN
C
C         read freestream Mach number
C
              CALL OptionsGetDouble(PETSC_NULL_CHARACTER,
     +        '-freestream_Mach_number',M_INFTY,IFLAG,IERR)
              IF( IFLAG .EQ. PETSC_FALSE )THEN
                  WRITE(ERRMSG(1:56),300)'-freestream_Mach_number'
                  NERR = 4
                  IOPT = 1
                  CALL SETERR(ERRMSG,56,NERR,IOPT)
              ENDIF
              IF(MY_PE.EQ.0)
     &        WRITE(IWUNIT,320)'Freestream Mach number','M_INFTY',
     +M_INFTY
          ENDIF
C
          IF(COMPRESSIBLE.AND.VISCOUS)THEN
C
C         read Freestream Static Temperature
C
              CALL OptionsGetDouble(PETSC_NULL_CHARACTER,
     +        '-freestream_temperature',TREF,IFLAG,IERR)
              IF( IFLAG .EQ. PETSC_FALSE )THEN
                  TREF = 273.d0
              ENDIF
              IF(MY_PE.EQ.0)
     &        WRITE(IWUNIT,320)'Freestream Temperature','TREF',TREF
          ENDIF
C
C         read flow angles
C
          NVALS = NDIM
          CALL OptionsGetDoubleArray(PETSC_NULL_CHARACTER,
     +    '-flow_angles',FLOWDIR,NVALS,IFLAG,IERR)
          IF( IFLAG .EQ. PETSC_FALSE )THEN
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
              CALL OptionsGetDouble(PETSC_NULL_CHARACTER,'-Reynolds',RE,
     +        IFLAG,IERR)
              IF( IFLAG .EQ. PETSC_FALSE )THEN
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
              Q_INFTY = ONE
              HREF = (TWO+GM1*M_INFTY*M_INFTY)/(TWO*GM1*M_INFTY*M_INFTY)
              U_INFTY(1) = ONE
              U_INFTY(2) = HREF
              U_INFTY(3) = FLOWDIR(1)
              U_INFTY(4) = FLOWDIR(2)
              U_INFTY(5) = FLOWDIR(3)
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
          CALL OptionsGetDouble(PETSC_NULL_CHARACTER,
     +    '-outlet_static_pressure',POUTLET,IFLAG,IERR)
          IF( IFLAG .EQ. PETSC_FALSE )THEN
C
C     If one does not specify the value of the outlet
C     static pressure, this is
C     computed using the freestream conditions
C
              IF(COMPRESSIBLE)THEN
                  POUTLET = PRESSC(NDIM,U_INFTY(1))
              ELSE
                  POUTLET = ONE
              ENDIF
          ENDIF
C
C
      ELSEIF( NONDIM .EQ. INTERNAL_FLOW )THEN
C
C
      NVALS = NOFVAR 
      CALL OptionsGetDoubleArray(PETSC_NULL_CHARACTER,
     +'-inlet_parameter_vector',ZIN,NVALS,IFLAG,IERR)
      CALL OptionsGetDoubleArray(PETSC_NULL_CHARACTER,
     +'-outlet_parameter_vector',ZOUT,NVALS,IFLAG,IERR)
C
          CALL OptionsGetDouble(PETSC_NULL_CHARACTER,
     +    '-inlet_total_temperature',TTOT,IFLAG,IERR)
          IF( IFLAG .EQ. PETSC_FALSE )THEN
              WRITE(ERRMSG(1:56),300)'-inlet_total_temperature'
              NERR = 4
              IOPT = 1
              CALL SETERR(ERRMSG,56,NERR,IOPT)
          ENDIF
          TREF=TTOT
          IF(MY_PE.EQ.0)
     &    WRITE(IWUNIT,320)'Inlet Total Temperature','TREF',TREF
C
          CALL OptionsGetDouble(PETSC_NULL_CHARACTER,
     +    '-inlet_total_pressure',PTOT,IFLAG,IERR)
          IF( IFLAG .EQ. PETSC_FALSE )THEN
              WRITE(ERRMSG(1:56),300)'-inlet_total_pressure'
              NERR = 4
              IOPT = 1
              CALL SETERR(ERRMSG,56,NERR,IOPT)
          ENDIF
          IF(MY_PE.EQ.0)
     &    WRITE(IWUNIT,320)'Inlet Total Pressure','PTOT',PTOT
C
          CALL OptionsGetDouble(PETSC_NULL_CHARACTER,
     +    '-outlet_static_pressure',POUTLET, IFLAG,IERR)
 
          IF( IFLAG .EQ. PETSC_FALSE )THEN
              WRITE(ERRMSG(1:56),300)'-outlet_static_pressure'
              NERR = 4
              IOPT = 1
              CALL SETERR(ERRMSG,56,NERR,IOPT)
          ENDIF
          IF(MY_PE.EQ.0)
     &    WRITE(IWUNIT,320)'Outlet Static Pressure','POUTLET',POUTLET
          POUTLET=POUTLET/PTOT
C
          CALL OptionsGetDouble(PETSC_NULL_CHARACTER,
     +    '-reference_length',CHORD,IFLAG,IERR)
          IF( IFLAG .EQ. PETSC_FALSE )THEN
              WRITE(ERRMSG(1:56),300)'-reference_length'
              NERR = 4
              IOPT = 1
              CALL SETERR(ERRMSG,56,NERR,IOPT)
          ENDIF
          IF(MY_PE.EQ.0)
     &    WRITE(IWUNIT,320)'Reference Length ','CHORD',CHORD
C
C         read flow angles
C
          NVALS = NDIM
          CALL OptionsGetDoubleArray(PETSC_NULL_CHARACTER,
     +    '-flow_angles',FLOWDIR,NVALS,IFLAG,IERR)
          IF( IFLAG .EQ. PETSC_FALSE )THEN
              WRITE(ERRMSG(1:56),300)'-flow_angles'
              NERR = 4
              IOPT = 1
              CALL SETERR(ERRMSG,56,NERR,IOPT)
          ENDIF
C
	  IF ( COMPRESSIBLE )THEN
C
C     dimensional reference values: to be used to compute the
C     Reynolds number
C
          UREF = SQRT(RGAS*TTOT)
          RTOT = PTOT/(RGAS*TTOT)
          HREF = RGAS*TTOT
          RE = RTOT*UREF*CHORD/VISCL(TTOT)
          REINV = ONE/RE
          IF(MY_PE.EQ.0)WRITE(IWUNIT,320)'Reynolds Number ','RE',RE
C
C     the inlet mach number is only required to initialise
C     the flowfield
C
          CALL OptionsGetDouble(PETSC_NULL_CHARACTER,
     +    '-inlet_mach_number',MACH1,IFLAG,IERR)
          IF( IFLAG .EQ. PETSC_FALSE )THEN
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
          ELSEIF( INCOMPRESSIBLE )THEN
C
C     dimensional reference values: to be used to compute the
C     Reynolds number
C
C         compute (constant) density from total
C
C
          RTOT = PTOT/(RGAS*TTOT)
          UREF = SQRT(PTOT/RTOT)
caldoite(6,*)'RTOT;UREF;VISCL ',RTOT,UREF,VISCL(TTOT)
          RE = RTOT*UREF*CHORD/VISCL(TTOT)
          REINV = ONE/RE
          IF(MY_PE.EQ.0)WRITE(IWUNIT,320)'Reynolds Number ','RE',RE
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
C
          ENDIF
C
      ENDIF
C
C     A few options follow that are only meaningful for
C     compressible flows
C
      IF(COMPRESSIBLE.AND.VISCOUS)THEN
          CALL OptionsGetDouble(PETSC_NULL_CHARACTER,'-Prandtl',
     +    PRANDTL, IFLAG,IERR)
C
C     Prandtl number Defaults to the value for air, if not
C             explicitely given
C
          IF( IFLAG .EQ. PETSC_FALSE )PRANDTL = 0.72d0
          IF(MY_PE.EQ.0)
     &    WRITE(IWUNIT,320)"Prandtl Number ",'PRANDTL',PRANDTL
C
C     test for isothermal walls
C
          CALL OptionsHasName(PETSC_NULL_CHARACTER,'-isothermal_wall',
     +    IFLAG,IERR)
          IF( IFLAG .EQ. PETSC_TRUE )THEN
              IADIA = 1
              CALL OptionsGetDouble(PETSC_NULL_CHARACTER,
     +        '-wall_temperature',TWALL,IFLAG,IERR)
              IF( IFLAG .EQ. PETSC_FALSE )THEN
                  WRITE(ERRMSG,FMT=300)'-wall_temperature'
                  NERR = 4
                  IOPT = 1
                  CALL SETERR(ERRMSG,72,NERR,IOPT)
              ENDIF
          IF(MY_PE.EQ.0)
     &    WRITE(IWUNIT,320)"Wall temperature ",'TWALL',TWALL
          ELSE
              IADIA = 0
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
  320 FORMAT(5X,A30,2X,'(',A10,') =',E12.6)
C
      RETURN
      END
