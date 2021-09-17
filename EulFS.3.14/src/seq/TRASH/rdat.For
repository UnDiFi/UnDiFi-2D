      SUBROUTINE READVAL(MESHFILE,NGHBFILE,STARFILE,BAKFILE,VISCTFILE,
     +                   WDISTFILE,PBCSFILE,PTOTFILE,CLFILE,SKINFILE,
     +                   DATADIR,NTURB)
C
      IMPLICIT NONE
C
C     $Id: rdat.F,v 1.59 2008/02/25 09:45:04 abonfi Exp abonfi $
C
C
#include "include/finclude/petsc.h"
#include "include/finclude/petscviewer.h"
C
      PetscTruth IFLAG
C
      INCLUDE 'paramt.h'
      include 'bnd.h'
      INCLUDE 'bctype.h'
C
C
      INCLUDE 'implicit.h'
      INCLUDE 'visco.com'
      INCLUDE 'constants'
      INCLUDE 'conv.com'
      INCLUDE 'stream.com'
      INCLUDE 'chorin.com'
      INCLUDE 'scalar.com'
      INCLUDE 'flags.com'
      INCLUDE 'bodyf.com'
      INCLUDE 'periodic.com'
      INCLUDE 'turb.com'
      INCLUDE 'bnd.com'
      INCLUDE 'ibc2.com'
      INCLUDE 'io.com'
      INCLUDE 'verbose.com'
C
C
      INTEGER NTURB
      CHARACTER*(*) MESHFILE,NGHBFILE,STARFILE,BAKFILE,DATADIR,
     +              WDISTFILE,VISCTFILE,PBCSFILE,PTOTFILE,CLFILE,
     2              SKINFILE
C
      INTEGER MY_PE
      COMMON/MPICOM/MY_PE
C
      CHARACTER*80 STRINGA
      CHARACTER*72 ERRMSG
      CHARACTER*3 EXT
      INTEGER I,IERR,NVALS,KSPACE
      INTEGER NERR,IOPT
C
      INTEGER I1MACH
C
C     Executable Statements ..
C
C
C     Number of turbulence eqns.
C
      NTURB = 0
      TURBULENT = .FALSE.
      COUPLED = .FALSE.
      SEGREGATED = .FALSE.
C     when using turbulent viscosity multiplied
C     by \sqrt{\rho} set to .TRUE.
      LTSQR = .FALSE.
C
      WRITE(EXT,FMT='(I3.3)')MY_PE+1
C
      WRITE(NOUT,200)
#ifdef MPI
      IF(MY_PE.EQ.0)WRITE(IWUNIT,200)
#endif
C
C     ********************************************
C     Check options for time stepping
C     ********************************************
C
      CALL PetscOptionsGetString(PETSC_NULL_CHARACTER,'-timestepping',
     +                      STRINGA,IFLAG,IERR)
      IF( IFLAG .EQ. PETSC_FALSE )THEN
C
C     implicit time-stepping is the default
C
          TIMEIMPL = .TRUE.
      ELSE
          IF(STRINGA(1:8).EQ.'implicit') THEN
              TIMEIMPL = .TRUE.
          ELSEIF(STRINGA(1:8).EQ.'explicit') THEN
              TIMEIMPL = .FALSE.
          ELSE
              WRITE(ERRMSG,310)STRINGA,'-timestepping'
              NERR = 4
              IOPT = 1
!             SETERRQ(NERR,ERRMSG,ierr)
              SETERRQ(NERR,'invalid option for timestepping',ierr)
              CALL SETERR(ERRMSG,72,NERR,IOPT)
          ENDIF
      ENDIF
      CALL PetscOptionsGetString(PETSC_NULL_CHARACTER,'-timestep',
     +                      STRINGA,IFLAG,IERR)
      IF( IFLAG .EQ. PETSC_FALSE )THEN
C
C     local time-stepping is the default
C
          IGLOB = 1
      ELSE
          IF(STRINGA(1:3).EQ.'loc') THEN
              IGLOB = 1
          ELSEIF(STRINGA(1:3).EQ.'glo') THEN
              IGLOB = 0
          ELSE
              WRITE(ERRMSG,310)STRINGA,'-timestep'
              NERR = 4
              IOPT = 1
              CALL SETERR(ERRMSG,72,NERR,IOPT)
          ENDIF
      ENDIF
C
      CALL PetscOptionsGetReal(PETSC_NULL_CHARACTER,
     +                      '-time_step',DELT,IFLAG,IERR)
      IF( IFLAG .EQ. PETSC_FALSE )THEN
           WRITE(ERRMSG,310)STRINGA,'-time_step'
           NERR = 4
           IOPT = 1
           CALL SETERR(ERRMSG,72,NERR,IOPT)
      ENDIF
C
C 
C     ********************************************
C     Check options for linearization
C     ********************************************
C
      IF(TIMEIMPL)THEN
          CALL PetscOptionsGetString(PETSC_NULL_CHARACTER,
     +                          '-linearization',STRINGA,IFLAG,IERR)
          IF( IFLAG .EQ. PETSC_FALSE )THEN
C
C         Picard linearization is the default
C
              PICARD = .TRUE.
              NEWTON = .FALSE.
          ELSE
              IF(STRINGA(1:6).EQ.'picard') THEN
                  PICARD = .TRUE.
                  NEWTON = .FALSE.
              ELSEIF(STRINGA(1:6).EQ.'newton') THEN
                  PICARD = .FALSE.
                  NEWTON = .TRUE.
              ELSE
                  WRITE(ERRMSG,310)STRINGA,'-linearization'
                  NERR = 4
                  IOPT = 1
                  CALL SETERR(ERRMSG,72,NERR,IOPT)
              ENDIF
          ENDIF
      ELSE
          PICARD = .FALSE.
          NEWTON = .FALSE.
      ENDIF
C
C
C
C     maximum number of NON linear iterations
C
      CALL PetscOptionsGetInt(PETSC_NULL_CHARACTER,'-itmax',ITMAX,
     +IFLAG,IERR)
      IF( IFLAG .EQ. PETSC_FALSE )ITMAX = 0
C
C     convergence will be written to STDOUT each ISTMP iterations
C
      CALL PetscOptionsGetInt(PETSC_NULL_CHARACTER,'-istmp',ISTMP,
     +IFLAG,IERR)
      IF( IFLAG .EQ. PETSC_FALSE )ISTMP = 1
C
C     solution will be damped to disk each IBAK iterations
C
      CALL PetscOptionsGetInt(PETSC_NULL_CHARACTER,'-ibak',IBAK,IFLAG,
     +IERR)
      IF( IFLAG .EQ. PETSC_FALSE )IBAK = 10
C
C     variable upon which convergence is tested
C     need to be careful: at this stage the code does not know
C     how many (NOFVAR) variables we are solving for:
C     must ensure 1 <= IVCNVG <= NOFVAR 
C
      CALL PetscOptionsGetInt(PETSC_NULL_CHARACTER,'-ivcnvg',IVCNVG,
     +                   IFLAG,IERR)
      IF( IFLAG .EQ. PETSC_FALSE )IVCNVG = 1
C
C     threshold for terminating NON linear iterations
C
      CALL PetscOptionsGetReal(PETSC_NULL_CHARACTER,'-tolerance',
     +                      TOLER,IFLAG,IERR)
      IF( IFLAG .EQ. PETSC_FALSE )TOLER = -20.
C
      CALL PetscOptionsGetReal(PETSC_NULL_CHARACTER,'-cfl',CFL(1),
     +                      IFLAG,IERR)
C
      CALL PetscOptionsGetReal(PETSC_NULL_CHARACTER,'-cflmax',
     +                      CFLMAX(1),IFLAG,IERR)
C
C     ********************************************
C     Boundary colouring
C     ********************************************
C
      NVALS = NCOLOR+1
      CALL PetscOptionsGetIntArray(PETSC_NULL_CHARACTER,'-colors',
     +                        ICOLOR,NVALS,IFLAG,IERR)
      IF( IFLAG .EQ. PETSC_FALSE )THEN
          WRITE(ERRMSG,300)'-colors'
          NERR = 4
          IOPT = 1
          CALL SETERR(ERRMSG,72,NERR,IOPT)
      ENDIF
C
C     ********************************************
C     type of flow analysis
C     ********************************************
C
      CALL PetscOptionsGetString(PETSC_NULL_CHARACTER,'-equation',
     +                      STRINGA,IFLAG,IERR)
      IF( IFLAG .EQ. PETSC_FALSE )THEN
          WRITE(ERRMSG,300)'-equation'
          NERR = 4
          IOPT = 1
          CALL SETERR(ERRMSG,72,NERR,IOPT)
          WRITE(I1MACH(4),320)
      ENDIF
C
      IF( STRINGA(1:19) .EQ. 'advection-diffusion' )THEN
          KAN = +1
          CALL CheckPetscOptionsScalar
      ELSEIF( STRINGA(1:9) .EQ. 'advection' )THEN
          KAN = -1
          CALL CheckPetscOptionsScalar
      ELSEIF( STRINGA(1:5) .EQ. 'Euler' )THEN
          CALL CheckPetscOptionsEuler
      ELSEIF( STRINGA(1:13) .EQ. 'Navier-Stokes' )THEN
          CALL CheckPetscOptionsNS
      ELSEIF( STRINGA(1:4) .EQ. 'RANS' )THEN
          CALL CheckPetscOptionsNS
          CALL CheckPetscOptionsTurb
          NTURB = 1
          TURBULENT = .TRUE.
      ELSE
          WRITE(ERRMSG,310)STRINGA,'-equation'
          NERR = 4
          IOPT = 1
          CALL SETERR(ERRMSG,72,NERR,IOPT)
          WRITE(I1MACH(4),320)
      ENDIF
C
      IF(TURBULENT)THEN
C
C     boundary conditions for the turbulent variable(s)
C
          NVALS = NCOLOR+1
          CALL PetscOptionsGetIntArray(PETSC_NULL_CHARACTER,
     +                            '-turb_colors',ICOLOR_T,NVALS,
     +                            IFLAG,IERR)
          IF( IFLAG .EQ. PETSC_FALSE )THEN
              WRITE(ERRMSG(1:44),300)'-turb_colors'
              NERR = 4
              IOPT = 1
              CALL SETERR(ERRMSG,72,NERR,IOPT)
          ENDIF
C
C     time-stepping strategy for the turbulent variable(s)
C
          CALL PetscOptionsGetString(PETSC_NULL_CHARACTER,
     +                     '-turb_ts_type',
     +                     STRINGA,IFLAG,IERR)
          IF( IFLAG .EQ. PETSC_FALSE )THEN
              WRITE(ERRMSG(1:44),300)'-turb_ts_type'
              NERR = 4
              IOPT = 1
              CALL SETERR(ERRMSG,72,NERR,IOPT)
          ELSE
              IF(STRINGA(1:7).EQ.'coupled') THEN
                  COUPLED = .TRUE.
                  SEGREGATED = .FALSE.
              ELSEIF(STRINGA(1:10).EQ.'segregated') THEN
                  COUPLED = .FALSE.
                  SEGREGATED = .TRUE.
              ELSE
                  WRITE(ERRMSG,310)STRINGA,'-turb_ts_type'
                  NERR = 4
                  IOPT = 1
                  CALL SETERR(ERRMSG,72,NERR,IOPT)
              ENDIF
          ENDIF
      ENDIF
C
      IF(ABS(KAN).EQ.2)THEN
          CALL PetscOptionsGetReal(PETSC_NULL_CHARACTER,
     +         '-artificial_sound_speed',BETA,IFLAG,IERR)
          IF( IFLAG .EQ. PETSC_FALSE )BETA = 1.d0
      ENDIF 
C
C     for the incompressible UNpreconditioned eqns. no scalar scheme
C     is required, however a reasonable value for ISCHEME is still
C     required, since this is used in subroutine SCLSCH()
C
      IF(ABS(KAN).EQ.2.AND.DECOMP.EQ.7)THEN
         ISCHEME = 1
         GOTO 70 
      ENDIF
C
C     ********************************************
C     choose scalar and matrix schemes 
C     ********************************************
C
      CALL PetscOptionsGetString(PETSC_NULL_CHARACTER,'-scalar_scheme',
     +                      STRINGA,IFLAG,IERR)
C
      IF( IFLAG .EQ. PETSC_FALSE )THEN 
          WRITE(ERRMSG(1:44),300)'-scalar_scheme'
          NERR = 4
          IOPT = 1
          CALL SETERR(ERRMSG,72,NERR,IOPT)
      ENDIF
C
C     be careful with schemes like: N,NL,etc. 
C     test first the longer ones
C
      IF( STRINGA(1:3).EQ.'NL2' )THEN 
          ISCHEME = 7
      ELSEIF( STRINGA(1:2).EQ.'NL' )THEN 
          ISCHEME = 4
      ELSEIF( STRINGA(1:2).EQ.'NS' )THEN 
          ISCHEME = 9
      ELSEIF( STRINGA(1:2).EQ.'LW' )THEN 
          ISCHEME = 8
      ELSEIF( STRINGA(1:1).EQ.'N' )THEN 
          ISCHEME = 1
      ELSEIF( STRINGA(1:3).EQ.'PSI' )THEN 
          ISCHEME = 2
      ELSEIF( STRINGA(1:4).EQ.'LDA2' )THEN 
          ISCHEME = 11
      ELSEIF( STRINGA(1:3).EQ.'LDA' )THEN 
          ISCHEME = 3
      ELSEIF( STRINGA(1:3).EQ.'FCT' )THEN 
          ISCHEME = 6
      ELSEIF( STRINGA(1:2).EQ.'FV' )THEN 
          ISCHEME = 5
      ELSEIF( STRINGA(1:4).EQ.'SUPG' )THEN 
          ISCHEME = 10
      ELSE
          WRITE(ERRMSG,FMT=310)STRINGA(1:22),'-scalar_scheme'
          NERR = 4
          IOPT = 1
          CALL SETERR(ERRMSG,72,NERR,IOPT)
      ENDIF
   70 CONTINUE
C
C     if it is a scalar problem, we do not have to care about
C        the matrix scheme
C
      IF(ABS(KAN).EQ.1)THEN
          JSCHEME = 1
          GOTO 140
      ENDIF
C
C     ********************************************
C     choose matrix scheme
C     ********************************************
C
      CALL PetscOptionsGetString(PETSC_NULL_CHARACTER,'-matrix_scheme',
     +                      STRINGA,IFLAG,IERR)
      IF( IFLAG .EQ. PETSC_FALSE )THEN 
          WRITE(ERRMSG(1:44),300)'-matrix_scheme'
          NERR = 4
          IOPT = 1
          CALL SETERR(ERRMSG,72,NERR,IOPT)
      ENDIF
C
C     be careful with schemes like: LDA,LDA2,etc. 
C     test first the longer ones
C
      IF    ( STRINGA(1:1).EQ.'N' )THEN 
          JSCHEME = 1
      ELSEIF( STRINGA(1:1).EQ.'B' )THEN 
          JSCHEME = 10
      ELSEIF( STRINGA(1:4).EQ.'LDAF' )THEN 
          JSCHEME = 2
      ELSEIF( STRINGA(1:4).EQ.'LDA2' )THEN 
          JSCHEME = 6
      ELSEIF( STRINGA(1:3).EQ.'LDA' )THEN 
          JSCHEME = 3
      ELSEIF( STRINGA(1:3).EQ.'CDS' )THEN 
          JSCHEME = 7
      ELSEIF( STRINGA(1:3).EQ.'LW2' )THEN 
          JSCHEME = 4
      ELSEIF( STRINGA(1:4).EQ.'LimN' )THEN
          JSCHEME = 9
      ELSEIF( STRINGA(1:2).EQ.'LW' )THEN 
          JSCHEME = 8
      ELSEIF( STRINGA(1:2).EQ.'FV' )THEN 
          JSCHEME = 5
      ELSE
          WRITE(ERRMSG,FMT=310)STRINGA(1:22),'-matrix_scheme'
          NERR = 4
          IOPT = 1
          CALL SETERR(ERRMSG,72,NERR,IOPT)
      ENDIF
C
C
  140 CONTINUE
C
C     PetscOptions that do NOT depend on the type of eqns. being solved 
C
C     ********************************************
C     debugging option
C     ********************************************
C
      CALL PetscOptionsHasName(PETSC_NULL_CHARACTER,'-check',
     +IFLAG,IERR)
      IF( IFLAG .EQ. PETSC_TRUE )THEN
          ICHECK = 1
      ELSE 
          ICHECK = 0
      ENDIF
C
C     *******************************************
C     meshfile(s)
C     *******************************************
C
      CALL PetscOptionsGetString(PETSC_NULL_CHARACTER,'-data_dir',
     +                      DATADIR,IFLAG,IERR)
      IF( IFLAG .EQ. PETSC_FALSE )THEN
          WRITE(ERRMSG(1:44),300)'-data_dir'
          NERR = 4
          IOPT = 1
          CALL SETERR(ERRMSG,72,NERR,IOPT)
      ENDIF
C
      KSPACE = INDEX(DATADIR,' ')
C
#ifdef MPI
C
      MESHFILE(1:KSPACE+10)  = DATADIR(1:KSPACE-1)//'mesh'//ext//'.dat'
      NGHBfILE(1:KSPACE+10)  = DATADIR(1:KSPACE-1)//'neig'//ext//'.dat'
      WDISTFILE(1:KSPACE+10) = DATADIR(1:KSPACE-1)//'wdis'//ext//'.dat'
      PBCSFILE(1:KSPACE+10) = DATADIR(1:KSPACE-1)//'pbcs'//ext//'.dat'
      PTOTFILE(1:KSPACE+10) = DATADIR(1:KSPACE-1)//'ptot'//ext//'.dat'
      CLFILE(1:KSPACE+10) = DATADIR(1:KSPACE-1)//'clin'//ext//'.dat'
C
C     file containing the flow variables and viscosity are kept under
C     the current (working) directory
C
      VISCTFILE(1:11) = 'vist'//ext//'.dat'
      BAKFILE(1:11)   = 'file'//ext//'.dat'
      STARFILE(1:11)  = 'file'//ext//'.dat'
      SKINFILE(1:11) = 'skin'//ext//'.dat'
C
#else
C
      MESHFILE(1:KSPACE+10) = DATADIR(1:KSPACE-1)//'file001.dat'
      NGHBfILE(1:KSPACE+10) = DATADIR(1:KSPACE-1)//'file002.dat'
      STARFILE(1:KSPACE+10) = DATADIR(1:KSPACE-1)//'file003.dat'
      PBCSFILE(1:KSPACE+10) = DATADIR(1:KSPACE-1)//'file004.dat'
      WDISTFILE(1:KSPACE+10) = DATADIR(1:KSPACE-1)//'file019.dat'
      PTOTFILE(1:KSPACE+10) = DATADIR(1:KSPACE-1)//'file005.dat'
      CLFILE(1:KSPACE+10) = DATADIR(1:KSPACE-1)//'file006.dat'
      BAKFILE(1:11) = 'file010.dat'
      VISCTFILE(1:11) = 'file015.dat'
      SKINFILE(1:11) = 'file016.dat'
#endif
C
C
C     ********************************************
C     restart from a previous solution
C     ********************************************
C
      IF(IABS(KAN).EQ.1)GOTO 60
C
      CALL PetscOptionsGetString(PETSC_NULL_CHARACTER,
     +                      '-restart_file',STRINGA,IFLAG,IERR)
C
      IF( IFLAG .EQ. PETSC_TRUE )THEN

      ISTART = 1
C
C     in the parallel case the restart file cannot be changed from 
C     file$np$.dat
C
#ifndef MPI
          STARFILE = STRINGA
#endif
      ELSE 
          ISTART = 0
      ENDIF
   60 CONTINUE
cxxxx
!         TURBULENT = .TRUE.
!         SEGREGATED = .FALSE.
!         COUPLED = .TRUE.
cxxxx
C
C 
C     ********************************************
C     Body forces
C     ********************************************
C
      NVALS = 3
      CALL PetscOptionsGetRealArray(PETSC_NULL_CHARACTER,
     +'-body_force',GRAV,NVALS,IFLAG,IERR)
          IF( IFLAG .EQ. PETSC_FALSE )THEN
              GRAV(1) = ZERO
              GRAV(2) = ZERO
              GRAV(3) = ZERO
          ENDIF
C
C     type of boundary conditions: inflow/outflow
C
      CALL PetscOptionsGetString(PETSC_NULL_CHARACTER,'-bc_type',
     +                      STRINGA,IFLAG,IERR)
      IF( IFLAG .EQ. PETSC_FALSE )THEN
          WRITE(ERRMSG(1:44),300)'-bc_type'
          NERR = 4
          IOPT = 1
          CALL SETERR(ERRMSG,44,NERR,IOPT)
      ELSE
          IF(    STRINGA(1:4).EQ.'weak')THEN 
             IBCTYPE = BC_TYPE_WEAK
          ELSEIF(STRINGA(1:6).EQ.'strong')THEN 
             IBCTYPE = BC_TYPE_STRONG
          ELSE 
             WRITE(ERRMSG,FMT=310)STRINGA(1:20),'-bc_type'
             NERR = 4
             IOPT = 1
             CALL SETERR(ERRMSG,72,NERR,IOPT)
          ENDIF
      ENDIF
C
C     type of periodic conditions: annular or linear
C
      CALL PetscOptionsGetInt(PETSC_NULL_CHARACTER,'-number_of_blades',
     +                   NBLADES,IFLAG,IERR)
      IF( IFLAG .EQ. PETSC_FALSE )THEN
          ANNULAR = .FALSE.
      ELSE
          ANNULAR = .TRUE.
          PITCH = TWO*PETSC_PI/REAL(NBLADES)
          COSALPHA = COS(PITCH)
          SINALPHA = SIN(PITCH)
      ENDIF
C
C     type of periodic conditions: annular or linear
C
      CALL PetscOptionsHasName(PETSC_NULL_CHARACTER,
     +'-radial_equilibrium',IFLAG,IERR)
      IF( IFLAG .EQ. PETSC_FALSE )THEN
          RAD_EQUI = .FALSE.
      ELSE
          RAD_EQUI = .TRUE.
          IF(.NOT.ANNULAR)THEN
             STOP 'Cannot use ISRE on linear cascades'
          ENDIF
      ENDIF
C
      RETURN
C
   10 FORMAT(A)
   66 FORMAT(7(I10))
C
  300 FORMAT('Need to specify option: ',A20)
  310 FORMAT('Invalid argument: ',A22,' for option ',A18)
C
  320 FORMAT(5X,'Valid options are:'/,
     +       15X,'advection'/, 
     +       15X,'advection-diffusion'/, 
     +       15X,'Euler'/, 
     +       15X,'Navier-Stokes'/, 
     +       15X,'RANS')
C
  325 FORMAT(5X,'Valid options are:'/,
     +       15X,'weak'/, 
     +       15X,'mirror')
  330 FORMAT(5X,'Boundary type ',I2,' associated with color ',I2,
     &' does NOT exist')
  200 FORMAT(/,
     115X,'#######                 #######  #####',/,
     215X,'#        #    #  #      #       #     #',/,
     315X,'#        #    #  #      #       #',/,
     415X,'#####    #    #  #      #####    #####',/,
     515X,'#        #    #  #      #             #',/,
     615X,'#        #    #  #      #       #     #',/,
     715X,'#######   ####   ###### #        #####',/)
      END
C
      SUBROUTINE CheckPetscOptionsScalar
C
      IMPLICIT NONE 
C
#include "include/finclude/petsc.h"
C
C
      PetscTruth IFLAG
C
      INCLUDE 'paramt.h'
      INCLUDE 'implicit.h'
      INCLUDE 'visco.com'
      INCLUDE 'flags.com'
      INCLUDE 'conv.com'
      INCLUDE 'scalar.com'
C
      INTEGER I1MACH
      INTEGER IERR
      INTEGER NERR,IOPT
      CHARACTER*80 STRINGA
      CHARACTER*72 ERRMSG
C
C     ********************************************
C     check for testcase
C
C     the convection speed is computed depending on the testcase,
C     similarly for the source term, if present 
C     ********************************************
C
      CALL PetscOptionsGetInt(PETSC_NULL_CHARACTER,'-testcase',ICASE,
     +                   IFLAG,IERR)
      IF(IFLAG.EQ.PETSC_FALSE)THEN
          WRITE(I1MACH(4),300)'-testcase'
          WRITE(I1MACH(4),330)
          WRITE(I1MACH(4),340)
          NERR = 4
          IOPT = 1
          CALL SETERR(ERRMSG,72,NERR,IOPT)
      ENDIF
C
C     ********************************************
C     force restart to read boundary/initial data
C     ********************************************
C
      ISTART = 1 
C
      IF(KAN.EQ.1)THEN
C
C     read diffusion parameter
C
          CALL PetscOptionsGetReal(PETSC_NULL_CHARACTER,'-Reynolds',
     +                          RE,IFLAG,IERR)
          IF( IFLAG .EQ. PETSC_FALSE )THEN
              WRITE(ERRMSG(1:44),300)'-Reynolds'
              NERR = 4
              IOPT = 1
              CALL SETERR(ERRMSG,72,NERR,IOPT)
          ENDIF
          REINV=1.d0/RE 
      ELSEIF(KAN.EQ.-1)THEN
          RE = 1.E+38
          REINV=0.d0
      ELSE
          WRITE(I1MACH(4),FMT=*)
          STOP
      ENDIF
C 
  300 FORMAT(5X,'Need to specify option: ',A20)
  310 FORMAT('Invalid argument: ',A22,' for option ',A18)
  320 FORMAT(5X,'in order to set up initial/boundary conditions')
  330 FORMAT(5X,'Implemented testcases in TWO dimensions are:'/,
     +       15X,'1. Linear convection a = (2.,1.)'/, 
     +       15X,'2. Burger equation   a = (u.,1.)'/, 
     +       15X,'3. Linear convection a = (1.,1.)'/, 
     +15X,'4. Linear convection a = (1.,1.) + source term 4*x - 2*y'/, 
     +15X,'5. Linear convection a = (1.,1.) + source term '/, 
     +15X,'6. Linear convection a = (cos(d),sin(d)) + source term')
C
  340 FORMAT(5X,'Implemented testcases in THREE dimensions are:'/,
     +       15X,'1. Linear convection a = (0.75,0.875,1.)'/, 
     +       15X,'2. Linear convection a = (z,.2,-x)')
C
      RETURN
      END
C
      SUBROUTINE CheckPetscOptionsEuler
C
      IMPLICIT NONE 
C
#include "include/finclude/petsc.h"
C
C
      PetscTruth IFLAG
C
      INCLUDE 'visco.com'
      INCLUDE 'flags.com'
      INCLUDE 'implicit.h'
C
      INTEGER I1MACH
      INTEGER IERR,NERR,IOPT
      CHARACTER*80 STRINGA
      CHARACTER*72 ERRMSG
      CHARACTER*1 ANSW
C
C
C     type of fluid: compressible OR incompressible
C
      CALL PetscOptionsGetString(PETSC_NULL_CHARACTER,'-fluid',
     +                      STRINGA,IFLAG,IERR)
      IF( IFLAG .EQ. PETSC_FALSE )THEN
          WRITE(ERRMSG,FMT=310)STRINGA(1:20),'-fluid'
          NERR = 4
          IOPT = 1
          CALL SETERR(ERRMSG,72,NERR,IOPT)
      ELSE
          IF(    STRINGA(1:12).EQ.'compressible')THEN 
             KAN = -4
          ELSEIF(STRINGA(1:14).EQ.'incompressible')THEN 
             KAN = -2
          ELSE 
             WRITE(ERRMSG,FMT=310)STRINGA(1:20),'-fluid'
             NERR = 4
             IOPT = 1
             CALL SETERR(ERRMSG,72,NERR,IOPT)
          ENDIF
      ENDIF
C
C     Test option for preconditioner
C
      CALL PetscOptionsGetString(PETSC_NULL_CHARACTER,
     +         '-preconditioning',STRINGA,IFLAG,IERR)
      IF( IFLAG .EQ. PETSC_FALSE )THEN
          WRITE(ERRMSG,300)'-preconditioning'
          NERR = 4
          IOPT = 1
          CALL SETERR(ERRMSG,72,NERR,IOPT)
      ELSE
          ANSW = STRINGA(1:1)
          IF(ANSW.EQ.'y'.OR.ANSW.EQ.'Y')THEN 
              IF    (KAN.EQ.-4)THEN
                  DECOMP = 3
              ELSEIF(KAN.EQ.-2)THEN
                  DECOMP = 9
              ENDIF 
          ELSEIF(ANSW.EQ.'h'.OR.ANSW.EQ.'H')THEN
!                 DECOMP = 4
                  DECOMP = 2
          ELSEIF(ANSW.EQ.'t'.OR.ANSW.EQ.'T')THEN
                   DECOMP = 5
          ELSEIF(ANSW.EQ.'n'.OR.ANSW.EQ.'N')THEN 
              IF    (KAN.EQ.-4)THEN
                  DECOMP = 7
              ELSEIF(KAN.EQ.-2)THEN
                  DECOMP = 8
              ENDIF 
          ELSE 
              WRITE(ERRMSG,310)ANSW,'-preconditioning'
              NERR = 4
              IOPT = 1
              CALL SETERR(ERRMSG,72,NERR,IOPT)
          ENDIF
      ENDIF
      RE = 0.d0
      REINV = 1.D+38
C 
C
  300 FORMAT('Need to specify option: ',A30)
  310 FORMAT('Invalid argument: ',A22,' for option ',A18)
C
      RETURN
      END 

      SUBROUTINE CheckPetscOptionsNS
C
      IMPLICIT NONE 
C
#include "include/finclude/petsc.h"
C
C
      PetscTruth IFLAG
C
      INCLUDE 'visco.com'
      INCLUDE 'flags.com'
      INCLUDE 'implicit.h'
C
      INTEGER I1MACH
      INTEGER IERR
      INTEGER NERR,IOPT
      CHARACTER*80 STRINGA
      CHARACTER*72 ERRMSG
      CHARACTER*1 ANSW
C
C     type of fluid: compressible OR incompressible
C
      CALL PetscOptionsGetString(PETSC_NULL_CHARACTER,'-fluid',
     +                      STRINGA,IFLAG,IERR)
      IF( IFLAG .EQ. PETSC_FALSE )THEN
          WRITE(ERRMSG(1:44),300)'-fluid'
          NERR = 4
          IOPT = 1
          CALL SETERR(ERRMSG,44,NERR,IOPT)
      ELSE
          IF(    STRINGA(1:12).EQ.'compressible')THEN 
             KAN = 4
             DECOMP = 7
          ELSEIF(STRINGA(1:14).EQ.'incompressible')THEN 
             KAN = 2
             DECOMP = 8
          ELSE 
             WRITE(ERRMSG,FMT=310)STRINGA(1:20),'-fluid'
             NERR = 4
             IOPT = 1
             CALL SETERR(ERRMSG,72,NERR,IOPT)
          ENDIF
      ENDIF
C
C     Test option for preconditioner
C
      CALL PetscOptionsGetString(PETSC_NULL_CHARACTER,
     +         '-preconditioning',STRINGA,IFLAG,IERR)
      IF( IFLAG .EQ. PETSC_FALSE )THEN
          WRITE(ERRMSG,300)'-preconditioning [y,Y,n,N]'
          NERR = 4
          IOPT = 1
          CALL SETERR(ERRMSG,72,NERR,IOPT)
      ELSE
          ANSW = STRINGA(1:1)
          IF(ANSW.EQ.'y'.OR.ANSW.EQ.'Y')THEN 
              IF    (KAN.EQ.4)THEN
                  DECOMP = 3
              ELSEIF(KAN.EQ.2)THEN
                  DECOMP = 9
              ENDIF 
          ELSEIF(ANSW.EQ.'h'.OR.ANSW.EQ.'H')THEN
                  DECOMP = 4
          ELSEIF(ANSW.EQ.'t'.OR.ANSW.EQ.'T')THEN
                   DECOMP = 5
          ELSEIF(ANSW.EQ.'n'.OR.ANSW.EQ.'N')THEN 
              IF    (KAN.EQ.4)THEN
                  DECOMP = 7
              ELSEIF(KAN.EQ.2)THEN
                  DECOMP = 8
              ENDIF 
          ELSE 
              WRITE(I1MACH(4),310)ANSW,'-preconditioning'
              STOP
          ENDIF
      ENDIF
      RETURN
  300 FORMAT('Need to specify option: ',A48)
  310 FORMAT('Invalid argument: ',A22,' for option ',A18)
      END 
C
      SUBROUTINE CheckPetscOptionsTurb
C
      IMPLICIT NONE 
C
#include "include/finclude/petsc.h"
C
C
      PetscTruth IFLAG
C
      INCLUDE 'paramt.h'
      INCLUDE 'implicit.h'
      INCLUDE 'visco.com'
      INCLUDE 'flags.com'
      INCLUDE 'conv.com'
      INCLUDE 'turb.com'
C
      INTEGER I1MACH
      INTEGER IERR
      INTEGER NERR,IOPT
      CHARACTER*80 STRINGA
      CHARACTER*72 ERRMSG
C
C
C****************************************************************
C     
C     Time stepping (turbulence eqn.)
C     
C****************************************************************
C
      CALL PetscOptionsGetReal(PETSC_NULL_CHARACTER,'-cfl_max',
     +                      CFLMAX(2),IFLAG,IERR)
C
      CALL PetscOptionsGetReal(PETSC_NULL_CHARACTER,'-cfl_start',
     +                      CFL(2),IFLAG,IERR)
C
      CALL PetscOptionsGetReal(PETSC_NULL_CHARACTER,'-cfl_ratio',
     +                      CFLRATIO,IFLAG,IERR)
C
C     ********************************************
C     turbulence modeling (Spalart-Allmaras model)
C     ********************************************
C
C     Trip term on/off
C
      CALL PetscOptionsHasName(PETSC_NULL_CHARACTER,'-trip_term',
     +IFLAG,IERR)
C
      IF( IFLAG .EQ. PETSC_TRUE )THEN
              TTFLAG = 1
C
          CALL PetscOptionsGetInt(PETSC_NULL_CHARACTER,'-trip_cell',
     +                       TTELEM,IFLAG,IERR)
          CALL PetscOptionsGetReal(PETSC_NULL_CHARACTER,
     +                       '-trip_spacing',TDXT,IFLAG,IERR)
C
      ELSE
              TTFLAG = 0
      ENDIF
C
C     number of sub iterations on the turbulence equation
C
      CALL PetscOptionsGetInt(PETSC_NULL_CHARACTER,'-turb_sub_iter',
     +                   NSUBIT,IFLAG,IERR)
C
      IF( IFLAG .EQ. PETSC_FALSE .OR. NSUBIT .LE. 0 )NSUBIT = 1
C 
  300 FORMAT(5X,'Need to specify option: ',A20)
  310 FORMAT('Invalid argument: ',A22,' for option ',A18)
      RETURN
      END
C
