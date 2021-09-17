      INTEGER       ISCHEME,JSCHEME,KSCHEME,ICHECK,KAN,DECOMP,IBCTYPE,
     &              SLIP_FREE_BC_TYPE
      LOGICAL TURBULENT,SEGREGATED,COUPLED,LTSQR,LDUMP(5),LREAD(2),
     &LAPLACE
C
C     $Id: flags.com,v 1.6 2013/06/04 14:38:31 abonfi Exp $
C
      COMMON/CFLAGS/ISCHEME,JSCHEME,KSCHEME,ICHECK,KAN,DECOMP,IBCTYPE,
     &SLIP_FREE_BC_TYPE,TURBULENT,SEGREGATED,COUPLED,LTSQR,LDUMP,LREAD,
     3LAPLACE
!23456789012345678901234567890123456789012345678901234567890123456789012
C
C	ISCHEME	scalar scheme
C	JSCHEME	system scheme
C	KSCHEME	scalar scheme for turbulence transport eqn
C	ICHECK  if <> 0 some tests are performed
C	KAN:	Type of analysis
C	DECOMP
C       IBCTYPE
C       SLIP_FREE_BC_TYPE
C       TURBULENT = .TRUE. when using a turbulence model
C       SEGREGATED
C       LDUMP(1) = .TRUE. when using -dump_nodal_residual
C       LDUMP(2) = .TRUE. when using -dump_pseudo_timestep
C       LDUMP(3) = .TRUE. when using -dump_jacobian_matrix
C       LDUMP(4) = .TRUE. when using -dump_boundary_fluxes
C       LDUMP(5) = .TRUE. when using -dump_integral
C       LREAD(1) = .TRUE. when reading inflow profiles
C       LREAD(2) = .TRUE. when a bflux000.dat file should be read
C
