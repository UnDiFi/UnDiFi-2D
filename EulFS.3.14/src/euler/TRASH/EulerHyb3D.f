      SUBROUTINE EulerHyb3D(IELEM,VCN,VCZ,NDIM,NOFVERT,NOFVAR,
     +                    NTURB,NODRES,TSTEP,STIFEL,VOLUME,PICARD,
     +                    ScalarScheme,MatrixScheme)
C
      IMPLICIT NONE 
C
C     Hybrid scheme
C
C     $Id: EulerHyb3D.f,v 1.4 2000/10/20 08:15:02 aldo Exp $
C
C
      INCLUDE 'paramt.h'
C
C    NEQMAX is the max. no. of equations (4 in 3D)
C           for the matrix scheme (solves for dp/ra,du,dv,dw)
C
      INTEGER IMACH
      INCLUDE 'constants'
C
C
      DOUBLE PRECISION DSTAK
      COMMON /CSTAK/DSTAK(1)
      INTEGER ISTAK(1)
      EQUIVALENCE(ISTAK(1),DSTAK(1))
C
C
      INCLUDE 'three'
      INCLUDE 'nloc'
      INCLUDE 'conv.com'
C
C
      INTEGER IELEM,NDIM,NOFVERT,NOFVAR,NTURB,NERR,IOPT
      INTEGER ICHPSI(4)
C
      INTRINSIC SQRT
C
C     ScalarScheme and MatrixScheme must be for the subsonic case.
      EXTERNAL ScalarScheme,MatrixScheme,PSI_scheme,LDASys_scheme
     +         ,N_scheme, NSys_scheme, LimNSys_scheme
C
C
      INTEGER IVAR,IVERT,JVERT,IPOIN
      DOUBLE PRECISION SUM,MACHMIN
      LOGICAL LFLAG,PICARD
      CHARACTER*72 ERRMSG
C
C
      DOUBLE PRECISION NODRES(*)
C
C     NODRES(1:NOFVAR,1:NOFVERT) is used to accumulate the
C         nodal residual in conserved variables and scatter
C         it to the RHS PETSc vector
C
C     TSTEP(1:NOFVERT) is used to accumulate the timestep
C         and then scatter it to the DT PETSc vector 
C
C
C
      DOUBLE PRECISION TSTEP(*)
C
C
C
      DOUBLE PRECISION VCZ(NOFVAR,NOFVERT),VCN(*),VOLUME,STIFEL(*)
C
      INTEGER I,L,LDA
C
C
      INTEGER IADDR
      IADDR(L,LDA,I) = (I-1)*LDA+L
C
C
C
C
C
caldo
      DO 20 JVERT = 1, NOFVERT
          IPOIN = ISTAK(IADDR(LCELNOD,NOFVERT,IELEM))
          ICHPSI(JVERT) = ISTAK(IADDR(LCHPSI,1,IPOIN))
   20 CONTINUE
caldo
C
c     Compute minimum Mach number:

      MACHMIN = ZERO

            DO  JVERT = 1 , NOFVERT

      IF(NDIM.EQ.3)
     &UAVG(5) = VCZ(5,JVERT)/VCZ(1,JVERT) ! z componenet of the velocity vector
      UAVG(4) = VCZ(4,JVERT)/VCZ(1,JVERT) ! y componenet of the velocity vector
      UAVG(3) = VCZ(3,JVERT)/VCZ(1,JVERT) ! x componenet of the velocity vector
      UAVG(2) = VCZ(2,JVERT)/VCZ(1,JVERT) ! Total Enthalpy
      UAVG(1) = VCZ(1,JVERT)*VCZ(1,JVERT) ! Density
C
      KINETIC = UAVG(3)*UAVG(3) + UAVG(4)*UAVG(4)
      IF( NDIM .EQ. 3 )KINETIC = KINETIC + UAVG(5)*UAVG(5)
      KINETIC = HALF * KINETIC
      ASQR = GM1 * ( UAVG(2) - KINETIC )
      IF( ASQR .LT. ZERO )THEN
         WRITE(ERRMSG(1:68),FMT=333)ASQR,IELEM
  333 FORMAT('EULERHYB3D',1X,'Negative averaged sound speed ',F7.3,
     &       ' in element ',I8)
         NERR = 6
         IOPT = 1
         CALL SETERR(ERRMSG(1:68),68,NERR,IOPT)
      ENDIF
      ABAR      = SQRT(ASQR)
c     MACHSQR   = DMAX1( TWO * KINETIC / ASQR, EPSQR_STAG) ! Barth & Deconinck
      MACHSQR   = TWO * KINETIC / ASQR
      MACH      = SQRT( MACHSQR )

      MACHMIN = dmin1(MACHMIN,MACH)
               
            END DO    

c     Compute average Mach number:

      DO 10 IVAR = 1 , NOFVAR
       SUM = ZERO
            DO 12 JVERT = 1 , NOFVERT
               SUM = SUM + VCZ( IVAR , JVERT )
   12       CONTINUE
      ZAVG(IVAR) = SUM / NOFVERT
   10 CONTINUE

      IF(NDIM.EQ.3)
     &UAVG(5) = ZAVG(5)/ZAVG(1) ! z componenet of the velocity vector
      UAVG(4) = ZAVG(4)/ZAVG(1) ! y componenet of the velocity vector
      UAVG(3) = ZAVG(3)/ZAVG(1) ! x componenet of the velocity vector
      UAVG(2) = ZAVG(2)/ZAVG(1) ! Total Enthalpy
      UAVG(1) = ZAVG(1)*ZAVG(1) ! Density
C
      KINETIC = UAVG(3)*UAVG(3) + UAVG(4)*UAVG(4)
      IF( NDIM .EQ. 3 )KINETIC = KINETIC + UAVG(5)*UAVG(5)
      KINETIC = HALF * KINETIC
      ASQR = GM1 * ( UAVG(2) - KINETIC )
      IF( ASQR .LT. ZERO )THEN
         WRITE(ERRMSG(1:68),FMT=333)ASQR,IELEM
         NERR = 6
         IOPT = 1
         CALL SETERR(ERRMSG(1:68),68,NERR,IOPT)
      ENDIF
      ABAR      = SQRT(ASQR)
c     MACHSQR   = DMAX1( TWO * KINETIC / ASQR, EPSQR_STAG) ! Barth & Deconinck
      MACHSQR   = TWO * KINETIC / ASQR
      MACH      = SQRT( MACHSQR )

      IMACH = ZERO
      DO JVERT = 1, NOFVERT
         IMACH = IMACH + ICHPSI(JVERT)
      END DO


      if (MACH.lt.ONE) then

c     CALL EulerIIbis(IELEM,VCN,VCZ,NDIM, 
      CALL EulerVII(IELEM,VCN,VCZ,NDIM, 
     +NOFVERT,NOFVAR,NTURB,NODRES,TSTEP,
     +STIFEL,VOLUME,PICARD, ScalarScheme,MatrixScheme)

      else

c     CALL EulerIIbis(IELEM,VCN,VCZ,NDIM, 
      CALL EulerVII(IELEM,VCN,VCZ,NDIM, 
     +NOFVERT,NOFVAR,NTURB,NODRES,TSTEP,
     +STIFEL,VOLUME,PICARD, PSI_scheme,LimNSys_scheme)

      end if

      END
C
