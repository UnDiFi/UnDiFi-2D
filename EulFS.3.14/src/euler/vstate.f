C
C     $Id: vstate.f,v 1.13 2013/04/30 07:42:54 abonfi Exp $
C
C     A bunch of subroutines with the same calling sequence:
C     subroutine ghost*(zin,zstar,work,ndim)
C
C     used to set the starred (virtual) state for inflow/outflow b.c.s
C
C     zin(1:NOFVAR) in input
C                   is the dependent variable in the boundary node
C                   un-changed on return
C     zstar(1:NOFVAR) upon return
C                     is the dependent variable in the ghost node
C
      subroutine ghost2vii(zin,zstar,work,ndim)
      IMPLICIT NONE
C
C     ghost state for subsonic outflow b.c. (compressible flow):
C     +) static pressure is taken from the b.c.
C     +) all other primitive variables are taken from the interior
C
      integer ndim
      double precision zin(*),zstar(*),work(*)
C
      INCLUDE'constants.h'
      INCLUDE 'pfcgas.com'
      double precision temp,POUTLET
C
      POUTLET = WORK(1)
C
C     compute in the boundary node
C
      ZSTAR(1) = ZIN(1)
      ZSTAR(3) = ZIN(3)
      ZSTAR(4) = ZIN(4)
      IF(NDIM.EQ.3) ZSTAR(5) = ZIN(5)
C
      TEMP       = ZIN(3)*ZIN(3) + ZIN(4)*ZIN(4)
      IF( NDIM .EQ. 3 )TEMP = TEMP + ZIN(5)*ZIN(5)
      TEMP = HALF * TEMP
      ZSTAR(2) = (GOGM1*POUTLET+TEMP)/ZSTAR(1)
C
      RETURN
      END
C
      subroutine ghost8vii(zin,zstar,work,ndim)
      IMPLICIT NONE
C
C     ghost state for subsonic inflow b.c. (compressible flow):
C     +) Mach number is taken from the interior
C     +) total pressure, total temperature and flow angles are imposed
C
      integer ndim
      double precision zin(*),zstar(*),work(*)
      double precision uavg(5)
C
      INCLUDE'constants.h'
      INCLUDE'ibc8.com'
      INCLUDE 'pfcgas.com'
      integer nerr,iopt
      double precision asqr,abar,mach,machsqr,kinetic,pres,dens,
     +temp,dens0
      character*67 errmsg
C
C     work(1:*):
C     p0/p_ref,t0/t_ref,unused,n_x,n_y,n_z
C
C     compute in the boundary node
C
      IF(NDIM.EQ.3)
     &UAVG(5) = ZIN(5)/ZIN(1) ! z component of the velocity vector
      UAVG(4) = ZIN(4)/ZIN(1) ! y component of the velocity vector
      UAVG(3) = ZIN(3)/ZIN(1) ! x component of the velocity vector
      UAVG(2) = ZIN(2)/ZIN(1) ! Total Enthalpy
      UAVG(1) = ZIN(1)*ZIN(1) ! Density
C
      KINETIC = UAVG(3)*UAVG(3) + UAVG(4)*UAVG(4)
      IF( NDIM .EQ. 3 )KINETIC = KINETIC + UAVG(5)*UAVG(5) 
      KINETIC = HALF * KINETIC
      ASQR = GM1 * ( UAVG(2) - KINETIC )
      IF( ASQR .LT. ZERO )THEN
         WRITE(ERRMSG(1:67),FMT=333)ASQR,-10
  333 FORMAT('GHOST8VII',1X,'Negative averaged sound speed ',F7.3,
     &       ' in element ',I8)
         NERR = 6
         IOPT = 1
         CALL SETERR(ERRMSG(1:67),67,NERR,IOPT)
      ENDIF
      ABAR      = SQRT(ASQR)
      MACHSQR   = TWO * KINETIC / ASQR
      MACH      = SQRT( MACHSQR )
C
      TEMP = (ONE+HALF*GM1*MACHSQR)
C
C     rho/rho^0 = (1+0.5*gm1*mach**2)**expn * p^0/T^0
C
C     work(2) = p^0
C     work(3) = T^0
C
      DENS0 = WORK(2)/WORK(3)
      DENS = ( TEMP**(-ONE/GM1) ) * DENS0
      PRES = ( ONE/TEMP**GOGM1 ) * WORK(2)
      ABAR = SQRT(GAM*PRES/DENS)
      ZSTAR(1) = SQRT(DENS)
      TEMP = ZSTAR(1) * ABAR * MACH
      ZSTAR(2) = ZSTAR(1) * GOGM1 * WORK(3)
      ZSTAR(3) = TEMP * WORK(4)
      ZSTAR(4) = TEMP * WORK(5)
      IF(NDIM.EQ.3)ZSTAR(5) = TEMP * WORK(6)
C
      RETURN
      END
C
      subroutine ghost8viii(zin,zstar,work,ndim)
      IMPLICIT NONE
C
C     ghost state for subsonic inflow b.c. (INcompressible flow):
C     +) velocity magnitude is taken from the interior
C     +) total pressure is imposed
C
      integer ndim
      double precision zin(*),zstar(*),work(*)
      integer nerr,iopt
      character*67 errmsg
C
      INCLUDE'constants.h'
      INCLUDE'ibc8.com'
      INTEGER MY_PE
      COMMON/MPICOM/MY_PE
      double precision temp
      INTEGER I1MACH
      EXTERNAL I1MACH
C
C     velocity magnitude is taken from the interior
C
      TEMP = ZIN(2)*ZIN(2) + ZIN(3)*ZIN(3)
      IF(NDIM.EQ.3) TEMP = TEMP +  ZIN(4)*ZIN(4)
C
      ZSTAR(1) = WORK(2) - HALF * TEMP
C
C     flow angles are imposed
C
      TEMP = SQRT(TEMP)
      ZSTAR(2) = TEMP * WORK(4)
      ZSTAR(3) = TEMP * WORK(5)
      IF(NDIM.EQ.3)ZSTAR(4) = TEMP * WORK(6)
C
      RETURN
      END
C  
C
      subroutine ghost2viii(zin,zstar,work,ndim)
      IMPLICIT NONE
C
C     ghost state for subsonic outflow b.c.:
C     +) velocity magnitude is taken from the interior
C     +) static pressure is imposed
C
      integer ndim
      double precision zin(*),zstar(*),work(*)
C
      INCLUDE'constants.h'
C
      ZSTAR(1) = WORK(1)
      ZSTAR(2) = ZIN(2)
      ZSTAR(3) = ZIN(3)
      IF(NDIM.EQ.3)ZSTAR(4) = ZIN(4)
C
      RETURN
      END
C 
      subroutine ghost2vii4Ar(zin,zstar,work,ndim)
      IMPLICIT NONE
C
C     ghost state for subsonic outflow b.c. (plasma flow):
C     +) static pressure is taken from the b.c.
C     +) all other primitive variables are taken from the interior
C
      integer ndim,isp
      double precision zin(*),zstar(*),work(*)
C
      INCLUDE 'constants.h'
      INCLUDE 'plasma.h'
      INCLUDE 'chem.h'
      INCLUDE 'commonv.inc'
      INCLUDE 'dofs.com'
      INCLUDE 'pfcgas.com'
      INclude 'streamplasma.com'
      double precision temp,POUTLET,Zrho,HFTOT
C
      POUTLET = WORK(1)
C
C     compute in the boundary node
C
      DO ISP = 1,NSP 
         ZSTAR(ISP) = ZIN(ISP) 
      ENDDO
      ZSTAR(IX) = ZIN(IX)
      ZSTAR(IY) = ZIN(IY)
      IF(NDIM.EQ.3) ZSTAR(IZ) = ZIN(IZ)
C
      TEMP       = ZIN(IX)*ZIN(IX) + ZIN(IY)*ZIN(IY)
      IF( NDIM .EQ. 3 )TEMP = TEMP + ZIN(IZ)*ZIN(IZ)
      TEMP = HALF * TEMP
C
      ZRHO = ZERO
      DO ISP = 1,NSP
         ZRHO = ZRHO + ZSTAR(ISP)
      ENDDO 
C
C     compute the formation enthalpy
C 
        HFTOT = ZERO
         DO ISP = 1,NSP
             HFTOT = HFTOT + ZIN(ISP) * HF0(ISP)
         ENDDO
         HFTOT = HFTOT / HREFP
C
      ZSTAR(IE) = (GOGM1*POUTLET+TEMP)/ZRHO + HFTOT
C
      RETURN
      END
C
      subroutine ghost8vii4Ar(zin,zstar,work,ndim)
      IMPLICIT NONE
C
C     ghost state for subsonic inflow b.c. (plasma flow):
C     +) Mach number is taken from the interior
C     +) total pressure, total temperature, mass concentrations and flow angles are imposed
C
      INCLUDE 'plasma.h'    
      INCLUDE 'chem.h'    
      INCLUDE 'commonv.inc' 
C
      integer ndim,nofvar,isp
      parameter(nofvar=NSP+4)
      double precision zin(*),zstar(*),work(*)
      double precision uavg(NOFVAR)
C
      INCLUDE 'constants.h'
      INCLUDE 'pfcgas.com'
      INCLUDE 'ibc8.com'
      INCLUDE 'dofs.com'
      INCLUDE 'streamplasma.com'
      INCLUDE 'ioplasma.com'
      integer nerr,iopt
      double precision asqr,abar,mach,machsqr,kinetic,pres,dens,
     +temp,dens0,zrho,zrstar,hftot
      character*67 errmsg
C
C     work(1:*):
C     p0/p_ref,t0/t_ref,unused,n_x,n_y,n_z
C
C     compute Zrho
C
      ZRHO=ZERO
      DO ISP=1,NSP
         ZRHO=ZRHO + ZIN(ISP)
      ENDDO
C
C     compute mixture formation enthalpy
C        
      HFTOT = ZERO
      DO ISP = 1,NSP
          HFTOT = HFTOT + ZIN(ISP) * HF0(ISP)
      ENDDO
      HFTOT = HFTOT / ZRHO / HREFP
C
C     compute in the boundary node
C
      IF(NDIM.EQ.3)
     &UAVG(IZ) = ZIN(IZ)/ZRHO ! z component of the velocity vector
      UAVG(IY) = ZIN(IY)/ZRHO ! y comiponent of the velocity vector
      UAVG(IX) = ZIN(IX)/ZRHO ! x component of the velocity vector
      UAVG(IE) = ZIN(IE)/ZRHO ! Total Enthalpy
      DO ISP=1,NSP
         UAVG(ISP) = ZIN(ISP)*ZRHO ! Single species density
      ENDDO
C
      KINETIC = UAVG(IX)*UAVG(IX) + UAVG(IY)*UAVG(IY)
      IF( NDIM .EQ. 3 )KINETIC = KINETIC + UAVG(IZ)*UAVG(IZ)
      KINETIC = HALF * KINETIC
      ASQR = GM1 * ( UAVG(IE) - KINETIC - HFTOT )
      IF( ASQR .LT. ZERO )THEN
         WRITE(ERRMSG(1:67),FMT=333)ASQR,-10
  333 FORMAT('GHOST8VII4Ar',1X,'Negative averaged sound speed ',F7.3,
     &       ' in element ',I8)
         NERR = 6
         IOPT = 1
         CALL SETERR(ERRMSG(1:67),67,NERR,IOPT)
      ENDIF
      ABAR      = SQRT(ASQR)
      MACHSQR   = TWO * KINETIC / ASQR
      MACH      = SQRT( MACHSQR )
C   
C
      TEMP = (ONE+HALF*GM1*MACHSQR) 
C
C     rho/rho^0 = (1+0.5*gm1*mach**2)**expn * p^0/T^0
C
C     work(2) = p^0
C     work(3) = T^0
C
C     DENS0 = WORK(2)/WORK(3) !for perfect gas RSTAR=1
      DENS0 = WORK(2)/WORK(3)/RMIXSTAR  !RMIXSTAR depends on mixture
C     
      DENS = ( TEMP**(-ONE/GM1) ) * DENS0
      PRES = ( ONE/TEMP**GOGM1 ) * WORK(2)
      ABAR = SQRT(GAM*PRES/DENS)
C
      ZRSTAR = SQRT(DENS)
C
Crpepe    Freestreem concentration are stored in ALPHA1
Crpepe    ALPHA1 = rho_i1/rho_1
C
      DO ISP=1,NSP
         ZSTAR(ISP) = ALPHA1(ISP)*ZRSTAR
      ENDDO
      TEMP = ZRSTAR * ABAR * MACH
      ZSTAR(IE) = ZRSTAR * (GOGM1 * RMIXSTAR * WORK(3) + HFTOT) !rpepe verificare equazioni
      ZSTAR(IX) = TEMP * WORK(4)
      ZSTAR(IY) = TEMP * WORK(5)
      IF(NDIM.EQ.3)ZSTAR(IZ) = TEMP * WORK(6)
C
      RETURN
      END

