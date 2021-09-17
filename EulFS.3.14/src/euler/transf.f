      SUBROUTINE parm_to_cons(ZROE,NDIM,NOFVAR,NPOIN,LFLAG,IFAIL)
C
C     $Id: transf.f,v 1.11 2013/04/30 07:42:54 abonfi Exp $
C
      IMPLICIT NONE
C
      INCLUDE 'constants.h'
      INCLUDE 'pfcgas.com'
C
C     .. Array Arguments ..
C
      INTEGER NDIM,NOFVAR,NPOIN,IFAIL
      LOGICAL LFLAG
      DOUBLE PRECISION ZROE(NOFVAR,NPOIN)
C
C     .. Local Scalars ..
C
      INTEGER INODE 
      DOUBLE PRECISION SUM
C
C     .. External Functions ..
C
C
C     .. Executable Statements ..
C
      IFAIL = 0 
      do 100 INODE = 1 , NPOIN
         IF(LFLAG)ZROE(NOFVAR,INODE)=ZROE(1,INODE) * ZROE(NOFVAR,INODE)
         SUM = ZROE(3,INODE)**2 + ZROE(4,INODE)**2
         IF( NDIM .EQ. 3 )SUM = SUM + ZROE(5,INODE)**2
         ZROE(2,INODE) = GINV * ZROE(1,INODE) * ZROE(2,INODE) +
     .   HALF * GM1OG * SUM
C
         IF(NDIM.EQ.3)ZROE(5,INODE) = ZROE(1,INODE) * ZROE(5,INODE)
         ZROE(4,INODE) = ZROE(1,INODE) * ZROE(4,INODE)
         ZROE(3,INODE) = ZROE(1,INODE) * ZROE(3,INODE)
         ZROE(1,INODE) = ZROE(1,INODE) * ZROE(1,INODE)
  100 CONTINUE
c
      RETURN
      END

      SUBROUTINE cons_to_parm(ZROE,NDIM,NOFVAR,NPOIN,LFLAG,IFAIL)
C
      IMPLICIT NONE
C
c
c This routine transforms the conservative vector into the
c paramter vector overwriting the array Z.
c
      INCLUDE'constants.h'
      INCLUDE 'pfcgas.com'
C
C     .. Array Arguments ..
C
      INTEGER NDIM,NOFVAR,NPOIN,IFAIL
      LOGICAL LFLAG
      DOUBLE PRECISION ZROE(NOFVAR,NPOIN)
C
C     .. Local Scalars ..
C
      INTEGER INODE
      DOUBLE PRECISION SUM
C
C     .. External Functions ..
C
C
C     .. Executable Statements ..
C
      IFAIL = 0 
      DO 100 INODE = 1 , NPOIN
C
         IF(ZROE(1,INODE) .LE. ZERO )THEN
           IFAIL = INODE
           RETURN
         ENDIF
C
         ZROE(1,INODE) = SQRT(ZROE(1,INODE))
         ZROE(3,INODE) = ZROE(3,INODE) / ZROE(1,INODE)
         ZROE(4,INODE) = ZROE(4,INODE) / ZROE(1,INODE)
         IF(NDIM .EQ. 3 )ZROE(5,INODE) = ZROE(5,INODE) / ZROE(1,INODE)
         SUM = ZROE(3,INODE)**2 + ZROE(4,INODE)**2
         IF(NDIM .EQ. 3 )SUM = SUM + ZROE(5,INODE)**2
         ZROE(2,INODE) = GAM *( ZROE(2,INODE) - HALF * GM1OG * SUM )
     & / ZROE(1,INODE)
         IF(LFLAG)ZROE(NOFVAR,INODE)=ZROE(NOFVAR,INODE)/ ZROE(1,INODE)
  100 CONTINUE
C
      RETURN
      END
C
C ------------------------------+------------------------------
C
      SUBROUTINE PARM2PRIM(NDIM,IELEM)
C
      IMPLICIT NONE
C
      INCLUDE 'paramt.h'
      INCLUDE 'constants.h'
      INCLUDE 'three.com'
      INCLUDE 'pfcgas.com'
C
C     .. Scalar Arguments ..
C
      INTEGER NDIM,IELEM
      INTEGER NERR,IOPT
      CHARACTER*72 ERRMSG
      integer k
C
C     .. Intrinsic Functions ..
C
      INTRINSIC SQRT
C
C     .. External Functions ..
C
C
C     .. Executable Statements ..
C
      IF(NDIM.EQ.3)
     &UAVG(5) = ZAVG(5)/ZAVG(1) ! z component of the velocity vector
      UAVG(4) = ZAVG(4)/ZAVG(1) ! y component of the velocity vector
      UAVG(3) = ZAVG(3)/ZAVG(1) ! x component of the velocity vector
      UAVG(2) = ZAVG(2)/ZAVG(1) ! Total Enthalpy
      UAVG(1) = ZAVG(1)*ZAVG(1) ! Density
C
      KINETIC = UAVG(3)*UAVG(3) + UAVG(4)*UAVG(4)
      IF( NDIM .EQ. 3 )KINETIC = KINETIC + UAVG(5)*UAVG(5) 
      KINETIC = HALF * KINETIC
      ASQR = GM1 * ( UAVG(2) - KINETIC )
      IF( ASQR .LT. ZERO )THEN
         WRITE(6,FMT=*)ASQR,IELEM
         WRITE(ERRMSG(1:67),FMT=333)ASQR,IELEM
  333 FORMAT('PARM2PRIM',1X,'Negative averaged sound speed ',E7.3,
     &       ' in element ',I8)
         WRITE(6,FMT=*)'d = ',ndim
         WRITE(6,FMT=*)'k = ',kinetic
         WRITE(6,FMT=*)(uavg(k),k=1,NMAX)
         NERR = 6
         IOPT = 1
         CALL SETERR(ERRMSG(1:67),67,NERR,IOPT)
      ENDIF
      ABAR      = SQRT(ASQR)
c     MACHSQR   = DMAX1( TWO * KINETIC / ASQR, EPSQR_STAG) ! Barth & Deconinck
      MACHSQR   = TWO * KINETIC / ASQR
      MACH      = SQRT( MACHSQR )

      QINV      = ONE/(MACH*ABAR)
C
      RETURN
      END
C
C
C ------------------------------+------------------------------
C
      SUBROUTINE PARM2PRIM4Ar(NDIM,IELEM)
C
      IMPLICIT NONE
C
      INCLUDE 'paramt.h'
      INCLUDE 'constants.h'
      INCLUDE 'plasma.h'
      INCLUDE 'chem.h'
      INCLUDE 'commonv.inc'
      INCLUDE 'pfcgas.com'
      INCLUDE 'dofs.com'
      INCLUDE 'streamplasma.com'
      INCLUDE 'three.com'
      INCLUDE 'four.com'
C
C     .. Scalar Arguments ..
      INTEGER NDIM,IELEM
      INTEGER NERR,IOPT,ISP
      CHARACTER*72 ERRMSG
      DOUBLE PRECISION DPDR
      INTEGER k,I
C
C     .. Intrinsic Functions ..
C
      INTRINSIC SQRT
C
C     .. External Functions ..
C
      DOUBLE PRECISION PIM,PIE,PIR
C
C     .. Executable Statements ..
C
C
      SQRTR = ZERO
      DO ISP = 1,NSP
         SQRTR = SQRTR + ZAVG(ISP) ! square root of total density
      ENDDO
      IF(NDIM.EQ.3)
     &UAVG(IZ) = ZAVG(IZ)/SQRTR ! z component of the velocity vector
      UAVG(IY) = ZAVG(IY)/SQRTR ! y component of the velocity vector
      UAVG(IX) = ZAVG(IX)/SQRTR ! x component of the velocity vector
      UAVG(IE) = ZAVG(IE)/SQRTR ! Total Enthalpy
      DO ISP = NSP,1,-1
         UAVG(ISP) = ZAVG(ISP)*SQRTR ! Density of the NSP species
      ENDDO
C
      KINETIC = UAVG(IX)*UAVG(IX) + UAVG(IY)*UAVG(IY)
      IF( NDIM .EQ. 3 )KINETIC = KINETIC + UAVG(IZ)*UAVG(IZ) 
      KINETIC = HALF * KINETIC
C   
C      Pressure derivatives ..      
C
!rpepe      CALL PRESSDERIV(ZAVG,CHI,KAPPA)
!      write(6,*)'k=',KAPPA     
      
      KAPPA = GM1
C
      DENS = SQRTR*SQRTR
      DENSINV = ONE/DENS
C
!       write(6,*)'chis=',CHI
!       write(6,*)'KINETIC=',KINETIC
!       write(6,*)'UAVG=',UAVG     
C  
      DO ISP = 1,NSP
         CHI(ISP) = -GM1*HF0(ISP)/HREFP
!         write(6,*)'chi(',isp,') --> ',chi(isp)         
         ALPHA(ISP) = UAVG(ISP) * DENSINV
         DR(ISP) = PIR(UAVG(IX),NDIM,CHI(ISP),KAPPA)
!         write(6,*)'dr(',isp,') --> ',dr(isp)
!         DR(ISP) = CHI(ISP) + KAPPA*KINETIC
!         write(6,*)'dr2(',isp,') --> ',dr(isp)          
!         write(6,*)'ALPHA(',isp,') --> ',ALPHA(isp)
      ENDDO    
C
      DE = PIE(KAPPA)
!      write(6,*)'de --> ',DE
C
      DO 14 I = 1 , NDIM
         DM(I) = PIM(UAVG(IE+I),KAPPA)         
!         write(6,*)'dm(',I,') --> ',DM(I)
!         DM(I) = -KAPPA*UAVG(IE+I)
!         write(6,*)'dm2(',I,') --> ',DM(I)
   14 CONTINUE
C
C     Sound speed ..
C      
      DPDR = ZERO
      DO ISP = 1 , NSP
        DPDR = DPDR + ALPHA(ISP)*DR(ISP) 
      ENDDO       
      ASQR = DPDR + DE * (UAVG(IE) - TWO * KINETIC) ! a^2 = \sum_{i=1}^{N_s} \alpha_i \Pi_{\rho_i}  + \Pi_{\rho E}  (H - \mathbf{u} \cdot \mathbf{u})       
!       write(*,*)'ASQR pressure derivatives=',ASQR
!       ASQR = GM1 * ( UAVG(IE) - KINETIC ) 
!      write(*,*)'ASQR Z=',ASQR      
!        pause
      IF( ASQR .LT. ZERO )THEN
         WRITE(6,FMT=*)ASQR,IELEM
         WRITE(ERRMSG(1:70),FMT=333)ASQR,IELEM
  333 FORMAT('PARM2PRIM4Ar',1X,'Negative averaged sound speed ',E7.3,
     &       ' in element ',I8)
         WRITE(6,FMT=*)'d = ',ndim
         WRITE(6,FMT=*)'DPDR = ',DPDR
         WRITE(6,FMT=*)'k = ',kinetic
         WRITE(6,FMT=*)(uavg(k),k=1,MAXNOFVAR)
         NERR = 6
         IOPT = 1
         CALL SETERR(ERRMSG(1:67),67,NERR,IOPT)
      ENDIF
      ABAR      = SQRT(ASQR)
c     MACHSQR   = DMAX1( TWO * KINETIC / ASQR, EPSQR_STAG) ! Barth & Deconinck
      MACHSQR   = TWO * KINETIC / ASQR
      MACH      = SQRT( MACHSQR )

      QINV      = ONE/(MACH*ABAR)
C
      RETURN
      END
C

      SUBROUTINE parm_to_cons4Ar(ZROE,NDIM,NOFVAR,NPOIN,LFLAG,IFAIL)
C
C     $Id: transf.f,v 1.11 2013/04/30 07:42:54 abonfi Exp $
C
      IMPLICIT NONE
C
      INCLUDE 'paramt.h'
      INCLUDE 'constants.h'
      INCLUDE 'plasma.h'
      INCLUDE 'chem.h'
      INCLUDE 'commonv.inc' 
      INCLUDE 'dofs.com'
      INCLUDE 'streamplasma.com'
      INCLUDE 'pfcgas.com'
C
C     .. Array Arguments ..
C
      INTEGER NDIM,NOFVAR,NPOIN,IFAIL
      LOGICAL LFLAG
      DOUBLE PRECISION ZROE(NOFVAR,NPOIN)
C
C     .. Local Scalars ..
C
      INTEGER INODE,ISP
      DOUBLE PRECISION KINE,SQRTR,HFTOT
C
C     .. External Functions ..
C
C
C     .. Executable Statements ..
C
      IFAIL = 0      
C
      do 100 INODE = 1 , NPOIN
         SQRTR = ZERO
         DO ISP = 1,NSP
            SQRTR = SQRTR + ZROE(ISP,INODE) ! \sqrt{\rho} = \sum_{s=1}^{NSP} \sqrt{\rho_s}
         ENDDO
         KINE = ZROE(IX,INODE)**2 + ZROE(IY,INODE)**2
         IF( NDIM .EQ. 3 )KINE = KINE + ZROE(IZ,INODE)**2
         KINE = HALF * KINE ! \rho (u^2+v^2+w^2)
C
C        FORMATION ENTHALPY
         HFTOT = ZERO
         DO ISP = 1,NSP
            HFTOT = HFTOT + ZROE(ISP,INODE) * HF0(ISP)
         ENDDO
         HFTOT = HFTOT * SQRTR / HREFP
C         write(6,*) "HFTOT=",HFTOT
C         write(6,*) "HREFP",HREFP
C         pause
C
C        pressure
!        DPDRZ = ZERO
!        DO ISP = 1,NSP
!           DPDRZ = DPDRZ + DR(ISP) * ZROE(ISP,INODE)
!        ENDDO
!        PRESS = DPDRZ + DE*ZROE(IE,INODE) + DM(1)*ZROE(IX,INODE)
!    &         + DM(2)*ZROE(IY,INODE)
!        IF (NDIM .EQ. 3) THEN
!           PRESS = PRESS + DM(3)*ZROE(IZ,INODE)  
!        ENDIF
!        PRESS = PRESS * SQRTR / (ONE + DE)
C         
!        ZROE(IE,INODE) = SQRTR * ZROE(IE,INODE) - PRESS    
!         ZROE(IE,INODE) = GINV*(SQRTR*ZROE(IE,INODE)+GM1*KINE) ! has become \rho E !rpepe manca Hf0
         ZROE(IE,INODE) = GINV*(SQRTR*ZROE(IE,INODE)+GM1*(KINE+HFTOT))
C
         IF(NDIM.EQ.3)ZROE(IZ,INODE) = SQRTR * ZROE(IZ,INODE) ! becomes \rho w
         ZROE(IY,INODE) = SQRTR * ZROE(IY,INODE) ! becomes \rho v
         ZROE(IX,INODE) = SQRTR * ZROE(IX,INODE) ! becomes \rho u
         DO ISP = IE-1,1,-1 ! loop over species
            ZROE(ISP,INODE) = SQRTR * ZROE(ISP,INODE) ! is now \rho_s
         ENDDO
  100 CONTINUE
c
      RETURN
      END

      SUBROUTINE cons_to_parm4Ar(ZROE,NDIM,NOFVAR,NPOIN,LFLAG,IFAIL)
C
      IMPLICIT NONE
C
c
c This routine transforms the conservative vector into the
c paramter vector overwriting the array Z.
c
      INCLUDE 'paramt.h' 
      INCLUDE 'constants.h'
      INCLUDE 'plasma.h'
      INCLUDE 'chem.h'
      INCLUDE 'commonv.inc'
      INCLUDE 'dofs.com'
      INCLUDE 'streamplasma.com'
      INCLUDE 'pfcgas.com'
C
C     .. Array Arguments ..
C
      INTEGER NDIM,NOFVAR,NPOIN,IFAIL
      LOGICAL LFLAG
      DOUBLE PRECISION ZROE(NOFVAR,NPOIN)
C
C     .. Local Scalars ..
C
      INTEGER INODE,ISP
      DOUBLE PRECISION KINE,SQRTR,RHO,PRESS,DPDRU,HFTOT
C
C     .. External Functions ..
C
C
C     .. Executable Statements ..
C
      IFAIL = 0 
      DO 100 INODE = 1 , NPOIN
C
         RHO = ZERO ! total density
         DO ISP = 1, NSP
            RHO = RHO + ZROE(ISP,INODE)
         ENDDO
         IF(RHO .LT. ZERO )THEN
            WRITE(6,*)'r = ',RHO,INODE,(ZROE(ISP,INODE),ISP=1,NOFVAR) 
            IFAIL = INODE
            RETURN
         ENDIF
C
         SQRTR = SQRT(RHO)
C
C        Formation Enthalpy
         HFTOT = ZERO
         DO ISP = 1,NSP
C             HFTOT = ZROE(ISP,INODE) * HF0(ISP)
             HFTOT = HFTOT + ZROE(ISP,INODE) * HF0(ISP)
         ENDDO
!         HFTOT = HFTOT * SQRTR / HREFP
         HFTOT = HFTOT/ HREFP
C
C        pressure
!        DPDRU = 0 
!        DO ISP = 1,NSP
!           DPDRU = DPDRU + DR(ISP) * ZROE(ISP,INODE)
!        ENDDO
!        PRESS = DPDRU + DE*ZROE(IE,INODE) + DM(1)*ZROE(IX,INODE)
!    &         + DM(2)*ZROE(IY,INODE)
!        IF (NDIM .EQ. 3) THEN
!           PRESS = PRESS + DM(3)*ZROE(IZ,INODE)  
!        ENDIF
C         
         ZROE(IX,INODE) = ZROE(IX,INODE) / SQRTR ! \rho u := \sqrt{\rho} u
         ZROE(IY,INODE) = ZROE(IY,INODE) / SQRTR ! \rho u := \sqrt{\rho} u
         IF(NDIM .EQ. 3 )ZROE(IZ,INODE) = ZROE(IZ,INODE) / SQRTR ! \rho w := \sqrt{\rho} w
         KINE = ZROE(IX,INODE)**2 + ZROE(IY,INODE)**2
         IF(NDIM .EQ. 3)KINE = KINE + ZROE(IZ,INODE)**2
         KINE = HALF * KINE  ! = \rho ( u^2 + v^2 + w^2 ) / 2
!        ZROE(IE,NODE) = (ZROE(IE,NODE) + PRESS) / SQRTR  !Da utilizzare per miscele di gas reali
!         ZROE(IE,INODE) = GAM *( ZROE(IE,INODE) - GM1OG * KINE ) / SQRTR !manca entalpia di formazione
         ZROE(IE,INODE) = GAM *(ZROE(IE,INODE) - GM1OG * (KINE 
     &                     + HFTOT)) / SQRTR   
         DO ISP = IE-1,1,-1
            ZROE(ISP,INODE) = ZROE(ISP,INODE) / SQRTR ! \rho_s := \rho_s / \sqrt{\rho_1}
         ENDDO
  100 CONTINUE
C
      RETURN
      END
