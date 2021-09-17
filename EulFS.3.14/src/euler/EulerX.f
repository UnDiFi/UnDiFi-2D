!> \copydetails EulerIX()
      SUBROUTINE EulerX(IELEM,VCN,VCB,VCZ,NDIM,NOFVERT,NOFVAR,
     &                    NTURB,NODRES,TSTEP,STIFEL,VOLUME,PICARD,
     &                    ScalarScheme,MatrixScheme)
C
      IMPLICIT NONE 
C
C     Merkle's Preconditioned equations in primitive variables ..
C
C     $Id: EulerX.f,v 1.9 2020/03/28 09:51:15 abonfi Exp $
C
C
C
      INCLUDE 'paramt.h'
      INCLUDE 'constants.h'
      INCLUDE 'bnd.h'
      INCLUDE 'time.h'
C
C     NEQMAX is the max. no. of equations (5 in 3D)
C            for the matrix scheme (solves for dp,du,dv,dw,dT)
C     MAXNOFEQN is the max. no. of mean flow equations (5 in 3D)
C
      INTEGER NEQMAX,LNNVV
      DOUBLE PRECISION TOLER
      PARAMETER (NEQMAX=5,TOLER=1.D-15)
      PARAMETER (LNNVV=MAX_NOFVAR_SQR*MAX_NOFVERT_SQR)
      INTEGER IWRK
      PARAMETER(IWRK=10)
      INTEGER MOVE(IWRK)
C
C
      INCLUDE 'dofs.com'
      INCLUDE 'bodyf.com'
      INCLUDE 'flags.com'
      INCLUDE 'merkle.com'
      INCLUDE 'three.com'
      INCLUDE 'time.com'
      INCLUDE 'transf.com'
      INCLUDE 'visco.com'
C
C     .. Scalar Arguments ..
      INTEGER IELEM,NDIM,NOFVERT,NOFVAR,NTURB
      DOUBLE PRECISION VOLUME(*)
      LOGICAL PICARD
C
C     .. Array Arguments ..
      DOUBLE PRECISION VCZ(*),VCN(*),STIFEL(*),NODRES(*),TSTEP(*),VCB(*)
C
C
C     .. External Arguments ..
      EXTERNAL ScalarScheme,MatrixScheme
C
C
C
C     .. Local Scalar ..
      INTEGER IVAR,IVERT,JVERT,IDIM,IADD,JADD,KADD,I,J
      INTEGER NORDER,ifail,M,N,MN,N4,JCOL
      DOUBLE PRECISION FLUCT,SCALRES,PEH,HELP,SUMK,AM2PGR,AM2VIS,DP,BASE
     &,HEIGHT,CUTOFF,ULOCAL
      LOGICAL LFLAG
      LOGICAL unitmat
C
C
C
C     NODRES(1:NOFVAR,1:NOFVERT) is used to accumulate the
C         nodal residual in conserved variables and scatter
C         it to the RHS PETSc vector
C
C     TSTEP(1:NOFVERT) is used to accumulate the timestep
C         and then scatter it to the DT PETSc vector 
C
C     PRIMV(1:NOFEQN,1:NOFVERT) is used to store the vector
C         of symmetrizing variables
C
C     DPRIMV(1:NOFEQN,1:NOFVERT) is used to store the change
C         in the vector of symmetrizing variables
C
      DOUBLE PRECISION PRIMV(MAXNOFEQN*MAXNOFVERT),PHI(NEQMAX)
      DOUBLE PRECISION DPRIMV(MAXNOFEQN*MAXNOFVERT),WKSP(MAXNOFEQN),
     1                 TAUX(MAXNOFEQN*MAXNOFVERT),SOURCE(MAXNOFVAR),
     2                 Jacobian(NEQMAX,NEQMAX*3),RESIDUAL(2*NEQMAX),
     3                 TEMPA(MAX_NOFVERT_SQR*MAXNOFEQN**2),
     4                 TEMPB(MAX_NOFVERT_SQR*MAXNOFEQN**2)
      DOUBLE PRECISION DVDU(NEQMAX*NEQMAX),
     6                 BETA(MAXNOFEQN*MAXNOFEQN*MAXNOFVERT),
     8                 STIFC(MAX_NOFVERT_SQR*NEQMAX*NEQMAX)
C
C     NOFEQN (= DIM+2) is actual no. of mean flow equations
C     NORDER (= DIM+2) is actual no. of equations being solved
C                      with the system scheme
C
      INTEGER NOFEQN
C
C     RESIDUAL[1:NORDER] stores the residual computed by
C                        the Matrix scheme as \sum K_j U_j
C     RESIDUAL[NORDER+1:2*NORDER]
C                        stores the residual computed by
C                        the Matrix scheme as \sum C_{ij} U_j
C     it is used just for debugging purposes, to be compared with
C     the residual computed as:
C     dF/dU * dU/dX + dG/dU * dU/dy + [ dH/dU * dU/dz ]
C
C
      DOUBLE PRECISION DDOT,DNRM2
      EXTERNAL DDOT,DNRM2
C
      EXTERNAL MatSplitNum,MatSplitX
C
      DATA SOURCE/MAXNOFVAR*ZERO/
C
C
C     Statement function
C
      INTEGER JADDR
      JADDR(IVERT,JVERT,N) = (((JVERT-1)*NOFVERT+IVERT-1)*N*N) + 1
C
      NOFEQN = NDIM + 2
      NORDER = NOFEQN
      N4 = NOFVAR*NOFVAR*NOFVERT*NOFVERT
C
C     Here we set the preconditioner Mach number M_p
C
C     the length scale h is set equal to the distance btw the inflow
C     and outflow points based on the velocity vector
C
C     h = (Vol * |u|)/(sum k_j^+)
C
      SUMK = ZERO
      DO IVERT = 1, NOFVERT
         HELP = DDOT(NDIM,VCN((IVERT-1)*NDIM+1),1,UAVG(IX),1)
         SUMK = SUMK + MAX(ZERO,HELP)
      ENDDO
      HELP = REAL(NDIM)
      SUMK = SUMK/HELP ! SUMK = (u \cdot n_j )/d
      ULOCAL = DNRM2(NDIM,UAVG(IX),1)
      HEIGHT = VOLUME(1) * ULOCAL / SUMK
      BASE = HELP*VOLUME(1)/HEIGHT 
      IF(NDIM.EQ.3)BASE = SQRT(4.d0*BASE/3.141593)
C
C     here we assume that the local (non-dimensional) kinematic viscosity is 1
C     or, equivalently, that \nu = 1/Re
C     hence:
C     Pe_h = (u * h)/\nu =  
C
      PEH = MIN(BASE,HEIGHT) * RE * UAVG(1) * ULOCAL / ONE ! 1.d0 should be replaced by viscosity
C
C     Compute the pressure gradient
C
!     DP = ZERO
!     DO I = 1,NDIM
!        WKSP(I) = ZAVG(2)*GRAD_PRIM(1,I)+ZAVG(1)*GRAD_PRIM(2,I)
!        J = I+2
!        WKSP(I) = WKSP(I)-ZAVG(J)*GRAD_PARM(J,I)
!        WKSP(I) = WKSP(I) * GM1OG
!        DP = DP + UAVG(J) * WKSP(I)
!     ENDDO
C
C     Compute the pressure difference along the velocity
C
!     DP = ABS(DP) * VOLUME(1) / SUMK
!     AM2PGR = DP/(UAVG(1)*ASQR)
      AM2PGR = ZERO
      AM2VIS = MACHSQR/(PEH*PEH)
!     AM2VIS = ZERO
      CUTOFF = MERKLE_CUTOFF**2
      AMPSQR = MIN(MAX(MACHSQR,AM2PGR,AM2VIS,CUTOFF),ONE)
!     IF(MAX(MACHSQR,AM2PGR,AM2VIS).LT.CUTOFF)THEN
!     write(6,*)IELEM,MAX(MACHSQR,AM2PGR,AM2VIS),CUTOFF
!     ENDIF
!     write(22,*)IELEM,DSQRT(MACHSQR),DSQRT(AM2PGR),DSQRT(AM2VIS),
!    &DSQRT(AMPSQR),PEH,height,base
C
C     The element stiffness matrix is initialized to 0.d0
C
      IF (PICARD) THEN
          CALL DINIT(N4,ZERO,STIFEL,1)
          CALL DINIT((NORDER*NOFVERT)**2,ZERO,STIFC,1)
      ENDIF
      IF (LTIME) THEN
         DTVOL = DELT/VOLUME(1)
         IF(DUALTS)THEN 
            CALL DINIT(NOFEQN*NOFEQN*NOFVERT,ZERO,BETA,1)
         ENDIF
      ENDIF
C
C     set local residual and timestep to zero (should maybe bring it
C     in the calling routine)
C
      CALL DINIT(NOFVERT*NOFEQN,ZERO,DPRIMV,1)
      CALL DINIT(NOFVERT*NORDER,ZERO,TAUX,1)
C
C --------------- Debugging code starts here ---------------
C
      IF (ICHECK.NE.0) THEN
C
C    Some initializations ....
C
          CALL CHECK(IELEM,NDIM,NOFEQN)
      ENDIF
C
C --------------- Debugging code ends here ---------------
C
C
C     The Jacobian Matrix of the subsystem is assembled and
C         the eigenvectors computed ..
C
      CALL Eigen_X(Jacobian,NEQMAX,dVdZ,dUdV,dVdU,NDIM,NOFEQN)
C
      IF (ICHECK.EQ.0) GOTO 7
C
C --------------- Debugging code starts here ---------------
C
C     Here we compute first
C
C     dPdx dPdy dPdz
C     dTdx dTdy dTdz
C     dudx dudy dudz
C     dvdx dvdy dvdz
C     dwdx dwdy dwdz
C
C
      CALL DGEMM('No Transpose','No Transpose',NOFEQN,NDIM,NOFEQN,
     +           ONE,dVdZ,NOFEQN,GRAD_PARM,MAXNOFVAR,ZERO,TEMPA,NOFEQN)
C
!     CALL R8Mat_Print('G',' ',Nofeqn,Ndim,TEMPA,Nofeqn,
!    +      'd(P,u,v,T) array ',IFAIL)
C
C     COMPUTES THE RESIDUAL/VOLUME as:
C     dF/dU * dU/dX + dG/dU * dU/dy + [ dH/dU * dU/dz ]
C     for debugging purposes
C
      CALL DINIT(NORDER,ZERO,PHI,1)
      DO 12 idim = 1,NDIM
          JCOL = (idim-1)*NEQMAX + 1
          CALL DGEMV('N',NORDER,NORDER,ONE,Jacobian(1,JCOL),NEQMAX,
     +               TEMPA((idim-1)*NOFEQN+1),1,ONE,PHI,1)

   12 CONTINUE
C
C --------------- Debugging code ends here ---------------
C
    7 CONTINUE
C
C     The quasi linear form we discretize is in
C     primitive variables PRIMV = (dp,du,dv,dw,dT)
C
C
C     primitive variables are computed through the relation
C     \tilde{V_p}_i = M Z_j
C     where M = \frac{\partial V_p}{\partial Z}
C
C     dp(1) dp(1) dp(1) dp(1)
C     dT(1) dT(2) dT(3) dT(4)
C     du(1)   ....    ..    ..    ..
C     dv(1)   ....    ..    ..    ..
C     dw(1)   ....    ..    ..    ..
C
      CALL DGEMM('No Transpose','No Transpose',NOFEQN,NOFVERT,NOFEQN,
     +           ONE,dVdZ,NOFEQN,VCZ,NOFVAR,ZERO,PRIMV,NOFEQN)
C
!     CALL R8Mat_Print('G',' ',Nofeqn,NOFeqn,dVdZ,Nofeqn,
!    +      'dVdZ array ',IFAIL)
!     CALL R8Mat_Print('G',' ',Nofvar,NOFvert,vcz,Nofvar,
!    +      'Z variables ',IFAIL)
!     CALL R8Mat_Print('G',' ',NOFEQN,Nofvert,PRIMV,NofEQN,
!    +      '(p,u,v,T) variables ',IFAIL)
!     pause
C
C     ---------- Matrix scheme ----------
C
      SOURCE(2) = -GRAV(1)*VOLUME(1)
      SOURCE(3) = -GRAV(2)*VOLUME(1)
      IF(NDIM.EQ.3)SOURCE(4) = -GRAV(3)*VOLUME(1)
C
C REM: calling with TAUX(1) will add new contributions
C      to those already computed when solving entropy
C
      CALL MatrixScheme(MatSplitX,PRIMV(1),DPRIMV(1),TAUX(1),
     +                  BETA,STIFEL,NORDER,NORDER,NOFVERT,VCN,NDIM,
     +                  Jacobian,NEQMAX,RESIDUAL,SOURCE,IELEM,PICARD)
!     write(6,*)ielem
!     CALL R8Mat_Print('G',' ',Nofeqn,Nofvert,dPRIMV,Nofeqn,
!    +      'residual in prim. vars. ',IFAIL)
!     CALL R8Mat_Print('G',' ',NOrder,NOrder*nofvert,BETA,
!    +      NOrder,'Beta in symm only last d+1 eqns',IFAIL)
C
C
C     dp(1) dp(2) dp(3) dp(4)
C     dT(1) dT(2) dT(3) dT(4)
C     du(1)   ....    ..    ..    ..
C     dv(1)   ....    ..    ..    ..
C     dw(1)   ....    ..    ..    ..
C
C
C     copy the timestep from TAUX into TSTEP (update)
C
C
C     compute the timestep
C
      DO IVERT = 1, NOFVERT ! loop over the vertices
         IADD = (IVERT-1)*NORDER
         HELP = ZERO 
C
C     sum over the dofs
C
         DO IVAR = 1,NORDER
            JADD = IADD + IVAR
            IF( CHAR_TIMESTEPPING )THEN
               HELP = MAX(HELP,TAUX(JADD))
            ELSE
               HELP = HELP + TAUX(JADD)
            ENDIF
         ENDDO
         IADD = (IVERT-1)*NOFVAR
         DO IVAR = 1,NORDER
            JADD = IADD + IVAR
            TSTEP(JADD) = TSTEP(JADD) + HELP
         ENDDO
      ENDDO
C
C Transform the nodal residual into conserved variables
C
      CALL DGEMM('No Transpose','No Transpose',NOFEQN,NOFVERT,NOFEQN,
     +           ONE,dUdV,NOFEQN,DPRIMV,NOFEQN,ZERO,NODRES,
     +           NOFVAR)
C
!     CALL R8Mat_Print('G',' ',NOFEQN,NOFEQN,dUdV,nofeqn,
!    +      'dUdV variables ',IFAIL)
!     CALL R8Mat_Print('G',' ',NOFvar,NOFvert,NODRES,nofvar,
!    +      'residual (before) ',IFAIL)
C     pause
C
C
!     CALL R8Mat_Print('G',' ',NOFEQN,NOFEQN*nofvert,BETA,
!    +      NOFEQN,'Beta after insertion',IFAIL)
!     pause
C
C --------------- Debugging code starts here ---------------
C
C    Checks the decomposition ..
C
      IF (ICHECK.NE.0) THEN
          CALL DSCAL(2*NORDER,ONE/VOLUME(1),RESIDUAL,1)
          LFLAG = .TRUE.
          DO 18 IVAR = 1,NORDER
              IF (DABS(PHI(IVAR)-RESIDUAL(IVAR)).GT.
     +            TOLER) LFLAG = .FALSE.
   18     CONTINUE
          IF (LFLAG .EQV. .FALSE.) THEN
              WRITE (6,99999) IELEM
              DO 22 IVAR = 1,NORDER
                  WRITE (6,*) PHI(IVAR),RESIDUAL(IVAR),
     +              DABS(PHI(IVAR)-RESIDUAL(IVAR))
   22         CONTINUE
C           PAUSE
C
          ENDIF
C
C     transforms the residual into conserved variables
C
          CALL DGEMV('N',NOFEQN,NORDER,ONE,dUdV,
     +               NOFEQN,RESIDUAL,1,ZERO,WKSP,1)
C
C     test the residual as computed by the "explicit" scheme
C     WKSP := dUdV * \Phi
C
          CALL TEST(DivFlux,WKSP,TOLER,IELEM,NOFEQN)
      ENDIF
C
C --------------- Debugging code ends here ---------------
C
      IF (PICARD) THEN
C
C     compute the transformation matrices from
C     conserved to parameter variables in the vertices
C
          DO 33 IVERT = 1,NOFVERT
              IADD = (IVERT-1)*NOFVAR+1
              JADD = (IVERT-1)*NOFEQN*NOFEQN+1
c
c a bit of care here: Z and dZdU are dimensioned NOFVAR in
c subr MatdZdU(), but this should not be a problem if called
c with NOFEQN
c
          CALL MatdZdU(VCZ(IADD),dZdU(JADD),NDIM,NOFEQN)
   33     CONTINUE
      ENDIF
C
      IF (LTIME) THEN
          DO IVERT = 1,NOFVERT
              IADD = (IVERT-1)*NOFVAR+1
              JADD = (IVERT-1)*NOFEQN*NOFEQN+1
              CALL PARM2CONS(VCZ(IADD),TEMPB(JADD),NOFEQN,NDIM)
          ENDDO
      ENDIF
C
!     CALL R8Mat_Print('G',' ',NOFVAR,nofvert,NODRES,
!    +      NOFVAR,'Nodal update in U before td',IFAIL)
!
      IF(LTIME)THEN 
C
         CALL DGEMM('No Transpose','No Transpose',NOFEQN,
     +               NOFEQN*NOFVERT,NOFEQN,ONE,dUdV,NOFEQN,
     2               BETA,NOFEQN,ZERO,TEMPA,NOFEQN) 
         DO IVERT = 1,NOFVERT
            JADD = (IVERT-1)*NOFEQN*NOFEQN+1
            CALL DGEMM('No Transpose','No Transpose',NOFEQN,
     +               NOFEQN,NOFEQN,ONE,TEMPA(JADD),NOFEQN,DVDU,NOFEQN,
     2               ZERO,BETA(JADD),NOFEQN) 
         ENDDO
C
         CALL UNSTEADY4(TEMPB,BETA,VCZ,NOFVAR,NODRES,STIFEL,VOLUME,
     &                  NOFEQN,NDIM,NOFVERT,PICARD)
      ENDIF
C
!     write(6,*)'ielem = ',ielem
!     CALL R8Mat_Print('G',' ',Nofeqn,nofeqn*nofvert,beta,Nofeqn,
!    +      'beta(U) ',IFAIL)
!     CALL R8Mat_Print('G',' ',NOFEQN,nofvert,NODRES,
!    +      NOFEQN,'Nodal update in U',IFAIL)
C
C
C     --------------- If explicit, return now  ---------------
C
      IF (.NOT.PICARD) RETURN
C
C     transform the element stiffness matrix into conserved variables
C
      IF( NOFVAR .EQ. NOFEQN )THEN
          CALL DGEMM('No Transpose','No Transpose',NOFVAR,
     +               NOFVAR*NOFVERT*NOFVERT,NOFVAR,MONE,dUdV,
     +               NOFVAR,STIFEL,NOFVAR,ZERO,TEMPA,NOFVAR)
      ELSE
          DO 13 JVERT = 1, NOFVERT
             DO 13 IVERT = 1, NOFVERT
                JADD = JADDR(IVERT,JVERT,NOFEQN)
                IADD = JADDR(IVERT,JVERT,NOFVAR)
                CALL DGEMM('No Transpose','No Transpose',NOFEQN,
     +               NOFEQN,NOFEQN,MONE,dUdV,
     +               NOFEQN,STIFEL(IADD),NOFVAR,ZERO,TEMPA(JADD),NOFEQN)
!     write(6,*)ivert,jvert
!     CALL R8Mat_Print('G',' ',NOFEQN,NOFEQN,tempa(jadd),NOFEQN,
!    +      'TEMPA ',IFAIL)
!     pause
   13     CONTINUE
      ENDIF
C
C
      CALL DGEMM('No Transpose','No Transpose',NOFEQN,
     +               NOFEQN*NOFVERT,NOFEQN,TWO,dVdZ,NOFEQN,dZdU,
     +               NOFEQN,ZERO,TEMPB,NOFEQN)
!     CALL R8Mat_Print('G',' ',NOFEQN,NOFEQN*nofvert,tempb,NOFEQN,
!    +      'TEMPB ',IFAIL)
!     pause
C
          DO 35 JVERT = 1,NOFVERT
                  JADD = (JVERT-1)*NOFEQN*NOFEQN + 1
          DO 35 IVERT = 1,NOFVERT
                  IADD = JADDR(IVERT,JVERT,NOFVAR)
                  KADD = JADDR(IVERT,JVERT,NOFEQN)
C
               CALL DGEMM('No Transpose','No Transpose',NOFEQN,
     +              NOFEQN,NOFEQN,ONE,TEMPA(KADD),NOFEQN,
     +              TEMPB(JADD),NOFEQN,ZERO,STIFEL(IADD),NOFVAR)
C
C     now STIFEL contains the convection stiffness matrix in
C     conserved variables
C
!     write(6,*)ivert,jvert
!     CALL R8Mat_Print('G',' ',Nofvar,Nofvar,stifel(iadd),Nofvar,
!    +      'C(i,j) ',IFAIL)
!
   35     CONTINUE
!     pause
C
C
      IF (ICHECK.EQ.0) RETURN
C
C --------------- Debugging code starts here ---------------
C
C     test the residual as computed by the "implicit" scheme
C     WKSP := dUdV(1,1) * \Phi
C
      CALL DGEMV('N',NOFEQN,NORDER,ONE,dUdV(1),
     +           NOFEQN,RESIDUAL(NORDER+1),1,ZERO,WKSP,1)
      CALL TEST(DivFlux,WKSP,TOLER,-IELEM,NOFEQN)
C
C --------------- Debugging code ends here ---------------
C
      RETURN

99999 FORMAT (5X,'Vector residual in Element ',I6,' EulerX')

      END
