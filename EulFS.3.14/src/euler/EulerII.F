!> \copydetails EulerIX()
      SUBROUTINE EulerII(IELEM,VCN,VCB,VCZ,NDIM,NOFVERT,NOFVAR,
     &NTURB,NODRES,TSTEP,STIFEL,VOLUME,PICARD,ScalarScheme,MatrixScheme)
C
C     .. Hyperbolic Elliptic splitting using the vanLeer, Lee, Roe
C     preconditioning matrix ..
C     EulerII is ONLY valid for supersonic 2D flows, where
C     the preconditioned equations can be reduced to 4 scalar
C     equations describing the advection of entropy, total enthalpy
C     and two acoustic waves
C
C     $Id: EulerII.F,v 1.26 2020/03/28 09:51:15 abonfi Exp $
C
      IMPLICIT NONE
C
C
      INCLUDE 'paramt.h'
      INCLUDE 'bnd.h'
      INCLUDE 'constants.h'
C
C	.. NEQMAX is the max. no. of equations
C	.. NORDER is actual no. of equations = DIM
C
      INTEGER NEQMAX,IWRK,NORDER,FRSTEQ
      PARAMETER(NEQMAX=4,IWRK=10,NORDER=4,FRSTEQ=1)
      INTEGER MOVE(IWRK)
      DOUBLE PRECISION TOLER
      PARAMETER(TOLER=1.d-14) 
C
      INCLUDE 'blkEulerII.com'
      INCLUDE 'three.com'
      INCLUDE 'time.com'
      INCLUDE 'transf.com'
      INCLUDE 'flags.com'
C
      INTEGER IELEM,NDIM,NOFVERT,NOFVAR,i1mach
C
      EXTERNAL scalarScheme,MatrixScheme
C
      INTEGER KWAVE,JCOL,IDIM,ifail,I,J,N,M,IDX,N4
      INTEGER IOFF,IADD,JADD,IVERT,JVERT,NTURB,IVAR
      DOUBLE PRECISION FLUCT,SCALRES,HELP
      LOGICAL PICARD
C
      DOUBLE PRECISION JacobianMatrix(NEQMAX,NEQMAX*3),
     2PHI(NEQMAX),CHARV(12),DCHARV(12),
     3WKSP(MAXNOFVAR),SOURCE(MAXNOFVAR),WORK(MAX_NOFVERT_SQR),
     4TAUX(MAXNOFVAR*MAXNOFVERT)
      DOUBLE PRECISION NODRES(*)
      DOUBLE PRECISION TSTEP(*)
      DOUBLE PRECISION VCZ(*),VCN(*),VCB(*),
     +VOLUME,STIFC(144),TEMPA(400),TEMPB(100),STIFEL(*)
C
      DOUBLE PRECISION DDOT
      EXTERNAL         DDOT
C
      DATA SOURCE/MAXNOFVAR*ZERO/
C
      IDX(I,J,N,M) = (((J-1)*M+I-1)*N*N)+1
C
      IF( NDIM .NE. 2 )THEN
          STOP 'Only works in 2 space dimensions'
      ENDIF 
      IF( MACHSQR .LT. ONE )THEN
          WRITE(6,*)' Subsonic Mach number ',MACH,' in cell ',IELEM
          STOP
      ENDIF
      IF( NOFVAR .NE. (NDIM+2) )THEN
          WRITE(6,*)' Model NOT suited for turbulent calculations'
          STOP
      ENDIF
C
      N4 = NOFVAR*NOFVAR*NOFVERT*NOFVERT
C
C     The element stiffness matrix is initialized to 0.d0
C
      IF (PICARD) THEN
          CALL DINIT(N4,ZERO,STIFEL,1)
          CALL DINIT(N4,ZERO,STIFC,1)
      ENDIF
C
      CALL DINIT(NOFVERT*NOFVAR,ZERO,DCHARV,1)
C
      CALL DINIT(NORDER*NOFVERT,ZERO,TAUX,1)
C
C
C	Sets a stream aligned frame ..
C
      CALL StreamAlignedFrame(NDIM)
C
C       The Jacobian Matrix of the subsystem is assembled and
C          the eigenvectors computed ..
C
      CALL Eigen_II(JacobianMatrix,NEQMAX,DVDZ,DUDV,NDIM,NOFVAR)
C
C	---------- System scheme ----------
C
      IF(ICHECK.EQ.0)GOTO 7
C
C	Some initializations ....
C
      CALL DINIT(NORDER,ZERO,PHI,1)
      CALL DINIT(NORDER,ZERO,WKSP,1)
      CALL CHECK(IELEM,NDIM,NOFVAR)
c
c     COMPUTES THE RESIDUAL/VOLUME : dF/dU * dU/dX + dG/dU * dU/dy + ...
c
      DO 12 idim = 1 , NDIM
         jcol = (idim-1) * NEQMAX + 1
         CALL DGEMV( 'N' , NORDER , NORDER , ONE , JacobianMatrix
     +   (1,jcol) , NEQMAX , GRAD_CHAR(FrstEq,idim) ,  1, ONE , PHI ,
     +    1)
   12 CONTINUE
C
    7 CONTINUE
C
C	---------- Upwind scheme ----------
C
      CALL DGEMM('Transpose','Transpose',NOFVERT,NOFVAR,NOFVAR,
     +           ONE,VCZ,NOFVAR,dVdZ,NOFVAR,ZERO,CHARV,NOFVERT)
C
      DO 18 KWAVE = 1, NOFVAR
          IOFF = (KWAVE-1)*NOFVERT+1
          CALL DINIT(NOFVERT*NOFVERT,ZERO,STIFC,1)
C
C     Each advection eqn. is solved using an upwind scheme
C
          CALL ScalarScheme( IELEM, VCN, R_SPEED(1,KWAVE), SCALRES,
     +       ZERO, CHARV(IOFF), TAUX(IOFF), DCHARV(IOFF), WORK, STIFC,
     +       NDIM, NOFVERT, PICARD )
C
C
          IF (PICARD) THEN
C
C     STIFC(i,j) is copied into STIFEL(KWAVE,KWAVE,i,j)
C
              IADD = KWAVE + (KWAVE-1)*NOFVAR
              CALL DCOPY(NOFVERT*NOFVERT,STIFC,1,STIFEL(IADD),
     +               NOFVAR*NOFVAR)
          ENDIF
C
C
         IF( ICHECK .NE. 0 )THEN
C
C	 Checks for the scalar residual ..
C
            FLUCT = DDOT(NDIM,R_SPEED(1,KWAVE),1,GRAD_CHAR(KWAVE,1),LDW)
C
            SCALRES = SCALRES/VOLUME
            CALL DAXPY(NOFVAR,SCALRES,dUdV((KWAVE-1)*NOFVAR+1),1,WKSP,1)
C
            IF( DABS( FLUCT - SCALRES ) .GT. TOLER)THEN
               WRITE(I1MACH(2),89999)IELEM,KWAVE
               WRITE(I1MACH(2),*)FLUCT,SCALRES
               PAUSE
            ENDIF
C
         ENDIF
C
   18 CONTINUE ! end loop over waves
C
C Transform into conserved variables
C
      CALL DGEMM('No Transpose','Transpose',NOFVAR,NOFVERT,NOFVAR,
     +           ONE,dUdV,NOFVAR,DCHARV(1),NOFVERT,ZERO,NODRES,NOFVAR) 
C
!     CALL R8Mat_Print('G',' ',NOFVERT,NORDER,TAUX,
!    +            NOFVERT,'Nodal timestep ',IFAIL)
C
C     compute the timestep: TAUX is transposed, see the table above for DSYMM
C
      DO IVERT = 1, NOFVERT ! loop over the vertices
         JADD = IVERT
         HELP = ZERO 
C
C     sum or maximum over the dofs
C
         DO IVAR = 1,NOFVAR ! loop over the NDIM+2 dofs
            IF( CHAR_TIMESTEPPING )THEN
               HELP = MAX(HELP,TAUX(JADD))
            ELSE
               HELP = HELP + TAUX(JADD)
            ENDIF
            JADD = JADD + NOFVERT
         ENDDO ! end loop over the dofs
         IADD = (IVERT-1)*NOFVAR
         DO IVAR = 1,NOFVAR
            JADD = IADD + IVAR
            TSTEP(JADD) = TSTEP(JADD) + HELP
         ENDDO
      ENDDO
C
      IF( ICHECK .NE. 0 )
     +CALL TEST( DivFlux , WKSP , 1.D-15, IELEM , NOFVAR )
C
C     --------------- If explicit, return now  ---------------
C
      IF (.NOT.PICARD) RETURN
C
C
C     compute the transformation matrices from
C     conserved to parameter variables in the vertices
C
          DO 33 IVERT = 1,NOFVERT
              JADD = (IVERT-1)*NOFVAR*NOFVAR+1
              CALL MatdZdU(VCZ((IVERT-1)*NOFVAR+1),dZdU(JADD),
     +                     NDIM,NOFVAR)
   33     CONTINUE
C
C     transform the element stiffness matrix into conserved variables
C
          CALL DGEMM('No Transpose','No Transpose',NOFVAR,
     +               NOFVAR*NOFVERT*NOFVERT,NOFVAR,MONE,dUdV,
     +               NOFVAR,STIFEL,NOFVAR,ZERO,TEMPA,NOFVAR)
          CALL DGEMM('No Transpose','No Transpose',NOFVAR,
     +               NOFVAR*NOFVERT,NOFVAR,TWO,dVdZ,NOFVAR,dZdU,
     +               NOFVAR,ZERO,TEMPB,NOFVAR)
C
corig             CALL DGEMM('No Transpose','No Transpose',NOFVAR,
corig+                       NOFVAR,NOFVAR,ONE,STIFEL(IADD),NOFVAR,
corig+                       dVdZ,NOFVAR,ZERO,TEMPB,NOFVAR)
C     sign changed
corig             CALL DGEMM('No Transpose','No Transpose',NOFVAR,
corig+                       NOFVAR,NOFVAR,MONE,dUdV,NOFVAR,TEMPB,NOFVAR,
corig+                       ZERO,TEMPA,NOFVAR)
C

              DO 35 JVERT = 1,NOFVERT
                  JADD = (JVERT-1)*NOFVAR*NOFVAR+1
                  DO 35 IVERT = 1,NOFVERT
                  IADD = IDX(IVERT,JVERT,NOFVAR,NOFVERT)
                  CALL DGEMM('No Transpose','No Transpose',NOFVAR,
     +                       NOFVAR,NOFVAR,ONE,TEMPA(IADD),NOFVAR,
     +                       TEMPB(JADD),NOFVAR,ZERO,STIFEL(IADD),
     +                       NOFVAR)
C
   35     CONTINUE
C
C
      RETURN
89999 FORMAT(5X,'Scalar residual in Element ',I6,' Wave # ',I1)
      END
