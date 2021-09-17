!> \details
!> @param[in] IELEM the current simplicial element
!> @param[in] VCN the NDIM cartesian components of the inward face normal, scaled by its measure in the NOFVERT vertices
!> @param[in] VCB the NDIM cartesian components of the grid velocity
!> @param[in] VCZ the NOFVAR dofs in the NOFVERT vertices
!> @param[in] NDIM the dimension of the space
!> @param[in] NOFVERT the nof vertices of the current simplex (=NDIM+1)
!> @param[in] NOFVAR the nof dependent variables
!> @param[in] NTURB the nof turbulent variables
!> @param[out] NODRES the NOFVAR components of the residual vector in the NOFVERT vertices
!> @param[out] TSTEP is the contribution of the current cell to the global array 
!> @param[out] STIFEL is the contribution of the current cell to the global jacobian matrix
!> @param[in] VOLUME is the area/volume of the current cell
!> @param[in] PICARD is .TRUE. if the (approximate) jacobian matrix has to be assembled
!> @param[in] ScalarScheme the subroutine implementing the FS scheme for scalar equations
!> @param[in] MatrixScheme the subroutine implementing the FS scheme for hyperbolic system
C
      SUBROUTINE EulerIX(IELEM,VCN,VCB,VCZ,NDIM,NOFVERT,NOFVAR,
     &                   NTURB,NODRES,TSTEP,STIFEL,VOLUME,PICARD,
     &                   ScalarScheme,MatrixScheme)
C
C     Hyperbolic Elliptic splitting for incompressible
C     flows using the Turkel preconditioning matrix ..
C
C     $Id: EulerIX.f,v 1.18 2013/05/03 09:55:22 abonfi Exp $
C
C
      IMPLICIT NONE
C
      INCLUDE 'paramt.h'
      INCLUDE 'constants.h'
      INCLUDE 'bnd.h'
C
C     MAXNORDER is the max. no. of equations (3 in 3D)
C     NORDER is actual no. of equations = DIM
C
      INTEGER   MAXNORDER
      PARAMETER(MAXNORDER=3)
C
      INCLUDE 'dofs.com'
      INCLUDE 'three.com'
      INCLUDE 'transf.com'
      INCLUDE 'flags.com'
      INCLUDE 'time.com'
C
      INTEGER  IELEM,NDIM,NOFVERT,NOFVAR
C
      DOUBLE PRECISION DNRM2
      EXTERNAL DNRM2
      EXTERNAL ScalarScheme,MatrixScheme
C
C
      INTEGER FrstEq,JCOL,IDIM,ifail,IWRK
      PARAMETER (IWRK=10,FRSTEQ=2)
      INTEGER NORDER,MOVE(IWRK)
      DOUBLE PRECISION SCALRES,HELP
      LOGICAL PICARD
C
      INTEGER IDX,I,J,N,M,N4,IVERT,JVERT,NTURB,IADDR,JADDR,NOFEQN,IVAR
      DOUBLE PRECISION VCZ(*),VCB(*),VCN(*),NODRES(*),
     &VOLUME,STIFEL(*),TSTEP(*)
      DOUBLE PRECISION STIFC((MAXNORDER*VMAX)**2),PHI(MAXNORDER),
     2CHARV(MAXNOFEQN*VMAX),Jacobian(MAXNORDER,MAXNORDER*3),
     3DCHARV(MAXNOFEQN*VMAX),RESIDUAL(2*MAXNORDER),WKSP(5),
     4TEMPA(MAX_NOFVAR_SQR*MAX_NOFVERT_SQR),TAUX(MAXNOFVAR*MAXNOFVERT),
     5SOURCE(MAXNOFVAR)
C
C
      EXTERNAL MatSplitNum,MatSplitIX
      DATA SOURCE/MAXNOFVAR*ZERO/
C
      IDX(I,J,N,M) = (((J-1)*M+I-1)*N*N)+1
C
      NOFEQN = NDIM+1
      NORDER = NDIM
      N4 = NOFVAR*NOFVAR*NOFVERT*NOFVERT
C
C     The element stiffness matrix is initialized to 0.d0
C
      IF (PICARD) THEN
          CALL DINIT(N4,ZERO,STIFEL,1)
          CALL DINIT((NORDER*NOFVERT)**2,ZERO,STIFC,1)
      ENDIF
C
C     set local residual and timestep to zero
C
      CALL DINIT(NOFVERT*NOFEQN,ZERO,DCHARV,1)
      CALL DINIT(NOFVERT*NOFEQN,ZERO,TAUX,1)
C
C     smthg I do not like very much, but.......for the time being
C
      UAVG(3)=ZAVG(2) ! why am I doing this?
      UAVG(4)=ZAVG(3)
      UAVG(5)=ZAVG(4)
      QINV = ONE/DNRM2(NDIM,ZAVG(IX),1)
C
C     Sets up a stream aligned frame ..
C
      CALL StreamAlignedFrame(NDIM)
C
C     The Jacobian Matrix of the subsystem is assembled and
C         the eigenvectors computed ..
C
      CALL Eigen_IX(Jacobian,MAXNORDER,DVDZ,DUDV,NDIM,NOFEQN)
C
C --------------- Debugging code ends here ---------------
C
      IF(ICHECK.EQ.0)GOTO 7
C
C     Some initializations ....
C
      CALL DINIT(NORDER,ZERO,PHI,1)
      CALL CHECK(IELEM,NDIM,NOFEQN)
c
c     COMPUTES THE RESIDUAL/VOLUME : dF/dU * dU/dX + dG/dU * dU/dy + ...
c
      DO 12 idim = 1 , NDIM
         jcol = (idim-1) * MAXNORDER + 1
         CALL DGEMV( 'N' , NORDER , NORDER , ONE , Jacobian
     +   (1,jcol) , MAXNORDER , GRAD_CHAR(FrstEq,idim) ,  1, ONE , PHI ,
     +    1)
   12 CONTINUE
C
C --------------- Debugging code ends here ---------------
C
    7 CONTINUE
C
C     Builds the nodal vector of characteristic variables (CHARV)
C     call these (H,P,Q,R)
C
      CALL DGEMM('Transpose','Transpose',NOFVERT,NOFEQN,NOFEQN,
     +           ONE,VCZ,NOFVAR,dVdZ,NOFEQN,ZERO,CHARV,NOFVERT)
C
C     The matrix CHARV now looks like:
C
C     dP(1)  qdr(1)  dp(1)  qds(1)
C      ...    ...    ...    ...
C      ...    ...    ...    ...
C           
C     dP(4)  qdr(4)  dp(4)  qds(4)
C
C     the matrix of the characteristic variables qdr,dp,qds
C     is now transposed 
C
      IADDR = NOFVERT+1 
      CALL TRANS(CHARV(IADDR),NOFVERT,NORDER,NOFVERT*NORDER,
     +           MOVE,IWRK,IFAIL)
      IF(IFAIL.NE.0)WRITE(6,*)'TRANS HAS RETURNED IFAIL = ',IFAIL
C
C     so that these are now stored 
C     (starting at CHARV(IADDR)=qdr(1)) as CHARV(1:NORDER,1:NOFVERT)
C
C      qdr(1)  ...   qdr(4)
C      dp(1)   ...   dp(4) 
C      qds(1)  ...   qds(4)
C
C ************************************************************
C The non commuting system is solved using a matrix Scheme
C ************************************************************
C
          CALL MatrixScheme(MatSplitIX,CHARV(IADDR),DCHARV(IADDR),
     +                      TAUX(NOFVERT+1),TEMPA,STIFC,NORDER,NORDER,
     +                      NOFVERT,VCN,NDIM,Jacobian,MAXNORDER,
     &                      RESIDUAL,SOURCE,IELEM,PICARD)
C
      IF (PICARD) THEN
          CALL MATINS(STIFEL,NOFVAR,STIFC,NORDER,NOFVERT,NOFVERT,1)
          CALL DINIT(NORDER*NORDER*NOFVERT*NOFVERT,ZERO,STIFC,1)
      ENDIF
C
C ************************************************************
C     Solve the total pressure transport eqn. 
C ************************************************************
C
      CALL ScalarScheme(IELEM,VCN,R_SPEED(1,1),SCALRES,ZERO,
     +                  CHARV(1),TAUX(1),DCHARV(1),TEMPA,STIFC,
     +                  NDIM,NOFVERT,PICARD)
C
C    Copy the convective jacobian (STIFC) into STIFEL and
C
      IF (PICARD) THEN
          CALL DCOPY(NOFVERT*NOFVERT,STIFC,1,STIFEL(1),NOFVAR*NOFVAR)
      ENDIF
C
C
C ************************************************************
C ************************************************************
C
C
C     at this stage the matrix DCHARV is as follows
C     same structure for TAUX
C
C  v       v  a  r  i  a  b  l  e
C  e  dP(1) 
C  r   ...   qdr(1)   ...   qdr(4)
C  t   ...   dp(1)    ...   dp(4) 
C  e  dP(4)  qds(1)   ...   qds(4)
C  x
C     the block involving the last three variables (qdr,dp,qds)
C     is now transposed
C     IADDR gives the location of qdr(1) in DCHARV
C
      IADDR = NOFVERT+1
      CALL TRANS(DCHARV(IADDR),NORDER,NOFVERT,NORDER*NOFVERT,MOVE,
     +           IWRK,IFAIL)
      IF(IFAIL.NE.0)THEN
         WRITE(6,*)'TRANS HAS RETURNED IFAIL = ',IFAIL
         CALL EXIT(IFAIL)
      ENDIF
C
C     so that now DCHARV looks like:
C
C     dP(1)  qdr(1)  dp(1)  qds(1)
C      ...    ...    ...    ...
C      ...    ...    ...    ...
C           
C     dP(4)  qdr(4)  dp(4)  qds(4)
C
C     transpose the timestep
C
      CALL TRANS(TAUX(IADDR),NORDER,NOFVERT,NOFVERT*NORDER,
     +           MOVE,IWRK,IFAIL)
      IF(IFAIL.NE.0)THEN
         WRITE(6,*)'TRANS HAS RETURNED IFAIL = ',IFAIL
         CALL EXIT(IFAIL)
      ENDIF
C
C     compute the timestep: TAUX is transposed, see the table above for DSYMM
C
      DO IVERT = 1, NOFVERT ! loop over the vertices
         JADDR = IVERT
         HELP = ZERO 
C
C     sum or maximum over the dofs
C
         DO IVAR = 1,NOFEQN ! loop over the NDIM+1 dofs
            IF( CHAR_TIMESTEPPING )THEN
               HELP = MAX(HELP,TAUX(JADDR))
            ELSE
               HELP = HELP + TAUX(JADDR)
            ENDIF
            JADDR = JADDR + NOFVERT
         ENDDO ! end loop over the dofs
         IADDR = (IVERT-1)*NOFVAR
         DO IVAR = 1,NOFEQN
            JADDR = IADDR + IVAR
            TSTEP(JADDR) = TSTEP(JADDR) + HELP
         ENDDO
      ENDDO ! loop over the vertices
C
C Transform the nodal residual into conserved variables
C           note that DCHARV is transposed during the MM product
C
      CALL DGEMM('No Transpose','Transpose',NOFEQN,NOFVERT,NOFEQN,
     +           ONE,dUdV,NOFEQN,DCHARV(1),NOFVERT,ZERO,NODRES,
     +           NOFVAR)
C
C     Checks the decomposition ..
C
      IF( ICHECK .NE. 0 )THEN
C
         CALL DINIT(NOFEQN,ZERO,WKSP,1)
         DO 9 IVERT = 1,NOFVERT
             IADDR = (IVERT-1)*NOFVAR+1
             CALL DAXPY(NOFEQN,-ONE/VOLUME,NODRES(IADDR),1,WKSP,1)
    9    CONTINUE 
C
          CALL TEST( DivFlux , WKSP , 1.D-15, IELEM , NOFEQN )
C
      ENDIF
C
C     --------------- If explicit, return now  ---------------
C
      IF (.NOT.PICARD) RETURN
C
C     Add the element stiffness matrix to the global stiffness matrix
C
C
C     transform the element stiffness matrix into conserved variables
C
      IF(NOFVAR.EQ.NOFEQN)THEN
          CALL DGEMM('No Transpose','No Transpose',NOFVAR,
     +               NOFVAR*NOFVERT*NOFVERT,NOFVAR,MONE,dUdV,
     +               NOFVAR,STIFEL,NOFVAR,ZERO,TEMPA,NOFVAR)
      ELSE
          DO 37 JVERT = 1,NOFVERT
             DO 37 IVERT = 1,NOFVERT
                IADDR = IDX(IVERT,JVERT,NOFVAR,NOFVERT)
                JADDR = IDX(IVERT,JVERT,NOFEQN,NOFVERT)
          CALL DGEMM('No Transpose','No Transpose',NOFEQN,
     +         NOFEQN,NOFEQN,MONE,dUdV,NOFEQN,
     +         STIFEL(IADDR),NOFVAR,ZERO,TEMPA(JADDR),NOFEQN)
   37     CONTINUE
      ENDIF
C

              DO 35 JVERT = 1,NOFVERT
                  DO 35 IVERT = 1,NOFVERT
                  IADDR = IDX(IVERT,JVERT,NOFVAR,NOFVERT)
                  JADDR = IDX(IVERT,JVERT,NOFEQN,NOFVERT)
                  CALL DGEMM('No Transpose','No Transpose',NOFEQN,
     +                       NOFEQN,NOFEQN,ONE,TEMPA(JADDR),NOFEQN,
     +                       dVdZ,NOFEQN,ZERO,STIFEL(IADDR),NOFVAR)
C
   35     CONTINUE
C
C
      RETURN
      END
