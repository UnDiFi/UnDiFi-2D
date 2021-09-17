!> \copydetails EulerIX()
      SUBROUTINE EulerIIbis(IELEM,VCN,VCB,VCZ,NDIM,NOFVERT,NOFVAR,
     &NTURB,NODRES,TSTEP,STIFEL,VOLUME,PICARD,ScalarScheme,MatrixScheme)
C
C     Hyperbolic Elliptic splitting using the vanLeer, Lee, Roe
C     preconditioning matrix ..
C
C     $Id: EulerIIbis.f,v 1.29 2020/03/28 09:51:15 abonfi Exp $
C
C
      IMPLICIT NONE
C
C
      INCLUDE 'paramt.h'
      INCLUDE 'constants.h'
      INCLUDE 'time.h'
      INCLUDE 'bnd.h'
      INCLUDE 'blkEulerII.com'
      INCLUDE 'time.com'
      INCLUDE 'three.com'
      INCLUDE 'transf.com'
      INCLUDE 'flags.com'
C
C
C     MAXNORDER is the max. no. of equations (3 in 3D)
C     NORDER is actual no. of equations = DIM
C
      INTEGER MAXNORDER
      PARAMETER(MAXNORDER=3)
      INTEGER IELEM,NDIM,NOFVERT,NOFVAR,NTURB
C
      EXTERNAL ScalarScheme,MatrixScheme
C
      INTEGER FrstEq,JCOL,IDIM,ifail,IWRK,IADDR,JADDR,KADDR,NOFEQN
      PARAMETER (IWRK=10,FRSTEQ=3)
      INTEGER NORDER,MOVE(IWRK),IADD,JADD
      DOUBLE PRECISION SCALRES(2)
      LOGICAL PICARD
C
C
      DOUBLE PRECISION Jacobian(MAXNORDER,MAXNORDER*3)
      DOUBLE PRECISION PHI(MAXNORDER)
      DOUBLE PRECISION NODRES(*),WKSP(MAXNOFEQN)
      DOUBLE PRECISION TSTEP(*),RESIDUAL(2*MAXNORDER)
      INTEGER IDX,I,J,N,M,N4,IVERT,JVERT,IVAR,ORDSQR
      integer lwork
      parameter(lwork=2*maxnofeqn) 
      DOUBLE PRECISION VCZ(*),VCN(*),VCB(*),STIFEL(*),
     +VOLUME(*),STIFC(144),CHARV(MAXNOFEQN*MAXNOFVERT),
     +DCHARV(MAXNOFEQN*MAXNOFVERT),TEMPC(MAXNOFVERT*MAXNOFEQN**2),
     +TEMPA(MAXNOFEQN**2*MAX_NOFVERT_SQR),SOURCE(MAXNOFVAR),
     +TEMPB(MAXNOFEQN*MAXNOFEQN*MAXNOFVERT),TAUX(MAXNOFEQN*MAXNOFVERT),
     +TEMPD(MAXNOFEQN*MAXNOFEQN),work(lwork)
      integer ipiv(lwork) 
      DOUBLE PRECISION TOLER,HELP
      PARAMETER(TOLER=1.D-15)
      LOGICAL LFLAG
      LOGICAL UNITMAt
C
C     STIFEL(NOFVAR,NOFVAR,NOFVERT,NOFVERT)
C     STIFC (NORDER,NORDER,NOFVERT,NOFVERT) for matrix schemes and
C     STIFC (NOFVERT,NOFVERT) for scalar schemes
C
C
      EXTERNAL MatSplitNum,MatSplitII
      DATA SOURCE/MAXNOFVAR*ZERO/
C
      IDX(I,J,N,M) = (((J-1)*M+I-1)*N*N)+1
C
      NORDER = NDIM
      ORDSQR = NORDER*NORDER
      NOFEQN = NDIM+2
      N4 = NOFVAR*NOFVAR*NOFVERT*NOFVERT
C
C     The element stiffness matrix is initialized to 0.d0
C
      IF (PICARD) THEN
          CALL DINIT(N4,ZERO,STIFEL,1)
          CALL DINIT((NORDER*NOFVERT)**2,ZERO,STIFC,1)
      ENDIF
C
      CALL DINIT(NOFVERT*NOFEQN,ZERO,DCHARV,1)
      CALL DINIT(NOFVERT*NOFEQN,ZERO,TAUX,1)
C
C     Sets up a stream aligned frame ..
C
      CALL StreamAlignedFrame(NDIM)
C
C     The Jacobian Matrix of the subsystem is assembled and
C         the eigenvectors computed ..
C
      CALL Eigen_IIbis(Jacobian,MAXNORDER,DVDZ,DUDV,NDIM,NOFEQN)
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
C     call these (S,H,P,Q,R)
C
      CALL DGEMM('Transpose','Transpose',NOFVERT,NOFEQN,NOFEQN,
     +           ONE,VCZ,NOFVAR,dVdZ,NOFEQN,ZERO,CHARV,NOFVERT)
C
C     The matrix CHARV now looks like:
C
C     dS(1) dH(1)  dP(1)  dQ(1)  dR(1)
C     ....   ...    ...    ...    ...
C     ....   ...    ...    ...    ...
C           
C     dS(4) dH(4)  dP(4)  dQ(4)  dR(4)
C
C     the matrix of the characteristic variables P,Q,R
C     is now transposed 
C
      IADDR = 2*NOFVERT+1 
      CALL TRANS(CHARV(IADDR),NOFVERT,NORDER,NOFVERT*NORDER,
     +           MOVE,IWRK,IFAIL)
      IF(IFAIL.NE.0)WRITE(6,*)'TRANS HAS RETURNED IFAIL = ',IFAIL
C
C     so that these are now stored 
C     (starting at CHARV(IADDR)=dP(1)) 
C     as CHARV(1:NORDER,1:NOFVERT)
C
C      dP(1)   ...   dP(4)
C      dQ(1)   ...   dQ(4) 
C      dR(1)   ...   dR(4)
C
C ************************************************************
C The non commuting system is solved using a matrix Scheme
C ************************************************************
C
caldo CALL R8Mat_Print('G',' ',NORDER,NOFVERT,CHARV(IADDR),NORDER,
caldo+            'W matrix',IFAIL)
caldo CALL R8Mat_Print('G',' ',NORDER,NOFVERT,DCHARV(IADDR),NORDER,
caldo+            'dW matrix',IFAIL)
      CALL MatrixScheme(MatSplitII,CHARV(IADDR),DCHARV(IADDR),
     +                  TAUX(IADDR),TEMPA,STIFC,NORDER,NORDER,NOFVERT,
     +                  VCN,NDIM,Jacobian,MAXNORDER,RESIDUAL,SOURCE,
     +                  IELEM,PICARD)
C
!           CALL R8Mat_Print('G',' ',NORDER,NORDER*nofvert,TEMPA,
!    +      NORDER,'Beta for P,Q,R ',IFAIL)
C
caldo CALL R8Mat_Print('G',' ',NORDER,NOFVERT,DCHARV(IADDR),NORDER,
caldo+            'dW matrix',IFAIL)
caldo   pause
C
C    Checks the cell residual in characteristic variables
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
C
          ENDIF
      ENDIF
C
      IF (LTIME) THEN
          CALL DINIT(NOFEQN*NOFEQN*NOFVERT,ZERO,TEMPC,1) ! set distr. matrices == 0
          CALL MATINS(TEMPC,NOFEQN,TEMPA,NORDER,NOFVERT,1,2)
          CALL DINIT(ORDSQR*NOFVERT*NOFVERT,ZERO,TEMPA,1) ! should be useless
!         CALL R8Mat_Print('G',' ',NOFEQN,NOFEQN*nofvert,TEMPC,
!    +                NOFEQN,'Beta (2) for P,Q,R ',IFAIL)
      ENDIF
      IF (PICARD) THEN
          CALL MATINS(STIFEL,NOFVAR,STIFC,NORDER,NOFVERT,NOFVERT,2)
          CALL DINIT(ORDSQR*NOFVERT*NOFVERT,ZERO,STIFC,1)
      ENDIF
C
C ************************************************************
C     Solve the entropy transport eqn. 
C ************************************************************
C
      CALL ScalarScheme(IELEM,VCN,R_SPEED(1,1),SCALRES(1),ZERO,
     +                  CHARV(1),TAUX(1),DCHARV(1),TEMPA,STIFC,
     +                  NDIM,NOFVERT,PICARD)
C
C     copy entries of the entropy distribution matrix
C     into the "full" distribution matrix
C
      IF (LTIME) THEN
!         CALL R8Mat_Print('G',' ',1,nofvert,TEMPA,
!    +                1,'Beta for S ',IFAIL)
          CALL DCOPY(NOFVERT,TEMPA,1,TEMPC(1),NOFEQN*NOFEQN)
          CALL DINIT(NOFVERT*NOFVERT,ZERO,TEMPA,1) ! should be useless
!         CALL R8Mat_Print('G',' ',NOFEQN,NOFEQN*nofvert,TEMPC,
!    +                NOFEQN,'Beta (2) for S,P,Q,R ',IFAIL)
      ENDIF
C
C    Copy the convective jacobian (STIFC) into STIFEL and
C    reset it to 0. since it will be reused by the scalar schemes
C
      IF (PICARD) THEN
          CALL DCOPY(NOFVERT*NOFVERT,STIFC,1,STIFEL(1),NOFVAR*NOFVAR)
          CALL DINIT(NOFVERT*NOFVERT,ZERO,STIFC,1)
      ENDIF
C
C ************************************************************
C     Solve the total enthalpy transport eqn. 
C ************************************************************
C
      IADDR = NOFVERT+1
      CALL ScalarScheme(IELEM,VCN,R_SPEED(1,2),SCALRES(2),ZERO,
     +                  CHARV(IADDR),TAUX(IADDR),DCHARV(IADDR),TEMPA,
     +                  STIFC,NDIM,NOFVERT,PICARD)
C
C     BETA(i,j) is copied into TEMPC(2,2,i,j)
C
      IF (LTIME) THEN
!         CALL R8Mat_Print('G',' ',1,nofvert,TEMPA,
!    +                1,'Beta for H ',IFAIL)
          IADDR = NOFEQN+2 
          CALL DCOPY(NOFVERT,TEMPA,1,TEMPC(IADDR),
     +               NOFEQN*NOFEQN)
!           CALL R8Mat_Print('G',' ',NOFEQN,NOFEQN*nofvert,TEMPC,
!    +      NOFEQN,'Beta (2) for H,S,P,Q,R ',IFAIL)
      ENDIF
      IF (PICARD) THEN
C
C     STIFC(i,j) is copied into STIFEL(2,2,i,j)
C
          IADDR = NOFVAR+2 
          CALL DCOPY(NOFVERT*NOFVERT,STIFC,1,STIFEL(IADDR),
     +               NOFVAR*NOFVAR)
      ENDIF
C
C     at this stage the matrix DCHARV is as follows
C
C  v       v  a  r  i  a  b  l  e
C  e  dS(1) dH(1) 
C  r  ....   ...   dP(1)   ...   dP(4)
C  t         ...   dQ(1)   ...   dQ(4) 
C  e  dS(4) dH(4)  dR(1)   ...   dR(4)
C  x
C     the block involving the last three variables (dP,dQ,dR)
C     is now transposed
C     IADDR gives the location of dP(1) in DCHARV
C
      IADDR = 2*NOFVERT+1
      CALL TRANS(DCHARV(IADDR),NORDER,NOFVERT,NORDER*NOFVERT,MOVE,
     +           IWRK,IFAIL)
      IF(IFAIL.NE.0)WRITE(6,*)'TRANS HAS RETURNED IFAIL = ',IFAIL
C
C     so that now DCHARV looks like:
C
C     dS(1) dH(1)  dP(1)  dQ(1)  dR(1)
C     ....   ...    ...    ...    ...
C     ....   ...    ...    ...    ...
C           
C     dS(4) dH(4)  dP(4)  dQ(4)  dR(4)
C
C     a little trick so that we can directly feed
C     TAUX to the scalar schemes:
C
      CALL TRANS(TAUX(IADDR),NORDER,NOFVERT,NORDER*NOFVERT,MOVE,IWRK,
     &IFAIL)
      IF(IFAIL.NE.0)WRITE(6,*)'TRANS HAS RETURNED IFAIL = ',IFAIL
!     CALL R8Mat_Print('G',' ',NOFVERT,NOFEQN,TAUX,
!    +      NOFVERT,'time step is :',IFAIL)
C
C     compute the timestep: TAUX is transposed
C
      DO IVERT = 1, NOFVERT ! loop over the vertices
         JADD = IVERT
         HELP = ZERO 
C
C     sum or maximum over the dofs
C
         DO IVAR = 1,NOFEQN ! loop over the NDIM+2 dofs
            IF( CHAR_TIMESTEPPING )THEN
               HELP = MAX(HELP,TAUX(JADD))
            ELSE ! standard (default) approach
               HELP = HELP + TAUX(JADD)
            ENDIF
            JADD = JADD + NOFVERT
         ENDDO ! end loop over the dofs
         IADD = (IVERT-1)*NOFVAR
         DO IVAR = 1,NOFEQN
            JADD = IADD + IVAR
            TSTEP(JADD) = TSTEP(JADD) + HELP
         ENDDO
      ENDDO ! end loop over the vertices
C
C Transform the nodal residual into conserved variables
C           note that DCHARV is transposed during the MM product
C
      CALL DGEMM('No Transpose','Transpose',NOFEQN,NOFVERT,NOFEQN,
     +           ONE,dUdV,NOFEQN,DCHARV(1),NOFVERT,ZERO,NODRES,
     +           NOFVAR)
C
C     Checks that the residual in conserved variables
C     equals the flux integral computed as dF/dZ Z_x +
C
      IF( ICHECK .NE. 0 )THEN
C
         CALL DINIT(NOFEQN,ZERO,WKSP,1)
         DO 9 IVERT = 1,NOFVERT
             IADDR = (IVERT-1)*NOFVAR+1
             CALL DAXPY(NOFEQN,MONE/VOLUME(1),NODRES(IADDR),1,WKSP,1)
    9    CONTINUE 
C
          CALL TEST( DivFlux , WKSP , TOLER, IELEM , NOFEQN )
C
      ENDIF
C
      IF(LTIME)THEN 
         DTVOL = DELT/VOLUME(1)
c
c     compute transformation matrix dW/dU
c
copt1    CALL CONS2PARM(ZAVG,dZdU,NDIM,NOFEQN)
copt1    CALL DGEMM('No Transpose','No Transpose',NOFEQN,
copt1+               NOFEQN,NOFEQN,ONE,dVdZ,NOFEQN,dZdU,NOFEQN,
copt12               ZERO,TEMPB,NOFEQN)
c
copt1 call dcopy(NOFEQN*NOFEQN,TEMPB,1,TEMPD,1)
copt1 CALL DGETRF(NOFEQN,NOFEQN,TEMPD,NOFEQN,IPIV,IFAIL)
copt1 CALL DGETRI(NOFEQN,TEMPD,NOFEQN,IPIV,WORK,LWORK,IFAIL)
caldo
caldo compute the inverse of M P^{-1} R
caldo
      call dcopy(NOFEQN*NOFEQN,dUdV,1,TEMPB,1)
      CALL DGETRF(NOFEQN,NOFEQN,TEMPB,NOFEQN,IPIV,IFAIL)
      CALL DGETRI(NOFEQN,TEMPB,NOFEQN,IPIV,WORK,LWORK,IFAIL)
caldo
caldo    Here we check that the betas sum up to the identity matrix
caldo
!        CALL DINIT(NOFEQN*NOFEQN*NOFVERT,ZERO,TEMPA,1)
!        DO IVERT = 1,NOFVERT ! si puo forse ridurre il nof di DGEMM !?!?!?!?!
!           JADD = (IVERT-1)*NOFEQN*NOFEQN+1
!           CALL DAXPY(NOFEQN*NOFEQN,ONE,TEMPC(JADD),1,TEMPA,1) 
!        ENDDO
!           IF(.NOT.UNITMAT(TEMPA,NOFEQN,NOFEQN,NOFEQN,1.D-12))THEN
!     CALL R8Mat_Print('G',' ',NOFEQN,NOFEQN,TEMPA,
!    +      NOFEQN,'should be I ',IFAIL)
!           pause
!           ENDIF
caldo
C
!           CALL R8Mat_Print('G',' ',NOFEQN,NOFEQN*nofvert,TEMPC,
!    +      NOFEQN,'Beta in char vars',IFAIL)
C
         DO IVERT = 1,NOFVERT ! si puo forse ridurre il nof di DGEMM !?!?!?!?!
            JADD = (IVERT-1)*NOFEQN*NOFEQN+1
            CALL DGEMM('No Transpose','No Transpose',NOFEQN,
     +               NOFEQN,NOFEQN,ONE,dUdV,NOFEQN,TEMPC(JADD),NOFEQN,
     2               ZERO,TEMPA(JADD),NOFEQN) 
copt1       CALL DGEMM('No Transpose','No Transpose',NOFEQN,
copt1+               NOFEQN,NOFEQN,ONE,TEMPD,NOFEQN,TEMPC(JADD),NOFEQN,
copt12               ZERO,TEMPA(JADD),NOFEQN) 
            CALL DGEMM('No Transpose','No Transpose',NOFEQN,
     +               NOFEQN,NOFEQN,ONE,TEMPA(JADD),NOFEQN,TEMPB,NOFEQN,
     2               ZERO,TEMPC(JADD),NOFEQN) 
         ENDDO
!        goto 392
!           CALL DGEMM('No Transpose','No Transpose',NOFEQN,
!    +               NOFEQN*NOFVERT,NOFEQN,ONE,dUdV,NOFEQN,TEMPC(JADD),
!    2               NOFEQN,ZERO,TEMPA(JADD),NOFEQN) 
!           CALL DGEMM('No Transpose','No Transpose',NOFEQN,
!    +               NOFEQN*NOFVERT,NOFEQN,ONE,TEMPA(JADD),NOFEQN,TEMPB,
!    2               NOFEQN,ZERO,TEMPC(JADD),NOFEQN) 
! 392 continue
C
!        CALL R8Mat_Print('G',' ',NOFEQN,NOFEQN*nofvert,TEMPC,
!    +               NOFEQN,'Beta in conserved vars',IFAIL)
C
C        compute nodal values of the dudz matrix
C
         DO IVERT = 1,NOFVERT
              IADD = (IVERT-1)*NOFVAR+1
              JADD = (IVERT-1)*NOFEQN*NOFEQN+1
              CALL PARM2CONS(VCZ(IADD),TEMPB(JADD),NOFEQN,NDIM)
         ENDDO
C
         IF(DUALTS)
     1   CALL UNSTEADY4(TEMPB,TEMPC,VCZ,NOFVAR,NODRES,STIFEL,VOLUME,
     &                  NOFEQN,NDIM,NOFVERT,PICARD)
      ENDIF
C
C     --------------- If explicit, return now  ---------------
C
      IF (.NOT.PICARD) RETURN
C
C     compute the transformation matrices from
C     conserved to parameter variables in the vertices
C
          DO 33 IVERT = 1,NOFVERT
              JADDR = (IVERT-1)*NOFEQN*NOFEQN+1
              CALL MatdZdU(VCZ((IVERT-1)*NOFVAR+1),dZdU(JADDR),
     +                     NDIM,NOFEQN)
   33     CONTINUE
C
C     transform the element stiffness matrix into conserved variables
C
      IF(NOFVAR.EQ.NOFEQN)THEN
          CALL DGEMM('No Transpose','No Transpose',NOFVAR,
     +               NOFVAR*NOFVERT*NOFVERT,NOFVAR,MONE,dUdV,
     +               NOFVAR,STIFEL,NOFVAR,ZERO,TEMPA,NOFVAR)
      ELSE
          DO 37 IVERT = 1,NOFVERT
          DO 37 JVERT = 1,NOFVERT
              IADDR = IDX(IVERT,JVERT,NOFVAR,NOFVERT)
              JADDR = IDX(IVERT,JVERT,NOFEQN,NOFVERT)
              CALL DGEMM('No Transpose','No Transpose',NOFEQN,
     +               NOFEQN,NOFEQN,MONE,dUdV,NOFEQN,
     +               STIFEL(IADDR),NOFVAR,ZERO,TEMPA(JADDR),NOFEQN)
   37     CONTINUE
      ENDIF 
          CALL DGEMM('No Transpose','No Transpose',NOFEQN,
     +               NOFEQN*NOFVERT,NOFEQN,TWO,dVdZ,NOFEQN,dZdU,
     +               NOFEQN,ZERO,TEMPB,NOFEQN)
C
C

              DO 35 JVERT = 1,NOFVERT
                  JADDR = (JVERT-1)*NOFEQN*NOFEQN+1
                  DO 35 IVERT = 1,NOFVERT
                  IADDR = IDX(IVERT,JVERT,NOFVAR,NOFVERT)
                  KADDR = IDX(IVERT,JVERT,NOFEQN,NOFVERT)
                  CALL DGEMM('No Transpose','No Transpose',NOFEQN,
     +                       NOFEQN,NOFEQN,ONE,TEMPA(KADDR),NOFEQN,
     +                       TEMPB(JADDR),NOFEQN,ZERO,STIFEL(IADDR),
     +                       NOFVAR)
!     write(6,*)ivert,jvert
!     CALL R8Mat_Print('G',' ',NOfvar,NOfvar,stifel(IADDR),NOfvar,
!    +            'Cij matrix',IFAIL)
!     pause
C
   35     CONTINUE
C
      IF(ICHECK.EQ.0)RETURN
C
C --------------- Debugging code starts here ---------------
C
C     test the residual as computed by the "implicit" scheme
C     WKSP := dUdV(1) * \phi_{entropy} 
C             dUdV(2) * \phi_{enthalpy} 
C           + dUdV(1,FrstEq) * \Phi
C
      CALL DINIT(NOFEQN,ZERO,WKSP,1)
      CALL DAXPY(NOFEQN,SCALRES(1)/VOLUME(1),dUdV,1,WKSP,1)
      CALL DAXPY(NOFEQN,SCALRES(2)/VOLUME(1),dUdV(NOFEQN+1),1,WKSP,1)
      CALL DGEMV('N',NOFEQN,NORDER,ONE,dUdV((FrstEq-1)*NOFEQN+1),
     +           NOFEQN,RESIDUAL(NORDER+1),1,ONE,WKSP,1)
      CALL TEST(DivFlux,WKSP,TOLER,-IELEM,NOFEQN)
C
C --------------- Debugging code ends here ---------------
C

89999 FORMAT (5X,'Scalar residual in Element ',I6,' Wave # ',I1)
99999 FORMAT (5X,'Vector residual in Element ',I6,' EulerIIbis')

C
      RETURN
      END
