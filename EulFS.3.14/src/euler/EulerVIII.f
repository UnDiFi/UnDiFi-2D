!> \copydetails EulerIX()
      SUBROUTINE EulerVIII(IELEM,VCN,VCB,VCZ,NDIM,NOFVERT,NOFVAR,
     &                     NTURB,NODRES,TSTEP,STIFEL,VOLUME,
     &                     MATRIX_ASSEMBLY,ScalarScheme,MatrixScheme)
C
C     $Id: EulerVIII.f,v 1.32 2013/05/03 09:55:22 abonfi Exp $
C
C
      IMPLICIT NONE 
C
C     The INCOMPRESSIBLE Euler eqns.
C        (using a pseudo compressibility formulation)
C        are solved in primitive variables (i.e. pressure-velocity)..
C
C
C
      INCLUDE 'paramt.h'
      INCLUDE 'constants.h'
      INCLUDE 'bnd.h'
      INCLUDE 'time.h'
      INCLUDE 'three.com'
      INCLUDE 'transf.com'
      INCLUDE 'flags.com'
      INCLUDE 'bodyf.com'
      INCLUDE 'time.com'
C
C     NEQMAX is the max. no. of equations for the system schemes
C     (3 in 2D,4 in 3D)
C     NOFVAR is actual no. of equations
C     NOFEQN is no. of mean flow eqns. = DIM+1
C
      INTEGER NEQMAX,NOFEQN,FrstEQN
      PARAMETER (NEQMAX=4,FRSTEQN=1)
C
C
      INTEGER IELEM,NDIM,NOFVERT,NOFVAR,NELEM
C
      LOGICAL MATRIX_ASSEMBLY
      LOGICAL LFLAG
C
      EXTERNAL ScalarScheme,MatrixScheme
C
      INTEGER IVAR,IVERT,JVERT,NTURB,I,JCOL,IDIM,IFAIL,N4,
     &IADD,JADD
C
      DOUBLE PRECISION NODRES(*),TSTEP(*)
      DOUBLE PRECISION VOLUME,HELP
      DOUBLE PRECISION VCZ(*),VCN(*),VCB(*),STIFEL(*)
      DOUBLE PRECISION STIFC(MAX_NOFVERT_SQR*NEQMAX*NEQMAX),
     1                 Jacobian(NEQMAX,NEQMAX*3),SOURCE(MAXNOFVAR),
     2                 PHI(NEQMAX),WKSP(5),RESIDUAL(2*NMAX),
     3                 BETA(NEQMAX*NEQMAX*MAXNOFVERT),
     4                 TAUX(MAXNOFVAR*MAXNOFVERT)
C
C     RESIDUAL[1:NOFVAR] stores the residual computed by
C                        the Matrix scheme as \sum K_j U_j
C     RESIDUAL[NOFVAR+1:2*NOFVAR]
C                        stores the residual computed by
C                        the Matrix scheme as \sum_C_{ij} U_j
C     it is used just for debugging purposes, to be compared with
C     the residual computed as:
C     dF/dU * dU/dX + dG/dU * dU/dy + [ dH/dU * dU/dz ]
C
C
      DATA SOURCE/MAXNOFVAR*ZERO/
C
      EXTERNAL MatSplitNum,MatSplitVIII
C
      NOFEQN = NDIM + 1
      N4 = NOFVAR*NOFVAR*NOFVERT*NOFVERT
C
C     The element stiffness matrix is initialized to 0.d0
C
      IF (MATRIX_ASSEMBLY) CALL DINIT(N4,ZERO,STIFC,1)
      CALL DINIT(NOFEQN*NOFVERT,ZERO,TAUX,1)
C
C --------------- Debugging code starts here ---------------
C
      IF (ICHECK.EQ.0) GOTO 7
C
C     Some initializations ....
C
          CALL DINIT(NEQMAX,ZERO,PHI,1)
          CALL DINIT(5,ZERO,WKSP,1)
          CALL CHECK(IELEM,NDIM,NOFEQN)
C
C     Subr. Eigen_VIII computes the jacobian matrix of the 
C           inviscid fluxes; this should also be computed
C           if a numerical decomposition of the matrix 
C           is required.
C           It is anyway required whenever ICHECK <> 0 
C
      CALL Eigen_VIII(Jacobian,NEQMAX,DVDZ,DUDV,NDIM,NOFEQN)
c
c     COMPUTES THE RESIDUAL/VOLUME : dF/dU * dU/dX + dG/dU * dU/dy + ...
c
      DO 12 idim = 1,NDIM
          jcol = (idim-1)*NEQMAX + 1
          CALL DGEMV('N',NOFVAR,NOFVAR,ONE,Jacobian(1,jcol),NEQMAX,
     +               GRAD_PARM(FrstEqn,idim),1,ONE,PHI,1)

   12 CONTINUE
C
C --------------- Debugging code ends here ---------------
C
    7 CONTINUE
C
      SOURCE(2) =-GRAV(1)*VOLUME
      SOURCE(3) =-GRAV(2)*VOLUME
      IF(NDIM.EQ.3)SOURCE(4) =-GRAV(3)*VOLUME
C
      IF(LTIME)DTVOL = DELT/VOLUME
C
C The system is solved using a Matrix Scheme
C
      CALL MatrixScheme(MatSplitVIII,VCZ,NODRES,TAUX,BETA,STIFC,NOFEQN,
     +                  NOFVAR,NOFVERT,VCN,NDIM,Jacobian,NEQMAX,
     +                  RESIDUAL,SOURCE,IELEM,MATRIX_ASSEMBLY)
C
C     compute the timestep
C
      DO IVERT = 1, NOFVERT ! loop over the vertices
         IADD = (IVERT-1)*NOFEQN
         HELP = ZERO 
C
C     sum over the dofs
C
         DO IVAR = 1,NOFEQN
            JADD = IADD + IVAR
            IF( CHAR_TIMESTEPPING )THEN
               HELP = MAX(HELP,TAUX(JADD))
            ELSE
               HELP = HELP + TAUX(JADD)
            ENDIF
         ENDDO
         IADD = (IVERT-1)*NOFVAR
         DO IVAR = 1,NOFEQN
            JADD = IADD + IVAR
            TSTEP(JADD) = TSTEP(JADD) + HELP
         ENDDO
      ENDDO
C
C
C Add the un-steady contribution
C
      IF(LTIME.AND.(NTIMLEVS.GT.1))THEN ! why the 2nd test ?
         CALL UNSTEADY2(Jacobian,BETA,VCZ,NOFVAR,NODRES,STIFC,NOFEQN,
     &                  NDIM,NOFVERT,MATRIX_ASSEMBLY)
      ENDIF
C
C --------------- Debugging code starts here ---------------
C
C     Checks the decomposition ..
C
      IF (ICHECK.NE.0) THEN
          CALL DSCAL(2*NOFVAR,ONE/VOLUME,RESIDUAL,1)
          LFLAG = .TRUE.
          DO 18 IVAR = 1,NOFEQN
              IF (DABS(PHI(IVAR)-RESIDUAL(IVAR)).GT.
     +            1.D-13) LFLAG = .FALSE.
   18     CONTINUE
          IF (LFLAG .EQV. .FALSE.) THEN
              WRITE (6,99999) IELEM
              DO 22 IVAR = 1,NOFVAR
                  WRITE (6,*) PHI(IVAR),RESIDUAL(IVAR),
     +              DABS(PHI(IVAR)-RESIDUAL(IVAR))
   22         CONTINUE
              PAUSE
          ENDIF
          CALL TEST(DivFlux,RESIDUAL,1.D-15,IELEM,NOFEQN)
      ENDIF
C
C --------------- Debugging code ends here ---------------
C
C
C
C --------------- If explicit, return now  --------------------
C
      IF (MATRIX_ASSEMBLY) THEN
C
C     Add the element stiffness matrix to the global stiffness matrix
C
           IF(NOFVAR.EQ.NOFEQN)THEN ! this is the Navier-Stokes case
              DO 33 I = 1,N4
                  STIFEL(I) = -STIFC(I)
   33         CONTINUE
           ELSE ! this is the RANS case
               CALL DSCAL(N4,MONE,STIFC,1)
               CALL MATINS(STIFEL,NOFVAR,STIFC,NOFEQN,NOFVERT,NOFVERT,0)
           ENDIF

      ENDIF
C
      IF (ICHECK.EQ.0) RETURN
C
C --------------- Debugging code starts here ---------------
C
C     test the residual as computed by the "implicit" scheme
C
      CALL TEST(DivFlux,RESIDUAL(NOFVAR+1),1.D-15,IELEM,NOFVAR)
C
C --------------- Debugging code ends here ---------------
C
      RETURN
99999 FORMAT (5X,'Vector residual in Element ',I6,' EulerVIII')
      END
