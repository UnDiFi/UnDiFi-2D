!> \par Purpose
!>
!> Here we compute the signals \c NODRES sent to the vertices of cell \c IELEM
!> when dealing with scalar convection equations:
!>
!> \f[
!> \frac{\partial u}{\partial t} + \mathbf{\lambda} \cdot \nabla u = f
!> \f]
!>
!> We compute the cell averaged velocity \f$\hat{\lambda}-\hat{\mathbf{b}}\f$ relative to the grid, where:
!>
!> \f[
!> \hat{\lambda} = \frac{1}{d+1} \sum_{j=1}^{d+1} u_j \quad\quad \hat{\mathbf{b}} = \frac{1}{d+1} \sum_{j=1}^{d+1} \mathbf{b}_j
!> \f]
!> and
!> \f$(\mathbf{b} = 0)\f$ on a fixed grid 
!>
!> We compute the source term \f${\tt SOURCE} = -\int_{T_e} f \mathrm{d}V\f$
!>
!> We compute the signals due to the inviscid contribution \f$\int_{T_e}\left(\hat{\lambda}-\hat{\mathbf{b}}\right)\cdot\nabla u\,\mathrm{d}V\f$ 
!>
!> We compute the signals due to the temporal contribution 
!>
!> @param[in] IELEM current element
!> @param[in] VCN \c NDIM Cartesian components of the \c NOFVERT inward normals for the current cell
!> @param[in] VCB \c NDIM Cartesian components of the \c NOFVERT grid velocity vectors
!> @param[in] VCZ \c dependent variable in the \c NOFVERT vertices of the current cell
!> @param[in] NDIM dimensionality of the space
!> @param[in] NOFVERT nof vertices of the current cell \c (=NDIM+1) dimensionality of the space
!> @param[in] NOFVAR nof dofs within each vertex \c (NOFVAR=1) for a scalar problem
!> @param[in] NDUMMY dummy integer used for compatibility with similar subroutines sharing the same calling sequence
!> @param[in,out] NODRES the nodal residual is updated with the signals scattered to the vertices of the current cell
!> @param[in,out] TSTEP the nodal timestep (actually \f$V_i/\Delta t_i\f$) is updated with contributions from the current cell
!> @param[in,out] STIFEL the elemental Jacobian matrix (\f$ C_{ij} = \frac{\partial R_i}{\partial u_j}\f$ is updated with contributions from the current cell, only if \c PICARD is set to \c .TRUE.
!> @param[in] VOLUME is the array with cell volumes at times: current, \c n+1, \c n, \c n-1
!> @param[in] PICARD should be set to \c TRUE when the Jacobian matrix has to be computed analytically
!> @param[in] SCALARSCHEME is the Fluctuation Splitting scheme to be used to discretize the convective term of the governing PDEs
!> @param[in] MATRIXSCHEME is unused and left for compatibility with similar subroutines sharing the same calling sequence
!> \author $Author: abonfi $
!> \version $Revision: 1.26 $
!> \date $Date: 2013/09/18 10:39:31 $
!>
!>
      SUBROUTINE SCALAR(IELEM,VCN,VCB,VCZ,NDIM,NOFVERT,NOFVAR,NDUMMY,
     +                  NODRES,TSTEP,STIFEL,VOLUME,PICARD, 
     +                  SCALARSCHEME,MATRIXSCHEME)
C
      IMPLICIT NONE
C
C     $Id: scalar.f,v 1.26 2013/09/18 10:39:31 abonfi Exp $
C
      INCLUDE 'paramt.h'
      INCLUDE 'constants.h'
      INCLUDE 'time.h'
C
C
C     FORTRAN stack
C
      DOUBLE PRECISION DSTAK(1)
      COMMON /CSTAK/ DSTAK
      INTEGER ISTAK(1)
      EQUIVALENCE(DSTAK(1),ISTAK(1))
C
C     ..
C     .. Common blocks ..
C     ..
C
      INCLUDE 'flags.com'
      INCLUDE 'io.com'
      INCLUDE 'nloc.com'
      INCLUDE 'three.com'
      INCLUDE 'time.com'
C
C
C     .. Parameters ..
C     ..
C     .. Scalar Arguments ..
      INTEGER NDIM,NOFVAR,NOFVERT,NDUMMY
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION TSTEP(*),NODRES(*),
     +STIFEL(NOFVERT,NOFVERT),VCN(*),VCZ(*),VCB(*),VOLUME(*)
C     ..
C     .. Subroutine Arguments ..
      EXTERNAL SCALARSCHEME,MATRIXSCHEME
C     ..
C     .. Arrays in Common ..
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION ERR,RESIDUAL,S,SOURCE,DIVB
      INTEGER IELEM,IFAIL,LOCA,LOCB,IOFF,I,J,IPOIN
      LOGICAL PICARD
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION LAMBDA(3),VCP(3*MAXNOFVERT),
     2STIFC(MAX_NOFVERT_SQR),BETA(MAXNOFVERT)
      INTEGER ICN(MAXNOFVERT)
C     ..
C     .. External Functions ..
      DOUBLE PRECISION DDOT,DNRM2,DIV
      EXTERNAL DDOT,DNRM2,DIV
C     ..
C     .. External Subroutines ..
      EXTERNAL ADVECT,DINIT
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,MAX,SIGN
C
C     .. the following is also included in SUPG_scheme.f
      DOUBLE PRECISION VOL
      COMMON/ABC/VOL
      DOUBLE PRECISION sumdiv
      COMMON/DEF/sumdiv
C
C     Sets residual and local timestep to zero
C
      CALL DINIT(NOFVERT,ZERO,TSTEP,1)
      CALL DINIT(NOFVERT,ZERO,NODRES,1)
C
      LOCA = LCELNOD+(IELEM-1)*NOFVERT -1
      DO 10 I = 1, NOFVERT
         IPOIN = ISTAK(LOCA+I)
         LOCB = LCORG + (IPOIN-1)*NDIM-1
         IOFF = (I-1)*NDIM
         DO 10 J = 1, NDIM 
         VCP(IOFF+J) = DSTAK(LOCB+J)
   10 CONTINUE
C
c     computes the advection vector and source term
c
      CALL ADVECT(IELEM,LAMBDA,VCP,VCZ,NDIM,NOFVERT,SOURCE)
c
      SOURCE = -SOURCE*VOLUME(1)
c 
c
      IF(LALE)THEN
         CALL DAXPY(NDIM,MONE,BAVG,1,LAMBDA,1) ! substracts the average grid velocity for ALE calculations
c
c     Not a very good solution, but.... for the time being 
c     When using the explicit LW scheme, DUALTS = .FALSE. and the term -ZAVG(1)*DIVB
c     should be added to the source term. This is done in UNSTEADY1 with all other
c     schemes
c
         IF(.NOT.DUALTS)THEN
            DIVB = DIV(NDIM,NOFVERT,VCN,VCB) ! computes the divergence of the grid velocity
            SOURCE = SOURCE - ZAVG(1)*DIVB ! adds the term: -\div \mathbf{b} <u>
         ENDIF
      ENDIF
c
c
      IF (ICHECK.NE.0) THEN
c
c Dots the adv. speed with the gradient
c
          S = DDOT(NDIM,LAMBDA,1,GRAD_PARM(1,1),NMAX)
C
      ENDIF
C
C     ugly,ugly,ugly: for the SUPG scheme..
C
      VOL = VOLUME(1)
C
C     ugly,ugly,ugly: for the LW scheme..
C
      IF(LTIME)DTVOL = DELT/VOLUME(1)
C
C     ugly,ugly,ugly
C
      CALL SCALARSCHEME(IELEM,VCN,LAMBDA,RESIDUAL,SOURCE,VCZ,TSTEP,
     +                  NODRES,BETA,STIFC,NDIM,NOFVERT,PICARD)
C
      IF(LTIME.AND.DUALTS)THEN 
          CALL UNSTEADY1(BETA,BETA,VCZ,1,NODRES,STIFC,VOLUME,
     2                   1,NDIM,NOFVERT,MMTYPE,PICARD)
      ENDIF
C
C
      IF (ICHECK.NE.0) THEN
              RESIDUAL = RESIDUAL/VOLUME(1)
              ERR = ABS(S-RESIDUAL)
              IF (ERR.GT.5.D-15) THEN
                  WRITE (NOUT,FMT=200) IELEM,S,RESIDUAL,ERR
                  PAUSE

              ENDIF

          ENDIF
C
C     CALL dcopy(nofvert*nofvert,stifc,1,stif_copy,1)
C     CALL daxpy(nofvert*nofvert,-1.d0,stifd,1,stif_copy,1)
C
C     Assembling the Stiffness matrix ...
C
      IF (.NOT.PICARD) RETURN
C
      CALL DSCAL(NOFVERT*NOFVERT,MONE,STIFC,1) 
      CALL DCOPY(NOFVERT*NOFVERT,STIFC,1,STIFEL,1) 
C
      RETURN

  200 FORMAT (5X,'Error on scalar residual in ELEM # ',I6,/,12X,'true',
     +       17X,'computed',14X,'error',/,3 (10X,D12.5))
C
      END
