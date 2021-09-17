!> \par Purpose
!>
!> Compute the integral of the flux \c FLUXF through a linear/triangular face
!> using Simpson's rule
!>
!> @param[in] NDIM dimension of the space
!> @param[in] NOFVERT nof vertices for a cell = NDIM+1
!> @param[in] NDOF nof degrees of freedom in vector ZROE, also leading dimension of FLXN
!> @param[in] ZROE Roe's parameter vector
!> @param[in] XYZDOT the cartesian components of the grid velocity \f$ \mathbf{b}\f$
!> @param[in] FACN the cartesian components of the scaled inward normal to the face
!> @param[out] FLXN the eulerian (inviscid) nodal fluxes: the first NDIM columns return the fluxes in the NDIM vertices of the face,
!>             the last one (NOFVERT) is used to store the face value computed according to Simpson's rule
!> @param[in] FLUXF is the subroutine used to compute the flux
!>            the calling sequence is: 
!> \verbatim
!> FLUXF((NDIM,FACN,ZROE(1,I),XYZDOT,FLXN(1,I))
!> \endverbatim
C
!> \author $Author: abonfi $
!> \version $Revision: 1.7 $
!> \date $Date: 2020/03/28 09:54:05 $
C
!> \bug In 3D it uses the midpoint rule, NOT Simpson's
      SUBROUTINE SIMPSON(NDIM,NOFVERT,NDOF,ZROE,XYZDOT,FACN,FLXN,FLUXF)
      IMPLICIT NONE
C
      INCLUDE 'paramt.h'
      INCLUDE 'constants.h'
C
      INTEGER NDIM,NOFVERT,NDOF
C
      DOUBLE PRECISION FACN(NDIM),ZROE(NDOF,NOFVERT),
     &                 XYZDOT(NDIM,NOFVERT),FLXN(NDOF,NOFVERT)
C
      DOUBLE PRECISION WC,WO,WD,HELP
      DOUBLE PRECISION DNRM2
C
      EXTERNAL FLUXF
C
      INTEGER I ,IFAIL
      LOGICAL VERBOSE
      PARAMETER(VERBOSE=.FALSE.)
!     PARAMETER(VERBOSE=.TRUE.)
C
      IF(NDIM.EQ.2)THEN
         WC = 4.d0/6.d0
         WO = ONE/6.d0
      ELSE
!        WRITE(6,*)"Simpson's flux in 3D should be checked"
!
! the rule should be Area/3*(f((1+2)/2)+f((2+3)/2)+f((3+1)/2)) i.e.
! compute the flux ad midpoints
!
! https://math.stackexchange.com/questions/1353599/quadrature-formula-on-triangle?answertab=active#tab-top
!
!
         WC = 0.d0
         WO = ONE/3.d0
      ENDIF
      WD = ONE/REAL(NOFVERT-1)
C
      CALL DINIT(NDOF,ZERO,ZROE(1,NOFVERT),1)
      CALL DINIT(NDOF,ZERO,FLXN(1,NOFVERT),1)
      CALL DINIT(NDIM,ZERO,XYZDOT(1,NOFVERT),1)
      DO I = 1,(NOFVERT-1)
            CALL DAXPY(NDOF,ONE,ZROE(1,I),1,ZROE(1,NOFVERT),1)
            CALL DAXPY(NDIM,ONE,XYZDOT(1,I),1,XYZDOT(1,NOFVERT),1) ! compute average grid velocity
            CALL FLUXF(NDIM,FACN,ZROE(1,I),XYZDOT(1,I),FLXN(1,I)) ! compute flux in the vertices of the face
      ENDDO
      CALL DSCAL(NDOF,WD,ZROE(1,NOFVERT),1) ! compute average Z in the edge/face
      CALL DSCAL(NDIM,WD,XYZDOT(1,NOFVERT),1) ! compute average grid velocity in the edge/face
      CALL FLUXF(NDIM,FACN,ZROE(1,NOFVERT),XYZDOT(1,NOFVERT),
     &           FLXN(1,NOFVERT)) ! compute flux in the center of gravity of the face
C
      IF(VERBOSE)THEN
         CALL R8Mat_Print('General',' ',NDOF,NOFVERT,FLXN,
     +            NDOF,'Nodal fluxes (vertex + gravity) ',IFAIL)
      ENDIF
C
      CALL DSCAL(NDOF,WC,FLXN(1,NOFVERT),1)
C
      IF(VERBOSE)THEN
         CALL R8Mat_Print('General',' ',NDOF,NOFVERT,ZROE,
     +            NDOF,'Nodal values incl. average ',IFAIL)
      ENDIF
C
      DO I = 1,(NOFVERT-1)
         CALL DAXPY(NDOF,WO,FLXN(1,I),1,FLXN(1,NOFVERT),1)
      ENDDO
!
      IF(VERBOSE)THEN
         HELP = DNRM2(NDIM,FACN,1)
         WRITE(6,*)'Face area is ',HELP
         CALL R8Mat_Print('General',' ',NDOF,NOFVERT,FLXN,
     +            NDOF,'Nodal fluxes in subr. Simpson ',IFAIL)
      ENDIF
C
      RETURN
      END 
