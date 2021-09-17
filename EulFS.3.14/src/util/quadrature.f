!> \par Purpose
!>
!> Compute the integral of the flux \c FLUXF through a triangular face
!> using a quadrature formula which is exact at least for quadratic polynomials
!> we use a three points rule
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
!> \version $Revision: 1.2 $
!> \date $Date: 2020/03/28 09:54:05 $
C
!> \bug Not sure the grid velociy is accounted for correctly
      SUBROUTINE QUADRATURE(NDIM,NOFVERT,NDOF,ZROE,XYZDOT,FACN,FLXN,
     &FLUXF)
      IMPLICIT NONE
C
      INCLUDE 'paramt.h'
      INCLUDE 'constants.h'
C
      INTEGER NDIM,NOFVERT,NDOF
      INTEGER NM1
C
      DOUBLE PRECISION FACN(NDIM),ZROE(NDOF,NOFVERT),
     &                 XYZDOT(NDIM,NOFVERT),FLXN(NDOF,NOFVERT)
C
      DOUBLE PRECISION XYW(4,3)
!     PARAMETER(XYW(1,1)=1.66666666666666657414808128123695d-01)
!    1          XYW(2,1)=1.66666666666666657414808128123695d-01,
!    2          XYW(3,1)=6.66666666666666740681534975010436d-01,
!    3          XYW(4,1)=3.33333333333333314829616256247391d-01,
!    4          XYW(1,2)=1.66666666666666657414808128123695d-01,
!    5          XYW(2,2)=6.66666666666666740681534975010436d-01,
!    6          XYW(3,2)=1.66666666666666629659232512494782d-01,
!    7          XYW(4,2)=3.33333333333333314829616256247391d-01,
!    8          XYW(1,3)=6.66666666666666740681534975010436d-01,
!    9          XYW(2,3)=1.66666666666666657414808128123695d-01,
!    &          XYW(3,3)=1.66666666666666601903656896865868d-01,
!    &          XYW(4,3)=3.33333333333333314829616256247391d-01)
      DOUBLE PRECISION HELP
           
C
      DOUBLE PRECISION DNRM2
      EXTERNAL FLUXF
C
      INTEGER I,J,IFAIL
      LOGICAL VERBOSE
      PARAMETER(VERBOSE=.FALSE.)
!     PARAMETER(VERBOSE=.TRUE.)
      XYW(1,1)=1.66666666666666657414808128123695d-01 ! area coordinate
      XYW(2,1)=1.66666666666666657414808128123695d-01 ! area coordinate
      XYW(3,1)=6.66666666666666740681534975010436d-01 ! area coordinate
      XYW(4,1)=3.33333333333333314829616256247391d-01 ! weight
      XYW(1,2)=1.66666666666666657414808128123695d-01 ! area coordinate
      XYW(2,2)=6.66666666666666740681534975010436d-01 ! area coordinate
      XYW(3,2)=1.66666666666666629659232512494782d-01 ! area coordinate
      XYW(4,2)=3.33333333333333314829616256247391d-01 ! weight
      XYW(1,3)=6.66666666666666740681534975010436d-01 ! area coordinate
      XYW(2,3)=1.66666666666666657414808128123695d-01 ! area coordinate
      XYW(3,3)=1.66666666666666601903656896865868d-01 ! area coordinate
      XYW(4,3)=3.33333333333333314829616256247391d-01 ! weight
C
      NM1 = (NOFVERT-1)
C
      DO I = 1,NM1
      CALL DINIT(NDOF,ZERO,ZROE(1,NOFVERT),1)
      CALL DINIT(NDIM,ZERO,XYZDOT(1,NOFVERT),1)
         DO J = 1,NM1
            CALL DAXPY(NDOF,XYW(I,J),ZROE(1,J),1,ZROE(1,NOFVERT),1)
            CALL DAXPY(NDIM,XYW(I,J),XYZDOT(1,J),1,XYZDOT(1,NOFVERT),1) ! compute average grid velocity
         ENDDO
         CALL FLUXF(NDIM,FACN,ZROE(1,NOFVERT),XYZDOT(1,NOFVERT),
     &              FLXN(1,I)) ! compute flux in vertex I of the face
      ENDDO
C
      IF(VERBOSE)THEN
         CALL R8Mat_Print('General',' ',NDOF,NOFVERT,FLXN,
     +            NDOF,'Nodal fluxes (vertex + gravity) ',IFAIL)
      ENDIF
C
      CALL DINIT(NDOF,ZERO,FLXN(1,NOFVERT),1)
      DO I = 1,NM1
         CALL DAXPY(NDOF,XYW(4,I),FLXN(1,I),1,FLXN(1,NOFVERT),1) ! use weights
      ENDDO
C
      IF(VERBOSE)THEN
         CALL R8Mat_Print('General',' ',NDOF,NOFVERT,ZROE,
     +            NDOF,'Nodal values incl. average ',IFAIL)
      ENDIF
C
!
      IF(VERBOSE)THEN
         HELP = DNRM2(NDIM,FACN,1)
         WRITE(6,*)'Face area is ',HELP
         CALL R8Mat_Print('General',' ',NDOF,NOFVERT,FLXN,
     +            NDOF,'Nodal fluxes in subr. Quadrature ',IFAIL)
      ENDIF
C
      RETURN
      END
