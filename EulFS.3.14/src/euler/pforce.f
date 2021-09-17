!> \brief \b PFORCE
C
*>    \par Purpose:
*>
*     ====================== 
CCCCCCC*>    \verbatim
*>    PFORCE computes the pressure force acting on a
*>    boundary face as the aritmetic mean of the face pressure
*>    values, ie. assuming linear variation of pressure;
*>    while this is o.k. for the INcompressible flow eqns., for
*>    COmpressible flows one should compute
*>    the EXACT integral of pressure which is a quadratic function
CCCCCCC*>    \endverbatim
*
C
!> @param[in] ICLR colour of the patch the boundary face belongs to 
!> @param[in] IVERT vertex (in local numbering) facing the boundary face   
!> @param[in] VCN cartesian components of the scaled inward normals of the current cell
!> @param[in] NDIM dimension of the space
!> @param[in] VCZ nofvar values of the dependent variable in the nofvert vertices of the current cell
!> @param[in] NOFVAR number of degrees of freedom
!> @param[in] NOFVERT number of vertices of the current simplicial element
!> @param[in] PRESSURE subroutine used to compute pressure from the set of dependent variables
!> @return the cartesian components of the pressure force acting on the current boundary face are stored in PRESF(1:d,iclr)
C
!> \author $Author: abonfi $
!> \version $Revision: 1.9 $
!> \date $Date: 2013/08/21 10:55:07 $
!> \see subroutine PFORCESimpson
C
      SUBROUTINE PFORCE(ICLR,IVERT,VCN,NDIM,VCZ,NOFVAR,NOFVERT,PRESSURE)
C
C     $Id: pforce.f,v 1.9 2013/08/21 10:55:07 abonfi Exp $
C
      IMPLICIT NONE
C
      INCLUDE 'constants.h'
      INCLUDE 'paramt.h'
      INCLUDE 'bnd.h'
      INCLUDE 'bnd.com'
C
      INTEGER ICLR,IVERT,NDIM,NOFVERT,NOFVAR
C
C
      DOUBLE PRECISION VCZ(NOFVAR,NOFVERT),VCN(NDIM,NOFVERT)
C
      DOUBLE PRECISION P
      INTEGER I,IV
C
      DOUBLE PRECISION PRESSURE
      INTEGER ICYCL
      EXTERNAL ICYCL,PRESSURE
C
C  .. Loop over the nodes of the face
C
      P = ZERO
      DO 1 I = 1,NOFVERT - 1
          IV = ICYCL(IVERT+I,NOFVERT)
          P = P + PRESSURE(NDIM,VCZ(1,IV))
    1 CONTINUE
C
      PRESF(1,ICLR) = PRESF(1,ICLR) - P*VCN(1,IVERT)/NDIM
      PRESF(2,ICLR) = PRESF(2,ICLR) - P*VCN(2,IVERT)/NDIM
      IF(NDIM.EQ.3)PRESF(3,ICLR) = PRESF(3,ICLR) - P*VCN(3,IVERT)/NDIM
C
      RETURN

      END
!> \brief \b PFORCESimpson
C
!>    \par Purpose:
!>
!>    This routine computes the pressure force acting on a
!>    boundary face using Simpson's rule which is exact
!>    when pressure is a quadratic function of the dependent variables
C     ====================== 
C
!> @param[in] ICLR colour of the patch the boundary face belongs to 
!> @param[in] IVERT vertex (in local numbering) facing the boundary face   
!> @param[in] VCN cartesian components of the scaled inward normals of the current cell
!> @param[in] NDIM dimension of the space
!> @param[in] VCZ nofvar values of the dependent variable in the nofvert vertices of the current cell
!> @param[in] NOFVAR number of degrees of freedom
!> @param[in] NOFVERT number of vertices of the current simplicial element
!> @param[in] PRESSURE subroutine used to compute pressure from the set of dependent variables
!> @return the cartesian components of the pressure force acting on the current boundary face are stored in PRESF(1:d,iclr)
C
!> \warning The 3D version is un-implemented yet
!> \author $Author: abonfi $
!> \version $Revision: 1.9 $
!> \date $Date: 2013/08/21 10:55:07 $
!> \see subroutine PFORCE
C
      SUBROUTINE PFORCESimpson(ICLR,IVERT,VCN,NDIM,VCZ,NOFVAR,NOFVERT,
     &PRESSURE)
C
C     $Id: pforce.f,v 1.9 2013/08/21 10:55:07 abonfi Exp $
C
      IMPLICIT NONE
C
      INCLUDE 'constants.h'
      INCLUDE 'paramt.h'
      INCLUDE 'bnd.h'
      INCLUDE 'bnd.com'
C
      INTEGER ICLR,IVERT,NDIM,NOFVERT,NOFVAR
C
C
C
      DOUBLE PRECISION VCZ(NOFVAR,NOFVERT),VCN(NDIM,NOFVERT)
      DOUBLE PRECISION ZG(MAXNOFVAR)
C
      DOUBLE PRECISION P
      INTEGER I,IV
C
      DOUBLE PRECISION PRESSURE
      INTEGER ICYCL
      EXTERNAL ICYCL,PRESSURE
C
C  .. Loop over the nodes of the face
C
      IF(NDIM.NE.2)THEN
         STOP 'Unimplemented feature in 3D'
      ENDIF
      CALL DINIT(NOFVAR,ZERO,ZG,1)
      P = ZERO
      DO 1 I = 1,NOFVERT - 1
          IV = ICYCL(IVERT+I,NOFVERT)
          P = P + PRESSURE(NDIM,VCZ(1,IV))
          CALL DAXPY(NOFVAR,ONE,VCZ(1,IV),1,ZG,1) ! compute average Z
    1 CONTINUE
      CALL DSCAL(NOFVAR,ONE/REAL(NDIM),ZG,1)
      P = (P+4.d0*PRESSURE(NDIM,ZG))/6.d0 ! use Simpson's rule
C
      PRESF(1,ICLR) = PRESF(1,ICLR) - P*VCN(1,IVERT)
      PRESF(2,ICLR) = PRESF(2,ICLR) - P*VCN(2,IVERT)
      IF(NDIM.EQ.3)PRESF(3,ICLR) = PRESF(3,ICLR) - P*VCN(3,IVERT)
C
      RETURN

      END
