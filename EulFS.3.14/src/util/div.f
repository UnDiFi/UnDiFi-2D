!> \par Purpose
!>
!> Compute the volume mutiplied divergence of a vector field.
!>
!> @param[in] NDIM dimension of the space
!> @param[in] NOFVERT nof vertices for a cell = NDIM+1
!> @param[in] VCN the cartesian components of the scaled inward normal to the face
!> @param[in] V the cartesian components of the grid velocity \f$ \mathbf{b}\f$ or any other vector field
!>
!>
!> \f{eqnarray*}{
!> \nabla \cdot \mathbf{b} \; |T_e| &=& - \sum_{\ell=1}^{d} \frac{1}{d+1} \left(\mathbf{b}_i+\mathbf{b}_j+\mathbf{b}_k\right) \cdot \mathbf{n}_{\ell} \\
!> &=& \frac{1}{d} \sum_{j=1}^{d+1} \mathbf{n}_j \cdot \mathbf{b}_j
!>
!> \f}
C
!> \author $Author: abonfi $
!> \version $Revision: 1.2 $
!> \date $Date: 2014/04/15 10:16:20 $
C
      DOUBLE PRECISION FUNCTION DIV(NDIM,NOFVERT,VCN,V)
      IMPLICIT NONE
C
C     computes the divergence of the vector field V
C     (defined in the vertices of a simplicial element)
C     VCN are the inward scaled normals to the vertices
C
C     $Id: div.f,v 1.2 2014/04/15 10:16:20 abonfi Exp $
C
      INCLUDE "constants.h"
C
C     VCN are the inward face normals, scaled by their measure
C     V is the vector field defined in the vertices of the mesh
C
      INTEGER NDIM,NOFVERT
      INTEGER I
      DOUBLE PRECISION VCN(NDIM,NOFVERT),V(NDIM,NOFVERT)
      DOUBLE PRECISION DDOT
      EXTERNAL DDOT
C
      DIV = ZERO
      DO I = 1,NOFVERT
         DIV = DIV + DDOT(NDIM,V(1,I),1,VCN(1,I),1)
      ENDDO
      DIV=DIV/REAL(NDIM)
      RETURN
      END
