!> \par Purpose
!>
!> retrieve data for cell \c IELEM from global arrays
!>
!> @param[in] IELEM is the current cell in local (per processor) numbering
!> @param[in] NELEM nof cells on the current PE
!> @param[in] ICELNOD Cell to node pointer: \c ICELNOD(i,j) gives the global node number of the i-th vertex of the j-th element
!> @param[in] ICELFAC Cell to face pointer: \c ICELFAC(i,j} gives the global face number of the face opposite the i-th vertex of the j-th element
!> @param[in] VOL area/volume of the simplicial elements (triangles,tetrahedra)
!> @param[in] ZROE stores the dependent \c NOFVAR variables within all meshpoints; ZROE(1:NOFVAR,1,*) stores the values at time level \c n+1,k ; \c ZROE(1:NOFVAR,1,*) those at tme level n and \c ZROE(1:NOFVAR,3,*) those at time level \c n-1; could actually feed any other vector using the same layout, such as, e.g. the nodal coordinates; however, there must be room enough in the unsteady case
!> @param[in] FACNOR Cartesian components of the normals to a face, multiplied by the face area
!> @param[in] XYZDOT Cartesian components of the nodal grid velocities
!> @param[in] NDIM dimension of the space
!> @param[in] NOFVERT nof boundary faces
!> @param[in] NOFVAR nof dofs within each meshpoint
!> @param[in] LDA second leading dimension of \c ZROE, should equal NPOIN+NGHOST+NPNOD
!> @param[out] ICN returns the vertices of the current simplex in 0-based indexing; to be used to insert values into PETSc vecs and/or mats 
!> @param[out] VCZ returns the \c NOFVAR dofs of the NOFVERT vertices of cell \c IELEM; VCZ(1:NOFVAR,*,1) stores the values at time level \c n+1,k ; \c VCZ(1:NOFVAR,*,2) those at tme level n and \c VCZ(1:NOFVAR,*,3) those at time level \c n-1
!> @param[out] VCN returns the \c NDIM Cartesian components of the NOFVERT faces of cell \c IELEM
!> @param[out] VCB returns the \c NDIM Cartesian components of the nodal grid velocities at time \c n+1/2 of the NOFVERT vertices of cell \c IELEM; only if \c LALE is \c .TRUE.
!> @param[out] VOLUME returns the volume of simplex \c IELEM 
!> \author $Author: abonfi $
!> \version $Revision: 1.10 $
!> \date $Date: 2020/03/28 09:45:57 $
!> \warning on periodic meshes we access the 2nd copy of the cell to node pointer, i.e. the original addressing
      SUBROUTINE CELPTR(IELEM,NELEM,ICELNOD,ICELFAC,VOL,ZROE,FACNOR,
     2                  XYZDOT,NDIM,NOFVERT,NOFVAR,LDA,ICN,VCZ,VCN,VCB,
     3                  VOLUME)
C
C     $Id: celptr.f,v 1.10 2020/03/28 09:45:57 abonfi Exp $
C
      IMPLICIT NONE
C
C
C     .. Parameters ..
      INCLUDE 'paramt.h'
      INCLUDE 'constants.h'
      INCLUDE 'time.h'
      INCLUDE 'time.com'
C     ..
C     .. Scalar Arguments ..
C     Input:
      INTEGER IELEM,NELEM,NDIM,NOFVAR,NOFVERT,LDA
C     Output:
      DOUBLE PRECISION VOLUME(*)
C     ..
C     .. Array Arguments ..
C     Input:
      DOUBLE PRECISION FACNOR(NDIM,*),ZROE(NOFVAR,LDA,*),XYZDOT(NDIM,*)
      INTEGER ICELFAC(NOFVERT,*),ICELNOD(NOFVERT,*)
C     Output:
      DOUBLE PRECISION VCN(NDIM,NOFVERT),VCZ(NOFVAR,NOFVERT,*),
     2                 VCB(NDIM,NOFVERT),VOL(*)
      INTEGER ICN(NOFVERT)
C     ..
C     .. Local Scalars ..
      INTEGER IDIM,IFACE,IPOIN,IVAR,IVERT,ILEV
C     ..
C     .. Intrinsic Functions ..
C
      INTRINSIC ABS,SIGN
C     ..
      DO 10 IVERT = 1,NOFVERT
C
C     in periodic grids, IPOIN is
C     the original addressing
C     (note that we access the 2nd copy of the
C     cell to node pointer)
C
          IPOIN = ICELNOD(IVERT,IELEM+NELEM)
C
C     zero indexing for PETSc
C     in periodic grids, ICN will keep 
C     the re-mapped addressing
C
          ICN(IVERT) = ICELNOD(IVERT,IELEM) - 1
C
          IFACE = ICELFAC(IVERT,IELEM)
C
          DO 8 ILEV = 1, NTIMLEVS
             DO 8 IVAR = 1,NOFVAR
                 VCZ(IVAR,IVERT,ILEV) = ZROE(IVAR,IPOIN,ILEV)
    8     CONTINUE
C
          IF(LALE)THEN
              DO 9 IDIM = 1,NDIM
                  VCB(IDIM,IVERT) = XYZDOT(IDIM,IPOIN)
    9         CONTINUE
          ENDIF
C
          IF (IFACE.GT.0) THEN
              DO 6 IDIM = 1,NDIM
                  VCN(IDIM,IVERT) = FACNOR(IDIM,IFACE)
    6         CONTINUE

          ELSE
              IFACE = -IFACE
              DO 7 IDIM = 1,NDIM
                  VCN(IDIM,IVERT) = -FACNOR(IDIM,IFACE)
    7         CONTINUE
          ENDIF
C
   10 CONTINUE
C
      VOLUME(1) = VOL(IELEM)
      IF(LTIME.AND.LALE.AND.(NTIMLEVS.EQ.3))THEN
         VOLUME(2) = VOL(IELEM+NELEM) ! time level n+1
         VOLUME(3) = VOL(IELEM+2*NELEM) ! time level n
         VOLUME(4) = VOL(IELEM+3*NELEM) ! time level n-1
      ENDIF
C
C     CALL R8Mat_Print('General',' ',NDIM,NOFVERT,VCN(1,1),NDIM,
C    +            'VCN',IFAce)
!     CALL R8Mat_Print('General',' ',NDIM,NOFVERT,VCB(1,1),NDIM,
!    +            'VCB',IFAce)
C     CALL R8Mat_Print('General',' ',NDIM,NOFVERT,VCP(1,1),NDIM,
C    +            'VCP',IFAce)
!     do idim = 1,ntimlevs
!      write(6,*)'Time lev within celptr is = ',idim
!     CALL R8Mat_Print('General',' ',NOFVAR,NOFVERT,VCZ(1,1,idim),
!    +            NOFVAR,'VCZ',IFAce)
!     enddo
C     PAUSE
C
      RETURN

      END
