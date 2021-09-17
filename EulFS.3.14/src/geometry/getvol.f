!> \par Purpose
!>
!> compute the area (2D) volume (3D) of a simplicial element
!> A formula for computing the volume of a simplicial element \c T
!> (triangle or tetrahedron) is the following:
!> \f[
!> |T| = \frac{1}{d^2} \sum_{j=1}^{d+1} \mathbf{x}_j \cdot \mathbf{n}_j.
!> \f]
!>
!>
!> Indeeed, using Gauss theorem one gets:
!>
!> \f{eqnarray*}{
!> d \, |T| &=& \int_{T} \nabla \cdot \mathbf{x} \, \mathrm{d}V =
!>            - \oint_{\partial T} \mathbf{x} \cdot d \mathbf{n}
!>           =  - \underline{\sum_{j=1}^{d+1}
!>                \frac{1}{d} \left( \sum_{\ell \neq j}
!>                \mathbf{x}_{\ell} \right)} \cdot \mathbf{n}_j \\\
!>          &=& - \frac{1}{d} \sum_{j=1}^{d+1} \left[ \mathbf{x}_j \cdot
!>                \underbrace{ \left( \sum_{\ell \neq j}
!>                \mathbf{n}_{\ell} \right) }_{ - \mathbf{n}_j } \right]
!>          = \frac{1}{d} \sum_{j=1}^{d+1} \mathbf{x}_j \cdot \mathbf{n}_j.
!> \f}
!> The underlined term in the above equation
!> is the arithmetic average of the values of the dependent variable (\c x)
!> over the face opposite vertex \c j.
!> It equals, without involving any
!> approximation, the surface integral.
!>
!> This is a consequence of the fact that dependent variable
!> is linear in \f$x\f$ and that
!> the faces over which integration is carried
!> are planar (in 2D sides they are are straight segments).
!>
!> @param[in] ICELNOD Cell to node pointer: \c ICELNOD(i,j) gives the global node number of the i-th vertex of the j-th element
!> @param[in] ICELFAC Cell to face pointer: \c ICELFAC(i,j} gives the global face number of the face opposite the i-th vertex of the ielem-th cell; If ICELFAC(i,ielem) > 0 the normal vector ( which is stored in VFACNOR(1:NDIM,IABS(ICELFAC(i,ielem))) ) points inside the ielem-th cell; outside if ICELFAC(i,ielem) < 0;  the storage used for the face (edge) normals is as follows: a real array FACNOR(1:NDIM,1:NFACE) stores all the face(edge) normals of the mesh and an integer pointer: -NFACE <= ICELFAC(1:NOFVERT,1:NELEM) <= NFACE gives, with its absolute value, the normal opposite a given vertex of a given element
!> @param[in] NOFVERT nof boundary faces
!> @param[in] NELEM nof boundary faces
!> @param[in] FACNOR Cartesian components of the normals to a face, multiplied by the face area; Face normals : FACNOR(i,iface) gives the i-th cartesian component of the iface-th face. This is the vector normal to a triangular face (edge segment in 2D) scaled by the face area (edge length in 2D)
!> @param[in] NFACE nof faces in the mesh
!> @param[in] VCORG Cartesian coordinates of the meshpoints
!> @param[in] NDIM dimension of the space
!> @param[in] NPOIN nof interior nodes in the mesh
!> @param[out] VOL area/volume of the simplicial elements (triangles,tetrahedra)
!> @param[in] TIME the time when the volumes are computed
!>
!> \author $Author: abonfi $
!> \version $Revision: 1.8 $
!> \date $Date: 2020/03/28 09:46:02 $
      SUBROUTINE CMPVOL(ICELNOD,ICELFAC,NOFVERT,NELEM,FACNOR,NFACE,
     1VCORG,NDIM,NPOIN,VOL,TIME)
C
      IMPLICIT NONE
C
      INCLUDE'constants.h'
      INCLUDE'io.com'
C
      INTEGER NOFVERT,NELEM,NDIM,NFACE,NPOIN
      DOUBLE PRECISION TIME
C
C     .. Array Arguments ..
C
      INTEGER ICELNOD(NOFVERT,NELEM),ICELFAC(NOFVERT,NELEM)
      DOUBLE PRECISION FACNOR(NDIM,NFACE),VCORG(NDIM,NPOIN),VOL(NELEM)
C
C     .. Local Scalars ..
C
      INTEGER I,J,K,IFAIL
      INTEGER IELEM,IPOIN,JPOIN,IFACE,IFREQ
      DOUBLE PRECISION SNORM,VOLUME,TVOLUME,VOLMIN,VOLMAX
C
C     .. Local Arrays ..
C
      DOUBLE PRECISION WKSP(12)
C
C     .. External Functions ..
C
      DOUBLE PRECISION DDOT
      EXTERNAL         DDOT
      INTEGER  ICYCL
      EXTERNAL ICYCL
C
C     .. External Subroutines ..
C
      EXTERNAL      DAXPY
C
C     .. Intrinsic Functions ..
C
      INTRINSIC MAX0,DMIN1,DMAX1
C
C     .. Executable Statements ..
C
C     This routine computes the volume of all elements
C     using the fact that int(div(x)) = NDIM*volume
C
C     1/NDIM*SUM_i(x(j)+x(k)+x(l))*n(i)=-int(div(x))
C     (i,j,k,l) cyclic permutation
C
C
      TVOLUME = ZERO
      VOLMIN = 1.d+38
      VOLMAX = ZERO
      IFREQ = MAX0( 1 , NELEM / 20 )
C
      WRITE(NOUT,2000)TIME
2000  FORMAT(//' COMPUTATION OF THE CELL VOLUMES at time t = ',F10.6,
     &/,31('=')/)
      WRITE(NOUT,"(10X,'VOLUMES',$)")
C
      DO 1 IELEM = 1 , NELEM
C
      IF((IELEM/IFREQ)*IFREQ .EQ. IELEM)WRITE(NOUT,111)
         VOLUME = ZERO
         DO 2 I =  1, NOFVERT
            IPOIN = ICELNOD(I,IELEM)
            IFACE = ICELFAC(I,IELEM)
            SNORM = SIGN( 1 , IFACE )
            IFACE = IABS( IFACE )
C
C     Compute the position vector of the
C     center of gravity of the face opposite node I
C
            DO 3 K =  1, NDIM
    3       WKSP(K) = ZERO
            DO 4 J =  1, NDIM
               JPOIN = ICELNOD(ICYCL(I+J,NOFVERT),IELEM)
    4       CALL DAXPY(NDIM,ONE/NDIM,VCORG(1,JPOIN),1,WKSP,1)
            VOLUME = VOLUME - SNORM*DDOT(NDIM,WKSP,1,FACNOR(1,IFACE),1)
    2    CONTINUE
         VOLUME = VOLUME/NDIM
         VOLMIN = DMIN1( VOLMIN , VOLUME )
         VOLMAX = DMAX1( VOLMAX , VOLUME )
         TVOLUME = TVOLUME + VOLUME
         VOL(IELEM) = VOLUME
         IF(VOLUME.LE.ZERO)THEN
               WRITE(NOUT,888)IELEM,VOLUME
               DO 32 I =  1, NOFVERT
                  IPOIN = ICELNOD(I,IELEM)
                  CALL DCOPY(NDIM,VCORG(1,IPOIN),1,WKSP((I-1)*NDIM+1),1)
   32 CONTINUE
                  CALL R8Mat_Print('General',' ',NDIM,NOFVERT,WKSP,
     +            NDIM,'Nodal coordinates ',IFAIL)
                  CALL I4Mat_Print('General',' ',1,NOFVERT,
     +            ICELNOD(1,IELEM),1,'Cell connectivity',IFAIL)
         ENDIF
    1 CONTINUE ! End loop over elements
      WRITE(NOUT,110)NELEM
      WRITE(NOUT,100)TVOLUME,VOLMIN,VOLMAX,VOLMAX/VOLMIN
C
C     I/O FORMATS
C
  100 FORMAT(/10X,'TOTAL VOLUME ......... ',E10.3/, 10X,
     +'MIN.  VOLUME ......... ',E10.3/, 10X,'MAX.  VOLUME ......... ',
     +E10.3/, 10X,'MAX/MIN RATIO ........ ',E10.3/)

  110 FORMAT(I8,/)
  111 FORMAT('.',$)
  888 FORMAT(5X,'WARNING! NEGATIVE OR ZERO VOLUME IN CELL # ',I6,' VOLUM
     1E= ',E14.6)
      RETURN
      END
