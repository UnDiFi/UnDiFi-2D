      SUBROUTINE ARATIO(IVERT,ICN,XYZ,FACNOR,SURFACE,NDIM,NOFVERT,AR,
     +                  HWALL)
C
C     This routine computes cell aspect ratio and wall distance
C     for solid wall, viscous boundary faces.
C     By "wall distance" we mean the distance from the wall face
C     of the node opposite that face
C
C     On entry:
C     --------
C
C     .. IVERT   is the local nodenumber of the vertex opposite
C                the boundary face
C     .. ICN     is the list of nodes of the current boundary element
C     .. XYZ     are the coordinates of the meshpoints
C     .. FACNOR  are the cartesian components of the boundary face
C     .. SURFACE is the surface measure of the boundary face
C     .. NDIM    is the space dimension
C     .. NOFVERT is the number of vertices
C
C     Upon return:
C     -----------
C
C     .. AR      is the aspect ratio of the current boundary element
C     .. HWALL   is the "wall distance" of the current boundary element
C
C
C
C
C
C
C     XYZG is the distance between the barycenter of the boundary
C          face and the node opposite the boundary face
C
C
C
C
C     .. Parameters ..
      DOUBLE PRECISION CNST
      PARAMETER (CNST=3.D0/1.4142135623730950488D0)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION AR,HWALL,SURFACE
      INTEGER IVERT,NDIM,NOFVERT
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION FACNOR(NDIM),XYZ(NDIM,*)
      INTEGER ICN(NOFVERT)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION TEMP
      INTEGER IDIM,IPOIN,IV
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION XYZG(3)
C     ..
C     .. External Functions ..
      DOUBLE PRECISION DDOT
      INTEGER ICYCL
      EXTERNAL DDOT,ICYCL
C     ..
      TEMP = 1.D0/ (NOFVERT-1)
C
      DO 3 IDIM = 1,NDIM
              XYZG(IDIM) = 0.D0
    3 CONTINUE
C
C     compute the vector joining the center of gravity of the
C     boundary face with the opposite vertex
C
      DO 2 IV = 1,NOFVERT - 1
          IPOIN = ICN(ICYCL(IVERT+IV,NOFVERT))
          DO 2 IDIM = 1,NDIM
              XYZG(IDIM) = XYZG(IDIM) + XYZ(IDIM,IPOIN)*TEMP
    2 CONTINUE
C
      IPOIN = ICN(IVERT)
C
      DO 1 IDIM = 1,NDIM
          XYZG(IDIM) = XYZ(IDIM,IPOIN)-XYZG(IDIM)
    1 CONTINUE
C
C     take the projection on the face normal
C
      HWALL = DDOT(NDIM,XYZG,1,FACNOR,1)/SURFACE
C
C     .N.B the following definition of AR only applies in 2D
C          where it equals one for equilateral elements
C
      AR = SURFACE*CNST/HWALL
C
      RETURN

      END
