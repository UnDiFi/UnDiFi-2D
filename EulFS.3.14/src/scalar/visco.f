      SUBROUTINE VISCO(IELEM,VCZ,NODRES,DT,NOFVAR,VCN,NDIM,NOFVERT,
     +                 VOLUME,STIFD,EPSILON,DUMMY,PICARD)
C
C
C     Subroutine to compute the diffusion term over a triangle/
C     tetrahedron for scalar problems. The diffusion coefficient
C     EPSILON is taken to be constant
C
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION EPSILON,VOLUME,DUMMY
      INTEGER IELEM,NDIM,NOFVERT
      LOGICAL PICARD
C
C     EPSILON 
C             On entry, EPSILON specifies the diffusion coefficient. 
C             Unchanged on exit.
C
C     VOLUME 
C             On entry, VOLUME specifies the area/volume of the current
C             element.  Unchanged on exit.
C
C     IELEM 
C             On entry, IELEM specifies the number of the current
C             element.  Unchanged on exit.
C
C     NDIM 
C             On entry, NDIM specifies the dimension of the space.
C             Unchanged on exit.
C
C     NOFVERT 
C             On entry, NOFVERT specifies the number of vertices
C             of the current element.  Unchanged on exit.
C
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION STIFD(NOFVERT,NOFVERT),DT(NOFVERT),
     +                 NODRES(NOFVERT),VCN(NDIM,NOFVERT),VCZ(NOFVERT)
C
C     STIFD    On exit, the array STIFD is overwritten with 
C             the element diffusion matrix
C
C     DT      On entry, DT contains the of the local timestep
C             in the nodes of the current element. 
C             On exit, the timestep due to the viscous terms is
C             added to.
C
C     NODRES  On entry, NODRES contains the residual
C             in the nodes of the current element.
C             On exit, the nodal residual due to the viscous terms is
C             added to.
C
C     VCN     On entry, VCN contains the components of the
C             inward normals to the edges/faces of the current element.
C             Unchanged on exit.
C
C     VCZ       On entry, VCZ contains the values of the dependent
C             variable in te nodes of the current element.
C             Unchanged on exit.
C
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION SUM,TEMPB
      INTEGER I,J
      DOUBLE PRECISION TMPIJ(4,4)
C     ..
C     .. External Functions ..
      DOUBLE PRECISION DDOT
      EXTERNAL DDOT
C     ..
C     ..
      TEMPB = EPSILON/ (NDIM*NDIM*VOLUME)
C
      DO 1 J = 1,NOFVERT
          DO 1 I = 1,NOFVERT
              IF (J.LE.I) THEN
                  TMPIJ(I,J) = TEMPB*DDOT(NDIM,VCN(1,I),1,VCN(1,J),1)

              ELSE
                  TMPIJ(I,J) = TMPIJ(J,I)
              ENDIF

    1 CONTINUE
C
      DO 2 I = 1,NOFVERT
          SUM = 0.D0
          DO 3 J = 1,NOFVERT
              SUM = SUM + TMPIJ(I,J)*VCZ(J)
              IF(PICARD)STIFD(I,J)=STIFD(I,J)+TMPIJ(I,J)
    3     CONTINUE
          NODRES(I) = NODRES(I) - SUM


          DT(I) = DT(I) + TMPIJ(I,I)
    2 CONTINUE
C
      RETURN

      END
