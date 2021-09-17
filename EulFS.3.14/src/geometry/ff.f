!> \par Purpose
!>
!> Compute the faces of the mesh; one face normal (its \c NDIM cartesian components) is stored for each edge/face shared by two elements i.e. the normals of the inner (shared) faces are NOT duplicated; the
!> pointer \c ICELFAC(i,j) is used to address the faces of a cell and know whether it is pointing inside or outside
!>
!> @param[in] ICELNOD Cell to node pointer: \c ICELNOD(i,j) gives the global node number of the i-th vertex of the j-th element
!> @param[in] ICELCEL Cell to cell pointer: \c ICELCEL(i,j) gives the element number that shares the face opposite the i-th vertex of the j-th element if 0 or > \c NELEM, that face is a boundary face
!> @param[in] ICELFAC Cell to face pointer: \c ICELFAC(i,j) gives the global face number of the face opposite the i-th vertex of the j-th element
!> @param[in] NOFVERT nof boundary faces
!> @param[in] NELEM nof boundary faces
!> @param[in] CORG Cartesian coordinates of the meshpoints
!> @param[in] NDIM dimension of the space
!> @param[in] NPOIN nof interior nodes in the mesh
!> @param[in,out] FACNOR Cartesian components of the normals to a face, multiplied by the face area
!> @param[in] NFACE nof faces in the mesh
!> @param[in] NBFAC nof boundary faces
!> @param[in] NBINT nof inter-processor faces
!>
!> \author $Author: abonfi $
!> \version $Revision: 1.11 $
!> \date $Date: 2020/03/28 09:45:59 $
      SUBROUTINE FF(ICELNOD,ICELCEL,ICELFAC,NOFVERT,NELEM,CORG,NDIM,
     +              NPOIN,FACNOR,NFACE,NBFAC,NBINT)
C
      IMPLICIT NONE
C
      INCLUDE 'constants.h'
      INCLUDE 'io.com'
C
C     This routines finds all faces of a 3D mesh
C
C     .. Scalar Arguments ..
      INTEGER NBFAC,NDIM,NELEM,NFACE,NOFVERT,NPOIN,NBINT
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION CORG(NDIM,NPOIN),FACNOR(NDIM,NFACE)
      INTEGER ICELCEL(NOFVERT,NELEM),ICELFAC(NOFVERT,NELEM),
     +        ICELNOD(NOFVERT,NELEM)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION DUMMY,S,W
      INTEGER I,IELEM,IFACE,IFAIL,IFREQ,II,IV,J,JV,JVERT,K,L,LL,N,
     +        NBFCHK,NEIGHB,NERR,IOPT
      CHARACTER*72 ERRMSG
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION WKSP(3),XJI(3),XJK(3),XJL(3)
C     ..
C     .. External Functions ..
      DOUBLE PRECISION DDOT,DNRM2
      INTEGER ICYCL,I1MACH
      EXTERNAL DDOT,DNRM2,ICYCL,I1MACH
C     ..
C     .. External Subroutines ..
      EXTERNAL CROSS_PROD,DAXPY,DSCAL
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,IABS,INT,ISIGN,MAX0,SIGN
C     ..
C     .. Data statements ..
C     ..
C
      IFACE = 0
      NBFCHK = 0
      IFREQ = MAX0(1,INT(NFACE)/20)
      WRITE (NOUT,FMT=113)
C
C Loop over the cells to find the face normals
C
      DO 10 IELEM = 1,NELEM
C
C Loop over the vertices of the cell
C
          DO 12 JVERT = 1,NOFVERT
C
              NEIGHB = ICELCEL(JVERT,IELEM)
              I = ICELNOD(JVERT,IELEM)
              J = ICELNOD(ICYCL(JVERT+1,NOFVERT),IELEM)
              K = ICELNOD(ICYCL(JVERT+2,NOFVERT),IELEM)
              L = ICELNOD(ICYCL(JVERT+3,NOFVERT),IELEM)
C
C If the neighbouring cell opposite to node "jvert" has an index
C greater than "IELEM" or the face opposite node "jvert" is a
C boundary face( NEIGHB <= 0 ), then : COMPUTE the normal
C
              IF (NEIGHB.GT.IELEM .OR. NEIGHB.LE.0) THEN
C
                  IFACE = IFACE + 1
C
                  IF ((IFACE/IFREQ)*IFREQ.EQ.IFACE) WRITE (NOUT,FMT=111)
C
C test that the space =NFACE allocated for storing the normals
C is enough
C
                  IF (IFACE.GT.NFACE) THEN
                      WRITE (I1MACH(4),FMT=405)
                      WRITE (ERRMSG(1:43),FMT=406)
                      NERR = 9
                      IOPT = 1
                      CALL SETERR(ERRMSG,43,NERR,IOPT)
                  ENDIF
C
                  ICELFAC(JVERT,IELEM) = IFACE
C
                  GOTO (200,300) NDIM - 1
C----2D case
  200             CONTINUE
                  DO 201 II = 1,NDIM
                      XJI(II) = CORG(II,I) - CORG(II,J)
                      XJK(II) = CORG(II,K) - CORG(II,J)
  201             CONTINUE
                  FACNOR(1,IFACE) = -XJK(2)
                  FACNOR(2,IFACE) = XJK(1)
                  GOTO 400
C---3D case
  300             CONTINUE
                  DO 301 II = 1,NDIM
                      XJI(II) = CORG(II,I) - CORG(II,J)
                      XJK(II) = CORG(II,K) - CORG(II,J)
                      XJL(II) = CORG(II,L) - CORG(II,J)
  301             CONTINUE
C
C Compute the "scaled" normal
C
                  CALL CROSS_PROD(XJK,XJL,FACNOR(1,IFACE))
                  CALL DSCAL(NDIM,HALF,FACNOR(1,IFACE),1)
C
  400             CONTINUE
C
C Checks that it is NOT 0.0E+00
C
                  IF (DNRM2(NDIM,FACNOR(1,IFACE),1).LE.1.E-10) THEN
                      WRITE (NOUT,FMT=210) IFACE,J,K,L
                  ENDIF
C
C Check whether the normal is inward or outward
C
                  W = DDOT(NDIM,XJI,1,FACNOR(1,IFACE),1)

                  DUMMY = SIGN(ONE,W)
                  ICELFAC(JVERT,IELEM) = INT(DUMMY)*ICELFAC(JVERT,IELEM)
C
C The neighbour of a boundary element is flagged 0 or with an
C element number > NELEM
C If the face is a boundary face, we want the normal to be pointing
C INSIDE the computational domain
C
                  IF (NEIGHB.LE.0 .OR. NEIGHB.GT.NELEM) THEN
                     NBFCHK = NBFCHK + 1
                     CALL DSCAL(NDIM,DUMMY,FACNOR(1,IFACE),1)
                     ICELFAC(JVERT,IELEM) = IABS(ICELFAC(JVERT,IELEM))
                  ENDIF
C
C
C If the neighbouring cell opposite to node "jvert" has an index
C smaller than "icell" : the normal already exist !
C
              ELSEIF (NEIGHB.LT.IELEM) THEN
C
C Find the only node of the neighbouring cell which does not belong
C to the present.
C Loop over the verices of the neighbouring cell
                  DO 31 JV = 1,NOFVERT
                      LL = ICELNOD(JV,NEIGHB)
C Loop over the verices of the current cell and get the only one
C which does not belong to the neighbouring cell
                      DO 33 IV = 1,NOFVERT
                          IF (ICELNOD(IV,IELEM).EQ.LL) GOTO 31
   33                 CONTINUE
                      ICELFAC(JVERT,IELEM) = -ICELFAC(JV,NEIGHB)
                      GOTO 12

   31             CONTINUE
                  WRITE (*,FMT=*) ' Smthg. has gone wrong in cell no. ',
     +              IELEM
                  WRITE (*,FMT=*) (ICELNOD(JV,IELEM),JV=1,NOFVERT)
              ENDIF
C
   12     CONTINUE
C
   10 CONTINUE
C
      WRITE (NOUT,FMT=115) IFACE
C
C     test number of faces
C
      IF (IFACE.NE.NFACE) THEN
          WRITE (I1MACH(4),FMT=205) NFACE,IFACE
          WRITE (ERRMSG(1:40),FMT=206)
          NERR = 10
          IOPT = 1
!         CALL SETERR(ERRMSG,40,NERR,IOPT)
      ENDIF
C
C     make some checks
C
C     Loop over all elements to check that
C     the scaled normals sum up to 0.E+00 for EACH element
C
      NERR = 0
      DO 20 IELEM = 1,NELEM
C
          DO 7 II = 1,NDIM
    7     WKSP(II) = ZERO
C
          DO 5 JV = 1,NOFVERT
              N = ICELFAC(JV,IELEM)
              S = ISIGN(1,N)
              N = IABS(N)
              CALL DAXPY(NDIM,S,FACNOR(1,N),1,WKSP,1)
    5     CONTINUE
C
C Check whether the scaled normals sum up to 0.0E+00
C
          S = DNRM2(NDIM,WKSP,1)
          IF (ABS(S).GT.1.E-6) THEN
              NERR = NERR + 1 
              WRITE (NOUT,FMT=999) IELEM,S
              DO JV = 1,NOFVERT 
                 N = ICELFAC(JV,IELEM)
                 S = ISIGN(1,N)
                 N = IABS(N)
                 WRITE(NOUT,*)N,S
                 WRITE(NOUT,*)(FACNOR(II,N),II=1,NDIM)
              ENDDO
              WRITE(NOUT,*)(WKSP(II),II=1,NDIM)
          ENDIF
C
   20 CONTINUE ! end loop over cells
       IF(NERR.NE.0)THEN
          CALL I4Mat_Print('General',' ',NOFVERT,NELEM,ICELFAC,
     +    NOFVERT,'ICELFAC array ',IFAIL)
          IOPT = 1
          WRITE (ERRMSG(1:43),FMT=208)
          CALL SETERR(ERRMSG,43,NERR,IOPT)
       ENDIF
C
C     test boundary faces
C
      IF ((NBFAC+NBINT).NE.NBFCHK) THEN
          WRITE (I1MACH(4),FMT=207) NBFAC,NBINT,NBFCHK
          WRITE (ERRMSG(1:43),FMT=208)
          NERR = 11
          IOPT = 1
          CALL SETERR(ERRMSG,43,NERR,IOPT)
      ENDIF
C
      RETURN
C
C     I/O FORMATS
C
  111 FORMAT ('.',$)
  113 FORMAT (10X,'COMPUTING NORMALS ',$)
  115 FORMAT (I8,/)
  207 FORMAT (5X,10 ('*'),' WARNING ',10 ('*'),/,5X,I7,
     +       ' BOUNDARY FACES(EDGES) WERE EXPECTED',/,5X,I7,
     +       ' INTER-PROC FACES(EDGES) HAVE BEEN FOUND',/,5X,I7,
     +       ' BOUNDARY FACES(EDGES) HAVE BEEN FOUND',/)
  208 FORMAT ('FF -- CHECK ON BOUNDARY FACES(EDGES) FAILED')
  205 FORMAT (5X,10 ('*'),' ERROR ',10 ('*'),/,5X,I7,
     +       ' FACES(EDGES) WERE EXPECTED',/,5X,I7,
     +       ' FACES(EDGES) HAVE BEEN FOUND',/)
  206 FORMAT ('FF -- TOO MUCH ROOM TO STORE THE NORMALS') 
  210 FORMAT (3X,'WARNING !! ZERO NORMAL FOR FACE ',I5,' NODES ',
     +       3 (3X,I6))
  405 FORMAT (/,25X,10 ('*'),' ERROR ',10 ('*'),/,5X,
     +       'THE ESTIMATED NUMBER OF FACES IS INSUFFICIENT',/,5X,
     +       'CHECK THE CONNECTIVITY OF THE MESH')
  406 FORMAT (' FF -- NOT ENOUGH ROOM TO STORE THE NORMALS')
  999 FORMAT (5X,'ERROR : Scaled normals in cell ',I6,' sum up to ',
     +       E9.4)

      END
