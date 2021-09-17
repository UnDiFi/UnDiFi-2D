      SUBROUTINE MatGetSizeSeq(ICELNOD,DEGREE,ICELCEL,NODCODE,
     +                    NDIM,NOFVERT,NELEM,NR,NNZR,IOUT,VERBOSE,IFAIL)
C
C     $Id: MatGetSizeSeq.f,v 1.4 2020/02/05 14:52:01 abonfi Exp $
C
C
C
C     .. This routine estimates the workspace required by the
C        stiffness matrix ..
C
C
C ICELNOD -- IN Integer ICELNOD(1:NOFVERT,1:NELEM)
C            Cell to Node pointer : ICELNOD(i,ielem) gives the
C            global node number of the i-th vertex of the ielem-th cell
C
C DEGREE  -- OUT Integer DEGREE(1:NR)
C            DEGREE(i) gives the number of nonzero (block) entries
C            for the i-th row of the stiffness matrix
C
C ICELCEL -- IN Integer ICELCEL(1:NOFVERT,1:NELEM)
C            Cell to Cell pointer : ICELCEL(i,ielem) gives the
C            global number of the element sharing with ielem
C            the face opposite the i-th vertex of the ielem-th cell
C            If ICELCEL(i,ielem) = 0 or ICELCEL(i,ielem) > NELEM
C            the element ielem is a boundary element and the face
C            opposite its i-th vertex is a boundary face
C
      IMPLICIT NONE
C
C     .. Scalar Arguments ..
      INTEGER NDIM,NELEM,NNZR,NOFVERT,NR,IFAIL
      LOGICAL VERBOSE
C     ..
C     .. Array Arguments ..
      INTEGER DEGREE(NR+1),ICELCEL(NOFVERT,NELEM),
     +        ICELNOD(NOFVERT,NELEM),NODCODE(NR)
C     ..
C     .. Local Scalars ..
      INTEGER IELEM,IOUT,IPOIN,IVERT,JELEM,JVERT,KVERT,MAXDEG,
     +        MINDEG,TEMP
C     ..
C     .. External Functions ..
      INTEGER ICYCL
      EXTERNAL ICYCL
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC FLOAT
C     ..
      WRITE(IOUT,FMT=305) 
      IFAIL = 0
C
C     DEGREE(i) gives the number of elements meeting in node i
C
      NNZR = 0
C
C     .. DEGREE(i) is initialized to 1 to account for the diagonal
C        element in the stiffness matrix ..
C
      DO 10 IPOIN = 1,NR + 1
          DEGREE(IPOIN) = 1
   10 CONTINUE
C
C     .. depending on the space dimension ..
C
      GOTO (2000,3000) NDIM - 1
C
 2000 CONTINUE
C
C     .. 2 dimensions ..
C
C     The number of nonzero entries in the IPOIN-th row equals
C     the number of elements meeting in node IPOIN
C      + 1 to account for the diagonal element
C      + 1 if the node is a boundary node
C
      DO 100 IELEM = 1,NELEM
          DO 100 IVERT = 1,NOFVERT
              IPOIN = ICELNOD(IVERT,IELEM)
              DEGREE(IPOIN) = DEGREE(IPOIN) + 1
              NNZR = NNZR + 1
  100 CONTINUE
C
      DO 20 IPOIN = 1,NR
          IF (NODCODE(IPOIN).NE.0) THEN
              DEGREE(IPOIN) = DEGREE(IPOIN) + 1
C
C     .. NNZR needs to be increased
C        to account for boundary points ..
C
              NNZR = NNZR + 1
          ENDIF

   20 CONTINUE
C
C     .. NNZR needs to be increased by NR to account for the
C        diagonal elements of the stiffness matrix
C        and by to account for boundary points ..
C
      NNZR = NNZR + NR
      GOTO 1000
C
 3000 CONTINUE
C
C     .. 3 dimensions ..
C
C     The number of nonzero entries in the IPOIN-th row ia
C     equal to:
C     the number of faces meeting in IPOIN -
C     the number of elements meeting in node IPOIN +
C     2 if an internal node
C     1 if a  boundary node
C     (+1 to account for the diagonal element of the matrix)
C
C     .. counts the number of faces meeting in node IPOIN
C        by looping over all faces of the mesh ..
C
      DO 90 IELEM = 1,NELEM
          DO 80 IVERT = 1,NOFVERT
              JELEM = ICELCEL(IVERT,IELEM)
C
C     .. if the face has already been encountered, skip to the next ..
C
              IF (JELEM.GT.IELEM .OR. JELEM.LE.0) THEN
                  DO 70 JVERT = 1,NDIM
                      KVERT = ICYCL(IVERT+JVERT,NOFVERT)
                      IPOIN = ICELNOD(KVERT,IELEM)
                      DEGREE(IPOIN) = DEGREE(IPOIN) + 1
   70             CONTINUE
              ENDIF

   80     CONTINUE
   90 CONTINUE
C
C     .. substracts the number of elements meeting in node IPOIN ..
C        by looping over all elements of the mesh ..
C
      DO 40 IELEM = 1,NELEM
          DO 40 IVERT = 1,NOFVERT
              IPOIN = ICELNOD(IVERT,IELEM)
              DEGREE(IPOIN) = DEGREE(IPOIN) - 1
   40 CONTINUE
C
C     add 2 for an internal node
C     add 1 for a  boundary node
C
      DO 60 IPOIN = 1,NR
          TEMP = DEGREE(IPOIN)
          IF (NODCODE(IPOIN).EQ.0) THEN
              TEMP = TEMP + 2

          ELSE
              TEMP = TEMP + 1
          ENDIF

          NNZR = NNZR + TEMP
          DEGREE(IPOIN) = TEMP
   60 CONTINUE
C
 1000 CONTINUE
C
      MAXDEG = 0
      MINDEG = 100000
C
      DO 50 IPOIN = 1,NR
          TEMP = DEGREE(IPOIN)
          IF (TEMP.GT.MAXDEG) MAXDEG = TEMP
          IF (TEMP.LT.MINDEG) MINDEG = TEMP
          IF (TEMP.EQ.1) THEN
              WRITE (IOUT,FMT=310) IPOIN,TEMP
              IFAIL = IFAIL + 1 
          ELSEIF (TEMP.LT.1) THEN
              WRITE (IOUT,FMT=310) IPOIN,TEMP
              IFAIL = IFAIL + 1 
          ENDIF

   50 CONTINUE
C
      IF (VERBOSE) WRITE (IOUT,FMT=320) MINDEG,
     +    MAXDEG,NNZR/FLOAT(NR),NNZR
C
  300 FORMAT (15X,'NODE ',I6,' DOES NOT BELONG TO ANY ELEMENT')
  305 FORMAT (/,/,5X,'SUBROUINE MatGetSizeSeq ',/,/)
  310 FORMAT ('MatGetSizeSeq ERROR: NODE ',I6,' HAS INVALID DEGREE ',I4)
  320 FORMAT (/,/,'MatGetSizeSeq: STIFFNESS MATRIX',/,' ',15 ('='),/,
     +       10X,'MIN/MAX/AVG VERTEX DEGREE : ',2 (I2,2X),F4.1,/,10X,
     +       'NNZR = ',I10,' NONZERO BLOCK ENTRIES',/)
C
C
C
      RETURN

      END
