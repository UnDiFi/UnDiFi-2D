C
C
      SUBROUTINE CSRFMT(NR,NELEM,NOFVERT,ICELNOD,IA,JA,IWK,JWK)
C
      IMPLICIT NONE
C
C     .. This routine creates the Compressed Sparse Row structure
C        of the stiffness matrix, i.e. the arrays:
C        ia(1:NR+1) and ja(1:NNZ)
C        It has been adapted from one of the routines of the
C        SparsKit package by Y. Saad ..
C
C
C
C
C ICELNOD -- Integer ICELNOD(1:NOFVERT,1:NELEM)
C            Cell to Node pointer : ICELNOD(i,ielem) gives the
C            global node number of the i-th vertex of the ielem-th cell
C
C
C JWK(1:NNZR) -- Integer workarray
C
C
C
C     .. Scalar Arguments ..
      INTEGER NELEM,NOFVERT,NR
C     ..
C     .. Array Arguments ..
      INTEGER IA(NR+1),ICELNOD(NOFVERT,NELEM),IWK(NR),JA(*),JWK(*)
C     ..
C     .. Local Scalars ..
      INTEGER IELEM,II,ILAST,IROWST,J,JJ,K,KA,KB,KSAV,KSAVN
C     ..
      KSAV = IA(1)
      IA(1) = 1
      DO 101 J = 2,NR + 1
          KSAVN = IA(J)
          IA(J) = IA(J-1) + KSAV
          IWK(J-1) = IA(J-1) - 1
  101 KSAV = KSAVN
      WRITE(6,*)'Elements in the matrix (NNZR) ',ia(NR+1)-ia(1)
c
c-----------------
c main loop
c-----------------
c
      DO 102 IELEM = 1,NELEM
c
c
          DO 120 KA = 1,NOFVERT
              II = ICELNOD(KA,IELEM)
c
c unpack row into jwk
c
              IROWST = IA(II)
              ILAST = IWK(II)
              DO 109 K = IROWST,ILAST
                  JWK(JA(K)) = K
  109         CONTINUE
c
              DO 108 KB = 1,NOFVERT
c
c column number = jj
c
                  JJ = ICELNOD(KB,IELEM)
c
                  K = JWK(JJ)
                  IF (K.EQ.0) THEN
                      ILAST = ILAST + 1
                      JWK(JJ) = ILAST
                      JA(ILAST) = JJ
                  ENDIF

  108         CONTINUE
c
c refresh jwk
c
              DO 119 K = IROWST,ILAST
                  JWK(JA(K)) = 0
  119         CONTINUE
              IWK(II) = ILAST
  120     CONTINUE
c
  102 CONTINUE
c
      WRITE(6,*)'Elements in the matrix (NNZR) ',ia(NR+1)-ia(1)
c
      RETURN

      END
