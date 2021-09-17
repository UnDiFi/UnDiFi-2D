      SUBROUTINE REORDR(N,A,INCX,ORD,WORK)

C
C     This routine is used to reorder the array A
C     according to the permutation array ORD
C     i.e. A(I) := A(ORD(I)) 
C     It should follow a previous call to SORTRX or QSORT
C
C     On entry
C     WORK must contain a copy of A
C     WORK(i) = A(j) i=1,N j = 1,INCX*N,INCX
C     A is allowed to have non unit stride
C
C

C     .. Scalar Arguments ..
      INTEGER INCX,N
C     ..
C     .. Array Arguments ..
      INTEGER A(*),ORD(N),WORK(N)
C     ..
C     .. Local Scalars ..
      INTEGER I,XADDR
C     ..
      XADDR = 1
C
      IF  ( INCX .LT. 0 )  THEN
          XADDR = (-N+1)*INCX + 1
      ENDIF
C
      DO 1 I = 1,N

          A(XADDR) = WORK(ORD(I))
          XADDR    = XADDR + INCX
    1 CONTINUE

      RETURN

      END
