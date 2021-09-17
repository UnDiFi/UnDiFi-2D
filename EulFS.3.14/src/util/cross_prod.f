c
      SUBROUTINE cross_prod(a,b,c)
C
      IMPLICIT NONE
C
C     .. Array Arguments ..
C
      REAL*8 a(*) , b(*) , c(*)
C
C     .. Executable Statements ..
C
      c(1) = a(2) * b(3) - a(3) * b(2)
      c(2) = a(3) * b(1) - a(1) * b(3)
      c(3) = a(1) * b(2) - a(2) * b(1)
c     do m = 1 , 3
c     c(m) = a( JCYCL(m+1) ) * b( JCYCL(m+2) ) -
c    .       a( JCYCL(m+2) ) * b( JCYCL(m+1) )
c     enddo
c
      return
      end
C
