      subroutine MY_GEMV(TRANS, M, N, ALPHA, A,  LDA,  X,  INCX,
     +                   BETA, Y, INCY )

      DOUBLE       PRECISION ALPHA, BETA

      INTEGER      INCX, INCY, LDA, M, N

      CHARACTER*1  TRANS

      DOUBLE       PRECISION A( LDA, * ), X( * ), Y( * )


      IF( M .NE. N )STOP 'M MUST = N in MY_GEMV'

      GOTO(2,3,4,5)M-1

      STOP 'M out of range in MY_GEMV'

    2 CONTINUE

      Y(1) = ALPHA*(A(1,1)*X(1)+A(1,2)*X(2)) + BETA*Y(1) 
      Y(2) = ALPHA*(A(2,1)*X(1)+A(2,2)*X(2)) + BETA*Y(2) 
      RETURN 

    3 CONTINUE

      Y(1) = ALPHA*(A(1,1)*X(1)+A(1,2)*X(2)+A(1,3)*X(3)) + BETA*Y(1) 
      Y(2) = ALPHA*(A(2,1)*X(1)+A(2,2)*X(2)+A(2,3)*X(3)) + BETA*Y(2) 
      Y(3) = ALPHA*(A(3,1)*X(1)+A(3,2)*X(2)+A(3,3)*X(3)) + BETA*Y(3) 
      RETURN 

    4 CONTINUE

      Y(1) = ALPHA*(A(1,1)*X(1)+A(1,2)*X(2)+A(1,3)*X(3)+A(1,4)*X(4))
     &+ BETA*Y(1) 
      Y(2) = ALPHA*(A(2,1)*X(1)+A(2,2)*X(2)+A(2,3)*X(3)+A(2,4)*X(4)) 
     &+ BETA*Y(2) 
      Y(3) = ALPHA*(A(3,1)*X(1)+A(3,2)*X(2)+A(3,3)*X(3)+A(3,4)*X(4)) 
     &+ BETA*Y(3) 
      Y(4) = ALPHA*(A(4,1)*X(1)+A(4,2)*X(2)+A(4,3)*X(3)+A(4,4)*X(4)) 
     &+ BETA*Y(4) 
      RETURN 
    5 CONTINUE

      Y(1) = ALPHA*(A(1,1)*X(1)+A(1,2)*X(2)+A(1,3)*X(3)+A(1,4)*X(4)+
     &              A(1,5)*X(5)) + BETA*Y(1) 
      Y(2) = ALPHA*(A(2,1)*X(1)+A(2,2)*X(2)+A(2,3)*X(3)+A(2,4)*X(4)+ 
     &              A(2,5)*X(5)) + BETA*Y(2) 
      Y(3) = ALPHA*(A(3,1)*X(1)+A(3,2)*X(2)+A(3,3)*X(3)+A(3,4)*X(4)+ 
     &              A(3,5)*X(5)) + BETA*Y(3) 
      Y(4) = ALPHA*(A(4,1)*X(1)+A(4,2)*X(2)+A(4,3)*X(3)+A(4,4)*X(4)+ 
     &              A(4,5)*X(5)) + BETA*Y(4) 
      Y(5) = ALPHA*(A(5,1)*X(1)+A(5,2)*X(2)+A(5,3)*X(3)+A(5,4)*X(4)+ 
     &              A(5,5)*X(5)) + BETA*Y(5) 
      RETURN
      END 
