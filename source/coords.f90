SUBROUTINE COORDS(X,Y,X1,X2,X3,Y1,Y2,Y3,R,IER)
  INTEGER IER
  DOUBLE PRECISION X,Y,X1,X2,X3,Y1,Y2,Y3,R(3)
!
!***********************************************************
!
!                                               ROBERT RENKA
!                                       OAK RIDGE NATL. LAB.
!                                             (615) 576-5139
!
!   THIS ROUTINE COMPUTES THE THREE BARYCENTRIC COORDINATES
! OF A POINT IN THE PLANE FOR A GIVEN TRIANGLE.
!
! INPUT PARAMETERS - X,Y - X AND Y COORDINATES OF THE POINT
!                          WHOSE BARYCENTRIC COORDINATES ARE
!                          DESIRED.
!
!      X1,X2,X3,Y1,Y2,Y3 - COORDINATES OF THE VERTICES OF
!                          THE TRIANGLE.
!
! INPUT PARAMETERS ARE NOT ALTERED BY THIS ROUTINE.
!
! OUTPUT PARAMETERS -  R - 3-VECTOR OF BARYCENTRIC COORDI-
!                          NATES UNLESS IER = 1.  NOTE THAT
!                          R(I) .LT. 0. IFF (X,Y) IS TO THE
!                          RIGHT OF THE VECTOR FROM VERTEX
!                          I+1 TO VERTEX I+2 (CYCLICAL
!                          ARITHMETIC).
!
!                    IER - ERROR INDICATOR
!                          IER = 0 IF NO ERRORS WERE
!                                  ENCOUNTERED.
!                          IER = 1 IF THE VERTICES OF THE
!                                  TRIANGLE ARE COLLINEAR.
!
! MODULES REFERENCED BY COORDS - NONE
!
!***********************************************************
!
  DOUBLE PRECISION U(3),V(3),AREA,XP,YP
!
! LOCAL PARAMETERS -
!
! U(K),V(K) = X AND Y COMPONENTS OF THE VECTOR REPRESENTING
!               THE SIDE OPPOSITE VERTEX K FOR K = 1,2,3.
! AREA =      TWICE THE AREA OF THE TRIANGLE.
! XP,YP =     X-X1, Y-Y1
!
  U(1) = X3 - X2
  U(2) = X1 - X3
  U(3) = X2 - X1
!
  V(1) = Y3 - Y2
  V(2) = Y1 - Y3
  V(3) = Y2 - Y1
!
! AREA = 3-1 X 3-2
!
  AREA = U(1)*V(2) - U(2)*V(1)
  IF (AREA.EQ.0.D0) GOTO 1
!
! R(1) = (2-3 X 2-(X,Y))/AREA, R(2) = (1-(X,Y) X 1-3)/AREA,
!   R(3) = (1-2 X 1-(X,Y))/AREA
!
  R(1) = (U(1)* (Y-Y2)-V(1)* (X-X2))/AREA
  XP = X - X1
  YP = Y - Y1
  R(2) = (U(2)*YP-V(2)*XP)/AREA
  R(3) = (U(3)*YP-V(3)*XP)/AREA
  IER = 0
  RETURN
!
! VERTICES ARE COLLINEAR
!
1 IER = 1
  RETURN

end subroutine COORDS
