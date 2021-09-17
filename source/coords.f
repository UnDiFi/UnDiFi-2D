      SUBROUTINE COORDS(X,Y,X1,X2,X3,Y1,Y2,Y3,R,IER)
      INTEGER IER
      DOUBLE PRECISION X,Y,X1,X2,X3,Y1,Y2,Y3,R(3)
C
C***********************************************************
C
C                                               ROBERT RENKA
C                                       OAK RIDGE NATL. LAB.
C                                             (615) 576-5139
C
C   THIS ROUTINE COMPUTES THE THREE BARYCENTRIC COORDINATES
C OF A POINT IN THE PLANE FOR A GIVEN TRIANGLE.
C
C INPUT PARAMETERS - X,Y - X AND Y COORDINATES OF THE POINT
C                          WHOSE BARYCENTRIC COORDINATES ARE
C                          DESIRED.
C
C      X1,X2,X3,Y1,Y2,Y3 - COORDINATES OF THE VERTICES OF
C                          THE TRIANGLE.
C
C INPUT PARAMETERS ARE NOT ALTERED BY THIS ROUTINE.
C
C OUTPUT PARAMETERS -  R - 3-VECTOR OF BARYCENTRIC COORDI-
C                          NATES UNLESS IER = 1.  NOTE THAT
C                          R(I) .LT. 0. IFF (X,Y) IS TO THE
C                          RIGHT OF THE VECTOR FROM VERTEX
C                          I+1 TO VERTEX I+2 (CYCLICAL
C                          ARITHMETIC).
C
C                    IER - ERROR INDICATOR
C                          IER = 0 IF NO ERRORS WERE
C                                  ENCOUNTERED.
C                          IER = 1 IF THE VERTICES OF THE
C                                  TRIANGLE ARE COLLINEAR.
C
C MODULES REFERENCED BY COORDS - NONE
C
C***********************************************************
C
      DOUBLE PRECISION U(3),V(3),AREA,XP,YP
C
C LOCAL PARAMETERS -
C
C U(K),V(K) = X AND Y COMPONENTS OF THE VECTOR REPRESENTING
C               THE SIDE OPPOSITE VERTEX K FOR K = 1,2,3.
C AREA =      TWICE THE AREA OF THE TRIANGLE.
C XP,YP =     X-X1, Y-Y1
C
      U(1) = X3 - X2
      U(2) = X1 - X3
      U(3) = X2 - X1
C
      V(1) = Y3 - Y2
      V(2) = Y1 - Y3
      V(3) = Y2 - Y1
C
C AREA = 3-1 X 3-2
C
      AREA = U(1)*V(2) - U(2)*V(1)
      IF (AREA.EQ.0.D0) GOTO 1
C
C R(1) = (2-3 X 2-(X,Y))/AREA, R(2) = (1-(X,Y) X 1-3)/AREA,
C   R(3) = (1-2 X 1-(X,Y))/AREA
C
      R(1) = (U(1)* (Y-Y2)-V(1)* (X-X2))/AREA
      XP = X - X1
      YP = Y - Y1
      R(2) = (U(2)*YP-V(2)*XP)/AREA
      R(3) = (U(3)*YP-V(3)*XP)/AREA
      IER = 0
      RETURN
C
C VERTICES ARE COLLINEAR
C
    1 IER = 1
      RETURN

      END
