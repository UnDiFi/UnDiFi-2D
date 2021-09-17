C
      INTEGER FUNCTION GR_LINE_CROSS (X1,Y1,X2,Y2,X3,Y3,X4,Y4,S,T)
C
      IMPLICIT NONE
      DOUBLE PRECISION X1,Y1,X2,Y2,X3,Y3,X4,Y4
C 
C..........
C	      Function to determine if the line segments
C	      (X1,Y1)-(X2,Y2) and (X3,Y3)-(X4,Y4) cross
C
C	      This is done by computing the parameterized intersection
C	      and seeing if the parameters for both lines are in the
C	      interval [0,1]
C	      For this problem use Cramer's Rule to solve the 
C	      equations
C	      (x2-x1) T + (x3-x4) S = (x3 - x1)
C	      (y2-y1) T + (y3-y4) S = (y3 - y1)
C
C	      If there is no solution, the lines are parallel so they
C	      do not cross anyway.
C	
C	      S and T are the parameters of the crossing on each line segment
C	      S goes from 0 to 1 as we go from 3 - 4
C	      T goes from 0 to 1 as we go from 1 - 2
C
C	      GR_LINE_CROSS = 1 if the lines cross, 0 otherwise
C..........
C
      DOUBLE PRECISION S,T,XX1,XX2,YY1,YY2,XXX,YYY,DET,EPSILON
      PARAMETER (EPSILON = 1.E-20)
C	
      GR_LINE_CROSS = 0
      S = 0.
      T = 0.
      XX1 = X2 - X1
      XX2 = X3 - X4
      XXX = X3 - X1
      YY1 = Y2 - Y1
      YY2 = Y3 - Y4
      YYY = Y3 - Y1
      DET = XX1*YY2 - XX2*YY1
      IF (ABS(DET) .LT. EPSILON) THEN
C         WRITE(*,*) 'DET too small'
          GO TO 1
      ENDIF
      T = (XXX*YY2 - XX2*YYY) / DET
      S = (XX1*YYY - XXX*YY1) / DET
C
C..........
C            Check to see if there is an intersection within the parameter
C            ranges [0,1) 
C..........
C
      IF (T .GE. 0 .AND. T .LE. 1) THEN
          IF (S .GT. 0.0 .AND. S .LE. 1.0 .AND. T.NE.1) THEN
              GR_LINE_CROSS = 1
          ELSE IF (S .EQ. 1.D-10)THEN
                GR_LINE_CROSS = 10000
          ENDIF
      ENDIF
C
 1    RETURN
      END
