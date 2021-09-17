      Subroutine Add_to_List( IV , NumItems , NewItem ) 
C
C       .. This routines adds to the list IV the element NewItem
C       and keeps the ascending order of the elements ..
C       .. If NewItem is already in the list, the control is
C       returned to the calling routine
C
      IMPLICIT NONE
C
C     .. Scalar Arguments ..
C
      INTEGER       NumItems,NewItem
C
C     .. Array Arguments ..
C
      INTEGER       IV(1)
C
C     .. Local Scalars ..
C
      INTEGER       Lower,Middle,Upper,N,I
C
C     .. External Subroutines ..
C
C
C     .. External Functions ..
C
      INTEGER       binsrch
      EXTERNAL       binsrch
C
C     .. Intrinsic Functions ..
C
C
C     .. Data Statements ..
C
      DATA Lower,Middle,Upper / 0,0,0 /
C
C     .. Executable Statements ..
C
C       .. If NewItem is > IV(NumItems) it is simply added
C       at the end of the list
C
      IF( NumItems .EQ. 0 .OR. (NewItem .GT. IV(NumItems)) )THEN
       NumItems = NumItems + 1
       IV(NumItems) = NewItem
      ELSE
       N = binsrch( IV , NumItems , NewItem , Lower , Middle , Upper )
C
       IF( N .NE. -1 )RETURN ! If already in the list: do NOTHING
C
C
C
       NumItems = NumItems + 1
       DO 10 i = NumItems , Upper+1 , -1
         IV(i) = IV(i-1)
   10       CONTINUE
              IV(Upper) = NewItem
      ENDIF
C
      RETURN
      END
C
C
C
      INTEGER FUNCTION binsrch( A , N , X , lower , middle , upper )
C
      INTEGER A(*)
      INTEGER lower,middle,upper,X
C
      lower = 1
      upper = N
C
    1 IF    ( lower .GT. upper )THEN
        binsrch = -1
        RETURN
      ELSEIF( lower .LT. upper )THEN
        middle = (lower+upper)/2
        IF( X .GT. A(middle) )THEN
          lower = middle + 1
        ELSE
          upper = middle
        ENDIF
        GOTO 1
      ELSE
        IF( X .EQ. A(lower) )THEN
          binsrch = lower
        ELSE
          binsrch = -1
        ENDIF
        RETURN
      ENDIF
      END
