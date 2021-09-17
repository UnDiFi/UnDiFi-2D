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
