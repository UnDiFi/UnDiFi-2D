      SUBROUTINE RPERM(P,N)
!===Generate a random permutation, P, of the first N integers.
!   (Shuffling, equivalent to sampling WITHOUT REPLACEMENT).
!   Adaptation of Knuth Volume 2, Algorithm 3.4.2P.
      INTEGER N,P(N), K,J,I,IPJ,ITEMP,M
      REAL U(100)
      DO I=1,N
        P(I)=I
      END DO
!---Generate up to 100 U(0,1) numbers at a time.
      DO I=1,N,100
        M=MIN(N-I+1,100)
        CALL RAND(U,M)
        DO J=1,M
          IPJ=I+J-1
          K=INT(U(J)*(N-IPJ+1))+IPJ
          ITEMP=P(IPJ)
          P(IPJ)=P(K)
          P(K)=ITEMP
        END DO
      END DO
      RETURN
      END
      SUBROUTINE IRAND(S,N,LOW,HI)
!===Generate a random integer sequence: S(1),S(2), ... ,S(N)
!   such that each element is in the closed interval <LOW,HI> and
!   sampled WITH REPLACEMENT.                            HDK, JUNE 1971.
      INTEGER N,S(N),LOW,HI,IX,I
      REAL U(1)
      DOUBLE PRECISION X
      DO I=1,N
        CALL RAND(U,1)
!---Use DP arithmetic to effect a more precise transformation.
        X=DBLE((HI+1)-LOW)*U(1) + DBLE(LOW)
        IX=X
        IF(X.LT.0 .AND. IX.NE.X) IX=X-1.D0
        S(I)=IX
      END DO
      RETURN
      END
      BLOCK DATA
!=======================================================================
!  Portable pseudo-random integer generator, especially for
!  microcomputers with a limitation of 16 bit integers. Translated from
!  original Pascal version(1) to Fortran 77 by H. D. Knoble, PSU.
!
!  Corrected typo (change 30308.0 to 30307.0) courtesy of Ed Wynn
!  and Dave Dodson, 3 May 2002.
!
!   The supporting paper is:
!   (1) B. Wichmann & D. Hill, "Building a Random-Number Generator",
!             BYTE, March, 1987, 12(3):127-128.
!
!   Also see the following related works:
!   (2) Wichmann, B.A. and Hill, I.D., "An Efficient and Portable",
!             Pseudo-random Number Generator", Applied Statistics,
!             Volume 31, 1982, pages 188-190.
!   (3) Haas, Alexander, "The Multiple Prime Random Number Generator",
!             ACM Transactions on Mathematical Software; December,
!             1987, 13(4):368-381.
!   (4) L'Ecuyer, Pierre, "Efficient and Portable Combined Random Number
!             Generators", Communications of the ACM; June, 1988,
!             31(6):742-749,774.
!
! Use...
!      CALL RAND(U,N)
!          To generate a sequence, U, of N Uniform(0,1) numbers.
!          Cycle length is ((30269-1)*(30307-1)*(30323-1))/4  or
!          6953607871644  > 6.95E+12.
!
!     To access the SEED vector in the calling program use statements:
!     INTEGER SEED(3)
!     COMMON/RANDOM/SEED
!
!  The common variable SEED is the array of three current seeds.
      INTEGER SEED(3)
      COMMON/RANDOM/SEED
      DATA SEED(1),SEED(2),SEED(3)/1,10000,3000/
      END
!=======================================================================
      SUBROUTINE RAND(U,N)
      INTEGER N,X,Y,Z
      REAL U(N),V
      COMMON/RANDOM/X,Y,Z
      IF(N.LE.0) RETURN
      DO  I=1,N
        X=171*MOD(X,177)-2*(X/177)
        IF(X.LT.0) X=X+30269
        Y=172*MOD(Y,176)-35*(Y/176)
        IF(Y.LT.0) Y=Y+30307
        Z=170*MOD(Z,178)-63*(Z/178)
        IF(Z.LT.0) Z=Z+30323
        V=X/30269.0 + Y/30307.0 + Z/30323.0
        U(I)=V-INT(V)
      END DO
      RETURN
      END
