      SUBROUTINE GAUSS(NPO)
C
C****************************************************************/
C* GAUSS.F                                                      */
C*--------------------------------------------------------------*/
C* FORMS ARRAYS OF COORDINATES AND WEIGHTS AT INTEGRATION       */
C* POINTS FOR TRIANGULAR ELEMENTS (HAMMER'S RULES)              */
C* (Courtesy of Jean Christophe)
C* Input                                                        */
C*      Npts :          number of integration points            */
C* "Output"                                                     */
C*      Gcoor :    triangular coordinates of integration        */
C*                 points                                       */
C*      Gwght :    weights at the integration points            */
C*
C*  Watch out : the 7 points weights and/or
C*              coordinates are probrably wrong.
C****************************************************************/
C
C
C
c
c
C     .. Parameters ..
      INTEGER IPMAX
      PARAMETER (IPMAX=7)
      DOUBLE PRECISION ONE,TWO,THREE,FOUR,FIVE,SIX,ZERO
      PARAMETER (ONE=1.D0,TWO=2.D0,THREE=3.D0,FOUR=4.D0,FIVE=5.D0,
     +          SIX=6.D0,ZERO=0.D0)
C     ..
C     .. Scalar Arguments ..
      INTEGER NPO
C     ..
C     .. Scalars in Common ..
      INTEGER NPTS
C     ..
C     .. Arrays in Common ..
      DOUBLE PRECISION GCOOR(IPMAX,3),GWGHT(IPMAX)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION C1,C2,C3,C4,C5,C6,C7,C8
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC DSQRT
C     ..
C     .. Common blocks ..
      COMMON /GAUSSBLOCK/GCOOR,GWGHT,NPTS
C     ..
      NPTS=NPO
      GOTO (10,999,30,40,999,999,70) NPTS

  999 WRITE (6,FMT=*) 'Npts out of range'
      STOP

   10 CONTINUE
      C1 = ONE/THREE
      GWGHT(1) = ONE/TWO
      GCOOR(1,1) = C1
      GCOOR(1,2) = C1
      GCOOR(1,3) = C1
      RETURN

   30 CONTINUE
      C1 = ONE/SIX
      C2 = ONE/TWO
      GWGHT(1) = C1
      GWGHT(2) = C1
      GWGHT(3) = C1
      GCOOR(1,1) = ZERO
      GCOOR(1,2) = C2
      GCOOR(1,3) = C2
      GCOOR(2,1) = C2
      GCOOR(2,2) = ZERO
      GCOOR(2,3) = C2
      GCOOR(3,1) = C2
      GCOOR(3,2) = C2
      GCOOR(3,3) = ZERO
      RETURN

   40 CONTINUE
      C1 = ONE/THREE
      C2 = -27.D0/96.D0
      C3 = 25.D0/96.D0
      C4 = THREE/FIVE
      C5 = ONE/FIVE
      GWGHT(1) = C2
      GWGHT(2) = C3
      GWGHT(3) = C3
      GWGHT(4) = C3
      GCOOR(1,1) = C1
      GCOOR(1,2) = C1
      GCOOR(1,3) = C1
      GCOOR(2,1) = C4
      GCOOR(2,2) = C5
      GCOOR(2,3) = C5
      GCOOR(3,1) = C5
      GCOOR(3,2) = C4
      GCOOR(3,3) = C5
      GCOOR(4,1) = C5
      GCOOR(4,2) = C5
      GCOOR(4,3) = C4
      RETURN

   70 CONTINUE
      C1 = 9.D0/80.D0
      C2 = (155.D0+DSQRT(15.D0))/2400.D0
      C3 = 31.D0/240.D0 - C2
      C4 = (SIX+DSQRT(15.D0))/21.D0
      C5 = FOUR/7.D0 - C4
      C6 = ONE/THREE
      C7 = ONE - TWO*C4
      C8 = ONE - TWO*C5
      GWGHT(1) = C1
      GWGHT(2) = C2
      GWGHT(3) = C2
      GWGHT(4) = C2
      GWGHT(5) = C3
      GWGHT(6) = C3
      GWGHT(7) = C3
      GCOOR(1,1) = C1
      GCOOR(1,2) = C1
      GCOOR(1,3) = C1
      GCOOR(2,1) = C7
      GCOOR(2,2) = C4
      GCOOR(2,3) = C4
      GCOOR(3,1) = C4
      GCOOR(3,2) = C7
      GCOOR(3,3) = C4
      GCOOR(4,1) = C4
      GCOOR(4,2) = C4
      GCOOR(4,3) = C7
      GCOOR(5,1) = C8
      GCOOR(5,2) = C5
      GCOOR(5,3) = C5
      GCOOR(6,1) = C5
      GCOOR(6,2) = C8
      GCOOR(6,3) = C5
      GCOOR(7,1) = C5
      GCOOR(7,2) = C5
      GCOOR(7,3) = C8
      RETURN

      END
