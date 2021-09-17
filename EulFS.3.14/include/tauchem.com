!
!     This common block contains chemical characteristic times
!
      DOUBLE PRECISION tau(NSP)
      COMMON/TAUCHEM/   tau
!
!
