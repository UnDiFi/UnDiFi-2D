      INTEGER MAXNVT
      PARAMETER(MAXNVT=8000)
      DOUBLE PRECISION tmpconv(MAXNVT)
      DOUBLE PRECISION tmpdif1(MAXNVT)
      DOUBLE PRECISION tmpdif2(MAXNVT)
      DOUBLE PRECISION tmpsou1(MAXNVT)
      DOUBLE PRECISION tmpsou2(MAXNVT)
      DOUBLE PRECISION tmpsou3(MAXNVT)
      DOUBLE PRECISION tmpdiff(MAXNVT)
      DOUBLE PRECISION tmpsum(MAXNVT)
      DOUBLE PRECISION tmpdum(MAXNVT)
      INTEGER KCN(4)
      COMMON/TDEBUG/tmpconv,tmpdif1,tmpdif2,tmpsou1,tmpsou2,
     +tmpsou3,tmpdiff,tmpsum,tmpdum,kcn
