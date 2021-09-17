C
C     $Id: bnd.h,v 1.6 2013/01/25 08:18:05 abonfi Exp $
C
      INTEGER   NBTYPE    ,NCOLOR   ,MBODIES
      PARAMETER(NBTYPE= 12,NCOLOR=50,MBODIES=NCOLOR)
C
C     NBTYPE is the number of boundary types currently implemented
C     NCOLOR is the number of colours
C     MBODIES is the max. number of bodies currently allowed
C
C     IF a NEW boundary type is add "src/blockdata.f" MUST
C     be accordingly modified
C
C
C     flags defining boundary conditions on boundary faces(edges) 
C
      INTEGER BC_TYPE_SUPS_INLET  ,BC_TYPE_SUBS_OUTLET ,
     &        BC_TYPE_SUPS_OUTLET ,BC_TYPE_SLIP_FREE   ,
     &        BC_TYPE_FAR_FIELD   ,BC_TYPE_NO_SLIP     ,
     &        BC_TYPE_PROFILE     ,BC_TYPE_SUBS_INLET  ,
     &        BC_TYPE_PERIODIC    ,BC_TYPE_X_SYMMETRY  ,
     &        BC_TYPE_Y_SYMMETRY  ,BC_TYPE_Z_SYMMETRY  ,
     &        BC_TYPE_PRESCRIBED_FLUX
      PARAMETER(BC_TYPE_SUPS_INLET  = 1,BC_TYPE_SUBS_OUTLET = 2,
     &          BC_TYPE_SUPS_OUTLET = 3,BC_TYPE_SLIP_FREE   = 4,
     &          BC_TYPE_FAR_FIELD   = 5,BC_TYPE_NO_SLIP     = 6,
     &          BC_TYPE_PROFILE     = 7,BC_TYPE_SUBS_INLET  = 8,
     &          BC_TYPE_PERIODIC    = 0,BC_TYPE_X_SYMMETRY  = 9,
     &          BC_TYPE_Y_SYMMETRY  =10,BC_TYPE_Z_SYMMETRY  =11,
     &          BC_TYPE_PRESCRIBED_FLUX  =12)
C
      INTEGER BC_TYPE_MIRROR,BC_TYPE_FLUX
      PARAMETER(BC_TYPE_MIRROR=102,BC_TYPE_FLUX=103)
C
C     BC_TYPE_MIRROR and BC_TYPE_WEAK are two different treatments for
C     the inviscid wall bndry condition BC_TYPE_SLIP_FREE
