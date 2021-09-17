      INTEGER NTRIP,LOCTIA,LOCTJA,LOCTRP,LOCDXT,LOCTST,LTTD,LTTI
C
C     NTRIP number of trip points 
C     LTNTRP pointer to the linear list of trip point indexes
C     LOCTIA,LOCTJA pointers to the linked-list storing 
C                   the elements sharing the trip points
C
      COMMON /TRIPCOM/ NTRIP,LOCTIA,LOCTJA,LOCTRP,LOCDXT,LOCTST,
     +LTTD,LTTI
C
C     CHARACTER*80 TRIPFILE
C
