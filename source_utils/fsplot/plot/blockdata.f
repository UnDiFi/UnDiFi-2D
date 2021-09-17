      BLOCK DATA
C
      IMPLICIT NONE
C
C BLOCK DATA for EulFS
C
      INCLUDE 'constants'
      INCLUDE 'bnd.h'
      INCLUDE 'bnd'
C     INCLUDE 'alloc'
      INTEGER LXX(10)
      COMMON/NLOC/LXX
C
      INTEGER*4	I,J
C
C============================================================
C       INITIALIZE LABELLED COMMONS
C============================================================
C
C
C
C---------- COMMON /ALLOC/
C
C     DATA IVA/1/,IVAMAX/1/,NTBL/5/
C
C---------- COMMON /NLOC/
C
      DATA LXX/10*1/
C
C.......... DEFINE HERE THE NUMBER OF INTEGERS CONTAINED IN A REAL
C
C               SINGLE PRECISION        NREEL.EQ.1
C               DOUBLE PRECISION        NREEL.EQ.2
C
C     DATA NREEL/2/
C
C---------- COMMON /bnd/
C
      DATA (MCOLOR(J),J=0,NCOLOR) / 51*0 /
      DATA (CBTYPE(J),J=0,NBTYPE) /
     &20HSUBSONIC INLET       ,
     120HSUPERSONIC INLET     ,
     220HSUBSONIC OUTLET      ,
     320HSUPERSONIC OUTLET    ,
     420HINVISCID WALL        ,
     520HFAR FIELD            /
C
      END
