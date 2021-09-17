      BLOCK DATA
C
      IMPLICIT NONE
C
C BLOCK DATA for EulFS
C
C
C     .. Parameters ..
C
      INCLUDE 'paramt.h'
      INCLUDE 'bnd.h'
      INCLUDE 'constants.h'
      INCLUDE 'implicit.h'
      INTEGER V3,VN,N2,N3,LENB
      PARAMETER(V3=MAXNOFVAR*3,VN=MAXNOFVAR*NMAX,N2=2*MAXNOFVAR,
     &N3=MAXTIMLEVS*MAXNOFVAR,LENB=(MAXNOFVERT-1)*(MAXNOFEQN**2))
C
C     .. Commons ..
C
      INCLUDE 'bnd.com'
      INCLUDE 'io.com'
      INCLUDE 'conv.com'
      INCLUDE 'transf.com'
      INCLUDE 'nboun.com'
      INCLUDE 'verbose.com'
      INCLUDE 'three.com'
      INTEGER LXX(30)
      COMMON/NLOC/LXX
C
C     .. Local Scalars ..
C
      INTEGER I,J
C
C     .. Executable Statements ..
C
C
C============================================================
C       INITIALIZE LABELLED COMMONS
C============================================================
C
C---------- COMMON /io.com/
C
      DATA IHST1,IHST2,IHST3,IHST4,ITIM1,IPROBE/1,2,3,116,4,55/
C
C---------- COMMON /VERBOSE/
C
      DATA IVERBOSE/0/
C
C---------- COMMON /CONV/
C
      DATA RESMAX,RESL2,DELMAX,DELL2,RESL20,RESMAX0,CFL,CFLMAX,TOLER,
     &OMEGA/N2*ZERO,N2*ZERO,N2*ZERO,N2*ZERO,8*ONE,ZERO,ONE/
      DATA INMAX,INDEL,NITER,ITMAX,IVCNVG,ISTMP,ISTART,IBAK
     &/N2*1,N2*1,0,5*1/
C
C     The right eigenvectors matrix is initialized to the identity matrix
C
      DATA (dUdV(J) , J=1,5   ) / ONE,ZERO,ZERO,ZERO,ZERO/
      DATA (dUdV(J) , J=6,10  ) / ZERO,ONE,ZERO,ZERO,ZERO/
      DATA (dUdV(J) , J=11,15 ) / ZERO,ZERO,ONE,ZERO,ZERO/
      DATA (dUdV(J) , J=16,20 ) / ZERO,ZERO,ZERO,ONE,ZERO/
      DATA (dUdV(J) , J=21,25 ) / ZERO,ZERO,ZERO,ZERO,ONE/
C
      DATA (dVdZ(J) , J=1,5   ) / ONE,ZERO,ZERO,ZERO,ZERO/
      DATA (dVdZ(J) , J=6,10  ) / ZERO,ONE,ZERO,ZERO,ZERO/
      DATA (dVdZ(J) , J=11,15 ) / ZERO,ZERO,ONE,ZERO,ZERO/
      DATA (dVdZ(J) , J=16,20 ) / ZERO,ZERO,ZERO,ONE,ZERO/
      DATA (dVdZ(J) , J=21,25 ) / ZERO,ZERO,ZERO,ZERO,ONE/
C
      DATA DZDU/LENB*ZERO/ 
C
C---------- COMMON /THREE/
C
      DATA ZAVG/N3*ZERO/
      DATA UAVG/NMAX*ZERO/
C
C---------- COMMON /FIX/
C
caldo DATA STAGFIX / ONE /
C
C---------- COMMON /NLOC/
C
      DATA LXX/30*1/
C
C---------- COMMON /bnd/
C
      DATA (IMUNIT(J),J=0,NCOLOR) / 11,12,13,14,15,16,17,18,19,20,21,22,
     1                              23,24,25,26,27,28,29,30,31,32,33,34,
     2                              35,36,37,38,39,40,41,42,43,44,45,46,
     3                              47,48,49,50,51,52,53,54,100,101,102,
     4                              103,104,105,106/
      DATA (IFUNIT(J),J=0,NCOLOR) / 99,98,97,96,95,94,93,92,91,90,89,88,
     1                              87,86,85,84,83,82,81,80,79,78,77,76,
     2                              75,74,73,72,71,70,69,68,67,66,65,64,
     4                              63,62,61,60,59,58,57,56,107,108,109,
     4                              110,111,112,113/
      DATA (MCOLOR(J),J=0,NCOLOR) / 51*0 /
      DATA (SCOLOR(J),J=0,NCOLOR) / 51*ZERO /
      DATA (CBTYPE(J),J=0,NBTYPE) /'PERIODIC            ',
     1'SUPERSONIC INLET    ','SUBSONIC OUTLET     ',
     2'SUPERSONIC OUTLET   ','INVISCID WALL       ',
     3'FAR FIELD           ','VISCOUS WALL        ',
     4'PRESCRIBED VELOCITY ','SUBSONIC INLET      ',
     5'X SYMMETRY          ','Y SYMMETRY          ',
     6'Z SYMMETRY          ','PRESCIBED FLUX      '/
      DATA IOALE/115/
C
C
C----6----- COMMON /impsol/
C
      DATA TIMEIMPL,PICARD,NEWTON/2*.TRUE.,.FALSE./
C
C---------- COMMON /nboun/
C
      DATA NPOIN1,NPOIN6,NPOIN7/3*0/
C
      END
