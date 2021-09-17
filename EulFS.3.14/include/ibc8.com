C
C     $Id: ibc8.com,v 1.2 2013/01/25 08:04:29 abonfi Exp $
C
      INTEGER nVarsInlet
      PARAMETER(nVarsInlet=6)
      INTEGER LKLIST,LVLIST,NLIST
caldo LOGICAL LREAD
      COMMON/COMIBC8/LKLIST,LVLIST,NLIST
caldo COMMON/COMLBC8/LREAD
C
C     common variables related to inflow boundary conditions
C     of type 8: subsonic inflow
C     inflow vertices with corresponding boundary conditions
C     are stored in:
C     file005.dat
C     pbcs$nnn$.dat
C     in the sequential/parallel case
C
C     LREAD(1)  == .TRUE. whenever there is a file to read from
C     LREAD(1) is now in flags.com
C     LKLIST POINTER to the list of nodes with ibc=8
C     LVLIST POINTER to the array of inlet bcs
C     NLIST  number of meshpoints with ibc=8
C
C
C     up to 0.11.8 the three entries in VLIST(1:*,*) are:
C
C     p/p0 t/t0 u/u_ref
C
C     starting with 0.11.9 the six entries in VLIST(1:*,*) are:
C
C     p/p0 t/t0 unused n_x n_y n_z
C
