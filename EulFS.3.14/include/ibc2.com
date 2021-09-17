C
C     $Id: ibc2.com,v 1.1 2013/01/25 08:09:32 abonfi Exp $
C
      INTEGER NCL,LCLDA
      LOGICAL RAD_EQUI,LCLTHERE
      PARAMETER (LCLDA=4)
      COMMON/COMIBC2/NCL
      COMMON/COMLBC2/RAD_EQUI,LCLTHERE
C
C     common variables related to outflow boundary conditions
C     of type 2: subsonic outflow
C     c-lines are stored in:
C     file006.dat
C     clin$nnn$.dat
C     in the sequential/parallel case
C
C     NCL number of c-lines, it is read from file006.dat
C     RAD_EQUI is set .TRUE. if the option
C     -radial_equilibrium is used
C
C     c-lines are stored in a linked list
C
C     A(LCLDA,*),JA(LCLDA,*),IA(1:NCL+1)
