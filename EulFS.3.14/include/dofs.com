      INTEGER IE,IX,IY,IZ
      COMMON/DOFCOM/IE,IX,IY,IZ
C
C     $Id: dofs.com,v 1.1 2013/04/30 07:02:02 abonfi Exp $
C
C     degree of freedom where the x,y,z momentum is stored
C
C     IE = 2 for   compressible flow (IABS(KAN)=4)
C     IE = -1 for incompressible flow (IABS(KAN)=4)
C     IX = 2 for incompressible flow (IABS(KAN)=2)
C     IX = NSP+2 for plasma     flow (IABS(KAN)=3)
C     IX = 3 for   compressible flow (IABS(KAN)=4)
