      SUBROUTINE BNDRY_ISET(NodalBcs, IDX, IBNDFAC, NBFAC, ICELNOD,
     &NOFVERT, NODCODE, NDIM, NOFVAR, NPOIN, NGHOST, NPNOD, NELEM )
C
C     $Id: bndry_iset.F,v 1.4 2020/04/23 09:46:19 abonfi Exp $
C
C     This routine creates PETSc Index Sets addressing
C          boundary meshpoints with different colours
C
C     Warnings:
C     --------
C     Nodes are adressed by local numbering 1-based indexing
C
C     NODCODE will address interior nodes as well ghost ones
C     so it must be allocated for NPOIN+NGHOST+NPNOD in the calling routine
C
C     Note that Index sets are NOT created, if empty
C
C     It might be BUGGY in the parallel, periodic case
C
C
C
#include "petsc/finclude/petscis.h"
      use petscis
      IMPLICIT NONE
      INCLUDE 'bnd.h'
      INCLUDE 'bctype.h'
C
C
C
      INCLUDE 'bnd.com'
      INCLUDE 'io.com'
      INCLUDE 'flags.com'
      INCLUDE 'visco.com'
      INCLUDE 'ibc8.com'
C
      IS, dimension(0:*) :: NodalBcs
C
      INTEGER MY_PE
      COMMON/MPICOM/MY_PE
C
C     Scalar arguments:
C
      INTEGER NPOIN,NDIM,NOFVAR,NOFVERT,NGHOST,NPNOD,NBFAC,NELEM
C
C     Array arguments:
C
      INTEGER NODCODE(NPOIN+NGHOST),IDX(NPOIN),IBNDFAC(3,NBFAC)
      INTEGER ICELNOD(NOFVERT,*)
C
C     Local scalars:
C
      INTEGER IPOIN,NI,IFAIL,IBC,IVERT,I,ICLR,IFACE,IELEM
      INTEGER ICYCL
C
C
      WRITE(NOUT,1000)MY_PE
C
      DO 50 IBC = 0, NCOLOR ! loop over all colours
C
         IF( MCOLOR(IBC) .EQ. 0 )GOTO 50 ! if empty, skip
C
C     +-----------------------------------------------------+
C     Test for 
C     +-----------------------------------------------------+
C
C
         CALL IINIT(NPOIN+NGHOST+NPNOD,0,NODCODE,1)
C
         NI = 0
         DO 1 IFACE = 1 , NBFAC ! loop over all boundary faces
C
            ICLR = (IBNDFAC(3,IFACE))
C
C        faces with no-slip b.c. or prescribed velocity profile
C
            IF( .NOT.(ICLR .EQ. IBC) )GOTO 1 ! skip current face
C
               NI = NI + 1
               IELEM = IBNDFAC(1,IFACE)
               IVERT = IBNDFAC(2,IFACE)
C
               DO 2 I =  1, NOFVERT - 1
C
                   IPOIN = ICELNOD( ICYCL( IVERT+I , NOFVERT ) , IELEM )
C
                   NODCODE(IPOIN) = IBC
C
    2         CONTINUE ! end loop over the vertices of the bndry face
    1    CONTINUE
C
         IF( NI .NE. MCOLOR(IBC) )THEN
             WRITE(NOUT,*)'Smthg. went wrong in bndry_iset ',NI,
     &MCOLOR(IBC)
             CALL EXIT(1)
         ENDIF 
C
         NI=0
         DO 230 IPOIN = 1, NPOIN+NGHOST+NPNOD ! this should be checked in the parallel and/or periodic case
            IF(NODCODE(IPOIN).EQ.IBC)THEN
                NI=NI+1
                IDX(NI) = IPOIN
            ENDIF
  230 CONTINUE
         CALL ISCreateGeneral(PETSC_COMM_SELF,NI,IDX,PETSC_COPY_VALUES,
     &                        NodalBcs(IBC),IFAIL)
         WRITE(NOUT,2100)IBC,NI
   50 CONTINUE ! loop over colours (i.e. patches)
C
C
C
      RETURN
 1000 FORMAT(/,5X,'NODAL BOUNDARY COLOURS ON PE # ',I4,/ 5X,38("="))
 
 2100 FORMAT(5X,' BOUNDARY PATCH COLOURED ',I2,' HAS ',I6,' NODES ')
      END
C
