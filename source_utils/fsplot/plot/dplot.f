C
      SUBROUTINE DPLOT(plotfile,POINT,ICELNOD,ZROE,BNDPNTR)
C
      IMPLICIT NONE
C
CWRITES a DPLOT (2D) format datafile
C
      include 'constants'
      include 'bnd.h'
      include 'bnd'
      include 'dim_flags'
      include 'int_flags'
      include 'mesh_i4'
      include 'IO'
      double precision fren
      common /freestream/ fren
C
C	.. Arguments ..
C
      DOUBLE PRECISION ZROE(NOFVAR,1),POINT(DIM,1)
      INTEGER ICELNOD(NOFVERT,1),BNDPNTR(3,1)
C
      INTEGER nCells,nNodes,nCorners,nBodies,i,j,IELEM,IVERT
      CHARACTER*(*) plotfile
C
      INTEGER  ICYCL
      EXTERNAL ICYCL
C
      EQUIVALENCE(nCells,NELEM)
      EQUIVALENCE(nNodes,NPOIN)
*     EQUIVALENCE(nCorners,NOFVERT)
      nCorners = NOFVERT
C
      DATA nBodies /  0/
C
C
      OPEN(UNIT=1,FILE=PLOTFILE,FORM='FORMATTED',STATUS='UNKNOWN')
C
C	UNSTRUCTURED DATA FILE
C
      write (1,"(a22)") "unstructured grid data"
c
c..Needs "unstr" or "UNSTR" as first five characters
c
      write (1,*) nCells
      do i=1,nCells
         write (1,*) nCorners,(ICELNOD(j,i),j=1,nCorners)
      end do
      write (1,*) nNodes
*     write (1,*) U1(inf),U2(inf),U3(inf),U4(inf)
*     write (1,*)(U_infty(I),I=1,4)
c
c..Freestream state values
c
      IF(NOFVAR.EQ.1)THEN
         write (1,*)ONE
         do i=1,nNodes
            write (1,*) (POINT(j,i),j=1,DIM),ZROE(1,i)
         end do
      ELSEIF(NOFVAR.EQ.3)THEN
         write (1,*)ONE,ONE,ONE,ONE
         do i=1,nNodes
            write (1,*) (POINT(j,i),j=1,DIM),1.,ZROE(1,i),ZROE
     +      (2,i),ZROE(3,i)
         end do
      ELSEIF(NOFVAR.EQ.4)THEN
         write (1,*)ONE,ONE,0.d0,fren
         do i=1,nNodes
            write (1,*) (POINT(j,i),j=1,DIM),ZROE(1,i),ZROE(3,i),ZROE
     +      (4,i),ZROE(2,i)
         end do
      ELSE
         WRITE(6,*)" Don't know what to do with NOFVAR = ",NOFVAR
         CALL EXIT(1)
      ENDIF
C
      write (1,*) nBodies
*     do j=1,nBodies
*     write (1,*) nBodyFaces(j)
*      do i=1,nBodyFaces(j)
*write (1,*) node1(i),node2(i)
*end do
*     end do
c
c..Two nodes for face on body numbered as above
c
      write (1,*) nBoundaryFaces
      do i=1,nBoundaryFaces
         IELEM = BndPntr(1,i)
         IVERT = BndPntr(2,i)
         write (1,*) (ICELNOD(ICYCL(IVERT+j,NOFVERT),IELEM),j=1,DIM)
         write (16,*)ielem,BndPntr(3,i),
     + (ICELNOD(ICYCL(IVERT+j,NOFVERT),IELEM),j=1,DIM)
      end do
c
c..Two nodes for face on boundary numbered as above
c
      close(1)
c
      WRITE(NOUT,66)plotfile
   66 FORMAT(5X,'Dplot file WRITTEN to ..',A60/)
      return
      end
