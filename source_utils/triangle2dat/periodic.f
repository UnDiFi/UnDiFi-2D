      subroutine periodic(XY,RV,icelnod,icelcel,irank,ndim,nofvert,
     +npoin,nelem,nbfac,iwork,npnod,ZROE,nofvar,ixy)
C
C    This routine renumbers a periodic mesh
C
C     NPNOD is half the # of periodic nodes
C     Nodes must be ordered as follows in the datafile
C     the lists of periodic nodes that have been deleted must appear last
C     and there must be a correspondance between these
C     two lists; the correspondance is established by the array MAP
C
C
C     +-------------------------+-------+-------+
C      ^                       ^       ^       ^ 
C      |                       |       |       |
C      |                       |       |     NPOIN
C      1                       |       |
C                              |   NPOIN-NPNOD
C                              | 
C                        NPOIN-2*NPNOD
C
C
      implicit none
C
C     .. Local Scalars ..
      INTEGER I,IA,IB,IC,IELEM,IFAIL,loc,ndim,nofvert,J,
     +        NBFAC,NELEM,NPNOD,NPOIN,nofvar,NGHOST,ixy,IXDRS
      CHARACTER*1 ANSW
      DOUBLE PRECISION PETSC_PI,ALPHA
      PARAMETER(PETSC_PI = 3.14159265358979323846264d0)
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION RV(npoin),XY(ndim,Npoin),ZROE(nofvar,*)
      INTEGER icelnod(nofvert,NELEM),IRANK(Npoin),
     +        icelcel(nofvert,NELEM),iwork(NPOIN)
      INTEGER INITXDR,IXDRCLOSE,IXDRIMAT,IXDRINT
C
      CALL IINIT(NPOIN,0,IWORK,1)
C
      OPEN (10,FILE='corresp')
C
C     1st col are the nodes to be deleted
C     entries in the 2nd col can be duplicated
C
      READ (10,FMT=*) NPNOD
      DO 12 I = 1,NPNOD
          READ (10,FMT=*) IA,IB
          IWORK(IA) = IB ! ia is mapped into ab
          write(6,*)ia,ib,xy(ixy,ia),xy(ixy,ib),' <-- these x/y coords s
     &hould be equal'
   12 CONTINUE
!     WRITE (6,FMT=*) 'ranking criterion',(iwork(i),i=1,npoin)
C     rank according to iwork so that flagged nodes will appear last
cnag  CALL M01DBF(IWORK,1,NPOIN,'Ascending',IRANK,IFAIL)
      CALL QSORTI(IRANK,NPOIN,IWORK)
cnag  WRITE (6,FMT=*) 'M01DBF (ranking) has returned IFAIL =',IFAIL
cnag  CALL M01ZBF(IRANK,1,NPOIN,IFAIL) ! CHECKS THE VALIDITY OF A PERMUTATION
cnag  WRITE (6,FMT=*) 'M01ZBF (validity) has returned IFAIL =',IFAIL
c
c     invert!
c
      CALL ICOPY(NPOIN,IRANK,1,IWORK(NPOIN+1),1) 
cnag  CALL M01ZAF(IWORK(NPOIN+1),1,NPOIN,IFAIL)
      CALL RNKIDX(IRANK(1),1,NPOIN,IFAIL)
      WRITE (6,FMT=*) 'RNKIDX has returned IFAIL =',IFAIL
c
c     WRITE (6,FMT=*) 'rank',(irank(i),i=1,npoin)
      OPEN(UNIT=13,FILE='rank.dat',FORM='formatted')
      do i = 1,npoin
      write(13,*)i,iwork(i),irank(i),iwork(i+npoin)
      enddo
      CLOSE(13)
cnag  CALL M01EBF(IWORK,1,NPOIN,IRANK,IFAIL) ! re-arranges according to a given rank
      CALL I4RANK(IWORK,1,NPOIN,IRANK,IFAIL) ! re-arranges according to a given rank

      WRITE (6,FMT=*) 'I4RANK (re-arranges) has returned IFAIL =',IFAIL
!     WRITE (6,FMT=*) 'work before',(iwork(i),i=1,npoin)
! update work, because the node addressed in iwork has changed position
      do i = 1,npoin
         ib = iwork(i)
         if( ib .EQ. 0 )then
! do nothing
         else
             iwork(i) = irank(ib)
         endif
      enddo
!     WRITE (6,FMT=*) 'work after',(iwork(i),i=1,npoin)
!
!     CALL M01ZBF(IRANK,1,NPOIN,IFAIL) ! CHECKS THE VALIDITY OF A PERMUTATION
!     iwork should be non-zero only for the periodic nodes 
C
c
      DO 18 IELEM = 1,NELEM
          DO 18 J = 1,nofvert
              IA = icelnod(J,IELEM)
              icelnod(J,IELEM) = IRANK(IA)
   18 CONTINUE
!     WRITE (6,FMT=*) 'rank',(irank(i),i=1,npoin)
cnag  CALL M01ZBF(IRANK,1,NPOIN,IFAIL) ! CHECKS THE VALIDITY OF A PERMUTATION
cnag  WRITE (6,FMT=*) 'M01ZBF (check) has returned IFAIL =',IFAIL
      DO 30 I = 1,ndim
          CALL DCOPY(NPOIN,XY(I,1),ndim,RV,1)
!         WRITE (6,FMT=*) 'rv',i,(rv(j),j=1,npoin)
cnag      CALL M01EAF(RV,1,NPOIN,IRANK,IFAIL)
          CALL R8RANK(RV,1,NPOIN,IRANK,IFAIL)
          WRITE (6,FMT=*) 'R8RANK (ranking coords) has returned IFAIL ='
     +,IFAIL
!         WRITE (6,FMT=*) 'rv',i,(rv(j),j=1,npoin)
          CALL DCOPY(NPOIN,RV,1,XY(I,1),ndim)
   30 CONTINUE
c
      DO 32 I = 1,nofvar
          CALL DCOPY(NPOIN,ZROE(I,1),nofvar,RV,1)
          CALL R8RANK(RV,1,NPOIN,IRANK,IFAIL)
          WRITE (6,FMT=*) 'R8RANK (ranking values) has returned IFAIL ='
     +,IFAIL
          CALL DCOPY(NPOIN,RV,1,ZROE(I,1),nofvar)
   32 CONTINUE
c
!     ic = 0
!     do ia = npoin-npnod+1,npoin
!        ib = iwork(ia)
!        write(6,*)'1',ia,ib,xy(2,ia),xy(2,ib)
!     enddo
C
C     write the mapping periodic ---> interior nodes
C
      OPEN(21,FILE='scratch') ! dump the mapping
      LOC = NPOIN - NPNOD
      DO 22 I = 1,NPNOD
          WRITE(21,*)IWORK(LOC+I)
   22 CONTINUE
      CLOSE(21)
      CALL IINIT(NPOIN,0,IWORK,1)
      OPEN(21,FILE='scratch') ! read the mapping
      DO 23 I = 1,NPNOD
          READ(21,*)IWORK(I)
   23 CONTINUE
      CLOSE(21)
c
      ic = 0
      do ia = npoin-npnod+1,npoin
      ic = ic+1
         ib = iwork(ic)
         write(6,*)'2',ia,xy(ixy,ia),ib,xy(ixy,ib),' <--- these x/y coor
     &ds should be equal'
      enddo
C
      IXDRS = INITXDR('file004.dat','w',.false.)
      IFAIL = IXDRINT(IXDRS,NPNOD)
      IFAIL = IXDRIMAT(IXDRS,NPNOD,IWORK)
      IFAIL = IXDRCLOSE(IXDRS)
C
      return
      END
