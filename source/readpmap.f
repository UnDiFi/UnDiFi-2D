!     in the main define lpmap among pointers and npnod(0:2) among variables
!
!     integer lbndfac(0:2),lcelcel(0:2),lcelnod(0:2),lcorg(0:2),
!    & lnodcod(0:2),lzroe(0:2),ledgptr(0:2),lnodptr(0:2),
!    & lshnor,lxyshold,lxyshnew,lwork, lpmap(0:2)
!
!      integer npnod(0:2)
!
!     after readmesh()
!
!     call readmap(npoin(0),npnod(0),lpmap(0))
!
!     when the whole vector is needed use:
!
!     call subroutine( ...,istak(lpmap(0)),....

      subroutine readpmap(nitems,npnod,lpmap)
      implicit none
      integer nitems,npnod,lpmap 

!<     nitems (input) the nof of gridpoints incls. both sets of periodic nodes
!<     npnod  (output) the nof periodic gridpoints in only one of the two sets
!<     lpmap  (output) a pointer to an integer array of length nitems
!<     for a given meshpoint i
!<     j = lpmap(i) <> 0 is the periodic meshpoint corresponding to i
!<     j = lpmap(i) = 0 means that i is not among the periodic meshpoints

      logical lflag
      integer i,j,n

!     .. arrays in common ..
      double precision dstak(1)
      integer istak(1)
      common/cstak/dstak

!     .. equivalences ..
      equivalence (dstak(1),istak(1))

!     .. external functions ..
      integer  istkgt
      external istkgt

!     read the file with the index of all nodes and in case 
!     the corresponding node on the periodic boundary
      inquire(file="pnodes0.dat",exist=lflag)
      if(.not.lflag)then
!     if the file is absent, build a table with all zero (i.e. with any periodic node)
         n=nitems
         lpmap = istkgt(n,2)
         do i = lpmap,lpmap+n-1
           istak(i)=0
         enddo
         return
      endif 
      open(13,file="pnodes0.dat") 
      read(13,*)n,npnod
      if(n.ne.nitems)then
           write(6,*)'the nof meshpoints in the dataset and in pnodes0.d
     &at do not match'
           call exit(13)
      endif
      lpmap = istkgt(n,2)
      do i = lpmap,lpmap+n-1
         read(13,*)j,istak(i)
      enddo
      close(13)
      return
      end
