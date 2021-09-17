      subroutine readmesh(lbndfac,lcelcel,lcelnod,lcorg,ledgptr,
     &           lnodcod,lnodptr,lzroe,nbfac,npoin,nelem,nhole,nbpoin,
     &           nvt,nedge,fname,lia,lja,liclr,nclr,fndbnds)

!     implicit none
      include 'paramt.h'

!     this subroutine reads a mesh written by the triangle code
!     the arrays are allocated into a 1d array (called dstak)
!     using pseudo pointers as implemented in the port library
!
!     there are various triangle file that can be read:
!     node neigh edge poly
!     their format is explained in triangle/triangle.ps
!
!     ndim  is the space dimension = 2
!     nvt = ndim+1 is the number of vertices
!
!     nshmax  : max number of shock
!     npshmax : max number of shock points for each shock.
!     neshmax : max number of shock element for each shocks

!     .. local scalars ..
      integer lbndfac,lcelcel,lcelnod,lcorg,ledgptr,lnodcod,lnodptr,
     &        lzroe,lia,lja,liclr,nclr
      integer nbfac,npoin,nelem,nhole,nvt,nedge
      integer ifail,k
      logical fndbnds

!     .. local arrays ..
      integer kspace,lenstr

!     .. external subroutines ..
      external dinit,iinit,istkin,istkrl

      character*(*) fname

!     .. arrays in common ..
      double precision dstak(1)
      integer istak(1)
      common/cstak/dstak

!     .. equivalences ..
      equivalence (dstak(1),istak(1))

!     .. external functions ..
      integer  istkgt
      external istkgt

!     open log file
      open(8,file='log/readmesh.log')

!     ********** node file ***************
      call rtri(istak(ledgptr),nedge,istak(lbndfac),nbfac,istak(lcelnod)
     .  ,istak(lcelcel),nvt,nelem,dstak(lcorg),dstak(lzroe),
     .  istak(lnodcod),npoin,"node",0,fname)

!     allocate nodal and shock coordinates  corg(1:ndim,1:npoin+nshmax*npshmax)
      lcorg = istkgt((npoin+2*nshmax*npshmax)*ndim,4)
      call dinit((npoin+2*nshmax*npshmax)*ndim,zero,dstak(lcorg),1)

!     allocate node flag nodcode(1:npoin+nshmax*npshmax)
      lnodcod = istkgt(npoin+2*nshmax*npshmax,2)
      call iinit(npoin+2*nshmax*npshmax,0,istak(lnodcod),1)

!     allocate nodal values zroe(1:ndof,1:npoin+nshmax*npshmax)
!     write(6,*)'allocating for ndof = ',ndof
      lzroe = istkgt((npoin+2*nshmax*npshmax)*ndof,4)
      call dinit((npoin+2*nshmax*npshmax)*ndof,zero,dstak(lzroe),1)

      call rtri(istak(ledgptr),nedge,istak(lbndfac),nbfac,istak(lcelnod)
     .  ,istak(lcelcel),nvt,nelem,dstak(lcorg),dstak(lzroe),
     .   istak(lnodcod),npoin,"node",1,fname)

      write(8,*)'rtri node  --> ok'

!     ********** ele file ***************
      call rtri(istak(ledgptr),nedge,istak(lbndfac),nbfac,istak(lcelnod)
     .  ,istak(lcelcel),nvt,nelem,dstak(lcorg),dstak(lzroe),
     .   istak(lnodcod),npoin,"ele",0,fname)

!     allocate cell to node pointer, then read icelnod(1:nvt,1:nelem)
!     gives the nvt vertices of a given cell
      lcelnod = istkgt(nelem*nvt,2)
      call iinit(nelem*nvt,0,istak(lcelnod),1)

      call rtri(istak(ledgptr),nedge,istak(lbndfac),nbfac,istak(lcelnod)
     .  ,istak(lcelcel),nvt,nelem,dstak(lcorg),dstak(lzroe),
     .   istak(lnodcod),npoin,"ele",1,fname)

      write(8,*)'rtri ele  --> ok'

!     ********** neigh file ***************
      call rtri(istak(ledgptr),nedge,istak(lbndfac),nbfac,istak(lcelnod)
     .  ,istak(lcelcel),nvt,nelem,dstak(lcorg),dstak(lzroe),
     .   istak(lnodcod),npoin,"neigh",0,fname)

!     allocate cell to cell pointer, then read icelcel(1:nvt,1:nelem)
!     gives the nvt neighbouring triangles of a given cell
!     icelcel(i,ielem) is the triangle sharing the edge opposite
!     the i-th vertex of triangle ielem
      lcelcel = istkgt(nelem*nvt,2)
      call iinit(nelem*nvt,0,istak(lcelcel),1)

      call rtri(istak(ledgptr),nedge,istak(lbndfac),nbfac,istak(lcelnod)
     .  ,istak(lcelcel),nvt,nelem,dstak(lcorg),dstak(lzroe),
     .   istak(lnodcod),npoin,"neigh",1,fname)

      write(8,*)'rtri neigh --> ok'

!     ********** edge file ***************
      call rtri(istak(ledgptr),nedge,istak(lbndfac),nbfac,istak(lcelnod)
     .  ,istak(lcelcel),nvt,nelem,dstak(lcorg),dstak(lzroe),
     .   istak(lnodcod),npoin,"edge",0,fname)

!     allocate edge to vertex pointer, then read
      ledgptr = istkgt(nedge*3,2)
      call iinit(nedge*3,0,istak(ledgptr),1)

      call rtri(istak(ledgptr),nedge,istak(lbndfac),nbfac,istak(lcelnod)
     .  ,istak(lcelcel),nvt,nelem,dstak(lcorg),dstak(lzroe),
     .   istak(lnodcod),npoin,"edge",1,fname)

!     call x04eaf('general',' ',nvt,nelem,istak(lcelnod),nvt,
!    +            'mesh connectivity',ifail)
!     call x04eaf('general',' ',nvt,nelem,istak(lcelcel),nvt,
!    +            'cell neighbours',ifail)
!     call x04caf('general',' ',ndim,npoin,dstak(lcorg),ndim,
!    +            'nodal coordinates',ifail)

      write(8,*)'rtri edge  --> ok'

!     ********** poly file ***************
      call rtri(istak(ledgptr),nedge,istak(lbndfac),nbfac,istak(lcelnod)
     .  ,istak(lcelcel),nvt,nelem,dstak(lcorg),dstak(lzroe),
     .   istak(lnodcod),npoin,"poly",0,fname)

!     boundary structure; ibndfac(1:3,1:nbfac+nshmax*neshmax)
!     ibndfac(1,ibfac) and ibndfac(2,ibfac) are the endpoints
!     of bndry segment ibfac; ibndfac(3,ibfac) is the colour
!     ibndfac stores not only boundary edges, but also shock edges
!     shock edges are coloured 10
      lbndfac = istkgt(3*(nbfac+2*nshmax*neshmax),2) ! leave room for duplicated nodes
      call iinit(3*(nbfac+2*nshmax*neshmax),0,istak(lbndfac),1)

      call rtri(istak(ledgptr),nedge,istak(lbndfac),nbfac,istak(lcelnod)
     .  ,istak(lcelcel),nvt,nelem,dstak(lcorg),dstak(lzroe),
     .   istak(lnodcod),npoin,"poly",1,fname)

!     call x04eaf('general',' ',3,nbfac,istak(lbndfac),3,
!    +            'bndry pointer',ifail)

      k = lenstr(fname)
      write(8,*)'done reading ',fname(1:k),' triangle files'

!     lbndfac,lnodcod are pointers for integer arrays returned by rtri;
!     lnodptr,lia,lja,liclr are pointers for integer arrays created within setbndrynodeptr
!     the number of colours within the mesh is nclr
!     if you need to use these pointers elsewhere, they should be included in a common
      if(fndbnds)then
      call setbndrynodeptr(lbndfac,lnodcod,nbfac,npoin,
     &nbpoin,lnodptr,lia,lja,liclr,nclr)

!     the following subroutine writes the sequence of gridpoints of each
!     boundary (patch) separately within a different file bndry??.dat
!     it also illustrates how to use the various pointers
      call check(istak(lia),istak(lja),istak(liclr),nclr,dstak(lcorg),
     &ndim)
      endif
!     call exit(-1)
      return
      end
