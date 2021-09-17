      subroutine wascii(coor2,ndim,nnodes,ibndfac,
     +nwall,tria,neighb,nofvert,ntria)
c
c     write mesh in the format for hydra
c
      implicit none
c
c  **************************************************************
c
      integer ndim,nnodes,ntria,nofvar,nwall,nofvert
      real*8 coor2(ndim,nnodes)
      integer tria(nofvert,ntria),ibndfac(3,nwall),neighb(ntria,nofvert)
      integer i,l,iunit,iel,ibc,jbc,ielem,ivert,j,i1,i2
      integer unit1
      integer icycl
      external icycl
c
c     ******************************************************************
c     Unstructured mesh
c     ******************************************************************
c
c
c scalars
c ----------
c ndim                 : dimension of the space
c nofvar               : number of degrees of freedom per meshpoint
c nnodes               : number of meshpoint in the mesh
c nwall                : number of boundary edges /faces
c nofvert              : number of vertices of each cell (=3,4)
c ntria                : number of cells in the mesh
c arrays
c ------
c coor2                : x,y coordinates
c q2                   : array of variables
c ntri                 : connectivity table
c ibndfac(1,i)         : element (cell) the i-th boundary face belongs to
c ibndfac(2,i)         : vertex in front of the the i-th cell
c ibndfac(3,i)         : face label
c neighb(i,iel)        : is the elemt that shares with iel
c                        the face opposide vertex i
c                        if that is a bndry face, then either 
c                        neighb(i,iel) > ntria OR neighb(i,iel) < 1
c
c
c   ***************************************************************
c
      write (6,FMT=*)' writing file for Mario'
      unit1=22
      open (unit1,file='bump.ascii',form='formatted')
c
c     write nodal coords
c
      write (unit1,FMT=*)nnodes
      do 3 i=1,nnodes
          write (unit1,FMT=*) (coor2(j,i),j=1,ndim)
    3 continue
c
c        write the connectivity table
c
c
c
      write (unit1,FMT=*) ntria
 
      do 203 iel=1,ntria
          write (unit1,FMT=*) (tria(i,iel),i=1,nofvert)
c
c NOT referenced
c    &    (neighb(i,iel),i=1,nofvert)
  203 continue
c
      write (unit1,FMT=*) nwall
c
c     boundary edges (faces)
c
      do 6 i=1,nwall
          ielem = ibndfac(1,i)
          ivert = ibndfac(2,i)
c
c     print out also the global coords 
c     of the vertices of the face
c
              write (unit1,FMT=*) (ibndfac(j,i),j=1,3),
     +        ( tria(icycl(ivert+j,nofvert),ielem),j=1,ndim )
c
    6 continue
c
      close(unit1)
c
      return
      end
