      subroutine hydra2(coor2,ndim,q2,nofvar,nnodes,ibndfac,wall,
     +nwall,tria,nofvert,ntria)
c
c     write mesh in the format for hydra
c
      implicit none
c
c  **************************************************************
c
      integer ndim,nnodes,ntria,nofvar,nwall,nofvert
      real*8 coor2(ndim,nnodes),q2(nofvar,nnodes),q1(6)
      real*8 ax_ch,x1,y1,x2,y2,xg,yg,anx,any,dx,dy
      real*8 ptot
      parameter(ptot=101956.)
c     parameter(ax_ch=58.85e-3)
      parameter(ax_ch=228.6e-3)
      real*8 cosalpha,sinalpha
      integer tria(nofvert,ntria),ibndfac(3,nwall),wall(2,nwall)
      integer i,l,iunit,iel,ibc,jbc,ielem,ivert,j,i1,i2,i3
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
c nwall                : number of boundary edges
c nofvert              : number of vertices of each cell (=3)
c ntria                : number of cells in the mesh
c arrays
c ------
c coor2                : x,y coordinates
c q2                   : array of variables
c ntri                 : connectivity table
c ibndfac(1,i)         : element (cell) the i-th boundary face belongs to
c ibndfac(2,i)         : vertex in front of the the i-th cell
c ibndfac(3,i)         : face label (same as mwl)
c
c
c   ***************************************************************
c
      write (6,FMT=*)' writing file for hydra'
      unit1=22
      open (unit1,file='hydra',form='formatted')
c
c     boundary edges
c
      do 6 i=1,nwall
          ielem = ibndfac(1,i)
          ivert = ibndfac(2,i)
          i1 = tria(icycl(ivert+1,nofvert),ielem)
          i2 = tria(icycl(ivert+2,nofvert),ielem)
          i3 = tria(ivert,ielem)
c
          x1 = coor2(1,i1)
          x2 = coor2(1,i2)
          y1 = coor2(2,i1)
          y2 = coor2(2,i2)
c
c se ho ben capito, la normale e` ottenuta routando in senso ANTIorario
c il vettore orientato che va dal primo al secondo vertice del bnd.
c
          anx = y2-y1
          any = x1-x2
c
c baricentro della faccia
c
          xg = 0.5*(x1+x2)
          yg = 0.5*(y1+y2)
c
c vettore che unisce il vertice opposto alla faccia al baricentro di q.ta
c
          dx = xg-coor2(1,i3)
          dy = yg-coor2(2,i3)
c
c if the dot product of the latter with the face normal is positive,
c then the order of the bnd nodes is o.k. otherwise swap
c
c     write (6,FMT=*)i,dx,dy
          if ( sign(1.d0,dx*anx+dy*any) .gt. 0.d0 )then
              wall(1,i) = i1
              wall(2,i) = i2
          else
              wall(2,i) = i1
              wall(1,i) = i2
          endif
c
c     write (6,FMT=*)(wall(j,i),j=1,2)
c
    6 continue
c
c  etichetta boundaries
c
      do 1 i=1,nwall
          ibc = ibndfac(3,i)
c
c converte i miei colori in quelli che vuole Sergio
c                      1 su inlet plane
c                      2 su outlet plane
c                      3 su lower periodic boundary
c                      4 su upper periodic boundary
c                      5 su blade surface (SS + PS)
C
C THERE ARE    80 BOUNDARY FACES COLOURED  1 inlet
C THERE ARE    80 BOUNDARY FACES COLOURED  2 outlet
C THERE ARE   160 BOUNDARY FACES COLOURED  3 periodic lower 
C THERE ARE   160 BOUNDARY FACES COLOURED  4 periodic upper
C THERE ARE   256 BOUNDARY FACES COLOURED  5 blade lower
C THERE ARE   256 BOUNDARY FACES COLOURED  6 blade upper
C
c
          if ( ibc .eq.  1)then
              jbc = 1
          elseif ( ibc .eq.  2)then
              jbc = 2
          elseif ( ibc .eq.  3)then
              jbc = 3
          elseif ( ibc .eq.  4)then
              jbc = 4
          elseif ( ibc .eq.  5)then
              jbc = 5
          elseif ( ibc .eq.  6)then
              jbc = 5
          else
	      write(6,*)'do not know what to do with ibc = ',ibc
              stop
          endif
          jbc = ibc
    1 continue
c
      write (unit1,FMT=*) nwall,ntria,nnodes
      write (unit1,FMT=*) ((wall(j,i),j=1,2),i=1,nwall)
      write (unit1,FMT=*) (ibndfac(3,j),j=1,nwall)
c
c        write the connectivity table
c
      write (unit1,*) ((tria(i,iel),i=1,nofvert),
     +    tria(1,iel),iel=1,ntria)
c
c     redimensionalise and write coords
c
      write (unit1,FMT=*) ((coor2(j,i)*ax_ch,j=1,ndim),1.d0,i=1,nnodes)
      q1(1) = 1.205
      q1(2) = 33.d0*cosd(51.5d0)
      q1(3) = 33.d0*sind(51.5d0)
      q1(4) = 0.d0
      q1(5) = 1.013e5
      q1(6) = .001
      write (unit1,FMT=*) (1.205,(33.d0*q2(j,i),j=2,3),0.d0,
     +ptot*q2(1,i),0.001,i=1,nnodes)
c
      close(unit1)
c
      return
      end
