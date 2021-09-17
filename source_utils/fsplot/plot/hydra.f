      subroutine hydra(coor2,ndim,q,q2,nofvar,nnodes,ibndfac,
     +nwall,tria,nofvert,ntria)
c
c     q(nofvar,*) stores the parameter vector (INPUT)
c
c     write mesh in the format for hydra
c
      implicit none
c
c  **************************************************************
c
      integer ndim,nnodes,ntria,nofvar,nwall,nofvert
      real*8 coor2(ndim,nnodes),q(nofvar,nnodes),q2(6,nnodes)
      real*8 ax_ch,x1,y1,x2,y2,xg,yg,anx,any,dx,dy,rho,ux,uy,pres
      real*8 zk
c     parameter(ax_ch=58.85e-3)
      parameter(ax_ch=228.6e-3,zk=1.4)
      integer tria(nofvert,ntria),ibndfac(3,nwall)
      integer i,l,iunit,iel,ibc,jbc,ielem,ivert,j,i1,i2
      integer unit1,ipde,npdes
      integer icycl
      external icycl
c
c     Sergio reference values
c
      real pref,uref,rref,lref
      parameter(pref=1.013e5,rref=1.226,lref=1.)
c
c     Aldo reference values
c
      real pref2,uref2,rref2,lref2,tref2
      parameter(pref2=124600.,tref2=318.,rref2=pref/(287.*tref2),
     +lref2=1.)
      uref=sqrt(pref/rref)
      uref2=sqrt(pref2/rref2)
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
      open (unit1,file='hydra56',form='formatted')
c
      write (unit1,FMT=*) nwall
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
c
          if ( ibc .eq.  1)then
              jbc = 3
          elseif ( ibc .eq.  2)then
              jbc = 2
          elseif ( ibc .eq.  3)then
              jbc = 4
          elseif ( ibc .eq.  4)then
              jbc = 1
          elseif ( ibc .eq.  5)then
              jbc = 5
          elseif ( ibc .eq.  6)then
!             jbc = 5
              jbc = 6
          else
              write(6,*)'do not know what to do with ibc = ',ibc
              stop
          endif
c         jbc = ibc
          write (unit1,FMT=*) jbc
    1 continue
c
c
c     boundary edges
c
      do 6 i=1,nwall
          ielem = ibndfac(1,i)
          ivert = ibndfac(2,i)
          i1 = tria(icycl(ivert+1,nofvert),ielem)
          i2 = tria(icycl(ivert+2,nofvert),ielem)
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
          dx = xg-coor2(1,ivert)
          dy = yg-coor2(2,ivert)
c
c if the dot product of the latter with the face normal is positive,
c then the order of the bnd nodes is o.k. otherwise swap
c
          if ( sign(1.d0,dx*anx+dy*any) .gt. 0.d0 )then
              write (unit1,FMT=*) i1,i2
          else
              write (unit1,FMT=*) i2,i1
          endif
c
c
    6 continue
c
c        write the connectivity table
c
      write (unit1,FMT=*) ntria
 
      do 203 iel=1,ntria
          write (unit1,*) (tria(i,iel),i=1,nofvert)
  203 continue
c
c     redimensionalise and write coords
c
      write (unit1,FMT=*)nnodes
      do 3 i=1,nnodes
          write (unit1,FMT=*) (coor2(j,i)*ax_ch,j=1,ndim),1.d0
    3 continue

      close(unit1)
c
c     recover (r,u,v,w,p) from Roe parameter vector
c
      do 4 i =1,nnodes
         pres = zk/(zk-1.)*(q(1,i)*q(2,i)-
     &   0.5*(q(3,i)**2+q(4,i)**2) )
         ux = q(1,i)*q(3,i)
         uy = q(1,i)*q(4,i)
         rho = q(1,i)*q(1,i)
         q2(1,i) = rho*rref2/rref
         q2(2,i) = ux*uref2/uref
         q2(3,i) = uy*uref2/uref
         q2(4,i) = 0.d0
         q2(5,i) = pres*pref2/pref
         q2(6,i) = 0.01
    4 continue
c
      npdes = 6
      open (unit=1,file='flow.aldo',form='formatted')
      write (1,*) nwall,ntria,nnodes
      write (1,*) ((q2(ipde,i),ipde=1,npdes),i=1,nnodes)
      close(unit=1)
c
      return
      end
