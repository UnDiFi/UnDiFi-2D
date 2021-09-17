      subroutine xa112d(icelnod,icelcel,nelem,ibndfac,nbfac,
     +coor,iwork,lwork,ndim,npoin,iclr,ixy)
c
c     assigns b.c.s on faces
c
c
      implicit none

      integer nelem,nbfac,ndim,npoin,lwork,ixy
      integer icelnod(3,*),icelcel(3,*),ibndfac(3,nbfac)
      integer iwork(lwork,2),ic(2),iclr(2)
      integer ielem,ivert,iface,icycl,i,ifail,istmp,ipoin,jpoin,j
      double precision coor(ndim,npoin)
      double precision xi,xj,yi,yj,di,dj
      double precision eps
      parameter(eps=1.e-9)
!     parameter(eps=5.e-5)
!     parameter(eps=1.e-4)
!     parameter(eps=1.e-2)
      integer n,k,ibc
c
c     loop over all elements
c
      istmp = max(1,nbfac/20)
c     write(6,*)(coor(k,npoin),k=1,4)
c
c     loop over the two periodic surfaces
c
      do 5 k = 1, 2
         write(6,*)'periodic surface labeled ',iclr(k)
         ic(k)=0
         do 1 iface = 1, nbfac
            ibc = ibndfac(3,iface)
!           if( (iface/istmp)*istmp .eq. iface)write(6,*)iface
            if( ibc .NE. iclr(k) )goto 1
            ielem = ibndfac(1,iface)
            ivert = ibndfac(2,iface)
c
c     loop over the vertices of each element
c
            do 3 i = 1, 2
               ic(k) = ic(k) + 1
               if(ic(k).GT.(lwork/2))then
                  write(6,*)'Not enough room in work array, now = ',
     &lwork
                  call exit(3)
               endif
c
c     consider boundary elements only
c
               iwork(ic(k),k) = icelnod(icycl(ivert+i,3),ielem)
c
    3       continue
    1    continue ! loop over boundary faces
    5 continue ! loop over periodic faces
c
c     check that the same number of nodes is on both surfaces
c
      if( ic(1) .NE. ic(2) )then
          write(6,*)'Non matching periodic nodes (1) ',(ic(k),k=1,2)
          call exit(1)
      else
          write(6,*)'Found ',ic(1),' periodic nodes before removing
     + duplicated entries'
      endif

      do k = 1,2
         n = ic(k)
         call sortsp(n,iwork(1,k),ic(k))
      enddo
c
c     check that the same number of nodes is on both surfaces
c
      if( ic(1) .NE. ic(2) )then
          write(6,*)'Non matching periodic nodes (2) ',(ic(k),k=1,2)
          call exit(1)
      else
          write(6,*)'Found ',ic(1),' periodic nodes after removing
     + duplicated entries'
      endif
c
      open(14,FILE='corresp')
      write(14,*)ic(1)
      ifail = 0
      do 18 i = 1,ic(1)
         ipoin = iwork(i,1)
         xi = coor(1,ipoin)
         yi = coor(2,ipoin)
         do 8 j = 1,ic(2)
            jpoin = iwork(j,2)
            xj = coor(1,jpoin)
            yj = coor(2,jpoin)
            if( ((ixy.EQ.1) .AND. (abs(xi-xj) .LE. EPS)) 
     &      .OR.((ixy.EQ.2) .AND. (abs(yi-yj) .LE. EPS)))then
               write(14,*)ipoin,jpoin
               goto 18
            endif
    8    continue
         write(6,*)'Cannot find matching node for ',ipoin
         ifail = ifail + 1
   18 continue
      close(14)
      if(ifail.NE.0)then
         write(6,*)'Found ',ifail,' non matching entries'
         call exit(ifail)
      else
         write(6,*)'Nodes on bndry coloured ',iclr(2),
     &   ' shall be removed'
         write(6,*)'Corresponding nodes have been written in file "corre
     &sp"'
      endif 
c
      return
      end
