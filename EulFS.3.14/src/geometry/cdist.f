      subroutine c_distij ( node_xyz, ndim, node_boundary, node_num, 
     &distij, icelnod, nofvert, nelem )
c
c     $Id: cdist.f,v 1.2 2013/07/17 10:30:42 abonfi Exp $
c
c     node_num = NPOIN + NGHOST + NPNOD 
c
      implicit none
c
      INCLUDE 'constants.h'
c
      integer ndim,node_num,nofvert,nelem
      integer node_boundary(node_num)
      integer icelnod(nofvert,nelem) 
      double precision node_xyz(ndim,node_num),distij(nelem)
      integer i,j,ielem,ivert,ipoin,idx
      double precision help,dist
      double precision g(3)
      logical verbose
      parameter(verbose=.FALSE.)
!
!      find the distance of the current node from the nearest boundary
!
!     write(6,*)'Calculating distance in ',node_num,' gridpoints'
      if(verbose)open(24,FILE="distance.txt") 
      do ielem = 1, nelem ! loop over all elements
!
!     pick up the coords of the centroid of cell ielem
!
         do i = 1, ndim
            g(i) = ZERO
         enddo
         do ivert = 1, nofvert ! loop over the vertices
            ipoin = icelnod(ivert,nelem+ielem) ! for periodic grids
            do i = 1, ndim 
               g(i) = g(i) + node_xyz(i,ipoin)
            enddo ! loop over dimensions
         enddo ! loop over vertices
         do i = 1, ndim
            g(i) = g(i)/real(nofvert)
         enddo
c
c the centroid is now in g
c
          help = 1.e+38
          idx = -1
          do j = 1, node_num ! loop over all nodes, including ghost and periodic
             if( node_boundary(j) .NE. 0 )then ! but pick up only the boundary ones
                dist = ZERO
                do i = 1,ndim
                   dist = dist + (node_xyz(i,j)-g(i))**2
                enddo
                dist = sqrt(dist)
                if( dist .LT. help )then
                    help = dist
                    idx = j
                endif 
                help = min(help,dist)
!               write(24,*)ielem,j,dist,help
             endif
          enddo ! end loop over all (boundary) nodes
          if(abs(help).LT.1.d-8)then
             distij(ielem) = 1.d-8
          else
             distij(ielem) = help
          endif
          if(verbose)write(24,FMT="(I6,3(1X,E12.4),1X,I5)")ielem,
     &(g(i),i=1,ndim),distij(ielem),idx
      enddo ! loop over ielem
      if(verbose)close(24)
      write(6,*)'done ! ' 
      do ielem = 1,nelem
         distij(ielem) = ONE/max(distij(ielem),0.01d0)
      enddo 
      return
      end
