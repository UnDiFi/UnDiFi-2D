      subroutine seg3d(head,nubo,nu,ns,nt,nnsg)

      implicit none
c
c     a code for computing the edges of a tetrahedral mesh;
c     taken (or stolen) from a code by B.Mohammadi
c     upon return:
c     nubo lists the forming nodes of all the edges
c     head returns the degree of each vertex

      integer ns,nt,nnsg
      integer is,jt,i,js1,js2,is1,is2,nseg,k,iseg
      integer head(ns),nubo(2,nnsg),nu(4,nt)
      integer nex(6),nor(6)

      data nex,nor/1,1,1,2,2,3,2,3,4,3,4,4/

      do 1 is = 1, ns
         head(is) = 0
    1 continue

      nseg = 0

      do 1000 jt = 1, nt
         do 500 k = 1,6
            js1 = nu(nor(k),jt)
            js2 = nu(nex(k),jt)
            is1 = nu(nor(k),jt)
            is2 = nu(nex(k),jt)
            if(is1.lt.is2) then
               is2 = nu(nor(k),jt)
               is1 = nu(nex(k),jt)
            endif
C     is IS2 the neighbor of is1 ?
C     head of the list of edge with the same first vertex
C     check the existance of the edge (is2,is1)
            is = head(is1)
  100 continue
            if(is.gt.0)then
               if( is2.eq.nubo(2,is) )then
c              we find the edge is1,is2
                   goto 200
               else
                   is = nubo(1,is)
c              try the next edge with the same first vertex
                   goto 100
               endif
            endif
c           new edge
            nseg = nseg + 1
            if(nseg.gt.nnsg) then
               write(6,*) 'seg3d:increase nnsg : ',nnsg
            endif
            is = nseg
c           link at the begin of list:
c                the vertices is1
            nubo(1,nseg)=head(is1)
            nubo(2,nseg)=is2 
            head(is1)=nseg
c
  200 continue
  500 continue
 1000 continue

      do 2600 is1 = 1,ns
	 iseg=head(is1)
 2500 continue
      if(iseg.gt.0)then
	 is=nubo(1,iseg)
	 nubo(1,iseg) = is1
	 iseg = is
	 goto 2500
      endif
 2600 continue
       
      do 3 is = 1, ns
         head(is) = 1
    3 continue


       do 5 i = 1, nseg
          is1 = nubo(1,i)
          is2 = nubo(2,i)
          head(is1)=head(is1)+1
          head(is2)=head(is2)+1
    5 continue
       return
       end
