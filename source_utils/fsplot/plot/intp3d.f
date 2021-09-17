      subroutine intp3d(idx,a,zroe,zintp,nofvar,nofvert)
c
c     $Id:$
c
c     interpolates data on a linear element
c     (triangle/tetrahedron) given the
c     area coordinates of the point belonging
c     to the element
c
c
c     data in input:
c     idx()  global index of the nofvert meshpoints
c            defining the simplex 
c     a      nofvert area coordinates of the interpolated point
c     zroe   dependent variables
c     nofvar leading dimension of zroe
c     output:
c     zintp     interpolated value(s) of the dependent variable  
c
      implicit none
      integer nofvar,ldz,nofvert
      double precision zroe(nofvar,*),zintp(*),a(*)
      integer idx(*)
      integer ip,i,j
      double precision q
      do j = 1,nofvar
         zintp(j) = 0.d0
      enddo
      do i = 1,nofvert
         ip = idx(i)
         q  = a(i)
         do j = 1,nofvar
            zintp(j) = zintp(j) + zroe(j,ip) * q
         enddo
      enddo
      end
