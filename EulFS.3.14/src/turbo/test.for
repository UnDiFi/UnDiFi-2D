      program test
      include 'turb.com'
      double precision r,dr,s,smin,smax,ds
      integer i
      double precision c2,c3
      parameter(c2=0.7d0,c3=0.9d0) 
      double precision tfv1,tfw
      npts = 300
      goto 7
      dr = 2.50/real(npts)
      do i = 0,npts
         r = 0.+i*dr
         write(6,*)r,tfv1(r),tfw(r)
      enddo
    7 continue 
      smin = -5.d0
      smax = 5.d0
      ds = (smax-smin)/real(npts)
      do i = 0,npts
         s = smin+i*ds
         if(s .LT. -c2)then
            r = 1.d0 + (c2*c2+c3*s)/((c3-2.d0*c2)-s)
         else
            r = 1.d0 + s
         endif
         write(6,*)s,r
      enddo
      end

