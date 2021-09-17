C
C  Subject 5.05: How do I find the intersection of a line and a plane?

C    If the plane is defined as:

C       a*x + b*y + c*z + d = 0

C   and the line is defined as:

C       x = x1 + (x2 - x1)*t = x1 + i*t
C       y = y1 + (y2 - y1)*t = y1 + j*t
C       z = z1 + (z2 - z1)*t = z1 + k*t

C   Then just substitute these into the plane equation. You end up
C   with:

C       t = - (a*x1 + b*y1 + c*z1 + d)/(a*i + b*j + c*k)

C   When the denominator is zero, the line is contained in the plane 
C   if the numerator is also zero (the point at t=0 satisfies the
C   plane equation), otherwise the line is parallel to the plane.
C---------------------------------------------------------------------
      double precision function plpt(a,b,c,d,x1,x2,y1,y2,z1,z2,ierr)
      implicit none
      integer ierr
      double precision a,b,c,d,x1,x2,y1,y2,z1,z2,i,j,k,t,denom
      double precision EPS,LARGE
      parameter (EPS=1.d-12,LARGE=1.d+38)
      i = x2 - x1
      j = y2 - y1
      k = z2 - z1
      denom = (a*i + b*j + c*k)
!     write(6,*)denom,EPS
      if( abs(denom) .LE. EPS )then
          plpt = LARGE
          ierr = 1
      else
          plpt = - (a*x1 + b*y1 + c*z1 + d)/denom
          ierr = 0
      endif
      return
      end
