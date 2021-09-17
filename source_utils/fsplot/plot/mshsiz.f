      subroutine mshsiz(icelnod,nofvert,nelem,coor,ndim)
      implicit none
      integer nofvert,nelem,ndim
      integer icelnod(nofvert,nelem)
      real*8 coor(ndim,*)
      real*8 t(12),pc(3)
      real*8 r1min,r1max,r2min,r2max,r1avg,r2avg,wg
      real*8 std1,std2,r1sum,r2sum,r1,r2,ar
      integer i,ipoin,ielem,ivert,j
      r1min = 1.d+38
      r2min = 1.d+38
      r1max = 0.d+00
      r2max = 0.d+00
      r1sum = 0.d0
      r2sum = 0.d0
      r1avg = 0.d0
      r2avg = 0.d0
      wg = 1.d0/real(nelem)
caldo open(120,FILE='fort.12')
caldo write(120,*)nelem,1
      do ielem = 1, nelem
         do ivert = 1,nofvert
            j = (ivert-1)*ndim
            ipoin = icelnod(ivert,ielem)
            do i= 1,ndim
               t(j+i) = coor(i,ipoin)
            enddo 
         enddo
         if(ndim.EQ.2)then
            call triangle_inradius_2d(t,r1)
            call triangle_circumcircle_2d(t,r2,pc)
!           call triangle_diameter_2d(t,r2)
!           r2 = r2*0.5d0
            ar = r2/(2.d0*r1)
caldo       write(120,*)ielem,ar

         elseif(ndim.EQ.3)then
            call tetrahedron_insphere_3d(t,r1,pc)
            call tetrahedron_circumsphere_3d(t,r2,pc)
         else
            write(6,*)'Smthg. wrong with ndim = ',ndim
            call exit(ndim)
         endif
         r1avg = r1avg + wg * r1
         r2avg = r2avg + wg * r2
         r1sum = r1sum + r1 * r1
         r2sum = r2sum + r2 * r2
         r1min = min(r1,r1min)
         r1max = max(r1,r1max)
         r2min = min(r2,r2min)
         r2max = max(r2,r2max)
      enddo ! ielem
      std1 = sqrt(r1sum * wg - r1avg*r1avg)/r1avg
      std2 = sqrt(r2sum * wg - r2avg*r2avg)/r2avg
      write(6,100)'Inradius     (min/avg/rel.std./max) : ',r1min,r1avg,
     &std1,r1max
      write(6,100)'Circumradius (min/avg/rel.std./max) : ',r2min,r2avg,
     &std2,r2max
caldo close(120)
  100 FORMAT(A,4(1X,E12.6))
      end
