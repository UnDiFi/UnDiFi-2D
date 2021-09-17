! Write the files *.node (of the backgound mesh) in triangle format

      subroutine wtri0(xy,
!    +                 ndim,
     +                 zroe,
!    +                 ndof,
     +                 nodcod,
     +                 npoin,
     +                 fname)

      implicit none
      include 'paramt.h'

!     .. scalar arguments ..
      integer npoin

!     .. array arguments ..
      character*(*) fname
      character     fwork*255
      double precision xy(ndim,*),zroe(ndof,*)
      integer nodcod(*)

!     .. local scalars ..
      integer ia,k,ipoin,ilist

!     .. external functions ..
      integer  icycl,lenstr
      external icycl,lenstr

!     .. intrinsic functions ..
      k = lenstr(fname)
      fwork(1:k+5) = fname(1:k)//".node"
      open(unit=19,file=fwork(1:k+5))
!     node files

!     * first line: <# of vertices> <dimension (must be 2)> <# of
!     * attributes> <# of boundary markers (0 or 1)>
!     * remaining lines: <vertex #> <x> <y> [attributes]
!     * [boundary marker]

       ilist=npoin

!     vertices must be numbered consecutively, starting from one or zero.

      write(19,*)ilist,ndim,ndof,1

      do  ipoin = 1, npoin
      if(nodcod(ipoin).eq.-1)then
             write(19,fmt=300)ipoin,(xy(ia,ipoin),ia=1,ndim),
     &       (zroe(ia,ipoin),ia=1,ndof),0
      elseif(nodcod(ipoin).eq.-2)then
             write(19,fmt=300)ipoin,(xy(ia,ipoin),ia=1,ndim),
     &       (zroe(ia,ipoin),ia=1,ndof),2
      else
             write(19,fmt=300)ipoin,(xy(ia,ipoin),ia=1,ndim),
     &       (zroe(ia,ipoin),ia=1,ndof),nodcod(ipoin)
      endif
      end do

      close(19)

! 300 format(i,6(1x,f),1x,i)
  300 format(i6,6(1x,e22.15),1x,i2)

      return
      end
