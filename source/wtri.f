! Write the files *.node *.poly (of modified meshes) in triangle format

      subroutine wtri(
     .              ibndfac,
     .              nbfac,
     .              nbfac_sh,
     .              icelnod,
     .              nofvert,
     .              xy,
     .              xysh,
     .              xyshu,
     .              xyshd,
     +              zroe,
     .              zroeshu,
     .              zroeshd,
     +              nodcod,
     +              nodcodsh,
     .              npoin,
     .              fname,
     +              nshocks,
     .              nshockpoints,
     .              nshocksegs,
     +              nphpoin)

!     write a node and a poly file in triangle fmt

      implicit none
      include 'paramt.h'

!     .. scalar arguments ..
      integer nofvert,npoin,nbfac,nbfac_sh,nphpoin,
     +nshocks,nshockpoints(*),nshocksegs(*)

!     .. array arguments ..
      character*(*) fname
      character     fwork*255
      double precision xy(ndim,*),
     +                 xysh(ndim,npshmax,*),
     +                 xyshu(ndim,npshmax,*),
     +                 xyshd(ndim,npshmax,*),
     +                 zroe(ndof,*),
     +                 zroeshu(ndof,npshmax,*),
     +                 zroeshd(ndof,npshmax,*)
      integer          ibndfac(3,*),
     +                 icelnod(nofvert,*),
     +                 nodcod(*),
     +                 nodcodsh(npshmax,*)

!     .. local scalars ..
      double precision x0,y0,dum1,dum2,dum3,dum4,dum,distmin3,distmin4,
     +                 distmin31,distmin41,
     +                 dum4x,dum2x,dum4y,dum2y,dum3x,dum3y,
     +                 dumx,dumy,
     +                 alpha,alpha3,alpha4
      integer ia,k,iface,ipoin,n,ibc,icheck,ilist,i,ish,
     +        ishplist(npshmax,nshmax),nholes,ihole,
     +        ipoinatdistmin3,ipoinatdistmin4,
     +        ipoinatdistmin31,ipoinatdistmin41,
     +        ish1,ish2,i1,i2,icount

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

! create the shock nodes list

       ilist=npoin+2*nshmax*npshmax

!      do ish=1,nshocks
!        do i=1,nshockpoints(ish)
!          ilist=ilist+1
!          ishplist(i,ish)=ilist
!        end do
!        do i=nshockpoints(ish)+1,2*nshockpoints(ish)
!          ilist=ilist+1
!          ishplist(i,ish)=ilist
!        end do
!      end do
!      write(*,*)'npoin',npoin
!      write(*,*)'ilist',ilist

!     vertices must be numbered consecutively, starting from one or zero.

      write(19,*)ilist,ndim,ndof,1

! write in the poly file all nodes. phantom nodes (the
! color of the phantom nodes is -1) are given the same
! coordinates as node 1, so that triangle will ignore them

      do  ipoin = 1, npoin
          if(nodcod(ipoin).ge.0)then
             write(19,fmt=300)ipoin,(xy(ia,ipoin),ia=1,ndim),
     &       (zroe(ia,ipoin),ia=1,ndof),nodcod(ipoin)

          elseif((nodcod(ipoin).eq.-1).or.(nodcod(ipoin).eq.-2))then
             write(19,fmt=300)ipoin,(xy(ia,1),ia=1,ndim),
     &       (zroe(ia,1),ia=1,ndof),-1
          endif
      end do

       icount=npoin
       do ish=1,nshmax
         do i=1,npshmax
           icount=icount+1
           if(nodcodsh(i,ish).eq.10)then
            write(19,fmt=300)icount,
     +      (xyshu(ia,i,ish),ia=1,ndim),
     +      (zroeshu(ia,i,ish),ia=1,ndof),
     +      nodcodsh(i,ish)
!          elseif(ish.eq.3.and.i.eq.2*nshockpoints(3)+1)then
!           write(19,fmt=300)npoin+(ish-1)*npshmax+i,
!    +      (xysh(ia,i,ish),ia=1,ndim),
!    +      (zroesh(ia,i,ish),ia=1,ndof),-99
!!   +      nodcodsh(i,ish)
           else
             write(19,fmt=300)icount,
     &       (xy(ia,1),ia=1,ndim),
     &       (zroe(ia,1),ia=1,ndof),-99
!    &       nodcodsh(i,ish)
           endif
         end do
      end do

       do ish=1,nshmax
         do i=1,npshmax
           icount=icount+1
           if(nodcodsh(i,ish).eq.10)then
            write(19,fmt=300)icount,
     +      (xyshd(ia,i,ish),ia=1,ndim),
     +      (zroeshd(ia,i,ish),ia=1,ndof),
     +      nodcodsh(i,ish)
!          elseif(ish.eq.3.and.i.eq.2*nshockpoints(3)+1)then
!           write(19,fmt=300)npoin+(ish-1)*npshmax+i,
!    +      (xysh(ia,i,ish),ia=1,ndim),
!    +      (zroesh(ia,i,ish),ia=1,ndof),-99
!!   +      nodcodsh(i,ish)
           else
             write(19,fmt=300)icount,
     &       (xy(ia,1),ia=1,ndim),
     &       (zroe(ia,1),ia=1,ndof),-99
!    &       nodcodsh(i,ish)
           endif
         end do
      end do

      close(19)

      fwork(1:k+5) = fname(1:k)//".poly"
      open(unit=19,file=fwork(1:k+5))
!     empty node list; this is just a poly file
      write(19,*)0,ndim,0,1
!     the number of boundary faces is:
!     nbfac(0) + 2*(nshockedges+1) + 4 -2
!     (-2) is linked at the two bndry edges of the background grid
!     which are deactivated
      icheck = 0
      do 80 iface = 1,nbfac_sh
         ibc = ibndfac(3,iface)
         if(ibc.gt.0) icheck = icheck+1
   80 continue

!     write(*,*)'nbfac:',nbfac
!     write(*,*)'nbfac_sh:',nbfac_sh
!     write(*,*)'icheck:',icheck

      write(19,*)icheck,1

      icheck = 0

      do 90 iface = 1,nbfac_sh
      ibc = ibndfac(3,iface)

! faces with ibc<0 are deactivated

      if(ibc.gt.0)then
         icheck = icheck+1
         write(19,fmt=*)icheck,(ibndfac(ia,iface),ia=1,3)
      endif
   90 continue

! compute the number of holes

       nholes=0
       do ish=1,nshocks
         nholes=nholes+nshockpoints(ish)-2
       end do

! warning: modified for q1d
! adds a hole point in the innermost circle
       write(19,*)nholes+naddholes
!      write(19,*)nholes

       ihole=0
       do ish=1,nshocks
         do i=2,nshockpoints(ish)-1
          ihole=ihole+1
          x0 = xysh(1,i,ish)
          y0 = xysh(2,i,ish)
          write(19,*)ihole,x0,y0
         enddo
       enddo
       do i=1, naddholes
        write(19,*)ihole+i,caddhole(1,i),caddhole(2,i)
       enddo

! triple point n.1
!      x0 = xysh(1,nshockpoints(1),1)
!      y0 = xysh(2,nshockpoints(1),1)
!      write(19,*)nholes+1,x0,y0
!
!!     here we assume that:
!!     1) shock points are stored btw
!!
!!     we skip the first and last (ia=2, nholes-1)
!!
!      k = (npoin-2*nholes) +1
!      n = (npoin-  nholes) +1
!      do ia = 2, nholes-1
!         k=k+1
!         n=n+1
!         x0 = 0.5d0*( coor(1,k) + coor(1,n) )
!         y0 = 0.5d0*( coor(2,k) + coor(2,n) )
!      write(19,*)ia-1,x0,y0
!      enddo
      close(19)
! 300 format(i7,6(1x,e12.6),1x,i3)
  300 format(i7,6(1x,e22.15e3),1x,i3)

      return
      end
