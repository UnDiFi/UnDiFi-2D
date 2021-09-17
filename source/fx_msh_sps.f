! This routine:
!
!     1) identifies the bndry segments cut by the shock lines
!     2) modifies the bndry structure (i.e. the pointer
!        IBNDFAC)
!        2a) to include the shock segments
!        2b) to account for the bndry segments of the background mesh
!            cut by the shock

      subroutine fx_msh_sps(
     .                ibndfac,
     .                nodcod,
     .                nbfac,
     .                nbfac_sh,
     .                nvt,
     .                nelem,
     .                xy,
     .                xysh,
     .                xyshu,
     .                xyshd,
     .                npoin,
     .                nshocks,
     .                nshockpoints,
     .                nshockedges,
     .                nspecpoints,
     .                typespecpoints,
     .                shinspps,
     .                ispclr)

      implicit none
      include 'paramt.h'

      integer nelem,npoin,nvt,nbfac,nbfac_sh,nbfac_new
      integer nshocks,nshockedges(nshmax),nshockpoints(nshmax)
      integer nspecpoints,shinspps(2,5,*),ispclr(5,*),nodcod(npoin)

!     .. array arguments ..
      double precision   xy(ndim,*),
     .                   xysh(ndim,npshmax,*),
     .                   xyshu(ndim,npshmax,*),
     .                   xyshd(ndim,npshmax,*)
      integer            ibndfac(3,*)

!     .. array arguments ..
!     character*(*) fname
      character*5 typespecpoints(*)

!     .. local scalars ..
      double precision x0,y0,help,s1,s2,x1,y1,x2,y2,x3,y3,x4,y4
      integer i,iedg1,iedg2,ifail,ipoin(4),i1,i2,ibc,
     +        ibfac,ish,ish1,ip1,isppnts,ish2,ip2,
     +        ilist,ibf,
     +        ishplistu(npshmax,nshmax),
     +        ishplistd(npshmax,nshmax),
     +        idum1,idum2,j1,j2,
     +        ishel1

!     .. external functions ..
      integer findbedg

!     open log file
      open(8,file='log/fx_msh_sps.log')

!     check the intersections of the shock lines
!
!     do i = 1, nbfac
!        ibc=ibndfac(3,i)
!        if(ibc.lt.0)ibndfac(3,i)=iabs(ibc)
!     enddo
!
!     call x04eaf('general',' ',3,nbfac,ibndfac,3,
!    +            'bndry pointer in chkbndry',ifail)
!     pause
!
!      do ish=1,nshocks
!      s(ish,1) = -100.d0
!      s(ish,2) = -100.d0
!      s(ish,3) = -100.d0
!      s(ish,4) = -100.d0
!      ipoin(1) = 1
!      x0 = xyshd(1,ipoin(1),ish)
!      y0 = xyshd(2,ipoin(1),ish)
!      iedg(ish,1) = findbedg(xy,ndim,ibndfac,nbfac,x0,y0,s(ish,1))
!      write(8,*)'s(1) ' ,s(ish,1),x0,y0,iedg(ish,1)
!      if( iedg(ish,1) .eq. -1 )then
!          write(8,*)'failed matching 1st shock point of the shock n.'
!     +,ish
!!         stop
!      else
!          write(8,*)'shockpoint (1) ' ,x0,y0,' falls within ',
!     &(ibndfac(i,iedg(ish,1)),i=1,2)
!      endif
!      ipoin(2) = ipoin(1)
!      x0 = xyshu(1,ipoin(2),ish)
!      y0 = xyshu(2,ipoin(2),ish)
!      iedg(ish,2) = findbedg(xy,ndim,ibndfac,nbfac,x0,y0,s(ish,2))
!      write(8,*)'shockpoint (2) ' ,x0,y0,' falls within ',
!     &(ibndfac(i,iedg(ish,2)),i=1,2)
!     write(8,*)'s(2) ' ,s(ish,2),x0,y0,iedg(ish,2)
!      if( iedg(ish,2) .eq. -1 )then
!       write(8,*)'failed matching 2nd shock point of the shock n.',ish
!!         stop
!      endif

!     check if the shock crosses a bndry point

!      do i =1,2
!      if( s(ish,i) .lt. 0.d0 .or. s(ish,i) .gt. 1.d0 )then
!          write(8,*)'s(',i,') out of bounds',s(ish,i)
!!         stop
!      endif
!      enddo
!      if( iedg(ish,2) .ne. iedg(ish,1) )then
!          write(8,*)'shock points (1) (2) not on the same bndry edge'
!          write(8,*)iedg(ish,1),iedg(ish,2)
!          write(8,*)s(ish,1),s(ish,2)
!          stop
!      endif

!     the other end of the shock line

!      ipoin(3) = nshockpoints(ish)
!      x0 = xyshd(1,ipoin(3),ish)
!      y0 = xyshd(2,ipoin(3),ish)
!      iedg(ish,3) = findbedg(xy,ndim,ibndfac,nbfac,x0,y0,s(ish,3))
!      write(8,*)'shockpoint (3) ' ,x0,y0,' falls within ',
!     &(ibndfac(i,iedg(ish,3)),i=1,2)
!      write(8,*)'s(3) ' ,s(ish,3),x0,y0,iedg(ish,3)
!      if( iedg(ish,3) .eq. -1 )then
!       write(8,*)'failed matching 3rd shock point of the shock n.',ish
!      endif
!      ipoin(4) = ipoin(3)
!      x0 = xyshu(1,ipoin(4),ish)
!      y0 = xyshu(2,ipoin(4),ish)
!      iedg(ish,4) = findbedg(xy,ndim,ibndfac,nbfac,x0,y0,s(ish,4))
!      write(8,*)'shockpoint (4) ' ,x0,y0,' falls within ',
!     &(ibndfac(i,iedg(ish,4)),i=1,2)
!      write(8,*)'s(4) ' ,s(ish,4),x0,y0,iedg(ish,4)
!      if( iedg(ish,4) .eq. -1 )then
!       write(8,*)'failed matching 4th shock point of the shock n.',ish
!      endif

!     check if the shock crosses a bndry point

!      do i =3,4
!      if( s(ish,i) .lt. 0.d0 .or. s(ish,i) .gt. 1.d0 )then
!          write(8,*)'s(',i,') out of bounds',s(ish,i)
!!         stop
!      endif
!      end do
!      if( iedg(ish,3) .ne. iedg(ish,4) )then
!          write(8,*)'shock points (2) (3) not on the same bndry edge'
!          write(8,*)iedg(ish,3),iedg(ish,4)
!         stop
!      endif
!      end do

*****************************************
! TODO: check how to modify this part

!      create the shock nodes list
       ilist=npoin
       do ish=1,nshmax
         do i=1,npshmax
           ilist=ilist+1
           ishplistu(i,ish)=ilist
         end do
       end do

       do ish=1,nshmax
         do i=1,npshmax
           ilist=ilist+1
           ishplistd(i,ish)=ilist
         end do
       end do

!     correction for the case of a shock splitted in two connetted parts
!     the first point of shock 2 coincides with the first point of shock 1+4
!     with this correction the first shock point of shock 1+4 disappeares from connectivity
!     ishplist(1,1+4)=ishplist(1,2)
!     ishplist(1+nshockpoints(1+4),1+4)=ishplist(1+nshockpoints(2),2)

! triple point 1
!     ishplistu(nshockpoints(3),3)=ishplistu(nshockpoints(1),1)
!     ishplistu(nshockpoints(2),2)=ishplistd(nshockpoints(1),1)
!     ishplistu(nshockpoints(4),4)=ishplistd(nshockpoints(2),2)
!     ishplistd(nshockpoints(4),4)=ishplistd(nshockpoints(3),3)

****************************************

!     create shock edges (downstream)
      ibfac=nbfac
      do ish=1, nshocks
       nshockedges(ish) = nshockpoints(ish) - 1
       do i = 1, nshockedges(ish)
         ibfac=ibfac+1
         ibndfac(1,ibfac)=ishplistd(i,ish)
         ibndfac(2,ibfac)=ishplistd(i+1,ish)
         ibndfac(3,ibfac)=10
       end do

       do i = 1, nshockedges(ish)
         ibfac=ibfac+1
         ibndfac(1,ibfac)=ishplistu(i,ish)
         ibndfac(2,ibfac)=ishplistu(i+1,ish)
         ibndfac(3,ibfac)=10
       end do
      end do

!     create the 2 extra shock edges needed to close the shock hole
!
!      ibfac=ibfac+1
!      ibndfac(1,ibfac)=ishplistd(nshockpoints(ish),ish)
!      ibndfac(2,ibfac)=ishplistu(nshockpoints(ish),ish)
!      ibndfac(3,ibfac)=10
!
!      ibfac=ibfac+1
!      ibndfac(1,ibfac)=ishplistd(1,ish)
!      ibndfac(2,ibfac)=ishplistu(1,ish)
!      ibndfac(3,ibfac)=10
!
!     nbfacnew=ibfac

      do isppnts=1,nspecpoints
        if(typespecpoints(isppnts).eq.'IPX'.or.
     +     typespecpoints(isppnts).eq.'IPY'.or.
     +     typespecpoints(isppnts).eq.'OPX'.or.
     +     typespecpoints(isppnts).eq.'OPY'.or.
     +     typespecpoints(isppnts).eq.'FWP'.or.
!    +     typespecpoints(isppnts).eq.'PC' .or.
     +     typespecpoints(isppnts).eq.'WPNRX'.or.
     +     typespecpoints(isppnts).eq.'WPNRY')then

          ish1 = shinspps(1,1,isppnts)
          i    = shinspps(2,1,isppnts)-1
          ip1  = 1+i*(nshockpoints(ish1)-1)

          x0 = xyshd(1,ip1,ish1)
          y0 = xyshd(2,ip1,ish1)

caldo
!         do ibf=1,nbfac
!          i1 = ibndfac(1,ibf)
!          i2 = ibndfac(2,ibf)
!          x1 = xy(1,i1)
!          y1 = xy(2,i1)
!          x2 = xy(1,i2)
!          y2 = xy(2,i2)
!         write(8,*)ibf,(y0-y1)*(x2-x1)-(x0-x1)*(y2-y1)
!         enddo
caldo

          iedg1 = findbedg(xy,ndim,ibndfac,nbfac,x0,y0,s1)
          write(8,*)'typespecpoints:',typespecpoints(isppnts)
          write(8,*)'s(1) ' ,s1,x0,y0,iedg1
          if( iedg1 .eq. -1 )then
           write(8,*)'failed matching 1st shock point of the shock n.'
     +                ,ish1
           stop
          else
          write(8,*)'shockpoint (1) ' ,x0,y0,' falls within ',
     &   (ibndfac(i,iedg1),i=1,2)
          endif
          x0 = xyshu(1,ip1,ish1)
          y0 = xyshu(2,ip1,ish1)

!         do ibf=1,nbfac
!          i1 = ibndfac(1,ibf)
!          i2 = ibndfac(2,ibf)
!          x1 = xy(1,i1)
!          y1 = xy(2,i1)
!          x2 = xy(1,i2)
!          y2 = xy(2,i2)
!         write(8,*)ibf,(y0-y1)*(x2-x1)-(x0-x1)*(y2-y1)
!         enddo

          iedg2 = findbedg(xy,ndim,ibndfac,nbfac,x0,y0,s2)

          write(8,*)'shockpoint (2) ' ,x0,y0,' falls within ',
     &   (ibndfac(i,iedg2),i=1,2)
         write(8,*)'s(2) ' ,s2,x0,y0,iedg2
         if( iedg2 .eq. -1 )then
          write(8,*)'failed matching 2nd shock point of the shock n.'
     +    ,ish1
          stop
         endif

!        check if the shock crosses a bndry point
         if( s1 .lt. 0.d0 .or. s1 .gt. 1.d0.or.
     +       s2 .lt. 0.d0 .or. s2 .gt. 1.d0 )then
          write(8,*)'s(',i,') out of bounds',s1,s2
          stop
        endif
        if( iedg2 .ne. iedg1 )then
          write(8,*)'shock points (1) (2) not on the same bndry edge'
          write(8,*)iedg1,iedg2
          write(8,*)s1,s2
          stop
        endif

!        split the existing edges of the background mesh
!        if the first point of the shock  is at the boundary
         write(8,*)'**********************'
         write(8,*)'shock:',ish1
         write(8,*)iedg1
         write(8,*)'**********************'
         if(iedg1.gt.0)then

            i   = iedg1
            i1  = ibndfac(1,i)
            i2  = ibndfac(2,i)
            ibc = ibndfac(3,i)

! ************************************************************

            if(nodcod(i1).lt.0.or.nodcod(i2).lt.0)then
!             write(*,*) 'stopped'
!             write(*,*)'nodcod(1):',nodcod(i1)
!             write(*,*)'nodcod(2):',nodcod(i2)
!             write(*,*)'s1:',s1
!             write(*,*)'s2:',s2
!             write(*,*)
            if(nodcod(i1).lt.0.d+0.and.s1.lt.s2)then
              ibndfac(1,i)=ishplistu(ip1,ish1)
              idum1=i1
              idum2=ishplistd(ip1,ish1)
            elseif(nodcod(i1).lt.0.d+0.and.s1.gt.s2)then
              ibndfac(1,i)=ishplistd(ip1,ish1)
              idum1=i1
              idum2=ishplistu(ip1,ish1)
            elseif(nodcod(i2).lt.0.d+0.and.s1.lt.s2)then
              ibndfac(2,i)=ishplistd(ip1,ish1)
              idum1=i2
              idum2=ishplistu(ip1,ish1)
            elseif(nodcod(i2).lt.0.d+0.and.s1.gt.s2)then
              ibndfac(2,i)=ishplistu(ip1,ish1)
              idum1=i2
              idum2=ishplistd(ip1,ish1)
             endif
             do ibf=1,nbfac
               if(ibndfac(1,ibf).eq.idum1)ibndfac(1,ibf)=idum2
               if(ibndfac(2,ibf).eq.idum1)ibndfac(2,ibf)=idum2
              enddo

            else
! *****************************************************
            ibndfac(3,i)=-ibc
            write(8,*)'removing background edge ' ,i,i1,i2

!          create 2 new edges at boundary
           ibfac=ibfac+1
           ibndfac(3,ibfac)=ibc
           if( s1 .lt. s2 )then
             ibndfac(1,ibfac)=i1
             ibndfac(2,ibfac)=ishplistd(ip1,ish1)
           else
            ibndfac(1,ibfac)=i1
            ibndfac(2,ibfac)=ishplistu(ip1,ish1)
           endif

           ibfac=ibfac+1
           ibndfac(3,ibfac)=ibc
           if( s1 .lt. s2 )then
             ibndfac(1,ibfac)=ishplistu(ip1,ish1)
             ibndfac(2,ibfac)=i2
           else
            ibndfac(1,ibfac)=ishplistd(ip1,ish1)
            ibndfac(2,ibfac)=i2
           endif
          endif
         endif

        elseif(typespecpoints(isppnts).eq.'PC') then

! point 1

          ish1 = shinspps(1,1,isppnts)
          i    = shinspps(2,1,isppnts)-1
          ip1  = 1+i*(nshockpoints(ish1)-1)

          x0 = xyshd(1,ip1,ish1)
          y0 = xyshd(2,ip1,ish1)
caldo
          do ibf=1,nbfac
           i1 = ibndfac(1,ibf)
           i2 = ibndfac(2,ibf)
           x1 = xy(1,i1)
           y1 = xy(2,i1)
           x2 = xy(1,i2)
           y2 = xy(2,i2)
          write(8,*)ibf,(y0-y1)*(x2-x1)-(x0-x1)*(y2-y1)
          enddo
caldo

          iedg1 = findbedg(xy,ndim,ibndfac,nbfac,x0,y0,s1)
          write(8,*)'typespecpoints:',typespecpoints(isppnts)
          write(8,*)'s(1) ' ,s1,x0,y0,iedg1
          if( iedg1 .eq. -1 )then
           write(8,*)'failed matching 1st shock point of the shock n.'
     +                ,ish1
           stop
          else
          write(8,*)'shockpoint (1) ' ,x0,y0,' falls within ',
     &   (ibndfac(i,iedg1),i=1,2)
          endif
          x0 = xyshu(1,ip1,ish1)
          y0 = xyshu(2,ip1,ish1)

          do ibf=1,nbfac
           i1 = ibndfac(1,ibf)
           i2 = ibndfac(2,ibf)
           x1 = xy(1,i1)
           y1 = xy(2,i1)
           x2 = xy(1,i2)
           y2 = xy(2,i2)
          write(8,*)ibf,(y0-y1)*(x2-x1)-(x0-x1)*(y2-y1)
          enddo

          iedg2 = findbedg(xy,ndim,ibndfac,nbfac,x0,y0,s2)

          write(8,*)'shockpoint (2) ' ,x0,y0,' falls within ',
     &   (ibndfac(i,iedg2),i=1,2)
         write(8,*)'s(2) ' ,s2,x0,y0,iedg2
         if( iedg2 .eq. -1 )then
          write(8,*)'failed matching 2nd shock point of the shock n.'
     +    ,ish1
          stop
         endif

!        check if the shock crosses a bndry point
         if( s1 .lt. 0.d0 .or. s1 .gt. 1.d0.or.
     +       s2 .lt. 0.d0 .or. s2 .gt. 1.d0 )then
          write(8,*)'s(',i,') out of bounds',s1,s2
          stop
        endif
        if( iedg2 .ne. iedg1 )then
          write(8,*)'shock points (1) (2) not on the same bndry edge'
          write(8,*)iedg1,iedg2
          write(8,*)s1,s2
          stop
        endif

!        split the existing edges of the background mesh
!        if the first point of the shock is at the boundary
         write(8,*)'**********************'
         write(8,*)'shock:',ish1
         write(8,*)iedg1
         write(8,*)'**********************'
         if(iedg1.gt.0)then

            i   = iedg1
            i1  = ibndfac(1,i)
            i2  = ibndfac(2,i)
            ibc = ibndfac(3,i)

! ************************************************************

            if(nodcod(i1).lt.0.or.nodcod(i2).lt.0)then
!             write(*,*) 'stopped!'
!             write(*,*)'nodcod(1):',nodcod(i1)
!             write(*,*)'nodcod(2):',nodcod(i2)
!             write(*,*)'s1:',s1
!             write(*,*)'s2:',s2
!             write(*,*)
            if(nodcod(i1).lt.0.d+0.and.s1.lt.s2)then
              ibndfac(1,i)=ishplistu(ip1,ish1)
              idum1=i1
              idum2=ishplistd(ip1,ish1)
            elseif(nodcod(i1).lt.0.d+0.and.s1.gt.s2)then
              ibndfac(1,i)=ishplistd(ip1,ish1)
              idum1=i1
              idum2=ishplistu(ip1,ish1)
            elseif(nodcod(i2).lt.0.d+0.and.s1.lt.s2)then
              ibndfac(2,i)=ishplistd(ip1,ish1)
              idum1=i2
              idum2=ishplistu(ip1,ish1)
            elseif(nodcod(i2).lt.0.d+0.and.s1.gt.s2)then
              ibndfac(2,i)=ishplistu(ip1,ish1)
              idum1=i2
              idum2=ishplistd(ip1,ish1)
             endif
             do ibf=1,nbfac
               if(ibndfac(1,ibf).eq.idum1)ibndfac(1,ibf)=idum2
               if(ibndfac(2,ibf).eq.idum1)ibndfac(2,ibf)=idum2
              enddo

            else

! *****************************************************

            ibndfac(3,i)=-ibc
            write(8,*)'removing background edge ' ,i,i1,i2

!          create 2 new edges at  boundary
           ibfac=ibfac+1
           ibndfac(3,ibfac)=ibc
           if( s1 .lt. s2 )then
             ibndfac(1,ibfac)=i1
             ibndfac(2,ibfac)=ishplistd(ip1,ish1)
           else
            ibndfac(1,ibfac)=i1
            ibndfac(2,ibfac)=ishplistu(ip1,ish1)
           endif

           ibfac=ibfac+1
           ibndfac(3,ibfac)=ibc
           if( s1 .lt. s2 )then
             ibndfac(1,ibfac)=ishplistu(ip1,ish1)
             ibndfac(2,ibfac)=i2
           else
            ibndfac(1,ibfac)=ishplistd(ip1,ish1)
            ibndfac(2,ibfac)=i2
           endif
          endif
         endif

! point 2

          ish1 = shinspps(1,2,isppnts)
          i    = shinspps(2,2,isppnts)-1
          ip1  = 1+i*(nshockpoints(ish1)-1)

          x0 = xyshd(1,ip1,ish1)
          y0 = xyshd(2,ip1,ish1)

caldo
          do ibf=1,nbfac
           i1 = ibndfac(1,ibf)
           i2 = ibndfac(2,ibf)
           x1 = xy(1,i1)
           y1 = xy(2,i1)
           x2 = xy(1,i2)
           y2 = xy(2,i2)
          write(8,*)ibf,(y0-y1)*(x2-x1)-(x0-x1)*(y2-y1)
          enddo
caldo

          iedg1 = findbedg(xy,ndim,ibndfac,nbfac,x0,y0,s1)
          write(8,*)'typespecpoints:',typespecpoints(isppnts)
          write(8,*)'s(1) ' ,s1,x0,y0,iedg1
          if( iedg1 .eq. -1 )then
           write(8,*)'failed matching 1st shock point of the shock n.'
     +                ,ish1
           stop
          else
          write(8,*)'shockpoint (1) ' ,x0,y0,' falls within ',
     &   (ibndfac(i,iedg1),i=1,2)
          endif
          x0 = xyshu(1,ip1,ish1)
          y0 = xyshu(2,ip1,ish1)

          do ibf=1,nbfac
           i1 = ibndfac(1,ibf)
           i2 = ibndfac(2,ibf)
           x1 = xy(1,i1)
           y1 = xy(2,i1)
           x2 = xy(1,i2)
           y2 = xy(2,i2)
          write(8,*)ibf,(y0-y1)*(x2-x1)-(x0-x1)*(y2-y1)
          enddo

          iedg2 = findbedg(xy,ndim,ibndfac,nbfac,x0,y0,s2)

          write(8,*)'shockpoint (2) ' ,x0,y0,' falls within ',
     &   (ibndfac(i,iedg2),i=1,2)
         write(8,*)'s(2) ' ,s2,x0,y0,iedg2
         if( iedg2 .eq. -1 )then
          write(8,*)'failed matching 2nd shock point of the shock n.'
     +    ,ish1
          stop
         endif

!        check if the shock crosses a bndry point
         if( s1 .lt. 0.d0 .or. s1 .gt. 1.d0.or.
     +       s2 .lt. 0.d0 .or. s2 .gt. 1.d0 )then
          write(8,*)'s(',i,') out of bounds',s1,s2
          stop
        endif
        if( iedg2 .ne. iedg1 )then
          write(8,*)'shock points (1) (2) not on the same bndry edge'
          write(8,*)iedg1,iedg2
          write(8,*)s1,s2
          stop
        endif

!        split the existing edges of the background mesh
!        if the first point of the shock is at the boundary
         write(8,*)'**********************'
         write(8,*)'shock:',ish1
         write(8,*)iedg1
         write(8,*)'**********************'
         if(iedg1.gt.0)then

            i   = iedg1
            i1  = ibndfac(1,i)
            i2  = ibndfac(2,i)
            ibc = ibndfac(3,i)

! ***********************************************************

            if(nodcod(i1).lt.0.or.nodcod(i2).lt.0)then
!             write(*,*) 'stopped!'
!             write(*,*)'nodcod(1):',nodcod(i1)
!             write(*,*)'nodcod(2):',nodcod(i2)
!             write(*,*)'s1:',s1
!             write(*,*)'s2:',s2
!             write(*,*)
            if(nodcod(i1).lt.0.d+0.and.s1.lt.s2)then
              ibndfac(1,i)=ishplistu(ip1,ish1)
              idum1=i1
              idum2=ishplistd(ip1,ish1)
            elseif(nodcod(i1).lt.0.d+0.and.s1.gt.s2)then
              ibndfac(1,i)=ishplistd(ip1,ish1)
              idum1=i1
              idum2=ishplistu(ip1,ish1)
            elseif(nodcod(i2).lt.0.d+0.and.s1.lt.s2)then
              ibndfac(2,i)=ishplistd(ip1,ish1)
              idum1=i2
              idum2=ishplistu(ip1,ish1)
            elseif(nodcod(i2).lt.0.d+0.and.s1.gt.s2)then
              ibndfac(2,i)=ishplistu(ip1,ish1)
              idum1=i2
              idum2=ishplistd(ip1,ish1)
             endif
             do ibf=1,nbfac
               if(ibndfac(1,ibf).eq.idum1)ibndfac(1,ibf)=idum2
               if(ibndfac(2,ibf).eq.idum1)ibndfac(2,ibf)=idum2
              enddo

            else

! *****************************************************

            ibndfac(3,i)=-ibc
            write(8,*)'removing background edge ' ,i,i1,i2

!          create 2 new edges at boundary
           ibfac=ibfac+1
           ibndfac(3,ibfac)=ibc
           if( s1 .lt. s2 )then
             ibndfac(1,ibfac)=i1
             ibndfac(2,ibfac)=ishplistd(ip1,ish1)
           else
            ibndfac(1,ibfac)=i1
            ibndfac(2,ibfac)=ishplistu(ip1,ish1)
           endif

           ibfac=ibfac+1
           ibndfac(3,ibfac)=ibc
           if( s1 .lt. s2 )then
             ibndfac(1,ibfac)=ishplistu(ip1,ish1)
             ibndfac(2,ibfac)=i2
           else
            ibndfac(1,ibfac)=ishplistd(ip1,ish1)
            ibndfac(2,ibfac)=i2
           endif
          endif
         endif

         elseif(typespecpoints(isppnts).eq.'TE')then

          ish1 = shinspps(1,1,isppnts)
          i    = shinspps(2,1,isppnts)-1
          ip1  = 1+i*(nshockpoints(ish1)-1)

          ish2 = shinspps(1,3,isppnts)
          i    = shinspps(2,3,isppnts)-1
          ip2  = 1+i*(nshockpoints(ish2)-1)

          x0 = xysh(1,ip1,ish1)
          y0 = xysh(2,ip1,ish1)
!         write(*,*)'x0,y0 1 ',x0,y0

          iedg1 = findbedg(xy,ndim,ibndfac,nbfac,x0,y0,s1)
          write(8,*)'typespecpoints:',typespecpoints(isppnts)
          write(8,*)'s(1) ' ,s1,x0,y0,iedg1
          if( iedg1 .eq. -1)then
           write(8,*)'failed matching 1st shock point of the shock n.'
     +               ,ish1
           stop
          else
           write(8,*)'shockpoint (1)' ,x0,y0,' falls within ',
     &        (ibndfac(i,iedg1),i=1,2)
          endif

          i   = iedg1
          i1  = ibndfac(1,i)
          i2  = ibndfac(2,i)
          ibc = ibndfac(3,i)

          ibndfac(3,i)=-ibc
           write(8,*)'removing background edge ' ,i,i1,i2

!           create a new edge at boundary
            ibfac=ibfac+1
            ibndfac(3,ibfac)=ibc
            if( nodcod(i1)  .lt. 0.0d0 )then
             ibndfac(1,ibfac)=ishplistu(ip1,ish1)
             ibndfac(2,ibfac)=i2
             j1=1
            elseif(nodcod(i2)  .lt. 0.0d0)then
             ibndfac(1,ibfac)=i1
             ibndfac(2,ibfac)=ishplistu(ip1,ish1)
             j1=2
            else
             write(*,*)'condition not considered'
             stop
            endif

          x0 = xysh(1,ip2,ish2)
          y0 = xysh(2,ip2,ish2)
!         write(*,*)'x0,y0 2 ',x0,y0

          iedg2 = findbedg(xy,ndim,ibndfac,nbfac,x0,y0,s1)
          write(8,*)'typespecpoints:',typespecpoints(isppnts)
          write(8,*)'s(2) ' ,s1,x0,y0,iedg2
          if( iedg2 .eq. -1)then
           write(8,*)'failed matching 1st shock point of the shock n.'
     +               ,ish2
           stop
          else
           write(8,*)'shockpoint (2)' ,x0,y0,' falls within ',
     &        (ibndfac(i,iedg2),i=1,2)
          endif

          i   = iedg2
          i1  = ibndfac(1,i)
          i2  = ibndfac(2,i)
          ibc = ibndfac(3,i)
!          write(*,*)'iedg2    :',iedg2
!          write(*,*)'nodcod(1):',nodcod(i1),i1
!          write(*,*)'nodcod(2):',nodcod(i2),i2
!          write(*,*)'ibc:      ',ibc

          ibndfac(3,i)=-ibc
           write(8,*)'removing background edge ' ,i,i1,i2

!           create a new edge at boundary
            ibfac=ibfac+1
            ibndfac(3,ibfac)=ibc
            if( nodcod(i1)  .lt. 0.0d0 )then
             ibndfac(1,ibfac)=ishplistu(ip2,ish2)
             ibndfac(2,ibfac)=i2
             j2=1
            elseif(nodcod(i2)  .lt. 0.0d0)then
             ibndfac(1,ibfac)=i1
             ibndfac(2,ibfac)=ishplistu(ip2,ish2)
             j2=2
            else
             write(*,*)'condition not considered'
             stop
            endif

!          check cross
           x1= xy(1,ibndfac(1,ibfac))
           y1= xy(2,ibndfac(1,ibfac))
           x2= xy(1,ibndfac(2,ibfac))
           y2= xy(2,ibndfac(2,ibfac))
           x3= xy(1,ibndfac(1,ibfac-1))
           y3= xy(2,ibndfac(1,ibfac-1))
           x4= xy(1,ibndfac(2,ibfac-1))
           y4= xy(2,ibndfac(2,ibfac-1))

           idum1=ishel1(x1, y1, x1, y1 ,x2, y2, x3, y3,
     +                       x4, y4)

!          write(*,*)'ishel1:',idum1
           if(idum1.eq.0.)then

            if( j1.eq.1 )then
             ibndfac(1,ibfac-1)=ishplistu(ip2,ish2)
            elseif(j1.eq.2)then
             ibndfac(2,ibfac-1)=ishplistu(ip2,ish2)
            else
             write(*,*)'condition not considered'
             stop
            endif

            if( j2.eq.1 )then
             ibndfac(1,ibfac)=ishplistu(ip1,ish1)
            elseif(j2.eq.2)then
             ibndfac(2,ibfac)=ishplistu(ip1,ish1)
            else
             write(*,*)'condition not considered'
             stop
            endif

           endif

!          i   = iedg2
!          i1  = ibndfac(1,i)
!          i2  = ibndfac(2,i)
!          ibc = ibndfac(3,i)
!          write(*,*)'iedg2    :',iedg2
!          write(*,*)'nodcod(1):',nodcod(i1),i1
!          write(*,*)'nodcod(2):',nodcod(i2),i2
!          write(*,*)'ibc:      ',ibc

c         pause

         elseif(typespecpoints(isppnts).eq.'RRX'.or.
     +          typespecpoints(isppnts).eq.'RR')then

          ish1 = shinspps(1,1,isppnts)
          i    = shinspps(2,1,isppnts)-1
          ip1  = 1+i*(nshockpoints(ish1)-1)

          ish2 = shinspps(1,2,isppnts)
          i    = shinspps(2,2,isppnts)-1
          ip2  = 1+i*(nshockpoints(ish2)-1)

          x0 = xyshd(1,ip2,ish2)
          y0 = xyshd(2,ip2,ish2)

          write(8,*)'typespecpoints:',typespecpoints(isppnts)
          iedg1 = findbedg(xy,ndim,ibndfac,nbfac,x0,y0,s1)
          write(8,*)'s(1) ' ,s1,x0,y0,iedg1
          if( iedg1 .eq. -1 )then
           write(8,*)'failed matching 1st shock point of the shock n.'
     +                ,ish2
           stop
          else
          write(8,*)'shockpoint (1) ' ,x0,y0,' falls within ',
     &   (ibndfac(i,iedg1),i=1,2)
          endif
          x0 = xyshu(1,ip1,ish1)
          y0 = xyshu(2,ip1,ish1)
          iedg2 = findbedg(xy,ndim,ibndfac,nbfac,x0,y0,s2)

          write(8,*)'shockpoint (2) ' ,x0,y0,' falls within ',
     &   (ibndfac(i,iedg2),i=1,2)
         write(8,*)'s(2) ' ,s2,x0,y0,iedg2
         if( iedg2 .eq. -1 )then
          write(8,*)'failed matching 2nd shock point of the shock n.'
     +    ,ish1
          stop
         endif

!     check if the shock crosses a bndry point

         if( s1 .lt. 0.d0 .or. s1 .gt. 1.d0. or.
     +       s2 .lt. 0.d0 .or. s2 .gt. 1.d0 )then
          write(8,*)'out of bounds',s1,s2
          stop
        endif
        if( iedg1 .ne. iedg2 )then
          write(8,*)'shock points (1) (2) not on the same bndry edge'
          write(8,*)iedg1,iedg2
          write(8,*)s1,s2
          stop
        endif

!        split the existing edges of the background mesh
!        if the first point of the shock is at the boundary
         write(8,*)'**********************'
         write(8,*)'shock:',ish1
         write(8,*)iedg1
         write(8,*)'**********************'
         if(iedg1.gt.0)then
            i   = iedg1
            i1  = ibndfac(1,i)
            i2  = ibndfac(2,i)
            ibc = ibndfac(3,i)

! ***********************************************************

            if(nodcod(i1).lt.0.or.nodcod(i2).lt.0)then
!             write(*,*) 'stopped!'
!             write(*,*)'nodcod(1):',nodcod(i1)
!             write(*,*)'nodcod(2):',nodcod(i2)
!             write(*,*)'s1:',s1
!             write(*,*)'s2:',s2
!             write(*,*)
            if(nodcod(i1).lt.0.d+0.and.s1.lt.s2)then
              ibndfac(1,i)=ishplistu(ip1,ish1)
              idum1=i1
              idum2=ishplistd(ip2,ish2)
            elseif(nodcod(i1).lt.0.d+0.and.s1.gt.s2)then
              ibndfac(1,i)=ishplistd(ip2,ish2)
              idum1=i1
              idum2=ishplistu(ip1,ish1)
            elseif(nodcod(i2).lt.0.d+0.and.s1.lt.s2)then
              ibndfac(2,i)=ishplistd(ip2,ish2)
              idum1=i2
              idum2=ishplistu(ip1,ish1)
            elseif(nodcod(i2).lt.0.d+0.and.s1.gt.s2)then
              ibndfac(2,i)=ishplistu(ip1,ish1)
              idum1=i2
              idum2=ishplistd(ip2,ish2)
             endif
             do ibf=1,nbfac
               if(ibndfac(1,ibf).eq.idum1)ibndfac(1,ibf)=idum2
               if(ibndfac(2,ibf).eq.idum1)ibndfac(2,ibf)=idum2
              enddo

            else

! *****************************************************

            ibndfac(3,i)=-ibc
            write(8,*)'removing background edge ' ,i,i1,i2

!          create 2 new edges at boundary
           ibfac=ibfac+1
           ibndfac(3,ibfac)=ibc
           if( s1 .lt. s2 )then
             ibndfac(1,ibfac)=i1
             ibndfac(2,ibfac)=ishplistd(ip2,ish2)
           else
            ibndfac(1,ibfac)=i1
            ibndfac(2,ibfac)=ishplistu(ip1,ish1)
           endif

           ibfac=ibfac+1
           ibndfac(3,ibfac)=ibc
           if( s1 .lt. s2 )then
             ibndfac(1,ibfac)=ishplistu(ip1,ish1)
             ibndfac(2,ibfac)=i2
           else
            ibndfac(1,ibfac)=ishplistd(ip2,ish2)
            ibndfac(2,ibfac)=i2
           endif
          endif
         endif

! nothing to do
         elseif(typespecpoints(isppnts).eq.'TP')then

         elseif(typespecpoints(isppnts).eq.'QP')then

         elseif(typespecpoints(isppnts).eq.'EP')then

         elseif(typespecpoints(isppnts).eq.'C')then

         elseif(typespecpoints(isppnts).eq.'SP')then

         else

         write(*,*)'condition not implemented!'
         write(*,*)
         stop
         endif

        enddo

!     split the existing edges of the background mesh
!     if the last point of the shock is at the boundary
!
!        if(iedg(ish,3).gt.0)then
!          i   = iedg(ish,3)
!          i1  = ibndfac(1,i)
!          i2  = ibndfac(2,i)
!          ibc = ibndfac(3,i)
!          ibndfac(3,i)=-ibc
!          write(8,*)'removing background edge ' ,i,i1,i2
!
!     create 2 new edges at  boundary
!
!          ibfac=ibfac+1
!          ibndfac(3,ibfac)=ibc
!          if( s(ish,3) .lt. s(ish,4) )then
!           ibndfac(1,ibfac)=i1
!           ibndfac(2,ibfac)=ishplistd(nshockpoints(ish),ish)
!          else
!           ibndfac(1,ibfac)=i1
!           ibndfac(2,ibfac)=ishplistu(nshockpoints(ish),ish)
!          endif
!
!          ibfac=ibfac+1
!          ibndfac(3,ibfac)=ibc
!          if( s(ish,3) .lt. s(ish,4) )then
!           ibndfac(1,ibfac)=ishplistu(nshockpoints(ish),ish)
!           ibndfac(2,ibfac)=i2
!          else
!           ibndfac(1,ibfac)=ishplistd(nshockpoints(ish),ish)
!           ibndfac(2,ibfac)=i2
!          endif
!
!       endif
!     end do

      nbfac_sh=ibfac

!     create 4 new edges near the triple point

      close(8)

      close(130)

      return
      end


      integer function findbedg(xy,ndim,ibndfac,nbfac,xsh,ysh,s)

!     finds the bndry edge (of the background mesh)
!     the shock point (xsh,ysh) belongs to

      implicit none
      integer ndim,nbfac
      double precision xy(ndim,*)
      integer ibndfac(3,*)
      double precision xsh,ysh,s
      double precision x,y,x1,x2,y1,y2
      integer i1,i2,ibc,ibfac
      double precision toler
      parameter (toler=0.2d-6)
      double precision tline
      tline(x,y) = (y-y1)*(x2-x1)-(x-x1)*(y2-y1)

      do 10 ibfac = 1, nbfac
          i1 = ibndfac(1,ibfac)
          i2 = ibndfac(2,ibfac)
          ibc = ibndfac(3,ibfac)

          if( ibc.lt.0 )goto 10
          x1 = xy(1,i1)
          y1 = xy(2,i1)
          x2 = xy(1,i2)
          y2 = xy(2,i2)

          if( abs(tline(xsh,ysh)).le.toler )then

              if( abs(y2-y1) .gt. abs(x2-x1) )then
                  s = (ysh-y1)/(y2-y1)
              else
                  s = (xsh-x1)/(x2-x1)

              endif
!             if( 0.d0 .le. s .and. s .le. 1.d0 )then
              if( 0.d0 .le. s .and. s .le. 1.d0 )then
                  findbedg = ibfac
                  return
              endif
          endif
   10 continue
      findbedg = -1

      return
      end
