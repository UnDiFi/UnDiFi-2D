! Read the file containing the information concerning shocks and discontinuities

      subroutine re_sdw_info(xysh,
     +                       zroeshu,
     +                       zroeshd,
     +                       zroeshuold,
     +                       zroeshdold,
     +                       nodcodsh,
     +                       coor,      !vale
     +                       ibndfac,   !vale
     +                       nbfac,     !vale
     +                       npoin,     !vale
     +                       nshocks,
     +                       nshockpoints,
     +                       nshockedges,
     +                       typesh,
     +                       nspecpoints,
     +                       typespecpoints,
     +                       shinspps,
     +                       ispclr)

      implicit none
      include 'paramt.h'
      include 'shock.com'

!     .. scalar arguments ..
      integer nshocks,nspecpoints,nshockedges(*),nshockpoints(*),
     +        isppnts,idummy,nshe
      character*1 typesh(*)
      character*5 typespecpoints(*)

!     .. array arguments ..
      double precision xysh(ndim,npshmax,*),
     &                 zroeshu(ndof,npshmax,*),
     &                 zroeshd(ndof,npshmax,*),
     &                 zroeshuold(ndof,npshmax,*),
     &                 zroeshdold(ndof,npshmax,*)
! vale
      integer inode, ibfac
      integer npoin, nbfac, ibndfac(3,*)
      double precision coor(ndim,*)
      double precision xywedge(2)
      logical foundbgnwedge
! vale
      integer nodcodsh(npshmax,*),
     +        shinspps(2,5,*),
     +        ispclr(5,*)

!     .. local scalars ..
      integer i,k,ish

!     open log file
      open(8,file='log/re_sdw_info.log')

!     initialize nodcodsh which is part of nodcod
!     if the code -99 is used this means no shock point
!     if the code 10 is used this means shock point
      do ish = 1,nshmax
        do  k = 1,npshmax
          nodcodsh(k,ish)=-99
        enddo
      enddo

!     open file outnew containing the information concerning shocks and discontinuities
      open(12,file='sh00.dat',status='old')
      write(8,*)' open file sh00.dat'
      read(12,*)nshocks
      write(8,*)' found n.',nshocks,'shocks/discontinuities'
      do ish=1,nshocks
        write(8,*)' shock/discontinuity n.',ish
        read(12,*)nshockpoints(ish),typesh(ish)
        write(8,*)' kind of discontinuity:',typesh(ish)
        write(8,*)' n. of points',nshockpoints(ish)
        write(8,*)'shock-point coordinates, downstream, upstream states'

        nshockedges(ish) = nshockpoints(ish)-1
        do k = 1, nshockpoints(ish)
          read(12,*)(xysh(i,k,ish),i=1,ndim),
     &              (zroeshd(i,k,ish),i=1,ndof),
     &              (zroeshu(i,k,ish),i=1,ndof)
          write(8,*)(xysh(i,k,ish),i=1,ndim),
     &              (zroeshd(i,k,ish),i=1,ndof),
     &              (zroeshu(i,k,ish),i=1,ndof)

          nodcodsh(k,ish)=10
        end do

        do k = 1,nshockpoints(ish)
          do i = 1,ndof
            zroeshdold(i,k,ish)=zroeshd(i,k,ish)
            zroeshuold(i,k,ish)=zroeshu(i,k,ish)
          enddo
        enddo

!       write(8,*)' in shock n.:', ish
!       write(8,*)' there are ', nshockpoints(ish),' shock vertices'
!       write(8,*)' there are ', nshockedges(ish),' shock edges'

      end do

      idummy=0

      read(12,*)nspecpoints
      write(8,*)nspecpoints
      do isppnts=1,nspecpoints
        read(12,*)typespecpoints(isppnts)
        write(8,*)typespecpoints(isppnts)
        if(typespecpoints(isppnts).eq.'TP')then                 ! internal special point: triple point
          nshe=4
          idummy=idummy+nshe

          do k = 1,nshe
            read(12,*)shinspps(1,k,isppnts),shinspps(2,k,isppnts)
            write(8,*)shinspps(1,k,isppnts),shinspps(2,k,isppnts)
          enddo

        elseif(typespecpoints(isppnts).eq.'QP')then             ! internal special point: quad point
          nshe=5
          idummy=idummy+nshe

          do k = 1,nshe
            read(12,*)shinspps(1,k,isppnts),shinspps(2,k,isppnts)
            write(8,*)shinspps(1,k,isppnts),shinspps(2,k,isppnts)
          enddo

        elseif(typespecpoints(isppnts).eq.'TE')then             ! trailing edge point:
          nshe=3
          idummy=idummy+nshe

          do k = 1,nshe
            read(12,*)shinspps(1,k,isppnts),shinspps(2,k,isppnts)
            write(8,*)shinspps(1,k,isppnts),shinspps(2,k,isppnts)
          enddo

        elseif(typespecpoints(isppnts).eq.'RRX')then            ! boundary special point: regular reflection along x-wall
          nshe=2                                                ! obsolte - use rr
          idummy=idummy+nshe

          do k = 1,nshe
            read(12,*)shinspps(1,k,isppnts),shinspps(2,k,isppnts)
            write(8,*)shinspps(1,k,isppnts),shinspps(2,k,isppnts)
          enddo

        elseif(typespecpoints(isppnts).eq.'RR')then             ! boundary special point: regular reflection along
                                                                ! a generic (also curved) wall
          nshe=2
          idummy=idummy+nshe

          do k = 1,nshe
            read(12,*)shinspps(1,k,isppnts),shinspps(2,k,isppnts),
     +                ispclr(k,isppnts)
            write(8,*)shinspps(1,k,isppnts),shinspps(2,k,isppnts),
     +                ispclr(k,isppnts)
          enddo

        elseif(typespecpoints(isppnts).eq.'WPNRX')then          ! boundary special point: wall point without reflection
          nshe=1                                                ! floating along x direction
          idummy=idummy+nshe                                    ! obsolete - use fwp

          do k = 1,nshe
            read(12,*)shinspps(1,k,isppnts),shinspps(2,k,isppnts)
            write(8,*)shinspps(1,k,isppnts),shinspps(2,k,isppnts)
          enddo

        elseif(typespecpoints(isppnts).eq.'WPNRY')then          ! boundary special point: wall point without reflection
          nshe=1                                                ! floating along y direction
          idummy=idummy+nshe                                    ! obsolete - use fwp

          do k = 1,nshe
            read(12,*)shinspps(1,k,isppnts),shinspps(2,k,isppnts)
            write(8,*)shinspps(1,k,isppnts),shinspps(2,k,isppnts)
          enddo

        elseif(typespecpoints(isppnts).eq.'IPX')then            ! boundary special point: inlet point
          nshe=1                                                ! floating along x direction
          idummy=idummy+nshe

          do k = 1,nshe
            read(12,*)shinspps(1,k,isppnts),shinspps(2,k,isppnts)
            write(8,*)shinspps(1,k,isppnts),shinspps(2,k,isppnts)
          enddo

        elseif(typespecpoints(isppnts).eq.'IPY')then            ! boundary special point: inlet point
          nshe=1                                                ! floating along y direction
          idummy=idummy+nshe

          do k = 1,nshe
            read(12,*)shinspps(1,k,isppnts),shinspps(2,k,isppnts)
            write(8,*)shinspps(1,k,isppnts),shinspps(2,k,isppnts)
          enddo

        elseif(typespecpoints(isppnts).eq.'OPX')then            ! boundary special point: outlet point
          nshe=1                                                ! floating along x direction
          idummy=idummy+nshe

          do k = 1,nshe
            read(12,*)shinspps(1,k,isppnts),shinspps(2,k,isppnts)
            write(8,*)shinspps(1,k,isppnts),shinspps(2,k,isppnts)
          enddo

        elseif(typespecpoints(isppnts).eq.'OPY')then            ! boundary special point: outlet point
          nshe=1                                                ! floating along y direction
          idummy=idummy+nshe

          do k = 1,nshe
            read(12,*)shinspps(1,k,isppnts),shinspps(2,k,isppnts)
            write(8,*)shinspps(1,k,isppnts),shinspps(2,k,isppnts)
          enddo

        elseif(typespecpoints(isppnts).eq.'EP')then             ! boundary special point: end point
          nshe=1
          idummy=idummy+nshe

          do k = 1,nshe
            read(12,*)shinspps(1,k,isppnts),shinspps(2,k,isppnts)
            write(8,*)shinspps(1,k,isppnts),shinspps(2,k,isppnts)
          enddo

        elseif(typespecpoints(isppnts).eq.'C')then              ! connection between two shocks
          nshe=2
          idummy=idummy+nshe

          do k = 1,nshe
            read(12,*)shinspps(1,k,isppnts),shinspps(2,k,isppnts)
            write(8,*)shinspps(1,k,isppnts),shinspps(2,k,isppnts)
          enddo

       elseif(typespecpoints(isppnts).eq.'FWP')then             ! boundary special point: floating wall point along a coloured wall
          nshe=1
          idummy=idummy+nshe

          do k = 1,nshe
            read(12,*)shinspps(1,k,isppnts),shinspps(2,k,isppnts),
     +                ispclr(k,isppnts)
            write(8,*)shinspps(1,k,isppnts),shinspps(2,k,isppnts),
     +                ispclr(k,isppnts)
          enddo

        elseif(typespecpoints(isppnts).eq.'PC')then             ! periodic connection between two shocks
          nshe=2
          idummy=idummy+nshe

          do k = 1,nshe
            read(12,*)shinspps(1,k,isppnts),shinspps(2,k,isppnts),
     +                ispclr(k,isppnts)
            write(8,*)shinspps(1,k,isppnts),shinspps(2,k,isppnts),
     +                ispclr(k,isppnts)
          enddo

        elseif(typespecpoints(isppnts).eq.'SP')then              ! boundary special point: start point
                                                                 ! (shock point originated from a characteristic coalescence)
          nshe=1
          idummy=idummy+nshe

          do k = 1,nshe
            read(12,*)shinspps(1,k,isppnts),shinspps(2,k,isppnts)
            write(8,*)shinspps(1,k,isppnts),shinspps(2,k,isppnts)
          enddo

        else
          write(8,*)'condition not implemented'
          write(*,*)'condition not implemented'
          stop
        endif
      enddo

!     check condition on special points
      if(idummy.ne.2*nshocks)then
        write(8,*)'wrong n. of conditions on special points '
        write(8,*)'n. imposed condition on special points:',nspecpoints
        write(8,*)'n. required conditions:',2*nshocks
        write(8,*)'n. imposed conditions:',idummy
        write(*,*)'wrong n. of conditions on special points '
        stop
      endif

      close(8)

      return
      end
