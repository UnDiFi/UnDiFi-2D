! Write the file containing the information concerning shocks discontinuities

      subroutine wrt_sdw_info(xysh,
     +                        zroeshu,
     +                        zroeshd,
     +                        nodcodsh,
     +                        nshocks,
     +                        nshockpoints,
     +                        nshockedges,
     +                        typesh,
     +                        nspecpoints,
     +                        typespecpoints,
     +                        shinspps,
     +                        ispclr,
     +                        ispclr1,
     +                        ispclr2)

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
     &                 zroeshd(ndof,npshmax,*)

      integer nodcodsh(npshmax,*),
     +        shinspps(2,5,*),
     +        ispclr(5,*),
     +        ispclr1(5,*),
     +        ispclr2(5,*)

!     .. local scalars ..
      integer i,k,ish

!     open log file
      open(8,file='log/wrt_sdw_info.log')

!     open file sh99 containing the informations concerning shocks and discontinuities
      open(12,file='sh99.dat')
      write(8,*)' open file sh99.dat'
      write(12,*)nshocks
      write(8,*)' n.',nshocks,'shocks/discontinuities'
      do ish=1,nshocks
        write(8,*)' shock/discontinuity n.',ish
        write(12,*)nshockpoints(ish),typesh(ish)
        write(8,*)' kind of discontinuity:',typesh(ish)
        write(8,*)' n. of points',nshockpoints(ish)

        do k = 1,nshockpoints(ish)
          write(12,*)(xysh(i,k,ish),i=1,ndim),
     &              (zroeshd(i,k,ish),i=1,ndof),
     &              (zroeshu(i,k,ish),i=1,ndof)

        end do

!       write(8,*)' in shock n.:', ish
!       write(8,*)' there are ', nshockpoints(ish),' shock vertices'
!       write(8,*)' there are ', nshockedges(ish),' shock edges'

      end do

      write(12,*)nspecpoints
      write(8,*) nspecpoints
      do isppnts = 1,nspecpoints
        write(12,*)typespecpoints(isppnts)
        write(8,*) typespecpoints(isppnts)
        if (typespecpoints(isppnts) .eq. 'TP') then         ! internal special point: triple point
          nshe = 4

          do k=1,nshe
            write(12,*)shinspps(1,k,isppnts),shinspps(2,k,isppnts)
            write(8,*) shinspps(1,k,isppnts),shinspps(2,k,isppnts)
          enddo

        elseif (typespecpoints(isppnts) .eq. 'QP') then     ! internal special point: quad point
          nshe = 5

          do k=1,nshe
           write(12,*)shinspps(1,k,isppnts),shinspps(2,k,isppnts)
           write(8,*)shinspps(1,k,isppnts),shinspps(2,k,isppnts)
          enddo

        elseif (typespecpoints(isppnts) .eq. 'TE') then     ! trailing edge point
          nshe = 3

          do k=1,nshe
           write(12,*)shinspps(1,k,isppnts),shinspps(2,k,isppnts)
           write(8,*)shinspps(1,k,isppnts),shinspps(2,k,isppnts)
          enddo

        elseif (typespecpoints(isppnts) .eq. 'RRX') then    ! boundary special point: regular reflection along
          nshe = 2

          do k=1,nshe
           write(12,*)shinspps(1,k,isppnts),shinspps(2,k,isppnts)
           write(8,*)shinspps(1,k,isppnts),shinspps(2,k,isppnts)
          enddo

        elseif (typespecpoints(isppnts) .eq. 'RR') then     ! boundary special point: regular reflection along a generic (also curved) wall
          nshe = 2
          idummy=idummy+nshe

          do k=1,nshe
           write(12,*)shinspps(1,k,isppnts),shinspps(2,k,isppnts),
     +               ispclr(k,isppnts)
           write(8,*)shinspps(1,k,isppnts),shinspps(2,k,isppnts),
     +               ispclr(k,isppnts)
          enddo

        elseif (typespecpoints(isppnts) .eq. 'WPNRX') then  ! boundary special point: wall point without reflection
          nshe = 1                                          ! floating along x direction

          do k=1,nshe                                   ! obsolete - use fwp
           write(12,*)shinspps(1,k,isppnts),shinspps(2,k,isppnts)
           write(8,*)shinspps(1,k,isppnts),shinspps(2,k,isppnts)
          enddo

        elseif (typespecpoints(isppnts) .eq. 'WPNRY') then  ! boundary special point: wall point without reflection
          nshe = 1                                          ! floating along y direction

          do k=1,nshe                                   ! obsolete - use fwp
           write(12,*)shinspps(1,k,isppnts),shinspps(2,k,isppnts)
           write(8,*)shinspps(1,k,isppnts),shinspps(2,k,isppnts)
          enddo

        elseif (typespecpoints(isppnts) .eq. 'FWP') then    ! boundary special point: floating wall point along a coloured wall
          nshe = 1

          do k=1,nshe
           write(12,*)shinspps(1,k,isppnts),shinspps(2,k,isppnts),
     +                ispclr(k,isppnts)
           write(8,*)shinspps(1,k,isppnts),shinspps(2,k,isppnts),
     +               ispclr(k,isppnts)
          enddo

        elseif (typespecpoints(isppnts) .eq. 'IPX') then    ! boundary special point: inlet point
          nshe = 1                                          ! floating along x direction

          do k=1,nshe
           write(12,*)shinspps(1,k,isppnts),shinspps(2,k,isppnts)
           write(8,*)shinspps(1,k,isppnts),shinspps(2,k,isppnts)
          enddo

        elseif (typespecpoints(isppnts) .eq. 'IPY') then    ! boundary special point: inlet point
          nshe = 1                                          ! floating along y direction

          do k=1,nshe
           write(12,*)shinspps(1,k,isppnts),shinspps(2,k,isppnts)
           write(8,*)shinspps(1,k,isppnts),shinspps(2,k,isppnts)
          enddo

        elseif (typespecpoints(isppnts) .eq. 'OPX') then    ! boundary special point: outlet point
          nshe = 1                                          ! floating along x direction

          do k=1,nshe
           write(12,*)shinspps(1,k,isppnts),shinspps(2,k,isppnts)
           write(8,*) shinspps(1,k,isppnts),shinspps(2,k,isppnts)
          enddo

        elseif (typespecpoints(isppnts) .eq. 'OPY') then    ! boundary special point: outlet point
          nshe = 1                                          ! floating along y direction

          do k=1,nshe
           write(12,*)shinspps(1,k,isppnts),shinspps(2,k,isppnts)
           write(8,*) shinspps(1,k,isppnts),shinspps(2,k,isppnts)
          enddo

        elseif (typespecpoints(isppnts) .eq. 'EP') then     ! boundary special point: end point
          nshe = 1

          do k=1,nshe
           write(12,*)shinspps(1,k,isppnts),shinspps(2,k,isppnts)
           write(8,*) shinspps(1,k,isppnts),shinspps(2,k,isppnts)
          enddo

!       watch out: conditions SPX not added because not tested yet.
        elseif (typespecpoints(isppnts) .eq. 'SP') then     ! boundary special point: sonic point inside the field
          nshe=1

          do k=1,nshe
           write(12,*)shinspps(1,k,isppnts),shinspps(2,k,isppnts)
           write(8,*)shinspps(1,k,isppnts),shinspps(2,k,isppnts)
          enddo

        elseif (typespecpoints(isppnts) .eq. 'C') then      !  connection point
          nshe=2

          do k=1,nshe
            write(12,*)shinspps(1,k,isppnts),shinspps(2,k,isppnts)
            write(8,*)shinspps(1,k,isppnts),shinspps(2,k,isppnts)
          enddo

       elseif (typespecpoints(isppnts) .eq. 'PC') then      ! periodic connection between two shocks
          nshe = 2

          do k = 1,nshe
            write(12,*)shinspps(1,k,isppnts),shinspps(2,k,isppnts),
     +                 ispclr(k,isppnts)
            write(8,*)shinspps(1,k,isppnts),shinspps(2,k,isppnts),
     +                ispclr(k,isppnts)
          enddo

        else
          write(8,*)'condition not implemented!'
          write(*,*)'condition not implemented!'
          stop
        endif
      enddo

      close(8)

      return
      end
