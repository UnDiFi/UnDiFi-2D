! Re-distribute the shock points in each shock with uniform spacing

      subroutine rd_dps_eq(
     +           xysh,
     +           zroeshu,
     +           zroeshd,
     +           nshocks,
     +           nshockpoints,
     +           nshockedges)

      implicit none
      include 'paramt.h'
      include 'shock.com'

!     .. scalar arguments ..
      integer iter,nshocks,nshockpoints(nshmax),nshockedges(nshmax)

!     .. array arguments ..
      double precision
     +                 xysh   (ndim, npshmax, *),
     +                 zroeshu(ndof, npshmax, *),
     +                 zroeshd(ndof, npshmax, *)

!     .. local scalar  ..
      double precision dum,dx,sh_edge_lgth,alpha,beta

!     .. array arguments ..
!     character*(*) fname
      double precision xysh_new(ndim,npshmax),
     +                 zroeshu_new(ndof, npshmax),
     +                 zroeshd_new(ndof, npshmax),
     +                 sh_absc(npshmax), sh_absc_new(npshmax)

!     .. local scalars ..
      double precision ds,dsj
      integer i,im,iv,ish,k,j,nshockpoints_new

!     open log file
      open(8,file='log/rd_dps_eq.log')

      do ish = 1,nshocks

        write(8,*)'shock n.:', ish

!       compute length of shock edge
        sh_absc(1)=0.0
        do iv = 1, nshockedges(ish)
          i=iv+1
          sh_edge_lgth=0.0d0
          do k=1,ndim
            sh_edge_lgth=sh_edge_lgth+
     .                   (xysh(k,iv,ish)-xysh(k,iv+1,ish))**2
          enddo
          sh_edge_lgth=sqrt(sh_edge_lgth)
          sh_absc(i)=sh_absc(i-1)+sh_edge_lgth
        enddo
        write(8,*)'# of points old distribution:',nshockpoints(ish)
!       enddo
!       write(*,*)

!       compute number of shock points of redistribution
        nshockpoints_new=sh_absc(nshockpoints(ish))/dxcell+1

!       check the new number of shock points
        if(nshockpoints_new.gt.npshmax)then
          write(6,*)'too many shock points! increase npshmax in
     +               paramt.f'
          call exit(1)
        endif

!       compute distribution step
        dx=sh_absc(nshockpoints(ish))/(nshockpoints_new-1)

        sh_absc_new(1)=0.0
        do iv = 1,nshockpoints_new
          i=iv+1
          sh_absc_new(i)=sh_absc_new(i-1)+dx
        enddo

        write(8,*)'number of points new distribution:',nshockpoints_new
!       do i = 1,nshockpoints_new
!         write(*,*)i,sh_absc_new(i)
!       enddo
!       write(*,*)

!       interpolation of new shock points and set the first and the last point
        do k = 1,ndim
          xysh_new(k,1)=xysh(k,1,ish)
        enddo
        do k = 1,ndof
          zroeshu_new(k, 1)=zroeshu(k, 1,ish)
          zroeshd_new(k, 1)=zroeshd(k, 1,ish)
        enddo

        do k = 1,ndim
          xysh_new(k,nshockpoints_new)=xysh(k,nshockpoints(ish),ish)
        enddo
        do k = 1,ndof
        zroeshu_new(k,nshockpoints_new)=zroeshu(k,nshockpoints(ish),ish)
        zroeshd_new(k,nshockpoints_new)=zroeshd(k,nshockpoints(ish),ish)
        enddo

!       computation of internal points
        do i=2,nshockpoints_new-1
          do j=2,nshockpoints(ish)
            if(sh_absc_new(i).ge.sh_absc(j-1).and.
     +        sh_absc_new(i).lt.sh_absc(j))then
              ds=sh_absc(j)-sh_absc(j-1)
              dsj=sh_absc_new(i)-sh_absc(j-1)
              alpha=dsj/ds
              beta=1.0d0-alpha
              do k=1,ndim
                xysh_new(k,i)=beta*xysh(k,j-1,ish)+alpha*xysh(k,j,ish)
              enddo
              do k=1,ndof
                zroeshu_new(k, i)=beta *zroeshu(k, j-1,ish)+
     +                            alpha*zroeshu(k, j  ,ish)
                zroeshd_new(k, i)=beta *zroeshd(k, j-1,ish)+
     +                            alpha*zroeshd(k, j  ,ish)
              enddo
            endif
          enddo
        enddo

!       re-write the state arrays and indices
        write(8,*)'new shock point coordinates'
        do i=1,nshockpoints_new
          do k=1,ndim
            xysh(k,i,ish)=xysh_new(k,i)
          enddo
          write(8,*)k,(xysh(k,i,ish),k=1,ndim)
          do k=1,ndof
            zroeshu(k, i,ish)=zroeshu_new(k, i)
            zroeshd(k, i,ish)=zroeshd_new(k, i)
          enddo

          nshockpoints(ish)=nshockpoints_new
          nshockedges(ish)=nshockpoints_new-1
        enddo
      enddo

      close(8)
      return
      end
