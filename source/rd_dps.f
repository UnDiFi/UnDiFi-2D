! Remove discontinuity points where they are too close and insert new
! discontinuity points where they are too far.

      subroutine rd_dps(
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
      double precision dum,length_rel_min,length_rel_max

!     .. array arguments ..
!     character*(*) fname
      double precision sh_edge_lgth(npshmax)

!     .. local scalars ..
      double precision dt
      integer i,im,iv,ish,k,ile_min,ile_max,np_c,np_i

!     open log file
      open(8,file='log/rd_dps.log')

      do ish = 1,nshocks
        ile_min=0
        length_rel_min=1.0
        ile_max=0
        length_rel_max=1.0

!       compute length of shock edge
        do iv = 1, nshockedges(ish)
          sh_edge_lgth(iv)=0.0d0
          do k=1,ndim
            sh_edge_lgth(iv)=sh_edge_lgth(iv)+
     .                       (xysh(k,iv,ish)-xysh(k,iv+1,ish))**2
          enddo
          sh_edge_lgth(iv)=sqrt(sh_edge_lgth(iv))
          dum=sh_edge_lgth(iv)/dxcell
          if(dum.lt.0.5)then
            if(length_rel_min.gt.dum)then
              length_rel_min=dum
              ile_min=iv
            endif
          endif
          if(dum.gt.1.5)then
            if(length_rel_max.lt.dum)then
              length_rel_max=dum
              ile_max=iv
            endif
          endif
        enddo

        if(nshockpoints(ish).le.3)then
          ile_min=0.0d0
          ile_max=0.0d0
        endif

        if(ile_min.ne.0)then
          write(8,*)'before'
          do iv = 1,nshockpoints(ish)
            write(8,*)iv,zroeshd(1,iv,ish),zroeshd(2,iv,ish)
          enddo

          np_c=ile_min
          if(sh_edge_lgth(ile_min-1).gt.sh_edge_lgth(ile_min+1))
     +      np_c=ile_min+1
          if(ile_min.eq.1)np_c=2
          if(ile_min.eq.nshockedges(ish))np_c=ile_min

          do iv=np_c+1,nshockpoints(ish)
            do k=1,ndim
              xysh(k,iv-1,ish)=xysh(k,iv,ish)
            enddo
          do k=1,ndof
            zroeshu(k,iv-1,ish)=zroeshu(k,iv,ish)
            zroeshd(k,iv-1,ish)=zroeshd(k,iv,ish)
          enddo
         enddo
         do k=1,ndim
          xysh(k,iv,ish)=0.0d+0
         enddo
         do k=1,ndof
          zroeshu(k,iv,ish)=0.0d+0
          zroeshd(k,iv,ish)=0.0d+0
         enddo

         nshockpoints(ish)=nshockpoints(ish)-1
         nshockedges(ish) =nshockedges (ish)-1

         write(8,*)'after'
         do iv=1,nshockpoints(ish)
          write(8,*)iv,zroeshd(1,iv,ish),zroeshd(2,iv,ish)
         enddo

        endif

        if(ile_max.ne.0)then
         write(8,*)'before'
         do iv=1,nshockpoints(ish)
          write(8,*)iv,zroeshd(1,iv,ish),zroeshd(2,iv,ish)
         enddo

        np_i=ile_max

         do iv=nshockpoints(ish),np_i+1,-1
          do k=1,ndim
           xysh(k,iv+1,ish)=xysh(k,iv,ish)
          enddo
          do k=1,ndof
           zroeshu(k,iv+1,ish)=zroeshu(k,iv,ish)
           zroeshd(k,iv+1,ish)=zroeshd(k,iv,ish)
          enddo
         enddo
         do k=1,ndim
          xysh(k,np_i+1,ish)=0.5*(xysh(k,np_i,ish)+xysh(k,np_i+2,ish))
         enddo
         do k=1,ndof
          zroeshu(k,np_i+1,ish)=0.5*(zroeshu(k,np_i,ish)+
     .                               zroeshu(k,np_i+2,ish))
          zroeshd(k,np_i+1,ish)=0.5*(zroeshd(k,np_i,ish)+
     .                               zroeshd(k,np_i+2,ish))
         enddo

         nshockpoints(ish)=nshockpoints(ish)+1
         nshockedges(ish) =nshockedges (ish)+1

         write(8,*)'after'
         do iv=1,nshockpoints(ish)
          write(8,*)iv,zroeshd(1,iv,ish),zroeshd(2,iv,ish)
         enddo

        endif

        enddo

      close(8)
      return
      end
