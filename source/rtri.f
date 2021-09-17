      subroutine rtri(iedgptr,
     +                nface,
     +                ibndfac,
     +                nbfac,
     +                icelnod,
     +                icelcel,
     +                nvt,
     +                nelem,
     +                coor,
     +                zroe,
     +                nodcode,
     +                npoin,
     +                ftype,
     +                imode,
     +                fname)

!      read data from:
!      a node  file if ftype = "node"
!      a poly  file if ftype = "poly"
!      a ele   file if ftype = "ele"
!      a neigh file if ftype = "neigh"
!      a edge  file if ftype = "edge"

!      only dimensions are read if imode == 0
!      the actual data are read if imode != 0

      implicit none
      include 'paramt.h'

!     .. scalar arguments ..
      integer nface,nelem,npoin,nvt,nbfac,imode

!     .. array arguments ..
      character*(*) fname,ftype
      character fwork*255
      double precision coor(ndim,npoin),zroe(ndof,npoin)
      integer ibndfac(3,nbfac),icelnod(nvt,nelem),nodcode(npoin),
     .icelcel(nvt,nelem),iedgptr(3,nface)

!     .. local scalars ..
      integer ia,k,ielem,iface,ipoin,idum,iattr,nhole,i,idum1,idum2

!     .. external functions ..
      integer  icycl,lenstr
      external icycl,lenstr

!     open log file
      open(88,file='log/rtri.log')

      idum = 0
      k = lenstr(fname)
      if (ftype(1:4).eq.'node') then
        fwork(1:k+5) = fname(1:k)//".node"
!       write(8,*)'file name:',fwork(1:k+5)
        open(unit=19,file=fwork(1:k+5),status='old')

!     node files (check triangle docs for explanations)
!
!     * first line: <# of vertices> <dimension (must be 2)> <# of
!     * attributes> <# of boundary markers (0 or 1)>
!     * remaining lines: <vertex #> <x> <y> [attributes]
!     * [boundary marker]
!
!     blank lines and comments prefixed by `#' may be placed
!     anywhere. vertices must be numbered consecutively, starting from one or zero.

        if(imode.eq.0)then
!         read(19,*)npoin,ndim,ndof,iattr
          read(19,*)npoin,idum1,idum2,iattr
          if(idum1.ne.ndim.or.idum2.ne.ndof)then
           write(88,*)' warning !!!!! '
           write(88,*)' warning !!!!! '
           write(88,*)' ndim:',idum1 ,' ndof:',idum2,' in file *.node'
           write(88,*)' ndim:',ndim ,' ndof:',ndof,' in param.h'
           write(88,*)' warning !!!!! '
           write(88,*)' warning !!!!! '
           stop
          endif

!         write(88,*)' there are ',npoin,' gridpoints in ',fwork(1:k+5)
!         write(88,*)' there are ',ndof,' degrees of freedom in ',
!    &    fwork(1:k+5)
          if(iattr.ne.1)
     .       stop '<# of attributes should be 1 in node file>'
          return
          if(ndof.ne.4)
     .       stop '<# of variables should be 4 in node file>'
          return
        else
          do 92 ipoin = 1,npoin
             read(19,fmt=*)i,(coor(ia,ipoin),ia=1,ndim),
     1       (zroe(ia,ipoin),ia=1,ndof), nodcode(ipoin)
!            write(88,*)i,(coor(ia,ipoin),ia=1,ndim)
             if(nodcode(ipoin).ne.0)idum = idum+1
   92     continue
          write(88,*)' there are ',idum,' bndry nodes'
          close(19)
        endif
      elseif(ftype(1:4).eq.'poly')then
        if(imode.eq.0)then
          fwork(1:k+5) = fname(1:k)//".poly"
          open(unit=19,file=fwork(1:k+5),status='old')

!     empty node list; this is just a poly file

          read(19,*)ia,idum,idum,idum
          if(ia.ne.0)
     &       stop 'node list should be empty in the poly file'
          read(19,*)nbfac,iattr

          if(iattr.ne.1)
     .       stop '<# of attributes should be 1 in poly file>'
          write(88,*)' there are ',nbfac,' polylines in ',fwork(1:k+5)
          return
        else

          do 90 i= 1,nbfac
            read(19,fmt=*)iface,(ibndfac(ia,i),ia=1,3)
   90     continue
!         read nof holes
          read(19,*)nhole
          close(19)
          return
        endif

!     read ele files
      elseif(ftype(1:3).eq.'ele')then
        if(imode.eq.0)then
          fwork(1:k+4) = fname(1:k)//".ele"
          open(unit=19,file=fwork(1:k+4),status='old')
          read(19,*)nelem,nvt,idum
          if(nvt.ne.3)stop '<nodes per triangle> must be 3 in .ele file'
          if(idum.ne.0)stop '<# of attributes> must be 0 in .ele file'
          write(88,*)' there are ',nelem,' triangles in ',fwork(1:k+4)
          return
        else
          do 94 ielem = 1,nelem
            read(19,fmt=*)idum,(icelnod(ia,ielem),ia=1,nvt)
   94     continue
          close(19)
          return
        endif

!     read neigh files
      elseif(ftype(1:5).eq.'neigh')then
        if(imode.eq.0)then
          fwork(1:k+6) = fname(1:k)//".neigh"
          open(unit=19,file=fwork(1:k+6),status='old')
          read(19,*)nelem,nvt
          if(nvt.ne.3)stop '<nodes per triangle> must be 3 in .ele file'
          write(88,*)' there are ',nelem,' triangles in ',fwork(1:k+6)
          return
        else
          do 96 ielem = 1,nelem
            read(19,fmt=*)idum,(icelcel(ia,ielem),ia=1,nvt)
   96     continue
          close(19)
          return
        endif

!     read edge files
      elseif(ftype(1:4).eq.'edge')then
        if(imode.eq.0)then
          fwork(1:k+5) = fname(1:k)//".edge"
          open(unit=19,file=fwork(1:k+5),status='old')
          read(19,*)nface,iattr
          write(88,*)' there are ',nface,' edges in ',fwork(1:k+5)
          if(iattr.ne.1)
     .       stop '<# of attributes should be 1 in edge file>'
          return
        else
             do 98 iface = 1,nface
             read(19,fmt=*)idum,(iedgptr(ia,iface),ia=1,3)
   98     continue
          close(19)
          return
        endif
      else
        stop 'invalid file type in subr rtri'
      endif ! on ftype

      close(88)

      return
      end
