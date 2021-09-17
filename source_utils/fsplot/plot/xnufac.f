      subroutine xnufac(nbfac,icelnod,nofvert,ibndfac,nufac,logfac)

      implicit none

      integer nbfac,nofvert
      integer icelnod(nofvert,*),ibndfac(3,nbfac),
     +nufac(nofvert-1,nbfac),logfac(nbfac)
      integer icycl
      integer icolr,jvert,j,iface,ielem,ivert

      do 1 iface = 1,nbfac
          ielem = ibndfac(1,iface)
          ivert = ibndfac(2,iface)
          icolr = ibndfac(3,iface)
          logfac(iface)=icolr+10
          do 3 j = 1,nofvert-1
             jvert = icycl(j+ivert,nofvert)
             nufac(j,iface) = icelnod(jvert,ielem)
    3    continue
    1 continue

      return
      end
