      subroutine wmesh(ndimax,idimax,jdimax,kdimax,nsmax,neltmax,
     &           nfacmax,npatchmax,iunit,
     &           ndim,nblock,idim,jdim,kdim,istruc,npatchi,npatchf,
     &           ns,nfac,nelt,ihmgel,ihmgfac,itypel,itypfac,
     &           coor1,coor2,nnperel,nuvoltyp,nuvol,logvol,
     &           nuvolp,nufacp,nufactyp,nufac,logfac, 
     &           ibegi,iendi,jbegi,jendi,kbegi,kendi,
     &           ibegf,iendf,jbegf,jendf,kbegf,kendf)
c
      real*4 coor1(ndimax,idimax,jdimax,kdimax),coor2(ndimax,nsmax)
      integer nnperel(7),iunit
      integer nuvoltyp(neltmax),nuvol(8*neltmax),logvol(neltmax)
      integer nufactyp(nfacmax),nufac(4*nfacmax),logfac(nfacmax)
      integer nuvolp(neltmax),nufacp(nfacmax)
      integer ibegi(npatchmax),iendi(npatchmax),
     &        jbegi(npatchmax),jendi(npatchmax),
     &        kbegi(npatchmax),kendi(npatchmax)
      integer ibegf(npatchmax),iendf(npatchmax),
     &        jbegf(npatchmax),jendf(npatchmax),
     &        kbegf(npatchmax),kendf(npatchmax)
c=============================================================================
c
c Definition of element or boundary element types : itypel,itypfac 
c The type (from 1 to 7) give implicitely the number of nodes per element 
c     itypX=1 ==> HE8 8-noded 3D hex-cube element 
      nnperel(1)=8
c     itypX=2 ==> P6 6-noded 3D prismatic element 
      nnperel(2)=6
c     itypX=3 ==> P5 5-noded 3D pyramidal element 
      nnperel(3)=5
c     itypX=4 ==> TE4 4-noded 3D tetrahedral element 
      nnperel(4)=4
c     itypX=5 ==> Q4 4-noded 2D quadrangular element 
      nnperel(5)=4
c     itypX=6 ==> T3 3-noded 2D triangular element 
      nnperel(6)=3
c     itypX=7 ==> S2 2-noded 1D edge element 
      nnperel(7)=2
c
c  2) Type of each block  : 
c  ========================
c     istruc=1 Structured Block
c     istruc=0 Unstructured Block
      WRITE(iunit) istruc
cform WRITE(6,*) istruc
c
c=====Case of structured block :
      if(istruc.eq.1) then
c
c   3a) Number of interior patches, number of boundary patches :
      WRITE(iunit) npatchi,npatchf 
c
c   3b) Number of points in X1,X2[,X3] directions : (In 2D, kdim=1) 
      WRITE(iunit) idim,jdim,kdim
      if(ndim.eq.2) kdim=1
c Number of nodes
      ns=idim*jdim*kdim
c
c   3c) Coordinates : 
      WRITE(iunit)((((coor1(l,i,j,k),l=1,ndim)
     &                ,i=1,idim),j=1,jdim),k=1,kdim)
c
c   3d) for each interior patch, indices min et max in each direction :
c  ( In 2D, kbegi=kendi=1 )
      do ip=1,npatchi
      WRITE(iunit)ibegi(ip),iendi(ip),jbegi(ip),jendi(ip)
     &           ,kbegi(ip),kendi(ip),logvol(ip)
      enddo
c
c   3e) for each boundary patch, indices min et max in each direction :
c  ( In 2D, kbegf=kendf=1 )
      do ip=1,npatchf
      WRITE(iunit)ibegf(ip),iendf(ip),jbegf(ip),jendf(ip)
     &           ,kbegf(ip),kendf(ip),logfac(ip)
      enddo
c
c=====end of istruc=1
      endif
c
c=====Case of unstructured block :
      if(istruc.eq.0) then
c 
c   4a) Numbers of nodes, elements, boundary elements :
      WRITE(iunit) ns,nelt,nfac
cform WRITE(6,*) ns,nelt,nfac
c
c   4b) Homegeneous element type, homegeneous boundary element type :
c     ihmgX=1 Homegeneous 
c     ihmgX=0 Non-Homegeneous 
      WRITE(iunit) ihmgel,ihmgfac
cform WRITE(6,*) ihmgel,ihmgfac
c
c   4c) Coordinates:
      WRITE(iunit) ((coor2(l,i),l=1,ndim),i=1,ns)
cform WRITE(6,*) ((coor2(l,i),l=1,ndim),i=1,ns)
c
c-----Case of homogeneous element type :
      if(ihmgel.eq.1) then
c
c   4d1) Type of elements : (see the definitions above)
      WRITE(iunit) itypel 
cform WRITE(6,*) itypel 
c
c   4d2) For each element, nodal connectivity and logic of element :
      iptr=0
      do iel=1,nelt
      WRITE(iunit)(nuvol(iptr+i),i=1,nnperel(itypel)),logvol(iel) 
cform WRITE(6,*)(nuvol(iptr+i),i=1,nnperel(itypel)),logvol(iel) 
      iptr=iptr+nnperel(itypel)
      enddo
c
c-----end of ihmgel=1
      endif
c   
c-----Case of non-homogeneous element type :
      if(ihmgel.eq.0) then
c
      iptr=0
      do iel=1,nelt
      nuvolp(iel)=iptr+1
      itypel=nuvoltyp(iel)
c   4e1) Type of elements : 
      WRITE(iunit) itypel
c   4e2) nodal connectivity and logic of element :
      WRITE(iunit)(nuvol(iptr+i),i=1,nnperel(itypel)),logvol(iel) 
      iptr=iptr+nnperel(itypel)
      enddo
c
c-----end of ihmgel=0
      endif
c   
c-----Case of homogeneous boundary element type :
      if(ihmgfac.eq.1) then
c
c   4f1) Type of elements : (see the definitions above)
      WRITE(iunit) itypfac 
cform WRITE(6,*) itypfac 
c
c   4f2) For each element, nodal connectivity and logic of element :
      iptr=0
      do ifac=1,nfac
      WRITE(iunit)(nufac(iptr+i),i=1,nnperel(itypfac)),logfac(ifac) 
cform WRITE(6,*)(nufac(iptr+i),i=1,nnperel(itypfac)),logfac(ifac) 
      iptr=iptr+nnperel(itypfac)
      enddo
c
c-----end of ihmgfac=1
      endif
c   
c-----Case of non-homogeneous boundary element type :
      if(ihmgfac.eq.0) then
c
      iptr=0
      do ifac=1,nfac
      nufacp(ifac)=iptr+1
      itypfac=nufactyp(ifac)
c   4g1) Type of elements : 
      WRITE(iunit) itypfac
c   4g2) nodal connectivity and logic of element :
      WRITE(iunit)(nufac(iptr+i),i=1,nnperel(itypfac)),logfac(ifac) 
      iptr=iptr+nnperel(itypfac)
      enddo
c
c-----end of ihmgfac=0
      endif
c   
c=====end of istruc=0
      endif
c
      return
      end
