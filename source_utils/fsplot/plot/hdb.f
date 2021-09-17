      subroutine whdb(coor2,ndim,q2,nnu,nnodes,nnufr,nnodesfr,
     +iconnec,nelem,nelemtypes,logfac,nfac)
c
      implicit none
c 
c NOTE In this example, there are only two variables (named q 
c      and qfr hereafter). q contains the set of variable in 
c      the whole field, it has nnu unknown for each grid point.
c      The variable qfr contains other variable on specific 
c      surface S of the mesh (inflow, outflow, wall,etc~...). 
c      It contains nnufr variables to be stored at the grid 
c      points of the surface S.
c       The unit number(s) below (unit1 and unit2) are arbitrary
c      and to be selected by the participant.
c
c  **************************************************************
c
      integer ndim,nnodes,nelem,nelemtypes,nnu,nfac,nnufr,nnodesfr
      real*8 coor2(ndim,nnodes)
      real*8 q2(nnu,nnodes)
      integer nnufrmax
      parameter (nnufrmax=6)
      real q2fr(4,nnufrmax)
      integer nnperelement(6),iconnec(*)
      integer itypel,itype,nblock,npatch,i,l,nu,iunit,nb,nn,iel,iptr,iv
      integer unit1,ifac
      integer logfac(3,*)
      integer icycl
      external icycl
c
c parameters
c ----------
c
c ndimax               : space dimension (1,2 or 3)
c idimax,jdimax,kdimax : maximum of indices in the i,j,k direction,
c                        structured meshes
c idifrmax,jdifrmax,
c kdifrmax             : maximum of indices in the ifr,jfr,kfr direction,
c                      : structured mesh for S
c nnodesmax            : maximum number of nodes, unstructured meshes
c nnodesfrmax          : maximum number of nodes on S, unstructured
c                        meshes
c nelemax              : maximum number of element, ustructured meshes
c nnumax,
c nnufrmax             : maximum number of variables
c nfacmax              : maximum number of elements on S
c npatchmax            : maximum number of patches (structured meshes)
c 
c arrays
c ------
c coor1,coor2          : x,y,z coordinates 
c q1,q2,q1fr,q2fr      : array of variables
c nnperlement          : number of nodes per element
c ietyp                : define for a given element its type
c iconnec              : connectivity table
c iconnecp(nelemax)    : connectivity table pointer
c iconnecpfac(nfacmax) : connectivity table pointer for the frontier
c nufac(6,nfacmax)     : connectivity table for the frontier
c logfac(nfacmax)      : logic type of the frontier elements
c curv                 : values of variables in file unit2
c npoint               : number of value on the curves curv(*,*)
c ibegfr,...kendfr     : extreme indices of surfaces where surface
c                        values are asked for.
c
c note : the dimension 6*nelemax for iconnec has to be related to the
c size of nnperlement.
c
c   ***************************************************************
c
      unit1=22
      open (unit1,file='soldata',form='unformatted')
c
      if(ndim.eq.2)then
         itypel = 3
      else
         itypel = 2
      endif
c
c     Number of blocks. Note that blocks should be either structured
c     or unstructured, but one database file can contain both
c     structured and unstructured blocks
c 
      nblock = 1
c
      write (unit1) nblock
c
      do 103 nb=1,nblock
c
c     Read type of block, and number of unknowns stored for this block 
c
c     itype = 0: Unstructured mesh
c             1: Structured mesh
c     nnu   = number of unknowns stored
c     nnufr = number of unknown stored on the frontier
c
      itype =0
      npatch =0
c
         write (unit1) itype,nnu,nnufr,npatch
         write (6,*)'itype,nnu,nnufr,npatch'
         write (6,*) itype,nnu,nnufr,npatch
c
c     ******************************************************************
c     Structured mesh
c     ******************************************************************
c
c
c     ******************************************************************
c     Unstructured mesh
c     ******************************************************************
c
c
c     nelemtypes = number of element types
c     nelem      = number of elements
c     nnodes     = number of nodes
c     nfac       = number of facets on the frontier
c     nnodesfr   = number of nodes on the frontier where data is stored
c                  they correspond to type of frontiers that are given
c                  by the next paragraph, problem per problem
c  remark : only on frontier is assumed to be of interest
c           same remarks as for structured
c
            write (unit1) nnodes,nelem,nelemtypes,nfac,nnodesfr,ndim
            write (6,*)'nnodes,nelem,nelemtypes,nfac,nnodesfr,ndim'
            write (6,*) nnodes,nelem,nelemtypes,nfac,nnodesfr,ndim
c
            write (unit1) ((sngl(coor2(l,i)),l=1,ndim),i=1,nnodes)
            write (unit1) ((sngl(q2(nu,i)),i=1,nnodes),nu=1,nnu)
C    MODIFICATION
C            write(unit1) ((q2fr(i,nu),i=1,nnodesfr),nu=1,nnufr)
c
c     Element names implicitly define the element characteristics:
c
c     HE8    8-noded (linear) 3D hex-cube cell      
c     TE4    4-noded (linear) 3D tetrahedral cell
c     T3     3-noded triangular 2D cell
c     Q4     4-noded quadrilateral 2D cell
c     P5     5-noded pyramid element
c     P6     6-noded prismatic element
c
c    itypel      = element type number in table 
c                 (order is  1=HE8, 2=TE4, 3=T3, 4=Q4, 5=P5, 6=P6)
c The ordering of the nodes are defined implicitely, so that
c the faces of the element itypel are known at priory (see figure)
c
           nnperelement(1)=8
           nnperelement(2)=4
           nnperelement(3)=3
           nnperelement(4)=4
           nnperelement(5)=5
           nnperelement(6)=6
c
c  write the connectivity table
c

            do 203 iel=1,nelem
               iptr = (iel-1)*nnperelement(itypel)
               write (unit1) itypel
               write (unit1) (iconnec(iptr+i),i=1,nnperelement(itypel))
  203       continue
c
c  logic of facet
c  
c   logfac= 1 : entrance facet
c   logfac=-1 : symmetry facet
c   logfac= 2 : outflow facet
c   logfac= 3 : freestream facet
c   logfac= 4 : wall nodes (no slip condition)
c   logfac= 5 : wall nodes (velocity=0)
c 
           do 204 ifac=1,nfac
            write(unit1) itypel
            write(unit1) logfac(3,ifac)
            iel=logfac(1,ifac)
            iv =logfac(3,ifac)
            iptr = (iel-1)*nnperelement(itypel)
c           write (unit1) (iconnec(iptr+i),i=1,nnperelement(itypel))
            write (unit1) 
     &      (iconnec(iptr+icycl(iv+i,nnperelement(itypel))),
     &      i=1,nnperelement(itypel)-1)
204        continue
c
c values on frontiers MODIFICATIONS
c
         do 205 ifac=1,nfac
          iel=logfac(1,ifac)
          iv =logfac(3,ifac)
          iptr = (iel-1)*nnperelement(itypel)
          do l = 1,nnperelement(itypel)-1
             i = iconnec(iptr+icycl(iv+l,nnperelement(itypel)-1))
             do nn = 1,nnufr
                q2fr(l,nn) = q2(nn,i)
             end do
          end do
          write(unit1)((q2fr(l,nn),l=1,nnperelement(itypel)-1),
     +                 nn=1,nnufr)
205       continue

c
  103 continue
      close(unit1)
c
      return
      end
