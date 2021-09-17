! Update data structure on the boundary to take into account for the presence of phantom points

      subroutine chg_bnd_ptr(
     .                   nodcode,
     .                   ibndptr,
     .                   nbfac,
     .                   inodptr,
     .                   npoin,
     &                   nbpoin)

!     This routine updates the boundary data structure to take into 
!     account possible phantom nodes on the boundary
!
!     ndim  is the space dimension =2
!     nvt = ndim+1 is the number of vertices
!
!     nbfac  boundary faces (shock segments excluded)
!     nbpoin boundary nodes (shock points   excluded)
      implicit none

      integer nbfac,npoin,nbpoin
      integer nodcode(*)
      integer ibndptr(3,nbfac),inodptr(nbpoin,3)

!     ibndptr(1,i) first node of the i-th boundary face
!     ibndptr(2,i) second node of the i-th boundary face
!     ibndptr(3,i) flag of the i-th boundary face; 
!                  if <0 the face has been deactivated
!
!     inodptr is a nodal pointer for boundary nodes
!
!     inodptr(j,1) index of the j-th boundary node
!     inodptr(j,2) first boundary face which the node belong to; 1 <= inodptr(j,2) <= nbfac
!     inodptr(j,3) first boundary face which the node belong to; 1 <= inodptr(j,3) <= nbfac
!
!     the vector ibndptr is read from file *.poly; 
!     the vector inodptr is built from it in setbndrynodeptr() recalled from readmesh()

!     .. local scalars ..
      integer ipoin,ipos,last,ifail,j,k,iface

!     .. local arrays ..
      integer inode(2)

!     .. external subroutines ..
      external dinit,iinit

!     open log file
      open(8,file='log/ChangeBndryPtr.log')

      write(8,*)'Subr ChangeBndryPtr; NBFAC was = ',NBFAC
      write(8,*)'Subr ChangeBndryPtr; NBPOIN was = ',NBPOIN
      write(8,*)'Subr ChangeBndryPtr; NPOIN was = ',NPOIN

      do 1 ipoin = 1, npoin

!     if the node has been deactivated (nodcode=-2) is searched in the vector inodptr
          if( nodcode(ipoin) .eq. -2 )then
              call binsrc(ipoin,inodptr(1,1),nbpoin,ipos,last)
              if(ipos.eq.0)then
                  write(8,*)'Entry NOT found for ',IPOIN
                  write(*,*)'Entry NOT found for ',IPOIN
                  stop
               endif
               write(8,*)'Removing node ',ipoin,
     &                   ' belongs to edges ',(inodptr(ipos,k),k=2,3)
               do 6 k = 2,3

                  iface = inodptr(ipos,k) 

                  if( ibndptr(1,iface) .ne. ipoin )then
                      inode(k-1) = ibndptr(1,iface)
                  else 
                      inode(k-1) = ibndptr(2,iface)
                  endif
    6          continue

!           the face is modified
            iface = inodptr(ipos,2)
            ibndptr(1,iface) = inode(1)
            ibndptr(2,iface) = inode(2)

            write(8,*)'Face ',iface,' has been updated with ',
     &(inode(k),k=1,2)
caldo
            iface = inodptr(ipos,3)
            write(8,*)'face ',iface,' has been removed'
            ibndptr(1,iface) = ibndptr(1,iface)
            ibndptr(2,iface) = ibndptr(2,iface)
            ibndptr(3,iface) =-ibndptr(3,iface)
         endif
    1 continue
      write(8,*)'Subr ChangeBndryPtr; NBFAC is now = ',NBFAC

!     call x04eaf('general',' ',3,nbfac,ibndptr,3,
!    +            'bndry pointer in changebndryptr',ifail)
!     pause
!     stop

      close(8)

      return
      end
