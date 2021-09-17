      subroutine renumber(xy,ndim,np,zroe,wksp,irank,ndof,
     &icelnod,nodcode,nofvert,nelem)
      implicit none
      integer ndim,nbkgrd,np,ndof,nelem,nofvert
      double precision xy(ndim,*),zroe(ndof,*),wksp(*) 
      integer icelnod(nofvert,nelem),irank(*),nodcode(*)
      integer nitems
      integer i,j,ia,ielem,ifail
C
C     reading rank
C
      OPEN(UNIT=23,FILE='rank.dat')
      do i = 1,np
         read(23,*)j,j,j,irank(i)
      enddo
      CLOSE(23)
cnag  CALL M01ZBF(IRANK,1,NP,IFAIL) ! CHECKS THE VALIDITY OF A PERMUTATION
cnag  WRITE (6,FMT=*) 'M01ZBF (check) has returned IFAIL =',IFAIL
c
      DO 18 IELEM = 1,NELEM
          DO 18 J = 1,nofvert
              IA = icelnod(J,IELEM)
              icelnod(J,IELEM) = IRANK(IA)
   18 CONTINUE
!     WRITE (6,FMT=*) 'rank',(irank(i),i=1,npoin)
      WRITE (6,FMT=*) 'ranking nodal coordinates'
      DO 30 I = 1,ndim
          CALL DCOPY(NP,XY(I,1),ndim,wksp,1)
!         WRITE (6,FMT=*) 'rv',i,(rv(j),j=1,npoin)
cnag      CALL M01EAF(wksp,1,NP,IRANK,IFAIL)
          CALL R8RANK(wksp,1,NP,IRANK,IFAIL)
          WRITE (6,FMT=*) 'R8RANK (ranking coords) has returned IFAIL ='
     +,IFAIL
!         WRITE (6,FMT=*) 'rv',i,(rv(j),j=1,npoin)
          CALL DCOPY(NP,wksp,1,XY(I,1),ndim)
   30 CONTINUE
c
      WRITE (6,FMT=*) 'ranking nodal values'
      DO 32 I = 1,ndof
          CALL DCOPY(NP,ZROE(I,1),ndof,wksp,1)
cnag      CALL M01EAF(wksp,1,NP,IRANK,IFAIL)
          CALL R8RANK(wksp,1,NP,IRANK,IFAIL)
          WRITE (6,FMT=*) 'R8RANK (ranking values) has returned IFAIL ='
     +,IFAIL
          CALL DCOPY(NP,wksp,1,ZROE(I,1),ndof)
   32 CONTINUE
      WRITE (6,FMT=*) 'ranking nodal values'
cnag  CALL M01EBF(nodcode,1,NP,IRANK,IFAIL)
      CALL I4RANK(nodcode,1,NP,IRANK,IFAIL)
          WRITE (6,FMT=*) 'I4RANK (ranking coords) has returned IFAIL ='
     +,IFAIL
C
      return
      end
