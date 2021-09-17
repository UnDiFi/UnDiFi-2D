      subroutine getdia2 (nrow,ncol,job,ja,ia,len,idiag,ioff)
      integer nrow, ncol, job, len, ioff, ia(*), ja(*), idiag(*)
c-----------------------------------------------------------------------
c this subroutine extracts a given diagonal from a matrix stored in csr 
c format. the output matrix may be transformed with the diagonal removed
c from it if desired (as indicated by job.) 
c----------------------------------------------------------------------- 
c our definition of a diagonal of matrix is a vector of length nrow
c (always) which contains the elements in rows 1 to nrow of
c the matrix that are contained in the diagonal offset by ioff
c with respect to the main diagonal. if the diagonal element
c falls outside the matrix then it is defined as a zero entry.
c thus the proper definition of diag(*) with offset ioff is 
c
c     diag(i) = a(i,ioff+i) i=1,2,...,nrow
c     with elements falling outside the matrix being defined as zero.
c 
c----------------------------------------------------------------------- 
c 
c on entry:
c---------- 
c
c nrow	= integer. the row dimension of the matrix a.
c ncol	= integer. the column dimension of the matrix a.
c job   = integer. job indicator.  if job = 0 then
c         the matrix a, ja, ia, is not altered on return.
c         if job.ne.0  then getdia will remove the entries
c         collected in diag from the original matrix.
c         this is done in place.
c
c a,ja,
c    ia = matrix stored in compressed sparse row a,ja,ia,format
c ioff  = integer,containing the offset of the wanted diagonal
c	  the diagonal extracted is the one corresponding to the
c	  entries a(i,j) with j-i = ioff.
c	  thus ioff = 0 means the main diagonal
c
c on return:
c----------- 
c len   = number of nonzero elements found in diag.
c         (len .le. min(nrow,ncol-ioff)-max(1,1-ioff) + 1 )
c
c diag  = real*8 array of length nrow containing the wanted diagonal.
c	  diag contains the diagonal (a(i,j),j-i = ioff ) as defined 
c         above. 
c
c idiag = integer array of  length len, containing the poisitions 
c         in the original arrays a and ja of the diagonal elements
c         collected in diag. a zero entry in idiag(i) means that 
c         there was no entry found in row i belonging to the diagonal.
c         
c a, ja,
c    ia = if job .ne. 0 the matrix is unchanged. otherwise the nonzero
c         diagonal entries collected in diag are removed from the 
c         matrix and therefore the arrays a, ja, ia will change.
c	  (the matrix a, ja, ia will contain len fewer elements) 
c 
c----------------------------------------------------------------------c
c     Y. Saad, sep. 21 1989 - modified and retested Feb 17, 1996.      c 
c----------------------------------------------------------------------c
c     local variables
      integer istart, max, iend, i, kold, k, kdiag, ko
c     
      istart = max(0,-ioff)
      iend = min(nrow,ncol-ioff)
      len = 0
      do 1 i=1,nrow
         idiag(i) = 0
 1    continue
c     
c     extract  diagonal elements
c     
      do 6 i=istart+1, iend
         do 51 k= ia(i),ia(i+1) -1
            if (ja(k)-i .eq. ioff) then
               idiag(i) = k
               len = len+1
               goto 6
            endif
 51      continue
 6    continue
      if (job .eq. 0 .or. len .eq.0) return
c
c     remove diagonal elements and rewind structure
c 
      ko = 0
      do  7 i=1, nrow 
         kold = ko
         kdiag = idiag(i) 
         do 71 k= ia(i), ia(i+1)-1 
            if (k .ne. kdiag) then
               ko = ko+1
               ja(ko) = ja(k)
            endif
 71      continue
         ia(i) = kold+1
 7    continue
c
c     redefine ia(nrow+1)
c
      ia(nrow+1) = ko+1
      return
c------------end-of-getdia2----------------------------------------------
c-----------------------------------------------------------------------
      end
