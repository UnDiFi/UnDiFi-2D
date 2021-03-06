      subroutine dump (i1,i2,values,a,ja,ia,iout)
      integer i1, i2, ia(*), ja(*), iout
      real*8 a(*) 
      logical values 
c-----------------------------------------------------------------------
c outputs rows i1 through i2 of a sparse matrix stored in CSR format 
c (or columns i1 through i2 of a matrix stored in CSC format) in a file, 
c one (column) row at a time in a nice readable format. 
c This is a simple routine which is useful for debugging. 
c-----------------------------------------------------------------------
c on entry:
c---------
c i1    = first row (column) to print out
c i2    = last row (column) to print out 
c values= logical. indicates whether or not to print real values.
c         if value = .false. only the pattern will be output.
c a,
c ja, 
c ia    =  matrix in CSR format (or CSC format) 
c iout  = logical unit number for output.
c---------- 
c the output file iout will have written in it the rows or columns 
c of the matrix in one of two possible formats (depending on the max 
c number of elements per row. The values are output with only 
c two digits of accuracy (D9.2). )
c-----------------------------------------------------------------------
c     local variables
      integer maxr, i, k1, k2 
c
c select mode horizontal or vertical 
c
        maxr = 0
        do 1 i=i1, i2
           maxr = max0(maxr,ia(i+1)-ia(i))
 1      continue
        
        if (maxr .le. 8) then
c
c able to do one row acros line
c
        do 2 i=i1, i2
           write(iout,100) i
	   k1=ia(i)
	   k2 = ia(i+1)-1
	   write (iout,101) (ja(k),k=k1,k2)
	   if (values) write (iout,102) (a(k),k=k1,k2)
 2      continue
      else 
c
c unable to one row acros line. do three items at a time
c across a line 
         do 3 i=i1, i2
            if (values) then
               write(iout,200) i
            else
               write(iout,203) i               
            endif
            k1=ia(i)
            k2 = ia(i+1)-1
            if (values) then
               write (iout,201) (ja(k),a(k),k=k1,k2)
            else
                write (iout,202) (ja(k),k=k1,k2)
             endif
 3       continue
      endif 
c
c formats :
c
 100  format (1h ,34(1h-),' row',i6,1x,34(1h-) )
 101  format(' col:',8(i5,6h     : ))
 102  format(' val:',8(D9.2,2h :) )
 200  format (1h ,30(1h-),' row',i3,1x,30(1h-),/
     *     3('  columns :    values  * ') )
c-------------xiiiiiihhhhhhddddddddd-*-
 201  format(3(1h ,i6,6h   :  ,D9.2,3h * ) )
 202  format(6(1h ,i5,6h  *    ) ) 
 203  format (1h ,30(1h-),' row',i3,1x,30(1h-),/
     *     3('  column  :  column   *') )
      return
c----end-of-dump--------------------------------------------------------
c-----------------------------------------------------------------------
      end
