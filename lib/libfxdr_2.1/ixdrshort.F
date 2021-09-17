c	-----------------------------
c	Returns < 0 if error, 0 if OK
c	-----------------------------
      	integer function ixdrshort (ixdrs, num)

	implicit none

      	integer ixdrs,ixdrsm1,num, ierr

      	ixdrsm1=ixdrs - 1

      	call cxdrshort( ixdrsm1, num, ierr )

	ixdrshort = ierr

      	return
      	end

