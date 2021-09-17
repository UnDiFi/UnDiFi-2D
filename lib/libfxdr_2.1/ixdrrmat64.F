c	-----------------------------
c	Returns < 0 if error, 0 if OK
c	-----------------------------
	integer function ixdrrmat64( ixdrs, nels, rval )

	implicit none

	integer   ixdrs, ixdrsm1, ierr, nels
	real	  rval(nels)

	ixdrsm1 = ixdrs - 1

	call cxdrrmat64( ixdrsm1, nels, rval, ierr )

	ixdrrmat64 = ierr

	return
	end


