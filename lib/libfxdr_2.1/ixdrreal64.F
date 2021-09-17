c	-----------------------------
c	Returns < 0 if error, 0 if OK
c	-----------------------------
	integer function ixdrreal64( ixdrs, rval )

	implicit none

	integer   ixdrs, ixdrsm1, ierr
	real	  rval

	ixdrsm1 = ixdrs - 1

	call cxdrreal64( ixdrsm1, rval, ierr )

	ixdrreal64 = ierr

	return
	end


