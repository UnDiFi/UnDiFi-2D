c	-----------------------------
c	Returns < 0 if error, 0 if OK
c	-----------------------------
	integer function ixdrreal( ixdrs, rval )

	implicit none

	integer   ixdrs, ixdrsm1, ierr
	real	  rval

	ixdrsm1 = ixdrs - 1

	call cxdrreal( ixdrsm1, rval, ierr )

	ixdrreal = ierr

	return
	end


