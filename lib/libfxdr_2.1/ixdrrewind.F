c	-----------------------------
c	Returns < 0 if error, 0 if OK
c	-----------------------------
	integer function ixdrrewind( ixdrs )

	implicit none

	integer   ixdrs, ixdrsm1, ierr

	ixdrsm1 = ixdrs - 1

	call cxdrrewind( ixdrsm1, ierr )

	ixdrrewind = ierr

	return
	end
