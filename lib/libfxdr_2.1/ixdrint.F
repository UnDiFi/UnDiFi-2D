c	-----------------------------
c	Returns < 0 if error, 0 if OK
c	-----------------------------
	integer function ixdrint( ixdrs, ival )

	implicit none

	integer   ixdrs, ixdrsm1, ierr
	integer	  ival

	ixdrsm1 = ixdrs - 1
	call cxdrint( ixdrsm1, ival, ierr )

	ixdrint = ierr

	return
	end


