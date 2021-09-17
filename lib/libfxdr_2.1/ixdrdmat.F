c	-----------------------------
c	Returns < 0 if error, 0 if OK
c	-----------------------------
	integer function ixdrdmat( ixdrs, nels, dval )

	implicit none

	integer   ixdrs, ixdrsm1, ierr, nels
	double precision  dval(nels)

	ixdrsm1 = ixdrs - 1
	call cxdrdmat( ixdrsm1, nels, dval, ierr )

	ixdrdmat = ierr

	return
	end


