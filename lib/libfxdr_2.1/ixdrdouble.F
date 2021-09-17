c	-----------------------------
c	Returns < 0 if error, 0 if OK
c	-----------------------------
	integer function ixdrdouble( ixdrs, dval )

	implicit none

	integer   ixdrs, ixdrsm1, ierr
	double precision  dval

	ixdrsm1 = ixdrs - 1

	call cxdrdouble( ixdrsm1, dval, ierr )

	ixdrdouble = ierr

	return
	end


