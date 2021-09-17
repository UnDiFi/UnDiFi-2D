c	-----------------------------
c	Returns < 0 if error, 0 if OK
c	-----------------------------
	integer function ixdrclose( ixdr )

	implicit none

	integer ixdr, ixdrm1, ierr

	ixdrm1 = ixdr - 1
	call cxdrclose( ixdrm1, ierr )

	ixdrclose = ierr

	return
	end
