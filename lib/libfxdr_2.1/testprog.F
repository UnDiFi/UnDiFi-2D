	program testprog

#ifdef cray
	real x1, x2
#else
	double precision x1, x2
#endif

	common /cb1/ x1, x2

	x1 = 3.0
	x2 = 4.0
	call sub( x1, x2)

	end

c========================================================

	subroutine sub( x1, x2 )

#ifdef cray
	real x1, x2
#else
	double precision x1, x2
#endif

	print *, x1, x2

	return
	end

