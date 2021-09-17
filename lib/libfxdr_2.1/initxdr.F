c	-------------------------------------------------
c	Returns < 0 on error, otherwise the xdrid used to
c	access the xdr file 
c
c	Args:
c		filename: character*(*): name of file to 
c			access
c		mode: character*1: either 'r', 'w', or 'a' 
c			(for read, write, and append)
c		returnonerror: logical: if .TRUE., then
c			routines return even if there is
c			an I/O error.  If .FALSE., then
c			routines halt on I/O error, 
c			printing out a useful message 
c			(this is like what fortran does,
c			so set to .FALSE. for most fortran-
c			like behavior)
c	-------------------------------------------------
	integer function initxdr( filename, mode, returnonerror )

	implicit none

	character*(*) 	filename
	character*1  	mode
	logical		returnonerror

	integer     	iretval, cxdrinit, iroe, imode, lenf, istrlen
	integer		ifxdrstrlen

	if( (mode .eq. 'r').or.(mode .eq. 'R') ) then 
		imode = 1
	else if( (mode .eq. 'w').or.(mode .eq. 'W') ) then
		imode = 2
	else if( (mode .eq. 'a').or.(mode .eq. 'A') ) then
		imode = 3
	else
		write(0,*) 'Error: fxdr library, initxdr called '
		write(0,*) 'with unknown mode (should be r, w, or a):',
     &				mode
		stop 'fxdr library: initxdr routine, bad mode'
	endif

	lenf    = len(filename)
	istrlen = ifxdrstrlen( filename )
	if( istrlen .lt. lenf ) then
		lenf = istrlen
	endif

	if( returnonerror ) then
		iroe = 1
	else
		iroe = 0
	endif

	iretval = cxdrinit( lenf, filename, imode, iroe )
	if( iretval .lt. 0 ) then
		initxdr = iretval
		return
	endif

c	-----------------------------------------------------------
c	Add one to the returned value because the C library returns 
c	starting with zero, which is awkward considering that 
c	unitialized variables often start with zero in Fortran.
c	-----------------------------------------------------------
	initxdr = iretval + 1

	return
	end

c========================================================================

	integer function ifxdrstrlen( s )

	character*(*) s

	i = 1
100	continue
		if( s(i:i) .eq. ' ' ) then
			goto 200
		endif
		i = i + 1
		if( i .gt. len(s) ) then 
			goto 200
		endif
		goto 100

200	continue
	ifxdrstrlen = i-1

	return
	end
	
