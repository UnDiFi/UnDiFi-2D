	program testfxdr

c	===========================================================
c	Test the FXDR library by exercising the calls, and
c	comparing the results to the known good .xdr file.
c
c	David Pierce   Scripps Inst. Oceanography / Climate Research
c	Aug 10, 1999
c	===========================================================

	implicit none

	include 'fxdr.inc'

	integer ien, jen
	parameter(ien=3,jen=2)

	integer iwrite, ierr, isign, j, i2, imatwrite, ixdrs, iread
	real*4 	rwrite, rsign, rxpon, fmatwrite, rread
#ifdef cray
	real	dwrite
#else
	double precision dwrite
#endif

	dimension imatwrite(ien)
	dimension fmatwrite(ien,jen)
#ifdef cray
	real dmatwrite(ien,jen), dsign, dxpon
#else
	double precision dmatwrite(ien,jen), dsign, dxpon
#endif
	character*26 stringwrite

	common /cb1/ iwrite, imatwrite
	common /cb2/ rwrite, fmatwrite
	common /cb3/ dwrite, dmatwrite
	common /cb4/ stringwrite

c	-----------------------------------------------
c	Initialize our variables to recognizable values
c	-----------------------------------------------
	iwrite = 987654
	rwrite = 3.1415926535e0
#ifdef cray
	dwrite = 1.234567890123456789e0
#else
	dwrite = 1.234567890123456789d0
#endif
	stringwrite = 'this is a FXDR test string'

c	-------------------------------------------------
c	Fill in the integer, single precision, and double
c	precision arrays with useful values.
c	-------------------------------------------------
	dsign = 1.0d0
	rsign = 1.0
	isign = 1
	rxpon = -20.0
	dxpon = -20.0d0
	do j=1, jen
	do i2=1, ien
		imatwrite(i2)   = isign*i2
		fmatwrite(i2,j) = rsign*(float(i2) + 
     &				float(j)/100.0)*10.0**rxpon
#ifdef cray
		dmatwrite(i2,j) = dsign*(float(i2)  + 
     &				float(j)/10000000000.0)*10.0**dxpon
#else
		dmatwrite(i2,j) = dsign*(dble(i2)  + 
     &				dble(j)/10000000000.0)*10.0**dxpon
#endif
		dsign           = -1.0d0*dsign
		rsign           = -1.0*rsign
		isign           = -isign 
		rxpon           = rxpon + 7.0
		dxpon           = dxpon + 7.0d0
	enddo
	enddo

	print *, 'here are the check values:'
	print *, fmatwrite
	print *, dmatwrite

c	-------------------------------------------------
c	Open the test XDR file.  Set it to halt on errors
c	(.FALSE.), so don't bother checking error codes 
c	below.
c	-------------------------------------------------
	print *, 'generating test XDR file test.xdr...'
	ixdrs = initxdr( 'test.xdr', 'w', .FALSE. )

c	------------------------------------------------------
c	Put the scalar values: integer, real, double precision
c	------------------------------------------------------
	ierr = ixdrint   ( ixdrs, iwrite )
	ierr = ixdrreal  ( ixdrs, rwrite )
#ifdef cray
	ierr = ixdrreal64( ixdrs, dwrite )
#else
	ierr = ixdrdouble( ixdrs, dwrite )
#endif

c	-----------------------------------------------------
c	Put the array values: integer, real, double precision
c	-----------------------------------------------------
	ierr = ixdrimat  ( ixdrs, ien,     imatwrite )
	ierr = ixdrrmat  ( ixdrs, ien*jen, fmatwrite )
#ifdef cray
	ierr = ixdrrmat64( ixdrs, ien*jen, dmatwrite )
#else
	ierr = ixdrdmat  ( ixdrs, ien*jen, dmatwrite )
#endif

	ierr = ixdrstring( ixdrs, stringwrite )

	ierr = ixdrclose( ixdrs )
	print *, 'successfully wrote the test file test.xdr'

	call compareto( 'test.xdr' )
	call compareto( 'test.orig.xdr' )

	print *
	print *, 'Testing the rewind function...'
	ixdrs = initxdr( 'test.xdr', 'r', .FALSE. )
	ierr = ixdrint   ( ixdrs,  iread )
	call checki   ( iwrite, iread )
	ierr = ixdrreal  ( ixdrs,  rread )
	call checkr   ( rwrite, rread )
	ierr = ixdrrewind( ixdrs )
	ierr = ixdrint   ( ixdrs,  iread )
	call checki   ( iwrite, iread )
	ierr = ixdrreal  ( ixdrs,  rread )
	call checkr   ( rwrite, rread )
	ierr = ixdrclose ( ixdrs )

c	------------------------------------------
c	Now check the error handling functionality
c	------------------------------------------
	print *
	print *, '---------------------------------------------'
	print *, 'Testing the error handling functionality ... '
	print *, 'these routine should NOT STOP!'
	print *, '---------------------------------------------'
	print *

c	----------------------------------------
c	Try to open a read-only file for writing
c	----------------------------------------
	print *, 'opening a read-only file for writing...'
	ixdrs = initxdr( 'test_read_only.xdr', 'w', .TRUE. )
	if( ixdrs .ge. 0 ) then
		print *, 'TEST FAILED!  No error returned when trying'
		print *, 'to open a read only file for writing!!'
		print *, 'Make SURE file test_read_only.xdr does NOT'
		print *, 'have write permissions!'
		print *, '--------- TEST FAILED ----------'
		stop
	endif

c	----------------------------------------
c	Try to open a no-access file for reading
c	----------------------------------------
	print *, 'opening a no-access file for reading...'
	ixdrs = initxdr( 'test_no_read.xdr', 'r', .TRUE. )
	if( ixdrs .ge. 0 ) then
		print *, 'TEST FAILED!  No error returned when trying'
		print *, 'to open a no-access file for reading!!'
		print *, 'Make SURE file test_no_read.xdr does NOT'
		print *, 'have read permissions!'
		print *, '--------- TEST FAILED ----------'
		stop
	endif

c	-------------------------------------------------------------
c	Write a single integer, rewind, then try to read TWO integers
c	-------------------------------------------------------------
	print *, 'trying to read integer past EOF...'
	ixdrs = initxdr( 'test.xdr', 'w', .TRUE. )
	ierr  = ixdrint( ixdrs, iwrite )
	ierr  = ixdrclose( ixdrs )
	ixdrs = initxdr( 'test.xdr', 'r', .TRUE. )
	ierr  = ixdrint( ixdrs, iread )
	if( ierr .lt. 0 ) then
		print *, 'TEST FAILED!  Anomalous error returned when'
		print *, 'trying to read a test integer from file'
		print *, 'test.xdr!'
		print *, '--------- TEST FAILED ----------'
		stop
	endif
	ierr  = ixdrint( ixdrs, iread )
	if( ierr .ne. FXDRERR_READERR ) then
		print *, 'TEST FAILED!  No error returned when trying'
		print *, 'to read a test integer past the end of file!'
		print *, '--------- TEST FAILED ----------'
		stop
	endif

c	-------------------------------------------------------
c	Write a single real, rewind, then try to read TWO reals
c	-------------------------------------------------------
	print *, 'trying to read real past EOF...'
	ixdrs = initxdr( 'test.xdr', 'w', .TRUE. )
	ierr  = ixdrreal( ixdrs, rwrite )
	ierr  = ixdrclose( ixdrs )
	ixdrs = initxdr( 'test.xdr', 'r', .TRUE. )
	ierr  = ixdrreal( ixdrs, rread )
	if( ierr .lt. 0 ) then
		print *, 'TEST FAILED!  Anomalous error returned when'
		print *, 'trying to read a test real from file'
		print *, 'test.xdr!'
		print *, '--------- TEST FAILED ----------'
		stop
	endif
	ierr  = ixdrint( ixdrs, iread )
	if( ierr .ne. FXDRERR_READERR ) then
		print *, 'TEST FAILED!  No error returned when trying'
		print *, 'to read a test real past the end of file!'
		print *, '--------- TEST FAILED ----------'
		stop
	endif

c	-----------------------------------------------------
c	Try to read wrong # of elements from an integer array
c	-----------------------------------------------------
	print *,'try to read wrong # of elements from integer array...'
	ixdrs = initxdr( 'test.xdr', 'w', .FALSE. )
	ierr = ixdrimat  ( ixdrs, ien,     imatwrite )
	ierr  = ixdrclose( ixdrs )
	ixdrs = initxdr( 'test.xdr', 'r', .TRUE. )
	ierr = ixdrimat  ( ixdrs, 2*ien,     imatwrite )
	if( ierr .ge. 0 ) then
		print *, 'TEST FAILED!  No error returned when trying'
		print *, 'to read incorrect # integer elements from'
		print *, 'an array!'
		print *, '--------- TEST FAILED ----------'
		stop
	endif

c	-----------------------------------------------------
c	Try to read wrong # of elements from an real array
c	-----------------------------------------------------
	print *, 'try to read wrong # of elements from a real array...'
	ixdrs = initxdr( 'test.xdr', 'w', .FALSE. )
	ierr = ixdrrmat  ( ixdrs, ien*jen, fmatwrite )
	ierr  = ixdrclose( ixdrs )
	ixdrs = initxdr( 'test.xdr', 'r', .TRUE. )
	ierr = ixdrrmat  ( ixdrs, ien, fmatwrite )
	if( ierr .ge. 0 ) then
		print *, 'TEST FAILED!  No error returned when trying'
		print *, 'to read incorrect # real elements from an'
		print *, 'array!'
		print *, 'ierr=', ierr
		print *, '--------- TEST FAILED ----------'
		stop
	endif

c	--------------------
c	Test append function
c	--------------------
	print *
	print *, '--------------------------'
	print *, 'Testing append function...'
	print *, '--------------------------'
	print *
c	... make test file with single integer ...
	ixdrs = initxdr( 'test.xdr', 'w', .FALSE. )
	ierr  = ixdrint( ixdrs, iwrite )
	ierr  = ixdrclose( ixdrs )
c	... read back same integer ...
	ixdrs = initxdr( 'test.xdr', 'r', .TRUE. )
	ierr  = ixdrint( ixdrs, iread )
	if( ierr .lt. 0 ) then
		print *, 'ERROR reading first int from test file!!'
		print *, '-------- TEST FAILED -------'
		stop
	endif
	call checki( iwrite, iread )
c	... read back ANOTHER real, should fail ...
	ierr  = ixdrreal( ixdrs, rread )
	if( ierr .ge. 0 ) then
		print *, 'ERROR: no error value returned when reading'
		print *, 'past end of file!!'
		print *, '-------- TEST FAILED -------'
		stop
	endif
	ierr  = ixdrclose( ixdrs )
c	... now try to append real to this file ...
	ixdrs = initxdr( 'test.xdr', 'a', .TRUE. )
	ierr  = ixdrreal( ixdrs, rwrite )
	if( ierr .lt. 0 ) then
		print *, 'ERROR trying to append real to test file!'
		print *, '-------- TEST FAILED -------'
		stop
	endif
	ierr  = ixdrclose( ixdrs )
c	... check to make sure appendation (?) went OK ...
	ixdrs = initxdr( 'test.xdr', 'r', .FALSE. )
	ierr  = ixdrint( ixdrs, iread )
	call checki( iwrite, iread )
	ierr  = ixdrreal( ixdrs, rread )
	call checkr( rwrite, rread )
	ierr  = ixdrclose( ixdrs )

	print *
	print *, '---------------------------------------------------'
	print *, 'All tests passed successfully: precision loss was'
	print *, 'never more than 1 part in 10**-10 for double'
	print *, 'precision nor more than 1 part in 10**-5 for single'
	print *, 'precision.'
	print *, '---------------------------------------------------'

	end

c=====================================================================

	subroutine compareto( filename )

	character*(*) filename

	parameter(ien=3,jen=2)

	integer iread, iwrite
	real*4 	rread, rwrite, fmatwrite, fmatread
#ifdef cray
	real 	dread, dwrite
#else
	double precision dread, dwrite
#endif

	common /cb1/ iwrite, imatwrite
	common /cb2/ rwrite, fmatwrite
	common /cb3/ dwrite, dmatwrite
	common /cb4/ stringwrite

	dimension imatwrite(ien), imatread(ien)
	dimension fmatwrite(ien,jen), fmatread(ien,jen)
#ifdef cray
	real	dmatwrite(ien,jen), dmatread(ien,jen)
#else
	double precision dmatwrite(ien,jen), dmatread(ien,jen)
#endif
	character*26 stringwrite, stringread

	print *
	print *, '--------------------------------------------'
	print *, 'checking file ', filename, ' for correctness.'
	print *, '--------------------------------------------'

c	-----------------------------------------
c	Open the file and try to read it back in.
c	-----------------------------------------
	ixdrs = initxdr( filename, 'r', .FALSE. )

c	-------------------------------------------------------
c	Read the scalar values: integer, real, double precision
c	-------------------------------------------------------
	ierr = ixdrint   ( ixdrs,  iread )
	call checki   ( iwrite, iread )
	ierr = ixdrreal  ( ixdrs,  rread )
	call checkr   ( rwrite, rread )
#ifdef cray
	ierr = ixdrreal64( ixdrs,  dread )
#else
	ierr = ixdrdouble( ixdrs,  dread )
#endif
	call checkd   ( dwrite, dread )

c	------------------------------------------------------
c	Read the array values: integer, real, double precision
c	------------------------------------------------------
	ierr = ixdrimat  ( ixdrs, ien,     imatread )
	do i=1, ien
		call checki( imatwrite(i), imatread(i) )
	enddo
	ierr = ixdrrmat  ( ixdrs, ien*jen, fmatread )
	do j=1, jen
	do i=1, ien
		call checkr( fmatwrite(i,j), fmatread(i,j) )
	enddo
	enddo
#ifdef cray
	ierr = ixdrrmat64( ixdrs, ien*jen, dmatread )
#else
	ierr = ixdrdmat  ( ixdrs, ien*jen, dmatread )
#endif
	do j=1, jen
	do i=1, ien
		call checkd( dmatwrite(i,j), dmatread(i,j) )
	enddo
	enddo

	ierr = ixdrstring( ixdrs, stringread )
	do i=1, len(stringread)
		if( stringread(i:i) .ne. stringwrite(i:i) ) then 
			print *, 'Error on string read! I expected: ', 
     &				stringwrite(i:i)
			print *, 'but got: ', stringread(i:i)
		endif
	enddo

	ierr = ixdrclose( ixdrs )
	
	return
	end

c====================================================================

	subroutine checkr( expected, got )

	real*4 expected, got

	if( got .ne. expected ) then
		diffnorm = abs((got-expected)/(.5*(got+expected)))
		if( diffnorm .gt. 1.0e-5 ) then
	print *, 'Fatal error!  On single precision read, I got ',
     &				got
			print *, 'but expected ', expected
	print *, 'The library failed its tests--please do not use it!'
			stop
		else
			print *, 'checking single precision...'
	print *, '===> precision loss of 1 part in 10**', 
     &				nint(log10(1.0/diffnorm))
		endif
	else
		print *, 'checking single precision...OK'
	endif

	return
	end

c====================================================================

	subroutine checkd( expected, got )

#ifdef cray
	real	expected, got, diffnorm
#else
	double precision expected, got, diffnorm
#endif

	if( got .ne. expected ) then
#ifdef cray
		diffnorm = abs((got-expected)/(.5d0*(got+expected)))
#else
		diffnorm = dabs((got-expected)/(.5d0*(got+expected)))
#endif
		if( diffnorm .gt. 1.0d-10 ) then
	print *, 'Fatal error!  On double precision read, I got ',
     &				got
			print *, 'but expected ', expected
	print *, 'The library failed its tests--please do not use it!'
			stop
		else
			print *, 'checking double precision...'
#ifdef cray
			iprecis = nint(log10(1.0/diffnorm))
#else
			iprecis = nint(dlog10(1.0/diffnorm))
#endif
	print *, '===> precision loss of 1 part in 10**', iprecis
		endif
	else
		print *, 'checking double precision...OK'
	endif

	return
	end

c====================================================================

	subroutine checki( expected, got )

	integer expected, got

	if( got .ne. expected ) then
		print *, 'Fatal error! On integer read, I got ', got
		print *, 'but expected ', expected
		print *, 'The library failed its tests--'
		print *, 'please do not use it!'
		stop
	endif
	print *, 'checking integer...OK'

	return
	end
