#include <stdio.h>
#include <rpc/rpc.h>
#include "cfxdr.h"

extern XDR_element xdrfile[MAX_N_XDR_FILES];

	void

#if (IFSTYLE==1)
	CXDRCLOSE( int *ixdrid, int *retval )

#elif (IFSTYLE==2)
	cxdrclose( ixdrid, retval )
	int	*ixdrid;
	int	*retval;

#else
	cxdrclose_( int *ixdrid, int *retval )
#endif
{
	XDR 	*xdrs;
	FILE	*f;

	*retval = 0;	/* No errors */

	xdrs = xdrfile[*ixdrid].xdrs;

	/* Get the file pointer, which was saved in the user's
	 * data area, so we can close the file after destroying
	 * the XDR handle.
	 */
	f = (FILE *)xdrs->x_public;

	xdr_destroy( xdrs );

	fclose( f );

	xdrfile[*ixdrid].filename = NULL;
}

