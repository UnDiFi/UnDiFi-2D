#include <stdio.h>
#include <rpc/rpc.h>
#include "cfxdr.h"

extern XDR_element xdrfile[MAX_N_XDR_FILES];

	void
#if (IFSTYLE==1)
	CXDRREWIND( int *ixdrid, int *retval )

#elif  (IFSTYLE==2)
	cxdrrewind( ixdrid, retval )
	int    	*ixdrid;
	int	*retval;

#else
	cxdrrewind_( int *ixdrid, int *retval )
#endif
{
	XDR 	*xdrs;
	u_int	position;

	*retval = 0;	/* No error */

	xdrs     = xdrfile[*ixdrid].xdrs;
	position = 0;

	if( ! xdr_setpos( xdrs, position ) ) {
		if( xdrfile[*ixdrid].return_on_error ) {
			*retval = FXDRERR_REWIND;
			return;
			}
		else
			{
			fprintf( stderr, "FXDR library error while trying to rewind ");
			fprintf( stderr, "file %s\n", 
				xdrfile[*ixdrid].filename );
			exit( -1 );
			}
		}
}
