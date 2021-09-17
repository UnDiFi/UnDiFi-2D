#include <stdio.h>
#include <rpc/rpc.h>
#include "cfxdr.h"

extern XDR_element xdrfile[MAX_N_XDR_FILES];

	void
#if defined(cray_twobytecptrs)
	CXDRSTRING( int *ixdrid, u_int *len, char *s, int *dummy, int *retval )

#elif (IFSTYLE==1)
	CXDRSTRING( int *ixdrid, u_int *len, char *s, int *retval )

#elif (IFSTYLE==2)
	cxdrstring( ixdrid, len, s, retval )
	int     *ixdrid;
	u_int	*len;
	char	*s;
	int	*retval;

#else
	cxdrstring_( int *ixdrid, u_int *len, char *s, int *retval )
#endif
{
	XDR 	*xdrs;

	*retval	= 0;	/* No error */

	xdrs = xdrfile[*ixdrid].xdrs;

	if( ! xdr_bytes( xdrs, &s, len, *len ) ) {
		if( xdrs->x_op == XDR_ENCODE ) {
			if( xdrfile[*ixdrid].return_on_error ) {
				*retval = FXDRERR_WRITEERR;
				return;
				}
			else
				{
				fprintf( stderr, "FXDR library error while trying to write a");
				fprintf( stderr, "string from file %s\n", 
					xdrfile[*ixdrid].filename );
				exit( -1 );
				}
			}
		else
			{
			if( xdrfile[*ixdrid].return_on_error ) {
				*retval = FXDRERR_READERR;
				return;
				}
			else
				{
				fprintf( stderr, "FXDR library error while trying to read a");
				fprintf( stderr, "string from file %s\n", 
					xdrfile[*ixdrid].filename );
				exit( -1 );
				}
			}
		}
}

