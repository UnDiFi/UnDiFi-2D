#include <stdio.h>
#include <rpc/rpc.h>
#include "cfxdr.h"

extern XDR_element xdrfile[MAX_N_XDR_FILES];

	void
#if (IFSTYLE==1)
	CXDRINT( int *ixdrid, int *i, int *retval )

#elif (IFSTYLE==2)
	cxdrint( ixdrid, i, retval )
	int    *ixdrid;
	int	*i;
	int	*retval;

#else
	cxdrint_( int *ixdrid, int *i, int *retval )
#endif
{
	XDR 	*xdrs;

	*retval = 0;	/* No error */

	xdrs = xdrfile[*ixdrid].xdrs;

	if( ! xdr_int( xdrs, i ) ) {
		if( xdrs->x_op == XDR_ENCODE ) {
			if( xdrfile[*ixdrid].return_on_error ) {
				*retval = FXDRERR_WRITEERR;
				return;
				}
			else
				{
				fprintf( stderr, "FXDR library error!  Call to write integer ");
				fprintf( stderr, "did not complete successfully!\n" );
				fprintf( stderr, "Error occured in file %s\n", 
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
				fprintf( stderr, "FXDR library error!  Call to read integer ");
				fprintf( stderr, "did not complete successfully!\n" );
				fprintf( stderr, "Error occured in file %s\n", 
					xdrfile[*ixdrid].filename );
				exit( -1 );
				}
			}
		}
}

