#include <stdio.h>
#include <rpc/rpc.h>
#include "cfxdr.h"

extern XDR_element xdrfile[MAX_N_XDR_FILES];

	void
#if (IFSTYLE==1)
	CXDRREAL( int *ixdrid, float *r, int *retval )

#elif (IFSTYLE==2)
	cxdrreal( ixdrid, r, retval )
	int    	*ixdrid;
	float	*r;
	int	*retval;

#else
	cxdrreal_( int *ixdrid, float *r, int *retval )
#endif
{
	XDR 	*xdrs;

	*retval = 0;	/* No error */

	xdrs = xdrfile[*ixdrid].xdrs;

	if( ! xdr_float( xdrs, r ) ) {
		if( xdrs->x_op == XDR_ENCODE ) {
			if( xdrfile[*ixdrid].return_on_error ) {
				*retval = FXDRERR_WRITEERR;
				return;
				}
			else
				{
				fprintf( stderr, "FXDR library error!  Call to write real ");
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
				fprintf( stderr, "FXDR library error!  Call to read real ");
				fprintf( stderr, "did not complete successfully!\n" );
				fprintf( stderr, "Error occured in file %s\n", 
					xdrfile[*ixdrid].filename );
				exit( -1 );
				}
			}
		}
}

