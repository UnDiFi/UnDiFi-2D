#include <stdio.h>
#include <rpc/rpc.h>
#include "cfxdr.h"

extern XDR_element xdrfile[MAX_N_XDR_FILES];

	void
#if (IFSTYLE==1)
	CXDRDOUBLE( int *ixdrid, double *d, int *retval )

#elif (IFSTYLE==2)
	cxdrdouble( ixdrid, d, retval )
	int    	*ixdrid;
	double	*d;
	int	*retval;

#else
	cxdrdouble_( int *ixdrid, double *d, int *retval )
#endif
{
	XDR 	*xdrs;

	*retval = 0;	/* No error */

	xdrs = xdrfile[*ixdrid].xdrs;

	if( ! xdr_double( xdrs, d )) {
		if( xdrs->x_op == XDR_ENCODE ) {
			if( xdrfile[*ixdrid].return_on_error ) {
				*retval = FXDRERR_WRITEERR;
				return;
				}
			else
				{
				fprintf( stderr, "FXDR library error while trying to write " );
				fprintf( stderr, "a double precision number from file %s\n", 
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
				fprintf( stderr, "FXDR library error while trying to read " );
				fprintf( stderr, "a double precision number to file %s\n", 
						xdrfile[*ixdrid].filename );
				exit( -1 );
				}
			}
		}
}
