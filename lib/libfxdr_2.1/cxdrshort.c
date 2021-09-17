/* Write/Read a short */

#include <stdio.h>
#include <rpc/rpc.h>
#include "cfxdr.h"

extern XDR_element xdrfile[MAX_N_XDR_FILES];

	void
#if (IFSTYLE==1)
	CXDRSHORT( int *ixdrid, int *i, int *retval )

#elif (IFSTYLE==2)
	cxdrshort( ixdrid, i, retval )
	int    	*ixdrid;
	int	*i;
	int	*retval;

#else
	cxdrshort_( int *ixdrid, int *i, int *retval )
#endif
{
	XDR 	*xdrs;
	short ii;

	*retval = 0;	/* No error */

	xdrs = xdrfile[*ixdrid].xdrs;

	if (xdrs->x_op == XDR_ENCODE) ii=(short) *i;

	if( ! xdr_short( xdrs, &ii ) ) {
		if( xdrs->x_op == XDR_ENCODE ) {
			if( xdrfile[*ixdrid].return_on_error ) {
				*retval = FXDRERR_WRITEERR;
				return;
				}
			else
				{
				fprintf( stderr, "FXDR library error!  Call to write short ");
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
				fprintf( stderr, "FXDR library error!  Call to read short ");
				fprintf( stderr, "did not complete successfully!\n" );
				fprintf( stderr, "Error occured in file %s\n", 
					xdrfile[*ixdrid].filename );
				exit( -1 );
				}
			}
		}

	if (xdrs->x_op == XDR_DECODE) *i=(int)ii;
}

