#include <stdio.h>
#include <rpc/rpc.h>
#include "cfxdr.h"

extern XDR_element xdrfile[MAX_N_XDR_FILES];

	void
#if (IFSTYLE==1)
	CXDRDMAT( int *ixdrid, int *nels, double *d, int *retval )

#elif (IFSTYLE==2)
	cxdrdmat( ixdrid, nels, d, retval )
	int     *ixdrid;
	int	*nels;
	double	*d;
	int	*retval;

#else
	cxdrdmat_( int *ixdrid, int *nels, double *d, int *retval )
#endif
{
	XDR 	*xdrs;
	u_int	actual_nels, ne;

	*retval = 0;	/* No error */

	xdrs = xdrfile[*ixdrid].xdrs;

	if( xdrs->x_op == XDR_ENCODE ) {
		if( *nels < 0 ) {
			if( xdrfile[*ixdrid].return_on_error ) {
				*retval = FXDRERR_WRNEGNELS;
				return;
				}
			else
				{
				fprintf( stderr, "FXDR library error in call to ");
				fprintf( stderr, "write double precision array to file %s:\n",
							xdrfile[*ixdrid].filename );
				fprintf( stderr, "negative number of elements specified!\n" );
				exit( -1 );
				}
			}
		ne = *nels;
		if( ! xdr_array( xdrs, (char **)&d, &ne, *nels, sizeof(double), xdr_double ) ) {
			if( xdrfile[*ixdrid].return_on_error ) {
				*retval = FXDRERR_WRITEERR;
				return;
				}
			else
				{
				fprintf( stderr, "FXDR library error in call to " );
				fprintf( stderr, "write double precision array to file %s:\n",
							xdrfile[*ixdrid].filename );
				fprintf( stderr, "write error (disk full?)\n" );
				exit( -1 );
				}
			}
		}
	else
		{
		if( ! xdr_array( xdrs, (char **)&d, &actual_nels, *nels, sizeof(double), xdr_double )) {
			if( xdrfile[*ixdrid].return_on_error ) {
				*retval = FXDRERR_READERR;
				return;
				}
			else
				{
				fprintf( stderr, "FXDR library error in call to " );
				fprintf( stderr, "read double precision array from file %s:\n",
							xdrfile[*ixdrid].filename );
				fprintf( stderr, "error reading the file\n" );
				exit( -1 );
				}
			}

		if( actual_nels != *nels ) {
			if( xdrfile[*ixdrid].return_on_error ) {
				*retval = FXDRERR_READWRONGNELS;
				return;
				}
			else
				{
				fprintf( stderr, "FXDR library error!  Asked for " );
				fprintf( stderr, "%d double precision elements, but actually read %d!\n",
						*nels, actual_nels );
				fprintf( stderr, "Error occured in file %s\n",
					xdrfile[*ixdrid].filename );
				exit( -1 );
				}
			}
		} 
}
