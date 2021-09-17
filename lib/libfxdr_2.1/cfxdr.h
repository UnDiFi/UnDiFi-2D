
#define 	MAX_N_XDR_FILES	32

typedef struct {
	char	*filename;
	int 	return_on_error;	/* 0=FALSE, 1=TRUE */
	XDR	*xdrs;
} XDR_element;

#define FXDRERR_WRNEGNELS	-10
#define FXDRERR_WRITEERR	-11
#define FXDRERR_READERR		-12
#define FXDRERR_READWRONGNELS	-13
#define FXDRERR_REWIND		-14

