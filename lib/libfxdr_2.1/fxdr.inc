	integer initxdr
	integer ixdrclose
	integer ixdrdmat
	integer ixdrdouble
	integer ixdrimat
	integer ixdrint
	integer ixdrreal
	integer ixdrreal64
	integer ixdrrewind
	integer ixdrrmat
	integer ixdrrmat64
	integer ixdrshort
	integer ixdrstring

	integer FXDRERR_WRNEGNELS, FXDRERR_WRITEERR
	integer FXDRERR_READERR, FXDRERR_READWRONGNELS
	integer FXDRERR_REWIND
	parameter(FXDRERR_WRNEGNELS=-10, FXDRERR_WRITEERR=-11)
	parameter(FXDRERR_READERR=-12, FXDRERR_READWRONGNELS=-13)
	parameter(FXDRERR_REWIND=-14)

