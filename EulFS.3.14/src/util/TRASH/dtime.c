/*
  From David Klein (klein@sunrise.huji.ac.za)
*/

#include <sys/times.h>
#include <unistd.h>

#include "f2c.h"

doublereal dtime_(tarray)
real *tarray;
{
  struct tms buf;
  time_t ticks;

    /* System generated locals */
    real ret_val;

    /* Parameter adjustments */
    --tarray;

    /* Function Body */
  ticks=times(&buf);
  tarray[1]=buf.tms_utime/100.0;
  ret_val=tarray[1];
  return ret_val;
}
