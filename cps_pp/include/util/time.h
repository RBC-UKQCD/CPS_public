/*!\file
  \brief Declaration of functions for timing and performance measurement.

  $Id: time.h,v 1.5 2004-09-02 16:58:15 zs Exp $
*/

#ifndef UTIL_TIME_H
#define UTIL_TIME_H

#include <config.h>
#include <util/data_types.h>
#include <sys/time.h>

CPS_START_NAMESPACE
/*! \defgroup profiling Timing and performance functions
  @{
*/

//! Gets the wall-clock time.
Float dclock(void);
//! Prints the FLOPS rate to stdout
Float print_flops(int nflops, Float time);
//! Prints the FLOPS rate to stdout
Float print_flops(char *cname, char *fname,int nflops, Float time);
//! Prints the FLOPS rate to stdout
Float print_flops(int nflops, struct timeval *start, struct timeval *end);
//! Prints the FLOPS rate to stdout
Float print_flops(char *cname, char *fname, int nflops, struct timeval *start, struct timeval *end);

/*! @} */

CPS_END_NAMESPACE
#endif
