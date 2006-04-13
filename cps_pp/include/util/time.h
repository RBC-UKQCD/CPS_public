/*!\file
  \brief Declaration of functions for timing and performance measurement.

  $Id: time.h,v 1.7 2006-04-13 19:12:20 chulwoo Exp $
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
Float print_time(char *cname, char *fname, Float time);

//! Prints the FLOPS rate to stdout
Float print_flops(unsigned long long nflops, Float time);
//! Prints the FLOPS rate to stdout
Float print_flops(char *cname, char *fname,unsigned long long nflops, Float time);
//! Prints the FLOPS rate to stdout
Float print_flops(unsigned long long nflops, struct timeval *start, struct timeval *end);
//! Prints the FLOPS rate to stdout
Float print_flops(char *cname, char *fname, unsigned long long nflops, struct timeval *start, struct timeval *end);

/*! @} */

CPS_END_NAMESPACE
#endif
