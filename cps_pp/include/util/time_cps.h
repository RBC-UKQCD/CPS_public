/*!\file
  \brief Declaration of functions for timing and performance measurement.

  $Id: time_cps.h,v 1.6 2012-03-27 05:02:40 chulwoo Exp $
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
Float print_time(const char *cname, const char *fname, Float time);

//! Prints the FLOPS rate to stdout
Float print_flops(double nflops, Float time);
//! Prints the FLOPS rate to stdout
//Float print_flops(const char cname[], const char fname[],unsigned long long nflops, Float time);
Float print_flops(const char cname[], const char fname[],double nflops, Float time);
//! Prints the FLOPS rate to stdout
Float print_flops(double nflops, struct timeval *start, struct timeval *end);
//! Prints the FLOPS rate to stdout
//Float print_flops(const char cname[], const char fname[], unsigned long long nflops, struct timeval *start, struct timeval *end);
Float print_flops(const char cname[], const char fname[], double nflops, struct timeval *start, struct timeval *end);

/*! @} */

CPS_END_NAMESPACE
#endif
