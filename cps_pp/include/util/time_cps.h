/*!\file
  \brief Declaration of functions for timing and performance measurement.

  $Id: time_cps.h,v 1.7 2013-04-05 17:46:30 chulwoo Exp $
*/

#ifndef UTIL_TIME_H
#define UTIL_TIME_H

#include <config.h>
#include <util/data_types.h>
#include <sys/time.h>
#include <time.h>

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


static double time_now, time_prev;
inline double time_elapse(){
  time_now = CPS_NAMESPACE::dclock();
  double elp  = time_now - time_prev;
  time_prev=time_now;
  return elp;
}

void print_asctime_();

#define print_asctime(msg,a ...) do {			\
    if(! UniqueID()){					\
      printf("asctime[%05d] " msg, UniqueID() ,##a); 	\
      print_asctime_(); }}				\
  while(0);



/*! @} */

CPS_END_NAMESPACE
#endif
