/*!\file
  \brief Declaration of functions for timing and performance measurement.

*/

#ifndef UTIL_TIME_H
#define UTIL_TIME_H

#include <config.h>
#include <util/data_types.h>
#include <comms/sysfunc_cps.h>
#include <sys/time.h>
#include <time.h>
#include <string>

CPS_START_NAMESPACE
/*! \defgroup profiling Timing and performance functions
  @{
*/


//! Gets the wall-clock time.

Float dclock(bool stop = false);
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

//CK: static timestamp class with optional output stream and 'depth' tagging of stamps enabling easy parsing of output
//    writes output to a file

class TimeStamp{
 protected:
  static int cur_depth;
  static FILE *stream;
  static bool enabled;
  static Float start;
  
  static void vstamp(const char *format, va_list args);
 public:
  static void set_file(const char *filename);
  static void reset();
  static void close_file();
  static void stamp(const char *format,...);
  static void incr_depth();
  static void decr_depth();
  
  static void stamp_incr(const char *format,...);
  static void stamp_decr(const char *format,...);

  static void incr_stamp(const char *format,...);
  static void decr_stamp(const char *format,...);

  static void start_func(const char* cls, const char *fnc);
  static void end_func(const char* cls, const char *fnc);
};

struct Elapsed{
  int hours;
  int mins;
  Float secs;

  Elapsed():hours(0),mins(0),secs(0){}
  Elapsed(const Float &delta);
  inline void print(const std::string &stamp = "Elapsed", FILE *to = stdout){
    if(!UniqueID()){
      fprintf(stdout,"%s %dh %dm %fs\n",stamp.c_str(),hours,mins,secs); fflush(stdout);
    }
  }
};

// Moving Timer to util/timer.h and merging with the existing version 
#if 0
class Timer{
 protected:
  static Float dtime_begin;
  static Float dtime_last;

 public:
  static void reset();
  //Time since last reset
  static Elapsed elapsed_time();
  //Time since last call to this function
  static Elapsed relative_time();
};
#endif
CPS_END_NAMESPACE
#endif
