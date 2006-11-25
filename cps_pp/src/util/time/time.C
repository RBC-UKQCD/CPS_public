/*!\file
  \brief Implementation of functions for timing and performance measurement.

  $Id: time.C,v 1.9 2006-11-25 19:10:50 chulwoo Exp $
*/

#include <config.h>
#include <sys/time.h>
#include <util/qcdio.h>
#include <util/data_types.h>

CPS_START_NAMESPACE

/*!
  \return The current time in seconds (accurate to the microsecond),
*/
Float dclock(void){
	struct timeval tp;
	gettimeofday(&tp, NULL);
	return((Float) tp.tv_sec + (Float) tp.tv_usec * 1e-6);
}

/*!
  The user counts the operations and calls \a gettimeofday before
  they start and after they finish to get the elapsed time.
  \param n_flops The number of floating point operations performed
  \param start The initial result of calling \a gettimeofday
  \param end The final result of calling \a gettimeofday  
  \return The FLOPS rate.
*/

Float print_flops(unsigned long long nflops, struct timeval *start, struct timeval *end){
	int sec = end->tv_sec - start->tv_sec; 
	int usec = end->tv_usec - start->tv_usec; 
	Float time = sec + 1.e-6*usec;
	if(!UniqueID())
	printf("%e flops /%e seconds = %e MFlops\n",(Float)nflops,time,(Float)nflops/(time*1.e6));
	return nflops/time;
}

/*!
  The user counts the operations and calls \a gettimeofday before
  they start and after they finish to get the elapsed time.
  The class name and method name are printed with the FLOPS rate.  
  \param cname The class name.
  \param fname The method name.
  \param n_flops The number of floating point operations performed
  \param start The initial result of calling \a gettimeofday
  \param end The final result of calling \a gettimeofday  
  \return The FLOPS rate.
*/
Float print_flops(char *cname, char *fname, unsigned long long nflops, struct timeval *start, struct timeval *end){
  if (!UniqueID()){
	printf("%s:%s: ",cname,fname);
  }
	return print_flops(nflops,start,end);
}

/*!
  The user counts the operations and calls dclock before they start and after
  they finish  to get the elapsed time.
  \param n_flops The number of floating point operations performed
  \param time The elapsed time in seconds.
  \return The FLOPS rate.
*/
Float print_time(const char *cname, const char *fname, Float time){
  if (!UniqueID())
	printf("%s::%s: %e seconds\n",cname,fname,time);
	return time;
}

Float print_flops(unsigned long long nflops, Float time){
	if(!UniqueID())
	printf("%e flops /%e seconds = %e MFlops\n",(Float)nflops,time,(Float)nflops/(time*1.e6));
	return nflops/time;
}

/*!
  The user counts the operations and calls dclock before they start and
  after they finish to get the elapsed time.
  The class name and method name are printed with the FLOPS rate.  
  \param cname The class name.
  \param fname The method name.
  \param n_flops The number of floating point operations performed.
  \param time The elapsed time in seconds.
  \return The FLOPS rate.
*/

Float print_flops(char *cname, char *fname, unsigned long long nflops, Float time){
	if(!UniqueID())
	printf("%s::%s: ",cname,fname);
	return print_flops(nflops,time);
}

CPS_END_NAMESPACE
