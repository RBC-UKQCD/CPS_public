#include <config.h>
#include <sys/time.h>
#include <stdio.h>
#include <util/data_types.h>
CPS_START_NAMESPACE
Float dclock(void){
	struct timeval tp;
	gettimeofday(&tp, NULL);
	return((Float) tp.tv_sec + (Float) tp.tv_usec * 1e-6);
}

Float print_flops(int nflops, struct timeval *start, struct timeval *end){
	int sec = end->tv_sec - start->tv_sec; 
	int usec = end->tv_usec - start->tv_usec; 
	Float time = sec + 1.e-6*usec;
	printf("%e flops /%e seconds = %e MFlops\n",(Float)nflops,time,(Float)nflops/(time*1.e6));
}

Float print_flops(int nflops, Float time){
	printf("%e flops /%e seconds = %e MFlops\n",(Float)nflops,time,(Float)nflops/(time*1.e6));
}

CPS_END_NAMESPACE
