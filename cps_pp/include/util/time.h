#ifndef UTIL_TIME_H
#define UTIL_TIME_H
#include <config.h>
#include <util/data_types.h>
#include <sys/time.h>
CPS_START_NAMESPACE
Float dclock(void);
Float print_flops(int nflops, Float time);
Float print_flops(int nflops, struct timeval *start, struct timeval *end);
CPS_END_NAMESPACE
#endif
