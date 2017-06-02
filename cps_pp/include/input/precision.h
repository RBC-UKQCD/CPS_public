//#include <stdint.h>
#ifndef _PRECISION_H_
#define _PRECISION_H_
typedef double Float;
typedef double IFloat;
//typedef uint64_t IntFlopCounter;
#ifdef HAVE_UINT64_T
typedef uint64_t Pointer;
#else
typedef uint32_t Pointer;
#endif
#endif


