#ifndef _INLINE_SSE_H_
#define _INLINE_SSE_H_

/* These assembly inline versions work only with the gnu C compiler on
   P3 or P4 processors */

#include "qla_types.h"

#include "inline_headers.h"

typedef struct
{
   unsigned int c1,c2,c3,c4;
} sse_mask __attribute__ ((aligned (16)));

static sse_mask _sse_sgn13 __attribute__ ((unused)) ={0x80000000, 0x00000000, 0x80000000, 0x00000000};
static sse_mask _sse_sgn24 __attribute__ ((unused)) ={0x00000000, 0x80000000, 0x00000000, 0x80000000};
static sse_mask _sse_sgn3 __attribute__  ((unused)) ={0x00000000, 0x00000000, 0x80000000, 0x00000000};
static sse_mask _sse_sgn4 __attribute__  ((unused)) ={0x00000000, 0x00000000, 0x00000000, 0x80000000};
static sse_mask _sse_sgn2 __attribute__  ((unused)) ={0x00000000, 0x80000000, 0x00000000, 0x00000000};

#endif
