#include "config.h"
#define PRECISION 8
#if ( SIZEOF_VOID_P == SIZEOF_LONG )
typedef long asmint;
#else 
#if ( SIZEOF_VOID_P == SIZEOF_INT )
typedef int asmint;
#else
#error Dont know how to represent integers
#endif
#endif


#if ( SIZEOF_DOUBLE == PRECISION )
typedef double fpoint;
#else
#if ( SIZEOF_FLOAT == PRECISION )
typedef float fpoint;
#else
#error Cannot represent PRECISION 
#endif
#endif

