#include <config.h>
#include <util/data_types.h>
#include <util/gjp.h>
#include <util/wilson.h>
#define BND_COMM
// use macro version or inlines
#define USE_MACROS
#define __RESTRICT
// The sizes of the block are assumed to be even numbers
#define Z_BLOCK 16
#define Y_BLOCK 16
#define X_DIR_ON
// 2.95(inline), 3.02 (no-inline)
#define Y_DIR_ON
// 3.78(inline), 4.04 (no-inline)
#define Z_DIR_ON
// 3.17(inline), 3.24 (no-inline) 
#define T_DIR_ON
// 3.44(inline), 4.70 (no-inline)

#define POS_DIR_ON
#define NEG_DIR_ON

CPS_START_NAMESPACE
/* wilson dlash kernel pieces */
#define  SSE_C_FLOAT float
//#define  SSE_C_FLOAT IFloat
#define SSE_ALIGNED 
//__attribute__ ((aligned (64)))
#ifdef SSE_TO_C
typedef struct{
SSE_C_FLOAT d[2]  SSE_ALIGNED ;
} M128D;
#else
#define M128D __m128d 
#endif

CPS_END_NAMESPACE
