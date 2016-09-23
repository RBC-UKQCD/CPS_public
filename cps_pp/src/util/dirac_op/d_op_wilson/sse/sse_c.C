#if (defined USE_SSE)||(defined SSE_TO_C)
#include <config.h>
#include <util/data_types.h>
#include <util/gjp.h>
#include <util/wilson.h>
#include "sse-defs.h"
#define BND_COMM
//#define USE_MACROS
#define __RESTRICT
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


//#define  SSE_C_FLOAT IFloat

CPS_START_NAMESPACE

#include "sse-subs.h"

#include "sse-blk-dag0.h"
#include "sse-blk-dag1.h"

#include "sse-bnd-dag0.h"
#include "sse-bnd-dag1.h"

CPS_END_NAMESPACE
#endif
