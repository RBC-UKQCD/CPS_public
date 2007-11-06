#ifndef _WFM_PAB_CONFIG_H_
#define _WFM_PAB_CONFIG_H_

/*---------------------------------------------------------------------------
 * Configuration of the library.
 *---------------------------------------------------------------------------
 */

/*
 * Precision control of the C++ library.
 * Make system must couple this to linking the appropriate assembler kernels
 *
 * Under the control of the preprocessor flags:
 *
 *precision:
 * USE_DOUBLE
 *
 *allocation:
 * USE_QALLOC
 * USE_MEMALIGN
 * USE_POSIX_MEMALIGN
 * USE_MALLOC
 * USE_GMALLOC
 * USE_QMPALLOC
 *
 *comms:
 * USE_COMMS_NONE
 * USE_COMMS_FAKE
 * USE_COMMS_SCU
 * USE_COMMS_QMP
 * USE_COMMS_MPI
 * USE_COMMS_GM
 *
 * These should be mutually exclusive (i.e. only one set) 
 * and have the subsidiary preprocessor flag:
 *
 * WFM_ALIGN_ARG
 *
 * Set to the appropriate alignment value. (defaults 128 bytes in this header)
 *
 */
#define USE_DOUBLE
#if TARGET == QCDOC
#define USE_CPS
#define USE_COMMS_SCU
#else
#undef USE_CPS
#define USE_COMMS_QMP
#define USE_MALLOC
#endif

#ifdef USE_DOUBLE
typedef double Float ;
#else 
typedef float Float ;
#endif

/*
 * Heap allocation control.
 *
 */


#ifdef USE_MALLOC
#define ALLOC(A) malloc(A)
#define SEND_ALLOC(A) malloc(A)
#define FREE(A)  free(A)
#endif

#ifdef USE_CPS
#include <qalloc.h>
#include <util/gjp.h>
#define ALLOC(A) qalloc((GJP.WfmAllocFlag()),A)
#define SEND_ALLOC(A) qalloc((GJP.WfmSendAllocFlag()),A)
#define FREE(A)  qfree(A)
#endif

#ifdef USE_QALLOC
#include <qalloc.h>
#define ALLOC(A) qalloc(QCOMMS|QFAST,A)
#define SEND_ALLOC(A) qalloc(QNONCACHE|QFAST,A)
#define FREE(A)  qfree(A)
#endif

#ifdef USE_MEMALIGN

#ifndef WFM_ALIGN_ARG
#define WFM_ALIGN_ARG (128)
#endif

#define ALLOC(A) memalign(WFM_ALIGN_ARG,A)
#define SEND_ALLOC(A) memalign(WFM_ALIGN_ARG,A)
#define FREE(A)  free(A)

#endif

#ifdef USE_POSIX_MEMALIGN
#ifndef WFM_ALIGN_ARG
#define WFM_ALIGN_ARG (128)
#endif
#define ALLOC(A) posix_memalign(WFM_ALIGN_ARG,A)
#define SEND_ALLOC(A) posix_memalign(WFM_ALIGN_ARG,A)
#define FREE(A)  free(A)
#endif

#ifdef USE_GMALLOC
#error No GMALLOC Support
#endif

#ifdef USE_QMPALLOC
#error No QMPALLOC Support
#endif

/*
 * Communications library control
 */

#endif
