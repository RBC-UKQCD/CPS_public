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
 * These are now set by including wfm_options.h
 */

#include "wfm_options.h"

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
/*On a regular paged machine bump up the alloc size to allow linear prefetch off the end without hitting a page boundary*/
#define ALLOC(A) malloc(A+256)
#define SEND_ALLOC(A) malloc(A+256)
#define FREE(A)  free(A)
#endif

#ifdef USE_QALLOC
#include <qalloc.h>
#define ALLOC(A) qalloc(QCOMMS|QFAST,A)
#define SEND_ALLOC(A) qalloc(QNONCACHE|QFAST,A)
#define FREE(A)  qfree(A)

#define TABLE_ALLOC(A) qalloc(0,A)
#define TABLE_FREE(A)  qfree(A)

#else

#define TABLE_ALLOC(A) ALLOC(A)
#define TABLE_FREE(A)  FREE(A)

#endif

#ifdef USE_MEMALIGN
#include<stdlib.h>
#ifndef WFM_ALIGN_ARG
#define WFM_ALIGN_ARG (128)
#endif

extern "C" { 
   void* memalign(size_t, size_t);
};
#define ALLOC(A) memalign(WFM_ALIGN_ARG,A+256)
#define SEND_ALLOC(A) memalign(WFM_ALIGN_ARG,A+256)
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
#ifdef USE_COMMS_QMP
#include <qmp.h>
#endif

#endif
