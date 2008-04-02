/******************************************************************************/
/*                                                                            */
/*                            SMP & VN MODE BARRIERS                          */
/*                                                         30.06.07 S.Krieg   */
/******************************************************************************/

#ifndef _NODEBARRIER_H_
#define _NODEBARRIER_H_

#ifdef SMP_MODE
#ifdef VN_MODE
#error Both SMP_MODE and VN_MODE specified.
#endif
#endif

#ifdef SMP_MODE
#include "SMP_NodeBarrier.h"
extern inline void NodeBarrier(void) { 
    SMP_NodeBarrier();
}
#define NODE_BARRIER SMP_NodeBarrier()
#else
#ifdef VN_MODE
#include "VN_NodeBarrier.h"
extern inline void NodeBarrier(void) {
    VN_NodeBarrier();
}
#define NODE_BARRIER VN_NodeBarrier()
#else
#error Barrier not using VN or SMP mode, but nothing else (incl DUAL mode) is defined yet.
#endif
#endif

/* for compatibility */
#define CHIP_BARRIER NODE_BARRIER

#endif
