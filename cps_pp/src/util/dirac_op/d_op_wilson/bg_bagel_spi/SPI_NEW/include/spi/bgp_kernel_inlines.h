/* begin_generated_IBM_copyright_prolog                             */
/*                                                                  */
/* This is an automatically generated copyright prolog.             */
/* After initializing,  DO NOT MODIFY OR MOVE                       */
/*  --------------------------------------------------------------- */
/*                                                                  */
/* (C) Copyright IBM Corp.  2007, 2007                              */
/* IBM CPL License                                                  */
/*                                                                  */
/*  --------------------------------------------------------------- */
/*                                                                  */
/* end_generated_IBM_copyright_prolog                               */
#ifndef _BGP_SPI_KERNEL_INLINES_H_ // Prevent multiple inclusion
#define _BGP_SPI_KERNEL_INLINES_H_

#include <common/namespace.h>

__BEGIN_DECLS
#include <bpcore/ppc450_inlines.h>
#include <cnk/bgp_SysCall_Extensions.h>
#include <cnk/bgp_SPRG_Usage.h>

#ifndef __INLINE__
#define __INLINE__ extern inline
#endif

/*! \brief LockBox allocate syscall definition
 * \param[in] lockid Indicates which counter ID is to be obtained.  Counter IDs vary from 0-1023
 * \param[in] numlocks The number of sequencial counter IDs that will be obtained
 * \param[out] ptr An array of pointers that will be filled in with the counter virtual addresses.
 * \param[in] flags Optional flags
 * \warning Must storage indicated by ptr must be large enough to whole numlocks*sizeof(uint32_t) bytes
 * \internal This is an internal syscall - do not use.
 * \see LockBox_AllocateCounter 
 * \see LockBox_AllocateMutex
 * \see LockBox_AllocateBarrier
 */
__INLINE__ int Kernel_AllocateLockBox(uint32_t lockid, uint32_t numlocks, uint32_t** ptr, uint32_t flags)
{
  int rc = 0; // this syscall returns RC in r3 and does not use errno
 
  asm __volatile__ ("li 0,%1;"
                    "mr 3,%2;"
                    "mr 4,%3;"
                    "mr 5,%4;"
                    "mr 6,%5;"
                     "sc;"
                    "mr %0, 3;"
                    : "=&r" (rc)  // early clobber
                    : "i" (_BGP_SYSCALL_NR_ALLOC_COUNTER),
                    "r" (lockid),
                    "r" (numlocks),
                    "r" (ptr),
                    "r" (flags)
                    : "r0", "r3", "r4", "r5", "r6", "cc", "memory" );

  return( rc );
}


/*! \brief Returns the physical processor ID of the running PPC450 core.
 *
 * \return Physical processor ID
 * \retval 0 Running on processor 0
 * \retval 1 Running on processor 1
 * \retval 2 Running on processor 2
 * \retval 3 Running on processor 3
 */
__INLINE__ int Kernel_PhysicalProcessorID( void )
{
  uint32_t dst2  = _bgp_mfspr( _BGP_SPRGRO_DST2  );

  return( dst2 & 0x3 );
}

__END_DECLS
  
  
  
#endif // Add nothing below this line

