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
/**
 * \file spi/lockbox_interface.h
 */

#ifndef	_BGP_SPI_LOCKBOX_INLINES_H_ // Prevent multiple inclusion
#define	_BGP_SPI_LOCKBOX_INLINES_H_

#include <common/namespace.h>

__BEGIN_DECLS

#include <common/linkage.h>
#include <common/bgp_bitnumbers.h>
#include <bpcore/bgp_types.h>
#include <bpcore/ppc450_core.h>
#include <bpcore/ppc450_inlines.h>
#include <spi/bgp_kernel_inlines.h>

#include <errno.h>
// op codes encoded in lockbox address bits 21:22
#define _BGP_LOCKBOX_OP_QUERY  _B2(22,0x0)        // addr[21:22] = 0b00 (RW)
#define _BGP_LOCKBOX_OP_INC    _B2(22,0x1)        // addr[21:22] = 0b01 (RO)
#define _BGP_LOCKBOX_OP_DEC    _B2(22,0x2)        // addr[21:22] = 0b10 (RO)
#define _BGP_LOCKBOX_OP_CLEAR  _B2(22,0x3)        // addr[21:22] = 0b11 (RO)

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
//__INLINE__ int Kernel_AllocateLockBox(uint32_t lockid, uint32_t numlocks, uint32_t** ptr, uint32_t flags)
//{
//  int rc = 0; // this syscall returns RC in r3 and does not use errno

//  asm __volatile__ ("li 0,%1;"
//                    "mr 3,%2;"
//                    "mr 4,%3;"
//                    "mr 5,%4;"
//                    "mr 6,%5;"
//                     "sc;"
//                    "mr %0, 3;"
//                    : "=&r" (rc)  // early clobber
//                    : "i" (_BGP_SYSCALL_NR_ALLOC_COUNTER),
//                    "r" (lockid),
//                    "r" (numlocks),
//                    "r" (ptr),
//                    "r" (flags)
//                    : "r0", "r3", "r4", "r5", "r6", "cc", "memory" );

//  return( rc );
//}

/*! \brief Assures that lock allocate is ordered across all CPUs partipating in the allocation */
#define LOCKBOX_ORDERED_ALLOC 0x100

typedef uint32_t* LockBox_Counter_t;  /*!< Counter ID definition */

/*! \brief Structure contains information required to perform a lockbox mutex.
 * This structure will be filled out by the LockBox_AllocateMutex() function
 * \see LockBox_AllocateMutex
 */
typedef LockBox_Counter_t LockBox_Mutex_t;

/*! \brief Structure contains information required to perform a lockbox barrier.
 * This structure will be filled out by the LockBox_AllocateBarrier() function
 * \see LockBox_AllocateBarrier
 */
typedef struct
{
  LockBox_Counter_t ctrl_lock;  /*!< Counter controls which lock phase we're using */
  LockBox_Counter_t lock[2];    /*!< Two counter locks to perform barrier synchronization */
  LockBox_Counter_t status[2];  /*!< Two counter locks to accumulate status */
  uint8_t mastercore;           /*!< Described which PPC450 is controlling the barrier synchronization */
  uint8_t numcores;             /*!< Number of cores participating in the barrier */
} LockBox_Barrier_t;

/*! \brief Queries and returns the value of the lockbox counter
 *
 * Lockbox counters can range from 0x00000000 to 0xffffffff
 * Counter values are unsigned.
 * \pre Counter must be allocated via the LockBox_AllocateCounter() function.
 *
 * \param[in] lba lockbox counter
 *
 */
// Return the current value of a lock without changing lock contents
__INLINE__ uint32_t LockBox_Query(LockBox_Counter_t lba)
{
  return (uint32_t)_bgp_In32((uint32_t*)((char *)lba + _BGP_LOCKBOX_OP_QUERY));
}

/*! \brief Sets the value of the counter.
 *
 * Value may range from 0x00000000 to 0xffffffff (unsigned)
 * \pre Counter must have been allocated via the LockBox_AllocateCounter() routine.
 * \warning LockBox_Write is not an atomic operation
 *
 * \param[in] lba   lockbox counter
 * \param[in] value new value for the lockbox
 *
 */
__INLINE__ void  LockBox_Write(LockBox_Counter_t lba, uint32_t value )
{
   _bgp_Out32((uint32_t*)((uint32_t *)((char*)lba + _BGP_LOCKBOX_OP_QUERY)), value);
}

/*! \brief Fetches and increments the specified counter
 *
 * Fetches the current value of the counter and returns it.  Meanwhile, the
 * fetch has a side effect of atomically incrementing the counter value.
 *
 * Value ranges from 0x00000000 to 0xffffffff (unsigned)
 * Counter must have been allocated via LockBox_AllocateCounter() function.
 *
 * \param[in] lba lockbox counter
 *
 */
__INLINE__ uint32_t LockBox_FetchAndInc(LockBox_Counter_t lba)
{
   return (uint32_t)_bgp_In32((uint32_t*)((char*)lba + _BGP_LOCKBOX_OP_INC));
}

/*! \brief Fetches and decrements the specified counter
 *
 * Fetches the current value of the counter and returns it.  Meanwhile, the
 * fetch has a side effect of atomically decrementing the counter value.
 *
 * Value ranges from 0x00000000 to 0xffffffff (unsigned)
 * Counter must have been allocated via LockBox_AllocateCounter() function.
 *
 * \param[in] lba lockbox counter
 *
 */
__INLINE__ uint32_t LockBox_FetchAndDec(LockBox_Counter_t lba)
{
   return (uint32_t)_bgp_In32( (uint32_t*)((char*)lba + _BGP_LOCKBOX_OP_DEC) );
}

/*! \brief Fetches and Clear the specified counter
 *
 * Fetches the current value of the counter and returns it.  Meanwhile, the
 * fetch has a side effect of atomically clearing the counter to 0x00000000.
 *
 * Counter must have been allocated via LockBox_AllocateCounter() function.
 *
 * \param[in] lba lockbox counter
 */
__INLINE__ uint32_t LockBox_FetchAndClear(LockBox_Counter_t lba)
{
  return (uint32_t)_bgp_In32( (uint32_t*)((char*)lba + _BGP_LOCKBOX_OP_CLEAR ));
}

/*! \brief Allocates a Lockbox Counter
 *
 * Allocates a lockbox counter.  Also has the ability to have multiple
 * cores rendevous on a particular counter ID.
 *
 * \param[in]  desired    The desired counter identifier.  Ranges from 0-1023
 * \param[out] ptr        The pointer to storage for the lockbox counter ptr
 * \param[in]  mastercore The PPC450 core ID responsible for coordinating the allocation
 * \param[in]  numcores   The number of PPC450 core's that will be rendevousing on the counter allocation
 * \param[in]  flags      Optional flags to LockBox_AllocateCounter
 *
 * \return Contains the error indication
 * \retval 0  Success
 * \retval -1 An error occurred, check errno
 * \todo Investigate better way of lock rendevous
 */
__INLINE__ int LockBox_AllocateCounter(uint32_t desired, LockBox_Counter_t* ptr, uint32_t mastercore, uint32_t numcores, uint32_t flags)
{
  if(numcores == 0)
    return EINVAL;
  if(numcores > 4)
    return EINVAL;
  return Kernel_AllocateLockBox(desired, 1, ptr, (mastercore<<4) | numcores | flags);
}

/*! \brief Allocates a Lockbox Mutex
 *
 * Allocates a lockbox-based mutex.  Also has the ability to have multiple
 * cores rendevous on a particular counter ID.
 *
 * \param[in]  desired    The desired counter identifier.  Ranges from 0-1023
 * \param[out] ptr        The pointer to storage for the LockBox_Mutex ptr
 * \param[in]  mastercore The PPC450 core ID responsible for coordinating the allocation
 * \param[in]  numcores   The number of PPC450 core's that will be rendevousing on the counter allocation
 * \param[in]  flags      Optional flags to LockBox_AllocateCounter
 *
 * \return Contains the error indication
 * \retval 0  Success
 * \retval -1 An error occurred, check errno
 * \todo Investigate better way of lock rendevous
 */
__INLINE__ int LockBox_AllocateMutex(uint32_t desired, LockBox_Mutex_t* ptr, uint32_t mastercore, uint32_t numcores, uint32_t flags)
{
  if(numcores == 0)
    return EINVAL;
  if(numcores > 4)
    return EINVAL;
  return Kernel_AllocateLockBox(desired, 1, ptr, (mastercore<<4) | numcores | flags);
}

/*! \brief Allocates a Lockbox Barrier
 *
 * Allocates a lockbox-based barrier.  Also has the ability to have multiple
 * cores rendevous on a particular counter ID.
 *
 * Note: Barrier's consume 5 counters
 *
 * \param[in] desired    The desired counter identifier.  Ranges from 0-1023
 * \param[out] ptr       The pointer to storage for the lockbox counter ptr
 * \param[in] mastercore The PPC450 core ID responsible for coordinating the allocation
 * \param[in] numcores   The number of PPC450 core's that will be rendevousing on the counter allocation
 * \param[in] flags      Optional flags to LockBox_AllocateCounter
 *
 * \return Contains the error indication
 * \retval 0  Success
 * \retval -1 An error occurred, check errno
 * \todo Investigate better way of lock rendevous
 */
__INLINE__ int LockBox_AllocateBarrier(uint32_t desired, LockBox_Barrier_t* ptr, uint32_t mastercore, uint32_t numcores, uint32_t flags)
{
  if(numcores == 0)
    return EINVAL;
  if(numcores > 4)
    return EINVAL;
  ptr->mastercore = mastercore;
  ptr->numcores = numcores;
  return Kernel_AllocateLockBox(desired, 5, &ptr->ctrl_lock, (mastercore<<4) | numcores | flags);
}


/*! \brief Locks a lockbox-based Mutex
 *
 * Locks a mutex.  If the mutex is already locked by another core,
 * LockBox_MutexLock will spin indefinitely waiting to acquire the
 * lock.
 *
 * \param[in] ptr Pointer to an allocated lockbox-based mutex
 *
 * \return Contains the error indication
 * \retval 0  Success
 * \retval -1 An error occurred
 */
__INLINE__ int LockBox_MutexLock(LockBox_Mutex_t mutex)
{
  while ( LockBox_FetchAndInc((LockBox_Counter_t)mutex) != 0 ) {}
  
  return 0;
}

/*! \brief Attempts to lock a lockbox-based Mutex
 *
 * Attempts to lock a mutex.  If the mutex is already locked by
 * another core, LockBox_MutexTryLock will return EAGAIN.
 *
 * \param[in] ptr Pointer to an allocated lockbox-based mutex
 *
 * \return Contains the error indication
 * \retval 0     Success
 * \retval EBUSY The lock was already taken
 */
__INLINE__ int LockBox_MutexTryLock(LockBox_Mutex_t mutex)
{
  uint32_t value;

  value = LockBox_FetchAndInc((LockBox_Counter_t)mutex);
  if(value == 0)
    {
      return 0;
    }
  return EBUSY;
}

/*! \brief Unlocks a lockbox-based Mutex
 *
 * Unlocks a locked mutex
 *
 * \param[in] ptr Pointer to an allocated lockbox-based mutex
 *
 * \return Contains the error indication
 * \retval 0  Success
 * \retval -1 An error occurred, check errno
 */
__INLINE__ int LockBox_MutexUnlock(LockBox_Mutex_t mutex)
{
  LockBox_FetchAndClear((LockBox_Counter_t)mutex);
  return 0;
}

/*! \brief Performs a lockbox-based barrier
 *
 * Will spin until the number of PPC450 cores that have entered the barrier is equal to ptr->numcores.
 *
 * It is a requirement that ptr->mastercore participates in the barrier.
 *
 * \param[in] ptr Pointer to an allocated lockbox-based barrier
 *
 * \return Contains the error indication
 * \retval 0  Success
 * \retval -1 An error occurred, check errno
 */
__INLINE__ int LockBox_Barrier(LockBox_Barrier_t* ptr)
{
  uint32_t lockup, value;
  _bgp_msync();
  lockup = LockBox_Query(ptr->ctrl_lock);
  LockBox_FetchAndInc(ptr->lock[lockup]);

  // Wait until all cores in group have checked in.
  // Can't just loop while value < ptr->numcores, because
  // the master could enter the barrier while our non-master
  // is looping here, clear the counter, and we would
  // loop forever waiting for the counter to reach ptr->numcores.
  // Thus, must check for zero as a loop termination condition too.
  do
  {
    value = LockBox_Query( ptr->lock[lockup] );
  } while ( (value > 0) && (value < ptr->numcores) );

  if (Kernel_PhysicalProcessorID() == ptr->mastercore)
    {
      if(lockup)
        LockBox_FetchAndDec(ptr->ctrl_lock);
      else
        LockBox_FetchAndInc(ptr->ctrl_lock);

      LockBox_FetchAndClear(ptr->status[lockup]);
      LockBox_FetchAndClear(ptr->lock[lockup]);
    }
  else
    {
      // wait until master releases the barrier by clearing the lock
      while( LockBox_Query(ptr->lock[lockup]) > 0 )
        {
        }
    }
  return 0;
}

/*! \brief Performs a lockbox-based barrier and accumulates status
 *
 * Will spin until the number of PPC450 cores that have entered the barrier is equal to ptr->numcores.
 *
 * It is a requirement that ptr->mastercore participates in the barrier.
 *
 * \param[in] ptr Pointer to an allocated lockbox-based barrier
 * \param[in] status Status bit that is accumulated with other cores status bits participating in the barrier
 *
 * Accumulated status is returned.
 *
 * \return Contains the error indication
 * \retval 0  Success
 * \retval -1 An error occurred, check errno
 *
 */
__INLINE__ int LockBox_BarrierStatus(LockBox_Barrier_t* ptr, uint32_t status)
{
  uint32_t lockup, value;
  uint32_t status_rc;
  _bgp_msync();

  lockup = LockBox_Query(ptr->ctrl_lock);
  if(status)
    {
      LockBox_FetchAndInc(ptr->status[lockup]);
    }
  LockBox_FetchAndInc(ptr->lock[lockup]);

  // Wait until all cores in group have checked in.
  // Can't just loop while value < ptr->numcores, because
  // the master could enter the barrier while our non-master
  // is looping here, clear the counter, and we would
  // loop forever waiting for the counter to reach ptr->numcores.
  // Thus, must check for zero as a loop termination condition too.
  do
  {
    value = LockBox_Query( ptr->lock[lockup] );
  } while ( (value > 0) && (value < ptr->numcores) );

  // all cores read the accumulated status
  status_rc = LockBox_Query(ptr->status[lockup]);

  // all cores indicate that they have read status by checking in again
  LockBox_FetchAndInc(ptr->lock[lockup]);

  // wait until all cores in group have checked in for the second time
  while( LockBox_Query(ptr->lock[lockup]) < (unsigned)(2 * ptr->numcores) )
    {
    }

  if (Kernel_PhysicalProcessorID() == ptr->mastercore)
    {
      if(lockup)
	LockBox_FetchAndDec(ptr->ctrl_lock);
      else
	LockBox_FetchAndInc(ptr->ctrl_lock);

      LockBox_FetchAndClear(ptr->status[lockup]);
      LockBox_FetchAndClear(ptr->lock[lockup]);
    }
  else
    {
      // wait until master releases the barrier by clearing the lock
      while( LockBox_Query(ptr->lock[lockup]) > 0 )
	{
	}
    }

  return( (int)status_rc );
}

__END_DECLS



#endif // Add nothing below this line
