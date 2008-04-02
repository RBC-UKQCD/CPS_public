/******************************************************************************/
/*                                                                            */
/*                              VN MODE BARRIER                               */
/*                                                         30.06.07 S.Krieg   */
/******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <spi/bgp_SPI.h>
#include <spi/lockbox_interface.h>
#include <errno.h>
#include "NodeBarrier.h"
//#include "Helpers.h"
//#include "Errors.h"
const int ERR_ALLOC_LOCKBOX_COUNTER = -4;
const int ERR_LOCKBOX_BARRIER =-5;

inline int KillJob(int rc){
  printf("Killed with return code %d\n",rc);
  exit(rc);
}

#define GET_PID Kernel_PhysicalProcessorID()

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
/* __INLINE__ int LockBox_AllocateBarrier(uint32_t desired, LockBox_Barrier_t* ptr, uint32_t mastercore, uint32_t numcores, uint32_t flags); */

static LockBox_Barrier_t mem_LockBox;
static LockBox_Barrier_t *ptr_LockBox = &mem_LockBox;
static uint32_t desired = 197;
const uint32_t masterCore = 0;
const uint32_t numCores = 4;
const uint32_t flags = LOCKBOX_ORDERED_ALLOC;

static int VN_AllocateBarrier(void){
    int returnCode = 1;
    returnCode = LockBox_AllocateBarrier( desired,     //uint32_t desired, 
					  ptr_LockBox, //LockBox_Barrier_t* ptr, 
					  masterCore,  //uint32_t mastercore, 
					  numCores,    //uint32_t numcores, 
					  flags);      //uint32_t flags)
    return returnCode; 
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
/* __INLINE__ int LockBox_Barrier(LockBox_Barrier_t* ptr) */
static int LockBoxAllocated = 0;

    
void VN_NodeBarrier(void){
    char *fname = "VN_Barrier";
    int Pid = GET_PID;
    extern int errno;
    
    if (!LockBoxAllocated){
	if (VN_AllocateBarrier()){
	    printf("%s: pid %i: ERROR fixed desired LockBox counter %i allocate failed. errno = %i\n",
		   fname, Pid, desired, errno);
	    fflush(stdout);
	    KillJob (ERR_ALLOC_LOCKBOX_COUNTER);
	}
	LockBoxAllocated = 1;
    }
    
    /* _bgp_msync(); */ /* this is included in the current lockbox spi */
    if (LockBox_Barrier(ptr_LockBox)){
	printf("%s: pid %i: ERROR during LockBox_Barrier, counter %i, errno = %i\n",
	       fname, Pid, desired, errno);
	fflush(stdout);
	KillJob (ERR_LOCKBOX_BARRIER);
    }
    return;
}
