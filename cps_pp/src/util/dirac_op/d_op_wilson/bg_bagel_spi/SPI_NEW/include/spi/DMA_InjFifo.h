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

#ifndef	_DMA_INJFIFO_H_ // Prevent multiple inclusion
#define	_DMA_INJFIFO_H_


/*!
 * \file spi/DMA_InjFifo.h
 *
 * \brief DMA SPI Injection Fifo Definitions and Inline Functions
 *
 * This include file contains inline functions that are used to interface with
 * BG/P DMA injection fifos at the lowest level.
 * Functions include
 * - initialize
 * - get fifo start, head, tail, end, size, free space, descriptor count
 * - set fifo head, tail, start PA, head PA, tail PA, end PA
 * - increment tail
 * - inject descriptor(s)
 * - query status: not empty, available, threshold crossed, activated,
 *   descriptor done.
 * - set status: clear threshold crossed, activate, deactivate
 *
 * Data structures are defined to manipulate the injection fifos:
 * - An injection fifo group structure defining a group of injection fifos
 * - Within the group are injection fifo structures
 * - Within each injection fifo structure is a software fifo structure
 * - Each software fifo structure points to its corresponding hardware
 *   fifo structure in the DMA SRAM
 * 
 * \verbatim Picture of data structures:
  
   ========DDR MEMORY===================|==========DMA SRAM MEMORY==========
   ------------------------------       |
   | DMA_InjFifoGroup_t         |       |
   |                            |       |     -----------------------------
   | status --------------------|-------|---->| DMA_InjFifoStatus_t       |
   | fifo[0..31]                |       |     -----------------------------
   |   ------------------------ |       | 
   |   | DMA_InjFifo_t        | |       |
   |   |                      | |       |
   | 0 |  ------------------- | |       |     -----------------------------
   |   |  | DMA_Fifo_t      |-|-|-------|---->| DMA_FifoHW_t              |
   |   |  ------------------- | |       |     ----------------------------- 
   |   ------------------------ |       |
   |             .              |       |
   |             .              |       |
   |             .              |       |
   |   ------------------------ |       | 
   |   | DMA_InjFifo_t        | |       |
   |   |                      | |       |
   |31 |  ------------------- | |       |     -----------------------------
   |   |  | DMA_Fifo_t      |-|-|-------|---->| DMA_FifoHW_t              |
   |   |  ------------------- | |       |     ----------------------------- 
   |   ------------------------ |       |
   ------------------------------       |
  
   \endverbatim
 *
 * Definitions:
 * - A fifo represents a contiguous block of DDR memory
 * - A fifo has a starting address and an ending address (defines the memory 
 *   block)
 * - An injection fifo is a series of 32-byte descriptors.  There is a count
 *   of the number of descriptors ever injected into this fifo.  It will never
 *   wrap in the expected lifetime of a job.
 * - Injection consists of copying a 32-byte descriptor into the next available 
 *   slot (pointed to by the tail), incrementing the tail pointer, and
 *   incrementing the descriptor count for the fifo.
 * - The DMA engine asynchronously processes descriptors, beginning with the 
 *   descriptor pointed to by head, and ending with the descriptor just prior
 *   to tail.
 * - There are injection (DMA InjFifo) and reception (DMA RecFifo) fifos 
 *   (separate interfaces)
 * - There are DMA_NUM_INJ_FIFO_GROUPS injection fifo groups
 * - There are DMA_NUM_INJ_FIFOS_PER_GROUP injection fifos per group
 * - Thus, there are DMA_NUM_INJ_FIFOS injection fifos per node
 * - There are DMA_NUM_REC_FIFO_GROUPS reception fifo groups
 * - There are DMA_NUM_REC_FIFOS_PER_GROUP reception fifos per group
 * - Thus, there are DMA_NUM_REC_FIFOS reception fifos per node
 * - A "shadow" refers to a copy of the elements of the fifo (start, end, head, 
 *   tail) that is maintained by these inline functions.  The shadows may be 
 *   used to calculate other values such as free space.  The shadows are updated
 *   by these inlines whenever the hardware fifo is read or written.
 *
 * \note These functions do not try to detect things that software shouldn't do, 
 *       like injecting a descriptor into a remote_get fifo, since the hardware 
 *       doesn't distinguish between remote get fifos and normal injection 
 *       fifos.  That sort of checking should be done in a higher level. 
 *
 * \note Memory consistency/coherency inside these inlines is achieved using
 *       mbar and msync.  
 * 
 *       MBAR is used to make sure that all writes to memory issued by the
 *       calling core have been accepted by the memory system before 
 *       continuing.  This guarantees that writes and reads to/from different 
 *       addresses to go in defined order.
 *       
 *       MBAR EXAMPLE 1:  When a store is done to DMA SRAM, it may not complete 
 *       for a period of time.  If a counter value is set, and then an injection
 *       fifo tail pointer is set, DMA may see the tail pointer update and begin 
 *       the operation before the counter value has been set.  Inserting an mbar 
 *       between the setting of the counter and the setting of the tail pointer 
 *       guarantees that the counter will be set before the tail pointer is 
 *       updated.
 *
 *       MBAR EXAMPLE 2:  A counter hits zero.  We process the hit-zero and write
 *       a "clear hit zero" to DMA SRAM, and then go read that counter's hit-zero
 *       status (different address).  The hit-zero status will still indicate 
 *       that it hit zero, even though we have already processed it, unless an 
 *       mbar is inserted between clearing the hit-zero and reading the hit-zero
 *       status.
 *   
 *       MBAR PHILOSOPHY:  After DMA SRAM is updated in the DMA inline functions,
 *       they always do at least an mbar (possibly an msync instead...see below).
 *
 *       MSYNC does what mbar does, plus ensures consistency across cores.  That 
 *       is, it waits for snoops (invalidations of L1 cache) on the other cores
 *       to complete before continuing.  This guarantees that all of the cores
 *       will see a consistent view of memory after the msync.
 *
 *       MSYNC EXAMPLE:  When a reception counter has hit zero, we assume the 
 *       DMA'd data is available to be read by any core.  However, old copies of 
 *       that data may still be in the L1 caches.  Inserting an msync after 
 *       detecting that a counter has hit zero guarantees that the old data has 
 *       been removed from the L1 caches.
 *
 *       MSYNC PHILOSOPHY:  After the inline functions detect that a counter has 
 *       hit zero, they always do an msync.
 *
 *       SPECULATIVE EXECUTION OF MSYNC:  There are cases where msync is done
 *       conditionally.  The CPU will begin execution of both sides of the
 *       condition before the result of the condition has been determined.
 *       Then, it will cancel the execution of one side once the result of the
 *       condition has been determined.  This speculation is unwanted when
 *       the first instruction on one side of the condition is msync because
 *       cancelling an msync is similar to executing the complete msync.
 *       To avoid this speculative execution of msync, we call 
 *       _bgp_msync_nonspeculative().  This will trick the CPU so it won't begin
 *       the msync until the result of the condition is known.
 *
 *       CALLER ADVICE:  Users of these functions should not need to do 
 *       mbar/msync themselves, unless they are doing something like the 
 *       following:  Read a counter and operate on the result when the counter
 *       hasn't reached zero.  The caller will need to perform an msync after
 *       reading the counter in order to ensure that snoops have completed
 *       on all CPUs before operating on the DMA'd data.
 *
 * \note General discussion on injection fifo interrupts.  Both the warning 
 * threshold crossed and full fifo interrupts...
 *
 * For remote gets, a fifo is considered available if it has at least 512 bytes
 * free (32 16B quads). An arriving remote get can be written if there are 512 
 * bytes free, but after that the available goes low and no further remote gets 
 * can be written to any fifo.  Furthermore, if any injection fifo has less than
 * 512 bytes free, the fifo becomes unavailable and any arriving remote get 
 * packet will cause an interrupt to fire and the rDMA will stop. 
 *
 * Specifically, if an injection fifo has less than 512 B (by either injecting 
 * or remote gets) the iDMA will continue to operate and the rDMA will continue
 * to operate until any remote get packet arrives to any fifo, at which point
 * an interrupt fires and the rDMA stops.  
 *
 * Note that these interrupts were put in for warnings of remote get fifos 
 * becoming nearly full. However the time between when the warning fires and the
 * condition is cleared may be long, reconfiguring an almost full remote get 
 * fifo is difficult, and recovery from full remote get injection fifos is very 
 * difficult.  Since software can prevent this, and since recovery is so 
 * difficult, we consider injection fifo threshold crossing interrupts and 
 * injection fifo full interrupts to be fatal. Thus there is no handler function
 * in the injection fifo allocation routine.  
 *
 * So software needs to manage injection and remote get fifo space so that there
 * are always at least 512 bytes of free space in every fifo. To accomplish 
 * this, software needs to guarantee it won't inject descriptors if doing so 
 * would trigger an interrupt or make the fifo unavailable.  
 *
 * This can be done by setting the interrupt threshold to 0 (interrupt fires if 
 * free space <= threshhold), and not injecting if after injection there are 
 * less than DMA_MIN_INJECT_SIZE_IN_QUADS (=32) slots.  Furthermore, remote 
 * get space should not be allocated if doing so might result in strictly less 
 * than DMA_MIN_INJECT_SIZE_IN_QUADS slots.
 *
 */

#include <common/namespace.h>
#include <memory.h>


__BEGIN_DECLS


/*!
 * \brief __INLINE__ definition
 * 
 * Option 1:
 * Make all functions be "static inline":
 * - They are inlined if the compiler can do it
 * - If the compiler does not inline it, a single copy of the function is
 *   placed in the translation unit (eg. xxx.c)for use within that unit.
 *   The function is not externalized for use by another unit...we want this
 *   so we don't end up with multiple units exporting the same function,
 *   which would result in linker errors.
 *
 * Option 2:
 * A GNU C model: Use "extern inline" in a common header (this one) and provide
 * a definition in a .c file somewhere, perhaps using macros to ensure that the
 * same code is used in each case. For instance, in the header file:
 *
   \verbatim
   #ifndef INLINE
   # define INLINE extern inline
   #endif
   INLINE int max(int a, int b) {
     return a > b ? a : b;
   }
   \endverbatim
 *
 * ...and in exactly one source file (in runtime/SPI), that is included in a
 * library...
 *
   \verbatim
   #define INLINE
   #include "header.h"
   \endverbatim
 * 
 * This allows inlining, where possible, but when not possible, only one 
 * instance of the function is in storage (in the library).
 */
#ifndef __INLINE__
#define __INLINE__ extern inline
#endif



#include <spi/DMA_Assert.h>
#include <spi/DMA_Fifo.h>
#include <spi/DMA_Descriptors.h>




/*!
 * \brief Number of Injection Fifo Groups
 */
#define DMA_NUM_INJ_FIFO_GROUPS     4


/*!
 * \brief Number of Injection Fifos per Group
 */
#define DMA_NUM_INJ_FIFOS_PER_GROUP 32


/*!
 * \brief Number of Injection Fifos (total)
 */
#define DMA_NUM_INJ_FIFOS (DMA_NUM_INJ_FIFO_GROUPS*DMA_NUM_INJ_FIFOS_PER_GROUP)


/*!
 * \brief Minimum Free Space Required After Injection
 *
 * This is the number of 16-byte quads that need to be free in a fifo after
 * injection of a descriptor.
 */
#define DMA_MIN_INJECT_SIZE_IN_QUADS 32


/*!
 * \brief Number of 16-byte quads in a fifo descriptor
 *
 */
#define DMA_FIFO_DESCRIPTOR_SIZE_IN_QUADS 2


/*!
 * \brief Number of bytes in a fifo descriptor
 *
 */
#define DMA_FIFO_DESCRIPTOR_SIZE_IN_BYTES DMA_FIFO_DESCRIPTOR_SIZE_IN_QUADS*16 


/*!
 * \brief Minimum size of a fifo, somewhat arbitrary
 */
#define DMA_MIN_INJ_FIFO_SIZE_IN_BYTES (256*4)


/*!
 * \brief Injection DMA Fifo Structure
 *
 * This structure contains a software DMA fifo structure (defined in DMA_Fifo.h)
 * and other fields that are specific to an injection fifo used by software.
 *
 * \todo   Some more careful thought should be given how to group these so as to
 *         get best memory system performance.
 *         eg.  Probably want to ALIGN_L3_CACHE the fifo_hw_ptr.
 *
 */
typedef struct DMA_InjFifo_t
{
  DMA_Fifo_t         dma_fifo;   /*!< Common software fifo structure          */
  unsigned short int fifo_id;    /*!< The fifo identifier (0 to                 
                                      DMA_NUM_INJ_FIFOS_PER_GROUP-1).         */

  unsigned long long desc_count; /*!< The number of descriptors that have       
                                      ever been injected into this fifo.      */

  unsigned int occupiedSize;     /*!< The number of 16B quads in the fifo that  
                                      are logically occupied.  This does not    
                                      include the DMA_MIN_INJECT_SIZE_IN_QUADS  
                                      that always remains logically occupied. */
  /*!
   * \note The following fields contain info about the fifo that affects the 
   *       DCR values configuring the fifo.
   */
  unsigned short int priority;   /*!< 0 = Normal priority, 1 = High priority.   
                                      The DMA uses this to determine which      
                                      injection fifo to serve next.             
                                      Reflected in DCR addresses                
                                      _BGP_DCR_iDMA_FIFO_PRIORITY(i), where i   
                                      is the group_id.  0xD32 - 0xD35.          
                                      Fifo j is high priority if bit j in the   
                                      DCR is 1, otherwise it is normal          
                                      priority.                               */

  unsigned short int local;      /*!< 0 = non-local, 1 = local.                 
                                      If 0, this fifo uses the torus and        
                                      ts_inj_map must be non-zero.              
                                      If 1, this fifo is used for tranfsers     
                                      local to the node only.                   
                                      Reflected in DCR addresses                
                                      _BGP_DCR_iDMA_LOCAL_COPY(i), where i      
                                      is the group_id.  0xD5C - 0xD5F.             
                                      Fifo j is for local transfers if bit j    
                                      in the DCR is 1, otherwise it is for      
                                      torus transfers.                        */

  unsigned char ts_inj_map;      /*!< 8 bit vector mask indicating which torus  
                                      fifos can be used by this DMA fifo.       
                                      Reflected in DCR addresses                
                                      _BGP_DCR_iDMA_TS_INJ_FIFO_MAP(k) where k  
                                      is the fifo_id.  0xD3C - 0xD5B.           
                                      Fifo k can inject in torus fifo j if      
                                      bit j of the k'th DCR byte is 1.        */
}
DMA_InjFifo_t;


/*!
 * \brief DMA Injection Fifo Status structure
 *
 * This structure maps the DMA SRAM for a particular group of 
 * DMA_NUM_INJ_FIFOS_PER_GROUP fifos.
 *
 */
typedef struct DMA_InjFifoStatus_t
{
  volatile unsigned  not_empty;              /*!< R bitmask, 1 bit/fifo:        
                                                    Injection FIFO not empty. */

  volatile unsigned  reserved_0;             /*!< HOLE                        */

  volatile unsigned  available;              /*!< R bitmask, 1 bit/fifo:        
                                                    Injection FIFO available. */

  volatile unsigned  reserved_1;             /*!< HOLE                        */

  volatile unsigned  threshold_crossed;      /*!< R bitmask, 1 bit/fifo:        
                                                    Threshold crossed.        */

  volatile unsigned  reserved_2;             /*!< HOLE                        */

  volatile unsigned  clear_threshold_crossed;/*!< W bitmask, 1 bit/fifo:        
                                                    Clear threshold crossed.  */

  volatile unsigned  reserved_3;             /*!< HOLE                        */

  volatile unsigned  activated;              /*!< R bitmask, 1 bit/fifo:        
                                                    Retrieve activated fifos. */

  volatile unsigned  activate;               /*!< W bitmask, 1 bit/fifo:        
                                                    Set "1" to activate fifo. */

  volatile unsigned  deactivate;             /*!< W bitmask, 1 bit/fifo:        
                                                    Set "1" to deactivate fifo*/
}
DMA_InjFifoStatus_t;


/*!
 * \brief DMA Injection Fifo Group Structure
 *
 * This structure defines a DMA InjFifo Group.  It points to a
 * DMA InjFifo Status structure, and contains DMA_NUM_INJ_FIFOS_PER_GROUP 
 * DMA InjFifo structures.
 *
 * It is passed into the DMA_InjFifoGroupAllocate system call.  
 * The system call sets up the requested fifos, and fills in this fifo group
 * structure, including the appropriate DMA InjFifo structures within it.
 *
 * It also contains permission bits to use the fifos, one bit per fifo.
 * When the permission bit is on, the corresponding fifo belongs to this
 * group and can be used.  Otherwise, the fifo should not be used as part
 * of this group.  These permission bits are used as follows:
 *   1. Inline functions will ASSERT when an attempt is made
 *      to use a fifo that is not part of this group.
 *   2. Inline functions will use the permission bits as a mask 
 *      to return status information only for fifos that are allocated
 *      to this group.
 *
 */
typedef struct DMA_InjFifoGroup_t
{
  DMA_InjFifoStatus_t *status_ptr;  /*!< Pointer to fifo status.              */

  DMA_InjFifo_t        fifos[DMA_NUM_INJ_FIFOS_PER_GROUP];/*!< Array            
                                         of fifo structures.  The i-th struct   
                                         is defined and usable only if          
                                         bit i of permissions = 1.            */   

  unsigned int permissions;         /*!< Permissions bit vector.  Bit i is 1    
                                         if permitted to use fifo i.  The fifo            
                                         is allocated to this group.          */

  unsigned int group_id;            /*!< The id of this group (0 to            
                                         DMA_NUM_INJ_FIFO_GROUPS-1).          */
}
DMA_InjFifoGroup_t;


/*!
 * \brief Remote Get Fifo Full Handler Function Prototype
 *
 * A function with this signature receives control when one or more remote 
 * get fifos have filled.  This function should do the following to help 
 * make space in the fifo(s):
 * 1. Determine if there are any remote get fifos full or nearly full.
 * 2. For each such fifo:
 *    1. Allocate a larger fifo
 *    2. Copy the descriptors from the old fifo to the new fifo
 *    3. Call DMA_InjFifoInitById() to register the new fifo with the DMA
 *    4. Call DMA_InjFifoSetTailById() to set the new fifo's tail pointer
 *    5. Free the old fifo
 *
 * A function of this type can be registered on DMA_InjFifoGroupAllocate().
 *
 * \param[in]  fg_ptr  Pointer to the fifo group associated with this fifo.
 * \param[in]  f_num   The fifo number that has filled.  This is
 *                     relative to the DMA fifo group
 *                     (0 to DMA_NUM_INJ_FIFOS_PER_GROUP-1).
 * \param[in]  handler_param  An opaque pointer provided by the caller who 
 *                            registered this handler.
 */
typedef void (*DMA_InjFifoRgetFifoFullHandler_t)(
					  DMA_InjFifoGroup_t *fg_ptr,
                                          int                 f_num,
                                          void               *handler_parm
					 );


/*!
 * \brief Remote Get Fifo Full Handler Table Entry
 *
 * This defines an entry in the Remote Get Fifo Full Handler Table.
 * It identifies the fifo group pointer associated with the full fifo,
 * and the pointer to the handler function to receive control to handle
 * the fifo full condition and the opaque pointer to be passed to the
 * handler function when it is called.  The core number of the core that
 * will process the condition is associated with each entry.  
 */ 
typedef struct DMA_InjFifoRgetFifoFullHandlerEntry_t
{
  DMA_InjFifoGroup_t               *fg_ptr;  /*!< Pointer to injection fifo group    */
  DMA_InjFifoRgetFifoFullHandler_t  handler; /*!< Pointer to handler function */
  void                             *handler_parm; /*!< Pointer to be passed to
                                                       the handler.           */
  uint32_t                          core_num;/*!< Core number of the core that
                                                  will process the condition. */
} DMA_InjFifoRgetFifoFullHandlerEntry_t;


/*!
 * 
 * \brief Remote Get Fifo Full Handler Table
 *
 * An array of entries, one per injection fifo.  Each entry specifies the fifo
 * group structure and the handler function that will receive control to
 * handle a remote get fifo full condition for fifos in that fifo group.
 */
extern DMA_InjFifoRgetFifoFullHandlerEntry_t DMA_RgetFifoFullHandlerTable[DMA_NUM_INJ_FIFOS];


/*!
 * \brief Remote Get Fifo Full Init Has Been Done Indicator
 *
 *  0 means the initialization has not been done.
 *  1 means the initialization has been done.
 */
extern int DMA_InjFifoRgetFifoFullInitHasBeenDone;


/*!
 * \brief Remote Get Fifo Full Initialization
 *
 * Initialize data structures and interrupt handlers to handle a remote get
 * fifo full condition.
 *
 * \param[in]  interruptGroup  The handle that identifies the remote get fifo
 *                             full interrupts (only one interrupt, in this
 *                             case, group 3, irq 24).
 * \param[in]  rget_barrier    A barrier that is used by the handler function
 *                             to synchronize all cores in the node as they
 *                             each handle the interrupt (it is a broadcasted
 *                             interrupt).
 */
void DMA_InjFifoRgetFifoFullInit( Kernel_InterruptGroup_t   interruptGroup,
				  LockBox_Barrier_t        *rget_barrier  );


/*!
 * \brief Query Free DMA InjFifos within a Group
 *
 * This function is a wrapper around a system call that returns a list of the 
 * free (available to be allocated) fifos within the specified group.
 *
 * \param[in]   grp            Group number being queried 
 *                             (0 to DMA_NUM_INJ_FIFOS_PER_GROUP-1)
 * \param[out]  num_fifos      Pointer to an int where the number of free
 *                             fifos in the specified group is returned
 * \param[out]  fifo_ids       Pointer to an array of num_fifos short ints where
 *                             the list of free fifos is returned.
 *                             Each short int is the fifo number 
 *                             (0 to DMA_NUM_INJ_FIFOS_PER_GROUP-1).
 *                             The caller must provide space for 
 *                             DMA_NUM_INJ_FIFOS_PER_GROUP ints,
 *                             in case the entire fifo group is free.
 *
 * \retval  0  Successful.  num_fifos and fifo_ids array set as described.
 * \retval  -1 Unsuccessful.  errno gives the reason.
 *
 */
__INLINE__ int  DMA_InjFifoGroupQueryFree(
					  int  grp,
					  int *num_fifos,
					  int *fifo_ids
					 )
{
  return Kernel_InjFifoGroupQueryFree( grp, 
				       (uint32_t*)num_fifos, 
				       (uint32_t*)fifo_ids);
}


/*!
 * \brief Allocate DMA InjFifos From A Group 
 *
 * This function is a wrapper around a system call that allocates specified 
 * DMA injection fifos from the specified group.  Parameters specify whether 
 * each fifo is high or normal priority, local or non-local, and which torus 
 * fifos it maps to.  A DMA_InjFifoGroup_t structure is returned for 
 * use in other inline functions to operate on the allocated fifos.
 * 
 * Refer to the interrupt discussion at the top of this include file to see why
 * there are no interrupt-related parameters.
 *
 * \param[in]   grp          Group number whose DMA injection fifos are being 
 *                           allocated (0 to DMA_NUM_INJ_FIFO_GROUPS-1)
 * \param[in]   num_fifos    Number of fifos to be allocated from the group 
 *                           (1 to DMA_NUM_INJ_FIFOS_PER_GROUP)
 * \param[in]   fifo_ids     Pointer to an array of num_fifos ints where
 *                           the list of fifos to be allocated is provided.
 *                           Each int is the fifo number 
 *                           (0 to DMA_NUM_INJ_FIFOS_PER_GROUP-1).
 * \param[in]   priorities   Pointer to an array of num_fifos short ints where
 *                           the list of priorities to be assigned to the fifos
 *                           is provided.  Each short int indicates the priority
 *                           to be assigned to each of the fifos identified in
 *                           the fifo_ids array (0 is normal, 1 is high priority).
 * \param[in]   locals       Pointer to an array of num_fifos short ints where
 *                           an indication is provided of whether each fifo will 
 *                           be used for local transfers (within the same node) 
 *                           or torus transfers.  Each short int indicates the 
 *                           local/non-local attribute to be assigned to each of 
 *                           the fifos identified in the fifo_ids array (0 is
 *                           non-local, 1 is local).  If 0, the corresponding 
 *                           array element in ts_inj_maps indicates which torus 
 *                           fifos can be injected.
 * \param[in]   ts_inj_maps  Pointer to an array of num_fifos chars where
 *                           the torus fifos that can be injected are specified 
 *                           for each fifo.  Each char specifies which of 
 *                           the 8 torus injection fifos can be injected when a 
 *                           descriptor is injected into the DMA injection fifo.
 *                           Must be non-zero when the corresponding "locals" 
 *                           is 0.
 *                           Bits 0-3 are for torus group 0.
 *                           Bits 4-7 are for torus group 1.
 *                           Bits 3 and 7 are the high priority fifos.
 * \param[in]   rget_handler Pointer to a function with prototype
 *                           DMA_InjFifoRgetFifoFullHandler_t that will handle
 *                           a remote get fifo full condition for fifos in this
 *                           fifo group.  If NULL is specified, the condition 
 *                           will not be handled.
 * \param[in]   rget_handler_parm   A pointer to opaque storage that will be 
 *                                  passed to the rget_handler.
 * \param[in]   rget_interruptGroup  A InterruptGroup_t that identifies the
 *                           group of interrupts that handle the remote get
 *                           fifo full condition.  It is only one interrupt:
 *                           group 3, irq 24.
 * \param[in]   rget_barrier A barrier that is used by the rget fifo full 
 *                           interrupt handler.  This barrier should be across
 *                           all processes on this compute node.
 * \param[out]  fg_ptr       Pointer to a structure that is filled in upon 
 *                           successful return for use in other inline functions 
 *                           to operate on the allocated fifos.
 *                           \li fifos - Array of fifo structures.  Structures
 *                                       for allocated fifos are initialized as 
 *                                       documented below.  Structures for
 *                                       fifos not allocated by this instance of
 *                                       this syscall are initialized to binary
 *                                       zeros.  Allocated fifos are enabled.
 *                           \li status_ptr  - Points to status area within the 
 *                                             DMA memory map.
 *                           \li permissions - Bits indicating which fifos were 
 *                                             allocated during this syscall.
 *                           \li group_id    - The id of this group.
 *
 * \retval  0  Successful.  Fifos allocated and fg_ptr structure filled in as 
 *                          described.
 * \retval  -1 Unsuccessful.  errno gives the reason.
 *
 * \return The group fifo structure pointed to by fg_ptr is completely 
 *         initialized as follows:
 *         - status_ptr points to the appropriate fifo group DMA memory map
 *         - fifo structures array.  Fifo structures for fifos not allocated
 *           during this syscall are initialized to binary zeros.  Fifo 
 *           structures for fifos allocated during this syscall are initialized:
 *             - fifo_hw_ptr points to the DMA memory map for this fifo.  The 
 *               hardware start, end, head, and tail are set to zero by the 
 *               kernel.
 *             - All other fields in the structure are set to zero by the kernel
 *               except priority, local, and ts_inj_map are set to reflect what 
 *               was requested in the priorities, locals, and ts_inj_maps 
 *               syscall parameters.
 *
 */
__INLINE__ int  DMA_InjFifoGroupAllocate(
			 int                               grp,
			 int                               num_fifos,
		         int                              *fifo_ids,
			 unsigned short int               *priorities,
			 unsigned short int               *locals,
			 unsigned char                    *ts_inj_maps,
                         DMA_InjFifoRgetFifoFullHandler_t  rget_handler,
                         void                             *rget_handler_parm,
                         Kernel_InterruptGroup_t           rget_interruptGroup,
                         LockBox_Barrier_t                *rget_barrier,
			 DMA_InjFifoGroup_t               *fg_ptr
					)
{
  int rc;
  int i, global_fifo_id;

  rc = Kernel_InjFifoGroupAllocate( grp, 
				    num_fifos, 
				    (uint32_t*)fifo_ids, 
				    (uint16_t*)priorities, 
				    (uint16_t*)locals, 
				    (uint8_t*)ts_inj_maps, 
				    (uint32_t*)fg_ptr);

  if ( rc == 0 )
  {
     /*
      * If a remote get fifo full handler has been provided, update the table
      * to indicate that this handler will handle full conditions on the fifos
      * just allocated.
      */
     if ( rget_handler )
     {
        /*
         * If rget handler init has not been done, do it:
         */
        if ( DMA_InjFifoRgetFifoFullInitHasBeenDone == 0 )
          DMA_InjFifoRgetFifoFullInit( rget_interruptGroup,
				       rget_barrier );

        for (i=0; i<num_fifos; i++)
        {
           global_fifo_id = (grp * DMA_NUM_INJ_FIFOS_PER_GROUP) + fifo_ids[i];
           DMA_RgetFifoFullHandlerTable[global_fifo_id].fg_ptr  = fg_ptr;
           DMA_RgetFifoFullHandlerTable[global_fifo_id].handler = rget_handler;
	   DMA_RgetFifoFullHandlerTable[global_fifo_id].handler_parm = 
                                                             rget_handler_parm;
           DMA_RgetFifoFullHandlerTable[global_fifo_id].core_num= 
                                                  Kernel_PhysicalProcessorID();
        }

	/*
	 * Indicate done with initialization.
	 */
	DMA_InjFifoRgetFifoFullInitHasBeenDone = 1;
     }
  }

  return(rc);
}


/*!
 * \brief Free DMA InjFifos From A Group 
 *
 * This function is a wrapper around a system call that frees DMA injection
 * counters from the specified group.
 *
 * \param[in]   grp          Group number whose DMA injection fifos are being 
 *                           freed (0 to DMA_NUM_INJ_FIFO_GROUPS-1)
 * \param[in]   num_fifos    Number of fifos to be freed from the group 
 *                           (1 to DMA_NUM_INJ_FIFOS_PER_GROUP)
 * \param[in]   fifo_ids     Pointer to an array of num_fifos ints where
 *                           the list of fifos to be freed is provided.
 *                           Each int is the fifo number (0 to num_fifos-1).
 * \param[in]   fg_ptr       Pointer to the structure previously filled in when 
 *                           these fifos were allocated.  Upon successful 
 *                           return, this structure is updated to reflect the 
 *                           freed fifos:
 *                           \li fifos - Structures for freed fifos zero'd.
 *                                       Freed fifos are disabled.
 *                           \li permissions - Bits cleared for each freed fifo.
 *
 * \retval  0  Successful.  Fifos freed and fg_ptr structure updated as described.
 * \retval  -1 Unsuccessful.  errno gives the reason.
 *
 * \note  This is a fatal error if any of the fifos are non empty and activated
 *
 */
__INLINE__ int  DMA_InjFifoGroupFree(
				     int                 grp,
				     int                 num_fifos,
				     int                *fifo_ids,
				     DMA_InjFifoGroup_t *fg_ptr
				    )
{
  return Kernel_InjFifoGroupFree( grp, 
				  num_fifos, 
				  (uint32_t*)fifo_ids, 
				  (uint32_t*)fg_ptr);
}




/*
 * -----------------------------------------------------------------------------
 * Calls to access the Fifo, given a pointer to the injection fifo structure
 * -----------------------------------------------------------------------------
 */




/*!
 * \brief Set DMA Injection Fifo Head
 *
 * Set a DMA injection fifo's "head", given an injection fifo structure
 *
 * \param[in]  f_ptr    Pointer to the injection fifo structure
 * \param[in]  va_head  Virtual address of the head to be set
 *
 * \return  None
 *
 * \post va_head is set in both the hardware and software fifo structures,
 *       and the fifo free space is recalculated.
 *
 * \note Normally, for an injection fifo, the dma manipulates the head, but in
 *       optimized persistant communications the core can do it if it is sure 
 *       the fifo is empty at the time this is called.
 */
__INLINE__ void DMA_InjFifoSetHead(  
				   DMA_InjFifo_t *f_ptr,
				   void        *va_head 
				  )
{
  SPI_assert( f_ptr != NULL );

  DMA_FifoSetHead( &f_ptr->dma_fifo,
		   va_head );
}


/*!
 * \brief Increment DMA Injection Fifo Tail
 *
 * Increment a DMA injection fifo's "tail", given an injection fifo structure
 *
 * \param[in]  f_ptr  Pointer to the injection fifo structure
 * \param[in]  incr   The number of quads (16 byte units) to increment the 
 *                    tail pointer by.  This value must be even (ie. descriptors
 *                    are 32 bytes).
 *
 * \retval  None
 *
 * \post va_tail is set in both the hardware and software fifo structures,
 *       the fifo free space is recalculated, and the fifo's descriptor count
 *       is incremented according to the incr.
 *
 * \note This function does not check if there is free space in the fifo 
 *       for this many quads.  It must be preceeded by a check of the 
 *       free space.
 */
__INLINE__ void DMA_InjFifoIncrementTail(  
					 DMA_InjFifo_t *f_ptr,
					 unsigned int   incr 
					)
{
  SPI_assert( f_ptr != NULL );
  SPI_assert( (incr & 0x1) == 0 );

  void *va_tail = DMA_FifoGetTailFromShadow( &f_ptr->dma_fifo );

  void *va_end  = DMA_FifoGetEndFromShadow( &f_ptr->dma_fifo );

  unsigned int incr_bytes = incr << 4;

  unsigned int bytes_to_end = (unsigned)va_end - (unsigned)va_tail;

  /* 
   * Note:  The following check must be >= instead of just >.  We never want
   *        the tail to be equal to the end so we can always copy a descriptor
   *        to the tail, safely.
   */
  if ( incr_bytes >= bytes_to_end )
    {
      va_tail = (char *)
	          ( (unsigned)DMA_FifoGetStartFromShadow( &f_ptr->dma_fifo ) + 
		    ( incr_bytes - bytes_to_end ) );
    }
  else
    {
      va_tail = (char *)( (unsigned)va_tail + incr_bytes );
    }

  DMA_FifoSetTail( &f_ptr->dma_fifo,
		   va_tail );
  
  f_ptr->desc_count += (incr >> 1);

}


/*!
 * \brief Get DMA Injection Fifo Descriptor Count
 *
 * Get a DMA injection fifo's "descriptor count", given an injection fifo 
 * structure
 *
 * \param[in]  f_ptr  Pointer to the injection fifo structure
 *
 * \retval  desc_count  The descriptor count for the specified fifo
 *
 */
__INLINE__ unsigned long long DMA_InjFifoGetDescriptorCount( 
						      DMA_InjFifo_t *f_ptr
						     )
{
  SPI_assert( f_ptr != NULL );

  return f_ptr->desc_count;
} 


/*!
 * \brief Is DMA Descriptor Done
 *
 * Return whether a specified descriptor is still in the specified injection
 * fifo (not done).  The descriptor is identified by the descriptor count
 * immediately after the descriptor was injected into the fifo (returned by
 * DMA_InjFifoIncrementTail().
 *
 * \param[in]  f_ptr       Pointer to the injection fifo structure
 * \param[in]  desc_count  The descriptor count immediately after the
 *                         descriptor in question was injected into 
 *                         the fifo.
 * \param[in]  update      0  Do not update the fifo's shadow information.
 *                         1  Update the fifo's shadow information.
 *                         It is a performance optimization to only update the
 *                         shadow information once for a group of descriptors
 *                         being processed.
 *
 * \retval  0  False.  The descriptor identified by desc_count is not done.
 *                     It is still in the fifo.
 * \retval  1  True.   The descriptor identified by desc_count is done.
 *                     It is no longer in the fifo.
 *
 */
__INLINE__ unsigned int DMA_InjFifoIsDescriptorDone( 
					    DMA_InjFifo_t      *f_ptr,
					    unsigned long long  desc_count,
					    unsigned int        update
					   )
{
  unsigned int num_desc_in_fifo;
  unsigned int free_space;
  DMA_Fifo_t  *fifo_ptr;

  SPI_assert( f_ptr != NULL );

  fifo_ptr = &(f_ptr->dma_fifo);

  /* If caller wants a fresh look in the fifo, update its free space. 
   * Otherwise, fetch the free space based on shadows.
   */
  if (update)
    free_space = DMA_FifoGetFreeSpace (fifo_ptr, 1, 0);
  else
    free_space = DMA_FifoGetFreeSpaceNoUpdateCalculation(fifo_ptr);

  /* Compute the desc_count of the oldest descriptor in the fifo (minus 1)
   * Note:  Each desc is a 32B unit and the below are 16B entities
   */
  num_desc_in_fifo = ( DMA_FifoGetSize(fifo_ptr) - free_space ) / 2;

  /* Determine if the specified desc_count is still in the fifo.
   * We take the current descriptor count for this fifo and subtract the
   * number of descriptors still in the fifo.  This is the descriptor count
   * of the oldest descriptor still remaining in the fifo (minus 1).  
   * We compare that with the caller's desc_count to determine if the
   * caller's descriptor is still in the fifo.
   */
  if (desc_count <= (unsigned) (DMA_InjFifoGetDescriptorCount(f_ptr) -
				num_desc_in_fifo) )
    return (1); /* Descriptor is done */
  else
    return (0); /* Descriptor is not done */

} 


/*!
 * \brief DMA Injection Fifo Reserve Descriptor Storage
 *
 * Reserve storage in a DMA injection fifo for a remote get descriptor, given
 * an injection fifo structure.
 *
 * \param[in]  f_ptr   Pointer to the injection fifo structure
 *
 * \retval  0   Successful.  There was enough space in the fifo and the
 *              storage was reserved.
 * \retval -1   Unsuccessful.  There was not enough space in the fifo.
 *
 * \note Internally, this increments the occupiedSize of the fifo by
 * DMA_FIFO_DESCRIPTOR_SIZE_IN_QUADS. 
 *
 */
__INLINE__ int DMA_InjFifoReserveDescriptorStorage( 
						   DMA_InjFifo_t *f_ptr
						  )
{
  SPI_assert( f_ptr != NULL );

  if ( (DMA_FifoGetSize(&f_ptr->dma_fifo) - f_ptr->occupiedSize) >=
       (DMA_MIN_INJECT_SIZE_IN_QUADS + DMA_FIFO_DESCRIPTOR_SIZE_IN_QUADS) ) {
     f_ptr->occupiedSize += DMA_FIFO_DESCRIPTOR_SIZE_IN_QUADS;
    return (0);
  }
  else {
     return (-1);
  }
} 


/*!
 * \brief DMA Injection Fifo Free Descriptor Storage Reservation
 *
 * Free a reservation for storage for a remote get descriptor in a DMA injection
 * fifo, previously reserved using DMA_InjFifoReserveDescriptorStorageById().
 *
 * \param[in]  f_ptr  Pointer to the injection fifo structure
 *
 * \return  None
 *
 * \note Internally, this decrements the occupiedSize of the fifo by
 * DMA_FIFO_DESCRIPTOR_SIZE_IN_QUADS. 
 *
 */
__INLINE__ void DMA_InjFifoFreeDescriptorStorageReservation( 
						 DMA_InjFifo_t *f_ptr
						)
{
  SPI_assert( f_ptr != NULL );
  SPI_assert( f_ptr->occupiedSize >= DMA_FIFO_DESCRIPTOR_SIZE_IN_QUADS );

  f_ptr->occupiedSize -= DMA_FIFO_DESCRIPTOR_SIZE_IN_QUADS;
} 


/*!
 * \brief Check If An Injection Fifo Has Space For Injection
 *
 * Check if an injection fifo has enough space for a single descriptor to be
 * injected.
 *
 * \param[in]  f_ptr  Pointer to the injection fifo structure
 *
 * \retval  hasSpace  An indicator of whether the fifo has space for a 
 *                    descriptor.
 *                    - 0 (false) means the fifo is full.
 *                    - 1 (true)  means the fifo has space.
 *
 */
__INLINE__ unsigned int DMA_InjFifoHasSpace(  
					    DMA_InjFifo_t       *f_ptr
					   )
{
  SPI_assert( f_ptr != NULL );

  unsigned int free_space;

  /* Get the free space in the fifo using the shadow value */
  free_space = DMA_FifoGetFreeSpace( &f_ptr->dma_fifo,
				     0, /* Use shadow head */
				     0);/* use shadow tail */

  /* 
   * If after injecting, (DMA_FIFO_DESCRIPTOR_SIZE_IN_QUADS is the amount we
   * are going to inject), there is still at least the minimum allowable free 
   * space left in the fifo, go ahead and inject.  We want at least 
   * DMA_MIN_INJECT_SIZE_IN_QUADS free space after injection.
   * 
   * Otherwise, read the hardware head pointer and recalculate the free space,
   * and check again.  Note:  We don't need to read the hardware tail
   * pointer because only software updates that, and we recalculate the
   * free space at that time.
   *
   * If there is still not enough room in the fifo, return 0, indicating that
   * the descriptor could not be injected.
   *
   */
  if ( free_space < DMA_MIN_INJECT_SIZE_IN_QUADS + 
                    DMA_FIFO_DESCRIPTOR_SIZE_IN_QUADS )
    {
      free_space = DMA_FifoGetFreeSpace( &f_ptr->dma_fifo,
					 1,  /* Use hardware head */
					 0); /* Use shadow tail   */

      if ( free_space < DMA_MIN_INJECT_SIZE_IN_QUADS + 
	                DMA_FIFO_DESCRIPTOR_SIZE_IN_QUADS ) return 0;
    }

  return 1; // There is space in the fifo.
}


/*!
 * \brief Inject a Descriptor into a DMA Injection Fifo Without Checking for 
 *        Space
 *
 * Inject a descriptor into a DMA injection fifo, given an injection fifo 
 * structure, without checking to see if there is enough space in the fifo.
 * It is assumed that the caller has already checked for enough space using
 * the DMA_InjFifoHasSpace() function.
 *
 * \param[in]  f_ptr  Pointer to the injection fifo structure
 * \param[in]  desc   A pointer to the descriptor to be injected.
 *                    Must be 16-byte aligned.
 *
 * \retval  numDescInjected  The number of descriptors injected.
 *                           - 1 means it was successfully injected.
 *
 * \see DMA_InjFifoHasSpace()
 */
__INLINE__ int DMA_InjFifoInjectDescriptorNoSpaceCheck(  
					   DMA_InjFifo_t       *f_ptr,
					   DMA_InjDescriptor_t *desc
					  )
{
  SPI_assert( f_ptr != NULL );
  SPI_assert( desc  != NULL );

  char *load_ptr, *store_ptr;

  /*
   * Copy the descriptor to the current va_tail of the fifo.
   * Msync to ensure the descriptor has been written to memory and the L1 caches
   * are in sync.
   * Move the tail past the descriptor so the DMA knows the descriptor is there.
   *   - handle wrapping
   *   - update free space
   *
   */

  if ( ( (unsigned)desc & 0xFFFFFFF0 ) == (unsigned)desc ) /* 16B aligned? */
    {
      load_ptr  = (char*)desc;
      store_ptr = (char*)DMA_FifoGetTailFromShadow( &f_ptr->dma_fifo );
      _bgp_QuadLoad ( load_ptr,     0 );
      _bgp_QuadLoad ( load_ptr+16,  1 );
      _bgp_QuadStore( store_ptr,    0 );
      _bgp_QuadStore( store_ptr+16, 1 );
    }
  else
    {
      memcpy( DMA_FifoGetTailFromShadow( &f_ptr->dma_fifo ),
	      desc,
	      DMA_FIFO_DESCRIPTOR_SIZE_IN_BYTES );
    }

  // _bgp_msync();  mbar is good enough
  _bgp_mbar(); 

  DMA_InjFifoIncrementTail( f_ptr,
			    DMA_FIFO_DESCRIPTOR_SIZE_IN_QUADS );

  return 1; /* Success */
}


/*!
 * \brief Inject a Descriptor into a DMA Injection Fifo
 *
 * Inject a descriptor into a DMA injection fifo, given an injection fifo 
 * structure
 *
 * \param[in]  f_ptr  Pointer to the injection fifo structure
 * \param[in]  desc   A pointer to the descriptor to be injected.
 *                    Must be 16-byte aligned.
 *
 * \retval  numDescInjected  The number of descriptors injected.
 *                           - 0 means it was not injected, most likely because
 *                             the fifo is full.
 *                           - 1 means it was successfully injected
 *
 */
__INLINE__ int DMA_InjFifoInjectDescriptor(  
					   DMA_InjFifo_t       *f_ptr,
					   DMA_InjDescriptor_t *desc
					  )
{
  SPI_assert( f_ptr != NULL );
  SPI_assert( desc  != NULL );

  unsigned int free_space;
  char *load_ptr, *store_ptr;

  /* Get the free space in the fifo using the shadow value */
  free_space = DMA_FifoGetFreeSpace( &f_ptr->dma_fifo,
				     0, /* Use shadow head */
				     0);/* use shadow tail */

  /* 
   * If after injecting, (DMA_FIFO_DESCRIPTOR_SIZE_IN_QUADS is the amount we
   * are going to inject), there is still at least the minimum allowable free 
   * space left in the fifo, go ahead and inject.  We want at least 
   * DMA_MIN_INJECT_SIZE_IN_QUADS free space after injection.
   * 
   * Otherwise, read the hardware head pointer and recalculate the free space,
   * and check again.  Note:  We don't need to read the hardware tail
   * pointer because only software updates that, and we recalculate the
   * free space at that time.
   *
   * If there is still not enough room in the fifo, return 0, indicating that
   * the descriptor could not be injected.
   *
   */
  if ( free_space < DMA_MIN_INJECT_SIZE_IN_QUADS + 
                    DMA_FIFO_DESCRIPTOR_SIZE_IN_QUADS )
    {
      free_space = DMA_FifoGetFreeSpace( &f_ptr->dma_fifo,
					 1,  /* Use hardware head */
					 0); /* Use shadow tail   */

      if ( free_space < DMA_MIN_INJECT_SIZE_IN_QUADS + 
	                DMA_FIFO_DESCRIPTOR_SIZE_IN_QUADS ) return 0;
    }

  /*
   * We have enough room in the fifo.
   * Copy the descriptor to the current va_tail of the fifo.
   * Msync to ensure the descriptor has been written to memory and the L1 caches
   * are in sync.
   * Move the tail past the descriptor so the DMA knows the descriptor is there.
   *   - handle wrapping
   *   - update free space
   *
   */

  if ( ( (unsigned)desc & 0xFFFFFFF0 ) == (unsigned)desc ) /* 16B aligned? */
    {
      load_ptr  = (char*)desc;
      store_ptr = (char*)DMA_FifoGetTailFromShadow( &f_ptr->dma_fifo );
      _bgp_QuadLoad ( load_ptr,     0 );
      _bgp_QuadLoad ( load_ptr+16,  1 );
      _bgp_QuadStore( store_ptr,    0 );
      _bgp_QuadStore( store_ptr+16, 1 );
    }
  else
    {
      memcpy( DMA_FifoGetTailFromShadow( &f_ptr->dma_fifo ),
	      desc,
	      DMA_FIFO_DESCRIPTOR_SIZE_IN_BYTES );
    }

  // _bgp_msync();  mbar is good enough
  _bgp_mbar(); 

  DMA_InjFifoIncrementTail( f_ptr,
			    DMA_FIFO_DESCRIPTOR_SIZE_IN_QUADS );

  return 1; /* Success */
}


/*!
 * \brief Inject Multiple Descriptors into a DMA Injection Fifo
 *
 * Inject multiple descriptors into a DMA injection fifo, given an injection fifo 
 * structure
 *
 * \param[in]  f_ptr     Pointer to the injection fifo structure
 * \param[in]  num_desc  Number of descriptors to be injected
 * \param[in]  desc      A pointer to an array of pointers to descriptors to be 
 *                       injected.  The descriptors must be 16-byte aligned.
 *
 * \retval  numDescInjected  The number of descriptors injected.
 *                           - less than num_desc means some were not injected, 
 *                             most likely because the fifo is full.
 *                           - num_desc means all were successfully injected
 *
 */
__INLINE__ int DMA_InjFifoInjectDescriptors(  
					    DMA_InjFifo_t        *f_ptr,
					    int                   num_desc,
					    DMA_InjDescriptor_t **desc
					   )
{
  unsigned int  free_space;
  unsigned int  total_space_needed_in_quads;
  void         *va_tail;
  void         *va_end;
  void         *va_start;
  char         *target;
  unsigned int  num_quads_to_inject, num_quads_remaining;
  int           i;
  char         *load_ptr, *store_ptr;

  SPI_assert( f_ptr != NULL );
  SPI_assert( desc  != NULL );
  SPI_assert( num_desc > 0  );

  /* Get the free space in the fifo using the shadow value */
  free_space = DMA_FifoGetFreeSpace( &f_ptr->dma_fifo,
				     0, /* Use shadow head */
				     0);/* Use shadow tail */

  total_space_needed_in_quads = num_desc * 
                                  DMA_FIFO_DESCRIPTOR_SIZE_IN_QUADS;

  /* 
   * If after injecting all descriptors (DMA_FIFO_DESCRIPTOR_SIZE_IN_QUADS 
   * per descriptor is the amount we are going to inject), there is still at 
   * least the minimum allowable free space left in the fifo, go ahead and 
   * inject.  We want at least DMA_MIN_INJECT_SIZE_IN_QUADS free space 
   * after injection.
   * 
   * Otherwise, read the hardware head pointer and recalculate the free space,
   * and check again.
   *
   * If there is still not enough room in the fifo for any descriptors, 
   * return 0.  Otherwise, continue and inject as many descriptors as possible.
   *
   */
  if ( free_space < DMA_MIN_INJECT_SIZE_IN_QUADS +
                      total_space_needed_in_quads )
    {
      free_space = DMA_FifoGetFreeSpace( &f_ptr->dma_fifo,
				         1,  /* Use hardware head */
					 0); /* Use shadow tail   */

      if ( free_space < DMA_MIN_INJECT_SIZE_IN_QUADS + 
	                  DMA_FIFO_DESCRIPTOR_SIZE_IN_QUADS ) return 0;
    }

  /*
   * We have enough room in the fifo for at least some descriptors.
   * Copy the descriptors (as many as will fit) to the current va_tail of the 
   * fifo.
   * Msync to ensure the descriptor has been written to memory and the L1 caches
   * are in sync.
   * Move the tail past the descriptor so the DMA knows the descriptor is there.
   *   - handle wrapping
   *   - update free space
   *
   */
  va_tail             = DMA_FifoGetTailFromShadow( &f_ptr->dma_fifo );
  va_start            = DMA_FifoGetStartFromShadow( &f_ptr->dma_fifo );
  va_end              = DMA_FifoGetEndFromShadow( &f_ptr->dma_fifo );
  target              = (char*)va_tail;

  if ( free_space < DMA_MIN_INJECT_SIZE_IN_QUADS + total_space_needed_in_quads ) {
     num_quads_to_inject = free_space - DMA_MIN_INJECT_SIZE_IN_QUADS;
  }
  else { 
    num_quads_to_inject = total_space_needed_in_quads;
  }
  num_quads_remaining = num_quads_to_inject;
  i                   = 0;

  while ( num_quads_remaining > 0 )
    {
      SPI_assert( desc[i] != NULL );

      if ( ( (unsigned)desc[i] & 0xFFFFFFF0 ) == (unsigned)desc[i] ) /* 16B aligned? */
	{
	  load_ptr  = (char*)desc[i];
	  store_ptr = (char*)target;
	  _bgp_QuadLoad ( load_ptr,     0 );
	  _bgp_QuadLoad ( load_ptr+16,  1 );
	  _bgp_QuadStore( store_ptr,    0 );
	  _bgp_QuadStore( store_ptr+16, 1 );
	}
      else
	{
	  memcpy( (char*)target,
		  desc[i],
		  DMA_FIFO_DESCRIPTOR_SIZE_IN_BYTES );
	}

      i++;
      num_quads_remaining -= DMA_FIFO_DESCRIPTOR_SIZE_IN_QUADS;
      target += DMA_FIFO_DESCRIPTOR_SIZE_IN_BYTES;
      if ( target >= (char*)va_end )
	target = (char*)va_start;
    }
  
  // _bgp_msync();  mbar good enough
  _bgp_mbar();

  DMA_InjFifoIncrementTail( f_ptr,
			    num_quads_to_inject );

  return i; /* Success */

}


/*!
 * \brief Get DMA Injection Fifo Group Number
 *
 * Get the DMA Injection Fifo Group number, given an injection fifo group 
 * structure.
 *
 * \param[in]  fg_ptr       Pointer to the structure previously filled in when the
 *                          group was allocated.
 *
 * \return  The DMA Injection Fifo Group number
 *
 */
__INLINE__ int DMA_InjFifoGetGroupNum(
				      const DMA_InjFifoGroup_t *fg_ptr
				     )
{
  SPI_assert( fg_ptr != NULL );

  return fg_ptr->group_id;
}

                  
/*!
 * \brief Get the "Not Empty" Status of an Injection Fifo Group
 *
 * Get the "Not Empty" status of the fifos that the specified fifo group has
 * permission to use.
 *
 * \param[in]  fg_ptr     Pointer to the injection fifo group structure
 *
 * \retval  notEmptyStatus  A 32-bit value, one bit per fifo.
 *                          Bit i is 1 if the specified fifo group has 
 *                          permission to use fifo i and fifo i is not
 *                          empty.
 *                          Bit i is 0 if the specified fifo group either
 *                          does not have permission to use fifo i, or fifo i
 *                          is empty.
 *
 */
__INLINE__ unsigned DMA_InjFifoGetNotEmpty( 
					   DMA_InjFifoGroup_t *fg_ptr
					  ) 
{
  SPI_assert( fg_ptr != NULL );
  SPI_assert( fg_ptr->status_ptr != NULL );

  return ( fg_ptr->status_ptr->not_empty & fg_ptr->permissions );

}

                  
/*!
 * \brief Get the "available" Status of an Injection Fifo Group
 *
 * Get the "available" status of the fifos that the specified fifo group has
 * permission to use.  "available" means that the fifo is enabled and
 * activated.
 *
 * \param[in]  fg_ptr     Pointer to the injection fifo group structure
 *
 * \retval  availableStatus  A 32-bit value, one bit per fifo.
 *                           Bit i is 1 if the specified fifo group has 
 *                           permission to use fifo i and fifo i is available
 *                           Bit i is 0 if the specified fifo group either
 *                           does not have permission to use fifo i, or fifo i
 *                           is not available.
 *
 * \note Normally, there should be a 1 in every position except those that
 *       the specified fifo group does not have permission to use.
 *
 */
__INLINE__ unsigned DMA_InjFifoGetAvailable( 
					    DMA_InjFifoGroup_t *fg_ptr
					   ) 
{
  SPI_assert( fg_ptr != NULL );
  SPI_assert( fg_ptr->status_ptr != NULL );

  return ( fg_ptr->status_ptr->available & fg_ptr->permissions );

}


/*!
 * \brief Get the "threshold crossed" Status of an Injection Fifo Group
 *
 * Get the "threshold crossed" status of the fifos that the specified fifo 
 * group has permission to use.
 *
 * \param[in]  fg_ptr     Pointer to the injection fifo group structure
 *
 * \retval  thresholdCrossedStatus  A 32-bit value, one bit per fifo.
 *                           Bit i is 1 if the specified fifo group has 
 *                           permission to use fifo i and fifo i has crossed
 *                           a threshold.
 *                           Bit i is 0 if the specified fifo group either
 *                           does not have permission to use fifo i, or fifo i
 *                           has not crossed a threshold.
 *
 * \note Normally, there should be a 0 in every position.
 *
 */
__INLINE__ unsigned DMA_InjFifoGetThresholdCrossed( 
						   DMA_InjFifoGroup_t *fg_ptr
						  ) 
{
  SPI_assert( fg_ptr != NULL );
  SPI_assert( fg_ptr->status_ptr != NULL );

  return ( fg_ptr->status_ptr->threshold_crossed & fg_ptr->permissions );

}


/*!
 * \brief Set the "clear threshold crossed" Status of an Injection Fifo Group
 *
 * Set the "clear threshold crossed" status of the fifos that the specified fifo
 * group has permission to use.
 *
 * \param[in]  fg_ptr  Pointer to the injection fifo group structure
 * \param[in]  clr     A 32-bit value, one bit per fifo.
 *                     Only bits that the specified fifo group has
 *                     permission to use should be set to 1.
 *                     Set bit i to 1 to clear the threshold crossed status
 *                     of fifo i.
 *
 * \return  None
 *
 */
__INLINE__ void DMA_InjFifoSetClearThresholdCrossed( 
						    DMA_InjFifoGroup_t *fg_ptr, 
						    unsigned int      clr
						   ) 
{
  SPI_assert( fg_ptr != NULL );
  SPI_assert( fg_ptr->status_ptr != NULL );
  SPI_assert( (clr & fg_ptr->permissions) == clr );

  fg_ptr->status_ptr->clear_threshold_crossed = clr;

}


/*!
 * \brief Get the "activated" Status of an Injection Fifo Group
 *
 * Get the "activated" status of the fifos that the specified fifo 
 * group has permission to use.
 *
 * \param[in]  fg_ptr     Pointer to the injection fifo group structure
 *
 * \retval  activatedStatus  A 32-bit value, one bit per fifo.
 *                           Bit i is 1 if the specified fifo group has 
 *                           permission to use fifo i and fifo i is activated
 *                           Bit i is 0 if the specified fifo group either
 *                           does not have permission to use fifo i, or fifo i
 *                           is not activated.
 *
 */
__INLINE__ unsigned DMA_InjFifoGetActivated( 
					    DMA_InjFifoGroup_t *fg_ptr
					   ) 
{
  SPI_assert( fg_ptr != NULL );
  SPI_assert( fg_ptr->status_ptr != NULL );

  return ( fg_ptr->status_ptr->activated & fg_ptr->permissions );

}


/*!
 * \brief Set the "activate" Status of an Injection Fifo Group
 *
 * Set the "activate" status of the fifos that the specified fifo
 * group has permission to use.
 *
 * \param[in]  fg_ptr  Pointer to the injection fifo group structure
 * \param[in]  act     A 32-bit value, one bit per fifo.
 *                     Only bits that the specified fifo group has
 *                     permission to use should be set to 1.
 *                     Set bit i to 1 to activate fifo i.
 *
 * \return  None
 *
 */
__INLINE__ void DMA_InjFifoSetActivate( 
				       DMA_InjFifoGroup_t *fg_ptr, 
				       unsigned int        act
				      ) 
{
  SPI_assert( fg_ptr != NULL );
  SPI_assert( fg_ptr->status_ptr != NULL );
  SPI_assert( (act & fg_ptr->permissions) == act );

  fg_ptr->status_ptr->activate = act;

}


/*!
 * \brief Set the "deactivate" Status of an Injection Fifo Group
 *
 * Set the "deactivate" status of the fifos that the specified fifo
 * group has permission to use.
 *
 * \param[in]  fg_ptr     Pointer to the injection fifo group structure
 * \param[in]  deact      A 32-bit value, one bit per fifo.
 *                        Only bits that the specified fifo group has
 *                        permission to use should be set to 1.
 *                        Set bit i to 1 to deactivate fifo i.
 *
 * \return  None
 *
 */
__INLINE__ void DMA_InjFifoSetDeactivate( 
					 DMA_InjFifoGroup_t *fg_ptr, 
					 unsigned int        deact
				      ) 
{
  SPI_assert( fg_ptr != NULL );
  SPI_assert( fg_ptr->status_ptr != NULL );
  SPI_assert( (deact & fg_ptr->permissions) == deact );

  fg_ptr->status_ptr->deactivate = deact;

}




/*
 * -----------------------------------------------------------------------------
 * Calls to access the Fifo, given a fifo group and a fifo ID
 * -----------------------------------------------------------------------------
 */




/*!
 * \brief DMA InjFifo Initialization By Id
 *
 * - For an allocated injection DMA fifo, initialize its start, head, tail, and 
 *   end.  
 * - Compute fifo size and free space.
 * - Initialize descriptor count.
 * - Activate the fifo.
 *
 * \param[in]  fg_ptr    Pointer to fifo group structure.
 * \param[in]  fifo_id   Id of the fifo to be initialized
 *                       (0 to DMA_NUM_INJ_FIFOS_PER_GROUP-1).
 * \param[in]  va_start  Virtual address of the start of the fifo.
 * \param[in]  va_head   Virtual address of the head of the fifo (typically 
 *                       equal to va_start).
 * \param[in]  va_end    Virtual address of the end of the fifo.
 *
 * \retval   0  Successful.
 * \retval  -1  Unsuccessful.  Error checks include
 *              - va_start < va_end
 *              - va_start <= va_head <= 
 *                  (va_end - DMA_FIFO_DESCRIPTOR_SIZE_IN_QUADS)
 *              - va_start and va_end are 32-byte aligned
 *              - fifo_size is larger than (DMA_MIN_INJECT_SIZE_IN_QUADS + 
 *                                          DMA_FIFO_DESCRIPTOR_SIZE_IN_QUADS)
 *
 */
__INLINE__ int   DMA_InjFifoInitById(
				     DMA_InjFifoGroup_t *fg_ptr,
				     int                 fifo_id,
				     void               *va_start,
				     void               *va_head,
				     void               *va_end
				    )
{
  int rc;

  SPI_assert( fg_ptr != NULL );
  SPI_assert( fifo_id >= 0 && fifo_id < DMA_NUM_INJ_FIFOS_PER_GROUP );
  SPI_assert( (fg_ptr->permissions & _BN(fifo_id)) != 0 );
  SPI_assert( va_start < va_end );
  SPI_assert( va_start <= va_head );
  SPI_assert( ((uint32_t) va_head)  <= ( ((uint32_t) va_end) - DMA_FIFO_DESCRIPTOR_SIZE_IN_BYTES) );
  SPI_assert( ( ( (uint32_t) va_start) & 0xFFFFFFE0)  == (uint32_t) va_start );
  SPI_assert( ( ( (uint32_t) va_end  ) & 0xFFFFFFE0)  == (uint32_t) va_end );
  SPI_assert( ( (unsigned)va_end - (unsigned)va_start ) >=
	   ( (DMA_MIN_INJECT_SIZE_IN_QUADS + DMA_FIFO_DESCRIPTOR_SIZE_IN_QUADS) * 16 ) );

  /* 
   * Initialize the fifo by invoking a system call.  This system call
   * deactivates the fifo, initializes the start, end, head, and tail,
   * and activates the fifo.
   */
    
   rc = Kernel_InjFifoInitById( 
			       (uint32_t*)fg_ptr,
			       fifo_id,
			       (uint32_t*)va_start,
			       (uint32_t*)va_head,
			       (uint32_t*) va_end
			      );
   
  if (rc) return rc;

  /* Initialize the descriptor count */
  fg_ptr->fifos[fifo_id].desc_count = 0;

  return 0;
}


/*!
 * \brief Get DMA InjFifo Start Pointer from the Shadow Using a Fifo Id
 *
 * Get a DMA injection fifo's start pointer, given a fifo group and fifo id.
 *
 * \param[in]  fg_ptr   Pointer to the fifo group structure
 * \param[in]  fifo_id  Id of the fifo within the group
 *                      (0 to DMA_NUM_INJ_FIFOS_PER_GROUP-1).
 *
 * \retval  va_start  The virtual address of the start of the fifo
 *
 */
__INLINE__ void * DMA_InjFifoGetStartFromShadowById(  
						    DMA_InjFifoGroup_t *fg_ptr,
						    int                 fifo_id
						   )
{
  SPI_assert( (fifo_id >= 0) && (fifo_id < DMA_NUM_INJ_FIFOS_PER_GROUP) );
  SPI_assert( fg_ptr != NULL );
  SPI_assert( ( fg_ptr->permissions & _BN(fifo_id) ) != 0 );

  return DMA_FifoGetStartFromShadow( &fg_ptr->fifos[fifo_id].dma_fifo );
} 


/*!
 * \brief Get DMA InjFifo Head Pointer Using a Fifo Id
 *
 * Get a DMA injection fifo's head pointer, given a fifo group and fifo id.
 *
 * \param[in]  fg_ptr   Pointer to the fifo group structure
 * \param[in]  fifo_id  Id of the fifo within the group
 *                      (0 to DMA_NUM_INJ_FIFOS_PER_GROUP-1).
 *
 * \retval  va_head  The virtual address of the head of the fifo.
 *
 */
__INLINE__ void * DMA_InjFifoGetHeadById(  
					 DMA_InjFifoGroup_t *fg_ptr,
					 int                 fifo_id
					)
{
  SPI_assert( (fifo_id >= 0) && (fifo_id < DMA_NUM_INJ_FIFOS_PER_GROUP) );
  SPI_assert( fg_ptr != NULL );
  SPI_assert( ( fg_ptr->permissions & _BN(fifo_id) ) != 0 );

  return DMA_FifoGetHead( &fg_ptr->fifos[fifo_id].dma_fifo );
} 


/*!
 * \brief Get DMA InjFifo Tail Pointer Using a Fifo Id
 *
 * Get a DMA injection fifo's tail pointer, given a fifo group and fifo id.
 *
 * \param[in]  fg_ptr   Pointer to the fifo group structure
 * \param[in]  fifo_id  Id of the fifo within the group
 *                      (0 to DMA_NUM_INJ_FIFOS_PER_GROUP-1).
 *
 * \retval  va_tail  The virtual address of the tail of the fifo
 *
 */
__INLINE__ void *DMA_InjFifoGetTailById(  
					DMA_InjFifoGroup_t *fg_ptr,
					int                 fifo_id
				       )
{
  SPI_assert( (fifo_id >= 0) && (fifo_id < DMA_NUM_INJ_FIFOS_PER_GROUP) );
  SPI_assert( fg_ptr != NULL );
  SPI_assert( ( fg_ptr->permissions & _BN(fifo_id) ) != 0 );

  return   DMA_FifoGetTail( &fg_ptr->fifos[fifo_id].dma_fifo );
} 


/*!
 * \brief Get DMA InjFifo End Pointer from the Shadow Using a Fifo Id
 *
 * Get a DMA injection fifo's end pointer, given a fifo group and fifo id.
 *
 * \param[in]  fg_ptr   Pointer to the fifo group structure
 * \param[in]  fifo_id  Id of the fifo within the group
 *                      (0 to DMA_NUM_INJ_FIFOS_PER_GROUP-1).
 *
 * \retval  va_end  The virtual address of the end of the fifo
 *
 */
__INLINE__ void *DMA_InjFifoGetEndById(  
				       DMA_InjFifoGroup_t *fg_ptr,
				       int                 fifo_id
				      )
{
  SPI_assert( (fifo_id >= 0) && (fifo_id < DMA_NUM_INJ_FIFOS_PER_GROUP) );
  SPI_assert( fg_ptr != NULL );
  SPI_assert( ( fg_ptr->permissions & _BN(fifo_id) ) != 0 );

  return   DMA_FifoGetEndFromShadow( &fg_ptr->fifos[fifo_id].dma_fifo );
} 


/*!
 * \brief Get DMA InjFifo Size Using a Fifo Id
 *
 * Get a DMA injection fifo's size, given a fifo group and fifo id.
 *
 * \param[in]  fg_ptr   Pointer to the fifo group structure
 * \param[in]  fifo_id  Id of the fifo within the group
 *                      (0 to DMA_NUM_INJ_FIFOS_PER_GROUP-1).
 *
 * \retval  size   The size of the DMA fifo, in units of 16B quads.
 *
 */
__INLINE__ unsigned int DMA_InjFifoGetSizeById(  
					       DMA_InjFifoGroup_t *fg_ptr,
					       int                 fifo_id
					      )
{
  SPI_assert( (fifo_id >= 0) && (fifo_id < DMA_NUM_INJ_FIFOS_PER_GROUP) );
  SPI_assert( fg_ptr != NULL );
  SPI_assert( ( fg_ptr->permissions & _BN(fifo_id) ) != 0 );

  return   DMA_FifoGetSize( &fg_ptr->fifos[fifo_id].dma_fifo );
}


/*!
 * \brief Get DMA InjFifo Free Space Using a Fifo Id
 *
 * Get a DMA injection fifo's free space, given a fifo group and fifo id.
 *
 * \param[in]  fg_ptr   Pointer to the fifo group structure
 * \param[in]  fifo_id  Id of the fifo within the group
 *                      (0 to DMA_NUM_INJ_FIFOS_PER_GROUP-1).
 * \param[in]  read_head  Indicates whether to read the head from the hardware
 *                        fifo before calculating the free space.
 *                          - 1 means to read the hardware head
 *                          - 0 means to use the current head shadow
 * \param[in]  read_tail  Indicates whether to read the tail from the hardware
 *                        fifo before calculating the free space.
 *                          - 1 means to read the hardware tail
 *                          - 0 means to use the current tail shadow
 *
 * \retval  freeSpace  The amount of free space in the fifo, in units of 
 *                     16B quads.
 *
 */
__INLINE__ unsigned int  DMA_InjFifoGetFreeSpaceById(  
					 DMA_InjFifoGroup_t *fg_ptr,
					 int                 fifo_id,
					 unsigned int        read_head, 
					 unsigned int        read_tail
					)
{
  SPI_assert( (fifo_id >= 0) && (fifo_id < DMA_NUM_INJ_FIFOS_PER_GROUP) );
  SPI_assert( fg_ptr != NULL );
  SPI_assert( ( fg_ptr->permissions & _BN(fifo_id) ) != 0 );

  return DMA_FifoGetFreeSpace( &fg_ptr->fifos[fifo_id].dma_fifo,
			       read_head, 
			       read_tail );
} 


/*!
 * \brief Set DMA InjFifo Head Pointer Using a Fifo Id
 *
 * Set a DMA injection fifo's head pointer, given a fifo group and fifo id.
 *
 * \param[in]  fg_ptr   Pointer to the fifo group structure
 * \param[in]  fifo_id  Id of the fifo within the group
 *                      (0 to DMA_NUM_INJ_FIFOS_PER_GROUP-1).
 * \param[in]  va_head  The virtual address of the head of the fifo.
 *
 * \return  None
 *
 */
__INLINE__ void DMA_InjFifoSetHeadById(  
				       DMA_InjFifoGroup_t *fg_ptr,
				       int                 fifo_id,
				       void               *va_head
				       )
{
  SPI_assert( (fifo_id >= 0) && (fifo_id < DMA_NUM_INJ_FIFOS_PER_GROUP) );
  SPI_assert( fg_ptr != NULL );
  SPI_assert( ( fg_ptr->permissions & _BN(fifo_id) ) != 0 );

  DMA_InjFifoSetHead( &fg_ptr->fifos[fifo_id],
		      va_head);
}


/*!
 * \brief Set DMA InjFifo Tail Pointer Using a Fifo Id
 *
 * Set a DMA injection fifo's tail pointer, given a fifo group and fifo id.
 *
 * \param[in]  fg_ptr   Pointer to the fifo group structure
 * \param[in]  fifo_id  Id of the fifo within the group
 *                      (0 to DMA_NUM_INJ_FIFOS_PER_GROUP-1).
 * \param[in]  va_tail  The virtual address of the tail of the fifo.
 *
 * \return  None
 *
 */
__INLINE__ void DMA_InjFifoSetTailById(  
				       DMA_InjFifoGroup_t *fg_ptr,
				       int                 fifo_id, 
				       void               *va_tail 
				      )
{
  SPI_assert( (fifo_id >= 0) && (fifo_id < DMA_NUM_INJ_FIFOS_PER_GROUP) );
  SPI_assert( fg_ptr != NULL );
  SPI_assert( ( fg_ptr->permissions & _BN(fifo_id) ) != 0 );

  DMA_FifoSetTail( &fg_ptr->fifos[fifo_id].dma_fifo,
		   va_tail);  
}


/*!
 * \brief Increment DMA InjFifo Tail Pointer Using a Fifo Id
 *
 * Increment a DMA injection fifo's tail pointer, given a fifo group and 
 * fifo id.
 *
 * \param[in]  fg_ptr   Pointer to the fifo group structure
 * \param[in]  fifo_id  Id of the fifo within the group
 *                      (0 to DMA_NUM_INJ_FIFOS_PER_GROUP-1).
 * \param[in]  incr     The number of quads (16 byte units) to increment the 
 *                      tail pointer by.
 *
 * \return  None
 *
 * \note This function does not check if there is free space in the fifo 
 *       for this many quads.  It must be preceeded by a check of the 
 *       free space.
*/
__INLINE__ void DMA_InjFifoIncrementTailById( 
					     DMA_InjFifoGroup_t *fg_ptr,
					     int                 fifo_id, 
					     unsigned int        incr
					    )
{
  SPI_assert( (fifo_id >= 0) && (fifo_id < DMA_NUM_INJ_FIFOS_PER_GROUP) );
  SPI_assert( fg_ptr != NULL );
  SPI_assert( ( fg_ptr->permissions & _BN(fifo_id) ) != 0 );

  DMA_InjFifoIncrementTail( &fg_ptr->fifos[fifo_id],
			    incr);
} 


/*!
 * \brief Get DMA InjFifo Descriptor Count Using a Fifo Id
 *
 * Get a DMA injection fifo's descriptor count, given a fifo group and 
 * fifo id.
 *
 * \param[in]  fg_ptr   Pointer to the fifo group structure
 * \param[in]  fifo_id  Id of the fifo within the group
 *                      (0 to DMA_NUM_INJ_FIFOS_PER_GROUP-1).
 *
 * \return  None
 *
 */
__INLINE__ unsigned long long DMA_InjFifoGetDescriptorCountById(  
					DMA_InjFifoGroup_t *fg_ptr,
					int                 fifo_id
				       )
{
  SPI_assert( (fifo_id >= 0) && (fifo_id < DMA_NUM_INJ_FIFOS_PER_GROUP) );
  SPI_assert( fg_ptr != NULL );
  SPI_assert( ( fg_ptr->permissions & _BN(fifo_id) ) != 0 );

  return  DMA_InjFifoGetDescriptorCount( &fg_ptr->fifos[fifo_id] );
} 


/*!
 * \brief Is DMA Descriptor Done Using a Fifo Id
 *
 * Return whether a specified descriptor is still in the specified injection
 * fifo (not done).  The descriptor is identified by the descriptor count
 * immediately after the descriptor was injected into the fifo (returned by
 * DMA_InjFifoIncrementTail().
 *
 * \param[in]  fg_ptr      Pointer to the fifo group structure
 * \param[in]  fifo_id     Id of the fifo within the group
 *                         (0 to DMA_NUM_INJ_FIFOS_PER_GROUP-1).
 * \param[in]  desc_count  The descriptor count immediately after the
 *                         descriptor in question was injected into 
 *                         the fifo.
 * \param[in]  update      0  Do not update the fifo's shadow information.
 *                         1  Update the fifo's shadow information.
 *                         It is a performance optimization to only update the
 *                         shadow information once for a group of descriptors
 *                         being processed.
 *
 * \retval  0  False.  The descriptor identified by desc_count is not done.
 *                     It is still in the fifo.
 * \retval  1  True.   The descriptor identified by desc_count is done.
 *                     It is no longer in the fifo.
 *
 */
__INLINE__ unsigned int DMA_InjFifoIsDescriptorDoneById( 
					    DMA_InjFifoGroup_t *fg_ptr,
					    int                 fifo_id,
					    unsigned long long  desc_count,
					    unsigned int        update
					   )
{
  SPI_assert( (fifo_id >= 0) && (fifo_id < DMA_NUM_INJ_FIFOS_PER_GROUP) );
  SPI_assert( fg_ptr != NULL );
  SPI_assert( ( fg_ptr->permissions & _BN(fifo_id) ) != 0 );

  return(DMA_InjFifoIsDescriptorDone( &fg_ptr->fifos[fifo_id],
				      desc_count,
				      update ) );
  
}


/*!
 * \brief DMA Injection Fifo Reserve Descriptor Storage Using a Fifo Id
 *
 * Reserve storage in a DMA injection fifo for a remote get descriptor, given
 * a fifo group and fifo id.
 *
 * \param[in]  fg_ptr   Pointer to the fifo group structure
 * \param[in]  fifo_id  Id of the fifo within the group
 *                      (0 to DMA_NUM_INJ_FIFOS_PER_GROUP-1).
 *
 * \retval  0   Successful.  There was enough space in the fifo and the
 *              storage was reserved.
 * \retval -1   Unsuccessful.  There was not enough space in the fifo.
 *
 * \note Internally, this increments the occupiedSize of the fifo by
 * DMA_FIFO_DESCRIPTOR_SIZE_IN_QUADS. 
 *
 */
__INLINE__ int DMA_InjFifoReserveDescriptorStorageById( 
					DMA_InjFifoGroup_t *fg_ptr,
					int                 fifo_id
				       )
{
  SPI_assert( (fifo_id >= 0) && (fifo_id < DMA_NUM_INJ_FIFOS_PER_GROUP) );
  SPI_assert( fg_ptr != NULL );
  SPI_assert( ( fg_ptr->permissions & _BN(fifo_id) ) != 0 );

  return ( DMA_InjFifoReserveDescriptorStorage( &fg_ptr->fifos[fifo_id] ) );
} 


/*!
 * \brief DMA Injection Fifo Free Descriptor Storage Reservation Using a Fifo Id
 *
 * Free a reservation for storage for a remote get descriptor in a DMA injection
 * fifo, previously reserved using DMA_InjFifoReserveDescriptorStorageById().
 *
 * \param[in]  fg_ptr   Pointer to the fifo group structure
 * \param[in]  fifo_id  Id of the fifo within the group
 *                      (0 to DMA_NUM_INJ_FIFOS_PER_GROUP-1).
 *
 * \return  None
 *
 * \note Internally, this decrements the occupiedSize of the fifo by
 * DMA_FIFO_DESCRIPTOR_SIZE_IN_QUADS. 
 *
 */
__INLINE__ void DMA_InjFifoFreeDescriptorStorageReservationById( 
					DMA_InjFifoGroup_t *fg_ptr,
					int                 fifo_id
						)
{
  SPI_assert( (fifo_id >= 0) && (fifo_id < DMA_NUM_INJ_FIFOS_PER_GROUP) );
  SPI_assert( fg_ptr != NULL );
  SPI_assert( ( fg_ptr->permissions & _BN(fifo_id) ) != 0 );

  DMA_InjFifoFreeDescriptorStorageReservation( &fg_ptr->fifos[fifo_id] );
  return;
} 


/*!
 * \brief Check If An Injection Fifo Has Space For Injection Using a Fifo Id
 *
 * Check if an injection fifo has enough space for a single descriptor to be
 * injected, given a fifo group and fifo id.
 *
 * \param[in]  fg_ptr   Pointer to the fifo group structure
 * \param[in]  fifo_id  Id of the fifo within the group
 *                      (0 to DMA_NUM_INJ_FIFOS_PER_GROUP-1).
 *
 * \retval  hasSpace  An indicator of whether the fifo has space for a 
 *                    descriptor.
 *                    - 0 (false) means the fifo is full.
 *                    - 1 (true)  means the fifo has space.
 *
 */
__INLINE__ unsigned int DMA_InjFifoHasSpaceById(  
						DMA_InjFifoGroup_t    *fg_ptr,
						int                    fifo_id
					       )
{
  SPI_assert( (fifo_id >= 0) && (fifo_id < DMA_NUM_INJ_FIFOS_PER_GROUP) );
  SPI_assert( fg_ptr != NULL );
  SPI_assert( ( fg_ptr->permissions & _BN(fifo_id) ) != 0 );

  return  DMA_InjFifoHasSpace( &fg_ptr->fifos[fifo_id] );
}


/*!
 * \brief Inject a Descriptor into a DMA Injection Fifo Without Checking for 
 *        Space Using a Fifo Id
 *
 * Inject a descriptor into a DMA injection fifo, given a fifo group and 
 * fifo id, without checking to see if there is enough space in the fifo.
 * It is assumed that the caller has already checked for enough space using
 * the DMA_InjFifoHasSpaceById() function.
 *
 * \param[in]  fg_ptr   Pointer to the fifo group structure
 * \param[in]  fifo_id  Id of the fifo within the group
 *                      (0 to DMA_NUM_INJ_FIFOS_PER_GROUP-1).
 * \param[in]  desc   A pointer to the descriptor to be injected.
 *                    Must be 16-byte aligned.
 *
 * \retval  numDescInjected  The number of descriptors injected.
 *                           - 1 means it was successfully injected.
 *
 * \see DMA_InjFifoHasSpaceById()
 */
__INLINE__ int DMA_InjFifoInjectDescriptorNoSpaceCheckById(  
					DMA_InjFifoGroup_t    *fg_ptr,
					int                    fifo_id,
					DMA_InjDescriptor_t   *desc
				       )
{
  SPI_assert( (fifo_id >= 0) && (fifo_id < DMA_NUM_INJ_FIFOS_PER_GROUP) );
  SPI_assert( fg_ptr != NULL );
  SPI_assert( ( fg_ptr->permissions & _BN(fifo_id) ) != 0 );

  return  DMA_InjFifoInjectDescriptorNoSpaceCheck( &fg_ptr->fifos[fifo_id],
						   desc );
}


/*!
 * \brief Inject Descriptor into a DMA InjFifo Using a Fifo Id
 *
 * Inject a descriptor into a DMA injection fifo, given a fifo group and 
 * fifo id.
 *
 * \param[in]  fg_ptr   Pointer to the fifo group structure
 * \param[in]  fifo_id  Id of the fifo within the group
 *                      (0 to DMA_NUM_INJ_FIFOS_PER_GROUP-1).
 * \param[in]  desc     Pointer to the descriptor to be injected into the fifo.
 *
 * \retval  numDescInjected  The number of descriptors injected.
 *                           - 0 means it was not injected, most likely because
 *                             the fifo is full.
 *                           - 1 means it was successfully injected
 *
 */
__INLINE__ int DMA_InjFifoInjectDescriptorById( 
					DMA_InjFifoGroup_t    *fg_ptr,
					int                    fifo_id,
					DMA_InjDescriptor_t   *desc  
				       )
{
  SPI_assert( (fifo_id >= 0) && (fifo_id < DMA_NUM_INJ_FIFOS_PER_GROUP) );
  SPI_assert( fg_ptr != NULL );
  SPI_assert( ( fg_ptr->permissions & _BN(fifo_id) ) != 0 );

  return  DMA_InjFifoInjectDescriptor( &fg_ptr->fifos[fifo_id],
				       desc );
}


/*!
 * \brief Inject Multiple Descriptors into a DMA InjFifo Using a Fifo Id
 *
 * Inject multiple descriptors into a DMA injection fifo, given a fifo group and
 * fifo id.
 *
 * \param[in]  fg_ptr    Pointer to the fifo group structure
 * \param[in]  fifo_id   Id of the fifo within the group
 *                       (0 to DMA_NUM_INJ_FIFOS_PER_GROUP-1).
 * \param[in]  num_desc  Number of descriptors to be injected 
 * \param[in]  desc      Pointer to an array of pointers to the descriptors to 
 *                       be injected into the fifo.
 *
 * \retval  numDescInjected  The number of descriptors injected.
 *                           - less than num_desc means that some were not 
 *                             injected, most likely because the fifo is full.
 *                           - equal to num_desc means that all were 
 *                             successfully injected.
 *
 */
__INLINE__ int DMA_InjFifoInjectDescriptorsById(  
					DMA_InjFifoGroup_t     *fg_ptr,
                                        int                     fifo_id,
                                        int                     num_desc,
                                        DMA_InjDescriptor_t   **desc
				       )
{
  SPI_assert( (fifo_id >= 0) && (fifo_id < DMA_NUM_INJ_FIFOS_PER_GROUP) );
  SPI_assert( fg_ptr != NULL );
  SPI_assert( ( fg_ptr->permissions & _BN(fifo_id) ) != 0 );

  return  DMA_InjFifoInjectDescriptors ( &fg_ptr->fifos[fifo_id],
					 num_desc,
					 desc );
}


/*!
 * \brief Get DMA InjFifo Not Empty Status Using a Fifo Id
 *
 * Get a DMA injection fifo's not empty status, given a fifo group and 
 * fifo id.
 *
 * \param[in]  fg_ptr   Pointer to the fifo group structure
 * \param[in]  fifo_id  Id of the fifo within the group
 *                      (0 to DMA_NUM_INJ_FIFOS_PER_GROUP-1).
 *
 * \return  32-bit status.  status bit "fifo_id" is 1 if the 
 *          fifo is not empty, 0 if empty.  That is, the return value
 *          is 0 if empty, non-zero if not empty.
 *
 */
__INLINE__ unsigned DMA_InjFifoGetNotEmptyById( 
					DMA_InjFifoGroup_t *fg_ptr,
					int                 fifo_id
				       ) 
{
  SPI_assert( (fifo_id >= 0) && (fifo_id < DMA_NUM_INJ_FIFOS_PER_GROUP) );
  SPI_assert( fg_ptr != NULL );
  SPI_assert( ( fg_ptr->permissions & _BN(fifo_id) ) != 0 );

  return ( DMA_InjFifoGetNotEmpty( fg_ptr ) & _BN(fifo_id) );
}


/*!
 * \brief Get DMA InjFifo Available Status Using a Fifo Id
 *
 * Get a DMA injection fifo's available status, given a fifo group and 
 * fifo id.
 *
 * \param[in]  fg_ptr   Pointer to the fifo group structure
 * \param[in]  fifo_id  Id of the fifo within the group
 *                      (0 to DMA_NUM_INJ_FIFOS_PER_GROUP-1).
 *
 * \return  32-bit status.  status bit fifo_id is 1 if the 
 *          fifo is available, 0 if empty.
 *
 * \note  "available" status means the fifo is enabled and activated.
 *
 */
__INLINE__ unsigned DMA_InjFifoGetAvailableById( 
					DMA_InjFifoGroup_t *fg_ptr,
					int                 fifo_id
				       ) 
{
  SPI_assert( (fifo_id >= 0) && (fifo_id < DMA_NUM_INJ_FIFOS_PER_GROUP) );
  SPI_assert( fg_ptr != NULL );
  SPI_assert( ( fg_ptr->permissions & _BN(fifo_id) ) != 0 );

  return ( DMA_InjFifoGetAvailable( fg_ptr ) & _BN(fifo_id) );
}


/*!
 * \brief Get DMA InjFifo Threshold Crossed Status Using a Fifo Id
 *
 * Get a DMA injection fifo's threshold crossed status, given a fifo group and
 * fifo id.
 *
 * \param[in]  fg_ptr   Pointer to the fifo group structure
 * \param[in]  fifo_id  Id of the fifo within the group
 *                      (0 to DMA_NUM_INJ_FIFOS_PER_GROUP-1).
 *
 * \return  32-bit status.  status bit fifo_id is 1 if the 
 *          fifo's threshold has been crossed, 0 if not.
 *
 * \note  This will always be zero.
 *
 */
__INLINE__ unsigned DMA_InjFifoGetThresholdCrossedById( 
					DMA_InjFifoGroup_t *fg_ptr,
					int                 fifo_id
				       ) 
{
  SPI_assert( (fifo_id >= 0) && (fifo_id < DMA_NUM_INJ_FIFOS_PER_GROUP) );
  SPI_assert( fg_ptr != NULL );
  SPI_assert( ( fg_ptr->permissions & _BN(fifo_id) ) != 0 );

  return ( DMA_InjFifoGetThresholdCrossed( fg_ptr ) & _BN(fifo_id) );
}


/*!
 * \brief Clear DMA InjFifo Threshold Crossed Status Using a Fifo Id
 *
 * Clear a DMA injection fifo's threshold crossed status, given a fifo group and
 * fifo id.
 *
 * \param[in]  fg_ptr   Pointer to the fifo group structure
 * \param[in]  fifo_id  Id of the fifo within the group
 *                      (0 to DMA_NUM_INJ_FIFOS_PER_GROUP-1).
 *
 * \return  None
 *
 */
__INLINE__ void DMA_InjFifoSetClearThresholdCrossedById( 
					DMA_InjFifoGroup_t *fg_ptr,
					int                 fifo_id
				       ) 
{
  SPI_assert( (fifo_id >= 0) && (fifo_id < DMA_NUM_INJ_FIFOS_PER_GROUP) );
  SPI_assert( fg_ptr != NULL );
  SPI_assert( ( fg_ptr->permissions & _BN(fifo_id) ) != 0 );

  DMA_InjFifoSetClearThresholdCrossed( fg_ptr,
				       _BN(fifo_id) );
}


/*!
 * \brief Get DMA InjFifo Activated Status Using a Fifo Id
 *
 * Get a DMA injection fifo's activated status, given a fifo group and 
 * fifo id.
 *
 * \param[in]  fg_ptr   Pointer to the fifo group structure
 * \param[in]  fifo_id  Id of the fifo within the group
 *                      (0 to DMA_NUM_INJ_FIFOS_PER_GROUP-1).
 *
 * \return  32-bit status.  status bit fifo_id is 1 if the 
 *          fifo is activated, 0 if empty.
 *
 */
__INLINE__ unsigned DMA_InjFifoGetActivatedById( 
				DMA_InjFifoGroup_t *fg_ptr,
				int                 fifo_id
			       ) 
{
  SPI_assert( (fifo_id >= 0) && (fifo_id < DMA_NUM_INJ_FIFOS_PER_GROUP) );
  SPI_assert( fg_ptr != NULL );
  SPI_assert( ( fg_ptr->permissions & _BN(fifo_id) ) != 0 );

  return ( DMA_InjFifoGetActivated( fg_ptr ) & _BN(fifo_id) );
}


/*!
 * \brief Activate DMA InjFifo Using a Fifo Id
 *
 * Activate a DMA injection fifo, given a fifo group and fifo id.
 *
 * \param[in]  fg_ptr   Pointer to the fifo group structure
 * \param[in]  fifo_id  Id of the fifo within the group
 *                      (0 to DMA_NUM_INJ_FIFOS_PER_GROUP-1).
 *
 * \return  None
 *
 */
__INLINE__ void DMA_InjFifoSetActivateById( 
				DMA_InjFifoGroup_t *fg_ptr,
				int                 fifo_id
			       ) 
{
  SPI_assert( (fifo_id >= 0) && (fifo_id < DMA_NUM_INJ_FIFOS_PER_GROUP) );
  SPI_assert( fg_ptr != NULL );
  SPI_assert( ( fg_ptr->permissions & _BN(fifo_id) ) != 0 );

  DMA_InjFifoSetActivate( fg_ptr,
			  _BN(fifo_id) );
}


/*!
 * \brief Deactivate DMA InjFifo Using a Fifo Id
 *
 * Deactivate a DMA injection fifo, given a fifo group and fifo id.
 *
 * \param[in]  fg_ptr   Pointer to the fifo group structure
 * \param[in]  fifo_id  Id of the fifo within the group
 *                      (0 to DMA_NUM_INJ_FIFOS_PER_GROUP-1).
 *
 * \return  None
 *
 */
__INLINE__ void DMA_InjFifoSetDeactivateById( 
					     DMA_InjFifoGroup_t *fg_ptr,
					     int                 fifo_id
					    ) 
{
  SPI_assert( (fifo_id >= 0) && (fifo_id < DMA_NUM_INJ_FIFOS_PER_GROUP) );
  SPI_assert( fg_ptr != NULL );
  SPI_assert( ( fg_ptr->permissions & _BN(fifo_id) ) != 0 );

  DMA_InjFifoSetDeactivate( fg_ptr,
			    _BN(fifo_id) );
}


__END_DECLS


#endif 
