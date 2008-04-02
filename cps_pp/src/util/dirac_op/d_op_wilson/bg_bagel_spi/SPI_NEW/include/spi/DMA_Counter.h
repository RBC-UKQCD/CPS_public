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

#ifndef	_DMA_COUNTER_H_ /* Prevent multiple inclusion */
#define	_DMA_COUNTER_H_


/*!
 * \file spi/DMA_Counter.h
 * 
 * \brief DMA SPI Counter Definitions and Inline Functions
 *
 * This include file contains inline functions that are used to interface with
 * BG/P DMA injection and reception counters at the lowest level.
 * Functions include
 * - set and get a counter's value and base address
 * - enable and disable a counter or group of counters
 * - query whether a counter or group of counters has hit zero
 * - clear a counter's or group of counters' hit-zero state
 * - set and get a reception counter's maximum address 
 *
 * Definitions:
 * - A counter is a 32-bit value containing the number of bytes being 
 *   transferred from/to memory
 * - Associated with a counter is a base address indicating where the data is
 *   being transferred from/to
 * - Associated with a reception counter is a max address bounding the DMA
 *   transfer.
 * - There are injection (iDMA) and reception (rDMA) counters
 * - There are DMA_NUM_COUNTERS iDMA counters and DMA_NUM_COUNTERS rDMA 
 *   counters
 * - A counter group consists of DMA_NUM_COUNTERS_PER_GROUP counters
 * - There are DMA_NUM_COUNTER_GROUPS iDMA counter groups and 
 *   DMA_NUM_COUNTER_GROUPS rDMA counter groups
 * - A subgroup consists of DMA_NUM_COUNTERS_PER_SUBGROUP counters.  This is 
 *   the unit of counter allocation.
 * - The highest-level counter inlines in this include file work with virtual
 *   addresses.  They are converted to physical addresses and placed into the
 *   counter.
 * - The counter's base and max addresses reside in the DMA memory map (DMA 
 *   SRAM).  The DMA_CounterHw_t structure, known as the hardware counter
 *   structure maps a single counter in this storage.  They are "shadowed" by 
 *   these inline functions to a DMA_Counter_t structure in DDR memory, 
 *   known as the software counter structure, and their associated virtual 
 *   address is also stored in that structure for easy retrieval.  The 
 *   physical addresses really don't have to reside in this shadow structure, 
 *   but it is faster to access them there than from the DMA's SRAM.
 * - The counter's base and max addresses are stored in the DMA SRAM as
 *   16B-aligned 4-bit shifted physical addresses.  That is, the 36-bit
 *   physical address is right shifted 4 bits, aligning it on a 16B boundary
 *   leaving 32 bits.  The following naming conventions are used to store
 *   addresses:
 *   - pa_xxxx: Physical address (32-bit, 16B-aligned 4-bit shifted)
 *   - va_xxxx: Virtual address (32 bits).
 * 
 * \verbatim Picture of data structures:
   
   ========DDR MEMORY===================|==========DMA SRAM MEMORY=============
   ------------------------------       |
   | DMA_CounterGroup_t         |       |
   |                            |       |     --------------------------------
   | status --------------------|-------|---->| DMA_CounterStatus_t          |
   | counter[0..63]             |       |     --------------------------------
   |   ------------------------ |       | 
   |   | DMA_Counter_t        | |       |     ----------------------------- 
   | 0 | (software counter)   | |       |     | DMA_CounterHw_t           |
   |   | counter_hw_ptr-------|-|-------|---->| (hardware counter)        |
   |   ------------------------ |       |     -----------------------------
   |             .              |       |
   |             .              |       |
   |             .              |       |
   |   ------------------------ |       | 
   |   | DMA_Counter_t        | |       |     -----------------------------
   |63 | (software counter)   | |       |     | DMA_CounterHw_t           |
   |   | counter_hw_ptr-------|-|-------|---->| (hardware counter)        |
   |   ------------------------ |       |     -----------------------------
   |             .              |       |
   ------------------------------       |
  
   \endverbatim
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
 */


#include <common/namespace.h>


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
 * \verbatim
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
 * \verbatim
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


#include <errno.h>
#include <bpcore/ppc450_inlines.h>  /* For bgp_msync_nonspeculative() */
#include <spi/bpcore_interface.h>   /* For _BGP_IC_DMA_NFT_G3_HIER_POS*/
#include <spi/kernel_interface.h>   /* For Kernel_Virtual2Physical()  */
#include <spi/DMA_Assert.h>


/*
 * ------------------------------------------------------------------------------
 * Definitions
 * ------------------------------------------------------------------------------
 */

/*!
 * \brief Number of DMA counter groups
 *
 * There are 4 counter groups.
 *
 */
#define DMA_NUM_COUNTER_GROUPS     4


/*!
 * \brief Number of DMA counters in a counter group
 *
 * There are 64 counters in a counter group.
 *
 */
#define DMA_NUM_COUNTERS_PER_GROUP 64


/*!
 * \brief Number of DMA counters in a counter subgroup
 *
 * There are 8 counters in a counter subgroup.
 *
 */
#define DMA_NUM_COUNTERS_PER_SUBGROUP 8


/*!
 * \brief Number of DMA counter subgroups in a group
 *
 * There are 8 subgroups in a counter group.
 *
 */
#define DMA_NUM_COUNTER_SUBGROUPS_PER_GROUP (DMA_NUM_COUNTERS_PER_GROUP / DMA_NUM_COUNTERS_PER_SUBGROUP)


/*!
 * \brief Number of DMA counter subgroups, in total, across all groups
 *
 * There are 32 subgroups in total.
 *
 */
#define DMA_NUM_COUNTER_SUBGROUPS (DMA_NUM_COUNTER_SUBGROUPS_PER_GROUP * DMA_NUM_COUNTER_GROUPS)


/*!
 * \brief Initial value for a DMA counter
 *
 * This value is somewhat arbitrary, but is chosen to be different from zero,
 * because zero means the counter has hit zero, and may cause false interupts.
 *
 */
#define DMA_NUM_COUNTERS ( DMA_NUM_COUNTER_GROUPS * DMA_NUM_COUNTERS_PER_GROUP) 


/*!
 * \brief Initial value for a DMA counter
 *
 * This value is somewhat arbitrary, but is chosen to be different from zero,
 * because zero means the counter has hit zero, and may cause false interupts.
 *
 */
#define DMA_COUNTER_INIT_VAL 0xFFFFFFFF


/*!
 * \brief Max Number of Cores Per Node
 *
 * This is the maximum number of cores that can run on a compute node.
 */
#define DMA_MAX_NUM_CORES 4


/*!
 * \brief Returns the word number that the specified counter is in
 *
 * \param[in]  counter_id  The ID of the counter (0 to 
 *                         DMA_NUM_COUNTERS_PER_GROUP-1)
 * 
 * \return The number of the word that the specified counter is in (0 or 1)
 *
 * Used as an index in the "enabled", "enable", "disable", "hit_zero", and
 * "clear_hit_zero" fields of the DMA_CounterStatus_t structure, and
 * the permissions field of the DMA_CounterGroup_t structure.
 *
 */
#define DMA_COUNTER_GROUP_WORD_ID(counter_id) ((counter_id)>>5)


/*!
 * \brief Returns the bit within the word that the specified counter is in
 *
 * \param[in]  counter_id  The ID of the counter (0 to 
 *                         DMA_NUM_COUNTERS_PER_GROUP-1)
 *
 * \return The bit position within the word that the specified counter is
 *         in (0-31)
 *
 * Used with the "enabled", "enable", "disable", "hit_zero", and
 * "clear_hit_zero" fields of the DMA_CounterStatus_t structure, and
 * the permissions" field of the DMA_CounterGroup_t structure.
 *
 */
#define DMA_COUNTER_GROUP_WORD_BIT_ID(counter_id) ((counter_id) & 0x0000001F)


/* 
 * -----------------------------------------------------------------------------
 * Structures
 * -----------------------------------------------------------------------------
 */

/*!
 * \brief Hardware DMA Counter
 *
 * This maps a DMA counter as it is in the DMA memory map (DMA SRAM).
 *
 */
typedef struct DMA_CounterHw_t
{
  volatile unsigned  counter;   /*!< RW Value of the counter                */
  volatile unsigned  increment; /*!< W Increment the counter by this value  */
  volatile unsigned  pa_base;   /*!< RW Base address of the counter, 32 bit
                                        16B-aligned 4-bit shifted address   */
  volatile unsigned  pa_max;    /*!< RW Maximum payload address (rDMA only),
				     16B-aligned 4-bit shifted address      */ 
} 
DMA_CounterHw_t;


/*!
 * \brief DMA Counter Hardware Status structure
 *
 * This structure maps the DMA SRAM for a particular group of 
 * DMA_NUM_COUNTERS_PER_GROUP counters.
 *
 * This is a common structure between iDMA and rDMA.
 * 
 * \see DMA_COUNTER_GROUP_WORD_ID
 * \see DMA_COUNTER_GROUP_WORD_BIT_ID
 * 
 */
typedef struct DMA_CounterStatus_t
{
  volatile unsigned enabled[2];        /*!< R bitmask (1 bit/counter):
                                              Counter is enabled (1=enabled)   */
  volatile unsigned enable[2];         /*!< W bitmask (1 bit/counter):         
                                              Counter enable: writing a 1 to   
                                              bit i enables counter i.  This   
                                              changes the corrresponding bit   
                                              in enabled.                      */
  volatile unsigned disable[2];        /*!< W bitmask (1 bit/counter):         
                                              Counter disble: writing a 1 to   
                                              bit i disbles counter i.  This   
                                              changes the corrresponding bit   
                                              in enabled.                      */
  volatile unsigned reserved[2];       /*!< HOLE                               */
  volatile unsigned hit_zero[2];       /*!< R bitmask (1 bit/counter):         
                                              Counter hit zero                 
                                              (1=counter hit zero)             */
  volatile unsigned clear_hit_zero[2]; /*!< W bitmask (1 bit/counter):         
                                              Clear counter hit zero: writing  
                                              a 1 to bit i clears the          
                                              corresponding bit in hit_zero.   */
  volatile unsigned grp_status;        /*!< R bitmask (1 bit/subgroup):        
                                              bit i is 1 if or-reduce over     
                                              sub-group i of the hit_zero bits 
                                              anded with the enable bits.      
                                              Note this includes info about    
                                              all DMA_NUM_COUNTERS counters,   
                                              not just those in this group.    */
} 
DMA_CounterStatus_t;


/*!
 * \brief Software DMA Counter Structure
 *
 * This structure provides a shadow (recent copy) of the hardware counter's 
 * base and max.  While accessing the actual hardware DMA counter's base and
 * max is equivalent, it is slower than accessing them from here.
 *
 * Additionally, it stores the corresponding virtual addresses, for easy
 * retrieval, since the hardware counter does not maintain the virtual
 * address.
 *
 * Finally, it contains a pointer to the corresponding hardware counter in 
 * DMA SRAM.
 *
 */
typedef struct DMA_Counter_t 
{
  void         *va_base;    /*!< Shadow virtual address of the base            */
  unsigned int  pa_base;    /*!< Shadow physical address of the base.          
                                 16B-aligned 4-bit shifted address.            */
  void         *va_max;     /*!< Shadow virtual address of the max (rDMA only) */  
  unsigned int  pa_max;     /*!< Shadow physical address of the max (rDMA only)
                                 16B-aligned 4-bit shifted address.            */
  DMA_CounterHw_t *counter_hw_ptr; /*!< Pointer to the hardware counter        */
}
ALIGN_L1D_CACHE  DMA_Counter_t; 
/*!
 * \todo  Re-think whether we need to align this structure on a L1 cache line boundary 
 *
 */


/*!
 * \enum DMA_Type_t
 * \brief DMA type (injection/reception) enum
 *
 */
typedef enum DMA_Type_t 
{
  DMA_Type_Injection = 0,  /*!< Injection type of DMA */
  DMA_Type_Reception = 1   /*!< Reception type of DMA */

} 
DMA_Type_t;


/*!
 * \brief DMA Counter Group Structure
 *
 * This structure defines a DMA Counter Group.  It is filled in by the kernel
 * during the DMA_CounterGroupAllocate system call.  It points to a
 * DMA Counter Status structure, and contains up to DMA_NUM_COUNTERS_PER_GROUP
 * software DMA Counter structures making up this group.  
 *
 * It also contains permission bits to use the counters, one bit per counter.
 * When the permission bit is on, the corresponding counter belongs to this
 * group and can be used.  Otherwise, the counter should not be used as part
 * of this group.  These permission bits are used as follows:
 *   1. Inline functions will ASSERT when an attempt is made
 *      to use a counter that is not part of this group.
 *   2. Inline functions will use the permission bits as a mask 
 *      to return status information only for the counters allocated
 *      to this group.
 * Use the DMA_COUNTER_GROUP_WORD_ID and DMA_COUNTER_GROUP_WORD_BIT_ID
 * macros to locate the appropriate "permitted_counters" bit.
 *
 * Allocations are done in subgroups (groups of DMA_NUM_COUNTERS_PER_SUBGROUP 
 * counters).  This structure contains a bit mask of the subgroups that belong 
 * to this group.
 * 
 * \see DMA_COUNTER_GROUP_WORD_ID
 * \see DMA_COUNTER_GROUP_WORD_BIT_ID
 *
 */
typedef struct DMA_CounterGroup_t
{

  DMA_CounterStatus_t *status_ptr;        /*!< Pointer to counter status       */
  unsigned int permissions[2];            /*!< Bit i is 1 if permitted to use  
                                               counter i, 0 otherwise. One bit 
                                               per counter,                    
                                               DMA_NUM_COUNTERS_PER_GROUP      
                                               counters.                       */
  unsigned int grp_permissions;           /*!< Bit i is 1 if permitted to use  
                                               subgroup i, 0 otherwise. One    
                                               bit per subgroup, 8 subgroups.  */
  unsigned int group_id;                  /*!< The id of this group (0 to        
					       DMA_NUM_COUNTER_GROUPS-1).      */
  DMA_Type_t type;                        /*!< The type of the DMA (injection  
                                               or reception)                   */
  DMA_Counter_t counter[DMA_NUM_COUNTERS_PER_GROUP]; /*!<                     
                                               Software Counter Structures.    
                                               i-th structure's hardware       
                                               pointer is non-NULL if          
                                               permissions[i]=1, NULL if       
                                               permissions[i]=0.               */
} 
DMA_CounterGroup_t;


/*!
 *
 * \brief Counter Application Segment
 *
 * A segment of user-addressible memory.
 * Each segment consists of a virtual address, physical address, and length 
 * defining a contiguous segment of storage that is accessible from the
 * application.
 */
typedef struct DMA_CounterAppSegment_t
{
  unsigned int length;     /*!< Length in bytes of the segment                */
  uint32_t     va_base;    /*!< Virtual address of the segment base           */
  uint32_t     pa_base;    /*!< Shifted physical address of the segment base  */
  uint32_t     va_max;     /*!< Virtual address of the last byte of segment   */
} DMA_CounterAppSegment_t;


/*!
 * 
 * \brief Counter Application Segments
 *
 * An array of application segments.  There are N application segments per core 
 * on a node.  Thus there are N * (number of cores on a node) application 
 * segments in this array.  The first group of segments in the array correspond
 * to core 0.  The second group, core 1, etc.  
 */
extern DMA_CounterAppSegment_t *DMA_CounterAppSegmentArray;


/*!
 * \brief Number of application segments for a core
 *     
 * The number of application segments is the same for all cores.
 */
extern uint32_t                  DMA_CounterNumAppSegments;


/*!
 * \brief The index of the last application segment accessed for a core.
 */
extern int                       DMA_CounterCachedAppSegmentIndex[DMA_MAX_NUM_CORES];


/*!
 * \brief The Minimum 4-bit Shifted Physical Address Accessible From User Mode
 */
extern uint32_t                  DMA_CounterMinPaAccessibleFromUserMode[DMA_MAX_NUM_CORES];

/*!
 *
 * \brief Initialize Counter Application Segments
 *
 * Initialize the array of application segments and the global pointer to it.
 * This identifies the memory regions that the application can access.
 * 
 * Also, initialize the minimum physical address accessible from user mode 
 * for each core.
 *
 * \retval  0            Success
 * \retval  errorNumber  Failure
 */
int DMA_CounterInitAppSegments(void);


/*!
 *
 * \brief Get Number of Counter Application Segments
 *
 * \returns  Number of application segments for a core.
 */
__INLINE__ uint32_t DMA_CounterGetNumAppSegments( void )
{
  return ( DMA_CounterNumAppSegments );
}


/*!
 *
 * \brief Get Pointer to Counter Application Segments
 *
 * \param[in]  Core number whose application segments pointer is to be
 *             returned.
 *
 * \returns  Pointer to application segments
 */
__INLINE__ DMA_CounterAppSegment_t * DMA_CounterGetAppSegments( unsigned int coreNum )
{
  SPI_assert ( coreNum >= 0 );
  SPI_assert ( coreNum <= DMA_MAX_NUM_CORES );
  
  unsigned int index = coreNum * DMA_CounterGetNumAppSegments();
  return ( & ( DMA_CounterAppSegmentArray [ index ] ) );
}


/*!
 *
 * \brief Get Virtual Addresses for the Min and Max Physical Addresses 
 *        for User Space
 *
 * Based on information in the DMA_CounterAppSegments array, return the
 * virtual addresses associated with the min and max physical addresses
 * allowed for user space.
 *
 * \param[out]  va_min  Pointer to a pointer.  Upon return, the pointer is 
 *                      set to the virtual address associated with the
 *                      minimum physical address allowed for user space.
 * \param[out]  va_max  Pointer to a pointer.  Upon return, the pointer is 
 *                      set to the virtual address associated with the
 *                      maximum physical address allowed for user space.
 *
 * If the DMA_CounterNumAppSegments array has not been initialized yet
 * (it is initialized in DMA_CounterGroupAllocate()), a value of 0 for the
 * min and 0xFFFFFFFF max is returned.
 */
__INLINE__ void DMA_CounterGetMinMaxVa(void ** va_min,
				       void ** va_max)
{
  /* Determine the core we are running on so the correct application 
   * segments are consulted 
   */
  unsigned int              coreNum         = Kernel_PhysicalProcessorID();
  DMA_CounterAppSegment_t * appSegmentArray = DMA_CounterGetAppSegments(coreNum);
  uint32_t                  numAppSegments  = DMA_CounterGetNumAppSegments();

  if ( appSegmentArray )
  {
    uint32_t minPaBase=0xFFFFFFFF, maxPa=0;
    uint32_t segmentPaBase, segmentPaMax;
    uint32_t i, minIndex=0, maxIndex=0;

    for (i=0; i<numAppSegments; i++)
    {
      segmentPaBase = appSegmentArray[i].pa_base;
      if ( segmentPaBase < minPaBase )
      {
	minPaBase = segmentPaBase;
	minIndex  = i;
      }

      segmentPaMax = appSegmentArray[i].pa_base + (appSegmentArray[i].length >> 4);
      if ( segmentPaMax > maxPa )
      {
	maxPa     = segmentPaMax;
	maxIndex  = i;
      }
    }

    *va_min = (void*)(appSegmentArray[minIndex].va_base);
    *va_max = (void*)(appSegmentArray[maxIndex].va_max);

/*        printf("coreNum=%d, va_min = 0x%08x, minIndex=%d, va_max = 0x%08x, maxIndex=%d, minPa=0x%08x maxPa=0x%08x\n",coreNum,(unsigned)*va_min, minIndex, (unsigned)*va_max, maxIndex, minPaBase, maxPa); */
/*        fflush(stdout); */
  }
  else
  {
    *va_min = (void*)0;
    *va_max = (void*)0xFFFFFFFF;
  }
}


/*!
 * \brief Convert a 32-bit virtual address to a 32-bit physical address
 *
 * This function is a wrapper around _bgp_Virtual2Physical(), only it combines
 * its 36-bit output into a 32-bit physical address by right-shifting it 4 bits.
 * Thus, the physical address returned corresponds to the input virtual address
 * rounded down to the next lowest 16-byte boundary.
 *
 * \param[in]   VA     32-bit virtual address to be converted
 * \param[in]   vsize  Size in bytes of virtual range
 * \param[out]  pPA    Pointer to 32-bit physical address.  The output physical
 *                     address is returned in the storage pointed to by pPA.
 *
 * \retval   0  Successful.  The output physical address is in *pPA
 * \retval  -1  Invalid Virtual Address for this process.  *pPA unmodified.
 * \retval  -2  The range from VA to (VA+vsize-1) is not physically 
 *              contiguous
 * \retval  -3  Virtual Address is in Scratch, but no Scratch, or not enough 
 *              Scratch, is enabled.  *pPA unmodified.
 * 
 */
__INLINE__ int Kernel_VaTo4bitShiftedPa(void     *VA,
					size_t    vsize,
					uint32_t *pPA )
{
  int rc;
  uint32_t ua_out, pa_out;

  SPI_assert( pPA != NULL );

  rc = Kernel_Virtual2Physical(VA,
			       vsize,
			       &ua_out,
			       &pa_out );

  if ( rc == 0 )
    {
      *pPA = (ua_out << 28) | (pa_out >> 4);
    }

  return (rc);
}


/*
 *------------------------------------------------------------------------------
 *
 * The following are inline function wrappers around system calls that
 * operate on DMA counters.
 *
 *------------------------------------------------------------------------------
 */


/*!
 * \brief Query Free DMA Counter Subgroups within a Group
 *
 * This function is a wrapper around a system call that returns a list of the 
 * free (available) subgroups within the specified group.
 *
 * \param[in]   type           Specifies whether this is an injection or 
 *                             reception counter group (DMA_Type_Injection
 *                             or DMA_Type_Reception)
 * \param[in]   grp            Group number being queried (0 to 
 *                             DMA_NUM_COUNTER_GROUPS-1)
 * \param[out]  num_subgroups  Pointer to an int where the number of free
 *                             subgroups in the specified group is returned
 * \param[out]  subgroups      Pointer to an array of num_subgroups ints where
 *                             the list of num_subgroups subgroups is returned.
 *                             Each int is the subgroup number 
 *                             (0 to DMA_NUM_COUNTERS_PER_SUBGROUP-1).  The 
 *                             caller must provide space for 
 *                             DMA_NUM_COUNTERS_PER_SUBGROUP ints, in case the
 *                             entire counter group is free.
 *
 * \retval  0  Successful.  num_subgroups and subgroups array set as described.
 * \retval  -1 Unsuccessful.  errno gives the reason.
 *
 * \note The kernel may need to synchronize with other cores performing
 *       allocate or free syscalls.
 *
 */
__INLINE__ int DMA_CounterGroupQueryFree(
					 DMA_Type_t  type,
					 int         grp,
					 int        *num_subgroups,
					 int        *subgroups
					)
{
  return Kernel_CounterGroupQueryFree( (uint32_t)type, 
				       grp, 
				       (uint32_t*)num_subgroups, 
				       (uint32_t*)subgroups);
}


/*!
 * \brief Allocate DMA Counters From A Group 
 *
 * This function is a wrapper around a system call that allocates DMA counters
 * from the specified group.  Counters may be allocated in subgroups of 
 * DMA_NUM_COUNTERS_PER_SUBGROUP counters.  Parameters specify how interrupts,
 * generated when a counter hits zero, are to be handled.  A  
 * DMA_CounterGroup_t structure is returned for use in other inline 
 * functions to operate on the allocated counters.
 *
 * \param[in]   type           Specifies whether this is an injection or 
 *                             reception counter group (DMA_Type_Injection
 *                             or DMA_Type_Reception)
 * \param[in]   grp            Group number whose counters are being allocated 
 *                             (0 to DMA_NUM_COUNTER_GROUPS-1)
 * \param[in]   num_subgroups  Number of subgroups to be allocated from the group
 *                             (1 to DMA_NUM_COUNTERS_PER_SUBGROUP)
 * \param[in]   subgroups      Pointer to an array of num_subgroups ints where
 *                             the list of subgroups to be allocated is provided.
 *                             Each int is the subgroup number 
 *                             (0 to num_subgroups-1).
 * \param[in]   target         The core that will receive the interrupt when a
 *                             counter in this allocation hits zero 
 *                             (0 to DMA_NUM_COUNTER_GROUPS-1)
 * \param[in]   handler        A pointer to the function to receive control in 
 *                             the I/O thread to handle the interrupt when a 
 *                             counter in this allocation hits zero.  This 
 *                             function must be coded to take 4 uint32_t 
 *                             parameters:
 *                             - A pointer to storage specific to this 
 *                               handler.  This is the handler_parm 
 *                               specified on this allocation function.
 *                             - Three unint32_t parameters that are not used.
 *                             If handler is NULL, hit-zero interrupts will not 
 *                             be enabled for these counters.
 * \param[in]   handler_parm   A pointer to storage that should be passed to the 
 *                             interrupt handling function (see handler 
 *                             parameter)
 * \param[in]   interruptGroup A InterruptGroup_t that identifies the
 *                             group of interrupts that the counters being
 *                             allocated will become part of.
 * \param[out]  cg_ptr         Pointer to a structure that is filled in upon 
 *                             successful return for use in other inline 
 *                             functions to operate on the allocated counters.
 *                             \li counter -     Array of software counter
 *                                               structures.  Each element
 *                                               points to the corresponding
 *                                               hardware counter in DMA SRAM.
 *                                               Pointers are null if not 
 *                                               allocated).
 *                                               Counters are initialized to 
 *                                               DMA_COUNTER_INIT_VAL, 
 *                                               disabled, their hit_zero bit
 *                                               is off, base and max are NULL.
 *                             \li status_ptr  - Points to status area within the
 *                                               DMA memory map.
 *                             \li permissions - Bits set for each allocated 
 *                                               counter
 *                             \li grp_permissions - Permissions for each 
 *                                                   subgroup
 *                             \li group_id    - The group number
 *                             \li type        - The type of DMA (injection or
 *                                               reception)
 *
 * \retval  0  Successful.  Counters allocated and cg_ptr structure filled in as 
 *                          described.
 * \retval  -1 Unsuccessful.  errno gives the reason.  Nothing has been 
 *                            allocated.
 *
 * \note The kernel may need to synchronize with other cores performing queries
 *       or frees.
 *
 */
__INLINE__ int  DMA_CounterGroupAllocate(
		        DMA_Type_t                type,
			int                       grp,
			int                       num_subgroups,
			int                      *subgroups,
			int                       target,
			Kernel_CommThreadHandler  handler,
			void                     *handler_parm,
			Kernel_InterruptGroup_t   interruptGroup,
			DMA_CounterGroup_t       *cg_ptr
		       )
{
  int rc;
  /*
   * Initialize the Counter Application Segment array and its global pointer if
   * it has not been initialized yet.
   */
  if ( DMA_CounterAppSegmentArray == NULL )
  {
    rc = DMA_CounterInitAppSegments();
    if (rc) return(rc);
  }

  /*
   * If an interrupt handler has been specified, invoke the system call
   * to configure the kernel to invoke the handler when the hit zero
   * interrupt fires.
   */

  if (handler)
  {
    /*
     * Calculate the IRQ to be one of 
     * - 0: inj counter hit zero vector 0
     * - 1: inj counter hit zero vector 1
     * - 2: inj counter hit zero vector 2
     * - 3: inj counter hit zero vector 3
     *
     * - 4: rec counter hit zero vector 0
     * - 5: rec counter hit zero vector 1
     * - 6: rec counter hit zero vector 2
     * - 7: rec counter hit zero vector 3
     * based on the counter type and the DMA group number.
     */
    unsigned irqInGroup = (type == DMA_Type_Injection) ? 0 + grp : 4 + grp;

    /*
     * Calculate an interrupt ID, which is the BIC interrupt group (3)
     * combined with the IRQ number.
     */
    int interruptID = Kernel_MkInterruptID(_BGP_IC_DMA_NFT_G3_HIER_POS, 
					   irqInGroup);
    
    /*
     * Calculate the opcode indicating
     * - the target core for interrupt
     * - to call the specified function when the interrupt fires
     * - to disable interrupts before calling the specified function
     * - to enable interrupts after callling the specified function
     */
    int opcode = ( COMMTHRD_OPCODE_CORE0 + target ) | 
                   COMMTHRD_OPCODE_CALLFUNC |
                   COMMTHRD_OPCODE_DISABLEINTONENTRY |
                   COMMTHRD_OPCODE_ENABLEINTONPOOF  ;
    
    /*
     * Configure this interrupt with the kernel.
     */
    rc = Kernel_SetCommThreadConfig(interruptID, 
				    opcode, 
				    (uint32_t*)interruptGroup,
				    handler,
				    (uint32_t)handler_parm,
				    (uint32_t)NULL,
				    (uint32_t)NULL,
				    (uint32_t)NULL);
    if (rc) return rc;
  }

  /*
   * Invoke the system call to allocate the counters.
   * This system call also sets up the DMA DCRs to interrupt when the
   * counters hit zero.
   */
  rc = Kernel_CounterGroupAllocate( (uint32_t)type, 
				    grp, 
				    num_subgroups, 
				    (uint32_t*)subgroups, 
				    target,
				    (uint32_t) NULL, /* Handler.        Not used */
				    (uint32_t*)NULL, /* Handler_parm.   Not used */
				    (uint32_t) NULL, /* InterruptGroup. Not used */
				    (uint32_t*)cg_ptr);
  return rc;
}


/*!
 * \brief Free DMA Counters From A Group 
 *
 * This function is a wrapper around a system call that frees DMA counters
 * from the specified group.  Counters may be freed in subgroups of 
 * DMA_NUM_COUNTERS_PER_SUBGROUP counters.
 *
 * \param[in]   grp            Group number whose counters are being freed 
 *                             (0 to DMA_NUM_COUNTER_GROUPS-1)
 * \param[in]   num_subgroups  Number of subgroups to be freed from the group 
 *                             (1-DMA_NUM_COUNTERS_PER_SUBGROUP)
 * \param[in]   subgroups      Pointer to an array of num_subgroups ints where
 *                             the list of subgroups to be freed is provided.
 *                             Each int is the subgroup number 
 *                             (0 to DMA_NUM_COUNTERS_PER_SUBGROUP-1).
 * \param[out]  cg_ptr         Pointer to the structure previously filled in when
 *                             these counters were allocated.  Upon successful 
 *                             return, this structure is updated to reflect the 
 *                             freed counters:
 *                             \li counter[]  -  Counter structures Pointers to 
 *                                               freed counters nulled.
 *                             \li permissions - Bits cleared for each freed 
 *                                               counter.
 *
 * \retval  0  Successful.  Counters freed and cg_ptr structure updated as 
 *                          described.
 * \retval  -1 Unsuccessful.  errno gives the reason.
 *
 * \note The kernel may need to synchronize with other cores performing allocates
 *       or queries.
 */
__INLINE__ int  DMA_CounterGroupFree(
				     int                 grp,
				     int                 num_subgroups,
				     int                *subgroups,
				     DMA_CounterGroup_t *cg_ptr
				    )
{
  return Kernel_CounterGroupFree( grp, 
				  num_subgroups, 
				  (uint32_t*)subgroups, 
				  (uint32_t*)cg_ptr);
}



/*!
 * \brief Enable or Disable Counter Overflow and Underflow Interrupts
 *
 * This function is a wrapper around a system call that enables or disables
 * the 4 counter overflow/underflow interrupts for all counters:
 * 1. Injection counter overflow
 * 2. Injection counter underflow
 * 3. Reception counter overflow
 * 4. Reception counter underflow
 *
 * \param[in]  enable  Specifies whether to enable or disable the interrupts
 *                     0 = Disable, 1 = Enable.
 *
 * \retval  0            Successful
 * \retval  error_value  An error value defined in the _BGP_RAS_DMA_ErrCodes
 *                       enum located in bgp/arch/include/common/bgp_ras.h
 *
 */
__INLINE__ int DMA_CounterInterruptControl(unsigned int enable)
{
  return Kernel_ChgCounterInterruptEnables( (uint32_t)enable );

}



/*
 * -----------------------------------------------------------------------------
 * The following inline functions operate directly on the Hardware DMA Counter.
 * Note that MSYNC and MBAR are not performed by these hardware functions...
 * it is up to the caller to perform them.
 *------------------------------------------------------------------------------
 */


/*!
 * \brief Set DMA Hardware Counter Value
 *
 * Set a DMA hardware counter's value, given a pointer to the hardware counter.
 *
 * \param[in]  c_hw   Pointer to the hardware counter structure
 * \param[in]  value  The value to be set into the counter
 *
 * \return None
 *
 * \note No MSYNC or MBAR is done in this function.  It is the responsibility
 *       of the caller to do it.
 *
 */
__INLINE__ void DMA_CounterSetValueHw( 
				      DMA_CounterHw_t *c_hw,
				      unsigned int     value
				     )
{
  SPI_assert( c_hw != NULL );

  c_hw->counter = value; 
}  


/*!
 * \brief Set DMA Hardware Counter Base
 *
 * Set a DMA hardware counter's base, given a pointer to the hardware counter.
 *
 * \param[in]  c_hw     Pointer to the hardware counter structure
 * \param[in]  pa_base  The base physical address to be associated with the 
 *                      counter.  16B-aligned 4-bit shifted physical address.
 *
 * \return None
 *
 * \note No MSYNC or MBAR is done in this function.  It is the responsibility
 *       of the caller to do it.
 *
 */
__INLINE__ void DMA_CounterSetBaseHw(  
				     DMA_CounterHw_t *c_hw,
				     unsigned int     pa_base
				     )
{
  SPI_assert( c_hw != NULL );

  c_hw->pa_base = pa_base;
}


/*!
 * \brief Increment DMA Hardware Counter Value
 *
 * Increment a DMA hardware counter's value, given a pointer to the hardware 
 * counter.
 *
 * \param[in]  c_hw  Pointer to the hardware counter structure
 * \param[in]  incr  The amount to increment the counter by
 *
 * \return None
 *
 * \note No MSYNC or MBAR is done in this function.  It is the responsibility
 *       of the caller to do it.
 *
 */
__INLINE__ void DMA_CounterIncrementHw(  
				       DMA_CounterHw_t *c_hw,
				       unsigned int     incr
				      )
{
  SPI_assert( c_hw != NULL );

  c_hw->increment = incr;
}


/*!
 * \brief Decrement DMA Hardware Counter Value
 *
 * Decrement a DMA hardware counter's value, given a pointer to the hardware 
 * counter.
 *
 * \param[in]  c_hw  Pointer to the hardware counter structure
 * \param[in]  decr  The amount to decrement the counter by
 *
 * \return None
 *
 * \note No MSYNC or MBAR is done in this function.  It is the responsibility
 *       of the caller to do it.
 *
 * \note The counter overflow interrupt will fire as a result of this operation.
 *       Consider disabling this interrupt.
 *
 */
__INLINE__ void DMA_CounterDecrementHw(  
				       DMA_CounterHw_t *c_hw,
				       unsigned int     decr
				      )
{
  SPI_assert( c_hw != NULL );

  /* Decrement the counter by incrementing with a large value, which will
   * cause the counter to wrap.
   */
  c_hw->increment = (0 - decr);
}


/*!
 * \brief Set Reception DMA Hardware Counter Max
 *
 * Set a reception DMA hardware counter's maximum payload address, given a 
 * pointer to the hardware counter.
 *
 * \param[in]  c_hw    Pointer to the hardware counter structure
 * \param[in]  pa_max  The max physical address to be associated with the 
 *                     counter.  16B-aligned 4-bit shifted physical address.
 *
 * \return None
 *
 * \pre The caller has ASSERTed that (c_hw != NULL)
 *
 * \note No MSYNC or MBAR is done in this function.  It is the responsibility
 *       of the caller to do it.
 *
 */
__INLINE__ void DMA_CounterSetMaxHw(  
				    DMA_CounterHw_t *c_hw,
				    unsigned int     pa_max
				   )
{
  c_hw->pa_max = pa_max;
}


/*!
 * \brief Set DMA Hardware Counter Value and Base
 *
 * Set a DMA hardware counter's value and base, given a pointer to the hardware 
 * counter.
 *
 * \param[in]  c_hw     Pointer to the hardware counter structure
 * \param[in]  value    The value to be set into the counter
 * \param[in]  pa_base  The base physical address to be associated with the 
 *                      counter.  16B-aligned 4-bit shifted physical address.
 *
 * \return None
 *
 * \note No MSYNC or MBAR is done in this function.  It is the responsibility
 *       of the caller to do it.
 *
 */
__INLINE__ void DMA_CounterSetValueBaseHw( 
					  DMA_CounterHw_t *c_hw,
					  unsigned int     value,
					  unsigned int     pa_base
					 )
{
  SPI_assert( c_hw != NULL );

  c_hw->counter = value; 
  c_hw->pa_base = pa_base; 

}


/*!
 * \brief Set Reception DMA Hardware Counter Value, Base, and Max
 *
 * Set a reception DMA hardware counter's value, base, and max, given a pointer
 * to the hardware counter.
 *
 * \param[in]  c_hw     Pointer to the hardware counter structure
 * \param[in]  value    The value to be set into the counter
 * \param[in]  pa_base  The base physical address to be associated with the 
 *                      counter.  16B-aligned 4-bit shifted physical address.
 * \param[in]  pa_max   The max physical address to be associated with the 
 *                      counter.  16B-aligned 4-bit shifted physical address.
 *
 * \return None
 *
 * \note No MSYNC or MBAR is done in this function.  It is the responsibility
 *       of the caller to do it.
 *
 */
__INLINE__ void DMA_CounterSetValueBaseMaxHw( 
					     DMA_CounterHw_t *c_hw,
					     unsigned int     value,
					     unsigned int     pa_base,
					     unsigned int     pa_max
					    )
{
  SPI_assert( c_hw != NULL );
  SPI_assert( pa_max >= pa_base);

  c_hw->counter = value;
  c_hw->pa_base = pa_base;
  c_hw->pa_max  = pa_max;
}


/*!
 * \brief Get DMA Hardware Counter Value
 *
 * Get a DMA hardware counter's value, given a pointer to the hardware counter.
 *
 * \param[in]  c_hw   Pointer to the hardware counter structure
 *
 * \retval  value  The current value of the counter
 *
 * \note No MSYNC or MBAR is done in this function.  It is the responsibility
 *       of the caller to do it.
 *
 */
__INLINE__ unsigned int DMA_CounterGetValueHw(  
					      const DMA_CounterHw_t *c_hw
					     )
{
  SPI_assert( c_hw != NULL );

  return( c_hw->counter );
}


/*!
 * \brief Get DMA Hardware Counter Base
 *
 * Get a DMA hardware counter's base, given a pointer to the hardware counter.
 *
 * \param[in]  c_hw     Pointer to the hardware counter structure
 *
 * \retval  pa_base  The base physical address associated with the counter.
 *                   16B-aligned 4-bit shifted physical address.
 *
 * \note No MSYNC or MBAR is done in this function.  It is the responsibility
 *       of the caller to do it.
 *
 */
__INLINE__ unsigned int DMA_CounterGetBaseHw(
					     const DMA_CounterHw_t *c_hw
					    )
{
  SPI_assert( c_hw != NULL );

  return( c_hw->pa_base );
}


/*!
 * \brief Get Reception DMA Hardware Counter Max
 *
 * Get a reception DMA hardware counter's max payload address, given a pointer 
 * to the hardware counter.
 *
 * \param[in]  c_hw     Pointer to the hardware counter structure
 *
 * \retval  pa_max  The max physical address associated with the counter.
 *                  16B-aligned 4-bit shifted physical address.
 *
 * \note No MSYNC or MBAR is done in this function.  It is the responsibility
 *       of the caller to do it.
 *
 */
__INLINE__ unsigned int DMA_CounterGetMaxHw(  
					    const DMA_CounterHw_t *c_hw
					   )
{
  SPI_assert( c_hw != NULL );

  return( c_hw->pa_max );
}




/*
 * -----------------------------------------------------------------------------
 * The following inline functions operate indirectly on a hardware DMA counter
 * through the Software DMA Counter structure.
 *------------------------------------------------------------------------------
 */




/*!
 * \brief Set DMA Counter Value
 *
 * Set a DMA counter's value, given a pointer to the software DMA counter
 * structure.
 *
 * \param[in]  c_sw   Pointer to the software counter structure
 * \param[in]  value  The value to be set into the counter
 *
 * \return None
 *
 */
__INLINE__ void DMA_CounterSetValue(  
				    DMA_Counter_t *c_sw, 
				    unsigned int   value
				   )
{
  SPI_assert( c_sw != NULL );

  DMA_CounterSetValueHw(c_sw->counter_hw_ptr,
			value);
  _bgp_mbar();    /* Make sure these writes have been accepted by the memory */
                  /* system before continuing                                */  
                         
}


/*!
 * \brief Set DMA Counter Base Address
 *
 * Set a DMA counter's base address, given a pointer to the software counter
 * structure.
 *
 * \param[in]  c_sw        Pointer to the software counter structure
 * \param[in]  va_base_in  The base virtual address to be associated with the 
 *                         counter.
 *
 * \retval  0   Success
 * \retval  -1  Failure.  errno contains the reason.  Most likely EFAULT due to
 *              the va_base_in being a bad virtual address.
 *
 * \post In the software counter structure, va_base and pa_base are set.
 *       In the hardware counter structure, pa_base is set.
 *
 * \note This function does an MBAR after setting the counter to ensure the
 *       writes have been accepted by the memory system before allowing other
 *       memory accesses to to occur.
 *
 * \note The va_base in the software counter structure is the va_base_in rounded
 *       down to the next lowest 16B-aligned address.  The pa_base is the 4-bit
 *       shifted version of va_base.
 *
 */
__INLINE__ int DMA_CounterSetBase(  
				  DMA_Counter_t *c_sw, 
				  void          *va_base_in
				 )
{
  int rc;

  SPI_assert( c_sw != NULL );

  /*
   * 16-B align the virtual address and store result in software counter
   * structure
   */
  c_sw->va_base = (char*)( (unsigned)va_base_in & 0xFFFFFFF0 );

  rc =  Kernel_VaTo4bitShiftedPa(c_sw->va_base, 
				 1,
				 &(c_sw->pa_base) ); 
  if ( rc != 0 ) 
    {
      errno = EFAULT;
      return (-1);
    }

  /*
   * Write physical address to the hardware counter 
   */
  DMA_CounterSetBaseHw(c_sw->counter_hw_ptr,
		       c_sw->pa_base); 

  _bgp_mbar();    /* Make sure these writes have been accepted by the memory */
                  /* system before continuing                                */

  return (0);

}


/*!
 * \brief Increment DMA Counter
 *
 * Increment a DMA counter's value, given a pointer to the software counter
 * structure.
 *
 * \param[in]  c_sw   Pointer to the software counter structure
 * \param[in]  incr   The amount to increment the counter by
 *
 * \return None
 *
 * \note This function does an MBAR after setting the counter to ensure the
 *       writes have been accepted by the memory system before allowing other
 *       memory accesses to to occur.
 *
 */
__INLINE__ void DMA_CounterIncrement(  
				     DMA_Counter_t *c_sw, 
				     unsigned int   incr
				    )
{
  SPI_assert( c_sw != NULL );

  DMA_CounterIncrementHw(c_sw->counter_hw_ptr,
			 incr); 

  _bgp_mbar();    /* Make sure these writes have been accepted by the memory */
                  /* system before continuing                                */
}


/*!
 * \brief Decrement DMA Counter
 *
 * Decrement a DMA counter's value, given a pointer to the software counter
 * structure.
 *
 * \param[in]  c_sw   Pointer to the software counter structure
 * \param[in]  decr   The amount to decrement the counter by
 *
 * \return None
 *
 * \note This function does an MBAR after setting the counter to ensure the
 *       writes have been accepted by the memory system before allowing other
 *       memory accesses to to occur.
 *
 */
__INLINE__ void DMA_CounterDecrement(  
				     DMA_Counter_t *c_sw, 
				     unsigned int   decr
				    )
{
  SPI_assert( c_sw != NULL );

  DMA_CounterDecrementHw(c_sw->counter_hw_ptr,
			 decr); 

  _bgp_mbar();    /* Make sure these writes have been accepted by the memory */
                  /* system before continuing                                */
}


/*!
 * \brief Set DMA Counter Max Address
 *
 * Set a DMA counter's max address, given a pointer to the software counter
 * structure.
 *
 * \param[in]  c_sw       Pointer to the software counter structure
 * \param[in]  va_max_in  The max virtual address to be associated with the 
 *                        counter.
 *
 * \retval  0   Success
 * \retval  -1  Failure.  errno contains the reason.  Most likely EFAULT due to
 *              the va_max_in being a bad virtual address.
 *
 * \post In the software counter structure, va_max and pa_max are set.
 *       In the hardware counter structure, pa_max is set.
 *
 * \note This function does an MBAR after setting the counter to ensure the
 *       writes have been accepted by the memory system before allowing other
 *       memory accesses to to occur.
 *
 * \note The va_max in the software counter structure is the va_max_in rounded 
 *       up to the next highest 16B-aligned address.  The pa_max is the 4-bit
 *       shifted version of va_max.
 *
 */
__INLINE__ int DMA_CounterSetMax(  
				 DMA_Counter_t *c_sw, 
				 void          *va_max_in
				)
{
  int rc;

  SPI_assert( c_sw != NULL );

  /*
   * Round up to 16B boundary and 16-B align the virtual address and store 
   * result in software counter structure.
   */
  c_sw->va_max = (char*) ( (unsigned)va_max_in & 0xFFFFFFF0 );
  if ( c_sw->va_max != va_max_in )  c_sw->va_max = (char*)c_sw->va_max + 0x00000010;

  /*
   * Get the 16B-aligned 4-bit shifted physical address from the virtual address.
   */
  rc = Kernel_VaTo4bitShiftedPa(c_sw->va_max, 
				1,
				&(c_sw->pa_max) ); 

  if ( rc != 0 ) 
    {
      errno = EFAULT;
      return (-1);
    }

  /*
   * Write physical address to the hardware counter 
   */
  DMA_CounterSetMaxHw(c_sw->counter_hw_ptr,
		      c_sw->pa_max); 

  _bgp_mbar();    /* Make sure these writes have been accepted by the memory */
                  /* system before continuing                                */      

  return (0);

}


/*!
 * \brief Set DMA Counter Value and Base Address
 *
 * Set a DMA counter's value and base address, given a pointer to the software
 * counter structure.
 *
 * \param[in]  c_sw        Pointer to the software counter structure
 * \param[in]  value       The value to be set into the counter
 * \param[in]  va_base_in  The base virtual address to be associated with the 
 *                         counter.
 *
 * \retval  0   Success
 * \retval  -1  Failure.  errno contains the reason.  Most likely EFAULT due to
 *              the va_base_in being a bad virtual address.
 *
 * \post In the software counter structure, va_base and pa_base are set.
 *       In the hardware counter structure, pa_base and value are set.
 *
 * \note This function does an MBAR after setting the counter to ensure the
 *       writes have been accepted by the memory system before allowing other
 *       memory accesses to to occur.
 *
 * \note The va_base in the software counter structure is the va_base_in rounded
 *       down to the next lowest 16B-aligned address.  The pa_base is the 4-bit
 *       shifted version of va_base.
 *                      
 */
__INLINE__ int DMA_CounterSetValueBase( 
				       DMA_Counter_t *c_sw,
				       unsigned int   value,
				       void          *va_base_in
				      )  
{
  int rc=0;

  SPI_assert( c_sw != NULL );

  /*
   * 16-B align the virtual address and store result in software counter
   * structure
   */
  c_sw->va_base = (char*) ( (unsigned)va_base_in & 0xFFFFFFF0 ); 

  /*
   * Get the 16B-aligned 4-bit shifted physical address from the virtual address.
   */
  rc = Kernel_VaTo4bitShiftedPa(c_sw->va_base, 
				1,
				&(c_sw->pa_base) ); 
  if ( rc != 0 ) 
    {
      errno = EFAULT;
      return (-1);
    }

  /*
   * Write the value and physical address to the hardware counter 
   */
  DMA_CounterSetValueBaseHw(c_sw->counter_hw_ptr,
			    value,
			    c_sw->pa_base );

  _bgp_mbar();    /* Make sure these writes have been accepted by the memory */
                  /* system before continuing                                */  

  return (0);
}


/*!
 * \brief Set DMA Counter Value, Base Address, and Max Address
 *
 * Set a reception DMA counter's value, base address, and max address, given a 
 * pointer to the software counter structure.
 *
 * \param[in]  c_sw        Pointer to the software counter structure
 * \param[in]  value       The value to be set into the counter
 * \param[in]  va_base_in  The base virtual address to be associated with the 
 *                         counter.
 * \param[in]  va_max_in   The max virtual address to be associated with the 
 *                         counter.
 *
 * \retval  0   Success
 * \retval  -1  Failure.  errno contains the reason.  Most likely EFAULT due to
 *              the va_base_in or va_max_in being a bad virtual address.
 *
 * \post In the software counter structure, va_base, pa_base, va_max, and pa_max
 *       are set.  In the hardware counter structure, pa_base, pa_max, and value 
 *       are set.
 *
 * \note This function does an MBAR after setting the counter to ensure the
 *       writes have been accepted by the memory system before allowing other
 *       memory accesses to to occur.
 *
 * \note The va_base in the software counter structure is the va_base_in rounded
 *       down to the next lowest 16B-aligned address.  The pa_base is the 4-bit
 *       shifted version of va_base.
 *
 * \note The va_max in the software counter structure is the va_max_in rounded 
 *       up to the next highest 16B-aligned address.  The pa_max is the 4-bit
 *       shifted version of va_max.
 *                      
 */
__INLINE__ int DMA_CounterSetValueBaseMax( 
					  DMA_Counter_t *c_sw,
					  unsigned int   value,
					  void          *va_base_in,
					  void          *va_max_in
					 )  
{
  int rc=0;
  void *va_base, *va_max;
  unsigned int pa_base, pa_max;

  SPI_assert( c_sw != NULL );

  /* 
   * Process the base address: 
   * - 16-B align the virtual address and store result in software counter
   *   structure
   * - Get the 16B-aligned 4-bit shifted physical address from the virtual 
   *   address.
   */
  va_base = c_sw->va_base = (char*) ( (unsigned)va_base_in & 0xFFFFFFF0 ); 

  rc = Kernel_VaTo4bitShiftedPa(va_base, 
				1,
				&pa_base ); 
  if ( rc != 0 ) 
    {
      errno = EFAULT;
      return (-1);
    }

  c_sw->pa_base = pa_base;

  /*
   * Process the max address:
   * - 16B align the virtual address and store result in software counter structure.
   *   Note: we can't round up or the address may be one byte out of range.
   * - Get the 16B-aligned 4-bit shifted physical address from the virtual 
   *   address.
   */
  va_max = (char*) ( (unsigned)va_max_in & 0xFFFFFFF0 );

  rc = Kernel_VaTo4bitShiftedPa(va_max, 
				1,
				&pa_max ); 
/*    printf("SetValueBaseMax: va_max_in=0x%08x, va_max=0x%08x, pa_max=0x%08x, rc=%d\n",(unsigned)va_max_in, (unsigned)va_max,pa_max,rc); */
/*    fflush(stdout); */
  if ( rc != 0 ) 
    {
      errno = EFAULT;
      return (-1);
    }

  c_sw->pa_max = pa_max;

  /*
   * Write the value, base, and max to the hardware counter 
   */
  DMA_CounterSetValueBaseMaxHw(c_sw->counter_hw_ptr,
			       value,
			       pa_base,
			       pa_max); 

  _bgp_mbar();    /* Make sure these writes have been accepted by the memory */
                  /* system before continuing                                */  

  return (0);
}


/*!
 * \brief Get DMA Counter Value
 *
 * Get a DMA counter's value, given a pointer to the software counter
 * structure.
 *
 * \param[in]  c_sw  Pointer to the software counter structure
 *
 * \retval  value  The value of the specified counter
 *
 * \note This function does an MSYNC after fetching the counter's value
 *       to ensure that the data that was just DMA'd is available to all
 *       cores.
 *
 */
__INLINE__ unsigned int DMA_CounterGetValue(  
					    const DMA_Counter_t *c_sw
					   )
{
  unsigned int val;

  SPI_assert( c_sw != NULL );

  val = DMA_CounterGetValueHw( c_sw->counter_hw_ptr );

  _bgp_msync();

  return val;  

}


/*!
 * \brief Get DMA Counter Value with No Msync
 *
 * Get a DMA counter's value, given a pointer to the software counter
 * structure.  No Msync is done.  It is up to the caller to do it,
 * if necessary.
 *
 * \param[in]  c_sw  Pointer to the software counter structure
 *
 * \retval  value  The value of the specified counter
 *
 */
__INLINE__ unsigned int DMA_CounterGetValueNoMsync(  
						   const DMA_Counter_t *c_sw
						  )
{
  unsigned int val;

  SPI_assert( c_sw != NULL );

  val = DMA_CounterGetValueHw( c_sw->counter_hw_ptr );

  return val;  

}


/*!
 * \brief Get DMA Base Address
 *
 * Get a DMA counter's base virtual address, given a pointer to the software
 * counter structure.
 *
 * \param[in]  c_sw  Pointer to the software counter structure
 *
 * \retval  va_base  The base virtual address associated with the specified 
 *                   counter
 *
 * \note This returns the shadow va_base directly out of the software counter
 *       structure.  This should correspond with the physical address in the
 *       hardware counter, but it is a rounded-down-to-the-previous-16B-boundary
 *       version of the actual base virtual address of the buffer the caller is
 *       working with.
 *
 */
__INLINE__ void * DMA_CounterGetBase(  
				     const DMA_Counter_t *c_sw
				    )
{
  SPI_assert( c_sw != NULL );

  return( c_sw->va_base );
}


/*!
 * \brief Get Reception DMA Max Address
 *
 * Get a reception DMA counter's max virtual address, given a pointer to 
 * the software counter structure.
 *
 * \param[in]  c_sw  Pointer to the software counter structure
 *
 * \retval  va_max  The max virtual address associated with the specified 
 *                  counter
 *
 * \note This returns the shadow va_max directly out of the software counter
 *       structure.  This should correspond with the physical address in the
 *       hardware counter, but it is a rounded-up-to-the-next-16B-boundary
 *       version of the actual max virtual address of the buffer the caller is
 *       working with.
 *
 */
__INLINE__ void *DMA_CounterGetMax(  
				   const DMA_Counter_t *c_sw
				  )
{
  SPI_assert( c_sw != NULL );

  return( c_sw->va_max );
}


/*!
 * \brief Get Offset from DMA Base Address
 *
 * Given a virtual address, get the offset from the base address associated with
 * a counter.
 *
 * \param[in]  c_sw  Pointer to the software counter structure
 * \param[in]  va    Virtual address whose offset from the counter's base is
 *                   to be returned.
 * \param[in]  length   The number of bytes in the buffer pointed to by va.
 * \param[in]  coreNum  The number of the core in which the virtual
 *                      address resides (0 to DMA_MAX_NUM_CORES).
 *
 * \retval  offset   The offset of the va from the counter's base.
 *
 * \note This uses the counter's physical base address and the application's 
 *       memory segments (see DMA_CounterAppSegment_t).
 *
 * \note It is assumed that if the coreNum is not our core, then the counter's
 *       base address (used in calculating the offset) is the smallest physical 
 *       address accessible from user space on coreNum
 *       (DMA_CounterMinPaAccessibleFromUserMode[coreNum]).
 *
 */
__INLINE__ unsigned int DMA_CounterGetOffsetFromBase(  
						     const DMA_Counter_t *c_sw,
						     void                *va,
						     unsigned int         length,
						     unsigned int         coreNum
						    )
{
  SPI_assert( c_sw != NULL );
  SPI_assert( va   != NULL );
  SPI_assert ( coreNum >= 0 );
  SPI_assert ( coreNum <= DMA_MAX_NUM_CORES );

  DMA_CounterAppSegment_t *appSegmentArray = DMA_CounterGetAppSegments( coreNum );
  uint32_t                 numAppSegments;
  uint32_t                 i;
  uint32_t                 segmentVaBase;
  uint32_t                 offset;
  uint32_t                 ourCoreNum = Kernel_PhysicalProcessorID();
  uint32_t                 counterPaBase;


  /* Determine which application segment the virtual address is in. */
  /* First, check the last app segment accessed.                    */
  i             = DMA_CounterCachedAppSegmentIndex[coreNum];
  segmentVaBase = appSegmentArray[i].va_base;
  if ( ! ( ((uint32_t)va >= segmentVaBase) &&
	   ((uint32_t)va - segmentVaBase < appSegmentArray[i].length) ) )
  {
    /* It is not the last app segment accessed.  Search them. */
    numAppSegments = DMA_CounterGetNumAppSegments();
    for (i=0; i<numAppSegments; i++)
    {
      segmentVaBase = appSegmentArray[i].va_base;
      if ( ((uint32_t)va >= segmentVaBase) &&
	   ((uint32_t)va - segmentVaBase < appSegmentArray[i].length) )
	break;
    }
    SPI_assert(i < numAppSegments);
    DMA_CounterCachedAppSegmentIndex[coreNum] = i;
  }

  /*
   * Make sure buffer fits in app segment.
   */
  if ( ( (uint32_t)va + length - 1 ) > appSegmentArray[i].va_max )
  {
    printf("DMA_CounterGetOffsetFromBase: Buffer 0x%08x of length %d is out of bounds.  Check length.\n",
	   (unsigned)va, length);
    SPI_assert(0);
  }
  
  /* 
   * If coreNum is our core, use the offset from our core's counter base to
   * calculate the DMA offset.
   * Otherwise, assume the counter base is the smallest physical address 
   * accessible from user space on coreNum and use that.
   */
  if ( ourCoreNum == coreNum )
    counterPaBase = c_sw->pa_base;
  else
    counterPaBase = DMA_CounterMinPaAccessibleFromUserMode[coreNum];

  /*
   * If the base physical address of the application segment found above is
   * greater than or equal to the counter's base physical address (typical
   * case), proceed with the calculation based on that.
   *
   * Otherwise, use a slightly different calculation (see else leg).
   */
  if ( appSegmentArray[i].pa_base >= counterPaBase )
  {
    /* 
     * Calculate the offset from the counter base:
     * - offset from app segment's virtual address base (va - segmentVaBase) + 
     * - segment's physical base (shifted) - counter's base (shifted) * 16
     */
    offset = 
      ((unsigned)va - segmentVaBase) +
      ( (appSegmentArray[i].pa_base - counterPaBase) << 4 );
    
/*     printf("GetOffsetFromBase:  va=0x%08x, length=%d, offset=0x%08x, index=%d, segmentVaBase=0x%08x, appSegmentArrayPaBase=0x%08x, counterBase=0x%08x\n",(unsigned)va, length, offset, i, */
/* 	   segmentVaBase, appSegmentArray[i].pa_base, counterPaBase); */
/*     fflush(stdout); */
  }
  /*
   * Handle the case where the counter's base exceeds the app segment's base.
   * This occurs when the counter's base is set to the base of the buffer
   * rather than the min base of all the app segments.
   */
  else
  {
    offset = 
      ((unsigned)va - segmentVaBase) -
      ( (counterPaBase - appSegmentArray[i].pa_base) << 4 );
    
/*     printf("GetOffsetFromBase2:  va=0x%08x, length=%d, offset=0x%08x, index=%d, segmentVaBase=0x%08x, appSegmentArrayPaBase=0x%08x, counterBase=0x%08x\n",(unsigned)va, length, offset, i, */
/* 	   segmentVaBase, appSegmentArray[i].pa_base, counterPaBase); */
/*     fflush(stdout); */
  }

  return ( offset );
}




/* 
 * ------------------------------------------------------------------------------
 *
 * The following functions access counters by specifying the group pointer and 
 * counter_id.
 *
 * ------------------------------------------------------------------------------
 */




/*!
 * \brief Set DMA Counter Value using a Counter ID
 *
 * Set a DMA counter's value, given a counter group structure and counter ID.
 *
 * \param[in]  cg_ptr      Pointer to the structure previously filled in when 
 *                         the counter was allocated
 * \param[in]  counter_id  Identifier of the counter being set 
 *                         (0 to DMA_NUM_COUNTERS_PER_GROUP-1)
 * \param[in]  value       The value to be set into the counter
 *
 * \return None
 *
 * \note This function does an MBAR after setting the counter to ensure the
 *       writes have been accepted by the memory system before allowing other
 *       memory accesses to to occur.
 *
 */
__INLINE__ void DMA_CounterSetValueById(
					DMA_CounterGroup_t *cg_ptr,
					int                 counter_id,      
					unsigned int        value
				       )
{
  SPI_assert( (counter_id >= 0) && (counter_id < DMA_NUM_COUNTERS_PER_GROUP) );
  SPI_assert( cg_ptr != NULL );
  SPI_assert( (cg_ptr->permissions[DMA_COUNTER_GROUP_WORD_ID(counter_id)] & 
	   _BN(DMA_COUNTER_GROUP_WORD_BIT_ID(counter_id))) != 0 );

  DMA_CounterSetValue( &cg_ptr->counter[counter_id],
		       value );

  // Note: it is assumed that the above function call performs an MBAR
}


/*!
 * \brief Set DMA Counter Base Address using a Counter ID
 *
 * Set a DMA counter's base address, given a counter group structure and 
 * counter ID.
 *
 * \param[in]  cg_ptr      Pointer to the structure previously filled in when
 *                         the counter was allocated.
 * \param[in]  counter_id  Identifier of the counter being set 
 *                         (0 to DMA_NUM_COUNTERS_PER_GROUP-1)
 * \param[in]  va_base_in  The base virtual address to be associated with the 
 *                         counter.
 *
 * \retval  0   Success
 * \retval  -1  Failure.  errno contains the reason.  Most likely EFAULT due to
 *              the va_base_in being a bad virtual address.
 *
 * \note This function does an MBAR after setting the counter to ensure the
 *       writes have been accepted by the memory system before allowing other
 *       memory accesses to to occur.
 *
 */
__INLINE__ int DMA_CounterSetBaseById(
				      DMA_CounterGroup_t *cg_ptr,
				      int                 counter_id,      
				      void               *va_base_in
				     )
{
  int rc;
  SPI_assert( (counter_id >= 0) && (counter_id < DMA_NUM_COUNTERS_PER_GROUP) );
  SPI_assert( cg_ptr != NULL );
  SPI_assert( (cg_ptr->permissions[DMA_COUNTER_GROUP_WORD_ID(counter_id)] & 
	   _BN(DMA_COUNTER_GROUP_WORD_BIT_ID(counter_id))) != 0 );

  rc = DMA_CounterSetBase( &cg_ptr->counter[counter_id],
			   va_base_in );

  // Note: it is assumed that the above function call performs an MBAR

  return rc;
}


/*!
 * \brief Increment DMA Counter using a Counter ID
 *
 * Increment a DMA counter's value, given a counter group structure and 
 * counter ID.
 *
 * \param[in]  cg_ptr      Pointer to the structure previously filled in when
 *                         the counter was allocated.
 * \param[in]  counter_id  Identifier of the counter being incremented
 *                         (0 to DMA_NUM_COUNTERS_PER_GROUP-1)
 * \param[in]  incr        The amount to increment the counter by
 *
 * \return None
 *
 * \note This function does an MBAR after setting the counter to ensure the
 *       writes have been accepted by the memory system before allowing other
 *       memory accesses to to occur.
 *
 */
__INLINE__ void DMA_CounterIncrementById(
					 DMA_CounterGroup_t *cg_ptr,
					 int                 counter_id,      
					 unsigned int        incr
					)
{
  SPI_assert( (counter_id >= 0) && (counter_id < DMA_NUM_COUNTERS_PER_GROUP) );
  SPI_assert( cg_ptr != NULL );
  SPI_assert( (cg_ptr->permissions[DMA_COUNTER_GROUP_WORD_ID(counter_id)] & 
	   _BN(DMA_COUNTER_GROUP_WORD_BIT_ID(counter_id))) != 0 );

  DMA_CounterIncrement( &cg_ptr->counter[counter_id],
			incr );

  // Note: it is assumed that the above function call performs an MBAR
}


/*!
 * \brief Decrement DMA Counter using a Counter ID
 *
 * Decrement a DMA counter's value, given a counter group structure and 
 * counter ID.
 *
 * \param[in]  cg_ptr      Pointer to the structure previously filled in when
 *                         the counter was allocated.
 * \param[in]  counter_id  Identifier of the counter being decremented
 *                         (0 to DMA_NUM_COUNTERS_PER_GROUP-1)
 * \param[in]  decr        The amount to decrement the counter by
 *
 * \return None
 *
 * \note This function does an MBAR after setting the counter to ensure the
 *       writes have been accepted by the memory system before allowing other
 *       memory accesses to to occur.
 *
 */
__INLINE__ void DMA_CounterDecrementById(
					 DMA_CounterGroup_t *cg_ptr,
					 int                 counter_id,      
					 unsigned int        decr
					)
{
  SPI_assert( (counter_id >= 0) && (counter_id < DMA_NUM_COUNTERS_PER_GROUP) );
  SPI_assert( cg_ptr != NULL );
  SPI_assert( (cg_ptr->permissions[DMA_COUNTER_GROUP_WORD_ID(counter_id)] & 
	   _BN(DMA_COUNTER_GROUP_WORD_BIT_ID(counter_id))) != 0 );

  DMA_CounterDecrement( &cg_ptr->counter[counter_id],
			decr );

  // Note: it is assumed that the above function call performs an MBAR
}


/*!
 * \brief Set Reception DMA Counter Max Address using a Counter ID
 *
 * Set a reception DMA counter's base address, given a counter group structure 
 * and counter ID.
 *
 * \param[in]  cg_ptr      Pointer to the structure previously filled in when
 *                         the counter was allocated.
 * \param[in]  counter_id  Identifier of the counter being set 
 *                         (0 to DMA_NUM_COUNTERS_PER_GROUP-1)
 * \param[in]  va_max_in   The max virtual address to be associated with the 
 *                         counter.
 *
 * \retval  0   Success
 * \retval  -1  Failure.  errno contains the reason.  Most likely EFAULT due to
 *              the va_max_in being a bad virtual address.
 *
 * \note This function does an MBAR after setting the counter to ensure the
 *       writes have been accepted by the memory system before allowing other
 *       memory accesses to to occur.
 *
 */
__INLINE__ int DMA_CounterSetMaxById(
				     DMA_CounterGroup_t *cg_ptr,
				     int                 counter_id,      
				     void               *va_max_in
				    ) 
{
  int rc;

  SPI_assert( (counter_id >= 0) && (counter_id < DMA_NUM_COUNTERS_PER_GROUP) );
  SPI_assert( cg_ptr != NULL );
  SPI_assert( (cg_ptr->permissions[DMA_COUNTER_GROUP_WORD_ID(counter_id)] & 
	   _BN(DMA_COUNTER_GROUP_WORD_BIT_ID(counter_id))) != 0 );

  rc = DMA_CounterSetMax( &cg_ptr->counter[counter_id],
			  va_max_in );

  // Note: it is assumed that the above function call performs an MBAR

  return rc;

}


/*!
 * \brief Get DMA Counter Value using a Counter ID
 *
 * Get a DMA counter's value, given a counter group structure and counter ID.
 *
 * \param[in]  cg_ptr      Pointer to the structure previously filled in when
 *                         the counter was allocated.
 * \param[in]  counter_id  Identifier of the counter 
 *                         (0 to DMA_NUM_COUNTERS_PER_GROUP-1)
 *
 * \retval  value  The value of the counter
 *
 */
__INLINE__ unsigned int DMA_CounterGetValueById(
					const DMA_CounterGroup_t *cg_ptr,
			                const int                 counter_id
				       )
{
  SPI_assert( (counter_id >= 0) && (counter_id < DMA_NUM_COUNTERS_PER_GROUP) );
  SPI_assert( cg_ptr != NULL );
  SPI_assert( (cg_ptr->permissions[DMA_COUNTER_GROUP_WORD_ID(counter_id)] & 
	   _BN(DMA_COUNTER_GROUP_WORD_BIT_ID(counter_id))) != 0 );

  return ( DMA_CounterGetValue( &cg_ptr->counter[counter_id] ) );
}


/*!
 * \brief Get DMA Counter Base Address using a Counter ID
 *
 * Get a DMA counter's base virtual address, given a counter group structure and 
 * counter ID.
 *
 * \param[in]  cg_ptr      Pointer to the structure previously filled in when the
 *                         counter was allocated.
 * \param[in]  counter_id  Identifier of the counter 
 *                         (0 to DMA_NUM_COUNTERS_PER_GROUP-1)
 *
 * \retval  va_base  The base virtual address associated with the specified 
 *                   counter
 *
 * \note This returns the shadow va_base directly out of the software counter
 *       structure.  This should correspond with the physical address in the
 *       hardware counter, but it is a rounded-down-to-the-previous-16B-boundary
 *       version of the actual base virtual address of the buffer the caller is
 *       working with.
 *
 */
__INLINE__ void * DMA_CounterGetBaseById(
				const DMA_CounterGroup_t *cg_ptr,
				int                       counter_id
			       )
{
  SPI_assert( (counter_id >= 0) && (counter_id < DMA_NUM_COUNTERS_PER_GROUP) );
  SPI_assert( cg_ptr != NULL );
  SPI_assert( (cg_ptr->permissions[DMA_COUNTER_GROUP_WORD_ID(counter_id)] & 
	   _BN(DMA_COUNTER_GROUP_WORD_BIT_ID(counter_id))) != 0 );

  return( DMA_CounterGetBase( &cg_ptr->counter[counter_id] ) );
}


/*!
 * \brief Get Offset from DMA Base Address using a Counter ID
 *
 * Given a virtual address, get the offset from the base address associated with
 * the specified counter.
 *
 * \param[in]  cg_ptr       Pointer to the structure previously filled in when
 *                          the counter was allocated
 * \param[in]  counter_id   Identifier of the counter
 *                          (0 to DMA_NUM_COUNTERS_PER_GROUP-1)
 * \param[in]  va           Virtual address whose offset from the counter's base
 *                          is to be returned.
 * \param[in]  length       The number of bytes in the buffer pointed to by va.
 * \param[in]  coreNum      The number of the core in which the virtual
 *                          address resides (0 to DMA_MAX_NUM_CORES).
 *
 * \retval  offset   The offset of the va from the counter's base.
 *
 * \note This works with the shadow va_base directly out of the software counter
 *       structure.  This should correspond with the physical address in the
 *       hardware counter, but it is a rounded-down-to-the-previous-16B-boundary
 *       version of the actual base virtual address of the buffer the caller is
 *       working with.
 *
 * \note No effort is given to flag the case where va is less than the base
 *       address.  In that case, (va - va_base) is returned, whatever that is.
 *
 */
__INLINE__ unsigned int DMA_CounterGetOffsetFromBaseById(  
				const DMA_CounterGroup_t *cg_ptr,
				int                       counter_id,
   			        void                     *va,
				unsigned int              length,
				unsigned int              coreNum
			       )
{
  SPI_assert( (counter_id >= 0) && (counter_id < DMA_NUM_COUNTERS_PER_GROUP) );
  SPI_assert( cg_ptr != NULL );
  SPI_assert( (cg_ptr->permissions[DMA_COUNTER_GROUP_WORD_ID(counter_id)] & 
	   _BN(DMA_COUNTER_GROUP_WORD_BIT_ID(counter_id))) != 0 );
/*   printf("Getting offset from counter %d for core %d\n",counter_id,coreNum); */
  return( DMA_CounterGetOffsetFromBase( &cg_ptr->counter[counter_id],
					va,
					length,
					coreNum ) );
}


/*!
 * \brief Get Reception DMA Counter Max Address Using a Counter ID
 *
 * Get a reception DMA counter's maximum virtual address, given a counter group 
 * structure and counter ID.
 *
 * \param[in]  cg_ptr      Pointer to the structure previously filled in when the
 *                         counter was allocated.
 * \param[in]  counter_id  Identifier of the counter 
 *                         (0 to DMA_NUM_COUNTERS_PER_GROUP-1)
 *
 * \retval  va_max  The virtual address of the max of the counter
 *
 * \note This returns the shadow va_max directly out of the software counter
 *       structure.  This should correspond with the physical address in the
 *       hardware counter, but it is a rounded-up-to-the-next-16B-boundary
 *       version of the actual max virtual address of the buffer the caller is
 *       working with.
 *
 */
__INLINE__ void * DMA_CounterGetMaxById(
				const DMA_CounterGroup_t *cg_ptr,
				const int                 counter_id
			       )
{
  SPI_assert( (counter_id >= 0) && (counter_id < DMA_NUM_COUNTERS_PER_GROUP) );
  SPI_assert( cg_ptr != NULL );
  SPI_assert( (cg_ptr->permissions[DMA_COUNTER_GROUP_WORD_ID(counter_id)] & 
	   _BN(DMA_COUNTER_GROUP_WORD_BIT_ID(counter_id))) != 0 );

  return ( DMA_CounterGetMax( &cg_ptr->counter[counter_id] ) );
}


/*!
 * \brief Set DMA Counter Value and Base Address using a Counter ID
 *
 * Set a DMA counter's value and base address, given a counter group structure
 * and counter ID.
 *
 * \param[in]  cg_ptr      Pointer to the structure previously filled in when the
 *                         counter was allocated.
 * \param[in]  counter_id  Identifier of the counter being set 
 *                         (0 to DMA_NUM_COUNTERS_PER_GROUP-1)
 * \param[in]  value       The value to be set into the counter
 * \param[in]  va_base_in  The base virtual address to be associated with the 
 *                         counter.
 *
 * \retval  0   Success
 * \retval  -1  Failure.  errno contains the reason.  Most likely EFAULT due to
 *              the va_base_in being a bad virtual address.
 *
 * \note This function does an MBAR after setting the counter to ensure the
 *       writes have been accepted by the memory system before allowing other
 *       memory accesses to to occur.
 *
 */
__INLINE__ int DMA_CounterSetValueBaseById( 
					   DMA_CounterGroup_t *cg_ptr,
					   int                 counter_id,
					   unsigned int        value,
					   void               *va_base_in
					  )
{
  int rc;
  SPI_assert( (counter_id >= 0) && (counter_id < DMA_NUM_COUNTERS_PER_GROUP) );
  SPI_assert( cg_ptr != NULL );
  SPI_assert( (cg_ptr->permissions[DMA_COUNTER_GROUP_WORD_ID(counter_id)] & 
	   _BN(DMA_COUNTER_GROUP_WORD_BIT_ID(counter_id))) != 0 );

  rc =  DMA_CounterSetValueBase( &cg_ptr->counter[counter_id],
				 value,
				 va_base_in );

  // Note: it is assumed that the above function call performs an MBAR

  return rc;
}




/*!
 * \brief Set Reception DMA Counter Value, Base Address, and Max Address  using
 *        a Counter ID
 *
 * Set a reception DMA counter's value, base address, and max address, given a 
 * counter group structure and counter ID.
 *
 * \param[in]  cg_ptr      Pointer to the structure previously filled in when the
 *                         counter was allocated.
 * \param[in]  counter_id  Identifier of the counter being set 
 *                         (0 to DMA_NUM_COUNTERS_PER_GROUP-1)
 * \param[in]  value       The value to be set into the counter
 * \param[in]  va_base_in  The base virtual address to be associated with the 
 *                         counter.
 * \param[in]  va_max_in   The max virtual address to be associated with the 
 *                         counter.
 *
 * \retval  0   Success
 * \retval  -1  Failure.  errno contains the reason.  Most likely EFAULT due to
 *              the va_base_in or va_max_in being a bad virtual address.
 *
 * \note This function does an MBAR after setting the counter to ensure the
 *       writes have been accepted by the memory system before allowing other
 *       memory accesses to to occur.
 *
 */
__INLINE__ int DMA_CounterSetValueBaseMaxById( 
				DMA_CounterGroup_t *cg_ptr,
				int                 counter_id,
				unsigned int        value,
				void               *va_base_in,
				void               *va_max_in
			       )
{
  int rc;
  SPI_assert( (counter_id >= 0) && (counter_id < DMA_NUM_COUNTERS_PER_GROUP) );
  SPI_assert( cg_ptr != NULL );
  SPI_assert( (cg_ptr->permissions[DMA_COUNTER_GROUP_WORD_ID(counter_id)] & 
	   _BN(DMA_COUNTER_GROUP_WORD_BIT_ID(counter_id))) != 0 );

  rc =  DMA_CounterSetValueBaseMax( &cg_ptr->counter[counter_id],
				    value,
				    va_base_in,
				    va_max_in );

  // Note: it is assumed that the above function call performs an MBAR

  return rc;

}


/*!
 * \brief Enable DMA Counter using a Counter ID
 *
 * Enable a DMA counter, given a counter group structure and counter ID.
 *
 * \param[in]  cg_ptr      Pointer to the structure previously filled in when the
 *                         counter was allocated.
 * \param[in]  counter_id  Identifier of the counter being enabled 
 *                         (0 to DMA_NUM_COUNTERS_PER_GROUP-1)
 *
 * \return None
 *
 * \note This function does an MBAR after setting the counter to ensure the
 *       writes have been accepted by the memory system before allowing other
 *       memory accesses to to occur.
 *
 */
__INLINE__ void DMA_CounterSetEnableById(  
				DMA_CounterGroup_t *cg_ptr, 
				int                 counter_id
			       )
{
  SPI_assert( (counter_id >= 0) && (counter_id < DMA_NUM_COUNTERS_PER_GROUP) );
  SPI_assert( cg_ptr != NULL );
  SPI_assert( (cg_ptr->permissions[DMA_COUNTER_GROUP_WORD_ID(counter_id)] & 
	   _BN(DMA_COUNTER_GROUP_WORD_BIT_ID(counter_id))) != 0 );
  SPI_assert( cg_ptr->status_ptr != 0);

  /* Enable the counter by writing 1 to the appropriate bit */
  int r =  DMA_COUNTER_GROUP_WORD_ID(counter_id);      
  int c =  DMA_COUNTER_GROUP_WORD_BIT_ID(counter_id);   
  cg_ptr->status_ptr->enable[r] = _BN(c);

  _bgp_mbar();    /* Make sure these writes have been accepted by the memory */
                  /* system before continuing                                */
}


/*!
 * \brief Disable DMA Counter using a Counter ID
 *
 * Disable a DMA counter, given a counter group structure and counter ID.
 *
 * \param[in]  cg_ptr      Pointer to the structure previously filled in when the
 *                         counter was allocated.
 * \param[in]  counter_id  Identifier of the counter being disabled 
 *                         (0 to DMA_NUM_COUNTERS_PER_GROUP-1)
 *
 * \return None
 *
 * \note This function does an MBAR after setting the counter to ensure the
 *       writes have been accepted by the memory system before allowing other
 *       memory accesses to to occur.
 *
 */
__INLINE__ void DMA_CounterSetDisableById(  
				DMA_CounterGroup_t *cg_ptr, 
				int                 counter_id
			       )
{
  SPI_assert( (counter_id >= 0) && (counter_id < DMA_NUM_COUNTERS_PER_GROUP) );
  SPI_assert( cg_ptr != NULL );
  SPI_assert( (cg_ptr->permissions[DMA_COUNTER_GROUP_WORD_ID(counter_id)] & 
	   _BN(DMA_COUNTER_GROUP_WORD_BIT_ID(counter_id))) != 0 );
  SPI_assert( cg_ptr->status_ptr != 0);

  /* Disable the counter by writing 1 to the appropriate bit */
  int r =  DMA_COUNTER_GROUP_WORD_ID(counter_id);      
  int c =  DMA_COUNTER_GROUP_WORD_BIT_ID(counter_id);  
  cg_ptr->status_ptr->disable[r] = _BN(c);

  _bgp_mbar();    /* Make sure these writes have been accepted by the memory */
                  /* system before continuing                                */
}


/*!
 * \brief Determine Whether a DMA Counter is Enabled using a Counter ID
 *
 * Determine whether a DMA counter is enabled, given a counter group structure 
 * and counter ID.
 *
 * \param[in]  cg_ptr      Pointer to the structure previously filled in when the
 *                         counter was allocated.
 * \param[in]  counter_id  Identifier of the counter being queried 
 *                         (0 to DMA_NUM_COUNTERS_PER_GROUP-1)
 *
 * \retval  0  The counter is disabled
 * \retval  1  The counter is enabled
 *
 */
__INLINE__ int DMA_CounterGetEnabledById(  
				const DMA_CounterGroup_t *cg_ptr, 
				int                       counter_id
			       )
{
  SPI_assert( (counter_id >= 0) && (counter_id < DMA_NUM_COUNTERS_PER_GROUP) );
  SPI_assert( cg_ptr != NULL );
  SPI_assert( (cg_ptr->permissions[DMA_COUNTER_GROUP_WORD_ID(counter_id)] & 
	   _BN(DMA_COUNTER_GROUP_WORD_BIT_ID(counter_id))) != 0 );
  SPI_assert( cg_ptr->status_ptr != 0);

  /* Return 0 or 1 if counter is disabled/enabled */
  int r =  DMA_COUNTER_GROUP_WORD_ID(counter_id);      
  int c =  DMA_COUNTER_GROUP_WORD_BIT_ID(counter_id);  
  if ( ( cg_ptr->status_ptr->enabled[r] & _BN(c) ) == 0 ) {return 0;}
  else { return 1;}
}


/*!
 * \brief Determine Whether a DMA Counter is Has Hit Zero using a Counter ID
 *
 * Determine whether a DMA counter has hit zero, given a counter group structure 
 * and counter ID.
 *
 * \param[in]  cg_ptr      Pointer to the structure previously filled in when the
 *                         counter was allocated.
 * \param[in]  counter_id  Identifier of the counter being queried 
 *                         (0 to DMA_NUM_COUNTERS_PER_GROUP-1)
 *
 * \retval  0  The counter has not hit zero
 * \retval  1  The counter has hit zero
 *
 * \note This function does an MSYNC after determining that the counter has hit
 *       zero to ensure that the data that was just DMA'd is available to all
 *       cores.  The msync is only done if this is a reception counter group,
 *       since there is nothing to sync for injection counters that have hit zero.
 *
 */
__INLINE__ int DMA_CounterGetHitZeroById(  
				const DMA_CounterGroup_t *cg_ptr, 
				int                       counter_id
			       )
{
  SPI_assert( (counter_id >= 0) && (counter_id < DMA_NUM_COUNTERS_PER_GROUP) );
  SPI_assert( cg_ptr != NULL );
  SPI_assert( (cg_ptr->permissions[DMA_COUNTER_GROUP_WORD_ID(counter_id)] & 
	   _BN(DMA_COUNTER_GROUP_WORD_BIT_ID(counter_id))) != 0 );
  SPI_assert( cg_ptr->status_ptr != 0);

  /* Return 0 or 1 if counter has hit zero */
  int r =  DMA_COUNTER_GROUP_WORD_ID(counter_id);      
  int c =  DMA_COUNTER_GROUP_WORD_BIT_ID(counter_id);  
  if ( ( cg_ptr->status_ptr->hit_zero[r] & _BN(c) ) == 0 ) {return 0;}
  else {
    /* By convention, we assume that if counter has hit zero, then it will be 
     * used.  This requires an msync to ensure snoops from the DMA arbiter 
     * have hit the cores.  That is, the data that was just DMA'd is available 
     * to all cores.  
     *
     * Furthermore, If we just put a _bgp_msync() here, it could get 
     * speculatively executed and withdrawn even if the counter hasn't hit zero,
     * so we call a special version of this function that ensures the speculation
     * does not occur.
     *
     * It only needs to be done for reception counters since there is nothing
     * to sync when sending data.
     */
    if ( cg_ptr->type == DMA_Type_Reception ) _bgp_msync_nonspeculative();    
    return 1;
  }
}


/*!
 * \brief Clear a DMA Counter's Hit Zero Status using a Counter ID
 *
 * Clear a DMA counter's "hit zero" status, given a counter group structure 
 * and counter ID.
 *
 * \param[in]  cg_ptr      Pointer to the structure previously filled in when the
 *                         counter was allocated.
 * \param[in]  counter_id  Identifier of the counter whose "hit zero" status is 
 *                         being cleared (0 to DMA_NUM_COUNTERS_PER_GROUP-1)
 *
 * \return None
 *
 * \note This function does an MBAR after setting the counter to ensure the
 *       writes have been accepted by the memory system before allowing other
 *       memory accesses to to occur.
 *
 */
__INLINE__ void DMA_CounterClearHitZeroById(   
					DMA_CounterGroup_t *cg_ptr, 
					int                 counter_id
				       )
{
  SPI_assert( (counter_id >= 0) && (counter_id < DMA_NUM_COUNTERS_PER_GROUP) );
  SPI_assert( cg_ptr != NULL );
  SPI_assert( (cg_ptr->permissions[DMA_COUNTER_GROUP_WORD_ID(counter_id)] & 
	   _BN(DMA_COUNTER_GROUP_WORD_BIT_ID(counter_id))) != 0 );
  SPI_assert( cg_ptr->status_ptr != 0);

  /* Clear the hit zero bit of a counter */
  int r =  DMA_COUNTER_GROUP_WORD_ID(counter_id);      
  int c =  DMA_COUNTER_GROUP_WORD_BIT_ID(counter_id);  
  cg_ptr->status_ptr->clear_hit_zero[r] = _BN(c) ;

  _bgp_mbar();    /* Make sure these writes have been accepted by the memory */
                  /* system before continuing                                */
}


/* 
 * ------------------------------------------------------------------------------
 *
 * The following functions manipulate or get the status of multiple counters
 *
 * ------------------------------------------------------------------------------
 */


/*!
 * \brief Enable Multiple DMA Counters
 *
 * Enable multiple DMA counters, given a counter group structure and mask.
 *
 * \param[in]  cg_ptr       Pointer to the structure previously filled in when the
 *                          counter was allocated.
 * \param[in]  reg          Identifies the "word" (0 or 1) of the counters
 *                          being manipulated.  This is the index into the
 *                          enable array.
 * \param[in]  counterBits  Identifies which counters in the "word" are being
 *                          manipulated.
 *
 * \return None
 *
 * \note This function does an MBAR after setting the counter to ensure the
 *       writes have been accepted by the memory system before allowing other
 *       memory accesses to to occur.
 *
 */
__INLINE__ void DMA_CounterSetEnable(
				     DMA_CounterGroup_t *cg_ptr, 
				     int                 reg, 
				     unsigned int        counterBits
				    )
{
  SPI_assert( cg_ptr != NULL );
  SPI_assert( ( ( reg == 0 ) || ( reg == 1 ) ) );
  SPI_assert( counterBits == (counterBits & cg_ptr->permissions[reg]) );
  SPI_assert( cg_ptr->status_ptr != 0);

  cg_ptr->status_ptr->enable[reg] = counterBits;

  _bgp_mbar();    /* Make sure these writes have been accepted by the memory */
                  /* system before continuing                                */
}


/*!
 * \brief Disable Multiple DMA Counters
 *
 * Disable multiple DMA counters, given a counter group structure and mask.
 *
 * \param[in]  cg_ptr       Pointer to the structure previously filled in when the
 *                          counter was allocated.
 * \param[in]  reg          Identifies the "word" (0 or 1) of the counters
 *                          being manipulated.  This is the index into the
 *                          disable array.
 * \param[in]  counterBits  Identifies which counters in the "word" are being
 *                          manipulated.
 *
 * \return None
 *
 * \note This function does an MBAR after setting the counter to ensure the
 *       writes have been accepted by the memory system before allowing other
 *       memory accesses to to occur.
 *
 */
__INLINE__ void DMA_CounterSetDisable(DMA_CounterGroup_t *cg_ptr, 
				      int                 reg, 
				      unsigned int        counterBits
				     )
{
  SPI_assert( cg_ptr != NULL );
  SPI_assert( ( ( reg == 0 ) || ( reg == 1 ) ) );
  SPI_assert( counterBits == (counterBits & cg_ptr->permissions[reg]) );
  SPI_assert( cg_ptr->status_ptr != 0);

  cg_ptr->status_ptr->disable[reg] = counterBits;

  _bgp_mbar();    /* Make sure these writes have been accepted by the memory */
                  /* system before continuing                                */
}


/*!
 * \brief Get Enabled DMA Counters
 *
 * Get the enabled status of DMA counters, given a counter group structure 
 * and "word".
 *
 * \param[in]  cg_ptr       Pointer to the structure previously filled in when the
 *                          counter was allocated.
 * \param[in]  reg          Identifies the "word" (0 or 1) of the counters
 *                          being queried.  This is the index into the
 *                          enabled array.
 *
 * \return  32 bit mask indicating which counters in the specified word are enabled.
 *          Only the counters that the caller has allocated will have their status
 *          returned.  The status for other counters will be 0.
 *
 */
__INLINE__ unsigned int DMA_CounterGetEnabled(
				const DMA_CounterGroup_t *cg_ptr, 
				int                       reg
			       )
{
  SPI_assert( ( ( cg_ptr != NULL ) &&
	    ( ( reg == 0 ) || ( reg == 1 ) ) ) );
  SPI_assert( cg_ptr->status_ptr != 0);

  return (cg_ptr->permissions[reg] & cg_ptr->status_ptr->enabled[reg]);
}


/*!
 * \brief Get Hit Zero Status of DMA Counters
 *
 * Get the "hit zero" status of DMA counters, given a counter group structure 
 * and "word".
 *
 * \param[in]  cg_ptr       Pointer to the structure previously filled in when the
 *                          counter was allocated.
 * \param[in]  reg          Identifies the "word" (0 or 1) of the counters
 *                          being queried.  This is the index into the
 *                          hit zero array.
 *
 * \return  32 bit mask indicating which counters in the specified word hit zero.
 *          Only the counters that the caller has allocated will have their status
 *          returned.  The status for other counters will be 0.
 *
 * \note This function does an MSYNC after determining that the counter has hit
 *       zero to ensure that the data that was just DMA'd is available to all
 *       cores.  The msync is only done if this is a reception counter group,
 *       since there is nothing to sync for injection counters that have hit zero.
 *
 */
__INLINE__ unsigned int DMA_CounterGetHitZero(
					      const DMA_CounterGroup_t *cg_ptr, 
					      int                       reg
					     )
{
  unsigned int x;

  SPI_assert( ( ( cg_ptr != NULL ) &&
	    ( ( reg == 0 ) || ( reg == 1 ) ) ) );
  SPI_assert( cg_ptr->status_ptr != 0);

  x =  cg_ptr->status_ptr->hit_zero[reg];

  if ( x != 0 ) {

    x &= cg_ptr->permissions[reg];

    if ( ( cg_ptr->type == DMA_Type_Reception ) &&
	 ( x != 0 ) )
      _bgp_msync_nonspeculative();   

  }

  return (x);
}


/*!
 * \brief Get Hit Zero Status of All DMA Counters In the Specified Group
 *
 * Get the "hit zero" status of all DMA counters in the group specified by the
 * counter group structure.
 *
 * \param[in]  cg_ptr       Pointer to the structure previously filled in when
 *                          the counter was allocated.
 * \param[in,out]  word0    Pointer to the first status word, for the first 32 
 *                          counters.
 * \param[in,out]  word1    Pointer to the second status word, for the second 32
 *                          counters.
 *
 * \return  word0 and word1 are set to the status of the counters.
 *          Only the counters that the caller has allocated will have their 
 *          status returned.  The status for other counters will be 0.
 *
 * \note This function does an MSYNC after determining that at least 1 counter 
 *       has hit zero to ensure that the data that was just DMA'd is available 
 *       to all cores.  The msync is only done if this is a reception counter 
 *       group, since there is nothing to sync for injection counters that have 
 *       hit zero.
 *
 */
__INLINE__ void DMA_CounterGetAllHitZero(
					 const DMA_CounterGroup_t *cg_ptr, 
					 unsigned int             *word0,
					 unsigned int             *word1
					)
{
  unsigned int x,y;

  SPI_assert( ( cg_ptr != NULL ) &&
	      ( word0  != NULL ) &&
	      ( word1  != NULL ) );
  SPI_assert( cg_ptr->status_ptr != 0 );

  x = cg_ptr->status_ptr->hit_zero[0];
  y = cg_ptr->status_ptr->hit_zero[1];

  if ( (x | y) != 0 ) {
    x &= cg_ptr->permissions[0];
    y &= cg_ptr->permissions[1];
    
    if ( ( cg_ptr->type == DMA_Type_Reception ) &&
	 ( (x | y) != 0 ) )
      _bgp_msync_nonspeculative(); 
  }

  *word0 = x;
  *word1 = y;

  return;
}


/*!
 * \brief Clear Hit Zero Status of DMA Counters
 *
 * Clear the "hit zero" status of DMA counters, given a counter group structure, 
 * a "word", and a mask of counters.
 *
 * \param[in]  cg_ptr       Pointer to the structure previously filled in when the
 *                          counter was allocated.
 * \param[in]  reg          Identifies the "word" (0 or 1) of the counters
 *                          being manipulated.  This is the index into the
 *                          clear_hit_zero array.
 * \param[in]  counterBits  Identifies which counters in the "word" are being
 *                          manipulated.
 *
 * \return  None
 *
 * \note This function does an MBAR after setting the counter to ensure the
 *       writes have been accepted by the memory system before allowing other
 *       memory accesses to to occur.
 *
 */
__INLINE__ void DMA_CounterGroupClearHitZero(
				DMA_CounterGroup_t *cg_ptr, 
				int                 reg,
				unsigned int        counterBits 
			       )
{
  SPI_assert( cg_ptr != NULL );
  SPI_assert( ( ( reg == 0 ) || ( reg == 1 ) ) );
  SPI_assert( counterBits == (counterBits & cg_ptr->permissions[reg]) );
  SPI_assert( cg_ptr->status_ptr != 0);

  cg_ptr->status_ptr->clear_hit_zero[reg] = counterBits;

  _bgp_mbar();    /* Make sure these writes have been accepted by the memory */
                  /* system before continuing                                */
}


/*!
 * \brief Get DMA Counter Group Status
 *
 * Get the DMA Counter Group Status, given a counter group structure.
 *
 * \param[in]  cg_ptr       Pointer to the structure previously filled in when the
 *                          counters were allocated.
 *
 * \return  32 bit mask indicating which subgroups have counters that are enabled and
 *          have hit zero.  Only the subgroups that the caller has allocated will have 
 *          their status returned.  The status for other subgroups will be 0.
 *
 * \note This function does an MSYNC after determining that the counter has hit
 *       zero to ensure that the data that was just DMA'd is available to all
 *       cores.  The msync is only done if this is a reception counter group,
 *       since there is nothing to sync for injection counters that have hit zero.
 *
 */
__INLINE__ unsigned int DMA_CounterGetGroupStatus(
				const DMA_CounterGroup_t *cg_ptr
			       )
{
  unsigned int x;

  SPI_assert( cg_ptr != NULL );
  SPI_assert( cg_ptr->status_ptr != 0);

  x = cg_ptr->status_ptr->grp_status;

  if ( x != 0 ) {

    x &= cg_ptr->grp_permissions;

    if (  ( cg_ptr->type == DMA_Type_Reception ) &&
	  ( x != 0 ) )
      _bgp_msync_nonspeculative();

  }

  return x;
}


/*!
 * \brief Get DMA Counter Group Number
 *
 * Get the DMA Counter Group number, given a counter group structure.
 *
 * \param[in]  cg_ptr       Pointer to the structure previously filled in when the
 *                          counters were allocated.
 *
 * \return  The DMA Counter Group number
 *
 */
__INLINE__ int DMA_CounterGetGroupNum(
				      const DMA_CounterGroup_t *cg_ptr
				     )
{
  SPI_assert( cg_ptr != NULL );

  return cg_ptr->group_id;
}


/*!
 * \brief Get DMA Counter Global Id
 *
 * Get the global Id of a DMA Counter, given a counter group structure and a counter Id.
 *
 * \param[in]  cg_ptr      Pointer to the structure previously filled in when the
 *                         counters were allocated.
 * \param[in]  counter_id  Identifier of the counter
 *
 * \return  The DMA Counter Global Id (0 to DMA_NUM_COUNTERS-1)
 *
 */
__INLINE__ int DMA_CounterGetGlobalId( 
				      const DMA_CounterGroup_t *cg_ptr,   
				      int                       counter_id
				     )
{
  SPI_assert( ( cg_ptr != NULL ) &&
	  ( counter_id >= 0 ) && 
	  ( counter_id < DMA_NUM_COUNTERS_PER_GROUP ) );

  return( counter_id + (DMA_NUM_COUNTERS_PER_GROUP * cg_ptr->group_id) );
}


/*!
 * \brief Get DMA Counter Local Id
 *
 * Get the local Id of a DMA Counter, given a counter group structure and a Global 
 * counter Id.
 *
 * \param[in]  counter_id  Global Identifier of the counter (0 to DMA_NUM_COUNTERS-1)
 *
 * \return  The DMA Counter Local Id (0 to DMA_NUM_COUNTERS_PER_GROUP-1)
 *
 */
__INLINE__ int DMA_CounterGetLocalId(  
				     int counter_id
				    )
{
  return( counter_id % DMA_NUM_COUNTERS_PER_GROUP );
}




__END_DECLS


#endif 
