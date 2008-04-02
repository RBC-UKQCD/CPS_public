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

#ifndef	_DMA_FIFO_H_ /* Prevent multiple inclusion */
#define	_DMA_FIFO_H_


/*!
 * \file spi/DMA_Fifo.h
 *
 * \brief DMA SPI Fifo Definitions and Inline Functions Common to Injection
 *        and Reception Fifos
 *
 * This include file contains data structures and inline functions that are 
 * common among injection and reception fifos.  The inlines are used to 
 * interface with the fifos at the lowest level.
 *
 * There are two levels of access:  hardware and software.  For direct 
 * hardware access, the DMA_FifoHW_t structure describes fields that reside
 * in the "hardware fifo" in DMA SRAM.  For normal software access, the 
 * DMA_Fifo_t structure contains a pointer to the hardware structure,
 * shadows (snapshot copies) of the fields in the hardware structure, and
 * size information calculated from the shadows.
 * 
 * \verbatim Picture of fifo structures
  
   ========DDR MEMORY===================|==========DMA SRAM MEMORY==========
   ------------------------------       |
   | DMA_Fifo_t                 |       |
   |                            |       |
   |   Software Fifo            |       |
   |                            |       |
   |                            |       |     -----------------------------
   |   fifo_hw_ptr--------------|-------|---->| DMA_FifoHW_t              |
   |                            |       |     |                           |
   |                            |       |     |   Hardware Fifo           |
   |   Shadow Pointers          |       |     -----------------------------
   |             .              |       |
   ------------------------------       |
  
   \endverbatim
 *
 * For normal messaging software, one should access the DMA using the
 * DMA_Fifo_t, DMA_InjFifo_t, or DMA_RecFifo_t structures since 
 * they maintain shadows.  This include file contains inline functions that 
 * operate on the DMA_Fifo_t for this purpose.  Functions include:
 * - get va_start, va_head, va_tail, va_end, fifo size, fifo free_space
 * - set va_head, va_tail
 * - update fifo free-space based upon current shadows
 *
 * However, for bringup or diagnostic software, there is a need for direct
 * access to the hardware fifos.  This include file contains functions that
 * operate on the DMA_FifoHW_t for this purpose.  Functions include:
 * - get pa_start, pa_head, pa_tail, pa_end
 * - set pa_start, pa_head, pa_tail, pa_end
 * While it probably doesn't make sense to have a stand-alone 
 * DMA_FifoSetStartPa() or DMA_FifoSetEndPa() since this dynamically 
 * messes up the fifo, causing unpredictable results.  But bringup or 
 * diagnostic software will need this (with dma disabled, or the fifo 
 * disabled).  Therefore we provide direct interfaces using physical 
 * addresses and no shadows (for speed).
 *
 * Definitions:
 * - A fifo represents a contiguous block of DDR memory
 * - A fifo has a starting address and an ending address (defines the memory
 *   block)
 * - An injection fifo is a series of 32-byte descriptors.
 * - Injection consists of copying a 32-byte descriptor into the next available
 *   slot (pointed to by the tail) and incrementing the tail pointer.  
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
#include <spi/kernel_interface.h>



/*!
 * \brief Number of fifo groups
 */
#define DMA_NUM_FIFO_GROUPS 4


/*!
 * \brief Hardware DMA Fifo
 *
 * This maps the hardware fifo (the DMA SRAM) for a fifo.  These fields are
 * common for injection and reception fifos.
 *
 * The fifo represents a physically contiguous block of memory.
 *
 */
typedef struct DMA_FifoHW_t
{
  volatile unsigned  pa_start; /*!< RW fifo start address.                      
                                       16B-aligned 4-bit shifted address.     */

  volatile unsigned  pa_end;   /*!< RW fifo end address.                        
                                       16B-aligned 4-bit shifted address.     */

  volatile unsigned  pa_head;  /*!< RW fifo head pointer.                       
                                       16B-aligned 4-bit shifted address.       
                                       Injection fifo head moved by DMA.        
                                       Reception fifo head moved by cores.      
                                       Remote get injection fifo head moved     
                                       by DMA.                                */

  volatile unsigned  pa_tail;  /*!< RW fifo tail pointer.                       
                                       16B-aligned 4-bit shifted address.       
                                       Injection fifo tail moved by cores.      
                                       Reception fifo tail moved by DMA.        
                                       Remote get injection fifo tail moved     
                                       by DMA.                                */
} 
DMA_FifoHW_t;


/*!
 * \brief Software DMA Fifo structure
 *
 * This structure contains a pointer to the hardware fifo, and other fields that
 * describe software's view of the fifo.  These fields are common for injection
 * and reception fifos.
 *
 * \todo   Some more careful thought should be given how to group these so as to
 *         get best memory system performance.
 *         eg.  Probably want to ALIGN_L3_CACHE the fifo_hw_ptr.
 *         Might want to have an assert to check that sizeof( DMA_Fifo_t)
 *         is 32.
 *         COMMENT:  I think below definition puts the entire structure in one
 *                   L1 line.
 */
typedef struct DMA_Fifo_t
{
  DMA_FifoHW_t *fifo_hw_ptr;     /*!< Pointer to hardware fifo.               */

  unsigned int free_space;       /*!< Shadow of how much free space is in the   
                                      fifo, in units of 16B quads.            */

  unsigned int fifo_size;        /*!< Shadow of how much total space is in the  
                                      fifo, in units of 16B quads.            */

  unsigned int pa_start;         /*!< Physical address of the start. (shadow)   
                                      16B-aligned 4-bit shifted address.        
                                      Enables simple calculation of va_head,    
                                      va_tail, and va_end.                    */  
  /*! 
   * \note  The following 4 fields are shadows of the hardware fifo.
   *        They should be in the same L1 cache line for performance.
   *        They are updated by the inline functions in this file upon each
   *        read or write to the fifo.
   */
  void *va_start;                /*!< Shadow of the virtual address start of    
                                      the fifo.  Must be 32B aligned.         */

  void *va_head;                 /*!< Shadow of the virtual address head of     
                                      the fifo.                               */  

  void *va_tail;                 /*!< Shadow of the virtual address tail of     
                                      the fifo.                               */       

  void *va_end;                  /*!< Shadow of the virtual address end  of     
                                      the fifo.  Must be 32B aligned.         */

}
/*! 
 * With above, there should be 8 fields x 4 bytes/field = 32 bytes in the 
 * structure.  Below ensures these 32 bytes are in the same cache line.
 */
ALIGN_L1D_CACHE DMA_Fifo_t;

/*
 *------------------------------------------------------------------------------
 * The following functions operate on fields in the hardware and software fifo 
 * structures.
 *------------------------------------------------------------------------------
 */


/*!
 * \brief Update DMA Fifo Free Space from the Shadow
 *
 * Force a recalculation of a DMA fifo's amount of free space, given a software 
 * fifo structure.
 *
 * \param[in]  f_ptr      Pointer to the software fifo structure
 *
 * \return  None
 *
 * \note  WARNING:  The calculation is based on the current shadow values of the
 *                  head and tail, not the actual hardware values.
 *
 */
__INLINE__ void DMA_FifoUpdateFreeSpaceFromShadow( 
						  DMA_Fifo_t *f_ptr
						 )
{
  SPI_assert( f_ptr           != NULL );
  SPI_assert( f_ptr->fifo_hw_ptr != NULL );

  /* 
   * Recompute the amount of free space in the fifo, given the current shadows.
   */
 
  if ( f_ptr->va_tail >= f_ptr->va_head) 
    { 
      f_ptr->free_space = f_ptr->fifo_size - 
	                    ( ( (unsigned)(f_ptr->va_tail) - 
			        (unsigned)(f_ptr->va_head) ) >> 4 );
    }
  else 
    {
      f_ptr->free_space = ( (unsigned)(f_ptr->va_head) - 
                            (unsigned)(f_ptr->va_tail) ) >> 4;
    }

}


/*!
 * \brief Get DMA Fifo Start Virtual Address from the Shadow
 *
 * Get a DMA fifo's "start" virtual address, given a software fifo structure
 *
 * \param[in]  f_ptr  Pointer to the software fifo structure
 *
 * \retval  va_start  The virtual address of the start of the fifo
 *
 * \note WARNING: This function does not read the DMA SRAM, but instead returns
 *                the current shadow va_start.  To actually issue a read to the
 *                DMA, use DMA_FifoGetStartPa().
 */
__INLINE__ void * DMA_FifoGetStartFromShadow(  
					     DMA_Fifo_t *f_ptr
					    )
{
  SPI_assert( f_ptr              != NULL );
  SPI_assert( f_ptr->fifo_hw_ptr != NULL );

  return  f_ptr->va_start;   
} 


/*!
 * \brief Get DMA Fifo Head Virtual Address
 *
 * Get a DMA fifo's "head" virtual address, given a software fifo structure
 *
 * \param[in]  f_ptr  Pointer to the software fifo structure
 *
 * \retval  va_head  The virtual address of the head of the fifo
 *
 * \post va_head is recalculated from the current hardware head, updated in
 *       the software fifo structure, and returned.  Additionally, the free
 *       space in the software fifo structure is updated.
 *
 */
__INLINE__ void * DMA_FifoGetHead(  
				  DMA_Fifo_t *f_ptr
				 )
{
  unsigned int val;

  SPI_assert( f_ptr              != NULL );
  SPI_assert( f_ptr->fifo_hw_ptr != NULL );

  /* Read the DMA to get the head.
   * Recompute va_head based upon the va_start and the current hardware head.
   * Update free_space.
   */
  
  val = f_ptr->fifo_hw_ptr->pa_head;

  f_ptr->va_head = (char*)( (unsigned)f_ptr->va_start + 
			    ( ( val - f_ptr->pa_start ) << 4 ) );

  DMA_FifoUpdateFreeSpaceFromShadow( f_ptr );

  return  f_ptr->va_head;

} 


/*!
 * \brief Get DMA Fifo Head Virtual Address Without Updating Free Space
 *
 * Get a DMA fifo's "head" virtual address, given a software fifo structure,
 * without updating the fifo's free space.  It is up to the caller to ensure
 * this update occurs later, if necessary.
 *
 * \param[in]  f_ptr  Pointer to the software fifo structure
 *
 * \retval  va_head  The virtual address of the head of the fifo
 *
 * \post va_head is recalculated from the current hardware head, updated in
 *       the software fifo structure, and returned.
 *
 */
__INLINE__ void * DMA_FifoGetHeadNoFreeSpaceUpdate(  
						   DMA_Fifo_t *f_ptr
						  )
{
  unsigned int val;

  SPI_assert( f_ptr              != NULL );
  SPI_assert( f_ptr->fifo_hw_ptr != NULL );

  /* Read the DMA to get the head.
   * Recompute va_head based upon the va_start and the current hardware head.
   */
  
  val = f_ptr->fifo_hw_ptr->pa_head;

  f_ptr->va_head = (char*)( (unsigned)f_ptr->va_start + 
			    ( ( val - f_ptr->pa_start ) << 4 ) );

  return  f_ptr->va_head;

} 


/*!
 * \brief Get DMA Fifo Tail Virtual Address
 *
 * Get a DMA fifo's "tail" virtual address, given a software fifo structure
 *
 * \param[in]  f_ptr  Pointer to the software fifo structure
 *
 * \retval  va_tail  The virtual address of the tail of the fifo
 *
 * \post va_tail is recalculated from the current hardware tail, updated in
 *       the software fifo structure, and returned.  Additionally, the free
 *       space in the software fifo structure is updated.
 *
 */
__INLINE__ void * DMA_FifoGetTail(  
				  DMA_Fifo_t *f_ptr
				 )
{
  unsigned int val;

  SPI_assert( f_ptr              != NULL );
  SPI_assert( f_ptr->fifo_hw_ptr != NULL );

  /* Read the DMA to get the tail.
   * Recompute va_tail based upon the va_start and the current hardware tail.
   * Update free_space.
   */
  
  val = f_ptr->fifo_hw_ptr->pa_tail;

  f_ptr->va_tail = (char*)( (unsigned)f_ptr->va_start + 
			    ( ( val - f_ptr->pa_start ) << 4 ) );

  DMA_FifoUpdateFreeSpaceFromShadow( f_ptr );

  return  f_ptr->va_tail;

 
} 


/*!
 * \brief Get DMA Fifo Tail Virtual Address Without Updating Free Space
 *
 * Get a DMA fifo's "tail" virtual address, given a software fifo structure,
 * without updating the fifo's free space.  It is up to the caller to 
 * invoke DMA_FifoUpdateFreeSpaceFromShadow() at a later time.
 *
 * \param[in]  f_ptr  Pointer to the software fifo structure
 *
 * \retval  va_tail  The virtual address of the tail of the fifo
 *
 * \post va_tail is recalculated from the current hardware tail, updated in
 *       the software fifo structure, and returned.  
 *
 */
__INLINE__ void * DMA_FifoGetTailNoFreeSpaceUpdate(  
						   DMA_Fifo_t *f_ptr
						  )
{
  unsigned int val;

  SPI_assert( f_ptr              != NULL );
  SPI_assert( f_ptr->fifo_hw_ptr != NULL );

  /* Read the DMA to get the tail.
   * Recompute va_tail based upon the va_start and the current hardware tail.
   */
  
  val = f_ptr->fifo_hw_ptr->pa_tail;

  f_ptr->va_tail = (char*)( (unsigned)f_ptr->va_start + 
			    ( ( val - f_ptr->pa_start ) << 4 ) );

  return  f_ptr->va_tail;
 
} 


/*!
 * \brief Get DMA Fifo Tail Virtual Address from Shadow
 *
 * Get a DMA fifo's "tail" virtual address, given a software fifo structure
 *
 * \param[in]  f_ptr  Pointer to the software fifo structure
 *
 * \retval  va_tail  The virtual address of the tail of the fifo
 *
 * \post va_tail is obtained from the shadow, NOT recalculated from the 
 *       current hardware tail.  The free space in the software fifo
 *       structure is NOT updated.
 *
 */
__INLINE__ void * DMA_FifoGetTailFromShadow(  
					    DMA_Fifo_t *f_ptr
					   )
{
  SPI_assert( f_ptr              != NULL );
  SPI_assert( f_ptr->fifo_hw_ptr != NULL );

  return  f_ptr->va_tail;
 
} 


/*!
 * \brief Get DMA Fifo End Virtual Address from the Shadow
 *
 * Get a DMA fifo's "end" virtual address, given a software fifo structure
 *
 * \param[in]  f_ptr  Pointer to the software fifo structure
 *
 * \retval  va_end  The virtual address of the end of the fifo
 *
 * \note WARNING: This function does not read the DMA SRAM, but instead returns
 *                the current shadow va_end.  To actually issue a read to the
 *                DMA, use DMA_FifoGetEndPa().
 */
__INLINE__ void * DMA_FifoGetEndFromShadow(  
					   DMA_Fifo_t *f_ptr
					  )
{
  SPI_assert( f_ptr              != NULL );
  SPI_assert( f_ptr->fifo_hw_ptr != NULL );

  return  f_ptr->va_end;
} 


/*!
 * \brief Get DMA Fifo Size
 *
 * Get a DMA fifo's size, given a software fifo structure
 *
 * \param[in]  f_ptr  Pointer to the software fifo structure
 *
 * \retval  size  The size of the DMA fifo, in units of 16B quads.
 *
 * \note WARNING: This function does not calculate the size based on the DMA
 *                SRAM's current start and end values, but instead returns the
 *                size that was calculated when the fifo was initialized.
 */
__INLINE__ unsigned int DMA_FifoGetSize(  
					DMA_Fifo_t *f_ptr
				       )
{
  SPI_assert( f_ptr              != NULL );
  SPI_assert( f_ptr->fifo_hw_ptr != NULL );

  return  f_ptr->fifo_size;
} 


/*!
 * \brief Get DMA Fifo Free Space With No Update Calculation
 *
 * Get a DMA fifo's amount of free space, given a software fifo structure.
 * Do not perform update calculations.
 *
 * \param[in]  f_ptr      Pointer to the software fifo structure
 *
 * \retval  freeSpace  The amount of free space in the fifo, in units of 
 *                     16B quads.
 */
__INLINE__ unsigned int DMA_FifoGetFreeSpaceNoUpdateCalculation(
                                             DMA_Fifo_t   *f_ptr
                                            )
{
  SPI_assert( f_ptr              != NULL );

  return f_ptr->free_space;
}


/*!
 * \brief Get DMA Fifo Free Space
 *
 * Get a DMA fifo's amount of free space, given a software fifo structure
 *
 * \param[in]  f_ptr      Pointer to the software fifo structure
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
 * \note If both read_head and read_tail are false, the amount of free space is
 *       calculated based on the current shadow values of head and tail.
 */
__INLINE__ unsigned int DMA_FifoGetFreeSpace(  
					     DMA_Fifo_t   *f_ptr, 
					     unsigned int  read_head, 
					     unsigned int  read_tail 
					    )
{
  SPI_assert( f_ptr              != NULL );
  SPI_assert( f_ptr->fifo_hw_ptr != NULL );
  SPI_assert( read_head == 1 || read_head == 0 );
  SPI_assert( read_tail == 1 || read_tail == 0 );

  /*
   * If both read_head and read_tail are 0, return the current shadow.
   * If read_head != 0, read the head of the fifo first and recompute free space.
   * If read_tail != 0, read the tail of the fifo first and recompute free space.
   */

  if ( (read_head == 0) && ( read_tail == 0) )  
    DMA_FifoUpdateFreeSpaceFromShadow( f_ptr);
  else
    {
      if ( read_head == 1) DMA_FifoGetHead(f_ptr);    /* This does an update    */
                                                      /* of the free space.     */
      if ( read_tail == 1) DMA_FifoGetTail(f_ptr);    /* This does an update    */
                                                      /* of the free space.     */
    }

  return f_ptr->free_space;
  
}


/*!
 * \brief Set DMA Fifo Head
 *
 * Set a DMA fifo's "head", given a software fifo structure
 *
 * \param[in]  f_ptr    Pointer to the software fifo structure
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
__INLINE__ void DMA_FifoSetHead(  
				DMA_Fifo_t *f_ptr,
				void       *va_head 
			       )
{
  unsigned int pa_head;

  SPI_assert( f_ptr              != NULL );
  SPI_assert( f_ptr->fifo_hw_ptr != NULL );
  SPI_assert( va_head >= f_ptr->va_start && 
	   va_head <  f_ptr->va_end );

  /* 
   * Calculate new pa_head based on the shadow pa_start and va_start.
   */
  pa_head = f_ptr->pa_start + ( ( (unsigned)va_head - 
				  (unsigned)f_ptr->va_start ) >> 4 );

  /*
   * Set the hardware head
   */
  f_ptr->fifo_hw_ptr->pa_head = pa_head;
  _bgp_mbar();

  /*
   * Update the software fifo structure's head and free space.
   */
  f_ptr->va_head             = va_head;

  DMA_FifoUpdateFreeSpaceFromShadow( f_ptr );

} 


/*!
 * \brief Set DMA Fifo Tail
 *
 * Set a DMA fifo's "tail", given a software fifo structure
 *
 * \param[in]  f_ptr    Pointer to the software fifo structure
 * \param[in]  va_tail  Virtual address of the tail to be set
 *
 * \return  None
 *
 * \post va_tail is set in both the hardware and software fifo structures,
 *       and the fifo free space is recalculated.
 *
 */
__INLINE__ void DMA_FifoSetTail(  
				DMA_Fifo_t *f_ptr,
				void       *va_tail 
			       )
{
  unsigned int pa_tail;

  SPI_assert( f_ptr              != NULL );
  SPI_assert( f_ptr->fifo_hw_ptr != NULL );
  SPI_assert( va_tail >= f_ptr->va_start && 
	   va_tail <  f_ptr->va_end );

  /* 
   * Calculate new pa_tail based on the shadow pa_start and va_start.
   */
  pa_tail = f_ptr->pa_start + ( ( (unsigned)va_tail - 
				  (unsigned)f_ptr->va_start ) >> 4 );

  /*
   * Set the hardware tail
   */
  f_ptr->fifo_hw_ptr->pa_tail = pa_tail;
  _bgp_mbar();

  /*
   * Update the software fifo structure's tail and free space.
   */
  f_ptr->va_tail             = va_tail;

  DMA_FifoUpdateFreeSpaceFromShadow( f_ptr );

} 




/*
 *------------------------------------------------------------------------------
 * The following functions operate directly on the hardware fifo.  Normally,
 * users should use the software fifo routines (previously defined), but for
 * bringup or diagnostics, it may be desirable to use these.
 *------------------------------------------------------------------------------
 */




/*!
 * \brief Set DMA Hardware Fifo Start
 *
 * Set a DMA fifo's "start", given a hardware fifo structure
 *
 * \param[in]  fifo_hw_ptr  Pointer to the hardware fifo structure
 * \param[in]  pa_start     Physical address of the start to be set.
 *                          16B-aligned 4-bit shifted physical address.
 *
 * \return  None
 *
 * \note This function does an MBAR after setting the fifo to ensure the
 *       writes have been accepted by the memory system before allowing other
 *       memory accesses to to occur.
 */
__INLINE__ void DMA_FifoSetStartPa(  
				   DMA_FifoHW_t *fifo_hw_ptr,
				   unsigned int  pa_start 
				  )
{
  SPI_assert( fifo_hw_ptr != NULL );

  fifo_hw_ptr->pa_start = pa_start;

  _bgp_mbar();

}
 

/*!
 * \brief Set DMA Hardware Fifo Head
 *
 * Set a DMA fifo's "head", given a hardware fifo structure
 *
 * \param[in]  fifo_hw_ptr  Pointer to the hardware fifo structure
 * \param[in]  pa_head      Physical address of the head to be set.
 *                          16B-aligned 4-bit shifted physical address.
 *
 * \return  None
 *
 * \note This function does an MBAR after setting the fifo to ensure the
 *       writes have been accepted by the memory system before allowing other
 *       memory accesses to to occur.
 */
__INLINE__ void DMA_FifoSetHeadPa(  
				  DMA_FifoHW_t *fifo_hw_ptr,
				  unsigned int  pa_head 
				 )
{
  SPI_assert( fifo_hw_ptr != NULL );

  fifo_hw_ptr->pa_head = pa_head;

  _bgp_mbar();  

}
 

/*!
 * \brief Set DMA Hardware Fifo Tail
 *
 * Set a DMA fifo's "tail", given a hardware fifo structure
 *
 * \param[in]  fifo_hw_ptr  Pointer to the hardware fifo structure
 * \param[in]  pa_tail      Physical address of the tail to be set.
 *                          16B-aligned 4-bit shifted physical address.
 *
 * \return  None
 *
 * \note This function does an MBAR after setting the fifo to ensure the
 *       writes have been accepted by the memory system before allowing other
 *       memory accesses to to occur.
 */
__INLINE__ void DMA_FifoSetTailPa(  
				  DMA_FifoHW_t *fifo_hw_ptr,
				  unsigned int  pa_tail 
				 )

{
  SPI_assert( fifo_hw_ptr != NULL );

  fifo_hw_ptr->pa_tail = pa_tail;

  _bgp_mbar();  

}
 

/*!
 * \brief Set DMA Hardware Fifo End
 *
 * Set a DMA fifo's "end", given a hardware fifo structure
 *
 * \param[in]  fifo_hw_ptr  Pointer to the hardware fifo structure
 * \param[in]  pa_end       Physical address of the end to be set.
 *                          16B-aligned 4-bit shifted physical address.
 *
 * \return  None
 *
 * \note This function does an MBAR after setting the fifo to ensure the
 *       writes have been accepted by the memory system before allowing other
 *       memory accesses to to occur.
 */
__INLINE__ void DMA_FifoSetEndPa(  
				 DMA_FifoHW_t *fifo_hw_ptr,
				 unsigned int  pa_end 
				)
{
  SPI_assert( fifo_hw_ptr != NULL );

  fifo_hw_ptr->pa_end = pa_end;

  _bgp_mbar();  

}


/*!
 * \brief Get DMA Hardware Fifo Start
 *
 * Get a DMA fifo's "start", given a hardware fifo structure
 *
 * \param[in]  fifo_hw_ptr  Pointer to the hardware fifo structure
 *
 * \retval  pa_start  Physical address of the fifo start.
 *                    16B-aligned 4-bit shifted physical address.
 *
 * \return  None
 *
 */
__INLINE__ unsigned int DMA_FifoGetStartPa(  
					   DMA_FifoHW_t *fifo_hw_ptr
					  )
{
  SPI_assert( fifo_hw_ptr != NULL );

  return fifo_hw_ptr->pa_start;
}


/*!
 * \brief Get DMA Hardware Fifo Head
 *
 * Get a DMA fifo's "head", given a hardware fifo structure
 *
 * \param[in]  fifo_hw_ptr  Pointer to the hardware fifo structure
 *
 * \retval  pa_head  Physical address of the fifo head.
 *                   16B-aligned 4-bit shifted physical address.
 *
 * \return  None
 *
 */
__INLINE__ unsigned int DMA_FifoGetHeadPa(  
					  DMA_FifoHW_t *fifo_hw_ptr
					 )
{
  SPI_assert( fifo_hw_ptr != NULL );

  return fifo_hw_ptr->pa_head;
}


/*!
 * \brief Get DMA Hardware Fifo Tail
 *
 * Get a DMA fifo's "tail", given a hardware fifo structure
 *
 * \param[in]  fifo_hw_ptr  Pointer to the hardware fifo structure
 *
 * \retval  pa_tail  Physical address of the fifo tail.
 *                   16B-aligned 4-bit shifted physical address.
 *
 * \return  None
 *
 */
__INLINE__ unsigned int DMA_FifoGetTailPa(  
					  DMA_FifoHW_t *fifo_hw_ptr
					 )
{
  SPI_assert( fifo_hw_ptr != NULL );

  return fifo_hw_ptr->pa_tail;
}


/*!
 * \brief Get DMA Hardware Fifo End
 *
 * Get a DMA fifo's "end", given a hardware fifo structure
 *
 * \param[in]  fifo_hw_ptr  Pointer to the hardware fifo structure
 *
 * \retval  pa_end  Physical address of the fifo end.
 *                  16B-aligned 4-bit shifted physical address.
 *
 * \return  None
 *
 */
__INLINE__ unsigned int DMA_FifoGetEndPa(
					 DMA_FifoHW_t *fifo_hw_ptr
					)
{
  SPI_assert( fifo_hw_ptr != NULL );

  return fifo_hw_ptr->pa_end;
}


__END_DECLS


#endif 
