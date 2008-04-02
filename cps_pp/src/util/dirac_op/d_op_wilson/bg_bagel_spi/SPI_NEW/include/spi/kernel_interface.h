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
 * \file spi/kernel_interface.h
 */

#ifndef _BGP_VIRT2PHYS_H_ // Prevent multiple inclusion
#define _BGP_VIRT2PHYS_H_



#include <common/namespace.h>

__BEGIN_DECLS

#include <common/linkage.h>
#include <bpcore/bgp_types.h>
#include <bpcore/ppc450_core.h>
#include <bpcore/ppc450_inlines.h>
#include <spi/bpcore_interface.h>
#include <spi/bgp_kernel_inlines.h>
#include <common/bgp_ras.h>
#include <cnk/bgp_VirtualMap.h>
#include <cnk/bgp_vmm.h>
#include <cnk/bgp_SPRG_Usage.h>
#include <cnk/bgp_SysCall_Extensions.h>
#include <fcntl.h>
#include <limits.h>
#include <string.h>
#include <errno.h>

#ifdef __CNK__
extern int _bgp_DMA_ClearFullReceptionFifo(void);
#endif

#if ((!defined(__CNK__)) && (!defined(__BL__)))
#include <pthread.h>
#endif

#ifndef __INLINE__
#define __INLINE__ extern inline
#endif

								

/*!
 * \brief Communication Thread interrupt handler function prototype
 *
 * \param[in] arg1 1st argument to commthread
 * \param[in] arg2 2nd argument to commthread
 * \param[in] arg3 3rd argument to commthread
 */
typedef void (*Kernel_CommThreadHandler)(uint32_t arg1, uint32_t arg2, uint32_t arg3, uint32_t arg4);

/*!
 * \brief Interrupt Group Prototype
 *
 * This data type is used to group interrupts of various devices together
 * so they can be enabled or disabled simultaneously.  A given interrupt user
 * (eg. messaging, QCD, etc) specifies a value of this data type when its
 * interrupt resources are allocated.  The kernel associates those resources
 * with the specified value so when this value is specified on the enable or
 * disable interupts system call, all of the interrupts in the group are
 * operated upon.  Examples of devices that can be grouped in this way include
 * DMA fifos, torus, tree, etc.
 *
 * \todo The kernel should provide interfaces to allocate a
 *       Kernel_InterruptGroup_t and deallocate it.
 */
typedef void * Kernel_InterruptGroup_t;




/*! \brief Returns the physical processor ID of the running PPC450 core.
 *
 * \return Physical processor ID
 * \retval 0 Running on processor 0
 * \retval 1 Running on processor 1
 * \retval 2 Running on processor 2
 * \retval 3 Running on processor 3
 */
//__INLINE__ int Kernel_PhysicalProcessorID( void )
//{
//  uint32_t dst2  = _bgp_mfspr( _BGP_SPRGRO_DST2  );

//  return( dst2 & 0x3 );
//}

/*! \brief Returns the number of Processes (Virtual Nodes) running on this Physical Node.
 *
 * \return Process Count
 * \retval 1 Running in Single Process "SMP Mode"
 * \retval 2 Running in "2 Virtual Node Mode"
 * \retval 3 Running in "3 Virtual Node Mode"
 * \retval 4 Running in "4 Virtual Node Mode"
 */
__INLINE__ int Kernel_ProcessCount( void )
{
  uint32_t shm  = _bgp_mfspr( _BGP_SPRGRO_SHMem  );

  return( (shm & 0x3) + 1 );
}

/*! \brief Returns the number of Processors (cores) running in this Process (Virtual Node)
 *
 * \return Processor Count
 * \retval 1 Single Processor in this Process (usually 4-VN Mode).
 * \retval 2 Two Processors in this Process (usually 2-VN Mode).
 * \retval 3 Three Processors in this Process.
 * \retval 4 Four Processors in this Process (usually SMP Mode).
 */
__INLINE__ int Kernel_ProcessorCount( void )
{
  uint32_t shm  = _bgp_mfspr( _BGP_SPRGRO_SHMem  );

  return( ((shm & 0xC) >> 2) + 1 );
}

__INLINE__ int Kernel_GetAppSegmentCount(uint32_t* count)
{
   _BGP_SprgShMem shm;
   
   shm.shmem = _bgp_mfspr(_BGP_SPRGRO_SHMem);
   if(shm.IsStaticMap)
   {
      if(Kernel_ProcessCount() == 1)
      {
	 *count = 3;  /* text/rodata, data, heap */
      }
      else
      {
	 *count = 4;  /* text/rodata, data, heap, shared (in dual/vn) */
      }
   }
   else
   {
      if(Kernel_ProcessCount() == 1)
      {
	 *count = 2;  /* text/rodata, data/heap */
      }
      else
      {
	 *count = 3;  /* text/rodata, data/heap, shared (in dual/vn) */
      }
   }
   return 0;
}

__INLINE__ int Kernel_GetAppSegmentMapping(uint32_t segmentID, uint32_t coreID, uint32_t* va, uint64_t* pa, uint32_t* length)
{
   int rc = 0;   
   _BGP_SprgShMem shm;   
   shm.shmem = _bgp_mfspr(_BGP_SPRGRO_SHMem);
   if((!shm.IsStaticMap)&&(segmentID > 1))
      segmentID++;
   
   asm __volatile__ ("li 0,%1;"
		     "mr 3,%2;"
		     "mr 4,%3;"
		     "mr 5,%4;"
		     "mr 6,%5;"
		     "mr 7,%6;"
		     "sc;"
		     "mr %0, 3;"
		     : "=&r" (rc)  // early clobber
		     : "i" (_BGP_SYSCALL_NR_GETAPPSEGMENTMAPPING),
		     "r" (segmentID),
		     "r" (coreID),
		     "r" (va),
		     "r" (pa),
		     "r" (length)
		     : "r0", "r3", "r4", "r5", "r6", "r7", "cc", "memory" );
   return rc;
}

/*! \brief Translate a 32bit Virtual Address to a 36bit Physical Address, returning separated upper and lower parts.
 *
 * \param[in] pVA   32bit virtual address in the calling process
 * \param[in] vsize size in bytes of the virtual range
 * \param[out] ua_out upper 4 physical address bits
 * \param[out] pa_out lower 32 physical address bits
 * \return Error condition for translation
 * \retval  0 Successful translation, with ua_out and pa_out filled in
 * \retval -1 Invalid Virtual Address for this process, ua_out and pa_out unmodified.
 * \retval -2 The range from vaddr to (vaddr+vsize) is not physically contiguous.
 * \retval -3 vaddr in Scratch, but no Scratch, or not enough Scratch, is enabled.
 * \retval -4 invalid parameter
 *
 *  \warning Supports only Text, Data, Stack, and (optional) eDRAM Scratch translation
 *  \warning CNK "pagesize" is 1MB.
 *  \warning Text and Data are virtually contiguous, but not necessarily physically contiguous.
 *  \todo Does not (currently) support > 4GB DDR space.
 *  \todo Does not (currently) support Shared Memory Area.
 */
__INLINE__ int Kernel_Virtual2Physical( void     *pVA,     // input: 32bit Virtual start address
                                    size_t   vsize,    // input: size in bytes of virtual range
                                    uint32_t *ua_out,  // output: upper  4 Physical Address bits
                                    uint32_t *pa_out ) // output: lower 32 Physical Address bits
{
   _BGP_SprgShMem shmem;

   shmem.shmem = _bgp_mfspr(_BGP_SPRGRO_SHMem);
   if(shmem.IsStaticMap)
   {
      uint32_t x;
      static int static_v2p_initialized;
      static uint32_t segcnt;

#define KERNEL_V2P_MAXSEGMENTS 5
      static uint32_t segva[KERNEL_V2P_MAXSEGMENTS];
      static uint64_t segpa[KERNEL_V2P_MAXSEGMENTS];
      static size_t   segsz[KERNEL_V2P_MAXSEGMENTS];
#undef KERNEL_V2P_MAXSEGMENTS
      if(static_v2p_initialized == 0)
      {
	 Kernel_GetAppSegmentCount(&segcnt);
	 for(x=0; x<segcnt; x++)
	 {
	    if(Kernel_GetAppSegmentMapping(x, Kernel_PhysicalProcessorID(), &segva[x], &segpa[x], &segsz[x]))
	       return -1;
	 }
	 static_v2p_initialized = 1;
      }
      for(x=0; x<segcnt; x++)
      {
	 if(((uint32_t)pVA >= segva[x]) && (segsz[x] > (uint32_t)pVA - segva[x] + vsize) && ((uint32_t)pVA + vsize > (uint32_t)pVA))
	 {
	    *ua_out = (uint32_t)((segpa[x] + ((uint32_t)pVA-segva[x])) >> 32);
	    *pa_out = (uint32_t)((segpa[x] + ((uint32_t)pVA-segva[x]))&0xffffffff);
	    return 0;
	 }
      }
      return -1;
   }

   uint32_t vaddr = (uint32_t)pVA;
   uint32_t texti = _bgp_mfspr( _BGP_SPRGRO_TextI );
   uint32_t datai = _bgp_mfspr( _BGP_SPRGRO_DataI );
   uint32_t dst2  = _bgp_mfspr( _BGP_SPRGRO_DST2  );
   uint32_t shm   = (_bgp_mfspr( _BGP_SPRGRO_SHMem ) & 0xFFFFFFC0);
   uint32_t text_v_start = (texti & 0xFFF00000);
   uint32_t data_v_start = (datai & 0xFFF00000); // text_v_limit is (data_v_start - 1)
   uint32_t text_ua      = ((texti & 0x000000C0) >> 6);
   uint32_t text_p_start = ((texti & 0x000FFF00) << 12);
   uint32_t data_ua      = ((datai & 0x000000C0) >> 6);
   uint32_t data_p_start = ((datai & 0x000FFF00) << 12);
   uint32_t data_v_size  = (dst2  & 0xFFF00000);
   uint32_t data_v_limit = (data_v_start + data_v_size + _BGP_VMM_PAGE_MASK);
   uint32_t vend    = (vaddr + vsize - 1);
   uint32_t vpage   = (vaddr & ~_BGP_VMM_PAGE_MASK);  // which 1MB page?
   uint32_t voffset = (vaddr & _BGP_VMM_PAGE_MASK);   // offset within 1MB page

   // printf("V2P: texti=0x%08x, datai=0x%08x, dst2=0x%08x\n", texti, datai, dst2 );
   // printf("V2P: vaddr=0x%08x, vend=0x%08x, text_v_start=0x%08x, data_v_limit=0x%08x\n",
   //         vaddr, vend, text_v_start, data_v_limit );

   // parm check
   if ( !vsize || !ua_out || !pa_out )
      return(-4);

   // range check: below text or off end of data, or in eDRAM Scratch
   if ( (vaddr < text_v_start) || (vend > data_v_limit) )
      {
      // Scratch?
      if ( vaddr >= _BGP_VA_SCRATCH )
         {
         uint32_t scratchMB   = ((dst2 & 0x00000078) << (20-3));
         uint32_t scratch_end = (_BGP_VA_SCRATCH + scratchMB);

         if ( !scratchMB || (vend > scratch_end) )
            return(-3);

         *ua_out = (uint32_t)_BGP_UA_SCRATCH;
         *pa_out = (vaddr & _BGP_VM_SCRATCH);
         return(0);
         }
      else if ( shm ) // Shared Memory? If any, always mapped V=R.
         {
         uint32_t shm_v_start = (shm & 0xFFF00000);
         uint32_t shm_v_end   = (shm_v_start + ((shm & 0x000FFF00) << 12));
         uint32_t shm_ua      = ((shm & 0x000000C0) >> 6);

         if ( (vaddr >= shm_v_start) && (vend <= shm_v_end) )
            {
            *ua_out = shm_ua;
            *pa_out = vaddr;
            return(0);
            }
         }

      return(-1);
      }

   // Text? (includes Read-Only Data)
   if ( vaddr < data_v_start )
      {
      // if range starts in Text but ends in Data, then discontiguous
      if ( vend >= data_v_start )
         return(-2);

      *ua_out = text_ua;
      *pa_out = (text_p_start + (vpage - text_v_start) + voffset);

      return(0);
      }

   // Data
   *ua_out = data_ua;
   *pa_out = (data_p_start + (vpage - data_v_start) + voffset);

   return(0);
}



/*! \brief Returns a copy of the node's personality
 *
 * \param[out] personality Location of personality structure that will be filled in by Kernel_GetPersonality
 * \param[in]  size Size, in bytes, that was allocated to hold the personality structure
 * \return Error indication
 * \retval  0 Success
 * \retval -1 Invalid parameters
 */
__INLINE__ int Kernel_GetPersonality(_BGP_Personality_t* personality, size_t size)
{
  int rc = 0; // this syscall returns RC in r3 and does not use errno
  asm __volatile__ ("li 0,%3;"
		    "mr 3,%1;"
		    "mr 4,%2;"
                     "sc;"
		    "mr %0, 3;"
                    : "=&r" (rc)  // early clobber
                    : "r" (personality),
		      "r" (size),
		      "i" (_BGP_SYSCALL_NR_GET_PERSONALITY)
                    : "r0", "r3", "r4", "cc", "memory" );

  return( rc );
}

/*! \brief Starts to checkpoint/restore the Kernel data structures for CNK
 *
 * \param[out] personality Location of personality structure that will be filled in by Kernel_GetPersonality
 * \param[in]  size Size, in bytes, that was allocated to hold the personality structure
 * \param[in]  int operation, The type of operation that the kernel needs to provide (e.g. CHECKPOINT_START, CHECKPOINT_RESTART,CHECKPOINT_COMPLETE)
 * \return Error indication
 * \retval  0 Success
 * \retval -1 Invalid parameters
 */
__INLINE__ int Kernel_checkpoint(int component, int operation, void *buffer, uint32_t size, uint32_t *actualSize, uint32_t*basePtr)
{
  int rc = 0; // this syscall returns RC in r3 and does not use errno

  asm __volatile__ ("li 0,%1;"
                    "mr 3,%2;"
                    "mr 4,%3;"
                    "mr 5,%4;"
                    "mr 6,%5;"
                    "mr 7,%6;"
                    "mr 8,%7;"
                    "sc;"
                    "mr %0, 3;"
                    : "=&r" (rc)  // early clobber
                    : "i" (_BGP_SYSCALL_NR_CHECKPOINT),
                    "r" (component),
                    "r" (operation),
                    "r" (buffer),
                    "r" (size),
                    "r" (actualSize),
                    "r" (basePtr)
                    : "r0", "r3", "r4", "r5", "r6", "r7", "r8", "cc", "memory" );

  return( rc );
}

/*! \brief Returns the contents of the running PPC450's processor version register.
 * \return Contents of PPC450 PVR register
 */
__INLINE__ int Kernel_GetProcessorVersion()
{
  int rc = 0; // this syscall returns RC in r3 and does not use errno

  asm __volatile__ ("li 0,%1;"
                     "sc;"		
		    "mr %0, 3;"

                    : "=&r" (rc)  // early clobber
                    : "i" (_BGP_SYSCALL_NR_GET_PERSONALITY)
                    : "r0", "r3", "cc", "memory" );

  return( rc );
}

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
//
//  asm __volatile__ ("li 0,%1;"
//		    "mr 3,%2;"
//		    "mr 4,%3;"
//		    "mr 5,%4;"
//		    "mr 6,%5;"
// /                    "sc;"
//		    "mr %0, 3;"
//                    : "=&r" (rc)  // early clobber
//                    : "i" (_BGP_SYSCALL_NR_ALLOC_COUNTER),
//		    "r" (lockid),
//		    "r" (numlocks),
//		    "r" (ptr),
//		    "r" (flags)
// /                   : "r0", "r3", "r4", "r5", "r6", "cc", "memory" );
//
//  return( rc );
//}

/*! \brief Converts a Rank into a XYZT Coordinate
 *
 * \param[in] rank Rank for the node
 * \param[out] xcoord X Coordinate for the specified node
 * \param[out] ycoord Y Coordinate for the specified node
 * \param[out] zcoord Z Coordinate for the specified node
 * \param[out] tcoord T Coordinate for the specified node
 * \return Error status
 * \retval 0 Success
 * \retval non-zero Error
 */
__INLINE__ int Kernel_Rank2Coord(uint32_t rank, uint32_t* xcoord, uint32_t* ycoord, uint32_t* zcoord, uint32_t* tcoord)
{
  int rc = 0; // this syscall returns RC in r3 and does not use errno

  asm __volatile__ ("li 0,%1;"
                    "mr 3,%2;"
                    "mr 4,%3;"
                    "mr 5,%4;"
                    "mr 6,%5;"
		    "mr 7,%6;"
                     "sc;"
		    "mr %0, 3;"
                    : "=&r" (rc)  // early clobber
                    : "i" (_BGP_SYSCALL_NR_RANK2COORD),
                    "r" (rank),
                    "r" (xcoord),
                    "r" (ycoord),
                    "r" (zcoord),
		    "r" (tcoord)
                    : "r0", "r3", "r4", "r5", "r6", "r7", "cc", "memory" );

  return( rc );
}

/*! \brief Convert a XYZT Coordinate into a Rank.  Also returns number of nodes
 * \param[in] xcoord X Coordinate used to specify the desired node
 * \param[in] ycoord Y Coordinate used to specify the desired node
 * \param[in] zcoord Z Coordinate used to specify the desired node
 * \param[in] tcoord T Coordinate used to specify the desired node
 * \param[out] rank Rank of the desired node
 * \param[out] numnodes Number of Nodes in the partition
 * \return Error indication
 * \retval 0 Success
 * \retval non-zero Error
 */

__INLINE__ int Kernel_Coord2Rank(uint32_t xcoord, uint32_t ycoord, uint32_t zcoord, uint32_t tcoord, uint32_t* rank, uint32_t* numnodes)
{
    int rc = 0; // this syscall returns RC in r3 and does not use errno

  asm __volatile__ ("li 0,%1;"
                    "mr 3,%2;"
                    "mr 4,%3;"
                    "mr 5,%4;"
                    "mr 6,%5;"
                    "mr 7,%6;"
		    "mr 8,%7;"
                     "sc;"
		    "mr %0, 3;"
                    : "=&r" (rc)  // early clobber
                    : "i" (_BGP_SYSCALL_NR_COORD2RANK),
                    "r" (xcoord),
                    "r" (ycoord),
                    "r" (zcoord),
                    "r" (tcoord),
                    "r" (rank),
		    "r" (numnodes)
                    : "r0", "r3", "r4", "r5", "r6", "r7", "r8", "cc", "memory" );

  return( rc );
}

/*! \brief Returns the Job ID
 * \return Contains the control system JobID
 */
__INLINE__ uint32_t Kernel_GetJobID()
{
  int rc = 0; // this syscall returns RC in r3 and does not use errno

  asm __volatile__ ("li 0,%1;"
                     "sc;"
		    "mr %0, 3;"
                    : "=&r" (rc)  // early clobber
                    : "i" (_BGP_SYSCALL_NR_GETJOBID)
                    : "r0", "r3", "cc", "memory" );

  return( rc );
}

/*! \brief Read from a privileged DCR
 * \param[in] dcrid Number of the DCR register
 * \param[out] value Contents of DCR register
 * \return Error indication
 * \retval  0 Success
 * \retval -1 Invalid DCR
 * \note Only selected previleged DCRs will be accessible via this system call.
 */
__INLINE__ uint32_t Kernel_ReadDCR(uint32_t dcrid, uint32_t* value)
{
    int rc = 0; // this syscall returns RC in r3 and does not use errno

  asm __volatile__ ("li 0,%1;"
                    "mr 3,%2;"
                    "mr 4,%3;"
                     "sc;"
		    "mr %0, 3;"
                    : "=&r" (rc)  // early clobber
                    : "i" (_BGP_SYSCALL_NR_READDCR),
		    "r" (dcrid),
		    "r" (value)
                    : "r0", "r3", "r4", "cc", "memory" );

  return( rc );
}

/*! \brief Write to a privileged DCR
 * \param[in] dcrid Number of the DCR register
 * \param[in] value Contents of DCR register
 * \return Error indication
 * \retval  0 Success
 * \retval -1 Invalid DCR
 * \note Only selected previleged DCRs will be accessible via this system call.
 */
__INLINE__ uint32_t Kernel_WriteDCR(uint32_t dcrid, uint32_t value)
{
    int rc = 0; // this syscall returns RC in r3 and does not use errno

  asm __volatile__ ("li 0,%1;"
                    "mr 3,%2;"
                    "mr 4,%3;"
                     "sc;"
		    "mr %0, 3;"
                    : "=&r" (rc)  // early clobber
                    : "i" (_BGP_SYSCALL_NR_WRITEDCR),
		    "r" (dcrid),
		    "r" (value)
                    : "r0", "r3", "r4", "cc", "memory" );

  return( rc );
}

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
 * \internal This function is not intended to be called directly
 * \see DMA_CounterGroupQueryFree()
 * \note The kernel may need to synchronize with other cores performing
 *       allocate or free syscalls.
 *
 */
__INLINE__ uint32_t Kernel_CounterGroupQueryFree(uint32_t type, uint32_t group, uint32_t* num_subgroups, uint32_t* subgroups)
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
                    : "i" (_BGP_SYSCALL_NR_COUNTERGRPQUERYFREE),
		    "r" (type),
		    "r" (group),
			"r" (num_subgroups),
			"r" (subgroups)
                    : "r0", "r3", "r4", "r5", "r6", "cc", "memory" );

  return( rc );
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
 * \internal This function is not intended to be called directly
 * \see DMA_CounterGroupAllocate()
 * \note The kernel may need to synchronize with other cores performing queries
 *       or frees.
 *
 */
__INLINE__ uint32_t Kernel_CounterGroupAllocate(uint32_t type, uint32_t group, uint32_t num_subgroups, uint32_t* subgroups, uint32_t target, uint32_t handler, uint32_t* handler_parm, uint32_t interruptGroup, uint32_t* cg_ptr)
{
  int rc = 0; // this syscall returns RC in r3 and does not use errno

  asm __volatile__ ("li 0,%1;"
                    "mr 3,%2;"
                    "mr 4,%3;"
                    "mr 5,%4;"
					"mr 6,%5;"
					"mr 7,%6;"
                    "mr 8,%7;"
                    "mr 9,%8;"
					"mr 10,%9;"
					"mr 11,%10;"
					"sc;"
		    "mr %0, 3;"
                    : "=&r" (rc)  // early clobber
                    : "i" (_BGP_SYSCALL_NR_COUNTERGRPALLOCATE),
		    "r" (type),
		    "r" (group),
			"r" (num_subgroups),
			"r" (subgroups),
			"r" (target),
			"r" (handler),
			"r" (handler_parm),
			"r" (interruptGroup),
			"r" (cg_ptr)
                    : "r0", "r3", "r4", "r5", "r6", "r7", "r8", "r9", "r10", "r11", "cc", "memory" );

  return( rc );
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
 * \internal This function is not intended to be called directly
 * \see DMA_CounterGroupFree()
 * \note The kernel may need to synchronize with other cores performing allocates
 *       or queries.
 */
__INLINE__ uint32_t Kernel_CounterGroupFree(uint32_t group, uint32_t num_subgroups, uint32_t* subgroups, uint32_t* cg_ptr)
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
                    : "i" (_BGP_SYSCALL_NR_COUNTERGRPFREE),
		    "r" (group),
			"r" (num_subgroups),
			"r" (subgroups),
			"r" (cg_ptr)
                    : "r0", "r3", "r4", "r5", "r6", "cc", "memory" );

  return( rc );
}


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
 * \param[out]  fifo_ids       Pointer to an array of num_fifos ints where
 *                             the list of free fifos is returned.
 *                             Each int is the fifo number
 *                             (0 to DMA_NUM_INJ_FIFOS_PER_GROUP-1).
 *                             The caller must provide space for
 *                             DMA_NUM_INJ_FIFOS_PER_GROUP ints,
 *                             in case the entire fifo group is free.
 *
 * \retval  0  Successful.  num_fifos and fifo_ids array set as described.
 * \retval  -1 Unsuccessful.  errno gives the reason.
 * \internal This function is not intended to be called directly
 * \see DMA_InjFifoGroupQueryFree()
 */
__INLINE__ uint32_t Kernel_InjFifoGroupQueryFree(uint32_t group, uint32_t* num_fifos, uint32_t* fifo_ids)
{
  int rc = 0; // this syscall returns RC in r3 and does not use errno

  asm __volatile__ ("li 0,%1;"
                    "mr 3,%2;"
                    "mr 4,%3;"
                    "mr 5,%4;"
                     "sc;"
		    "mr %0, 3;"
                    : "=&r" (rc)  // early clobber
                    : "i" (_BGP_SYSCALL_NR_INJFIFOGRPQUERYFREE),
		    "r" (group),
			"r" (num_fifos),
			"r" (fifo_ids)
                    : "r0", "r3", "r4", "r5", "cc", "memory" );

  return( rc );
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
 *                           Each int is the fifo number (0 to num_fifos-1).
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
 * \param[in]   ts_inj_maps  Pointer to an array of num_fifos short ints where
 *                           the torus fifos that can be injected are specified
 *                           for each fifo.  Each short int specifies which of
 *                           the 8 torus injection fifos can be injected when a
 *                           descriptor is injected into the DMA injection fifo.
 *                           Must be non-zero when the corresponding "locals"
 *                           is 0.
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
 * \internal This function is not intended to be called directly
 * \see DMA_InjFifoGroupAllocate()
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
__INLINE__ uint32_t Kernel_InjFifoGroupAllocate(uint32_t group, uint32_t num_fifos, uint32_t* fifo_ids, uint16_t* priorities, uint16_t* locals, uint8_t* ts_inj_maps, uint32_t* fg_ptr)
{
  int rc = 0; // this syscall returns RC in r3 and does not use errno

  asm __volatile__ ("li 0,%1;"
                    "mr 3,%2;"
                    "mr 4,%3;"
                    "mr 5,%4;"
					"mr 6,%5;"
					"mr 7,%6;"
                    "mr 8,%7;"
                    "mr 9,%8;"
		    "sc;"
		    "mr %0, 3;"
                    : "=&r" (rc)  // early clobber
                    : "i" (_BGP_SYSCALL_NR_INJFIFOGRPALLOCATE),
		    "r" (group),
			"r" (num_fifos),
			"r" (fifo_ids),
			"r" (priorities),
			"r" (locals),
			"r" (ts_inj_maps),
			"r" (fg_ptr)
                    : "r0", "r3", "r4", "r5", "r6", "r7", "r8", "r9", "cc", "memory" );

  return( rc );
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
 * \internal This function is not intended to be called directly
 * \see DMA_InjFifoGroupFree()
 * \note  This is a fatal error if any of the fifos are non empty and activated
 *
 */
__INLINE__ uint32_t Kernel_InjFifoGroupFree(uint32_t group, uint32_t num_fifos, uint32_t* fifo_ids, uint32_t* fg_ptr)
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
                    : "i" (_BGP_SYSCALL_NR_INJFIFOGRPFREE),
		    "r" (group),
			"r" (num_fifos),
			"r" (fifo_ids),
			"r" (fg_ptr)
                    : "r0", "r3", "r4", "r5", "r6", "cc", "memory" );

  return( rc );
}

/*!
 * \brief DMA InjFifo Initialization By Id
 *
 * - For an allocated injection DMA fifo, initialize its start, head, tail, and
 *   end.
 * - Compute fifo size and free space.
 * - Initialize wrap count.
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
__INLINE__ uint32_t Kernel_InjFifoInitById(uint32_t* fg_ptr,
				     int  fifo_id,
				     uint32_t* va_start,
				     uint32_t* va_head,
				     uint32_t* va_end)
{
    	int rc = 0; // this syscall returns RC in r3 and does not use errno

	asm __volatile__ ("li 0,%1;"
                    "mr 3,%2;"
                    "mr 4,%3;"
                    "mr 5,%4;"
					"mr 6,%5;"
					"mr 7,%6;"
					"sc;"
		    "mr %0, 3;"
                    : "=&r" (rc)  // early clobber
                    : "i" (_BGP_SYSCALL_NR_INJFIFOINITID),
		    "r" (fg_ptr),
			"r" (fifo_id),
			"r" (va_start),
			"r" (va_head),
			"r" (va_end)
                    : "r0", "r3", "r4", "r5", "r6", "r7", "cc", "memory" );

  return( rc );
}


/*!
 * \brief Set DMA Reception Fifo Map
 *
 * This function is a wrapper around a system call that
 * - Sets DCRs establishing the map between the hardware torus fifos and the
 *   DMA reception fifos that are to receive the packets from those hardware
 *   torus fifos.
 * - Sets DCRs establishing the DMA reception fifos that are to receive
 *   local transfer packets.
 * - Sets the DCRs establishing the type (0 or 1) of each reception fifo.
 * - Sets the DCRs establishing the threshold for type 0 and 1 reception fifos.
 * - Leaves all of the fifos that are used in a "disabled" state.
 *   DMA_RecFifoInitById() initializes and enables the fifos.
 *
 * \param[in]  rec_map  Reception Fifo Map structure, defining the mapping.
 *
 * \retval  0            Successful
 * \retval  error_value  An error value defined in the _BGP_RAS_DMA_ErrCodes
 *                       enum located in bgp/arch/include/common/bgp_ras.h
 *
 * \internal This is an internal syscall
 * \see DMA_RecFifoSetMap
 * \note  This function should be called once per job, after DMA_ResetRelease().
 *        It may be called by any core, but once a core has called it, other
 *        calls by that same core or any other core will fail.
 *
 * \note  During job init, the kernel sets up the DCR clear masks for each
 *        reception fifo group (DCRs 0xD68 - 0xD6C) such that a write to clear
 *        a fifo in group g only clears group g.
 *
 */
__INLINE__ int Kernel_RecFifoSetMap(uint32_t* rec_map)
{
  int rc = 0; // this syscall returns RC in r3 and does not use errno

  asm __volatile__ ("li 0,%1;"
                    "mr 3,%2;"
					"sc;"
		    "mr %0, 3;"
                    : "=&r" (rc)  // early clobber
                    : "i" (_BGP_SYSCALL_NR_RECFIFOSETMAP),
		    "r" (rec_map)
                    : "r0", "r3", "cc", "memory" );
  return( rc );
}

/*!
 * \brief Get DMA Reception Fifo Map
 *
 * This function is a wrapper around a system call that returns a DMA
 * reception fifo map structure, filled in according to the DCRs.
 *
 * \param[in,out]  rec_map  A pointer to a Reception Fifo Map structure
 *                          that will be filled-in upon return.
 *
 * \retval  0            Successful
 * \retval  error_value  An error value defined in the _BGP_RAS_DMA_ErrCodes
 *                       enum located in bgp/arch/include/common/bgp_ras.h
 *
 */
__INLINE__ int Kernel_RecFifoGetMap(uint32_t* rec_map)
{
  int rc = 0; // this syscall returns RC in r3 and does not use errno

  asm __volatile__ ("li 0,%1;"
                    "mr 3,%2;"
					"sc;"
		    "mr %0, 3;"
                    : "=&r" (rc)  // early clobber
                    : "i" (_BGP_SYSCALL_NR_RECFIFOGETMAP),
		    "r" (rec_map)
                    : "r0", "r3", "cc", "memory" );
  return( rc );
}

/*!
 * \brief Get DMA Reception Fifo Group
 *
 * This is a wrapper around a System Call. This function returns THE
 * one-and-only pointer to the fifo group structure, with the entries all
 * filled in from info in the DCRs.  If called multiple times with the same
 * group, it will always return the same pointer, and the system call will
 * not be invoked again.
 *
 * It must be called AFTER DMA_RecFifoSetMap().
 *
 * By convention, the same "target" is used for normal and header fifo
 * interrupts (could be changed).  In addition, by convention, interrupts for
 * fifos in group g come out of the DMA as non-fatal irq bit 28+g,
 * ie, only fifos in group g can cause the "type g" threshold interrupts.
 *
 * \param[in]  grp      The group number (0 through DMA_NUM_REC_FIFO_GROUPS).
 * \param[in]  target   The core that will receive the interrupt when a
 *                      fifo in this group reaches its threshold
 *                      (0 to DMA_NUM_REC_FIFO_GROUPS-1).
 *                      Ignored on subsequent call with the same group.
 * \param[in]  normal_handler  A pointer to the function to receive control in
 *                             the I/O thread to handle the interrupt when a 
 *                             normal fifo in this group reaches its threshold.
 *                             This function must be coded to take 4 uint32_t
 *                             parameters:
 *                             - A pointer to storage specific to this 
 *                               handler.  This is the normal_handler_parm
 *                               specified on this function call.
 *                             - 3 uint32_t parameters that are not used.
 *                             If normal_handler is NULL, threshold interrupts
 *                             are not delivered for normal fifos in this group.
 *                             Ignored on subsequent call with the same group.
 * \param[in]  normal_handler_parm   A pointer to storage that should be passed 
 *                                   to the normal interrupt handling function 
 *                                   (see normal_handler parameter).
 *                                   Ignored on subsequent call with the same 
 *                                   group.
 * \param[in]  header_handler  ** This parameter is deprecated.  Specify NULL.**
 *                             A pointer to the function to receive control in
 *                             the I/O thread to handle the interrupt when a 
 *                             header fifo in this group reaches its threshold.
 *                             This function must be coded to take 2 parameters:
 *                               void* A pointer to storage specific to this 
 *                                     handler.  This is the header_handler_parm
 *                                     specified on this function call.
 *                               int   The global fifo ID of the fifo that hit
 *                                     its threshold (0 through 
 *                                     NUM_DMA_REC_FIFOS-1).
 *                             If header_handler is NULL, threshold interrupts
 *                             are not delivered for header fifos in this group.
 *                             Ignored on subsequent call with the same group.
 * \param[in]  header_handler_parm   ** This parameter is deprecated.  Specify
 *                                      NULL. **
 *                                   A pointer to storage that should be passed 
 *                                   to the header interrupt handling function 
 *                                   (see header_handler parameter).
 *                                   Ignored on subsequent call with the same 
 *                                   group.
 * \param[in]  interruptGroup  A InterruptGroup_t that identifies the
 *                             group of interrupts that the fifos in this group
 *                             will become part of.
 *                             Ignored on subsequent call with the same group.
 *
 * \return  RecFifoGroupStruct  Pointer to a DMA Reception Fifo Group structure
 *                              that reflects the fifos that are being used in
 *                              this group.  This same structure is shared by
 *                              all users of this reception fifo group.
 *                              NULL is returned if an error occurs.
 *
 * \note  The following comments from Phil about the internals of the syscall:
 *   - error checks
 *     - 0 <= group_id < 4
 *     - the start of the fifo group is a valid virtual address (tlb mapped)?
 *   - disable the rDMA
 *   - call _BGP_rDMA_Fifo_Get_Map to get the DCR mapping information
 *   - loop through the map to determine how many and which fifos in this group
 *     are used, including headers
 *   - filling in the addresses of used fifos
 *     - In particular, any pointer to any fifo in the group that is not used
 *       will have a null pointer
 *   - furthermore,
 *     - write starting values to all used fifos
 *     - make sure all interrupts are cleared
 *     - enable rDMA
 *
 */
__INLINE__ int Kernel_RecFifoGetFifoGroup(
			uint32_t*		      	  fifogroup,
			int                               grp,
			int                               target,
			Kernel_CommThreadHandler          normal_handler,
			void                             *normal_handler_parm,
			Kernel_CommThreadHandler          header_handler,
			void                             *header_handler_parm,
			Kernel_InterruptGroup_t           interruptGroup
		       )
{
	int rc = 0; // this syscall returns RC in r3 and does not use errno

	asm __volatile__ ("li 0,%1;"
                    "mr 3,%2;"
                    "mr 4,%3;"
                    "mr 5,%4;"
					"mr 6,%5;"
					"mr 7,%6;"
					"mr 8,%7;"
					"mr 9,%8;"
					"mr 10,%9;"
					"sc;"
		    "mr %0, 3;"
                    : "=&r" (rc)  // early clobber
                    : "i" (_BGP_SYSCALL_NR_RECGETFIFOGROUP),
			  "r" (fifogroup),
			  "r" (grp),
			"r" (target),
			"r" (normal_handler),
			"r" (normal_handler_parm),
			"r" (header_handler),
			"r" (header_handler_parm),
			"r" (interruptGroup)
                    : "r0", "r3", "r4", "r5", "r6", "r7", "r8", "r9", "r10", "cc", "memory" );

  return( rc );
}

/*!
 * \brief DMA RecFifo Initialization By Id
 *
 * - For a DMA reception fifo, initialize its start, head, tail, and end.
 * - Compute fifo size and free space.
 *
 * \param[in]  fg_ptr    Pointer to fifo group structure.
 * \param[in]  fifo_id   Id of the fifo to be initialized
 *                       (0 to DMA_NUM_REC_FIFOS_PER_GROUP-1).
 * \param[in]  va_start  Virtual address of the start of the fifo.
 * \param[in]  va_head   Virtual address of the head of the fifo (typically
 *                       equal to va_start).
 * \param[in]  va_end    Virtual address of the end of the fifo.
 *
 * \retval   0  Successful.
 * \retval  -1  Unsuccessful.  Error checks include
 *              - va_start <  va_end
 *              - va_start <= va_head < va_end
 *              - va_start and va_end are 32-byte aligned
 *              - fifo_size >= DMA_MIN_REC_FIFO_SIZE_IN_BYTES
 *
 */
__INLINE__ int Kernel_RecFifoInitById(
				   uint32_t*		  fg_ptr,
				   int                fifo_id,
				   void               *va_start,
				   void               *va_head,
				   void               *va_end
				  )
{
	int rc = 0; // this syscall returns RC in r3 and does not use errno

	asm __volatile__ ("li 0,%1;"
                    "mr 3,%2;"
                    "mr 4,%3;"
                    "mr 5,%4;"
					"mr 6,%5;"
					"mr 7,%6;"
					"sc;"
		    "mr %0, 3;"
                    : "=&r" (rc)  // early clobber
                    : "i" (_BGP_SYSCALL_NR_RECFIFOINITID),
		    "r" (fg_ptr),
			"r" (fifo_id),
			"r" (va_start),
			"r" (va_head),
			"r" (va_end)
                    : "r0", "r3", "r4", "r5", "r6", "r7", "cc", "memory" );

  return( rc );
}

 /*! 
  * \brief Injects a binary (RAW) RAS message to the control system
  * 
  * Ships a RAS message of the given facility, unit, errcode, and packed data to the control system.  No checking is done on the
  * correctness of the data.  Can be used to simulate a RAS message for testing purposes.
  * 
  * \param[in] facility High level component detecting the condition. (e.g., _bgp_fac_kernel, _bgp_fac_application, _bgp_fac_diags)
  * \param[in] unit Unit generating the RAS event.  (e.g., _bgp_unit_ppc450, _bgp_unit_snoop)
  * \param[in] err_code Error code for RAS event (e.g., _bgp_err_ppc450_l1d_dpe0)
  * \param[in] numwords Number of 32-bit integers in the packed binary array
  * \param[in] array Pointer to the array of packed binary data.
  * 
  * Restriction.  There is currently a limit of eight 32-bit words of packed binary data.
  * 
  * \internal This function is intended for testing purposes only.  It should not be used in a production system as it could introduce false RAS messages.
 */

__INLINE__ int Kernel_InjectRAWRAS(
				   _BGP_Facility facility,
				   _BGP_RAS_Units unit,
				   uint16_t err_code,
				   int numwords,
				   const uint32_t* array)
{
  int rc = 0; // this syscall returns RC in r3 and does not use errno                                                                           

  asm __volatile__ ("li 0,%1;"
                    "mr 3,%2;"
                    "mr 4,%3;"
                    "mr 5,%4;"
		    "mr 6,%5;"
		    "mr 7,%6;"
		    "sc;"
                    "mr %0, 3;"
                    : "=&r" (rc)  // early clobber                                                                                                    
                    : "i" (_BGP_SYSCALL_NR_RAWRASINJECT),
                    "r" (facility),
		    "r" (unit),
		    "r" (err_code),
		    "r" (numwords),
		    "r" (array)
                    : "r0", "r3", "r4", "r5", "r6", "r7", "cc", "memory" );

  return( rc );
}

 /*! 
  * \brief Injects a ASCII (Textual) RAS message to the control system
  * 
  * Ships a RAS message of the given facility, unit, errcode, and an ASCII string to the control system.  No checking is done on the
  * correctness of the facility or unit.  Can be used to simulate a RAS message for testing purposes.
  * 
  * \param[in] facility High level component detecting the condition. (e.g., _bgp_fac_kernel, _bgp_fac_application, _bgp_fac_diags)
  * \param[in] unit Unit generating the RAS event.  (e.g., _bgp_unit_ppc450, _bgp_unit_snoop)
  * \param[in] err_code Error code for RAS event (e.g., _bgp_err_ppc450_l1d_dpe0)
  * \param[in] text Pointer to a string of null-terminated ASCII characters
  * 
  * \internal This function is intended for testing purposes only.  It should not be used in a production system as it could introduce false RAS messages.
 */
__INLINE__ int Kernel_InjectASCIIRAS(
				   _BGP_Facility facility,
				   _BGP_RAS_Units unit,
				   uint16_t err_code,
				   const uint8_t* text)
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
                    : "i" (_BGP_SYSCALL_NR_ASCIIRASINJECT),
                    "r" (facility),
		    "r" (unit),
		    "r" (err_code),
		    "r" (text)
                    : "r0", "r3", "r4", "r5", "r6", "cc", "memory" );

  return( rc );
}



/*!
 * \brief Enables/Disables the counter overflow/underflow interrupts
 *
 * This function is a wrapper around a system call that can enable or disable the 4 counter overflow/underflow interrupts
 *
 * \param[in]  enable/disable boolean
 *
 * \retval  0            Successful
 * \retval  error_value  An error value defined in the _BGP_RAS_DMA_ErrCodes
 *                       enum located in bgp/arch/include/common/bgp_ras.h
 *
 */
__INLINE__ int Kernel_ChgCounterInterruptEnables(uint32_t enable)
{
  int rc = 0; // this syscall returns RC in r3 and does not use errno

  asm __volatile__ ("li 0,%1;"
                    "mr 3,%2;"
                    "sc;"
                    "mr %0, 3;"
                    : "=&r" (rc)  // early clobber
                    : "i" (_BGP_SYSCALL_NR_CHGDMACTRINTERRUPT),
                    "r" (enable)
                    : "r0", "r3", "cc", "memory" );
  return( rc );
}


/*!
 * \brief Clears the Full Reception FIFO (DD1 workaround)
 *
 * This function exists to reset the DMA reception fifos - it is a workaround for DD1 only.  It should not be needed in DD2.
 *
 * \retval  0            Successful
 * \retval  error_value  An error value defined in the _BGP_RAS_DMA_ErrCodes
 *                       enum located in bgp/arch/include/common/bgp_ras.h
 *
 */
__INLINE__ int Kernel_ClearFullReceptionFifo()
{
  int rc = 0; // this syscall returns RC in r3 and does not use errno

#ifdef __CNK__

  rc =  _bgp_DMA_ClearFullReceptionFifo();


#else
  asm __volatile__ ("li 0,%1;"
                    "sc;"
                    "mr %0, 3;"
                    : "=&r" (rc)  // early clobber
                    : "i" (_BGP_SYSCALL_NR_CLEARFULLRECPFIFO)
                    : "r0", "r3", "cc", "memory" );
#endif
  return( rc );
}

#include <spi/lockbox_interface.h>

#if ((!defined(__CNK__)) && (!defined(__BL__)))
/*! \brief Creates a pthread with a commthread attribute
 *
 * \note CNK restriction:  1 CommThread per core is allowed
 * \note In Dual or VNM, each process must allocate its own commthreads
 * \note CommThreads are pinned per core.  (e.g., in SMP mode, this SPI must be called 4 times to create enough CommThreads for each processor)
 * \warning non-portable pthread API
 * \param[in] thread pthread_t structure
 * \param[in] attr   pthread_attr_t structure
 * \param[in] start_routine function pointer of the thread's main()
 * \param[in] arg    1st argument to the pthread
 * \return Error condition from pthread_create
 * \retval 0 success
 * \retval -1 error, check errno
 */
__INLINE__ int pthread_create_CommThread_np( pthread_t *thread,
                                          pthread_attr_t *attr,
                                          void *(*start_routine)(void *),
                                          void *arg )
{
  uint32_t usprg0 = _bgp_mfspr( SPRN_USPRG0 ); // save orig usprg0

  _bgp_mtspr( SPRN_USPRG0, _BGP_COMMTHREAD_MAGIC );

  int rc = pthread_create( thread, attr, start_routine, arg );
  _bgp_mtspr( SPRN_USPRG0, usprg0 );  // restore orig usprg0

  return( rc );
}
#endif

/*! \brief Causes a commthread to disappear from the runqueue
 *
 *  \note Kernel does not guarantee that the instruction pointer, stack pointer, and register state are preserved across a poof.  
 *  \note TLS data is preserved across a poof
 *  \note This SPI is only executable on a comm. thread.
 *  \warning non-portable pthread API
 *  \return error indication
 *  \retval success Does not return.  Thread has "poofed"
 *  \retval -1 Calling thread is not a CommThread, so cannot poof
 */
__INLINE__ int pthread_poof_np( void )
{
    int rc = 0; // this syscall returns RC in r3 and does not use errno

  asm __volatile__ ("li 0,%1;"
                     "sc;"
		    "mr %0, 3;"
		    : "=&r" (rc)  // early clobber
		    : "i" (_BGP_SYSCALL_NR_PTHREAD_POOF)
		    : "r0", "r3", "cc", "memory" );

  return( rc );
}



/*! \defgroup COMMTHRD_OPCODES CommThread Opcodes
 *  \{
 * \note Only 1 interrupt route can be specified per opcode
 * \note CallFunc, DisableIntOnEntry, EnableIntOnPoof can be specified in any combination
 * \note Current support requires that DisableIntOnEntry and EnableIntOnPoof be specified
 */
#define COMMTHRD_OPCODE_DISABLE            0x00 //!< Interrupt route - Not routed / interrupt disabled
#define COMMTHRD_OPCODE_CORE0              0x01 //!< Interrupt route - Dispatched on core0
#define COMMTHRD_OPCODE_CORE1              0x02 //!< Interrupt route - Dispatched on core1
#define COMMTHRD_OPCODE_CORE2              0x03 //!< Interrupt route - Dispatched on core2
#define COMMTHRD_OPCODE_CORE3              0x04 //!< Interrupt route - Dispatched on core3
#define COMMTHRD_OPCODE_BCAST              0x05 //!< Interrupt route - Dispatched on all cores
#define COMMTHRD_OPCODE_ROUTEMASK          0x0F //!< Interrupt route mask
#define COMMTHRD_OPCODE_CALLFUNC           0x10 //!< The provided function will be called on the comm. thread
#define COMMTHRD_OPCODE_DISABLEINTONENTRY  0x20 //!< Interrupts using cntrid will be disabled when comm. thread is invoked
#define COMMTHRD_OPCODE_ENABLEINTONPOOF    0x40 //!< Interrupts using cntrid will be enabled when comm. thread poofs
/*!
 * \}
 */

/*! \brief Generates an InterruptID value
 * \param[in] group group of the interrupt.  range 0-9.
 * \param[in] irq_in_group irq within the group.  range 0-31.
 * \return Composite value able to be passed to Kernel_SetCommThreadConfig
 * \see Kernel_SetCommThreadConfig
 */
#define Kernel_MkInterruptID(group, irq_in_group) ((group<<5)|(irq_in_group&0x1f))

/*! 
 * \brief Sets kernel data structures needed to dispatch a communications thread
 * 
 * Each interrupt on BGP can be used to launch a communications thread.  Since access to the
 * interrupt controller is privileged, the function exposes some interrupt control to the 
 * user application.
 * \pre Counter must have been allocated via the LockBox_AllocateCounter() routine.
 * \pre It is recommended that Kernel_DisableInteruptClass() be called twice on the counter
 *      to ensure that the interrupt is disabled until all interrupts for the counter
 *      have been configured.
 * \pre All 
 * \post After the last call to Kernel_SetCommThreadConfig for the counter, invoke
 *       Kernel_EnableInterruptClass() and Kernel_HardEnableInterruptClass() on 
 *       that counter to enable the interrupts for that class.
 * \see LockBox_AllocateCounter
 * \see Kernel_DisableInterruptClass
 * \see Kernel_EnableInterruptClass
 * \see Kernel_HardEnableInterruptClass
 * \note An interrupt can only belong to 1 interrupt class (a.k.a., lockbox counter)
 * \note The effects of this function span the entire node regardless of SMP, Dual, or VNM settings. 
 * \note Kernel may prevent changing interrupt settings for certain InterruptID values.  
 * \note If an interrupt fires on a core without a comm. thread, results are not guaranteed.  
 * \return Completion status of the command.
 * \retval 0 no error occurred
 * \retval EINVAL invalid parameter
 * \param[in] InterruptID  Identifies a unique interrupt line.  \see Kernel_MkInterruptID
 * \param[in] opcode       Specifies what operation to perform when the interrupt occurs. Valid \ref COMMTHRD_OPCODES
 * \param[in] cntrid       ID of the lockbox counter used for interrupt enable/disable control
 * \param[in] funcptr      Function pointer that will be invoked when the interrupt fires
 * \param[in] arg1         1st argument to the funcptr when the interrupt fires
 * \param[in] arg2         2nd argument to the funcptr when the interrupt fires
 * \param[in] arg3         3rd argument to the funcptr when the interrupt fires
 * 
 */
__INLINE__ int Kernel_SetCommThreadConfig(int InterruptID, int opcode, LockBox_Counter_t cntrid, 
					  Kernel_CommThreadHandler funcptr, 
					  uint32_t arg1, uint32_t arg2, uint32_t arg3, uint32_t arg4)
{
   int rc = 0;
   asm __volatile__ ("li 0,%1;"
		     "mr 3, %2;"
		     "mr 4, %3;"
		     "mr 5, %4;"
		     "mr 6, %5;"
		     "mr 7, %6;"
		     "mr 8, %7;"
		     "mr 9, %8;"
		     "mr 10, %9;"
		     "sc;"
		     "mr %0, 3;"
		     : "=&r" (rc) // early clobber
		     : "i" (_BGP_SYSCALL_NR_SETCOMMTHREADCONFIG),
		     "r" (InterruptID),
		     "r" (opcode),
		     "r" (cntrid),
		     "r" (funcptr),
		     "r" (arg1),
		     "r" (arg2),
		     "r" (arg3),
		     "r" (arg4)
		    : "r0", "r3", "r4", "r5", "r6", "r7", "r8", "r9", "r10", "cc", "memory" );
  return rc;
}

/*! 
 * \brief Returns the kernel data structures that were specified to dispatch communication thread
 * 
 * Each interrupt on BGP can be used to launch a communications thread.  Since access to the
 * interrupt controller is privileged, the function exposes some interrupt control to the 
 * user application.
 *
 * \param[in] InterruptID  Identifies a unique interrupt line.  
 * \param[out] opcode    Storage for opcode value.  Specifies which core receives the interrupt.  It also controls whether the interrupt disables a class of interrupts.  Valid \ref COMMTHRD_OPCODES
 * \param[out] cntrid       Storage for ID of the lockbox counter used for interrupt enable/disable control
 * \param[out] funcptr      Storage for Function pointer that will be invoked when the interrupt fires
 * \param[out] arg1         Storage for 1st argument to the funcptr when the interrupt fires
 * \param[out] arg2         Storage for 2nd argument to the funcptr when the interrupt fires
 * \param[out] arg3         Storage for 3rd argument to the funcptr when the interrupt fires
 * \return Completion status of the command.
 * \retval 0 no error occurred
 * \retval EINVAL invalid parameter
 * 
 */
__INLINE__ int Kernel_GetCommThreadConfig(int InterruptID, int* opcode, LockBox_Counter_t* cntrid, 
					  Kernel_CommThreadHandler* funcptr, 
					  uint32_t* arg1, uint32_t* arg2, uint32_t* arg3, uint32_t* arg4)
{
  int rc = 0;
  asm __volatile__ ("li 0,%1;"
                    "mr 3, %2;"
                    "mr 4, %3;"
                    "mr 5, %4;"
                    "mr 6, %5;"
                    "mr 7, %6;"
                    "mr 8, %7;"
                    "mr 9, %8;"
		    "mr 10, %9;"
                    "sc;"
                    "mr %0, 3;"
                    : "=&r" (rc) // early clobber
                    : "i" (_BGP_SYSCALL_NR_GETCOMMTHREADCONFIG),
		    "r" (InterruptID),
		    "r" (opcode),
		    "r" (cntrid),
		    "r" (funcptr),
		    "r" (arg1),
		    "r" (arg2),
		    "r" (arg3),
		    "r" (arg4)
                    : "r0", "r3", "r4", "r5", "r6", "r7", "r8", "r9", "r10", "cc", "memory" );
  return rc;
}

/*! \brief Flush interrupt enable/disable state
 *
 * For each interrupt that has a lockbox counter associated with it, this SPI will 
 * update the interrupt controller to match the state specified by the lockbox counter.
 * \note The effects of this function span the entire node regardless of SMP, Dual, or VNM settings. 
 * \note Kernel is responsible for updating the interrupt controller to match all lockbox counters
 * 
 * \return Completion status of the command.
 * \retval 0 no error occurred
 */
__INLINE__ int Kernel_FlushInterruptState()
{
   int rc;
   asm __volatile__ ("li 0,%1;"
		     "sc;"
		     "mr %0, 3;"
		     : "=&r" (rc)  // early clobber
		     : "i" (_BGP_SYSCALL_NR_FLUSHINTSTATE)
		     : "r0", "r3", "cc", "memory" );
   return rc;
}

/*! \brief Indicates that the kernel should disable the interrupt
 *
 * Updates the interrupt class's lockbox to indicate that the kernel should disable the interrupt.
 * Kernel will disable the interrupt at its leisure, but it should ensure that no communications thread
 * is invoked for that interrupt class.  
 *
 * The lockbox values have the following meanings:
 * 0: Interrupts for this classid are enabled
 * 1: Interrupts for this classid are logically disabled.  
 *    If an interrupt occurs, the kernel will hard-disable them and ignore the interrupt.
 * 2: Interrupts for this classid are hard-disabled.  The interrupt will not disturb the core.
 *
 * \note The effects of this function span the entire node regardless of SMP, Dual, or VNM settings. 
 * \note Do not disable an already disabled interrupt class.
 * \note A disabled interrupt class is disabled for all 4 cores, regardless of mode.
 * \param[in] classid An allocated lockbox that is being used to control a set of interrupt enable/disable lines
 * 
 */
__INLINE__ uint32_t Kernel_DisableInterruptClass(LockBox_Counter_t classid)
{
  return ( LockBox_FetchAndInc(classid) );
}

/*! \brief Indicates that the kernel should enable the interrupt
 * 
 * Updates the interrupt class's lockbox to indicate that the kernel should leave this interrupt enabled.
 * This does not hard-enable the interrupts for this classid (see Kernel_HardEnableInterruptClass).
 *
 * The lockbox values have the following meanings:
 * 0: Interrupts for this classid are enabled
 * 1: Interrupts for this classid are logically disabled.  
 *    If an interrupt occurs, the kernel will hard-disable them and ignore the interrupt.
 * 2: Interrupts for this classid are hard-disabled.  The interrupt will not disturb the core.
 *
 * \note The effects of this function span the entire node regardless of SMP, Dual, or VNM settings. 
 * \note The kernel is responsible for incrementing the lockbox counter when the interrupt is hard-disabled. 
 * \note There is potential race condition that must be avoided in the kernel.  The kernel will need to Query the lockbox when an interrupt occurs, and if it is non-zero, then increment it (another core could enable the interrupt class between those 2 events).  One solution is to always FetchAndInc, but that may lead to an extranous (but rare) FlushInterruptState() call, followed by a FetchAndDec if zero.  There are fancier solutions as well.
 * \param[in] classid An allocated lockbox that is being used to control a set of interrupt enable/disable lines
 * 
 */
__INLINE__ uint32_t Kernel_EnableInterruptClass(LockBox_Counter_t classid)
{
  return ( LockBox_FetchAndDec(classid) );
}

/*! \brief Indicates that the kernel should hard enable the interrupt
 * 
 * Updates the interrupt class's lockbox to indicate that the kernel has hard-enabled this interrupt.
 * If the kernel has actually disabled the interrupt, this SPI will enable the interrupt by using the
 * Kernel_FlushInterruptState() SPI.  
 *
 * The lockbox values have the following meanings:
 * 0: Interrupts for this classid are enabled
 * 1: Interrupts for this classid are logically disabled.  
 *    If an interrupt occurs, the kernel will hard-disable them and ignore the interrupt.
 * 2: Interrupts for this classid are hard-disabled.  The interrupt will not disturb the core.
 *
 * \note The effects of this function span the entire node regardless of SMP, Dual, or VNM settings. 
 * \note The kernel is responsible for incrementing the lockbox counter when the interrupt is disabled. 
 * \note There is potential race condition that must be avoided in the kernel.  The kernel will need to Query the lockbox when an interrupt occurs, and if it is non-zero, then increment it (another core could enable the interrupt class between those 2 events).  One solution is to always FetchAndInc, but that may lead to an extranous (but rare) FlushInterruptState() call, followed by a FetchAndDec if zero.  There are fancier solutions as well.
 * \param[in] classid An allocated lockbox that is being used to control a set of interrupt enable/disable lines
 * 
 */
__INLINE__ void Kernel_HardEnableInterruptClass(LockBox_Counter_t classid)
{
  LockBox_FetchAndDec(classid);
  Kernel_FlushInterruptState();
}

/*! \brief Delivers an interrupt to the cores specified in the mask
 * \param[in] coremask Bitmask describing which processor cores will receive the interrupt.  Processor 0 is the least significant bit (1<<0 in C parlance).  Processor 3 is 1<<3.  Any combination of processors can be interrupted.  
 * \note It is possible to interrupt yourself.  
 */
__INLINE__ int Kernel_DeliverCommSignal(uint32_t ipiset, uint32_t coremask)
{
   int rc = 0;
   asm __volatile__ ("li 0,%1;"
		     "mr 3, %2;"
		     "mr 4, %3;"
		     "sc;"
                     "mr %0, 3;"
		     : "=&r" (rc)  // early clobber
		     : "i" (_BGP_SYSCALL_NR_DELIVERCOMMSIGNAL),
		     "r" (ipiset),
		     "r" (coremask)
		     : "r0", "r3", "r4", "cc", "memory" );
  return rc;
}

/*!
 * \brief Suspends/Resumes a core
 *
 * \param[in]  target core ID
 * \param[in]  suspend Boolean.  TRUE if core is to be suspended
 *
 * \retval  0            Successful
 * \retval  error_value  An error value defined in the _BGP_RAS_DMA_ErrCodes
 *                       enum located in bgp/arch/include/common/bgp_ras.h
 * \note In a threaded application, use care to avoid suspending a thread containing a lock needed by the active thread.  (e.g., if the other core is performing a printf, it may have the glibc io subsystem locked with a mutex.  If that happens, the main thread may deadlock if it also happens to call printf)
 *
 */
__INLINE__ int Kernel_ChangeCoreEnables(uint32_t target_core, uint32_t suspend)
{
  int rc = 0; // this syscall returns RC in r3 and does not use errno

  asm __volatile__ ("li 0,%1;"
                    "mr 3,%2;"
                    "mr 4,%3;"
                    "sc;"
                    "mr %0, 3;"
                    : "=&r" (rc)  // early clobber
                    : "i" (_BGP_SYSCALL_NR_CHGCOREENABLES),
                    "r" (target_core),
		    "r" (suspend)
                    : "r0", "r3", "cc", "memory" );
  return( rc );
}

/*! \brief Persistent Shared Memory interface to application. Currently, simlpy a wrapper to open(2),
 *         with a prefix of /dev/persist
 */
__INLINE__ int persist_open( char *name, int oflag, mode_t mode )
{
   char pathName[PATH_MAX];
   strcpy(pathName, "/dev/persist/");
   strncat(pathName, name, PATH_MAX - strlen("/dev/persist/") - 1);
   return open(pathName, oflag, mode);
}

#if SPI_DEPRECATED

//! \see Kernel_PhysicalProcessorID
#define BGP_PhysicalProcessorID  Kernel_PhysicalProcessorID

//! \see Kernel_Virtual2Physical
#define _bgp_Virtual2Physical    Kernel_Virtual2Physical

//! \see Kernel_GetPersonality
#define rts_get_personality(p,s)    Kernel_GetPersonality(p,s)

//! \see Kernel_PhysicalProcessorID
#define rts_get_processor_id()      Kernel_PhysicalProcessorID()

//! \see Kernel_GetProcessorVersion
#define rts_get_processor_version() Kernel_GetProcessorVersion()
#endif

__END_DECLS




#endif // Add nothing below this line
