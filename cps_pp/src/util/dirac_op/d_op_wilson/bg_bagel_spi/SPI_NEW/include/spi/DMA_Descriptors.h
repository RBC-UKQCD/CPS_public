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

#ifndef _DMA_DESCRIPTORS_H_ /* Prevent multiple inclusion */
#define _DMA_DESCRIPTORS_H_

/*!
 * \file spi/DMA_Descriptors.h
 *
 * \brief DMA SPI Descriptor Definitions and Inline Functions
 *
 * This header file contains the definition of the DMA_InjDescriptor_t, which is
 * put into the tail of an injection fifo to initiate a DMA transfer.
 *
 * The following defines the terms used for describing the various kinds of
 * descriptors:
 * - "Torus" means the transfer is between nodes.
 * - "Local" means the transfer is within the same node.
 * - "Direct-put" means the data is put directly into the destination node's
 *   memory.
 * - "MemFifo" means the packets are put into the destination node's reception
 *   fifo.
 * - "Remote-get" means the packet payload contains an injection descriptor
 *   to be injected into the destination node's injection fifo.
 * - Prefetch-only" means the payload is just pre-fetched into L3.  It is not
 *   transferred to the destination node.
 *
 * The following are the functions provided for creating injection descriptors:
 * - DMA_TorusDirectPutDescriptor
 * - DMA_LocalDirectPutDescriptor
 * - DMA_LocalPrefetchOnlyDescriptor
 * - DMA_TorusRemoteGetDescriptor
 * - DMA_LocalRemoteGetDescriptor
 * - DMA_TorusMemFifoDescriptor
 * - DMA_LocalMemFifoDescriptor
 * - DMA_TorusDirectPutBcastDescriptor
 * - DMA_TorusMemFifoBcastDescriptor
 *
 *
 * There are also functions for setting or changing specific values in the
 * injection descriptors.
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




#include <bpcore/bgp_types.h>
#include <common/alignment.h>
#include <common/bgp_bitnumbers.h>
#include <spi/DMA_Packet.h>
#include <spi/DMA_Assert.h>




/*!
 * \brief Packet Header - Checksum Skip Bytes
 *
 * Default number of 2 byte units to skip from the top of a packet before
 * including the packet bytes into the running checksum of the torus
 * injection fifo where this packet is injected.
 *
 * 8 corresponds to skipping 16 bytes, which is the DMA packet header size
 * (hardware header + software header).
 */
#define DMA_CSUM_SKIP  8


/*!
 * \brief Packet Header - Checksum Skip Packet
 *
 * Default value for the torus injection checksum skip packet bit.
 *   - 0 includes the packet (excluding the portion designated by DMA_CSUM_SKIP)
 *     in the checksum.
 *   - 1 excludes the entire packet from the checksum.
 */
#define DMA_CSUM_BIT   0




/*!
 * \brief DMA Injection Descriptor Structure
 *
 */
typedef struct DMA_InjDescriptor_t
{
  union {
    unsigned word1;                 /*!< For accessing fields as 32-bit word  */

    struct {
      unsigned rsvd0          : 24; /*!< 3 bytes: unused                      */

      unsigned rsvd1          :  6; /*!< Bits 0-5: unused flags               */

      unsigned prefetch_only  :  1; /*!< Bit 6: prefetch only, on local       
                                         memcopy:                             
                                         0 = Data is both read and written,   
                                         1 = Data is only read.               
                                         This bit is ignored for torus        
                                         packets.                             */

      unsigned local_memcopy  :  1; /*!< Bit 7: local memory copy bit:        
                                         0 = The message is a torus message,  
                                         1 = The message is a local copy.     */
    };
  };

  union {
    unsigned word2;                 /*!< For accessing fields as 32-bit word  */

    struct {
      unsigned rsvd2          : 24; /*!< 3 bytes: unused                      */

      unsigned idma_counterId :  8; /*!< 1 byte: Injection Counter Id.        */
    };
  };

  unsigned base_offset        : 32; /*!< 4 bytes: pointer to base address of  
                                         message payload.  This gets added to 
                                         the base address associated with the 
                                         idma_counterId injection counter.    */

  unsigned msg_length         : 32; /*!< 4 bytes: message length (in bytes)   */

  DMA_PacketHeader_t hwHdr;         /*!< DMA Hardware Packet Header           */

}
DMA_InjDescriptor_t ALIGN_QUADWORD;
/*!
 * \todo Change to ALIGN_L1D_CACHE when it works.
 *
 */


/*!
 * \brief Static Info from Personality
 *
 * The following structure defines information from the personality.
 * It is intended to be static so, once the info is retrieved from
 * the personality, it does not need to be retrieved again (it is a
 * system call to retrieve personality info).
 *
 */
typedef struct DMA_PersonalityInfo_t
{
  unsigned int personalityRetrieved; /*!< 0 = Personality Info not            
                                              retrieved into this             
                                              structure yet.                  
                                          1 = Personality Info in this        
                                              structure is valid.             */
  uint8_t      nodeXCoordinate;      /*!< X coord of the calling node.        */
  uint8_t      nodeYCoordinate;      /*!< Y coord of the calling node.        */
  uint8_t      nodeZCoordinate;      /*!< Z coord of the calling node.        */
  uint8_t      xNodes;               /*!< X dimension of the block.           */
  uint8_t      yNodes;               /*!< Y dimension of the block.           */
  uint8_t      zNodes;               /*!< Z dimension of the block.           */
}
DMA_PersonalityInfo_t;


/*!
 * \brief Create a DMA Descriptor For a Torus Direct Put Message
 *
 * A torus direct put message is one that is sent to another node and its data
 * is directly put into memory by the DMA on the destination node...it does
 * not go into a reception fifo.
 *
 * A torus direct-put DMA descriptor contains the following:
 *
 * - 16 bytes of control information:
 *   - prefetch_only   = 0
 *   - local_memcopy   = 0
 *   - idma_counterId  = Injection counter ID associated with the data being
 *                       sent.  This counter contains the base address of the
 *                       message and the message length.  Set based on caller's
 *                       inj_ctr_grp_id and inj_ctr_id.
 *   - base_offset     = Message offset (from the injection counter's base
 *                       address).  Set to caller's send_offset.
 *   - msg_length      = Message length.  Set to caller's msg_len.
 *
 * - 8 byte torus hardware header
 *   - CSum_Skip       = DMA_CSUM_SKIP.
 *   - Sk              = DMA_CSUM_BIT.
 *   - Hint            = Set to caller's "hints".
 *   - Pid0, Pid1      = Set based on caller's "recv_ctr_grp_id" (see note).
 *   - Chunks          = Set to largest size consistent with msg_len.
 *   - Dm              = 1 (Indicates a direct-put packet).
 *   - Dynamic         = Set based on caller's "vc".
 *   - VC              = Set to caller's "vc".
 *   - X,Y,Z           = Set to caller's "x", "y", "z".
 *
 * - 8 byte software header (initial values used by iDMA).
 *   - Put_Offset      = Destination message offset (from the reception
 *                       counter's base address).  Set to caller's recv_offset.
 *   - rDMA_Counter    = Reception counter ID.  This counter is located on the
 *                       destination node and contains the base address of the
 *                       message and the message length.  Set based on caller's
 *                       recv_ctr_grp_id and recv_ctr_id.
 *   - Payload_Bytes   = Number of valid bytes in the payload.  Set by iDMA.
 *   - Flags           = Pacing     = 0.
 *                       Remote-Get = 0.
 *   - iDMA_Fifo_ID    = 0 (not used).
 *   - Func_Id         = 0 (not used).
 *
 * This function creates the above descriptor.
 *
 * \param[in,out]  desc             Pointer to the storage where the descriptor
 *                                  will be created.
 * \param[in]      x                The destination's x coordinate (8 bits).
 * \param[in]      y                The destination's y coordinate (8 bits).
 * \param[in]      z                The destination's z coordinate (8 bits).
 * \param[in]      hints            Hint bits for torus routing (6 bits).
 *                                  Each bit corresponds to x+, x-, y+, y-,
 *                                  z+, z-.  If a bit is set, it indicates that
 *                                  the packet wants to travel along the
 *                                  corresponding direction.  If all bits are
 *                                  zero, the hardware calculates the hint bits.
 *                                  Both of x+ and x- cannot be set at the same
 *                                  time...same with y and z.
 * \param[in]      vc               The virtual channel that the packet must go
 *                                  into if it fails to win the bypass
 *                                  arbitration in the receiving node.
 *                                  - 0 = Virtual channel dynamic 0
 *                                  - 1 = Virtual channel dynamic 1
 *                                  - 2 = Virtual channel deterministic bubble
 *                                  - 3 = Virtual channel deterministic priority
 * \param[in]      inj_ctr_grp_id   Injection counter group ID
 *                                  (0 to DMA_NUM_COUNTER_GROUPS-1).
 * \param[in]      inj_ctr_id       Injection counter ID (within the inj counter
 *                                  group) (0 to DMA_NUM_COUNTERS_PER_GROUP-1).
 * \param[in]      send_offset      Offset of the send payload from the pa_base
 *                                  associated with the specified injection
 *                                  counter.
 * \param[in]      recv_ctr_grp_id  Reception counter group ID
 *                                  (0 to DMA_NUM_COUNTER_GROUPS-1).
 * \param[in]      recv_ctr_id      Reception counter ID (within the recv counter
 *                                  group) (0 to DMA_NUM_COUNTERS_PER_GROUP-1).
 * \param[in]      recv_offset      Offset of the payload from the pa_base
 *                                  associated with the specified reception
 *                                  counter.
 * \param[in]      msg_len          Total message length (in bytes).
 *
 * \retval  0         Success
 * \retval  non-zero  Failure
 *
 * \note By default, all payload bytes are included in the torus injection
 *       checksum.  In the first byte of the torus hardware packet header,
 *       this corresponds to setting CSum_Skip = 0x8 (16 bytes) and Sk=0.
 *       The defaults can be changed by changing DMA_CSUM_SKIP and
 *       DMA_CSUM_BIT in this include file.
 *
 * \note By default, the packet size is set to the largest value consistent
 *       with the message size.  For example,
 *       - if msg_len >= 209, there will be 8 32-byte chunks in each packet,
 *         with the possible exception of the last packet, which could contain
 *         fewer chunks (209... of payload + 16 header).
 *       - if 177 <= msg_len < 208, there will be 7 chunks in the packet, etc.
 *
 * \note By default, for direct-put DMA messages, the pid0 and pid1 bits in the
 *       torus hardware packet header are determined by the recv_ctr_grp_id:
 *       - if recv_ctr_grp_id = 0 => (pid0,pid1) = (0,0)
 *       - if recv_ctr_grp_id = 1 => (pid0,pid1) = (0,1)
 *       - if recv_ctr_grp_id = 2 => (pid0,pid1) = (1,0)
 *       - if recv_ctr_grp_id = 3 => (pid0,pid1) = (1,1)
 *       Pid0 determines into which physical torus fifo group on the destination
 *       node the packet is put, prior to the dma receiving it.  Other than that,
 *       the only use for the pid bits is for debug, ie, if headers are being
 *       saved.
 */
int  DMA_TorusDirectPutDescriptor(
				  DMA_InjDescriptor_t *desc,
				  unsigned int         x,
				  unsigned int         y,
				  unsigned int         z,
				  unsigned int         hints,
				  unsigned int         vc,
				  unsigned int         inj_ctr_grp_id,
				  unsigned int         inj_ctr_id,
				  unsigned int         send_offset,
				  unsigned int         recv_ctr_grp_id,
				  unsigned int         recv_ctr_id,
				  unsigned int         recv_offset,
				  unsigned int         msg_len
				 );


/*!
 * \brief Create a DMA Descriptor For a Local Direct Put Message
 *
 * A local direct put message is one that is targeted within the same node, and
 * its data is directly put into memory by the DMA...it does not go into a
 * reception fifo.  This is essentially a memcpy via DMA.
 *
 * A local direct-put DMA descriptor contains the following:
 *
 * - 16 bytes of control information:
 *   - prefetch_only   = 0
 *   - local_memcopy   = 1
 *   - idma_counterId  = Injection counter ID associated with the data being
 *                       sent.  This counter contains the base address of the
 *                       message and the message length.  Set based on caller's
 *                       inj_ctr_grp_id and inj_ctr_id.
 *   - base_offset     = Message offset (from the injection counter's base
 *                       address).  Set to caller's send_offset.
 *   - msg_length      = Message length.  Set to caller's msg_len.
 *
 * - 8 byte torus hardware header
 *   - CSum_Skip       = 0 (not used).
 *   - Sk              = 0 (not used).
 *   - Hint            = 0 (not used).
 *   - Pid0, Pid1      = Set based on caller's "recv_ctr_grp_id".
 *   - Chunks          = Set to largest size consistent with msg_len.
 *   - Dm              = 1 (Indicates a direct-put packet).
 *   - Dynamic         = 0 (not used).
 *   - VC              = 0 (not used).
 *   - X,Y,Z           = 0 (not used).
 *
 * - 8 byte software header (initial values used by iDMA).
 *   - Put_Offset      = Destination message offset (from the reception
 *                       counter's base address).  Set to caller's recv_offset.
 *   - rDMA_Counter    = Reception counter ID.  This counter is located on the
 *                       destination node and contains the base address of the
 *                       message and the message length..  Set based on caller's
 *                       recv_ctr_grp_id and recv_ctr_id.
 *   - Payload_Bytes   = Number of valid bytes in the payload.  Set by iDMA.
 *   - Flags           = Pacing     = 0.
 *                       Remote-Get = 0.
 *   - iDMA_Fifo_ID    = 0 (not used).
 *   - Func_Id         = 0 (not used).
 *
 * This function creates the above descriptor.
 *
 * \param[in,out]  desc             Pointer to the storage where the descriptor
 *                                  will be created.
 * \param[in]      inj_ctr_grp_id   Injection counter group ID
 *                                  (0 to DMA_NUM_COUNTER_GROUPS-1).
 * \param[in]      inj_ctr_id       Injection counter ID (within the inj counter
 *                                  group) (0 to DMA_NUM_COUNTERS_PER_GROUP-1).
 * \param[in]      send_offset      Offset of the send payload from the pa_base
 *                                  associated with the specified injection
 *                                  counter.
 * \param[in]      recv_ctr_grp_id  Reception counter group ID
 *                                  (0 to DMA_NUM_COUNTER_GROUPS-1).
 * \param[in]      recv_ctr_id      Reception counter ID (within the recv counter
 *                                  group) (0 to DMA_NUM_COUNTERS_PER_GROUP-1).
 * \param[in]      recv_offset      Offset of the payload from the pa_base
 *                                  associated with the specified reception
 *                                  counter.
 * \param[in]      msg_len          Total message length (in bytes).
 *
 * \retval  0         Success
 * \retval  non-zero  Failure
 *
 * \note By default, the packet size is set to the largest value consistent
 *       with the message size.  For example,
 *       - if msg_len >= 209, there will be 8 32-byte chunks in each packet,
 *         with the possible exception of the last packet, which could contain
 *         fewer chunks (209... of payload + 16 header).
 *       - if 177 <= msg_len < 208, there will be 7 chunks in the packet, etc.
 *
 * \note By default, for direct-put DMA messages, the pid0 and pid1 bits in the
 *       torus hardware packet header are determined by the recv_ctr_grp_id:
 *       - if recv_ctr_grp_id = 0 => (pid0,pid1) = (0,0)
 *       - if recv_ctr_grp_id = 1 => (pid0,pid1) = (0,1)
 *       - if recv_ctr_grp_id = 2 => (pid0,pid1) = (1,0)
 *       - if recv_ctr_grp_id = 3 => (pid0,pid1) = (1,1)
 *       The only use for the pid bits is for debug, ie, if headers are
 *       being saved.
 */
int  DMA_LocalDirectPutDescriptor(
				  DMA_InjDescriptor_t *desc,
				  unsigned int         inj_ctr_grp_id,
				  unsigned int         inj_ctr_id,
				  unsigned int         send_offset,
				  unsigned int         recv_ctr_grp_id,
				  unsigned int         recv_ctr_id,
				  unsigned int         recv_offset,
				  unsigned int         msg_len
				 );


/*!
 * \brief Create a DMA Descriptor For a Local L3 Prefetch Only Message
 *
 * A local prefetch is one in which the DMA simply prefetches the send buffer
 * into L3.
 *
 * A local prefetch DMA descriptor contains the following:
 *
 * - 16 bytes of control information:
 *   - prefetch_only   = 1
 *   - local_memcopy   = 1
 *   - idma_counterId  = Injection counter ID associated with the message being
 *                       prefetched.  This counter contains the base address of
 *                       the message and the message length.  Set based on caller's
 *                       inj_ctr_grp_id and inj_ctr_id.
 *   - base_offset     = Message offset (from the injection counter's base
 *                       address).  Set to caller's send_offset.
 *   - msg_length      = Message length.  Set to caller's msg_len.
 *
 * - 8 byte torus hardware header
 *   - CSum_Skip       = 0 (not used).
 *   - Sk              = 0 (not used).
 *   - Hint            = 0 (not used).
 *   - Pid0, Pid1      = 0 (not used).
 *   - Chunks          = Set to largest size consistent with msg_len.
 *   - Dm              = 1 (Indicates a DMA packet).
 *   - Dynamic         = 0 (not used).
 *   - VC              = 0 (not used).
 *   - X,Y,Z           = 0 (not used).
 *
 * - 8 byte software header (initial values used by iDMA).
 *   - Put_Offset      = 0 (not used).
 *   - rDMA_Counter    = 0 (not used).
 *   - Payload_Bytes   = 0 (not used).
 *   - Flags           = Pacing     = 0.
 *                       Remote-Get = 0.
 *   - iDMA_Fifo_ID    = 0 (not used).
 *   - Func_Id         = 0 (not used).
 *
 * This function creates the above descriptor.
 *
 * \param[in,out]  desc             Pointer to the storage where the descriptor
 *                                  will be created.
 * \param[in]      inj_ctr_grp_id   Injection counter group ID
 *                                  (0 to DMA_NUM_COUNTER_GROUPS-1).
 * \param[in]      inj_ctr_id       Injection counter ID (within the inj counter
 *                                  group) (0 to DMA_NUM_COUNTERS_PER_GROUP-1).
 * \param[in]      send_offset      Offset of the send payload from the pa_base
 *                                  associated with the specified injection
 *                                  counter.
 * \param[in]      msg_len          Total message length (in bytes).
 *
 * \retval  0         Success
 * \retval  non-zero  Failure
 *
 * \note By default, the packet size is set to the largest value consistent
 *       with the message size.  For example,
 *       - if msg_len >= 209, there will be 8 32-byte chunks in each packet,
 *         with the possible exception of the last packet, which could contain
 *         fewer chunks (209... of payload + 16 header).
 *       - if 177 <= msg_len < 208, there will be 7 chunks in the packet, etc.
 *
 */
int  DMA_LocalPrefetchOnlyDescriptor(
				     DMA_InjDescriptor_t *desc,
				     unsigned int         inj_ctr_grp_id,
				     unsigned int         inj_ctr_id,
				     unsigned int         send_offset,
				     unsigned int         msg_len
				    );


/*!
 * \brief Create a DMA Descriptor For a Torus Remote-Get Message
 *
 * A torus remote-get message is one that is sent to another node and its data
 * is directly put by the DMA into an injection fifo on the destination
 * node...it does not go into a reception fifo.  Therefore, the payload of this
 * message is one (or more) descriptors for another message that is to be sent
 * back to the originating node.
 *
 * By default, we assume that the payload of this remote get packet is a single
 * descriptor.  Thus, Chunks = (2)-1 (64 byte packet) and msg_length = 32.
 * For remote gets whose payload is greater than 1 descriptor, the caller can
 * change the packet Chunks and msg_length after this function builds the
 * default descriptor.
 *
 * It is also assumed that the payload is NOT checksummed, since it is not
 * always reproducible.  Things like idma_counterId and base_offset may be
 * different on another run, making checksumming inconsistent.
 *
 * A torus remote-get DMA descriptor contains the following:
 *
 * - 16 bytes of control information:
 *   - prefetch_only   = 0
 *   - local_memcopy   = 0
 *   - idma_counterId  = Injection counter ID associated with the data being
 *                       sent.  This counter contains the base address of the
 *                       message and the message length.  Set based on caller's
 *                       inj_ctr_grp_id and inj_ctr_id.
 *   - base_offset     = Message offset (from the injection counter's base
 *                       address).  Set to caller's send_offset.
 *   - msg_length      = 32.
 *
 * - 8 byte torus hardware header
 *   - CSum_Skip       = 0 (not used because Sk is 1).
 *   - Sk              = 1 (do not checksum this packet).
 *   - Hint            = Set to caller's "hints".
 *   - Pid0, Pid1      = Set based on caller's "recv_inj_fifo_id" (see note).
 *   - Chunks          = Set to (2)-1 = 1.
 *   - Dm              = 1 (Indicates a DMA packet).
 *   - Dynamic         = Set based on caller's "vc".
 *   - VC              = Set to caller's "vc".
 *   - X,Y,Z           = Set to caller's "x", "y", "z".
 *
 * - 8 byte software header (initial values used by iDMA).
 *   - Put_Offset      = 0 (not used).
 *   - rDMA_Counter    = 0 (not used).
 *   - Payload_Bytes   = Number of valid bytes in the payload.  Set by iDMA.
 *   - Flags           = Pacing     = 0.
 *                       Remote-Get = 1.
 *   - iDMA_Fifo_ID    = Injection fifo ID where the payload will be injected.
 *                       Set based on caller's recv_inj_ctr_grp_id and
 *                       recv_inj_ctr_id.
 *   - Func_Id         = 0 (not used).
 *
 * This function creates the above descriptor.
 *
 * \param[in,out]  desc             Pointer to the storage where the descriptor
 *                                  will be created.
 * \param[in]      x                The destination's x coordinate (8 bits).
 * \param[in]      y                The destination's y coordinate (8 bits).
 * \param[in]      z                The destination's z coordinate (8 bits).
 * \param[in]      hints            Hint bits for torus routing (6 bits).
 *                                  Each bit corresponds to x+, x-, y+, y-,
 *                                  z+, z-.  If a bit is set, it indicates that
 *                                  the packet wants to travel along the
 *                                  corresponding direction.  If all bits are
 *                                  zero, the hardware calculates the hint bits.
 *                                  Both of x+ and x- cannot be set at the same
 *                                  time...same with y and z.
 * \param[in]      vc               The virtual channel that the packet must go
 *                                  into if it fails to win the bypass
 *                                  arbitration in the receiving node.
 *                                  - 0 = Virtual channel dynamic 0
 *                                  - 1 = Virtual channel dynamic 1
 *                                  - 2 = Virtual channel deterministic bubble
 *                                  - 3 = Virtual channel deterministic priority
 * \param[in]      inj_ctr_grp_id   Injection counter group ID
 *                                  (0 to DMA_NUM_COUNTER_GROUPS-1).
 * \param[in]      inj_ctr_id       Injection counter ID (within the inj counter
 *                                  group) (0 to DMA_NUM_COUNTERS_PER_GROUP-1).
 * \param[in]      send_offset      Offset of the send payload from the pa_base
 *                                  associated with the specified injection
 *                                  counter.
 * \param[in]      recv_inj_fifo_grp_id  Injection fifo group ID where payload
 *                                       will be injected on destination node
 *                                       (0 to DMA_NUM_INJ_FIFO_GROUPS-1).
 * \param[in]      recv_inj_fifo_id      Injection fifo ID (within the
 *                                       recv_inj_fifo_grp_id group)
 *                                       (0 to DMA_NUM_INJ_FIFOS_PER_GROUP-1).
 *
 * \retval  0         Success
 * \retval  non-zero  Failure
 *
 * \note By default, for remote-get DMA messages, the pid0 and pid1 bits in the
 *       torus hardware packet header are determined by the recv_inj_fifo_grp_id:
 *       - if recv_inj_fifo_grp_id = 0 => (pid0,pid1) = (0,0)
 *       - if recv_inj_fifo_grp_id = 1 => (pid0,pid1) = (0,1)
 *       - if recv_inj_fifo_grp_id = 2 => (pid0,pid1) = (1,0)
 *       - if recv_inj_fifo_grp_id = 3 => (pid0,pid1) = (1,1)
 *       Pid0 determines into which physical torus fifo group on the destination
 *       node the packet is put, prior to the dma receiving it.  Other than that,
 *       the only use for the pid bits is for debug, ie, if headers are being
 *       saved.
 */
int  DMA_TorusRemoteGetDescriptor(
				  DMA_InjDescriptor_t *desc,
				  unsigned int         x,
				  unsigned int         y,
				  unsigned int         z,
				  unsigned int         hints,
				  unsigned int         vc,
				  unsigned int         inj_ctr_grp_id,
				  unsigned int         inj_ctr_id,
				  unsigned int         send_offset,
				  unsigned int         recv_inj_fifo_grp_id,
				  unsigned int         recv_inj_fifo_id
				 );
					

/*!
 * \brief Create a DMA Descriptor For a Local Remote-Get Message
 *
 * A local remote-get message is one whose data is directly put by the DMA into
 * an injection fifo on the local node...it does not go into a reception fifo.
 * Therefore, the payload of this message is one (or more) descriptors for
 * another message that is to be injected on the local node.
 *
 * By default, we assume that the payload of this remote get packet is a single
 * descriptor.  Thus, Chunks = (2)-1 (64 byte packet) and msg_length = 32.
 * For remote gets whose payload is greater than 1 descriptor, the caller can
 * change the packet Chunks and msg_length after this function builds the
 * default descriptor.
 *
 * A local remote-get DMA descriptor contains the following:
 *
 * - 16 bytes of control information:
 *   - prefetch_only   = 0
 *   - local_memcopy   = 1
 *   - idma_counterId  = Injection counter ID associated with the data being
 *                       sent.  This counter contains the base address of the
 *                       message and the message length.  Set based on caller's
 *                       inj_ctr_grp_id and inj_ctr_id.
 *   - base_offset     = Message offset (from the injection counter's base
 *                       address).  Set to caller's send_offset.
 *   - msg_length      = 32.
 *
 * - 8 byte torus hardware header
 *   - CSum_Skip       = 0 (not used).
 *   - Sk              = 0 (not used).
 *   - Hint            = 0 (Set to caller's "hints".
 *   - Pid0, Pid1      = Set based on caller's "recv_inj_fifo_id" (see note).
 *   - Chunks          = Set to (2)-1 = 1.
 *   - Dm              = 1 (Indicates a DMA packet).
 *   - Dynamic         = 0 (not used).
 *   - VC              = 0 (not used).
 *   - X,Y,Z           = 0 (not used).
 *
 * - 8 byte software header (initial values used by iDMA).
 *   - Put_Offset      = 0 (not used).
 *   - rDMA_Counter    = 0 (not used).
 *   - Payload_Bytes   = Number of valid bytes in the payload.  Set by iDMA.
 *   - Flags           = Pacing     = 0.
 *                       Remote-Get = 1.
 *   - iDMA_Fifo_ID    = Injection fifo ID where the payload will be injected.
 *                       Set based on caller's inj_ctr_grp_id and inj_ctr_id.
 *   - Func_Id         = 0 (not used).
 *
 * This function creates the above descriptor.
 *
 * \param[in,out]  desc             Pointer to the storage where the descriptor
 *                                  will be created.
 * \param[in]      inj_ctr_grp_id   Injection counter group ID
 *                                  (0 to DMA_NUM_COUNTER_GROUPS-1).
 * \param[in]      inj_ctr_id       Injection counter ID (within the inj counter
 *                                  group) (0 to DMA_NUM_COUNTERS_PER_GROUP-1).
 * \param[in]      send_offset      Offset of the send payload from the pa_base
 *                                  associated with the specified injection
 *                                  counter.
 * \param[in]      recv_inj_fifo_grp_id  Injection fifo group ID where payload
 *                                       will be injected on local node
 *                                       (0 to DMA_NUM_INJ_FIFO_GROUPS-1).
 * \param[in]      recv_inj_fifo_id      Injection fifo ID (within the
 *                                       recv_inj_fifo_grp_id group)
 *                                       (0 to DMA_NUM_INJ_FIFOS_PER_GROUP-1).
 *
 * \retval  0         Success
 * \retval  non-zero  Failure
 *
 * \note By default, for remote-get DMA messages, the pid0 and pid1 bits in the
 *       hardware packet header are determined by the recv_inj_fifo_grp_id:
 *       - if recv_inj_fifo_grp_id = 0 => (pid0,pid1) = (0,0)
 *       - if recv_inj_fifo_grp_id = 1 => (pid0,pid1) = (0,1)
 *       - if recv_inj_fifo_grp_id = 2 => (pid0,pid1) = (1,0)
 *       - if recv_inj_fifo_grp_id = 3 => (pid0,pid1) = (1,1)
 *
 */
int  DMA_LocalRemoteGetDescriptor(
				  DMA_InjDescriptor_t *desc,
				  unsigned int         inj_ctr_grp_id,
				  unsigned int         inj_ctr_id,
				  unsigned int         send_offset,
				  unsigned int         recv_inj_fifo_grp_id,
				  unsigned int         recv_inj_fifo_id
				 );

					
/*!
 * \brief Create a DMA Descriptor For a Torus Memory Fifo Message
 *
 * A torus memory fifo message is one that is sent to another node and its data
 * is put into a reception memory fifo by the DMA on the destination node.
 *
 * A torus memory fifo DMA descriptor contains the following:
 *
 * - 16 bytes of control information:
 *   - prefetch_only   = 0
 *   - local_memcopy   = 0
 *   - idma_counterId  = Injection counter ID associated with the data being
 *                       sent.  This counter contains the base address of the
 *                       message and the message length.  Set based on caller's
 *                       inj_ctr_grp_id and inj_ctr_id.
 *   - base_offset     = Message offset (from the injection counter's base
 *                       address).  Set to caller's send_offset.
 *   - msg_length      = Message length.  Set to caller's msg_len.
 *
 * - 8 byte torus hardware header
 *   - CSum_Skip       = DMA_CSUM_SKIP.
 *   - Sk              = DMA_CSUM_BIT.
 *   - Hint            = Set to caller's "hints".
 *   - Pid0, Pid1      = Set based on caller's "recv_fifo_grp_id" (see note).
 *   - Chunks          = Set to largest size consistent with msg_len.
 *   - Dm              = 0 (Indicates a memory fifo packet).
 *   - Dynamic         = Set based on caller's "vc".
 *   - VC              = Set to caller's "vc".
 *   - X,Y,Z           = Set to caller's "x", "y", "z".
 *
 * - 8 byte software header (initial values used by iDMA).
 *   - Put_Offset      = 0 (not used).
 *   - rDMA_Counter    = 0 (not used).
 *   - Payload_Bytes   = 0 (not used).
 *   - Flags           = Pacing     = 0.
 *                       Remote-Get = 0.
 *   - iDMA_Fifo_ID    = 0 (not used).
 *   - SW_Arg          = User-defined 24 bits.  Set to caller's sw_arg.
 *   - Func_Id         = The registration ID of a function to receive control
 *                       on the destination node to process the packet.
 *                       Set to caller's function_id.
 *
 * This function creates the above descriptor.
 *
 * \param[in,out]  desc             Pointer to the storage where the descriptor
 *                                  will be created.
 * \param[in]      x                The destination's x coordinate (8 bits).
 * \param[in]      y                The destination's y coordinate (8 bits).
 * \param[in]      z                The destination's z coordinate (8 bits).
 * \param[in]      recv_fifo_grp_id Reception fifo group ID
 *                                  (0 to DMA_NUM_REC_FIFO_GROUPS-1).
 * \param[in]      hints            Hint bits for torus routing (6 bits).
 *                                  Each bit corresponds to x+, x-, y+, y-,
 *                                  z+, z-.  If a bit is set, it indicates that
 *                                  the packet wants to travel along the
 *                                  corresponding direction.  If all bits are
 *                                  zero, the hardware calculates the hint bits.
 *                                  Both of x+ and x- cannot be set at the same
 *                                  time...same with y and z.
 * \param[in]      vc               The virtual channel that the packet must go
 *                                  into if it fails to win the bypass
 *                                  arbitration in the receiving node.
 *                                  - 0 = Virtual channel dynamic 0
 *                                  - 1 = Virtual channel dynamic 1
 *                                  - 2 = Virtual channel deterministic bubble
 *                                  - 3 = Virtual channel deterministic priority
 * \param[in]      sw_arg           User-defined 24 bits to be placed into the
 *                                  packets (bits 8-31).
 * \param[in]      function_id      Function id (8 bit registration ID) of the
 *                                  function to receive control on the
 *                                  destination node to process packets for this
 *                                  message.
 * \param[in]      inj_ctr_grp_id   Injection counter group ID
 *                                  (0 to DMA_NUM_COUNTER_GROUPS-1).
 * \param[in]      inj_ctr_id       Injection counter ID (within the inj counter
 *                                  group) (0 to DMA_NUM_COUNTERS_PER_GROUP-1).
 * \param[in]      send_offset      Offset of the send payload from the pa_base
 *                                  associated with the specified injection
 *                                  counter.
 * \param[in]      msg_len          Total message length (in bytes).
 *
 * \retval  0         Success
 * \retval  non-zero  Failure
 *
 * \note By default, all payload bytes are included in the torus injection
 *       checksum.  In the first byte of the torus hardware packet header,
 *       this corresponds to setting CSum_Skip = 0x8 (16 bytes) and Sk=0.
 *       The defaults can be changed by changing DMA_CSUM_SKIP and
 *       DMA_CSUM_BIT in this include file.
 *
 * \note By default, the packet size is set to the largest value consistent
 *       with the message size.  For example,
 *       - if msg_len >= 209, there will be 8 32-byte chunks in each packet,
 *         with the possible exception of the last packet, which could contain
 *         fewer chunks (209... of payload + 16 header).
 *       - if 177 <= msg_len < 208, there will be 7 chunks in the packet, etc.
 *
 * \note By default, for DMA messages, the pid0 and pid1 bits in the
 *       torus hardware packet header are determined by the recv_fifo_grp_id:
 *       - if recv_fifo_grp_id = 0 => (pid0,pid1) = (0,0)
 *       - if recv_fifo_grp_id = 1 => (pid0,pid1) = (0,1)
 *       - if recv_fifo_grp_id = 2 => (pid0,pid1) = (1,0)
 *       - if recv_fifo_grp_id = 3 => (pid0,pid1) = (1,1)
 *       Pid0 determines into which physical torus fifo group on the destination
 *       node the packet is put, prior to the dma receiving it.  Other than that,
 *       the only use for the pid bits is for debug, ie, if headers are being
 *       saved.
*/
int  DMA_TorusMemFifoDescriptor(
				DMA_InjDescriptor_t *desc,
				unsigned int         x,
				unsigned int         y,
				unsigned int         z,
				unsigned int         recv_fifo_grp_id,
				unsigned int         hints,
				unsigned int         vc,
				unsigned int         sw_arg,
				unsigned int         function_id,
				unsigned int         inj_ctr_grp_id,
				unsigned int         inj_ctr_id,
				unsigned int         send_offset,
				unsigned int         msg_len
			       );


/*!
 * \brief Create a DMA Descriptor For a Local Memory Fifo Message
 *
 * A local memory fifo message is one whose data is put into a reception
 * memory fifo on the same node by the DMA.
 *
 * A local memory fifo DMA descriptor contains the following:
 *
 * - 16 bytes of control information:
 *   - prefetch_only   = 0
 *   - local_memcopy   = 0
 *   - idma_counterId  = Injection counter ID associated with the data being
 *                       sent.  This counter contains the base address of the
 *                       message and the message length.  Set based on caller's
 *                       inj_ctr_grp_id and inj_ctr_id.
 *   - base_offset     = Message offset (from the injection counter's base
 *                       address).  Set to caller's send_offset.
 *   - msg_length      = Message length.  Set to caller's msg_len.
 *
 * - 8 byte torus hardware header
 *   - CSum_Skip       = 0 (not used).
 *   - Sk              = 0 (not used).
 *   - Hint            = 0 (not used).
 *   - Pid0, Pid1      = Set based on caller's "recv_fifo_grp_id" (see note).
 *   - Chunks          = Set to largest size consistent with msg_len.
 *   - Dm              = 0 (Indicates a memory fifo packet).
 *   - Dynamic         = 0 (not used).
 *   - VC              = 0 (not used).
 *   - X,Y,Z           = 0 (not used).
 *
 * - 8 byte software header (initial values used by iDMA).
 *   - Put_Offset      = 0 (not used).
 *   - rDMA_Counter    = 0 (not used).
 *   - Payload_Bytes   = 0 (not used).
 *   - Flags           = Pacing     = 0.
 *                       Remote-Get = 0.
 *   - iDMA_Fifo_ID    = 0 (not used).
 *   - SW_Arg          = User-defined 24 bits.  Set to caller's sw_arg.
 *   - Func_Id         = The registration ID of a function to receive control
 *                       on this local node to process the packet.
 *                       Set to caller's function_id.
 *
 * This function creates the above descriptor.
 *
 * \param[in,out]  desc             Pointer to the storage where the descriptor
 *                                  will be created.
 * \param[in]      recv_fifo_grp_id Reception fifo group ID
 *                                  (0 to DMA_NUM_REC_FIFO_GROUPS-1).
 * \param[in]      sw_arg           User-defined 24 bits to be placed into the
 *                                  packets (bits 8-31).
 * \param[in]      function_id      Function id (8 bit registration ID) of the
 *                                  function to receive control on this
 *                                  local node to process packets for this
 *                                  message.
 * \param[in]      inj_ctr_grp_id   Injection counter group ID
 *                                  (0 to DMA_NUM_COUNTER_GROUPS-1).
 * \param[in]      inj_ctr_id       Injection counter ID (within the inj counter
 *                                  group) (0 to DMA_NUM_COUNTERS_PER_GROUP-1).
 * \param[in]      send_offset      Offset of the send payload from the pa_base
 *                                  associated with the specified injection
 *                                  counter.
 * \param[in]      msg_len          Total message length (in bytes).
 *
 * \retval  0         Success
 * \retval  non-zero  Failure
 *
 * \note By default, the packet size is set to the largest value consistent
 *       with the message size.  For example,
 *       - if msg_len >= 209, there will be 8 32-byte chunks in each packet,
 *         with the possible exception of the last packet, which could contain
 *         fewer chunks (209... of payload + 16 header).
 *       - if 177 <= msg_len < 208, there will be 7 chunks in the packet, etc.
 *
 * \note By default, for direct-put DMA messages, the pid0 and pid1 bits in the
 *       torus hardware packet header are determined by the recv_fifo_grp_id:
 *       - if recv_fifo_grp_id = 0 => (pid0,pid1) = (0,0)
 *       - if recv_fifo_grp_id = 1 => (pid0,pid1) = (0,1)
 *       - if recv_fifo_grp_id = 2 => (pid0,pid1) = (1,0)
 *       - if recv_fifo_grp_id = 3 => (pid0,pid1) = (1,1)
*/
int  DMA_LocalMemFifoDescriptor(
				DMA_InjDescriptor_t *desc,
				unsigned int         recv_fifo_grp_id,
				unsigned int         sw_arg,
				unsigned int         function_id,
				unsigned int         inj_ctr_grp_id,
				unsigned int         inj_ctr_id,
				unsigned int         send_offset,
				unsigned int         msg_len
			       );


/*!
 * \brief Create a DMA Descriptor For a Torus Direct Put Broadcast Message
 *
 * A torus direct put broadcast message is one that is sent to all of the nodes
 * in a specified direction along a specified line, its data
 * is directly put into memory on the nodes along that line by the DMA on those
 * nodes...it does not go into a reception fifo.  Only one hint bit can be
 * specified, dictating the direction (plus or minus) and line (x, y, or z).
 *
 * By default, the packet is included in the checksum.  Retransmitted packets
 * should not be included in the checksum.
 *
 * By default, the deterministic bubble normal virtual channel is used.
 *
 * A torus direct-put broadcast DMA descriptor contains the following:
 *
 * - 16 bytes of control information:
 *   - prefetch_only   = 0
 *   - local_memcopy   = 0
 *   - idma_counterId  = Injection counter ID associated with the data being
 *                       sent.  This counter contains the base address of the
 *                       message and the message length.  Set based on caller's
 *                       inj_ctr_grp_id and inj_ctr_id.
 *   - base_offset     = Message offset (from the injection counter's base
 *                       address).  Set to caller's send_offset.
 *   - msg_length      = Message length.  Set to caller's msg_len.
 *
 * - 8 byte torus hardware header
 *   - CSum_Skip       = DMA_CSUM_SKIP.
 *   - Sk              = DMA_CSUM_BIT.
 *   - Hint            = Set to caller's "hints".
 *   - Pid0, Pid1      = Set based on caller's "recv_ctr_grp_id" (see note).
 *   - Chunks          = Set to largest size consistent with msg_len.
 *   - Dm              = 1 (Indicates a direct-put packet).
 *   - Dynamic         = 0 (Deterministic).
 *   - VC              = Virtual Channel: Deterministic Bubble Normal.
 *   - X,Y,Z           = Set according to the hints:
 *                       Two of the directions are set to this node's
 *                       coordinates (no movement in those directions).
 *                       One direction is set to the dest specified
 *                       by the caller.
 *
 * - 8 byte software header (initial values used by iDMA).
 *   - Put_Offset      = Destination message offset (from the reception
 *                       counter's base address).  Set to caller's recv_offset.
 *   - rDMA_Counter    = Reception counter ID.  This counter is located on the
 *                       destination node and contains the base address of the
 *                       message and the message length.  Set based on caller's
 *                       recv_ctr_grp_id and recv_ctr_id.
 *   - Payload_Bytes   = Number of valid bytes in the payload.  Set by iDMA.
 *   - Flags           = Pacing     = 0.
 *                       Remote-Get = 0.
 *   - iDMA_Fifo_ID    = 0 (not used).
 *   - Func_Id         = 0 (not used).
 *
 * This function creates the above descriptor.
 *
 * \param[in,out]  desc             Pointer to the storage where the descriptor
 *                                  will be created.
 * \param[in]      dest             The final torus destination coordinate
 *                                  along the line specified by the hints.
 *                                  Should not exceed the number of nodes in
 *                                  the direction of travel.
 * \param[in]      hints            Hint bits for torus routing (6 bits).
 *                                  Each bit corresponds to x+, x-, y+, y-,
 *                                  z+, z-.  If a bit is set, it indicates that
 *                                  the packet wants to travel along the
 *                                  corresponding direction.  If all bits are
 *                                  zero, the hardware calculates the hint bits.
 *                                  Only one bit may be specified.
 * \param[in]      inj_ctr_grp_id   Injection counter group ID
 *                                  (0 to DMA_NUM_COUNTER_GROUPS-1).
 * \param[in]      inj_ctr_id       Injection counter ID (within the inj counter
 *                                  group) (0 to DMA_NUM_COUNTERS_PER_GROUP-1).
 * \param[in]      send_offset      Offset of the send payload from the pa_base
 *                                  associated with the specified injection
 *                                  counter.
 * \param[in]      recv_ctr_grp_id  Reception counter group ID
 *                                  (0 to DMA_NUM_COUNTER_GROUPS-1).
 * \param[in]      recv_ctr_id      Reception counter ID (within the recv counter
 *                                  group) (0 to DMA_NUM_COUNTERS_PER_GROUP-1).
 * \param[in]      recv_offset      Offset of the payload from the pa_base
 *                                  associated with the specified reception
 *                                  counter.
 * \param[in]      msg_len          Total message length (in bytes).
 *
 * \retval  0         Success
 * \retval  non-zero  Failure
 *
 * \note By default, all payload bytes are included in the torus injection
 *       checksum.  In the first byte of the torus hardware packet header,
 *       this corresponds to setting CSum_Skip = 0x8 (16 bytes) and Sk=0.
 *       The defaults can be changed by changing DMA_CSUM_SKIP and
 *       DMA_CSUM_BIT in this include file.
 *
 * \note By default, the packet size is set to the largest value consistent
 *       with the message size.  For example,
 *       - if msg_len >= 209, there will be 8 32-byte chunks in each packet,
 *         with the possible exception of the last packet, which could contain
 *         fewer chunks (209... of payload + 16 header).
 *       - if 177 <= msg_len < 208, there will be 7 chunks in the packet, etc.
 *
 * \note By default, for direct-put DMA messages, the pid0 and pid1 bits in the
 *       torus hardware packet header are determined by the recv_ctr_grp_id:
 *       - if recv_ctr_grp_id = 0 => (pid0,pid1) = (0,0)
 *       - if recv_ctr_grp_id = 1 => (pid0,pid1) = (0,1)
 *       - if recv_ctr_grp_id = 2 => (pid0,pid1) = (1,0)
 *       - if recv_ctr_grp_id = 3 => (pid0,pid1) = (1,1)
 *       Pid0 determines into which physical torus fifo group on the destination
 *       node the packet is put, prior to the dma receiving it.  Other than that,
 *       the only use for the pid bits is for debug, ie, if headers are being
 *       saved.
*/
int  DMA_TorusDirectPutBcastDescriptor(
				       DMA_InjDescriptor_t *desc,
				       unsigned int         dest,
				       unsigned int         hints,
				       unsigned int         inj_ctr_grp_id,
				       unsigned int         inj_ctr_id,
				       unsigned int         send_offset,
				       unsigned int         recv_ctr_grp_id,
				       unsigned int         recv_ctr_id,
				       unsigned int         recv_offset,
				       unsigned int         msg_len
				      );


/*!
 * \brief Create a DMA Descriptor For a Torus Memory Fifo Broadcast Message
 *
 * A torus memory fifo broadcast message is one that is sent to all of the nodes
 * in a specified direction along a specified line, its data is
 * put into a reception memory fifo by the DMA on the destination nodes along
 * that line.  Only one hint bit can be specified, dictating the direction
 * (plus or minus) and line (x, y, or z).
 *
 * By default, the packet is included in the checksum.  Retransmitted packets
 * should not be included in the checksum.
 *
 * By default, the deterministic bubble normal virtual channel is used.
 *
 * A torus memory fifo broadcast DMA descriptor contains the following:
 *
 * - 16 bytes of control information:
 *   - prefetch_only   = 0
 *   - local_memcopy   = 0
 *   - idma_counterId  = Injection counter ID associated with the data being
 *                       sent.  This counter contains the base address of the
 *                       message and the message length.  Set based on caller's
 *                       inj_ctr_grp_id and inj_ctr_id.
 *   - base_offset     = Message offset (from the injection counter's base
 *                       address).  Set to caller's send_offset.
 *   - msg_length      = Message length.  Set to caller's msg_len.
 *
 * - 8 byte torus hardware header
 *   - CSum_Skip       = DMA_CSUM_SKIP.
 *   - Sk              = DMA_CSUM_BIT.
 *   - Hint            = Set to caller's "hints".
 *   - Pid0, Pid1      = Set based on caller's "recv_fifo_grp_id" (see note).
 *   - Chunks          = Set to largest size consistent with msg_len.
 *   - Dm              = 0 (Indicates a memory fifo packet).
 *   - Dynamic         = 0 (Deterministic).
 *   - VC              = Virtual Channel: Deterministic Bubble Normal.
 *   - X,Y,Z           = Set according to the hints:
 *                       Two of the directions are set to this node's
 *                       coordinates (no movement in those directions).
 *                       One direction is set to the dest specified
 *                       by the caller.
 *
 * - 8 byte software header (initial values used by iDMA).
 *   - Put_Offset      = 0 (not used).
 *   - rDMA_Counter    = 0 (not used).
 *   - Payload_Bytes   = 0 (not used).
 *   - Flags           = Pacing     = 0.
 *                       Remote-Get = 0.
 *   - iDMA_Fifo_ID    = 0 (not used).
 *   - SW_Arg          = User-defined 24 bits.  Set to caller's sw_arg.
 *   - Func_Id         = The registration ID of a function to receive control
 *                       on the destination node to process the packet.
 *                       Set to caller's function_id.
 *
 * This function creates the above descriptor.
 *
 * \param[in,out]  desc             Pointer to the storage where the descriptor
 *                                  will be created.
 * \param[in]      dest             The final torus destination coordinate
 *                                  along the line specified by the hints.
 *                                  Should not exceed the number of nodes in
 *                                  the direction of travel.
 * \param[in]      recv_fifo_grp_id Reception fifo group ID
 *                                  (0 to DMA_NUM_REC_FIFO_GROUPS-1).
 * \param[in]      hints            Hint bits for torus routing (6 bits).
 *                                  Each bit corresponds to x+, x-, y+, y-,
 *                                  z+, z-.  If a bit is set, it indicates that
 *                                  the packet wants to travel along the
 *                                  corresponding direction.  If all bits are
 *                                  zero, the hardware calculates the hint bits.
 *                                  Only one bit may be specified.
 * \param[in]      sw_arg           User-defined 24 bits to be placed into the
 *                                  packets (bits 8-31).
 * \param[in]      function_id      Function id (8 bit registration ID) of the
 *                                  function to receive control on the
 *                                  destination node to process packets for this
 *                                  message.
 * \param[in]      inj_ctr_grp_id   Injection counter group ID
 *                                  (0 to DMA_NUM_COUNTER_GROUPS-1).
 * \param[in]      inj_ctr_id       Injection counter ID (within the inj counter
 *                                  group) (0 to DMA_NUM_COUNTERS_PER_GROUP-1).
 * \param[in]      send_offset      Offset of the send payload from the pa_base
 *                                  associated with the specified injection
 *                                  counter.
 * \param[in]      msg_len          Total message length (in bytes).
 *
 * \retval  0         Success
 * \retval  non-zero  Failure
 *
 * \note By default, all payload bytes are included in the torus injection
 *       checksum.  In the first byte of the torus hardware packet header,
 *       this corresponds to setting CSum_Skip = 0x8 (16 bytes) and Sk=0.
 *       The defaults can be changed by changing DMA_CSUM_SKIP and
 *       DMA_CSUM_BIT in this include file.
 *
 * \note By default, the packet size is set to the largest value consistent
 *       with the message size.  For example,
 *       - if msg_len >= 209, there will be 8 32-byte chunks in each packet,
 *         with the possible exception of the last packet, which could contain
 *         fewer chunks (209... of payload + 16 header).
 *       - if 177 <= msg_len < 208, there will be 7 chunks in the packet, etc.
 *
 * \note By default, for direct-put DMA messages, the pid0 and pid1 bits in the
 *       torus hardware packet header are determined by the recv_fifo_grp_id:
 *       - if recv_fifo_grp_id = 0 => (pid0,pid1) = (0,0)
 *       - if recv_fifo_grp_id = 1 => (pid0,pid1) = (0,1)
 *       - if recv_fifo_grp_id = 2 => (pid0,pid1) = (1,0)
 *       - if recv_fifo_grp_id = 3 => (pid0,pid1) = (1,1)
 *       Pid0 determines into which physical torus fifo group on the destination
 *       node the packet is put, prior to the dma receiving it.  Other than that,
 *       the only use for the pid bits is for debug, ie, if headers are being
 *       saved.
*/
int  DMA_TorusMemFifoBcastDescriptor(
				     DMA_InjDescriptor_t *desc,
				     unsigned int         dest,
				     unsigned int         recv_fifo_grp_id,
				     unsigned int         hints,
				     unsigned int         sw_arg,
				     unsigned int         function_id,
				     unsigned int         inj_ctr_grp_id,
				     unsigned int         inj_ctr_id,
				     unsigned int         send_offset,
				     unsigned int         msg_len
				    );


/*!
 * \brief Set or Change the Hint Bits in a Fifo Descriptor
 *
 * \param[in,out]  desc   Pointer to descriptor to be set or changed.
 * \param[in]      hints  Hint bits to be set.
 *
 * \return None
 *
 */
__INLINE__ void DMA_SetHints(
			     DMA_InjDescriptor_t *desc,
			     unsigned int         hints
			    )
{
  SPI_assert( desc != NULL );
  desc->hwHdr.Hint = hints;

}


/*!
 * \brief Set or Change the Virtual Channel and Dynamic Bit in a Descriptor
 *
 * \param[in,out]  desc  Pointer to descriptor to be set or changed.
 * \param[in]      vc    Input virtual channel
 *                       - 0 = Virtual channel dynamic 0
 *                       - 1 = Virtual channel dynamic 1
 *                       - 2 = Virtual channel deterministic bubble
 *                       - 3 = Virtual channel deterministic priority
 *
 * \return None
 *
 * \post The Dynamic bit is set according to the specified virtual channel.
 *
 */
__INLINE__ void DMA_SetVc(
			  DMA_InjDescriptor_t *desc,
			  unsigned int         vc
			 )
{
  SPI_assert( desc != NULL );

  switch(vc) {
   case DMA_PACKET_VC_D0:
   case DMA_PACKET_VC_D1:
     desc->hwHdr.Dynamic =1;
     break;

   case DMA_PACKET_VC_BN:
   case DMA_PACKET_VC_BP:
    desc->hwHdr.Dynamic =0;
    break;

   default:
     SPI_assert(0);
  }
  desc->hwHdr.VC = vc;

}


/*!
 * \brief Set Descriptor Pid Bits
 *
 * Given a pointer to the descriptor and the receive-side counter group number,
 * set the Pid0 and Pid1 bits in the torus hardware header portion of the
 * descriptor.
 *
 * \param[in]  desc  Pointer to injection descriptor
 * \param[in]  g     Reception-side counter group number
 *                   (0 through DMA_NUM_COUNTER_GROUPS).
 *
 * \return None
 *
 */
__INLINE__ void DMA_SetDescriptorPids(
				      DMA_InjDescriptor_t *desc,
				      unsigned int         g
				     )
{
  /* Set the pid bits according to the group id g */
  desc->hwHdr.Pid0 = _GN(g,30);
  desc->hwHdr.Pid1 = _GN(g,31);
/* ---------------------------------
  The above code performs the following:

  switch(g) {
  case 0:
    desc->hwHdr.Pid0      = 0;
    desc->hwHdr.Pid1      = 0;
    break;

  case 1
    desc->hwHdr.Pid0      = 0;
    desc->hwHdr.Pid1      = 1;
    break;

  case 2
    desc->hwHdr.Pid0      = 1;
    desc->hwHdr.Pid1      = 0;
    break;

  case 3
    desc->hwHdr.Pid0      = 1;
    desc->hwHdr.Pid1      = 1;
    break;

  default:
    SPI_assert(0);

  }
  --------------------------------- */
}


/*!
 * \brief Set or Change the Number of Chunks in a Fifo Descriptor
 *
 * \param[in,out]  desc           Pointer to the descriptor to be set or
 *                                changed.
 * \param[in]      packet_chunks  Number of 32B chunks in the packet
 *                                (1 through 8).
 *
 * \return None
 *
 */
__INLINE__ void DMA_SetChunks(
			      DMA_InjDescriptor_t *desc,
			      unsigned int         packet_chunks
			     )
{
  SPI_assert( desc != NULL );
  SPI_assert( packet_chunks >=1);
  SPI_assert( packet_chunks <=8);
  desc->hwHdr.Chunks = (packet_chunks-1) ;
}


/*!
 * \brief Set or Change the Message Length in a Fifo Descriptor
 *
 * \param[in,out]  desc     Pointer to the descriptor to be set or changed.
 * \param[in]      msg_len  Number of bytes in the payload of the message.
 *
 * \return None
 *
 */
__INLINE__ void DMA_SetMessageLength(
				     DMA_InjDescriptor_t *desc,
				     unsigned int msg_len
				    )
{
  SPI_assert( desc != NULL );

  desc->msg_length= msg_len;
}


/*!
 * \brief Change the Checksum Characteristics in a Fifo Descriptor
 *
 * \param[in,out]  desc       Pointer to the descriptor to be changed.
 * \param[in]      csum_skip  The number of 2-bytes to skip in the checksum
 *                            (7 bits).
 * \param[in]      skip       The checksum skip attribute:
 *                            0 = The packets participates in the injection
 *                                checksum.
 *                            1 = The packet does not participate in the
 *                                injection checksum.
 *
 * \return None
 *
 */
__INLINE__ void DMA_SetInjCsum(
			       DMA_InjDescriptor_t *desc,
			       unsigned int         csum_skip,
			       unsigned int         skip
			      )
{
  SPI_assert( desc != NULL );
  SPI_assert( skip <=1 );

  desc->hwHdr.CSum_Skip = csum_skip;
  desc->hwHdr.Sk        = skip;

}


/*!
 * \brief Determine the Number of Packet Chunks for the First Packet of a
 *        Message
 *
 * Compute the best (largest) packet size in units of 32B chunks given the
 * message length.
 *
 * \param[in]  msg_len  Message length
 *
 * \retval  numPacketChunks  Number of 32B chunks needed in the first packet
 *                           of a message whose length is msg_len.
 *                           This will be a number from 1 through 8.
 * \retval  0                This is considered an error, resulting from a
 *                           msg_len = 0.  The DMA must send at least 1 byte.
 */
__INLINE__ int  DMA_PacketChunks(
				 unsigned msg_len
				)
{
  /* Do most common case first */
  if (msg_len > 208) return 8;

  /* Error case...the DMA must send at least one byte of data */
  SPI_assert( msg_len > 0);

  /* Basically add in the packet header and round to 32B multiple */
  int chunks = ( msg_len - 1 + sizeof(DMA_PacketHeader_t) ) / 32;
  return (1+chunks);

}


/*!
 * \brief Zero Out All Fields a Descriptor
 *
 * \param[in]  desc  Pointer to descriptor to be zero'd.
 *
 * \post The descriptor is zero'd.
 *
 */
__INLINE__ void  DMA_ZeroOutDescriptor(
				       DMA_InjDescriptor_t *desc
				      )
{
  /*
   * Possible optimizations:
   * There are 32 bytes in the descriptor and it should be L1 aligned.
   * SPI_assert(( desc & 0x000000FF) == 0); // check alignment, not needed if can't
   *                                    // easily use double hummer.
   * _bgp_dcache_zero_line(desc);       //Not allowed with SWOA
   * Should be a better way to do this.
   */

  SPI_assert( desc != NULL );

  int *addr = (int *) desc ;

  /* Generates 8 stw's */
  addr[0] = 0;
  addr[1] = 0;
  addr[2] = 0;
  addr[3] = 0;
  addr[4] = 0;
  addr[5] = 0;
  addr[6] = 0;
  addr[7] = 0;

}



/*!
 * \brief Update the Offset and Length in a Descriptor
 *
 * \param[in]  desc    Pointer to descriptor to be updated.
 * \param[in]  offset  The new offset value.
 * \param[in]  length  The new length value.
 *
 * \post The descriptor is updated.
 *
 */
__INLINE__ void DMA_DescriptorUpdateOffsetLength (DMA_InjDescriptor_t *desc,
						  unsigned offset, 
						  unsigned length) 
{
  desc->base_offset = offset;
  desc->msg_length  = length;
}



/*!
 * \brief Set the Put Offset in a Descriptor
 *
 * This sets the "put_offset" field of the software packet header in the
 * provided descriptor.  This field is placed into the packet header by
 * the DMA.  In the first packet, this field is placed into the packet
 * unchanged.  In each subsequent packet, the DMA adds to this field 
 * the number of payload bytes from the previous packet.
 *
 * \param[in]  desc    Pointer to descriptor.
 * \param[in]  offset  The offset value to be set.
 *
 * \post The Put Offset in the descriptor is set.
 *
 */
__INLINE__ void DMA_DescriptorSetPutOffset (DMA_InjDescriptor_t *desc,
					    unsigned offset)
{
  desc->hwHdr.Put_Offset = offset;
}

__END_DECLS

#endif 
