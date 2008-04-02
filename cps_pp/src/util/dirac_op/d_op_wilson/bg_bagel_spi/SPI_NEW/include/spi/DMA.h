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

#ifndef _DMA_H_ // Prevent multiple inclusion
#define _DMA_H_

/*!
 * \file spi/DMA.h
 *
 * \brief DMA Master Include File for all DMA System Programming Intefaces (SPIs)
 * 
 * The DMA (Direct Memory Access) has an injection and reception engine.
 *
 * The Injection DMA (iDMA) reads data from memory and writes it into the torus 
 * hardware injection fifos.  This activity is defined by a "descriptor".
 * Descriptors are injected into in-memory injection fifos that are processed
 * by the iDMA.  Each descriptor describes a particular memory-to-torus 
 * transfer.  It may involve a small amount of data, or a very large block of
 * data.  When there is space in a torus hardware injection fifo, the
 * iDMA selects the next descriptor in an in-memory injection fifo,
 * manufactures packets, and injects the packets into the torus hardware 
 * injection fifo.  Each in-memory injection fifo is mapped to one or more of
 * the torus hardware injection fifos.  Each compute node may have up to 128
 * in-memory injection fifos, typically 32 per core.
 *
 * A descriptor designates one of two types of transfers: memory-fifo or direct.
 * A memory fifo transfer means that the packets arriving at the destination
 * node are placed into in-memory reception fifos for processing by the core.
 * The core usually copies the packet payload to its final memory destination.
 * A direct transfer means that the packets arriving at the destination are
 * self-describing such that their payload is placed directly into their final
 * memory destination by the DMA.
 * 
 * The Reception DMA (rDMA) reads packets from the torus hardware reception
 * fifos and processes them as memory-fifo or direct transfers depending on the
 * Dm bit setting in the packet header.  For memory-fifo transfers, the rDMA 
 * uses the Pid1 bit in the packet together with the torus hardware reception 
 * fifo Id where the packet arrived to determine which of the (up to 28) 
 * in-memory reception fifos to copy the packet into.  The reception fifos are 
 * polled by the cores, and registered functions are called to process the 
 * packets based on the Func_Id in the packet header.  Direct transfers are 
 * processed by the rDMA.
 *
 * Both the iDMA and rDMA use DMA counters to count the number of payload bytes
 * injected by the iDMA and received by the rDMA.  The cores set up the counters
 * and poll them to determine when an injection or reception is complete.
 * There are 256 iDMA and 256 rDMA counters.  Each counter contains the
 * following:
 * - A base address.  The descriptor and/or packet header have an offset.  
 *   The combination of the base address plus the offset identify the first 
 *   byte of the transfer.
 * - A max address.  This is used by the rDMA as an upper bound on the memory
 *   that it will write into.  Ignored by the iDMA.
 * - A value.  This is typically initialized by the core to the number of bytes
 *   to be transferred.  The DMA decrements it as it transfers the data.
 * - An increment value.  This is used by the core to increment or decrement
 *   the value.
 * There is always an iDMA counter associated with each descriptor that is
 * injected.  There is always an rDMA counter associated with each direct
 * transfer.  There is no rDMA counter associated with a memory-fifo transfer.
 *
 * The DMA has a remote-get capability.  Node A sends a remote-get packet
 * to node B containing one or more descriptors.  The descriptors describe
 * data on node B that is to be sent to node A.  The rDMA on node B receives
 * the packet and injects the descriptors into one of node B's injection
 * fifos.  The data on node B described by the descriptors is received by
 * node A.  From node A's perspective, it did a remote-get of node B's data.
 * Only the DMA on node B is involved (no core on node B is involved).
 *
 * The DMA has local-copy capability where data is transferred between cores 
 * on the same node.  This is accomplished using a local-copy bit in the
 * descriptor, and designating a local-copy in-memory reception fifo.  Local
 * copies may be memory-fifo or direct.
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


#include <spi/DMA_Assert.h>       /* DMA Assert definitions         */
#include <spi/DMA_Packet.h>       /* DMA Packet definitions         */
#include <spi/DMA_Counter.h>      /* DMA Counter definitions        */
#include <spi/DMA_Fifo.h>         /* DMA Common Fifo definitions    */
#include <spi/DMA_InjFifo.h>      /* DMA Injection Fifo definitions */
#include <spi/DMA_RecFifo.h>      /* DMA Reception Fifo definitions */
#include <spi/DMA_Descriptors.h>  /* DMA Descriptor definitions     */


__END_DECLS


#endif 
