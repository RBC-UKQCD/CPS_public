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
#ifndef _BGP_TORUS_SPI_H_ // Prevent multiple inclusion
#define _BGP_TORUS_SPI_H_



//
// Torus inlines and interfaces available in user space and in kernel.
//

#include <common/namespace.h>

__BEGIN_DECLS

#include <common/linkage.h>
#include <common/alignment.h>
#include <bpcore/bgp_types.h>
#include <bpcore/bgp_torus_packet.h>
#include <bpcore/ppc450_inlines.h>

// Virtual Addresses seen by software to access the status registers
#define _BGP_TORUS_FIFO_STATUS0 (_BGP_VA_TORUS0 | _BGP_MEM_TORUS_STATUS0_OFFSET)
#define _BGP_TORUS_FIFO_STATUS1 (_BGP_VA_TORUS1 | _BGP_MEM_TORUS_STATUS0_OFFSET)
#define _BGP_MEM_TORUS_STATUS0_OFFSET           0x7000 // offset within group of *this* group's status
#define _BGP_MEM_TORUS_STATUS1_OFFSET           0x7010 // offset within group of *other* group's status


//
// Torus Fifo Status:
//   Each byte contains the number of 32 byte chunks of data in a fifo.
//   Each Group (0 and 1) of fifos has its own status register.
typedef struct T_BGP_TorusFifoStatus
                {
                uint8_t r0,r1,r2,r3,r4,r5, // Normal Priority Reception Fifos
                        rH,                // High   Priority Reception Fifo
                        i0,i1,i2,          // Normal Priority Injection Fifos
                        iH,                // High   Priority Injection Fifo
                        pad[5];            // Pad to 16 bytes
                }
                ALIGN_QUADWORD _BGP_TorusFifoStatus;

typedef _QuadWord_t _BGP_TorusFifo_t;
typedef _QuadWord_t _BGP_TorusPacketPayload_t;

// Obtain status for Torus Fifo Group 0.
//   Status is returned at the 'stat' memory address.
extern inline void _bgp_TorusGetStatus0( _BGP_TorusFifoStatus *stat )
{
   _bgp_QuadMove( _BGP_TORUS_FIFO_STATUS0, stat, 0 );
}

// Obtain status for Torus Fifo Group 1
//   Status is returned at the 'stat' memory address.
extern inline void _bgp_TorusGetStatus1( _BGP_TorusFifoStatus *stat )
{
   _bgp_QuadMove( _BGP_TORUS_FIFO_STATUS1, stat, 0 );
}


// Non-DMA Mode Torus Packet Header Creation
//  Defaults used:  Point-to-Point,
//                  Normal Priority,
//                  Dynamic on VCD0,
//                  Supplied Hint Bits (if 0, hardware calculates them),
//                  Checksum will skip only Hardware Header.
//  Hdr address must be aligned on 16Byte boundary, i.e. quadword aligned.
extern inline _BGP_TorusPacketHeader_t *_bgp_TorusMakeHeader(
                     _BGP_TorusPacketHeader_t *hdr, // Filled in on return
                     int dest_x,                    // destination coordinates
                     int dest_y,
                     int dest_z,
                     int dest_f,                    // destination Fifo Group
                     int hints,                     // hint bits (required if TS loopback
                     uint32_t sw0,                  // Software Header Word 0
                     uint32_t sw1  )                // Software Header Word 1
{
   // fill out hardware header word 0
   hdr->hwh0.CSum_Skip = 4;                            // skip 8 byte hardware header
   hdr->hwh0.Sk        = 0;                            // use CSum_Skip (ie, don't skip entire packet)
   hdr->hwh0.Hint      = hints;
   hdr->hwh0.Dp        = 0;                            // Point to Point Packet
   hdr->hwh0.Pid0      = dest_f;                       // Destination Fifo Group
   hdr->hwh0.Dm        = 0;                            // Non-DMA, i.e., FIFO Mode.
   hdr->hwh0.Chunks    = 0;                            // init chunks, updated during send
   hdr->hwh0.Pid1      = 0;                            // N/A in Non-DMA mode
   hdr->hwh0.Dynamic   = _BGP_TORUS_PKT_DYNAMIC;       // use Dynamic Routing
   hdr->hwh0.VC        = _BGP_TORUS_PKT_VC_D0;         // inject into VC Dynamic 0
   hdr->hwh0.X         = dest_x;                       // Destination X coord

   // fill out hardware header word 1
   hdr->hwh1.Y         = dest_y;                       // Destination Y coord
   hdr->hwh1.Z         = dest_z;                       // Destination Z coord

   // fill out software header word 0
   hdr->hwh2.word      = sw0;                          // Software part of Hardware Header
   hdr->hwh3.word      = sw1;                          // Software part of Hardware Header

   return( hdr );
}


// Torus FIFO-Mode Packet Injection
extern inline void _bgp_TorusInjectPacket(
           _BGP_TorusFifo_t          *fifo,     // Injection FIFO to use (Must have space available!)
           _BGP_TorusPacketHeader_t  *header,   // Packet Header
           _BGP_TorusPacketPayload_t *payload,  // Packet payload
           int                       bytes    ) // Bytes in packet payload
{
   int chunks = ((bytes + 15) >> 5);
   uint8_t *pdata = (uint8_t *)payload;

   // update chunks in the header
   header->hwh0.Chunks = chunks;

   // compute total quads to store in addition to the header
   int quads = ((chunks << 1) + 1);

   // load header
   _bgp_QuadLoad( header, 0 );

   register int u = 16;
   int n4 = (quads >> 2);

   // store header
   _bgp_QuadStore( fifo, 0 );

   while( n4-- )
      {
      asm volatile( "lfpdx   1,0,%0;"
                    "lfpdux  2,%0,%1;"
                    "lfpdux  3,%0,%1;"
                    "lfpdux  4,%0,%1;"
                    "stfpdx  1,0,%2;"
                    "stfpdx  2,0,%2;"
                    "stfpdx  3,0,%2;"
                    "stfpdx  4,0,%2;"
                    : "+b" (pdata)
                    : "b"  (u),
                      "b"  (fifo),
                      "0"  (pdata)
                    : "memory",
                      "fr1", "fr2", "fr3", "fr4" );

      pdata += 16;
      }

   if ( quads & 2 )
      {
      asm volatile( "lfpdx   1,0,%0;"
                    "lfpdux  2,%0,%1;"
                    "stfpdx  1,0,%2;"
                    "stfpdx  2,0,%2;"
                    : "+b" (pdata)
                    : "b"  (u),
                      "b"  (fifo),
                      "0"  (pdata)
                    : "memory",
                      "fr1", "fr2" );
      pdata += 16;
      }

   if ( quads & 1 )
      {
      _bgp_QuadMove( pdata, fifo, 1 );
      }
}


// Torus FIFO-mode Packet Reception
//   Returns bytes in payload received.
extern inline int _bgp_TorusReceivePacket(
           _BGP_TorusFifo_t          *fifo,     // Injection FIFO to use (Must have space available!)
           _BGP_TorusPacketHeader_t  *header,   // Packet Header received here
           _BGP_TorusPacketPayload_t *payload ) // Packet payload received here (rounded up to chunks)
{
   int chunks, quads;
   uint8_t *dst = (uint8_t *)payload;

   // read first chunk, includes header and 1st quad of payload
   asm volatile( "lfpdx   1,0,%0;"
                 "lfpdx   2,0,%0;"
                 "stfpdx  1,0,%1;"
                 "stfpdx  2,0,%2;"
                 :
                 : "b" (fifo),
                   "b" (header),
                   "b" (payload)
                 : "memory",
                   "fr1", "fr2" );

   chunks = header->hwh0.Chunks; // we already read 0th chunk, so this is chunks left in fifo

   quads  = (chunks << 1); // 2 quads per chunk

   register int u = 16;
   int n4 = (quads >> 2);

   while( n4-- )
      {
      asm volatile( "lfpdx   1,0,%2;"
                    "lfpdx   2,0,%2;"
                    "lfpdx   3,0,%2;"
                    "lfpdx   4,0,%2;"
                    "stfpdux 1,%0,%1;" /* this pre-inc is correct! */
                    "stfpdux 2,%0,%1;"
                    "stfpdux 3,%0,%1;"
                    "stfpdux 4,%0,%1;"
                    : "+b" (dst)
                    : "r" (u),
                      "r" (fifo),
                      "0" (dst)
                    : "memory",
                      "fr1", "fr2", "fr3", "fr4" );
      }

   if ( quads & 2 )
      {
      asm volatile( "lfpdx   1,0,%2;"
                    "lfpdx   2,0,%2;"
                    "stfpdux 1,%0,%1;"
                    "stfpdux 2,%0,%1;"
                    : "+b" (dst)
                    : "r" (u),
                      "r" (fifo),
                      "0" (dst)
                    : "memory",
                      "fr1", "fr2" );
      }

   // there cannot be an odd number of quads
   return( (quads * 16) + 16 );
}

__END_DECLS



#endif // Add nothing below this line
