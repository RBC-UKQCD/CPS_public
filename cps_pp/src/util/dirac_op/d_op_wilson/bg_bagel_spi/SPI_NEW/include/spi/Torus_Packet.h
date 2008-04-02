/* begin_generated_IBM_copyright_prolog                             */
/*                                                                  */
/* This is an automatically generated copyright prolog.             */
/* After initializing,  DO NOT MODIFY OR MOVE                       */
/* ---------------------------------------------------------------- */
/* IBM Confidential                                                 */
/*                                                                  */
/* Licensed Machine Code Source Materials                           */
/*                                                                  */
/* Product(s):                                                      */
/*     Blue Gene/P Licensed Machine Code                            */
/*                                                                  */
/* (C) Copyright IBM Corp.  2007, 2007                              */
/*                                                                  */
/* The Source code for this program is not published  or otherwise  */
/* divested of its trade secrets,  irrespective of what has been    */
/* deposited with the U.S. Copyright Office.                        */
/*  --------------------------------------------------------------- */
/*                                                                  */
/* end_generated_IBM_copyright_prolog                               */

#ifndef _TORUS_PACKET_SPI_H_ // Prevent multiple inclusion
#define _TORUS_PACKET_SPI_H_

//
// Torus Packet Interface (non-DMA Mode) interfaces.
//

#include <common/namespace.h>

__BEGIN_DECLS

#include <common/linkage.h>
#include <common/alignment.h>
#include <bpcore/bgp_types.h>
#include <bpcore/bgp_MemoryMap.h>
#include <bpcore/bgp_torus_memmap.h>
#include <bpcore/bgp_torus_packet.h>
#include <bpcore/ppc450_inlines.h>


/*!
 * \brief Torus Fifo Status Data Structure
 *
 *   Each byte contains the number of 32 byte chunks of data in a fifo.
 *   Each Group (0 and 1) of fifos has its own status register.
 *
 * \warning An injection fifo should be considered "full" when its status byte is greater than or equal to 24.
 * \warning This structure when stored in memory must be 16-Byte aligned.
 */
typedef struct Torus_FifoStatus_t
                {
                uint8_t r0,r1,r2,r3,r4,r5, // Normal Priority Reception Fifos
                        rH,                // High   Priority Reception Fifo
                        i0,i1,i2,          // Normal Priority Injection Fifos
                        iH,                // High   Priority Injection Fifo
                        pad[5];            // Pad to 16 bytes
                }
                ALIGN_QUADWORD Torus_FifoStatus_t;

/*!
 * \brief Torus Fifo Data Type
 *
 *   Torus Fifos may only be access via quadword loads and stores.
 */
typedef _QuadWord_t Torus_Fifo_t;

/*!
 * \brief Torus Packet Hardware-Header Constants
 *
 */
// Torus Packet Header Hint Bits (for use in structure bit-fields)
#define TORUS_PKT_HINT_XP  (0x20)   //! Packet can travel in the X+ direction
#define TORUS_PKT_HINT_XM  (0x10)   //! Packet can travel in the X- direction
#define TORUS_PKT_HINT_YP  (0x08)   //! Packet can travel in the Y+ direction
#define TORUS_PKT_HINT_YM  (0x04)   //! Packet can travel in the Y- direction
#define TORUS_PKT_HINT_ZP  (0x02)   //! Packet can travel in the Z+ direction
#define TORUS_PKT_HINT_ZM  (0x01)   //! Packet can travel in the Z- direction

// Torus Packet Header Virtual Circuits
#define TORUS_PKT_VC_D0 (0)         //! Packet enters network on Dynamic 0 Virtual Circuit
#define TORUS_PKT_VC_D1 (1)         //! Packet enters network on Dynamic 1 Virtual Circuit
#define TORUS_PKT_VC_BN (2)         //! Packet enters network on Bubble Normal Virtual Circuit
#define TORUS_PKT_VC_BP (3)         //! Packet enters network on Bubble Priority Virtual Circuit

// Torus Packet Header Dynamic Routing
#define TORUS_PKT_DETERMINSTIC (0)  //! Packet is routed determinstically.
#define TORUS_PKT_DYNAMIC      (1)  //! Packet is routed dynamically.

// Torus Packet Header Routing
#define TORUS_PKT_POINT2POINT  (0)  //! Packet is Point-to-Point.
#define TORUS_PKT_CLASS        (1)  //! Packet is Deposited at all nodes along route.


/*!
 * \brief Torus Packet Hardware-Header Word 0
 *
 */
typedef struct Torus_HdrHW0_t    // first 32bit word of a Torus Hardware Header
  {
  union {
      uint32_t word;                // access fields as 32-bit word
      struct {
           unsigned CSum_Skip : 7;  // Number of shorts (2B) to skip in CheckSum.
           unsigned Sk        : 1;  // 0=use CSum_Skip, 1=Skip entire pkt.

           unsigned Hint      : 6;  // Hint Bits
           unsigned Dp        : 1;  // Deposit Bit for Class Routed MultiCast
           unsigned Pid0      : 1;  // Destination Fifo Group MSb

           unsigned Chunks    : 3;  // Size in Chunks of 32B (0 for 1 chunk, ... , 7 for 8 chunks)
           unsigned Pid1      : 1;  // Destination Fifo Group LSb
           unsigned Dm        : 1;  // 1=DMA Mode, 0=Fifo Mode
           unsigned Dynamic   : 1;  // 1=Dynamic Routing, 0=Deterministic Routing.
           unsigned VC        : 2;  // Virtual Channel (0=D0,1=D1,2=BN,3=BP)

           unsigned X         : 8;  // Destination X Physical Coordinate
          };
        };
  }
  Torus_HdrHW0_t;


/*!
 * \brief Torus Packet Hardware-Header Word 1
 *
 */
typedef struct Torus_HdrHW1_t    // second 32bit word of a Torus Hardware Header
  {
  union {
        uint32_t word;              // access fields as 32-bit word
        struct {
             unsigned Y       : 8;  // Destination Y Physical Coordinate
             unsigned Z       : 8;  // Destination Z Physical Coordinate
             unsigned Resvd0  : 8;  // Reserved (pkt crc)
             unsigned Resvd1  : 8;  // Reserved (pkt crc)
             };
        };
  }
  Torus_HdrHW1_t;


/*!
 * \brief Torus Packet Hardware-Header Word 2
 *
 */
typedef union Torus_HdrHW2_t     // third 32bit word of a Torus Hardware Header (DMA Mode)
  {
  uint32_t word;
  uint32_t Put_Offset;
  }
  Torus_HdrHW2_t;


/*!
 * \brief Torus Packet Hardware-Header Word 3
 *
 */
typedef struct Torus_HdrHW3_t    // fourth 32bit word of a Torus Hardware Header (DMA Mode)
  {
  union {
        uint32_t word;                      // access fields as 32-bit word
        struct {                            // not used in non-DMA mode
             unsigned rDMA_Counter  : 8;
             unsigned Payload_Bytes : 8;
             unsigned Flags         : 8;    // Flags[6]=Pacing, Flags[7]=Remote-Get
             unsigned iDMA_Fifo_ID  : 8;
             };
        };
  }
  Torus_HdrHW3_t;


/*!
 * \brief Torus Packet Hardware-Header
 *
 */
typedef struct Torus_PacketHeader_t
   {
   Torus_Hdr_HW0_t hwh0;      //! \see Torus_Hdr_HW0_t
   Torus_Hdr_HW1_t hwh1;      //! \see Torus_Hdr_HW1_t
   Torus_Hdr_HW2_t hwh2;      //! In non-DMA Mode, this is available for software use.
   Torus_Hdr_HW3_t hwh3;      //! In non-DMA Mode, this is available for software use.
   }
   ALIGN_QUADWORD Torus_PacketHeader_t;


/*!
 * \brief Torus Packet Payload Data Type
 *
 *   Torus Fifos may only be access via quadword loads and stores.
 *   To facilitate rapid data movement to/from fifos, payloads are
 *   stored in memory as a 16-Byte aligned array of quadwords.
 */
typedef _QuadWord_t Torus_PacketPayload_t;


/*!
 * \brief Obtain status for Torus Fifo Group 0
 *   Status is returned at the 'stat' memory address, which must be 16 Byte aligned.
 */
__INLINE__ void Torus_GetStatus0( Torus_FifoStatus_t *stat )
{
   _bgp_QuadMove( _BGP_TORUS_FIFO_STATUS0, stat, 0 );
}


/*!
 * \brief Obtain status for Torus Fifo Group 1
 *   Status is returned at the 'stat' memory address, which must be 16 Byte aligned.
 */
__INLINE__ void Torus_GetStatus1( Torus_FifoStatus_t *stat )
{
   _bgp_QuadMove( _BGP_TORUS_FIFO_STATUS1, stat, 0 );
}


/*
 * \brief Non-DMA Mode Torus Packet Header Creation
 *
 *  Defaults used:  Point-to-Point,
 *                  Normal Priority,
 *                  Dynamic on VCD0,
 *                  Supplied Hint Bits (if 0, hardware calculates them),
 *                  Checksum will skip only Hardware Header.
 *
 *  It is recommended that headers be pre-created for regular and repeated
 *   communication patterns.  This function returns the address of the header
 *   so that it may optionally be used to create headers on-the-fly passed to
 *   the packet injection routines.
 *
 *  \warning Header address must be aligned on 16Byte boundary, i.e. quadword aligned.
 */
__INLINE__ Torus_PacketHeader_t *Torus_MakeHeader(
                     Torus_PacketHeader_t *hdr,    // Filled in on return
                     int dest_x,                   // destination node X dimension physical coordinate
                     int dest_y,                   // destination node Y dimension physical coordinate
                     int dest_z,                   // destination node Z dimension physical coordinate
                     int dest_f,                   // destination Fifo Group: 0 or 1.
                     int hints,                    // hint bits: required to be non-zero if Torus loopback
                     uint32_t sw0,                 // Software Header Word 0
                     uint32_t sw1  )               // Software Header Word 1
{
   // fill out hardware header word 0
   hdr->hwh0.CSum_Skip = 4;                        // skip 8 byte hardware header
   hdr->hwh0.Sk        = 0;                        // use CSum_Skip (ie, don't skip entire packet)
   hdr->hwh0.Hint      = hints;                    // Hint Bits to apply.  May be 0
   hdr->hwh0.Dp        = 0;                        // Point to Point Packet
   hdr->hwh0.Pid0      = dest_f;                   // Destination Fifo Group
   hdr->hwh0.Dm        = 0;                        // Non-DMA, i.e., FIFO Mode.
   hdr->hwh0.Chunks    = 0;                        // init chunks, updated during send
   hdr->hwh0.Pid1      = 0;                        // N/A in Non-DMA mode
   hdr->hwh0.Dynamic   = TORUS_PKT_DYNAMIC;        // use Dynamic Routing
   hdr->hwh0.VC        = TORUS_PKT_VC_D0;          // inject into VC Dynamic 0
   hdr->hwh0.X         = dest_x;                   // Destination X coord

   // fill out hardware header word 1
   hdr->hwh1.Y         = dest_y;                   // Destination Y coord
   hdr->hwh1.Z         = dest_z;                   // Destination Z coord

   // fill out software header word 0
   hdr->hwh2.word      = sw0;                      // Software part of Hardware Header Word 0
   hdr->hwh3.word      = sw1;                      // Software part of Hardware Header Word 1

   return( hdr );
}


//
// Torus FIFO-Mode Packet Injection
//   Caller must have pre-checked space in the injection fifo.
//   Bytes will be rounded up to fill the chunk.
__INLINE__ void Torus_InjectPacket(
                 Torus_Fifo_t          *fifo,     // Injection FIFO to use
                 Torus_PacketHeader_t  *header,   // Packet Header
                 Torus_PacketPayload_t *payload,  // Packet payload
                 int                   bytes    ) // Bytes in packet payload
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
//   Caller must have checked that a packet is available in the reception fifo.
__INLINE__ int Torus_ReceivePacket(
           Torus_Fifo_t          *fifo,     // Injection FIFO to use
           Torus_PacketHeader_t  *header,   // Packet Header received here
           Torus_PacketPayload_t *payload ) // Packet payload received here (rounded up to chunks)
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
                    : "b" (u),
                      "b" (fifo),
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
                    : "b" (u),
                      "b" (fifo),
                      "0" (dst)
                    : "memory",
                      "fr1", "fr2" );
      }

   // there cannot be an odd number of quads
   return( (quads * 16) + 16 );
}

__END_DECLS



#endif // Add nothing below this line
