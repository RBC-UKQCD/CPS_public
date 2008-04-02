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

#ifndef	_DMA_PACKET_H_ /* Prevent multiple inclusion */
#define	_DMA_PACKET_H_


/*!
 * \file spi/DMA_Packet.h
 *
 * \brief DMA SPI Packet Definitions
 *
 */


#include <common/namespace.h>



__BEGIN_DECLS


/*!
 * \brief Hint Bit: Packet wants to travel in the X plus direction.
 */
#define DMA_PACKET_HINT_XP  (0x20)


/*!
 * \brief Hint Bit: Packet wants to travel in the X minus direction.
 */
#define DMA_PACKET_HINT_XM  (0x10)


/*!
 * \brief Hint Bit: Packet wants to travel in the Y plus direction.
 */
#define DMA_PACKET_HINT_YP  (0x08)


/*!
 * \brief Hint Bit: Packet wants to travel in the Y minus direction.
 */
#define DMA_PACKET_HINT_YM  (0x04)


/*!
 * \brief Hint Bit: Packet wants to travel in the Z plus direction.
 */
#define DMA_PACKET_HINT_ZP  (0x02)


/*!
 * \brief Hint Bit: Packet wants to travel in the Z minus direction.
 */
#define DMA_PACKET_HINT_ZM  (0x01)


/*!
 * \brief Virtual Channel Bits:  Dynamic 0.
 */
#define DMA_PACKET_VC_D0 (0)


/*!
 * \brief Virtual Channel Bits:  Dynamic 1.
 */
#define DMA_PACKET_VC_D1 (1)


/*!
 * \brief Virtual Channel Bits:  Deterministic Bubble Normal.
 */
#define DMA_PACKET_VC_BN (2)


/*!
 * \brief Virtual Channel Bits:  Deterministic Bubble Priority.
 */
#define DMA_PACKET_VC_BP (3)


/*!
 * \brief Dynamic Routing Bit:  Follows deterministic Routing.
 */
#define DMA_PACKET_DETERMINSTIC (0)


/*!
 * \brief Dynamic Routing Bit:  Follows dynamic Routing.
 */
#define DMA_PACKET_DYNAMIC      (1)


/*!
 * \brief Torus Hardware Packet Header Constants for Routing: Point to Point.
 */
#define DMA_PACKET_POINT2POINT  (0)


/*!
 * \brief Torus Hardware Packet Header Constants for Routing: Class.
 */
#define DMA_PACKET_CLASS        (1)


/*!
 * \brief Torus DMA Hardware Packet Header.
 *
 * There are two sections of the packet header:  The hardware header
 * and the software header.
 *
 * The same 8-byte hardware header as was used on Blue Gene/L is used
 * for Blue Gene/P, except that two bits that were previously unused
 * will be used as follows:
 *
 * - The Pid bit on Blue Gene/L indicates the logical destination group.
 *   This determines the reception fifo group a packet ends up in.  
 *   This bit is now called Pid0.  The new Pid1 bit expands the logical 
 *   destination group from two to four.  This corresponds to the increase 
 *   in cores from two to four.
 *
 * - The new Dm bit indicates the DMA mode:  Memory fifo or direct.
 *   In memory fifo mode, the DMA receives packets from the torus fifos into
 *   reception fifos located in memory.  Then the core copies the payload 
 *   from the memory fifo to its final destination.  In direct mode, the DMA 
 *   moves the packet payload directly from the torus fifos to its final 
 *   destination.
 *
 * The 8-byte software header was used by the software on Blue Gene/L for
 * its own purposes.  On Blue Gene/P, parts of it are used by the DMA, 
 * depending on the type of DMA transfer being used.  The usage of the fields
 * in the software header is as follows for the typical types of DMA transfers:
 *
 * - In memory fifo mode, 
 *   - The first 4 bytes of the software header contain the "put offset".  
 *     This is the offset from the injection counter's base address, in bytes, 
 *     of the memory being transferred in this packet.
 *   - The last 4 bytes of the software header is for use by software.
 *
 * - In direct put mode, 
 *   - The first 4 bytes of the software header contain the "put offset".  
 *     This is the offset from the reception counter's base address, in bytes, 
 *     of the memory where the payload in this packet is placed.  
 *   - The fifth byte of the software header is the reception counter Id.  
 *   - The sixth byte of the software header is the number of valid bytes of
 *     payload in this packet.
 *   - The seventh byte of the software header contains DMA flags. Specifically,
 *     the remote-get flag is 0.
 *   - The last byte of the software header is for use by software.
 *
 * - In remote get mode, the payload contains one or more injection descriptors
 *   describing data to be transferred by the DMA.  When the DMA receives this
 *   packet, it injects the descriptors into injection fifos to perform the
 *   specified data transfer.
 *   - The first 5 bytes of the software header are for use by software.
 *   - The sixth byte of the software header is the number of valid bytes of
 *     payload in this packet.  This will be a multiple of 32, since the payload
 *     consists of one or more 32 byte DMA descriptors.
 *   - The seventh byte of the software header contains DMA flags. Specifically,
 *     the remote-get flag is 1.
 *   - The eighth byte of the software header is the injection fifo Id where
 *     the descriptors in the payload will be injected.
 *
 */
typedef struct DMA_PacketHeader_t
{
  union {
      unsigned word0;             /*!< First 4 bytes of packet header.        */
    
      struct {
	unsigned CSum_Skip : 7;   /*!< Number of 2 byte units to skip from
				       the top of a packet before including
                                       the packet bytes into the running
                                       checksum of the torus injection fifo
                                       where this packet is injected.
				  */

	  unsigned Sk        : 1; /*!< Torus injection checksum skip packet
	                               bit.
	                               - 0 includes the packet (excluding the
                                         portion designated by DMA_CSUM_SKIP)
                                         in the checksum.
                                       - 1 excludes the entire packet from
                                         the checksum.
				  */

	  unsigned Hint      : 6; /*!< Hint bits for torus routing (6 bits).
	                               Each bit corresponds to x+, x-, y+, y-,
	                               z+, z-.  If a bit is set, it indicates
                                       that the packet wants to travel along
                                       the corresponding direction.  If all
                                       bits are zero, the hardware calculates
                                       the hint bits.  Both x+ and x- cannot
                                       be set at the same time...same with y
                                       and z.
				  */

	  unsigned Dp        : 1; /*!< Deposit Bit for Class Routed MultiCast.
	                               If this bit is set to 1, then as the
                                       packet travels along a straight line
                                       to its destination, it also deposits
	                               a copy of itself into each node as it
                                       goes through.  This feature must be
                                       used only if the packet is set to
	                               travel along a straight line.
				  */

	  unsigned Pid0      : 1; /*!< Destination Fifo Group Most Significant
                                       Bit.  (Pid0,Pid1) specifies which of 4
                                       reception fifo groups that this packet
                                       is destined for.
				  */

	  unsigned Chunks    : 3; /*!< Size in Chunks of 32B (0 for 1 chunk,
                                       ... , 7 for 8 chunks).
				  */				       

          unsigned Pid1      : 1; /*!< Destination Fifo Group Least
                                       significant bit.  Refer to Pid0.
				  */

  	  unsigned Dm        : 1; /*!< 1=DMA Mode, 0=Fifo Mode.               */

	  unsigned Dynamic   : 1; /*!< 1=Dynamic Routing,
				       0=Deterministic Routing.
				  */

	  unsigned VC        : 2; /*!< Virtual Channel
        	                       - 0=Dynamic 0
                                       - 1=Dynamic 1
                                       - 2=Deterministic Bubble Normal
                                       - 3=Deterministic Bubble Priority
				  */

	  unsigned X         : 8; /*!< Destination X Physical Coordinate.     */

      }; /* End of individual fields in Word 0 */

  }; /* End of Word 0 */


  union {

        unsigned word1;           /*!< Second 4 bytes of packet header.       */

        struct {
	  unsigned Y         : 8; /*!< Destination Y Physical Coordinate.     */

	  unsigned Z         : 8; /*!< Destination Z Physical Coordinate.     */

	  unsigned Resvd0    : 8; /*!< Reserved (pkt crc).                    */

	  unsigned Resvd1    : 8; /*!< Reserved (pkt crc).                    */

        }; /* End of individual fields in Word 1 */

  }; /* End of Word 1 */


  union {

    unsigned word2;               /*!< Third 4 bytes of packet header.        */

    unsigned Put_Offset;          /*!< For a memory fifo packet, gives a
                                       unique ID to each packet in a long
                                       message.  Derived from the put offset
                                       of the torus packet header in the
                                       descriptor, and updated for each
                                       packet.
                                       For a direct-put packet, the rDMA
                                       writes the first payload byte to this
                                       offset plus the base address
                                       corresponding to the rDMA counter ID.
				  */

    unsigned Single_Packet_Parameter; /*!< For a single memory fifo packet,
					   this is essentially unused space
					   that can be used to pass a 
					   parameter to the target node.
				      */
  }; /* End of Word 2 */


  union {

      unsigned word3;             /*!< Fourth 4 bytes of packet header.       */
      
      struct {
	unsigned rDMA_Counter  : 8; /*!< For a direct-put packet, this is the
                                         number of the rDMA counter associated
	                                 with this packet.
				    */

	unsigned Payload_Bytes : 8; /*!< For a direct-put packet, this is the
                                         number of valid bytes in the payload.
                                         This is set by the iDMA, based on the
                 	                 message length in the injection
                                         descriptor.
				    */

	unsigned Flags         : 8; /*!< Flags[6]=Pacing, Flags[7]=Remote-Get.*/

	unsigned iDMA_Fifo_ID  : 8; /*!< For a remote-get packet, this is the
	                                 iDMA fifo ID to be injected during
	                                 remote-get processing.
				    */
      };
  
      struct {                      /*   For memory fifo packets...           */

	unsigned SW_Arg   : 24;     /*!< User-defined.                        */

	unsigned Func_Id  : 8 ;     /*!< Function ID for dispatching receiver
                                         functions from Polling reception
	                                 fifos.
				    */
      };

  }; /* End of Word 3 */

}
ALIGN_QUADWORD DMA_PacketHeader_t;




__END_DECLS


#endif 
