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

/*! \file spi/bgp_collective_spi.h
 *  \brief BG Collective SPI
 */
#ifndef	_BGP_COLLECTIVE_SPI_H_ 
#define	_BGP_COLLECTIVE_SPI_H_

#include <common/namespace.h>

__BEGIN_DECLS

#include <common/linkage.h>
#include <common/alignment.h>
#include <bpcore/bgp_types.h>
#include <bpcore/bgp_collective.h>
#include <bpcore/bgp_collective_inlines.h>
#include <sys/types.h>

typedef _BGP_TreeHwHdr CollHwHdr;
/*!
  Introduction to the Collective Network System Programming Interface        
  --------------------------------------------------------------------------- 
  All tree packets consist of a 32-bit header and a 256-byte payload. 
  Point-to-point (P2P) tree packets, which are primarily used for communication 
  between compute and I/O nodes, contain a hardware header, and software header 
  at the beginning of the payload.  Collective tree packets contain just the 
  hardware header and a payload. 
  --------------------------------------------------------------------------- 

  --------------------------------------------------------------------------- 
  Diagnostic Injection Checksum Routines                                 
  --------------------------------------------------------------------------- 
  These routines assume that partial packets are always rounded to quadword 
  granularity and padded with zeros. 
*/

/*! Define to cause CollectiveSend* to maintain the software injection checksums */
/* As per Issue 2804, only do this when running as a cnk diagnostic test */
#ifdef __CNK__
#define _SPI_TREE_DO_INJ_CSUM
#endif

/*@{*/
/*! Use these when calling CollectiveInjectionCsumGet: */
#define _BGP_TREE_CSUM_VC0_HEADER   (0)
#define _BGP_TREE_CSUM_VC0_PAYLOAD  (1)
#define _BGP_TREE_CSUM_VC1_HEADER   (2)
#define _BGP_TREE_CSUM_VC1_PAYLOAD  (3)
/*@}*/

#define CollectiveConfigureClass     _bgp_Collective_ConfigureClass
#define CollectiveMakeCollectiveHdr  _bgp_TreeMakeCollectiveHdr
#define CollectiveMakePtpHdr         _bgp_TreeMakePtpHdr
#define CollectiveMakeSendHdr        _bgp_TreeMakeSendHdr
#define CollectiveMakeBcastHdr       _bgp_TreeMakeBcastHdr
#define CollectiveMakeAluHdr         _bgp_TreeMakeAluHdr
#define CollectiveGetStatus          _bgp_TreeGetStatus
#define CollectiveInvertClassMsb     _bgp_TreeInvertClassMsb
#define CollectiveWrapChannel        _bgp_TreeWrapChannel
#define CollectiveUnwrapChannel      _bgp_TreeUnwrapChannel
#define CollectiveRawSendHeader      _bgp_TreeRawSendHeader
#define CollectiveRawReceiveHeader   _bgp_TreeRawReceiveHeader
#define CollectiveRawSendPacket0     _bgp_TreeRawSendPacket0
#define CollectiveRawSendPacket      _bgp_TreeRawSendPacket
#define CollectiveRawReceivePacket   _bgp_TreeRawReceivePacket
#define CollectiveRawReceivePacketNoHdr   _bgp_TreeRawReceivePacketNoHdr
#define CollectiveRawReceivePacketNoHdrNoStore   _bgp_TreeRawReceivePacketNoHdrNoStore

/*@{*/
/*!
 * \brief Collective Reception Fifo Receive Function Prototype
 *
 * Receive functions must be registered through the 
 * Collective_RecFifoRegisterRecvFunction interface, which assigns a registration ID
 * to the function.  When the polling functions process a packet in a 
 * reception fifo, the appropriate receive function for that packet is
 * called with a pointer to the packet header, pointer to the payload, and
 * length of the payload.  
 *
 * \param[in]  vc              Virtual Channel packet has arrived on, For user space code must use vc=0
 * \param[in]  recv_func_parm  Pointer to storage specific to this receive
 *                             function.  This pointer was specified when the 
 *                             receive function was registered with the kernel,
 *                             and is passed to the receive function
 *                             unchanged.
 * \param[in]  hdr_ptr         Pointer to the first bytes of payload
 * \param[in]  hdr_bytes       Number of bytes in the first bytes of payload
 * \param[in]  payload_bytes   Number of bytes to read from the tree
 *
 * \retval   0                 No errors found while processing the payload.
 * \retval   negative_number   Errors found while processing the payload.
 */
typedef int  (*CollectiveRecvP2PFunction_t)(int   vc,
					    void  *recv_func_parm,
					    void  *hdr_ptr,
					    int    hdr_bytes,
					    int    payload_bytes
					    );

typedef int  (*CollectiveRecvCollFunction_t)(int   vc,
					     void  *recv_func_parm
					     );
/*@}*/	

/*@{*/		
/*!
 * \brief Register a Reception Fifo Receive Function
 *
 * Register a specified receive function to handle packets having a specific
 * "registration ID".  It returns a registration ID that is to be used 
 * For collective packets, the packet header has unused tag bits that can be used
 * for dispatch
 * For P2P packets, the bits are not available, so the first quad is used for dispatch.
 *
 * \param[in]  recv_func          Pointer to the receive function.
 * \param[in]  recv_func_parm     Arbitrary pointer to be passed to the
 *                                recv_func when it is called.
 * \param[in]  payload_bytes(p2p only) How much to deliver in payload, the rest is info
 * \retval   0            This is a registration ID if is_error_function=0 and
 *                        is_header_fifo=0.  Otherwise, it indicates success.
 *           1-255        This is a registration ID.  Successful.
 *           negative     Failure.  This is a negative errno value.
 */
int CollectiveRegisterRecvP2PFunction(CollectiveRecvP2PFunction_t  recv_func,
				      void                        *recv_func_parm,
				      int                          payload_bytes
				      );

int CollectiveRegisterRecvCollFunction(CollectiveRecvCollFunction_t  recv_func,
				       void                         *recv_func_parm
				       );
/*@}*/


/*! Initialize header and payload injection checksums. */
void CollectiveInjectionCsumInit(int      vc, 
				 uint32_t hdr_init, 
				 uint32_t pyld_init);

/*! Add header and payload to injection checksums. */
void CollectiveInjectionCsumAdd(int       vc, 
				uint32_t  hdr, 
				uint16_t *pyld, 
				int       bytes);

/*! Return a checksum. */
uint32_t CollectiveInjectionCsumGet(int offset);

/*! Print all calculated checksums. */
void CollectiveInjectionCsumPrint(char *preprint, 
				  int   print_header);


/*---------------------------------------------------------------------------*
 *     Low-Level Collective Network Packet Send Routines                     *
 *---------------------------------------------------------------------------*/

/*! The following routines should only be called when there is fifo room to send. 
  That is, they will inject unconditionally and over-run if necessary. 
  Note: For on-the-fly header creation, call one of the above header creation 
  routines as a parameter to this function. 
  Use for aligned payload of 16 quadwords. 
  User space code must use virtual channel 0 (vc=0)
*/
/*! Use for aligned payload of 256 byte packets. */


/*! Use for aligned payload of 1-16 quadwords. */
void CollectiveRawSendPacketPartial(int            vc,       /*! Virtual channel (0 or 1) */
				    CollHwHdr      *hdrHW,   /*! Previously created hardware header (any type) */
				    void           *pyld,    /*! Source address of payload (16-byte aligned) */
				    int            bytes );  /*! Payload bytes (multiple of 16, up to 256) */

/*! Use for unaligned payload or bytes = 1-256. */
int CollectiveRawSendPacketUnaligned(int            vc,       /*! Virtual channel (0 or 1) */
				     CollHwHdr      *hdrHW,   /*! Previously created hardware header (any type) */
				     void           *pyld,    /*! Source address of payload */
				     int            bytes );  /*! Payload bytes (1-256) */

/*! Use for aligned payload of 1-16 quadwords. */
int CollectiveRawSendPayloadPartial(int            vc,       /*! Virtual channel (0 or 1) */
				    void           *pyld,    /*! Source address of payload (16-byte aligned) */
				    int            bytes );  /*! Payload bytes (multiple of 16, up to 256) */

/*! Use for unaligned payload or bytes = 1-256. */
int CollectiveRawSendPayloadUnaligned(int            vc,       /*! Virtual channel (0 or 1) */
				      void           *pyld,    /*! Source address of payload */
				      int            bytes );  /*! Payload bytes (1-256) */

/*! Inject a single header */

/*! Use for unaligned payload or bytes = 1-256. This function will not fill out a payload, but will fill to quadword. */
int CollectiveRawInjectBytesUnaligned(int            vc,       /*! Virtual channel (0 or 1) */
				      void           *pyld,    /*! Source address of payload */
				      int            bytes );  /*! Payload bytes (1-256) */

/*! The following routines will block indefinitely until there is room in the specified
  injection fifo. This implies that they poll the status register.
  !! WARNING: If the payload of a packet send routine is greater than a single packet's 
  !! worth of data, then every packet will use the specified header. 
  Use for aligned payload that is a multiple of 16 bytes in length. 
*/
int CollectiveSend(int            vc,       /*! Virtual channel (0 or 1) */
		   CollHwHdr      *hdrHW,   /*! Previously created hardware header (any type) */
		   void           *pyld,    /*! Source address of payload (16-byte aligned) */
		   int            bytes );  /*! Payload bytes (multiple of 16) */

/*! Use for any size payload at any alignment. */
int CollectiveSendUnaligned(int            vc,       /*! Virtual channel (0 or 1) */
			    CollHwHdr      *hdrHW,   /*! Previously created hardware header (any type) */
			    void           *pyld,    /*! Source address of payload */
			    int            bytes );  /*! Payload bytes */


/*---------------------------------------------------------------------------*
 *     Low-Level Collective Network Packet Receive Routines                  *
 *---------------------------------------------------------------------------*/

/*! The following routines should only be called when there is something to receive. 
    That is, they will receive unconditionally and under-run if necessary. 
    Use for aligned payload of 16 quadwords. 
*/

/*! Use for aligned payload of 1-16 quadwords. */
void CollectiveRawReceivePacketPartial(int            vc,       /*! Virtual channel (0 or 1) */
				       CollHwHdr      *hdrHW,   /*! Hardware header buffer */
				       void           *pyld,    /*! Payload buffer (16-byte aligned) */
				       int            bytes );  /*! Payload bytes (multiple of 16, up to 256) */

/*! Use for unaligned payload or bytes = 1-256. */
int CollectiveRawReceivePacketUnaligned(int            vc,       /*! Virtual channel (0 or 1) */
					CollHwHdr      *hdrHW,   /*! Hardware header buffer */
					void           *pyld,    /*! Payload buffer */
					int            bytes );  /*! Payload bytes (1-256) */

/*! Use for aligned payload of 1-16 quadwords. */
int CollectiveRawReceivePayloadPartial(int            vc,       /*! Virtual channel (0 or 1) */
				       void           *pyld,    /*! Payload buffer (16-byte aligned) */
				       int            bytes );  /*! Payload bytes (multiple of 16, up to 256) */

/*! Use for unaligned payload or bytes = 1-256. */
int CollectiveRawReceivePayloadUnaligned(int            vc,       /*! Virtual channel (0 or 1) */
					 void           *pyld,    /*! Payload buffer */
					 int            bytes );  /*! Payload bytes (1-256) */

/*! Receive a single header */

/*! Use for unaligned payload or bytes = 1-256. This function will not dump the remainder of a payload, but will 
  dump the remainder of a quadword. 
*/
int CollectiveRawReceiveBytesUnaligned(int            vc,       /*! Virtual channel (0 or 1) */
				       void           *pyld,    /*! Source address of payload */
				       int            bytes );  /*! Payload bytes (1-256) */

/*! The following routines will block indefinitely until there is data in the specified 
  receive fifo. This implies that they poll the status register. 
  Use for aligned payload that is a multiple of 16 bytes in length. 
*/
int CollectiveReceive(int            vc,       /*! Virtual channel (0 or 1) */
		      CollHwHdr      *hdrHW,   /*! Hardware header buffer */
		      void           *pyld,    /*! Payload buffer (16-byte aligned) */
		      int            bytes );  /*! Payload bytes (multiple of 16) */

/*! Use for any size payload at any alignment. */
int CollectiveReceiveUnaligned(int            vc,       /*! Virtual channel (0 or 1) */
			       CollHwHdr      *hdrHW,   /*! Hardware header buffer */
			       void           *pyld,    /*! Payload buffer */
			       int            bytes );  /*! Payload bytes */


extern inline uint32_t _issue3193_In32( uint32_t *vaddr )
{
	volatile uint32_t *va = (volatile uint32_t *)vaddr;

	uint32_t v = *va;
	/*
	 * Warning! Both 'mbar' and 'isync' are required, in that order.
	 * They also must be place *after* the load instruction.
	 * Otherwise there is a timing (instruction path length) issue
	 * whereby the results of the "In32" may not be present in the
	 * destination when it is first used.
	 *
	 * Note, 'msync' will not fix the problem.
	 */
	_bgp_mbar();
	_bgp_isync();

	return( v );
}

static inline int CollectiveFifoStatus(int vc, 
			 unsigned *rechcount,
			 unsigned *recdcount,
			 unsigned *injhcount,
			 unsigned *injdcount)
{
  _BGP_TreeFifoStatus stat;
  if ( vc )
      stat.status_word = _issue3193_In32( (uint32_t *)_BGP_TR1_S1 );
  else
      stat.status_word = _issue3193_In32( (uint32_t *)_BGP_TR0_S0 );

  *rechcount = stat.RecHdrCount;
  *recdcount = stat.RecPyldCount;

  *injhcount = stat.InjHdrCount;
  *injdcount = stat.InjPyldCount;
  return 0;
}

typedef struct P2PTblEntry_t
{
  CollectiveRecvP2PFunction_t   P2PrecvFunction;
  void                         *recv_func_parm;
  int                           payload_bytes;
} P2PTblEntry_t;

typedef struct CollTblEntry_t
{
  CollectiveRecvCollFunction_t  CollrecvFunction;
  void                         *recv_func_parm;
} CollTblEntry_t;

extern P2PTblEntry_t  P2PrecvFunctions[256];
extern CollTblEntry_t CollrecvFunctions[256];

static inline int CollectiveRecFifoPoll(int vc, 
					unsigned rechcount, 
					unsigned recdcount)
{
  int packets = 0;
  int done    = 0;
  volatile _BGP_TreeHwHdr hdr;
  void                   *hdrptr;

  if(vc == 0)
    hdrptr = (void *)_BGP_TR0_HR;
  else
    hdrptr = (void *)_BGP_TR1_HR;
  

  while (rechcount>0 && recdcount>0 && !done)
    {
      /* Read Packet Header */
      hdr.CtvHdr.word = _bgp_In32((uint32_t*)hdrptr);
      if(hdr.CtvHdr.Ptp == 0)
	{
	  /* Collective Packet Dispatch */
	  unsigned protonum = hdr.CtvHdr.Tag;
	  CollectiveRecvCollFunction_t       recv_func = 
	    CollrecvFunctions[protonum].CollrecvFunction;
	  void                         *recv_func_parm = 
	    CollrecvFunctions[protonum].recv_func_parm;
	  
	  done = recv_func(vc, recv_func_parm);
	  rechcount -=1;
	  recdcount -=16;
	  packets++;
	}
      else
	{
	  /* P2P Packet Dispatch */
	  _BGP_TreePayload pkt;
	  CollectiveRawReceivePayloadPartial(vc,&pkt,16);
	  unsigned protonum = pkt.u32[0];
	  CollectiveRecvP2PFunction_t        recv_func = 
	    P2PrecvFunctions[protonum].P2PrecvFunction;
	  void                         *recv_func_parm = 
	    P2PrecvFunctions[protonum].recv_func_parm;
	  int                           payload_bytes  =
	    P2PrecvFunctions[protonum].payload_bytes;

	  int hdr_bytes = 256 - payload_bytes;
	  if(hdr_bytes != 16)
	    CollectiveRawReceivePayloadPartial(vc,
					       &pkt.q[1],
					       hdr_bytes-16);
	  done = recv_func(vc,
			   recv_func_parm, 
			   (void*)&pkt,
			   hdr_bytes,
			   payload_bytes);
	  rechcount -=1;
	  recdcount -=16;
	  packets++;
	}
    }
  return packets;
}




__END_DECLS



#endif 
