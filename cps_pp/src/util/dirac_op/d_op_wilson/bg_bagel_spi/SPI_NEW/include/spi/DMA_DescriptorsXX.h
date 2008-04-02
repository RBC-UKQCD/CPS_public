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


#ifndef _DMA_DESCRIPTORSXX_H_ /* Prevent multiple inclusion */
#define _DMA_DESCRIPTORSXX_H_

#include <spi/DMA_Packet.h>
#include <spi/DMA_DescriptorBaseXX.h>
#include <bpcore/ppc450_inlines.h>

////////////////////////////////////////////////////////////////////////////////
///
/// \file spi/DMA_DescriptorsXX.h
///
/// \brief C++ DMA SPI Descriptor Classes
///
/// \see DMA_Descriptors.h for C interfaces
///
////////////////////////////////////////////////////////////////////////////////

/// \brief Indicator that setPreFetchOnly() has been called
#define DMA_DESC_SET_PREFETCHONLY         0x00000001

/// \brief Indicator that setLocalMemCopy() has been called
#define DMA_DESC_SET_LOCALMEMCOPY         0x00000002

/// \brief Indicator that setBaseOffset() has been called
#define DMA_DESC_SET_BASEOFFSET           0x00000004

/// \brief Indicator that setMsgLength() has been called
#define DMA_DESC_SET_MSGLENGTH            0x00000008

/// \brief Indicator that setDest() has been called
#define DMA_DESC_SET_DEST                 0x00000010

/// \brief Indicator that setHints() has been called
#define DMA_DESC_SET_HINTS                0x00000020

/// \brief Indicator that setVirtualChannel() has been called
#define DMA_DESC_SET_VIRTUALCHANNEL       0x00000040

/// \brief Indicator that setInjCounterInfo() has been called
#define DMA_DESC_SET_INJCOUNTERINFO       0x00000080

/// \brief Indicator that setRecFifoGroupID() has been called
#define DMA_DESC_SET_RECFIFOGROUPID       0x00000100

/// \brief Indicator that setRecCounterInfo() has been called
#define DMA_DESC_SET_RECCOUNTERINFO       0x00000200

/// \brief Indicator that setSwArg() has been called
#define DMA_DESC_SET_SWARG                0x00000400

/// \brief Indicator that setFunctionID() has been called
#define DMA_DESC_SET_FUNCTIONID           0x00000800

/// \brief Indicator that setFifoNum() has been called
#define DMA_DESC_SET_FIFONUM              0x00001000

/// \brief Indicator that setRemoteGetInjFifoInfo() has been called
#define DMA_DESC_SET_REMOTEGETINJFIFOINFO 0x00002000

/// \brief Indicator that setPutOffset() has been called
#define DMA_DESC_SET_PUTOFFSET            0x00004000

////////////////////////////////////////////////////////////////////////////////
///
/// \brief DMA_Descriptor Class
///
/// This class contains a DMA_InjDescriptor_t which will ultimately be
/// injected into a DMA injection fifo.  There are setter and getter functions
/// to operate on the descriptor's data members.
///
/// This class also contains some additional data members that assist with
/// managing the descriptor.
///
/// This class is intended to be used as a base class for the various types
/// of descriptors.  Thus, users of descriptors should instantiate those
/// derived classes, not this base class.  Derived classes are also defined in
/// this include file.
///
/// The derived classes document which setter functions must be called to set
/// all of the necessary fields into the descriptor.  There are also optional
/// setters that may be called.  When the descriptor is ready to be injected,
/// the user should invoke checkIfReadyToInject() to verify that all required
/// fields have been set.
///
////////////////////////////////////////////////////////////////////////////////

class DMA_Descriptor
{
 public:

  //////////////////////////////////////////////////////////////////////////////
  ///
  /// \brief DMA_Descriptor constructor
  ///
  /// The caller must specify indicators (defined above) of the functions that
  /// must be called to set all of the necessary fields into the descriptor.
  ///
  /// \param[in]  fieldsThatMustBeSet  Indicators of the fields that must be
  ///                                  set before the descriptor can be
  ///                                  injected.
  /// \param[in]  desc                 Pointer to the descriptor that is
  ///                                  managed by this class.
  ///
  //////////////////////////////////////////////////////////////////////////////

  DMA_Descriptor(unsigned int fieldsThatMustBeSet,
		 DMA_DescriptorBase *desc=NULL) :
  _desc(desc),
  _fieldsSet(fieldsThatMustBeSet),
  _fifo(-1),
  _numInjected(0),
  _cb_done(NULL),
  _clientdata(NULL),
  _nextPtr(NULL),
  _prevPtr(NULL)
  { }

  //////////////////////////////////////////////////////////////////////////////
  ///
  /// \brief Default DMA_Descriptor Constructor
  ///
  /// Not to be directly used
  ///
  //////////////////////////////////////////////////////////////////////////////

  DMA_Descriptor() {}

  //////////////////////////////////////////////////////////////////////////////
  ///
  /// \brief Set Descriptor Pointer
  ///
  /// Set the pointer to the descriptor managed by this class
  ///
  //////////////////////////////////////////////////////////////////////////////

  inline void setDescriptorPtr(DMA_DescriptorBase *desc)
  {
    _desc = desc;
  }

  //////////////////////////////////////////////////////////////////////////////
  ///
  /// \brief Get Descriptor Pointer
  ///
  /// Get the pointer to the descriptor managed by this class
  ///
  //////////////////////////////////////////////////////////////////////////////

  inline DMA_DescriptorBase * getDescriptorPtr()
  {
    return (_desc);
  }

  //////////////////////////////////////////////////////////////////////////////
  ///
  /// \brief Set the Checksum limit
  ///
  /// \param[in]      csum_skip  The number of 2-bytes to skip in the checksum
  ///                            (7 bits).
  /// \param[in]      skip       The checksum skip attribute:
  ///                            0 = The packets participates in the injection
  ///                                checksum.
  ///                            1 = The packet does not participate in the
  ///                                injection checksum.
  ///
  //////////////////////////////////////////////////////////////////////////////

  inline void setInjCsum(unsigned csum_skip,
                         unsigned skip)
  {
    SPI_assert( skip <=1 );

    _desc->hwHdr.CSum_Skip = csum_skip;
    _desc->hwHdr.Sk        = skip;
  }

  //////////////////////////////////////////////////////////////////////////////
  ///
  /// \brief Set the "Prefetch Only" indicator.
  ///
  /// The "prefetch only" indicator tells the DMA that the data should not be
  /// transferred.  It should only be fetched into L3 cache.
  ///
  //////////////////////////////////////////////////////////////////////////////

  inline void setPrefetchOnly()
  {
    indicateFieldSet(DMA_DESC_SET_PREFETCHONLY);

    _desc->setPrefetchOnly();
  }

  //////////////////////////////////////////////////////////////////////////////
  ///
  /// \brief Set the "Local Memory Copy" indicator.
  ///
  /// The "local memory copy" indicator tells the DMA that the data will be
  /// transferred to the same node.  Thus, it is simplified to a memory
  /// transfer.
  ///
  //////////////////////////////////////////////////////////////////////////////

  inline void setLocalMemcopy()
  {
    indicateFieldSet(DMA_DESC_SET_LOCALMEMCOPY);

    _desc->setLocalMemcopy();
  }

  //////////////////////////////////////////////////////////////////////////////
  ///
  /// \brief Set the Base Offset Field of the Descriptor
  ///
  /// The base offset field is the offset of the data buffer from the base
  /// offset in the injection counter associated with this descriptor.
  /// It tells the DMA where the data is located.
  ///
  /// \param[in]  offset  The offset in bytes
  ///
  //////////////////////////////////////////////////////////////////////////////

  inline void setBaseOffset(unsigned int offset)
  {
    indicateFieldSet(DMA_DESC_SET_BASEOFFSET);

    _desc->setBaseOffset(offset);
  }

  //////////////////////////////////////////////////////////////////////////////
  ///
  /// \brief Set the Message Length Field of the Descriptor
  ///
  /// The message length field is the length in bytes of the data to be
  /// transferred.  This field is set into the descriptor, and the number of
  /// 32-byte chunks in the first packet is calculated and set into the packet
  /// header in the descriptor.
  ///
  /// \param[in]  len  The length in bytes
  ///
  //////////////////////////////////////////////////////////////////////////////

  inline void setMsgLength(unsigned int len)
  {
    indicateFieldSet(DMA_DESC_SET_MSGLENGTH);

    _desc->setMsgLength(len);

  }

  //////////////////////////////////////////////////////////////////////////////
  ///
  /// \brief Get the Message Length
  ///
  /// The message length field in the descriptor is returned.
  ///
  /// \retval  msgLength  The length of the message in bytes.
  ///
  //////////////////////////////////////////////////////////////////////////////

  inline unsigned int getMsgLength() const
  {
    return( _desc->getMsgLength() );
  }

  //////////////////////////////////////////////////////////////////////////////
  ///
  /// \brief  Set the Destination Fields in the Descriptor
  ///
  /// The x, y, and z coordinates of the destination are set into the
  /// descriptor.
  ///
  /// \param[in]  x  The x coordinate of the destination
  /// \param[in]  y  The y coordinate of the destination
  /// \param[in]  z  The z coordinate of the destination
  ///
  //////////////////////////////////////////////////////////////////////////////

  inline void setDest(unsigned int x,
		      unsigned int y,
		      unsigned int z)
  {
    indicateFieldSet(DMA_DESC_SET_DEST);

    _desc->setDest(x,y,z);
  }

  //////////////////////////////////////////////////////////////////////////////
  ///
  /// \brief Get the X Coordinate of the Destination
  ///
  /// \retval  X  The x coordinate of the destination
  ///
  //////////////////////////////////////////////////////////////////////////////

  inline unsigned int x() const
  {
    return( _desc->x() );
  }

  //////////////////////////////////////////////////////////////////////////////
  ///
  /// \brief Get the Y Coordinate of the Destination
  ///
  /// \retval  Y  The y coordinate of the destination
  ///
  //////////////////////////////////////////////////////////////////////////////

  inline unsigned int y() const
  {
    return( _desc->y() );
  }

  //////////////////////////////////////////////////////////////////////////////
  ///
  /// \brief Get the Z Coordinate of the Destination
  ///
  /// \retval  Z  The z coordinate of the destination
  ///
  //////////////////////////////////////////////////////////////////////////////

  inline unsigned int z() const
  {
    return( _desc->z() );
  }

  //////////////////////////////////////////////////////////////////////////////
  ///
  /// \brief Set the Hint Bits Into the Descriptor's Hardware Packet Header
  ///
  /// \param[in]  hints  The 6 hint bits, located in the low-order bits of this
  ///                    field.
  ///
  //////////////////////////////////////////////////////////////////////////////

  inline void setHints(unsigned int hints)
  {
    indicateFieldSet(DMA_DESC_SET_HINTS);

    _desc->setHints(hints);
  }

  //////////////////////////////////////////////////////////////////////////////
  ///
  /// \brief
  ///
  //////////////////////////////////////////////////////////////////////////////

  inline void setVirtualChannel(unsigned int vc)
  {
    indicateFieldSet(DMA_DESC_SET_VIRTUALCHANNEL);

    _desc->setVirtualChannel(vc);
  }

  //////////////////////////////////////////////////////////////////////////////
  ///
  /// \brief
  ///
  //////////////////////////////////////////////////////////////////////////////

  inline void setInjCounterInfo(unsigned int groupID,
				unsigned int counterID)
  {
    indicateFieldSet(DMA_DESC_SET_INJCOUNTERINFO);

    _desc->setInjCounterInfo(groupID,
			     counterID);
  }

  //////////////////////////////////////////////////////////////////////////////
  ///
  /// \brief Get Injection Counter Number
  ///
  /// \retval  counterNumber  The global counter number
  ///
  //////////////////////////////////////////////////////////////////////////////

  inline int getInjCounterNum()
  {
    return ( _desc->getInjCounterNum() );
  }
  
  //////////////////////////////////////////////////////////////////////////////
  ///
  /// \brief
  ///
  //////////////////////////////////////////////////////////////////////////////

  inline void setRecFifoGroupID(unsigned int groupID)
  {
    indicateFieldSet(DMA_DESC_SET_RECFIFOGROUPID);

    _desc->setRecFifoGroupID(groupID);
  }

  //////////////////////////////////////////////////////////////////////////////
  ///
  /// \brief
  ///
  //////////////////////////////////////////////////////////////////////////////

  inline void setRecCounterInfo(unsigned int groupID,
				unsigned int counterID)
  {
    indicateFieldSet(DMA_DESC_SET_RECCOUNTERINFO);

    _desc->setRecCounterInfo(groupID,
			     counterID);
  }

  //////////////////////////////////////////////////////////////////////////////
  ///
  /// \brief
  ///
  //////////////////////////////////////////////////////////////////////////////

  inline void setSwArg(unsigned int swArg)
  {
    indicateFieldSet(DMA_DESC_SET_SWARG);

    _desc->setSwArg(swArg);
  }

  //////////////////////////////////////////////////////////////////////////////
  ///
  /// \brief
  ///
  //////////////////////////////////////////////////////////////////////////////

  inline unsigned int getSwArg() const
  {
    return( _desc->getSwArg() );
  }

  //////////////////////////////////////////////////////////////////////////////
  ///
  /// \brief
  ///
  //////////////////////////////////////////////////////////////////////////////

  inline void setFunctionID(unsigned int functionID)
  {
    indicateFieldSet(DMA_DESC_SET_FUNCTIONID);

    _desc->setFunctionID(functionID);
  }

  //////////////////////////////////////////////////////////////////////////////
  ///
  /// \brief
  ///
  //////////////////////////////////////////////////////////////////////////////

  inline void setFifoNum(int fnum)
  {
    indicateFieldSet(DMA_DESC_SET_FIFONUM);

    this->_fifo = fnum;
  }

  //////////////////////////////////////////////////////////////////////////////
  ///
  /// \brief
  ///
  //////////////////////////////////////////////////////////////////////////////

  inline short getFifoNum() const
  {
    return( this->_fifo );
  }

  //////////////////////////////////////////////////////////////////////////////
  ///
  /// \brief
  ///
  //////////////////////////////////////////////////////////////////////////////

  inline void setSinglePacketParameter(unsigned int parameter)
  {
    _desc->setSinglePacketParameter(parameter);
  }

  //////////////////////////////////////////////////////////////////////////////
  ///
  /// \brief
  ///
  //////////////////////////////////////////////////////////////////////////////

  inline unsigned int getSinglePacketParameter() const
  {
    return( _desc->getSinglePacketParameter() );
  }

  //////////////////////////////////////////////////////////////////////////////
  ///
  /// \brief
  ///
  //////////////////////////////////////////////////////////////////////////////

  inline void setDMAFifoStatus(unsigned long long status)
  {
    this->_numInjected = status;
  }

  //////////////////////////////////////////////////////////////////////////////
  ///
  /// \brief
  ///
  //////////////////////////////////////////////////////////////////////////////

  inline unsigned long long getNumInjected() const
  {
    return(this->_numInjected);
  }

  //////////////////////////////////////////////////////////////////////////////
  ///
  /// \brief Set Remote Get Injection Fifo Information
  ///
  /// This sets the injection fifo group number and fifo number where the
  /// DMA that receives the remote get packet will inject the payload.
  ///
  /// \param[in] groupID  The injection fifo group number
  /// \param[in] fifoNum  The injection fifo number
  ///
  //////////////////////////////////////////////////////////////////////////////

  inline void setRemoteGetInjFifoInfo(unsigned int groupID,
				      int          fifoNum)
  {
    indicateFieldSet(DMA_DESC_SET_REMOTEGETINJFIFOINFO);

    _desc->setRemoteGetInjFifoInfo(groupID,
				   fifoNum);
  }

  //////////////////////////////////////////////////////////////////////////////
  ///
  /// \brief Set the putOffset field of the descriptor
  ///
  /// This is the offset to the buffer being sent.
  ///
  /// \param[in]  putOffset  The put offset
  ///
  //////////////////////////////////////////////////////////////////////////////

  inline void setPutOffset(unsigned int putOffset)
  {
    indicateFieldSet(DMA_DESC_SET_PUTOFFSET);

    _desc->setPutOffset(putOffset);
  }

  //////////////////////////////////////////////////////////////////////////////
  ///
  /// \brief Set the descriptor callback information
  ///
  /// \param[in]  cb_done     Callback function
  /// \param[in]  clientdata  Pointer to opaque data for the callback function
  ///
  //////////////////////////////////////////////////////////////////////////////

  inline void setDoneCallback ( void (*cb_done)(void *),
			        void * clientdata )
  {
    _cb_done    = cb_done;
    _clientdata = clientdata;
  }

  //////////////////////////////////////////////////////////////////////////////
  ///
  /// \brief Return Whether a Callback is Desired When this Descriptor has
  ///        been Processed by the DMA
  ///
  /// \retval  0  No callback is desired
  /// \retval  1  A callback is desired
  ///
  //////////////////////////////////////////////////////////////////////////////

  inline unsigned int isCallbackDesired ( )
  {
    return ( _cb_done != NULL );
  }

  //////////////////////////////////////////////////////////////////////////////
  ///
  /// \brief Perform the callback
  ///
  //////////////////////////////////////////////////////////////////////////////

  inline void performCallback ( )
  {
    this->_cb_done ( this->_clientdata );
  }

  //////////////////////////////////////////////////////////////////////////////
  ///
  /// \brief Set Next Pointer
  ///
  /// Set the "next" pointer in this descriptor to point to the next descriptor
  /// in the list.
  ///
  /// \param[in]  nextPtr  Next pointer to be set
  ///
  //////////////////////////////////////////////////////////////////////////////

  inline void setNextPtr ( DMA_Descriptor *nextPtr )
  {
    this->_nextPtr = nextPtr;
  }

  //////////////////////////////////////////////////////////////////////////////
  ///
  /// \brief Get Next Pointer
  ///
  /// Get the "next" pointer in this descriptor
  ///
  /// \retval  nextPtr  Pointer to next descriptor in the list
  ///
  //////////////////////////////////////////////////////////////////////////////

  inline DMA_Descriptor * getNextPtr ( )
  {
    return ( this->_nextPtr );
  }

  //////////////////////////////////////////////////////////////////////////////
  ///
  /// \brief Set the Deposit Bit (for a Broadcast)
  ///
  /// This function is to be called on a DMA_MemoryFifoDescriptor or a
  /// DMA_DirectPutDescriptor to designate that a broadcast is to be done.
  /// The packets will be sent out across the torus along one line.
  /// This requires that the X, Y, Z parameters be specified such that two
  /// of the coordinates match the origin node, and the third is different,
  /// designating the direction of the travel.  The hints parameter must also
  /// designate the direction of the travel...only the one hint bit should be
  /// specified.
  ///
  //////////////////////////////////////////////////////////////////////////////

  inline void setDepositBit()
  {
    _desc->setDepositBit();
  }

  //////////////////////////////////////////////////////////////////////////////
  ///
  /// \brief
  ///
  //////////////////////////////////////////////////////////////////////////////

  inline void checkIfReadyToInject() const
  {
#ifndef ASSERT_ABORT
#ifndef ASSERT_PROD
    // If all of the required fields have been set, _fieldsSet will
    // be zero.  If not, assert.
    if (_fieldsSet != 0)
    {
      printf("DMA_DescriptorsXX.h: Fields Not Set = 0x%08x\n",_fieldsSet);
      SPI_assert(0);
    }
#endif
#endif
  }

  //////////////////////////////////////////////////////////////////////////////
  ///
  /// \brief Dump Descriptor
  ///
  //////////////////////////////////////////////////////////////////////////////

  inline void dump()
  {
    printf("Descriptor fieldsSet = 0x%08x\n",_fieldsSet);
    _desc->dump();
  }

 protected:

  ///
  /// \brief _desc       Pointer to descriptor managed by this class
  ///
  DMA_DescriptorBase * _desc;

  ///
  /// \brief _fieldsSet  Bits indicating which fields in the descriptor must be
  ///                    set before the descriptor can be injected into a fifo.
  ///                    The constructor of the derived descriptor class sets
  ///                    the bits, and the corresponding functions turn their
  ///                    bit off.  They should all be off when the descriptor
  ///                    is ready to inject.  Call checkIfReadyToInject()
  ///                    prior to injecting to verify.
  ///
  /// \see checkIfReadyToInject()
  unsigned int       _fieldsSet;


  ///
  /// \brief _fifo  The fifo number (relative to the fifo group) where the
  ///               descriptor will be injected.
  /// \see setFifoNum()
  /// \see getFifoNum()
  ///
  short              _fifo;

  ///
  /// \brief _numInjected  The sequence number assigned to this descriptor when
  ///                      it is injected.  This sequence number grows
  ///                      sequentially (won't wrap in the lifetime of a job).
  ///                      It is recalled and used later to determine if the
  ///                      descriptor has been processed by the DMA (it has left
  ///                      the fifo).
  ///
  /// \see setDMAFifoStatus()
  /// \see getNumInjected()
  ///
  unsigned long long _numInjected;

  ///
  /// \brief _cb_done  Callback when this descriptor has completed.
  ///                  This may be NULL when counters are used for completion.
  ///                  Normally, a descriptor that uses a shared counter will
  ///                  specify this so "descriptor" notification will take place.
  ///                  A descriptor that uses an exclusive counter will leave
  ///                  this NULL because "counter" notification will take place.
  ///
  void (*_cb_done)(void *);

  ///
  /// \brief  _clientdata  A pointer to opaque data to be passed to _cb_done.
  ///
  void *_clientdata;

  ///
  /// \brief  nextPtr  Pointer to next descriptor on the queue
  ///
  DMA_Descriptor *_nextPtr;

  ///
  /// \brief  prevPtr  Pointer to previous descriptor on the queue
  ///
  DMA_Descriptor *_prevPtr;


 private:

  inline void indicateFieldSet(unsigned int field)
  {
#ifndef ASSERT_ABORT
#ifndef ASSERT_PROD
    _fieldsSet &= ~field;
#endif
#endif
  }

}; // End: class DMA_Descriptor


////////////////////////////////////////////////////////////////////////////////
///
/// \brief DMA Memory Fifo Descriptor Class
///
/// This class is derived from the DMA_Descriptor class.  It establishes the
/// list of fields that must be set before the descriptor can be injected.
///
/// It also defaults the descriptor's hardware packet header such that the
/// packets associated with this descriptor will be routed deterministically,
/// so that all of the packets arrive at the destination in order.  This can be
/// overridden using setVirtualChannel().
///
/// It also hard-codes the descriptor's hardware packet header such that the
/// entire packet except the header is checksum'd.
///
/// The following functions must be called before injecting the descriptor:
///
/// DMA_MemoryFifoDescriptor() (constructor)
/// setBaseOffset()
/// setMsgLength()
/// setDest()
/// setInjCounterInfo()
/// setRecFifoGroupID()
/// setSwArg()
/// setFunctionID()
/// setFifoNum()
/// setVirtualChannel(DMA_PACKET_VC_D0) (optional, to set dynamic routing)
/// setLocalMemcopy() (optional, to set local copy (non-torus) mode)
/// setPreFetchOnly() (optional, when in local copy mode, prefetch only, but
///                    do not copy)
/// setDepositBit()   (optional, for a broadcast)
/// checkIfReadyToInject() (last)
///
////////////////////////////////////////////////////////////////////////////////

class DMA_MemoryFifoDescriptor : public DMA_Descriptor
{
 public:

  DMA_MemoryFifoDescriptor( DMA_DescriptorBase * desc ) :
    DMA_Descriptor(DMA_DESC_SET_BASEOFFSET     |
		   DMA_DESC_SET_MSGLENGTH      |
		   DMA_DESC_SET_DEST           |
		   DMA_DESC_SET_INJCOUNTERINFO |
		   DMA_DESC_SET_RECFIFOGROUPID |
		   DMA_DESC_SET_SWARG          |
		   DMA_DESC_SET_FUNCTIONID     |
		   DMA_DESC_SET_FIFONUM,
		   desc)
  {
    setVirtualChannel(DMA_PACKET_VC_BN); // Default to deterministic routing.

    _desc->hwHdr.CSum_Skip = DMA_CSUM_SKIP;    // Checksum all but header
    _desc->hwHdr.Sk        = DMA_CSUM_BIT;     // Checksum entire packet
  }

  //////////////////////////////////////////////////////////////////////////////
  ///
  /// \brief Default DMA_MemoryFifoDescriptor Constructor
  ///
  /// Not to be directly used
  ///
  //////////////////////////////////////////////////////////////////////////////

  DMA_MemoryFifoDescriptor() {}

  //////////////////////////////////////////////////////////////////////////////
  ///
  /// \brief Clone the Descriptor
  ///
  /// This is a clone, copying "this" object to a specified Descriptor object.
  ///
  /// \param[in]  obj  Reference to the target Descriptor object to copy to.
  ///
  //////////////////////////////////////////////////////////////////////////////

  inline void clone (DMA_MemoryFifoDescriptor & obj)
  {
    SPI_assert( sizeof(DMA_MemoryFifoDescriptor) == 40 );

    // Move first 32 bytes
    _bgp_QuadMove( this, &obj, 0 );
    _bgp_QuadMove( ((char*)this)+16, ((char*)&obj)+16, 1);

    // Move last 8 bytes
    *(((unsigned*)&obj)+8) = *(((unsigned*)this)+8);
    *(((unsigned*)&obj)+9) = *(((unsigned*)this)+9);
  }

 protected:


 private:


}; // End: class DMA_MemoryFifoDescriptor


////////////////////////////////////////////////////////////////////////////////
///
/// \brief DMA Direct Put Descriptor Class
///
/// A Direct Put Descriptor is injected by the core into an injection fifo,
/// or it may be the payload of a remote get descriptor such that it is
/// injected into a remote get injection fifo by the DMA.  This descriptor
/// contains information about the origin side and the target side.  This
/// allows the data to be put directly into the target buffer by the DMA.
///
/// This class is derived from the DMA_Descriptor class.  It establishes the
/// list of fields that must be set before the descriptor can be injected.
///
/// It defaults the descriptor's hardware packet header such that the
/// packets associated with this descriptor will be routed dynamically
/// for performance.  Packets may arrive at the destination out of order.
/// This can be overridden using setVirtualChannel().
///
/// It also hard-codes the descriptor's hardware packet header such that the
/// entire packet except the header is checksum'd.
///
/// It hard-codes the packet to be a DMA packet (not a memory fifo packet).
///
/// The following functions must be called before injecting the descriptor:
///
/// DMA_DirectPutDescriptor() (constructor)
/// setBaseOffset()
/// setMsgLength()
/// setDest()
/// setInjCounterInfo()
/// setRecFifoGroupID()
/// setPutOffset()
/// setFifoNum()
/// setVirtualChannel(DMA_PACKET_VC_BN) (optional, to set deterministic routing)
/// setLocalMemcopy() (optional, to set local copy (non-torus) mode)
/// setPreFetchOnly() (optional, when in local copy mode, prefetch only, but
///                    do not copy)
/// setDepositBit()   (optional, for a broadcast)
/// checkIfReadyToInject() (last)
///
////////////////////////////////////////////////////////////////////////////////

class DMA_DirectPutDescriptor : public DMA_Descriptor
{
 public:
  DMA_DirectPutDescriptor( DMA_DescriptorBase * desc ) :
    DMA_Descriptor(DMA_DESC_SET_BASEOFFSET     |
		   DMA_DESC_SET_MSGLENGTH      |
		   DMA_DESC_SET_DEST           |
		   DMA_DESC_SET_INJCOUNTERINFO |
		   DMA_DESC_SET_RECFIFOGROUPID |
		   DMA_DESC_SET_PUTOFFSET      |
		   DMA_DESC_SET_FIFONUM,
		   desc)
  {
    setVirtualChannel(DMA_PACKET_VC_D0); // Default to dynamic routing.

    _desc->hwHdr.CSum_Skip = DMA_CSUM_SKIP;    // Checksum all but header
    _desc->hwHdr.Sk        = DMA_CSUM_BIT;     // Checksum entire packet
    _desc->hwHdr.Dm        = 1;                // 1=DMA Mode, 0=Fifo Mode
  }

  //////////////////////////////////////////////////////////////////////////////
  ///
  /// \brief Default DMA_DirectPutDescriptor Constructor
  ///
  /// Not to be directly used
  ///
  //////////////////////////////////////////////////////////////////////////////

  DMA_DirectPutDescriptor() {}

  //////////////////////////////////////////////////////////////////////////////
  ///
  /// \brief Clone the Descriptor
  ///
  /// This is a clone, copying "this" object to a specified Descriptor object.
  ///
  /// \param[in]  obj  Reference to the target Descriptor object to copy to.
  ///
  //////////////////////////////////////////////////////////////////////////////

  inline void clone (DMA_DirectPutDescriptor & obj)
  {
    SPI_assert( sizeof(DMA_DirectPutDescriptor) == 40 );

    // Move first 32 bytes
    _bgp_QuadMove( this, &obj, 0 );
    _bgp_QuadMove( ((char*)this)+16, ((char*)&obj)+16, 1);

    // Move last 8 bytes
    *(((unsigned*)&obj)+8) = *(((unsigned*)this)+8);
    *(((unsigned*)&obj)+9) = *(((unsigned*)this)+9);
  }

 protected:


 private:

}; // End: class DMA_DirectPutDescriptor


////////////////////////////////////////////////////////////////////////////////
///
/// \brief DMA Remote Get Descriptor Class
///
/// A Remote Get Descriptor is injected by the core into an injection fifo.
/// This descriptor refers to a payload that is itself a descriptor that
/// will be injected into a remote get injection fifo on the node that
/// receives the remote get packet.
///
/// This class is derived from the DMA_Descriptor class.  It establishes the
/// list of fields that must be set before the descriptor can be injected.
///
/// It defaults the descriptor's hardware packet header such that the
/// packets associated with this descriptor will be routed deterministically,
/// so that all of the packets arrive at the destination in order.  This can be
/// overridden using setVirtualChannel().
///
/// It hard-codes the descriptor's hardware packet header such that the
/// entire packet is not checksum'd.
///
/// It hard-codes the size of the packet to be 64 bytes.
///
/// It hard-codes the packet to be a DMA packet (not a memory fifo packet).
///
/// It hard-codes the packet to be a remote get packet.
///
/// The following functions must be called before injecting the descriptor:
///
/// DMA_RemoteGetDescriptor() (constructor)
/// setBaseOffset()
/// setMsgLength()
/// setDest()
/// setInjCounterInfo()
/// setRecFifoGroupID()
/// setRemoteGetInjFifoInfo()
/// setLocalMemcopy() (optional, to set local copy (non-torus) mode)
/// setFifoNum()
/// setVirtualChannel(DMA_PACKET_VC_D0) (optional, to set dynamic routing)
/// checkIfReadyToInject() (last)
////////////////////////////////////////////////////////////////////////////////

class DMA_RemoteGetDescriptor : public DMA_Descriptor
{
 public:

  DMA_RemoteGetDescriptor( DMA_DescriptorBase * desc ) :
    DMA_Descriptor(DMA_DESC_SET_BASEOFFSET           |
		   DMA_DESC_SET_MSGLENGTH            |
		   DMA_DESC_SET_DEST                 |
		   DMA_DESC_SET_INJCOUNTERINFO       |
		   DMA_DESC_SET_RECFIFOGROUPID       |
		   DMA_DESC_SET_REMOTEGETINJFIFOINFO |
		   DMA_DESC_SET_FIFONUM,
		   desc)
  {
    setVirtualChannel(DMA_PACKET_VC_BN); // Default to deterministic routing
    _desc->hwHdr.Sk        = 1;           // Don't checksum this packet
    _desc->hwHdr.Chunks    = 1;           // Size in Chunks of 32B 1 => 64 bytes
    _desc->hwHdr.Dm        = 1;           // 1=DMA Mode, 0=Fifo Mode
    _desc->hwHdr.Flags     = 0x01;        // Flags[7]=Remote-Get
  }

  //////////////////////////////////////////////////////////////////////////////
  ///
  /// \brief Default DMA_RemoteGetDescriptor Constructor
  ///
  /// Not to be directly used
  ///
  //////////////////////////////////////////////////////////////////////////////

  DMA_RemoteGetDescriptor() {}

  //////////////////////////////////////////////////////////////////////////////
  ///
  /// \brief Clone the Descriptor
  ///
  /// This is a clone, copying "this" object to a specified Descriptor object.
  ///
  /// \param[in]  obj  Reference to the target Descriptor object to copy to.
  ///
  //////////////////////////////////////////////////////////////////////////////

  inline void clone (DMA_RemoteGetDescriptor & obj)
  {
    SPI_assert( sizeof(DMA_RemoteGetDescriptor) == 40 );

    // Move first 32 bytes
    _bgp_QuadMove( this, &obj, 0 );
    _bgp_QuadMove( ((char*)this)+16, ((char*)&obj)+16, 1);

    // Move last 8 bytes
    *(((unsigned*)&obj)+8) = *(((unsigned*)this)+8);
    *(((unsigned*)&obj)+9) = *(((unsigned*)this)+9);
  }

 protected:


 private:


}; // End: class DMA_RemoteGetDescriptor



#endif
