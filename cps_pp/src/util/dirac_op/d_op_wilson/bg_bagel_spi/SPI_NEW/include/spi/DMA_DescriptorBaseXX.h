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


#ifndef _DMA_DESCRIPTORBASEXX_H_ /* Prevent multiple inclusion */
#define _DMA_DESCRIPTORBASEXX_H_

#include <spi/DMA_Packet.h>
#include <spi/DMA_Descriptors.h>

////////////////////////////////////////////////////////////////////////////////
///
/// \file spi/DMA_DescriptorBaseXX.h
///
/// \brief C++ DMA SPI Base Descriptor Classes
///
/// \see DMA_Descriptors.h for C interfaces
///
////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
///
/// \brief DMA_Descriptor Class
///
/// This class references a DMA_InjDescriptor_t which will ultimately be
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

class DMA_DescriptorBase : public DMA_InjDescriptor_t 
{
 public:
  
  //////////////////////////////////////////////////////////////////////////////
  ///
  /// \brief DMA_Descriptor constructor
  ///
  /// Zeros out the descriptor
  ///
  //////////////////////////////////////////////////////////////////////////////

  DMA_DescriptorBase()
  {
    DMA_ZeroOutDescriptor(this);              // Clear descriptor.
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
    this->prefetch_only = 1;
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
    this->local_memcopy = 1;
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
    this->base_offset = offset;
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
    int c;
    
    this->msg_length = len;
    
    c = DMA_PacketChunks(len); // Calculate number of 32 byte chunks in 
                               // the first packet.
    SPI_assert( c != 0 );
    this->hwHdr.Chunks = c - 1;// Packet header has 0 for 1 chunk, ... ,
                               // 7 for 8 chunks).
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
    return(this->msg_length);
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
    this->hwHdr.X = x;
    this->hwHdr.Y = y;
    this->hwHdr.Z = z;
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
    return(this->hwHdr.X);
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
    return(this->hwHdr.Y);
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
    return(this->hwHdr.Z);
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
    SPI_assert( (hints & 0x0000003F) == hints );

    this->hwHdr.Hint = hints;
  }
  
  //////////////////////////////////////////////////////////////////////////////
  ///
  /// \brief
  ///
  //////////////////////////////////////////////////////////////////////////////

  inline void setVirtualChannel(unsigned int vc)
  {
    SPI_assert( vc <= 3 );

    DMA_SetVc( this,
	       vc );
  }
  
  //////////////////////////////////////////////////////////////////////////////
  ///
  /// \brief
  ///
  //////////////////////////////////////////////////////////////////////////////

  inline void setInjCounterInfo(unsigned int groupID,
				unsigned int counterID)
  {
    SPI_assert( groupID   < DMA_NUM_COUNTER_GROUPS );
    SPI_assert( counterID < DMA_NUM_COUNTERS_PER_GROUP );

    this->idma_counterId =  
      counterID + groupID * DMA_NUM_COUNTERS_PER_GROUP;
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
    return (this->idma_counterId);
  }
  
  //////////////////////////////////////////////////////////////////////////////
  ///
  /// \brief
  ///
  //////////////////////////////////////////////////////////////////////////////

  inline void setRecFifoGroupID(unsigned int groupID)
  {
    SPI_assert( groupID < DMA_NUM_REC_FIFO_GROUPS );
    
    DMA_SetDescriptorPids(this,
			  groupID);
  }
  
  //////////////////////////////////////////////////////////////////////////////
  ///
  /// \brief
  ///
  //////////////////////////////////////////////////////////////////////////////

  inline void setRecCounterInfo(unsigned int groupID,
				unsigned int counterID)
  {
    DMA_SetDescriptorPids(this,
			  groupID);
    this->hwHdr.rDMA_Counter = 
    counterID + groupID * DMA_NUM_COUNTERS_PER_GROUP;
  }
  
  //////////////////////////////////////////////////////////////////////////////
  ///
  /// \brief
  ///
  //////////////////////////////////////////////////////////////////////////////

  inline void setSwArg(unsigned int swArg)
  {
    this->hwHdr.SW_Arg = swArg;
  }
  
  //////////////////////////////////////////////////////////////////////////////
  ///
  /// \brief
  ///
  //////////////////////////////////////////////////////////////////////////////

  inline unsigned int getSwArg() const
  {
    return(this->hwHdr.SW_Arg);
  }
  
  //////////////////////////////////////////////////////////////////////////////
  ///
  /// \brief
  ///
  //////////////////////////////////////////////////////////////////////////////

  inline void setFunctionID(unsigned int functionID)
  {
    SPI_assert(functionID < 256);

    this->hwHdr.Func_Id = functionID;
  }
    
  //////////////////////////////////////////////////////////////////////////////
  ///
  /// \brief
  ///
  //////////////////////////////////////////////////////////////////////////////

  inline void setSinglePacketParameter(unsigned int parameter)
  {
    this->hwHdr.Single_Packet_Parameter = parameter;
  }
  
  //////////////////////////////////////////////////////////////////////////////
  ///
  /// \brief
  ///
  //////////////////////////////////////////////////////////////////////////////

  inline unsigned int getSinglePacketParameter() const
  {
    return(this->hwHdr.Single_Packet_Parameter);
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
    SPI_assert( groupID < DMA_NUM_INJ_FIFO_GROUPS );
    SPI_assert( fifoNum < DMA_NUM_INJ_FIFOS_PER_GROUP );

    this->hwHdr.iDMA_Fifo_ID =  
      fifoNum + groupID * DMA_NUM_INJ_FIFOS_PER_GROUP;
  }
   
  //////////////////////////////////////////////////////////////////////////////
  ///
  /// \brief
  ///
  //////////////////////////////////////////////////////////////////////////////

  inline void setPutOffset(unsigned int putOffset)
  {
    this->hwHdr.Put_Offset = putOffset;
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
    this->hwHdr.Dp = 1;
  }

  //////////////////////////////////////////////////////////////////////////////
  ///
  /// \brief Clone the Descriptor
  ///
  /// This is a clone, copying "this" object to a specified DescriptorBase
  /// object.
  ///
  /// \param[in]  obj  Reference to the target DescriptorBase object to copy to.
  ///
  //////////////////////////////////////////////////////////////////////////////

  inline void clone (DMA_DescriptorBase & obj)
  {
    SPI_assert( sizeof(DMA_DescriptorBase) == 32 );

    _bgp_QuadMove( this, &obj, 0 );
    _bgp_QuadMove( ((char*)this)+16, ((char*)&obj)+16, 1);
  }
   
  //////////////////////////////////////////////////////////////////////////////
  ///
  /// \brief Dump this descriptor
  ///
  //////////////////////////////////////////////////////////////////////////////
  
  inline void dump()
  {
    printf("%08x: %08x %08x %08x %08x %08x %08x %08x %08x\n",
	   (unsigned)this,
	   (unsigned)(*(((unsigned*)this)+0)),
	   (unsigned)(*(((unsigned*)this)+1)),
	   (unsigned)(*(((unsigned*)this)+2)),
	   (unsigned)(*(((unsigned*)this)+3)),
	   (unsigned)(*(((unsigned*)this)+4)),
	   (unsigned)(*(((unsigned*)this)+5)),
	   (unsigned)(*(((unsigned*)this)+6)),
	   (unsigned)(*(((unsigned*)this)+7)));
  }

}; // End: class DMA_DescriptorBase




#endif
