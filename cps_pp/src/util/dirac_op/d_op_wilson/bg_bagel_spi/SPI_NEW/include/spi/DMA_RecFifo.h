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


#ifndef _DMA_RECFIFO_H_ /* Prevent multiple inclusion */
#define _DMA_RECFIFO_H_


/*! 
 * \file spi/DMA_RecFifo.h
 *
 * \brief DMA SPI Reception Fifo Definitions and Inline Functions
 * 
 * There are 
 * - 6 normal-priority torus hardware fifos
 * - 1 high-priority torus hardware fifo
 * - 1 local memcpy reception fifo
 * If we assigned a reception fifo to each, there would be 8.  These are called
 * "normal reception fifos".
 *
 * There is one reception fifo used to store packet headers only (for debug).
 * This is called the "header reception fifo".
 *
 * The hardware packet header's (Pid0, Pid1) bits specify up to four processors.
 * There is one set of "normal" and one "header" reception fifos per processor,
 * called a "reception fifo group".  Thus, there are 4 groups.
 *
 */


#include <common/namespace.h>
#include <common/bgp_ras.h>


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




#include <spi/DMA_Fifo.h>
#include <spi/DMA_Assert.h>
#include <spi/DMA_Packet.h>




/*!
 * \brief Number of Normal (non-header) Reception Fifos Per Group
 * 
 * These will have fifo IDs 0 through 7 in a group.
 * 
 */
#define  DMA_NUM_NORMAL_REC_FIFOS_PER_GROUP 8


/*!
 * \brief Number of Header Reception Fifos Per Group
 *
 */
#define  DMA_NUM_HEADER_REC_FIFOS_PER_GROUP 1


/*!
 * \brief Fifo ID of the Header Reception Fifo in a group.
 *
 * This will be fifo ID 8 in a group.
 */
#define  DMA_HEADER_REC_FIFO_ID  (DMA_NUM_NORMAL_REC_FIFOS_PER_GROUP)


/*!
 * \brief Number of Reception Fifos Per Group
 *
 * 8 Normal + 1 Header
 */
#define  DMA_NUM_REC_FIFOS_PER_GROUP (DMA_NUM_NORMAL_REC_FIFOS_PER_GROUP + \
                                      DMA_NUM_HEADER_REC_FIFOS_PER_GROUP)


/*!
 * \brief Number of Reception Fifo Groups
 *
 * One group for each processor, identified by (Pid0,Pid1) in the packet header.
 */
#define  DMA_NUM_REC_FIFO_GROUPS 4


/*!
 * \brief Total Number of Normal Reception Fifos
 */
#define  DMA_NUM_NORMAL_REC_FIFOS (DMA_NUM_REC_FIFO_GROUPS * \
                                   DMA_NUM_NORMAL_REC_FIFOS_PER_GROUP)

/*!
 * \brief Total Number of Header Reception Fifos
 */
#define  DMA_NUM_HEADER_REC_FIFOS (DMA_NUM_REC_FIFO_GROUPS * \
                                   DMA_NUM_HEADER_REC_FIFOS_PER_GROUP)


/*!
 * \brief The maximum number of packets that can be processed by a polling
 *        function before it must update the fifo's hardware head.  The
 *        polling function can keep track of the head in the va_head shadow
 *        and only update the hardware head when this limit is reached to
 *        reduce overhead.
 */
#define  DMA_MAX_NUM_PACKETS_BEFORE_MOVING_HEAD 32


/*!
 * \brief Minimum Reception Fifo Size in bytes
 *
 * It is important that this size be enough to hold more packets than 
 * DMA_MAX_NUM_PACKETS_BEFORE_MOVING_HEAD.  Otherwise, the polling
 * function may deadlock with the DMA (the DMA considers the fifo full,
 * but we have actually processed all of the packets).
 * Additionally, we add 512 bytes to this, since the DMA will only fill
 * the fifo to fifoSize - 512 bytes.
 * 
 *
 */
#define DMA_MIN_REC_FIFO_SIZE_IN_BYTES (512 + (256 * DMA_MAX_NUM_PACKETS_BEFORE_MOVING_HEAD))

/*!
 * \brief DMA Reception Fifo Map Structure
 *
 * This structure defines the basic reception fifo configuration.
 * It is common across all reception fifo groups.
 *
 * Example 1:
 * In SMP mode you might have only two reception fifos:
 * - One for normal-priority torus and local transfers, and 
 * - one for high-priority torus transfers.
 * In this case, only one fifo group would be needed.
 *
 * Example 2:
 * In virtual node mode, you might have two reception fifos per group (as
 * described in Example 1), and 4 groups, one for each processor.
 * All packets with the same (pid0,pid1) bits use the same group.
*/
typedef struct DMA_RecFifoMap_t
{
  unsigned short int save_headers; /*!< Flag that specifies whether header      
                                        will be used to store packet headers.   
                                        - 0 => Do not store headers               
                                        - 1 => Store headers (debug mode)     */

  unsigned int fifo_types[ DMA_NUM_NORMAL_REC_FIFOS ]; /*!< The type of each    
                                        normal rec fifo.  If entry i is         
                                        - 0, rec fifo i is type 0                 
                                        - 1, rec fifo i is type 1
                 
                                        For type i fifo, threshold interrupt    
                                        fires if fifo free space <=             
                                        threshold[i], in units of 16B quads.  */

  unsigned int hdr_fifo_types[ DMA_NUM_HEADER_REC_FIFOS ]; /*!< The type of     
                                        each header reception fifo.  If entry   
                                        i is                                    
                                        - 0, header rec fifo i is type 0          
                                        - 1, header rec fifo i is type 1          

                                        For type i fifo, threshold interrupt    
                                        fires if fifo free space <=             
                                        threshold[i], in units of 16B quads.  */

  unsigned int threshold[2];       /*!< For type i fifo, threshold interrupt    
                                        fires if fifo free space <=             
                                        threshold[i], in units of 16B quads.  */

  unsigned char ts_rec_map[4][8];  /*!< Torus Reception Map.                    
              This array contains the rec fifo ID into which torus              
              packets are deposited that originate from                         
              - each hardware torus fifo (x+,x-,y+,y-,z+,z-) (0 through 5)      
              - high-priority hardware torus fifo (6)                           
              - a local transfer (7)                                            

              for each group (0 through 3).                                     
                                                                                
              For ts_rec_map[i][j],                                             
              i is the rec fifo group ID, as defined by (pid0,pid1) pair.       
              j identifies the hardware torus fifo (0-5 for                     
              x+,x-,y+,y-,z+,z- repectively), high-priority                     
              torus fifo (6), and local transfer (7).                           
              The value in each arrary element must be a global fifo ID         
              (between 0 and DMA_NUM_NORMAL_REC_FIFOS-1).  That is, the value   
              identifies the normal rec fifo that will receive packets          
              originating at i,j.                                               
              Note that the global fifo ID must be 0-7 for group 0,             
              8-15 for group 1, 16-23 for group 2, and 24-31 for group 3.       
                                                                                
              This affects DCRS 0xd60 to 0xd67 as defined by the following      
              - _BGP_DCR_rDMA_TS_REC_FIFO_MAP_G0_PID00_XY  (_BGP_DCR_DMA+0x60)  
              - _BGP_DCR_rDMA_TS_REC_FIFO_MAP_G0_PID00_ZHL (_BGP_DCR_DMA+0x61)  
              - _BGP_DCR_rDMA_TS_REC_FIFO_MAP_G0_PID01_XY  (_BGP_DCR_DMA+0x62)  
              - _BGP_DCR_rDMA_TS_REC_FIFO_MAP_G0_PID01_ZHL (_BGP_DCR_DMA+0x63)  
              - _BGP_DCR_rDMA_TS_REC_FIFO_MAP_G1_PID10_XY  (_BGP_DCR_DMA+0x64)  
              - _BGP_DCR_rDMA_TS_REC_FIFO_MAP_G1_PID10_ZHL (_BGP_DCR_DMA+0x65)  
              - _BGP_DCR_rDMA_TS_REC_FIFO_MAP_G1_PID11_XY  (_BGP_DCR_DMA+0x66)  
              - _BGP_DCR_rDMA_TS_REC_FIFO_MAP_G1_PID11_ZHL (_BGP_DCR_DMA+0x67)  
                                                                                
              e.g., for (pid0,pid1) = (0,1)                                     
              - ts_rec_map[1][0] = fifo id for torus x+ receiver packets          
              - ts_rec_map[1][1] = fifo id for torus x- receiver packets          
              - ts_rec_map[1][2] = fifo id for torus y+ receiver packets          
              - ts_rec_map[1][3] = fifo id for torus y- receiver packets          
              - ts_rec_map[1][4] = fifo id for torus z+ receiver packets          
              - ts_rec_map[1][5] = fifo id for torus z- receiver packets          
              - ts_rec_map[1][6] = fifo id for torus high priority packets        
              - ts_rec_map[1][7] = fifo id for local transfer packets         */

} DMA_RecFifoMap_t;


/*!
 * \brief DMA Reception Fifo Status Structure
 * 
 * Defines the DMA SRAM Status Area for Reception Fifos
 *
 */
typedef struct DMA_RecFifoStatus_t
{
  volatile unsigned  not_empty[2]; /*!< R bit mask, 1 bit/FIFO:                 
                                        Reception FIFO not empty status.        
                                        Not_empty[0] contains 1 bit for each    
                                        of the 32 normal fifos.                 
                                        Each bit corresponds to a               
                                        global fifo ID.                         
                                        Not_empty[1] :                          
                                        - bit 7  for group 0 header fifo,       
                                        - bit 15 for group 1 header fifo,       
                                        - bit 23 for group 2 header fifo,       
                                        - bit 31 for group 3 header fifo.     */

  volatile unsigned  available[2]; /*!< R bitmask, 1 bit/FIFO:                  
                                        Reception FIFO available status.        
                                        Bits are as above for available[0]      
                                        and available[1].                     */

  volatile unsigned  threshold_crossed[2]; /*!< R bitmask, 1 bit/FIFO:          
                                                Threshold crossed status.       
                                                Bits are as above for           
                                                threshhold_crossed[0] and       
                                                threshhold_crossed[1].        */
   
  volatile unsigned  clear_threshold_crossed[2]; /*!< W bitmask, 1 bit/FIFO:    
                                        Clear Threshold Crossed Status.         
                                        Bits are as above for                   
                                        clear_threshold_crossed[0] and          
                                        clear_threshold_crossed[1].           */
}
DMA_RecFifoStatus_t;


/*!
 * \brief Returns the word number that the specified reception fifo is in
 *
 * \param[in]  global_rec_fifo_id  The global ID of the reception fifo
 *                                 (0 to DMA_NUM_REC_FIFOS-1).
 * 
 * \return The number of the word that the specified fifo is in (0 or 1).
 *
 * Used as an index in the "not_empty", "available", "threshold_crossed", and
 * "clear_threshold_crossed" fields of the DMA_RecFifoStatus_t structure.
 *
 */
#define DMA_REC_FIFO_GROUP_WORD_ID(global_rec_fifo_id) \
                                      ( ((global_rec_fifo_id)>>5) & _BN(31) )


/*!
 * \brief Reception DMA Fifo Structure
 *
 * This structure contains a software DMA fifo structure (defined in DMA_Fifo.h)
 * and other fields that are specific to a reception fifo used by software.
 */
typedef struct DMA_RecFifo_t
{
  DMA_Fifo_t         dma_fifo;       /*!< Common software fifo structure      */

  unsigned char      global_fifo_id; /*!< Global fifo ID:                       
                                          - 0 to DMA_NUM_NORMAL_REC_FIFOS-1     
                                            for normal fifos,                   
                                          - 32-35 for header fifos.           */
  /*!
   * \note The following field contains info about the fifo that reflects the 
   *       DCR values configuring the fifo.
   */

  unsigned char type;                /*!< 0 or 1 for type of fifo.            */

  /*!
   * \note The following field is used by the reception fifo polling functions.
   *       It counts the number of packets processed since the fifo's hardware
   *       head was last updated.  When DMA_MAX_NUM_PACKETS_BEFORE_MOVING_HEAD
   *       is reached, the hardware head is moved and this counter is reset
   *       to zero.  This helps to minimize the number of times the hardware
   *       head is updated, which can be an expensive operation.
   */
  unsigned int  num_packets_processed_since_moving_fifo_head; /*!< Tracks when
								to move the 
								fifo head.    */
} 
DMA_RecFifo_t;


/*!
 * \brief DMA Reception Fifo Group Structure
 *
 * This structure defines a DMA Reception Fifo Group.  It points to a
 * Reception Fifo Status structure, and contains DMA_NUM_REC_FIFOS_PER_GROUP
 * Reception Fifo structures.
 *
 * It is returned from the DMA_RecFifoGetFifoGroup system call wrapper function
 * defined in this header file.  This same structure must be used by all users 
 * of reception fifos in this group because the fifos will contain packets 
 * destined for these different users, and this structure contains shadows of
 * the hardware fifos in the DMA SRAM that must be maintained as the fifos 
 * change.  This common structure is located in static storage 
 * declared in DMA_RecFifo.c.
 *
 */
typedef struct DMA_RecFifoGroup_t
{
  int group_id;         /*!< Group ID (0 through DMA_NUM_REC_FIFO_GROUPS-1)   */

  int num_normal_fifos; /*!< Number of normal fifos used in this group          
                             (0 through DMA_NUM_NORMAL_REC_FIFOS_PER_GROUP)   */

  int num_hdr_fifos;    /*!< Number of header fifos used in this group          
                             - 0 - headers not being saved,                       
                             - 1 - headers being saved.                       */

  unsigned  mask;       /*!< All reads to the status for this group are         
                             masked by this, so you only see results for        
                             this group.                                        
                             - Group 0: 0xFF000000                                
                             - Group 1: 0x00FF0000                                
                             - Group 2: 0x0000FF00                                
                             - Group 3: 0x000000FF                            */

  void *unused1;        /*!< Unused space                                     */
 
  DMA_RecFifoStatus_t *status_ptr; /*!< Pointer to the status, in DMA SRAM.   */

  DMA_RecFifo_t        fifos[DMA_NUM_REC_FIFOS_PER_GROUP]; /*!< Rec Fifos.      
                  Indexes 0 through DMA_NUM_NORMAL_REC_FIFOS_PER_GROUP-1 are    
                  the normal fifos in the group.                                
                  Index DMA_HEADER_REC_FIFO_ID is the header fifo in the        
                  group.                                                        
                  Note:  fifos[0] may not be local fifo number 0 in the group.*/
}
ALIGN_L1D_CACHE  DMA_RecFifoGroup_t;


/*!
 * \brief DMA Reception Fifo Receive Function Prototype
 *
 * Receive functions must be registered through the 
 * DMA_RecFifoRegisterRecvFunction interface, which assigns a registration ID
 * to the function.  When the polling functions process a packet in a 
 * reception fifo, the appropriate receive function for that packet is
 * called with a pointer to the packet header, pointer to the payload, and
 * length of the payload.  The packet header is always 16 bytes of
 * contiguous storage, in the fifo.  Because the fifo is a circular buffer,
 * the payload of a packet may wrap from the end of the fifo to the beginning.
 * For large fifos, this happens infrequently.  To make it easier for 
 * user/messaging code, the poll function will always return a starting payload
 * address and number of bytes so that the receive function can treat the packet
 * as contiguous storage in memory.  If the packet does not wrap, the starting
 * payload address will be a pointer to the appropriate address in the fifo.
 * If the packet does wrap, the poll function will copy bytes from the fifo to 
 * a contiguous buffer (on the stack) and call the receive function with a 
 * payload pointer pointing to this temporary buffer.  In either case, when the
 * receive function returns, user code cannot assume that the payload buffer is 
 * permanent, i.e., after return, it may be overwritten by either the DMA or 
 * the poll function.  To keep a copy of the packet, the receive function would 
 * have to copy it to some other location.
 *
 * \param[in]  f_ptr           Pointer to the reception fifo.
 * \param[in]  packet_ptr      Pointer to the packet header (== va_head).
 *                             This is quad-aligned for optimized copying.
 * \param[in]  recv_func_parm  Pointer to storage specific to this receive
 *                             function.  This pointer was specified when the 
 *                             receive function was registered with the kernel,
 *                             and is passed to the receive function
 *                             unchanged.
 * \param[in]  payload_ptr     Pointer to the beginning of the payload.
 *                             This is quad-aligned for optimized copying.
 * \param[in]  payload_bytes   Number of bytes in the payload
 *
 * \retval   0                 No errors found while processing the payload.
 * \retval   negative_number   Errors found while processing the payload.
 */
typedef int  (*DMA_RecFifoRecvFunction_t)(
					  DMA_RecFifo_t      *f_ptr,
					  DMA_PacketHeader_t *packet_ptr,
					  void               *recv_func_parm,
					  char               *payload_ptr,
					  int                 payload_bytes
					 );


/*!
 * \brief DMA Reception Fifo Default Error Receive Function
 *
 * This is the default function that will handle packets having an 
 * unregistered registration ID.
 * 
 * Receive functions must be registered through the 
 * DMA_RecFifoRegisterRecvFunction interface, which assigns a registration ID
 * to the function.  When the polling functions process a packet in a 
 * reception fifo that has a registration ID that does not have a corresponding
 * receive function registered, this error receive function is
 * called with a pointer to the packet header, pointer to the payload, and
 * length of the payload.  The packet header is always be 16 bytes of
 * contiguous storage, in the fifo.  Because the fifo is a circular buffer,
 * the payload of a packet may wrap from the end of the fifo to the beginning.
 * For large fifos, this happens infrequently.  To make it easier for 
 * user/messaging code, the poll function will always return a starting payload
 * address and number of bytes so that the receive function can treat the packet
 * as contiguous storage in memory.  If the packet does not wrap, the starting
 * payload address will be a pointer to the appropriate address in the fifo.
 * If the packet does wrap, the poll function will copy bytes from the fifo to 
 * a contiguous buffer (on the stack) and call the receive function with a 
 * payload pointer pointing to this temporary buffer.  In either case, when the
 * receive function returns, user code cannot assume that the payload buffer is 
 * permanent, i.e., after return, it may be overwritten by either the DMA or 
 * the poll function.  To keep a copy of the packet, the receive function would 
 * have to copy it to some other location.
 *
 * \param[in]  f_ptr           Pointer to the reception fifo.
 * \param[in]  packet_ptr      Pointer to the packet header (== va_head).
 *                             This is quad-aligned for optimized copying.
 * \param[in]  recv_func_parm  Pointer to storage specific to this receive
 *                             function.  This pointer was specified when the 
 *                             receive function was registered with the kernel,
 *                             and is passed to the receive function
 *                             unchanged.
 * \param[in]  payload_ptr     Pointer to the beginning of the payload.
 *                             This is quad-aligned for optimized copying.
 * \param[in]  payload_bytes   Number of bytes in the payload
 *
 * \retval  -1  An unregistered packet was just processed.  This is considered
 *              an error.
 */
int  DMA_RecFifoDefaultErrorRecvFunction(
					 DMA_RecFifo_t      *f_ptr,
					 DMA_PacketHeader_t *packet_ptr,
					 void               *recv_func_parm,
					 char               *payload_ptr,
					 int                 payload_bytes
					);


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
 *                       enum located in bgp/arch/include/common/bgp_ras.h.
 *                       _bgp_err_dma_rfifo_map_twice means the mapping has
 *                       already been set.
 *
 * \note  This function should be called once per job, after DMA_ResetRelease().
 *        It may be called by any core, but once a core has called it, other
 *        calls by that same core or any other core will fail.
 *
 * \note  During job init, the kernel sets up the DCR clear masks for each 
 *        reception fifo group (DCRs 0xD68 - 0xD6C) such that a write to clear 
 *        a fifo in group g only clears group g.
 * 
 */
__INLINE__ int DMA_RecFifoSetMap( 
				 DMA_RecFifoMap_t * rec_map
				)
{
	int rc;
	rc = Kernel_RecFifoSetMap((uint32_t*)rec_map);
	return rc;
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
__INLINE__ int DMA_RecFifoGetMap(  
				 DMA_RecFifoMap_t  *rec_map
				)
{
	int rc;
	rc = Kernel_RecFifoGetMap((uint32_t*)rec_map);
	return rc;
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
 */
DMA_RecFifoGroup_t *
DMA_RecFifoGetFifoGroup( 
			int                       grp,
			int                       target,
			Kernel_CommThreadHandler  normal_handler,
			void                     *normal_handler_parm,
			Kernel_CommThreadHandler  header_handler,
			void                     *header_handler_parm,
			Kernel_InterruptGroup_t   interruptGroup
		       );




/*
 * -----------------------------------------------------------------------------
 * Calls to access the Fifo, given a reception fifo structure
 * -----------------------------------------------------------------------------
 */




/*!
 * \brief Increment DMA Reception Fifo Head
 *
 * Increment a DMA reception fifo's "head", given a reception fifo structure
 *
 * \param[in]  f_ptr  Pointer to the reception fifo structure
 * \param[in]  incr   The number of quads (16 byte units) to increment the 
 *                    head pointer by.
 *
 * \return  None
 *
 * \post va_head is set in both the hardware and software fifo structures,
 *       and the fifo free space is recalculated.
 *
 */
__INLINE__ void DMA_RecFifoIncrementHead(  
					 DMA_RecFifo_t *f_ptr,
					 unsigned int   incr 
					)
{
  SPI_assert( f_ptr != NULL );

  void *va_head = DMA_FifoGetHeadNoFreeSpaceUpdate( &f_ptr->dma_fifo );

  void *va_end  = DMA_FifoGetEndFromShadow( &f_ptr->dma_fifo );

  unsigned int incr_bytes = incr << 4;

  unsigned int bytes_to_end = (unsigned)va_end - (unsigned)va_head;

  /* 
   * Note:  The following check must be >= instead of just >.  We never want
   *        the head to be equal to the end so we can always copy a quad
   *        from the head, safely.
   */
  if ( incr_bytes >= bytes_to_end )
    {
      va_head = (char *)
	          ( (unsigned)DMA_FifoGetStartFromShadow( &f_ptr->dma_fifo ) + 
		    ( incr_bytes - bytes_to_end ) );
    }
  else
    {
      va_head = (char *)( (unsigned)va_head + incr_bytes );
    }

  /* Set the head and update the fifo's free space */
  DMA_FifoSetHead( &f_ptr->dma_fifo,
		   va_head );
}

                  
/*!
 * \brief Get the "Not Empty" Status of a Reception Fifo Group
 *
 * Get the "Not Empty" status of the reception fifos that are being used in the
 * specified "not empty" word.
 *
 * \param[in]  fg_ptr     Pointer to the reception fifo group structure
 * \param[in]  word       The word (0 or 1) of the "not empty" status to be
 *                        returned.
 *
 * \retval  notEmptyStatus  A 32-bit value:
 *                          If "word" is 0, bit i is 1 if normal rec fifo i is
 *                          in use and is not empty.
 *                          If "word" is 1, bits 7, 15, 23, 31 are 1 if header
 *                          rec fifos for groups 1, 2, 3, 4 respectively are in
 *                          use and not empty.
 *
 */
__INLINE__ unsigned DMA_RecFifoGetNotEmpty( 
					   DMA_RecFifoGroup_t *fg_ptr,
					   int                 word
					  ) 
{
  SPI_assert( fg_ptr != NULL );
  SPI_assert( fg_ptr->status_ptr != NULL );
  SPI_assert( (word == 0) || (word == 1) );

  //  printf("RecFifoGetNotEmpty: group=%d, status addr=0x%08x, not_empty=0x%08x, mask=0x%08x, RecFifoHwAddr=0x%08x, RecFifo4PaTail=0x%08x, PaHead=0x%08x\n",
  //	 fg_ptr->group_id, (unsigned)(&(fg_ptr->status_ptr->not_empty[word])),
  //	 fg_ptr->status_ptr->not_empty[word], fg_ptr->mask,
  //	 (unsigned)(fg_ptr->fifos[4].dma_fifo.fifo_hw_ptr),
  // 	 fg_ptr->fifos[4].dma_fifo.fifo_hw_ptr->pa_tail,
  //	 fg_ptr->fifos[4].dma_fifo.fifo_hw_ptr->pa_head);

  return ( fg_ptr->status_ptr->not_empty[word] & fg_ptr->mask );

}

                  
/*!
 * \brief Get the "Available" Status of a Reception Fifo Group
 *
 * Get the "available" status of the reception fifos that are being used in the
 * specified "available" word.
 *
 * \param[in]  fg_ptr     Pointer to the reception fifo group structure
 * \param[in]  word       The word (0 or 1) of the "available" status to be
 *                        returned.
 *
 * \retval  availableStatus  A 32-bit value:
 *                           If "word" is 0, bit i is 1 if normal rec fifo i is
 *                           in use and is available.
 *                           If "word" is 1, bits 7, 15, 23, 31 are 1 if header
 *                           rec fifos for groups 1, 2, 3, 4 respectively are in
 *                           use and available.
 *
 */
__INLINE__ unsigned DMA_RecFifoGetAvailable( 
					    DMA_RecFifoGroup_t *fg_ptr,
					    int                 word
					   ) 
{
  SPI_assert( fg_ptr != NULL );
  SPI_assert( fg_ptr->status_ptr != NULL );
  SPI_assert( (word == 0) || (word == 1) );

  return ( fg_ptr->status_ptr->available[word] & fg_ptr->mask );  
}


/*!
 * \brief Get the "Threshold Crossed" Status of a Reception Fifo Group
 *
 * Get the "threshold crossed" status of the reception fifos that are being 
 * used in the specified "threshold crossed" word.
 *
 * \param[in]  fg_ptr     Pointer to the reception fifo group structure
 * \param[in]  word       The word (0 or 1) of the "threshold crossed" status 
 *                        to be returned.
 *
 * \retval  thresholdCrossedStatus  A 32-bit value:
 *                           If "word" is 0, bit i is 1 if normal rec fifo i is
 *                           in use and its threshold has been crossed.
 *                           If "word" is 1, bits 7, 15, 23, 31 are 1 if header
 *                           rec fifos for groups 1, 2, 3, 4 respectively are in
 *                           use and their threshold has been crossed.
 *
 */
__INLINE__ unsigned DMA_RecFifoGetThresholdCrossed( 
					    DMA_RecFifoGroup_t *fg_ptr,
					    int                 word
					   ) 
{
  SPI_assert( fg_ptr != NULL );
  SPI_assert( fg_ptr->status_ptr != NULL );
  SPI_assert( (word == 0) || (word == 1) );

  return ( fg_ptr->status_ptr->threshold_crossed[word] & fg_ptr->mask );  
}


/*!
 * \brief Set the "Clear Threshold Crossed" Status of Specified Fifos in a 
 *        Reception Fifo Group
 *
 * Set the "clear threshold crossed" status of the specified reception fifos 
 * in the specified "clear threshold crossed" word.
 *
 * \param[in]  fg_ptr     Pointer to the reception fifo group structure
 * \param[in]  clr        32-bit value, specifying which fifos are to have
 *                        their "clear threshold crossed" status set.
 *                        If "word" is 0, bit i is 1 if normal rec fifo i is
 *                        to have its "clear threshold crossed" status set.
 *                        If "word" is 1, one of bits 7, 15, 23, 31 is 1 if 
 *                        header fifo for group 1, 2, 3, 4 respectively is to
 *                        have its "clear threshold crossed" status set.
 *                        Fifos that are not in the group will not have their
 *                        status set.
 * \param[in]  word       The word (0 or 1) of the "clear threshold crossed" 
 *                        status to be set.
 *
 * \return  None
 *
 * \note This function does an MBAR after setting the status to ensure the
 *       writes have been accepted by the memory system before allowing other
 *       memory accesses to to occur.
*/
__INLINE__ void DMA_RecFifoSetClearThresholdCrossed( 
					    DMA_RecFifoGroup_t *fg_ptr,
					    unsigned int        clr,
					    int                 word
					   ) 
{
  SPI_assert( fg_ptr != NULL );
  SPI_assert( fg_ptr->status_ptr != NULL );
  SPI_assert( (word == 0) || (word == 1) );

  fg_ptr->status_ptr->clear_threshold_crossed[word] = clr & fg_ptr->mask; 

  _bgp_mbar();

}


/*
 * -----------------------------------------------------------------------------
 * Calls to access the Fifo, given a fifo group and a fifo ID
 * -----------------------------------------------------------------------------
 */




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
 * \retval  0            Successful.
 * \retval  error_value  An error value defined in the _BGP_RAS_DMA_ErrCodes 
 *                       enum located in bgp/arch/include/common/bgp_ras.h.
 *                       _bgp_err_dma_rfifo_map_twice means this fifo has
 *                       already been initialized
 *
 */
__INLINE__ int DMA_RecFifoInitById(
				   DMA_RecFifoGroup_t *fg_ptr,
				   int                 fifo_id,
				   void               *va_start,
				   void               *va_head,
				   void               *va_end
				  )
{
  int rc;

  SPI_assert( fg_ptr != NULL );
  SPI_assert( (fifo_id >= 0) && (fifo_id < DMA_NUM_REC_FIFOS_PER_GROUP) );
  SPI_assert( ( (uint32_t) va_end - (uint32_t)va_start ) >=
	                                       DMA_MIN_REC_FIFO_SIZE_IN_BYTES );

  /* 
   * Initialize the fifo by invoking a system call.
   */

  rc = Kernel_RecFifoInitById( 
			      (uint32_t*)fg_ptr,
			      fifo_id,
			      va_start,
			      va_head,
			      va_end);

  return rc;
}


/*!
 * \brief Get DMA RecFifo Start Pointer from the Shadow Using a Fifo Id
 *
 * Get a DMA reception fifo's start pointer, given a fifo group and fifo id.
 *
 * \param[in]  fg_ptr   Pointer to the fifo group structure
 * \param[in]  fifo_id  Id of the fifo within the group
 *                      (0 to DMA_NUM_REC_FIFOS_PER_GROUP-1).
 *
 * \retval  va_start  The virtual address of the start of the fifo
 *
 */
__INLINE__ void * DMA_RecFifoGetStartById(  
					  DMA_RecFifoGroup_t *fg_ptr, 
					  int                 fifo_id
					 )
{
  SPI_assert( (fifo_id >= 0) && (fifo_id < DMA_NUM_REC_FIFOS_PER_GROUP) );
  SPI_assert( fg_ptr != NULL );

  return DMA_FifoGetStartFromShadow( &fg_ptr->fifos[fifo_id].dma_fifo );
}


/*!
 * \brief Get DMA RecFifo Head Pointer Using a Fifo Id
 *
 * Get a DMA reception fifo's head pointer, given a fifo group and fifo id.
 *
 * \param[in]  fg_ptr   Pointer to the fifo group structure
 * \param[in]  fifo_id  Id of the fifo within the group
 *                      (0 to DMA_NUM_REC_FIFOS_PER_GROUP-1).
 *
 * \retval  va_head  The virtual address of the head of the fifo
 *
 */
__INLINE__ void * DMA_RecFifoGetHeadById(  
					 DMA_RecFifoGroup_t *fg_ptr,
					 int                 fifo_id
					)
{
  SPI_assert( (fifo_id >= 0) && (fifo_id < DMA_NUM_REC_FIFOS_PER_GROUP) );
  SPI_assert( fg_ptr != NULL );

  return DMA_FifoGetHead( &fg_ptr->fifos[fifo_id].dma_fifo );
} 


/*!
 * \brief Get DMA RecFifo Tail Pointer Using a Fifo Id
 *
 * Get a DMA reception fifo's tail pointer, given a fifo group and fifo id.
 *
 * \param[in]  fg_ptr   Pointer to the fifo group structure
 * \param[in]  fifo_id  Id of the fifo within the group
 *                      (0 to DMA_NUM_REC_FIFOS_PER_GROUP-1).
 *
 * \retval  va_tail  The virtual address of the tail of the fifo
 *
 */
__INLINE__ void *DMA_RecFifoGetTailById(  
					DMA_RecFifoGroup_t *fg_ptr,
					int                 fifo_id
				       )
{
  SPI_assert( (fifo_id >= 0) && (fifo_id < DMA_NUM_REC_FIFOS_PER_GROUP) );
  SPI_assert( fg_ptr != NULL );

  return   DMA_FifoGetTail( &fg_ptr->fifos[fifo_id].dma_fifo );
} 


/*!
 * \brief Get DMA RecFifo End Pointer from the Shadow Using a Fifo Id
 *
 * Get a DMA reception fifo's end pointer, given a fifo group and fifo id.
 *
 * \param[in]  fg_ptr   Pointer to the fifo group structure
 * \param[in]  fifo_id  Id of the fifo within the group
 *                      (0 to DMA_NUM_REC_FIFOS_PER_GROUP-1).
 *
 * \retval  va_end  The virtual address of the end of the fifo
 *
 */
__INLINE__ void *DMA_RecFifoGetEndById(  
				       DMA_RecFifoGroup_t *fg_ptr,
				       int                 fifo_id
				      )
{
  SPI_assert( (fifo_id >= 0) && (fifo_id < DMA_NUM_REC_FIFOS_PER_GROUP) );
  SPI_assert( fg_ptr != NULL );

  return   DMA_FifoGetEndFromShadow( &fg_ptr->fifos[fifo_id].dma_fifo );
} 


/*!
 * \brief Get DMA RecFifo Size Using a Fifo Id
 *
 * Get a DMA reception fifo's size, given a fifo group and fifo id.
 *
 * \param[in]  fg_ptr   Pointer to the fifo group structure
 * \param[in]  fifo_id  Id of the fifo within the group
 *                      (0 to DMA_NUM_REC_FIFOS_PER_GROUP-1).
 *
 * \retval  size   The size of the DMA fifo, in units of 16B quads.
 *
 */
__INLINE__ unsigned int DMA_RecFifoGetSizeById(  
					       DMA_RecFifoGroup_t *fg_ptr,
					       int                 fifo_id
					      )
{
  SPI_assert( (fifo_id >= 0) && (fifo_id < DMA_NUM_REC_FIFOS_PER_GROUP) );
  SPI_assert( fg_ptr != NULL );

  return   DMA_FifoGetSize( &fg_ptr->fifos[fifo_id].dma_fifo );
}


/*!
 * \brief Get DMA RecFifo Free Space Using a Fifo Id
 *
 * Get a DMA reception fifo's free space, given a fifo group and fifo id.
 *
 * \param[in]  fg_ptr   Pointer to the fifo group structure
 * \param[in]  fifo_id  Id of the fifo within the group
 *                      (0 to DMA_NUM_REC_FIFOS_PER_GROUP-1).
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
 */
__INLINE__ unsigned int  DMA_RecFifoGetFreeSpaceById(  
					 DMA_RecFifoGroup_t *fg_ptr,
					 int                 fifo_id,
					 unsigned int        read_head, 
					 unsigned int        read_tail
					)
{
  SPI_assert( (fifo_id >= 0) && (fifo_id < DMA_NUM_REC_FIFOS_PER_GROUP) );
  SPI_assert( fg_ptr != NULL );

  return DMA_FifoGetFreeSpace( &fg_ptr->fifos[fifo_id].dma_fifo,
			       read_head, 
			       read_tail );
} 


/*!
 * \brief Set DMA RecFifo Head Pointer Using a Fifo Id
 *
 * Set a DMA reception fifo's head pointer, given a fifo group and fifo id.
 *
 * \param[in]  fg_ptr   Pointer to the fifo group structure
 * \param[in]  fifo_id  Id of the fifo within the group
 *                      (0 to DMA_NUM_REC_FIFOS_PER_GROUP-1).
 * \param[in]  va_head  The virtual address of the head of the fifo.
 *
 * \return  None
 *
 */
__INLINE__ void DMA_RecFifoSetHeadById(  
				       DMA_RecFifoGroup_t *fg_ptr,
				       int                 fifo_id,
				       void               *va_head
				       )
{
  SPI_assert( (fifo_id >= 0) && (fifo_id < DMA_NUM_REC_FIFOS_PER_GROUP) );
  SPI_assert( fg_ptr != NULL );

  DMA_FifoSetHead( &fg_ptr->fifos[fifo_id].dma_fifo,
		   va_head);
}


/*!
 * \brief Set DMA RecFifo Tail Pointer Using a Fifo Id
 *
 * Set a DMA reception fifo's tail pointer, given a fifo group and fifo id.
 *
 * \param[in]  fg_ptr   Pointer to the fifo group structure
 * \param[in]  fifo_id  Id of the fifo within the group
 *                      (0 to DMA_NUM_REC_FIFOS_PER_GROUP-1).
 * \param[in]  va_tail  The virtual address of the tail of the fifo.
 *
 * \return  None
 *
 */
__INLINE__ void DMA_RecFifoSetTailById(  
				       DMA_RecFifoGroup_t *fg_ptr,
				       int                 fifo_id, 
				       void               *va_tail 
				      )
{
  SPI_assert( (fifo_id >= 0) && (fifo_id < DMA_NUM_REC_FIFOS_PER_GROUP) );
  SPI_assert( fg_ptr != NULL );

  DMA_FifoSetTail( &fg_ptr->fifos[fifo_id].dma_fifo,
		   va_tail);  
}


/*!
 * \brief Increment DMA RecFifo Head Pointer Using a Fifo Id
 *
 * Increment a DMA reception fifo's head pointer, given a fifo group and 
 * fifo id.
 *
 * \param[in]  fg_ptr   Pointer to the fifo group structure
 * \param[in]  fifo_id  Id of the fifo within the group
 *                      (0 to DMA_NUM_REC_FIFOS_PER_GROUP-1).
 * \param[in]  incr     The number of quads (16 byte units) to increment the 
 *                      head pointer by.
 *
 * \return  None
 *
*/
__INLINE__ void DMA_RecFifoIncrementHeadById( 
					     DMA_RecFifoGroup_t *fg_ptr,
					     int                 fifo_id, 
					     unsigned int        incr
					    )
{
  SPI_assert( (fifo_id >= 0) && (fifo_id < DMA_NUM_REC_FIFOS_PER_GROUP) );
  SPI_assert( fg_ptr != NULL );

  DMA_RecFifoIncrementHead( &fg_ptr->fifos[fifo_id],
			    incr);
} 


/*!
 * \brief Get DMA RecFifo Not Empty Status Using a Fifo Id
 *
 * Get a DMA reception fifo's not empty status, given a fifo group and 
 * fifo id.
 *
 * \param[in]  fg_ptr   Pointer to the fifo group structure
 * \param[in]  fifo_id  Id of the fifo within the group
 *                      (0 to DMA_NUM_REC_FIFOS_PER_GROUP-1).
 *
 * \retval  0         The specified fifo is empty.
 *          non-zero  The specified fifo is not empty.
 *
 */
__INLINE__ unsigned DMA_RecFifoGetNotEmptyById( 
					DMA_RecFifoGroup_t *fg_ptr,
					int                 fifo_id
				       ) 
{
  SPI_assert( (fifo_id >= 0) && (fifo_id < DMA_NUM_REC_FIFOS_PER_GROUP) );
  SPI_assert( fg_ptr != NULL );

  int word = DMA_REC_FIFO_GROUP_WORD_ID(fg_ptr->fifos[fifo_id].global_fifo_id);

  unsigned status; 
  status = DMA_RecFifoGetNotEmpty(fg_ptr,
				  word);
  if ( word == 0 ) 
    {
      /* If normal fifo, mask with the correct bit number. */
      status = status & _BN(fg_ptr->fifos[fifo_id].global_fifo_id); 
    }
  /* For header fifo, don't need additional mask because the status word was
   * already masked by the 8 bits for this group, leaving only the 1 bit for
   * the header fifo.
   */

  return status;

}


/*!
 * \brief Get DMA RecFifo Available Status Using a Fifo Id
 *
 * Get a DMA reception fifo's available status, given a fifo group and 
 * fifo id.
 *
 * \param[in]  fg_ptr   Pointer to the fifo group structure
 * \param[in]  fifo_id  Id of the fifo within the group
 *                      (0 to DMA_NUM_REC_FIFOS_PER_GROUP-1).
 *
 * \retval  0         The specified fifo is not available.
 *          non-zero  The specified fifo is available.
 *
 */
__INLINE__ unsigned DMA_RecFifoGetAvailableById( 
					DMA_RecFifoGroup_t *fg_ptr,
					int                 fifo_id
				       ) 
{
  SPI_assert( (fifo_id >= 0) && (fifo_id < DMA_NUM_REC_FIFOS_PER_GROUP) );
  SPI_assert( fg_ptr != NULL );

  int word = DMA_REC_FIFO_GROUP_WORD_ID(fg_ptr->fifos[fifo_id].global_fifo_id);

  unsigned status; 
  status = DMA_RecFifoGetAvailable(fg_ptr,
				   word);
  if ( word == 0 ) 
    {
      /* If normal fifo, mask with the correct bit number. */
      status = status & _BN(fg_ptr->fifos[fifo_id].global_fifo_id); 
    }
  /* For header fifo, don't need additional mask because the status word was
   * already masked by the 8 bits for this group, leaving only the 1 bit for
   * the header fifo.
   */

  return status;

}


/*!
 * \brief Get DMA RecFifo Threshold Crossed Status Using a Fifo Id
 *
 * Get a DMA reception fifo's threshold crossed status, given a fifo group and 
 * fifo id.
 *
 * \param[in]  fg_ptr   Pointer to the fifo group structure
 * \param[in]  fifo_id  Id of the fifo within the group
 *                      (0 to DMA_NUM_REC_FIFOS_PER_GROUP-1).
 *
 * \retval  0         The specified fifo has not had its threshold crossed
 *          non-zero  The specified fifo has had its threshold crossed
 *
 */
__INLINE__ unsigned DMA_RecFifoGetThresholdCrossedById( 
					DMA_RecFifoGroup_t *fg_ptr,
					int                 fifo_id
				       ) 
{
  SPI_assert( (fifo_id >= 0) && (fifo_id < DMA_NUM_REC_FIFOS_PER_GROUP) );
  SPI_assert( fg_ptr != NULL );

  int word = DMA_REC_FIFO_GROUP_WORD_ID(fg_ptr->fifos[fifo_id].global_fifo_id);

  unsigned status; 
  status = DMA_RecFifoGetThresholdCrossed(fg_ptr,
					  word);
  if ( word == 0 ) 
    {
      /* If normal fifo, mask with the correct bit number. */
      status = status & _BN(fg_ptr->fifos[fifo_id].global_fifo_id); 
    }
  /* For header fifo, don't need additional mask because the status word was
   * already masked by the 8 bits for this group, leaving only the 1 bit for
   * the header fifo.
   */

  return status;

}


/*!
 * \brief Set DMA RecFifo Clear Threshold Crossed Status Using a Fifo Id
 *
 * Set a DMA reception fifo's "clear threshold crossed" status, given a fifo 
 * group and fifo id.
 *
 * \param[in]  fg_ptr   Pointer to the fifo group structure
 * \param[in]  fifo_id  Id of the fifo within the group
 *                      (0 to DMA_NUM_REC_FIFOS_PER_GROUP-1).
 *
 * \return  None
 *
 */
__INLINE__ void DMA_RecFifoSetClearThresholdCrossedById( 
					   DMA_RecFifoGroup_t *fg_ptr,
					   int                 fifo_id
					  )
{
  unsigned int clr;
  int          word;

  SPI_assert( (fifo_id >= 0) && (fifo_id < DMA_NUM_REC_FIFOS_PER_GROUP) );
  SPI_assert( fg_ptr != NULL );

  word = DMA_REC_FIFO_GROUP_WORD_ID(fg_ptr->fifos[fifo_id].global_fifo_id); 

  if ( word == 0 ) 
    {
      /* If normal fifo, mask with the correct bit number so we only specify the
       * bit corresponding to this normal fifo.
       * Note:  The fg_ptr->mask shouldn't be necessary, but it is a bit safer.
       */
      clr = ( _BN(fg_ptr->fifos[fifo_id].global_fifo_id) & fg_ptr->mask); 
    }
  else 
    {
      /* If header fifo, it is ok to just clear all of the mask bits for this
       * group, since only 1 bit is used inside the DMA.
       */
      clr = fg_ptr->mask; 
    }

  DMA_RecFifoSetClearThresholdCrossed(fg_ptr,
				      clr,
				      word); /* Write to the DMA SRAM */
}
 

/*!
 * \brief Register a Reception Fifo Receive Function
 *
 * Register a specified receive function to handle packets having a specific
 * "registration ID".  It returns a registration ID (0-255) that is to be used 
 * in the packet header Func_Id field, such that packets that arrive in a 
 * reception fifo will result in the corresponding receive function being called
 * when that fifo is processed by a polling or interrupt handler function.
 *
 * \param[in]  recv_func          Pointer to the receive function.
 * \param[in]  recv_func_parm     Arbitrary pointer to be passed to the
 *                                recv_func when it is called.
 * \param[in]  is_error_function  If non-zero, this is the receiver function
 *                                to be called if a packet contains an invalid
 *                                (unregistered) registration ID.  The return
 *                                value from this function is zero, indicating
 *                                success, not indicating a registration ID.
 *                                A default function is provided if one is not
 *                                registered.  If there is already a non-default
 *                                error receive function registered, -EBUSY is 
 *                                returned.
 * \param[in]  is_header_fifo     Indicates whether the fifo is normal or 
 *                                header:
 *                                - 0 is normal.  The return code is the
 *                                  registration ID.
 *                                - 1 is header.  The return code is 0,
 *                                  indicating success, because packets in
 *                                  header fifos are direct-put packets, and
 *                                  hence have no registration ID.
 *                                If there is already a header receive function
 *                                registered, -EBUSY is returned.
 *
 * If both is_error_function and is_header_fifo are 1, -EINVAL is returned.
 *
 * \retval   0            This is a registration ID if is_error_function=0 and
 *                        is_header_fifo=0.  Otherwise, it indicates success.
 *           1-255        This is a registration ID.  Successful.
 *           negative     Failure.  This is a negative errno value.
 */
int DMA_RecFifoRegisterRecvFunction( 
			        DMA_RecFifoRecvFunction_t  recv_func,
				void                      *recv_func_parm,
				int                        is_error_function,
				int                        is_header_fifo
			       );


/*! 
 * \brief Poll Normal Reception Fifos
 *
 * Poll the "normal" reception fifos in the specified fifo group, removing one 
 * packet after another from the fifos, dispatching the appropriate receive 
 * function for each packet, until one of the following occurs:
 * 1.  Total_packets packets are received
 * 2.  All the fifos are empty
 * 3.  A receive function returns a non-zero value
 * 4.  The last packet removed from a fifo has an invalid registration id.  The
 *     error receive function will have been called, but polling ends.
 *     The invalid packet is counted as a processed packet, and the return
 *     code from the error receive function is returned.
 *
 * Polling occurs in a round-robin fashion through the array of normal fifos in 
 * the group, beginning with array index starting_index. If a fifo has a packet,
 * the appropriate receive function is called.  Upon return, the packet is 
 * removed from the fifo (the fifo head is moved past the packet).
 *
 * After processing packets_per_fifo packets in a fifo (or emptying that fifo),
 * the next fifo in the group is processed.  When the last index in the fifo
 * array is processed, processing continues with the first fifo in the array.
 * Multiple loops through the array of fifos in the group may occur.
 *
 * The receive functions must be registered through the 
 * DMA_RecFifoRegisterRecvFunction interface.  The receive function is
 * called with a pointer to the packet header, pointer to the payload, and
 * length of the payload.  The packet header is always be 16 bytes of
 * contiguous storage, in the fifo.  Because the fifo is a circular buffer,
 * the payload of a packet may wrap from the end of the fifo to the beginning.
 * For large fifos, this happens infrequently.  To make it easier for 
 * user/messaging code, the poll function will always return a starting payload
 * address and number of bytes so that the receive function can treat the packet
 * as contiguous storage in memory.  If the packet does not wrap, the starting
 * payload address will be a pointer to the appropriate address in the fifo.
 * If the packet does wrap, the poll function will copy bytes from the fifo to 
 * a contiguous buffer (on the stack) and call the receive function with a 
 * payload pointer pointing to this temporary buffer.  In either case, when the
 * receive function returns, user code cannot assume that the payload buffer is 
 * permanent, i.e., after return, it may be overwritten by either the DMA or 
 * the poll function.  To keep a copy of the packet, the receive function would 
 * have to copy it to some other location.  The packet header and payload are
 * 16-byte aligned for optimized copying.
 *
 * \param[in]  total_packets     The maximum number of packets that will be
 *                               processed.
 * \param[in]  packets_per_fifo  The maximum number of packets that will be
 *                               processed in a given fifo before switching
 *                               to the next fifo.
 * \param[in]  starting_index    The fifos in the fifo group are maintained
 *                               in an array.  This is the array index of the
 *                               first fifo to be processed (0 through
 *                               DMA_NUM_NORMAL_REC_FIFOS_PER_GROUP-1).
 * \param[in]  num_empty_passes  The number of passes over the normal fifos 
 *                               while they are empty that this function 
 *                               should tolerate before giving up and 
 *                               returning.  This is an optimization
 *                               to catch late arriving packets.
 * \param[in]  not_empty_poll_delay  The number of pclks to delay between polls
 *                                   of the not-empty status when the fifos are
 *                                   empty.
 * \param[in]  fg_ptr            Pointer to the fifo group.
 * \param[out] next_fifo_index   Pointer to an int where the recommended 
 *                               starting_index for the next call is returned.
 *
 * \retval  num_packets_received  The number of packets received and processed.
 *                                next_fifo_index is set.
 * \retval  negative_value        The return code from the receive function that
 *                                caused polling to end.  next_fifo_index is
 *                                set.
 *
 * \pre  The caller is responsible for disabling interrupts before invoking this
 *       function.  
 * \todo By setting fg_ptr->interrupt_lock? or by calling
 *       the system call to disable a class of interrupts?
 *
 * \note next_fifo_index is set to the index of the fifo that had the last
 *       packet received if all packets_per_fifo packets were not received from
 *       that fifo.  However, if all packets_per_fifo packets were received
 *       from that fifo, the index of the next fifo will be returned.
 *
 */
int DMA_RecFifoPollNormalFifos(int                 total_packets,    
			       int                 packets_per_fifo,
			       int                 starting_index,
			       int                 num_empty_passes,
			       int                 not_empty_poll_delay,
			       DMA_RecFifoGroup_t *fg_ptr,
			       int                *next_fifo_index
			      );


/*! 
 * \brief Simple Poll Normal Reception Fifos
 *
 * Poll the "normal" reception fifos in the specified fifo group, removing one 
 * packet after another from the fifos, dispatching the appropriate receive 
 * function for each packet, until one of the following occurs:
 * 1.  All packets in all of the fifos have been received.
 * 2.  A receive function returns a non-zero value.
 * 3.  The last packet removed from a fifo has an invalid registration id.  The
 *     error receive function will have been called, but polling ends.
 *     The invalid packet is counted as a processed packet, and the return
 *     code from the error receive function is returned.
 * 4.  There have been fruitfulPollLimit polls attempted (summed across all 
 *     fifos).
 *
 * Polling occurs in a round-robin fashion through the array of normal fifos in 
 * the group.  If a fifo has a packet, the appropriate receive function is 
 * called.  Upon return, the packet is removed from the fifo (the fifo head is 
 * moved past the packet).
 *
 * After processing all of the packets in a fifo (or emptying that fifo),
 * the next fifo in the group is processed.  When the last index in the fifo
 * array is processed, processing continues with the first fifo in the array.
 * Multiple loops through the array of fifos in the group may occur until all
 * fifos are empty or fruitfulPollLimit polls have been completed.
 *
 * It is risky to set the fruitfulPollLimit to zero, allowing this function to
 * poll indefinitely as long as there are packets to be processed.  This may
 * starve the node in a scenario where other nodes send "polling" packets to
 * our node, and our node never gets a chance to do anything else except
 * process those polling packets.
 *
 * The receive functions must be registered through the 
 * DMA_RecFifoRegisterRecvFunction interface.  The receive function is
 * called with a pointer to the packet header, pointer to the payload, and
 * length of the payload.  The packet header is always be 16 bytes of
 * contiguous storage, in the fifo.  Because the fifo is a circular buffer,
 * the payload of a packet may wrap from the end of the fifo to the beginning.
 * For large fifos, this happens infrequently.  To make it easier for 
 * user/messaging code, the poll function will always return a starting payload
 * address and number of bytes so that the receive function can treat the packet
 * as contiguous storage in memory.  If the packet does not wrap, the starting
 * payload address will be a pointer to the appropriate address in the fifo.
 * If the packet does wrap, the poll function will copy bytes from the fifo to 
 * a contiguous buffer (on the stack) and call the receive function with a 
 * payload pointer pointing to this temporary buffer.  In either case, when the
 * receive function returns, user code cannot assume that the payload buffer is 
 * permanent, i.e., after return, it may be overwritten by either the DMA or 
 * the poll function.  To keep a copy of the packet, the receive function would 
 * have to copy it to some other location.  The packet header and payload are
 * 16-byte aligned for optimized copying.
 *
 * \param[in]  fg_ptr             Pointer to the fifo group.
 * \param[in]  fruitfulPollLimit  The limit on the number of fruitful polls that
 *                                will be attempted (summed across all fifos).  
 *                                If the limit is reached, this function 
 *                                returns.  A value of zero means there is no
 *                                limit imposed.  A fruitful poll is one where
 *                                at least one packet has arrived in the fifo
 *                                since the last poll.
 *
 * \retval  num_packets_received  The number of packets received and processed.

 * \retval  negative_value        The return code from the receive function that
 *                                caused polling to end.
 *
 * \pre  The caller is responsible for disabling interrupts before invoking this
 *       function.  
 *
 */
int DMA_RecFifoSimplePollNormalFifos( DMA_RecFifoGroup_t *fg_ptr,
				      int                 fruitfulPollLimit);

/*! 
 * \brief Poll Normal Reception Fifo Given a Fifo Group and Fifo ID
 *
 * Poll the specified "normal" reception fifo in the specified fifo group, 
 * removing one packet after another from the fifo, dispatching the appropriate 
 * receive function for each packet, until one of the following occurs:
 * 1.  num_packets packets are received
 * 2.  The specified fifo is empty
 * 3.  A receive function returns a non-zero value
 * 4.  The last packet removed from the fifo has an invalid registration id. The
 *     error receive function will have been called, but polling ends.
 *     The invalid packet is counted as a processed packet, and the return
 *     code from the error receive function is returned.
 *
 * If the specified fifo has a packet, the appropriate receive function is 
 * called.  Upon return, the packet is removed from the fifo (the fifo head is 
 * moved past the packet).
 *
 * After processing num_packets packets in the fifo (or emptying that fifo),
 * the function returns the number of packets processed *
 * The receive functions must be registered through the 
 * DMA_RecFifoRegisterRecvFunction interface.  The receive function is
 * called with a pointer to the packet header, pointer to the payload, and
 * length of the payload.  The packet header is always be 16 bytes of
 * contiguous storage, in the fifo.  Because the fifo is a circular buffer,
 * the payload of a packet may wrap from the end of the fifo to the beginning.
 * For large fifos, this happens infrequently.  To make it easier for 
 * user/messaging code, the poll function will always return a starting payload
 * address and number of bytes so that the receive function can treat the packet
 * as contiguous storage in memory.  If the packet does not wrap, the starting
 * payload address will be a pointer to the appropriate address in the fifo.
 * If the packet does wrap, the poll function will copy bytes from the fifo to 
 * a contiguous buffer (on the stack) and call the receive function with a 
 * payload pointer pointing to this temporary buffer.  In either case, when the
 * receive function returns, user code cannot assume that the payload buffer is 
 * permanent, i.e., after return, it may be overwritten by either the DMA or 
 * the poll function.  To keep a copy of the packet, the receive function would 
 * have to copy it to some other location.  The packet header and payload are
 * 16-byte aligned for optimized copying.
 *
 * \param[in]  num_packets       The maximum number of packets that will be
 *                               processed.
 * \param[in]  fifo_id           The ID of the fifo to be polled.
 *                               (0 through 
 *                               DMA_NUM_NORMAL_REC_FIFOS_PER_GROUP-1).
 * \param[in]  num_empty_passes  The number of passes over the fifo
 *                               while it is empty that this function 
 *                               should tolerate before giving up and 
 *                               returning.  This is an optimization
 *                               to catch late arriving packets.
 * \param[in]  not_empty_poll_delay  The number of pclks to delay between polls
 *                                   of the not-empty status when the fifos are
 *                                   empty.
 * \param[in]  fg_ptr            Pointer to the fifo group.
 *
 * \retval  num_packets_received  The number of packets received and processed.
 * \retval  negative_value        The return code from the receive function that
 *                                caused polling to end.
 *
 * \pre  The caller is responsible for disabling interrupts before invoking this
 *       function.  
 * \todo By setting fg_ptr->interrupt_lock? or by calling
 *       the system call to disable a class of interrupts?
 *
 */
int DMA_RecFifoPollNormalFifoById( int                 num_packets,
				   int                 fifo_id,  
				   int                 num_empty_passes,
				   int                 not_empty_poll_delay,
				   DMA_RecFifoGroup_t *fg_ptr
				 );


/*! 
 * \brief Simple Poll Normal Reception Fifo Given a Fifo Group and Fifo ID
 *
 * Poll the specified "normal" reception fifo in the specified fifo group, 
 * removing one packet after another from the fifo, dispatching the appropriate 
 * receive function for each packet, until one of the following occurs:
 * 1.  All packets in the fifo have been received.
 * 2.  The specified fifo is empty.
 * 3.  A receive function returns a non-zero value.
 * 4.  The last packet removed from the fifo has an invalid registration id. The
 *     error receive function will have been called, but polling ends.
 *     The invalid packet is counted as a processed packet, and the return
 *     code from the error receive function is returned.
 * 5.  There have been fruitfulPollLimit polls attempted.
 *
 * If the specified fifo has a packet, the appropriate receive function is 
 * called.  Upon return, the packet is removed from the fifo (the fifo head is 
 * moved past the packet).
 *
 * After processing all of the packets in the fifo (emptying that fifo),
 * or the fruitfulPollLimit has been reached, the function returns the number 
 * of packets processed.
 *
 * It is risky to set the fruitfulPollLimit to zero, allowing this function to
 * poll indefinitely as long as there are packets to be processed.  This may
 * starve the node in a scenario where other nodes send "polling" packets to
 * our node, and our node never gets a chance to do anything else except
 * process those polling packets.
 *
 * The receive functions must be registered through the 
 * DMA_RecFifoRegisterRecvFunction interface.  The receive function is
 * called with a pointer to the packet header, pointer to the payload, and
 * length of the payload.  The packet header is always be 16 bytes of
 * contiguous storage, in the fifo.  Because the fifo is a circular buffer,
 * the payload of a packet may wrap from the end of the fifo to the beginning.
 * For large fifos, this happens infrequently.  To make it easier for 
 * user/messaging code, the poll function will always pass a starting payload
 * address and number of bytes so that the receive function can treat the packet
 * as contiguous storage in memory.  If the packet does not wrap, the starting
 * payload address will be a pointer to the appropriate address in the fifo.
 * If the packet does wrap, the poll function will copy bytes from the fifo to 
 * a contiguous buffer (on the stack) and call the receive function with a 
 * payload pointer pointing to this temporary buffer.  In either case, when the
 * receive function returns, user code cannot assume that the payload buffer is 
 * permanent, i.e., after return, it may be overwritten by either the DMA or 
 * the poll function.  To keep a copy of the packet, the receive function has 
 * to copy it to some other location.  The packet header and payload are
 * 16-byte aligned for optimized copying.
 *
 * \param[in]  fifo_id           The ID of the fifo to be polled.
 *                               (0 through 
 *                               DMA_NUM_NORMAL_REC_FIFOS_PER_GROUP-1).
 * \param[in]  fg_ptr            Pointer to the fifo group.
 * \param[in]  fruitfulPollLimit  The limit on the number of fruitful polls that
 *                                will be attempted.
 *                                If the limit is reached, this function 
 *                                returns.  A value of zero means there is no
 *                                limit imposed.  A fruitful poll is one where
 *                                at least one packet has arrived in the fifo
 *                                since the last poll.
 *
 * \retval  num_packets_received  The number of packets received and processed.
 * \retval  negative_value        The return code from the receive function that
 *                                caused polling to end.
 *
 * \pre  The caller is responsible for disabling interrupts before invoking this
 *       function.  
 *
 */
int DMA_RecFifoSimplePollNormalFifoById( int                 fifo_id,  
					 DMA_RecFifoGroup_t *fg_ptr,
					 int                 fruitfulPollLimit
				       );



/*! 
 * \brief Poll Header Reception Fifo Given a Fifo Group
 *
 * Poll the "header" reception fifo in the specified fifo group, 
 * removing one packet after another from the fifo, dispatching the appropriate 
 * receive function for each packet, until one of the following occurs:
 * 1.  Total_packets packets are received
 * 2.  The specified fifo is empty
 * 3.  A receive function returns a non-zero value
 *
 * If the header fifo has a packet, the appropriate receive function is 
 * called.  Upon return, the packet is removed from the fifo (the fifo head is 
 * moved past the packet).
 *
 * After processing num_packets packets in the fifo (or emptying that fifo),
 * the function returns the number of packets processed.
 *
 * The receive function must be registered through the 
 * DMA_RecFifoRegisterRecvFunction interface.  The receive function is
 * called with a pointer to the packet header. The packet header is always 
 * 16 bytes of contiguous storage, in the fifo.  When the
 * receive function returns, user code cannot assume that the buffer is 
 * permanent, i.e., after return, it may be overwritten by either the DMA or 
 * the poll function.  To keep a copy of the packet, the receive function would 
 * have to copy it to some other location.  The packet header is 16-byte aligned
 * for optimized copying.
 *
 * \param[in]  num_packets       The maximum number of packets that will be
 *                               processed.
 * \param[in]  num_empty_passes  The number of passes over the fifo
 *                               while it is empty that this function 
 *                               should tolerate before giving up and 
 *                               returning.  This is an optimization
 *                               to catch late arriving packets.
 * \param[in]  not_empty_poll_delay  The number of pclks to delay between polls
 *                                   of the not-empty status when the fifos are
 *                                   empty.
 * \param[in]  fg_ptr            Pointer to the fifo group.
 *
 * \retval  num_packets_received  The number of packets received and processed.
 * \retval  negative_value        The return code from the receive function that
 *                                caused polling to end.
 *
 * \pre  The caller is responsible for disabling interrupts before invoking this
 *       function.  
 *
 */
int DMA_RecFifoPollHeaderFifo( int                 num_packets, 
			       int                 num_empty_passes,
			       int                 not_empty_poll_delay,
			       DMA_RecFifoGroup_t *fg_ptr
			     );
   


__END_DECLS


#endif 
