#include<config.h>
CPS_START_NAMESPACE
//-------------------------------------------------------------------
/*!\file
  \brief  Declarations of network communications routines for BG/L
*/

/*******************************************************************
*
* bgl_net.h
*
* The minimum amount of definitions needed to write a basic layer
* to send/receive torus packets.
*
* Several of the definitions and functions below were copied from 
* (in alphabetical order) G. Almasi, J. Castanos and M. Giampapa
*
*******************************************************************/

#ifndef INCLUDED_BGL_NET_H
#define INCLUDED_BGL_NET_H

CPS_END_NAMESPACE

#include <util/data_types.h>
#include <sys/bgl/bgl_sys_all.h>

CPS_START_NAMESPACE


//------------------------------------------------------------------
// typedefs
//------------------------------------------------------------------

typedef struct BGLCPSTorusStatus_tag {
  unsigned R0 :8;	/* chunks in receive FIFO 0 */
  unsigned R1 :8;	/* chunks in receive FIFO 1 */
  unsigned R2 :8;	/* chunks in receive FIFO 2 */
  unsigned R3 :8;	/* chunks in receive FIFO 3 */
  unsigned R4 :8;	/* chunks in receive FIFO 4 */
  unsigned R5 :8;	/* chunks in receive FIFO 5 */
  unsigned RH :8;	/* chunks in high priority receive FIFO */
  unsigned IN0:8;	/* chunks in send FIFO 0 */
  unsigned IN1:8;	/* chunks in send FIFO 1 */
  unsigned IN2:8;	/* chunks in send FIFO 2 */
  unsigned INH:8;	/* chunks in high priority send FIFO */
  unsigned dm1:8;	/* unused */
  unsigned dm2:32;	/* also unused */

} BGLCPSTorusStatus __attribute__((aligned(BGL_QUAD_ALIGNSIZE)));


//------------------------------------------------------------------
// Torus packet header 16 bytes (= quadword)
//------------------------------------------------------------------
typedef struct BGLCPSTorusPacketHardHeader_tag {
  unsigned type        :8;	/* BYTE 0 (also used as checksum hint) */
  unsigned hintXPlus   :1;	/* BYTE 1 bit 7 : X+ hint */
  unsigned hintXMinus  :1;	/* BYTE 1 bit 6 : X- hint */
  unsigned hintYPlus   :1;	/* BYTE 1 bit 5 : Y+ hint */
  unsigned hintYMinus  :1;	/* BYTE 1 bit 4 : Y- hint */
  unsigned hintZPlus   :1;	/* BYTE 1 bit 3 : Z+ hint */
  unsigned hintZMinus  :1;	/* BYTE 1 bit 2 : Z- hint */
  unsigned deposit     :1;	/* BYTE 1 bit 1 : deposit bit */
  unsigned fifoID      :1;	/* BYTE 1 bit 0 : dest. FIFO number */
  unsigned size        :3;	/* BYTE 2 bits 5 */
  unsigned avail       :2;	/* BYTE 2 bits 3 */
  unsigned dynamicRoute:1;	/* BYTE 2 bit 2 */
  unsigned vc          :2;	/* BYTE 2 bits 0 */
  unsigned destX       :8;	/* BYTE 3 : destination coordinate X */
  unsigned destY       :8;	/* BYTE 4 : destination coordinate Y */
  unsigned destZ       :8;	/* BYTE 5 : destination coordinate Z */
  unsigned seqno       :8;	/* BYTE 6 : sequence number */
  unsigned chksum      :8;	/* BYTE 7 : checksum */
} BGLCPSTorusPacketHardHeader __attribute__((aligned(BGL_QUAD_ALIGNSIZE)));

typedef struct _BGLCPSTorusPacketSoftHeader_tag
{
  unsigned val1;
  unsigned val2;
} BGLCPSTorusPacketSoftHeader __attribute__((aligned(BGL_DWORD_ALIGNSIZE)));

typedef struct _BGLCPSTorusPacketHeader_tag
{
  BGLCPSTorusPacketHardHeader hh;
  BGLCPSTorusPacketSoftHeader sh;
} BGLCPSTorusPacketHeader __attribute__((aligned(BGL_QUAD_ALIGNSIZE)));



//------------------------------------------------------------------
// Function defs
//------------------------------------------------------------------
void BGLCPSTorusPacketHeader_Init (BGLCPSTorusPacketHeader * h,
				   int hxp, int hxm, 
				   int hyp, int hym, 
				   int hzp, int hzm, 
				   int x, int y, int z, 
				   int fifo_id, int size);

void BGLCPSTorus_send(int dir, int size, IFloat *data);

void BGLCPSTorus_recv(int dir, int size, IFloat *data);

void BGLCPSTorus_send_spinor(int dir, BGLQuad *data);

void BGLCPSTorus_recv_spinor(int dir, BGLQuad *data);

void BGLCPSTorusPacketHeader_InitFill (void);

void BGLCPSGrid_InitFill (void);

void BGLCPSTorusMFifo_Init (void);

void BGLCPSVarious_Init (void);

void BGLCPSNet_Init (void);


#endif

CPS_END_NAMESPACE
