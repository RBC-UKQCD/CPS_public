**--------------------------------------------------------------------
**  CVS keywords
**
**  $Author: mcneile $
**  $Date: 2003-06-22 13:34:46 $
**  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_wilson/qcdsp_nos/wfm_comm_backward.asm,v 1.1.1.1 2003-06-22 13:34:46 mcneile Exp $
**  $Id: wfm_comm_backward.asm,v 1.1.1.1 2003-06-22 13:34:46 mcneile Exp $
**  $Name: not supported by cvs2svn $
**  $Locker:  $
**  $Log: not supported by cvs2svn $
**  Revision 1.2  2001/06/19 18:13:01  anj
**  Serious ANSIfication.  Plus, degenerate double64.h files removed.
**  Next version will contain the new nga/include/double64.h.  Also,
**  Makefile.gnutests has been modified to work properly, propagating the
**  choice of C++ compiler and flags all the way down the directory tree.
**  The mpi_scu code has been added under phys/nga, and partially
**  plumbed in.
**
**  Everything has newer dates, due to the way in which this first alteration was handled.
**
**  Anj.
**
**  Revision 1.2  2001/05/25 06:16:07  cvs
**  Added CVS keywords to phys_v4_0_0_preCVS
**
**  $RCSfile: wfm_comm_backward.asm,v $
**  $Revision: 1.1.1.1 $
**  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_wilson/qcdsp_nos/wfm_comm_backward.asm,v $
**  $State: Exp $
**
**--------------------------------------------------------------------
****************************************************************************************
*---------------------------------------------------------------------------------------
*
* _wfm_comm_backward
*
* This routine performs the backward communications ( send to -mu, receive from +mu ).
* The maximum number of words before access to DRAM is relinquished 
* and the initialization of the SCU DMA registers with the communication 
* parameters must be set before this routine is called. This is done 
* by calling the routine scu_init.
*
* If STAND_ALONE = 1 it can be called from C as:
*
* wfm_comm_backward(float *af0,                 ; af for mu = 0 
*	            float *af1,                 ; af for mu = 1 
*	      	    float *af2,                 ; af for mu = 2 
*   	            float *af3,                 ; af for mu = 3 
*	            Wilson *wilson_p);          ; Wilson struct.
*
* mu = {0,1,2,3} <-> {x,y,z,t}
*
* If STAND_ALONE = 0 it can not be called from C. Instead the following is expected
* to be set up before the routine is called:
*
* 1) all the .ref in the STAND_ALONE = 0  case below must be defined
* 2) The following registers must be loaded as follows:
*    AR0 <-- Wilson pointer
*
*---------------------------------------------------------------------------------------
****************************************************************************************

	.version	30

	.include	"../../include/wilson.hasm"
	.include	"wfm_macros.hasm"
	.include	"wfm_nga_struct.hasm"

	.def	_wfm_comm_backward

*---------------------------------------------------------------------------------------
* References
*---------------------------------------------------------------------------------------
	.ref _wfm_copy_backward

	.ref	af0
	.ref	af1
	.ref	af2
	.ref	af3
	.ref	scu_b
	.ref	_wfm_wire_map
	.ref _gjp_local_axis
	
*---------------------------------------------------------------------------------------
* definitions
*---------------------------------------------------------------------------------------
FP		.set	AR3			; use FP for arguments only
WILSON_AD	.set	AR0
SCUB		.set	AR1
ATMP		.set	AR2
SCU		.set	AR4
WIRE_MAP	.set	AR5
LOCAL_AXIS	.set	AR6	
WIRE		.set	R1
AXIS_FLAG	.set	R2
TMP		.set	R0

*---------------------------------------------------------------------------------------
* allocations
*--------------------------------------------------------------------------------------
	.text
wmap	.word	_wfm_wire_map
laxis	.word	_gjp_local_axis
	
****************************************************************************************
* _wfm_comm_backward
****************************************************************************************
	.text
_wfm_comm_backward:

*---------------------------------------------------------------------------------------
* If STAND_ALONE = 1: C-callable
*---------------------------------------------------------------------------------------
	.if	STAND_ALONE = 1
	PUSH    FP
	LDI     SP,FP
*  Save all registers that are important to C
	PUSH    R4
	PUSH    R5
	RND	R6
	PUSHF   R6
	RND	R7
	PUSHF   R7
	PUSH    AR4
	PUSH    AR5
	PUSH    AR6
	PUSH    AR7
	PUSH    FP              ; Local frame pointer
	PUSH    DP
*  Load arguments from stack to registers
	LDP	@af0
	LDI     *-FP(2), 	TMP			; Address of destination
	STI	TMP,		@af0
	LDI     *-FP(3), 	TMP			; Address of destination
	STI	TMP,		@af1
	LDI     *-FP(4), 	TMP			; Address of destination
	STI	TMP,		@af2
	LDI     *-FP(5), 	TMP			; Address of destination
	STI	TMP,		@af3
	LDI     *-FP(6), 	WILSON_AD		; Wilson structure
	.endif

*---------------------------------------------------------------------------------------
* Routine starts here
*---------------------------------------------------------------------------------------
	PUSH    DP
	LDP	@wmap
	LDI	@wmap,				WIRE_MAP		; wire map 

	LDP	@laxis
	LDI	@laxis,				LOCAL_AXIS		; local axis 

	LDP	@scu_b
	LDI	@scu_b,				SCUB			; SCU base address

	LDI	*+SCUB(scu.poll_all), 		TMP			; debug

	LDP	@af0

* Send in the -mu direction
* Send base-addresses = 0 + af[mu];
	LDI	*+LOCAL_AXIS(0),		AXIS_FLAG
	BNZ	SNDX
	LDI	*+WIRE_MAP(1),			WIRE
	ADDI	WIRE,	SCUB,			SCU
	LDI	@af0,				TMP
	STI	TMP,				*+SCU(scu.send0)	; -X
SNDX:	NOP
	
	LDI	*+LOCAL_AXIS(1),		AXIS_FLAG
	BNZ	SNDY
	LDI	*+WIRE_MAP(3),			WIRE
	ADDI	WIRE,	SCUB,			SCU
	LDI	@af1,				TMP
	STI	TMP,				*+SCU(scu.send0)	; -Y
SNDY:	NOP

	LDI	*+LOCAL_AXIS(2),		AXIS_FLAG
	BNZ	SNDZ
	LDI	*+WIRE_MAP(5),			WIRE
	ADDI	WIRE,	SCUB,			SCU
	LDI	@af2,				TMP
	STI	TMP,				*+SCU(scu.send0)	; -Z
SNDZ:	NOP

	LDI	*+LOCAL_AXIS(3),		AXIS_FLAG
	BNZ	SNDT
	LDI	*+WIRE_MAP(7),			WIRE
	ADDI	WIRE,	SCUB,			SCU
	LDI	@af3,				TMP
	STI	TMP,				*+SCU(scu.send0)	; -T
SNDT:	NOP

* Receive from the +mu direction
* Receive base-addresses = wilson_p->comm_offset[mu] + af[mu];
	LDI	*+LOCAL_AXIS(0),			AXIS_FLAG
	BNZ	RCVX
	LDI	*+WIRE_MAP(0),				WIRE
	ADDI	WIRE,	SCUB,				SCU
	LDI	@af0,					TMP
	ADDI	*+WILSON_AD(Wilson.comm_offset+0),	TMP
	STI	TMP,				*+SCU(scu.receive0)	; +X
RCVX:	NOP
		
	LDI	*+LOCAL_AXIS(1),			AXIS_FLAG
	BNZ	RCVY
	LDI	*+WIRE_MAP(2),				WIRE
	ADDI	WIRE,	SCUB,				SCU
	LDI	@af1,					TMP
	ADDI	*+WILSON_AD(Wilson.comm_offset+1),	TMP
	STI	TMP,				*+SCU(scu.receive0)	; +Y
RCVY:	NOP

	LDI	*+LOCAL_AXIS(2),			AXIS_FLAG
	BNZ	RCVZ
	LDI	*+WIRE_MAP(4),				WIRE
	ADDI	WIRE,	SCUB,				SCU
	LDI	@af2,					TMP
	ADDI	*+WILSON_AD(Wilson.comm_offset+2),	TMP
	STI	TMP,				*+SCU(scu.receive0)	; +Z
RCVZ:	NOP

	LDI	*+LOCAL_AXIS(3),			AXIS_FLAG
	BNZ	RCVT
	LDI	*+WIRE_MAP(6),				WIRE
	ADDI	WIRE,	SCUB,				SCU
	LDI	@af3,					TMP
	ADDI	*+WILSON_AD(Wilson.comm_offset+3),	TMP
	STI 	TMP,				*+SCU(scu.receive0)	; +T
RCVT:	NOP


* If needed call the wfm_copy_backward routine to do
* any local communications by copying.
	LDI	*+LOCAL_AXIS(5),	AXIS_FLAG
	BZ	NCP
	PUSH	WILSON_AD
	LDI	@af3,	TMP
	PUSH	TMP
	LDI	@af2,	TMP
	PUSH	TMP
	LDI	@af1,	TMP
	PUSH	TMP
	LDI	@af0,	TMP
	PUSH	TMP
	CALL	_wfm_copy_backward
	SUBI	5,	SP
NCP:	NOP
			
	LDI	*+SCU(scu.poll_all), 		TMP			; debug

	POP	DP
*-------------------------------------------------------------------------
*
*  Restore the registers before returning to the C program                
*
	.if	STAND_ALONE = 1
	POP     DP
	POP     FP              ;  This is our frame pointer
	POP     AR7
	POP     AR6
	POP     AR5
	POP     AR4
	POPF    R7
	POPF    R6
	POP     R5
	POP     R4

	POP     FP              ;  This is our caller's frame pointer
	.endif

	RETS
	.end
