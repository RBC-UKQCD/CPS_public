**--------------------------------------------------------------------
**  CVS keywords
**
**  $Author: zs $
**  $Date: 2004-08-18 11:57:56 $
**  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_wilson/qcdsp/wfm_comm_forward.asm,v 1.4 2004-08-18 11:57:56 zs Exp $
**  $Id: wfm_comm_forward.asm,v 1.4 2004-08-18 11:57:56 zs Exp $
**  $Name: not supported by cvs2svn $
**  $Locker:  $
**  $Revision: 1.4 $
**  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_wilson/qcdsp/wfm_comm_forward.asm,v $
**  $State: Exp $
**
**--------------------------------------------------------------------
****************************************************************************************
*---------------------------------------------------------------------------------------
*
* _wfm_comm_forward
*
* This routine performs the forward communications ( send to +mu, receive from -mu ).
* The maximum number of words before access to DRAM is relinquished 
* and the initialization of the SCU DMA registers with the communication 
* parameters must be set before this routine is called. This is done 
* by calling the routine scu_init.
*
* If STAND_ALONE = 1 it can be called from C as:
*
* wfm_comm_forward(float *ab0,                 ; ab for mu = 0 
*	            float *ab1,                 ; ab for mu = 1 
*	            float *ab2,                 ; ab for mu = 2 
*   	            float *ab3,                 ; ab for mu = 3 
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

	.def _wfm_comm_forward

*---------------------------------------------------------------------------------------
* References
*---------------------------------------------------------------------------------------
	.ref	ab0
	.ref	ab1
	.ref	ab2
	.ref	ab3
	.ref	_WilsonSCUCommForward
        .ref	_wfm_copy_forward
        .ref	_gjp_local_axis
		
*---------------------------------------------------------------------------------------
* definitions
*---------------------------------------------------------------------------------------
FP		.set	AR3			; use FP for arguments only
WILSON_AD	.set	AR0
SCU		.set	AR1
ATMP		.set	AR2
TMP		.set	R0
LOCAL_AXIS      .set    AR6
AXIS_FLAG       .set    R2

*---------------------------------------------------------------------------------------
* allocations
*--------------------------------------------------------------------------------------
        .text
laxis   .word   _gjp_local_axis
	 
		
****************************************************************************************
* _wfm_comm_forward
****************************************************************************************
	.text
_wfm_comm_forward:

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
	LDP	@ab0
	LDI     *-FP(2), 	TMP			; Address of destination
	STI	TMP,		@ab0
	LDI     *-FP(3), 	TMP			; Address of destination
	STI	TMP,		@ab1
	LDI     *-FP(4), 	TMP			; Address of destination
	STI	TMP,		@ab2
	LDI     *-FP(5), 	TMP			; Address of destination
	STI	TMP,		@ab3
	LDI     *-FP(6), 	WILSON_AD		; Wilson structure
	.endif

*---------------------------------------------------------------------------------------
* Routine starts here
*---------------------------------------------------------------------------------------
	PUSH    DP
	
	PUSH	WILSON_AD
	
	LDP	@ab3

* Receive from the -mu direction
* Receive base-addresses = 0 + ab[mu];
	LDI	@ab3,				TMP
	PUSH	TMP						; -T
	LDI	@ab2,				TMP
	PUSH	TMP						; -Z
	LDI	@ab1,				TMP
	PUSH	TMP						; -Y
	LDI	@ab0,				TMP
	PUSH	TMP						; -X

* Send in  the +mu direction
* Send base-addresses = wilson_p->comm_offset[mu] + ab[mu];
	LDI	@ab3,					TMP
	ADDI	*+WILSON_AD(Wilson.comm_offset+3),	TMP
	PUSH	TMP						; +T
	LDI	@ab2,					TMP
	ADDI	*+WILSON_AD(Wilson.comm_offset+2),	TMP
	PUSH	TMP						; +Z
	LDI	@ab1,					TMP
	ADDI	*+WILSON_AD(Wilson.comm_offset+1),	TMP
	PUSH	TMP						; +Y
	LDI	@ab0,					TMP
	ADDI	*+WILSON_AD(Wilson.comm_offset+0),	TMP
	PUSH	TMP						; +X

	CALL	_WilsonSCUCommForward
	
	SUBI	8,	SP

	POP	WILSON_AD
	
        LDP     @laxis
        LDI     @laxis,                         LOCAL_AXIS;  local axis
		
	LDP	@ab0
* If needed call the wfm_copy_forward routine to do
* any local communications by copying.
        LDI     *+LOCAL_AXIS(5),        AXIS_FLAG
        BZ      NCP
        PUSH    WILSON_AD
        LDI     @ab3,   TMP
        PUSH    TMP
        LDI     @ab2,   TMP
        PUSH    TMP
        LDI     @ab1,   TMP
        PUSH    TMP
        LDI     @ab0,   TMP
        PUSH    TMP
	CALL    _wfm_copy_forward
        SUBI    5,      SP
NCP:    NOP
		
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

