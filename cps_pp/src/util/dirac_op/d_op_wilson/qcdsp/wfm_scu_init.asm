**--------------------------------------------------------------------
**  CVS keywords
**
**  $Author: zs $
**  $Date: 2004-08-18 11:57:56 $
**  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_wilson/qcdsp/wfm_scu_init.asm,v 1.4 2004-08-18 11:57:56 zs Exp $
**  $Id: wfm_scu_init.asm,v 1.4 2004-08-18 11:57:56 zs Exp $
**  $Name: not supported by cvs2svn $
**  $Locker:  $
**  $Revision: 1.4 $
**  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_wilson/qcdsp/wfm_scu_init.asm,v $
**  $State: Exp $
**
**--------------------------------------------------------------------
****************************************************************************************
*---------------------------------------------------------------------------------------
*
* _wfm_scu_init
*
* This routine initializes the stride block-length and number-of-blocks in
* the Serial Communications Unit. Since these values depend on the lattice
* size only, and since the SCU does not change them during operation, they
* only need to be set once.
*
* If STAND_ALONE = 1 it can be called from C as:
*
* wfm_scu_init(Wilson *wilson_p);          ; Wilson struct.
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
	.include	"wfm_nga_struct.hasm"

	.def _wfm_scu_init

*---------------------------------------------------------------------------------------
* References
*---------------------------------------------------------------------------------------
	.ref	_WilsonSCUSetDMA

*---------------------------------------------------------------------------------------
* definitions
*--------------------------------------------------------------------------------------
FP		.set	AR3			; use FP for arguments only
WILSON_AD	.set	AR0
TMP		.set	R0

****************************************************************************************
* _wfm_scu_init
****************************************************************************************
	.sect	".text"
_wfm_scu_init:

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
	LDI     *-FP(2), 	WILSON_AD		; Wilson structure
	.endif

*---------------------------------------------------------------------------------------
* Routine starts here
*---------------------------------------------------------------------------------------
	PUSH    DP


* Set up communication parameters by initializing the SCU DMA registers. 
* wilson_p->comm[mu] contains a 32 bit word created by packaging
* the stride, number of blocks, and block length: 
* numblk[10bits]-blklen[10bits]-stride[12bits]
	LDI	*+WILSON_AD(Wilson.comm+3),	TMP
	PUSH	TMP
	LDI	*+WILSON_AD(Wilson.comm+2),	TMP
	PUSH	TMP
	LDI	*+WILSON_AD(Wilson.comm+1),	TMP
	PUSH	TMP
	LDI	*+WILSON_AD(Wilson.comm+0),	TMP
	PUSH	TMP

	CALL	_WilsonSCUSetDMA

	SUBI	4,		SP


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
