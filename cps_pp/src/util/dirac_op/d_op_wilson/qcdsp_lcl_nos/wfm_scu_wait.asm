**--------------------------------------------------------------------
**  CVS keywords
**
**  $Author: chulwoo $
**  $Date: 2004-06-04 21:14:09 $
**  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_wilson/qcdsp_lcl_nos/wfm_scu_wait.asm,v 1.3 2004-06-04 21:14:09 chulwoo Exp $
**  $Id: wfm_scu_wait.asm,v 1.3 2004-06-04 21:14:09 chulwoo Exp $
**  $Name: not supported by cvs2svn $
**  $Locker:  $
**  $Revision: 1.3 $
**  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_wilson/qcdsp_lcl_nos/wfm_scu_wait.asm,v $
**  $State: Exp $
**
**--------------------------------------------------------------------
****************************************************************************************
*---------------------------------------------------------------------------------------
*
* _wfm_scu_wait
*
* This routine waits until the Serial Communications Unit has completed 
* all communications
*
* If STAND_ALONE = 1 it can be called from C as:
*
* wfm_scu_wait();
*
* If STAND_ALONE = 0 it can not be called from C. Instead the following is expected
* to be set up before the routine is called:
*
* 1) all the .ref in the STAND_ALONE = 0  case below must be defined
*
*---------------------------------------------------------------------------------------
****************************************************************************************

	.version	30

	.include	"wfm_nga_struct.hasm"

	.def _wfm_scu_wait

*---------------------------------------------------------------------------------------
* References
*---------------------------------------------------------------------------------------
	.ref	scu_b

*---------------------------------------------------------------------------------------
* definitions
*--------------------------------------------------------------------------------------
FP	.set	AR3			; use FP for arguments only
SCU	.set	AR1
TMP	.set	R0

****************************************************************************************
* _wfm_scu_wait
****************************************************************************************
	.sect 	"T:wfm0"
_wfm_scu_wait:

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
	LDP	@scu_b
	.endif

*---------------------------------------------------------------------------------------
* Routine starts here
*---------------------------------------------------------------------------------------
	PUSH    DP
	LDP	@scu_b

	LDI	@scu_b,				SCU	; SCU base address
poll?	LDI	*+SCU(scu.poll_all), 		TMP
	CMPI	0, TMP
	BNZ	poll?

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

