**--------------------------------------------------------------------
**  CVS keywords
**
**  $Author: chulwoo $
**  $Date: 2004-06-04 21:14:10 $
**  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_wilson/qcdsp_nos/wfm_scu_wait.asm,v 1.3 2004-06-04 21:14:10 chulwoo Exp $
**  $Id: wfm_scu_wait.asm,v 1.3 2004-06-04 21:14:10 chulwoo Exp $
**  $Name: not supported by cvs2svn $
**  $Locker:  $
**  $Revision: 1.3 $
**  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_wilson/qcdsp_nos/wfm_scu_wait.asm,v $
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
	.ref	_wfm_wire_map
	.ref _gjp_local_axis

*---------------------------------------------------------------------------------------
* definitions
*--------------------------------------------------------------------------------------
FP		.set	AR3			; use FP for arguments only
SCUB		.set	AR1
SCU		.set	AR2
WIRE_MAP	.set	AR4
LOCAL_AXIS	.set	AR5	
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
* _wfm_scu_wait
****************************************************************************************
	.text
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

	PUSH    DP
*---------------------------------------------------------------------------------------
* Routine starts here
*---------------------------------------------------------------------------------------
	LDP	@wmap
	LDI	@wmap,				WIRE_MAP	; wire map 

	LDP	@laxis
	LDI	@laxis,				LOCAL_AXIS	; local axis 

	LDP	@scu_b
	LDI	@scu_b,				SCUB		; SCU base address

* Poll each direction.
*-------------------------------------------------------------------------
	LDI	*+LOCAL_AXIS(0),		AXIS_FLAG
	BNZ	CHKX
	LDI	*+WIRE_MAP(0),			WIRE
	ADDI	WIRE,	SCUB,			SCU
PLXP:	LDI	*+SCU(scu.poll_wire0),		TMP ; +X
	CMPI	0, TMP
	BNZ	PLXP
	LDI	*+WIRE_MAP(1),			WIRE
	ADDI	WIRE,	SCUB,			SCU
PLXM:	LDI	*+SCU(scu.poll_wire0),		TMP ; -X
	CMPI	0, TMP
	BNZ	PLXM	
CHKX:	NOP

	LDI	*+LOCAL_AXIS(1),		AXIS_FLAG
	BNZ	CHKY
	LDI	*+WIRE_MAP(2),			WIRE
	ADDI	WIRE,	SCUB,			SCU
PLYP:	LDI	*+SCU(scu.poll_wire0),		TMP ; +Y
	CMPI	0, TMP
	BNZ	PLYP
	LDI	*+WIRE_MAP(3),			WIRE
	ADDI	WIRE,	SCUB,			SCU
PLYM:	LDI	*+SCU(scu.poll_wire0),		TMP ; -Y
	CMPI	0, TMP
	BNZ	PLYM	
CHKY:	NOP

	LDI	*+LOCAL_AXIS(2),		AXIS_FLAG
	BNZ	CHKZ
	LDI	*+WIRE_MAP(4),			WIRE
	ADDI	WIRE,	SCUB,			SCU
PLZP:	LDI	*+SCU(scu.poll_wire0),		TMP ; +Z
	CMPI	0, TMP
	BNZ	PLZP
	LDI	*+WIRE_MAP(5),			WIRE
	ADDI	WIRE,	SCUB,			SCU
PLZM:	LDI	*+SCU(scu.poll_wire0),		TMP ; -Z
	CMPI	0, TMP
	BNZ	PLZM	
CHKZ:	NOP

	LDI	*+LOCAL_AXIS(3),		AXIS_FLAG
	BNZ	CHKT
	LDI	*+WIRE_MAP(6),			WIRE
	ADDI	WIRE,	SCUB,			SCU
PLTP:	LDI	*+SCU(scu.poll_wire0),		TMP ; +T
	CMPI	0, TMP
	BNZ	PLTP
	LDI	*+WIRE_MAP(7),			WIRE
	ADDI	WIRE,	SCUB,			SCU
PLTM:	LDI	*+SCU(scu.poll_wire0),		TMP ; -T
	CMPI	0, TMP
	BNZ	PLTM	
CHKT:	NOP

	
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

