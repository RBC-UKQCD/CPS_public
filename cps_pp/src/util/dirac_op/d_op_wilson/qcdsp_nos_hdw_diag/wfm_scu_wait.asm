**--------------------------------------------------------------------
**  CVS keywords
**
**  $Author: zs $
**  $Date: 2004-08-18 11:58:00 $
**  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_wilson/qcdsp_nos_hdw_diag/wfm_scu_wait.asm,v 1.4 2004-08-18 11:58:00 zs Exp $
**  $Id: wfm_scu_wait.asm,v 1.4 2004-08-18 11:58:00 zs Exp $
**  $Name: not supported by cvs2svn $
**  $Locker:  $
**  $Revision: 1.4 $
**  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_wilson/qcdsp_nos_hdw_diag/wfm_scu_wait.asm,v $
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
	.ref	_gjp_local_axis
	.ref	_wfm_scu_diag
	.ref	_wfm_max_scu_poll
	.ref	_wilson_scu_error

*---------------------------------------------------------------------------------------
* definitions
*--------------------------------------------------------------------------------------
FP		.set	AR3			; use FP for arguments only
SCUB		.set	AR1
SCU		.set	AR2
WIRE_MAP	.set	AR4
LOCAL_AXIS	.set	AR5	
SCU_DIAG	.set	AR6
SCU_DIAG_BASE	.set	AR7
TMP		.set	R0
WIRE		.set	R1
AXIS_FLAG	.set	R2
POLL_COUNT	.set	R3
MAX_SCU_POLL	.set	R4
FAILED_FLAG	.set	R5
	
*---------------------------------------------------------------------------------------
* allocations
*--------------------------------------------------------------------------------------
	.text
wmap		.word	_wfm_wire_map
wscu_diag	.word	_wfm_scu_diag
laxis		.word	_gjp_local_axis


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

	LDP	@wscu_diag
	LDI	@wscu_diag,			SCU_DIAG_BASE	; scu diagnostic array 

	LDP	@_wfm_max_scu_poll
	LDI	@_wfm_max_scu_poll,		MAX_SCU_POLL    ; maximum scu poll

	LDP	@scu_b
	LDI	@scu_b,				SCUB		; SCU base address



* Set failed flag	
*-------------------------------------------------------------------------
	LDI	0,				FAILED_FLAG


	
* Poll each direction.
*-------------------------------------------------------------------------

	LDI	*+LOCAL_AXIS(0),		AXIS_FLAG
	BNZ	CHKX
*---------------------------------------------------------------------- +X
	LDI	*+WIRE_MAP(0),			WIRE
	ADDI	WIRE,	SCUB,			SCU
	ADDI	WIRE,	SCU_DIAG_BASE,		SCU_DIAG
	STI	MAX_SCU_POLL,			*+SCU_DIAG(3)	
PLXP:	LDI	*+SCU_DIAG(3),			POLL_COUNT
	BZ	FXP
	SUBI	1,				POLL_COUNT
	STI	POLL_COUNT,			*+SCU_DIAG(3)	
	LDI	*+SCU(scu.poll_wire0),		TMP ; +X
	CMPI	0, TMP	
	BNZ	PLXP
	BU	CNTXP
FXP:	LDI	1,				FAILED_FLAG
CNTXP:	NOP
*---------------------------------------------------------------------- -X
	LDI	*+WIRE_MAP(1),			WIRE
	ADDI	WIRE,	SCUB,			SCU
	ADDI	WIRE,	SCU_DIAG_BASE,		SCU_DIAG
	STI	MAX_SCU_POLL,			*+SCU_DIAG(3)	
PLXM:	LDI	*+SCU_DIAG(3),			POLL_COUNT
	BZ	FXM
	SUBI	1,				POLL_COUNT
	STI	POLL_COUNT,			*+SCU_DIAG(3)	
	LDI	*+SCU(scu.poll_wire0),		TMP ; -X
	CMPI	0, TMP
	BNZ	PLXM	
	BU	CNTXM
FXM:	LDI	1,				FAILED_FLAG
CNTXM:	NOP

CHKX:	NOP


	
	
	LDI	*+LOCAL_AXIS(1),		AXIS_FLAG
	BNZ	CHKY
*---------------------------------------------------------------------- +Y
	LDI	*+WIRE_MAP(2),			WIRE
	ADDI	WIRE,	SCUB,			SCU
	ADDI	WIRE,	SCU_DIAG_BASE,		SCU_DIAG
	STI	MAX_SCU_POLL,			*+SCU_DIAG(3)	
PLYP:	LDI	*+SCU_DIAG,			POLL_COUNT
	BZ	FYP
	SUBI	1,				POLL_COUNT
	STI	POLL_COUNT,			*+SCU_DIAG(3)	
	LDI	*+SCU(scu.poll_wire0),		TMP ; +Y
	CMPI	0, TMP
	BNZ	PLYP
	BU	CNTYP
FYP:	LDI	1,				FAILED_FLAG
CNTYP:	NOP
*---------------------------------------------------------------------- -Y
	LDI	*+WIRE_MAP(3),			WIRE
	ADDI	WIRE,	SCUB,			SCU
	ADDI	WIRE,	SCU_DIAG_BASE,		SCU_DIAG
	STI	MAX_SCU_POLL,			*+SCU_DIAG(3)	
PLYM:	LDI	*+SCU_DIAG(3),			POLL_COUNT
	BZ	FYM
	SUBI	1,				POLL_COUNT
	STI	POLL_COUNT,			*+SCU_DIAG(3)	
	LDI	*+SCU(scu.poll_wire0),		TMP ; -Y
	CMPI	0, TMP
	BNZ	PLYM	
	BU	CNTYM
FYM:	LDI	1,				FAILED_FLAG
CNTYM:	NOP

CHKY:	NOP


	
	
	LDI	*+LOCAL_AXIS(2),		AXIS_FLAG
	BNZ	CHKZ
*---------------------------------------------------------------------- +Z
	LDI	*+WIRE_MAP(4),			WIRE
	ADDI	WIRE,	SCUB,			SCU
	ADDI	WIRE,	SCU_DIAG_BASE,		SCU_DIAG
	STI	MAX_SCU_POLL,			*+SCU_DIAG(3)	
PLZP:	LDI	*+SCU_DIAG(3),			POLL_COUNT
	BZ	FZP
	SUBI	1,				POLL_COUNT
	STI	POLL_COUNT,			*+SCU_DIAG(3)	
	LDI	*+SCU(scu.poll_wire0),		TMP ; +Z
	CMPI	0, TMP
	BNZ	PLZP
	BU	CNTZP
FZP:	LDI	1,				FAILED_FLAG
CNTZP:	NOP

*---------------------------------------------------------------------- -Z
	LDI	*+WIRE_MAP(5),			WIRE
	ADDI	WIRE,	SCUB,			SCU
	ADDI	WIRE,	SCU_DIAG_BASE,		SCU_DIAG
	STI	MAX_SCU_POLL,			*+SCU_DIAG(3)	
PLZM:	LDI	*+SCU_DIAG(3),			POLL_COUNT
	BZ	FZM
	SUBI	1,				POLL_COUNT
	STI	POLL_COUNT,			*+SCU_DIAG(3)	
	LDI	*+SCU(scu.poll_wire0),		TMP ; -Z
	CMPI	0, TMP
	BNZ	PLZM	
	BU	CNTZM
FZM:	LDI	1,				FAILED_FLAG
CNTZM:	NOP

CHKZ:	NOP


	
	
	
	LDI	*+LOCAL_AXIS(3),		AXIS_FLAG
	BNZ	CHKT
*---------------------------------------------------------------------- +T
	LDI	*+WIRE_MAP(6),			WIRE
	ADDI	WIRE,	SCUB,			SCU
	ADDI	WIRE,	SCU_DIAG_BASE,		SCU_DIAG
	STI	MAX_SCU_POLL,			*+SCU_DIAG(3)	
PLTP:	LDI	*+SCU_DIAG(3),			POLL_COUNT
	BZ	FTP
	SUBI	1,				POLL_COUNT
	STI	POLL_COUNT,			*+SCU_DIAG(3)	
	LDI	*+SCU(scu.poll_wire0),		TMP ; +T
	CMPI	0, TMP
	BNZ	PLTP
	BU	CNTTP
FTP:	LDI	1,				FAILED_FLAG
CNTTP:	NOP

*---------------------------------------------------------------------- -T
	LDI	*+WIRE_MAP(7),			WIRE
	ADDI	WIRE,	SCUB,			SCU
	ADDI	WIRE,	SCU_DIAG_BASE,		SCU_DIAG
	STI	MAX_SCU_POLL,			*+SCU_DIAG(3)	
PLTM:	LDI	*+SCU_DIAG(3),			POLL_COUNT
	BZ	FTM
	SUBI	1,				POLL_COUNT
	STI	POLL_COUNT,			*+SCU_DIAG(3)	
	LDI	*+SCU(scu.poll_wire0),		TMP ; -T
	CMPI	0, TMP
	BNZ	PLTM	
	BU	CNTTM
FTM:	LDI	1,				FAILED_FLAG
CNTTM:	NOP

CHKT:	NOP
	

* Copy the scu status register 813040 to SCU_DIAG(2)	
*-------------------------------------------------------------------------
	LDI	*+SCUB(scu.poll_all),		TMP
	STI	TMP,				*+SCU_DIAG_BASE(2)
	
* Copy the scu error status register 813060 - 813067 to 
* SCU_DIAG_BASE(11) - 18
*-------------------------------------------------------------------------
	LDI	*+SCUB(scu.error0),		TMP
	STI	TMP,				*+SCU_DIAG_BASE(11)
	LDI	*+SCUB(scu.error1),		TMP
	STI	TMP,				*+SCU_DIAG_BASE(12)
	LDI	*+SCUB(scu.error2),		TMP
	STI	TMP,				*+SCU_DIAG_BASE(13)
	LDI	*+SCUB(scu.error3),		TMP
	STI	TMP,				*+SCU_DIAG_BASE(14)
	LDI	*+SCUB(scu.error4),		TMP
	STI	TMP,				*+SCU_DIAG_BASE(15)
	LDI	*+SCUB(scu.error5),		TMP
	STI	TMP,				*+SCU_DIAG_BASE(16)
	LDI	*+SCUB(scu.error6),		TMP
	STI	TMP,				*+SCU_DIAG_BASE(17)
	LDI	*+SCUB(scu.error7),		TMP
	STI	TMP,				*+SCU_DIAG_BASE(18)
	

* Branch to scu_wait_error if any communication failed
*-------------------------------------------------------------------------
	CMPI	0,				FAILED_FLAG
	BNZ	scu_wait_error
	

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

*-------------------------------------------------------------------------
* Exit with error
*-------------------------------------------------------------------------
scu_wait_error:	
	CALL	_wilson_scu_error
	RETS
	
