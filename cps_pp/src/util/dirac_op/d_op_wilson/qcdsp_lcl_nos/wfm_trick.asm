**--------------------------------------------------------------------
**  CVS keywords
**
**  $Author: chulwoo $
**  $Date: 2004-06-04 21:14:09 $
**  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_wilson/qcdsp_lcl_nos/wfm_trick.asm,v 1.3 2004-06-04 21:14:09 chulwoo Exp $
**  $Id: wfm_trick.asm,v 1.3 2004-06-04 21:14:09 chulwoo Exp $
**  $Name: not supported by cvs2svn $
**  $Locker:  $
**  $Revision: 1.3 $
**  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_wilson/qcdsp_lcl_nos/wfm_trick.asm,v $
**  $State: Exp $
**
**--------------------------------------------------------------------
*****************************************************************************************
*---------------------------------------------------------------------------------------
*
* _wfm_trick:
*
* If STAND_ALONE = 1 it can be called from C as:
*
*  	   wfm_trick(chi, 
*		     float *ab0, float *ab1, float *ab2, float *ab3,
*		     float sign, Wilson *wilson_p, int cb);
*
* where:
*
* chi      --> 4 component spinors (not padded)
* ab       --> 2 component spproj'ed padded spinors for the backward direction
* sign     --> sign in spin project
* wilson_p --> pointer to Wilson struct
* cb       --> checkerboard
*
* If STAND_ALONE = 0 it can not be called from C. Instead the following is expected
* to be set up before the routine is called:
*
* 1) all the .ref in the STAND_ALONE = 0  case below must be defined
* 2) The following registers must be loaded as follows:
*    AR0  = wilson_p
*    AR4  = chi
*    AR7  = cb
*    R6  = sign
*
*---------------------------------------------------------------------------------------
*****************************************************************************************

	.version	30

	.include	"../../include/wilson.hasm"
	.include	"wfm_macros.hasm"

	.def	 _wfm_trick

	.ref	_wfm_trick_segment_2
	.ref	_wfm_trick_segment_3

*---------------------------------------------------------------------------------------
* References
*---------------------------------------------------------------------------------------
	.ref	direct
	.ref	bank0
	.ref	bank1a
	.ref	bank1b
	.ref	bank2a
	.ref	bank2b
	.ref	bank3a
	.ref	bank3b
	.ref	bank4a
	.ref	bank4b
	.ref	cb_cntrl_b	
	.ref	ab0
	.ref	ab1
	.ref	ab2
	.ref	ab3
	.ref	ab0_t
	.ref	ab1_t
	.ref	ab2_t
	.ref	ab3_t

*---------------------------------------------------------------------------------------
* definitions
*---------------------------------------------------------------------------------------
FP	.set	AR3			; use FP for arguments only

PTRB	.set	AR7
LCYZT	.set	R7
LCC	.set	AR6
SIGN	.set	R6
CBB	.set	AR0
WILSON_AD .set	AR0

AB0	.set	AR0
AB1	.set	AR1
AB2	.set	AR2
AB3	.set	AR3

CHI	.set	AR4
CHI_ST	.set	AR5

CHI_R0	.set	R0
CHI_I0	.set	R1
CHI_R1	.set	R2
CHI_I1	.set	R3
CHI_R2	.set	R4
CHI_I2	.set	R5

TMP	.set	R0

*---------------------------------------------------------------------------------------
* Local data and pointers
*---------------------------------------------------------------------------------------
*  Set up circular buffer modes.
 	.if	CBUF=1
	.sect	"T:wfm0"
	DEF_CB_MODE	cb_reg0,  1,  1,  1,  1,  1, 0, 0, 1, 1, 0
	DEF_CB_MODE	cb_reg1,  6,  1,  6,  6,  6, 1, 0, 1, 1, 1
	.endif

****************************************************************************************
* _wfm_trick
****************************************************************************************
	.sect 	"T:wfm0"
_wfm_trick:

*---------------------------------------------------------------------------------------
* If STAND_ALONE = 1: C-callable
*---------------------------------------------------------------------------------------
	.if 	STAND_ALONE = 1
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
	LDI     *-FP(2), 	CHI  				; Address of destination
	LDI     *-FP(3), 	TMP 				; Address of AB0
	STI	TMP,		@ab0
	LDI     *-FP(4), 	TMP 				; Address of AB1
	STI	TMP,		@ab1
	LDI     *-FP(5),	TMP 				; Address of AB2
	STI	TMP,		@ab2
	LDI     *-FP(6), 	TMP 				; Address of AB3
	STI	TMP,		@ab3
	LDF     *-FP(7),	SIGN				; Sign in trick
	LDI     *-FP(8), 	WILSON_AD			; Wilson structure
	LDI     *-FP(9), 	PTRB				; Checkerboard of AB 
	.endif

*---------------------------------------------------------------------------------------
* Routine starts here
*---------------------------------------------------------------------------------------
	PUSH    DP
	LDP	@ab0

* Various addresses used for padding
	LSH	2,				PTRB		; shift left by 2 and construct
	OR	3,				PTRB		; the pointer into the PTRB array
	MPYI	*+WILSON_AD(Wilson.offset),	PTRB		; calculate and load the
	ADDI	*+WILSON_AD(Wilson.ptr),	PTRB		; base address of the PTRB array
	LDI	*+WILSON_AD(Wilson.yztmax),	LCYZT		; load the Y,Z,T loop counter

* Initialize address to store CHI in DRAM 
	LDI	CHI,				CHI_ST		

	.if	CBUF=1
* Prepare addresses for Circular Buffer access.
* CHI, AB0, AB1, AB2, AB3 accesses are in mode: 6, 1, 6, 6, 6, 1, 0, 1, 1, 1
* All other accesses are in mode	      : 1, 1, 1, 1, 1, 0, 0, 1, 1, 0
* The sub-banks must alternate to start the next process. 
* Stores must be through bank0 = 0h. CHI_ST is in bank0.
* Because clrbit=1 after a store a new process can be started 
* even in the same bank as before the store.
*  Set up circular buffer.
	LDI     @cb_cntrl_b, 	CBB			;  CBB has cbctrl address
	LDI     @cb_reg0, 	TMP
	STI     TMP, 		*+CBB(0)		;  Set flags for cntrl reg 0
	LDI     @cb_reg1,	TMP
	STI     TMP,		*+CBB(1)		;  Set flags for cntrl reg 1
* Put variables in the appropriate banks
	LDI	@ab0,		TMP
	ADDI	@bank1a,	TMP
	STI	TMP,		@ab0_t
	LDI	@ab1,		TMP
	ADDI	@bank1b,	TMP
	STI	TMP,		@ab1_t
	LDI	@ab2,		TMP
	ADDI	@bank1a,	TMP
	STI	TMP,		@ab2_t
	LDI	@ab3,		TMP
	ADDI	@bank1b,	TMP
	STI	TMP,		@ab3_t
	ADDI	@bank1a,	CHI
	.else
	LDI	@ab0,		TMP
	STI	TMP,		@ab0_t
	LDI	@ab1,		TMP
	STI	TMP,		@ab1_t
	LDI	@ab2,		TMP
	STI	TMP,		@ab2_t
	LDI	@ab3,		TMP
	STI	TMP,		@ab3_t
	.endif


*
* Start main loop
*
LYZT: 								; loop over Y, Z, T

	LDI	@ab0_t,		TMP
	ADDI	*PTRB++(1),	TMP,		AB0		; add displacement to AB0
	LDI	@ab1_t,		TMP
	ADDI	*PTRB++(1),	TMP,		AB1		; add displacement to AB1
	LDI	@ab2_t,		TMP
	ADDI	*PTRB++(1),	TMP,		AB2		; add displacement to AB2
	LDI	@ab3_t,		TMP
	ADDI	*PTRB++(1),	TMP,		AB3		; add displacement to AB3
	LDI	*PTRB++(1),	RC				; set the X loop counter

*
* Start X loop
*
	RPTB	LX						; loop over X

******************************************************************************
*			trick
******************************************************************************
*
* Trick AB into CHI
*

* First and second spinor components

	LDI	1,		LCC				; initialize the loop counter
LCCM:
	LDF	*AB0++(1),	CHI_R0				; + AB0
	LDF	*AB0++(1),	CHI_I0				; + AB0
	LDF	*AB0++(1),	CHI_R1				; + AB0
	LDF	*AB0++(1),	CHI_I1				; + AB0
	LDF	*AB0++(1),	CHI_R2				; + AB0
	LDF	*AB0++(1),	CHI_I2				; + AB0

	ADDF	*AB1++(1),	CHI_R0				; + AB1
	ADDF	*AB1++(1),	CHI_I0				; + AB1
	ADDF	*AB1++(1),	CHI_R1				; + AB1
	ADDF	*AB1++(1),	CHI_I1				; + AB1
	ADDF	*AB1++(1),	CHI_R2				; + AB1
	ADDF	*AB1++(1),	CHI_I2				; + AB1

	ADDF	*AB2++(1),	CHI_R0				; + AB2
	ADDF	*AB2++(1),	CHI_I0				; + AB2
	ADDF	*AB2++(1),	CHI_R1				; + AB2
	ADDF	*AB2++(1),	CHI_I1				; + AB2
	ADDF	*AB2++(1),	CHI_R2				; + AB2
	ADDF	*AB2++(1),	CHI_I2				; + AB2

	ADDF	*AB3++(1),	CHI_R0				; + AB3
	ADDF	*AB3++(1),	CHI_I0				; + AB3
	ADDF	*AB3++(1),	CHI_R1				; + AB3
	ADDF	*AB3++(1),	CHI_I1				; + AB3
	ADDF	*AB3++(1),	CHI_R2				; + AB3
	ADDF	*AB3++(1),	CHI_I2				; + AB3

	ADDF	*CHI++(1),	CHI_R0				; + CHI
	ADDF	*CHI++(1),	CHI_I0				; + CHI
	ADDF	*CHI++(1),	CHI_R1				; + CHI
	ADDF	*CHI++(1),	CHI_I1				; + CHI
	ADDF	*CHI++(1),	CHI_R2				; + CHI
	ADDF	*CHI++(1),	CHI_I2				; + CHI

	.if	ROUND=1
	RND	CHI_R0
	.endif
	STF	CHI_R0,		*CHI_ST++(1)			; -> CHI_R0  DRAM
	.if	ROUND=1
	RND	CHI_I0
	.endif
	STF	CHI_I0,		*CHI_ST++(1)			; -> CHI_I0  DRAM
	.if	ROUND=1
	RND	CHI_R1
	.endif
	STF	CHI_R1,		*CHI_ST++(1)			; -> CHI_R1  DRAM
	.if	ROUND=0
	DBD	LCC,		LCCM				; branch to LCCM (if round=0)
	.endif
	.if	ROUND=1
	RND	CHI_I1
	.endif
	STF	CHI_I1,		*CHI_ST++(1)			; -> CHI_I1  DRAM
	.if	ROUND=1
	RND	CHI_R2
	DBD	LCC,		LCCM				; branch to LCCM (if round=1)
	.endif
	STF	CHI_R2,		*CHI_ST++(1)			; -> CHI_R2  DRAM
	.if	ROUND=1
	RND	CHI_I2
	.endif
	STF	CHI_I2,		*CHI_ST++(1)			; -> CHI_I2  DRAM

*---------------------------------------------------------------; end of LCCM loop

	SUBI	6,		AB0
	SUBI	6,		AB1
	SUBI	12,		AB2
	SUBI	12,		AB3

* Third spinor component
	NEGF	*AB0++(1),	CHI_I0				; - AB0
	LDF	*AB0++(1),	CHI_R0				; + AB0
	NEGF	*AB0++(1),	CHI_I1				; - AB0
	LDF	*AB0++(1),	CHI_R1				; + AB0
	NEGF	*AB0++(1),	CHI_I2				; - AB0
	LDF	*AB0--(11),	CHI_R2				; + AB0

	CALL	_wfm_trick_segment_2			; trick segment 2

	MPYF	SIGN,		CHI_R0				; * SIGN
	MPYF	SIGN,		CHI_I0				; * SIGN
	MPYF	SIGN,		CHI_R1				; * SIGN
	MPYF	SIGN,		CHI_I1				; * SIGN
	MPYF	SIGN,		CHI_R2				; * SIGN
	MPYF	SIGN,		CHI_I2				; * SIGN

	ADDF	*CHI++(1),	CHI_R0				; + CHI
	ADDF	*CHI++(1),	CHI_I0				; + CHI
	ADDF	*CHI++(1),	CHI_R1				; + CHI
	ADDF	*CHI++(1),	CHI_I1				; + CHI
	ADDF	*CHI++(1),	CHI_R2				; + CHI
	ADDF	*CHI++(1),	CHI_I2				; + CHI

	.if	ROUND=1
	RND	CHI_R0
	.endif
	STF	CHI_R0,		*CHI_ST++(1)			; -> CHI_R0  DRAM
	.if	ROUND=1
	RND	CHI_I0
	.endif
	STF	CHI_I0,		*CHI_ST++(1)			; -> CHI_I0  DRAM
	.if	ROUND=1
	RND	CHI_R1
	.endif
	STF	CHI_R1,		*CHI_ST++(1)			; -> CHI_R1  DRAM
	.if	ROUND=1
	RND	CHI_I1
	.endif
	STF	CHI_I1,		*CHI_ST++(1)			; -> CHI_I1  DRAM
	.if	ROUND=1
	RND	CHI_R2
	.endif
	STF	CHI_R2,		*CHI_ST++(1)			; -> CHI_R2  DRAM
	.if	ROUND=1
	RND	CHI_I2
	.endif
	STF	CHI_I2,		*CHI_ST++(1)			; -> CHI_I2  DRAM

* Fourth spinor component
	NEGF	*AB0++(1),	CHI_I0				; - AB0
	LDF	*AB0++(1),	CHI_R0				; + AB0
	NEGF	*AB0++(1),	CHI_I1				; - AB0
	LDF	*AB0++(1),	CHI_R1				; + AB0
	NEGF	*AB0++(1),	CHI_I2				; - AB0
	LDF	*AB0++(7),	CHI_R2				; + AB0

	CALL	_wfm_trick_segment_3			; trick segment 3

	MPYF	SIGN,		CHI_R0				; * SIGN
	MPYF	SIGN,		CHI_I0				; * SIGN
	MPYF	SIGN,		CHI_R1				; * SIGN
	MPYF	SIGN,		CHI_I1				; * SIGN
	MPYF	SIGN,		CHI_R2				; * SIGN
	MPYF	SIGN,		CHI_I2				; * SIGN

	ADDF	*CHI++(1),	CHI_R0				; + CHI
	ADDF	*CHI++(1),	CHI_I0				; + CHI
	ADDF	*CHI++(1),	CHI_R1				; + CHI
	ADDF	*CHI++(1),	CHI_I1				; + CHI
	ADDF	*CHI++(1),	CHI_R2				; + CHI
	ADDF	*CHI++(1),	CHI_I2				; + CHI

	.if	ROUND=1
	RND	CHI_R0
	.endif
	STF	CHI_R0,		*CHI_ST++(1)			; -> CHI_R0  DRAM
	.if	ROUND=1
	RND	CHI_I0
	.endif
	STF	CHI_I0,		*CHI_ST++(1)			; -> CHI_I0  DRAM
	.if	ROUND=1
	RND	CHI_R1
	.endif
	STF	CHI_R1,		*CHI_ST++(1)			; -> CHI_R1  DRAM
	.if	ROUND=1
	RND	CHI_I1
	.endif
	STF	CHI_I1,		*CHI_ST++(1)			; -> CHI_I1  DRAM
	.if	ROUND=1
	RND	CHI_R2
	.endif
	STF	CHI_R2,		*CHI_ST++(1)			; -> CHI_R2  DRAM
	.if	ROUND=1
	RND	CHI_I2
	.endif
LX:	STF	CHI_I2,		*CHI_ST++(1)			; -> CHI_I2  DRAM
*---------------------------------------------------------------; end of X loop

	SUBI	1,		LCYZT
	BNN	LYZT
*---------------------------------------------------------------; end of Y, Z, T loop

	POP	DP

*-------------------------------------------------------------------------
*
*  Restore the registers before returning to the C program                
*
	.if 	STAND_ALONE = 1
	POP     DP
	POP     FP              				;  This is our frame pointer
	POP     AR7
	POP     AR6
	POP     AR5
	POP     AR4
	POPF    R7
	POPF    R6
	POP     R5
	POP     R4

	POP     FP              				;  This is our caller's frame pointer
	.endif

	RETS
	.end
