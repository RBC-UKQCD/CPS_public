**--------------------------------------------------------------------
**  CVS keywords
**
**  $Author: mcneile $
**  $Date: 2003-06-22 13:34:46 $
**  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_wilson/qcdsp_nos_hdw_diag/wfm_trick_kappa.asm,v 1.1.1.1 2003-06-22 13:34:46 mcneile Exp $
**  $Id: wfm_trick_kappa.asm,v 1.1.1.1 2003-06-22 13:34:46 mcneile Exp $
**  $Name: not supported by cvs2svn $
**  $Locker:  $
**  $Log: not supported by cvs2svn $
**  Revision 1.2  2001/06/19 18:13:11  anj
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
**  Revision 1.2  2001/05/25 06:16:08  cvs
**  Added CVS keywords to phys_v4_0_0_preCVS
**
**  $RCSfile: wfm_trick_kappa.asm,v $
**  $Revision: 1.1.1.1 $
**  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_wilson/qcdsp_nos_hdw_diag/wfm_trick_kappa.asm,v $
**  $State: Exp $
**
**--------------------------------------------------------------------
*****************************************************************************************
*---------------------------------------------------------------------------------------
*
* _wfm_trick_kappa:
*
* If STAND_ALONE = 1 it can be called from C as:
*
*  	   wfm_trick_kappa(float *chi, float *phi, 
*			   float *ab0, float *ab1, float *ab2, float *ab3,
*		 	   float -kappa_sq, float 1/kappa_sq, 
*			   Wilson *wilson_p, int cb);
*
* where:
*
* chi        --> 4 component spinor result (not padded)
* phi        --> 4 component spinor intermediate (input, not padded)
* ab         --> 2 component spproj'ed padded spinors for the backward direction
* -kappa_sq  --> -kappa^2, kappa is Wilson's parameter
* 1/kappa_sq --> 1/kappa^2, kappa is Wilson's parameter
* sign       --> sign in spin project
* wilson_p   --> pointer to Wilson struct
* cb         --> checkerboard
*
* If STAND_ALONE = 0 it can not be called from C. Instead the following is expected
* to be set up before the routine is called:
*
* 1) all the .ref in the STAND_ALONE = 0  case below must be defined
* 2) The following registers must be loaded as follows:
*    AR0  = wilson_p
*    AR4  = chi
*    AR5  = phi
*    AR7  = cb
*
*---------------------------------------------------------------------------------------
*****************************************************************************************

	.version	30

	.include	"../../include/wilson.hasm"
	.include	"wfm_macros.hasm"

	.def _wfm_trick_kappa

	.ref	_wfm_trick_segment_1
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
	.ref	af0
	.ref	af1
	.ref	af2
	.ref	af3
	.ref	ab0
	.ref	ab1
	.ref	ab2
	.ref	ab3
	.ref	parameters
	.ref	kappa_sq
	.ref	i_kappa_sq
	.ref	m_kappa_sq
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
WILSON_AD .set	AR0
CBB	.set	AR0

AB0	.set	AR0
AB1	.set	AR1
AB2	.set	AR2
AB3	.set	AR3

CHI	.set	AR4
OCHI	.set	AR5
CHI_ST	.set	AR6

LC	.set	AR7	

CHI_R0	.set	R0
CHI_I0	.set	R1
CHI_R1	.set	R2
CHI_I1	.set	R3
CHI_R2	.set	R4
CHI_I2	.set	R5
K2	.set	R6
IK2	.set	R7

TMP	.set	R0

*---------------------------------------------------------------------------------------
* Local data and pointers
*---------------------------------------------------------------------------------------
*  Set up circular buffer modes.
 	.if	CBUF=1
	.sect	"T:wfm0"
*  DEF_CB_MODE	.macro	mode_name, wcntset, gnteeset, lookcnt, highrqset, hysterhigh, 
*				   wcnt_enable, kick_enable, clrbit, waitgntee, zflag
	DEF_CB_MODE	cb_reg0,  1,  1,  1,  1,  1, 0, 0, 1, 1, 0
	DEF_CB_MODE	cb_reg1,  6,  1,  6,  6,  6, 1, 0, 1, 1, 1
	.endif

	.sect 	"T:wfm1"
* reserve memory to temporarily save address of PTRB
ptrb	.space	1
lcyzt	.space	1

****************************************************************************************
* _wfm_trick_kappa
****************************************************************************************
	.sect 	"T:wfm1"
_wfm_trick_kappa

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
	LDI     *-FP(2), 	CHI  			; Address of CHI
	LDI     *-FP(3), 	OCHI  			; Address of old Chi
	LDI     *-FP(4), 	TMP 			; Address of AB0
	STI	TMP,		@ab0
	LDI     *-FP(5), 	TMP 			; Address of AB1
	STI	TMP,		@ab1
	LDI     *-FP(6),	TMP 			; Address of AB2
	STI	TMP,		@ab2
	LDI     *-FP(7), 	TMP 			; Address of AB3
	STI	TMP,		@ab3
	LDF     *-FP(8),	TMP			; -kappa_sq
	NEGF	TMP,		TMP
	RND	TMP
	STF	TMP,		@kappa_sq		; kappa_sq
	LDF     *-FP(9),	TMP
	STF	TMP,		@i_kappa_sq		; 1 / kappa_sq 
	LDI     *-FP(10), 	WILSON_AD		; Wilson structure
	LDI     *-FP(11), 	PTRB			; Checkerboard of AB 
	.endif

*---------------------------------------------------------------------------------------
* Routine starts here
*---------------------------------------------------------------------------------------
	PUSH    DP
	LDP	@ab0

	LDI	CHI,		CHI_ST			; CHI = CHI_ST = address to store CHI
	LDF	@kappa_sq,	K2

* Various addresses used for padding
* For trick:
	LSH	2,				PTRB	; shift left by 2 and construct
	OR	3,				PTRB	; the pointer into the PTRB array
	MPYI	*+WILSON_AD(Wilson.offset),	PTRB	; calculate and load the
	ADDI	*+WILSON_AD(Wilson.ptr),	PTRB	; base address of the PTRB array
	LDI	*+WILSON_AD(Wilson.yztmax),	LCYZT	; load the Y,Z,T loop counter
	STI	PTRB,		@ptrb			; save address of PTRB
*
* Note: DP same for all of cram and is set above
*

	.if	CBUF=1
* Prepare addresses for Circular Buffer access.
* OCHI, CHI, AB0, AB1, AB2, AB3 accesses are in mode: 6, 1, 6, 6, 6, 1, 0, 1, 1, 1
* All other accesses are in mode		    : 1, 1, 1, 1, 1, 0, 0, 1, 1, 0
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
	ADDI	@bank1a,	OCHI
	ADDI	@bank1b,	CHI
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
LYZT: 							; loop over Y, Z, T
	STI	LCYZT,		@lcyzt			; save loop counter
	LDI	@ptrb,		PTRB			; retrieve address of PTRB
	LDI	@ab0_t,		TMP
	ADDI	*PTRB++(1),	TMP,	AB0		; add displacement to AB0
	LDI	@ab1_t,		TMP
	ADDI	*PTRB++(1),	TMP,	AB1		; add displacement to AB1
	LDI	@ab2_t,		TMP
	ADDI	*PTRB++(1),	TMP,	AB2		; add displacement to AB2
	LDI	@ab3_t,		TMP
	ADDI	*PTRB++(1),	TMP,	AB3		; add displacement to AB3
	LDI	*PTRB++(1),	RC			; set the X loop counter
	STI	PTRB,		@ptrb			; save address of PTRB
	LDF	@i_kappa_sq,	IK2			; load 1/ k^2
	
*
* Start X loop
*
	RPTB	LX					; loop over X
******************************************************************************
*			trick
******************************************************************************
*
* Trick AB into CHI
*

* First and second spinor components
	LDI	1,		LC
LCC:
	CALL	_wfm_trick_segment_1			; trick segment 1

	.if	ROUND=1
	RND	CHI_R0
	.endif
	MPYF	K2,		CHI_R0,		CHI_R0		; * kappa^2
	.if	ROUND=1
	RND	CHI_R0
	.endif
	STF	CHI_R0,		*CHI_ST++(1)			; -> CHI_R0
	.if	ROUND=1
	RND	CHI_I0
	.endif
	MPYF	K2,		CHI_I0,		CHI_I0		; * kappa^2
	.if	ROUND=1
	RND	CHI_I0
	.endif
	STF	CHI_I0,		*CHI_ST++(1)			; -> CHI_I0
	.if	ROUND=1
	RND	CHI_R1
	.endif
	MPYF	K2,		CHI_R1,		CHI_R1		; * kappa^2
	.if	ROUND=1
	RND	CHI_R1
	.endif
	STF	CHI_R1,		*CHI_ST++(1)			; -> CHI_R1
	.if	ROUND=1
	RND	CHI_I1
	.endif
	MPYF	K2,		CHI_I1,		CHI_I1		; * kappa^2
	.if	ROUND=1
	RND	CHI_I1
	.endif
	STF	CHI_I1,		*CHI_ST++(1)			; -> CHI_I1
	.if	ROUND=1
	RND	CHI_R2
	.endif
	MPYF	K2,		CHI_R2,		CHI_R2		; * kappa^2
	.if	ROUND=1
	RND	CHI_R2
	.endif


	.if	ROUND=0
	DBD	LC,		LCC				; branch to LCC (if round=0)
	.endif


	STF	CHI_R2,		*CHI_ST++(1)			; -> CHI_R2

	.if	ROUND=1
	RND	CHI_I2
	DBD	LC,		LCC				; branch to LCC (if round=1)
	.endif

	MPYF	K2,		CHI_I2,		CHI_I2		; * kappa^2
	.if	ROUND=1
	RND	CHI_I2
	.endif
	STF	CHI_I2,		*CHI_ST++(1)			; -> CHI_I2

*---------------------------------------------------------------; end of LCC loop

	SUBI	6,	AB0
	SUBI	6,	AB1
	SUBI	12,	AB2
	SUBI	12,	AB3

* Third spinor component
	MPYF	*OCHI++(1),	IK2,	CHI_R0		; OCHI / kappa_sq
	MPYF	*OCHI++(1),	IK2,	CHI_I0		; OCHI / kappa_sq
	MPYF	*OCHI++(1),	IK2,	CHI_R1		; OCHI / kappa_sq
	MPYF	*OCHI++(1),	IK2,	CHI_I1		; OCHI / kappa_sq
	MPYF	*OCHI++(1),	IK2,	CHI_R2		; OCHI / kappa_sq
	MPYF	*OCHI++(1),	IK2,	CHI_I2		; OCHI / kappa_sq
	
	SUBF	*CHI++(1),	CHI_R0				; - CHI
	SUBF	*CHI++(1),	CHI_I0				; - CHI
	SUBF	*CHI++(1),	CHI_R1				; - CHI
	SUBF	*CHI++(1),	CHI_I1				; - CHI
	SUBF	*CHI++(1),	CHI_R2				; - CHI
	SUBF	*CHI++(1),	CHI_I2				; - CHI

	SUBF	*AB0++(1),	CHI_I0				; - AB0
	ADDF	*AB0++(1),	CHI_R0				; + AB0
	SUBF	*AB0++(1),	CHI_I1				; - AB0
	ADDF	*AB0++(1),	CHI_R1				; + AB0
	SUBF	*AB0++(1),	CHI_I2				; - AB0
	ADDF	*AB0--(11),	CHI_R2				; + AB0

	CALL	_wfm_trick_segment_2			; trick segment 2

	.if	ROUND=1
	RND	CHI_R0
	.endif
	MPYF	K2,		CHI_R0,		CHI_R0		; * kappa^2
	.if	ROUND=1
	RND	CHI_R0
	.endif
	STF	CHI_R0,		*CHI_ST++(1)			; -> CHI_R0
	.if	ROUND=1
	RND	CHI_I0
	.endif
	MPYF	K2,		CHI_I0,		CHI_I0		; * kappa^2
	.if	ROUND=1
	RND	CHI_I0
	.endif
	STF	CHI_I0,		*CHI_ST++(1)			; -> CHI_I0
	.if	ROUND=1
	RND	CHI_R1
	.endif
	MPYF	K2,		CHI_R1,		CHI_R1		; * kappa^2
	.if	ROUND=1
	RND	CHI_R1
	.endif
	STF	CHI_R1,		*CHI_ST++(1)			; -> CHI_R1
	.if	ROUND=1
	RND	CHI_I1
	.endif
	MPYF	K2,		CHI_I1,		CHI_I1		; * kappa^2
	.if	ROUND=1
	RND	CHI_I1
	.endif
	STF	CHI_I1,		*CHI_ST++(1)			; -> CHI_I1
	.if	ROUND=1
	RND	CHI_R2
	.endif
	MPYF	K2,		CHI_R2,		CHI_R2		; * kappa^2
	.if	ROUND=1
	RND	CHI_R2
	.endif
	STF	CHI_R2,		*CHI_ST++(1)			; -> CHI_R2
	.if	ROUND=1
	RND	CHI_I2
	.endif
	MPYF	K2,		CHI_I2,		CHI_I2		; * kappa^2
	.if	ROUND=1
	RND	CHI_I2
	.endif
	STF	CHI_I2,		*CHI_ST++(1)			; -> CHI_I2

* Fourth spinor component
	MPYF	*OCHI++(1),	IK2,	CHI_R0		; OCHI / kappa_sq
	MPYF	*OCHI++(1),	IK2,	CHI_I0		; OCHI / kappa_sq
	MPYF	*OCHI++(1),	IK2,	CHI_R1		; OCHI / kappa_sq
	MPYF	*OCHI++(1),	IK2,	CHI_I1		; OCHI / kappa_sq
	MPYF	*OCHI++(1),	IK2,	CHI_R2		; OCHI / kappa_sq
	MPYF	*OCHI++(1),	IK2,	CHI_I2		; OCHI / kappa_sq
	
	SUBF	*CHI++(1),	CHI_R0				; - CHI
	SUBF	*CHI++(1),	CHI_I0				; - CHI
	SUBF	*CHI++(1),	CHI_R1				; - CHI
	SUBF	*CHI++(1),	CHI_I1				; - CHI
	SUBF	*CHI++(1),	CHI_R2				; - CHI
	SUBF	*CHI++(1),	CHI_I2				; - CHI

	SUBF	*AB0++(1),	CHI_I0				; - AB0
	ADDF	*AB0++(1),	CHI_R0				; + AB0
	SUBF	*AB0++(1),	CHI_I1				; - AB0
	ADDF	*AB0++(1),	CHI_R1				; + AB0
	SUBF	*AB0++(1),	CHI_I2				; - AB0
	ADDF	*AB0++(7),	CHI_R2				; + AB0

	CALL	_wfm_trick_segment_3			; trick segment 3

	.if	ROUND=1
	RND	CHI_R0
	.endif
	MPYF	K2,		CHI_R0,		CHI_R0		; * kappa^2
	.if	ROUND=1
	RND	CHI_R0
	.endif
	STF	CHI_R0,		*CHI_ST++(1)			; -> CHI_R0
	.if	ROUND=1
	RND	CHI_I0
	.endif
	MPYF	K2,		CHI_I0,		CHI_I0		; * kappa^2
	.if	ROUND=1
	RND	CHI_I0
	.endif
	STF	CHI_I0,		*CHI_ST++(1)			; -> CHI_I0
	.if	ROUND=1
	RND	CHI_R1
	.endif
	MPYF	K2,		CHI_R1,		CHI_R1		; * kappa^2
	.if	ROUND=1
	RND	CHI_R1
	.endif
	STF	CHI_R1,		*CHI_ST++(1)			; -> CHI_R1
	.if	ROUND=1
	RND	CHI_I1
	.endif
	MPYF	K2,		CHI_I1,		CHI_I1		; * kappa^2
	.if	ROUND=1
	RND	CHI_I1
	.endif
	STF	CHI_I1,		*CHI_ST++(1)			; -> CHI_I1
	.if	ROUND=1
	RND	CHI_R2
	.endif
	MPYF	K2,		CHI_R2,		CHI_R2		; * kappa^2
	.if	ROUND=1
	RND	CHI_R2
	.endif
	STF	CHI_R2,		*CHI_ST++(1)			; -> CHI_R2
	.if	ROUND=1
	RND	CHI_I2
	.endif
	MPYF	K2,		CHI_I2,		CHI_I2		; * kappa^2
	.if	ROUND=1
	RND	CHI_I2
	.endif
LX:	STF	CHI_I2,		*CHI_ST++(1)			; -> CHI_I2
								; end of X-loop
	LDI	@lcyzt,		LCYZT				; retrieve YZT loop counter
	SUBI	1,		LCYZT
	BNN	LYZT						; end of Y, Z, T loop

	POP	DP
*-------------------------------------------------------------------------

*
*  Restore the registers before returning to the C program                
*
	.if 	STAND_ALONE = 1
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

