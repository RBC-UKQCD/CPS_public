**--------------------------------------------------------------------
**  CVS keywords
**
**  $Author: mcneile $
**  $Date: 2003-06-22 13:34:46 $
**  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_wilson/qcdsp_nos/wfm_trick_spproj.asm,v 1.1.1.1 2003-06-22 13:34:46 mcneile Exp $
**  $Id: wfm_trick_spproj.asm,v 1.1.1.1 2003-06-22 13:34:46 mcneile Exp $
**  $Name: not supported by cvs2svn $
**  $Locker:  $
**  $Log: not supported by cvs2svn $
**  Revision 1.2  2001/06/19 18:13:05  anj
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
**  $RCSfile: wfm_trick_spproj.asm,v $
**  $Revision: 1.1.1.1 $
**  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_wilson/qcdsp_nos/wfm_trick_spproj.asm,v $
**  $State: Exp $
**
**--------------------------------------------------------------------
*****************************************************************************************
*---------------------------------------------------------------------------------------
*
* _wfm_trick_spproj:
*
* If STAND_ALONE = 1 it can be called from C as:
*
*  	   wfm_trick_spproj(chi, 
*			    float *af0, float *af1, float *af2, float *af3,
*			    float *ab0, float *ab1, float *ab2, float *ab3,
*		 	     float sign, Wilson *wilson_p, int cb, int skip_spproj);
*
* where:
*
* chi      --> 4 component spinors (not padded)
* af       --> 2 component spproj'ed padded spinors for the forward direction
* ab       --> 2 component spproj'ed padded spinors for the backward direction
* sign     --> sign in spin project
* wilson_p --> pointer to Wilson struct
* cb       --> checkerboard
* skip_spproj --> if =1 it skips the spproj section of the routine 
*                 (i.e. it only executes the trick part).
*                 if =0 the fill trick-spproj routine is executed.
*
* If STAND_ALONE = 0 it can not be called from C. Instead the following is expected
* to be set up before the routine is called:
*
* 1) all the .ref in the STAND_ALONE = 0  case below must be defined
* 2) The following registers must be loaded as follows:
*    AR0  = wilson_p
*    AR4  = chi
*    AR7  = cb
*    R6   = sign
*    R7   = skip_spproj
*
*---------------------------------------------------------------------------------------
*****************************************************************************************

	.version	30

	.include	"../../include/wilson.hasm"
	.include	"wfm_macros.hasm"

	.def	 _wfm_trick_spproj

	.ref	_wfm_trick_segment_1
	.ref	_wfm_trick_segment_2
	.ref	_wfm_trick_segment_3
	.ref	_wfm_spproj_segment

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
	.ref	tpsi0_p
	.ref	ab0_t
	.ref	ab1_t
	.ref	ab2_t
	.ref	ab3_t

*---------------------------------------------------------------------------------------
* definitions
*---------------------------------------------------------------------------------------
FP	.set	AR3			; use FP for arguments only

PTRB	.set	AR7
PTRF	.set	AR5
LCYZT	.set	R7
LCC	.set	AR7
SIGN	.set	R6
SSPR	.set	R7
CBB	.set	AR0
WILSON_AD .set	AR0

AB0	.set	AR0
AB1	.set	AR1
AB2	.set	AR2
AB3	.set	AR3

AF0	.set	AR0
AF1	.set	AR1
AF2	.set	AR2
AF3	.set	AR3

CHI	.set	AR4
PSI	.set	AR6
CHI_ST	.set	AR5

CHI_R0	.set	R0
CHI_I0	.set	R1
CHI_R1	.set	R2
CHI_I1	.set	R3
CHI_R2	.set	R4
CHI_I2	.set	R5

TMP	.set	R0

P2R	.set	R1
P2I	.set	R2
P3R	.set	R3
P3I	.set	R4

*---------------------------------------------------------------------------------------
* Local data and pointers
*---------------------------------------------------------------------------------------
*  Set up circular buffer modes.
 	.if	CBUF=1
	.sect	"T:wfm0"
	DEF_CB_MODE	cb_reg0,  1,  1,  1,  1,  1, 0, 0, 1, 1, 0
	DEF_CB_MODE	cb_reg1,  6,  1,  6,  6,  6, 1, 0, 1, 1, 1
	.endif

	.sect	"T:wfm1"
* Reserve space for AB addresses.  
tab0	.space	1
tab1	.space	1
tab2	.space	1
tab3	.space	1

* Reserve space for AF addresses.  
taf0	.space	1
taf1	.space	1
taf2	.space	1
taf3	.space	1

* Reserve space for PTRB address
ptrb	.space	1
ptrf	.space	1

* Reserve space for CHI_STore address
chi_st	.space	1

* Reserve space for skip_spproj
skip_spproj	.space	1

****************************************************************************************
* _wfm_trick_spproj
****************************************************************************************
	.sect 	"T:wfm0"
_wfm_trick_spproj:

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
	LDI     *-FP(3), 	TMP 				; Address of AF0
	STI	TMP,		@af0
	LDI     *-FP(4), 	TMP 				; Address of AF1
	STI	TMP,		@af1
	LDI     *-FP(5),	TMP 				; Address of AF2
	STI	TMP,		@af2
	LDI     *-FP(6), 	TMP 				; Address of AF3
	STI	TMP,		@af3
	LDI     *-FP(7), 	TMP 				; Address of AB0
	STI	TMP,		@ab0
	LDI     *-FP(8), 	TMP 				; Address of AB1
	STI	TMP,		@ab1
	LDI     *-FP(9),	TMP 				; Address of AB2
	STI	TMP,		@ab2
	LDI     *-FP(10), 	TMP 				; Address of AB3
	STI	TMP,		@ab3
	LDF     *-FP(11),	SIGN				; Sign in trick
	LDI     *-FP(12), 	WILSON_AD			; Wilson structure
	LDI     *-FP(13), 	PTRB				; Checkerboard of AB 
	LDI     *-FP(14), 	SSPR				; skip_spproj 
	.endif

*---------------------------------------------------------------------------------------
* Routine starts here
*---------------------------------------------------------------------------------------
	PUSH    DP
	LDP	@ab0

	STI	SSPR,		@skip_spproj			; save the skip_spproj

	LDI	1, 		TMP
	SUBI 	PTRB,		TMP,		PTRF		; Checkerboard of AF 
	STI	CHI,		@chi_st				; Address to store CHI in DRAM

* Various addresses used for padding
* For trick:
	LSH	2,				PTRB		; shift left by 2 and construct
	OR	3,				PTRB		; the pointer into the PTRB array
	MPYI	*+WILSON_AD(Wilson.offset),	PTRB		; calculate and load the
	ADDI	*+WILSON_AD(Wilson.ptr),	PTRB		; base address of the PTRB array
* For spproj:
	LSH	2,				PTRF		; shift left by 2 and construct
;;;	OR	0,				PTRF		; the pointer into the PTRF array
	MPYI	*+WILSON_AD(Wilson.offset),	PTRF		; calculate and load the
	ADDI	*+WILSON_AD(Wilson.ptr),	PTRF		; base address of the PTRF array
	LDI	*+WILSON_AD(Wilson.yztmax),	LCYZT		; load the Y,Z,T loop counter
	STI	PTRB,				@ptrb   	; save register
	STI	PTRF,				@ptrf   	; save register

* Initialize address for the CRAM PSI
	LDI	@tpsi0_p,	PSI
*
* Initialize addresses for TRICK.
* AB holds the address of the first element of  ab 
*
* Note: DP same for all of cram and is set above
*
	LDI	6,		IR0				; set first stride of spproj
	LDI	13,		IR1				; set second stride of spproj

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

	LDI	@ptrb,		PTRB				; retrieve PTRB
	LDI	@ab0_t,		TMP
	ADDI	*PTRB++(1),	TMP,		AB0		; add displacement to AB0
	LDI	@ab1_t,		TMP
	ADDI	*PTRB++(1),	TMP,		AB1		; add displacement to AB1
	LDI	@ab2_t,		TMP
	ADDI	*PTRB++(1),	TMP,		AB2		; add displacement to AB2
	LDI	@ab3_t,		TMP
	ADDI	*PTRB++(1),	TMP,		AB3		; add displacement to AB3
	LDI	*PTRB++(1),	RC				; set the X loop counter
	STI	PTRB,		@ptrb				; save register	

	STI	AB0,		@tab0
	STI	AB1,		@tab1
	STI	AB2,		@tab2
	STI	AB3,		@tab3

	LDI	@ptrf,		PTRF				; retrieve PTRF
	LDI	@af0,		TMP
	ADDI	*PTRF++(1),	TMP,		AF0		; add displacement to AF0
	LDI	@af1,		TMP
	ADDI	*PTRF++(1),	TMP,		AF1		; add displacement to AF1
	LDI	@af2,		TMP
	ADDI	*PTRF++(1),	TMP,		AF2		; add displacement to AF2
	LDI	@af3,		TMP
	ADDI	*PTRF++(1),	TMP,		AF3		; add displacement to AF3
	LDI	*PTRF++(1),	RC				; set the X loop counter
	STI	PTRF,		@ptrf				; save register	

	STI	AF0,		@taf0
	STI	AF1,		@taf1
	STI	AF2,		@taf2
	STI	AF3,		@taf3

	LDI	@chi_st,	CHI_ST				; Initialize address to store CHI in DRAM 

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
	LDI	@tab0,		AB0
	LDI	@tab1,		AB1
	LDI	@tab2,		AB2
	LDI	@tab3,		AB3

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
	STF	CHI_R0,		*PSI++(1)			; -> CHI_R0  CRAM 
||	STF	CHI_R0,		*CHI_ST++(1)			; -> CHI_R0  DRAM
	.if	ROUND=1
	RND	CHI_I0
	.endif
	STF	CHI_I0,		*PSI++(1)			; -> CHI_I0  CRAM 
||	STF	CHI_I0,		*CHI_ST++(1)			; -> CHI_I0  DRAM
	.if	ROUND=1
	RND	CHI_R1
	.endif
	STF	CHI_R1,		*PSI++(1)			; -> CHI_R1  CRAM 
||	STF	CHI_R1,		*CHI_ST++(1)			; -> CHI_R1  DRAM
	.if	ROUND=0
	DBD	LCC,		LCCM				; branch to LCCM (if round=0)
	.endif
	.if	ROUND=1
	RND	CHI_I1
	.endif
	STF	CHI_I1,		*PSI++(1)			; -> CHI_I1  CRAM 
||	STF	CHI_I1,		*CHI_ST++(1)			; -> CHI_I1  DRAM
	.if	ROUND=1
	RND	CHI_R2
	DBD	LCC,		LCCM				; branch to LCCM (if round=1)
	.endif
	STF	CHI_R2,		*PSI++(1)			; -> CHI_R2  CRAM 
||	STF	CHI_R2,		*CHI_ST++(1)			; -> CHI_R2  DRAM
	.if	ROUND=1
	RND	CHI_I2
	.endif
	STF	CHI_I2,		*PSI++(1)			; -> CHI_I2  CRAM 
||	STF	CHI_I2,		*CHI_ST++(1)			; -> CHI_I2  DRAM

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
	STF	CHI_R0,		*PSI++(1)			; -> CHI_R0  CRAM 
||	STF	CHI_R0,		*CHI_ST++(1)			; -> CHI_R0  DRAM
	.if	ROUND=1
	RND	CHI_I0
	.endif
	STF	CHI_I0,		*PSI++(1)			; -> CHI_I0  CRAM 
||	STF	CHI_I0,		*CHI_ST++(1)			; -> CHI_I0  DRAM
	.if	ROUND=1
	RND	CHI_R1
	.endif
	STF	CHI_R1,		*PSI++(1)			; -> CHI_R1  CRAM 
||	STF	CHI_R1,		*CHI_ST++(1)			; -> CHI_R1  DRAM
	.if	ROUND=1
	RND	CHI_I1
	.endif
	STF	CHI_I1,		*PSI++(1)			; -> CHI_I1  CRAM 
||	STF	CHI_I1,		*CHI_ST++(1)			; -> CHI_I1  DRAM
	.if	ROUND=1
	RND	CHI_R2
	.endif
	STF	CHI_R2,		*PSI++(1)			; -> CHI_R2  CRAM 
||	STF	CHI_R2,		*CHI_ST++(1)			; -> CHI_R2  DRAM
	.if	ROUND=1
	RND	CHI_I2
	.endif
	STF	CHI_I2,		*PSI++(1)			; -> CHI_I2  CRAM 
||	STF	CHI_I2,		*CHI_ST++(1)			; -> CHI_I2  DRAM

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
	STF	CHI_R0,		*PSI++(1)			; -> CHI_R0  CRAM 
||	STF	CHI_R0,		*CHI_ST++(1)			; -> CHI_R0  DRAM
	.if	ROUND=1
	RND	CHI_I0
	.endif
	STF	CHI_I0,		*PSI++(1)			; -> CHI_I0  CRAM 
||	STF	CHI_I0,		*CHI_ST++(1)			; -> CHI_I0  DRAM
	.if	ROUND=1
	RND	CHI_R1
	.endif
	STF	CHI_R1,		*PSI++(1)			; -> CHI_R1  CRAM 
||	STF	CHI_R1,		*CHI_ST++(1)			; -> CHI_R1  DRAM
	.if	ROUND=1
	RND	CHI_I1
	.endif
	STF	CHI_I1,		*PSI++(1)			; -> CHI_I1  CRAM 
||	STF	CHI_I1,		*CHI_ST++(1)			; -> CHI_I1  DRAM
	.if	ROUND=1
	RND	CHI_R2
	.endif
	STF	CHI_R2,		*PSI++(1)			; -> CHI_R2  CRAM 
||	STF	CHI_R2,		*CHI_ST++(1)			; -> CHI_R2  DRAM
	.if	ROUND=1
	RND	CHI_I2
	.endif
	STF	CHI_I2,		*PSI++(1)			; -> CHI_I2  CRAM 
||	STF	CHI_I2,		*CHI_ST++(1)			; -> CHI_I2  DRAM

	SUBI	12,		PSI				; Reset PSI address

	STI	AB0,		@tab0
	STI	AB1,		@tab1
	STI	AB2,		@tab2
	STI	AB3,		@tab3

	LDI	@skip_spproj,	TMP				; check if you 
	BNZ	LX						; should skip spproj

******************************************************************************
*			spproj
******************************************************************************

* Spin project (spproj) PSI into AF

	LDI	@taf0,		AF0
	LDI	@taf1,		AF1
	LDI	@taf2,		AF2
	LDI	@taf3,		AF3

	LDI	2,		LCC				; initialize color loop counter

LCSP:								; color loop
	MPYF	SIGN,		*PSI++(IR0),	P2R
	MPYF	SIGN,		*PSI++,		P3R
	MPYF	SIGN,		*PSI--(IR0),	P3I
	MPYF	SIGN,		*PSI--(IR1),	P2I

	CALL	_wfm_spproj_segment				; spproj segment

	.if	ROUND=0
	DBD	LCC,		LCSP				; branch to LCSP (if round=0)
	.endif
	.if	ROUND=1
	RND	TMP
	.endif
	SUBF	P3I,		*PSI--(IR0),	TMP		; TMP  = P1I - P3I
||	STF	TMP,		*AF3++				; AF3 += TMP
	.if	ROUND=1
	RND	TMP
	DBD	LCC,		LCSP				; branch to LCSP (if round=1)
	.endif
	SUBF	P2I,		*PSI++(IR1),	TMP		; TMP  = P0I - P2I
||	STF	TMP,		*AF3--(IR0)			; AF3 += TMP
	.if	ROUND=1
	RND	TMP
	.endif
	STF	TMP,		*AF3++			
*---------------------------------------------------------------; end of LCSP color loop
								; end color loop
	ADDI	IR0,		AF0
	ADDI	IR0,		AF1
	ADDI	IR0,		AF2
	ADDI	IR0,		AF3
	SUBI	6,		PSI				; Reset PSI address (partially)

	STI	AF0,		@taf0
	STI	AF1,		@taf1
	STI	AF2,		@taf2
	STI	AF3,		@taf3

LX:	SUBI	12,		PSI				; Reset PSI address

*---------------------------------------------------------------; end of X loop

	STI	CHI_ST,		@chi_st				; Save CHI_STore address
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
