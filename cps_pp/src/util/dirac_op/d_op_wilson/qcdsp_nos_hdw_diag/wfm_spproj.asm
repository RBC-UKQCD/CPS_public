**--------------------------------------------------------------------
**  CVS keywords
**
**  $Author: mcneile $
**  $Date: 2003-06-22 13:34:46 $
**  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_wilson/qcdsp_nos_hdw_diag/wfm_spproj.asm,v 1.1.1.1 2003-06-22 13:34:46 mcneile Exp $
**  $Id: wfm_spproj.asm,v 1.1.1.1 2003-06-22 13:34:46 mcneile Exp $
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
**  $RCSfile: wfm_spproj.asm,v $
**  $Revision: 1.1.1.1 $
**  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_wilson/qcdsp_nos_hdw_diag/wfm_spproj.asm,v $
**  $State: Exp $
**
**--------------------------------------------------------------------
****************************************************************************************
*---------------------------------------------------------------------------------------
* _wfm_spproj
*
* It performs the spin projection by equating the top 2 spin components of
* [(1 +/- gamma_mu)/2 psi] to af_mu. The bottom two components can be reconstructed 
* using the trick routines.
*
* If STAND_ALONE = 1 it can be called from C as:
*
* wfm_spproj(float *af0,                 ; af for mu = 0 
*	      float *af1,                 ; af for mu = 1 
*	      float *af2,                 ; af for mu = 2 
*   	      float *af3,                 ; af for mu = 3 
* 	      float *psi,		  ; 4 comp. spinor
*	      float sign, 		  ; the +/- sign
*	      Wilson *wilson_p,           ; Wilson struct.
*	      int cb);		          ; checkerboard 0/1 = even/odd
*
* If STAND_ALONE = 0 it can not be called from C. Instead the following is expected
* to be set up before the routine is called:
*
* 1) all the .ref in the STAND_ALONE = 0  case below must be defined
* 2) The following registers must be loaded as follows:
*    AR6 <-- psi pointer
*    R5  <-- sign 
*    AR0 <-- Wilson pointer
*    AR1 <-- cb
*
*---------------------------------------------------------------------------------------
****************************************************************************************
	.version	30

	.include	"../../include/wilson.hasm"
	.include	"wfm_macros.hasm"

	.def	 _wfm_spproj

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

*---------------------------------------------------------------------------------------
* definitions
*---------------------------------------------------------------------------------------
FP	.set	AR3
WILSON_AD .set	AR0
CBB	.set	AR0
PTR	.set	AR1
AF0	.set	AR0
AF1	.set	AR2
AF2	.set	AR3
AF3	.set	AR4
PSI	.set	AR6
LCC	.set	AR7
SGN	.set	AR5
TMP	.set	R0
PSI0R	.set	R0
PSI0I	.set	R1
PSI1R	.set	R2
PSI1I	.set	R3
PSI2R	.set	R4
PSI2I	.set	R5
PSI3R	.set	R6
PSI3I	.set	R7
SIGN	.set	R5
LCYZT	.set	R7

*---------------------------------------------------------------------------------------
*  Set up circular buffer modes.
*---------------------------------------------------------------------------------------
 	.if	CBUF=1
	.sect	"T:wfm1"
*  DEF_CB_MODE	.macro	mode_name, wcntset, gnteeset, lookcnt, highrqset, hysterhigh, 
*				   wcnt_enable, kick_enable, clrbit, waitgntee, zflag
	DEF_CB_MODE	cb_reg0,  1,  1,  1,  1,  1, 0, 0, 1, 1, 0
	DEF_CB_MODE	cb_reg1, 24,  5, 24, 24, 24, 1, 1, 0, 1, 1
	.endif

*---------------------------------------------------------------------------------------
* Reserve space for local data
*---------------------------------------------------------------------------------------
	.sect 	"T:wfm0"
lcyzt	.space	1
psi0r	.space	1
sgn	.word	sign
sign	.space	1



****************************************************************************************
* _wfm_spproj
****************************************************************************************
	.sect 	"T:wfm1"
_wfm_spproj:

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
	LDI     *-FP(6), 	PSI 			; Address of PSI
	LDF     *-FP(7),	SIGN			; Sign in spin project (spproj)
	LDI     *-FP(8), 	WILSON_AD		; Wilson structure
	LDI     *-FP(9), 	PTR			; Checkerboard
	.endif

*---------------------------------------------------------------------------------------
* Routine starts here
*---------------------------------------------------------------------------------------
	PUSH    DP
	LDP	@af0

	LDI	5,		IR0
	LDI	17,		IR1
	LDI	2,		LCC			; initialize color loop counter
	LDI	@sgn,		SGN			; set up sign so it can be accessed
	.if	ROUND=1
	RND	SIGN
	.endif
	STF	SIGN,		@sign			; indirectly

* Various addresses used for padding. Note, OR'ing with 0 is a nop
	LSH	2,				PTR	; shift left by 2 and construct
;;;	OR	0,				PTR	; the pointer into the PTR array
	MPYI	*+WILSON_AD(Wilson.offset),	PTR	; calculate and load the
	ADDI	*+WILSON_AD(Wilson.ptr),	PTR	; base address of the PTR array
	LDI	*+WILSON_AD(Wilson.yztmax),	LCYZT	; load the Y,Z,T loop counter

*
* AF holds the address of the first element of  af 
*
* Note: DP same for all of cram and is set above
*
 	.if	CBUF=1
* Prepare addresses for Circular Buffer access.
* PSI accesses are in mode       :  24,  5, 24, 24, 24, 1, 0, 1, 1, 1
* All other accesses are in mode :   1,  1,  1,  1,  1, 0, 0, 1, 1, 0
* The sub-banks must alternate to start the next process. 
* Stores must be through bank0 = 0h. AF0, AF1, AF2, AF3 are in bank0.
* The PSI accesses use kickstart which initiates a new process even
* from the same bank. Another method would be not to use kickstart
* but set the clrbt to 1 so that the indermidiate stores will invalidate
* the current process. This however will result to extra cycles every time
* around the color loop.
*
*  Set up circular buffer.
	LDI     @cb_cntrl_b, 	CBB			;  AR0 has cbctrl address
	LDI     @cb_reg0, 	TMP
	STI     TMP, 		*+CBB(0)		;  Set flags for cntrl reg 0
	LDI     @cb_reg1,	TMP
	STI     TMP,		*+CBB(1)		;  Set flags for cntrl reg 1
* Put variables in the appropriate banks
 	ADDI	@bank1a,	PSI			; PSI is controled from CB register 1
 	.endif

*
* Start main loop
*
LYZT: 							; loop over Y, Z, T
	STI	LCYZT,		@lcyzt
	LDI	@af0,		TMP
	ADDI	*PTR++(1),	TMP,	AF0		; add displacement to AF0
	LDI	@af1,		TMP
	ADDI	*PTR++(1),	TMP,	AF1		; add displacement to AF1
	LDI	@af2,		TMP
	ADDI	*PTR++(1),	TMP,	AF2		; add displacement to AF2
	LDI	@af3,		TMP
	ADDI	*PTR++(1),	TMP,	AF3		; add displacement to AF3
	LDI	*PTR++(1),	RC			; set the X loop counter
*
* Start X loop
*
	RPTB	LX					; loop over X
*
* Spin project (spproj) PSI into AF
*
* Color Loop
	LDF	*PSI,		TMP			; kick start circular buffer
LCTR:							; color loop
	LDF	*PSI++(1),	PSI0R
	LDF	*PSI++(5),	PSI0I
	LDF	*PSI++(1),	PSI1R
	LDF	*PSI++(5),	PSI1I
	MPYF	*PSI++(1),	*SGN,		PSI2R
	MPYF	*PSI++(IR0),	*SGN,		PSI2I
	MPYF	*PSI++(1),	*SGN,		PSI3R
	MPYF	*PSI--(IR1),	*SGN,		PSI3I

	.if	ROUND=1
	RND	PSI0R
	.endif
	STF	PSI0R,		@psi0r

	ADDF	PSI3I,		PSI0R,		TMP
	.if	ROUND=1
	RND	TMP
	.endif
	STF	TMP,		*AF0++(1)
	SUBF	PSI3R,		PSI0I,		TMP
	.if	ROUND=1
	RND	TMP
	.endif
	STF	TMP,		*AF0++(5)
	ADDF	PSI2I,		PSI1R,		TMP
	.if	ROUND=1
	RND	TMP
	.endif
	STF	TMP,		*AF0++(1)
	SUBF	PSI2R,		PSI1I,		TMP
	.if	ROUND=1
	RND	TMP
	.endif
	STF	TMP,		*AF0--(5)
	LDF	@psi0r,		PSI0R

	ADDF	PSI3R,		PSI0R,		TMP
	.if	ROUND=1
	RND	TMP
	.endif
	STF	TMP,		*AF1++(1)
	ADDF	PSI3I,		PSI0I,		TMP
	.if	ROUND=1
	RND	TMP
	.endif
	STF	TMP,		*AF1++(5)
	SUBF	PSI2R,		PSI1R,		TMP
	.if	ROUND=1
	RND	TMP
	.endif
	STF	TMP,		*AF1++(1)
	SUBF	PSI2I,		PSI1I,		TMP
	.if	ROUND=1
	RND	TMP
	.endif
	STF	TMP,		*AF1--(5)
	LDF	@psi0r,		PSI0R

	ADDF	PSI2I,		PSI0R,		TMP
	.if	ROUND=1
	RND	TMP
	.endif
	STF	TMP,		*AF2++(1)
	SUBF	PSI2R,		PSI0I,		TMP
	.if	ROUND=1
	RND	TMP
	.endif
	STF	TMP,		*AF2++(5)
	SUBF	PSI3I,		PSI1R,		TMP
	.if	ROUND=1
	RND	TMP
	.endif
	STF	TMP,		*AF2++(1)
	ADDF	PSI3R,		PSI1I,		TMP
	.if	ROUND=1
	RND	TMP
	.endif
	STF	TMP,		*AF2--(5)
	LDF	@psi0r,		PSI0R

	SUBF	PSI2R,		PSI0R,		TMP
	.if	ROUND=1
	RND	TMP
	.endif
	STF	TMP,		*AF3++(1)
	SUBF	PSI2I,		PSI0I,		TMP
	.if	ROUND=1
	RND	TMP
	.endif
	STF	TMP,		*AF3++(5)
	SUBF	PSI3R,		PSI1R,		TMP

	.if	ROUND=0
	DBD	LCC,		LCTR			; branch to LCTR (if round=0)
	.endif

	.if	ROUND=1
	RND	TMP
	.endif
	STF	TMP,		*AF3++(1)

	.if	ROUND=1
	DBD	LCC,		LCTR			; branch to LCTR (if round=1)
	.endif

	SUBF	PSI3I,		PSI1I,		TMP
	.if	ROUND=1
	RND	TMP
	.endif
	STF	TMP,		*AF3--(5)

*-------------------------------------------------------; end LCTR loop



	LDI	2,	LCC				; initialize color loop counter
	ADDI	18,	PSI
	ADDI	6,	AF0
	ADDI	6,	AF1
	ADDI	6,	AF2
LX:	ADDI	6,	AF3
						; end of X-loop
	LDI	@lcyzt,		LCYZT
	SUBI	1,		LCYZT
	BNN	LYZT				; end of Y, Z, T loop


	POP    DP
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
