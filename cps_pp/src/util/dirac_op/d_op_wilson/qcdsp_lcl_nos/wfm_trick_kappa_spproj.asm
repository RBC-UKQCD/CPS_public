**--------------------------------------------------------------------
**  CVS keywords
**
**  $Author: mcneile $
**  $Date: 2003-06-22 13:34:46 $
**  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_wilson/qcdsp_lcl_nos/wfm_trick_kappa_spproj.asm,v 1.1.1.1 2003-06-22 13:34:46 mcneile Exp $
**  $Id: wfm_trick_kappa_spproj.asm,v 1.1.1.1 2003-06-22 13:34:46 mcneile Exp $
**  $Name: not supported by cvs2svn $
**  $Locker:  $
**  $Log: not supported by cvs2svn $
**  Revision 1.2  2001/06/19 18:12:59  anj
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
**  $RCSfile: wfm_trick_kappa_spproj.asm,v $
**  $Revision: 1.1.1.1 $
**  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_wilson/qcdsp_lcl_nos/wfm_trick_kappa_spproj.asm,v $
**  $State: Exp $
**
**--------------------------------------------------------------------
*****************************************************************************************
*---------------------------------------------------------------------------------------
*
* _wfm_trick_kappa_spproj:
*
* If STAND_ALONE = 1 it can be called from C as:
*
*   wfm_trick_kappa_spproj(float *phi, 
*			   float *af0, float *af1, float *af2, float *af3,
*                          float *psi, 
*			   float *ab0, float *ab1, float *ab2, float *ab3,
*			   float *mp_sq_p,
*		 	   float -kappa_sq, float 1/kappa_sq, 
*			   Wilson *wilson_p, int cb, int skip_spproj);
*
* where:
*
* phi        --> 4 component spinor (not padded)
* af         --> 2 component spproj'ed padded spinors for the forward direction
* psi        --> 4 component spinor (not padded)
* ab         --> 2 component spproj'ed padded spinors for the backward direction
* mp_sq_p    --> pointer where the sum of [M * psi]^2 is stored 
* -kappa_sq  --> -kappa^2, kappa is Wilson's parameter
* 1/kappa_sq --> 1/kappa^2, kappa is Wilson's parameter
* wilson_p   --> pointer to Wilson struct
* cb         --> checkerboard
* skip_spproj --> if =1 it skips the sum and spproj section of the routine 
*                 (i.e. it only executes the trick-kappa part).
*                 if =0 the fill trick-kappa-spproj routine is executed.
*
*
* If STAND_ALONE = 0 it can not be called from C. Instead the following is expected
* to be set up before the routine is called:
*
* 1) all the .ref in the STAND_ALONE = 0  case below must be defined
* 2) The following registers must be loaded as follows:
*    AR0  = wilson_p
*    AR4  = phi
*    AR5  = psi
*    AR7  = cb
*    R7   = skip_spproj
*
*---------------------------------------------------------------------------------------
*****************************************************************************************

	.version	30

	.include	"../../include/wilson.hasm"
	.include	"wfm_macros.hasm"

	.def _wfm_trick_kappa_spproj

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
	.ref	mp_sq_p
 	.ref	af0
 	.ref	af1
 	.ref	af2
 	.ref	af3
  	.ref	ab0
 	.ref	ab1
 	.ref	ab2
 	.ref	ab3
  	.ref	tpsi1_p
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

LCYZT	.set	R7			; to be preserved by the core
LCC	.set	AR7
SSPR	.set	R7

K2	.set	R6
IK2	.set	R7

PTRF	.set	AR6
PTRB	.set	AR7
WILSON_AD .set	AR0
CBB	.set	AR0

AB0	.set	AR0
AB1	.set	AR1
AB2	.set	AR2
AB3	.set	AR3

AF0	.set	AR0
AF1	.set	AR1
AF2	.set	AR2
AF3	.set	AR3

CHI	.set	AR4
OCHI	.set	AR5
PSI	.set	AR6
CHI_ST	.set	AR7

CHI_R0	.set	R0
CHI_I0	.set	R1
CHI_R1	.set	R2
CHI_I1	.set	R3
CHI_R2	.set	R4
CHI_I2	.set	R5

PSI_SQ	.set	R0
MP_SQ	.set	R2
ACCUM	.set	R5


TMP1	.set	R4

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
*  DEF_CB_MODE	.macro	mode_name, wcntset, gnteeset, lookcnt, highrqset, hysterhigh, 
*				   wcnt_enable, kick_enable, clrbit, waitgntee, zflag
	DEF_CB_MODE	cb_reg0,  1,  1,  1,  1,  1, 0, 0, 1, 1, 0
	DEF_CB_MODE	cb_reg1,  6,  1,  6,  6,  6, 1, 0, 1, 1, 1
	.endif

	.sect 	"T:wfm1" 
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

ptrb	.space	1
ptrf	.space	1

* Reserve space for chi store address
chi_str	.space	1

* Reserve space for the YZT loop counter
lcyzt	.space	1
* Reserve space for the component loop counter
lc	.space	1

* Reserve space for skip_spproj
skip_spproj	.space	1


****************************************************************************************
* _wfm_trick_kappa_spproj
****************************************************************************************
	.sect 	"T:wfm0"
_wfm_trick_kappa_spproj

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
	LDI     *-FP(2), 	CHI  			; Address of destination
	LDI     *-FP(3), 	TMP 			; Address of AF0
	STI	TMP,		@af0
	LDI     *-FP(4), 	TMP 			; Address of AF1
	STI	TMP,		@af1
	LDI     *-FP(5),	TMP 			; Address of AF2
	STI	TMP,		@af2
	LDI     *-FP(6), 	TMP 			; Address of AF3
	STI	TMP,		@af3
	LDI     *-FP(7), 	OCHI  			; Address of old Chi
	LDI     *-FP(8), 	TMP 			; Address of AB0
	STI	TMP,		@ab0
	LDI     *-FP(9), 	TMP 			; Address of AB1
	STI	TMP,		@ab1
	LDI     *-FP(10),	TMP 			; Address of AB2
	STI	TMP,		@ab2
	LDI     *-FP(11), 	TMP 			; Address of AB3
	STI	TMP,		@ab3
	LDI     *-FP(12), 	TMP 			; Address of AB3
	STI	TMP,		@mp_sq_p
	LDF     *-FP(13),	TMP			; Pointer to M*Psi^2
	STF	TMP,		@m_kappa_sq		; -kappa_sq
	NEGF	TMP,		TMP
	RND	TMP
	STF	TMP,		@kappa_sq		; kappa_sq
	LDF     *-FP(14),	TMP			;
	STF	TMP,		@i_kappa_sq		; 1 / kappa_sq 
	LDI     *-FP(15), 	WILSON_AD		; Wilson structure
	LDI     *-FP(16), 	PTRB			; Checkerboard of AB 
	LDI     *-FP(17), 	SSPR			; skip_spproj 
	.endif

*---------------------------------------------------------------------------------------
* Routine starts here
*---------------------------------------------------------------------------------------
	PUSH    DP
	LDP	@ab0

	STI	SSPR,		@skip_spproj		; save the skip_spproj

	LDI	1, 		TMP
	SUBI	PTRB,		TMP,		PTRF	; Checkerboard of AF 
	LDF	@kappa_sq,	K2			; put k^2 to a register
	STI	CHI,		@chi_str		; address to stotre chi

* Initialize mp_sq to 0
	LDF	0.0,		MP_SQ
	LDI	@mp_sq_p,	AR3
	STF	MP_SQ,		*AR3

* Various addresses used for padding
* For trick:
	LSH	2,				PTRB	; shift left by 2 and construct
	OR	3,				PTRB	; the pointer into the PTRB array
	MPYI	*+WILSON_AD(Wilson.offset),	PTRB	; calculate and load the
	ADDI	*+WILSON_AD(Wilson.ptr),	PTRB	; base address of the PTRB array
* For spproj:
	LSH	2,				PTRF	; shift left by 2 and construct
;;;	OR	0,				PTRF	; the pointer into the PTRF array
	MPYI	*+WILSON_AD(Wilson.offset),	PTRF	; calculate and load the
	ADDI	*+WILSON_AD(Wilson.ptr),	PTRF	; base address of the PTRF array
	STI	PTRB,				@ptrb   ; save register
	STI	PTRF,				@ptrf   ; save register
	LDI	*+WILSON_AD(Wilson.yztmax),	LCYZT	; load the Y,Z,T loop counter

*
* Initialize addresses for TRICK.
* AB holds the address of the first element of  ab 
*
* Note: DP same for all of cram and is set above
*
	LDI	6,		IR0			; set first stride of spproj
	LDI	13,		IR1			; set second stride ; of spproj

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
	LDI	@ptrb,		PTRB			; retrieve PTRB
	LDI	@ab0_t,		TMP
	ADDI	*PTRB++(1),	TMP,	AB0		; add displacement to AB0
	LDI	@ab1_t,		TMP
	ADDI	*PTRB++(1),	TMP,	AB1		; add displacement to AB1
	LDI	@ab2_t,		TMP
	ADDI	*PTRB++(1),	TMP,	AB2		; add displacement to AB2
	LDI	@ab3_t,		TMP
	ADDI	*PTRB++(1),	TMP,	AB3		; add displacement to AB3
	LDI	*PTRB++(1),	RC			; set the X loop counter
	STI	PTRB,		@ptrb			; save register	

	STI	AB0,		@tab0			; save AB0 address
	STI	AB1,		@tab1			; save AB1 address
	STI	AB2,		@tab2			; save AB2 address
	STI	AB3,		@tab3			; save AB3 address

	LDI	@ptrf,		PTRF			; retrieve PTRF
	LDI	@af0,		TMP
	ADDI	*PTRF++(1),	TMP,	AF0		; add displacement to AF0
	LDI	@af1,		TMP
	ADDI	*PTRF++(1),	TMP,	AF1		; add displacement to AF1
	LDI	@af2,		TMP
	ADDI	*PTRF++(1),	TMP,	AF2		; add displacement to AF2
	LDI	@af3,		TMP
	ADDI	*PTRF++(1),	TMP,	AF3		; add displacement to AF3
	LDI	*PTRF++(1),	RC			; set the X loop counter
	STI	PTRF,		@ptrf			; save register	

	STI	AF0,		@taf0			; save AF0 address
	STI	AF1,		@taf1			; save AF1 address
	STI	AF2,		@taf2			; save AF2 address
	STI	AF3,		@taf3			; save AF3 address

	LDI	@tpsi1_p,	PSI			; Initialize address for the CRAM PSI

	LDF	@i_kappa_sq,	IK2				; put 1/k^2 to a register

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

	LDI	@chi_str,	CHI_ST				; Initialize address to store CHI in DRAM 
	LDI	@tab0,		AB0				; Retrive AB0 address
	LDI	@tab1,		AB1				; Retrive AB1 address
	LDI	@tab2,		AB2				; Retrive AB2 address
	LDI	@tab3,		AB3				; Retrive AB3 address


* First and second spinor components
	LDI	0,		TMP
	STI	TMP,		@lc
LCCM:
	CALL	_wfm_trick_segment_1				; trick segment 1

	.if	ROUND=1
	RND	CHI_R0
	.endif
	MPYF	K2,		CHI_R0				; * kappa^2
	.if	ROUND=1
	RND	CHI_R0
	.endif
	STF	CHI_R0,		*PSI++(1)			; -> CHI_R0  CRAM 
||	STF	CHI_R0,		*CHI_ST++(1)			; -> CHI_R0  DRAM
	.if	ROUND=1
	RND	CHI_I0
	.endif
	MPYF	K2,		CHI_I0				; * kappa^2
	.if	ROUND=1
	RND	CHI_I0
	.endif
	STF	CHI_I0,		*PSI++(1)			; -> CHI_I0  CRAM 
||	STF	CHI_I0,		*CHI_ST++(1)			; -> CHI_I0  DRAM
	.if	ROUND=1
	RND	CHI_R1
	.endif
	MPYF	K2,		CHI_R1				; * kappa^2
	.if	ROUND=1
	RND	CHI_R1
	.endif
	STF	CHI_R1,		*PSI++(1)			; -> CHI_R1  CRAM 
||	STF	CHI_R1,		*CHI_ST++(1)			; -> CHI_R1  DRAM
	.if	ROUND=1
	RND	CHI_I1
	.endif
	MPYF	K2,		CHI_I1				; * kappa^2
	.if	ROUND=1
	RND	CHI_I1
	.endif
	STF	CHI_I1,		*PSI++(1)			; -> CHI_I1  CRAM 
||	STF	CHI_I1,		*CHI_ST++(1)			; -> CHI_I1  DRAM
	.if	ROUND=1
	RND	CHI_R2
	.endif
	MPYF	K2,		CHI_R2				; * kappa^2
	.if	ROUND=1
	RND	CHI_R2
	.endif
	STF	CHI_R2,		*PSI++(1)			; -> CHI_R2  CRAM 
||	STF	CHI_R2,		*CHI_ST++(1)			; -> CHI_R2  DRAM
	.if	ROUND=1
	RND	CHI_I2
	.endif
	MPYF	K2,		CHI_I2				; * kappa^2

	.if	ROUND=0
	LDI	@lc,		TMP				
	BNND	LCCM						; branch to LCCM (if round=0)
	.endif

	.if	ROUND=1
	RND	CHI_I2
	LDI	@lc,		TMP				
	BNND	LCCM						; branch to LCCM (if round=1)
	.endif

	STF	CHI_I2,		*PSI++(1)			; -> CHI_I2  CRAM 
||	STF	CHI_I2,		*CHI_ST++(1)			; -> CHI_I2  DRAM
	SUBI	1,		TMP
	STI	TMP,		@lc
*---------------------------------------------------------------; end of LCCM loop

	SUBI	6,		AB0
	SUBI	6,		AB1
	SUBI	12,		AB2
	SUBI	12,		AB3

* Third spinor component
	NEGF	K2,		K2				;  K2 = -kappa_sq
	RND	K2
	NEGF	IK2,		IK2				; IK2 = -1/kappa_sq
	RND	IK2
	MPYF	IK2,		*OCHI++(1),	CHI_R0		; -OCHI / kappa_sq
	MPYF	IK2,		*OCHI++(1),	CHI_I0		; -OCHI / kappa_sq
	MPYF	IK2,		*OCHI++(1),	CHI_R1		; -OCHI / kappa_sq
	MPYF	IK2,		*OCHI++(1),	CHI_I1		; -OCHI / kappa_sq
	MPYF	IK2,		*OCHI++(1),	CHI_R2		; -OCHI / kappa_sq
	MPYF	IK2,		*OCHI++(1),	CHI_I2		; -OCHI / kappa_sq
	
	ADDF	*CHI++(1),	CHI_R0				; + CHI
	ADDF	*CHI++(1),	CHI_I0				; + CHI
	ADDF	*CHI++(1),	CHI_R1				; + CHI
	ADDF	*CHI++(1),	CHI_I1				; + CHI
	ADDF	*CHI++(1),	CHI_R2				; + CHI
	ADDF	*CHI++(1),	CHI_I2				; + CHI

	SUBF	*AB0++(1),	CHI_I0				; + AB0
	ADDF	*AB0++(1),	CHI_R0				; - AB0
	SUBF	*AB0++(1),	CHI_I1				; + AB0
	ADDF	*AB0++(1),	CHI_R1				; - AB0
	SUBF	*AB0++(1),	CHI_I2				; + AB0
	ADDF	*AB0--(11),	CHI_R2				; - AB0

	CALL	_wfm_trick_segment_2				; trick segment 2

	.if	ROUND=1
	RND	CHI_R0
	.endif
	MPYF	K2,		CHI_R0				; * kappa^2
	.if	ROUND=1
	RND	CHI_R0
	.endif
	STF	CHI_R0,		*PSI++(1)			; -> CHI_R0  CRAM 
||	STF	CHI_R0,		*CHI_ST++(1)			; -> CHI_R0  DRAM
	.if	ROUND=1
	RND	CHI_I0
	.endif
	MPYF	K2,		CHI_I0				; * kappa^2
	.if	ROUND=1
	RND	CHI_I0
	.endif
	STF	CHI_I0,		*PSI++(1)			; -> CHI_I0  CRAM 
||	STF	CHI_I0,		*CHI_ST++(1)			; -> CHI_I0  DRAM
	.if	ROUND=1
	RND	CHI_R1
	.endif
	MPYF	K2,		CHI_R1				; * kappa^2
	.if	ROUND=1
	RND	CHI_R1
	.endif
	STF	CHI_R1,		*PSI++(1)			; -> CHI_R1  CRAM 
||	STF	CHI_R1,		*CHI_ST++(1)			; -> CHI_R1  DRAM
	.if	ROUND=1
	RND	CHI_I1
	.endif
	MPYF	K2,		CHI_I1				; * kappa^2
	.if	ROUND=1
	RND	CHI_I1
	.endif
	STF	CHI_I1,		*PSI++(1)			; -> CHI_I1  CRAM 
||	STF	CHI_I1,		*CHI_ST++(1)			; -> CHI_I1  DRAM
	.if	ROUND=1
	RND	CHI_R2
	.endif
	MPYF	K2,		CHI_R2				; * kappa^2
	.if	ROUND=1
	RND	CHI_R2
	.endif
	STF	CHI_R2,		*PSI++(1)			; -> CHI_R2  CRAM 
||	STF	CHI_R2,		*CHI_ST++(1)			; -> CHI_R2  DRAM
	.if	ROUND=1
	RND	CHI_I2
	.endif
	MPYF	K2,		CHI_I2				; * kappa^2
	.if	ROUND=1
	RND	CHI_I2
	.endif
	STF	CHI_I2,		*PSI++(1)			; -> CHI_I2  CRAM 
||	STF	CHI_I2,		*CHI_ST++(1)			; -> CHI_I2  DRAM

* Fourth spinor component
	MPYF	IK2,		*OCHI++(1),	CHI_R0		; -OCHI / kappa_sq
	MPYF	IK2,		*OCHI++(1),	CHI_I0		; -OCHI / kappa_sq
	MPYF	IK2,		*OCHI++(1),	CHI_R1		; -OCHI / kappa_sq
	MPYF	IK2,		*OCHI++(1),	CHI_I1		; -OCHI / kappa_sq
	MPYF	IK2,		*OCHI++(1),	CHI_R2		; -OCHI / kappa_sq
	MPYF	IK2,		*OCHI++(1),	CHI_I2		; -OCHI / kappa_sq
	
	ADDF	*CHI++(1),	CHI_R0				; + CHI
	ADDF	*CHI++(1),	CHI_I0				; + CHI
	ADDF	*CHI++(1),	CHI_R1				; + CHI
	ADDF	*CHI++(1),	CHI_I1				; + CHI
	ADDF	*CHI++(1),	CHI_R2				; + CHI
	ADDF	*CHI++(1),	CHI_I2				; + CHI

	SUBF	*AB0++(1),	CHI_I0				; + AB0
	ADDF	*AB0++(1),	CHI_R0				; - AB0
	SUBF	*AB0++(1),	CHI_I1				; + AB0
	ADDF	*AB0++(1),	CHI_R1				; - AB0
	SUBF	*AB0++(1),	CHI_I2				; + AB0
	ADDF	*AB0++(7),	CHI_R2				; - AB0

	CALL	_wfm_trick_segment_3				; trick segment 3

	.if	ROUND=1
	RND	CHI_R0
	.endif
	MPYF	K2,		CHI_R0				; * kappa^2
	.if	ROUND=1
	RND	CHI_R0
	.endif
	STF	CHI_R0,		*PSI++(1)			; -> CHI_R0  CRAM 
||	STF	CHI_R0,		*CHI_ST++(1)			; -> CHI_R0  DRAM
	.if	ROUND=1
	RND	CHI_I0
	.endif
	MPYF	K2,		CHI_I0				; * kappa^2
	.if	ROUND=1
	RND	CHI_I0
	.endif
	STF	CHI_I0,		*PSI++(1)			; -> CHI_I0  CRAM 
||	STF	CHI_I0,		*CHI_ST++(1)			; -> CHI_I0  DRAM
	.if	ROUND=1
	RND	CHI_R1
	.endif
	MPYF	K2,		CHI_R1				; * kappa^2
	.if	ROUND=1
	RND	CHI_R1
	.endif
	STF	CHI_R1,		*PSI++(1)			; -> CHI_R1  CRAM 
||	STF	CHI_R1,		*CHI_ST++(1)			; -> CHI_R1  DRAM
	.if	ROUND=1
	RND	CHI_I1
	.endif
	MPYF	K2,		CHI_I1				; * kappa^2
	.if	ROUND=1
	RND	CHI_I1
	.endif
	STF	CHI_I1,		*PSI++(1)			; -> CHI_I1  CRAM 
||	STF	CHI_I1,		*CHI_ST++(1)			; -> CHI_I1  DRAM
	.if	ROUND=1
	RND	CHI_R2
	.endif
	MPYF	K2,		CHI_R2				; * kappa^2
	.if	ROUND=1
	RND	CHI_R2
	.endif
	STF	CHI_R2,		*PSI++(1)			; -> CHI_R2  CRAM 
||	STF	CHI_R2,		*CHI_ST++(1)			; -> CHI_R2  DRAM
	.if	ROUND=1
	RND	CHI_I2
	.endif
	MPYF	K2,		CHI_I2				; * kappa^2
	.if	ROUND=1
	RND	CHI_I2
	.endif
	STF	CHI_I2,		*PSI++(1)			; -> CHI_I2  CRAM 
||	STF	CHI_I2,		*CHI_ST++(1)			; -> CHI_I2  DRAM

	SUBI	24,		PSI				; Reset PSI address
	STI	CHI_ST,		@chi_str			; Save CHI_STore address
	STI	AB0,		@tab0				; Save AB0 address
	STI	AB1,		@tab1				; Save AB1 address
	STI	AB2,		@tab2				; Save AB2 address
	STI	AB3,		@tab3				; Save AB3 address
	LDF	@kappa_sq,	K2				; K2 = kappa_sq
	LDF	@i_kappa_sq,	IK2				; IK2 = 1/kappa_sq

	LDI	@skip_spproj,	TMP				; check if you should
	BNZ	LX						; skip spproj 
*---------------------------------------------------------------; and sum

* On processor sum used to calculate the global sum needed by the CG
* Square PSI and accumulate onto mp_sq 
	LDI	8,		ACCUM
	LDI	@mp_sq_p,	AR3
	LDF	*AR3,		MP_SQ

	MPYF3	*PSI,		*PSI++(1),	PSI_SQ		; Initial multiply

accum:
	SUBI	1,		ACCUM
	BNZD	accum

	MPYF3	*PSI,		*PSI++(1),	PSI_SQ
||	ADDF3	MP_SQ,		PSI_SQ,		MP_SQ
	MPYF3	*PSI,		*PSI++(1),	PSI_SQ
||	ADDF3	MP_SQ,		PSI_SQ,		MP_SQ
	MPYF3	*PSI,		*PSI++(1),	PSI_SQ
||	ADDF3	MP_SQ,		PSI_SQ,		MP_SQ


	RND	MP_SQ						; Round the sum
	STF	MP_SQ,		*AR3				; Save intermediate sum

	SUBI	13,		PSI				; Reset PSI address
*---------------------------------------------------------------; end sum

******************************************************************************
*			spproj
******************************************************************************
* Spin project (spproj) PSI into AF
	LDI	@taf0,		AF0				; Retrive AF0 address
	LDI	@taf1,		AF1				; Retrive AF1 address
	LDI	@taf2,		AF2				; Retrive AF2 address
	LDI	@taf3,		AF3				; Retrive AF3 address
	LDI	2,		LCC				; initialize color loop counter
LCSP:								; color loop
	NEGF	*PSI++(IR0),	P2R
	NEGF	*PSI++,		P3R
	NEGF	*PSI--(IR0),	P3I
	NEGF	*PSI--(IR1),	P2I

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

	ADDI	IR0,		AF0				; point to next site
	ADDI	IR0,		AF1				; point to next site
	ADDI	IR0,		AF2				; point to next site
	ADDI	IR0,		AF3				; point to next site
	SUBI	12,		PSI				; Reset PSI address

	STI	AF0,		@taf0				; Save AF0 address
	STI	AF1,		@taf1				; Save AF1 address
	STI	AF2,		@taf2				; Save AF2 address
	STI	AF3,		@taf3				; Save AF3 address

LX:	NOP
*---------------------------------------------------------------; end of X-loop

	LDI	@lcyzt,		LCYZT				; retrieve loop counter
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

