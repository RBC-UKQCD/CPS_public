**--------------------------------------------------------------------
**  CVS keywords
**
**  $Author: chulwoo $
**  $Date: 2004-06-04 21:14:10 $
**  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_wilson/qcdsp_nos_hdw_diag/wfm_cmat_spproj.asm,v 1.3 2004-06-04 21:14:10 chulwoo Exp $
**  $Id: wfm_cmat_spproj.asm,v 1.3 2004-06-04 21:14:10 chulwoo Exp $
**  $Name: not supported by cvs2svn $
**  $Locker:  $
**  $Revision: 1.3 $
**  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_wilson/qcdsp_nos_hdw_diag/wfm_cmat_spproj.asm,v $
**  $State: Exp $
**
**--------------------------------------------------------------------
*****************************************************************************************
*---------------------------------------------------------------------------------------
*
* _wfm_cmat_spproj:
*
* If STAND_ALONE = 1 it can be called from C as:
*
* wfm_cmat_spproj(float *ab0, float *ab1, float *ab2, float *ab3,
*                 float *u, float *chi, float sign,
*                 Wilson *wilson_p, int cb)
*
* where:
*
* ab       --> 2 component spproj'ed padded spinors for the backward direction
* u        --> gauge fields (not padded)
* chi      --> 4 component spinors (not padded)
* sign     --> sign in spin project
* wilson_p --> pointer to Wilson struct
* cb       --> checkerboard
*
* If STAND_ALONE = 0 it can not be called from C. Instead the following is expected
* to be set up before the routine is called:
*
* 1) all the .ref in the STAND_ALONE = 0  case below must be defined
* 2) The following registers must be loaded as follows:
*    AR1  = cb
*    AR2  = wilson_p
*    AR6  = chi
*    AR7  = u
*     R6  = sign
*
*---------------------------------------------------------------------------------------
*****************************************************************************************

	.version	30

	.include	"../../include/wilson.hasm"
	.include	"wfm_macros.hasm"

	.def 	_wfm_cmat_spproj
	.def	_wfm_cmat_spproj_count

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
	.ref	tas1
	.ref	tas1_p
	.ref	c_u0_p
	.ref	c_u1_p
	.ref	c_u0
	.ref	c_u1

*---------------------------------------------------------------------------------------
* definitions
*---------------------------------------------------------------------------------------
WILSON_AD .set	AR2
CBB	.set	AR0
PTR	.set	AR1
MUPTR	.set	AR4
LCYZT	.set	AR5
AB	.set	AR3
AB0	.set	R7
AB1	.set	R7
AB2	.set	R7
AB3	.set	R7

FP	.set	AR3			; use FP for arguments only
U	.set	AR7			; stored in dma source address
TAS1	.set	AR2
C_U	.set	AR0
PSI	.set	AR6
DMABASE	.set	AR3

REAL	.set	R2
IMAG	.set	R3
SIGN	.set	R6
ZERO	.set	R5
COUNT	.set	R3
CODE_ADDR .set	R7
CACHE_CMAT .set	R1

LC	.set	AR7
P2R	.set	R1
P2I	.set	R2
P3R	.set 	R3	
P3I	.set	R4
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
	DEF_CB_MODE	cb_reg1, 72,  1, 10,  5, 10, 1, 0, 0, 1, 1
	DEF_CB_MODE	cb_reg2, 24, 12, 24,  8, 12, 1, 1, 1, 1, 1
	.endif

	.sect	"T:wfm0"
_wfm_cmat_spproj_count	.space	1

	.sect 	"T:wfm0" 
* Reserve space for ab base addresses.
muptr_p	.word	muptr
muptr	.space	4  

* Reserve space for DMA'd U in block 1.
u_p	.space	1
c_u_t	.space	1

* Space for code addresses
cache_cmat_p	.word	_cache_cmat
cmat_cont_p	.word	cmat_cont

* DMA stuff
dmabase	.word	0808000h		; base address for DMA control regs
dmagc	.set	0			; offset of DMA global ctrl reg
dmasrc	.set	4			; offset of DMA src addr reg
dmadest	.set	6			; offset of DMA dest addr reg
dmatc	.set	8			; offset of DMA transfer counter reg
dmagccf	.set	010001010011b		; START | INCSRC | INCDST | TC

* Cache stuff
cf	.word	0010000000000b		; cache freeze mask 
ce	.word	0100000000000b		; cache enable mask
cc	.word	1000000000000b		; cache clear mask

****************************************************************************************
* _wfm_cmat_spproj
****************************************************************************************
	.sect 	"T:wfm0"
_wfm_cmat_spproj:

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
	LDP	@ab0					; initialize DP for cram
	LDI     *-FP(2), 	TMP  			; Address of destination
	STI	TMP,		@ab0
	LDI     *-FP(3), 	TMP  			; Address of destination
	STI	TMP,		@ab1
	LDI     *-FP(4), 	TMP  			; Address of destination
	STI	TMP,		@ab2
	LDI     *-FP(5), 	TMP  			; Address of destination
	STI	TMP,		@ab3
	LDI     *-FP(6), 	U   			; Address of U
	LDI     *-FP(7), 	PSI 			; Address of PSI
	LDF     *-FP(8),	SIGN			; Sign in spin project (spproj)
	LDI     *-FP(9), 	WILSON_AD		; Wilson structure
	LDI     *-FP(10), 	PTR			; Checkerboard
	.endif

*---------------------------------------------------------------------------------------
* Routine starts here
*---------------------------------------------------------------------------------------
	PUSH	ST					; save status register
	PUSH    DP
	LDP	@ab0					; initialize DP for cram

	.if	CBUF=1
* Prepare addresses for Circular Buffer access.
* U   accesses are in mode   : 72,  1, 10,  5, 10, 1, 0, 0, 1, 1
* PSI accesses are in mode   : 24, 12, 24,  8, 12, 1, 1, 1, 1, 1
* All other accesses in mode :  1,  1,  1,  1,  1, 0, 0, 1, 1, 0
* Stores are through bank0 = 0h.
*  Set up circular buffer.
	LDI     @cb_cntrl_b, 	CBB			;  CBB has cbctrl address
	LDI     @cb_reg0, 	TMP
	STI     TMP, 		*+CBB(0)		;  Set flags for cntrl reg 0
	LDI     @cb_reg1,	TMP
	STI     TMP,		*+CBB(1)		;  Set flags for cntrl reg 1
	LDI     @cb_reg2,	TMP
	STI     TMP,		*+CBB(2)		;  Set flags for cntrl reg 2
* Put variables in the appropriate banks
	ADDI	@bank1a,	U
	ADDI	@bank2a,	PSI
	.endif

*
* Enable instruction cache
*
	ANDN	@cf, ST					; clear CF -- cache not frozen
	OR	@cc, ST					; set CC -- cache initially cleared
	OR	@ce, ST					; set CE -- cache enabled
;;	ANDN	@ce, ST					; set CE -- cache disabled

*
* Initialize DMA
*
	LDI	@dmabase, 	DMABASE			; DMABASE -> DMA control regs
	XOR	TMP,	  	TMP			; clear TMP
	STI	TMP,	  	*+DMABASE(dmagc)	; halt any current DMA operation

	LDI	U,		TMP
	ADDI	U_SIZE,		TMP			; next site of U
	STI	TMP,	  	*+DMABASE(dmasrc)	; DMA src = u[72]
	XOR	TMP,	 	TMP
	STI	TMP, 		*+DMABASE(dmatc)	; DMA TC = 0
	LDI	dmagccf,	TMP
	STI	TMP, 		*+DMABASE(dmagc)	; load DMA config and start DMA

	XOR	TMP,		TMP
	STI	TMP,		@_wfm_cmat_spproj_count	; DMA spin count

*
* Initialize addresses for CMAT_SPPROJ.
* U holds the address of the first element of   u
* TAS1 holds the address of the first element of  tas1 
* AB holds the address of the first element of  ab 
*
* Note: DP same for all of cram and is set above
*
	STI	U,		@u_p			; not used yet
	LDI	@c_u0_p,	C_U
	LDI	@c_u1_p,	TMP
	XOR	C_U,		TMP			; toggle for C_U
	STI	TMP,		@c_u_t	
	LDI	6,		IR1			; stride in spproj and cmat
	LDF	0.0,		ZERO			; Constant. Could be moved.

* Various addresses used for padding. 
	LSH	2,				PTR	; shift left by 2 and construct
	OR	2,				PTR	; the pointer into the PTR array
	MPYI	*+WILSON_AD(Wilson.offset),	PTR	; calculate and load the
	ADDI	*+WILSON_AD(Wilson.ptr),	PTR	; base address of the PTR array

	LDI	*+WILSON_AD(Wilson.yztmax),	LCYZT	; load the Y,Z,T loop counter

	LDI	@muptr_p,			MUPTR	; get base address of MUPTR
							; The auxiliary register MUPTR
							; addresses four words of
							; memory that contain the base
							; addresses of AB0,AB1,AB2,AB3

*
* Load the first U field directly
*
	LDF	*U++,		TMP
	RPTS	U_SIZE-2
	LDF	*U++,		TMP
||	STF	TMP,		*C_U++
       	STF	TMP,		*C_U--(U_SIZE-1)	; store and reset

* Toggle the C_U block to point to the next block
	XOR	@c_u_t,		C_U

*
* Start main loop
*
LYZT: 							; loop over Y, Z, T

	LDI	@ab0,		AB0
	ADDI	*PTR++(1),	AB0,		TMP	; add displacement to AB0
	LDI	@ab1,		AB1
	ADDI	*PTR++(1),	AB1,		TMP	; add displacement to AB1
||	STI	TMP,		*MUPTR++(1)		; store base address of AB0
	LDI	@ab2,		AB2
	ADDI	*PTR++(1),	AB2,		TMP	; add displacement to AB2
||	STI	TMP,		*MUPTR++(1)		; store base address of AB1
	LDI	@ab3,		AB3
	ADDI	*PTR++(1),	AB3,		TMP	; add displacement to AB3
||	STI	TMP,		*MUPTR++(1)		; store base address of AB2
	STI	TMP,		*MUPTR--(3)		; store base address of AB3
	LDI	*PTR++(1),	RC			; set the X loop counter

*
* Start X loop
*
	RPTB	LX					; loop over X

*
* Spin project (spproj) PSI into TAS1
*
	LDF	*PSI++(12),	TMP			; prepare first address of psi
	LDI	@tas1_p,	TAS1			; reload TAS1
	LDI	13,		IR0			; second stride
	LDI	2,		LC			; set loop counter


spproj:							; loop over colors
	MPYF	SIGN,		*PSI++(IR1),	P2R
	MPYF	SIGN,		*PSI++(1),	P3R
	MPYF	SIGN,		*PSI--(IR1),	P3I
	MPYF	SIGN,		*PSI--(IR0),	P2I
	SUBF	P3I,		*PSI++(IR1),	TMP
	.if	ROUND=1
	RND	TMP
	.endif
	SUBF	P2I,		*PSI--(IR1),	TMP
||	STF	TMP,		*TAS1++(IR1)
	.if	ROUND=1
	RND	TMP
	.endif
	SUBF	P3R,		*PSI++(IR1),	TMP
||	STF	TMP,		*TAS1++(IR1)
	.if	ROUND=1
	RND	TMP
	.endif
	ADDF	*PSI--(IR1),	P2R,		TMP
||	STF	TMP,		*TAS1++(IR1)
	.if	ROUND=1
	RND	TMP
	.endif
	SUBF	P2I,		*PSI++(IR1),	TMP
||	STF	TMP,		*TAS1++(IR1)
	.if	ROUND=1
	RND	TMP
	.endif
	ADDF	*PSI--(IR1),	P3I,		TMP
||	STF	TMP,		*TAS1++(IR1)
	.if	ROUND=1
	RND	TMP
	.endif
	ADDF	*PSI++(IR1),	P2R,		TMP
||	STF	TMP,		*TAS1++(IR1)
	.if	ROUND=1
	RND	TMP
	.endif
	ADDF	*PSI++(1),	P3R,		TMP
||	STF	TMP,		*TAS1++(IR1)
	.if	ROUND=1
	RND	TMP
	.endif
	ADDF	*PSI--(IR1),	P3I,		TMP
||	STF	TMP,		*TAS1++(1)
	.if	ROUND=1
	RND	TMP
	.endif
	ADDF	*PSI++(IR1),	P2I,		TMP
||	STF	TMP,		*TAS1--(IR1)
	.if	ROUND=1
	RND	TMP
	.endif
	SUBF	P3R,		*PSI--(IR1),	TMP
||	STF	TMP,		*TAS1--(IR1)
	.if	ROUND=1
	RND	TMP
	.endif
	ADDF	*PSI++(IR1),	P2R,		TMP
||	STF	TMP,		*TAS1--(IR1)
	.if	ROUND=1
	RND	TMP
	.endif
	ADDF	*PSI--(IR1),	P2I,		TMP
||	STF	TMP,		*TAS1--(IR1)
	.if	ROUND=1
	RND	TMP
	.endif
	SUBF	P3I,		*PSI++(IR1),	TMP
||	STF	TMP,		*TAS1--(IR1)
	.if	ROUND=1
	RND	TMP
	.endif

	.if	ROUND=0
	DBD	LC,		spproj			; delayed branch (if no round)
	.endif

	ADDF	*PSI--(IR1),	P2R,		TMP
||	STF	TMP,		*TAS1--(IR1)
	.if	ROUND=1
	RND	TMP
	.endif

	.if	ROUND=1
	DBD	LC,		spproj			; delayed branch (if round)
	.endif

	ADDF	*PSI++(IR0),	P3R,		TMP
||	STF	TMP,		*TAS1--(IR1)
	.if	ROUND=1
	RND	TMP
	.endif
	STF	TMP,		*TAS1++(1)
*-------------------------------------------------------; end of spproj loop

	ADDI	6,		PSI			; point to next site
	NOP	*TAS1--(6)				; point to beginning of block

* Gauge field successfully loaded.
* DMA load next U field. 
* NOTE: no error if read beyond end of field, but will not be used
	LDI	@dmabase,	DMABASE			; DMABASE -> DMA control regs
	STI	C_U, 		*+DMABASE(dmadest)	; DMA dest = c_u[0]
	LDI	U_SIZE,	 	TMP
	STI	TMP, 		*+DMABASE(dmatc)	; DMA TC = U_SIZE, start dma
;***	LDI	dmagccf,	TMP
;***	STI	TMP, 		*+DMABASE(dmagc)	; load DMA config and start DMA

* Switch C_U to other buffer
	XOR	@c_u_t,		C_U

*
* Start multiply 
* REAL,IMAG : accumulators for real/imaginary part
* TMP   : temporary storage for multiply
* ZERO  : register holding floating zero
*
* Loop over 4 directions
* Multiply by TAS1 by C_U^dag for 2 spin components
*
	LDI	2,		IR0
	LDI	3,		LC

	MPYF3	*C_U++(IR0),	*TAS1++(IR0),	TMP	; TMP  = Ur00 * Ar00

* Jump to CMAT routine
	LDI	@cache_cmat_p,	CACHE_CMAT		; dram cmat code
	LDI	@cmat_cont_p,	CODE_ADDR		; return address
	B	CACHE_CMAT				; jump to dram
cmat_cont:

	.if	ROUND=1					; if round less code fits in cache
	
	MPYF3	*C_U--(IR0),	*TAS1--(IR0),	TMP	; TMP = Ui21 * Ar21
||	ADDF3	TMP,	   	IMAG,		IMAG	; IMAG += TMP
	MPYF3	*C_U--(IR0),	*TAS1--(IR0),	TMP	; TMP = Ui11 * Ar11
||	SUBF3	TMP,	   	IMAG,		IMAG	; IMAG -= TMP
	MPYF3	*C_U,		*TAS1++,	TMP	; TMP = Ui01 * Ar11
||	SUBF3	TMP,	   	IMAG,		IMAG	; IMAG -= TMP

* Br11
	MPYF3	*C_U++(IR0),	*TAS1++(IR0),	TMP	; TMP = Ui01 * Ai01
||	SUBF3	TMP,	   	IMAG,		IMAG	; IMAG -= TMP
	MPYF3	*C_U++(IR0),	*TAS1++(IR0),	TMP	; TMP = Ui11 * Ai11
||	ADDF3	TMP, 	   	ZERO,		REAL	; REAL  = TMP

	MPYF3	*C_U--,		*TAS1--,	TMP	; TMP = Ui21 * Ai21
||	ADDF3	TMP,	   	REAL,		REAL	; REAL += TMP
	MPYF3	*C_U--(IR0),	*TAS1--(IR0),	TMP	; TMP = Ur21 * Ar21
||	ADDF3	TMP,	   	REAL,		REAL	; REAL += TMP
	MPYF3	*C_U--(IR0),	*TAS1--(IR0),	TMP	; TMP = Ur11 * Ar11
||	ADDF3	TMP,	   	REAL,		REAL	; REAL += TMP

	.endif

	MPYF3	*C_U--(IR1),	*TAS1++,	TMP	; TMP = Ur01 * Ar01
||	ADDF3	TMP,	   	REAL,		REAL	; REAL += TMP

* Bi01
	MPYF3	*C_U++(IR0),	*TAS1++(IR0),	TMP	; TMP = Ur00 * Ai01
||	ADDF3	TMP,	   	REAL,		REAL	; REAL += TMP

	.if	ROUND=1
	RND	REAL
	RND	IMAG
	.endif
	STF	REAL,	   	*AB--(IR0)		; Br11 = REAL
||	STF	IMAG,	   	*+AB			; Bi11 = IMAG

	MPYF3	*C_U++(IR0),	*TAS1++(IR0),	TMP	; TMP = Ur10 * Ai11
||	ADDF3	TMP, 	   	ZERO,		IMAG	; IMAG  = TMP
	MPYF3	*C_U++,		*TAS1--,	TMP	; TMP = Ur20 * Ai21
||	ADDF3	TMP,	   	IMAG,		IMAG	; IMAG += TMP
	MPYF3	*C_U--(IR0),	*TAS1--(IR0),	TMP	; TMP = Ui20 * Ar21
||	ADDF3	TMP,	   	IMAG,		IMAG	; IMAG += TMP
	MPYF3	*C_U--(IR0),	*TAS1--(IR0),	TMP	; TMP = Ui10 * Ar11
||	SUBF3	TMP,	   	IMAG,		IMAG	; IMAG -= TMP
	MPYF3	*C_U,		*TAS1++,	TMP	; TMP = Ui00 * Ar01
||	SUBF3	TMP,	   	IMAG,		IMAG	; IMAG -= TMP

* Br01
	MPYF3	*C_U++(IR0),	*TAS1++(IR0),	TMP	; TMP = Ui00 * Ai01
||	SUBF3	TMP,	   	IMAG,		IMAG	; IMAG -= TMP
	MPYF3	*C_U++(IR0),	*TAS1++(IR0),	TMP	; TMP = Ui10 * Ai11
||	ADDF3	TMP, 	   	ZERO,		REAL	; REAL  = TMP
	MPYF3	*C_U--,		*TAS1--,	TMP	; TMP = Ui20 * Ai21
||	ADDF3	TMP,	   	REAL,		REAL	; REAL += TMP
	MPYF3	*C_U--(IR0),	*TAS1--(IR0),	TMP	; TMP = Ur20 * Ar21
||	ADDF3	TMP,	   	REAL,		REAL	; REAL += TMP
	MPYF3	*C_U--(IR0),	*TAS1--(IR0),	TMP	; TMP = Ur10 * Ar11
||	ADDF3	TMP,	   	REAL,		REAL	; REAL += TMP
	MPYF3	*C_U,		*TAS1++(IR1),	TMP	; TMP = Ur00 * Ar01
||	ADDF3	TMP,	   	REAL,		REAL	; REAL += TMP

	NOP	*C_U++(18)				; jump to next dir

	.if	ROUND=0
	DBD	LC,		CACHE_CMAT		; branch back (if round)
	.endif

	MPYF3	*C_U++(IR0),	*TAS1++(IR0),	TMP	; TMP = Ur00 * Ar00
||	ADDF3	TMP,	   	REAL,		REAL	; REAL += TMP

	.if	ROUND=1
	RND	REAL
	DBD	LC,		CACHE_CMAT		; branch back (if no round)
	RND	IMAG
	.endif

	STF	REAL,	   	*AB++(IR1)		; Br01 = REAL
||	STF	IMAG,	   	*+AB			; Bi01 = IMAG

* Store modified address of  AB(mu). AB points to the would be next site
	STI	AB,		*MUPTR++(1)
*-------------------------------------------------------; end of cmat loop

* Reset C_U to the beginning of the block. The reset is IR0+U_SIZE
	NOP	*C_U--(74)

* Insure U field is loaded 
	LDI	@dmabase,	DMABASE			; DMABASE -> DMA control regs
	LDI	@_wfm_cmat_spproj_count,	COUNT
dmatest:
	ADDI	1,		COUNT
	LDI	*+DMABASE(dmatc), TMP			; inspect value of counter
	BNZ	dmatest					; spin until DMA is complete

	SUBI	1,		COUNT
	STI	COUNT,		@_wfm_cmat_spproj_count

* End of loops
LX:	NOP	*MUPTR--(4)				; end of X-loop

	DBD	LCYZT,		LYZT			; end of Y, Z, T loop
	NOP
	NOP
	NOP

	POP	DP
	POP	ST

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


;-------------------------------------------------------------------------
	.text

	.align

_cache_cmat:

* Load address of destination
	LDI	*MUPTR,		AB			; set the address

* Br00
	MPYF3	*C_U++(IR0),	*TAS1++(IR0),	TMP	; TMP = Ur10 * Ar10
||	ADDF3	TMP, 	   	ZERO,		REAL	; REAL += TMP
	MPYF3	*C_U++,		*TAS1++,	TMP	; TMP = Ur20 * Ar20
||	ADDF3	TMP, 	   	REAL,		REAL	; REAL += TMP
	MPYF3	*C_U--(IR0),	*TAS1--(IR0),	TMP	; TMP = Ui20 * Ai20
||	ADDF3	TMP,	   	REAL,		REAL	; REAL += TMP
	MPYF3	*C_U--(IR0),	*TAS1--(IR0),	TMP	; TMP = Ui10 * Ai10
||	ADDF3	TMP,	   	REAL,		REAL	; REAL += TMP
	MPYF3	*C_U,		*TAS1--,	TMP	; TMP = Ui00 * Ai00
||	ADDF3	TMP,	   	REAL,		REAL	; REAL += TMP

* Bi00
	MPYF3	*C_U++(IR0),	*TAS1++(IR0),	TMP	; TMP = Ui00 * Ar00
||	ADDF3	TMP,	   	REAL,		REAL	; REAL += TMP
	MPYF3	*C_U++(IR0),	*TAS1++(IR0),	TMP	; TMP = Ui10 * Ar10
||	ADDF3	TMP, 	   	ZERO,		IMAG	; IMAG  = TMP
	MPYF3	*C_U--,		*TAS1++,	TMP	; TMP = Ui20 * Ar20
||	ADDF3	TMP,	   	IMAG,		IMAG	; IMAG += TMP
	MPYF3	*C_U--(IR0),	*TAS1--(IR0),	TMP	; TMP = Ur20 * Ai20
||	ADDF3	TMP,	   	IMAG,		IMAG	; IMAG += TMP
	MPYF3	*C_U--(IR0),	*TAS1--(IR0),	TMP	; TMP = Ur10 * Ai10
||	SUBF3	IMAG,		TMP,		IMAG	; IMAG  = TMP - IMAG
	MPYF3	*C_U++(IR1),	*TAS1--,	TMP	; TMP = Ur00 * Ai00
||	ADDF3	TMP,		IMAG,		IMAG	; IMAG += TMP

* Br10
	MPYF3	*C_U++(IR0),	*TAS1++(IR0),	TMP	; TMP = Ur01 * Ar00
||	ADDF3	TMP, 		IMAG,		IMAG	; IMAG += TMP

	.if	ROUND=1
	RND	REAL
	RND	IMAG
	.endif
	STF	REAL,	   	*AB++(IR0)		; Br00 = REAL
||	STF	IMAG,	   	*+AB			; Bi00 = IMAG

	MPYF3	*C_U++(IR0),	*TAS1++(IR0),	TMP	; TMP = Ur11 * Ar10
||	ADDF3	TMP, 	   	ZERO,		REAL	; REAL  = TMP
	MPYF3	*C_U++,		*TAS1++,	TMP	; TMP = Ur21 * Ar20
||	ADDF3	TMP, 	   	REAL,		REAL	; REAL += TMP
	MPYF3	*C_U--(IR0),	*TAS1--(IR0),	TMP	; TMP = Ui21 * Ai20
||	ADDF3	TMP,	   	REAL,		REAL	; REAL += TMP
	MPYF3	*C_U--(IR0),	*TAS1--(IR0),	TMP	; TMP = Ui11 * Ai10
||	ADDF3	TMP,	   	REAL,		REAL	; REAL += TMP
	MPYF3	*C_U,		*TAS1--,	TMP	; TMP = Ui01 * Ai00
||	ADDF3	TMP,	   	REAL,		REAL	; REAL += TMP

* Bi10
	MPYF3	*C_U++(IR0),	*TAS1++(IR0),	TMP	; TMP = Ui01 * Ar00
||	ADDF3	TMP,	   	REAL,		REAL	; REAL += TMP
	MPYF3	*C_U++(IR0),	*TAS1++(IR0),	TMP	; TMP = Ui11 * Ar10
||	ADDF3	TMP, 	   	ZERO,		IMAG	; IMAG  = TMP
	MPYF3	*C_U--,		*TAS1++,	TMP	; TMP = Ui21 * Ar20
||	ADDF3	TMP,	   	IMAG,		IMAG	; IMAG += TMP
	MPYF3	*C_U--(IR0),	*TAS1--(IR0),	TMP	; TMP = Ur21 * Ai20
||	ADDF3	TMP,	   	IMAG,		IMAG	; IMAG += TMP
	MPYF3	*C_U--(IR0),	*TAS1--(IR0),	TMP	; TMP = Ur11 * Ai10
||	SUBF3	IMAG,	   	TMP,		IMAG	; IMAG  = IMAG - TMP
	MPYF3	*C_U++(IR1),	*TAS1--,	TMP	; TMP = Ur01 * Ai00
||	ADDF3	TMP,	   	IMAG,		IMAG	; IMAG += TMP

* Br20
	MPYF3	*C_U++(IR0),	*TAS1++(IR0),	TMP	; TMP = Ur02 * Ar00
||	ADDF3	TMP, 	   	IMAG,		IMAG	; IMAG += TMP

	.if	ROUND=1
	RND	REAL
	RND	IMAG
	.endif
	STF	REAL,	   	*AB++(IR0)		; Br10 = REAL
||	STF	IMAG,	   	*+AB			; Bi10 = IMAG

	MPYF3	*C_U++(IR0),	*TAS1++(IR0),	TMP	; TMP = Ur12 * Ar10
||	ADDF3	TMP, 	   	ZERO,		REAL	; REAL  = TMP
	MPYF3	*C_U++,		*TAS1++,	TMP	; TMP = Ur22 * Ar20
||	ADDF3	TMP, 	   	REAL,		REAL	; REAL += TMP
	MPYF3	*C_U--(IR0),	*TAS1--(IR0),	TMP	; TMP = Ui22 * Ai20
||	ADDF3	TMP,	   	REAL,		REAL	; REAL += TMP
	MPYF3	*C_U--(IR0),	*TAS1--(IR0),	TMP	; TMP = Ui21 * Ai10
||	ADDF3	TMP,	   	REAL,		REAL	; REAL += TMP
	MPYF3	*C_U,		*TAS1--,	TMP	; TMP = Ui02 * Ai00
||	ADDF3	TMP,	   	REAL,		REAL	; REAL += TMP

* Bi20
	MPYF3	*C_U++(IR0),	*TAS1++(IR0),	TMP	; TMP = Ui02 * Ar00
||	ADDF3	TMP,	   	REAL,		REAL	; REAL += TMP
	MPYF3	*C_U++(IR0),	*TAS1++(IR0),	TMP	; TMP = Ui12 * Ar10
||	ADDF3	TMP, 	   	ZERO,		IMAG	; IMAG  = TMP
	MPYF3	*C_U--,		*TAS1++,	TMP	; TMP = Ui22 * Ar20
||	ADDF3	TMP,	   	IMAG,		IMAG	; IMAG += TMP
	MPYF3	*C_U--(IR0),	*TAS1--(IR0),	TMP	; TMP = Ur22 * Ai20
||	ADDF3	TMP,	   	IMAG,		IMAG	; IMAG += TMP
	MPYF3	*C_U--(IR0),	*TAS1--(IR0),	TMP	; TMP = Ur12 * Ai10
||	SUBF3	IMAG,	   	TMP,		IMAG	; IMAG  = TMP - IMAG
	MPYF3	*C_U,		*TAS1++(IR1),	TMP	; TMP = Ur02 * Ai00
||	ADDF3	TMP,	   	IMAG,		IMAG	; IMAG += TMP

* Reverse direction
* Bi21
	MPYF3	*C_U++(IR0),	*TAS1++(IR0),	TMP	; TMP = Ur02 * Ai01
||	ADDF3	TMP,	   	IMAG,		IMAG	; IMAG += TMP

	.if	ROUND=1
	RND	REAL
	RND	IMAG
	.endif
	STF	REAL,	   	*AB++(IR1)		; Br20 = REAL
||	STF	IMAG,	   	*+AB			; Bi20 = IMAG

	MPYF3	*C_U++(IR0),	*TAS1++(IR0),	TMP	; TMP = Ur12 * Ai11
||	ADDF3	TMP, 	   	ZERO,		IMAG	; IMAG  = TMP
	MPYF3	*C_U++,		*TAS1--,	TMP	; TMP = Ur22 * Ai21
||	ADDF3	TMP,	   	IMAG,		IMAG	; IMAG += TMP
	MPYF3	*C_U--(IR0),	*TAS1--(IR0),	TMP	; TMP = Ui22 * Ar21
||	ADDF3	TMP,	   	IMAG,		IMAG	; IMAG += TMP
	MPYF3	*C_U--(IR0),	*TAS1--(IR0),	TMP	; TMP = Ui12 * Ar11
||	SUBF3	TMP,	   	IMAG,		IMAG	; IMAG -= TMP
	MPYF3	*C_U,		*TAS1++,	TMP	; TMP = Ui02 * Ar11
||	SUBF3	TMP,	   	IMAG,		IMAG	; IMAG -= TMP

* Br21
	MPYF3	*C_U++(IR0),	*TAS1++(IR0),	TMP	; TMP = Ui02 * Ai01
||	SUBF3	TMP,	   	IMAG,		IMAG	; IMAG -= TMP
	MPYF3	*C_U++(IR0),	*TAS1++(IR0),	TMP	; TMP = Ui12 * Ai11
||	ADDF3	TMP, 	   	ZERO,		REAL	; REAL  = TMP
	MPYF3	*C_U--,		*TAS1--,	TMP	; TMP = Ui22 * Ai21
||	ADDF3	TMP,	   	REAL,		REAL	; REAL += TMP
	MPYF3	*C_U--(IR0),	*TAS1--(IR0),	TMP	; TMP = Ur22 * Ar21
||	ADDF3	TMP,	   	REAL,		REAL	; REAL += TMP
	MPYF3	*C_U--(IR0),	*TAS1--(IR0),	TMP	; TMP = Ur12 * Ar11
||	ADDF3	TMP,	   	REAL,		REAL	; REAL += TMP
	MPYF3	*C_U--(IR1),	*TAS1++,	TMP	; TMP = Ur02 * Ar01
||	ADDF3	TMP,	   	REAL,		REAL	; REAL += TMP

* Bi11
	MPYF3	*C_U++(IR0),	*TAS1++(IR0),	TMP	; TMP = Ur01 * Ai01
||	ADDF3	TMP,	   	REAL,		REAL	; REAL += TMP

	.if	ROUND=1
	RND	REAL
	RND	IMAG
* Branch delayed back to cram (if round = 1 less code fits in cache)
	BD	CODE_ADDR
	.endif

	STF	REAL,	   	*AB--(IR0)		; Br21 = REAL
||	STF	IMAG,	   	*+AB			; Bi21 = IMAG

	MPYF3	*C_U++(IR0),	*TAS1++(IR0),	TMP	; TMP = Ur11 * Ai11
||	ADDF3	TMP, 	   	ZERO,		IMAG	; IMAG  = TMP
	MPYF3	*C_U++,		*TAS1--,	TMP	; TMP = Ur21 * Ai21
||	ADDF3	TMP,	   	IMAG,		IMAG	; IMAG += TMP

	.if	ROUND=0

	MPYF3	*C_U--(IR0),	*TAS1--(IR0),	TMP	; TMP = Ui21 * Ar21
||	ADDF3	TMP,	   	IMAG,		IMAG	; IMAG += TMP
	MPYF3	*C_U--(IR0),	*TAS1--(IR0),	TMP	; TMP = Ui11 * Ar11
||	SUBF3	TMP,	   	IMAG,		IMAG	; IMAG -= TMP
	MPYF3	*C_U,		*TAS1++,	TMP	; TMP = Ui01 * Ar11
||	SUBF3	TMP,	   	IMAG,		IMAG	; IMAG -= TMP

* Br11
	MPYF3	*C_U++(IR0),	*TAS1++(IR0),	TMP	; TMP = Ui01 * Ai01
||	SUBF3	TMP,	   	IMAG,		IMAG	; IMAG -= TMP
	MPYF3	*C_U++(IR0),	*TAS1++(IR0),	TMP	; TMP = Ui11 * Ai11
||	ADDF3	TMP, 	   	ZERO,		REAL	; REAL  = TMP

* Branch delayed back to cram (if round = 0 more code fits in cache)
	BD	CODE_ADDR

	MPYF3	*C_U--,		*TAS1--,	TMP	; TMP = Ui21 * Ai21
||	ADDF3	TMP,	   	REAL,		REAL	; REAL += TMP
	MPYF3	*C_U--(IR0),	*TAS1--(IR0),	TMP	; TMP = Ur21 * Ar21
||	ADDF3	TMP,	   	REAL,		REAL	; REAL += TMP
	MPYF3	*C_U--(IR0),	*TAS1--(IR0),	TMP	; TMP = Ur11 * Ar11
||	ADDF3	TMP,	   	REAL,		REAL	; REAL += TMP

	.endif



	.end

