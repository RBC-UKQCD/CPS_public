**--------------------------------------------------------------------
**  CVS keywords
**
**  $Author: chulwoo $
**  $Date: 2004-06-04 21:14:11 $
**  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_wilson/qcdsp_nos_hdw_diag/wfm_mat_trick.asm,v 1.3 2004-06-04 21:14:11 chulwoo Exp $
**  $Id: wfm_mat_trick.asm,v 1.3 2004-06-04 21:14:11 chulwoo Exp $
**  $Name: not supported by cvs2svn $
**  $Locker:  $
**  $Revision: 1.3 $
**  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_wilson/qcdsp_nos_hdw_diag/wfm_mat_trick.asm,v $
**  $State: Exp $
**
**--------------------------------------------------------------------
*****************************************************************************************
*---------------------------------------------------------------------------------------
*
* _wfm_mat_trick:
*
* If STAND_ALONE = 1 it can be called from C as:
*
* sublattice_mat_trick(float *chi, float *u, 
*                      float *af0, float *af1, float *af2, float *af3,
*                      float sign, Wilson *wilson_p, int shftptr)
*
* where:
*
* chi      --> 4 component spinors (not padded)
* u        --> gauge fields (not padded)
* af       --> 2 component spproj'ed padded spinors for the backward direction
* sign     --> sign in spin project
* wilson_p --> pointer to Wilson structure
* cb       --> chekerboard
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

	.def 	_wfm_mat_trick
	.def	_wfm_mat_trick_count

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
	.ref	c_u0_p
	.ref	c_u1_p
	.ref	c_u0
	.ref	c_u1
	.ref	tas0
	.ref	tas0_p
	.ref	tas1
	.ref	tas1_p

*---------------------------------------------------------------------------------------
* definitions
*---------------------------------------------------------------------------------------
WILSON_AD .set	AR2
CBB	.set	AR0
PTR	.set	AR1
MUPTR	.set	AR5
LCYZT	.set	R7
AF	.set	AR3
AF0	.set	R1
AF1	.set	R1
AF2	.set	R1
AF3	.set	R1

FP	.set	AR3					; use FP for arguments only
U	.set	AR7					; temporary
TAS0	.set	AR2
TAS1	.set	AR4
AX	.set	AR4
C_U	.set	AR0
CHI	.set	AR6
DMABASE	.set	AR3					; temporary

SIGN_P	.set	AR3					; temporary

REAL	.set	R2
IMAG	.set	R3
SIGN	.set	R6
ZERO	.set	R5
COUNT	.set	R3
CODE_ADDR .set	R4
CACHE_MAT .set	R1

LC	.set	AR7					; temporary
CHIR0	.set	R1
CHIR1	.set	R2
CHIR2	.set 	R3	
CHIR3	.set	R4
CHII0	.set	R0
CHII1	.set	R3
CHII2	.set 	R2	
CHII3	.set	R1
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
	DEF_CB_MODE	cb_reg2, 12,  1, 12, 12, 12, 1, 0, 1, 1, 1
	.endif

	.sect	"T:wfm0"
_wfm_mat_trick_count	.space	1

	.sect 	"T:wfm0" 
* Researve space for the sign
sign_p	.word	sign
sign	.space	2

* Reserve space for temporary AF.  
af_t	.word	af0_t
af0_t	.space	1
af1_t	.space	1
af2_t	.space	1
af3_t	.space	1

* Reserve space for temporary AF to be used to set the circular buffer access.  
af0_tmp	.word	1
af1_tmp	.word	1
af2_tmp	.word	1
af3_tmp	.word	1

* Reserve space for DMA'd U in block 1.
u_p	.space	1
c_u_t	.space	1

* Space for code addresses
cache_mat_p	.word	_cache_mat
cont_p	.word	cont

* DMA stuff
dmabase	.word	0808000h				; base address for DMA control regs
dmagc	.set	0					; offset of DMA global ctrl reg
dmasrc	.set	4					; offset of DMA src addr reg
dmadest	.set	6					; offset of DMA dest addr reg
dmatc	.set	8					; offset of DMA transfer counter reg
dmagccf	.set	010001010011b				; START | INCSRC | INCDST | TC

* Cache stuff
cf	.word	0010000000000b				; cache freeze mask 
ce	.word	0100000000000b				; cache enable mask
cc	.word	1000000000000b				; cache clear mask

****************************************************************************************
* _wfm_mat_trick
****************************************************************************************
	.sect 	"T:wfm0"
_wfm_mat_trick:

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
	LDI     *-FP(2), 	CHI  			; Address of destination
	LDI     *-FP(3), 	U   			; Address of U
	LDI     *-FP(4), 	TMP 			; Address of AF0
	STI	TMP,		@af0
	LDI     *-FP(5), 	TMP 			; Address of AF1
	STI	TMP,		@af1
	LDI     *-FP(6),	TMP 			; Address of AF2
	STI	TMP,		@af2
	LDI     *-FP(7), 	TMP 			; Address of AF3
	STI	TMP,		@af3
	LDF     *-FP(8),	SIGN			; Sign in trick
	LDI     *-FP(9), 	WILSON_AD		; Wilson structure
	LDI     *-FP(10), 	PTR			; Checkerboard
	.endif

*---------------------------------------------------------------------------------------
* Routine starts here
*---------------------------------------------------------------------------------------
	PUSH	ST					; save status register
	PUSH    DP
	LDP	@af0

	.if	CBUF=1
* Prepare addresses for Circular Buffer access.
* U   accesses are in mode   		  : 72,  1, 10,  5, 10, 1, 0, 0, 1, 1
* AF0, AF1, AF2, AF3 accesses are in mode : 12,  1, 12, 12, 12, 1, 0, 1, 1, 1
* All other accesses are in mode	  :  1,  1,  1,  1,  1, 0, 0, 1, 1, 0
* The sub-banks of AF's must alternate to start the next process. 
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
	LDI	@af0,		TMP
	ADDI	@bank2a,	TMP
	STI	TMP,		@af0_tmp
	LDI	@af1,		TMP
	ADDI	@bank2b,	TMP
	STI	TMP,		@af1_tmp
	LDI	@af2,		TMP
	ADDI	@bank2a,	TMP
	STI	TMP,		@af2_tmp
	LDI	@af3,		TMP
	ADDI	@bank2b,	TMP
	STI	TMP,		@af3_tmp
	.else
	LDI	@af0,		TMP
	STI	TMP,		@af0_tmp
	LDI	@af1,		TMP
	STI	TMP,		@af1_tmp
	LDI	@af2,		TMP
	STI	TMP,		@af2_tmp
	LDI	@af3,		TMP
	STI	TMP,		@af3_tmp
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
	STI	TMP,		@_wfm_mat_trick_count	; DMA spin count

*
* Initialize addresses for MAT_TRICK.
* U holds the address of the first element of   u
* TAS0 holds the address of the first element of  tas0 
* AF holds the address of the first element of  af 
*
* Note: DP same for all of cram and is set above
*
	STI	U,		@u_p			; not used yet
	LDI	@c_u0_p,	C_U
	LDI	@c_u1_p,	TMP
	XOR	C_U,		TMP			; toggle for C_U
	STI	TMP,		@c_u_t	
	LDI	6,		IR1			; stride in trick and mat
	LDF	0.0,		ZERO			; Constant
	LDI	@sign_p,	SIGN_P
	.if	ROUND=1
	RND	SIGN
	.endif
	STF	SIGN,		*SIGN_P++		; save  sign
        NEGF    SIGN   ,        TMP                     ; -sign
	.if	ROUND=1
	RND	TMP
	.endif
	STF	TMP,		*SIGN_P--		; save -sign

* Various addresses used for padding
	LSH	2,				PTR	; shift left by 2 and construct
	OR	1,				PTR	; the pointer into the PTR array
	MPYI	*+WILSON_AD(Wilson.offset),	PTR	; calculate and load the
	ADDI	*+WILSON_AD(Wilson.ptr),	PTR	; base address of the PTR array
	LDI	*+WILSON_AD(Wilson.yztmax),	LCYZT	; load the Y,Z,T loop counter
	LDI	@af_t,				MUPTR	; get base address of af0_t
							; The auxiliary register MUPTR
							; addresses four words of
							; memory that contain the base
							; addresses of AF0,AF1,AF2,AF3
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

	LDI	@af0_tmp,	AF0
	ADDI	*PTR++(1),	AF0,		TMP	; add displacement to AF0
	LDI	@af1_tmp,	AF1
	ADDI	*PTR++(1),	AF1,		TMP	; add displacement to AF1
||	STI	TMP,		*MUPTR++(1)		; store base address of AF0
	LDI	@af2_tmp,	AF2
	ADDI	*PTR++(1),	AF2,		TMP	; add displacement to AF2
||	STI	TMP,		*MUPTR++(1)		; store base address of AF1
	LDI	@af3_tmp,	AF3
	ADDI	*PTR++(1),	AF3,		TMP	; add displacement to AF3
||	STI	TMP,		*MUPTR++(1)		; store base address of AF2
	STI	TMP,		*MUPTR--(3)		; store base address of AF3
	LDI	*PTR++(1),	RC			; set the X loop counter

*
* Start X loop
*
	RPTB	LX					; loop over X
*
* Preload the AF spinors into cram pointed by TAS1
*
	LDI	@tas1_p,	TAS1			; source spinors for MAT

* Load address of destination  AF[0]
	LDI	@af0_t,		AF			; set the address
	LDF	*AF++,		TMP			; AF[0]

	LDF	*AF++,		TMP			; AF[1]
||	STF	TMP,		*TAS1++			; TAS1[0]
	LDF	*AF++,		TMP			; AF[2]
||	STF	TMP,		*TAS1++			; TAS1[1]
	LDF	*AF++,		TMP			; AF[3]
||	STF	TMP,		*TAS1++			; TAS1[2]
	LDF	*AF++,		TMP			; AF[4]
||	STF	TMP,		*TAS1++			; TAS1[3]
	LDF	*AF++,		TMP			; AF[5]
||	STF	TMP,		*TAS1++			; TAS1[4]
	LDF	*AF++,		TMP			; AF[6]
||	STF	TMP,		*TAS1++			; TAS1[5]
	LDF	*AF++,		TMP			; AF[7]
||	STF	TMP,		*TAS1++			; TAS1[6]
	LDF	*AF++,		TMP			; AF[8]
||	STF	TMP,		*TAS1++			; TAS1[7]
	LDF	*AF++,		TMP			; AF[9]
||	STF	TMP,		*TAS1++			; TAS1[8]
	LDF	*AF++,		TMP			; AF[10]
||	STF	TMP,		*TAS1++			; TAS1[9]
	LDF	*AF++,		TMP			; AF[11]
||	STF	TMP,		*TAS1++			; TAS1[10]

	STF	TMP,		*TAS1++			; TAS1[11]

* Store modified address of  AF[0]. AF points to the would be next site
	STI	AF,		@af0_t

* Load address of destination  AF[1]
	LDI	@af1_t,		AF			; set the address
	LDF	*AF++,		TMP			; AF[0]

	LDF	*AF++,		TMP			; AF[1]
||	STF	TMP,		*TAS1++			; TAS1[0]
	LDF	*AF++,		TMP			; AF[2]
||	STF	TMP,		*TAS1++			; TAS1[1]
	LDF	*AF++,		TMP			; AF[3]
||	STF	TMP,		*TAS1++			; TAS1[2]
	LDF	*AF++,		TMP			; AF[4]
||	STF	TMP,		*TAS1++			; TAS1[3]
	LDF	*AF++,		TMP			; AF[5]
||	STF	TMP,		*TAS1++			; TAS1[4]
	LDF	*AF++,		TMP			; AF[6]
||	STF	TMP,		*TAS1++			; TAS1[5]
	LDF	*AF++,		TMP			; AF[7]
||	STF	TMP,		*TAS1++			; TAS1[6]
	LDF	*AF++,		TMP			; AF[8]
||	STF	TMP,		*TAS1++			; TAS1[7]
	LDF	*AF++,		TMP			; AF[9]
||	STF	TMP,		*TAS1++			; TAS1[8]
	LDF	*AF++,		TMP			; AF[10]
||	STF	TMP,		*TAS1++			; TAS1[9]
	LDF	*AF++,		TMP			; AF[11]
||	STF	TMP,		*TAS1++			; TAS1[10]

	STF	TMP,		*TAS1++			; TAS1[11]

* Store modified address of  AF[1]. AF points to the would be next site
	STI	AF,		@af1_t

* Load address of destination  AF[2]
	LDI	@af2_t,		AF			; set the address
	LDF	*AF++,		TMP			; AF[0]

	LDF	*AF++,		TMP			; AF[1]
||	STF	TMP,		*TAS1++			; TAS1[0]
	LDF	*AF++,		TMP			; AF[2]
||	STF	TMP,		*TAS1++			; TAS1[1]
	LDF	*AF++,		TMP			; AF[3]
||	STF	TMP,		*TAS1++			; TAS1[2]
	LDF	*AF++,		TMP			; AF[4]
||	STF	TMP,		*TAS1++			; TAS1[3]
	LDF	*AF++,		TMP			; AF[5]
||	STF	TMP,		*TAS1++			; TAS1[4]
	LDF	*AF++,		TMP			; AF[6]
||	STF	TMP,		*TAS1++			; TAS1[5]
	LDF	*AF++,		TMP			; AF[7]
||	STF	TMP,		*TAS1++			; TAS1[6]
	LDF	*AF++,		TMP			; AF[8]
||	STF	TMP,		*TAS1++			; TAS1[7]
	LDF	*AF++,		TMP			; AF[9]
||	STF	TMP,		*TAS1++			; TAS1[8]
	LDF	*AF++,		TMP			; AF[10]
||	STF	TMP,		*TAS1++			; TAS1[9]
	LDF	*AF++,		TMP			; AF[11]
||	STF	TMP,		*TAS1++			; TAS1[10]

	STF	TMP,		*TAS1++			; TAS1[11]

* Store modified address of  AF[2]. AF points to the would be next site
	STI	AF,		@af2_t

* Load address of destination  AF[3]
	LDI	@af3_t,		AF			; set the address
	LDF	*AF++,		TMP			; AF[0]

	LDF	*AF++,		TMP			; AF[1]
||	STF	TMP,		*TAS1++			; TAS1[0]
	LDF	*AF++,		TMP			; AF[2]
||	STF	TMP,		*TAS1++			; TAS1[1]
	LDF	*AF++,		TMP			; AF[3]
||	STF	TMP,		*TAS1++			; TAS1[2]
	LDF	*AF++,		TMP			; AF[4]
||	STF	TMP,		*TAS1++			; TAS1[3]
	LDF	*AF++,		TMP			; AF[5]
||	STF	TMP,		*TAS1++			; TAS1[4]
	LDF	*AF++,		TMP			; AF[6]
||	STF	TMP,		*TAS1++			; TAS1[5]
	LDF	*AF++,		TMP			; AF[7]
||	STF	TMP,		*TAS1++			; TAS1[6]
	LDF	*AF++,		TMP			; AF[8]
||	STF	TMP,		*TAS1++			; TAS1[7]
	LDF	*AF++,		TMP			; AF[9]
||	STF	TMP,		*TAS1++			; TAS1[8]
	LDF	*AF++,		TMP			; AF[10]
||	STF	TMP,		*TAS1++			; TAS1[9]
	LDF	*AF++,		TMP			; AF[11]
||	STF	TMP,		*TAS1++			; TAS1[10]

	STF	TMP,		*TAS1--(AF_SIZE-1)	; reset TAS1 to beginning of block

* Store modified address of  AF[3]. AF points to the would be next site
	STI	AF,		@af3_t

*
* Gauge field successfully loaded.
* DMA load next U field. 
* NOTE: no error if read beyond end of field, but will not be used
*
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
* R2,R3	: accumulators for real/imaginary part
* R0	: temporary storage for multiply
* R6    : register holding floating zero
*
* Loop over 4 directions
* Multiply by TAS0 by C_U^dag for 2 spin components
*
	LDI	@tas0_p,	TAS0			; reload TAS0
	LDI	2,		IR0
	LDI	3,		LC

	MPYF3	*C_U++(IR1), 	*TAS1++(IR0),	TMP	; TMP  = Ur00 * Ar00

* Jump to MAT routine
	LDI	@cache_mat_p,	CACHE_MAT		; dram mat code
	LDI	@cont_p,	CODE_ADDR		; return address
	B	CACHE_MAT				; jump to dram

cont:
	.if	ROUND=1					; if round=1 cache holds less code
	MPYF3	*C_U--(IR1), 	*TAS1--(IR0),	TMP	; TMP  = Ui11 * Ar11
||	ADDF3	TMP,	   	R3,		R3	; R3 += TMP
	MPYF3	*C_U,	   	*TAS1++,	TMP	; TMP  = Ui10 * Ar11
||	ADDF3	TMP,	   	R3,		R3	; R3 += TMP

* Br11
	MPYF3	*C_U++(IR1), 	*TAS1++(IR0),	TMP	; TMP  = Ui10 * Ai01
||	ADDF3	TMP,	   	R3,		R3	; R3 += TMP
	MPYF3	*C_U++(IR1), 	*TAS1++(IR0),	TMP	; TMP  = Ui11 * Ai11
||	ADDF3	TMP, 	   	ZERO,		R2	; R2  = TMP
	MPYF3	*C_U--,      	*TAS1--,	TMP	; TMP  = Ui12 * Ai21
||	ADDF3	TMP,	   	R2,		R2	; R2 += TMP
	MPYF3	*C_U--(IR1), 	*TAS1--(IR0),	TMP	; TMP  = Ur12 * Ar21
||	ADDF3	TMP,	   	R2,		R2	; R2 += TMP
	MPYF3	*C_U--(IR1), 	*TAS1--(IR0),	TMP	; TMP  = Ur11 * Ar11
||	SUBF3	R2,	   	TMP,		R2	; R2  = TMP - R2
	MPYF3	*C_U--(IR0), 	*TAS1++,	TMP	; TMP  = Ur10 * Ar01
||	ADDF3	TMP,	   	R2,		R2	; R2 += TMP
	.endif

* Bi01
	MPYF3	*C_U++(IR1), 	*TAS1++(IR0),	TMP	; TMP  = Ur00 * Ai01
||	ADDF3	TMP,	   	R2,		R2	; R2 += TMP
	.if	ROUND=1
	RND	R2
	RND	R3
	.endif
	STF	R2,	   	*TAS0--(IR0)		; Br11 = R2
||	STF	R3,	   	*+TAS0			; Bi11 = R3
	MPYF3	*C_U++(IR1), 	*TAS1++(IR0),	TMP	; TMP  = Ur01 * Ai11
||	ADDF3	TMP, 	   	ZERO,		R3	; R3  = TMP
	MPYF3	*C_U++,      	*TAS1--,	TMP	; TMP  = Ur02 * Ai21
||	ADDF3	TMP,	   	R3,		R3	; R3 += TMP
	MPYF3	*C_U--(IR1), 	*TAS1--(IR0),	TMP	; TMP  = Ui02 * Ar21
||	ADDF3	TMP,	   	R3,		R3	; R3 += TMP
	MPYF3	*C_U--(IR1), 	*TAS1--(IR0),	TMP	; TMP  = Ui01 * Ar11
||	ADDF3	TMP,	   	R3,		R3	; R3 += TMP
	MPYF3	*C_U,	   	*TAS1++,	TMP	; TMP  = Ui00 * Ar01
||	ADDF3	TMP,	   	R3,		R3	; R3 += TMP

* Br01
	MPYF3	*C_U++(IR1), 	*TAS1++(IR0),	TMP	; TMP  = Ui00 * Ai01
||	ADDF3	TMP,	  	 R3,		R3	; R3 += TMP
	MPYF3	*C_U++(IR1), 	*TAS1++(IR0),	TMP	; TMP  = Ui01 * Ai11
||	ADDF3	TMP, 	   	ZERO,		R2	; R2  = TMP
	MPYF3	*C_U--,      	*TAS1--,	TMP	; TMP  = Ui02 * Ai21
||	ADDF3	TMP,	   	R2,		R2	; R2 += TMP
	MPYF3	*C_U--(IR1), 	*TAS1--(IR0),	TMP	; TMP  = Ur02 * Ar21
||	ADDF3	TMP,	   	R2,		R2	; R2 += TMP
	MPYF3	*C_U--(IR1), 	*TAS1--(IR0),	TMP	; TMP  = Ur01 * Ar11
||	SUBF3	R2,	   	TMP,		R2	; R2  = TMP - R2
	MPYF3	*C_U, 	   	*TAS1++(IR1),	TMP	; TMP  = Ur00 * Ar01
||	ADDF3	TMP,	   	R2,		R2	; R2 += TMP

	.if	ROUND=0
	DBD	LC,		CACHE_MAT		; branch back (if round=0)
	.endif

	NOP	*C_U++(18)				; jump to next dir
	MPYF3	*C_U++(IR1), 	*TAS1++(IR0),	TMP	; TMP  = Ur00 * Ar00
||	ADDF3	TMP,	   	R2,		R2	; R2 += TMP

	.if	ROUND=1
	DBD	LC,		CACHE_MAT		; branch back (if round=1)
	RND	R2
	RND	R3
	.endif

	STF	R2,	   	*TAS0++(IR1)		; Br01 = R2
||	STF	R3,	   	*+TAS0			; Bi01 = R3


;-------------------------------------------------------; end of mat loop
* Reset C_U to the beginning of the block. The reset is IR1+U_SIZE
* Reset TAS0 to the beginning of the block. The reset is AF_SIZE
* Reset TAS1 to the beginning of the block. The reset is IR0+AF_SIZE
	NOP	*C_U--(78)
	NOP	*TAS0--(AF_SIZE)

*
* Trick TAS0 into CHI
*
	LDI	@sign_p,	SIGN_P

	LDI	11,		IR0			; second stride

	LDI	2,		LC			; set loop counter
	LDI	TAS0,		AX
	ADDI	36,		AX			; AX = &TAS0(0,0,0,3)

trick:							;loop over colors
	ADDF3	*TAS0++(IR1),	*AX++(IR1),	CHIR0	; CHIr0  = Ar00 + Ar03
	ADDF3	*TAS0++(IR1),	*AX++,		CHIR1	; CHIr1  = Ar10 + Ar13
	ADDF	*TAS0++(IR1),	CHIR0			; CHIr0 += Ar01
	ADDF	*TAS0++(IR1),	CHIR1			; CHIr1 += Ar11
	ADDF	*TAS0++(IR1),	CHIR0			; CHIr0 += Ar02
	.if	ROUND=1
	RND	CHIR0
	.endif
	ADDF3	*TAS0++,	CHIR1,		CHIR1	; CHIr1 += Ar12
||	STF	CHIR0,		*CHI++(IR1)		; CHI(0,0) = CHIr0

	ADDF3	*TAS0--(IR1),	*AX--(IR1),	CHII1	; CHIi1  = Ai12 + Ai13
	ADDF3	*TAS0--(IR1),	*AX--,		CHII0	; CHIi0  = Ai02 + Ai03
	.if	ROUND=1
	RND	CHIR1
	.endif
	ADDF	*TAS0--(IR1),	CHII1,		CHII1	; CHIi1 += Ai11
||	STF	CHIR1,		*CHI++			; CHI(0,1) = CHIr1
	ADDF	*TAS0--(IR1),	CHII0			; CHIi0 += Ai01
	ADDF	*TAS0--(IR1),	CHII1			; CHIi1 += Ai10
	.if	ROUND=1
	RND	CHII1
	.endif
	ADDF3	*TAS0++(IR1),	CHII0,		CHII0	; CHIi0 += Ai00
||	STF	CHII1,		*CHI--(IR1)		; CHI(1,1) = CHIi1

	ADDF3	*TAS0--(IR1),	*AX++(IR1),	CHIR2	; CHIr2  = Ai10 + Ar03
	ADDF3	*TAS0--,	*AX++,		CHIR3	; CHIr3  = Ai00 + Ar13
	SUBF3	*AX--(IR1),	*TAS0++(IR1),	CHII3	; CHIi3  = Ar00 - Ai13
	SUBF3	*AX++,		*TAS0++(IR1),	CHII2	; CHIi2  = Ar10 - Ai03
	.if	ROUND=1
	RND	CHII0
	.endif
	SUBF3	CHIR3,		*TAS0,		CHIR3	; CHIr3  = Ar01 - CHIr3
||	STF	CHII0,		*CHI++(IR0)		; CHI(1,0) = CHIi0
	ADDF	*+TAS0(6),	CHIR2			; CHIr2 += Ar11
	ADDF	*+TAS0(13),	CHIR2			; CHIr2 += Ai02  &

	.if	ROUND=1
	RND	CHIR2
	.endif
	MPYF3	*+SIGN_P,	CHIR2,		CHIR2	; CHIR2  = CHIr2 * -sign
	ADDF	*+TAS0(19),	CHIR3			; CHIr3 += Ai12  &
	.if	ROUND=1
	RND	CHIR2
	RND	CHIR3
	.endif
	MPYF3	*SIGN_P,	CHIR3,		CHIR3	; CHIR3  = CHIr3 * sign
||	STF	CHIR2,		*CHI++(IR1)		; CHI(0,2) = CHIr2
	SUBF	*+TAS0(18),	CHII3			; CHIi3 -= Ar12
	ADDF	*+TAS0(1),	CHII3			; CHIi3 += Ai01  &
	.if	ROUND=1
	RND	CHIR3
	RND	CHII3
	.endif
	MPYF3	*SIGN_P,	CHII3,		CHII3	; CHIi3  = CHIi3 * sign
||	STF	CHIR3,		*CHI++			; CHI(0,3) = CHIr3
	ADDF	*+TAS0(12),	CHII2			; CHIi2 += Ar02
	SUBF	*+TAS0(7),	CHII2			; CHIi2 -= Ai11  &

	.if	ROUND=0	
	DBD	LC,		trick			; branch to trick (if round=0)
	.endif	

	.if	ROUND=1
	RND	CHII3
	RND	CHII2
	.endif
	MPYF3	*SIGN_P,	CHII2,		CHII2	; CHIi2  = CHIi2 * sign
||	STF	CHII3,		*CHI--(IR1)		; CHI(1,3) = CHIi3

	.if	ROUND=1
	DBD	LC,		trick			; branch to trick (if round=1)
	.endif	

	NOP	*TAS0--(10)
	.if	ROUND=1
	RND	CHII2
	.endif
	STF	CHII2,		*CHI--(IR0)		; CHI(1,2) = CHIi2

*-------------------------------------------------------; end of trick loop

	NOP	*TAS0--(6)				; point to beginning of block
	ADDI	18,		CHI			; CHI to next site

* Insure U field is loaded 
	LDI	@dmabase,	DMABASE			; DMABASE -> DMA control regs
	LDI	@_wfm_mat_trick_count,	COUNT
dmatest:
	ADDI	1,		COUNT
	LDI	*+DMABASE(dmatc), TMP			; inspect value of counter
	BNZ	dmatest					; spin until DMA is complete

	SUBI	1,		COUNT
	STI	COUNT,		@_wfm_mat_trick_count

* End of loops
LX:	NOP						; end of X-loop
	SUBI	1,		LCYZT
	BNN	LYZT					; end of Y, Z, T loop

*------------------------------------------------------------------------	[-
* Clear last DMA load
	LDI	@dmabase,	DMABASE			; DMABASE -> DMA control regs
	XOR	TMP,		TMP			; clear TMP
	STI	TMP, 		*+DMABASE(dmagc)	; halt any current DMA operation

	POP	DP
	POP	ST
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


;-------------------------------------------------------------------------
	.text

	.align

_cache_mat:

* Br00
	MPYF3	*C_U++(IR1),	*TAS1++(IR0),	TMP	; TMP  = Ur01 * Ar10
||	ADDF3	TMP, 	   	ZERO,		R2	; R2 += TMP
	MPYF3	*C_U++,      	*TAS1++,	TMP	; TMP  = Ur02 * Ar20
||	ADDF3	TMP, 	   	R2,		R2	; R2 += TMP
	MPYF3	*C_U--(IR1), 	*TAS1--(IR0),	TMP	; TMP  = Ui02 * Ai20
||	ADDF3	TMP,	   	R2,		R2	; R2 += TMP
	MPYF3	*C_U--(IR1), 	*TAS1--(IR0),	TMP	; TMP  = Ui01 * Ai10
||	SUBF3	TMP,	   	R2,		R2	; R2 -= TMP
	MPYF3	*C_U, 	   	*TAS1--,	TMP	; TMP  = Ui00 * Ai00
||	SUBF3	TMP,	   	R2,		R2	; R2 -= TMP

* Bi00
	MPYF3	*C_U++(IR1), 	*TAS1++(IR0),	TMP	; TMP  = Ui00 * Ar00
||	SUBF3	TMP,	   	R2,		R2	; R2 -= TMP
	MPYF3	*C_U++(IR1), 	*TAS1++(IR0),	TMP	; TMP  = Ui01 * Ar10
||	ADDF3	TMP, 	   	ZERO,		R3	; R3  = TMP
	MPYF3	*C_U--,      	*TAS1++,	TMP	; TMP  = Ui02 * Ar20
||	ADDF3	TMP,	   	R3,		R3	; R3 += TMP
	MPYF3	*C_U--(IR1), 	*TAS1--(IR0),	TMP	; TMP  = Ur02 * Ai20
||	ADDF3	TMP,	   	R3,		R3	; R3 += TMP
	MPYF3	*C_U--(IR1), 	*TAS1--(IR0),	TMP	; TMP  = Ur01 * Ai10
||	ADDF3	TMP,	  	R3,		R3	; R3 += TMP
	MPYF3	*C_U++(IR0), 	*TAS1--,	TMP	; TMP  = Ur00 * Ai00
||	ADDF3	TMP,	   	R3,		R3	; R3 += TMP

* Br10
	MPYF3	*C_U++(IR1), 	*TAS1++(IR0),	TMP	; TMP  = Ur10 * Ar00
||	ADDF3	TMP, 	   	R3,		R3	; R3 += TMP

	.if	ROUND=1
	RND	R2
	RND	R3
	.endif
	STF	R2,	   	*TAS0++(IR0)		; BTMP0 = R2
||	STF	R3,	   	*+TAS0			; Bi00 = R3

	MPYF3	*C_U++(IR1), 	*TAS1++(IR0),	TMP	; TMP  = Ur11 * Ar10
||	ADDF3	TMP, 	   	ZERO,		R2	; R2 += TMP
	MPYF3	*C_U++,      	*TAS1++,	TMP	; TMP  = Ur12 * Ar20
||	ADDF3	TMP, 	   	R2,		R2	; R2 += TMP
	MPYF3	*C_U--(IR1), 	*TAS1--(IR0),	TMP	; TMP  = Ui12 * Ai20
||	ADDF3	TMP,	   	R2,		R2	; R2 += TMP
	MPYF3	*C_U--(IR1), 	*TAS1--(IR0),	TMP	; TMP  = Ui11 * Ai10
||	SUBF3	TMP,	   	R2,		R2	; R2 -= TMP
	MPYF3	*C_U, 	   	*TAS1--,	TMP	; TMP  = Ui10 * Ai00
||	SUBF3	TMP,	   	R2,		R2	; R2 -= TMP

* Bi10
	MPYF3	*C_U++(IR1), 	*TAS1++(IR0),	TMP	; TMP  = Ui10 * Ar00
||	SUBF3	TMP,	   	R2,		R2	; R2 -= TMP
	MPYF3	*C_U++(IR1), 	*TAS1++(IR0),	TMP	; TMP  = Ui11 * Ar10
||	ADDF3	TMP, 	   	ZERO,		R3	; R3  = TMP
	MPYF3	*C_U--,      	*TAS1++,	TMP	; TMP  = Ui12 * Ar20
||	ADDF3	TMP,	   	R3,		R3	; R3 += TMP
	MPYF3	*C_U--(IR1), 	*TAS1--(IR0),	TMP	; TMP  = Ur12 * Ai20
||	ADDF3	TMP,	   	R3,		R3	; R3 += TMP
	MPYF3	*C_U--(IR1), 	*TAS1--(IR0),	TMP	; TMP  = Ur11 * Ai10
||	ADDF3	TMP,	   	R3,		R3	; R3 += TMP
	MPYF3	*C_U++(IR0), 	*TAS1--,	TMP	; TMP  = Ur10 * Ai00
||	ADDF3	TMP,	   	R3,		R3	; R3 += TMP

* Br20
	MPYF3	*C_U++(IR1), 	*TAS1++(IR0),	TMP	; TMP  = Ur20 * Ar00
||	ADDF3	TMP, 	   	R3,		R3	; R3 += TMP
	.if	ROUND=1
	RND	R2
	RND	R3
	.endif
	STF	R2,	   	*TAS0++(IR0)		; Br10 = R2
||	STF	R3,	   	*+TAS0			; Bi10 = R3
	MPYF3	*C_U++(IR1), 	*TAS1++(IR0),	TMP	; TMP  = Ur21 * Ar10
||	ADDF3	TMP, 	   	ZERO,		R2	; R2  = TMP
	MPYF3	*C_U++,      	*TAS1++,	TMP	; TMP  = Ur22 * Ar20
||	ADDF3	TMP, 	   	R2,		R2	; R2 += TMP
	MPYF3	*C_U--(IR1), 	*TAS1--(IR0),	TMP	; TMP  = Ui22 * Ai20
||	ADDF3	TMP,	   	R2,		R2	; R2 += TMP
	MPYF3	*C_U--(IR1), 	*TAS1--(IR0),	TMP	; TMP  = Ui21 * Ai10
||	SUBF3	TMP,	   	R2,		R2	; R2 -= TMP
	MPYF3	*C_U, 	   	*TAS1--,	TMP	; TMP  = Ui20 * Ai00
||	SUBF3	TMP,	   	R2,		R2	; R2 -= TMP

* Bi20
	MPYF3	*C_U++(IR1), 	*TAS1++(IR0),	TMP	; TMP  = Ui20 * Ar00
||	SUBF3	TMP,	   	R2,		R2	; R2 -= TMP
	MPYF3	*C_U++(IR1), 	*TAS1++(IR0),	TMP	; TMP  = Ui21 * Ar10
||	ADDF3	TMP, 	   	ZERO,		R3	; R3  = TMP
	MPYF3	*C_U--,      	*TAS1++,	TMP	; TMP  = Ui22 * Ar20
||	ADDF3	TMP,	   	R3,		R3	; R3 += TMP
	MPYF3	*C_U--(IR1), 	*TAS1--(IR0),	TMP	; TMP  = Ur22 * Ai20
||	ADDF3	TMP,	   	R3,		R3	; R3 += TMP
	MPYF3	*C_U--(IR1), 	*TAS1--(IR0),	TMP	; TMP  = Ur21 * Ai10
||	ADDF3	TMP,	   	R3,		R3	; R3 += TMP
	MPYF3	*C_U, 	   	*TAS1++(IR1),	TMP	; TMP  = Ur20 * Ai00
||	ADDF3	TMP,	   	R3,		R3	; R3 += TMP

* Reverse direction
* Bi21
	MPYF3	*C_U++(IR1), 	*TAS1++(IR0),	TMP	; TMP  = Ur20 * Ai01
||	ADDF3	TMP,	   	R3,		R3	; R3 += TMP
	.if	ROUND=1
	RND	R2
	RND	R3
	.endif
	STF	R2,	   	*TAS0++(IR1)		; Br20 = R2
||	STF	R3,	   	*+TAS0			; Bi20 = R3
	MPYF3	*C_U++(IR1), 	*TAS1++(IR0),	TMP	; TMP  = Ur21 * Ai11
||	ADDF3	TMP, 	   	ZERO,		R3	; R3  = TMP
	MPYF3	*C_U++,      	*TAS1--,	TMP	; TMP  = Ur22 * Ai21
||	ADDF3	TMP,	   	R3,		R3	; R3 += TMP
	MPYF3	*C_U--(IR1), 	*TAS1--(IR0),	TMP	; TMP  = Ui22 * Ar21
||	ADDF3	TMP,	   	R3,		R3	; R3 += TMP
	MPYF3	*C_U--(IR1), 	*TAS1--(IR0),	TMP	; TMP  = Ui21 * Ar11
||	ADDF3	TMP,	   	R3,		R3	; R3 += TMP
	MPYF3	*C_U,	   	*TAS1++,	TMP	; TMP  = Ui20 * Ar11
||	ADDF3	TMP,	   	R3,		R3	; R3 += TMP

* Br21
	MPYF3	*C_U++(IR1), 	*TAS1++(IR0),	TMP	; TMP  = Ui20 * Ai01
||	ADDF3	TMP,	   	R3,		R3	; R3 += TMP
	MPYF3	*C_U++(IR1), 	*TAS1++(IR0),	TMP	; TMP  = Ui21 * Ai11
||	ADDF3	TMP, 	   	ZERO,		R2	; R2  = TMP
	MPYF3	*C_U--,      	*TAS1--,	TMP	; TMP  = Ui22 * Ai21
||	ADDF3	TMP,	   	R2,		R2	; R2 += TMP
	MPYF3	*C_U--(IR1), 	*TAS1--(IR0),	TMP	; TMP  = Ur22 * Ar21
||	ADDF3	TMP,	   	R2,		R2	; R2 += TMP
	MPYF3	*C_U--(IR1), 	*TAS1--(IR0),	TMP	; TMP  = Ur21 * Ar11
||	SUBF3	R2,	   	TMP,		R2	; R2  = TMP - R2
	MPYF3	*C_U--(IR0), 	*TAS1++,	TMP	; TMP  = Ur20 * Ar01
||	ADDF3	TMP,	   	R2,		R2	; R2 += TMP

* Bi11
	MPYF3	*C_U++(IR1), 	*TAS1++(IR0),	TMP	; TMP  = Ur10 * Ai01
||	ADDF3	TMP,	   	R2,		R2	; R2 += TMP
	.if	ROUND=1
	RND	R2
	RND	R3
	.endif
	STF	R2,	   	*TAS0--(IR0)		; Br21 = R2
||	STF	R3,	   	*+TAS0			; Bi21 = R3

	.if	ROUND=1					; if round=1 cache holds less code
* Branch delayed back to cram
	BD	CODE_ADDR
	.endif

	MPYF3	*C_U++(IR1), 	*TAS1++(IR0),	TMP	; TMP  = Ur11 * Ai11
||	ADDF3	TMP, 	   	ZERO,		R3	; R3  = TMP
	MPYF3	*C_U++,      	*TAS1--,	TMP	; TMP  = Ur12 * Ai21
||	ADDF3	TMP,	   	R3,		R3	; R3 += TMP
	MPYF3	*C_U--(IR1), 	*TAS1--(IR0),	TMP	; TMP  = Ui12 * Ar21
||	ADDF3	TMP,	   	R3,		R3	; R3 += TMP

	.if	ROUND=0					; if round=0 cache holds more code
	MPYF3	*C_U--(IR1), 	*TAS1--(IR0),	TMP	; TMP  = Ui11 * Ar11
||	ADDF3	TMP,	   	R3,		R3	; R3 += TMP
	MPYF3	*C_U,	   	*TAS1++,	TMP	; TMP  = Ui10 * Ar11
||	ADDF3	TMP,	   	R3,		R3	; R3 += TMP

* Br11
	MPYF3	*C_U++(IR1), 	*TAS1++(IR0),	TMP	; TMP  = Ui10 * Ai01
||	ADDF3	TMP,	   	R3,		R3	; R3 += TMP
	MPYF3	*C_U++(IR1), 	*TAS1++(IR0),	TMP	; TMP  = Ui11 * Ai11
||	ADDF3	TMP, 	   	ZERO,		R2	; R2  = TMP
	MPYF3	*C_U--,      	*TAS1--,	TMP	; TMP  = Ui12 * Ai21
||	ADDF3	TMP,	   	R2,		R2	; R2 += TMP

* Branch delayed back to cram
	BD	CODE_ADDR

	MPYF3	*C_U--(IR1), 	*TAS1--(IR0),	TMP	; TMP  = Ur12 * Ar21
||	ADDF3	TMP,	   	R2,		R2	; R2 += TMP
	MPYF3	*C_U--(IR1), 	*TAS1--(IR0),	TMP	; TMP  = Ur11 * Ar11
||	SUBF3	R2,	   	TMP,		R2	; R2  = TMP - R2
	MPYF3	*C_U--(IR0), 	*TAS1++,	TMP	; TMP  = Ur10 * Ar01
||	ADDF3	TMP,	   	R2,		R2	; R2 += TMP
	.endif

	.end


