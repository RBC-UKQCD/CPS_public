**--------------------------------------------------------------------
**  CVS keywords
**
**  $Author: chulwoo $
**  $Date: 2004-06-04 21:14:05 $
**  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_clover/qcdsp/clover_mat_mlt_asm.asm,v 1.3 2004-06-04 21:14:05 chulwoo Exp $
**  $Id: clover_mat_mlt_asm.asm,v 1.3 2004-06-04 21:14:05 chulwoo Exp $
**  $Name: not supported by cvs2svn $
**  $Locker:  $
**  $Revision: 1.3 $
**  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_clover/qcdsp/clover_mat_mlt_asm.asm,v $
**  $State: Exp $
**
**--------------------------------------------------------------------

;=============================================================================
; clover_mat_mlt_asm.asm                             last modified on 11/18/98
;=============================================================================
; Implements void clover_mat_mlt_asm(float *out, const float *mat, 
;				     const float* in, int num_sites);
; Circular buffer:           2 and 4 for mat and vector respectively
; Chip ram:                  72 words for the clover matrix
; Registers:                 R0, R1, R2, R3, R4, AR0, AR1, AR2, AR3,
;                            AR4, DP, RS, RE, RC, BK(?)
; Stack frame            : full (frame pointer in AR3)                     
;                            *FP       old FP                                
;                            *-FP(1)   return address                        
;                            *-FP(2)   Y
;                            *-FP(3)   A
;                            *-FP(4)   X
;                            *-FP(5)   n
;=============================================================================
		.version 30

FP		.set            AR3
Y		.set		AR0
A		.set		AR1
X		.set		AR2
A_DRAM		.set		AR4
N_SITES		.set		R4
MEM_RC		.set		70		; 72 (MAT_SIZE) - 2
OFFSET		.set		71
CACHE_FREEZE	.set		400h
CACHE_ENABLE	.set		800h
CACHE_CLEAR	.set		1000h

		.ref	_clover_cram_scratch_addr        
		.def	_clover_mat_mlt_asm

		.text
_clover_mat_mlt_asm:
	PUSH	FP
	LDI	SP, FP

	PUSH	A_DRAM
	PUSH	N_SITES

; enable circular buffer  (use AR0, R0, R1 before messing up A,Y,X...)
	LDP	cbuf3_ctrl,	DP
	LDIU	@cbuf3_ctrl,	AR0
	LDIU	@cbuf2_mode,	R0
	LDIU	@cbuf4_mode,	R1
	STI	R0,		*--AR0
  ||	STI	R1,		*++AR0

; should not use AR0,AR1,AR2,AR3,AR4,R4 for other purpose from now on.
	LDIU	*-FP(2),	Y
	LDIU	*-FP(3),	A_DRAM
	LDIU	*-FP(4),	X
	LDIU	*-FP(5),	N_SITES

; use circular buffer
	LDP	BANK2_BASE,	DP
	LDIU	@BANK2_BASE,	R1
	ADDI	R1,		A_DRAM
	LDP	BANK4_BASE,	DP
	LDIU	@BANK4_BASE,	R0
	ADDI	R0,		X

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; DOES NOT SEEM TO MATTER AT ALL
; cache cleared
;	LDIU	CACHE_CLEAR,	R0
;	OR	R0,		ST
; cache enabled
;	LDIU	CACHE_ENABLE,	R1
;	OR	R1,		ST
; cache enabled and not frozen
;	LDIU	CACHE_FREEZE,	R0
;	ANDN	R0,		ST
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; IR0 = 5, A = CRAM_SCRATCH
	LDI	5,		IR0
	LDP	_clover_cram_scratch_addr,	DP	
	LDIU	@_clover_cram_scratch_addr,	A
	LDI	MEM_RC,		R2

Loop_over_sites:
; move A to cram
	LDI	*A_DRAM++,	R0
	RPTS	R2
	LDI	*A_DRAM++,	R0
    ||  STI	R0,		*A++
	STI	R0,		*A--(OFFSET)

; IR0 = 5, IR1 = 25
	LDI	1,		RC
	RPTB	clover_mat_mlt_end
	LDI     25,		IR1	

; R0,R1,R2,R3	are used for parallel instructions
; R2  : real part of Y     
; R3  : imaginary part of Y
; R0, R1  :   temporary storage for R2 and R3 respectively

; row 0
;    Y[0] = A[0]*X[0] + A[1]*X[2] + A[2]*X[3] + A[4]*X[4] + A[5]*X[5] 
;         + A[9]*X[6] + A[10]*X[7] + A[16]*X[8] + A[17]*X[9]
;         + A[25]*X[10] + A[26]*X[11]; 
;    Y[1] = A[0]*X[1] + A[1]*X[3] - A[2]*X[2] + A[4]*X[5] - A[5]*X[4] 
;         + A[9]*X[7] - A[10]*X[6] + A[16]*X[9] - A[17]*X[8] 
;         + A[25]*X[11] - A[26]*X[10];

	MPYF3	*A,	*X++,	R2	;R2 = A[0] * X[0]; 
	MPYF3	*A++,	*X++,	R3	;R3 = A[0] * X[1]; 
	MPYF3	*A,	*X,	R0	;R0 = A[1] * X[2];
	MPYF3	*A++,	*+X,	R1	;R1 = A[1] * X[3];
    ||	ADDF3	R0,	R2,	R2	;R2+= A[1] * X[2];
	MPYF3	*A,	*X++,	R1	;R1 = A[2] * X[2];
    ||	ADDF3	R1,	R3,	R3	;R3+= A[1] * X[3];
	MPYF3	*A++,	*X++,	R0	;R0 = A[2] * X[3]
    ||	SUBF3	R1,	R3,	R3	;R3-= A[2] * X[2];
	MPYF3	*++A,	*X,	R0	;R0 = A[4] * X[4];
    ||	ADDF3	R0,	R2,	R2	;R2+= A[2] * X[3]
	MPYF3	*A++,	*+X,	R1	;R1 = A[4] * X[5];
    ||	ADDF3	R0,	R2,	R2	;R2+= A[4] * X[4];
	MPYF3	*A,	*X++,	R1	;R1 = A[5] * X[4];
    ||	ADDF3	R1,	R3,	R3	;R3+= A[4] * X[5];
	MPYF3	*A--,	*X++,	R0	;R0 = A[5] * X[5];
    ||	SUBF3	R1,	R3,	R3	;R3-= A[5] * X[4];
	MPYF3	*++A(IR0),*X,	R0	;R0 = A[9] * X[6];
    ||	ADDF3	R0,	R2,	R2	;R2+= A[5] * X[5];
	MPYF3	*A++,	*+X,	R1	;R1 = A[9] * X[7]
    ||	ADDF3	R0,	R2,	R2	;R2+= A[9] * X[6]
	MPYF3	*A,	*X++,	R1	;R1 = A[10] * X[6]
    ||	ADDF3	R1,	R3,	R3	;R3+= A[9] * X[7]
	MPYF3	*A++,	*X++,	R0	;R0 = A[10] * X[7]
    ||	SUBF3	R1,	R3,	R3	;R3-= A[10] * X[6]
	MPYF3	*++A(IR0),*X,	R0	;R0 = A[16] * X[8]
    ||	ADDF3	R0,	R2,	R2	;R2+= A[10] * X[7]
	MPYF3	*+A,	*X++,	R1	;R1 = A[17] * X[8]
    ||	ADDF3	R0,	R2,	R2	;R2+= A[16] * X[8]
	MPYF3	*+A,	*X,	R0	;R0 = A[17] * X[9]
    ||	SUBF3	R1,	R3,	R3	;R3-= A[17] * X[8]
	MPYF3	*A++(IR0),*X++,	R1	;R1 = A[16] * X[9]
    ||	ADDF3	R0,	R2,	R2	;R2+= A[17] * X[9]
	MPYF3	*++A(IR0),*X,	R1	;R1 = A[26] * X[10]
    ||	ADDF3	R1,	R3,	R3	;R3+= A[16] * X[9]
	MPYF3	*--A,	*X++,	R0	;R0 = A[25] * X[10]
    ||	SUBF3	R1,	R3,	R3	;R3-= A[26] * X[10]
	MPYF3	*A++,	*X,	R1	;R1 = A[25] * X[11]
    ||	ADDF3	R0,	R2,	R2	;R2+= A[25] * X[10]
	MPYF3	*A--(IR1),*X--(IR0),R0	;R0 = A[26] * X[11]
    ||	ADDF3	R1,	R3,	R3	;R3+= A[25] * X[11]     done!
	RND	R3			;                       ready

; row 1
;    Y[2] = A[1]*X[0] - A[2]*X[1] + A[3]*X[2] + A[6]*X[4] + A[7]*X[5] 
;         + A[11]*X[6] + A[12]*X[7] + A[18]*X[8] + A[19]*X[9] 
;         + A[27]*X[10] + A[28]*X[11];
;    Y[3] = A[1]*X[1] + A[2]*X[0] + A[3]*X[3] + A[6]*X[5] - A[7]*X[4] 
;         + A[11]*X[7] - A[12]*X[6] + A[18]*X[9] - A[19]*X[8] 
;         + A[27]*X[11] - A[28]*X[10];

	MPYF3	*A,	*--X(IR0),R1	;R1 = A[1] * X[1]
    ||	ADDF3	R0,	R2,	R2	;R2+= A[26] * X[11]     done!
        RND	R2			;                       ready
	LDF	*-X,	R0		;R0 = X[0]
    ||  STF	R2,	*Y++		;                       Y[0] = real
	MPYF3	*A++,	R0,	R2	;R2 = A[1] * X[0]
    ||  STF	R3,	*Y++		;                       Y[1] = imag
	MPYF3	*A,	R0,	R3	;R3 = A[2] * X[0]
	MPYF3	*A++,	*X++,	R0	;R0 = A[2] * X[1]
    ||	ADDF3	R1,	R3,	R3	;R3 = A[1] * X[1] + A[2] * X[0]
	MPYF3	*A,	*X++,	R0	;R0 = A[3] * X[2]
    ||	SUBF3	R0,	R2,	R2	;R2 = A[1]*X[0] - A[2]*X[1]
	MPYF3	*A--,	*X++,	R1	;R1 = A[3] * X[3]
    ||	ADDF3	R0,	R2,	R2	;R2+= A[3] * X[2]
	MPYF3	*++A(IR0),*X,	R1	;R1 = A[7] * X[4]
    ||	ADDF3	R1,	R3,	R3	;R3+= A[3] * X[3]
	MPYF3	*--A,	*X++,	R0	;R0 = A[6] * X[4]
    ||	SUBF3	R1,	R3,	R3	;R3-= A[7] * X[4]
	MPYF3	*A++,	*X,	R1	;R1 = A[6] * X[5]
    ||	ADDF3	R0,	R2,	R2	;R2+= A[6] * X[4]
	MPYF3	*A++(IR0),*X++,	R0	;R0 = A[7] * X[5]
    ||	ADDF3	R1,	R3,	R3	;R3+= A[6] * X[5]
	MPYF3	*A--,	*X,	R1	;R1 = A[12] * X[6]
    ||	ADDF3	R0,	R2,	R2	;R2+= A[7] * X[5]
	MPYF3	*A,	*X++,	R0	;R0 = A[11] * X[6]
    ||	SUBF3	R1,	R3,	R3	;R3-= A[12] * X[6]
 	MPYF3	*A++,	*X,	R1	;R1 = A[11] * X[7]
    ||	ADDF3	R0,	R2,	R2	;R2+= A[11] * X[6]
	MPYF3	*A++(IR0),*X++,	R0	;R0 = A[12] * X[7]
    ||	ADDF3	R1,	R3,	R3	;R3+= A[11] * X[7]
	MPYF3	*++A,	*X,	R0	;R0 = A[18] * X[8]
    ||	ADDF3	R0,	R2,	R2	;R2+= A[12] * X[7]
 	MPYF3	*++A,	*X++,	R1	;R1 = A[19] * X[8]
    ||	ADDF3	R0,	R2,	R2	;R2+= A[18] * X[8]
	MPYF3	*A--,	*X,	R0	;R0 = A[19] * X[9]
    ||	SUBF3	R1,	R3,	R3	;R3-= A[19] * X[8]
	MPYF3	*A++(IR0),*X++,	R1	;R1 = A[18] * X[9]
    ||	ADDF3	R0,	R2,	R2	;R2+= A[19] * X[9]
	MPYF3	*++A(IR0),*X,	R1	;R1 = A[28] * X[10]
    ||	ADDF3	R1,	R3,	R3	;R3+= A[18] * X[9]
	MPYF3	*--A,	*X++,	R0	;R0 = A[27] * X[10]
    ||	SUBF3	R1,	R3,	R3	;R3-= A[28] * X[10]
	MPYF3	*A++,	*X,	R1	;R1 = A[27] * X[11]
    ||	ADDF3	R0,	R2,	R2	;R2+= A[27] * X[10]
	MPYF3	*A--(IR1),*X--(IR0),R0	;R0 = A[28] * X[11]
    ||	ADDF3	R1,	R3,	R3	;R3+= A[27] * X[11]     done!
	RND	R3			;                       ready
	
; row 2
;    Y[4] = A[4]*X[0] - A[5]*X[1] + A[6]*X[2] - A[7]*X[3] 
;         + A[8]*X[4] + A[13]*X[6] + A[14]*X[7] + A[20]*X[8] 
;         + A[21]*X[9] + A[29]*X[10] + A[30]*X[11];
;    Y[5] = A[4]*X[1] + A[5]*X[0] + A[6]*X[3] + A[7]*X[2] 
;         + A[8]*X[5] + A[13]*X[7] - A[14]*X[6] + A[20]*X[9] 
;         - A[21]*X[8] + A[29]*X[11] - A[30]*X[10];


	MPYF3	*++A,	*--X(IR0),R1	;R1 = A[4] * X[1]
    ||	ADDF3	R0,	R2,	R2	;R2+= A[28] * X[11]     done!
        RND	R2			;                       ready
	LDF	*-X,	R0		;R0 = X[0]
    ||  STF	R2,	*Y++		;                       Y[2] = real
	MPYF3	*A++,	R0,	R2	;R2 = A[4] * X[0]
    ||  STF	R3,	*Y++		;                       Y[3] = imag
	MPYF3	*A,	R0,	R3	;R3 = A[5] * X[0]
	MPYF3	*A++,	*X++,	R0	;R0 = A[5] * X[1]
    ||	ADDF3	R1,	R3,	R3	;R3 = A[4] * X[1] + A[5] * X[0]
	MPYF3	*A,	*X,	R0	;R0 = A[6] * X[2]
    ||	SUBF3	R0,	R2,	R2	;R2 = A[4] * X[0] - A[5] * X[1]
	MPYF3	*A++,	*+X,	R1	;R1 = A[6] * X[3]
    ||	ADDF3	R0,	R2,	R2	;R2+= A[6] * X[2]
	MPYF3	*A,	*X++,	R1	;R1 = A[7] * X[2]
    ||	ADDF3	R1,	R3,	R3	;R3+= A[6] * X[3]
	MPYF3	*A++,	*X++,	R0	;R0 = A[7] * X[3]
    ||	ADDF3	R1,	R3,	R3	;R3+= A[7] * X[2]
	MPYF3	*A,	*X++,	R0	;R0 = A[8] * X[4]
    ||	SUBF3	R0,	R2,	R2	;R2-= A[7] * X[3]
	MPYF3	*A++(IR0),*X++,	R1	;R1 = A[8] * X[5]
    ||	ADDF3	R0,	R2,	R2	;R2+= A[8] * X[4]
	MPYF3	*A,	*X,	R0	;R0 = A[13] * X[6]
    ||	ADDF3	R1,	R3,	R3	;R3+= A[8] * X[5]
	MPYF3	*A++,	*+X,	R1	;R1 = A[13] * X[7]
    ||	ADDF3	R0,	R2,	R2	;R2+= A[13] * X[6]
	MPYF3	*A,	*X++,	R1	;R1 = A[14] * X[6]
    ||	ADDF3	R1,	R3,	R3	;R3+= A[13] * X[7]
	MPYF3	*A++(IR0),*X++,	R0	;R0 = A[14] * X[7]
    ||	SUBF3	R1,	R3,	R3	;R3-= A[14] * X[6]
	MPYF3	*++A,	*X,	R0	;R0 = A[20] * X[8]
    ||	ADDF3	R0,	R2,	R2	;R2+= A[14] * X[7]
	MPYF3	*+A,	*X++,	R1	;R1 = A[21] * X[8]
    ||	ADDF3	R0,	R2,	R2	;R2+= A[20] * X[8]
	MPYF3	*+A,	*X,	R0	;R0 = A[21] * X[9]
    ||	SUBF3	R1,	R3,	R3	;R3-= A[21] * X[8]
	MPYF3	*A++(IR0),*X++,	R1	;R1 = A[20] * X[9]
    ||	ADDF3	R0,	R2,	R2	;R2+= A[21] * X[9]
	MPYF3	*++A(IR0),*X++,	R1	;R1 = A[30] * X[10]
    ||	ADDF3	R1,	R3,	R3	;R3+= A[20] * X[9]
	MPYF3	*A--,	*X,	R0	;R0 = A[30] * X[11]
    ||	SUBF3	R1,	R3,	R3	;R3-= A[30] * X[10]
	MPYF3	*A,	*X--,	R1	;R1 = A[29] * X[11]
    ||	ADDF3	R0,	R2,	R2	;R2+= A[30] * X[11]
	MPYF3	*A--(IR1),*X--(IR0),R0	;R0 = A[29] * X[10]
    ||	ADDF3	R1,	R3,	R3	;R3+= A[29] * X[11]     
	RND	R3			;                       done!

; row 3
;    Y[6] = A[9]*X[0] - A[10]*X[1] + A[11]*X[2] - A[12]*X[3] 
;         + A[13]*X[4] - A[14]*X[5] + A[15]*X[6] + A[22]*X[8] 
;         + A[23]*X[9] + A[31]*X[10] + A[32]*X[11];
;    Y[7] = A[9]*X[1] + A[10]*X[0] + A[11]*X[3] + A[12]*X[2] 
;         + A[13]*X[5] + A[14]*X[4] + A[15]*X[7] + A[22]*X[9] 
;         - A[23]*X[8] + A[31]*X[11] - A[32]*X[10];


	MPYF3	*++A(IR0),*--X(IR0),R0	;R0 = A[9] * X[0]
    ||	ADDF3	R0,	R2,	R2	;R2+= A[29] * X[10]     
        RND	R2			;                       done!
	LDF	*X++,	R1		;R1 = X[0]
    ||  STF	R2,	*Y++		;                       Y[4] = real
	MPYF3	*+A,	R1,	R3	;R3 = A[10] * X[0]
    ||  STF	R3,	*Y++		;                       Y[5] = imag
	MPYF3	*+A,	*X,	R2	;R2 = A[10] * X[1]
	MPYF3	*A++,	*X++,	R1	;R1 = A[9] * X[1]
    ||	SUBF3	R2,	R0,	R2	;R2 = A[9] * X[0] - A[10] * X[1]
; IR0 = 5, IR1 = 7
        LDI     7,	IR1
	MPYF3	*++A,	*X,	R0	;R0 = A[11] * X[2]
    ||	ADDF3	R1,	R3,	R3	;R3 = A[10] * X[0] + A[9] * X[1]
	MPYF3	*A++,	*+X,	R1	;R1 = A[11] * X[3]
    ||	ADDF3	R0,	R2,	R2	;R2+= A[11] * X[2]
	MPYF3	*A,	*X++,	R1	;R1 = A[12] * X[2]
    ||	ADDF3	R1,	R3,	R3	;R3+= A[11] * X[3]
	MPYF3	*A++,	*X++,	R0	;R0 = A[12] * X[3]
    ||	ADDF3	R1,	R3,	R3	;R3+= A[12] * X[2]
	MPYF3	*A,	*X,	R0	;R0 = A[13] * X[4]
    ||	SUBF3	R0,	R2,	R2	;R2-= A[12] * X[3]
	MPYF3	*A++,	*+X,	R1	;R1 = A[13] * X[5]
    ||	ADDF3	R0,	R2,	R2	;R2+= A[13] * X[4]
	MPYF3	*A,	*X++,	R1	;R1 = A[14] * X[4]
    ||	ADDF3	R1,	R3,	R3	;R3+= A[13] * X[5]
	MPYF3	*A++,	*X++,	R0	;R0 = A[14] * X[5]
    ||	ADDF3	R1,	R3,	R3	;R3+= A[14] * X[4]
	MPYF3	*A,	*X++,	R0	;R0 = A[15] * X[6]
    ||	SUBF3	R0,	R2,	R2	;R2-= A[14] * X[5]
	MPYF3	*A++(IR1),*X++,	R1	;R1 = A[15] * X[7]
    ||	ADDF3	R0,	R2,	R2	;R2+= A[15] * X[6]
	MPYF3	*A,	*X,	R0	;R0 = A[22] * X[8]
    ||	ADDF3	R1,	R3,	R3	;R3+= A[15] * X[7]
	MPYF3	*A++,	*+X,	R1	;R1 = A[22] * X[9]
    ||	ADDF3	R0,	R2,	R2	;R2+= A[22] * X[8]
	MPYF3	*A,	*X++,	R1	;R1 = A[23] * X[8]
    ||	ADDF3	R1,	R3,	R3	;R3+= A[22] * X[9]
	MPYF3	*A++(IR1),*X++,	R0	;R0 = A[23] * X[9]
    ||	SUBF3	R1,	R3,	R3	;R3-= A[23] * X[8]
; IR0 = 5, IR1 = 10
        LDI     10,	IR1
	MPYF3	*++A,	*X,	R0	;R0 = A[31] * X[10]
    ||	ADDF3	R0,	R2,	R2	;R2+= A[23] * X[9]
	MPYF3	*A++,	*+X,	R1	;R1 = A[31] * X[11]
    ||	ADDF3	R0,	R2,	R2	;R2+= A[31] * X[10]
	MPYF3	*A,	*X++,	R1	;R1 = A[32] * X[10]
    ||	ADDF3	R1,	R3,	R3	;R3+= A[31] * X[11]
	MPYF3	*A--(IR1),*X--(IR1),R0	;R0 = A[32] * X[11]
    ||	SUBF3	R1,	R3,	R3	;R3-= A[32] * X[10]     done!
	RND	R3			;                       ready

; row 4
;    Y[8] = A[16]*X[0] - A[17]*X[1] + A[18]*X[2] - A[19]*X[3] 
;         + A[20]*X[4] - A[21]*X[5] + A[22]*X[6] - A[23]*X[7] 
;         + A[24]*X[8] + A[33]*X[10] + A[34]*X[11];
;    Y[9] = A[16]*X[1] + A[17]*X[0] + A[18]*X[3] + A[19]*X[2] 
;         + A[20]*X[5] + A[21]*X[4] + A[22]*X[7] + A[23]*X[6] 
;         + A[24]*X[9] + A[33]*X[11] - A[34]*X[10];

	MPYF3	*--A(IR0),*X,	R0	;R0 = A[17] * X[1]
    ||	ADDF3	R0,	R2,	R2	;R2+= A[32] * X[11]     done!
        RND	R2			;                       ready
	LDF	*-X,	R1		;R1 = X[0]
    ||  STF	R2,	*Y++		;                       Y[6] = real
	MPYF3	*A--,	R1,	R3	;R3 = A[17] * X[0]
    ||  STF	R3,	*Y++		;                       Y[7] = imag
	MPYF3	*A,	R1,	R2	;R2 = A[16] * X[0]
	MPYF3	*A++,	*X++,	R1	;R1 = A[16] * X[1]
    ||	SUBF3	R0,	R2,	R2	;R2 = A[16] * X[0] - A[17] * X[1]
	MPYF3	*++A,	*X,	R0	;R0 = A[18] * X[2]
    ||	ADDF3	R1,	R3,	R3	;R3 = A[16] * X[1] + A[17] * X[0]
	MPYF3	*A++,	*+X,	R1	;R1 = A[18] * X[3]
    ||	ADDF3	R0,	R2,	R2	;R2+= A[18] * X[2]
	MPYF3	*A,	*X++,	R1	;R1 = A[19] * X[2]
    ||	ADDF3	R1,	R3,	R3	;R3+= A[18] * X[3]
	MPYF3	*A++,	*X++,	R0	;R0 = A[19] * X[3]
    ||	ADDF3	R1,	R3,	R3	;R3+= A[19] * X[2]
	MPYF3	*A,	*X,	R0	;R0 = A[20] * X[4]
    ||	SUBF3	R0,	R2,	R2	;R2-= A[19] * X[3]
	MPYF3	*A++,	*+X,	R1	;R1 = A[20] * X[5]
    ||	ADDF3	R0,	R2,	R2	;R2+= A[20] * X[4]
	MPYF3	*A,	*X++,	R1	;R1 = A[21] * X[4]
    ||	ADDF3	R1,	R3,	R3	;R3+= A[20] * X[5]
	MPYF3	*A++,	*X++,	R0	;R0 = A[21] * X[5]
    ||	ADDF3	R1,	R3,	R3	;R3+= A[21] * X[4]
	MPYF3	*A,	*X,	R0	;R0 = A[22] * X[6]
    ||	SUBF3	R0,	R2,	R2	;R2-= A[21] * X[5]
	MPYF3	*A++,	*+X,	R1	;R1 = A[22] * X[7]
    ||	ADDF3	R0,	R2,	R2	;R2+= A[22] * X[6]
	MPYF3	*A,	*X++,	R1	;R1 = A[23] * X[6]
    ||	ADDF3	R1,	R3,	R3	;R3+= A[22] * X[7]
	MPYF3	*A++,	*X++,	R0	;R0 = A[23] * X[7]
    ||	ADDF3	R1,	R3,	R3	;R3+= A[23] * X[6]
	MPYF3	*A,	*X++,	R0	;R0 = A[24] * X[8]
    ||	SUBF3	R0,	R2,	R2	;R2-= A[23] * X[7]
	MPYF3	*A++(IR1),*X++,	R1	;R1 = A[24] * X[9]
    ||	ADDF3	R0,	R2,	R2	;R2+= A[24] * X[8]
	MPYF3	*--A,	*X,	R0	;R0 = A[33] * X[10]
    ||	ADDF3	R1,	R3,	R3	;R3+= A[23] * X[6]
	MPYF3	*A++,	*+X,	R1	;R1 = A[33] * X[11]
    ||	ADDF3	R0,	R2,	R2	;R2+= A[33] * X[10]
	MPYF3	*A,	*X++,	R1	;R1 = A[34] * X[10]
    ||	ADDF3	R1,	R3,	R3	;R3+= A[33] * X[11]
	MPYF3	*A--(IR1),*X--(IR1),R0	;R0 = A[34] * X[11]
    ||	SUBF3	R1,	R3,	R3	;R3-= A[34] * X[10]
	RND	R3			;			done!

; row 5
;    Y[10] = A[25]*X[0] - A[26]*X[1] + A[27]*X[2] - A[28]*X[3] 
;          + A[29]*X[4] - A[30]*X[5] + A[31]*X[6] - A[32]*X[7] 
;          + A[33]*X[8] - A[34]*X[9] + A[35]*X[10];
;    Y[11] = A[25]*X[1] + A[26]*X[0] + A[27]*X[3] + A[28]*X[2] 
;          + A[29]*X[5] + A[30]*X[4] + A[31]*X[7] + A[32]*X[6] 
;          + A[33]*X[9] + A[34]*X[8] + A[35]*X[11];

	MPYF3	*++A,	*X,	R1	;R1 = A[25] * X[1]
    ||	ADDF3	R0,	R2,	R2	;R2+= A[32] * X[11]     done!
        RND	R2			;                       ready

	LDF	*-X,	R0		;R0 = X[0]
    ||  STF	R2,	*Y++		;                       Y[8] = real
	MPYF3	*+A,	R0,	R3	;R3 = A[26] * X[0]
    ||  STF	R3,	*Y++		;                       Y[9] = imag
	MPYF3	*A++,	R0,	R2	;R2 = A[25] * X[0]
	MPYF3	*A++,	*X++,	R0	;R0 = A[26] * X[1]
    ||	ADDF3	R1,	R3,	R3	;R3 = A[25] * X[1] + A[26] * X[0]
	MPYF3	*A,	*X,	R0	;R0 = A[27] * X[2]
    ||	SUBF3	R0,	R2,	R2	;R2 = A[25] * X[0] - A[26] * X[1]
	MPYF3	*A++,	*+X,	R1	;R1 = A[27] * X[3]
    ||	ADDF3	R0,	R2,	R2	;R2+= A[27] * X[2]
	MPYF3	*A,	*X++,	R1	;R1 = A[28] * X[2]
    ||	ADDF3	R1,	R3,	R3	;R3+= A[27] * X[3]
	MPYF3	*A++,	*X++,	R0	;R0 = A[28] * X[3]
    ||	ADDF3	R1,	R3,	R3	;R3+= A[28] * X[2]
	MPYF3	*A,	*X,	R0	;R0 = A[29] * X[4]
    ||	SUBF3	R0,	R2,	R2	;R2-= A[28] * X[3]
	MPYF3	*A++,	*+X,	R1	;R1 = A[29] * X[5]
    ||	ADDF3	R0,	R2,	R2	;R2+= A[29] * X[4]
	MPYF3	*A,	*X++,	R1	;R1 = A[30] * X[4]
    ||	ADDF3	R1,	R3,	R3	;R3+= A[29] * X[5]
	MPYF3	*A++,	*X++,	R0	;R0 = A[30] * X[5]
    ||	ADDF3	R1,	R3,	R3	;R3+= A[30] * X[4]
	MPYF3	*A,	*X,	R0	;R0 = A[31] * X[6]
    ||	SUBF3	R0,	R2,	R2	;R2-= A[30] * X[5]
	MPYF3	*A++,	*+X,	R1	;R1 = A[31] * X[7]
    ||	ADDF3	R0,	R2,	R2	;R2+= A[30] * X[6]
	MPYF3	*A,	*X++,	R1	;R1 = A[32] * X[6]
    ||	ADDF3	R1,	R3,	R3	;R3+= A[31] * X[7]
	MPYF3	*A++,	*X++,	R0	;R0 = A[32] * X[7]
    ||	ADDF3	R1,	R3,	R3	;R3+= A[32] * X[6]
	MPYF3	*A,	*X,	R0	;R0 = A[33] * X[8]
    ||	SUBF3	R0,	R2,	R2	;R2-= A[32] * X[7]
	MPYF3	*A++,	*+X,	R1	;R1 = A[33] * X[9]
    ||	ADDF3	R0,	R2,	R2	;R2+= A[33] * X[8]
	MPYF3	*A,	*X++,	R1	;R1 = A[34] * X[8]
    ||	ADDF3	R1,	R3,	R3	;R3+= A[33] * X[9]
	MPYF3	*A++,	*X++,	R0	;R0 = A[34] * X[9]
    ||	ADDF3	R1,	R3,	R3	;R3+= A[34] * X[8]
	MPYF3	*A,	*X++,	R0	;R0 = A[35] * X[10]
    ||	SUBF3	R0,	R2,	R2	;R2-= A[34] * X[9]
	MPYF3	*A++,	*X++,	R1	;R1 = A[35] * X[11]
    ||	ADDF3	R0,	R2,	R2	;R2+= A[35] * X[10]
	ADDF3	R1,	R3,	R3	;R3+= A[35] * X[11] 
	RND	R2			;			done!
        RND	R3			;                       ready
	STF	R2,	*Y++		;                       Y[10] = real
clover_mat_mlt_end:
	STF	R3,	*Y++		;                       Y[11] = imag

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; DOES NOT SEEM TO MATTER AT ALL
; freeze cache
;	LDIU	CACHE_FREEZE,	R0
;	OR	R0,		ST
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	SUBI	1,	N_SITES
	BNZD	Loop_over_sites
	LDP	_clover_cram_scratch_addr,	DP	
	LDIU	@_clover_cram_scratch_addr,	A
	LDI	MEM_RC,		R2

	POP	N_SITES
	POP	A_DRAM
	POP	FP
	RETSU

		.data
cbuf3_ctrl	.word		815803h
cbuf2_mode	.word		0cca52148h
cbuf4_mode	.word		0ce318118h
BANK2_BASE	.word		0a00000h
BANK4_BASE	.word		0c00000h

	.end


