**--------------------------------------------------------------------
**  CVS keywords
**
**  $Author: zs $
**  $Date: 2004-08-18 11:57:38 $
**  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/alg/alg_ghb/qcdsp/cmhb_kern.asm,v 1.4 2004-08-18 11:57:38 zs Exp $
**  $Id: cmhb_kern.asm,v 1.4 2004-08-18 11:57:38 zs Exp $
**  $Name: not supported by cvs2svn $
**  $Locker:  $
**  $Revision: 1.4 $
**  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/alg/alg_ghb/qcdsp/cmhb_kern.asm,v $
**  $State: Exp $
**
**--------------------------------------------------------------------
*****************************************************************
*  You must be using assign_cram.asm if you are using this	*
* routine.  Do NOT attempt to use this routine with the		*
* necessary variables defined in C, it will not work!		*
*****************************************************************

	.asg	1,	ROUND

*  Make the name of this routine global
	.def	_cmhb_kernel__FP6rfloatT1 ; publish run time address

*	.def	src_update	; publish first load time address
*	.def	end_update	; publish last load time address

*  Some symbols which might be useful
	.ref	_sqrt
	.ref	_cos
	.ref	_bad_sqrt
	.ref	_log_cram
	.ref	_grand
	.ref	_m_multiply3
	.ref	_m_equal
	.ref	_m_identity
	.ref	_exit
	.ref	_m_print
	.ref	_iprint
	.ref	_fprint
	.ref	_fprint2
	.ref	_xprint

*  External variables
	.ref	_el_seed_p
	.ref	_core_iscratch
	.ref	_core_fscratch

	.asg	AR3,	FP
	.asg	AR5,	seed_for_make_rand
	.asg	AR6,	iscratch
	.asg	AR7,	fscratch

	.asg	*+AR0(0),	bman0
	.asg	*+AR0(1),	bman1
	.asg	*+AR0(IR0),	bman2
	.asg	*+AR0(IR1),	bman3
	.asg	*+AR1(0),	rinv0
	.asg	*+AR1(1),	rinv1
	.asg	*+AR1(IR0),	rinv2
	.asg	*+AR1(IR1),	rinv3
	.asg	*+AR5(0),	twoxtwo_alpha0
	.asg	*+AR5(1),	twoxtwo_alpha1
	.asg	*+AR5(IR0),	twoxtwo_alpha2
	.asg	*+AR5(IR1),	twoxtwo_alpha3

* define the offsets in scratch which locate particular
* integer variables. All such offsets are 
* designated by the "I_" prefix.
	.asg	0,	I_subblock
	.asg	1,	I_link_p
	.asg	2,	I_sigma_p

* define the offsets in scratch which locate particular
* floating point variables. All such offsets are 
* designated by the "F_" prefix.
	.asg	0,	F_u_sigma00r
	.asg	1,	F_u_sigma01r
	.asg	2,	F_u_sigma02r
	.asg	3,	F_u_sigma10r
	.asg	4,	F_u_sigma11r
	.asg	5,	F_u_sigma12r
	.asg	6,	F_u_sigma20r
	.asg	7,	F_u_sigma21r
	.asg	8,	F_u_sigma22r
	.asg	9,	F_u_sigma00i
	.asg	10,	F_u_sigma01i
	.asg	11,	F_u_sigma02i
	.asg	12,	F_u_sigma10i
	.asg	13,	F_u_sigma11i
	.asg	14,	F_u_sigma12i
	.asg	15,	F_u_sigma20i
	.asg	16,	F_u_sigma21i
	.asg	17,	F_u_sigma22i
	.asg	18,	F_r0
	.asg	19,	F_r1
	.asg	20,	F_r2
	.asg	21,	F_r3
	.asg	22,	F_k
	.asg	23,	F_lambda
	.asg	24,	F_lambda_inv
	.asg	25,	F_lambda_inv_square
	.asg	26,	F_lambda_exp_inv
	.asg	27,	F_lambda_distribution
	.asg	28,	F_z
	.asg	29,	F_log_z
	.asg	30,	F_bman0
	.asg	31,	F_bman1
	.asg	32,	F_bman2
	.asg	33,	F_bman3
	.asg	34,	F_twoxtwo_alpha0
	.asg	35,	F_twoxtwo_alpha1
	.asg	36,	F_twoxtwo_alpha2
	.asg	37,	F_twoxtwo_alpha3
	.asg	38,	F_temp00r
	.asg	56,	F_temp22i
	.asg	57,	F_gamma
	.asg	58,	F_scale

*	.sect	".update"
	.text
*	.label	src_update		;	label load time address
_cmhb_kernel__FP6rfloatT1
*  Push the old value of AR3 onto the stack
	PUSH	FP
*  Load the Stack Pointer into FP
	LDI	SP, FP
*  Save all registers that are important to C
	PUSH	DP

        PUSH    R4
        PUSH    R5
        PUSHF   R6
        PUSHF   R7
        PUSH    AR4
        PUSH    AR5
        PUSH    AR6
        PUSH    AR7
*        PUSH    FP              ; Local frame pointer
*	PUSH	DP

*	LDI	99, R0		;	DEBUG
*	PUSH	R0		;	DEBUG
*	CALL	_iprint		;	DEBUG
*	SUBI	1, SP		;	DEBUG
_initialization
	LDP	@_core_iscratch
	LDI	@_core_iscratch, iscratch
*	LDP	@_core_fscratch
	LDI	@_core_fscratch, fscratch
	
	LDI	0, R0
	STI	R0, *+iscratch(I_subblock)

	LDI	*-FP(3), R0	; link
	STI	R0, *+iscratch(I_link_p)

	LDI	*-FP(2), R0	; sigma
	STI	R0, *+iscratch(I_sigma_p)
	
_make_usigma
*****************************************************
*  In this section the u-sigma matrix is created
* by multiplying the two arguments to cmhb_kernel()
* togeather.  This is done by an explicit call to
* m_multiply3 instead  of an in-lined code segment
* both to save memory in cram, and because of the
* assumption that m_multiply3 is in the cache.
*****************************************************
*	CALL	_mus
*	LDI	*+iscratch(I_subblock), R0
*	CMPI	0, R0
*	BNZ	_extract_subblock_one
*	BU	_extract_subblock_zero
*****************************************************

*	PUSH	FP
	LDI	*+iscratch(I_sigma_p), R0
	PUSH	R0
	LDI	*+iscratch(I_link_p), R0
	PUSH	R0
	PUSH	fscratch
	CALL	_m_multiply3
	SUBI	3, SP
*	POP	FP

	LDI	4, IR0		;  Load these two here since they
	LDI	3, IR1		; are trampled bu _m_multiply3

	LDI	fscratch, AR0
	LDI	fscratch, AR4
	ADDI	F_u_sigma00r, AR0
	ADDI	F_r0, AR4

	LDI	*+iscratch(I_subblock), R0
*	CMPI	0, R0
	BNZD	_extract_subblock_one
	LDI	AR0, AR1
	ADDI	9, AR1
	LDF	0.5, R5


_extract_subblock_zero
*****************************************************
*  Here the 4-vector form of the SU2 subblock located
* in the (0,0) corner of the matrix u-sigma is 
* extracted, it is not a SU2 object yet, that will
* only be true after some rescaling in the
* _process_subblock section.
*****************************************************
*	CALL	_es0
*	BU	_process_subblock
*****************************************************

	ADDF3	*+AR0(IR0), *+AR0(0), R2
	ADDF3	*+AR1(IR1), *+AR1(1), R3
||	MPYF3	R5, R2, R0
	SUBF3	*+AR0(IR1), *+AR0(1), R2
||	MPYF3	R5, R3, R1
	SUBF3	*+AR1(IR0), *+AR1(0), R3
	MPYF	R5, R2

	.if	ROUND
	MPYF	R5, R3
	RND	R0
	RND	R1
	RND	R2
	RND	R3
	BD	_process_subblock
	LDI	2, IR0
	.else
	MPYF	R5, R3
	BD	_process_subblock
	LDI	2, IR0
	.endif
	STF	R0, *+AR4(0)
||	STF	R1, *+AR4(1)
	STF	R2, *+AR4(IR0)
||	STF	R3, *+AR4(IR1)
	
_extract_subblock_one
*****************************************************
*  Here the 4-vector form of the SU2 subblock located
* in the (0,0) corner of the matrix u-sigma is 
* extracted, it is not a SU2 object yet, that will
* only be true after some rescaling in the
* _process_subblock section.
*****************************************************
*	CALL	_es1
*	BU	_process_subblock
*****************************************************
	ADDI	4, AR0
	ADDI	4, AR1

	ADDF3	*+AR0(IR0), *+AR0(0), R2
	ADDF3	*+AR1(IR1), *+AR1(1), R3
||	MPYF3	R5, R2, R0
	SUBF3	*+AR0(IR1), *+AR0(1), R2
||	MPYF3	R5, R3, R1
	SUBF3	*+AR1(IR0), *+AR1(0), R3
	MPYF	R5, R2

	.if	ROUND
	MPYF	R5, R3
	RND	R0
	RND	R1
	RND	R2
	RND	R3
	BD	_process_subblock
	LDI	2, IR0
	.else
	MPYF	R5, R3
	LDI	2, IR0
	.endif

	STF	R0, *+AR4(0)
||	STF	R1, *+AR4(1)
	STF	R2, *+AR4(IR0)
||	STF	R3, *+AR4(IR1)

_process_subblock
*****************************************************
*  Here the twoxtwo subblock extracted from u-sigma is
* rescaled into a SU2 object.  Thereafter it is
* used to generate the boltzman distributed matrix
* which will be used for the update.
*****************************************************
*	CALL	_ps
*	LDI	*+iscratch(I_subblock), R0
*	BNZ	_make_alpha1
*	BU	_make_alpha0
*****************************************************
*  calculate the normalization constant "k"

	MPYF3	*+AR4(0), *+AR4(0), R2		; *+fscratch(F_r0)^2
	MPYF3	*+AR4(1), *+AR4(1), R3		; *+fscratch(F_r1)^2
	MPYF3	*+AR4(IR0), *+AR4(IR0), R0	; *+fscratch(F_r2)^2
||	ADDF3	R2, R3, R2
	MPYF3	*+AR4(IR1), *+AR4(IR1), R1	; *+fscratch(F_r3)^2
||	ADDF3	R2, R0, R2
	ADDF3	R2, R1, R2
	PUSHF	R2
	CALL	_sqrt_cram
	SUBI	1,SP

	.if	ROUND
	RND	R0
	.endif
	
	STF	R0, *+fscratch(F_k)	

	LDF	R0, R1
	LDF	1., R0
	CALL	DIV_F30_cram_c

	LDI	AR4, AR5
	MPYF3	*AR4++, R0, R2
	.if	ROUND
	RND	R2
	.endif
	MPYF3	*AR4++, R0, R2
||	STF	R2, *AR5++
	.if	ROUND
	RND	R2
	.endif
	MPYF3	*AR4++, R0, R2
||	STF	R2, *AR5++
	.if	ROUND
	RND	R2
	.endif
	MPYF3	*AR4++, R0, R2
||	STF	R2, *AR5++
	.if	ROUND
	RND	R2
	.endif
	STF	R2, *AR5

***********************************************
* Begin here the Kennedy Pendleton algorithm. *
***********************************************
kp:
	PUSH	DP
	LDP	@_el_seed_p
	LDI	@_el_seed_p, seed_for_make_rand
	POP	DP

_kp_generate_a_zero

*	PUSH	DP	;	DEBUGGING
*	LDF	*+fscratch(F_k),R0	;	DEBUGGING
*	PUSHF	R0	;	DEBUGGING
*	CALL	_fprint	;	DEBUGGING
*	SUBI	1, SP	;	DEBUGGING
*	POP	DP	;	DEBUGGING

*	LDP	@_make_rand	;	DEBUG
*	LDI	@_make_rand, R0	;	DEBUG
*	PUSH	R0		;	DEBUG
*	CALL	_iprint		;	DEBUG
*	SUBI	1, SP		;	DEBUG
*	LDP	@_el_seed_p	;	DEBUG

	CALL	_grand
	LDF	R0, R7			;  R7 = R''
	ABSF	R7
	CALL	_grand
	ABSF	R0
	PUSHF	R0
	CALL	_log_cram
	SUBI	1, SP
	LDF	R0, R6			;  R6 = log(R')
	CALL	_grand
	ABSF	R0
	PUSHF	R0
	CALL	_log_cram
	SUBI	1, SP
	LDF	R0, R5			;  R5 = log(R)
					;  el_seed_p on stack

	LDF	*+fscratch(F_gamma), R1 ;  R1 = 2 beta / 3
	MPYF	*+fscratch(F_k), R1	;  alpha = (2/3) beta k

	LDF	1., R0
	CALL	DIV_F30_cram_c		;  R0 = 1/alpha
	MPYF	R0, R6			;  R6 = -X' = log(R')/alpha
	MPYF	R0, R5			;  R5 = -X  = log(R )/alpha
*
	PUSH	DP
	LDP	@UC_PI
	MPYF	@UC_PI, R7
	POP	DP

	MPYF	2., R7
	PUSHF	R7
	CALL	_cos_cram		;  R0 = cos( 2*PI*R'')
	SUBI	1, SP
	MPYF3	R0, R0, R1		;  R1 = C = cos( 2*PI*R'')^2

	MPYF	R5, R1			;  R1 = -A = (-X)C
	ADDF	R6, R1			;  R1 = -d = -A - X'

	LDF	R1, R7			;  R7 = -d
	CALL	_grand			;  R0 = R'''
	ABSF	R0
	MPYF	R0, R0			;  R0 = (R''')^2
	MPYF	0.5, R7
	SUBF	R7, R0
	SUBF	1., R0			;  R0 = (R''')^2 + (d/2) - 1
*****
*	PUSH	DP			;	DEBUGGING
*	PUSHF	R7			;	DEBUGGING
*	PUSHF	R0			;	DEBUGGING
*	CALL	_fprint			;	DEBUGGING
*	POPF	R0			;	DEBUGGING
*	POPF	R7			;	DEBUGGING
*	POP	DP			;	DEBUGGING
*	LDF	R0, R1			;	DEBUGGING
*****

	BNN	_kp_generate_a_zero	;  Branch if R''' is too big.
	MPYF	2., R7			;  R7 = -d
	ADDF	1., R7			;  R7 = 1 - d
	.if	ROUND
	RND	R7
	.endif
	STF	R7, *+fscratch(F_bman0)	
*****
*	PUSH	DP			;	DEBUGGING
*	PUSHF	R7			;	DEBUGGING
*	CALL	_fprint			;	DEBUGGING
*	POPF	R7			;	DEBUGGING
*	POP	DP			;	DEBUGGING
*****
	
*********************************************
* End here the Kennedy Pendleton algorithm. *
*********************************************

*  Use the trigonometric algorithm to generate points on 
* the unit sphere:
trig:	MPYF3	R7, R7, R1	;	R1 = b0^2
	SUBRF	1., R1		;	R1 = 1-b0^2
	PUSHF	R1
	CALL	_sqrt_cram	;	R0 = sqrt( 1-b0^2 )
	SUBI	1, SP
	.if	ROUND
	RND	R0
	.endif
	STF	R0, *+fscratch(F_scale)

*	PUSHF	R0		;	DEBUG
*	CALL	_fprint		;	DEBUG
*	POPF	R0		;	DEBUG

	CALL	_grand
	LDF	R0, R6		;	R6 = cos(theta) e [-1,1)

	CALL	_grand
	LDF	R0, R7
	MPYF	@UC_PI, R7		;	R7 = phi e [-pi,pi)

* Now R6 = cos(theta), R7 = phi

*	;	Print Phi
*	PUSH	DP
*	PUSH	FP
*	PUSHF	R6		;	DEBUG
*	PUSHF	R7		;	DEBUG
*	CALL	_fprint		;	DEBUG
*	POPF	R7		;	DEBUG
*	POPF	R6		;	DEBUG
*	POP	FP
*	POP	DP

*	;	Print cos(Theta)
*	PUSH	DP
*	PUSH	FP
*	PUSHF	R7		;	DEBUG
*	PUSHF	R6		;	DEBUG
*	CALL	_fprint		;	DEBUG
*	POPF	R6		;	DEBUG
*	POPF	R7		;	DEBUG
*	POP	FP
*	POP	DP

	;	Calculate bman3 and store it
	LDF	R6, R0
	MPYF	*+fscratch(F_scale), R0
	.if	ROUND
	RND	R0
	.endif
	STF	R0, *+fscratch(F_bman3)
	
	;	cos(Theta) is no longer needed, replace it
	;	with sin(Theta).
	MPYF3	R6, R6, R0
	SUBRF	1., R0
	ABSF	R0	
	PUSHF	R0
	CALL	_sqrt_cram
	SUBI	1, SP
	LDF	R0, R6		;	R6 = sin(Theta)

*	;	Print sin(Theta)
*	PUSHF	R7		;	DEBUG
*	PUSHF	R6		;	DEBUG
*	CALL	_fprint		;	DEBUG
*	POPF	R6		;	DEBUG
*	POPF	R7		;	DEBUG

	;	Calculate bman1 and store it.
	PUSHF	R7
	CALL	_cos_cram
	SUBI	1, SP		;	R0 = cos(Phi)
	MPYF3	R0, R6, R1	;	R1 = sin(Theta)*cos(Phi)
	MPYF	*+fscratch(F_scale), R1
	.if	ROUND
	RND	R1
	.endif
	STF	R1, *+fscratch(F_bman1)

*	;	Print cos(Phi)
*	PUSHF	R7		;	DEBUG
*	PUSHF	R6		;	DEBUG
*	PUSHF	R0		;	DEBUG
*	CALL	_fprint		;	DEBUG
*	POPF	R0		;	DEBUG
*	POPF	R6		;	DEBUG
*	POPF	R7		;	DEBUG

	;	cos(Phi) is no longer needed, replace it
	;	with sin(Phi)
	MPYF3	R0, R0, R2
	SUBRF	1., R2
	ABSF	R2	
	PUSHF	R2
	CALL	_sqrt_cram
	SUBI	1, SP		;	R0 = sin(Phi)
	CMPF	0., R7		;	Compare phi to 0 
	BNN	noneg
	NEGF	R0		;	If (Phi < 0) sin(Phi) *= -1;
noneg:

*	;	Print sin(Phi)
*	PUSHF	R7		;	DEBUG
*	PUSHF	R6		;	DEBUG
*	PUSHF	R0		;	DEBUG
*	CALL	_fprint		;	DEBUG
*	POPF	R0		;	DEBUG
*	POPF	R6		;	DEBUG
*	POPF	R7		;	DEBUG

	;	Calculate bman2 and store it.
	MPYF3	R6, R0, R1	;	R1 = sin(Theta)*sin(Phi)
	MPYF	*+fscratch(F_scale), R1
	.if	ROUND
	RND	R1
	.endif
	STF	R1, *+fscratch(F_bman2)

*	LDF	*+fscratch(F_bman0), R2;DEBUG
*	PUSHF	R2		;	DEBUG
*	CALL	_fprint		;	DEBUG
*	SUBI	1, SP		;	DEBUG
*	LDF	*+fscratch(F_bman1), R2;DEBUG
*	PUSHF	R2		;	DEBUG
*	CALL	_fprint		;	DEBUG
*	SUBI	1, SP		;	DEBUG
*	LDF	*+fscratch(F_bman2), R2;DEBUG
*	PUSHF	R2		;	DEBUG
*	CALL	_fprint		;	DEBUG
*	SUBI	1, SP		;	DEBUG
*	LDF	*+fscratch(F_bman3), R2;DEBUG
*	PUSHF	R2		;	DEBUG
*	CALL	_fprint		;	DEBUG
*	SUBI	1, SP		;	DEBUG

_make_2x2_alpha
*****************************************************
*  Now that B-> is in scratch the next thing to do is to calculate
* the inverse of the two by two matrix r->. Since the determinant
* is one, the inverse is just the transpose of the cofactor matrix.
*****************************************************

	LDI	F_r0, R4
	LDI	fscratch, AR0
	ADDI3	R4, AR0, AR1

	LDF	*+AR1(0), R0 
	STF	R0, *+AR1(0) 	
||	NEGF	*+AR1(1), R1 
	STF	R1, *+AR1(1) 	
||	NEGF	*+AR1(IR0), R2 
	STF	R2, *+AR1(IR0) 	
||	NEGF	*+AR1(IR1), R3 
	STF	R3, *+AR1(IR1)	

*  The matrix alpha is specified by the twoxtwo subblock which
* is given by the product of B and r_inv, call this subblock "a".
*
*	Set up things for the multiplication so that the
*	maximal number of parallel instructions can be used.
*	rinv0             = AR1 = scratch+F_r0
*	bman0             = AR0 = scratch+F_bman
*	twoxtwo_alpha0    = AR5 = fscratch+F_twoxtwo_alpha

	ADDI	F_bman0, AR0
	LDI	fscratch, AR5
	ADDI	F_twoxtwo_alpha0, AR5

*	begin the multiplication
	MPYF3	bman0, rinv0, R0
	MPYF3	bman3, rinv3, R1
	MPYF3	bman2, rinv2, R0
||	SUBF3	R1, R0, R2
	MPYF3	bman1, rinv1, R0
||	SUBF3	R0, R2, R2

	MPYF3	bman3, rinv2, R0
||	SUBF3	R0, R2, R2
	MPYF3	bman0, rinv1, R1
	.if	ROUND
	RND	R2
	.endif
	STF	R2, twoxtwo_alpha0
	MPYF3	bman2, rinv3, R0
||	ADDF3	R0, R1, R2
	MPYF3	bman1, rinv0, R0
||	SUBF3	R0, R2, R2

	MPYF3	bman0, rinv2, R0
||	ADDF3	R0, R2, R2
	MPYF3	bman3, rinv1, R1
	.if	ROUND
	RND	R2
	.endif
	STF	R2, twoxtwo_alpha1
	MPYF3	bman2, rinv0, R0
||	SUBF3	R1, R0, R2
	MPYF3	bman1, rinv3, R0
||	ADDF3	R0, R2, R2

	MPYF3	bman0, rinv3, R0
||	ADDF3	R0, R2, R2
	MPYF3	bman3, rinv0, R1
	.if	ROUND
	RND	R2
	.endif
	STF	R2, twoxtwo_alpha2
	MPYF3	bman1, rinv2, R0
||	ADDF3	R1, R0, R2
	MPYF3	bman2, rinv1, R0
||	SUBF3	R0, R2, R2

	ADDF3	R0, R2, R2
	.if	ROUND
	RND	R2
	.endif
	STF	R2, twoxtwo_alpha3

*	LDF	*+fscratch(F_twoxtwo_alpha0), R0;	DEBUGGING
*	PUSHF	R0				;	DEBUGGING
*	CALL	_fprint2			;	DEBUGGING
*	POPF	R0				;	DEBUGGING
*	LDF	*+fscratch(F_twoxtwo_alpha1), R0;	DEBUGGING
*	PUSHF	R0				;	DEBUGGING
*	CALL	_fprint2			;	DEBUGGING
*	POPF	R0				;	DEBUGGING
*	LDF	*+fscratch(F_twoxtwo_alpha2), R0;	DEBUGGING
*	PUSHF	R0				;	DEBUGGING
*	CALL	_fprint2			;	DEBUGGING
*	POPF	R0				;	DEBUGGING
*	LDF	*+fscratch(F_twoxtwo_alpha3), R0;	DEBUGGING
*	PUSHF	R0				;	DEBUGGING
*	CALL	_fprint2			;	DEBUGGING
*	POPF	R0				;	DEBUGGING

*  Now the tho by two form of alpha (the update matrix)
* has been constructed and is is time to update the
* link matrix by multiplying it by a block diagonal
* matrix containing alpha (2x2) and 1 (1x1) as it's 'blocks'.

*  Set up for the multiplication:
*	AR0 = &link (Re)
*	AR1 = &link (Im)
*	AR2 = &temp (Re)
*	AR4 = &temp (Im)
	LDI	*+iscratch(I_link_p), AR0
	LDI	AR0, AR1
	ADDI	9, AR1
	LDI	fscratch, AR2

*  If this is subblock 0, set up some registers and jump
* over some additional stuff to get started.
	LDI	*+iscratch(I_subblock), R0
	BZD	r_update0

	ADDI	F_temp00r, AR2
	LDI	AR2, AR4
	ADDI	9, AR4

r_update1
*  If subblock 1 is to be updated, then some of the pointers
* need to be incremented further.
	ADDI	3, AR0
	ADDI	3, AR1
	ADDI	3, AR2
	ADDI	3, AR4


r_update0
*  All comments are apropriate to subblock 0, if this is
* the second pass through, all of the comments are not 
* exactly correct.
* Calculate t00r
	MPYF3	twoxtwo_alpha0, *AR0, R0	;  R0 = a00r*u00r
	MPYF3	twoxtwo_alpha2, *++AR0(IR1), R1	;  R1 = a01r*u10r
	MPYF3	twoxtwo_alpha3, *AR1, R0	;  R0 = a00i*u00i
||	ADDF3	R0, R1, R2			;  R2 = a00r*u00r + a01r*u10r
	MPYF3	twoxtwo_alpha1, *++AR1(IR1), R0	;  R0 = a01i*u10i
||	SUBF3	R0, R2, R2			;  R2 = a00r*u00r + a01r*u10r
						;     - a00i*u00i

* Calculate t01r/t02r
	SUBI	1, AR2
	LDI	1, RC
	RPTB	block00
	MPYF3	twoxtwo_alpha0, *--AR0(IR0), R0	;  R0 = a00r*u01r
||	SUBF3	R0, R2, R2			;  R2 = a00r*u00r + a01r*u10r
						;     - a00i*u00i - a01i*u10i
	MPYF3	twoxtwo_alpha2, *++AR0(IR1), R1	;  R1 = a01r*u11r
	.if	ROUND
	RND	R2
	.endif
	STF	R2, *++AR2(1)			;  t00r = a00r*u00r + a01r*u10r
	MPYF3	twoxtwo_alpha3, *--AR1(IR0), R0	;  R0 = a00i*u01i
||	ADDF3	R0, R1, R2			;  R2 = a00r*u01r + a01r*u11r
block00
	MPYF3	twoxtwo_alpha1, *++AR1(IR1), R0 ;  R0 = a01i*u11i
||	SUBF3	R0, R2, R2			;  R2 = a00r*u01r + a01r*u11r
						;     - a00i*u01i

* Calculate t10r
	SUBI	5, AR0
	SUBI	5, AR1
	MPYF3	twoxtwo_alpha2, *AR0, R0	;  R0 =-a10r*u00r
||	SUBF3	R0, R2, R2			;  R2 = a00r*u02r + a01r*u12r
						;     - a00i*u02i - a01i*u12i
	MPYF3	twoxtwo_alpha0, *++AR0(IR1), R1	;  R1 = a11r*u10r
	.if	ROUND
	RND	R2
	.endif
	STF	R2, *++AR2(1)			;  t02r = a00r*u02r + a01r*u12r
						;       - a00i*u02i - a01i*u12i
	MPYF3	twoxtwo_alpha1, *AR1, R0	;  R0 = a10i*u00i
||	SUBF3	R0, R1, R2			;  R2 = a10r*u00r + a11r*u10r
	MPYF3	twoxtwo_alpha3, *++AR1(IR1), R0 ;  R0 =-a11i*u10i
||	SUBF3	R0, R2, R2			;  R2 = a10r*u00r + a11r*u10r
						;     - a10i*u00i
	
* Calculate t11r/t12r
	LDI	1, RC
	RPTB	block01
	MPYF3	twoxtwo_alpha2, *--AR0(IR0), R0	;  R0 =-a10r*u01r
||	ADDF3	R0, R2, R2			;  R2 = a10r*u00r + a11r*u10r
						;     - a10i*u00i - a11i*u10i
	MPYF3	twoxtwo_alpha0, *++AR0(IR1), R1	;  R1 = a11r*u11r
	.if	ROUND
	RND	R2
	.endif
	STF	R2, *++AR2(1)			;  t10r = a10r*u00r + a11r*u10r
						;       - a10i*u00i - a11i*u10i
	MPYF3	twoxtwo_alpha1, *--AR1(IR0), R0	;  R0 = a10i*u01i
||	SUBF3	R0, R1, R2			;  R2 = a10r*u01r + a11r*u11r
block01
	MPYF3	twoxtwo_alpha3, *++AR1(IR1), R0 ;  R0 =-a11i*u11i
||	SUBF3	R0, R2, R2			;  R2 = a10r*u01r + a11r*u11r
						;     - a10i*u01i
	
	ADDF3	R0, R2, R2			;  R2 = a10r*u01r + a11r*u11r
						;     - a10i*u01i - a11i*u11i
	.if	ROUND
	RND	R2
	.endif
	STF	R2, *++AR2(1)			;  t12r = a10r*u01r + a11r*u11r
						;     - a10i*u01i - a11i*u11i
	ADDI	3, AR2

	LDP	@_core_iscratch
	LDI	@_core_iscratch, iscratch
	LDI	*+iscratch(I_subblock), R0
	BZD	i_update0
	LDI	*+iscratch(I_link_p), AR0
	LDI	AR0, AR1
	ADDI	9, AR1

i_update1
	ADDI	3, AR0
	ADDI	3, AR1

i_update0
* Calculate t00i
	MPYF3	twoxtwo_alpha3, *AR0, R0	;  R0 = a00i*u00r
	MPYF3	twoxtwo_alpha1, *++AR0(IR1), R1	;  R1 = a01i*u10r
	MPYF3	twoxtwo_alpha0, *AR1, R0	;  R0 = a00r*u00i
||	ADDF3	R0, R1, R2			;  R2 = a00i*u00r + a01i*u10r
	MPYF3	twoxtwo_alpha2, *++AR1(IR1), R0	;  R0 = a01r*u10i
||	ADDF3	R0, R2, R2			;  R2 = a00i*u00r + a01i*u10r
						;     + a00r*u00i

* Calculate t01i
	MPYF3	twoxtwo_alpha3, *--AR0(IR0), R0	;  R0 = a00i*u01r
||	ADDF3	R0, R2, R2			;  R2 = a00i*u00r + a01i*u10r
						;     + a00r*u00i + a01r*u10i
	MPYF3	twoxtwo_alpha1, *++AR0(IR1), R1	;  R1 = a01i*u11r
	.if	ROUND
	RND	R2
	.endif
	STF	R2, *++AR2(1)			;  t00i = a00i*u00r + a01i*u10r
						;       + a00r*u00i + a01r*u10i
	MPYF3	twoxtwo_alpha0, *--AR1(IR0), R0	;  R0 = a00r*u01i
||	ADDF3	R0, R1, R2			;  R2 = a00i*u01r + a01i*u11r
	MPYF3	twoxtwo_alpha2, *++AR1(IR1), R0 ;  R0 = a01r*u11i
||	ADDF3	R0, R2, R2			;  R2 = a00i*u01r + a01i*u11r
						;     + a00r*u01i

* Calculate t02i
	MPYF3	twoxtwo_alpha3, *--AR0(IR0), R0	;  R0 = a00i*u02r
||	ADDF3	R0, R2, R2			;  R2 = a00i*u01r + a01i*u11r
						;     + a00r*u01i + a01r*u11i
	MPYF3	twoxtwo_alpha1, *++AR0(IR1), R1	;  R1 = a01i*u12r
	.if	ROUND
	RND	R2
	.endif
	STF	R2, *++AR2(1)			;  t01i = a00i*u01r + a01i*u11r
						;       + a00r*u01i + a01r*u11i
	MPYF3	twoxtwo_alpha0, *--AR1(IR0), R0	;  R0 = a00r*u02i
||	ADDF3	R0, R1, R2			;  R2 = a00i*u02r + a01i*u12r
	MPYF3	twoxtwo_alpha2, *++AR1(IR1), R0 ;  R0 = a01r*u12i
||	ADDF3	R0, R2, R2			;  R2 = a00i*u02r + a01i*u12r
						;     + a00r*u02i

* Calculate t10i
	SUBI	5, AR0
	SUBI	5, AR1
	MPYF3	twoxtwo_alpha1, *AR0, R0	;  R0 = a10i*u00r
||	ADDF3	R0, R2, R2			;  R2 = a00i*u02r + a01i*u12r
						;     + a00r*u02i + a01r*u12i
	MPYF3	twoxtwo_alpha3, *++AR0(IR1), R1	;  R1 =-a11i*u10r
	.if	ROUND
	RND	R2
	.endif
	STF	R2, *++AR2(1)			;  t02i = a00i*u02r + a01i*u12r
						;       + a00r*u02i + a01r*u12i
	MPYF3	twoxtwo_alpha2, *AR1, R0	;  R0 =-a10r*u00i
||	SUBF3	R1, R0, R2			;  R2 = a10i*u00r + a11i*u10r
	MPYF3	twoxtwo_alpha0, *++AR1(IR1), R0 ;  R0 = a11r*u10i
||	SUBF3	R0, R2, R2			;  R2 = a10i*u00r + a11i*u10r
						;     + a10r*u00i
	
* Calculate t11i
	MPYF3	twoxtwo_alpha1, *--AR0(IR0), R0	;  R0 = a10i*u01r
||	ADDF3	R0, R2, R2			;  R2 = a10i*u00r + a11i*u10r
						;     + a10r*u00i + a11r*u10i
	MPYF3	twoxtwo_alpha3, *++AR0(IR1), R1	;  R1 =-a11i*u11r
	.if	ROUND
	RND	R2
	.endif
	STF	R2, *++AR2(1)			;  t11i = a10i*u00r + a11i*u10r
						;       + a10r*u00i + a11r*u10i
	MPYF3	twoxtwo_alpha2, *--AR1(IR0), R0	;  R0 =-a10r*u01i
||	SUBF3	R1, R0, R2			;  R2 = a10i*u01r + a11i*u11r
	MPYF3	twoxtwo_alpha0, *++AR1(IR1), R0 ;  R0 = a11r*u11i
||	SUBF3	R0, R2, R2			;  R2 = a10i*u01r + a11i*u11r
						;     + a10r*u01i
	
* Calculate t12i
	MPYF3	twoxtwo_alpha1, *--AR0(IR0), R0	;  R0 = a10i*u02r
||	ADDF3	R0, R2, R2			;  R2 = a10i*u01r + a11i*u11r
						;     + a10r*u01i + a11r*u11i
	MPYF3	twoxtwo_alpha3, *++AR0(IR1), R1	;  R1 =-a11i*u12r
	.if	ROUND
	RND	R2
	.endif
	STF	R2, *++AR2(1)			;  t12i = a10i*u01r + a11i*u11r
						;       + a10r*u01i + a11r*u11i
	MPYF3	twoxtwo_alpha2, *--AR1(IR0), R0	;  R0 =-a10r*u02i
||	SUBF3	R1, R0, R2			;  R2 = a10i*u02r + a11i*u12r
	MPYF3	twoxtwo_alpha0, *++AR1(IR1), R0 ;  R0 = a11r*u12i
||	SUBF3	R0, R2, R2			;  R2 = a10i*u02r + a11i*u12r
						;     + a10r*u02i
	ADDF3	R0, R2, R2			;  R2 = a10i*u01r + a11i*u11r
						;     + a10r*u01i + a11r*u11i
	.if	ROUND
	RND	R2
	.endif
	STF	R2, *++AR2(1)			;  t20i = a10i*u01r + a11i*u11r
						;       + a10r*u01i + a11r*u11i
*  Copy this matrix from temp, writing over the origional U 
* with the updated U.
	LDI	@_core_iscratch, iscratch
	LDI	*+iscratch(I_subblock), R0
	BZD	copy0
	LDI	*+iscratch(I_link_p), AR0
	LDI	@_core_fscratch, AR1
	ADDI	F_temp00r, AR1
copy1
	ADDI	3, AR0
	ADDI	3, AR1

copy0
	LDF	*AR1++, R0

	RPTS	5
	LDF	*AR1++, R0
||	STF	R0, *AR0++
	ADDI	2, AR1
	ADDI	3, AR0
	LDF	*AR1++, R0
	RPTS	5
	LDF	*AR1++, R0
||	STF	R0, *AR0++

*	LDI	*+iscratch(I_link_p), R0;	DEBUGGING
*	PUSH	R0			;	DEBUGGING
*	CALL	_m_print		;	DEBUGGING
*	POP	R0			;	DEBUGGING

*  Increment the subblock counter, and branch back to do subblock1 if
* the counter's current value is 0.  Otherwise, just continue on and
* close up.
	LDI	*+iscratch(I_subblock), R0
	BZD	_make_usigma
	ADDI	1, R0
	STI	R0, *+iscratch(I_subblock)
	NOP

_clean_up_and_leave
*	CALL	_cu
*  Restore the registers before returning to the C program 
*
*        POP     DP
*        POP     FP              ;  This is our frame pointer
        POP     AR7
        POP     AR6
        POP     AR5
        POP     AR4
        POPF    R7
        POPF    R6
        POP     R5
        POP     R4

	POP	DP
        POP     FP              ;  This is our caller's frame pointer

        RETS

UC_PI	.float	3.14159265358979323846

*******************************************************************
	.def	faststack_p_c
faststack_p_c	.word	faststack_c	; label run time address
faststack_c	.space	27 ; 27
*******************************************************************

_sqrt_cram:
	LDI     SP,AR0
	LDF     *-AR0(1),R2   ; x
 
	BGT     pos0           
	LDFLE   0.0,R0        ; if x == 0, return 0
	RETSEQ

	.ref	_errno        ;*Reference this myself.
	LDP     @_errno       ; in case big model
	LDI     1 ,R1         ; indicate error
	STI     R1,@_errno
	CALL	_bad_sqrt     ;*This is my own error case.
*	RETS

pos0:    LDF     R2,R3        ; save x
	MPYF	2.0,R2        ; add a rounding bit in exponent
	PUSHF	R2            ; push x as float
	POP	R1            ; pop as int
	ASH	-25,R1        ; e = exponent(x) / 2
;
; determine initial estimate .25 * 2**(-e/2)
;
	NEGI	R1            ; negate exponent 
	ASH	24,R1         ; shift into place
	PUSH	R1            ; push as int 
	POPF	R1            ; pop as float
	MPYF	0.25,R2       ; remove rounding bit
;
; iterate 5 times
;
	LDI	4, RC
	RPTB	sqloop0
	RND     R1
	MPYF	R1,R2,R0      ; R0 = x[4] * (v/2)
	MPYF	R1,R0         ; R0 = (v/2) * x[4] * x[4]
	SUBRF	1.5,R0	      ; R0 = 1.5 - (v/2) * x[4] * x[4]
sqloop0:	MPYF	R0,R1	      ; x[5] = x[4] * (1.5 - v/2 * x[4] * x[4])
;
	RND     R1,R2
	MPYF    R3,R2,R0      ; sqrt(x) = x * sqrt(1/x)

	RETS


******************************************************

	.def	DIV_F30_cram_c
DIV_F30_cram_c:  
        POP     AR1         ; Pop return address
	PUSH    R2          ; Save R2: integer part
	PUSHF   R2          ; Save R2: floating point part
	PUSHF   R0          ; Save u (dividend) 
	LDI     R1,AR0      ; Save mantissa of v to remember sign
        ABSF    R1          ; The algorithm uses v = |v|.
;
;   Extract the exponent of v.
;
        PUSHF   R1
        POP     R2
        ASH     -24,R2      ; The 8 LSBs of R2 contain the exponent of v.
;
; A few comments on boundary conditions.  If e = -128, then v = 0.  The
; following x[0] calculation yields R2 = --128 - 1 = 127 and the algorithm will
; overflow and saturate since x[0] is large.  This seems reasonable.  If e =
; 127, the R2 = -127 - 1 = -128.  Thus x[0] = 0 and this will cause the
; algorithm to yield zero.  Since the mantissa of v is always between 1 and 2,
; this is also reasonable.  As a result, boundary conditions are handled
; automatically in a reasonable fashion.
;
;   x[0] formation given the exponent of v.
;
        NEGI    R2
        SUBI    1,R2            ; Now we have -e-1, the exponent of x[0].
        ASH     24,R2
        PUSH    R2
        POPF    R2              ; Now R2 = x[0] = 1.0 * 2**(-e-1).
;
; Now the iterations begin.
;
	LDI	3, RC
        MPYF    R2,R1,R0        ; R0 = v * x[0]
        SUBRF   2.0,R0          ; R0 = 2.0 - v * x[0]
        MPYF    R0,R2           ; R2 = x[1] = x[0] * (2.0 - v * x[0])
 
        MPYF    R2,R1,R0        ; R0 = v * x[1]
        SUBRF   2.0,R0          ; R0 = 2.0 - v * x[1]
        MPYF    R0,R2           ; R2 = x[2] = x[1] * (2.0 - v * x[1])
 
        MPYF    R2,R1,R0        ; R0 = v * x[2]
        SUBRF   2.0,R0          ; R0 = 2.0 - v * x[2]
        MPYF    R0,R2           ; R2 = x[3] = x[2] * (2.0 - v * x[2])
* 
        MPYF    R2,R1,R0        ; R0 = v * x[3]
        SUBRF   2.0,R0          ; R0 = 2.0 - v * x[3]
        MPYF    R0,R2           ; R2 = x[4] = x[3] * (2.0 - v * x[3])
 
        RND     R2              ; This minimizes error in the LSBs.
;
; For the last iteration we use the formulation:
; x[5] = (x[4] * (1.0 - (v * x[4]))) + x[4]
;
        MPYF    R2,R1,R0        ; R0 = v * x[4] = 1.0..01.. => 1
        SUBRF   1.0,R0          ; R0 = 1.0 - v * x[4] = 0.0..01... => 0
        MPYF    R2,R0           ; R0 = x[4] * (1.0 - v * x[4])
        ADDF    R0,R2,R1        ; R0 = x[5] = (x[4]*(1.0-(v*x[4])))+x[4]
;
; R1 contains 1/v.  Multiply by u to make result.
;
        RND     R1              ; Round since this is follow by a MPYF.
	POPF    R0              ; Pop u
	MPYF    R1,R0           ; Result = u * (1/v)
;
; Branch (delayed) return.  Use delay slots to negate the result if v < 0.
;
	POPF    R2              ; Restore R2: floating point part
	POP     R2              ; Restore R2: integer part

	BD      AR1             ; Delayed branch to return
        NEGF    R0,R1           ; R1 = -(1/|v|)
        CMPI    0,AR0           ; See if v was negative
        LDFN    R1,R0           ; If v < 0, then R1 = -R1
***     B       AR1             ; BRANCH OCCURS (RETURN)

*****************************************************
*	.def _cos_cram
;
; Initialization: get arguement, setup data pointers, and save registers
;
_cos_cram:
	POP	AR2		      ;return address -> AR2
	LDI     SP,AR0      ;setup frame pointer
        LDF     *-AR0(0),R2           ;x -> R2
	LDP	COS_ADR 	      ;save data page
        ABSF    R2                    ;Y = absolute value of x
        LDI     @COS_ADR,AR0
;
; Compute the result
;
        ADDF    *AR0++,R2,R1          ;Y += Pi / 2
        MPYF    *AR0++,R1             ;Y / PI -> R1
        FIX     R1,R0                 ;N = integer R1
        FLOAT   R0                    ;XN = float N
        SUBF    R0,R1                 ;R1 - XN -> R1
	CMPF	0.5,R1		      ;compare R1 to 0.5
	LDFNN	1.0,R1		      ;if R1 >= 0.5, 1 -> R1
	LDFN	0.0,R1		      ;if R1 < 0.5, 0 -> R1
        ADDF    R0,R1                 ;R2 + R1 to round R2
        FIX     R1,R0                 ;XN = integer R1
        TSTB    1,R0                  ;logical AND XN and 1
	LDINZ	-1,AR1		      ;if XN is odd, - 1 -> R3
	LDIZ	1,AR1		      ;if XN is even, 1 -> R3
	SUBF	0.5,R1		      ;sign *= R3
        MPYF    *AR0++,R1,R3          ;3.140625 * XN -> R3
        SUBF    R3,R2                 ;Y - R3 -> R3
        MPYF    *AR0++,R1,R3          ;9.67653589796e-4 * XN -> R2
        SUBF    R3,R2,R0              ;f = Y - XN * PI
	MPYF	R0,R0,R3	      ;g = f * f
        MPYF    *AR0++,R3,R2          ;0.2601903036e-5 * g -> R2
        ADDF    *AR0++,R2             ;-0.1980741872e-3 + R2 -> R2
	MPYF	R3,R2		      ;R2 * g -> R2
        ADDF    *AR0++,R2             ;0.8333025739e-2 + R2 -> R2
	MPYF	R3,R2		      ;R2 * g -> R2
        ADDF    *AR0++,R2             ;-0.1666665668 + R2 -> R2
	MPYF	R3,R2		      ;g * R2 -> R2
	MPYF	R0,R2		      ;f * R2 -> R2
	BD	AR2		      ;return from routine
	ADDF	R0,R2		      ;result = f + R2
	FLOAT	AR1,R1		      ;sign -> R1
	MPYF	R2,R1,R0	      ;result *= sign

***********************************************************************
*  DEFINE CONSTANTS FOR COSINE
***********************************************************************
COS_ADR:        .word   COS

COS             .float  1.570796326794896
		.float	0.318309886183790
		.float	3.140625
		.float	9.67653589793e-4
		.float	0.2601903036e-5
		.float -0.1980741872e-3
		.float	0.8333025139e-2
		.float -0.1666665668

*****************************************************

*
*	.label	end_update	; label end of load time address.

	.end

