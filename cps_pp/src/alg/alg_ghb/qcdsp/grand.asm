**--------------------------------------------------------------------
**  CVS keywords
**
**  $Author: chulwoo $
**  $Date: 2004-06-04 21:13:59 $
**  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/alg/alg_ghb/qcdsp/grand.asm,v 1.3 2004-06-04 21:13:59 chulwoo Exp $
**  $Id: grand.asm,v 1.3 2004-06-04 21:13:59 chulwoo Exp $
**  $Name: not supported by cvs2svn $
**  $Locker:  $
**  $Revision: 1.3 $
**  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/alg/alg_ghb/qcdsp/grand.asm,v $
**  $State: Exp $
**
**--------------------------------------------------------------------

	.asg	AR3, FP			; FP is AR3

	.def	_grand
	.ref	_fprint2
	.ref	_el_seed_p

*	.sect	".grand"
	.text
_grand:
	PUSH	DP

	LDP	@_el_seed_p		; set DP to _el_seed_p data page
	LDI	@_el_seed_p, AR0	; set AR0 to _el_seed_p

	LDI	*+AR0(1), R1
	LDI	*+AR0(0), R0
	ADDI	*+AR0(4), R0

	SUBRI	R1, R0
	LDIU	0, AR2			; new carry if no borrow 
	BNC	store_carry

	LDIU	1, AR2			; new carry if borrow
	SUBI	18, R0			; mod 2^32 - 18

store_carry:
	LDI	*+AR0(2), RC
	STI	R1, *+AR0(0)		; _el_seed11 = _el_seed12
	STI	AR2, *+AR0(4)
	STI	RC, *+AR0(1)		; _el_seed12 = _el_seed13
	STI	R0, *+AR0(2)		; _el_seed13 = s

*****
* next section based on MPY_I30 using:
*  69069*n = (3533*(n>>16) + 1*(n&0xFFFF))<<16 + 3533*(n&0xFFFF)
*****

	LDI	*+AR0(3), R1
	LDI	R1, AR1			; n
	LSH	-16, R1			; 16 MSBs of n
	AND	0FFFFh, AR1		; 16 LSBs of n
	LDI	AR1, R2			; 1*(n&0xFFFF)
	MPYI	3533, R1		; 3533*(n>>16)
	MPYI	3533, AR1		; 3533*(n&0xFFFF)
	ADDI	R2, R1
	LSH	16, R1			;
	ADDI	AR1, R1			;final result n = 69069*n

	LDI	*+AR0(5), R2
	ADDI	R2, R1			; n = 69069*n + lconst
	STI	R1, *+AR0(3)		; _el_seed21 = n
	ADDI	R1, R0			; combine generators
	XOR	R1, R1
	PUSH	R1
	POPF	R1
*	ABSI	R0, R0			; rem for [-1,1]
	OR	R0, R1
	NORM	R1, R0
	RND	R0, R0

	POP	DP
	RETS

	.end
