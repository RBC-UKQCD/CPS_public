**--------------------------------------------------------------------
**  CVS keywords
**
**  $Author: zs $
**  $Date: 2004-08-18 11:57:39 $
**  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/alg/alg_ghb/qcdsp/metro_kern.asm,v 1.4 2004-08-18 11:57:39 zs Exp $
**  $Id: metro_kern.asm,v 1.4 2004-08-18 11:57:39 zs Exp $
**  $Name: not supported by cvs2svn $
**  $Locker:  $
**  $Revision: 1.4 $
**  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/alg/alg_ghb/qcdsp/metro_kern.asm,v $
**  $State: Exp $
**
**--------------------------------------------------------------------
******************************************************
*    TMS320C30 C COMPILER     Version 4.50
******************************************************
	.version	30
FP	.set		AR3

	.asg	1,	ROUND

*  Make the name of this routine global
	.def	_metropolis_kernel__FP6rfloatT1; publish run time address
*	.def	src_update	; publish first load time address
*	.def	end_update	; publish last load time address

*  Some symbols which might be useful
*	.ref	_sqrt
	.ref	_exp
*	.ref	INV_F30

	.ref	_grand
	.ref	_m_multiply3
	.ref	_m_equal
	.ref	_m_identity

*  External variables
	.ref	_core_iscratch
	.ref	_core_fscratch

******************************************************
* FUNCTION DEF : _metropolis_kernel__FP6rfloatT1
******************************************************
*	.sect	".update"
	.text
*	.label	src_update	;	label load time address
_metropolis_kernel__FP6rfloatT1
	PUSH	FP
	LDI	SP,FP
	ADDI	3,SP
	PUSH	R4
	PUSH	R5
	PUSHF	R6
	PUSHF	R7
	PUSH	AR4
	.globl	_grand
*
* R2	assigned to temp var  C$2
* R2	assigned to temp var  C$4
* R2	assigned to temp var  C$7
* R2	assigned to variable  die_roll
* R2	assigned to variable  norm
* R2	assigned to variable  select
* R4	assigned to parameter U
* R5	assigned to parameter sigma
* R6	assigned to temp var  K$19
* R6	assigned to temp var  K$42
* R7	assigned to variable  a
* R7	assigned to variable  accept_probability
* R7	assigned to variable  old_action
* AR2	assigned to temp var  C$1
* AR2	assigned to temp var  C$3
* AR2	assigned to temp var  C$6
* AR2	assigned to temp var  C$8
* AR4	assigned to temp var  C$5
* AR4	assigned to variable  v
*
	LDI	*-FP(2),R5
	LDI	*-FP(3),R4
*** 11	-----------------------    core_fscratch[58] = 2.0e-1F;
	LDP	@_core_fscratch
	LDI	@_core_fscratch,AR0
	LDP	@METRO_CONST+0
	LDF	@METRO_CONST+0,R0
	STF	R0,*+AR0(58)
*** 12	-----------------------    *core_iscratch = 10;
	LDP	@_core_iscratch
	LDI	@_core_iscratch,AR1
	LDI	10,R1
	STI	R1,*AR1
*** 14	-----------------------    core_iscratch[1] = 0;
	LDI	0,R2
	STI	R2,*+AR1(1)
*** 14	-----------------------    C$8 = core_iscratch;
	LDP	@_core_iscratch
	LDI	@_core_iscratch,AR2
*** 14	-----------------------    if ( C$8[1] >= *C$8 ) goto g11;
	CMPI	*AR2,*+AR2(1)
	BGE	EPI0_1
***  	-----------------------    K$42 = K$19 = 0.0F;
	LDF	0.0,R6
L3:
***	-----------------------g3:
*** 21	-----------------------    v = core_fscratch;
	LDP	@_core_fscratch
	LDI	@_core_fscratch,AR4
*** 23	-----------------------    a = core_fscratch[58]*grand();
	CALL	_grand
	LDP	@_core_fscratch
	LDI	@_core_fscratch,AR0
	MPYF	*+AR0(58),R0
	LDF	R0,R7
*** 23	-----------------------    a = (a > K$19) ? -a : a;
	CMPF	R6,R0
	BLE	LL3
	NEGF	R0,R7
LL3:
*** 23	-----------------------    b = grand();
	CALL	_grand
	STF	R0,*+FP(1)
*** 25	-----------------------    c = grand();
	CALL	_grand
	STF	R0,*+FP(2)
*** 26	-----------------------    d = grand();
	CALL	_grand
	STF	R0,*+FP(3)
*** 23	-----------------------    a += C$7 = 1.0F;
	LDF	1.0,R2
	ADDF	R2,R7
*** 28	-----------------------    norm = sqrt((C$7-a*a)/(b*b+c*c+d*d));
	MPYF	*+FP(1),*+FP(1),R0
	LDF	*+FP(2),R1
	MPYF	R1,R1,R3
	ADDF	R3,R0
	LDF	*+FP(3),R3
	MPYF	R3,R3
	ADDF	R3,R0
	CALL	INV_F30_cram
	RND	R0
	MPYF	R7,R7,R1
	SUBF	R1,R2,R1
	MPYF	R1,R0
	PUSHF	R0
	CALL	_sqrt_cram
	SUBI	1,SP
	LDF	R0,R2
*** 29	-----------------------    b *= norm;
	MPYF	R0,*+FP(1),R1
	.if ROUND
	RND	R1
	.endif
	STF	R1,*+FP(1)
*** 29	-----------------------    c *= norm;
	MPYF	*+FP(2),R0
	.if ROUND
	RND	R0
	.endif
	STF	R0,*+FP(2)
*** 29	-----------------------    d *= norm;
	MPYF	*+FP(3),R2
	.if ROUND
	RND	R2
	.endif
	STF	R2,*+FP(3)
*** 31	-----------------------    select = grand();
	CALL	_grand
	LDF	R0,R2
*** 33	-----------------------    if ( select < 3.3333333e-1F ) goto g7;
	LDP	@METRO_CONST+1
	CMPF	@METRO_CONST+1,R0
	BLT	L7
*** 42	-----------------------    if ( select > 6.6666666e-1F ) goto g6;
	LDP	@METRO_CONST+2
	CMPF	@METRO_CONST+2,R2
	BGT	L6
*** 53	-----------------------    *v++ = a;
	STF	R7,*AR4++
*** 53	-----------------------    *v++ = K$42;
	STF	R6,*AR4++
*** 53	-----------------------    *v++ = c;
	LDF	*+FP(2),R1
	STF	R1,*AR4++
*** 54	-----------------------    *v++ = K$42;
	STF	R6,*AR4++
*** 54	-----------------------    *v++ = 1.0F;
	LDF	1.0,R2
	STF	R2,*AR4++
*** 54	-----------------------    *v++ = K$42;
	STF	R6,*AR4++
*** 55	-----------------------    *v++ = -c;
	NEGF	R1,R3
	STF	R3,*AR4++
*** 55	-----------------------    *v++ = K$42;
	STF	R6,*AR4++
*** 55	-----------------------    *v++ = a;
	STF	R7,*AR4++
*** 57	-----------------------    *v++ = b;
	LDF	*+FP(1),R3
	STF	R3,*AR4++
*** 57	-----------------------    *v++ = K$42;
	STF	R6,*AR4++
*** 57	-----------------------    *v++ = d;
	LDF	*+FP(3),R3
	STF	R3,*AR4++
*** 58	-----------------------    *v++ = K$42;
	STF	R6,*AR4++
*** 58	-----------------------    *v++ = K$42;
	STF	R6,*AR4++
*** 58	-----------------------    *v++ = K$42;
	STF	R6,*AR4++
*** 59	-----------------------    *v++ = d;
	STF	R3,*AR4++
*** 59	-----------------------    *v++ = K$42;
	BD	L8
	STF	R6,*AR4++
*** 59	-----------------------    *v++ = -b;
	NEGF	*+FP(1),R3
	STF	R3,*AR4++
*** 59	-----------------------    goto g8;
***	B	L8	;BRANCH OCCURS
L6:
***	-----------------------g6:
*** 44	-----------------------    *v++ = 1.0F;
	LDF	1.0,R1
	STF	R1,*AR4++
*** 44	-----------------------    *v++ = K$42;
	STF	R6,*AR4++
*** 44	-----------------------    *v++ = K$42;
	STF	R6,*AR4++
*** 45	-----------------------    *v++ = K$42;
	STF	R6,*AR4++
*** 45	-----------------------    *v++ = a;
	STF	R7,*AR4++
*** 45	-----------------------    *v++ = c;
	LDF	*+FP(2),R2
	STF	R2,*AR4++
*** 46	-----------------------    *v++ = K$42;
	STF	R6,*AR4++
*** 46	-----------------------    *v++ = -c;
	NEGF	R2,R3
	STF	R3,*AR4++
*** 46	-----------------------    *v++ = a;
	STF	R7,*AR4++
*** 48	-----------------------    *v++ = K$42;
	STF	R6,*AR4++
*** 48	-----------------------    *v++ = K$42;
	STF	R6,*AR4++
*** 48	-----------------------    *v++ = K$42;
	STF	R6,*AR4++
*** 49	-----------------------    *v++ = K$42;
	STF	R6,*AR4++
*** 49	-----------------------    *v++ = b;
	LDF	*+FP(1),R3
	STF	R3,*AR4++
*** 49	-----------------------    *v++ = d;
	LDF	*+FP(3),R3
	STF	R3,*AR4++
*** 50	-----------------------    *v++ = K$42;
	STF	R6,*AR4++
*** 50	-----------------------    *v++ = d;
	BD	L8
	STF	R3,*AR4++
*** 50	-----------------------    *v++ = -b;
	NEGF	*+FP(1),R3
	STF	R3,*AR4++
***  	-----------------------    goto g8;
***	B	L8	;BRANCH OCCURS
L7:
***	-----------------------g7:
*** 34	-----------------------    *v++ = a;
	STF	R7,*AR4++
*** 34	-----------------------    *v++ = c;
	LDF	*+FP(2),R1
	STF	R1,*AR4++
*** 34	-----------------------    *v++ = K$42;
	STF	R6,*AR4++
*** 35	-----------------------    *v++ = -c;
	NEGF	R1,R2
	STF	R2,*AR4++
*** 35	-----------------------    *v++ = a;
	STF	R7,*AR4++
*** 35	-----------------------    *v++ = K$42;
	STF	R6,*AR4++
*** 36	-----------------------    *v++ = K$42;
	STF	R6,*AR4++
*** 36	-----------------------    *v++ = K$42;
	STF	R6,*AR4++
*** 36	-----------------------    *v++ = 1.0F;
	LDF	1.0,R2
	STF	R2,*AR4++
*** 38	-----------------------    *v++ = b;
	LDF	*+FP(1),R3
	STF	R3,*AR4++
*** 38	-----------------------    *v++ = d;
	LDF	*+FP(3),R3
	STF	R3,*AR4++
*** 38	-----------------------    *v++ = K$42;
	STF	R6,*AR4++
*** 39	-----------------------    *v++ = d;
	STF	R3,*AR4++
*** 39	-----------------------    *v++ = -b;
	NEGF	*+FP(1),R3
	STF	R3,*AR4++
*** 39	-----------------------    *v++ = K$42;
	STF	R6,*AR4++
*** 40	-----------------------    *v++ = K$42;
	STF	R6,*AR4++
*** 40	-----------------------    *v++ = K$42;
	STF	R6,*AR4++
*** 40	-----------------------    *v++ = K$42;
	STF	R6,*AR4++
L8:
***	-----------------------g8:
*** 64	-----------------------    m_multiply3(core_fscratch+18, sigma, U);
	PUSH	R4
	PUSH	R5
	LDP	@_core_fscratch
	LDI	@_core_fscratch,R1
	ADDI	18,R1
	PUSH	R1
	CALL	_m_multiply3
	SUBI	3,SP
*** 65	-----------------------    C$6 = core_fscratch;
	LDP	@_core_fscratch
	LDI	@_core_fscratch,AR2
*** 65	-----------------------    C$5 = &C$6[18];
	LDI	AR2,AR4
	ADDI	18,AR4
*** 65	-----------------------    old_action = C$6[57]*(9.0F-(*C$5+C$6[22]+C$6[26])*5.0e-1F);
	LDF	*+AR2(22),R7
	ADDF	*AR4,R7
	ADDF	*+AR2(26),R7
	MPYF	5.0e-1,R7
	SUBRF	9.0,R7
	MPYF	*+AR2(57),R7
*** 67	-----------------------    m_multiply3(C$5, C$6, U);
	PUSH	R4
	PUSH	AR2
	PUSH	AR4
	CALL	_m_multiply3
	SUBI	3,SP
*** 68	-----------------------    C$4 = core_fscratch;
	LDP	@_core_fscratch
	LDI	@_core_fscratch,R2
*** 68	-----------------------    m_multiply3(C$4+36, sigma, C$4+18);
	LDI	R2,R0
	ADDI	18,R0
	PUSH	R0
	PUSH	R5
	ADDI	36,R2
	PUSH	R2
	CALL	_m_multiply3
	SUBI	3,SP
*** 71	-----------------------    C$3 = core_fscratch;
	LDP	@_core_fscratch
	LDI	@_core_fscratch,AR2
*** 71	-----------------------    accept_probability = exp(old_action-C$3[57]*(9.0F-(C$3[36]+C$3[40]+C$3[44])*5.0e-1F));
	LDF	*+AR2(36),R0
	ADDF	*+AR2(40),R0
	ADDF	*+AR2(44),R0
	MPYF	5.0e-1,R0
	SUBRF	9.0,R0
	MPYF	*+AR2(57),R0
	SUBF	R0,R7,R0
	PUSHF	R0
	CALL	_exp_cram
	SUBI	1,SP
	LDF	R0,R7
*** 73	-----------------------    die_roll = grand();
	CALL	_grand
	LDF	R0,R2
*** 74	-----------------------    die_roll = (die_roll < K$19) ? -die_roll : die_roll;
	CMPF	R6,R0
	BGE	LL5
	NEGF	R0,R2
LL5:
*** 74	-----------------------    if ( die_roll >= accept_probability ) goto g10;
	CMPF	R7,R2
	BGE	L10
*** 77	-----------------------    C$2 = core_fscratch;
	LDP	@_core_fscratch
	LDI	@_core_fscratch,R2
*** 77	-----------------------    m_multiply3(C$2+18, C$2, U);
	PUSH	R4
	PUSH	R2
	ADDI	18,R2
	PUSH	R2
	CALL	_m_multiply3
	SUBI	3,SP
*** 78	-----------------------    m_equal(U, core_fscratch+18);
	LDP	@_core_fscratch
	LDI	@_core_fscratch,R0
	ADDI	18,R0
	PUSH	R0
	PUSH	R4
	CALL	_m_equal
	SUBI	2,SP
L10:
***	-----------------------g10:
*** 14	-----------------------    ++(core_iscratch[1]);
	LDP	@_core_iscratch
	LDI	@_core_iscratch,AR0
	LDI	1,R0
	ADDI	R0,*+AR0(1),R1
	STI	R1,*+AR0(1)
*** 14	-----------------------    C$1 = core_iscratch;
	LDP	@_core_iscratch
	LDI	@_core_iscratch,AR2
*** 14	-----------------------    if ( C$1[1] < *C$1 ) goto g3;
	CMPI	*AR2,*+AR2(1)
	BLT	L3
***	-----------------------g11:
***  	-----------------------    return;
EPI0_1:
	.if ROUND
	RND	R0
	.endif
	LDI	*-FP(1),R1
	LDI	*FP,FP
	POP	AR4
	POPF	R7
	POPF	R6
	BD	R1
	POP	R5
	POP	R4
	SUBI	5,SP
***	B	R1	;BRANCH OCCURS
******************************************************
* DEFINE METRO_CONSTANTS                             *
******************************************************

METRO_CONST:
	.float	2.0e-1           ;0
	.float	3.3333333e-1     ;1
	.float	6.6666666e-1     ;2

*******************************************************************
	.def	faststack_p_m
faststack_p_m	.word	faststack_m	; label run time address
faststack_m	.space	60 ; 27
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
	.ref	_bad_sqrt
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

INV_F30_cram:  
	POP     AR1         ; Pop return address
	PUSH    R2          ; Save R2: integer part
	PUSHF   R2          ; Save R2: floating point part
	LDI     R0,AR0      ; Save mantissa of v to remember sign
        ABSF    R0          ; The algorithm uses v = |v|.
;
;   Extract the exponent of v.
;
        PUSHF   R0
        POP     R1
        ASH     -24,R1      ; The 8 LSBs of R1 contain the exponent of v.
;
; A few comments on boundary conditions.  If e = -128, then v = 0.  The
; following x[0] calculation yields R1 = --128 - 1 = 127 and the algorithm will
; overflow and saturate since x[0] is large.  This seems reasonable.  If e =
; 127, the R1 = -127 - 1 = -128.  Thus x[0] = 0 and this will cause the
; algorithm to yield zero.  Since the mantissa of v is always between 1 and 2,
; this is also reasonable.  As a result, boundary conditions are handled
; automatically in a reasonable fashion.
;
;   x[0] formation given the exponent of v.
;
        NEGI    R1
        SUBI    1,R1            ; Now we have -e-1, the exponent of x[0].
        ASH     24,R1
        PUSH    R1
        POPF    R1              ; Now R1 = x[0] = 1.0 * 2**(-e-1).
;
; Now the iterations begin.
;
        MPYF    R1,R0,R2        ; R2 = v * x[0]
        SUBRF   2.0,R2          ; R2 = 2.0 - v * x[0]
        MPYF    R2,R1           ; R1 = x[1] = x[0] * (2.0 - v * x[0])
 
        MPYF    R1,R0,R2        ; R2 = v * x[1]
        SUBRF   2.0,R2          ; R2 = 2.0 - v * x[1]
        MPYF    R2,R1           ; R1 = x[2] = x[1] * (2.0 - v * x[1])
 
        MPYF    R1,R0,R2        ; R2 = v * x[2]
        SUBRF   2.0,R2          ; R2 = 2.0 - v * x[2]
        MPYF    R2,R1           ; R1 = x[3] = x[2] * (2.0 - v * x[2])
 
        MPYF    R1,R0,R2        ; R2 = v * x[3]
        SUBRF   2.0,R2          ; R2 = 2.0 - v * x[3]
        MPYF    R2,R1           ; R1 = x[4] = x[3] * (2.0 - v * x[3])
 
        RND     R1              ; This minimizes error in the LSBs.
;
; For the last iteration we use the formulation:
; x[5] = (x[4] * (1.0 - (v * x[4]))) + x[4]
;
        MPYF    R1,R0,R2        ; R2 = v * x[4] = 1.0..01.. => 1
        SUBRF   1.0,R2          ; R2 = 1.0 - v * x[4] = 0.0..01... => 0
        MPYF    R1,R2           ; R2 = x[4] * (1.0 - v * x[4])
        ADDF    R2,R1,R0        ; R0 = x[5] = (x[4]*(1.0-(v*x[4])))+x[4]
;
; Return (delayed). Use delay slots to negate the result if v < 0.
;
	POPF    R2              ; Restore R2: floating point part
	POP     R2              ; Restore R2: integer part

	BD      AR1             ; delayed branch to return
        NEGF    R0,R1           ; R1 = -(1/v)
        CMPI    0,AR0           ; See if v was negative
        LDFN    R1,R0           ; If v < 0, then R1 = -R1
***     B       AR1             ; BRANCH OCCURS (RETURN)

*****************************************************************

        .globl _errno

ERANGE  .set    2

;
; Initialization: get arguements, setup data pointers, and save registers
;
_exp_cram:
	LDI	SP,AR0		 	;setup stack pointer
	LDF	*-AR0(1),R1		;x -> R1
*        .if .REGPARM == 0
*	SAFESP	"LDI     SP,AR0"      ;setup stack pointer
*        LDF     *-AR0(1),R1           ;x -> R1
*        .else
*        LDF     R2,R1                 ;x -> R1
*        .endif
	LDP	METRO_EXP_ADR		;load data page
	LDI	@METRO_EXP_ADR,AR0	;load data address in AR0
	PUSH	R4		      ;save integer value in R4
        PUSHF   R6                    ;save float value of R6
        CMPF    *+AR0(6),R1           ;if x >= MAXX
        BGED    EPI0_2                ;check for range error
;
; Compute R
;
	LDF	*AR0++,R4	      ;1 / ln(2) -> R4
        MPYF    R1,R4                 ;x / ln(2) -> R4
        ADDF    0.5,R4                ;round before truncating

        FIX     R4,R6                 ;integer part of R4 -> R6
        FLOAT   R6,R4                 ;make R4 into a floating point number
	LDF	*AR0++,R2	      ;0.693359375 -> R2
	MPYF	R4,R2		      ;R4 * R2 -> R2
	LDF	*AR0++,R3	      ;-2.1219444005469 -> R3
	MPYF	R4,R3		      ;R4 * R3 -> R3
	SUBF	R2,R1		      ;R1 - R2 -> R1
	SUBF	R3,R1		      ;g = x - (int)x * ln(2)
	MPYF	R1,R1,R2	      ;z = g * g
	LDF	*AR0++,R0	      ;0.41602886568e-2 -> R0
	MPYF	R2,R0		      ;z * R0 -> R0
	ADDF	*AR0++,R0	      ;P = 0.2499999995 + R0
	MPYF	R1,R0		      ;g * P -> R0
	LDF	*AR0++,R1	      ;0.49987178778e-1 -> R1
	MPYF	R2,R1		      ;z * R1 -> R1
	ADDF	0.5,R1		      ;Q = 0.5 + R1
	SUBF	R0,R1		      ;Q - g * P -> R1
				      ;g * P / ( Q - g * P) -> R0
	CALL    DIV_F30_cram_m
	ADDF	0.5,R0		      ;R = 0.5 + R0
        ADDI    1,R6                  ;(int)x + 1 -> R6
        PUSHF   R0                    ;convert R0 from floating -
        POP     R1                    ;point to integer format to extract -
        ASH     -24,R1                ;exponent: store in R1
        ADDI    R6,R1                 ;Add n + exponent of x
        ASH     24,R1                 ;restore the exponent to integer field -
        PUSH    R1                    ;convert back to floating point #
        POPF    R1                    ;store into R1
        LDE     R1,R0                 ;Load exponent of R1 into result
        POPF    R6                    ;restore R6
	POP	R4		      ;restore R4
        RETS
;
; Error service module
;
EPI0_2:
        LDP     _errno
        LDF     *+AR0(6),R0	      ;return HUGE_VAL
        POPF    R6                    ;restore R6
        POP     R4                    ;restore R4
        LDI     ERANGE,R1
        STI     R1,@_errno            ;set error flag = 2(range error)
        RETS
***********************************************************************
*  DEFINE CONSTANTS
***********************************************************************
METRO_EXP_ADR:       .word   METRO_EXP

METRO_EXP             .float  1.44269504088896
		.float	0.693359375
		.float -2.1219444005469e-4
		.float	0.41602886268e-2
		.float	0.249999995
		.float	0.49987178778e-1
                .float 88.72283906                     ; MAXX
                .float 3.4028235e+38                   ; HUGE_VAL

*******************************************************************

	.def	DIV_F30_cram_m
DIV_F30_cram_m:  
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

*	.label	end_update	; label end of load time address.

	.end
