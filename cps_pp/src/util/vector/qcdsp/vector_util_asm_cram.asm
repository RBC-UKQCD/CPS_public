**--------------------------------------------------------------------
**  CVS keywords
**
**  $Author: zs $
**  $Date: 2004-08-18 11:58:09 $
**  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/vector/qcdsp/vector_util_asm_cram.asm,v 1.4 2004-08-18 11:58:09 zs Exp $
**  $Id: vector_util_asm_cram.asm,v 1.4 2004-08-18 11:58:09 zs Exp $
**  $Name: not supported by cvs2svn $
**  $Locker:  $
**  $Revision: 1.4 $
**  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/vector/qcdsp/vector_util_asm_cram.asm,v $
**  $State: Exp $
**
**--------------------------------------------------------------------
*
*  vector_util_asm_cram.asm
*  To be relocated to chip ram
*

	.version        30

	.global		_moveMem
	.global		_dotProduct
	.global		_fTimesV1PlusV2


	.text

;=====================================================================
; moveMem(float *b, float *a, int len)

	.def	_moveMem

_moveMem:
	POP	BK

	POP	AR0
	POP	AR1
	POP	RC
	ADDI	3,SP


	LDI	*AR1++,	R0

	SUBI	4,RC
	RPTB	loop

loop:	LDI	*AR1++, R0
    ||	STI	R0,*AR0++


	BUD	BK
	LDI	*AR1++, R0
    ||	STI	R0,*AR0++
	LDI	*AR1++, R0
    ||	STI	R0,*AR0++
	STI	R0,*AR0


;=====================================================================
; float dotProduct(const float *a, const float *b, int n)	     |
;								     |
;   Total words of code    : 18					     |
;   Volatile registers used: R0,R1,R2,AR0,AR1,DP,BK,RS,RE,RC	     |
;   Registers for locals   : R2		holds sum		     |
;   Parameters             : AR0	holds a			     |
;			     AR1	holds b			     |
;			     R0		holds n			     |
;   Stack frame            : quick (AR3 points to some old frame)    |
;=====================================================================
	.def	_dotProduct

_dotProduct:
	POP	BK

	POP	AR0
	POP	AR1
	POP	RC
	ADDI	3,SP

	.word	40628000h	;	LDFU	0.0,R2
	SUBI	3,RC
	MPYF3	*AR1++(1),*AR0++(1),R1
	RPTB	L53

; begin loop 1-1
L53:	MPYF3	*AR1++(1),*AR0++(1),R1
     || ADDF3	R1,R2,R2
; end loop 1-1

	BUD	BK
	MPYF3	*AR1++(1),*AR0++(1),R1
     || ADDF3	R1,R2,R2
	ADDF	R1,R2
	RND	R2,R0


;=====================================================================
; fTimesV1PlusV2(float *a, float b, const float *c,		     |
;	const float *d, int n)					     |
;								     |
;   Total words of code    : 16					     |
;   Volatile registers used: R0,R1,R2,AR0,AR1,AR2,DP,BK,RS,RE,RC     |
;   Parameters             : AR0	holds a			     |
;			     R0		holds b			     |
;			     AR1	holds c			     |
;			     AR2	holds d			     |
;			     R1		holds n			     |
;   Stack frame            : quick (AR3 points to some old frame)    |
;=====================================================================
	.def	_fTimesV1PlusV2

_fTimesV1PlusV2:
	POP	BK

	POP	AR0		; a
	POPF	R1		; b
	POP	AR1		; c
	POP	AR2		; d
	POP	RC		; n
	ADDI	5,SP

	MPYF3	*AR1++(1),R1,R3

	SUBI	2,RC
	RPTB	L101

; begin loop 5-1
	ADDF3	*AR2++(1),R3,R0
	RND	R0, R2
L101:	MPYF3	*AR1++(1),R1,R3
    ||	STF	R2,*AR0++(1)
; end loop 5-1

	BUD	BK
	ADDF3	*AR2,R3,R0
	RND	R0, R2
    	STF	R2,*AR0
