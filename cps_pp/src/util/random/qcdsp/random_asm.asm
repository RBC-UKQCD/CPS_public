**--------------------------------------------------------------------
**  CVS keywords
**
**  $Author: zs $
**  $Date: 2004-08-18 11:58:07 $
**  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/random/qcdsp/random_asm.asm,v 1.4 2004-08-18 11:58:07 zs Exp $
**  $Id: random_asm.asm,v 1.4 2004-08-18 11:58:07 zs Exp $
**  $Name: not supported by cvs2svn $
**  $Locker:  $
**  $Revision: 1.4 $
**  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/random/qcdsp/random_asm.asm,v $
**  $State: Exp $
**
**--------------------------------------------------------------------
	.version 30
;=====================================================================
; This file is   : random_asm.asm				     |
; generated from : random_asm					     |
; on (M/D/Y H:M) : 08/26/1997 05:32				     |
;--------------------------------------------------------------------|


	.global	_ma__15RandomGenerator		; !C COMMON!
	.global	_inext__15RandomGenerator	; !C COMMON!
	.global	_inextp__15RandomGenerator	; !C COMMON!
	.def	_Rand__15RandomGeneratorFv
	.def	_Rand__15RandomGeneratorFv$LAJ
;=====================================================================
; LOCAL OPTIONS							     |
;   Calling convention     : TI C stack parameters		     |
;								     |
; GENERATED CODE PROPERTIES					     |
;   Total words of code    : 54					     |
;   Volatile registers used: R0,R1,AR0,AR1,DP,BK,RC		     |
;   Registers for locals   : RC		holds mj		     |
;   Parameters             : Never Used holds TEMP4		     |
;   Stack frame            : quick (AR3 points to some old frame)    |
;=====================================================================

	.sect	"T:_Rand__15RandomGeneratorFv"

_Rand__15RandomGeneratorFv$LAJ:
	PUSH	BK
_Rand__15RandomGeneratorFv:
	PUSH	DP

	LDP	_inext__15RandomGenerator,AR0
	LSH	16,AR0
	OR	_inext__15RandomGenerator,AR0
	LDIU	*AR0,R0
	ADDI	1,R0
	CMPI	55,R0
	BNED	L40
	LDP	_inext__15RandomGenerator,DP
	LDIU	0,R1
	STI	R0,@_inext__15RandomGenerator
	LDP	_inext__15RandomGenerator,DP
	STI	R1,@_inext__15RandomGenerator
L40:	LDP	_inextp__15RandomGenerator,AR0
	LSH	16,AR0
	OR	_inextp__15RandomGenerator,AR0
	LDIU	*AR0,R0
	ADDI	1,R0
	CMPI	55,R0
	BNED	L41
	LDP	_inextp__15RandomGenerator,DP
	STI	R0,@_inextp__15RandomGenerator
	LDIEQ	(_inextp__15RandomGenerator)>>16,DP
	LDIU	0,R1
	STI	R1,@_inextp__15RandomGenerator
L41:	LDP	_inext__15RandomGenerator,AR0
	LSH	16,AR0
	OR	_inext__15RandomGenerator,AR0
	LDP	_ma__15RandomGenerator,DP
	LSH	16,DP
	LDIU	*AR0,AR0
	OR	_ma__15RandomGenerator,DP
	ADDI	DP,AR0
	LDP	_inextp__15RandomGenerator,AR1
	LSH	16,AR1
	OR	_inextp__15RandomGenerator,AR1
	LDP	_ma__15RandomGenerator,DP
	LSH	16,DP
	LDIU	*AR1,AR1
	OR	_ma__15RandomGenerator,DP
	ADDI	DP,AR1
	SUBI3	*AR1,*AR0,RC
	CMPI	0,RC
	BGED	L42
	LDIU	15258,R0
	LSH	16,R0
	OR	51712,R0
	ADDI	R0,RC
L42:	STI	RC,*AR0

	LDP	FAC,DP
	LDF	@FAC,R0

	POP	DP
	POP	BK

	BUD	BK
	FLOAT	RC,R1
	MPYF	R0,R1
	RND	R1,R0



	.def	_Rand__22UniformRandomGeneratorFv
	.def	_Rand__22UniformRandomGeneratorFv$LAJ
;=====================================================================
; LOCAL OPTIONS							     |
;   Calling convention     : TI C stack parameters		     |
;								     |
; GENERATED CODE PROPERTIES					     |
;   Total words of code    : 13					     |
;   Volatile registers used: R0,R1,AR0,AR1,AR2,DP,BK,RC		     |
;   Parameters             : AR0	holds TEMP6		     |
;   Stack frame            : quick (AR3 points to some old frame)    |
;=====================================================================

	.sect	"T:_Rand__22UniformRandomGeneratorFv"

_Rand__22UniformRandomGeneratorFv$LAJ:
	PUSH	BK
_Rand__22UniformRandomGeneratorFv:
	LDIU	SP,AR0
	LDIU	*-AR0(1),AR0

	LDIU	AR0,AR2

	PUSH	AR2
	CALL	_Rand__15RandomGeneratorFv
	SUBI	1,SP

	LDFU	*+AR2(2),R1

	SUBF	*+AR2(1),R1
	RND	R1

	POP	BK

	BUD	BK
	MPYF	R0,R1
	ADDF3	R1,*+AR2(1),R0
	RND	R0




	.def	_Rand__23GaussianRandomGeneratorFv
	.ref	DEFALT
	.ref	ARTALOG32
	.ref	ARTDIVF32UZ
	.ref	ARTSQRT32
	.def	_Rand__23GaussianRandomGeneratorFv$LAJ
;=====================================================================
; LOCAL OPTIONS							     |
;   Calling convention     : TI C stack parameters		     |
;								     |
; GENERATED CODE PROPERTIES					     |
;   Total words of code    : 56					     |
;   Registers for locals   : R6		holds v2		     |
;			     R7		holds rsq, fac		     |
;   Parameters             : *-AR3(2)   holds TEMP8		     |
;   Stack frame            : full (frame pointer in AR3)	     |
;=====================================================================

	.sect	"T:_Rand__23GaussianRandomGeneratorFv"

_Rand__23GaussianRandomGeneratorFv$LAJ:
	PUSH	BK
_Rand__23GaussianRandomGeneratorFv:
	PUSH	AR3
	LDIU	SP,AR3
	ADDI	1,SP
	PUSHF	R6
	PUSHF	R7

	PUSH	DP

	LDIU	*-AR3(2),AR2
	LDI	*+AR2(2),R0
	BNE	L34

; begin loop 1-1
L36:

	PUSH	AR2
	CALL	_Rand__15RandomGeneratorFv

	ADDF	R1,R1,R0
	.word	17e00000h	;	SUBF	0.1e1,R0
	RND	R0
	STF	R0,*+AR3(1)

	PUSH	AR2
	CALL	_Rand__15RandomGeneratorFv
	SUBI	2,SP

	ADDF	R1,R1,R0
	.word	40610000h	;	LDFU	0.1e1,R1
	SUBF3	R1,R0,R6
	LDFU	*+AR3(1),R0
	MPYF	R0,R0
	RND	R6
	MPYF3	R6,R6,R1
	ADDF3	R1,R0,R7
	.word	4670000h	;	CMPF	0.1e1,R7
	BGE	L36
	.word	4678000h	;	CMPF	0.0,R7
	BEQ	L36
; end loop 1-1


	RND	R7,R0
	LDP	DEFALT,DP
	CALL	ARTALOG32
	RND	R0

	LDFU	*+AR2(1),R1
	.word	0a610800h	;	MPYF	-0.2e1,R1
	RND	R1
	MPYF	R1,R0
	RND	R0

	RND	R7,R1
	LDP	DEFALT,DP
	CALL	ARTDIVF32UZ
	RND	R0


	LDP	DEFALT,DP
	CALL	ARTSQRT32
	RND	R0,R7

	MPYF	*+AR3(1),R0
	RND	R0
	STF	R0,*+AR2(3)
	LDIU	1,R1
	BRD	L54
	STI	R1,*+AR2(2)
	MPYF3	R7,R6,R0
	RND	R0


L34:	LDIU	0,R0
	STI	R0,*+AR2(2)
	LDFU	*+AR2(3),R0
L54:

	POP	DP

	LDFU	*+AR3(3),R7
	LDFU	*+AR3(2),R6
	LDIU	*-AR3(1),BK
	BUD	BK
	LDIU	AR3,SP
	LDIU	*AR3,AR3
	SUBI	2,SP







	.data
GLOBAL:	.word	0
	.word	0
	.word	0
	.word	0
	.word	0
	.word	_Rand__22UniformRandomGeneratorFv
	.word	0
	.word	0
	.word	0


	.word	0
	.word	0
	.word	0
	.word	0
	.word	0
	.word	_Rand__23GaussianRandomGeneratorFv
	.word	0
	.word	0
	.word	0


	.word	0
	.word	0
	.word	0
	.word	0
	.word	0
	.word	_Rand__15RandomGeneratorFv
	.word	0
	.word	0
	.word	0


	.data
OWNX:	.word	0e209705fh	;0.1e-8


	.space	2
	.sect	"o:DEFAULT"
DEFAULT:	.word	0
	.def	___vtbl__15RandomGenerator
___vtbl__15RandomGenerator	.set	GLOBAL+18
	.def	___vtbl__22UniformRandomGenerator
___vtbl__22UniformRandomGenerator	.set	GLOBAL
	.def	___vtbl__23GaussianRandomGenerator
___vtbl__23GaussianRandomGenerator	.set	GLOBAL+9

MBIG	.set	OWNX+1

MZ	.set	OWNX+2

FAC	.set	OWNX




	.end

