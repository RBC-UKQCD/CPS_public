#include<config.h>
CPS_START_NAMESPACE
**--------------------------------------------------------------------
**  CVS keywords
**
**  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/tests/various/pmalloc/main.asm,v $
**  $State: Exp $
**
**--------------------------------------------------------------------
	.version 30
;=============================================================================
; This file is   : main.asm						     |
; generated from : main							     |
; using          : Tartan C/C++ Compiler for the TMS320C3x/4x, v2.1.1	     |
; on (M/D/Y H:M) : 10/08/1997 12:54					     |
;----------------------------------------------------------------------------|
;									     |
; GLOBAL ATTRIBUTES AND ASSERTIONS					     |
;   Wait states            : Code = 0, Stack = 0, Data page = 0, Heap = 0    |
;   Max loop iterations    : 2**23 (no huge loops)			     |
;   Code Span              : 24-bits of address				     |
;									     |
; GLOBAL OPTIONS							     |
;   Optimization heuristics: default speed/space tradeoff		     |
;   Debug information      : off					     |
;   Profiling              : off					     |
;   Max RPTS allowed       : 2**31					     |
;   Data Page              : No_Page					     |
;   Data Page is in ROM    : false					     |
;   *ARn same mem as *ARm? : true					     |
;									     |
;=============================================================================

;
	.ref	_pmalloc__Fi
	.ref	_pfree__FPv
	.ref	_pclear__Fv
	.ref	_printf
	.def	_main
	.ref	__main
	.ref	TEMP5
	.ref	TEMP6
	.ref	TEMP11
	.ref	TEMP16
	.def	_main$LAJ
;=============================================================================
; LOCAL OPTIONS								     |
;   Calling convention     : TI C stack parameters			     |
;									     |
; GENERATED CODE PROPERTIES						     |
;   Total words of code    : 64						     |
;   Volatile registers used: R0,R1,R2,R3,R6,R7,AR0,AR1,AR2,DP,IR0,IR1,BK,ST, |
;			     IE,IF,IOF,RS,RE,RC				     |
;   Registers for locals   : AR7	holds p1			     |
;   Parameters             : *-AR3(2)   holds TEMP3			     |
;			     *-AR3(3)   holds TEMP4			     |
;   Stack frame            : full (frame pointer in AR3)		     |
;=============================================================================

	.sect	"T:_main"

_main$LAJ:
	PUSH	BK
_main:
	PUSH	AR3
	LDIU	SP,AR3
	ADDI	10,SP
	PUSH	AR6
	PUSH	AR7


	CALL	__main

	LDIU	AR3,AR7
	ADDI	1,AR7
	LDIU	9,AR6

	LDIU	1,R0

; begin loop 1-1
L40:	PUSH	R0
	CALL	_pmalloc__Fi

	STI	R0,*AR7++(1)

	PUSH	R0
	LDP	DEF3,R1
	LSH	16,R1
	OR	DEF3,R1
	PUSH	R1
	CALL	_printf

	DBUD	AR6,L40
	SUBI	3,SP
	LDIU	1,R0
	NOP	
; end loop 1-1


	LDIU	*+AR3(2),R0
	PUSH	R0
	CALL	_pfree__FPv


	LDIU	*+AR3(1),R0
	PUSH	R0
	CALL	_pfree__FPv


	LDIU	5000,R0
	PUSH	R0
	CALL	_pmalloc__Fi


	PUSH	R0
	LDP	DEF6,R0
	LSH	16,R0
	OR	DEF6,R0
	PUSH	R0
	CALL	_printf
	SUBI	5,SP


	CALL	_pclear__Fv

	LDIU	14,AR6

	LDIU	97,R0

; begin loop 2-1
L41:	PUSH	R0
	CALL	_pmalloc__Fi
	LDIU	R0,AR7


	PUSH	R0
	LDP	DEF3,R1
	LSH	16,R1
	OR	DEF3,R1
	PUSH	R1
	CALL	_printf

	DBUD	AR6,L41
	SUBI	3,SP
	LDIU	97,R0
	NOP	
; end loop 2-1

	LDIU	0,R0

	LDIU	*+AR3(12),AR7
	LDIU	*+AR3(11),AR6
	LDIU	*-AR3(1),BK
	BUD	BK
	LDIU	AR3,SP
	LDIU	*AR3,AR3
	SUBI	2,SP

;
	.def	_CWInitmain$C
	.def	_CWInitmain$C$LAJ
;=============================================================================
; LOCAL OPTIONS								     |
;   Calling convention     : Tartan register parameters			     |
;									     |
; GENERATED CODE PROPERTIES						     |
;   Total words of code    : 3						     |
;   Volatile registers used: R0,DP,BK					     |
;   Stack frame            : quick (AR3 points to some old frame)	     |
;=============================================================================

	.sect	"T:RTN2"

_CWInitmain$C:
	POP	BK
_CWInitmain$C$LAJ:

	LDIU	0,R0

	BU	BK

;


	.sect	".init"
DEF3:	.byte	37,120	; %x
	.byte	10,0	;   

DEF6:	.byte	37,120	; %x
	.byte	32,46	;  .
	.byte	46,46	; ..
	.byte	46,10	; . 
	.byte	0	;  

;=============================================================================
; Total words of code = 67						     |
; Total words of data = 13						     |
;=============================================================================

	.end

CPS_END_NAMESPACE
