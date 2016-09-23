#include<config.h>
CPS_START_NAMESPACE
**--------------------------------------------------------------------
**  CVS keywords
**
**  $Source: /space/cvs/cps/cps++/tests/various/overlay/test.asm,v $
**  $State: Exp $
**
**--------------------------------------------------------------------
	.version	30

	.def	_test
	.def	tmpxxx

	.ref 	_out_buf

FP	.set	AR3

	.sect	"T:wfm1"
tmpxxx	.word	tmp1
tmp1	.word	1


	.sect	"T:wfm1"
_test:

*---------------------------------------------------------------------------------------
	PUSH	FP
	LDI	SP, FP
*  Save all registers that are important to C
	PUSH    R4
	PUSH    R5
	PUSHF   R6
	PUSHF   R7
	PUSH    AR4
	PUSH    AR5
	PUSH    AR6
	PUSH    AR7
	PUSH    FP              ; Local frame pointer
	PUSH    DP
	PUSH	ST		; Save the status register

	LDI	1010,	R0

	LDP	_out_buf,	AR7
	LSH	16,AR7
	OR	_out_buf,	AR7

	LDI	444h,		R0
	STI	R0,		*AR7++
	LDI	333h,		R0
	STI	R0,		*AR7++

	LDP	@tmpxxx
	LDI	@tmpxxx,	R0
	STI	R0,		*AR7++


*---------------------------------------------------------------------------------------
*  Restore the registers before returning to the C program                
*---------------------------------------------------------------------------------------
	POP	ST		; Restore the status register
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
	RETS
	.end



CPS_END_NAMESPACE
