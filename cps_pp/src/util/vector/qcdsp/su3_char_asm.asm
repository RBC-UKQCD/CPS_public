**--------------------------------------------------------------------
**  CVS keywords
**
**  $Author: mcneile $
**  $Date: 2003-06-22 13:34:46 $
**  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/vector/qcdsp/su3_char_asm.asm,v 1.1.1.1 2003-06-22 13:34:46 mcneile Exp $
**  $Id: su3_char_asm.asm,v 1.1.1.1 2003-06-22 13:34:46 mcneile Exp $
**  $Name: not supported by cvs2svn $
**  $Locker:  $
**  $Log: not supported by cvs2svn $
**  Revision 1.2  2001/06/19 18:13:39  anj
**  Serious ANSIfication.  Plus, degenerate double64.h files removed.
**  Next version will contain the new nga/include/double64.h.  Also,
**  Makefile.gnutests has been modified to work properly, propagating the
**  choice of C++ compiler and flags all the way down the directory tree.
**  The mpi_scu code has been added under phys/nga, and partially
**  plumbed in.
**
**  Everything has newer dates, due to the way in which this first alteration was handled.
**
**  Anj.
**
**  Revision 1.2  2001/05/25 06:16:11  cvs
**  Added CVS keywords to phys_v4_0_0_preCVS
**
**  $RCSfile: su3_char_asm.asm,v $
**  $Revision: 1.1.1.1 $
**  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/vector/qcdsp/su3_char_asm.asm,v $
**  $State: Exp $
**
**--------------------------------------------------------------------
	.version 30
;=============================================================================
; This file is   : su3_char.asm						     |
; generated from : su3_char						     |
; using          : Tartan C/C++ Compiler for the TMS320C3x/4x, v2.1.1	     |
; on (M/D/Y H:M) : 02/24/1998 15:42					     |
;----------------------------------------------------------------------------|
;									     |
; GLOBAL ATTRIBUTES AND ASSERTIONS					     |
;   Wait states            : Code = 0, Stack = 0, Data page = 0, Heap = 0    |
;   Max loop iterations    : 2**23 (no huge loops)			     |
;   Code Span              : 24-bits of address				     |
;									     |
; GLOBAL OPTIONS							     |
;   Optimization heuristics: space					     |
;   Debug information      : on						     |
;   Profiling              : off					     |
;   Max RPTS allowed       : 2**31					     |
;   Data Page              : No_Page					     |
;   Data Page is in ROM    : false					     |
;   *ARn same mem as *ARm? : true					     |
;									     |
;=============================================================================

;
	.def	_reChar6__FPf
	.def	_reChar6__FPf$LAJ
;=============================================================================
; LOCAL OPTIONS								     |
;   Calling convention     : TI C stack parameters			     |
;									     |
; GENERATED CODE PROPERTIES						     |
;   Total words of code    : 53						     |
;   Volatile registers used: R0,R1,R2,R3,R4,R5,R6,AR0,AR1,AR2,DP,BK,RE,RC    |
;   Parameters             : AR0	holds p				     |
;   Stack frame            : full (frame pointer in AR3)		     |
;=============================================================================

	.sect	"T:RTN1"

_reChar6__FPf:
	POP	BK
_reChar6__FPf$LAJ:
	ADDI	1,SP
				; 2 cycle register conflict pipeline delay
	PUSH	AR3
	LDIU	SP,AR3
	LDIU	R5,RC
	LDIU	R4,RE
	PUSHF	R6
	LDIU	*-AR3(2),AR0
				; end of routine prologue

	LDIU	AR0,AR2		; AR0 -> AR2
	LDIU	AR0,AR1		; AR0 -> AR1
				; 2 cycle register conflict pipeline delay

				; source line 25
	LDFU	*+AR0(8),R1	; *+AR0(8) -> R1
	ADDF3	R1,*AR0,R2	; indirect_thru(p) + R1 -> R2
	LDFU	*+AR0(16),R3	; *+AR0(16) -> R3
	ADDF	R3,R2		; R2 + R3 -> R2
	RND	R2,R2		; GRF
	MPYF	*AR0,R2		; R2 * indirect_thru(p) -> R2
	LDFU	*+AR0(9),R4	; *+AR0(9) -> R4
	ADDF3	R4,*+AR0(1),R5	; *+AR0(1) + R4 -> R5
	LDFU	*+AR0(17),R6	; *+AR0(17) -> R6
	ADDF	R6,R5		; R5 + R6 -> R5
	RND	R5,R5		; GRF
	MPYF	*++AR1(1),R5	; R5 * *++AR1(1) -> R5
	SUBF	R5,R2		; R2 - R5 -> R2
	ADDF3	R3,R1,R5	; R1 + R3 -> R5
	RND	R5,R5		; GRF
	MPYF	R5,R1		; R1 * R5 -> R1
	ADDF	R2,R1		; R1 + R2 -> R1
	ADDF3	R6,R4,R2	; R4 + R6 -> R2
	RND	R2,R2		; GRF
	MPYF	R4,R2		; R2 * R4 -> R2
	SUBF	R2,R1		; R1 - R2 -> R1
	MPYF	R3,R3		; R3 * R3 -> R3
	ADDF	R3,R1		; R1 + R3 -> R1
	ADDI	6,AR2		; AR2 + 6 -> AR2
	MPYF	R6,R6		; R6 * R6 -> R6
	SUBF	R6,R1		; R1 - R6 -> R1
	MPYF3	*AR2++(1),*++AR1(1),R2
				; *++AR1(1) * *AR2++(1) -> R2
	ADDF	R2,R1		; R1 + R2 -> R1
	MPYF3	*AR2,*++AR1(1),R2
				; *++AR1(1) * *AR2 -> R2
	ADDI	5,AR2		; AR2 + 5 -> AR2
	SUBF	R2,R1		; R1 - R2 -> R1
				; 1 cycle register conflict pipeline delay
	MPYF3	*AR2++(1),*++AR1(1),R2
				; *++AR1(1) * *AR2++(1) -> R2
	ADDF	R2,R1		; R1 + R2 -> R1
	MPYF3	*AR2++(1),*++AR1(1),R2
				; *++AR1(1) * *AR2++(1) -> R2
	ADDI	5,AR1		; AR1 + 5 -> AR1
	SUBF	R2,R1		; R1 - R2 -> R1
				; 1 cycle register conflict pipeline delay
	MPYF3	*AR2++(1),*AR1++(1),R2
				; *AR1++(1) * *AR2++(1) -> R2
	ADDF	R2,R1		; R1 + R2 -> R1
	MPYF3	*AR2,*AR1,R2	; *AR1 * *AR2 -> R2
	SUBF3	R2,R1,R0	; R1 - R2 -> my_return_place

				; source line 13

				; begin of routine epilogue
	LDFU	*+AR3(1),R6
	LDIU	RE,R4
	LDIU	RC,R5
	BUD	BK
	LDIU	AR3,SP
	LDIU	*AR3,AR3
	SUBI	2,SP

;
	.def	_imChar6__FPf
	.def	_imChar6__FPf$LAJ
;=============================================================================
; LOCAL OPTIONS								     |
;   Calling convention     : TI C stack parameters			     |
;									     |
; GENERATED CODE PROPERTIES						     |
;   Total words of code    : 51						     |
;   Volatile registers used: R0,R1,R2,R3,R4,R5,R6,AR0,AR1,AR2,DP,BK,RE,RC    |
;   Parameters             : AR0	holds p				     |
;   Stack frame            : full (frame pointer in AR3)		     |
;=============================================================================

	.sect	"T:RTN2"

_imChar6__FPf:
	POP	BK
_imChar6__FPf$LAJ:
	ADDI	1,SP
				; 2 cycle register conflict pipeline delay
	PUSH	AR3
	LDIU	SP,AR3
	LDIU	R5,RC
	LDIU	R4,RE
	PUSHF	R6
	LDIU	*-AR3(2),AR0
				; end of routine prologue

	LDIU	AR0,AR2		; AR0 -> AR2
	LDIU	AR0,AR1		; AR0 -> AR1
				; 2 cycle register conflict pipeline delay

				; source line 48
	LDFU	*+AR0(9),R1	; *+AR0(9) -> R1
	ADDF3	R1,*+AR0(1),R2	; *+AR0(1) + R1 -> R2
	LDFU	*+AR0(17),R3	; *+AR0(17) -> R3
	ADDF	R3,R2		; R2 + R3 -> R2
	RND	R2,R2		; GRF
	MPYF	*AR0,R2		; R2 * indirect_thru(p) -> R2
	LDFU	*++AR1(8),R4	; *++AR1(8) -> R4
	ADDF3	R4,*AR0,R5	; indirect_thru(p) + R4 -> R5
	LDFU	*+AR0(16),R6	; *+AR0(16) -> R6
	ADDF	R6,R5		; R5 + R6 -> R5
	RND	R5,R5		; GRF
	MPYF	*++AR2(1),R5	; R5 * *++AR2(1) -> R5
	ADDF	R5,R2		; R2 + R5 -> R2
	ADDF3	R3,R1,R5	; R1 + R3 -> R5
	RND	R5,R5		; GRF
	MPYF	R4,R5		; R5 * R4 -> R5
	ADDF	R5,R2		; R2 + R5 -> R2
	ADDF	R6,R4		; R4 + R6 -> R4
	RND	R4,R4		; GRF
	MPYF	R4,R1		; R1 * R4 -> R1
	ADDF	R2,R1		; R1 + R2 -> R1
	MPYF	R6,R3		; R3 * R6 -> R3
	ADDF	R3,R1		; R1 + R3 -> R1
	ADDF	R3,R1		; R1 + R3 -> R1
	MPYF3	*--AR1(1),*++AR2(1),R2
				; *++AR2(1) * *--AR1(1) -> R2
	ADDF	R2,R1		; R1 + R2 -> R1
	MPYF3	*--AR1(1),*++AR2(1),R2
				; *++AR2(1) * *--AR1(1) -> R2
	ADDI	7,AR1		; AR1 + 7 -> AR1
	ADDF	R2,R1		; R1 + R2 -> R1
				; 1 cycle register conflict pipeline delay
	MPYF3	*AR1--(1),*++AR2(1),R2
				; *++AR2(1) * *AR1--(1) -> R2
	ADDF	R2,R1		; R1 + R2 -> R1
	MPYF3	*AR1--(1),*++AR2(1),R2
				; *++AR2(1) * *AR1--(1) -> R2
	ADDI	10,AR2		; AR2 + 10 -> AR2
	ADDF	R2,R1		; R1 + R2 -> R1
				; 1 cycle register conflict pipeline delay
	MPYF3	*AR2--(1),*-AR1(1),R2
				; *-AR1(1) * *AR2--(1) -> R2
	ADDF	R2,R1		; R1 + R2 -> R1
	MPYF3	*AR2,*AR1,R2	; *AR1 * *AR2 -> R2
	ADDF3	R2,R1,R0	; R1 + R2 -> my_return_place

				; source line 36

				; begin of routine epilogue
	LDFU	*+AR3(1),R6
	LDIU	RE,R4
	LDIU	RC,R5
	BUD	BK
	LDIU	AR3,SP
	LDIU	*AR3,AR3
	SUBI	2,SP

;
	.def	_reChar8__FPf
	.def	_reChar8__FPf$LAJ
;=============================================================================
; LOCAL OPTIONS								     |
;   Calling convention     : TI C stack parameters			     |
;									     |
; GENERATED CODE PROPERTIES						     |
;   Total words of code    : 73						     |
;   Volatile registers used: R0,R1,R2,R3,R4,R5,R6,R7,AR0,AR1,AR2,DP,BK,RE,RC |
;   Parameters             : AR0	holds p				     |
;   Stack frame            : full (frame pointer in AR3)		     |
;=============================================================================

	.sect	"T:RTN3"

_reChar8__FPf:
	POP	BK
_reChar8__FPf$LAJ:
	ADDI	1,SP
				; 2 cycle register conflict pipeline delay
	PUSH	AR3
	LDIU	SP,AR3
	LDIU	R5,RC
	LDIU	R4,RE
	PUSHF	R6
	PUSHF	R7
	LDIU	*-AR3(2),AR0
				; end of routine prologue

	LDIU	AR0,AR2		; AR0 -> AR2
	LDIU	AR0,AR1		; AR0 -> AR1
	LDP	f2o3,DP		; f2o3 -> DP
				; 1 cycle register conflict pipeline delay

				; source line 60
	LDFU	*++AR1(8),R1	; *++AR1(8) -> R1
	LDFU	*++AR2(16),R2	; *++AR2(16) -> R2
	ADDF3	R2,R1,R3	; R1 + R2 -> R3
	RND	R3,R3		; GRF
	MPYF	*AR0,R3		; R3 * indirect_thru(p) -> R3
	LDF	*++AR1(1),R4	; *++AR1(1) -> R4
				; in parallel with
     || LDF	*++AR2(1),R5	; *++AR2(1) -> R5
	ADDF3	R5,R4,R6	; R4 + R5 -> R6
	RND	R6,R6		; GRF
	MPYF	*--AR1(8),R6	; R6 * *--AR1(8) -> R6
	ADDF	R6,R3		; R3 + R6 -> R3
	MPYF3	R2,R1,R6	; R1 * R2 -> R6
	ADDF	R6,R3		; R3 + R6 -> R3
	MPYF3	R5,R4,R6	; R4 * R5 -> R6
	ADDF	R6,R3		; R3 + R6 -> R3
	ADDF	R3,R3		; R3 + R3 -> R3
	MPYF3	*+AR0(1),*+AR0(1),R6
				; *+AR0(1) * *+AR0(1) -> R6
	MPYF3	*AR0,*AR0,R7	; indirect_thru(p) * indirect_thru(p) -> R7
	ADDF	R7,R6		; R6 + R7 -> R6
	MPYF	R1,R1		; R1 * R1 -> R1
	ADDF	R6,R1		; R1 + R6 -> R1
	MPYF	R4,R4		; R4 * R4 -> R4
	ADDF	R4,R1		; R1 + R4 -> R1
	MPYF	R2,R2		; R2 * R2 -> R2
	ADDF	R2,R1		; R1 + R2 -> R1
	MPYF	R5,R5		; R5 * R5 -> R5
	ADDF	R5,R1		; R1 + R5 -> R1
	RND	R1,R1		; GRF
	MPYF	@f2o3,R1	; R1 * f2o3 -> R1
	ADDF	R3,R1		; R1 + R3 -> R1
	MPYF3	*+AR1(1),*++AR1(1),R2
				; *++AR1(1) * *+AR1(1) -> R2
	MPYF3	*+AR1(1),*++AR1(1),R3
				; *++AR1(1) * *+AR1(1) -> R3
	ADDF	R3,R2		; R2 + R3 -> R2
	MPYF3	*+AR1(1),*++AR1(1),R3
				; *++AR1(1) * *+AR1(1) -> R3
	ADDF	R3,R2		; R2 + R3 -> R2
	MPYF3	*+AR1(1),*++AR1(1),R3
				; *++AR1(1) * *+AR1(1) -> R3
	ADDF	R3,R2		; R2 + R3 -> R2
	MPYF3	*+AR1(1),*++AR1(1),R3
				; *++AR1(1) * *+AR1(1) -> R3
	ADDF	R3,R2		; R2 + R3 -> R2
	MPYF3	*+AR1(1),*++AR1(1),R3
				; *++AR1(1) * *+AR1(1) -> R3
	ADDI	3,AR1		; AR1 + 3 -> AR1
	LDP	f1o3,DP		; f1o3 -> DP
	ADDF	R3,R2		; R2 + R3 -> R2
	MPYF3	*AR1,*AR1++(1),R3
				; *AR1++(1) * *AR1 -> R3
	ADDF	R3,R2		; R2 + R3 -> R2
	MPYF3	*AR1,*AR1++(1),R3
				; *AR1++(1) * *AR1 -> R3
	ADDF	R3,R2		; R2 + R3 -> R2
	MPYF3	*AR1,*AR1++(1),R3
				; *AR1++(1) * *AR1 -> R3
	ADDF	R3,R2		; R2 + R3 -> R2
	MPYF3	*AR1,*AR1++(1),R3
				; *AR1++(1) * *AR1 -> R3
	ADDF	R3,R2		; R2 + R3 -> R2
	MPYF3	*AR1,*AR1++(1),R3
				; *AR1++(1) * *AR1 -> R3
	ADDF	R3,R2		; R2 + R3 -> R2
	MPYF3	*AR1,*AR1,R3	; *AR1 * *AR1 -> R3
	ADDF	R3,R2		; R2 + R3 -> R2
	RND	R2,R2		; GRF
	MPYF	@f1o3,R2	; R2 * f1o3 -> R2
	SUBF3	R2,R1,R0	; R1 - R2 -> my_return_place

				; source line 58

				; begin of routine epilogue
	LDFU	*+AR3(2),R7
	LDFU	*+AR3(1),R6
	LDIU	RE,R4
	LDIU	RC,R5
	BUD	BK
	LDIU	AR3,SP
	LDIU	*AR3,AR3
	SUBI	2,SP

;
	.def	_reChar10__FPf
	.def	_reChar10__FPf$LAJ
;=============================================================================
; LOCAL OPTIONS								     |
;   Calling convention     : TI C stack parameters			     |
;									     |
; GENERATED CODE PROPERTIES						     |
;   Total words of code    : 3						     |
;   Volatile registers used: R0,DP,BK					     |
;   Parameters             : Never Used holds p				     |
;   Stack frame            : quick (AR3 points to some old frame)	     |
;=============================================================================

	.sect	"T:RTN4"

_reChar10__FPf$LAJ:
	PUSH	BK
_reChar10__FPf:
				; end of routine prologue


				; source line 76
	LDFU	0.0,R0		; 0.0 -> my_return_place

				; source line 74

				; begin of routine epilogue
	RETSU	

;
	.def	_imChar10__FPf
	.def	_imChar10__FPf$LAJ
;=============================================================================
; LOCAL OPTIONS								     |
;   Calling convention     : TI C stack parameters			     |
;									     |
; GENERATED CODE PROPERTIES						     |
;   Total words of code    : 3						     |
;   Volatile registers used: R0,DP,BK					     |
;   Parameters             : Never Used holds p				     |
;   Stack frame            : quick (AR3 points to some old frame)	     |
;=============================================================================

	.sect	"T:RTN5"

_imChar10__FPf$LAJ:
	PUSH	BK
_imChar10__FPf:
				; end of routine prologue


				; source line 81
	LDFU	0.0,R0		; 0.0 -> my_return_place

				; source line 79

				; begin of routine epilogue
	RETSU	

;
	.def	_CWInitsu3_char$C
	.def	_CWInitsu3_char$C$LAJ
;=============================================================================
; LOCAL OPTIONS								     |
;   Calling convention     : Tartan register parameters			     |
;									     |
; GENERATED CODE PROPERTIES						     |
;   Total words of code    : 3						     |
;   Volatile registers used: R0,DP,BK					     |
;   Stack frame            : quick (AR3 points to some old frame)	     |
;=============================================================================

	.sect	"T:RTN6"

_CWInitsu3_char$C:
	POP	BK
_CWInitsu3_char$C$LAJ:
				; end of routine prologue

	LDIU	0,R0		; 0 -> my_return_place

				; begin of routine epilogue
	BU	BK

;

	.data
OWN:	.word	0ff2aaaabh	;0.66666666666666666666

	.word	0fe2aaaabh	;0.33333333333333333333


f1o3	.set	OWN+1

f2o3	.set	OWN






;=============================================================================
; Total words of code = 186						     |
; Total words of data = 2						     |
;=============================================================================

	.end

