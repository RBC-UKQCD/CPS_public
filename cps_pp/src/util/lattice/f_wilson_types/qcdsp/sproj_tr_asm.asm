**--------------------------------------------------------------------
**  CVS keywords
**
**  $Author: chulwoo $
**  $Date: 2004-06-04 21:14:13 $
**  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/lattice/f_wilson_types/qcdsp/sproj_tr_asm.asm,v 1.3 2004-06-04 21:14:13 chulwoo Exp $
**  $Id: sproj_tr_asm.asm,v 1.3 2004-06-04 21:14:13 chulwoo Exp $
**  $Name: not supported by cvs2svn $
**  $Locker:  $
**  $Revision: 1.3 $
**  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/lattice/f_wilson_types/qcdsp/sproj_tr_asm.asm,v $
**  $State: Exp $
**
**--------------------------------------------------------------------
	.version 30
;=============================================================================
; This file is   : sproj_tr.asm						     |
; generated from : sproj_tr						     |
; using          : Tartan C/C++ Compiler for the TMS320C3x/4x, v2.1.1	     |
; on (M/D/Y H:M) : 01/22/1998 14:13					     |
;----------------------------------------------------------------------------|
;									     |
; GLOBAL ATTRIBUTES AND ASSERTIONS					     |
;   Wait states            : Code = 0, Stack = 0, Data page = 0, Heap = 0    |
;   Max loop iterations    : 2**23 (no huge loops)			     |
;   Code Span              : 24-bits of address				     |
;									     |
; GLOBAL OPTIONS							     |
;   Optimization heuristics: space					     |
;   Debug information      : off					     |
;   Profiling              : off					     |
;   Max RPTS allowed       : 2**31					     |
;   Data Page              : No_Page					     |
;   Data Page is in ROM    : false					     |
;   *ARn same mem as *ARm? : true					     |
;									     |
;=============================================================================

;
	.def	_sprojTrXp__FPfN21iN24
	.def	_sprojTrXp__FPfN21iN24$LAJ
;=============================================================================
; LOCAL OPTIONS								     |
;   Calling convention     : TI C stack parameters			     |
;									     |
; GENERATED CODE PROPERTIES						     |
;   Total words of code    : 123					     |
;   Volatile registers used: R0,R1,R2,R3,R4,R5,R6,R7,AR0,AR1,AR2,DP,IR0,IR1, |
;			     BK,RS,RE,RC				     |
;   Registers for locals   : AR0	holds tw			     |
;			     AR5	holds tv			     |
;			     AR6	holds tf			     |
;   Parameters             : *-AR3(2)   holds f				     |
;			     *-AR3(3)   holds v				     |
;			     *-AR3(4)   holds w				     |
;			     *-AR3(5)   holds num_blk			     |
;			     *-AR3(6)   holds v_stride			     |
;			     *-AR3(7)   holds w_stride			     |
;   Stack frame            : full (frame pointer in AR3)		     |
;=============================================================================

	.sect	"T:RTN1"

_sprojTrXp__FPfN21iN24$LAJ:
	PUSH	BK
_sprojTrXp__FPfN21iN24:
	PUSH	AR3
	LDIU	SP,AR3
	ADDI	3,SP
	PUSH	R4
	PUSH	R5
	RND	R6, R6		; GRF
	PUSHF	R6
	RND	R7, R7		; GRF
	PUSHF	R7
	PUSH	AR5
	PUSH	AR6

	LDIU	17,RC
	LDP	DEF1,DP
	LDIU	*-AR3(2),IR1
	LDIU	IR1,AR6
	.word	40608000h	;	LDFU	0.0,R0
	RPTB	L116

; begin loop 1-1
L116:	STF	R0,*AR6++(1)
; end loop 1-1

	LDIU	*-AR3(5),R0
	CMPI	1,R0
	BLT	L4
	BR	L115

; begin loop 2-1
L3:	LDIU	IR1,AR6
	LDIU	*-AR3(3),AR5
	BRD	L114
	LDIU	2,R1
	LDIU	*-AR3(4),R0
	LDIU	R0,AR0

; begin loop 3-2
L5:	LDIU	2,RC
	RPTB	L117

; begin loop 4-3
	LDIU	AR5,AR2
	LDIU	AR6,AR1
	ADDI	1,AR6
	LDFU	*+AR5(19),R1
	SUBF3	R1,*AR5,R2
	RND	R2, R2		; GRF
	MPYF	*AR0,R2
	LDFU	*+AR5(18),R3
	ADDF3	R3,*+AR5(1),R4
	RND	R4, R4		; GRF
	MPYF	*+AR0(1),R4
	ADDF	R4,R2
	LDFU	*+AR5(6),R4
	LDFU	*++AR2(13),R5
	SUBF3	R5,R4,R6
	RND	R6, R6		; GRF
	MPYF	*+AR0(6),R6
	ADDF	R6,R2
	LDFU	*+AR5(7),R6
	ADDF3	*-AR2(1),R6,R7
	RND	R7, R7		; GRF
	MPYF	*+AR0(7),R7
	ADDF	R7,R2
	ADDF	*+AR5(12),R6
	RND	R6, R6		; GRF
	MPYF	*+AR0(12),R6
	ADDF	R6,R2
	SUBRF	R5,R4
	RND	R4, R4		; GRF
	MPYF	*+AR0(13),R4
	ADDF	R4,R2
	ADDF	*+AR5(1),R3
	RND	R3, R3		; GRF
	MPYF	*+AR0(18),R3
	ADDF	R3,R2
	SUBF	*AR5,R1
	RND	R1, R1		; GRF
	MPYF	*+AR0(19),R1
	ADDF	R2,R1
	ADDF	*AR1,R1
	RND	R1, R1		; GRF
	STF	R1,*AR1
	LDIU	AR6,AR1
	ADDI	1,AR6
	LDFU	*+AR5(18),R1
	ADDF3	R1,*+AR5(1),R2
	RND	R2, R2		; GRF
	MPYF	*AR0++(1),R2
	LDFU	*+AR5(19),R3
	SUBF3	R3,*AR5,R4
	RND	R4, R4		; GRF
	MPYF	*AR0++(5),R4
	SUBF	R4,R2
	LDFU	*+AR5(7),R4
	LDFU	*+AR5(12),R5
	ADDF3	R5,R4,R6
	RND	R6, R6		; GRF
	MPYF	*AR0++(1),R6
	ADDF	R6,R2
	LDFU	*+AR5(6),R6
	SUBF3	*AR2,R6,R7
	RND	R7, R7		; GRF
	MPYF	*AR0++(5),R7
	SUBF	R7,R2
	SUBRF	*+AR5(13),R6
	RND	R6, R6		; GRF
	MPYF	*AR0++(1),R6
	ADDF	R6,R2
	ADDF	R5,R4
	RND	R4, R4		; GRF
	MPYF	*AR0++(5),R4
	SUBF	R4,R2
	SUBF	*AR5,R3
	RND	R3, R3		; GRF
	MPYF	*AR0++(1),R3
	ADDF	R3,R2
	ADDF	*+AR5(1),R1
	RND	R1, R1		; GRF
	MPYF	*AR0--(17),R1
	SUBRF	R2,R1
	ADDF	*AR1,R1
	RND	R1, R1		; GRF
L117:	STF	R1,*AR1
; end loop 4-3

	LDIU	R0,AR0
	ADDI	2,AR5
	LDIU	*+AR3(2),R1
	SUBI	1,R1
L114:	STI	R1,*+AR3(2)
	CMPI	0,R1
	BGE	L5
; end loop 3-2

	LDIU	*-AR3(3),IR0
	ADDI	*-AR3(6),IR0
	ADDI	24,IR0
	STI	IR0,*-AR3(3)
	ADDI	*-AR3(7),R0
	ADDI	24,R0
	STI	R0,*-AR3(4)
	LDIU	*+AR3(1),R0
L115:	SUBI	1,R0
	STI	R0,*+AR3(1)
	CMPI	0,R0
	BGE	L3
; end loop 2-1

L4:	LDIU	*-AR3(5),R0

	LDIU	*+AR3(9),AR6
	LDIU	*+AR3(8),AR5
	LDFU	*+AR3(7),R7
	LDFU	*+AR3(6),R6
	LDIU	*+AR3(5),R5
	LDIU	*+AR3(4),R4
	LDIU	AR3,SP
	POP	AR3
	RETSU	

;
	.def	_sprojTrXm__FPfN21iN24
	.def	_sprojTrXm__FPfN21iN24$LAJ
;=============================================================================
; LOCAL OPTIONS								     |
;   Calling convention     : TI C stack parameters			     |
;									     |
; GENERATED CODE PROPERTIES						     |
;   Total words of code    : 123					     |
;   Volatile registers used: R0,R1,R2,R3,R4,R5,R6,R7,AR0,AR1,AR2,DP,IR0,IR1, |
;			     BK,RS,RE,RC				     |
;   Registers for locals   : AR0	holds tw			     |
;			     AR5	holds tv			     |
;			     AR6	holds tf			     |
;   Parameters             : *-AR3(2)   holds f				     |
;			     *-AR3(3)   holds v				     |
;			     *-AR3(4)   holds w				     |
;			     *-AR3(5)   holds num_blk			     |
;			     *-AR3(6)   holds v_stride			     |
;			     *-AR3(7)   holds w_stride			     |
;   Stack frame            : full (frame pointer in AR3)		     |
;=============================================================================

	.sect	"T:RTN2"

_sprojTrXm__FPfN21iN24$LAJ:
	PUSH	BK
_sprojTrXm__FPfN21iN24:
	PUSH	AR3
	LDIU	SP,AR3
	ADDI	3,SP
	PUSH	R4
	PUSH	R5
	RND	R6, R6		; GRF
	PUSHF	R6
	RND	R7, R7		; GRF
	PUSHF	R7
	PUSH	AR5
	PUSH	AR6

	LDIU	17,RC
	LDP	DEF1,DP
	LDIU	*-AR3(2),IR1
	LDIU	IR1,AR6
	.word	40608000h	;	LDFU	0.0,R0
	RPTB	L154

; begin loop 5-1
L154:	STF	R0,*AR6++(1)
; end loop 5-1

	LDIU	*-AR3(5),R0
	CMPI	1,R0
	BLT	L14
	BR	L153

; begin loop 6-1
L13:	LDIU	IR1,AR6
	LDIU	*-AR3(3),AR5
	BRD	L152
	LDIU	2,R1
	LDIU	*-AR3(4),R0
	LDIU	R0,AR0

; begin loop 7-2
L15:	LDIU	2,RC
	RPTB	L155

; begin loop 8-3
	LDIU	AR5,AR2
	LDIU	AR6,AR1
	ADDI	1,AR6
	LDFU	*+AR5(19),R1
	ADDF3	R1,*AR5,R2
	RND	R2, R2		; GRF
	MPYF	*AR0,R2
	LDFU	*+AR5(18),R3
	SUBF3	R3,*+AR5(1),R4
	RND	R4, R4		; GRF
	MPYF	*+AR0(1),R4
	ADDF	R4,R2
	LDFU	*+AR5(6),R4
	LDFU	*++AR2(13),R5
	ADDF3	R5,R4,R6
	RND	R6, R6		; GRF
	MPYF	*+AR0(6),R6
	ADDF	R6,R2
	LDFU	*+AR5(7),R6
	SUBF3	*-AR2(1),R6,R7
	RND	R7, R7		; GRF
	MPYF	*+AR0(7),R7
	ADDF	R7,R2
	SUBRF	*+AR5(12),R6
	RND	R6, R6		; GRF
	MPYF	*+AR0(12),R6
	ADDF	R6,R2
	ADDF	R5,R4
	RND	R4, R4		; GRF
	MPYF	*+AR0(13),R4
	ADDF	R4,R2
	SUBF	*+AR5(1),R3
	RND	R3, R3		; GRF
	MPYF	*+AR0(18),R3
	ADDF	R3,R2
	ADDF	*AR5,R1
	RND	R1, R1		; GRF
	MPYF	*+AR0(19),R1
	ADDF	R2,R1
	ADDF	*AR1,R1
	RND	R1, R1		; GRF
	STF	R1,*AR1
	LDIU	AR6,AR1
	ADDI	1,AR6
	LDFU	*+AR5(18),R1
	SUBF3	R1,*+AR5(1),R2
	RND	R2, R2		; GRF
	MPYF	*AR0++(1),R2
	LDFU	*+AR5(19),R3
	ADDF3	R3,*AR5,R4
	RND	R4, R4		; GRF
	MPYF	*AR0++(5),R4
	SUBF	R4,R2
	LDFU	*+AR5(7),R4
	LDFU	*+AR5(12),R5
	SUBF3	R5,R4,R6
	RND	R6, R6		; GRF
	MPYF	*AR0++(1),R6
	ADDF	R6,R2
	LDFU	*+AR5(6),R6
	ADDF3	*AR2,R6,R7
	RND	R7, R7		; GRF
	MPYF	*AR0++(5),R7
	SUBF	R7,R2
	ADDF	*+AR5(13),R6
	RND	R6, R6		; GRF
	MPYF	*AR0++(1),R6
	ADDF	R6,R2
	SUBRF	R5,R4
	RND	R4, R4		; GRF
	MPYF	*AR0++(5),R4
	SUBF	R4,R2
	ADDF	*AR5,R3
	RND	R3, R3		; GRF
	MPYF	*AR0++(1),R3
	ADDF	R3,R2
	SUBF	*+AR5(1),R1
	RND	R1, R1		; GRF
	MPYF	*AR0--(17),R1
	SUBRF	R2,R1
	ADDF	*AR1,R1
	RND	R1, R1		; GRF
L155:	STF	R1,*AR1
; end loop 8-3

	LDIU	R0,AR0
	ADDI	2,AR5
	LDIU	*+AR3(2),R1
	SUBI	1,R1
L152:	STI	R1,*+AR3(2)
	CMPI	0,R1
	BGE	L15
; end loop 7-2

	LDIU	*-AR3(3),IR0
	ADDI	*-AR3(6),IR0
	ADDI	24,IR0
	STI	IR0,*-AR3(3)
	ADDI	*-AR3(7),R0
	ADDI	24,R0
	STI	R0,*-AR3(4)
	LDIU	*+AR3(1),R0
L153:	SUBI	1,R0
	STI	R0,*+AR3(1)
	CMPI	0,R0
	BGE	L13
; end loop 6-1

L14:	LDIU	*-AR3(5),R0

	LDIU	*+AR3(9),AR6
	LDIU	*+AR3(8),AR5
	LDFU	*+AR3(7),R7
	LDFU	*+AR3(6),R6
	LDIU	*+AR3(5),R5
	LDIU	*+AR3(4),R4
	LDIU	AR3,SP
	POP	AR3
	RETSU	

;
	.def	_sprojTrYp__FPfN21iN24
	.def	_sprojTrYp__FPfN21iN24$LAJ
;=============================================================================
; LOCAL OPTIONS								     |
;   Calling convention     : TI C stack parameters			     |
;									     |
; GENERATED CODE PROPERTIES						     |
;   Total words of code    : 123					     |
;   Volatile registers used: R0,R1,R2,R3,R4,R5,R6,R7,AR0,AR1,AR2,DP,IR0,IR1, |
;			     BK,RS,RE,RC				     |
;   Registers for locals   : AR0	holds tw			     |
;			     AR5	holds tv			     |
;			     AR6	holds tf			     |
;   Parameters             : *-AR3(2)   holds f				     |
;			     *-AR3(3)   holds v				     |
;			     *-AR3(4)   holds w				     |
;			     *-AR3(5)   holds num_blk			     |
;			     *-AR3(6)   holds v_stride			     |
;			     *-AR3(7)   holds w_stride			     |
;   Stack frame            : full (frame pointer in AR3)		     |
;=============================================================================

	.sect	"T:RTN3"

_sprojTrYp__FPfN21iN24$LAJ:
	PUSH	BK
_sprojTrYp__FPfN21iN24:
	PUSH	AR3
	LDIU	SP,AR3
	ADDI	3,SP
	PUSH	R4
	PUSH	R5
	RND	R6, R6		; GRF
	PUSHF	R6
	RND	R7, R7		; GRF
	PUSHF	R7
	PUSH	AR5
	PUSH	AR6

	LDIU	17,RC
	LDP	DEF1,DP
	LDIU	*-AR3(2),IR1
	LDIU	IR1,AR6
	.word	40608000h	;	LDFU	0.0,R0
	RPTB	L192

; begin loop 9-1
L192:	STF	R0,*AR6++(1)
; end loop 9-1

	LDIU	*-AR3(5),R0
	CMPI	1,R0
	BLT	L24
	BR	L191

; begin loop 10-1
L23:	LDIU	IR1,AR6
	LDIU	*-AR3(3),AR5
	BRD	L190
	LDIU	2,R1
	LDIU	*-AR3(4),R0
	LDIU	R0,AR0

; begin loop 11-2
L25:	LDIU	2,RC
	RPTB	L193

; begin loop 12-3
	LDIU	AR5,AR2
	LDIU	AR6,AR1
	ADDI	1,AR6
	LDFU	*+AR5(18),R1
	SUBF3	R1,*AR5,R2
	RND	R2, R2		; GRF
	MPYF	*AR0,R2
	LDFU	*+AR5(19),R3
	SUBF3	R3,*+AR5(1),R4
	RND	R4, R4		; GRF
	MPYF	*+AR0(1),R4
	ADDF	R4,R2
	LDFU	*+AR5(6),R4
	LDFU	*++AR2(12),R5
	ADDF3	R5,R4,R6
	RND	R6, R6		; GRF
	MPYF	*+AR0(6),R6
	ADDF	R6,R2
	LDFU	*+AR5(7),R6
	ADDF3	*+AR2(1),R6,R7
	RND	R7, R7		; GRF
	MPYF	*+AR0(7),R7
	ADDF	R7,R2
	ADDF	R5,R4
	RND	R4, R4		; GRF
	MPYF	*+AR0(12),R4
	ADDF	R4,R2
	ADDF	*+AR5(13),R6
	RND	R6, R6		; GRF
	MPYF	*+AR0(13),R6
	ADDF	R6,R2
	SUBF	*AR5,R1
	RND	R1, R1		; GRF
	MPYF	*+AR0(18),R1
	ADDF	R2,R1
	SUBF	*+AR5(1),R3
	RND	R3, R3		; GRF
	MPYF	*+AR0(19),R3
	ADDF	R3,R1
	ADDF	*AR1,R1
	RND	R1, R1		; GRF
	STF	R1,*AR1
	LDIU	AR6,AR1
	ADDI	1,AR6
	LDFU	*+AR5(19),R1
	SUBF3	R1,*+AR5(1),R2
	RND	R2, R2		; GRF
	MPYF	*AR0++(1),R2
	LDFU	*+AR5(18),R3
	SUBF3	R3,*AR5,R4
	RND	R4, R4		; GRF
	MPYF	*AR0++(5),R4
	SUBF	R4,R2
	LDFU	*+AR5(7),R4
	LDFU	*+AR5(13),R5
	ADDF3	R5,R4,R6
	RND	R6, R6		; GRF
	MPYF	*AR0++(1),R6
	ADDF	R6,R2
	LDFU	*+AR5(6),R6
	ADDF3	*AR2,R6,R7
	RND	R7, R7		; GRF
	MPYF	*AR0++(5),R7
	SUBF	R7,R2
	ADDF	R5,R4
	RND	R4, R4		; GRF
	MPYF	*AR0++(1),R4
	ADDF	R4,R2
	ADDF	*+AR5(12),R6
	RND	R6, R6		; GRF
	MPYF	*AR0++(5),R6
	SUBF	R6,R2
	SUBF	*+AR5(1),R1
	RND	R1, R1		; GRF
	MPYF	*AR0++(1),R1
	ADDF	R2,R1
	SUBF	*AR5,R3
	RND	R3, R3		; GRF
	MPYF	*AR0--(17),R3
	SUBF	R3,R1
	ADDF	*AR1,R1
	RND	R1, R1		; GRF
L193:	STF	R1,*AR1
; end loop 12-3

	LDIU	R0,AR0
	ADDI	2,AR5
	LDIU	*+AR3(2),R1
	SUBI	1,R1
L190:	STI	R1,*+AR3(2)
	CMPI	0,R1
	BGE	L25
; end loop 11-2

	LDIU	*-AR3(3),IR0
	ADDI	*-AR3(6),IR0
	ADDI	24,IR0
	STI	IR0,*-AR3(3)
	ADDI	*-AR3(7),R0
	ADDI	24,R0
	STI	R0,*-AR3(4)
	LDIU	*+AR3(1),R0
L191:	SUBI	1,R0
	STI	R0,*+AR3(1)
	CMPI	0,R0
	BGE	L23
; end loop 10-1

L24:	LDIU	*-AR3(5),R0

	LDIU	*+AR3(9),AR6
	LDIU	*+AR3(8),AR5
	LDFU	*+AR3(7),R7
	LDFU	*+AR3(6),R6
	LDIU	*+AR3(5),R5
	LDIU	*+AR3(4),R4
	LDIU	AR3,SP
	POP	AR3
	RETSU	

;
	.def	_sprojTrYm__FPfN21iN24
	.def	_sprojTrYm__FPfN21iN24$LAJ
;=============================================================================
; LOCAL OPTIONS								     |
;   Calling convention     : TI C stack parameters			     |
;									     |
; GENERATED CODE PROPERTIES						     |
;   Total words of code    : 123					     |
;   Volatile registers used: R0,R1,R2,R3,R4,R5,R6,R7,AR0,AR1,AR2,DP,IR0,IR1, |
;			     BK,RS,RE,RC				     |
;   Registers for locals   : AR0	holds tw			     |
;			     AR5	holds tv			     |
;			     AR6	holds tf			     |
;   Parameters             : *-AR3(2)   holds f				     |
;			     *-AR3(3)   holds v				     |
;			     *-AR3(4)   holds w				     |
;			     *-AR3(5)   holds num_blk			     |
;			     *-AR3(6)   holds v_stride			     |
;			     *-AR3(7)   holds w_stride			     |
;   Stack frame            : full (frame pointer in AR3)		     |
;=============================================================================

	.sect	"T:RTN4"

_sprojTrYm__FPfN21iN24$LAJ:
	PUSH	BK
_sprojTrYm__FPfN21iN24:
	PUSH	AR3
	LDIU	SP,AR3
	ADDI	3,SP
	PUSH	R4
	PUSH	R5
	RND	R6, R6		; GRF
	PUSHF	R6
	RND	R7, R7		; GRF
	PUSHF	R7
	PUSH	AR5
	PUSH	AR6

	LDIU	17,RC
	LDP	DEF1,DP
	LDIU	*-AR3(2),IR1
	LDIU	IR1,AR6
	.word	40608000h	;	LDFU	0.0,R0
	RPTB	L230

; begin loop 13-1
L230:	STF	R0,*AR6++(1)
; end loop 13-1

	LDIU	*-AR3(5),R0
	CMPI	1,R0
	BLT	L34
	BR	L229

; begin loop 14-1
L33:	LDIU	IR1,AR6
	LDIU	*-AR3(3),AR5
	BRD	L228
	LDIU	2,R1
	LDIU	*-AR3(4),R0
	LDIU	R0,AR0

; begin loop 15-2
L35:	LDIU	2,RC
	RPTB	L231

; begin loop 16-3
	LDIU	AR5,AR2
	LDIU	AR6,AR1
	ADDI	1,AR6
	LDFU	*+AR5(18),R1
	ADDF3	R1,*AR5,R2
	RND	R2, R2		; GRF
	MPYF	*AR0,R2
	LDFU	*+AR5(19),R3
	ADDF3	R3,*+AR5(1),R4
	RND	R4, R4		; GRF
	MPYF	*+AR0(1),R4
	ADDF	R4,R2
	LDFU	*+AR5(6),R4
	LDFU	*++AR2(12),R5
	SUBF3	R5,R4,R6
	RND	R6, R6		; GRF
	MPYF	*+AR0(6),R6
	ADDF	R6,R2
	LDFU	*+AR5(7),R6
	SUBF3	*+AR2(1),R6,R7
	RND	R7, R7		; GRF
	MPYF	*+AR0(7),R7
	ADDF	R7,R2
	SUBRF	R5,R4
	RND	R4, R4		; GRF
	MPYF	*+AR0(12),R4
	ADDF	R4,R2
	SUBRF	*+AR5(13),R6
	RND	R6, R6		; GRF
	MPYF	*+AR0(13),R6
	ADDF	R6,R2
	ADDF	*AR5,R1
	RND	R1, R1		; GRF
	MPYF	*+AR0(18),R1
	ADDF	R2,R1
	ADDF	*+AR5(1),R3
	RND	R3, R3		; GRF
	MPYF	*+AR0(19),R3
	ADDF	R3,R1
	ADDF	*AR1,R1
	RND	R1, R1		; GRF
	STF	R1,*AR1
	LDIU	AR6,AR1
	ADDI	1,AR6
	LDFU	*+AR5(19),R1
	ADDF3	R1,*+AR5(1),R2
	RND	R2, R2		; GRF
	MPYF	*AR0++(1),R2
	LDFU	*+AR5(18),R3
	ADDF3	R3,*AR5,R4
	RND	R4, R4		; GRF
	MPYF	*AR0++(5),R4
	SUBF	R4,R2
	LDFU	*+AR5(7),R4
	LDFU	*+AR5(13),R5
	SUBF3	R5,R4,R6
	RND	R6, R6		; GRF
	MPYF	*AR0++(1),R6
	ADDF	R6,R2
	LDFU	*+AR5(6),R6
	SUBF3	*AR2,R6,R7
	RND	R7, R7		; GRF
	MPYF	*AR0++(5),R7
	SUBF	R7,R2
	SUBRF	R5,R4
	RND	R4, R4		; GRF
	MPYF	*AR0++(1),R4
	ADDF	R4,R2
	SUBRF	*+AR5(12),R6
	RND	R6, R6		; GRF
	MPYF	*AR0++(5),R6
	SUBF	R6,R2
	ADDF	*+AR5(1),R1
	RND	R1, R1		; GRF
	MPYF	*AR0++(1),R1
	ADDF	R2,R1
	ADDF	*AR5,R3
	RND	R3, R3		; GRF
	MPYF	*AR0--(17),R3
	SUBF	R3,R1
	ADDF	*AR1,R1
	RND	R1, R1		; GRF
L231:	STF	R1,*AR1
; end loop 16-3

	LDIU	R0,AR0
	ADDI	2,AR5
	LDIU	*+AR3(2),R1
	SUBI	1,R1
L228:	STI	R1,*+AR3(2)
	CMPI	0,R1
	BGE	L35
; end loop 15-2

	LDIU	*-AR3(3),IR0
	ADDI	*-AR3(6),IR0
	ADDI	24,IR0
	STI	IR0,*-AR3(3)
	ADDI	*-AR3(7),R0
	ADDI	24,R0
	STI	R0,*-AR3(4)
	LDIU	*+AR3(1),R0
L229:	SUBI	1,R0
	STI	R0,*+AR3(1)
	CMPI	0,R0
	BGE	L33
; end loop 14-1

L34:	LDIU	*-AR3(5),R0

	LDIU	*+AR3(9),AR6
	LDIU	*+AR3(8),AR5
	LDFU	*+AR3(7),R7
	LDFU	*+AR3(6),R6
	LDIU	*+AR3(5),R5
	LDIU	*+AR3(4),R4
	LDIU	AR3,SP
	POP	AR3
	RETSU	

;
	.def	_sprojTrZp__FPfN21iN24
	.def	_sprojTrZp__FPfN21iN24$LAJ
;=============================================================================
; LOCAL OPTIONS								     |
;   Calling convention     : TI C stack parameters			     |
;									     |
; GENERATED CODE PROPERTIES						     |
;   Total words of code    : 123					     |
;   Volatile registers used: R0,R1,R2,R3,R4,R5,R6,R7,AR0,AR1,AR2,DP,IR0,IR1, |
;			     BK,RS,RE,RC				     |
;   Registers for locals   : AR0	holds tw			     |
;			     AR5	holds tv			     |
;			     AR6	holds tf			     |
;   Parameters             : *-AR3(2)   holds f				     |
;			     *-AR3(3)   holds v				     |
;			     *-AR3(4)   holds w				     |
;			     *-AR3(5)   holds num_blk			     |
;			     *-AR3(6)   holds v_stride			     |
;			     *-AR3(7)   holds w_stride			     |
;   Stack frame            : full (frame pointer in AR3)		     |
;=============================================================================

	.sect	"T:RTN5"

_sprojTrZp__FPfN21iN24$LAJ:
	PUSH	BK
_sprojTrZp__FPfN21iN24:
	PUSH	AR3
	LDIU	SP,AR3
	ADDI	3,SP
	PUSH	R4
	PUSH	R5
	RND	R6, R6		; GRF
	PUSHF	R6
	RND	R7, R7		; GRF
	PUSHF	R7
	PUSH	AR5
	PUSH	AR6

	LDIU	17,RC
	LDP	DEF1,DP
	LDIU	*-AR3(2),IR1
	LDIU	IR1,AR6
	.word	40608000h	;	LDFU	0.0,R0
	RPTB	L268

; begin loop 17-1
L268:	STF	R0,*AR6++(1)
; end loop 17-1

	LDIU	*-AR3(5),R0
	CMPI	1,R0
	BLT	L44
	BR	L267

; begin loop 18-1
L43:	LDIU	IR1,AR6
	LDIU	*-AR3(3),AR5
	BRD	L266
	LDIU	2,R1
	LDIU	*-AR3(4),R0
	LDIU	R0,AR0

; begin loop 19-2
L45:	LDIU	2,RC
	RPTB	L269

; begin loop 20-3
	LDIU	AR5,AR2
	LDIU	AR6,AR1
	ADDI	1,AR6
	LDFU	*+AR5(13),R1
	SUBF3	R1,*AR5,R2
	RND	R2, R2		; GRF
	MPYF	*AR0,R2
	LDFU	*+AR5(12),R3
	ADDF3	R3,*+AR5(1),R4
	RND	R4, R4		; GRF
	MPYF	*+AR0(1),R4
	ADDF	R4,R2
	LDFU	*+AR5(6),R4
	LDFU	*++AR2(19),R5
	ADDF3	R5,R4,R6
	RND	R6, R6		; GRF
	MPYF	*+AR0(6),R6
	ADDF	R6,R2
	LDFU	*+AR5(7),R6
	SUBF3	*-AR2(1),R6,R7
	RND	R7, R7		; GRF
	MPYF	*+AR0(7),R7
	ADDF	R7,R2
	ADDF	*+AR5(1),R3
	RND	R3, R3		; GRF
	MPYF	*+AR0(12),R3
	ADDF	R3,R2
	SUBF	*AR5,R1
	RND	R1, R1		; GRF
	MPYF	*+AR0(13),R1
	ADDF	R2,R1
	SUBRF	*+AR5(18),R6
	RND	R6, R6		; GRF
	MPYF	*+AR0(18),R6
	ADDF	R6,R1
	ADDF	R5,R4
	RND	R4, R4		; GRF
	MPYF	*+AR0(19),R4
	ADDF	R4,R1
	ADDF	*AR1,R1
	RND	R1, R1		; GRF
	STF	R1,*AR1
	LDIU	AR6,AR1
	ADDI	1,AR6
	LDFU	*+AR5(12),R1
	ADDF3	R1,*+AR5(1),R2
	RND	R2, R2		; GRF
	MPYF	*AR0++(1),R2
	LDFU	*+AR5(13),R3
	SUBF3	R3,*AR5,R4
	RND	R4, R4		; GRF
	MPYF	*AR0++(5),R4
	SUBF	R4,R2
	LDFU	*+AR5(7),R4
	LDFU	*+AR5(18),R5
	SUBF3	R5,R4,R6
	RND	R6, R6		; GRF
	MPYF	*AR0++(1),R6
	ADDF	R6,R2
	LDFU	*+AR5(6),R6
	ADDF3	*AR2,R6,R7
	RND	R7, R7		; GRF
	MPYF	*AR0++(5),R7
	SUBF	R7,R2
	SUBF	*AR5,R3
	RND	R3, R3		; GRF
	MPYF	*AR0++(1),R3
	ADDF	R3,R2
	ADDF	*+AR5(1),R1
	RND	R1, R1		; GRF
	MPYF	*AR0++(5),R1
	SUBRF	R2,R1
	ADDF	*+AR5(19),R6
	RND	R6, R6		; GRF
	MPYF	*AR0++(1),R6
	ADDF	R6,R1
	SUBRF	R5,R4
	RND	R4, R4		; GRF
	MPYF	*AR0--(17),R4
	SUBF	R4,R1
	ADDF	*AR1,R1
	RND	R1, R1		; GRF
L269:	STF	R1,*AR1
; end loop 20-3

	LDIU	R0,AR0
	ADDI	2,AR5
	LDIU	*+AR3(2),R1
	SUBI	1,R1
L266:	STI	R1,*+AR3(2)
	CMPI	0,R1
	BGE	L45
; end loop 19-2

	LDIU	*-AR3(3),IR0
	ADDI	*-AR3(6),IR0
	ADDI	24,IR0
	STI	IR0,*-AR3(3)
	ADDI	*-AR3(7),R0
	ADDI	24,R0
	STI	R0,*-AR3(4)
	LDIU	*+AR3(1),R0
L267:	SUBI	1,R0
	STI	R0,*+AR3(1)
	CMPI	0,R0
	BGE	L43
; end loop 18-1

L44:	LDIU	*-AR3(5),R0

	LDIU	*+AR3(9),AR6
	LDIU	*+AR3(8),AR5
	LDFU	*+AR3(7),R7
	LDFU	*+AR3(6),R6
	LDIU	*+AR3(5),R5
	LDIU	*+AR3(4),R4
	LDIU	AR3,SP
	POP	AR3
	RETSU	

;
	.def	_sprojTrZm__FPfN21iN24
	.def	_sprojTrZm__FPfN21iN24$LAJ
;=============================================================================
; LOCAL OPTIONS								     |
;   Calling convention     : TI C stack parameters			     |
;									     |
; GENERATED CODE PROPERTIES						     |
;   Total words of code    : 123					     |
;   Volatile registers used: R0,R1,R2,R3,R4,R5,R6,R7,AR0,AR1,AR2,DP,IR0,IR1, |
;			     BK,RS,RE,RC				     |
;   Registers for locals   : AR0	holds tw			     |
;			     AR5	holds tv			     |
;			     AR6	holds tf			     |
;   Parameters             : *-AR3(2)   holds f				     |
;			     *-AR3(3)   holds v				     |
;			     *-AR3(4)   holds w				     |
;			     *-AR3(5)   holds num_blk			     |
;			     *-AR3(6)   holds v_stride			     |
;			     *-AR3(7)   holds w_stride			     |
;   Stack frame            : full (frame pointer in AR3)		     |
;=============================================================================

	.sect	"T:RTN6"

_sprojTrZm__FPfN21iN24$LAJ:
	PUSH	BK
_sprojTrZm__FPfN21iN24:
	PUSH	AR3
	LDIU	SP,AR3
	ADDI	3,SP
	PUSH	R4
	PUSH	R5
	RND	R6, R6		; GRF
	PUSHF	R6
	RND	R7, R7		; GRF
	PUSHF	R7
	PUSH	AR5
	PUSH	AR6

	LDIU	17,RC
	LDP	DEF1,DP
	LDIU	*-AR3(2),IR1
	LDIU	IR1,AR6
	.word	40608000h	;	LDFU	0.0,R0
	RPTB	L306

; begin loop 21-1
L306:	STF	R0,*AR6++(1)
; end loop 21-1

	LDIU	*-AR3(5),R0
	CMPI	1,R0
	BLT	L54
	BR	L305

; begin loop 22-1
L53:	LDIU	IR1,AR6
	LDIU	*-AR3(3),AR5
	BRD	L304
	LDIU	2,R1
	LDIU	*-AR3(4),R0
	LDIU	R0,AR0

; begin loop 23-2
L55:	LDIU	2,RC
	RPTB	L307

; begin loop 24-3
	LDIU	AR5,AR2
	LDIU	AR6,AR1
	ADDI	1,AR6
	LDFU	*+AR5(13),R1
	ADDF3	R1,*AR5,R2
	RND	R2, R2		; GRF
	MPYF	*AR0,R2
	LDFU	*+AR5(12),R3
	SUBF3	R3,*+AR5(1),R4
	RND	R4, R4		; GRF
	MPYF	*+AR0(1),R4
	ADDF	R4,R2
	LDFU	*+AR5(6),R4
	LDFU	*++AR2(19),R5
	SUBF3	R5,R4,R6
	RND	R6, R6		; GRF
	MPYF	*+AR0(6),R6
	ADDF	R6,R2
	LDFU	*+AR5(7),R6
	ADDF3	*-AR2(1),R6,R7
	RND	R7, R7		; GRF
	MPYF	*+AR0(7),R7
	ADDF	R7,R2
	SUBF	*+AR5(1),R3
	RND	R3, R3		; GRF
	MPYF	*+AR0(12),R3
	ADDF	R3,R2
	ADDF	*AR5,R1
	RND	R1, R1		; GRF
	MPYF	*+AR0(13),R1
	ADDF	R2,R1
	ADDF	*+AR5(18),R6
	RND	R6, R6		; GRF
	MPYF	*+AR0(18),R6
	ADDF	R6,R1
	SUBRF	R5,R4
	RND	R4, R4		; GRF
	MPYF	*+AR0(19),R4
	ADDF	R4,R1
	ADDF	*AR1,R1
	RND	R1, R1		; GRF
	STF	R1,*AR1
	LDIU	AR6,AR1
	ADDI	1,AR6
	LDFU	*+AR5(12),R1
	SUBF3	R1,*+AR5(1),R2
	RND	R2, R2		; GRF
	MPYF	*AR0++(1),R2
	LDFU	*+AR5(13),R3
	ADDF3	R3,*AR5,R4
	RND	R4, R4		; GRF
	MPYF	*AR0++(5),R4
	SUBF	R4,R2
	LDFU	*+AR5(7),R4
	LDFU	*+AR5(18),R5
	ADDF3	R5,R4,R6
	RND	R6, R6		; GRF
	MPYF	*AR0++(1),R6
	ADDF	R6,R2
	LDFU	*+AR5(6),R6
	SUBF3	*AR2,R6,R7
	RND	R7, R7		; GRF
	MPYF	*AR0++(5),R7
	SUBF	R7,R2
	ADDF	*AR5,R3
	RND	R3, R3		; GRF
	MPYF	*AR0++(1),R3
	ADDF	R3,R2
	SUBF	*+AR5(1),R1
	RND	R1, R1		; GRF
	MPYF	*AR0++(5),R1
	SUBRF	R2,R1
	SUBRF	*+AR5(19),R6
	RND	R6, R6		; GRF
	MPYF	*AR0++(1),R6
	ADDF	R6,R1
	ADDF	R5,R4
	RND	R4, R4		; GRF
	MPYF	*AR0--(17),R4
	SUBF	R4,R1
	ADDF	*AR1,R1
	RND	R1, R1
L307:	STF	R1,*AR1
; end loop 24-3

	LDIU	R0,AR0
	ADDI	2,AR5
	LDIU	*+AR3(2),R1
	SUBI	1,R1
L304:	STI	R1,*+AR3(2)
	CMPI	0,R1
	BGE	L55
; end loop 23-2

	LDIU	*-AR3(3),IR0
	ADDI	*-AR3(6),IR0
	ADDI	24,IR0
	STI	IR0,*-AR3(3)
	ADDI	*-AR3(7),R0
	ADDI	24,R0
	STI	R0,*-AR3(4)
	LDIU	*+AR3(1),R0
L305:	SUBI	1,R0
	STI	R0,*+AR3(1)
	CMPI	0,R0
	BGE	L53
; end loop 22-1

L54:	LDIU	*-AR3(5),R0

	LDIU	*+AR3(9),AR6
	LDIU	*+AR3(8),AR5
	LDFU	*+AR3(7),R7
	LDFU	*+AR3(6),R6
	LDIU	*+AR3(5),R5
	LDIU	*+AR3(4),R4
	LDIU	AR3,SP
	POP	AR3
	RETSU	

;
	.def	_sprojTrTp__FPfN21iN24
	.def	_sprojTrTp__FPfN21iN24$LAJ
;=============================================================================
; LOCAL OPTIONS								     |
;   Calling convention     : TI C stack parameters			     |
;									     |
; GENERATED CODE PROPERTIES						     |
;   Total words of code    : 123					     |
;   Volatile registers used: R0,R1,R2,R3,R4,R5,R6,R7,AR0,AR1,AR2,DP,IR0,IR1, |
;			     BK,RS,RE,RC				     |
;   Registers for locals   : AR0	holds tw			     |
;			     AR5	holds tv			     |
;			     AR6	holds tf			     |
;   Parameters             : *-AR3(2)   holds f				     |
;			     *-AR3(3)   holds v				     |
;			     *-AR3(4)   holds w				     |
;			     *-AR3(5)   holds num_blk			     |
;			     *-AR3(6)   holds v_stride			     |
;			     *-AR3(7)   holds w_stride			     |
;   Stack frame            : full (frame pointer in AR3)		     |
;=============================================================================

	.sect	"T:RTN7"

_sprojTrTp__FPfN21iN24$LAJ:
	PUSH	BK
_sprojTrTp__FPfN21iN24:
	PUSH	AR3
	LDIU	SP,AR3
	ADDI	3,SP
	PUSH	R4
	PUSH	R5
	RND	R6, R6		; GRF
	PUSHF	R6
	RND	R7, R7		; GRF
	PUSHF	R7
	PUSH	AR5
	PUSH	AR6

	LDIU	17,RC
	LDP	DEF1,DP
	LDIU	*-AR3(2),IR1
	LDIU	IR1,AR6
	.word	40608000h	;	LDFU	0.0,R0
	RPTB	L344

; begin loop 25-1
L344:	STF	R0,*AR6++(1)
; end loop 25-1

	LDIU	*-AR3(5),R0
	CMPI	1,R0
	BLT	L64
	BR	L343

; begin loop 26-1
L63:	LDIU	IR1,AR6
	LDIU	*-AR3(3),AR5
	BRD	L342
	LDIU	2,R1
	LDIU	*-AR3(4),R0
	LDIU	R0,AR0

; begin loop 27-2
L65:	LDIU	2,RC
	RPTB	L345

; begin loop 28-3
	LDIU	AR5,AR2
	LDIU	AR6,AR1
	ADDI	1,AR6
	LDFU	*+AR5(12),R1
	ADDF3	R1,*AR5,R2
	RND	R2, R2		; GRF
	MPYF	*AR0,R2
	LDFU	*+AR5(13),R3
	ADDF3	R3,*+AR5(1),R4
	RND	R4, R4		; GRF
	MPYF	*+AR0(1),R4
	ADDF	R4,R2
	LDFU	*+AR5(6),R4
	LDFU	*++AR2(18),R5
	ADDF3	R5,R4,R6
	RND	R6, R6		; GRF
	MPYF	*+AR0(6),R6
	ADDF	R6,R2
	LDFU	*+AR5(7),R6
	ADDF3	*+AR2(1),R6,R7
	RND	R7, R7		; GRF
	MPYF	*+AR0(7),R7
	ADDF	R7,R2
	ADDF	*AR5,R1
	RND	R1, R1		; GRF
	MPYF	*+AR0(12),R1
	ADDF	R2,R1
	ADDF	*+AR5(1),R3
	RND	R3, R3		; GRF
	MPYF	*+AR0(13),R3
	ADDF	R3,R1
	ADDF	R5,R4
	RND	R4, R4		; GRF
	MPYF	*+AR0(18),R4
	ADDF	R4,R1
	ADDF	*+AR5(19),R6
	RND	R6, R6		; GRF
	MPYF	*+AR0(19),R6
	ADDF	R6,R1
	ADDF	*AR1,R1
	RND	R1, R1		; GRF
	STF	R1,*AR1
	LDIU	AR6,AR1
	ADDI	1,AR6
	LDFU	*+AR5(13),R1
	ADDF3	R1,*+AR5(1),R2
	RND	R2, R2		; GRF
	MPYF	*AR0++(1),R2
	LDFU	*+AR5(12),R3
	ADDF3	R3,*AR5,R4
	RND	R4, R4		; GRF
	MPYF	*AR0++(5),R4
	SUBF	R4,R2
	LDFU	*+AR5(7),R4
	LDFU	*+AR5(19),R5
	ADDF3	R5,R4,R6
	RND	R6, R6		; GRF
	MPYF	*AR0++(1),R6
	ADDF	R6,R2
	LDFU	*+AR5(6),R6
	ADDF3	*AR2,R6,R7
	RND	R7, R7		; GRF
	MPYF	*AR0++(5),R7
	SUBF	R7,R2
	ADDF	*+AR5(1),R1
	RND	R1, R1		; GRF
	MPYF	*AR0++(1),R1
	ADDF	R2,R1
	ADDF	*AR5,R3
	RND	R3, R3		; GRF
	MPYF	*AR0++(5),R3
	SUBF	R3,R1
	ADDF	R5,R4
	RND	R4, R4		; GRF
	MPYF	*AR0++(1),R4
	ADDF	R4,R1
	ADDF	*+AR5(18),R6
	RND	R6, R6		; GRF
	MPYF	*AR0--(17),R6
	SUBF	R6,R1
	ADDF	*AR1,R1
	RND	R1, R1		; GRF
L345:	STF	R1,*AR1
; end loop 28-3

	LDIU	R0,AR0
	ADDI	2,AR5
	LDIU	*+AR3(2),R1
	SUBI	1,R1
L342:	STI	R1,*+AR3(2)
	CMPI	0,R1
	BGE	L65
; end loop 27-2

	LDIU	*-AR3(3),IR0
	ADDI	*-AR3(6),IR0
	ADDI	24,IR0
	STI	IR0,*-AR3(3)
	ADDI	*-AR3(7),R0
	ADDI	24,R0
	STI	R0,*-AR3(4)
	LDIU	*+AR3(1),R0
L343:	SUBI	1,R0
	STI	R0,*+AR3(1)
	CMPI	0,R0
	BGE	L63
; end loop 26-1

L64:	LDIU	*-AR3(5),R0

	LDIU	*+AR3(9),AR6
	LDIU	*+AR3(8),AR5
	LDFU	*+AR3(7),R7
	LDFU	*+AR3(6),R6
	LDIU	*+AR3(5),R5
	LDIU	*+AR3(4),R4
	LDIU	AR3,SP
	POP	AR3
	RETSU	

;
	.def	_sprojTrTm__FPfN21iN24
	.def	_sprojTrTm__FPfN21iN24$LAJ
;=============================================================================
; LOCAL OPTIONS								     |
;   Calling convention     : TI C stack parameters			     |
;									     |
; GENERATED CODE PROPERTIES						     |
;   Total words of code    : 123					     |
;   Volatile registers used: R0,R1,R2,R3,R4,R5,R6,R7,AR0,AR1,AR2,DP,IR0,IR1, |
;			     BK,RS,RE,RC				     |
;   Registers for locals   : AR0	holds tw			     |
;			     AR5	holds tv			     |
;			     AR6	holds tf			     |
;   Parameters             : *-AR3(2)   holds f				     |
;			     *-AR3(3)   holds v				     |
;			     *-AR3(4)   holds w				     |
;			     *-AR3(5)   holds num_blk			     |
;			     *-AR3(6)   holds v_stride			     |
;			     *-AR3(7)   holds w_stride			     |
;   Stack frame            : full (frame pointer in AR3)		     |
;=============================================================================

	.sect	"T:RTN8"

_sprojTrTm__FPfN21iN24$LAJ:
	PUSH	BK
_sprojTrTm__FPfN21iN24:
	PUSH	AR3
	LDIU	SP,AR3
	ADDI	3,SP
	PUSH	R4
	PUSH	R5
	RND	R6, R6		; GRF
	PUSHF	R6
	RND	R7, R7		; GRF
	PUSHF	R7
	PUSH	AR5
	PUSH	AR6

	LDIU	17,RC
	LDP	DEF1,DP
	LDIU	*-AR3(2),IR1
	LDIU	IR1,AR6
	.word	40608000h	;	LDFU	0.0,R0
	RPTB	L382

; begin loop 29-1
L382:	STF	R0,*AR6++(1)
; end loop 29-1

	LDIU	*-AR3(5),R0
	CMPI	1,R0
	BLT	L74
	BR	L381

; begin loop 30-1
L73:	LDIU	IR1,AR6
	LDIU	*-AR3(3),AR5
	BRD	L380
	LDIU	2,R1
	LDIU	*-AR3(4),R0
	LDIU	R0,AR0

; begin loop 31-2
L75:	LDIU	2,RC
	RPTB	L383

; begin loop 32-3
	LDIU	AR5,AR2
	LDIU	AR6,AR1
	ADDI	1,AR6
	LDFU	*+AR5(12),R1
	SUBF3	R1,*AR5,R2
	RND	R2, R2		; GRF
	MPYF	*AR0,R2
	LDFU	*+AR5(13),R3
	SUBF3	R3,*+AR5(1),R4
	RND	R4, R4		; GRF
	MPYF	*+AR0(1),R4
	ADDF	R4,R2
	LDFU	*+AR5(6),R4
	LDFU	*++AR2(18),R5
	SUBF3	R5,R4,R6
	RND	R6, R6		; GRF
	MPYF	*+AR0(6),R6
	ADDF	R6,R2
	LDFU	*+AR5(7),R6
	SUBF3	*+AR2(1),R6,R7
	RND	R7, R7		; GRF
	MPYF	*+AR0(7),R7
	ADDF	R7,R2
	SUBF	*AR5,R1
	RND	R1, R1		; GRF
	MPYF	*+AR0(12),R1
	ADDF	R2,R1
	SUBF	*+AR5(1),R3
	RND	R3, R3		; GRF
	MPYF	*+AR0(13),R3
	ADDF	R3,R1
	SUBRF	R5,R4
	RND	R4, R4		; GRF
	MPYF	*+AR0(18),R4
	ADDF	R4,R1
	SUBRF	*+AR5(19),R6
	RND	R6, R6		; GRF
	MPYF	*+AR0(19),R6
	ADDF	R6,R1
	ADDF	*AR1,R1
	RND	R1, R1		; GRF
	STF	R1,*AR1
	LDIU	AR6,AR1
	ADDI	1,AR6
	LDFU	*+AR5(13),R1
	SUBF3	R1,*+AR5(1),R2
	RND	R2, R2		; GRF
	MPYF	*AR0++(1),R2
	LDFU	*+AR5(12),R3
	SUBF3	R3,*AR5,R4
	RND	R4, R4		; GRF
	MPYF	*AR0++(5),R4
	SUBF	R4,R2
	LDFU	*+AR5(7),R4
	LDFU	*+AR5(19),R5
	SUBF3	R5,R4,R6
	RND	R6, R6		; GRF
	MPYF	*AR0++(1),R6
	ADDF	R6,R2
	LDFU	*+AR5(6),R6
	SUBF3	*AR2,R6,R7
	RND	R7, R7		; GRF
	MPYF	*AR0++(5),R7
	SUBF	R7,R2
	SUBF	*+AR5(1),R1
	RND	R1, R1		; GRF
	MPYF	*AR0++(1),R1
	ADDF	R2,R1
	SUBF	*AR5,R3
	RND	R3, R3		; GRF
	MPYF	*AR0++(5),R3
	SUBF	R3,R1
	SUBRF	R5,R4
	RND	R4, R4		; GRF
	MPYF	*AR0++(1),R4
	ADDF	R4,R1
	SUBRF	*+AR5(18),R6
	RND	R6, R6		; GRF
	MPYF	*AR0--(17),R6
	SUBF	R6,R1
	ADDF	*AR1,R1
	RND	R1, R1		; GRF
L383:	STF	R1,*AR1
; end loop 32-3

	LDIU	R0,AR0
	ADDI	2,AR5
	LDIU	*+AR3(2),R1
	SUBI	1,R1
L380:	STI	R1,*+AR3(2)
	CMPI	0,R1
	BGE	L75
; end loop 31-2

	LDIU	*-AR3(3),IR0
	ADDI	*-AR3(6),IR0
	ADDI	24,IR0
	STI	IR0,*-AR3(3)
	ADDI	*-AR3(7),R0
	ADDI	24,R0
	STI	R0,*-AR3(4)
	LDIU	*+AR3(1),R0
L381:	SUBI	1,R0
	STI	R0,*+AR3(1)
	CMPI	0,R0
	BGE	L73
; end loop 30-1

L74:	LDIU	*-AR3(5),R0

	LDIU	*+AR3(9),AR6
	LDIU	*+AR3(8),AR5
	LDFU	*+AR3(7),R7
	LDFU	*+AR3(6),R6
	LDIU	*+AR3(5),R5
	LDIU	*+AR3(4),R4
	LDIU	AR3,SP
	POP	AR3
	RETSU	

;
	.def	_CWInitsproj_tr$C
	.def	_CWInitsproj_tr$C$LAJ
;=============================================================================
; LOCAL OPTIONS								     |
;   Calling convention     : Tartan register parameters			     |
;									     |
; GENERATED CODE PROPERTIES						     |
;   Total words of code    : 3						     |
;   Volatile registers used: R0,DP,BK					     |
;   Stack frame            : quick (AR3 points to some old frame)	     |
;=============================================================================

	.sect	"T:RTN9"

_CWInitsproj_tr$C:
	POP	BK
_CWInitsproj_tr$C$LAJ:

	LDIU	0,R0

	BU	BK

;









	.sect	".const"
DEF1:	.word	80000000h	;0.0
;=============================================================================
; Total words of code = 987						     |
; Total words of data = 1						     |
;=============================================================================

	.end

