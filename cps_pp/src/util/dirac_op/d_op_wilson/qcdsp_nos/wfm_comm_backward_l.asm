**--------------------------------------------------------------------
**  CVS keywords
**
**  $Author: zs $
**  $Date: 2004-08-18 11:57:59 $
**  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_wilson/qcdsp_nos/wfm_comm_backward_l.asm,v 1.4 2004-08-18 11:57:59 zs Exp $
**  $Id: wfm_comm_backward_l.asm,v 1.4 2004-08-18 11:57:59 zs Exp $
**  $Name: not supported by cvs2svn $
**  $Locker:  $
**  $Revision: 1.4 $
**  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_wilson/qcdsp_nos/wfm_comm_backward_l.asm,v $
**  $State: Exp $
**
**--------------------------------------------------------------------
	.version 30
;=============================================================================
; This file is   : wfm_comm_backward_l.asm				     |
; generated from : wfm_comm_backward_l					     |
; using          : Tartan C/C++ Compiler for the TMS320C3x/4x, v2.1.1	     |
; on (M/D/Y H:M) : 11/02/1998 18:34					     |
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
;   Data Page              : DEFALT					     |
;   Data Page is in ROM    : false					     |
;   *ARn same mem as *ARm? : true					     |
;									     |
;=============================================================================

;
	.def	_wfm_comm_backward_l
	.def	_wfm_comm_backward_l$LAJ
;=============================================================================
; LOCAL OPTIONS								     |
;   Calling convention     : TI C stack parameters			     |
;									     |
; GENERATED CODE PROPERTIES						     |
;   Total words of code    : 57						     |
;   Volatile registers used: R0,R1,R2,R3,AR0,AR1,AR2,IR0,IR1,BK,RS,RE,RC     |
;   Registers for locals   : AR1	holds send_ad			     |
;			     AR2	holds receive_ad		     |
;			     RE		holds i				     |
;			     RC		holds j				     |
;   Parameters             : *-AR3(2)   holds af0			     |
;			     *-AR3(3)   holds af1			     |
;			     *-AR3(4)   holds af2			     |
;			     *-AR3(5)   holds af3			     |
;			     *-AR3(6)   holds wilson_p			     |
;   Stack frame            : full (frame pointer in AR3)		     |
;=============================================================================

	.sect	"T:RTN1"

_wfm_comm_backward_l$LAJ:
	PUSH	BK
_wfm_comm_backward_l:
	PUSH	AR3
	LDIU	SP,AR3
	ADDI	5,SP
	PUSH	AR4
	PUSH	AR5

	LDIU	*-AR3(2),R0
	STI	R0,*+AR3(2)
	LDIU	*-AR3(3),R1
	STI	R1,*+AR3(3)
	LDIU	*-AR3(4),R2
	STI	R2,*+AR3(4)
	LDIU	*-AR3(5),R3
	STI	R3,*+AR3(5)
	LDIU	0,R0
	STI	R0,*+AR3(1)
	LDIU	*-AR3(6),AR5
	LDIU	AR3,AR4
	ADDI	2,AR4
	LDIU	3,IR0
	LDIU	1,BK

; begin loop 1-1
L16:	LDIU	*-AR3(6),AR0
	ADDI	*+AR3(1),AR0
	LDI	*+AR0(15),R0
	BLED	L17
	LDIU	*AR4,AR1
	ADDI3	*+AR0(IR0),AR1,AR2
	LDIU	0,RE

; begin loop 2-2
L14:	LDIU	0,RC
	LDI	*+AR0(11),R0
	BLE	L18

; begin loop 3-3
L12:	LDIU	*AR1++(1),R0
	STI	R0,*AR2++(1)
	ADDI	BK,RC
	CMPI	*+AR0(11),RC
	BLT	L12
; end loop 3-3

L18:	LDIU	*+AR5(7),IR1
	ADDI3	AR1,IR1,RS
	ADDI	BK,RE
	CMPI	*+AR0(15),RE
	BLTD	L14
	SUBI3	BK,RS,AR1
	ADDI	AR2,IR1
	SUBI3	BK,IR1,AR2
; end loop 2-2

L17:	ADDI3	BK,*+AR3(1),R2
	CMPI	4,R2
	BLTD	L16
	ADDI	BK,AR5
	ADDI	BK,AR4
	STI	R2,*+AR3(1)
; end loop 1-1


	LDIU	*+AR3(7),AR5
	LDIU	*+AR3(6),AR4
	LDIU	*-AR3(1),BK
	BUD	BK
	LDIU	AR3,SP
	LDIU	*AR3,AR3
	SUBI	2,SP

;
	.def	_CWInitwfm_comm_backward_l$C
	.def	_CWInitwfm_comm_backward_l$C$LAJ
;=============================================================================
; LOCAL OPTIONS								     |
;   Calling convention     : Tartan register parameters			     |
;									     |
; GENERATED CODE PROPERTIES						     |
;   Total words of code    : 3						     |
;   Volatile registers used: R0,BK					     |
;   Stack frame            : quick (AR3 points to some old frame)	     |
;=============================================================================

	.sect	"T:RTN2"

_CWInitwfm_comm_backward_l$C:
	POP	BK
_CWInitwfm_comm_backward_l$C$LAJ:

	LDIU	0,R0

	BU	BK

;
	.sect	"o:DEFALT"
DEFALT:	.word	0


;=============================================================================
; Total words of code = 60						     |
; Total words of data = 0						     |
;=============================================================================

	.end

