**--------------------------------------------------------------------
**  CVS keywords
**
**  $Author: mcneile $
**  $Date: 2003-06-22 13:34:46 $
**  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_stag/qcdsp/dirac_serial.asm,v 1.1.1.1 2003-06-22 13:34:46 mcneile Exp $
**  $Id: dirac_serial.asm,v 1.1.1.1 2003-06-22 13:34:46 mcneile Exp $
**  $Name: not supported by cvs2svn $
**  $Locker:  $
**  $Log: not supported by cvs2svn $
**  Revision 1.2  2001/06/19 18:12:46  anj
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
**  Revision 1.2  2001/05/25 06:16:06  cvs
**  Added CVS keywords to phys_v4_0_0_preCVS
**
**  $RCSfile: dirac_serial.asm,v $
**  $Revision: 1.1.1.1 $
**  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_stag/qcdsp/dirac_serial.asm,v $
**  $State: Exp $
**
**--------------------------------------------------------------------
	.version 30
;=============================================================================
; This file is   : dirac_serial.asm					     |
; generated from : dirac_serial						     |
; using          : Tartan C/C++ Compiler for the TMS320C3x/4x, v2.1.1	     |
; on (M/D/Y H:M) : 03/08/1999 14:01					     |
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
	.ref	_SCUTrans__FPP9SCUDirArgi
	.ref	_SCUTrans__FP9SCUDirArgPUii
	.ref	_SCUTransComplete__Fv
	.def	_dirac_SCU__FPP9SCUDirArgT1PfiPiPPiT6N25
	.ref	_Addr__9SCUDirArgFPv
	.def	_dirac_SCU__FPP9SCUDirArgT1PfiPiPPiT6N25$LAJ
;=============================================================================
; LOCAL OPTIONS								     |
;   Calling convention     : TI C stack parameters			     |
;									     |
; GENERATED CODE PROPERTIES						     |
;   Total words of code    : 121					     |
;   Volatile registers used: R0,R1,R2,R3,R6,R7,AR0,AR1,AR2,DP,IR0,IR1,BK,ST, |
;			     IE,IF,IOF,RS,RE,RC				     |
;   Parameters             : *-AR3(2)   holds Xarg			     |
;			     *-AR3(3)   holds Rarg			     |
;			     *-AR3(4)   holds a				     |
;			     *-AR3(5)   holds a_odd			     |
;			     *-AR3(6)   holds Xoffset			     |
;			     *-AR3(7)   holds ToffsetP			     |
;			     *-AR3(8)   holds ToffsetM			     |
;			     *-AR3(9)   holds countP			     |
;			     *-AR3(10)  holds countM			     |
;   Stack frame            : full (frame pointer in AR3)		     |
;=============================================================================

	.sect	"T:RTN1"

_dirac_SCU__FPP9SCUDirArgT1PfiPiPPiT6N25$LAJ:
	PUSH	BK
_dirac_SCU__FPP9SCUDirArgT1PfiPiPPiT6N25:
	PUSH	AR3
	LDIU	SP,AR3
	ADDI	2,SP
	PUSH	R7
	PUSH	AR4
	PUSH	AR5
	PUSH	AR6
	PUSH	AR7

	LDIU	*-AR3(6),AR5
	LDIU	AR5,AR7
	ADDI	1,AR7
	LDIU	*-AR3(2),AR4
	ADDI	1,AR4
	STI	AR4,*+AR3(1)
	LDIU	AR4,AR6

	LDIU	*-AR3(4),AR0
	ADDI	*AR7++(1),AR0
	PUSH	AR0
	LDIU	*AR6++(1),R0
	PUSH	R0
	CALL	_Addr__9SCUDirArgFPv


	LDIU	*-AR3(4),AR0
	ADDI	*AR7++(1),AR0
	PUSH	AR0
	LDIU	*AR6++(1),R0
	PUSH	R0
	CALL	_Addr__9SCUDirArgFPv


	LDIU	*-AR3(4),AR0
	ADDI	*AR7++(1),AR0
	PUSH	AR0
	LDIU	*AR6++(1),R0
	PUSH	R0
	CALL	_Addr__9SCUDirArgFPv


	LDIU	*-AR3(4),R0
	PUSH	R0
	LDIU	*-AR3(2),AR0
	LDIU	*AR0,R1
	PUSH	R1
	CALL	_Addr__9SCUDirArgFPv


	LDIU	3,R0
	PUSH	R0
	LDIU	*+AR3(1),R1
	PUSH	R1
	CALL	_SCUTrans__FPP9SCUDirArgi


	LDIU	4,R0
	PUSH	R0
	ADDI	*-AR3(3),R0
	PUSH	R0
	CALL	_SCUTrans__FPP9SCUDirArgi


	LDIU	*-AR3(9),AR0
	ADDI	*-AR3(5),AR0
	LDIU	*AR0,R0
	PUSH	R0
	LDIU	*-AR3(7),AR0
	ADDI	*-AR3(5),AR0
	LDIU	*AR0,R0
	PUSH	R0
	LDIU	*-AR3(2),AR0
	LDIU	*AR0,R0
	PUSH	R0
	CALL	_SCUTrans__FP9SCUDirArgPUii
	SUBI	15,SP


	CALL	_SCUTransComplete__Fv

	LDIU	AR5,AR6
	ADDI	5,AR6
	LDIU	*-AR3(2),AR5
	ADDI	5,AR5
	LDIU	AR5,AR7
	LDIU	2,AR4
	STI	AR4,*+AR3(1)

; begin loop 2-1
L21:

	LDIU	*-AR3(4),AR0
	ADDI	*AR6++(1),AR0
	PUSH	AR0
	LDIU	*AR7++(1),R0
	PUSH	R0
	CALL	_Addr__9SCUDirArgFPv

	LDI	*+AR3(1),R0
	BNED	L21
	SUBI	1,R0
	STI	R0,*+AR3(1)
	SUBI	2,SP
; end loop 2-1


	LDIU	*-AR3(4),R0
	PUSH	R0
	LDIU	*-AR3(2),AR0
	LDIU	*+AR0(4),R1
	PUSH	R1
	CALL	_Addr__9SCUDirArgFPv


	LDIU	3,R0
	PUSH	R0
	PUSH	AR5
	CALL	_SCUTrans__FPP9SCUDirArgi


	LDIU	4,R0
	PUSH	R0
	LDIU	*-AR3(3),R1
	PUSH	R1
	CALL	_SCUTrans__FPP9SCUDirArgi


	LDIU	*-AR3(10),AR0
	ADDI	*-AR3(5),AR0
	LDIU	*AR0,R0
	PUSH	R0
	LDIU	*-AR3(8),AR0
	ADDI	*-AR3(5),AR0
	LDIU	*AR0,R0
	PUSH	R0
	LDIU	*-AR3(2),AR0
	LDIU	*+AR0(4),R0
	PUSH	R0
	CALL	_SCUTrans__FP9SCUDirArgPUii
	SUBI	9,SP


	CALL	_SCUTransComplete__Fv


	LDIU	*+AR3(7),AR7
	LDIU	*+AR3(6),AR6
	LDIU	*+AR3(5),AR5
	LDIU	*+AR3(4),AR4
	LDIU	*+AR3(3),R7
	LDIU	*-AR3(1),BK
	BUD	BK
	LDIU	AR3,SP
	LDIU	*AR3,AR3
	SUBI	2,SP

;
	.ref	_cmv
	.def	_dirac_cmv__FiPPfPfN31T3
	.def	_dirac_cmv__FiPPfPfN31T3$LAJ
;=============================================================================
; LOCAL OPTIONS								     |
;   Calling convention     : TI C stack parameters			     |
;									     |
; GENERATED CODE PROPERTIES						     |
;   Total words of code    : 65						     |
;   Volatile registers used: R0,R1,R2,R3,R4,R5,R6,R7,AR0,AR1,AR2,DP,IR0,IR1, |
;			     BK,ST,IE,IF,IOF,RS,RE,RC			     |
;   Registers for locals   : AR0	holds fp0			     |
;			     R5		holds i				     |
;			     AR6	holds dest			     |
;			     AR7	holds fp1			     |
;   Parameters             : *-AR3(2)   holds sites			     |
;			     *-AR3(3)   holds chi			     |
;			     *-AR3(4)   holds u				     |
;			     *-AR3(5)   holds a				     |
;			     *-AR3(6)   holds b				     |
;			     *-AR3(7)   holds add_flag			     |
;			     *-AR3(8)   holds cmv_vector		     |
;   Stack frame            : full (frame pointer in AR3)		     |
;=============================================================================

	.sect	"T:RTN2"

	.def	_cmv_vector
_cmv_vector
	.space	30h

;====================================================================
;  Copy function arguments to _cmv_stack so code runs faster.
;====================================================================
	.def	_cmv_stack
	.space	10h
_cmv_stack
	.space	10h

Y	.set		AR0
U	.set		AR2
X	.set		AR1

_dirac_cmv__FiPPfPfN31T3$LAJ:
	PUSH	BK
_dirac_cmv__FiPPfPfN31T3:
	PUSH	AR3
	LDIU	SP,AR3
	LDIU	SP,IR1
	PUSH	R4
	PUSH	R5
	PUSH	AR4
	PUSH	AR5
	PUSH	AR6
	PUSH	AR7

;  Copy argument stack to CRAM

	LDI	80h, AR2
	LSH	16, AR2
	OR	_cmv_stack, AR2

	LDI	10, R0
	SUBI3	R0, AR3, AR1
	SUBI3	R0, AR2, AR0

	LDI	*AR1++, R0
	RPTS	20
	LDI	*AR1++, R0
     ||	STI	R0, *AR0++
	STI	R0, *AR0++

	LDI	AR2, AR3	;  set AR3 to cmv_stack

;  End copy argument stack

;  Set IR0 for cmv.  Also add flag

	LDI	2,	IR0
	FLOAT	*-AR3(7),R4	;  add flag

	LDI	*-AR3(2),R0
	BLED	L69
	LDIU	*-AR3(4),U	;  U pointer
	LDIU	0,R5
	LDIU	*+AR3(6),AR7
	LDIU	*-AR3(3),AR5
	ADDI	R5,AR5

; begin loop 5-1
L71:	LDIU	7,RC

	LDIU	*+AR5(8),R0
	ADDI	*-AR3(6),R0
	LDI	R0, Y

	LDIU	256,R1
	LSH	16,R1

	LDIU	*-AR3(8),AR6

	RPTB	L72

; begin loop 6-2
	LDIU	*AR5++(1),AR7
	TSTB	AR7,R1
	BEQ	L60

	ADDI	*-AR3(5),AR7
L60:	LDI	*AR7++(1),R0
	LDI	*AR7++(1),R0
     ||	STI	R0,*AR6++(1)
	LDI	*AR7++(1),R0
     ||	STI	R0,*AR6++(1)
	LDI	*AR7++(1),R0
     ||	STI	R0,*AR6++(1)
	LDI	*AR7++(1),R0
     ||	STI	R0,*AR6++(1)
	LDI	*AR7++(1),R0
     ||	STI	R0,*AR6++(1)
L72:	STI	R0,*AR6++(1)
; end loop 6-2

************************************************************************
*       cmv.asm
*	2/14/99  RDM
*
*	derived from dirac.asm, written by Dong Chen.
*
*	cmv( (FLOAT *) y, (FLOAT *) U, (FLOAT * )x, float add_flag )
*
************************************************************************

;  First row
	LDIU	*-AR3(8), X

	LDF	0,	R1
	LDI	23,	RC

	MPYF3   *Y,	R4,	R2 ;  R4 = 1.0 for add, 0.0 otherwise
	MPYF3   *+Y,	R4,	R3 ;  R4 = 1.0 for add, 0.0 otherwise

	RPTB	loop0

	MPYF3	*U,	*X++,	R0	;R0 = U[0] * X1[0]
    ||	ADDF3	R1, 	R3,	R3

	MPYF3	*U++,	*X,	R1	;R1 = U[0] * X1[1]
    ||	ADDF3	R0, 	R2,	R2

	MPYF3	*U,	*X--,	R0	;R0 = U[1] * X1[1]
    ||	ADDF3	R1, 	R3,	R3

loop0
	MPYF3	*U++,	*X++(IR0),R1	;R1 = U[1] * X1[0]
    ||	SUBF3	R0, 	R2,	R2

	RND	R2,	R2
	STF	R2,	*Y++
      	ADDF3	R1, 	R3,	R3
	RND	R3,	R3
	STF	R3,	*Y++

;  Second row

	LDIU	*-AR3(8), X

	LDF	0,	R1
	LDI	23,	RC

	MPYF3   *Y,	R4,	R2 ;  R4 = 1.0 for add, 0.0 otherwise
	MPYF3   *+Y,	R4,	R3 ;  R4 = 1.0 for add, 0.0 otherwise

	RPTB	loop1

	MPYF3	*U,	*X++,	R0	;R0 = U[0] * X1[0]
    ||	ADDF3	R1, 	R3,	R3

	MPYF3	*U++,	*X,	R1	;R1 = U[0] * X1[1]
    ||	ADDF3	R0, 	R2,	R2

	MPYF3	*U,	*X--,	R0	;R0 = U[1] * X1[1]
    ||	ADDF3	R1, 	R3,	R3

loop1
	MPYF3	*U++,	*X++(IR0),R1	;R1 = U[1] * X1[0]
    ||	SUBF3	R0, 	R2,	R2

	RND	R2,	R2
	STF	R2,	*Y++
      	ADDF3	R1, 	R3,	R3
	RND	R3,	R3
	STF	R3,	*Y++

;  Thrid row

	LDIU	*-AR3(8), X

	LDF	0,	R1
	LDI	23,	RC

	MPYF3   *Y,	R4,	R2 ;  R4 = 1.0 for add, 0.0 otherwise
	MPYF3   *+Y,	R4,	R3 ;  R4 = 1.0 for add, 0.0 otherwise

	RPTB	loop2

	MPYF3	*U,	*X++,	R0	;R0 = U[0] * X1[0]
    ||	ADDF3	R1, 	R3,	R3

	MPYF3	*U++,	*X,	R1	;R1 = U[0] * X1[1]
    ||	ADDF3	R0, 	R2,	R2

	MPYF3	*U,	*X--,	R0	;R0 = U[1] * X1[1]
    ||	ADDF3	R1, 	R3,	R3

loop2
	MPYF3	*U++,	*X++(IR0),R1	;R1 = U[1] * X1[0]
    ||	SUBF3	R0, 	R2,	R2

	RND	R2,	R2
	STF	R2,	*Y++
      	ADDF3	R1, 	R3,	R3
	RND	R3,	R3
	STF	R3,	*Y++

;  End cmv insertion
******************************************************

	ADDI	9,R5
	CMPI	*-AR3(2),R5
	BLTD	L71
;	SUBI	4,SP
	LDIU	*-AR3(3),AR5
	NOP
	ADDI	R5,AR5
; end loop 5-1


	LDIU	*+AR3(6),AR7
L69:	LDIU	*+AR3(5),AR6
	LDIU	*+AR3(4),AR5
	LDIU	*+AR3(3),AR4
	LDIU	*+AR3(2),R5
	LDIU	*+AR3(1),R4
	LDIU	*-AR3(1),BK
	LDIU	IR1,AR3
	BUD	BK
	LDIU	AR3,SP
	LDIU	*AR3,AR3
	SUBI	2,SP

;
	.def	_CWInitdirac_serial$C
	.def	_CWInitdirac_serial$C$LAJ
;=============================================================================
; LOCAL OPTIONS								     |
;   Calling convention     : Tartan register parameters			     |
;									     |
; GENERATED CODE PROPERTIES						     |
;   Total words of code    : 3						     |
;   Volatile registers used: R0,DP,BK					     |
;   Stack frame            : quick (AR3 points to some old frame)	     |
;=============================================================================

	.sect	"T:RTN3"

_CWInitdirac_serial$C:
	POP	BK
_CWInitdirac_serial$C$LAJ:

	LDIU	0,R0

	BU	BK

;



;=============================================================================
; Total words of code = 216						     |
; Total words of data = 0						     |
;=============================================================================

	.end


