**--------------------------------------------------------------------
**  CVS keywords
**
**  $Author: mcneile $
**  $Date: 2003-06-22 13:34:46 $
**  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/rcomplex/qcdsp/rcomplex_asm.asm,v 1.1.1.1 2003-06-22 13:34:46 mcneile Exp $
**  $Id: rcomplex_asm.asm,v 1.1.1.1 2003-06-22 13:34:46 mcneile Exp $
**  $Name: not supported by cvs2svn $
**  $Locker:  $
**  $Log: not supported by cvs2svn $
**  Revision 1.2  2001/06/19 18:13:36  anj
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
**  Revision 1.2  2001/05/25 06:16:10  cvs
**  Added CVS keywords to phys_v4_0_0_preCVS
**
**  $RCSfile: rcomplex_asm.asm,v $
**  $Revision: 1.1.1.1 $
**  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/rcomplex/qcdsp/rcomplex_asm.asm,v $
**  $State: Exp $
**
**--------------------------------------------------------------------
	.version 30
;=============================================================================
; This file is   : rcomplex.asm						     |
; generated from : rcomplex						     |
; using          : Tartan C/C++ Compiler for the TMS320C3x/4x, v2.1.1	     |
; on (M/D/Y H:M) : 10/30/1997 15:43					     |
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
	.def	_norm__8RcomplexCFv
	.def	_norm__8RcomplexCFv$LAJ
;=============================================================================
; LOCAL OPTIONS								     |
;   Calling convention     : TI C stack parameters			     |
;									     |
; GENERATED CODE PROPERTIES						     |
;   Total words of code    : 9						     |
;   Volatile registers used: R0,R1,R2,AR0,DP,BK				     |
;   Parameters             : AR0	holds TEMP1			     |
;   Stack frame            : quick (AR3 points to some old frame)	     |
;=============================================================================

	.sect	"T:_norm__8RcomplexCFv"

_norm__8RcomplexCFv:
	POP	BK
_norm__8RcomplexCFv$LAJ:
	POP	AR0
	ADDI	1,SP

	LDFU	*AR0,R1
	MPYF	R1,R1

	BUD	BK
	LDFU	*+AR0(1),R2
	MPYF	R2,R2
	ADDF3	R2,R1,R0
	RND	R0;

;
	.def	___apl__8RcomplexFRC8Rcomplex
	.def	___apl__8RcomplexFRC8Rcomplex$LAJ
;=============================================================================
; LOCAL OPTIONS								     |
;   Calling convention     : TI C stack parameters			     |
;									     |
; GENERATED CODE PROPERTIES						     |
;   Total words of code    : 10						     |
;   Volatile registers used: R0,R1,AR0,AR1,DP,BK			     |
;   Parameters             : AR0	holds TEMP3			     |
;			     AR1	holds a				     |
;   Stack frame            : quick (AR3 points to some old frame)	     |
;=============================================================================

	.sect	"T:___apl__8RcomplexFRC8Rcomplex"

___apl__8RcomplexFRC8Rcomplex:
	POP	BK
___apl__8RcomplexFRC8Rcomplex$LAJ:
	POP	AR0
	POP	AR1
	ADDI	2,SP

	ADDF3	*AR1,*AR0,R0
	RND	R0
	STF	R0,*AR0

	ADDF3	*+AR1(1),*+AR0(1),R1
	BUD	BK
	RND	R1
	STF	R1,*+AR0(1)
	LDIU	AR0,R0

;
	.def	___apl__8RcomplexFf
	.def	___apl__8RcomplexFf$LAJ
;=============================================================================
; LOCAL OPTIONS								     |
;   Calling convention     : TI C stack parameters			     |
;									     |
; GENERATED CODE PROPERTIES						     |
;   Total words of code    : 8						     |
;   Volatile registers used: R0,AR0,DP,BK				     |
;   Parameters             : AR0	holds TEMP4			     |
;			     R0		holds a				     |
;   Stack frame            : quick (AR3 points to some old frame)	     |
;=============================================================================

	.sect	"T:___apl__8RcomplexFf"

___apl__8RcomplexFf:
	POP	BK
___apl__8RcomplexFf$LAJ:
	POP	AR0
	POPF	R0
	ADDI	2,SP


	ADDF	*AR0,R0
	BUD	BK
	RND	R0
	STF	R0,*AR0
	LDIU	AR0,R0

;
	.def	___ami__8RcomplexFRC8Rcomplex
	.def	___ami__8RcomplexFRC8Rcomplex$LAJ
;=============================================================================
; LOCAL OPTIONS								     |
;   Calling convention     : TI C stack parameters			     |
;									     |
; GENERATED CODE PROPERTIES						     |
;   Total words of code    : 10						     |
;   Volatile registers used: R0,R1,AR0,AR1,DP,BK			     |
;   Parameters             : AR0	holds TEMP5			     |
;			     AR1	holds a				     |
;   Stack frame            : quick (AR3 points to some old frame)	     |
;=============================================================================

	.sect	"T:___ami__8RcomplexFRC8Rcomplex"

___ami__8RcomplexFRC8Rcomplex:
	POP	BK
___ami__8RcomplexFRC8Rcomplex$LAJ:
	POP	AR0
	POP	AR1
	ADDI	2,SP

	SUBF3	*AR1,*AR0,R0
	RND	R0
	STF	R0,*AR0

	SUBF3	*+AR1(1),*+AR0(1),R1
	BUD	BK
	RND	R1
	STF	R1,*+AR0(1)
	LDIU	AR0,R0

;
	.def	___ami__8RcomplexFf
	.def	___ami__8RcomplexFf$LAJ
;=============================================================================
; LOCAL OPTIONS								     |
;   Calling convention     : TI C stack parameters			     |
;									     |
; GENERATED CODE PROPERTIES						     |
;   Total words of code    : 8						     |
;   Volatile registers used: R0,AR0,DP,BK				     |
;   Parameters             : AR0	holds TEMP6			     |
;			     R0		holds a				     |
;   Stack frame            : quick (AR3 points to some old frame)	     |
;=============================================================================

	.sect	"T:___ami__8RcomplexFf"

___ami__8RcomplexFf:
	POP	BK
___ami__8RcomplexFf$LAJ:
	POP	AR0
	POPF	R0
	ADDI	2,SP


	SUBRF	*AR0,R0
	BUD	BK
	RND	R0
	STF	R0,*AR0
	LDIU	AR0,R0

;
	.def	___amu__8RcomplexFRC8Rcomplex
	.def	___amu__8RcomplexFRC8Rcomplex$LAJ
;=============================================================================
; LOCAL OPTIONS								     |
;   Calling convention     : TI C stack parameters			     |
;									     |
; GENERATED CODE PROPERTIES						     |
;   Total words of code    : 26						     |
;   Volatile registers used: R0,R1,R2,R3,R4,R5,R7,AR0,AR1,DP,BK,RE,RC	     |
;   Registers for locals   : R7		holds sre			     |
;   Parameters             : AR0	holds TEMP7			     |
;			     AR1	holds a				     |
;   Stack frame            : full (frame pointer in AR3)		     |
;=============================================================================

	.sect	"T:___amu__8RcomplexFRC8Rcomplex"

___amu__8RcomplexFRC8Rcomplex:
	POP	BK
___amu__8RcomplexFRC8Rcomplex$LAJ:
	ADDI	1,SP
	PUSH	AR3
	LDIU	SP,AR3
	LDIU	R5,RC
	LDIU	R4,RE
	PUSHF	R7
	LDIU	*-AR3(2),AR0
	LDIU	*-AR3(3),AR1

	LDF	*AR0,R0
     || LDF	*AR1,R1
	MPYF3	R1,R0,R2
	LDF	*+AR0(1),R3
     || LDF	*+AR1(1),R4
	MPYF3	R4,R3,R5
	SUBF3	R5,R2,R7
	MPYF	R4,R0
	MPYF	R3,R1
	ADDF	R1,R0
	RND	R7
	RND	R0
	STF	R7,*AR0
     || STF	R0,*+AR0(1)
	LDIU	AR0,R0

	LDIU	RE,R4
	LDFU	*+AR3(1),R7
	LDIU	RC,R5
	BUD	BK
	LDIU	AR3,SP
	LDIU	*AR3,AR3
	SUBI	2,SP

;
	.def	___amu__8RcomplexFf
	.def	___amu__8RcomplexFf$LAJ
;=============================================================================
; LOCAL OPTIONS								     |
;   Calling convention     : TI C stack parameters			     |
;									     |
; GENERATED CODE PROPERTIES						     |
;   Total words of code    : 9						     |
;   Volatile registers used: R0,R1,AR0,DP,BK				     |
;   Parameters             : AR0	holds TEMP8			     |
;			     R0		holds a				     |
;   Stack frame            : quick (AR3 points to some old frame)	     |
;=============================================================================

	.sect	"T:___amu__8RcomplexFf"

___amu__8RcomplexFf:
	POP	BK
___amu__8RcomplexFf$LAJ:
	POP	AR0
	POPF	R0
	ADDI	2,SP

	MPYF3	R0,*AR0,R1
	RND	R1
	MPYF3	*+AR0(1),R0,R0
     || STF	R1,*AR0

	BUD	BK
	RND	R0
	STF	R0,*+AR0(1)
	LDIU	AR0,R0

;
	.def	___adv__8RcomplexFRC8Rcomplex
        .ref    DEFALT
	.ref	ARTINVERSE32Z
	.def	___adv__8RcomplexFRC8Rcomplex$LAJ
;=============================================================================
; LOCAL OPTIONS								     |
;   Calling convention     : TI C stack parameters			     |
;									     |
; GENERATED CODE PROPERTIES						     |
;   Total words of code    : 34						     |
;   Volatile registers used: R0,R1,R2,R3,R5,R6,R7,AR0,AR2,DP,BK,RC	     |
;   Registers for locals   : R6		holds sre			     |
;			     R7		holds norm2_a_1			     |
;   Parameters             : *-AR3(2)   holds TEMP9			     |
;			     *-AR3(3)   holds a				     |
;   Stack frame            : full (frame pointer in AR3)		     |
;=============================================================================

	.sect	"T:___adv__8RcomplexFRC8Rcomplex"

___adv__8RcomplexFRC8Rcomplex$LAJ:
	PUSH	BK
___adv__8RcomplexFRC8Rcomplex:
	PUSH	AR3
	LDIU	SP,AR3
	LDIU	R5,RC
	PUSHF	R6
	PUSHF	R7

	LDIU	*-AR3(3),AR2
	MPYF3	*AR2,*AR2,R5
	MPYF3	*+AR2(1),*+AR2(1),R0
	ADDF	R5,R0

	LDP	DEFALT,DP
	CALL	ARTINVERSE32Z
	LDFU	R0,R7
	RND	R7

	LDIU	*-AR3(2),AR0
	LDFU	*AR0,R1
	MPYF3	*AR2,R1,R0
	LDFU	*+AR0(1),R2
	MPYF3	*+AR2(1),R2,R3
	ADDF	R3,R0
	RND	R0
	MPYF3	R7,R0,R6
	MPYF	*AR2,R2
	MPYF	*+AR2(1),R1
	SUBRF	R2,R1
	RND	R1
	MPYF	R7,R1
	RND	R6
	RND	R1
	STF	R6,*AR0
     || STF	R1,*+AR0(1)
	LDIU	AR0,R0

	LDIU	RC,R5
	LDFU	*+AR3(2),R7
	LDFU	*+AR3(1),R6
	LDIU	*-AR3(1),BK
	BUD	BK
	LDIU	AR3,SP
	LDIU	*AR3,AR3
	SUBI	2,SP

;
	.def	___adv__8RcomplexFf
	.ref	ARTDIVF32UZ
	.def	___adv__8RcomplexFf$LAJ
;=============================================================================
; LOCAL OPTIONS								     |
;   Calling convention     : TI C stack parameters			     |
;									     |
; GENERATED CODE PROPERTIES						     |
;   Total words of code    : 23						     |
;   Volatile registers used: R0,R1,R2,R3,AR0,AR2,DP,BK			     |
;   Parameters             : AR0	holds TEMP10			     |
;			     R0		holds a				     |
;   Stack frame            : full (frame pointer in AR3)		     |
;=============================================================================

	.sect	"T:___adv__8RcomplexFf"

___adv__8RcomplexFf$LAJ:
	PUSH	BK
___adv__8RcomplexFf:
	PUSH	AR3
	LDIU	SP,AR3
	LDIU	*-AR3(2),AR0
	LDFU	*-AR3(3),R0

	ADDI	1,SP
	LDIU	AR0,AR2
	STF	R0,*+AR3(1)

	LDFU	R0,R1
	LDP	DEFALT,DP
	LDFU	*AR2,R0
	CALL	ARTDIVF32UZ
	RND	R0
	STF	R0,*AR2


	LDP	DEFALT,DP
	LDF	*+AR2(1),R0
     || LDF	*+AR3(1),R1
	CALL	ARTDIVF32UZ
	RND	R0
	STF	R0,*+AR2(1)

	LDIU	AR2,R0

	LDIU	*-AR3(1),BK
	BUD	BK
	LDIU	AR3,SP
	LDIU	*AR3,AR3
	SUBI	2,SP


;=============================================================================
; friend function in Rcomplex class. NEGF operation cause rounding
; error
;=============================================================================
        .def    _conj__FRC8Rcomplex
        .ref    ___ct__8RcomplexFfT1
        .def    _conj__FRC8Rcomplex$LAJ
 
                .sect   "T:_conj__FRC8Rcomplex"
 
_conj__FRC8Rcomplex$LAJ:
        PUSH    BK
_conj__FRC8Rcomplex:
        LDIU    SP,AR0
        LDIU    *-AR0(2),AR1
        LDIU    *-AR0(1),AR0
 
 
        NEGF    *+AR1(1),R0
        RND     R0
        PUSHF   R0
        LDFU    *AR1,R0
        PUSHF   R0
        PUSH    AR0
        CALL    ___ct__8RcomplexFfT1
        SUBI    3,SP
 
 
        RETSU

;=============================================================================
; friend function of Rcomplex: operator negate
;=============================================================================
        .def    ___mi__FRC8Rcomplex
        .ref    ___ct__8RcomplexFfT1
        .def    ___mi__FRC8Rcomplex$LAJ

        .sect   "T:___mi__FRC8Rcomplex"
 
___mi__FRC8Rcomplex$LAJ:
        PUSH    BK
___mi__FRC8Rcomplex:
        LDIU    SP,AR0
        LDIU    *-AR0(2),AR1
        LDIU    *-AR0(1),AR0
 
 
        NEGF    *+AR1(1),R0
        RND     R0
        PUSHF   R0
        NEGF    *AR1,R0
        RND     R0
        PUSHF   R0
        PUSH    AR0
        CALL    ___ct__8RcomplexFfT1
        SUBI    3,SP
 
 
        RETSU  

;=============================================================================
; Total words of data = 0						     |
;=============================================================================

	.end

