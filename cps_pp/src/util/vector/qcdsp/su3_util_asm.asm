**--------------------------------------------------------------------
**  CVS keywords
**
**  $Author: mcneile $
**  $Date: 2003-06-22 13:34:46 $
**  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/vector/qcdsp/su3_util_asm.asm,v 1.1.1.1 2003-06-22 13:34:46 mcneile Exp $
**  $Id: su3_util_asm.asm,v 1.1.1.1 2003-06-22 13:34:46 mcneile Exp $
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
**  $RCSfile: su3_util_asm.asm,v $
**  $Revision: 1.1.1.1 $
**  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/vector/qcdsp/su3_util_asm.asm,v $
**  $State: Exp $
**
**--------------------------------------------------------------------
;
;  su3_util_asm.asm
;
;
;
	.version 30




	.def	_Dagger__6MatrixFPCf
	.def	_Dagger__6MatrixFPCf$LAJ
;=====================================================================
; LOCAL OPTIONS							     |
;   Calling convention     : TI C stack parameters		     |
;								     |
; GENERATED CODE PROPERTIES					     |
;   Total words of code    : 41					     |
;   Volatile registers used: R0,R1,R2,R3,AR0,AR1,DP,BK		     |
;   Parameters             : AR0	holds TEMP22		     |
;			     AR1	holds a			     |
;   Stack frame            : quick (AR3 points to some old frame)    |
;=====================================================================

	.sect	"T:_Dagger__6MatrixFPCf"

_Dagger__6MatrixFPCf:
	POP	BK
_Dagger__6MatrixFPCf$LAJ:
	POP	AR0
	POP	AR1
	ADDI	2,SP

	LDI	5, IR0
	LDI	11, IR1

	LDF	*AR1++,R0		; a[0], AR1->a[1]
	NEGF	*AR1++,R1		; a[1], AR1->a[2]
    ||	STF	R0,*AR0++		; AR0 = a[0]
	RND	R1
	LDF	*AR1++,R0		; a[2]
    ||	STF	R1,*AR0++(IR0)
	NEGF	*AR1++,R1		; a[3]
    ||	STF	R0,*AR0++
	RND	R1
	LDF	*AR1++,R0		; a[4]
    ||	STF	R1,*AR0++(IR0)
	NEGF	*AR1++,R1		; a[5]
    ||	STF	R0,*AR0++
	RND	R1
	LDF	*AR1++,R0		; a[6]
    ||	STF	R1,*AR0--(IR1)
	NEGF	*AR1++,R1		; a[7]
    ||	STF	R0,*AR0++
	RND	R1
	LDF	*AR1++,R0		; a[8]
    ||	STF	R1,*AR0++(IR0)
	NEGF	*AR1++,R1		; a[9]
    ||	STF	R0,*AR0++
	RND	R1
	LDF	*AR1++,R0		; a[10]
    ||	STF	R1,*AR0++(IR0)
	NEGF	*AR1++,R1		; a[11]
    ||	STF	R0,*AR0++
	RND	R1
	LDF	*AR1++,R0		; a[12]
    ||	STF	R1,*AR0--(IR1)
	NEGF	*AR1++,R1		; a[13]
    ||	STF	R0,*AR0++
	RND	R1
	LDF	*AR1++,R0		; a[14]
    ||	STF	R1,*AR0++(IR0)
	NEGF	*AR1++,R1		; a[15]
    ||	STF	R0,*AR0++
	RND	R1
	LDF	*AR1++,R0		; a[16]
    ||	STF	R1,*AR0++(IR0)

	BUD	BK
	NEGF	*AR1,R1			; a[17]
    ||	STF	R0,*AR0++
	RND	R1
	STF	R1,*AR0



;=====================================================================
;								     |
;   Total words of code    : 41					     |
;   Registers for locals   : R7		holds c			     |
;   Parameters             : *-AR3(2)   holds TEMP71		     |
;			     *-AR3(3)   holds dag		     |
;   Stack frame            : full (frame pointer in AR3)	     |
;=====================================================================
	.ref	_vecTimesEquFloat
	.ref	_vecMinusEquVec
	.def	_TrLessAntiHermMatrix__6MatrixFRC6Matrix
	.def	_TrLessAntiHermMatrix__6MatrixFRC6Matrix$LAJ

	.sect	"T:_TrLessAntiHermMatrix__6MatrixFRC6Matrix"

_TrLessAntiHermMatrix__6MatrixFRC6Matrix$LAJ:
	PUSH	BK
_TrLessAntiHermMatrix__6MatrixFRC6Matrix:
	PUSH	AR3
	LDIU	SP,AR3
	ADDI	3,SP

	LDIU	*-AR3(3),R0
	STI	R0,*+AR3(1)
	LDIU	*-AR3(2),AR2
	STI	AR2,*+AR3(3)

	LDIU	18,R1
	PUSH	R1
	PUSH	R0
	PUSH	AR2
	CALL	_vecMinusEquVec


	LDIU	18,R0
	PUSH	R0
	.word	4061f000h	;	LDFU	0.5,R1
	PUSHF	R1
	PUSH	AR2
	CALL	_vecTimesEquFloat
	LDIU	AR2,AR1
	LDIU	AR2,AR0
	ADDI	9,AR0

	PUSH	DP

	LDP	inv3,DP
	LDFU	*+AR2(1),R0
	ADDF3	*AR0,R0,R1
	ADDF	*++AR1(17),R1
	LDFU	@inv3,R2
	RND	R1
	MPYF	R1,R2
	SUBF	R2,R0
	RND	R0
	SUBF3	R2,*AR0,R1
     || STF	R0,*+AR2(1)
	RND	R1
	SUBF3	R2,*AR1,R2
     || STF	R1,*AR0
	RND	R2
	STF	R2,*+AR2(17)

	POP	DP

	LDIU	*-AR3(1),BK
	BUD	BK
	LDIU	AR3,SP
	LDIU	*AR3,AR3
	SUBI	2,SP


;=====================================================================
;								     |
;   Total words of code    : 69					     |
;   Registers for locals   : R3		holds j			     |
;			     R4		holds m			     |
;			     AR5	holds p			     |
;			     AR6	holds p2		     |
;			     RS		holds n			     |
;			     RE		holds m2		     |
;   Parameters             : *-AR3(2)   holds TEMP72		     |
;			     *-AR3(3)   holds v1		     |
;			     *-AR3(4)   holds v2		     |
;   Stack frame            : full (frame pointer in AR3)	     |
;=====================================================================
	.def	_Cross2__6MatrixFRC6VectorT1
	.def	_Cross2__6MatrixFRC6VectorT1$LAJ

	.sect	"T:_Cross2__6MatrixFRC6VectorT1"

_Cross2__6MatrixFRC6VectorT1$LAJ:
	PUSH	BK
_Cross2__6MatrixFRC6VectorT1:
	PUSH	AR3
	LDIU	SP,AR3
	ADDI	5,SP
	PUSH	R4
	PUSH	AR5
	PUSH	AR6

	LDIU	0,R1
	LDIU	*-AR3(3),R0
	STI	R0,*+AR3(4)
	LDIU	R1,R0
	CMPI	2,R0
	BGTD	L347
	LDIU	*-AR3(4),AR6
	LDIU	*-AR3(2),AR5
	STI	R1,*+AR3(2)
	LDIU	1,BK
	LDIU	*+AR3(2),AR0
	ADDI	AR0,AR0

; begin loop 1-1
L352:	ADDI3	BK,AR0,R0
	STI	AR0,*+AR3(3)
	STI	R0,*+AR3(5)
	LDIU	0,R3
	LDIU	*+AR3(2),R1
	LDIU	R1,R2
	ADDI	R1,R1
	ADDI	R1,R2
	BRD	L342
	LDIU	*+AR3(4),AR1
	ADDI	AR1,AR0
	ADDI	*+AR3(5),AR1

; begin loop 2-2
L350:	ADDI3	IR0,R4,RS
	ADDI3	R3,R3,RE
	ADDI	RE,IR0
	ADDI	BK,R3
	STI	IR0,*+AR3(1)
	ADDI3	AR6,RE,AR2
	MPYF3	*AR2,*AR0,R0
	MPYF3	*+AR6(IR0),*AR1,R1
	LDIU	R4,IR1
	ADDF	R1,R0
	ADDF	R0,R0
	RND	R0
	STF	R0,*+AR5(IR1)
	MPYF3	*AR2,*AR1,R0
	MPYF3	*+AR6(IR0),*AR0,R1
	ADDI3	AR5,RS,AR2
	SUBF	R1,R0
	ADDF	R0,R0
	RND	R0
	STF	R0,*AR2
L342:	CMPI	2,R3
	BLED	L350
	ADDI3	R3,R2,R4
	ADDI	R4,R4
	LDIU	BK,IR0
; end loop 2-2

	LDIU	*+AR3(2),R0
	ADDI	BK,R0
	CMPI	2,R0
	BLED	L352
	STI	R0,*+AR3(2)
	LDIU	R0,AR0
	ADDI	AR0,AR0
; end loop 1-1


L347:	LDIU	*+AR3(8),AR6
	LDIU	*+AR3(7),AR5
	LDIU	*+AR3(6),R4
	LDIU	*-AR3(1),BK
	BUD	BK
	LDIU	AR3,SP
	LDIU	*AR3,AR3
	SUBI	2,SP


;=====================================================================
;								     |
;   Total words of code    : 45					     |
;   Volatile registers used: R0,R1,R2,R3,AR0,AR1,AR2,DP,BK,RC	     |
;   Registers for locals   : R3		holds s3		     |
;			     AR2	holds p			     |
;   Parameters             : AR0	holds TEMP73		     |
;			     AR1	holds a			     |
;   Stack frame            : quick (AR3 points to some old frame)    |
;=====================================================================
	.def	_AntiHermMatrix__6MatrixFPCf
	.def	_AntiHermMatrix__6MatrixFPCf$LAJ

	.sect	"T:_AntiHermMatrix__6MatrixFPCf"

_AntiHermMatrix__6MatrixFPCf:
	POP	BK
_AntiHermMatrix__6MatrixFPCf$LAJ:
	POP	AR0
	POP	AR1
	ADDI	2,SP

	PUSH	DP

	LDP	DEF11,DP
	LDFU	*+AR1(7),R3
	LDIU	AR0,AR2
	MPYF	@DEF11,R3
	.word	40608000h	;	LDFU	0.0,R0
	STF	R0,*AR2
	STF	R0,*+AR2(8)
	STF	R0,*+AR2(16)
	LDFU	*+AR1(2),R1
	ADDF	R3,R1
	RND	R1
	STF	R1,*+AR2(1)
	NEGF	*+AR1(2),R2
	ADDF	R3,R2
	RND	R2
	STF	R2,*+AR2(9)
	RND	R3
	.word	0a630800h	;	MPYF	-0.2e1,R3
	RND	R3
	STF	R3,*+AR2(17)
	LDIU	*+AR1(1),RC
	STI	RC,*+AR2(2)
	LDIU	*AR1,R1
	STI	R1,*+AR2(3)
	LDIU	*+AR1(4),R2
	STI	R2,*+AR2(4)
	LDIU	*+AR1(3),RC
	STI	RC,*+AR2(5)
	NEGF	*+AR1(1),R1
	RND	R1
	STF	R1,*+AR2(6)
	LDIU	*AR1,R2
	STI	R2,*+AR2(7)
	LDIU	*+AR1(6),RC
	STI	RC,*+AR2(10)
	LDIU	*+AR1(5),R1
	STI	R1,*+AR2(11)
	NEGF	*+AR1(4),R2
	RND	R2
	STF	R2,*+AR2(12)
	LDIU	*+AR1(3),RC
	STI	RC,*+AR2(13)
	NEGF	*+AR1(6),R1
	RND	R1

	POP	DP

	BUD	BK
	STF	R1,*+AR2(14)
	LDIU	*+AR1(5),R2
	STI	R2,*+AR2(15)


;=====================================================================
;								     |
;   Total words of code    : 113				     |
;   Registers for locals   : R6		holds im0		     |
;			     R7		holds re0		     |
;			     AR6	holds q			     |
;			     AR7	holds p			     |
;   Parameters             : *-AR3(2)   holds TEMP76		     |
;			     *-AR3(3)   holds c			     |
;   Stack frame            : full (frame pointer in AR3)	     |
;=====================================================================
	.def	_Det__6MatrixCFPf
	.def	_Normalize__6VectorFv
	.def	_Orthogonalize__6VectorFRC6Vector
	.ref	DEFALT

	.sect	"T:_Det__6MatrixCFPf"

_Det__6MatrixCFPf$LAJ:
	PUSH	BK
_Det__6MatrixCFPf:
	PUSH	AR3
	LDIU	SP,AR3
	ADDI	7,SP
	LDIU	AR7,RC
	LDIU	AR4,RE
	LDIU	R5,RS
	LDIU	R4,IR1
	PUSHF	R6
	PUSHF	R7

	LDIU	AR3,AR4
	LDIU	AR3,AR0
	LDIU	*-AR3(2),AR7
	LDIU	AR7,AR2
	LDIU	AR7,AR1
	LDFU	*+AR7(8),R0
	LDFU	*++AR1(16),R1
	MPYF3	R1,R0,R2
	LDFU	*+AR7(9),R3
	MPYF3	*+AR1(1),R3,R4
	SUBF	R4,R2
	LDFU	*++AR2(10),R4
	MPYF	*+AR7(14),R4
	SUBF	R4,R2
	MPYF3	*-AR1(1),*+AR2(1),R4
	ADDF3	R4,R2,R7
	MPYF3	*+AR1(1),R0,R2
	MPYF3	R1,R3,R4
	ADDF	R4,R2
	MPYF3	*--AR1(1),*AR2++(1),R4
	SUBF	R4,R2
	MPYF3	*--AR1(1),*AR2--(1),R6
	SUBRF	R2,R6
	LDFU	*+AR7(12),R2
	MPYF3	R2,*AR2,R4
	LDFU	*+AR7(13),R5
	STF	R5,*++AR0(5)
	MPYF	*+AR7(11),R5
	SUBF	R5,R4
	LDFU	*+AR7(6),R5
	STF	R5,*+AR3(6)
	MPYF	R1,R5
	SUBF	R5,R4
	LDFU	*+AR7(7),R5
	STF	R5,*+AR3(7)
	MPYF	*++AR1(3),R5
	ADDF	R5,R4
	RND	R4
	STF	R4,*+AR3(1)
	MPYF3	*AR0++(1),*AR2++(1),R5
	MPYF3	R2,*AR2,R4
	ADDF	R5,R4
	MPYF3	*AR1,*AR0,R5
	MPYF	*+AR3(7),R1
	ADDI	-3,AR1
	SUBF	R5,R4
	SUBRF	R4,R1
	RND	R1
	STF	R1,*+AR3(2)
	MPYF3	*AR1++(1),*AR0++(1),R1
	MPYF3	*AR1,*AR0--(1),R4
	SUBF	R4,R1
	MPYF3	R2,R0,R4
	SUBF	R4,R1
	MPYF3	*-AR0(1),R3,R4
	ADDF	R4,R1
	RND	R1
	STF	R1,*++AR4(3)
	MPYF	R3,R2
	MPYF3	*AR1--(1),*AR0++(1),R1
	MPYF3	*AR1,*AR0,R4
	ADDF	R4,R1
	MPYF	*+AR3(5),R0
	SUBRF	R1,R0
	SUBF	R2,R0
	RND	R0
	STF	R0,*+AR3(4)
	LDIU	*-AR3(3),AR0
	ADDI	-12,AR1
	MPYF3	R7,*AR7,R0
	MPYF3	R6,*+AR7(1),R1
	SUBF	R1,R0
	LDFU	*+AR3(1),R1
	MPYF3	R1,*AR1++(1),R2
	ADDF	R2,R0
	LDFU	*+AR3(2),R2
	MPYF3	R2,*AR1++(1),R3
	SUBF	R3,R0
	LDFU	*+AR3(3),R3
	MPYF3	R3,*AR1++(1),R4
	ADDF	R4,R0
	MPYF3	*++AR4(1),*AR1,R4
	SUBF	R4,R0
	RND	R0
	STF	R0,*AR0
	MPYF	*AR7,R6
	MPYF	*+AR7(1),R7
	ADDF	R7,R6
	MPYF	*+AR7(2),R2
	ADDF	R6,R2
	MPYF	*+AR7(3),R1
	ADDF	R2,R1
	LDFU	*+AR7(4),R2
	MPYF	*+AR3(4),R2
	ADDF	R2,R1
	MPYF	*+AR7(5),R3
	ADDF	R3,R1
	RND	R1
	STF	R1,*+AR0(1)

	LDIU	IR1,R4
	LDIU	RS,R5
	LDFU	*+AR3(9),R7
	LDFU	*+AR3(8),R6
	LDIU	RE,AR4
	LDIU	RC,AR7
	LDIU	*-AR3(1),BK
	BUD	BK
	LDIU	AR3,SP
	LDIU	*AR3,AR3
	SUBI	2,SP


;=====================================================================
;								     |
;   Total words of code    : 8					     |
;   Volatile registers used: R0,R1,AR0,AR2,DP,BK		     |
;   Registers for locals   : AR2	holds p			     |
;   Parameters             : AR0	holds TEMP77		     |
;   Stack frame            : quick (AR3 points to some old frame)    |
;=====================================================================
	.def	_ReTr__6MatrixCFv
	.def	_ReTr__6MatrixCFv$LAJ

	.sect	"T:_ReTr__6MatrixCFv"

_ReTr__6MatrixCFv:
	POP	BK
_ReTr__6MatrixCFv$LAJ:
	POP	AR0
	ADDI	1,SP
	LDFU	*AR0,R0

	BUD	BK
	ADDF	*+AR0(8),R0
	ADDF	*+AR0(16),R0
	RND	R0


;=====================================================================
;								     |
;   Total words of code    : 31					     |
;   Volatile registers used: R0,R1,R2,AR0,AR2,DP,BK		     |
;   Registers for locals   : AR2	holds p			     |
;   Parameters             : AR0	holds TEMP78		     |
;   Stack frame            : quick (AR3 points to some old frame)    |
;=====================================================================
	.def	_NegHalfTrSquare__6MatrixCFv
	.def	_NegHalfTrSquare__6MatrixCFv$LAJ

	.sect	"T:_NegHalfTrSquare__6MatrixCFv"

_NegHalfTrSquare__6MatrixCFv:
	POP	BK
_NegHalfTrSquare__6MatrixCFv$LAJ:
	POP	AR0
	ADDI	1,SP

	LDIU	AR0,AR2
	LDFU	*+AR2(3),R1
	MPYF	R1,R1
	LDFU	*+AR2(2),R2
	MPYF	R2,R2
	ADDF	R2,R1
	LDFU	*+AR2(1),R2
	SUBF	*+AR2(9),R2
	RND	R2
	MPYF	R2,R2
	RND	R2
	.word	0a62e000h	;	MPYF	0.25,R2
	ADDF	R2,R1
	LDFU	*+AR2(5),R2
	MPYF	R2,R2
	ADDF	R2,R1
	LDFU	*+AR2(4),R2
	MPYF	R2,R2
	ADDF	R2,R1
	LDFU	*+AR2(11),R2
	MPYF	R2,R2
	ADDF	R2,R1
	LDFU	*+AR2(10),R2
	MPYF	R2,R2
	ADDF	R2,R1
	LDFU	*+AR2(17),R2
	MPYF	R2,R2
	RND	R2

	BUD	BK
	.word	0a62f400h	;	MPYF	0.75,R2
	ADDF3	R2,R1,R0
	RND	R0


;=====================================================================
;								     |
;   Total words of code    : 34					     |
;   Registers for locals   : R7		holds norm		     |
;   Parameters             : *-AR3(2)   holds TEMP87		     |
;   Stack frame            : full (frame pointer in AR3)	     |
;=====================================================================
	.ref	ARTSQRT32
	.ref	ARTINVERSE32Z
	.ref	_dotProduct

	.sect	"T:_Normalize__6VectorFv"

_Normalize__6VectorFv$LAJ:
	PUSH	BK
_Normalize__6VectorFv:
	PUSH	AR3
	LDIU	SP,AR3
	PUSHF	R7

	PUSH	DP


	LDIU	6,R0
	PUSH	R0
	LDIU	*-AR3(2),R1
	PUSH	R1
	PUSH	R1
	CALL	_dotProduct
	LDFU	R0,R7

	.word	4670000h	;	CMPF	0.1e1,R7
	BEQD	L425
	SUBI	3,SP
	LDFU	*+AR3(1),R7
	LDIU	*-AR3(1),BK

	LDP	DEFALT,DP
	CALL	ARTSQRT32
	RND	R0
	LDFU	R0,R7


	LDP	DEFALT,DP
	CALL	ARTINVERSE32Z
	RND	R0
	LDFU	R0,R7


	LDIU	6,R1
	PUSH	R1
	PUSHF	R7
	LDIU	*-AR3(2),R2
	PUSH	R2
	CALL	_vecTimesEquFloat

	POP	DP

	LDFU	*+AR3(1),R7
	LDIU	*-AR3(1),BK
L425:	BUD	BK
	LDIU	AR3,SP
	LDIU	*AR3,AR3
	SUBI	2,SP


;=====================================================================
; LOCAL OPTIONS							     |
;   Calling convention     : TI C stack parameters		     |
;								     |
; GENERATED CODE PROPERTIES					     |
;   Total words of code    : 80					     |
;   Registers for locals   : R5		holds i			     |
;			     R6		holds re		     |
;			     R7		holds im		     |
;			     AR5	holds n, n		     |
;			     RS		holds q1		     |
;			     RE		holds q			     |
;   Parameters             : AR0	holds TEMP88		     |
;			     AR1	holds p1		     |
;   Stack frame            : full (frame pointer in AR3)	     |
;=====================================================================

	.sect	"T:_Orthogonalize__6VectorFRC6Vector"

_Orthogonalize__6VectorFRC6Vector$LAJ:
	PUSH	BK
_Orthogonalize__6VectorFRC6Vector:
	PUSH	AR3
	LDIU	SP,AR3
	ADDI	2,SP
	LDIU	AR5,IR1
	LDIU	R5,IR0
	PUSH	R4
	PUSHF	R6
	PUSHF	R7
	LDIU	*-AR3(2),AR0
	LDIU	*-AR3(3),AR1

	.word	40668000h	;	LDFU	0.0,R6
	.word	40678000h	;	LDFU	0.0,R7
	LDIU	AR0,RE
	LDIU	AR1,RS
	LDIU	0,R5
	LDIU	1,BK

; begin loop 6-1
L270:	ADDI3	R5,R5,AR0
	ADDI	BK,R5
	STI	AR0,*+AR3(1)
	ADDI3	BK,AR0,AR5
	ADDI3	RS,AR0,AR1
	LDFU	*AR1,R0
	ADDI	RE,AR0
	LDFU	*AR0,R1
	ADDI3	RS,AR5,AR0
	MPYF3	R1,R0,R2
	LDFU	*AR0,R3
	ADDI3	RE,AR5,AR0
	LDFU	*AR0,R4
	STF	R4,*+AR3(2)
	MPYF	R3,R4
	ADDF	R4,R2
	ADDF	R2,R6
	MPYF	*+AR3(2),R0
	CMPI	3,R5
	BLTD	L270
	MPYF	R3,R1
	SUBF	R1,R0
	ADDF	R0,R7
; end loop 6-1

	.word	4668000h	;	CMPF	0.0,R6
	BNED	L431
	LDIU	0,R5
	LDIU	1,RC
	ADDI3	R5,R5,AR0
	CMPF	R7,R6
	BEQD	L432
	LDIU	0,R5
	LDIU	1,RC
	LDFEQ	*+AR3(5),R7

; begin loop 7-1
L275:	ADDI3	R5,R5,AR0
L431:	ADDI	RC,R5
	STI	AR0,*+AR3(1)
	ADDI3	RC,AR0,AR5
	ADDI3	RE,AR0,AR1
	ADDI	RS,AR0
	RND	R6
	MPYF3	*AR0,R6,R0
	ADDI3	RS,AR5,AR2
	RND	R7
	MPYF3	*AR2,R7,R1
	SUBF	R1,R0
	SUBRF	*AR1,R0
	RND	R0
	STF	R0,*AR1
	ADDI3	RE,AR5,AR1
	MPYF3	*AR2,R6,R0
	MPYF3	*AR0,R7,R1
	ADDF	R1,R0
	CMPI	3,R5
	BLTD	L275
	SUBRF	*AR1,R0
	RND	R0
	STF	R0,*AR1
; end loop 7-1


	LDFU	*+AR3(5),R7
L432:	LDFU	*+AR3(4),R6
	LDIU	*+AR3(3),R4
	LDIU	IR1,AR5
	LDIU	IR0,R5
	LDIU	*-AR3(1),BK
	BUD	BK
	LDIU	AR3,SP
	LDIU	*AR3,AR3
	SUBI	2,SP


	.data
OWN:	.word	0fe2aaaabh	;0.33333333333333333333
inv3	.set	OWN

	.sect	".const"
DEF11:	.word	0ff13cd36h	;0.57735


	.end

