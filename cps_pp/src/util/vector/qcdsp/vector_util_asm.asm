**--------------------------------------------------------------------
**  CVS keywords
**
**  $Author: zs $
**  $Date: 2004-08-18 11:58:09 $
**  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/vector/qcdsp/vector_util_asm.asm,v 1.4 2004-08-18 11:58:09 zs Exp $
**  $Id: vector_util_asm.asm,v 1.4 2004-08-18 11:58:09 zs Exp $
**  $Name: not supported by cvs2svn $
**  $Locker:  $
**  $Revision: 1.4 $
**  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/vector/qcdsp/vector_util_asm.asm,v $
**  $State: Exp $
**
**--------------------------------------------------------------------
*
*  vector_util_asm.asm
*
*

	.version        30

	.global		_uDotXEqual
	.global		_mDotMEqual
	.global		_mDotMPlus
	.global		_vecTimesEquFloat
	.global		_vecAddEquVec
	.global		_vecMinusEquVec
	.global		_vecNegative
	.global		_fTimesV1MinusV2
	.global		_oneMinusfTimesMatrix


	.text


;=====================================================================
; LOCAL OPTIONS							     |
;   Calling convention     : TI C stack parameters		     |
;								     |
; GENERATED CODE PROPERTIES					     |
;   Total words of code    : 96					     |
;   Volatile registers used: R0,R1,R2,R3,AR0,AR1,AR2,DP,BK,RE,RC     |
;   Parameters             : AR0	holds y			     |
;			     AR1	holds u			     |
;			     AR2	holds x			     |
;   Stack frame            : full (frame pointer in AR3)	     |
;=====================================================================

Y	.set	AR0
U	.set	AR1
X	.set	AR2

;=====================================================================
*
* Y = U . X	(Complex 3x3 matrix dot 3 vector)
* 
* R0,R1,R2,R3	are used for parallel instructions
* R0,R1	: temporary storage
* R2  : real part of Y
* R3  : imaginary part of Y
*

	.def	_uDotXEqual

_uDotXEqual:
	POP	BK

	ADDI	1,SP
	PUSH	AR3
	LDIU	SP,AR3
	LDIU	*-AR3(2),Y
	LDIU	*-AR3(3),U
	LDIU	*-AR3(4),X
	LDI     5, IR0		; VECT_LEN - 1


*
* Y = U . X	(Complex 3x3 matrix dot 3 vector)
* 
* R0,R1,R2,R3	are used for parallel instructions
* R0,R1	: temporary storage
* R2  : real part of Y
* R3  : imaginary part of Y
*
* first row of U
	MPYF3	*U,	*X,	R2	;R2 = U[0] * X[0]
	MPYF3	*U++,	*++X,	R3	;R3 = U[0] * X[1]
	MPYF3	*U,	*X,	R0	;R0 = U[1] * X[1]
	MPYF3	*U++,	*-X,	R1
    ||	SUBF3	R0, 	R2,	R2	;R1 = U[1] * X[0] || R2 -= R0
	MPYF3	*U,	*++X,	R0
    ||	ADDF3	R1, 	R3,	R3	;R0 = U[2] * X[2] || R3 += R1
	MPYF3	*U++, *++X,	R1
    ||	ADDF3	R0, 	R2,	R2	;R1 = U[2] * X[3] || R2 += R0
	MPYF3	*U,	*X,	R0
    ||	ADDF3	R1,	R3,	R3	;R0 = U[3] * X[3] || R3 += R1
	MPYF3	*U++,	*-X,	R1
    ||	SUBF3	R0,	R2,	R2	;R1 = U[3] * X[2] || R2 -= R0
	MPYF3	*U,	*++X,	R0
    ||	ADDF3	R1, 	R3,	R3	;R0 = U[4] * X[4] || R3 += R1
	MPYF3	*U++, *++X,	R1
    ||	ADDF3	R0, 	R2,	R2	;R1 = U[4] * X[5] || R2 += R0
	MPYF3	*U,	*X,	R0
    ||	ADDF3	R1,	R3,	R3	;R0 = U[5] * X[5] || R3 += R1
	MPYF3	*U++,	*-X,	R1
    ||	SUBF3	R0,	R2,	R2	;R1 = U[5] * X[4] || R2 -= R0

* second row of U
	MPYF3	*U,	*--X(IR0),R0
    ||	ADDF3	R1,	R3,	R3	;R0 = U[6] * X[0] || R3 += R1
	RND	R2
	LDF	*++X,	R2
    ||	STF	R2,	*Y++		;R2 = X[1] 	   || Y[0] = R2
	RND	R3
	MPYF3	*U++,	R2,	R3
    ||	STF	R3,	*Y++		;R3 = U[6] * X[1] || Y[1] = R3
	MPYF3	*U,	*X,	R2	;R2 = U[7] * X[1] || R3 += R1
	MPYF3	*U++,	*-X,	R1
    ||	SUBF3	R2, 	R0,	R2	;R1 = U[7] * X[0] || R2 = R0-R2
	MPYF3	*U,	*++X,	R0
    ||	ADDF3	R1, 	R3,	R3	;R0 = U[8] * X[2] || R3 += R1
	MPYF3	*U++, *++X,	R1
    ||	ADDF3	R0, 	R2,	R2	;R1 = U[8] * X[3] || R2 += R0
	MPYF3	*U,	*X,	R0
    ||	ADDF3	R1,	R3,	R3	;R0 = U[9] * X[3] || R3 += R1
	MPYF3	*U++,	*-X,	R1
    ||	SUBF3	R0,	R2,	R2	;R1 = U[9] * X[2] || R2 -= R0
	MPYF3	*U,	*++X,	R0
    ||	ADDF3	R1, 	R3,	R3	;R0 = U[10] * X[4] || R3 += R1
	MPYF3	*U++,	*++X,	R1
    ||	ADDF3	R0, 	R2,	R2	;R1 = U[10] * X[5] || R2 += R0
	MPYF3	*U,	*X,	R0
    ||	ADDF3	R1,	R3,	R3	;R0 = U[11] * X[5] || R3 += R1
	MPYF3	*U++,	*-X,	R1
    ||	SUBF3	R0,	R2,	R2	;R1 = U[11] * X[4] || R2 -= R0

* third row of U
	MPYF3	*U,	*--X(IR0), R0
    ||	ADDF3	R1,	R3,	R3	;R0 = U[12] * X[0] || R3 += R1
	RND	R2
	LDF	*++X,	R2		;
    ||	STF	R2,	*Y++		;R2 = X[1] 	    || Y[2] = R2
	RND	R3
	MPYF3	*U++,	R2,	R3
    ||	STF	R3,	*Y++		;R3 = U[12] * X[1] || Y[3] = R3
	MPYF3	*U,	*X,	R2	;R2 = U[13] * X[1] || R3 += R1
	MPYF3	*U++,  *-X,	R1
    ||	SUBF3	R2, 	R0,	R2	;R1 = U[13] * X[0] || R2 = R0-R2
	MPYF3	*U,	*++X,	R0
    ||	ADDF3	R1, 	R3,	R3	;R0 = U[14] * X[2] || R3 += R1
	MPYF3	*U++,	*++X,	R1
    ||	ADDF3	R0, 	R2,	R2	;R1 = U[14] * X[3] || R2 += R0
	MPYF3	*U,	*X,	R0
    ||	ADDF3	R1,	R3,	R3	;R0 = U[15] * X[3] || R3 += R1
	MPYF3	*U++,	*-X,	R1
    ||	SUBF3	R0,	R2,	R2	;R1 = U[15] * X[2] || R2 -= R0
	MPYF3	*U,	*++X,	R0
    ||	ADDF3	R1, 	R3,	R3	;R0 = U[16] * X[4] || R3 += R1
	MPYF3	*U++,	*++X,	R1
    ||	ADDF3	R0, 	R2,	R2	;R1 = U[16] * X[5] || R2 += R0
	MPYF3	*U,	*X,	R0
    ||	ADDF3	R1,	R3,	R3	;R0 = U[17] * X[5] || R3 += R1
	MPYF3	*U,     *-X,    R1
    ||	SUBF3	R0,	R2,	R2	;R1 = U[17] * X[4] || R2 -= R0

	ADDF3	R1,	R3, R3		;R3 += R1
	RND	R2
	RND	R3
    	STF	R2,	*Y		;Y[4] = R2
    ||	STF	R3,	*+Y		;Y[5] = R3

	BUD	BK
	LDIU	AR3,SP
	LDIU	*AR3,AR3
	SUBI	2,SP

;=====================================================================
; LOCAL OPTIONS							     |
;   Calling convention     : TI C stack parameters		     |
;								     |
; GENERATED CODE PROPERTIES					     |
;   Total words of code    : 321				     |
;   Volatile registers used: R0,R1,R2,R3,AR0,AR1,AR2,DP,BK,RE,RC     |
;   Parameters             : AR0	holds C			     |
;			     AR1	holds A			     |
;			     AR2	holds B			     |
;   Stack frame            : full (frame pointer in AR3)	     |
;=====================================================================
C	.set	AR0
A	.set	AR1
B	.set	AR2
;=====================================================================
	.def	_mDotMEqual

_mDotMEqual:
	POP	BK

	ADDI	1,SP
	PUSH	AR3
	LDIU	SP,AR3
	LDIU	AR5, RC
	LDIU	*-AR3(2),C
	LDIU	AR4, RE
	LDIU	*-AR3(3),A
	LDIU	*-AR3(4),B

	LDIU    5, IR0
	LDIU    10, IR1

	LDIU	2,   AR5		; loop counter
	LDIU	2,   AR4		; loop counter

mdme_LOOP0:
	MPYF3	*A,	*B,	R0	;R0 = A[0]*B[0]
	MPYF3	*A++,	*++B,	R3	;R3 = A[0]*B[1]

mdme_LOOP1:
	LDFU	R0,	R2
	MPYF3	*A,	*B,	R0	;R0 = A[1]*B[1]
	MPYF3	*A++,	*-B,	R1	; A=1,B=1 ==> A=2,B=1
    ||	SUBF3	R0, 	R2,	R2	;R1 = A[1]*B[0] || R2 -= R0
	MPYF3	*A, *++B(IR0),	R0	; A=2,B=6
    ||	ADDF3	R1, 	R3,	R3	;R0 = A[2]*B[6] || R3 += R1
	MPYF3	*A++, *++B,	R1	; A=2, B=7 ==> A=3,B=7
    ||	ADDF3	R0, 	R2,	R2	;R1 = A[2]*B[7] || R2 += R0
	MPYF3	*A,	*B,	R0	; A=3,B=7
    ||	ADDF3	R1,	R3,	R3	;R0 = A[3]*B[7] || R3 += R1
	MPYF3	*A++,	*-B,	R1	; A=3,B=7 ==> A=4,B=7
    ||	SUBF3	R0,	R2,	R2	;R1 = A[3]*B[6] || R2 -= R0
	MPYF3	*A, *++B(IR0),	R0	; A=4,B=12
    ||	ADDF3	R1, 	R3,	R3	;R0 = A[4]*B[12] || R3 += R1
	MPYF3	*A++, *++B,	R1	; A=4,B=13 ==> A=5,B=13
    ||	ADDF3	R0, 	R2,	R2	;R1 = A[4]*B[13] || R2 += R0
	MPYF3	*A,	*B,	R0	; A=5,B=13
    ||	ADDF3	R1,	R3,	R3	;R0 = A[5]*B[13] || R3 += R1
	MPYF3	*A,	*-B,	R1	; A=5,B=13
    ||	SUBF3	R0,	R2,	R2	;R1 = U[5]*B[12] || R2 -= R0

	RND	R2
	LDF	*--B(IR1),   R2		; B=3
    ||  STF	R2,     *C++		;R2 = B[3] || C[0] = R2
	DBUD	AR4,    mdme_LOOP1
	MPYF3	*--A(IR0),*-B,  R0	; A=0,B=3
    ||	ADDF3	R1,	R3,	R3	;R0 = A[0]*B[2] || R3 += R1
	RND	R3
	MPYF3	*A++,	R2,	R3	; A=0 ==> A=1,B=3
    ||  STF	R3,	*C++		;R3 = A[0]*B[3] || C[1] = R3
* Branch here

	DBUD	AR5,	mdme_LOOP0
	ADDI	5,	A
	SUBI	7,	B		; B = 0
	LDIU	2,   	AR4		; loop counter
* Branch here


	LDIU    RE,AR4
	LDIU    RC,AR5
	BUD	BK
	LDIU	AR3,SP
	LDIU	*AR3,AR3
	SUBI	2,SP
;=====================================================================
	.def	_mDotMPlus

_mDotMPlus:
	POP	BK

	ADDI	1,SP
	PUSH	AR3
	LDIU	SP,AR3
	LDIU	AR5, RC
	LDIU	*-AR3(2),C
	LDIU	AR4, RE
	LDIU	*-AR3(3),A
	LDIU	*-AR3(4),B

	LDIU    5, 	IR0
	LDIU    11, 	IR1

	LDIU	2,   AR5		; loop counter
	LDIU	2,   AR4		; loop counter

mdmp_LOOP0:
	MPYF3	*A,	*B,	R0	;R0 = A[0]*B[0]

mdmp_LOOP1:
	LDF	*C++,	R2		;R2 = C[0]
    ||	LDF	*+C,	R3		;R3 = C[1]

	MPYF3	*A++,	*++B,	R1	;R1 = A[0]*B[1]
    ||  ADDF3   R0,     R2,     R2      ;R2 += R0
	MPYF3	*A,	*B,	R0	;R0 = A[1]*B[1]
    ||	ADDF3	R1, 	R3,	R3	;R3 += R1
	MPYF3	*A++,	*-B,	R1	; A=1,B=1 ==> A=2,B=1
    ||	SUBF3	R0, 	R2,	R2	;R1 = A[1]*B[0] || R2 -= R0
	MPYF3	*A, *++B(IR0),	R0	; A=2,B=6
    ||	ADDF3	R1, 	R3,	R3	;R0 = A[2]*B[6] || R3 += R1
	MPYF3	*A++, *++B,	R1	; A=2, B=7 ==> A=3,B=7
    ||	ADDF3	R0, 	R2,	R2	;R1 = A[2]*B[7] || R2 += R0
	MPYF3	*A,	*B,	R0	; A=3,B=7
    ||	ADDF3	R1,	R3,	R3	;R0 = A[3]*B[7] || R3 += R1
	MPYF3	*A++,	*-B,	R1	; A=3,B=7 ==> A=4,B=7
    ||	SUBF3	R0,	R2,	R2	;R1 = A[3]*B[6] || R2 -= R0
	MPYF3	*A, *++B(IR0),	R0	; A=4,B=12
    ||	ADDF3	R1, 	R3,	R3	;R0 = A[4]*B[12] || R3 += R1
	MPYF3	*A++, *++B,	R1	; A=4,B=13 ==> A=5,B=13
    ||	ADDF3	R0, 	R2,	R2	;R1 = A[4]*B[13] || R2 += R0
	MPYF3	*A,	*B,	R0	; A=5,B=13
    ||	ADDF3	R1,	R3,	R3	;R0 = A[5]*B[13] || R3 += R1

	MPYF3	*A,	*-B,	R1	; A=5,B=13
    ||	SUBF3	R0,	R2,	R2	;R1 = U[5]*B[12] || R2 -= R0
	MPYF3   *--A(IR0),*--B(IR1),R0  ; A=0,B=2
    ||  ADDF3   R1,     R3,     R3      ;R0 = A[0]*B[2] || R3 += R1
	DBUD	AR4,	mdmp_LOOP1
	RND	R2
	RND	R3
    	STF	R2,	*-C		;C[0] = R2
    ||	STF	R3,	*C++		;C[1] = R3
* Branch here

	DBUD	AR5,	mdmp_LOOP0
	ADDI	6,	A
	SUBI	6,	B		; B = 0
	LDIU	2,   	AR4		; loop counter
* Branch here


	LDIU    RE,AR4
	LDIU    RC,AR5
	BUD	BK
	LDIU	AR3,SP
	LDIU	*AR3,AR3
	SUBI	2,SP


;=====================================================================
; void vecTimesEquFloat(float *a, float c, int n)		     |
;								     |
;   Total words of code    : 15					     |
;   Volatile registers used: R0,R1,R2,AR0,AR2,DP,BK,RS,RE,RC	     |
;   Parameters             : AR0	holds a			     |
;			     R0		holds c			     |
;			     R1		holds n			     |
;   Stack frame            : quick (AR3 points to some old frame)    |
;=====================================================================
	.def	_vecTimesEquFloat

_vecTimesEquFloat:
	POP	BK

	POP	AR0		; a
	POPF	R0		; b
	POP	RC		; n
	ADDI	3,SP

	LDI	AR0,AR1
	MPYF3	*AR1++,R0,R1

	SUBI	2,RC
	RPTB	L65

; begin loop 2-1
	RND	R1,R2
L65:	MPYF3	*AR1++,R0,R1
    ||	STF	R2,*AR0++
; end loop 2-1

	BUD	BK
	RND	R1,R2
    	STF	R2,*AR0
	NOP


;=====================================================================
; void vecAddEquVec(float *a, const float *b, int n)		     |
;								     |
;   Total words of code    : 15					     |
;   Volatile registers used: R0,R1,AR0,AR1,AR2,DP,BK,RS,RE,RC	     |
;   Parameters             : AR0	holds a			     |
;			     AR1	holds b			     |
;			     R0		holds n			     |
;   Stack frame            : quick (AR3 points to some old frame)    |
;=====================================================================
	.def	_vecAddEquVec

_vecAddEquVec:
	POP	BK

	POP	AR0
	POP	AR1
	POP	RC
	ADDI	3,SP

	SUBI	2,RC
	RPTB	L77

; begin loop 3-1
	ADDF3	*AR1++(1),*AR0,R1
	RND	R1,R0
L77:	STF	R0,*AR0++
; end loop 3-1

	BUD	BK
	ADDF3	*AR1,*AR0,R1
	RND	R1,R0
	STF	R0,*AR0


;=====================================================================
; void vecMinusEquVec(float *a, const float *b, int n)		     |
;								     |
;   Total words of code    : 15					     |
;   Volatile registers used: R0,R1,AR0,AR1,AR2,DP,BK,RS,RE,RC	     |
;   Parameters             : AR0	holds a			     |
;			     AR1	holds b			     |
;			     R0		holds n			     |
;   Stack frame            : quick (AR3 points to some old frame)    |
;=====================================================================
	.def	_vecMinusEquVec

_vecMinusEquVec:
	POP	BK

	POP	AR0
	POP	AR1
	POP	RC
	ADDI	3,SP

	SUBI	2,RC
	RPTB	L89

; begin loop 4-1
	SUBF3	*AR1++(1),*AR0,R1
	RND	R1,R0
L89:	STF	R0,*AR0++
; end loop 4-1

	BUD	BK
	SUBF3	*AR1,*AR0,R1
	RND	R1,R0
	STF	R0,*AR0

;=====================================================================
	.def	_vecNegative

_vecNegative:
	POP	BK

	POP	AR0
	POP	AR1
	POP	RC
	ADDI	3,SP

	NEGF	*AR1++(1),R1

	SUBI	2,RC
	RPTB	L91

; begin loop 4-1
	RND	R1,R0
L91:	NEGF	*AR1++(1),R1
    ||	STF	R0,*AR0++
; end loop 4-1

	BUD	BK
	RND	R1,R0
	STF	R0,*AR0
	NOP



;=====================================================================
; void fTimesV1MinusV2(float *a, float b, const float *c,	     |
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
	.def	_fTimesV1MinusV2

_fTimesV1MinusV2:
	POP	BK

	POP	AR0		; a
	POPF	R1		; b
	POP	AR1		; c
	POP	AR2		; d
	POP	RC		; n
	ADDI	5,SP

	MPYF3	*AR1++(1),R1,R3

	SUBI	2,RC
	RPTB	L113

; begin loop 6-1
	SUBF3	*AR2++(1),R3,R0
	RND	R0,R2
L113:	MPYF3	*AR1++(1),R1,R3
    ||	STF	R2,*AR0++(1)
; end loop 6-1

	BUD	BK
	SUBF3	*AR2,R3,R0
	RND	R0,R2
    	STF	R2,*AR0



;=====================================================================
; void oneMinusfTimesMatrix(float *a, float b, const float *c, int n)|
;								     |
;   Total words of code    : 25					     |
;   Volatile registers used: R0,R1,R2,AR0,AR1,AR2,DP,BK,RS,RE,RC     |
;   Registers for locals   : AR2	holds p			     |
;   Parameters             : AR0	holds a			     |
;			     R0		holds b			     |
;			     AR1	holds c			     |
;			     R1		holds n			     |
;   Stack frame            : quick (AR3 points to some old frame)    |
;=====================================================================
	.def	_oneMinusfTimesMatrix


_oneMinusfTimesMatrix:
	POP	BK

	POP	AR0
	POPF	R0
	POP	AR1
	POP	RC
	ADDI	4,SP

	LDIU	AR0,AR2
	NEGF	R0,R0
	RND	R0
	MPYF3	*AR1++(1),R0,R2
	SUBI	2,RC
	RPTB	L121

; begin loop 1-1
	RND	R2,R1
L121:	MPYF3	*AR1++(1),R0,R2
    ||	STF	R1,*AR2++(1)
; end loop 1-1
	RND     R2,R1
	STF     R1,*AR2++(1)

	LDFU	*AR0,R0
	.word	1e00000h	;	ADDF	0.1e1,R0
	RND	R0
	STF	R0,*AR0
	LDFU	*+AR0(8),R1
	.word	1e10000h	;	ADDF	0.1e1,R1
	RND	R1
	STF	R1,*+AR0(8)
	LDFU	*+AR0(16),R2

	BUD	BK
	.word	1e20000h	;	ADDF	0.1e1,R2
	RND	R2
	STF	R2,*+AR0(16)



	
;=====================================================================
;=====================================================================
; The routines below are the direct compiler output. RND's have
; been inserted by hand. However, no hand-made optimizations
; have been performed. This is quick and dirty and should be
; rewritten.
;=====================================================================
;=====================================================================


	.def	_compDotProduct
	.def	_compDotProduct$LAJ
;=============================================================================
; LOCAL OPTIONS								     |
;   Calling convention     : TI C stack parameters			     |
;									     |
; GENERATED CODE PROPERTIES						     |
;   Total words of code    : 43						     |
;   Volatile registers used: R0,R1,AR0,AR1,AR2,BK,RS,RE,RC		     |
;   Registers for locals   : RC		holds i				     |
;   Parameters             : *-AR3(2)   holds c_r			     |
;			     *-AR3(3)   holds c_i			     |
;			     *-AR3(4)   holds a				     |
;			     *-AR3(5)   holds b				     |
;			     *-AR3(6)   holds len			     |
;   Stack frame            : full (frame pointer in AR3)		     |
;=============================================================================


_compDotProduct$LAJ:
	PUSH	BK
_compDotProduct:
	PUSH	AR3
	LDIU	SP,AR3
	LDIU	AR5,RE
	LDIU	AR4,RS

	LDIU	*-AR3(2),AR0
	LDIU	*-AR3(3),AR1
	LDIU	0,RC
	CMPI	*-AR3(6),RC
	BGED	L16
	.word	40608000h	;	LDFU	0.0,R0
	RND	R0
	STF	R0,*AR1
	STF	R0,*AR0

; begin loop 1-1
L9:	LDIU	*-AR3(4),AR2
	LDIU	*-AR3(5),AR4
	MPYF3	*AR4,*AR2,R0
	MPYF3	*+AR4(1),*+AR2(1),R1
	ADDF	R1,R0
	ADDF	*AR0,R0
	RND	R0
	STF	R0,*AR0
	LDIU	*-AR3(5),AR5
	MPYF3	*+AR5(1),*AR2,R0
	LDIU	*-AR3(4),AR2
	MPYF3	*AR4,*+AR2(1),R1
	SUBF	R1,R0
	ADDF	*AR1,R0
	RND	R0
	STF	R0,*AR1
	ADDI	2,RC
	LDIU	*-AR3(4),R1
	ADDI	2,R1
	STI	R1,*-AR3(4)
	CMPI	*-AR3(6),RC
	BLTD	L9
	LDIU	*-AR3(5),R0
	ADDI	2,R0
	STI	R0,*-AR3(5)
; end loop 1-1


L16:	LDIU	RS,AR4
	LDIU	RE,AR5
	LDIU	*-AR3(1),BK
	BUD	BK
	LDIU	AR3,SP
	LDIU	*AR3,AR3
	SUBI	2,SP



	.def	_cTimesV1PlusV2
	.def	_cTimesV1PlusV2$LAJ
;=============================================================================
; LOCAL OPTIONS								     |
;   Calling convention     : TI C stack parameters			     |
;									     |
; GENERATED CODE PROPERTIES						     |
;   Total words of code    : 33						     |
;   Volatile registers used: R0,R1,R2,R3,R4,AR0,AR1,AR2,BK,RE,RC	     |
;   Registers for locals   : RC		holds i				     |
;   Parameters             : AR0	holds a				     |
;			     R0		holds re			     |
;			     R1		holds im			     |
;			     AR1	holds c				     |
;			     AR2	holds d				     |
;			     R2		holds len			     |
;   Stack frame            : full (frame pointer in AR3)		     |
;=============================================================================

_cTimesV1PlusV2:
	POP	BK
_cTimesV1PlusV2$LAJ:
	ADDI	1,SP
	PUSH	AR3
	LDIU	SP,AR3
	LDIU	R4,RE
	LDIU	*-AR3(2),AR0
	LDFU	*-AR3(3),R0
	LDIU	*-AR3(7),R2

	LDIU	0,RC
	CMPI	RC,R2
	BLED	L21
	LDFU	*-AR3(4),R1
	LDIU	*-AR3(5),AR1
	LDIU	*-AR3(6),AR2

; begin loop 2-1
L11:	ADDI	2,RC
	MPYF3	*AR1,R0,R3
	MPYF3	*+AR1(1),R1,R4
	SUBF	R4,R3
	ADDF	*AR2++(1),R3
	RND	R3
	STF	R3,*AR0++(1)
	MPYF3	*+AR1(1),R0,R3
	MPYF3	*AR1,R1,R4
	ADDF	R4,R3
	ADDF	*AR2++(1),R3
	ADDI	2,AR1
	CMPI	RC,R2
	BGTD	L11
	RND	R3
	STF	R3,*AR0++(1)
	NOP
; end loop 2-1


L21:	LDIU	RE,R4
	BUD	BK
	LDIU	AR3,SP
	LDIU	*AR3,AR3
	SUBI	2,SP
