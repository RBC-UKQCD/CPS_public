*
*  wmat.asm
*
*

	.version        30

	.global		_wmatMultwmat
	.global		_Tracewmatwmat

	.text

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
	.def	_wmatMultwmat

_wmatMultwmat:
	POP	BK

	ADDI	1,SP
	PUSH	AR3
	LDIU	SP,AR3
	LDIU	AR5, RC
	LDIU	*-AR3(2),C
	LDIU	AR4, RE
	LDIU	*-AR3(3),A
	LDIU	*-AR3(4),B

	LDIU    23, IR0
	LDIU    262, IR1

	LDIU	11,   AR5		; loop counter
	LDIU	11,   AR4		; loop counter

mdme_LOOP0:
	MPYF3	*A,	*B,	R0	;R0 = A[0]*B[0]
	MPYF3	*A++,	*++B,	R3	;R3 = A[0]*B[1]

mdme_LOOP1:
	LDFU	R0,	R2
	MPYF3	*A,	*B,	R0	;R0 = A[1]*B[1]
	MPYF3	*A++,	*-B,	R1	; A=1,B=1 ==> A=2,B=1
    ||	SUBF3	R0, 	R2,	R2	;R1 = A[1]*B[0] || R2 -= R0
	MPYF3	*A, *++B(IR0),	R0	; A=2,B=24
    ||	ADDF3	R1, 	R3,	R3	;R0 = A[2]*B[6] || R3 += R1
	MPYF3	*A++, *++B,	R1	; A=2,B=25 ==> A=3,B=7
    ||	ADDF3	R0, 	R2,	R2	;R1 = A[2]*B[7] || R2 += R0
	MPYF3	*A,	*B,	R0	; A=3,B=25
    ||	ADDF3	R1,	R3,	R3	;R0 = A[3]*B[7] || R3 += R1
	MPYF3	*A++,	*-B,	R1	; A=3,B=25 ==> A=4,B=7
    ||	SUBF3	R0,	R2,	R2	;R1 = A[3]*B[6] || R2 -= R0
	MPYF3	*A, *++B(IR0),	R0	; A=4,B=48
    ||	ADDF3	R1, 	R3,	R3	;R0 = A[4]*B[12] || R3 += R1
	MPYF3	*A++, *++B,	R1	; A=4,B=49 ==> A=5,B=13
    ||	ADDF3	R0, 	R2,	R2	;R1 = A[4]*B[13] || R2 += R0
	MPYF3	*A,	*B,	R0	; A=5,B=49
    ||	ADDF3	R1,	R3,	R3	;R0 = A[5]*B[13] || R3 += R1
	MPYF3	*A++,	*-B,	R1	; A=5,B=49
    ||	SUBF3	R0,	R2,	R2	;R1 = A[5]*B[12] || R2 -= R0
	MPYF3	*A, *++B(IR0),	R0	; A=6;B=72
    ||	ADDF3	R1, 	R3,	R3	;R0 = A[4]*B[12] || R3 += R1
	MPYF3	*A++, *++B,	R1	; A=6,B=73 ==> A=5,B=13
    ||	ADDF3	R0, 	R2,	R2	;R1 = A[4]*B[13] || R2 += R0
	MPYF3	*A,	*B,	R0	; A=7,B=73
    ||	ADDF3	R1,	R3,	R3	;R0 = A[5]*B[13] || R3 += R1
	MPYF3	*A++,	*-B,	R1	; A=7,B=73
    ||	SUBF3	R0,	R2,	R2	;R1 = A[5]*B[12] || R2 -= R0
	MPYF3	*A, *++B(IR0),	R0	; A=8;B=96
    ||	ADDF3	R1, 	R3,	R3	;R0 = A[4]*B[12] || R3 += R1
	MPYF3	*A++, *++B,	R1	; A=8,B=97 ==> A=5,B=13
    ||	ADDF3	R0, 	R2,	R2	;R1 = A[4]*B[13] || R2 += R0
	MPYF3	*A,	*B,	R0	; A=9,B=97
    ||	ADDF3	R1,	R3,	R3	;R0 = A[5]*B[13] || R3 += R1
	MPYF3	*A++,	*-B,	R1	; A=9,B=97
    ||	SUBF3	R0,	R2,	R2	;R1 = A[5]*B[12] || R2 -= R0
	MPYF3	*A, *++B(IR0),	R0	; A=10;B=120
    ||	ADDF3	R1, 	R3,	R3	;R0 = A[4]*B[12] || R3 += R1
	MPYF3	*A++, *++B,	R1	; A=10,B=121 ==> A=5,B=13
    ||	ADDF3	R0, 	R2,	R2	;R1 = A[4]*B[13] || R2 += R0
	MPYF3	*A,	*B,	R0	; A=11,B=121
    ||	ADDF3	R1,	R3,	R3	;R0 = A[5]*B[13] || R3 += R1
	MPYF3	*A++,	*-B,	R1	; A=11,B=121
    ||	SUBF3	R0,	R2,	R2	;R1 = A[5]*B[12] || R2 -= R0
	MPYF3	*A, *++B(IR0),	R0	; A=12;B=144
    ||	ADDF3	R1, 	R3,	R3	;R0 = A[4]*B[12] || R3 += R1
	MPYF3	*A++, *++B,	R1	; A=12,B=145 ==> A=5,B=13
    ||	ADDF3	R0, 	R2,	R2	;R1 = A[4]*B[13] || R2 += R0
	MPYF3	*A,	*B,	R0	; A=13,B=145
    ||	ADDF3	R1,	R3,	R3	;R0 = A[5]*B[13] || R3 += R1
	MPYF3	*A++,	*-B,	R1	; A=13,B=145
    ||	SUBF3	R0,	R2,	R2	;R1 = A[5]*B[12] || R2 -= R0
	MPYF3	*A, *++B(IR0),	R0	; A=14;B=168
    ||	ADDF3	R1, 	R3,	R3	;R0 = A[4]*B[12] || R3 += R1
	MPYF3	*A++, *++B,	R1	; A=14,B=169 ==> A=5,B=13
    ||	ADDF3	R0, 	R2,	R2	;R1 = A[4]*B[13] || R2 += R0
	MPYF3	*A,	*B,	R0	; A=15,B=169
    ||	ADDF3	R1,	R3,	R3	;R0 = A[5]*B[13] || R3 += R1
	MPYF3	*A++,	*-B,	R1	; A=15,B=169
    ||	SUBF3	R0,	R2,	R2	;R1 = A[5]*B[12] || R2 -= R0
	MPYF3	*A, *++B(IR0),	R0	; A=16;B=192
    ||	ADDF3	R1, 	R3,	R3	;R0 = A[4]*B[12] || R3 += R1
	MPYF3	*A++, *++B,	R1	; A=16,B=193 ==> A=5,B=13
    ||	ADDF3	R0, 	R2,	R2	;R1 = A[4]*B[13] || R2 += R0
	MPYF3	*A,	*B,	R0	; A=17,B=193
    ||	ADDF3	R1,	R3,	R3	;R0 = A[5]*B[13] || R3 += R1
	MPYF3	*A++,	*-B,	R1	; A=17,B=193
    ||	SUBF3	R0,	R2,	R2	;R1 = A[5]*B[12] || R2 -= R0
	MPYF3	*A, *++B(IR0),	R0	; A=18;B=216
    ||	ADDF3	R1, 	R3,	R3	;R0 = A[4]*B[12] || R3 += R1
	MPYF3	*A++, *++B,	R1	; A=18,B=217 ==> A=5,B=13
    ||	ADDF3	R0, 	R2,	R2	;R1 = A[4]*B[13] || R2 += R0
	MPYF3	*A,	*B,	R0	; A=19,B=217
    ||	ADDF3	R1,	R3,	R3	;R0 = A[5]*B[13] || R3 += R1
	MPYF3	*A++,	*-B,	R1	; A=19,B=217
    ||	SUBF3	R0,	R2,	R2	;R1 = A[5]*B[12] || R2 -= R0
	MPYF3	*A, *++B(IR0),	R0	; A=20;B=240
    ||	ADDF3	R1, 	R3,	R3	;R0 = A[4]*B[12] || R3 += R1
	MPYF3	*A++, *++B,	R1	; A=20,B=241 ==> A=5,B=13
    ||	ADDF3	R0, 	R2,	R2	;R1 = A[4]*B[13] || R2 += R0
	MPYF3	*A,	*B,	R0	; A=21,B=241
    ||	ADDF3	R1,	R3,	R3	;R0 = A[5]*B[13] || R3 += R1
	MPYF3	*A++,	*-B,	R1	; A=21,B=241
    ||	SUBF3	R0,	R2,	R2	;R1 = A[5]*B[12] || R2 -= R0
	MPYF3	*A, *++B(IR0),	R0	; A=22;B=264
    ||	ADDF3	R1, 	R3,	R3	;R0 = A[4]*B[12] || R3 += R1
	MPYF3	*A++, *++B,	R1	; A=22,B=265 ==> A=5,B=13
    ||	ADDF3	R0, 	R2,	R2	;R1 = A[4]*B[13] || R2 += R0
	MPYF3	*A,	*B,	R0	; A=23,B=265
    ||	ADDF3	R1,	R3,	R3	;R0 = A[5]*B[13] || R3 += R1
	MPYF3	*A,	*-B,	R1	; A=23,B=265
    ||	SUBF3	R0,	R2,	R2	;R1 = A[5]*B[12] || R2 -= R0

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
	ADDI	23,	A
	SUBI	25,	B		; B = 0
	LDIU	11,   	AR4		; loop counter
* Branch here


	LDIU    RE,AR4
	LDIU    RC,AR5
	BUD	BK
	LDIU	AR3,SP
	LDIU	*AR3,AR3
	SUBI	2,SP
;=====================================================================

	.def	_Tracewmatwmat

_Tracewmatwmat:
	POP	BK

	ADDI	1,SP
	PUSH	AR3
	LDIU	SP,AR3
	LDIU	*-AR3(2),C
	LDIU	AR4, RE
	LDIU	*-AR3(3),A
	LDIU	*-AR3(4),B

	LDIU    23, IR0
	LDIU    263, IR1

	LDIU	11,   AR4		; loop counter

	PUSH	R4
	PUSH	R5

	LDF	*C++,	R4		; initialize real sum
	LDF	*C--,	R5		; initialize imag sum

; multiply Nth row of A by Nth column of B.
; accumulate real part of trace in R4
; accumulate imag part of trace in R5

tr_LOOP0:
	MPYF3	*A,	*B,	R0	;R0 = A[0]*B[0]
	MPYF3	*A++,	*++B,	R3	;R3 = A[0]*B[1]
	LDFU	R0,	R2		;R2 = A[0]*B[0]
	MPYF3	*A,	*B,	R0	;R0 = A[1]*B[1]
	MPYF3	*A++,	*-B,	R1	; A=1,B=1 ==> A=2,B=1
    ||	SUBF3	R0, 	R2,	R2	;R1 = A[1]*B[0] || R2 -= R0
	MPYF3	*A, *++B(IR0),	R0	; A=2,B=24
    ||	ADDF3	R1, 	R3,	R3	;R0 = A[2]*B[6] || R3 += R1
	MPYF3	*A++, *++B,	R1	; A=2,B=25 ==> A=3,B=7
    ||	ADDF3	R0, 	R2,	R2	;R1 = A[2]*B[7] || R2 += R0
	MPYF3	*A,	*B,	R0	; A=3,B=25
    ||	ADDF3	R1,	R3,	R3	;R0 = A[3]*B[7] || R3 += R1
	MPYF3	*A++,	*-B,	R1	; A=3,B=25 ==> A=4,B=7
    ||	SUBF3	R0,	R2,	R2	;R1 = A[3]*B[6] || R2 -= R0
	MPYF3	*A, *++B(IR0),	R0	; A=4,B=48
    ||	ADDF3	R1, 	R3,	R3	;R0 = A[4]*B[12] || R3 += R1
	MPYF3	*A++, *++B,	R1	; A=4,B=49 ==> A=5,B=13
    ||	ADDF3	R0, 	R2,	R2	;R1 = A[4]*B[13] || R2 += R0
	MPYF3	*A,	*B,	R0	; A=5,B=49
    ||	ADDF3	R1,	R3,	R3	;R0 = A[5]*B[13] || R3 += R1
	MPYF3	*A++,	*-B,	R1	; A=5,B=49
    ||	SUBF3	R0,	R2,	R2	;R1 = A[5]*B[12] || R2 -= R0
	MPYF3	*A, *++B(IR0),	R0	; A=6;B=72
    ||	ADDF3	R1, 	R3,	R3	;R0 = A[6]*B[12] || R3 += R1
	MPYF3	*A++, *++B,	R1	; A=6,B=73 ==> A=5,B=13
    ||	ADDF3	R0, 	R2,	R2	;R1 = A[6]*B[13] || R2 += R0
	MPYF3	*A,	*B,	R0	; A=7,B=73
    ||	ADDF3	R1,	R3,	R3	;R0 = A[7]*B[13] || R3 += R1
	MPYF3	*A++,	*-B,	R1	; A=7,B=73
    ||	SUBF3	R0,	R2,	R2	;R1 = A[7]*B[12] || R2 -= R0
	MPYF3	*A, *++B(IR0),	R0	; A=8;B=96
    ||	ADDF3	R1, 	R3,	R3	;R0 = A[8]*B[12] || R3 += R1
	MPYF3	*A++, *++B,	R1	; A=8,B=97 ==> A=5,B=13
    ||	ADDF3	R0, 	R2,	R2	;R1 = A[8]*B[13] || R2 += R0
	MPYF3	*A,	*B,	R0	; A=9,B=97
    ||	ADDF3	R1,	R3,	R3	;R0 = A[9]*B[13] || R3 += R1
	MPYF3	*A++,	*-B,	R1	; A=9,B=97
    ||	SUBF3	R0,	R2,	R2	;R1 = A[9]*B[12] || R2 -= R0
	MPYF3	*A, *++B(IR0),	R0	; A=10;B=120
    ||	ADDF3	R1, 	R3,	R3	;R0 = A[10]*B[12] || R3 += R1
	MPYF3	*A++, *++B,	R1	; A=10,B=121 ==> A=5,B=13
    ||	ADDF3	R0, 	R2,	R2	;R1 = A[10]*B[13] || R2 += R0
	MPYF3	*A,	*B,	R0	; A=11,B=121
    ||	ADDF3	R1,	R3,	R3	;R0 = A[11]*B[13] || R3 += R1
	MPYF3	*A++,	*-B,	R1	; A=11,B=121
    ||	SUBF3	R0,	R2,	R2	;R1 = A[11]*B[12] || R2 -= R0
	MPYF3	*A, *++B(IR0),	R0	; A=12;B=144
    ||	ADDF3	R1, 	R3,	R3	;R0 = A[12]*B[12] || R3 += R1
	MPYF3	*A++, *++B,	R1	; A=12,B=145 ==> A=5,B=13
    ||	ADDF3	R0, 	R2,	R2	;R1 = A[12]*B[13] || R2 += R0
	MPYF3	*A,	*B,	R0	; A=13,B=145
    ||	ADDF3	R1,	R3,	R3	;R0 = A[13]*B[13] || R3 += R1
	MPYF3	*A++,	*-B,	R1	; A=13,B=145
    ||	SUBF3	R0,	R2,	R2	;R1 = A[13]*B[12] || R2 -= R0
	MPYF3	*A, *++B(IR0),	R0	; A=14;B=168
    ||	ADDF3	R1, 	R3,	R3	;R0 = A[14]*B[12] || R3 += R1
	MPYF3	*A++, *++B,	R1	; A=14,B=169 ==> A=5,B=13
    ||	ADDF3	R0, 	R2,	R2	;R1 = A[14]*B[13] || R2 += R0
	MPYF3	*A,	*B,	R0	; A=15,B=169
    ||	ADDF3	R1,	R3,	R3	;R0 = A[15]*B[13] || R3 += R1
	MPYF3	*A++,	*-B,	R1	; A=15,B=169
    ||	SUBF3	R0,	R2,	R2	;R1 = A[15]*B[12] || R2 -= R0
	MPYF3	*A, *++B(IR0),	R0	; A=16;B=192
    ||	ADDF3	R1, 	R3,	R3	;R0 = A[16]*B[12] || R3 += R1
	MPYF3	*A++, *++B,	R1	; A=16,B=193 ==> A=5,B=13
    ||	ADDF3	R0, 	R2,	R2	;R1 = A[16]*B[13] || R2 += R0
	MPYF3	*A,	*B,	R0	; A=17,B=193
    ||	ADDF3	R1,	R3,	R3	;R0 = A[17]*B[13] || R3 += R1
	MPYF3	*A++,	*-B,	R1	; A=17,B=193
    ||	SUBF3	R0,	R2,	R2	;R1 = A[17]*B[12] || R2 -= R0
	MPYF3	*A, *++B(IR0),	R0	; A=18;B=216
    ||	ADDF3	R1, 	R3,	R3	;R0 = A[18]*B[12] || R3 += R1
	MPYF3	*A++, *++B,	R1	; A=18,B=217 ==> A=5,B=13
    ||	ADDF3	R0, 	R2,	R2	;R1 = A[18]*B[13] || R2 += R0
	MPYF3	*A,	*B,	R0	; A=19,B=217
    ||	ADDF3	R1,	R3,	R3	;R0 = A[19]*B[13] || R3 += R1
	MPYF3	*A++,	*-B,	R1	; A=19,B=217
    ||	SUBF3	R0,	R2,	R2	;R1 = A[19]*B[12] || R2 -= R0
	MPYF3	*A, *++B(IR0),	R0	; A=20;B=240
    ||	ADDF3	R1, 	R3,	R3	;R0 = A[20]*B[12] || R3 += R1
	MPYF3	*A++, *++B,	R1	; A=20,B=241 ==> A=5,B=13
    ||	ADDF3	R0, 	R2,	R2	;R1 = A[20]*B[13] || R2 += R0
	MPYF3	*A,	*B,	R0	; A=21,B=241
    ||	ADDF3	R1,	R3,	R3	;R0 = A[21]*B[13] || R3 += R1
	MPYF3	*A++,	*-B,	R1	; A=21,B=241
    ||	SUBF3	R0,	R2,	R2	;R1 = A[21]*B[12] || R2 -= R0
	MPYF3	*A, *++B(IR0),	R0	; A=22;B=264
    ||	ADDF3	R1, 	R3,	R3	;R0 = A[22]*B[12] || R3 += R1
	MPYF3	*A++, *++B,	R1	; A=22,B=265 ==> A=5,B=13
    ||	ADDF3	R0, 	R2,	R2	;R1 = A[22]*B[13] || R2 += R0
	MPYF3	*A,	*B,	R0	; A=23,B=265
    ||	ADDF3	R1,	R3,	R3	;R0 = A[23]*B[13] || R3 += R1
	MPYF3	*A++,	*-B,	R1	; A=23,B=265
    ||	SUBF3	R0,	R2,	R2	;R1 = A[23]*B[12] || R2 -= R0
    	ADDF3	R1,	R3,	R3	;		     R3 += R1


	DBUD	AR4,	tr_LOOP0
	ADDF	R2,	R4
	ADDF	R3,	R5
	LDFU	*--B(IR1),	R2	; B is now top of next col.
					; A is start of next row.
* Branch here

	RND	R4
	RND	R5

	STF	R4,	*C++
	STF	R5,	*C

	POP	R5
	POP	R4

	LDIU    RE,AR4
	BUD	BK
	LDIU	AR3,SP
	LDIU	*AR3,AR3
	SUBI	2,SP
;=====================================================================


