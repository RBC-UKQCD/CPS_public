/* 
 * BAGEL machine generated output.  
 * Copyright UKQCD Collaboration  
 * This software is provided for NON-COMMERCIAL use only,  
 * and may not be redistributed without permission.  
 * It is provided as is, and is not guaranteed fit for any purpose.
 * Written by Peter Boyle
 */  
.text
.align 2
.globl QCDOC_ChRecon_p_su3
QCDOC_ChRecon_p_su3:

	or    %r9 , %r3 , %r3
	la   %r1,	-1408(%r1)
	or    %r11 , %r4 , %r4
	or    %r23 , %r5 , %r5
	addi  %r22 , %r1 , 864
	or    %r24 , %r6 , %r6
	addi  %r13 , %r22 , 64
	lwz   %r24,	0(%r24)
	or.  %r24 , %r24 , %r24
	bf gt,   lab0 
	addi  %r20 , %r11 , 96
	li   %r25,0
	addi  %r15 , %r22 , 160
	li   %r26,32
	addi  %r16 , %r22 , 256
	li   %r27,64
	addi  %r17 , %r22 , 352
	dcbt 	%r25,%r13
	addi  %r18 , %r22 , 448
	lfd   0,	0(%r11)
	addi  %r12 , %r15 , 0
	lfd   6,	16(%r11)
	lfd   1,	8(%r11)
	dcbt 	%r26,%r13
	lfd   7,	24(%r11)
	lfd   2,	48(%r11)
	lfd   20,	64(%r11)
	dcbt 	%r27,%r13
	lfd   3,	56(%r11)
	dcbt 	%r25,%r23
	lfd   21,	72(%r11)
	dcbt 	%r26,%r23
	lfd   8,	0(%r23)
	dcbt 	%r27,%r23
	lfd   9,	8(%r23)
	dcbt 	%r25,%r11
	lfd   10,	16(%r23)
	dcbt 	%r26,%r11
	lfd   11,	24(%r23)
	dcbt 	%r27,%r11
	lfd   12,	32(%r23)
	dcbt 	%r25,%r20
	lfd   13,	40(%r23)
	dcbt 	%r26,%r20
	lfd   14,	48(%r23)
	lfd   15,	56(%r23)
	mtctr 	%r24
	b        lab1 
.align 4
lab1:
	li   %r19,3
	mtctr 	%r19
	b        lab2 
.align 4
lab2:
	addi  %r20 , %r11 , 144
	stfd  26,	48(%r13)
	fmul   24 , 0 , 8
	lfd   4,	96(%r11)
	fmul   25 , 0 , 9
	stfd  27,	56(%r13)
	fmul   26 , 0 , 10
	lfd   22,	112(%r11)
	fmul   27 , 0 , 11
	stfd  28,	64(%r13)
	or    %r10 , %r11 , %r11
	lfd   5,	104(%r11)
	fmul   28 , 6 , 8
	stfd  29,	72(%r13)
	fnmsub 24 , 1 , 9 , 24
	lfd   23,	120(%r11)
	fmul   29 , 6 , 9
	stfd  30,	80(%r13)
	fmadd  25 , 1 , 8 , 25
	lfd   16,	64(%r23)
	fnmsub 26 , 1 , 11 , 26
	stfd  31,	88(%r13)
	fmadd  27 , 1 , 10 , 27
	dcbt 	%r25,%r20
	fnmsub 28 , 7 , 9 , 28
	addi  %r11 , %r11 , 144
	fmadd  29 , 7 , 8 , 29
	lfd   17,	72(%r23)
	fmadd  24 , 2 , 12 , 24
	lfd   18,	80(%r23)
	fmadd  25 , 2 , 13 , 25
	lfd   19,	88(%r23)
	fmadd  26 , 2 , 14 , 26
	addi  %r13 , %r12 , 0
	fmadd  27 , 2 , 15 , 27
	dcbt 	%r26,%r20
	fmadd  28 , 20 , 12 , 28
	addi  %r23 , %r23 , 128
	fmadd  29 , 20 , 13 , 29
	lfd   0,	32(%r10)
	fnmsub 24 , 3 , 13 , 24
	lfd   1,	40(%r10)
	fmadd  25 , 3 , 12 , 25
	la   %r19,	-1(%r19)
	fnmsub 26 , 3 , 15 , 26
	fmadd  27 , 3 , 14 , 27
	dcbt 	%r27,%r20
	fnmsub 28 , 21 , 13 , 28
	fmadd  29 , 21 , 12 , 29
	addi  %r20 , %r20 , 96
	fmadd  24 , 4 , 16 , 24
	fmadd  25 , 4 , 17 , 25
	fmadd  2 , 4 , 18 , 26
	fmadd  3 , 4 , 19 , 27
	dcbt 	%r25,%r23
	fmadd  28 , 22 , 16 , 28
	fmadd  29 , 22 , 17 , 29
	lfd   4,	128(%r10)
	fnmsub 24 , 5 , 17 , 24
	fmadd  25 , 5 , 16 , 25
	fnmsub 2 , 5 , 19 , 2
	fmadd  3 , 5 , 18 , 3
	dcbt 	%r26,%r23
	fnmsub 28 , 23 , 17 , 28
	fmadd  29 , 23 , 16 , 29
	stfd  24,	0(%r12)
	fmul   30 , 6 , 10
	lfd   5,	136(%r10)
	fmul   31 , 6 , 11
	stfd  25,	8(%r12)
	fmul   24 , 0 , 8
	fmul   25 , 0 , 9
	stfd  2,	16(%r12)
	fmul   26 , 0 , 10
	lfd   6,	16(%r11)
	fmul   27 , 0 , 11
	stfd  3,	24(%r12)
	fnmsub 30 , 7 , 11 , 30
	dcbt 	%r27,%r23
	fmadd  31 , 7 , 10 , 31
	fnmsub 24 , 1 , 9 , 24
	lfd   2,	80(%r10)
	fmadd  25 , 1 , 8 , 25
	lfd   3,	88(%r10)
	fnmsub 26 , 1 , 11 , 26
	stfd  28,	32(%r12)
	fmadd  27 , 1 , 10 , 27
	lfd   0,	0(%r11)
	fmadd  9 , 20 , 14 , 30
	stfd  29,	40(%r12)
	fmadd  8 , 20 , 15 , 31
	dcbt 	%r25,%r20
	addi  %r12 , %r12 , 96
	fmadd  24 , 2 , 12 , 24
	fmadd  25 , 2 , 13 , 25
	lfd   1,	8(%r11)
	fmadd  11 , 2 , 14 , 26
	lfd   7,	24(%r11)
	fmadd  10 , 2 , 15 , 27
	lfd   20,	64(%r11)
	fnmsub 9 , 21 , 15 , 9
	dcbt 	%r26,%r20
	fmadd  8 , 21 , 14 , 8
	fnmsub 24 , 3 , 13 , 24
	fmadd  25 , 3 , 12 , 25
	lfd   2,	48(%r11)
	fnmsub 11 , 3 , 15 , 11
	lfd   21,	72(%r11)
	fmadd  10 , 3 , 14 , 10
	lfd   12,	32(%r23)
	fmadd  26 , 22 , 18 , 9
	lfd   13,	40(%r23)
	fmadd  27 , 22 , 19 , 8
	lfd   3,	56(%r11)
	fmadd  28 , 4 , 16 , 24
	lfd   9,	8(%r23)
	fmadd  29 , 4 , 17 , 25
	lfd   8,	0(%r23)
	fmadd  30 , 4 , 18 , 11
	lfd   14,	48(%r23)
	fmadd  31 , 4 , 19 , 10
	lfd   15,	56(%r23)
	fnmsub 26 , 23 , 19 , 26
	lfd   11,	24(%r23)
	fmadd  27 , 23 , 18 , 27
	lfd   10,	16(%r23)
	fnmsub 28 , 5 , 17 , 28
	fmadd  29 , 5 , 16 , 29
	fnmsub 30 , 5 , 19 , 30
	fmadd  31 , 5 , 18 , 31
	bdnz     lab2 
	addi  %r20 , %r11 , 144
	stfd  26,	48(%r13)
	fmul   24 , 0 , 8
	lfd   4,	96(%r11)
	fmul   25 , 0 , 9
	stfd  27,	56(%r13)
	fmul   26 , 0 , 10
	lfd   22,	112(%r11)
	fmul   27 , 0 , 11
	stfd  28,	64(%r13)
	or    %r10 , %r11 , %r11
	lfd   5,	104(%r11)
	fmul   28 , 6 , 8
	stfd  29,	72(%r13)
	fnmsub 24 , 1 , 9 , 24
	lfd   23,	120(%r11)
	fmul   29 , 6 , 9
	stfd  30,	80(%r13)
	fmadd  25 , 1 , 8 , 25
	lfd   16,	64(%r23)
	fnmsub 26 , 1 , 11 , 26
	stfd  31,	88(%r13)
	fmadd  27 , 1 , 10 , 27
	dcbt 	%r25,%r20
	fnmsub 28 , 7 , 9 , 28
	addi  %r11 , %r11 , 144
	fmadd  29 , 7 , 8 , 29
	lfd   17,	72(%r23)
	fmadd  24 , 2 , 12 , 24
	lfd   18,	80(%r23)
	fmadd  25 , 2 , 13 , 25
	lfd   19,	88(%r23)
	fmadd  26 , 2 , 14 , 26
	addi  %r13 , %r12 , 0
	fmadd  27 , 2 , 15 , 27
	dcbt 	%r26,%r20
	fmadd  28 , 20 , 12 , 28
	addi  %r23 , %r23 , 128
	fmadd  29 , 20 , 13 , 29
	lfd   0,	32(%r10)
	fnmsub 24 , 3 , 13 , 24
	lfd   1,	40(%r10)
	fmadd  25 , 3 , 12 , 25
	la   %r24,	-1(%r24)
	fnmsub 26 , 3 , 15 , 26
	fmadd  27 , 3 , 14 , 27
	dcbt 	%r27,%r20
	fnmsub 28 , 21 , 13 , 28
	fmadd  29 , 21 , 12 , 29
	addi  %r20 , %r20 , 96
	fmadd  24 , 4 , 16 , 24
	fmadd  25 , 4 , 17 , 25
	fmadd  2 , 4 , 18 , 26
	fmadd  3 , 4 , 19 , 27
	dcbt 	%r25,%r23
	fmadd  28 , 22 , 16 , 28
	fmadd  29 , 22 , 17 , 29
	lfd   4,	128(%r10)
	fnmsub 24 , 5 , 17 , 24
	fmadd  25 , 5 , 16 , 25
	fnmsub 2 , 5 , 19 , 2
	fmadd  3 , 5 , 18 , 3
	dcbt 	%r26,%r23
	fnmsub 28 , 23 , 17 , 28
	fmadd  29 , 23 , 16 , 29
	stfd  24,	0(%r12)
	fmul   30 , 6 , 10
	lfd   5,	136(%r10)
	fmul   31 , 6 , 11
	stfd  25,	8(%r12)
	fmul   24 , 0 , 8
	fmul   25 , 0 , 9
	stfd  2,	16(%r12)
	fmul   26 , 0 , 10
	fmul   27 , 0 , 11
	stfd  3,	24(%r12)
	fnmsub 30 , 7 , 11 , 30
	dcbt 	%r27,%r23
	fmadd  31 , 7 , 10 , 31
	fnmsub 24 , 1 , 9 , 24
	lfd   2,	80(%r10)
	fmadd  25 , 1 , 8 , 25
	lfd   3,	88(%r10)
	fnmsub 26 , 1 , 11 , 26
	stfd  28,	32(%r12)
	fmadd  27 , 1 , 10 , 27
	lfd   0,	0(%r15)
	fmadd  9 , 20 , 14 , 30
	stfd  29,	40(%r12)
	fmadd  8 , 20 , 15 , 31
	dcbt 	%r25,%r20
	addi  %r12 , %r12 , 96
	fmadd  24 , 2 , 12 , 24
	fmadd  25 , 2 , 13 , 25
	lfd   1,	8(%r15)
	fmadd  11 , 2 , 14 , 26
	lfd   20,	32(%r16)
	fmadd  10 , 2 , 15 , 27
	addi  %r12 , %r15 , 0
	fnmsub 9 , 21 , 15 , 9
	dcbt 	%r26,%r20
	fmadd  8 , 21 , 14 , 8
	fnmsub 24 , 3 , 13 , 24
	fmadd  25 , 3 , 12 , 25
	lfd   2,	16(%r15)
	fnmsub 11 , 3 , 15 , 11
	lfd   13,	8(%r17)
	fmadd  10 , 3 , 14 , 10
	lfd   12,	0(%r17)
	fmadd  26 , 22 , 18 , 9
	lfd   21,	40(%r16)
	fmadd  27 , 22 , 19 , 8
	lfd   3,	24(%r15)
	fmadd  28 , 4 , 16 , 24
	lfd   9,	8(%r16)
	fmadd  29 , 4 , 17 , 25
	lfd   8,	0(%r16)
	fmadd  30 , 4 , 18 , 11
	lfd   22,	48(%r16)
	fmadd  31 , 4 , 19 , 10
	lfd   14,	16(%r17)
	fnmsub 26 , 23 , 19 , 26
	lfd   11,	24(%r16)
	fmadd  27 , 23 , 18 , 27
	lfd   10,	16(%r16)
	fnmsub 28 , 5 , 17 , 28
	fmadd  29 , 5 , 16 , 29
	lfd   15,	24(%r17)
	fnmsub 30 , 5 , 19 , 30
	lfd   17,	8(%r18)
	fmadd  31 , 5 , 18 , 31
	stfd  26,	48(%r13)
	fsub   6 , 1 , 8
	lfd   16,	0(%r18)
	fsub   5 , 11 , 2
	stfd  27,	56(%r13)
	fadd   4 , 3 , 10
	lfd   18,	16(%r18)
	fadd   7 , 0 , 9
	stfd  28,	64(%r13)
	fadd   1 , 1 , 9
	lfd   19,	24(%r18)
	fadd   0 , 0 , 8
	stfd  29,	72(%r13)
	fadd   2 , 2 , 10
	lfd   23,	56(%r16)
	fadd   3 , 3 , 11
	stfd  30,	80(%r13)
	fneg   7 , 7
	lfd   8,	32(%r17)
	fadd   25 , 1 , 13
	stfd  31,	88(%r13)
	fadd   24 , 0 , 12
	lfd   9,	40(%r17)
	fadd   26 , 2 , 14
	lfd   1,	40(%r15)
	fadd   27 , 3 , 15
	lfd   0,	32(%r15)
	fadd   28 , 4 , 13
	lfd   2,	48(%r15)
	fsub   29 , 5 , 12
	lfd   3,	56(%r15)
	fsub   30 , 6 , 15
	lfd   10,	48(%r17)
	fadd   31 , 7 , 14
	lfd   11,	56(%r17)
	fadd   24 , 24 , 16
	lfd   12,	32(%r18)
	fadd   25 , 25 , 17
	lfd   13,	40(%r18)
	fadd   26 , 26 , 18
	lfd   14,	48(%r18)
	fadd   27 , 27 , 19
	lfd   15,	56(%r18)
	fadd   28 , 28 , 16
	addi  %r13 , %r9 , 96
	fadd   29 , 29 , 17
	stfd  24,	0(%r9)
	fadd   30 , 30 , 18
	lfd   16,	64(%r16)
	fadd   31 , 31 , 19
	stfd  25,	8(%r9)
	fadd   4 , 3 , 22
	lfd   17,	72(%r16)
	fsub   5 , 23 , 2
	stfd  26,	16(%r9)
	fsub   6 , 1 , 20
	lfd   18,	80(%r16)
	fadd   7 , 0 , 21
	stfd  27,	24(%r9)
	fadd   1 , 1 , 21
	lfd   19,	88(%r16)
	fadd   0 , 0 , 20
	stfd  28,	32(%r9)
	fadd   2 , 2 , 22
	lfd   21,	72(%r17)
	fadd   3 , 3 , 23
	stfd  29,	40(%r9)
	fneg   7 , 7
	lfd   20,	64(%r17)
	fadd   25 , 1 , 9
	stfd  30,	48(%r9)
	fadd   24 , 0 , 8
	lfd   22,	80(%r17)
	fadd   26 , 2 , 10
	stfd  31,	56(%r9)
	fadd   27 , 3 , 11
	lfd   0,	64(%r15)
	fadd   28 , 4 , 9
	lfd   1,	72(%r15)
	fsub   29 , 5 , 8
	lfd   2,	80(%r15)
	fsub   30 , 6 , 11
	lfd   3,	88(%r15)
	fadd   31 , 7 , 10
	lfd   23,	88(%r17)
	fadd   24 , 24 , 12
	lfd   8,	64(%r18)
	fadd   25 , 25 , 13
	lfd   9,	72(%r18)
	fadd   26 , 26 , 14
	lfd   10,	80(%r18)
	fadd   27 , 27 , 15
	lfd   11,	88(%r18)
	fadd   28 , 28 , 12
	fadd   29 , 29 , 13
	stfd  24,	64(%r9)
	fadd   30 , 30 , 14
	lfd   12,	32(%r23)
	fadd   31 , 31 , 15
	stfd  25,	72(%r9)
	fadd   4 , 3 , 18
	lfd   13,	40(%r23)
	fsub   5 , 19 , 2
	stfd  26,	80(%r9)
	fsub   6 , 1 , 16
	lfd   14,	48(%r23)
	fadd   7 , 0 , 17
	stfd  27,	88(%r9)
	fadd   1 , 1 , 17
	lfd   15,	56(%r23)
	fadd   0 , 0 , 16
	stfd  28,	96(%r9)
	fadd   2 , 2 , 18
	fadd   3 , 3 , 19
	stfd  29,	104(%r9)
	fneg   7 , 7
	fadd   25 , 1 , 21
	stfd  30,	112(%r9)
	fadd   24 , 0 , 20
	fadd   26 , 2 , 22
	stfd  31,	120(%r9)
	fadd   27 , 3 , 23
	lfd   0,	0(%r11)
	fadd   28 , 4 , 21
	lfd   1,	8(%r11)
	fsub   29 , 5 , 20
	lfd   2,	48(%r11)
	fsub   30 , 6 , 23
	lfd   3,	56(%r11)
	fadd   31 , 7 , 22
	lfd   20,	64(%r11)
	fadd   24 , 24 , 8
	lfd   6,	16(%r11)
	fadd   25 , 25 , 9
	lfd   7,	24(%r11)
	fadd   26 , 26 , 10
	lfd   21,	72(%r11)
	fadd   27 , 27 , 11
	fadd   28 , 28 , 8
	fadd   29 , 29 , 9
	stfd  24,	128(%r9)
	fadd   30 , 30 , 10
	lfd   8,	0(%r23)
	fadd   31 , 31 , 11
	stfd  25,	136(%r9)
	lfd   9,	8(%r23)
	addi  %r9 , %r9 , 192
	lfd   10,	16(%r23)
	lfd   11,	24(%r23)
# Instruction order barrier inserted
	or.  %r24 , %r24 , %r24
	bt gt,   lab1 
	stfd  26,	48(%r13)
	stfd  27,	56(%r13)
	stfd  28,	64(%r13)
	stfd  29,	72(%r13)
	stfd  30,	80(%r13)
	stfd  31,	88(%r13)
lab0:
	la   %r1,	1408(%r1)
	blr
