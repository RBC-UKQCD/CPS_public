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
.globl QCDOC_Clover_mat_mlt
QCDOC_Clover_mat_mlt:

	or    %r10 , %r3 , %r3
	la   %r1,	-736(%r1)
	or    %r11 , %r4 , %r4
	or    %r12 , %r5 , %r5
	or    %r13 , %r6 , %r6
	lwz   %r9,	0(%r10)
	or.  %r9 , %r9 , %r9
	bf gt,   lab0 
	addi  %r14 , %r13 , 0
	li   %r15,96
	li   %r16,128
	li   %r17,160
	li   %r18,288
	li   %r19,320
	li   %r20,352
	li   %r21,384
	li   %r22,416
	li   %r23,448
	li   %r24,480
	li   %r25,512
	li   %r26,544
	lfd   12,	0(%r12)
	lfd   13,	8(%r12)
	lfd   14,	16(%r12)
	lfd   15,	24(%r12)
	lfd   16,	32(%r12)
	lfd   17,	40(%r12)
	lfd   24,	0(%r11)
	lfd   25,	8(%r11)
	lfd   26,	16(%r11)
	lfd   27,	24(%r11)
	mtctr 	%r9
	b        lab1 
.align 4
lab1:
	fmul   0 , 24 , 12
	lfd   28,	32(%r11)
	fmul   1 , 24 , 13
	lfd   29,	40(%r11)
	fmul   2 , 25 , 12
	stfd  7,	56(%r14)
	fmul   3 , 25 , 13
	lfd   24,	48(%r11)
	lfd   18,	48(%r12)
	fmadd  0 , 25 , 14 , 0
	lfd   19,	56(%r12)
	fmadd  1 , 25 , 15 , 1
	stfd  8,	64(%r14)
	fnmsub 2 , 26 , 13 , 2
	stfd  9,	72(%r14)
	fmadd  3 , 26 , 12 , 3
	lfd   25,	56(%r11)
	fmul   4 , 28 , 12
	stfd  10,	80(%r14)
	fmadd  0 , 26 , 15 , 0
	stfd  11,	88(%r14)
	fnmsub 1 , 26 , 14 , 1
	lfd   20,	64(%r12)
	fmadd  2 , 27 , 14 , 2
	lfd   21,	72(%r12)
	fmadd  3 , 27 , 15 , 3
	lfd   26,	64(%r11)
	fmul   5 , 28 , 13
	lfd   22,	80(%r12)
	fnmsub 4 , 29 , 13 , 4
	lfd   27,	72(%r11)
	fmadd  0 , 28 , 16 , 0
	lfd   23,	88(%r12)
	fmadd  1 , 28 , 17 , 1
	addi  %r14 , %r13 , 0
	fmadd  2 , 24 , 16 , 2
	dcbt 	%r15,%r12
	fmadd  5 , 29 , 12 , 5
	lfd   28,	80(%r11)
	fmadd  4 , 24 , 14 , 4
	la   %r13,	96(%r13)
	fmadd  0 , 29 , 17 , 0
	la   %r9,	-1(%r9)
	fnmsub 1 , 29 , 16 , 1
	dcbt 	%r16,%r12
	fmadd  3 , 24 , 17 , 3
	fmadd  5 , 24 , 15 , 5
	lfd   29,	88(%r11)
	fnmsub 4 , 25 , 15 , 4
	fmadd  2 , 25 , 17 , 2
	lfd   24,	96(%r11)
	fmul   6 , 27 , 12
	dcbt 	%r17,%r12
	fnmsub 3 , 25 , 16 , 3
	fmadd  5 , 25 , 14 , 5
	la   %r12,	96(%r12)
	fmadd  4 , 26 , 16 , 4
	fmul   7 , 27 , 13
	lfd   25,	104(%r11)
	fnmsub 6 , 28 , 13 , 6
	dcbt 	%r18,%r11
	fmadd  0 , 27 , 18 , 0
	fmadd  5 , 26 , 17 , 5
	fmadd  1 , 27 , 19 , 1
	fmadd  7 , 28 , 12 , 7
	lfd   26,	112(%r11)
	fmadd  6 , 29 , 14 , 6
	lfd   27,	120(%r11)
	fmadd  0 , 28 , 19 , 0
	dcbt 	%r19,%r11
	fmadd  2 , 29 , 18 , 2
	fnmsub 1 , 28 , 18 , 1
	fmadd  7 , 29 , 15 , 7
	fnmsub 6 , 24 , 15 , 6
	lfd   28,	192(%r11)
	fmadd  3 , 29 , 19 , 3
	dcbt 	%r20,%r11
	fmadd  2 , 24 , 19 , 2
	fmadd  4 , 25 , 18 , 4
	lfd   29,	280(%r11)
	fmadd  7 , 24 , 14 , 7
	fmadd  6 , 25 , 16 , 6
	dcbt 	%r21,%r11
	fnmsub 3 , 24 , 18 , 3
	fmadd  5 , 25 , 19 , 5
	fmadd  4 , 26 , 19 , 4
	lfd   24,	128(%r11)
	fmadd  7 , 25 , 17 , 7
	dcbt 	%r22,%r11
	fnmsub 6 , 26 , 17 , 6
	fmul   8 , 28 , 20
	lfd   25,	136(%r11)
	fnmsub 5 , 26 , 18 , 5
	fmul   9 , 28 , 21
	dcbt 	%r23,%r11
	fmadd  7 , 26 , 16 , 7
	fmadd  6 , 27 , 18 , 6
	lfd   28,	144(%r11)
	fmul   10 , 29 , 22
	lfd   26,	200(%r11)
	fmul   11 , 29 , 23
	dcbt 	%r24,%r11
	fmadd  8 , 24 , 12 , 8
	fmadd  7 , 27 , 19 , 7
	lfd   29,	152(%r11)
	fmadd  9 , 24 , 13 , 9
	fmadd  0 , 24 , 20 , 0
	lfd   27,	208(%r11)
	fmadd  1 , 24 , 21 , 1
	dcbt 	%r25,%r11
	fnmsub 8 , 25 , 13 , 8
	fmadd  10 , 26 , 12 , 10
	lfd   24,	216(%r11)
	fmadd  9 , 25 , 12 , 9
	fmadd  0 , 25 , 21 , 0
	dcbt 	%r26,%r11
	fnmsub 1 , 25 , 20 , 1
	fmadd  11 , 26 , 13 , 11
	fnmsub 10 , 27 , 13 , 10
	lfd   25,	224(%r11)
	fmadd  8 , 28 , 14 , 8
	fmadd  0 , 26 , 22 , 0
	lfd   13,	8(%r12)
	fmadd  1 , 26 , 23 , 1
	fmadd  11 , 27 , 12 , 11
	fmadd  9 , 28 , 15 , 9
	lfd   26,	160(%r11)
	fnmsub 8 , 29 , 15 , 8
	lfd   12,	0(%r12)
	fmadd  0 , 27 , 23 , 0
	fnmsub 1 , 27 , 22 , 1
	fmadd  2 , 28 , 20 , 2
	fmadd  9 , 29 , 14 , 9
	lfd   27,	168(%r11)
	fmadd  3 , 28 , 21 , 3
	fmadd  10 , 24 , 14 , 10
	stfd  0,	0(%r14)
	fmadd  11 , 24 , 15 , 11
	lfd   28,	232(%r11)
	fmadd  2 , 29 , 21 , 2
	stfd  1,	8(%r14)
	fmadd  8 , 26 , 16 , 8
	fnmsub 3 , 29 , 20 , 3
	fnmsub 10 , 25 , 15 , 10
	fmadd  11 , 25 , 14 , 11
	lfd   29,	240(%r11)
	fmadd  2 , 24 , 22 , 2
	lfd   15,	24(%r12)
	fmadd  9 , 26 , 17 , 9
	lfd   14,	16(%r12)
	fmadd  3 , 24 , 23 , 3
	fnmsub 8 , 27 , 17 , 8
	fmadd  4 , 26 , 20 , 4
	lfd   24,	176(%r11)
	fmadd  2 , 25 , 23 , 2
	fmadd  9 , 27 , 16 , 9
	fnmsub 3 , 25 , 22 , 3
	fmadd  5 , 26 , 21 , 5
	fmadd  4 , 27 , 21 , 4
	lfd   25,	184(%r11)
	fmadd  10 , 28 , 16 , 10
	lfd   26,	248(%r11)
	fmadd  11 , 28 , 17 , 11
	stfd  2,	16(%r14)
	fmadd  8 , 24 , 18 , 8
	stfd  3,	24(%r14)
	fnmsub 5 , 27 , 20 , 5
	fmadd  4 , 28 , 22 , 4
	fnmsub 10 , 29 , 17 , 10
	lfd   27,	256(%r11)
	fmadd  11 , 29 , 16 , 11
	fmadd  9 , 24 , 19 , 9
	lfd   17,	40(%r12)
	fmadd  5 , 28 , 23 , 5
	lfd   16,	32(%r12)
	fmadd  4 , 29 , 23 , 4
	fnmsub 8 , 25 , 19 , 8
	lfd   28,	264(%r11)
	fmadd  6 , 24 , 20 , 6
	fmadd  9 , 25 , 18 , 9
	fnmsub 5 , 29 , 22 , 5
	fmadd  7 , 24 , 21 , 7
	stfd  4,	32(%r14)
	fmadd  10 , 26 , 18 , 10
	lfd   29,	272(%r11)
	fmadd  6 , 25 , 21 , 6
	fmadd  11 , 26 , 19 , 11
	la   %r11,	288(%r11)
	fmadd  8 , 28 , 22 , 8
	stfd  5,	40(%r14)
	fnmsub 7 , 25 , 20 , 7
	lfd   24,	0(%r11)
	fnmsub 10 , 27 , 19 , 10
	fmadd  6 , 26 , 22 , 6
	lfd   25,	8(%r11)
	fmadd  11 , 27 , 18 , 11
	fmadd  9 , 28 , 23 , 9
	fmadd  7 , 26 , 23 , 7
	fmadd  10 , 28 , 20 , 10
	fmadd  6 , 27 , 23 , 6
	lfd   26,	16(%r11)
	fmadd  11 , 28 , 21 , 11
	fmadd  8 , 29 , 23 , 8
	fnmsub 7 , 27 , 22 , 7
	fnmsub 10 , 29 , 21 , 10
	fnmsub 9 , 29 , 22 , 9
	stfd  6,	48(%r14)
	fmadd  11 , 29 , 20 , 11
	lfd   27,	24(%r11)
	bdnz     lab1 
	stfd  7,	56(%r14)
	stfd  8,	64(%r14)
	stfd  9,	72(%r14)
	stfd  10,	80(%r14)
	stfd  11,	88(%r14)
lab0:
	la   %r1,	736(%r1)
	blr
