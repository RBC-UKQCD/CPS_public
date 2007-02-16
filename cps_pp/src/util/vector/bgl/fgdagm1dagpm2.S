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
.globl fgdagm1dagpm2
fgdagm1dagpm2:

	or    %r11 , %r3 , %r3
	la   %r1,	-736(%r1)
	or    %r12 , %r5 , %r5
	or    %r13 , %r6 , %r6
	stw   %r14,	368(%r1)
	stw   %r15,	376(%r1)
	stw   %r16,	384(%r1)
	or    %r15 , %r7 , %r7
	stw   %r17,	392(%r1)
	or    %r16 , %r4 , %r4
	stw   %r18,	400(%r1)
	stw   %r19,	408(%r1)
	or    %r18 , %r8 , %r8
	stw   %r20,	416(%r1)
	stw   %r21,	424(%r1)
	stw   %r22,	432(%r1)
	stw   %r23,	440(%r1)
	stw   %r24,	448(%r1)
	stw   %r25,	456(%r1)
	stw   %r26,	464(%r1)
	stw   %r27,	472(%r1)
	stw   %r28,	480(%r1)
	stw   %r29,	488(%r1)
	stw   %r30,	496(%r1)
	stw   %r31,	504(%r1)
	stfd  14,	112(%r1)
	stfd  15,	120(%r1)
	stfd  16,	128(%r1)
	stfd  17,	136(%r1)
	stfd  18,	144(%r1)
	stfd  19,	152(%r1)
	stfd  20,	160(%r1)
	stfd  21,	168(%r1)
	stfd  22,	176(%r1)
	stfd  23,	184(%r1)
	stfd  24,	192(%r1)
	stfd  25,	200(%r1)
	stfd  26,	208(%r1)
	stfd  27,	216(%r1)
	stfd  28,	224(%r1)
	stfd  29,	232(%r1)
	stfd  30,	240(%r1)
	stfd  31,	248(%r1)
	lwz   %r17,	0(%r18)
	lfd   30,	0(%r16)
	or.  %r17 , %r17 , %r17
	bf gt,  fgdagm1dagpm2_lab0
	or    %r14 , %r11 , %r11
	li   %r19,0
	li   %r20,32
	li   %r21,64
	li   %r22,96
	li   %r23,128
	li   %r24,160
	li   %r25,192
	li   %r26,224
	li   %r27,256
	lfd   0,	0(%r12)
	lfd   1,	8(%r12)
	lfd   2,	48(%r12)
	lfd   3,	56(%r12)
	lfd   24,	0(%r15)
	lfd   25,	8(%r15)
	lfd   6,	0(%r13)
	fmul   0 , 30 , 0
	lfd   7,	8(%r13)
	fmul   1 , 30 , 1
	lfd   8,	48(%r13)
	fmul   2 , 30 , 2
	lfd   9,	56(%r13)
	fmul   3 , 30 , 3
	lfd   10,	96(%r13)
	lfd   11,	104(%r13)
	mtctr 	%r17
	b       fgdagm1dagpm2_lab1
.align 4
fgdagm1dagpm2_lab1:
	fmadd  24 , 0 , 6 , 24
	lfd   26,	16(%r15)
	fnmsub 25 , 0 , 7 , 25
	lfd   27,	24(%r15)
	or    %r14 , %r11 , %r11
	lfd   12,	16(%r13)
	lfd   13,	24(%r13)
	lfd   14,	64(%r13)
	fnmsub 24 , 1 , 7 , 24
	lfd   15,	72(%r13)
	fnmsub 25 , 1 , 6 , 25
	lfd   16,	112(%r13)
	fmadd  26 , 0 , 8 , 26
	lfd   17,	120(%r13)
	fnmsub 27 , 0 , 9 , 27
	lfd   18,	32(%r13)
	lfd   19,	40(%r13)
	fmadd  24 , 2 , 12 , 24
	lfd   20,	80(%r13)
	fnmsub 25 , 2 , 13 , 25
	lfd   21,	88(%r13)
	fnmsub 26 , 1 , 9 , 26
	lfd   22,	128(%r13)
	fnmsub 27 , 1 , 8 , 27
	lfd   23,	136(%r13)
	lfd   4,	96(%r12)
	fnmsub 24 , 3 , 13 , 24
	lfd   5,	104(%r12)
	fnmsub 25 , 3 , 12 , 25
	lfd   28,	32(%r15)
	fmadd  26 , 2 , 14 , 26
	lfd   29,	40(%r15)
	fnmsub 27 , 2 , 15 , 27
	li   %r10,576
	la   %r17,	-1(%r17)
	mulli  %r10 , %r10 , 1
	fmul   4 , 30 , 4
	add   %r10 , %r10 , %r12
	fmul   5 , 30 , 5
	dcbt 	%r19,%r10
	fmadd  28 , 0 , 10 , 28
	fnmsub 29 , 0 , 11 , 29
	fnmsub 26 , 3 , 15 , 26
	fnmsub 27 , 3 , 14 , 27
	lfd   0,	16(%r12)
	fmadd  24 , 4 , 18 , 24
	dcbt 	%r20,%r10
	fnmsub 28 , 1 , 11 , 28
	fnmsub 29 , 1 , 10 , 29
	fnmsub 25 , 4 , 19 , 25
	fmadd  26 , 4 , 20 , 26
	lfd   1,	24(%r12)
	fnmsub 27 , 4 , 21 , 27
	dcbt 	%r21,%r10
	fmadd  28 , 2 , 16 , 28
	fnmsub 29 , 2 , 17 , 29
	fnmsub 24 , 5 , 19 , 24
	fnmsub 25 , 5 , 18 , 25
	lfd   2,	64(%r12)
	fnmsub 26 , 5 , 21 , 26
	dcbt 	%r22,%r10
	fnmsub 28 , 3 , 17 , 28
	fnmsub 29 , 3 , 16 , 29
	fnmsub 27 , 5 , 20 , 27
	stfd  24,	0(%r11)
	fmul   0 , 30 , 0
	lfd   3,	72(%r12)
	fmul   1 , 30 , 1
	fmadd  28 , 4 , 22 , 28
	stfd  25,	8(%r11)
	fnmsub 29 , 4 , 23 , 29
	stfd  26,	16(%r11)
	fmul   2 , 30 , 2
	lfd   4,	112(%r12)
	stfd  27,	24(%r11)
	fnmsub 28 , 5 , 23 , 28
	fnmsub 29 , 5 , 22 , 29
	lfd   24,	48(%r15)
	fmul   3 , 30 , 3
	lfd   25,	56(%r15)
	lfd   5,	120(%r12)
	fmul   4 , 30 , 4
	stfd  28,	32(%r11)
	stfd  29,	40(%r11)
	fmadd  24 , 0 , 6 , 24
	lfd   26,	64(%r15)
	fnmsub 25 , 0 , 7 , 25
	lfd   27,	72(%r15)
	fmul   5 , 30 , 5
	lfd   28,	80(%r15)
	lfd   29,	88(%r15)
	fnmsub 24 , 1 , 7 , 24
	dcbt 	%r23,%r10
	fnmsub 25 , 1 , 6 , 25
	fmadd  26 , 0 , 8 , 26
	fnmsub 27 , 0 , 9 , 27
	dcbt 	%r23,%r13
	fmadd  28 , 0 , 10 , 28
	fnmsub 29 , 0 , 11 , 29
	fnmsub 26 , 1 , 9 , 26
	fnmsub 27 , 1 , 8 , 27
	dcbt 	%r24,%r13
	fmadd  24 , 2 , 12 , 24
	fnmsub 28 , 1 , 11 , 28
	fnmsub 29 , 1 , 10 , 29
	lfd   0,	32(%r12)
	fnmsub 25 , 2 , 13 , 25
	dcbt 	%r25,%r13
	fmadd  26 , 2 , 14 , 26
	fnmsub 27 , 2 , 15 , 27
	fmadd  28 , 2 , 16 , 28
	lfd   1,	40(%r12)
	fnmsub 29 , 2 , 17 , 29
	dcbt 	%r26,%r13
	fnmsub 24 , 3 , 13 , 24
	fnmsub 25 , 3 , 12 , 25
	fnmsub 26 , 3 , 15 , 26
	lfd   2,	80(%r12)
	fnmsub 27 , 3 , 14 , 27
	dcbt 	%r27,%r13
	fnmsub 28 , 3 , 17 , 28
	fnmsub 29 , 3 , 16 , 29
	fmadd  24 , 4 , 18 , 24
	la   %r13,	144(%r13)
	fnmsub 25 , 4 , 19 , 25
	dcbt 	%r23,%r15
	fmadd  26 , 4 , 20 , 26
	fnmsub 27 , 4 , 21 , 27
	fmadd  28 , 4 , 22 , 28
	lfd   3,	88(%r12)
	fnmsub 29 , 4 , 23 , 29
	dcbt 	%r24,%r15
	fnmsub 24 , 5 , 19 , 24
	fnmsub 25 , 5 , 18 , 25
	fnmsub 26 , 5 , 21 , 26
	lfd   4,	128(%r12)
	fnmsub 27 , 5 , 20 , 27
	dcbt 	%r25,%r15
	fnmsub 28 , 5 , 23 , 28
	fnmsub 29 , 5 , 22 , 29
	fmul   0 , 30 , 0
	stfd  24,	48(%r11)
	fmul   1 , 30 , 1
	dcbt 	%r26,%r15
	fmul   2 , 30 , 2
	fmul   3 , 30 , 3
	fmul   4 , 30 , 4
	lfd   5,	136(%r12)
	dcbt 	%r27,%r15
	stfd  25,	56(%r11)
	stfd  26,	64(%r11)
	la   %r12,	576(%r12)
	fmul   5 , 30 , 5
	stfd  27,	72(%r11)
	stfd  28,	80(%r11)
	stfd  29,	88(%r11)
	lfd   24,	96(%r15)
	lfd   25,	104(%r15)
	lfd   26,	112(%r15)
	lfd   27,	120(%r15)
	lfd   28,	128(%r15)
	fmadd  24 , 0 , 6 , 24
	lfd   29,	136(%r15)
	fnmsub 25 , 0 , 7 , 25
	fmadd  26 , 0 , 8 , 26
	la   %r15,	144(%r15)
	fnmsub 27 , 0 , 9 , 27
	fmadd  28 , 0 , 10 , 28
	fnmsub 24 , 1 , 7 , 24
	fnmsub 29 , 0 , 11 , 29
	fnmsub 25 , 1 , 6 , 25
	lfd   7,	8(%r13)
	fnmsub 26 , 1 , 9 , 26
	lfd   0,	0(%r12)
	fnmsub 27 , 1 , 8 , 27
	lfd   6,	0(%r13)
	fnmsub 28 , 1 , 11 , 28
	fnmsub 29 , 1 , 10 , 29
	lfd   8,	48(%r13)
	fmadd  24 , 2 , 12 , 24
	lfd   9,	56(%r13)
	fnmsub 25 , 2 , 13 , 25
	lfd   1,	8(%r12)
	fmadd  26 , 2 , 14 , 26
	fnmsub 27 , 2 , 15 , 27
	lfd   10,	96(%r13)
	fmadd  28 , 2 , 16 , 28
	lfd   11,	104(%r13)
	fnmsub 29 , 2 , 17 , 29
	fnmsub 24 , 3 , 13 , 24
	fnmsub 25 , 3 , 12 , 25
	lfd   2,	48(%r12)
	fnmsub 26 , 3 , 15 , 26
	fnmsub 27 , 3 , 14 , 27
	fnmsub 28 , 3 , 17 , 28
	fnmsub 29 , 3 , 16 , 29
	fmadd  24 , 4 , 18 , 24
	fnmsub 25 , 4 , 19 , 25
	lfd   3,	56(%r12)
	fmadd  26 , 4 , 20 , 26
	fnmsub 27 , 4 , 21 , 27
	fmadd  28 , 4 , 22 , 28
	fnmsub 29 , 4 , 23 , 29
	fnmsub 24 , 5 , 19 , 24
	fnmsub 25 , 5 , 18 , 25
	fnmsub 26 , 5 , 21 , 26
	fnmsub 27 , 5 , 20 , 27
	fnmsub 28 , 5 , 23 , 28
	fnmsub 29 , 5 , 22 , 29
	stfd  24,	96(%r11)
	fmul   0 , 30 , 0
	stfd  25,	104(%r11)
	fmul   1 , 30 , 1
	fmul   2 , 30 , 2
	stfd  26,	112(%r11)
	fmul   3 , 30 , 3
	stfd  27,	120(%r11)
	stfd  28,	128(%r11)
	stfd  29,	136(%r11)
	la   %r11,	144(%r11)
	lfd   24,	0(%r15)
	lfd   25,	8(%r15)
	bdnz    fgdagm1dagpm2_lab1
	stfd  26,	112(%r14)
	stfd  27,	120(%r14)
	stfd  28,	128(%r14)
	stfd  29,	136(%r14)
fgdagm1dagpm2_lab0:
	lwz   %r14,	368(%r1)
	lwz   %r15,	376(%r1)
	lwz   %r16,	384(%r1)
	lwz   %r17,	392(%r1)
	lwz   %r18,	400(%r1)
	lwz   %r19,	408(%r1)
	lwz   %r20,	416(%r1)
	lwz   %r21,	424(%r1)
	lwz   %r22,	432(%r1)
	lwz   %r23,	440(%r1)
	lwz   %r24,	448(%r1)
	lwz   %r25,	456(%r1)
	lwz   %r26,	464(%r1)
	lwz   %r27,	472(%r1)
	lwz   %r28,	480(%r1)
	lwz   %r29,	488(%r1)
	lwz   %r30,	496(%r1)
	lwz   %r31,	504(%r1)
	lfd   14,	112(%r1)
	lfd   15,	120(%r1)
	lfd   16,	128(%r1)
	lfd   17,	136(%r1)
	lfd   18,	144(%r1)
	lfd   19,	152(%r1)
	lfd   20,	160(%r1)
	lfd   21,	168(%r1)
	lfd   22,	176(%r1)
	lfd   23,	184(%r1)
	lfd   24,	192(%r1)
	lfd   25,	200(%r1)
	lfd   26,	208(%r1)
	lfd   27,	216(%r1)
	lfd   28,	224(%r1)
	lfd   29,	232(%r1)
	lfd   30,	240(%r1)
	lfd   31,	248(%r1)
	la   %r1,	736(%r1)
	blr
