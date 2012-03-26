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
.globl m1m2_lookup
m1m2_lookup:

	or    %r12 , %r3 , %r3
	la   %r1,	-736(%r1)
	or    %r13 , %r4 , %r4
	std   %r14,	368(%r1)
	std   %r15,	376(%r1)
	or    %r14 , %r5 , %r5
	std   %r16,	384(%r1)
	std   %r17,	392(%r1)
	or    %r16 , %r6 , %r6
	std   %r18,	400(%r1)
	or    %r17 , %r8 , %r8
	std   %r19,	408(%r1)
	or    %r18 , %r9 , %r9
	std   %r20,	416(%r1)
	or    %r19 , %r7 , %r7
	std   %r21,	424(%r1)
	or.  %r16 , %r16 , %r16
	std   %r22,	432(%r1)
	std   %r23,	440(%r1)
	std   %r24,	448(%r1)
	std   %r25,	456(%r1)
	std   %r26,	464(%r1)
	std   %r27,	472(%r1)
	std   %r28,	480(%r1)
	std   %r29,	488(%r1)
	std   %r30,	496(%r1)
	std   %r31,	504(%r1)
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
	bf gt,  m1m2_lookup_lab0
	or    %r20 , %r13 , %r13
	li   %r21,0
	or    %r26 , %r14 , %r14
	li   %r22,32
	or    %r27 , %r12 , %r12
	li   %r23,64
	li   %r24,96
	li   %r25,128
	lwz   %r11,	0(%r17)
	mtctr 	%r16
	mulli  %r11 , %r11 , 144
	add   %r13 , %r20 , %r11
	lfd   0,	0(%r13)
	lwz   %r11,	0(%r18)
	lfd   2,	16(%r13)
	mulli  %r11 , %r11 , 144
	lfd   1,	8(%r13)
	add   %r14 , %r26 , %r11
	lfd   3,	24(%r13)
	lfd   6,	0(%r14)
	lwz   %r11,	0(%r19)
	lfd   7,	8(%r14)
	mulli  %r11 , %r11 , 144
	lfd   8,	16(%r14)
	add   %r12 , %r27 , %r11
	lfd   9,	24(%r14)
	or    %r15 , %r12 , %r12
	lfd   10,	32(%r14)
	lfd   11,	40(%r14)
	b       m1m2_lookup_lab1
.align 4
m1m2_lookup_lab1:
	fmul   24 , 0 , 6
	stfd  26,	112(%r15)
	fmul   25 , 0 , 7
	lfd   12,	48(%r14)
	fmul   26 , 0 , 8
	stfd  27,	120(%r15)
	fmul   28 , 0 , 10
	lfd   13,	56(%r14)
	fmul   27 , 0 , 9
	lfd   14,	64(%r14)
	fmul   29 , 0 , 11
	lfd   15,	72(%r14)
	fnmsub 24 , 1 , 7 , 24
	lfd   16,	80(%r14)
	fmadd  25 , 1 , 6 , 25
	lfd   17,	88(%r14)
	fnmsub 26 , 1 , 9 , 26
	lfd   18,	96(%r14)
	fmadd  27 , 1 , 8 , 27
	lfd   19,	104(%r14)
	fnmsub 28 , 1 , 11 , 28
	lfd   20,	112(%r14)
	fmadd  29 , 1 , 10 , 29
	lfd   21,	120(%r14)
	fmadd  24 , 2 , 12 , 24
	lfd   22,	128(%r14)
	fmadd  25 , 2 , 13 , 25
	lfd   23,	136(%r14)
	fmadd  26 , 2 , 14 , 26
	lfd   4,	32(%r13)
	fmadd  27 , 2 , 15 , 27
	lfd   5,	40(%r13)
	fmadd  28 , 2 , 16 , 28
	stfd  30,	128(%r15)
	fmadd  29 , 2 , 17 , 29
	lwz   %r11,	4(%r17)
	fnmsub 24 , 3 , 13 , 24
	stfd  31,	136(%r15)
	fmadd  25 , 3 , 12 , 25
	lfd   0,	48(%r13)
	fnmsub 26 , 3 , 15 , 26
	mulli  %r11 , %r11 , 144
	fmadd  27 , 3 , 14 , 27
	lfd   1,	56(%r13)
	fnmsub 28 , 3 , 17 , 28
	lfd   2,	64(%r13)
	fmadd  29 , 3 , 16 , 29
	add   %r11 , %r20 , %r11
	fmadd  24 , 4 , 18 , 24
	or.  %r10 , %r10 , %r10
	fmadd  25 , 4 , 19 , 25
	lfd   3,	72(%r13)
	fmadd  26 , 4 , 20 , 26
	fmadd  27 , 4 , 21 , 27
	fmadd  30 , 4 , 22 , 28
	fmadd  31 , 4 , 23 , 29
	fnmsub 24 , 5 , 19 , 24
	fmadd  25 , 5 , 18 , 25
	lfd   4,	80(%r13)
	fnmsub 26 , 5 , 21 , 26
	fmadd  27 , 5 , 20 , 27
	fnmsub 30 , 5 , 23 , 30
	fmadd  31 , 5 , 22 , 31
	stfd  24,	0(%r12)
	stfd  25,	8(%r12)
	stfd  26,	16(%r12)
	stfd  27,	24(%r12)
	stfd  30,	32(%r12)
	stfd  31,	40(%r12)
	lfd   5,	88(%r13)
	bf lt,  m1m2_lookup_lab2
	or    %r11 , %r1 , %r1
m1m2_lookup_lab2:
	fmul   24 , 0 , 6
	dcbt 	%r21,%r11
	fmul   25 , 0 , 7
	or.  %r10 , %r10 , %r10
	fmul   26 , 0 , 8
	fmul   27 , 0 , 9
	fmul   28 , 0 , 10
	dcbt 	%r22,%r11
	fmul   29 , 0 , 11
	fnmsub 24 , 1 , 7 , 24
	fmadd  25 , 1 , 6 , 25
	fnmsub 26 , 1 , 9 , 26
	dcbt 	%r23,%r11
	fmadd  27 , 1 , 8 , 27
	fnmsub 28 , 1 , 11 , 28
	fmadd  29 , 1 , 10 , 29
	fmadd  24 , 2 , 12 , 24
	dcbt 	%r24,%r11
	fmadd  25 , 2 , 13 , 25
	fmadd  26 , 2 , 14 , 26
	fmadd  27 , 2 , 15 , 27
	fmadd  28 , 2 , 16 , 28
	dcbt 	%r25,%r11
	fmadd  29 , 2 , 17 , 29
	fnmsub 24 , 3 , 13 , 24
	fmadd  25 , 3 , 12 , 25
	lwz   %r11,	4(%r18)
	fnmsub 26 , 3 , 15 , 26
	dcbt 	%r22,%r17
	fmadd  27 , 3 , 14 , 27
	mulli  %r11 , %r11 , 144
	fnmsub 28 , 3 , 17 , 28
	add   %r11 , %r26 , %r11
	fmadd  29 , 3 , 16 , 29
	bf lt,  m1m2_lookup_lab3
	or    %r11 , %r1 , %r1
m1m2_lookup_lab3:
	fmadd  24 , 4 , 18 , 24
	dcbt 	%r21,%r11
	fmadd  25 , 4 , 19 , 25
	or    %r15 , %r12 , %r12
	fmadd  26 , 4 , 20 , 26
	fmadd  27 , 4 , 21 , 27
	lfd   0,	96(%r13)
	fmadd  30 , 4 , 22 , 28
	dcbt 	%r22,%r11
	fmadd  31 , 4 , 23 , 29
	fnmsub 24 , 5 , 19 , 24
	fmadd  25 , 5 , 18 , 25
	lfd   1,	104(%r13)
	fnmsub 26 , 5 , 21 , 26
	dcbt 	%r23,%r11
	fmadd  27 , 5 , 20 , 27
	fnmsub 30 , 5 , 23 , 30
	fmadd  31 , 5 , 22 , 31
	stfd  24,	48(%r12)
	fmul   28 , 0 , 10
	dcbt 	%r24,%r11
	fmul   24 , 0 , 6
	fmul   29 , 0 , 11
	stfd  25,	56(%r12)
	dcbt 	%r25,%r11
	fmul   25 , 0 , 7
	fnmsub 24 , 1 , 7 , 24
	fnmsub 28 , 1 , 11 , 28
	stfd  26,	64(%r12)
	fmadd  29 , 1 , 10 , 29
	dcbt 	%r22,%r18
	fmul   26 , 0 , 8
	fmadd  25 , 1 , 6 , 25
	la   %r18,	4(%r18)
	stfd  27,	72(%r12)
	lwz   %r11,	0(%r18)
	fmul   27 , 0 , 9
	lfd   2,	112(%r13)
	mulli  %r11 , %r11 , 144
	lfd   3,	120(%r13)
	add   %r14 , %r26 , %r11
	la   %r17,	4(%r17)
	fnmsub 26 , 1 , 9 , 26
	lfd   4,	128(%r13)
	lfd   5,	136(%r13)
	fmadd  27 , 1 , 8 , 27
	lwz   %r11,	0(%r17)
	stfd  30,	80(%r12)
	mulli  %r11 , %r11 , 144
	stfd  31,	88(%r12)
	add   %r13 , %r20 , %r11
	la   %r19,	4(%r19)
	fmadd  24 , 2 , 12 , 24
	lfd   0,	0(%r13)
	fmadd  25 , 2 , 13 , 25
	lwz   %r11,	0(%r19)
	fmadd  26 , 2 , 14 , 26
	lfd   1,	8(%r13)
	fmadd  27 , 2 , 15 , 27
	mulli  %r11 , %r11 , 144
	fmadd  28 , 2 , 16 , 28
	add   %r12 , %r27 , %r11
	fmadd  29 , 2 , 17 , 29
	lfd   6,	0(%r14)
	fnmsub 24 , 3 , 13 , 24
	lfd   7,	8(%r14)
	fmadd  25 , 3 , 12 , 25
	lfd   2,	16(%r13)
	fnmsub 26 , 3 , 15 , 26
	la   %r16,	-1(%r16)
	fmadd  27 , 3 , 14 , 27
	lfd   8,	16(%r14)
	fnmsub 28 , 3 , 17 , 28
	lfd   9,	24(%r14)
	fmadd  29 , 3 , 16 , 29
	fmadd  24 , 4 , 18 , 24
	lfd   10,	32(%r14)
	fmadd  25 , 4 , 19 , 25
	lfd   3,	24(%r13)
	fmadd  26 , 4 , 20 , 26
	fmadd  27 , 4 , 21 , 27
	lfd   11,	40(%r14)
	fmadd  30 , 4 , 22 , 28
	fmadd  31 , 4 , 23 , 29
	fnmsub 24 , 5 , 19 , 24
	fmadd  25 , 5 , 18 , 25
	fnmsub 26 , 5 , 21 , 26
	fmadd  27 , 5 , 20 , 27
	fnmsub 30 , 5 , 23 , 30
	fmadd  31 , 5 , 22 , 31
	stfd  24,	96(%r15)
	stfd  25,	104(%r15)
	bdnz    m1m2_lookup_lab1
	stfd  26,	112(%r15)
	stfd  27,	120(%r15)
	stfd  30,	128(%r15)
	stfd  31,	136(%r15)
m1m2_lookup_lab0:
	ld   %r14,	368(%r1)
	ld   %r15,	376(%r1)
	ld   %r16,	384(%r1)
	ld   %r17,	392(%r1)
	ld   %r18,	400(%r1)
	ld   %r19,	408(%r1)
	ld   %r20,	416(%r1)
	ld   %r21,	424(%r1)
	ld   %r22,	432(%r1)
	ld   %r23,	440(%r1)
	ld   %r24,	448(%r1)
	ld   %r25,	456(%r1)
	ld   %r26,	464(%r1)
	ld   %r27,	472(%r1)
	ld   %r28,	480(%r1)
	ld   %r29,	488(%r1)
	ld   %r30,	496(%r1)
	ld   %r31,	504(%r1)
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
