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
.globl copy_gauge
copy_gauge:

	or    %r10 , %r3 , %r3
	la   %r1,	-736(%r1)
	or    %r9 , %r4 , %r4
	or    %r12 , %r5 , %r5
	or    %r13 , %r6 , %r6
	std   %r14,	368(%r1)
	std   %r15,	376(%r1)
	std   %r16,	384(%r1)
	std   %r17,	392(%r1)
	std   %r18,	400(%r1)
	std   %r19,	408(%r1)
	std   %r20,	416(%r1)
	std   %r21,	424(%r1)
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
	lwz   %r11,	0(%r12)
	or.  %r11 , %r11 , %r11
	bf gt,  copy_gauge_lab0
	or    %r14 , %r10 , %r10
	li   %r15,0
	li   %r16,32
	li   %r17,64
	li   %r18,96
	li   %r19,128
	ld   %r8,	0(%r13)
	li   %r20,160
	mulli  %r8 , %r8 , 144
	li   %r21,192
	add   %r10 , %r14 , %r8
	li   %r22,224
	li   %r23,256
	mtctr 	%r11
	b       copy_gauge_lab1
.align 4
copy_gauge_lab1:
	lfd   0,	8(%r9)
	lfd   1,	16(%r9)
	lfd   2,	24(%r9)
	lfd   3,	32(%r9)
	lfd   4,	40(%r9)
	lfd   5,	48(%r9)
	dcbt 	%r19,%r9
	stfd  0,	0(%r10)
	dcbt 	%r20,%r9
	stfd  1,	8(%r10)
	dcbt 	%r21,%r9
	stfd  2,	16(%r10)
	dcbt 	%r22,%r9
	stfd  3,	24(%r10)
	dcbt 	%r23,%r9
	stfd  4,	32(%r10)
	lfd   6,	56(%r9)
	stfd  5,	40(%r10)
	lfd   7,	64(%r9)
	lfd   8,	72(%r9)
	lfd   9,	80(%r9)
	lfd   10,	88(%r9)
	lfd   11,	96(%r9)
	stfd  6,	48(%r10)
	lfd   12,	104(%r9)
	stfd  7,	56(%r10)
	lfd   13,	112(%r9)
	stfd  8,	64(%r10)
	lfd   14,	120(%r9)
	stfd  9,	72(%r10)
	lfd   15,	128(%r9)
	stfd  10,	80(%r10)
	lfd   16,	136(%r9)
	stfd  11,	88(%r10)
	lfd   17,	144(%r9)
	stfd  12,	96(%r10)
	la   %r13,	8(%r13)
	stfd  13,	104(%r10)
	ld   %r8,	0(%r13)
	stfd  14,	112(%r10)
	mulli  %r8 , %r8 , 144
	la   %r9,	152(%r9)
	stfd  15,	120(%r10)
	la   %r11,	-1(%r11)
	stfd  16,	128(%r10)
	stfd  17,	136(%r10)
	add   %r10 , %r14 , %r8
	bdnz    copy_gauge_lab1
copy_gauge_lab0:
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
