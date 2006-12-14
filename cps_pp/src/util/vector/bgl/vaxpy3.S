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
.globl vaxpy3
vaxpy3:

	or    %r12 , %r3 , %r3
	la   %r1,	-736(%r1)
	or    %r10 , %r5 , %r5
	or    %r11 , %r6 , %r6
	or    %r13 , %r7 , %r7
	stw   %r14,	368(%r1)
	or.  %r13 , %r13 , %r13
	stw   %r15,	376(%r1)
	or    %r14 , %r4 , %r4
	stw   %r16,	384(%r1)
	stw   %r17,	392(%r1)
	stw   %r18,	400(%r1)
	stw   %r19,	408(%r1)
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
	lfd   0,	0(%r14)
	li   %r16,32
	li   %r17,64
	bf gt,  vaxpy3_lab0
	lfd   1,	0(%r10)
	lfd   7,	0(%r11)
	lfd   2,	8(%r10)
	lfd   8,	8(%r11)
	lfd   3,	16(%r10)
	lfd   9,	16(%r11)
	lfd   4,	24(%r10)
	lfd   10,	24(%r11)
	mtctr 	%r13
	b       vaxpy3_lab1
.align 4
vaxpy3_lab1:
	addi  %r15 , %r12 , 0
	lfd   5,	32(%r10)
	fmadd  13 , 0 , 1 , 7
	lfd   11,	32(%r11)
	fmadd  14 , 0 , 2 , 8
	lfd   6,	40(%r10)
	fmadd  15 , 0 , 3 , 9
	lfd   12,	40(%r11)
	fmadd  16 , 0 , 4 , 10
	la   %r10,	48(%r10)
	la   %r11,	48(%r11)
	stfd  13,	0(%r15)
	stfd  14,	8(%r15)
	fmadd  17 , 0 , 5 , 11
	stfd  15,	16(%r15)
	stfd  16,	24(%r15)
	fmadd  18 , 0 , 6 , 12
	la   %r12,	48(%r12)
	dcbt 	%r16,%r10
	lfd   1,	0(%r10)
	stfd  17,	32(%r15)
	lfd   7,	0(%r11)
	stfd  18,	40(%r15)
	lfd   2,	8(%r10)
	dcbt 	%r17,%r10
	lfd   8,	8(%r11)
	lfd   3,	16(%r10)
	lfd   9,	16(%r11)
	lfd   4,	24(%r10)
	lfd   10,	24(%r11)
	dcbt 	%r16,%r11
	la   %r13,	-1(%r13)
	dcbt 	%r17,%r11
	bdnz    vaxpy3_lab1
vaxpy3_lab0:
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
