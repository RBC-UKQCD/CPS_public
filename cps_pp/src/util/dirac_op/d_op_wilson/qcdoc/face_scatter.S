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
.globl face_scatter
face_scatter:

	or    %r11 , %r3 , %r3
	la   %r1,	-864(%r1)
	or    %r12 , %r4 , %r4
	or    %r13 , %r5 , %r5
	stw   %r14,	368(%r1)
	stw   %r15,	376(%r1)
	or    %r14 , %r6 , %r6
	stw   %r16,	384(%r1)
	or.  %r14 , %r14 , %r14
	stw   %r17,	392(%r1)
	stw   %r18,	400(%r1)
	stw   %r19,	408(%r1)
	stw   %r20,	416(%r1)
	stw   %r21,	424(%r1)
	stw   %r22,	432(%r1)
	stw   %r23,	440(%r1)
	or    %r22 , %r11 , %r11
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
	li   %r15,0
	li   %r16,32
	li   %r17,64
	li   %r18,128
	li   %r19,160
	li   %r20,192
	dcbt 	%r15,%r12
	li   %r21,96
	lwz   %r10,	0(%r13)
	mulli  %r10 , %r10 , 128
	dcbt 	%r16,%r12
	add   %r11 , %r22 , %r10
	dcbt 	%r17,%r12
	bf gt,  face_scatter_lab0
	mtctr 	%r14
	b       face_scatter_lab1
.align 4
face_scatter_lab1:
	lfd   0,	0(%r12)
	lfd   1,	8(%r12)
	lfd   2,	16(%r12)
	lfd   3,	24(%r12)
	lfd   4,	32(%r12)
	lfd   5,	40(%r12)
	dcbt 	%r18,%r12
	lfd   6,	48(%r12)
	stfd  0,	0(%r11)
	lfd   7,	56(%r12)
	stfd  1,	8(%r11)
	lfd   8,	64(%r12)
	stfd  2,	16(%r11)
	lfd   9,	72(%r12)
	dcbt 	%r19,%r12
	lfd   10,	80(%r12)
	stfd  3,	24(%r11)
	lfd   11,	88(%r12)
	stfd  4,	32(%r11)
	la   %r13,	4(%r13)
	stfd  5,	40(%r11)
	lwz   %r10,	0(%r13)
	dcbt 	%r20,%r12
	mulli  %r10 , %r10 , 128
	stfd  6,	48(%r11)
	la   %r12,	128(%r12)
	stfd  7,	56(%r11)
	la   %r14,	-1(%r14)
	stfd  8,	64(%r11)
	stfd  9,	72(%r11)
	stfd  10,	80(%r11)
	stfd  11,	88(%r11)
	add   %r11 , %r22 , %r10
	bdnz    face_scatter_lab1
face_scatter_lab0:
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
	la   %r1,	864(%r1)
	blr
