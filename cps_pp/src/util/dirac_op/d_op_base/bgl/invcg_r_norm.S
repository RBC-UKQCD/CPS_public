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
.globl invcg_r_norm
invcg_r_norm:

	or    %r13 , %r3 , %r3
	la   %r1,	-736(%r1)
	or    %r11 , %r5 , %r5
	or    %r12 , %r6 , %r6
	stw   %r14,	368(%r1)
	stw   %r15,	376(%r1)
	or    %r14 , %r7 , %r7
	stw   %r16,	384(%r1)
	or    %r15 , %r4 , %r4
	stw   %r17,	392(%r1)
	stw   %r18,	400(%r1)
	or    %r17 , %r8 , %r8
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
	lfd   0,	0(%r15)
	li   %r18,0
	stw   %r18,	0(%r17)
	stw   %r18,	4(%r17)
invcg_r_norm_lab0:
	or.  %r14 , %r14 , %r14
	lfd   19,	0(%r17)
	lfd   20,	0(%r17)
	lfd   21,	0(%r17)
	lfd   22,	0(%r17)
	lfd   23,	0(%r17)
	lfd   24,	0(%r17)
	li   %r19,96
	li   %r20,128
	li   %r21,160
	li   %r22,0
	li   %r23,32
	li   %r24,64
	bf gt,  invcg_r_norm_lab1
	lfd   1,	0(%r11)
	lfd   7,	0(%r12)
	lfd   2,	8(%r11)
	lfd   8,	8(%r12)
	lfd   3,	16(%r11)
	lfd   9,	16(%r12)
	lfd   4,	24(%r11)
	lfd   10,	24(%r12)
	mtctr 	%r14
	b       invcg_r_norm_lab2
.align 4
invcg_r_norm_lab2:
	addi  %r16 , %r13 , 0
	lfd   5,	32(%r11)
	fmadd  13 , 0 , 1 , 7
	lfd   11,	32(%r12)
	fmadd  14 , 0 , 2 , 8
	lfd   6,	40(%r11)
	fmadd  15 , 0 , 3 , 9
	lfd   12,	40(%r12)
	fmadd  16 , 0 , 4 , 10
	lfd   1,	48(%r11)
	lfd   7,	48(%r12)
	fmadd  19 , 13 , 13 , 19
	stfd  13,	0(%r16)
	fmadd  20 , 14 , 14 , 20
	stfd  14,	8(%r16)
	fmadd  21 , 15 , 15 , 21
	stfd  15,	16(%r16)
	fmadd  22 , 16 , 16 , 22
	stfd  16,	24(%r16)
	fmadd  17 , 0 , 5 , 11
	lfd   2,	56(%r11)
	fmadd  18 , 0 , 6 , 12
	lfd   8,	56(%r12)
	fmadd  13 , 0 , 1 , 7
	lfd   3,	64(%r11)
	lfd   9,	64(%r12)
	lfd   4,	72(%r11)
	fmadd  23 , 17 , 17 , 23
	stfd  17,	32(%r16)
	fmadd  24 , 18 , 18 , 24
	stfd  18,	40(%r16)
	fmadd  19 , 13 , 13 , 19
	lfd   10,	72(%r12)
	fmadd  14 , 0 , 2 , 8
	stfd  13,	48(%r16)
	lfd   5,	80(%r11)
	fmadd  15 , 0 , 3 , 9
	lfd   11,	80(%r12)
	lfd   6,	88(%r11)
	lfd   12,	88(%r12)
	fmadd  20 , 14 , 14 , 20
	stfd  14,	56(%r16)
	fmadd  16 , 0 , 4 , 10
	la   %r11,	96(%r11)
	fmadd  21 , 15 , 15 , 21
	stfd  15,	64(%r16)
	la   %r12,	192(%r12)
	fmadd  17 , 0 , 5 , 11
	la   %r13,	192(%r13)
	dcbt 	%r19,%r11
	fmadd  22 , 16 , 16 , 22
	stfd  16,	72(%r16)
	fmadd  18 , 0 , 6 , 12
	li   %r10,192
	lfd   1,	0(%r11)
	fmadd  23 , 17 , 17 , 23
	stfd  17,	80(%r16)
	mulli  %r10 , %r10 , 1
	lfd   7,	0(%r12)
	add   %r10 , %r10 , %r12
	dcbt 	%r20,%r11
	fmadd  24 , 18 , 18 , 24
	stfd  18,	88(%r16)
	lfd   2,	8(%r11)
	lfd   8,	8(%r12)
	lfd   3,	16(%r11)
	lfd   9,	16(%r12)
	dcbt 	%r21,%r11
	lfd   4,	24(%r11)
	lfd   10,	24(%r12)
	la   %r14,	-1(%r14)
	dcbt 	%r22,%r10
	dcbt 	%r23,%r10
	dcbt 	%r24,%r10
	bdnz    invcg_r_norm_lab2
	fadd   19 , 19 , 20
	fadd   21 , 21 , 22
	fadd   23 , 23 , 24
	fadd   19 , 19 , 21
	fadd   19 , 19 , 23
	stfd  19,	0(%r17)
invcg_r_norm_lab1:
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
