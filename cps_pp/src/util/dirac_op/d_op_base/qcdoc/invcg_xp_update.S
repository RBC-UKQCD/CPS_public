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
.globl invcg_xp_update
invcg_xp_update:

	or    %r12 , %r7 , %r7
	la   %r1,	-736(%r1)
	or    %r13 , %r8 , %r8
	stw   %r14,	368(%r1)
	stw   %r15,	376(%r1)
	or    %r14 , %r3 , %r3
	stw   %r16,	384(%r1)
	or    %r15 , %r4 , %r4
	stw   %r17,	392(%r1)
	or    %r16 , %r9 , %r9
	stw   %r18,	400(%r1)
	or    %r17 , %r5 , %r5
	stw   %r19,	408(%r1)
	or    %r18 , %r6 , %r6
	stw   %r20,	416(%r1)
	or.  %r16 , %r16 , %r16
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
	lfd   0,	0(%r17)
	lfd   1,	0(%r18)
	li   %r21,96
	li   %r22,128
	li   %r23,160
	li   %r24,192
	li   %r25,224
	li   %r26,256
	li   %r27,288
	li   %r28,320
	li   %r29,352
	li   %r30,0
	li   %r31,32
	li   %r3,64
	bf gt,  invcg_xp_update_lab0
	lfd   2,	0(%r12)
	lfd   8,	0(%r13)
	lfd   14,	96(%r13)
	lfd   3,	8(%r12)
	lfd   9,	8(%r13)
	lfd   15,	104(%r13)
	lfd   4,	16(%r12)
	lfd   10,	16(%r13)
	lfd   16,	112(%r13)
	lfd   5,	24(%r12)
	lfd   11,	24(%r13)
	lfd   17,	120(%r13)
	mtctr 	%r16
	b       invcg_xp_update_lab1
.align 4
invcg_xp_update_lab1:
	addi  %r19 , %r14 , 0
	addi  %r20 , %r15 , 0
	fmadd  20 , 0 , 2 , 8
	lfd   6,	32(%r12)
	fmadd  26 , 1 , 2 , 14
	lfd   12,	32(%r13)
	fmadd  21 , 0 , 3 , 9
	lfd   18,	128(%r13)
	fmadd  27 , 1 , 3 , 15
	lfd   7,	40(%r12)
	fmadd  22 , 0 , 4 , 10
	lfd   13,	40(%r13)
	fmadd  28 , 1 , 4 , 16
	stfd  20,	0(%r19)
	fmadd  23 , 0 , 5 , 11
	lfd   19,	136(%r13)
	fmadd  29 , 1 , 5 , 17
	stfd  26,	0(%r20)
	fmadd  24 , 0 , 6 , 12
	lfd   2,	48(%r12)
	fmadd  30 , 1 , 6 , 18
	stfd  21,	8(%r19)
	lfd   8,	48(%r13)
	fmadd  25 , 0 , 7 , 13
	stfd  27,	8(%r20)
	lfd   14,	144(%r13)
	fmadd  31 , 1 , 7 , 19
	stfd  22,	16(%r19)
	lfd   3,	56(%r12)
	stfd  28,	16(%r20)
	lfd   9,	56(%r13)
	fmadd  20 , 0 , 2 , 8
	stfd  23,	24(%r19)
	lfd   15,	152(%r13)
	fmadd  26 , 1 , 2 , 14
	stfd  29,	24(%r20)
	lfd   4,	64(%r12)
	stfd  24,	32(%r19)
	lfd   10,	64(%r13)
	fmadd  21 , 0 , 3 , 9
	stfd  30,	32(%r20)
	lfd   16,	160(%r13)
	fmadd  27 , 1 , 3 , 15
	stfd  25,	40(%r19)
	lfd   5,	72(%r12)
	stfd  31,	40(%r20)
	lfd   11,	72(%r13)
	fmadd  22 , 0 , 4 , 10
	lfd   17,	168(%r13)
	stfd  20,	48(%r19)
	fmadd  28 , 1 , 4 , 16
	lfd   6,	80(%r12)
	stfd  26,	48(%r20)
	lfd   12,	80(%r13)
	stfd  21,	56(%r19)
	fmadd  23 , 0 , 5 , 11
	lfd   18,	176(%r13)
	fmadd  29 , 1 , 5 , 17
	stfd  27,	56(%r20)
	lfd   7,	88(%r12)
	stfd  22,	64(%r19)
	lfd   13,	88(%r13)
	fmadd  24 , 0 , 6 , 12
	stfd  28,	64(%r20)
	lfd   19,	184(%r13)
	fmadd  30 , 1 , 6 , 18
	stfd  23,	72(%r19)
	la   %r12,	96(%r12)
	stfd  29,	72(%r20)
	dcbt 	%r21,%r12
	fmadd  25 , 0 , 7 , 13
	stfd  24,	80(%r19)
	la   %r13,	192(%r13)
	fmadd  31 , 1 , 7 , 19
	stfd  30,	80(%r20)
	la   %r14,	192(%r14)
	la   %r15,	96(%r15)
	stfd  25,	88(%r19)
	dcbt 	%r22,%r12
	stfd  31,	88(%r20)
	lfd   2,	0(%r12)
	lfd   8,	0(%r13)
	lfd   14,	96(%r13)
	lfd   3,	8(%r12)
	dcbt 	%r23,%r12
	lfd   9,	8(%r13)
	lfd   15,	104(%r13)
	lfd   4,	16(%r12)
	lfd   10,	16(%r13)
	lfd   16,	112(%r13)
	dcbt 	%r24,%r13
	lfd   5,	24(%r12)
	lfd   11,	24(%r13)
	lfd   17,	120(%r13)
	la   %r16,	-1(%r16)
	dcbt 	%r25,%r13
	dcbt 	%r26,%r13
	dcbt 	%r27,%r13
	dcbt 	%r28,%r13
	dcbt 	%r29,%r13
	bdnz    invcg_xp_update_lab1
invcg_xp_update_lab0:
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
