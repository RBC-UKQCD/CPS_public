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
.globl FtV1pV2Skip_asm
FtV1pV2Skip_asm:

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
	li   %r16,0
	li   %r17,32
	li   %r18,64
	bf gt,  FtV1pV2Skip_asm_lab0
	lfd   1,	0(%r10)
	lfd   13,	0(%r11)
	lfd   2,	8(%r10)
	lfd   14,	8(%r11)
	lfd   3,	16(%r10)
	lfd   15,	16(%r11)
	lfd   4,	24(%r10)
	lfd   16,	24(%r11)
	lfd   5,	32(%r10)
	lfd   17,	32(%r11)
	lfd   6,	40(%r10)
	lfd   18,	40(%r11)
	mtctr 	%r13
	b       FtV1pV2Skip_asm_lab1
.align 4
FtV1pV2Skip_asm_lab1:
	fmadd  25 , 0 , 1 , 13
	li   %r9,192
	fmadd  26 , 0 , 2 , 14
	lfd   7,	48(%r10)
	lfd   19,	48(%r11)
	add   %r9 , %r9 , %r11
	lfd   8,	56(%r10)
	fmadd  27 , 0 , 3 , 15
	dcbt 	%r16,%r9
	fmadd  28 , 0 , 4 , 16
	lfd   20,	56(%r11)
	fmadd  29 , 0 , 5 , 17
	lfd   9,	64(%r10)
	fmadd  30 , 0 , 6 , 18
	lfd   21,	64(%r11)
	lfd   10,	72(%r10)
	lfd   22,	72(%r11)
	dcbt 	%r17,%r9
	lfd   11,	80(%r10)
	lfd   23,	80(%r11)
	lfd   12,	88(%r10)
	lfd   24,	88(%r11)
	stfd  25,	0(%r12)
	dcbt 	%r18,%r9
	fmadd  25 , 0 , 7 , 19
	stfd  26,	8(%r12)
	stfd  27,	16(%r12)
	fmadd  26 , 0 , 8 , 20
	la   %r10,	192(%r10)
	fmadd  27 , 0 , 9 , 21
	stfd  28,	24(%r12)
	stfd  29,	32(%r12)
	fmadd  28 , 0 , 10 , 22
	la   %r11,	192(%r11)
	fmadd  29 , 0 , 11 , 23
	stfd  30,	40(%r12)
	stfd  25,	48(%r12)
	fmadd  30 , 0 , 12 , 24
	li   %r19,192
	stfd  26,	56(%r12)
	stfd  27,	64(%r12)
	add   %r19 , %r19 , %r10
	lfd   1,	0(%r10)
	stfd  28,	72(%r12)
	stfd  29,	80(%r12)
	dcbt 	%r16,%r19
	stfd  30,	88(%r12)
	lfd   13,	0(%r11)
	la   %r12,	192(%r12)
	lfd   2,	8(%r10)
	lfd   14,	8(%r11)
	dcbt 	%r17,%r19
	lfd   3,	16(%r10)
	lfd   15,	16(%r11)
	lfd   4,	24(%r10)
	lfd   16,	24(%r11)
	lfd   5,	32(%r10)
	dcbt 	%r18,%r19
	lfd   17,	32(%r11)
	lfd   6,	40(%r10)
	lfd   18,	40(%r11)
	la   %r13,	-1(%r13)
	bdnz    FtV1pV2Skip_asm_lab1
FtV1pV2Skip_asm_lab0:
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
