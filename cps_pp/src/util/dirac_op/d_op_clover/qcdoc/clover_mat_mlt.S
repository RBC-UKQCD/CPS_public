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
.globl clover_mat_mlt_asm
clover_mat_mlt_asm:

	or    %r12 , %r3 , %r3
	la   %r1,	-736(%r1)
	or    %r10 , %r4 , %r4
	or    %r11 , %r5 , %r5
	or    %r9 , %r6 , %r6
	stw   %r14,	368(%r1)
	add   %r9 , %r9 , %r9
	stw   %r15,	376(%r1)
	or.  %r9 , %r9 , %r9
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
	bf gt,  clover_mat_mlt_asm_lab0
	addi  %r13 , %r12 , 0
	li   %r14,96
	li   %r15,128
	li   %r16,160
	li   %r17,288
	li   %r18,320
	li   %r19,352
	li   %r20,384
	li   %r21,416
	li   %r22,448
	li   %r23,480
	li   %r24,512
	li   %r25,544
	lfd   12,	0(%r11)
	lfd   13,	8(%r11)
	lfd   14,	16(%r11)
	lfd   15,	24(%r11)
	lfd   16,	32(%r11)
	lfd   17,	40(%r11)
	lfd   18,	48(%r11)
	lfd   24,	0(%r10)
	lfd   25,	8(%r10)
	lfd   26,	16(%r10)
	lfd   27,	24(%r10)
	mtctr 	%r9
	b       clover_mat_mlt_asm_lab1
.align 4
clover_mat_mlt_asm_lab1:
	fmul   0 , 24 , 12
	lfd   28,	32(%r10)
	fmul   1 , 24 , 13
	lfd   29,	40(%r10)
	fmul   2 , 25 , 12
	lfd   19,	56(%r11)
	fmul   3 , 25 , 13
	lfd   24,	48(%r10)
	stfd  8,	64(%r13)
	fmadd  0 , 25 , 14 , 0
	stfd  9,	72(%r13)
	fmadd  1 , 25 , 15 , 1
	stfd  10,	80(%r13)
	fnmsub 2 , 26 , 13 , 2
	stfd  11,	88(%r13)
	fmadd  3 , 26 , 12 , 3
	lfd   25,	56(%r10)
	fmul   4 , 28 , 12
	lfd   20,	64(%r11)
	fmadd  0 , 26 , 15 , 0
	lfd   21,	72(%r11)
	fnmsub 1 , 26 , 14 , 1
	lfd   22,	80(%r11)
	fmadd  2 , 27 , 14 , 2
	lfd   23,	88(%r11)
	fmadd  3 , 27 , 15 , 3
	lfd   26,	64(%r10)
	fmul   5 , 28 , 13
	addi  %r13 , %r12 , 0
	fnmsub 4 , 29 , 13 , 4
	lfd   27,	72(%r10)
	fmadd  0 , 28 , 16 , 0
	dcbt 	%r14,%r11
	fmadd  1 , 28 , 17 , 1
	fmadd  2 , 24 , 16 , 2
	la   %r12,	96(%r12)
	fmadd  5 , 29 , 12 , 5
	lfd   28,	80(%r10)
	fmadd  4 , 24 , 14 , 4
	la   %r9,	-1(%r9)
	fmadd  0 , 29 , 17 , 0
	fnmsub 1 , 29 , 16 , 1
	fmadd  3 , 24 , 17 , 3
	fmadd  5 , 24 , 15 , 5
	lfd   29,	88(%r10)
	fnmsub 4 , 25 , 15 , 4
	fmadd  2 , 25 , 17 , 2
	lfd   24,	96(%r10)
	fmul   6 , 27 , 12
	dcbt 	%r15,%r11
	fnmsub 3 , 25 , 16 , 3
	fmadd  5 , 25 , 14 , 5
	fmadd  4 , 26 , 16 , 4
	fmul   7 , 27 , 13
	lfd   25,	104(%r10)
	fnmsub 6 , 28 , 13 , 6
	fmadd  0 , 27 , 18 , 0
	fmadd  5 , 26 , 17 , 5
	dcbt 	%r16,%r11
	fmadd  1 , 27 , 19 , 1
	fmadd  7 , 28 , 12 , 7
	lfd   26,	112(%r10)
	fmadd  6 , 29 , 14 , 6
	lfd   27,	120(%r10)
	fmadd  0 , 28 , 19 , 0
	la   %r11,	96(%r11)
	fmadd  2 , 29 , 18 , 2
	fnmsub 1 , 28 , 18 , 1
	fmadd  7 , 29 , 15 , 7
	fnmsub 6 , 24 , 15 , 6
	lfd   28,	192(%r10)
	fmadd  3 , 29 , 19 , 3
	dcbt 	%r17,%r10
	fmadd  2 , 24 , 19 , 2
	fmadd  4 , 25 , 18 , 4
	lfd   29,	280(%r10)
	fmadd  7 , 24 , 14 , 7
	fmadd  6 , 25 , 16 , 6
	fnmsub 3 , 24 , 18 , 3
	fmadd  5 , 25 , 19 , 5
	fmadd  4 , 26 , 19 , 4
	lfd   24,	128(%r10)
	fmadd  7 , 25 , 17 , 7
	dcbt 	%r18,%r10
	fnmsub 6 , 26 , 17 , 6
	fmul   8 , 28 , 20
	lfd   25,	136(%r10)
	fnmsub 5 , 26 , 18 , 5
	fmul   9 , 28 , 21
	fmadd  7 , 26 , 16 , 7
	fmadd  6 , 27 , 18 , 6
	lfd   28,	144(%r10)
	fmul   10 , 29 , 22
	lfd   26,	200(%r10)
	fmul   11 , 29 , 23
	dcbt 	%r19,%r10
	fmadd  8 , 24 , 12 , 8
	fmadd  7 , 27 , 19 , 7
	lfd   29,	152(%r10)
	fmadd  9 , 24 , 13 , 9
	fmadd  0 , 24 , 20 , 0
	lfd   27,	208(%r10)
	fmadd  1 , 24 , 21 , 1
	fnmsub 8 , 25 , 13 , 8
	fmadd  10 , 26 , 12 , 10
	lfd   24,	216(%r10)
	fmadd  9 , 25 , 12 , 9
	dcbt 	%r20,%r10
	fmadd  0 , 25 , 21 , 0
	fnmsub 1 , 25 , 20 , 1
	fmadd  11 , 26 , 13 , 11
	fnmsub 10 , 27 , 13 , 10
	lfd   25,	224(%r10)
	fmadd  8 , 28 , 14 , 8
	fmadd  0 , 26 , 22 , 0
	lfd   13,	8(%r11)
	fmadd  1 , 26 , 23 , 1
	dcbt 	%r21,%r10
	fmadd  11 , 27 , 12 , 11
	fmadd  9 , 28 , 15 , 9
	lfd   26,	160(%r10)
	fnmsub 8 , 29 , 15 , 8
	lfd   12,	0(%r11)
	fmadd  0 , 27 , 23 , 0
	fnmsub 1 , 27 , 22 , 1
	fmadd  2 , 28 , 20 , 2
	fmadd  9 , 29 , 14 , 9
	lfd   27,	168(%r10)
	fmadd  3 , 28 , 21 , 3
	fmadd  10 , 24 , 14 , 10
	stfd  0,	0(%r13)
	fmadd  11 , 24 , 15 , 11
	lfd   28,	232(%r10)
	fmadd  2 , 29 , 21 , 2
	stfd  1,	8(%r13)
	fmadd  8 , 26 , 16 , 8
	dcbt 	%r22,%r10
	fnmsub 3 , 29 , 20 , 3
	fnmsub 10 , 25 , 15 , 10
	fmadd  11 , 25 , 14 , 11
	lfd   29,	240(%r10)
	fmadd  2 , 24 , 22 , 2
	lfd   15,	24(%r11)
	fmadd  9 , 26 , 17 , 9
	lfd   14,	16(%r11)
	fmadd  3 , 24 , 23 , 3
	fnmsub 8 , 27 , 17 , 8
	fmadd  4 , 26 , 20 , 4
	lfd   24,	176(%r10)
	fmadd  2 , 25 , 23 , 2
	dcbt 	%r23,%r10
	fmadd  9 , 27 , 16 , 9
	fnmsub 3 , 25 , 22 , 3
	fmadd  5 , 26 , 21 , 5
	fmadd  4 , 27 , 21 , 4
	lfd   25,	184(%r10)
	fmadd  10 , 28 , 16 , 10
	lfd   26,	248(%r10)
	fmadd  11 , 28 , 17 , 11
	stfd  2,	16(%r13)
	fmadd  8 , 24 , 18 , 8
	stfd  3,	24(%r13)
	fnmsub 5 , 27 , 20 , 5
	dcbt 	%r24,%r10
	fmadd  4 , 28 , 22 , 4
	fnmsub 10 , 29 , 17 , 10
	lfd   27,	256(%r10)
	fmadd  11 , 29 , 16 , 11
	fmadd  9 , 24 , 19 , 9
	lfd   17,	40(%r11)
	fmadd  5 , 28 , 23 , 5
	lfd   16,	32(%r11)
	fmadd  4 , 29 , 23 , 4
	fnmsub 8 , 25 , 19 , 8
	lfd   28,	264(%r10)
	fmadd  6 , 24 , 20 , 6
	dcbt 	%r25,%r10
	fmadd  9 , 25 , 18 , 9
	fnmsub 5 , 29 , 22 , 5
	fmadd  7 , 24 , 21 , 7
	stfd  4,	32(%r13)
	fmadd  10 , 26 , 18 , 10
	lfd   29,	272(%r10)
	fmadd  6 , 25 , 21 , 6
	fmadd  11 , 26 , 19 , 11
	la   %r10,	288(%r10)
	fmadd  8 , 28 , 22 , 8
	stfd  5,	40(%r13)
	fnmsub 7 , 25 , 20 , 7
	lfd   24,	0(%r10)
	fnmsub 10 , 27 , 19 , 10
	fmadd  6 , 26 , 22 , 6
	lfd   25,	8(%r10)
	fmadd  11 , 27 , 18 , 11
	fmadd  9 , 28 , 23 , 9
	fmadd  7 , 26 , 23 , 7
	lfd   18,	48(%r11)
	fmadd  10 , 28 , 20 , 10
	fmadd  6 , 27 , 23 , 6
	lfd   26,	16(%r10)
	fmadd  11 , 28 , 21 , 11
	fmadd  8 , 29 , 23 , 8
	fnmsub 7 , 27 , 22 , 7
	fnmsub 10 , 29 , 21 , 10
	fnmsub 9 , 29 , 22 , 9
	stfd  6,	48(%r13)
	fmadd  11 , 29 , 20 , 11
	lfd   27,	24(%r10)
	stfd  7,	56(%r13)
	bdnz    clover_mat_mlt_asm_lab1
	stfd  8,	64(%r13)
	stfd  9,	72(%r13)
	stfd  10,	80(%r13)
	stfd  11,	88(%r13)
clover_mat_mlt_asm_lab0:
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
