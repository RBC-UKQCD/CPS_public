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
.globl dec_hsu3
dec_hsu3:

	or    %r9 , %r3 , %r3
	la   %r1,	-1280(%r1)
	or    %r10 , %r4 , %r4
	addi  %r13 , %r1 , 736
	stw   %r14,	368(%r1)
	stw   %r15,	376(%r1)
	or    %r14 , %r5 , %r5
	stw   %r16,	384(%r1)
	stw   %r17,	392(%r1)
	stw   %r18,	400(%r1)
	stw   %r19,	408(%r1)
	stw   %r20,	416(%r1)
	stw   %r21,	424(%r1)
	stw   %r22,	432(%r1)
	or    %r21 , %r6 , %r6
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
	li   %r22,0
	li   %r23,32
	li   %r24,64
	li   %r25,96
	li   %r26,128
	lwz   %r14,	0(%r14)
	or.  %r14 , %r14 , %r14
	bf gt,  dec_hsu3_lab0
	addi  %r13 , %r1 , 736
	lwz   %r12,	0(%r21)
	lfd   0,	0(%r9)
	lfd   1,	8(%r9)
	lfd   6,	48(%r9)
	lfd   7,	56(%r9)
	lfd   12,	96(%r9)
	lfd   13,	104(%r9)
	lfd   18,	144(%r9)
	lfd   19,	152(%r9)
	lfd   2,	16(%r9)
	lfd   3,	24(%r9)
	lfd   8,	64(%r9)
	lfd   9,	72(%r9)
	lfd   14,	112(%r9)
	lfd   15,	120(%r9)
	lfd   20,	160(%r9)
	lfd   21,	168(%r9)
	mtctr 	%r14
	b       dec_hsu3_lab1
.align 4
dec_hsu3_lab1:
	addi  %r15 , %r13 , 160
	lfd   4,	32(%r9)
	addi  %r18 , %r12 , 0
	la   %r21,	4(%r21)
	addi  %r17 , %r13 , 256
	lfd   5,	40(%r9)
	fadd   24 , 0 , 19
	lwz   %r12,	0(%r21)
	fsub   25 , 1 , 18
	lfd   10,	80(%r9)
	addi  %r16 , %r12 , 0
	la   %r21,	4(%r21)
	fsub   26 , 0 , 19
	lfd   11,	88(%r9)
	fadd   27 , 1 , 18
	lwz   %r12,	0(%r21)
	fadd   28 , 6 , 13
	lfd   16,	128(%r9)
	fsub   29 , 7 , 12
	stfd  24,	0(%r18)
	fsub   30 , 6 , 13
	lfd   17,	136(%r9)
	fadd   31 , 7 , 12
	stfd  25,	8(%r18)
	fadd   24 , 2 , 21
	lfd   22,	176(%r9)
	fsub   25 , 3 , 20
	stfd  26,	0(%r15)
	addi  %r11 , %r13 , 160
	lfd   23,	184(%r9)
	fsub   26 , 2 , 21
	stfd  27,	8(%r15)
	stfd  28,	16(%r18)
	fadd   27 , 3 , 20
	stfd  29,	24(%r18)
	fadd   28 , 8 , 15
	stfd  30,	16(%r15)
	fsub   29 , 9 , 14
	stfd  31,	24(%r15)
	fsub   30 , 8 , 15
	stfd  24,	32(%r18)
	fadd   31 , 9 , 14
	stfd  25,	40(%r18)
	fadd   24 , 4 , 23
	stfd  26,	32(%r15)
	fsub   25 , 5 , 22
	stfd  27,	40(%r15)
	fsub   26 , 4 , 23
	stfd  28,	48(%r18)
	fadd   27 , 5 , 22
	stfd  29,	56(%r18)
	fadd   28 , 10 , 17
	stfd  30,	48(%r15)
	fsub   29 , 11 , 16
	stfd  31,	56(%r15)
	fsub   30 , 10 , 17
	stfd  24,	64(%r18)
	fadd   31 , 11 , 16
	stfd  25,	72(%r18)
	fadd   24 , 0 , 18
	stfd  26,	64(%r15)
	fadd   25 , 1 , 19
	stfd  27,	72(%r15)
	fsub   26 , 0 , 18
	stfd  28,	80(%r18)
	fsub   27 , 1 , 19
	stfd  29,	88(%r18)
	fsub   28 , 6 , 12
	stfd  30,	80(%r15)
	fsub   29 , 7 , 13
	stfd  31,	88(%r15)
	fadd   30 , 6 , 12
	stfd  24,	0(%r16)
	fadd   31 , 7 , 13
	stfd  25,	8(%r16)
	fadd   24 , 2 , 20
	stfd  26,	0(%r17)
	fadd   25 , 3 , 21
	stfd  27,	8(%r17)
	fsub   26 , 2 , 20
	stfd  28,	16(%r16)
	fsub   27 , 3 , 21
	stfd  29,	24(%r16)
	fsub   28 , 8 , 14
	stfd  30,	16(%r17)
	fsub   29 , 9 , 15
	stfd  31,	24(%r17)
	fadd   30 , 8 , 14
	stfd  24,	32(%r16)
	fadd   31 , 9 , 15
	stfd  25,	40(%r16)
	fadd   24 , 4 , 22
	stfd  26,	32(%r17)
	fadd   25 , 5 , 23
	stfd  27,	40(%r17)
	fsub   26 , 4 , 22
	stfd  28,	48(%r16)
	fsub   27 , 5 , 23
	stfd  29,	56(%r16)
	fsub   28 , 10 , 16
	stfd  30,	48(%r17)
	fsub   29 , 11 , 17
	stfd  31,	56(%r17)
	fadd   30 , 10 , 16
	stfd  24,	64(%r16)
	fadd   31 , 11 , 17
	stfd  25,	72(%r16)
	addi  %r15 , %r13 , 352
	stfd  26,	64(%r17)
	addi  %r18 , %r12 , 0
	stfd  27,	72(%r17)
	fadd   24 , 0 , 13
	stfd  28,	80(%r16)
	fsub   25 , 1 , 12
	stfd  29,	88(%r16)
	fsub   26 , 0 , 13
	stfd  30,	80(%r17)
	fadd   27 , 1 , 12
	stfd  31,	88(%r17)
	fsub   28 , 6 , 19
	dcbt 	%r22,%r10
	fadd   29 , 7 , 18
	la   %r21,	4(%r21)
	fadd   30 , 6 , 19
	stfd  24,	0(%r18)
	fsub   31 , 7 , 18
	lwz   %r12,	0(%r21)
	fadd   24 , 2 , 15
	stfd  25,	8(%r18)
	dcbt 	%r23,%r10
	fsub   25 , 3 , 14
	stfd  26,	0(%r15)
	stfd  27,	8(%r15)
	fsub   26 , 2 , 15
	stfd  28,	16(%r18)
	fadd   27 , 3 , 14
	stfd  29,	24(%r18)
	fsub   28 , 8 , 21
	dcbt 	%r24,%r10
	fadd   29 , 9 , 20
	stfd  30,	16(%r15)
	stfd  31,	24(%r15)
	fadd   30 , 8 , 21
	stfd  24,	32(%r18)
	fsub   31 , 9 , 20
	stfd  25,	40(%r18)
	fadd   24 , 4 , 17
	dcbt 	%r25,%r10
	fsub   25 , 5 , 16
	stfd  26,	32(%r15)
	stfd  27,	40(%r15)
	fsub   26 , 4 , 17
	stfd  28,	48(%r18)
	fadd   27 , 5 , 16
	stfd  29,	56(%r18)
	fsub   28 , 10 , 23
	dcbt 	%r26,%r10
	fadd   29 , 11 , 22
	stfd  30,	48(%r15)
	stfd  31,	56(%r15)
	fadd   30 , 10 , 23
	stfd  24,	64(%r18)
	fsub   31 , 11 , 22
	stfd  25,	72(%r18)
	fsub   24 , 0 , 12
	stfd  26,	64(%r15)
	fsub   25 , 1 , 13
	stfd  27,	72(%r15)
	fadd   26 , 0 , 12
	stfd  28,	80(%r18)
	fadd   27 , 1 , 13
	stfd  29,	88(%r18)
	fsub   28 , 6 , 18
	stfd  30,	80(%r15)
	addi  %r18 , %r12 , 0
	stfd  31,	88(%r15)
	fsub   29 , 7 , 19
	la   %r21,	4(%r21)
	addi  %r15 , %r13 , 448
	stfd  24,	0(%r18)
	fadd   30 , 6 , 18
	lwz   %r12,	0(%r21)
	fadd   31 , 7 , 19
	stfd  25,	8(%r18)
	fsub   24 , 2 , 14
	stfd  26,	0(%r15)
	fsub   25 , 3 , 15
	stfd  27,	8(%r15)
	fadd   26 , 2 , 14
	stfd  28,	16(%r18)
	fadd   27 , 3 , 15
	stfd  29,	24(%r18)
	addi  %r16 , %r18 , 0
	stfd  30,	16(%r15)
	stfd  31,	24(%r15)
	fadd   30 , 8 , 20
	li   %r20,3
	fadd   31 , 9 , 21
	stfd  24,	32(%r18)
	lfd   0,	0(%r11)
	fadd   24 , 4 , 16
	stfd  25,	40(%r18)
	lfd   2,	16(%r11)
	fadd   25 , 5 , 17
	stfd  26,	32(%r15)
	lfd   1,	8(%r11)
	fadd   26 , 10 , 22
	stfd  27,	40(%r15)
	lfd   3,	24(%r11)
	fadd   27 , 11 , 23
	stfd  30,	48(%r15)
	lfd   12,	0(%r10)
	fsub   30 , 8 , 20
	stfd  31,	56(%r15)
	lfd   18,	48(%r10)
	fsub   31 , 9 , 21
	stfd  24,	64(%r15)
	lfd   13,	8(%r10)
	fsub   24 , 4 , 16
	stfd  25,	72(%r15)
	lfd   19,	56(%r10)
	fsub   25 , 5 , 17
	stfd  26,	80(%r15)
	lfd   4,	32(%r11)
	fsub   26 , 10 , 22
	stfd  27,	88(%r15)
	lfd   6,	48(%r11)
	fsub   27 , 11 , 23
	lfd   5,	40(%r11)
	lfd   7,	56(%r11)
	mtctr 	%r20
	b       dec_hsu3_lab2
.align 4
dec_hsu3_lab2:
	fmul   28 , 18 , 0
	stfd  30,	48(%r16)
	fmul   29 , 18 , 1
	lfd   14,	16(%r10)
	addi  %r19 , %r10 , 144
	stfd  31,	56(%r16)
	fmul   30 , 18 , 2
	lfd   20,	64(%r10)
	fmul   31 , 18 , 3
	stfd  24,	64(%r16)
	fmadd  28 , 19 , 1 , 28
	lfd   15,	24(%r10)
	fmul   24 , 12 , 0
	stfd  25,	72(%r16)
	fnmsub 29 , 19 , 0 , 29
	lfd   21,	72(%r10)
	fmul   25 , 12 , 1
	stfd  26,	80(%r16)
	fmadd  30 , 19 , 3 , 30
	lfd   8,	64(%r11)
	fmul   26 , 12 , 2
	stfd  27,	88(%r16)
	fmadd  24 , 13 , 1 , 24
	lfd   10,	80(%r11)
	fmul   27 , 12 , 3
	addi  %r16 , %r12 , 0
	fnmsub 25 , 13 , 0 , 25
	lfd   9,	72(%r11)
	fmadd  28 , 20 , 4 , 28
	la   %r21,	4(%r21)
	fmadd  26 , 13 , 3 , 26
	lfd   11,	88(%r11)
	fmadd  24 , 14 , 4 , 24
	la   %r20,	-1(%r20)
	fnmsub 27 , 13 , 2 , 27
	lfd   16,	32(%r10)
	fmadd  25 , 14 , 5 , 25
	addi  %r11 , %r11 , 96
	fmadd  29 , 20 , 5 , 29
	lfd   22,	80(%r10)
	fmadd  26 , 14 , 6 , 26
	fmadd  24 , 15 , 5 , 24
	lfd   17,	40(%r10)
	fmadd  27 , 14 , 7 , 27
	fnmsub 25 , 15 , 4 , 25
	lfd   23,	88(%r10)
	fmadd  28 , 21 , 5 , 28
	fmadd  26 , 15 , 7 , 26
	lfd   12,	96(%r10)
	fnmsub 29 , 21 , 4 , 29
	fnmsub 27 , 15 , 6 , 27
	lfd   13,	104(%r10)
	fmadd  24 , 16 , 8 , 24
	fmadd  25 , 16 , 9 , 25
	lfd   14,	112(%r10)
	fmadd  26 , 16 , 10 , 26
	fmadd  28 , 22 , 8 , 28
	lfd   15,	120(%r10)
	fmadd  27 , 16 , 11 , 27
	dcbt 	%r22,%r19
	fmadd  29 , 22 , 9 , 29
	fmadd  24 , 17 , 9 , 24
	lfd   16,	128(%r10)
	fnmsub 25 , 17 , 8 , 25
	fmadd  26 , 17 , 11 , 26
	fnmsub 27 , 17 , 10 , 27
	dcbt 	%r23,%r19
	fmadd  28 , 23 , 9 , 28
	fnmsub 29 , 23 , 8 , 29
	stfd  24,	0(%r12)
	fnmsub 31 , 19 , 2 , 31
	lfd   17,	136(%r10)
	fmul   24 , 12 , 0
	stfd  25,	8(%r12)
	addi  %r10 , %r10 , 144
	fmul   25 , 12 , 1
	stfd  26,	16(%r12)
	fmadd  30 , 20 , 6 , 30
	lfd   18,	48(%r10)
	fmul   26 , 12 , 2
	stfd  27,	24(%r12)
	fmadd  24 , 13 , 1 , 24
	lfd   19,	56(%r10)
	fmul   27 , 12 , 3
	stfd  28,	32(%r12)
	fnmsub 25 , 13 , 0 , 25
	lfd   1,	8(%r11)
	fmadd  31 , 20 , 7 , 31
	stfd  29,	40(%r12)
	fmadd  26 , 13 , 3 , 26
	dcbt 	%r24,%r19
	fmadd  24 , 14 , 4 , 24
	addi  %r19 , %r19 , 96
	lwz   %r12,	0(%r21)
	fnmsub 27 , 13 , 2 , 27
	lfd   0,	0(%r11)
	fmadd  25 , 14 , 5 , 25
	fmadd  26 , 14 , 6 , 26
	dcbt 	%r22,%r19
	fmadd  30 , 21 , 7 , 30
	fnmsub 31 , 21 , 6 , 31
	lfd   2,	16(%r11)
	fmadd  27 , 14 , 7 , 27
	fmadd  24 , 15 , 5 , 24
	lfd   3,	24(%r11)
	fnmsub 25 , 15 , 4 , 25
	dcbt 	%r23,%r19
	fmadd  26 , 15 , 7 , 26
	fmadd  30 , 22 , 10 , 30
	lfd   12,	0(%r10)
	fnmsub 27 , 15 , 6 , 27
	fmadd  31 , 22 , 11 , 31
	lfd   13,	8(%r10)
	fmadd  24 , 16 , 8 , 24
	fmadd  25 , 16 , 9 , 25
	lfd   4,	32(%r11)
	fmadd  26 , 16 , 10 , 26
	fmadd  27 , 16 , 11 , 27
	lfd   6,	48(%r11)
	fmadd  30 , 23 , 11 , 30
	fnmsub 31 , 23 , 10 , 31
	lfd   5,	40(%r11)
	fmadd  24 , 17 , 9 , 24
	fnmsub 25 , 17 , 8 , 25
	lfd   7,	56(%r11)
	fmadd  26 , 17 , 11 , 26
	fnmsub 27 , 17 , 10 , 27
	bdnz    dec_hsu3_lab2
dec_hsu3_lab3:
	addi  %r9 , %r9 , 192
	stfd  30,	48(%r16)
	addi  %r19 , %r9 , 96
	lfd   14,	16(%r10)
	fmul   28 , 18 , 0
	stfd  31,	56(%r16)
	fmul   29 , 18 , 1
	lfd   20,	64(%r10)
	fmul   30 , 18 , 2
	stfd  24,	64(%r16)
	fmul   31 , 18 , 3
	lfd   15,	24(%r10)
	fmul   24 , 12 , 0
	stfd  25,	72(%r16)
	fmadd  28 , 19 , 1 , 28
	lfd   21,	72(%r10)
	fmul   25 , 12 , 1
	stfd  26,	80(%r16)
	fnmsub 29 , 19 , 0 , 29
	dcbt 	%r22,%r9
	fmul   26 , 12 , 2
	stfd  27,	88(%r16)
	fmadd  24 , 13 , 1 , 24
	lfd   8,	64(%r11)
	fmul   27 , 12 , 3
	addi  %r16 , %r12 , 0
	fnmsub 25 , 13 , 0 , 25
	lfd   10,	80(%r11)
	fmadd  28 , 20 , 4 , 28
	dcbt 	%r23,%r9
	fmadd  26 , 13 , 3 , 26
	lfd   9,	72(%r11)
	fmadd  24 , 14 , 4 , 24
	la   %r21,	4(%r21)
	fnmsub 27 , 13 , 2 , 27
	lfd   11,	88(%r11)
	fmadd  25 , 14 , 5 , 25
	fmadd  29 , 20 , 5 , 29
	dcbt 	%r24,%r9
	fmadd  26 , 14 , 6 , 26
	lfd   16,	32(%r10)
	fmadd  24 , 15 , 5 , 24
	addi  %r11 , %r11 , 96
	fmadd  27 , 14 , 7 , 27
	lfd   22,	80(%r10)
	fnmsub 25 , 15 , 4 , 25
	fmadd  28 , 21 , 5 , 28
	dcbt 	%r22,%r19
	fmadd  26 , 15 , 7 , 26
	lfd   17,	40(%r10)
	fnmsub 29 , 21 , 4 , 29
	fnmsub 27 , 15 , 6 , 27
	lfd   23,	88(%r10)
	fmadd  24 , 16 , 8 , 24
	fmadd  25 , 16 , 9 , 25
	dcbt 	%r23,%r19
	fmadd  26 , 16 , 10 , 26
	lfd   12,	96(%r10)
	fmadd  28 , 22 , 8 , 28
	fmadd  27 , 16 , 11 , 27
	lfd   13,	104(%r10)
	fmadd  29 , 22 , 9 , 29
	fmadd  24 , 17 , 9 , 24
	dcbt 	%r24,%r19
	fnmsub 25 , 17 , 8 , 25
	lfd   14,	112(%r10)
	fmadd  26 , 17 , 11 , 26
	fnmsub 27 , 17 , 10 , 27
	lfd   15,	120(%r10)
	fmadd  28 , 23 , 9 , 28
	fnmsub 29 , 23 , 8 , 29
	stfd  24,	0(%r12)
	fmadd  30 , 19 , 3 , 30
	lfd   16,	128(%r10)
	fmul   24 , 12 , 0
	stfd  25,	8(%r12)
	fnmsub 31 , 19 , 2 , 31
	lfd   17,	136(%r10)
	fmul   25 , 12 , 1
	stfd  26,	16(%r12)
	addi  %r10 , %r10 , 144
	fmul   26 , 12 , 2
	stfd  27,	24(%r12)
	fmadd  24 , 13 , 1 , 24
	fmul   27 , 12 , 3
	stfd  28,	32(%r12)
	fnmsub 25 , 13 , 0 , 25
	fmadd  30 , 20 , 6 , 30
	stfd  29,	40(%r12)
	fmadd  26 , 13 , 3 , 26
	fmadd  31 , 20 , 7 , 31
	lwz   %r12,	0(%r21)
	fnmsub 27 , 13 , 2 , 27
	fmadd  24 , 14 , 4 , 24
	fmadd  25 , 14 , 5 , 25
	fmadd  26 , 14 , 6 , 26
	fmadd  30 , 21 , 7 , 30
	fmadd  27 , 14 , 7 , 27
	fnmsub 31 , 21 , 6 , 31
	fmadd  24 , 15 , 5 , 24
	fnmsub 25 , 15 , 4 , 25
	fmadd  26 , 15 , 7 , 26
	fnmsub 27 , 15 , 6 , 27
	fmadd  30 , 22 , 10 , 30
	fmadd  31 , 22 , 11 , 31
	fmadd  24 , 16 , 8 , 24
	fmadd  25 , 16 , 9 , 25
	fmadd  26 , 16 , 10 , 26
	fmadd  27 , 16 , 11 , 27
	fmadd  30 , 23 , 11 , 30
	fnmsub 31 , 23 , 10 , 31
	fmadd  24 , 17 , 9 , 24
	fnmsub 25 , 17 , 8 , 25
	fmadd  26 , 17 , 11 , 26
	fnmsub 27 , 17 , 10 , 27
dec_hsu3_lab4:
	stfd  30,	48(%r16)
	lfd   0,	0(%r9)
	stfd  31,	56(%r16)
	lfd   1,	8(%r9)
	stfd  24,	64(%r16)
	lfd   6,	48(%r9)
	stfd  25,	72(%r16)
	lfd   7,	56(%r9)
	stfd  26,	80(%r16)
	lfd   12,	96(%r9)
	stfd  27,	88(%r16)
	lfd   13,	104(%r9)
	lfd   18,	144(%r9)
	lfd   19,	152(%r9)
	lfd   2,	16(%r9)
	lfd   3,	24(%r9)
	lfd   8,	64(%r9)
	lfd   9,	72(%r9)
	lfd   14,	112(%r9)
	lfd   15,	120(%r9)
	lfd   20,	160(%r9)
	la   %r14,	-1(%r14)
	lfd   21,	168(%r9)
# Instruction order barrier inserted
	or.  %r14 , %r14 , %r14
	bt gt,  dec_hsu3_lab1
dec_hsu3_lab0:
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
	la   %r1,	1280(%r1)
	blr
