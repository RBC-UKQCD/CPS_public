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
.globl rec_su3t
rec_su3t:

	or    %r9 , %r3 , %r3
	la   %r1,	-1408(%r1)
	or    %r10 , %r4 , %r4
	or    %r12 , %r5 , %r5
	stw   %r14,	368(%r1)
	stw   %r15,	376(%r1)
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
	addi  %r26 , %r1 , 864
	stw   %r28,	480(%r1)
	or    %r27 , %r6 , %r6
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
	lwz   %r27,	0(%r27)
	or.  %r27 , %r27 , %r27
	bf gt,  rec_su3_lab0
	li   %r28,0
	li   %r29,32
	li   %r30,64
	lfd   0,	0(%r10)
	lfd   6,	16(%r10)
	lfd   1,	8(%r10)
	lfd   7,	24(%r10)
	lfd   2,	48(%r10)
	lfd   20,	64(%r10)
	lfd   3,	56(%r10)
	lfd   21,	72(%r10)
	lfd   8,	0(%r12)
	lfd   9,	8(%r12)
	lfd   10,	16(%r12)
	lfd   11,	24(%r12)
	lfd   16,	32(%r12)
	lfd   17,	40(%r12)
	lfd   18,	48(%r12)
	lfd   19,	56(%r12)
	mtctr 	%r27
	mftbl %r8 
rec_su3_lab1:
	addi  %r11 , %r10 , 0
	dcbt 	%r28,%r12
	addi  %r25 , %r11 , 96
	addi  %r13 , %r26 , 160
	addi  %r24 , %r26 , 64
	dcbt 	%r29,%r12
	li   %r23,3
	dcbt 	%r28,%r11
	lfd   8,	0(%r12)
	dcbt 	%r29,%r11
	lfd   9,	8(%r12)
	dcbt 	%r30,%r11
	lfd   10,	16(%r12)
	dcbt 	%r30,%r12
	lfd   11,	24(%r12)
	dcbt 	%r28,%r25
	lfd   16,	32(%r12)
	dcbt 	%r29,%r25
	lfd   17,	40(%r12)
	lfd   18,	48(%r12)
	lfd   19,	56(%r12)
	mtctr 	%r23
rec_su3_lab2:
	addi  %r14 , %r12 , 128
	stfd  30,	48(%r24)
	addi  %r11 , %r10 , 144
	lfd   12,	64(%r12)
	addi  %r25 , %r11 , 96
	stfd  31,	56(%r24)
	fmul   28 , 6 , 8
	lfd   13,	72(%r12)
	fmul   29 , 6 , 9
	stfd  24,	64(%r24)
	fmul   30 , 6 , 10
	lfd   14,	80(%r12)
	fmul   24 , 0 , 8
	stfd  25,	72(%r24)
	fmul   31 , 6 , 11
	lfd   15,	88(%r12)
	fmul   25 , 0 , 9
	stfd  26,	80(%r24)
	fnmsub 28 , 7 , 9 , 28
	lfd   4,	96(%r10)
	fmul   26 , 0 , 10
	stfd  27,	88(%r24)
	fnmsub 24 , 1 , 9 , 24
	lfd   22,	112(%r10)
	fmul   27 , 0 , 11
	lfd   5,	104(%r10)
	fmadd  25 , 1 , 8 , 25
	dcbt 	%r28,%r14
	fmadd  29 , 7 , 8 , 29
	lfd   23,	120(%r10)
	fnmsub 26 , 1 , 11 , 26
	lfd   0,	32(%r10)
	fmadd  24 , 2 , 16 , 24
	addi  %r12 , %r12 , 128
	fmadd  27 , 1 , 10 , 27
	addi  %r24 , %r13 , 0
	fmadd  25 , 2 , 17 , 25
	la   %r23,	-1(%r23)
	fmadd  28 , 20 , 16 , 28
	dcbt 	%r29,%r14
	fmadd  26 , 2 , 18 , 26
	lfd   1,	40(%r10)
	fmadd  29 , 20 , 17 , 29
	fmadd  27 , 2 , 19 , 27
	fnmsub 24 , 3 , 17 , 24
	fmadd  25 , 3 , 16 , 25
	lfd   2,	80(%r10)
	fnmsub 26 , 3 , 19 , 26
	dcbt 	%r30,%r14
	fnmsub 28 , 21 , 17 , 28
	fmadd  27 , 3 , 18 , 27
	fmadd  29 , 21 , 16 , 29
	fmadd  24 , 4 , 12 , 24
	lfd   3,	88(%r10)
	fmadd  25 , 4 , 13 , 25
	fmadd  26 , 4 , 14 , 26
	dcbt 	%r28,%r11
	fmadd  27 , 4 , 15 , 27
	fmadd  28 , 22 , 12 , 28
	fmadd  29 , 22 , 13 , 29
	lfd   4,	128(%r10)
	fnmsub 24 , 5 , 13 , 24
	fmadd  25 , 5 , 12 , 25
	fnmsub 26 , 5 , 15 , 26
	dcbt 	%r29,%r11
	fmadd  27 , 5 , 14 , 27
	fnmsub 28 , 23 , 13 , 28
	fmadd  29 , 23 , 12 , 29
	stfd  24,	0(%r13)
	fnmsub 30 , 7 , 11 , 30
	lfd   5,	136(%r10)
	fmul   24 , 0 , 8
	stfd  25,	8(%r13)
	addi  %r10 , %r10 , 144
	dcbt 	%r30,%r11
	fmul   25 , 0 , 9
	stfd  26,	16(%r13)
	fmadd  31 , 7 , 10 , 31
	lfd   6,	16(%r10)
	fmul   26 , 0 , 10
	stfd  27,	24(%r13)
	fnmsub 24 , 1 , 9 , 24
	lfd   7,	24(%r10)
	fmul   27 , 0 , 11
	stfd  28,	32(%r13)
	fmadd  25 , 1 , 8 , 25
	dcbt 	%r28,%r25
	fmadd  30 , 20 , 18 , 30
	stfd  29,	40(%r13)
	fnmsub 26 , 1 , 11 , 26
	lfd   0,	0(%r10)
	fmadd  31 , 20 , 19 , 31
	addi  %r13 , %r13 , 96
	fmadd  27 , 1 , 10 , 27
	lfd   8,	0(%r12)
	fmadd  24 , 2 , 16 , 24
	lfd   20,	64(%r10)
	fmadd  25 , 2 , 17 , 25
	dcbt 	%r29,%r25
	fmadd  26 , 2 , 18 , 26
	lfd   1,	8(%r10)
	fnmsub 30 , 21 , 19 , 30
	lfd   9,	8(%r12)
	fmadd  27 , 2 , 19 , 27
/*	mftbl %r8  357,397 cycles (93,147 dropped) (355,355)*/
	fmadd  31 , 21 , 18 , 31
	lfd   10,	16(%r12)
	fnmsub 24 , 3 , 17 , 24
	lfd   2,	48(%r10)
	fmadd  25 , 3 , 16 , 25
	lfd   21,	72(%r10)
	fnmsub 26 , 3 , 19 , 26
	fmadd  27 , 3 , 18 , 27
	lfd   11,	24(%r12)
	fmadd  30 , 22 , 14 , 30
	lfd   16,	32(%r12)
	fmadd  31 , 22 , 15 , 31
	lfd   3,	56(%r10)
	fmadd  24 , 4 , 12 , 24
/*	mftbl %r8 343,388 cycles (93,147 dropped)*/
	fmadd  25 , 4 , 13 , 25
	lfd   17,	40(%r12)
	fmadd  26 , 4 , 14 , 26
	lfd   18,	48(%r12)
	fmadd  27 , 4 , 15 , 27
	fnmsub 30 , 23 , 15 , 30
	lfd   19,	56(%r12)
	fmadd  31 , 23 , 14 , 31
	fnmsub 24 , 5 , 13 , 24
	fmadd  25 , 5 , 12 , 25
	fnmsub 26 , 5 , 15 , 26
	fmadd  27 , 5 , 14 , 27
	bdnz    rec_su3_lab2
	nop
rec_su3_lab3:
	nop
   /*	mftbl %r31 */
	addi  %r14 , %r12 , 128
	stfd  30,	48(%r24)
	fmul   28 , 6 , 8
	stfd  31,	56(%r24)
	fmul   29 , 6 , 9
	lfd   12,	64(%r12)
	fmul   30 , 6 , 10
	stfd  24,	64(%r24)
	fmul   31 , 6 , 11
	stfd  25,	72(%r24)
	fmul   24 , 0 , 8
	lfd   13,	72(%r12)
	fmul   25 , 0 , 9
	stfd  26,	80(%r24)
	fnmsub 28 , 7 , 9 , 28
	stfd  27,	88(%r24)
	fmul   26 , 0 , 10
	lfd   14,	80(%r12)
	fmul   27 , 0 , 11
	lfd   15,	88(%r12)
	fnmsub 24 , 1 , 9 , 24
	addi  %r24 , %r13 , 0
	fmadd  25 , 1 , 8 , 25
	lfd   4,	96(%r10)
	fmadd  29 , 7 , 8 , 29
	lfd   22,	112(%r10)
	fnmsub 26 , 1 , 11 , 26
	addi  %r12 , %r12 , 128
	fmadd  27 , 1 , 10 , 27
	lfd   5,	104(%r10)
	fmadd  24 , 2 , 16 , 24
	lfd   23,	120(%r10)
	fmadd  25 , 2 , 17 , 25
	dcbt 	%r28,%r14
	fmadd  28 , 20 , 16 , 28
	addi  %r19 , %r26 , 160
	fmadd  26 , 2 , 18 , 26
	lfd   0,	32(%r10)
	fmadd  27 , 2 , 19 , 27
	lfd   1,	40(%r10)
	fmadd  29 , 20 , 17 , 29
	addi  %r15 , %r12 , 0
	fnmsub 24 , 3 , 17 , 24
	dcbt 	%r29,%r14
	fmadd  25 , 3 , 16 , 25
	addi  %r15 , %r12 , 0
	fnmsub 26 , 3 , 19 , 26
	lfd   2,	80(%r10)
	fmadd  27 , 3 , 18 , 27
	addi  %r15 , %r12 , 0
	fnmsub 28 , 21 , 17 , 28
	addi  %r20 , %r26 , 256
	fmadd  29 , 21 , 16 , 29
	dcbt 	%r30,%r14
	fmadd  24 , 4 , 12 , 24
	addi  %r15 , %r12 , 0
	addi  %r14 , %r14 , 128
	lfd   3,	88(%r10)
	fmadd  25 , 4 , 13 , 25
	addi  %r16 , %r12 , 128
	fmadd  26 , 4 , 14 , 26
	addi  %r21 , %r26 , 352
	fmadd  27 , 4 , 15 , 27
	dcbt 	%r28,%r14
	fmadd  28 , 22 , 12 , 28
	addi  %r22 , %r26 , 448
	fmadd  29 , 22 , 13 , 29
	lfd   4,	128(%r10)
	fnmsub 24 , 5 , 13 , 24
	addi  %r17 , %r16 , 128
	fmadd  25 , 5 , 12 , 25
	addi  %r18 , %r17 , 128
	fnmsub 26 , 5 , 15 , 26
	dcbt 	%r29,%r14
	fmadd  27 , 5 , 14 , 27
	fnmsub 28 , 23 , 13 , 28
	la   %r27,	-1(%r27)
	fmadd  29 , 23 , 12 , 29
	stfd  24,	0(%r13)
	fnmsub 30 , 7 , 11 , 30
	stfd  25,	8(%r13)
	fmul   24 , 0 , 8
	dcbt 	%r30,%r14
	fmul   25 , 0 , 9
	addi  %r14 , %r14 , 128
	stfd  26,	16(%r13)
	fmadd  31 , 7 , 10 , 31
	stfd  27,	24(%r13)
	fmul   26 , 0 , 10
	lfd   5,	136(%r10)
	fmul   27 , 0 , 11
	dcbt 	%r28,%r14
	addi  %r10 , %r10 , 144
	fnmsub 24 , 1 , 9 , 24
	stfd  28,	32(%r13)
	fmadd  25 , 1 , 8 , 25
	stfd  29,	40(%r13)
	fnmsub 26 , 1 , 11 , 26
	lfd   9,	8(%r19)
	fmadd  27 , 1 , 10 , 27
	dcbt 	%r29,%r14
	fmadd  30 , 20 , 18 , 30
	addi  %r13 , %r13 , 96
	fmadd  31 , 20 , 19 , 31
	lfd   8,	0(%r19)
	fmadd  24 , 2 , 16 , 24
	lfd   10,	16(%r19)
	fmadd  25 , 2 , 17 , 25
	fmadd  26 , 2 , 18 , 26
	dcbt 	%r30,%r14
	fmadd  27 , 2 , 19 , 27
	addi  %r14 , %r14 , 128
	lfd   11,	24(%r19)
	fnmsub 30 , 21 , 19 , 30
	fmadd  31 , 21 , 18 , 31
	lfd   20,	0(%r16)
	fnmsub 24 , 3 , 17 , 24
	dcbt 	%r28,%r14
	fmadd  25 , 3 , 16 , 25
	fnmsub 26 , 3 , 19 , 26
	lfd   17,	8(%r12)
	fmadd  27 , 3 , 18 , 27
	lfd   16,	0(%r12)
	fmadd  30 , 22 , 14 , 30
	/*mftbl %r8 244,289 cycles (2,47 dropped)*/
	fmadd  31 , 22 , 15 , 31
	dcbt 	%r29,%r14
	fmadd  24 , 4 , 12 , 24
	fmadd  25 , 4 , 13 , 25
	lfd   18,	16(%r12)
	fmadd  26 , 4 , 14 , 26
	lfd   19,	24(%r12)
	fmadd  27 , 4 , 15 , 27
	fnmsub 30 , 23 , 15 , 30
	dcbt 	%r30,%r14
	fmadd  31 , 23 , 14 , 31
	/*mftbl %r8 236,279 cycles (2,47 dropped)*/
	addi  %r12 , %r18 , 128
	fnmsub 24 , 5 , 13 , 24
	lfd   21,	8(%r16)
	fmadd  25 , 5 , 12 , 25
	fnmsub 26 , 5 , 15 , 26
	lfd   13,	8(%r20)
	/*mftbl %r8 231,277 cycles (2,47 dropped)*/
	fmadd  27 , 5 , 14 , 27
	stfd  30,	48(%r24)
	/* mftbl %r8 230,275 cycles (2,47 dropped)*/
	fadd   0 , 16 , 8
	/*mftbl %r8 229 cycles*/
	stfd  31,	56(%r24)
	/* mftbl %r8  228 cycles*/
	fadd   1 , 17 , 9
	stfd  24,	64(%r24)
	fadd   2 , 18 , 10
	stfd  25,	72(%r24)
	fadd   3 , 19 , 11
	/* mftbl %r8 223 cycles*/
	stfd  26,	80(%r24)
	fsub   4 , 19 , 11
	stfd  27,	88(%r24)
	fsub   5 , 10 , 18
	lfd   12,	0(%r20)
	fsub   6 , 17 , 9
	fsub   7 , 8 , 16
	lfd   14,	16(%r20)
	fadd   0 , 0 , 20
/*	mftbl %r8 218,219 cycles (6 dropped)*/
	fadd   1 , 1 , 21
	lfd   22,	16(%r16)
	fsub   6 , 6 , 20
	lfd   15,	24(%r20)
	fsub   7 , 7 , 21
	fadd   0 , 0 , 12
	lfd   23,	24(%r16)
	fadd   1 , 1 , 13
	lfd   8,	0(%r21)
	fadd   2 , 2 , 22
	fadd   4 , 4 , 22
	lfd   16,	0(%r17)
	fadd   6 , 6 , 12
	fadd   7 , 7 , 13
	lfd   9,	8(%r21)
	fadd   3 , 3 , 23
	fadd   2 , 2 , 14
	lfd   17,	8(%r17)
/*	mftbl %r8 204 cycles */
	fadd   5 , 5 , 23
	fsub   4 , 4 , 14
	lfd   10,	16(%r21)
	fadd   0 , 0 , 16
	fadd   3 , 3 , 15
	lfd   18,	16(%r17)
	fsub   5 , 5 , 15
	lfd   11,	24(%r21)
	fadd   1 , 1 , 17
	fadd   0 , 0 , 8
	lfd   19,	24(%r17)
	fadd   4 , 4 , 17
	lfd   12,	0(%r22)
	fadd   2 , 2 , 18
	fadd   1 , 1 , 9
	lfd   20,	0(%r18)
	fsub   5 , 5 , 16
	fsub   4 , 4 , 9
	lfd   13,	8(%r22)
	fadd   3 , 3 , 19
	fadd   2 , 2 , 10
	lfd   21,	8(%r18)
	fsub   6 , 6 , 19
	fadd   5 , 5 , 8
	lfd   14,	16(%r22)
	fadd   7 , 7 , 18
	fadd   3 , 3 , 11
	lfd   22,	16(%r18)
	fadd   0 , 0 , 20
	fadd   6 , 6 , 11
	lfd   15,	24(%r22)
	fadd   1 , 1 , 21
	fsub   7 , 7 , 10
	lfd   23,	24(%r18)
	fadd   4 , 4 , 20
	fadd   5 , 5 , 21
	lfd   8,	32(%r19)
	fadd   2 , 2 , 22
	fadd   6 , 6 , 22
	lfd   16,	32(%r15)
	fadd   0 , 0 , 12
	fadd   1 , 1 , 13
	lfd   9,	40(%r19)
	fadd   3 , 3 , 23
	fadd   7 , 7 , 23
	lfd   17,	40(%r15)
	fadd   2 , 2 , 14
	fsub   4 , 4 , 12
	lfd   10,	48(%r19)
	fsub   5 , 5 , 13
	stfd  0,	0(%r9)
	fadd   3 , 3 , 15
	lfd   18,	48(%r15)
	fsub   6 , 6 , 14
	stfd  1,	8(%r9)
	fsub   7 , 7 , 15
	lfd   11,	56(%r19)
	fadd   0 , 16 , 8
	stfd  2,	48(%r9)
	fadd   1 , 17 , 9
	lfd   19,	56(%r15)
	stfd  3,	56(%r9)
	stfd  4,	96(%r9)
	fadd   2 , 18 , 10
	stfd  5,	104(%r9)
	stfd  6,	144(%r9)
	fsub   5 , 10 , 18
	stfd  7,	152(%r9)
	fsub   6 , 17 , 9
	lfd   12,	32(%r20)
	fadd   3 , 19 , 11
	fsub   4 , 19 , 11
	lfd   20,	32(%r16)
	fsub   7 , 8 , 16
	lfd   13,	40(%r20)
	lfd   21,	40(%r16)
	lfd   14,	48(%r20)
	fadd   0 , 0 , 20
	fsub   6 , 6 , 20
	lfd   22,	48(%r16)
	lfd   15,	56(%r20)
	fadd   1 , 1 , 21
	fadd   0 , 0 , 12
	lfd   23,	56(%r16)
	fsub   7 , 7 , 21
	fadd   6 , 6 , 12
	lfd   8,	32(%r21)
	fadd   2 , 2 , 22
	fadd   1 , 1 , 13
	lfd   16,	32(%r17)
	fadd   4 , 4 , 22
	fadd   7 , 7 , 13
	lfd   9,	40(%r21)
	fadd   3 , 3 , 23
	fadd   2 , 2 , 14
	lfd   17,	40(%r17)
	fadd   5 , 5 , 23
	fsub   4 , 4 , 14
	lfd   10,	48(%r21)
	fadd   0 , 0 , 16
	fadd   3 , 3 , 15
	lfd   18,	48(%r17)
	fsub   5 , 5 , 15
	lfd   11,	56(%r21)
	fadd   1 , 1 , 17
	fadd   0 , 0 , 8
	lfd   19,	56(%r17)
	fadd   4 , 4 , 17
	lfd   12,	32(%r22)
	fadd   2 , 2 , 18
	fadd   1 , 1 , 9
	lfd   20,	32(%r18)
	fsub   5 , 5 , 16
	fsub   4 , 4 , 9
	lfd   13,	40(%r22)
	fadd   3 , 3 , 19
	fadd   2 , 2 , 10
	lfd   21,	40(%r18)
	fsub   6 , 6 , 19
	fadd   5 , 5 , 8
	lfd   14,	48(%r22)
	fadd   7 , 7 , 18
	fadd   3 , 3 , 11
	lfd   22,	48(%r18)
	fadd   0 , 0 , 20
	fadd   6 , 6 , 11
	lfd   15,	56(%r22)
	fadd   1 , 1 , 21
	fsub   7 , 7 , 10
	lfd   23,	56(%r18)
	fadd   4 , 4 , 20
	fadd   5 , 5 , 21
	lfd   8,	64(%r19)
	fadd   2 , 2 , 22
	fadd   6 , 6 , 22
	lfd   16,	64(%r15)
	fadd   0 , 0 , 12
	fadd   1 , 1 , 13
	lfd   9,	72(%r19)
	fadd   3 , 3 , 23
	fadd   7 , 7 , 23
	lfd   17,	72(%r15)
	fadd   2 , 2 , 14
	fsub   4 , 4 , 12
	lfd   10,	80(%r19)
	fsub   5 , 5 , 13
	stfd  0,	16(%r9)
	fadd   3 , 3 , 15
	lfd   18,	80(%r15)
	fsub   6 , 6 , 14
	stfd  1,	24(%r9)
	fsub   7 , 7 , 15
	lfd   11,	88(%r19)
	fadd   0 , 16 , 8
	stfd  2,	64(%r9)
	fadd   1 , 17 , 9
	lfd   19,	88(%r15)
	stfd  3,	72(%r9)
	stfd  4,	112(%r9)
	fadd   2 , 18 , 10
	stfd  5,	120(%r9)
	stfd  6,	160(%r9)
	fsub   5 , 10 , 18
	stfd  7,	168(%r9)
	fsub   6 , 17 , 9
	lfd   12,	64(%r20)
	fadd   3 , 19 , 11
	fsub   4 , 19 , 11
	lfd   20,	64(%r16)
	fsub   7 , 8 , 16
	lfd   13,	72(%r20)
	lfd   21,	72(%r16)
	lfd   14,	80(%r20)
	fadd   0 , 0 , 20
	fsub   6 , 6 , 20
	lfd   22,	80(%r16)
	lfd   15,	88(%r20)
	fadd   1 , 1 , 21
	fadd   0 , 0 , 12
	lfd   23,	88(%r16)
	fsub   7 , 7 , 21
	fadd   6 , 6 , 12
	lfd   8,	64(%r21)
	fadd   2 , 2 , 22
	fadd   1 , 1 , 13
	lfd   16,	64(%r17)
	fadd   4 , 4 , 22
	fadd   7 , 7 , 13
	lfd   9,	72(%r21)
	fadd   3 , 3 , 23
	fadd   2 , 2 , 14
	lfd   17,	72(%r17)
	fadd   5 , 5 , 23
	fsub   4 , 4 , 14
	lfd   10,	80(%r21)
	fadd   0 , 0 , 16
	fadd   3 , 3 , 15
	lfd   18,	80(%r17)
	fsub   5 , 5 , 15
	lfd   11,	88(%r21)
	fadd   1 , 1 , 17
	fadd   0 , 0 , 8
	lfd   19,	88(%r17)
	fadd   4 , 4 , 17
	lfd   12,	64(%r22)
	fadd   2 , 2 , 18
	fadd   1 , 1 , 9
	lfd   20,	64(%r18)
	fsub   5 , 5 , 16
	fsub   4 , 4 , 9
	lfd   13,	72(%r22)
	fadd   3 , 3 , 19
	fadd   2 , 2 , 10
	lfd   21,	72(%r18)
	fsub   6 , 6 , 19
	fadd   5 , 5 , 8
	lfd   14,	80(%r22)
	fadd   7 , 7 , 18
	fadd   3 , 3 , 11
	lfd   22,	80(%r18)
	fadd   0 , 0 , 20
	fadd   6 , 6 , 11
	lfd   15,	88(%r22)
	fadd   1 , 1 , 21
	fsub   7 , 7 , 10
	lfd   23,	88(%r18)
	fadd   4 , 4 , 20
	fadd   5 , 5 , 21
	fadd   2 , 2 , 22
	lfd   20,	64(%r10)
	fadd   6 , 6 , 22
	fadd   0 , 0 , 12
	lfd   21,	72(%r10)
	fadd   1 , 1 , 13
	fadd   3 , 3 , 23
	fadd   7 , 7 , 23
	fadd   2 , 2 , 14
	fsub   4 , 4 , 12
	stfd  0,	32(%r9)
	fsub   5 , 5 , 13
	stfd  1,	40(%r9)
	fadd   3 , 3 , 15
	lfd   0,	0(%r10)
	fsub   6 , 6 , 14
	fsub   7 , 7 , 15
	stfd  2,	80(%r9)
	stfd  4,	128(%r9)
	stfd  5,	136(%r9)
	stfd  3,	88(%r9)
	stfd  6,	176(%r9)
	stfd  7,	184(%r9)
	lfd   6,	16(%r10)
	addi  %r9 , %r9 , 192
	lfd   1,	8(%r10)
	lfd   7,	24(%r10)
	lfd   2,	48(%r10)
	lfd   3,	56(%r10)
# Instruction order barrier inserted
	or.  %r27 , %r27 , %r27
	bt gt,  rec_su3_lab1
rec_su3_lab0:
	mftbl %r31
	stw %r8,0(%r7)
	stw %r31,4(%r7)
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
	la   %r1,	1408(%r1)
	blr
