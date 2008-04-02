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
.globl rec_su3_dag
rec_su3_dag:

	or    %r9 , %r3 , %r3
	la   %r1,	-1280(%r1)
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
	addi  %r26 , %r1 , 736
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
	bf gt,  rec_su3_dag_lab0
	addi  %r24 , %r26 , 64
	li   %r28,0
	li   %r29,32
	li   %r30,64
	li   %r31,96
	li   %r3,128
	lfd   4,	0(%r10)
	lfd   10,	16(%r10)
	lfd   5,	8(%r10)
	lfd   11,	24(%r10)
	lfd   18,	48(%r10)
	lfd   16,	64(%r10)
	lfd   19,	56(%r10)
	lfd   17,	72(%r10)
	dcbt 	%r28,%r12
	lfd   24,	0(%r12)
	lfd   25,	8(%r12)
	lfd   26,	16(%r12)
	lfd   27,	24(%r12)
	dcbt 	%r29,%r12
	lfd   28,	32(%r12)
	lfd   29,	40(%r12)
	lfd   30,	48(%r12)
	lfd   31,	56(%r12)
	dcbt 	%r30,%r12
	mtctr 	%r27
	b       rec_su3_dag_lab1
.align 4
rec_su3_dag_lab1:
	addi  %r13 , %r26 , 160
	li   %r23,3
	mtctr 	%r23
	b       rec_su3_dag_lab2
.align 4
rec_su3_dag_lab2:
	addi  %r14 , %r12 , 128
	stfd  6,	48(%r24)
	addi  %r11 , %r10 , 144
	lfd   0,	64(%r12)
	fmul   20 , 10 , 24
	stfd  7,	56(%r24)
	fmul   21 , 10 , 25
	lfd   1,	72(%r12)
	fmul   6 , 10 , 26
	stfd  14,	64(%r24)
	fmul   7 , 10 , 27
	lfd   2,	80(%r12)
	fmul   14 , 4 , 24
	stfd  15,	72(%r24)
	fnmsub 20 , 11 , 25 , 20
	lfd   3,	88(%r12)
	fmul   15 , 4 , 25
	stfd  22,	80(%r24)
	fmadd  21 , 11 , 24 , 21
	lfd   8,	96(%r10)
	fmul   22 , 4 , 26
	stfd  23,	88(%r24)
	fnmsub 14 , 5 , 25 , 14
	lfd   12,	112(%r10)
	fmul   23 , 4 , 27
	lfd   9,	104(%r10)
	fmadd  15 , 5 , 24 , 15
	addi  %r12 , %r12 , 128
	fmadd  20 , 16 , 28 , 20
	lfd   13,	120(%r10)
	fnmsub 22 , 5 , 27 , 22
	dcbt 	%r28,%r14
	fmadd  14 , 18 , 28 , 14
	addi  %r24 , %r13 , 0
	fmadd  23 , 5 , 26 , 23
	lfd   4,	32(%r10)
	fmadd  15 , 18 , 29 , 15
	la   %r23,	-1(%r23)
	fmadd  21 , 16 , 29 , 21
	lfd   5,	40(%r10)
	fmadd  22 , 18 , 30 , 22
	fnmsub 14 , 19 , 29 , 14
	dcbt 	%r29,%r14
	fmadd  23 , 18 , 31 , 23
	fmadd  15 , 19 , 28 , 15
	fnmsub 20 , 17 , 29 , 20
	lfd   18,	80(%r10)
	fnmsub 22 , 19 , 31 , 22
	fmadd  21 , 17 , 28 , 21
	fmadd  23 , 19 , 30 , 23
	dcbt 	%r30,%r14
	fmadd  14 , 8 , 0 , 14
	fmadd  15 , 8 , 1 , 15
	lfd   19,	88(%r10)
	fmadd  22 , 8 , 2 , 22
	fmadd  20 , 12 , 0 , 20
	fmadd  23 , 8 , 3 , 23
	fmadd  21 , 12 , 1 , 21
	dcbt 	%r28,%r11
	fnmsub 14 , 9 , 1 , 14
	fmadd  15 , 9 , 0 , 15
	lfd   8,	128(%r10)
	fnmsub 22 , 9 , 3 , 22
	fmadd  23 , 9 , 2 , 23
	fnmsub 20 , 13 , 1 , 20
	fmadd  21 , 13 , 0 , 21
	dcbt 	%r29,%r11
	fnmsub 6 , 11 , 27 , 6
	fmadd  7 , 11 , 26 , 7
	stfd  14,	0(%r13)
	lfd   9,	136(%r10)
	fmul   14 , 4 , 24
	stfd  15,	8(%r13)
	addi  %r10 , %r10 , 144
	fmul   15 , 4 , 25
	dcbt 	%r30,%r11
	fmadd  6 , 16 , 30 , 6
	fmadd  7 , 16 , 31 , 7
	stfd  22,	16(%r13)
	fnmsub 14 , 5 , 25 , 14
	lfd   10,	16(%r10)
	fmul   22 , 4 , 26
	stfd  23,	24(%r13)
	fmadd  15 , 5 , 24 , 15
	lfd   11,	24(%r10)
	fmul   23 , 4 , 27
	dcbt 	%r31,%r11
	fnmsub 6 , 17 , 31 , 6
	fmadd  14 , 18 , 28 , 14
	stfd  20,	32(%r13)
	fnmsub 22 , 5 , 27 , 22
	lfd   4,	0(%r10)
	fmadd  15 , 18 , 29 , 15
	stfd  21,	40(%r13)
	fmadd  23 , 5 , 26 , 23
	lfd   16,	64(%r10)
	fmadd  7 , 17 , 30 , 7
	dcbt 	%r3,%r11
	fnmsub 14 , 19 , 29 , 14
	addi  %r13 , %r13 , 96
	fmadd  22 , 18 , 30 , 22
	lfd   5,	8(%r10)
	fmadd  15 , 19 , 28 , 15
	lfd   17,	72(%r10)
	fmadd  23 , 18 , 31 , 23
	fmadd  6 , 12 , 2 , 6
	lfd   24,	0(%r12)
	fmadd  7 , 12 , 3 , 7
	lfd   18,	48(%r10)
	fnmsub 22 , 19 , 31 , 22
	fmadd  14 , 8 , 0 , 14
	lfd   25,	8(%r12)
	fmadd  23 , 19 , 30 , 23
	lfd   26,	16(%r12)
	fmadd  15 , 8 , 1 , 15
	fnmsub 6 , 13 , 3 , 6
	lfd   19,	56(%r10)
	fmadd  22 , 8 , 2 , 22
	lfd   27,	24(%r12)
	fmadd  7 , 13 , 2 , 7
	fmadd  23 , 8 , 3 , 23
	lfd   28,	32(%r12)
	fnmsub 14 , 9 , 1 , 14
	lfd   29,	40(%r12)
	fmadd  15 , 9 , 0 , 15
	fnmsub 22 , 9 , 3 , 22
	lfd   30,	48(%r12)
	lfd   31,	56(%r12)
	fmadd  23 , 9 , 2 , 23
	bdnz    rec_su3_dag_lab2
rec_su3_dag_lab3:
	addi  %r14 , %r12 , 128
	stfd  6,	48(%r24)
	fmul   20 , 10 , 24
	stfd  7,	56(%r24)
	fmul   21 , 10 , 25
	lfd   0,	64(%r12)
	fmul   6 , 10 , 26
	stfd  14,	64(%r24)
	fmul   7 , 10 , 27
	stfd  15,	72(%r24)
	fmul   14 , 4 , 24
	lfd   1,	72(%r12)
	fmul   15 , 4 , 25
	stfd  22,	80(%r24)
	fnmsub 20 , 11 , 25 , 20
	stfd  23,	88(%r24)
	fmul   22 , 4 , 26
	lfd   2,	80(%r12)
	fmul   23 , 4 , 27
	lfd   3,	88(%r12)
	fnmsub 14 , 5 , 25 , 14
	addi  %r24 , %r13 , 0
	fmadd  15 , 5 , 24 , 15
	lfd   8,	96(%r10)
	fmadd  21 , 11 , 24 , 21
	lfd   12,	112(%r10)
	fnmsub 22 , 5 , 27 , 22
	addi  %r12 , %r12 , 128
	fmadd  23 , 5 , 26 , 23
	lfd   9,	104(%r10)
	fmadd  14 , 18 , 28 , 14
	lfd   13,	120(%r10)
	fmadd  15 , 18 , 29 , 15
	dcbt 	%r28,%r14
	fmadd  20 , 16 , 28 , 20
	addi  %r19 , %r26 , 160
	fmadd  22 , 18 , 30 , 22
	lfd   4,	32(%r10)
	fmadd  23 , 18 , 31 , 23
	lfd   5,	40(%r10)
	fmadd  21 , 16 , 29 , 21
	dcbt 	%r29,%r14
	fnmsub 14 , 19 , 29 , 14
	addi  %r15 , %r12 , 0
	fmadd  15 , 19 , 28 , 15
	lfd   18,	80(%r10)
	fnmsub 22 , 19 , 31 , 22
	addi  %r15 , %r12 , 0
	fmadd  23 , 19 , 30 , 23
	dcbt 	%r30,%r14
	fnmsub 20 , 17 , 29 , 20
	addi  %r15 , %r12 , 0
	addi  %r14 , %r14 , 128
	lfd   19,	88(%r10)
	fmadd  21 , 17 , 28 , 21
	addi  %r15 , %r12 , 0
	fmadd  14 , 8 , 0 , 14
	dcbt 	%r28,%r14
	fmadd  15 , 8 , 1 , 15
	addi  %r20 , %r26 , 256
	fmadd  22 , 8 , 2 , 22
	addi  %r16 , %r12 , 128
	fmadd  23 , 8 , 3 , 23
	addi  %r21 , %r26 , 352
	fmadd  20 , 12 , 0 , 20
	dcbt 	%r29,%r14
	fmadd  21 , 12 , 1 , 21
	addi  %r22 , %r26 , 448
	fnmsub 14 , 9 , 1 , 14
	lfd   8,	128(%r10)
	fmadd  15 , 9 , 0 , 15
	addi  %r17 , %r16 , 128
	fnmsub 22 , 9 , 3 , 22
	dcbt 	%r30,%r14
	fmadd  23 , 9 , 2 , 23
	addi  %r18 , %r17 , 128
	addi  %r14 , %r14 , 128
	la   %r27,	-1(%r27)
	fnmsub 20 , 13 , 1 , 20
	stfd  14,	0(%r13)
	fmadd  21 , 13 , 0 , 21
	dcbt 	%r28,%r14
	fmul   14 , 4 , 24
	fnmsub 6 , 11 , 27 , 6
	stfd  15,	8(%r13)
	fmadd  7 , 11 , 26 , 7
	stfd  22,	16(%r13)
	fmul   15 , 4 , 25
	dcbt 	%r29,%r14
	fmul   22 , 4 , 26
	fnmsub 14 , 5 , 25 , 14
	stfd  23,	24(%r13)
	fmadd  6 , 16 , 30 , 6
	stfd  20,	32(%r13)
	fmul   23 , 4 , 27
	dcbt 	%r30,%r14
	fmadd  15 , 5 , 24 , 15
	addi  %r14 , %r14 , 128
	stfd  21,	40(%r13)
	fnmsub 22 , 5 , 27 , 22
	lfd   9,	136(%r10)
	fmadd  7 , 16 , 31 , 7
	dcbt 	%r28,%r14
	addi  %r10 , %r10 , 144
	addi  %r13 , %r13 , 96
	fmadd  23 , 5 , 26 , 23
	lfd   24,	0(%r19)
	fmadd  14 , 18 , 28 , 14
	lfd   25,	8(%r19)
	fmadd  15 , 18 , 29 , 15
	dcbt 	%r29,%r14
	fmadd  22 , 18 , 30 , 22
	addi  %r11 , %r10 , 0
	fnmsub 6 , 17 , 31 , 6
	lfd   26,	16(%r19)
	fmadd  23 , 18 , 31 , 23
	lfd   27,	24(%r19)
	fmadd  7 , 17 , 30 , 7
	dcbt 	%r30,%r14
	fnmsub 14 , 19 , 29 , 14
	fmadd  15 , 19 , 28 , 15
	lfd   16,	0(%r20)
	fnmsub 22 , 19 , 31 , 22
	lfd   29,	8(%r12)
	fmadd  23 , 19 , 30 , 23
	lfd   28,	0(%r12)
	fmadd  6 , 12 , 2 , 6
	dcbt 	%r29,%r11
	fmadd  7 , 12 , 3 , 7
	lfd   30,	16(%r12)
	fmadd  14 , 8 , 0 , 14
	lfd   31,	24(%r12)
	fmadd  15 , 8 , 1 , 15
	dcbt 	%r30,%r11
	fmadd  22 , 8 , 2 , 22
	lfd   20,	0(%r16)
	fmadd  23 , 8 , 3 , 23
	dcbt 	%r31,%r11
	fnmsub 6 , 13 , 3 , 6
	lfd   17,	8(%r20)
	fmadd  7 , 13 , 2 , 7
	dcbt 	%r3,%r11
	fnmsub 14 , 9 , 1 , 14
	lfd   21,	8(%r16)
	fmadd  15 , 9 , 0 , 15
	fnmsub 22 , 9 , 3 , 22
	lfd   18,	16(%r20)
	fmadd  23 , 9 , 2 , 23
	stfd  6,	48(%r24)
	fadd   0 , 28 , 24
	stfd  7,	56(%r24)
	fadd   1 , 29 , 25
	lfd   19,	24(%r20)
	fadd   2 , 30 , 26
	stfd  14,	64(%r24)
	fadd   3 , 31 , 27
	stfd  15,	72(%r24)
	fsub   4 , 27 , 31
	fsub   5 , 30 , 26
	stfd  22,	80(%r24)
	fsub   6 , 25 , 29
	stfd  23,	88(%r24)
	fsub   7 , 28 , 24
	lfd   22,	16(%r16)
	fadd   0 , 0 , 20
	addi  %r24 , %r9 , 96
	fadd   1 , 1 , 21
	lfd   23,	24(%r16)
	fadd   6 , 6 , 20
	lfd   24,	0(%r21)
	fadd   7 , 7 , 21
	fadd   2 , 2 , 22
	lfd   28,	0(%r17)
	fadd   0 , 0 , 16
	lfd   20,	0(%r18)
	fadd   3 , 3 , 23
	lfd   25,	8(%r21)
	fadd   1 , 1 , 17
	fsub   4 , 4 , 22
	lfd   29,	8(%r17)
	fadd   2 , 2 , 18
	lfd   21,	8(%r18)
	fsub   5 , 5 , 23
	lfd   26,	16(%r21)
	fadd   3 , 3 , 19
	fsub   6 , 6 , 16
	lfd   30,	16(%r17)
	fadd   4 , 4 , 18
	lfd   22,	16(%r18)
	fsub   7 , 7 , 17
	lfd   27,	24(%r21)
	fadd   5 , 5 , 19
	lfd   16,	0(%r22)
	fadd   0 , 0 , 28
	lfd   31,	24(%r17)
	fadd   1 , 1 , 29
	fadd   2 , 2 , 30
	lfd   17,	8(%r22)
	fsub   4 , 4 , 29
	lfd   18,	16(%r22)
	fadd   5 , 5 , 28
	fadd   0 , 0 , 24
	lfd   19,	24(%r22)
	fadd   3 , 3 , 31
	lfd   23,	24(%r18)
	fadd   1 , 1 , 25
	fadd   2 , 2 , 26
	lfd   28,	32(%r15)
	fadd   4 , 4 , 25
	lfd   29,	40(%r15)
	fsub   5 , 5 , 24
	fadd   3 , 3 , 27
	lfd   25,	40(%r19)
	fadd   6 , 6 , 31
	lfd   24,	32(%r19)
	fsub   7 , 7 , 30
	fadd   0 , 0 , 20
	lfd   31,	56(%r15)
	fadd   1 , 1 , 21
	lfd   30,	48(%r15)
	fadd   2 , 2 , 22
	fsub   6 , 6 , 27
	fadd   7 , 7 , 26
	fadd   3 , 3 , 23
	lfd   27,	56(%r19)
	fsub   4 , 4 , 20
	lfd   26,	48(%r19)
	fsub   5 , 5 , 21
	fsub   6 , 6 , 22
	lfd   20,	32(%r16)
	fsub   7 , 7 , 23
	fadd   0 , 0 , 16
	lfd   21,	40(%r16)
	fadd   1 , 1 , 17
	fadd   2 , 2 , 18
	lfd   22,	48(%r16)
	fadd   3 , 3 , 19
	fadd   4 , 4 , 16
	lfd   23,	56(%r16)
	fadd   5 , 5 , 17
	stfd  0,	0(%r9)
	fadd   6 , 6 , 18
	stfd  1,	8(%r9)
	fadd   7 , 7 , 19
	lfd   16,	32(%r20)
	fadd   8 , 28 , 24
	fadd   9 , 29 , 25
	lfd   17,	40(%r20)
	fadd   10 , 30 , 26
	fadd   11 , 31 , 27
	lfd   18,	48(%r20)
	fsub   12 , 27 , 31
	fsub   13 , 30 , 26
	lfd   19,	56(%r20)
	fsub   14 , 25 , 29
	fsub   15 , 28 , 24
	lfd   26,	48(%r21)
	fadd   8 , 8 , 20
	lfd   25,	40(%r21)
	fadd   9 , 9 , 21
	lfd   24,	32(%r21)
	fadd   10 , 10 , 22
	fadd   11 , 11 , 23
	lfd   28,	32(%r17)
	fsub   12 , 12 , 22
	fadd   8 , 8 , 16
	lfd   29,	40(%r17)
	fadd   9 , 9 , 17
	lfd   22,	48(%r18)
	fadd   10 , 10 , 18
	lfd   30,	48(%r17)
	fadd   11 , 11 , 19
	fsub   13 , 13 , 23
	lfd   27,	56(%r21)
	fadd   12 , 12 , 18
	fadd   14 , 14 , 20
	lfd   31,	56(%r17)
	fadd   15 , 15 , 21
	lfd   18,	48(%r22)
	fadd   8 , 8 , 28
	lfd   20,	32(%r18)
	fadd   13 , 13 , 19
	lfd   21,	40(%r18)
	fadd   9 , 9 , 29
	fsub   14 , 14 , 16
	lfd   19,	56(%r22)
	fsub   15 , 15 , 17
	lfd   23,	56(%r18)
	fadd   10 , 10 , 30
	lfd   16,	32(%r22)
	fadd   11 , 11 , 31
	lfd   17,	40(%r22)
	fadd   8 , 8 , 24
	fadd   9 , 9 , 25
	lfd   0,	64(%r20)
	fsub   12 , 12 , 29
	fadd   10 , 10 , 26
	lfd   1,	72(%r20)
	fadd   11 , 11 , 27
	lfd   29,	72(%r15)
	fadd   13 , 13 , 28
	fadd   14 , 14 , 31
	fadd   12 , 12 , 25
	lfd   28,	64(%r15)
	fsub   15 , 15 , 30
	lfd   31,	88(%r15)
	fadd   8 , 8 , 20
	lfd   25,	72(%r19)
	fsub   13 , 13 , 24
	lfd   30,	80(%r15)
	fsub   14 , 14 , 27
	fadd   9 , 9 , 21
	lfd   24,	64(%r19)
	fadd   15 , 15 , 26
	lfd   27,	88(%r19)
	fadd   10 , 10 , 22
	fadd   11 , 11 , 23
	lfd   26,	80(%r19)
	fsub   12 , 12 , 20
	fsub   13 , 13 , 21
	fsub   14 , 14 , 22
	fsub   15 , 15 , 23
	fadd   8 , 8 , 16
	fadd   9 , 9 , 17
	fadd   10 , 10 , 18
	fadd   11 , 11 , 19
	fadd   12 , 12 , 16
	fadd   13 , 13 , 17
	stfd  8,	16(%r9)
	fadd   14 , 14 , 18
	stfd  9,	24(%r9)
	fadd   15 , 15 , 19
	lfd   8,	80(%r20)
	fadd   16 , 28 , 24
	fadd   17 , 29 , 25
	lfd   9,	88(%r20)
	fadd   18 , 30 , 26
	fadd   19 , 31 , 27
	fsub   20 , 27 , 31
	fsub   21 , 30 , 26
	fsub   22 , 25 , 29
	lfd   31,	88(%r16)
	fsub   23 , 28 , 24
	lfd   30,	80(%r16)
	lfd   29,	72(%r16)
	lfd   28,	64(%r16)
	lfd   24,	64(%r21)
	fadd   19 , 19 , 31
	fadd   18 , 18 , 30
	lfd   25,	72(%r21)
	fadd   17 , 17 , 29
	fadd   16 , 16 , 28
	lfd   26,	80(%r21)
	fsub   20 , 20 , 30
	fadd   19 , 19 , 9
	lfd   27,	88(%r21)
	fadd   18 , 18 , 8
	lfd   30,	80(%r17)
	fadd   17 , 17 , 1
	fadd   16 , 16 , 0
	fsub   21 , 21 , 31
	fadd   20 , 20 , 8
	fadd   22 , 22 , 28
	lfd   31,	88(%r17)
	fadd   23 , 23 , 29
	lfd   8,	80(%r22)
	fadd   18 , 18 , 30
	lfd   28,	64(%r17)
	fadd   21 , 21 , 9
	lfd   29,	72(%r17)
	fsub   22 , 22 , 0
	lfd   9,	88(%r22)
	fsub   23 , 23 , 1
	fadd   19 , 19 , 31
	lfd   0,	64(%r22)
	fadd   16 , 16 , 28
	lfd   1,	72(%r22)
	fadd   18 , 18 , 26
	fadd   17 , 17 , 29
	fsub   20 , 20 , 29
	fadd   19 , 19 , 27
	fadd   16 , 16 , 24
	lfd   29,	72(%r18)
	fadd   21 , 21 , 28
	fadd   17 , 17 , 25
	fadd   20 , 20 , 25
	lfd   28,	64(%r18)
	fadd   22 , 22 , 31
	fsub   23 , 23 , 30
	fsub   21 , 21 , 24
	lfd   31,	88(%r18)
	fadd   17 , 17 , 29
	lfd   30,	80(%r18)
	fsub   22 , 22 , 27
	addi  %r18 , %r18 , 0
	fadd   23 , 23 , 26
	addi  %r12 , %r18 , 128
	fadd   16 , 16 , 28
	dcbt 	%r28,%r12
	fadd   19 , 19 , 31
	dcbt 	%r29,%r12
	fadd   18 , 18 , 30
	dcbt 	%r30,%r12
	fsub   20 , 20 , 28
	lfd   24,	0(%r12)
	fsub   21 , 21 , 29
	lfd   25,	8(%r12)
	fsub   22 , 22 , 30
	fsub   23 , 23 , 31
	lfd   26,	16(%r12)
	fadd   16 , 16 , 0
	lfd   27,	24(%r12)
	fadd   17 , 17 , 1
	fadd   18 , 18 , 8
	lfd   28,	32(%r12)
	fadd   19 , 19 , 9
	lfd   29,	40(%r12)
	fadd   20 , 20 , 0
	fadd   21 , 21 , 1
	stfd  16,	32(%r9)
	fadd   22 , 22 , 8
	stfd  17,	40(%r9)
	fadd   23 , 23 , 9
	lfd   16,	64(%r10)
	stfd  2,	48(%r9)
	stfd  3,	56(%r9)
	lfd   17,	72(%r10)
	stfd  10,	64(%r9)
	stfd  11,	72(%r9)
	lfd   10,	16(%r10)
	stfd  18,	80(%r9)
	stfd  19,	88(%r9)
	lfd   11,	24(%r10)
	stfd  4,	96(%r9)
	stfd  5,	104(%r9)
	lfd   4,	0(%r10)
	stfd  12,	112(%r9)
	stfd  13,	120(%r9)
	lfd   5,	8(%r10)
	stfd  20,	128(%r9)
	stfd  21,	136(%r9)
	lfd   18,	48(%r10)
	addi  %r9 , %r9 , 192
	lfd   19,	56(%r10)
	lfd   30,	48(%r12)
	lfd   31,	56(%r12)
# Instruction order barrier inserted
	or.  %r27 , %r27 , %r27
	bt gt,  rec_su3_dag_lab1
	stfd  6,	48(%r24)
	stfd  7,	56(%r24)
	stfd  14,	64(%r24)
	stfd  15,	72(%r24)
	stfd  22,	80(%r24)
	stfd  23,	88(%r24)
rec_su3_dag_lab0:
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
