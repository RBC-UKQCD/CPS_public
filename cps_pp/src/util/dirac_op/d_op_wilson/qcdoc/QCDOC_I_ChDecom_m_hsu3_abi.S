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
.globl QCDOC_ChDecom_m_hsu3
QCDOC_ChDecom_m_hsu3:

	or    %r9 , %r3 , %r3
	la   %r1,	-1408(%r1)
	or    %r11 , %r4 , %r4
	li   %r26,0
	addi  %r14 , %r1 , 864
	li   %r27,32
	or    %r16 , %r5 , %r5
	li   %r28,64
	or    %r19 , %r6 , %r6
	lwz   %r16,	0(%r16)
	dcbt 	%r26,%r19
	or.  %r16 , %r16 , %r16
	bf gt,   lab0 
	addi  %r24 , %r9 , 96
	lwz   %r13,	0(%r19)
	addi  %r20 , %r14 , 160
	dcbt 	%r26,%r9
	addi  %r21 , %r14 , 256
	dcbt 	%r27,%r9
	addi  %r22 , %r14 , 352
	dcbt 	%r28,%r9
	addi  %r23 , %r14 , 448
	dcbt 	%r26,%r24
	dcbt 	%r27,%r24
	dcbt 	%r28,%r24
	dcbt 	%r26,%r11
	dcbt 	%r27,%r11
	dcbt 	%r28,%r11
	lfd   0,	0(%r9)
	lfd   1,	8(%r9)
	lfd   2,	16(%r9)
	lfd   3,	24(%r9)
	lfd   4,	32(%r9)
	lfd   5,	40(%r9)
	lfd   6,	48(%r9)
	lfd   7,	56(%r9)
	lfd   8,	64(%r9)
	lfd   9,	72(%r9)
	lfd   10,	80(%r9)
	lfd   11,	88(%r9)
	lfd   12,	96(%r9)
	lfd   13,	104(%r9)
	lfd   14,	112(%r9)
	lfd   15,	120(%r9)
	lfd   16,	128(%r9)
	lfd   17,	136(%r9)
	lfd   18,	144(%r9)
	lfd   19,	152(%r9)
	lfd   20,	160(%r9)
	lfd   21,	168(%r9)
	lfd   22,	176(%r9)
	lfd   23,	184(%r9)
	stw   %r20,	0(%r14)
	stw   %r21,	4(%r14)
	stw   %r22,	8(%r14)
	stw   %r23,	12(%r14)
	mtctr 	%r16
	addi  %r14 , %r1 , 864
	b        lab1 
.align 4
lab1:
	fsub   24 , 12 , 7
	dcbt 	%r27,%r19
	fadd   25 , 13 , 6
	addi  %r9 , %r9 , 192
	fsub   26 , 18 , 1
	stw   %r9,	16(%r14)
	fadd   27 , 19 , 0
	addi  %r15 , %r14 , 0
	fsub   28 , 14 , 9
	li   %r18,3
	fadd   29 , 15 , 8
	stfd  24,	0(%r20)
	fsub   30 , 20 , 3
	stfd  25,	8(%r20)
	fadd   31 , 21 , 2
	stfd  26,	16(%r20)
	fsub   24 , 16 , 11
	stfd  27,	24(%r20)
	fadd   25 , 17 , 10
	stfd  28,	32(%r20)
	fsub   26 , 22 , 5
	stfd  29,	40(%r20)
	fadd   27 , 23 , 4
	stfd  30,	48(%r20)
	fsub   28 , 12 , 6
	stfd  31,	56(%r20)
	fsub   29 , 13 , 7
	stfd  24,	64(%r20)
	fadd   30 , 18 , 0
	stfd  25,	72(%r20)
	fadd   31 , 19 , 1
	stfd  26,	80(%r20)
	fsub   24 , 14 , 8
	stfd  27,	88(%r20)
	fsub   25 , 15 , 9
	stfd  28,	0(%r21)
	fadd   26 , 20 , 2
	stfd  29,	8(%r21)
	fadd   27 , 21 , 3
	stfd  30,	16(%r21)
	fsub   28 , 16 , 10
	stfd  31,	24(%r21)
	fsub   29 , 17 , 11
	stfd  24,	32(%r21)
	fadd   30 , 22 , 4
	stfd  25,	40(%r21)
	fadd   31 , 23 , 5
	stfd  26,	48(%r21)
	fsub   24 , 12 , 1
	stfd  27,	56(%r21)
	fadd   25 , 13 , 0
	stfd  28,	64(%r21)
	fadd   26 , 18 , 7
	stfd  29,	72(%r21)
	fsub   27 , 19 , 6
	stfd  30,	80(%r21)
	fsub   28 , 14 , 3
	stfd  31,	88(%r21)
	fadd   29 , 15 , 2
	lwz   %r12,	0(%r15)
	fadd   30 , 20 , 9
	stfd  24,	0(%r22)
	fsub   31 , 21 , 8
	addi  %r13 , %r23 , 0
	fsub   24 , 16 , 5
	stfd  25,	8(%r22)
	addi  %r19 , %r19 , -4
	addi  %r24 , %r11 , 96
	fadd   25 , 17 , 4
	stfd  26,	16(%r22)
	mtctr 	%r18
	fadd   26 , 22 , 11
	stfd  27,	24(%r22)
	fsub   27 , 23 , 10
	stfd  28,	32(%r22)
	fsub   28 , 12 , 0
	stfd  29,	40(%r22)
	fsub   29 , 13 , 1
	stfd  30,	48(%r22)
	lfd   12,	0(%r11)
	fsub   30 , 18 , 6
	stfd  31,	56(%r22)
	lfd   13,	8(%r11)
	fsub   31 , 19 , 7
	stfd  24,	64(%r22)
	lfd   18,	48(%r11)
	fsub   24 , 14 , 2
	stfd  25,	72(%r22)
	lfd   19,	56(%r11)
	fsub   25 , 15 , 3
	stfd  26,	80(%r22)
	lfd   0,	0(%r12)
	fsub   26 , 20 , 8
	stfd  27,	88(%r22)
	lfd   1,	8(%r12)
	fsub   27 , 21 , 9
	stfd  28,	0(%r23)
	lfd   2,	16(%r12)
	fsub   28 , 16 , 4
	stfd  29,	8(%r23)
	lfd   3,	24(%r12)
	fsub   29 , 17 , 5
	stfd  30,	16(%r23)
	lfd   4,	32(%r12)
	fsub   30 , 22 , 10
	stfd  31,	24(%r23)
	lfd   5,	40(%r12)
	fsub   31 , 23 , 11
	stfd  24,	32(%r23)
	lfd   6,	48(%r12)
	stfd  25,	40(%r23)
	lfd   7,	56(%r12)
	stfd  26,	48(%r23)
	stfd  27,	56(%r23)
	b        lab2 
.align 4
lab2:
	fmul   24 , 12 , 0
	stfd  27,	56(%r13)
	fmul   25 , 12 , 1
	la   %r19,	4(%r19)
	fmul   26 , 12 , 2
	stfd  28,	64(%r13)
	fmul   27 , 12 , 3
	lfd   14,	16(%r11)
	fmul   28 , 18 , 0
	stfd  29,	72(%r13)
	fmadd  24 , 13 , 1 , 24
	lfd   20,	64(%r11)
	fmul   29 , 18 , 1
	stfd  30,	80(%r13)
	fnmsub 25 , 13 , 0 , 25
	lfd   15,	24(%r11)
	fmadd  26 , 13 , 3 , 26
	stfd  31,	88(%r13)
	fnmsub 27 , 13 , 2 , 27
	or    %r10 , %r11 , %r11
	fmadd  28 , 19 , 1 , 28
	lwz   %r13,	0(%r19)
	fnmsub 29 , 19 , 0 , 29
	lfd   21,	72(%r11)
	fmadd  24 , 14 , 4 , 24
	lfd   16,	32(%r11)
	fmadd  25 , 14 , 5 , 25
	lfd   22,	80(%r11)
	fmadd  26 , 14 , 6 , 26
	lfd   17,	40(%r11)
	fmadd  27 , 14 , 7 , 27
	lfd   23,	88(%r11)
	fmadd  28 , 20 , 4 , 28
	dcbt 	%r26,%r24
	fmadd  29 , 20 , 5 , 29
	addi  %r11 , %r11 , 144
	fmadd  24 , 15 , 5 , 24
	fnmsub 25 , 15 , 4 , 25
	lfd   8,	64(%r12)
	fmadd  26 , 15 , 7 , 26
	lfd   9,	72(%r12)
	fnmsub 27 , 15 , 6 , 27
	lfd   10,	80(%r12)
	fmadd  28 , 21 , 5 , 28
	dcbt 	%r27,%r24
	fnmsub 29 , 21 , 4 , 29
	fmul   30 , 18 , 2
	fmul   31 , 18 , 3
	lfd   11,	88(%r12)
	fmadd  24 , 16 , 8 , 24
	lfd   12,	96(%r10)
	fmadd  25 , 16 , 9 , 25
	lwz   %r12,	4(%r15)
	fmadd  14 , 16 , 10 , 26
	dcbt 	%r28,%r24
	fmadd  28 , 22 , 8 , 28
	addi  %r15 , %r15 , 4
	fmadd  29 , 22 , 9 , 29
	addi  %r24 , %r24 , 96
	fmadd  24 , 17 , 9 , 24
	lfd   13,	104(%r10)
	fmadd  15 , 16 , 11 , 27
	lfd   18,	48(%r11)
	fnmsub 25 , 17 , 8 , 25
	la   %r18,	-1(%r18)
	fmadd  14 , 17 , 11 , 14
	lfd   16,	128(%r10)
	fmadd  28 , 23 , 9 , 28
	fnmsub 29 , 23 , 8 , 29
	stfd  24,	0(%r13)
	fnmsub 15 , 17 , 10 , 15
	fmul   24 , 12 , 0
	stfd  25,	8(%r13)
	fmul   26 , 12 , 2
	lfd   17,	136(%r10)
	fmul   25 , 12 , 1
	stfd  14,	16(%r13)
	fmul   27 , 12 , 3
	fmadd  30 , 19 , 3 , 30
	stfd  15,	24(%r13)
	fnmsub 31 , 19 , 2 , 31
	lfd   14,	112(%r10)
	fmadd  24 , 13 , 1 , 24
	stfd  28,	32(%r13)
	fnmsub 25 , 13 , 0 , 25
	lfd   15,	120(%r10)
	fmadd  26 , 13 , 3 , 26
	stfd  29,	40(%r13)
	fnmsub 27 , 13 , 2 , 27
	dcbt 	%r26,%r24
	fmadd  1 , 20 , 6 , 30
	fmadd  0 , 20 , 7 , 31
	fmadd  24 , 14 , 4 , 24
	lfd   12,	0(%r11)
	fmadd  25 , 14 , 5 , 25
	lfd   13,	8(%r11)
	fmadd  3 , 14 , 6 , 26
	lfd   19,	56(%r11)
	fmadd  2 , 14 , 7 , 27
	dcbt 	%r27,%r24
	fmadd  1 , 21 , 7 , 1
	fnmsub 0 , 21 , 6 , 0
	fmadd  24 , 15 , 5 , 24
	fnmsub 25 , 15 , 4 , 25
	fmadd  3 , 15 , 7 , 3
	lfd   5,	40(%r12)
	fnmsub 2 , 15 , 6 , 2
	dcbt 	%r28,%r24
	fmadd  26 , 22 , 10 , 1
	fmadd  27 , 22 , 11 , 0
	addi  %r24 , %r24 , 96
	fmadd  28 , 16 , 8 , 24
	lfd   1,	8(%r12)
	fmadd  29 , 16 , 9 , 25
	lfd   0,	0(%r12)
	fmadd  30 , 16 , 10 , 3
	lfd   4,	32(%r12)
	fmadd  31 , 16 , 11 , 2
	lfd   6,	48(%r12)
	fmadd  26 , 23 , 11 , 26
	lfd   3,	24(%r12)
	fnmsub 27 , 23 , 10 , 27
	lfd   2,	16(%r12)
	fmadd  28 , 17 , 9 , 28
	fnmsub 29 , 17 , 8 , 29
	lfd   7,	56(%r12)
	fmadd  30 , 17 , 11 , 30
	fnmsub 31 , 17 , 10 , 31
	stfd  26,	48(%r13)
	bdnz     lab2 
	fmul   24 , 12 , 0
	stfd  27,	56(%r13)
	fmul   25 , 12 , 1
	la   %r19,	4(%r19)
	fmul   26 , 12 , 2
	stfd  28,	64(%r13)
	fmul   27 , 12 , 3
	lfd   14,	16(%r11)
	fmul   28 , 18 , 0
	stfd  29,	72(%r13)
	fmadd  24 , 13 , 1 , 24
	lfd   20,	64(%r11)
	fmul   29 , 18 , 1
	stfd  30,	80(%r13)
	fnmsub 25 , 13 , 0 , 25
	lfd   15,	24(%r11)
	fmadd  26 , 13 , 3 , 26
	stfd  31,	88(%r13)
	fnmsub 27 , 13 , 2 , 27
	or    %r10 , %r11 , %r11
	fmadd  28 , 19 , 1 , 28
	lwz   %r13,	0(%r19)
	fnmsub 29 , 19 , 0 , 29
	lfd   21,	72(%r11)
	fmadd  24 , 14 , 4 , 24
	lfd   16,	32(%r11)
	fmadd  25 , 14 , 5 , 25
	lfd   22,	80(%r11)
	fmadd  26 , 14 , 6 , 26
	lfd   17,	40(%r11)
	fmadd  27 , 14 , 7 , 27
	lfd   23,	88(%r11)
	fmadd  28 , 20 , 4 , 28
	dcbt 	%r26,%r9
	fmadd  29 , 20 , 5 , 29
	addi  %r11 , %r11 , 144
	fmadd  24 , 15 , 5 , 24
	addi  %r24 , %r9 , 96
	fnmsub 25 , 15 , 4 , 25
	lfd   8,	64(%r12)
	fmadd  26 , 15 , 7 , 26
	lfd   9,	72(%r12)
	fnmsub 27 , 15 , 6 , 27
	lfd   10,	80(%r12)
	fmadd  28 , 21 , 5 , 28
	dcbt 	%r27,%r9
	fnmsub 29 , 21 , 4 , 29
	fmul   30 , 18 , 2
	la   %r19,	4(%r19)
	fmul   31 , 18 , 3
	lfd   11,	88(%r12)
	fmadd  24 , 16 , 8 , 24
	lfd   12,	96(%r10)
	fmadd  25 , 16 , 9 , 25
	lwz   %r12,	4(%r15)
	fmadd  14 , 16 , 10 , 26
	dcbt 	%r28,%r9
	fmadd  28 , 22 , 8 , 28
	addi  %r15 , %r15 , 4
	fmadd  29 , 22 , 9 , 29
	addi  %r9 , %r9 , 0
	fmadd  24 , 17 , 9 , 24
	lfd   13,	104(%r10)
	fmadd  15 , 16 , 11 , 27
	lfd   18,	144(%r9)
	fnmsub 25 , 17 , 8 , 25
	la   %r16,	-1(%r16)
	fmadd  14 , 17 , 11 , 14
	lfd   16,	128(%r10)
	fmadd  28 , 23 , 9 , 28
	fnmsub 29 , 23 , 8 , 29
	stfd  24,	0(%r13)
	fnmsub 15 , 17 , 10 , 15
	fmul   24 , 12 , 0
	stfd  25,	8(%r13)
	fmul   26 , 12 , 2
	lfd   17,	136(%r10)
	fmul   25 , 12 , 1
	stfd  14,	16(%r13)
	fmul   27 , 12 , 3
	fmadd  30 , 19 , 3 , 30
	stfd  15,	24(%r13)
	fnmsub 31 , 19 , 2 , 31
	lfd   14,	112(%r10)
	fmadd  24 , 13 , 1 , 24
	stfd  28,	32(%r13)
	fnmsub 25 , 13 , 0 , 25
	lfd   15,	120(%r10)
	fmadd  26 , 13 , 3 , 26
	stfd  29,	40(%r13)
	fnmsub 27 , 13 , 2 , 27
	dcbt 	%r26,%r24
	fmadd  1 , 20 , 6 , 30
	fmadd  0 , 20 , 7 , 31
	lfd   12,	96(%r9)
	fmadd  24 , 14 , 4 , 24
	lfd   13,	104(%r9)
	fmadd  25 , 14 , 5 , 25
	lfd   19,	152(%r9)
	fmadd  3 , 14 , 6 , 26
	lfd   20,	160(%r9)
	fmadd  2 , 14 , 7 , 27
	dcbt 	%r27,%r24
	fmadd  1 , 21 , 7 , 1
	fnmsub 0 , 21 , 6 , 0
	lfd   14,	112(%r9)
	fmadd  24 , 15 , 5 , 24
	fnmsub 25 , 15 , 4 , 25
	lfd   21,	168(%r9)
	fmadd  3 , 15 , 7 , 3
	lfd   5,	40(%r9)
	fnmsub 2 , 15 , 6 , 2
	dcbt 	%r28,%r24
	fmadd  26 , 22 , 10 , 1
	fmadd  27 , 22 , 11 , 0
	lfd   4,	32(%r9)
	fmadd  28 , 16 , 8 , 24
	lfd   1,	8(%r9)
	fmadd  29 , 16 , 9 , 25
	lfd   0,	0(%r9)
	fmadd  30 , 16 , 10 , 3
	lfd   6,	48(%r9)
	fmadd  31 , 16 , 11 , 2
	lfd   7,	56(%r9)
	fmadd  26 , 23 , 11 , 26
	lfd   3,	24(%r9)
	fnmsub 27 , 23 , 10 , 27
	lfd   2,	16(%r9)
	fmadd  28 , 17 , 9 , 28
	lfd   15,	120(%r9)
	fnmsub 29 , 17 , 8 , 29
	lfd   16,	128(%r9)
	fmadd  30 , 17 , 11 , 30
	lfd   9,	72(%r9)
	fnmsub 31 , 17 , 10 , 31
	stfd  26,	48(%r13)
	addi  %r24 , %r11 , 0
	lfd   8,	64(%r9)
	stfd  27,	56(%r13)
	lfd   10,	80(%r9)
	stfd  28,	64(%r13)
	lfd   11,	88(%r9)
	stfd  29,	72(%r13)
	lfd   17,	136(%r9)
	stfd  30,	80(%r13)
	lfd   22,	176(%r9)
	stfd  31,	88(%r13)
	lfd   23,	184(%r9)
	lwz   %r13,	0(%r19)
# Instruction order barrier inserted
	or.  %r16 , %r16 , %r16
	bt gt,   lab1 
lab0:
	la   %r1,	1408(%r1)
	blr
