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
.globl QCDOC_ChRecon_p_add
QCDOC_ChRecon_p_add:

	or    %r9 , %r3 , %r3
	la   %r1,	-1056(%r1)
	or    %r18 , %r4 , %r4
	or    %r19 , %r5 , %r5
	lwz   %r19,	0(%r19)
	or.  %r19 , %r19 , %r19
	bf gt,   lab0 
	addi  %r10 , %r1 , 864
	li   %r20,0
	addi  %r13 , %r18 , 0
	li   %r21,32
	addi  %r14 , %r13 , 128
	li   %r22,64
	addi  %r15 , %r14 , 128
	dcbt 	%r20,%r13
	addi  %r16 , %r15 , 128
	dcbt 	%r21,%r13
	addi  %r18 , %r16 , 128
	addi  %r17 , %r16 , 0
	dcbt 	%r22,%r13
	dcbt 	%r20,%r14
	dcbt 	%r21,%r14
	dcbt 	%r22,%r14
	dcbt 	%r20,%r15
	dcbt 	%r21,%r15
	dcbt 	%r22,%r15
	lfd   0,	0(%r13)
	lfd   1,	8(%r13)
	lfd   2,	16(%r13)
	lfd   3,	24(%r13)
	lfd   8,	0(%r14)
	lfd   9,	8(%r14)
	lfd   10,	16(%r14)
	lfd   11,	24(%r14)
	lfd   12,	0(%r15)
	lfd   13,	8(%r15)
	lfd   14,	16(%r15)
	lfd   15,	24(%r15)
	mtctr 	%r19
	b        lab1 
.align 4
lab1:
	fsub   6 , 1 , 8
	stfd  28,	128(%r10)
	fadd   7 , 0 , 9
	lfd   16,	0(%r16)
	fadd   4 , 3 , 10
	stfd  29,	136(%r10)
	fsub   5 , 11 , 2
	lfd   17,	8(%r16)
	fadd   0 , 0 , 8
	stfd  30,	176(%r10)
	fadd   1 , 1 , 9
	lfd   18,	16(%r16)
	fneg   7 , 7
	stfd  31,	184(%r10)
	fadd   2 , 2 , 10
	lfd   19,	24(%r16)
	fadd   3 , 3 , 11
	dcbt 	%r20,%r17
	fadd   24 , 0 , 12
	lfd   20,	96(%r9)
	fadd   25 , 1 , 13
	lfd   21,	104(%r9)
	fadd   28 , 4 , 13
	lfd   22,	144(%r9)
	fadd   26 , 2 , 14
	lfd   23,	152(%r9)
	fadd   27 , 3 , 15
	dcbt 	%r21,%r17
	fsub   29 , 5 , 12
	lfd   8,	0(%r9)
	fsub   30 , 6 , 15
	lfd   9,	8(%r9)
	fadd   31 , 7 , 14
	lfd   10,	16(%r9)
	fadd   24 , 24 , 16
	lfd   11,	24(%r9)
	fadd   25 , 25 , 17
	dcbt 	%r22,%r17
	fadd   26 , 26 , 18
	lfd   12,	32(%r9)
	addi  %r17 , %r17 , 128
	lfd   13,	40(%r9)
	fadd   27 , 27 , 19
	lfd   14,	48(%r9)
	fadd   28 , 28 , 16
	lfd   15,	56(%r9)
	fadd   29 , 29 , 17
	dcbt 	%r20,%r17
	fadd   30 , 30 , 18
	addi  %r9 , %r9 , 0
	fadd   31 , 31 , 19
	lfd   0,	32(%r13)
	fadd   24 , 24 , 8
	lfd   1,	40(%r13)
	fadd   25 , 25 , 9
	lfd   2,	48(%r13)
	fadd   26 , 26 , 10
	dcbt 	%r21,%r17
	fadd   27 , 27 , 11
	lfd   3,	56(%r13)
	lfd   16,	32(%r14)
	fadd   28 , 28 , 12
	stfd  24,	0(%r9)
	fadd   29 , 29 , 13
	lfd   17,	40(%r14)
	fadd   30 , 30 , 14
	dcbt 	%r22,%r17
	fadd   31 , 31 , 15
	stfd  25,	8(%r9)
	addi  %r17 , %r17 , 128
	lfd   18,	48(%r14)
	stfd  26,	48(%r9)
	lfd   19,	56(%r14)
	stfd  27,	56(%r9)
	lfd   8,	32(%r15)
	stfd  28,	96(%r9)
	lfd   9,	40(%r15)
	fsub   6 , 1 , 16
	stfd  29,	104(%r9)
	lfd   10,	48(%r15)
	fadd   7 , 0 , 17
	stfd  30,	144(%r9)
	fadd   1 , 1 , 17
	lfd   11,	56(%r15)
	fadd   0 , 0 , 16
	stfd  31,	152(%r9)
	fadd   4 , 3 , 18
	lfd   12,	32(%r16)
	lfd   13,	40(%r16)
	fsub   5 , 19 , 2
	lfd   14,	48(%r16)
	fneg   7 , 7
	lfd   15,	56(%r16)
	fadd   2 , 2 , 18
	lfd   16,	64(%r9)
	fadd   3 , 3 , 19
	lfd   17,	72(%r9)
	fadd   24 , 0 , 8
	lfd   18,	80(%r9)
	fadd   25 , 1 , 9
	lfd   19,	88(%r9)
	fadd   28 , 4 , 9
	lfd   0,	64(%r13)
	fadd   26 , 2 , 10
	lfd   1,	72(%r13)
	fadd   27 , 3 , 11
	lfd   9,	72(%r16)
	fsub   29 , 5 , 8
	lfd   2,	80(%r13)
	fsub   30 , 6 , 11
	lfd   3,	88(%r13)
	fadd   31 , 7 , 10
	lfd   8,	64(%r16)
	fadd   24 , 24 , 12
	lfd   11,	120(%r9)
	fadd   25 , 25 , 13
	lfd   10,	112(%r9)
	fadd   26 , 26 , 14
	dcbt 	%r20,%r17
	fadd   27 , 27 , 15
	addi  %r13 , %r18 , 0
	fadd   28 , 28 , 12
	la   %r19,	-1(%r19)
	fadd   29 , 29 , 13
	fadd   30 , 30 , 14
	lfd   12,	64(%r14)
	fadd   31 , 31 , 15
	lfd   13,	72(%r14)
	fadd   24 , 24 , 16
	lfd   14,	80(%r14)
	fadd   25 , 25 , 17
	lfd   15,	88(%r14)
	fadd   26 , 26 , 18
	lfd   16,	64(%r15)
	fadd   27 , 27 , 19
	lfd   17,	72(%r15)
	fadd   28 , 28 , 20
	lfd   18,	80(%r15)
	fadd   29 , 29 , 21
	lfd   19,	88(%r15)
	fadd   30 , 30 , 10
	lfd   20,	160(%r9)
	fadd   31 , 31 , 11
	lfd   21,	168(%r9)
	addi  %r14 , %r13 , 128
	lfd   10,	80(%r16)
	addi  %r9 , %r9 , 0
	lfd   11,	88(%r16)
	fsub   6 , 1 , 12
	stfd  24,	16(%r9)
	fadd   7 , 0 , 13
	dcbt 	%r21,%r17
	fadd   4 , 3 , 14
	stfd  25,	24(%r9)
	fsub   5 , 15 , 2
	addi  %r15 , %r14 , 128
	fadd   0 , 0 , 12
	stfd  26,	64(%r9)
	fadd   1 , 1 , 13
	addi  %r16 , %r15 , 128
	fneg   7 , 7
	stfd  27,	72(%r9)
	fadd   2 , 2 , 14
	dcbt 	%r22,%r17
	fadd   3 , 3 , 15
	stfd  28,	112(%r9)
	addi  %r17 , %r17 , 128
	lfd   12,	128(%r9)
	fadd   24 , 0 , 16
	stfd  29,	120(%r9)
	fadd   25 , 1 , 17
	lfd   13,	136(%r9)
	fadd   26 , 2 , 18
	stfd  30,	160(%r9)
	fadd   27 , 3 , 19
	dcbt 	%r20,%r17
	fadd   28 , 4 , 17
	stfd  31,	168(%r9)
	fsub   29 , 5 , 16
	lfd   14,	176(%r9)
	fsub   30 , 6 , 19
	lfd   15,	184(%r9)
	fadd   31 , 7 , 18
	lfd   0,	0(%r13)
	fadd   24 , 24 , 8
	dcbt 	%r21,%r17
	fadd   25 , 25 , 9
	addi  %r9 , %r9 , 0
	fadd   26 , 26 , 10
	addi  %r10 , %r9 , 0
	fadd   27 , 27 , 11
	lfd   1,	8(%r13)
	fadd   28 , 28 , 8
	lfd   2,	16(%r13)
	fadd   29 , 29 , 9
	dcbt 	%r22,%r17
	fadd   30 , 30 , 10
	lfd   3,	24(%r13)
	addi  %r17 , %r17 , 128
	lfd   8,	0(%r14)
	fadd   31 , 31 , 11
	lfd   9,	8(%r14)
	fadd   24 , 24 , 12
	lfd   10,	16(%r14)
	fadd   25 , 25 , 13
	lfd   11,	24(%r14)
	fadd   26 , 26 , 22
	lfd   12,	0(%r15)
	fadd   27 , 27 , 23
	lfd   13,	8(%r15)
	fadd   28 , 28 , 20
	addi  %r18 , %r16 , 128
	fadd   29 , 29 , 21
	stfd  24,	32(%r10)
	fadd   30 , 30 , 14
	addi  %r9 , %r9 , 192
	fadd   31 , 31 , 15
	stfd  25,	40(%r10)
	lfd   14,	16(%r15)
	stfd  26,	80(%r10)
	lfd   15,	24(%r15)
	stfd  27,	88(%r10)
	bdnz     lab1 
	stfd  28,	128(%r10)
	stfd  29,	136(%r10)
	stfd  30,	176(%r10)
	stfd  31,	184(%r10)
lab0:
	la   %r1,	1056(%r1)
	blr
