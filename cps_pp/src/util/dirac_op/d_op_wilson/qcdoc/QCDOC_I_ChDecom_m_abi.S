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
.globl QCDOC_ChDecom_m
QCDOC_ChDecom_m:

	or    %r9 , %r3 , %r3
	la   %r1,	-928(%r1)
	or    %r16 , %r4 , %r4
	li   %r26,0
	addi  %r14 , %r1 , 864
	li   %r27,32
	or    %r19 , %r5 , %r5
	li   %r28,64
	dcbt 	%r26,%r19
	lwz   %r16,	0(%r16)
	or.  %r16 , %r16 , %r16
	bf gt,   lab0 
	addi  %r24 , %r9 , 96
	lwz   %r13,	0(%r19)
	addi  %r14 , %r1 , 864
	dcbt 	%r26,%r9
	dcbt 	%r27,%r9
	dcbt 	%r28,%r9
	dcbt 	%r26,%r24
	dcbt 	%r27,%r24
	dcbt 	%r28,%r24
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
	mtctr 	%r16
	b        lab1 
.align 4
lab1:
	addi  %r20 , %r13 , 0
	la   %r19,	4(%r19)
	fsub   24 , 12 , 7
	addi  %r9 , %r9 , 192
	fadd   25 , 13 , 6
	lwz   %r13,	0(%r19)
	fsub   26 , 18 , 1
	dcbt 	%r27,%r19
	fadd   27 , 19 , 0
	stw   %r9,	16(%r14)
	fsub   28 , 14 , 9
	addi  %r21 , %r13 , 0
	fadd   29 , 15 , 8
	stfd  24,	0(%r20)
	fsub   30 , 20 , 3
	la   %r19,	4(%r19)
	fadd   31 , 21 , 2
	stfd  25,	8(%r20)
	fsub   24 , 16 , 11
	lwz   %r13,	0(%r19)
	fadd   25 , 17 , 10
	stfd  26,	16(%r20)
	addi  %r22 , %r13 , 0
	la   %r19,	4(%r19)
	fsub   26 , 22 , 5
	stfd  27,	24(%r20)
	addi  %r24 , %r9 , 192
	lwz   %r13,	0(%r19)
	fadd   27 , 23 , 4
	stfd  28,	32(%r20)
	addi  %r23 , %r13 , 0
	la   %r19,	4(%r19)
	fsub   28 , 12 , 6
	stfd  29,	40(%r20)
	lwz   %r13,	0(%r19)
	fsub   29 , 13 , 7
	stfd  30,	48(%r20)
	dcbt 	%r26,%r24
	fadd   30 , 18 , 0
	stfd  31,	56(%r20)
	la   %r16,	-1(%r16)
	fadd   31 , 19 , 1
	stfd  24,	64(%r20)
	fsub   24 , 14 , 8
	stfd  25,	72(%r20)
	fsub   25 , 15 , 9
	stfd  26,	80(%r20)
	dcbt 	%r27,%r24
	fadd   26 , 20 , 2
	stfd  27,	88(%r20)
	fadd   27 , 21 , 3
	stfd  28,	0(%r21)
	fsub   28 , 16 , 10
	stfd  29,	8(%r21)
	fsub   29 , 17 , 11
	stfd  30,	16(%r21)
	dcbt 	%r28,%r24
	fadd   30 , 22 , 4
	stfd  31,	24(%r21)
	addi  %r24 , %r24 , 96
	fadd   31 , 23 , 5
	stfd  24,	32(%r21)
	fsub   24 , 12 , 1
	stfd  25,	40(%r21)
	fadd   25 , 13 , 0
	stfd  26,	48(%r21)
	dcbt 	%r26,%r24
	fadd   26 , 18 , 7
	stfd  27,	56(%r21)
	fsub   27 , 19 , 6
	stfd  28,	64(%r21)
	fsub   28 , 14 , 3
	stfd  29,	72(%r21)
	fadd   29 , 15 , 2
	stfd  30,	80(%r21)
	dcbt 	%r27,%r24
	fadd   30 , 20 , 9
	stfd  31,	88(%r21)
	fsub   31 , 21 , 8
	stfd  24,	0(%r22)
	fsub   24 , 16 , 5
	stfd  25,	8(%r22)
	fadd   25 , 17 , 4
	stfd  26,	16(%r22)
	dcbt 	%r28,%r24
	fadd   26 , 22 , 11
	stfd  27,	24(%r22)
	fsub   27 , 23 , 10
	stfd  28,	32(%r22)
	fsub   28 , 12 , 0
	stfd  29,	40(%r22)
	fsub   29 , 13 , 1
	stfd  30,	48(%r22)
	lfd   0,	0(%r9)
	fsub   30 , 18 , 6
	stfd  31,	56(%r22)
	lfd   1,	8(%r9)
	fsub   31 , 19 , 7
	stfd  24,	64(%r22)
	lfd   6,	48(%r9)
	fsub   24 , 14 , 2
	stfd  25,	72(%r22)
	lfd   7,	56(%r9)
	fsub   25 , 15 , 3
	stfd  26,	80(%r22)
	lfd   2,	16(%r9)
	fsub   26 , 20 , 8
	stfd  27,	88(%r22)
	lfd   3,	24(%r9)
	fsub   27 , 21 , 9
	stfd  28,	0(%r23)
	lfd   8,	64(%r9)
	fsub   28 , 16 , 4
	stfd  29,	8(%r23)
	lfd   9,	72(%r9)
	fsub   29 , 17 , 5
	stfd  30,	16(%r23)
	lfd   4,	32(%r9)
	fsub   30 , 22 , 10
	stfd  31,	24(%r23)
	lfd   5,	40(%r9)
	fsub   31 , 23 , 11
	stfd  24,	32(%r23)
	lfd   10,	80(%r9)
	stfd  25,	40(%r23)
	lfd   11,	88(%r9)
	stfd  26,	48(%r23)
	lfd   12,	96(%r9)
	stfd  27,	56(%r23)
	lfd   13,	104(%r9)
	stfd  28,	64(%r23)
	lfd   14,	112(%r9)
	stfd  29,	72(%r23)
	lfd   15,	120(%r9)
	stfd  30,	80(%r23)
	lfd   16,	128(%r9)
	stfd  31,	88(%r23)
	lfd   17,	136(%r9)
	lfd   18,	144(%r9)
	lfd   19,	152(%r9)
	lfd   20,	160(%r9)
	lfd   21,	168(%r9)
	lfd   22,	176(%r9)
	lfd   23,	184(%r9)
	bdnz     lab1 
lab0:
	la   %r1,	928(%r1)
	blr
