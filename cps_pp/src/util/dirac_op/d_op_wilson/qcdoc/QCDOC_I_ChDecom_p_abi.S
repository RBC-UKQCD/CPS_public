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
.globl QCDOC_ChDecom_p
QCDOC_ChDecom_p:

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
	fsub   24 , 0 , 19
	addi  %r9 , %r9 , 192
	fadd   25 , 1 , 18
	lwz   %r13,	0(%r19)
	fsub   26 , 6 , 13
	dcbt 	%r27,%r19
	fadd   27 , 7 , 12
	stw   %r9,	16(%r14)
	fsub   28 , 2 , 21
	addi  %r21 , %r13 , 0
	fadd   29 , 3 , 20
	stfd  24,	0(%r20)
	fsub   30 , 8 , 15
	la   %r19,	4(%r19)
	fadd   31 , 9 , 14
	stfd  25,	8(%r20)
	fsub   24 , 4 , 23
	lwz   %r13,	0(%r19)
	fadd   25 , 5 , 22
	stfd  26,	16(%r20)
	addi  %r22 , %r13 , 0
	la   %r19,	4(%r19)
	fsub   26 , 10 , 17
	stfd  27,	24(%r20)
	addi  %r24 , %r9 , 192
	lwz   %r13,	0(%r19)
	fadd   27 , 11 , 16
	stfd  28,	32(%r20)
	addi  %r23 , %r13 , 0
	la   %r19,	4(%r19)
	fsub   28 , 0 , 18
	stfd  29,	40(%r20)
	lwz   %r13,	0(%r19)
	fsub   29 , 1 , 19
	stfd  30,	48(%r20)
	dcbt 	%r26,%r24
	fadd   30 , 6 , 12
	stfd  31,	56(%r20)
	la   %r16,	-1(%r16)
	fadd   31 , 7 , 13
	stfd  24,	64(%r20)
	fsub   24 , 2 , 20
	stfd  25,	72(%r20)
	fsub   25 , 3 , 21
	stfd  26,	80(%r20)
	dcbt 	%r27,%r24
	fadd   26 , 8 , 14
	stfd  27,	88(%r20)
	fadd   27 , 9 , 15
	stfd  28,	0(%r21)
	fsub   28 , 4 , 22
	stfd  29,	8(%r21)
	fsub   29 , 5 , 23
	stfd  30,	16(%r21)
	dcbt 	%r28,%r24
	fadd   30 , 10 , 16
	stfd  31,	24(%r21)
	addi  %r24 , %r24 , 96
	fadd   31 , 11 , 17
	stfd  24,	32(%r21)
	fsub   24 , 0 , 13
	stfd  25,	40(%r21)
	fadd   25 , 1 , 12
	stfd  26,	48(%r21)
	dcbt 	%r26,%r24
	fadd   26 , 6 , 19
	stfd  27,	56(%r21)
	fsub   27 , 7 , 18
	stfd  28,	64(%r21)
	fsub   28 , 2 , 15
	stfd  29,	72(%r21)
	fadd   29 , 3 , 14
	stfd  30,	80(%r21)
	dcbt 	%r27,%r24
	fadd   30 , 8 , 21
	stfd  31,	88(%r21)
	fsub   31 , 9 , 20
	stfd  24,	0(%r22)
	fsub   24 , 4 , 17
	stfd  25,	8(%r22)
	fadd   25 , 5 , 16
	stfd  26,	16(%r22)
	dcbt 	%r28,%r24
	fadd   26 , 10 , 23
	stfd  27,	24(%r22)
	fsub   27 , 11 , 22
	stfd  28,	32(%r22)
	fadd   28 , 0 , 12
	stfd  29,	40(%r22)
	fadd   29 , 1 , 13
	stfd  30,	48(%r22)
	lfd   0,	0(%r9)
	fadd   30 , 6 , 18
	stfd  31,	56(%r22)
	lfd   1,	8(%r9)
	fadd   31 , 7 , 19
	stfd  24,	64(%r22)
	lfd   6,	48(%r9)
	fadd   24 , 2 , 14
	stfd  25,	72(%r22)
	lfd   7,	56(%r9)
	fadd   25 , 3 , 15
	stfd  26,	80(%r22)
	lfd   2,	16(%r9)
	fadd   26 , 8 , 20
	stfd  27,	88(%r22)
	lfd   3,	24(%r9)
	fadd   27 , 9 , 21
	stfd  28,	0(%r23)
	lfd   8,	64(%r9)
	fadd   28 , 4 , 16
	stfd  29,	8(%r23)
	lfd   9,	72(%r9)
	fadd   29 , 5 , 17
	stfd  30,	16(%r23)
	lfd   4,	32(%r9)
	fadd   30 , 10 , 22
	stfd  31,	24(%r23)
	lfd   5,	40(%r9)
	fadd   31 , 11 , 23
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
