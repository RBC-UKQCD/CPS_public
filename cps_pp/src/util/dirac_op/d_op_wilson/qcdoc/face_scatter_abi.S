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
.globl face_scatter
face_scatter:

	or    %r11 , %r3 , %r3
	la   %r1,	-864(%r1)
	or    %r12 , %r4 , %r4
	li   %r15,0
	or    %r13 , %r5 , %r5
	li   %r16,32
	or    %r14 , %r6 , %r6
	li   %r17,64
	or    %r22 , %r11 , %r11
	li   %r18,128
	or.  %r14 , %r14 , %r14
	li   %r19,160
	li   %r20,192
	dcbt 	%r15,%r12
	li   %r21,96
	lwz   %r10,	0(%r13)
	mulli  %r10 , %r10 , 128
	dcbt 	%r16,%r12
	add   %r11 , %r22 , %r10
	dcbt 	%r17,%r12
	bf gt,   lab0 
	mtctr 	%r14
	b        lab1 
.align 4
lab1:
	lfd   0,	0(%r12)
	lfd   1,	8(%r12)
	lfd   2,	16(%r12)
	lfd   3,	24(%r12)
	lfd   4,	32(%r12)
	lfd   5,	40(%r12)
	dcbt 	%r18,%r12
	lfd   6,	48(%r12)
	stfd  0,	0(%r11)
	lfd   7,	56(%r12)
	stfd  1,	8(%r11)
	lfd   8,	64(%r12)
	stfd  2,	16(%r11)
	lfd   9,	72(%r12)
	dcbt 	%r19,%r12
	lfd   10,	80(%r12)
	stfd  3,	24(%r11)
	lfd   11,	88(%r12)
	stfd  4,	32(%r11)
	la   %r13,	4(%r13)
	stfd  5,	40(%r11)
	lwz   %r10,	0(%r13)
	dcbt 	%r20,%r12
	mulli  %r10 , %r10 , 128
	stfd  6,	48(%r11)
	la   %r12,	128(%r12)
	stfd  7,	56(%r11)
	la   %r14,	-1(%r14)
	stfd  8,	64(%r11)
	stfd  9,	72(%r11)
	stfd  10,	80(%r11)
	stfd  11,	88(%r11)
	add   %r11 , %r22 , %r10
	bdnz     lab1 
lab0:
	la   %r1,	864(%r1)
	blr
