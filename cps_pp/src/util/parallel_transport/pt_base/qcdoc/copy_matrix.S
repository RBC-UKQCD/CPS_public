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
.globl copy_matrix
copy_matrix:

	or    %r11 , %r3 , %r3
	la   %r1,	-736(%r1)
	or    %r10 , %r4 , %r4
	or    %r13 , %r5 , %r5
	stw   %r14,	368(%r1)
	stw   %r15,	376(%r1)
	or    %r14 , %r7 , %r7
	stw   %r16,	384(%r1)
	or    %r15 , %r6 , %r6
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
	lwz   %r12,	0(%r13)
	or.  %r12 , %r12 , %r12
	bf gt,  copy_matrix_lab0
	or    %r16 , %r10 , %r10
	li   %r17,0
	or    %r22 , %r11 , %r11
	li   %r18,32
	li   %r19,64
	li   %r20,96
	li   %r21,128
	lwz   %r9,	0(%r14)
	mtctr 	%r12
	mulli  %r9 , %r9 , 144
	add   %r10 , %r16 , %r9
	lwz   %r9,	0(%r15)
	mulli  %r9 , %r9 , 144
	add   %r11 , %r22 , %r9
	b       copy_matrix_lab1
.align 4
copy_matrix_lab1:
	or.  %r8 , %r8 , %r8
	lfd   0,	0(%r10)
	lfd   1,	8(%r10)
	lfd   2,	16(%r10)
	lfd   3,	24(%r10)
	lfd   4,	32(%r10)
	lfd   5,	40(%r10)
	lwz   %r9,	4(%r14)
	mulli  %r9 , %r9 , 144
	add   %r9 , %r16 , %r9
	bf lt,  copy_matrix_lab2
	or    %r9 , %r1 , %r1
copy_matrix_lab2:
	dcbt 	%r17,%r9
	stfd  0,	0(%r11)
	dcbt 	%r18,%r9
	stfd  1,	8(%r11)
	dcbt 	%r19,%r9
	stfd  2,	16(%r11)
	dcbt 	%r20,%r9
	stfd  3,	24(%r11)
	dcbt 	%r21,%r9
	stfd  4,	32(%r11)
	dcbt 	%r18,%r14
	stfd  5,	40(%r11)
	lfd   6,	48(%r10)
	lfd   7,	56(%r10)
	lfd   8,	64(%r10)
	lfd   9,	72(%r10)
	lfd   10,	80(%r10)
	lfd   11,	88(%r10)
	lfd   12,	96(%r10)
	stfd  6,	48(%r11)
	lfd   13,	104(%r10)
	stfd  7,	56(%r11)
	lfd   14,	112(%r10)
	stfd  8,	64(%r11)
	lfd   15,	120(%r10)
	stfd  9,	72(%r11)
	lfd   16,	128(%r10)
	stfd  10,	80(%r11)
	lfd   17,	136(%r10)
	stfd  11,	88(%r11)
	la   %r15,	4(%r15)
	stfd  12,	96(%r11)
	lwz   %r9,	0(%r15)
	stfd  13,	104(%r11)
	mulli  %r9 , %r9 , 144
	la   %r14,	4(%r14)
	stfd  14,	112(%r11)
	la   %r12,	-1(%r12)
	stfd  15,	120(%r11)
	stfd  16,	128(%r11)
	stfd  17,	136(%r11)
	add   %r11 , %r22 , %r9
	lwz   %r9,	0(%r14)
	mulli  %r9 , %r9 , 144
	add   %r10 , %r16 , %r9
	bdnz    copy_matrix_lab1
copy_matrix_lab0:
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
