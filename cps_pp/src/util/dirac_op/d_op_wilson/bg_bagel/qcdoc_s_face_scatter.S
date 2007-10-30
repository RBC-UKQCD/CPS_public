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
.globl qcdoc_s_face_scatter
qcdoc_s_face_scatter:

	or    %r11 , %r3 , %r3
	la   %r1,	-736(%r1)
	or    %r12 , %r4 , %r4
	or    %r13 , %r5 , %r5
	stw   %r14,	368(%r1)
	stw   %r15,	376(%r1)
	or    %r14 , %r6 , %r6
	stw   %r16,	384(%r1)
	or.  %r14 , %r14 , %r14
	stw   %r17,	392(%r1)
	stw   %r18,	400(%r1)
	stw   %r19,	408(%r1)
	stw   %r20,	416(%r1)
	stw   %r21,	424(%r1)
	stw   %r22,	432(%r1)
	or    %r21 , %r11 , %r11
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
	li   %r15,0
	li   %r16,32
	li   %r17,64
	li   %r18,128
	li   %r19,160
	li   %r20,192
	dcbt 	%r15,%r12
	lwz   %r10,	0(%r13)
	li   %r22,8
	mulli  %r10 , %r10 , 64
	li   %r23,16
	add   %r11 , %r21 , %r10
	dcbt 	%r16,%r12
	li   %r24,24
	li   %r25,40
	dcbt 	%r17,%r12
	bf gt,  qcdoc_s_face_scatter_lab0
	mtctr 	%r14
	b       qcdoc_s_face_scatter_lab1
.align 4
qcdoc_s_face_scatter_lab1:
	lfpsx   0,%r15,%r12
	lfpsx   1,%r22,%r12
	dcbt 	%r15,%r12
	lfpsx   2,%r23,%r12
	lfpsx   3,%r24,%r12
	lfpsx   4,%r16,%r12
	stfpsx  0,%r15,%r11
	stfpsx  1,%r22,%r11
	lfpsx   5,%r25,%r12
	stfpsx  2,%r23,%r11
	dcbt 	%r16,%r12
	stfpsx  3,%r24,%r11
	stfpsx  4,%r16,%r11
	la   %r13,	4(%r13)
	stfpsx  5,%r25,%r11
	lwz   %r10,	0(%r13)
	la   %r14,	-1(%r14)
	mulli  %r10 , %r10 , 64
	add   %r11 , %r21 , %r10
	dcbt 	%r17,%r12
	la   %r12,	64(%r12)
	bdnz    qcdoc_s_face_scatter_lab1
qcdoc_s_face_scatter_lab0:
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
