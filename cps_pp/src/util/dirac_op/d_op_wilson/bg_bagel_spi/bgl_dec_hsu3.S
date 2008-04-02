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
.globl bgl_dec_hsu3
bgl_dec_hsu3:

	or    %r10 , %r3 , %r3
	la   %r1,	-736(%r1)
	or    %r11 , %r4 , %r4
	or    %r13 , %r5 , %r5
	stw   %r14,	368(%r1)
	stw   %r15,	376(%r1)
	or    %r14 , %r6 , %r6
	stw   %r16,	384(%r1)
	or    %r15 , %r7 , %r7
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
	li   %r18,0
	stfd  15,	120(%r1)
	li   %r19,16
	stfd  16,	128(%r1)
	li   %r20,32
	stfd  17,	136(%r1)
	li   %r21,48
	stfd  18,	144(%r1)
	li   %r22,64
	stfd  19,	152(%r1)
	li   %r23,80
	stfd  20,	160(%r1)
	li   %r24,96
	stfd  21,	168(%r1)
	li   %r25,112
	stfd  22,	176(%r1)
	li   %r26,128
	stfd  23,	184(%r1)
	li   %r27,144
	stfd  24,	192(%r1)
	li   %r28,160
	stfd  25,	200(%r1)
	li   %r29,176
	stfd  26,	208(%r1)
	lwz   %r13,	0(%r13)
	stfd  27,	216(%r1)
	or.  %r13 , %r13 , %r13
	stfd  28,	224(%r1)
	stfd  29,	232(%r1)
	stfd  30,	240(%r1)
	stfd  31,	248(%r1)
	bf gt,  bgl_dec_hsu3_lab0
	lwz   %r12,	0(%r14)
	lfpdx   0,%r18,%r10
	lfpdx   3,%r21,%r10
	lfpdx   6,%r24,%r10
	lfpdx   9,%r27,%r10
	lfpdx   1,%r19,%r10
	lfpdx   4,%r22,%r10
	lfpdx   7,%r25,%r10
	lfpdx   10,%r28,%r10
	lfpdx   2,%r20,%r10
	lfpdx   5,%r23,%r10
	lfpdx   8,%r26,%r10
	lfpdx   11,%r29,%r10
	lfpdx   12,%r18,%r15
	mtctr 	%r13
	b       bgl_dec_hsu3_lab1
.align 4
bgl_dec_hsu3_lab1:
	addi  %r16 , %r11 , 0
	la   %r14,	4(%r14)
	fxcxnpma  26 , 12 , 9 , 0
	lfpdx   13,%r18,%r11
	fxcxnpma  27 , 12 , 10 , 1
	lfpdx   16,%r21,%r11
	fxcxnpma  28 , 12 , 11 , 2
	lfpdx   14,%r19,%r11
	fxcxnsma  20 , 12 , 9 , 0
	lfpdx   17,%r22,%r11
	fxcxnsma  21 , 12 , 10 , 1
	lfpdx   15,%r20,%r11
	fxcxnsma  22 , 12 , 11 , 2
	lfpdx   18,%r23,%r11
	fxcxnpma  29 , 12 , 6 , 3
	lfpdx   19,%r24,%r11
	fxcxnpma  30 , 12 , 7 , 4
	addi  %r16 , %r16 , 144
	fxcxnpma  31 , 12 , 8 , 5
	stfpdx  20,%r18,%r12
	fxcxnsma  23 , 12 , 6 , 3
	addi  %r10 , %r10 , 192
	fxcxnsma  24 , 12 , 7 , 4
	stfpdx  21,%r19,%r12
	fxcxnsma  25 , 12 , 8 , 5
	la   %r13,	-1(%r13)
	fxpmul  20 , 13 , 26
	stfpdx  22,%r20,%r12
	fxpmul  21 , 16 , 26
	lfpdx   12,%r25,%r11
	fxpmul  22 , 19 , 26
	stfpdx  23,%r21,%r12
	fxpmul  23 , 13 , 29
	stfpdx  24,%r22,%r12
	fxcxnsma  20 , 13 , 26 , 20
	fxcxnsma  21 , 16 , 26 , 21
	stfpdx  25,%r23,%r12
	fxpmul  24 , 16 , 29
	fxpmul  25 , 19 , 29
	lwz   %r12,	0(%r14)
	fxcxnsma  23 , 13 , 29 , 23
	addi  %r17 , %r12 , 0
	la   %r14,	4(%r14)
	fxcpmadd  20 , 14 , 27 , 20
	lfpdx   13,%r26,%r11
	fxcpmadd  21 , 17 , 27 , 21
	lwz   %r12,	0(%r14)
	fxcxnsma  24 , 16 , 29 , 24
	addi  %r11 , %r11 , 144
	fxcpmadd  23 , 14 , 30 , 23
	dcbt 	%r20,%r16
	fxcxnsma  22 , 19 , 26 , 22
	fxcxnsma  20 , 14 , 27 , 20
	la   %r14,	4(%r14)
	fxcxnsma  21 , 17 , 27 , 21
	lfpdx   16,%r19,%r11
	fxcxnsma  25 , 19 , 29 , 25
	dcbt 	%r22,%r16
	fxcxnsma  23 , 14 , 30 , 23
	fxcpmadd  24 , 17 , 30 , 24
	lfpdx   19,%r22,%r11
	fxcpmadd  20 , 15 , 28 , 20
	lfpdx   14,%r18,%r15
	fxcpmadd  21 , 18 , 28 , 21
	dcbt 	%r24,%r16
	fxcpmadd  22 , 12 , 27 , 22
	fxcpmadd  23 , 15 , 31 , 23
	lfpdx   14,%r25,%r11
	fxcpmadd  25 , 12 , 30 , 25
	fxcxnsma  20 , 15 , 28 , 20
	dcbt 	%r26,%r16
	fxcxnsma  21 , 18 , 28 , 21
	fxcxnsma  24 , 17 , 30 , 24
	addi  %r16 , %r16 , 144
	fxcxnsma  23 , 15 , 31 , 23
	fxcxnsma  22 , 12 , 27 , 22
	lfpdx   17,%r20,%r11
	fxcxnsma  25 , 12 , 30 , 25
	stfpdx  20,%r18,%r17
	fpsub  26 , 0 , 9
	lfpdx   15,%r18,%r11
	fxcpmadd  24 , 18 , 31 , 24
	stfpdx  21,%r19,%r17
	fpsub  27 , 1 , 10
	lfpdx   12,%r23,%r11
	fxcpmadd  22 , 13 , 28 , 22
	stfpdx  23,%r21,%r17
	fxcpmadd  25 , 13 , 31 , 25
	dcbt 	%r20,%r16
	fpadd  20 , 0 , 9
	fxcxnsma  24 , 18 , 31 , 24
	fpadd  21 , 1 , 10
	fxcxnsma  22 , 13 , 28 , 22
	lfpdx   18,%r21,%r11
	fxcxnsma  25 , 13 , 31 , 25
	fpsub  28 , 2 , 11
	stfpdx  20,%r18,%r12
	fpadd  29 , 3 , 6
	stfpdx  24,%r22,%r17
	fpadd  30 , 4 , 7
	lfpdx   13,%r24,%r11
	fpadd  31 , 5 , 8
	stfpdx  22,%r20,%r17
	fpsub  23 , 3 , 6
	fpadd  22 , 2 , 11
	stfpdx  25,%r23,%r17
	fpsub  24 , 4 , 7
	fpsub  25 , 5 , 8
	stfpdx  21,%r19,%r12
	fxpmul  20 , 15 , 26
	fxpmul  21 , 18 , 26
	stfpdx  23,%r21,%r12
	stfpdx  22,%r20,%r12
	fxpmul  23 , 15 , 29
	fxpmul  22 , 13 , 26
	stfpdx  24,%r22,%r12
	fxcxnsma  20 , 15 , 26 , 20
	fxcxnsma  21 , 18 , 26 , 21
	stfpdx  25,%r23,%r12
	fxpmul  24 , 18 , 29
	fxcxnsma  23 , 15 , 29 , 23
	lwz   %r12,	0(%r14)
	fxpmul  25 , 13 , 29
	addi  %r17 , %r12 , 0
	lfpdx   15,%r26,%r11
	fxcpmadd  20 , 16 , 27 , 20
	la   %r14,	4(%r14)
	fxcpmadd  21 , 19 , 27 , 21
	addi  %r11 , %r11 , 144
	fxcpmadd  23 , 16 , 30 , 23
	lwz   %r12,	0(%r14)
	fxcxnsma  24 , 18 , 29 , 24
	dcbt 	%r22,%r16
	fxcxnsma  22 , 13 , 26 , 22
	fxcxnsma  20 , 16 , 27 , 20
	la   %r14,	4(%r14)
	fxcxnsma  21 , 19 , 27 , 21
	lfpdx   18,%r19,%r11
	fxcxnsma  23 , 16 , 30 , 23
	dcbt 	%r24,%r16
	fxcxnsma  25 , 13 , 29 , 25
	fxcpmadd  24 , 19 , 30 , 24
	lfpdx   16,%r18,%r15
	fxcpmadd  20 , 17 , 28 , 20
	lfpdx   13,%r22,%r11
	fxcpmadd  21 , 12 , 28 , 21
	dcbt 	%r26,%r16
	fxcpmadd  23 , 17 , 31 , 23
	fxcpmadd  22 , 14 , 27 , 22
	addi  %r16 , %r16 , 144
	fxcpmadd  25 , 14 , 30 , 25
	fxcxnsma  20 , 17 , 28 , 20
	dcbt 	%r20,%r16
	fxcxnsma  21 , 12 , 28 , 21
	fxcxnsma  23 , 17 , 31 , 23
	fxcxnsma  24 , 19 , 30 , 24
	fxcxnsma  22 , 14 , 27 , 22
	lfpdx   17,%r18,%r11
	fxcxnsma  25 , 14 , 30 , 25
	stfpdx  20,%r18,%r17
	fxcxnpma  26 , 16 , 6 , 0
	lfpdx   19,%r20,%r11
	fxcxnpma  27 , 16 , 7 , 1
	stfpdx  21,%r19,%r17
	fxcpmadd  24 , 12 , 31 , 24
	lfpdx   14,%r23,%r11
	fxcpmadd  22 , 15 , 28 , 22
	stfpdx  23,%r21,%r17
	fxcpmadd  25 , 15 , 31 , 25
	dcbt 	%r22,%r16
	fxcxnsma  20 , 16 , 6 , 0
	fxcxnsma  21 , 16 , 7 , 1
	fxcxnsma  24 , 12 , 31 , 24
	fxcxnsma  22 , 15 , 28 , 22
	fxcxnsma  25 , 15 , 31 , 25
	lfpdx   12,%r21,%r11
	fxcxnpma  28 , 16 , 8 , 2
	stfpdx  20,%r18,%r12
	fxcxnsma  29 , 16 , 9 , 3
	lfpdx   15,%r24,%r11
	fxcxnsma  30 , 16 , 10 , 4
	stfpdx  24,%r22,%r17
	fxcxnsma  31 , 16 , 11 , 5
	stfpdx  22,%r20,%r17
	fxcxnpma  23 , 16 , 9 , 3
	fxcxnsma  22 , 16 , 8 , 2
	stfpdx  25,%r23,%r17
	fxcxnpma  24 , 16 , 10 , 4
	fxcxnpma  25 , 16 , 11 , 5
	stfpdx  21,%r19,%r12
	fxpmul  20 , 17 , 26
	fxpmul  21 , 12 , 26
	stfpdx  23,%r21,%r12
	stfpdx  22,%r20,%r12
	fxpmul  23 , 17 , 29
	lfpdx   16,%r25,%r11
	fxpmul  22 , 15 , 26
	stfpdx  24,%r22,%r12
	fxcxnsma  20 , 17 , 26 , 20
	fxcxnsma  21 , 12 , 26 , 21
	stfpdx  25,%r23,%r12
	fxpmul  24 , 12 , 29
	fxcxnsma  23 , 17 , 29 , 23
	lwz   %r12,	0(%r14)
	fxpmul  25 , 15 , 29
	addi  %r17 , %r12 , 0
	lfpdx   17,%r26,%r11
	fxcpmadd  20 , 18 , 27 , 20
	la   %r14,	4(%r14)
	fxcpmadd  21 , 13 , 27 , 21
	addi  %r11 , %r11 , 144
	fxcpmadd  23 , 18 , 30 , 23
	lwz   %r12,	0(%r14)
	fxcxnsma  24 , 12 , 29 , 24
	dcbt 	%r24,%r16
	fxcxnsma  22 , 15 , 26 , 22
	fxcxnsma  20 , 18 , 27 , 20
	la   %r14,	4(%r14)
	fxcxnsma  21 , 13 , 27 , 21
	lfpdx   12,%r19,%r11
	fxcxnsma  23 , 18 , 30 , 23
	dcbt 	%r26,%r16
	fxcxnsma  25 , 15 , 29 , 25
	fxcpmadd  24 , 13 , 30 , 24
	lfpdx   18,%r18,%r15
	fxcpmadd  20 , 19 , 28 , 20
	lfpdx   15,%r22,%r11
	fxcpmadd  21 , 14 , 28 , 21
	lfpdx   18,%r25,%r11
	fxcpmadd  23 , 19 , 31 , 23
	addi  %r16 , %r16 , 144
	fxcpmadd  22 , 16 , 27 , 22
	dcbt 	%r20,%r16
	fxcpmadd  25 , 16 , 30 , 25
	fxcxnsma  20 , 19 , 28 , 20
	fxcxnsma  21 , 14 , 28 , 21
	fxcxnsma  23 , 19 , 31 , 23
	dcbt 	%r22,%r16
	fxcxnsma  24 , 13 , 30 , 24
	fxcxnsma  22 , 16 , 27 , 22
	lfpdx   19,%r18,%r11
	fxcxnsma  25 , 16 , 30 , 25
	stfpdx  20,%r18,%r17
	fpadd  26 , 0 , 6
	lfpdx   13,%r20,%r11
	fpadd  27 , 1 , 7
	stfpdx  21,%r19,%r17
	fxcpmadd  24 , 14 , 31 , 24
	lfpdx   16,%r23,%r11
	fxcpmadd  22 , 17 , 28 , 22
	stfpdx  23,%r21,%r17
	fxcpmadd  25 , 17 , 31 , 25
	dcbt 	%r24,%r16
	fpsub  20 , 0 , 6
	fpsub  21 , 1 , 7
	fxcxnsma  24 , 14 , 31 , 24
	lfpdx   0,%r18,%r10
	fxcxnsma  22 , 17 , 28 , 22
	lfpdx   6,%r24,%r10
	fxcxnsma  25 , 17 , 31 , 25
	lfpdx   14,%r21,%r11
	fpadd  28 , 2 , 8
	stfpdx  20,%r18,%r12
	fpadd  29 , 3 , 9
	lfpdx   17,%r24,%r11
	fpadd  30 , 4 , 10
	stfpdx  24,%r22,%r17
	fpadd  31 , 5 , 11
	stfpdx  22,%r20,%r17
	fpsub  23 , 3 , 9
	lfpdx   1,%r19,%r10
	fpsub  22 , 2 , 8
	stfpdx  25,%r23,%r17
	fpsub  24 , 4 , 10
	lfpdx   3,%r21,%r10
	fpsub  25 , 5 , 11
	stfpdx  21,%r19,%r12
	fxpmul  20 , 19 , 26
	lfpdx   9,%r27,%r10
	fxpmul  21 , 14 , 26
	stfpdx  23,%r21,%r12
	stfpdx  22,%r20,%r12
	fxpmul  23 , 19 , 29
	lfpdx   4,%r22,%r10
	fxpmul  22 , 17 , 26
	stfpdx  24,%r22,%r12
	fxcxnsma  20 , 19 , 26 , 20
	lfpdx   7,%r25,%r10
	fxcxnsma  21 , 14 , 26 , 21
	stfpdx  25,%r23,%r12
	fxpmul  24 , 14 , 29
	lfpdx   10,%r28,%r10
	fxcxnsma  23 , 19 , 29 , 23
	lwz   %r12,	0(%r14)
	fxpmul  25 , 17 , 29
	lfpdx   2,%r20,%r10
	addi  %r17 , %r12 , 0
	lfpdx   19,%r26,%r11
	fxcpmadd  20 , 12 , 27 , 20
	la   %r14,	4(%r14)
	fxcpmadd  21 , 15 , 27 , 21
	addi  %r11 , %r11 , 144
	fxcpmadd  23 , 12 , 30 , 23
	lwz   %r12,	0(%r14)
	fxcxnsma  24 , 14 , 29 , 24
	dcbt 	%r26,%r16
	fxcxnsma  22 , 17 , 26 , 22
	fxcxnsma  20 , 12 , 27 , 20
	lfpdx   5,%r23,%r10
	fxcxnsma  21 , 15 , 27 , 21
	lfpdx   8,%r26,%r10
	fxcxnsma  23 , 12 , 30 , 23
	lfpdx   11,%r29,%r10
	fxcxnsma  25 , 17 , 29 , 25
	fxcpmadd  24 , 15 , 30 , 24
	lfpdx   12,%r18,%r15
	fxcpmadd  20 , 13 , 28 , 20
	fxcpmadd  21 , 16 , 28 , 21
	fxcpmadd  23 , 13 , 31 , 23
	fxcpmadd  22 , 18 , 27 , 22
	fxcpmadd  25 , 18 , 30 , 25
	fxcxnsma  20 , 13 , 28 , 20
	fxcxnsma  21 , 16 , 28 , 21
	fxcxnsma  23 , 13 , 31 , 23
	fxcxnsma  24 , 15 , 30 , 24
	fxcxnsma  22 , 18 , 27 , 22
	fxcxnsma  25 , 18 , 30 , 25
	stfpdx  20,%r18,%r17
	stfpdx  21,%r19,%r17
	fxcpmadd  24 , 16 , 31 , 24
	fxcpmadd  22 , 19 , 28 , 22
	stfpdx  23,%r21,%r17
	fxcpmadd  25 , 19 , 31 , 25
	fxcxnsma  24 , 16 , 31 , 24
	fxcxnsma  22 , 19 , 28 , 22
	fxcxnsma  25 , 19 , 31 , 25
	stfpdx  24,%r22,%r17
	stfpdx  22,%r20,%r17
	stfpdx  25,%r23,%r17
	bdnz    bgl_dec_hsu3_lab1
bgl_dec_hsu3_lab0:
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
