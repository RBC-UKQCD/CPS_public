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
.globl s_bgl_rec_su3
s_bgl_rec_su3:

	or    %r10 , %r3 , %r3
	la   %r1,	-736(%r1)
	or    %r11 , %r4 , %r4
	or    %r12 , %r5 , %r5
	or    %r13 , %r6 , %r6
	stw   %r14,	368(%r1)
	stw   %r15,	376(%r1)
	or    %r14 , %r7 , %r7
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
	li   %r16,16
	li   %r17,32
	li   %r18,48
	li   %r19,64
	li   %r20,80
	li   %r21,96
	li   %r22,112
	li   %r23,128
	li   %r24,144
	li   %r25,160
	li   %r26,176
	li   %r27,8
	li   %r28,24
	li   %r29,40
	lwz   %r13,	0(%r13)
	or.  %r13 , %r13 , %r13
	bf gt,  s_bgl_rec_su3_lab0
	mtctr 	%r13
	b       s_bgl_rec_su3_lab1
.align 4
s_bgl_rec_su3_lab1:
	lfpdx   24,%r15,%r11
	lfpsx   12,%r15,%r12
	lfpsx   15,%r28,%r12
	lfpdx   27,%r16,%r11
	lfpdx   25,%r18,%r11
	lfpsx   13,%r27,%r12
	lfpsx   16,%r17,%r12
	fxpmul  18 , 24 , 12
	lfpdx   28,%r19,%r11
	fxpmul  21 , 24 , 15
	lfpdx   26,%r21,%r11
	fxpmul  19 , 27 , 12
	lfpsx   14,%r16,%r12
	fxpmul  22 , 27 , 15
	lfpsx   17,%r29,%r12
	lfpdx   29,%r22,%r11
	fxcxnpma  18 , 24 , 12 , 18
	lfpdx   30,%r17,%r11
	fxcxnpma  21 , 24 , 15 , 21
	lfpdx   31,%r20,%r11
	fxcxnpma  19 , 27 , 12 , 19
	addi  %r12 , %r12 , 48
	fxcxnpma  22 , 27 , 15 , 22
	lfpdx   24,%r23,%r11
	dcbt 	%r17,%r12
	fxcpmadd  18 , 25 , 13 , 18
	addi  %r11 , %r11 , 144
	fxcpmadd  21 , 25 , 16 , 21
	lfpdx   27,%r18,%r11
	fxcpmadd  19 , 28 , 13 , 19
	la   %r13,	-1(%r13)
	fxpmul  20 , 30 , 12
	dcbt 	%r19,%r12
	fxpmul  23 , 30 , 15
	fxcxnpma  18 , 25 , 13 , 18
	fxcxnpma  21 , 25 , 16 , 21
	fxcxnpma  19 , 28 , 13 , 19
	dcbt 	%r17,%r11
	fxcxnpma  20 , 30 , 12 , 20
	lfpdx   25,%r15,%r14
	fxcxnpma  23 , 30 , 15 , 23
	fxcpmadd  18 , 26 , 14 , 18
	lfpsx   12,%r15,%r12
	fxcpmadd  21 , 26 , 17 , 21
	dcbt 	%r19,%r11
	fxcpmadd  19 , 29 , 14 , 19
	lfpsx   15,%r28,%r12
	fxcpmadd  22 , 28 , 16 , 22
	lfpdx   30,%r19,%r11
	fxcpmadd  20 , 31 , 13 , 20
	lfpdx   25,%r20,%r11
	fxcxnpma  18 , 26 , 14 , 18
	dcbt 	%r21,%r11
	fxcxnpma  21 , 26 , 17 , 21
	fxcxnpma  19 , 29 , 14 , 19
	fxcpmadd  23 , 31 , 16 , 23
	lfpdx   26,%r15,%r11
	fxcxnpma  22 , 28 , 16 , 22
	dcbt 	%r23,%r11
	fxcxnpma  20 , 31 , 13 , 20
	fpadd  0 , 12 , 18
	lfpdx   28,%r21,%r11
	fpadd  3 , 15 , 21
	lfpsx   13,%r27,%r12
	fxcxnpma  23 , 31 , 16 , 23
	dcbt 	%r25,%r11
	fxcpmadd  22 , 29 , 17 , 22
	fxcpmadd  20 , 24 , 14 , 20
	lfpsx   16,%r17,%r12
	fpsub  6 , 15 , 21
	lfpdx   31,%r22,%r11
	fpsub  9 , 12 , 18
	fxcpmadd  23 , 24 , 17 , 23
	fxcxnpma  22 , 29 , 17 , 22
	fxcxnpma  20 , 24 , 14 , 20
	fpadd  1 , 13 , 19
	lfpdx   29,%r16,%r11
	fpsub  10 , 13 , 19
	lfpsx   14,%r16,%r12
	fxcxnpma  23 , 24 , 17 , 23
	fpadd  4 , 16 , 22
	fpsub  7 , 16 , 22
	lfpsx   17,%r29,%r12
	lfpdx   24,%r17,%r11
	addi  %r12 , %r12 , 48
	fpadd  2 , 14 , 20
	dcbt 	%r17,%r12
	fpsub  11 , 14 , 20
	lfpsx   12,%r15,%r12
	lfpsx   15,%r28,%r12
	fpadd  5 , 17 , 23
	lfpsx   13,%r27,%r12
	fpsub  8 , 17 , 23
	dcbt 	%r19,%r12
	lfpsx   16,%r17,%r12
	lfpsx   14,%r16,%r12
	fxpmul  18 , 26 , 12
	lfpsx   17,%r29,%r12
	fxpmul  21 , 26 , 15
	fxpmul  19 , 29 , 12
	addi  %r12 , %r12 , 48
	fxpmul  22 , 29 , 15
	dcbt 	%r17,%r12
	fxpmul  20 , 24 , 12
	fxcxnpma  18 , 26 , 12 , 18
	fxcxnpma  21 , 26 , 15 , 21
	fxcxnpma  19 , 29 , 12 , 19
	dcbt 	%r19,%r12
	fxpmul  23 , 24 , 15
	lfpdx   26,%r23,%r11
	fxcxnpma  22 , 29 , 15 , 22
	fxcpmadd  18 , 27 , 13 , 18
	addi  %r11 , %r11 , 144
	fxcpmadd  21 , 27 , 16 , 21
	dcbt 	%r17,%r11
	fxcpmadd  19 , 30 , 13 , 19
	fxcxnpma  20 , 24 , 12 , 20
	fxcxnpma  23 , 24 , 15 , 23
	fxcxnpma  18 , 27 , 13 , 18
	dcbt 	%r19,%r11
	fxcxnpma  21 , 27 , 16 , 21
	lfpsx   12,%r15,%r12
	fxcxnpma  19 , 30 , 13 , 19
	lfpsx   15,%r28,%r12
	fxcpmadd  22 , 30 , 16 , 22
	lfpdx   27,%r15,%r14
	fxcpmadd  20 , 25 , 13 , 20
	dcbt 	%r21,%r11
	fxcpmadd  18 , 28 , 14 , 18
	fxcpmadd  21 , 28 , 17 , 21
	fxcpmadd  19 , 31 , 14 , 19
	fxcpmadd  23 , 25 , 16 , 23
	dcbt 	%r23,%r11
	fxcxnpma  22 , 30 , 16 , 22
	fxcxnpma  18 , 28 , 14 , 18
	fxcxnpma  21 , 28 , 17 , 21
	fxcxnpma  19 , 31 , 14 , 19
	dcbt 	%r25,%r11
	fxcxnpma  20 , 25 , 13 , 20
	fxcxnpma  23 , 25 , 16 , 23
	fxcpmadd  22 , 31 , 17 , 22
	lfpsx   13,%r27,%r12
	fpadd  0 , 0 , 12
	lfpsx   16,%r17,%r12
	fpadd  3 , 3 , 15
	fxcpmadd  20 , 26 , 14 , 20
	fxcpmadd  23 , 26 , 17 , 23
	fxcxnpma  22 , 31 , 17 , 22
	fpadd  1 , 1 , 13
	fpadd  4 , 4 , 16
	fxcxnpma  20 , 26 , 14 , 20
	fxcxnpma  23 , 26 , 17 , 23
	fpadd  0 , 0 , 18
	lfpsx   14,%r16,%r12
	fpadd  1 , 1 , 19
	lfpsx   17,%r29,%r12
	fpadd  3 , 3 , 21
	lfpdx   26,%r15,%r11
	addi  %r12 , %r12 , 48
	fpadd  4 , 4 , 22
	dcbt 	%r17,%r12
	fxcxnsma  28 , 27 , 6 , 15
	fpadd  2 , 2 , 14
	fpadd  5 , 5 , 17
	lfpsx   15,%r28,%r12
	fxcxnsma  29 , 27 , 7 , 16
	dcbt 	%r19,%r12
	fxcxnsma  30 , 27 , 8 , 17
	fpsub  6 , 28 , 21
	lfpsx   16,%r17,%r12
	fpadd  2 , 2 , 20
	lfpsx   17,%r29,%r12
	fpadd  5 , 5 , 23
	lfpdx   28,%r21,%r11
	fpsub  7 , 29 , 22
	fpsub  8 , 30 , 23
	fxcxnsma  31 , 27 , 9 , 18
	lfpdx   29,%r16,%r11
	fxcxnsma  24 , 27 , 10 , 19
	lfpdx   30,%r19,%r11
	fxcxnsma  25 , 27 , 11 , 20
	fxpmul  21 , 26 , 15
	lfpdx   27,%r18,%r11
	fpsub  9 , 31 , 12
	fpsub  10 , 24 , 13
	fpsub  11 , 25 , 14
	lfpsx   12,%r15,%r12
	fxcxnpma  21 , 26 , 15 , 21
	lfpsx   13,%r27,%r12
	fxpmul  22 , 29 , 15
	lfpsx   14,%r16,%r12
	lfpdx   31,%r22,%r11
	addi  %r12 , %r12 , 48
	lfpdx   24,%r17,%r11
	lfpdx   25,%r20,%r11
	fxpmul  18 , 26 , 12
	dcbt 	%r17,%r12
	fxpmul  19 , 29 , 12
	fxcpmadd  21 , 27 , 16 , 21
	fxcxnpma  22 , 29 , 15 , 22
	fxpmul  20 , 24 , 12
	dcbt 	%r19,%r12
	fxcxnpma  18 , 26 , 12 , 18
	fxcxnpma  19 , 29 , 12 , 19
	fxcxnpma  21 , 27 , 16 , 21
	lfpdx   26,%r23,%r11
	fxpmul  23 , 24 , 15
	fxcxnpma  20 , 24 , 12 , 20
	addi  %r11 , %r11 , 144
	fxcpmadd  18 , 27 , 13 , 18
	dcbt 	%r17,%r11
	fxcpmadd  19 , 30 , 13 , 19
	lfpsx   12,%r15,%r12
	fxcpmadd  21 , 28 , 17 , 21
	lfpdx   29,%r18,%r11
	fxcxnpma  23 , 24 , 15 , 23
	fxcpmadd  22 , 30 , 16 , 22
	dcbt 	%r19,%r11
	fxcxnpma  18 , 27 , 13 , 18
	lfpsx   15,%r28,%r12
	fxcxnpma  19 , 30 , 13 , 19
	lfpdx   24,%r19,%r11
	fxcxnpma  21 , 28 , 17 , 21
	lfpdx   27,%r15,%r14
	fxcpmadd  20 , 25 , 13 , 20
	dcbt 	%r21,%r11
	fxcpmadd  23 , 25 , 16 , 23
	fxcpmadd  18 , 28 , 14 , 18
	fxcpmadd  19 , 31 , 14 , 19
	fxcxnpma  22 , 30 , 16 , 22
	dcbt 	%r23,%r11
	fxcxnpma  20 , 25 , 13 , 20
	fxcxnpma  23 , 25 , 16 , 23
	lfpdx   30,%r21,%r11
	fxcxnpma  18 , 28 , 14 , 18
	lfpsx   13,%r27,%r12
	fxcxnpma  19 , 31 , 14 , 19
	dcbt 	%r25,%r11
	fxcpmadd  22 , 31 , 17 , 22
	lfpsx   16,%r17,%r12
	fxcpmadd  20 , 26 , 14 , 20
	lfpdx   28,%r15,%r11
	fxcpmadd  23 , 26 , 17 , 23
	lfpdx   25,%r22,%r11
	fpadd  0 , 0 , 12
	fpadd  1 , 1 , 13
	fxcxnpma  22 , 31 , 17 , 22
	fxcxnpma  20 , 26 , 14 , 20
	fxcxnpma  23 , 26 , 17 , 23
	lfpdx   31,%r16,%r11
	fpadd  3 , 3 , 15
	lfpsx   14,%r16,%r12
	fpadd  4 , 4 , 16
	lfpsx   17,%r29,%r12
	fpadd  0 , 0 , 18
	lfpdx   26,%r17,%r11
	addi  %r12 , %r12 , 48
	fpadd  1 , 1 , 19
	dcbt 	%r17,%r12
	fpadd  3 , 3 , 21
	fpadd  2 , 2 , 14
	fpadd  5 , 5 , 17
	fpadd  4 , 4 , 22
	dcbt 	%r19,%r12
	fxcxnsma  6 , 27 , 12 , 6
	fxcxnsma  7 , 27 , 13 , 7
	fpadd  2 , 2 , 20
	lfpsx   12,%r15,%r12
	fpadd  5 , 5 , 23
	lfpsx   13,%r27,%r12
	fxcxnsma  8 , 27 , 14 , 8
	fxcxnpma  9 , 27 , 15 , 9
	fxcxnpma  10 , 27 , 16 , 10
	lfpsx   14,%r16,%r12
	fxcxnpma  11 , 27 , 17 , 11
	lfpsx   15,%r28,%r12
	fxcxnpma  6 , 27 , 18 , 6
	lfpsx   16,%r17,%r12
	fxcxnpma  7 , 27 , 19 , 7
	lfpsx   17,%r29,%r12
	fxcxnpma  8 , 27 , 20 , 8
	fxcxnsma  9 , 27 , 21 , 9
	addi  %r12 , %r12 , 48
	fxcxnsma  10 , 27 , 22 , 10
	dcbt 	%r17,%r12
	fxcxnsma  11 , 27 , 23 , 11
	fxpmul  18 , 28 , 12
	fxpmul  21 , 28 , 15
	lfpdx   27,%r20,%r11
	fxpmul  19 , 31 , 12
	dcbt 	%r19,%r12
	fxpmul  22 , 31 , 15
	fxpmul  20 , 26 , 12
	fxcxnpma  18 , 28 , 12 , 18
	fxcxnpma  21 , 28 , 15 , 21
	fxcxnpma  19 , 31 , 12 , 19
	fxpmul  23 , 26 , 15
	lfpdx   28,%r23,%r11
	fxcxnpma  22 , 31 , 15 , 22
	fxcpmadd  18 , 29 , 13 , 18
	addi  %r11 , %r11 , 144
	fxcpmadd  21 , 29 , 16 , 21
	dcbt 	%r17,%r11
	fxcpmadd  19 , 24 , 13 , 19
	fxcxnpma  20 , 26 , 12 , 20
	fxcxnpma  23 , 26 , 15 , 23
	fxcxnpma  18 , 29 , 13 , 18
	dcbt 	%r19,%r11
	fxcxnpma  21 , 29 , 16 , 21
	lfpsx   12,%r15,%r12
	fxcxnpma  19 , 24 , 13 , 19
	lfpsx   15,%r28,%r12
	fxcpmadd  22 , 24 , 16 , 22
	lfpdx   29,%r15,%r14
	fxcpmadd  20 , 27 , 13 , 20
	dcbt 	%r21,%r11
	fxcpmadd  18 , 30 , 14 , 18
	fxcpmadd  21 , 30 , 17 , 21
	fxcpmadd  19 , 25 , 14 , 19
	fxcpmadd  23 , 27 , 16 , 23
	dcbt 	%r23,%r11
	fxcxnpma  22 , 24 , 16 , 22
	fxcxnpma  18 , 30 , 14 , 18
	fxcxnpma  21 , 30 , 17 , 21
	fxcxnpma  19 , 25 , 14 , 19
	dcbt 	%r25,%r11
	fxcxnpma  20 , 27 , 13 , 20
	fxcxnpma  23 , 27 , 16 , 23
	fxcpmadd  22 , 25 , 17 , 22
	lfpsx   13,%r27,%r12
	fpadd  0 , 0 , 12
	lfpsx   16,%r17,%r12
	fpadd  3 , 3 , 15
	fxcpmadd  20 , 28 , 14 , 20
	fxcpmadd  23 , 28 , 17 , 23
	fxcxnpma  22 , 25 , 17 , 22
	fpadd  1 , 1 , 13
	fpadd  4 , 4 , 16
	fxcxnpma  20 , 28 , 14 , 20
	fxcxnpma  23 , 28 , 17 , 23
	fpadd  0 , 0 , 18
	lfpsx   14,%r16,%r12
	fpadd  1 , 1 , 19
	lfpsx   17,%r29,%r12
	fpadd  3 , 3 , 21
	addi  %r12 , %r12 , 48
	fpadd  4 , 4 , 22
	dcbt 	%r17,%r12
	fpadd  6 , 6 , 12
	stfpdx  0,%r15,%r10
	fpadd  2 , 2 , 14
	stfpdx  1,%r16,%r10
	fpadd  5 , 5 , 17
	stfpdx  3,%r18,%r10
	fpadd  7 , 7 , 13
	dcbt 	%r19,%r12
	fpadd  8 , 8 , 14
	stfpdx  4,%r19,%r10
	fpadd  9 , 9 , 15
	fpadd  2 , 2 , 20
	fpadd  5 , 5 , 23
	fpadd  10 , 10 , 16
	fpadd  11 , 11 , 17
	fpsub  6 , 6 , 18
	fpsub  7 , 7 , 19
	stfpdx  2,%r17,%r10
	fpsub  8 , 8 , 20
	stfpdx  5,%r20,%r10
	fpsub  9 , 9 , 21
	fpsub  10 , 10 , 22
	fpsub  11 , 11 , 23
	stfpdx  6,%r21,%r10
	stfpdx  7,%r22,%r10
	stfpdx  8,%r23,%r10
	stfpdx  9,%r24,%r10
	stfpdx  10,%r25,%r10
	stfpdx  11,%r26,%r10
	addi  %r10 , %r10 , 192
	dcbt 	%r17,%r10
	dcbt 	%r19,%r10
	dcbt 	%r21,%r10
	dcbt 	%r23,%r10
	dcbt 	%r25,%r10
	bdnz    s_bgl_rec_su3_lab1
s_bgl_rec_su3_lab0:
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
