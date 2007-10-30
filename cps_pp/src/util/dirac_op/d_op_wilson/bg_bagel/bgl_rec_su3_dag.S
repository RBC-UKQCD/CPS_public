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
.globl bgl_rec_su3_dag
bgl_rec_su3_dag:

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
	lwz   %r13,	0(%r13)
	or.  %r13 , %r13 , %r13
	bf gt,  bgl_rec_su3_dag_lab0
	mtctr 	%r13
	b       bgl_rec_su3_dag_lab1
.align 4
bgl_rec_su3_dag_lab1:
	lfpdx   24,%r15,%r11
	lfpdx   18,%r15,%r12
	lfpdx   21,%r18,%r12
	lfpdx   27,%r16,%r11
	lfpdx   25,%r18,%r11
	lfpdx   19,%r16,%r12
	lfpdx   22,%r19,%r12
	fxpmul  12 , 24 , 18
	lfpdx   28,%r19,%r11
	fxpmul  15 , 24 , 21
	lfpdx   26,%r21,%r11
	fxpmul  13 , 27 , 18
	lfpdx   20,%r17,%r12
	fxpmul  16 , 27 , 21
	lfpdx   23,%r20,%r12
	lfpdx   29,%r22,%r11
	fxcxnpma  12 , 24 , 18 , 12
	lfpdx   30,%r17,%r11
	fxcxnpma  15 , 24 , 21 , 15
	lfpdx   31,%r20,%r11
	fxcxnpma  13 , 27 , 18 , 13
	addi  %r12 , %r12 , 96
	fxcxnpma  16 , 27 , 21 , 16
	lfpdx   24,%r23,%r11
	dcbt 	%r17,%r12
	fxcpmadd  12 , 25 , 19 , 12
	addi  %r11 , %r11 , 144
	fxcpmadd  15 , 25 , 22 , 15
	lfpdx   27,%r18,%r11
	fxcpmadd  13 , 28 , 19 , 13
	la   %r13,	-1(%r13)
	fxpmul  14 , 30 , 18
	dcbt 	%r19,%r12
	fxpmul  17 , 30 , 21
	fxcxnpma  12 , 25 , 19 , 12
	fxcxnpma  15 , 25 , 22 , 15
	fxcxnpma  13 , 28 , 19 , 13
	dcbt 	%r21,%r12
	fxcxnpma  14 , 30 , 18 , 14
	lfpdx   25,%r15,%r14
	fxcxnpma  17 , 30 , 21 , 17
	fxcpmadd  12 , 26 , 20 , 12
	lfpdx   18,%r15,%r12
	fxcpmadd  15 , 26 , 23 , 15
	dcbt 	%r17,%r11
	fxcpmadd  13 , 29 , 20 , 13
	lfpdx   21,%r18,%r12
	fxcpmadd  16 , 28 , 22 , 16
	lfpdx   30,%r19,%r11
	fxcpmadd  14 , 31 , 19 , 14
	lfpdx   25,%r20,%r11
	fxcxnpma  12 , 26 , 20 , 12
	dcbt 	%r19,%r11
	fxcxnpma  15 , 26 , 23 , 15
	fxcxnpma  13 , 29 , 20 , 13
	fxcpmadd  17 , 31 , 22 , 17
	lfpdx   26,%r15,%r11
	fxcxnpma  16 , 28 , 22 , 16
	dcbt 	%r21,%r11
	fxcxnpma  14 , 31 , 19 , 14
	fpadd  0 , 12 , 18
	lfpdx   28,%r21,%r11
	fpadd  3 , 15 , 21
	lfpdx   19,%r16,%r12
	fxcxnpma  17 , 31 , 22 , 17
	dcbt 	%r23,%r11
	fxcpmadd  16 , 29 , 23 , 16
	fxcpmadd  14 , 24 , 20 , 14
	lfpdx   22,%r19,%r12
	fpsub  6 , 15 , 21
	lfpdx   31,%r22,%r11
	fpsub  9 , 12 , 18
	dcbt 	%r25,%r11
	fxcpmadd  17 , 24 , 23 , 17
	fxcxnpma  16 , 29 , 23 , 16
	fxcxnpma  14 , 24 , 20 , 14
	fpadd  1 , 13 , 19
	lfpdx   29,%r16,%r11
	fpsub  10 , 13 , 19
	lfpdx   20,%r17,%r12
	fxcxnpma  17 , 24 , 23 , 17
	fpadd  4 , 16 , 22
	fpsub  7 , 16 , 22
	lfpdx   23,%r20,%r12
	lfpdx   24,%r17,%r11
	addi  %r12 , %r12 , 96
	fpadd  2 , 14 , 20
	dcbt 	%r17,%r12
	fpsub  11 , 14 , 20
	lfpdx   18,%r15,%r12
	lfpdx   21,%r18,%r12
	fpadd  5 , 17 , 23
	lfpdx   19,%r16,%r12
	fpsub  8 , 17 , 23
	dcbt 	%r19,%r12
	lfpdx   22,%r19,%r12
	lfpdx   20,%r17,%r12
	fxpmul  12 , 26 , 18
	lfpdx   23,%r20,%r12
	fxpmul  15 , 26 , 21
	dcbt 	%r21,%r12
	fxpmul  13 , 29 , 18
	fxpmul  16 , 29 , 21
	addi  %r12 , %r12 , 96
	fxpmul  14 , 24 , 18
	fxcxnpma  12 , 26 , 18 , 12
	dcbt 	%r17,%r12
	fxcxnpma  15 , 26 , 21 , 15
	fxcxnpma  13 , 29 , 18 , 13
	fxpmul  17 , 24 , 21
	lfpdx   26,%r23,%r11
	fxcxnpma  16 , 29 , 21 , 16
	dcbt 	%r19,%r12
	fxcpmadd  12 , 27 , 19 , 12
	addi  %r11 , %r11 , 144
	fxcpmadd  15 , 27 , 22 , 15
	fxcpmadd  13 , 30 , 19 , 13
	fxcxnpma  14 , 24 , 18 , 14
	dcbt 	%r21,%r12
	fxcxnpma  17 , 24 , 21 , 17
	fxcxnpma  12 , 27 , 19 , 12
	lfpdx   18,%r15,%r12
	fxcxnpma  15 , 27 , 22 , 15
	lfpdx   21,%r18,%r12
	fxcxnpma  13 , 30 , 19 , 13
	dcbt 	%r17,%r11
	fxcpmadd  16 , 30 , 22 , 16
	lfpdx   27,%r15,%r14
	fxcpmadd  14 , 25 , 19 , 14
	fxcpmadd  12 , 28 , 20 , 12
	fxcpmadd  15 , 28 , 23 , 15
	dcbt 	%r19,%r11
	fxcpmadd  13 , 31 , 20 , 13
	fxcpmadd  17 , 25 , 22 , 17
	fxcxnpma  16 , 30 , 22 , 16
	fxcxnpma  12 , 28 , 20 , 12
	dcbt 	%r21,%r11
	fxcxnpma  15 , 28 , 23 , 15
	fxcxnpma  13 , 31 , 20 , 13
	fxcxnpma  14 , 25 , 19 , 14
	fxcxnpma  17 , 25 , 22 , 17
	dcbt 	%r23,%r11
	fxcpmadd  16 , 31 , 23 , 16
	lfpdx   19,%r16,%r12
	fpadd  0 , 0 , 12
	lfpdx   22,%r19,%r12
	fpadd  1 , 1 , 13
	fxcpmadd  14 , 26 , 20 , 14
	dcbt 	%r25,%r11
	fxcpmadd  17 , 26 , 23 , 17
	fxcxnpma  16 , 31 , 23 , 16
	fpadd  3 , 3 , 15
	fpadd  0 , 0 , 18
	fxcxnpma  14 , 26 , 20 , 14
	fxcxnpma  17 , 26 , 23 , 17
	fpadd  4 , 4 , 16
	lfpdx   20,%r17,%r12
	fpadd  1 , 1 , 19
	lfpdx   23,%r20,%r12
	fpadd  3 , 3 , 21
	lfpdx   26,%r15,%r11
	addi  %r12 , %r12 , 96
	fpadd  2 , 2 , 14
	dcbt 	%r17,%r12
	fpadd  5 , 5 , 17
	fpadd  4 , 4 , 22
	fxcxnsma  28 , 27 , 6 , 15
	fxcxnsma  29 , 27 , 7 , 16
	dcbt 	%r19,%r12
	fpadd  2 , 2 , 20
	fpadd  5 , 5 , 23
	fxcxnsma  30 , 27 , 8 , 17
	fpsub  6 , 28 , 21
	dcbt 	%r21,%r12
	fpsub  7 , 29 , 22
	fxcxnsma  31 , 27 , 9 , 18
	lfpdx   21,%r18,%r12
	fxcxnsma  24 , 27 , 10 , 19
	lfpdx   29,%r16,%r11
	fpsub  8 , 30 , 23
	lfpdx   18,%r15,%r12
	fxcxnsma  25 , 27 , 11 , 20
	lfpdx   19,%r16,%r12
	lfpdx   22,%r19,%r12
	fpsub  9 , 31 , 12
	lfpdx   27,%r18,%r11
	fpsub  10 , 24 , 13
	lfpdx   30,%r19,%r11
	fxpmul  15 , 26 , 21
	lfpdx   28,%r21,%r11
	fpsub  11 , 25 , 14
	lfpdx   20,%r17,%r12
	fxpmul  12 , 26 , 18
	lfpdx   23,%r20,%r12
	fxpmul  13 , 29 , 18
	lfpdx   31,%r22,%r11
	fxpmul  16 , 29 , 21
	lfpdx   24,%r17,%r11
	fxcxnpma  15 , 26 , 21 , 15
	lfpdx   25,%r20,%r11
	addi  %r12 , %r12 , 96
	fxcxnpma  12 , 26 , 18 , 12
	dcbt 	%r17,%r12
	fxcxnpma  13 , 29 , 18 , 13
	fxcxnpma  16 , 29 , 21 , 16
	lfpdx   26,%r23,%r11
	fxcpmadd  15 , 27 , 22 , 15
	fxpmul  14 , 24 , 18
	dcbt 	%r19,%r12
	fxcpmadd  12 , 27 , 19 , 12
	addi  %r11 , %r11 , 144
	fxcpmadd  13 , 30 , 19 , 13
	lfpdx   29,%r18,%r11
	fxpmul  17 , 24 , 21
	fxcxnpma  15 , 27 , 22 , 15
	dcbt 	%r21,%r12
	fxcxnpma  14 , 24 , 18 , 14
	fxcxnpma  12 , 27 , 19 , 12
	fxcxnpma  13 , 30 , 19 , 13
	lfpdx   18,%r15,%r12
	fxcxnpma  17 , 24 , 21 , 17
	dcbt 	%r17,%r11
	fxcpmadd  15 , 28 , 23 , 15
	lfpdx   27,%r15,%r14
	fxcpmadd  16 , 30 , 22 , 16
	lfpdx   21,%r18,%r12
	fxcpmadd  12 , 28 , 20 , 12
	lfpdx   24,%r19,%r11
	fxcpmadd  13 , 31 , 20 , 13
	dcbt 	%r19,%r11
	fxcpmadd  14 , 25 , 19 , 14
	fxcxnpma  15 , 28 , 23 , 15
	fxcpmadd  17 , 25 , 22 , 17
	fxcxnpma  12 , 28 , 20 , 12
	dcbt 	%r21,%r11
	fxcxnpma  13 , 31 , 20 , 13
	fxcxnpma  16 , 30 , 22 , 16
	lfpdx   28,%r15,%r11
	fxcxnpma  14 , 25 , 19 , 14
	fxcxnpma  17 , 25 , 22 , 17
	dcbt 	%r23,%r11
	fpadd  0 , 0 , 12
	lfpdx   19,%r16,%r12
	fpadd  1 , 1 , 13
	lfpdx   22,%r19,%r12
	fxcpmadd  16 , 31 , 23 , 16
	lfpdx   30,%r21,%r11
	fxcpmadd  14 , 26 , 20 , 14
	dcbt 	%r25,%r11
	fxcpmadd  17 , 26 , 23 , 17
	lfpdx   25,%r22,%r11
	fpadd  3 , 3 , 15
	fpadd  0 , 0 , 18
	fxcxnpma  16 , 31 , 23 , 16
	fxcxnpma  14 , 26 , 20 , 14
	fxcxnpma  17 , 26 , 23 , 17
	lfpdx   31,%r16,%r11
	fpadd  1 , 1 , 19
	lfpdx   20,%r17,%r12
	fpadd  3 , 3 , 21
	lfpdx   23,%r20,%r12
	fpadd  4 , 4 , 16
	lfpdx   26,%r17,%r11
	addi  %r12 , %r12 , 96
	fpadd  2 , 2 , 14
	dcbt 	%r17,%r12
	fpadd  5 , 5 , 17
	fxcxnsma  6 , 27 , 12 , 6
	fpadd  4 , 4 , 22
	fxcxnsma  7 , 27 , 13 , 7
	dcbt 	%r19,%r12
	fpadd  2 , 2 , 20
	fpadd  5 , 5 , 23
	fxcxnsma  8 , 27 , 14 , 8
	fxcxnpma  9 , 27 , 15 , 9
	dcbt 	%r21,%r12
	fxcxnpma  10 , 27 , 16 , 10
	fxcxnpma  11 , 27 , 17 , 11
	fxcxnpma  6 , 27 , 18 , 6
	fxcxnpma  7 , 27 , 19 , 7
	fxcxnpma  8 , 27 , 20 , 8
	lfpdx   18,%r15,%r12
	fxcxnsma  9 , 27 , 21 , 9
	lfpdx   19,%r16,%r12
	fxcxnsma  10 , 27 , 22 , 10
	lfpdx   20,%r17,%r12
	fxcxnsma  11 , 27 , 23 , 11
	lfpdx   21,%r18,%r12
	lfpdx   22,%r19,%r12
	lfpdx   23,%r20,%r12
	fxpmul  12 , 28 , 18
	lfpdx   27,%r20,%r11
	fxpmul  13 , 31 , 18
	addi  %r12 , %r12 , 96
	fxpmul  14 , 26 , 18
	dcbt 	%r17,%r12
	fxpmul  15 , 28 , 21
	fxpmul  16 , 31 , 21
	fxcxnpma  12 , 28 , 18 , 12
	fxcxnpma  13 , 31 , 18 , 13
	dcbt 	%r19,%r12
	fxpmul  17 , 26 , 21
	fxcxnpma  15 , 28 , 21 , 15
	fxcxnpma  16 , 31 , 21 , 16
	fxcpmadd  12 , 29 , 19 , 12
	lfpdx   28,%r23,%r11
	fxcpmadd  13 , 24 , 19 , 13
	dcbt 	%r21,%r12
	fxcxnpma  14 , 26 , 18 , 14
	addi  %r11 , %r11 , 144
	fxcpmadd  15 , 29 , 22 , 15
	fxcxnpma  17 , 26 , 21 , 17
	lfpdx   18,%r15,%r12
	fxcxnpma  12 , 29 , 19 , 12
	dcbt 	%r17,%r11
	fxcxnpma  13 , 24 , 19 , 13
	lfpdx   21,%r18,%r12
	fxcpmadd  16 , 24 , 22 , 16
	fxcxnpma  15 , 29 , 22 , 15
	fxcpmadd  14 , 27 , 19 , 14
	dcbt 	%r19,%r11
	fxcpmadd  12 , 30 , 20 , 12
	lfpdx   29,%r15,%r14
	fxcpmadd  13 , 25 , 20 , 13
	fxcpmadd  17 , 27 , 22 , 17
	fxcpmadd  15 , 30 , 23 , 15
	dcbt 	%r21,%r11
	fxcxnpma  16 , 24 , 22 , 16
	fxcxnpma  12 , 30 , 20 , 12
	fxcxnpma  13 , 25 , 20 , 13
	fxcxnpma  14 , 27 , 19 , 14
	dcbt 	%r23,%r11
	fxcxnpma  15 , 30 , 23 , 15
	fxcxnpma  17 , 27 , 22 , 17
	lfpdx   19,%r16,%r12
	fxcpmadd  16 , 25 , 23 , 16
	fpadd  0 , 0 , 12
	dcbt 	%r25,%r11
	fxcpmadd  14 , 28 , 20 , 14
	lfpdx   22,%r19,%r12
	fpadd  1 , 1 , 13
	fxcpmadd  17 , 28 , 23 , 17
	fxcxnpma  16 , 25 , 23 , 16
	fpadd  3 , 3 , 15
	fxcxnpma  14 , 28 , 20 , 14
	fpadd  0 , 0 , 18
	fxcxnpma  17 , 28 , 23 , 17
	lfpdx   20,%r17,%r12
	fpadd  4 , 4 , 16
	fpadd  1 , 1 , 19
	lfpdx   23,%r20,%r12
	fpadd  2 , 2 , 14
	addi  %r12 , %r12 , 96
	stfpdx  0,%r15,%r10
	fpadd  5 , 5 , 17
	dcbt 	%r17,%r12
	fpadd  3 , 3 , 21
	fpadd  4 , 4 , 22
	stfpdx  1,%r16,%r10
	fpadd  2 , 2 , 20
	fpadd  6 , 6 , 12
	dcbt 	%r19,%r12
	fpadd  5 , 5 , 23
	fpadd  7 , 7 , 13
	stfpdx  3,%r18,%r10
	fpadd  8 , 8 , 14
	stfpdx  4,%r19,%r10
	fpadd  9 , 9 , 15
	dcbt 	%r21,%r12
	fpadd  10 , 10 , 16
	stfpdx  2,%r17,%r10
	fpadd  11 , 11 , 17
	stfpdx  5,%r20,%r10
	fpsub  6 , 6 , 18
	fpsub  7 , 7 , 19
	fpsub  8 , 8 , 20
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
	bdnz    bgl_rec_su3_dag_lab1
bgl_rec_su3_dag_lab0:
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
