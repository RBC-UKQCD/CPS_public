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
.globl s_bgl_dec_hsu3_dag
s_bgl_dec_hsu3_dag:

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
	li   %r30,8
	stfd  27,	216(%r1)
	li   %r31,24
	stfd  28,	224(%r1)
	li   %r3,40
	stfd  29,	232(%r1)
	lwz   %r13,	0(%r13)
	stfd  30,	240(%r1)
	or.  %r13 , %r13 , %r13
	stfd  31,	248(%r1)
	bf gt,  s_bgl_dec_hsu3_dag_lab0
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
	b       s_bgl_dec_hsu3_dag_lab1
.align 4
s_bgl_dec_hsu3_dag_lab1:
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
	stfpsx  26,%r18,%r12
	fxcxnpma  29 , 12 , 6 , 3
	lfpdx   18,%r23,%r11
	fxcxnpma  30 , 12 , 7 , 4
	stfpsx  27,%r30,%r12
	fxcxnpma  31 , 12 , 8 , 5
	lfpdx   19,%r24,%r11
	fxcxnsma  23 , 12 , 6 , 3
	stfpsx  28,%r19,%r12
	fxcxnsma  24 , 12 , 7 , 4
	addi  %r16 , %r16 , 144
	fxcxnsma  25 , 12 , 8 , 5
	stfpsx  29,%r31,%r12
	fxpmul  26 , 13 , 20
	addi  %r10 , %r10 , 192
	fxpmul  27 , 16 , 20
	stfpsx  30,%r20,%r12
	fxpmul  29 , 13 , 23
	lfpdx   12,%r25,%r11
	fxpmul  30 , 16 , 23
	stfpsx  31,%r3,%r12
	fxpmul  28 , 19 , 20
	la   %r13,	-1(%r13)
	fxcxnsma  26 , 13 , 20 , 26
	lwz   %r12,	0(%r14)
	fxcxnsma  27 , 16 , 20 , 27
	fxcxnsma  29 , 13 , 23 , 29
	la   %r14,	4(%r14)
	fxpmul  31 , 19 , 23
	fxcxnsma  30 , 16 , 23 , 30
	lfpdx   13,%r26,%r11
	fxcpmadd  26 , 14 , 21 , 26
	dcbt 	%r20,%r16
	fxcpmadd  27 , 17 , 21 , 27
	addi  %r11 , %r11 , 144
	fxcpmadd  29 , 14 , 24 , 29
	lfpdx   16,%r19,%r11
	fxcxnsma  28 , 19 , 20 , 28
	fxcxnsma  31 , 19 , 23 , 31
	dcbt 	%r22,%r16
	fxcxnsma  26 , 14 , 21 , 26
	fxcxnsma  27 , 17 , 21 , 27
	lfpdx   19,%r22,%r11
	fxcxnsma  29 , 14 , 24 , 29
	fxcpmadd  30 , 17 , 24 , 30
	dcbt 	%r24,%r16
	fxcpmadd  28 , 12 , 21 , 28
	fxcpmadd  26 , 15 , 22 , 26
	lfpdx   14,%r18,%r15
	fxcpmadd  27 , 18 , 22 , 27
	fxcpmadd  29 , 15 , 25 , 29
	dcbt 	%r26,%r16
	fxcpmadd  31 , 12 , 24 , 31
	fxcxnsma  30 , 17 , 24 , 30
	lfpdx   14,%r25,%r11
	fxcxnsma  26 , 15 , 22 , 26
	addi  %r16 , %r16 , 144
	fxcxnsma  27 , 18 , 22 , 27
	lfpdx   17,%r20,%r11
	fxcxnsma  29 , 15 , 25 , 29
	dcbt 	%r20,%r16
	fxcxnsma  28 , 12 , 21 , 28
	fxcxnsma  31 , 12 , 24 , 31
	lfpdx   15,%r18,%r11
	fxcpmadd  30 , 18 , 25 , 30
	stfpsx  26,%r18,%r12
	fpadd  20 , 0 , 9
	lfpdx   12,%r23,%r11
	fpsub  26 , 0 , 9
	stfpsx  27,%r30,%r12
	fxcpmadd  28 , 13 , 22 , 28
	fxcpmadd  31 , 13 , 25 , 31
	stfpsx  29,%r31,%r12
	fxcxnsma  30 , 18 , 25 , 30
	dcbt 	%r22,%r16
	fpsub  27 , 1 , 10
	fpadd  21 , 1 , 10
	lfpdx   18,%r21,%r11
	fxcxnsma  28 , 13 , 22 , 28
	fxcxnsma  31 , 13 , 25 , 31
	fpadd  22 , 2 , 11
	stfpsx  30,%r20,%r12
	fpadd  29 , 3 , 6
	lfpdx   13,%r24,%r11
	fpadd  30 , 4 , 7
	fpsub  23 , 3 , 6
	stfpsx  28,%r19,%r12
	fpsub  24 , 4 , 7
	fpsub  28 , 2 , 11
	stfpsx  31,%r3,%r12
	fpsub  25 , 5 , 8
	fpadd  31 , 5 , 8
	lwz   %r12,	0(%r14)
	stfpsx  26,%r18,%r12
	la   %r14,	4(%r14)
	fxpmul  26 , 15 , 20
	stfpsx  27,%r30,%r12
	fxpmul  27 , 18 , 20
	stfpsx  28,%r19,%r12
	fxpmul  28 , 13 , 20
	stfpsx  29,%r31,%r12
	fxcxnsma  26 , 15 , 20 , 26
	fxpmul  29 , 15 , 23
	stfpsx  30,%r20,%r12
	fxcxnsma  27 , 18 , 20 , 27
	fxpmul  30 , 18 , 23
	stfpsx  31,%r3,%r12
	fxcxnsma  28 , 13 , 20 , 28
	fxcpmadd  26 , 16 , 21 , 26
	lwz   %r12,	0(%r14)
	fxcxnsma  29 , 15 , 23 , 29
	fxcpmadd  27 , 19 , 21 , 27
	la   %r14,	4(%r14)
	fxpmul  31 , 13 , 23
	lfpdx   15,%r26,%r11
	fxcxnsma  30 , 18 , 23 , 30
	dcbt 	%r24,%r16
	fxcxnsma  26 , 16 , 21 , 26
	addi  %r11 , %r11 , 144
	fxcpmadd  29 , 16 , 24 , 29
	lfpdx   18,%r19,%r11
	fxcxnsma  27 , 19 , 21 , 27
	fxcxnsma  31 , 13 , 23 , 31
	dcbt 	%r26,%r16
	fxcpmadd  30 , 19 , 24 , 30
	fxcpmadd  26 , 17 , 22 , 26
	lfpdx   13,%r22,%r11
	fxcxnsma  29 , 16 , 24 , 29
	addi  %r16 , %r16 , 144
	fxcpmadd  27 , 12 , 22 , 27
	fxcpmadd  28 , 14 , 21 , 28
	lfpdx   16,%r18,%r15
	fxcpmadd  31 , 14 , 24 , 31
	dcbt 	%r20,%r16
	fxcxnsma  26 , 17 , 22 , 26
	fxcpmadd  29 , 17 , 25 , 29
	fxcxnsma  27 , 12 , 22 , 27
	fxcxnsma  30 , 19 , 24 , 30
	dcbt 	%r22,%r16
	fxcxnsma  28 , 14 , 21 , 28
	fxcxnsma  31 , 14 , 24 , 31
	stfpsx  26,%r18,%r12
	fxcxnsma  29 , 17 , 25 , 29
	lfpdx   19,%r20,%r11
	fxcxnpma  26 , 16 , 6 , 0
	stfpsx  27,%r30,%r12
	fxcpmadd  30 , 12 , 25 , 30
	lfpdx   17,%r18,%r11
	fxcpmadd  28 , 15 , 22 , 28
	lfpdx   14,%r23,%r11
	fxcpmadd  31 , 15 , 25 , 31
	fxcxnpma  27 , 16 , 7 , 1
	stfpsx  29,%r31,%r12
	fxcxnsma  20 , 16 , 6 , 0
	dcbt 	%r24,%r16
	fxcxnsma  30 , 12 , 25 , 30
	fxcxnsma  28 , 15 , 22 , 28
	fxcxnsma  31 , 15 , 25 , 31
	lfpdx   12,%r21,%r11
	fxcxnsma  21 , 16 , 7 , 1
	fxcxnsma  22 , 16 , 8 , 2
	lfpdx   15,%r24,%r11
	fxcxnsma  29 , 16 , 9 , 3
	stfpsx  30,%r20,%r12
	fxcxnpma  23 , 16 , 9 , 3
	stfpsx  28,%r19,%r12
	fxcxnsma  30 , 16 , 10 , 4
	fxcxnpma  28 , 16 , 8 , 2
	stfpsx  31,%r3,%r12
	fxcxnpma  24 , 16 , 10 , 4
	fxcxnsma  31 , 16 , 11 , 5
	lwz   %r12,	0(%r14)
	fxcxnpma  25 , 16 , 11 , 5
	stfpsx  26,%r18,%r12
	la   %r14,	4(%r14)
	fxpmul  26 , 17 , 20
	stfpsx  27,%r30,%r12
	lfpdx   16,%r25,%r11
	fxpmul  27 , 12 , 20
	stfpsx  28,%r19,%r12
	fxpmul  28 , 15 , 20
	stfpsx  29,%r31,%r12
	fxcxnsma  26 , 17 , 20 , 26
	fxpmul  29 , 17 , 23
	stfpsx  30,%r20,%r12
	fxcxnsma  27 , 12 , 20 , 27
	fxpmul  30 , 12 , 23
	stfpsx  31,%r3,%r12
	fxcxnsma  28 , 15 , 20 , 28
	fxcpmadd  26 , 18 , 21 , 26
	lwz   %r12,	0(%r14)
	fxcxnsma  29 , 17 , 23 , 29
	fxcpmadd  27 , 13 , 21 , 27
	la   %r14,	4(%r14)
	fxpmul  31 , 15 , 23
	lfpdx   17,%r26,%r11
	fxcxnsma  30 , 12 , 23 , 30
	dcbt 	%r26,%r16
	fxcxnsma  26 , 18 , 21 , 26
	addi  %r11 , %r11 , 144
	fxcpmadd  29 , 18 , 24 , 29
	lfpdx   12,%r19,%r11
	fxcxnsma  27 , 13 , 21 , 27
	addi  %r16 , %r16 , 144
	fxcxnsma  31 , 15 , 23 , 31
	dcbt 	%r20,%r16
	fxcpmadd  30 , 13 , 24 , 30
	fxcpmadd  26 , 19 , 22 , 26
	lfpdx   15,%r22,%r11
	fxcxnsma  29 , 18 , 24 , 29
	fxcpmadd  27 , 14 , 22 , 27
	fxcpmadd  28 , 16 , 21 , 28
	lfpdx   18,%r18,%r15
	fxcpmadd  31 , 16 , 24 , 31
	fxcxnsma  26 , 19 , 22 , 26
	lfpdx   18,%r25,%r11
	fxcpmadd  29 , 19 , 25 , 29
	dcbt 	%r22,%r16
	fxcxnsma  27 , 14 , 22 , 27
	fxcxnsma  30 , 13 , 24 , 30
	fxcxnsma  28 , 16 , 21 , 28
	fxcxnsma  31 , 16 , 24 , 31
	stfpsx  26,%r18,%r12
	fxcxnsma  29 , 19 , 25 , 29
	lfpdx   13,%r20,%r11
	fpadd  26 , 0 , 6
	stfpsx  27,%r30,%r12
	fxcpmadd  30 , 14 , 25 , 30
	lfpdx   19,%r18,%r11
	fxcpmadd  28 , 17 , 22 , 28
	lfpdx   16,%r23,%r11
	fxcpmadd  31 , 17 , 25 , 31
	fpadd  27 , 1 , 7
	stfpsx  29,%r31,%r12
	fpsub  20 , 0 , 6
	dcbt 	%r24,%r16
	fxcxnsma  30 , 14 , 25 , 30
	fxcxnsma  28 , 17 , 22 , 28
	lfpdx   0,%r18,%r10
	fxcxnsma  31 , 17 , 25 , 31
	lfpdx   14,%r21,%r11
	fpsub  21 , 1 , 7
	lfpdx   6,%r24,%r10
	fpsub  22 , 2 , 8
	lfpdx   17,%r24,%r11
	fpadd  29 , 3 , 9
	stfpsx  30,%r20,%r12
	fpsub  23 , 3 , 9
	stfpsx  28,%r19,%r12
	fpadd  30 , 4 , 10
	lfpdx   1,%r19,%r10
	fpadd  28 , 2 , 8
	stfpsx  31,%r3,%r12
	fpsub  24 , 4 , 10
	lfpdx   3,%r21,%r10
	fpadd  31 , 5 , 11
	lwz   %r12,	0(%r14)
	fpsub  25 , 5 , 11
	lfpdx   9,%r27,%r10
	stfpsx  26,%r18,%r12
	la   %r14,	4(%r14)
	fxpmul  26 , 19 , 20
	stfpsx  27,%r30,%r12
	lfpdx   4,%r22,%r10
	fxpmul  27 , 14 , 20
	stfpsx  28,%r19,%r12
	lfpdx   7,%r25,%r10
	fxpmul  28 , 17 , 20
	stfpsx  29,%r31,%r12
	fxcxnsma  26 , 19 , 20 , 26
	lfpdx   10,%r28,%r10
	fxpmul  29 , 19 , 23
	stfpsx  30,%r20,%r12
	fxcxnsma  27 , 14 , 20 , 27
	lfpdx   2,%r20,%r10
	fxpmul  30 , 14 , 23
	stfpsx  31,%r3,%r12
	fxcxnsma  28 , 17 , 20 , 28
	lfpdx   5,%r23,%r10
	fxcpmadd  26 , 12 , 21 , 26
	lwz   %r12,	0(%r14)
	fxcxnsma  29 , 19 , 23 , 29
	lfpdx   8,%r26,%r10
	fxcpmadd  27 , 15 , 21 , 27
	la   %r14,	4(%r14)
	fxpmul  31 , 17 , 23
	lfpdx   19,%r26,%r11
	fxcxnsma  30 , 14 , 23 , 30
	dcbt 	%r26,%r16
	fxcxnsma  26 , 12 , 21 , 26
	addi  %r11 , %r11 , 144
	fxcpmadd  29 , 12 , 24 , 29
	lfpdx   11,%r29,%r10
	fxcxnsma  27 , 15 , 21 , 27
	fxcxnsma  31 , 17 , 23 , 31
	fxcpmadd  30 , 15 , 24 , 30
	fxcpmadd  26 , 13 , 22 , 26
	fxcxnsma  29 , 12 , 24 , 29
	fxcpmadd  27 , 16 , 22 , 27
	fxcpmadd  28 , 18 , 21 , 28
	lfpdx   12,%r18,%r15
	fxcpmadd  31 , 18 , 24 , 31
	fxcxnsma  26 , 13 , 22 , 26
	fxcpmadd  29 , 13 , 25 , 29
	fxcxnsma  27 , 16 , 22 , 27
	fxcxnsma  30 , 15 , 24 , 30
	fxcxnsma  28 , 18 , 21 , 28
	fxcxnsma  31 , 18 , 24 , 31
	stfpsx  26,%r18,%r12
	fxcxnsma  29 , 13 , 25 , 29
	stfpsx  27,%r30,%r12
	fxcpmadd  30 , 16 , 25 , 30
	fxcpmadd  28 , 19 , 22 , 28
	fxcpmadd  31 , 19 , 25 , 31
	stfpsx  29,%r31,%r12
	fxcxnsma  30 , 16 , 25 , 30
	fxcxnsma  28 , 19 , 22 , 28
	fxcxnsma  31 , 19 , 25 , 31
	stfpsx  30,%r20,%r12
	stfpsx  28,%r19,%r12
	stfpsx  31,%r3,%r12
	lwz   %r12,	0(%r14)
	bdnz    s_bgl_dec_hsu3_dag_lab1
s_bgl_dec_hsu3_dag_lab0:
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
