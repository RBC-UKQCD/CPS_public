	.section	".text"
	.align 2
	.globl pt_cmv_dag_pad
	.type	pt_cmv_dag_pad,@function
pt_cmv_dag_pad:
	la   %r1,	-864(%r1)
	stw	%r14,	368(%r1)
	stw	%r15,	376(%r1)
	stw	%r16,	384(%r1)
	stw	%r17,	392(%r1)
	stw	%r18,	400(%r1)
	stw	%r19,	408(%r1)
	stw	%r20,	416(%r1)
	stw	%r21,	424(%r1)
	stw	%r22,	432(%r1)
	stw	%r23,	440(%r1)
	stw	%r24,	448(%r1)
	stw	%r25,	456(%r1)
	stw	%r26,	464(%r1)
	stw	%r27,	472(%r1)
	stw	%r28,	480(%r1)
	stw	%r29,	488(%r1)
	stw	%r30,	496(%r1)
	stw	%r31,	504(%r1)
	stfd	14,	112(%r1)
	stfd	15,	120(%r1)
	stfd	16,	128(%r1)
	stfd	17,	136(%r1)
	stfd	18,	144(%r1)
	stfd	19,	152(%r1)
	stfd	20,	160(%r1)
	stfd	21,	168(%r1)
	stfd	22,	176(%r1)
	stfd	23,	184(%r1)
	stfd	24,	192(%r1)
	stfd	25,	200(%r1)
	stfd	26,	208(%r1)
	stfd	27,	216(%r1)
	stfd	28,	224(%r1)
	stfd	29,	232(%r1)
	stfd	30,	240(%r1)
	stfd	31,	248(%r1)
	li   8,	0
	li   9,	16
	li   10,	32
	li   11,	48
	li   12,	64
	li   13,	80
	li   14,	96
	li   15,	112
	li   16,	128
	li   17,	160
	dcbt	 %r8,	 %r4
	dcbt	 %r10,	 %r4
	lwz	%r19, 0(%r4)
	lwz	%r23, 4(%r4)
	lwz	%r27, 8(%r4)
	dcbt	 %r12,	 %r4
	lwz	%r20, 12(%r4)
	lwz	%r24, 16(%r4)
	lwz	%r28, 20(%r4)
	add    %r27, %r27, %r5
	add    %r19, %r19, %r6
	add    %r23, %r23, %r7
	add    %r20, %r20, %r6
	add    %r24, %r24, %r7
	add    %r28, %r28, %r5
	dcbt	 %r14,	 %r4
	or.	 3,	3,	3
	mr	 %r26,	 %r24
	bf gt,	pt_cmv_dag_pad_lab0
	lfpdx	9, %r27, %r8
	lfpdx	10, %r27, %r9
	lfpdx	11, %r27, %r10
	add    %r18, %r5, %r8
	lfpdx	0, %r19, %r8
	mr	 %r26,	 %r24
	subic.	 3,	3,	1
	bf gt,	pt_cmv_dag_pad_lab2
	b  pt_cmv_dag_pad_lab1
pt_cmv_dag_pad_lab1:
	lwz	%r21, 24(%r4)
	fxpmul   3,  9,  0
	lwz	%r25, 28(%r4)
	dcbt	 %r16,	 %r30
	lwz	%r29, 32(%r4)
	fxpmul   4, 10,  0
	add    %r21, %r21, %r6
	lfpdx	12, %r27, %r11
	add    %r25, %r25, %r7
	fxpmul   5, 11,  0
	add    %r29, %r29, %r5

	fxcxnsma   3,  9,  0,  3
	stfpdx	6, %r24, %r8
	fxcxnsma   4, 10,  0,  4

	fxcxnsma   5, 11,  0,  5
	lfpdx	13, %r27, %r12
	stfpdx	7, %r24, %r9
	lfpdx	14, %r27, %r13
	stfpdx	8, %r24, %r10
	lfpdx	1, %r19, %r9
	stfpdx	8, %r24, %r10
	lfpdx	15, %r27, %r14
	fxcpmadd   3, 12,  1,  3
	dcbt	 %r8,	 %r29
	fxcpmadd   4, 13,  1,  4
	lfpdx	16, %r27, %r15
	fxcpmadd   5, 14,  1,  5

	fxcxnsma   3, 12,  1,  3
	lfpdx	17, %r27, %r16
	fxcxnsma   4, 13,  1,  4

	fxcxnsma   5, 14,  1,  5
	lfpdx	2, %r19, %r10
	dcbt	 %r10,	 %r29
	lfpdx	9, %r28, %r8
	mr	 %r24,	 %r26
	fxcpmadd   3, 15,  2,  3

	fxcpmadd   4, 16,  2,  4
	lfpdx	10, %r28, %r9
	fxcpmadd   5, 17,  2,  5

	fxcxnsma   3, 15,  2,  3
	lfpdx	11, %r28, %r10
	fxcxnsma   4, 16,  2,  4
	dcbt	 %r12,	 %r29
	fxcxnsma   5, 17,  2,  5
	lfpdx	0, %r20, %r8

	mr	 %r19,	 %r21
	mr	 %r27,	 %r29
	subic.	 3,	3,	1
	lwz	%r22, 36(%r4)
	fxpmul   6,  9,  0
	lwz	%r26, 40(%r4)
	lfpdx	12, %r28, %r11
	lwz	%r30, 44(%r4)
	fxpmul   7, 10,  0
	add    %r22, %r22, %r6

	add    %r26, %r26, %r7
	fxpmul   8, 11,  0
	add    %r30, %r30, %r5
	stfpdx	3, %r23, %r8
	fxcxnsma   6,  9,  0,  6
	dcbt	 %r8,	 %r21
	fxcxnsma   7, 10,  0,  7
	lfpdx	13, %r28, %r12
	fxcxnsma   8, 11,  0,  8
	dcbt	 %r10,	 %r21
	stfpdx	4, %r23, %r9

	lfpdx	14, %r28, %r13
	dcbt	 %r14,	 %r29
	stfpdx	5, %r23, %r10

	lfpdx	1, %r20, %r9

	stfpdx	5, %r23, %r10
	dcbt	 %r16,	 %r29


	dcbt	 %r16,	 %r4
	dcbt	 %r17,	 %r4

	dcbt	 %r8,	 %r30
	fxcpmadd   6, 12,  1,  6
	lfpdx	15, %r28, %r14
	fxcpmadd   7, 13,  1,  7

	fxcpmadd   8, 14,  1,  8
	lfpdx	16, %r28, %r15
	fxcxnsma   6, 12,  1,  6

	fxcxnsma   7, 13,  1,  7
	lfpdx	17, %r28, %r16
	fxcxnsma   8, 14,  1,  8
	dcbt	 %r10,	 %r30
	lfpdx	2, %r20, %r10

	lfpdx	9, %r29, %r8
	mr	 %r23,	 %r25
	fxcpmadd   6, 15,  2,  6
	addi %r4, %r4, 24

	mr	 %r20,	 %r22
	fxcpmadd   7, 16,  2,  7
	mr	 %r28,	 %r30
	lfpdx	10, %r29, %r9

	fxcpmadd   8, 17,  2,  8

	dcbt	 %r12,	 %r30

	fxcxnsma   6, 15,  2,  6

	lfpdx	11, %r29, %r10

	fxcxnsma   7, 16,  2,  7



	fxcxnsma   8, 17,  2,  8
	lfpdx	0, %r21, %r8

	dcbt	 %r8,	 %r22
	dcbt	 %r10,	 %r22

	dcbt	 %r14,	 %r30


	dcbt	 %r16,	 %r30


	bf gt,	pt_cmv_dag_pad_lab2
	b  pt_cmv_dag_pad_lab1
pt_cmv_dag_pad_lab2:
	fxpmul   3,  9,  0
	lfpdx	12, %r27, %r11
	fxpmul   4, 10,  0
	stfpdx	6, %r24, %r8
	fxpmul   5, 11,  0
	lfpdx	13, %r27, %r12
	fxcxnsma   3,  9,  0,  3
	stfpdx	7, %r24, %r9
	fxcxnsma   4, 10,  0,  4
	lfpdx	14, %r27, %r13
	fxcxnsma   5, 11,  0,  5
	stfpdx	8, %r24, %r10
	lfpdx	1, %r19, %r9
	stfpdx	8, %r24, %r10
	lfpdx	15, %r27, %r14
	mr	 %r20,	 %r22
	fxcpmadd   3, 12,  1,  3
	mr	 %r28,	 %r30
	lfpdx	16, %r27, %r15
	fxcpmadd   4, 13,  1,  4
	lfpdx	17, %r27, %r16
	fxcpmadd   5, 14,  1,  5
	lfpdx	2, %r19, %r10
	fxcxnsma   3, 12,  1,  3
	fxcxnsma   4, 13,  1,  4
	fxcxnsma   5, 14,  1,  5
	fxcpmadd   3, 15,  2,  3
	mr	 %r24,	 %r26
	lfpdx	9, %r28, %r8
	fxcpmadd   4, 16,  2,  4
	lfpdx	10, %r28, %r9
	fxcpmadd   5, 17,  2,  5
	lfpdx	11, %r28, %r10
	fxcxnsma   3, 15,  2,  3
	lfpdx	0, %r20, %r8
	fxcxnsma   4, 16,  2,  4
	fxcxnsma   5, 17,  2,  5
	lfpdx	12, %r28, %r11
	fxpmul   6,  9,  0
	stfpdx	3, %r23, %r8
	fxpmul   7, 10,  0
	lfpdx	13, %r28, %r12
	fxpmul   8, 11,  0
	stfpdx	4, %r23, %r9
	fxcxnsma   6,  9,  0,  6
	lfpdx	14, %r28, %r13
	fxcxnsma   7, 10,  0,  7
	stfpdx	5, %r23, %r10
	fxcxnsma   8, 11,  0,  8
	lfpdx	1, %r20, %r9
	stfpdx	5, %r23, %r10
	lfpdx	15, %r28, %r14
	fxcpmadd   6, 12,  1,  6
	lfpdx	16, %r28, %r15
	fxcpmadd   7, 13,  1,  7
	lfpdx	17, %r28, %r16
	fxcpmadd   8, 14,  1,  8
	lfpdx	2, %r20, %r10
	fxcxnsma   6, 12,  1,  6
	fxcxnsma   7, 13,  1,  7
	fxcxnsma   8, 14,  1,  8
	fxcpmadd   6, 15,  2,  6

	fxcpmadd   7, 16,  2,  7

	fxcpmadd   8, 17,  2,  8

	fxcxnsma   6, 15,  2,  6

	fxcxnsma   7, 16,  2,  7

	fxcxnsma   8, 17,  2,  8


	stfpdx	6, %r24, %r8
	stfpdx	7, %r24, %r9
	stfpdx	8, %r24, %r10
	stfpdx	8, %r24, %r10
pt_cmv_dag_pad_lab0:
	lwz	%r14,	368(%r1)
	lwz	%r15,	376(%r1)
	lwz	%r16,	384(%r1)
	lwz	%r17,	392(%r1)
	lwz	%r18,	400(%r1)
	lwz	%r19,	408(%r1)
	lwz	%r20,	416(%r1)
	lwz	%r21,	424(%r1)
	lwz	%r22,	432(%r1)
	lwz	%r23,	440(%r1)
	lwz	%r24,	448(%r1)
	lwz	%r25,	456(%r1)
	lwz	%r26,	464(%r1)
	lwz	%r27,	472(%r1)
	lwz	%r28,	480(%r1)
	lwz	%r29,	488(%r1)
	lwz	%r30,	496(%r1)
	lwz	%r31,	504(%r1)
	lfd	14,	112(%r1)
	lfd	15,	120(%r1)
	lfd	16,	128(%r1)
	lfd	17,	136(%r1)
	lfd	18,	144(%r1)
	lfd	19,	152(%r1)
	lfd	20,	160(%r1)
	lfd	21,	168(%r1)
	lfd	22,	176(%r1)
	lfd	23,	184(%r1)
	lfd	24,	192(%r1)
	lfd	25,	200(%r1)
	lfd	26,	208(%r1)
	lfd	27,	216(%r1)
	lfd	28,	224(%r1)
	lfd	29,	232(%r1)
	lfd	30,	240(%r1)
	lfd	31,	248(%r1)
	la   %r1,	864(%r1)
	blr
